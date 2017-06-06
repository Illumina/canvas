using Illumina.Common;
using Illumina.Common.FileSystem;
using Illumina.SecondaryAnalysis.VariantCalling;
using Isas.Framework.DataTypes;
using Isas.Framework.Logging;
using Isas.Framework.WorkManagement;
using Isas.Ploidy;
using Isas.SequencingFiles;
using Isas.SequencingFiles.Vcf;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Canvas.Wrapper
{
    [Obsolete("We plan to switch to providing (and consuming) the ploidy information in vcf format, as we do for Strelka and Manta")]
    public class CanvasPloidyBedCreator
    {
        private const string ExpectedSexChromosomeKaryotypeHeader = "##ExpectedSexChromosomeKaryotype";
        private const string PredictedSexChromosomeKaryotypeHeader = "##PredictedSexChromosomeKaryotype";
        public const string ReferenceSexChromosomeKaryotype = "##ReferenceSexChromosomeKaryotype";
        private readonly IWorkManager _workManager;
        private readonly ILogger _logger;
        private readonly PloidyCorrector _ploidyFixer;

        public CanvasPloidyBedCreator(ILogger logger, IWorkManager workManager, PloidyCorrector ploidyFixer)
        {
            _logger = logger;
            _workManager = workManager;
            _ploidyFixer = ploidyFixer;
        }

        /// <summary>
        ///  Write out the ploidy bed file if ploidy information is available from the vcf header
        /// Only create the normal XX or XY ploidy bed file so that Canvas can properly classify any abnormalities as variant. 
        /// If ploidy Y is > 1 produce the XY ploidy bed file, otherwise produce the XX ploidy bed file
        /// </summary>
        public IFileLocation CreateGermlinePloidyBed(SexPloidyInfo sexPloidy, GenomeMetadata genomeMetadata, IDirectoryLocation sampleSandbox)
        {
            if (sexPloidy.PloidyY > 0)
            {
                sexPloidy = SexPloidyInfo.NormalMale;
            }
            else
            {
                sexPloidy = SexPloidyInfo.NormalFemale;
            }
            return CreatePloidyBed(sexPloidy, genomeMetadata, sampleSandbox);
        }

        public IFileLocation CreatePloidyBed(SexPloidyInfo sexPloidy, GenomeMetadata genomeMetadata, IDirectoryLocation sampleSandbox)
        {
            var ploidyInfo = new SamplePloidyInfo();
            ploidyInfo.ProvidedPloidy = sexPloidy;
            IFileLocation ploidyBed = sampleSandbox.GetFileLocation("ploidy.bed.gz");
            var sexKaryotype = PloidyCorrector.PrettyPrintPloidy(sexPloidy.PloidyX, sexPloidy.PloidyY);
            _logger.Info($"Creating sex chromosome ploidy bed file at '{ploidyBed}' for sex karyotype '{sexKaryotype}'");
            string headerLine = $"{PloidyCorrector.ReferenceSexChromosomeKaryotype}={sexKaryotype}";
            this.WritePloidyBedFile(ploidyInfo, genomeMetadata, _ploidyFixer.GetParRegions(genomeMetadata),
                ploidyBed.FullName, headerLine, ploidy => true);
            return ploidyBed;
        }

        private static void WritePloidyBedLine(BgzipOrStreamWriter writer, string chr, long start, long end, int ploidy, Func<int, bool> includePloidy)
        {
            if (includePloidy(ploidy))
                writer.WriteLine($"{chr}\t{start}\t{end}\tploidy\t{ploidy}");
        }

        private void WritePloidyBedFile(
            SamplePloidyInfo ploidyInfo,
            GenomeMetadata genome,
            Dictionary<string, List<Interval>> parIntervals,
            BgzipOrStreamWriter writer, Func<int, bool> includePloidy)
        {
            var sampleInfo = new SampleInfo("none", "none");
            var ploidies = new SampleSet<SamplePloidyInfo> { [sampleInfo] = ploidyInfo }.SelectData(info => info.GetChromosomePloidySettings(genome));
            InternalWritePloidyLine writePloidyLine =
                (chr, start, end, samplePloidy) => WritePloidyBedLine(writer, chr, start, end, samplePloidy.Single().Value.Value, includePloidy);
            WritePloidyFile(ploidies, genome, parIntervals, writePloidyLine);
        }

        [Obsolete("Use WritePloidyVcf")]
        private void WritePloidyBedFile(
            SamplePloidyInfo ploidyInfo,
            GenomeMetadata genome,
            Dictionary<string, List<Interval>> parIntervals,
            string bedPath, string headerLine, Func<int, bool> includePloidy)
        {
            using (BgzipOrStreamWriter writer = new BgzipOrStreamWriter(bedPath))
            {
                if (!string.IsNullOrEmpty(headerLine)) writer.WriteLine(headerLine); // Include ##ReferenceSexChromosomeKaryotype line!
                WritePloidyBedFile(ploidyInfo, genome, parIntervals, writer, includePloidy);
            }
            //_tabixWrapper.BuildTabixIndex(new FileLocation(bedPath), TabixFileType.Bed);
        }

        private int GetChromosomeIndex(GenomeMetadata genome, string chr)
        {
            return genome.Sequences.FindIndex(sequence => sequence.Name == chr);
        }

        // TODO: remove this when removing WritePloidyLine delegate since in the future we will only write vcf lines
        private delegate void InternalWritePloidyLine(string chromosome, long zeroBasedStartInclusive, long zeroBasedEndExclusive, SampleSet<int?> ploidy);

        private bool DifferentPloidies(SampleSet<int?> ploidiesA, SampleSet<int?> ploidiesB)
        {
            return ploidiesA.Except(ploidiesB).Any();
        }

        private static int? CalculateParPloidy(Dictionary<string, int> chromosomePloidySettings, GenomeMetadata genome)
        {
            GenomeMetadata.SequenceMetadata chrX;
            if (!genome.TryGetChromosomeX(out chrX))
                return null;
            GenomeMetadata.SequenceMetadata chrY;
            if (!genome.TryGetChromosomeY(out chrY))
                return null;
            if (!chromosomePloidySettings.ContainsKey(chrX.Name) &&
                !chromosomePloidySettings.ContainsKey(chrY.Name))
            {
                return null;
            }
            int parPloidyX;
            int parPloidyY;
            if (!chromosomePloidySettings.TryGetValue(chrX.Name, out parPloidyX) ||
                !chromosomePloidySettings.TryGetValue(chrY.Name, out parPloidyY))
                throw new IlluminaException("Both ploidy X and ploidy Y must be specified. Only one was specified.");
            return parPloidyX + parPloidyY;
        }

        private void WritePloidyFile(
            SampleSet<Dictionary<string, int>> plodiesPerChromsome,
            GenomeMetadata genome,
            Dictionary<string, List<Interval>> parIntervals,
            InternalWritePloidyLine writePloidyLine)
        {
            var parPloidies = plodiesPerChromsome.SelectData(chromosomePloidy => CalculateParPloidy(chromosomePloidy, genome));
            var chromosomeKeys = plodiesPerChromsome.SampleData.SelectMany(chromosomePloidy => chromosomePloidy.Keys).Distinct().ToList();
            chromosomeKeys.Sort((strA, strB) => GetChromosomeIndex(genome, strA).CompareTo(GetChromosomeIndex(genome, strB)));
            foreach (string chromosome in chromosomeKeys)
            {
                var chromosomePloidies = plodiesPerChromsome.SelectData(chromosomePloidy => chromosomePloidy.ContainsKey(chromosome) ? chromosomePloidy[chromosome] : (int?)null);
                GenomeMetadata.SequenceMetadata sequence = genome.GetSequence(chromosome);

                // possibly multiple entries per chromosome to handle the PAR intervals
                int currentStart = 0;
                SampleSet<int?> currentPloidies = chromosomePloidies;
                var chromosomeParIntervals = parIntervals.ContainsKey(chromosome) ? parIntervals[chromosome] : new List<Interval>();
                foreach (var parInterval in chromosomeParIntervals)
                {
                    if (DifferentPloidies(currentPloidies, parPloidies))
                    {
                        if (parInterval.OneBasedStart > currentStart)
                            writePloidyLine(chromosome, currentStart, parInterval.OneBasedStart, currentPloidies);
                        currentStart = parInterval.OneBasedStart;
                        currentPloidies = parPloidies;
                    }
                    if (DifferentPloidies(currentPloidies, chromosomePloidies))
                    {
                        writePloidyLine(chromosome, currentStart, parInterval.OneBasedEnd, currentPloidies);
                        currentStart = parInterval.OneBasedEnd;
                        currentPloidies = chromosomePloidies;
                    }
                }
                if (currentStart < sequence.Length)
                    writePloidyLine(chromosome, currentStart, sequence.Length, chromosomePloidies);
            }
        }

        private static SexPloidyInfo GetSexPloidy(string sexChromosomeKaryotype)
        {
            var karyotypeLower = sexChromosomeKaryotype.ToLowerInvariant();
            var ploidyX = karyotypeLower.Count(letter => letter == 'x');
            var ploidyY = karyotypeLower.Count(letter => letter == 'y');
            return new SexPloidyInfo(ploidyX, ploidyY);
        }
        public IFileLocation GeneratePloidyBedFileFromSexChromosomeKaryotype(GenomeMetadata genome, string sexChromosomeKaryotype, IDirectoryLocation sampleSandbox)
        {
            var sexPloidy = GetSexPloidy(sexChromosomeKaryotype);
            return CreatePloidyBed(sexPloidy, genome, sampleSandbox);
        }
    }
}