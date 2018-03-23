using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using CanvasCommon;
using Illumina.Common;
using Isas.SequencingFiles;
using Isas.SequencingFiles.Vcf;
using Illumina.Common.FileSystem;
using Illumina.SecondaryAnalysis.VariantCalling;
using Isas.ClassicBioinfoTools.Tabix;
using Isas.Framework.DataTypes;
using Isas.Framework.Logging;
using Isas.Framework.Settings;
using Isas.Framework.WorkManagement;
using Isas.Framework.WorkManagement.CommandBuilding;

namespace EvaluateCNV
{
    public class CNInterval
    {
        public string Chromosome { get; }

        /// <summary>
        /// 0-based inclusive
        /// </summary>
        public int Start;
        /// <summary>
        /// 0-based exclusive
        /// </summary>
        public int End;
        public int Cn;
        public int? ReferenceCopyNumber; // set based on ploidy from GT fields or command line parameter
        public int BasesCovered;
        public int BasesExcluded;
        public int BasesCalledCorrectly;

        public int BasesNotCalled => Length - BasesExcluded - BasesCalledCorrectly - BasesCalledIncorrectly;

        public CNInterval(string chromosome, int start, int end, int cn)
        {
            Chromosome = chromosome;
            Start = start;
            End = end;
            Cn = cn;
        }
        public CNInterval(string chromosome)
        {
            Chromosome = chromosome;
        }

        public int Length => End - Start;

        public int BasesCalledIncorrectly;

        public override string ToString()
        {
            return $"{Chromosome}:{Start + 1}-{End}";
        }

        public void InitializeInterval()
        {
            BasesCovered = 0;
            BasesExcluded = 0;
            BasesCalledCorrectly = 0;
            BasesCalledIncorrectly = 0;
        }
    }

    public class CnvCall
    {
        public int Length => End - Start;

        public string Chr { get; }

        public int Start { get; }

        public int End { get; }

        public int CN { get; }

        public int? RefPloidy { get; }

        public bool PassFilter { get; }

        public string AltAllele { get; }

        public CnvCall(string chr, int start, int end, int cn, int? refPloidy, bool passFilter, string altAllele)
        {
            Chr = chr;
            Start = start;
            End = end;
            CN = cn;
            RefPloidy = refPloidy;
            PassFilter = passFilter;
            AltAllele = altAllele;
        }

        public bool IsAltVariant => CN != RefPloidy;

        public override string ToString()
        {
            return $"{Chr}:{Start}-{End} CN={CN}";
        }
    }

    public class CNVChecker
    {
        #region Members
        public Dictionary<string, List<CNInterval>> RegionsOfInterest = new Dictionary<string, List<CNInterval>>();
        public Dictionary<string, List<CNInterval>> ExcludeIntervals;
        private readonly IPloidyCorrector _ploidyCorrector;
        public double? DQscoreThreshold { get; }

        public CNVChecker(double? dQscoreThreshold, Dictionary<string, List<CNInterval>> excludeIntervals, IPloidyCorrector ploidyCorrector)
        {
            DQscoreThreshold = dQscoreThreshold;
            ExcludeIntervals = excludeIntervals;
            _ploidyCorrector = ploidyCorrector;
        }
        #endregion


        /// <summary>
        /// Load known CN data from a .bed file.  File lines have fields:
        /// chromosome, start, end, chromcountA, chromcountB
        /// So, copy number is the sum of the last 2 fields, major chromosome count is the max of the last 2 fields.
        /// </summary>
        /// <param name="oracleBedPath"></param>
        /// <param name="getCn"></param>
        /// <param name="heterogeneityFraction"></param>
        protected static Dictionary<string, List<CNInterval>> LoadIntervalsFromBed(string oracleBedPath, bool getCn, double heterogeneityFraction)
        {
            bool stripChr = false;
            int count = 0;
            long totalBases = 0;
            Dictionary<string, List<CNInterval>> bedIntervals = new Dictionary<string, List<CNInterval>>();
            using (FileStream stream = new FileStream(oracleBedPath, FileMode.Open, FileAccess.Read))
            using (StreamReader reader = new StreamReader(stream))
            {
                while (true)
                {
                    string fileLine = reader.ReadLine();
                    if (fileLine == null) break;
                    if (fileLine.Length == 0 || fileLine[0] == '#') continue;
                    string[] bits = fileLine.TrimEnd('\t').Split('\t');
                    if (bits.Length < 3) continue;
                    string chromosome = bits[0];
                    if (stripChr) chromosome = chromosome.Replace("chr", "");
                    if (!bedIntervals.ContainsKey(chromosome)) bedIntervals[chromosome] = new List<CNInterval>();
                    CNInterval interval = new CNInterval(chromosome);
                    interval.Start = int.Parse(bits[1]);
                    interval.End = int.Parse(bits[2]);
                    if (getCn) // bits.Length >= 5)
                    {
                        if (heterogeneityFraction < 1 && bits.Length > 5 && int.Parse(bits[3]) == 1 && int.Parse(bits[4]) == 1)
                            if (heterogeneityFraction > double.Parse(bits[5]))
                                continue;
                        interval.Cn = int.Parse(bits[3]) + int.Parse(bits[4]);
                    }
                    totalBases += interval.Length;
                    bedIntervals[chromosome].Add(interval);
                    count++;
                }
            }
            Console.WriteLine(">>>Loaded {0} CN intervals ({1} bases)", count, totalBases);
            return bedIntervals;
        }

        protected static Dictionary<string, List<CNInterval>> LoadKnownCNVCF(string oracleVcfPath)
        {
            bool stripChr = false;
            var knownCn = new Dictionary<string, List<CNInterval>>();
            // Load our "oracle" of known copy numbers:
            int count = 0;
            using (GzipReader reader = new GzipReader(oracleVcfPath))
            {
                while (true)
                {
                    string fileLine = reader.ReadLine();
                    if (fileLine == null) break;
                    if (fileLine.Length == 0 || fileLine[0] == '#') continue;
                    string[] bits = fileLine.Split('\t');
                    string chromosome = bits[0];
                    if (stripChr) chromosome = chromosome.Replace("chr", "");
                    if (!knownCn.ContainsKey(chromosome)) knownCn[chromosome] = new List<CNInterval>();
                    CNInterval interval = new CNInterval(chromosome)
                    {
                        Start = int.Parse(bits[1]),
                        Cn = -1
                    };
                    string[] infoBits = bits[7].Split(';');
                    foreach (string subBit in infoBits)
                    {
                        if (subBit.StartsWith("CN="))
                        {
                            float tempCn = float.Parse(subBit.Substring(3));
                            if (subBit.EndsWith(".5"))
                            {
                                interval.Cn = (int)Math.Round(tempCn + 0.1); // round X.5 up to X+1
                            }
                            else
                            {
                                interval.Cn = (int)Math.Round(tempCn); // Round off
                            }
                        }
                        if (subBit.StartsWith("END="))
                        {
                            interval.End = int.Parse(subBit.Substring(4));
                        }
                    }
                    // Parse CN from Canvas output:
                    if (bits.Length > 8)
                    {
                        string[] subBits = bits[8].Split(':');
                        string[] subBits2 = bits[9].Split(':');
                        for (int subBitIndex = 0; subBitIndex < subBits.Length; subBitIndex++)
                        {
                            if (subBits[subBitIndex] == "CN")
                            {
                                interval.Cn = int.Parse(subBits2[subBitIndex]);
                            }
                        }
                    }
                    if (interval.End == 0 || interval.Cn < 0)
                    {
                        Console.WriteLine("Error - bogus record!");
                        Console.WriteLine(fileLine);
                    }
                    else
                    {
                        knownCn[chromosome].Add(interval);
                        count++;
                    }
                }
            }
            Console.WriteLine(">>>Loaded {0} known-CN intervals", count);
            return knownCn;
        }

        protected void LoadRegionsOfInterest(string bedPath)
        {
            if (string.IsNullOrEmpty(bedPath)) return;
            if (!File.Exists(bedPath))
            {
                throw new ArgumentException(string.Format("* Error: ROI bed file not found at '{0}'", bedPath));
            }
            RegionsOfInterest = LoadIntervalsFromBed(bedPath, false, 1.0);
            var keys = RegionsOfInterest.Keys.ToList();
            foreach (string key in keys)
            {
                RegionsOfInterest[string.Format("chr{0}", key)] = RegionsOfInterest[key];
            }
        }

        protected static Dictionary<string, List<CNInterval>> LoadKnownCn(string oraclePath, double heterogeneityFraction)
        {
            if (!File.Exists(oraclePath))
            {
                throw new ArgumentException(string.Format("* Error: Truth vcf not found at '{0}'", oraclePath));
            }
            if (oraclePath.EndsWith(".bed"))
                return LoadIntervalsFromBed(oraclePath, true, heterogeneityFraction);
            var knownCn = LoadKnownCNVCF(oraclePath);
            SummarizeTruthSetStatistics(knownCn);
            return knownCn;
        }

        public void InitializeIntervalMetrics(Dictionary<string, List<CNInterval>> knownCN)
        {
            foreach (var chromosomeIntervals in knownCN.Values)
                foreach (var interval in chromosomeIntervals)
                    interval.InitializeInterval();
        }

        protected static void SummarizeTruthSetStatistics(Dictionary<string, List<CNInterval>> knownCn)
        {
            List<long> eventSizes = new List<long>();
            double meanEventSize = 0;
            int countUnder10KB = 0;
            int count10kb50kb = 0;
            int count50kb500kb = 0;
            int count500kbplus = 0;
            foreach (string key in knownCn.Keys)
            {
                foreach (CNInterval interval in knownCn[key])
                {
                    if (interval.Cn == 2) continue;
                    long length = interval.Length;
                    eventSizes.Add(length);
                    meanEventSize += length;
                    if (length <= 10000)
                    {
                        countUnder10KB++;
                    }
                    else if (length <= 50000)
                    {
                        count10kb50kb++;
                    }
                    else if (length <= 500000)
                    {
                        count50kb500kb++;
                    }
                    else
                    {
                        count500kbplus++;
                    }
                }
            }
            eventSizes.Sort();
            meanEventSize /= eventSizes.Count;
            Console.WriteLine("Known CN: Mean length {0:F4}", meanEventSize);
            Console.WriteLine("up to 10kb: {0}", countUnder10KB);
            Console.WriteLine("10kb-50kb: {0}", count10kb50kb);
            Console.WriteLine("50kb-500kb: {0}", count50kb500kb);
            Console.WriteLine("500kb+: {0}", count500kbplus);
            if (eventSizes.Count > 0)
                Console.WriteLine("Median size: {0}", eventSizes[eventSizes.Count / 2]);
        }

        protected static int GetCopyNumber(VcfVariant variant, out int end)
        {
            int CN = -1;
            end = -1;
            if (variant.GenotypeColumns != null && variant.GenotypeColumns.Count > 0)
            {
                Dictionary<string, string> genotype = variant.GenotypeColumns[variant.GenotypeColumns.Count - 1];
                if (genotype.ContainsKey("CN"))
                {
                    CN = int.Parse(genotype["CN"]);
                }
                if (genotype.ContainsKey("END"))
                {
                    end = int.Parse(genotype["END"]);
                }
            }
            if (variant.InfoFields.ContainsKey("END"))
            {
                end = int.Parse(variant.InfoFields["END"]);
            }
            if (variant.InfoFields.ContainsKey("CN"))
            {
                CN = int.Parse(variant.InfoFields["CN"]);
            }

            return CN;
        }

        private static int? TryGetRefPloidy(VcfVariant variant)
        {
            var genotypeColumn = variant.GenotypeColumns[variant.GenotypeColumns.Count - 1];
            if (!genotypeColumn.TryGetValue("GT", out var genotype)) return null;

            // genotype == "." could be reference copy number 0 (e.g. chrY for female sample) or 1
            // since we don't know return null. NB: if this causes a problem downstream the reference copy number should really be encoded in the truth file!
            if (genotype == ".") return null;

            return genotype.Split('/', '|').Length;
        }


        public void CountExcludedBasesInTruthSetIntervals(Dictionary<string, List<CNInterval>> knownCn)
        {
            foreach (string key in knownCn.Keys)
            {
                if (!ExcludeIntervals.ContainsKey(key)) continue;
                foreach (CNInterval interval in knownCn[key])
                {
                    foreach (CNInterval excludeInterval in ExcludeIntervals[key])
                    {
                        int overlapStart = Math.Max(interval.Start, excludeInterval.Start);
                        int overlapEnd = Math.Min(interval.End, excludeInterval.End);
                        if (overlapStart >= overlapEnd) continue;
                        interval.BasesExcluded += overlapEnd - overlapStart;
                        //Console.WriteLine("Interval {0}:{1}-{2} excludes {3} bases due to overlap with excluded interval {4}:{5}-{6}",
                        //    key, interval.Start, interval.End, overlapEnd - overlapStart,
                        //    key, excludeInterval.Start, excludeInterval.End);
                    }
                }
            }
        }

        private static List<string> GetVcfHeaderLines(IFileLocation vcfPath)
        {
            using (var reader = new VcfReader(vcfPath.FullName, false))
            {
                return reader.HeaderLines;
            }
        }

        private static void HandleHeaderLine(TextWriter writer, List<string> headerLines, string headerKey,
            Action<TextWriter, string> processValue)
        {
            if (headerLines.Any(stringToCheck => stringToCheck.Contains(headerKey)))
                processValue(writer, headerLines.Find(stringToCheck => stringToCheck.Contains(headerKey)));
        }

        private static void LogPurity(TextWriter writer, string value)
        {
            double purity = double.Parse(value.Split("=")[1]);
            writer.WriteLine($"Purity\t{purity}");
        }

        private static void LogPloidy(TextWriter writer, string value)
        {
            double ploidy = double.Parse(value.Split("=")[1]);
            writer.WriteLine($"Ploidy\t{ploidy}");
        }

        public void HandleVcfHeaderInfo(TextWriter outputWriter, IFileLocation vcfPath)
        {
            var headerLines = GetVcfHeaderLines(vcfPath);
            HandleHeaderLine(outputWriter, headerLines, "EstimatedTumorPurity", LogPurity);
            HandleHeaderLine(outputWriter, headerLines, "OverallPloidy", LogPloidy);
        }

        public static Dictionary<string, List<CnvCall>> GetCnvCallsFromVcf(string vcfPath, double? dQscoreThreshold)
        {
            var calls = new Dictionary<string, List<CnvCall>>();
            using (VcfReader reader = new VcfReader(vcfPath, false))
            {
                if (dQscoreThreshold.HasValue)
                {
                    var match = reader.HeaderLines.FirstOrDefault(stringToCheck => stringToCheck.Contains("DQ"));
                    if (match == null)
                        throw new ArgumentException($"File {vcfPath} does not contain DQ INFO field.");
                }

                foreach (VcfVariant variant in reader.GetVariants())
                {
                    if (!calls.ContainsKey(variant.ReferenceName)) calls[variant.ReferenceName] = new List<CnvCall>();
                    int end;
                    int cn = GetCopyNumber(variant, out end);
                    int? refPloidy = TryGetRefPloidy(variant);
                    var passFilter = variant.Filters == "PASS";
                    if (dQscoreThreshold.HasValue)
                    {
                        var genotypeColumn = variant.GenotypeColumns.Single();
                        if (!genotypeColumn.Keys.Contains("DQ"))
                            continue;
                        if (variant.Identifier.Contains("REF"))
                            continue;
                        if (genotypeColumn["DQ"] == ".")
                            continue;
                        if (Double.Parse(genotypeColumn["DQ"]) < dQscoreThreshold.Value)
                            continue;
                    }
                    calls[variant.ReferenceName].Add(new CnvCall(variant.ReferenceName, variant.ReferencePosition,
                        end, cn, refPloidy, passFilter, variant.VariantAlleles.First()));
                }
            }
            return calls;
        }

        public static void Evaluate(string truthSetPath, string cnvCallsPath, string excludedBed, string outputPath, EvaluateCnvOptions options)
        {
            double heterogeneityFraction = options.HeterogeneityFraction;
            var knownCn = LoadKnownCn(truthSetPath, heterogeneityFraction);
            knownCn = knownCn.SelectValues(
                truthEntries => truthEntries.Where(truthEntry => truthEntry.Length >= options.MinEntrySize).ToList());
            var calls = GetCnvCallsFromVcf(cnvCallsPath, options.DQscoreThreshold);
            calls = calls.SelectValues(
                chromosomeCalls => chromosomeCalls.Where(call => call.Length >= options.MinEntrySize).ToList());

            // LoadRegionsOfInterest(options.RoiBed?.FullName);
            var excludeIntervals = new Dictionary<string, List<CNInterval>>();
            if (!string.IsNullOrEmpty(excludedBed))
            {
                var excludeIntervalsTmp = LoadIntervalsFromBed(excludedBed, false, 1.0);
                List<string> keys = excludeIntervalsTmp.Keys.ToList();
                foreach (string key in keys)
                {
                    string chr = key;
                    if (!calls.ContainsKey(chr)) chr = key.Replace("chr", "");
                    if (!calls.ContainsKey(chr)) chr = "chr" + key;
                    if (!calls.ContainsKey(chr))
                    {
                        Console.WriteLine($"Error: Skipping exclude intervals for chromosome {key} with no truth data." +
                                          $"Check that chromosome names are spelled correctly for exclude intervals");
                        continue;
                    }
                    excludeIntervals[chr] = excludeIntervalsTmp[key];
                }
            }
            Console.WriteLine("TruthSet\t{0}", truthSetPath);
            Console.WriteLine("CNVCalls\t{0}", cnvCallsPath);

            bool includePassingOnly = Path.GetFileName(cnvCallsPath).ToLower().Contains("vcf");
            var logger = new Logger(new[] { Console.Out }, new[] { Console.Error });
            var settings = IsasConfigurationSettings.GetConfigSettings();
            var output = new DirectoryLocation(outputPath);
            var workerDirectory = new DirectoryLocation(Isas.Framework.Utilities.Utilities.GetAssemblyFolder(typeof(CNVChecker)));
            var commandManager = new CommandManager(new ExecutableProcessor(settings, logger, workerDirectory));
            WorkDoerFactory.RunWithWorkDoer(logger, settings, output, workDoer =>
            {
                var tabixWrapper = TabixWrapperFactory.GetTabixWrapper(logger, workDoer, commandManager);
                var ploidyCorrector = new PloidyCorrector(logger, workDoer,
                    new PloidyEstimator(logger, workDoer, null, false, commandManager), tabixWrapper, false);
                var checker = new CNVChecker(options.DQscoreThreshold, excludeIntervals, ploidyCorrector);
                if (options.PloidyInfo.SexPloidyInfo != null)
                {
                    Console.WriteLine($">>>Getting reference ploidy from provided ploidy information and PAR bed file '{options.PloidyInfo.ParBed}'");

                    var ploidy = checker.GetPloidy(options.PloidyInfo, output);
                    var referencePloidy = LoadReferencePloidy(options.PloidyInfo.SexPloidyInfo, options.PloidyInfo.ParBed);
                    knownCn = GetKnownCopyNumberWithReferencePloidy(referencePloidy, knownCn);
                    calls = GetCallsWithRefPloidy(calls, ploidy);
                }
                var cnvEvaluator = new CnvEvaluator(checker);

                if (checker.DQscoreThreshold.HasValue && !Path.GetFileName(cnvCallsPath).ToLower().Contains("vcf"))
                    throw new ArgumentException("CNV.vcf must be in a vcf format when --dqscore option is used");
                cnvEvaluator.ComputeAccuracy(knownCn, cnvCallsPath, outputPath, includePassingOnly, options, calls);
                if (includePassingOnly)
                    cnvEvaluator.ComputeAccuracy(knownCn, cnvCallsPath, outputPath, false, options, calls);
                Console.WriteLine(">>>Done - results written to {0}", outputPath);
            });
        }

        private static Dictionary<string, List<CNInterval>> GetKnownCopyNumberWithReferencePloidy(ReferencePloidy referencePloidy, Dictionary<string, List<CNInterval>> knownCn)
        {
            return knownCn.SelectValues(knownCnIntervals =>
                GetKnownCopyNumberWithReferencePloidy(knownCnIntervals, referencePloidy).ToList());
        }

        private static IEnumerable<CNInterval> GetKnownCopyNumberWithReferencePloidy(List<CNInterval> knownCnIntervals, ReferencePloidy referencePloidy)
        {
            foreach (var knownCnInterval in knownCnIntervals)
            {
                var interval = new ReferenceInterval(knownCnInterval.Chromosome, new Interval(knownCnInterval.Start + 1, knownCnInterval.End));
                var refPloidy = referencePloidy.GetSingleReferencePloidy(interval);
                yield return new CNInterval(knownCnInterval.Chromosome,
                    knownCnInterval.Start, knownCnInterval.End, knownCnInterval.Cn)
                {
                    ReferenceCopyNumber = refPloidy
                };
            }
        }

        private PloidyInfo GetPloidy((SexPloidyInfo SexPloidyInfo, IFileLocation ParBed) ploidyInfo, IDirectoryLocation output)
        {
            var vcf = new Vcf(output.GetFileLocation("ploidy.vcf.gz"));
            var genome = new ReferenceGenome(ploidyInfo.ParBed.Directory.Parent).GenomeMetadata;
            var sample = new SampleSet<SexPloidyInfo>
            {
                {new SampleInfo("EvaluateCNVSample", "EvaluateCNVSample"), ploidyInfo.SexPloidyInfo}
            };
            _ploidyCorrector.WritePloidyVcfFile(vcf, sample, genome);
            return PloidyInfo.LoadPloidyFromVcfFileNoSampleId(vcf.VcfFile.FullName);
        }

        private static ReferencePloidy LoadReferencePloidy(SexPloidyInfo sexPloidyInfo, IFileLocation parBed)
        {
            var genome = new ReferenceGenome(parBed.Directory.Parent).GenomeMetadata;
            string ploidyVcfString;
            var sampleId = "EvaluateCNVSample";
            using (var stringWriter = new StringWriter())
            using (var writer = new BgzipOrStreamWriter(stringWriter))
            {
                var sample = new SampleSet<SexPloidyInfo>
                {
                    {new SampleInfo(sampleId, sampleId), sexPloidyInfo}
                };
                PloidyCorrector.WritePloidyVcfFile(writer, sample, genome, parBed);
                ploidyVcfString = stringWriter.ToString();
            }

            using (var stringReader = new StringReader(ploidyVcfString))
            using (var reader = new VcfReader(stringReader))
            {
                return ReferencePloidy.Load(reader, new SampleId(sampleId));
            }
        }

        private static Dictionary<string, List<CnvCall>> GetCallsWithRefPloidy(
            Dictionary<string, List<CnvCall>> allCalls, PloidyInfo ploidy)
        {
            return allCalls.Select(kvp =>
                    (kvp.Key, GetCallsWithRefPloidy(kvp.Value, ploidy.PloidyByChromosome.ContainsKey(kvp.Key) ? ploidy.PloidyByChromosome[kvp.Key] : new List<PloidyInterval>()).ToList()))
                .ToDictionary();
        }

        private static IEnumerable<CnvCall> GetCallsWithRefPloidy(List<CnvCall> calls, List<PloidyInterval> ploidyRegions)
        {
            foreach (var call in calls)
            {
                // any truth entries not in ploidy regions are considered diploid
                int refPloidy = TryGetPloidyForCall(call, ploidyRegions, out var ploidyRegion) ? ploidyRegion.Ploidy : 2;
                yield return new CnvCall(call.Chr, call.Start, call.End, call.CN, refPloidy, call.PassFilter, call.AltAllele);
            }
        }

        private static bool TryGetPloidyForCall(CnvCall call, List<PloidyInterval> ploidyRegions,
            out PloidyInterval ploidyRegion)
        {
            foreach (var currentPloidyRegion in ploidyRegions)
            {
                // interval must be completely contained within the ploidy region	
                if (call.Start >= currentPloidyRegion.Start && call.End <= currentPloidyRegion.End)
                {

                    ploidyRegion = currentPloidyRegion;
                    return true;
                }
            }
            // majority of the variant falls into ploidy region
            foreach (var currentPloidyRegion in ploidyRegions)
            {
                // for now assign variant that overlaps several ploidy regions to the one with the highest overlap 
                int overlapStart = Math.Max(call.Start, currentPloidyRegion.Start);
                int overlapEnd = Math.Min(call.End, currentPloidyRegion.End);
                if (overlapStart >= overlapEnd) continue;
                double overlapRatio = (overlapEnd - overlapStart) / (double)call.Length;
                if (overlapRatio > 0.5 && call.RefPloidy.HasValue)
                {
                    ploidyRegion = currentPloidyRegion;
                    return true;
                }
            }
            ploidyRegion = null;
            var refPloidy = 2;
            if (call.RefPloidy.HasValue && call.RefPloidy.Value != refPloidy)
                throw new IlluminaException(
                    $"call '{call}' had unexpected reference ploidy '{call.RefPloidy}'. From provided ploidy information expected '{refPloidy}'");
            return false;
        }
    }
}
