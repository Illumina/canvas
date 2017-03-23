using System;
using System.IO;
using System.Collections;
using System.Diagnostics;
using System.Collections.Generic;
using System.Threading;
using ProtoBuf;
using System.Linq;
using System.Threading.Tasks;
using Illumina.Common.FileSystem;
using Isas.Manifests.NexteraManifest;
using Isas.SequencingFiles;
using CanvasCommon;

namespace CanvasBin
{

    public class SampleHitArrays
    {
        private readonly Dictionary<string, HitArray> _observedAlignments;

        public SampleHitArrays(Dictionary<string, HitArray> observedAlignments)
        {
            _observedAlignments = observedAlignments;
        }

        /// <summary>
        /// Returns ratios of observed and possible alignments for each chromosome
        /// </summary>
        public List<double> GetRates(Dictionary<string, BitArray> possibleAlignments, NexteraManifest manifest = null)
        {
            List<double> rates = new List<double>();

            Dictionary<string, List<NexteraManifest.ManifestRegion>> manifestRegionsByChrom = null;
            if (manifest != null)
            {
                manifestRegionsByChrom = manifest.GetManifestRegionsByChromosome();
            }

            var tasks = new List<Action>();
            foreach (string chr in possibleAlignments.Keys)
            {
                // We don't want to include the sex chromosomes because they may not be copy number 2
                if (!GenomeMetadata.SequenceMetadata.IsAutosome(chr))
                    continue;
                HitArray observed = _observedAlignments[chr];
                BitArray possible = possibleAlignments[chr];
                List<NexteraManifest.ManifestRegion> regions = null;
                if (manifestRegionsByChrom != null)
                {
                    if (!manifestRegionsByChrom.ContainsKey(chr)) { continue; }
                    regions = manifestRegionsByChrom[chr];
                }

                tasks.Add(() =>
                {
                    int numberObserved = observed.CountSetBits(regions);
                    int numberPossible = CanvasBin.CountSetBits(possible, regions);

                    double rate = numberObserved / (double)numberPossible;

                    lock (rates)
                    {
                        rates.Add(rate);
                    }

                });
            }
            Parallel.ForEach(tasks, task => task.Invoke());
            return rates;
        }

        public int GetBinSize(int countsPerBin, Dictionary<string, BitArray> possibleAlignments, NexteraManifest manifest)
        {
            List<double> rates = GetRates(possibleAlignments, manifest);
            return GetBinSize(countsPerBin, rates);
        }

        public static int GetBinSize(int countsPerBin, List<double> rates)
        {
            double medianRate = Utilities.Median(rates);
            return (int)(countsPerBin / medianRate);
        }
    }

    public class MultiSampleHitArrays
    {
        private readonly List<SampleHitArrays> _sampleHitArrays;
        public MultiSampleHitArrays(List<SampleHitArrays> sampleHitArrays)
        {
            _sampleHitArrays = sampleHitArrays;
        }

        public List<double> GetRates(Dictionary<string, BitArray> possibleAlignments, NexteraManifest manifest)
        {
            List<double> rates = new List<double>();

            foreach (SampleHitArrays sampleHitArrays in _sampleHitArrays)
            {
                rates.AddRange(sampleHitArrays.GetRates(possibleAlignments, manifest));
            }
            return rates;
        }

        public int GetBinSize(int countsPerBin, Dictionary<string, BitArray> possibleAlignments, NexteraManifest manifest)
        {
            var rates = GetRates(possibleAlignments, manifest);
            return SampleHitArrays.GetBinSize(countsPerBin, rates);
        }
    }

    public class CanvasBin
    {
        private static readonly int numberOfGCbins = 101;

        /// <summary>
        /// Counts the number of 'on' bits in a BitArray.
        /// </summary>
        /// <param name="bits">BitArray to count the bits of</param>
        /// <returns>number of 'on' bits</returns>
        public static int CountSetBits(BitArray bits, List<NexteraManifest.ManifestRegion> regions)
        {
            if (regions == null) { return CountSetBits(bits); }

            int total = 0;
            int i = -1;
            foreach (var region in regions)
            {
                if (i < region.Start) // avoid overlapping targeted regions
                {
                    i = region.Start - 1; // i is 0-based; manifest coordinates are 1-based.
                }
                for (; i < bits.Length && i < region.End; i++)
                {
                    if (bits[i]) { total++; }
                }
            }

            return total;
        }

        /// <summary>
        /// Counts the number of 'on' bits in a BitArray.
        /// </summary>
        /// <param name="bits">BitArray to count the bits of</param>
        /// <returns>number of 'on' bits</returns>
        static int CountSetBits(BitArray bits)
        {
            int total = 0;

            for (int i = 0; i < bits.Length; i++)
            {
                if (bits[i])
                    total++;
            }
            return total;
        }

        /// <summary>
        /// Estimate mean fragment length of read pairs.
        /// </summary>
        /// <param name="c">Dictionary of fragment length arrays for each chromosome.</param>
        /// <returns>Mean fragment length.</returns>
        static Int16 MeanFragmentSize(Dictionary<string, Int16[]> fragmentLengths)
        {
            List<Int16> fragmentMeanLengths = new List<Int16>();

            foreach (KeyValuePair<string, Int16[]> kvp in fragmentLengths)
            {
                fragmentMeanLengths.Add(CanvasCommon.Utilities.NonZeroMean(kvp.Value));
            }
            return CanvasCommon.Utilities.NonZeroMean(fragmentMeanLengths.ToArray());

        }

        /// <summary>
        /// Scale fragment sizes to a 0 .. 255 byte range
        /// </summary>
        /// <param name="fragmentLength">Int fragment length.</param>
        /// <param name="max">Maximal input fragment length.</param>
        /// <param name="min">Minimal input  fragment length.</param>
        /// <returns>Fragment length scaled to 0 .. 255 byte range.</returns>
        static byte ScaleIntToByte(int fragmentLength, int max, int min)
        {
            double result = (fragmentLength - min) * 255 / (double)(max - min);
            if (result > 255) return 255;
            if (result < 0) return 0;
            return (byte)result;
        }

        /// <summary>
        /// Sets up two Dictionaries holding BitArrays, one BitArray for each chromosome in a fasta file. One bit for each nucleotide.
        /// </summary>
        /// <param name="fastaFile">Fasta file containing uniquemer-marked reference genome.</param>
        /// <param name="possibleAlignments">Stores which alignments are possible (perfect and unique).</param>
        /// <param name="observedAlignments">Stores observed alignments from a sample.</param>
        /// <param name="fragmentLengths">Stores fragment length (Int16).</param>
        static void InitializeAlignmentArrays(string fastaFile, string chromosome, CanvasCoverageMode coverageMode, IDictionary<string, BitArray> possibleAlignments, IDictionary<string, HitArray> observedAlignments, IDictionary<string, Int16[]> fragmentLengths)
        {
            string referenceBases = FastaLoader.LoadFastaSequence(fastaFile, chromosome);

            BitArray possible = new BitArray(referenceBases.Length);
            possibleAlignments[chromosome] = possible;
            observedAlignments[chromosome] = new HitArray(referenceBases.Length);
            if (coverageMode == CanvasCoverageMode.GCContentWeighted)
                fragmentLengths[chromosome] = new Int16[referenceBases.Length];
            else
                fragmentLengths[chromosome] = new Int16[0];
            // Mark which k-mers in the fasta file are unique. These are indicated by upper-case letters.
            for (int i = 0; i < referenceBases.Length; i++)
            {
                if (char.IsUpper(referenceBases[i]))
                    possible[i] = true;
            }
        }

        /// <summary>
        /// Reads in a bam file and marks within the BitArrays which genomic mers are present.
        /// </summary>
        /// <param name="bamFile">bam file read alignments from.</param>
        /// <param name="observedAlignments">Dictioanry of BitArrays, one for each chromosome, to store the alignments in.</param>
        static void LoadObservedAlignmentsBAM(string bamFile, bool isPairedEnd, string chromosome, CanvasCoverageMode coverageMode, HitArray observed, Int16[] fragmentLengths)
        {
            // Sanity check: The .bai file must exist, in order for us to seek to our target chromosome!
            string indexPath = bamFile + ".bai";
            if (!File.Exists(indexPath))
            {
                throw new Exception(string.Format("Fatal error: Bam index not found at {0}", indexPath));
            }

            using (BamReader reader = new BamReader(bamFile))
            {
                int desiredRefIndex = -1;
                desiredRefIndex = reader.GetReferenceIndex(chromosome);
                if (desiredRefIndex == -1)
                {
                    throw new Illumina.Common.IlluminaException(
                        string.Format("Unable to retrieve the reference sequence index for {0} in {1}.", chromosome,
                        bamFile));
                }
                bool result = reader.Jump(desiredRefIndex, 0);
                if (!result)
                {
                    // Note: This is not necessarily an error, it just means that there *are* no reads for this chromosome in this 
                    // .bam file.  That is not uncommon e.g. for truseq amplicon.
                    return;
                }
                int readCount = 0;
                int keptReadCount = 0;
                string header = reader.GetHeader();
                BamAlignment alignment = new BamAlignment();
                while (reader.GetNextAlignment(ref alignment, true))
                {
                    readCount++;

                    // Flag check - Require reads to be aligned, passing filter, non-duplicate:
                    if (!alignment.IsMapped()) continue;
                    if (alignment.IsFailedQC()) continue;
                    if (alignment.IsDuplicate()) continue;
                    if (alignment.IsReverseStrand()) continue;
                    if (!alignment.IsMainAlignment()) continue;

                    // Require the alignment to start with 35 bases of non-indel:
                    if (alignment.CigarData[0].Type != 'M' || alignment.CigarData[0].Length < 35) continue;

                    if (isPairedEnd && !alignment.IsProperPair()) continue;

                    int refID = alignment.RefID;

                    // quit if the current reference index is different from the desired reference index
                    if (refID != desiredRefIndex)
                        break;

                    if (refID == -1)
                        continue;

                    keptReadCount++;
                    if (coverageMode == CanvasCoverageMode.Binary)
                    {
                        observed.Data[alignment.Position] = 1;
                    }
                    else
                    {
                        observed.Set(alignment.Position);
                    }
                    // store fragment size, make sure it's within Int16 range and is positive (simplification for now)
                    if (coverageMode == CanvasCoverageMode.GCContentWeighted)
                        fragmentLengths[alignment.Position] = Convert.ToInt16(Math.Max(Math.Min(Int16.MaxValue, alignment.FragmentLength), 0));
                }
                Console.WriteLine("Kept {0} of {1} total reads", keptReadCount, readCount);
            }
        }

        /// <summary>
        /// Calculates how many possible alignments corresponds to the desired number of observed alignments per bin.
        /// </summary>
        /// <param name="countsPerBin">Desired number of observed alignments per bin.</param>
        /// <param name="possibleAlignments">BitArrays of possible alignments (unique mers).</param>
        /// <param name="observedAlignments">BitArrays storing the observed alignments.</param>
        /// <returns>Median alignment rate observed on the autosomes.</returns>
        static int CalculateNumberOfPossibleAlignmentsPerBin(int countsPerBin, Dictionary<string, BitArray> possibleAlignments,
            Dictionary<string, HitArray> observedAlignments, NexteraManifest manifest = null)
        {
            var sampleHitArrays = new SampleHitArrays(observedAlignments);
            return sampleHitArrays.GetBinSize(countsPerBin, possibleAlignments, manifest);
        }


        /// <summary>
        /// Calculates how many possible alignments corresponds to the desired number of observed alignments per bin.
        /// </summary>
        /// <param name="countsPerBin">Desired number of observed alignments per bin.</param>
        /// <param name="possibleAlignments">BitArrays of possible alignments (unique mers).</param>
        /// <param name="observedAlignments">BitArrays storing the observed alignments.</param>
        /// <returns>Median alignment rate observed on the autosomes.</returns>
        static int CalculateNumberOfPossibleMultiSampleAlignmentsPerBin(int countsPerBin, List<IntermediateDataHolder> intermediateData, NexteraManifest manifest = null)
        {
            List<SampleHitArrays> allSampleHitArrays = intermediateData.Select(data => new SampleHitArrays(data.ObservedAlignments)).ToList();
            var multiSampleHitArrays = new MultiSampleHitArrays(allSampleHitArrays);
            return multiSampleHitArrays.GetBinSize(countsPerBin, intermediateData.First().PossibleAlignments, manifest);
        }


        class BinState
        {
            public int NucleotideCount;
            public int GCCount;
            public int PossibleCount;
            public int ObservedCount;
            public int StartPosition = -1;

            public void Reset()
            {
                NucleotideCount = 0;
                GCCount = 0;
                PossibleCount = 0;
                ObservedCount = 0;
                StartPosition = -1;
            }
        }


        /// <summary>
        /// Computes fragment-based GC normalization correction factor 
        /// </summary>
        /// <returns>An array of observed vs expected GC counts.</returns>
        static float[] ComputeObservedVsExpectedGC(Dictionary<string, HitArray> observedAlignments,
            Dictionary<string, byte[]> readGCContent, NexteraManifest manifest,
            bool debugGC, string outFile)
        {

            Dictionary<string, List<NexteraManifest.ManifestRegion>> regionsByChrom = null;
            if (manifest != null)
            {
                regionsByChrom = manifest.GetManifestRegionsByChromosome();
            }

            long[] expectedReadCountsByGC = new long[numberOfGCbins];
            long[] observedReadCountsByGC = new long[numberOfGCbins];
            foreach (KeyValuePair<string, byte[]> chromosomeReadGCContent in readGCContent)
            {
                string chr = chromosomeReadGCContent.Key;
                if (!observedAlignments.ContainsKey(chr)) { continue; }

                if (manifest == null) // look at the entire genome
                {
                    for (int i = 0; i < chromosomeReadGCContent.Value.Length; i++)
                    {
                        expectedReadCountsByGC[chromosomeReadGCContent.Value[i]]++;
                        observedReadCountsByGC[chromosomeReadGCContent.Value[i]] += observedAlignments[chr].Data[i];
                    }
                }
                else // look at only the targeted regions
                {
                    if (!regionsByChrom.ContainsKey(chr)) { continue; }
                    int i = -1;
                    foreach (var region in regionsByChrom[chr])
                    {
                        if (i < region.Start) // avoid overlapping targeted regions
                        {
                            i = region.Start - 1; // i is 0-based; manifest coordinates are 1-based.
                        }
                        for (; i < chromosomeReadGCContent.Value.Length && i < region.End; i++)
                        {
                            expectedReadCountsByGC[chromosomeReadGCContent.Value[i]]++;
                            observedReadCountsByGC[chromosomeReadGCContent.Value[i]] += observedAlignments[chr].Data[i];
                        }
                    }
                }
            }

            // calculate ratio of observed to expected read counts for each read GC bin
            float[] observedVsExpectedGC = new float[numberOfGCbins];
            for (int i = 0; i < numberOfGCbins; i++)
                observedVsExpectedGC[i] = 1;
            long sumObserved = 0;
            long sumExpected = 0;
            foreach (long gcContent in observedReadCountsByGC)
                sumObserved += gcContent;
            foreach (long gcContent in expectedReadCountsByGC)
                sumExpected += gcContent;
            for (int binIndex = 0; binIndex < numberOfGCbins; binIndex++)
            {
                if (expectedReadCountsByGC[binIndex] == 0)
                    expectedReadCountsByGC[binIndex] = 1;
                if (observedReadCountsByGC[binIndex] == 0)
                    observedReadCountsByGC[binIndex] = 1;
                observedVsExpectedGC[binIndex] = ((float)observedReadCountsByGC[binIndex] / (float)expectedReadCountsByGC[binIndex]) * ((float)sumExpected / (float)sumObserved);
            }

            if (debugGC)
            {
                using (GzipWriter writer = new GzipWriter(outFile + ".gcstat"))
                {
                    for (int binIndex = 0; binIndex < numberOfGCbins; binIndex++)
                    {
                        writer.WriteLine(string.Format("{0}\t{1}\t{2}", expectedReadCountsByGC[binIndex], observedReadCountsByGC[binIndex], observedVsExpectedGC[binIndex]));
                    }
                }
            }
            return observedVsExpectedGC;
        }

        /// <summary>
        /// Bin alignments.
        /// </summary>
        /// <param name="referenceFile">Reference fasta file.</param>
        /// <param name="binSize">Desired number of alignments per bin.</param>
        /// <param name="possibleAlignments">BitArrays of possible alignments.</param>
        /// <param name="observedAlignments">BitArrays of observed alignments.</param>
        /// <param name="predefinedBins">Pre-defined bins. null if not available.</param>
        /// <returns>A list of bins.</returns>
        static List<SampleGenomicBin> BinCounts(string referenceFile, int binSize, CanvasCoverageMode coverageMode, NexteraManifest manifest,
            Dictionary<string, BitArray> possibleAlignments,
            Dictionary<string, HitArray> observedAlignments,
            Dictionary<string, Int16[]> fragmentLengths,
            Dictionary<string, List<SampleGenomicBin>> predefinedBins,
            string outFile)
        {
            bool debugGCCorrection = false; // write value of GC bins and correction factor
            Dictionary<string, GenericRead> fastaEntries = new Dictionary<string, GenericRead>();
            List<string> chromosomes = new List<string>();
            Int16 meanFragmentSize = 0;
            Int16 meanFragmentCutoff = 3;
            if (coverageMode == CanvasCoverageMode.GCContentWeighted)
            {
                meanFragmentSize = MeanFragmentSize(fragmentLengths);
                if (meanFragmentSize <= 0)
                {
                    throw new Exception("CNV input observedAlignments error - Unable to determine fragment size, likely because of very low genome coverage.");
                }
            }

            using (FastaReader reader = new FastaReader(referenceFile))
            {
                GenericRead fastaEntry = new GenericRead();

                // Loop through each chromosome in the reference.
                while (reader.GetNextEntry(ref fastaEntry))
                {
                    chromosomes.Add(fastaEntry.Name);
                    fastaEntries[fastaEntry.Name] = fastaEntry;
                    fastaEntry = new GenericRead();
                }
            }

            // calculate GC content of the forward read at every position along the genome  
            Dictionary<string, byte[]> readGCContent = new Dictionary<string, byte[]>();
            if (coverageMode == CanvasCoverageMode.GCContentWeighted)
            {
                byte gcCap = (byte)numberOfGCbins;
                var normalizationTasks = new List<Action>();
                foreach (KeyValuePair<string, Int16[]> fragmentLengthsKVP in fragmentLengths)
                {
                    string chr = fragmentLengthsKVP.Key;
                    GenericRead fastaEntry = fastaEntries[chr];

                    normalizationTasks.Add(() =>
                    {
                        // contains GC content of the forward read at every position for current chr
                        byte[] gcContent = new byte[fastaEntry.Bases.Length];

                        uint gcCounter = 0;

                        // Iteratively calculate GC content of "reads" using fasta genome reference
                        for (int pos = 0; pos < fastaEntry.Bases.Length - meanFragmentSize * meanFragmentCutoff - 1; pos++)
                        {
                            Int16 currentFragment = 0;

                            if (fragmentLengthsKVP.Value[pos] == 0)
                                currentFragment = meanFragmentSize;
                            else
                                currentFragment = Convert.ToInt16(Math.Min(fragmentLengthsKVP.Value[pos], meanFragmentSize * meanFragmentCutoff));
                            for (int i = pos; i < pos + currentFragment; i++)
                            {
                                switch (fastaEntry.Bases[i])
                                {
                                    case 'C':
                                    case 'c':
                                    case 'G':
                                    case 'g':
                                        gcCounter++;
                                        break;
                                    default:
                                        break;
                                }
                            }
                            gcContent[pos] = (byte)Math.Min(100 * gcCounter / currentFragment, gcCap);
                            gcCounter = 0;
                        }
                        lock (readGCContent)
                        {
                            readGCContent[chr] = gcContent;
                        }
                    });
                }

                Console.WriteLine("{0} Launching normalization tasks.", DateTime.Now);
                Console.Out.Flush();
                Parallel.ForEach(normalizationTasks, task => task.Invoke());
                Console.WriteLine("{0} Normalization tasks complete.", DateTime.Now);
                Console.Out.Flush();
            }

            // populate observed and expected read GC bin vectors
            float[] observedVsExpectedGC = new float[0];
            if (coverageMode == CanvasCoverageMode.GCContentWeighted)
                observedVsExpectedGC = ComputeObservedVsExpectedGC(observedAlignments, readGCContent, manifest, debugGCCorrection, outFile);

            Dictionary<string, List<SampleGenomicBin>> perChromosomeBins = new Dictionary<string, List<SampleGenomicBin>>();
            var binningTasks = new List<Action>();
            foreach (KeyValuePair<string, GenericRead> fastaEntryKVP in fastaEntries)
            {
                string chr = fastaEntryKVP.Key;
                if (!possibleAlignments.ContainsKey(chr)) continue;
                if (predefinedBins != null && !predefinedBins.ContainsKey(chr)) continue;

                BinTaskArguments args = new BinTaskArguments();
                args.FastaEntry = fastaEntryKVP.Value;
                args.Chromosome = chr;
                args.PossibleAlignments = possibleAlignments[chr];
                args.ObservedAlignments = observedAlignments[chr];
                args.CoverageMode = coverageMode;
                perChromosomeBins[chr] = predefinedBins == null ? new List<SampleGenomicBin>() : predefinedBins[chr];
                args.Bins = perChromosomeBins[chr];
                args.BinSize = binSize;
                if (coverageMode == CanvasCoverageMode.GCContentWeighted)
                    args.ReadGCContent = readGCContent[chr];
                else
                    args.ReadGCContent = null;
                args.ObservedVsExpectedGC = observedVsExpectedGC;
                binningTasks.Add(() => { BinCountsForChromosome(args); });
            }
            Console.WriteLine("{0} Launch BinCountsForChromosome jobs...", DateTime.Now);
            Console.Out.WriteLine();
            Parallel.ForEach(binningTasks, task => task.Invoke());
            Console.WriteLine("{0} Completed BinCountsForChromosome jobs.", DateTime.Now);
            Console.Out.WriteLine();

            List<SampleGenomicBin> finalBins = new List<SampleGenomicBin>();
            foreach (string chr in chromosomes)
            {
                if (!perChromosomeBins.ContainsKey(chr)) continue;
                finalBins.AddRange(perChromosomeBins[chr]);
            }
            return finalBins;
        }

        class BinTaskArguments
        {
            public GenericRead FastaEntry;
            public string Chromosome;
            public BitArray PossibleAlignments;
            public HitArray ObservedAlignments;
            public CanvasCommon.CanvasCoverageMode CoverageMode;
            public List<SampleGenomicBin> Bins;
            public int BinSize;
            public byte[] ReadGCContent;
            public float[] ObservedVsExpectedGC;
        }

        /// <summary>
        /// Populate the list of GenomicBin objects for this chromosome.  
        /// </summary>
        static void BinCountsForChromosome(BinTaskArguments arguments)
        {
            List<SampleGenomicBin> bins = arguments.Bins;
            bool usePredefinedBins = bins.Any();
            int predefinedBinIndex = 0;
            GenericRead fastaEntry = arguments.FastaEntry; //fastaEntryKVP.Value;
            BinState currentBin = new BinState();
            string chr = arguments.Chromosome;
            BitArray possibleAlignments = arguments.PossibleAlignments;
            HitArray observedAlignments = arguments.ObservedAlignments;
            CanvasCoverageMode coverageMode = arguments.CoverageMode;
            int pos = usePredefinedBins ? bins[predefinedBinIndex].Start : 0;

            // Skip past leading Ns
            while (fastaEntry.Bases[pos].Equals('n'))
                pos++;
            List<float> binPositions = new List<float>();
            List<int> binObservations = new List<int>();
            for (; pos < fastaEntry.Bases.Length; pos++)
            {
                // Sets the start of the bin
                if (currentBin.StartPosition == -1)
                    currentBin.StartPosition = pos;

                if (!fastaEntry.Bases[pos].Equals("n"))
                    currentBin.NucleotideCount++;

                switch (fastaEntry.Bases[pos])
                {
                    case 'C':
                    case 'c':
                    case 'G':
                    case 'g':
                        currentBin.GCCount++;
                        break;
                }

                if (possibleAlignments[pos])
                {
                    currentBin.PossibleCount++;
                    currentBin.ObservedCount += observedAlignments.Data[pos];
                    binObservations.Add(observedAlignments.Data[pos]);
                    if (coverageMode == CanvasCoverageMode.GCContentWeighted)
                        binPositions.Add(arguments.ObservedVsExpectedGC[arguments.ReadGCContent[pos]]);
                }

                // We've seen the desired number of possible alignment positions.
                if ((!usePredefinedBins && currentBin.PossibleCount == arguments.BinSize)
                    || (usePredefinedBins && pos == bins[predefinedBinIndex].Stop - 1))
                {
                    if (coverageMode == CanvasCoverageMode.TruncatedDynamicRange) // Truncated dynamic range
                    {
                        currentBin.ObservedCount = 0;
                        foreach (int Value in binObservations)
                        {
                            currentBin.ObservedCount += Math.Min(10, Value);
                        }
                    }
                    if (coverageMode == CanvasCoverageMode.GCContentWeighted) // read GC content weighted 
                    {
                        currentBin.ObservedCount = 0;
                        float tmpObservedCount = 0;
                        for (int i = 0; i < binObservations.Count; i++)
                        {
                            tmpObservedCount += Math.Min(10, (float)binObservations[i] / binPositions[i]);
                        }
                        currentBin.ObservedCount = (int)Math.Round(tmpObservedCount);

                    }

                    int gc = (int)(100f * currentBin.GCCount / currentBin.NucleotideCount); // Note: Use 100f to avoid int overflow for a very large bin

                    if (usePredefinedBins)
                    {
                        bins[predefinedBinIndex].GenomicBin.GC = gc;
                        bins[predefinedBinIndex].Count = currentBin.ObservedCount;
                        predefinedBinIndex++;
                        if (predefinedBinIndex >= bins.Count) { break; } // we have processed all the bins
                        pos = bins[predefinedBinIndex].Start - 1; // jump to right before the next predefined bin
                    }
                    else
                    {
                        // Note the pos + 1 to make the first three conform to bed specification
                        SampleGenomicBin bin = new SampleGenomicBin(chr, currentBin.StartPosition, pos + 1, gc, currentBin.ObservedCount);
                        bins.Add(bin);
                    }

                    // Reset all relevant variables
                    currentBin.Reset();
                    binObservations.Clear();
                    binPositions.Clear();
                }
            }
        }

        /// <summary>
        /// Remove possible alignment positions if they intersect a supplied bed file.
        /// </summary>
        /// <param name="filterFile">Name of bed file to use for filtering.</param>
        /// <param name="tags">BitArrays of possible alignment positions.</param>
        static void ExcludeTagsOverlappingFilterFile(string filterFile, IDictionary<string, BitArray> tags)
        {
            using (FileStream stream = new FileStream(filterFile, FileMode.Open, FileAccess.Read))
            using (StreamReader reader = new StreamReader(stream))
            {
                string row;

                while ((row = reader.ReadLine()) != null)
                {

                    string[] fields = row.Split('\t');

                    string chr = fields[0];
                    int start = Convert.ToInt32(fields[1]);
                    int stop = Convert.ToInt32(fields[2]);

                    if (!tags.ContainsKey(chr))
                        continue;

                    for (int i = start; i < stop; i++)
                        tags[chr][i] = false;
                }
            }

        }

        /// <summary>
        /// Remove any observed alignment if it wasn't 'possible'.
        /// </summary>
        /// <param name="observedAlignments">BitArrays of observed alignment positions.</param>
        /// <param name="possibleAlignments">BitArrays of possible alignment positions.</param>
        static void ScreenObservedTags(IDictionary<string, HitArray> observedAlignments, IDictionary<string, BitArray> possibleAlignments)
        {

            foreach (string chr in possibleAlignments.Keys)
            {
                if (!observedAlignments.ContainsKey(chr))
                    continue;
                HitArray observed = observedAlignments[chr];
                BitArray possible = possibleAlignments[chr];
                for (int i = 0; i < possible.Length; i++)
                {
                    if (!possible[i])
                    {
                        observed.Data[i] = 0;
                    }
                }
            }
        }

        /// <summary>
        /// Deserialize CanvasBin object in multiple threads 
        /// </summary>
        /// <param name="inputFile">inputFile with per-chromosome CanvasBin objects.</param>
        /// <param name="possibleAlignments">Stores which alignments are possible (perfect and unique).</param>
        /// <param name="observedAlignments">Stores observed alignments from a sample.</param>
        /// <param name="fragmentLengths">Stores fragment length in byte format.</param>
        public static void DeserializeCanvasData(IFileLocation inputFile, Dictionary<string, BitArray> possibleAlignments,
            Dictionary<string, HitArray> observedAlignments, Dictionary<string, Int16[]> fragmentLengths,
            Object semaphore, CanvasCoverageMode coverageMode)
        {
            IntermediateData data = null;
            using (FileStream stream = new FileStream(inputFile.FullName, FileMode.Open, FileAccess.Read, FileShare.Read))
            {
                Stopwatch watch = new Stopwatch();
                watch.Start();
                data = ProtoBuf.Serializer.Deserialize<IntermediateData>(stream);
                watch.Stop();
                Console.WriteLine("File: {0}", inputFile);
                Console.WriteLine("Time elapsed: {0}", watch.Elapsed);
            }
            Dictionary<string, BitArray> tempPossibleAlignments;
            Dictionary<string, HitArray> tempObservedAlignments;
            Dictionary<string, Int16[]> tempFragmentLengths;
            data.GetData(out tempPossibleAlignments, out tempObservedAlignments, out tempFragmentLengths);
            lock (semaphore)
            {
                foreach (KeyValuePair<string, BitArray> kvp in tempPossibleAlignments)
                {
                    possibleAlignments.Add(kvp.Key, kvp.Value);
                }
                foreach (KeyValuePair<string, HitArray> kvp in tempObservedAlignments)
                {
                    observedAlignments.Add(kvp.Key, kvp.Value);
                }
                if (coverageMode == CanvasCoverageMode.GCContentWeighted)
                {
                    foreach (KeyValuePair<string, Int16[]> kvp in tempFragmentLengths)
                    {
                        fragmentLengths.Add(kvp.Key, kvp.Value);
                    }
                }

            }
        }


        private static int BinOneGenomicInterval(CanvasBinParameters parameters,
            Dictionary<string, BitArray> possibleAlignments,
            Dictionary<string, HitArray> observedAlignments,
            Dictionary<string, Int16[]> fragmentLengths)
        {
            InitializeAlignmentArrays(parameters.referenceFile, parameters.chromosome, parameters.coverageMode, possibleAlignments, observedAlignments, fragmentLengths);
            Console.WriteLine("{0} Initialized alignment arrays", DateTime.Now);
            LoadObservedAlignmentsBAM(parameters.bamFile, parameters.isPairedEnd, parameters.chromosome, parameters.coverageMode, observedAlignments[parameters.chromosome], fragmentLengths[parameters.chromosome]);
            Console.WriteLine("{0} Loaded observed alignments", DateTime.Now);

            // Filter on BED file.
            if (parameters.filterFile != null)
                ExcludeTagsOverlappingFilterFile(parameters.filterFile, possibleAlignments);

            // Make sure we don't have an 'impossible' observed alignment.
            ScreenObservedTags(observedAlignments, possibleAlignments);

            Console.WriteLine("{0} Serialize intermediate observedAlignments", DateTime.Now);
            //output binary intermediate file
            IntermediateData data = new IntermediateData(possibleAlignments, observedAlignments, fragmentLengths, parameters.coverageMode);
            Directory.CreateDirectory(Path.GetDirectoryName(parameters.outFile));
            using (FileStream stream = new FileStream(parameters.outFile, FileMode.Create, FileAccess.Write))
            {
                ProtoBuf.Serializer.Serialize<IntermediateData>(stream, data);
            }
            Console.WriteLine("{0} Intermediate observedAlignments serialized", DateTime.Now);
            return 0;
        }

        public class IntermediateDataHolder
        {
            public string SampleName;
            public IReadOnlyList<IFileLocation> InputFiles;
            public Dictionary<string, BitArray> PossibleAlignments;
            public Dictionary<string, HitArray> ObservedAlignments;
            public Dictionary<string, short[]> FragmentLengths;

        }

        /// <summary>
        /// Read BAM file, calculate and serialize observedAlignments and possibleAlignments to IntermediateData file
        /// </summary>
        public static void CalculateSampleHits(CanvasBinParameters parameters)
        {
            Dictionary<string, BitArray> possibleAlignments = new Dictionary<string, BitArray>();

            // Will hold a bunch of HitArrays, one for each chromosome.
            // Each one's length corresponds to the length of the chromosome it represents.
            // A position will be marked with the number of times the mer starting at that position 
            // is observed in the SAM file.
            Dictionary<string, HitArray> observedAlignments = new Dictionary<string, HitArray>();

            // Will hold a bunch of byte arrays, one for each chromosome.
            // Each one's length corresponds to the length of the chromosome it represents.
            // A value at a given index will represents fragment length of the read starting at that index
            Dictionary<string, Int16[]> fragmentLengths = new Dictionary<string, Int16[]>();
            BinOneGenomicInterval(parameters, possibleAlignments, observedAlignments, fragmentLengths);
        }

        public static void CalculateSingleSampleBinSize(CanvasBinParameters parameters,
            Dictionary<string, BitArray> possibleAlignments,
            Dictionary<string, HitArray> observedAlignments,
            Dictionary<string, Int16[]> fragmentLengths, NexteraManifest manifest)
        {
            var inputFiles = new List<IFileLocation>(parameters.intermediatePaths);
            Console.WriteLine("Start deserialization:");
            Console.Out.Flush();
            LoadIntermediateData(parameters, inputFiles, possibleAlignments, observedAlignments, fragmentLengths);
            Console.WriteLine("{0} Deserialization complete", DateTime.Now);
            Console.Out.Flush();
            if (parameters.binSize == -1)
                // Turn the desired # of alignments per bin into the number of possible alignments expected per bin.
                parameters.binSize = CalculateNumberOfPossibleAlignmentsPerBin(parameters.countsPerBin, possibleAlignments, observedAlignments,
                    manifest: manifest);
            Console.WriteLine("{0} binSize", parameters.binSize);
        }

        public static int CalculateMultiSampleBinSize(Dictionary<IFileLocation, IReadOnlyList<IFileLocation>> intermediateDataPaths,
            int binSize, int countsPerBin, CanvasCoverageMode coverageMode)
        {
            List<IntermediateDataHolder> intermediateData = new List<IntermediateDataHolder>();
            foreach (var intermediateDataPath in intermediateDataPaths.Values)
            {
                IntermediateDataHolder singleSample = new IntermediateDataHolder();
                singleSample.InputFiles = intermediateDataPath;
                singleSample.PossibleAlignments = new Dictionary<string, BitArray>();
                singleSample.ObservedAlignments = new Dictionary<string, HitArray>();
                singleSample.FragmentLengths = new Dictionary<string, Int16[]>();
                intermediateData.Add(singleSample);
            }
            Console.WriteLine("Start deserialization:");
            Console.Out.Flush();
            LoadIntermediateMultiSampleData(coverageMode, intermediateData);
            Console.WriteLine("{0} Deserialization complete", DateTime.Now);
            Console.Out.Flush();
            // Turn the desired # of alignments per bin into the number of possible alignments expected per bin.
            if (binSize == -1)
                binSize = CalculateNumberOfPossibleMultiSampleAlignmentsPerBin(countsPerBin, intermediateData);
            Console.WriteLine("{0} binSize", binSize);
            return binSize;
        }

        public static void CalculateSingleSampleBins(CanvasBinParameters parameters,
            Dictionary<string, BitArray> possibleAlignments,
            Dictionary<string, HitArray> observedAlignments,
            Dictionary<string, Int16[]> fragmentLengths, NexteraManifest manifest)
        {
            Dictionary<string, List<SampleGenomicBin>> predefinedBins = null;
            if (parameters.predefinedBinsFile != null)
            {
                // Read predefined bins
                predefinedBins = Utilities.LoadBedFile(parameters.predefinedBinsFile);
            }
            List<SampleGenomicBin> bins = BinCounts(parameters.referenceFile, parameters.binSize, parameters.coverageMode, manifest,
            possibleAlignments, observedAlignments, fragmentLengths, predefinedBins, parameters.outFile);
            // Output!
            Console.WriteLine("{0} Output binned counts:", DateTime.Now);
            CanvasIO.WriteToTextFile(parameters.outFile, bins);
            Console.WriteLine("{0} Output complete", DateTime.Now);
            Console.Out.Flush();
        }

        public static void CalculateMultiSampleBins(CanvasBinParameters parameters,
            List<IntermediateDataHolder> intermediateData, NexteraManifest manifest)
        {
            Dictionary<string, List<SampleGenomicBin>> predefinedBins = null;
            if (parameters.predefinedBinsFile != null)
            {
                // Read predefined bins
                predefinedBins = Utilities.LoadBedFile(parameters.predefinedBinsFile);
            }
            foreach (IntermediateDataHolder singleSample in intermediateData)
            {
                List<SampleGenomicBin> bins = BinCounts(parameters.referenceFile, parameters.binSize,
                    parameters.coverageMode, manifest,
                    singleSample.PossibleAlignments, singleSample.ObservedAlignments, singleSample.FragmentLengths,
                    predefinedBins, parameters.outFile);
                // Output!
                string outFile = Path.Combine(Path.GetDirectoryName(parameters.outFile) + singleSample.SampleName + ".binned");
                Console.WriteLine("{0} Output binned counts:", DateTime.Now);
                CanvasIO.WriteToTextFile(outFile, bins);
                Console.WriteLine("{0} Output complete", DateTime.Now);
                Console.Out.Flush();
            }
        }

        /// <summary>
        /// Compute bin size for multiple samples, compute and write CanvasBin data to a file separately for each sample
        /// </summary>
        public static void RunSingleSample(CanvasBinParameters parameters)
        {
            // Will hold information on intermediate Canvas data
            Dictionary<string, BitArray> possibleAlignments = new Dictionary<string, BitArray>();
            Dictionary<string, HitArray> observedAlignments = new Dictionary<string, HitArray>();
            Dictionary<string, Int16[]> fragmentLengths = new Dictionary<string, Int16[]>();
            // Load manifest
            NexteraManifest manifest = parameters.manifestFile == null ? null : new NexteraManifest(parameters.manifestFile, null, Console.WriteLine);
            // Derive bins
            CalculateSingleSampleBinSize(parameters, possibleAlignments, observedAlignments, fragmentLengths, manifest);
            if (parameters.binSizeOnly)
            {
                // Write bin size to file
                System.IO.File.WriteAllText(parameters.outFile + ".binsize", "" + parameters.binSize);
                return;
            }
            CalculateSingleSampleBins(parameters, possibleAlignments, observedAlignments, fragmentLengths, manifest);
        }

        /// <summary>
        /// Compute bin size for multiple samples, compute and write CanvasBin data to a file separately for each sample
        /// </summary>
        public static void RunMultiSample(CanvasBinParameters parameters)
        {
            // Will hold information on intermediate Canvas data
            List<IntermediateDataHolder> intermediateData = new List<IntermediateDataHolder>();
            // Load manifest
            NexteraManifest manifest = parameters.manifestFile == null ? null : new NexteraManifest(parameters.manifestFile, null, Console.WriteLine);
            // Derive bins
            CalculateMultiSampleBins(parameters, intermediateData, manifest);
        }

        /// <summary>
        /// Implements the Canvas binning algorithm
        /// </summary>
        public static int Run(CanvasBinParameters parameters)
        {
            // CanvasBin BAM reading step
            if (parameters.intermediatePaths.Count == 0 && parameters.inJson == null)
                CalculateSampleHits(parameters);
            else if (parameters.intermediatePaths.Count != 0 && parameters.inJson != null)
                throw new ArgumentException("-i/--infile or -j/--injson are mutually exclusive arguments");
            // MultiSample Canvas bin merge
            else if (parameters.inJson != null)
                RunMultiSample(parameters);
            // Single sample Canvas bin merge
            else
                RunSingleSample(parameters);
            return 0;
        }

        private static void LoadIntermediateData(CanvasBinParameters parameters, IReadOnlyList<IFileLocation> readOnlyInputFiles, Dictionary<string, BitArray> possibleAlignments,
            Dictionary<string, HitArray> observedAlignments, Dictionary<string, short[]> fragmentLengths)
        {
            Object semaphore = new object(); // control access to possibleAlignments, observedAlignments, fragmentLengths
            // int processorCoreCount = Environment.ProcessorCount;
            int processorCoreCount = 1; // Limit # of deserialization threads to avoid (rare) protobuf issue.
            List<Thread> threads = new List<Thread>();
            var inputFiles = readOnlyInputFiles.ToList();
            while (threads.Count > 0 || inputFiles.Count > 0)
            {
                // Remove defunct threads:
                threads.RemoveAll(t => !t.IsAlive);
                if (threads.Count == processorCoreCount)
                {
                    Thread.Sleep(1000);
                    continue;
                }
                while (inputFiles.Count > 0 && threads.Count < processorCoreCount)
                {
                    var inputFile = inputFiles.First();
                    ThreadStart threadDelegate =
                        new ThreadStart(
                            () =>
                                DeserializeCanvasData(inputFile, possibleAlignments, observedAlignments, fragmentLengths,
                                    semaphore, parameters.coverageMode));
                    Thread newThread = new Thread(threadDelegate);
                    threads.Add(newThread);
                    newThread.Name = "CanvasBin " + inputFiles[0];
                    Console.WriteLine(newThread.Name);
                    newThread.Start();
                    inputFiles.RemoveAt(0);
                }
            }
        }

        private static void LoadIntermediateMultiSampleData(CanvasCoverageMode coverageMode, List<IntermediateDataHolder> intermediateData)
        {
            Object semaphore = new object(); // control access to possibleAlignments, observedAlignments, fragmentLengths
            // int processorCoreCount = Environment.ProcessorCount;
            int processorCoreCount = 1; // Limit # of deserialization threads to avoid (rare) protobuf issue.
            List<Thread> threads = new List<Thread>();
            Console.WriteLine("Start deserialization:");
            Console.Out.Flush();
            foreach (IntermediateDataHolder singleSample in intermediateData)
            {
                var inputFiles = singleSample.InputFiles.ToList();
                while (threads.Count > 0 || inputFiles.Count > 0)
                {
                    // Remove defunct threads:
                    threads.RemoveAll(t => !t.IsAlive);
                    if (threads.Count == processorCoreCount)
                    {
                        Thread.Sleep(1000);
                        continue;
                    }
                    while (inputFiles.Count > 0 && threads.Count < processorCoreCount)
                    {
                        var inputFile = inputFiles.First();
                        ThreadStart threadDelegate =
                            new ThreadStart(
                                () =>
                                    DeserializeCanvasData(inputFile, singleSample.PossibleAlignments, singleSample.ObservedAlignments,
                                    singleSample.FragmentLengths, semaphore, coverageMode));
                        Thread newThread = new Thread(threadDelegate);
                        threads.Add(newThread);
                        newThread.Name = "CanvasBin " + inputFiles[0];
                        Console.WriteLine(newThread.Name);
                        newThread.Start();
                        inputFiles.RemoveAt(0);
                    }
                }
            }
        }

        [ProtoContract]
        private class IntermediateData
        {
            #region Members
            [ProtoMember(1)]
            public Dictionary<string, byte[]> PossibleAlignments = new Dictionary<string, byte[]>();
            [ProtoMember(2)]
            public Dictionary<string, byte[]> ObservedAlignments = new Dictionary<string, byte[]>();
            [ProtoMember(3)]
            public Dictionary<string, int> BitsInLastBytePossibleAlignments = new Dictionary<string, int>();
            [ProtoMember(4)]
            public Dictionary<string, Int16[]> FragmentLengths = new Dictionary<string, Int16[]>();
            #endregion

            public IntermediateData(Dictionary<string, BitArray> possibleAlignments, Dictionary<string, HitArray> observedAlignments, Dictionary<string, Int16[]> fragmentLengths, CanvasCoverageMode coverageMode)
            {
                foreach (KeyValuePair<string, BitArray> kvp in possibleAlignments)
                {
                    int bitsInLastByte = kvp.Value.Length % 8;
                    byte[] bytes = new byte[kvp.Value.Length / 8 + (bitsInLastByte == 0 ? 0 : 1)];
#if DotNetCore
                    BitArray array = kvp.Value;
                    int byteIndex = 0;
                    int bytePosition = 0;
                    for (int index = 0; index < array.Length; index++)
                    {
                        bytes[byteIndex] *= 2;
                        if (array[index]) bytes[byteIndex]++;
                        bytePosition = (bytePosition + 1) % 8;
                        if (bytePosition == 0) byteIndex++;
                    }
#else
                    kvp.Value.CopyTo(bytes, 0);
#endif
                    this.PossibleAlignments[kvp.Key] = bytes;
                    BitsInLastBytePossibleAlignments[kvp.Key] = bitsInLastByte;
                }

                foreach (KeyValuePair<string, HitArray> kvp in observedAlignments)
                {
                    this.ObservedAlignments[kvp.Key] = kvp.Value.Data;
                }

                if (coverageMode == CanvasCoverageMode.GCContentWeighted)
                {
                    foreach (KeyValuePair<string, Int16[]> kvp in fragmentLengths)
                    {
                        this.FragmentLengths[kvp.Key] = kvp.Value;
                    }
                }

            }

            public IntermediateData() { }

            public void GetData(out Dictionary<string, BitArray> possibleAlignments, out Dictionary<string, HitArray> observedAlignments, out Dictionary<string, Int16[]> fragmentLengths)
            {
                possibleAlignments = Convert(this.PossibleAlignments, this.BitsInLastBytePossibleAlignments);
                observedAlignments = new Dictionary<string, HitArray>();
                foreach (string key in this.ObservedAlignments.Keys)
                {
                    observedAlignments[key] = new HitArray(this.ObservedAlignments[key]);
                }
                fragmentLengths = new Dictionary<string, Int16[]>();
                foreach (string key in this.FragmentLengths.Keys)
                {
                    fragmentLengths[key] = this.FragmentLengths[key];
                }
            }

            private static Dictionary<string, BitArray> Convert(Dictionary<string, byte[]> alignments,
                Dictionary<string, int> bitsInLastByteAlignments)
            {
                Dictionary<string, BitArray> tempAlignments = new Dictionary<string, BitArray>();
                foreach (KeyValuePair<string, byte[]> kvp in alignments)
                {
                    int numberBitsInLastByte = bitsInLastByteAlignments[kvp.Key];
                    BitArray bits = null;
                    if (numberBitsInLastByte > 0)
                    {
                        int numBitsFromFullBytes = 8 * (kvp.Value.Length - 1);
                        bits = new BitArray(numBitsFromFullBytes + numberBitsInLastByte);
                        BitArrayCopyFrom(bits, kvp.Value, numBitsFromFullBytes);
                        byte lastByte = kvp.Value[kvp.Value.Length - 1];
                        for (int bitIndexLastByte = 0; bitIndexLastByte < numberBitsInLastByte; bitIndexLastByte++)
                        {
                            bool bit = (lastByte & (1 << bitIndexLastByte)) != 0;
                            bits.Set(numBitsFromFullBytes + bitIndexLastByte, bit);
                        }
                    }
                    else
                    {
                        bits = new BitArray(kvp.Value);
                    }

                    tempAlignments[kvp.Key] = bits;
                }
                return tempAlignments;
            }

            public static void BitArrayCopyFrom(BitArray bits, byte[] bytes, int numBitsToCopy)
            {
                BitArray tempBits = new BitArray(bytes);
                if (numBitsToCopy > tempBits.Length)
                    throw new ArgumentException("Source byte array does not contain enough bits for copying");
                if (numBitsToCopy > bits.Length)
                    throw new ArgumentException("Destination BitArray must be at least as large as number of bits to copy");
                for (int bitIndex = 0; bitIndex < numBitsToCopy; bitIndex++)
                {
                    bits.Set(bitIndex, tempBits.Get(bitIndex));
                }
            }
        }
    }
}
