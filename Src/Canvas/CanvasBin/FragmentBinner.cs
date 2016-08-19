using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;

using CanvasCommon;
using Isas.SequencingFiles;
using Isas.Shared.DataTypes;
using Isas.Shared.Utilities.FileSystem;

namespace CanvasBin
{
    public class FragmentBinner
    {
        private CanvasBinParameters parameters;

        public FragmentBinner(CanvasBinParameters parameters)
        {
            this.parameters = parameters;
        }

        /// <summary>
        /// Performs fragment binning.
        /// </summary>
        /// <returns></returns>
        public int Bin()
        {
            if (parameters.predefinedBinsFile == null)
            {
                throw new ApplicationException("Predefined bins in BED is required for fragment binning.");
            }
            if (!parameters.isPairedEnd) // Janus-SRS-189
            {
                throw new ApplicationException("Paired-end reads are required for fragment binning.");
            }

            Dictionary<string, List<GenomicBin>> predefinedBins = Utilities.LoadBedFile(parameters.predefinedBinsFile, gcIndex: 3);
            List<string> chromosomes = GetChromosomesInBam(); // used to order chromosomes

            if (!Utilities.IsSubset(predefinedBins.Keys, chromosomes))
            {
                throw new ApplicationException(
                    String.Format("Not all chromosomes in {0} are found in {1}.", parameters.predefinedBinsFile, parameters.bamFile));
            }

            // Count fragments by chromosome
            List<ThreadStart> binningThreads = new List<ThreadStart>();
            List<BinTask> tasks = new List<BinTask>();
            foreach (string chrom in chromosomes)
            {
                if (!predefinedBins.ContainsKey(chrom)) { continue; }
                BinTask task = new BinTask(parameters.referenceFile, chrom, parameters.bamFile, predefinedBins[chrom]);
                tasks.Add(task);
                binningThreads.Add(new ThreadStart(() => { task.DoIt(); }));
            }

            Console.WriteLine("Launch fragment binning jobs...");
            Console.Out.WriteLine();
            Parallel.ForEach(binningThreads, t => { t.Invoke(); });
            Console.WriteLine("Completed fragment binning jobs.");
            Console.Out.WriteLine();

            long usableFragmentCount = tasks.Select(t => t.UsableFragmentCount).Sum();
            if (usableFragmentCount == 0)
            {
                throw new ApplicationException(String.Format("No passing-filter fragments overlapping bins are found in {0}", parameters.bamFile));
            }

            // Aggregate bins
            List<GenomicBin> finalBins = new List<GenomicBin>();
            foreach (string chrom in chromosomes)
            {
                if (!predefinedBins.ContainsKey(chrom)) { continue; }
                finalBins.AddRange(predefinedBins[chrom]);
            }

            // Output!
            CanvasIO.WriteToTextFile(parameters.outFile, finalBins);

            return 0;
        }

        /// <summary>
        /// Gets chromosomes in parameters.bamFile.
        /// </summary>
        /// <returns></returns>
        private List<string> GetChromosomesInBam()
        {
            using (BamReader reader = new BamReader(parameters.bamFile))
            {
                return reader.GetReferenceNames();
            }
        }

        public class BinTask
        {
            public string FastaFile { get; private set; }
            public string Chromosome { get; private set; }
            public Bam Bam { get; private set; }
            public List<GenomicBin> Bins { get; private set; }
            public long UsableFragmentCount { get { return usableFragmentCount; } } // Passing-filter fragments overlapping bins
            private long usableFragmentCount;

            /// <summary>
            /// 
            /// </summary>
            /// <param name="fastaEntry">FASTA entry for the chromosome</param>
            /// <param name="chrom">chromosome</param>
            /// <param name="bamFile">path to BAM</param>
            /// <param name="bins">predefined bins</param>
            public BinTask(string fastaFile, string chrom, string bamFile, List<GenomicBin> bins)
            {
                FastaFile = fastaFile;
                Chromosome = chrom;
                Bam = new Bam(new FileLocation(bamFile));
                Bins = bins;
            }

            /// <summary>
            /// Performs fragment binning for the chromosome.
            /// </summary>
            public void DoIt()
            {
                // Initialize bins by setting Count = 0
                InitializeBins();

                // Populate GC for each bin if not available
                if (!IsBinGCAvailable()) { PopulateBinGC(); }

                // Bin fragments
                binFragments();
            }

            /// <summary>
            /// Sets bin counts to 0.
            /// </summary>
            private void InitializeBins()
            {
                Console.WriteLine("Initializing bin counts to 0 for {0}...", Chromosome);
                foreach (GenomicBin bin in Bins) { bin.Count = 0; }
            }

            /// <summary>
            /// Is %GC available for all bins?
            /// </summary>
            /// <returns></returns>
            private bool IsBinGCAvailable()
            {
                foreach (GenomicBin bin in Bins)
                {
                    if (bin.GC < 0)
                    {
                        return false;
                    }
                }

                return true;
            }

            /// <summary>
            /// Calculates %GC for bins on the chromosome.
            /// </summary>
            private void PopulateBinGC()
            {
                Console.WriteLine("Calculating %GC for each bin on {0}...", Chromosome);
                string referenceBases = FastaLoader.LoadFastaSequence(FastaFile, Chromosome);
                foreach (GenomicBin bin in Bins)
                {
                    double ntCount = 0;
                    double gcCount = 0;
                    for (int pos = bin.Start; pos < bin.Stop; pos++)
                    {
                        if (referenceBases[pos].Equals('n')) { continue; }
                        ntCount++;
                        if (Utilities.IsGC(referenceBases[pos])) { gcCount++; }
                    }
                    int gc = ntCount > 0 ? (int)(100 * gcCount / ntCount) : 0;
                    bin.GC = gc;
                }
            }

            /// <summary>
            /// Bins fragments.
            /// </summary>
            private void binFragments()
            {
                // Sanity check: The BAM index file must exist, in order for us to seek to our target chromosome!
                if (!Bam.Index.Exists)
                {
                    throw new Exception(string.Format("Fatal error: Bam index not found at {0}", Bam.Index.FullName));
                }

                long pairedAlignmentCount = 0; // keep track of paired alignments
                usableFragmentCount = 0;
                using (BamReader reader = new BamReader(Bam.BamFile.FullName))
                {
                    int desiredRefIndex = -1;
                    desiredRefIndex = reader.GetReferenceIndex(Chromosome);
                    if (desiredRefIndex == -1)
                    {
                        throw new ApplicationException(
                            string.Format("Unable to retrieve the reference sequence index for {0} in {1}.", Chromosome, Bam.BamFile.FullName));
                    }
                    bool result = reader.Jump(desiredRefIndex, 0);
                    if (!result)
                    {
                        // Note: This is not necessarily an error, it just means that there *are* no reads for this chromosome in this 
                        // .bam file.  That is not uncommon e.g. for truseq amplicon.
                        return;
                    }

                    Dictionary<string, int> readNameToBinIndex = new Dictionary<string, int>();
                    HashSet<string> samePositionReadNames = new HashSet<string>();
                    int binIndexStart = 0;
                    int prevPosition = -1;
                    BamAlignment alignment = new BamAlignment();
                    while (reader.GetNextAlignment(ref alignment, true))
                    {
                        int refID = alignment.RefID;

                        // quit if the current reference index is different from the desired reference index
                        if (refID != desiredRefIndex)
                            break;

                        if (refID == -1)
                            continue;

                        if (alignment.Position < prevPosition) // Make sure the BAM is properly sorted
                        {
                            throw new ApplicationException(
                                string.Format("The alignment on {0} are not properly sorted in {1}: {2}", Chromosome, Bam.BamFile.FullName, alignment.Name));
                        }
                        prevPosition = alignment.Position;

                        if (alignment.IsPaired()) { pairedAlignmentCount++; }

                        BinOneAlignment(alignment, FragmentBinnerConstants.MappingQualityThreshold, readNameToBinIndex,
                            samePositionReadNames, ref usableFragmentCount, Bins, ref binIndexStart);
                    }
                }
                if (pairedAlignmentCount == 0)
                {
                    throw new ApplicationException(string.Format("No paired alignments found for {0} in {1}", Chromosome, Bam.BamFile.FullName));
                }
            }

            /// <summary>
            /// Bins the fragment identified by alignment. Increases bin count if the first read of a pair passes all the filters.
            /// Decreases bin count if the second read of a pair does not pass all the filters.
            /// </summary>
            /// <param name="alignment"></param>
            /// <param name="qualityThreshold">minimum mapping quality</param>
            /// <param name="readNameToBinIndex">Dictionary of read name to bin index</param>
            /// <param name="usableFragmentCount">number of usable fragments</param>
            /// <param name="bins">predefined bins</param>
            /// <param name="binIndexStart">bin index from which to start searching for the best bin</param>
            public static void BinOneAlignment(BamAlignment alignment, uint qualityThreshold, Dictionary<string, int> readNameToBinIndex,
                HashSet<string> samePositionReadNames, ref long usableFragmentCount, List<GenomicBin> bins, ref int binIndexStart)
            {
                if (!alignment.IsMapped()) { return; }
                if (!alignment.IsMateMapped()) { return; }
                if (!alignment.IsPrimaryAlignment()) { return; }
                if (!(alignment.IsPaired() && alignment.IsProperPair())) { return; }

                bool duplicateFailedQCLowQuality = IsDuplicateFailedQCLowQuality(alignment, qualityThreshold);

                // Check whether we have binned the fragment using the mate
                if (readNameToBinIndex.ContainsKey(alignment.Name))
                {
                    // Undo binning when one of the reads is a duplicate, fails QC or has low mapping quality
                    if (duplicateFailedQCLowQuality)
                    {
                        usableFragmentCount--;
                        bins[readNameToBinIndex[alignment.Name]].Count--;
                    }
                    readNameToBinIndex.Remove(alignment.Name); // clean up
                    return;
                }
                if (duplicateFailedQCLowQuality) { return; }

                if (alignment.RefID != alignment.MateRefID) { return; } // does this ever happen?

                if (IsRightMostInPair(alignment)) { return; } // look at only one read of the pair
                // handle the case where alignment.Position == alignment.MatePosition
                if (alignment.Position == alignment.MatePosition)
                {
                    if (samePositionReadNames.Contains(alignment.Name))
                    {
                        samePositionReadNames.Remove(alignment.Name);
                        return;
                    }
                    samePositionReadNames.Add(alignment.Name);
                }
                if (alignment.FragmentLength == 0) { return; } // Janus-SRS-190: 0 when the information is unavailable

                // Try to bin the fragment
                int fragmentStart = alignment.Position; // 0-based, inclusive
                int fragmentStop = alignment.Position + alignment.FragmentLength; // 0-based, exclusive
                while (binIndexStart < bins.Count && bins[binIndexStart].Stop <= fragmentStart) // Bins[binIndexStart] on the left of the fragment
                {
                    binIndexStart++;
                }
                if (binIndexStart >= bins.Count) { return; } // all the remaining fragments are on the right of the last bin

                // now Bins[binIndexStart].Stop > fragmentStart
                int bestBinIndex = FindBestBin(bins, binIndexStart, fragmentStart, fragmentStop);
                if (bestBinIndex >= 0) // Bin the fragment
                {
                    usableFragmentCount++;
                    bins[bestBinIndex].Count++;
                    readNameToBinIndex[alignment.Name] = bestBinIndex;
                }
            }

            /// <summary>
            /// Checks if any of the conditions is true:
            /// 1. The read is a duplicate,
            /// 2. The read failed QC,
            /// 3. The read is of low mapping quality.
            /// </summary>
            /// <param name="alignment"></param>
            /// <returns></returns>
            public static bool IsDuplicateFailedQCLowQuality(BamAlignment alignment, uint qualityThreshold)
            {
                if (alignment.IsDuplicate()) { return true; }
                if (alignment.IsFailedQC()) { return true; }
                if (alignment.MapQuality == FragmentBinnerConstants.MappingQualityNotAvailable
                    || alignment.MapQuality < qualityThreshold)
                {
                    return true;
                }

                return false;
            }

            /// <summary>
            /// Is the read the right-most one (by genomic position) in a pair?
            /// </summary>
            /// <param name="alignment"></param>
            /// <returns></returns>
            public static bool IsRightMostInPair(BamAlignment alignment)
            {
                return alignment.Position > alignment.MatePosition;
            }

            /// <summary>
            /// Starting from the bin at binIndexStart, increments the bin index and searches for the bin
            /// that overlaps the fragment the most. In case of a tie, the bin encountered first is returned.
            /// </summary>
            /// <param name="binIndexStart"></param>
            /// <param name="fragmentStart"></param>
            /// <param name="fragmentStop"></param>
            /// <returns></returns>
            public static int FindBestBin(List<GenomicBin> bins, int binIndexStart, int fragmentStart, int fragmentStop)
            {
                int bestBinIndex = -1;
                int bestOverlap = 0;
                for (int binIndex = binIndexStart; binIndex < bins.Count; binIndex++)
                {
                    int overlapStart = Math.Max(bins[binIndex].Start, fragmentStart);
                    int overlapStop = Math.Min(bins[binIndex].Stop, fragmentStop);
                    int overlap = overlapStop - overlapStart;
                    if (overlap <= 0) { break; }
                    if (overlap > bestOverlap)
                    {
                        bestOverlap = overlap;
                        bestBinIndex = binIndex;
                    }
                }

                return bestBinIndex;
            }
        }
    }
}
