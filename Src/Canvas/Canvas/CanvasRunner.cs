using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using System.Text;
using System.Threading.Tasks;
using Canvas;
using CanvasBin;
using CanvasPartition;
using CanvasCommon;
using Illumina.Common;
using Isas.SequencingFiles;
using Isas.Shared.Checkpointing;
using Isas.Shared.DataTypes;
using Isas.Shared.Utilities;
using Isas.Shared.Utilities.FileSystem;
using Utilities = Isas.Shared.Utilities.Utilities;

namespace Illumina.SecondaryAnalysis
{
    /// <summary>
    /// Run Canvas tools to generate CNV calls:
    /// </summary>
    public class CanvasRunner
    {
        #region Members

        private readonly string _canvasFolder;
        private readonly CanvasCoverageMode _coverageMode = CanvasCoverageMode.TruncatedDynamicRange;
        private readonly int _countsPerBin;
        private readonly ILogger _logger;
        private readonly IWorkManager _workManager;
        private readonly ICheckpointRunnerAsync _checkpointRunner;
        private readonly bool _isSomatic;
        private readonly Dictionary<string, string> _customParameters = new Dictionary<string, string>();
        #endregion

        public CanvasRunner(ILogger logger, IWorkManager workManager, ICheckpointRunnerAsync checkpointRunner, bool isSomatic, CanvasCoverageMode coverageMode,
            int countsPerBin, Dictionary<string, string> customParameters = null)
        {
            _logger = logger;
            _workManager = workManager;
            _checkpointRunner = checkpointRunner;
            _isSomatic = isSomatic;
            _canvasFolder = Path.Combine(Utilities.GetAssemblyFolder(typeof(CanvasRunner)));
            _coverageMode = coverageMode;
            _countsPerBin = countsPerBin;
            if (customParameters != null)
            {
                _customParameters = new Dictionary<string, string>(customParameters, StringComparer.InvariantCultureIgnoreCase);
            }
        }

        private string SmallestFile(List<string> paths)
        {
            long minFileSize = long.MaxValue;
            string smallestBamPath = null;
            foreach (string path in paths)
            {
                long fileSize = (new FileInfo(path)).Length;
                if (smallestBamPath == null || minFileSize > fileSize)
                {
                    smallestBamPath = path;
                    minFileSize = fileSize;
                }
            }
            return smallestBamPath;
        }

        private int GetBinSize(CanvasCallset callset, string bamPath, List<string> intermediateDataPaths,
            string canvasReferencePath, string canvasBedPath)
        {
            string canvasBinPath = Path.Combine(_canvasFolder, "CanvasBin.exe");
            string executablePath = canvasBinPath;
            if (CrossPlatform.IsThisMono())
                executablePath = Utilities.GetMonoPath();

            StringBuilder commandLine = new StringBuilder();
            if (CrossPlatform.IsThisMono())
            {
                commandLine.AppendFormat("{0} ", canvasBinPath);
            }
            commandLine.AppendFormat("-b \"{0}\" ", bamPath);
            commandLine.AppendFormat("-p "); // Paired-end input mode (Isaac or BWA output)
            commandLine.AppendFormat("-r \"{0}\" ", canvasReferencePath);

            foreach (string path in intermediateDataPaths)
            {
                commandLine.AppendFormat("-i \"{0}\" ", path);
            }

            commandLine.AppendFormat("-y "); // bin size only

            if (callset.IsEnrichment) // manifest
            {
                if (!File.Exists(callset.TempManifestPath)) { NexteraManifestUtils.WriteNexteraManifests(callset.Manifest, callset.TempManifestPath); }
                commandLine.AppendFormat("-t \"{0}\" ", callset.TempManifestPath);
            }

            string outputStub = Path.Combine(Path.GetDirectoryName(callset.BinSizePath), Path.GetFileNameWithoutExtension(callset.BinSizePath));
            commandLine.AppendFormat("-f \"{0}\" -d {1} -o \"{2}\"", canvasBedPath, _countsPerBin, outputStub);

            UnitOfWork binJob = new UnitOfWork()
            {
                ExecutablePath = executablePath,
                LoggingFolder = _workManager.LoggingFolder.FullName,
                LoggingStub = Path.GetFileNameWithoutExtension(callset.BinSizePath),
                CommandLine = commandLine.ToString()
            };
            if (_customParameters.ContainsKey("CanvasBin"))
            {
                binJob.CommandLine = Utilities.MergeCommandLineOptions(binJob.CommandLine, _customParameters["CanvasBin"], true);
            }
            _workManager.DoWorkSingleThread(binJob);

            int binSize;
            using (StreamReader reader = new StreamReader(callset.BinSizePath))
            {
                binSize = int.Parse(reader.ReadLine());
            }

            return binSize;
        }

        /// <summary>
        /// Invoke CanvasBin.  Return null if this fails and we need to abort CNV calling for this sample.
        /// </summary>
        protected IFileLocation InvokeCanvasBin(CanvasCallset callset, string canvasReferencePath, string canvasBedPath, string ploidyBedPath)
        {
            StringBuilder commandLine = new StringBuilder();
            string canvasBinPath = Path.Combine(_canvasFolder, "CanvasBin.exe");
            string executablePath = canvasBinPath;
            if (CrossPlatform.IsThisMono())
                executablePath = Utilities.GetMonoPath();

            //use bam as input
            if (callset.Bam == null)
            {
                Console.WriteLine("Input bam file not seen for sample {0}_{1} - no CNV calls", callset.SampleName, callset.Id);
                return null;
            }
            List<string> bamPaths = new List<string>();
            bamPaths.Add(callset.Bam.BamFile.FullName);
            if (!(callset.IsEnrichment && callset.Manifest.CanvasControlAvailable)) // do not add normal BAMs if Canvas Control is available
            {
                bamPaths.AddRange(callset.NormalBamPaths.Select(bam => bam.BamFile.FullName));
            }

            // enrichment options
            if (callset.IsEnrichment) // manifest
            {
                if (!File.Exists(callset.TempManifestPath))
                {
                    NexteraManifestUtils.WriteNexteraManifests(callset.Manifest, callset.TempManifestPath);
                }
            }
          
            // read bams 
            var intermediateDataPathsByBamPath = ReadBams(callset.GenomeMetadata, callset.Bam.IsPairedEnd, new List<string>(){callset.Id}, callset.TempFolder,
                canvasReferencePath, canvasBedPath, bamPaths, commandLine, canvasBinPath, executablePath, callset.TempManifestPath);

            // get bin size (of the smallest BAM) if normal BAMs are given
            var intermediateDataPathsByBamPathCopy = (from x in intermediateDataPathsByBamPath
                                                      select x).ToDictionary(x => x.Key, x => x.Value.Select(y => y).ToList()); // deep dictionary copy
            int binSize = -1;
            if (bamPaths.Count > 1)
            {
                string smallestBamPath = SmallestFile(bamPaths);
                binSize = GetBinSize(callset, smallestBamPath, intermediateDataPathsByBamPathCopy[smallestBamPath],
                    canvasReferencePath, canvasBedPath);
            }
            else if (callset.IsEnrichment && callset.Manifest.CanvasControlAvailable)
            {
                binSize = callset.Manifest.CanvasBinSize.Value;
            }

            // derive Canvas bins
            var bamToBinned = BamToBinned(callset.TempFolder, callset.Bam.IsPairedEnd, new List<string>() {callset.Id}, canvasReferencePath, canvasBedPath, bamPaths, commandLine, canvasBinPath, binSize, intermediateDataPathsByBamPath, executablePath);

            string tumorBinnedPath = bamToBinned[callset.Bam.BamFile.FullName]; // binned tumor sample
            string outputPath = tumorBinnedPath;
            if (callset.NormalBamPaths.Any() || (callset.IsEnrichment && callset.Manifest.CanvasControlAvailable))
            {
                outputPath = InvokeCanvasNormalize(callset, tumorBinnedPath, bamToBinned, ploidyBedPath);
            }
            return new FileLocation(outputPath);
        }

        /// <summary>
        /// Invoke CanvasBin.  Return null if this fails and we need to abort CNV calling for this sample.
        /// </summary>
        protected List<string> InvokeCanvasBin(SmallPedigreeCallset callset, string canvasReferencePath, string canvasBedPath)
        {
            StringBuilder commandLine = new StringBuilder();
            string canvasBinPath = Path.Combine(_canvasFolder, "CanvasBin.exe");
            string executablePath = canvasBinPath;
            if (CrossPlatform.IsThisMono())
                executablePath = Utilities.GetMonoPath();

            //use bam as input
            List<string> bamPaths = new List<string>();
            foreach (Bam bam in callset.BamPaths)
            {
                if (bam.BamFile.FullName == null || !File.Exists(bam.BamFile.FullName))
                {
                    Console.WriteLine("Input bam file not seen {0}", bam.BamFile.FullName);
                    return null;
                }
                bamPaths.Add(bam.BamFile.FullName);
            }

            // read bams 
            var intermediateDataPathsByBamPath = ReadBams(callset.GenomeMetadata, true, callset.SampleNames, callset.TempFolder,
                canvasReferencePath, canvasBedPath, bamPaths, commandLine, canvasBinPath, executablePath, null);

            // get bin size (of the smallest BAM) if normal BAMs are given
            int binSize = -1;
            var intermediateDataPathsByBamPathCopy = (from x in intermediateDataPathsByBamPath
                                                      select x).ToDictionary(x => x.Key, x => x.Value.Select(y=>y).ToList()); // deep dictionary copy

            if (bamPaths.Count > 1)
            {
                string smallestBamPath = SmallestFile(bamPaths);
                CanvasBin.CanvasBin canvasBin = new CanvasBin.CanvasBin();
                binSize = canvasBin.CalculateMultiSampleBinSize(intermediateDataPathsByBamPathCopy, 
                    binSize, _countsPerBin, CanvasCommon.CanvasCoverageMode.TruncatedDynamicRange);
            }

            // derive Canvas bins
            var bamToBinned = BamToBinned(callset.TempFolder, true, callset.SampleNames, canvasReferencePath, canvasBedPath, bamPaths, commandLine, canvasBinPath, binSize, intermediateDataPathsByBamPath, executablePath);
            return bamToBinned.Values.ToList();
        }

        private Dictionary<string, string> BamToBinned(string tempFolder, bool isPairedEnd, IEnumerable<string> Id, string canvasReferencePath, string canvasBedPath, List<string> bamPaths,
            StringBuilder commandLine, string canvasBinPath, int binSize, Dictionary<string, List<string>> intermediateDataPathsByBamPaths,
            string executablePath)
        {
            Dictionary<string, string> bamToBinned = new Dictionary<string, string>();
            List<UnitOfWork> finalBinJobs = new List<UnitOfWork>();
            int bamIdx = 0;
            foreach (List<string> intermediateDataPathsByBamPath in intermediateDataPathsByBamPaths.Values)
            { 
                string bamPath = bamPaths[bamIdx];
                // finish up CanvasBin step by merging intermediate data and finally binning                
                string binnedPath = Path.Combine(tempFolder, string.Format("{0}_{1}.binned", Id.ToList()[bamIdx], bamIdx));
                bamToBinned[bamPath] = binnedPath;
                commandLine.Clear();
                if (CrossPlatform.IsThisMono())
                {
                    commandLine.AppendFormat("{0} ", canvasBinPath);
                }
                commandLine.AppendFormat("-b \"{0}\" ", bamPath);
                if (isPairedEnd) commandLine.AppendFormat("-p ");

                commandLine.AppendFormat("-r \"{0}\" ", canvasReferencePath);
                commandLine.AppendFormat("-f \"{0}\" -d {1} -o \"{2}\" ", canvasBedPath, _countsPerBin, binnedPath);
                if (binSize != -1)
                {
                    commandLine.AppendFormat("-z \"{0}\" ", binSize);
                }

                foreach (string path in intermediateDataPathsByBamPath)
                {
                    commandLine.AppendFormat("-i \"{0}\" ", path);
                    Console.WriteLine("path: {0}", path);
                }

                commandLine.AppendFormat("-m {0} ", _coverageMode);

                UnitOfWork finalBinJob = new UnitOfWork()
                {
                    ExecutablePath = executablePath,
                    LoggingFolder = _workManager.LoggingFolder.FullName,
                    LoggingStub = Path.GetFileName(binnedPath),
                    CommandLine = commandLine.ToString()
                };
                if (_customParameters.ContainsKey("CanvasBin"))
                {
                    finalBinJob.CommandLine = Utilities.MergeCommandLineOptions(finalBinJob.CommandLine,
                        _customParameters["CanvasBin"], true);
                }
                finalBinJobs.Add(finalBinJob);
                bamIdx++;
            }
            _workManager.DoWorkParallel(finalBinJobs, new TaskResourceRequirements(8, 25));
                // CanvasBin itself is multi-threaded
            return bamToBinned;
        }

        private Dictionary<string, List<string>> ReadBams(GenomeMetadata genomeInfo, bool isPairedEnd, 
            IEnumerable<string> Id, string tempFolder, string canvasReferencePath, string canvasBedPath, List<string> bamPaths,
            StringBuilder commandLine, string canvasBinPath, string executablePath, string tempManifestPath = null)
        {
            GenomeMetadata genomeMetadata = genomeInfo;
            List<UnitOfWork> binJobs = new List<UnitOfWork>();

            Dictionary<string, List<string>> intermediateDataPathsByBamPath = new Dictionary<string, List<string>>();
            for (int bamIndex = 0; bamIndex < bamPaths.Count; bamIndex++)
            {
                intermediateDataPathsByBamPath[Id.ToList()[bamIndex]] = new List<string>();
                foreach (
                    GenomeMetadata.SequenceMetadata sequenceMetadata in
                        genomeMetadata.Sequences.OrderByDescending(sequence => sequence.Length))
                {
                    // Only invoke CanvasBin for autosomes + allosomes;
                    // don't invoke it for mitochondrial chromosome or extra contigs or decoys
                    if (sequenceMetadata.Type != GenomeMetadata.SequenceType.Allosome && !sequenceMetadata.IsAutosome())
                        continue;

                    string bamPath = bamPaths[bamIndex];
                    commandLine.Clear();
                    if (CrossPlatform.IsThisMono())
                    {
                        commandLine.AppendFormat("{0} ", canvasBinPath);
                    }
                    commandLine.AppendFormat("-b \"{0}\" ", bamPath);
                    if (isPairedEnd) commandLine.AppendFormat("-p ");
                    commandLine.AppendFormat("-r \"{0}\" ", canvasReferencePath);
                    commandLine.AppendFormat("-c {0} ", sequenceMetadata.Name);
                    commandLine.AppendFormat("-m {0} ", _coverageMode);

                    string intermediateDataPath = Path.Combine(tempFolder, string.Format("{0}_{1}_{2}.dat",
                        Id.ToList()[bamIndex], bamIndex, sequenceMetadata.Name));
                    intermediateDataPathsByBamPath[Id.ToList()[bamIndex]].Add(intermediateDataPath);
                    commandLine.AppendFormat("-f \"{0}\" -d {1} -o \"{2}\" ", canvasBedPath, _countsPerBin, intermediateDataPath);
                    if (tempManifestPath != null)
                        commandLine.AppendFormat("-t \"{0}\" ", tempManifestPath);

                    UnitOfWork binJob = new UnitOfWork()
                    {
                        ExecutablePath = executablePath,
                        LoggingFolder = _workManager.LoggingFolder.FullName,
                        LoggingStub = Path.GetFileName(intermediateDataPath),
                        CommandLine = commandLine.ToString()
                    };
                    if (_customParameters.ContainsKey("CanvasBin"))
                    {
                        binJob.CommandLine = Utilities.MergeCommandLineOptions(binJob.CommandLine,
                            _customParameters["CanvasBin"], true);
                    }
                    binJobs.Add(binJob);
                }
            }
            _workManager.DoWorkParallelThreads(binJobs);
            return intermediateDataPathsByBamPath;
        }


        /// <summary>
        /// Invoke CanvasNormalize.
        /// </summary
        /// <param name="callset"></param>
        /// <returns>path to the bin ratio bed file</returns>
        protected string InvokeCanvasNormalize(CanvasCallset callset, string tumorBinnedPath, Dictionary<string, string> bamToBinned,
            string ploidyBedPath, string mode = "weightedaverage")
        {
            string ratioBinnedPath = Path.Combine(callset.TempFolder, string.Format("{0}.ratio.binned", callset.Id));

            string canvasNormalizePath = Path.Combine(_canvasFolder, "CanvasNormalize.exe");
            string executablePath = canvasNormalizePath;
            if (CrossPlatform.IsThisMono())
                executablePath = Utilities.GetMonoPath();

            StringBuilder commandLine = new StringBuilder();
            if (CrossPlatform.IsThisMono())
            {
                commandLine.AppendFormat("{0} ", canvasNormalizePath);
            }

            commandLine.AppendFormat("-t {0} ", tumorBinnedPath.WrapWithShellQuote()); // tumor bed

            if (callset.IsEnrichment && callset.Manifest.CanvasControlAvailable)
            {
                commandLine.AppendFormat("-n {0} ", callset.Manifest.CanvasControlBinnedPath.WrapWithShellQuote()); // normal bed
            }
            else
            {
                foreach (string normalBinnedPath in callset.NormalBamPaths.Select(path => bamToBinned[path.BamFile.FullName]))
                {
                    commandLine.AppendFormat("-n {0} ", normalBinnedPath.WrapWithShellQuote()); // normal bed
                }
            }

            commandLine.AppendFormat("-w {0} ", callset.NormalBinnedPath.WrapWithShellQuote()); // weighted average normal bed

            commandLine.AppendFormat("-o {0} ", ratioBinnedPath.WrapWithShellQuote()); // ratio bed

            if (callset.IsEnrichment) // manifest
            {
                if (!File.Exists(callset.TempManifestPath)) { NexteraManifestUtils.WriteNexteraManifests(callset.Manifest, callset.TempManifestPath); }
                commandLine.AppendFormat("-f {0} ", callset.TempManifestPath.WrapWithShellQuote());
            }

            commandLine.AppendFormat("-m {0} ", mode.WrapWithShellQuote());

            if (!string.IsNullOrEmpty(ploidyBedPath))
            {
                commandLine.AppendFormat("-p {0} ", ploidyBedPath.WrapWithShellQuote());
            }

            UnitOfWork normalizeJob = new UnitOfWork()
            {
                ExecutablePath = executablePath,
                LoggingFolder = _workManager.LoggingFolder.FullName,
                LoggingStub = Path.GetFileName(ratioBinnedPath),
                CommandLine = commandLine.ToString()
            };
            if (_customParameters.ContainsKey("CanvasNormalize"))
            {
                normalizeJob.CommandLine = Utilities.MergeCommandLineOptions(normalizeJob.CommandLine, _customParameters["CanvasNormalize"], true);
            }
            _workManager.DoWorkSingleThread(normalizeJob);

            return ratioBinnedPath;
        }

        /// <summary>
        /// Intersect bins with the targeted regions defined in callset.Manifest.
        /// Assumes that the targeted regions don't intersect, the bins are sorted by genomic location and the bins don't intersect.
        /// </summary>
        /// <param name="callset"></param>
        /// <param name="partitionedPath">Output of CanvasPartition. Bins are assumed to be sorted</param>
        /// <returns></returns>
        private IFileLocation IntersectBinsWithTargetedRegions(CanvasCallset callset, IFileLocation partitionedPath)
        {
            if (!partitionedPath.Exists) { return partitionedPath; }
            var rawPartitionedPath = partitionedPath.AppendName(".raw");
            if (rawPartitionedPath.Exists) { rawPartitionedPath.Delete(); }
            partitionedPath.MoveTo(rawPartitionedPath);

            //callset.Manifest
            Dictionary<string, List<NexteraManifest.ManifestRegion>> manifestRegionsByChrom = callset.Manifest.GetManifestRegionsByChromosome();

            // CanvasPartition output file is in the BED format
            //   start: 0-based, inclusive
            //   end: 0-based, exclusive
            // Manifest
            //   start: 1-based, inclusive
            //   end: 1-based, inclusive
            using (GzipReader reader = new GzipReader(rawPartitionedPath.FullName))
            using (GzipWriter writer = new GzipWriter(partitionedPath.FullName))
            {
                string currentChrom = null;
                int manifestRegionIdx = 0;
                string line;
                string[] toks;
                while ((line = reader.ReadLine()) != null)
                {
                    toks = line.Split('\t');
                    string chrom = toks[0];
                    int start = int.Parse(toks[1]) + 1; // 1-based, inclusive
                    int end = int.Parse(toks[2]); // 1-based, inclusive
                    if (chrom != currentChrom)
                    {
                        currentChrom = chrom;
                        manifestRegionIdx = 0;
                    }
                    if (!manifestRegionsByChrom.ContainsKey(currentChrom)) { continue; }
                    while (manifestRegionIdx < manifestRegionsByChrom[currentChrom].Count
                        && manifestRegionsByChrom[currentChrom][manifestRegionIdx].End < start) // |- manifest region -| |- bin -|
                    {
                        manifestRegionIdx++;
                    }
                    if (manifestRegionIdx >= manifestRegionsByChrom[currentChrom].Count || // |- last manifest region -| |- bin -|
                        end < manifestRegionsByChrom[currentChrom][manifestRegionIdx].Start) // |- bin -| |- manifest region -|
                    {
                        continue; // skip bin
                    }

                    // |- bin -|
                    //       |- manifest region -|
                    while (manifestRegionIdx < manifestRegionsByChrom[currentChrom].Count &&
                        end >= manifestRegionsByChrom[currentChrom][manifestRegionIdx].Start)
                    {
                        // calculate intersection
                        int intersectionStart = Math.Max(start, manifestRegionsByChrom[currentChrom][manifestRegionIdx].Start); // 1-based, inclusive
                        int intersectionEnd = Math.Min(end, manifestRegionsByChrom[currentChrom][manifestRegionIdx].End); // 1-based, inclusive
                                                                                                                          // start/end in BED format
                        toks[1] = String.Format("{0}", intersectionStart - 1); // 0-based, inclusive
                        toks[2] = String.Format("{0}", intersectionEnd); // 0-based, exclusive

                        // write intersected bin
                        writer.WriteLine(String.Join("\t", toks));

                        manifestRegionIdx++;
                    }
                }
            }

            return partitionedPath;
        }

        /// <summary>
        /// Invoke CanvasSNV on SmallPedigreeCallset callsets.  Return null if this fails and we need to abort CNV calling for this sample.
        /// </summary>
        protected void InvokeCanvasSnv(SmallPedigreeCallset callsets)
        {
            foreach (CanvasCallset callset in callsets.Callset)
                InvokeCanvasSnv(callset);
        }

        /// <summary>
        /// Invoke CanvasSNV.  Return null if this fails and we need to abort CNV calling for this sample.
        /// </summary>
        protected void InvokeCanvasSnv(CanvasCallset callset)
        {
            List<UnitOfWork> jobList = new List<UnitOfWork>();
            List<string> outputPaths = new List<string>();
            GenomeMetadata genomeMetadata = callset.GenomeMetadata;

            string tumorBamPath = callset.Bam.BamFile.FullName;
            string normalVcfPath = callset.NormalVcfPath.FullName;
            foreach (GenomeMetadata.SequenceMetadata chromosome in genomeMetadata.Sequences)
            {
                // Only invoke for autosomes + allosomes;
                // don't invoke it for mitochondrial chromosome or extra contigs or decoys
                if (chromosome.Type != GenomeMetadata.SequenceType.Allosome && !chromosome.IsAutosome())
                    continue;

                UnitOfWork job = new UnitOfWork();
                job.ExecutablePath = Path.Combine(_canvasFolder, "CanvasSNV.exe");
                if (CrossPlatform.IsThisMono())
                {
                    job.CommandLine = job.ExecutablePath;
                    job.ExecutablePath = Utilities.GetMonoPath();
                }

                string outputPath = Path.Combine(callset.TempFolder, string.Format("{0}-{1}.SNV.txt.gz", chromosome.Name, callset.Id));
                outputPaths.Add(outputPath);
                job.CommandLine += $" {chromosome.Name} {normalVcfPath} {tumorBamPath} {outputPath}";
                if (_customParameters.ContainsKey("CanvasSNV"))
                {
                    job.CommandLine = Utilities.MergeCommandLineOptions(job.CommandLine, _customParameters["CanvasSNV"], true);
                }
                job.LoggingFolder = _workManager.LoggingFolder.FullName;
                job.LoggingStub = string.Format("CanvasSNV-{0}-{1}", callset.Id, chromosome.Name);
                jobList.Add(job);
            }
            Console.WriteLine("Invoking {0} processor jobs...", jobList.Count);

            // Invoke CanvasSNV jobs:
            Console.WriteLine(">>>CanvasSNV start...");
            _workManager.DoWorkParallelThreads(jobList);
            Console.WriteLine(">>>CanvasSNV complete!");

            // Concatenate CanvasSNV results:
            ConcatenateCanvasSNVResults(callset.VfSummaryPath, outputPaths);
            ConcatenateCanvasSNVBafResults(callset.VfSummaryBafPath, outputPaths.Select(path => path + ".baf"));
        }

        protected void ConcatenateCanvasSNVResults(string vfSummaryPath, IEnumerable<string> outputPaths)
        {
            using (GzipWriter writer = new GzipWriter(vfSummaryPath))
            {
                bool headerWritten = false;
                foreach (string outputPath in outputPaths)
                {
                    using (GzipReader reader = new GzipReader(outputPath))
                    {
                        while (true)
                        {
                            string fileLine = reader.ReadLine();
                            if (fileLine == null) break;
                            if (fileLine.Length > 0 && fileLine[0] == '#')
                            {
                                if (headerWritten) continue;
                                headerWritten = true;
                            }
                            writer.WriteLine(fileLine);
                        }
                    }
                }
            }
        }

        protected void ConcatenateCanvasSNVBafResults(string vfSummaryBafPath, IEnumerable<string> outputBafPaths)
        {
            using (StreamWriter writer = new StreamWriter(vfSummaryBafPath))
            {
                string headers = null;
                foreach (string outputPath in outputBafPaths)
                {
                    using (StreamReader reader = new StreamReader(outputPath))
                    {
                        string fileLine = reader.ReadLine(); // headers
                        if (headers == null)
                        {
                            headers = fileLine;
                            writer.WriteLine(headers);
                        }
                        while (true)
                        {
                            fileLine = reader.ReadLine();
                            if (fileLine == null) break;
                            writer.WriteLine(fileLine);
                        }
                    }
                }
            }
        }

        public class CanvasCleanOutput
        {
            public IFileLocation CleanedPath { get; set; }
            public IFileLocation FfpePath { get; set; }

            public CanvasCleanOutput(IFileLocation cleanedPath, IFileLocation ffpePath)
            {
                CleanedPath = cleanedPath;
                FfpePath = ffpePath;
            }
        }

        public void CallSample(CanvasCallset callset)
        {
            Task.Run(() => CallSampleInternal(callset)).GetAwaiter().GetResult();
        }

        public void CallPedigree(SmallPedigreeCallset callset)
        {
            Task.Run(() => CallSampleInternal(callset)).GetAwaiter().GetResult();
        }


        /// <summary>
        /// Germline workflow:
        /// - Run CanvasBin, CanvasClean, CanvasPartition, CanvasDiploidCaller
        /// 
        /// Somatic workflow:
        /// - Run CanvasBin, CanvasClean, CanvasPartition, CanvasSNV, CanvasSomaticCaller
        /// </summary>
        private async Task CallSampleInternal(CanvasCallset callset)
        {
            Directory.CreateDirectory(callset.TempFolder);
            string canvasReferencePath = callset.KmerFasta.FullName;
            string canvasBedPath = callset.FilterBed.FullName;
            if (!File.Exists(canvasReferencePath))
            {
                throw new ApplicationException(string.Format("Error: Missing reference fasta file required for CNV calling at '{0}'", canvasReferencePath));
            }
            if (!File.Exists(canvasBedPath))
            {
                throw new ApplicationException(string.Format("Error: Missing filter bed file required for CNV calling at '{0}'", canvasBedPath));
            }

            // CanvasSNV
            var canvasSnvTask = _checkpointRunner.RunCheckpointAsync("CanvasSNV", () => InvokeCanvasSnv(callset));

            // Prepare ploidy file:
            string ploidyBedPath = callset.PloidyBed?.FullName;

            // CanvasBin:
            var binnedPath = _checkpointRunner.RunCheckpoint("CanvasBin", () => InvokeCanvasBin(callset, canvasReferencePath, canvasBedPath, ploidyBedPath));
            if (binnedPath == null) return;

            // CanvasClean:
            var canvasCleanOutput = _checkpointRunner.RunCheckpoint("CanvasClean", () => InvokeCanvasClean(callset, binnedPath));

            // CanvasPartition:
            var partitionedPath = _checkpointRunner.RunCheckpoint("CanvasPartition", () => InvokeCanvasPartition(callset, canvasCleanOutput.CleanedPath, canvasBedPath));

            // Intersect bins with manifest
            if (callset.IsEnrichment)
            {
                partitionedPath = _checkpointRunner.RunCheckpoint("Intersect bins with manifest",
                    () => IntersectBinsWithTargetedRegions(callset, partitionedPath));
            }

            await canvasSnvTask;

            // Variant calling
            _checkpointRunner.RunCheckpoint("Variant calling", () =>
            {
                if (_isSomatic)
                {
                    RunSomaticCalling(partitionedPath, callset, canvasBedPath, ploidyBedPath, canvasCleanOutput.FfpePath);
                }
                else
                {
                    RunGermlineCalling(partitionedPath, callset, ploidyBedPath);
                }
            });
        }


        private async Task CallSampleInternal(SmallPedigreeCallset callset)
        {
            Directory.CreateDirectory(callset.TempFolder);
            foreach (CanvasCallset singleSampleCallset in callset.Callset)
                Directory.CreateDirectory(singleSampleCallset.TempFolder);

            string canvasReferencePath = callset.KmerFasta.FullName;
            string canvasBedPath = callset.FilterBed.FullName;
            if (!File.Exists(canvasReferencePath))
            {
                throw new ApplicationException(string.Format("Error: Missing reference fasta file required for CNV calling at '{0}'", canvasReferencePath));
            }
            if (!File.Exists(canvasBedPath))
            {
                throw new ApplicationException(string.Format("Error: Missing filter bed file required for CNV calling at '{0}'", canvasBedPath));
            }

            // Prepare ploidy file:
            List<string> ploidyBedPaths = callset.Callset.Select(x=>x.PloidyBed?.FullName).ToList();

            // CanvasBin:
            var binnedPaths = _checkpointRunner.RunCheckpoint("CanvasBin", () => InvokeCanvasBin(callset, canvasReferencePath, canvasBedPath));
            if (binnedPaths == null) return;

            // CanvasClean:
            var canvasCleanOutput = _checkpointRunner.RunCheckpoint("CanvasClean", () => InvokeCanvasClean(callset, binnedPaths));

            // CanvasPartition:
            var partitionedPaths = _checkpointRunner.RunCheckpoint("CanvasPartition", () => InvokeCanvasPartition(callset, canvasCleanOutput, canvasBedPath));

            // CanvasSNV
            var canvasSnvTask = _checkpointRunner.RunCheckpointAsync("CanvasSNV", () => InvokeCanvasSnv(callset));

            // Variant calling
            await canvasSnvTask;
            RunGermlineCalling(partitionedPaths, callset, ploidyBedPaths);
        }

        private List<IFileLocation> WriteMergedCanvasPartition(List<IFileLocation> partitionedPaths, List<string> tempFolders, List<string> sampleNames)
        {
            Dictionary<string, List<GenomicBin>> multisamplePartitions = CanvasCommon.Utilities.LoadMultiSamplePartiotionedBedFile(partitionedPaths);
            List <IFileLocation> outPaths = new List<IFileLocation>();
            int count = 0;
            string mergedOutPath = Path.Combine(Directory.GetParent(tempFolders[0]).ToString(), "merged.partitioned");

            using (StreamWriter writer = new StreamWriter(mergedOutPath))
            {
                foreach (string chr in multisamplePartitions.Keys)
                {
                    foreach (GenomicBin genomicBin in multisamplePartitions[chr])
                    {
                        writer.Write(string.Format($"{genomicBin.Chromosome}\t{genomicBin.Start}\t{genomicBin.Stop}"));
                        for (int i = 0; i < sampleNames.Count; i++)
                            writer.Write(string.Format($"\t{genomicBin.CountBins.Count[i]}\t{genomicBin.CountBins.SegmentId[i]}"));
                        writer.Write("\n");
                    }
                }
            }

            foreach (string sampleName in sampleNames)
            {
                string outPath = Path.Combine(tempFolders[count], string.Format("{0}_merged.partitioned", sampleName));
                using (GzipWriter writer = new GzipWriter(outPath))
                {
                    foreach (string chr in multisamplePartitions.Keys)
                    {
                        foreach (GenomicBin genomicBin in multisamplePartitions[chr])
                        {
                            int segmentNum = genomicBin.CountBins.SegmentId.Max().Value;
                            writer.WriteLine(string.Format("{0}\t{1}\t{2}\t{3}\t{4}", genomicBin.Chromosome,
                                genomicBin.Start, genomicBin.Stop, genomicBin.CountBins.Count[count], segmentNum));
                        }
                    }
                }
                count++;
                outPaths.Add(new FileLocation(outPath));
            }
            return outPaths;
        }


        private List<IFileLocation> InvokeCanvasPartition(SmallPedigreeCallset callsets, List<IFileLocation> cleanedPaths, string canvasBedPath)
        {
            List<IFileLocation> partitionedPaths = new List<IFileLocation>();
            if (callsets.Callset.Count != cleanedPaths.Count)
                throw new Exception($"Number of output CanvasClean files {cleanedPaths.Count} is not equal to the number of Canvas callsets {callsets.Callset.Count}");
            for (int i = 0; i < callsets.Callset.Count; i++)
                partitionedPaths.Add(InvokeCanvasPartition(callsets.Callset[i], cleanedPaths[i], canvasBedPath));
            List<string> sampleNames = callsets.Callset.Select(x => x.SampleName).ToList();
            List<string> tmpFolders = callsets.Callset.Select(x => x.TempFolder).ToList();
            return WriteMergedCanvasPartition(partitionedPaths, tmpFolders, sampleNames);
        }

        private IFileLocation InvokeCanvasPartition(CanvasCallset callset, IFileLocation cleanedPath, string canvasBedPath)
        {
            StringBuilder commandLine = new StringBuilder();
            string executablePath = Path.Combine(_canvasFolder, "CanvasPartition.exe");
            if (CrossPlatform.IsThisMono())
            {
                commandLine.AppendFormat("{0} ", executablePath);
                executablePath = Utilities.GetMonoPath();
            }
            commandLine.AppendFormat("-i \"{0}\" ", cleanedPath);
            commandLine.AppendFormat("-b \"{0}\" ", canvasBedPath);
            string partitionedPath = Path.Combine(callset.TempFolder, string.Format("{0}.partitioned", callset.Id));
            commandLine.AppendFormat("-o \"{0}\" ", partitionedPath);
            if (!_isSomatic)
                commandLine.AppendFormat(" -g");

            UnitOfWork partitionJob = new UnitOfWork()
            {
                ExecutablePath = executablePath,
                LoggingFolder = _workManager.LoggingFolder.FullName,
                LoggingStub = Path.GetFileName(partitionedPath),
                CommandLine = commandLine.ToString()
            };
            if (_customParameters.ContainsKey("CanvasPartition"))
            {
                partitionJob.CommandLine = Utilities.MergeCommandLineOptions(partitionJob.CommandLine, _customParameters["CanvasPartition"], true);
            }
            _workManager.DoWorkSingleThread(partitionJob);
            return new FileLocation(partitionedPath);
        }

        /// <summary>
        /// Invoke CanvasClean on SmallPedigreeCallset callsets. 
        /// </summary>
        protected List<IFileLocation> InvokeCanvasClean(SmallPedigreeCallset callsets, List<string> binnedPaths)
        {
            List<IFileLocation> cleanedPaths = new List<IFileLocation>();
            if (callsets.Callset.Count != binnedPaths.Count)
                throw new Exception($"Number of output CanvasBin files {binnedPaths.Count} is not equal to the number of Canvas callsets {callsets.Callset.Count}");
            for (int i = 0; i < callsets.Callset.Count; i++)
            {
                IFileLocation binnedPath = new FileLocation(binnedPaths[i]);
                cleanedPaths.Add(InvokeCanvasClean(callsets.Callset[i], binnedPath).CleanedPath);
            }
            return cleanedPaths;
        }

        private CanvasCleanOutput InvokeCanvasClean(CanvasCallset callset, IFileLocation binnedPath)
        {
            StringBuilder commandLine = new StringBuilder();
            commandLine.Length = 0;
            string executablePath = Path.Combine(_canvasFolder, "CanvasClean.exe");
            if (CrossPlatform.IsThisMono())
            {
                commandLine.AppendFormat("{0} ", executablePath);
                executablePath = Utilities.GetMonoPath();
            }
            commandLine.AppendFormat("-i \"{0}\" ", binnedPath);
            var tempFolder = new DirectoryLocation(callset.TempFolder);
            var cleanedPath = tempFolder.GetFileLocation($"{callset.Id}.cleaned");
            commandLine.AppendFormat("-o \"{0}\" ", cleanedPath);
            commandLine.AppendFormat("-g");

            IFileLocation ffpePath = null;

            // TruSight Cancer has 1,737 targeted regions. The cut-off 2000 is somewhat arbitrary.
            // TruSight One has 62,309 targeted regions.
            // Nextera Rapid Capture v1.1 has 411,513 targeted regions.
            if (!callset.IsEnrichment || callset.Manifest.Regions.Count > 2000)
            {
                ffpePath = tempFolder.GetFileLocation("FilterRegions.txt");
                commandLine.AppendFormat(" -s -r -f \"{0}\"", ffpePath);
            }
            if (callset.IsEnrichment) // manifest
            {
                if (!File.Exists(callset.TempManifestPath))
                {
                    NexteraManifestUtils.WriteNexteraManifests(callset.Manifest, callset.TempManifestPath);
                }
                commandLine.AppendFormat(" -t \"{0}\"", callset.TempManifestPath);
            }
            UnitOfWork cleanJob = new UnitOfWork()
            {
                ExecutablePath = executablePath,
                LoggingFolder = _workManager.LoggingFolder.FullName,
                LoggingStub = cleanedPath.Name,
                CommandLine = commandLine.ToString()
            };
            if (_customParameters.ContainsKey("CanvasClean"))
            {
                cleanJob.CommandLine = Utilities.MergeCommandLineOptions(cleanJob.CommandLine, _customParameters["CanvasClean"], true);
            }
            _workManager.DoWorkSingleThread(cleanJob);

            var canvasCleanOutput = new CanvasCleanOutput(cleanedPath, ffpePath);
            return canvasCleanOutput;
        }

        protected void RunSomaticCalling(IFileLocation partitionedPath, CanvasCallset callset, string canvasBedPath,
            string ploidyBedPath, IFileLocation ffpePath)
        {

            // get somatic SNV output:
            string somaticSnvPath = callset.SomaticVcfPath?.FullName;

            // Prepare and run CanvasSomaticCaller job:
            UnitOfWork callerJob = new UnitOfWork();
            var cnvVcfPath = callset.OutputVcfPath;
            callerJob.ExecutablePath = Path.Combine(this._canvasFolder, "CanvasSomaticCaller.exe");
            if (CrossPlatform.IsThisMono())
            {
                callerJob.CommandLine = callerJob.ExecutablePath;
                callerJob.ExecutablePath = Utilities.GetMonoPath();
            }
            callerJob.CommandLine += string.Format(" -v {0}", callset.VfSummaryPath);
            callerJob.CommandLine += string.Format(" -i {0}", partitionedPath);
            callerJob.CommandLine += string.Format(" -o {0}", cnvVcfPath);
            callerJob.CommandLine += string.Format(" -b {0}", canvasBedPath);
            if (!string.IsNullOrEmpty(ploidyBedPath))
                callerJob.CommandLine += string.Format(" -p {0}", ploidyBedPath);
            callerJob.CommandLine += string.Format(" -n {0}", callset.SampleName);
            if (callset.IsEnrichment)
                callerJob.CommandLine += " -e";
            if (callset.IsDbSnpVcf) // a dbSNP VCF file is used in place of the normal VCF file
                callerJob.CommandLine += " -d";
            // get localSD metric:
            if (ffpePath != null)
            {
                // Sanity-check: CanvasClean does not always write this file. 
                // If it's not present, just carry on:
                if (ffpePath.Exists)
                {
                    callerJob.CommandLine += string.Format(" -f \"{0}\"", ffpePath);
                }
                else
                {
                    _logger.Info("Note: SD file not found at '{0}'", ffpePath);
                }
            }

            if (!string.IsNullOrEmpty(somaticSnvPath))
                callerJob.CommandLine += string.Format(" -s {0}", somaticSnvPath);
            callerJob.CommandLine += string.Format(" -r \"{0}\" ", callset.WholeGenomeFastaFolder);
            if (_customParameters.ContainsKey("CanvasSomaticCaller"))
            {
                callerJob.CommandLine = Utilities.MergeCommandLineOptions(callerJob.CommandLine, _customParameters["CanvasSomaticCaller"], true);
            }
            callerJob.LoggingFolder = _workManager.LoggingFolder.FullName;
            callerJob.LoggingStub = string.Format("SomaticCNV-{0}", callset.Id);
            _workManager.DoWorkSingleThread(callerJob);
        }

        protected void RunGermlineCalling(List<IFileLocation> partitionedPaths, SmallPedigreeCallset callsets, List<string> ploidyBedPaths)
        {
            List<CanvasCleanOutput> cleanedPaths = new List<CanvasCleanOutput>();
            if (callsets.Callset.Count != partitionedPaths.Count)
                throw new Exception($"Number of output CanvasPartition files {partitionedPaths.Count} is not equal to the number of Canvas callsets {callsets.Callset.Count}");
            for (int i = 0; i < callsets.Callset.Count; i++)
            {
                IFileLocation partitionedPath = partitionedPaths[i];
                RunGermlineCalling(partitionedPath, callsets.Callset[i],  ploidyBedPaths[i]);
            }
        }

        protected void RunGermlineCalling(IFileLocation partitionedPath, CanvasCallset callset, string ploidyBedPath)
        {
            StringBuilder commandLine = new StringBuilder();
            ////////////////////////////////////////////////////////
            // CanvasDiploidCaller:
            commandLine.Length = 0;
            string executablePath = Path.Combine(_canvasFolder, "CanvasDiploidCaller.exe");
            if (CrossPlatform.IsThisMono())
            {
                commandLine.AppendFormat("{0} ", executablePath);
                executablePath = Utilities.GetMonoPath();
            }
            commandLine.AppendFormat("-i \"{0}\" ", partitionedPath);
            commandLine.AppendFormat("-v \"{0}\" ", callset.VfSummaryPath);
            var cnvVcfPath = callset.OutputVcfPath;
            commandLine.AppendFormat("-o \"{0}\" ", cnvVcfPath);
            commandLine.AppendFormat("-n \"{0}\" ", callset.SampleName);
            commandLine.AppendFormat("-r \"{0}\" ", callset.WholeGenomeFastaFolder);
            if (!string.IsNullOrEmpty(ploidyBedPath))
            {
                commandLine.AppendFormat("-p \"{0}\" ", ploidyBedPath);
            }
            if (callset.IsDbSnpVcf) // a dbSNP VCF file is used in place of the normal VCF file
                commandLine.AppendFormat("-d ");
            UnitOfWork callJob = new UnitOfWork()
            {
                ExecutablePath = executablePath,
                LoggingFolder = _workManager.LoggingFolder.FullName,
                LoggingStub = cnvVcfPath.Name,
                CommandLine = commandLine.ToString()
            };
            if (_customParameters.ContainsKey("CanvasDiploidCaller"))
            {
                callJob.CommandLine = Utilities.MergeCommandLineOptions(callJob.CommandLine, _customParameters["CanvasDiploidCaller"], true);
            }
            _workManager.DoWorkSingleThread(callJob);
        }
    }
}
