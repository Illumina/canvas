using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Canvas.SmallPedigree;
using Canvas.Visualization;
using CanvasCommon;
using CanvasCommon.Visualization;
using Illumina.Common;
using Illumina.Common.FileSystem;
using Isas.ClassicBioinfoTools.Tabix;
using Isas.Framework.Checkpointing;
using Isas.Framework.Checkpointing.Legacy;
using Isas.Framework.Logging;
using Isas.Framework.Settings;
using Isas.Framework.Utilities;
using Isas.Framework.WorkManagement;
using Isas.Framework.WorkManagement.CommandBuilding;
using Isas.Manifests.NexteraManifest;
using Isas.SequencingFiles;

namespace Canvas
{
    /// <summary>
    /// Run Canvas tools to generate CNV calls:
    /// </summary>
    public class CanvasRunner
    {
        #region Members
        private readonly string _canvasFolder;
        private readonly CanvasCoverageMode _coverageMode = CanvasCoverageMode.TruncatedDynamicRange;
        private readonly CanvasNormalizeMode _normalizeMode = CanvasNormalizeMode.WeightedAverage;
        private readonly int _countsPerBin;
        public ILogger Logger { get; }
        private readonly IWorkDoer _workDoer;
        private readonly ICheckpointRunner _checkpointRunner;
        private readonly bool _isSomatic;
        private readonly Dictionary<string, string> _customParameters = new Dictionary<string, string>();
        private readonly IFileLocation _runtimeExecutable;
        private readonly Func<string, ICommandFactory> _runtimeCommandPrefix; // Path to either mono or dotnet
        private readonly IBAlleleBedGraphWriter _bAlleleBedGraphWriter;
        #endregion

        public CanvasRunner(ILogger logger, IWorkDoer workDoer,
            ICheckpointRunner checkpointRunner, IFileLocation runtimeExecutable,
            Func<string, ICommandFactory> runtimeCommandPrefix, bool isSomatic, CanvasCoverageMode coverageMode,
            int countsPerBin, IBAlleleBedGraphWriter bAlleleBedGraphWriter,
            Dictionary<string, string> customParameters, string canvasFolder)
        {
            Logger = logger;
            _workDoer = workDoer;
            _checkpointRunner = checkpointRunner;
            _isSomatic = isSomatic;
            _canvasFolder = canvasFolder;
            _coverageMode = coverageMode;
            _countsPerBin = countsPerBin;
            _bAlleleBedGraphWriter = bAlleleBedGraphWriter;
            _runtimeExecutable = runtimeExecutable;
            _runtimeCommandPrefix = runtimeCommandPrefix;
            if (customParameters != null)
            {
                _customParameters = new Dictionary<string, string>(customParameters, StringComparer.OrdinalIgnoreCase);
                UpdateCoverageMode(ref _coverageMode);
                UpdateNormalizeMode(ref _normalizeMode);
            }
        }

        private void UpdateCoverageMode(ref CanvasCoverageMode mode)
        {
            if (_customParameters.ContainsKey("CanvasBin"))
            {
                string beforeFirstOption;
                var options = Canvas.CommandOptionsUtilities.GetCommandOptions(_customParameters["CanvasBin"], out beforeFirstOption);
                foreach (var option in options)
                {
                    if (option.Key != "-m" && option.Key != "--mode")
                        continue;
                    mode = CanvasCommon.Utilities.ParseCanvasCoverageMode(option.Value.TrimStart('=').Trim());
                }
                // remove mode from custom parameters
                _customParameters["CanvasBin"] = Canvas.CommandOptionsUtilities.MergeCommandLineOptions(_customParameters["CanvasBin"], "#m #mode");
            }
        }

        private void UpdateNormalizeMode(ref CanvasNormalizeMode mode)
        {
            if (_customParameters.ContainsKey("CanvasNormalize"))
            {
                string beforeFirstOption;
                var options = Canvas.CommandOptionsUtilities.GetCommandOptions(_customParameters["CanvasNormalize"], out beforeFirstOption);
                foreach (var option in options)
                {
                    if (option.Key != "-m" && option.Key != "--mode")
                        continue;
                    mode = CanvasCommon.Utilities.ParseCanvasNormalizeMode(option.Value.TrimStart('=').Trim());
                }
                // remove mode from custom parameters
                _customParameters["CanvasNormalize"] = Canvas.CommandOptionsUtilities.MergeCommandLineOptions(_customParameters["CanvasNormalize"], "#m #mode");
            }
        }

        private IFileLocation SmallestFile(List<IFileLocation> paths)
        {
            long minFileSize = long.MaxValue;
            IFileLocation smallestBamPath = null;
            foreach (var path in paths)
            {
                long fileSize = path.Length;
                if (smallestBamPath == null || minFileSize > fileSize)
                {
                    smallestBamPath = path;
                    minFileSize = fileSize;
                }
            }
            return smallestBamPath;
        }

        /// <summary>
        /// Get the executable path, and stub command-line, for a canvas executable.
        /// .NET core: GetExecutablePath("CanvasBin") -> /path/to/dotnet, /path/toCanvasBin.dll
        /// .NET 4.x windows: GetExecutablePath("CanvasBin") -> /path/to/Canvas.exe, ''
        /// .NET 4.x mono: GetExecutablePath("CanvasBin") -> /path/to/mono, /path/toCanvasBin.exe
        /// </summary>
        internal string GetExecutablePath(string canvasExecutableStub, StringBuilder commandLineBuilder)
        {
            commandLineBuilder.Append(Path.Combine(_canvasFolder, canvasExecutableStub, string.Format("{0}.dll", canvasExecutableStub)));
            commandLineBuilder.Append(" ");
            return _runtimeExecutable.FullName;
        }

        private int GetBinSize(CanvasCallset callset, IFileLocation bamPath, IReadOnlyList<IFileLocation> intermediateDataPaths,
            string canvasReferencePath, string canvasBedPath)
        {
            var command = _runtimeCommandPrefix("CanvasBin").GetNewCommandBuilder();

            command.AddArgument("-b", bamPath);
            command.AddArgument("-p"); // Paired-end input mode (Isaac or BWA output)
            command.AddArgument("-r", canvasReferencePath);

            foreach (var path in intermediateDataPaths)
            {
                command.AddArgument("-i", path);
            }

            command.AddArgument("-y"); // bin size only

            if (callset.IsEnrichment) // manifest
            {
                if (!File.Exists(callset.TempManifestPath)) { NexteraManifestUtils.WriteNexteraManifests(callset.Manifest, callset.TempManifestPath); }
                command.AddArgument("-t", callset.TempManifestPath);
            }

            IFileLocation binSizeFile = new FileLocation(callset.SingleSampleCallset.BinSizePath);

            var outputFileNameStub = Path.GetFileNameWithoutExtension(binSizeFile.FullName);
            var outputStub = binSizeFile.Directory.GetFileLocation(outputFileNameStub);
            command.AddArgument("-f", canvasBedPath);
            command.AddArgument("-d", _countsPerBin);
            command.AddArgument("-o", outputStub);

            var result = _workDoer.DoWork(
                WorkResourceRequest.CreateLinear(1, 5, 23, 1),
                (workResources, jobLauncher) => new FileLocation(callset.SingleSampleCallset.BinSizePath)).Await();

            int binSize;
            using (var stream = result.OpenRead())
            using (StreamReader reader = new StreamReader(stream))
            {
                binSize = int.Parse(reader.ReadLine());
            }

            return binSize;
        }

        /// <summary>
        /// Invoke CanvasBin.  Return null if this fails and we need to abort CNV calling for this sample.
        /// </summary>
        protected IFileLocation InvokeCanvasBin(CanvasCallset callset, string canvasReferencePath, string canvasBedPath, string ploidyVcfPath)
        {
            // use bam as input
            if (callset.SingleSampleCallset.Bam == null)
            {
                Console.WriteLine($"Input bam file not seen for sample {callset.SingleSampleCallset.SampleName}- no CNV calls");
                return null;
            }

            if (_coverageMode == CanvasCoverageMode.Fragment)
            {
                return InvokeCanvasBinFragment(callset, canvasReferencePath, canvasBedPath, ploidyVcfPath);
            }
            else
            {
                return InvokeCanvasBin35Mers(callset, canvasReferencePath, canvasBedPath, ploidyVcfPath);
            }
        }

        /// <summary>
        /// Invoke CanvasBin in a mode that is not the Fragment mode.  Return null if this fails and we need to abort CNV calling for this sample.
        /// </summary>
        protected IFileLocation InvokeCanvasBin35Mers(CanvasCallset callset, string canvasReferencePath, string canvasBedPath, string ploidyVcfPath)
        {
            StringBuilder commandLine = new StringBuilder();

            //use bam as input
            if (callset.SingleSampleCallset.Bam == null)
            {
                Console.WriteLine("Input bam file not seen for sample {0}_{1} - no CNV calls", callset.SingleSampleCallset.SampleName, callset.SingleSampleCallset.SampleName);
                return null;
            }
            var bamPaths = new List<IFileLocation> { callset.SingleSampleCallset.Bam.BamFile };
            if (!(callset.IsEnrichment && callset.Manifest.CanvasControlAvailable)) // do not add normal BAMs if Canvas Control is available
            {
                bamPaths.AddRange(callset.NormalBamPaths.Select(bam => bam.BamFile));
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
            var intermediateDataPathsByBamPath = GetIntermediateBinnedFilesByBamPath(
                callset.AnalysisDetails.GenomeMetadata, callset.SingleSampleCallset.Bam.IsPairedEnd, 
                new List<string>() { callset.SingleSampleCallset.SampleName }, callset.AnalysisDetails.TempDirectory,
                canvasReferencePath, canvasBedPath, bamPaths, callset.IsEnrichment ? callset.TempManifestPath : null);

            int binSize = -1;
            if (bamPaths.Count > 1)
            {
                var smallestBamPath = SmallestFile(bamPaths);
                binSize = GetBinSize(callset, smallestBamPath, intermediateDataPathsByBamPath[smallestBamPath],
                    canvasReferencePath, canvasBedPath);
            }
            else if (callset.IsEnrichment && callset.Manifest.CanvasControlAvailable)
            {
                if (callset.Manifest.CanvasBinSize != null) binSize = callset.Manifest.CanvasBinSize.Value;
            }

            // derive Canvas bins
            var bamToBinned = BamToBinned(callset.SingleSampleCallset.SampleOutputFolder, callset.SingleSampleCallset.Bam.IsPairedEnd, 
                new List<string>() { callset.SingleSampleCallset.SampleName }, canvasReferencePath, canvasBedPath, bamPaths, binSize, intermediateDataPathsByBamPath);

            var tumorBinnedPath = bamToBinned[callset.SingleSampleCallset.Bam.BamFile]; // binned tumor sample
            var outputPath = tumorBinnedPath;
            if (callset.NormalBamPaths.Any() || (callset.IsEnrichment && callset.Manifest.CanvasControlAvailable))
            {
                outputPath = InvokeCanvasNormalize(callset, tumorBinnedPath, bamToBinned, ploidyVcfPath);
            }
            return outputPath;
        }

        /// <summary>
        /// Invoke CanvasBin.  Return null if this fails and we need to abort CNV calling for this sample.
        /// </summary>
        protected List<IFileLocation> InvokeCanvasBin(SmallPedigreeCallset callset, string canvasReferencePath, string canvasBedPath)
        {
            //use bam as input
            var bamPaths = callset.PedigreeSample.Select(sample => sample.Sample.Bam.BamFile).ToList();

            var sampleNames = callset.PedigreeSample.Select(x => x.Sample.SampleName).ToList();
            // read bams 
            var intermediateDataPathsByBamPath = GetIntermediateBinnedFilesByBamPath(callset.AnalysisDetails.GenomeMetadata, true, sampleNames, callset.AnalysisDetails.TempDirectory,
                canvasReferencePath, canvasBedPath, bamPaths);

            int binSize = -1;
            if (bamPaths.Count > 1)
            {
                binSize = CanvasBin.CanvasBin.CalculateMultiSampleBinSize(intermediateDataPathsByBamPath,
                    binSize, _countsPerBin, CanvasCoverageMode.TruncatedDynamicRange);
            }

            // derive Canvas bins
            var bamToBinned = BamToBinned(callset.AnalysisDetails.TempDirectory, true, sampleNames, canvasReferencePath, canvasBedPath, bamPaths, binSize, intermediateDataPathsByBamPath);
            return bamToBinned.Values.ToList();
        }

        private Dictionary<IFileLocation, IFileLocation> BamToBinned(
            IDirectoryLocation tempFolder, bool isPairedEnd, List<string> sampleIds,
            string canvasReferencePath, string canvasBedPath, List<IFileLocation> bamPaths,
            int binSize, Dictionary<IFileLocation, IReadOnlyList<IFileLocation>> intermediateDataPathsByBam)
        {
            List<WorkToDo<(IFileLocation, IFileLocation)>> finalBinWork = Enumerable.Range(0, bamPaths.Count).SelectWork(
                WorkResourceRequest.CreateExact(8, 25),
                (bamIdx, resources, jobLauncher) =>
                {
                    var sampleId = sampleIds[bamIdx];
                    var bamPath = bamPaths[bamIdx];

                    var intermediateDataPaths = intermediateDataPathsByBam[bamPath];
                    // finish up CanvasBin step by merging intermediate data and finally binning 
                    var binnedPath = tempFolder.GetFileLocation($"{sampleId}_{bamIdx}.binned");
                    var commandLine = new StringBuilder();
                    string executablePath = GetExecutablePath("CanvasBin", commandLine);

                    commandLine.Append($"-b \"{bamPath}\" ");
                    if (isPairedEnd) commandLine.AppendFormat("-p ");

                    commandLine.Append($"-r \"{canvasReferencePath}\" ");
                    commandLine.Append($"-f \"{canvasBedPath}\" -d {_countsPerBin} -o \"{binnedPath}\" ");
                    if (binSize != -1)
                    {
                        commandLine.Append($"-z {binSize} ");
                    }

                    foreach (var path in intermediateDataPaths)
                    {
                        commandLine.Append($"-i \"{path}\" ");
                    }

                    commandLine.Append($"-m {_coverageMode} ");

                    var command = commandLine.ToString();
                    if (_customParameters.ContainsKey("CanvasBin"))
                    {
                        command = Canvas.CommandOptionsUtilities.MergeCommandLineOptions(
                            command, _customParameters["CanvasBin"], true);
                    }

                    var job = new JobInfo(executablePath, command, binnedPath.Name);
                    jobLauncher.LaunchJob(job);

                    return (bamPath, binnedPath);
                }).ToList();

            var pathTups = _workDoer.DoWork(finalBinWork).Await();

            return pathTups.ToDictionary();
        }

        private Dictionary<IFileLocation, IReadOnlyList<IFileLocation>> GetIntermediateBinnedFilesByBamPath(GenomeMetadata genomeMetadata, bool isPairedEnd,
            List<string> sampleIds, IDirectoryLocation tempFolder, string canvasReferencePath, string canvasBedPath, List<IFileLocation> bamPaths, string tempManifestPath = null)
        {
            List<WorkToDo<(IFileLocation, IFileLocation)>> binJobs = Enumerable
                .Range(0, bamPaths.Count)
                .SelectMany((bamIndex) =>
                {
                    var sampleId = sampleIds[bamIndex];
                    var bamPath = bamPaths[bamIndex];
                    var work = genomeMetadata
                        .Contigs()
                        .OrderByDescending(sequence => sequence.Length)
                        // Only invoke CanvasBin for autosomes + allosomes;
                        // don't invoke it for mitochondrial chromosome or extra contigs or decoys
                        .Where(contig => contig.Type == GenomeMetadata.SequenceType.Allosome || contig.IsAutosome())
                        .SelectWork(
                            WorkResourceRequest.CreateExact(1, 10),
                            (sequenceMetadata, _, jobLauncher) =>
                            {
                                var commandLine = new StringBuilder();
                                string executablePath = GetExecutablePath("CanvasBin", commandLine);

                                commandLine.AppendFormat("-b \"{0}\" ", bamPath);
                                if (isPairedEnd) commandLine.AppendFormat("-p ");
                                commandLine.AppendFormat("-r \"{0}\" ", canvasReferencePath);
                                commandLine.AppendFormat("-c {0} ", sequenceMetadata.Name);
                                commandLine.AppendFormat("-m {0} ", _coverageMode);

                                var intermediateDataPath = tempFolder.GetFileLocation($"{sampleId}_{bamIndex}_{sequenceMetadata.Name}.dat");
                                commandLine.AppendFormat("-f \"{0}\" -d {1} -o \"{2}\" ", canvasBedPath, _countsPerBin, intermediateDataPath);
                                if (tempManifestPath != null)
                                    commandLine.AppendFormat("-t \"{0}\" ", tempManifestPath);

                                var command = commandLine.ToString();
                                if (_customParameters.ContainsKey("CanvasBin"))
                                {
                                    command = Canvas.CommandOptionsUtilities.MergeCommandLineOptions(
                                        command, _customParameters["CanvasBin"], true);
                                }
                                var job = new JobInfo(executablePath, command, intermediateDataPath.Name);
                                jobLauncher.LaunchJob(job);
                                return (bamPath, intermediateDataPath);
                            });
                    return work;
                }).ToList();

            var intermediateDataPathsByBamPath = new Dictionary<IFileLocation, List<IFileLocation>>();

            var tuples = _workDoer.DoWork(binJobs).Await();
            foreach ((var bamPath, var dataPath) in tuples)
                if (intermediateDataPathsByBamPath.TryGetValue(bamPath, out var list))
                    list.Add(dataPath);
                else
                    (intermediateDataPathsByBamPath[bamPath] = new List<IFileLocation>()).Add(dataPath);

            return intermediateDataPathsByBamPath.ToDictionary(kvp => kvp.Key, kvp => (IReadOnlyList<IFileLocation>)kvp.Value);
        }

        /// <summary>
        /// Invoke CanvasBin in the Fragment mode.  Return null if this fails and we need to abort CNV calling for this sample.
        /// </summary>
        protected IFileLocation InvokeCanvasBinFragment(CanvasCallset callset, string canvasReferencePath, string canvasBedPath, string ploidyVcfPath)
        {
            // require predefined bins
            string predefinedBinsPath = GetPredefinedBinsPath();
            if (string.IsNullOrEmpty(predefinedBinsPath))
            {
                Console.WriteLine("Predefined bins are required to run CanvasBin in the Fragment mode");
                return null;
            }

            var bamPaths = new List<IFileLocation>();
            bool isPairedEnd = true;
            bamPaths.Add(callset.SingleSampleCallset.Bam.BamFile);
            isPairedEnd = isPairedEnd && callset.SingleSampleCallset.Bam.IsPairedEnd;
            if (!(callset.IsEnrichment && callset.Manifest.CanvasControlAvailable)) // do not add normal BAMs if Canvas Control is available
            {
                bamPaths.AddRange(callset.NormalBamPaths.Select(bam => bam.BamFile));
                isPairedEnd = isPairedEnd && callset.NormalBamPaths.All(bam => bam.IsPairedEnd);
            }

            // require paired-end reads
            if (!isPairedEnd)
            {
                Console.WriteLine("Paired-end reads are required to run CanvasBin in the Fragment mode");
                return null;
            }

            var work = Enumerable.Range(0, bamPaths.Count).SelectWork(
                WorkResourceRequest.CreateExact(8, 25), // CanvasBin itself is multi-threaded
                (bamIndex, resources, jobLauncher) => 
                {
                    var bamPath = bamPaths[bamIndex];
                    var binnedPath = callset.SingleSampleCallset.SampleOutputFolder.GetFileLocation($"{callset.SingleSampleCallset.SampleName}_{bamIndex}.binned");

                    var commandLine = new StringBuilder();
                    string executablePath = GetExecutablePath("CanvasBin", commandLine);

                    commandLine.AppendFormat(" -p -b \"{0}\" ", bamPath);
                    commandLine.AppendFormat("-r \"{0}\" ", canvasReferencePath);
                    commandLine.AppendFormat("-m {0} ", _coverageMode);
                    commandLine.AppendFormat("-f \"{0}\" -o \"{1}\" ", canvasBedPath, binnedPath);
                    commandLine.AppendFormat("-n {0} ", predefinedBinsPath); // assumes that predefinedBinsPath has been properly quoted

                    var command = commandLine.ToString();
                    if (_customParameters.ContainsKey("CanvasBin"))
                    {
                        command = Canvas.CommandOptionsUtilities.MergeCommandLineOptions(
                            command, _customParameters["CanvasBin"], true);
                    }
                    var job = new JobInfo(executablePath, command, binnedPath.Name);
                    jobLauncher.LaunchJob(job);

                    return (bamPath, binnedPath);
                }).ToList();

            var bamToBinned = _workDoer.DoWork(work).Await().ToDictionary(); 

            return NormalizeCoverage(callset, bamToBinned, ploidyVcfPath);
        }

        private string GetPredefinedBinsPath()
        {
            string path = null;
            if (_customParameters.ContainsKey("CanvasBin"))
            {
                string beforeFirstOption;
                var options = Canvas.CommandOptionsUtilities.GetCommandOptions(_customParameters["CanvasBin"], out beforeFirstOption);
                foreach (var option in options)
                {
                    if (option.Key != "-n" && option.Key != "--bins")
                        continue;
                    path = option.Value.TrimStart('=').Trim();
                }
                // remove bins from custom parameters
                _customParameters["CanvasBin"] = Canvas.CommandOptionsUtilities.MergeCommandLineOptions(_customParameters["CanvasBin"], "#n #bins");
            }
            return path;
        }

        protected IFileLocation NormalizeCoverage(CanvasCallset callset, Dictionary<IFileLocation, IFileLocation> bamToBinned, string ploidyVcfPath)
        {
            var tumorBinnedPath = bamToBinned[callset.SingleSampleCallset.Bam.BamFile]; // binned tumor sample
            var outputPath = tumorBinnedPath;
            if (callset.NormalBamPaths.Any() ||
                (callset.IsEnrichment && (callset.Manifest.CanvasControlAvailable)) ||
                _normalizeMode == CanvasNormalizeMode.PCA)
            {
                outputPath = InvokeCanvasNormalize(callset, tumorBinnedPath, bamToBinned, ploidyVcfPath);
            }

            return outputPath;
        }

        /// <summary>
        /// Invoke CanvasNormalize.
        /// </summary
        /// <param name="callset"></param>
        /// <returns>path to the bin ratio bed file</returns>
        protected IFileLocation InvokeCanvasNormalize(CanvasCallset callset, IFileLocation tumorBinnedPath, Dictionary<IFileLocation, IFileLocation> bamToBinned,
            string ploidyVcfPath)
        {
            var ratioBinnedPath = callset.SingleSampleCallset.SampleOutputFolder.GetFileLocation($"{callset.SingleSampleCallset.SampleName}.ratio.binned");
            StringBuilder commandLine = new StringBuilder();
            string executablePath = GetExecutablePath("CanvasNormalize", commandLine);

            commandLine.AppendFormat("-t \"{0}\" ", tumorBinnedPath); // tumor bed

            if ((callset.IsEnrichment && callset.Manifest.CanvasControlAvailable) ||
                _normalizeMode == CanvasNormalizeMode.PCA)
            {
                commandLine.AppendFormat("-n \"{0}\" ", callset.Manifest.CanvasControlBinnedPath); // normal bed
            }
            else
            {
                foreach (var normalBinnedPath in callset.NormalBamPaths.Select(path => bamToBinned[path.BamFile]))
                {
                    commandLine.AppendFormat("-n \"{0}\" ", normalBinnedPath); // normal bed
                }
            }

            commandLine.AppendFormat("-w \"{0}\" ", callset.SingleSampleCallset.NormalBinnedPath); // weighted average normal bed

            commandLine.AppendFormat("-o \"{0}\" ", ratioBinnedPath); // ratio bed

            if (callset.IsEnrichment) // manifest
            {
                if (!File.Exists(callset.TempManifestPath)) { NexteraManifestUtils.WriteNexteraManifests(callset.Manifest, callset.TempManifestPath); }
                commandLine.AppendFormat("-f \"{0}\" ", callset.TempManifestPath);
            }

            commandLine.AppendFormat("-m {0} ", _normalizeMode);

            if (!string.IsNullOrEmpty(ploidyVcfPath))
            {
                commandLine.AppendFormat("-p \"{0}\" ", ploidyVcfPath);
            }

            var command = commandLine.ToString();
            if (_customParameters.ContainsKey("CanvasNormalize"))
            {
                command = Canvas.CommandOptionsUtilities.MergeCommandLineOptions(
                    command, _customParameters["CanvasNormalize"], true);
            }

            var job = new JobInfo(executablePath, command, ratioBinnedPath.Name);
            return _workDoer.DoWork(WorkResourceRequest.CreateExact(1, 8), job, ratioBinnedPath).Await();
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
            foreach (PedigreeSample callset in callsets.PedigreeSample)
            {
                var canvasCallset = new CanvasCallset(callset.Sample, callsets.AnalysisDetails, null, null, null);
                InvokeCanvasSnv(canvasCallset, false, callset.Sample.SampleName);
            }
        }

        /// <summary>
        /// Invoke CanvasSNV.  Return null if this fails and we need to abort CNV calling for this sample.
        /// </summary>
        protected IFileLocation InvokeCanvasSnv(CanvasCallset callset, bool isSomatic = false, string sampleName = null)
        {
            GenomeMetadata genomeMetadata = callset.AnalysisDetails.GenomeMetadata;
            string bamPath = callset.SingleSampleCallset.Bam.BamFile.FullName;
            string normalVcfPath = callset.SingleSampleCallset.NormalVcfPath.FullName;

            var work = genomeMetadata.Contigs()
                .Where(chromosome => chromosome.Type == GenomeMetadata.SequenceType.Allosome || chromosome.IsAutosome())
                .SelectWork(
                    WorkResourceRequest.CreateExact(1, 10), 
                    (chromosome, resources, jobLauncher) => 
                    {
                        StringBuilder commandLine = new StringBuilder();
                        var executablePath = GetExecutablePath("CanvasSNV", commandLine);
                        IFileLocation outputPath = callset.SingleSampleCallset.SampleOutputFolder.GetFileLocation($"{chromosome.Name}-{callset.SingleSampleCallset.SampleName}.SNV.txt.gz");
                        commandLine.Append($" -c {chromosome.Name} -v {normalVcfPath} -b {bamPath} -o {outputPath}");
                        if (!sampleName.IsNullOrEmpty())
                            commandLine.Append($" -n {sampleName}");
                        if (callset.SingleSampleCallset.IsDbSnpVcf)
                            commandLine.Append(" --isDbSnpVcf");
                        if (isSomatic)
                            commandLine.Append(" --isSomatic");
                        var command = commandLine.ToString();
                        if (_customParameters.ContainsKey("CanvasSNV"))
                        {
                            command = Canvas.CommandOptionsUtilities.MergeCommandLineOptions(
                                command, _customParameters["CanvasSNV"], true);
                        }

                        var job = new JobInfo(executablePath, command, $"CanvasSNV-'{callset.SingleSampleCallset.SampleName}'-'{chromosome.Name}'");
                        jobLauncher.LaunchJob(job);
                        return outputPath;
                    }).ToList();

            Console.WriteLine($"Invoking {work.Count} processor jobs...for sample {callset.SingleSampleCallset.SampleName}");

            // Invoke CanvasSNV jobs:
            Console.WriteLine($"CanvasSNV start for sample {callset.SingleSampleCallset.SampleName}");
            var outputPaths = _workDoer.DoWork(work).Await();
            Console.WriteLine($"CanvasSNV complete for sample {callset.SingleSampleCallset.SampleName}");

            // Concatenate CanvasSNV results:
            ConcatenateCanvasSNVResults(callset.SingleSampleCallset.VfSummaryPath, outputPaths);
            ConcatenateCanvasSNVBafResults(callset.SingleSampleCallset.VfSummaryBafPath, outputPaths.Select(path => path + ".baf"));

            var bafFile = new FileLocation(callset.SingleSampleCallset.VfSummaryBafPath);
            var fileConversionMessage = $"converting '{bafFile}' to '{callset.SingleSampleCallset.BAlleleCoverageBedGraph.FileLocation}'";
            Logger.Info($"Begin {fileConversionMessage}");
            var benchmark = new Benchmark();
            _bAlleleBedGraphWriter.Write(bafFile, callset.SingleSampleCallset.BAlleleCoverageBedGraph);
            Logger.Info($"Finished {fileConversionMessage}. Elapsed time: {benchmark.GetElapsedTime()}");
            return new FileLocation(callset.SingleSampleCallset.VfSummaryPath);
        }

        protected void ConcatenateCanvasSNVResults(string vfSummaryPath, IEnumerable<IFileLocation> outputPaths)
        {
            using (GzipWriter writer = new GzipWriter(vfSummaryPath))
            {
                bool headerWritten = false;
                foreach (var outputPath in outputPaths)
                {
                    using (GzipReader reader = new GzipReader(outputPath.FullName))
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
            using (FileStream outStream = new FileStream(vfSummaryBafPath, FileMode.Create, FileAccess.Write))
            using (StreamWriter writer = new StreamWriter(outStream))
            {
                string headers = null;
                foreach (string outputPath in outputBafPaths)
                {
                    using (FileStream inStream = new FileStream(outputPath, FileMode.Open, FileAccess.Read))
                    using (StreamReader reader = new StreamReader(inStream))
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
            Logger.Info($"Normal Vcf path: {callset.SingleSampleCallset.NormalVcfPath}");
            Directory.CreateDirectory(callset.SingleSampleCallset.SampleOutputFolder.FullName);
            string canvasReferencePath = callset.AnalysisDetails.KmerFasta.FullName;
            string canvasBedPath = callset.AnalysisDetails.FilterBed.FullName;
            if (!File.Exists(canvasReferencePath))
            {
                throw new Illumina.Common.IlluminaException(string.Format("Error: Missing reference fasta file required for CNV calling at '{0}'", canvasReferencePath));
            }
            if (!File.Exists(canvasBedPath))
            {
                throw new Illumina.Common.IlluminaException(string.Format("Error: Missing filter bed file required for CNV calling at '{0}'", canvasBedPath));
            }

            // CanvasSNV
            var canvasSnvTask = _checkpointRunner.RunCheckpointAsync<IFileLocation>("CanvasSNV", () => _isSomatic ?
                InvokeCanvasSnv(callset, isSomatic: _isSomatic) : InvokeCanvasSnv(callset));

            // Prepare ploidy file:
            string ploidyVcfPath = callset.AnalysisDetails.PloidyVcf?.FullName;

            // CanvasBin:
            var binnedPath = _checkpointRunner.RunCheckpoint("CanvasBin", () => InvokeCanvasBin(callset, canvasReferencePath, canvasBedPath, ploidyVcfPath));
            if (binnedPath == null) return;

            // CanvasClean:
            var canvasCleanOutput = _checkpointRunner.RunCheckpoint("CanvasClean", () => InvokeCanvasClean(callset, binnedPath));

            var canvasSnvPath = await canvasSnvTask;
            // CanvasPartition:
            var partitionedPath = _checkpointRunner.RunCheckpoint("CanvasPartition", () => InvokeCanvasPartition(callset, canvasCleanOutput.CleanedPath, canvasBedPath, canvasSnvPath));

            // Intersect bins with manifest
            if (callset.IsEnrichment)
            {
                partitionedPath = _checkpointRunner.RunCheckpoint("Intersect bins with manifest",
                    () => IntersectBinsWithTargetedRegions(callset, partitionedPath));
            }

            // Variant calling
            _checkpointRunner.RunCheckpoint("Variant calling", () =>
            {
                if (_isSomatic)
                {
                    RunSomaticCalling(partitionedPath, callset, canvasBedPath, ploidyVcfPath, canvasCleanOutput.FfpePath, canvasSnvPath);
                }
                else
                {
                    RunGermlineCalling(partitionedPath, callset, ploidyVcfPath, canvasSnvPath);
                }
            });
        }


        private async Task CallSampleInternal(SmallPedigreeCallset callset)
        {
            callset.AnalysisDetails.TempDirectory.Create();
            foreach (var pedigreeSample in callset.PedigreeSample)
                Directory.CreateDirectory(pedigreeSample.Sample.SampleOutputFolder.FullName);

            string canvasReferencePath = callset.AnalysisDetails.KmerFasta.FullName;
            string canvasBedPath = callset.AnalysisDetails.FilterBed.FullName;
            string commonCnvsBed = null;
            if (callset.AnalysisDetails.CommonCnvsBed != null)
                commonCnvsBed = callset.AnalysisDetails.CommonCnvsBed.FullName;
            if (!File.Exists(canvasReferencePath))
            {
                throw new Illumina.Common.IlluminaException(
                    $"Error: Missing reference fasta file required for CNV calling at '{canvasReferencePath}'");
            }
            if (!File.Exists(canvasBedPath))
            {
                throw new Illumina.Common.IlluminaException(
                    $"Error: Missing filter bed file required for CNV calling at '{canvasBedPath}'");
            }

            // CanvasSNV
            var canvasSnvTask = _checkpointRunner.RunCheckpointAsync("CanvasSNV", () => InvokeCanvasSnv(callset));

            // CanvasBin:
            var binnedPaths = _checkpointRunner.RunCheckpoint("CanvasBin", () => InvokeCanvasBin(callset, canvasReferencePath, canvasBedPath));
            if (binnedPaths == null) return;

            // CanvasClean:
            var canvasCleanOutput = _checkpointRunner.RunCheckpoint("CanvasClean", () => InvokeCanvasClean(callset, binnedPaths));

            // CanvasPartition:
            var partitionedPaths = _checkpointRunner.RunCheckpoint("CanvasPartition", () => InvokeCanvasPartitionMultisample(callset, canvasCleanOutput, canvasBedPath, commonCnvsBed));

            // Variant calling
            await canvasSnvTask;
            _checkpointRunner.RunCheckpoint("Variant calling", () =>
            {
                RunSmallPedigreeCalling(partitionedPaths, callset);
            });
        }

        private void NormalizeCanvasClean(List<IFileLocation> cleanedPaths, string tempFolder)
        {
            Dictionary<string, List<MultiSampleGenomicBin>> normalizedCanvasClean = CanvasCommon.Utilities.MergeMultiSampleCleanedBedFile(cleanedPaths);
            int fileCounter = 0;
            foreach (IFileLocation cleanedPath in cleanedPaths)
            {
                using (GzipWriter writer = new GzipWriter(cleanedPath.FullName))
                {
                    foreach (string chr in normalizedCanvasClean.Keys)
                    {
                        foreach (MultiSampleGenomicBin genomicBin in normalizedCanvasClean[chr])
                        {
                            string outLine = string.Format($"{genomicBin.Bin.Chromosome}\t{genomicBin.Bin.Interval.Start}\t{genomicBin.Bin.Interval.End}");
                            outLine += string.Format($"\t{genomicBin.Counts[fileCounter]}");
                            writer.WriteLine(outLine);
                        }
                    }
                }
                fileCounter++;
            }
        }

        private List<IFileLocation> InvokeCanvasPartitionMultisample(SmallPedigreeCallset callsets, List<IFileLocation> cleanedPaths, string canvasBedPath, string commonCnvsBed)
        {
            NormalizeCanvasClean(cleanedPaths, callsets.AnalysisDetails.TempDirectory.FullName);
            StringBuilder commandLine = new StringBuilder();
            string executablePath = GetExecutablePath("CanvasPartition", commandLine);

            foreach (IFileLocation cleanedPath in cleanedPaths)
                commandLine.AppendFormat("-i \"{0}\" ", cleanedPath);

            commandLine.AppendFormat("-b \"{0}\" ", canvasBedPath);

            if (!commonCnvsBed.IsNullOrEmpty())
                commandLine.AppendFormat("-c \"{0}\" ", commonCnvsBed);

            List<IFileLocation> partitionedPaths = new List<IFileLocation>();
            foreach (var pedigreeSample in callsets.PedigreeSample)
            {
                IFileLocation partitionedPath = pedigreeSample.Sample.PartitionedPath;
                partitionedPaths.Add(partitionedPath);
                commandLine.AppendFormat("-o \"{0}\" ", partitionedPath);
            }
            commandLine.AppendFormat("-r \"{0}\" ", callsets.AnalysisDetails.WholeGenomeFastaFolder);
            commandLine.AppendFormat("-m HMM");

            var command = commandLine.ToString();
            if (_customParameters.ContainsKey("CanvasPartition"))
            {
                command = Canvas.CommandOptionsUtilities.MergeCommandLineOptions
                    (command, _customParameters["CanvasPartition"], true);
            }
            var job = new JobInfo(executablePath, command, partitionedPaths.First().Name);
            return _workDoer.DoWork(WorkResourceRequest.CreateExact(1,8), job, partitionedPaths).Await();
        }

        private IFileLocation InvokeCanvasPartition(CanvasCallset callset, IFileLocation cleanedPath, string canvasBedPath, IFileLocation canvasSnvPath)
        {
            StringBuilder commandLine = new StringBuilder();
            string executablePath = GetExecutablePath("CanvasPartition", commandLine);
            commandLine.Append($" -v {canvasSnvPath} ");
            commandLine.AppendFormat("-i \"{0}\" ", cleanedPath);
            commandLine.AppendFormat("-b \"{0}\" ", canvasBedPath);
            var partitionedPath = callset.SingleSampleCallset.PartitionedPath;
            commandLine.AppendFormat("-o \"{0}\" ", partitionedPath);
            commandLine.Append($" -r \"{callset.AnalysisDetails.WholeGenomeFastaFolder}\" ");
            if (!_isSomatic)
                commandLine.AppendFormat(" -g");
            else
            {
                if (!callset.IsEnrichment || callset.Manifest.Regions.Count > 2000)
                {
                    var tempFolder = new DirectoryLocation(callset.SingleSampleCallset.SampleOutputFolder.FullName);
                    var ffpePath = tempFolder.GetFileLocation("FilterRegions.txt");
                    commandLine.AppendFormat("-f \"{0}\" ", ffpePath);
                }
            }

            var command = commandLine.ToString();
            if (_customParameters.ContainsKey("CanvasPartition"))
            {
                command = Canvas.CommandOptionsUtilities.MergeCommandLineOptions(
                    command, _customParameters["CanvasPartition"], true);
            }
            var job = new JobInfo(executablePath, command, partitionedPath.Name);

            return _workDoer.DoWork(WorkResourceRequest.CreateExact(1, 8), job, partitionedPath).Await();
        }

        /// <summary>
        /// Invoke CanvasClean on SmallPedigreeCallset callsets. 
        /// </summary>
        protected List<IFileLocation> InvokeCanvasClean(SmallPedigreeCallset callsets, List<IFileLocation> binnedPaths)
        {
            List<IFileLocation> cleanedPaths = new List<IFileLocation>();
            if (callsets.PedigreeSample.Count != binnedPaths.Count)
                throw new Exception($"Number of output CanvasBin files {binnedPaths.Count} is not equal to the number of Canvas callsets {callsets.PedigreeSample.Count}");
            for (int i = 0; i < callsets.PedigreeSample.Count; i++)
            {
                IFileLocation binnedPath = binnedPaths[i];

                var canvasCallset = new CanvasCallset(callsets.PedigreeSample[i].Sample, callsets.AnalysisDetails, null, null, null);
                Console.WriteLine($"Created callset");

                cleanedPaths.Add(InvokeCanvasClean(canvasCallset, binnedPath).CleanedPath);
                Console.WriteLine($"Run InvokeCanvasClean");

            }
            return cleanedPaths;
        }

        private CanvasCleanOutput InvokeCanvasClean(CanvasCallset callset, IFileLocation binnedPath)
        {
            StringBuilder commandLine = new StringBuilder();
            string executablePath = GetExecutablePath("CanvasClean", commandLine);

            commandLine.AppendFormat("-i \"{0}\" ", binnedPath);
            var tempFolder = callset.SingleSampleCallset.SampleOutputFolder;
            var cleanedPath = tempFolder.GetFileLocation($"{callset.SingleSampleCallset.SampleName}.cleaned");
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

            var command = commandLine.ToString();
            if (_customParameters.ContainsKey("CanvasClean"))
            {
                command = Canvas.CommandOptionsUtilities.MergeCommandLineOptions(
                    command, _customParameters["CanvasClean"], true);
            }
            var job = new JobInfo(executablePath, command, cleanedPath.Name);
            return _workDoer.DoWork(WorkResourceRequest.CreateExact(1, 8), job, new CanvasCleanOutput(cleanedPath, ffpePath)).Await();
        }

        protected void RunSomaticCalling(IFileLocation partitionedPath, CanvasCallset callset, string canvasBedPath,
            string ploidyVcfPath, IFileLocation ffpePath, IFileLocation canvasSnvPath)
        {

            // get somatic SNV output:
            string somaticSnvPath = callset.SomaticVcfPath?.FullName;

            // Prepare and run CanvasSomaticCaller job:
            var cnvVcfPath = callset.SingleSampleCallset.OutputVcfPath;
            StringBuilder commandLine = new StringBuilder();
            var executablePath = GetExecutablePath("CanvasSomaticCaller", commandLine);

            commandLine.Append($" -v {canvasSnvPath}");
            commandLine.Append($" -i {partitionedPath}");
            commandLine.Append($" -o {cnvVcfPath}");
            commandLine.Append($" -b {canvasBedPath}");
            if (!string.IsNullOrEmpty(ploidyVcfPath))
                commandLine.Append($" -p {ploidyVcfPath}");
            commandLine.Append($" -n {callset.SingleSampleCallset.SampleName}");
            if (callset.IsEnrichment)
                commandLine.Append(" -e");
            if (callset.SingleSampleCallset.IsDbSnpVcf) // a dbSNP VCF file is used in place of the normal VCF file
                commandLine.Append(" -d");
            // get localSD metric:
            if (ffpePath != null)
            {
                // Sanity-check: CanvasClean does not always write this file. 
                // If it's not present, just carry on:
                if (ffpePath.Exists)
                {
                    commandLine.Append($" -f \"{ffpePath}\"");
                }
                else
                {
                    Logger.Info("Note: SD file not found at '{0}'", ffpePath);
                }
            }

            if (!string.IsNullOrEmpty(somaticSnvPath))
                commandLine.Append($" -s {somaticSnvPath}");
            commandLine.Append($" -r \"{callset.AnalysisDetails.WholeGenomeFastaFolder}\" ");
            var command = commandLine.ToString();
            if (_customParameters.ContainsKey("CanvasSomaticCaller"))
            {
                command = Canvas.CommandOptionsUtilities.MergeCommandLineOptions(
                    command, _customParameters["CanvasSomaticCaller"], true);
            }
            var job = new JobInfo(executablePath, command, $"SomaticCNV-{callset.SingleSampleCallset.SampleName}");
            _workDoer.DoWork(WorkResourceRequest.CreateExact(1, 8), job, true).Await();
        }

        protected void RunSmallPedigreeCalling(List<IFileLocation> partitionedPaths, SmallPedigreeCallset callsets)
        {

            if (callsets.PedigreeSample.Count != partitionedPaths.Count)
                throw new Exception($"Number of output CanvasPartition files {partitionedPaths.Count} is not equal to the number of Canvas callsets {callsets.PedigreeSample.Count}");

            // CanvasSmallPedigreeCaller:
            StringBuilder commandLine = new StringBuilder();
            string executablePath = GetExecutablePath("CanvasPedigreeCaller", commandLine);

            foreach (IFileLocation partitionedPath in partitionedPaths)
                commandLine.AppendFormat("-i \"{0}\" ", partitionedPath);

            foreach (var callset in callsets.PedigreeSample)
            {
                commandLine.AppendFormat("-v \"{0}\" ", callset.Sample.VfSummaryPath);
                commandLine.AppendFormat("-n \"{0}\" ", callset.Sample.SampleName);
                commandLine.AppendFormat("-t \"{0}\" ", callset.SampleType.ToString());
            }
            var outDir = callsets.AnalysisDetails.OutputFolder;
            commandLine.Append($"-o \"{outDir}\" ");
            commandLine.AppendFormat("-r \"{0}\" ", callsets.AnalysisDetails.WholeGenomeFastaFolder);

            if (callsets.AnalysisDetails.PloidyVcf != null)
                commandLine.AppendFormat("-p \"{0}\" ", callsets.AnalysisDetails.PloidyVcf);

            var command = commandLine.ToString();
            if (_customParameters.ContainsKey("CanvasPedigreeCaller"))
            {
                command = Canvas.CommandOptionsUtilities.MergeCommandLineOptions(
                    command, _customParameters["CanvasPedigreeCaller"], true);
            }
            var job = new JobInfo(executablePath, command, "CanvasPedigreeCaller");
            _workDoer.DoWork(WorkResourceRequest.CreateExact(1, 8), job, true).Await();
        }

        private static string WritePedigreeFile(SmallPedigreeCallset callsets)
        {
            string outFile = Path.Combine(callsets.AnalysisDetails.OutputFolder.FullName, "pedigree.ped");
            string motherSampleName = callsets.PedigreeSample.Where(x => x.SampleType == SampleType.Mother).Select(x => x.Sample.SampleName).SingleOrDefault() ?? "0";
            string fatherSampleName = callsets.PedigreeSample.Where(x => x.SampleType == SampleType.Father).Select(x => x.Sample.SampleName).SingleOrDefault() ?? "0";
            using (FileStream stream = new FileStream(outFile, FileMode.Create, FileAccess.Write))
            using (StreamWriter writer = new StreamWriter(stream))
            {
                foreach (PedigreeSample callset in callsets.PedigreeSample)
                    if (callset.SampleType == SampleType.Mother || callset.SampleType == SampleType.Father)
                        writer.WriteLine($"1\t{callset.Sample.SampleName}\t0\t0\t0\t0");
                    else
                    {
                        string phenotype = callset.SampleType == SampleType.Proband ? "affected" : "0";
                        writer.WriteLine($"1\t{callset.Sample.SampleName}\t{fatherSampleName}\t{motherSampleName}\t0\t{phenotype}");
                    }
            }
            return outFile;
        }

        protected void RunGermlineCalling(IFileLocation partitionedPath, CanvasCallset callset, string ploidyVcfPath, IFileLocation canvasSnvPath)
        {
            StringBuilder commandLine = new StringBuilder();
            ////////////////////////////////////////////////////////
            // CanvasDiploidCaller:
            string executablePath = GetExecutablePath("CanvasDiploidCaller", commandLine);

            commandLine.AppendFormat("-i \"{0}\" ", partitionedPath);
            commandLine.AppendFormat("-v \"{0}\" ", canvasSnvPath);
            var cnvVcfPath = callset.SingleSampleCallset.OutputVcfPath;
            commandLine.AppendFormat("-o \"{0}\" ", cnvVcfPath);
            commandLine.AppendFormat("-n \"{0}\" ", callset.SingleSampleCallset.SampleName);
            commandLine.AppendFormat("-r \"{0}\" ", callset.AnalysisDetails.WholeGenomeFastaFolder);
            if (!string.IsNullOrEmpty(ploidyVcfPath))
            {
                commandLine.AppendFormat("-p \"{0}\" ", ploidyVcfPath);
            }
            if (callset.SingleSampleCallset.IsDbSnpVcf) // a dbSNP VCF file is used in place of the normal VCF file
                commandLine.AppendFormat("-d ");
            
            var command = commandLine.ToString();
            if (_customParameters.ContainsKey("CanvasDiploidCaller"))
            {
                command = Canvas.CommandOptionsUtilities.MergeCommandLineOptions(
                    command, _customParameters["CanvasDiploidCaller"], true);
            }
            var job = new JobInfo(executablePath, command, cnvVcfPath.Name);
            _workDoer.DoWork(WorkResourceRequest.CreateExact(1, 8), job, true).Await();
        }
    }
}
