using System;
using System.Collections.Generic;
using Canvas.CommandLineParsing;
using Illumina.SecondaryAnalysis;
using Isas.Shared;
using Isas.Shared.Checkpointing;
using Isas.Shared.Utilities;
using Newtonsoft.Json;

namespace Canvas
{
    public interface IModeLauncher
    {
        int Launch();
    }

    public interface IModeRunner
    {
        CommonOptions CommonOptions { get; }
        void Run(ILogger logger, ICheckpointRunnerAsync checkpointRunner, IWorkManager workManager);
    }

    public class ModeLauncher : IModeLauncher
    {
        private readonly IModeRunner _modeRunner;
        private readonly IEnumerable<string> _args;
        private readonly string _version;
        private readonly string _mode;

        public ModeLauncher(IModeRunner modeRunner, IEnumerable<string> args, string version, string mode)
        {
            _modeRunner = modeRunner;
            _args = args;
            _version = version;
            _mode = mode;
        }

        public int Launch()
        {
            CommonOptions commonOptions = _modeRunner.CommonOptions;
            IDirectoryLocation outFolder = commonOptions.OutputDirectory;
            var log = outFolder.GetFileLocation("CanvasLog.txt");
            var error = outFolder.GetFileLocation("CanvasError.txt");
            using (ILogger logger = new Logger(log, error))
            {
                try
                {
                    logger.Info($"Running Canvas {_mode} {_version}");
                    logger.Info($"Command-line arguments: {string.Join(" ", _args)}");
                    using (ICheckpointRunnerAsync checkpointRunner =
                        GetCheckpointRunner(
                            logger,
                            outFolder,
                            commonOptions.StartCheckpoint,
                            commonOptions.StopCheckpoint,
                            commonOptions.WholeGenomeFasta))
                    {
                        IDirectoryLocation loggingFolder = outFolder.CreateSubdirectory("Logging");
                        IsasConfiguration config = IsasConfiguration.GetConfiguration();
                        IWorkManager workManager = new LocalWorkManager(logger, loggingFolder, 0, config.MaximumMemoryGB, config.MaximumHoursPerProcess);
                        _modeRunner.Run(logger, checkpointRunner, workManager);
                    }
                }
                catch (StopCheckpointFoundException) { }
                catch (Exception e)
                {
                    logger.Error($"Canvas workflow error: {e}");
                    return -1;
                }
            }
            return 0;
        }

        private ICheckpointRunnerAsync GetCheckpointRunner(ILogger logger, IDirectoryLocation outputDirectory, string startCheckpoint, string stopCheckpoint, IDirectoryLocation wholeGenomeFastaFolder)
        {

            var parentDirectories = new Dictionary<string, IDirectoryLocation>
            {
                {"Output", outputDirectory},
            };
            // Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta
            IDirectoryLocation genomeRoot = wholeGenomeFastaFolder?.Parent?.Parent?.Parent?.Parent?.Parent;
            if (genomeRoot != null) parentDirectories.Add("Genome", genomeRoot);

            JsonConverter[] converters = { new FileSystemLocationConverter(parentDirectories) };

            ICheckpointSerializerAsync serializer = new CheckpointJsonSerializerAsync(CheckpointManagerFactory.GetCheckpointFolder(outputDirectory), logger, converters);
            return CheckpointRunnerAsync.Create(serializer, logger, outputDirectory, startCheckpoint, stopCheckpoint, true);
        }
    }

    public class NullModeLauncher : IModeLauncher
    {
        public int Launch()
        {
            return 0;
        }
    }
}