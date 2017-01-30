using System;
using System.Collections.Generic;
using Canvas.CommandLineParsing;
using Illumina.Common.FileSystem;
using Isas.Framework.Checkpointing;
using Isas.Framework.Logging;
using Isas.Framework.Settings;
using Isas.Framework.Utilities;
using Isas.Framework.WorkManagement;
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
        void Run(ILogger logger, ICheckpointRunner checkpointRunner, IWorkManager workManager);
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
                // Get work manager
                IDirectoryLocation loggingFolder = outFolder.CreateSubdirectory("Logging");
                IsasConfiguration config = IsasConfiguration.GetConfiguration();
                IWorkManager workManager = new LocalWorkManager(logger, loggingFolder, 0, config.MaximumMemoryGB, config.MaximumHoursPerProcess);
                try
                {
                    logger.Info($"Running Canvas {_mode} {_version}");
                    logger.Info($"Command-line arguments: {string.Join(" ", _args)}");
                    // Manager factory
                    CheckpointManagerFactory checkpointManagerFactory = new CheckpointManagerFactory(logger, commonOptions.StartCheckpoint, commonOptions.StopCheckpoint);
                    IDirectoryLocation genomeRoot = commonOptions.WholeGenomeFasta?.Parent?.Parent?.Parent?.Parent?.Parent;
                    CheckpointSerializerFactory serializerFactory = new CheckpointSerializerFactory(logger, outFolder, genomeRoot);

                    checkpointManagerFactory.RunWithArgument(checkpointManager =>
                    {
                        var basicCheckpointerFactory = new SerializingCheckpointerFactory(logger, checkpointManager, serializerFactory.GetJsonSerializer());
                        SandboxCheckpointerFactory checkpointerFactory = new SandboxCheckpointerFactory(basicCheckpointerFactory, outFolder, true);
                        checkpointerFactory.RunWithArgument(checkpointer =>
                        {
                            _modeRunner.Run(logger, checkpointer, workManager);
                        });
                    });
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
    }

    public class NullModeLauncher : IModeLauncher
    {
        public int Launch()
        {
            return 0;
        }
    }
}