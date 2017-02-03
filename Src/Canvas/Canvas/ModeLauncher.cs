using System;
using System.Collections.Generic;
using Canvas.CommandLineParsing;
using Illumina.Common.FileSystem;
using Isas.Framework.Checkpointing;
using Isas.Framework.FrameworkFactory;
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
            IsasConfiguration config = IsasConfiguration.GetConfiguration();
            IDirectoryLocation genomeRoot = commonOptions.WholeGenomeFasta?.Parent?.Parent?.Parent?.Parent?.Parent;
            int returnValue = 0;
            IsasFrameworkFactory.RunWithIsasFramework(outFolder, log, error, commonOptions.StartCheckpoint, commonOptions.StopCheckpoint, 0,
                config.MaximumMemoryGB, config.MaximumHoursPerProcess, false, genomeRoot,
                frameworkServices =>
                {
                    var logger = frameworkServices.Logger;
                    try
                    {
                        frameworkServices.Logger.Info($"Running Canvas {_mode} {_version}");
                        logger.Info($"Command-line arguments: {string.Join(" ", _args)}");
                        _modeRunner.Run(logger, frameworkServices.Checkpointer, frameworkServices.WorkManager);
                    }
                    catch (Exception e)
                    {
                        logger.Error($"Canvas workflow error: {e}");
                        returnValue = -1;
                    }
                    returnValue = 0;
                });
            return returnValue;
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
