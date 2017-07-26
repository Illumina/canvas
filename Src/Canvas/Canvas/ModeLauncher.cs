using Canvas.CommandLineParsing;
using Illumina.Common.FileSystem;
using Isas.Framework.Checkpointing;
using Isas.Framework.FrameworkFactory;
using Isas.Framework.Logging;
using Isas.Framework.Settings;
using Isas.Framework.WorkManagement;
using Illumina.Common;
using System;
using System.IO;
using System.Linq;
using System.Collections.Generic;

namespace Canvas
{
    public interface IModeLauncher
    {
        int Launch();
    }

    public interface IModeRunner
    {
        void Run(ILogger logger, ICheckpointRunner checkpointRunner, IWorkManager workManager, IFileLocation runtimeExecutable);
    }

    public class ModeLauncher : IModeLauncher
    {
        private readonly FrameworkServices _frameworkServices;
        private readonly CommonOptions _commonOptions;
        private readonly IModeRunner _modeRunner;
        private readonly IEnumerable<string> _args;
        private readonly string _version;
        private readonly string _mode;

        public ModeLauncher(FrameworkServices frameworkServices, CommonOptions commonOptions, IModeRunner modeRunner, IEnumerable<string> args, string version, string mode)
        {
            _frameworkServices = frameworkServices;
            _commonOptions = commonOptions;
            _modeRunner = modeRunner;
            _args = args;
            _version = version;
            _mode = mode;
        }
        public int Launch()
        {
            var logger = _frameworkServices.Logger;
            try
            {
                var executableProcessor = new ExecutableProcessor(new SettingsProcessor(), logger);
                var runtimeExecutable = new FileLocation(executableProcessor.GetEnvironmentExecutablePath("dotnet"));
                _frameworkServices.Logger.Info($"Running Canvas {_mode} {_version}");
                logger.Info($"Command-line arguments: {string.Join(" ", _args)}");
                _modeRunner.Run(logger, _frameworkServices.Checkpointer, _frameworkServices.WorkManager, runtimeExecutable);
                return 0;
            }
            catch (Exception e)
            {
                logger.Error($"Canvas workflow error: {e}");
                return -1;
            }
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
