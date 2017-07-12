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
        CommonOptions CommonOptions { get; }
        void Run(ILogger logger, ICheckpointRunner checkpointRunner, IWorkManager workManager, IFileLocation runtimeExecutable);
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
            // From $Root/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta to $Root
            IDirectoryLocation genomeRoot = commonOptions.WholeGenomeFasta?.Parent?.Parent?.Parent?.Parent?.Parent;
            int returnValue = 0;
            IsasFrameworkFactory.RunWithIsasFramework(outFolder, log, error, commonOptions.StartCheckpoint, commonOptions.StopCheckpoint, 0,
                config.MaximumMemoryGB, config.MaximumHoursPerProcess, false, genomeRoot,
                frameworkServices =>
                {
                    var logger = frameworkServices.Logger;
                    try
                    {
                        var workerDirectory = new DirectoryLocation(Isas.Framework.Utilities.Utilities.GetAssemblyFolder(typeof(ModeLauncher)));
                        var executableProcessor = new ExecutableProcessor((ISampleSettings)(new SettingsProcessor()), logger, workerDirectory);
                        var runtimeExecutable = new FileLocation(executableProcessor.GetEnvironmentExecutablePath("dotnet"));
                        frameworkServices.Logger.Info($"Running Canvas {_mode} {_version}");
                        logger.Info($"Command-line arguments: {string.Join(" ", _args)}");
                        _modeRunner.Run(logger, frameworkServices.Checkpointer, frameworkServices.WorkManager, runtimeExecutable);
                        returnValue = 0;
                    }
                    catch (Exception e)
                    {
                        logger.Error($"Canvas workflow error: {e}");
                        returnValue = -1;
                    }
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
