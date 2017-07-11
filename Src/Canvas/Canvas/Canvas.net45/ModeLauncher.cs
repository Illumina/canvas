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
                        var executableProcessor = new NullExecutableProcessor(logger, workerDirectory);
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

    public class NullExecutableProcessor : IExecutableProcessor
    {
        private ILogger _logger;
        protected readonly List<IDirectoryLocation> _execDirs;

        public NullExecutableProcessor(ILogger log, params IDirectoryLocation[] binaryDirs) 
        {
            _logger = log;
            _execDirs = binaryDirs.ToList();
        }

        public string GetExecutableParameters(string executableFileNameWithoutExtension)
        {
            return null;
        }

        public string GetExecutable(string executableFileNameWithoutExtension, params string[] relativeBinDirectory)
        {
            return GetExecutableFileLocation(executableFileNameWithoutExtension, relativeBinDirectory).FullName;
        }

        public string GetEnvironmentExecutablePath(string executableFileNameWithoutExtension)
        {
            IFileLocation path;

            // Now check the environment paths
            var paths = Environment.GetEnvironmentVariable("PATH").Split(Path.PathSeparator).Select(dir => new DirectoryLocation(dir));
            foreach (var candidatePath in paths)
            {
                var fullCandidatePath =
                    GetPlatformDependentExecutablePath(candidatePath.GetFileLocation(executableFileNameWithoutExtension));

                if (fullCandidatePath.Exists)
                {
                    return fullCandidatePath.FullName;
                }
            }
            throw new Exception("Could not find executable '" + executableFileNameWithoutExtension + "' in environment path");
        }

        public bool TryGetCustomExecutableFileLocation(string executableFileNameWithoutExtension, out IFileLocation path)
        {
            path = null;
            return false;
        }
        public bool TryGetCustomExecutablePath(string executableFileNameWithoutExtension, out string path)
        {
            path = null;
            return false;
        }


        public IFileLocation GetExecutableFileLocation(string executableFileNameWithoutExtension, params string[] executableRelativePath)
        {
            IFileLocation path;

            IDirectoryLocation SearchDir(IDirectoryLocation execDir) => executableRelativePath.Aggregate(execDir,
                (dir, relativePath) => dir.GetDirectoryLocation(relativePath), result => result);

            var execDirs = _execDirs.Select(SearchDir).ToList();
            foreach (var execDir in execDirs)
            {
                path = GetPlatformDependentExecutablePath(
                    execDir.GetFileLocation(executableFileNameWithoutExtension));
                if (path.Exists) return path;
            }
            string searchedDirs = string.Join(", ", execDirs);
            throw new Exception($"Could not find executable {executableFileNameWithoutExtension} in {searchedDirs}");
        }

        private IFileLocation GetPlatformDependentExecutablePath(IFileLocation executablePathWithoutExtension)
        {
            var path = executablePathWithoutExtension.AppendName(".exe");
            if (path.Exists && !CrossPlatform.IsThisLinux())
            {
                return path;
            }
            path = executablePathWithoutExtension.AppendName(".dll");
            if (path.Exists)
            {
                return path;
            }
            return executablePathWithoutExtension; // fallback
        }

    }
    
}
