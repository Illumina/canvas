using System;
using System.Collections.Generic;
using Canvas.CommandLineParsing;
using Illumina.Common;
using Illumina.Common.FileSystem;
using Isas.Framework.Checkpointing;
using Isas.Framework.DataTypes;
using Isas.Framework.FrameworkFactory;
using Isas.Framework.Logging;
using Isas.Framework.Settings;
using Isas.Framework.WorkManagement;

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
            IDirectoryLocation genomeRoot = commonOptions.WholeGenomeFasta?.Parent?.Parent?.Parent?.Parent?.Parent;
            int returnValue = 0;
            IsasFrameworkFactory.RunWithIsasFramework(outFolder, log, error, commonOptions.StartCheckpoint, commonOptions.StopCheckpoint, 0,
                config.MaximumMemoryGB, config.MaximumHoursPerProcess, false, genomeRoot,
                frameworkServices =>
                {
                    var logger = frameworkServices.Logger;
                    try
                    {
                        var executableProcessor = new ExecutableProcessor(new NullSampleSettings(), logger);
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

    public class NullSampleSettings : ISampleSettings
    {
        /// <summary>
        /// Returns the value associated with the given key in the Header section.
        /// </summary>
        /// <param name="key">The key to look up.</param>
        /// <returns>The value associated with the key, or null if key not present.</returns>
        public string GetHeader(string key)
        {
            return null;
        }

        /// <summary>
        /// Returns the value associated with the given key in the Settings section.
        /// </summary>
        /// <param name="key">The key to look up.</param>
        /// <returns>The value associated with the key, or null if key not present.</returns>
        public string GetSetting(string key)
        {
            return null;
        }

        /// <summary>Returns all the setting keys.</summary>
        /// <returns></returns>
        public IEnumerable<string> GetSettingKeys()
        {
            return null;
        }

        /// <summary>
        /// </summary>
        /// <param name="key"></param>
        /// <returns>SampleSet&lt;List&lt;string&gt;&gt;, or SampleSet&lt;List&lt;(string)null&gt;&gt; if key not present.</returns>
        public SampleSet<List<string>> GetDataColumn(string key)
        {
            throw new NotImplementedException();
        }

        public SampleSet<SampleInfo> GetSamples()
        {
            throw new NotImplementedException();
        }

        public string GetManifest(string key)
        {
            throw new NotImplementedException();
        }

        public IEnumerable<string> GetManifestKeys()
        {
            throw new NotImplementedException();
        }

        public List<string> GetSection(string sectionName)
        {
            throw new NotImplementedException();
        }

        public void CheckUnusedEntries()
        {
            throw new NotImplementedException();
        }

        public string SampleSheetPath => null;
    }
}
