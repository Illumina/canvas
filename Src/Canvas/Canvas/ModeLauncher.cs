﻿using Illumina.Common.FileSystem;
using Isas.Framework.Checkpointing;
using Isas.Framework.Logging;
using Isas.Framework.Settings;
using Isas.Framework.WorkManagement;
using System;
using System.Collections.Generic;
using Isas.Framework.WorkManagement.CommandBuilding;

namespace Canvas
{
    public interface IModeLauncher
    {
        int Launch();
    }

    public interface IModeRunner
    {
        void Run(ILogger logger, ICheckpointRunner checkpointRunner, IWorkManager workManager, IWorkDoer workDoer, IFileLocation runtimeExecutable, Func<string, ICommandFactory> runtimeCommandPrefix);
    }

    public class ModeLauncher : IModeLauncher
    {
        private readonly ILogger _logger;
        private readonly ISettings _settings;
        private readonly ICheckpointRunner _checkpointRunner;
        private readonly IWorkManager _workManager;
        private readonly IWorkDoer _workDoer;
        private readonly IModeRunner _modeRunner;
        private readonly IEnumerable<string> _args;
        private readonly string _version;
        private readonly string _mode;

        public ModeLauncher(ILogger logger, ISettings settings, ICheckpointRunner checkpointRunner, IWorkManager workManager, IWorkDoer workDoer, IModeRunner modeRunner, IEnumerable<string> args, string version, string mode)
        {
            _logger = logger;
            _settings = settings;
            _checkpointRunner = checkpointRunner;
            _workManager = workManager;
            _workDoer = workDoer;
            _modeRunner = modeRunner;
            _args = args;
            _version = version;
            _mode = mode;
        }

        public int Launch()
        {
            try
            {
                var executableProcessor = new ExecutableProcessor(_settings, _logger);
                var runtimeExecutable = new FileLocation(executableProcessor.GetEnvironmentExecutablePath("dotnet"));

                ICommandFactory GetDotnetCommand(string component)
                {

                    return new CommandManager(executableProcessor).GetDotnetCommandFactory($"{component}.dll",component);
                }
                _logger.Info($"Running Canvas {_mode} {_version}");
                _logger.Info($"Command-line arguments: {string.Join(" ", _args)}");
                _modeRunner.Run(_logger, _checkpointRunner, _workManager, _workDoer, runtimeExecutable, GetDotnetCommand);
                return 0;
            }
            catch (Exception e)
            {
                _logger.Error($"Canvas workflow error: {e}");
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
