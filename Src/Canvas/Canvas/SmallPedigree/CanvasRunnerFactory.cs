using System;
using System.Collections.Generic;
using Canvas.Visualization;
using CanvasCommon;
using CanvasCommon.Visualization;
using Illumina.Common.FileSystem;
using Isas.ClassicBioinfoTools.Tabix;
using Isas.Framework.Checkpointing;
using Isas.Framework.Logging;
using Isas.Framework.Settings;
using Isas.Framework.WorkManagement;
using Isas.Framework.WorkManagement.CommandBuilding;

namespace Canvas.SmallPedigree
{
    public class CanvasRunnerFactory
    {
        private readonly ILogger _logger;
        private readonly ICheckpointRunner _checkpointRunner;
        private readonly IWorkManager _workManager;
        private readonly IWorkDoer _workDoer;
        private readonly IFileLocation _runtimeExecutable;
        private readonly Func<string, ICommandFactory> _runtimeCommandPrefix;

        public CanvasRunnerFactory(ILogger logger, ICheckpointRunner checkpointRunner, IWorkManager workManager,
            IWorkDoer workDoer, IFileLocation runtimeExecutable, Func<string, ICommandFactory> runtimeCommandPrefix)
        {
            _logger = logger;
            _checkpointRunner = checkpointRunner;
            _workManager = workManager;
            _workDoer = workDoer;
            _runtimeExecutable = runtimeExecutable;
            _runtimeCommandPrefix = runtimeCommandPrefix;
        }

        public CanvasRunner Create(bool isSomatic, CanvasCoverageMode coverageMode,
            int countsPerBin, Dictionary<string, string> customParameters)
        {
            var settings = IsasConfigurationSettings.GetConfigSettings();
            var canvasFolder = new DirectoryLocation(Isas.Framework.Utilities.Utilities.GetAssemblyFolder(typeof(CanvasRunner)));

            var commandManager = new CommandManager(new ExecutableProcessor(settings, _logger, canvasFolder));
            var tabixWrapper = TabixWrapperFactory.GetTabixWrapper(_logger, _workDoer, commandManager);
            var bAlleleBedGraphWriter = new BAlleleBedGraphWriter(new BgzfBedGraphWriter(new RoundingBedGraphWriter(new BedGraphWriterFacade(), 4), tabixWrapper));

            return new CanvasRunner(_logger, _workManager, _workDoer, _checkpointRunner, _runtimeExecutable, _runtimeCommandPrefix, isSomatic, coverageMode, countsPerBin, bAlleleBedGraphWriter, customParameters, canvasFolder.FullName);
        }
    }
}