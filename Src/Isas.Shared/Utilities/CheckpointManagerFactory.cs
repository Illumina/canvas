using System.Collections.Generic;
using Illumina.SecondaryAnalysis;
using Illumina.SecondaryAnalysis.Workflow;
using Isas.Shared;
using Newtonsoft.Json;

namespace Isas.Shared
{
    public class CheckpointManagerFactory
    {
        private readonly ILogger _logger;
        private readonly IDirectoryLocation _genomesRoot;
        private readonly bool _retainTemps;

        public CheckpointManagerFactory(ILogger logger, IDirectoryLocation genomesRoot, bool retainTemps)
        {
            _logger = logger;
            _genomesRoot = genomesRoot;
            _retainTemps = retainTemps;
        }

        public ICheckpointRunner GetCheckpointManager(IDirectoryLocation analysisFolder,
            string startingCheckpointName, string stopCheckpointName)
        {
            // Checkpoint Manager
            var parentDirectories = new Dictionary<string, IDirectoryLocation>
            {
                {"AnalysisFolder", analysisFolder},
                {"GenomesRoot", _genomesRoot}
            };
            JsonConverter[] converters =
            {
                new FileSystemLocationConverter(parentDirectories),
                new EnumerableConverter()
            };

            var settings = new JsonSerializerSettings
            {
                TypeNameHandling = TypeNameHandling.Auto,
                MissingMemberHandling = MissingMemberHandling.Error,
                Converters = converters,
                ContractResolver = new WritablePropertiesOnlyResolver(),
                ConstructorHandling = ConstructorHandling.AllowNonPublicDefaultConstructor,
            };
            ICheckpointSerializer serializer = new CheckpointJsonSerializer(settings, IsasFilePaths.GetCheckpointFolder(analysisFolder), _logger);

            //todo: remove this increment logic once the legacy checkpointer is no longer used
            int tempStartCheckpoint;
            if (int.TryParse(startingCheckpointName, out tempStartCheckpoint))
            {
                tempStartCheckpoint++;
                startingCheckpointName = tempStartCheckpoint.ToString();
            }
            ICheckpointManager manager = new CheckpointManager(_logger, startingCheckpointName, stopCheckpointName);
            return new CheckpointRunner(_logger, analysisFolder, manager, serializer, _retainTemps);
        }
    }
}
