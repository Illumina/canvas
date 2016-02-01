using System.Collections.Generic;
using Illumina.SecondaryAnalysis;
using Newtonsoft.Json;

namespace Isas.Shared.Checkpointing
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
            ICheckpointSerializer serializer = new CheckpointJsonSerializer(IsasFilePaths.GetCheckpointFolder(analysisFolder), _logger, GetJsonConverters(analysisFolder));
            ICheckpointManager manager = new CheckpointManager(_logger, GetStartingCheckpointName(startingCheckpointName), stopCheckpointName);
            return new CheckpointRunner(_logger, analysisFolder, manager, serializer, _retainTemps);
        }

        private static string GetStartingCheckpointName(string startingCheckpointName)
        {
            //todo: remove this increment logic once the legacy checkpointer is no longer used
            int tempStartCheckpoint;
            if (int.TryParse(startingCheckpointName, out tempStartCheckpoint))
            {
                tempStartCheckpoint++;
                startingCheckpointName = tempStartCheckpoint.ToString();
            }
            return startingCheckpointName;
        }

        private JsonConverter[] GetJsonConverters(IDirectoryLocation analysisFolder)
        {
            var parentDirectories = new Dictionary<string, IDirectoryLocation>
            {
                {
                    "AnalysisFolder", analysisFolder
                },
                {
                    "GenomesRoot", _genomesRoot
                }
            };
            JsonConverter[] converters =
            {
                new FileSystemLocationConverter(parentDirectories),
                new EnumerableConverter()
            };
            return converters;
        }

        public ICheckpointRunner GetAsyncCheckpointManager(IDirectoryLocation analysisFolder,
            string startingCheckpointName, string stopCheckpointName)
        {
            var serializer = new CheckpointJsonSerializerAsync(IsasFilePaths.GetCheckpointFolder(analysisFolder), _logger, GetJsonConverters(analysisFolder));
            return CheckpointRunnerAsync.Create(serializer, _logger, analysisFolder, GetStartingCheckpointName(startingCheckpointName), stopCheckpointName, _retainTemps);
        }

        public static JsonSerializerSettings GetJsonSerializerSettings(params JsonConverter[] converters)
        {
            return new JsonSerializerSettings
            {
                TypeNameHandling = TypeNameHandling.Auto,
                MissingMemberHandling = MissingMemberHandling.Error,
                Converters = converters,
                ContractResolver = new WritablePropertiesOnlyResolver(),
                ConstructorHandling = ConstructorHandling.AllowNonPublicDefaultConstructor,
            };
        }
    }
}
