using System;
using System.Collections.Generic;
using Illumina.SecondaryAnalysis;
using Newtonsoft.Json;
using SampleSettingsProcessing;

namespace Isas.Shared.Checkpointing
{
    public class CheckpointManagerFactory
    {
        private readonly ISampleSettings _settings;
        private readonly ILogger _logger;
        private readonly IDirectoryLocation _genomesRoot;
        private readonly bool _retainTemps;
        public static string RunStepsSerially = "RunStepsSerially";
        public CheckpointManagerFactory(ISampleSettings settings, ILogger logger, IDirectoryLocation genomesRoot, bool retainTemps)
        {
            _settings = settings;
            _logger = logger;
            _genomesRoot = genomesRoot;
            _retainTemps = retainTemps;
        }

        public ICheckpointRunner GetCheckpointManager(IDirectoryLocation analysisFolder,
            string startingCheckpointName, string stopCheckpointName)
        {
            // Checkpoint Manager
            ICheckpointSerializer serializer = new CheckpointJsonSerializer(GetCheckpointFolder(analysisFolder), _logger, GetJsonConverters(analysisFolder));
            ICheckpointManager manager = new CheckpointManager(_logger, GetStartingCheckpointName(startingCheckpointName), stopCheckpointName);
            return new CheckpointRunner(_logger, analysisFolder, manager, serializer, _retainTemps);
        }

        private static string GetStartingCheckpointName(string startingCheckpointName)
        {
            //todo: remove this increment logic once the legacy checkpointer is no longer used
            int tempStartCheckpoint;
            if (Int32.TryParse(startingCheckpointName, out tempStartCheckpoint))
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

        public ICheckpointRunnerAsync GetAsyncCheckpointManager(IDirectoryLocation analysisFolder,
            string startingCheckpointName, string stopCheckpointName, bool defaultRunStepsSerially = false)
        {
            bool runStepsSerially = _settings.GetBooleanSetting(RunStepsSerially, defaultRunStepsSerially);
            string message = $"Workflow steps will {(runStepsSerially ? "not" : "")} run in parallel";
            _logger.Info(message);
            if (runStepsSerially != defaultRunStepsSerially)
            {
                _logger.Warn($"Overriding default value for setting '{RunStepsSerially}'. {message}");
            }
            var serializer = new CheckpointJsonSerializerAsync(GetCheckpointFolder(analysisFolder), _logger, GetJsonConverters(analysisFolder));
            return CheckpointRunnerAsync.Create(serializer, _logger, analysisFolder, GetStartingCheckpointName(startingCheckpointName), stopCheckpointName, _retainTemps, runStepsSerially);
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

        public static IDirectoryLocation GetCheckpointFolder(IDirectoryLocation analysisFolder)
        {
            return analysisFolder.CreateSubdirectory("Checkpoints");
        }
    }
}
