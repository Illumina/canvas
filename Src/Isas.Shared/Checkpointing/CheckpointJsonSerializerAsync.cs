using System;
using System.IO;
using Illumina.SecondaryAnalysis;
using ILMNcommon.Common;
using Newtonsoft.Json;

namespace Isas.Shared.Checkpointing
{
    public interface ICheckpointSerializerAsync
    {
        bool CanLoad(Checkpoint checkpoint);
        void Save<T>(Checkpoint checkpoint, T output);
        T Load<T>(Checkpoint checkpoint);
        ICheckpointSerializerAsync GetChildSerializer(IDirectoryLocation childRepository);
    }

    public class CheckpointJsonSerializerAsync : ICheckpointSerializerAsync
    {
        private readonly JsonSerializerSettings _settings;
        private readonly IDirectoryLocation _repository;
        private readonly ILogger _logger;
        private CheckpointJsonSerializerAsync(JsonSerializerSettings settings, IDirectoryLocation repository, ILogger logger)
        {
            _settings = settings;
            _repository = repository;
            _logger = logger;
        }

        public CheckpointJsonSerializerAsync(IDirectoryLocation repository, ILogger logger, params JsonConverter[] converters) :
            this(CheckpointManagerFactory.GetJsonSerializerSettings(converters), repository, logger)
        {
        }

        public bool CanLoad(Checkpoint checkpoint)
        {
            return GetCheckpointFile(checkpoint).Exists;
        }

        public void Save<T>(Checkpoint checkpoint, T output)
        {
            IFileLocation path = GetCheckpointFile(checkpoint);
            _logger.Info("Saving checkpoint results to {0}", path);
            Serialize(path, output);
        }

        public T Load<T>(Checkpoint checkpoint)
        {
            IFileLocation path = GetCheckpointFile(checkpoint);
            if (!path.Exists)
                throw new ApplicationException($"Unable to load checkpoint results. Missing expected checkpoint file at {path}");
            _logger.Info("Loading checkpoint results from {0}", path);
            return Deserialize<T>(path);
        }

        private IFileLocation GetCheckpointFile(Checkpoint checkpoint)
        {
            return _repository.GetFileLocation($"{checkpoint.Number:00}-{checkpoint.Name.RemoveWhiteSpace()}.json");
        }

        private void Serialize<T>(IFileLocation path, T obj)
        {
            path.Directory.Create();
            using (StreamWriter writer = new StreamWriter(path.OpenWrite()))
                writer.Write(JsonConvert.SerializeObject(obj, Formatting.Indented, _settings));
        }

        private T Deserialize<T>(IFileLocation path)
        {
            using (StreamReader reader = path.OpenText())
                return JsonConvert.DeserializeObject<T>(reader.ReadToEnd(), _settings);
        }

        public ICheckpointSerializerAsync GetChildSerializer(IDirectoryLocation childRepository)
        {
            return new CheckpointJsonSerializerAsync(_settings, childRepository, _logger);
        }
    }
}