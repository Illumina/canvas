using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using Illumina.SecondaryAnalysis;
using Newtonsoft.Json;
using Newtonsoft.Json.Serialization;

namespace Isas.Shared.Checkpointing
{
    public interface ICheckpointSerializer
    {
        void Save<T>(string checkpoint, T input);
        T Load<T>(string checkpoint);
        bool Exists(string checkpoint);
    }

    public class CheckpointJsonSerializer : ICheckpointSerializer
    {
        private readonly JsonSerializerSettings _settings;
        private readonly IDirectoryLocation _repository;
        private readonly ILogger _logger;

        private CheckpointJsonSerializer(JsonSerializerSettings settings, IDirectoryLocation repository, ILogger logger)
        {
            _logger = logger;
            _repository = repository;
            if (!_repository.Exists)
                _repository.Create();
            _settings = settings;
        }
        public CheckpointJsonSerializer(IDirectoryLocation repository, ILogger logger, params JsonConverter[] converters) :
            this(CheckpointManagerFactory.GetJsonSerializerSettings(converters), repository, logger)
        {
        }

        public void Save<T>(string checkpoint, T input)
        {
            IFileLocation path = GetCurrentCheckpointFile(checkpoint);
            _logger.Info("Saving checkpoint results to {0}", path);
            Serialize(path, input);
        }

        public T Load<T>(string checkpoint)
        {
            IFileLocation path = GetCurrentCheckpointFile(checkpoint);
            if (!Exists(checkpoint))
                throw new ApplicationException(string.Format("Unable to load checkpoint results. Missing expected checkpoint file at {0}", path));
            _logger.Info("Loading checkpoint results from {0}", path);
            return Deserialize<T>(path);
        }

        public bool Exists(string checkpoint)
        {
            IFileLocation path = GetCurrentCheckpointFile(checkpoint);
            if (!path.Exists) return false;
            return true;
        }

        private void Serialize<T>(IFileLocation path, T obj)
        {
            using (StreamWriter writer = new StreamWriter(path.OpenWrite()))
                writer.Write(JsonConvert.SerializeObject(obj, Formatting.Indented, _settings));
        }

        private T Deserialize<T>(IFileLocation path)
        {
            using (StreamReader reader = path.OpenText())
                return JsonConvert.DeserializeObject<T>(reader.ReadToEnd(), _settings);
        }

        private IFileLocation GetCurrentCheckpointFile(string checkpoint)
        {
            return _repository.GetFileLocation(checkpoint + ".json");
        }
    }

    /// <summary>
    /// Serialization settings:
    ///     Include:
    ///         - public/private auto properties with setters
    ///         - public/private fields that are not the backing fields for auto properties with setters
    ///         - private backing fields for auto properties without setters (C# 6.0 only)
    ///     Exclude:
    ///         - private backing fields for auto properties with setters (we serialize the property instead) 
    ///         - properties without setters
    /// </summary>
    public class WritablePropertiesOnlyResolver : DefaultContractResolver
    {
        protected override IList<JsonProperty> CreateProperties(Type type, MemberSerialization memberSerialization)
        {
            var props = type.GetProperties(BindingFlags.Public | BindingFlags.NonPublic | BindingFlags.Instance)
                            .Where(IsAutoPropertyWithSetter)
                            .Select(property => base.CreateProperty(property, memberSerialization));

            var fields = type.GetFields(BindingFlags.Public | BindingFlags.NonPublic | BindingFlags.Instance)
                             .Where(field => !IsBackingFieldForAutoPropertyWithSetter(field))
                             .Select(field => base.CreateProperty(field, memberSerialization));

            var jsonProperties = props.Union(fields).ToList();
            jsonProperties.ForEach(p => { p.Writable = true; p.Readable = true; });
            return jsonProperties;
        }

        private static bool IsAutoPropertyWithSetter(PropertyInfo property)
        {
            return property.DeclaringType.GetFields(BindingFlags.Public | BindingFlags.NonPublic | BindingFlags.Instance)
                                     .Any(field => IsBackingFieldForAutoPropertyWithSetter(field, property));
        }

        private static bool IsBackingFieldForAutoPropertyWithSetter(FieldInfo field)
        {

            return field.DeclaringType.GetProperties(BindingFlags.Public | BindingFlags.NonPublic | BindingFlags.Instance)
                                      .Any(property => IsBackingFieldForAutoPropertyWithSetter(field, property));
        }

        private static bool IsBackingFieldForAutoPropertyWithSetter(FieldInfo field, PropertyInfo property)
        {
            Type fieldType = field.DeclaringType;
            Type propertyType = property.DeclaringType;
            if (fieldType != propertyType) return false;
            bool hasSetter = property.GetSetMethod(true) != null;
            return hasSetter && field.Name.Contains("<" + property.Name + ">");
        }
    }

    public interface IAbsolutePathInstantiator
    {
        object InstantiateAbsolutePath(string path, Type objectType);
    }

    public class RealFileAbsolutePathInstantiator : IAbsolutePathInstantiator
    {
        public object InstantiateAbsolutePath(string path, Type objectType)
        {
            if (typeof(IFileLocation).IsAssignableFrom(objectType))
            {
                return new FileLocation(path);
            }
            if (typeof(IDirectoryLocation).IsAssignableFrom(objectType))
            {
                return new DirectoryLocation(path);
            }
            throw new ApplicationException("Tried to deserialize something that wasn't IFileLocation or IDirectoryLocation.");
        }
    }
    public class FileSystemLocationConverter : JsonConverter
    {
        private static readonly Type[] Types = { typeof(IFileLocation), typeof(IDirectoryLocation) };
        private readonly Dictionary<string, IDirectoryLocation> _parentDirectories;
        private readonly IAbsolutePathInstantiator _absolutePathInstantiator;
        /// <summary>
        /// Constructor for custom FileSystemLocationConverter
        /// </summary>
        /// <param name="parentDirectories">List of parent directories to use when serializing to store absolute paths as relative paths.
        /// If a location starts with one of the parent directories, the parent directory will be removed from the path and replaced with a prefix: "${key}"
        /// For example if a parent directory with key "AnalysisFolder" and value "/path/to/AnalysisFolder" exists, then the path: "/path/to/AnalysisFolder/myfile"
        /// will get serialized as "${AnalysisFolder}/myfile</param>
        public FileSystemLocationConverter(Dictionary<string, IDirectoryLocation> parentDirectories = null,
            IAbsolutePathInstantiator absPath = null)
        {
            _parentDirectories = parentDirectories ?? new Dictionary<string, IDirectoryLocation>();
            _absolutePathInstantiator = absPath ?? new RealFileAbsolutePathInstantiator();
        }

        public override bool CanConvert(Type objectType)
        {
            return Types.Any(t => t.IsAssignableFrom(objectType));
        }

        public override object ReadJson(JsonReader reader, Type objectType, object existingValue, JsonSerializer serializer)
        {
            // These are serialized as (string) FullName - create new objects from these strings.
            var path = (string)serializer.Deserialize(reader, typeof(string));
            if (path == null)
                return null;
            return ConvertRelativeToAbsolute(path, objectType);
        }

        public override void WriteJson(JsonWriter writer, object value, JsonSerializer serializer)
        {
            // Serialize FileLocation and DirectoryLocation as their FullName (string)
            string fileSysPath;
            if (typeof(IFileLocation).IsAssignableFrom(value.GetType()))
                fileSysPath = ((IFileLocation)value).FullName;
            else if (typeof(IDirectoryLocation).IsAssignableFrom(value.GetType()))
                fileSysPath = ((IDirectoryLocation)value).FullName;
            else
                throw new ApplicationException("Tried to deserialize something that wasn't IFileLocation or IDirectoryLocation.");
            serializer.Serialize(writer, ConvertAbsoluteToRelative(fileSysPath));
        }

        private object ConvertRelativeToAbsolute(string path, Type objectType)
        {
            if (Path.IsPathRooted(path))
                return _absolutePathInstantiator.InstantiateAbsolutePath(path, objectType);
            foreach (var kvp in _parentDirectories)
            {
                string key = kvp.Key;
                IDirectoryLocation parent = kvp.Value;
                string pathPrefix = GetPathPrefix(key);
                if (path.StartsWith(pathPrefix))
                {
                    string relativePath = path.Substring(pathPrefix.Length).Trim('\\', '/');
                    if (typeof(IFileLocation).IsAssignableFrom(objectType))
                    {
                        return parent.GetFileLocation(relativePath);
                    }
                    if (typeof(IDirectoryLocation).IsAssignableFrom(objectType))
                    {
                        return parent.GetDirectoryLocation(relativePath);
                    }
                    throw new ApplicationException("Tried to deserialize something that wasn't IFileLocation or IDirectoryLocation.");
                }
            }
            throw new ApplicationException($"Could not find parent directory to create relative path {path}.");
        }

        private static string GetPathPrefix(string key)
        {
            return "${" + key + "}";
        }

        private string ConvertAbsoluteToRelative(string path)
        {
            foreach (var kvp in _parentDirectories)
            {
                string key = kvp.Key;
                IDirectoryLocation parent = kvp.Value;
                string parentPath = parent.FullName;
                if (path.StartsWith(parentPath))
                {
                    string relativePath = path.Substring(parentPath.Length).Trim('\\', '/');
                    return GetPathPrefix(key) + Path.DirectorySeparatorChar + relativePath;
                }
            }
            return path;
        }
    }

    public class EnumerableConverter : JsonConverter
    {
        public override bool CanConvert(Type objectType)
        {
            return typeof(IEnumerable<dynamic>).IsAssignableFrom(objectType)
                && !typeof(List<dynamic>).IsAssignableFrom(objectType);
        }

        public override object ReadJson(JsonReader reader, Type objectType, object existingValue, JsonSerializer serializer)
        {
            throw new NotImplementedException();
        }

        public override void WriteJson(JsonWriter writer, object value, JsonSerializer serializer)
        {
            var enumer = value as IEnumerable<dynamic>;
            serializer.Serialize(writer, enumer.ToList());
        }

        public override bool CanRead => false;
    }
}
