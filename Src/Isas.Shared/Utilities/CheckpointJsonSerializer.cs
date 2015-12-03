using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using Isas.Shared;
using Newtonsoft.Json;
using Newtonsoft.Json.Serialization;

namespace Illumina.SecondaryAnalysis.Workflow
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

        public CheckpointJsonSerializer(JsonSerializerSettings settings, IDirectoryLocation repository, ILogger logger)
        {
            _logger = logger;
            _repository = repository;
            if (!_repository.Exists)
                _repository.Create();
            _settings = settings;
        }

        public void Save<T>(string checkpoint, T input)
        {
            IFileLocation path = GetCurrentCheckpointFile(checkpoint);
            _logger.Info("Saving checkpoint results to {0}", path);
            Serialize(path.FullName, input);
        }

        public T Load<T>(string checkpoint)
        {
            IFileLocation path = GetCurrentCheckpointFile(checkpoint);
            if (!Exists(checkpoint))
                throw new ApplicationException(string.Format("Unable to load checkpoint results. Missing expected checkpoint file at {0}", path));
            _logger.Info("Loading checkpoint results from {0}", path);
            return Deserialize<T>(path.FullName);
        }

        public bool Exists(string checkpoint)
        {
            IFileLocation path = GetCurrentCheckpointFile(checkpoint);
            if (!path.Exists) return false;
            return true;
        }

        private void Serialize<T>(string path, T obj)
        {
            File.WriteAllText(path, JsonConvert.SerializeObject(obj, Formatting.Indented, _settings));
        }

        private T Deserialize<T>(string path)
        {
            return JsonConvert.DeserializeObject<T>(File.ReadAllText(path), _settings);
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

    public class FileSystemLocationConverter : JsonConverter
    {
        private static readonly Type[] Types = { typeof(IFileLocation), typeof(IDirectoryLocation) };
        private readonly Dictionary<string, IDirectoryLocation> _parentDirectories;

        /// <summary>
        /// Constructor for custom FileSystemLocationConverter
        /// </summary>
        /// <param name="parentDirectories">List of parent directories to use when serializing to store absolute paths as relative paths.
        /// If a location starts with one of the parent directories, the parent directory will be removed from the path and replaced with a prefix: "${key}"
        /// For example if a parent directory with key "AnalysisFolder" and value "/path/to/AnalysisFolder" exists, then the path: "/path/to/AnalysisFolder/myfile"
        /// will get serialized as "${AnalysisFolder}/myfile</param>
        public FileSystemLocationConverter(Dictionary<string, IDirectoryLocation> parentDirectories = null)
        {
            _parentDirectories = parentDirectories ?? new Dictionary<string, IDirectoryLocation>();
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
            path = ConvertRelativeToAbsolute(path);
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

        public override void WriteJson(JsonWriter writer, object value, JsonSerializer serializer)
        {
            // Serialize FileLocation and DirectoryLocation as their FullName (string)
            var fileSysObj = value as FileSystemLocationBase;
            serializer.Serialize(writer, ConvertAbsoluteToRelative(fileSysObj.FullName));
        }

        private string ConvertRelativeToAbsolute(string path)
        {
            if (Path.IsPathRooted(path))
                return path;
            foreach (var kvp in _parentDirectories)
            {
                string key = kvp.Key;
                IDirectoryLocation parent = kvp.Value;
                string pathPrefix = GetPathPrefix(key);
                if (path.StartsWith(pathPrefix))
                {
                    string relativePath = path.Substring(pathPrefix.Length).Trim('\\', '/');
                    return Path.GetFullPath(Path.Combine(parent.FullName, relativePath));
                }
            }
            throw new ApplicationException(string.Format("Could not find parent directory for relative path {0}."));
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
}
