using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Isas.Shared;
using SequencingFiles;

namespace SampleSettingsProcessing
{
    public static class SampleSettingsExtensions
    {
        public static string GetStringSetting(this ISampleSettings processor, string setting, string _default)
        {
            return processor.GetSetting(setting) ?? _default;
        }

        public static IFileLocation GetFileLocationSetting(this ISampleSettings processor, string setting, IFileLocation _default)
        {
            string value = processor.GetStringSetting(setting, null);
            return value == null ? _default : new FileLocation(value);
        }

        public static bool HasSetting(this ISampleSettings processor, string key)
        {
            return processor.GetSetting(key) != null;
        }

        public static int GetIntegerSetting(this ISampleSettings processor, string setting, int _default)
        {
            string value = processor.GetSetting(setting);
            if (value != null)
            {
                return int.Parse(value);
            }
            else
            {
                return _default;
            }
        }

        public static int? GetPositiveNullableIntegerSetting(this ISampleSettings processor, string setting, int? _default)
        {
            int? value = processor.GetNullableIntegerSetting(setting, _default);
            if (value < 0)
            {
                throw new ApplicationException($"{setting} cannot be negative. Use default: {_default}");
            }
            return value;
        }

        /// <summary>
        /// Parse an optional integer setting from the Samplesheet
        /// </summary>
        /// <param name="processor">Samplesheet info</param>
        /// <param name="setting">Name of the setting</param>
        /// <returns>Integer value of setting or null if not present</returns>
        public static int? GetNullableIntegerSetting(this ISampleSettings processor, string setting, int? _default = null)
        {
            string value = processor.GetSetting(setting);
            if (value != null)
                return int.Parse(value);
            else
                return _default;
        }

        public static float GetFloatSetting(this ISampleSettings processor, string setting, float _default)
        {
            string value = processor.GetSetting(setting);
            if (value != null)
            {
                return float.Parse(value);
            }
            else
            {
                return _default;
            }
        }

        /// <summary>
        /// Parse an optional float setting from the Samplesheet
        /// </summary>
        /// <param name="processor">Samplesheet info</param>
        /// <param name="setting">Name of the setting</param>
        /// <returns>Integer value of setting or null if not present</returns>
        public static float? GetNullableFloatSetting(this ISampleSettings processor, string setting)
        {
            string value = processor.GetSetting(setting);
            if (value != null)
            {
                return float.Parse(value);
            }
            else
            {
                return null;
            }
        }

        public static bool? ParseBooleanSetting(this string value)
        {
            var trueStrings = new HashSet<string>(StringComparer.OrdinalIgnoreCase) { "true", "t", "yes", "y", "1" };
            var falseStrings = new HashSet<string>(StringComparer.OrdinalIgnoreCase) { "false", "f", "no", "n", "0" };
            bool isTrue = trueStrings.Contains(value);
            bool isFalse = falseStrings.Contains(value);
            // If it's not set to anything we recognize, return null
            if ((!isTrue && !isFalse) || (isTrue && isFalse))
                return null;
            else
                return isTrue;
        }

        public static bool GetBooleanSetting(this ISampleSettings processor, string setting, bool _default)
        {
            string value = processor.GetSetting(setting);
            if (value == null) return _default;
            bool? boolvalue = ParseBooleanSetting(value);

            // If it's not set to anything we recognize, report an error!
            if (!boolvalue.HasValue)
            {
                throw new Exception($"Error: Sample sheet setting '{setting} = {value}' not understood; use true or false (or 1 or 0) for this setting");
            }
            return boolvalue.Value;
        }

        public static bool? GetNullableBooleanSetting(this ISampleSettings processor, string setting, bool? _default)
        {
            return setting.ParseBooleanSetting() ?? _default;
        }

        /// <summary>
        /// Get the data column corresponding to the list of keys.
        /// Checks each key in the list until one is present in the sample data section.
        /// </summary>
        /// <param name="keys">Possible names of the desired data column</param>
        /// <returns>SampleSet<List<string>>, or SampleSet<List<(string)null>> if key not present.</returns>
        public static SampleSet<List<string>> GetDataColumn(this ISampleSettings processor, params string[] keys)
        {
            SampleSet<List<string>> value = null;
            foreach (string key in keys)
            {
                value = processor.GetDataColumn(key);
                // Check that the "strings" being returned aren't null
                //  which means that the column doesn't exist.
                if (value.First().Value.First() != null) return value;
            }
            return value;
        }

        /// <summary>
        /// Like GetDataColumn, but returns null if none of the keys are found.
        /// </summary>
        /// <param name="processor">The instance of the settings processor this method extends</param>
        /// <param name="keys">Possible names of the desired data column</param>
        /// <returns></returns>
        public static SampleSet<List<string>> GetDataColumnOrNull(this ISampleSettings processor, params string[] keys)
        {
            SampleSet<List<string>> value = null;
            foreach (string key in keys)
            {
                value = processor.GetDataColumn(key);
                // Check that the "strings" being returned aren't null
                //  which means that the column doesn't exist.
                if (value.First().Value.First() != null) return value;
            }
            return null;
        }

        /// <summary>
        /// Get the data column and check it for consistency across the sample IDs.
        /// Compares the string values for all entries in the column belonging to the same
        /// sample ID. If these values are not identical, throws an exception.
        /// Returns a SampleSet of string, as opposed to a SampleSet of List of string.
        /// </summary>
        /// <param name="processor"></param>
        /// <param name="requireData">If true, also ensures that all entries are not empty.</param>
        /// <param name="keys">Possible names of the column.</param>
        /// <returns>A SampleSet of string</returns>
        public static SampleSet<string> GetConsistentDataColumn(this ISampleSettings processor, bool requireData,
            params string[] keys)
        {
            SampleSet<List<string>> column = processor.GetDataColumn(keys);
            if (column.First().Value.First() == null) return column.SelectData(s => (string)null);

            string columnName = string.Join(" | ", keys);
            foreach (var kvp in column)
            {
                // Check that all values for this sample are the same
                if (kvp.Value.Distinct().Count() > 1)
                    throw new ApplicationException(string.Format("Sample {0} has inconsistent values in column {1}",
                        kvp.Key.Id, columnName));

                if (requireData && string.IsNullOrWhiteSpace(kvp.Value.First()))
                    throw new ApplicationException(string.Format("Sample {0} is missing a value for column {1}",
                        kvp.Key.Id, columnName));
            }
            return column.SelectData(s => s.First());
        }

        /// <summary>
        /// Get a data column in the sample settings, check the values for consistency,
        /// and convert the values to booleans.
        /// Returns null if the column does not exist.
        /// </summary>
        /// <param name="processor"></param>
        /// <param name="keys">Possible names of the column</param>
        /// <returns>SampleSet of bool</returns>
        public static SampleSet<bool> GetBooleanDataColumn(this ISampleSettings processor, params string[] keys)
        {
            if (!processor.HasDataColumn(keys))
                return null;
            SampleSet<List<string>> column = processor.GetDataColumn(keys);
            SampleSet<bool> boolcolumn = new SampleSet<bool>();
            foreach (SampleInfo sample in column.Keys)
            {
                foreach (string value in column[sample])
                {
                    bool? boolvalue = ParseBooleanSetting(value);
                    if (!boolvalue.HasValue)
                    {
                        string keystring = string.Join(" OR ", keys);
                        throw new Exception(string.Format("Errror: Sample data column: '{0}', ID: '{1}', value: '{2}' not understood; use true or false (or 1 or 0) for this setting",
                            keystring, sample.Id, value));
                    }
                    else
                    {
                        // Make sure this value matches others for the same sample ID
                        if (boolcolumn.ContainsKey(sample) && boolcolumn[sample] != boolvalue.Value)
                        {
                            string keystring = string.Join(" OR ", keys);
                            throw new Exception(string.Format("Error: Inconsistent data for sample '{0}' in column '{1}'",
                                sample.Id, keystring));
                        }
                        boolcolumn.Add(sample, boolvalue.Value);
                    }
                }
            }
            return boolcolumn;
        }

        /// <summary>
        /// A convenience method to improve code readability
        /// </summary>
        public static bool HasDataColumn(this ISampleSettings processor, params string[] keys)
        {
            return processor.GetDataColumn(keys).First().Value.First() != null;
        }

        /// <summary>
        /// Returns a SampleSet of SampleInfo.
        /// Assumes that requesting the column "sample_id" will at least return a SampleSet.
        /// </summary>
        /// <param name="processor"></param>
        /// <returns>A SampleSet of the SampleInfo's in the settings.</returns>
        public static SampleSet<SampleInfo> GetSampleInfos(this ISampleSettings processor)
        {
            return processor.GetDataColumn("sample_id").SelectData((info, d) => info);
        }

        /// <summary>
        /// Returns a SampleSet of nulls.
        /// </summary>
        public static SampleSet<T> GetNullSet<T>(this ISampleSettings processor) where T : class
        {
            return processor.GetSampleInfos().SelectData<T>(info => null);
        }
    }

    public class SampleSettingsProcessor : ISampleSettings
    {
        // Sample columns with info about the experimental setup; Ignored by Isas.
        static List<string> _knownColumns = new List<string>()
        {
            "Sample_Plate", "Sample_Well", "Sample_Project", "I7_Index_ID", "I5_Index_ID",
            "Description"
        };
        public string SampleSheetPath { get; private set; }
        private ILogger logger;
        private readonly Func<string, string[]> ParseCommaDelimitedLine;

        private Dictionary<string, List<string>> Sections;
        private Dictionary<string, string> Header;
        private Dictionary<string, string> Settings;
        private Dictionary<string, SampleSet<List<string>>> Data;
        private Dictionary<string, string> Manifests;

        private ISet<string> _accessedSettings;
        private ISet<string> _accessedColumns;
        private ISet<string> _accessedManifests;

        public SampleSettingsProcessor(string sampleSheetPath, ILogger logger)
        {
            SampleSheetPath = sampleSheetPath;
            this.logger = logger;

            ParseCommaDelimitedLine = s => CSVReader.ParseCommaDelimitedLine(s).Select(bit => bit.Trim()).ToArray();

            Sections = null;

            try
            {
                using (var reader = new StreamReader(sampleSheetPath))
                {
                    Sections = ParseSections(reader);
                }
            }
            catch (Exception ex)
            {
                try
                {
                    // When Excel has a file open, it locks it, prohibiting even reading the file.
                    // In these cases, however, we can make a copy of the file and read the copy
                    // So, just in case Excel is being a control freak, we'll try the copy trick
                    // before giving up. 
                    logger.Error("{0} may be locked - attempting to read a copy", sampleSheetPath);
                    string newSampleSheetPath = Path.Combine(
                        IsasConfiguration.GetConfiguration().TempFolder, sampleSheetPath + ".readable_copy");
                    Utilities.SafeDelete(newSampleSheetPath);
                    File.Copy(sampleSheetPath, newSampleSheetPath);
                    using (var reader = new StreamReader(newSampleSheetPath))
                    {
                        Sections = ParseSections(reader);
                    }
                    File.Delete(newSampleSheetPath);
                }
                catch
                {
                    // Swallow the 2nd error and just report the original error
                    logger.Error("Unable to open sample sheet {0}\n{1}", sampleSheetPath, ex);
                    throw ex;
                }
            }

            _accessedSettings = new HashSet<string>();
            _accessedColumns = new HashSet<string>();
            _accessedManifests = new HashSet<string>();

            Header = ParseHeader(Sections);
            Settings = ParseSettings(Sections);
            Data = ParseData(Sections);
            Manifests = ParseManifests(Sections);
        }

        private Dictionary<string, string> ParseHeader(Dictionary<string, List<string>> Sections)
        {
            // Header section begins with [Header]
            Dictionary<string, string> header = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
            foreach (string line in Sections["Header"])
            {
                var tokens = ParseCommaDelimitedLine(line);
                header[tokens[0]] = tokens[1];
            }
            return header;
        }

        private Dictionary<string, string> ParseSettings(Dictionary<string, List<string>> Sections)
        {
            // Settings section begins with [Settings]
            Dictionary<string, string> settings = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
            if (!Sections.ContainsKey("Settings")) return settings;
            foreach (var line in Sections["Settings"])
            {
                var tokens = ParseCommaDelimitedLine(line);
                if (tokens[0] == "") continue;
                if (tokens.Length < 2)
                    throw new ArgumentException(string.Format("Missing value for sample sheet setting {0}", tokens[0]));
                settings[tokens[0]] = tokens[1];
            }
            return settings;
        }

        private Dictionary<string, SampleSet<List<string>>> ParseData(Dictionary<string, List<string>> Sections)
        {
            // Data section begins with [Data]
            Dictionary<string, SampleSet<List<string>>> data = new Dictionary<string, SampleSet<List<string>>>(StringComparer.OrdinalIgnoreCase);
            if (!Sections.ContainsKey("Data") || !Sections["Data"].Any())
                return data;

            ISet<SampleInfo> sampleInfos = new HashSet<SampleInfo>();

            // First row contains column names
            Dictionary<string, int> Columns = ParseColumnNames(Sections["Data"][0]);
            foreach (var key in Columns.Keys)
            {
                data[key] = new SampleSet<List<string>>();
            }

            int SampleNumber = 1;
            foreach (string line in Sections["Data"].Skip(1))
            {
                var tokens = ParseCommaDelimitedLine(line);
                if (tokens.All(string.IsNullOrWhiteSpace)) continue;

                // Parse the sample info
                SampleInfo info = ParseSampleInfo(tokens, Columns);
                if (sampleInfos.Contains(info))
                {
                    // We've already seen this sample ID, so get the existing entry (which has a sample number assigned)
                    SampleInfo tempInfo = sampleInfos.First(infoA => infoA.Equals(info));

                    // Make sure the ID case matches
                    if (info.Id != tempInfo.Id)
                        throw new FileLoadException(string.Format("Found Sample ID {0} and {1} that differ only in case", info.Id, tempInfo.Id));

                    // Make sure the SampleInfo.Name matches a previous SampleInfo
                    if (info.Name != tempInfo.Name)
                        throw new FileLoadException(string.Format("SampleName must be the same for all rows with SampleID {0}", info.Id));
                    info = tempInfo;
                }
                else
                {
                    // Haven't seen this sample ID before, so add it to the set and assign it a sample number
                    info = new SampleInfo(info.Id, info.Name, SampleNumber++);
                    sampleInfos.Add(info);
                }

                // Add the sample information for each column
                foreach (var kvp in Columns) // key is column name, value is position in row
                {
                    if (!data[kvp.Key].ContainsKey(info))
                    {
                        // Haven't seen this sample ID for this column yet - initialize the list
                        data[kvp.Key][info] = new List<string>();
                    }
                    if (tokens.Length <= kvp.Value)
                    {
                        // header has more columns than this data row. Assume missing columns are null 
                        data[kvp.Key][info].Add(null);
                    }
                    else
                    {
                        data[kvp.Key][info].Add(tokens[kvp.Value]);
                    }
                }
            }
            return data;
        }

        private Dictionary<string, int> ParseColumnNames(string line)
        {
            Dictionary<string, int> Columns = new Dictionary<string, int>(StringComparer.OrdinalIgnoreCase);
            var tokens = ParseCommaDelimitedLine(line);
            for (int ii = 0; ii < tokens.Count(); ii++)
            {
                var column = tokens[ii];
                // Only add the column to the dictionary if it has a value
                //  and is not already present in the dictionary.
                if (!String.IsNullOrWhiteSpace(column))
                {
                    if (Columns.ContainsKey(column))
                    {
                        throw new Exception(String.Format("Duplicate column headers found in Data section: {0}", column));
                    }
                    Columns[column] = ii;
                }
            }
            return Columns;
        }

        private SampleInfo ParseSampleInfo(string[] tokens, Dictionary<string, int> columns)
        {
            // Sample info needs ID, Name, and Number
            const string sampleId = "SampleId";
            const string sample_Id = "Sample_Id";
            string id = null;
            if (!columns.ContainsKey(sampleId) && !columns.ContainsKey(sample_Id))
                throw new FileLoadException(string.Format("{0} or {1} must be specified", sampleId, sample_Id));
            if (columns.ContainsKey(sampleId) && columns.ContainsKey(sample_Id))
                throw new FileLoadException(string.Format("Cannot specify both {0} and {1}", sampleId, sample_Id));
            if (columns.ContainsKey(sampleId))
            {
                id = tokens[columns[sampleId]];
                _accessedColumns.Add(sampleId);
            }
            if (columns.ContainsKey(sample_Id))
            {
                id = tokens[columns[sample_Id]];
                _accessedColumns.Add(sample_Id);
            }
            if (string.IsNullOrEmpty(id))
                throw new FileLoadException(string.Format("{0} or {1} must be specified", sampleId, sample_Id));

            const string sampleName = "SampleName";
            const string sample_Name = "Sample_Name";
            if (columns.ContainsKey(sampleName) && columns.ContainsKey(sample_Name))
                throw new FileLoadException(string.Format("Cannot specify both {0} and {1}", sampleName, sample_Name));
            string name = id;
            if (columns.ContainsKey(sampleName) && !string.IsNullOrEmpty(tokens[columns[sampleName]]))
            {
                name = tokens[columns[sampleName]];
                _accessedColumns.Add(sampleName);
            }
            if (columns.ContainsKey(sample_Name) && !string.IsNullOrEmpty(tokens[columns[sample_Name]]))
            {
                name = tokens[columns[sample_Name]];
                _accessedColumns.Add(sample_Name);
            }

            return new SampleInfo(id, name);
        }

        private Dictionary<string, string> ParseManifests(Dictionary<string, List<string>> Sections)
        {
            // Manifests section begins with [Manifests]
            Dictionary<string, string> manifests = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
            if (!Sections.ContainsKey("Manifests")) return manifests;
            foreach (var line in Sections["Manifests"])
            {
                var tokens = ParseCommaDelimitedLine(line);
                manifests[tokens[0]] = tokens[1];
            }
            return manifests;
        }

        private Dictionary<string, List<string>> ParseSections(StreamReader reader)
        {
            // Read through file
            // Section labels become the keys
            // Lines under that section are added to the list
            Dictionary<string, List<string>> Sections = new Dictionary<string, List<string>>(StringComparer.OrdinalIgnoreCase);

            string line;
            string sectionLabel = null;
            while ((line = reader.ReadLine()) != null)
            {
                string[] fields = ParseCommaDelimitedLine(line);
                if (fields.All(string.IsNullOrEmpty)) continue;
                if (line.StartsWith("["))
                {
                    // This is a section label
                    sectionLabel = ParseSectionName(fields[0]);
                    Sections[sectionLabel] = new List<string>();
                }
                else
                {
                    if (string.IsNullOrEmpty(sectionLabel))
                        throw new FileLoadException(string.Format("Found invalid line: {0}", line));
                    // Add this line to the current section
                    Sections[sectionLabel].Add(line);
                }
            }
            return Sections;
        }

        private string ParseSectionName(string sectionName)
        {
            // return string between [ ]
            int start = sectionName.IndexOf('[');
            if (start < 0)
                throw new ApplicationException($"{sectionName} is not a Section name.");
            int end = sectionName.IndexOf(']', start);
            if (end < 0)
                throw new ApplicationException($"{sectionName} is not a Section name.");
            return sectionName.Substring(start + 1, end - start - 1);
        }

        private TValue GetValueOrDefault<TKey, TValue>(IDictionary<TKey, TValue> dictionary, TKey key, TValue value)
        {
            if (dictionary.ContainsKey(key))
            {
                return dictionary[key];
            }
            else
            {
                return value;
            }
        }

        public string GetHeader(string key)
        {
            return GetValueOrDefault(Header, key, null);
        }

        public string GetSetting(string key)
        {
            _accessedSettings.Add(key);
            return GetValueOrDefault(Settings, key, null);
        }

        public SampleSet<List<string>> GetDataColumn(string key)
        {
            _accessedColumns.Add(key);
            // default return value is a null string for each row in the sample sheet
            var defaultSampleSet = Data.First().Value.SelectData(list => list.Select(s => (string)null).ToList());
            return GetValueOrDefault(Data, key, defaultSampleSet);
        }

        public string GetManifest(string key)
        {
            _accessedManifests.Add(key);
            return GetValueOrDefault(Manifests, key, null);
        }

        /// <summary>
        /// Warn about settings in the Samplesheet or Manifest that are not being used by Isas
        /// </summary>
        public void CheckUnusedEntries()
        {
            // Check Settings
            var unusedSettings = CompareSets(_accessedSettings, Settings.Keys);
            LogUnused(unusedSettings, "Settings");

            // Check Data Columns
            var unusedColumns = CompareSets(_accessedColumns.Union(_knownColumns), Data.Keys);
            LogUnused(unusedColumns, "Data Columns");

            // Check Manifests
            var unusedManifests = CompareSets(_accessedManifests, Manifests.Keys);
            LogUnused(unusedManifests, "Manifests");
        }

        // Return the keys that were provided by not accessed
        private IEnumerable<string> CompareSets(IEnumerable<string> accessed, IEnumerable<string> provided)
        {
            return provided.Except(accessed, StringComparer.OrdinalIgnoreCase);
        }

        private void LogUnused(IEnumerable<string> unused, string category)
        {
            if (unused.Count() == 0) return;

            logger.Warn("The following {0} were unused:", category);
            foreach (var entry in unused)
            {
                logger.Warn("\t{0}", entry);
            }
        }

        public SampleSet<SampleInfo> GetSamples()
        {
            var samples = new SampleSet<SampleInfo>();
            foreach (SampleInfo sample in Data.First().Value.SampleInfos)
            {
                samples[sample] = sample;
            }
            return samples;
        }

        public List<string> GetSection(string sectionName)
        {
            return GetValueOrDefault(Sections, sectionName, null);
        }
    }
}
