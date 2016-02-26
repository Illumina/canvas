using Isas.Shared;
using System.Collections.Generic;

namespace SampleSettingsProcessing
{
    public interface ISampleSettings
    {
        // Someday we hope to get rid of this method, but we need it for now to make old and refactored workflows play nicely together
        string SampleSheetPath { get; }

        /// <summary>
        /// Returns the value associated with the given key in the Header section.
        /// </summary>
        /// <param name="key">The key to look up.</param>
        /// <returns>The value associated with the key, or null if key not present.</returns>
        string GetHeader(string key);

        /// <summary>
        /// Returns the value associated with the given key in the Settings section.
        /// </summary>
        /// <param name="key">The key to look up.</param>
        /// <returns>The value associated with the key, or null if key not present.</returns>
        string GetSetting(string key);

        /// <summary>
        /// 
        /// </summary>
        /// <param name="key"></param>
        /// <returns>SampleSet<List<string>>, or SampleSet<List<(string)null>> if key not present.</returns>
        SampleSet<List<string>> GetDataColumn(string key);

        SampleSet<SampleInfo> GetSamples();

        string GetManifest(string key);

        List<string> GetSection(string sectionName);

        void CheckUnusedEntries();

    }


}
