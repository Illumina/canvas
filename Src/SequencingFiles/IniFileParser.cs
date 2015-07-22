using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;

namespace SequencingFiles
{

	public class IniFileParser
	{
		private Dictionary<string, Dictionary<string, string>> sectionKeyPairs = new Dictionary<string, Dictionary<string, string>>(StringComparer.InvariantCultureIgnoreCase);
		private Dictionary<string, string> rootKeyPairs = new Dictionary<string, string>(StringComparer.InvariantCultureIgnoreCase);
		private readonly string _iniFilePath;
		private static List<string> _skipLineStart = new List<string> { "#" };
		/// <summary>
		/// Opens the INI file at the given path if it exists and enumerates the values in the IniParser.
		/// white space lines and lines starting with '#' are skipped
		/// </summary>
		/// <param name="iniPath">Full path to INI file.</param>
		public IniFileParser(string iniPath)
		{
			_iniFilePath = iniPath;
			if (!File.Exists(iniPath))
				return;

			using (TextReader iniFile = new StreamReader(_iniFilePath))
			{
				//start with root section
				string currentSectionName = null;
				while (true)
				{
					string strLine = iniFile.ReadLine();
					if (strLine == null) break;

					try
					{
						strLine = strLine.Trim();
						if (string.IsNullOrEmpty(strLine)) continue;
						if (_skipLineStart.Where(s => strLine.StartsWith(s)).Any()) continue;
						if (strLine.StartsWith("[") && strLine.EndsWith("]"))
						{
							currentSectionName = strLine.Substring(1, strLine.Length - 2);
							AddSection(ref currentSectionName);
							continue;
						}
						// return everything after the first '=' as one string
						string[] keyPair = strLine.Split(new char[] { '=' }, 2);
						keyPair = Array.ConvertAll(keyPair, p => p.Trim());

						string settingName = keyPair[0];
						string settingValue = null;
						if (keyPair.Length > 1)
							settingValue = keyPair[1];
						ValidateSettingName(ref settingName);
						ValidateSettingValue(ref settingValue);

						// if we have a duplicate setting it better have the same value!
						Dictionary<string, string> section;
						if (currentSectionName == null)
							section = rootKeyPairs;
						else
							section = sectionKeyPairs[currentSectionName];
						if (section.ContainsKey(settingName))
						{
							string oldSettingValue = section[settingName];
							if (settingValue != oldSettingValue)
								throw new ArgumentException(string.Format("Setting '{0}' already exists. Old value: '{1}'  New value: '{2}'", settingName, oldSettingValue ?? "null", settingValue ?? "null"));
						}

						// finally, just add the setting
						if (currentSectionName == null)
							AddRootSetting(settingName, settingValue);
						else
							AddSetting(currentSectionName, settingName, settingValue);
					}
					catch (ArgumentException e)
					{
						throw new ApplicationException(string.Format("Failed to parse ini file {0} in section '{1}'. Bad line: {2}", iniPath, currentSectionName ?? "ROOT", strLine), e);
					}
				}
			}
		}

		/// <summary>
		/// Returns the value for the given setting at the top level.
		/// </summary>
		/// <param name="settingName">Key name.</param>
		public List<string> GetRootSettings()
		{
			return rootKeyPairs.Keys.ToList();
		}

		/// <summary>
		/// Returns the value for the given setting at the top level.
		/// </summary>
		/// <param name="settingName">Key name.</param>
		public string GetRootSetting(string settingName)
		{
			ValidateSettingName(ref settingName);
			if (rootKeyPairs.ContainsKey(settingName))
				return rootKeyPairs[settingName];
			return null;
		}

		/// <summary>
		/// Returns the value for the given section, key pair.
		/// </summary>
		/// <param name="sectionName">Section name.</param>
		/// <param name="settingName">Key name.</param>
		public string GetSetting(string sectionName, string settingName)
		{
			ValidateSectionName(ref sectionName);
			ValidateSettingName(ref settingName);
			if (sectionKeyPairs.ContainsKey(sectionName) && sectionKeyPairs[sectionName].ContainsKey(settingName))
				return sectionKeyPairs[sectionName][settingName];
			return null;
		}

		/// <summary>
		/// Return a list of all the setting keys for a given section
		/// </summary>
		/// <param name="sectionName">Section to enum.</param>
		public List<string> GetSettings(string sectionName)
		{
			ValidateSectionName(ref sectionName);
			if (sectionKeyPairs.ContainsKey(sectionName))
				return sectionKeyPairs[sectionName].Keys.ToList();
			return null;
		}

		/// <summary>
		/// User cannot tell the difference between a setting that has a null value
		/// and a setting that is just missing. Using this they can
		/// </summary>
		/// <param name="sectionName">Section to enum.</param>
		public bool HasSetting(string sectionName, string settingName)
		{
			ValidateSectionName(ref sectionName);
			ValidateSettingName(ref settingName);
			if (sectionKeyPairs.ContainsKey(sectionName) && sectionKeyPairs[sectionName].ContainsKey(settingName))
				return true;
			return false;
		}

		/// <summary>
		/// User cannot tell the difference between a setting that has a null value
		/// and a setting that is just missing. Using this they can
		/// </summary>
		/// <param name="sectionName">Section to enum.</param>
		public bool HasRootSetting(string settingName)
		{
			ValidateSettingName(ref settingName);
			if (rootKeyPairs.ContainsKey(settingName))
				return true;
			return false;
		}

		/// <summary>
		/// Adds or replaces a setting to the table to be saved.
		/// </summary>
		/// <param name="sectionName">Section to add under.</param>
		/// <param name="settingName">Key name to add.</param>
		/// <param name="settingValue">Value of key.</param>
		public void AddSetting(string sectionName, string settingName, string settingValue)
		{
			AddSection(ref sectionName);
			ValidateSettingName(ref settingName);
			sectionKeyPairs[sectionName][settingName] = settingValue;
		}

		private void ValidateSectionName(ref string sectionName)
		{
			if (string.IsNullOrWhiteSpace(sectionName))
				throw new ArgumentException("Error: ini file section name cannot be null or white space");

			sectionName = sectionName.Trim();
		}

		private void ValidateSettingName(ref string settingName)
		{
			if (string.IsNullOrWhiteSpace(settingName))
				throw new ArgumentException("Error: ini file setting name cannot be null or white space");

			settingName = settingName.Trim();

			if (settingName.Contains('='))
				throw new ArgumentException(string.Format("Error: setting name {0} cannot contain '='", settingName));

			foreach (string skipStart in _skipLineStart)
			{
				if (settingName.StartsWith(skipStart))
					throw new ArgumentException(string.Format("Error: illegal setting name {0} for ini file. Setting name cannot start with '{0}' as these lines will be skipped when parsing.", settingName, skipStart));
			}
		}

		private void ValidateSettingValue(ref string value)
		{
			if (!string.IsNullOrEmpty(value))
				value = value.Trim();
		}

		/// <summary>
		/// Sometimes we want to add sections without any settings
		/// </summary>
		/// <param name="sectionName">Section to add.</param>
		public void AddSection(string sectionName)
		{
			AddSection(ref sectionName);
		}

		private void AddSection(ref string sectionName)
		{
			ValidateSectionName(ref sectionName);
			if (!sectionKeyPairs.ContainsKey(sectionName))
				sectionKeyPairs[sectionName] = new Dictionary<string, string>(StringComparer.InvariantCultureIgnoreCase);
		}

		/// <summary>
		/// Sometimes we want to remove entire sections
		/// </summary>
		/// <param name="sectionName">Section to add.</param>
		public void DeleteSection(string sectionName)
		{
			ValidateSectionName(ref sectionName);
			if (sectionKeyPairs.ContainsKey(sectionName))
				sectionKeyPairs.Remove(sectionName);
		}

		/// <summary>
		/// Removes all top level settings
		/// all settings under sections are retained
		/// </summary>
		public void DeleteRootSettings()
		{
			rootKeyPairs.Clear();
		}

		/// <summary>
		/// Utility function to get all section starting with a given string
		/// </summary>
		/// <param name="sectionName">Section to add.</param>
		public List<string> GetSectionsStartingWith(string start)
		{
			return sectionKeyPairs.Where(kvp => kvp.Key.StartsWith(start)).Select(kvp => kvp.Key).ToList();
		}

		/// <summary>
		/// Adds or replaces a setting to the table to be saved. 
		/// Setting is at the top level, not associated with a section
		/// </summary>
		/// <param name="settingName">Key name to add.</param>
		/// <param name="settingValue">Value to add.</param>
		public void AddRootSetting(string settingName, string settingValue = null)
		{
			ValidateSettingName(ref settingName);
			ValidateSettingValue(ref settingValue);
			rootKeyPairs[settingName] = settingValue;
		}

		/// <summary>
		/// Deletes a setting in the table 
		/// Setting is at the top level, not associated with a section
		/// </summary>
		/// <param name="settingName">Key name to delete.</param>
		public void DeleteRootSetting(string settingName)
		{
			ValidateSettingName(ref settingName);
			rootKeyPairs.Remove(settingName);
		}

		/// <summary>
		/// Remove a setting.
		/// </summary>
		/// <param name="sectionName">Section containing setting to remove.</param>
		/// <param name="settingName">Setting to remove.</param>
		public void DeleteSetting(string sectionName, string settingName)
		{
			ValidateSectionName(ref sectionName);
			ValidateSettingName(ref settingName);
			if (sectionKeyPairs.ContainsKey(sectionName))
				sectionKeyPairs[sectionName].Remove(settingName);
		}

		/// <summary>
		/// Save settings to new file.
		/// </summary>
		/// <param name="newFilePath">New file path.</param>
		public void SaveSettings(string newFilePath)
		{
			Directory.CreateDirectory(Path.GetDirectoryName(newFilePath));
			using (StreamWriter writer = new StreamWriter(newFilePath))
			{
				//write out root settings
				foreach (KeyValuePair<string, string> setting in rootKeyPairs)
				{
					WriteSettingLine(writer, setting);
				}

				//write out section settings
				foreach (KeyValuePair<string, Dictionary<string, string>> section in sectionKeyPairs)
				{
					writer.WriteLine("[{0}]", section.Key);
					foreach (KeyValuePair<string, string> setting in section.Value)
					{
						WriteSettingLine(writer, setting);
					}
				}
			}
		}

		private void WriteSettingLine(StreamWriter writer, KeyValuePair<string, string> setting)
		{
			if (setting.Value == null)
				writer.WriteLine("{0} =", setting.Key);
			else
				writer.WriteLine("{0} = {1}", setting.Key, setting.Value);
		}

		/// <summary>
		/// Save settings back to ini file.
		/// </summary>
		public void SaveSettings()
		{
			SaveSettings(_iniFilePath);
		}

		/// <summary>
		/// Return a list of all the section keys
		/// </summary>
		public List<string> GetSections()
		{
			return sectionKeyPairs.Keys.ToList();
		}

		/// <summary>
		/// are there any settings to write out or would we just write out an empty ini file?
		/// </summary>
		public bool Empty()
		{
			return true;
		}
	}
}
