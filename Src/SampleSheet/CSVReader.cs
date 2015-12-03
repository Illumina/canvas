using System;
using System.Collections.Generic;
using System.IO;
using System.Reflection;

namespace SequencingFiles
{
	/// <summary>
	///     Utility class to handle the reading of Comma Separated Values (CSV) files.
	/// </summary>
	public class CSVReader
	{

		//Parse CSV files and return a list of lists
		static public List<List<string>> ParseCSVReportFile(string csvFile, bool hasHeader = false)
		{
			StreamReader reader = new StreamReader(File.OpenRead(csvFile));
			List<List<string>> csvData = new List<List<string>>();
			bool gotHeader = false;

			while (!reader.EndOfStream)
			{
				string line = reader.ReadLine();
				if (string.IsNullOrEmpty(line) || line[0].Equals('#')) continue; // skip blank line or comment
				if (hasHeader && !gotHeader)
				{
					gotHeader = true;
					continue;
				}
				string[] values = CSVReader.ParseCommaDelimitedLine(line);
				if (values.Length > 1) //ignore lines with no data or only one entry (headers, etc.)
				{
					if (values[0].EndsWith(":"))
					{
						values[0] = values[0].Substring(0, values[0].Length - 1);
					}
					csvData.Add(new List<string>(values));
				}
			}

			reader.Close();

			// Return the result
			return csvData;
		}

		/// <summary>
		///     Parse a csv file line.  Normally FileLine.split(',') gives us what we want, but we must handle
		///     embedded commas and double-quotes.  Details: http://creativyst.com/Doc/Articles/CSV/CSV01.htm
		/// </summary>
		public static string[] ParseCommaDelimitedLine(string fileLine)
		{
			List<string> bits = new List<string>();
			bool quoted = false;
			string currentBit = "";
			for (int position = 0; position < fileLine.Length; position++)
			{
				char character = fileLine[position];
				switch (character)
				{
					case '\r': // ignore newlines!
					case '\n':
						break;
					case '"':
						if (quoted)
						{
							quoted = false;
						}
						else
						{
							quoted = true;
							if (position > 0 && fileLine[position - 1] == '"')
							{
								currentBit += "\"";
							}
						}
						break;
					case ',':
						if (position > 0 && fileLine[position - 1] == ',')
						// look behind for another comma--if so, this is an empty field
						{
							bits.Add("");
							currentBit = "";
						}
						else if (!quoted)
						{
							bits.Add(currentBit.Trim()); // Trim to avoid white space
							currentBit = "";
						}
						else
						{
							currentBit += ",";
						}
						break;
					default:
						currentBit += character;
						break;
				}
			}
			bits.Add(currentBit.Trim());
			return bits.ToArray();
		}
	}
}