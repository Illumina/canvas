using System;
using System.ComponentModel;
using System.Globalization;
using System.Text.RegularExpressions;

namespace Illumina.SecondaryAnalysis
{
	/// <summary>
	/// Provides Id and Name
	/// Workflow/Wrapper input/output classes will typically want to derive from SampleInfo
	/// </summary>
	[TypeConverter(typeof(SampleInfoConverter))]
	public class SampleInfo
	{
		public readonly string Id;
		public readonly string Name;
		public readonly int? Number;

		public SampleInfo(string id, string name, int? number = null)
		{
			Regex invalidCharacters = new Regex(@"[^A-Za-z0-9_-]");
			// Match anything that's not a letter, number, dash or underscore.  (Also match weird unicode characters that are still considered "letters")
			if (invalidCharacters.IsMatch(id))
				throw new ApplicationException(string.Format("Illegal characters in sample ID '{0}' - the ID must consist of letters, numbers, dashes and underscores.", id));
			if (invalidCharacters.IsMatch(name))
				throw new ApplicationException(string.Format("Illegal characters in sample name '{0}' - the name must consist of letters, numbers, dashes and underscores.", name));

			Id = id;
			Name = name;
			Number = number;
		}

		public override bool Equals(object obj)
		{
			SampleInfo sampleInfo = obj as SampleInfo;
			if (sampleInfo == null) return false;
			return string.Equals(Id, sampleInfo.Id, StringComparison.OrdinalIgnoreCase);
		}

		public override int GetHashCode()
		{
			return Id.ToLowerInvariant().GetHashCode();
		}

		public override string ToString()
		{
			return string.Format(@"(Id: {0}, Name: {1}, Number: {2})", Id, Name, Number.HasValue ? Number.Value.ToString() : "null");
		}
	}

	public class SampleInfoConverter : TypeConverter
	{
		public override bool CanConvertFrom(ITypeDescriptorContext context, Type source)
		{
			if (source == typeof(string))
				return true;
			return base.CanConvertFrom(context, source);
		}

		public override object ConvertFrom(ITypeDescriptorContext context, CultureInfo culture, object value)
		{
			string str = value as string;
			if (str == null)
				return base.ConvertFrom(context, culture, value);


			// "( id, name, number )"
			Regex regex = new Regex(@"\(\s*Id:(?<id>.*),\s*Name:(?<name>.*),\s*Number:(?<number>.*)\)");
			Match match = regex.Match(str);
			if (String.IsNullOrEmpty(match.Value))
				throw new ApplicationException(String.Format(@"Unrecognized SampleInfo string: ""{0}"". Expected format is ""(Id: SampleID, Name: SampleName, Number: SampleNumberOrNull)""", str));

			GroupCollection groups = match.Groups;
			string id = groups["id"].Value.Trim();
			string name = groups["name"].Value.Trim();
			string numberStr = groups["number"].Value.Trim();
			if (numberStr.Equals("null", StringComparison.OrdinalIgnoreCase))
				return new SampleInfo(id, name);

			int number;
			if (!Int32.TryParse(numberStr, out number))
				throw new ApplicationException(String.Format("Could not parse sample number from SampleInfo string: {0}. Value must be an integer or null", str));

			return new SampleInfo(id, name, number);
		}
	}
}
