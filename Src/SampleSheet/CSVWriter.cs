using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace SequencingFiles
{
	/// <summary>
	///     A helper class to write csv files
	/// </summary>
	public static class CSVWriter
	{
		public static string GetLine(params string[] values)
		{
			return string.Join(",", values.Select(Escape));
		}

		public static string GetLine(List<string> values)
		{
			return GetLine(values.ToArray());
		}

		public static string Escape(string s)
		{
			if (s.Contains(QUOTE))
				s = s.Replace(QUOTE, ESCAPED_QUOTE);

			if (s.IndexOfAny(CHARACTERS_THAT_MUST_BE_QUOTED) > -1)
				s = QUOTE + s + QUOTE;

			return s;
		}

		private const string QUOTE = "\"";
		private const string ESCAPED_QUOTE = "\"\"";
		private static char[] CHARACTERS_THAT_MUST_BE_QUOTED = { ',', '"', '\n' };
	}
}
