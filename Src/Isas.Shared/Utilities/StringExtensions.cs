using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace Isas.Shared.Utilities
{
    public static class StringExtensions
    {
        private static readonly Regex Whitespace = new Regex(@"\s+");

        public static bool IsInt(this string s)
        {
            int val;
            return int.TryParse(s, out val);
        }

        public static bool IsNullOrEmpty(this string value)
        {
            return string.IsNullOrEmpty(value);
        }

        public static bool IsNullOrWhiteSpace(this string value)
        {
            return string.IsNullOrWhiteSpace(value);
        }

        public static string RemoveWhiteSpace(this string value)
        {
            return Whitespace.Replace(value, "");
        }

        /// <summary>
        /// unlike string.Replace only one match if found will be replaced.
        /// match must occur at the start of the string
        /// </summary>
        /// <param name="value"></param>
        /// <param name="oldPrefix"></param>
        /// <param name="newPrefix"></param>
        /// <returns></returns>
        public static string ReplaceStart(this string value, string oldPrefix, string newPrefix)
        {
            if (!value.StartsWith(oldPrefix)) return value;
            return newPrefix + value.Substring(oldPrefix.Length);
        }

        /// <summary>
        /// unlike string.Replace only one match if found will be replaced.
        /// match must occur at the end of the string
        /// </summary>
        /// <param name="value"></param>
        /// <param name="oldSuffix"></param>
        /// <param name="newSuffix"></param>
        /// <returns></returns>
        public static string ReplaceEnd(this string value, string oldSuffix, string newSuffix)
        {
            if (!value.EndsWith(oldSuffix)) return value;
            return value.Substring(0, value.Length - oldSuffix.Length) + newSuffix;
        }

        /// <summary>
        /// returns the minimum string length from a sequence
        /// null strings are treated as zero length strings
        /// if the sequence is empty returns 0
        /// </summary>
        /// <param name="strings"></param>
        /// <returns></returns>
        public static int MinLength(this IEnumerable<string> strings)
        {
            if (!strings.Any()) return 0;
            if (strings.Any(item => item == null)) return 0;
            return strings.Select(item => item.Length).Min();
        }

        public static string[] Split(this string value, string separator)
        {
            return value.Split(new[] { separator }, StringSplitOptions.None);
        }
    }
}
