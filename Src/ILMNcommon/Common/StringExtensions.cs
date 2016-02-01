using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;

namespace ILMNcommon.Common
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

        public static int CommonPrefixLength(this string a, string b)
        {
            int commonPrefixLength = 0;
            for (; commonPrefixLength < a.Length && commonPrefixLength < b.Length; commonPrefixLength++)
            {
                if (a[commonPrefixLength] != b[commonPrefixLength]) break;
            }
            return commonPrefixLength;
        }

        public static int CommonSuffixLength(this string a, string b)
        {
            int commonSuffixLength = 0;
            for (; commonSuffixLength < a.Length && commonSuffixLength < b.Length; commonSuffixLength++)
            {
                if (a[a.Length - 1 - commonSuffixLength] !=
                    b[b.Length - 1 - commonSuffixLength])
                    break;
            }
            return commonSuffixLength;
        }

        /// <summary>
        /// Find the shortest repeating substring in this string. For example, ATATAT --> AT and ATA --> ATA
        /// The first part of this algorithm is the table building part of the Knuth-Morris-Pratt algorithm which finds
        /// the widest border (same string at beginning and end) for each prefix. The length of the widest border for the entire string
        /// is then used to determine the length of the shortest repeat.
        /// </summary>
        public static string ShortestRepeatingSubstring(this string input)
        {
            var widestBorder = new int[input.Length];
            foreach (var i in Enumerable.Range(1, input.Length - 1))
            {
                var k = widestBorder[i - 1];
                while (true)
                {
                    if (input[i] == input[k])
                    {
                        widestBorder[i] = k + 1;
                        break;
                    }
                    else if (k == 0)
                    {
                        widestBorder[i] = 0;
                        break;
                    }
                    else
                    {
                        k = widestBorder[k - 1];
                    }
                }
            }
            var shortestRepeatLength = input.Length - widestBorder.Last();
            if (input.Length % shortestRepeatLength != 0)
                return input;

            return input.Substring(0, shortestRepeatLength);
        }
    }
}
