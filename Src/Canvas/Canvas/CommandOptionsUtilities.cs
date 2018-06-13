using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;

namespace Canvas
{
    public static class CommandOptionsUtilities
    {
        /// <summary>
        ///     Takes an input command line and merges in arbitrary options. Returns the updated command line
        ///<para>- If the option already exists in the command line, then the value is overridden.</para>
        ///<para>- If the option does not exist in the command line, then the option and value are added after any other updated options</para>
        ///<para>- If there are no other updated options then we insert at end or beginning which is controlled by insertAtEnd parameter</para>
        ///<para>- Additional option #foo means: remove option foo if it's in the command-line now, otherwise do nothing.</para>
        ///<para>  (Remove "-foo" or "--foo" or "-foo bar" or "--foo bar")</para>
        ///<para>"-param -4","--param -4" "--param=-4" and "-p-4" are all supported. </para>
        /// </summary>
        public static string MergeCommandLineOptions(string commandLineWithOptions, string moreOptions, bool insertAtEnd = false)
        {
            if (string.IsNullOrEmpty(moreOptions))
                return commandLineWithOptions;

            //parse the existing command line and the new options to generate list of option key/value pairs for each
            string beforeFirstOption;
            List<KeyValuePair<string, string>> updatedOptions = GetCommandOptions(commandLineWithOptions, out beforeFirstOption);
            string beforeNewOptions;
            List<KeyValuePair<string, string>> newOptions = GetCommandOptions(moreOptions, out beforeNewOptions);

            //validate the new options
            if (!string.IsNullOrEmpty(beforeNewOptions.Trim()))
                throw new ArgumentException(string.Format("Unknown options format {0}", moreOptions));

            //special case, the last option value may include additional arguments to the command that are not part of that option
            //if we are updating this last option we can try to infer the number of white space delimited values the option has (usually 0 or 1) from the new values
            string afterLastOption = "";
            foreach (KeyValuePair<string, string> option in newOptions)
            {
                if (updatedOptions.Any() && option.Key == updatedOptions.Last().Key)
                {
                    string values = GetValues(updatedOptions.Last().Value, option.Value, out afterLastOption);
                    updatedOptions[updatedOptions.Count - 1] = new KeyValuePair<string, string>(option.Key, values);
                    break;
                }
            }

            //another special case, if we find a option to update make sure we do that first so we have a point of reference for inserting new options
            //inserting new options after this reference point has a better chance of being correct than arbitrarily choosing the beginning or the end of the command for insertion
            List<KeyValuePair<string, string>> prioritizedNewOptions = new List<KeyValuePair<string, string>>(newOptions);
            foreach (KeyValuePair<string, string> option in newOptions)
            {
                if (updatedOptions.Exists(kvp => option.Key == kvp.Key))
                {
                    prioritizedNewOptions.Remove(option);
                    prioritizedNewOptions.Insert(0, option);
                    break;
                }
            }

            //now we finally do the updating
            UpdateCommandOptions(updatedOptions, prioritizedNewOptions, insertAtEnd);

            //build up the final command that includes the update options
            StringBuilder updatedCommand = new StringBuilder(beforeFirstOption);
            foreach (KeyValuePair<string, string> option in updatedOptions)
            {
                updatedCommand.AppendWithSpace(option.Key);
                updatedCommand.Append(option.Value);
            }
            updatedCommand.AppendWithSpace(afterLastOption);
            return updatedCommand.ToString();
        }

        public static void UpdateCommandOptions(List<KeyValuePair<string, string>> originalOptions, List<KeyValuePair<string, string>> newOptions, bool insertAtEnd)
        {
            int lastUpdatedOptionIndex = -1;
            foreach (KeyValuePair<string, string> newOption in newOptions)
            {
                // Handle option removals:
                if (newOption.Key.StartsWith("#"))
                {
                    List<KeyValuePair<string, string>> removals = new List<KeyValuePair<string, string>>();
                    for (int index = 0; index < originalOptions.Count; index++)
                    {
                        var option = originalOptions[index];
                        if (option.Key.TrimStart('-') == newOption.Key.Substring(1))
                        {
                            removals.Add(option);
                            if (lastUpdatedOptionIndex >= index) lastUpdatedOptionIndex--;
                        }
                    }
                    foreach (var removal in removals) originalOptions.Remove(removal);
                    continue;
                }

                // Handle option overrides and insertions:
                int newOptionIndex = originalOptions.FindIndex(kvp => kvp.Key == newOption.Key);
                if (newOptionIndex != -1)
                {
                    originalOptions[newOptionIndex] = newOption;
                }
                else if (lastUpdatedOptionIndex != -1)
                    originalOptions.Insert(lastUpdatedOptionIndex + 1, newOption);
                else if (insertAtEnd)
                    originalOptions.Add(newOption);
                else
                    originalOptions.Insert(0, newOption);

                lastUpdatedOptionIndex = originalOptions.IndexOf(newOption);
            }
        }

        //count the number of values in the expectedValues string and grab the corresponding number of values from the currentValues string.
        //anything left over goes in the remainingValues string
        public static string GetValues(string currentValues, string expectedValues, out string remainingValues)
        {
            //a value can start at the beginning of the string with an optional = character and space --> \A=?\s*)
            //or simply by having space --> \s
            //the actual value is:
            //       a double quoted string containing any non double quote characters or any escaped double quote characters --> [""]([^""]|(?<=\\)[""])+[""]
            //     or
            //       a single quoted string containing any non single quote characters or any escaped single quote characters --> [']([^']|(?<=\\)['])+[']
            //     or
            //       a sequence of one or more non-whitespace characters --> \S+
            //finally the value must be followed by more space or the end of the string --> \s|\z
            Regex optionRegex = new Regex(@"((\A=?\s*|\s)\s*)([""]([^""]|(?<=\\)[""])+[""]|[']([^']|(?<=\\)['])+[']|\S+)(?=\s|\z)");
            int expectedValuesCount = optionRegex.Matches(expectedValues).Count;
            MatchCollection matches = optionRegex.Matches(currentValues);
            Match lastValue = matches[expectedValuesCount - 1];
            int indexRemainingValues = lastValue.Index + lastValue.Length;
            string values = currentValues.Substring(0, indexRemainingValues);
            remainingValues = currentValues.Substring(indexRemainingValues);
            return values;
        }

        //split the command into option key/value pairs and separately save anything that comes before the first option 
        public static List<KeyValuePair<string, string>> GetCommandOptions(string command, out string beforeFirstOption)
        {
            //this regex is the meat of this feature. We want to capture all valid options from the command
            //each option must happen after white space or at the beginning of the command
            //the option must start with a hyphen and then may contain more hyphens only if they are followed by a non-digit
            //the end of the option is identified by an = character or whitespace or a hyphen followed by a digit or the end of the string 
            Regex optionRegex = new Regex(@"(?<=\A|\s)(((--?)|(#))[a-zA-Z](?:\-(?=\D)|[_a-zA-Z])*)(?=[=\s]|\-?\d|\z)");
            List<KeyValuePair<string, string>> options = new List<KeyValuePair<string, string>>();
            Match m = optionRegex.Match(command);
            if (m.Success)
                beforeFirstOption = command.Substring(0, m.Index);
            else
                beforeFirstOption = command;
            while (m.Success)
            {
                int valueStartIndex = m.Index + m.Length;
                Match next = m.NextMatch();
                string value;
                if (next.Success)
                    value = command.Substring(valueStartIndex, next.Index - valueStartIndex);
                else
                    value = command.Substring(valueStartIndex);
                options.Add(new KeyValuePair<string, string>(m.Value, value));
                m = m.NextMatch();
            }
            return options;
        }

        //add white space if necessary
        private static void AppendWithSpace(this StringBuilder start, string end)
        {
            if (start.Length > 0 && !Char.IsWhiteSpace(start[start.Length - 1]) && !string.IsNullOrEmpty(end) && !Char.IsWhiteSpace(end[0]))
                start.Append(" ");
            start.Append(end);
        }
    }
}
