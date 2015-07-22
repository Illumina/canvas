using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Runtime.InteropServices;
using System.Text;

namespace Illumina.Common
{
    /// <summary>
    /// Misc functions for manipulating exceptions
    /// </summary>
    public class ExceptionUtils
    {
        public static string GetMessagesRecursiveAsSingleString(Exception ex, string separator = null)
        {
            separator = separator ?? (";" + Environment.NewLine);

            string[] errors = GetMessagesRecursive(ex, true);
            StringBuilder sb = new StringBuilder(10000);
            sb.AppendFormat(errors[0]);

            for (int i = 1; i < errors.Length; i++)
            {
                sb.Append(separator);
                sb.AppendFormat(errors[i]);
            }

            return sb.ToString();
        }

        /// <summary>
        /// Recurses all the inner messages of "ex", and returns the "Message" property
        /// of each as an element of the returned array. The "Message" property is enriched
        /// with additional information for certain exception types.
        /// The order of the returned elements is such that the 0 element is the message from the originating top-level 
        /// exception and 
        /// the Length-1 element is the lowest level enclosing exception (GetBaseException()).
        /// Note that is the same order returned by 
        /// GetExceptionTypesRecursive() and GetStackTraceRecursive().
        /// The number of entries may differ from the aforesaid functions if some of the messages 
        /// have embedded newlines, since the messages are split by newline.
        /// The only caveat here is if AggregateException objects encountered,
        /// in which case the order of the messages may be somewhat interleaved.
        /// </summary>
        public static string[] GetMessagesRecursive(Exception ex)
        {
            return GetMessagesRecursive(ex, false);
        }

        public static string[] GetMessagesRecursive(Exception exception, bool ignoreInternalExceptions)
        {
            List<string> msgs = new List<string>();

            IEnumerable<Exception> exceptions = GetExceptionsRecursive(exception);

            foreach (Exception ex in exceptions)
            {
                string msg = ex.Message;

#if SplitAtNewlines
                string[] lines = SplitAtNewlines(msg);
                foreach (string line in lines)
                    msgs.Add(line);
#else
                msg = RemoveNewlines(msg);

                if (!ignoreInternalExceptions)
                    msgs.Add(msg);
#endif

                if (ex is ArgumentException)
                {
                    string paramName = (ex as ArgumentException).ParamName;

                    if (!string.IsNullOrEmpty(paramName) || !ignoreInternalExceptions)
                        msgs.Add(String.Format("Invalid parameter name is \"{0}\"", paramName));
                }

                if (ex is ExternalException)
                {
                    int errorCode = ((ExternalException)(ex)).ErrorCode;

                    msgs.Add(String.Format("Error code is \"{0}\"", errorCode));
                }

                if (ex is ReflectionTypeLoadException)
                    msgs.AddRange((ex as ReflectionTypeLoadException).LoaderExceptions.Select(loadEx => RemoveNewlines(loadEx.Message)));
            }
            return msgs.ToArray();
        }

        public static void DebugPrintException(Exception ex)
        {
            string[] errorMessages = GetMessagesRecursive(ex);
            foreach (string t in errorMessages)
            {
                Console.WriteLine("* {0}", t);
            }
        }

        /// <summary>
        /// Returns a list of all the exceptions represented by the exception argument.  The list has the top-level 
        /// exception in position 0 (the exception argument to this function), and progressively does
        /// a depth-first traversal down into inner exceptions, adding to list as it goes.  Breadth-wise traversal occurs for 
        /// System.AggregateException; exceptions of those types
        /// are not included in the list, only their children.
        /// </summary>
        public static IEnumerable<Exception> GetExceptionsRecursive(Exception ex)
        {
            List<Exception> exceptions = new List<Exception>();
            GetExceptionsRecursiveImpl(ex, exceptions);
            return exceptions;
        }

        private static void GetExceptionsRecursiveImpl(Exception ex, List<Exception> exceptions)
        {
            if (ex is AggregateException)
            {
                AggregateException aggregate = ex as AggregateException;

                foreach (Exception thisEx in aggregate.InnerExceptions)
                    GetExceptionsRecursiveImpl(thisEx, exceptions);
            }
            else
            {
                if (ex != null)
                {
                    exceptions.Add(ex);

                    if (ex.InnerException != null)
                        GetExceptionsRecursiveImpl(ex.InnerException, exceptions);
                }
            }
        }



        /// <summary>
        /// Recurses all the inner messages of "ex", and determines the exception type of each 
        /// as an element of the returned array.
        /// This includes iteration of all exceptions within any System.AggregateException objects encountered.
        /// Return true if "ex" or any of the inner exceptions is of type "type", otherwise return false.
        /// </summary>
        public static bool AnyChildOfExceptionIsType(Exception ex, Type type)
        {
            Type[] types = GetExceptionTypesRecursive(ex);
            return types.Any(type1 => type1 == type);
        }

        /// <summary>
        /// Recurses all the inner messages of "ex", and returns the exception type of each 
        /// as an element of the returned array.
        /// The order of the returned strings is identical to the order returned by GetStackTraceRecursive(),
        /// and the number of elements in the returned array is always exactly equal to the number
        /// of exceptions recursively contained within "ex", including itself.
        /// This includes iteration of all exceptions within any System.AggregateException objects encountered.
        /// </summary>
        public static Type[] GetExceptionTypesRecursive(Exception ex)
        {
            List<Type> types = new List<Type>();
            GetExceptionTypesRecursiveImpl(ex, types);
            return types.ToArray();
        }

        private static void GetExceptionTypesRecursiveImpl(Exception ex, List<Type> types)
        {
            if (ex is AggregateException)
            {
                AggregateException multi = ex as AggregateException;

                foreach (Exception thisEx in multi.InnerExceptions)
                    GetExceptionTypesRecursiveImpl(thisEx, types);
            }
            else
            {
                types.Add(ex.GetType());

                if (ex.InnerException != null)
                    GetExceptionTypesRecursiveImpl(ex.InnerException, types);
            }
        }


        /// <summary>
        /// Returns a list of strings describing the stack trace of "ex.  All inner exceptions are recursed,
        /// and appended to this list as well.  Additionally this includes iteration of all exceptions within 
        /// any System.AggregateException objects encountered.
        /// The strings are intended to be output one per line. The
        /// order of the returned stack is such that the 0 element is the top-level function call, the Length-1
        /// element is the lowest level function that threw the original exception.  The number of elements 
        /// returned may not be the same as the actual stack depth of the inner-most exception since
        /// the StackTrace property is virtual and can be overridden in Exception-derived classes.
        /// </summary>
        public static void GetStackTraceRecursive(
            Exception ex,
            bool abbreviateText,
            out string[] functionsRet,
            out string[] linesRet,
            out string[] filesRet)
        {
            List<string> functions = new List<string>();
            List<string> lines = new List<string>();
            List<string> files = new List<string>();

            GetStackTraceRecursiveImpl(ex, abbreviateText, functions, lines, files);

            functionsRet = functions.ToArray();
            linesRet = lines.ToArray();
            filesRet = files.ToArray();
        }


        private static void GetStackTraceRecursiveImpl(
            Exception ex,
            bool abbreviateText,
            List<string> functionsRet,
            List<string> functionslinesRet,
            List<string> functionsfilesRet
            )
        {
            if (ex is AggregateException)
            {
                AggregateException multi = ex as AggregateException;

                foreach (Exception thisEx in multi.InnerExceptions)
                    GetStackTraceRecursiveImpl(thisEx, abbreviateText, functionsRet, functionslinesRet, functionsfilesRet);
            }
            else
            {
                //if (functionsRet.Count > 0)
                //{
                //    // Stick in blank lines between the stack traces of seperate exceptions for better readability
                //    functionsRet.Add("");
                //    functionslinesRet.Add("");
                //    functionsfilesRet.Add("");
                //}

                AppendStackTrace(ex, abbreviateText, functionsRet, functionslinesRet, functionsfilesRet);

                // Show inner exception stack trace last, so that the lowest-level stack frame comes last 
                // if the exception has been re-thrown, re-re-thrown, and re-re-re-thrown
                if (ex.InnerException != null)
                    GetStackTraceRecursiveImpl(ex.InnerException, abbreviateText, functionsRet, functionslinesRet, functionsfilesRet);

            }
        }


        private static void AppendStackTrace(
            Exception ex,
            bool abbreviateText,
            List<string> functionsRet,
            List<string> functionslinesRet,
            List<string> functionsfilesRet
            )
        {
            string[] stackEntries = SplitAtNewlines(ex.StackTrace);

            // Output lines from high-level to low-level:
            for (int idx = stackEntries.Length - 1; idx >= 0; idx--)
            {
                string stackEntry = stackEntries[idx];

                // Ignores if it's an empty string
                if (stackEntry.Equals(""))
                    continue;

                string functionName;
                string lineNumberString;
                string filename;

                ParseStackEntry(stackEntry, abbreviateText, out functionName, out lineNumberString, out filename);

                functionsRet.Add(functionName);
                functionslinesRet.Add(lineNumberString);
                functionsfilesRet.Add(filename);
            }
        }
        private static string CreateStackEntry(string func, string line, string file)
        {
            return string.IsNullOrEmpty(func)
                ? string.Empty
                : string.Format(" at {0} in {2}:line {1}", func, line, file);
        }

        private static void ParseStackEntry(string stackEntry, bool abbreviateText, out string func, out string line, out string file)
        {
            // Strip off the preceding "at" in a way that is immune to minor formatting
            // changes due to leading whitespace.
            const string atString = "at ";
            int  atIndex = stackEntry.IndexOf(atString, StringComparison.Ordinal);
            string cleanStackEntry = (atIndex == -1) ? stackEntry : stackEntry.Substring(atIndex + atString.Length);

            // Strip off and remember the function name
            const string inString = " in ";
            int inIndex = cleanStackEntry.IndexOf(inString, StringComparison.Ordinal);
            // For CLR functions, there is no file and line number - in thoses cases, inString is not present
            func = (inIndex == -1) ? cleanStackEntry : cleanStackEntry.Substring(0, inIndex);

            // This now contains only the file and line number
            string remainder = (inIndex == -1) ? "" : cleanStackEntry.Substring(inIndex + inString.Length);

            // Split apart into file and line number
            const string lineString = ":line ";
            int lineIndex = remainder.LastIndexOf(lineString, StringComparison.Ordinal);
            file = (lineIndex == -1) ? "" : remainder.Substring(0, lineIndex);
            line = (lineIndex == -1) ? "" : remainder.Substring(lineIndex + lineString.Length);

            if (abbreviateText)
            {
                // Rip out func args, and just leave empty parens ()
                int openingParenIdx = func.IndexOf('(');
                func = (openingParenIdx == -1) || (openingParenIdx >= func.Length - 1)
                    ? func
                    : func.Substring(0, openingParenIdx + 1) + ")";

                // Remove full directory path from file, and just return name
                file = Path.GetFileName(file);
            }
        }

        /// <summary>
        /// Returns a list of the individual lines of the StackTrace property of "ex", including all inner exceptions which are recursed,
        /// and appended to this list.  The strings include the function name, line number, and file exactly as the
        /// StackTrace property provides.
        /// </summary>
        public static string[] GetStackTraceRecursive(Exception ex, bool abbreviateText)
        {
            string[] funcs, lines, files;
            GetStackTraceRecursive(ex, abbreviateText, out funcs, out lines, out files);

            return funcs.Select((func, idx) => CreateStackEntry(funcs[idx], lines[idx], files[idx])).ToArray();
        }

        /// <summary>
        /// Returns a string array formatted with all error messages recursed from all inner exceptions,
        /// followed by a stack trace recursed similarly.  Useful for error displays where you need the srings
        /// split up indivudually such as a listbox.
        /// </summary>
        public static string[] GetComprehensiveErrorStringArray(Exception ex, bool abbreviateStackText, bool blankLineBeforeStack)
        {
            System.Collections.ArrayList strings = new System.Collections.ArrayList();

            IEnumerable<Exception> exceptions = GetExceptionsRecursive(ex);
            string[] stackEntries = GetStackTraceRecursive(ex, abbreviateStackText);

            foreach (Exception e in exceptions)
                strings.Add(string.Format("{0}: {1}", e.GetType().Name, e.Message));

            if (blankLineBeforeStack)
                strings.Add("");

            foreach (string stack in stackEntries)
                strings.Add(stack);

            return (string[])strings.ToArray(typeof(string));
        }


        /// <summary>
        /// Returns a single string formatted with all error messages recursed from all inner exceptions,
        /// and a stack trace recursed similarly.  Useful for simple error displays such as message box.
        /// </summary>
        public static string GetComprehensiveErrorString(Exception ex, bool abbreviateStackText, bool blankLineBeforeStack)
        {
            string output = "";
            string[] strings = GetComprehensiveErrorStringArray(ex, abbreviateStackText, blankLineBeforeStack);

            foreach (string str in strings)
            {
                if (output.Length > 0)
                    output += Environment.NewLine;
                output += str;
            }

            return output;
        }

        public static string GetSimpleErrorStringForDisplay(Exception ex)
        {
            string[] msgs = GetMessagesRecursive(ex);
            Array.Reverse(msgs);
            StringBuilder sb = new StringBuilder("An exception occurred. Please see the log for details. The messages are:" + Environment.NewLine + Environment.NewLine);

            foreach (string str in msgs)
            {
                sb.AppendLine(str);
            }

            return sb.ToString();
        }



        private static string[] SplitAtNewlines(string str)
        {
            return string.IsNullOrEmpty(str)
                ? new string[0]
                : str.Split(new[] { Environment.NewLine }, StringSplitOptions.None);
        }

        private static string RemoveNewlines(string str)
        {
            return string.IsNullOrEmpty(str)
                ? string.Empty
                : str.Replace(Environment.NewLine, " ");
        }
    }
}
