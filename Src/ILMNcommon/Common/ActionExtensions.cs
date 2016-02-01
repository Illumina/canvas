using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ILMNcommon.Common
{
    public static class ActionExtensions
    {
        public static bool Try(this Action a, out Exception e)
        {
            e = null;
            try
            {
                a.Invoke();
                return true;
            }
            catch (Exception ex)
            {
                e = ex;
                return false;
            }
        }

        public static bool Try(this Action a)
        {
            Exception e;
            return a.Try(out e);
        }

        /// <summary>
        /// Run the action and wrap any exception with a newly created exception 
        /// </summary>
        /// <param name="a">The action to run</param>
        /// <param name="createFromNestedException">The function for creating the new exception from the nested exception</param>
        public static void WrapAnyExceptionWith(this Action a, Func<Exception, Exception> createFromNestedException)
        {
            Exception e;
            if (!a.Try(out e))
                throw createFromNestedException(e);
        }
    }
}
