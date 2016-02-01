using System.Collections.Generic;
using System.Linq;

namespace ILMNcommon.Common
{
    public static class IDictionaryExtensions
    {
        /// <summary>
        /// Perform an action on each item. LINQ missed this one for some reason
        /// </summary>
        public static IDictionary<TKey, TValue> SelectKeys<TKey, TValue>(this IDictionary<TKey, TValue> dict, IEnumerable<TKey> keys)
        {
            return keys.Distinct().ToDictionary(key => key, key => dict[key]);
        }

        /// <summary>
        /// Perform an action on each item. LINQ missed this one for some reason
        /// </summary>
        public static Dictionary<TKey, TValue> ExcludeKeys<TKey, TValue>(this IDictionary<TKey, TValue> dict, IEnumerable<TKey> keys)
        {
            return dict.Keys.Where(key => !keys.Contains(key)).ToDictionary(key => key, key => dict[key]);
        }
    }
}
