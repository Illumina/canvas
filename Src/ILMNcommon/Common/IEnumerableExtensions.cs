using System;
using System.Collections.Generic;
using System.Linq;

namespace ILMNcommon.Common
{
    public static class IEnumerableExtensions
    {
        /// <summary>
        /// Perform an action on each item. LINQ missed this one for some reason
        /// </summary>
        public static void ForEach<T>(this IEnumerable<T> value, Action<T> action)
        {
            foreach (T item in value)
            {
                action(item);
            }
        }

        public static bool IsNullOrEmpty<T>(this IEnumerable<T> items)
        {
            return items == null || items.Empty();
        }

        public static bool Empty<T>(this IEnumerable<T> items)
        {
            return !items.Any();
        }

        /// <summary>
        /// return true if predicate is either
        ///     true for all items
        ///            OR
        ///     false for all items
        /// </summary>
        public static bool AllOrNone<T>(this IEnumerable<T> items, Func<T, bool> predicate)
        {
            return items.All(predicate) || items.All(item => !predicate(item));
        }

        public static bool AllAreNullOrNoneAreNull(params object[] items)
        {
            return items.AllOrNone(item => item != null);
        }

        public static bool AnyNull(params object[] items)
        {
            return items.Any(item => item == null);
        }

        public static T MaxBy<T>(this IEnumerable<T> items, Func<T, int> transform)
        {
            T itemWithMaxValue = items.FirstOrDefault();
            int maxValue = itemWithMaxValue == null ? int.MinValue : transform(itemWithMaxValue);
            foreach (var item in items)
            {
                int value = transform(item);
                if (value <= maxValue) continue;
                itemWithMaxValue = item;
                maxValue = value;

            }
            return itemWithMaxValue;
        }

        public static IEnumerable<Tuple<T, T>> GetPairs<T>(this IEnumerable<T> items)
        {
            // IEnumerable doesn't guarantee the same order for multiple enumerations, but IList does
            IList<T> itemsList = items.ToList();
            return itemsList.Select((value, index) => new { value, index })
                            .SelectMany(x => itemsList.Skip(x.index + 1),
                                        (x, y) => Tuple.Create(x.value, y));
        }

        public static IEnumerable<T> ToSingleItemEnumerable<T>(this T item)
        {
            yield return item;
        }
    }
}
