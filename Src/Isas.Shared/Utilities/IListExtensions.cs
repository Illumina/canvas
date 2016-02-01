using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Illumina.SecondaryAnalysis
{
	public static class IListExtensions
	{
		public static T RemoveLast<T>(this IList<T> list)
		{
			T result = list.Last();
			list.RemoveAt(list.Count() - 1);
			return result;
		}
	}
}
