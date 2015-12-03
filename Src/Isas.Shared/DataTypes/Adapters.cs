using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Illumina.SecondaryAnalysis
{
	public class Adapters
	{
		public readonly List<string> Read1;
		public readonly List<string> Read2;

		public Adapters(List<string> read1, List<string> read2)
		{
			Read1 = read1;
			Read2 = read2;
		}

		public IEnumerable<string> All()
		{
			return Read1.Union(Read2);
		}
	}
}
