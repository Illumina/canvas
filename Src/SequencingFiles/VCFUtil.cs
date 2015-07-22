using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Text;
using Extensions;
namespace Extensions
{
	public static class Extension
	{
		public static void WriteTab(this SequencingFiles.BgzipOrStreamWriter w, string s)
		{
			w.Write("{0}\t", s);
		}
	}
}

// This module contains a simplified general-purpose VCF writer

namespace SequencingFiles.VCFUtil
{
	/// <summary>
	///     generalized variant represented by each line of a VCF record
	/// </summary>
	public class VCFLocus : IComparable
	{
		public List<string> Alt;
		public string ChromName;
		public List<string> Filter;
		public string Id;
		public uint Pos;
		public string Ref;

		public VCFLocus()
		{
			Alt = new List<string>();
			Filter = new List<string>();
			_info = new List<KeyValuePair<string, string>>();
			_sampleValues = new List<KeyValuePair<string, string[]>>();
		}

		int IComparable.CompareTo(object obj)
		{
			VCFLocus rhs = (VCFLocus)obj;
			if (rhs == null) return -1;
			int res = String.Compare(ChromName, rhs.ChromName);
			if (res != 0) return res;
			if (Pos > rhs.Pos) return 1;
			if (Pos < rhs.Pos) return -1;
			return 0;
		}

		#region private

		private readonly List<KeyValuePair<string, string>> _info;
		private readonly List<KeyValuePair<string, string[]>> _sampleValues;
		private string _qual;

		private static string StringCheck(string s)
		{
			if (String.IsNullOrEmpty(s)) return ".";
			return s;
		}

		private void WriteStringFieldWithTab(BgzipOrStreamWriter w, string s)
		{
			w.WriteTab(StringCheck(s));
		}

		private void WriteStringFieldNoTab(BgzipOrStreamWriter w, string s)
		{
			w.Write(StringCheck(s));
		}

		private string GetAltString()
		{
			return String.Join(",", Alt.ToArray());
		}

		private string GetFilterString()
		{
			if (Filter.Count == 0) return "PASS";
			return String.Join(",", Filter.ToArray());
		}

		private string GetInfoString()
		{
			StringBuilder sb = new StringBuilder();
			foreach (KeyValuePair<string, string> kv in _info)
			{
				if (kv.Key == null) continue;
				if (sb.Length > 0) sb.Append(";");
				sb.Append(kv.Key);
				if (kv.Value == null) continue;
				sb.Append('=' + kv.Value);
			}
			return sb.ToString();
		}

		private string GetFormatString()
		{
			StringBuilder sb = new StringBuilder();
			foreach (KeyValuePair<string, string[]> kv in _sampleValues)
			{
				if (sb.Length > 0) sb.Append(':');
				sb.Append(kv.Key);
			}
			return sb.ToString();
		}

		private string GetSampleString(int sampleIndex)
		{
			StringBuilder sb = new StringBuilder();
			foreach (KeyValuePair<string, string[]> kv in _sampleValues)
			{
				if (sb.Length > 0) sb.Append(':');
				if (sampleIndex < kv.Value.Length)
				{
					sb.Append(kv.Value[sampleIndex]);
				}
				else
				{
					sb.Append('.');
				}
			}
			return sb.ToString();
		}

		#endregion

		public void Write(BgzipOrStreamWriter w, int nSamples)
		{
			SanityCheck(nSamples);

			WriteStringFieldWithTab(w, ChromName);
			w.WriteTab(Pos.ToString());
			WriteStringFieldWithTab(w, Id);
			WriteStringFieldWithTab(w, Ref);
			WriteStringFieldWithTab(w, GetAltString());
			WriteStringFieldWithTab(w, _qual);
			WriteStringFieldWithTab(w, GetFilterString());
			WriteStringFieldWithTab(w, GetInfoString());
			WriteStringFieldWithTab(w, GetFormatString());
			for (int sampleIndex = 0; sampleIndex < nSamples; ++sampleIndex)
			{
				if (sampleIndex == (nSamples - 1))
					WriteStringFieldNoTab(w, GetSampleString(sampleIndex));
				else
					WriteStringFieldWithTab(w, GetSampleString(sampleIndex));
			}
			w.WriteLine();
		}

		// Qual is protected to preserve/enforce floating point representation:
		public string Qual()
		{
			return _qual;
		}

		public void SetQual(string value)
		{
			double test;
			if (!Double.TryParse(value, out test))
			{
				throw new Exception(String.Format("Can't set VCF qual field to '{0}'", value));
			}
			_qual = value;
		}

		public void AddInfo(string key, string val)
		{
			_info.Add(new KeyValuePair<string, string>(key, val));
		}

		public void AddInfo(string key)
		{
			AddInfo(key, null);
		}

		public void AddSampleValue(string key, string[] val)
		{
			_sampleValues.Add(new KeyValuePair<string, string[]>(key, val));
		}

		public void AddSampleValue(string key, string val)
		{
			AddSampleValue(key, new[] { val });
		}

		public void SanityCheck(int nSamples)
		{
			foreach (KeyValuePair<string, string[]> kv in _sampleValues)
			{
				Debug.Assert(kv.Key != null);
				Debug.Assert(kv.Value.Length <= nSamples);
			}
		}
	}


	public class VCFOutStreamer
	{
		// TODO: handle vcf version specializations by inheritance -- for now all cases are vcf 4.1
		public const string VCFVersion = "4.1";

		private readonly string[] _sampleIDs;
		private readonly List<KeyValuePair<string, string>> _headerList;
		private readonly string[] _sampleNames;
		private readonly BgzipOrStreamWriter _writer;
		private bool _isHeaderWritten;

		public VCFOutStreamer(BgzipOrStreamWriter writer, string[] sampleNames, string[] sampleIds)
		{
			_sampleNames = sampleNames;
			_sampleIDs = sampleIds;
			_writer = writer;
			_headerList = new List<KeyValuePair<string, string>>();
		}

		public VCFOutStreamer(BgzipOrStreamWriter writer, string sampleName, string sampleId)
			: this(writer, new[] { sampleName }, new[] { sampleId })
		{
		}

		public VCFOutStreamer(BgzipOrStreamWriter writer)
			: this(writer, "SAMPLE", "SAMPLE")
		{
		}

		public void AddHeaderKeyVal(string key, string val)
		{
			// TODO: check for repeated keys:
			_headerList.Add(new KeyValuePair<string, string>(key, val));
		}

		public void WriteHeader()
		{
			Debug.Assert(!_isHeaderWritten);
			if (_isHeaderWritten) return;

			// write some custom fields 'by hand':
			_writer.WriteLine("##fileformat=VCFv" + VCFVersion);
			_writer.WriteLine("##fileDate=" + string.Format("{0:yyyyMMdd}", DateTime.Now));

			foreach (KeyValuePair<string, string> item in _headerList)
			{
				_writer.WriteLine("##{0}={1}", item.Key, item.Value);
			}

			_writer.Write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
			foreach (string id in _sampleIDs)
			{
				_writer.Write("\t{0}", id);
			}
			_writer.WriteLine();
			_isHeaderWritten = true;
		}

		public void WriteLocus(VCFLocus locus)
		{
			if (locus == null) return;
			locus.Write(_writer, _sampleNames.Length);
		}
	}

#if FOO
    // roughly how I'd like to interact with this class:

    x = new VCFOutStreamer(FileStream,SampleNames);

    x = AddHeaderKeyVal();
    x = AddHeadLine();
    x.writeHeader;
    rec = new VCFRecord(chrom,pos);
    x.Write(rec);
#endif
}