using Isas.Shared;
using System.Collections.Generic;

namespace Isas.Shared
{
	public interface IMoveable<out T>
	{
		T Move(FileNamingConvention newName);
	}

	/// <summary>
	/// Naming convention for files identified by a SampleInfo
	/// </summary>
	/// <param name="sample">The associated SampleInfo</param>
	/// <param name="filename">The name of the original file: i.e. IFileLocation.Name</param>
	/// <returns>An IFileLocation named by this naming convention.</returns>
	public delegate IFileLocation SampleNamingConvention(SampleInfo sample, string filename);

	public delegate IFileLocation SampleStubNamingConvention(SampleInfo sample);

	/// <summary>
	/// Naming convention for files identified by order in an enumeration (e.g. IEnumerable<Fastq>)
	/// </summary>
	/// <param name="index">The index of the element in the enumeration</param>
	/// <param name="filename">The name of the original file: i.e. IFileLocation.Name</param>
	/// <returns>An IFileLocation named by this naming convention.</returns>
	public delegate IFileLocation EnumerableNamingConvention(int index, string filename);

	/// <summary>
	/// Naming convention for files identified by both a SampleInfo and the order in an enumeration
	/// (e.g. SampleSet<IEnumberable<Fastq>>)
	/// </summary>
	/// <param name="sample">The associated SampleInfo</param>
	/// <param name="index">The index of the element in the enumeration</param>
	/// <param name="filename">The name of the original file: i.e. IFileLocation.Name</param>
	/// <returns>An IFileLocation named by this naming convention.</returns>
	public delegate IFileLocation SampleEnumerableNamingConvention(SampleInfo sample, int index, string filename);

	/// <summary>
	/// Basic naming convention
	/// </summary>
	/// <param name="filename">The name of the original file: i.e. IFileLocation.Name</param>
	/// <returns>An IFileLocation named by this naming convention.</returns>
	public delegate IFileLocation FileNamingConvention(string filename);

	public class NamingConventions
	{
		/// <summary>
		/// Strips <paramref name="sample"/>.Id and <paramref name="sample"/>.Name,
		/// if present, from the filename before returning the new naming stub.
		/// This helps to avoid repeatedly preppending identifiers to the filenames
		/// during subsequent moves.
		/// </summary>
		/// <param name="sample">The SampleInfo associated with the file.</param>
		/// <param name="filename">The name of the file (i.e. IFileLocation.Name)</param>
		/// <returns>A new version of <paramref name="filename"/> without duplicate <paramref name="sample"/> information.</returns>
		public static string GetStandardSampleStub(SampleInfo sample, string filename)
		{
			// Make sure filename doesn't have ID and Name in it already
			// We don't want to keep preppending ID_Name everytime we move something.
			if (filename.StartsWith(sample.Id + "_"))
			{
				filename = filename.Substring(sample.Id.Length + 1);
			}
			if (filename.StartsWith(sample.Name + "_"))
			{
				filename = filename.Substring(sample.Name.Length + 1);
			}
			return string.Format("{0}_{1}_{2}", sample.Id, sample.Name, filename);
		}

		/// <summary>
		/// Strips <paramref name="sample"/>.Id , <paramref name="sample"/>.Name,
		/// and <paramref name="index"/> information,
		/// if present, from the filename before returning the new naming stub.
		/// This helps to avoid repeatedly preppending identifiers to the filenames
		/// during subsequent moves.
		/// </summary>
		/// <param name="sample">The SampleInfo associated with the file.</param>
		/// <param name="index">The index associated with the file.</param>
		/// <param name="filename">The name of the file (i.e. IFileLocation.Name)</param>
		/// <returns>A new version of <paramref name="filename"/> without duplicate <paramref name="sample"/> 
		/// and <paramref name="index"/> information.</returns>
		public static string GetStandardSampleEnumerableStub(SampleInfo sample, int index, string filename)
		{
			// Make sure filename doesn't have ID and Name in it already
			// We don't want to keep preppending ID_Name everytime we move something.
			if (filename.StartsWith(sample.Id + "_"))
			{
				filename = filename.Substring(sample.Id.Length + 1);
			}
			if (filename.StartsWith(sample.Name + "_"))
			{
				filename = filename.Substring(sample.Name.Length + 1);
			}
			if (filename.StartsWith(index.ToString() + "_"))
			{
				filename = filename.Substring(index.ToString().Length + 1);
			}
			return string.Format("{0}_{1}_{2}_{3}", sample.Id, sample.Name, index, filename);
		}

		/// <summary>
		/// OriginalDir/OrignalName.Ext => <paramref name="targetDir"/>/SampleID_SampleName_OriginalName.Ext
		/// </summary>
		/// <param name="targetDir">The target directory to be used in the naming convention.</param>
		/// <returns>A SampleNamingConvention function.</returns>
		public static SampleNamingConvention GetStandardSampleNamingConvention(IDirectoryLocation targetDir)
		{
			return delegate(SampleInfo sample, string filename)
			{
				return targetDir.GetFileLocation(GetStandardSampleStub(sample, filename));
			};
		}

		/// <summary>
		/// OriginalDir/OrignalName.Ext => <paramref name="targetDir"/>/SampleID/SampleID_SampleName_OriginalName.Ext
		/// </summary>
		/// <param name="targetDir">The target directory to be used in the naming convention.</param>
		/// <returns>A SampleNamingConvention function.</returns>
		public static SampleNamingConvention GetNestedSampleNamingConvention(IDirectoryLocation targetDir)
		{
			return delegate(SampleInfo sample, string filename)
			{
				return targetDir.CreateSubdirectory(sample.Id)
								.GetFileLocation(GetStandardSampleStub(sample, filename));
			};
		}

		/// <summary>
		/// OriginalDir/OrignalName.Ext => <paramref name="targetDir"/>/SampleID_SampleName_#_OriginalName.Ext
		/// The # comes from the order the objects are enumerated.
		/// </summary>
		/// <param name="targetDir">The target directory to be used in the naming convention.</param>
		/// <returns>A SampleEnumerableNamingConvention function.</returns>
		public static SampleEnumerableNamingConvention GetStandardSampleEnumerableNamingConvention(IDirectoryLocation targetDir)
		{
			return delegate(SampleInfo sample, int index, string filename)
			{
				return targetDir.GetFileLocation(GetStandardSampleEnumerableStub(sample, index, filename));
			};
		}

		/// <summary>
		/// OriginalDir/OrignalName.Ext => <paramref name="targetDir"/>/SampleID/SampleID_SampleName_#_OriginalName.Ext
		/// The # comes from the order the objects are enumerated.
		/// </summary>
		/// <param name="targetDir">The target directory to be used in the naming convention.</param>
		/// <returns>A SampleEnumerableNamingConvention function.</returns>
		public static SampleEnumerableNamingConvention GetNestedSampleEnumerableNamingConvention(IDirectoryLocation targetDir)
		{
			return delegate(SampleInfo sample, int index, string filename)
			{
				return targetDir.CreateSubdirectory(sample.Id)
								.GetFileLocation(GetStandardSampleEnumerableStub(sample, index, filename));
			};
		}
	}

	public static class SampleSetMove
	{
		/// <summary>
		/// Move all members of the SampleSet according to the naming convention
		/// </summary>
		/// <typeparam name="Tmove">Must implement IMoveable<<typeparamref name="Tmove"/>></typeparam>
		/// <param name="samples">The instance of the SampleSet this method extends</param>
		/// <param name="sampleNaming">The sample naming convention used to rename the samples</param>
		/// <returns>A SampleSet of <typeparamref name="Tmove"/></returns>
		public static SampleSet<Tmove> Move<Tmove>(this SampleSet<Tmove> samples,
			SampleNamingConvention sampleNaming) where Tmove : IMoveable<Tmove>
		{
			SampleSet<Tmove> movedSamples = new SampleSet<Tmove>();
			foreach (var sample in samples)
			{
				FileNamingConvention naming = delegate(string filename)
				{
					return sampleNaming(sample.Key, filename);
				};
				movedSamples.Add(sample.Key, sample.Value.Move(naming));
			}
			return movedSamples;
		}

		/// <summary>
		/// Move all IEnumerable collections of IMoveable objects in the SampleSet according to the naming convention
		/// </summary>
		/// <typeparam name="Tmove">Must implement IMoveable<<typeparamref name="Tmove"/>></typeparam>
		/// <param name="samples">The instance of the SampleSet this method extends</param>
		/// <param name="sampleNaming">The sample naming convention used to rename the samples</param>
		/// <returns>A SampleSet of IEnumerable of <typeparamref name="Tmove"/></returns>
		public static SampleSet<IEnumerable<Tmove>> Move<Tmove>(this SampleSet<IEnumerable<Tmove>> samples,
			SampleEnumerableNamingConvention sampleNaming) where Tmove : IMoveable<Tmove>
		{
			SampleSet<IEnumerable<Tmove>> movedSamples = new SampleSet<IEnumerable<Tmove>>();
			foreach (var sampleList in samples)
			{
				// Create an EnumerableNamingConvention by plugging the SampleInfo (sampleList.Key)
				//  into the SampleEnumerableNamingConvention.
				EnumerableNamingConvention listNaming = delegate(int index, string filename)
				{
					return sampleNaming(sampleList.Key, index, filename);
				};
				IEnumerable<Tmove> newSamples = sampleList.Value.Move(listNaming);
				movedSamples.Add(sampleList.Key, newSamples);
			}
			return movedSamples;
		}

		/// <summary>
		/// Move an IEnumerable collection of IMoveable objects according to the naming convention
		/// </summary>
		/// <typeparam name="Tmove">Must implement IMoveable<<typeparamref name="Tmove"/>></typeparam>
		/// <param name="samples">The instance of the SampleSet this method extends</param>
		/// <param name="sampleNaming">The sample naming convention used to rename the samples</param>
		/// <returns>An IEnumerable of <typeparamref name="Tmove"/></returns>
		public static IEnumerable<Tmove> Move<Tmove>(this IEnumerable<Tmove> samples,
			EnumerableNamingConvention namingConvention) where Tmove : IMoveable<Tmove>
		{
			List<Tmove> newSamples = new List<Tmove>();
			foreach (IMoveable<Tmove> sample in samples)
			{
				FileNamingConvention getFileName = delegate(string filename)
				{
					return namingConvention(newSamples.Count, filename);
				};
				newSamples.Add(sample.Move(getFileName));
			}
			return newSamples;
		}
	}

	public static class FileLocationMove
	{
		/// <summary>
		/// Move this IFileLocation according to the <paramref name="namingConvention"/> and add a symlink
		/// in the original location.
		/// </summary>
		/// <param name="source">This IFileLocation</param>
		/// <param name="namingConvention">The NamingConvention that retuns the destination IFileLocation</param>
		/// <returns>The destination IFileLocation</returns>
		public static IFileLocation Move(this IFileLocation source, FileNamingConvention namingConvention)
		{
			return source.MoveAndLink(namingConvention(source.Name));
		}

		/// <summary>
		/// Move this IDirectoryLocation according to the <paramref name="namingConvention"/> and add a symlink
		/// in the original location.
		/// </summary>
		/// <param name="source">This IDirectoryLocation</param>
		/// <param name="namingConvention">The NamingConvention that retuns the destination IFileLocation. 
		/// While the return type of <paramref name="namingConvention"/> is an IFileLocation, this is converted
		/// to an IDirectoryLocation.</param>
		/// <returns>The destination IDirectoryLocation</returns>
		public static IDirectoryLocation Move(this IDirectoryLocation source, FileNamingConvention namingConvention)
		{
			// FileNamingConvention returns a file, but we'll use a hack to make it a directory
			IFileLocation fakeDestination = namingConvention(source.Name);
			IDirectoryLocation realDestination = fakeDestination.Directory.CreateSubdirectory(fakeDestination.Name);
			return source.MoveAndLink(realDestination);
		}

		/// <summary>
		/// Move all IFileLocations in the SampleSet according to the naming convention
		/// </summary>
		/// <param name="samples">The instance of the SampleSet this method extends</param>
		/// <param name="sampleNaming">The sample naming convention used to rename the samples</param>
		/// todo: make IFileLocation implement IMoveable and then remove this method
		/// <returns>A SampleSet of IFileLocation</returns>
		public static SampleSet<IFileLocation> Move(this SampleSet<IFileLocation> samples,
			SampleNamingConvention sampleNaming)
		{
			SampleSet<IFileLocation> movedSamples = new SampleSet<IFileLocation>();
			foreach (var sample in samples)
			{
				FileNamingConvention naming = delegate(string filename)
				{
					return sampleNaming(sample.Key, filename);
				};
				movedSamples.Add(sample.Key, sample.Value.Move(naming));
			}
			return movedSamples;
		}

		/// <summary>
		/// Move all IDirectoryLocation in the SampleSet according to the naming convention
		/// </summary>
		/// <param name="samples">The instance of the SampleSet this method extends</param>
		/// <param name="sampleNaming">The sample naming convention used to rename the samples</param>
		/// todo: make IDirectoryLocation implement IMoveable and then remove this method
		/// <returns>A SampleSet of IDirectoryLocation</returns>
		public static SampleSet<IDirectoryLocation> Move(this SampleSet<IDirectoryLocation> samples,
			SampleNamingConvention sampleNaming)
		{
			SampleSet<IDirectoryLocation> movedSamples = new SampleSet<IDirectoryLocation>();
			foreach (var sample in samples)
			{
				FileNamingConvention naming = delegate(string filename)
				{
					return sampleNaming(sample.Key, filename);
				};
				movedSamples.Add(sample.Key, sample.Value.Move(naming));
			}
			return movedSamples;
		}
	}
}
