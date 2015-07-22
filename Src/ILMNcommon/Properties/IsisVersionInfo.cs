using System.Reflection;


// Version information for an assembly consists of the following four values:
//
//      Major Version
//      Minor Version 
//      Build Number
//      Revision
//
// You can specify all the values or you can default the Build and Revision Numbers 
// by using the '*' as shown below:
// [assembly: AssemblyVersion("1.0.*")]
[assembly: AssemblyCompany(Illumina.Shared.Version.VersionInfo.AssemblyCompany)]
[assembly: AssemblyCopyright(Illumina.Shared.Version.VersionInfo.AssemblyCopyright)]

namespace Illumina.Shared.Version
{
	internal class VersionInfo
	{
		public const string VersionString = "2.6.28";
		public const string BuildDateString = "7/13/2015 3:05:42 PM";
		public const string AssemblyCompany = CompanyName;
		public const string AssemblyCopyright = "Copyright © Illumina 2015";
		/// <summary>
		///     Important that this is the same as the RegistryKeyName
		///     because .NET will use the assembly attribute AssemblyCompany
		///     for locating all registry keys automatically.
		/// </summary>
		public const string CompanyName = RegistryKeyName;

		/// <summary>
		///     NOTE: Do not append "Inc." or anything to this - it is also used
		///     as the .NET default registry key for this app.
		/// </summary>
		public const string RegistryKeyName = "Illumina";

	}
}

