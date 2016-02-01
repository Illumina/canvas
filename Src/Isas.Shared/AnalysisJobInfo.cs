using System;
using System.IO;
using System.Text;
using System.Xml;

namespace Isas.Shared
{
    // universal statistics across all job-types
	// ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
	public class AnalysisJobInfo : WorkflowStatisticsCore
    {
        public string AnalysisSoftwareVersion;
        public string RTAOutputFolder;
        public string ReportPath;
		// ReSharper restore InconsistentNaming

		public const string FileName = "CompletedJobInfo.xml";
		private static readonly object ConfigLock = new object();
	
		public AnalysisJobInfo()
        {
            StartTime = DateTime.Now;
        }

        public AnalysisJobInfo(string runFolder, string softwareVersion, string analysisFolder)
            : this(runFolder, softwareVersion)
        {
            AnalysisFolder = analysisFolder;
        }


        public AnalysisJobInfo(string runFolder, string softwareVersion)
            : this()
        {
            AnalysisSoftwareVersion = softwareVersion;
            RunFolder = runFolder;
            if (Directory.Exists(runFolder))
            {
                string intensitiesDir = Path.Combine(Path.Combine(runFolder, "Data"), "Intensities");
                string rtaConfigFile = Path.Combine(intensitiesDir, "RTAConfiguration.xml");
                if (File.Exists(rtaConfigFile))
                {
                    try
                    {
                        XmlDocument document = new XmlDocument();
                        lock (ConfigLock)
                        {
                            document.Load(rtaConfigFile);
                        }
                        XmlNodeList nodeList = document.GetElementsByTagName("OutputDirectory");
                        if (nodeList != null && nodeList.Count > 0)
                        {
                            RTAOutputFolder = nodeList[0].InnerText;
                        }
                    }
                    catch (IOException)
                    {
                        RTAOutputFolder = "Error";
                    }
                }
            }
        }

        public override string ToString()
        {
            StringBuilder builder = new StringBuilder();

            builder.Append(RunFolder != null ? string.Format("[{0}]", Path.GetFileName(RunFolder)) : null);

            if (Workflow != null)
            {
                builder.Append(string.Format(" *** Workflow: {0}", Workflow));
            }
            return builder.ToString();
        }

        public string GetRunFolderName()
        {
            return new DirectoryInfo(RunFolder).Name.ToUpper();
        }

        public string GetFullRunFolderPath()
        {
            return new DirectoryInfo(RunFolder).FullName;
        }

        public override bool Equals(object obj)
        {
            if (obj is AnalysisJobInfo)
            {
                AnalysisJobInfo jobInfo = (AnalysisJobInfo) obj;
                return (jobInfo.ToString().Equals(ToString(), StringComparison.OrdinalIgnoreCase)
                    && jobInfo.AnalysisFolder.Equals(AnalysisFolder, StringComparison.OrdinalIgnoreCase));
            }

            return false;
        }

        public override int GetHashCode()
        {
            return ToString().GetHashCode();
        }
    }
}