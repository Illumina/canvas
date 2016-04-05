using System;
using System.IO;
using System.Text.RegularExpressions;

namespace Isas.Shared
{
    /// <summary>
    /// Encapsulates aligner stats .xml file output by Isaac, where we can parse metrics (e.g. # of adapter bases) not obtainable from 
    /// the .bam file.  Each instance is particular to a file location (i.e. when the file moves, you get a new IsaacStatsXML file)
    /// </summary>
    public class IsaacStatsXML : IMoveable<IsaacStatsXML>, IMoveableResult<IsaacStatsXML>
    {
        public readonly IFileLocation XmlPath;

        // used by json.net to bypass existence check
        private IsaacStatsXML() { }

        public IsaacStatsXML(IFileLocation path)
        {
            this.XmlPath = path;
        }

        public IsaacStatsXML Move(FileNamingConvention newName)
        {
            IFileLocation newXML = XmlPath.Move(newName);
            return new IsaacStatsXML(newXML);
        }

        public IsaacStatsXML Move(IFileLocation newPath)
        {
            IsaacStatsXML movedXML = new IsaacStatsXML(newPath);
            return movedXML;
        }

        public IsaacStatsXML Move(IsaacStatsXML destination)
        {
            XmlPath.MoveAndLink(destination.XmlPath);
            return destination;
        }

        public IsaacStatsXML MoveWithStub(IFileLocation newStub, Action<IFileLocation, IFileLocation> move)
        {
            IFileLocation newPath = newStub.AppendName(".bam.xml");
            return Move(newPath, move);
        }


        public IsaacStatsXML Move(IFileLocation newPath, Action<IFileLocation, IFileLocation> move)
        {
            move(XmlPath, newPath);
            IsaacStatsXML movedXML = new IsaacStatsXML(newPath);
            return movedXML;
        }
        public void Delete()
        {
            XmlPath.Delete();
        }

        public long? ParseAdapterBases()
        {
            long? result = null;

            if (!File.Exists(XmlPath.FullName))
            {
                // Return null if the .bam.xml file does not exist (e.g. if we resumed a workflow)
                Console.Error.Write("Warning: xml file not found at '{0}'", XmlPath.FullName);
                return result;
            }
            Regex adapterBasesRegex = new Regex(@"<AdapterBases>(\d+)</AdapterBases>");
            try
            {
                result = 0;
                using (StreamReader reader = new StreamReader(XmlPath.FullName))
                {
                    while (true)
                    {
                        string fileLine = reader.ReadLine();
                        if (fileLine == null) break;
                        Match match = adapterBasesRegex.Match(fileLine);
                        if (match != null && match.Success)
                        {
                            result += long.Parse(match.Groups[1].Value);
                        }
                    }
                }
            }
            catch (Exception ex)
            {
                Console.Error.WriteLine("Exception parsing stats file at '{0}'", XmlPath.FullName);
                Console.Error.WriteLine("Details: {0}", ex.ToString());
                result = null;
            }
            return result;
        }

    }
}