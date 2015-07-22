using System;
using System.IO;

namespace SequencingFiles
{
    public class IntensityStats
    {
        public float Dev;
        public float Max;
        public float Mean;
        public float Med;
        public float Min;
        public float P90;

        public void Save(string filePath)
        {
            using (TextWriter writer = new StreamWriter(filePath))
            {
                writer.WriteLine(Mean);
                writer.WriteLine(Med);
                writer.WriteLine(Min);
                writer.WriteLine(Max);
                writer.WriteLine(P90);
                writer.WriteLine(Dev);
            }
        }

        public static IntensityStats Load(string fullFileName)
        {
            IntensityStats stats = new IntensityStats();
            using (TextReader reader = new StreamReader(fullFileName))
            {
                stats = new IntensityStats
                            {
                                Mean = Convert.ToSingle(reader.ReadLine()),
                                Med = Convert.ToSingle(reader.ReadLine()),
                                Min = Convert.ToSingle(reader.ReadLine()),
                                Max = Convert.ToSingle(reader.ReadLine()),
                                P90 = Convert.ToSingle(reader.ReadLine()),
                                Dev = Convert.ToSingle(reader.ReadLine())
                            };
            }
            return stats;
        }
    }
}