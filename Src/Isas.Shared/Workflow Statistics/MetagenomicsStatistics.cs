using System.Collections.Generic;
using System.IO;
using System.Xml.Serialization;

namespace Isas.Shared
{
    // ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
    public class MetagenomicsStatistics
    {
        public long[] ClassifiedReadCount = new long[MetagenomicsConstants.NumTaxonomicLevels]; // by level
        public long ClusterCount;
        public long ClusterCountPF;
        public string SampleID;
        public string SampleName;
        public int SampleNumber;
        // ReSharper restore InconsistentNaming

        public void Serialize(string saveFileName)
        {
            XmlSerializer serializer = new XmlSerializer(typeof(MetagenomicsStatistics));

            using (FileStream stream = new FileStream(saveFileName, FileMode.Create))
            {
                serializer.Serialize(stream, this);
            }
        }
    }

    // ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
    public class StatisticsMetagenomics : WorkflowStatistics
    {
        public MetagenomicsStatistics[] StatisticsBySample;
        // ReSharper restore InconsistentNaming

        public StatisticsMetagenomics()
        {
            //these are not used so set them to null so that they aren't serialized
            OverallSamples = null;
            PairedEndByGenome = null; // by genome
            StatsSamples = null;
        }
    }
}