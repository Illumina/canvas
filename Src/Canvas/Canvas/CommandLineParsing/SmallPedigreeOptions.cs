using System.Collections.Generic;
using Illumina.Common.FileSystem;

namespace Canvas.CommandLineParsing
{
    public class SmallPedigreeOptions
    {
        public List<SmallPedigreeSampleOptions> Samples { get; }
        public IFileLocation CommonCnvsBed { get; }
        public IFileLocation BAlleleSites { get; }
        public bool IsPopulationBAlleleSites { get; }
        public IFileLocation MultiSamplePloidyVcf { get; }


        public SmallPedigreeOptions(List<SmallPedigreeSampleOptions> samples, IFileLocation commonCnvsBed, IFileLocation bAlleleSites, bool isPopulationBAlleleSites, IFileLocation multiSamplePloidyVcf)
        {
            Samples = samples;
            CommonCnvsBed = commonCnvsBed;
            BAlleleSites = bAlleleSites;
            MultiSamplePloidyVcf = multiSamplePloidyVcf;
            IsPopulationBAlleleSites = isPopulationBAlleleSites;
        }
    }
}