using Illumina.Common.FileSystem;
using Isas.Framework.Settings;

namespace Canvas.Wrapper
{
    public class DbSnpVcfProcessor
    {
        private readonly ISettings _sampleSettings;

        [System.Obsolete("Pass in ISettings")]
        public DbSnpVcfProcessor(ISampleSettings sampleSettings) : this((ISettings)sampleSettings)
        {
        }
        public DbSnpVcfProcessor(ISettings sampleSettings)
        {
            _sampleSettings = sampleSettings;
        }


        public IFileLocation GetDbSnpVcfPath()
        {
            return _sampleSettings.GetSetting(DbSnpVcfPathSetting);
        }

        public Setting<IFileLocation> DbSnpVcfPathSetting = SampleSettings.CreateSetting<IFileLocation>(
            "DBSnpVcfPath",
            "Path to a vcf file containing SNP sites to use for calculating B-allele frequencies during CNV calling",
            null,
            null,
            path => new FileLocation(path));
    }
}
