using Illumina.Common.FileSystem;
using Isas.Framework.Settings;

namespace Canvas.Wrapper
{
    public class DbSnpVcfProcessor
    {
        private readonly ISampleSettings _sampleSettings;

        public DbSnpVcfProcessor(ISampleSettings sampleSettings)
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
