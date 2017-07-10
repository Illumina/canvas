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
            string dbSnpVcfPath = _sampleSettings.GetStringSetting("dbsnpvcfpath", null);
            IFileLocation dbSnpVcf = null;
            if (dbSnpVcfPath != null)
                dbSnpVcf = new FileLocation(dbSnpVcfPath);
            return dbSnpVcf;
        }
    }
}
