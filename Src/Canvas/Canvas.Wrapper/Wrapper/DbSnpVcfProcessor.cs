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
            System.Console.WriteLine("%%% GetDbSnpVcfPath!");
            foreach (var item in _sampleSettings.GetSettingKeys())
            {
                System.Console.WriteLine("{0} - {1}", item, _sampleSettings.GetStringSetting(item, "default"));
            }
            string dbSnpVcfPath = _sampleSettings.GetStringSetting("dbsnpvcfpath", null);
            IFileLocation dbSnpVcf = null;
            if (dbSnpVcfPath != null)
                dbSnpVcf = new FileLocation(dbSnpVcfPath);
            System.Console.WriteLine("Return dbSNPVCF - '{0}'", dbSnpVcf);
            return dbSnpVcf;
        }
    }
}
