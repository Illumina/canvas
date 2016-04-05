using SampleSettingsProcessing;

namespace Isas.Shared
{
    public class WorkManagerFactory
    {
        public static LocalWorkManager GetLocalWorkManager(ILogger logger, IDirectoryLocation analysisFolder,
            ISampleSettings settings, IsasConfiguration config)
        {
            int maxThreadCount = settings.GetIntegerSetting(SampleSheetUtils.MaximumThreadCount, 0);
            float maxMemGB = settings.GetFloatSetting(SampleSheetUtils.MaximumMemoryGB, config.MaximumMemoryGB);
            var loggingFolder = analysisFolder.CreateSubdirectory("Logging");

            return new LocalWorkManager(logger, loggingFolder, maxThreadCount, maxMemGB, config.MaximumHoursPerProcess);
        }

        public static SGEWorkManager GetSGEWorkManager(ILogger logger, IDirectoryLocation analysisFolder,
            ISampleSettings settings, IsasConfiguration config)
        {
            LocalWorkManager local = GetLocalWorkManager(logger, analysisFolder, settings, config);
            return new SGEWorkManager(logger, local,
                config.ClusterQueueName, config.ClusterParallelEnvironmentName);
        }

        public static IWorkManager GetWorkManager(ILogger logger, IDirectoryLocation analysisFolder,
            ISampleSettings settings, IsasConfiguration config)
        {
            bool useSGE = settings.GetBooleanSetting(SampleSheetUtils.UseSGE, config.UseCluster);
            config.UseCluster = useSGE;

            if (config.UseCluster)
            {
                return GetSGEWorkManager(logger, analysisFolder, settings, config);
            }
            else
            {
                return GetLocalWorkManager(logger, analysisFolder, settings, config);
            }
        }
    }
}
