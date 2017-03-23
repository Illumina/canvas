using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Illumina.Common.FileSystem;
using Isas.Framework.DataTypes;
using Isas.Framework.WorkManagement;
using Isas.SequencingFiles;

namespace Canvas.Wrapper
{
    public interface ICanvasSingleSampleInputCommandLineBuilder
    {
        StringBuilder GetSingleSampleCommandLine(string sampleId, Bam bam, GenomeMetadata genomeMetadata, IDirectoryLocation sampleSandbox);
        StringBuilder GetCustomParameters(Dictionary<string, string> moreCustomParameters = null);
        StringBuilder MergeCustomCanvasParameters(StringBuilder commandLine);
    }

    /// <summary>
    /// Run Canvas on enrichment data to generate CNV calls:
    /// </summary>
    public class CanvasSingleSampleInputCommandLineBuilder : ICanvasSingleSampleInputCommandLineBuilder
    {
        private readonly ICanvasAnnotationFileProvider _annotationFileProvider;
        private readonly Dictionary<string, string> _customParameters;
        private readonly string _customCanvasParameters;

        public CanvasSingleSampleInputCommandLineBuilder(ICanvasAnnotationFileProvider annotationFileProvider, Dictionary<string, string> customParameters, string customCanvasParameters)
        {
            _annotationFileProvider = annotationFileProvider;
            _customParameters = customParameters;
            _customCanvasParameters = customCanvasParameters;
        }

        public StringBuilder GetSingleSampleCommandLine(string sampleId, Bam bam, GenomeMetadata genomeMetadata, IDirectoryLocation sampleSandbox)
        {
            StringBuilder commandLine = new StringBuilder();

            commandLine.Append($" --bam \"{bam.BamFile}\"");
            commandLine.Append($" --sample-name \"{sampleId}\"");
            IFileLocation kmerFasta = _annotationFileProvider.GetKmerFasta(genomeMetadata);
            commandLine.Append($" --reference \"{kmerFasta}\"");
            IDirectoryLocation wholeGenomeFasta = new FileLocation(genomeMetadata.Sequences.First().FastaPath).Directory;
            commandLine.Append($" --genome-folder \"{wholeGenomeFasta}\"");
            IFileLocation filterBed = _annotationFileProvider.GetFilterBed(genomeMetadata);
            commandLine.Append($" --filter-bed \"{filterBed}\"");
            commandLine.Append($" --output \"{sampleSandbox}\"");

            return commandLine;
        }

        public StringBuilder GetCustomParameters(Dictionary<string, string> moreCustomParameters = null)
        {
            moreCustomParameters = moreCustomParameters ?? new Dictionary<string, string>();
            var allModules = new HashSet<string>(_customParameters.Keys, StringComparer.OrdinalIgnoreCase);
            allModules.UnionWith(moreCustomParameters.Keys);
            var commandLine = new StringBuilder();
            foreach (var module in allModules)
            {
                var options = "";
                if (_customParameters.ContainsKey(module))
                    options = _customParameters[module];
                if (moreCustomParameters.ContainsKey(module))
                    options = Isas.Framework.Settings.CommandOptionsUtilities.MergeCommandLineOptions(options, moreCustomParameters[module]);
                string customParameters = module + "," + options;
                commandLine.Append($" --custom-parameters \"{customParameters}\"");
            }
            return commandLine;
        }

        public StringBuilder MergeCustomCanvasParameters(StringBuilder commandLine)
        {
            return new StringBuilder(Isas.Framework.Settings.CommandOptionsUtilities.MergeCommandLineOptions(commandLine.ToString(), _customCanvasParameters));
        }
    }
}