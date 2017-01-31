using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using Illumina.Common.FileSystem;
using Illumina.SecondaryAnalysis.VariantCalling;
using Isas.Framework.DataTypes;
using Isas.Manifests.NexteraManifest;
using Isas.SequencingFiles;

namespace Canvas.Wrapper
{
    public interface ICanvasEnrichmentInput : ICanvasCheckpointInput
    {
        NexteraManifest NexteraManifest { get; }
    }

    public class CanvasEnrichmentInput : ICanvasEnrichmentInput
    {
        public Bam Bam { get; }
        public GenomeMetadata GenomeMetadata { get; }
        public IEnumerable<Bam> NormalBamPaths { get; }
        public NexteraManifest NexteraManifest { get; }
        public CanvasEnrichmentPrecomputedControl PrecomputedControl { get; }
        public SamplePloidyInfo PloidyInfo { get; }
        public IFileLocation PredefinedBinsFile { get; }
        public CanvasPcaModels PcaModels { get; }
        public bool IsCanvasNormalizePcaMode => PcaModels != null && (PcaModels.MaleModelFile != null || PcaModels.FemaleModelFile != null);

        public CanvasEnrichmentInput(Bam bam, GenomeMetadata genomeMetadata,
            IEnumerable<Bam> controlBamPaths,
            NexteraManifest nexteraManifest,
            CanvasEnrichmentPrecomputedControl precomputedControl,
            SamplePloidyInfo ploidyInfo,
            IFileLocation predefinedBinsFile,
            CanvasPcaModels pcaModels)
        {
            Bam = bam;
            GenomeMetadata = genomeMetadata;
            NexteraManifest = nexteraManifest;
            PrecomputedControl = precomputedControl;
            NormalBamPaths = new ReadOnlyCollection<Bam>(controlBamPaths.ToList());
            PloidyInfo = ploidyInfo;
            PredefinedBinsFile = predefinedBinsFile;
            PcaModels = pcaModels;
        }
    }

    public class CanvasEnrichmentPrecomputedControl
    {
        public int BinSize { get; }
        public IFileLocation BinnedPath { get; }
        public string SexChromosomeKaryotype { get; }

        public CanvasEnrichmentPrecomputedControl(IFileLocation binnedPath, int binSize, string sexChromosomeKaryotype)
        {
            BinSize = binSize;
            BinnedPath = binnedPath;
            SexChromosomeKaryotype = sexChromosomeKaryotype;
        }
    }

    public class CanvasPcaModels
    {
        public readonly IFileLocation MaleModelFile;
        public readonly IFileLocation FemaleModelFile;

        public CanvasPcaModels(IFileLocation maleModelFile, IFileLocation femaleModelFile)
        {
            MaleModelFile = maleModelFile;
            FemaleModelFile = femaleModelFile;
        }
    }
}