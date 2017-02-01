using System.Collections.Generic;
using Illumina.Common.FileSystem;
using Isas.Framework.DataTypes;
using Isas.Manifests.NexteraManifest;

namespace Canvas.Wrapper
{
    public class CanvasEnrichmentInputCreator<TCanvasEnrichmentInput> where TCanvasEnrichmentInput : ICanvasEnrichmentInput
    {
        public Dictionary<NexteraManifest, IFileLocation> WriteManifests(SampleSet<TCanvasEnrichmentInput> inputs, IDirectoryLocation sandbox)
        {
            var manifests = new Dictionary<NexteraManifest, IFileLocation>();
            foreach (var input in inputs.SampleData)
            {
                var manifest = input.NexteraManifest;
                if (manifests.ContainsKey(manifest)) continue;
                manifests[manifest] = WriteManifest(manifest, sandbox);
            }
            return manifests;
        }

        private IFileLocation WriteManifest(NexteraManifest manifest, IDirectoryLocation sandbox)
        {
            var path = sandbox.GetFileLocation(manifest.Name);
            NexteraManifestUtils.WriteNexteraManifests(manifest, path.FullName);
            return path;
        }

        public Dictionary<NexteraManifest, IFileLocation> CreateDbSnpVcfForManifests(SampleSet<TCanvasEnrichmentInput> inputs, IDirectoryLocation sandBox, ICanvasAnnotationFileProvider annotationFileProvider)
        {
            var dbSnpVcfs = new Dictionary<NexteraManifest, IFileLocation>();
            foreach (var input in inputs.SampleData)
            {
                var manifest = input.NexteraManifest;
                if (dbSnpVcfs.ContainsKey(manifest)) continue;
                var fullDbSnpVcf = annotationFileProvider.GetDbSnpVcf(input.GenomeMetadata);
                dbSnpVcfs[manifest] = CreateDbSnpVcfForManifest(fullDbSnpVcf, manifest, sandBox);
            }
            return dbSnpVcfs;
        }

        private IFileLocation CreateDbSnpVcfForManifest(IFileLocation fullDbSnpVcf, NexteraManifest manifest, IDirectoryLocation sandBox)
        {
            IFileLocation targetedDbSnpVcf = sandBox.GetFileLocation($"{manifest.Name}_{fullDbSnpVcf.Name}");
            Isas.Manifests.NexteraManifest.VcfUtilities.IntersectVcfWithManifest(fullDbSnpVcf.FullName, targetedDbSnpVcf.FullName, manifest);
            return targetedDbSnpVcf;
        }
    }
}