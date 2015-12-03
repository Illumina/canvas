using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace SequencingFiles
{
    public enum NovelVariantAnnotator
    {
        None = 0,
        VEP,
        Nirvana
    }

    /// <summary>
    /// Parameters to pass in to our gvcf annotator/splitter.
    /// </summary>
    public class AnnotationParameters
    {
        public string[] InputPaths; // Input gvcf (or variants-only vcf) paths.
        public string[] OutputPaths; // Output paths.
        public string[] OutputVCFPaths; // Optional.  If provided, we output JUST the variants to this file.
        public string AnnotationDatabase;
        public bool VerboseLogging;
        public string TranscriptSource = "refseq";
        public string AnnotationDatabaseVersion = "72.5";
        public bool SkipMito = false;
        public string RefMinorAlleleBedPath = null;
        public string GenomeXMLPath;
        public NovelVariantAnnotator NovelAnnotator;
        public bool GenerateAntFile = true;
    }
}
