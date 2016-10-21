using System;
using System.Collections.Generic;
using System.Text;
using System.IO;
using System.Linq;
using System.Xml.Schema;
using CanvasCommon;

namespace CanvasPedigreeCaller
{
    class CanvasPedigreeCaller
    {
        #region Members
        // Static:
        static private int MaximumCopyNumber = 10;

        // Data:
        List<SegmentPloidy> AllPloidies;
        double DiploidCoverage = 0;

        // Parameters:
        public static double CoverageWeighting = 0.6;
        private double CoverageWeightingFactor; // Computed from CoverageWeighting
        public bool IsDbsnpVcf = false;
        static protected int MinimumVariantFrequenciesForInformativeSegment = 50;
        protected int MedianHetSnpsDistance = 463; // based on NA12878 VFResults.txt.gz file
        CopyNumberOracle CNOracle = null;
        public QualityScoreParameters germlineScoreParameters;
        public int QualityFilterThreshold { get; set; } = 10;

        // File paths:
        public string TempFolder;
        CoverageModel Model;
        Pedigree Pedigree;

        #endregion


        internal int CallVariants(List<string> variantFrequencyFiles, List<string> segmentFiles, string outFile, List<string> ploidyBedPaths, string referenceFolder, List<string> sampleNames, string truthDataPath)
        {
            // load files
            // initialize data structures and classes
            int fileCounter = 0;
            List<PedigreeMember> pedigreeMembers = new List<PedigreeMember>();
            foreach (string sampleName in sampleNames)
            {
                PedigreeMember pedigreeMember = new PedigreeMember();
                pedigreeMember.Name = sampleName;
                pedigreeMember.Segments = CanvasSegment.ReadSegments(segmentFiles[fileCounter]);
                pedigreeMember.MeanCoverage = CanvasIO.LoadVariantFrequencies(variantFrequencyFiles[fileCounter], pedigreeMember.Segments);
                pedigreeMember.Ploidy = PloidyInfo.LoadPloidyFromBedFile(ploidyBedPaths[fileCounter]);
                pedigreeMembers.Add(pedigreeMember);
            }
            Pedigree = new Pedigree(pedigreeMembers);
            Pedigree parents = Pedigree.GetSubset(PedigreeMember.PedigreeMemberType.Parent);
            Pedigree offspring = Pedigree.GetSubset(PedigreeMember.PedigreeMemberType.Parent);
            likelihood = parents.CopyNumberLikelihood();
            offspring.CopyNumberLikelihood(likelihood);
            return 0;
        }

        public static int AggregateVariantCoverage(ref List<CanvasSegment> segments)
        {
            var variantCoverage = segments.SelectMany(segment => segment.VariantTotalCoverage).ToList();
            return variantCoverage.Any() ? Utilities.Median(variantCoverage) : 0;
        }


    }
}
