using System;
using System.Collections.Generic;
using System.Text;
using System.IO;
using System.Linq;
using System.Linq.Expressions;
using System.Xml.Schema;
using CanvasCommon;
using MathNet.Numerics.Distributions;

namespace CanvasPedigreeCaller
{
    class CanvasPedigreeCaller
    {
        #region Members
        // Static:
        static private int MaximumCopyNumber = 10;
        public int QualityFilterThreshold { get; set; } = 10;
        #endregion


        internal int CallVariants(List<string> variantFrequencyFiles, List<string> segmentFiles, string outDir, List<string> ploidyBedPaths, string referenceFolder, List<string> sampleNames, string truthDataPath)
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

            Pedigree<PedigreeMember> pedigree= new Pedigree<PedigreeMember>(pedigreeMembers[0], true);
            pedigree.AddParent(pedigreeMembers[1]);
            pedigree.AddChildren(pedigreeMembers[2]);

            var numberOfSegments = pedigreeMembers[0].Segments.Count;
            var segmentIndex = 0;

            List<PedigreeMember> parents = pedigree.GetParents();
            List<PedigreeMember> offsprings = pedigree.GetChildren();
            double[][] transitionMatrix = GetTransitionMatrix(numberOfSegments);

            foreach (PedigreeMember parent in parents)
                parent.CnModel = new CopyNumberModel(numberOfSegments, parent.MeanCoverage, parent.Variance);
            foreach (PedigreeMember offspring in offsprings)
                offspring.CnModel = new CopyNumberModel(numberOfSegments, offspring.MeanCoverage, offspring.Variance);
            List<List<int>> copyNumbers = new List<List<int>>();
            List<double> lls = new List<double>();

            while (segmentIndex < numberOfSegments)
            {
                var segmentCn = new List<int>();
                var ll = MaximalLikelyhood(parents, offsprings, segmentIndex, ref segmentCn, transitionMatrix);
                copyNumbers.Add(segmentCn);
                lls.Add(ll);
                segmentIndex++;
            }

            for (int pedigreeMemberIndex = 0; pedigreeMemberIndex < pedigreeMembers.Count; pedigreeMemberIndex++)
            {             
                for (segmentIndex = 0; segmentIndex < pedigreeMembers[pedigreeMemberIndex].Segments.Count; segmentIndex++)
                {
                    pedigreeMembers[pedigreeMemberIndex].Segments[segmentIndex].CopyNumber = copyNumbers[pedigreeMemberIndex][segmentIndex];
                }
                var outFile = Path.Combine(outDir, pedigreeMembers[pedigreeMemberIndex].Name + "CNV.vcf");
                CanvasSegment.WriteCoveragePlotData(pedigreeMembers[pedigreeMemberIndex].Segments, pedigreeMembers[pedigreeMemberIndex].MeanCoverage * 2.0, pedigreeMembers[pedigreeMemberIndex].Ploidy, outFile, referenceFolder);
            }
            return 0;
        }


        public static int AggregateVariantCoverage(ref List<CanvasSegment> segments)
        {
            var variantCoverage = segments.SelectMany(segment => segment.VariantTotalCoverage).ToList();
            return variantCoverage.Any() ? Utilities.Median(variantCoverage) : 0;
        }

        double[][] GetTransitionMatrix(int numCnStates)
        {
            double[][] transitionMatrix = CanvasCommon.Utilities.MatrixCreate(numCnStates, numCnStates);
            for (int cn = 0; cn < numCnStates; cn++)
            {
                var gtLikelyhood = new Poisson(Math.Round(cn / 2.0));
                for (int gt = 0; gt < cn; gt++)
                    transitionMatrix[cn][gt] = gtLikelyhood.Probability(gt);
            }
            return transitionMatrix;
        }

        double MaximalLikelyhood(List<PedigreeMember> parents, List<PedigreeMember> children, int segmentPosition, ref List<int> bestCnStates, double[][] transitionMatrix)
        {
            double maximalLikelyhood = Double.MinValue;
            double marginals = 0;
            List<int> currentCnStates = new List<int>();
            var parent1Likelihood = parents[0].CnModel.GetCnLikelihood(parents[0].GetCoverage(segmentPosition));
            var parent2Likelihood = parents[1].CnModel.GetCnLikelihood(parents[1].GetCoverage(segmentPosition));
            for (int cn1 = 0; cn1 < parent1Likelihood.Count; cn1++)
            {
                currentCnStates.Add(cn1);
                double currentLikelihood = parent1Likelihood[cn1];
                for (int cn2 = 0; cn2 < parent2Likelihood.Count; cn1++)
                {
                    currentCnStates.Add(cn2);
                    currentLikelihood *= parent2Likelihood[cn2];
                    for (int gt1 = 0; gt1 < cn1; gt1++)
                    {
                        for (int gt2 = 0; gt2 < cn2; gt2++)
                        {
                            foreach (PedigreeMember child in children)
                            {
                                currentCnStates.Add(gt1 + gt2);
                                currentLikelihood *= parent1Likelihood[cn1] * parent2Likelihood[cn2] * transitionMatrix[cn1][gt1] * transitionMatrix[cn2][gt2] *
                                    child.CnModel.GetCnLikelihood(child.GetCoverage(segmentPosition))[gt1 + gt2];
                            }
                            if (currentLikelihood > maximalLikelyhood)
                            {
                                maximalLikelyhood = currentLikelihood;
                                bestCnStates = currentCnStates;
                            }
                            currentCnStates = new List<int>();
                            marginals += maximalLikelyhood;
                        }
                    }
                }
            }
            return maximalLikelyhood / marginals;
        }
    }
}
