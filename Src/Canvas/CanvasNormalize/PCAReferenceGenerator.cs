using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using CanvasCommon;
using Isas.SequencingFiles;
using Isas.Shared.Utilities.FileSystem;

namespace CanvasNormalize
{
    public class PCAReferenceGenerator : IReferenceGenerator
    {
        private readonly IFileLocation _sampleBinnedFile;
        private readonly PCAModel _model;
        private readonly NexteraManifest _manifest;
        private readonly double _minBinCount;
        private readonly double _maxBinCount;
        private readonly RawRatioCalculator _ratioCalculator;

        public PCAReferenceGenerator(IFileLocation sampleBinnedFile, IFileLocation pcaModelFile,
            NexteraManifest manifest, double minBinCount = 1, double maxBinCount = double.PositiveInfinity)
        {
            _sampleBinnedFile = sampleBinnedFile;
            _model = new PCAModel(pcaModelFile);
            _manifest = manifest;
            _minBinCount = minBinCount;
            _maxBinCount = maxBinCount;
            _ratioCalculator = new RawRatioCalculator(manifest, _minBinCount, _maxBinCount);
        }

        public void Run(IFileLocation outputFile)
        {
            List<SampleGenomicBin> sampleBins = CanvasIO.ReadFromTextFile(_sampleBinnedFile.FullName);
            VerifyBinOrder(sampleBins);

            // set bin count to 1 if less than 1
            foreach (var bin in sampleBins)
                bin.Count = Math.Max(1, bin.Count);

            // center the sample
            var centeredSampleVector = Enumerable.Zip(sampleBins, _model.Mu, (bin, mu) => (double)bin.Count - mu.Count).ToArray();

            // project onto the axes
            var projectedSampleVector = CanvasCommon.Utilities.Project(centeredSampleVector, _model.Axes);

            // undo centering and set bin count to 1 if less than 1
            var referenceVector = Enumerable.Zip(_model.Mu, projectedSampleVector, (bin, count) => Math.Max(1, bin.Count + count));

            // write temporary reference count file
            var tempReferenceFile = new FileLocation(Path.GetTempFileName());
            var tempReferenceBins = Enumerable.Zip(sampleBins, referenceVector,
                (bin, count) => new SampleGenomicBin(bin.GenomicBin.Chromosome, bin.Start, bin.Stop, bin.GenomicBin.GC, (float)count));
            CanvasIO.WriteToTextFile(tempReferenceFile.FullName, tempReferenceBins);

            // calcualte median ratio
            var ratios = new BinCounts(_ratioCalculator.Run(_sampleBinnedFile, tempReferenceFile), _manifest);
            double medianRatio = ratios.OnTargetMedianBinCount;

            // delete temporary reference count file
            Isas.Shared.Utilities.Utilities.SafeDelete(tempReferenceFile.FullName);

            // multiply reference counts by the median ratio
            var referenceBins = Enumerable.Zip(sampleBins, referenceVector,
                (bin, count) => new SampleGenomicBin(bin.GenomicBin.Chromosome, bin.Start, bin.Stop, bin.GenomicBin.GC, (float)(count * medianRatio)));

            // write reference count file
            CanvasIO.WriteToTextFile(outputFile.FullName, referenceBins);
        }

        /// <summary>
        /// 
        /// </summary>
        protected void VerifyBinOrder(IEnumerable<SampleGenomicBin> bins)
        {
            bool mismatch = Enumerable.Zip(bins, _model.Mu, (bin1, bin2) => bin1.IsSameBin(bin2))
                .SkipWhile(sameBin => sameBin).TakeWhile(sameBin => !sameBin).Any();

            if (mismatch)
                throw new ApplicationException("Bins must be in the same order as those in the model file.");
        }

        class PCAModel
        {
            public readonly SampleGenomicBin[] Mu;
            public readonly double[][] Axes;

            public PCAModel(IFileLocation pcaModelFile)
            {
                List<SampleGenomicBin> mu;
                List<double[]> axes;
                LoadModel(pcaModelFile, out mu, out axes);
                Mu = mu.ToArray();
                Axes = axes.Select(a => a.ToArray()).ToArray();
            }

            public PCAModel(SampleGenomicBin[] mu, double[][] axes)
            {
                Mu = mu;
                Axes = axes;
            }

            private static void LoadModel(IFileLocation modelFile, out List<SampleGenomicBin> mu, out List<double[]> axes)
            {
                mu = new List<SampleGenomicBin>();
                axes = new List<double[]>();
                List<List<double>> tempAxes = new List<List<double>>();

                using (GzipReader reader = new GzipReader(modelFile.FullName))
                {
                    string line = reader.ReadLine();
                    for (int i = 0; i < line.Split('\t').Length - 4; i++) // initialize axes
                        tempAxes.Add(new List<double>());

                    while (line != null)
                    {
                        string[] toks = line.Split('\t');
                        string chrom = toks[0];
                        int start = int.Parse(toks[1]);
                        int stop = int.Parse(toks[2]);
                        float mean = float.Parse(toks[3]);
                        mu.Add(new SampleGenomicBin(chrom, start, stop, -1, mean));
                        for (int i = 0; i < tempAxes.Count; i++)
                        {
                            tempAxes[i].Add(double.Parse(toks[i + 4]));
                        }
                        line = reader.ReadLine();
                    }
                }

                foreach (var axis in tempAxes)
                {
                    axes.Add(CanvasCommon.Utilities.NormalizeBy2Norm(axis.ToArray()));
                }

                if (!AreOrthogonal(axes))
                    throw new ApplicationException(String.Format("Axes are not orthogonal to each other in {0}.",
                        modelFile.FullName));
            }

            private static bool AreOrthogonal(List<double[]> axes)
            {
                for (int i = 0; i < axes.Count; i++)
                {
                    for (int j = i + 1; j < axes.Count; j++)
                    {
                        if (!CanvasCommon.Utilities.AreOrthogonal(axes[i], axes[j]))
                            return false;
                    }
                }
                return true;
            }
        }
    }
}
