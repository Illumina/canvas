using System;
using System.Collections.Generic;
using System.IO;
using CanvasBin;
using CanvasCommon;
using Isas.SequencingFiles;
using Illumina.Common;
using Xunit;

namespace CanvasTest
{
    public class TestCanvasBin
    {
        public void TestBinOneAlignment(int pos1, int pos2)
        {
            uint qualityThreshold = 3;
            Dictionary<string, int> readNameToBinIndex = new Dictionary<string, int>();
            HashSet<string> samePositionReadNames = new HashSet<string>();
            long usableFragmentCount = 0;
            List<SampleGenomicBin> bins = new List<SampleGenomicBin>()
            {
                new SampleGenomicBin("chr1", 100, 200, 50, 0)
            };
            int binIndexStart = 0;

            BamAlignment alignment1 = new BamAlignment();
            BamAlignment alignment2 = new BamAlignment();
            alignment1.Name = alignment2.Name = "ReadName";
            alignment1.AlignmentFlag = 0x1 | 0x2;
            alignment2.AlignmentFlag = 0x1 | 0x2;
            alignment1.Position = pos1;
            alignment1.MatePosition = pos2;
            alignment1.FragmentLength = 100;
            alignment2.Position = pos2;
            alignment2.MatePosition = pos1;
            alignment2.FragmentLength = -100;
            alignment1.MapQuality = 10;
            alignment2.MapQuality = 10;

            // Both reads pass filters
            FragmentBinner.BinTask.BinOneAlignment(alignment1, qualityThreshold, readNameToBinIndex, samePositionReadNames,
                ref usableFragmentCount, bins, ref binIndexStart);
            FragmentBinner.BinTask.BinOneAlignment(alignment2, qualityThreshold, readNameToBinIndex, samePositionReadNames,
                ref usableFragmentCount, bins, ref binIndexStart);
            Assert.Equal(bins[0].Count, 1);

            // First read passes filters
            bins[0].Count = 0; // reset bin count
            alignment2.MapQuality = 2; // below quality threshold of 3
            FragmentBinner.BinTask.BinOneAlignment(alignment1, qualityThreshold, readNameToBinIndex, samePositionReadNames,
                ref usableFragmentCount, bins, ref binIndexStart);
            FragmentBinner.BinTask.BinOneAlignment(alignment2, qualityThreshold, readNameToBinIndex, samePositionReadNames,
                ref usableFragmentCount, bins, ref binIndexStart);
            Assert.Equal(bins[0].Count, 0);

            // Second read passes filters
            bins[0].Count = 0; // reset bin count
            alignment1.MapQuality = 2; // below quality threshold of 3
            alignment2.MapQuality = 10;
            FragmentBinner.BinTask.BinOneAlignment(alignment1, qualityThreshold, readNameToBinIndex, samePositionReadNames,
                ref usableFragmentCount, bins, ref binIndexStart);
            FragmentBinner.BinTask.BinOneAlignment(alignment2, qualityThreshold, readNameToBinIndex, samePositionReadNames,
                ref usableFragmentCount, bins, ref binIndexStart);
            Assert.Equal(bins[0].Count, 0);

            // Both fail filters
            bins[0].Count = 0; // reset bin count
            alignment1.MapQuality = 2; // below quality threshold of 3
            alignment2.MapQuality = 2; // below quality threshold of 3
            FragmentBinner.BinTask.BinOneAlignment(alignment1, qualityThreshold, readNameToBinIndex, samePositionReadNames,
                ref usableFragmentCount, bins, ref binIndexStart);
            FragmentBinner.BinTask.BinOneAlignment(alignment2, qualityThreshold, readNameToBinIndex, samePositionReadNames,
                ref usableFragmentCount, bins, ref binIndexStart);
            Assert.Equal(bins[0].Count, 0);
        }

        [Fact]
        public void TestBinOneAlignment()
        {
            TestBinOneAlignment(100, 120);
        }

        [Fact]
        public void TestBinOneAlignmentSamePosition()
        {
            TestBinOneAlignment(100, 100);
        }

        [Fact]
        public void TestBinSingleEndBam()
        {
            string assemblyFolder = Isas.Framework.Utilities.Utilities.GetAssemblyFolder(typeof(TestCanvasBin));
            string dataFolder = Path.Combine(assemblyFolder, "Data");
            string bedPath = Path.Combine(dataFolder, "bins_chrM.bed");
            string bamPath = Path.Combine(dataFolder, "single-end.bam");
            Dictionary<string, List<SampleGenomicBin>> bins = CanvasCommon.Utilities.LoadBedFile(bedPath, gcIndex: 3);
            string chrom = "chrM";
            FragmentBinner.BinTask binTask = new FragmentBinner.BinTask(null, chrom, bamPath, bins[chrom]);
            bool exceptionCaught = false;
            try
            {
                binTask.DoIt();
            }
            catch (IlluminaException e)
            {
                if (e.Message.Contains("No paired alignments found"))
                    exceptionCaught = true;
            }
            Assert.True(exceptionCaught);
        }

        [Fact]
        public void TestAllChromsInBedAreInBam()
        {
            CanvasBinParameters parameters = new CanvasBinParameters();
            string assemblyFolder = Isas.Framework.Utilities.Utilities.GetAssemblyFolder(typeof(TestCanvasBin));
            string dataFolder = Path.Combine(assemblyFolder, "Data");
            parameters.predefinedBinsFile = Path.Combine(dataFolder, "bins_chrU.bed");
            parameters.bamFile = Path.Combine(dataFolder, "single-end.bam");
            parameters.isPairedEnd = true;

            FragmentBinner fragmentBinner = new FragmentBinner(parameters);
            bool exceptionCaught = false;
            try
            {
                fragmentBinner.Bin();
            }
            catch (IlluminaException e)
            {
                if (e.Message.Contains(String.Format("Not all chromosomes in {0} are found in {1}.", parameters.predefinedBinsFile, parameters.bamFile)))
                    exceptionCaught = true;
            }
            Assert.True(exceptionCaught);
        }
    }
}
