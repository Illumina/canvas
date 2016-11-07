using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using CanvasCommon;

namespace CanvasPedigreeCaller
{
    public class PedigreeMember
    {
        public enum Kinship
        {
            Parent, Offspring
        }
        public List<CanvasSegment> Segments = new List<CanvasSegment>();
        public double MeanCoverage { get; set; }
        public double MeanMafCoverage { get; set; }
        public string Name { get; set; }
        public List<string> Parents { get; set; }
        public List<string> Offspring { get; set; }
        public PloidyInfo Ploidy { get; set; }
        public double Variance { get; internal set; }
        public double MafVariance { get; internal set; }

        public CopyNumberModel CnModel { get; set; }
        public Kinship Kin { get; set; }

        public double GetCoverage(int segmentIndex)
        {
            return Segments[segmentIndex].MedianCount;
        }
        public Tuple<int,int> GetAlleleCounts(int segmentIndex)
        {
            double allele1 = Segments[segmentIndex].Alleles.Counts.Select(x=>x.Item1).Average();
            double allele2 = Segments[segmentIndex].Alleles.Counts.Select(x => x.Item2).Average();
            return new Tuple<int, int>(Convert.ToInt32(allele1), Convert.ToInt32(allele2));
        }
    }
}