using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using CanvasCommon;
using Illumina.Common;
using Isas.Framework.DataTypes;

namespace CanvasPedigreeCaller
{
    public class PedigreeMember
    {
        public PedigreeMember(SampleId id, Kinship kin)
        {
            Kin = kin;
            Id = id;
        }
        public enum Kinship
        {
            Other, Parent, Proband
        }
        public List<CanvasSegmentsSet> SegmentSets = new List<CanvasSegmentsSet>(); 
        public SampleId Id { get; }
        public string Name => Id.ToString();
        public Kinship Kin { get; }


        public double GetCoverage(int setPosition, int segmentPosition, SegmentsSet segmentsSet, int numberOfTrimmedBins)
        {
            return SegmentSets[setPosition].GetSet(segmentsSet)[segmentPosition].MedianCount;
        }

        public List<Tuple<int, int>> GetAlleleCounts(int setPosition, int segmentPosition, SegmentsSet segmentsSet)
        {
            return SegmentSets[setPosition].GetSet(segmentsSet)[segmentPosition].Balleles.GetAlleleCounts();
        }

        public List<CanvasSegment> GetCanvasSegments()
        {
            return SegmentSets.SelectMany(set => set.GetSet(set.SelectedSet)).ToList();
        }
    }
}