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


        public double GetCoverage(CanvasSegmentIndex segmentIndex, int numberOfTrimmedBins)
        {
            return GetCanvasSegment(segmentIndex).MedianCount;
        }

        public List<Tuple<int, int>> GetAlleleCounts(CanvasSegmentIndex segmentIndex)
        {
            return GetCanvasSegment(segmentIndex).Balleles.GetAlleleCounts();
        }

        public List<CanvasSegment> GetCanvasSegments()
        {
            return SegmentSets.SelectMany(set => set.GetSet(set.SelectedSet)).ToList();
        }
        public CanvasSegment GetCanvasSegment(CanvasSegmentIndex segmentIndex)
        {
            return SegmentSets[segmentIndex.SetPosition].GetSet(segmentIndex.Set)[segmentIndex.SegmentPosition];
        }
    }
}