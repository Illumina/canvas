using System.Collections.Generic;
using ProtoBuf;

namespace Isas.Shared
{
    [ProtoContract(SkipConstructor = true)]
    public class ProbeSet
    {
        #region Members
        [ProtoMember(1, AsReference = true)]
        public List<ProbeSetTarget> Targets = new List<ProbeSetTarget>();
        [ProtoMember(2)]
        public string LocusID { get; set; } // assay ID in targeted RNA-Seq
        [ProtoMember(3)]
        public string Species { get; set; }
        [ProtoMember(4)]
        public string BuildID { get; set; }
        [ProtoMember(5)]
        public string SubmittedSequence { get; set; }
        [ProtoMember(6)]
        public string Chromosome { get; set; }
        [ProtoMember(7)]
        public int StartPosition { get; set; }
        [ProtoMember(8)]
        public int EndPosition { get; set; }
        [ProtoMember(9)]
        public bool ReverseStrand { get; set; }
        [ProtoMember(10)]
        public bool SubmittedSequenceReverseStrand { get; set; }
        [ProtoMember(11)]
        public string Read1Tag { get; set; } // upstream 
        [ProtoMember(12)]
        public string Read2Tag { get; set; } // downstream
        #endregion
    }
}