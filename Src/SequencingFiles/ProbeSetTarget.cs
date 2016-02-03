using ProtoBuf;
using SequencingFiles;

namespace Isas.Shared
{
    [ProtoContract(SkipConstructor = true)]
    public class ProbeSetTarget
    {
        #region Members
        [ProtoMember(1)]
        public string AmpliconID;
        [ProtoMember(2)]
        public string AssayID;
        [ProtoMember(3)]
        public int ConcatenatedOffset;
        [ProtoMember(4)]
        public string FullAmpliconSequence;
        [ProtoMember(5)]
        public string FullReverseSequence;
        [ProtoMember(6)]
        public string GeneName;
        [ProtoMember(7, AsReference = true)]
        public AmpliconManifest Manifest;
        [ProtoMember(8)]
        public string Name;
        [ProtoMember(9)]
        public CigarAlignment ReferenceAlignment; // optional: used for targets which recreate an indel relative to the reference
        [ProtoMember(10)]
        public int ReferenceIndex; // used for BAM I/O
        [ProtoMember(11)]
        public string Chromosome { get; set; }
        [ProtoMember(12)]
        public int StartPosition { get; set; } // 1-based
        [ProtoMember(13)]
        public int EndPosition { get; set; } // 1-based, inclusive!
        [ProtoMember(14)]
        public bool ReverseStrand { get; set; }
        [ProtoMember(15)]
        public string Sequence { get; set; }
        [ProtoMember(16)]
        public string ReverseSequence { get; set; }
        [ProtoMember(17, AsReference = true)]
        public ProbeSet ProbeA { get; set; }
        [ProtoMember(18, AsReference = true)]
        public ProbeSet ProbeB { get; set; }
        [ProtoMember(19)]
        public int Index { get; set; }
        [ProtoMember(20)]
        public string Build { get; set; }
        [ProtoMember(21)]
        public string Species { get; set; }
        [ProtoMember(22)]
        public string SoftClipReadDirection { get; set; }
        [ProtoMember(23)]
        public int SoftClipUpstream { get; set; }
        [ProtoMember(24)]
        public int SoftClipDownstream { get; set; }
        #endregion
    }
}