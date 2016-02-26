
namespace Isas.Shared
{
    public class BclRunFolder
    {
        public readonly IDirectoryLocation RunFolder;
        public readonly Adapters Adapters;
        public readonly ReadStructure ReadStructure;
        public readonly string FlowcellId;
        public readonly int FlowcellNumber;

        // one dictionary entry per lane for this runfolder
        // if an entry has an empty set of tiles, then all tiles from that lane will be used
        public readonly TileSelection TileSelection;
        private readonly BclFormat _bclFormat;

        public BclRunFolder(IDirectoryLocation runFolder, string flowcellId, int flowcellNumber, Adapters adapters, ReadStructure readStructure, TileSelection tileSelection, BclFormat bclFormat)
        {
            RunFolder = runFolder;
            Adapters = adapters;
            ReadStructure = readStructure;
            TileSelection = tileSelection;
            _bclFormat = bclFormat;
            FlowcellId = flowcellId;
            FlowcellNumber = flowcellNumber;
        }

        public override bool Equals(object obj)
        {
            BclRunFolder runFolder = obj as BclRunFolder;
            if (runFolder == null)
                return false;
            return RunFolder.Equals(runFolder.RunFolder);
        }

        public override int GetHashCode()
        {
            return RunFolder.FullName.GetHashCode();
        }
    }

    public enum BclFormat
    {
        BGZF, // nova .bcl.bgzf
        GZ, // HiSeq .bcl.gz
        UNCOMPRESSED //MiSeq .bcl
    }
}
