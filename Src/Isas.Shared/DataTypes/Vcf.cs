using System;
using System.Collections.Generic;
using System.Linq;
using System.Security.Cryptography;
using System.Text;
using System.Threading.Tasks;

namespace Isas.Shared.DataTypes
{
    public class Vcf
    {
        private const string Extension = ".vcf.gz";
        public IFileLocation VcfFile { get; }

        public IFileLocation TabixIndex => VcfFile.AppendName(".tbi");
        public Vcf() { }
        public Vcf(IFileLocation vcfFile)
        {
            if (vcfFile == null)
                throw new ArgumentNullException(nameof(vcfFile));
            if (!vcfFile.Name.EndsWith(Extension))
                throw new ArgumentException($"Invalid vcf path {vcfFile}. Path must end with .vcf.gz");
            VcfFile = vcfFile;
        }

        public Vcf MoveWithStub(IFileLocation newStub, Action<IFileLocation, IFileLocation> move)
        {
            var newVcf = GetVcfFromStub(newStub);
            Move(newVcf, move);
            return newVcf;
        }

        public Vcf Move(IFileLocation newPath, Action<IFileLocation, IFileLocation> move)
        {
            var newVcf = new Vcf(newPath);
            Move(newVcf, move);
            return newVcf;
        }

        public void Move(Vcf destination, Action<IFileLocation, IFileLocation> move)
        {
            move(VcfFile, destination.VcfFile);
            if (TabixIndex.Exists)
                move(TabixIndex, destination.TabixIndex);
        }

        public void Delete()
        {
            VcfFile.Delete();
            TabixIndex.Delete();
        }

        public static Vcf GetVcfFromStub(IFileLocation stub)
        {
            return new Vcf(stub.AppendName(Extension));
        }

        public IFileLocation GetStub()
        {
            return VcfFile.ReplaceEnd(Extension, "");
        }
    }
}
