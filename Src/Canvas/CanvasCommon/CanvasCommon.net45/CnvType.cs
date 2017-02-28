using System;

namespace CanvasCommon
{
    public enum CnvType
    {
        Gain,
        Loss,
        LossOfHeterozygosity,
        Reference,
        ComplexCnv,
    }

    public static class CnvTypeExtensions
    {
        public static string ToVcfId(this CnvType cnvType)
        {
            switch (cnvType)
            {
                case CnvType.Gain:
                    return "GAIN";
                case CnvType.Loss:
                    return "LOSS";
                case CnvType.Reference:
                    return "REF";
                case CnvType.LossOfHeterozygosity:
                    return "LOH";
                case CnvType.ComplexCnv:
                    return "COMPLEXCNV";                   
                default:
                    throw new Illumina.Common.IlluminaException($"Unsupported CNV type: {cnvType}");
            }
        }

        public static string ToSvType(this CnvType cnvType)
        {
            switch (cnvType)
            {
                case CnvType.Gain:
                case CnvType.Loss:
                case CnvType.ComplexCnv:
                    return "CNV";
                case CnvType.LossOfHeterozygosity:
                    return "LOH";
                default:
                    throw new Illumina.Common.IlluminaException($"SVTYPE field is unsupported for CNV type: {cnvType}");
            }
        }

        public static string ToAltId(this CnvType cnvType)
        {
            switch (cnvType)
            {
                case CnvType.Gain:
                case CnvType.Loss:
                case CnvType.LossOfHeterozygosity:
                case CnvType.ComplexCnv:
                    return "<CNV>";
                case CnvType.Reference:
                    return ".";
                default:
                    throw new Illumina.Common.IlluminaException($"ALT is unsupported for CNV type: {cnvType}");
            }
        }
    }
}