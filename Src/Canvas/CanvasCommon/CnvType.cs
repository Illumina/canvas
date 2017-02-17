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
                default:
                    throw new ApplicationException($"Unsupported CNV type: {cnvType}");
            }
        }

        public static string ToSvType(this CnvType cnvType)
        {
            switch (cnvType)
            {
                case CnvType.Gain:
                case CnvType.Loss:
                    return "CNV";
                case CnvType.LossOfHeterozygosity:
                    return "LOH";
                default:
                    throw new ApplicationException($"SVTYPE field is unsupported for CNV type: {cnvType}");
            }
        }

        public static string ToAltId(this CnvType cnvType)
        {
            switch (cnvType)
            {
                case CnvType.Gain:
                case CnvType.Loss:
                case CnvType.LossOfHeterozygosity:
                    return "<CNV>";
                case CnvType.Reference:
                    return ".";
                default:
                    throw new ApplicationException($"ALT is unsupported for CNV type: {cnvType}");
            }
        }
    }
}