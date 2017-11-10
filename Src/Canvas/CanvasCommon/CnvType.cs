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
        public const string CnvTag = "CNV";
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
                    return CnvTag;
                case CnvType.LossOfHeterozygosity:
                    return "LOH";
                default:
                    throw new Illumina.Common.IlluminaException($"SVTYPE field is unsupported for CNV type: {cnvType}");
            }
        }

        public static bool IsReference(CnvType cnvType) => cnvType.Equals(CnvType.Reference);

    }
}