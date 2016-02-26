namespace SequencingFiles.Vcf
{
    public static class VcfCommon
    {
        #region header titles

        public const string ChromosomeHeader = "#CHROM";
        public const string FilterTag = "##FILTER=";
        public const string FormatTag = "##FORMAT=";
        public const string InfoTag = "##INFO=";
        public const string SourceTag = "##source=";

        #endregion

        #region storage classes

        private enum NumberClasses
        {
            None,
            Period,
            Single,
            Allele,
            Genotype
        }

        private enum StorageClasses
        {
            Float,
            Integer,
            String
        }

        #endregion

        #region column indexes

        public const int ChromIndex = 0;
        public const int PosIndex = 1;
        public const int IDIndex = 2;
        public const int RefIndex = 3;
        public const int AltIndex = 4;
        public const int QualIndex = 5;
        public const int FilterIndex = 6;
        public const int InfoIndex = 7;
        public const int FormatIndex = 8;
        public const int GenotypeIndex = 9;

        #endregion

        private const int MinNumColumns = 9;

        public static class GenotypeColumnFields
        {
            public const string Genotype = "GT";
        }

        public static class InfoColumnFields
        {
            public const string End = "End";
        }
    }
}