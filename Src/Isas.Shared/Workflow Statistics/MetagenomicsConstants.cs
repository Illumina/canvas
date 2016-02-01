namespace Isas.Shared
{
    public static class MetagenomicsConstants
    {
        public const int BootstrapIterationCount = 20;
        public const int NumTaxonomicLevels      = 8;
        public const int SeedSize                = 8;

        public const int KingdomLevel = 0;
        public const int PhylumLevel  = 1;
        public const int ClassLevel   = 2;
        public const int OrderLevel   = 3;
        public const int FamilyLevel  = 4;
        public const int GenusLevel   = 5;
        public const int SpeciesLevel = 6;
        public const int SubtypeLevel = 7;

        public static string[] TaxonomicNames = new[] { "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Subtype" };

        public const string UnknownTaxonomy        = "(unknown)";
        public const int MinSeedsForClassification = 42;
    }
}
