namespace CanvasPedigreeCaller
{
    class PedigreeCallerParameters
    {
        public int MaximumCopyNumber { get; set; } = 5;
        public int MaxAlleleNumber { get; set; } = 3;
        public int DefaultAlleleDensityThreshold { get; set; } = 1000;
        public double MaxQscore { get; set; } = 60.0;
        public int DefaultPerSegmentAlleleMaxCounts { get; set; } = 100;
        public int DefaultAlleleCountThreshold { get; set; } = 4;
        public int MaxNumOffspringGenotypes { get; set; } = 500;
        public double DeNovoRate { get; set; } = 0.00001;
        public int MinimumCallSize { get; set; } = 1000;
    }
}