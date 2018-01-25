using System;
using System.Collections.Generic;
using System.Text;

namespace EvaluateCNV
{
    public class MetricsCalculator
    {
        // Compute overall stats:

        public long TotalBases { get; }
        public long TotalBasesRight { get; }
        public long TotalBasesRightDirection { get; }
        public long IsGainBases { get; }
        public long CallGainBases { get; }
        public long IsGainBasesCorrect { get; }
        public long IsGainBasesCorrectDirection { get; }
        public long IsLossBases { get; }
        public long CallLossBases { get; }
        public long IsLossBasesCorrect { get; }
        public long IsLossBasesCorrectDirection { get; }
        public long RoiBases { get; }
        public long RoiBasesCorrect { get; }
        public long RoiBasesCorrectDirection { get; }

        private double FractionalPrecision => (IsGainBasesCorrect + IsLossBasesCorrect) / (double)(CallGainBases + CallLossBases);
        private double FractionalRecall => (IsGainBasesCorrect + IsLossBasesCorrect) / (double)(IsGainBases + IsLossBases);
        public double Precision => FractionalPrecision * 100;
        public double Recall => FractionalRecall * 100;
        public double F1Score => 2 * FractionalPrecision * FractionalRecall / (FractionalPrecision + FractionalRecall);
        public double Accuracy => TotalBasesRight / (double)TotalBases;
        // SK: I felt the direction based performance metrices make more sense
        public double DirectionAccuracy => 100 * TotalBasesRightDirection / (double)TotalBases;

        public double DirectionRecall => 100 * (IsGainBasesCorrectDirection + IsLossBasesCorrectDirection) / (double)(IsGainBases + IsLossBases);
        public double DirectionPrecision => 100 * (IsGainBasesCorrectDirection + IsLossBasesCorrectDirection) / (double)(CallGainBases + CallLossBases);
        public double GainRecall => 100 * IsGainBasesCorrect / (double)IsGainBases;
        public double GainDirectionRecall => 100 * IsGainBasesCorrectDirection / (double)IsGainBases;
        public double GainPrecision => 100 * IsGainBasesCorrect / (double)CallGainBases;
        public double GainDirectionPrecision => 100 * IsGainBasesCorrectDirection / (double)CallGainBases;
        public double LossRecall => 100 * IsLossBasesCorrect / (double)IsLossBases;
        public double LossDirectionRecall => 100 * IsLossBasesCorrectDirection / (double)IsLossBases;
        public double LossPrecision => 100 * IsLossBasesCorrect / (double)CallLossBases;
        public double LossDirectionPrecision => 100 * IsLossBasesCorrectDirection / (double)CallLossBases;
        public double ROIAccuracy => 100 * RoiBasesCorrect / (double)RoiBases;
        public double ROIDirectionAccuracy => 100 * RoiBasesCorrectDirection / (double)RoiBases;

        public MetricsCalculator(long totalBases, long totalBasesRight, long totalBasesRightDirection, long isGainBases, long callGainBases,
                long isGainBasesCorrect, long isGainBasesCorrectDirection, long isLossBases, long callLossBases, long isLossBasesCorrect,
                long isLossBasesCorrectDirection, long roiBases, long roiBasesCorrect, long roiBasesCorrectDirection)
        {
            TotalBases = totalBases;
            TotalBasesRight = totalBasesRight;
            TotalBasesRightDirection = totalBasesRightDirection;

            IsGainBases = isGainBases;
            CallGainBases = callGainBases;
            IsGainBasesCorrect = isGainBasesCorrect;
            IsGainBasesCorrectDirection = isGainBasesCorrectDirection;

            IsLossBases = isLossBases;
            CallLossBases = callLossBases;
            IsLossBasesCorrect = isLossBasesCorrect;
            IsLossBasesCorrectDirection = isLossBasesCorrectDirection;

            RoiBases = roiBases;
            RoiBasesCorrect = roiBasesCorrect;
            RoiBasesCorrectDirection = roiBasesCorrectDirection;
        }

        public static MetricsCalculator CalculateMetrics(BaseCounter baseCounter, int maxCn, int maxPloidy)
        {
            long totalBases = 0;
            long totalBasesRight = 0;
            long totalBasesRightDirection = 0;

            long isGainBases = 0;
            long callGainBases = 0;
            long isGainBasesCorrect = 0;
            long isGainBasesCorrectDirection = 0;

            long isLossBases = 0;
            long callLossBases = 0;
            long isLossBasesCorrect = 0;
            long isLossBasesCorrectDirection = 0;

            for (var ploidy = 0; ploidy <= maxPloidy; ploidy++)
            {
                for (var trueCn = 0; trueCn <= maxCn; trueCn++)
                {
                    for (var callCn = 0; callCn <= maxCn; callCn++)
                    {
                        long bases = baseCounter.BaseCount[trueCn, callCn, ploidy];
                        totalBases += bases;
                        if (trueCn == callCn) totalBasesRight += bases;
                        if (trueCn < ploidy && callCn < ploidy || trueCn == ploidy && callCn == ploidy || trueCn > ploidy && callCn > ploidy)
                            totalBasesRightDirection += bases;
                        if (trueCn < ploidy) isLossBases += bases;
                        if (trueCn > ploidy) isGainBases += bases;
                        if (callCn < ploidy) callLossBases += bases;
                        if (callCn > ploidy) callGainBases += bases;
                        if (trueCn == callCn && trueCn < ploidy) isLossBasesCorrect += bases;
                        if (trueCn == callCn && trueCn > ploidy) isGainBasesCorrect += bases;
                        if (trueCn > ploidy && callCn > ploidy) isGainBasesCorrectDirection += bases;
                        if (trueCn < ploidy && callCn < ploidy) isLossBasesCorrectDirection += bases;
                    }
                }
            }


            // Compute ROI stats:
            long roiBases = 0;
            long roiBasesCorrect = 0;
            long roiBasesCorrectDirection = 0;
            if (baseCounter.RoiBaseCount != null)
            {
                for (var ploidy = 0; ploidy <= maxPloidy; ploidy++)
                {
                    for (var trueCn = 0; trueCn <= maxCn; trueCn++)
                    {
                        for (var callCn = 0; callCn <= maxCn; callCn++)
                        {
                            long bases = baseCounter.RoiBaseCount[trueCn, callCn, ploidy];
                            roiBases += bases;
                            if (trueCn == callCn) roiBasesCorrect += bases;
                            if (trueCn < ploidy && callCn < ploidy || trueCn == ploidy && callCn == ploidy || trueCn > ploidy && callCn > ploidy)
                                roiBasesCorrectDirection += bases;
                        }
                    }
                }
            }
            return new MetricsCalculator(totalBases, totalBasesRight, totalBasesRightDirection, isGainBases, callGainBases,
            isGainBasesCorrect, isGainBasesCorrectDirection, isLossBases, callLossBases, isLossBasesCorrect,
            isLossBasesCorrectDirection, roiBases, roiBasesCorrect, roiBasesCorrectDirection);
        }
    }
}
