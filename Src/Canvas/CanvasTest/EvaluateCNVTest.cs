
using System;
using System.Collections.Generic;
using CanvasCommon;
using EvaluateCNV;
using Xunit;
using CNInterval = EvaluateCNV.CNInterval;

namespace CanvasTest
{
    public class EvaluateCNVTest
    {
        [Fact]
        public void TestAllosomes()
        {
            var cnvEvaluator = new CnvEvaluator(new CNVChecker(null));

            var baseCounter = new BaseCounter(5, 0, 4999);
            const string chr = "1";
            var calls = new Dictionary<string, List<CnvCall>>
            {
                [chr] = new List<CnvCall>
                {
                    new CnvCall(chr, start: 1, end: 1000, cn: 2, refPloidy: 1, passFilter: true, altAllele: "."),
                    new CnvCall(chr, start: 2001, end: 3000, cn: 1, refPloidy: 2, passFilter: true, altAllele: "."),
                    new CnvCall(chr, start: 4001, end: 5000, cn: 2, refPloidy: 1, passFilter: true, altAllele: "."),
                    new CnvCall(chr, start: 6001, end: 7000, cn: 2, refPloidy: 2, passFilter: true, altAllele: ".")
                }
            };

            var knownCN = new Dictionary<string, List<CNInterval>>
            {
                [chr] = new List<CNInterval>
                {
                    new CNInterval(chr, start: 1, end: 1000, cn: 2, referenceCopyNumber: 1),
                    new CNInterval(chr, start: 2001, end: 3000, cn: 1, referenceCopyNumber: 2),
                    new CNInterval(chr, start: 3001, end: 4000, cn: 1, referenceCopyNumber: 2),
                    new CNInterval(chr, start: 4001, end: 5000, cn: 1, referenceCopyNumber: 1),
                    new CNInterval(chr, start: 6001, end: 7000, cn: 2, referenceCopyNumber: 2)
                }
            };
            var metrics = cnvEvaluator.CalculateMetrics(knownCN, calls, baseCounter, optionsSkipDiploid: false, includePassingOnly: true);
            Assert.Equal(100, Convert.ToInt32(metrics.Recall));
        }

        [Fact]
        public void TestFalseNegatives()
        {
            var cnvEvaluator = new CnvEvaluator(new CNVChecker(null));

            var baseCounter = new BaseCounter(5, 0, 4999);
            const string chr = "1";
            var calls = new Dictionary<string, List<CnvCall>>
            {
                [chr] = new List<CnvCall>
                {
                    new CnvCall(chr, start: 1, end: 1000, cn: 2, refPloidy: 1, passFilter: false, altAllele: "."),
                    new CnvCall(chr, start: 2001, end: 3000, cn: 1, refPloidy: 2, passFilter: true, altAllele: "."),
                    new CnvCall(chr, start: 4001, end: 5000, cn: 2, refPloidy: 1, passFilter: true, altAllele: "."),
                    new CnvCall(chr, start: 6001, end: 7000, cn: 2, refPloidy: 2, passFilter: true, altAllele: ".")
                }
            };

            var knownCN = new Dictionary<string, List<CNInterval>>
            {
                [chr] = new List<CNInterval>
                {
                    new CNInterval(chr, start: 1, end: 1000, cn: 2, referenceCopyNumber: 1),
                    new CNInterval(chr, start: 2001, end: 3000, cn: 1, referenceCopyNumber: 2),
                    new CNInterval(chr, start: 3001, end: 4000, cn: 1, referenceCopyNumber: 2),
                    new CNInterval(chr, start: 4001, end: 5000, cn: 1, referenceCopyNumber: 1),
                    new CNInterval(chr, start: 6001, end: 7000, cn: 2, referenceCopyNumber: 2)
                }
            };
            var metrics = cnvEvaluator.CalculateMetrics(knownCN, calls, baseCounter, optionsSkipDiploid: false, includePassingOnly: true);
            Assert.Equal(33, Convert.ToInt32(metrics.Recall));
        }
    }
}
