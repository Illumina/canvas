using Microsoft.VisualStudio.TestTools.UnitTesting;
using CanvasCommon;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CanvasTest
{
    [TestClass()]
    public class TestUtilities
    {
        [TestMethod()]
        public void TestGoldenSectionSearch()
        {
            double res = Utilities.GoldenSectionSearch(x => x * x, -5, 5);
            Assert.IsTrue(Math.Abs(res) < 0.001);

            res = Utilities.GoldenSectionSearch(x => x * x, 0, 5);
            Assert.IsTrue(Math.Abs(res) < 0.001);

            res = Utilities.GoldenSectionSearch(x => x * x, -5, 0);
            Assert.IsTrue(Math.Abs(res) < 0.001);
        }
    }
}