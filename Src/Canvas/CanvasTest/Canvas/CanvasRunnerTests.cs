using System.IO;
using System.Runtime.CompilerServices;
using System.Text;
using Canvas;
using CanvasCommon;
using Illumina.Common.FileSystem;
using Isas.Framework.Checkpointing;
using Isas.Framework.Logging;
using Isas.Framework.WorkManagement;
using Xunit;
using NSubstitute;

namespace CanvasTest.Canvas
{
    public sealed class CanvasRunnerTests
    {
        [Fact]
        public void GetExecutablePath_test()
        {
            var logger = Substitute.For<ILogger>();
            var workManager = Substitute.For<IWorkManager>();
            var checkpointRunner = Substitute.For<ICheckpointRunner>();
            string dotnetPath = Path.Combine("path", "to", "executable");
            var runtimeExecutable = new FileLocation(dotnetPath);
            bool isSomatic = true;
            var coverageMode = new CanvasCoverageMode();
            int countsPerBin = 0;
            string canvasFolder = Path.Combine(Isas.Framework.Utilities.Utilities.GetAssemblyFolder(typeof(CanvasRunner)));
            var canvasRunner = new CanvasRunner(logger, workManager, checkpointRunner, runtimeExecutable, isSomatic, coverageMode, countsPerBin, null, canvasFolder);
            string prefix = "something before ";
            var commandLineBuilder = new StringBuilder(prefix);
            string canvasExecutableStub = "CanvasBin";
            string dllName = canvasExecutableStub + ".dll";
            string fullName = canvasRunner.GetExecutablePath(canvasExecutableStub, commandLineBuilder);

            Assert.Equal(Path.Combine(canvasFolder, dotnetPath), fullName);
            Assert.Equal(prefix + Path.Combine(canvasFolder, canvasExecutableStub, dllName) + " ", commandLineBuilder.ToString());
        }
    }
}
