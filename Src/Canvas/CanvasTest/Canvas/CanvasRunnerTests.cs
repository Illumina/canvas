using System;
using System.IO;
using System.Runtime.CompilerServices;
using System.Text;
using Canvas;
using Canvas.Visualization;
using CanvasCommon;
using Illumina.Common.FileSystem;
using Isas.Framework.Checkpointing;
using Isas.Framework.Logging;
using Isas.Framework.WorkManagement;
using Isas.Framework.WorkManagement.CommandBuilding;
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
            var workDoer = Substitute.For<IWorkDoer>();
            var checkpointRunner = Substitute.For<ICheckpointRunner>();
            Func<string, ICommandFactory> runtimePrefix = component => Substitute.For<ICommandFactory>();
            string dotnetPath = @"C:\path\to\dotnet.exe";
            var runtimeExecutable = new FileLocation(dotnetPath);
            bool isSomatic = true;
            var coverageMode = new CanvasCoverageMode();
            int countsPerBin = 0;
            string canvasFolder = @"C:\path\to\Canvas\"; 
            var bAlleleBedGraphWriter = Substitute.For<IBAlleleBedGraphWriter>();
            var canvasRunner = new CanvasRunner(logger, workDoer, checkpointRunner, runtimeExecutable, runtimePrefix, isSomatic, coverageMode, countsPerBin, bAlleleBedGraphWriter, null, canvasFolder);
            string prefix = "something before ";
            var commandLineBuilder = new StringBuilder(prefix);
            string canvasExecutableStub = "CanvasBin";
            string fullName = canvasRunner.GetExecutablePath(canvasExecutableStub, commandLineBuilder);

            Assert.Equal(@"C:\path\to\dotnet.exe", fullName);
            Assert.Equal(@"something before C:\path\to\Canvas\CanvasBin\CanvasBin.dll ", commandLineBuilder.ToString());
        }
    }
}
