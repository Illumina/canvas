using Canvas.CommandLineParsing;
using CanvasCommon.CommandLineParsing.OptionProcessing;
using Isas.Framework.UnitTestUtilities;
using Xunit;

namespace CanvasTest
{
    public class SerializationTests
    {
        private readonly TemporaryDirectoryFixture _tempDirectory;
        private readonly SerializationTester _tester;

        public SerializationTests()
        {
            _tempDirectory = new TemporaryDirectoryFixture();
            _tester = new SerializationTester(_tempDirectory);
        }

        ~SerializationTests()
        {
            _tempDirectory.Dispose();
        }

        [Fact]
        public void CanSerializeFailedParsingResult()
        {
            var parsingResult = ParsingResult<SmallPedigreeInput>.FailedResult("error");

            var serializedResult = _tester.RoundTrip(parsingResult);

            Assert.Equal("error", serializedResult.ErrorMessage);
        }

        [Fact]
        public void CanSerializeSuccessfulParsingResult()
        {
            var parsingResult = ParsingResult<string>.SuccessfulResult("success");

            var serializedResult = _tester.RoundTrip(parsingResult);

            Assert.Equal("success", serializedResult.Result);
        }
    }
}