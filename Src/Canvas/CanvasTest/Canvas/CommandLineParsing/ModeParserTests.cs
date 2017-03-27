using System.IO;
using Canvas.CommandLineParsing;
using CanvasCommon.CommandLineParsing.CoreOptionTypes;
using CanvasCommon.CommandLineParsing.OptionProcessing;
using Illumina.Common.FileSystem;
using Ploeh.AutoFixture.Xunit2;
using Xunit;

namespace CanvasTest.Canvas.CommandLineParsing
{
    public class ModeParserTests
    {
        private const string Version = "Version";
        private const string Copyright = "Copyright";

        [Theory]
        [InlineAutoNSubstituteData("Error: no mode specified")]
        [InlineAutoNSubstituteData("Available modes:")]
        [InlineAutoNSubstituteData("ModeName - ModeDescription", "ModeName", "ModeDescription")]
        [InlineAutoNSubstituteData("Options:")]
        [InlineAutoNSubstituteData("-h, --help")]
        [InlineAutoNSubstituteData("-v, --version")]
        public void Parse_NoArguments_DisplaysError(string messageToDisplay, string name, string description, StringWriter standardWriter,
            StringWriter errorWriter)
        {
            // arrange
            GermlineWgsModeParser germlineWgsModeParser = new GermlineWgsModeParser(name, description);
            MainParser parser = GetMainParser(germlineWgsModeParser);
            string[] noArgs =
            {
            };

            // act
            var result = parser.Parse(noArgs, standardWriter, errorWriter);
            string errorOutput = errorWriter.ToString();

            // assert
            Assert.False(result.Success);
            Assert.Contains(messageToDisplay, errorOutput);
            Assert.Empty(standardWriter.ToString());
        }

        private static MainParser GetMainParser(GermlineWgsModeParser germlineWgsModeParser)
        {

            return new MainParser(Version, Copyright, germlineWgsModeParser);
        }

        [Theory]
        [InlineAutoNSubstituteData("Available modes:")]
        [InlineAutoNSubstituteData("ModeName - ModeDescription", "ModeName", "ModeDescription")]
        [InlineAutoNSubstituteData("Options:")]
        [InlineAutoNSubstituteData("-h, --help")]
        [InlineAutoNSubstituteData("-v, --version")]
        public void Parse_WithHelpArgument_ReturnsSuccecssAndDisplaysHelp(string messageToDisplay, string name, string description, StringWriter standardWriter,
            StringWriter errorWriter)
        {
            // arrange
            GermlineWgsModeParser germlineWgsModeParser = new GermlineWgsModeParser(name, description);
            MainParser parser = GetMainParser(germlineWgsModeParser);
            string[] args =
            {
                "-h"
            };

            // act
            var result = parser.Parse(args, standardWriter, errorWriter);
            string standardOutput = standardWriter.ToString();

            // assert
            Assert.True(result.Success);
            Assert.Contains(messageToDisplay, standardOutput);
            Assert.Empty(errorWriter.ToString());
        }

        [Theory]
        [InlineAutoNSubstituteData("Error: found unexpected arguments '--unknown-option'")]
        [InlineAutoNSubstituteData("Available modes:")]
        [InlineAutoNSubstituteData("ModeName - ModeDescription", "ModeName", "ModeDescription")]
        [InlineAutoNSubstituteData("Options:")]
        [InlineAutoNSubstituteData("-h, --help")]
        [InlineAutoNSubstituteData("-v, --version")]
        public void Parse_WithHelpArgumentAndUnkownArgument_DisplaysError(string messageToDisplay, string name, string description, StringWriter standardWriter,
            StringWriter errorWriter)
        {
            // arrange
            GermlineWgsModeParser germlineWgsModeParser = new GermlineWgsModeParser(name, description);
            MainParser parser = GetMainParser(germlineWgsModeParser);
            string[] args =
            {
                "-h", "--unknown-option"
            };

            // act
            var result = parser.Parse(args, standardWriter, errorWriter);
            string errorOutput = errorWriter.ToString();

            // assert
            Assert.False(result.Success);
            Assert.Contains(messageToDisplay, errorOutput);
            Assert.Empty(standardWriter.ToString());
        }

        [Theory]
        [InlineAutoNSubstituteData("required")]
        public void Parse_ModeWithMissingRequiredArgument_DisplaysError(string messageToDisplay, string name, string description, StringWriter standardWriter,
            StringWriter errorWriter)
        {
            // arrange
            GermlineWgsModeParser germlineWgsModeParser = new GermlineWgsModeParser("WGS", "Run Canvas from WGS data");
            MainParser parser = GetMainParser(germlineWgsModeParser);
            string[] modeArgs =
            {
                "WGS"
            };

            // act
            parser.Parse(modeArgs, standardWriter, errorWriter);
            string errorOutput = errorWriter.ToString();

            // assert
            Assert.Contains(messageToDisplay, errorOutput);
            Assert.Empty(standardWriter.ToString());
        }

        [Theory]
        [InlineAutoNSubstituteData("required")]
        public void Parse_ModeWithVersion_ReturnsSuccecssAndDisplaysVersion(string messageToDisplay, string name, string description,
            StringWriter standardWriter, StringWriter errorWriter)
        {
            // arrange
            GermlineWgsModeParser germlineWgsModeParser = new GermlineWgsModeParser("WGS", "Run Canvas from WGS data");
            MainParser parser = GetMainParser(germlineWgsModeParser);
            string[] modeArgs =
            {
                "WGS", "-v"
            };

            // act
            var result = parser.Parse(modeArgs, standardWriter, errorWriter);
            string output = standardWriter.ToString().Trim();

            // assert
            Assert.True(result.Success);
            Assert.Equal(Version, output);
            Assert.Empty(errorWriter.ToString());
        }

        [Fact]
        public void ParseOptionInfo_WithStringInputArgument_ReturnsStringInputArgument()
        {
            // arrange
            string inputArgument = "input";
            StringOption valueOption = StringOption.Create("value option", "value");
            string[] stringInputArgument =
            {
                "--value", inputArgument
            };
            // act
            string parsedResult = valueOption.Parse(stringInputArgument).Result;

            // assert
            Assert.Equal(inputArgument, parsedResult);
        }

        [Fact]
        public void ParseMultiStringOption_WithMultipleStringInputArguments_ReturnsMultipleStrings()
        {
            // arrange
            string inputArgument1 = "input";
            string inputArgument2 = "input2";
            var multiStringOption = new MultiValueOption<string>(StringOption.Create("value option", "value"));
            string[] stringInputArgument =
            {
                "--value", inputArgument1, "--value", inputArgument2
            };

            // act
            var result = multiStringOption.Parse(stringInputArgument);

            // assert
            Assert.True(result.Success);
            Assert.Equal(inputArgument1, result.Result[0]);
            Assert.Equal(inputArgument2, result.Result[1]);
        }

        [Fact]
        public void Description_RequiredMultiStringOption_ShouldIncludeMultipleOccurrencesAllowedAndRequired()
        {
            // arrange
            var multiStringOption = new MultiValueOption<string>(StringOption.CreateRequired("value option", "value"));

            // act
            var description = multiStringOption.Info.Description;

            // assert
            Assert.Equal("value option. Option can be specified multiple times. (required)", description);
        }

        [Fact]
        public void ParseOptionInfo_WithMissingInputArgument_ReturnsFailedResult()
        {
            // arrange
            StringOption valueOption = StringOption.Create("value option", "value");
            string[] stringInputArgument =
            {
                "--value"
            };

            // act
            ParsingResult<string> result = valueOption.Parse(stringInputArgument);

            // assert
            Assert.False(result.Success);
            Assert.Contains("required", result.ErrorMessage);
        }

        [Theory]
        [AutoData]
        public void ParseCommonOptions_WithRequiredArguments_ReturnsSuccessfulResult(TemporaryDirectoryFixture tempDirectory)
        {
            // arrange
            Option<CommonOptions> commonOptionsParser = new CommonOptionsParser();
            var kmerFasta = tempDirectory.CreateFile("kmerv2.fa");
            var filterBed = tempDirectory.GetFileLocation("filter.bed").Touch();
            var output = tempDirectory.CreateSubdirectory("output");
            var genome = tempDirectory.CreateSubdirectory("WholeGenomeFasta");
            string[] stringInputArgument =
            {
                "-r", kmerFasta.ToString(), "-o", output.ToString(), "-g", genome.ToString(), "--filter-bed", filterBed.ToString()
            };

            // act
            ParsingResult<CommonOptions> result = commonOptionsParser.Parse(stringInputArgument);

            // assert
            Assert.Equal("", result.ErrorMessage);
            Assert.True(result.Success);
            Assert.Equal(kmerFasta, result.Result.KmerFasta);
            Assert.Equal(output, result.Result.OutputDirectory);
            Assert.Equal(genome, result.Result.WholeGenomeFasta);
        }

        [Theory]
        [AutoData]
        public void ParseCommonOptions_KmerFastaDoesntExist_DisplaysError(TemporaryDirectoryFixture tempDirectory)
        {
            // arrange
            Option<CommonOptions> commonOptionsParser = new CommonOptionsParser();
            var kmerFasta = tempDirectory.GetFileLocation("kmerv2.fa");
            var filterBed = tempDirectory.GetFileLocation("filter.bed").Touch();
            var output = tempDirectory.CreateSubdirectory("output");
            var genome = tempDirectory.CreateSubdirectory("WholeGenomeFasta");
            string[] stringInputArgument =
            {
                "-r", kmerFasta.ToString(), "-o", output.ToString(), "-g", genome.ToString(), "--filter-bed", filterBed.ToString()
            };

            // act
            ParsingResult<CommonOptions> result = commonOptionsParser.Parse(stringInputArgument);

            // assert
            Assert.False(result.Success);
            Assert.Contains("kmerv2.fa", result.ErrorMessage);
            Assert.Contains("does not exist", result.ErrorMessage);
        }

        [Theory]
        [AutoData]
        public void GermlineWgsParse_WithRequiredArguments_ReturnsSuccessfulCallsetResult(string name, string description, StringWriter standardWriter,
            StringWriter errorWriter, TemporaryDirectoryFixture tempDirectory)
        {
            // arrange
            GermlineWgsModeParser germlineWgsModeParser = new GermlineWgsModeParser(name, description);
            MainParser parser = GetMainParser(germlineWgsModeParser);
            string[] args =
            {
                "-h"
            };

            // act
            var result = parser.Parse(args, standardWriter, errorWriter);

            // assert
            Assert.True(result.Success);
        }

        [Theory]
        [AutoData]
        public void ParseMultiFileOption_WithMultipleFileArguments_ReturnsListOfFileLocations(TemporaryDirectoryFixture tempDirectory)
        {
            // arrange
            MultiValueOption<IFileLocation> multiFileOption = new MultiValueOption<IFileLocation>(FileOption.CreateRequired("multiple files", "file"));
            var file1 = tempDirectory.CreateFile("file1");
            var file2 = tempDirectory.CreateFile("file2");
            string[] args =
            {
                "--file", file1.ToString(), "--file", file2.ToString()
            };

            // act
            var result = multiFileOption.Parse(args);

            // assert
            Assert.Equal("", result.ErrorMessage);
            Assert.True(result.Success);
            Assert.Equal(file1, result.Result[0]);
            Assert.Equal(file2, result.Result[1]);
        }

        [Theory]
        [AutoData]
        public void ParseExclusiveOption_WithOnlyOneOption_ReturnsOneValue(TemporaryDirectoryFixture tempDirectory)
        {
            FileOption option1 = FileOption.CreateRequired("file1 option", "file1");
            FileOption option2 = FileOption.CreateRequired("file2 option", "file2");
            // arrange
            ExclusiveFileOption multiFileOption = ExclusiveFileOption.CreateRequired(option1, option2);
            var file1 = tempDirectory.CreateFile("file1");
            string[] args =
            {
                "--file1", file1.ToString()
            };

            // act
            var result = multiFileOption.Parse(args);

            // assert
            Assert.Equal("", result.ErrorMessage);
            Assert.True(result.Success);
            Assert.Equal(file1, result.Result.Result);
            Assert.Equal(option1, result.Result.MatchedOption);
        }

        [Theory]
        [AutoData]
        public void ParseExclusiveOption_WithOnlyTwoOption_ReturnsFailedParseResult(TemporaryDirectoryFixture tempDirectory)
        {
            FileOption option1 = FileOption.CreateRequired("file1 option", "file1");
            FileOption option2 = FileOption.CreateRequired("file2 option", "file2");
            // arrange
            ExclusiveFileOption multiFileOption = ExclusiveFileOption.CreateRequired(option1, option2);
            var file1 = tempDirectory.CreateFile("file1");
            string[] args =
            {
                "--file1", file1.ToString(), "--file2", file1.ToString()
            };

            // act
            var result = multiFileOption.Parse(args);

            // assert
            Assert.Contains("not both", result.ErrorMessage);
            Assert.False(result.Success);
        }

        [Theory]
        [AutoData]
        public void ParseRequiredExclusiveOption_WithNeitherOptionSpecified_ReturnsFailedParseResult(TemporaryDirectoryFixture tempDirectory)
        {
            FileOption option1 = FileOption.CreateRequired("file1 option", "file1");
            FileOption option2 = FileOption.CreateRequired("file2 option", "file2");
            // arrange
            ExclusiveFileOption multiFileOption = ExclusiveFileOption.CreateRequired(option1, option2);
            string[] args =
            {

            };

            // act
            var result = multiFileOption.Parse(args);

            // assert
            Assert.Contains("must be specified", result.ErrorMessage);
            Assert.False(result.Success);
        }

        [Theory]
        [AutoData]
        public void ParseDictionaryOption_WithMultipleKeyValueArguments_ReturnsDictionary(TemporaryDirectoryFixture tempDirectory)
        {
            // arrange
            DictionaryOption dictOption = DictionaryOption.Create("dictionary", "kvp");
            string key1 = "key1";
            string value1 = "value1";
            string key2 = "key2";
            string value2 = "value2";
            string[] args =
            {
                "--kvp", $"{key1}, {value1}","--kvp",  $"{key2}, {value2}"
            };

            // act
            var result = dictOption.Parse(args);

            // assert
            Assert.Equal("", result.ErrorMessage);
            Assert.True(result.Success);
            Assert.Equal(value1, result.Result[key1]);
            Assert.Equal(value2, result.Result[key2]);
        }

        [Theory]
        [AutoData]
        public void ParseDictionaryOption_WithKeyOnlyArgument_ReturnsFailedResult(TemporaryDirectoryFixture tempDirectory)
        {
            // arrange
            DictionaryOption dictOption = DictionaryOption.Create("dictionary", "kvp");
            string key = "key1";
            string[] args =
            {
                "--kvp", key
            };

            // act
            var result = dictOption.Parse(args);

            // assert
            Assert.False(result.Success);
            Assert.Contains("Error", result.ErrorMessage);
            Assert.Contains("format", result.ErrorMessage);
        }

        [Fact]
        public void ParseValueOption_ReturnsSuccessfulResult()
        {
            // arrange
            uint inputArgument = 100;
            ValueOption<uint> valueOption = ValueOption<uint>.Create("value option", "value");
            string[] stringInputArgument =
            {
                "--value", inputArgument.ToString()
            };
            // act
            uint parsedResult = valueOption.Parse(stringInputArgument).Result;

            // assert
            Assert.Equal(inputArgument, parsedResult);
        }

        [Fact]
        public void ParseNullableValueOption_ReturnsSuccessfulResult()
        {
            // arrange
            ValueOption<uint?> valueOption = ValueOption<uint?>.Create("value option", "value");

            // act
            uint? parsedResult = valueOption.Parse("").Result;

            // assert
            Assert.False(parsedResult.HasValue);
        }
    }
}
