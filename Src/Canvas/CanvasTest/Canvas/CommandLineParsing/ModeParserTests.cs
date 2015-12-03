using System.Collections.Generic;
using System.IO;
using System.Linq;
using Canvas;
using Canvas.CommandLineParsing;
using Canvas.CommandLineParsing.CoreOptionTypes;
using Isas.Shared;
using Ploeh.AutoFixture.Xunit2;
using UnitTests;
using Xunit;

namespace CanvasTest.Canvas.CommandLineParsing
{
    public class ModeParserTests
    {

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
            MainParser parser = new MainParser(germlineWgsModeParser);
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
            MainParser parser = new MainParser(germlineWgsModeParser);
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
            MainParser parser = new MainParser(germlineWgsModeParser);
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
            MainParser parser = new MainParser(germlineWgsModeParser);
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
            MainParser parser = new MainParser(germlineWgsModeParser);
            string[] modeArgs =
            {
                "WGS", "-v"
            };

            // act
            var result = parser.Parse(modeArgs, standardWriter, errorWriter);
            string output = standardWriter.ToString();

            // assert
            Assert.True(result.Success);
            string version = typeof(Program).Assembly.GetName().Version.ToString();
            Assert.Contains(version, output);
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
            var kmerFasta = tempDirectory.CreateFile("kmer.fa");
            var bAlleleVcf = tempDirectory.GetFileLocation("ballele.vcf").Touch();
            var filterBed = tempDirectory.GetFileLocation("filter.bed").Touch();
            var output = tempDirectory.CreateSubdirectory("output");
            var genome = tempDirectory.CreateSubdirectory("WholeGenomeFasta");
            string[] stringInputArgument =
            {
                "-r", kmerFasta.ToString(), "-o", output.ToString(), "-g", genome.ToString(), "--b-allele-vcf", bAlleleVcf.ToString(), "--filter-bed", filterBed.ToString(), "--sample-name", "SampleName"
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
            var kmerFasta = tempDirectory.GetFileLocation("kmer.fa");
            var bAlleleVcf = tempDirectory.GetFileLocation("ballele.vcf").Touch();
            var filterBed = tempDirectory.GetFileLocation("filter.bed").Touch();
            var output = tempDirectory.CreateSubdirectory("output");
            var genome = tempDirectory.CreateSubdirectory("WholeGenomeFasta");
            string[] stringInputArgument =
            {
                "-r", kmerFasta.ToString(), "-o", output.ToString(), "-g", genome.ToString(), "--b-allele-vcf", bAlleleVcf.ToString(), "--filter-bed", filterBed.ToString(), "--sample-name", "SampleName"
            };

            // act
            ParsingResult<CommonOptions> result = commonOptionsParser.Parse(stringInputArgument);

            // assert
            Assert.False(result.Success);
            Assert.Contains("kmer.fa", result.ErrorMessage);
            Assert.Contains("does not exist", result.ErrorMessage);
        }

        [Theory]
        [AutoData]
        public void GermlineWgsParse_WithRequiredArguments_ReturnsSuccessfulCallsetResult(string name, string description, StringWriter standardWriter,
            StringWriter errorWriter, TemporaryDirectoryFixture tempDirectory)
        {
            // arrange
            GermlineWgsModeParser germlineWgsModeParser = new GermlineWgsModeParser(name, description);
            MainParser parser = new MainParser(germlineWgsModeParser);
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
    }
}
