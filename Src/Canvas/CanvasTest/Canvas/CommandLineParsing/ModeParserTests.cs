using System;
using System.IO;
using Canvas.CommandLineParsing;
using CanvasCommon.CommandLineParsing.CoreOptionTypes;
using CanvasCommon.CommandLineParsing.OptionProcessing;
using Illumina.Common.FileSystem;
using Ploeh.AutoFixture.Xunit2;
using Xunit;
using System.Linq;

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

        private static MainParser GetMainParser(ModeParser germlineWgsModeParser)
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
        public void ParsePositionalOption_WithMultipleStringInputArguments_ReturnsMultipleStrings()
        {
            // arrange
            string inputArgument1 = "input";
            string inputArgument2 = "input2";
            string inputArgument3 = "input3";
            Func<string, string, string, ParsingResult<Tuple<string, string, string>>> parse = (value1, value2, value3) => ParsingResult<Tuple<string, string, string>>.SuccessfulResult(Tuple.Create(value1, value2, value3));
            var option1 = StringOption.CreateRequired("value1 option", "value1");
            var option2 = StringOption.CreateRequired("value2 option", "value2");
            var option3 = StringOption.CreateRequired("value3 option", "value2");
            var multiStringOption = new PositionalOption<string, string, string, Tuple<string, string, string>>(parse, true, option1, option2, option3, "positional-option");
            string[] stringInputArgument =
            {
                "--positional-option", inputArgument1, inputArgument2, inputArgument3
            };

            // act
            var result = multiStringOption.Parse(stringInputArgument);

            // assert
            Assert.True(result.Success);
            Assert.Equal(inputArgument1, result.Result[0].Item1);
            Assert.Equal(inputArgument2, result.Result[0].Item2);
            Assert.Equal(inputArgument3, result.Result[0].Item3);
        }

        [Fact]
        public void ParsePositionalOption_WithDefaultInputArguments_ReturnsMultipleStrings()
        {
            // arrange
            string inputArgument1 = "input";
            string inputArgument2 = "input2";
            Func<string, string, string, ParsingResult<Tuple<string, string, string>>> parse = (value1, value2, value3) => ParsingResult<Tuple<string, string, string>>.SuccessfulResult(Tuple.Create(value1, value2, value3));
            var option1 = StringOption.CreateRequired("value1 option", "value1");
            var option2 = StringOption.CreateRequired("value2 option", "value2");
            var option3 = StringOption.Create("value3 option", "value2");
            var multiStringOption = new PositionalOption<string, string, string, Tuple<string, string, string>>(parse, true, option1, option2, option3, "positional-option");
            string[] stringInputArgument =
            {
                "--positional-option", inputArgument1, inputArgument2
            };

            // act
            var result = multiStringOption.Parse(stringInputArgument);

            // assert
            Assert.True(result.Success);
            Assert.Equal(inputArgument1, result.Result[0].Item1);
            Assert.Equal(inputArgument2, result.Result[0].Item2);
            Assert.Equal(null, result.Result[0].Item3);
        }

        [Fact]
        public void ParsePositionalOption_WithMissingInputArguments_ReturnsFailedResult()
        {
            // arrange
            string inputArgument1 = "input";
            string inputArgument2 = "input2";
            Func<string, string, string, ParsingResult<Tuple<string, string, string>>> parse = (value1, value2, value3) => ParsingResult<Tuple<string, string, string>>.SuccessfulResult(Tuple.Create(value1, value2, value3));
            var option1 = StringOption.CreateRequired("value1 option", "value1");
            var option2 = StringOption.CreateRequired("value2 option", "value2");
            var option3 = StringOption.Create("value3 option", "value2");
            var multiStringOption = new PositionalOption<string, string, string, Tuple<string, string, string>>(parse, true, option1, option2, option3, "positional-option");
            string[] stringInputArgument =
            {
                "--positional-option", inputArgument1
            };

            // act
            var result = multiStringOption.Parse(stringInputArgument);

            // assert
            Assert.False(result.Success);
        }

        [Fact]
        public void ParseRequiredPositionalOption_WithNoInputArguments_ReturnsFailedResult()
        {
            Func<string, string, string, ParsingResult<Tuple<string, string, string>>> parse = (value1, value2, value3) => ParsingResult<Tuple<string, string, string>>.SuccessfulResult(Tuple.Create(value1, value2, value3));
            var option1 = StringOption.Create("value1 option", "value1");
            var option2 = StringOption.Create("value2 option", "value2");
            var option3 = StringOption.Create("value3 option", "value2");
            var multiStringOption = new PositionalOption<string, string, string, Tuple<string, string, string>>(parse, true, option1, option2, option3, "positional-option");
            string[] stringInputArgument =
            {
            };

            // act
            var result = multiStringOption.Parse(stringInputArgument);

            // assert
            Assert.False(result.Success);
        }

        [Fact]
        public void ParsePositionalOption_WithTooManyInputArguments_ReturnsFailedResult()
        {
            Func<string, string, string, ParsingResult<Tuple<string, string, string>>> parse = (value1, value2, value3) => ParsingResult<Tuple<string, string, string>>.SuccessfulResult(Tuple.Create(value1, value2, value3));
            var option1 = StringOption.Create("value1 option", "value1");
            var option2 = StringOption.Create("value2 option", "value2");
            var option3 = StringOption.Create("value3 option", "value2");
            var multiStringOption = new PositionalOption<string, string, string, Tuple<string, string, string>>(parse, true, option1, option2, option3, "positional-option");
            OptionCollection<Tuple<Tuple<string, string, string>, string>> collection = new OptionCollection<Tuple<Tuple<string, string, string>, string>>();
            collection.Add(multiStringOption);
            collection.Add(option1);
            string[] stringInputArgument =
            {
                "--positional-option", "1", "2", "3", "4", "--value1", "value"
            };

            // act
            var result = collection.Parse(stringInputArgument);

            // assert
            ParsingResult<Tuple<Tuple<string, string, string>, string>> failedResult;
            Assert.False(result.Validate(out failedResult));
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

        [Fact]
        public void ParseSampleType_ReturnsSuccessfulResult()
        {
            string[] stringInputArgument =
            {
                $"--{SmallPedigreeOptionsParser.SampleType.Info.Name}", "mother"
            };

            var result = SmallPedigreeOptionsParser.SampleType.Parse(stringInputArgument);

            Assert.True(result.Success);
            Assert.Equal(SampleType.Mother, result.Result);
        }

        [Fact]
        public void ParseBam_ReturnsSuccessfulResultWithDefaultValues()
        {
            string assemblyFolder = Isas.Framework.Utilities.Utilities.GetAssemblyFolder(typeof(ModeParserTests));
            var dataFolder = new DirectoryLocation(assemblyFolder).GetDirectoryLocation("Data");
            var bamPath = dataFolder.GetFileLocation("Tiny_COLO829BL_S1.bam");
            string[] stringInputArgument =
            {
                $"--{SmallPedigreeOptionsParser.Bams.Info.Name}", bamPath.FullName, $"--{SmallPedigreeOptionsParser.Bams.Info.Name}", bamPath.FullName, "father", "sampleID"
            };

            var result = SmallPedigreeOptionsParser.Bams.Parse(stringInputArgument);
            Assert.True(result.Success);

            var sampleResult = result.Result.First();
            Assert.Equal(bamPath, sampleResult.Bam);
            Assert.Equal(SampleType.Other, sampleResult.SampleType);
            Assert.Equal("COLO829BL", sampleResult.SampleName);

            sampleResult = result.Result[1];
            Assert.Equal(bamPath, sampleResult.Bam);
            Assert.Equal(SampleType.Father, sampleResult.SampleType);
            Assert.Equal("sampleID", sampleResult.SampleName);
        }

        [Fact]
        public void ParsePedigreeSample_ReturnsSuccessfulResultWithDefaultValues()
        {
            string assemblyFolder = Isas.Framework.Utilities.Utilities.GetAssemblyFolder(typeof(ModeParserTests));
            var dataFolder = new DirectoryLocation(assemblyFolder).GetDirectoryLocation("Data");
            var bamPath = dataFolder.GetFileLocation("Tiny_COLO829BL_S1.bam");
            string[] stringInputArgument =
            {
                $"--{SmallPedigreeOptionsParser.Bams.Info.Name}", bamPath.FullName, "--sample-b-allele-vcf", bamPath.FullName
            };

            var result = new SmallPedigreeOptionsParser().Parse(stringInputArgument);
            Assert.True(result.Success);
        }

        [Fact]
        public void Parse_SpwWgsMode_ReturnsSuccess()
        {
            // arrange
            var germlineWgsModeParser = new SmallPedigreeModeParser("SmallPedigree-WGS", "run spw");
            MainParser parser = GetMainParser(germlineWgsModeParser);
            string[] args =
            {
                "SmallPedigree-WGS",
                "--bam", "/illumina/scratch/STE_DataAnalysis/repo/SmallPedigreeWorkflow/Chr15Bam/NA12877-PcrFree_S1.REF_chr15.bam", "father", "NA12877-PcrFree",
                "--bam", "/illumina/scratch/STE_DataAnalysis/repo/SmallPedigreeWorkflow/Chr15Bam/NA12878-PcrFree_S2.REF_chr15.bam", "mother", "NA12878-PcrFree",
                "--bam", "/illumina/scratch/STE_DataAnalysis/repo/SmallPedigreeWorkflow/Chr15Bam/NA12882-PcrFree_S3.REF_chr15.bam", "proband", "NA12882-PcrFree",
                "--reference", "/illumina/development/Isas/Genomes/Homo_sapiens/NCBI/GRCh38Decoy/Annotation/Canvas/kmerv2.fa",
                "--genome-folder", "/illumina/development/Isas/Genomes/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta",
                "--filter-bed", "/illumina/development/Isas/Genomes/Homo_sapiens/NCBI/GRCh38Decoy/Annotation/Canvas/filter13.bed",
                "--output", "/illumina/scratch/bioinfoSD/eroller/SPW_NewCanvas/Temp_06-DetectCNV",
                "--sample-b-allele-vcf", "/illumina/scratch/bioinfoSD/eroller/SPW_NewCanvas/Pedigree.vcf.gz",
                "--ploidy-vcf", "/illumina/scratch/bioinfoSD/eroller/SPW_NewCanvas/Temp_06-DetectCNV/ploidy.vcf.gz",
                "--custom-parameters", "CanvasPartition,--commoncnvs /illumina/development/Isas/Genomes/Homo_sapiens/NCBI/GRCh38Decoy/Annotation/Canvas/commoncnvs.bed"
            };

            // act
            using (var standardWriter = new StringWriter())
            using (var errorWriter = new StringWriter())
            {
                var result = parser.Parse(args, standardWriter, errorWriter);
                // assert
                Assert.False(result.Success);
                Assert.Contains("does not exist", result.ErrorMessage);
            }
        }

    }
}
