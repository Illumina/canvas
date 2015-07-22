using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;

namespace SequencingFiles
{
    public class CFReport
    {
        public enum CFReportType
        {
            CfDiagnostic,
            CfCarrierScreening
        }

        //These emun-strings get turned directly to string for the report.  We cant change them without changing the report.
        public enum CfVariantType
        {
            PolyTGPolyT,    //special cftr poly TG
            Del,            //special cftr deleltions
            Div,            //a normal smallindel
            SNV             //snp
        }

        public enum CfZygosity
        {
            //NA, //TJD: not used?
            Het,
            Hom,
            Reference  //TJD: formerly "Nd"
        };

        //nd is not detected

        public enum ControlType
        {
            Positive,
            Negative,
            Regular
        };

        public const string SampleDiagnosticColumnHeaders = "Variant Type\t" +
                                                            "Coordinate\t" +
                                                            "Chromosome\t" +
                                                            "Frequency\t" +
                                                            "Depth\t" +
                                                            "cDNA Name (HGVS)\t" +
                                                            "Protein Name (HGVS)\t" +
                                                            "dbSNP ID\t" +
                                                            "Reference\t" +
                                                            "Result\t" +
                                                            "Interpretation\t";

        public const string SampleCarrierColumnHeaders = "Mutations (Common Name)\t" +
                                                         "Mutation Type\t" +
                                                         "dbSNP rsID\t" +
                                                         "CFTR gene region\t" +
                                                         "Genomic location\t" +
                                                         "cDNA Name (HGVS)\t" +
                                                         "Protein Name (HGVS)\t" +
                                                         "Result\t" +
                                                         "Interpretation";

        // Positive control and regular sample: Display PASS if call rate is greater than or equal to 99% and FAIL otherwise
        public const float PositiveCallRateThreshold = 99f;

        public const float DiagnosticVariantCallRateThreshold = 99f;
        public const float CarrierScreenVariantCallRateThreshold = 50f;

        // Negative control sample: Display PASS if call rate is less than or equal to 2% and FAIL otherwise
        public const float NegativeCallRateThreshold = 10f;

        public const string CoordinatesNotCalledToken = "Coordinates Not Called:";
        public const string SampleIdToken = "Sample ID";
        public const string SampleNameToken = "Sample Name";
        public const string SampleCallRateToken = "Sample Call Rate";
        public const string SampleControlToken = "Control";
        public const string SampleCommentToken = "Comment";
        public const string PerformanceToken = "Performance";
        public const string RunIdToken = "Run ID";
        public const string RunDateToken = "Run Date";
        public const string SampleFailedToken = "Sample Failed";
        public const string AnalysisVerToken = "Analysis Version";

        public CFReport(CFReportType reportType, string reportPath)
            : this(reportType, reportPath, null, null, null)
        {
        }

        public CFReport(CFReportType reportType, string reportPath, string runId, string runDate, string analysisVer)
        {
            CfSamples = new Dictionary<string, CFSample>();
            ReportType = reportType;
            ReportPath = reportPath;
            ReportRunId = runId;
            ReportRunDate = runDate;
            AnalysisVer = analysisVer;
            if (ReportType == CFReportType.CfDiagnostic)
                VariantsDiagnosticResultList = new List<VariantCFDx>();
            else
                VariantsScreenResultList = new List<VariantScreenResults>();
        }

        public Dictionary<string, CFSample> CfSamples { get; set; }
        public List<VariantScreenResults> VariantsScreenResultList { get; set; }
        public List<VariantCFDx> VariantsDiagnosticResultList { get; set; }
        public CFReportType ReportType { get; set; }
        public string ReportRunId { get; set; }
        public string ReportRunDate { get; set; }
        public string ReportPath { get; set; }
        public string AnalysisVer { get; set; }

        public static CFReport LoadDiagnosticReport(string reportPath)
        {
            CFReport diagnosticReport = new CFReport(CFReportType.CfDiagnostic, reportPath);
            diagnosticReport.GetCfReportResults();
            return diagnosticReport;
        }

        public static CFReport LoadCarrierReport(string reportPath)
        {
            CFReport carrierReport = new CFReport(CFReportType.CfCarrierScreening, reportPath);
            carrierReport.GetCfReportResults();
            return carrierReport;
        }

        //public static string ConvertCfVariantTypeToString(CfVariantType variantType)
        //{
        //    switch (variantType)
        //    {
        //        case CfVariantType.PolyTGPolyT:
        //            return "PolyTGPolyT";
        //        case CfVariantType.Del:
        //            return "DEL";
        //        case CfVariantType.Div:
        //            return "DIV";
        //        case CfVariantType.SNV:
        //            return "SNV";
        //        default:
        //            return "SNV";
        //    }
        //}

        public static CfVariantType ConvertStringToCfVariantType(string variantTypeString)
        {
            switch (variantTypeString.ToLower())
            {
                case "polytg/polyt":
                case "polytgpolyt":
                    return CfVariantType.PolyTGPolyT;
                case "del":
                    return CfVariantType.Del;
                case "div":
                    return CfVariantType.Div;
                case "snv":
                    return CfVariantType.SNV;
                default:
                    return CfVariantType.SNV;
            }
        }

        public void SaveReportToFile()
        {
            string reportTitle = ReportType == CFReportType.CfDiagnostic
                                     ? "CF Clinical Sequencing Assay"
                                     : "CF 139-Variant Assay";

            string polyTGReference = "TGTGTGTGTGTGTGTGTGTGTGTTTTTTT";
            string shortennedPolyTG = "(TG)11(T)7";

            StringBuilder strBuilderReport = new StringBuilder();

            strBuilderReport.AppendLine(String.Format("Test\t{0}", reportTitle));
            strBuilderReport.AppendLine(String.Format("For In Vitro Diagnostic Use."));
            strBuilderReport.AppendLine(String.Format("{0}\t{1}", RunIdToken, ReportRunId));
            strBuilderReport.AppendLine(String.Format("{0}\t{1}", RunDateToken, ReportRunDate));
            strBuilderReport.AppendLine(String.Format("{0}\t{1}", AnalysisVerToken, AnalysisVer));

            //foreach (string key in CfSamples.Keys.ToList().OrderBy(sampleId=>sampleId))
            foreach (string key in CfSamples.Keys.ToList())
            {
                CFSample sample = CfSamples[key];
                strBuilderReport.AppendLine();
                strBuilderReport.AppendLine();
                strBuilderReport.AppendLine(String.Format("{0}\t{1}", SampleIdToken, key));
                strBuilderReport.AppendLine(String.Format("{0}\t{1}", SampleNameToken, sample.SampleName));
                strBuilderReport.AppendLine(string.Format("{0}\t{1}", SampleControlToken, sample.Control));
                strBuilderReport.AppendLine(string.Format("{0}\t{1}", SampleCommentToken, sample.Comment));

                // performance
                // Negative control sample: Display PASS if call rate is less than or equal to the threshold and FAIL otherwise
                // Positive control and regular sample: Display PASS if call rate is greater than or equal to 99% and FAIL otherwise
                bool isPass = false;
                if (sample.Control != null &&
                    sample.Control.Equals(ControlType.Negative.ToString(), StringComparison.InvariantCultureIgnoreCase))
                    isPass = sample.CallRate <= NegativeCallRateThreshold;
                else
                    isPass = sample.CallRate >= PositiveCallRateThreshold;
                strBuilderReport.AppendLine(string.Format("{0}\t{1}", PerformanceToken, isPass ? "Pass" : "Fail"));

                strBuilderReport.AppendLine(String.Format("{0}\t{1:0.00}%", SampleCallRateToken, sample.CallRate));
                strBuilderReport.AppendLine();

                if (ReportType == CFReportType.CfDiagnostic)
                {
                    strBuilderReport.AppendLine(SampleDiagnosticColumnHeaders);

                    if (sample.CallRate >= DiagnosticVariantCallRateThreshold)
                    {
                        foreach (VariantCFDx variantCFDx in VariantsDiagnosticResultList)
                        {
                            if (!variantCFDx.SampleID.Equals(key))
                                continue;

                            if (variantCFDx.VariantStatsForReport.Type == CfVariantType.PolyTGPolyT)
                            {
                                //in this case, we dont report the frequency or depth
                                variantCFDx.VariantStatsForReport.Depth = -1;
                                variantCFDx.VariantStatsForReport.Frequency = -1;
                            }

                            string frequency = variantCFDx.VariantStatsForReport.Frequency <= 0
                                                   ? "N/A"
                                                   : variantCFDx.VariantStatsForReport.Frequency.ToString("0.0000000");
                            string depth = String.Format("{0}",
                                                         variantCFDx.VariantStatsForReport.Depth <= 0
                                                             ? "N/A"
                                                             : variantCFDx.VariantStatsForReport.Depth.ToString(
                                                                 CultureInfo.InvariantCulture));

                            //variantType, coordinate, chromosome, frequency, depth, dbSnpId, reference,result, qMaxScore
                            string variantLine =
                                String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}",
                                              variantCFDx.VariantStatsForReport.Type,
                                              variantCFDx.VariantStatsForReport.Position,
                                              variantCFDx.VariantStatsForReport.ChrName,
                                              frequency,
                                              depth,
                                              variantCFDx.VariantStatsForReport.CDNAName,
                                              variantCFDx.VariantStatsForReport.ProteinName,
                                              variantCFDx.VariantStatsForReport.DbSNP,
                                              variantCFDx.VariantStatsForReport.ReferenceGenotype.Equals(polyTGReference) ?
                                                shortennedPolyTG : variantCFDx.VariantStatsForReport.ReferenceGenotype,
                                              variantCFDx.VariantStatsForReport.Call,
                                              variantCFDx.VariantStatsForReport.Interpretation);
                            strBuilderReport.AppendLine(variantLine);
                        }
                        strBuilderReport.AppendLine(String.Format("Coordinates Not Called:\t{0}", sample.CoordNotCalled));
                    }
                    else
                    {
                        strBuilderReport.AppendLine(SampleFailedToken);
                    }
                }
                else //carrier report
                {
                    strBuilderReport.AppendLine(SampleCarrierColumnHeaders);

                    if (sample.CallRate >= CarrierScreenVariantCallRateThreshold)
                    {
                        //yzhao: collect all variants in the sample, apply clinical reporting logic, then write to file
                        List<VariantScreenResults> sampleVariantScreenResults = new List<VariantScreenResults>();
                        bool hasI506V = false;
                        bool hasI507V = false;
                        bool hasF508C = false;

                        foreach (VariantScreenResults variantScreenResults in VariantsScreenResultList)
                        {
                            if (variantScreenResults.SampleID.Equals(key))
                            {
                                sampleVariantScreenResults.Add(variantScreenResults);
                                if (variantScreenResults.MutationName.Equals("I506V"))
                                    hasI506V = true;
                                if (variantScreenResults.MutationName.Equals("I507V"))
                                    hasI507V = true;
                                if (variantScreenResults.MutationName.Equals("F508C"))
                                    hasF508C = true;
                            }
                        }

                        bool conditionalVariantAbsent = false;
                        if (sampleVariantScreenResults.Count > 0)
                        {
                            foreach (VariantScreenResults variantScreenResults in sampleVariantScreenResults)
                            {
                                // Not used right now, per FDA request
                                //string confidenceValue = variantScreenResults.ConfidenceValue == 0
                                //                             ? String.Empty
                                //                             : variantScreenResults.ConfidenceValue.ToString(
                                //                                 CultureInfo.InvariantCulture);

                                string variantLine = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}",
                                                                   variantScreenResults.MutationName,
                                                                   variantScreenResults.Type,
                                                                   variantScreenResults.RsID,
                                                                   variantScreenResults.Region,
                                                                   variantScreenResults.GenomicLocation,
                                                                   variantScreenResults.CDNAName,
                                                                   variantScreenResults.ProteinName,
                                                                   variantScreenResults.Result,
                                                                   variantScreenResults.Interpretation
                                                                  );
                                strBuilderReport.AppendLine(variantLine);

                                if ((variantScreenResults.MutationName.Equals("F508del") && variantScreenResults.Result.Equals("HOM")) ||
                                    (variantScreenResults.MutationName.Equals("I507del") && variantScreenResults.Result.Equals("HOM")))
                                    if (hasI506V == false && hasI507V == false && hasF508C == false)
                                        conditionalVariantAbsent = true;
                            }
                        }
                        else
                        {
                            strBuilderReport.AppendLine("No panel variants detected");
                        }

                        if (conditionalVariantAbsent == true)
                            strBuilderReport.AppendLine("I506V, I507V, F508C alleles are not present in the sample");
                        //foreach (VariantScreenResults variantScreenResults in VariantsScreenResultList)
                        //{
                        //    if (!variantScreenResults.SampleID.Equals(key))
                        //        continue;

                        //    string confidenceValue = variantScreenResults.ConfidenceValue == 0
                        //                                 ? String.Empty
                        //                                 : variantScreenResults.ConfidenceValue.ToString(
                        //                                     CultureInfo.InvariantCulture);

                        //    string variantLine = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}",
                        //                                       variantScreenResults.MutationName,
                        //                                       variantScreenResults.Type,
                        //                                       variantScreenResults.RsID,
                        //                                       variantScreenResults.Region,
                        //                                       variantScreenResults.GenomicLocation,
                        //                                       variantScreenResults.CDNAName,
                        //                                       variantScreenResults.ProteinName,
                        //                                       variantScreenResults.Result,
                        //                                       confidenceValue,
                        //                                       variantScreenResults.Interpretation
                        //        );
                        //    strBuilderReport.AppendLine(variantLine);
                        //}
                    }
                    else
                    {
                        strBuilderReport.AppendLine(SampleFailedToken);
                    }
                }
                strBuilderReport.AppendLine();
            }

            // If nothing is wrong, write to file

            FileInfo fi;
            string tempFileName = string.Empty;
            DateTime orgCreationTime = DateTime.Now;
            if (File.Exists(ReportPath))
            {
                fi = new FileInfo(ReportPath);
                tempFileName = ReportPath + "." + DateTime.Now.ToString("MMddyy-HHmmss");
                orgCreationTime = fi.CreationTime;
                File.Copy(ReportPath, tempFileName); // copy file as backup
            }

            using (StreamWriter writer = new StreamWriter(ReportPath))
            {
                try
                {
                    writer.Write(strBuilderReport.ToString());
                    writer.Close();
                }
                catch (Exception ex)
                {
                    // In case anything went wrong, reverse the sample sheet back
                    if (File.Exists(ReportPath))
                    {
                        fi = new FileInfo(ReportPath);
                        DateTime creationTime = fi.CreationTime;
                        if (DateTime.Compare(orgCreationTime, creationTime) != 0)
                            File.Copy(tempFileName, ReportPath);
                    }
                    throw new Exception(
                        String.Format(
                            "Error saving Report file {0}.  Ensure existing Report is closed. Error Details: {1}",
                            ReportPath, ex.Message), ex);
                }
            }
        }




        private void GetCfReportResults()
        {
            if (File.Exists(ReportPath))
            {
                int numClinicalDataCols = 10, numCf139DataCols = 9;

                using (StreamReader reader = new StreamReader(ReportPath))
                {
                    string line;
                    string currentSample = String.Empty;

                    while ((line = reader.ReadLine()) != null)
                    {
                        string[] tokens = line.Split('\t');
                        if (tokens.Length >= 2 && tokens[0].Equals(SampleIdToken))
                        {
                            currentSample = tokens[1];

                            if (!CfSamples.ContainsKey(currentSample))
                                CfSamples.Add(currentSample, new CFSample());
                        }
                        else if (tokens.Length >= 2 && tokens[0].Equals(RunIdToken))
                        {
                            ReportRunId = tokens[1];
                        }
                        else if (tokens.Length >= 2 && tokens[0].Equals(RunDateToken))
                        {
                            ReportRunDate = tokens[1];
                        }
                        else if (tokens.Length >= 2 && tokens[0].Equals(SampleControlToken) &&
                                 !String.IsNullOrEmpty(currentSample))
                        {
                            CfSamples[currentSample].Control = tokens[1];
                        }
                        else if (tokens.Length >= 2 && tokens[0].Equals(SampleCommentToken) &&
                                 !String.IsNullOrEmpty(currentSample))
                        {
                            CfSamples[currentSample].Comment = tokens[1];
                        }
                        else if (tokens.Length >= 2 && tokens[0].Equals(PerformanceToken) &&
                                 !String.IsNullOrEmpty(currentSample))
                        {
                            CfSamples[currentSample].Performance = tokens[1];
                        }
                        else if (tokens.Length >= 2 && tokens[0].Equals(SampleNameToken) &&
                                 !String.IsNullOrEmpty(currentSample))
                        {
                            CfSamples[currentSample].SampleName = tokens[1];
                        }
                        else if (tokens.Length >= 2 && tokens[0].Equals(SampleCallRateToken) &&
                                 !String.IsNullOrEmpty(currentSample))
                        {
                            float callRate = Single.Parse(tokens[1].Replace("%", ""));
                            CfSamples[currentSample].CallRate = callRate;
                        }
                        else if (tokens.Length >= numCf139DataCols && ReportType == CFReportType.CfCarrierScreening
                                 && !String.IsNullOrEmpty(currentSample) && VariantsScreenResultList != null)
                        {
                            int mutationTypeCol, dbSnpIDCol, geneRegionCol, genomicLocationCol, cDNANameCol;
                            int proteinNameCol, resultCol, confValCol, interpretationCol;

                            GetCf139ColIndices(tokens, out mutationTypeCol, out dbSnpIDCol, out geneRegionCol,
                                out genomicLocationCol, out cDNANameCol, out proteinNameCol, out resultCol,
                                out confValCol, out interpretationCol);

                            // keep reading this "Mutations' block until seeing the new sample,
                            // break the inside while loop
                            while ((line = reader.ReadLine()) != null)
                            {
                                string[] variantTokens = line.Split('\t');
                                if (variantTokens[0].Equals(SampleIdToken, StringComparison.InvariantCultureIgnoreCase))
                                {
                                    currentSample = variantTokens[1];

                                    if (!CfSamples.ContainsKey(currentSample))
                                        CfSamples.Add(currentSample, new CFSample());

                                    break; // end of "Mutation" block of this sample
                                }
                                if (variantTokens[0].Equals(SampleFailedToken, StringComparison.InvariantCultureIgnoreCase))
                                {
                                    break; // end of "Mutation" block of this sample
                                }
                                if (variantTokens.Length >= numCf139DataCols)
                                {
                                    int confidence = -1;
                                    if (confValCol != -1)
                                        Int32.TryParse(variantTokens[confValCol], out confidence);

                                    VariantScreenResults variant = new VariantScreenResults
                                    {
                                        SampleID = currentSample,
                                        MutationName = variantTokens[0],
                                        Type = mutationTypeCol != -1 ? variantTokens[mutationTypeCol] : string.Empty,
                                        RsID = dbSnpIDCol != -1 ? variantTokens[dbSnpIDCol] : string.Empty,
                                        Region = geneRegionCol != -1 ? variantTokens[geneRegionCol] : string.Empty,
                                        GenomicLocation = genomicLocationCol != -1 ? variantTokens[genomicLocationCol] : string.Empty,
                                        CDNAName = cDNANameCol != -1 ? variantTokens[cDNANameCol] : string.Empty,
                                        ProteinName = proteinNameCol != -1 ? variantTokens[proteinNameCol] : string.Empty,
                                        Result = resultCol != -1 ? variantTokens[resultCol] : string.Empty,
                                        Interpretation = interpretationCol != -1 ? variantTokens[interpretationCol] : string.Empty,
                                        ConfidenceValue = confidence
                                    };

                                    VariantsScreenResultList.Add(variant);
                                }
                            }
                        }
                        else if (tokens.Length >= numClinicalDataCols && ReportType == CFReportType.CfDiagnostic
                            && !String.IsNullOrEmpty(currentSample))
                        {
                            int corCol, proteinNameCol, interpretationCol, confValCol, resultCol;
                            int refCol, dbSNPCol, cDNANameCol, depthCol, freqCol, chrCol;
                            // Is this the only way to detecte the col header indices? It has to do for each and all samples?
                            GetDxColIndices(tokens, out corCol, out chrCol, out freqCol,
                                            out depthCol, out cDNANameCol,
                                            out dbSNPCol, out refCol, out resultCol,
                                            out confValCol, out interpretationCol,
                                            out proteinNameCol);

                            // keep reading this "Variant Type' block until seeing "Coordinates Not Called:", 
                            // break the inside while loop                            
                            while ((line = reader.ReadLine()) != null)
                            {
                                string[] variantTokens = line.Split('\t');

                                if (variantTokens[0].Equals(CoordinatesNotCalledToken,
                                                            StringComparison.InvariantCultureIgnoreCase))
                                {
                                    CfSamples[currentSample].CoordNotCalled = variantTokens[1];
                                    break; // end of the "Variant Type" block of this sample
                                }
                                if (variantTokens[0].Equals(SampleFailedToken, StringComparison.InvariantCultureIgnoreCase))
                                {
                                    break;
                                }
                                if (variantTokens.Length >= numClinicalDataCols)
                                //Only process lines with 11+ tokens
                                {
                                    CfVariantType type = ConvertStringToCfVariantType(
                                            (variantTokens[0].ToLower()));

                                    int pos = -1, score = -1, depth = -1;
                                    // In case of SNP and Indel, the position is just a number, in case of PolyTG/PolyT,
                                    // position is the first number in thh range (number1 - number2)
                                    if (corCol != -1)
                                    {
                                        if (variantTokens[corCol].IndexOf('-') >= 0)
                                            Int32.TryParse(variantTokens[corCol].Split('-')[0].Trim(), out pos);
                                        else
                                            Int32.TryParse(variantTokens[corCol], out pos);
                                    }

                                    if (depthCol != -1)
                                        Int32.TryParse(variantTokens[depthCol], out depth);

                                    float freq = 0f;
                                    if (freqCol != -1)
                                        Single.TryParse(variantTokens[freqCol], out freq);

                                    if (confValCol != -1)
                                        Int32.TryParse(variantTokens[confValCol], out score);

                                    VariantCFDx variant = new VariantCFDx
                                    {
                                        SampleID = currentSample,
                                        VariantStatsForReport =
                                            new VariantStatsForReport
                                            {
                                                DbSNP = (dbSNPCol != -1) ? variantTokens[dbSNPCol] : string.Empty,
                                                CDNAName = (cDNANameCol != -1) ? variantTokens[cDNANameCol] : string.Empty,
                                                Type = type,
                                                Position = pos,
                                                ReferenceGenotype = (refCol != -1) ? variantTokens[refCol] : string.Empty,
                                                Score = score,
                                                Call = (resultCol != -1) ? variantTokens[resultCol] : string.Empty,
                                                Filter = "PASS",
                                                Depth = depth,
                                                Frequency = freq,
                                                ChrName = (chrCol != -1) ? variantTokens[chrCol] : string.Empty,
                                                //according to Joe, CF Diagnostic - the chr is always 7, other chr is off-target 
                                                Interpretation = (interpretationCol != -1) ? variantTokens[interpretationCol]
                                                    : string.Empty,
                                                ProteinName = (proteinNameCol != -1) ? variantTokens[proteinNameCol] : "N/A"
                                            }
                                    };

                                    // If  sample's callRate >= call rate threshold, add to the list
                                    // The CFDiagnostic report already took care of this
                                    VariantsDiagnosticResultList.Add(variant);
                                }
                            }
                        }
                    }
                }
            }
        }

        /// <summary>
        ///     Get the index of the column
        /// </summary>
        /// <param name="tokens"></param>
        /// <param name="corCol"></param>
        /// <param name="chrCol"></param>
        /// <param name="freqCol"></param>
        /// <param name="depthCol"></param>
        /// <param name="cDNANameCol"></param>
        /// <param name="dbSNPCol"></param>
        /// <param name="refCol"></param>
        /// <param name="resultCol"></param>
        /// <param name="confValCol"></param>
        /// <param name="interpretationCol"></param>
        /// <param name="proteinNameCol"></param>
        private static void GetDxColIndices(string[] tokens, out int corCol, out int chrCol, out int freqCol,
                                            out int depthCol, out int cDNANameCol, out int dbSNPCol, out int refCol,
                                            out int resultCol,
                                            out int confValCol, out int interpretationCol, out int proteinNameCol)
        {
            corCol = -1;
            chrCol = -1;
            freqCol = -1;
            depthCol = -1;
            cDNANameCol = -1;
            dbSNPCol = -1;
            refCol = -1;
            resultCol = -1;
            confValCol = -1;
            interpretationCol = -1;
            proteinNameCol = -1;
            for (int i = 1; i < tokens.Length; i++)
            {
                switch (tokens[i].ToLower())
                {
                    case "coordinate":
                        corCol = i;
                        break;
                    case "chromosome":
                        chrCol = i;
                        break;
                    case "frequency":
                        freqCol = i;
                        break;
                    case "depth":
                        depthCol = i;
                        break;
                    case "cdna name (hgvs)":
                        cDNANameCol = i;
                        break;
                    case "dbsnp id":
                        dbSNPCol = i;
                        break;
                    case "reference":
                        refCol = i;
                        break;
                    case "result":
                        resultCol = i;
                        break;
                    case "confidence value":
                        confValCol = i;
                        break;
                    case "interpretation":
                        interpretationCol = i;
                        break;
                    case "protein name (hgvs)":
                        proteinNameCol = i;
                        break;
                }
            }
        }

        private static void GetCf139ColIndices(string[] tokens, out int mutationTypeCol, out int dbSnpIDCol, out int geneRegionCol,
                                           out int genomicLocationCol, out int cDNANameCol, out int proteinNameCol, out int resultCol,
                                           out int confValCol, out int interpretationCol)
        {
            mutationTypeCol = -1;
            dbSnpIDCol = -1;
            geneRegionCol = -1;
            genomicLocationCol = -1;
            cDNANameCol = -1;
            proteinNameCol = -1;
            resultCol = -1;
            confValCol = -1;
            interpretationCol = -1;

            for (int i = 1; i < tokens.Length; i++)
            {
                switch (tokens[i].ToLower())
                {
                    case "mutation type":
                        mutationTypeCol = i;
                        break;
                    case "dbsnp rsid":
                        dbSnpIDCol = i;
                        break;
                    case "cftr gene region":
                        geneRegionCol = i;
                        break;
                    case "genomic location":
                        genomicLocationCol = i;
                        break;
                    case "cdna name (hgvs)":
                        cDNANameCol = i;
                        break;
                    case "protein name (hgvs)":
                        proteinNameCol = i;
                        break;
                    case "result":
                        resultCol = i;
                        break;
                    case "confidence value":
                        confValCol = i;
                        break;
                    case "interpretation":
                        interpretationCol = i;
                        break;
                }
            }
        }

        public class CFSample
        {
            public float CallRate { get; set; }
            public string CoordNotCalled { get; set; }
            public string SampleName { get; set; }
            public string Control { get; set; } // can be positive, negative, or empty
            public string Comment { get; set; }
            public string Performance { get; set; }
        }

        public class VariantCFDx
        {
            public string SampleID { get; set; }
            public VariantStatsForReport VariantStatsForReport { get; set; }
        }

        public class VariantScreenResults
        {
            public string SampleID { get; set; }
            public string MutationName { get; set; }
            public string Type { get; set; }
            public string RsID { get; set; }
            public string Region { get; set; }
            public string GenomicLocation { get; set; }
            public string CDNAName { get; set; }
            public string ProteinName { get; set; }
            public string Result { get; set; }
            public int ConfidenceValue { get; set; }
            public string Interpretation { get; set; }
        }

        public class VariantStatsForReport
        {
            public int Position { get; set; }
            public int Score { get; set; }
            public CfVariantType Type { get; set; }
            public string Call { get; set; }
            public float Frequency { get; set; }
            public int Depth { get; set; }
            public string Filter { get; set; }
            public string DbSNP { get; set; }
            public string ReferenceGenotype { get; set; }
            public CfZygosity CfZygosity { get; set; }
            public string PolySum { get; set; }
            public string AltBase { get; set; }
            public string ChrName { get; set; }
            public string CDNAName { get; set; }
            public string Interpretation { get; set; }

            public int ChangeInNumBases
            {
                get { return (AltBase.Length - ReferenceGenotype.Length); }
            }

            public string ProteinName { get; set; }
        }
    }
}