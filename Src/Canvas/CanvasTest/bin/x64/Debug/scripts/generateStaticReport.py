from __future__ import division
import argparse
import logging
import sys
import os
import xml.etree.ElementTree as ElementTree 
import jinja2
from tempfile import mktemp
import shutil
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment-directory", dest="alignmentDirectory", help="Directory containing alignment results", required=True)
    parser.add_argument("--template-directory", dest="templateDirectory", help="Directory containing the report template", required=True)
    parser.add_argument("--output-html", dest="outputHTML", help="File to store HTML output", required=True)

    args = parser.parse_args()
    args.outputHTML = os.path.abspath(args.outputHTML)
    return args


def GenerateTableHTML(header, rows):
        env = jinja2.Environment(autoescape = True, loader = jinja2.FileSystemLoader(args.templateDirectory))
        template = env.get_template("html/table.html")
        return template.render(header = header, rows = rows)

def GenerateBodyHTML(HTML_Dictionary, template = "html/body.html"):
        print(args.templateDirectory)
        env = jinja2.Environment(loader = jinja2.FileSystemLoader(args.templateDirectory))
        template = env.get_template(template)
        return template.render(HTML_Dictionary)

def ReadDataSeries(filename):
        result = []
        for line in open(filename):
                result.append(line.rstrip().split(","))
        return result


def safe_int(dic, key):
        try:
                val = int(dic[key])
        except:
                return "N/A"
        return "{0:,d}".format(val)

def safe_percent(dic, key, denom):
        try:
                val = float(dic[key]) / float(denom) * 100
        except:
                return "N/A"
        return "{0:0.2f}".format(val)

def safe_float(dic, key):
        try:
                val = float(dic[key])
        except:
                return "N/A"
        return "{0:0.2f}".format(val)


def ReadMetricsFile(metricsFile):
        result = {}
        for line in open(metricsFile):
                splat = line.rstrip().split(",")
                result[splat[0]] = splat[1]
        return result

def main():
        global args
        args = parse_args()

        metricDirectories = {} 
        samplesDirectory = os.path.join(args.alignmentDirectory, "samples")
        if not os.path.isdir(samplesDirectory):
                logging.error("No samples directory found in summary directory")
                sys.exit(-1)

        for sampleDirectory in os.listdir(samplesDirectory):
                sampleDirectory = os.path.join(samplesDirectory, sampleDirectory)
                replicatesDirectory = os.path.join(sampleDirectory, "replicates")
                if not os.path.isdir(replicatesDirectory):
                        logging.error("Sample directory " + sampleDirectory + " does not have any replicates")
                        sys.exit(-1)
                for replicateDirectory in os.listdir(replicatesDirectory):
                        metricDirectory = os.path.join(replicatesDirectory, replicateDirectory, "metrics")
                        if os.path.exists(metricDirectory):
                            metricDirectories[(os.path.basename(sampleDirectory), os.path.basename(replicateDirectory))] = metricDirectory
                        else:
                            logging.error("Did not find metrics file at " + metricDirectory)
                            sys.exit(-1)

        global tempDirectory
        tempDirectory = os.path.dirname(args.outputHTML)

        HTML_Dictionary = {}

        #
        # Generate alignment summary table
        #
        totalPairs = dict()
        header = ["Sample", "Replicate", "Read Count (Pairs)", "Uniquely mapped (%)", "Multi-mapped (%)", "Average mapped length (Pairs)", 
                  "Mismatch rate per base (%)", "Total Number of Splices", "Alignments spliced (%)"]
        rows = []
        for (sample, replicate) in sorted(metricDirectories.keys()):
                metricDirectory = metricDirectories[(sample, replicate)]
                metricFile = os.path.join(metricDirectory, "alignment.csv")
                metricDictionary = ReadMetricsFile(metricFile)
                rowData = [sample, replicate]
                totalPairs[(sample, replicate)] = int(metricDictionary["Number of reads (pairs)"])
                rowData.append(safe_int(metricDictionary,"Number of reads (pairs)"))
                rowData.append(safe_float(metricDictionary,"Uniquely mapped reads (pairs) (%)"))
                rowData.append(safe_float(metricDictionary,"Multi-mapped reads (%)"))
                rowData.append(safe_float(metricDictionary,"Average mapped length (pairs)"))
                rowData.append(safe_float(metricDictionary,"Mismatch rate per base (%)"))
                rowData.append(safe_int(metricDictionary,"Number of Splices"))
                rowData.append(safe_float(metricDictionary,"Spliced Alignments (%)"))
                rows.append(rowData)
        HTML_Dictionary["AlignmentTable"] = GenerateTableHTML(header, rows)

        #
        # Generate alignment summary table
        #
        header = ["Sample", "Replicate", "Assigned exonic reads (% of Pairs)", "Unassigned (non-exonic)", "Ambiguous", "Number of genes"]
        rows = []
        for (sample, replicate) in sorted(metricDirectories.keys()):
                metricDirectory = metricDirectories[(sample, replicate)]
                metricFile = os.path.join(metricDirectory, "counts.csv")
                metricDictionary = ReadMetricsFile(metricFile)
                rowData = [sample, replicate]
                nPairs = totalPairs.get((sample, replicate), None)
                rowData.append(safe_percent(metricDictionary,"Assigned reads", nPairs))
                rowData.append(safe_percent(metricDictionary,"Unassigned reads", nPairs))
                rowData.append(safe_percent(metricDictionary,"Ambiguous reads", nPairs))
                rowData.append(safe_int(metricDictionary,"Number of genes"))
                rows.append(rowData)
        HTML_Dictionary["CountTable"] = GenerateTableHTML(header, rows)

        seqCounts = dict()
        for repl in sorted(metricDirectories.keys()):
            metricFile = os.path.join(metricDirectories[repl], "seqCounts.csv")
            seqCounts[repl] = ReadMetricsFile(metricFile) if os.path.exists(metricFile) else dict()
        seqs = sorted(set(s for s in seqCounts[repl] for repl in seqCounts))
        header = ["Sample", "Replicate"] + [s + " (% of Reads)" for s in seqs]
        rows = []
        for (sample, replicate) in sorted(seqCounts):
            rowData = [sample, replicate]
            nPairs = totalPairs.get((sample, replicate), None)
            rowData += [safe_percent(seqCounts[(sample, replicate)], seq, nPairs) for seq in seqs]
            rows.append(rowData)
        HTML_Dictionary["SeqCountTable"] = GenerateTableHTML(header, rows)

        #
        # Generate differential expression summary
        #
        differentialRoot = os.path.join(args.alignmentDirectory, "differential")
        differentialDirectories = {}
        if os.path.isdir(differentialRoot):
                for file in os.listdir(differentialRoot):
                    if "_vs_" in file:
                        (g1, _, g2) = file.partition("_vs_")
                        differentialDirectories[(g1,g2)] = os.path.join(differentialRoot, file)

        header = ["Control", "Comparison", "Total Genes", "Tested Genes", "Differential Genes", "Links"]
        rows = []
        for cmps, diffDir in differentialDirectories.items():
                metricFile = os.path.join(diffDir, "deseq.metrics.csv")
                metricDictionary = ReadMetricsFile(metricFile)
                rowData = []
                rowData.append(cmps[0])
                rowData.append(cmps[1])
                rowData.append(safe_int(metricDictionary,"genes"))
                rowData.append(safe_int(metricDictionary,"tested"))
                rowData.append(safe_int(metricDictionary,"differential"))
                links = [] 
                try:
                        outputDir = os.path.dirname(args.outputHTML)
                        prefix = cmps[0] + "_vs_" + cmps[1]
                        outputFile = prefix + "_genes_browser.html"
                        dataFile = prefix + "_genes.json"
                        maFile =  prefix + "_MA.pdf"
                        heatmapFile = prefix + "_heatmap.pdf"
                        dispFile = prefix + "_disp.pdf"
                        #if os.path.isfile(targetFile): os.remove(targetFile)
                        shutil.copy(os.path.join(diffDir, "deseq.res.json"), os.path.join(outputDir, dataFile))
                        shutil.copy(os.path.join(diffDir, "deseq.ma.pdf"), os.path.join(outputDir, maFile))
                        if os.path.isfile(os.path.join(diffDir, "deseq.heatmap.pdf")):
                            shutil.copy(os.path.join(diffDir, "deseq.heatmap.pdf"), os.path.join(outputDir, heatmapFile))
                        shutil.copy(os.path.join(diffDir, "deseq.disp.pdf"), os.path.join(outputDir, dispFile))

                        handle = open(os.path.join(os.path.dirname(args.outputHTML), outputFile), "w")
                        handle.write(GenerateBodyHTML({"control" : cmps[0], "comparison" : cmps[1], "prefix" : prefix}, "html/gene_browser.html"))
                        handle.close()
                        links.append("<a href=\"" + outputFile + "\">Browse Genes</a>")
                except Exception as e:
                    logging.warning("Failed to create gene results browser page: " + str(e)) 

                rowData.append(jinja2.Markup(",".join(links)))
                rows.append(rowData)

        if len(rows) > 0:
                HTML_Dictionary["DifferentialTable"] = "<h1> Differential Expression </h1>\n" + GenerateTableHTML(header, rows)
        else:
                HTML_Dictionary["DifferentialTable"] = ""


        handle = open(args.outputHTML, "w")
        handle.write(GenerateBodyHTML(HTML_Dictionary))
        handle.close()

        sampleReplicatePairs = ["_".join(key) for key in sorted(metricDirectories.keys())]

        heatmap = os.path.join(args.alignmentDirectory, "differential", "global", "deseq.corr.png")
        if os.path.exists(heatmap):
            shutil.copyfile(heatmap, os.path.join(os.path.dirname(args.outputHTML), "correlation.png"))
        pca = os.path.join(args.alignmentDirectory, "differential", "global", "pca.png")
        if os.path.exists(pca):
            shutil.copyfile(pca, os.path.join(os.path.dirname(args.outputHTML), "pca.png"))
        helpFile = os.path.join(args.templateDirectory, "README.html")
        if os.path.exists(helpFile):
            shutil.copyfile(helpFile, os.path.join(os.path.dirname(args.outputHTML), "README.html"))

if __name__ == "__main__":
        main()
