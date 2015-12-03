import argparse
import os
import sys
import logging
import numpy as np
from math import log,isnan

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


DELIM = "\t"

def safe_float(value):
        try:
                return float(value)
        except:
                return 0.0

def ReadTableData(inputFile, keyFields, outputFields):
        result = {}
        inputHandle = open(inputFile)
        header = inputHandle.readline().rstrip().split(DELIM)
        for keyField in keyFields:
                if keyField not in header:
                        raise Exception("Input file " + inputFile + " does not contain header field " + keyField)
        for outputField in outputFields:
                if outputField not in header:
                        raise Exception("Input file " + inputFile + " does not contain header field " + outputField)
        keyIndices = [header.index(keyField) for keyField in keyFields]
        outputIndices = [header.index(outputField) for outputField in outputFields]

        for line in inputHandle:
                content = line.rstrip().split(DELIM)
                for idx in keyIndices + outputIndices:
                        if idx >= len(content):
                                raise Exception("Improperly formatted line in input file " + inputFile + "(" + line + ")")
                key = tuple([content[keyIdx] for keyIdx in keyIndices])
                output = [content[outputIdx] for outputIdx in outputIndices]
                if key in result:
                        raise Exception("Repeated key (" + DELIM.join(key) + ") in input file " + inputFile)
                result[key] = output
        inputHandle.close() 
        return result


def PruneData(inputData, FPKM_Threshold, cvThreshold):
        result = {}
        for key in inputData:
                (FPKM, FPKM_conf_lo, FPKM_conf_hi) = map(safe_float,inputData[key])
                if FPKM < FPKM_Threshold or FPKM <= 0:
                        continue
                if abs(FPKM_conf_hi - FPKM_conf_lo)/FPKM > cvThreshold:
                        continue
                result[key] = (FPKM, FPKM_conf_lo, FPKM_conf_hi)
        return result

def CalculateCorrelation(FPKM_File1, FPKM_File2, cvThreshold):
        data1 = ReadTableData(FPKM_File1, ["tracking_id", "locus",], ["FPKM", "FPKM_conf_lo", "FPKM_conf_hi"])
        data2 = ReadTableData(FPKM_File2, ["tracking_id", "locus",], ["FPKM", "FPKM_conf_lo", "FPKM_conf_hi"])

        data1 = PruneData(data1, 1.0, cvThreshold)
        data2 = PruneData(data2, 1.0, cvThreshold)

        array1 = []
        array2 = []

        for key in data1:
                if key in data2:
                        array1.append(log(data1[key][0]))
                        array2.append(log(data2[key][0]))

        if len(array1) < 50:
                return 0.0 
        return np.corrcoef(array1, array2)[0,1] 

parser = argparse.ArgumentParser(description="Plot correlation matrix with dendrogram")
parser.add_argument("--correlation-path", dest = "correlationPath", required = True, help = "Location of the correlation matrix file in CSV")
parser.add_argument("--grid-line", dest="gridLine", default=False, action="store_true", help="Plot grid lines")
options = parser.parse_args()

#
# validate the inputs
#
options.correlationPath = os.path.abspath(options.correlationPath)
if not os.path.exists(options.correlationPath):
        sys.stderr.write("Correlation matrix file does not exist\n")
        sys.exit(-1)

# create logger
logger = logging.getLogger('plot rna correlation')
logger.setLevel(logging.DEBUG)

# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)

# create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# add formatter to ch
ch.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)

#
# Read the correlation matrix
#
correlations = {}
inp = open(options.correlationPath, 'r')
colNames = []
for colName in inp.readline().rstrip('\r\n').split(','):
        if len(colName) > 1 and colName[0] == '"' and colName[-1] == '"':
                colName = colName[1:-1]
        if not colName: continue
        colNames.append(colName)
for line in inp.readlines():
        cells = line.rstrip('\r\n').split(',')
        rowName = cells[0]
        if len(rowName) > 1 and rowName[0] == '"' and rowName[-1] == '"':
                rowName = rowName[1:-1]
        cells = cells[1:]
        for i in xrange(len(cells)):
                colName = colNames[i]
                replicatePair = tuple(sorted([rowName, colName]))
                if replicatePair not in correlations:
                        correlations[replicatePair] = float(cells[i])
inp.close()
replicates = colNames

print correlations

from scipy.cluster.hierarchy import dendrogram, linkage

distanceArray = []
for replicateIdx in xrange(len(replicates)):
        for replicateIdy in xrange(replicateIdx + 1, len(replicates)):
                distanceArray.append(1 - correlations[tuple(sorted([replicates[replicateIdx], replicates[replicateIdy]]))]) 

#
# Figure out the ordering of leaves
# in the final clustering
# (this will determine how we actually plot
# the heat map)
#

#
# We don't need to do any re-ordering for
# the single sample case
#
if len(replicates) > 1:
        distanceArray = np.array(distanceArray)
        linkage_matrix = linkage(distanceArray, 'average')
        clustering_result = dendrogram(linkage_matrix, no_plot = True)

        replicates = [replicates[idx] for idx in clustering_result["leaves"]]
        replicates.reverse()


numberReplicates = len(replicates)


#
# Figure out the dimensions of the
# plot
#
MAX_DISPLAY_CHARS = 20 
MAX_DISPLAY_CHARS = min(max(map(len, replicates)),MAX_DISPLAY_CHARS)
H_PAD = 0.25
V_PAD = 0.25

legendHeight = 0.20
legendWidth = 1.8 

scatterHeight = 0.25*numberReplicates
scatterWidth =  0.25*numberReplicates
textWidth = 0.10 * MAX_DISPLAY_CHARS 

figureWidth = 2*scatterWidth + textWidth + 4*H_PAD
figureHeight = scatterHeight + 3 * V_PAD + legendHeight 

figure = plt.figure(figsize = (figureWidth, figureHeight), facecolor = 'white')
plt.axes([H_PAD/figureWidth,(2 * V_PAD + legendHeight)/figureHeight,scatterWidth/figureWidth,scatterHeight/figureHeight])

#
# Generate the heat map
#
correlationArray = np.zeros([len(replicates), len(replicates)])
for replicate in replicates:
        replicateIdx = replicates.index(replicate)
        correlationArray[replicateIdx,replicateIdx] = 1.0
for replicatePair in correlations:
        if not isnan(correlations[replicatePair]):
                idx1 = replicates.index(replicatePair[0])
                idx2 = replicates.index(replicatePair[1])
                correlationArray[idx1,idx2] = correlations[replicatePair]
                correlationArray[idx2,idx1] = correlations[replicatePair]

mymap = matplotlib.colors.LinearSegmentedColormap(
        "mycustom",
        {"red" :  [(0.0, 1.0, 1.0),
                   (0.33, 0.0, 0.0),
                   (0.66, 1.0, 1.0),
                   (1.0, 1.0, 1.0)],
        "green" :[(0.0, 1.0, 1.0),
                  (0.33, 0.5, 0.5),
                  (0.66, 1.0, 1.0),
                  (1.0, 0.0, 0.0)],
        "blue" : [(0.0, 1.0, 1.0),
                  (0.33, 0.0, 0.0),
                  (0.66, 0.0, 0.0),
                  (1.0, 0.0, 0.0)]})

image = plt.imshow(correlationArray, cmap=mymap, interpolation = "none", aspect='auto')
if options.gridLine:
	plt.xticks([float(i) - 0.5 for i in xrange(1, len(replicates))],
			   ['' for i in xrange(1, len(replicates))])
	plt.yticks([float(i) - 0.5 for i in xrange(1, len(replicates))],
			   ['' for i in xrange(1, len(replicates))])
	plt.grid(b=True, linestyle='-')
else:
	plt.xticks([])
	plt.yticks([])

#
# List the names of the replicates
#
ax = plt.axes([(scatterWidth + 2 * H_PAD)/figureWidth,(2 * V_PAD + legendHeight)/figureHeight,textWidth/figureWidth,scatterHeight / figureHeight])
for replicateIdx in xrange(len(replicates)):
        stepSize = 1/float(len(replicates))
        yPos = 1.0 - replicateIdx * stepSize - (stepSize)/2
        xPos = 0.0
        textValue = replicates[replicateIdx]
        if len(textValue) > MAX_DISPLAY_CHARS:
                textValue = textValue[0:(MAX_DISPLAY_CHARS-3)] + "..."
        ax.text(xPos, yPos, textValue, verticalalignment = "center")
ax.set_frame_on(False)
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)

#figure.text(0.5, 0.95, "A Title", horizontalalignment="center")

#
# Show the dendrogram
#
if len(replicates) > 1:
        ax = plt.axes([(scatterWidth + textWidth + 3 * H_PAD)/figureWidth, (2 * V_PAD + legendHeight)/figureHeight,scatterWidth/figureWidth,scatterHeight / figureHeight])
        dendrogram(linkage_matrix,
                        orientation = "left",
                        color_threshold = 0)
               
ax.set_frame_on(False)
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False) 

norm = matplotlib.colors.Normalize(-1.0, 1.0)
image.set_norm(norm)

colorBarAxes = plt.axes([(figureWidth/2.0 - legendWidth/2.0)/figureWidth, V_PAD/figureHeight, legendWidth/figureWidth, legendHeight/figureHeight])
colorBarAxes.set_title("correlation", fontsize = 10)
figure.colorbar(image, colorBarAxes, orientation='horizontal', ticks = [-1.0,1.0])
plotFilePrefix = os.path.splitext(options.correlationPath)[0]

maxPixels = 20000
dpi = min(maxPixels / max(figure.get_size_inches()), 600)
for format in ['png', 'pdf']:
        figure.savefig(plotFilePrefix + '.' + format, dpi=dpi, format=format)
