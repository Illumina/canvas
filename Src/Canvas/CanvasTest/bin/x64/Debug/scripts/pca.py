from sys import argv
from pcalib import *

rlogPath = argv[1]
samplesPath = argv[2]
outputFolder = argv[3]

# Read samples.txt
sampleID2Group = ReadSampleID2Group(samplesPath, groupColumnIdx=2)

#
# PCA
#
sampleIDs, genes, rlogMatrix = ReadRLog(rlogPath)

PerformPCA(sampleIDs, genes, rlogMatrix, sampleID2Group, outputFolder)