######################################################################
## 
## This file contains functions used by the Whole Genome RNA and
## RNA Express workflows for PCA.
##
######################################################################
import numpy
import sys
import os
import csv

# Reads mapping from sample ID to sample group (sample name)
# The file, samples.txt, is written by DeseqWrapper.WriteParamFile() or WholeGenomeRNAStatistics.WriteSampleInformation().
# The sample group column is given by name or by position or assumed to be the last column.
def ReadSampleID2Group(path, hasHeaders=True, groupColumnIdx=-1, groupColumnName=None):
    sampleID2Group = {}
    inp = open(path, 'r')
    if hasHeaders:
        headers = inp.readline().lower().rstrip('\r\n').split('\t') # skip headers
        if groupColumnName is not None:
            groupColumnIdx = headers.index(groupColumnName)
    for line in inp.readlines():
        cells = line.rstrip('\r\n').split('\t')
        sampleID2Group[cells[0]] = cells[groupColumnIdx]
    inp.close()
    return sampleID2Group

# Parses genes.rlog.csv written by the global R script
def ReadRLog(path):
    genes = []
    mat = []
    with open(path, 'r') as inp:
        csvReader = csv.reader(inp)
        headers = csvReader.next()
        for row in csvReader:
            try:
                rlog = map(float, row[1:])
            except:
                continue
            genes.append(row[0])
            mat.append(rlog)
    sampleIDs = headers[1:]
    return sampleIDs, genes, numpy.array(mat)

# Parses filtered.merged.genes.fpkm_tracking written by cuffnorm
def ReadFPKM(path):
    with open(path, 'r') as inp:
        header = inp.next().rstrip('\r\n').split('\t')
        sampleIDs = [sampleID.replace('_FPKM', '') for sampleID in header[2:]]
        genes = []
        mat = []
        for line in inp:
            cells = line.rstrip('\r\n').split('\t')
            if not cells: break
            try:
                row = map(float, cells[2:])
            except:
                continue
            genes.append(cells[0])
            mat.append(row)
    return sampleIDs, genes, numpy.array(mat)

# Writes PCA results in CSV
def WritePCA(path, pcs, percentVariance, sampleIDs, sampleGroups):
    out = open(path, 'w')
    out.write('"","%s"\n' % ('","'.join(sampleIDs), ))
    out.write('"SampleGroup","%s"\n' % ('","'.join(sampleGroups), ))
    nPCs = pcs.shape[1]
    for i in xrange(1, nPCs + 1):
        pc = pcs[:, -i]
        p = percentVariance[-i]
        out.write('"PC%s (%.2f%%)","%s"\n' % (i, p, '","'.join(map(str, pc))))
    out.close()

# Performs PCA on the input expression matrix
def PerformPCA(sampleIDs, genes, expMatrix, sampleID2Group, outputFolder):
    if len(sampleIDs) < 2:
        sys.stderr.write('Must have at least 2 samples to perform PCA.\n')
        sys.exit(1)
    sampleGroups = [sampleID2Group[sampleID] for sampleID in sampleIDs]
    mu = numpy.array([numpy.mean(expMatrix, axis=1)]).T
    centeredExpMatrix = numpy.subtract(expMatrix, mu)
    # The covariance matrix of genes: numpy.dot(centeredExpMatrix, centeredExpMatrix.T) / nSamples
    # Decomposing the covariance matrix gives the axes (the eigen vectors).
    # Projection of the samples onto the axes is needed to get the principal components.
    # Decomposing numpy.dot(centeredExpMatrix.T, centeredExpMatrix) / nSamples gives the principal components
    # as eigen vectors without the need of projection.
    nSamples = float(len(sampleIDs))
    sigma = numpy.dot(centeredExpMatrix.T, centeredExpMatrix) / nSamples # not the covariance matrix of genes
    res = numpy.linalg.eigh(sigma)
    mask = res[0] > 0
    nComponents = sum(mask)
    if nComponents < 2:
        sys.stderr.write('Not enough principal components to generate a 2-D scatter plot.\n')
    else:
        pcs = numpy.multiply(res[1][:, mask], numpy.sqrt(res[0][mask]))
        percentVariance = res[0][mask] / sum(res[0][mask]) * 100
        WritePCA(os.path.join(outputFolder, 'pca.csv'), pcs, percentVariance, sampleIDs, sampleGroups)
