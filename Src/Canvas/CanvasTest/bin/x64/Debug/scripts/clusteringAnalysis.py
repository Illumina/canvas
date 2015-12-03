from pandas import DataFrame, read_csv
from numpy import log, corrcoef, array, transpose
from scipy.cluster.hierarchy import dendrogram, linkage
import sys
import argparse
from tempfile import mkstemp
import os
DELIM="\t"
TEXT_COLUMNS = 3

parser = argparse.ArgumentParser("Perform clustering analysis on TREx hits data")


parser.add_argument("--input-file", dest = "input_file", required = True, help = "Input file containing raw hits (e.g., TargetHitsPerSample_M1.tsv)")
parser.add_argument("--output-dir", dest = "output_dir", required = True, help = "Directory to hold result files")
parser.add_argument("--output-suffix", dest = "output_suffix", required = True, help = "Suffix for naming output files")


args = parser.parse_args()


try:
    open(args.input_file)
except:
    sys.stderr.write("Unable to read input file: " + args.input_file + "\n")
    sys.exit(-1)

try:
    (handle, path_name) = mkstemp(dir = args.output_dir)
    os.close(handle)
    os.remove(path_name)
except:
    sys.stderr.write("Unable to write to output directory: " + args.output_dir + "\n")
    sys.exit(-1)


# read in data from file
sample2group = {}
with open(args.input_file) as hits_handle:
    header = hits_handle.readline().rstrip().split(DELIM)
    groups = hits_handle.readline().rstrip().split(DELIM)
    for group_index in xrange(3, len(header)):
        sample2group[header[group_index]] = groups[group_index]

samples = sample2group.keys()

data_frame = read_csv(args.input_file, sep = DELIM, skiprows = 2, names = header)

# apply minimum count threshold and log-transform
numerical_data = data_frame.ix[:,TEXT_COLUMNS:]
samples = map(str, numerical_data.columns)


numerical_data = numerical_data.clip_lower(1).apply(log)

# Perform median normalization per-sample
numerical_data = numerical_data.subtract(numerical_data.median(axis = 0), axis = 1)

# Restore data to original data frame
data_frame.ix[:,TEXT_COLUMNS:] = numerical_data

if len(samples) > 1 and len(data_frame["Amplicon ID"]) > 2:
    try:
        linkage_matrix = linkage(transpose(array(numerical_data)), method = 'average', metric = 'correlation')

        # may need to do clipping on distance matrix for small floating point errors
        for idx in xrange(len(linkage_matrix)):
            linkage_matrix[idx][2] = max(0, linkage_matrix[idx][2])
            
        clustering_result = dendrogram(linkage_matrix, no_plot = True)

        leaf_index2original_index = [0] * len(samples)
        for idx in xrange(len(clustering_result["leaves"])):
            leaf_index2original_index[clustering_result["leaves"][idx]] = idx

        # for simplicity, remap samples into optimal leaf ordering
        # pretend we had original specified the samples in this ordering
        samples = [samples[idx] for idx in clustering_result["leaves"]]

        with open(os.path.join(args.output_dir, "SampleDendrogram_" + args.output_suffix + ".csv"), "w") as dendrogram_handle:
            for row in linkage_matrix:
                result = []
                if row[0] < len(samples):
                    result.append(int(leaf_index2original_index[int(row[0])]))
                else:
                    result.append(int(row[0]))
                    
                if row[1] < len(samples):
                    result.append(int(leaf_index2original_index[int(row[1])]))
                else:
                    result.append(int(row[1]))
                result.append(row[2])
                result.append(row[3])
                            
                dendrogram_handle.write(",".join(map(str, result)) + "\n")

        # re-order the data columsn to match the order of the leaves
        data_frame = data_frame.reindex_axis(list(data_frame.columns.values[:TEXT_COLUMNS]) + samples, axis = 1)
    except:
        sys.stderr.write("Warning: failed to generate sample clustering information\n")

    
if len(samples) > 2 and len(data_frame["Amplicon ID"]) > 1:
    try:
        linkage_matrix = linkage(array(numerical_data), method = 'average', metric = 'correlation')

        # may need to do clipping on distance matrix for small floating point errors
        for idx in xrange(len(linkage_matrix)):
            linkage_matrix[idx][2] = max(0, linkage_matrix[idx][2])
            
        clustering_result = dendrogram(linkage_matrix, no_plot = True)

        leaf_index2original_index = [0] * len(data_frame["Amplicon ID"])
        for idx in xrange(len(clustering_result["leaves"])):
            leaf_index2original_index[clustering_result["leaves"][idx]] = idx

        # for simplicity, remap amplicons into optimal leaf ordering
        # pretend we had original specified the amplicons in this ordering
        amplicons = [data_frame["Amplicon ID"][idx] for idx in clustering_result["leaves"]]

        with open(os.path.join(args.output_dir, "AmpliconDendrogram_" + args.output_suffix + ".csv"), "w") as dendrogram_handle:
            for row in linkage_matrix:
                result = []
                if row[0] < len(data_frame["Amplicon ID"]):
                    result.append(int(leaf_index2original_index[int(row[0])]))
                else:
                    result.append(int(row[0]))
                    
                if row[1] < len(data_frame["Amplicon ID"]):
                    result.append(int(leaf_index2original_index[int(row[1])]))
                else:
                    result.append(int(row[1]))
                result.append(row[2])
                result.append(row[3])
                            
                dendrogram_handle.write(",".join(map(str, result)) + "\n")
        # re-order the data rows to match the order of the leaves
        data_frame = data_frame.reindex_axis(clustering_result["leaves"], axis = 0)
    except:
        sys.stderr.write("Warning: failed to generate amplicon clustering information\n")


# Calculate MAD normalization per-row
# median location shift
data_frame.ix[:,TEXT_COLUMNS:] = data_frame.ix[:,3:].subtract(data_frame.ix[:,TEXT_COLUMNS:].median(axis = 1), axis = 0)

# mad scaling
data_frame.ix[:,TEXT_COLUMNS:] = data_frame.ix[:,3:].divide( data_frame.ix[:,TEXT_COLUMNS:].mad(axis = 1), axis = 0)


    
# output the sorted data
data_frame.to_csv(os.path.join(args.output_dir, "NormalizedHitsPerSample_" + args.output_suffix + ".tsv"), sep = DELIM, index = False)




