Utility scripts 
-------
This folder contains various utility scripts designed to assist parsing of the Canvas output. 
The current usage includes:

``` Rscript plotSomaticModelMatrix.r [Canvas TempCNV* directory] [output directory] ```
– reads in and plots values of the best Canvas somatic model and observed coverage and allele ratios from CNVModeling.txt file in Canvas temporary output folder.

```optimizeSomaticCanvasModel.py -i TrainingSamples.txt -e CanvasPath -o OutputDir -c SomaticCanvasConfig -v EbaluateCNV Path```
– optimizes Canvas somatic model on the training data using CNV calling accuracy produced by EvaluateCNV as benchmarking set. Run ```optimizeSomaticCanvasModel.py -h``` for more information.

```Disclaimer: files in this folder are not part of the main Canvas distribution and might not be synchronized to the new Canvas releases. ```
