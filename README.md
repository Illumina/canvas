Canvas Copy Number Variant Caller
=================================

Canvas is a tool for calling copy number variants (CNVs) from human DNA sequencing data.  It can work either with germline data, or paired tumor/normal samples.  Its primary input is aligned reads (in .bam format), and its primary output is a report (in a .vcf file) giving the copy number status of the genome.  

Canvas is used as the copy number caller in the Isaac Whole Genome Sequencing workflow in BaseSpace (https://basespace.illumina.com), and in HiSeq Analysis Software (HAS) (http://support.illumina.com/sequencing/sequencing_software/hiseq-analysis-software.html).  

Canvas is written in C# and runs either under a recent version of Mono (e.g. 3.10.0), .NET 4.5.1 and .NET Core 1.1.

For more information about Canvas and the algorithms it uses see the [software design document][SDD].

[SDD]:SoftwareDesignDocument.pdf

License
-------

Copyright (c) 2013-2017 Illumina, Inc. All rights reserved.

This software is provided under the terms and conditions of the GNU GENERAL PUBLIC LICENSE Version 3

You should have received a copy of the GNU GENERAL PUBLIC LICENSE Version 3 along with this program. If not, see https://github.com/illumina/licenses/.

Canvas includes several third party packages provided under other open source licenses, please see [COPYRIGHT.txt](COPYRIGHT.txt) for additional details.

[COPYRIGHT.txt]:COPYRIGHT.txt
 
Build instructions
------------------

### Binaries:
It is recommended to start from one of the [binary distributions on the Canvas releases page][releases] if a suitable version is available.  Executables can be run either under .NET Core or mono. .NET Core environment is recommended as it provides higher speed and lower RAM usage.

[releases]:https://github.com/Illumina/canvas/releases

### Source code organization:
Canvas consists of several projects all built from one solution file (Src/Canvas/Canvas/Canvas.sln).  The main Canvas project is a command line tool for launching the various workflows. Additionally, there are projects for each Canvas module - e.g. CanvasBin counts coverage for each bin, CanvasSomaticCaller makes CNV calls for tumor/normal data - as well as some shared libraries with utility functions (math functions, file I/O for various formats, etc.)  

### Compiling from source
Open the solution file (Canvas.sln) using Visual Studio 2013 and above, and build the main solution configuration (x64 + Release).  The managed code can be run on a Windows system or on a Linux system using Mono or .NET Core.  The FileCompression library (unmanaged c++ code) can be rebuilt from source under Linux, or the prebuilt binary libFileCompression.so can be used.

### Operating System Guidelines

#### Linux
Canvas was tested under the following Linux distributions:
- CentOS 5, 6 (Mono 3.10.0, Mono 4.0.2, .NET Core 1.1)
- Ubuntu 14.04 (Mono 4.0.2, .NET Core 1.1)
- Ubuntu 16.04.4 (Mono 4.2.3, .NET Core 1.1)

Other Linux distributions and other recent Mono versions are likely to work as well but have not been explicitly tested.

#### Windows
Canvas is known to run on Windows 7 or Windows 8 systems using .NET 4.5.1

Run instructions
------------------

Canvas can be run on a variety of sequencing inputs. See the help information from the Canvas.exe command line executable for the supported workflows and required input files:

$Canvas.exe --help  (or Canvas.dll --help under .NET Core)
Canvas Canvas 1.25.0.49+master Copyright © Illumina 2017-03-23
Usage: Canvas.exe [MODE] [OPTIONS]+  
Available modes:  
- Germline-WGS - CNV calling of a germline sample from whole genome sequencing data
- Somatic-Enrichment - CNV calling of a somatic sample from targeted sequencing data
- Somatic-WGS - CNV calling of a somatic sample from whole genome sequencing data
- Tumor-normal-enrichment - CNV calling of a tumor/normal pair from targeted sequencing data
- SmallPedigree-WGS - CNV calling of a small pedigree from whole genome sequencing data
       
Options:  
  -h, --help                 show this message and exit  
  -v, --version              print version and exit  

#### Reference genome
The required input files for Human reference genome builds GRCh37, hg19, and GRCh38 can be downloaded from S3 http://canvas-cnv-public.s3.amazonaws.com/. You can use the AWS CLI to download from S3:

http://docs.aws.amazon.com/cli/latest/userguide/installing.html
http://docs.aws.amazon.com/cli/latest/reference/s3/sync.html

When using a custom reference genome the equivalent files need to be created. Use the FlagUniqueKmers project to generate the annotated fasta file (kmer.fa) for a custom reference genome. 

## Installation
The easiest way to install Canvas is to use the latest pre-copiled binaries from [releases]:https://github.com/Illumina/canvas/releases (just download and uncopress). 
### .Net Core 
To install .NET Core follow the instructions at https://www.microsoft.com/net/core#linuxubuntu
The .Net core runtime can be download here:
https://www.microsoft.com/net/download/core#/runtime
### Mono
Exact installation of mono environment depends on OS, below is an installation example for Ubuntu:
#### Compiling mono from source
```
mkdir mono-4.0.2_source
wget http://download.mono-project.com/sources/mono/mono-4.0.2.5.tar.bz2
tar xf mono-4.0.2.5.tar.bz2
cd mono-4.0.2
mkdir /home/ubuntu/mono-4.0.2
./configure --prefix=/home/ubuntu/mono-4.0.2
sudo apt-get update
sudo apt-get install gcc
sudo apt-get install g++
sudo apt-get install gettext
sudo apt-get install automake
sudo apt-get install libtool
./autogen.sh --prefix=/home/ubuntu/mono-4.0.2 --with-large-heap=yes --enable-parallel-mark --with-sgen=yes
```
#### Installing binaries
```
sudo apt-get install mono-runtime
sudo apt-get install mono-complete
```

## DEMO (SmallPedigree-WGS workflow) 
Here we provide an example on how to run Canvas SPW (Small Pedigree Workflow) on a simulated trio (bam files of 60x coverage) and then using EvaluateCNV (under Tools) to estimate performance metrics. This demo will work with the *Canvas release v1.25* and above. Amazon AWS m4.4xlarge instance was used to create this demo. It is recommended that the amount of RAM per core is 4G. More information on input options and output formats can be found on the [canvas wiki][wiki] and [software design document][SDD].

[wiki]:https://github.com/Illumina/canvas/wiki

#### Data and binaries
1. Install .Net Core https://www.microsoft.com/net/core#linuxubuntu and download Canvas binary (CanvasDIR)
2. Add BaseSpace project https://basespace.illumina.com/s/f1ganFhSPsBo with simulation bams to your account (you might need to register first). 
3. Install BaseMount and load the canvas-spw project 
```
sudo bash -c "$(curl -L https://basemount.basespace.illumina.com/install/)"
mkdir /tmp/BaseSpace
basemount --scopes="Create Global, Browse Global, Create Projects, Read Global" /tmp/BaseSpace
cd /tmp/BaseSpace
```
4. This should show the following folders under canvas-spw/AppResults 
```
- bams = simulated trio bams of 60x coverage aligned with Isaac
- canvasdata = hg19 genome reference files for running Canvas (can also be downloaded from S3 http://canvas-cnv-public.s3.amazonaws.com/) 
- snvvcf = SNV vcf files to accompany bams (joint germline CNV calls using Strelka2 https://github.com/Illumina/strelka)
- simdata = bed files with simulated inherited and de novo variants
```

#### Running Canvas
1. In this example we are accessing files through basemount (Canvas should be run as user rather than sudo root). Files could also be copied to a local drive and run from there.
2. Issue the following command (output directory - /tmp/gHapMixDemo)
```
dotnet /CanvasDIR/Canvas.dll SmallPedigree-WGS --bam=/basespace/Projects/canvas/AppResults/bams/Files/father.bam --bam=/basespace/Projects/canvas/AppResults/bams/Files/mother.bam --bam=/basespace/Projects/canvas/AppResults/bams/Files/child1.bam --mother=mother --father=father --proband=child1 -r /basespace/Projects/canvas/AppResults/canvasdata/Files/kmer.fa -g /basespace/Projects/canvas/AppResults/canvasdata/Files/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta --sample-b-allele-vcf /basespace/Projects/canvas/AppResults/snvvcf/Files/Pedigree.vcf.gz -f /basespace/Projects/canvas/AppResults/canvasdata/Files/filter13.bed -o /tmp/gHapMixDemo --ploidy-vcf="/basespace/Projects/canvas/AppResults/snvvcf/Files/MultiSamplePloidy.vcf"
```
3. The runtime will depend on the number of available CPUs and whereas bam files were copied to a local drive. The run on a bare Amazon m4.4xlarge instance (16 CPUs and 64G RAM) with network I/O took 03h34m. Results are available as VCF files: either a multi-sample VCF under gHapMixDemo or single-sample equivalents under gHapMixDemo/TempCNV folders. Here we will use EvaluateCNV tool supplied with Canvas distribution to calculate various performance metrics for inherited and de novo CNVs.

#### Running EvaluateCNV
1. First, we can run EvaluateCNV to produce recall and precision metrics for inherited Canvas CNV calls using truth variant files. 
```
zcat /tmp/gHapMixDemo/TempCNV_child1/CNV.vcf.gz | grep -v ":REF:" > /tmp/gHapMixDemo/TempCNV_child1/CNV.vcf (remove REF calls)
/CanvasDIR/Tools/EvaluateCNV/EvaluateCNV.dll /ihart/BaseSpace/Projects/CanvasSPW/AppResults/simdata/Files/child1_truth.bed /tmp/gHapMixDemo/TempCNV_child1/CNV.vcf /CanvasDIR/Tools/EvaluateCNV/generic.cnaqc.excluded_regions.bed inheritedCNVs.txt 

This gives us for PASS variants: 
Recall  97.46
Precision       93.85
```

2. Next, we run a similar command but using the de novo variant truth file and a -q 20 argument to extract variants with DQ20.
```
/CanvasDIR/Tools/EvaluateCNV/EvaluateCNV.dll /ihart/BaseSpace/Projects/CanvasSPW/AppResults/simdata/Files/child1_truth.bed /tmp/gHapMixDemo/TempCNV_child1/CNV.vcf.gz /CanvasDIR/Tools/EvaluateCNV/generic.cnaqc.excluded_regions.bed -q 20 denovoCNVs.txt  

This gives us for PASS variants:
Recall  97.98
Precision       96.51
```


## DEMO (Tumor-normal-enrichment workflow)
This demo will run Canvas on exome data for HCC2218 breast carcinoma cell lines and compare results with previously curated ground truth set. The demo presumes mono runtime and that binary files were installed to WORKDIR/canvas/canvas-1.3.4_x64/. 

#### Data 
To download demo data, add BaseSpace project https://basespace.illumina.com/s/DcPnOqHmtPNB to your account (you might need to register first). The actual files can then be downloaded from the following subdirectories:
https://basespace.illumina.com/analyses/30697313/files/28317292?projectId=26760736
https://basespace.illumina.com/analyses/30697313/files/28296383?projectId=26760736
In addition to manual download, a command line basemount (https://basemount.basespace.illumina.com ) can be used for file transfer. To install basemount run
```
sudo bash -c "$(curl -L https://basemount.basespace.illumina.com/install/)"
mkdir /tmp/BaseSpace
basemount  /tmp/BaseSpace
cd /tmp/BaseSpace
```
BaseSpace files are now available under your current directory. To run demo, transfer the following files into WORKDIR/testing/files/
```
“Projects/HiSeq 2500 RR: NRC Exome (HCC1187 & HCC2218)/AppResults/HCC1187BL/Files/HCC1187BL_S1.vcf" (germline vcf)
"Projects/HiSeq 2500 RR: NRC Exome (HCC1187 & HCC2218)/AppResults/HCC2218C/Files/HCC2218C_S1.bam" (somatic bam)
"Projects/HiSeq 2500 RR: NRC Exome (HCC1187 & HCC2218)/AppResults/HCC2218C/Files/HCC2218C_S1.bam.bai"
"Projects/HiSeq 2500 RR: NRC Exome (HCC1187 & HCC2218)/AppResults/HCC2218BL/Files/HCC2218BL_S1.bam" (normal bam)
"Projects/HiSeq 2500 RR: NRC Exome (HCC1187 & HCC2218)/AppResults/HCC2218BL/Files/HCC2218BL_S1.bam.bai"
“Projects/HiSeq 2500 RR:  NRC\ Exome\ (HCC1187 & HCC2218)/AppSessions/Isaac Enrichment 11|24|2015 9:23:23/AppResults.28295376.HCC1187BL/Files/Additional Files/NexteraRapidCapture_Exome_TargetedRegions_v1.2Used.txt” (targeted regions)
```
#### Genome reference files  
Download hg19 genome reference files from S3 (http://canvas-cnv-public.s3.amazonaws.com/) into WORKDIR/testing/hg19/.

#### Running demo
With all files copied and installed, we are now ready to run Canvas. This demo will use Tumor-normal-enrichment workflow that runs on Nextera exome data.  Execute the command below. 
```
/home/ubuntu/mono-4.0.2/bin/mono $WORKDIR/canvas/canvas-1.3.4_x64/Canvas.exe Tumor-normal-enrichment -b $WORKDIR/testing/files/HCC2218C_S1.bam --normal-bam=$WORKDIR/testing/files/HCC2218BL_S1.bam --reference=$WORKDIR/testing/hg19/kmer.fa --manifest=$WORKDIR/testing/files/NexteraRapidCapture_Exome_TargetedRegions_v1.2Used.txt -g $WORKDIR/testing/hg19/ -n HCC2218C -f $WORKDIR/testing/hg19/filter13.bed -o $WORKDIR/testing/HCC2218_v2 --b-allele-vcf=$WORKDIR/testing/files/HCC2218BL_S1.vcf --custom-parameters=CanvasBin,-m=TruncatedDynamicRange
```
CNV.vcf.gz files will be saved to HCC2218_v2 output directory. Depending on the number of available CPUs, the demo will take from few minutes to under an hour to complete.

#### Inspecting results 
Now we can test Canvas performance by using a set of previously curated HCC2218 copy number calls from whole-genome data (HCC2218Truth.vcf) and a set of repetitive or ambiguous regions (HCC2218.cnaqc.excluded_regions.bed), which are available in the TruthSets directory in S3 http://canvas-cnv-public.s3.amazonaws.com/.  The evaluation is accomplished by using EvaluateCNV; the latest binary distribution for the tool can be found in [releases]:https://github.com/Illumina/canvas/releases.  

EvaluateCNV usage info:
```
EvaluateCNV $TruthSetPath $CNV.vcf $ExcludedRegionsBed $OutputPath  [$RegionOfInterestBed]
```
In our case, given that truth files location in WORKDIR/tools/EvaluateCNV, the command is:
```
mono $WORKDIR/tools/EvaluateCNV/EvaluateCNV.exe WORKDIR/TruthSets/HCC2218Truth.vcf $WORKDIR/testing/HCC2218/CNV.vcf.gz 
$WORKDIR/TruthSets/HCC2218.cnaqc.excluded_regions.bed $WORKDIR/testing/HCC2218/EvaluateCNV.txt
```
This will save evaluation data into $WORKDIR/testing/HCC2218/EvaluateCNV.txt.
Inspecting it suggests that Canvas performed quite well in calling somatic CNV variants in HCC2218, below is an extract from the file (results obtained using Canvas 1.3.4 with the command line shown above, other versions and main/custom parameters might alter performance metrics)
```
Accuracy        92.0255
DirectionAccuracy       93.1368
Recall  88.0894
DirectionRecall 92.0237
Precision       81.3032
DirectionPrecision      84.9345
```

