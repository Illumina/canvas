# Canvas Software Copyright (c) 2013-2015 Illumina, Inc. All rights reserved.

This software is provided under the terms and conditions of the GNU GENERAL PUBLIC LICENSE Version 3

You should have received a copy of the GNU GENERAL PUBLIC LICENSE Version 3 along with this program. If not, see https://github.com/illumina/licenses/.

# Overview

Canvas is a tool for calling copy number variants (CNVs) from human DNA sequencing data.  It can work either with germline data, or paired tumor/normal samples.  Its primary input is aligned reads (in .bam format), and its primary output is a report (in a .vcf file) giving the copy number status of the genome.  

Canvas is used as the copy number caller in the Isaac Whole Genome Sequencing workflow in BaseSpace (https://basespace.illumina.com), and in HiSeq Analysis Software (HAS) (http://support.illumina.com/sequencing/sequencing_software/hiseq-analysis-software.html).  

For more information on Canvas, see the document Docs/CanvasSoftwareDesignDescription.pdf for a description of Canvas and the algorithms it uses.

# Source code

Canvas consists of serveral components all built from one solution file (Canvas.sln).  There are several executables - e.g. CanvasBin counts coverage for each bin, CanvasSomaticCaller makes CNV calls for tumor/normal data - as well as some shared libraries with utility functions (math functions, file I/O for various formats, etc.)
