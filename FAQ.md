Frequently asked questions
=================================

#### I've downloaded and installed the Canvas software and am trying to run the somatic WGS workflow.  The software calls for two files (one optional and one required).  One is the somatic small variant vcf file, and the other is the b-allele vcf file.  How do you generate these files?  
Canvas Somatic-WGS uses the following variant files: 

--somatic-vcf = SNVs and small indels in the tumor sample – used to estimate purity of the tumor sample when there are not enough copy number changes. This option has no effect on called CNVs. 

--b-allele-vcf = SNVs and small indels in the normal sample – used to determine which germline SNPs are heterozygous. These are later utilized to compute allelic ratios that are part of the Canvas model for predicting purity and ploidy. 

There are two options to generate file for --b-allele-vcf parameter 

1)	One can download dbsnp.vcf file (for the relevant genome reference) from https://illumina.app.box.com/CanvasPublic/ and let Canvas calculate allelic ratios 

2)	Or use germline SNV callers, such as the one available in samtools http://www.htslib.org/workflow/#mapping_to_variant - the output will need to be in vcf (or vcf.gz) format 

--somatic-vcf is required in Canvas v1.3.5 and below and optional in later versions. To generate the somatic vcf a joint tumor-normal SNV and small indel caller will need to be run on tumor and normal bams. For example Strelka https://sites.google.com/site/strelkasomaticvariantcaller/ . One could omit this parameter in Canvas versions above v1.3.5 or provide an empty VCF (header only) in Canvas v1.3.5 and below.


