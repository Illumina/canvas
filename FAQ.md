Frequently asked questions
=================================

#### I've downloaded and installed the Canvas software and trying to run the tumor-normal and somatic WGS workflow.  The software calls for two required files.  One is the somatic small variant vcf file, and the other is the b-allele vcf file.  How do you generate these files?  
Canvas Somatic-WGS uses the following vriant files: 

--somatic-vcf = SNVs and small indels in the control (“cancer”) sample – used to predict purity in the tumor sample that have no copy number changes 

--b-allele-vcf = SNVs and small indels in the reference (“normal”) sample – used to determine which germline SNPs are heterozygous. These are later utilized to compute allelic ratios that are part of the Canvas model for predicting purity and ploidy. 

There are two options to generate file for --b-allele-vcf parameter 

1)	One can download dbsnp.vcf file (for the relevant genome reference)  from https://illumina.app.box.com/CanvasPublic/ and let Canvas calculate allelic ratios 

2)	Or use germline SNV callers, such as the one available in samtools http://www.htslib.org/workflow/#mapping_to_variant - the output will need to be in vcf (or vcf.gz) format 

--somatic-vcf is a required in Canvas v1.3.5 and below (and optional in the versions above). To generate them a joint tumor-normal SNVs and small indels caller will need to be run on “cancer” and “normal” bams, for example Strelka https://sites.google.com/site/strelkasomaticvariantcaller/ . One could also omit this parameter in Canvas above v1.3.5 or provide a dummy VCF (header only) in Canvas v1.3.5 and below.


