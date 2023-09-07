# stepStone v2.0
A pipeline for identification of chromothripsis breakpoints and cancer rearrangements

### Download and Compile:

    $ git clone  https://github.com/wtsi-hpag/stepStone.git 
    $ cd stepStone 
    $ bash install.sh
		
If everything compiled successfully you must see the final comment: 
		"Congrats: installation successful!"		

#### External packages
The genome aligner BWA-mem2 (https://github.com/bwa-mem2/bwa-mem2), minimap2 (https://github.com/lh3/minimap2) and samtools (http://www.htslib.org) are downloaded and compiled by stepStone.

### Run the pipelines
Program: stepStone - a pipeline to identify chromothripsis breakpoints and rearrangements
Version: 2.0

Usage: ./stepStone command [options]			\

Commands:                                               \
-- breakpoints		Detect breakpoints              \
-- plot			Plot depth of coverage          \
-- align		Align reads to a reference      \

===Detect breakpoints with aligned, and name sorted sam,bam or cram files:     \

           $ /full/path/to/stepStone/src/stepStone breakpoint -nodes 60        \
		-data ccs/ont/ont-NLR/ngs-HiC/ngs-10X/ngs-SSR                  \
		-bam /myspace/desk/test-sorted.bam <output_breakpoints.vcf>    \
		--- Note: the sam/bam/cram file has to be Name sorted! ---     \
		--- If not read name sorted, please do the following:  ---     \
		--- samtools sort -@ 60 -n your.bam new.bam ---                \

===Detect breakpoints with fasta/fastq long read files:                        \ 

           $ /full/path/to/stepStone/src/stepStone breakpoint -nodes 60 -data ccs/ont 	\
		-reads input_long.fasta/q <Input_Reference> <breakpoints.vcf>  		\

===Plot depth of coverage for all data types:                                  		\ 

           $ /full/path/to/stepStone/src/stepStone plot -bam /myspace/test-sorted.bam -sample cancer \

===Align reads to a reference for all data types:                              			\

           $ /full/path/to/stepStone/src/stepStone align -data ccs/ont/ont-NLR 			\
		-reads input_long.fasta/q <Input_Reference> <Output_sorted_bam>			\

           $ /full/path/to/stepStone/src/stepStone align -data ngs-HiC         			\
		-fq1 read_1.fq.gz -fg2 read_2.fq.gz <Input_Reference> <Output_sorted_bam>  	\

           $ /full/path/to/stepStone/src/stepStone align -data ngs-10X                     	\
		-fq1 read_1.fq.gz -fg2 read_2.fq.gz <Input_Reference> <Output_sorted_bam>  	\

           $ /full/path/to/stepStone/src/stepStone align -data ngs-SSR                     	\
		-fq1 read_1.fq.gz -fg2 read_2.fq.gz <Input_Reference> <Output_sorted_bam>  	\

		-nodes     60     - Number of CPUs requested                 			\
		-data     ccs     - PacBio Hifi                 				\
		-data     ont     - Oxford Nanopore Q20 or Q30                 			\
		-data     ont-NLR - Oxford Nanopore normal long reads (NLR)    			\
		-data     ngs-HiC - Hi-C reads                                 			\
		-data     ngs-10X - 10X reads                                 			\
		-data     ngs-SSR - Standard short reads such as Illumina data 			\


           $ /full/path/to/stepStone/src/stepStone -nodes <nodes>  \

# stepStone
# stepStone
