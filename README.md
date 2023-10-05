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

#### Run the pipelines
Program: stepStone - a pipeline to identify chromothripsis breakpoints and rearrangements
Version: 2.0

Usage: ./stepStone command [options]			\

Commands:                                               \
-- breakpoints		Detect breakpoints              \
-- plot			Plot depth of coverage          \
-- align		Align reads to a reference      \


#### Read Alignments 

           $ /full/path/to/stepStone/src/stepStone align           \

===Align reads to a reference for all data types:                                                                                   \
===Output a coordinate sorted bam file:                                                                                             \
           $ ./stepStone align -nodes 60 -data ccs -reads input_long.fasta/q <Input_Reference> <Output_sorted_bam>            \
           $ ./stepStone align -nodes 60 -data ont -reads input_long.fasta/q <Input_Reference> <Output_sorted_bam>            \
           $ ./stepStone align -nodes 60 -data ont-NLR -reads input_long.fasta/q <Input_Reference> <Output_sorted_bam>            \
           $ ./stepStone align -nodes 60 -data ngs-HiC -reads input_long.fasta/q <Input_Reference> <Output_sorted_bam>            \
           $ ./stepStone align -nodes 60 -data ngs-10X -reads input_long.fasta/q <Input_Reference> <Output_sorted_bam>            \
           $ ./stepStone align -nodes 60 -data ngs-SSR -reads input_long.fasta/q <Input_Reference> <Output_sorted_bam>            \
      	 	-nodes    60      - Number of CPUs requested                                                                        \
      		-data     ccs     - PacBio Hifi                                                                                     \
		-data     ont     - Oxford Nanopore Q20 or Q30                                                                      \
		-data     ont-NLR - Oxford Nanopore normal long reads (NLR)                                                         \
		-data     ngs-HiC - Hi-C reads                                                                                      \
      		-data     ngs-10X - 10X reads                                                                                       \
		-data     ngs-SSR - Standard short reads such as Illumina data                                                      \

#### Breakpoint Detection

           $ /full/path/to/stepStone/src/stepStone breakpoint           \

===Detect breakpoints with aligned, and name sorted sam,bam or cram files:
           $ stepStone breakpoint -data ccs -bam my-sorted.bam <output_breakpoints.vcf>           \
           $ stepStone breakpoint -data ont -bam my-sorted.bam <output_breakpoints.vcf>           \
           $ stepStone breakpoint -data ngs-HiC -bam my-sorted.bam <output_breakpoints.vcf>           \
           $ stepStone breakpoint -data ngs-10x -bam my-sorted.bam <output_breakpoints.vcf>           \
           $ stepStone breakpoint -data ngs-SSR -bam my-sorted.bam <output_breakpoints.vcf>           \
      		-data     ccs     - PacBio Hifi
		-data     ont     - Oxford Nanopore Q20 or Q30
		-data     ngs-HiC - Hi-C reads
      		-data     ngs-10X - 10X reads
		-data     ngs-SSR - Standard short reads such as Illumina data
	--- Note: the sam/bam/cram file has to be Name sorted! ---
	--- If not read name sorted, please do the following:  ---
	--- samtools sort -@ 60 -n your.bam new.bam ---

===Detect breakpoints with fasta/fastq long read files:
	Usage: ./stepStone breakpoint -nodes 60 -data ccs/ont -reads input_long.fasta/q <Input_Reference> <breakpoints.vcf>

#### Coverage Plots 

           $ /full/path/to/stepStone/src/stepStone plot           \

===Plot depth of coverage for all data types:
===Input a coordinate sorted bam file
===Output a tmp directory containing coverage images for 23 chromosomes chr{1,22,X}
	Usage: ./stepStone plot -bam /myspace/desk/test-sorted.bam -sample cancer

===Without noise reduction:
	Usage: ./stepStone plot -bam /myspace/desk/test-sorted.bam -sample cancer -denoise 0


