/****************************************************************************
 ****************************************************************************
 *                                                                          *
 *  Copyright (C) 2023  Genome Research Ltd.                                *
 *                                                                          *
 *  Author: Zemin Ning (zn1@sanger.ac.uk)                                   *
 *                                                                          *
 *  This file is part of steppingStone pipeline                             *
 *                                                                          *
 *  steppingStone is a free software: you can redistribute it and/or modify *
 *  it under the terms of the GNU General Public License as published by the*
 *  Free Software Foundation, either version 3 of the License, or (at your  *
 *  option) any later version.                                              *
 *                                                                          *
 *  This program is distributed in the hope that it will be useful, but     *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of              *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU        *
 *  General Public License for more details.                                *
 *                                                                          *
 *  You should have received a copy of the GNU General Public License along *
 *  with this program.  If not, see <http://www.gnu.org/licenses/>.         *
 *                                                                          *
 ****************************************************************************
 ****************************************************************************/
/****************************************************************************/

 
#include <math.h>
#include <values.h>
#include <stdio.h>
#include <netinet/in.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include "fasta.h"

#define GT '>'
#define GT4 (((((GT<<8)+GT)<<8)+GT)<<8)+GT

#define ENDS_EXTRA 0
#define PADCHAR '-'
#define MAX_N_BRG 50000 
#define MAX_N_ROW 40000 
#define Max_N_NameBase 400 
#define Max_N_NameBase2 400 
#define Max_N_Pair 100
static char bindir[2000];
static char tmpdir[2000];
static char **S_Name;
static int *insert_siz,*insert_dev,*core_list,*ctg_list,*ctg_head,*read2contig;
static int *readIndex;

/* Default parameters   */
static int n_group=0;
static int seq_len = 40000;
static int file_tag = 0;
static int run_align = 1;
static int plot_flag = 0;
static int NLR_flag = 0;
static int flag_align = 0;
static int flag_plot = 0;
static int window_size = 100;
static int flag_breakpoint = 0;
static int min_len = 3000;
static int bam_flag = 0;
static int data_flag = 1;
static int platform_tag = 1;
static int denoise_flag = 3;
static int y_hight = 180;
static int num_chr = 23;
static int n_cover = 5;

void
RunSystemCommand(char *cmd)
{
    int ret;
    if ((ret = system(cmd)) != 0)
    {
        fprintf(stderr, "Error running command: %s\n", cmd);
        exit(EXIT_FAILURE);
    }
}


int main(int argc, char **argv)
{
    int i,nSeq,args;
    char *st,*ed;
    char **cmatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch);
    void ArraySort_String2(int n,char **Pair_Name,int *brr);
    fasta *seq;
    FILE *fp,*namef,*namef2;
    int size_file;
    int m_score,n_nodes,n_debug,num_sigma;
    void Memory_Allocate(int arr);
    char tempa[2000],tempc[2000],syscmd[2000],workdir[2000],commands[2000],window_word[100];
    char file_tarseq[2000],file_Stones[2000],file_snpfrq[2000],file_datas[2000],file_bambam[2000],file_ctspec[2000];
    char file_lgread[2000],samname[500],bamname[500],toolname[500],fastname[500],datname[500];
    char file_read1[2000],file_read2[2000],fq1name[500],fq2name[500];
    char datatype[10],tagname[10],sample[50];
    int systemRet = system (syscmd);
    int systemChd = chdir(tmpdir);
    pid_t pid;

    seq=NULL;
    
    if(argc < 2)
    {
         printf("Program: stepStone - a pipeline to identify chromothripsis breakpoints and rearrangements\n");
         printf("Version: 2.0\n");
         printf("\n");
         
         printf("Usage: %s <command> [options]\n",argv[0]);
         printf("\n");
         printf("Commands:\n");
         printf("-- breakpoints		Detect breakpoints\n");
         printf("-- plot			Plot depth of coverage\n");
         printf("-- align		Align reads to a reference\n");
         printf("\n");
	 
         printf("===Detect breakpoints with aligned, and name sorted sam,bam or cram files:\n");
         printf("	Usage: %s breakpoint -nodes 60 -data ccs/ont/ont-NLR/ngs-HiC/ngs-10X/ngs-SSR -bam /myspace/desk/test-sorted.bam <output_breakpoints.vcf>\n",argv[0]);
         printf("	--- Note: the sam/bam/cram file has to be Name sorted! ---\n");
         printf("	--- If not read name sorted, please do the following:  ---\n");
         printf("	--- samtools sort -@ 60 -n your.bam new.bam ---\n\n");

         printf("===Detect breakpoints with fasta/fastq long read files:\n");
         printf("	Usage: %s breakpoint -nodes 60 -data ccs/ont -reads input_long.fasta/q <Input_Reference> <breakpoints.vcf>\n",argv[0]);
         printf("\n");
         
         printf("===Plot depth of coverage for all data types:\n");
         printf("	Usage: %s plot -nodes 60 -bam /myspace/desk/test-sorted.bam -sample cancer\n",argv[0]);
         printf("\n");
         
         printf("===Align reads to a reference for all data types:\n");
         printf("	Usage: %s align -nodes 60 -data ccs/ont/ont-NLR -reads input_long.fasta/q <Input_Reference> <Output_sorted_bam>\n",argv[0]);
         printf("	Usage: %s align -nodes 60 -data ngs-HiC -fq1 read_1.fq.gz -fq2 read_2.fq.gz <Input_Reference> <Output_sorted_bam>\n",argv[0]);
         printf("	Usage: %s align -nodes 60 -data ngs-10X -fq1 read_1.fq.gz -fq2 read_2.fq.gz <Input_Reference> <Output_sorted_bam>\n",argv[0]);
         printf("	Usage: %s align -nodes 60 -data ngs-SSR -fq1 read_1.fq.gz -fq2 read_2.fq.gz <Input_Reference> <Output_sorted_bam>\n",argv[0]);
         printf("      	 	-nodes    60      - Number of CPUs requested\n");
         printf("      		-data     ccs     - PacBio Hifi\n");
         printf("		-data     ont     - Oxford Nanopore Q20 or Q30\n");
         printf("		-data     ont-NLR - Oxford Nanopore normal long reads (NLR)\n");
         printf("		-data     ngs-HiC - Hi-C reads\n");
         printf("      		-data     ngs-10X - 10X reads\n");
         printf("		-data     ngs-SSR - Standard short reads such as Illumina data\n");
         exit(1);
    }

    m_score = 50;
    n_nodes = 30;
    n_debug = 1;

    strcpy(toolname,"minimap2");
    strcpy(sample,"UNKNOW");
    memset(tagname,'\0',10);
    strcpy(datatype,"ccs");
    nSeq=0;
    args=1;
    for(i=1;i<argc;i++)
    {
       if(!strcmp(argv[i],"-nodes"))
       {
         sscanf(argv[++i],"%d",&n_nodes);
         args=args+2;
       }
       else if(!strcmp(argv[i],"align"))
       {
         flag_align = 1;
         args=args+1;
       }
       else if(!strcmp(argv[i],"plot"))
       {
         flag_plot = 1;
         args=args+1;
       }
       else if(!strcmp(argv[i],"breakpoint"))
       {
         flag_breakpoint = 1;
         args=args+1;
       }
       else if(!strcmp(argv[i],"-data"))
       {
         run_align = 0;
         sscanf(argv[++i],"%s",datatype);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-fq1"))
       {
         sscanf(argv[++i],"%s",fq1name);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-fq2"))
       {
         sscanf(argv[++i],"%s",fq2name);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-sample"))
       {
         sscanf(argv[++i],"%s",sample); 
         args=args+2;
       }       
       else if(!strcmp(argv[i],"-bam"))
       {
         run_align = 0;
         bam_flag = 1;
         sscanf(argv[++i],"%s",datname);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-window"))
       {
         sscanf(argv[++i],"%d",&window_size);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-chrnum"))
       {
         sscanf(argv[++i],"%d",&num_chr);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-platform"))
       {
         run_align = 1;
         sscanf(argv[++i],"%d",&platform_tag);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-denoise"))
       {
         sscanf(argv[++i],"%d",&denoise_flag);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-hight"))
       {
         sscanf(argv[++i],"%d",&y_hight);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-fasta"))
       {
         run_align = 1;
         sscanf(argv[++i],"%s",fastname);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-reads"))
       {
         run_align = 1;
         sscanf(argv[++i],"%s",fastname);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-help"))
       {
         printf("Program: stepStone - a pipeline to identify chromothripsis breakpoints and rearrangements\n");
         printf("Version: 2.0\n");
         printf("\n");
         
         printf("Usage: %s <command> [options]\n",argv[0]);
         printf("\n");
         printf("Commands:\n");
         printf("-- breakpoints		Detect breakpoints\n");
         printf("-- plot			Plot depth of coverage\n");
         printf("-- align		Align reads to a reference\n");
         printf("\n");
	 
         printf("===Detect breakpoints with aligned, and name sorted sam,bam or cram files:\n");
         printf("	Usage: %s breakpoint -nodes 60 -data ccs/ont/ngs-HiC/ngs-10X/ngs-SSR -bam /myspace/desk/test-sorted.bam <output_breakpoints.vcf>\n",argv[0]);
         printf("	--- Note: the sam/bam/cram file has to be Name sorted! ---\n");
         printf("	--- If not read name sorted, please do the following:  ---\n");
         printf("	--- samtools sort -@ 60 -n your.bam new.bam ---\n\n");

         printf("===Detect breakpoints with fasta/fastq long read files:\n");
         printf("	Usage: %s breakpoint -nodes 60 -data ccs/ont -reads input_long.fasta/q <Input_Reference> <breakpoints.vcf>\n",argv[0]);
         printf("\n");
         
         printf("===Plot depth of coverage for all data types:\n");
         printf("	Usage: %s plot -bam /myspace/desk/test-sorted.bam -sample cancer\n",argv[0]);
         printf("\n");
         
         printf("===Align reads to a reference for all data types:\n");
         printf("	Usage: %s align -nodes 60 -data ccs/ont/ont-NLR -reads input_long.fasta/q <Input_Reference> <Output_sorted_bam>\n",argv[0]);
         printf("	Usage: %s align -nodes 60 -data ngs-HiC -fq1 read_1.fq.gz -fq2 read_2.fq.gz <Input_Reference> <Output_sorted_bam>\n",argv[0]);
         printf("	Usage: %s align -nodes 60 -data ngs-10X -fq1 read_1.fq.gz -fq2 read_2.fq.gz <Input_Reference> <Output_sorted_bam>\n",argv[0]);
         printf("	Usage: %s align -nodes 60 -data ngs-SSR -fq1 read_1.fq.gz -fq2 read_2.fq.gz <Input_Reference> <Output_sorted_bam>\n",argv[0]);
         printf("      	 	-nodes    60      - Number of CPUs requested\n");
         printf("      		-data     ccs     - PacBio Hifi\n");
         printf("		-data     ont     - Oxford Nanopore Q20 or Q30\n");
         printf("		-data     ont-NLR - Oxford Nanopore normal long reads (NLR)\n");
         printf("		-data     ngs-HiC - Hi-C reads\n");
         printf("      		-data     ngs-10X - 10X reads\n");
         printf("		-data     ngs-SSR - Standard short reads such as Illumina data\n");
         exit(1);
       }
       else if(!strcmp(argv[i],"-debug"))
       {
         sscanf(argv[++i],"%d",&n_debug);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-file"))
       {
         sscanf(argv[++i],"%d",&file_tag); 
         args=args+2;
       }
    }

    if((argc == 2)&&(flag_align == 1))
    {
         printf("===Align reads to a reference for all data types:\n");
         printf("===Output a coordinate sorted bam file:\n");
         printf("	Usage: %s align -nodes 60 -data ccs -reads input_long.fasta/q <Input_Reference> <Output_sorted_bam>\n",argv[0]);
         printf("	Usage: %s align -nodes 60 -data ont -reads input_long.fasta/q <Input_Reference> <Output_sorted_bam>\n",argv[0]);
         printf("	Usage: %s align -nodes 60 -data ont-NLR -reads input_long.fasta/q <Input_Reference> <Output_sorted_bam>\n",argv[0]);
         printf("	Usage: %s align -nodes 60 -data ngs-HiC -fq1 read_1.fq.gz -fq2 read_2.fq.gz <Input_Reference> <Output_sorted_bam>\n",argv[0]);
         printf("	Usage: %s align -nodes 60 -data ngs-10X -fq1 read_1.fq.gz -fq2 read_2.fq.gz <Input_Reference> <Output_sorted_bam>\n",argv[0]);
         printf("	Usage: %s align -nodes 60 -data ngs-SSR -fq1 read_1.fq.gz -fq2 read_2.fq.gz <Input_Reference> <Output_sorted_bam>\n",argv[0]);
         printf("      	 	-nodes    60      - Number of CPUs requested\n");
         printf("      		-data     ccs     - PacBio Hifi\n");
         printf("		-data     ont     - Oxford Nanopore Q20 or Q30\n");
         printf("		-data     ont-NLR - Oxford Nanopore normal long reads (NLR)\n");
         printf("		-data     ngs-HiC - Hi-C reads\n");
         printf("      		-data     ngs-10X - 10X reads\n");
         printf("		-data     ngs-SSR - Standard short reads such as Illumina data\n");
         exit(1);
    }
    if((argc == 2)&&(flag_breakpoint == 1))
    {
         printf("===Detect breakpoints with aligned, and name sorted sam,bam or cram files:\n");
         printf("	Usage: %s breakpoint -data ccs -bam /myspace/desk/test-sorted.bam <output_breakpoints.vcf>\n",argv[0]);
         printf("	Usage: %s breakpoint -data ont -bam /myspace/desk/test-sorted.bam <output_breakpoints.vcf>\n",argv[0]);
         printf("	Usage: %s breakpoint -data ngs-HiC -bam /myspace/desk/test-sorted.bam <output_breakpoints.vcf>\n",argv[0]);
         printf("	Usage: %s breakpoint -data ngs-10X -bam /myspace/desk/test-sorted.bam <output_breakpoints.vcf>\n",argv[0]);
         printf("	Usage: %s breakpoint -data ngs-SSR -bam /myspace/desk/test-sorted.bam <output_breakpoints.vcf>\n",argv[0]);
         printf("      		-data     ccs     - PacBio Hifi\n");
         printf("		-data     ont     - Oxford Nanopore Q20 or Q30\n");
         printf("		-data     ngs-HiC - Hi-C reads\n");
         printf("      		-data     ngs-10X - 10X reads\n");
         printf("		-data     ngs-SSR - Standard short reads such as Illumina data\n");
         printf("	--- Note: the sam/bam/cram file has to be Name sorted! ---\n");
         printf("	--- If not read name sorted, please do the following:  ---\n");
         printf("	--- samtools sort -@ 60 -n your.bam new.bam ---\n\n");

         printf("===Detect breakpoints with fasta/fastq long read files:\n");
         printf("	Usage: %s breakpoint -nodes 60 -data ccs/ont -reads input_long.fasta/q <Input_Reference> <breakpoints.vcf>\n",argv[0]);
         printf("\n");
         exit(1);
    }
    if((argc == 2)&&(flag_plot == 1))
    {
         printf("===Plot depth of coverage for all data types:\n");
         printf("===Input a coordinate sorted bam file\n");
         printf("===Output a tmp directory containing coverage images for 23 chromosomes chr{1,22,X}\n");
         printf("	Usage: %s plot -bam /myspace/desk/test-sorted.bam -sample cancer-XXX -hight 180 -window 100 -denoise 1\n",argv[0]);
         printf("      		  -sample:   Sample name\n");
         printf("      		  -hight:    Maximum value in Y axis (read depth)\n");
         printf("      		  -window:   Window size to display chromosome coordinates\n");
         printf("      		  -chrnum:   Number of chromosomes {23}\n");
         printf("      		  -denoise:  Noise reduction option {0,1,2,3}\n");
         printf("      		     (0) No noise reduction option\n");
         printf("      		     (1) Average the data points after filtering high and low points within the window size\n");
         printf("      		     (2) Fit the datasets into copy numbers after filtering/smoothing\n");
         printf("      		     (3) Dimensionless copy numbers in logscale after filtering/smoothing\n");
         printf("\n");
         printf("===Without noise reduction:\n");
         printf("	Usage: %s plot -bam /myspace/desk/test-sorted.bam -sample cancer -denoise 0\n",argv[0]);
         printf("\n");
         exit(1);
    }
    pid = getpid();
    memset(tempa,'\0',2000);
    if (!getcwd(tempa, sizeof(tempa)))
    {
      exit(1);
    } 
    memset(tmpdir,'\0',2000);
    memset(workdir,'\0',2000);
    sprintf(tmpdir,"%s/",tempa);
    sprintf(workdir,"%s/tmp_rununik_%d/",tempa,pid);

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"mkdir %s",workdir);

    RunSystemCommand(syscmd);

    if(chdir(workdir) == -1)
    {
      printf("System command error: chdir\n");
      exit(EXIT_FAILURE);
    }
     
    st = argv[0];
    ed = strrchr(argv[0],'/');
    memset(tempc,'\0',2000);
    strncpy(tempc,argv[0],ed-st);
    memset(bindir,'\0',2000);
    sprintf(bindir,"%s/step-bin",tempc);

    if((strcmp(datatype,"ccs") == 0))
      data_flag = 1;
    else if((strcmp(datatype,"ont") == 0))
      data_flag = 2;
    else if((strcmp(datatype,"ngs-HiC") == 0))
      data_flag = 3;
    else if((strcmp(datatype,"ngs-10X") == 0))
      data_flag = 4;
    else if((strcmp(datatype,"ngs-SSR") == 0))
      data_flag = 5;
    else if((strcmp(datatype,"ont-NLR") == 0))
    {
      data_flag = 2;
      NLR_flag = 1;
    }
    else 
    {
      printf("Define the datatype: ccs, ont, ngs-10X etc!\n");
      exit(1);
    }

    if(flag_align == 1)
    {
      memset(file_bambam,'\0',2000);
      memset(file_tarseq,'\0',2000);
      memset(file_lgread,'\0',2000);
      memset(file_Stones,'\0',2000);
      memset(file_read1,'\0',2000);
      memset(file_read2,'\0',2000);
      sprintf(file_read1,"%s/%s",tempa,fq1name);
      sprintf(file_read2,"%s/%s",tempa,fq2name);
      sprintf(file_bambam,"%s/%s",tempa,datname);
      if(bam_flag == 1)
        sprintf(file_Stones,"%s/%s",tempa,argv[args]);
      else
      {
        sprintf(file_tarseq,"%s/%s",tempa,argv[args]);
        sprintf(file_Stones,"%s/%s",tempa,argv[args+1]);
      }
      if((data_flag == 4)||(data_flag == 5))
      { 
        memset(syscmd,'\0',2000);
        printf("Comms: %s\n",syscmd);
        printf("Ref file: %s\n",file_tarseq);
        sprintf(syscmd,"%s/bwa-mem2 index %s > try.out",bindir,file_tarseq);
        RunSystemCommand(syscmd);

        memset(syscmd,'\0',2000);
        sprintf(syscmd,"%s/bwa-mem2 mem -t %d %s %s %s> align.sam",bindir,n_nodes,file_tarseq,file_read1,file_read2);
        RunSystemCommand(syscmd);
      }
      else if(data_flag == 3)
      {
        memset(syscmd,'\0',2000);
        sprintf(syscmd,"%s/bwa-mem2 index %s > try.out",bindir,file_tarseq);
        RunSystemCommand(syscmd);

        memset(syscmd,'\0',2000);
        sprintf(syscmd,"%s/bwa-mem2 mem -t %d -5SPM %s %s %s > align.sam",bindir,n_nodes,file_tarseq,file_read1,file_read2);
        RunSystemCommand(syscmd);
      }
      else if(data_flag == 1)
      { 
        sprintf(file_lgread,"%s/%s",tempa,fastname);
        memset(syscmd,'\0',2000);
        sprintf(syscmd,"%s/minimap2 -a -x asm10 -t %d %s %s > align.sam",bindir,n_nodes,file_tarseq,file_lgread);
	RunSystemCommand(syscmd);
      }
      else if(data_flag == 2)
      {
        sprintf(file_lgread,"%s/%s",tempa,fastname);
        ed = strrchr(file_lgread,'.');
        memset(tagname,'\0',10);
        if(ed != NULL)
        {
          strcpy(tagname,ed);
          if((strncmp(tagname,".gz",3))==0)
          {
            memset(syscmd,'\0',2000);
            sprintf(syscmd,"%s/pigz -dc %s > reads_ont.fq ",bindir,file_lgread);
            RunSystemCommand(syscmd);

            memset(syscmd,'\0',2000);
            sprintf(syscmd,"%s/step_shred -rlength 8000 reads_ont.fq reads_ont.shred > try.out",bindir);
            RunSystemCommand(syscmd);

            memset(syscmd,'\0',2000);
            sprintf(syscmd,"rm -rf reads_ont.fq");
            RunSystemCommand(syscmd);
          }
          else
          {
            memset(syscmd,'\0',2000);
            sprintf(syscmd,"%s/step_shred -rlength 8000 %s reads_ont.shred ",bindir,file_lgread);
            RunSystemCommand(syscmd);
          }
        }
        else
        {
            memset(syscmd,'\0',2000);
            sprintf(syscmd,"%s/step_shred -rlength 8000 %s reads_ont.shred ",bindir,file_lgread);
            RunSystemCommand(syscmd);
        }	      
        memset(syscmd,'\0',2000);
	if(NLR_flag == 0)
          sprintf(syscmd,"%s/minimap2 -a -x asm10 -t %d %s reads_ont.shred > align.sam",bindir,n_nodes,file_tarseq);
	else
          sprintf(syscmd,"%s/minimap2 -a -x map-ont -t %d %s reads_ont.shred > align.sam",bindir,n_nodes,file_tarseq);
	RunSystemCommand(syscmd);
      }
      memset(syscmd,'\0',2000);
      sprintf(syscmd,"%s/samtools view -@ %d -bS align.sam > Sorted_names.bam",bindir,n_nodes);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"rm -rf align.sam");
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"%s/samtools fixmate -@ %d -m Sorted_names.bam Fixmate.bam > try.out",bindir,n_nodes);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"rm -rf Sorted_names.bam");
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"%s/samtools sort -@ %d -o Sorted.bam Fixmate.bam > try.out",bindir,n_nodes);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"rm -rf Fixmate.bam");
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"%s/samtools markdup -@ %d -r -s Sorted.bam Dupmarked.bam > try.out",bindir,n_nodes);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"rm -rf Sorted.bam");
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"mv Dupmarked.bam %s",file_Stones);
      RunSystemCommand(syscmd);
      return EXIT_SUCCESS;
    }


    memset(file_tarseq,'\0',2000);
    memset(file_lgread,'\0',2000);
    memset(file_Stones,'\0',2000);
    memset(file_snpfrq,'\0',2000);
    memset(file_ctspec,'\0',2000);
    memset(file_bambam,'\0',2000);
    memset(file_datas,'\0',2000);

    if(flag_plot == 1)
    {
      if(bam_flag == 1)
      {
        sprintf(file_bambam,"%s/%s",tempa,datname);
        printf("Bam file: %s %s\n",datname,file_bambam);
        if((namef = fopen(datname,"r")) == NULL)
        {
          printf("File not in the working directory!\n");
          if((namef = fopen(file_bambam,"r")) == NULL)
          {
            printf("BAM File %s not found and please copy it to your working directory!\n",datname);
            exit(1);
          }
          else
          {
            printf("Input bam file: %s\n",file_bambam);
	    memset(window_word,'\0',100);
	    sprintf(window_word,"($%s%d==0)","2%",window_size);
            memset(syscmd,'\0',2000);
            sprintf(syscmd,"%s/samtools depth -@ %d %s | awk '%s%s' > depth-t2t.dat",bindir,n_nodes,file_bambam,window_word,"{print $1,$2,$3}");
            printf("Command: %d %s\n",window_size,syscmd);
            RunSystemCommand(syscmd);
          }
        }
        else
        {
          printf("BAM name: %s\n",datname);
	  memset(window_word,'\0',100);
	  sprintf(window_word,"($%s%d==0)","2%",window_size);
          memset(syscmd,'\0',2000);
          sprintf(syscmd,"%s/samtools depth -@ %d %s | awk '%s%s' > depth-t2t.dat",bindir,n_nodes,datname,window_word,"{print $1,$2,$3}");
            printf("Command: %d %s\n",window_size,syscmd);
          RunSystemCommand(syscmd);
        }

/*      Plot coverage profiles chr by chr */
  	memset(commands,'\0',2000);
        sprintf(commands,"%s/step_depthPlot",bindir);

        memset(syscmd,'\0',2000);
        sprintf(syscmd,"%s/step_commsPlot -command %s -sample %s -denoise %d -hight %d -chrnum %d sh.plots",bindir,commands,sample,denoise_flag,y_hight,num_chr);
        RunSystemCommand(syscmd);

        memset(syscmd,'\0',2000);
        sprintf(syscmd,"sort -k 3,3n depth-t2t.dat | awk '%s' > coverage.dat","{print $2,$3}");
        RunSystemCommand(syscmd);

        memset(syscmd,'\0',2000);
        sprintf(syscmd,"%s/step_coverage coverage.dat coverage.freq",bindir);
        RunSystemCommand(syscmd);

/*      Plot coverage frequency chr by chr */
  	memset(commands,'\0',2000);
        sprintf(commands,"%s/step_freqPlot",bindir);

        memset(syscmd,'\0',2000);
        sprintf(syscmd,"%s/step_commsPlot -command %s -sample %s -chrnum %d sh.plot-freq",bindir,commands,sample,num_chr);
        RunSystemCommand(syscmd);

        memset(syscmd,'\0',2000);
        sprintf(syscmd,"bash sh.plots");
        RunSystemCommand(syscmd);

        memset(syscmd,'\0',2000);
        sprintf(syscmd,"bash sh.plot-freq");
        RunSystemCommand(syscmd);
      }
      else
      {
        printf("BAM File was not found!\n");
        exit(1);
      }

      return EXIT_SUCCESS;
    }


    if(bam_flag != 1)
    {
      sprintf(file_lgread,"%s/%s",tempa,fastname);
      sprintf(file_tarseq,"%s/%s",tempa,argv[args]);
      sprintf(file_Stones,"%s/%s",tempa,argv[args+1]);
//    sprintf(file_snpfrq,"%s/%s.snps",tempa,argv[args+1]);
//    sprintf(file_ctspec,"%s/%s.spec",tempa,argv[args+1]);

      if((namef = fopen(file_tarseq,"r")) == NULL)
      {
        printf("File not in the working directory!\n");
        if((namef = fopen(argv[args],"r")) == NULL)
        {
          printf("File %s not found and please copy it to your working directory!\n",argv[args]);
          exit(1);
        }
        else
        {
          memset(file_tarseq,'\0',2000);
          strcpy(file_tarseq,argv[args]);
          printf("Input target assembly file1: %s\n",file_tarseq);
        }
      }
      else
      {
        printf("Input target assembly file2: %s\n",file_tarseq);
      } 

      if((namef = fopen(file_lgread,"r")) == NULL)
      {
        printf("File not in the working directory!\n");
        if((namef = fopen(fastname,"r")) == NULL)
        {
          printf("File %s not found and please copy it to your working directory!\n",argv[args+1]);
          exit(1);
        }
        else
        {
          memset(file_lgread,'\0',2000);
          strcpy(file_lgread,fastname);
          printf("Input reference file: %s\n",file_lgread);
        }
      }
      else
      {
        printf("Input fasta file: %s\n",file_lgread);
      }

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"%s/step_fastq -name tarseq -len 20000 %s tarseq.fastq tarseq.tag > try.out",bindir,file_tarseq);
      RunSystemCommand(syscmd);
    
      if(data_flag == 2)
      {
	ed = strrchr(file_lgread,'.');
	memset(tagname,'\0',10);
	if(ed != NULL)
        {
          strcpy(tagname,ed);
          if((strncmp(tagname,".gz",3))==0)
	  {
            memset(syscmd,'\0',2000);
            sprintf(syscmd,"%s/pigz -dc %s > reads_ont.fq ",bindir,file_lgread);
            RunSystemCommand(syscmd);

            memset(syscmd,'\0',2000);
            sprintf(syscmd,"%s/step_shred -rlength 8000 reads_ont.fq reads_ont.shred > try.out",bindir);
            RunSystemCommand(syscmd);

            memset(syscmd,'\0',2000);
            sprintf(syscmd,"rm -rf reads_ont.fq");
            RunSystemCommand(syscmd);
	  }
          else
	  {
            memset(syscmd,'\0',2000);
            sprintf(syscmd,"%s/step_shred -rlength 8000 %s reads_ont.shred ",bindir,file_lgread);
            RunSystemCommand(syscmd);
	  }
        }
        else
        {
            memset(syscmd,'\0',2000);
            sprintf(syscmd,"%s/step_shred -rlength 8000 %s reads_ont.shred ",bindir,file_lgread);
            RunSystemCommand(syscmd);
        }
      }
    }
    else
    {
      memset(file_Stones,'\0',2000);
      sprintf(file_Stones,"%s/%s",tempa,argv[args]);
      printf("VCF file: %s %s\n",file_Stones,argv[args]);
    }
   
    if(run_align)
    {
      if((strcmp(toolname,"bwa") == 0))
      { 
        memset(syscmd,'\0',2000);
        sprintf(syscmd,"%s/bwa index tarseq.fastq > try.out",bindir);
        RunSystemCommand(syscmd);

        memset(syscmd,'\0',2000);
        sprintf(syscmd,"%s/bwa index tarseq.fastq > try.out",bindir);
        RunSystemCommand(syscmd);

        memset(syscmd,'\0',2000);
        sprintf(syscmd,"%s/bwa mem -t %d tarseq.fastq %s | egrep tarseq_ | awk '%s' | egrep -v '^@' > align.dat0",bindir,n_nodes,file_lgread,"{print $1,$2,$3,$4,$5,$6,$10}");
        RunSystemCommand(syscmd);
      }
      else if((strcmp(toolname,"minimap2") == 0))
      { 
        memset(syscmd,'\0',2000);
	if(data_flag == 2)
	{
          sprintf(syscmd,"%s/minimap2 -a -x asm10 -t %d tarseq.fastq reads_ont.shred | egrep tarseq_ | awk '%s' | egrep -v '^@' > align.dat0",bindir,n_nodes,"($5>0){print $1,$2,$3,$4,$5,$6,$10}");
	}
	else
	{
	  if(platform_tag == 1)
            sprintf(syscmd,"%s/minimap2 -a -x asm10 -t %d tarseq.fastq %s | egrep tarseq_ | awk '%s' | egrep -v '^@' > align.dat0",bindir,n_nodes,file_lgread,"($5>0){print $1,$2,$3,$4,$5,$6,$10}");
	  else if(platform_tag == 2)
            sprintf(syscmd,"%s/minimap2 -a -x map-ont -t %d tarseq.fastq %s | egrep tarseq_ | awk '%s' | egrep -v '^@' > align.dat0",bindir,n_nodes,file_lgread,"($5>0){print $1,$2,$3,$4,$5,$6,$10}");
	}
        RunSystemCommand(syscmd);
      }
      else if((strcmp(toolname,"minimap3") == 0))
      { 
        memset(syscmd,'\0',2000);
        sprintf(syscmd,"%s/minimap3 -a -x asm20 -t %d tarseq.fastq %s | egrep tarseq_ | awk '%s' | egrep -v '^@' > align.dat0",bindir,n_nodes,file_lgread,"{print $1,$2,$3,$4,$5,$6,$10}");
        RunSystemCommand(syscmd);
      }
      else if((strcmp(toolname,"smalt") == 0))
      {
        memset(syscmd,'\0',2000);
        sprintf(syscmd,"%s/smalt index -k 13 -k 13 hash_genome tarseq.fastq > try.out",bindir);
        RunSystemCommand(syscmd);

        memset(syscmd,'\0',2000);
        sprintf(syscmd,"%s/smalt map -n %d -m 100 -f samsoft -O hash_genome %s | egrep tarseq_ | awk '%s' | egrep -v '^@' > align.dat0",bindir,n_nodes,file_lgread,"($2<100){print $1,$2,$3,$4,$5,$6,$10}");
        RunSystemCommand(syscmd);
      }
      else
      {
        printf("Give an alignment tool! \n");
        exit(1);
      }
    }
    else
    {
      if(bam_flag == 1)
      {
        sprintf(file_bambam,"%s/%s",tempa,datname);
        printf("Bam file: %s %s\n",datname,file_bambam);
        if((namef = fopen(datname,"r")) == NULL)
        {
          printf("File not in the working directory!\n");
          if((namef = fopen(file_bambam,"r")) == NULL)
          {
            printf("BAM File %s not found and please copy it to your working directory!\n",datname);
            exit(1);
          }
          else
          {
            printf("Input bam file: %s\n",file_bambam);
            memset(syscmd,'\0',2000);
            sprintf(syscmd,"%s/samtools view -@ %d %s | awk '%s' > align.dat0",bindir,n_nodes,file_bambam,"(length($6)>5)&&($5>0){print $1,$2,$3,$4,$5,$6,$10}");
            RunSystemCommand(syscmd);
          }
        }
        else
        {
          printf("Input target assembly file2: %s\n",file_bambam);
          memset(syscmd,'\0',2000);
          sprintf(syscmd,"%s/samtools view -@ %d %s | awk '%s' > align.dat0",bindir,n_nodes,datname,"(length($6)>5)&&($5>0){print $1,$2,$3,$4,$5,$6,$10}");
          RunSystemCommand(syscmd);
        }
      }
      else
      {
        memset(syscmd,'\0',2000);
        sprintf(syscmd,"cp %s align.dat0",datname);
        RunSystemCommand(syscmd);
      } 
    } 

    memset(syscmd,'\0',2000);
    printf("%s/step_number align.dat0 align.dat > try.out",bindir);
    if(data_flag != 3)
      sprintf(syscmd,"%s/step_number align.dat0 align.dat > try.out",bindir);
    else
      sprintf(syscmd,"%s/step_short-number align.dat0 align.dat > try.out",bindir);
    RunSystemCommand(syscmd);

    memset(syscmd,'\0',2000);
    printf("%s/stepBreakPoint align.dat break.dat > try.out",bindir);
    if(data_flag != 3)
      sprintf(syscmd,"%s/stepBreakPoint align.dat break.dat > try.out",bindir);
    else
      sprintf(syscmd,"%s/stepBreakPoint-S align.dat break.dat > try.out",bindir);
    RunSystemCommand(syscmd);
    
    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/step_cleanProcess break.dat break.clean > try.out",bindir);
    RunSystemCommand(syscmd);

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"sort -k 3,3n -k 5,5n -k 4,4 break.clean > break.sort");
    RunSystemCommand(syscmd);

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/step_breakProcess break.sort break-point.dat > try.out",bindir);
    RunSystemCommand(syscmd);

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"sort -k 3,3n -k 4,4 -k 5,5n break-point.dat > break-point.sort");
    RunSystemCommand(syscmd);

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/step_breakProcess break-point.sort stones-point.dat > try.out",bindir);
    RunSystemCommand(syscmd);

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"sort -k 3,3n -k 5,5n -k 4,4 stones-point.dat > stones-point.sort");
    RunSystemCommand(syscmd);

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/step_breakSort stones-point.sort stones.dat > try.out",bindir);
    RunSystemCommand(syscmd);

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"sort -k 3,3 -k 4,4n -k 5,5n stones.dat > stones.sort");
    RunSystemCommand(syscmd);

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/step_linkStones stones.sort stones.link > try.out",bindir);
    RunSystemCommand(syscmd);

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"cat stones.link | awk '{print $1,$4,$5,$2,$3,$6,$7,$8,$9,$11,$10,$12,$11,$15,$14}' > stones.link2");
    RunSystemCommand(syscmd);

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"cat stones.link stones.link2 | sort -k 2,2n -k 4,4n -k 3,3n -k 5,5n > stones-all.link");
    RunSystemCommand(syscmd);

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/step_edgeStones stones-all.link stones.edge > try.out",bindir);
    RunSystemCommand(syscmd);

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"cp stones.edge %s",file_Stones);
    RunSystemCommand(syscmd);

    if(n_debug == 0)
    {
      memset(syscmd,'\0',2000);
      sprintf(syscmd,"rm -rf * > /dev/null");
    
      RunSystemCommand(syscmd);
      if(chdir(tmpdir) == -1)
      {
        printf("System command error: chdir\n");
        exit(EXIT_FAILURE);
      }

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"rm -rf %s > /dev/null",workdir);
    
      RunSystemCommand(syscmd);
    }
    return EXIT_SUCCESS;

}
/* end of the main */



#define SWAP(a,b) temp=(a);(a)=b;(b)=temp;

/*   Subroutine to sort an array arr[0,...,n-1] into ascending order while
     making the corresponding reaarangement of the array brr[0,...,n-1]
     by the use of Quicksort (Sedgwick, R. 1978, Communications o fthe ACM,
     vol. 21, pp. 847-857) also see Numerical Recipes in C                  */  

/* =============================== */
void ArraySort_Long(int n, B64_long *arr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,NSTACK=50,istack[NSTACK];
     B64_long a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
             }
             arr[i+1]=a;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}


/* =============================== */
void ArraySort_Int(int n, int *arr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,NSTACK=50,istack[NSTACK];
     int a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
             }
             arr[i+1]=a;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}


/* =============================== */
void ArraySort_Mix(int n, B64_long *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     B64_long a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}

/* =============================== */
void ArraySort_Int2(int n, int *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     int a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}

/*   function to sort an array into a decreasing order:  a>b>c>....    */  
/* =============================== */
void ArraySort2_Int2(int n, int *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     int a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]>=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);

          if(arr[m]<arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
          }

          if(arr[m+1]<arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
          }

          if(arr[m]<arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          for(;;)
          {
             do i++; while (arr[i]>a);
             do j--; while (arr[j]<a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}

/* =============================== */
void ArraySort_Mix3(int n, B64_long *arr, int *brr, int *crr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,c,NSTACK=50,istack[NSTACK];
     B64_long a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             c=crr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
                crr[i+1]=crr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
             crr[i+1]=c;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);
          SWAP(crr[k],crr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
            SWAP(crr[m],crr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
            SWAP(crr[m+1],crr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
            SWAP(crr[m],crr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          c=crr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
             SWAP(crr[i],crr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          crr[m+1]=crr[j];
          crr[j]=c;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}


/*   to swap the string arrays           */
/* ============================================= */
void s_swap(char **Pair_Name, int i, int j)
/* ============================================= */
{
     char temp[Max_N_NameBase];

     strcpy(temp,Pair_Name[i]);
     strcpy(Pair_Name[i],Pair_Name[j]);
     strcpy(Pair_Name[j],temp);
}


/*   to sort the string array in order          */
/* ============================================= */
void ArraySort_String(int n, char **Pair_Name, int *brr)
/* ============================================= */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     int temp,MIN=7;
     char p[Max_N_NameBase];

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             strcpy(p,Pair_Name[j]);
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(strcmp(Pair_Name[i],p)<=0) break;
                strcpy(Pair_Name[i+1],Pair_Name[i]);
                brr[i+1]=brr[i];
             }
             strcpy(Pair_Name[i+1],p);
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          s_swap(Pair_Name,k,m+1);
          SWAP(brr[k],brr[m+1]);

          if(strcmp(Pair_Name[m],Pair_Name[ir])>0)
          {
            s_swap(Pair_Name,m,ir);
            SWAP(brr[m],brr[ir]);
          }

          if(strcmp(Pair_Name[m+1],Pair_Name[ir])>0)
          {
            s_swap(Pair_Name,m+1,ir);
            SWAP(brr[m+1],brr[ir]);
          }

          if(strcmp(Pair_Name[m],Pair_Name[m+1])>0)
          {
            s_swap(Pair_Name,m,m+1);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          strcpy(p,Pair_Name[m+1]);
          b=brr[m+1];
          for(;;)
          {
             do i++; while (strcmp(Pair_Name[i],p)<0);
             do j--; while (strcmp(Pair_Name[j],p)>0);
             if(j<i) break;
             s_swap(Pair_Name,i,j);
             SWAP(brr[i],brr[j]);
          }
          strcpy(Pair_Name[m+1],Pair_Name[j]);
          strcpy(Pair_Name[j],p);
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}

/*   to swap the string arrays           */
/* ============================================= */
void s_swap2(char **Pair_Name, int i, int j)
/* ============================================= */
{
     char temp[Max_N_NameBase2];

     strcpy(temp,Pair_Name[i]);
     strcpy(Pair_Name[i],Pair_Name[j]);
     strcpy(Pair_Name[j],temp);
}


/*   to sort the string array in order          */
/* ============================================= */
void ArraySort_String2(int n, char **Pair_Name, int *brr)
/* ============================================= */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     int temp,MIN=7;
     char p[Max_N_NameBase2];

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             strcpy(p,Pair_Name[j]);
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(strcmp(Pair_Name[i],p)<=0) break;
                strcpy(Pair_Name[i+1],Pair_Name[i]);
                brr[i+1]=brr[i];
             }
             strcpy(Pair_Name[i+1],p);
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          s_swap2(Pair_Name,k,m+1);
          SWAP(brr[k],brr[m+1]);

          if(strcmp(Pair_Name[m],Pair_Name[ir])>0)
          {
            s_swap2(Pair_Name,m,ir);
            SWAP(brr[m],brr[ir]);
          }

          if(strcmp(Pair_Name[m+1],Pair_Name[ir])>0)
          {
            s_swap2(Pair_Name,m+1,ir);
            SWAP(brr[m+1],brr[ir]);
          }

          if(strcmp(Pair_Name[m],Pair_Name[m+1])>0)
          {
            s_swap2(Pair_Name,m,m+1);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          strcpy(p,Pair_Name[m+1]);
          b=brr[m+1];
          for(;;)
          {
             do i++; while (strcmp(Pair_Name[i],p)<0);
             do j--; while (strcmp(Pair_Name[j],p)>0);
             if(j<i) break;
             s_swap2(Pair_Name,i,j);
             SWAP(brr[i],brr[j]);
          }
          strcpy(Pair_Name[m+1],Pair_Name[j]);
          strcpy(Pair_Name[j],p);
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}


/* creat an int matrix with subscript ange m[nrl...nrh][ncl...nch]  */
int     **imatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch)
{
        B64_long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        int  **m;

        /* allocate pointers to rows        */
        if((m=(int **)calloc(nrow,sizeof(int*)))==NULL)
        {
           printf("error imatrix: calloc error No. 1 \n");
           return(NULL);
        }
        m+=0;
        m-=nrl;

        /* allocate rows and set pointers to them        */
        if((m[nrl]=(int *)calloc(nrow*ncol,sizeof(int)))==NULL)
        {
           printf("error imatrix: calloc error No. 2 \n");
           return(NULL);
        }
        m[nrl]+=0;
        m[nrl]-=nrl;

        for(i=nrl+1;i<=nrh;i++)
           m[i]=m[i-1]+ncol;
        /* return pointer to array of pointers to rows   */
        return m;
}

/* creat char matrix with subscript ange cm[nrl...nrh][ncl...nch]  */
char    **cmatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch)
{
        B64_long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        char **cm;

        /* allocate pointers to rows        */
        if((cm=(char **)calloc(nrow,sizeof(char*)))==NULL)
        {
           printf("error cmatrix: calloc error No. 1 \n");
           return(NULL);
        }
        cm+=0;
        cm-=nrl;

        /* allocate rows and set pointers to them        */
        if((cm[nrl]=(char *)calloc(nrow*ncol,sizeof(char)))==NULL)
        {
           printf("error cmatrix: calloc error No. 2 \n");
           return(NULL);
        }
        cm[nrl]+=0;
        cm[nrl]-=nrl;

        for(i=nrl+1;i<=nrh;i++)
           cm[i]=cm[i-1]+ncol;
        /* return pointer to array of pointers to rows   */
        return cm;
}
