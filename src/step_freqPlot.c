/****************************************************************************
 ****************************************************************************
 *                                                                          *
 *  Copyright (C) 2018  Genome Research Ltd.                                *
 *                                                                          *
 *  Author: Zemin Ning (zn1@sanger.ac.uk)                                   *
 *                                                                          *
 *  This file is part of stepStone pipeline.                                *
 *                                                                          *
 *  ScaffHiC is a free software: you can redistribute it and/or modify it   *
 *  under the terms of the GNU General Public License as published by the   *
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

#define GT '>'
#define GT4 (((((GT<<8)+GT)<<8)+GT)<<8)+GT

#define ENDS_EXTRA 0
#define PADCHAR '-'
#define MAX_N_BRG 50000 
#define MAX_N_ROW 50000 
#define Max_N_NameBase 60 
#define Max_N_Pair 100
static char **S_Name,**R_Name,**R_Name2,**T_Name,**cellname;
static int *hit_locus,*hit_index,*hit_cover,*hit_tagcode,*ctg_index,*ctg_cover,*ctg_offset;

/* SSAS default parameters   */
static int IMOD=0;
static int n_type=0;
static int barreads=5;
static int file_flag=2;
static int tiles_flag=0;
static int denoise_flag=3;
static int y_hight=180;
static int edge_flag=0;
static int nContig=0;
static int G_Size = 0;
static int n_lenn = 14;
static int max_len = 100000;
typedef struct
{
       int foffset;
       int fsindex;
} SIO;


static char rc_char[500000];
static char chromo[50];
static char sample[50];


int main(int argc, char **argv)
{
    FILE *namef;
    int i,j,nSeq,args;
    int n_contig,n_reads,n_readsMaxctg,nseq;
    void decodeReadpair(int nSeq);
    void HashFasta_Head(int i, int nSeq);
    void HashFasta_Table(int i, int nSeq);
    void Assemble_SM(int arr,int brr);
    void Mapping_Process(char **argv,int args,int nSeq);
    void Memory_Allocate(int arr);
    char line[2000]={0},temp[60],cc[60],RC[5],chrname[200],*st,*ed;
    char **cmatrix(long nrl,long nrh,long ncl,long nch);

    if(argc < 2)
    {
      printf("Usage: %s [-chr chr1] [-sample OES103] <SAM_depth_file> <sh.coverplot>\n",argv[0]);

      exit(1);
    }

    nSeq=0;
    args=1;
    for(i=1;i<argc;i++)
    {
       if(!strcmp(argv[i],"-mod"))
       {
         sscanf(argv[++i],"%d",&IMOD); 
         args=args+2;
       }
       else if(!strcmp(argv[i],"-chr"))
       {
         sscanf(argv[++i],"%s",chromo); 
         args=args+2;
       }
       else if(!strcmp(argv[i],"-sample"))
       {
         sscanf(argv[++i],"%s",sample); 
         args=args+2;
       }
       else if(!strcmp(argv[i],"-tile"))
       {
         sscanf(argv[++i],"%d",&tiles_flag);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-reads"))
       {
         sscanf(argv[++i],"%d",&barreads);
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
       else if(!strcmp(argv[i],"-max"))
       {
         sscanf(argv[++i],"%d",&max_len);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-file"))
       {
         sscanf(argv[++i],"%d",&file_flag);
         args=args+2;
       }
    }

    nseq=0;
    if((namef = fopen(argv[args],"r")) == NULL)
    {
      printf("ERROR main:: args \n");
      exit(1);
    }
    while(!feof(namef))
    {
      if(fgets(line,2000,namef) == NULL)
      {
//        printf("fgets command error:\n);
      }
      if(feof(namef)) break;
      nseq++;
    }
    fclose(namef); 
   
    if((ctg_cover = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - hit_locus2\n");
      exit(1);
    }
    if((ctg_offset = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - ctg_left\n");
      exit(1);
    }
    if((hit_index = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - ctg_left\n");
      exit(1);
    }
    if((hit_cover = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - ctg_left\n");
      exit(1);
    }

    S_Name=cmatrix(0,nseq+10,0,8);
    if((namef = fopen(argv[args],"r")) == NULL)
    {
      printf("ERROR main:: args \n");
      exit(1);
    }
/*  read the alignment files         */
    i=0;
    while(fscanf(namef,"%s %d %d",chrname,&ctg_offset[i],&ctg_cover[i])!=EOF)
    {
        if((strncmp(chrname,"chr",3))==0)
           strcpy(S_Name[i],chrname);
        else
        {
           sprintf(S_Name[i],"%s%s","chr",chrname);
        }
        i++;
    }
    fclose(namef);
 
    printf("Plot frequency: %s %s\n",chromo,sample);
    nContig = nseq;

    n_reads=i;
//    Readname_match(seq,argv,args,n_reads,nRead);
    Mapping_Process(argv,args,n_reads);
//    Read_Pairs(argv,args,seq,n_reads);

    return EXIT_SUCCESS;

}
/* end of the main */

/*   subroutine to sort out read pairs    */
/* =============================== */
void Mapping_Process(char **argv,int args,int nSeq)
/* =============================== */
{
     int i,j,k,m,n,n_chr,n_oando;
     int num_hits,num_pairs,num_blocks,num_trans,num_sends,num_reads;
     int stopflag,*chr_trans,*chr_rdhit,*chr_inchr;
     int offset,offset1,offset2,block_len,block_ttt;
     float rate1,rate2,rate3;
     FILE *namef;
     void Frequency_Process(char *plotname,int nSeq);
     char outputname[100];
     void ArraySort_Int2(int n, int *arr, int *brr);
     char **cmatrix(long nrl,long nrh,long ncl,long nch);
     void ArraySort_String(int n,char **Pair_Name,int *brr);
     char KKK0[100],KKK00[100],KKK1[100],KKK2[100],KKK22[100],KKK3[100],KKK4[100],KKK5[100],KKK6[100],KKK7[100];
     char *st,*ed,Sam_name[30],syscmd[2000],Chr_name[30];

     n_oando = 1;
     strcpy(KKK00,"set logscale x");     
     strcpy(KKK0,"set logscale y");     
     strcpy(KKK1,"set terminal svg");     
     strcpy(KKK2,"set style line 1 lt 1 lw 3 pt 3 linecolor rgb \\\"red\\\"");     
     strcpy(KKK22,"set style line 2 lt 1 lw 3 pt 3 linecolor rgb \\\"blue\\\"");     
     strcpy(KKK4,"set xlabel \\\"Coverage depth\\\"");     
     strcpy(KKK3,"set ylabel \\\"Percentage of bases\\\"");     
     strcpy(KKK5,"plot [ 1 to 1000 ] [ 0.0001 to 16.0 ] ");
//     strcpy(KKK7,"[ 1 to 300 ] ");
     strcpy(KKK6,"with lines ls 1,");     
     strcpy(KKK7,"with lines ls 2");     
     num_hits =0;
     k = 0;
     offset = 0;

     memset(Sam_name,'\0',30);
     strcpy(Sam_name,sample);
     memset(Chr_name,'\0',30);
     strcpy(Chr_name,chromo);

     sprintf(outputname,"%s.freq",chromo);
     for(i=0;i<(nSeq-1);i++)
     {
        stopflag=0;
        j=i+1;
        num_blocks = 0;
        while((j<nSeq)&&(stopflag==0))
        {
          if(strcmp(S_Name[j],chromo)==0)
          {
            j++;
          }
          else
            stopflag=1;
        }
        num_hits = j-i;
        num_reads++;
        if((num_hits > 2)) 
        {
          for(k=i;k<j;k++)
	  {
	     hit_cover[k-i] = ctg_cover[k];
             hit_index[k-i] = k-i;
	  }
          ArraySort_Int2(num_hits,hit_cover,hit_index);
          Frequency_Process(outputname,num_hits);
        }
        else
        {
        }
	num_hits = j-i;
        i=j-1;
     }
     
     if((namef = fopen(argv[args+1],"w")) == NULL)
     {
       printf("ERROR main:: args \n");
       exit(1);
     }

     G_Size = 1000;
     fprintf(namef,"#!/bin/sh\n");
     fprintf(namef,"\n"); 
     fprintf(namef,"function plotcmd\n");
     fprintf(namef,"\{\n");
     fprintf(namef,"printf \"%s\\n\"\n",KKK00);
     fprintf(namef,"printf \"%s\\n\"\n",KKK0);
     fprintf(namef,"printf \"%s\\n\"\n",KKK1);
     fprintf(namef,"printf \"%s\\n\"\n",KKK2);
     fprintf(namef,"printf \"%s\\n\"\n",KKK22);
     fprintf(namef,"printf \"%s\\n\"\n",KKK3);
     fprintf(namef,"printf \"%s\\n\"\n",KKK4);
//     fprintf(namef,"printf \"%s \\\"%s%s\\\" title \\\"%s %s \\\" %s\"\n",KKK5,chromo,".freq",Sam_name,Chr_name,KKK6);
//     fprintf(namef,"printf \"%s \\\"%s%s\\\" title \\\"%s %s \\\" %s\"\n \\\"%s\\\" title \\\"%s %s \\\" %s\"\n",KKK5,chromo,".freq",Sam_name,Chr_name,KKK6,"coverage.freq",Sam_name,"All-data",KKK7);
     fprintf(namef,"printf \"%s \\\"%s%s\\\" title \\\"%s %s \\\" %s\\\"%s\\\" title \\\"%s %s \\\" %s\"\n",KKK5,chromo,".freq",Sam_name,Chr_name,KKK6,"coverage.freq",Sam_name,"Whole-Genome",KKK7);
//     fprintf(namef,"printf \"%s \\\"%s%s\\\" title \\\"%s %s \\\" %s\ \\\"%s\\\" title \\\"%s %s \\\" %s\"\n",KKK5,chromo,".freq",Sam_name,Chr_name,KKK6,"coverage.freq",Sam_name,"All-data",KKK7);
     fprintf(namef,"}\n");
     fprintf(namef,"plotcmd | gnuplot > data.svg\n");
     fprintf(namef,"inkscape -z --export-text-to-path --export-pdf data.pdf data.svg\n");
     fprintf(namef,"gs -r600 -dNOPAUSE -dBATCH -sDEVICE=png256 -sOutputFile=%s%s-freq.png data.pdf\n",Sam_name,Chr_name);
     fprintf(namef,"\n"); 
     fprintf(namef,"\n"); 
     fclose(namef);
      
     memset(syscmd,'\0',2000);
     sprintf(syscmd,"bash %s > try.out",argv[args+1]); 
     if(system(syscmd) == -1)
     {
        printf("System command error:\n");
     }
}

/*   subroutine to calculate frequency   */
/* =============================== */
void Frequency_Process(char *plotname,int nSeq)
/* =============================== */
{
     int i=0,j=0,k,len=0,num_steps;
     int BAR = 0,nstep = 0,stopflag;
     char line[100],tempc1[100];
     FILE *namef,*namef2;
     long num_base,base;
     double rate;

    if((namef = fopen(plotname,"w")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }
    num_steps = 0;
    BAR = 1;
    nstep = 1;
    j = 0;
    for(i=0;i<nSeq;i++)
    {
       if((hit_cover[j]==0))
       {
         j++;
       }
    }
    if(j>=1)
    {
      rate = j;
      rate = rate/nSeq;
      rate = rate*100;
      printf("frequency: %d %lf\n",1,rate);
      fprintf(namef,"%d %lf\n",1,rate);
    }
    for(i=0;i<nSeq;i++)
    {
/*     search reads with an index < i     */
/*     search reads with an index > i     */
       stopflag=0;
       j=i+1;
       base = hit_cover[i];
       while((j<nSeq)&&(stopflag==0))
       {
         if((hit_cover[j]<(BAR+nstep))&&(hit_cover[i]>=BAR))
         {
           base = base + hit_cover[j];
           j++;
         }
         else
           stopflag=1;
       }
       if((j-i)>1)
       {
         rate = (j-i);
         rate = rate/nSeq;
         rate = rate*100;
         fprintf(namef,"%d %lf\n",BAR+nstep,rate);
         BAR = BAR+nstep;
         num_steps++;
       }
       else if((j-i)==1)
       {
         rate = 100;
         rate = rate/nSeq;
         num_steps++;
       }
       i=j-1;
     }
     fclose(namef);
}

#define SWAP(a,b) temp=(a);(a)=b;(b)=temp;

/*   Subroutine to sort an array arr[0,...,n-1] into ascending order while
     making the corresponding reaarangement of the array brr[0,...,n-1]
     by the use of Quicksort (Sedgwick, R. 1978, Communications o fthe ACM,
     vol. 21, pp. 847-857) also see Numerical Recipes in C                  */  

/* =============================== */
void ArraySort_Long(int n, long *arr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,NSTACK=50,istack[NSTACK];
     long a,temp,MIN=7;

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
void ArraySort_Mix(int n, long *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     long a,temp,MIN=7;

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
void ArraySort_Mix3(int n, long *arr, int *brr, int *crr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,c,NSTACK=50,istack[NSTACK];
     long a,temp,MIN=7;

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


/* creat an int matrix with subscript ange m[nrl...nrh][ncl...nch]  */
int     **imatrix(long nrl,long nrh,long ncl,long nch)
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
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
char    **cmatrix(long nrl,long nrh,long ncl,long nch)
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
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

