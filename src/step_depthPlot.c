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
static int *ctg_locus,*hit_index,*hit_score,*ctg_locnoi,*ctg_cover,*ctg_offset;

/* SSAS default parameters   */
static int IMOD=0;
static int n_type=0;
static int barreads=5;
static int file_flag=2;
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
       else if(!strcmp(argv[i],"-hight"))
       {
         sscanf(argv[++i],"%d",&y_hight);
         edge_flag=1;
         args=args+2;
       }
       else if(!strcmp(argv[i],"-denoise"))
       {
         sscanf(argv[++i],"%d",&denoise_flag);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-reads"))
       {
         sscanf(argv[++i],"%d",&barreads);
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
    if((ctg_locus = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - ctg_left\n");
      exit(1);
    }
    if((ctg_locnoi = (int *)calloc(nseq,sizeof(int))) == NULL)
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
//    printf("Plot: %s %s\n",S_Name[i],chrname);
	}
        i++;
    }
    fclose(namef);
 
    printf("Plot: %s %s\n",chromo,sample);
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
     int i,j,k,kkk,m,n_chr,n_oando,win_size,step_n;
     int num_hits,num_pairs,wgs_depth,num_trans,num_sends,num_reads;
     int stopflag,*chr_trans,*chr_rdhit,*chr_inchr;
     int offset,offset1,offset2,block_len,block_ttt,de_noise[1001],de_index[1001];
     long int wgs_all;
     float pca_sumx,pca_sumy,xdeg,*pca_x,*pca_y;
     FILE *namef,*namef2;
     int de_max,de_min,id_max,id_min,sum_locs,*pca_locus,*pca_index;
     void ArraySort_Mix(int n, long *arr, int *brr);
     char **DBname,RC[2],outputname[100],outpcaname[100];
     void ArraySort_Int2(int n, int *arr, int *brr);
     void ArraySort_Float2(int n, float *arr, int *brr);
     char **cmatrix(long nrl,long nrh,long ncl,long nch);
     void ArraySort_String(int n,char **Pair_Name,int *brr);
     char KKK1[100],KKK10[100],KKK2[100],KKK3[100],KKK4[100],KKK5[100],KKK6[100],KKK7[100];
     char *st,*ed,Sam_name[30],syscmd[2000],Chr_name[30];

     n_oando = 1;

     if((pca_index = (int *)calloc(10000,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_left\n");
       exit(1);
     }
     if((pca_locus = (int *)calloc(10000,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_left\n");
       exit(1);
     }
     if((pca_x = (float *)calloc(10000,sizeof(float))) == NULL)
     {
       printf("fmate: calloc - ctg_left\n");
       exit(1);
     }
     if((pca_y = (float *)calloc(10000,sizeof(float))) == NULL)
     {
       printf("fmate: calloc - ctg_left\n");
       exit(1);
     }

     if(denoise_flag == 3)
     {
       strcpy(KKK10,"set logscale y"); 
       y_hight = 10;    
     }
     strcpy(KKK1,"set terminal svg");     
     if(denoise_flag <= 1)
       strcpy(KKK2,"set style line 1 lt 1 lw 3 pt 7 linecolor rgb \\\"black\\\"");     
     else
       strcpy(KKK2,"set style line 1 lt 1 lw 2 pt 1 linecolor rgb \\\"black\\\"");     
     strcpy(KKK3,"set xlabel \\\"Chromosome coordinate\\\"");     
     if(denoise_flag == 3)
       strcpy(KKK4,"set ylabel \\\"Copy numbers (dimensionless)\\\"");     
     else
       strcpy(KKK4,"set ylabel \\\"Base coverage\\\"");     
     strcpy(KKK5,"plot [ 1 to");
//     strcpy(KKK7,"[ 1 to 150 ] ");
     if(denoise_flag == 3)
       strcpy(KKK7,"[ 0.6 to 6.0 ] ");
     else
       sprintf(KKK7,"%s%d%s","[ 1 to ",y_hight," ] ");
//     strcpy(KKK7,"[ 1 to 300 ] ");
     strcpy(KKK6,"with lines ls 1");     
     num_hits =0;
     k = 0;
     offset = 0;

     memset(Sam_name,'\0',30);
     strcpy(Sam_name,sample);
     memset(Chr_name,'\0',30);
     strcpy(Chr_name,chromo);

     sprintf(outputname,"%s.dat",chromo);
     sprintf(outpcaname,"%s.pca",chromo);
     if((namef = fopen(outputname,"w")) == NULL)
     {
       printf("ERROR main:: args \n");
       exit(1);
     }
     if((namef2 = fopen(outpcaname,"w")) == NULL)
     {
       printf("ERROR main:: args \n");
       exit(1);
     }

     wgs_depth = 0;
     wgs_all = 0;
     for(i=0;i<nSeq;i++)
        wgs_all = wgs_all + ctg_cover[i];

     if(nSeq > 0)
       wgs_depth = wgs_all/nSeq;
     else
     {
       printf("No input data! \n");
       exit(1);
     }

     wgs_depth = wgs_depth*0.5;
     printf("Depth: %d \n",wgs_depth);

     for(i=0;i<(nSeq-1);i++)
     {
        stopflag=0;
        j=i+1;
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
        if((num_hits > 20)) 
        {
            G_Size = 0;
	    win_size = 0;
            for(k=(i+1);k<j;k++)
            {
               if(ctg_offset[k] > G_Size)
                 G_Size = ctg_offset[k];
            }
	    win_size = G_Size/num_hits;
	    win_size = win_size/10;
	    win_size = win_size*10;
	    step_n = num_hits/9000;
	    kkk = 0;
	    xdeg = 3.14159/180;
	    if(denoise_flag == 0)
	    {
	      for(k=(i+1);k<j;k++)
                 fprintf(namef,"%d %d\n",ctg_offset[k],ctg_cover[k]);
//               fprintf(namef,"%d %d\n",ctg_offset[k],ctg_cover[k]);
	    }
	    else if(denoise_flag == 1)
	    {
	      for(k=(i+1);k<(j-1001);k++)
              {
                 memset(de_noise,'\0',4004);
                 memset(de_index,'\0',4004);
                 de_max = 0;
                 de_min = 999999999;
                 id_max = 0;
                 id_min = 0;
                 sum_locs = 0;
                 for(m=0;m<1000;m++)
                 {
                    de_noise[m] = ctg_cover[k+m];
                    de_index[m] = m;
                 }
                 ArraySort_Int2(1000,de_noise,de_index);
                 for(m=400;m<600;m++)
                 {
                    sum_locs = de_noise[m]+sum_locs;
                 }
                 ctg_locnoi[k] = sum_locs/200;
              }

              for(k=(i+1);k<(j-1001);k++)
              {
//               fprintf(namef,"%d %d\n",ctg_offset[k],ctg_cover[k]);
                 fprintf(namef,"%d %d\n",ctg_offset[k],ctg_locnoi[k]);
              }
	    }
	    else if(denoise_flag == 3)
	    {
	      for(k=(i+1);k<(j-1001);k++)
              {
                 float cover_ave = 0.0; 
                 memset(de_noise,'\0',4004);
                 memset(de_index,'\0',4004);
                 de_max = 0;
                 de_min = 999999999;
                 id_max = 0;
                 id_min = 0;
                 sum_locs = 0;
                 for(m=0;m<1000;m++)
                 {
                    de_noise[m] = ctg_cover[k+m];
                    de_index[m] = m;
                 }
                 ArraySort_Int2(1000,de_noise,de_index);
                 for(m=400;m<600;m++)
                 {
                    sum_locs = de_noise[m]+sum_locs;
                 }
		 if(wgs_depth == 0)
		   cover_ave = 0.0;
		 else
		 {
		   cover_ave = sum_locs/200;
		   cover_ave = cover_ave/wgs_depth;
		 }
		 if(cover_ave < 6.0)
                   ctg_locus[k] = sum_locs/200;
		 else 
                   ctg_locus[k] = 2*wgs_depth;
              }

	      for(k=(i+1);k<(j-1001);k++)
              {
                 memset(de_noise,'\0',4004);
                 memset(de_index,'\0',4004);
                 de_max = 0;
                 de_min = 999999999;
                 id_max = 0;
                 id_min = 0;
                 sum_locs = 0;
                 for(m=0;m<1000;m++)
                 {
                    de_noise[m] = ctg_locus[k+m];
                    de_index[m] = m;
                 }
                 ArraySort_Int2(1000,de_noise,de_index);
                 for(m=450;m<550;m++)
                 {
                    sum_locs = de_noise[m]+sum_locs;
                 }
                 ctg_locnoi[k] = sum_locs/100;
              }

              for(k=(i+1);k<(j-1001);k++)
              {
//               fprintf(namef,"%d %d\n",ctg_offset[k],ctg_cover[k]);
                 float cover_ave = 0.0; 
		 if(wgs_depth == 0)
		   cover_ave = 0.0;
		 else
		 {
		   cover_ave = ctg_locnoi[k];
		   cover_ave = cover_ave/wgs_depth;
		 }
                 fprintf(namef,"%d %f\n",ctg_offset[k],cover_ave);
		 if(((k-i)%step_n==0)&&(kkk<9000))
		 {
		   kkk++;
		   xdeg = 3.14159/180;
		   xdeg = xdeg*kkk;
		   xdeg = xdeg/100;
//		   pca_x[kkk-1] = cover_ave*cos(xdeg);
//		   pca_y[kkk-1] = cover_ave*sin(xdeg);
		   pca_x[kkk-1] = cover_ave;
		   pca_locus[kkk-1] = ctg_offset[k];
//		   pca_y[kkk-1] = (2.0*wgs_depth)-cover_ave;
//                   printf("PCA: %s %d %d %f %f %d %f || %f %f\n",chromo,k-i,kkk,xdeg,sin(xdeg),ctg_offset[k],cover_ave,cover_ave*sin(xdeg),cover_ave*cos(xdeg));
		 }
              }
/*	      for(k=0;k<9000;k++)
	      {
		 if(pca_x[k] < 7.0)
		   pca_xx[k] = pca_x[k]*1000;
		 else
		   pca_xx[k] = 2000;
		 pca_index[k] = k;
	      }   */
	      pca_sumx = 0;
	      pca_sumy = 0;
	      for(k=0;k<9000;k++)
	      {
		 if(pca_x[k] < 7.0)
		   pca_sumx = pca_sumx + pca_x[k];
//		 pca_sumy = pca_sumy + pca_y[k]*pca_y[k];
	      }
	      pca_sumx = pca_sumx/9000;
	      for(k=0;k<9000;k++)
              {
		 pca_y[k] = (pca_x[k]-1.9)*(pca_x[k]-1.9);
//		 pca_y[k] = (pca_x[k]-pca_sumx)*(pca_x[k]-pca_sumx);
//		 pca_y[k] = sqrt(pca_y[k]);
		 if(pca_y[k] < 1.0)
                   pca_sumy = pca_sumy + pca_y[k];
              }
	      pca_sumy = pca_sumy/9000;
	      pca_sumy = sqrt(pca_sumy);
              fprintf(namef2,"PCA: %s %s %f %f\n",Sam_name,chromo,pca_sumx,pca_sumy);
              printf("PCA: %s %s %f %f\n",Sam_name,chromo,pca_sumx,pca_sumy);
//	      for(k=0;k<9000;k++)
//                 printf("CO: %d %d %f %f %f\n",k,pca_locus[k],pca_sumx,pca_x[k],pca_y[k]);
	    }
	    else if(denoise_flag == 2)
	    {
              int cover_ave = wgs_depth;
              int ave = 0;
              int edg1 = cover_ave/3;
              int edg2 = cover_ave*1.5;
              int edg3 = cover_ave*2.5;
              int edg4 = cover_ave*3.5;
              int edg5 = cover_ave*4.5;
              int edg6 = cover_ave*5.5;
              int s_len;

	      for(k=(i+1);k<(j-1001);k++)
              {
	         memset(de_noise,'\0',4004);
	         memset(de_index,'\0',4004);
                 de_max = 0;
	         de_min = 999999999;
	         id_max = 0;
	         id_min = 0;
	         sum_locs = 0;
    	         for(m=0;m<1000;m++)
	         {
		    de_noise[m] = ctg_cover[k+m];
	            de_index[m] = m;	  
	         }
                 ArraySort_Int2(1000,de_noise,de_index);
	         for(m=400;m<600;m++)
	         {
	            sum_locs = de_noise[m]+sum_locs;
	         }
//	         ctg_locnoi[k] = sum_locs/200;
	         s_len = sum_locs/200;


		 if((s_len >= 15)&&(s_len < edg2))
                 {
                   ave = cover_ave;
                   s_len = ave + (s_len-ave)*0.5;
                 }
                 else if((s_len >= edg2)&&(s_len < edg3))
                 {
                   ave = 2*cover_ave;
                   s_len = ave + (s_len-ave)*0.4;
                 }
                 else if((s_len >= edg3)&&(s_len < edg4))
                 {
                   ave = 3*cover_ave;
                   s_len = ave + (s_len-ave)*0.3;
                 }
                 else if((s_len >= edg4)&&(s_len < edg5))
                 {
                   ave = 4*cover_ave;
                   s_len = ave + (s_len-ave)*0.2;
                 }
                 else if((s_len >= edg5)&&(s_len < edg6))
                 {
                   ave = 5*cover_ave;
                   s_len = ave + (s_len-ave)*0.15;
                 }
                 else if((s_len >= edg6)&&(s_len < 260))
                 {
                   ave = 6*cover_ave;
                   s_len = ave + (s_len-ave)*0.1;
                 }

		 ctg_locnoi[k] = s_len;
	      }

    	      for(k=(i+1);k<(j-1001);k++)
              {
//               fprintf(namef,"%d %d\n",ctg_offset[k],ctg_cover[k]);
                 fprintf(namef,"%d %d\n",ctg_offset[k],ctg_locnoi[k]);
              }
	    }
        }
        else
        {
        }
	num_hits = j-i;
        i=j-1;
     }
     fclose(namef);
     fclose(namef2);

//     G_Size = 250000000;
     printf("Max: %d \n",G_Size);
     if((namef = fopen(argv[args+1],"w")) == NULL)
     {
       printf("ERROR main:: args \n");
       exit(1);
     }

     fprintf(namef,"#!/bin/sh\n");
     fprintf(namef,"\n"); 
     fprintf(namef,"function plotcmd\n");
     fprintf(namef,"\{\n");
     if(denoise_flag == 3)
     {
       fprintf(namef,"printf \"%s\\n\"\n",KKK10);
     }
     fprintf(namef,"printf \"%s\\n\"\n",KKK1);
     fprintf(namef,"printf \"%s\\n\"\n",KKK2);
     fprintf(namef,"printf \"%s\\n\"\n",KKK3);
     fprintf(namef,"printf \"%s\\n\"\n",KKK4);
     fprintf(namef,"printf \"%s %d ] %s \\\"%s%s\\\" title \\\"%s %s \\\" %s\"\n",KKK5,G_Size,KKK7,chromo,".dat",Sam_name,Chr_name,KKK6);
     fprintf(namef,"}\n");
     fprintf(namef,"plotcmd | gnuplot > data.svg\n");
     fprintf(namef,"inkscape -z --export-text-to-path --export-pdf data.pdf data.svg\n");
     fprintf(namef,"gs -r600 -dNOPAUSE -dBATCH -sDEVICE=png256 -sOutputFile=%s%s.png data.pdf\n",Sam_name,Chr_name);
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

/* =============================== */
void ArraySort_Float2(int n, float *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     int MIN=7;
     float a,temp;

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

