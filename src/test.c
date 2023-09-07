/****************************************************************************
 ****************************************************************************
 *                                                                          *
 *  Copyright (C) 2017  Genome Research Ltd.                                *
 *                                                                          *
 *  Author: Zemin Ning (zn1@sanger.ac.uk)                                   *
 *                                                                          *
 *  This file is part of covidPileup pipeline.                              *
 *                                                                          *
 *  covidPileup is a free software: you can redistribute it and/or modify it*
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
#include "fasta.h"

#define GT '>'
#define GT4 (((((GT<<8)+GT)<<8)+GT)<<8)+GT

#define ENDS_EXTRA 0
#define PADCHAR '-'
#define MAX_N_BRG 50000 
#define MAX_N_ROW 50000 
#define Max_N_NameBase 100
static char **S_Name,**R_Name,**T_Name,**M_Name;
static int *snp_locus,*hit_length;

/* SSAS default parameters   */
static int IMOD=0;
static int n_type=0;
static int barreads=10;
static int file_flag=2;
static int n_cover=5;
static int edge_set=2000;
static int edge_flag=0;
static int nContig=0;
static int max_len = 100000;
typedef struct
{
       int foffset;
       int fsindex;
} SIO;

fasta *expt;


int main(int argc, char **argv)
{
    FILE *namef,*namef2;
    int i,j,nSeq,args,n_pileup,g_size,n_SNPs,n_nation;
    int *cod_list,*cod_head,**snp_list,**snp_head;
    int n,*snp_index,*hit_index,*snp_sample,i_country;
    int n_contig,n_reads,num_uniqs,num_allSNPs,nseq;
    fasta *seq;
    void Memory_Allocate(int arr);
    void ArraySort_String(int n,char **Pair_Name,int *brr);
    char line[2000]={0},tempc1[60],cc[60],RC[2],rdname[100],countryname[100],nationname[100],*st,*ed;
    char **cmatrix(long nrl,long nrh,long ncl,long nch);
    int **imatrix(long nrl,long nrh,long ncl,long nch);

    seq=NULL;
    if(argc < 2)
    {
      printf("Usage: %s -cover 5 -country UK <input_countryname_file> <input_SNP_file> <output_countrySNP_unique> <output_countrySNP_all>\n",argv[0]);

      exit(1);
    }

    strcpy(countryname,"UK");
    memset(nationname,'\0',100);

    nSeq=0;
    args=1;
    for(i=1;i<argc;i++)
    {
       if(!strcmp(argv[i],"-mod"))
       {
         sscanf(argv[++i],"%d",&IMOD); 
         args=args+2;
       }
       else if(!strcmp(argv[i],"-cover"))
       {
         sscanf(argv[++i],"%d",&n_cover);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-country"))
       {
         sscanf(argv[++i],"%s",countryname);
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

    if((snp_sample = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - snp_index\n");
      exit(1);
    }

    R_Name=cmatrix(0,nseq+10,0,Max_N_NameBase);

    n_nation = nseq;

    if((namef = fopen(argv[args],"r")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }

/*  read the alignment files         */
    i=0;
    while(fscanf(namef,"%s %d %s",tempc1,&snp_sample[i],R_Name[i])!=EOF)
    {
        i++;
    }
    fclose(namef);

    i_country = 305;
    if(strcmp(countryname,"UK")==0)
      strcpy(nationname,"England");
    else if(strcmp(countryname,"CHINA")==0)
      strcpy(nationname,"China");
    else if(strcmp(countryname,"EU")==0)
      strcpy(nationname,"Belgium");
    else
      strcpy(nationname,countryname);

    for(i=0;i<n_nation;i++)
    {
       
       if(strcmp(R_Name[i],nationname)==0)
         i_country = i;
    }
    if((namef = fopen(argv[args+1],"r")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }
    n_SNPs = 0;
    while(!feof(namef))
    {
      if(fgets(line,2000,namef) == NULL)
      {
//        printf("fgets command error:\n);
      }
      if(feof(namef)) break;
      n_SNPs++;
    }
    fclose(namef); 

    g_size = 40000;
    if((snp_index = (int *)calloc(n_SNPs,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - snp_index\n");
      exit(1);
    }
    if((hit_index = (int *)calloc(n_SNPs,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - snp_index\n");
      exit(1);
    }
    if((snp_locus = (int *)calloc(n_SNPs,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - snp_index\n");
      exit(1);
    }
    if((cod_list = (int *)calloc(g_size,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - snp_index\n");
      exit(1);
    }
    if((cod_head = (int *)calloc(g_size,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - snp_index\n");
      exit(1);
    }

    S_Name   = cmatrix(0,n_SNPs+10,0,Max_N_NameBase);
    snp_list = imatrix(0,g_size,0,n_nation);
    snp_head = imatrix(0,g_size,0,n_nation);

    for(i=0;i<g_size;i++)
    {
       cod_list[i] = 0;
    }

    if((namef = fopen(argv[args+1],"r")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }

    printf("Here \n");
/*  read the alignment files         */
    i=0;
    while(fscanf(namef,"%s %d %s %d %s %s %s",tempc1,&snp_index[i],tempc1,&snp_locus[i],S_Name[i],tempc1,tempc1)!=EOF)
    {
        st = S_Name[i];
        ed = strchr(S_Name[i],'/');
        memset(rdname,'\0',100);
        strncpy(rdname,S_Name[i],ed-st);

             printf("Nation: %d %s %s\n",i,rdname,R_Name[0]);
        if(strcmp(countryname,"UK")==0)
        {
          if((strcmp(rdname,"England")==0)||(strcmp(rdname,"Scotland")==0)||(strcmp(rdname,"Wales")==0)||(strcmp(rdname,"NorthernIreland")==0))
          {
            snp_list[snp_locus[i]][i_country]++;
            cod_list[snp_locus[i]]++;
//          printf("Nation: %d %s %s %d\n",j,rdname,S_Name[i],snp_locus[i]);
          }
          else
          {
            for(j=0;j<n_nation;j++)
            {
               if(strcmp(rdname,R_Name[j])==0)
               {
                 snp_list[snp_locus[i]][j]++;
                 cod_list[snp_locus[i]]++;
               printf("Nation: %d %s %s %d\n",j,rdname,S_Name[i],snp_locus[i]);
               }
            }
          }
        }
        else if(strcmp(countryname,"CHINA")==0)
        {
          int idd = 0;
          if(strcmp(rdname,"Anhui")==0)
            idd = 2;
          else if(strcmp(rdname,"Beijing")==0)
            idd = 6;
          else if(strcmp(rdname,"China")==0)
            idd = 13;
          else if(strcmp(rdname,"Chongqing")==0)
            idd = 14;
          else if(strcmp(rdname,"Fujian")==0)
            idd = 28;
          else if(strcmp(rdname,"Fuyang")==0)
            idd = 29;
          else if(strcmp(rdname,"Ganzhou")==0)
            idd = 31;
          else if(strcmp(rdname,"Guangdong")==0)
            idd = 36;
          else if(strcmp(rdname,"Guangzhou")==0)
            idd = 37;
          else if(strcmp(rdname,"Hangzhou")==0)
            idd = 38;
          else if(strcmp(rdname,"Hefei")==0)
            idd = 39;
          else if(strcmp(rdname,"Henan")==0)
            idd = 40;
          else if(strcmp(rdname,"Jian")==0)
            idd = 51;
          else if(strcmp(rdname,"Jiangsu")==0)
            idd = 52;
          else if(strcmp(rdname,"Jiangxi")==0)
            idd = 53;
          else if(strcmp(rdname,"Jingzhou")==0)
            idd = 54;
          else if(strcmp(rdname,"Jiujiang")==0)
            idd = 55;
          else if(strcmp(rdname,"Lishui")==0)
            idd = 60;
          else if(strcmp(rdname,"NanChang")==0)
            idd = 65;
          else if(strcmp(rdname,"Nanchang")==0)
            idd = 66;
          else if(strcmp(rdname,"Pingxiang")==0)
            idd = 79;
          else if(strcmp(rdname,"Shandong")==0)
            idd = 88;
          else if(strcmp(rdname,"Shanghai")==0)
            idd = 89;
          else if(strcmp(rdname,"Shangrao")==0)
            idd = 90;
          else if(strcmp(rdname,"Shenzhen")==0)
            idd = 91;
          else if(strcmp(rdname,"Sichuan")==0)
            idd = 92;
          else if(strcmp(rdname,"Tianmem")==0)
            idd = 104;
          else if(strcmp(rdname,"Wuhan")==0)
            idd = 110;
          else if(strcmp(rdname,"Wuhan-Hu-1")==0)
            idd = 111;
          else if(strcmp(rdname,"Xinyu")==0)
            idd = 112;
          else if(strcmp(rdname,"Yunnan")==0)
            idd = 113;
          else if(strcmp(rdname,"Zhejiang")==0)
            idd = 114;
          else
            idd = -1;
   
          if(idd > 0)
          {
            snp_list[snp_locus[i]][i_country]++;
            cod_list[snp_locus[i]]++;
//          printf("Nation: %d %s %s %d\n",j,rdname,S_Name[i],snp_locus[i]);
          }
          else
          {
            for(j=0;j<n_nation;j++)
            {
               if(strcmp(rdname,R_Name[j])==0)
               {
                 snp_list[snp_locus[i]][j]++;
                 cod_list[snp_locus[i]]++;
//               printf("Nation: %d %s %s %d\n",j,rdname,S_Name[i],snp_locus[i]);
               }
            }
          }
        }
        else if(strcmp(countryname,"EU")==0)
        {
          int idd = 0;
          if(strcmp(rdname,"Austria")==0)
            idd = 5;
          else if(strcmp(rdname,"Belgium")==0)
            idd = 9;
          else if(strcmp(rdname,"Denmark")==0)
            idd = 23;
          else if(strcmp(rdname,"France")==0)
            idd = 30;
          else if(strcmp(rdname,"Germany")==0)
            idd = 36;
          else if(strcmp(rdname,"Italy")==0)
            idd = 52;
          else if(strcmp(rdname,"Luxembourg")==0)
            idd = 66;
          else if(strcmp(rdname,"Netherlands")==0)
            idd = 73;
          else if(strcmp(rdname,"Portugal")==0)
            idd = 85;
          else if(strcmp(rdname,"Spain")==0)
            idd = 103;
          else
            idd = -1;
   
          if(idd > 0)
          {
            snp_list[snp_locus[i]][i_country]++;
            cod_list[snp_locus[i]]++;
//          printf("Nation: %d %s %s %d %d\n",j,rdname,S_Name[i],snp_locus[i],i_country);
          }
          else
          {
            for(j=0;j<n_nation;j++)
            {
               if(strcmp(rdname,R_Name[j])==0)
               {
                 snp_list[snp_locus[i]][j]++;
                 cod_list[snp_locus[i]]++;
//               printf("Nation: %d %s %s %d\n",j,rdname,S_Name[i],snp_locus[i]);
               }
            }
          }
        }
        else
        {
          for(j=0;j<n_nation;j++)
          {
             if(strcmp(rdname,R_Name[j])==0)
             {
               snp_list[snp_locus[i]][j]++;
               cod_list[snp_locus[i]]++;
//             printf("Nation: %d %s %s %d\n",j,rdname,S_Name[i],snp_locus[i]);
             }
          }
        }
        i++;
    }
    fclose(namef);

      printf("www: %d \n",nseq);
    cod_head[0] = 0;
    for(i=1;i<g_size;i++)
       cod_head[i] = cod_head[i-1] + cod_list[i-1];


    for(i=0;i<30000;i++)
       printf("Cod: %d %d %d \n",i,cod_list[i],cod_head[i]);
    n_pileup = 0;
    for(i=0;i<g_size;i++)
    {
       if(cod_list[i] > 0)
         n_pileup++;
       snp_head[i][0] = 0;
       for(j=1;j<n_nation;j++)
       {
          snp_head[i][j] = snp_head[i][j-1] + snp_list[i][j-1]; 
       }
    }

    for(i=0;i<g_size;i++)
    {
       cod_list[i] = 0;
       for(j=0;j<n_nation;j++)
          snp_list[i][j] = 0;
    }

    if((namef = fopen(argv[args+1],"r")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }

    printf("Here 2 \n");
/*  read the alignment files         */
    i=0;
    while(fscanf(namef,"%s %d %s %d %s %s %s",tempc1,&snp_index[i],tempc1,&snp_locus[i],S_Name[i],tempc1,tempc1)!=EOF)
    {
        st = S_Name[i];
        ed = strchr(S_Name[i],'/');
        memset(rdname,'\0',100);
        strncpy(rdname,S_Name[i],ed-st);

        if(strcmp(countryname,"UK")==0)
        {
          if((strcmp(rdname,"England")==0)||(strcmp(rdname,"Scotland")==0)||(strcmp(rdname,"Wales")==0)||(strcmp(rdname,"NorthernIreland")==0))
          {
             printf("Nation: %d %s %s %d %d %d\n",i,rdname,S_Name[i],snp_index[i],snp_locus[i],i_country);
            int idt = cod_head[snp_locus[i]]+snp_head[snp_locus[i]][i_country]+snp_list[snp_locus[i]][i_country];
          printf("Nation1: %d %s %s %d || %d %d %d\n",idt,rdname,S_Name[i],snp_locus[i],cod_head[snp_locus[i]],snp_head[snp_locus[i]][i_country],snp_list[snp_locus[i]][i_country]);
            hit_index[idt] = i;
          printf("Nation2: %d %s %s %d\n",j,rdname,S_Name[i],snp_locus[i]);
            snp_list[snp_locus[i]][i_country]++;
            cod_list[snp_locus[i]]++;
          printf("Nation3: %d %s %s %d\n",j,rdname,S_Name[i],snp_locus[i]);
          }
          else
          {
            for(j=0;j<n_nation;j++)
            {
               if(strcmp(rdname,R_Name[j])==0)
               {
                 snp_list[snp_locus[i]][j]++;
                 cod_list[snp_locus[i]]++;
//               printf("Nation: %d %s %s %d\n",j,rdname,S_Name[i],snp_locus[i]);
               }
            }
          }
        }
        else if(strcmp(countryname,"CHINA")==0)
        {
          int idd = 0;
          if(strcmp(rdname,"Anhui")==0)
            idd = 2;
          else if(strcmp(rdname,"Beijing")==0)
            idd = 6;
          else if(strcmp(rdname,"China")==0)
            idd = 13;
          else if(strcmp(rdname,"Chongqing")==0)
            idd = 14;
          else if(strcmp(rdname,"Fujian")==0)
            idd = 28;
          else if(strcmp(rdname,"Fuyang")==0)
            idd = 29;
          else if(strcmp(rdname,"Ganzhou")==0)
            idd = 31;
          else if(strcmp(rdname,"Guangdong")==0)
            idd = 36;
          else if(strcmp(rdname,"Guangzhou")==0)
            idd = 37;
          else if(strcmp(rdname,"Hangzhou")==0)
            idd = 38;
          else if(strcmp(rdname,"Hefei")==0)
            idd = 39;
          else if(strcmp(rdname,"Henan")==0)
            idd = 40;
          else if(strcmp(rdname,"Jian")==0)
            idd = 51;
          else if(strcmp(rdname,"Jiangsu")==0)
            idd = 52;
          else if(strcmp(rdname,"Jiangxi")==0)
            idd = 53;
          else if(strcmp(rdname,"Jingzhou")==0)
            idd = 54;
          else if(strcmp(rdname,"Jiujiang")==0)
            idd = 55;
          else if(strcmp(rdname,"Lishui")==0)
            idd = 60;
          else if(strcmp(rdname,"NanChang")==0)
            idd = 65;
          else if(strcmp(rdname,"Nanchang")==0)
            idd = 66;
          else if(strcmp(rdname,"Pingxiang")==0)
            idd = 79;
          else if(strcmp(rdname,"Shandong")==0)
            idd = 88;
          else if(strcmp(rdname,"Shanghai")==0)
            idd = 89;
          else if(strcmp(rdname,"Shangrao")==0)
            idd = 90;
          else if(strcmp(rdname,"Shenzhen")==0)
            idd = 91;
          else if(strcmp(rdname,"Sichuan")==0)
            idd = 92;
          else if(strcmp(rdname,"Tianmem")==0)
            idd = 104;
          else if(strcmp(rdname,"Wuhan")==0)
            idd = 110;
          else if(strcmp(rdname,"Wuhan-Hu-1")==0)
            idd = 111;
          else if(strcmp(rdname,"Xinyu")==0)
            idd = 112;
          else if(strcmp(rdname,"Yunnan")==0)
            idd = 113;
          else if(strcmp(rdname,"Zhejiang")==0)
            idd = 114;
          else
            idd = -1;
   
          if(idd > 0)
          {
            int idt = cod_head[snp_locus[i]]+snp_head[snp_locus[i]][i_country]+snp_list[snp_locus[i]][i_country];
            hit_index[idt] = i;
            snp_list[snp_locus[i]][i_country]++;
            cod_list[snp_locus[i]]++;
//          printf("Nation: %d %s %s %d\n",j,rdname,S_Name[i],snp_locus[i]);
          }
          else
          {
            for(j=0;j<n_nation;j++)
            {
               if(strcmp(rdname,R_Name[j])==0)
               {
                 snp_list[snp_locus[i]][j]++;
                 cod_list[snp_locus[i]]++;
//               printf("Nation: %d %s %s %d\n",j,rdname,S_Name[i],snp_locus[i]);
               }
            }
          }
        }
        else if(strcmp(countryname,"EU")==0)
        {
          int idd = 0;
          if(strcmp(rdname,"Austria")==0)
            idd = 5;
          else if(strcmp(rdname,"Belgium")==0)
            idd = 9;
          else if(strcmp(rdname,"Denmark")==0)
            idd = 23;
          else if(strcmp(rdname,"France")==0)
            idd = 30;
          else if(strcmp(rdname,"Germany")==0)
            idd = 36;
          else if(strcmp(rdname,"Italy")==0)
            idd = 52;
          else if(strcmp(rdname,"Luxembourg")==0)
            idd = 66;
          else if(strcmp(rdname,"Netherlands")==0)
            idd = 73;
          else if(strcmp(rdname,"Portugal")==0)
            idd = 85;
          else if(strcmp(rdname,"Spain")==0)
            idd = 103;
          else
            idd = -1;
   
          if(idd > 0)
          {
            int idt = cod_head[snp_locus[i]]+snp_head[snp_locus[i]][i_country]+snp_list[snp_locus[i]][i_country];
            hit_index[idt] = i;
            snp_list[snp_locus[i]][i_country]++;
            cod_list[snp_locus[i]]++;
//            printf("Nation: %d %s %s %d %d\n",j,rdname,S_Name[i],snp_locus[i],i_country);
          }
          else
          {
            for(j=0;j<n_nation;j++)
            {
               if(strcmp(rdname,R_Name[j])==0)
               {
//                 int idt = cod_head[snp_locus[i]]+snp_head[snp_locus[i]][j]+snp_list[snp_locus[i]][j];
//                 hit_index[idt] = i;
                 snp_list[snp_locus[i]][j]++;
                 cod_list[snp_locus[i]]++;
//               printf("EU: %d %s %s %d\n",j,rdname,S_Name[i],snp_locus[i]);
               }
            }
          }
        }
        else
        {
          for(j=0;j<n_nation;j++)
          {
             if(strcmp(rdname,R_Name[j])==0)
             {
               int idt = cod_head[snp_locus[i]]+snp_head[snp_locus[i]][j]+snp_list[snp_locus[i]][j];
               hit_index[idt] = i;
               snp_list[snp_locus[i]][j]++;
               cod_list[snp_locus[i]]++;
//             printf("Nation: %d %s %s %d\n",j,rdname,S_Name[i],snp_locus[i]);
             }
          }
        }
        i++;
    }
    fclose(namef);


/*
    for(j=0;j<n_nation;j++)
    {
       if(j==i_country)
       {
         for(i=0;i<snp_list[241][j];i++)
         {
            int idd = cod_head[241]+snp_head[241][j]+i;
            int idi = hit_index[idd];
            printf("xxx: %d %s %s %d %d %d\n",j,R_Name[j],S_Name[idi],snp_locus[idi],snp_list[241][j],idi);
         }
       }
    }
             */

    if((namef = fopen(argv[args+2],"w")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }
    if((namef2 = fopen(argv[args+3],"w")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }
    num_uniqs = 0;
    num_allSNPs = 0;
    for(i=0;i<g_size;i++)
    {
       if(cod_list[i] >= n_cover)
       {
         int num_sum = 0;
         for(j=0;j<n_nation;j++)
         {
            if(j != i_country)
              num_sum = num_sum + snp_list[i][j];
         }

         if((num_sum == 0)&&(snp_list[i][i_country] > 0))
         {
           fprintf(namef,"Uniq_%s SNP: %d %d %d\n",countryname,i,cod_list[i],snp_list[i][i_country]);
           num_uniqs++;
           for(j=0;j<cod_list[i];j++)
           {
              int idd = cod_head[i]+j;
              fprintf(namef,"Sample: %d %s %d %d %d\n",j,S_Name[idd],snp_index[idd],snp_locus[idd],snp_locus[cod_head[i]]);
           }
         }

         if((snp_list[i][i_country] >= n_cover))
         {
           
           fprintf(namef2,"All_%s SNP: %d %d %d\n",countryname,i,cod_list[i],snp_list[i][i_country]);
           num_allSNPs++;
           for(j=0;j<snp_list[i][i_country];j++)
           {
              int idd = cod_head[i]+snp_head[i][i_country]+j;
              int idi = hit_index[idd];
              fprintf(namef2,"Sample: %d %s %d %d %d\n",j,S_Name[idi],snp_index[idi],snp_locus[idi],idi);
           }
         }
       }
    }


    fprintf(namef,"==========================================\n");
    fprintf(namef,"Number of Total  SNPs: %d all countries\n",n_pileup);
    fprintf(namef,"Number of Unique SNPs: %d in %s!\n",num_uniqs,countryname);
    fprintf(namef,"==========================================\n");
    fclose(namef);
    fprintf(namef2,"Number of SNPs: %d in %s!\n",num_allSNPs,countryname);
    fclose(namef2);

    return EXIT_SUCCESS;

}
/* end of the main */


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

