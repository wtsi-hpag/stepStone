/****************************************************************************
 ****************************************************************************
 *                                                                          *
 *  Copyright (C) 2020  Genome Research Ltd.                                *
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

#include <values.h>
#include <stdio.h>
#include <netinet/in.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <ctype.h>

static int IMOD=0;
static int gap_flag = 1;

int main(int argc, char **argv)
{
    int i=0,j=0,k,args=0,num_steps,nSeq,nRead,num_samples;
    int *s_len,BAR = 0,nstep = 0,stopflag;
    char line[500],tempc1[100];
    FILE *namef;
    long num_base,base;
    double rate;

    if(argc < 2)
    {
      printf("Usage: %s -gap 0 <input frequence file> <output frequency file>\n",argv[0]);
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
       else if(!strcmp(argv[i],"-gap"))
       {
         sscanf(argv[++i],"%d",&gap_flag);
         args=args+2;
       }
    }

    if((namef = fopen(argv[args],"r")) == NULL)
    {
      printf("ERROR main:: reads group file: %s \n",argv[args]);
      exit(1);
    }

    nRead = 0;
    while(!feof(namef))
    {
      if(fgets(line,500,namef) == NULL)
        printf("Data_input_file_problem!%s\n",argv[args]);
      if(feof(namef)) break;
      nRead++;
    }
    fclose(namef); 

    if((s_len= (int *)calloc(nRead,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: Align_Process - hit_rcdex\n");
      exit(1);
    }

    if((namef = fopen(argv[args],"r")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }

/*  read the alignment files         */
    i=0;
    num_base = 0;
    num_samples = 0;
    while(fscanf(namef,"%s %d",tempc1,&s_len[i])!=EOF)
    {
       int idd = atoi(tempc1);
       if(idd > num_samples)
	 num_samples = idd;
//         printf("s: %d %d %s %s\n",i,s_len[i],tempc1,argv[1]);
      num_base = num_base + s_len[i];
      i++;
    }
    fclose(namef);

    if((namef = fopen(argv[args+1],"w")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }
    nSeq = i;
    BAR = 1;
    nstep = 1;
    j = 0;
    num_samples = num_samples + 1;
    for(i=0;i<nSeq;i++)
    {
/*     search reads with an index < i     */
/*     search reads with an index > i     */
       stopflag=0;
       j=i+1;
       base = s_len[i]; 
//         printf("www: %d %d %d %d\n",i,s_len[i],s_len[j],BAR);
       while((j<nSeq)&&(stopflag==0))
       {
         if(s_len[j]==s_len[i])
         {
           base = base + s_len[j];
           j++;
         }
         else
           stopflag=1;
       }
       if((j-i)>=1)
       {
         rate = (j-i);
//         rate = rate/num_base;
         if(gap_flag == 0)
           rate = rate/num_samples;
         else
           rate = rate/nSeq;
         rate = rate;
         printf("%d %f %d %d\n",s_len[i],rate,j-i,num_samples);
         fprintf(namef,"%d %d\n",s_len[i],j-i);
         BAR = BAR+nstep;
         num_steps++;
       }
       i=j-1;
     }
     fclose(namef);
}


