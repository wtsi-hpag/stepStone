#include <math.h>
#include <values.h>
#include <stdio.h>
#include <netinet/in.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <ctype.h>


static int *nn;

int main(int argc, char **argv)
{
    int i=0,j=0,k,len=0,num_steps,args,max_len;
    int *s_len,*s_hits,BAR = 0,nstep = 0,stopflag;
    char line[2000],tempc1[100];
    FILE *namef,*namef2;
    long num_base,base,nSeq;
    double rate,lenlog,frelog;
    int num_hits,i_hits,*hit_list,*hit_head;

    if(argc < 2)
    {
      printf("Usage: %s <align.size> <length_frequency>\n",argv[0]);
      exit(1);
    } 
    args=1;
    nSeq=0;
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
      nSeq++;
    }
    fclose(namef);

    if((s_len = (int *)calloc(nSeq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - ctg_left\n");
      exit(1);
    }
    if((hit_list = (int *)calloc(nSeq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - ctg_left\n");
      exit(1);
    }
    if((hit_head = (int *)calloc(nSeq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - ctg_left\n");
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
    max_len = 0;
    while(fscanf(namef,"%s %d",tempc1,&s_len[i])!=EOF)
    {
      if(s_len[i] > max_len)
        max_len = s_len[i];
//         printf("s: %d %d %s %s\n",i,s_len[i],tempc1,argv[1]);
//      num_base = num_base + s_len[i];
      i++;
    }
    fclose(namef);

     for(i=0;i<(nSeq-1);i++)
     {
        stopflag=0;
        j=i+1;
        while((j<nSeq)&&(stopflag==0))
        {
          if(s_len[i]==s_len[j])
          {
            j++;
          }
          else
            stopflag=1;
        }
        num_hits = j-i;
        hit_list[s_len[i]] = j-i;
         printf("sss: %d %d %d\n",i,s_len[i],hit_list[s_len[i]]);
        i = j-1;
    }

    exit(1);

    num_hits = max_len/100;
    num_hits = num_hits+1;
 
    printf("hits: %d %d %ld\n",num_hits,max_len,nSeq);
    if((s_hits = (int *)calloc(num_hits,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - ctg_left\n");
      exit(1);
    }

    for(i=0;i<nSeq;i++)
    {
       if(s_len[i] == 0)
         s_hits[0]++;
       else
       {
         i_hits = (s_len[i]-1)/100; 
         s_hits[i_hits+1]++;
       }
    }

    if((namef2 = fopen(argv[args+1],"w")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }

    for(i=0;i<num_hits;i++)
    {
       rate = s_hits[i];
       rate = rate/nSeq;
       rate = rate*100;
       lenlog = log(s_hits[i]);
       frelog = -log(rate);
       if((lenlog >= 10.0)&&(lenlog <= 15.0))
         fprintf(namef2,"%d %f %f %f\n",i*100+1,rate,lenlog,frelog);
//       if(s_hits[i] >= 10)
//         fprintf(namef2,"%d %f\n",i*100+1,rate);
    }
    return EXIT_SUCCESS;
}


