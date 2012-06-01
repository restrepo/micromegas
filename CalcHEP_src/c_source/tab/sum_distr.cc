#include <stdio.h>
#include "histogram.h"

#define const
#include"num_out.h"
#undef const

int  nvar_ext=0;
char * varName_ext[1];
double  va_ext[1];
int nfunc_ext=0; 

int nin_ext, nout_ext;
char* pinf_ext(int nsub,int num , double * mass){ return NULL;}

 
int main(int np, char ** par)
{
  int nf,nl,n;
  char fn[30],process[200];
  FILE*f;
  
  if(np!=3) 
  { printf(" This  routine is intended to sum  distributions produced\n" 
           " in  calcHEP numerical sessions (files dist_#).\n"
           " The first and the last session numbers should be submited as parameters\n");
    return 1;
  }

  if(sscanf(par[1],"%d",&nf)!=1){printf("First parameter must be a number\n");return 2;}
  if(sscanf(par[2],"%d",&nl)!=1){printf("Last parameter must be a number\n"); return 2;}

  if(nf>nl) {printf("wrong order of parameters\n"); return 2;}

  for(n=nf;n<=nl;n++)
  { 
    sprintf(fn,"distr_%d",n);
    f=fopen(fn,"r"); 
    if(!f) continue;
    if(add_hist(f,process)) return 3;
    fclose(f);
  } 
  sprintf(fn,"distr_%d_%d",nf,nl);
  f=fopen(fn,"w");
  wrt_hist2(f,process);
  fclose(f);

  return 0;
}
