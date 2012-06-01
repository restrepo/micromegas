/*
 Copyright (C) 1997 by Alexander Pukhov, email: pukhov@theory.npi.msu.su 
*/
#include"n_compil.h"
#include"out_ext.h"
#include"read_func.h"
#include"procvar.h"
#include"getmem.h"
#include"pvars.h"
int gwidth=0;
int rwidth=0;

const int quadra=0;

int pinf_(int nsub,int nprtcl,char * pname,long double * pmass )
{  
  if  (nsub<0 ||nsub>nprc_||nprtcl<0||nprtcl>nin_+nout_) return 1;
  if(pname)strcpy(pname,allcanal[nsub-1].prtclnames[nprtcl-1]);
  if(pmass) *pmass=ABS(*(allcanal[nsub-1]. prtclmasses[nprtcl-1]));
  return 0;
}



int vinf_(int numvar,char * name,double * val)
{
  int k;
  int l=0;

  if(numvar==0 && nin_==2) 
  { if(name)  strcpy(name,"Sqrt(S)");
    if(val) *val=sqrts;
    return 0;
  }  

  if(numvar<0||numvar>nvar_+nfunc_) return 1;

  for (k=1;k<=nmodelvar; k++) 
  if(vararr[k].used && !modelvars[k].func  ) if(++l == numvar) goto exi;

  for (k=1;k<=nmodelvar; k++) 
  if(vararr[k].used && modelvars[k].func  ) if(++l == numvar) goto exi;

  return 1;
  exi:
  if(name)  strcpy(name,vararr[k].alias);
  if(val)   *val=vararr[k].tmpvalue;
  return 0;
}


int asgn_(int numvar,double newval)
{
int k;
int l=0;


  if(numvar==0 && nin_==2) {sqrts=newval; return 0;}

  if(numvar <= 0 || numvar> nvar_ ) return 1;
  
  for (k=1;k<=nmodelvar; k++)
  if(vararr[k].used && !modelvars[k].func) if(++l == numvar)
  { vararr[k].tmpvalue=newval; return 0; }

}



static int  rd_3_(char* s, double * p)
{
   int k;
   
   if(strcmp(s,"i")==0 || strcmp(s,"GG")==0) return 0;

   if (isdigit(*s))
   {  sscanf(s,"%lf",p); return 1; }
                     
   k= modelVarPos(s);
   if(k>0) { *p=vararr[k].tmpvalue;  return 1;}
   return 0;
}



int  calcFunc(void)
{  unsigned  k; 
   int err;

   for (k =1; k <= nmodelvar; k++) if(vararr[k].used && modelvars[k].func )
   {   
 
     err=calcExpression(modelvars[k].func, rd_3_, &vararr[k].tmpvalue ); 
     if (err)  return err;
   }
   return 0;   
} 
