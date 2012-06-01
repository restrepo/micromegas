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


int asgn_(int numvar,double newval) { }
int  calcFunc(void){}




int main(void)
{

 double be1=0,be2=0;

 fillCutArray();
 fillRegArray();

 imkmom(be1,be2, &ndim,&cmMom);


do{

      L0=L; for (i =dim-1;i>=0; i--) {Kg[i]=L0%Ng[i]; L0=L0/Ng[i];}
      seedXX(drandXXstate);
      for (i =0;i<dim;i++) xlocal[i]=drandXX();
       Local2Global(xlocal,x,&f,NULL);

       mkmom(x, &factor_0);
 
       lorrot(rapidity,nin_+nout_);
  
  
  
  } while(1)

}
