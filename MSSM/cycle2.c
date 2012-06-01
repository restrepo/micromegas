/* 
    This is a test program which should reproduce numbers 
        presented in Table 2.
*/

#define RGE suspect
     /* choose 'suspect','isajet','softSusy','spheno'*/

#define SUGRA

/*====== Modules ===============
   Keys to switch on 
   various modules of micrOMEGAs
================================*/


#include"../sources/micromegas.h"
#include"lib/pmodel.h"

#define SUGRAMODEL_(A) A ## SUGRA
#define SUGRAMODEL(A) SUGRAMODEL_(A)

#define AMSBMODEL_(A) A ## AMSB
#define AMSBMODEL(A) AMSBMODEL_(A)

#define EWSBMODEL_(A) A ## EwsbMSSM
#define EWSBMODEL(A) EWSBMODEL_(A)

#define PRINTRGE_(A) printf(" Spectrum calculator is %s\n", #A)
#define PRINTRGE(A)  PRINTRGE_(A)

int main(int argc,char** argv)
{  int err;
   char mess[10];
   double dMb,dMd;
   double (*LF)(double,double,double,double);         

struct { char *lbl; double m0,mhf,A0,tb,sgn;}
testPoints[4]=
{
{"BP",   70, 250,-300,10, 1},
{"KP", 2500, 550, -80,40,-1},
{"IP",  180, 350,   0,35, 1},
{"NUH", 250, 530,   0,30, 1}
};


int I_, K;

/* to save RGE input/output files uncomment the next line */
/*delFiles(0);*/

PRINTRGE(RGE);

for(I_=0;I_<5;I_++) 
{    


if(I_==4) 
{
   printf("================================================ MSSM1 ======\n");
   err=readVarMSSM("mssm1.par"); 
   if(err==0)err=EWSBMODEL(RGE)();
}else
{
  double m0,mhf,a0,tb;
  double gMG1, gMG2, gMG3,  gAl, gAt, gAb,  sgn, gMHu,  gMHd,
         gMl2, gMl3, gMr2, gMr3, gMq2, gMq3, gMu2, gMu3, gMd2, gMd3;
         
  printf("\n========= mSUGRA scenario =====\n");
  printf("================================================ %s ======\n",testPoints[I_].lbl);

  m0=testPoints[I_].m0;
  mhf=testPoints[I_].mhf;
  a0=testPoints[I_].A0;
  tb=testPoints[I_].tb;
  sgn=testPoints[I_].sgn;
  
/*==== simulation of mSUGRA =====*/
  gMG1=mhf, gMG2=mhf,gMG3=mhf;
  gAl=a0,   gAt=a0,  gAb=a0;  gMHu=m0,  gMHd=m0;
  gMl2=m0,  gMl3=m0, gMr2=m0, gMr3=m0;
  gMq2=m0,  gMq3=m0, gMu2=m0, gMd2=m0, gMu3=m0, gMd3=m0;

  if(I_==3){
/*==== Non universal SUGRA====*/
  gMHu=2.*m0,  gMHd=-0.6*m0;
}


  err= SUGRAMODEL(RGE) (tb,  
    gMG1, gMG2, gMG3,  gAl,  gAt, gAb,  sgn, gMHu, gMHd,
    gMl2, gMl3, gMr2, gMr3, gMq2,  gMq3, gMu2, gMu3, gMd2, gMd3);

  if(err>0){ printf(" non fatal problems in RGE\n");err=0;}  
}

  if(err) 
  { printf("Problem with RGE solution or spectrum calculation\n");
     exit(1);
  }
  err=sortOddParticles(mess);
{
  int fast=1;
  double Beps=1.E-5;
  double Omega,Xf;   
  Omega=darkOmega(&Xf,fast,Beps);
  printf("Omega=%.2e\n",Omega);
}
  if(err) { printf("Can't calculate %s\n",mess); return 1;}

  dMb=findValW("dMb");
  dMd=findValW("dMd");

printf("dMb=%.2E   dMd=dMs=%.2E \n", dMb, dMd);
  
for(K=1;K<=5;K++)
{

switch(K)
{
  case 1: printf("================= Tree level \n"); 
          QCDcorrections=0; 
          assignValW("dMb",0.);  
          assignValW("dMd",0.);
          assignValW("dMs",0.); break;
  case 2: printf("================= QCD        \n"); 
          QCDcorrections=1; 
          assignValW("dMb",0.);
          assignValW("dMd",0.);
          assignValW("dMs",0.); 
             break;
  case 3: printf("================= dMb        \n"); 
           QCDcorrections=0;    
           assignValW("dMb",dMb);
           assignValW("dMd",dMd);
           assignValW("dMs",dMd);
           break;
  case 4: printf("================= QCD+dMb    \n"); 
           QCDcorrections=1;  
           assignValW("dMb",dMb);
           assignValW("dMd",dMd);
           assignValW("dMs",dMd);
           
           break;
  case 5: printf("================= box+QCD+dMb\n"); 
          QCDcorrections=1;  
           assignValW("dMb",dMb);
           assignValW("dMd",dMd);
           assignValW("dMs",dMd);          
          break;
}


{ double pA0[2],pA5[2],nA0[2],nA5[2];
  double Nmass=0.939; /*nucleon mass*/
  double SCcoeff;        

//printf("\n==== Calculation of WIMP-nucleons amplitudes  =====\n");   

if(K==5)  LF=FeScLoop; else LF=NULL;

  nucleonAmplitudes(LF, pA0,pA5,nA0,nA5);
//    printf("WIMP-nucleon amplitudes:\n");
//    printf("proton:  SI  %.3E  SD  %.3E\n",pA0[0],pA5[0]);
//    printf("neutron: SI  %.3E  SD  %.3E\n",nA0[0],nA5[0]); 

  SCcoeff=4/M_PI*3.8937966E8*pow(Nmass*Mcdm/(Nmass+ Mcdm),2.);
//    printf("WIMP-nucleon cross sections:\n");
    
    printf(" proton  SI %.3E  neutron SI %.3E\n",SCcoeff*pA0[0]*pA0[0],SCcoeff*nA0[0]*nA0[0]);
//    printf(" neutron SI %.3E  SD %.3E\n",SCcoeff*nA0[0]*nA0[0],SCcoeff*nA5[0]*nA5[0]);
}
  

}
}
  return 0;
}
