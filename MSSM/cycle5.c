/*
    This is a test program which should reproduce numbers
          presented in Table 5.
*/
                                                                      
                  
#define RGE suspect
     /* choose 'suspect','isajet','softSusy','spheno'*/

#define SUGRA
/*#define AMSB  */


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
         

struct { char *lbl; double m0,mhf,A0,tb,sgn;}
testPoints[9]=
{
{"AP",  130, 600,   0, 5, 1},
{"BP",   70, 250,-300,10, 1},
{"CP",   90, 400,   0,10, 1},
{"DP",  120, 500,-400,10,-1},
{"IP",  180, 350,   0,35, 1},
{"KP", 2500, 550, -80,40,-1},
{"MP", 1100,1100,   0,50, 1},
{"NUG",1620, 300,   0,10, 1},
{"NUH", 250, 530,   0,30, 1}
};

struct { char *lbl; double s_0, s_pi;}
scalarData[3]={
{"A",35,55},
{"B",35,70},
{"C",40,55}
};

struct { char *lbl; double q[3];}
vectorData[2]={
{"A'",{-0.48, 0.78, -0.15 }},
{"B'",{-0.427,0.842,-0.085}}
};


int I_,J;

PRINTRGE(RGE);

for(I_=0;I_<9;I_++) 
{  
  printf("================================================ %s ======\n",testPoints[I_].lbl);

{
  double m0,mhf,a0,tb;
  double gMG1, gMG2, gMG3,  gAl, gAt, gAb,  sgn, gMHu,  gMHd,
         gMl2, gMl3, gMr2, gMr3, gMq2, gMq3, gMu2, gMu3, gMd2, gMd3;
         
  printf("\n========= mSUGRA scenario =====\n");

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

   if(I_==7){
/*==== Non universal SUGRA====*/
  gMG1=mhf, gMG2=mhf,gMG3=0.7*mhf;
}

   
   if(I_==8){
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
  if(err) { printf("Can't calculate %s\n",mess); return 1;}  

{
  printf("\n=== MASSES OF HIGGS AND SUSY PARTICLES: ===\n");
  printHiggs(stdout);
  printMasses(stdout,1);
}


{ int fast=1;
  double Beps=1.E-5, cut=0.01;
  double Omega,Xf;   
  printf("\n==== Calculation of relic density =====\n");  
  if(strcmp(mess,"~o1")) printf(" ~o1 is not LSP\n"); 
                    else o1Contents(stdout);
  Omega=darkOmega(&Xf,fast,Beps);
  printf("Xf=%.2e Omega=%.2e\n",Xf,Omega);
//  printChannels(Xf,cut,Beps,1,stdout);
}


for(J=0;J<3;J++)
{
  double   ffS0P[3],ffS0N[3];
  double pA0[2],pA5[2],nA0[2],nA5[2];
  double Nmass=0.939; /*nucleon mass*/
  double SCcoeff;        

  printf("=================================== Scalar form factors %s\n",  scalarData[J].lbl);

      
  calcScalarFF(0.553,18.9,scalarData[J].s_pi,scalarData[J].s_0);
  nucleonAmplitudes(FeScLoop, pA0,pA5,nA0,nA5);

  SCcoeff=4/M_PI*3.8937966E8*pow(Nmass*Mcdm/(Nmass+ Mcdm),2.);
    
    printf(" proton  SI %.3E  ",SCcoeff*pA0[0]*pA0[0]);
    printf(" neutron SI %.3E \n",SCcoeff*nA0[0]*nA0[0]);
}  

for(J=0;J<2;J++)
{ int i;
  double   ffV0P[3],ffV0N[3];
  double pA0[2],pA5[2],nA0[2],nA5[2];
  double Nmass=0.939; /*nucleon mass*/
  double SCcoeff;        

  printf("=================================== Vector form factors %s\n",  vectorData[J].lbl);

//  for(i=0;i<3;i++) ffV0P[i]=vectorData[J].q[i];
  pVectorFFNu=pVectorFFPd=vectorData[J].q[0];
  pVectorFFNd=pVectorFFPu=vectorData[J].q[1];
  pVectorFFNs=pVectorFFPs=vectorData[J].q[2];

      
//  setProtonFF(NULL,ffV0P, NULL);
//  setNeutronFF(NULL,ffV0N,NULL);

  nucleonAmplitudes(FeScLoop, pA0,pA5,nA0,nA5);

  SCcoeff=4/M_PI*3.8937966E8*pow(Nmass*Mcdm/(Nmass+ Mcdm),2.);
    
    printf(" proton  SD %.3E  ",3*SCcoeff*pA5[0]*pA5[0]);
    printf(" neutron SD %.3E\n",3*SCcoeff*nA5[0]*nA5[0]);
}  

{ int i;
  double   ffS0P[3],ffS0N[3];
  double   ffV0P[3],ffV0N[3];  

  calcScalarFF(0.553,18.9,scalarData[0].s_pi,scalarData[0].s_0);

  pVectorFFNu=pVectorFFPd=vectorData[1].q[0];
  pVectorFFNd=pVectorFFPu=vectorData[1].q[1];
  pVectorFFNs=pVectorFFPs=vectorData[1].q[2];
        
  
//  for(i=0;i<3;i++) ffV0P[i]=vectorData[1].q[i];
//  ffV0N[0]= ffV0P[1];  ffV0N[1]= ffV0P[0]; ffV0N[2]= ffV0P[2];

//  setProtonFF(ffS0P,ffV0P, NULL);
//  setNeutronFF(ffS0N,ffV0N,NULL);

}



{ double dNdE[200];
  double nEvents;
  double rho=0.3; /* DM density GeV/sm^3 */

printf("\n======== Direct Detection ========\n");    
/*  SetFermi(1.2,0,0.5); */

  nEvents=nucleusRecoil(Maxwell,73,Z_Ge,J_Ge73,S00Ge73,S01Ge73,S11Ge73,FeScLoop,dNdE);

  printf("73Ge: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));                                   

  nEvents=nucleusRecoil(Maxwell,131,Z_Xe,J_Xe131,S00Xe131,S01Xe131,S11Xe131,FeScLoop,dNdE);

  printf("131Xe: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));                                   

}


}


                         
  return 0;
}
