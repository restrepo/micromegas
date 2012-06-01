#include<math.h>
#include"../sources/micromegas.h"

struct 
{ double A, Z, J, Sp,Sn; double (*S00)(double);double (*S01)(double);double (*S11)(double);}

sdInfo[19]={
{  19,Z_F  ,J_F19 ,  Sp_F19,   Sn_F19    ,S00F19    ,S01F19    ,S11F19   },
{  23,Z_Na ,J_Na23,  Sp_Na23,  Sn_Na23   ,S00Na23   ,S01Na23   ,S11Na23  },
{  27,Z_Al ,J_Al27,  Sp_Al27,  Sn_Al27   ,S00Al27   ,S01Al27   ,S11Al27  },
{  29,Z_Si ,J_Si29,  Sp_Si29,  Sn_Si29   ,S00Si29   ,S01Si29   ,S11Si29  },
{  39,Z_K  ,J_K39 ,  Sp_K39,   Sn_K39    ,S00K39    ,S01K39    ,S11K39   },
{  73,Z_Ge ,J_Ge73,  Sp_Ge73,  Sn_Ge73   ,S00Ge73   ,S01Ge73   ,S11Ge73  },
{  93,Z_Nb ,J_Nb93,  Sp_Nb93,  Sn_Nb93   ,S00Nb93   ,S01Nb93   ,S11Nb93  },
{  125,Z_Te ,J_Te125,Sp_Te125, Sn_Te125  ,S00Te125  ,S01Te125  ,S11Te125 },
{  127,Z_I  ,J_I127, Sp_I127,  Sn_I127   ,S00I127   ,S01I127   ,S11I127  },
{  129,Z_Xe, J_Xe129,Sp_Xe129, Sn_Xe129  ,S00Xe129  ,S01Xe129  ,S11Xe129 },
{  131,Z_Xe ,J_Xe131,Sp_Xe131, Sn_Xe131  ,S00Xe131  ,S01Xe131  ,S11Xe131 },

{  23,Z_Na ,J_Na23,  Sp_Na23, Sn_Na23  ,S00Na23A   ,S01Na23A   ,S11Na23A  },
{  29,Z_Si ,J_Si29,  Sp_Si29, Sn_Si29  ,S00Si29A   ,S01Si29A   ,S11Si29A  },
{  73,Z_Ge ,J_Ge73,  Sp_Ge73, Sn_Ge73  ,S00Ge73A   ,S01Ge73A   ,S11Ge73A  },
{  125,Z_Te ,J_Te125,Sp_Te125,Sn_Te125 ,S00Te125A  ,S01Te125A  ,S11Te125A },
{  127,Z_I  ,J_I127, Sp_I127, Sn_I127 , S00I127A   ,S01I127A   ,S11I127A  },
{  129,Z_Xe, J_Xe129,Sp_Xe129,Sn_Xe129 , S00Xe129A  ,S01Xe129A  ,S11Xe129A },
{  131,Z_Xe ,J_Xe131,Sp_Xe131,Sn_Xe131 ,S00Xe131A  ,S01Xe131A  ,S11Xe131A },
{  131,Z_Xe ,J_Xe131,Sp_Xe131,Sn_Xe131 ,S00Xe131B  ,S01Xe131B  ,S11Xe131B }

};

void restoreSpSn(double J,double S00,double S11,double S01,double *Sp,double*Sn)
{

double C= (2*J+1)*(J+1)/(4*M_PI*J);
S00/=C;
S11/=C;
S01/=2*C;

*Sp=sqrt(fabs(((S00+S11)/2+S01)/2));
*Sn=sqrt(fabs(((S00+S11)/2-S01)/2));


if((S00-S11)/4<0) *Sp*=-1;

}



#define RGE suspect
     /* choose 'suspect','isajet','softSusy','spheno'*/

/*=========   SUSY scenario  ==========
  One can define SUGRA or AMSB, 
  comment both for low scale input 
=======================================*/ 
#define SUGRA
/*#define AMSB  */

/*====== Modules ===============
   Keys to switch on 
   various modules of micrOMEGAs
================================*/

#define MASSES_INFO      
      /* Display information about SUSY and Higgs masses 
      */
#define CONSTRAINTS     
      /* Display  deltarho, B_>sgamma, Bs->mumu, gmuon and
         check LEP mass limits 
      */ 
#define OMEGA            
      /* Calculate relic density and display contribution of
         individual channels 
      */
#define INDIRECT_DETECTION  
      /* Compute spectra of gamma/positron/neutrinos
         for DM annihilation; calculate <sigma*v> and
         integrate gamma signal over DM galactic squared
         density for given line of sight.  
      */
//#define RESET_FORMFACTORS
      /* Modify default nucleus form factors, 
         DM velocity distribution,
         A-dependence of Fermi-dencity
      */     
#define WIMP_NUCLEON     
      /* Calculate amplitudes and cross-sections for 
         WIMP-mucleon collisions 
      */  
#define WIMP_NUCLEUS      
      /* Calculate number of events for 1kg*day 
         and recoil energy distibution for various nuclei
      */
#define CROSS_SECTIONS 
      /* Calculate cross sections and widths for 
         reactions specified by the user
      */
/*===== end of Modules  ======*/

/*===== Options ========*/
#define SHOWPLOTS
     /* Display  graphical plots on the screen */ 

/*===== End of DEFINE  settings ===== */


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
   int i;
  double E0=15E-6;
  double Sp,Sn;
  
  double Ap=-1.14;
  double An=1;
  double Aplus=Ap+An, Aminus=Ap-An;
  double rho=0.3;   
   
/* to save RGE input/output files uncomment the next line */
/*delFiles(0);*/
     
#ifdef SUGRA
{
  double m0,mhf,a0,tb;
  double gMG1, gMG2, gMG3,  gAl, gAt, gAb,  sgn, gMHu,  gMHd,
         gMl2, gMl3, gMr2, gMr3, gMq2, gMq3, gMu2, gMu3, gMd2, gMd3;
         
  printf("\n========= mSUGRA scenario =====\n");
  PRINTRGE(RGE);

  if(argc<5) 
  { 
    printf(" This program needs 4 parameters:\n"
           "   m0      common scalar mass at GUT scale\n"
           "   mhf     common gaugino mass at GUT scale\n"
           "   a0      trilinear soft breaking parameter at GUT scale\n"
           "   tb      tan(beta) \n");
    printf(" Auxiliary parameters are:\n"
           "   sgn     +/-1,  sign of Higgsino mass term (default 1)\n"    
           "   Mtp     top quark pole mass\n"
           "   MbMb    Mb(Mb) scale independent b-quark mass\n"
           "   alfSMZ  strong coupling at MZ\n");
    printf("Example: ./main 70 250 -300 10\n");                        
      exit(1); 
  } else  
  {  double Mtp,MbMb,alfSMZ;
     sscanf(argv[1],"%lf",&m0);
     sscanf(argv[2],"%lf",&mhf);
     sscanf(argv[3],"%lf",&a0);
     sscanf(argv[4],"%lf",&tb);
     if(argc>5)sscanf(argv[5],"%lf",&sgn); else sgn=1;
     if(argc>6){ sscanf(argv[6],"%lf",&Mtp);    assignValW("Mtp",Mtp);      }
     if(argc>7){ sscanf(argv[7],"%lf",&MbMb);   assignValW("MbMb",MbMb);    }
     if(argc>8){ sscanf(argv[8],"%lf",&alfSMZ); assignValW("alfSMZ",alfSMZ);}
  }

/*==== simulation of mSUGRA =====*/
  gMG1=mhf, gMG2=mhf,gMG3=mhf;
  gAl=a0,   gAt=a0,  gAb=a0;  gMHu=m0,  gMHd=m0;
  gMl2=m0,  gMl3=m0, gMr2=m0, gMr3=m0;
  gMq2=m0,  gMq3=m0, gMu2=m0, gMd2=m0, gMu3=m0, gMd3=m0;

  err= SUGRAMODEL(RGE) (tb,  
    gMG1, gMG2, gMG3,  gAl,  gAt, gAb,  sgn, gMHu, gMHd,
    gMl2, gMl3, gMr2, gMr3, gMq2,  gMq3, gMu2, gMu3, gMd2, gMd3);

  if(err>0){ printf(" non fatal problems in RGE\n");err=0;}  
}
#endif

  if(err) 
  { printf("Problem with RGE solution or spectrum calculation\n");
     exit(1);
  }

  err=sortOddParticles(mess);
  if(err) { printf("Can't calculate %s\n",mess); return 1;}

/*to print input parameters or model parameters in SLHA format
 uncomment correspondingly*/
/* 
  printVar(stdout);  
  writeLesH("slha.out"); 
*/

  
  for(i=0;i<18;i++)
  { double p=sqrt(2*sdInfo[i].A*0.939*E0);  /*/0.197327*; */
    double J=sdInfo[i].J,  A=sdInfo[i].A, Z=sdInfo[i].Z;
    
    double coeff;
    double R;
    double N1,N2;
    double dNdE[250];
    double Emin=5,Emax=50;
    double Sp,Sn;
    
    if(i<11) {Sp=sdInfo[i].Sp; Sn=sdInfo[i].Sn;}
    else{restoreSpSn(J,sdInfo[i].S00(0.),sdInfo[i].S11(0.),sdInfo[i].S01(0.),&Sp,&Sn);}
    
    nucleusRecoil(Maxwell,A,Z,J,sdInfo[i].S00,sdInfo[i].S01,sdInfo[i].S11,0,dNdE);
    N1=cutRecoilResult(dNdE,Emin,Emax);
    
    nucleusRecoil(Maxwell,A,Z,0,sdInfo[i].S00,sdInfo[i].S01,sdInfo[i].S11,0,dNdE);
    N1-=cutRecoilResult(dNdE,Emin,Emax);
    
    nucleusRecoil0(Maxwell,A,Z,J,Sp,Sn,0,dNdE);
    N2=cutRecoilResult(dNdE,Emin,Emax);
    nucleusRecoil0(Maxwell,A,Z,0,Sp,Sn,0,dNdE);
    N2-=cutRecoilResult(dNdE,Emin,Emax);

     printf("For Tab 3  %.0f %.2E  %.2E\n", A, N1,N2); 
  }

}

