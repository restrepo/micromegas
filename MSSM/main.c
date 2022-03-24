/*======  Spectrum calculator  ========= 
   Choose RGE from the list below. SuSpect is included 
   in micrOMEGAs, to use another code define the path 
   to the corresponding package in lib/Makefile
=====================================*/ 

#define RGE  spheno

     /* choose 'suspect','softSusy','spheno', 'tree' */

/*=========   SUSY scenario  ==========
  One can define SUGRA, AMSB, EWSB (for low scale input). 
  By default the program reads SLHA data file 
=======================================*/
//#define SUGRA 
//#define SUGRANUH
//#define AMSB 

#define EWSB 

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
//#define CheckMassMatrix      
//#define HIGGSBOUNDS 
//#define HIGGSSIGNALS
//#define SUPERISO
//#define LILITH
//#define SMODELS
//#define MONOJET       

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
//#define LoopGAMMA
      /* Calculate discrete  photon spectrum caused by annihilation of 
         neutralinos into two photons and Z-photon
      */ 
      
//#define RESET_FORMFACTORS
      /* Modify default nucleus form factors, 
         DM velocity distribution,
         A-dependence of Fermi-dencity
      */     
#define CDM_NUCLEON 
      /* Calculate amplitudes and cross-sections for 
         CDM-mucleon collisions 
      */  
//#define TEST_Direct_Detection 
      /* 
        Compare analytical formula for DD against micrOMEGAS calculation.
        As well as compare tree level and box improved approaches.
       */      
#define CDM_NUCLEUS
     // Calculate  exclusion rate for direct detection experiments Xenon1T and DarkSide50

#define NEUTRINO 
 /*  Neutrino signal of DM annihilation in Sun and Earth */
 
#define DECAYS 
      /* Calculate decay widths and branchings  */      
//#define CROSS_SECTIONS 
      /* Calculate cross sections of reactions specified by the user */

/*===== end of Modules  ======*/

/*===== Options ========*/
//#define SHOWPLOTS 
     /* Display  graphical plots on the screen */ 

#define CLEAN    to clean intermediate files

/*===== End of DEFINE  settings ===== */

#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"
#include"lib/pmodel.h"

#define SUGRAMODEL_(A) A ## SUGRA
#define SUGRAMODEL(A) SUGRAMODEL_(A)

#define SUGRANUHMODEL_(A) A ## SUGRAnuh
#define SUGRANUHMODEL(A) SUGRANUHMODEL_(A)

#define AMSBMODEL_(A) A ## AMSB
#define AMSBMODEL(A) AMSBMODEL_(A)

#define EWSBMODEL_(A) A ## EwsbMSSM
#define EWSBMODEL(A) EWSBMODEL_(A)

#define PRINTRGE_(A) printf(" Spectrum calculator is %s\n", #A)
#define PRINTRGE(A)  PRINTRGE_(A)



int main(int argc,char** argv)
{  int err;
   char cdmName[10];
   int spin2, charge3,cdim;

// sysTimeLim=1000; 
   ForceUG=0;   /* to Force Unitary Gauge assign 1 */
   //useSLHAwidth=0;
//  nPROCSS=0; /* to switch off multiprocessor calculations */
//   useSLHAwidth=1;

/*
   if you would like to work with superIso
    setenv("superIso","./superiso_v3.1",1);  
*/

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
/*    printf("Example: ./main 70 250 -300 10\n");  */
      printf("Example: ./main 120 500 -350 10 1 173.1 \n");
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
}
#elif defined(SUGRANUH)
{
  double m0,mhf,a0,tb;
  double gMG1, gMG2, gMG3,  gAl, gAt, gAb,  gMl2, gMl3, gMr2, gMr3, gMq2, gMq3, gMu2, gMu3, gMd2, gMd3,mu,MA;
         
  printf("\n========= mSUGRA non-universal Higgs scenario =====\n");
  PRINTRGE(RGE);

  if(argc<7) 
  { 
    printf(" This program needs 6 parameters:\n"
           "   m0      common scalar mass at GUT scale\n"
           "   mhf     common gaugino mass at GUT scale\n"
           "   a0      trilinear soft breaking parameter at GUT scale\n"
           "   tb      tan(beta) \n" 
           "   mu      mu(EWSB)\n"
           "   MA      mass of pseudoscalar Higgs\n");     
    printf(" Auxiliary parameters are:\n"
           "   Mtp     top quark pole mass\n"
           "   MbMb    Mb(Mb) scale independent b-quark mass\n"
           "   alfSMZ  strong coupling at MZ\n");
/*    printf("Example: ./main 70 250 -300 10\n");  */
      printf("Example: ./main 120 500 -350 10 680 760  \n");
      exit(1); 
  } else  
  {  double Mtp,MbMb,alfSMZ;
     sscanf(argv[1],"%lf",&m0);
     sscanf(argv[2],"%lf",&mhf);
     sscanf(argv[3],"%lf",&a0);
     sscanf(argv[4],"%lf",&tb);
     sscanf(argv[5],"%lf",&mu);
     sscanf(argv[6],"%lf",&MA); 
     if(argc>7){ sscanf(argv[7],"%lf",&Mtp);    assignValW("Mtp",Mtp);      }
     if(argc>8){ sscanf(argv[8],"%lf",&MbMb);   assignValW("MbMb",MbMb);    }
     if(argc>9){ sscanf(argv[9],"%lf",&alfSMZ); assignValW("alfSMZ",alfSMZ);}
  }

/*==== simulation of mSUGRA =====*/
  gMG1=mhf, gMG2=mhf,gMG3=mhf;
  gAl=a0,   gAt=a0,  gAb=a0;
  gMl2=m0,  gMl3=m0, gMr2=m0, gMr3=m0;
  gMq2=m0,  gMq3=m0, gMu2=m0, gMd2=m0, gMu3=m0, gMd3=m0;

  err= SUGRANUHMODEL(RGE) (tb,gMG1,gMG2,gMG3,gAl,gAt,gAb,gMl2,gMl3,gMr2,gMr3,gMq2,gMq3,gMu2,gMu3,gMd2,gMd3,mu,MA); 
}
#elif defined(AMSB)
{
  double m0,m32,sgn,tb;

  printf("\n========= AMSB scenario =====\n");
  PRINTRGE(RGE);
  if(argc<4) 
  { 
    printf(" This program needs 3 parameters:\n"
           "   m0      common scalar mass at GUT scale\n"
           "   m3/2    gravitino mass\n"
           "   tb      tan(beta) \n");
    printf(" Auxiliary parameters are:\n"
           "   sgn     +/-1,  sign of Higgsino mass term (default 1)\n"    
           "   Mtp     top quark pole mass\n"
           "   MbMb    Mb(Mb) scale independent b-quark mass\n"
           "   alfSMZ  strong coupling at MZ\n");
   printf("Example: ./main 450  60000 10\n");                                                                          
   exit(1); 
  } else  
  {  double Mtp,MbMb,alfSMZ;
     sscanf(argv[1],"%lf",&m0);
     sscanf(argv[2],"%lf",&m32);
     sscanf(argv[3],"%lf",&tb);
     if(argc>4)sscanf(argv[4],"%lf",&sgn); else sgn=1;
     if(argc>5){ sscanf(argv[5],"%lf",&Mtp);    assignValW("Mtp",Mtp);      }
     if(argc>6){ sscanf(argv[6],"%lf",&MbMb);   assignValW("MbMb",MbMb);    }
     if(argc>7){ sscanf(argv[7],"%lf",&alfSMZ); assignValW("alfSMZ",alfSMZ);}
  }

  err= AMSBMODEL(RGE)(m0,m32,tb,sgn);
 
}
#elif defined(EWSB)
{ 
   printf("\n========= EWSB scale input =========\n");
   PRINTRGE(RGE);

   if(argc <2) 
   {  printf("The program needs one argument:the name of file with MSSM parameters.\n"
            "Example: ./main mssm1.par \n");
      exit(1);
   }  
   
   printf("Initial file  \"%s\"\n",argv[1]);
     
   err=readVarMSSM(argv[1]);
          
   if(err==-1)     { printf("Can not open the file\n"); exit(2);}
   else if(err>0)  { printf("Wrong file contents at line %d\n",err);exit(3);}

   err=EWSBMODEL(RGE)();
}
#else
{
   printf("\n========= SLHA file input =========\n");

   if(argc <2) 
   {  printf("The program needs one argument:the name of SLHA input file.\n"
            "Example: ./main suspect2_lha.out \n");
      exit(1);
   }  
   
   printf("Initial file  \"%s\"\n",argv[1]);
   err=lesHinput(argv[1]);
   if(err) exit(2);
}
#endif
          
    if(err==-1)     { printf("Can not open the file\n"); exit(2);}
    else if(err>0)  { printf("Wrong file contents at line %d\n",err);exit(3);}
  
  { int nw;
    printf("Warnings from spectrum calculator:\n");
    nw=slhaWarnings(stdout);
    if(nw==0) printf(" .....none\n");
  } 

  if(err) exit(1);
/*  
  initQCD5(0.1184,1.27,4.23, 173.07);
  printf("Mb_pole=%E Mb_pole_1=%E alpha/pi=%E \n", Mbp(), 4.23*(1+4/3*alphaQCD(4.23)/M_PI), alphaQCD(4.23));
  exit(0);  
*/  

  err=sortOddParticles(cdmName);

  if(err) { printf("Can't calculate %s\n",cdmName); return 1;}

//  err=treeMSSM();  


  qNumbers(cdmName,&spin2, &charge3, &cdim);
  printf("\nDark matter candidate is '%s' with spin=%d/2  mass=%.2E\n",
  cdmName,       spin2, Mcdm); 
  
//  if(charge3) { printf("Dark Matter has electric charge %d/3\n",charge3); exit(1);}
  if(cdim!=1) { printf("Dark Matter is a color particle\n"); exit(1);}
  if(strcmp(cdmName,"~o1")) printf(" ~o1 is not CDM\n"); 
                              else o1Contents(stdout);

                             
#ifdef MASSES_INFO
{
  printf("\n=== MASSES OF HIGGS AND SUSY PARTICLES: ===\n");
  printHiggs(stdout);
  printMasses(stdout,1);  
}
#endif


#ifdef CONSTRAINTS
{ double SMbsg,dmunu,csLim;
  printf("\n\n==== Physical Constraints: =====\n"); 
  printf("deltartho=%.2E\n",deltarho());
  printf("gmuon=%.2E\n", gmuon());
  printf("bsgnlo=%.2E ", bsgnlo(&SMbsg)); printf("( SM %.2E )\n",SMbsg);

  printf("bsmumu=%.2E\n", bsmumu());
  printf("btaunu=%.2E\n", btaunu());

  printf("dtaunu=%.2E  ", dtaunu(&dmunu)); printf("dmunu=%.2E\n", dmunu);   
  printf("Rl23=%.3E\n", Rl23());
  if(Zinvisible()) printf("Excluded by Z->invisible\n");
  if(LspNlsp_LEP(&csLim)) printf("Excluded by LEP  by e+,e- -> DM q qbar. Cross section =%.2E [pb] \n",csLim);

  if(masslimits()==0) printf("MassLimits OK\n");
  
  if(blockExists("SPhenoLowEnergy"))
  {
    printf("\n SPheno  low energy observables\n");
    printf("  BR(b -> s gamma)                         %.3E \n", slhaVal("SPhenoLowEnergy",0.,  1,   1   ));  
    printf("  BR(b -> s mu+ mu-)                       %.3E \n", slhaVal("SPhenoLowEnergy",0.,  1,   2   ));  
    printf("  BR(b -> s nu nu)                         %.3E \n", slhaVal("SPhenoLowEnergy",0.,  1,   3   ));  
    printf("  BR(Bd -> e+ e-)                          %.3E \n", slhaVal("SPhenoLowEnergy",0.,  1,   4   ));  
    printf("  BR(Bd -> mu+ mu-)                        %.3E \n", slhaVal("SPhenoLowEnergy",0.,  1,   5   ));  
    printf("  BR(Bd -> tau+ tau-)                      %.3E \n", slhaVal("SPhenoLowEnergy",0.,  1,   6   ));  
    printf("  BR(Bs -> e+ e-)                          %.3E \n", slhaVal("SPhenoLowEnergy",0.,  1,   7   ));  
    printf("  BR(Bs -> mu+ mu-)                        %.3E \n", slhaVal("SPhenoLowEnergy",0.,  1,   8   ));  
    printf("  BR(Bs -> tau+ tau-)                      %.3E \n", slhaVal("SPhenoLowEnergy",0.,  1,   9   ));  
    printf("  BR(B_u -> tau nu)                        %.3E \n", slhaVal("SPhenoLowEnergy",0.,  1,  10   ));  
    printf("  BR(B_u -> tau nu)/BR(B_u -> tau nu)_SM   %.3E \n", slhaVal("SPhenoLowEnergy",0.,  1,  11   ));  
    printf("  |Delta(M_Bd)| [ps^-1]                    %.3E \n", slhaVal("SPhenoLowEnergy",0.,  1,  12   ));  
    printf("  |Delta(M_Bs)| [ps^-1]                    %.3E \n", slhaVal("SPhenoLowEnergy",0.,  1,  13   ));  
    printf("  epsilon_K                                %.3E \n", slhaVal("SPhenoLowEnergy",0.,  1,  16   ));  
    printf("  Delta(M_K)                               %.3E \n", slhaVal("SPhenoLowEnergy",0.,  1,  17   ));  
    printf("  BR(K^0 -> pi^0 nu nu)                    %.3E \n", slhaVal("SPhenoLowEnergy",0.,  1,  18   ));  
    printf("  BR(K^+ -> pi^+ nu nu)                    %.3E \n", slhaVal("SPhenoLowEnergy",0.,  1,  19   ));  
    printf("  Delta(g-2)_electron/2                    %.3E \n", slhaVal("SPhenoLowEnergy",0.,  1,  20   ));  
    printf("  Delta(g-2)_muon/2                        %.3E \n", slhaVal("SPhenoLowEnergy",0.,  1,  21   ));  
    printf("  Delta(g-2)_tau/2                         %.3E \n", slhaVal("SPhenoLowEnergy",0.,  1,  22   ));  
  }
 
}
#endif

#ifdef SUPERISO
{
  int err= callSuperIsoSLHA();
  if(err==0)
  { printf("\nSuperIso Flavour MSSM and ( SM)  observables :\n");
    printf("  BR(b->s gamma)                     %.3E  (%.3E)\n", slhaValFormat("FOBS",0.,"    5    1  %lf    0     2     3    22        "), slhaValFormat("FOBSSM",0.,"    5    1  %lf    0     2     3    22        ")); 
    printf("  Delta0(B->K* gamma)                %.3E  (%.3E)\n", slhaValFormat("FOBS",0.,"  521    4  %lf    0     2   313    22        "), slhaValFormat("FOBSSM",0.,"  521    4  %lf    0     2   313    22        ")); 
    printf("  BR(B_s->mu+ mu-)                   %.3E  (%.3E)\n", slhaValFormat("FOBS",0.,"  531    1  %lf    0     2    13   -13        "), slhaValFormat("FOBSSM",0.,"  531    1  %lf    0     2    13   -13        ")); 
    printf("  BR(B_u->tau nu)                    %.3E  (%.3E)\n", slhaValFormat("FOBS",0.,"  521    1  %lf    0     2   -15    16        "), slhaValFormat("FOBSSM",0.,"  521    1  %lf    0     2   -15    16        ")); 
    printf("  R(B_u->tau nu)                     %.3E  (%.3E)\n", slhaValFormat("FOBS",0.,"  521    2  %lf    0     2   -15    16        "), slhaValFormat("FOBSSM",0.,"  521    2  %lf    0     2   -15    16        ")); 
    printf("  BR(D_s->tau nu)                    %.3E  (%.3E)\n", slhaValFormat("FOBS",0.,"  431    1  %lf    0     2   -15    16        "), slhaValFormat("FOBSSM",0.,"  431    1  %lf    0     2   -15    16        ")); 
    printf("  BR(D_s->mu nu)                     %.3E  (%.3E)\n", slhaValFormat("FOBS",0.,"  431    1  %lf    0     2   -13    14        "), slhaValFormat("FOBSSM",0.,"  431    1  %lf    0     2   -13    14        ")); 
    printf("  BR(B+->D0 tau nu)                  %.3E  (%.3E)\n", slhaValFormat("FOBS",0.,"  521    1  %lf    0     3   421   -15    16  "), slhaValFormat("FOBSSM",0.,"  521    1  %lf    0     3   421   -15    16  ")); 
    printf("  BR(B+->D0 tau nu)/BR(B+-> D0 e nu) %.3E  (%.3E)\n", slhaValFormat("FOBS",0.,"  521   11  %lf    0     3   421   -15    16  "), slhaValFormat("FOBSSM",0.,"  521   11  %lf    0     3   421   -15    16  ")); 
    printf("  BR(K->mu nu)/BR(pi->mu nu)         %.3E  (%.3E)\n", slhaValFormat("FOBS",0.,"  321   11  %lf    0     2   -13    14        "), slhaValFormat("FOBSSM",0.,"  321   11  %lf    0     2   -13    14        ")); 
    printf("  R_mu23                             %.3E  (%.3E)\n", slhaValFormat("FOBS",0.,"  321   12  %lf    0     2   -13    14        "), slhaValFormat("FOBSSM",0.,"  321   12  %lf    0     2   -13    14        ")); 
  }
}
#endif

#ifdef CheckMassMatrix
{
  double MZ,SW,CW,sb,cb;
  MZ=findValW("MZ"); SW=findValW("SW"); CW=findValW("CW"); sb=findValW("sb"); cb=findValW("cb");
  printf("\n    Neutralino Mass Matrix\n");  
  printf(" i,j     Zki*Mk*Zkj    tree level\n");  
  printf(" 1 1    %10.3E     %10.3E\n", findValW("nmm11"), findValW("MG1"));
  printf(" 1 2    %10.3E     %10.3E\n", findValW("nmm12"), 0.);
  printf(" 1 3    %10.3E     %10.3E\n", findValW("nmm13"), -MZ*cb*SW);
  printf(" 1 4    %10.3E     %10.3E\n", findValW("nmm14"),  MZ*sb*SW);
  printf(" 2 2    %10.3E     %10.3E\n", findValW("nmm22"),  findValW("MG2"));
  printf(" 2 3    %10.3E     %10.3E\n", findValW("nmm23"),  MZ*cb*CW);
  printf(" 2 4    %10.3E     %10.3E\n", findValW("nmm24"), -MZ*sb*CW);
  printf(" 3 3    %10.3E     %10.3E\n", findValW("nmm33"),  0. );  
  printf(" 3 4    %10.3E     %10.3E\n", findValW("nmm34"), -findValW("mu"));

}
#endif


#if defined(HIGGSBOUNDS) || defined(HIGGSSIGNALS)
{  int NH0=3, NHch=1; // number of neutral and charged Higgs particles.
   int HB_id[3]={0,0,0},HB_result[3];
   double  HB_obsratio[3],HS_observ=-1,HS_chi2, HS_pval;
   char HB_chan[3][100]={""}, HB_version[50], HS_version[50]; 
   NH0=hbBlocksMO("HB.in",&NHch); 
   system("echo 'BLOCK DMASS\n 25  2  '>> HB.in");
#include "../include/hBandS.inc"
#ifdef HIGGSBOUNDS
   printf("HiggsBounds(%s)\n", HB_version);
   for(int i=0;i<3;i++) if(HB_id[i]) printf("  id= %d  result = %d  obsratio=%.2E  channel= %s \n", HB_id[i],HB_result[i],HB_obsratio[i],HB_chan[i]);
#endif 
#ifdef HIGGSSIGNALS
   if(HS_observ>=0)
   {
     printf("HiggsSignals(%s)\n",HS_version); 
     printf("  Nobservables=%.0f chi^2 = %.2E pval= %.2E\n",HS_observ,HS_chi2, HS_pval);
   }
#endif   
}
#endif

#ifdef LILITH
{  double m2logL, m2logL_reference=0,pvalue;
   int exp_ndf,n_par=0,ndf;
   char  Lilith_version[50];
   if( LilithMO("Lilith_in.xml"))
   {        
#include "../include/Lilith.inc"
      if(ndf)
      {
        printf("LILITH(DB%s):  -2*log(L): %.2f; -2*log(L_reference): %.2f; ndf: %d; p-value: %.2E \n",
        Lilith_version,m2logL,m2logL_reference,ndf,pvalue);
      }  
   } else printf("LILITH: there is no Higgs candidate\n");
}     
#endif


#ifdef SMODELS
{ int status=0, smodelsOK=0; 
  double Rvalue, Rexpected, SmoLsig, SmoLmax, SmoLSM;
  char analysis[50]={},topology[100]={},smodelsInfo[100];
  int LHCrun=LHC8|LHC13;  //  LHC8  - 8TeV; LHC13  - 13TeV;   

  printf("\n\n=====  LHC constraints with SModelS  =====\n\n");

#include "../include/SMODELS.inc" // SLHA interface with SModelS

  printf("SModelS %s \n",smodelsInfo);
  if(smodelsOK) 
  { printf(" highest r-value = %.2E",Rvalue); 
    if(Rvalue>0) 
    { printf(" from %s, topology: %s ",analysis,topology);
      if(Rexpected>0) 
      { printf("\n expected r = %.2E ",Rexpected);
        if(SmoLsig>0) 
        { printf("\n -2log (L_signal, L_max, L_SM) = %.2E %.2E %.2E", 
                  -2*log(SmoLsig),-2*log(SmoLmax),-2*log(SmoLSM)); }
      }
    }  
    if(status==1) { printf("\n excluded by SMS results"); }
    else if(status==0) printf("\n not excluded"); 
    else if(status==-1) printf("\n not not tested by results in SModelS database"); 
    printf("\n");
  } else system("cat smodels.err"); // problem: see smodels.err
}   

#endif 

#ifdef MONOJET
{ double CL=monoJet();
  printf(" Monojet signal exclusion CL is %.3e\n", CL);
}  
#endif




#ifdef OMEGA
{ int fast=1;
  double Beps=1.E-5, cut=0.01;
  double Omega,Xf=25; 
  
// to exclude processes with virtual W/Z in DM   annihilation      
    VZdecay=0; VWdecay=0; cleanDecayTable(); 
    
// to include processes with virtual W/Z  also  in co-annihilation 
//   VZdecay=2; VWdecay=2; cleanDecayTable(); 
    
  printf("\n==== Calculation of relic density =====\n");  

//  sortOddParticles(cdmName);

   Omega=darkOmega(&Xf,fast,Beps,&err);
   printf("Xf=%.2e Omega=%.2e\n",Xf,Omega);

   if(Omega>0)printChannels(Xf,cut,Beps,1,stdout);
/*   
   Omega=darkOmega2(fast,Beps);
   displayPlot("Y","T",Tend,Tstart, 0,2,"Y",0,YF,NULL,"Y1",0,Y1F,NULL);
   printf("Omega2=%e\n",Omega);
*/

// direct access for annihilation channels 

/*
if(omegaCh){
  int i; 
  for(i=0; omegaCh[i].weight>0  ;i++)
  printf(" %.2E %s %s -> %s %s\n", omegaCh[i].weight, omegaCh[i].prtcl[0],
  omegaCh[i].prtcl[1],omegaCh[i].prtcl[2],omegaCh[i].prtcl[3]); 
}  
*/

// to restore default switches  


    VZdecay=1; VWdecay=1; cleanDecayTable();

}
#endif
 

#ifdef INDIRECT_DETECTION
{ 
  int err,i;
  double Emin=1,SMmev=320;/*Energy cut in GeV and solar potential in MV*/
  double  sigmaV;
  char txt[100];
  double SpA[NZ],SpE[NZ],SpP[NZ];
  double FluxA[NZ],FluxE[NZ],FluxP[NZ];
  double SpNe[NZ],SpNm[NZ],SpNl[NZ];  
//  double * SpNe=NULL,*SpNm=NULL,*SpNl=NULL;
  double Etest=Mcdm/2;
 
/* default DarkSUSY parameters */

/*
    K_dif=0.036;
    L_dif=4;  
    Delta_dif=0.6; 
    Vc_dif=10;
    Rdisk=30;
    SMmev=320;
*/                        
  
printf("\n==== Indirect detection =======\n");  

  sigmaV=calcSpectrum(1+2+4,SpA,SpE,SpP,SpNe,SpNm,SpNl ,&err);
    /* Returns sigma*v in cm^3/sec.     SpX - calculated spectra of annihilation.
       Use SpectdNdE(E, SpX) to calculate energy distribution in  1/GeV units.
       
       First parameter 1-includes W/Z polarization
                       2-includes gammas for 2->2+gamma
                       4-print cross sections             
    */


  { 
     double fi=0.1,dfi=M_PI/180.; /* angle of sight and 1/2 of cone angle in [rad] */ 
                                                   /* dfi corresponds to solid angle 1.E-3sr */                                             
     printf("\nPhoton flux  for angle of sight f=%.2f[rad]\n"
     "and spherical region described by cone with angle %.4f[rad]\n",fi,2*dfi);
     gammaFluxTab(fi,dfi, sigmaV, SpA, FluxA);

     printf("Photon flux = %.2E[cm^2 s GeV]^{-1} for E=%.1f[GeV]\n",SpectdNdE(Etest, FluxA), Etest);

#ifdef SHOWPLOTS
     sprintf(txt,"Photon flux for angle of sight %.2f[rad] and cone angle %.2f[rad]",fi,2*dfi);
     displayPlot(txt,"E[GeV]",Emin,Mcdm,0,1,"",0,SpectdNdE,FluxA);
#endif
  }

  { 
    posiFluxTab(Emin, sigmaV, SpE, FluxE);
    if(SMmev>0)  solarModulation(SMmev,0.0005,FluxE,FluxE);
#ifdef SHOWPLOTS     
    displayPlot("positron flux [cm^2 s sr GeV]^{-1}","E[GeV]",Emin,Mcdm,0,1,"",0,SpectdNdE,FluxE);
#endif
    printf("\nPositron flux  =  %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
    SpectdNdE(Etest, FluxE),  Etest); 
  }
  
  {
    pbarFluxTab(Emin, sigmaV, SpP,  FluxP); 
    
    if(SMmev>0)  solarModulation(SMmev,1,FluxP,FluxP);     
#ifdef SHOWPLOTS    
     displayPlot("antiproton flux [cm^2 s sr GeV]^{-1}","E[GeV]",Emin,Mcdm,0,1,"",0,SpectdNdE,FluxP);
#endif
    printf("\nAntiproton flux  =  %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
    SpectdNdE(Etest, FluxP),  Etest);     
  }
}  
#endif

#ifdef LoopGAMMA
{    double vcs_gz,vcs_gg;
     double fi=0.,dfi=M_PI/180.; /* fi angle of sight[rad], dfi  1/2 of cone angle in [rad] */
                                 /* dfi corresponds to solid angle  pi*(1-cos(dfi)) [sr] */                             
     if(loopGamma(&vcs_gz,&vcs_gg)==0)
     {
         printf("\nGamma  ray lines:\n");
         printf("E=%.2E[GeV]  vcs(Z,A)= %.2E[cm^3/s], flux=%.2E[cm^2 s]^{-1}\n",Mcdm-91.19*91.19/4/Mcdm,vcs_gz,
                               gammaFlux(fi,dfi,vcs_gz));  
         printf("E=%.2E[GeV]  vcs(A,A)= %.2E[cm^3/s], flux=%.2E[cm^2 s]^{-1}\n",Mcdm,vcs_gg, 
                             2*gammaFlux(fi,dfi,vcs_gg));
     }

}     
#endif     



#ifdef RESET_FORMFACTORS
{
/* 
   The user has approach to form factors  which specifies quark contents 
   of  proton and nucleon via global parametes like
      <Type>FF<Nucleon><q>
   where <Type> can be "Scalar", "pVector", and "Sigma"; 
         <Nucleon>     "P" or "N" for proton and neutron
         <q>            "d", "u","s"

   calcScalarQuarkFF( Mu/Md, Ms/Md, sigmaPiN[MeV], sigmaS[MeV])  
   calculates and rewrites Scalar form factors
*/

  printf("protonFF (default) d %E, u %E, s %E\n",ScalarFFPd, ScalarFFPu,ScalarFFPs);                               
  printf("neutronFF(default) d %E, u %E, s %E\n",ScalarFFNd, ScalarFFNu,ScalarFFNs);


//                    To restore default form factors of  version 2  call 
     calcScalarQuarkFF(0.553,18.9,55.,243.5);

  printf("protonFF (new)     d %E, u %E, s %E\n",ScalarFFPd, ScalarFFPu,ScalarFFPs);                               
  printf("neutronFF(new)     d %E, u %E, s %E\n",ScalarFFNd, ScalarFFNu,ScalarFFNs);

 
//                    To restore default form factors of  current version  call 
//  calcScalarQuarkFF(0.56,20.2,34,42);

}
#endif

#ifdef CDM_NUCLEON
{ double pA0[2],pA5[2],nA0[2],nA5[2];
  double Nmass=0.939; /*nucleon mass*/
  double SCcoeff;        
  double csSIp,csSIn,csSDp,csSDn;
  int sI,sD; 
printf("\n==== Calculation of CDM-nucleons amplitudes  =====\n");   
#ifdef TEST_Direct_Detection
printf("         TREE LEVEL\n");

    MSSMDDtest(0, pA0,pA5,nA0,nA5);
    printf("Analitic formulae\n");
    printf(" proton:  SI %.3E  SD  %.3E\n",pA0[0],pA5[0]);
    printf(" neutron: SI %.3E  SD  %.3E\n",nA0[0],nA5[0]); 
    nucleonAmplitudes(CDM1, pA0,pA5,nA0,nA5);
    printf("%s-nucleon micrOMEGAs amplitudes:\n",CDM1);
    printf("proton:  SI  %.3E  SD  %.3E\n",pA0[0],pA5[0]);
    printf("neutron: SI  %.3E  SD  %.3E\n",nA0[0],nA5[0]); 

printf("         BOX DIAGRAMS\n");  

    MSSMDDtest(1, pA0,pA5,nA0,nA5);
    printf("Analitic formulae\n");
    printf(" proton:  SI %.3E  SD  %.3E\n",pA0[0],pA5[0]);
    printf(" neutron: SI %.3E  SD  %.3E\n",nA0[0],nA5[0]); 
    
#endif

    nucleonAmplitudes(CDM1,pA0,pA5,nA0,nA5);
    printf("%s-nucleon micrOMEGAs amplitudes:\n",CDM1);
    printf("proton:  SI  %.3E  SD  %.3E\n",pA0[0],pA5[0]);
    printf("neutron: SI  %.3E  SD  %.3E\n",nA0[0],nA5[0]); 

    SCcoeff=4/M_PI*3.8937966E8*pow(Nmass*Mcdm/(Nmass+ Mcdm),2.);
    csSIp=  SCcoeff*pA0[0]*pA0[0];  
    csSDp=3*SCcoeff*pA5[0]*pA5[0];  
    csSIn=  SCcoeff*nA0[0]*nA0[0];  
    csSDn=3*SCcoeff*nA5[0]*nA5[0];  
                    
    printf("\n==== %s-nucleon cross sections[pb] ====\n",CDM1);
    printf(" proton  SI %.3E  SD %.3E\n",csSIp,csSDp);
    printf(" neutron SI %.3E  SD %.3E\n",csSIn,csSDn);

    if(pA0[0]*nA0[0]<0) sI=-1; else sI=1;
    if(pA5[0]*nA5[0]<0) sD=-1; else sD=1;
    char*expName="";
    double pval=DD_pvalCS(AllDDexp, Maxwell, csSIp, sI*csSIn,csSDp,sD*csSDn, &expName);
    if(pval<0.1 )  printf("Excluded by %s [CDM_NUCLEON] %.1f%% \n", expName, 100*(1-pval)); 
    else printf("Not excluded by DD experiments  at 90%% level \n");          
}
#endif
  
#ifdef CDM_NUCLEUS
{ char* expName; 
  printf("\n===== Direct detection exclusion:======\n");
  double pval=DD_pval(AllDDexp, Maxwell, &expName);
       if(pval<0.1 )  printf("Excluded by %s [CDM_NUCLEUS]  %.1f%%\n", expName, 100*(1-pval)); 
  else printf("Not excluded by DD experiments  at 90%% level \n");  
}
#endif 

#ifdef NEUTRINO
{ double nu[NZ], nu_bar[NZ],mu[NZ];
  int forSun=1;
  double Emin=1;

WIMPSIM=0;
 
  printf("\n===============Neutrino Telescope=======  for  "); 
  if(forSun) printf("Sun\n"); else printf("Earth\n");  

  err=neutrinoFlux(Maxwell,forSun, nu,nu_bar);
#ifdef SHOWPLOTS
  displayPlot("neutrino fluxes [1/Year/km^2/GeV]","E[GeV]",Emin,Mcdm,0, 2,"dnu/dE",0,SpectdNdE,nu,"dnu_bar/dE",0,SpectdNdE,nu_bar);
#endif

printf(" E>%.1E GeV neutrino/anti-neutrino fluxes   %.2E/%.2E [1/Year/km^2]\n",Emin,
          spectrInfo(Emin,nu,NULL), spectrInfo(Emin,nu_bar,NULL));  
//  ICE CUBE
if(forSun)printf("IceCube22 exclusion confidence level = %.2E%%\n", 100*exLevIC22(nu,nu_bar,NULL));
  
/* Upward events */
  
  muonUpward(nu,nu_bar, mu);

#ifdef SHOWPLOTS
  displayPlot("Upward muons[1/Year/km^2/GeV]","E",Emin,Mcdm/2, 0,1,"mu",0,SpectdNdE,mu);  
#endif

  printf(" E>%.1E GeV Upward muon flux    %.2E [1/Year/km^2]\n",Emin,spectrInfo(Emin,mu,NULL));
  
/* Contained events */
  muonContained(nu,nu_bar,1., mu);
#ifdef SHOWPLOTS  
  displayPlot("Contained  muons[1/Year/km^3/GeV]","E",Emin,Mcdm,0,1,"",0,SpectdNdE,mu); 
#endif
  printf(" E>%.1E GeV Contained muon flux %.2E [1/Year/km^3]\n",Emin,spectrInfo(Emin,mu,NULL)); 
}        
#endif 


#ifdef DECAYS
{  
  txtList L;
   double width,br;
   char * pname;
   printf("\n================= Decays ==============\n");

   pname = "h";
   width=pWidth(pname,&L);
   printf("\n%s :   total width=%.2E \n",pname,width);
//   printTxtList(L,stdout);  
   printPartialWidth(width,L,stdout);

   pname = "H3";
   width=pWidth(pname,&L);
   printf("\n%s :   total width=%.2E \n",pname,width);
//   printTxtList(L,stdout);          
   printPartialWidth(width,L,stdout);

   pname = "H";
   width=pWidth(pname,&L);
   printf("\n%s :   total width=%.2E \n",pname,width);
//   printTxtList(L,stdout);   
   printPartialWidth(width,L,stdout);

   pname = "~1+";
   width=pWidth(pname,&L);
   printf("\n%s :   total width=%.2E \n",pname,width);
//   printTxtList(L,stdout);   
   printPartialWidth(width,L,stdout);

   printf("Example of 1->3 decay:\n"); 
   numout*cc=newProcess("~o2->~o1,e,E");
   int err;
   printf("width(~o2->~o1,e,E)=%e\n", pWidthCC(cc,&err));
   
}
#endif


#ifdef CROSS_SECTIONS
{
  char* next,next_;
  double nextM;
    
  next=nextOdd(1,&nextM); 
  if(next && nextM<1000)  
  { 
     double cs, Pcm=6500, Qren, Qfact, pTmin=0;
     int nf=3;
     char*next_=antiParticle(next);
     Qren=Qfact=nextM; 
 
     printf("\npp > nextOdd  at sqrt(s)=%.2E GeV\n",2*Pcm);  
  
     Qren=Qfact;
     cs=hCollider(Pcm,1,nf,Qren, Qfact, next,next_,pTmin,1);
     printf("Production of 'next' odd particle: cs(pp-> %s,%s)=%.2E[pb]\n",next,next_, cs);
  }  
}
#endif

#ifdef CLEAN
  system("rm -f suspect2_lha.in suspect2_lha.out suspect2.out");
  system("rm -f LesHouches.in Messages.out SPheno.spc");
  system("rm -f LesHin LesHout");
  system("rm -f  nngg.*  output.flha ");
//  system("rm -f HB.* HS.* hb.* hs.*  debug_channels.txt debug_predratio.txt  Key.dat");
  system("rm -f Lilith_*   particles.py*");
//  system("rm -f  smodels.*");  
#endif 

  killPlots();
  return 0;
}