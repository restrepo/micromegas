/*====== Modules ===============
   Keys to switch on
   various modules of micrOMEGAs
================================*/

#define MASSES_INFO
  /* Display information about mass spectrum  */

#define CONSTRAINTS

//#define MONOJET
//#define HIGGSBOUNDS
//#define HIGGSSIGNALS
//#define LILITH
//#define SMODELS

//#define OMEGA  /* Calculate relic density and display contribution of  individual channels */
#define FREEZEIN      /*  Calculate relic density in Freeze-in scenario  */


//#define INDIRECT_DETECTION
  /* Compute spectra of gamma/positron/antiprotons/neutrinos for DM annihilation;
     Calculate <sigma*v>;
     Integrate gamma signal over DM galactic squared density for given line
     of sight;
     Calculate galactic propagation of positrons and antiprotons.
  */

//#define RESET_FORMFACTORS
  /* Modify default nucleus form factors,
    DM velocity distribution,
    A-dependence of Fermi-dencity
  */
//#define CDM_NUCLEON
  /* Calculate amplitudes and cross-sections for  CDM-mucleon collisions */

//#define CDM_NUCLEUS
  /* Calculate number of events for 1kg*day and recoil energy distibution
      for various nuclei
  */
//#define NEUTRINO
 /*  Neutrino signal of DM annihilation in Sun and Earth */

#define DECAYS

//#define CROSS_SECTIONS

/*===== end of Modules  ======*/

/*===== Options ========*/
/*#define SHOWPLOTS*/
     /* Display  graphical plots on the screen */
#define CLEAN
/*===== End of DEFINE  settings ===== */


#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"
#include"lib/pmodel.h"

int main(int argc,char** argv)
{  int err;
   char cdmName[10];
   int spin2, charge3,cdim;

  ForceUG=1;  /* to Force Unitary Gauge assign 1 */

  VZdecay=0; VWdecay=0;

  if(argc==1)
  {
      printf(" Correct usage:  ./main  <file with parameters> \n");
      printf("Example: ./main data1.par\n");
      exit(1);
  }

  err=readVar(argv[1]);

  if(err==-1)     {printf("Can not open the file\n"); exit(1);}
  else if(err>0)  { printf("Wrong file contents at line %d\n",err);exit(1);}



  err=sortOddParticles(cdmName);
  if(err) { printf("Can't calculate %s\n",cdmName); return 1;}





/*
  if(CDM1)
  {
     qNumbers(CDM1, &spin2, &charge3, &cdim);
     printf("\nDark matter candidate is '%s' with spin=%d/2 mass=%.2E\n",CDM1,  spin2,Mcdm1);
     if(charge3) printf("Dark Matter has electric charge %d/3\n",charge3);
     if(cdim!=1) printf("Dark Matter is a color particle\n");
  }
  if(CDM2)
  {
     qNumbers(CDM2, &spin2, &charge3, &cdim);
     printf("\nDark matter candidate is '%s' with spin=%d/2 mass=%.2E\n",CDM2,spin2,Mcdm2);
     if(charge3) printf("Dark Matter has electric charge %d/3\n",charge3);
     if(cdim!=1) printf("Dark Matter is a color particle\n");
  }
*/

#ifdef MASSES_INFO
{
  printf("\n=== MASSES OF HIGGS AND ODD PARTICLES: ===\n");
  printHiggs(stdout);
  printMasses(stdout,1);
}
#endif

#ifdef CONSTRAINTS
{ double csLim;
  if(Zinvisible()) printf("Excluded by Z->invizible\n");
  if(LspNlsp_LEP(&csLim)) printf("LEP excluded by e+,e- -> DM q q-\\bar  Cross Section= %.2E pb\n",csLim);
}
#endif

#ifdef MONOJET
{ double CL=monoJet();
  printf(" Monojet signal exclusion CL is %.3e\n", CL);
}
#endif

#if defined(HIGGSBOUNDS) || defined(HIGGSSIGNALS)
{  int NH0,NHch;  // number of neutral and charged Higgs particles.
   double HB_result,HB_obsratio,HS_observ,HS_chi2, HS_pval;
   char HB_chan[100]={""}, HB_version[50], HS_version[50];
   NH0=hbBlocksMO("HB.in",&NHch);
   system("echo 'BLOCK DMASS\n 25  2  '>> HB.in");
#include "../include/hBandS.inc"
#ifdef HIGGSBOUNDS
   printf("HB(%s): result=%.0f  obsratio=%.2E  channel= %s \n", HB_version,HB_result,HB_obsratio,HB_chan);
#endif
#ifdef HIGGSSIGNALS
   printf("HS(%s): Nobservables=%.0f chi^2 = %.2E pval= %.2E\n",HS_version,HS_observ,HS_chi2, HS_pval);
#endif
}
#endif

#ifdef LILITH
{  double m2logL, m2logL_reference=0,pvalue;
   int exp_ndf,n_par=0,ndf;
   char call_lilith[100], Lilith_version[20];
   if(LilithMO("Lilith_in.xml"))
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
{  int result=0;
   double Rvalue=0;
   char analysis[30]={},topology[30]={};
   int LHCrun=LHC8|LHC13;  //  LHC8  - 8TeV; LHC13  - 13TeV; 
#include "../include/SMODELS.inc"
}
#endif

#ifdef OMEGA
{ int fast=1;
  double Beps=1.E-4, cut=0.01;
  double Omega;  
  int i,err; 
  printf("\n==== Calculation of relic density =====\n");   

  if(CDM1 && CDM2) 
  {
  
    Omega= darkOmega2(fast,Beps);
    printf("Omega_1h^2=%.2E\n", Omega*(1-fracCDM2));
    printf("Omega_2h^2=%.2E\n", Omega*fracCDM2);
  } else
  {  double Xf;
     Omega=darkOmega(&Xf,fast,Beps,&err);
     printf("Xf=%.2e Omega=%.2e\n",Xf,Omega);
     if(Omega>0)printChannels(Xf,cut,Beps,1,stdout);
  }
}

#endif



#ifdef FREEZEIN
{ 
  double omegaFi;
  int err;

  printf("\n==== Calculation of relic density =====\n");
  double TR=1E6;
  toFeebleList("~x1");
  omegaFi=darkOmegaFi(TR,&err);

  if(err) printf("darkOmegaFi error code = %d\n", err);
  else  printf("Freeze-in Omega h^2 =%.3E\n", omegaFi);
  printChannelsFi(0,0,stdout);
}

#endif



#ifdef INDIRECT_DETECTION
{
  int err,i;
  double Emin=1,/* Energy cut  in GeV   */  sigmaV;
  double vcs_gz,vcs_gg;
  char txt[100];
  double SpA[NZ],SpE[NZ],SpP[NZ];
  double FluxA[NZ],FluxE[NZ],FluxP[NZ];
  double * SpNe=NULL,*SpNm=NULL,*SpNl=NULL;
  double Etest=Mcdm/2;

printf("\n==== Indirect detection =======\n");

  sigmaV=calcSpectrum(1+2+4,SpA,SpE,SpP,SpNe,SpNm,SpNl ,&err);
    /* Returns sigma*v in cm^3/sec.     SpX - calculated spectra of annihilation.
       Use SpectdNdE(E, SpX) to calculate energy distribution in  1/GeV units.

       First parameter 1-includes W/Z polarization
                       2-includes gammas for 2->2+gamma
                       4-print cross sections
    */



  {
     double fi=0.1,dfi=0.05; /* angle of sight and 1/2 of cone angle in [rad] */

     gammaFluxTab(fi,dfi, sigmaV, SpA,  FluxA);
     printf("Photon flux  for angle of sight f=%.2f[rad]\n"
     "and spherical region described by cone with angle %.2f[rad]\n",fi,2*dfi);
#ifdef SHOWPLOTS
     sprintf(txt,"Photon flux[cm^2 s GeV]^{-1} at f=%.2f[rad], cone angle %.2f[rad]",fi,2*dfi);
     displaySpectrum(txt,Emin,Mcdm,FluxA);
#endif
     printf("Photon flux = %.2E[cm^2 s GeV]^{-1} for E=%.1f[GeV]\n",SpectdNdE(Etest, FluxA), Etest);
  }

  {
    posiFluxTab(Emin, sigmaV, SpE,  FluxE);
#ifdef SHOWPLOTS
    displaySpectrum("positron flux [cm^2 s sr GeV]^{-1}" ,Emin,Mcdm,FluxE);
#endif
    printf("Positron flux  =  %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
    SpectdNdE(Etest, FluxE),  Etest);
  }

  {
    pbarFluxTab(Emin, sigmaV, SpP,  FluxP  );
#ifdef SHOWPLOTS
     displaySpectrum("antiproton flux [cm^2 s sr GeV]^{-1}" ,Emin, Mcdm,FluxP);
#endif
    printf("Antiproton flux  =  %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
    SpectdNdE(Etest, FluxP),  Etest);
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

   calcScalarFF( Mu/Md, Ms/Md, sigmaPiN[MeV], sigma0[MeV])
   calculates and rewrites Scalar form factors
*/
  printf("\n======== RESET_FORMFACTORS ======\n");

  printf("protonFF (default) d %.2E, u %.2E, s %.2E\n",ScalarFFPd, ScalarFFPu,ScalarFFPs);
  printf("neutronFF(default) d %.2E, u %.2E, s %.2E\n",ScalarFFNd, ScalarFFNu,ScalarFFNs);

  calcScalarFF(0.553,18.9,70.,35.);

  printf("protonFF (new)     d %.2E, u %.2E, s %.2E\n",ScalarFFPd, ScalarFFPu,ScalarFFPs);
  printf("neutronFF(new)     d %.2E, u %.2E, s %.2E\n",ScalarFFNd, ScalarFFNu,ScalarFFNs);

}
#endif

#ifdef CDM_NUCLEON
{ double pA0[2],pA5[2],nA0[2],nA5[2];
  double Nmass=0.939; /*nucleon mass*/
  double SCcoeff;

printf("\n==== Calculation of CDM-nucleons amplitudes  =====\n");

  if(CDM1)
  {
    nucleonAmplitudes(CDM1, pA0,pA5,nA0,nA5);
    printf("CDM[antiCDM]-nucleon micrOMEGAs amplitudes for %s \n",CDM1);
    printf("proton:  SI  %.3E [%.3E]  SD  %.3E [%.3E]\n",pA0[0], pA0[1],  pA5[0], pA5[1] );
    printf("neutron: SI  %.3E [%.3E]  SD  %.3E [%.3E]\n",nA0[0], nA0[1],  nA5[0], nA5[1] );

  SCcoeff=4/M_PI*3.8937966E8*pow(Nmass*Mcdm/(Nmass+ Mcdm),2.);
    printf("CDM[antiCDM]-nucleon cross sections[pb]:\n");
    printf(" proton  SI %.3E [%.3E] SD %.3E [%.3E]\n",
       SCcoeff*pA0[0]*pA0[0],SCcoeff*pA0[1]*pA0[1],3*SCcoeff*pA5[0]*pA5[0],3*SCcoeff*pA5[1]*pA5[1]);
    printf(" neutron SI %.3E [%.3E] SD %.3E [%.3E]\n",
       SCcoeff*nA0[0]*nA0[0],SCcoeff*nA0[1]*nA0[1],3*SCcoeff*nA5[0]*nA5[0],3*SCcoeff*nA5[1]*nA5[1]);
  }
  if(CDM2)
  {
    nucleonAmplitudes(CDM2, pA0,pA5,nA0,nA5);
    printf("CDM[antiCDM]-nucleon micrOMEGAs amplitudes for %s \n",CDM2);
    printf("proton:  SI  %.3E [%.3E]  SD  %.3E [%.3E]\n",pA0[0], pA0[1],  pA5[0], pA5[1] );
    printf("neutron: SI  %.3E [%.3E]  SD  %.3E [%.3E]\n",nA0[0], nA0[1],  nA5[0], nA5[1] );

  SCcoeff=4/M_PI*3.8937966E8*pow(Nmass*Mcdm/(Nmass+ Mcdm),2.);
    printf("CDM[antiCDM]-nucleon cross sections[pb]:\n");
    printf(" proton  SI %.3E [%.3E] SD %.3E [%.3E]\n",
       SCcoeff*pA0[0]*pA0[0],SCcoeff*pA0[1]*pA0[1],3*SCcoeff*pA5[0]*pA5[0],3*SCcoeff*pA5[1]*pA5[1]);
    printf(" neutron SI %.3E [%.3E] SD %.3E [%.3E]\n",
       SCcoeff*nA0[0]*nA0[0],SCcoeff*nA0[1]*nA0[1],3*SCcoeff*nA5[0]*nA5[0],3*SCcoeff*nA5[1]*nA5[1]);
  }
}
#endif

#ifdef CDM_NUCLEUS
{ double dNdE[300];
  double nEvents;

printf("\n======== Direct Detection ========\n");

  nEvents=nucleusRecoil(Maxwell,73,Z_Ge,J_Ge73,SxxGe73,NULL,dNdE);

  printf("73Ge: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));

#ifdef SHOWPLOTS
    displayRecoilPlot(dNdE,"Distribution of recoil energy of 73Ge",0,199);
#endif

  nEvents=nucleusRecoil(Maxwell,131,Z_Xe,J_Xe131,SxxXe131,NULL,dNdE);

  printf("131Xe: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));
#ifdef SHOWPLOTS
    displayRecoilPlot(dNdE,"Distribution of recoil energy of 131Xe",0,199);
#endif

  nEvents=nucleusRecoil(Maxwell,23,Z_Na,J_Na23,SxxNa23,NULL,dNdE);

  printf("23Na: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));
#ifdef SHOWPLOTS
    displayRecoilPlot(dNdE,"Distribution of recoil energy of 23Na",0,199);
#endif

  nEvents=nucleusRecoil(Maxwell,127,Z_I,J_I127,SxxI127,NULL,dNdE);

  printf("I127: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));
#ifdef SHOWPLOTS
    displayRecoilPlot(dNdE,"Distribution of recoil energy of 127I",0,199);
#endif

}
#endif

#ifdef NEUTRINO
if(!CDM1 || !CDM2)
{ double nu[NZ], nu_bar[NZ],mu[NZ];
  double Ntot;
  int forSun=1;
  double Emin=1;

 printf("\n===============Neutrino Telescope=======  for  ");
 if(forSun) printf("Sun\n"); else printf("Earth\n");

  err=neutrinoFlux(Maxwell,forSun, nu,nu_bar);
#ifdef SHOWPLOTS

printf(  "nu[0]=%E nu_bar[0]=%E\n", nu[0], nu_bar[0]);

  displaySpectrum("nu flux from Sun [1/Year/km^2/GeV]",Emin,Mcdm,nu);
  displaySpectrum("nu-bar from Sun [1/Year/km^2/GeV]",Emin,Mcdm,nu_bar);
#endif
{
    printf(" E>%.1E GeV neutrino flux       %.2E [1/Year/km^2] \n",Emin,spectrInfo(Emin,nu,NULL));
    printf(" E>%.1E GeV anti-neutrino flux  %.2E [1/Year/km^2]\n",Emin,spectrInfo(Emin,nu_bar,NULL));
}

/* Upward events */

  muonUpward(nu,nu_bar, mu);
#ifdef SHOWPLOTS
  displaySpectrum("Upward muons[1/Year/km^2/GeV]",1,Mcdm/2,mu);
#endif
    printf(" E>%.1E GeV Upward muon flux    %.2E [1/Year/km^2]\n",Emin,spectrInfo(Emin,mu,NULL));

/* Contained events */
  muonContained(nu,nu_bar,1., mu);
#ifdef SHOWPLOTS
  displaySpectrum("Contained  muons[1/Year/km^3/GeV]",Emin,Mcdm,mu);
#endif
  printf(" E>%.1E GeV Contained muon flux %.2E [1/Year/km^3]\n",Emin,spectrInfo(Emin/Mcdm,mu,NULL));
}
#endif


#ifdef DECAYS
{ char*  pname = pdg2name(25);
  txtList L;
  double width;
  if(pname)
  {
    width=pWidth(pname,&L);
    printf("\n%s :   total width=%E \n and Branchings:\n",pname,width);
    printTxtList(L,stdout);
  }

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
  system("rm -f HB.* HB.* hb.* hs.*  debug_channels.txt debug_predratio.txt  Key.dat");
  system("rm -f Lilith_*   particles.py*");
  system("rm -f   smodels.in  smodels.log  smodels.out  summary.*");
#endif



  killPlots();
  return 0;
}
