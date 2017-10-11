/*====== Modules ===============
   Keys to switch on 
   various modules of micrOMEGAs  
================================*/

#define MASSES_INFO      
  /* Display information about mass spectrum  */

#define CONSTRAINTS     
      /* Display  deltarho, B_>sgamma, Bs->mumu, gmuon, Z1->inv,
         check LEP mass limits and Zprime limits
      */ 

#define OMEGA            
  /* Calculate relic density and display contribution of  individual channels */

#define INDIRECT_DETECTION  
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

//#define HIGGSBOUNDS  "../Packages/HiggsBounds-4.2.1"
//#define HIGGSSIGNALS "../Packages/HiggsSignals-1.4.0"
  
#define LILITH "../Packages/Lilith-1.1.2"

//#define SMODELS "../Packages/smodels-v1.0.3-micromegas"

#define CDM_NUCLEON     
  /* Calculate amplitudes and cross-sections for  CDM-mucleon collisions */  

//#define CDM_NUCLEUS
  /* Calculate number of events for 1kg*day and recoil energy distibution
      for various nuclei
  */
  
//#define NEUTRINO
 /*  Neutrino signal of DM annihilation in Sun and Earth */
 
#define DECAYS 
      /* Calculate decay widths and branchings  */

//#define CROSS_SECTIONS */
      /* Calculate cross sections and widths for 
         reactions specified by the user
      */

#define CLEAN    to clean intermediate files
  
/*===== end of Modules  ======*/

/*===== Options ========*/
/*#define SHOWPLOTS*/
     /* Display  graphical plots on the screen */ 

/*===== End of DEFINE  settings ===== */


#include"../sources/micromegas.h"
#include"../sources/micromegas_aux.h"
#include"lib/pmodel.h"

#define min(a,b) (a<=b?a:b)

int main(int argc,char** argv)
{  int err;
   char cdmName[10];
   int spin2, charge3,cdim;
   double Omega,bsg,DMd,DMs,bsmumu,btaunu,Damu,bxismu1,bxismu2,Delrho;

  ForceUG=1;  /* to Force Unitary Gauge assign 1 */

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

   
  err=UMSSMTools();      
  slhaWarnings(stdout);
   
#ifdef MASSES_INFO
{
  printf("\n=== MASSES OF HIGGS AND ODD PARTICLES: ===\n");
  printHiggs(stdout);
  printMasses(stdout,1);
}
#endif

  
  if(err) { printf("Can't calculate %s\n",cdmName); return 1;}
  
   qNumbers(cdmName, &spin2, &charge3, &cdim);
   printf("\nDark matter candidate is '%s' with spin=%d/2\n",
    cdmName,       spin2); 
   if(charge3) { printf("Dark Matter has electric charge %d/3\n",charge3); exit(1);}
   if(cdim!=1) { printf("Dark Matter is a color particle\n"); exit(1);}
   txtList L;
   double MNLSP,NLSPMass;
   qNumbers(nextOdd(1, &NLSPMass), &spin2, &charge3, &cdim);
   MNLSP=pMass(nextOdd(1, &NLSPMass));
   printf("NLSP is '%s' with spin=%d/2, electric charge %d/3 and the mass difference with the LSP is %f GeV\n",nextOdd(1, &NLSPMass),spin2,charge3,MNLSP-Mcdm);
   if(cdim!=1)printf("NLSP is a color particle\n");
   if (pWidth(nextOdd(1, &NLSPMass),&L)==0){printf("WARNING: the NLSP '%s' is stable\n",nextOdd(1, &NLSPMass));}


   if(strcmp(cdmName,"~o1")) printf("~o1 is not CDM\n"); 
   else o1Contents(stdout);
  
//Get the corrected Higgs branching ratios from UMSSMTools :
slhaRead("UMSSM_decay.dat",0);  
  
#ifdef OMEGA
{
  int fast=1;
  double Beps=1.E-5, cut=0.01;
  double Xf;   
  
// to include processes with virtual W/Z in DM   annihilation      
   VZdecay=0; VWdecay=0; cleanDecayTable(); 

// to include processes with virtual W/Z  also  in co-annihilation 
//   VZdecay=2; VWdecay=2; cleanDecayTable(); 
  
  printf("\n==== Calculation of relic density =====\n"); 
  sortOddParticles(cdmName);
  Omega=darkOmega(&Xf,fast,Beps);
  printf("Xf=%.2e Omega=%.2e\n",Xf,Omega);
  printChannels(Xf,cut,Beps,1,stdout);

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

#ifdef CONSTRAINTS
{
slhaRead("UMSSM_spectr.dat",0);
printf("\n\n==================================\n");
printf("==== Low energy observables: =====\n");
printf("==================================\n");


bsg=slhaVal("LOWEN",0.,1,1);
printf("\n BR(B -> Xs gamma) = %.3e",slhaVal("LOWEN",0.,1,1));
DMd=slhaVal("LOWEN",0.,1,2);
printf("\n Delta M_d = %.3e ps^-1",slhaVal("LOWEN",0.,1,2));
DMs=slhaVal("LOWEN",0.,1,3);
printf("\n Delta M_s = %.3e ps^-1",slhaVal("LOWEN",0.,1,3));
bsmumu=slhaVal("LOWEN",0.,1,4);
printf("\n BR(Bs->mu+mu-) = %.3e",slhaVal("LOWEN",0.,1,4));
btaunu=slhaVal("LOWEN",0.,1,5);
printf("\n BR(B+ -> tau+ + nu_tau) = %.3e",slhaVal("LOWEN",0.,1,5));
Damu=slhaVal("LOWEN",0.,1,6);
printf("\n BSM contr. to Delta a_mu = %.3e",slhaVal("LOWEN",0.,1,6));
bxismu1=slhaVal("LOWEN",0.,1,7);
printf("\n BR(B-->X_s mu+ mu) for low M_{l+l-}^2 = %.3e",slhaVal("LOWEN",0.,1,7));
bxismu2=slhaVal("LOWEN",0.,1,8);
printf("\n BR(B-->X_s mu+ mu) for high M_{l+l-}^2 = %.3e",slhaVal("LOWEN",0.,1,8));
printf("\n deltartho = %.3e\n",deltarho());
Delrho=deltarho();
if(Zinv()==0) {printf("Limit on the invisible Z1-width OK\n");}
if(masslimits()==0) printf("LEP limits OK\n");
if(Zprimelimits()==0) {printf("LHC limits on new Zprime OK\n");}
}
#endif


#ifdef HIGGSBOUNDS

printf("\n\n=======================\n");
printf("==== HIGGSBOUNDS: =====\n");
printf("=======================\n");

   if(access(HIGGSBOUNDS "/HiggsBounds",X_OK )) system( "cd " HIGGSBOUNDS "; ./configure; make ");
   system("cat UMSSM_spectr.dat UMSSM_decay.dat > HB.in");
   HBblocks("HB.in");
  System("%s/HiggsBounds  LandH SLHA 4 1 HB.in HB.out >hb.stdout",HIGGSBOUNDS);
  slhaRead("HB.out",1+4);
   printf("HB result= %.0E  obsratio=%.2E\n",slhaVal("HiggsBoundsResults",0.,2,1,2), slhaVal("HiggsBoundsResults",0.,2,1,3)  );
   { char hbInfo[100];
    if(0==slhaSTRFormat("HiggsBoundsResults","1 5 ||%[^|]||",hbInfo)) printf("Channel: %s\n",hbInfo);
   } 

#endif



#ifdef HIGGSSIGNALS
#define DataSet " latestresults "
#define Method  " peak " 
//#define  Method " mass "
//#define  Method " both "
#define PDF  " 2 "  // Gaussian
//#define PDF " 1 "  // box 
//#define PDF " 3 "  // box+Gaussia
#define dMh " 2 "
printf("\n\n========================\n");
printf("==== HIGGSSIGNALS: =====\n");
printf("========================\n");

   if(access(HIGGSSIGNALS "/HiggsSignals",X_OK )) system( "cd " HIGGSSIGNALS "; ./configure; make ");
     system("rm -f HS.in HS.out");
     system("cat UMSSM_spectr.dat UMSSM_decay.dat > HS.in");
     HBblocks("HS.in");
     system("echo 'BLOCK DMASS\n 25 " dMh " '>> HS.in");
     System(HIGGSSIGNALS "/HiggsSignals" DataSet Method  PDF " SLHA 4 1 HS.in > hs.stdout");
     System("grep -A 10000  HiggsSignalsResults HS.in > HS.out");
     slhaRead("HS.out",1+4);
     printf("  Number of observables %.0f\n",slhaVal("HiggsSignalsResults",0.,1,7));
     printf("  total chi^2= %.1E\n",slhaVal("HiggsSignalsResults",0.,1,12));
     printf("  HS p-value = %.1E\n", slhaVal("HiggsSignalsResults",0.,1,13));
#undef dMh
#undef PDF
#undef Method
#undef DataSet
     
#endif

#ifdef LILITH
   if(LiLithF("Lilith_in.xml"))
   {  double  like; 
      int exp_ndf;
      system("python " LILITH "/run_lilith.py  Lilith_in.xml  -s -r  Lilith_out.slha");
      slhaRead("Lilith_out.slha", 1);
      like = slhaVal("LilithResults",0.,1,0);
      exp_ndf = slhaVal("LilithResults",0.,1,1);
      printf("LILITH:  -2*log(L): %f; exp ndf: %d \n", like,exp_ndf );
   } else printf("LILITH: there is no Higgs candidate\n");
     
#endif

#ifdef SMODELS
{  int res;

   smodels(4000.,5, 0.1, "smodels.in",0);
   system("make -C " SMODELS); 
   system(SMODELS "/runTools.py xseccomputer -p -N -O -f smodels.in");
   system(SMODELS "/runSModelS.py -f smodels.in -s smodels.res -particles ./  > smodels.out "); 
   slhaRead("smodels.res", 1);
   res=slhaVal("SModelS_Exclusion",0.,2,0,0); 
   switch(res)
   { case -1: printf("SMODELS: no channels for testing\n");break;
     case  0: printf("SMODELS: not excluded\n");break; 
     case  1:  printf("SMODELS: excluded\n");break;
   }  
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
  printf("sigmav=%.2E[cm^3/s]\n",sigmaV);  


  if(SpA)
  { 
     double fi=0.1,dfi=0.05; /* angle of sight and 1/2 of cone angle in [rad] */ 

     gammaFluxTab(fi,dfi, sigmaV, SpA,  FluxA);     
     printf("Photon flux  for angle of sight f=%.2f[rad]\n"
     "and spherical region described by cone with angle %.2f[rad]\n",fi,2*dfi);
#ifdef SHOWPLOTS
     sprintf(txt,"Photon flux[cm^2 s GeV]^{1} at f=%.2f[rad], cone angle %.2f[rad]",fi,2*dfi);
     displaySpectrum(FluxA,txt,Emin,Mcdm,1);
#endif
     printf("Photon flux = %.2E[cm^2 s GeV]^{-1} for E=%.1f[GeV]\n",SpectdNdE(Etest, FluxA), Etest);       
  }

  if(SpE)
  { 
    posiFluxTab(Emin, sigmaV, SpE,  FluxE);
#ifdef SHOWPLOTS     
    displaySpectrum(FluxE,"positron flux [cm^2 s sr GeV]^{-1}" ,Emin,Mcdm, 1);
#endif
    printf("Positron flux  =  %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
    SpectdNdE(Etest, FluxE),  Etest);           
  }
  
  if(SpP)
  { 
    pbarFluxTab(Emin, sigmaV, SpP,  FluxP  ); 
#ifdef SHOWPLOTS    
     displaySpectrum(FluxP,"antiproton flux [cm^2 s sr GeV]^{-1}" ,Emin, Mcdm,1);
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

  calcScalarQuarkFF(0.553,18.9,45.5,26.);

  printf("protonFF (new)     d %.2E, u %.2E, s %.2E\n",ScalarFFPd, ScalarFFPu,ScalarFFPs);                               
  printf("neutronFF(new)     d %.2E, u %.2E, s %.2E\n",ScalarFFNd, ScalarFFNu,ScalarFFNs);


/* Option to change parameters of DM velocity  distribution  */   
//   SetfMaxwell(220.,600.);
/* 
    dN  ~  exp(-v^2/arg1^2)*Theta(v-arg2)  d^3v     
    Earth velocity with respect to Galaxy defined by 'Vearth' parameter.
    All parameters are  in [km/s] units.       
*/


}
#endif

#ifdef CDM_NUCLEON
{ double pA0[2],pA5[2],nA0[2],nA5[2];
  double Nmass=0.939; /*nucleon mass*/
  double SCcoeff;        

printf("\n==== Calculation of CDM-nucleons amplitudes  =====\n");   

    nucleonAmplitudes(CDM1,pA0,pA5,nA0,nA5);
    printf("CDM[antiCDM]-nucleon micrOMEGAs amplitudes:\n");
    printf("proton:  SI  %.3E [%.3E]  SD  %.3E [%.3E]\n",pA0[0], pA0[1],  pA5[0], pA5[1] );
    printf("neutron: SI  %.3E [%.3E]  SD  %.3E [%.3E]\n",nA0[0], nA0[1],  nA5[0], nA5[1] ); 

  SCcoeff=4/M_PI*3.8937966E8*pow(Nmass*Mcdm/(Nmass+ Mcdm),2.);
    printf("CDM[antiCDM]-nucleon cross sections[pb]:\n");
    printf(" proton  SI %.3E [%.3E] SD %.3E [%.3E]\n",
       SCcoeff*pA0[0]*pA0[0],SCcoeff*pA0[1]*pA0[1],3*SCcoeff*pA5[0]*pA5[0],3*SCcoeff*pA5[1]*pA5[1]);
    printf(" neutron SI %.3E [%.3E] SD %.3E [%.3E]\n",
       SCcoeff*nA0[0]*nA0[0],SCcoeff*nA0[1]*nA0[1],3*SCcoeff*nA5[0]*nA5[0],3*SCcoeff*nA5[1]*nA5[1]);

}
#endif

  
#ifdef CDM_NUCLEUS
{ double dNdE[300];
  double nEvents;

printf("\n======== Direct Detection ========\n");    

  nEvents=nucleusRecoil(Maxwell,73,Z_Ge,J_Ge73,SxxGe73,dNdE);

  printf("73Ge: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));
                                                                                                         
#ifdef SHOWPLOTS
    displayRecoilPlot(dNdE,"Distribution of recoil energy of 73Ge",0,199);
#endif

  nEvents=nucleusRecoil(Maxwell,131,Z_Xe,J_Xe131,SxxXe131,dNdE);

  printf("131Xe: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));                                   
#ifdef SHOWPLOTS
    displayRecoilPlot(dNdE,"Distribution of recoil energy of 131Xe",0,199);
#endif

  nEvents=nucleusRecoil(Maxwell,23,Z_Na,J_Na23,SxxNa23,dNdE);

  printf("23Na: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));                                   
#ifdef SHOWPLOTS
    displayRecoilPlot(dNdE,"Distribution of recoil energy of 23Na",0,199);
#endif

  nEvents=nucleusRecoil(Maxwell,127,Z_I,J_I127,SxxI127,dNdE);

  printf("I127: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));                                   
#ifdef SHOWPLOTS
    displayRecoilPlot(dNdE,"Distribution of recoil energy of 127I",0,199);
#endif
  
}
#endif 


#ifdef NEUTRINO
{ double nu[NZ], nu_bar[NZ],mu[NZ];
  int forSun=1;
  double Emin=1;
  
 printf("\n===============Neutrino Telescope=======  for  "); 
 if(forSun) printf("Sun\n"); else printf("Earth\n");  

  err=neutrinoFlux(Maxwell,forSun, nu,nu_bar);
#ifdef SHOWPLOTS
   displaySpectra("neutrino fluxes [1/Year/km^2/GeV]",Emin,Mcdm,2,nu,"nu",nu_bar,"nu_bar");
#endif

   printf(" E>%.1E GeV neutrino/anti-neutrin fluxes   %.2E/%.2E [1/Year/km^2]\n",Emin,
          spectrInfo(Emin,nu,NULL), spectrInfo(Emin,nu_bar,NULL));
 
//ICE CUBE  
if(forSun) printf("IceCube22 exclusion confidence level = %.2E%%\n", 100*exLevIC22(nu,nu_bar,NULL));
  
/* Upward events */
 
Emin=0.1;  
  muonUpward(nu,nu_bar, mu);
#ifdef SHOWPLOTS  
  displaySpectrum("Upward muons[1/Year/km^2/GeV]",Emin,Mcdm/2,mu);
#endif
  printf(" E>%.1E GeV Upward muon flux    %.2E [1/Year/km^2]\n",Emin,spectrInfo(Emin,mu,NULL));
   
/* Contained events */
  muonContained(nu,nu_bar,1., mu);
#ifdef SHOWPLOTS  
  displaySpectrum("Contained  muons[1/Year/km^3/GeV]",Emin,Mcdm,mu); 
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

   pname = "h1";
   width=pWidth(pname,&L);
   printf("\n%s :   total width=%.2E \n and Branchings:\n",pname,width);
   printTxtList(L,stdout);

   pname = "~o2";
   width=pWidth(pname,&L);
   printf("\n%s :   total width=%.2E \n and Branchings:\n",pname,width);
   printTxtList(L,stdout);            
}
#endif

#ifdef CROSS_SECTIONS
{
  double cs, Pcm=4000, Qren,Qfact,pTmin=200;
  int nf=3;
  
  Qfact=pMass(CDM1);
  Qren=pTmin;

  printf("pp -> DM,DM +jet(pt>%.2E GeV)  at %.2E GeV\n",pTmin,Pcm);  
  
  cs=hCollider(Pcm,1,nf,Qren, Qfact, CDM1,aCDM1,pTmin,1);
  printf("cs(pp->~o1,~o2)=%.2E[pb]\n",cs);
 
}
#endif 

#ifdef CLEAN
  killPlots();
  printf("\n==== Delete intermediate files =====\n");    
  system("rm -f Key.dat SM_decay.dat HB.* hb.stdout HS.* hs.stdout UMSSM* debug* Lilith* smodels.* summary.*  particles.py*");
  printf("\n====================================\n");    
#endif 

return 0;
}
