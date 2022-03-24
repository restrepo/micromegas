/*====== Modules ===============
   Keys to switch on 
   various modules of micrOMEGAs  
================================*/
      
#define MASSES_INFO      
  /* Display information about mass spectrum  */

#define CONSTRAINTS 

//#define HIGGSBOUNDS
//#define HIGGSSIGNALS 
#define LILITH
//#define SMODELS

#define OMEGA            
  /* Calculate relic density and display contribution of  individual channels */
#define INDIRECT_DETECTION  
  /* Compute spectra of gamma/positron/antiprotons/neutrinos for DM annihilation; 
     Calculate <sigma*v>;
     Integrate gamma signal over DM galactic squared density for given line 
     of sight; 
     Calculate galactic propagation of positrons and antiprotons.      
  */

/*#define RESET_FORMFACTORS*/
  /* Modify default nucleus form factors, 
    DM velocity distribution,
    A-dependence of Fermi-dencity
  */     
#define CDM_NUCLEON     
  /* Calculate amplitudes and cross-sections for  CDM-mucleon collisions */  

#define CDM_NUCLEUS      
  // Calculate  exclusion rate for direct detection experiments Xenon1T, DarkSide50, CRESST, and PICO-60
        
            
#define NEUTRINO 

//#define DECAYS
//#define CROSS_SECTIONS 
#define CLEAN
  
/*===== end of Modules  ======*/

/*===== Options ========*/
/*#define SHOWPLOTS*/
     /* Display  graphical plots on the screen */ 

/*===== End of DEFINE  settings ===== */

#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"
#include"lib/pmodel.h"


int main(int argc,char** argv)
{  int err;
   char cdmName[10];
   int spin2, charge3,cdim;
   
  ForceUG=0 ;  /* to Force Unitary Gauge assign 1 */
  //useSLHAwidth=0;
  if(argc==1)
  { 
      printf(" Correct usage:  ./main  <file with parameters> \n");
      printf("Example: ./main data1.par\n");
      exit(1);
  }
                               
  err=readVar(argv[1]);
  
  if(err==-1)     { printf("Can not open the file\n"); exit(1);}
  else if(err>0)  { printf("Wrong file contents at line %d\n",err);exit(1);}
           
  

  err=sortOddParticles(cdmName);

  if(err) { printf("Can't calculate %s\n",cdmName); return 1;}
  
   qNumbers(cdmName, &spin2, &charge3, &cdim);
   printf("\nDark matter candidate is '%s' with spin=%d/2\n",
    cdmName,       spin2); 
   if(charge3) { printf("Dark Matter has electric charge %d/3\n",charge3); exit(1);}
   if(cdim!=1) { printf("Dark Matter is a color particle\n"); exit(1);}
#ifdef MASSES_INFO
{
  printf("\n=== MASSES OF HIGGS AND ODD PARTICLES: ===\n");
  printHiggs(stdout);
  printMasses(stdout,1);   
}
#endif

#ifdef CONSTRAINTS
{ double csLim;  
  if(Zinvisible()) printf("Excluded by Z->invisible\n");
  if(LspNlsp_LEP(&csLim)) printf("LEP excluded by e+,e- -> DM q qbar  Cross section= %.2E pb\n",csLim);

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
   char Lilith_version[50];

   if(LilithMO("Lilith_in.xml"))
   {        
#include "../include/Lilith.inc"
      if(ndf)
      {  printf("LILITH(DB%s):  -2*log(L): %.2f; -2*log(L_reference): %.2f; ndf: %d; p-value: %.2E \n", 
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



#ifdef OMEGA
{ int fast=1;
  double Beps=1.E-5, cut=0.01;
  double Omega,Xf;   
  int i;
  printf("\n==== Calculation of relic density =====\n");  
  
// to exclude processes with virtual W/Z in DM   annihilation      
   VZdecay=0; VWdecay=1; cleanDecayTable(); 

// to include processes with virtual W/Z  also  in co-annihilation 
//   VZdecay=2; VWdecay=2; cleanDecayTable(); 
  
  Omega=darkOmega(&Xf,fast,Beps,&err);
  printf("Xf=%.2e Omega=%.2e\n",Xf,Omega);
  if(Omega>0)printChannels(Xf,cut,Beps,1,stdout);   

  VZdecay=1; VWdecay=1; cleanDecayTable();  // restore default
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

  sigmaV=calcSpectrum(4,SpA,SpE,SpP,SpNe,SpNm,SpNl ,&err);
    /* Returns sigma*v in cm^3/sec.     SpX - calculated spectra of annihilation.
       Use SpectdNdE(E, SpX) to calculate energy distribution in  1/GeV units.
       
       First parameter 1-includes W/Z polarization
                       2-includes gammas for 2->2+gamma
                       4-print cross sections             
    */
  printf("sigmav=%.2E[cm^3/s]\n",sigmaV);  


  { 
     double fi=0.1,dfi=0.05; /* angle of sight and 1/2 of cone angle in [rad] */ 

     gammaFluxTab(fi,dfi, sigmaV, SpA,  FluxA);     
     printf("Photon flux  for angle of sight f=%.2f[rad]\n"
     "and spherical region described by cone with angle %.2f[rad]\n",fi,2*dfi);
#ifdef SHOWPLOTS
     sprintf(txt,"Photon flux[cm^2 s GeV]^{1} at f=%.2f[rad], cone angle %.2f[rad]",fi,2*dfi);
     displayPlot(txt,"E[GeV]",Emin,Mcdm,0,1,"flux",0,SpectdNdE,FluxA);
#endif
     printf("Photon flux = %.2E[cm^2 s GeV]^{-1} for E=%.1f[GeV]\n",SpectdNdE(Etest, FluxA), Etest);       
  }

  { 
    posiFluxTab(Emin, sigmaV, SpE,  FluxE);
#ifdef SHOWPLOTS     
    displayPlot("positron flux [cm^2 s sr GeV]^{-1}" ,"E[GeV]",Emin,Mcdm,0,1,"flux",0,SpectdNdE,FluxE);
#endif
    printf("Positron flux  =  %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
    SpectdNdE(Etest, FluxE),  Etest);           
  }
  
  { 
    pbarFluxTab(Emin, sigmaV, SpP,  FluxP  ); 
#ifdef SHOWPLOTS    
     displayPlot("antiproton flux [cm^2 s sr GeV]^{-1}" ,"E[GeV]",Emin, Mcdm,0,1,"flux",0,SpectdNdE,FluxP);
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

   calcScalarQuarkFF( Mu/Md, Ms/Md, sigmaPiN[MeV], sigmaS[MeV])  
   calculates and rewrites Scalar form factors
*/
  printf("\n======== RESET_FORMFACTORS ======\n");
 
  printf("protonFF (default) d %.2E, u %.2E, s %.2E\n",ScalarFFPd, ScalarFFPu,ScalarFFPs);                               
  printf("neutronFF(default) d %.2E, u %.2E, s %.2E\n",ScalarFFNd, ScalarFFNu,ScalarFFNs);

//                    To restore default form factors of  version 2  call 
     calcScalarQuarkFF(0.553,18.9,55.,243.5);


  printf("protonFF (new)     d %.2E, u %.2E, s %.2E\n",ScalarFFPd, ScalarFFPu,ScalarFFPs);                               
  printf("neutronFF(new)     d %.2E, u %.2E, s %.2E\n",ScalarFFNd, ScalarFFNu,ScalarFFNs);

//                    To restore default form factors of current version  call 
//  calcScalarQuarkFF(0.56,20.2,34,42);


}
#endif

#ifdef CDM_NUCLEON
{ double pA0[2],pA5[2],nA0[2],nA5[2];
  double Nmass=0.939; /*nucleon mass*/
  double SCcoeff;
   double csSIp,csSIn,csSDp,csSDn, csSIp_,csSIn_,csSDp_,csSDn_;        

printf("\n==== Calculation of CDM-nucleons amplitudes  =====\n");   
  nucleonAmplitudes(CDM1, pA0,pA5,nA0,nA5);
  
    printf("%s[%s]-nucleon micrOMEGAs amplitudes\n",CDM1,aCDM1? aCDM1:CDM1 );
    printf("proton:  SI  %.3E [%.3E]  SD  %.3E [%.3E]\n",pA0[0], pA0[1],  pA5[0], pA5[1] );
    printf("neutron: SI  %.3E [%.3E]  SD  %.3E [%.3E]\n",nA0[0], nA0[1],  nA5[0], nA5[1] ); 

    SCcoeff=4/M_PI*3.8937966E8*pow(Nmass*Mcdm1/(Nmass+ Mcdm1),2.);
    csSIp=  SCcoeff*pA0[0]*pA0[0];  csSIp_=  SCcoeff*pA0[1]*pA0[1];
    csSDp=3*SCcoeff*pA5[0]*pA5[0];  csSDp_=3*SCcoeff*pA5[1]*pA5[1];
    csSIn=  SCcoeff*nA0[0]*nA0[0];  csSIn_=  SCcoeff*nA0[1]*nA0[1];
    csSDn=3*SCcoeff*nA5[0]*nA5[0];  csSDn_=3*SCcoeff*nA5[1]*nA5[1];
    
    printf("%s[%s]-nucleon cross sections[pb] :\n",CDM1,aCDM1);
    printf(" proton  SI %.3E [%.3E] SD %.3E [%.3E]\n", csSIp,csSIp_,csSDp,csSDp_);
    printf(" neutron SI %.3E [%.3E] SD %.3E [%.3E]\n", csSIn,csSIn_,csSDn,csSDn_);           
}
#endif
  
#ifdef CDM_NUCLEUS
{ char * expName; 
  printf("\n===== Direct detection exclusion:======\n");
  double pval=DD_pval(AllDDexp, Maxwell, &expName);
//printf("pval=%E  XENON1T_90(Mcdm)*1E36=%E expNmae=%p \n",pval,XENON1T_90(Mcdm)*1E36,expName);  
  if(pval<0.1 )  printf("Excluded by %s  %.1f%%\n", expName, 100*(1-pval)); 
  else printf("Not excluded by DD experiments  at 90%% level \n"); 
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
  displayPlot("neutrino fluxes [1/Year/km^2/GeV]","E[GeV]",Emin,Mcdm,0,2,"nu",0,SpectdNdE,nu,"nu_bar",0,SpectdNdE,nu_bar);
#endif
  printf(" E>%.1E GeV neutrino/anti-neutrino fluxes   %.2E/%.2E [1/Year/km^2]\n",Emin,
          spectrInfo(Emin,nu,NULL), spectrInfo(Emin,nu_bar,NULL));

//  ICE CUBE
  if(forSun) printf("IceCube22 exclusion confidence level = %.2E%%\n", 100*exLevIC22(nu,nu_bar,NULL));

/* Upward events */

  muonUpward(nu,nu_bar, mu);
#ifdef SHOWPLOTS  
  displayPlot("Upward muons[1/Year/km^2/GeV]","E[GeV]",1,Mcdm/2,0,1,"mu",0,SpectdNdE,mu);
#endif
  printf(" E>%.1E GeV Upward muon flux    %.3E [1/Year/km^2]\n",Emin,spectrInfo(Emin,mu,NULL));

/* Contained events */
  muonContained(nu,nu_bar,1., mu);
#ifdef SHOWPLOTS  
  displayPlot("Contained  muons[1/Year/km^3/GeV]","E[GeV]",Emin,Mcdm,0,1,"mu",0,SpectdNdE,mu);
#endif
  printf(" E>%.1E GeV Contained muon flux %.3E [1/Year/km^3]\n",Emin,spectrInfo(Emin,mu,NULL));
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
   printf("\n%s :   total width=%.3E \n and Branchings:\n",pname,width);
   printTxtList(L,stdout);

   pname = "~x2";
   width=pWidth(pname,&L);
   printf("\n%s :   total width=%.3E \n and Branchings:\n",pname,width);
   printTxtList(L,stdout);

   pname = "~H+";
   width=pWidth(pname,&L);
   printf("\n%s :   total width=%.3E \n and Branchings:\n",pname,width);
   printTxtList(L,stdout);

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
  system("rm -f nngg.in nngg.out"); 
  system("rm -f HB.* HS.* hb.* hs.*  debug_channels.txt debug_predratio.txt  Key.dat");
  system("rm -f Lilith_*   particles.py*");
//  system("rm -f   smodels.*");  
#endif 
  killPlots();
  return 0;
}
