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

//#define HIGGSBOUNDS 
//#define HIGGSSIGNALS
#define LILITH 
//#define SMODELS 


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
      /* Modify default nucleus form factors */     
#define CDM_NUCLEON 
      /* Calculate amplitudes and cross-sections for 
         CDM-mucleon collisions 
      */  
//#define CDM_NUCLEUS     
      // Calculate  exclusion rate for direct detection experiments Xenon1T and DarkSide50
           
            
#define NEUTRINO // neutrino telescope signals 
      
//#define DECAYS
    /* Calculate decay widths and branchings  */
      
//#define CROSS_SECTIONS 
      /* Calculate cross sections and widths for 
         reactions specified by the user
      */
/*===== end of Modules  ======*/

/*===== Options ========*/
/*#define SHOWPLOTS */
     /* Display  graphical plots on the screen */ 

//#define CLEAN   to clean intermidiate files 

/*===== End of DEFINE  settings ===== */


#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"
#include"lib/pmodel.h"


int main(int argc,char** argv)
{  int err;
   char cdmName[10];
   int spin2, charge3,cdim;

   ForceUG=0;    /* to Force Unitary Gauge assign 1 */
   //useSLHAwidth=0;
  if(argc==1)
  { 
      printf(" Correct usage:  ./main  <file with NMSSM parameters> \n");
      printf("Example: ./main data1.par\n");
      exit(1);
  }
                               
/*  err=readVar(argv[1]);*/
  err=readVarCPVMSSM(argv[1]);

  
  if(err==-1)     {printf("Can not open the file\n"); exit(1);}
  else if(err>0)  { printf("Wrong file contents at line %d\n",err);exit(1);}

  if(err) 
  { printf("Problem with RGE solution or spectrum calculation\n");
     exit(1);
  }

  err=sortOddParticles(cdmName);  
// to load CPsuperH decay tables  uncomment 
//slhaRead("cpsuperh2_slha.out",1+8);


  if(err) { printf("Can't calculate %s\n",cdmName); return 1;}
   qNumbers(cdmName,&spin2, &charge3, &cdim);
   printf("\nDark matter candidate is '%s' with spin=%d/2 \n",
    cdmName,       spin2); 
   if(charge3) { printf("Dark Matter has electric charge %d/3\n",charge3); exit(1);}
   if(cdim!=1) { printf("Dark Matter is a color particle\n"); exit(1);}
   if(strcmp(cdmName,"~o1")) printf(" ~o1 is not CDM\n"); 
                    else o1Contents(stdout);

// slhaRead("cpsuperh2_slha.out",1+8);



// If you would like to use widths calculated by NMSSMTools, 
// uncomment the line below

//   slhaRead("cpsuperh2_slha.out",1+8);


#ifdef MASSES_INFO
{
  printf("\n=== MASSES OF HIGGS AND SUSY PARTICLES: ===\n");
  printHiggs(stdout);
  printMasses(stdout,1);
}
#endif

#ifdef CONSTRAINTS
{ double rat,csLim;
  FILE *f;
  printf("\n==== Physical Constraints: =====\n");  
  printf("EDM(Thallium)=%.2E\n",               slhaValFormat("FDIPOLE",0.,"1000812050  1 1 %lf")  ) ;
  printf("EDM(Mercury)=%.2E",                  slhaValFormat("FDIPOLE",0.,"1000801990  1 1 %lf"));
  printf("  # %s\n" ,slhaComment) ;
//  printf("EDM(electron)=%.2E\n",               slhaValFormat("dipoleMOM",0.,"1000812030     1 %lf") ) ;

  printf("EDM(muon)=%.2E\n",                   slhaValFormat("FDIPOLE",0.,"13     1 1 %lf"));

  printf("SUSY contr. to muon Mag.Mom.=%.2E\n",slhaValFormat("FDIPOLE",0.,"13     2 1 %lf"));   

  printf("Br(Bs->mu+,mu-)=%.2E\n",             slhaValFormat("FOBS",0.,"531 1 1  %lf  %*f  2  13 -13"));
  printf("Br(Bd->tau+,tau-)=%.2E\n",           slhaValFormat("FOBS",0.,"511 1 1  %lf  %*f  2 -15 15"));
  printf("Br(B->Xs,gamma)=%.2E\n",             slhaValFormat("FOBS",0.,"5   1 1  %lf  %*f  2   3 22"));
  printf("CP_asymm(B->Xs,gamma)=%.2E[%%]\n",   slhaValFormat("FOBS",0.,"5   3 1  %lf  %*f  2   3 22"));
  printf("Br(Bu->tau,nu)/Br(SM)=%.2E\n",       slhaValFormat("FOBS",0.,"521 2 1  %lf  %*f  2 -15 16"));

//  printf("SUSY mass diff, Bd-bar{Bd} =%.2E[1/ps]\n", slhaValFormat("deltaM",0.,"511 1 %lf"));
//  printf("SUSY mass diff, Bs-bar{Bs} =%.2E[1/ps]\n", slhaValFormat("deltaM",0.,"531 1 %lf"));  

  if(Zinvisible()) printf("Excluded by Z->invisible\n");
  if(LspNlsp_LEP(NULL)) printf("LEP excluded by e+,e- -> DM q qbar Cross Section=%2E pb\n",csLim);
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
// LilithMO("Lilith_in.xml");  
   if(LilithMO("Lilith_in.xml"))
   {        
#include "../include/Lilith.inc"
      if(ndf)
      {
        printf("LILITH(DB%s):  -2*log(L): %.2f; -2*log(L_reference): %.2f; ndf: %d; p-value: %.2E \n", 
        Lilith_version,m2logL,m2logL_reference,ndf,pvalue);
      }  
   }  else printf("LILITH: there is no Higgs candidate\n");
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
  double Beps=1.E-4, cut=0.01;
  double Omega,Xf; 
  
// to exclude processes with virtual W/Z in DM   annihilation      
   VZdecay=0; VWdecay=0; cleanDecayTable(); 

// to include processes with virtual W/Z  also  in co-annihilation 
//   VZdecay=2; VWdecay=2; cleanDecayTable(); 
     
  printf("\n==== Calculation of relic density =====\n");  
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

  sigmaV=calcSpectrum(1+4,SpA,SpE,SpP,SpNe,SpNm,SpNl ,&err);
    /* Returns sigma*v in cm^3/sec.     SpX - calculated spectra of annihilation.
       Use SpectdNdE(E, SpX) to calculate energy distribution in  1/GeV units.
       
       First parameter 1-includes W/Z polarization
                       2-includes gammas for 2->2+gamma
                       4-print cross sections             
    */
  printf("sigmav=%.2E[cm^3/s]\n",sigmaV);  


  { 
     double fi=0.1,dfi=0.05; /* angle of sight and 1/2 of cone angle in [rad] */ 

     gammaFluxTab(fi,dfi, sigmaV, SpA, FluxA);     
    printf("Photon flux for angle of sight f=%.2f[rad]\n"
     "and spherical region described by cone with angle %.2f[rad]\n",fi,2*dfi);

#ifdef SHOWPLOTS
     sprintf(txt,"Photon flux[cm^2 s GeV]^{1} at f=%.2f[rad], cone angle %.2f[rad] ",fi,2*dfi);
     displayPlot(txt,"E[GeV]",Emin,Mcdm,0,1,"flux",0,SpectdNdE,FluxA);
#endif
     printf("Photon flux = %.2E[cm^2 s GeV]^{-1} for E=%.1f[GeV]\n",SpectdNdE(Etest, FluxA),Etest);       
  }

#ifdef LoopGAMMA
{
     double fi=0.1,dfi=0.05; /* angle of sight and 1/2 of cone angle in [rad] */ 
     double vcs_gz,vcs_gg;
     if(loopGamma(&vcs_gz,&vcs_gg)==0)
     {
         printf("Gamma  ray lines:\n");
         printf("E=%.2E[GeV]  vcs(Z,A)= %.2E[cm^3/s], flux=%.2E[cm^2 s]^{-1}\n",Mcdm-91.19*91.19/4/Mcdm,vcs_gz,
                               gammaFlux(fi,dfi,vcs_gz));  
         printf("E=%.2E[GeV]  vcs(A,A)= %.2E[cm^3/s], flux=%.2E[cm^2 s]^{-1}\n",Mcdm,vcs_gg, 
                             2*gammaFlux(fi,dfi,vcs_gg));
     }
}     
#endif     

  { 
    posiFluxTab(Emin, sigmaV, SpE, FluxE);
#ifdef SHOWPLOTS     
    displayPlot("positron flux [cm^2 s sr GeV]^{-1}","E[GeV]",Emin,Mcdm,0,1,"flux",0,SpectdNdE,FluxE);
#endif
    printf("Positron flux  =  %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
    SpectdNdE(Etest, FluxE),  Etest);           
  }
  
  { 
    pbarFluxTab(Emin, sigmaV, SpP,  FluxP  ); 
#ifdef SHOWPLOTS    
     displayPlot("antiproton flux [cm^2 s sr GeV]^{-1}" ,"E[GeV]",Emin,Mcdm,0,1,"flux",0,SpectdNdE, FluxP);
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

  printf("protonFF (default) d %E, u %E, s %E\n",ScalarFFPd, ScalarFFPu,ScalarFFPs);                               
  printf("neutronFF(default) d %E, u %E, s %E\n",ScalarFFNd, ScalarFFNu,ScalarFFNs);
  

//        To restore default form factors of  version 2  call 
    calcScalarQuarkFF(0.553,18.9,55.,243.5);
 
  printf("protonFF (new)     d %E, u %E, s %E\n",ScalarFFPd, ScalarFFPu,ScalarFFPs);                               
  printf("neutronFF(new)     d %E, u %E, s %E\n",ScalarFFNd, ScalarFFNu,ScalarFFNs);
//                    To restore default form factors of current version  call 
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
    nucleonAmplitudes(CDM1, pA0,pA5,nA0,nA5);
    printf("CDM-nucleon micrOMEGAs amplitudes:\n");
    printf("proton:  SI  %.3E  SD  %.3E\n",pA0[0],pA5[0]);
    printf("neutron: SI  %.3E  SD  %.3E\n",nA0[0],nA5[0]); 
    SCcoeff=4/M_PI*3.8937966E8*pow(Nmass*Mcdm/(Nmass+ Mcdm),2.);

    csSIp=  SCcoeff*pA0[0]*pA0[0];
    csSDp=3*SCcoeff*pA5[0]*pA5[0];
    csSIn=  SCcoeff*nA0[0]*nA0[0];
    csSDn=3*SCcoeff*nA5[0]*nA5[0];

    
    printf("CDM-nucleon cross sections[pb]:\n");

    printf(" proton  SI %.3E  SD %.3E \n", csSIp,csSDp);
    printf(" neutron SI %.3E  SD %.3E \n", csSIn,csSDn);
}
#endif
  
#ifdef CDM_NUCLEUS
{ char* expName; 
  printf("\n===== Direct detection exclusion:======\n");
  double pval=DD_pval(AllDDexp, Maxwell, &expName);
  if(pval<0.1 )  printf("Excluded by %s  %.1f%%\n", expName, 100*(1-pval)); 
  else printf("Not excluded by DD experiments  at 90%% level \n"); 
}
#endif 



#ifdef NEUTRINO
{ double nu[NZ], nu_bar[NZ],mu[NZ];
  double Ntot;
  int forSun=1;
  double Emin=1;

 printf("\n===============Neutrino Telescope=======  for  ");
 if(forSun) printf("Sun\n"); else printf("Earth\n");
  err=neutrinoFlux(Maxwell,forSun, nu,nu_bar);
#ifdef SHOWPLOTS
    displayPlot("neutrino fluxes [1/Year/km^2/GeV]","E[GeV]",Emin,Mcdm,0, 2,"dnu/dE",0,SpectdNdE,nu,"dnu_bar/dE",0,SpectdNdE,nu_bar)
#endif

  printf(" E>%.1E GeV neutrino/anti-neutrino fluxes   %.2E/%.2E [1/Year/km^2]\n",Emin,
              spectrInfo(Emin,nu,NULL), spectrInfo(Emin,nu_bar,NULL));  

// ICE CUBE 22
if(forSun) printf("IceCube22 exclusion confidence level = %.2E%%\n", 100*exLevIC22(nu,nu_bar,NULL));

/* Upward events */
  muonUpward(nu,nu_bar, mu);
#ifdef SHOWPLOTS  
  displayPlot("Upward muons[1/Year/km^2/GeV]","E[GeV]",1,Mcdm/2,0,1,"flux",0,SpectdNdE,mu);
#endif

  printf(" E>%.1E GeV Upward muon flux    %.3E [1/Year/km^2]\n",Emin,spectrInfo(Emin,mu,NULL));

/* Contained events */
  muonContained(nu,nu_bar,1., mu);
#ifdef SHOWPLOTS  
  displayPlot("Contained  muons[1/Year/km^3/GeV]","E[GeV]",Emin,Mcdm,0,1,"flux",0,SpectdNdE,mu);
#endif
    printf(" E>%.1E GeV Contained muon flux %.3E [1/Year/km^3]\n",Emin,spectrInfo(Emin,mu,NULL));
}

#endif

#ifdef DECAYS
{  
  txtList L;
   int dim;
   double width,br;
   char * pname;

   pname = "h1";
    width=pWidth(pname,&L);
    printf("Width(%s)=%.3E \n  Branchings:\n",pname,width);
    printTxtList(L,stdout);

   pname = "l";
    width=pWidth(pname,&L);
    printf("Width(%s)=%.3E \n and Branchings:\n",pname,width);
    printTxtList(L,stdout);
    printf("Br(e,Ne,nl)= %.3E\n",findBr(L,"e,Ne,nl"));

   pname = "~o2";
    width=pWidth(pname,&L);
    printf("Width(%s)=%.3E \n and Branchings:\n",pname,width);
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
 
  system(" rm -f CPsuperH.in cpsuperh2_slha.out CPsuperH.out");
  system(" rm -f nngg.in nngg.out");
  system("rm -f HB.* HS.* hb.* hs.*  debug_channels.txt debug_predratio.txt  Key.dat");
  system("rm -f Lilith_*   particles.py*");
//  system("rm -f  smodels.*");  
#endif

  killPlots();
  return 0;
}

