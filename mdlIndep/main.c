/*====== Modules ===============
   Keys to switch on 
   various modules of micrOMEGAs  
================================*/
      

#define INDIRECT_DETECTION  
/* Compute spectra of gamma/positron/antiprotons/neutrinos for DM annihilation; 
   Explicit input parameters are annihilation cross section vcs and channel fractions.
   Integrate gamma signal over DM galactic squared density for given line  of sight; 
   Calculate galactic propagation of positrons and antiprotons.      
*/
  
#define CDM_NUCLEUS
/* Calculate number of events for 1kg*day and recoil energy distibution for various nuclei.
   Input parameters are  Mcdm, spin dependent and spin independent cross sections on nucleons
*/
#define NEUTRINO 
/*  Neutrino signal of DM annihilation in Sun and Earth */
  
/*===== end of Modules  ======*/

/*===== Options ========*/
//#define SHOWPLOTS
#define CLEAN  to clean intermediate files

/*===== End of DEFINE  settings ===== */

#include"../sources/micromegas.h"
#include"../sources/micromegas_aux.h"

int main(int argc,char** argv)
{  int err,n,i;

// data corresponding to MSSM/mssmh.dat     
   double csSIp=5.489E-09,  csSIn=5.758E-09, csSDp=5.807E-05, csSDn=-4.503E-05;  //[pb] 
   double vcs = 0.310; // [pb] 
#define nCH 5   /* annihilation channels */ 
   int    pdgCH[nCH] ={   6,   5,  15,  23,   24};    
   double fracCH[nCH]={0.59,0.12,0.20,0.03, 0.06};
  
   Mcdm=189.0 ;
   

#ifdef INDIRECT_DETECTION
{ 
  double sigmaV;
  char txt[100];
  double SpA[NZ],SpE[NZ],SpP[NZ];
  double FluxA[NZ],FluxE[NZ],FluxP[NZ], buff[NZ];
  double SMmev=320;  /* solar potential in MV */
  double Etest=Mcdm/2;
   
printf("\n==== Indirect detection =======\n");  

  sigmaV=vcs*2.9979E-26; 
  printf("sigmav=%.2E[cm^3/s]\n",sigmaV);  

  SpA[0]=SpE[0]=SpP[0]=Mcdm;
  for(i=1;i<NZ;i++) { SpA[i]=SpE[i]=SpP[i]=0;} 

  for(n=0;n<nCH;n++) if(fracCH[n]>0)
  {    
    basicSpectra(Mcdm,pdgCH[n],0, buff);  for(i=1;i<NZ;i++)  SpA[i]+=buff[i]*fracCH[n];
    basicSpectra(Mcdm,pdgCH[n],1, buff);  for(i=1;i<NZ;i++)  SpE[i]+=buff[i]*fracCH[n];
    basicSpectra(Mcdm,pdgCH[n],2, buff);  for(i=1;i<NZ;i++)  SpP[i]+=buff[i]*fracCH[n];
  }  

  {  double Emin=1; /* Energy cut  in GeV   */
     double fi=0.0,dfi=M_PI/180.; /* angle of sight and 1/2 of cone angle in [rad] */ 

     gammaFluxTab(fi,dfi, sigmaV, SpA,  FluxA);     
     printf("Photon flux  for angle of sight f=%.2f[rad]\n"
     "and spherical region described by cone with angle %.2f[rad]\n",fi,2*dfi);
#ifdef SHOWPLOTS
     sprintf(txt,"Photon flux[cm^2 s GeV]^{1} at f=%.2f[rad], cone angle %.2f[rad]",fi,2*dfi);
     displaySpectrum(txt,Emin,Mcdm,FluxA);
#endif
     printf("Photon flux = %.2E[cm^2 s GeV]^{-1} for E=%.1f[GeV]\n",SpectdNdE(Etest, FluxA), Etest);       
  }

  { double Emin=1.;
    posiFluxTab(Emin, sigmaV, SpE,  FluxE);    
    if(SMmev>0)  solarModulation(SMmev,0.0005,FluxE,FluxE);

#ifdef SHOWPLOTS     
    displaySpectrum("positron flux [cm^2 s sr GeV]^{-1}" ,Emin,Mcdm,FluxE);
#endif
    printf("Positron flux  =  %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
    SpectdNdE(Etest, FluxE),  Etest);           
  }
  
  { double Emin=1;
    pbarFluxTab(Emin, sigmaV, SpP,  FluxP  ); 
    if(SMmev>0)  solarModulation(SMmev,1,FluxP,FluxP);
#ifdef SHOWPLOTS    
     displaySpectrum("antiproton flux [cm^2 s sr GeV]^{-1}" ,Emin, Mcdm,FluxP);
#endif
    printf("Antiproton flux  =  %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
    SpectdNdE(Etest, FluxP),  Etest);             
  }
}  
#endif

  
#ifdef CDM_NUCLEUS
{ double dNdE[300];
  double nEvents;

  printf("\n======== Direct Detection ========\n");    

  nEvents=nucleusRecoilAux(Maxwell,73,Z_Ge,J_Ge73,SxxGe73, csSIp,csSIn, csSDp,csSDn,dNdE);
                              
  printf("73Ge: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));
                                                                                                         
#ifdef SHOWPLOTS
    displayRecoilPlot(dNdE,"Distribution of recoil energy of 73Ge",0,199);
#endif

  nEvents=nucleusRecoilAux(Maxwell,131,Z_Xe,J_Xe131,SxxXe131,csSIp,csSIn, csSDp,csSDn,dNdE);

  printf("131Xe: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));                                   
#ifdef SHOWPLOTS
    displayRecoilPlot(dNdE,"Distribution of recoil energy of 131Xe",0,199);
#endif

  nEvents=nucleusRecoilAux(Maxwell,23,Z_Na,J_Na23,SxxNa23,csSIp,csSIn, csSDp,csSDn,dNdE);

  printf("23Na: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));                                   
#ifdef SHOWPLOTS
    displayRecoilPlot(dNdE,"Distribution of recoil energy of 23Na",0,199);
#endif

  nEvents=nucleusRecoilAux(Maxwell,127,Z_I,J_I127,SxxI127,csSIp,csSIn, csSDp,csSDn,dNdE);

  printf("I127: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));                                   
#ifdef SHOWPLOTS
    displayRecoilPlot(dNdE,"Distribution of recoil energy of 127I",0,199);
#endif
  
}
#endif 

#ifdef NEUTRINO
{ double nu[NZ], nu_bar[NZ],mu[NZ],buff[NZ];
  double Crate,R,Prop;
  int forSun=1;
  double Emin=1;
  int i,err;
  double yrs=31556925.2;
  
  printf("\n===============Neutrino Telescope=======  for  "); 
  if(forSun) printf("Sun\n"); else printf("Earth\n"); 
  
  Crate=captureAux(Maxwell,forSun, Mcdm,csSIp,csSIn,csSDp,csSDn);
    
  nu[0]=nu_bar[0]=Mcdm;
  for(i=1;i<NZ;i++){ nu[i]=nu_bar[i]=0;}

  for(n=0;n<nCH;n++) if(fracCH[n]>0)
  { double bn[NZ],bn_[NZ];   
    basicNuSpectra(forSun,Mcdm, pdgCH[n], 0, bn, bn_);
    for(i=1;i<NZ;i++) { nu[i] +=bn[i]*fracCH[n]; nu_bar[i] +=bn_[i]*fracCH[n];}
  }  


// neutrino flux at Earth
  if(forSun) R=150E6; else  R=6378.1;  // distance in km 
  Prop=yrs/(4*M_PI*R*R);        // in Year*km^2

  for(i=1;i<NZ;i++) { nu[i]*= 0.5*Crate*Prop;  nu_bar[i]*=0.5*Crate*Prop; }

  if(forSun) printf("IceCube22 exclusion confidence level = %.2E%%\n", exLevIC22(nu,nu_bar,NULL));
      
#ifdef SHOWPLOTS
  displaySpectra("neutrino fluxes [1/Year/km^2/GeV]",Emin,Mcdm,2, nu,"nu",nu_bar,"nu_bar");
#endif
{ 
    printf(" E>%.1E GeV neutrino flux       %.3E [1/Year/km^2] \n",Emin,spectrInfo(Emin,nu,NULL));
    printf(" E>%.1E GeV anti-neutrino flux  %.3E [1/Year/km^2]\n", Emin,spectrInfo(Emin,nu_bar,NULL));  
} 


if(forSun)
{
 double  nu[NZ],nuB[NZ]; // arrays for neutrino and antineutrino spectra 
 double  csSDp[6]={1.E-40,9.51E-41,1.04E-40, 2.67E-40,2.06e-39,6.0e-39};
                         // cross section TAB_1 4-th collumn  
 double  Mdm[6]={100,250,500,1000,3000,5000}; // DM masses 
 double exLev=0.9; // 90% exclusion confidence level 
 double Nsig;      // for number of signal events 
             
 WIMPSIM=1;
 
 printf("W channel:\n");
 
 for(int k=0;k<6;k++)   
 {  
    double Crate=captureAux(Maxwell,forSun,Mdm[k],0,0,csSDp[k]*1E36,0); // DM capture rate [s]
    double gamma=Crate*Prop/2;
    double f_90; // improving factor for 90% confidence level exclusion 
    printf("Dark matter mass is %.2f\n",Mdm[k]);
    basicNuSpectra(forSun,Mdm[k],24/*W*/, 0, nu, nu_bar);              // neurtuno spectra 
    for(int i=1;i<NZ;i++) { nu[i]*=gamma; nu_bar[i]*=gamma;}           // neutrino fluxes 
    if(k)
    {  IC22events(nu,nu_bar, 10, &Nsig,NULL, NULL);
       printf("     csSDp=%.2E[cm^2] ==>  Nsignal=%.2f (TAB_1 column 1)\n",csSDp[k],Nsig);
    } 
    f_90=fluxFactorIC22( exLev,nu,nu_bar);
    printf("      90%%  exclusion limit for csSDp is %.2E[cm^2]\n", csSDp[k]*f_90);
    for(int i=1;i<NZ;i++) {nu[i]*=f_90; nu_bar[i]*=f_90;}
    IC22events(nu,nu_bar, 10,&Nsig,NULL,NULL);
    printf("      90%% exclusion limit for number of signal events: %.2f \n",Nsig/1.2);

#ifdef SIGNAL_SIMULATION     
    { 
      double xx[15]={},dxx[15];
      double x_=0;
      for(i=0;i<1000;i++) 
      {  double f_90_=fluxFactorIC22_random( exLev,nu,nu_bar);      
         int n=f_90_*5;
         x_+=f_90_;
         if(n<15) xx[n]++;
      } 
      x_/=1000;  
      for(i=0;i<15;i++)  dxx[i]=sqrt(xx[i]);
      for(i=0;i<15;i++)  { xx[i]/=1000;  dxx[i]/=1000;}
      char s[1000];
      sprintf(s,"Mcdm=%.0f csSDp_exp=%.2E  <csSDp>=%.2E \n",Mdm[k], csSDp[k]*f_90,x_*csSDp[k]*f_90);
      displayPlot(s,0,csSDp[k]*f_90*3,"csSDp",15,1,xx,dxx,"N");
    } 
#endif     
 }
 
 
}  

  
/* Upward events */
  
  muonUpward(nu,nu_bar, mu);
#ifdef SHOWPLOTS  
  displaySpectrum("Upward muons[1/Year/km^2/GeV]",1,Mcdm/2,mu);
#endif

  printf(" E>%.1E GeV Upward muon flux    %.3E [1/Year/km^2]\n",Emin,spectrInfo(Emin,mu,NULL));  
  
/* Contained events */
  muonContained(nu,nu_bar,1., mu);
#ifdef SHOWPLOTS  
  displaySpectrum("Contained  muons[1/Year/km^3/GeV]",Emin,Mcdm,mu); 
#endif
  { 
    double Emin=1; //GeV
    printf(" E>%.1E GeV Contained muon flux %.3E [1/Year/km^3]\n",Emin,spectrInfo(Emin,mu,NULL));
  }  
}        
#endif
 


  killPlots();

#ifdef CLEAN


#endif

  return 0;
}
