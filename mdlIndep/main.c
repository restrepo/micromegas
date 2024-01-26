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
  
//#define NEUTRINO 
/*  Neutrino signal of DM annihilation in Sun and Earth */
  
/*===== end of Modules  ======*/

/*===== Options ========*/
//#define SHOWPLOTS
#define CLEAN  to clean intermediate files

/*===== End of DEFINE  settings ===== */

#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"


int main(int argc,char** argv)
{ int err,n,i;


{ Mcdm=500;
  double tabWA[NZ], tabWe[NZ], tabWp[NZ],
         tabBA[NZ], tabBe[NZ], tabBp[NZ];
  basicSpectra(Mcdm ,24, 0 , tabWA);
  basicSpectra(Mcdm ,24, 1 , tabWe);
  basicSpectra(Mcdm ,24, 2 , tabWp);
  displayPlot("W->EdN/dE","E",0,Mcdm,0,3,"gamma",0, eSpectdNdE,tabWA
                                        ,"positron",0,eSpectdNdE,tabWe
                                   ,"a-proton",0,eSpectdNdE,tabWp);

  basicSpectra(Mcdm ,6, 0 , tabBA);
  basicSpectra(Mcdm ,6, 1 , tabBe);
  basicSpectra(Mcdm ,6, 2 , tabBp);
  displayPlot("b->EdN/dE","E",0,Mcdm,0,3,"gamma",0, eSpectdNdE,tabBA
                                        ,"positron",0,eSpectdNdE,tabBe
                                   ,"a-proton",0,eSpectdNdE,tabBp);


  exit(0);
}      

{
   Mcdm=5;
   int pdg=13;
   double tab[10*NZ], tabA[10*NZ] ;
   double sum=0,sumA=0;

   for(int i=0;i<6;i++)
   {
      basicSpectra(Mcdm ,pdg, i , tab);
      basicSpectraA(Mcdm,pdg, i , tabA);
      char mess[30];
      sprintf(mess,"EdNdE for %s", outNames[i]);
      displayPlot(mess,"E", 1E-3 /*Mcdm/100*/,Mcdm*1.05,1,2
                                                  ,  "old",0,eSpectdNdE,tab
                                                    ,"Australia",0,eSpectdNdE,tabA);
      double dSum=simpson_arg(eSpectdNdE,tabA,1E-4, Mcdm, 1E-3,NULL);
      dSum/=Mcdm;
      if(i==0) sumA+=dSum; else sumA+=2*dSum;
      dSum=simpson_arg(eSpectdNdE,tab,1E-4, Mcdm, 1E-3,NULL);
      dSum/=Mcdm;
     if(i==0) sum+=dSum; else sum+=2*dSum;
  }
  printf(" Energy conservation sum=%E sumA=%E\n",  sum,sumA);

  exit(0); 
}

// data corresponding to MSSM/mssmh.dat     
// cross sections of DM-nucleon interaction 
   double csSIp=5.489E-09,  csSIn=5.758E-09, csSDp=5.807E-05, csSDn=-4.503E-05;  //[pb] 
// cross section of DM annihilation in halo
   double vcs = 0.310; // [pb] 
#define nCH 5   /* annihilation channels and their  fractions*/ 
   int    pdgCH[nCH] ={   6,   5,  15,  23,   24};    
   double fracCH[nCH]={0.59,0.12,0.20,0.03, 0.06};

   sortOddParticles(NULL);
// ======= millicharge   

/*
{     
 double  M=300, m1=150,m2=150, w1=1,w2=1;
 
w1=0;
w2=20;
for(w2=0.0001,w1=w2;w2<1;w2+=0.001) printf("w2=%E  decayPcmW=%E  decayPcm'=%E \n",w2,(double)decayPcmW(M,m1,m2,w1,w2,0)/sqrt(w2), (double)decayPcm(M,m1-w2/5,m2-w2/5)/sqrt(w2) );
 

// printf(" decayPcm=%E, decayPcmW=%E, decayPcm'=%E\n",(double)decayPcm(M,m1,m2), (double)decayPcmW(M,m1,m2,w1,w2,0), (double)decayPcm(M,m1-w1/20,m2-w2/20));
 exit(0);  

} 
*/
// ===   

// Mass of Dark Matter    
   Mcdm=189.0;

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

// Calculation on photon, positron and antiproton spectra: summation over channels 
  for(n=0;n<nCH;n++) if(fracCH[n]>0)
  {    
    basicSpectra(Mcdm,pdgCH[n],0, buff);  for(i=1;i<NZ;i++)  SpA[i]+=buff[i]*fracCH[n];
    basicSpectra(Mcdm,pdgCH[n],1, buff);  for(i=1;i<NZ;i++)  SpE[i]+=buff[i]*fracCH[n];
    basicSpectra(Mcdm,pdgCH[n],2, buff);  for(i=1;i<NZ;i++)  SpP[i]+=buff[i]*fracCH[n];
  }  

// Photon flux
  {  double Emin=1; /* Energy cut  in GeV   */
     double fi=0.0,dfi=M_PI/180.; /* angle of sight and 1/2 of cone angle in [rad] */ 

     gammaFluxTab(fi,dfi, sigmaV, SpA,  FluxA);     
     printf("Photon flux  for angle of sight f=%.2f[rad]\n"
     "and sph[Berical region described by cone with angle %.2f[rad]\n",fi,2*dfi);
#ifdef SHOWPLOTS
     sprintf(txt,"Photon flux[cm^2 s GeV]^{1} at f=%.2f[rad], cone angle %.2f[rad]",fi,2*dfi);
     displayPlot(txt,"E[GeV]",Emin,Mcdm,0,1,"flux",0,SpectdNdE,FluxA);
#endif
     printf("Photon flux = %.2E[cm^2 s GeV]^{-1} for E=%.1f[GeV]\n",SpectdNdE(Etest, FluxA), Etest);       
  }
  
// Positron flux
  { double Emin=1.;
    posiFluxTab(Emin, sigmaV, SpE,  FluxE);    
    if(SMmev>0)  solarModulation(SMmev,0.0005,FluxE,FluxE);
#ifdef SHOWPLOTS     
    displayPlot("positron flux [cm^2 s sr GeV]^{-1}" ,"E[GeV]",Emin,Mcdm,0,1,"flux",0,SpectdNdE,FluxE);
#endif
    printf("Positron flux  =  %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
    SpectdNdE(Etest, FluxE),  Etest);           
  }

// antiproton flux
  { double Emin=1;
    pbarFluxTab(Emin, sigmaV, SpP,  FluxP  ); 
    if(SMmev>0)  solarModulation(SMmev,1,FluxP,FluxP);
#ifdef SHOWPLOTS    
     displayPlot("antiproton flux [cm^2 s sr GeV]^{-1}" ,"E[GeV]",Emin, Mcdm,0,1,"flux",0,SpectdNdE,FluxP);
#endif
    printf("Antiproton flux  =  %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
    SpectdNdE(Etest, FluxP),  Etest);             
  }
}  
#endif


// Neutrino telescope 

#ifdef NEUTRINO
{ double nu[NZ], nu_bar[NZ],mu[NZ],buff[NZ];
  double Crate,R,Prop;
  int forSun=1;
  double Emin=1;
  int i,err;
  double yrs=31556925.2;
  
  printf("\n===============Neutrino Telescope=======  for  "); 
  if(forSun) printf("Sun\n"); else printf("Earth\n"); 

// Capture rate  
  Crate=captureCS(Maxwell,forSun, Mcdm,csSIp,csSIn,csSDp,csSDn);
    
  nu[0]=nu_bar[0]=Mcdm;
  for(i=1;i<NZ;i++){ nu[i]=nu_bar[i]=0;}

// Calculation of neutrino, antineutrino spectra: sumation over chanels
  for(n=0;n<nCH;n++) if(fracCH[n]>0)
  { double bn[NZ],bn_[NZ];   
    basicNuSpectra(forSun,Mcdm, pdgCH[n], 0, bn, bn_);
    for(i=1;i<NZ;i++) { nu[i] +=bn[i]*fracCH[n]; nu_bar[i] +=bn_[i]*fracCH[n];}
  }  


// propagation: neutrino flux at Earth
  if(forSun) R=150E6; else  R=6378.1;  // distance in km 
  Prop=yrs/(4*M_PI*R*R);        // in Year*km^2

  for(i=1;i<NZ;i++) { nu[i]*= 0.5*Crate*Prop;  nu_bar[i]*=0.5*Crate*Prop; }

//  Display  neutrino flux  
#ifdef SHOWPLOTS
  displayPlot("neutrino fluxes [1/Year/km^2/GeV]","E[GeV]",Emin,Mcdm,0,2,"nu",0,SpectdNdE,  nu, "nu_bar",0,SpectdNdE,nu_bar);
#endif
{ 
    printf(" E>%.1E GeV neutrino flux       %.3E [1/Year/km^2] \n",Emin,spectrInfo(Emin,nu,NULL));
    printf(" E>%.1E GeV anti-neutrino flux  %.3E [1/Year/km^2]\n", Emin,spectrInfo(Emin,nu_bar,NULL));  
} 

//Exclusion level 

  if(forSun) printf("IceCube22 exclusion confidence level = %.2E%%\n", exLevIC22(nu,nu_bar,NULL));
}      
#endif

#ifdef CLEAN


#endif

  killPlots();

  return 0;
}
