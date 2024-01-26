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

extern double complex HggF(double tau);

double hgg(double tau) { return cabs(HggF(tau));}


static double  fe(double E) //  Fig.3   1506.03811
{
  double ln10E[40]= { 3.6990, 3.9315, 4.1640, 4.3965, 4.6291, 4.8616, 5.0941, 5.3267, 5.5592, 5.7917, 6.0242, 6.2568, 6.4893, 6.7218, 6.9543, 7.1869, 7.4194, 7.6519, 7.8844, 8.1170, 8.3495, 8.5820, 8.8145, 9.0471, 9.2796, 9.5121, 9.7446, 9.9772, 10.2097, 10.4422, 10.6747, 10.9073, 11.1398, 11.3723, 11.6048, 11.8374, 12.0699, 12.3024, 12.5349, 12.7675};
  double  e[40]=    { 0.3084, 0.3072, 0.3044, 0.2983, 0.2944, 0.2819, 0.2593, 0.2283, 0.1904, 0.1527, 0.1770, 0.3680, 0.6173, 0.8228, 0.9298, 0.9756, 0.9932, 0.9889, 0.9288, 0.8256, 0.7032, 0.5589, 0.4296, 0.3952, 0.4660, 0.5714, 0.6388, 0.6243,  0.5593,  0.4868,  0.4502,  0.4357,  0.4083,  0.3994,  0.4049,  0.4064,  0.4070,  0.4016,  0.4036,  0.4056};
  return polint3(log10(E), 40,ln10E,e);
}


static double fA(double E)  //  Fig.3  1506.03811
{
  double ln10E[40]= { 3.6990, 3.9315, 4.1640, 4.3965, 4.6291, 4.8616, 5.0941, 5.3267, 5.5592, 5.7917, 6.0242, 6.2568, 6.4893, 6.7218, 6.9543, 7.1869, 7.4194, 7.6519, 7.8844, 8.1170, 8.3495, 8.5820, 8.8145, 9.0471, 9.2796, 9.5121, 9.7446, 9.9772, 10.2097, 10.4422, 10.6747, 10.9073, 11.1398, 11.3723, 11.6048, 11.8374, 12.0699, 12.3024, 12.5349, 12.7675};
  double A[40]=     { 0.8767, 0.7697, 0.7021, 0.6842, 0.6963, 0.6748, 0.5851, 0.4591, 0.3256, 0.2161, 0.1448, 0.1631, 0.2924, 0.4681, 0.5957, 0.6781, 0.7180, 0.7260, 0.7195, 0.6858, 0.6288, 0.5462, 0.4469, 0.3591, 0.3184, 0.3378, 0.3905, 0.4228,  0.4212,  0.3902,  0.3307,  0.3901,  0.4332,  0.4148,  0.4029,  0.4047,  0.4069,  0.4017,  0.4037,  0.4056};
  return polint3(log10(E), 40,ln10E,A);
}

static double intEe(double E, double *tab) { return eSpectdNdE(E,tab)*fe(E*1E9);}
static double intEA(double E, double *tab) { return eSpectdNdE(E,tab)*fA(E*1E9);} 


double   Planck(double vSigma, double *Sg,double*Se)  // 1506.03811
{ int err;
  double Egmin=5E-6, Eemin=5E-6; 
  if(Eemin<1E-7*Se[0]) Eemin=1E-7*Se[0];
  double r=simpson_arg(intEA,Sg, Egmin,Sg[0],1E-3,NULL)/2 + simpson_arg(intEe,Se,Eemin,Se[0], 1E-3,NULL);
  double nEff=0;
  for(int i=1;i<=Ncdm;i++) nEff+=fracCDM[i]/McdmN[i];  
  r*=vSigma*nEff*nEff;
  return r/3.2E-28;  
}



int main(int argc,char** argv)
{ int err,n,i;

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
