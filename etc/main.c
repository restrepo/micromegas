
/*====== Modules ===============
   Keys to switch on
  various modules of micrOMEGAs
================================*/

#define MASSES_INFO      
      /* Display information about masses of odd sector  
      */
#define CONSTRAINTS     
      /* This module is empty yet! 
      */ 
#define OMEGA            
      /* Calculate relic density and display contribution of
         individual channels 
      */
/*#define INDIRECT_DETECTION  */
      /* Compute spectra of gamma/positron/neutrinos
         for DM annihilation; calculate <sigma*v> and
         integrate gamma signal over DM galactic squared
         density for given line of sight.  
      */
/*#define RESET_FORMFACTORS*/
      /* Modify default nucleus form factors, 
         DM velocity distribution,
         A-dependence of Fermi-dencity
      */     
#define CDM_NUCLEON    
      /* Calculate amplitudes and cross-sections for 
         CDM-mucleon collisions 
      */  
/*#define CDM_NUCLEUS      */
      /* Calculate number of events for 1kg*day 
         and recoil energy distibution for various nuclei
      */
#define CROSS_SECTIONS 
      /* Calculate cross sections and widths for 
         reactions specified by the user
      */
/*===== end of Modules  ======*/

/*===== Options ========*/
/*#define SHOWPLOTS*/
     /* Display  graphical plots on the screen */ 

/*===== End of DEFINE  settings ===== */

#include"../sources/micromegas.h"
#include"lib/pmodel.h"

double dN_th(double z)
{  double x=exp(z);
   if(4*Mcdm*Mcdm*(1-x)/91.2/91.2 < exp(1)) return 0;
   return 1./137./M_PI*(2-2*x)*(log(4*Mcdm*Mcdm*(1-x)/91.2/91.2)-1);
}

#define LnXmin (-10*M_LN10)
static double zInterp_old(double zz, double * tab)  /* zz=log(E/Nmass) */
{  
   double z,dz,r;
   int j0;   
   
   if(zz>0) return 0;
   z=zz/LnXmin;  /*  z_j=0.5/NZ +j/NZ  */
//   if(NZ*z-0.5 <0) return tab[0];
   j0=NZ*z-0.5;
   if(j0<0) j0=0;
   if(j0>=NZ) return tab[NZ-1];
   dz=z*NZ-0.5-j0;
printf("z=%E j0=%d\n",zz,j0);   
   r=(1-dz)*tab[j0]+dz*tab[j0+1];
   if(r<0)r=0;
   return r; 
}

static double WA_rest(double z)
{
  double SpA_rest[NZ]={3.800E-04,1.194E-03,2.182E-03,4.528E-03,6.645E-03,9.793E-03,1.468E-02,2.152E-02,3.135E-02,4.322E-02,5.921E-02,7.913E-02,1.036E-01,1.355E-01,1.756E-01,2.185E-01,2.709E-01,3.330E-01,4.002E-01,4.806E-01,5.699E-01,6.726E-01,7.847E-01,9.048E-01,1.037E+00,1.187E+00,1.350E+00,1.522E+00,1.709E+00,1.905E+00,2.095E+00,2.328E+00,2.545E+00,2.770E+00,3.012E+00,3.269E+00,3.513E+00,3.792E+00,4.037E+00,4.300E+00,4.574E+00,4.822E+00,5.073E+00,5.324E+00,5.565E+00,5.798E+00,6.017E+00,6.227E+00,6.427E+00,6.598E+00,6.759E+00,6.891E+00,7.020E+00,7.131E+00,7.213E+00,7.241E+00,7.270E+00,7.270E+00,7.254E+00,7.211E+00,7.140E+00,7.044E+00,6.944E+00,6.822E+00,6.658E+00,6.489E+00,6.311E+00,6.114E+00,5.898E+00,5.672E+00,5.457E+00,5.221E+00,5.003E+00,4.778E+00,4.512E+00,4.283E+00,4.050E+00,3.822E+00,3.596E+00,3.373E+00,3.167E+00,2.959E+00,2.759E+00,2.584E+00,2.396E+00,2.226E+00,2.058E+00,1.903E+00,1.765E+00,1.630E+00,1.503E+00,1.389E+00,1.275E+00,1.179E+00,1.080E+00,9.881E-01,9.112E-01,8.392E-01,7.706E-01,7.069E-01,6.487E-01,5.974E-01,5.416E-01,5.015E-01,4.602E-01,4.211E-01,3.832E-01,3.506E-01,3.235E-01,2.965E-01,2.724E-01,2.522E-01,2.299E-01,2.102E-01,1.961E-01,1.810E-01,1.627E-01,1.526E-01,1.398E-01,1.303E-01,1.175E-01,1.094E-01,1.007E-01,9.311E-02,8.789E-02,7.980E-02,7.608E-02,7.083E-02,6.498E-02,5.894E-02,5.626E-02,5.131E-02,4.928E-02,4.651E-02,4.168E-02,4.064E-02,3.851E-02,3.499E-02,3.397E-02,3.170E-02,3.026E-02,2.934E-02,2.698E-02,2.603E-02,2.389E-02,2.408E-02,2.251E-02,2.050E-02,1.948E-02,1.911E-02,1.811E-02,1.808E-02,1.689E-02,1.655E-02,1.630E-02,1.523E-02,1.467E-02,1.334E-02,1.393E-02,1.314E-02,1.306E-02,1.278E-02,1.170E-02,1.166E-02,1.161E-02,1.147E-02,1.062E-02,1.033E-02,9.630E-03,1.005E-02,9.544E-03,8.957E-03,9.012E-03,8.708E-03,8.838E-03,8.838E-03,8.034E-03,8.089E-03,7.937E-03,7.513E-03,7.470E-03,7.231E-03,7.665E-03,7.198E-03,7.079E-03,6.775E-03,6.862E-03,6.894E-03,6.048E-03,6.558E-03,6.699E-03,5.863E-03,6.004E-03,5.809E-03,5.678E-03,5.505E-03,5.461E-03,5.830E-03,4.929E-03,5.461E-03,4.799E-03,5.049E-03,4.690E-03,4.918E-03,4.734E-03,5.005E-03,4.614E-03,4.647E-03,4.300E-03,4.158E-03,4.093E-03,3.876E-03,3.930E-03,4.072E-03,3.453E-03,3.757E-03,3.572E-03,3.268E-03,3.659E-03,3.659E-03,3.225E-03,3.235E-03,3.094E-03,3.073E-03,2.790E-03,2.801E-03,3.040E-03,2.866E-03,2.736E-03,2.660E-03,2.399E-03,2.617E-03,2.747E-03,2.638E-03,2.323E-03,2.660E-03,2.443E-03,2.063E-03,2.095E-03,2.074E-03,2.052E-03,2.063E-03,1.857E-03,1.933E-03,2.052E-03,1.759E-03,1.748E-03,1.759E-03,1.629E-03,1.639E-03};
  return zInterp_old(z,SpA_rest); 
}
static double WA_full(double z)
{ return WA_rest(z)+ dN_th(z);}
 
static double WA_pyth(double z)
{
  double SpA_pyth[NZ]={4.017E-04,1.455E-03,2.410E-03,4.093E-03,6.884E-03,1.028E-02,1.489E-02,2.264E-02,3.211E-02,4.311E-02,5.947E-02,7.892E-02,1.053E-01,1.338E-01,1.728E-01,2.156E-01,2.688E-01,3.318E-01,4.010E-01,4.827E-01,5.726E-01,6.722E-01,7.817E-01,9.044E-01,1.035E+00,1.190E+00,1.347E+00,1.509E+00,1.702E+00,1.893E+00,2.102E+00,2.312E+00,2.539E+00,2.774E+00,3.014E+00,3.263E+00,3.523E+00,3.786E+00,4.050E+00,4.307E+00,4.559E+00,4.826E+00,5.071E+00,5.321E+00,5.556E+00,5.791E+00,6.016E+00,6.231E+00,6.409E+00,6.597E+00,6.756E+00,6.901E+00,7.023E+00,7.122E+00,7.178E+00,7.240E+00,7.273E+00,7.275E+00,7.258E+00,7.215E+00,7.147E+00,7.045E+00,6.946E+00,6.818E+00,6.656E+00,6.493E+00,6.308E+00,6.116E+00,5.897E+00,5.681E+00,5.465E+00,5.238E+00,4.999E+00,4.758E+00,4.520E+00,4.285E+00,4.049E+00,3.829E+00,3.595E+00,3.382E+00,3.165E+00,2.960E+00,2.760E+00,2.564E+00,2.395E+00,2.220E+00,2.057E+00,1.915E+00,1.768E+00,1.633E+00,1.511E+00,1.390E+00,1.284E+00,1.177E+00,1.084E+00,9.931E-01,9.156E-01,8.394E-01,7.714E-01,7.092E-01,6.482E-01,5.943E-01,5.462E-01,4.987E-01,4.613E-01,4.217E-01,3.818E-01,3.508E-01,3.250E-01,3.003E-01,2.714E-01,2.517E-01,2.319E-01,2.121E-01,1.942E-01,1.796E-01,1.638E-01,1.502E-01,1.389E-01,1.307E-01,1.174E-01,1.116E-01,1.006E-01,9.493E-02,8.663E-02,7.951E-02,7.336E-02,6.851E-02,6.441E-02,5.919E-02,5.641E-02,5.309E-02,4.912E-02,4.539E-02,4.253E-02,3.948E-02,3.764E-02,3.549E-02,3.404E-02,3.190E-02,3.047E-02,2.866E-02,2.708E-02,2.528E-02,2.467E-02,2.250E-02,2.191E-02,2.138E-02,1.942E-02,1.967E-02,1.809E-02,1.667E-02,1.717E-02,1.635E-02,1.506E-02,1.500E-02,1.491E-02,1.438E-02,1.380E-02,1.278E-02,1.289E-02,1.274E-02,1.254E-02,1.175E-02,1.116E-02,1.117E-02,1.068E-02,1.033E-02,1.008E-02,9.663E-03,9.457E-03,9.554E-03,9.131E-03,9.153E-03,8.936E-03,8.838E-03,8.469E-03,8.338E-03,8.219E-03,8.013E-03,8.067E-03,7.416E-03,6.992E-03,7.242E-03,6.525E-03,6.666E-03,6.797E-03,6.656E-03,6.297E-03,6.514E-03,6.590E-03,5.950E-03,6.113E-03,5.961E-03,5.852E-03,5.559E-03,5.928E-03,5.537E-03,5.331E-03,5.494E-03,5.581E-03,4.886E-03,4.897E-03,4.462E-03,4.625E-03,4.690E-03,4.191E-03,4.528E-03,4.766E-03,4.289E-03,4.039E-03,4.028E-03,4.115E-03,3.844E-03,3.778E-03,3.583E-03,3.626E-03,3.409E-03,3.485E-03,3.442E-03,2.953E-03,2.997E-03,3.083E-03,2.573E-03,2.964E-03,3.018E-03,3.094E-03,2.790E-03,2.703E-03,2.530E-03,2.638E-03,2.910E-03,2.530E-03,2.269E-03,2.117E-03,2.074E-03,1.976E-03,2.030E-03,2.171E-03,1.911E-03,1.846E-03,2.041E-03,1.954E-03,1.629E-03,1.705E-03,1.661E-03,1.629E-03,1.672E-03,1.748E-03,1.650E-03};
  return zInterp_old(z,SpA_pyth);
}

int main(int argc,char** argv)
{  int err;
   char cdmName[10];
   int spin2, charge3,cdim;

/* to save RGE input/output files uncomment the next line */
/*delFiles(0);*/

  if(argc==1)
  { 
      printf(" Correct usage:  ./omg <file with parameters> \n");
      exit(1);
  }
                               

  err=readVar(argv[1]);
/*   err=readVarRHNM(argv[1]);*/
  if(err==-1)     {printf("Can not open the file\n"); exit(1);}
  else if(err>0)  { printf("Wrong file contents at line %d\n",err);exit(1);}


  err=sortOddParticles(cdmName);
  if(err) { printf("Can't calculate %s\n",cdmName); return 1;}

/*to print input parameters or model in SLHA format uncomment correspondingly*/
/* 
  printVar(stdout);  
  writeLesH("slha.out"); 
*/

#ifdef MASSES_INFO
{
  printf("\n=== MASSES OF PARTICLES OF ODD SECTOR: ===\n");
  printMasses(stdout,1);
}
#endif

#ifdef CONSTRAINTS
  printf("\n================= CONSTRAINTS =================\n");
#endif

#ifdef OMEGA
{ int fast=1;
  double Beps=1.E-2, cut=0.01;
  double Omega,Xf;   
  printf("\n==== Calculation of relic density =====\n");  

  Omega=darkOmega(&Xf,fast,Beps);
  printf("Xf=%.2e Omega=%.2e\n",Xf,Omega);
  printChannels(Xf,cut,Beps,1,stdout);
}
#endif


#ifdef INDIRECT_DETECTION
{ 
  int err,i;
  double Emin,Ntot,Xtot,Xsum=0,sigmaV;
  char txt[100];
  double SpA[NZ],SpE[NZ],SpP[NZ];
  double FluxA[NZ],FluxE[NZ],FluxP[NZ];
/*   double *SpA=NULL, *SpE=NULL; */
  double * SpNe=NULL,*SpNm=NULL,*SpNl=NULL;
/* double SpNe[NZ],SpNm[NZ],SpNl[NZ]; */

printf("\n==== Indirect detection =======\n");  

  Emin=1;  /* Energy cut  in GeV   */

  sigmaV=calcSpectrum(2,SpA,SpE,SpP,SpNe,SpNm,SpNl ,&err);  
             /* Returns sigma*v in cm^3/sec.  
                SpX - calculated spectra of annihilation.
                Use SpectradNdE(E, SpX) to calculate energy distribution
                in  1/GeV units.

             */
  printf("sigmav=%E[cm^3/s]\n",sigmaV);  

  if(SpA)
  {  double fi=0.1,dfi=0.05; /* angle of sight in radians */
     double HF;
     spectrInfo(Emin/Mcdm,SpA, &Ntot,&Xtot);
     printf(" Initial gamma spectrum: Ntot=%.2E Energy fraction =%.2E\n",Ntot, Xtot/2);
     displaySpectrum(SpA,"initial Gamma",Emin,Mcdm,1);

     Xsum+=Xtot/2;
     gammaFluxTab(fi,dfi, sigmaV, SpA, FluxA);     
#ifdef SHOWPLOTS
     sprintf(txt,"photon  flux[1/cm^2/s/GeV] at f=%.2f rad",fi);
     displaySpectrum(FluxA,txt,Emin,Mcdm,1); 
#endif     
  }

  if(SpE)
  { extern double kernel_test(double);
    spectrInfo(Emin/Mcdm,SpE, &Ntot,&Xtot);
    printf(" Initial positron spectrum: Ntot=%.2E Energy fraction =%.2E\n",Ntot, Xtot/2);
    displaySpectrum(SpE,"initial electrons",Emin,Mcdm,0); 
    Xsum+=Xtot;
    posiFluxTab(Emin, sigmaV, SpE, FluxE);
#ifdef SHOWPLOTS     
    displaySpectrum(FluxE,"positron flux [1/cm^2/s/sr/GeV]" ,Emin,Mcdm,1);
#endif
  }
  
  if(SpP)
  { 
    spectrInfo(Emin/Mcdm,SpP, &Ntot,&Xtot);
    printf(" Initial a-protons spectrum: Ntot=%.2E Energy fraction =%.2E\n",Ntot, Xtot/2);
    Xsum+=Xtot;
displaySpectrum(SpP,"initial antiproton flux " ,Emin,Mcdm,1);    
    pbarFluxTab(Emin, sigmaV, SpP, FluxP  ); 
#ifdef SHOWPLOTS    
     displaySpectrum(FluxP,"antiproton flux [1/cm^2/s/sr/GeV]" ,Emin,Mcdm,1);
#endif  
  }
  
  if(SpNe)
  { 
    spectrInfo(Emin/Mcdm,SpNe, &Ntot,&Xtot);
    printf(" Initial neutrino-e spectrum: Ntot=%.2E Energy fraction =%.2E\n",Ntot, Xtot/2);
    Xsum+=Xtot;
  }
  if(SpNm)
  { 
    spectrInfo(Emin/Mcdm,SpNm, &Ntot,&Xtot);
    printf(" Initial neutrino-e spectrum: Ntot=%.2E Energy fraction =%.2E\n",Ntot, Xtot/2);
    Xsum+=Xtot;
  }
  if(SpNl)
  { 
    spectrInfo(Emin/Mcdm,SpNl, &Ntot,&Xtot);
    printf(" Initial neutrino-e spectrum: Ntot=%.2E Energy fraction =%.2E\n",Ntot, Xtot/2);
    Xsum+=Xtot;
  }
  
  if(SpA && SpE && SpP && SpNe && SpNm && SpNl) printf("Check of energy conservation:"
            " Xsum=%.2E (sould be close to 1) \n",Xsum); 

}
#endif

#ifdef RESET_FORMFACTORS
{
/* 
   The default nucleon form factors can be completely or partially modified 
   by setProtonFF and setNeutronFF. For scalar form factors, one can first call
   getScalarFF( Mu/Md, Ms/Md, sigmaPiN[MeV], sigma0[MeV], protonFF,neutronFF)  
   or set the new coefficients by directly assigning numerical values.
*/
{ double   ffS0P[3]={0.033,0.023,0.26},
           ffS0N[3]={0.042,0.018,0.26},
           ffV5P[3]={-0.427, 0.842,-0.085},
           ffV5N[3]={ 0.842,-0.427,-0.085}; 

  printf("\n=========== Redefinition of form factors  =========\n");         
      
  getScalarFF(0.553,18.9,55.,35.,ffS0P, ffS0N);
  printf("protonFF  d %E, u %E, s %E\n",ffS0P[0],ffS0P[1],ffS0P[2]);                               
  printf("neutronFF d %E, u %E, s %E\n",ffS0N[0],ffS0N[1],ffS0N[2]);

/* Use NULL argument if there is no need for reassignment */
  setProtonFF(ffS0P,ffV5P, NULL);
  setNeutronFF(ffS0N,ffV5N,NULL);
}

/* Option to change parameters of DM velocity  distribution 
*/   
   SetfMaxwell(220.,244.4,600.);
     /* arg1- defines DM velocity distribution in Galaxy rest frame:
            ~exp(-v^2/arg1^2)d^3v
        arg2- Earth velocity with respect to Galaxy
        arg3- Maximal DM velocity in Sun orbit with respect to Galaxy.
        All parameters are  in [km/s] units.
     */
/* In case DM has velocity distribution close to delta-function 
   the DM velocity V[km/s] can be defined by
*/          
   SetfDelta(350.);

/* To reset parameters of Fermi nucleus distribution  */
   SetFermi(1.23,-0.6,0.52);
/*  with half-density radius for Fermi distribution: 
          c=arg1*A^(1/3) + arg2
    and arg3 is the surface thickness.
    All parameter in [fm].      
*/
}
#endif


#ifdef CDM_NUCLEON
{ double pA0[2],pA5[2],nA0[2],nA5[2];
  double Nmass=0.939; /*nucleon mass*/
  double SCcoeff;        
  double dpA0[2],dnA0[2];
 
printf("\n==== Calculation of CDM-nucleons amplitudes  =====\n");   

  nucleonAmplitudes(NULL, dpA0,pA5,dnA0,nA5);
printf("====OFF/On======\n");  
  nucleonAmplitudes(NULL, pA0,pA5,nA0,nA5);
  dpA0[0]-=pA0[0];
  dnA0[0]-=nA0[0];  
   
    printf("%s -nucleon amplitudes:\n",cdmName);
    printf("proton:  SI  %.3E  SD  %.3E\n",pA0[0],pA5[0]);
    printf("neutron: SI  %.3E  SD  %.3E\n",nA0[0],nA5[0]); 

  SCcoeff=4/M_PI*3.8937966E8*pow(Nmass*Mcdm/(Nmass+ Mcdm),2.);
    printf("%s-nucleon cross sections:\n",cdmName);
    
    printf(" proton  SI %.3E  SD %.3E\n",SCcoeff*pA0[0]*pA0[0],3*SCcoeff*pA5[0]*pA5[0]);
    printf(" neutron SI %.3E  SD %.3E\n",SCcoeff*nA0[0]*nA0[0],3*SCcoeff*nA5[0]*nA5[0]);

 printf(" twist-2 CS proton   %.3E   neutron %.3E \n",
 SCcoeff*dpA0[0]*dpA0[0], SCcoeff*dnA0[0]*dnA0[0]);
 
    printf("anti-%s -nucleon amplitudes:\n",cdmName);
    printf("proton:  SI  %.3E  SD  %.3E\n",pA0[1],pA5[1]);
    printf("neutron: SI  %.3E  SD  %.3E\n",nA0[1],nA5[1]); 

  SCcoeff=4/M_PI*3.8937966E8*pow(Nmass*Mcdm/(Nmass+ Mcdm),2.);
    printf("anti-%s-nucleon cross sections:\n",cdmName);
    
    printf(" proton  SI %.3E  SD %.3E\n",SCcoeff*pA0[1]*pA0[1],3*SCcoeff*pA5[1]*pA5[1]);
    printf(" neutron SI %.3E  SD %.3E\n",SCcoeff*nA0[1]*nA0[1],3*SCcoeff*nA5[1]*nA5[1]);

}
#endif
  
#ifdef CDM_NUCLEUS
{ double dNdE[200];
  double nEvents;
  double rho=0.3; /* DM density GeV/sm^3 */
printf("\n=========== Direct Detection ===============\n");


  nEvents=nucleusRecoil(rho,fDvMaxwell,73,Z_Ge,J_Ge73,S00Ge73,S01Ge73,S11Ge73,NULL,dNdE);
      /* See '../sources/micromegas.h' for description of arguments 
     
        Instead of Maxwell (DvMaxwell) one can use 'fDvDelta' Delta-function 
        velocity distribution.
      */

  printf("73Ge: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));                                   
#ifdef SHOWPLOTS
    displayRecoilPlot(dNdE,"Distribution of recoil energy of 73Ge",0,199);
#endif

  nEvents=nucleusRecoil(rho,fDvMaxwell,131,Z_Xe,J_Xe131,S00Xe131,S01Xe131,S11Xe131,NULL,dNdE);

  printf("131Xe: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));                                   
#ifdef SHOWPLOTS
    displayRecoilPlot(dNdE,"Distribution of recoil energy of 131Xe",0,199);
#endif

}
#endif 

#ifdef CROSS_SECTIONS
{
  double Pcm=500;
  numout* cc;
  double cosmin=-0.99, cosmax=0.99;
  double v=0.002;

printf("\n====== Calculation of widths and cross sections ====\n");  
  decay2Info("Z",stdout);
  decay2Info("H",stdout);

/*  Helicity[0]=0.45;
  Helicity[1]=-0.45;
  printf("Process e,E->2*x at Pcm=%.3E GeV\n",Pcm);
  cc=newProcess("e%,E%->2*x","eE_2x");
  if(cc)
  { int ntot,l;
    char * name[4];
    procInfo1(cc,&ntot,NULL,NULL);
    for(l=1;l<=ntot; l++)
    { int err;
      double cs;
      procInfo2(cc,l,name,NULL);
      printf("%3s,%3s -> %3s %3s  ",name[0],name[1],name[2],name[3]);
      cs= cs22(cc,l,Pcm,cosmin,cosmax,&err);
      if(err) printf("Error\n");
      else if(cs==0.) printf("Zero\n");
      else printf("%.2E [pb]\n",cs); 
    }
  } 
*/
  printf("\n CDM annihilation at V_rel=%.2E\n",v);
 
  cc=newProcess("",wimpAnnLib());
  assignValW("Q",2*Mcdm);
  if(cc)
  { int ntot,l;
    char * name[4];
    double mass[4];
    procInfo1(cc,&ntot,NULL,NULL);
    for(l=1;l<=ntot; l++)
    { int err;
      double cs;
      procInfo2(cc,l,name,mass);
      if(l==1) { Pcm=mass[0]*v/2; printf("(Pcm=%.2E)\n",Pcm);}
      printf("%3s,%3s -> %3s %3s  ",name[0],name[1],name[2],name[3]);
      cs= cs22(cc,l,Pcm,-1.,1.,&err);
      if(err) printf("Error\n");
      else if(cs==0.) printf("Zero\n");
      else printf("%.2E [pb] ( sigma*v=%.2E [cm^3/sec] )  \n",cs,cs*v*2.9979E-26); 
    }
  }
}

#endif
{ double SpA_tree[NZ];
  double SpA_rest[NZ]={3.800E-04,1.194E-03,2.182E-03,4.528E-03,6.645E-03,9.793E-03,1.468E-02,2.152E-02,3.135E-02,4.322E-02,5.921E-02,7.913E-02,1.036E-01,1.355E-01,1.756E-01,2.185E-01,2.709E-01,3.330E-01,4.002E-01,4.806E-01,5.699E-01,6.726E-01,7.847E-01,9.048E-01,1.037E+00,1.187E+00,1.350E+00,1.522E+00,1.709E+00,1.905E+00,2.095E+00,2.328E+00,2.545E+00,2.770E+00,3.012E+00,3.269E+00,3.513E+00,3.792E+00,4.037E+00,4.300E+00,4.574E+00,4.822E+00,5.073E+00,5.324E+00,5.565E+00,5.798E+00,6.017E+00,6.227E+00,6.427E+00,6.598E+00,6.759E+00,6.891E+00,7.020E+00,7.131E+00,7.213E+00,7.241E+00,7.270E+00,7.270E+00,7.254E+00,7.211E+00,7.140E+00,7.044E+00,6.944E+00,6.822E+00,6.658E+00,6.489E+00,6.311E+00,6.114E+00,5.898E+00,5.672E+00,5.457E+00,5.221E+00,5.003E+00,4.778E+00,4.512E+00,4.283E+00,4.050E+00,3.822E+00,3.596E+00,3.373E+00,3.167E+00,2.959E+00,2.759E+00,2.584E+00,2.396E+00,2.226E+00,2.058E+00,1.903E+00,1.765E+00,1.630E+00,1.503E+00,1.389E+00,1.275E+00,1.179E+00,1.080E+00,9.881E-01,9.112E-01,8.392E-01,7.706E-01,7.069E-01,6.487E-01,5.974E-01,5.416E-01,5.015E-01,4.602E-01,4.211E-01,3.832E-01,3.506E-01,3.235E-01,2.965E-01,2.724E-01,2.522E-01,2.299E-01,2.102E-01,1.961E-01,1.810E-01,1.627E-01,1.526E-01,1.398E-01,1.303E-01,1.175E-01,1.094E-01,1.007E-01,9.311E-02,8.789E-02,7.980E-02,7.608E-02,7.083E-02,6.498E-02,5.894E-02,5.626E-02,5.131E-02,4.928E-02,4.651E-02,4.168E-02,4.064E-02,3.851E-02,3.499E-02,3.397E-02,3.170E-02,3.026E-02,2.934E-02,2.698E-02,2.603E-02,2.389E-02,2.408E-02,2.251E-02,2.050E-02,1.948E-02,1.911E-02,1.811E-02,1.808E-02,1.689E-02,1.655E-02,1.630E-02,1.523E-02,1.467E-02,1.334E-02,1.393E-02,1.314E-02,1.306E-02,1.278E-02,1.170E-02,1.166E-02,1.161E-02,1.147E-02,1.062E-02,1.033E-02,9.630E-03,1.005E-02,9.544E-03,8.957E-03,9.012E-03,8.708E-03,8.838E-03,8.838E-03,8.034E-03,8.089E-03,7.937E-03,7.513E-03,7.470E-03,7.231E-03,7.665E-03,7.198E-03,7.079E-03,6.775E-03,6.862E-03,6.894E-03,6.048E-03,6.558E-03,6.699E-03,5.863E-03,6.004E-03,5.809E-03,5.678E-03,5.505E-03,5.461E-03,5.830E-03,4.929E-03,5.461E-03,4.799E-03,5.049E-03,4.690E-03,4.918E-03,4.734E-03,5.005E-03,4.614E-03,4.647E-03,4.300E-03,4.158E-03,4.093E-03,3.876E-03,3.930E-03,4.072E-03,3.453E-03,3.757E-03,3.572E-03,3.268E-03,3.659E-03,3.659E-03,3.225E-03,3.235E-03,3.094E-03,3.073E-03,2.790E-03,2.801E-03,3.040E-03,2.866E-03,2.736E-03,2.660E-03,2.399E-03,2.617E-03,2.747E-03,2.638E-03,2.323E-03,2.660E-03,2.443E-03,2.063E-03,2.095E-03,2.074E-03,2.052E-03,2.063E-03,1.857E-03,1.933E-03,2.052E-03,1.759E-03,1.748E-03,1.759E-03,1.629E-03,1.639E-03};
  double SpA_pyth[NZ]={4.017E-04,1.455E-03,2.410E-03,4.093E-03,6.884E-03,1.028E-02,1.489E-02,2.264E-02,3.211E-02,4.311E-02,5.947E-02,7.892E-02,1.053E-01,1.338E-01,1.728E-01,2.156E-01,2.688E-01,3.318E-01,4.010E-01,4.827E-01,5.726E-01,6.722E-01,7.817E-01,9.044E-01,1.035E+00,1.190E+00,1.347E+00,1.509E+00,1.702E+00,1.893E+00,2.102E+00,2.312E+00,2.539E+00,2.774E+00,3.014E+00,3.263E+00,3.523E+00,3.786E+00,4.050E+00,4.307E+00,4.559E+00,4.826E+00,5.071E+00,5.321E+00,5.556E+00,5.791E+00,6.016E+00,6.231E+00,6.409E+00,6.597E+00,6.756E+00,6.901E+00,7.023E+00,7.122E+00,7.178E+00,7.240E+00,7.273E+00,7.275E+00,7.258E+00,7.215E+00,7.147E+00,7.045E+00,6.946E+00,6.818E+00,6.656E+00,6.493E+00,6.308E+00,6.116E+00,5.897E+00,5.681E+00,5.465E+00,5.238E+00,4.999E+00,4.758E+00,4.520E+00,4.285E+00,4.049E+00,3.829E+00,3.595E+00,3.382E+00,3.165E+00,2.960E+00,2.760E+00,2.564E+00,2.395E+00,2.220E+00,2.057E+00,1.915E+00,1.768E+00,1.633E+00,1.511E+00,1.390E+00,1.284E+00,1.177E+00,1.084E+00,9.931E-01,9.156E-01,8.394E-01,7.714E-01,7.092E-01,6.482E-01,5.943E-01,5.462E-01,4.987E-01,4.613E-01,4.217E-01,3.818E-01,3.508E-01,3.250E-01,3.003E-01,2.714E-01,2.517E-01,2.319E-01,2.121E-01,1.942E-01,1.796E-01,1.638E-01,1.502E-01,1.389E-01,1.307E-01,1.174E-01,1.116E-01,1.006E-01,9.493E-02,8.663E-02,7.951E-02,7.336E-02,6.851E-02,6.441E-02,5.919E-02,5.641E-02,5.309E-02,4.912E-02,4.539E-02,4.253E-02,3.948E-02,3.764E-02,3.549E-02,3.404E-02,3.190E-02,3.047E-02,2.866E-02,2.708E-02,2.528E-02,2.467E-02,2.250E-02,2.191E-02,2.138E-02,1.942E-02,1.967E-02,1.809E-02,1.667E-02,1.717E-02,1.635E-02,1.506E-02,1.500E-02,1.491E-02,1.438E-02,1.380E-02,1.278E-02,1.289E-02,1.274E-02,1.254E-02,1.175E-02,1.116E-02,1.117E-02,1.068E-02,1.033E-02,1.008E-02,9.663E-03,9.457E-03,9.554E-03,9.131E-03,9.153E-03,8.936E-03,8.838E-03,8.469E-03,8.338E-03,8.219E-03,8.013E-03,8.067E-03,7.416E-03,6.992E-03,7.242E-03,6.525E-03,6.666E-03,6.797E-03,6.656E-03,6.297E-03,6.514E-03,6.590E-03,5.950E-03,6.113E-03,5.961E-03,5.852E-03,5.559E-03,5.928E-03,5.537E-03,5.331E-03,5.494E-03,5.581E-03,4.886E-03,4.897E-03,4.462E-03,4.625E-03,4.690E-03,4.191E-03,4.528E-03,4.766E-03,4.289E-03,4.039E-03,4.028E-03,4.115E-03,3.844E-03,3.778E-03,3.583E-03,3.626E-03,3.409E-03,3.485E-03,3.442E-03,2.953E-03,2.997E-03,3.083E-03,2.573E-03,2.964E-03,3.018E-03,3.094E-03,2.790E-03,2.703E-03,2.530E-03,2.638E-03,2.910E-03,2.530E-03,2.269E-03,2.117E-03,2.074E-03,1.976E-03,2.030E-03,2.171E-03,1.911E-03,1.846E-03,2.041E-03,1.954E-03,1.629E-03,1.705E-03,1.661E-03,1.629E-03,1.672E-03,1.748E-03,1.650E-03};


  Mcdm=5000;
  
  
  
//  basicSpectra(Mcdm, 24,0,SpA);
//  displaySpectrum(SpA_rest,"W_rest",1.E-1*Mcdm,Mcdm,  0);
//  displaySpectrum(SpA_pyth,"W_pyth",1.E-1*Mcdm,Mcdm,  0); 
//  displayFunc(dN_th,-2.,-0.001,"th");
  
//  basicSpectra(Mcdm, 24+'T',0,SpA);
//  displaySpectrum(SpA,"gammas from WT",1.,Mcdm,  1);

displayFunc(WA_rest,-10,0,"in rest");
displayFunc(WA_pyth,-10,0,"pyth");
//displayFunc(dN_th,  -0,5,"th");
displayFunc(WA_full,-10,0,"full");
}
  killPlots(); 
  return 0;
}

