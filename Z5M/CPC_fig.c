/*====== Modules ===============
   Keys to switch on
   various modules of micrOMEGAs
================================*/

#define MASSES_INFO
  /* Display information about mass spectrum  */

#define FREEZEIN  /*  Calculate relic density in Freeze-in scenario  */

#define DECAYS

#define CLEAN
/*===== End of DEFINE  settings ===== */


#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"
#include"lib/pmodel.h"


int main(int argc,char** argv)
{  int err;
   char cdmName[10];
   int spin2, charge3,cdim;

  ForceUG=0;  /* to Force Unitary Gauge assign 1 */
  //useSLHAwidth=0;
  VZdecay=0; VWdecay=0; cleanDecayTable();
  loadHeffGeff("GG.thg");

  if(argc==1)
  {
      printf(" Correct usage:  ./main  <file with parameters> \n");
      printf("Example: ./main data1.par\n");
      exit(1);
  };

  err=readVar(argv[1]);

  if(err==-1)     {printf("Can not open the file\n"); exit(1);}
  else if(err>0)  { printf("Wrong file contents at line %d\n",err);exit(1);}



  err=sortOddParticles(cdmName);
  if(err) { printf("Can't calculate %s\n",cdmName); return 1;}

  if(CDM[1])
  {
     qNumbers(CDM[1], &spin2, &charge3, &cdim);
     printf("\nDark matter candidate is '%s' with spin=%d/2 mass=%.2E\n",CDM[1],  spin2,McdmN[1]);
     if(charge3) printf("Dark Matter has electric charge %d/3\n",charge3);
     if(cdim!=1) printf("Dark Matter is a color particle\n");
  }
  if(CDM[2])
  {
     qNumbers(CDM[2], &spin2, &charge3, &cdim);
     printf("\nDark matter candidate is '%s' with spin=%d/2 mass=%.2E\n",CDM[2],spin2,McdmN[2]);
     if(charge3) printf("Dark Matter has electric charge %d/3\n",charge3);
     if(cdim!=1) printf("Dark Matter is a color particle\n");
  }


#ifdef MASSES_INFO
{
  printf("\n=== MASSES OF HIGGS AND ODD PARTICLES: ===\n");
  printHiggs(stdout);
  printMasses(stdout,1);
}
#endif


#ifdef FREEZEIN
{
  double TR=1E4;
  double omegaFi;
  int fast=0;
  double Y[2];

VZdecay=0; VWdecay=0; cleanDecayTable();

 double Xf,Omega;
   double Beps=1.E-6, cut=0.01;
   int i,err;
   printf("\n==== Calculation of relic density =====\n");

   Y[0]=0;
   Y[1]=YdmNEq(TR,"2");
   Omega=darkOmegaNTR(TR,Y,fast,Beps,&err);
   printf("Omega1 (N)=%.2E\n",Omega*fracCDM[1]);
   printf("Omega2 (N)=%.2E\n",Omega*fracCDM[2]);

   toFeebleList(CDM[1]);
/*   Omega=darkOmega(&Xf,fast,Beps,&err);
   printf("Omega2 (one-component)  =%.2E\n", Omega);
*/
  Omega=darkOmegaFi(TR,CDM[1],&err);
  printf("Omega_FI=%.2E\n",Omega);
  printChannelsFi(0,0,stdout);

   displayPlot("Y ","T",   Tend,TR,1,3
      , "Y1",0,YdmN, "1"
      , "Y2",0,YdmN, "2"
  //    , "Yom1",0,YF,NULL
      , "YFi",0,YFi,NULL
    );

  }
#endif





#ifdef CLEAN
  system("rm -f HB.* HB.* hb.* hs.*  debug_channels.txt debug_predratio.txt  Key.dat");
  system("rm -f Lilith_*   particles.py*");
//  system("rm -f   smodels.*");
#endif



  killPlots();
  return 0;
}
