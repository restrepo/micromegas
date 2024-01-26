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
  
#define OMEGA       /*  Calculate Freeze out relic density and display contribution of  individual channels */

#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"
#include"lib/pmodel.h"


int main(int argc,char** argv)
{  int err;
   char cdmName[10];
   int spin2, charge3,cdim;
  ForceUG=0;  /* to Force Unitary Gauge assign 1 */
  //useSLHAwidth=0;
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

  
#ifdef MASSES_INFO
{
  printf("\n=== MASSES OF HIGGS AND ODD PARTICLES: ===\n");
  printHiggs(stdout);
  printMasses(stdout,1);
}
#endif



#ifdef OMEGA
{ int fast=0;
  double Beps=1.E-6, cut=0.01;
  double Omega;  
  int i,err; 
   VZdecay=0; VWdecay=0; cleanDecayTable();
  printf("\n==== Calculation of relic density =====\n");   

     sortOddParticles(NULL);      
     printThermalSets();        //! new
     Tend=1E-3;
     Omega=darkOmega2(fast,Beps,&err);  
     printf("Omega2=%.4e\n",Omega);
     double Y[2];    
     err=darkOmegaN(Y,Beps,&err);    //! new
printf("err=%d\n", err);     
     if(err==0)for(int i=1;i<=Ncdm;i++)  printf(" Omega[%d]=%.3E ", i, Y[i-1]*2.742E8*McdmN[i]);
printf("\n");        

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
