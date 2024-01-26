#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"

static double hEffIdeal(double T) 
{
 
  double sum=  2+16+6*h1eff(80/T,-1)+ 3*h1eff(90/T,-1)          // bosons
              +7./8.*2*3                                        // neutrinos 
              +7./8.*2*2*2 + 4*h1eff(1.7/T,1)                   // charged leptons 
              +12*( 3*7./8.+ h1eff(1.4/T,1) +h1eff(5/T,1) + h1eff(175/T,1)) // quarks 
            ;   
  return sum;      
}

int main(int argc,char** argv)
{
  char * gh[7] = {"../sources/hgEff/std_thg.tab", "../sources/hgEff/std_thg_old.tab",   "../sources/hgEff/HP_A_thg.tab", "../sources/hgEff/HP_B_thg.tab","../sources/hgEff/HP_B2_thg.tab","../sources/hgEff/HP_B3_thg.tab","../sources/hgEff/HP_C_thg.tab"};  

  for(int i=1;i<7;i++) 
  { loadHeffGeff(gh[i]);
    displayPlot(gh[i],"T",0.01, 1000, 1,4 ,"gEff",0, gEff,NULL
                                          ,"geff2",0, gEff2,NULL
                                          ,"hEff",0, hEff,NULL
                                          ,"hEffIdeal",0, hEffIdeal  ); 
//     displayPlot(gh[i],"T",0.001, 10,0,1, "dlogh/dlogT", 0, hEffLnDiff,NULL);
  }

  double Year= 31558149.8; // is seconds
  printf(" Age of Univerce=%e\n", HubbleTime(0.1,T2_73K)/Year);
  double T_CMB=3000./300./38.681740*1E-9; //GeV
  printf(" Age of CMB Time=%E\n", HubbleTime(10,T_CMB)/Year);
  printf(" teq=%.2E yr\n", HubbleTime(10,Tmreq )/Year);
  double T1=1E-4, pm1=0.5;
  printf("Free Streaming (p/m=%.2, T1=%.2) = %.2E[Mpc]\n",pm1,T1, freeStreaming(pm1,T1, Tmreq));
  return 0;
}


