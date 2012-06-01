#include"../sources/micromegas.h"

int main(int narg, char** args)
{
  double xiP=1,xiN=1; /* WIMP-proton and WIMP-neutron SpinDependent couplings */   
/* Testing SD form factors */
  PlotSS(S11Ge73,73, "S11Ge73" ,50);
/* 
  PlotSS(S01Nb93,93, "S01Nb93" ,150);
  PlotSS(S11Nb93,93, "S11Nb93" ,150);
*/
  PlotSS(S11Ge73A,73, "S11Ge73A" ,150);
/*  
  PlotSS(S00F19,19, "S01F19" ,150);
  PlotSS(S11F19,19, "S11F19" ,150);
*/

/* testing  SD wimp-nuclei  scattering */

/*  
  Plot3SS(xiP,xiN,S00Si29,S01Si29,S11Si29,29,0.5, "Si29", 15);    

  Plot3SS0(xiP,xiN,Sp_Si29,Sn_Si29,29,0.5, "Si29-Gauss", 15);

  Plot3SS(xiP,xiN,S00I127,S01I127,S11I127,127,J_I127, "J127", 15);    

  Plot3SS0(xiP,xiN,Sp_I127,Sn_I127,127,J_I127, "I127-Gauss", 15);

  Plot3SS(xiP,xiN,S00Nb93,S01Nb93,S11Nb93,93,J_Nb93, "Nb93", 150);    

  Plot3SS0(xiP,xiN,Sp_Nb93,Sn_Nb93,       93,J_Nb93, "Nb93-Gauss", 150);
*/

  Plot3SS(xiP,xiN,S00Xe131,S01Xe131,S11Xe131,131,J_Xe131, "Xe131", 15);    

  Plot3SS0(xiP,xiN,Sp_Xe131,Sn_Xe131,131,J_Xe131, "Xe131-Gauss", 15);
  
  

}
