#include"../sources/micromegas.h"

int main(void)
{ double dNdE[200]; 
  double Mwimp, csSIp,csSIn,csSDp,csSDn,Leff;

  printf("Number of events for different nuclei ( /kg/day)\n");

  Mwimp=50;  csSIp=0; csSIn=0;  csSDp=0.01; csSDn=0.01;
  
printf("Parameters:\n Mwimp=%.2E csSIp=%.2E, csSIn=%.2E,csSDp=%.2E,csSDp=%.2E\n",
                           csSIp, csSIn,csSDp,csSDn);


 SetfMaxwell(244.4,700.);

printf("H1    %.2E\n",nucleusRecoil0Aux(Maxwell,  1,Z_H,J_H1,Sp_H1,Sn_H1,                   csSIp,csSIn,csSDp,csSDn,dNdE)); 
printf("He3   %.2E\n",nucleusRecoil0Aux(Maxwell,  3,Z_He,J_He3,Sp_He3,Sn_He3,               csSIp,csSIn,csSDp,csSDn,dNdE)); 
printf("F19   %.2E\n",nucleusRecoilAux(Maxwell,  19,Z_F,J_F19,S00F19,S01F19,S11F19,         csSIp,csSIn,csSDp,csSDn,dNdE)); 
printf("Na23  %.2E\n",nucleusRecoilAux(Maxwell,  23,Z_Na,J_Na23,S00Na23,S01Na23,S11Na23,    csSIp,csSIn,csSDp,csSDn,dNdE)); 
printf("Al27  %.2E\n",nucleusRecoilAux(Maxwell,  27,Z_Al,J_Al27,S00Al27,S01Al27,S11Al27,    csSIp,csSIn,csSDp,csSDn,dNdE)); 
printf("Si29  %.2E\n",nucleusRecoilAux(Maxwell,  29,Z_Si,J_Si29,S00Si29,S01Si29,S11Si29,    csSIp,csSIn,csSDp,csSDn,dNdE)); 
printf("K39   %.2E\n",nucleusRecoilAux(Maxwell,  39,Z_K,J_K39,S00K39,S01K39,S11K39,         csSIp,csSIn,csSDp,csSDn,dNdE)); 
printf("Ge73  %.2E\n",nucleusRecoilAux(Maxwell,  73,Z_Ge,J_Ge73,S00Ge73,S01Ge73,S11Ge73,    csSIp,csSIn,csSDp,csSDn,dNdE)); 
printf("Nb93  %.2E\n",nucleusRecoilAux(Maxwell,  93,Z_Nb,J_Nb93,S00Nb93,S01Nb93,S11Nb93,    csSIp,csSIn,csSDp,csSDn,dNdE)); 
printf("Te125 %.2E\n",nucleusRecoilAux(Maxwell, 125,Z_Te,J_Te125,S00Te125,S01Te125,S11Te125,csSIp,csSIn,csSDp,csSDn,dNdE));  
printf("I127  %.2E\n",nucleusRecoilAux(Maxwell, 127,Z_I,J_I127,S00I127,S01I127,S11I127,     csSIp,csSIn,csSDp,csSDn,dNdE));
printf("Xe129 %.2E\n",nucleusRecoilAux(Maxwell, 129,Z_Xe,J_Xe129,S00Xe129,S01Xe129,S11Xe129,csSIp,csSIn,csSDp,csSDn,dNdE));  
printf("Xe131 %.2E\n",nucleusRecoilAux(Maxwell, 131,Z_Xe,J_Xe131,S00Xe131,S01Xe131,S11Xe131,csSIp,csSIn,csSDp,csSDn,dNdE));  
printf("Cs133 %.2E\n",nucleusRecoil0Aux(Maxwell,133,Z_Cs,J_Cs133,Sp_Cs133,Sn_Cs133,         csSIp,csSIn,csSDp,csSDn,dNdE));  
}    
