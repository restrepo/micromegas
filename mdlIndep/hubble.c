// 1407.6387  1309.5258 (fig 6)  2303.17021 (fig.3)
#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"

int main(void)
{
  double Year= 31558149.8; // is seconds
  printf(" Age of Universe=%e\n", HubbleTime(0.1,T2_73K)/Year);
  double T_CMB=3000./300./38.681740*1E-9; //GeV
  printf(" Time of CMB =%E\n", HubbleTime(10,T_CMB)/Year);
  printf(" teq=%.2E yr\n", HubbleTime(10,Tmreq )/Year);
  double T1=1E-4, pm1=0.5;
  printf("Free Streaming (p/m=%.2, T1=%.2) = %.2E[Mpc]\n",pm1,T1, freeStreaming(pm1,T1, Tmreq));
  return 0;

}
