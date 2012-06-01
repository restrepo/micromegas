#ifndef __DECAY__
#define __DECAY__

extern  void decay_0(int nvin, double amm1, double amm2, double *factor, double *Emax);
extern  void decay_1(int nvpole, double *hsum, double *hdif);
extern  void decay_2(int nvout1, double *parcos);
extern  void decay_3(int nvath, double parcos, double parfi, int nvout1,int nvout2);
extern  void wrmom(int n);

#endif
