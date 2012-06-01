#ifndef __NUM_IN_
#define __NUM_IN_

#include<stdlib.h>
#include<string.h> 
#include<math.h>

extern double alpha_2(double);
typedef double (DNN)(double *,int *);
typedef double (FNN)(void);

/*extern  double sqrMom(int, char*, double*); */
extern  double Fmax;
extern  double DP[]; 
extern  double Helicity[2];
extern  double HelicityN[2];
extern  double *Q0,*Q1,*Q2;
extern  int    CalcConst;
extern  int    indx_(int k,int l);
extern  void sprod_(int ntot, double * momenta);
extern  int    prepDen(int nden, int nin, double BWrange2,   
                       double * dmass,double * dwidth, char **q,double * mom);
                       
#endif
