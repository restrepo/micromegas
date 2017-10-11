#ifndef __POLINT__
#define __POLINT__

extern double polint3(double x, int dim,  double *xa, double *ya);
extern double polint1(double x, int n,  double *xa, double *ya);
extern double polint1Exp(double x, int n,  double *xa, double *ya);
extern double  polintN(double x, int n,  double *xa, double *ya);
#endif
