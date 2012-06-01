#include<stdio.h>
#ifndef __4_VERTOR__
#define __4_VERTOR__

extern double vsqrt(double a);
extern double vdot4(int i, int j);
extern void   vsum4(int i, int j, int k, int isg);
extern void   vnull4(int i);
extern void   eps4(int n1, int n2, int n3, int n4);
extern void   pvFill(double mass, double * mom, int pos); 
extern void   lvtonv(char * lv, int nin, int nvpos);
extern double pvect[400];

extern void incomkin(double m1, double m2, double p1, double p2,  
           double *sqrt_S_,  double *Pcm_, double * rapidity_);
           
#endif
