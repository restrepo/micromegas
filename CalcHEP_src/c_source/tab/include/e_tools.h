
#ifndef _E_TOOLS_
#define _E_TOOLS_ 

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
 
#define FORMAT " %17.10E"
#define MAXNP 20
#define NAMELEN 20
#define STRSIZE 1000

#define ENERGY(m,p) sqrt((m)*(m) + *(p)**(p) + *(p+1)**(p+1)+ *(p+2)**(p+2))

extern int getNames(FILE * flow);
extern int getMasses(FILE * flow);
extern int getNinNout(FILE * flow);
extern void  boost(double * n, double *p);
extern void findBoost(double * p, double *n); 
extern int decay2(double M, double * p1, double * p2);

#endif
