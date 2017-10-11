#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<unistd.h>

/*=========================================
  Experimental constrains on the UMSSM in C
  =========================================*/
extern double deltarho(void);
extern int    masslimits(void);
extern int    Zinv(void);
extern int    HBblocks(char * fname);
extern int    LiLithF(char*fname);

/*============
  Zprime limit
  ============*/
extern double sBthtE6Fixed(double xx, int jj);
extern double tE6vssBorMZ2lim(double xx, int jj);
extern double sBthForAnytE6MZ2(double xx, double tE6);
extern int Zprimelimits(void);

/*==========
  UMSSMTools
  ==========*/
extern int    assignValFunc(char * name, double val);
extern int    UMSSMTools(void);
extern int    read_prmU(int n);
extern double UparC(int r);
extern double uparctof_(int *r);

/*=====
  Other
  =====*/
extern void o1Contents(FILE * f);
extern double randpar(char * xx, double start, double amp, int pwr, int sign, int logmod, int expval);
extern double SignParam (double x);
