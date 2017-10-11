#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<unistd.h>

/*======================================
  Higgs, LEP limits and deltarho
  ======================================*/
extern double deltarho_(void);
extern int    masslimits_(void);
extern int    hbBlocksMDL(char*fname,int*nHch);
extern int    LilithMDL(char*fname);

#define deltarho   deltarho_
#define masslimits masslimits_

/*==========
  UMSSMTools
  ==========*/
extern int    assignValFunc(char * name, double val);
extern int    umssmtools(int PDG_LSP);
extern int    read_prmU(int n);
extern double UparC(int r);
extern double uparctof_(int *r);

extern double bsg_(double *M, double*P);
extern double deltamd_(double *M, double*P);
extern double deltams_(double *M, double*P);
extern double bsmumu_(double *M, double*P);
extern double btaunu_(double *M, double*P);
extern double gmuon_(double *M, double*P);
extern double bxislllow_(double *M, double*P);
extern double bxisllhigh_(double *M, double*P);

extern double bdg_(double *M, double *P);
extern double bdmumu_(double *M, double *P);
extern double bxisnunu_(double *M, double *P);
extern double bpkpnunu_(double *M, double *P);
extern double bksnunu_(double *M, double *P);
extern double rdtaul_(double *M, double *P);
extern double rdstaul_(double *M, double *P);

extern double kppipnunu_(double *M, double *P);
extern double klpi0nunu_(double *M, double *P);
extern double deltamk_(double *M, double *P);
extern double epsk_(double *M, double *P);

#define bsg        bsg_
#define deltamd    deltamd_
#define deltams    deltams_
#define bsmumu     bsmumu_
#define btaunu     btaunu_
#define gmuon      gmuon_
#define bxislllow  bxislllow_
#define bxisllhigh bxisllhigh_

#define bdg        bdg_
#define bdmumu     bdmumu_
#define bxisnunu   bxisnunu_
#define bpkpnunu   bpkpnunu_
#define bksnunu    bksnunu_
#define rdtaul     rdtaul_
#define rdstaul    rdstaul_

#define kppipnunu  kppipnunu_
#define klpi0nunu  klpi0nunu_
#define deltamk    deltamk_
#define epsk       epsk_

/*=====================
  Les Houches interface
  =====================*/
extern int  readLesH(char*fname, int SM );
extern int  lesHinput(char * fname);
extern void FillVal(int mode);

/*=====
  Other
  =====*/
extern void   o1Contents(FILE * f);
extern double randpar(char * xx, double start, double amp, int pwr, int sign, int logmod, int expval);
extern double SignParam (double x);
