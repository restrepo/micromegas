#ifndef __DYNAMIC_ME__
#define __DYNAMIC_ME__

#include<stdio.h>

#ifdef __cplusplus
extern "C" {
#endif 

#ifndef  __MICROMEGAS__
#include "../../../include/num_out.h"

typedef struct numout
{
  void * handle;
  double ** link;
  double *Q,*SC,*GG;
  int init;
  CalcHEP_interface * interface; 
} numout;

#endif

#ifndef  __MICROMEGAS__


/*======= Subprocesses ===========*/
  typedef struct txtListStr
  {  struct txtListStr * next;
     char  *txt;
  } txtListStr;

  typedef txtListStr * txtList;

#endif

extern char  * libDir;
extern char  * modelDir;
extern char  * compDir;
extern char  * calchepDir;
extern int   modelNum;

extern int prepareWorkPlace(void);
extern int    cleanWorkPlace(void);
extern numout*getMEcode(int twidth,int UG,char*Process,
                          char * excludeVirtual,char*excludeOut,char*lib);

extern txtList  makeDecayList(char * pname, int nx);
extern txtList  makeProcList(char ** InNames, char** OutNames, int nx);
extern void cleanTxtList(txtList L);
extern void printTxtList(txtList L, FILE *f);

extern double*varAddress(char * name);
/*===========================================================*/
typedef struct{ double width; int dim; txtList pdList[2];}  decayTableStr;

extern double   (*sqme)(int nsub,double *pvect, int * err_code);
extern double   decayPcm(double am0,  double  am1,  double  am2);
extern int      pTabPos(char * name);
extern double   pMass(char * name);
extern int      ForceUG;
extern int      procInfo1(numout*cc, int *nsub, int * nin, int *nout);
extern int      procInfo2(numout*cc, int nsub,char**name,  double*mass);
extern long     pNum(char * name);
extern void     massFilter(double M, txtList * List);
extern void     gammaGluFilter(txtList * List);
extern void     process2Lib(char * process,char * lib);
extern          decayTableStr* decayTable;
extern void     cleanDecayTable(void);
extern void     pname2lib(char*pname, char * libname);
extern double   decay2Info(char * pname, FILE* f);
extern double   pWidth(char *name, txtList * LL,int *dim);
extern double   findBr(txtList L, char * pattern);
extern double   pWidth2(numout * cc, int nsub);
extern char *   pdg2name(int pdg);
extern int  qNumbers(char *pname, int *spin2, int * charge3, int * cdim);

#ifdef __cplusplus
}
#endif 


#endif
