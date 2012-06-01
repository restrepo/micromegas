
#include <stdio.h>
#include <stdio.h>
#include <math.h>   

#ifndef noCMPLX
#include <complex.h>
#define Complex double complex
#endif 
        
extern int slhaRead(char *fname,int mode);
extern double slhaVal(char * Block, double Q, int nKey, ...);
extern int slhaValExists(char * Block, int nKey, ...);
extern double slhaWidth(int pNum);
extern int slhaWrite(char *fname);
extern int slhaWarnings(FILE*f);
extern int slhaDecayExists(int pNum);
extern double slhaBranch(int pNum,int N, int * nCh);
extern int initDiagonal(void);
extern int rDiagonal(int nDim,...); 
extern int rDiagonalA(int nDim,...);
extern double MassArray(int id,  int i);
extern double MixMatrix(int id, int i,int j); 
extern double MixMatrixU(int id, int i,int j);
extern int cDiagonalH(int Dim,...);
extern int cDiagonalA(int Dim,...);
extern int cDiagonalS(int Dim,...);
#ifndef noCMPLX
extern Complex cMixMatrix(int id,int i,int j);
extern Complex cMixMatrixU(int id,int i,int j);
#endif
extern int System(char * format, ...);
extern int openAppend(char * fileName);
extern int aPrintF(char * format,...);

extern double initQCD(double,double,double,double);
extern double MbEff(double);
extern double MtEff(double);
extern double McEff(double);
extern double alphaQCD(double);

/* To avoid avto-prototyping  

extern double  sqrt(double);
extern double  sin(double);
extern double  cos(double);
extern double  tan(double);
extern double  asin(double);
extern double  acos(double);
extern double  atan(double);
extern double  exp(double);
extern double  log(double); 
extern double  pow(double,double);
extern double  fabs(double);
extern double  atan2(double,double);
extern double  log10(double);
extern double  sinh(double);
extern double  cosh(double);
extern double  tanh(double);
extern double  asinh(double);
extern double  acosh(double);
extern double  atanh(double);

extern double  creal(Complex);
extern double  cimag(Complex);
extern double  carg(Complex);
extern double  cabs(Complex);
extern Complex conj(Complex);
extern Complex cacos(Complex);
extern Complex casin(Complex);
extern Complex catan(Complex);
extern Complex ccos(Complex);
extern Complex csin(Complex);
extern Complex ctan(Complex);
extern Complex cacosh(Complex);
extern Complex casinh(Complex);
extern Complex catanh(Complex);
extern Complex ccosh(Complex);
extern Complex csinh(Complex);
extern Complex ctanh(Complex);
extern Complex cexp(Complex);
extern Complex clog(Complex);
extern Complex clog10(Complex);
extern Complex cpow(Complex,Complex);
extern Complex csqrt(Complex);
extern Complex cproj(Complex);

extern int printf(char*, ...); 

extern double  if(double,double,double);


*/
