#include"SLHAplus.h"
#include"aList.h"

/* transforms routines with variable  number  of  arguments
   into set of routines with fixed number of arguments. It is needed for 
   Fortran  and C++ ( in case of Complex arguments) versions of SLHAplus.
   fixArg.c uses substitutions  defined in aList.h 
*/                  

int rDiagonal2(aList3(REAL))   {return rDiagonal(2,aList3());}
int rDiagonal3(aList6(REAL))   {return rDiagonal(3,aList6());}
int rDiagonal4(aList10(REAL))  {return rDiagonal(4,aList10());}
int rDiagonal5(aList15(REAL))  {return rDiagonal(5,aList15());}

int rDiagonalA2(aList4(REAL))  {return rDiagonalA(2,aList4());}
int rDiagonalA3(aList9(REAL))  {return rDiagonalA(3,aList9());}
int rDiagonalA4(aList16(REAL)) {return rDiagonalA(4,aList16());}
int rDiagonalA5(aList25(REAL)) {return rDiagonalA(5,aList25());}

int cDiagonalH2(aList3(COMPLEX)) {return cDiagonalH(2,aList3());}
int cDiagonalH3(aList6(COMPLEX)) {return cDiagonalH(3,aList6());}
int cDiagonalH4(aList10(COMPLEX)){return cDiagonalH(4,aList10());}
int cDiagonalH5(aList15(COMPLEX)){return cDiagonalH(5,aList15());}

int cDiagonalS2(aList3(COMPLEX)) {return cDiagonalS(2,aList3());}
int cDiagonalS3(aList6(COMPLEX)) {return cDiagonalS(3,aList6());}
int cDiagonalS4(aList10(COMPLEX)){return cDiagonalS(4,aList10());}
int cDiagonalS5(aList15(COMPLEX)){return cDiagonalS(5,aList15());}

int cDiagonalA2(aList4(COMPLEX)) {return cDiagonalA(2,aList4());}
int cDiagonalA3(aList9(COMPLEX)) {return cDiagonalA(3,aList9());}
int cDiagonalA4(aList16(COMPLEX)){return cDiagonalA(4,aList16());}
int cDiagonalA5(aList25(COMPLEX)){return cDiagonalA(5,aList25());}

int System1(char*format){ return System(format);}
int System2(char*format,char*path){ return System(format,path);}

int aPrintF0(char * format) { return aPrintF(format);}
int aPrintF1(char*format,double x1) { return aPrintF(format,x1);}
int aPrintF2(char*format,double x1,double x2){return aPrintF(format,x1,x2);}
int aPrintF3(char * format, double x1,double x2,double x3)
    { return aPrintF(format,x1,x2,x3);}
int aPrintF4(char * format, double x1,double x2,double x3,double x4)
    { return aPrintF(format,x1,x2,x3,x4);}
int aPrintF5(char * format, double x1,double x2,double x3,double x4,double x5)
    { return aPrintF(format,x1,x2,x3,x4,x5);}
int aPrintF6(char * format, double x1,double x2,double x3,double x4,double x5,double x6)
    { return aPrintF(format,x1,x2,x3,x4,x5,x6);}
int aPrintF7(char * format, double x1,double x2,double x3,double x4,double x5,double x6,double x7)
    { return aPrintF(format,x1,x2,x3,x4,x5,x6,x7);}
int aPrintF8(char * format, double x1,double x2,double x3,double x4,double x5,double x6,double x7,double x8)
    { return aPrintF(format,x1,x2,x3,x4,x5,x6,x7,x8);}
int aPrintF9(char * format, double x1,double x2,double x3,double x4,double x5,double x6,double x7,double x8,double x9)
    { return aPrintF(format,x1,x2,x3,x4,x5,x6,x7,x8,x9);}

double dbl(REAL a){ return  (double)(a);}
double slhaVal0(char *block, double scale){                             return   slhaVal(block,scale,0);}
double slhaVal1(char *block, double scale,int i1){                      return   slhaVal(block,scale,1,i1);}
double slhaVal2(char *block, double scale,int i1,int i2){               return   slhaVal(block,scale,2,i1,i2);}
double slhaVal3(char *block, double scale,int i1,int i2,int i3){        return   slhaVal(block,scale,3,i1,i2,i3);}
double slhaVal4(char *block, double scale,int i1,int i2,int i3,int i4){ return   slhaVal(block,scale,4,i1,i2,i3,i4);}

int    slhaValExists0(char *block){                                     return   slhaValExists(block,0) ;}
int    slhaValExists1(char *block,int i1){                              return   slhaValExists(block,1,i1);}
int    slhaValExists2(char *block,int i1,int i2){                       return   slhaValExists(block,2,i1,i2);}
int    slhaValExists3(char *block,int i1,int i2,int i3){                return   slhaValExists(block,3,i1,i2,i3);}
int    slhaValExists4(char *block,int i1,int i2,int i3,int i4){         return   slhaValExists(block,4,i1,i2,i3,i4);}


    
#include"delList.h"
