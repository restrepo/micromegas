#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../CalcHEP_src/include/extern.h"
#include "../../CalcHEP_src/include/V_and_P.h"
#include "autoprot.h"
extern int  FError;
/*  Special model functions  */

int nModelParticles=23;
ModelPrtclsStr ModelPrtcls[23]=
{
  {"A","A", 22, "0","0",2,1,0}
, {"Z","Z", 23, "MZ","wZ",2,1,0}
, {"G","G", 21, "0","0",2,8,0}
, {"W+","W-", 24, "MW","wW",2,1,3}
, {"n1","N1", 12, "0","0",1,1,0}
, {"e1","E1", 11, "0","0",1,1,-3}
, {"n2","N2", 14, "0","0",1,1,0}
, {"e2","E2", 13, "Mm","0",1,1,-3}
, {"n3","N3", 16, "0","0",1,1,0}
, {"e3","E3", 15, "Mt","0",1,1,-3}
, {"u","U", 2, "Mup","0",1,3,2}
, {"d","D", 1, "Mdo","0",1,3,-1}
, {"c","C", 4, "Mc","0",1,3,2}
, {"s","S", 3, "Ms","0",1,3,-1}
, {"t","T", 6, "Mtop","wtop",1,3,2}
, {"b","B", 5, "Mb","0",1,3,-1}
, {"H","H", 37, "MH","wH",0,1,0}
, {"~H+","~H-", 37, "MH1","whc",0,1,3}
, {"~A0","~A0", 36, "MA0","wA0",0,1,0}
, {"~H0","~H0", 35, "MH0","wH0",0,1,0}
, {"~n1","~n1", 9900012, "MN1","wn1",1,1,0}
, {"~n2","~n2", 9900014, "MN2","wn2",1,1,0}
, {"~n3","~n3", 9900016, "MN3","wn3",1,1,0}
};
int nModelVars=51;
int nModelFunc=15;
int LastVar=51;
char*varNames[66]={
 "EE","SW","s12","s23","s13","MZ","wZ","wW","Mm","Mt"
,"Mup","Mdo","Mc","Ms","Mtop","wtop","Mb","MH","wH","MH1"
,"whc","MA0","wA0","MH0","wH0","La2","LaL","h11R","h12R","h13R"
,"h21R","h22R","h23R","h31R","h32R","h33R","h11I","h12I","h13I","h21I"
,"h22I","h23I","h31I","h32I","h33I","MN1","wn1","MN2","wn2","MN3"
,"wn3","Sqrt2","CW","c12","c23","c13","Vud","Vus","Vub","Vcd"
,"Vcs","Vcb","Vtd","Vts","Vtb","MW"};
double varValues[66]={
   3.133300E-01,  4.740000E-01,  2.210000E-01,  4.000000E-02,  3.500000E-03,  9.118700E+01,  2.502000E+00,  2.094000E+00,  1.057000E-01,  1.777000E+00
,  3.000000E-03,  6.000000E-03,  8.290000E-01,  2.000000E-01,  1.700000E+02,  1.442000E+00,  3.010000E+00,  1.300000E+02,  1.461000E+00,  2.800000E+02
,  1.000000E+00,  2.800000E+02,  1.000000E+00,  2.800000E+02,  1.000000E+00,  5.000000E-01, -4.000000E-01,  1.000000E+00,  1.000000E+00,  1.000000E+00
,  1.000000E+00,  1.000000E+00,  1.000000E+00,  1.000000E+00,  1.000000E+00,  1.000000E+00,  1.000000E+00,  1.000000E+00,  1.000000E+00,  1.000000E+00
,  1.000000E+00,  1.000000E+00,  1.000000E+00,  1.000000E+00,  1.000000E+00,  1.100000E+02,  0.000000E+00,  2.000000E+02,  0.000000E+00,  3.000000E+02
,  0.000000E+00};
int calcMainFunc(void)
{
   int i;
   static double * VV=NULL;
   static int iQ=-1;
   static int cErr=0;
   double *V=varValues;
   FError=0;
   if(VV&&!cErr)
   { for(i=0;i<nModelVars;i++) if(i!=iQ && VV[i]!=V[i]) break;
     if(i==nModelVars)     return 0;
   }
  cErr=1;
   V[51]=sqrt(2);
 LastVar=51;    if(!finite(V[51]) || FError) return 51;
   V[52]=sqrt(1-pow(V[1],2));
 LastVar=52;    if(!finite(V[52]) || FError) return 52;
   V[53]=sqrt(1-pow(V[2],2));
 LastVar=53;    if(!finite(V[53]) || FError) return 53;
   V[54]=sqrt(1-pow(V[3],2));
 LastVar=54;    if(!finite(V[54]) || FError) return 54;
   V[55]=sqrt(1-pow(V[4],2));
 LastVar=55;    if(!finite(V[55]) || FError) return 55;
   V[56]=V[53]*V[55];
 LastVar=56; 
   V[57]=V[2]*V[55];
 LastVar=57; 
   V[58]=V[4];
 LastVar=58; 
   V[59]=-V[53]*V[3]*V[4]-V[2]*V[54];
 LastVar=59; 
   V[60]=V[53]*V[54]-V[2]*V[3]*V[4];
 LastVar=60; 
   V[61]=V[3]*V[55];
 LastVar=61; 
   V[62]=V[2]*V[3]-V[53]*V[54]*V[4];
 LastVar=62; 
   V[63]=-V[2]*V[54]*V[4]-V[53]*V[3];
 LastVar=63; 
   V[64]=V[54]*V[55];
 LastVar=64; 
   V[65]=V[5]*V[52];
 LastVar=65; 
   if(VV==NULL) 
   {  VV=malloc(sizeof(double)*nModelVars);
      for(i=0;i<nModelVars;i++) if(strcmp(varNames[i],"Q")==0) iQ=i;
   }
   for(i=0;i<nModelVars;i++) VV[i]=V[i];
   cErr=0;
return 0;
}
