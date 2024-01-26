 // vertex:  C  c  h 
/*******************************
*    CalcHEP  3.8.9*
*******************************/
#include<stdlib.h>
#include<math.h>
#include<complex.h>
#include"nType.h"
static char*varName[5]={"zero"
,"EE"
,"SW"
,"MW"
,"Mc"};
 static double V[5];
static int vertexCoeff(double * coeff_out)
{
 double R=(double)(-V[4]*V[1])/(double)(+2*V[3]*V[2]);
 coeff_out[0]=R*(+1);
 return 0;
}
 static char *SymbVert[1]={"1"};
typedef struct  lVert
{
    void * handle;
    void ** link;
    int init;
    int GGpower; int nVar; char **varNames; double *varValues; int nTerms; char **SymbVert;  int (*calcCoeff)(double*);
}  lVert;
 extern lVert vert_a13p13p17;
 lVert vert_a13p13p17={NULL,NULL,0,0,4, varName+1, V, 1, SymbVert, vertexCoeff};
