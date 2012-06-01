/*
Copyright (C) 1997, Slava Ilyin
*/

#include<math.h>
#include"4_vector.h"

double pvect[400];  /*was [4][100] */

void lvtonv(char *lv, int nin, int nv)
{
  int i,n;
  vnull4(nv);
  for(i=0;n=lv[i] ;i++) if(n>nin) vsum4(nv,n,nv,1); else vsum4(nv,n,nv,-1);
} 

/* ****************************************** */
/*    Scalar product of two 4-vectors:     * */
/* ****************************************** */
double vdot4(int i, int j)
{
   i=4*i-4;
   j=4*j-4;
   return   pvect[i]*pvect[j] - pvect[i+1]*pvect[j+1] 
          - pvect[i+2]*pvect[j+2] - pvect[i+3]*pvect[j+3];
} 

/* ******************************************************* */
/*       SUM or Difference of two 4-vectors:            *  */
/* ISG=1  sum   and  ISG=-1   difference                *  */
/*             P(I) + ISG*P(J)=>P(K)                    *  */
/* ******************************************************* */

void vsum4(int i, int j, int k, int isg)
{   int l;
    i = 4*i - 4;
    j = 4*j - 4;
    k = 4*k - 4;
    if (isg == 1) {for(l=0;l<4;l++) pvect[k+l] = pvect[i+l] + pvect[j+l];}
    else          {for(l=0;l<4;l++) pvect[k+l] = pvect[i+l] - pvect[j+l];}
} /* vsum4_ */

/* ****************************************** */
/* NULLification of a 4-vector:  P(I) => 0 * */
/* ****************************************** */

void vnull4(int  i)
{
  int k;      
  for(k=4*i-4; k<4*i; k++)  pvect[k]=0;
} 

/* rest of this file is written by A.Pukhov */

void eps4(int n1, int  n2, int n3, int n4)
{
    double a10, a11, a12, a13, a20, a21, a22, a23, a30, a31, a32, 
	    a33, d1021, d1022, d1023, d1122, d1123, d1223;

n1*=4; n2*=4; n3*=4; n4*=4;

    a10 = pvect[n1 - 4];
    a11 = pvect[n1 - 3];
    a12 = pvect[n1 - 2];
    a13 = pvect[n1 - 1];
    a20 = pvect[n2 - 4];
    a21 = pvect[n2 - 3];
    a22 = pvect[n2 - 2];
    a23 = pvect[n2 - 1];
    a30 = pvect[n3 - 4];
    a31 = pvect[n3 - 3];
    a32 = pvect[n3 - 2];
    a33 = pvect[n3 - 1];
/*                               A10  A20  A30  X0 */
/*                               A11  A21  A31  X1 */
/*                               A12  A22  A32  X2 */
/*                               A13  A23  A33  X3 */
    d1021 = a10 * a21 - a20 * a11;
    d1022 = a10 * a22 - a20 * a12;
    d1023 = a10 * a23 - a20 * a13;
    d1122 = a11 * a22 - a21 * a12;
    d1123 = a11 * a23 - a21 * a13;
    d1223 = a12 * a23 - a22 * a13;
    pvect[n4 - 4] =   a31 * d1223 - a32 * d1123 + a33 * d1122;
    pvect[n4 - 3] =   a30 * d1223 - a32 * d1023 + a33 * d1022;
    pvect[n4 - 2] = -(a30 * d1123 - a31 * d1023 + a33 * d1021);
    pvect[n4 - 1] =   a30 * d1122 - a31 * d1022 + a32 * d1021;
} /* eps4_ */

void pvFill(double mass, double * mom4, int pos)
{

  int i,i0=4*(pos-1);
  pvect[i0]=mass;
  pvect[i0]*=pvect[i0];
  mass=mass*mass;

  for(i=1;i<4;i++) {pvect[i0+i]=mom4[i]; pvect[i0]+=pvect[i0+i]*pvect[i0+i];}
  pvect[i0]=sqrt(pvect[i0]);  
}



void incomkin(double m1, double m2, double p1, double p2, 
           double *sqrt_S_,  double *Pcm_, double * rapidity_)
{
  double sqrt_S, Pcm,rapidity;
  double e1=sqrt(m1*m1+p1*p1);
  double e2=sqrt(m2*m2+p2*p2);
   
  sqrt_S=(e1+e2)*(e1+e2)-(p1-p2)*(p1-p2);
  
  rapidity= atanh((p1-p2)/(e1+e2));

  Pcm=p1*cosh(rapidity)-e1*sinh(rapidity);

  if(sqrt_S_) *sqrt_S_=sqrt(sqrt_S);
  if(Pcm_) *Pcm_=Pcm; 
  if(rapidity_) *rapidity_=rapidity;
}


