#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//#include"../SLHAplus/SLHAplus.h"
#include"../../CalcHEP_src/c_source/SLHAplus/include/SLHAplus.h"
#include<complex.h>
#define Complex double complex


// Initialisation of functions :

double sgn (double x);
double min (double a, double b);
double max (double a, double b);

// Functions needed for cubic interpolation :

int  leftXN(int n,int dim,  double * xa, double x);
double  polintN(double x, int n,  double *xa, double *ya);
double polint(double x, int n,  double *xa, double *ya);

// Cubic interpolation for alfSQS :

double polint4alfSQS(double QSUSY);


// Definition of functions :

double sgn (double x)
{
  if (x>=0) return 1;
  else return -1;
}

double min (double a, double b)
{
  double mini;
  if (a>b) mini=b;
  else mini=a;
  return mini;
}


double max (double a, double b)
{
  double maxi;
  if (a>b) maxi=a;
  else maxi=b;
  return maxi;
}

int  leftXN(int n,int dim,  double * xa, double x)
{  int k1,k2,k3;

   k1=n/2;                         
   k2=dim-(n+1)/2-1;

   if(xa[0]< xa[dim-1])
   {              
     if(x<=xa[k1]) return 0;
     if(x>=xa[k2]) return dim-n;                   
     while(k2-k1>1)                
     { k3=(k1+k2)/2;               
       if(xa[k3]>x)k2=k3; else k1=k3;
     }
   } else 
   {  
     if(x>=xa[k1]) return 0;
     if(x<=xa[k2]) return dim-n;
     while(k2-k1>1)                
     { k3=(k1+k2)/2;               
       if(xa[k3]<x)k2=k3; else k1=k3;
     }
   }
   return k1+1-n/2;
}

double  polintN(double x, int n,  double *xa, double *ya)
{  double z[20];
   int i,m;
   for(i=0;i<n;i++) z[i]=ya[i];
   for(m=1;m<n;m++) for(i=0;i<n-m;i++)
   z[i]=(z[i]*(xa[i+m]-x) - z[i+1]*(xa[i]-x))/(xa[i+m]-xa[i]);
   return z[0];
}

double polint(double x, int n,  double *xa, double *ya)
{ int shift=leftXN(4,n, xa, x);
   return polintN(x,4,xa+shift, ya+shift);
}

double polint4alfSQS(double QSUSY)
{
	double scale[100]={100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,2600,2700,2800,2900,3000,3100,3200,3300,3400,3500,3600,3700,3800,3900,4000,4100,4200,4300,4400,4500,4600,4700,4800,4900,5000,5100,5200,5300,5400,5500,5600,5700,5800,5900,6000,6100,6200,6300,6400,6500,6600,6700,6800,6900,7000,7100,7200,7300,7400,7500,7600,7700,7800,7900,8000,8100,8200,8300,8400,8500,8600,8700,8800,8900,9000,9100,9200,9300,9400,9500,9600,9700,9800,9900,10000}; // scale corresponding to sqrt(QSTSB)
	double alfsca[100]={0.11677,0.10604,0.10105,0.09779,0.09541,0.09354,0.09203,0.09075,0.08965,0.08870,0.08785,0.08709,0.08640,0.08577,0.08520,0.08467,0.08417,0.08371,0.08328,0.08288,0.08250,0.08214,0.08180,0.08148,0.08117,0.08087,0.08059,0.08033,0.08007,0.07982,0.07958,0.07936,0.07913,0.07892,0.07872,0.07852,0.07833,0.07814,0.07796,0.07779,0.07762,0.07745,0.07729,0.07713,0.07698,0.07683,0.07669,0.07655,0.07641,0.07628,0.07615,0.07602,0.07589,0.07577,0.07565,0.07553,0.07542,0.07531,0.07520,0.07509,0.07498,0.07488,0.07478,0.07468,0.07458,0.07448,0.07439,0.07429,0.07420,0.07411,0.07402,0.07394,0.07385,0.07377,0.07368,0.07360,0.07352,0.07344,0.07337,0.07329,0.07321,0.07314,0.07306,0.07299,0.07292,0.07285,0.07278,0.07271,0.07264,0.07258,0.07251,0.07244,0.07238,0.07232,0.07225,0.07219,0.07213,0.07207,0.07201,0.07195};
	return polint(QSUSY,100,scale,alfsca); // corresponding alphas(sqrt(QSTSB))
}
