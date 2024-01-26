/*====== Modules ===============
   Keys to switch on 
   various modules of micrOMEGAs  
================================*/

#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"

// based on Eq.5.2,5.4 0902.0688
static long double aa,bb;

static void SommerfeldDeriv( double x, double *phi, double *dphi)
{
   dphi[0]=phi[1];
   if(x==0) dphi[1]=-phi[0];
   else dphi[1]=-2/x*phi[1] - (2*aa*expl(-bb*x)/x+1)*phi[0];
}

double Sommerfeld(double a, double b) 
{  aa=a, bb=b;
   double xstart=0; 
   if(b*xstart > 0.01) xstart=0.01/b;
   if(a*xstart > 0.01) xstart=0.01/a;
   double x=xstart;
   double phi[2];
   double c=(2*a*b-1)/6;
   phi[0]=1-a*x + c*x*x;
   phi[1]= -a+2*c*x;

   double phiArr[300];
   int n=0;   
   double r0=0;
   double stest[4];
   double dx=M_PI/2;
   for(; ;x+=dx)
   { 
     
     odeint(phi,2, x, x+dx, 1E-5, 
             0.1, SommerfeldDeriv);
     for(int i=0;i<3;i++) stest[i]=stest[i+1];
     stest[3]=phi[0]*phi[0]*(x+dx)*(x+dx);              
     double s=0.5*(stest[0]+stest[1]+stest[2]+stest[3]);
     
     if(x>6 && fabs(stest[0]+stest[1]-s)<0.001*s && fabs(stest[1]+stest[2]-s)<0.001*s && fabs(stest[2]+stest[3]-s)<0.001*s
     && fabs(stest[0]-stest[2])< 0.01*(stest[0]+stest[2]) && fabs(stest[1]-stest[3])< 0.01*(stest[1]+stest[3])   )
     { 
//       printf("x=%e stest= %E %E %E %E\n",x,stest[0],stest[1],stest[2],stest[3]);
        return 1/(stest[2]+stest[3]);
     }
   } 

   return 0;
}

