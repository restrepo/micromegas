#include"../sources/micromegas.h"

static double r0=1;
double rhoGauss(double r)
{ 
  if(r> 5*r0) return 0;
  return  Mcdm*Mcdm/pow(M_PI,3./2.)*exp(-r*r/(r0*r0))/r0/r0/r0;
}
 
double fluxGauss(double E,double dSigmadE)
{ double mp=0.938;
  double p=sqrt((E+mp)*(E+mp)-mp*mp);
  double Kdif=K_dif*pow(p,Delta_dif);
  double C=dSigmadE*0.9460729E24/(32*M_PI*M_PI*Kdif);
  int n;
  double sum=1/Rsun;
  double dsum=sum;

 for(n=1;  dsum>sum*1.E-3 ;n=n+2)
  { dsum=2*(1/sqrt(Rsun*Rsun +4*L_dif*L_dif*n*n)-1/sqrt(Rsun*Rsun +4*L_dif*L_dif*(n+1)*(n+1)));
    sum-=dsum;
  }
    sum-= 1/(2*L_dif*n);
  return C*sum;
}


int main(int argc,char** argv)
{
 double Etest=20, dSigmadE=1.E-30;
 
 Vc_dif=0;
 Gtot_style=0;
 rhoDM=1;
 setHaloProfiles(noClamps, rhoGauss);
 printf("Test:pbarFlux= %E  analitic= %E\n", 
       pbarFlux(Etest,dSigmadE), fluxGauss(Etest,dSigmadE));
 return 0;       
}
