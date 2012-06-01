#include<math.h>
#include<stdlib.h>

#include"../sources/micromegas.h"
#include"../sources/micromegas_aux.h"

static double MQ,MSQ,MNE;

extern double alphaQCD(double);

int main(int argc,char** argv)
{ double EE=0.31223;
  double SW=0.481;
  double CW=sqrt(1-SW*SW);
  double alpha_qcd=0.12;
  double MN=0.939; 
  double G,A,SCcoeff,CS;
  double mq,msq,mne;
  double b1s,b2s,ApB,g,Ap,An,DD,Scale;
  double Mc=1.4,Mb=4.5,Mt=175;  
  int i; 

  if(argc<3) 
  { printf("The program needs 2 arguments: mne,msq\n");
    return 1;
  }
  

  sscanf(argv[1],"%lf",&mne);
  sscanf(argv[2],"%lf",&msq);
Scale=(msq-mne)*(msq+mne)/msq/exp(0.5);
 

  SCcoeff=MN*mne/(MN+mne); SCcoeff*=4/M_PI*3.8937966E8*SCcoeff;
  DD=msq*msq-mne*mne;  
  printf("Tree level approach\n");
  for(i=1;i<=5;i++)
  { int in= i==1? 2:(i==2?1:i);
    ApB=(2-(i&1))*EE/3/CW; ApB*=ApB;
    g=-0.25*ApB/(DD*DD);
    Ap=-1.5*g*mne*MN*parton_x(i,Scale);
    An=-1.5*g*mne*MN*parton_x(in,Scale);

  printf("%d-quark twist-2 cross sections: proton %.2E  neutron %.2E\n",
     i,  Ap*Ap*SCcoeff,An*An*SCcoeff);
  }
  printf("Loop consideration\n");
  
  for(i=4;i<=6;i++)
  {   mq= i==4? 1.3:(i==5? 4.2:171);
    ApB=(2-(i&1))*EE/3/CW; ApB*=ApB;
    b1s =ApB*mne *LintIk(4,msq,mq,mne);	
    b2s =ApB     *LintIk(5,msq,mq,mne)/4.;

    Ap= 1.5*mne*MN*parton_alpha(mq)/(12*M_PI)*(b2s+mne*(b1s)/2)*parton_x(21,mq);
  printf("%d-quark twist-2 amplitude= %E cross sections: proton %.2E \n",
     i, Ap, Ap*Ap*SCcoeff);   
  }
  

printf(" Some tests:\n");
for(i=4;i<=5;i++)
{
  if(i==4) printf("  for c quark\n"); else  printf("  for b quark\n"); 
  mq= (i==4? 1.3:4.2);
  printf("parton_alpha/alphaQCD(%E)=%E\n",mq,parton_alpha(mq)/alphaQCD(mq));

  printf("alpha(mq)/alpha(middle point)=%E\n",parton_alpha(mq)/parton_alpha(sqrt(Scale*mq)));
                            
  printf("check RG solution via middle point  %E  %E\n", parton_x(21,mq)*parton_alpha(sqrt(mq*Scale))/3/M_PI*log(Scale/mq),
                            parton_x(i,Scale));
  printf( "check (I5 +2*mne*I4)/( log approximation)= %E\n", 
   (LintIk(5,msq,mq,mne)+2*mne*LintIk(4,msq,mq,mne))/(4/DD/DD*log(Scale/mq)));
                                                                                    
}  
  
  
   
      
  return 0;
}
