
#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"

static int N=0;
double p=0.8;
double w=1E-5;
double c=0.5;
double A=1;
//double fun(double x) { N++; if(x>0.557)  return sin(x); else return 0;}

double funBW(double x)   { N++; return 1/( x*x +w*w);   } 
double funPole(double x) { N++; return 1/pow(fabs(x),p); }
double funPoleArg(double x, void*p) { N++; double *pp=p; return  pow(fabs(x),*pp);};
double funRand(double x) { N++; return  0.5+ drand48(); }
double funStep(double x) { N++; if(x<0) return 1; else return 0;}
double funExp(double x)  { N++; return exp(A*x);        }
double funSin(double x)  { N++; return sin(x);          }

int main(int argc,char** argv)
{ int err,n,i;

  double a,b,eps, ans, ansSymb;

printf("\n");
N=0; 
  a=-1,b=1,eps=1E-3; 
  printf("  exp(%.2E*x) : [%.2E %.2E] eps=%.2E\n",A,a,b,eps);  
  ans=simpson(funExp,a,b,eps,NULL);
  ansSymb=(exp(A*b)-exp(A*a))/A;
  printf("       err=%.2E N=%d  \n", ansSymb/ans-1,N); 

N=0;  
  a=-10,b=20,eps=1E-3;
  printf("  sin(x)          : [%.2E %.2E]   eps=%.2E\n",a,b,eps); 
  ans=simpson(funSin,a,b,eps,NULL);
  ansSymb=cos(a)-cos(b);     
  printf("        err=%.2E N=%d\n\n", ans/ansSymb-1,N ); 
      
N=0;  
  a=-2, b=1, eps=1E-3;
  printf("  1/(x^2 +(%.2E)^2      : [%.2E %.2E]   eps=%.2E\n", w,a,b,eps);   
  ans=simpson(funBW,a,b,eps,NULL);
  ansSymb=atan(b/w)/w - atan(a/w)/w;
  printf("        err=%.2E N=%d\n\n", ans/ansSymb-1,N );  

N=0;
  a=-1,b=1,eps=1E-3;
  printf("  1/|x|^%.2E     : [%.2E %.2E]   eps=%.2E\n", p,a,b,eps);
  ans=simpson(funPole,a,b,eps,NULL);  
  ansSymb=(pow(fabs(a),1-p) + pow(fabs(b),1-p))/(1-p); 
  printf("          err=%.2E  N=%d\n", ans/ansSymb-1,N); 

N=0; 
  a=-1,b=1,eps=1E-3;
  printf(" Theta(x)      : [%.2E %.2E]   eps=%.2E\n", a,b,eps);
  ans=simpson(funStep,a,b,eps,NULL); 
  ansSymb=b;
  printf("          err=%.2E  N=%d\n", ans/ansSymb-1,N);

N=0;
  a=0,b=1,eps=1E-3;
  printf(" drand48(x) +%.2E       : [%.2E %.2E]   eps=%.2E\n",c, a,b,eps);
  ans=simpson(funRand,a,b,eps,NULL); 
  ansSymb=(b-a)*(c+0.5);
  printf("          err=%.2E  N=%d\n", ans/ansSymb-1,N);
}

