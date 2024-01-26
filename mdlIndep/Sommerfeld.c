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

double Sommerfeld(double a, double b)    // 0902.0688 Eq.5.2-5.4
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



double mu=100, Mr=1000, alpha=0.1, vv=0.01;
// mu - mass of mediator; Mr - reduced mass ~Mdm/2; alpha - Yukawa couplingh; vv-relative velosity

typedef struct { double mu;    //  mass of mediator 
                 double Mr;    //  reduced mass 
                 double alpha; //  Yukawa coupling
                 double vr;    //  relative velosity
               } SommPar;

double Sv(double v, void * par_)
{
  SommPar *par=par_;
  double a=par->alpha/v, b=par->mu/par->Mr/v;
  return Sommerfeld(a, b);
}

double Sfrac(double f, void * par_)
{
  SommPar *par=par_;
  double a=par->alpha/f/par->vr, b=par->mu/par->Mr/par->vr;
  return Sommerfeld(a, b)*f*f;
}

double Smu(double mu, void * par_)
{
  SommPar *par=par_;
  double a=par->alpha/par->vr, b=mu/par->Mr/par->vr;
  return Sommerfeld(a, b);
}

double SMcdm(double M, void * par_)
{
  SommPar *par=par_;
  double a=par->alpha/par->vr, b=par->mu/(M/2)/par->vr;
  return Sommerfeld(a, b);
}


    
double Scoul(double v,void*par_) 
{  SommPar *par=par_;
   double a=2*par->alpha/v;  
   return M_PI*a/(1-exp(-M_PI*a));
}


double SHulthen(void*par_)
{
   SommPar *par=par_;
   double eps_v=par->vr/par->alpha, eps_f=0.5*par->mu/par->alpha/par->Mr;
   
   return M_PI/eps_v*sinh(12*eps_v/(M_PI*eps_f))/
   ( coshl(12*eps_v/(M_PI*eps_f)) - cosl(2*sqrtl(6/eps_f)- powl(6*eps_v/M_PI/eps_f,2)) ); 
}

double Shulthen_mu(double mu, void*par_)
{   SommPar *par=par_;
    par->mu=mu;
    return SHulthen(par_);
}

double Shulthen_v(double v, void*par_)
{   SommPar *par=par_;
    par->vr=v;
    return SHulthen(par_);
}


int main(int argc,char** argv)
{  int err,n,i;

//   SommPar std={100, 1000, 0.1, 1E-2};
   SommPar std={1000,10000, 1, 3E-1};
   SommPar Coul=std, Pole=std, antiPole=std;


std.alpha=0.001;
std.mu=120;
std.Mr=10000;
std.vr=1E-3;

displayPlot("Sfrac^2","frac", 1E-4,1, 1,1,"S(alpha/frac)*frac^2",0,Sfrac,&std);

//displayPlot("S(v)","v", 1E-4,1, 1,1,"Sv",0,Sv,&std);

//displayPlot("S(v)","M", 1000,100000, 1,1,"SMcdm",0,SMcdm,&std);

exit(0); 
   char mess[100]; 
   sprintf(mess,"Expected resonance mu=%.2E",6*alpha*(2*Mr)/M_PI/M_PI);
   displayPlot(mess,"mu", 1,100000, 1,2,"Smu",0,Smu,&Pole, "hulthen", 0, Shulthen_mu,&Pole);


//displayPlot(mess,"v", 1E-3,1, 0,2,"Smu",0,Sv,&Pole, "hulthen", 0, Shulthen_v,&Pole);

//displayPlot("Comparison","v", 1E-4,100, 1,1 , "Coulomb",          0, Sv,&std);

exit(0);
   Coul.mu=0; // Coulomb
   Pole.mu=120;
   antiPole.mu=60;
   
   displayPlot("Comparison","v", 1E-4,1, 1,4 , "Coulomb",          0, Scoul,&Coul
                                             , "Sommerfeld mu=0",  0, Sv,   &Coul
                                             , "Sommerfeld mu=120",0, Sv,   &Pole
                                             , "Sommerfeld mu=60", 0, Sv,   &antiPole);
  

/*
   mu=6*alpha*(2*Mr)/M_PI/M_PI; // 1005.4678 Eq.7 
   sprintf(mess,"mu=%.1E ( resonane)",mu); 
   displayPlot(mess,"v", 1E-4,1, 1,2,"Sv",0,Sv,NULL,"Coulomb",0,Scoul,NULL);

   mu=6*alpha*(2*Mr)/M_PI/M_PI/2;  // anti-resonance
   sprintf(mess,"mu=%.1E (anti-resonane)",mu);  
   displayPlot(mess,"v", 1E-4,1, 1,2,"Sv",0,Sv,NULL,"Coulomb",0,Scoul,NULL);
*/            
   
   return 0;
}  
