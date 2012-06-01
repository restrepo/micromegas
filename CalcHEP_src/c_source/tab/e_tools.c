#include "e_tools.h"

#include"interface.h"

/*
#define const
#include"num_out.h"
#undef const
*/

/*
int nin_ext, nout_ext;
*/
double sum;
static  char   p_names[MAXNP][20];
static  long   p_codes[MAXNP];
static  double p_masses[MAXNP];


static char* pinf_ext(int nsub,int num , double * mass,long * N) 
{  
  if(nsub>1 || num>nin_int+nout_int) return NULL;
  if(mass) *mass=p_masses[num-1];
  if(N) *N=p_codes[num-1];
  return  p_names[num-1]; 
}


void  boost(double * n, double *p)
{

   double M=p[0]; /* particle mass */
   double shY=n[0]; /* h-sine of rapidity */
   double chY=sqrt(1+ shY*shY);
   double f= ENERGY(M,p+1)*shY + (p[1]*n[1]+p[2]*n[2]+p[3]*n[3])*(chY-1);

   int i;
   for(i=1;i<=3;i++) p[i]+=n[i]*f;
} 
  

void findBoost(double * p, double *n)
{
  double p2=p[1]*p[1]+p[2]*p[2]+p[3]*p[3];
  double M=p[0];  
  int i;

  if(p2 == 0) {for(i=0;i<=3;i++) n[i]=0; return;}  
  p2=sqrt(p2);
  n[0]=p2/M;
  for(i=1;i<=3;i++) n[i]=p[i]/p2;
}

int getNinNout(FILE* flow)
{ if(2!=fscanf(flow,"%d -> %d",&nin_int,&nout_int)) return 1;
  return 0;
}  

int getMasses(FILE * flow)
{
  int i; 
  for(i=0;1==fscanf(flow,"%lf",p_masses+i);i++);  
  if(i!=nin_int+nout_int) return 1; 
  return 0;
}

int getNames(FILE * flow)
{
  int i;
  for(i=0;i<nin_int+nout_int;i++)
  {  
    if(2!=fscanf(flow,"%ld(%[^)]%*c",p_codes+i, p_names+i)) return 1;  
    if(i==nin_int-1) fscanf(flow," -> ");
  }
  pinf_int=&pinf_ext;
  return 0;
}

/*
int decay2(double M, double * p1, double * p2)
{
double m1=p1[0];
double m2=p2[0];
double cosX,sinX, phi;
double P;
int i;

  if (m1+m2 >= M) return 1;

  P= sqrt( (M-m1-m2)*(M+m1+m2)*(M-m1+m2)*(M-m2+m1) )/(2*M);

  cosX=1-2*drand48();
  sinX=sqrt(1-cosX*cosX);

  phi= 2*M_PI*drand48();

  p1[3]=cosX;
  p1[2]=sinX*cos(phi);
  p1[1]=sinX*sin(phi);

  for(i=1;i<=3;i++) { p1[i]=P*p1[i]; p2[i]=-p1[i]; }
  return 0;
}

*/
