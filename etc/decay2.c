
#include"../sources/micromegas.h"
#include"../sources/micromegas_aux.h"
#include<math.h>

static double m0,m1,m2,w1,w2,z1;

static double Z1(double m){ return  atan((m*m-m1*m1)/(m1*w1));}
static double Z2(double m){ return  atan((m*m-m2*m2)/(m2*w2));}

static double M1(double z){ return  sqrt(fabs(m1*(m1+w1*tan(z))));}
static double M2(double z){ return  sqrt(fabs(m2*(m2+w2*tan(z))));}

static double z2_int(double z)
{ if(m0<=M1(z1)+M2(z)) return 0; else return  decayPcm(m0,M1(z1),M2(z)); }

static double z1_int(double z)
{ z1=z; return simpson(z2_int,Z2(0.),Z2(m0-M1(z1)), 1.E-3);}


int main(int argc, char ** argv)
{ double r;

  if(argc !=6) 
  { 
     printf(" 5 parameters are needes: m0, m1, m2, w1,w2\n");  
     exit(10);
  }

  if(1!=sscanf(argv[1],"%lf",&m0)){ printf("can not read 1 parameter\n"); exit(1);}
  if(1!=sscanf(argv[2],"%lf",&m1)){ printf("can not read 1 parameter\n"); exit(2);};
  if(1!=sscanf(argv[3],"%lf",&m2)){ printf("can not read 1 parameter\n"); exit(3);};
  if(1!=sscanf(argv[4],"%lf",&w1)){ printf("can not read 1 parameter\n"); exit(4);};
  if(1!=sscanf(argv[5],"%lf",&w2)){ printf("can not read 1 parameter\n"); exit(5);};


  if(m0>m1+m2) printf("Pcm=       %E\n",decayPcm(m0,m1,m2));
  r=simpson(z1_int,Z1(0.),Z1(m0),1.E-3)/(M_PI*M_PI);

  printf("integral = %E\n",r);
  exit(0);
}   
