#include<math.h>
#include<stdio.h>

typedef struct CalcHEP_interface
{
/* variables */
  int nvar; int nfunc; char ** varName; double* va;
/* precesses */
  int nin; int nout; int nprc; char* (*pinf)(int, int , double*,long *);
/* constraints calculation */
  int (*calcFunc)(void);
/* width service */
  double*BWrange;  int*twidth;  int*gtwidth; int*gswidth;
/* matrix element calculation */
  double (*sqme)(int , double*, int*);
/* information about propagators */
  char * (*den_info)(int, int, int *, int*);
/* color chains  service */
  void (*build_cb)(int);  int *cb_pow;  int *cb_nc;int**cb_chains; 
  double**cb_coeff;void (*destroy_cb)(void);  
} CalcHEP_interface;

 
extern CalcHEP_interface interface_out; /* generated code */
CalcHEP_interface * code;                 /* internal code pointer */ 

static int assignVal(char * name, double val)
{
  int i;
  for(i=1;i<=code->nvar+code->nfunc;i++)
  { 
    if(strcmp(name,code->varName[i])) continue;
    if(i<=code->nvar) {code->va[i]=val; return 0;} else return 1; 
  }
    return 2;
}
                      
                      

int main(void)
{ int err;

  double pvec[16];
  double *p1=pvec, *p2=pvec+4, *p3=pvec+8, *p4=pvec+12;

  double m1,m2,m3,m4;

  double Pin,Pout;
  double totcoef, sqrt_S, S, lambda12, lambda34,ms,md;
  double cos_fi, cos_step;
  double sigmaTot=0;
  int i;

  code=&interface_out;
/*an example how to change parametgers; if name of parameter is wrong, err=1 */ 
/*  err=assignVal("SW", 0.483);  */

  err=code->calcFunc();
  if(err>0) 
  { printf("Can not calculate constraints '%s'\n",code->varName[err]);
    exit(1);
  }      
/* find masses */ 
  code->pinf(1,1,&m1,NULL);code->pinf(1,2,&m2,NULL);  
  code->pinf(1,3,&m3,NULL);code->pinf(1,4,&m4,NULL);

/* Initial cms  momentum */

  Pin=100;

/* 2-2 kinematics */

  sqrt_S=sqrt(m1*m1+Pin*Pin) +sqrt(m2*m2+Pin*Pin);
  S=sqrt_S*sqrt_S;  
  lambda12=2*sqrt_S*Pin;
  ms = m3+m4; if (ms >= sqrt_S) return 1;
  md = m3 -m4;
  lambda34 = sqrt((S - ms*ms) * (S - md*md));
  totcoef = 3.8937966E8 * lambda34 /(32.0 * M_PI * lambda12 * S);        
  Pout=lambda34/(2*sqrt_S);

/* fill momenta of particles */
  
  for(i=0;i<16;i++) pvec[i]=0;
  
  p1[0]=  sqrt(Pin*Pin + m1*m1);
  p1[3]=  Pin;
  p2[0]=  sqrt(Pin*Pin + m2*m2);
  p2[3]= -Pin;
  p3[0]=  sqrt(Pout*Pout + m3*m3);
  p4[0]=  sqrt(Pout*Pout + m4*m4);
  
  cos_step=0.1; /* step to calculate dsigma/dcos */
  sigmaTot=0.;  /* total cross section */
/* cycle */
  for(cos_fi=1-cos_step/2; cos_fi>-1; cos_fi -= cos_step)
  { double  DsigmaDcos;
    double  sin_fi=sqrt((1-cos_fi)*(1+cos_fi));
    int err=0;
    p3[3]=Pout*cos_fi; p4[3]=-p3[3];
    p3[2]=Pout*sin_fi; p4[2]=-p3[2];
    DsigmaDcos=totcoef*code->sqme(1,pvec,&err); 
 
    sigmaTot += DsigmaDcos*cos_step;
    
    printf("DsigmaDcos=%E\n",DsigmaDcos);
  }
  printf("sigmaTot=%E\n",sigmaTot);

  return 0;
}

