#include "../include/micromegas_aux.h"
#define smallDelta 1E-10
//#define ERROR_PLOT

int init22_par(par_22*arg, numout*cc, int nsub)
{  int i;
   arg->cc=cc;
   arg->nsub=nsub;
   int ntot,nin,nout; 
   if(!cc) return 1;
   procInfo1(cc,&ntot,&nin,&nout);
   if(nsub>ntot || nin!=2 || nout!=2) return 2;    
   for(i=0;i<4;i++) arg->pdg[i]=arg->cc->interface->pinfAux(arg->nsub,1+i,arg->spin2+i,NULL,NULL,(arg->ndf)+i);
//   printf("arg->ndf = %d %d %d %d\n", arg->ndf[0], arg->ndf[1], arg->ndf[2], arg->ndf[3]);
   arg->T=0;
   arg->err=0;
   return 0;
}


static double statFactor(par_22*arg)
{ double c=1;
  for(int i=0;i<4;i++)if(arg->eta[i])
  { double p=0; for(int j=0;j<3;j++) p+=arg->pvect[j+1+4*i]*arg->n[j]; 
    double E=arg->pvect[4*i]*arg->ch - p*arg->sh;
    c/=1-arg->eta[i]*exp(-E/arg->T); 
  } return c;  
} 


void mass22_par(par_22*arg)
{
   int i;
   numout*cc=arg->cc; 
   int nsub=arg->nsub;
   for(i=0;i<4;i++) cc->interface->pinf(nsub,1+i,arg->pmass+i,NULL);
     
   
   arg->s13=arg->s14=0;
   for(int n=1;;n++)
   {  
     int m,w,pnum;
     char*s=cc->interface->den_info(nsub,n,&m,&w,&pnum);
     if(!s) break;
     if(s[0]==1 && s[1]==3)
     { 
       double mass=cc->interface->va[m];      
       if(arg->s13==0)
       { arg->s13=1- 2*((arg->spin2[0]+arg->spin2[2])&1);
         arg->M13=mass;  
       } else if(arg->M13>mass) arg->M13=mass;
     }
     if(s[0]==1 && s[1]==4)
     { 
       double mass=cc->interface->va[m];
       if(arg->s14==0)
       { arg->s14=1- 2*((arg->spin2[0]+arg->spin2[3])&1);
         arg->M14=mass;  
       } else if(arg->M14>mass) arg->M14=mass;
     }
   }
   arg->sqrtSmin=arg->pmass[0]+arg->pmass[1];
   if(arg->sqrtSmin < arg->pmass[2]+arg->pmass[3]) arg->sqrtSmin=arg->pmass[2]+arg->pmass[3];
//   arg->s13=arg->s14=0;      
}


int  kin22_par(par_22*arg, REAL sqrtS,double GG)
{  
   REAL ms,md,lambda,Pin,Pout,m[5];
   int i;
   numout*cc=arg->cc;
   arg->err-=arg->err&ErrE;

   int nsub=arg->nsub;
//printf("kin22_par\n");   
   for(i=0;i<4;i++) m[i+1]=arg->pmass[i];
//printf("masses again %E %E %E %E\n", m[1],m[2],m[3],m[4]);
   ms= m[1]+m[2];
   md= m[1]-m[2];
   if(sqrtS<=ms) { arg->err+=ErrE; return 1;} 
   lambda = sqrt((sqrtS-ms)*(sqrtS+ms)*(sqrtS-md)*(sqrtS+md));
   Pin=arg->PcmIn = lambda/(2*sqrtS);

   ms= m[3]+m[4];
   md= m[3]-m[4];
   if(sqrtS<=ms)  { arg->err+=ErrE; return 1;} 
   lambda = sqrt((sqrtS-ms)*(sqrtS+ms)*(sqrtS-md)*(sqrtS+md));
   Pout=arg->PcmOut = lambda/(2*sqrtS);

   for(i=0;i<16;i++) arg->pvect[i]=0;   
   
   for(i=1;i<5;i++) m[i]*=m[i];
   
   arg->pvect[3] = Pin;
   arg->pvect[7] =-Pin;
   arg->pvect[0] = sqrt(Pin*Pin + m[1]);
   arg->pvect[4] = sqrt(Pin*Pin + m[2]);
      
   arg->pvect[8] = sqrt(Pout*Pout + m[3]);
   arg->pvect[12]= sqrt(Pout*Pout + m[4]);


   double dE,dP;
   if(arg->s13)
   { dE=(m[1]-m[3]-m[2]+m[4])/2/sqrtS; 
     dE=arg->pvect[0]-arg->pvect[8];
     dP=Pin-Pout; 
     arg->e13=(arg->M13*arg->M13+(dP-dE)*(dP+dE))/2/Pin/Pout;
     if(arg->e13<=0) {/* printf("pole in t__channel\n");*/ arg->err=arg->err|Errt; return 1;}
   }
   if(arg->s14)
   {  
     dE=(m[1]-m[4]-m[2]+m[3])/2/sqrtS;
     dP=Pin-Pout;
//     dP=(m[1]+m[2]-m[3]-m[4])/2 - (m[1]-m[2])*(m[1]-m[2]) -(m[3]-m[4])*(m[3]-m[4])/4/sqrtS/sqrtS;
//     dP/=(Pin+Pout);
     arg->e14=(arg->M14*arg->M14-dE*dE+dP*dP )/2/Pin/Pout;
     if(arg->e14<=0) {/* printf("pole in u_channel\n");*/ arg->err=arg->err|Errt; return 1;}
   }    
   arg->totFactor=Pout /(32.0*M_PI*Pin*sqrtS*sqrtS);
   arg->GG=GG;
   return 0;
}


static double  dSqme_dCosR_arg(double  dCos,void*arg_)
{
   double  r;
   par_22*arg=arg_;
   REAL cosf=(REAL)1 -(REAL)dCos;
   REAL sinf=sqrt(fabs((REAL)dCos*(1+cosf)));
   int err_code=0;

//printf("arg->PcmOut=%E\n", arg->PcmOut);   
   arg->pvect[11]= arg->PcmOut*cosf;
   arg->pvect[15]=-arg->pvect[11];
   arg->pvect[10]= arg->PcmOut*sinf;
   arg->pvect[14]=-arg->pvect[10];

//for(int i=0;i<16;i++) printf("pvec %d= %e\n",i, (double)arg->pvect[i]);      
   r = arg->cc->interface->sqme(arg->nsub,arg->GG,arg->pvect,NULL,&err_code);
//if(err_code)printf("sqme=%e err_codes=%d\n", r,err_code);   
   if(arg->T>0) r*=statFactor(arg);
   if(err_code)  arg->err= arg->err | err_code ;
   return r;
}

static double  dSqme_dCosL_arg(double  dCos,void*arg_)
{
   double  r;
   par_22*arg=arg_;
   REAL cosf=-(REAL)1 +(REAL)dCos;
   REAL sinf=sqrt(fabs((REAL)dCos*(1-cosf)));
   int err_code=0;

//printf("arg->PcmOut=%E\n", arg->PcmOut);   
   arg->pvect[11]= arg->PcmOut*cosf;
   arg->pvect[15]=-arg->pvect[11];
   arg->pvect[10]= arg->PcmOut*sinf;
   arg->pvect[14]=-arg->pvect[10];

//for(int i=0;i<16;i++) printf("pvec %d= %e\n",i, (double)arg->pvect[i]);      
   r = arg->cc->interface->sqme(arg->nsub,arg->GG,arg->pvect,NULL,&err_code);
//if(err_code)printf("sqme=%e err_codes=%d\n", r,err_code);   
   if(arg->T>0) r*=statFactor(arg);
   if(err_code)  arg->err= arg->err | err_code ;
   return r;
}




static double  dSqme_dKsi13B(double ksi,void*arg_)
{
   double  r;
   par_22*arg=arg_;
   double e=  arg->e13;
   if(ksi>=1/e) return 0;
   REAL sinf=sin(2*asin(sqrt(0.5*(1/ksi-e)))); 
   if(!isfinite(sinf)) sinf=0;
   REAL cosf=sqrt(fabs((1-sinf)*(1+sinf)));
   int err_code=0;

//printf("arg->PcmOut=%E\n", arg->PcmOut);   
   arg->pvect[10]=arg->PcmOut*sinf;  // p3[2]
   arg->pvect[11]=arg->PcmOut*cosf;  // p3[3]
   arg->pvect[14]=-arg->pvect[10];
   arg->pvect[15]=-arg->pvect[11];   
      
   r = arg->cc->interface->sqme(arg->nsub,arg->GG,arg->pvect,NULL,&err_code)/ksi/ksi ;
   if(arg->T>0) r*=statFactor(arg);
   arg->sqme_err=err_code;
   return r;
}

static double  dSqme_dKsi13F(double ksi,void*arg_)
{
   double  r;
   par_22*arg=arg_;
   double e=  arg->e13;

   if( ksi>=-log(e)) return 0; 
   REAL sinf=sin(2*asin(sqrt(0.5*fabs(exp(-ksi)-e))));
   if(!isfinite(sinf)) sinf=0; 
   REAL cosf=sqrt(fabs((1-sinf)*(1+sinf)));
   int err_code=0;

//printf("arg->PcmOut=%E\n", arg->PcmOut);   
   arg->pvect[10]=arg->PcmOut*sinf;  // p3[2]
   arg->pvect[11]=arg->PcmOut*cosf;  // p3[3]
   arg->pvect[14]=-arg->pvect[10];
   arg->pvect[15]=-arg->pvect[11];   
      
   r = arg->cc->interface->sqme(arg->nsub,arg->GG,arg->pvect,NULL,&err_code)*exp(-ksi) ;
   if(arg->T>0) r*=statFactor(arg);
   if(!isfinite(r)){
    for(int i=0;i<16;i++) printf(" %e\n",  (double)arg->pvect[i]);
    printf("!isfinite(r) err_code=%d  exp(-ksi)=%E ksi=%E Pin=%E Pout=%e   \n",err_code,(double)exp(-ksi), (double)(ksi), arg->PcmIn,arg->PcmOut);
   } 
   arg->sqme_err=err_code;
   return r;

}


static double  dSqme_dKsi14B(double ksi,void*arg_)
{
   double  r;
   par_22*arg=arg_;

   double e=  arg->e14;

   if(ksi>=1/e) return 0;
   
   REAL sinf=sin(2*asin(sqrt(0.5*(1/ksi-e)))); 
   if(!isfinite(sinf)) sinf=0;
   REAL cosf=sqrt(fabs((1-sinf)*(1+sinf)));
   int err_code=0;

// printf("arg->PcmOut=%E\n", arg->PcmOut);   
   arg->pvect[14]=arg->PcmOut*sinf;   // p4[2]
   arg->pvect[15]=arg->PcmOut*cosf;   // p4[3]   
   arg->pvect[10]=-arg->pvect[14];    // p3[2]
   arg->pvect[11]=-arg->pvect[15];    // p3[3]      
   r = arg->cc->interface->sqme(arg->nsub,arg->GG,arg->pvect,NULL,&err_code)/ksi/ksi ;
   if(arg->T>0) r*=statFactor(arg);
   arg->sqme_err=err_code;
   return r;
}

static double  dSqme_dKsi14F(double ksi,void*arg_)
{
   double  r;
   par_22*arg=arg_;
   double e=  arg->e14;
   if(ksi>=-log(e)) return 0;
   REAL sinf=sin(2*asin(sqrt(0.5*fabs(exp(-ksi)-e))));
   if(!isfinite(sinf)) sinf=0; 
   REAL cosf=sqrt(fabs((1-sinf)*(1+sinf)));
   int err_code=0;

//printf("arg->PcmOut=%E\n", arg->PcmOut);   
   arg->pvect[14]=arg->PcmOut*sinf;   // p4[2]
   arg->pvect[15]=arg->PcmOut*cosf;   // p4[3]   
   arg->pvect[10]=- arg->pvect[14];   // p3[2]
   arg->pvect[11]=-arg->pvect[15];    // p3[3]      
   r = arg->cc->interface->sqme(arg->nsub,arg->GG,arg->pvect,NULL,&err_code)*exp(-ksi);
   if(arg->T>0) r*=statFactor(arg);
   arg->sqme_err=err_code ;
   return r;
}

static double  poleSqme(double x0,double x1,double x2,double delta,double f0, double f1,double f2)
{
//  printf("x= %e %e %e\n", x0,x1,x2);
//  printf("f= %e %e %e\n", f0,f1,f2);
//  printf("delta=%e\n", delta);  
  if(f2-f1==0) {  return f0*x0;}    
  double r=(f1-f0)/(f2-f0),z1=(x0+delta)/(x1+delta),z2=(x0+delta)/(x2+delta);
//printf("r=%E (%E %E  %E) \n", r, (1/z1-1)/(1/z2-1),   (z1-1)/(z2-1),  (z1*z1-1)/(z2*z2-1));   
  double p,a,b;
  if(r>=(1/z1-1)/(1/z2-1)) p=-1;
  else if(r<=(z1*z1-1)/(z2*z2-1)) p=2;
  else 
  {      double p1=-1,p2=2;
     double a1=(pow(z1,p1)-1)/(pow(z2,p1)-1)-r , a2=(pow(z1,p2)-1)/(pow(z2,p2)-1)-r;
     if(a1>a2) { a=a1;a1=a2;a2=a; p=p1;p1=p2;p2=p;}

// printf("solve a1=%E a2=%E\n",a1,a2);
     
     while(fabs(p1-p2)>0.001) 
     { p=0.5*(p1+p2);
       a=(pow(z1,p)-1)/(pow(z2,p)-1)-r; 
       if(a<0) { a1=a;p1=p;} else { a2=a;p2=p;} 
     }   
  }   
//p=1;
//printf("p=%E\n", p);
  
  if(fabs(1-p)<0.1)
  {  p=1;
     b=(f0-f1)/(1/(x0+delta)-1/(x1+delta));
     a=f0-b/(x0+delta);
//printf("f0=%e f1=%E x0=%e x1=%e\n", f0,f1,x0,x1);     
//printf("dI=%E a=%E b=%e p=%e\n",a*x0+ b*log(x0/delta+1),a,b,p);
     return   a*x0+ b*log(x0/delta+1);
  } 
                
  b=(f0-f1)/(pow(x0+delta,-p)-pow(x1+delta,-p)); 
  a=f0-b*pow(x1+delta,-p);
//printf("dI=%E a=%E b=%e p=%e\n", a*x0+b/(1-p)*(pow(x0+delta,1-p) - pow(delta,1-p)),a,b,p);  
  return a*x0+b/(1-p)*(pow(x0+delta,1-p) - pow(delta,1-p));

}  

// cross section [pb]= sqmeInt*3.8937966E8*arg.PcmOut/(32.0*M_PI*arg.PcmIn*sqrtS*sqrtS))

//static double mone(double x,void*arg) { return dSqme_dCos_arg(-1+x,arg);}

double sqmeInt(par_22*arg, double eps)
{  
     
  double res1=0;
  int i;
  int err=0;
  double x0=0,f0,f1,f2,in;  
//  if(arg->s13)printf("arg->e13=%e\n", arg->e13);
  int err_tmp=arg->err;
  if(arg->s13!=0 && arg->e13 < smallDelta)
  { 
     x0=smallDelta;
     for(;;)
     { arg->err=0;
       f0=dSqme_dCosR_arg(x0,arg);
       f1=dSqme_dCosR_arg(x0/2,arg);
       f2=dSqme_dCosR_arg(x0/4,arg);
       if(arg->err&ErrC && arg->err&ErrD) { arg->err=err_tmp|ErrP;  return 0;}
       if(arg->err&ErrD){x0*=2; continue;}
       if(arg->err&ErrC){x0/=2; continue;}      
       in=poleSqme(x0,x0/2,x0/4,arg->e13,f0,f1,f2);
       break;
     }     
     res1+=in;
  }
  
  arg->err=err_tmp;
  res1+= simpson_arg(dSqme_dCosR_arg,arg,x0,1,eps,&err); 
  if(err)  arg->err=arg->err|ErrIA;
/*
  if(err)
  {
   printf("simpson error code %d. File 22_par.c line 363\n",err);
#ifdef ERROR_PLOT  
    char txt[100];
    sprintf(txt,"Process %s,%s -> %s %s", arg->cc->interface->pinf(1,1,NULL,NULL),
                                          arg->cc->interface->pinf(1,2,NULL,NULL),
                                          arg->cc->interface->pinf(1,3,NULL,NULL),
                                          arg->cc->interface->pinf(1,4,NULL,NULL));
    displayPlot("dSqme_dCos","1-cos",x0,1,0,1,"",0,dSqme_dCosR_arg,arg);
    exit(0);
#endif    
  }   
*/  
  
  if(arg->pdg[0]==arg->pdg[1]) return 2*res1;

  err_tmp=arg->err;
  x0=0;    
  
  if(arg->s14!=0 && arg->e14< smallDelta )
  {
    x0=smallDelta;
    for(;;)
    { arg->err=0;
      f0=dSqme_dCosL_arg(x0,arg);  
      f1=dSqme_dCosL_arg(x0/2,arg);      
      f2=dSqme_dCosL_arg(x0/4,arg);  
      
      if(arg->err&ErrC && arg->err&ErrD) { arg->err=err_tmp|ErrP;  return 0;}
      if(arg->err&ErrD){x0*=2; continue;}
      if(arg->err&ErrC){x0/=2; continue;}
                                                
      in=poleSqme(x0,x0/2,x0/4,arg->e14,f0,f1,f2);
      break; 
    }        
    res1+=in;
  }                                                                     
  arg->err=err_tmp;
  res1+=simpson_arg(dSqme_dCosL_arg,arg,x0,1,eps,&err); 
  if(err) arg->err=arg->err|ErrIA;   
/*
  if(err)
  {
   printf("simpson error code %d. File 22_par.c line 401\n",err);
#ifdef ERROR_PLOT  
    char txt[100];
    sprintf(txt,"Process %s,%s -> %s %s", arg->cc->interface->pinf(1,1,NULL,NULL),
                                          arg->cc->interface->pinf(1,2,NULL,NULL),
                                          arg->cc->interface->pinf(1,3,NULL,NULL),
                                          arg->cc->interface->pinf(1,4,NULL,NULL));
    displayPlot("dSqme_dCos","cos+1",x0,1,0,1,"",0,dSqme_dCosL_arg,arg);
    exit(0);
#endif    
  }
*/     
  return res1;
} 
 