/*
 Copyright (C) 1997,2006, Alexander Pukhov 
*/
#include <math.h>
#include <stdio.h>
#include"1d_integration.h"

#define nErrMax 10

static double const X2[2]={2.113249E-01,7.886751E-01 };
static double const F2[2]={5.000000E-01,5.000000E-01 };
static double const X3[3]={1.127017E-01,5.000000E-01 ,8.872983E-01 };
static double const F3[3]={2.777778E-01,4.444444E-01 ,2.777778E-01 };
static double const X4[4]={6.943185E-02,3.300095E-01 ,6.699905E-01 ,9.305682E-01 };
static double const F4[4]={1.739274E-01,3.260726E-01 ,3.260726E-01 ,1.739274E-01 };
static double const X5[5]={4.691008E-02,2.307653E-01 ,5.000000E-01 ,7.692347E-01 ,9.530899E-01 };
static double const F5[5]={1.184634E-01,2.393143E-01 ,2.844445E-01 ,2.393143E-01 ,1.184634E-01 };
static double const X6[6]={3.376523E-02,1.693953E-01 ,3.806904E-01 ,6.193096E-01 ,8.306047E-01 ,9.662348E-01 };
static double const F6[6]={8.566223E-02,1.803808E-01 ,2.339570E-01 ,2.339570E-01 ,1.803808E-01 ,8.566225E-02 };
static double const X7[7]={2.544604E-02,1.292344E-01 ,2.970774E-01 ,5.000000E-01 ,7.029226E-01 ,8.707656E-01 ,9.745540E-01 };
static double const F7[7]={6.474248E-02,1.398527E-01 ,1.909150E-01 ,2.089796E-01 ,1.909150E-01 ,1.398527E-01 ,6.474248E-02 };



double gauss( double (*func)(double),double a,double b, int n)
{
        
  double ans=0;
  if(n<1) n=1;
  if(n>7) { printf(" 7 is a miximum number of points for Gauss integration (call with %d)\n",n); n=7;} 
  switch(n)
  {  int i;
    case 1: ans=(b-a)*func((a+b)/2);  break;
    case 2:
      for(i=0;i<n;i++) ans+=F2[i]*func(a+ (b-a)*X2[i]); break;
    case 3: 
      for(i=0;i<n;i++) ans+=F3[i]*func(a+ (b-a)*X3[i]); break;
    case 4:
      for(i=0;i<n;i++) ans+=F4[i]*func(a+ (b-a)*X4[i]); break;
    case 5:
      for(i=0;i<n;i++) ans+=F5[i]*func(a+ (b-a)*X5[i]); break;
    case 6:
      for(i=0;i<n;i++) ans+=F6[i]*func(a+ (b-a)*X6[i]); break;
    case 7:
      for(i=0;i<n;i++) ans+=F7[i]*func(a+ (b-a)*X7[i]); break;      
    default: 
      return 0;
  }
  return ans*(b-a);                       
 }



static void r_gauss( double(*func)(double),double a,double b, 
double eps, double * aEps, double * ans, double * aAns, int* N, int depth, int * err)
{
  int i,n;

  double s1,s2,s2a,s3,s3a,e_err,d=b-a;
  
  if(*N<0)     { *err=2; return;}
  if(depth>50) { *err=3;  printf("gauss345: depth>50 for [%e %e]\n", a,b);     return;}
  
  for(n=0,s1=0;n<3;n++) s1+=F3[n]*func(a+ d*X3[n]); s1*=d; *N-=3;
  for(n=0,s2=0,s2a=0;n<4;n++) {double  f=F4[n]*func(a+ d*X4[n]); s2+=f;s2a+=fabs(f);} 
  s2*=d; s2a*=fabs(d);*N-=4;
 
  if(!isfinite(s1) || ! isfinite(s2)) { *err=1; ; return;}  

  e_err=eps*s2a;
 
  if( fabs(s1-s2) <= 30*e_err)
  { 
    for(n=0,s3=0,s3a=0;n<5;n++) { double f=F5[n]*func(a+ d*X5[n]); s3+=f; s3a+=fabs(f); } 
    if(!isfinite(s3)) {  *err=1; return;} 
    s3*=d; s3a*=fabs(d); *N-=5;
    if(fabs(s3-s2) <= e_err) 
    { *ans+=s3;
      *aAns+=s3a;
      return;
    }
  }    

  if(fabs(s1-s2) <= 0.1*(*aEps)) 
  {  *ans+=s2;
     *aAns+=s2a;   
     *aEps -= fabs(s1-s2);
     return;
  }
         
  r_gauss(func,a,(a+b)/2,eps,aEps,ans,aAns,N,depth+1,err);
  if(*err) return;
  r_gauss(func,(a+b)/2,b,eps,aEps,ans,aAns,N,depth+1,err);
}   

double gauss345( double (*func)(double),double a,double b, double eps,int * err_code)
{
  double aEps; /* absolute error  */
  int n,k,err=0;	

  if(a==b) return 0;
  if(err_code) *err_code=0;

  for(n=0,aEps=0;n<4;n++) aEps+=F4[n]*fabs(func(a+ (b-a)*X4[n]));
  
  if(!isfinite(aEps)) { if(err_code) *err_code=1; else printf("gauss345: NaN in integrand\n"); 
                        return 0; 
                      }  
  if(aEps==0.)        return 0;
                      
  eps=eps/2;
  aEps = eps*aEps*fabs(b-a);


  for(k=0;;k++)
  {  double ans=0., aAns=0., aEps0=aEps;
     int N=50000*pow(2., (-log10(eps)-2)/2.);
     r_gauss(func,a,b,eps,&aEps,&ans,&aAns,&N,0,&err);
//printf("k=%d  aEps0=%E aEps=%E aAns=%E \n",k, aEps0, aEps, aAns);
     if(err) { if(err_code) *err_code=err; else 
               switch(err)
               { case 1: printf("gauss345: NaN in integrand\n"); break;
                 case 2: printf("gauss345: Too many points need for integration\n"); break;
                 case 3: printf("gauss345: Too deep recursion\n"); break;
               }
               return gauss(func,a,b,7);                                       
             }           
     if(aEps0-aEps < eps*aAns || k)   return ans; 
     aEps=aAns*eps;
  }
}

double gauss_arg( double (*func)(double,void*),void*par,double a,double b,  int n)
{
                
  double ans=0;
  if(n<1) n=1;
  if(n>7) { n=7; printf(" 7 is a miximum number of points for Gauss integration\n");} 
  switch(n)
  {  int i;
    case 1: ans=(b-a)*func((a+b)/2,par);  break;
    case 2:
      for(i=0;i<n;i++) ans+=F2[i]*func(a+ (b-a)*X2[i],par); break;
    case 3: 
      for(i=0;i<n;i++) ans+=F3[i]*func(a+ (b-a)*X3[i],par); break;
    case 4:
      for(i=0;i<n;i++) ans+=F4[i]*func(a+ (b-a)*X4[i],par); break;
    case 5:
      for(i=0;i<n;i++) ans+=F5[i]*func(a+ (b-a)*X5[i],par); break;
    case 6:
      for(i=0;i<n;i++) ans+=F6[i]*func(a+ (b-a)*X6[i],par); break;
    case 7:
      for(i=0;i<n;i++) ans+=F7[i]*func(a+ (b-a)*X7[i],par); break;      
  }
  return ans*(b-a);                       
 }


static void r_simpson( double(*func)(double),double * f,double a,double b, 
double eps, double * aEps, double * ans, double * aAns, int *N, int depth, int*nErr)
{
  double f1[5];
  int i;
  double s1,s2,s3,e_err;

  if(*N<0)     { *nErr=2; return;}
  if(depth>50) { *nErr=3; return;}
   
//printf("a=%E b=%E d=%d eps=%E \n",a,b,depth,eps);


  s1=(f[0]+4*f[4]+f[8])/6;  
  s2=(f[0]+4*f[2]+2*f[4]+4*f[6]+f[8])/12;
  s3=(f[0]+4*f[1]+2*f[2]+4*f[3]+2*f[4]+4*f[5]+2*f[6]+4*f[7]+f[8])/24;

  if(!isfinite(s3)){ *nErr=1;  return;}
  
  e_err=eps*fabs(s3);
  i=0;
  if( ( fabs(s3-s2) <= e_err && fabs(s3-s1) <= 16*e_err)) i=1; else
  if( fabs((s3-s2)*(b-a)) <= 0.1*(*aEps) && fabs((s3-s1)*(b-a)) <= 1.6*(*aEps)) 
  { i=1;  *aEps -= fabs((s3-s2)*(b-a));}


// 2^(-50) = 9E-16 

 int plost= (eps< 1E-14*(fabs(a)+fabs(b))/fabs(b-a));

 
  if(i)
  { *ans+=s3*(b-a);
    *aAns+=(fabs(f[0])+4*fabs(f[2])+2*fabs(f[4])+4*fabs(f[6])+fabs(f[8]))
          *fabs(b-a)/12;
    return;
  }        
  
  for(i=0;i<5;i++) f1[i]=f[4+i];
  for(i=8;i>0;i-=2)f[i]=f[i/2];

  for(i=1;i<8;i+=2) f[i]=(*func)(a+i*(b-a)/16);  *N-=4;

  r_simpson(func,f,a,(a+b)/2,eps,aEps,ans,aAns,N,depth+1,nErr);
  if(*nErr) return;
  for(i=0;i<5;i++) f[2*i]=f1[i];
  for(i=1;i<8;i+=2) f[i]=(*func)((a+b)/2+i*(b-a)/16); *N-=4;
  r_simpson(func, f,(a+b)/2,b,eps,aEps,ans,aAns,N,depth+1,nErr);
}

double simpson( double (*func)(double),double a,double b, double  eps, int *err)
{
  double f[9];
  double aEps; /* absolute error  */
  int i,j;	
  int N=50000*pow(2., (-log10(eps)-2)/2.);
  if(err) *err=0;
  int nErr=0;
  aEps=0;
  if(a==b) return 0;
  for(i=0;i<9;i++) { f[i]=(*func)(a+i*(b-a)/8); aEps +=fabs(f[i]); }
  if(aEps==0.)  return 0;
  eps=eps/2;
  aEps = eps*aEps*fabs(b-a)/9;
//printf("SIMPSON: esp=%E aEps=%E\n", eps,aEps);
  
  for(j=0;;j++)
  {  double ans=0., aAns=0.;
     r_simpson(func,f,a,b,eps,&aEps,&ans,&aAns,&N,0,&nErr);
     if(nErr) 
     {  if(err)  *err=nErr; else 
        switch(nErr)
        { case 1: printf("simpson: NaN in integrand\n"); break;
          case 2: printf("simpson: Too many points need for integration\n"); break;
          case 3: printf("simpson: Too deep recursion\n"); break;
        }   
        return gauss(func,a,b,7);                                             
     }
     
     if(5*aAns*eps > aEps || j>=2 ) return ans;
     for(i=0;i<9;i++)  f[i]=(*func)(a+i*(b-a)/8); N-=9;
     aEps=aAns*eps;
  }  
}


static void r_simpson_arg( double(*func)(double,void*par),double * f,double a,double b, void*par, 
                           double eps, double * aEps, double * ans, double * aAns, int*N, int depth,int*nErr)
{
  double f1[5];
  int i;
  double s1,s2,s3,e_err;

  if(*N<0)     {*nErr=2; return;}
  if(depth>50) {*nErr=3; /* printf("a=%E\n",a);*/ return;}
  s1=(f[0]+4*f[4]+f[8])/6;
  s2=(f[0]+4*f[2]+2*f[4]+4*f[6]+f[8])/12;
  s3=(f[0]+4*f[1]+2*f[2]+4*f[3]+2*f[4]+4*f[5]+2*f[6]+4*f[7]+f[8])/24;

  if(!isfinite(s3)){ *nErr=1; return;}
  
  e_err=eps*fabs(s3);
  i=0;
  if( ( fabs(s3-s2) <= e_err && fabs(s3-s1) <= 16*e_err)) i=1; else
  if( fabs((s3-s2)*(b-a)) <= 0.1*(*aEps) && fabs((s3-s1)*(b-a)) <= 1.6*(*aEps)) 
  { i=1;  *aEps -= fabs((s3-s2)*(b-a));}


  if(i)
  { *ans+=s3*(b-a);
    *aAns+=(fabs(f[0])+4*fabs(f[2])+2*fabs(f[4])+4*fabs(f[6])+fabs(f[8]))
          *fabs(b-a)/12;
     return;       
  }        
  
  for(i=0;i<5;i++) f1[i]=f[4+i];
  for(i=8;i>0;i-=2)f[i]=f[i/2];

  for(i=1;i<8;i+=2) f[i]=(*func)(a+i*(b-a)/16,par); *N-=4;

  r_simpson_arg(func,f,a,(a+b)/2,par,eps,aEps,ans,aAns,N,depth+1,nErr);
  if(*nErr) return;
  for(i=0;i<5;i++) f[2*i]=f1[i];
  for(i=1;i<8;i+=2) f[i]=(*func)((a+b)/2+i*(b-a)/16,par); *N-=4;
  r_simpson_arg(func, f,(a+b)/2,b,par,eps,aEps,ans,aAns,N,depth+1,nErr);
}

double simpson_arg( double (*func)(double,void*),void*par,double a,double b, double  eps, int *err)
{
  double f[9];
  double aEps; /* absolute error  */
  int i,j;	

  if(err) *err=0;

  aEps=0;
  if(a==b) return 0;
  for(i=0;i<9;i++) { f[i]=(*func)(a+i*(b-a)/8,par); aEps +=fabs(f[i]); }
  if(aEps==0.) { if(err) *err=0;  return 0;}
  eps=eps/2;
  aEps = eps*aEps*fabs(b-a)/9;
  int nErr=0;
  int N=50000*pow(2., (-log10(eps)-2)/2.);
  for(j=0;;j++)
  {  double ans=0., aAns=0.; 
     r_simpson_arg(func,f,a,b,par,eps,&aEps,&ans,&aAns,&N,0,&nErr);
     if(nErr)
     {  if(err)  *err=nErr; else
        switch(nErr)
        { case 1: printf("simpson_arg: NaN in integrand\n"); break;
          case 2: printf("simpson_arg: Too many points need for integration\n"); break;
          case 3: printf("simpson_arg: Too deep recursion\n"); break;
        }                                        
        return gauss_arg(func,par,a,b,7);
     }
                                                                        
     if(5*aAns*eps > aEps || j>=2 ) return ans;
     if(!isfinite(aAns)) return aAns;
     for(i=0;i<9;i++)  f[i]=(*func)(a+i*(b-a)/8,par);
     aEps=aAns*eps;
  }  
}


static void r_gauss_arg( double(*func)(double,void*),void*par, double a,double b, 
double eps, double * aEps, double * ans, double * aAns,int* N, int depth, int*err)
{
  int n;

  double s1,s2,s2a,s3,s3a,e_err,d=b-a;

  if(*N<0)      { *err=2;  return; }
  if(depth>50)  { *err=3;  return; }

  for(n=0,s1=0;n<3;n++) s1+=F3[n]*func(a+ d*X3[n],par); s1*=d; *N-=3;
  for(n=0,s2=0,s2a=0;n<4;n++) {double f=F4[n]*func(a+ d*X4[n],par); s2+=f; s2a+=fabs(f);} 
  s2*=d; s2a*=fabs(d);*N-=4;
  
  if(!isfinite(s1) || !isfinite(s2)) { *err=1; return;}
 
  e_err=eps*fabs(s2);
 
  if( fabs(s1-s2) <= 30*e_err)
  { 
    for(n=0,s3=0,s3a=0;n<5;n++) {double f=F5[n]*func(a+ d*X5[n],par); s3+=f; s3a+=fabs(f);} 
    if(!isfinite(s3)) { *err=1; return;} 
    s3*=d; s3a*=d; *N-=5;
    if(fabs(s3-s2) <= e_err) 
    { *ans+=s3;
      *aAns+=s3a;
      return;
    }
  }    

  if(fabs(s1-s2) <= 0.1*(*aEps)) 
  {  *ans+=s2;
     *aAns+=s2a;   
     *aEps -= fabs(s1-s2);
     return;
  }
  
  r_gauss_arg(func,par,a,(a+b)/2,eps,aEps,ans,aAns,N,depth+1,err);
  if(*err)  return;
  r_gauss_arg(func,par,(a+b)/2,b,eps,aEps,ans,aAns,N,depth+1,err);
}   

double gauss345_arg( double (*func)(double,void*),void*par,double a,double b, double eps,int * err_code)
{
  double d,s1,s2,s2a,s3,s3a,aEps; /* absolute error  */
  int k,n,err=0;	

  if(a==b) return 0;
  d=b-a;
  for(n=0,s1=0;n<3;n++) s1+=F3[n]*func(a+ d*X3[n],par); s1*=d;
  for(n=0,s2=0,s2a=0;n<4;n++) {double f=F4[n]*func(a+ d*X4[n],par); s2+=f; s2a+=fabs(f);}
       s2*=d; s2a*=fabs(d);
       
  if(err_code) *err_code=0;
  
  

  if(!isfinite(s1) || !isfinite(s2)  ) 
  { if(err_code) *err_code=1; else printf("gauss345_arg: NaN in integrand\n"); 
    return 0;
  }
       
  aEps= s2a*eps;
  if( fabs(s1-s2) <= 30*aEps)
  { 
    for(n=0,s3=0,s3a=0;n<5;n++) {double f=F5[n]*func(a+ d*X5[n],par); s3+=f; s3a+=fabs(f);} 
    if(!isfinite(s3)) 
    {  if(err_code) *err_code=1; 
       else printf("gauss345_arg: NaN in integrand\n"); 
       return 0;
    } 
    s3*=d; s3a*=fabs(d);
    aEps=eps*s3a;
    if(fabs(s3-s2) < aEps)  return s3; 
  }    

  aEps/=2;
  eps=eps/2;

  for(k=0;;k++)
  {  double ans=0., aAns=0., aEps0=aEps;
     int N=50000*pow(2., (-log10(eps)-2)/2.);
     r_gauss_arg(func,par,a,b,eps,&aEps,&ans,&aAns,&N,0,&err);
//printf("k=%d  aEps0=%E aEps=%E aAns=%E \n",k, aEps0, aEps, aAns);
     if(err) 
     {   if(err_code) *err_code=err;  else 
         switch(err)
         { case 1: printf("gauss345_arg: NaN in integrand\n");  break;
           case 2: printf("gauss345_arg: Too many points need for integration\n"); break;
           case 3: printf("gauss345_arg: Too deep recursion\n"); break;
         }
         return gauss_arg(func,par,a,b,7);
     }     
     if(aEps0-aEps < 1.5*eps*aAns || k)    return ans;
     aEps=aAns*eps;
  } 
}

