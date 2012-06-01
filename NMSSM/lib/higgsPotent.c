#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include"../../sources/micromegas.h"
#include"../../sources/micromegas_aux.h"
#include"pmodel.h"
#include"pmodel_aux.h"

static int nCall;

/* SASHA */
static double  min1(double xx1, double xx2, double (*FQ)(double))
{  double x1=xx1,x2=xx2;
   double x3,x4,f1,f2,f3,f4;

   nCall=0; 
   f1=FQ(x1);
   f2=FQ(x2); 
 
   if(f1<f2) { x3=x1;f3=f1;f1=f2;x1=x2;f2=f3;x2=x3;}
   for(;;) 
   {  x3=x2+ (x2-x1);
      f3=FQ(x3); if(nCall>1000) return x2;
      if(f3<f2) { x1=x2;x2=x3; f1=f2;f2=f3;}  else  break;
   }

   for(;;)
   {
      double d1,d2,d3,ff4,xx4;
      if(f1>f3) {f4=f3;x4=x3;f3=f1;x3=x1;x1=x4;f1=f4;}
            
      if(x1==x2|| x2==x3) return x2;
      d1=f1/(x2-x1)/(x3-x1);
      d2=f2/(x1-x2)/(x3-x2);   
      d3=f3/(x1-x3)/(x2-x3); 
      xx4=(d1*(x2+x3)+d2*(x1+x3)+d3*(x1+x2))/(2*(d1+d2+d3));
      ff4=d1*(x2-xx4)*(x3-xx4)+d2*(x1-xx4)*(x3-xx4)+d3*(x1-xx4)*(x2-xx4);
      
           if(fabs(x1-x2)>10*fabs(x2-x3)) x4=(x1+9*x2)/10;
      else if(fabs(x2-x3)>10*fabs(x1-x2)) x4=(x3+9*x2)/10; 
      else x4=xx4;
         
      if(x4==x2) return x2;
      f4=FQ(x4);  if(nCall>1000) return x2;
      if( (x4-x2)*(x1-x2) >0)
      { if(f4<f2) {x3=x2; f3=f2; x2=x4;f2=f4;} else {x1=x4;f1=f4;}}       
      else 
      { if(f4<f2) {x1=x2; f1=f2; x2=x4;f2=f4;} else {x3=x4;f3=f4;}} 
         
      if( (x1==x2|| x2==x3) ||
          (fabs(x1-x3)<0.001)||
          (fabs(f2-ff4) < 0.001*fabs(ff4)) ||
          (fabs(x1-x3)<0.1 &&  fabs(f2-ff4) < 0.01*fabs(ff4))
        ) return x2;
   } 
}


static struct { double  tb, mu,Lambda,Kappa,aLambda, aKappa,g2,vev,mw; 
              } extpar;

static double la[10], las[10],lass;
static double Pa[2][2], Zh[3][3];
static double mh[3],mha[2],mhc;
static double orimm[3][3],orimmE[3][3];

static double khi2(double hls)
{
  double tb=extpar.tb;
  double hl=extpar.Lambda;
  double hk=extpar.Kappa;
  double hks=extpar.aKappa;
  double xvev; /*=extpar.mu/extpar.Lambda;*/  
  double g2=extpar.g2;
  double cb=1.0/sqrt(1.0+tb*tb);
  double sb=tb*cb;
  double vev=extpar.vev; 
  double sum;

  double neum[3][3];
  double mw=extpar.mw;
  int i1,i2,k;
  
  for(i1=0;i1<10;i1++) {la[i1]=0;las[i1]=0;}
  lass=0;

  xvev=-(orimm[0][1] -hl*vev*hls)/(2*hk*hl*vev); 
  extpar.mu=xvev*hl;
  

  neum[0][0]= hl*xvev*(hls+hk*xvev)/(sb*cb);

  la[5] =(orimm[0][0] -neum[0][0])/(-2.0*vev*vev);
  
  hks=-(orimm[1][1]-4*vev*vev*hk*hl*cb*sb-vev*vev*hl*hls/xvev*cb*sb)/(3*hk*xvev);
  extpar.aKappa=hks;    

   la[4]=-(mhc*mhc-mw*mw-xvev*hl*(hls+hk*xvev)/cb/sb+vev*vev*(hl*hl+la[5])
  )/(vev*vev);

  neum[0][0]= g2*vev*vev*cb*cb+hl*hls*tb*xvev+hk*hl*tb*xvev*xvev;  
  
  neum[1][1]= g2*vev*vev*sb*sb+hl*hls/tb*xvev+hk*hl/tb*xvev*xvev;              
  neum[2][2]= 4*hk*hk*xvev*xvev+vev*vev*hl*hls/xvev*cb*sb+hk*hks*xvev;
  neum[0][1]=neum[1][0]= -g2*vev*vev*cb*sb+2.0*vev*vev*cb*sb*hl*hl
               -xvev*hl*(hls+xvev*hk) +2*vev*vev*cb*sb*(la[4]+la[5]);
  neum[0][2]=neum[2][0]=vev*(2.0*hl*hl*xvev*cb-2*hk*hl*xvev*sb-hl*hls*sb);
  neum[1][2]=neum[2][1]=vev*(2.0*hl*hl*xvev*sb-2*hk*hl*xvev*cb-hl*hls*cb);



  la[1] =(orimmE[0][0]-neum[0][0])/(2.0*vev*vev*cb*cb);
  la[2] =(orimmE[1][1]-neum[1][1])/(2.0*vev*vev*sb*sb); 
  la[3] =(orimmE[0][1]-neum[0][1])/(2.0*vev*vev*sb*cb);
  las[1]=(orimmE[0][2]-neum[0][2])/(2.0*vev*xvev*cb);
  las[2]=(orimmE[1][2]-neum[1][2])/(2.0*vev*xvev*sb);
  lass  =(orimmE[2][2]-neum[2][2])/(2.0*xvev*xvev);

  sum=   (la[1])*(la[1]) + (la[2])*(la[2]) 
          + la[3]*la[3] +la[4]*la[4] + la[5]*la[5] 
          + las[1]*las[1] + las[2]*las[2] + lass*lass;
  nCall++;        
  return sum;
}


double higgspotent(double nothing)
{ int i,j,k;
  double ee,sw,cw,lmax;
  double MZ,alfEMZ,tb;
  char name[10]; 
  double Q;

  MZ=findValW("MZ");
  alfEMZ=findValW("alfEMZ");
  tb=findValW("tb");
  
  for(i=1;i<=2;i++) for(j=1;j<=2;j++) 
  { sprintf(name,"Pa%d%d",i,j); Pa[i-1][j-1]=findValW(name);}

  for(i=1;i<=3;i++) for(j=1;j<=3;j++) 
  { sprintf(name,"Zh%d%d",i,j); Zh[i-1][j-1]=findValW(name);}
{
double * pos[6]={ mha,mha+1,mh,mh+1,mh+2,&mhc};
int code[6]={36,46,25,35,45,37};
for(i=0;i<6;i++)
 if(slhaValExists("MO_HIGGS",1,code[i])) *(pos[i])=slhaVal("MO_HIGGS",0.,1,code[i]);
 else  *(pos[i])= slhaVal("MASS",0.,1,code[i]);
}
 
  for(i=0;i<2;i++) for(j=0;j<2;j++)  for(k=0,orimm[i][j]=0 ;k<2;k++) 
  orimm[i][j]+=Pa[k][i]*Pa[k][j]*mha[k]*mha[k];
  for(i=0;i<3;i++) for(j=0;j<3;j++)  for(k=0,orimmE[i][j]=0;k<3;k++)
  orimmE[i][j]+= Zh[k][i]*Zh[k][j]*mh[k]*mh[k];

  extpar.tb      = findValW("tb");
  extpar.Lambda  = findValW("Lambda");
  extpar.Kappa   = findValW("Kappa");
  extpar.aLambda = slhaVal("MO_HIGGS",0.,1,3);  
  extpar.vev     =findValW("vev");
/*  extpar.mu     = slhaVal("HMIX",0.,1,1); */
  ee=sqrt(4*M_PI*alfEMZ);
  sw=sin(asin(extpar.vev*sqrt(2.0)*ee/MZ)/2);
  cw=sqrt(1.0-sw*sw);
  extpar.mw=MZ*cw; 
  extpar.g2=ee*ee/(sw*cw)/(sw*cw)/2; 
  
  { double A,B,C,D,x1,x2,al1,al2;
    double hl,hk,tb,g2,vev,cb;
    hl=extpar.Lambda;
    hk=extpar.Kappa;
    tb=extpar.tb;
    g2=extpar.g2;
    vev=extpar.vev;
    cb=1.0/sqrt(1.0+tb*tb);
    A=3*hl*hk*tb;
    B=tb*orimm[0][1]/vev;
    C=g2*vev*vev*cb*cb-orimmE[0][0];
    D=B*B-4*A*C;
/*printf("A=%E B=%E C=%E D=%E\n",A,B,C,D);*/
    if(D<0)D=0; else D=sqrt(D);
  
    x1=(-B+D)/(2*A);
    x2=(-B-D)/(2*A);
    al1=(x1*(2*hk*hl*vev)+orimm[0][1])/(hl*vev);
    al2=(x2*(2*hk*hl*vev)+orimm[0][1])/(hl*vev);
    if(fabs(al1-extpar.aLambda) < fabs(al2-extpar.aLambda)) extpar.aLambda=al1; else extpar.aLambda=al2; 
  }
  
  extpar.aLambda=min1(extpar.aLambda, extpar.aLambda+1, khi2);
  khi2(extpar.aLambda);
  
  lmax=0;
  
  for(i=0;i<10;i++)
  {  if(lmax<fabs(la[i]))lmax=fabs(la[i]);
     if(lmax<fabs(las[i]))lmax=fabs(las[i]);
  }
  if(lmax<fabs(lass))lmax=fabs(lass);
  if(nCall>2000) return -1;
  else   return lmax;
}

double La(double ix)  { int i=ix+0.0001;return la[i];}
double Las(double ix) { int i=ix+0.0001;return las[i];}
double Lass(double ix){ int i=ix+0.0001;return lass;}

double mu(double nothing)      { return extpar.mu;}
double aLambda(double nothing) { return extpar.aLambda;}
double aKappa(double nothing)  { return extpar.aKappa;}
 
