#include <sys/utsname.h>
#include "micromegas.h"
#include "micromegas_aux.h"
#include "micromegas_f.h"


double sWidth=0.01;
extern int  WIDTH_FOR_OMEGA;

static int neg_cs_flag;

static int NC=0;

static char ** inP;
static int  *  inAP;
static int  *  inG;
static int  *  inNum;

static numout ** code22; 
static int *inC;    /* combinatoric coefficients  NCxNC*/

static double *inDelta; /* inDelta[i]= (inMass[i]-M)/M; */
static double *inG_;    /* inG[i]*pow(1+inDelta[i],1.5);*/
static double **inMassAddress;
static double *inMass;  /* masses */
static int *sort;

aChannel* omegaCh=NULL;
aChannel* vSigmaTCh=NULL;

static int LSP;

double M1=0,M2=0;
static double DeltaXf;

static double pmass[4];

#define XSTEP 1.1
static double eps=0.001; /* precision of integration */

static int nsub=1;       /*  subprocess number  */

static double MassCut;

static double xf_;   /* to pass the Xf argument   */

#define MPlank 1.22E19 /*GeV*/
#define IMPROVE

static long PDGnum[4];
static  double sigma_simpson(double PcmIn)
{ double r;
  if(kin22(PcmIn,pmass)) return 0;
  _nsub_=nsub;

  r= simpson(dSigma_dCos,-1.,1.,0.3*eps);
  
  if(r<0) { neg_cs_flag=1;r=0;}
  
#ifdef IMPROVE
  improveCrossSection(PDGnum[0],PDGnum[1],PDGnum[2],PDGnum[3],PcmIn,&r);
#endif
 
/*
printf("PcmIN=%E sigma=%E\n",  PcmIn, r*3.8937966E8); 
*/
  return r;
}



static  double sigma_gauss(double PcmIn)
{ double r;
  if(kin22(PcmIn,pmass)) return 0.; 
  _nsub_=nsub;
  r=gauss(dSigma_dCos,-1.,1.,5);
  if(r<0) { neg_cs_flag=1;r=0;}

#ifdef IMPROVE
  improveCrossSection(PDGnum[0],PDGnum[1],PDGnum[2],PDGnum[3],PcmIn,&r);
#endif
  
  return r;
}  

static double (*sigma)(double)= sigma_gauss;

static double geff(double x)
{ double sum=0; int l;

  for(l=0;l<NC;l++)
  { int k=sort[l];
    double A=x*inDelta[k];
    if(A>15 || Mcdm +inMass[k] > MassCut) return sum;
    sum+=inG_[k]*exp(-A)*K2pol(1/(x+A));
  }
  return sum;
}


static void termod(double t, double * sqrt_gStar, double * heff)
{
  const  double table[276][3]=
  {
#include"geff_heff.tab"
  };

  if(t>= table[0][0])
  { if(sqrt_gStar) *sqrt_gStar=table[0][1];
    if(heff)             *heff=table[0][2];
  } else if(t<=0 )
  { if(sqrt_gStar) *sqrt_gStar=table[275][1];
    if(heff)             *heff=table[275][2];
  } else
  {
    int hi=0, lo=275, c;
    double tlo,thi;
    double slo,shi;

    for(c=138; lo-hi >1; c=(lo+hi)/2) if(table[c][0]< t) lo=c; else hi=c;

    tlo=table[lo][0];
    thi=table[hi][0];

    if(sqrt_gStar)
    {
       slo=table[lo][1];
       shi=table[hi][1];
       *sqrt_gStar=( (thi-t)*slo - (tlo-t)*shi)/(thi-tlo);
    }
    if(heff)
    {  slo=table[lo][2];
       shi=table[hi][2];
       *heff=( (thi-t)*slo - (tlo-t)*shi)/(thi-tlo);
    }
  }
}

static double y_pass;

static double weight_integrand(double v)
{  double x,gf;
   double sqrt_gStar;

   if(v==0.) return 0;
   x=xf_-3*log(v)/(y_pass-2);
   gf=geff(x);

   termod(Mcdm/x,&sqrt_gStar,NULL);
   return K1pol(1/(x*y_pass))*3*v*v*sqrt_gStar/(sqrt(x)*gf*gf*(y_pass-2));
}

static double weightBuff_x[1000];
static double weightBuff_y[1000];
static int inBuff=0;

static double weight(double y)
{ int i;
  double w;
  for(i=0;i<inBuff;i++) if(y==weightBuff_x[i]) return weightBuff_y[i];
  y_pass=y;
  w=  simpson(weight_integrand,0.,1.,0.3*eps);
  if(inBuff<1000){weightBuff_x[inBuff]=y; weightBuff_y[inBuff++]=w;}
  return w;
}

static int exi;

static double s_integrand( double u)
{  double z,y,sv_tot,w;
   double Xf_1;
   double ms,md,sqrtS,PcmIn;
   
   if(u==0. || u==1.) return 0.;

   z=1-u*u;
   y=2 +(DeltaXf - 3*log(z))/xf_;
   sqrtS=Mcdm*y;

   ms = M1 + M2;  if(ms>=sqrtS)  return 0;
   md = M1 - M2;
   PcmIn = sqrt((sqrtS-ms)*(sqrtS+ms)*(sqrtS-md)*(sqrtS+md))/(2*sqrtS);
   sv_tot=sigma(PcmIn);         

   if(exi) { w=weight(y); Xf_1=1;} else {w=K1pol(1/(xf_*y)); Xf_1=xf_;}

   return sqrt(2*Xf_1*y/M_PI)*y*(PcmIn*PcmIn/(Mcdm*Mcdm))*sv_tot*w*6*u*z*z;
}

static int Npow;

static double s_pow_integrand(double u)
{
   double ms,md,sqrtS;
   double z,y,sv_tot,pp,w,Xf_1;
   double PcmIn;
   
   if(u==0. || u==1.) return 0.;

   z=1-u*u;
   y=2+(DeltaXf - 3*log(z))/xf_;

   ms = M1 + M2;
   md = M1 - M2;
   sqrtS=(Mcdm*y);
   PcmIn = sqrt((sqrtS-ms)*(sqrtS+ms)*(sqrtS-md)*(sqrtS+md))/(2*sqrtS);
   pp= PcmIn* PcmIn;

   switch(Npow)
   { case 0: sv_tot=1; break;
     case 1: sv_tot= pp; break;
     case 2: sv_tot= pp*pp; break;
     case 3: sv_tot= pp*pp*pp; break;
   }

   if(exi){ w=weight(y); Xf_1=1;} else { w=K1pol(1/(xf_*y)); Xf_1=xf_;}
   return sqrt(Xf_1*y/(2*M_PI))*y*(PcmIn/Mcdm)*sv_tot*w*6*u*z*z;
}

static double m2u(double m) {return sqrt(1-exp((DeltaXf+xf_*(2-m/Mcdm))/3));}

typedef struct gridStr
{  int n;
   double ul[100];
   double ur[100];
   int pow[100];
}  gridStr;

static double u_max;

static gridStr   makeGrid(double mp,double wp)
{
  gridStr grd;

  int n=0,j;
  int pow_[6]={3,3,4,4,3,3};
  double c[5]={-8,-3,0,3,8};

  grd.ul[0]=0.;
  for(j=0;j<5;j++) if(mp+c[j]*wp>Mcdm*(2+DeltaXf/xf_))
  {  grd.ur[n]=m2u(mp+c[j]*wp);
     grd.pow[n]=pow_[j];
     grd.ul[n+1]=grd.ur[n];
     if( grd.ur[n]>u_max) { grd.ur[n]=u_max;  grd.n=n+1; return grd;}
     n++;
  }
  grd.ur[n]=u_max;
  grd.pow[n]=pow_[5];
  grd.n=n+1;
  return grd;
}

#ifdef DEBUG
static void printGrid(gridStr * grd)
{ int i;
  printf("~~~~~~~~~~~~~\n");
  for(i=0;i<grd->n;i++) printf("%E %E %d\n",grd->ul[i],grd->ur[i],grd->pow[i]);
  printf("~~~~~~~~~~~~~\n");
}
#endif

static gridStr  crossGrids(gridStr * grid1, gridStr * grid2)
{ gridStr grid;
  int n=0,i0=0,i1=0,i;
  grid.ul[0]=0.;
  while(i0<grid1->n && i1<grid2->n)
  { double d0= grid1->pow[i0]/(grid1->ur[i0]-grid1->ul[i0]);
    double d1= grid2->pow[i1]/(grid2->ur[i1]-grid2->ul[i1]);
    double d = ( d0>d1? d0:d1);
    int m=(grid1->pow[i0] > grid2->pow[i1]? grid1->pow[i0]:grid2->pow[i1]);

    if(grid1->ur[i0] < grid2->ur[i1]) { grid.ur[n]=grid1->ur[i0++];}
    else                              { grid.ur[n]=grid2->ur[i1++];}


    grid.pow[n]=0.999+d*(grid.ur[n]-grid.ul[n]);

    if(grid.pow[n]>m) grid.pow[n]=m;
    if(grid.pow[n]==1) grid.pow[n]=2;

    n++;
    grid.ul[n]=grid.ur[n-1];
  }
  grid.n=n;
  for(i=0;i<grid.n;i++) if(grid.ur[i]-grid.ul[i]>0.4
                         && grid.pow[i]<4)  grid.pow[i]=4;
  return grid;
}


static int  new_code(int k1,int k2);

static int testSubprocesses(void)
{
 static int first=1;
 int err,k1,k2,i,j;
 double *Q;
 double *GG;
 if(first)
 {
    first=0;
    if(createTableOddPrtcls())
    { printf("The model contains uncoupled odd patricles\n"); exit(10);}

    for(i=0,NC=0;i<Nodd;i++,NC++) 
        if(strcmp(OddPrtcls[i].name,OddPrtcls[i].aname))NC++;
        
    inP=(char**)malloc(NC*sizeof(char*));
    inAP=(int*)malloc(NC*sizeof(int));
    inG=(int*)malloc(NC*sizeof(int));
    inDelta=(double*)malloc(NC*sizeof(double)); 
    inG_=(double*)malloc(NC*sizeof(double));
    inMassAddress=(double**)malloc(NC*sizeof(double*));
    inMass=(double*)malloc(NC*sizeof(double));
    inNum= (int*)malloc(NC*sizeof(int));
    sort=(int*)malloc(NC*sizeof(int));

    code22=(numout**)malloc(NC*NC*sizeof(numout*));
    inC=(int*)malloc(NC*NC*sizeof(int)); 

      
    for(i=0,j=0;i<Nodd;i++)
    {  
       inP[j]=OddPrtcls[i].name;
       inNum[j]=OddPrtcls[i].NPDG;
       inG[j]=(OddPrtcls[i].spin2+1)*OddPrtcls[i].cdim;
       if(strcmp(OddPrtcls[i].name,OddPrtcls[i].aname))
       {
         inAP[j]=j+1;
         j++;
         inP[j]=OddPrtcls[i].aname;
         inG[j]=inG[j-1];
         inAP[j]=j-1;
         inNum[j]=-OddPrtcls[i].NPDG;
       } else inAP[j]=j;
       j++;
    }

    for(i=0;i<NC;i++) sort[i]=i;
    for(k1=0;k1<NC;k1++) for(k2=0;k2<NC;k2++) inC[k1*NC+k2]=-1;
    for(k1=0;k1<NC;k1++) for(k2=0;k2<NC;k2++) if(inC[k1*NC+k2]==-1)
    {  int kk1=inAP[k1];
       int kk2=inAP[k2];
       inC[k1*NC+k2]=1;
       if(inC[k2*NC+k1]==-1)   {inC[k2*NC+k1]=0;   inC[k1*NC+k2]++;}
       if(inC[kk1*NC+kk2]==-1) {inC[kk1*NC+kk2]=0; inC[k1*NC+k2]++;}
       if(inC[kk2*NC+kk1]==-1) {inC[kk2*NC+kk1]=0; inC[k1*NC+k2]++;}
    }

    for(k1=0;k1<NC;k1++) for(k2=0;k2<NC;k2++) code22[k1*NC+k2]=NULL;

    for(i=0,j=0;i<Nodd;i++)
    {
       inMassAddress[j]=varAddress(OddPrtcls[i].mass);
       if(!inMassAddress[j]) 
       { if(strcmp(OddPrtcls[i].mass ,"0")==0)
         { printf("Error: odd particle '%s' has zero mass.\n",OddPrtcls[i].name);
           exit(5);
         }  
         printf(" Model is not self-consistent:\n "
                " Mass identifier '%s' for particle '%s' is absent  among parameetrs\n",OddPrtcls[i].mass, OddPrtcls[i].name);
         exit(5);
       }

       if(strcmp(OddPrtcls[i].name,OddPrtcls[i].aname))
       {
         j++;
         inMassAddress[j]=inMassAddress[j-1];
       }
       j++;
    }
  }

  for(Q=NULL,GG=NULL,i=0;i<nModelVars;i++)  
  { if(strcmp(varNames[i],"Q")==0) Q=varValues+i;
    else if(strcmp(varNames[i],"GG")==0) GG=varValues+i;
  }  
  if(Q) *Q=100;
 
 err=calcMainFunc();
 if(err>0) return err;

 Mcdm=fabs(*(inMassAddress[0]));
 for(i=0;i<NC;i++) 
 { inMass[i]=fabs(*(inMassAddress[i]));
   if(Mcdm>inMass[i]) Mcdm=inMass[i];
 }

 if(Q) 
 { *Q=2*Mcdm;
    assignVal("Q",2*Mcdm);
    err=calcMainFunc();
    if(err>0) return err;
 }
 if(GG) *GG=parton_alpha(2*Mcdm/3.);
            
 for(i=0; i<NC-1;)
 { int i1=i+1;
   if(inMass[sort[i]] > inMass[sort[i1]])
   { int c=sort[i]; sort[i]=sort[i1]; sort[i1]=c;
     if(i) i--; else i++;
   } else i++;
 }

 LSP=sort[0];
 Mcdm=inMass[LSP];

 for(i=0;i<NC;i++)
 { inDelta[i]= (inMass[i]-Mcdm)/Mcdm;
   inG_[i]=inG[i]*pow(1+inDelta[i],1.5);
 }

  for(k1=0;k1<NC;k1++)  for(k2=0;k2<NC;k2++) if(code22[k1*NC+k2]) code22[k1*NC+k2]->init=0;
  cleanDecayTable();
return 0;
}


double aWidth(char * pName)
{  int dim;
   txtList LL;  
   return pWidth(pName,&LL,&dim);
}


 
int sortOddParticles(char * lsp)
{ int i,err;

  if(!modelNum)
  {
    int i,k,L;
    struct utsname buff;

    L=strlen(WORK);
    modelDir=malloc(L+15);  sprintf(modelDir,"%s/models",WORK);
    modelNum=1;
  
    calchepDir=malloc(strlen(WORK)+15);strcpy(calchepDir,WORK);
    for(i=L-1,k=2;i;i--) 
     {char ch=calchepDir[i]; calchepDir[i]=0; if(ch=='/') k--; if(k==0) break;}  
    strcat(calchepDir,"/CalcHEP_src");
    
    uname(&buff);
    compDir=malloc(strlen(WORK)+strlen(buff.nodename)+25);  
    strcpy(compDir,WORK);
    sprintf(compDir+strlen(compDir),"/_%s_%d_",buff.nodename,getpid());
      
    libDir=malloc(L+15); sprintf(libDir,"%s/so_generated",WORK);  
  }
  
  if(omegaCh) {free(omegaCh); omegaCh=NULL;}
  if(vSigmaTCh) {free(vSigmaTCh); vSigmaTCh=NULL;}
  if(vSigmaCh)  {free(vSigmaCh);  vSigmaCh=NULL; } 
  err=testSubprocesses();
  if(err)
  {
    strcpy(lsp,varNames[err]);
    return err;
  }

/*  
  if(Mcdm<0.1) 
  { sprintf(lsp,"Mcdm(%s)<0.1",inP[LSP]);
    return -1;
  } 
*/  
  if(sWidth>0) for(i=0;i< Nodd;i++) assignVal(OddPrtcls[i].width,sWidth*Mcdm);
  if(lsp) strcpy(lsp,inP[LSP]);
  return 0;
}

int  wimpPos(void) {return abs(pTabPos(inP[LSP]));}

char * wimpAnnLib(void)
{ static char out[20];
  int n=abs(pTabPos(inP[LSP]));
  if(strcmp(ModelPrtcls[n-1].name,ModelPrtcls[n-1].aname))
  sprintf(out,"omg_p%da%d",n,n); else sprintf(out,"omg_p%dp%d",n,n);
  return out;
}  

char * txtListOddParticles(void)
{ 
  static char * out=NULL;
  int i,len;

  if(out) return out;
  for(i=0,len=0;i<Nodd;i++) 
  { 
    len+=strlen(OddPrtcls[i].name)+1;
    if(strcmp(OddPrtcls[i].name,OddPrtcls[i].aname))
    len+=strlen(OddPrtcls[i].aname)+1;
  }
  out=malloc(len); out[0]=0;
  for(i=0;i<Nodd;i++) 
  { 
    if(i) strcat(out,",");
    strcat(out,OddPrtcls[i].name);
    if(strcmp(OddPrtcls[i].name,OddPrtcls[i].aname))
    { strcat(out,",");
      strcat(out,OddPrtcls[i].aname);
    } 
  }
  return out;
}

static int  new_code(int k1,int k2)
{
   char lib[40];
   char process[40];
   char lib1[12],lib2[12];
   numout*cc;
  
   pname2lib(inP[k1],lib1);
   pname2lib(inP[k2],lib2); 
    
   sprintf(lib,"omg_%s%s",lib1,lib2);
   sprintf(process,"%s,%s->2*x",inP[k1],inP[k2]);
   cc=getMEcode(0,ForceUG,process,NULL,txtListOddParticles(),lib);
   if(cc) 
   { (*cc->interface->twidth)=1;
        code22[k1*NC+k2]=cc;
   } else inC[k1*NC+k2]=0;

   return 0;
}   
   

static void gaussC2(double * c, double * x, double * f)
{
  double  A[2][2];
  int i,j;
  double det;
  double B[2];
    
  for(i=0;i<2;i++)
  { int l=1; for(j=0;j<2;j++) {A[i][j]=l*c[i+j]; l=-l;}  
     B[i]=l*c[2+i];
  }
  
  det=A[0][0]*A[1][1] - A[0][1]*A[1][0];
  
  f[0]= ( B[0]*A[1][1]-B[1]*A[0][1])/det;
  f[1]= (-B[0]*A[1][0]+B[1]*A[0][0])/det;
  
  det=sqrt(f[0]+f[1]*f[1]/4.);
   
  x[0]= -f[1]/2.-det;
  x[1]= -f[1]/2.+det;
 
  for(i=0;i<2;i++) { B[i]=c[i]; A[0][i]=1; }
  for(j=0;j<2;j++)   A[1][j]=A[0][j]*x[j];

  det= A[0][0]*A[1][1] - A[0][1]*A[1][0];
  
  f[0]= ( B[0]*A[1][1]-B[1]*A[0][1])/det;
  f[1]= (-B[0]*A[1][0]+B[1]*A[0][0])/det;
} 

static double aRate(double X, int average,int Fast, aChannel ** wPrc,int *NPrc)
{
  double Sum=0.;
  int i,l1,l2;
  int nPrc=0;
  char* pname[5];
  gridStr grid,grid1;
  double MassCutOut=MassCut+Mcdm*log(100.)/X;
  double Msmall,Mlarge;

  int nPrcTot=0;
  if(MassCutOut<Mcdm*(2+10/X)) MassCutOut=Mcdm*(2+10/X); 

  WIDTH_FOR_OMEGA=1;
  xf_=X;
  exi=average;

  if(wPrc) *wPrc=NULL;

  for(l1=0;l1<NC;l1++)
  { int k1=sort[l1]; if(Mcdm+inMass[k1]>MassCut) break;
  for(l2=0;l2<NC;l2++)
  {
    double Sumkk=0.;
    double x[2],f[2];
    double factor;
    int k2=sort[l2];
    CalcHEP_interface * CI;

    if(inMass[k1]+inMass[k2] > MassCut) break;

    if(inC[k1*NC+k2]<=0) continue;
    if(code22[k1*NC+k2]==NULL) new_code(k1,k2);
    if(inC[k1*NC+k2]<=0) continue;


    if(!code22[k1*NC+k2]->init)
    { numout * cd=code22[k1*NC+k2];
      CalcHEP_interface *cdi=cd->interface;
      for(i=1;i<=cdi->nvar;i++) if(cd->link[i]) cdi->va[i]=*(cd->link[i]);     
      if(  cdi->calcFunc()>0 ) {FError=1; WIDTH_FOR_OMEGA=0;  return -1;}
      cd->init=1;
    }

    if(wPrc)
    {  nPrcTot+=code22[k1*NC+k2]->interface->nprc;
       *wPrc=(aChannel*)realloc(*wPrc,sizeof(aChannel)*(nPrcTot+1));
    }

    sqme=code22[k1*NC+k2]->interface->sqme;
    DeltaXf=(inDelta[k1]+inDelta[k2])*X;
    inBuff=0;

    M1=inMass[k1];
    M2=inMass[k2];

    Msmall=M1>M2? M1-Mcdm*(1-sWidth): M2-Mcdm*(1-sWidth);
    Mlarge=M1>M2? M2+Mcdm*(1-sWidth): M1+Mcdm*(1-sWidth);

    u_max=m2u(MassCutOut);
    if(Fast)
    { 
      if(Fast==1)
      {  double c[4];

         for(Npow=0;Npow<4;Npow++) c[Npow]=simpson(s_pow_integrand, 0. ,1. ,1.E-4);
         gaussC2(c,x,f);
         for(i=0;i<2;i++){ x[i]=sqrt(x[i]); f[i]*=2*x[i]/Mcdm;}
      }else 
      {
         double c[2];
         for(Npow=0;Npow<2;Npow++) c[Npow]=simpson(s_pow_integrand, 0. ,1. ,1.E-4);
         x[0]= sqrt(c[1]/c[0]);
         f[0]= c[0]*2*x[0]/Mcdm;
      }
    }
    factor=inC[k1*NC+k2]*inG[k1]*inG[k2]*exp(-DeltaXf);
    CI=code22[k1*NC+k2]->interface;
    for(nsub=1; nsub<= CI->nprc;nsub++,nPrc++)
    { double u_min=0.;
      double a=0;

      for(i=0;i<4;i++)  pname[i]=CI->pinf(nsub,i+1,pmass+i,PDGnum+i);
      if(wPrc) 
      { (*wPrc)[nPrc].weight=0;
        for(i=0;i<4;i++) (*wPrc)[nPrc].prtcl[i]=pname[i];
      }
                     
      if(pmass[2]+pmass[3]>MassCutOut) continue;

      if( (pmass[2]>Mlarge && pmass[3]<Msmall)
        ||(pmass[3]>Mlarge && pmass[2]<Msmall))
           { *(CI->twidth)=1; *(CI->gtwidth)=1;}
      else { *(CI->twidth)=0; *(CI->gtwidth)=0;}
      *(CI->gswidth)=0;
                             
      if(pmass[2]+pmass[3] > pmass[0]+pmass[1])
      { double smin=pmass[2]+pmass[3];
        if((pmass[0]!=M1 || pmass[1]!=M2)&&(pmass[0]!=M2 || pmass[1]!=M1))
        { double ms=pmass[0]+pmass[1];
          double md=pmass[0]-pmass[1];
          double Pcm=sqrt((smin-ms)*(smin+ms)*(smin-md)*(smin+md))/(2*smin);
          smin=sqrt(M1*M1+Pcm*Pcm)+sqrt(M2*M2+Pcm*Pcm);
        }
        u_min=m2u(smin); 
      }else  u_min=0;
      
repeat:
      neg_cs_flag=0;

      if(!Fast) a=simpson(s_integrand,u_min,1.,eps); 
      else if(Fast!=1) a=f[0]*sigma(x[0]); else
      {
          int isPole=0;
          char * s;
          int m,w,n;
          double mass,width;

          for(n=1;(s=code22[k1*NC+k2]->interface->den_info(nsub,n,&m,&w));n++)
          if(m && w && strcmp(s,"\1\2")==0 )
          { mass=fabs(code22[k1*NC+k2]->interface->va[m]);
            width=code22[k1*NC+k2]->interface->va[w];
            if(mass<MassCutOut && mass+8*width > pmass[0]+pmass[1]
                            && mass+8*width > pmass[2]+pmass[3])
            { if((pmass[0]!=M1 || pmass[1]!=M2)&&(pmass[0]!=M2 || pmass[1]!=M1))
              { double ms=pmass[0]+pmass[1];
                double md=pmass[0]-pmass[1];
                double Pcm=sqrt((mass-ms)*(mass+ms)*(mass-md)*(mass+md))/(2*mass);
                mass=sqrt(M1*M1+Pcm*Pcm)+sqrt(M2*M2+Pcm*Pcm);
              }
              grid1=makeGrid(mass,width);
              if(isPole) grid=crossGrids(&grid,&grid1); else grid=grid1;
              isPole++;
            }
          }
          if(isPole==0)
          {  grid.n=1;
             grid.ul[0]=u_min;
             grid.ur[0]=u_max;
             grid.pow[0]=3;
          }

          if(grid.n==1 && pmass[0]+pmass[1]> 1.1*(pmass[2]+pmass[3]))
                a=f[0]*sigma(x[0])+f[1]*sigma(x[1]);
          else for(i=0;i<grid.n;i++)if(u_min<=grid.ur[i])
          {  
             double ul= u_min<grid.ul[i]? grid.ul[i]:u_min;
             double da=gauss(s_integrand,ul,grid.ur[i],grid.pow[i]);
             a+=da;             
          }
      }
      if(neg_cs_flag && *(CI->gswidth)==0)
      { *(CI->gswidth)=1;
         goto  repeat;
      }   
/*
printf("X=%.2E (%d) %.3E %s %s %s %s\n",X,average, a, pname[0],pname[1],pname[2],pname[3]);
*/
      Sumkk+=a;
      if(wPrc) (*wPrc)[nPrc].weight = a*factor;
    }
    Sum+=factor*Sumkk;
/*
printf("Sum=%E\n",Sum);
*/
  }
  }
  if(wPrc) 
  { for(i=0; i<nPrc;i++)  (*wPrc)[i].weight/=Sum;
    for(i=0;i<nPrc-1;)
    {  if((*wPrc)[i].weight >= (*wPrc)[i+1].weight) i++; 
       else
       {  aChannel buff;
          buff=(*wPrc)[i+1];(*wPrc)[i+1]=(*wPrc)[i];(*wPrc)[i]=buff;
          if(i)i--;else i++;
       }
    }          
    if(NPrc) *NPrc=nPrc; 
    (*wPrc)[nPrc].weight=0; for(i=0;i<4;i++) (*wPrc)[nPrc].prtcl[i]=NULL; 
  }  
  if(!average) { double gf=geff(X);  Sum/=gf*gf;}
/*
exit(1);
*/
  WIDTH_FOR_OMEGA=0;
  return Sum;
}


double vSigma(double T,double Beps ,int Fast)
{
    double X=Mcdm/T;
    assignVal("Q",2*Mcdm); 
    MassCut=Mcdm*(2-log(Beps)/X);    
    return  3.8937966E8*aRate(X, 0 ,Fast,&vSigmaTCh,NULL);
}


static double Yeq(double X)
{  double heff;
   termod(Mcdm/X,NULL,&heff);
   return (45/(4*M_PI*M_PI*M_PI*M_PI))*X*X*
                    geff(X)*sqrt(M_PI/(2*X))*exp(-X)/heff;
}
                          

struct {double*data;double xtop; int pow,size;} vSigmaGrid={NULL,0,0,0}; 

static void checkSgridUpdate(void)
{
  if(vSigmaGrid.pow==vSigmaGrid.size)
  { vSigmaGrid.size+=20;
    vSigmaGrid.data=(double*)realloc(vSigmaGrid.data,sizeof(double)*vSigmaGrid.size);
  }       
}

static double vSigmaI(double X, double Beps, int fast)
{ double XX;
  int i;
  if(vSigmaGrid.pow==0)
  { checkSgridUpdate();
    vSigmaGrid.pow=1;
    vSigmaGrid.xtop=X;
    MassCut=Mcdm*(2-log(Beps)/X);
    vSigmaGrid.data[0]= aRate(X,0,fast,NULL,NULL); 
    return vSigmaGrid.data[0];
  }

  while(X<= vSigmaGrid.xtop)
  { XX=vSigmaGrid.xtop/XSTEP;
    checkSgridUpdate();
    for(i=vSigmaGrid.pow;i;i--) vSigmaGrid.data[i]=vSigmaGrid.data[i-1];
    vSigmaGrid.xtop=XX;
    MassCut=Mcdm*(2-log(Beps)/XX);  
    vSigmaGrid.data[0]=aRate(XX,0,fast,NULL,NULL); 
    vSigmaGrid.pow++; 
  }

  for(XX=vSigmaGrid.xtop*pow(XSTEP,vSigmaGrid.pow-1);X>= XX;)
  { XX*=XSTEP;
    checkSgridUpdate();
    MassCut=Mcdm*(2-log(Beps)/XX);
    vSigmaGrid.data[vSigmaGrid.pow]=aRate(XX,0,fast,NULL,NULL);
    vSigmaGrid.pow++;
  }

  { double X0,X1,X2,sigmav0,sigmav1,sigmav2;
    i=log(X/vSigmaGrid.xtop)/log(XSTEP);
    if(i<0)i=0; 
    if(i>vSigmaGrid.pow-2) i=vSigmaGrid.pow-2;
    X1=vSigmaGrid.xtop*pow(XSTEP,i);   X2=X1*XSTEP;
    sigmav1=log(vSigmaGrid.data[i]);
    sigmav2=log(vSigmaGrid.data[i+1]);
    
    if(vSigmaGrid.pow==2)
    {
       return exp((sigmav1*(X-X2)-sigmav2*(X-X1))/(X1-X2));
    }
    if(i>0) { X0=X1/XSTEP; sigmav0=log(vSigmaGrid.data[i-1]); }
    else    { X0=X2*XSTEP; sigmav0=log(vSigmaGrid.data[i+2]); }
    return exp(  sigmav0*(X-X1)*(X-X2)/(X0-X1)/(X0-X2)
                +sigmav1*(X-X0)*(X-X2)/(X1-X0)/(X1-X2)
                +sigmav2*(X-X0)*(X-X1)/(X2-X0)/(X2-X1)
              );   
  }
}


static double dY(double X, double Beps,double fast)
{ double d, dYdX,Yeq0X, sqrt_gStar, vSig,cFactor=sqrt(M_PI/45)*MPlank*Mcdm;
  double epsY;
  d=X*0.001; MassCut=Mcdm*(2-log(Beps)/X); dYdX=(Yeq(X+d)-Yeq(X-d))/(2*d);
  Yeq0X=Yeq(X);
  epsY=deltaY/Yeq0X;
  Yeq0X/=X;
  termod(Mcdm/X,&sqrt_gStar,NULL);
  vSig=vSigmaI(X,Beps, fast);
  if(vSig==0){ FError=1; return 0;}
  return -dYdX/(vSig*cFactor*sqrt_gStar*Yeq0X*Yeq0X)/sqrt(1+epsY*epsY);
} 


static double darkOmega1(double * Xf,double Z1,double dZ1,int Fast,double Beps)
{
  double X = *Xf;
  double CCX=(Z1-1)*(Z1+1);
  double dCCX=(Z1-1+dZ1)*(Z1+1+dZ1)-CCX;
  double ddY;
  double dCC1,dCC2,X1,X2;
    
  sigma= (Fast)?sigma_gauss:sigma_simpson;

  if(Beps>=1) Beps=0.999;
  vSigmaGrid.pow=0;
  
  ddY=dY(X,Beps,Fast); if(FError || ddY==0)  return -1;
  if(fabs(CCX-ddY)<dCCX) { *Xf=X; MassCut=Mcdm*(2-log(Beps)/X); return Yeq(X)*sqrt(1+ddY);} 
   
  dCC1=dCC2=ddY-CCX; ;X1=X2=X; 
  while(dCC2>0) 
  {  
     X1=X2;
     dCC1=dCC2;
     X2=X2/XSTEP;
     X=X2;
     dCC2=-CCX+dY(X,Beps,Fast);
  }
             
  while (dCC1<0)
  {  
     X2=X1;
     dCC2=dCC1;
     X1=X1*XSTEP;
     X=X1;
     dCC1=-CCX+dY(X,Beps,Fast); 
  }
  for(;;)
  { double dCC;
    if(fabs(dCC1)<dCCX) 
      {*Xf=X1; MassCut=Mcdm*(2-log(Beps)/X1); return Yeq(X1)*sqrt(1+CCX+dCC1);}
    if(fabs(dCC2)<dCCX || fabs(X1-X2)<0.0001*X1) 
      {*Xf=X2; MassCut=Mcdm*(2-log(Beps)/X2); return Yeq(X2)*sqrt(1+CCX+dCC2);}
    X=0.5*(X1+X2); 
    dCC=-CCX+dY(X,Beps,Fast);
    if(dCC>0) {dCC1=dCC;X1=X;}  else {dCC2=dCC;X2=X;} 
  }
}

static double Beps_;
static int Fast_;

static void XderivLn(double x, double *Y, double *dYdx)
{
  double y=Y[0];
  double yeq, sqrt_gStar;
  
  termod(Mcdm/x,&sqrt_gStar,NULL);
  MassCut=Mcdm*(2-log(Beps_)/x); yeq=Yeq(x);
  if(y<yeq) *dYdx=0; else 
  {
    double epsY=deltaY/y; 
    *dYdx=-(Mcdm/x/x)*MPlank*sqrt_gStar*sqrt(M_PI/45)*vSigmaI(x,Beps_,Fast_)*(y*y-yeq*yeq)*sqrt(1+epsY*epsY);
  }
}


double darkOmegaFO(double * Xf_, int Fast, double Beps)
{
  double Yf,Yi;
  double  Z1=2.5;
  double  dZ1=0.05;
  double x;
  double Xf=25;
  
  if(omegaCh) {free(omegaCh); omegaCh=NULL;}
    
  if(Xf_) *Xf_=Xf; 
  assignVal("Q",2*Mcdm); 
  if(Beps>=1) Beps=0.999;
  Yf=  darkOmega1(&Xf, Z1, dZ1,Fast, Beps);
  if(Yf<0||FError) { return -1;}
  x=Xf;
  
  Yi=1/( (Mcdm/x)*sqrt(M_PI/45)*MPlank*aRate(x, 1,Fast, NULL,NULL) );
  if(!finite(Yi)||FError)  { return -1;}
  if(Xf_) *Xf_=Xf; 
  return  2.742E8*Mcdm/(1/Yf +  1/Yi); /* 2.828-old 2.755-new 2.742 next-new */
}


double darkOmega(double * Xf, int Fast, double Beps)
{
  double Yt,Yi,Xt=25;
  double Z1=1.1;
  double Z2=10; 
  int i;
  double Xf1;
  
  if(Xf)*Xf=Xt;
  assignVal("Q",2*Mcdm);
  if(Beps>=1) Beps=0.999;
  Beps_=Beps; Fast_=Fast;
  
  if(Z1<=1) Z1=1.1;

  Yt=  darkOmega1(&Xt, Z1, (Z1-1)/5,Fast, Beps);
  if(Yt<0||FError) {return -1;}
  Xf1=Xt;
  for(i=0; ;i++)
  { double X2=vSigmaGrid.xtop*pow(XSTEP,i+1);
    double y;

    if(Yt>=Z2*Yeq(Xt))  break;
    
    if(Xt>X2*0.999999) continue; 
    y=Yt;
    if(odeint(&y,1 , Xt , X2 , 1.E-3, (X2-Xt)/2, &XderivLn)){ printf("problem in solving diff.equation\n");      return -1;} 
    Yt=y;
    Xt=X2;
  }
  if(Xf) *Xf=0.5*(Xf1+Xt);
  Yi=1/( (Mcdm/Xt)*sqrt(M_PI/45)*MPlank*aRate(Xt,1,Fast,NULL,NULL));
  if(!finite(Yi)||FError)  return -1;
  if(deltaY==0.)
  { dmAsymm=0;
    return  2.742E8*Mcdm/(1/Yt  +  1/Yi); /* 2.828-old 2.755-new,2.742 -newnew */
  } else
  {  double a,f,z0,Y0;
     a=fabs(deltaY);
     f= (sqrt(Yt*Yt+a*a)-a)/
        (sqrt(Yt*Yt+a*a)+a)*exp(-2*a/Yi);
     z0=sqrt(f)*2*a/(1-f);
     Y0=sqrt(z0*z0+a*a);  
     dmAsymm=log((Y0+deltaY)/(Y0-deltaY));
     return 2.742E8*Mcdm*Y0;
  }   
}


double printChannels(double Xf ,double cut, double Beps, int prcn, FILE * f)
{ int i,nPrc,nform=log10(1/cut)-2;
  aChannel *wPrc;
  double Sum,s;
  
  if(omegaCh) {free(omegaCh); omegaCh=NULL;}
  
  MassCut=Mcdm*(2-log(Beps)/Xf);
  Sum=aRate(Xf, 1,1,&omegaCh,&nPrc)*(Mcdm/Xf)*sqrt(M_PI/45)*MPlank/(2.742E8*Mcdm);

  if(FError)     { return -1;}
  if(wPrc==NULL) { return  0;}
  if(nform<0)nform=0;
   
  if(f)
  {  int j;
     fprintf(f,"\nChannels which contribute to 1/(omega) more than %G%%.\n",100*cut );
     if(prcn) fprintf(f,"Relative contrubutions in %% are displyed\n");
        else  fprintf(f,"Absolut  contrubutions  are  displyed\n");
     for(i=0,s=0;i<nPrc;i++)  if(fabs(omegaCh[i].weight)>=cut)
     {  s+=omegaCh[i].weight;
        if(prcn)
        { if(cut <0.000001) fprintf(f,"%.1E%% ",100*omegaCh[i].weight);
          else              fprintf(f,"%*.*f%% ",nform+3,nform,
                                        100*omegaCh[i].weight);
        } else fprintf(f,"%.1E ",Sum*omegaCh[i].weight); 
        for(j=0;j<4;j++)
        {
           fprintf(f,"%s ",omegaCh[i].prtcl[j]);
           if(j==1) fprintf(f,"->");
           if(j==3) fprintf(f,"\n");
        }
     }
  }
  return 1/Sum;
}

static int strcmp_(char * n1, char *n2) { if( n1[0]=='*' &&  n1[1]==0) return 0; return strcmp(n1, n2);}

double oneChannel(double Xf,double Beps,char*n1,char*n2,char*n3,char*n4)
{ int j,nPrc;
  aChannel *wPrc;
  double Sum,res;
  
  MassCut=Mcdm*(2-log(Beps)/Xf);
  Sum=aRate(Xf, 1,1,&wPrc,&nPrc)*(Mcdm/Xf)*sqrt(M_PI/45)*MPlank/(2.742E8*Mcdm);
  
  if(FError)     { return -1;}
  if(wPrc==NULL) { return  0;}  

  for(res=0,j=0;j<nPrc;j++) 
  if( ( (strcmp_(n1,wPrc[j].prtcl[0])==0 && strcmp_(n2,wPrc[j].prtcl[1])==0) ||
        (strcmp_(n2,wPrc[j].prtcl[0])==0 && strcmp_(n1,wPrc[j].prtcl[1])==0)
      ) &&
      ( (strcmp_(n3,wPrc[j].prtcl[2])==0 && strcmp_(n4,wPrc[j].prtcl[3])==0) ||
        (strcmp_(n4,wPrc[j].prtcl[2])==0 && strcmp_(n3,wPrc[j].prtcl[3])==0)
      )
    )  {res+=wPrc[j].weight;} 
      
  free(wPrc); 
  return res;
}

void wimpannlib_(char * f_name, int len)
{
  char *c_name=wimpAnnLib();
  cName2f(c_name, f_name,len);
} 

int omegach_(int *i, double *w, char* name10)
{ int k,j;
  if(!omegaCh || *i<1)  return 0;
  for(k=0;k<*i;k++) if(omegaCh[k].weight==0) return 0;
  k--;
  *w=omegaCh[k].weight;
  for(j=0;j<40;j++) name10[j]=' '; 
  for(j=0;j<4;j++) strncpy(name10+j*10,omegaCh[k].prtcl[j],strlen(omegaCh[k].prtcl[j]));
  return 1;
}
