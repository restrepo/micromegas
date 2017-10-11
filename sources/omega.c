#include <sys/utsname.h>
#include "micromegas.h"
#include "micromegas_aux.h"
#include "micromegas_f.h"
#include "../CalcHEP_src/include/rootDir.h" 


char* CDM1=NULL, *CDM2=NULL,*aCDM1=NULL,*aCDM2=NULL;
aChannel* omegaCh=NULL;
aChannel* vSigmaTCh=NULL;
REAL *Qaddress=NULL;

/*
static int NT;
static double *XX,*YY;
static void loadINTER(int N, double *x, double *y)
{ NT=N;XX=x;YY=y;}
double INTER(double x) { return polint4(x,NT,XX,YY);} 
*/

extern double cs23(numout*cc, int nsub, double Pcm, int ii3);


typedef  struct{ int virt,i3;double br,w[2];numout*cc23; int nTab; double *pcmTab; double *csTab;} processAuxRec;

typedef  processAuxRec* processAux;
static   processAux* code22Aux0, *code22Aux1,*code22Aux2;

static   processAux AUX;

double sWidth=0.01;

extern int  WIDTH_FOR_OMEGA;

static int neg_cs_flag;

static int NC=0;

static char ** inP;
static int  *  inAP;
static int  *  inG;
static int  *  inNum;

static numout ** code22_0;
static numout ** code22_1;
static numout ** code22_2;

static int *inC0;    /* combinatoric coefficients  NCxNC*/
static int *inC1; 
static int *inC2;
static numout*cc23=NULL;

static REAL **inMassAddress;
static double *inMass;  /* masses */
static int *sort;
static int gaussInt=1;


static int LSP;

static double M1=0,M2=0;

static REAL pmass[6];
static int pdg[4];

#define XSTEP 1.1
static double eps=0.001; /* precision of integration */

static double MassCut;

static double s3f_;   /* to pass the Xf argument   */

static double T_;

static int Tdim=0;
static double *t_=NULL, *heff_=NULL, *geff_=NULL, *s3_=NULL;

static int Z4ch( char *name)
{  if(name[0]!='~') return 0;
   if(name[1]!='~') return 1;
   return 2;
}


int loadHeffGeff(char*fname)
{  double *tabs[3];
   int nRec,nCol,i;
   nRec=readTable(fname,&nCol,tabs);
   if(nRec<=0) return nRec;
   if(nCol!=3)
   { for(i=0;i<nCol;i++) free(tabs[i]);
     return 0;
   }
   if(t_)free(t_);  if(heff_)free(heff_);   if(geff_)free(geff_);  if(s3_)free(s3_);
   t_=tabs[0];
   heff_=tabs[1];
   geff_=tabs[2];
   s3_=malloc(nRec*sizeof(double));
   for(i=0;i<nRec;i++) s3_[i]=t_[i]*pow(2*M_PI*M_PI/45.*heff_[i],0.3333333333333333);  
   Tdim=nRec;
   return nRec;
}

double gEff(double T) 
{
  if(Tdim==0) 
  {  char*fName=malloc(strlen(micrO)+40);
     int err;
     sprintf(fName,"%s/sources/data/%s",micrO,"std_thg.tab");       
     err=loadHeffGeff(fName);    
     free(fName);
  }

  if(T< t_[0]) T=t_[0];
  if(T> t_[Tdim-1]) T=t_[Tdim-1];
  return polint2(T,Tdim,t_,geff_);
}

double hEff(double T) 
{
  if(Tdim==0) 
  {  char*fName=malloc(strlen(micrO)+40);
     int err;
     sprintf(fName,"%s/sources/data/%s",micrO,"std_thg.tab");       
     err=loadHeffGeff(fName);    
     free(fName);
  }
  if(T< t_[0]) T=t_[0];
  if(T> t_[Tdim-1]) T=t_[Tdim-1];
  return polint2(T,Tdim,t_,heff_);
} 

#define MPlank 1.22E19 /*GeV*/
#define IMPROVE


static double sigma23(double PcmIn)
{  int l,l_;
   double r;
   double brV1,MV1, MV2,wV1,wV2;
   
   brV1=AUX[nsub22].br; 
   r=cs23(cc23,1,PcmIn,AUX[nsub22].i3)/brV1/3.8937966E8;
   
   l=AUX[nsub22].virt;
   l_=5-l;  
   if(AUX[nsub22].w[l_-2])
   {  double m1,m2,sqrtS;
      m1=pmass[0];
      m2=pmass[1];
      MV1=pmass[l];
      MV2=pmass[l_];
      wV1=AUX[nsub22].w[l-2];
      wV2=AUX[nsub22].w[l_-2];
      sqrtS=sqrt(PcmIn*PcmIn+m1*m1) + sqrt(PcmIn*PcmIn+m2*m2); 
      if(wV1*wV2>0) r*=decayPcmW(sqrtS,MV1,MV2,wV1,wV2,0)/decayPcmW(sqrtS,MV1,MV2,wV1,0,0);    
   }   

   if(r<1.E-15) r=1.E-15;
   return log(r*PcmIn);
}   


static  double sigma(double PcmIn)
{ double r; 

  if(AUX[nsub22].nTab>0 && PcmIn<=AUX[nsub22].pcmTab[AUX[nsub22].nTab-1])
  {  if(PcmIn<AUX[nsub22].pcmTab[0]) r= 0; else 
     { 
        r=exp(polint4(PcmIn,AUX[nsub22].nTab,AUX[nsub22].pcmTab,AUX[nsub22].csTab))/PcmIn;
     }  
  } 
  else
  {  if(kin22(PcmIn,pmass)) return 0.; 
     if(gaussInt) r=gauss(dSigma_dCos,-1.,1.,5); else r=simpson(dSigma_dCos,-1.,1.,0.3*eps);;
     if(r<0) { neg_cs_flag=1;r=0;}

     if((VZdecay||VWdecay) && (AUX[nsub22].w[0] || AUX[nsub22].w[1] )) 
     { double f; 
       double sqrtS=sqrt(PcmIn*PcmIn+pmass[0]*pmass[0]) + sqrt(PcmIn*PcmIn+pmass[1]*pmass[1]); 
       if( sqrtS-pmass[2]-pmass[3] < 15*(AUX[nsub22].w[0]+AUX[nsub22].w[1]))   
         f=decayPcmW(sqrtS,pmass[2],pmass[3],AUX[nsub22].w[0],AUX[nsub22].w[1],5)/decayPcm(sqrtS,pmass[2],pmass[3]);     
       else   
       {   f=1; 
           if(AUX[nsub22].w[0]) f-= AUX[nsub22].w[0]/pmass[2]/M_PI;
           if(AUX[nsub22].w[1]) f-= AUX[nsub22].w[1]/pmass[3]/M_PI;
       }
       r*=f;
     }
  }  
#ifdef IMPROVE
  improveCrossSection(pdg[0],pdg[1],pdg[2],pdg[3],PcmIn,&r);
#endif  
  return r;
}  



static double geffDM(double T)
{ double sum=0; int l;

  for(l=0;l<NC;l++)
  { int k=sort[l];
    double A=Mcdm/T*(inMass[k]-Mcdm)/Mcdm   ;
    if(A>15 || Mcdm +inMass[k] > MassCut) return sum;
    sum+=inG[k]*pow(inMass[k]/Mcdm,1.5)*exp(-A)*K2pol(1/(Mcdm/T+A));
  }
  return sum;
}

static double y_pass;


static double weight_integrand(double s3)
{  double x,gf;
   double T,heff,geff;

   if(s3==0.) return 0;
   T=polint2(s3,Tdim,s3_,t_);
   heff=polint2(s3,Tdim,s3_,heff_);
   geff=polint2(s3,Tdim,s3_,geff_);
   gf=geffDM(T);
   return K1pol(T/(Mcdm*y_pass))*exp((1/T-1/T_)*(2-y_pass)*Mcdm)*sqrt(Mcdm/T)*heff/sqrt(geff)/(gf*gf*s3);
}


static double weightBuff_x[1000];
static double weightBuff_y[1000];
static int inBuff=0;


static double weight(double y)
{ int i;
  double w;
  for(i=0;i<inBuff;i++) if(y==weightBuff_x[i]) { return weightBuff_y[i]; }
  y_pass=y;
  w=  simpson(weight_integrand,0.,s3f_,0.3*eps);
  if(inBuff<1000){weightBuff_x[inBuff]=y; weightBuff_y[inBuff++]=w;}
  return w;
}

static int exi;

static double s_integrand( double u)
{  double z,y,sv_tot,w;
   double Xf_1;
   double ms,md,sqrtS,PcmIn,res0;
   
   if(u==0. || u==1.) return 0.;

   z=1-u*u;
   sqrtS=M1+M2-3*T_*log(z);
   y=sqrtS/Mcdm;
   ms = M1 + M2;  if(ms>=sqrtS)  return 0;
   md = M1 - M2;
   PcmIn = sqrt((sqrtS-ms)*(sqrtS+ms)*(sqrtS-md)*(sqrtS+md))/(2*sqrtS);
   sv_tot=sigma(PcmIn);         

   res0=sqrt(2*y/M_PI)*y*(PcmIn*PcmIn/(Mcdm*Mcdm))*sv_tot*6*u*z*z;
   
   if(exi) { return res0*weight(sqrtS/Mcdm); } else return  res0*K1pol(T_/sqrtS)*sqrt(Mcdm/T_);
}

static int Npow;

static double s_pow_integrand(double u)
{
   double ms,md,sqrtS;
   double z,y,sv_tot,pp,w,Xf_1;
   double PcmIn,res0;
   
   if(u==0. || u==1.) return 0.;

   z=1-u*u;
   sqrtS=M1+M2-3*T_*log(z);
   y=sqrtS/Mcdm;
   ms = M1 + M2;
   md = M1 - M2;
   PcmIn = sqrt((sqrtS-ms)*(sqrtS+ms)*(sqrtS-md)*(sqrtS+md))/(2*sqrtS);
   pp= PcmIn* PcmIn;

   switch(Npow)
   { case 0: sv_tot=1; break;
     case 1: sv_tot= pp; break;
     case 2: sv_tot= pp*pp; break;
     case 3: sv_tot= pp*pp*pp; break;
   }
   res0=sqrt(y/(2*M_PI))*y*(PcmIn/Mcdm)*sv_tot*6*u*z*z;
   if(exi)  return res0*weight(y);  else  return res0*K1pol(T_/sqrtS)*sqrt(Mcdm/T_);
}

static double m2u(double m) {return sqrt(1-exp(((M1+M2 -m)/T_)/3));}

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
  for(j=0;j<5;j++) if(mp+c[j]*wp>M1+M2)
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

static int testSubprocesses(void)
{
  static int first=1;
  int err,k1,k2,i,j;
  CDM1=CDM2=NULL;
  Mcdm=Mcdm1=Mcdm2=0;
   
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
    inMassAddress=(REAL**)malloc(NC*sizeof(REAL*));
    inMass=(double*)malloc(NC*sizeof(double));
    inNum= (int*)malloc(NC*sizeof(int));
    sort=(int*)malloc(NC*sizeof(int));

    code22_0 = (numout**)malloc(NC*NC*sizeof(numout*));
    code22_1 = (numout**)malloc(NC*NC*sizeof(numout*));
    code22_2 = (numout**)malloc(NC*NC*sizeof(numout*));
            
    inC0=(int*)malloc(NC*NC*sizeof(int)); 
    inC1=(int*)malloc(NC*NC*sizeof(int));
    inC2=(int*)malloc(NC*NC*sizeof(int));
            
    code22Aux0=(processAux*) malloc(NC*NC*sizeof(processAux));
    code22Aux1=(processAux*) malloc(NC*NC*sizeof(processAux));
    code22Aux2=(processAux*) malloc(NC*NC*sizeof(processAux));
    
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
    for(k1=0;k1<NC;k1++) for(k2=0;k2<NC;k2++) inC0[k1*NC+k2]=-1;
    for(k1=0;k1<NC;k1++) for(k2=0;k2<NC;k2++) if(inC0[k1*NC+k2]==-1)
    {  int kk1=inAP[k1];
       int kk2=inAP[k2];
       inC0[k1*NC+k2]=1;
       if(inC0[k2*NC+k1]==-1)   {inC0[k2*NC+k1]=0;   inC0[k1*NC+k2]++;}
       if(inC0[kk1*NC+kk2]==-1) {inC0[kk1*NC+kk2]=0; inC0[k1*NC+k2]++;}
       if(inC0[kk2*NC+kk1]==-1) {inC0[kk2*NC+kk1]=0; inC0[k1*NC+k2]++;}
    }

    for(k1=0;k1<NC;k1++) for(k2=0;k2<NC;k2++) 
    { int L=k1*NC+k2;
       code22_0[L]=NULL;
       code22_1[L]=NULL;

       code22_2[L]=NULL;
       code22Aux0[L]=NULL;
       code22Aux1[L]=NULL;
       code22Aux2[L]=NULL;
       inC1[L]=inC0[L];
       inC2[L]=inC0[L];
    }    

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
    for(Qaddress=NULL,i=0;i<nModelVars;i++) if(strcmp(varNames[i],"Q")==0) Qaddress=varValues+i;
  }

  cleanDecayTable(); 
  if(Qaddress) *Qaddress=100;
 
  err=calcMainFunc();
  if(err>0) return err;
   
  if(Nodd==0) { printf("No odd particles in the model\n"); return -1; }

  Mcdm=fabs(*(inMassAddress[0]));
  for(i=0;i<NC;i++) 
  { inMass[i]=fabs(*(inMassAddress[i]));
    if(Mcdm>inMass[i]) Mcdm=inMass[i];
  }
  
  if(Qaddress) 
  { *Qaddress=2*Mcdm;
     err=calcMainFunc();
     if(err>0) return err;
  }
            
  GGscale=2*Mcdm/3;

  for(i=0; i<NC-1;)
   {  int i1=i+1;
      if(inMass[sort[i]] > inMass[sort[i1]])
      { int c=sort[i]; sort[i]=sort[i1]; sort[i1]=c;
        if(i) i--; else i++;
      } else i++;
   }

 
  for(i=0;i<NC;i++) 
  {
    if(Z4ch(inP[i])==1) { if(!CDM1) { Mcdm1=inMass[i]; CDM1=inP[i];} else if(Mcdm1>inMass[i]) { Mcdm1=inMass[i];CDM1=inP[i];} }
    if(Z4ch(inP[i])==2) { if(!CDM2) { Mcdm2=inMass[i]; CDM2=inP[i];} else if(Mcdm2>inMass[i]) { Mcdm2=inMass[i];CDM2=inP[i];} }
  }
  
  if(CDM1 && CDM2) for(i=0;i<Nodd;i++) if(Z4ch(OddPrtcls[i].name) != Z4ch(OddPrtcls[i].aname))
  { if(Mcdm1>Mcdm2) { Mcdm1=Mcdm2; CDM1=CDM2;} 
    CDM2=NULL;
    Mcdm2=0;
  }  
  if(CDM1){for(i=0;i<Nodd;i++) if(OddPrtcls[i].name==CDM1) { aCDM1=OddPrtcls[i].aname; break;}}
  else aCDM1=NULL;
  if(CDM2){for(i=0;i<Nodd;i++) if(OddPrtcls[i].name==CDM2) { aCDM2=OddPrtcls[i].aname; break;}}
  else aCDM2=NULL;
  
            
  if(CDM1){strcpy(CDM1_,CDM1); i=strlen(CDM1); } else i=0;
  for(;i<20;i++) CDM1_[i]=' ';
  if(CDM2){strcpy(CDM2_,CDM2); i=strlen(CDM2); } else i=0;
  for(;i<20;i++) CDM2_[i]=' ';

  LSP=sort[0];
  Mcdm=inMass[LSP];

  for(k1=0;k1<NC;k1++)  for(k2=0;k2<NC;k2++) 
  {  if(code22_0[k1*NC+k2]) code22_0[k1*NC+k2]->init=0;
     if(code22_1[k1*NC+k2]) code22_1[k1*NC+k2]->init=0;
     if(code22_2[k1*NC+k2]) code22_2[k1*NC+k2]->init=0;
  }             
  cleanDecayTable();

  for(k1=0;k1<NC;k1++) for(k2=0;k2<NC;k2++) if( code22_0[k1*NC+k2])
  { int  nprc=code22_0[k1*NC+k2]->interface->nprc;
    processAux prc=code22Aux0[k1*NC+k2];
    int n;
    for(n=0;n<=nprc;n++)
    { 
        prc[n].w[0]=0;
        prc[n].w[1]=0;
        prc[n].virt=0;
        prc[n].i3=0;  
        prc[n].br=0;  
        prc[n].cc23=NULL;
        prc[n].virt=0; 
        
        prc[n].nTab=0;   
        if(prc[n].pcmTab) { free(prc[n].pcmTab); prc[n].pcmTab=NULL;}
        if(prc[n].csTab)  { free(prc[n].csTab);  prc[n].csTab=NULL; }  
    }     
  }

  for(k1=0;k1<NC;k1++) for(k2=0;k2<NC;k2++) if( code22_1[k1*NC+k2])
  { int  nprc=code22_1[k1*NC+k2]->interface->nprc;
    processAux prc=code22Aux1[k1*NC+k2];
    int n;
    for(n=0;n<=nprc;n++)
    { 
        prc[n].w[0]=0;
        prc[n].w[1]=0;
        prc[n].virt=0;
        prc[n].i3=0;  
        prc[n].br=0;  
        prc[n].cc23=NULL;
        prc[n].virt=0; 
        
        prc[n].nTab=0;   
        if(prc[n].pcmTab) { free(prc[n].pcmTab); prc[n].pcmTab=NULL;}
        if(prc[n].csTab)  { free(prc[n].csTab);  prc[n].csTab=NULL; }  
    }     
  }

  
  return 0;
}

/*
double aWidth(char * pName)
{  txtList LL;  
   return pWidth(pName,&LL);
}
*/
 
int sortOddParticles(char * lsp)
{ int i,err;

  if(Tdim==0) 
  {  char*fName=malloc(strlen(micrO)+40);
     sprintf(fName,"%s/sources/data/%s",micrO,"std_thg.tab");       
     err=loadHeffGeff(fName);    
     free(fName);
  }

  if(!modelNum)
  {
    int i,k,L;
    struct utsname buff;

    L=strlen(WORK);
    modelDir=malloc(L+15);  sprintf(modelDir,"%s/models",WORK);
    modelNum=1;
  
    calchepDir=malloc(strlen(rootDir)+1);strcpy(calchepDir,rootDir);
    uname(&buff);
    compDir=malloc(strlen(WORK)+strlen(buff.nodename)+25);  
    strcpy(compDir,WORK);
    sprintf(compDir+strlen(compDir),"/_%s_%d_",buff.nodename,getpid());
      
    libDir=malloc(L+15); sprintf(libDir,"%s/so_generated",WORK);  
  }
  
  if(omegaCh)   {free(omegaCh);   omegaCh=NULL;}
  if(vSigmaTCh) {free(vSigmaTCh); vSigmaTCh=NULL;}
  if(vSigmaCh)  {free(vSigmaCh);  vSigmaCh=NULL; } 
  err=testSubprocesses();
  if(err)
  { 
    if(err>0) {strcpy(lsp,varNames[err]); printf("can not calculate parameter %s\n",varNames[err]);}
    else strcpy(lsp,"Nodd=0");
    printf("sortOddparticles err=%d\n",err); 
    return err;
  }

  if(lsp) strcpy(lsp,inP[LSP]);
  return 0;
}


char * OddParticles(int mode)
{ 
  static char * out[3]={NULL,NULL,NULL};
  int i,len;
  if(mode<0||mode>2) return NULL;
  
  if(out[mode]) return out[mode];

  for(i=0,len=0;i<Nodd;i++) 
  { if(mode==0 || Z4ch(OddPrtcls[i].name)==mode) 
    len+=strlen(OddPrtcls[i].name)+1;
    if(strcmp(OddPrtcls[i].name,OddPrtcls[i].aname))
    len+=strlen(OddPrtcls[i].aname)+1;
  }
  out[mode]=malloc(len); out[mode][0]=0;

  for(i=0;i<Nodd;i++) if(mode==0 || Z4ch(OddPrtcls[i].name)==mode)
  { 
    if(out[mode][0]) strcat(out[mode],","); 
    strcat(out[mode],OddPrtcls[i].name);
    if(strcmp(OddPrtcls[i].name,OddPrtcls[i].aname))
    { strcat(out[mode],",");
      strcat(out[mode],OddPrtcls[i].aname);
    } 
  }
  return out[mode];
}


char * EvenParticles(void)
{ 
  static char * out=NULL;
  int i,len;

  if(out) return out;
  for(i=0,len=0;i<nModelParticles;i++) if(ModelPrtcls[i].name[0]!='~')
  {
    len+=strlen(ModelPrtcls[i].name)+1;
    if(strcmp(ModelPrtcls[i].name,ModelPrtcls[i].aname))
    len+=strlen(ModelPrtcls[i].aname)+1;
  }
  out=malloc(len); out[0]=0;
  for(i=0;i<nModelParticles;i++)if(ModelPrtcls[i].name[0]!='~') 
  { 
    if(i) strcat(out,",");
    strcat(out,ModelPrtcls[i].name);
    if(strcmp(ModelPrtcls[i].name,ModelPrtcls[i].aname))
    { strcat(out,",");
      strcat(out,ModelPrtcls[i].aname);
    } 
  }
  return out;
}


static int  new_code(int k1,int k2, int ch)
{
   char lib[40];
   char process[400];
   char lib1[12],lib2[12];
   numout*cc;
  
   pname2lib(inP[k1],lib1);
   pname2lib(inP[k2],lib2); 
   sprintf(lib,"omg_%s%s",lib1,lib2);
   switch(ch)
   { case  0: sprintf(process,"%s,%s->AllEven,1*x{%s",inP[k1],inP[k2],EvenParticles()); break;
     case  1: sprintf(process,"%s,%s->AllOdd1,AllOdd1{%s",inP[k1],inP[k2],OddParticles(1));
             strcat(lib,"_1"); break;
     case -1: sprintf(process,"%s,%s->AllOdd2,AllOdd2{%s",inP[k1],inP[k2],OddParticles(2));
             strcat(lib,"_2"); break;
     case  2: sprintf(process,"%s,%s->AllOdd1,AllOdd2{%s{%s",inP[k1],inP[k2],OddParticles(1),OddParticles(2));
                  strcat(lib,"_12"); break;                        
   } 
   cc=getMEcode(0,ForceUG,process,NULL,NULL,lib);
   if(cc) 
   {  int nprc,n;
      processAux prc;
      *(cc->interface->twidth)=1;
      switch(ch)
      { case  0:  code22_0[k1*NC+k2]=cc; break;
        case  1: 
        case -1:  code22_1[k1*NC+k2]=cc; break;
        case  2:  code22_2[k1*NC+k2]=cc; break;
      }  
      nprc=cc->interface->nprc;
      prc=(processAux)malloc((nprc+1)*sizeof(processAuxRec));
      
      switch(ch)
      { case  0: code22Aux0[k1*NC+k2]=prc; break;
        case  1: 
        case -1: code22Aux1[k1*NC+k2]=prc; break;
        case  2: code22Aux2[k1*NC+k2]=prc; break;
      }  
      for(n=0;n<=nprc;n++)
      { 
        prc[n].w[0]=0;
        prc[n].w[1]=0;
        prc[n].virt=0;
        prc[n].i3=0; 
        prc[n].br=0; 
        prc[n].cc23=NULL;
        prc[n].nTab=0;   
        prc[n].pcmTab=NULL;
        prc[n].csTab=NULL;                          
      }   
   } else 
   {  switch(ch)
      { case  0:  inC0[k1*NC+k2]=0; break;
        case  1:  
        case -1:  inC1[k1*NC+k2]=0; break;
        case  2:  inC2[k1*NC+k2]=0; break;
      }
   }   
   return 0;
}


static int Ntab=0;
static double*Ttab=NULL;


static double *vs1100T = NULL;
static double *vs1120T = NULL;
static double *vs1122T = NULL;
static double *vs1210T = NULL;
static double *vs2200T = NULL;
static double *vs2211T = NULL;

static double *vs1110T = NULL;
static double *vs2220T = NULL;
static double *vs1112T = NULL;
static double *vs1222T = NULL;
static double *vs1220T = NULL;
static double *vs2210T = NULL;
static double *vs2221T = NULL;
static double *vs1211T = NULL;

static double *Y1T=NULL;
static double *Y2T=NULL;

double vs1120F(double T){ return  3.8937966E8*polint2Exp(T,Ntab,Ttab,  vs1120T);}
double vs2200F(double T){ return  3.8937966E8*polint2Exp(T,Ntab,Ttab,  vs2200T);}
double vs1100F(double T){ return  3.8937966E8*polint2Exp(T,Ntab,Ttab,  vs1100T);}
double vs1210F(double T){ return  3.8937966E8*polint2Exp(T,Ntab,Ttab,  vs1210T);}
double vs1122F(double T){ return  3.8937966E8*polint2Exp(T,Ntab,Ttab,  vs1122T);}
double vs2211F(double T){ return  3.8937966E8*polint2Exp(T,Ntab,Ttab,  vs2211T);}

double vs1110F(double T){ return  3.8937966E8*polint2Exp(T,Ntab,Ttab,  vs1110T);}
double vs2220F(double T){ return  3.8937966E8*polint2Exp(T,Ntab,Ttab,  vs2220T);}
double vs1112F(double T){ return  3.8937966E8*polint2Exp(T,Ntab,Ttab,  vs1112T);}
double vs1222F(double T){ return  3.8937966E8*polint2Exp(T,Ntab,Ttab,  vs1222T);}
double vs1220F(double T){ return  3.8937966E8*polint2Exp(T,Ntab,Ttab,  vs1220T);}
double vs2210F(double T){ return  3.8937966E8*polint2Exp(T,Ntab,Ttab,  vs2210T);}
double vs2221F(double T){ return  3.8937966E8*polint2Exp(T,Ntab,Ttab,  vs2221T);}
double vs1211F(double T){ return  3.8937966E8*polint2Exp(T,Ntab,Ttab,  vs1211T);}


double dY1F(double T){ return polint2Exp(T,Ntab,Ttab, Y1T) ;}
double dY2F(double T){ return polint2Exp(T,Ntab,Ttab, Y2T) ;}

double Y1F(double T){ return  dY1F(T)+Yeq1(T);}
double Y2F(double T){ return  dY2F(T)+Yeq2(T);}



   
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

static double aRate(double X, int average,int Fast, double * alpha, aChannel ** wPrc,int *NPrc)
{
  double Sum=0.;
  double Sum1=0;
  int i,l1,l2;
  int nPrc=0;
  char* pname[5];
  gridStr grid,grid1;  
  double MassCutOut=MassCut+Mcdm*log(100.)/X;
  double Msmall,Mlarge;

  int nPrcTot=0;
  if(MassCutOut<Mcdm*(2+10/X)) MassCutOut=Mcdm*(2+10/X); 
  WIDTH_FOR_OMEGA=1;

  T_=Mcdm/X;
  s3f_=polint2(T_,Tdim,t_,s3_);
  exi=average;

  if(wPrc) *wPrc=NULL;

  for(l1=0;l1<NC;l1++)
  { int k1=sort[l1]; if(Mcdm+inMass[k1]>MassCut) break;
  for(l2=0;l2<NC;l2++)
  {
    double Sumkk=0.;
    double Sum1kk=0;
    double x[2],f[2];
    double factor;
    int kk,k2=sort[l2];
    CalcHEP_interface * CI;

    if(inMass[k1]+inMass[k2] > MassCut) break;

    if(inC0[k1*NC+k2]<=0) continue;
    if(code22_0[k1*NC+k2]==NULL) new_code(k1,k2,0);
    if(inC0[k1*NC+k2]<=0) continue;

    if(!code22_0[k1*NC+k2]->init)
    {
      if(Qaddress && *Qaddress!=inMass[k1]+inMass[k2]) 
      { *Qaddress=inMass[k1]+inMass[k2];
         calcMainFunc();
      }       
      if(passParameters(code22_0[k1*NC+k2])) {FError=1; WIDTH_FOR_OMEGA=0;  return -1;}
      code22_0[k1*NC+k2]->init=1;
    }

    if(wPrc)
    {  nPrcTot+=code22_0[k1*NC+k2]->interface->nprc;
       *wPrc=(aChannel*)realloc(*wPrc,sizeof(aChannel)*(nPrcTot+1));
    }

    sqme22=code22_0[k1*NC+k2]->interface->sqme;
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
    factor=inC0[k1*NC+k2]*inG[k1]*inG[k2]*exp(-(M1+M2 -2*Mcdm)/T_);
    CI=code22_0[k1*NC+k2]->interface;
    AUX=code22Aux0[k1*NC+k2];
    for(nsub22=1; nsub22<= CI->nprc;nsub22++,nPrc++)
    { double u_min=0.,smin;
      double a=0;
      double K=0;
      for(i=0;i<4;i++) pname[i]=CI->pinf(nsub22,i+1,pmass+i,pdg+i);  
      if(pmass[0]<Mcdm/2 || pmass[1]<Mcdm/2) continue; 
      if(wPrc) 
      { (*wPrc)[nPrc].weight=0;
        for(i=0;i<4;i++) (*wPrc)[nPrc].prtcl[i]=pname[i];
      }
      smin=pmass[2]+pmass[3];
      cc23=NULL;
      
      if(VZdecay||VWdecay)
      {  int l,l_,nVV;        

         if(!AUX[nsub22].virt )  for(l=2;l<4;l++) if(pdg[l]==21 ||pdg[l]==22) { AUX[nsub22].virt=-1; break;}
         
         if(!AUX[nsub22].virt)
         {  int vd[4]={0,0,0,0};
            int c_a =  (pmass[0]>Mcdm) || (pmass[1]>Mcdm);

            if(c_a){ for(l=2;l<4;l++) if((pdg[l]==23 && VZdecay>1)   || (abs(pdg[l])==24 && VWdecay>1)) vd[l]=1;} 
            else    for(l=2;l<4;l++)
            { 
            
              if((pdg[l]==23 && VZdecay)     || (abs(pdg[l])==24 && VWdecay)) vd[l]=1;
            } 

            for(l=2;l<4;l++) if(vd[l]) break; 
            if(l<4)
            {  l_=5-l; 
               if(vd[l_])
               { nVV=2;
                 if(pmass[l_]>pmass[l]) { l=l_; l_=5-l;}
               } else nVV=1; 
               AUX[nsub22].virt=l;  
               AUX[nsub22].w[l-2]=pWidth(pname[l],NULL);
               if(abs(pdg[l_])>16 && pmass[l_]> 2) AUX[nsub22].w[l_-2]=pWidth(pname[l_],NULL);
               if(AUX[nsub22].w[l_-2] < 0.1) AUX[nsub22].w[l_-2]=0;
            } else  AUX[nsub22].virt=-1;
         }        
         if(AUX[nsub22].virt>0)
         {  l=AUX[nsub22].virt;
            l_=5-l; 
            if(pmass[0]+pmass[1] < smin+ 4*AUX[nsub22].w[l-2]   && pmass[l_]< MassCutOut)
            { 
              if(AUX[nsub22].cc23) cc23=AUX[nsub22].cc23; else
              {  double  brV1,wV1;
                 int i3W;
                 AUX[nsub22].cc23=xVtoxll(2,2,pname,pdg, l, &wV1, &brV1);
                 if(pdg[l]==pdg[l_]) brV1*=2;
                 AUX[nsub22].br=brV1;
                 cc23=AUX[nsub22].cc23; 
                 if(cc23)
                 {   double Pcm0,PcmMax;
                    *(cc23->interface->BWrange)=10000; 
                    *(cc23->interface->gswidth)=1;  
                    for(i3W=2;i3W<5;i3W++) if(strcmp(cc23->interface->pinf(1,i3W+1,NULL,NULL),pname[l_])==0) break; 
                    AUX[nsub22].i3=i3W;   
                    PcmMax=decayPcm(pmass[2]+pmass[3]+10*AUX[nsub22].w[l-2], pmass[0],pmass[1]);
                    if( pmass[0]+pmass[1]>=pmass[l_]) Pcm0=Mcdm*0.01; else 
                    Pcm0=1.01*decayPcm(pmass[l_], pmass[0],pmass[1]); 
                    buildInterpolation(sigma23, Pcm0,PcmMax, 0.02,1.E-5, &(AUX[nsub22].nTab), &(AUX[nsub22].pcmTab), &(AUX[nsub22].csTab));
/*
    printf("nTab=%d\n", AUX[nsub22].nTab);
   displayFunc(sigma23,Pcm0,PcmMax,"sigma23");
   loadINTER(AUX[nsub22].nTab, AUX[nsub22].pcmTab,AUX[nsub22].csTab); 
                displayFunc(INTER,Pcm0,PcmMax,"INTER"); 
*/                
                       
                 }
              }    
              if(cc23){ smin=pmass[l_]; smin=pmass[0]+pmass[1]+0.1;}
            } 
         }
      }
//if(cc23)  printf("23  %s %s -> %s %s\n", pname[0],pname[1],pname[2],pname[3]);
    
     
//if(abs(pdg[2])!=24 && abs(pdg[3])!=24) continue; 
      if(smin>=MassCutOut) continue; 
      if(cc23==NULL) 
      {                               
         if( (pmass[2]>Mlarge && pmass[3]<Msmall)
           ||(pmass[3]>Mlarge && pmass[2]<Msmall))
            { *(CI->twidth)=1; *(CI->gtwidth)=1;} else { *(CI->twidth)=0; *(CI->gtwidth)=0;}
      } 
      *(CI->gswidth)=0;
                             
      if(smin > pmass[0]+pmass[1])
      { 
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

//     if(abs(pdg[2])==24 && abs(pdg[3])==24 && average) 
//     {   displayFunc(s_integrand,u_min,1,"s_integrand");
//         printf("int=%E\n",simpson(s_integrand,u_min,1.,eps)); 
//     } 


      if(!Fast) a=simpson(s_integrand,u_min,1.,eps); 
      else if(Fast!=1) a=f[0]*sigma(x[0]); else
      {
          int isPole=0;
          char * s;
          int m,w,n;
          double mass,width;

          for(n=1;(s=code22_0[k1*NC+k2]->interface->den_info(nsub22,n,&m,&w));n++)
          if(m && w && strcmp(s,"\1\2")==0 )
          { mass=fabs(code22_0[k1*NC+k2]->interface->va[m]);
            width=code22_0[k1*NC+k2]->interface->va[w];
            if(mass<MassCutOut && mass+8*width > pmass[0]+pmass[1]
                            && mass+8*width > smin)
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
          if(cc23)
          {  double mass,width;
             mass=pmass[2]+pmass[3];
             width= AUX[nsub22].w[0]+ AUX[nsub22].w[1];
             if(mass-width>M1+M2)
             { double u=m2u(mass-1*width);
                if(u<0.96 && u>0.4)
                { 
                  grid1.n=2;
                  grid1.ul[0]=0;
                  grid1.ur[0]=grid1.ul[1]=u;
                  grid1.pow[0]=3;
                  grid1.ur[1]=1; 
                  if(u<0.6) grid1.pow[1]=5; else grid1.pow[1]=3;
                  if(isPole) grid=crossGrids(&grid,&grid1); else grid=grid1;
                  isPole++;                  
                }  
             }
          } 
          if(isPole==0)
          {  grid.n=1;
             grid.ul[0]=u_min;
             grid.ur[0]=u_max;
             grid.pow[0]=5;
          }
//for(i=0;i<grid.n;i++) printf(" (%E %E) ",grid.ul[i],grid.ur[i]); printf("\n");
/*          if(grid.n==1 && pmass[0]+pmass[1]> 1.1*(smin))
                a=f[0]*sigma(x[0])+f[1]*sigma(x[1]);
          else
*/          
           for(i=0;i<grid.n;i++)if(u_min<grid.ur[i])
          {  
             double ul= u_min<grid.ul[i]? grid.ul[i]:u_min;
             double da;
//                printf("gauss %E %E %d\n", ul, grid.ur[i], grid.pow[i]);
                da=gauss(s_integrand,ul,grid.ur[i],grid.pow[i]);
//printf(" da=%E (%E %E %d)\n", da,ul,grid.ur[i],grid.pow[i]);                
             a+=da;             
          }
      }
      if(neg_cs_flag && *(CI->gswidth)==0)
      { *(CI->gswidth)=1;
         goto  repeat;
      }   

// printf("X=%.2E (%d) %.3E %s %s %s %s\n",X,average, a, pname[0],pname[1],pname[2],pname[3]);


      for(kk=2;kk<4;kk++) if(pmass[kk]>2*Mcdm && pname[kk][0]!='~')
      {  txtList LL;
         double BrSm=1;
         pWidth(pname[kk],&LL);
         for(;LL;LL=LL->next)
         { double b;
           char proc[40];
           sscanf(LL->txt,"%lf %[^\n]",&b,proc);
           if( strchr(proc,'~'))BrSm-=b;
         }
         a*=BrSm;
      }
      
      if(pname[2][0]=='~' || pname[3][0]=='~' ) { a/=2; Sum1kk+=a;}  
      Sumkk+=a;
      if(wPrc) (*wPrc)[nPrc].weight = a*factor;
    }
    Sum+=factor*Sumkk;
    Sum1+=factor*Sum1kk;
    
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
  if(!average) { double gf=geffDM(Mcdm/X);  Sum/=gf*gf; Sum1/=gf*gf;   }
/*
exit(1);
*/
  WIDTH_FOR_OMEGA=0;
  if(alpha) {  *alpha=Sum1/Sum;  /*  printf("ALPHA=%E\n",*alpha);*/    }    
  return Sum;
}


double vSigma(double T,double Beps ,int Fast)
{
    double X=Mcdm/T;
    double res;
    if(assignVal("Q",2*Mcdm+T)==0) calcMainFunc();
    GGscale=(2*Mcdm+T)/3; 
    MassCut=Mcdm*(2-log(Beps)/X);    
    res= 3.8937966E8*aRate(X, 0 ,Fast,NULL,&vSigmaTCh,NULL);
    return res;
}


double Yeq(double T)
{  double heff;
   double X=Mcdm/T;
   heff=polint2(T,Tdim,t_,heff_);
   return (45/(4*M_PI*M_PI*M_PI*M_PI))*X*X*geffDM(T)*sqrt(M_PI/(2*X))*exp(-X)/heff;
}
                          

struct {double*data; double *alpha; double xtop; int pow,size;} vSigmaGrid={NULL,NULL,0,0,0}; 

static void checkSgridUpdate(void)
{
  if(vSigmaGrid.pow==vSigmaGrid.size)
  { vSigmaGrid.size+=20;
    vSigmaGrid.data=(double*)realloc(vSigmaGrid.data,sizeof(double)*vSigmaGrid.size);
    vSigmaGrid.alpha=(double*)realloc(vSigmaGrid.alpha,sizeof(double)*vSigmaGrid.size);
  }       
}

static double vSigmaI(double T, double Beps, int fast,double * alpha_)
{ double XX,alpha;
  int i,n;
  double X=Mcdm/T;
  if(vSigmaGrid.pow==0)
  { checkSgridUpdate();
    vSigmaGrid.pow=1;
    vSigmaGrid.xtop=X;
    MassCut=Mcdm*(2-log(Beps)/X);
    vSigmaGrid.data[0]= aRate(X,0,fast,&alpha,NULL,NULL);
    vSigmaGrid.alpha[0]=alpha;
    if(alpha_) *alpha_=alpha;     
    return vSigmaGrid.data[0];
  }
  
  n=log(X/vSigmaGrid.xtop)/log(XSTEP); 
  while(n<0)
  { XX=vSigmaGrid.xtop/XSTEP;
    checkSgridUpdate();
    for(i=vSigmaGrid.pow;i;i--)
    { vSigmaGrid.data[i]=vSigmaGrid.data[i-1];
      vSigmaGrid.alpha[i]=vSigmaGrid.alpha[i-1];
    }
    vSigmaGrid.xtop=XX;
    MassCut=Mcdm*(2-log(Beps)/XX);  
    vSigmaGrid.data[0]=aRate(XX,0,fast,&alpha,NULL,NULL);
    vSigmaGrid.alpha[0]=alpha; 
    vSigmaGrid.pow++;
    n++; 
  }
  while(n+2>vSigmaGrid.pow-1)
  { 
    XX=vSigmaGrid.xtop* pow(XSTEP,vSigmaGrid.pow)  ;
    checkSgridUpdate();
    MassCut=Mcdm*(2-log(Beps)/XX);
    vSigmaGrid.data[vSigmaGrid.pow]=aRate(XX,0,fast,&alpha,NULL,NULL);
    vSigmaGrid.alpha[vSigmaGrid.pow]=alpha;
    vSigmaGrid.pow++;
  }

  { double X0,X1,X2,X3,sigmav0,sigmav1,sigmav2,sigmav3,alpha0,alpha1,alpha2,alpha3;
    i=log(X/vSigmaGrid.xtop)/log(XSTEP);
    if(i<0)i=0; 
    if(i>vSigmaGrid.pow-2) i=vSigmaGrid.pow-2;
    X0=vSigmaGrid.xtop*pow(XSTEP,n-1); X1=X0*XSTEP;  X2=X1*XSTEP; X3=X2*XSTEP; 

    sigmav1=log(vSigmaGrid.data[n]);   alpha1=vSigmaGrid.alpha[n];
    sigmav2=log(vSigmaGrid.data[n+1]); alpha2=vSigmaGrid.alpha[n+1];
    sigmav3=log(vSigmaGrid.data[n+2]); alpha3=vSigmaGrid.alpha[n+2];
    X=log(X);X0=log(X0); X1=log(X1); X2=log(X2); X3=log(X3);
    
    if(n==0)
    {
   
    if(alpha_)
    { if(alpha1==0) *alpha_=0; else
      *alpha_=    alpha1*       (X-X2)*(X-X3)/        (X1-X2)/(X1-X3)
                 +alpha2*(X-X1)*       (X-X3)/(X2-X1)/        (X2-X3)
                 +alpha3*(X-X1)*(X-X2)       /(X3-X1)/(X3-X2)        ;
    }
    return exp( 
   +sigmav1*       (X-X2)*(X-X3)/        (X1-X2)/(X1-X3) 
   +sigmav2*(X-X1)*       (X-X3)/(X2-X1)/        (X2-X3) 
   +sigmav3*(X-X1)*(X-X2)       /(X3-X1)/(X3-X2)        );                           
    }
    sigmav0=log(vSigmaGrid.data[n-1]); alpha0=vSigmaGrid.alpha[n-1]; 
    
    if(alpha_)
    { if(alpha1==0) *alpha_=0; 
      else   *alpha_=  alpha0*       (X-X1)*(X-X2)*(X-X3)/        (X0-X1)/(X0-X2)/(X0-X3)
                       +alpha1*(X-X0)*       (X-X2)*(X-X3)/(X1-X0)/        (X1-X2)/(X1-X3)
                       +alpha2*(X-X0)*(X-X1)*       (X-X3)/(X2-X0)/(X2-X1)/        (X2-X3)
                       +alpha3*(X-X0)*(X-X1)*(X-X2)       /(X3-X0)/(X3-X1)/(X3-X2)        ;
    }
    return exp( 
    sigmav0*       (X-X1)*(X-X2)*(X-X3)/        (X0-X1)/(X0-X2)/(X0-X3)
   +sigmav1*(X-X0)*       (X-X2)*(X-X3)/(X1-X0)/        (X1-X2)/(X1-X3) 
   +sigmav2*(X-X0)*(X-X1)*       (X-X3)/(X2-X0)/(X2-X1)/        (X2-X3) 
   +sigmav3*(X-X0)*(X-X1)*(X-X2)       /(X3-X0)/(X3-X1)/(X3-X2)        );                        
  }
}


static double dY(double s3, double Beps,double fast)  
{ double d, dlnYds3,Yeq0X, sqrt_gStar, vSig,res;;
  double epsY,alpha;
  double T,heff,geff;
  T=polint2(s3,Tdim,s3_,t_);   
  heff=polint2(s3,Tdim,s3_,heff_);
  geff=polint2(s3,Tdim,s3_,geff_);
  MassCut=2*Mcdm-T*log(Beps);
  d=0.001*s3;  dlnYds3=( log(Yeq(polint2(s3+d,Tdim,s3_,t_)))- log(Yeq(polint2(s3-d,Tdim,s3_,t_))) )/(2*d);

  epsY=deltaY/Yeq(T);

//  sqrt_gStar=polint2(Mcdm/X,Tdim,t_,sqrt_gstar_);

  vSig=vSigmaI(T,Beps, fast,&alpha);
  if(vSig <=0) return 10;
  if(vSig==0){ FError=1; return 0;}
  res= dlnYds3/(pow(2*M_PI*M_PI/45.*heff,0.66666666)/sqrt(8*M_PI/3.*M_PI*M_PI/30.*geff)*vSig*MPlank
  *(1-alpha/2)*sqrt(1+epsY*epsY))/Yeq(T);
  res=fabs(res);
  if(res>10) return 10;
  return res;
} 


static double darkOmega1(double * Xf,double Z1,double dZ1,int Fast,double Beps)
{
  double X = *Xf;
  double CCX=(Z1-1)*(Z1+1);
  double dCCX=(Z1-1+dZ1)*(Z1+1+dZ1)-CCX;
  double ddY;
  double dCC1,dCC2,X1,X2;
    
  gaussInt= Fast? 1 : 0;

  if(Beps>=1) Beps=0.999;
  vSigmaGrid.pow=0;
  
  ddY=dY(polint2(Mcdm/X,Tdim,t_,s3_) ,Beps,Fast); 
  if(FError || ddY==0)  return -1;
  if(fabs(CCX-ddY)<dCCX) 
  { *Xf=X; MassCut=Mcdm*(2-log(Beps)/X); 
    return Yeq(Mcdm/X)*sqrt(1+ddY);
  } 
   
  dCC1=dCC2=ddY-CCX; ;X1=X2=X; 
  while(dCC2>0) 
  {  
     X1=X2;
     dCC1=dCC2;
     X2=X2/XSTEP;
     X=X2;
     dCC2=-CCX+dY(polint2(Mcdm/X,Tdim,t_,s3_),Beps,Fast);
     if(Mcdm/X>1.E5) return -1;
  }
             
  while (dCC1<0)
  {  
     X2=X1;
     dCC2=dCC1;
     X1=X1*XSTEP;
     X=X1;
     dCC1=-CCX+dY(polint2(Mcdm/X,Tdim,t_,s3_),Beps,Fast); 
  }
  for(;;)
  { double dCC;
    if(fabs(dCC1)<dCCX) 
      {*Xf=X1; MassCut=Mcdm*(2-log(Beps)/X1); return Yeq(Mcdm/X1)*sqrt(1+CCX+dCC1);}
    if(fabs(dCC2)<dCCX || fabs(X1-X2)<0.0001*X1) 
      {*Xf=X2; MassCut=Mcdm*(2-log(Beps)/X2); return Yeq(Mcdm/X2)*sqrt(1+CCX+dCC2);}
    X=0.5*(X1+X2); 
    dCC=-CCX+dY(polint2(Mcdm/X,Tdim,t_,s3_),Beps,Fast);
    if(dCC>0) {dCC1=dCC;X1=X;}  else {dCC2=dCC;X2=X;} 
  }
}

static double Beps_;
int Fast_=1;

static void XderivLn(double s3, double *Y, double *dYdx)
{
  double y=Y[0];
  double yeq, sqrt_gStar;
  double T,heff,geff;
  
//  s3=polint2(T,Tdim,t_,s3_);  
  
  heff=polint2(s3,Tdim,s3_,heff_);
  geff=polint2(s3,Tdim,s3_,geff_);
    T=polint2(s3,Tdim,s3_,t_);  
//  sqrt_gStar=polint2(T,Tdim,t_,sqrt_gstar_);
  
  MassCut=2*Mcdm -T*log(Beps_); yeq=Yeq(T);
  if(y<yeq) *dYdx=0; else 
  { double vSig;
    double alpha; 
    double epsY=deltaY/y; 
    vSig=vSigmaI(T,Beps_,Fast_,&alpha);
//printf("T=%E alpha=%E\n", Mcdm/x, alpha);     
    *dYdx=MPlank
    *pow(2*M_PI*M_PI/45.*heff,0.666666666666)/sqrt(8*M_PI/3.*M_PI*M_PI/30.*geff)
//    *sqrt_gStar*sqrt(M_PI/45)
    *vSig*(y*y-(1-alpha)*yeq*yeq-alpha*y*yeq)*sqrt(1+epsY*epsY);
//printf(" T=%E  y=%E   yeq=%E  epsY=%E  alpha=%E \n",T, y,  yeq, epsY, alpha);   
  }
}


double darkOmegaFO(double * Xf_, int Fast, double Beps)
{
  double Yf,Yi;
  double Z1=2.5;
  double dZ1=0.05;
  double Xf=25;

  if(CDM1==NULL) fracCDM2=1; else
  if(CDM2==NULL) fracCDM2=0; else 
  if(Mcdm1<Mcdm2) fracCDM2=0; else fracCDM2=1;
  
  if(omegaCh) {free(omegaCh); omegaCh=NULL;}
    
  if(Xf_) *Xf_=Xf; 
  if(assignVal("Q",2*Mcdm)==0) calcMainFunc();
  GGscale=2*Mcdm/3;   
  if(Beps>=1) Beps=0.999;
  
  Yf=  darkOmega1(&Xf, Z1, dZ1,Fast, Beps);
  if(FError||Xf<1||Yf<=0) {  return -1;}
 
  Yi=1/( (Mcdm/Xf)*sqrt(M_PI/45)*MPlank*aRate(Xf, 1,Fast,NULL, NULL,NULL) );

  if(!finite(Yi)||FError)  {  return -1;}
  if(Xf_) *Xf_=Xf; 

  
  return  2.742E8*Mcdm/(1/Yf +  1/Yi); /* 2.828-old 2.755-new 2.742 next-new */
}



static double *Ytab=NULL;

double YF(double T){ return polint2Exp(T,Ntab,Ttab, Ytab) ;}


double darkOmega(double * Xf, int Fast, double Beps)
{
  double Yt,Yi,Xt=25;
  double Z1=1.1;
  double Z2=10; 
  int i;
  double Xf1;
  int Nt=25;
  
  Ytab=realloc(Ytab,sizeof(double)*Nt);
  Ttab=realloc(Ttab,sizeof(double)*Nt);
  Ntab=0;

  if(CDM1==NULL) fracCDM2=1; else
  if(CDM2==NULL) fracCDM2=0; else 
  if(Mcdm1<Mcdm2) fracCDM2=0; else fracCDM2=1;

  if(Xf)*Xf=Xt;
  if(assignVal("Q",2*Mcdm)==0) calcMainFunc() ;
  GGscale=2*Mcdm/3;
  if(Beps>=1) Beps=0.999;
  Beps_=Beps; Fast_=Fast;
  
  if(Z1<=1) Z1=1.1;

  Yt=  darkOmega1(&Xt, Z1, (Z1-1)/5,Fast, Beps);

  if(Yt<0||FError) { return -1;}
  Xf1=Xt;
  
  Tstart=Mcdm/Xt;
//printf("Tstart=%e\n",Tstart);

//printf("Yt=%E deltaY=%E\n", Yt,deltaY);
  
  if(Yt<fabs(deltaY)*1.E-15)
  {  
     if(deltaY>0) dmAsymm=1;  else dmAsymm=-1;
     if(Xf) *Xf=Xt;   
     return 2.742E8*Mcdm*deltaY;  
  }   
  
  Ntab=1;
  Ttab[0]=Tstart;
  Ytab[0]=Yt;
  
  for(i=0; ;i++)
  { double X2=vSigmaGrid.xtop*pow(XSTEP,i+1);
    double y,yeq,alpha;
    double s3_t,s3_2;
    yeq=Yeq(Mcdm/Xt);
    alpha=vSigmaGrid.alpha[i];    

    if(Xt>X2*0.999999) continue; 

    
    if(Yt*Yt>=Z2*Z2*( alpha*Yt*yeq+(1-alpha)*yeq*yeq))  break;
    if(Yt<fabs(deltaY*1E-15))
    {  if(Xf) *Xf=Xt;
       Tend=Mcdm/Xt;
       if(deltaY>0) dmAsymm=1; else dmAsymm=-1;   
       return 2.742E8*Mcdm*deltaY;
    }
    
    Tend=Mcdm/X2;                                 
    y=Yt;
    s3_t=polint2(Mcdm/Xt,Tdim,t_,s3_);
    s3_2=polint2(Mcdm/X2,Tdim,t_,s3_); 
//    if(odeint(&y,1 ,Mcdm/Xt , Mcdm/X2 , 1.E-3, (Mcdm/Xt-Mcdm/X2 )/2, &XderivLn)){ printf("problem in solving diff.equation\n"); return -1;}   
    if(odeint(&y,1 ,s3_t , s3_2 , 1.E-3, (s3_2-s3_t)/2, &XderivLn)){ printf("problem in solving diff.equation\n"); return -1;}
    Yt=y;
    Xt=X2;
    if(Ntab>=Nt)
    { Nt+=20;
      Ytab=realloc(Ytab,sizeof(double)*Nt);
      Ttab=realloc(Ttab,sizeof(double)*Nt);
    }      
    Ytab[Ntab]=Yt;
    Ttab[Ntab]=Tend;
    Ntab++;
  }
//printf("Ntab=%d\n",Ntab);
//for(int i=0;i<Ntab;i++) printf("Y(%.2E)=%.2E\n",Ttab[i], Ytab[i]); 
  
  if(Xf) *Xf=0.5*(Xf1+Xt);
  Yi=1/( (Mcdm/Xt)*sqrt(M_PI/45)*MPlank*aRate(Xt,1,Fast,NULL,NULL,NULL));
  if(!finite(Yi)||FError)  return -1;
  if(deltaY==0.)
  { dmAsymm=0;
    return  2.742E8*Mcdm/(1/Yt  +  1/Yi); /* 2.828-old 2.755-new,2.742 -newnew */
  } else
  {  double a,f,z0,Y0;
     a=fabs(deltaY);
//if(deltaY*1E-5 > Yt ) f= Yt*Yt/2/a*exp(-2*a/Yi);
//else  
     if(Yt<a*1.E-5)  f=Yt*Yt/4/a; else f= (sqrt(Yt*Yt+a*a)-a)/(sqrt(Yt*Yt+a*a)+a);
                        
      f *= exp(-2*a/Yi);
      z0=sqrt(f)*2*a/(1-f);
      Y0=sqrt(z0*z0+a*a);
      dmAsymm=deltaY/Y0;     

     return 2.742E8*Mcdm*Y0;
  }   
}

double printChannels(double Xf ,double cut, double Beps, int prcn, FILE * f)
{ int i,nPrc,nform=log10(1/cut)-2;
  double Sum,s;
  
  if(omegaCh) {free(omegaCh); omegaCh=NULL;}

  MassCut=Mcdm*(2-log(Beps)/Xf);
  Sum=aRate(Xf, 1,1,NULL,&omegaCh,&nPrc)*(Mcdm/Xf)*sqrt(M_PI/45)*MPlank/(2.742E8*Mcdm);
  if(FError)     { return -1;}
  if(nform<0)nform=0;
   
  if(f)
  {  int j;
     fprintf(f,"# Channels which contribute to 1/(omega) more than %G%%.\n",100*cut );
     if(prcn) fprintf(f,"# Relative contributions in %% are displayed\n");
        else  fprintf(f,"# Absolute contributions are displayed\n");
     for(i=0,s=0;i<nPrc;i++)  if(fabs(omegaCh[i].weight)>=cut)
     {  s+=omegaCh[i].weight;
        if(prcn)
        { if(cut <0.000001) fprintf(f,"  %.1E%% ",100*omegaCh[i].weight);
          else              fprintf(f,"  %*.*f%% ",nform+3,nform,
                                        100*omegaCh[i].weight);
        } else fprintf(f,"  %.1E ",Sum*omegaCh[i].weight); 
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
  Sum=aRate(Xf, 1,1,NULL,&wPrc,&nPrc)*(Mcdm/Xf)*sqrt(M_PI/45)*MPlank/(2.742E8*Mcdm);
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

/*================= Z4 ==================*/


static double geff1_(double T)
{ 
   double sum=0,t; int l;
   for(l=0;l<NC;l++) 
   { int k=sort[l];
     if(Z4ch(inP[k])==1) 
     { double bsk2; 
       double M=inMass[k];
       t=T/M;
       if(t<0.1) bsk2=K2pol(t)*exp((Mcdm1-M)/T)*sqrt(M_PI*t/2);
        else     bsk2=bessK2(1/t)*exp(Mcdm1/T);
       sum+=inG[k]*M*M*bsk2;
     }
   }      
   return sum;
}

static double geff2_(double T)
{ 
   double sum=0,t; int l;
   for(l=0;l<NC;l++) 
   { int k=sort[l];
     if(Z4ch(inP[k])==2) 
     { double bsk2; 
       double M=inMass[k];
       t=T/M;
       if(t<0.1) bsk2=K2pol(t)*exp(-1/t+Mcdm2/T)*sqrt(M_PI*t/2);
        else     bsk2=bessK2(1/t)*exp(Mcdm2/T);
       sum+=inG[k]*M*M*bsk2;
     }
   }  
   return sum;
}


static double McdmSum;
double Beps=1.E-4;

static double s_integrandT(double  sqrtS )
{  double sv_tot,t,bsk1;
   double ms,md,PcmIn;
   double res;
   
   ms = M1 + M2; 
   if(ms>=sqrtS)  return 0;
   md = M1 - M2;
   PcmIn = sqrt((sqrtS-ms)*(sqrtS+ms)*(sqrtS-md)*(sqrtS+md))/(2*sqrtS);
   sv_tot=sigma(PcmIn);         
   t=T_/sqrtS; 
   if(t<0.1) bsk1=K1pol(t)*exp(-1/t+McdmSum/T_)*sqrt(M_PI*t/2);
   else      bsk1=bessK1(sqrtS/T_)*exp(McdmSum/T_);
      
   res= sqrtS*sqrtS*(PcmIn*PcmIn)*sv_tot*bsk1/T_;
   return res;
}


static double u_integrand( double u)
{  double z,y,sv_tot,w;
   double Xf_1;
   double ms,md,sqrtS,PcmIn,res0;
   
   if(u==0. || u==1.) return 0.;

   z=1-u*u;
   sqrtS=M1+M2-3*T_*log(z);
   
   return s_integrandT(sqrtS )*6*T_*u/z;
   
}



static  int aRate4(double T, 
  double*vs1100_,double*vs2200_,
  double*vs1110_,double*vs1120_,double*vs1210_,double*vs1220_,double*vs2210_,double*vs2220_,
  double*vs1112_,double*vs1122_,double*vs1222_,double*vs1211_,double*vs2211_,double*vs2221_)
{
  int i,l1,l2,N12;
  char* pname[5];
  double X=Mcdm/T;
 
  double Msmall,Mlarge;
  double vs1100=0, vs2200=0, vs1110=0,vs2220=0,vs1120=0,vs1210=0,vs1122=0,vs2211=0,
         vs1112=0, vs1222=0, vs1220=0,vs2210=0,vs2221=0,vs1211=0;
     
  gridStr grid,grid1; 
   
  WIDTH_FOR_OMEGA=1;
  T_=T;

for(N12=0;N12<=2;N12++)
{ 
  if( N12==1 && (!CDM1 || !CDM2) ) break;
  for(l1=0;l1<NC;l1++)
  { int k1=sort[l1];
  for(l2=0;l2<NC;l2++)
  {
    double factor;
    int k2=sort[l2];
    numout * code=NULL;
    CalcHEP_interface * CI;
    double MassCutOut;
    int*inC;
    
    switch(N12)
    { case 0: inC=inC0; break;
      case 1: inC=inC1; break;
      case 2: inC=inC2; break;
    }  

    if(inC[k1*NC+k2]<=0) continue;

    MassCut=0;    
    if(Z4ch(inP[k1])==1) MassCut+=Mcdm1; else  MassCut+=Mcdm2;
    if(Z4ch(inP[k2])==1) MassCut+=Mcdm1; else  MassCut+=Mcdm2;
    McdmSum=MassCut;
    if(Beps>0)    MassCut-=T*log(Beps); else  MassCut+=1.E20;
    
    MassCutOut=MassCut+T*log(10000.); 
    
    if(inMass[k1]+inMass[k2] > MassCut)
    {
      continue;
    }
    switch(N12)
    { case 0:  
         if(code22_0[k1*NC+k2]==NULL) new_code(k1,k2,0);
         code=code22_0[k1*NC+k2];  break;  
      case 1:                   
       if(Mcdm1 > Mcdm2)
       { if(Z4ch(inP[k1])==2 &&  Z4ch(inP[k2])==2) continue;
         if(code22_1[k1*NC+k2]==NULL) new_code(k1,k2,-1);
         code=code22_1[k1*NC+k2];
       } else 
       {  if(Z4ch(inP[k1])==1 && Z4ch(inP[k2])==1) continue;
          if(code22_1[k1*NC+k2]==NULL) new_code(k1,k2,1);
          code=code22_1[k1*NC+k2];
       } break;
      case 2:
        if(Mcdm1>Mcdm2)
        {   if(Z4ch(inP[k1])==1 && Z4ch(inP[k2])==1) { if(code22_2[k1*NC+k2]==NULL) new_code(k1,k2,2); code=code22_2[k1*NC+k2];}
            else  continue;
        } else 
        {   if(Z4ch(inP[k1])==2 && Z4ch(inP[k2])==2) { if(code22_2[k1*NC+k2]==NULL) new_code(k1,k2,2); code=code22_2[k1*NC+k2];}   
            else  continue;
        }
        break; 
    }
    
    if(!code) continue;
    if(!code->init)
    { 
      CalcHEP_interface *cdi=code->interface;
      for(i=1;i<=cdi->nvar;i++) if(code->link[i]) cdi->va[i]=*(code->link[i]);     
      if(  cdi->calcFunc()>0 ) {FError=1; WIDTH_FOR_OMEGA=0;  return -1;}
      code->init=1;
    }
    
    sqme22=code->interface->sqme;
    
    inBuff=0;

    M1=inMass[k1];
    M2=inMass[k2];

    Msmall=M1>M2? M1-Mcdm*(1-sWidth): M2-Mcdm*(1-sWidth);
    Mlarge=M1>M2? M2+Mcdm*(1-sWidth): M1+Mcdm*(1-sWidth);

    u_max=m2u(MassCutOut);

    factor=inC[k1*NC+k2]*inG[k1]*inG[k2];

    CI=code->interface;
    switch(N12)
    { case 0:  AUX=code22Aux0[k1*NC+k2]; break;
      case 1:  AUX=code22Aux1[k1*NC+k2]; break;
      case 2:  AUX=code22Aux0[k1*NC+k2]; break;
    }                      
    for(nsub22=1; nsub22<= CI->nprc;nsub22++)
    { double smin;
      double a=0;

      int z4[4];
      for(i=0;i<4;i++)  pname[i]=CI->pinf(nsub22,i+1,pmass+i,pdg+i);

      if(pmass[0]==0) continue; // for the case of absence of process
      for(i=0;i<4;i++) z4[i]=Z4ch(pname[i]);
      smin=pmass[2]+pmass[3];      
      cc23=NULL;
      
      if(N12==0 && (VZdecay||VWdecay))
      {  int l,l_,nVV;        
         if(!AUX[nsub22].virt )  for(l=2;l<4;l++) if(pdg[l]==21 ||pdg[l]==22) { AUX[nsub22].virt=-1; break;}
         
         if(!AUX[nsub22].virt)
         {  int vd[4]={0,0,0,0};
            int c_a =  (pmass[0]>Mcdm) || (pmass[1]>Mcdm);

            if(c_a){ for(l=2;l<4;l++) if((pdg[l]==23 && VZdecay>1)   || (abs(pdg[l])==24 && VWdecay>1)) vd[l]=1;} 
            else    for(l=2;l<4;l++)
            { 
            
              if((pdg[l]==23 && VZdecay)     || (abs(pdg[l])==24 && VWdecay)) vd[l]=1;
            } 

            for(l=2;l<4;l++) if(vd[l]) break; 
            if(l<4)
            {  l_=5-l; 
               if(vd[l_])
               { nVV=2;
                 if(pmass[l_]>pmass[l]) { l=l_; l_=5-l;}
               } else nVV=1; 
               AUX[nsub22].virt=l;  
               AUX[nsub22].w[l-2]=pWidth(pname[l],NULL);
               if(abs(pdg[l_])>16 && pmass[l_]> 2) AUX[nsub22].w[l_-2]=pWidth(pname[l_],NULL);
               if(AUX[nsub22].w[l_-2] < 0.1) AUX[nsub22].w[l_-2]=0;
            } else  AUX[nsub22].virt=-1;
         }        
         if(AUX[nsub22].virt>0)
         {  l=AUX[nsub22].virt;
            l_=5-l; 
            if(pmass[0]+pmass[1] < smin+ 4*AUX[nsub22].w[l-2]   && pmass[l_]< MassCutOut)
            { 
              if(AUX[nsub22].cc23) cc23=AUX[nsub22].cc23; else
              {  double  brV1,wV1;
                 int i3W;
                 AUX[nsub22].cc23=xVtoxll(2,2,pname,pdg, l, &wV1, &brV1);
                 if(pdg[l]==pdg[l_]) brV1*=2;
                 AUX[nsub22].br=brV1;
                 cc23=AUX[nsub22].cc23; 
                 if(cc23)
                 {  double Pcm0,PcmMax; 
                    *(cc23->interface->BWrange)=10000; 
                    *(cc23->interface->gswidth)=0;
                    for(i3W=2;i3W<5;i3W++) if(strcmp(cc23->interface->pinf(1,i3W+1,NULL,NULL),pname[l_])==0) break; 
                    AUX[nsub22].i3=i3W;        
                    PcmMax=decayPcm(pmass[2]+pmass[3]+10*AUX[nsub22].w[l-2], pmass[0],pmass[1]);
                    if( pmass[0]+pmass[1]>=pmass[l_]) Pcm0=Mcdm*0.01; else 
                    Pcm0=1.01*decayPcm(pmass[l_], pmass[0],pmass[1]);  
                    buildInterpolation(sigma23, Pcm0,PcmMax, 0.02,1E-5,&(AUX[nsub22].nTab), &(AUX[nsub22].pcmTab), &(AUX[nsub22].csTab));
/*
    printf("nTab=%d\n", AUX[nsub22].nTab);
    displayFunc(sigma23,Pcm0,PcmMax,"sigma23");
    loadINTER(AUX[nsub22].nTab, AUX[nsub22].pcmTab,AUX[nsub22].csTab); 
    displayFunc(INTER,Pcm0,PcmMax,"INTER"); 
*/    
                 }
              }    
              if(cc23) smin=pmass[l_];
            } 
         }
      }


      if(smin>=MassCutOut) continue; 
      if(cc23==NULL) 
      {
        if( (pmass[2]>Mlarge && pmass[3]<Msmall)
          ||(pmass[3]>Mlarge && pmass[2]<Msmall))
             { *(CI->twidth)=1; *(CI->gtwidth)=1;}
        else { *(CI->twidth)=0; *(CI->gtwidth)=0;}
      *(CI->gswidth)=0;
      }
                                   
      if(pmass[0]+pmass[1]> smin) smin=pmass[0]+pmass[1];
repeat:
      neg_cs_flag=0;
      T_=T; 
      if(!Fast_) a=simpson(u_integrand,m2u(smin),m2u(MassCutOut) ,eps); else
      {
          int isPole=0;
          char * s;
          int m,w,n;
          double mass,width;

          for(n=1;(s=code->interface->den_info(nsub22,n,&m,&w));n++)
          if(m && w && strcmp(s,"\1\2")==0 )
          { mass=fabs(code->interface->va[m]);
            width=code->interface->va[w];
            if(mass<MassCutOut && mass+8*width > pmass[0]+pmass[1]
                            && mass+8*width > smin)
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
          if(cc23)
          {  double mass,width;
             mass=pmass[2]+pmass[3];
             width= AUX[nsub22].w[0]+ AUX[nsub22].w[1];
             if(mass-width>M1+M2)
             { grid1=makeGrid(mass,width);
               if(isPole) grid=crossGrids(&grid,&grid1); else grid=grid1;
               isPole++;                  
             }                  
          } 
          if(isPole==0)
          {  grid.n=1;
             grid.ul[0]=m2u(smin);
             grid.ur[0]=u_max;
             grid.pow[0]=5;
          }
          for(i=0;i<grid.n;i++)if(m2u(smin)<grid.ur[i])
          {  
             double ul= m2u(smin)<grid.ul[i]? grid.ul[i]: m2u(smin);
             double da,da_;
             char txt[100];
                da=gauss(u_integrand,ul,grid.ur[i],grid.pow[i]);                
             a+=da;             
          }
      }
       
      if(neg_cs_flag && *(CI->gswidth)==0)
      { *(CI->gswidth)=1;
         goto  repeat;
      }   
      

      { double br[2][5]={0}; 
//br[i][0]-SM, br[i][1]- cdm1, br[i][2]-cmd2, bri[i][3]-2*cdm1, br[i][4]-2*cdm2      
        int l1,l2,kk;   
        for(kk=2;kk<4;kk++)
        { 
          br[kk-2][z4[kk]]=1;
          if(pmass[kk]> 2*Mcdm && z4[kk]!=1 && abs(pdg[kk])!=24 && pdg[kk]!=23 ) 
          { txtList LL;
            pWidth(pname[kk],&LL);
            for(;LL;LL=LL->next)
            { double b;
              char proc[40];
              sscanf(LL->txt,"%lf %[^\n]",&b,proc);
              if(z4[kk]==0)
              {  if(strstr(proc,"~~"))     { br[kk-2][0]-=b; br[kk-2][4]+=b;}
                 else if(strstr(proc,"~")) { br[kk-2][0]-=b; br[kk-2][3]+=b;}
              } else  if(strstr(proc,"~~")==0 && strstr(proc,"~")){ br[kk-2][2]-=b; br[kk-2][3]+=b;}

            }  
          }
        }

        for(l1=0;l1<5;l1++) for(l2=0;l2<5;l2++)
        { double b=br[0][l1]*br[1][l2];
          int NN=0;
          if(b) 
          { if(l1<3 && l2<3) NN=z4[0]+3*z4[1]+9*l1+27*l2;
            else if((l1==0 && l2==3) || (l2==0 && l1==3) ) NN=z4[0]+3*z4[1]+9*1+27*1;
            else if((l1==0 && l2==4) || (l2==0 && l1==4) ) NN=z4[0]+3*z4[1]+9*2+27*2;
            else if((l1==1 && l2==3) || (l2==1 && l1==3) ) NN=0;    //  111               
            else if((l1==1 && l2==4) || (l2==1 && l1==4) ) NN=0;    //  122
            else if((l1==2 && l2==3) || (l2==2 && l1==3) ) NN=0;    //  211 
            else if((l1==2 && l2==4) || (l2==2 && l1==4) ) NN=0;    //  222 
            else if((l1==3 && l2==3)                     ) NN=0;    //  1111
            else if((l1==3 && l2==4) || (l2==3 && l1==4) ) NN=0;    //  1122
            else if((l1==4 && l2==4)                     ) NN=0;    //  2222 
              
            a*=factor*b;
            switch(NN)
            { case 1+3:       vs1100+=a; break;
              case 2+6:       vs2200+=a; break;
              case 1+3+9: 
              case 1+3+27:    vs1110+=a; break;
              case 2+6+18:
              case 2+6+54:    vs2220+=a; break;
              case 1+3+18:
              case 1+3+54:    vs1120+=a; break;
              case 1+6+9 :
              case 1+6+27:
              case 2+3+9 :
              case 2+3+27:    vs1210+=a; break;
              case 1+3+18+54: vs1122+=a; break;
              case 2+6+9+27 : vs2211+=a; break;
//====================              
              case 1+3+ 9+54:
              case 1+3+18+27: vs1112+=a; break;
              case 1+6+18+54:
              case 2+3+18+54: vs1222+=a; break;
              case 1+6+18:
              case 2+3+18:
              case 1+6+54:
              case 2+3+54:    vs1220+=a; break;
              case 2+6+9 :
              case 2+6+27:    vs2210+=a; break;
              case 2+6+9+54:
              case 2+6+18+27: vs2221+=a; break;
              case 1+6+9+27:   
              case 2+3+9+27:  vs1211+=a; break;
              
  //          default:    printf("unexpected type of process: %d %d => %d %d\n",z4[0],z4[1],z4[2],z4[3]);
            }   
          }                  
        }
      }
    }
  }
  }
  }
  { double g1=1,g2=1;
    if(CDM1) g1=geff1_(T);
    if(CDM2) g2=geff2_(T);
    
    *vs1100_=vs1100/(g1*g1);
    *vs1120_=vs1120/(g1*g1);
    *vs1210_=vs1210/(g1*g2);
    *vs2200_=vs2200/(g2*g2);
    *vs1110_=vs1110/(g1*g1);
    *vs2220_=vs2220/(g2*g2);
    *vs2210_=vs2210/(g2*g2);
//==================    

    *vs1220_=vs1220/(g1*g2);
    *vs2210_=vs2210/(g2*g2);

    if(CDM1&&CDM2) 
    { double C=exp(-fabs(Mcdm1/T -Mcdm2/T));
      if(Mcdm1>Mcdm2)
      { 
        *vs1122_=vs1122/(g1*g1);
        *vs2211_=C*C*vs1122/(g2*g2);
        *vs1112_=vs1112/(g1*g1);
        *vs1211_=C*vs1112/(g1*g2);    
        *vs1222_=vs1222/(g1*g2);
        *vs2221_=C*vs1222/(g2*g2);
      } else 
      {
        *vs1122_=C*C*vs2211/(g1*g1);
        *vs2211_=vs2211/(g2*g2);
        *vs1112_=C*vs1211/(g1*g1);
        *vs1211_=vs1211/(g1*g2);    
        *vs1222_=C*vs2221/(g1*g2);
        *vs2221_=vs2221/(g2*g2);
      }                  
    } else *vs1122_=*vs2211_=*vs1112_=*vs1211_=*vs1222_=*vs2221_=0;
  }
  WIDTH_FOR_OMEGA=0;
  
  return 0;
}


double Yeq1(double T)
{  double heff,s;
   s=2*M_PI*M_PI*T*T*T*hEff(T)/45;
   return  (T/(2*M_PI*M_PI*s))*geff1_(T)*exp(-Mcdm1/T);
}


double Yeq2(double T)
{  double s;
   s=2*M_PI*M_PI*T*T*T*hEff(T)/45;
   return  (T/(2*M_PI*M_PI*s))*geff2_(T)*exp(-Mcdm2/T); 
}


static double Y1SQ_Y2(double T)
{ 
  double s,X1,X2,g1_,g2_,res;
  X1=Mcdm1/T;
  X2=Mcdm2/T;
  if(X2-2*X1>500) return 0;
  s=2*M_PI*M_PI*T*T*T*hEff(T)/45;
  
  g1_=geff1_(T);
  g2_=geff2_(T);
      
  res = g1_/g2_*exp(X2-2*X1)*T*g1_/(2*M_PI*M_PI*s);
//  if(!finite(res)) { printf("T=%E  g1_=%E g2_=%E X1=%E X2=%E\n", T,g1_,g2_,X1,X2); exit(0);}
  return res;
}                          

static double Y2SQ_Y1(double T)
{ 
  double s,X1,X2,g1_,g2_,res;
  X1=Mcdm1/T;
  X2=Mcdm2/T;
  if(X1-2*X2>500) return 0;
  s=2*M_PI*M_PI*T*T*T*hEff(T)/45;
  
  g1_=geff1_(T);
  g2_=geff2_(T);
      
  res = g1_/g2_*exp(X1-2*X2)*T*g1_/(2*M_PI*M_PI*s);
//  if(!finite(res)) { printf("T=%E  g1_=%E g2_=%E X1=%E X2=%E\n", T,g1_,g2_,X1,X2); exit(0);}
  return res;
}                          




static int DMEQ0(double T, 
  double vs1100,double vs2200,
  double vs1110,double vs1120,double vs1210,double vs1220,double vs2210,double vs2220, 
  double vs1112,double vs1122,double vs1222,double vs1211,double vs2211,double vs2221,
  double *C, double *L, double *Q)
{
  double y1,y2,y1_, y2_,coef;
  int i;
  
  if(CDM1)  y1=Yeq1(T); else y1=0;
  if(CDM2)  y2=Yeq2(T); else y2=0;
  if(CDM1&&CDM2) { y1_=Y1SQ_Y2(T); y2_=Y2SQ_Y1(T);}  else { y1_=0; y2_=0;}

//vs1120*=0.0;

/* 
  if(Mcdm1 >Mcdm2)z1122= vs1122*(y1*y1-y2*y2*y1_y2_Q(T)); else  z1122=-vs2211*(y2*y2-y1*y1*y2_y1_Q(T));
  dYdT[0]=vs1100*(y1*y1-y1_*y1_)+     vs1120*(y1*y1-y2*y1q_y2) + z1122;
  dYdT[1]=vs2200*(y2*y2-y2_*y2_)- 0.5*vs1120*(y1*y1-y2*y1q_y2) - z1122 + 0.5*vs1210*y1*(y2-y2_) ;
*/


L[0]= y1*(2*(vs1100+vs1120+vs1122)+0.5*(vs1110+vs1112))+ 0.5*y2*(vs1220+vs1222)+0.5*y2_*vs2210; //A1_1
if(!isfinite(L[0])) { printf("y1=%e y2=%e y2_=%e\n", y1,y2,y2_);}
L[3]= y2*(2*(vs2200+vs2210+vs2211)+0.5*(vs2220+vs2221))+ 0.5*y1*(vs1210+vs1211)+0.5*y1_*vs1120; //A2_2



L[1]= -0.5*y1*vs1222 - y2*2*vs2211-0.5*y1*vs1211 -y1_*vs1120 -y2*vs2210;                //A1_2
//printf(" y1*y1*vs1122/(y2*y2*vs2211)=%E  y2*2*vs2211=%E  y1_*vs1120=%E\n",y1*y1*vs1122/(y2*y2*vs2211)  ,    y2*2*vs2211,y1_*vs1120);
L[2]= -0.5*y2*vs1211 - y1*2*vs1122-0.5*y2*vs1222 -y2_*vs2210 -y1*vs1120;                //A2_1

//printf("L[2]=%e = %e + %e  vs2210=%E   \n", L[2], - y1_*2*vs1122, -vs2210*y2q_y1, vs2210);
 
Q[0] =  vs1100+vs1120+vs1122+0.5*(vs1110+vs1112); //Q1_11
Q[5] =  vs2200+vs2210+vs2211+0.5*(vs2220+vs2221); //Q2_22

Q[1] =  0.5*(vs1222+vs1220-vs1211);  //Q1_12
Q[4] =  0.5*(vs1211+vs1210-vs1222);  //Q2_12

Q[2] = -vs2211 -0.5*vs2221 -0.5*vs2210;          //Q1_22
Q[3] = -vs1122 -0.5*vs1112 -0.5*vs1120;          //Q2_11

#ifdef OLD 
 
  if(Mcdm1 >Mcdm2){
                    Q[0]=vs1122;       Q[1]=0;     Q[2]=-vs1122*y1q_y2;
                    L[0]=2*vs1122*y1_; L[1]=-2*vs1122*y2_*y1q_y2;    
                  } else  
                  {
                    Q[0]=vs2211*y2q_y1; Q[1]=0; Q[2]=-vs2211;
                    L[0]=2*vs2211*y1_*y2q_y1; L[1]=-2*vs2211*y2_; 
                  }
                  
for(i=0;i<3;i++) Q[3+i]=-Q[i];
for(i=0;i<2;i++) L[2+i]=-L[i];                    
                  
    L[0]+=2*vs1100*y1_ +2*vs1120*y1_;

if(vs1120>0) L[1]+=-vs1120*y1q_y2;

if(!finite(L[1])){ printf("stop2 vs1120=%E y1q_y2=%E \n",vs1120,y1q_y2); exit(1);}
    Q[0]+=vs1100+vs1120;                                         
    L[2+0]+=-vs1120*y1_;
    
    L[2+1]+=2*vs2200*y2_ +0.5*vs1210*y1_;
    if(vs1120>0) L[2+1]+=0.5*vs1120*y1q_y2;

    Q[3+0]+=-0.5*vs1120;
    Q[3+1]+=0.5*vs1210;
    Q[3+2]+=vs2200;                                            

#endif 

  { double dT=T/100;
    double heff_=0.5*(hEff(T+dT)-hEff(T-dT))/dT;
    double sqrt_gStar=(hEff(T)+ T/3*heff_)/sqrt(gEff(T));
    coef=sqrt(M_PI/45)*MPlank*sqrt_gStar;
  }  
  for(i=0;i<4;i++) L[i]*=coef;
  for(i=0;i<6;i++) Q[i]*=coef;
  
  
  if(CDM1)  C[0]=(Yeq1(T*1.01)-Yeq1(T/1.01))/(2*log(1.01)*T); else C[0]=0;
  if(CDM2)  C[1]=(Yeq2(T*1.01)-Yeq2(T/1.01))/(2*log(1.01)*T); else C[1]=0;

{ int i;
for(i=0;i<2;i++) if(!isfinite(C[i])) {printf("T=%E, C[%d]=%E\n",T,i,C[i]); exit(1);}
for(i=0;i<4;i++) if(!isfinite(L[i])) {printf("T=%E, L[%d]=%E\n",T,i,L[i]); exit(1);}
for(i=0;i<6;i++) if(!isfinite(Q[i])) {printf("T=%E, Q[%d]=%E\n",T,i,Q[i]); exit(1);}
}
  return 0;
}

static int DMEQ(double T,double *C, double *L, double *Q)
{
  double  vs1100,vs2200,vs1110,vs2220,vs1120,vs1122,vs1210,vs2211,vs1112,vs1222,vs1220,vs2210,vs2221,vs1211;
        aRate4(T,&vs1100,&vs2200,
    &vs1110,&vs1120,&vs1210,&vs1220,&vs2210,&vs2220,
    &vs1112,&vs1122,&vs1222,&vs1211,&vs2211,&vs2221
        );
  return  DMEQ0(T, vs1100,vs2200,
    vs1110,vs1120,vs1210,vs1220,vs2210,vs2220,
    vs1112,vs1122,vs1222,vs1211,vs2211,vs2221, 
         C,L,Q);
}

static int dYstart(double T, double * dy, double * Lmin,double *Lmax)
{
  double C[2],L[4],Q[6],D;
  DMEQ(T, C, L,Q);

  if(!CDM2) { dy[0]=C[0]/L[0]; dy[1]=0; return 0;}
  if(!CDM1) { dy[0]=0; dy[1]=C[1]/L[3]; return 0;}
    
  D=L[0]*L[3]-L[1]*L[2];
  if(dy)
  { 
    dy[0]= ( C[0]*L[3]-C[1]*L[1])/D;
    dy[1]= (-C[0]*L[2]+C[1]*L[0])/D;
  }  
{ 
  double s,d;
  s=0.5*(L[0]+L[3]);
  d= s*s+L[1]*L[2]-L[0]*L[3];
  if(Lmin) *Lmin=s-sqrt(d);
  if(Lmax) *Lmax=s+sqrt(d);
}  
  return 0;
}




static void TderivZ4(double T, double *Y, double *dYdT)
{
  double C[2], L[4], Q[6],dy1,dy2;
   
  DMEQ(T, C, L, Q);
  
  dy1=Y[0];
  dy2=Y[1];
  
  dYdT[0]= -C[0] + L[0]*dy1 + L[1]*dy2 + Q[0]*dy1*dy1 + Q[1]*dy1*dy2 + Q[2]*dy2*dy2;
  dYdT[1]= -C[1] + L[2]*dy1 + L[3]*dy2 + Q[3]*dy1*dy1 + Q[4]*dy1*dy2 + Q[5]*dy2*dy2;
}


static void stiffDerives(double T, double*Y,double*f,double h,double*dfdx,double*dfdy)
{

 double C[2], L[4], Q[6],dy1,dy2;
 double dT=-0.001*T;
 int n=2;   
 DMEQ(T, C, L, Q);
        
 dy1=Y[0];
 dy2=Y[1];
              
  f[0]= -C[0] + L[0]*dy1 + L[1]*dy2 + Q[0]*dy1*dy1 + Q[1]*dy1*dy2 + Q[2]*dy2*dy2;
  f[1]= -C[1] + L[2]*dy1 + L[3]*dy2 + Q[3]*dy1*dy1 + Q[4]*dy1*dy2 + Q[5]*dy2*dy2;
 
  if(dfdy)
  { dfdy[0*n+0]= L[0]+2*Q[0]*dy1+Q[1]*dy2;               
    dfdy[0*n+1]= L[1]+2*Q[2]*dy2+Q[1]*dy1;
    dfdy[1*n+0]= L[2]+2*Q[3]*dy1+Q[4]*dy2;
    dfdy[1*n+1]= L[3]+2*Q[5]*dy2+Q[4]*dy1;
  }

  if(dfdx)
  { DMEQ(T+dT, C, L, Q);
    dfdx[0]= -C[0] + L[0]*dy1 + L[1]*dy2 + Q[0]*dy1*dy1 + Q[1]*dy1*dy2 + Q[2]*dy2*dy2;
    dfdx[1]= -C[1] + L[2]*dy1 + L[3]*dy2 + Q[3]*dy1*dy1 + Q[4]*dy1*dy2 + Q[5]*dy2*dy2;
    dfdx[0]-=f[0]; dfdx[0]/=dT;
    dfdx[1]-=f[1]; dfdx[1]/=dT;
  }
}

static void TabDmEq(double step, int show)
{
  int i,N;
  double T;

//printf("TabDmEq:  Tstart=%E Tend=%E\n",Tstart,Tend);  
  N=log(Tstart/Tend)/log(step)+2;
  if(N!=Ntab)
  {                 
    Ttab    = realloc(Ttab,   N*sizeof(double));
    vs1100T = realloc(vs1100T,N*sizeof(double));
    vs1120T = realloc(vs1120T,N*sizeof(double));
    vs1122T = realloc(vs1122T,N*sizeof(double));
    vs1210T = realloc(vs1210T,N*sizeof(double));
    vs2200T = realloc(vs2200T,N*sizeof(double));
    vs1110T = realloc(vs1110T,N*sizeof(double));
    vs2220T = realloc(vs2220T,N*sizeof(double));
    vs2211T = realloc(vs2211T,N*sizeof(double));
    vs1112T = realloc(vs1112T,N*sizeof(double)); 
    vs1222T = realloc(vs1222T,N*sizeof(double));
    vs1220T = realloc(vs1220T,N*sizeof(double));
    vs2210T = realloc(vs2210T,N*sizeof(double));
    vs2221T = realloc(vs2221T,N*sizeof(double));
    vs1211T = realloc(vs1211T,N*sizeof(double));
    Y1T     = realloc(Y1T,    N*sizeof(double));
    Y2T     = realloc(Y2T,    N*sizeof(double)); 
    Ntab=N;
  }  
  for(T=Tstart,i=0;i<N;i++)
  { 
     Ttab[i]=T;
     Y1T[i]=0;
     Y2T[i]=0;
     aRate4(T, vs1100T+i,vs2200T+i,
        vs1110T+i,vs1120T+i,vs1210T+i,vs1220T+i,vs2210T+i,vs2220T+i,
        vs1112T+i,vs1122T+i,vs1222T+i,vs1211T+i,vs2211T+i,vs2221T+i);
     T/=step;
//printf("vs1100T[%d]=%e\n",i,vs1100T[i]);          
  }
  if(show)
  {
    if(CDM1 && CDM2)
    {  displayFunc(vs1120F,Tend,Tstart,"vs1120");
       displayFunc(vs1210F,Tend,Tstart,"vs1210");
       if(Mcdm1<Mcdm2) displayFunc(vs2211F,Tend,Tstart,"vs2211");
       else            displayFunc(vs1122F,Tend,Tstart,"vs1122");
    }
    if(CDM2) displayFunc(vs2200F,Tend,Tstart,"vs2200");
    if(CDM1) displayFunc(vs1100F,Tend,Tstart,"vs1100");
  }  
}

static int DMEQtab(double T,double *C, double *L, double *Q)
{
  double vs1100,vs2200,vs1110,vs2220,vs1120,vs1122,vs1210,vs2211,vs1112,vs1222,vs1220,vs2210,vs2221,vs1211;
  vs1100 = polint2Exp(T,Ntab,Ttab, vs1100T);  if(vs1100<0) vs1100 = 0;
  vs1120 = polint2Exp(T,Ntab,Ttab, vs1120T);  if(vs1120<0) vs1120 = 0;
  vs1122 = polint2Exp(T,Ntab,Ttab, vs1122T);  if(vs1122<0) vs1122 = 0;
  vs1210 = polint2Exp(T,Ntab,Ttab, vs1210T);  if(vs1210<0) vs1210 = 0;
  vs2200 = polint2Exp(T,Ntab,Ttab, vs2200T);  if(vs2200<0) vs2200 = 0;
  vs2211 = polint2Exp(T,Ntab,Ttab, vs2211T);  if(vs2211<0) vs2211 = 0;
  vs1110 = polint2Exp(T,Ntab,Ttab, vs1110T);  if(vs1110<0) vs1110 = 0;
  vs2220 = polint2Exp(T,Ntab,Ttab, vs2220T);  if(vs2220<0) vs2220 = 0;
  vs1112 = polint2Exp(T,Ntab,Ttab, vs1112T);  if(vs1112<0) vs1112 = 0; 
  vs1222 = polint2Exp(T,Ntab,Ttab, vs1222T);  if(vs1222<0) vs1222 = 0;
  vs1220 = polint2Exp(T,Ntab,Ttab, vs1220T);  if(vs1220<0) vs1220 = 0;
  vs2210 = polint2Exp(T,Ntab,Ttab, vs2210T);  if(vs2210<0) vs2210 = 0;
  vs2221 = polint2Exp(T,Ntab,Ttab, vs2221T);  if(vs2221<0) vs2221 = 0;
  vs1211 = polint2Exp(T,Ntab,Ttab, vs1211T);  if(vs1211<0) vs1211 = 0;
  
// printf("vs1100=%E,vs1120=%E, vs1122=%E, vs1210=%E, vs2200=%E, vs2211=%E\n",vs1100,vs1120, vs1122, vs1210, vs2200, vs2211);
  return   DMEQ0(T,vs1100,vs2200,
    vs1110,vs1120,vs1210,vs1220,vs2210,vs2220, 
    vs1112,vs1122,vs1222, vs1211,vs2211,vs2221, 
    C,L,Q);
}


static void TderivZ4tab(double T, double *Y, double *dYdT)
{
  double C[2], L[4], Q[6],dy1,dy2;

  DMEQtab(T,C,L,Q);
  
  dy1=Y[0];
  dy2=Y[1];
// printf("T=%E C= %E %E L= %E %E %E %E  Q=%E %E %E %E %E %E \n",T,C[0],C[1],L[0],L[1],L[2],L[3],Q[0],Q[1],Q[2],Q[3],Q[4],Q[5] );      
  dYdT[0]= -C[0] + L[0]*dy1 + L[1]*dy2 + Q[0]*dy1*dy1 + Q[1]*dy1*dy2 + Q[2]*dy2*dy2;
  dYdT[1]= -C[1] + L[2]*dy1 + L[3]*dy2 + Q[3]*dy1*dy1 + Q[4]*dy1*dy2 + Q[5]*dy2*dy2;
}


static void TderivZ4tab2(double T, double *Y, double *dYdT)
{
  double C[2], L[4], Q[6],dy1,dy2;

  DMEQtab(T,C,L,Q);  
  dy2=Y[0];

if(!CDM1)dy1=0;else  dy1=( C[0] - L[1]*dy2 - Q[2]*dy2*dy2 )/(L[0]   + Q[1]*dy2);
   
  dYdT[0]= -C[1] + L[2]*dy1 + L[3]*dy2 + Q[3]*dy1*dy1 + Q[4]*dy1*dy2 + Q[5]*dy2*dy2;
}

static void TderivZ4tab1(double T, double *Y, double *dYdT)
{
  double C[2], L[4], Q[6],dy1,dy2;

  DMEQtab(T,C,L,Q);
  
  dy1=Y[0];

if(!CDM2) dy2=0; else { 
double D;
dy2= (C[1] - L[2]*dy1  - Q[3]*dy1*dy1)/(L[3]+Q[4]*dy1); 
D=(L[3]+Q[4]*dy1)/Q[5]/2;

dy2= -D + sqrt( (C[1]-L[2]*dy1-Q[3]*dy1*dy1)/Q[5] +D*D);

}
  dYdT[0]= -C[0] + L[0]*dy1 + L[1]*dy2 + Q[0]*dy1*dy1 + Q[1]*dy1*dy2 + Q[2]*dy2*dy2;
}




double darkOmega2( double fast, double Beps0)
{
  Fast_=fast;
  Beps=Beps0;
  int info=0;
  double Y[2],YY[2],T;
  double Lmin,Lmax;
  double step=1.1;
  double ips=0.01,ips_=0.005;  
  int i,err,N; 
  
  Tend=1.E-3;

  dmAsymm=0;
  if(!CDM1) Tstart= Mcdm2/20; else if(!CDM2) Tstart= Mcdm1/20;
  else { if(Mcdm1>Mcdm2) Tstart=Mcdm1/20; else Tstart= Mcdm2/20;} 
  
  dYstart(Tstart,Y,&Lmin,&Lmax);
  
  if(!CDM1)
  {     
     while(fabs(Y[1])>0.01 *Yeq2(Tstart)) { Tstart*=1.05; dYstart(Tstart,Y,&Lmin,&Lmax);}
     while(fabs(Y[1])<0.005*Yeq2(Tstart)) { Tstart/=1.05; dYstart(Tstart,Y,&Lmin,&Lmax);}
  } else if(!CDM2)
  {
      while(fabs(Y[0])>0.01 *Yeq1(Tstart)) { Tstart*=1.05; dYstart(Tstart,Y,&Lmin,&Lmax);}
      while(fabs(Y[0])<0.005*Yeq1(Tstart)) { Tstart/=1.05; dYstart(Tstart,Y,&Lmin,&Lmax);} 
  } else
  {           

     while( Lmin<100/Tstart  ) { Tstart*=1.05; dYstart(Tstart,Y,&Lmin,&Lmax);} 
     while( Lmin>200/Tstart  ) { Tstart/=1.05; dYstart(Tstart,Y,&Lmin,&Lmax);} 

//     while(fabs(Y[0])>ips*Yeq1(Tstart) || fabs(Y[1])>ips*Yeq2(Tstart)) { Tstart*=1.05; dYstart(Tstart,Y,&Lmin,&Lmax);} 
//     while(fabs(Y[0])<ips_*Yeq1(Tstart)&&fabs(Y[1])<ips_*Yeq2(Tstart)) { Tstart/=1.05; dYstart(Tstart,Y,&Lmin,&Lmax);} 
  }
  
//  if(info) 
  printf("Tstart=%.2E  dY1=%.2E(%.2E) dY2=%.2E(%.2E)  Lmin=%.2E Lmax=%.2E \n",Tstart, Y[0], Y[0]/Yeq1(Tstart),Y[1],Y[1]/Yeq2(Tstart),Lmin,Lmax);  

//Tstart=Tstart*pow(step,10);
//dYstart(Tstart,Y,&Lmin,&Lmax);

  TabDmEq(step,0);

//Y[0]*=5;
//Y[1]*=5;

  

  Y1T[0]=Y[0];
  Y1T[1]=Y[1]; 
  
 
  if(!CDM1)
  {
    for(i=1,err=0,T=Tstart; T> Tend && !err;i++ )
    { double T2=T/step;
      if(T2<Tend) T2=Tend;
      if(info) printf(" CMD2 : T1=%.2E Y=%.2E ",T,Y[1]);
      err=odeint(Y+1,1 , T ,T2 , 1.E-3, (T-T2) , &TderivZ4tab2);
      if(err) { printf(" error in odeint\n");  return  -1;}
      if(info)  printf(" CMD2 : T2=%.2E Y=%.2E\n",T2,Y[1]);
      Y1T[i]=0;      
      Y2T[i]=Y[1];
      T=T2;
    }
    fracCDM2=1;
    return Y[1]*2.742E8*Mcdm2; 
  } else if(!CDM2) 
  { 
  
    for(i=1,err=0,T=Tstart; T> Tend && !err;i++ )
    { double T2=T/step;
      if(T2<Tend) T2=Tend;
      if(info) printf(" CMD1 : T1=%.2E Y=%.2E ",T,Y[0]);
      err=odeint(Y,1 , T ,T2 , 1.E-3, (T-T2) , &TderivZ4tab1);
      if(err) { printf(" error in odeint\n");  return  -1;}
      if(info)  printf(" CMD1 : T2=%.2E Y=%.2E\n",T2,Y[0]);
      Y1T[i]=Y[0];      
      Y2T[i]=0;
      T=T2;
    }
    fracCDM2=0;
    return Y[0]*2.742E8*Mcdm1;
  }  else  
  { double h=0.01*Tstart*(1-1/step);
    for(i=1,err=0,T=Tstart; T> Tend && !err;i++ )
    { double T2=T/step;
      double Yscal[2]={1,1};
      if(T2<Tend) T2=Tend;
      
      
//      if(fast) err= RungeKuta2minus(Y,T,T2,info); else err=odeint(Y,2 , T ,T2 , 1.E-3, (T-T2) , TderivZ4tab);    
      err=stiff(i==1,T,T2,2,Y, Yscal,1.E-3, &h, stiffDerives);
//    err=stifbs(i==1,T, T2, 2, Y, Yscal, 1.E-3, &h,stiffDerives);
             
      if(err) { printf(" error in stiff\n");  return  -1;}
      Y1T[i]=Y[0];      
      Y2T[i]=Y[1];      
      T=T2;
      
/*
      if(i>15) 
      { int k;
        for(k=1;k<10;k++)
        { 
          if(fabs(Y1T[i-k]-Y1T[i]*Y1T[i-10]*(Ttab[i-10]-Ttab[i])/(Y1T[i-10]*(Ttab[i-10]-Ttab[i-k])-Y1T[i]*(Ttab[i]-Ttab[i-k])))> 0.001*fabs(Y[0]) ) break;
          if(fabs(Y2T[i-k]-Y2T[i]*Y2T[i-10]*(Ttab[i-10]-Ttab[i])/(Y2T[i-10]*(Ttab[i-10]-Ttab[i-k])-Y2T[i]*(Ttab[i]-Ttab[i-k])))> 0.001*fabs(Y[1]) ) break;
        }
        if(k==10) { Tend=T;
                    printf("i=%d Tstart=%E Tend=%E\n",i, Tstart,Tend);
        Y[0]=Y1T[i]*Y1T[i-10]*(Ttab[i-10]-Ttab[i])/(Y1T[i-10]*(Ttab[i-10])-Y1T[i]*(Ttab[i]));
        Y[1]=Y2T[i]*Y2T[i-10]*(Ttab[i-10]-Ttab[i])/(Y2T[i-10]*(Ttab[i-10])-Y2T[i]*(Ttab[i]));
          break;}
      }
*/      
        
//printf("T=%.5E h=%.2E Y={%.3E %.3E}\n",T2,h,Y[0],Y[1]);
    }
//i--; 
//   Y[0]=Y1T[i]*Y1T[i-10]*(Ttab[i-10]-Ttab[i])/(Y1T[i-10]*(Ttab[i-10])-Y1T[i]*(Ttab[i]));
//   Y[1]=Y2T[i]*Y2T[i-10]*(Ttab[i-10]-Ttab[i])/(Y2T[i-10]*(Ttab[i-10])-Y2T[i]*(Ttab[i]));
                 
    fracCDM2=Y[1]*Mcdm2/( Y[0]*Mcdm1 +Y[1]*Mcdm2);

    return  Y[0]*2.742E8*Mcdm1+ Y[1]*2.742E8*Mcdm2;
  }
}

