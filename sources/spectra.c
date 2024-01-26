#include"micromegas.h"
#include"micromegas_aux.h"
#include"micromegas_f.h"

int SpectraFlag=0;

#define NEn 36  /* Number of Energy points for interpolation */
#define Nin 21 
#define Nout 6
#define QCUT 25 //2.5
#define X1CUT (2.E-2)

#define V0 (vRot/299792.*1.732)  /*  1.732=sqrt(3)  */

//#define V0 1E-2
 
//#define DISPLAY_SPECTRA 
/*#define ECTEST*/

static double Mcdm0;

static char* errMess=NULL;



aChannel* vSigmaCh=NULL;
static int nAnCh=0;
static double GG=1.23;


static float phidiff[Nin][Nout][NEn][NZ];


double Zi(int i) { return log(1.E-7)*pow((double)(i-1)/(double)(NZ),1.5) ;}
static  int Iz(double z) { return NZ*pow(z/log(1.E-7),1./1.5)+1; }

void  fillSpect(double (*dNdE)(double ), double Emax, double * SpectAr)
{ int i;
  SpectAr[0]=Emax;
  for(i=1;i<NZ;i++) 
  { double E=exp(Zi(i))*Emax;  
    SpectAr[i]=dNdE(E)*E;
  }
}



// copied from CalcHEP/c_source/dynamicME/kin4.c
static double kinematic_1_3(REAL *pmass, int i3, double m12, double xcos, REAL * P)
{ 
  double factor;
  REAL pout,chY,shY,xsin, E1,P12,P13,E2,P22,P23, m0,m1,m2,m3;
  int i,i1,i2;
  
  for(i=1;i<4;i++)if(i3!=i) {i1=i; break;}
  for(i++;i<4;i++)if(i3!=i) {i2=i; break;}
  
  m0=pmass[0];
  m1=pmass[i1];
  m2=pmass[i2];
  m3=pmass[i3];

  if(m12<=m1+m2) return 0;
  for(i=0;i<16;i++) P[i]=0;

  P[0]=m0; 
  factor=1/(64*M_PI*M_PI*M_PI*m0*m0);
  
  pout=decayPcm(m0,m12,m3);
  if(!pout) return 0;  
  P[i3*4]=Sqrt(pout*pout+m3*m3); P[i3*4+3]=-pout; 

  factor*=pout;  
  
  shY=pout/m12;
  chY=Sqrt(1+shY*shY);  
  pout=decayPcm(m12,m1,m2);
  if(!pout) return 0;
  factor*=pout;
  xsin=Sqrt(1-xcos*xcos);
  E1=Sqrt(m1*m1+pout*pout);    E2=Sqrt(m2*m2+pout*pout);
  P13=xcos*pout;               P23=-P13;
  P12=xsin*pout;               P22=-P12;
  
  P[4*i1]  =chY*E1 + shY*P13;  P[4*i2]  =chY*E2 + shY*P23;
  P[4*i1+3]=shY*E1 + chY*P13;  P[4*i2+3]=shY*E2 + chY*P23;
  P[4*i1+2]=P12;               P[4*i2+2]=P22;
  
  return factor;
}

typedef struct
{  
   numout*cc;
   REAL mass[4];
   int i3;
   int nsub;
   double m12;
   double GG;
   int err;
}   argFor13Int;


static double dWidth13dCos(double xcos, void*arg )
{  REAL pvect[16];
   REAL mass[4];
   argFor13Int*arg13=arg;
   
   double factor=kinematic_1_3(arg13->mass,arg13->i3,arg13->m12,xcos, pvect);
   if(factor==0) return 0;
   return   factor*(arg13->cc->interface->sqme)(arg13->nsub,arg13->GG, pvect,NULL,&(arg13->err));
}

static double dWidth13dM(double m12, void*arg )
{
   ((argFor13Int*)arg)->m12=m12;
   int err;
   return simpson_arg(dWidth13dCos, arg, -1, 1, 1E-2,&err);
}



/* v*cs22 at v=0 */  
double  vcs22(numout * cc,int nsub,int * err)
{
   int i;
   double pcm,r;
   REAL pmass[4], pvect[16];
   
//printf("Mb4a=%E\n", findValW("Mb"));    
   for(i=1;i<=cc->interface->nvar;i++) if(cc->link[i]) cc->interface->va[i]=*(cc->link[i]);
//printf("Mbb=%E Q=%E\n", findValW("Mb"),  findValW("Q")); 
   if( cc->interface->calcFunc()>0 ) {*err=4;  return 0;}
//printf("Mbb_=%E Q=%E\n", findValW("Mb"),  findValW("Q"));
   *(cc->interface->gtwidth)=0;
   *(cc->interface->twidth)=0;
   *(cc->interface->gswidth)=0;
//printf("Mb4c=%E\n", findValW("Mb"));    
   for(i=0;i<4;i++)  cc->interface->pinf(nsub,1+i,pmass+i,NULL);   
   *err=0;
   if(pmass[0]+pmass[1] <= pmass[2]+pmass[3]) return 0;
   for(i=0;i<16;i++) pvect[i]=0;
//printf("Mb4d=%E\n", findValW("Mb")); 
   pcm= decayPcm(pmass[0]+pmass[1],pmass[2],pmass[3]);
   for(i=0;i<2; i++) pvect[4*i]=pmass[i];
   for(i=2;i<4; i++) pvect[4*i]=Sqrt(pmass[i]*pmass[i] +pcm*pcm);
   pvect[8+3]=pcm;
   pvect[12+3]=-pcm;
   r=cc->interface->sqme(nsub,GG,pvect,NULL,err);
   return 3.8937966E8*r*pcm/(16*M_PI*pmass[0]*pmass[1]*(pmass[0]+pmass[1]));  
}
 
/* New 2->3 */ 
static REAL pmass[5],pvect[20];
static int code[5];
static int iA,ix,iX;
static numout* cc23;
static double X0=0.01, dSigmadE_x0, Egamma,X1=0.8 ,dSigmadE_x1,dSigmadE_x1_e;
static int PrintOn=0;
static double dSigmadCos23(double csfi)
{
  REAL pcm1,pcm2,ms,md,chY,shY,M;
  double  r;
  int i,err_code;

  for(i=0;i<20;i++) pvect[i]=0;

  pvect[0]=pmass[0];
  pvect[4]=pmass[1];

  pvect[4*iA]= Egamma;
  pvect[4*iA+3] = -Egamma;  

  ms=pmass[ix]+pmass[iX];
  md=pmass[iX]-pmass[ix];
  pcm1=Egamma;

  M=4*pmass[0]*(pmass[0]-Egamma);
  if(M<=ms*ms) return 0;
  
  M= Sqrt(M);
  
  pcm2=Sqrt((M*M-ms*ms)*(M*M-md*md))/M/2;
  
  pvect[4*ix]=Sqrt(pmass[ix]*pmass[ix]+pcm2*pcm2);
  pvect[4*ix+3]=pcm2*csfi;
  pvect[4*ix+2]=pcm2*Sqrt(1-csfi*csfi);
  
  pvect[4*iX]= Sqrt(pmass[iX]*pmass[iX]+pcm2*pcm2);
  pvect[4*iX+3]=-pvect[4*ix+3];
  pvect[4*iX+2]=-pvect[4*ix+2];  

  chY=Sqrt(1+pcm1*pcm1/M/M);
  shY= pcm1/M;  
  
  { double p0,p3;
    p0=pvect[4*ix], p3=pvect[4*ix+3];
    pvect[4*ix]=  chY*p0 + shY*p3;
    pvect[4*ix+3]=shY*p0 + chY*p3;

    p0=pvect[4*iX]; p3=pvect[4*iX+3];
    pvect[4*iX]=  chY*p0 + shY*p3;
    pvect[4*iX+3]=shY*p0 + chY*p3;
  }
  
/*
  for(i=0;i<4;i++)
  { int j; 
    double s=0;
    s=pvect[i]+pvect[4+i]-pvect[8+i]-pvect[12+i]-pvect[16+i] ;
    printf("s(%d)=%Egamma\n",i,s);
  }              
*/  
err_code=0;
  r=(cc23->interface->sqme)(1,GG,pvect,NULL,&err_code)*pcm1*pcm2/(128*M*pmass[0]*pmass[0]*M_PI*M_PI*M_PI);
if(r<0) r=0;
  if(err_code) return 0;
//printf("csfi=%E r=%e\n", csfi,r);  
  return r*3.8937966E8;
}

static double dSigmadFi23(double fi)
{if(fi==0 || fi==M_PI) return 0; return dSigmadCos23(cos(fi))*sin(fi);}

static double dSigmadE(double E)
{ double r;
   Egamma=E;
/*
   if(pmass[ix]>1.E-3*Mcdm0 && pmass[iX]>1.E-3*Mcdm0  ) 
                       r= simpson(dSigmadFi23,0,M_PI,1.E-4,NULL);
   else { printf("dSigmadCos23\n"); //displayFunc("dSigmadCos23",dSigmadCos23, -1 , 1,0,"dSigmadCos23");
                 r= simpson(dSigmadCos23,-0.9,0.9,1.E-4,NULL); printf("ok\n");  
*/
//                   r= gauss(dSigmadCos23,-1,1,7); }
   r= gauss(dSigmadCos23,-1,1,7);                 
//   r= simpson(dSigmadFi23,0,M_PI,1.E-4,NULL);
   return r;
}

static int addGamma(int pdg)
{ pdg=abs(pdg); 
  if(pdg<=6) return 0;
  if(pdg>=11 && pdg<=15) return 0;
  if(pdg==81||pdg==83) return 0;
  return 1;
}  

static double dSigmadERest(double E)
{
  double res, eps,ms= (pmass[ix]+pmass[iX]);
  int l,nn[2]={ix,iX};
  
  if(E+sqrt(E*E + ms*ms)>=1.999999*Mcdm0) return 0;
  res= dSigmadE(E); if(res==0.) return 0;
  eps=ms/2/Mcdm0;

//  if( !(addGamma(code[iX]) && addGamma(code[ix])))
  { double subtract=0,norm=0,x=E/Mcdm0;
    double csmax,pcm,pcm0;
    
     if(pmass[ix]>1.E-3*Mcdm0 && pmass[iX]>1.E-3*Mcdm0 ) csmax=1; else csmax=0.999;
     pcm=decayPcm(2*Mcdm0*Sqrt(1-x),pmass[ix],pmass[iX]);
     for(l=0;l<2;l++)  
     {  
       double kappa;
       if(!addGamma(code[nn[l]]))
       { kappa=sqrt(1/(1+pow(pmass[nn[l]]/pcm,2)));
         subtract += log((1+kappa*csmax)/(1-kappa*csmax))*(1-x+x*x/2 -eps*eps/2)-kappa*(1-x)*csmax;
       }  
       pcm0=decayPcm(2*Mcdm0,pmass[ix],pmass[iX]);
       kappa=sqrt(1/(1+pow(pmass[nn[l]]/pcm0,2)));
       norm+=        log((1+csmax*kappa)/(1-csmax*kappa))*(1 -eps*eps/2)-kappa*csmax;
     }
     res-=subtract/norm*dSigmadE_x0/x; 
  }
  if(res<0) return 0;
  
  return res;
}

static double Xe;
static double FE(double cs)
{ 
  Egamma=2*Mcdm0*(1-Xe)/(1+cs);
  if(Egamma<X0*Mcdm0) return 0;
  if(Egamma/Mcdm0 -(1-Xe) < -1.E-4) return 0;
  return dSigmadCos23(cs)*2/(1+cs);
}

static double FEfi(double fi)
{ 
  if(fi< 1.E-4  || fi > M_PI - 1.E-4) return 0;
  return FE(cos(fi))*sin(fi);
}

static double dSigmadEe(double E)
{
   double xe,x1=0.96,x2=0.98;
   xe=Xe=E/Mcdm0; 
   if(xe> x2+1.E-4)
   { 
     return   ((x2-xe)*dSigmadEe(x1*Mcdm0) + (xe-x1)*dSigmadEe(x2*Mcdm0))/(x2-x1);
   }           
   if(Xe>X0 )
   { double csMin=(2*(1-Xe) -1)*1.0001;
     double csMax=(2*(1-Xe)/X0 -1);
     double r;
/*     csMax=1;     

        if(csMin<-0.98) csMin=-0.98;
        if(csMax> 0.98) csMax=0.98;
        return simpson(FE,csMin,csMax,1.E-3,NULL);

        if(csMin<-1) csMin=-1;
        if(csMax> 1) csMax=1;
        return gauss(FE,csMin,csMax,7);     
*/
     if(csMin<-1) csMin=-1;
     if(csMax> 1) csMax=1;
     csMin+=0.001; // !!!!!!!! Problem in  case of zero mass of electron. Reason is not clear 
int err;
     r= simpson(FEfi,acos(csMax),acos(csMin),1.E-4,NULL);
//if(err) displayPlot("dSigmadEe", "fi",acos(csMax),acos(csMin),0,1,"FEfi",0,FEfi,NULL);
     return r;
   } else return 0;
}
#ifdef DISPLAY_SPECTRA
static double EdSigmadEe(double E) {return E*dSigmadEe(E);}
static double EdSigmadERest(double E) {return E*dSigmadERest(E);}
#endif

static double Spectra22A(char*name1,char*name2,double ** Spectra,txtList plusA)
{
  double vcs=0;
  int i,l;
  txtList cRec;
  int spin2Dm;
  
  for(l=0;l<6;l++) if(Spectra[l]) { Spectra[l][0]=Mcdm0; for(i=1;i<NZ;i++) Spectra[l][i]=0;}
  if(Spectra[0]==NULL && Spectra[1]==NULL) return 0;

if(Spectra[0])   
  for(l=0;l<nModelParticles;l++)
  {
     if(ModelPrtcls[l].NPDG==22) 
     { outNames[0]=ModelPrtcls[l].name;
       break;
     }
  } 
  if(l==nModelParticles) return 0;
  
  for(cRec=plusA;cRec;cRec=cRec->next)
  { 
     char lib[50]="";
     char*N[5];
     process2Lib(cRec->txt,lib);
     cc23=getMEcode(0,0,cRec->txt,NULL,NULL,lib);
     if(cc23==NULL) continue;
     if(passParameters(cc23)) continue; 
     *(cc23->interface->gtwidth)=0;
     *(cc23->interface->twidth)=0;
     *(cc23->interface->gswidth)=0;

     for(i=0;i<5;i++) N[i]=cc23->interface->pinf(1,i+1,pmass+i,code+i);
    
     cc23->interface->pinfAux(1,1,&spin2Dm,NULL,NULL,NULL);
     
     for(ix=0,i=2;i<5;i++) if(code[i]==22) iA=i; else if(!ix)ix=i;else iX=i;
     if(Spectra[0] && dSigmadE(X1*Mcdm0) > X1CUT*dSigmadE_x1 )
     {  double x2Sum; 
        dSigmadE_x0=3*X0*(dSigmadE(X0*Mcdm0)+dSigmadE(3*X0*Mcdm0)-2*dSigmadE(2*X0*Mcdm0));
#ifdef DISPLAY_SPECTRA
     {  char buff[100];
        sprintf(buff,"d(vSigma(%s))/dE(A) [pb/GeV]",cRec->txt);      
        displayFunc(dSigmadERest, Mcdm0*X0, Mcdm0,buff);
     }
#endif 
       x2Sum=0; 
       for(i=1;i<NZ;i++) 
       { double dSigmaDz=0,x=exp(Zi(i));
         if(x>X0)
         { dSigmaDz=x*Mcdm0*dSigmadERest(x*Mcdm0);
           Spectra[0][i]+=dSigmaDz;
           x2Sum+=(Zi(i)-Zi(i+1))*dSigmaDz;
         }else  Spectra[0][i]+=dSigmaDz;
       }
       { 
         vcs+=x2Sum;
         vSigmaCh=realloc(vSigmaCh, (nAnCh+2)*sizeof(aChannel));
         vSigmaCh[nAnCh].weight=x2Sum;   
         { int j; for(j=0;j<5;j++) vSigmaCh[nAnCh].prtcl[j]=N[j]; }
         nAnCh++;                     
       }
     }
     if(Spectra[1] && abs(code[iX])==11 && code[ix]+code[iX]==0  /*&& pmass[ix]==0*/ &&   // !!! was blocked for Me!=0. Reason not clear.  
       code[0]==code[1] && spin2Dm==1 && dSigmadEe(X1*Mcdm0) > X1CUT*dSigmadE_x1_e  ) 
     {
        
#ifdef DISPLAY_SPECTRA 
  displayFunc(dSigmadEe, Mcdm0*X0, Mcdm0," electron spectrum");
{ double csA,csE,xcsA,xcsE; 
  csA=simpson(dSigmadE,Mcdm0*X0,Mcdm0,1.E-3,NULL);
  csE=simpson(dSigmadEe,Mcdm0*X0,Mcdm0,1.E-3,NULL);
  xcsA=simpson(EdSigmadERest,Mcdm0*X0,Mcdm0,1.E-3,NULL);
  xcsE=simpson(EdSigmadEe,Mcdm0*X0,Mcdm0,1.E-3,NULL);
/*  
  printf("vcs(A)= %E  vcs(E)= %E\n",csA,csE);
  printf("energy  fraction lost %E\n", (2*Mcdm0- xcsA/csA-2*xcsE/csA)/Mcdm0);  
*/  
}        
#endif
       for(i=1;i<NZ;i++)
       {
           double Ee=Mcdm0*exp(Zi(i));
           if(Ee>X0*Mcdm0) Spectra[1][i]+=Ee*dSigmadEe(Ee);
       }
     }         
  }  
  return vcs;
}

/* ------------------------------------------------- */


static int readSpectra(void)
{ static int rdOk=0;
  int k,l,i,n;
  FILE *f;
  char * buff;
  char *fnames[Nin]={"gg.dat","dd.dat","uu.dat","ss.dat","cc.dat",
                     "bb.dat","tt.dat","ee.dat","mm.dat","ll.dat",
                     "zz.dat","zz_t.dat","zz_l.dat",
                     "ww.dat","ww_t.dat","ww_l.dat",
                     "pi0pi0.dat","pippim.dat","klkl.dat","ksks.dat","kpkm.dat"};

  if(rdOk) return 0;
  buff=malloc(strlen(micrO)+100);

  for(n=0;n<Nin;n++)
  {  sprintf(buff,"%s/Data/data/%s",micrO,fnames[n]);
     f=  fopen(buff,"r");
//     fclose(f); f=NULL;
     
     if(f==NULL) { free(buff);return 1;}    
     for(i=0;i<6;i++)
     { char ppp[20];
       fscanf(f,"%s",ppp);
       for(k=0;k<NEn;k++)
       {
         for(l=0;l<NZ;l++) if(1!=fscanf(f,"%f",phidiff[n][i][k]+l)) break;
         fscanf(f,"%*s"); 
         for(;l<NZ;l++)phidiff[n][i][k][l]=0;
       }
     }
//printf("ok before fclose\n");     
     fclose(f);
  }
    
  
  free(buff);
  rdOk=1;   

#ifdef ECTEST
printf("Energy conservation test\n");
for(n=0;n<Nin;n++) for(k=0;k<NEn;k++)
{
    double mi[NEn]={0.1,0.15,0.2,0.3,0.5,0.8,1.1,1.5,2,5,10,25,50,80.3,85,91.2,92,95,100,110,120,125,130,140,150,176,200,250,350,500,750,1000,1500,2000,3000,5000};

  double x[6];
  for(n=0;n<Nin;n++) for(k=0;k<NEn;k++)
  { double sum=0;


  for(l=0;l<6;l++) x[l]=0;
  for(l=0;l<6;l++) if(phidiff[n][l][k])
  { 
     for(i=0;i<NZ;i++)
     { double xx= pow(1.E-7, pow((i+0.5)/NZ,1.5));
       double dx;
       if(i==0) dx=-0.5*log(1.E-7)*pow((double)(1.)/NZ,1.5);
       else dx= 0.5*log(1.E-7)*( pow((double)(i-1)/NZ,1.5) - pow((double)(i+1)/NZ,1.5) );
       if(l==2) xx+=1/mi[k];
       x[l]+=phidiff[n][l][k][i]*dx*xx;   
     }
  }  
  sum=x[0]+2*(x[1]+x[2]+x[3]+x[4]+x[5]);
   
  if(sum && (sum>2.1 || sum<1.9))   
   printf("Channel = %s Energy=%.1E  %.2e + %.2e + %.2e + %.2e + %.2e + %.2e = %.2E\n",fnames[n],mi[k], x[0],x[1],x[2],x[3],x[4],x[5],sum);
  }
}  
   
#endif  

  
  return 0;
}


double zInterp(double zz, double * tab) 
{  
   double dz,r;
   int j0;   
   if(zz>0) return 0;
   
   j0=Iz(zz); 
   if(j0<1) j0=1;
   if(j0>=NZ-1) return tab[NZ-1];
   
   dz= (zz-Zi(j0))/(Zi(j0+1)-Zi(j0));
   r=(1-dz)*tab[j0]+dz*tab[j0+1];
   if(r<0)r=0;
   return r; 
}

//char *fnames[Nin]={"gg.dat","dd.dat","uu.dat","ss.dat","cc.dat","bb.dat","tt.dat","ee.dat","mm.dat","ll.dat","zz.dat","zz_t.dat","zz_l.dat","ww.dat","ww_t.dat","ww_l.dat"};
//                      0        1        2        3        4        5        6        7        8        9        10         11         12       13         14         15

static void mInterp(double Nmass,  int  CHin,int  CHout, double*tab)
{  
   double mi[NEn]={0.1,0.15,0.2,0.3,0.5,0.8,1.1,1.5,2 ,5 ,10,25,50,80.3,85,91.2,92,95,100,110,120,125,130,140,150,176,200,250,350,500,750,1000,1500,2000,3000,5000};
//                  0   1    2   3   4   5   6   7  8  9  10 11 12   13 14  15  16 17  18  19  20  21  22  23  24  25
   int l,i0;
   double c0,c1;
   float *p0,*p1;
   for(i0=0; i0<NEn && Nmass>=mi[i0] ;i0++);
   if(i0) i0--;

   int i00=-1;
   switch(CHin)  /* below pole mass */
   { 
     case  6:                      if(i0<25) i00=25; break; /* t m<17   */
     case  7: case 9:              if(i0<8)  i00=8;  break; /* e and tau starts M=2 */
     case 13: case 14: case 15:    if(i0<13) i00=13; break;   /* W m<80.3 */
     case 10: case 11: case 12:    if(i0<15) i00=15; break;  /* Z m<91.2 */
   } 

   if(i00>=0) {for(l=1;l<NZ;l++) tab[l]=phidiff[CHin][CHout][i00][l-1]; tab[0]=Nmass; return; }

   p0=phidiff[CHin][CHout][i0];
   p1=phidiff[CHin][CHout][i0+1];
   if(i0==NEn-1 || Nmass<=2) for(l=1;l<NZ;l++) tab[l]= p0[l-1];
   else
   {
     c1=(Nmass*Nmass -mi[i0]*mi[i0])/(mi[i0+1]*mi[i0+1] - mi[i0]*mi[i0]);
     c0=1-c1;
     for(l=1;l<NZ;l++) tab[l]= c0*p0[l-1]+c1*p1[l-1];
   }
   tab[0]=Nmass;
}



static double zInterpE(double x, double *tab ){ return    zInterp(x,tab)*exp(x);}  


double  spectrInfo(double Emin,double*tab,double*Etot)
{
  if(Etot) *Etot=0;
  if(Emin>=tab[0]) return 0;  else 
  { int i1,i2;
    double Xmin=Emin/tab[0],zmin,zmax=0;
    if(Xmin<1.22E-7) Xmin=1.22E-7;
    zmin=log(Xmin); 
    for(i1=Iz(zmin); i1>1 && tab[i1]==0;i1--) continue;
    if(i1<NZ-1) i1++;
    if(zmin<Zi(i1)) zmin=Zi(i1)+1E-6;

    
    for(i2=1; tab[i2]==0 && i2<NZ-1 ;i2++ );
    if(i2>1) i2--;
    if(zmax>Zi(i2)) zmax=Zi(i2)-1E-6;

    if(zmax<zmin) return 0;    
  
    if(Etot)*Etot=tab[0]*simpson_arg(zInterpE,tab, zmin, zmax,1.E-4,NULL);       
    return simpson_arg(zInterp,tab, zmin, zmax,1.E-4,NULL);
  }  
}


double spectrInt(double Emin,double Emax, double * tab)
{ 
  double M=tab[0], zmin,zmax;
  int i1,i2;

  if(Emin<M*exp(Zi(NZ-1))) zmin=Zi(NZ-1); else if( Emin>=M) zmin=0; else zmin=log(Emin/M);
  if(Emax<M*exp(Zi(NZ-1))) zmax=Zi(NZ-1); else if( Emax>=M) zmax=0; else zmax=log(Emax/M);
  
  if(zmin>=zmax) return 0;
  
  for(i1=Iz(zmin) ;i1>1 && tab[i1]==0;i1--) continue;
  if(i1<NZ-1) i1++;
  if(zmin<Zi(i1)) zmin=Zi(i1)+1E-6;
 
  i2=Iz(zmax)+1; if(i2>NZ-1) i2=NZ-1;
  for( ;i2<NZ-1 && tab[i2]==0;i2--) continue;
  if(i2<2) i2--;
  if(zmax>Zi(i2)) zmax=Zi(i2)-1E-6;
  
  return simpson_arg(zInterp,tab, zmin, zmax,1.E-4,NULL);
}


void  spectrMult( double *spect, double(*func)(double))
{ 
  int i; 
  double M=spect[0];
  for(i=1;i<NZ;i++)
  { double E=M*exp(Zi(i));
    spect[i]*=func(E);
  }  
}


/*
typedef struct
{ double m;
  double * tab;
} boostParStr;




static double FUNB(double e, void *B)
{ 
  boostParStr * par=B;
  double r;

  if(e<=0) return 0;
  r=zInterp(log(e/par->tab[0]),par->tab)/(e*sqrt(e*(e+2*par->m)));
  return r;
}

static double FUNB_(double e, void *B)
{ 
  boostParStr * par=B;
  double r;

  if(e<=0) return 0;
  r=zInterp(log(e/par->tab[0]),par->tab);
  return r;
}

static double FUNB__(double e, void *B)
{ 
  boostParStr * par=B;
  double r;

  if(e<=0) return 0;
  r=1/(e*sqrt(e*(e+2*par->m)));
  return r;
}

*/


void boost(double Y, double M0, double mx, double*tab)
{ double chY=cosh(Y), shY=sinh(Y);
  int l,k;
  double tab_out[NZ];  

  double tab_in[NZ];
  
  tab_in[0]=tab[0];
  for(int i=1;i<NZ;i++)
  { double e=tab[0]*exp(Zi(i));
    tab_in[i]=tab[i]/(sqrt(e*(e+2*mx)));
  }  
/*
static int N=0;
N++;  
if(N==11) displayPlot("boost","e",tab[0]/1E7, tab[0],1,1,"",0,SpectdNdE,tab);
*/
  if(Y<0.01) return;
  
//m12=1.073733E-03 pcm=5.249327e-02 Y=4.582795e+00  tab[0]=5.368667e-04 cosh(Y)*tab[0]=2.625212E-02   M0=5.250000E-02 M=7.046556E-02

  double e=tab[0]+mx, p=sqrt(e*e-mx*mx), M=e*chY+p*shY-mx;
//  printf("Y=%e  tab[0]=%e cosh(Y)*tab[0]=%E   M0=%E M=%E\n",Y,tab[0],cosh(Y)*tab[0], M0,M);
  if(M0<M) M0=M;     
 
  double emin=tab[0]*exp(Zi(NZ-1));
  int kmin=NZ-1;
  for(; tab[kmin]==0 && kmin>0; kmin--) continue;
  if(kmin==0) { tab[0]=M0;  return;}
  if(kmin<NZ-1) kmin++;  // tab[kmin]==0, tab[kmin+1]!=0;
  emin=tab[0]*exp(Zi(kmin));
   
  for(l=1;l<NZ;l++)
  { double e=M0*exp(Zi(l)),e1,e2;
    if(mx==0) { e1=e*exp(-Y); e2=e*exp(Y);}  
    else 
    {
      if(e>mx) 
      {
        double YY=acosh(1+e/mx);
        e1=mx*(cosh(Y-YY)-1);
        e2=mx*(cosh(Y+YY)-1);
//printf("Y =%E YY=%E  e1=%E e2=%E\n", Y,YY,e1,e2);        
      }else 
      {
        double p=sqrt(e*(e+2*mx)); 
        e1=chY*(e+mx)-shY*p-mx;
        e2=chY*(e+mx)+shY*p-mx;
      }
    }
    if(e1>=tab[0]) {tab_out[l]=0; continue;}
    if(e2>tab[0]) e2=tab[0];
    
//    if(e1 < emin) k=NZ-1; else k=Iz(log(e1/tab[0]))+1;
//    if(k>kmin) k=kmin;
//    e1=tab[0]*exp(Zi(k));
      if(e1<emin) e1=emin;
      if(e2<=emin)  {tab_out[l]=0; continue;}  
        
//    if(e1>=tab[0]) {tab_out[l]=0; continue;}    
//    if(e2>tab[0]) e2=tab[0];

    int err=0,err_=0;
    if (e1>=e2) tab_out[l]=0; else tab_out[l]=e*simpson_arg(zInterp,tab_in, log(e1/tab[0]),log(e2/tab[0]),1.E-3,&err_)/2/shY; 
  }
  
  for(l=1;l<NZ;l++) tab[l]=tab_out[l];
  tab[0]=M0;
//  printf("tab[0]=%E M0=%E\n",tab[0],M0);
}


static double outMass[Nout]= {0.,0.511E-3,0.939,0.,0.,0.};
char* outNames[Nout]={"gamma","e+","p-","nu_e","nu_mu","nu_tau"};


static int basicSpectraS(double Mass, int pdgN, int outN, double * tab)
{ 
//  printf("basicSpectra: pdg=%d outN=%d\n", pdgN, outN);
  tab[0]=Mass;
  for(int i=1;i<NZ;i++) tab[i]=0;
  int N=abs(pdgN);
      
   
  int inP=-1;

  switch(N)
  { case 21:inP=0; break;  /*glu*/ 
    case 1: case 2: case 3: case 4: case 5: case 6:  inP=N; break; /* d,u,s,c,b,t*/
    case 11:case 13: case 15: inP=(N+1)/2+1; break; /*e,m,l*/  
    case 23:    inP=10; break;  /*z*/
    case 23+'T':inP=11; break;  
    case 23+'L':inP=12; break;     
    case 24:    inP=13; break;  /*w*/
    case 24+'T':inP=14; break;
    case 24+'L':inP=15; break;
    case 111   :inP=16; break;
    case 211   :inP=17; break;
    case 130   :inP=18; break; // K0_L
    case 310   :inP=19; break; // k0_S
    case 321   :inP=20; break;
    case 12:case 14: case 16: case 22:
    {  if( (N==22 && outN==0) ||  (N-12 == 2*(outN-3)))
       {
         tab[1]=2/(-Zi(3));
         tab[2]=2/(-Zi(3))*( 1-Zi(2)/Zi(3));       
       }
       if(N==22) {tab[1]*=2;tab[2]*=2;}
       return 0;
    }
  }
  if(inP==-1) { tab[0]=Mass; for(int i=1;i<NZ;i++) tab[i]=0;    return 1;}
  readSpectra(); 
  mInterp(Mass,inP,outN,tab);
//for(int i=0;i<NZ;i++) printf("tab[%d]=%E\n",i,tab[i]);
  
/*
  if(pdgN==11)
  { 
     if(outN==0 && Mass<2)
     {  double me=5.11E-4;
        if(Mass<me) for(int i=1;i<NZ;i++) tab[i]=0; else 
        { 
          for(int i=1;i<NZ;i++)
          {
             double E=Mass*exp(Zi(i));
             tab[i]=2*E*FSRdNdE(E,Mass,me,1,1);    
          }  
        } 
     }
  }
*/      
  return 0;
}

#define Nch   25
#define NM    62
#define NX    179


static int basicSpectraPPPC(double Mass, int pdgN, int outN, double * tab)
{
   tab[0]=Mass;
   for(int i=1;i<NZ;i++) tab[i]=0;
   static int init=0;
   static double dmMass[NM], x[NX];

   static float phidiffPPPC[Nch][Nout][NM][NX]; 

   if(init==0)
   {
//printf("READ\n");
      char* src[Nout]=
      {
         "AtProduction_gammas.dat", "AtProduction_positrons.dat", "AtProduction_antiprotons.dat"
        , "AtProduction_neutrinos_e.dat", "AtProduction_neutrinos_mu.dat", "AtProduction_neutrinos_tau.dat"
//        , "AtProduction_antideuterons.dat"            
      };
      for(int i=0;i<Nch;i++) for(int j=0;j<Nout;j++) for(int k=0;k<NM;k++) for(int l=0;l<NX;l++) phidiffPPPC[i][j][k][l]=0;

      char* buff;
      buff=malloc(strlen(micrO)+100);

      for(int k=0;k<Nout;k++) 
      { 
         sprintf(buff,"%s/Data/data/PPPC/%s",micrO,src[k]); 
         FILE*F=fopen(buff,"r"); 

//printf("F=%p\n",F);
         fscanf(F,"%*[^\n]\n");
         for(int m=0;m<NM;m++) for(int ix=0;ix<NX;ix++)
         {  fscanf(F,"%lf %lf",dmMass+m,x+ix);
//             printf("M[%d]=%E x[%d]=%E\n",m,dmMass[m],ix,pow(10,x[ix]));     
            for(int ich=0;ich<Nch;ich++) fscanf(F,"%f",&(phidiffPPPC[ich][k][m][ix]));
            fscanf(F,"%*[^\n]\n");
         }
         fclose(F);
               
         for(int m=0;m<NM;m++) for(int ix=0;ix<NX;ix++)for(int ich=0;ich<Nch;ich++) phidiffPPPC[ich][k][m][ix]/=log(10);
      }
//      exit(0);
      free(buff);
      for(int ix=0;ix<NX;ix++)  x[ix]=pow(10,x[ix]);           
      init=1;
   } 

//eL,eR, e,mL,mR, m,lL,lR, l, q, c, b, t,WL,WT, W,ZL,ZT, Z, g, A, h, nu_e,nu_m,nu_l
//0   1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20  21  22  23   24 

  int ch=-1;
  switch(pdgN) 
  { case 1  :  case 2 : case 3: 
                  ch=9 ; break;
    case 4      : ch=10; break;
    case 5      : ch=11; break;
    case 6      : ch=12; break;
    case 11+'L' : ch=0 ; break;
    case 11+'R' : ch=1 ; break;
    case 11     : ch=2 ; break;
    case 13+'L' : ch=3 ; break;
    case 13+'R' : ch=4 ; break;
    case 13     : ch=5 ; break;
    case 15+'L' : ch=6 ; break;
    case 15+'R' : ch=7 ; break;
    case 15     : ch=8 ; break;
    case 21     : ch=19; break;
    case 22     : ch=20; break;
    case 24+'L' : ch=13; break;
    case 24+'T' : ch=14; break;
    case 24     : ch=15; break;
    case 23+'L' : ch=16; break;
    case 23+'T' : ch=17; break;
    case 23     : ch=18; break;
    case 25     : ch=21; break;
    case 12     : ch=22; break;
    case 14     : ch=23; break;
    case 16     : ch=24; break;
  }
  if(ch<0) return 1;
  int m;
  for(m=0;m<NM && Mass>dmMass[m];m++) continue;
  m--;
  double stab[NX];
  if(m==NM-1) for(int ix=0;ix<NX;ix++) stab[ix]=phidiffPPPC[ch][outN][m][ix]; 
  else
  {  double a=(Mass-dmMass[m])/(dmMass[m+1]-dmMass[m]);
     for(int ix=0;ix<NX;ix++)  stab[ix]=(1-a)*phidiffPPPC[ch][outN][m][ix]+a*phidiffPPPC[ch][outN][m+1][ix];   
  }

  for(int i=1;i<NZ;i++) tab[i]=polint1(exp(Zi(i)), NX, x, stab);
//for(int i=0;i<NX;i++) printf("x=%e f=%e\n", x[i],stab[i]);  
   
  return 0;
}

#undef Nch  // 25
#undef NM   // 62
#undef NX   // 179   




#define Nch   14
#define NM    60
#define NX    180   


static int basicSpectraA(double Mass, int pdgN, int outN, double * tab)
{

   tab[0]=Mass;
   for(int i=1;i<NZ;i++) tab[i]=0;
   static int init=0;
   static double dmMass[NM], x[NX];

   static float phidiffA[Nch][Nout][NM][NX]; 

   if(init==0)
   {
//printf("READ\n");
      char* src[Nout]={"AtProduction-Ga.dat",  "AtProduction-Positrons.dat", "AtProduction-AntiP.dat",  
                    "AtProduction-Nuel.dat","AtProduction-Numu.dat",      "AtProduction-Nuta.dat"};
      for(int i=0;i<Nch;i++) for(int j=0;j<Nout;j++) for(int k=0;k<NM;k++) for(int l=0;l<NX;l++) phidiffA[i][j][k][l]=0;

      char* buff;
      buff=malloc(strlen(micrO)+100);

      for(int k=0;k<Nout;k++) 
      { 
         sprintf(buff,"%s/Data/data/woUncertainty/%s",micrO,src[k]); 
         FILE*F=fopen(buff,"r"); 

//printf("F=%p\n",F);
         fscanf(F,"%*[^\n]\n");
         for(int m=0;m<NM;m++) for(int ix=0;ix<NX;ix++)
         {  fscanf(F,"%lf %lf",dmMass+m,x+ix); /*printf("M=%E x=%E\n",dmMass[m],x[ix]);*/     for(int ich=0;ich<Nch;ich++){    fscanf(F,"%f",&(phidiffA[ich][k][m][ix]));}
         }
            for(int m=0;m<NM;m++) for(int ix=0;ix<NX;ix++)for(int ich=0;ich<Nch;ich++) phidiffA[ich][k][m][ix]/=log(10);
            fclose(F);   
      }
      free(buff);
      init=1;
   } 


  int ch=-1;
  switch(pdgN) 
  { case 1  : ch=4 ; break;
    case 2  : ch=3 ; break;
    case 3  : ch=5 ; break;
    case 4  : ch=6 ; break;
    case 5  : ch=7 ; break;
    case 6  : ch=8 ; break;
    case 11 : ch=0 ; break;
    case 13 : ch=1 ; break;
    case 15 : ch=2 ; break;
    case 21 : ch=12; break;
    case 22 : ch=9 ; break;
    case 23 : ch=11; break;
    case 24 : ch=10; break;
    case 25 : ch=13; break;
    case 12:case 14: case 16:
    {  if(pdgN-12 == 2*(outN-3))
       {
         tab[1]=2/(-Zi(3));
         tab[2]=2/(-Zi(3))*( 1-Zi(2)/Zi(3));       
       }
       return 0;
    }
  }
   
  if(ch<0) return 1;
  int m;
  for(m=0;m<NM && Mass>dmMass[m];m++) continue;
  m--;
  double stab[NX];
  if(m==NM-1) for(int ix=0;ix<NX;ix++) stab[ix]=phidiffA[ch][outN][m][ix]; 
  else
  {  double a=(Mass-dmMass[m])/(dmMass[m+1]-dmMass[m]);
     for(int ix=0;ix<NX;ix++)  stab[ix]=(1-a)*phidiffA[ch][outN][m][ix]+a*phidiffA[ch][outN][m+1][ix];   
  }

  for(int i=1;i<NZ;i++) tab[i]=polint1(exp(Zi(i)), NX, x, stab);

    
  if( pdgN==22 && outN==0 )
   {
         tab[1]=4/(-Zi(3));
         tab[2]=4/(-Zi(3))*( 1-Zi(2)/Zi(3));       
   }
  
  return 0;
}

#undef Nch  // 14
#undef NM   // 60
#undef NX   // 180   

int basicSpectra(double Mass, int pdgN, int outN, double * tab)
{ 
   if(SpectraFlag==0) return basicSpectraS(Mass,pdgN,outN,tab);
   if(SpectraFlag==1) return basicSpectraA(Mass,pdgN,outN,tab);
   if(SpectraFlag==2) return basicSpectraPPPC(Mass,pdgN,outN,tab);   
}



/*static*/ void getSpectrum2(int wPol, double M, char*n1,char*n2,int outP, double *tab);

void  decaySpectrum(char*pName,double M, int outP, double*tabD)
{ 
//  double M=pMass(pName); 
  double mx=outMass[outP];
  
  tabD[0]=M/2;
  for(int i=1;i<NZ;i++)tabD[i]=0;
/*          
  if( (CDM1 &&(strcmp(CDM1,pName)==0 ||strcmp(aCDM1,pName)==0))
    ||(CDM2 &&(strcmp(CDM2,pName)==0 ||strcmp(aCDM2,pName)==0))) return; 
*/    
  txtList L;   
  double w=pWidth(pName,&L);
  if(w==0) 
  { char  mess[100]; sprintf(mess," Can not decay '%s'; ",pName);
    addErrorMess(&errMess,mess);
    return;    
  }     
  
//double E,Es=0;
  double brs=0;          
  for(;L;L=L->next)
  { char p[5][20];
    double tab_p[NZ];
    double br;
    int n=sscanf(L->txt,"%lf %s -> %[^, ], %[^,], %[^,], %[^, ]",&br, p[0],p[1],p[2],p[3],p[4]);    
    n-=2;
//printf("%s %E\n",L->txt,br);    
    if(n==4) continue;
    else if(n==2) 
    { getSpectrum2(0,M,p[1],p[2],outP, tab_p);    
//      spectrInfo(1E-10,tab_p,&E);
//      Es+=E*br;
//      printf("E=%e Es=%E\n", E,Es);
    }  
    else if(n==3) 
    {  char process[50];
       sprintf(process,"%s->%s,%s,%s",p[0],p[1],p[2],p[3]); 
// printf("process=%s\n",process);       
       numout*cc13=newProcess(process);
       if(!cc13) { printf("Wrong decay channel passed via SLHA: %s\n", process); continue;}
       passParameters(cc13);
       argFor13Int arg13;
       arg13.cc=cc13;
       arg13.nsub=1;
       for(int i=0;i<4;i++) arg13.cc->interface->pinf(1,i+1,arg13.mass+i,NULL); 
       arg13.GG =sqrt(4*M_PI*alphaQCD(arg13.mass[0]));
                   
       tab_p[0]=tabD[0];
       for(int i=1;i<NZ;i++) tab_p[i]=0;
       
       int i1,i2,i3;
       i3=1;
       
       for(int i=2;i<4;i++)
       {
          int pdg=abs(pNum(p[i]));
          switch(pdg)
          { case 22: if(outP==0) continue; else break;
            case 11: if(outP==1) continue; else break;
            case 12: if(outP==3) continue; else break;
            case 14: if(outP==4) continue; else break;
            case 16: if(outP==5) continue; else break; 
          }  
          if(arg13.mass[i]>arg13.mass[i3]) i3=i;
       }  
//i3=1;    
//printf("p[i3]=%s\n",p[i3]);   
       for(i1=1;;i1++) if(i1!=i3) break;
       for(i2=i1+1;;i2++) if(i2!=i3) break; 
       double mMin=0;
       for(int i=1; i<4;i++) if(i!=i3) {mMin+=arg13.mass[i]; if(abs(pNum(p[i]))<3) mMin+=0.14;} 
       if(mMin<2*mx) mMin=2*mx;
       double mMax=arg13.mass[0]-arg13.mass[i3]; if(abs(pNum(p[i3]))<3) mMax-=0.14;
       arg13.i3=i3;
       if(mMin>=mMax) continue;
       double wt=0;
       int nInt=20;       
       for(int k=0;k<nInt;k++)
       { 
         double m12=mMin+ (k+0.5)*(mMax-mMin)/nInt;
         double dw=dWidth13dM(m12, &arg13)*(mMax-mMin)/nInt;
          
         wt+=dw;
         double pcm=decayPcm(arg13.mass[0], arg13.mass[i3], m12);            
         double tab12[NZ],tabD[NZ];

         getSpectrum2(0,m12,p[i1],p[i2],outP, tab12);
   
         double Y=asinh(pcm/m12);
         boost(Y, M/2, outMass[outP], tab12);

         double m=arg13.mass[i3];
         double e=sqrt(pcm*pcm+m*m);
         if(0==basicSpectra(e,abs(pNum(p[i3])),outP,tabD)) for(int i=1;i<NZ;i++) tabD[i]/=2; 
         else 
         {  
            decaySpectrum(p[i3],m, outP, tabD);
            Y=asinh(pcm/m);
            boost(Y, M/2, outMass[outP], tabD);
         }  
           addSpectrum(tab12, tabD);
           for(int i=1;i<NZ;i++) tab12[i]*=dw;  
           addSpectrum(tab_p, tab12);    
       }
//printf("wt=%e\n", wt);       
       for(int i=1;i<NZ;i++) tab_p[i]/=wt;
    }
    for(int i=1;i<NZ;i++) tab_p[i]*=br; 
    {
       double buff[NZ];
       for(int i=0;i<NZ;i++) buff[i]=tabD[i];
       addSpectrum(tabD,tab_p);
       
    }    
    brs+=br;
  }
}


void getSpectrum2(int pol, double M,char*n1,char*n2,int outP, double *tab)
{
  double m1=pMass(n1);
  double m2=pMass(n2);
  int N1=pNum(n1);
  int N2=pNum(n2);  
//printf("getSpectrum2: %s %s -> %d\n", n1,n2,outP);
  int i;
  int inP=-1;
  int N;

  if(M<=m1+m2) pol=0;
  if(pol!='T' && pol!='L') pol=0;
  tab[0]=M/2;
  for(i=1;i<NZ;i++) tab[i]=0;
   
  if(abs(N1)==abs(N2)) 
  {    
     N=abs(N1);
     if(N==23 || N==24)  N+=pol; 
     if(0==basicSpectra(M/2, N, outP, tab)) return;
  }  
         
    char* nn[2];
    double mm[2];
    double E[2]; 
    double pcm;
    int k;

    nn[0]=n1;
    nn[1]=n2;

    mm[0]=m1;
    mm[1]=m2;
    double tabAux[NZ];
    if(M>m1+m2) pcm=decayPcm(M,m1,m2);
    else 
    { pcm=0;
      if(abs(N1)==abs(N2)) { mm[0]=M/2; mm[1]=M/2;} else
      {
        if(N1==23 || abs(N1)==24) mm[0]=M-m2; else mm[1]=M-m1;
      }  
    } 
    E[0]=sqrt(mm[0]*mm[0]+pcm*pcm);
    E[1]=sqrt(mm[1]*mm[1]+pcm*pcm);
    
    for(i=1;i<NZ;i++) tab[i]=0;
    
    for(k=0;k<2;k++)
    { N=pNum(nn[k]);
      if(N==23||N==24) N+=pol;
      if(0==basicSpectra(E[k],N,outP,tabAux))  
      { 
         for(int i=1;i<NZ-1;i++) tabAux[i]/=2;
          addSpectrum(tab, tabAux);
      }     
      else 
      {  if(mm[k]==0)
         {
           char  mess[100]; sprintf(mess, "No hadronization  for massless '%s'  ",nn[k]);
           addErrorMess(&errMess,mess);
         }  
         else
         {
          decaySpectrum(nn[k],mm[k],outP, tabAux);
          double Y=asinh(pcm/mm[k]);
          boost(Y, M/2, outMass[outP], tabAux); 
          for(i=1;i<NZ;i++)tab[i]+=tabAux[i];
        }   
      } 
    }  
}


void addSpectrum(double *Spect, double * toAdd)
{
  double m1=Spect[0];
  double m2=toAdd[0];
  double buff[NZ];
  int i;
  
  if(fabs(m1-m2)<5E-3*(fabs(m1)+fabs(m2))) { for(int i=1; i<NZ; i++) Spect[i]+=toAdd[i]; return; }
  

  if(m1>m2) for(i=NZ-1;i>0;i--) 
  { 
     double E= m1*exp(Zi(i));
     if(E>m2) return; 
     Spect[i]+=  SpectdNdE(E,toAdd)*E;  
  } else  
  { for(i=0;i<NZ;i++) buff[i]=Spect[i];
    for(i=0;i<NZ;i++) Spect[i]=toAdd[i];
    for(i=NZ-1;i>0;i--) 
    {  
       double E= m2*exp(Zi(i)); 
       if(E>m1) return;
       Spect[i]+=  SpectdNdE(E,buff)*E;
    } 
  }
} 




static double calcSpectrum0(char *name1,char*name2, int key, double **Spectra, txtList*plusA)
{
  int i,k,l;
  double vcsSum=0; 
  int ntot,err;
  double * v_cs=NULL;
  char * photonName=NULL;

  numout * libPtr=NULL;
  
  Mcdm0=0.5*(pMass(name1)+pMass(name2));  
  for(l=0;l<6;l++) if(Spectra[l]) { Spectra[l][0]=Mcdm0; for(i=1;i<NZ;i++) Spectra[l][i]=0;}  
  if(plusA){ photonName=pdg2name(22); if(!photonName) plusA=NULL; } 
  
  if(Qaddress && *Qaddress!=2*Mcdm0) { *Qaddress=2*Mcdm0; calcMainFunc(); }   

//========   Dm,Dm->Even,1x ==============
    {
       char name1L[10],name2L[10], lib[20];
       char*even=EvenParticles();
       char* process=malloc(100+strlen(even));
       pname2lib(name1,name1L);
       pname2lib(name2,name2L);
       sprintf(lib,"omg_%s%s",name1L,name2L); 
       sprintf(process,"%s,%s->AllEven,1*x{%s",name1,name2,even);    
       libPtr=getMEcode(0,ForceUG,process,NULL,NULL,lib); 
       free(process); 
       if(libPtr) 
       {  passParameters(libPtr); 
          (*libPtr->interface->twidth)=0;
          procInfo1(libPtr,&ntot,NULL,NULL); 
          v_cs=malloc(sizeof(double)*ntot);
       }
       else 
       {  ntot=0;
          v_cs=NULL;
       }   
     }    
  
  for(k=0;k<ntot;k++)
  { REAL m[4];
    double wV,br;
    char *N[4];
    int pdg[4];
    int l,l_;
    for(i=0;i<4;i++) N[i]=libPtr->interface->pinf(k+1,i+1,m+i,pdg+i);
    cc23=NULL;
    v_cs[k]=0;
    if(VZdecay||VWdecay)
    {  int nVV;
       int vd[4]={0,0,0,0};
       for(l=2;l<4;l++) if((pdg[l]==23&&VZdecay) || (abs(pdg[l])==24&&VWdecay)) vd[l]=1;
            
       for(l=2;l<4;l++) if(vd[l]) break;
       if(l<3)
       {  l_=5-l; 
          if(vd[l_])
          { nVV=2;
            if(m[l_]>m[l]) { l=l_; l_=5-l;}
          } else nVV=1;
          
          if(m[0]+m[1] >  m[l_] +20  && m[0] + m[1] <  m[2]+m[3] + 4*nVV)
           cc23=xVtoxll(2,2,N,pdg,l,&wV,&br);                
       }
    }
    if(cc23)
    { int i3W;  
      double  r,m1,v0=0.001;
      passParameters(cc23);
      for(i3W=2;i3W<5;i3W++) if(strcmp(cc23->interface->pinf(1,i3W+1,NULL,NULL),N[l_])==0) break;
      { int err;
        r=v0*cs23(cc23,1,v0*Mcdm0/2,i3W,&err)/br;
        if(err) printf("error in simpson spectra.c line 902\n");
      }  
      if(pdg[l_]==23 || abs(pdg[l_])==24)
      { double wV2;
        
        wV2=pWidth(N[l_],NULL);
        r*=decayPcmW(2*Mcdm0,m[l],m[l_],wV,wV2,0)/decayPcmW(2*Mcdm0,m[l],m[l_],wV,0.,0);
        if(pdg[l]==pdg[l_]) r/=2;
      }
      v_cs[k]=r; 
      vcsSum+=r;                     
    }
    else  if(m[2]+m[3]< m[0]+m[1])
    {  err=0;
#ifdef V0
      double pcm=V0*m[0]*m[1]/(m[0]+m[1]);
      double csGeV=cs22(libPtr,k+1,pcm,-1.,1.,&err)/3.8937966E8;
      improveCrossSection(pdg[0],pdg[1],pdg[2],pdg[3],pcm,&csGeV);    
      v_cs[k]=V0*csGeV*3.8937966E8;  
#else 
      v_cs[k]= vcs22(libPtr,k+1,&err);
#endif 
/*
printf("%s %s => %s %s  %E \n", libPtr->interface->pinf(k+1,1,NULL,NULL)
                              , libPtr->interface->pinf(k+1,2,NULL,NULL)
                              , libPtr->interface->pinf(k+1,3,NULL,NULL)
                              , libPtr->interface->pinf(k+1,4,NULL,NULL)
                              , v_cs[k]    );
*/
     if(v_cs[k]<0) v_cs[k]=0; 
      vcsSum+=v_cs[k];
    } else v_cs[k]=-1;
  }
   
  for(k=0;k<ntot ;k++) if(v_cs[k]>=0)
  { char * N[4];
    REAL m[4];
    int l,charge3[2],spin2[2],pdg[2];
    int PlusAok=0;

    procInfo2(libPtr,k+1,N,m);
    
    for(l=0;l<2;l++)  pdg[l]=qNumbers(N[2+l],spin2+l,charge3+l,NULL);
    if(Spectra[0] && plusA && (charge3[0] || charge3[1])&& (m[2]+m[3]< m[0]+m[1])) 
    {
       double m1=m[2], m2=m[3], Eg=X1*Mcdm0;
       double kappa=4*Mcdm0*(Mcdm0-Eg), ms=m1+m2, md=m1-m2;
       if(ms*ms<kappa)
       {  double dp=(Mcdm0-Eg/2)*sqrt((1-ms*ms/kappa)*(1-md*md/kappa));
          double p1= Eg/2*(1+ms*md/kappa)-dp, E1=sqrt(p1*p1+m1*m1);
          double p2= Eg/2*(1-ms*md/kappa)-dp, E2=sqrt(p2*p2+m2*m2);
             
          double Q1=Mcdm0*Mcdm0+m1*m1-2*Mcdm0*E1;
          double Q2=Mcdm0*Mcdm0+m2*m2-2*Mcdm0*E2;             

          double  m_min=10*Mcdm0;
          int n,m,w;
          char *s;
/*             
printf("energy conservation:0=%E=%E\n",  (E1+Eg+sqrt(pow(p2+2*dp,2)+m2*m2))/2/Mcdm0-1   ,
                                                   (E2+Eg+sqrt(pow(p1+2*dp,2)+m1*m1))/2/Mcdm0-1 );
*/                      
          for(n=1;(s=libPtr->interface->den_info(k+1,n,&m,&w,NULL));n++)
          {   double mass=0; if(m) mass=fabs( libPtr->interface->va[m]);
                if( ( strcmp(s,"\1\3")==0  || strcmp(s,"\1\4")==0) && m  && m_min> mass) m_min=mass;
          }

          if(  m_min*m_min -Q1  < QCUT*Mcdm0*Mcdm0*abs(charge3[0])/3.  
            || m_min*m_min -Q2  < QCUT*Mcdm0*Mcdm0*abs(charge3[1])/3. 
            || (abs(pdg[0])==24 && abs(pdg[1])==24  && Mcdm0 > 500 && v_cs[k]>1.E-3*vcsSum ))
          { txtList new22A=malloc(sizeof(txtListStr));
               new22A->next=*plusA; 
               new22A->txt=malloc(50);
               *plusA=new22A;
               sprintf(new22A->txt ,"%s,%s->%s,%s,%s",N[0],N[1],N[2],N[3],photonName); 
               PlusAok=1; 
          }
       } 
    }
            
    if(v_cs[k]>0) 
    {  double tab2[NZ]; 
       int N3=pdg[0], N4=pdg[1];

//       if(PrintOn )
//       { char txt[100];
//         sprintf(txt,"%s,%s -> %s %s", N[0],N[1],N[2],N[3]);
//         printf("  %-20.20s  %.2E\n",txt,v_cs[k]*2.9979E-26);
//       }
       
       vSigmaCh=realloc(vSigmaCh, (nAnCh+2)*sizeof(aChannel));
       vSigmaCh[nAnCh].weight=v_cs[k];
       { int j; 
         for(j=0;j<4;j++) vSigmaCh[nAnCh].prtcl[j]=N[j];
         vSigmaCh[nAnCh].prtcl[4]=NULL;
       }
       nAnCh++;
       { double lng=0;
         int noP=-1;
         if(key&1 && (abs(pdg[0])==23 || abs(pdg[0])==24|| abs(pdg[1])==23 || abs(pdg[1])==24)) noP=vPolar(N, &lng);
//if((key&1) && noP>=0 ) printf(" %s %s %s %s noP=%d lng=%E \n",  N[0],N[1],N[2],N[3],noP,lng   );
         
         for(l=0;l<6;l++) if(Spectra[l])
         {
           if(noP) getSpectrum2(0,m[0]+m[1],N[2],N[3],l,tab2);
           else
           {
             double tabT[NZ],tabL[NZ]; 
             getSpectrum2('T',m[0]+m[1],N[2],N[3],l,tabT);
             getSpectrum2('L',m[0]+m[1],N[2],N[3],l,tabL);   
             for(int i=1;i<NZ;i++) tab2[i]=(1-lng)*tabT[i]+lng*tabL[i];
             tab2[0]=tabT[0];
           }
           for(i=1;i<NZ;i++) Spectra[l][i]+=tab2[i]*v_cs[k];
         }
       }
#ifdef ADDFSR

       if(Spectra[0] && charge3[0])
       {
          for(l=0;l<2;l++) if(addGamma(pdg[l])&& m[2+l]!=0.) for(i=1;i<NZ;i++)
          {  double x=exp(Zi(i));
             if(2*Mcdm0*sqrt(1-x) > m[2]+m[3])
             {  double pcm,kappa,one_kappa,f;
                pcm=decayPcm(2*Mcdm0*sqrt(1-x), m[2],m[3]); 
                kappa=sqrt(1/(1+m[2+l]*m[2+l]/pcm/pcm)); 
                if( m[2+l]/pcm > 1.E-2) one_kappa=1-kappa; 
                                   else one_kappa=0.5*m[2+l]*m[2+l]/pcm/pcm; 
               f=(1./137.)*charge3[l]*charge3[l]/9/M_PI*log((1+kappa)/one_kappa);
          
               if(spin2[l]&1) Spectra[0][i]+=f*(1-x*(1-x*0.5))* v_cs[k];
               else           Spectra[0][i]+=f*(1-x          )* v_cs[k];
             }  
          }
       }
#endif       
    }
  }
  if(v_cs)free(v_cs);  

  if(Ncdm >1)
  {
      char*inNames[3]={name1,name2,NULL}, *outNames[3];
       
      outNames[0]=OddParticles(0);
      outNames[1]=outNames[0];
      outNames[2]=NULL;
      
      txtList proc= makeProcList(inNames, outNames,0);  
      for(txtList p=proc;p;p=p->next)
      {    
         char*c=strstr(p->txt,"->");
         char N3[20], N4[20]; 
         sscanf(c+2,"%[^,],%s\n",N3,N4);
         trim(N3);trim(N4);   
         if(pMass(name1)+pMass(name2)<= pMass(N3)+pMass(N4)) continue;
         int dm3=0,dm4=0;
         for(int i=1;i<=Ncdm;i++) if(strcmp(CDM[i],N3)==0 || strcmp(CDM[i],antiParticle(N3))==0)   dm3=1;
         for(int i=1;i<=Ncdm;i++) if(strcmp(CDM[i],N4)==0 || strcmp(CDM[i],antiParticle(N4))==0)   dm4=1; 
//         if(dm3 && dm4) continue;
         numout *cc=newProcess(p->txt);
         passParameters(cc);
         REAL m[4];
         int pdg[4];
         char*N[4];
         for(i=0;i<4;i++) N[i]=cc->interface->pinf(1,i+1,m+i,pdg+i);
         double pcm=V0*m[0]*m[1]/(m[0]+m[1]);
         double csGeV=cs22(cc,1,pcm,-1.,1.,&err)/3.8937966E8;
         improveCrossSection(pdg[0],pdg[1],pdg[2],pdg[3],pcm,&csGeV);
         if(csGeV>0) 
         {  double vcs=V0*csGeV*3.8937966E8; 
            for(l=0;l<6;l++) if(Spectra[l])
            { double tab2[NZ];
              getSpectrum2(0,m[0]+m[1],N3,N4,l,tab2);
              for(i=1;i<NZ;i++) Spectra[l][i]+=tab2[i]*vcs;
            }
            vcsSum+=vcs;

            vSigmaCh=realloc(vSigmaCh, (nAnCh+2)*sizeof(aChannel));
            vSigmaCh[nAnCh].weight=vcs;
            { int j; 
              for(j=0;j<4;j++) vSigmaCh[nAnCh].prtcl[j]=N[j];
              vSigmaCh[nAnCh].prtcl[4]=NULL;
            }
            nAnCh++; 
         }   
      }
      cleanTxtList(proc);
  }  
  return  vcsSum;
}


double calcSpectrum(int key, double *Sg,double*Se, double*Sp, double*Sne,double*Snm,double*Snl, int *errcode)
{ int n,i,j,l,err;
  char  lop[20];
  double *Spectra[6],*Spectra_[6];
  double  buff[6][NZ];
  char * name, *aname;
  txtList plusA=NULL;
  txtList*plusAptr;

  if(SpectraFlag!=0 && key&1) { printf("If SpectraFlag==1, the  W/Z polarization is not taken into account \n"); key-=1;} 

  double vcs=0;
  
  if(errMess) { free(errMess); errMess=NULL;}
  
//   int checkGam[4][4]; ???

/*
  err=readSpectra(); 
  if(err) { printf("calcSpectrum: Can not read data files for spectra\n");
            if(errcode) *errcode=-1; 
            return 0;
          }
*/

  Mcdm0=0;
  for(int i=1;i<=Ncdm;i++) if(McdmN[i]>Mcdm0) Mcdm0=McdmN[i];
  
  Spectra[0]=Sg; Spectra[1]=Se; Spectra[2]=Sp,Spectra[3]=Sne,Spectra[4]=Snm; Spectra[5]=Snl;
  for(l=0;l<6;l++) if(Spectra[l])   Spectra_[l]=buff[l];  else Spectra_[l]=NULL;
  for(l=0;l<6;l++) if(Spectra[l]) { Spectra[l][0]=Mcdm0; for(i=1;i<NZ;i++) Spectra[l][i]=0;} 

  nAnCh=0;
  vSigmaCh=realloc(vSigmaCh, 2*sizeof(aChannel));
  if(errcode) *errcode=0;
    
  if(key&4) { PrintOn=1; printf("    Channel          vcs[cm^3/s]\n");} else PrintOn=0;  

  double Meff=0;
  for(int i=1;i<=Ncdm;i++) Meff+=fracCDM[i]/McdmN[i];
  Meff=1/Meff;
  
  for(int k1=1;k1<=Ncdm;k1++) for (int k2=k1;k2<=Ncdm;k2++) 
  {  double w=fracCDM[k1]*fracCDM[k2]*(Meff/McdmN[k1])*(Meff/McdmN[k2]) ;
     if(!w) continue;
     if(k1!=k2) w*=2;
     char * p1[2]={CDM[k1],antiParticle(CDM[k1])};
     char * p2[2]={CDM[k2],antiParticle(CDM[k2])};
     double w1[2]={ 0.5*(1+dmAsymm),0.5*(1-dmAsymm)};
     double w2[2]={ 0.5*(1+dmAsymm),0.5*(1-dmAsymm)}; 
     if(p1[0]==p1[1]) { w1[0]+=w1[1]; w1[1]=0;}
     if(p2[0]==p2[1]) { w2[0]+=w2[1]; w2[1]=0;}
     double ww[2][2]={ { w1[0]*w2[0], w1[0]*w2[1]}, {w1[1]*w2[0],w1[1]*w2[1]}};
     ww[0][1]+=ww[1][0]; ww[1][0]=0;
     if(k1==k2) { ww[0][0]+=ww[1][1]; ww[1][1]=0;}
     if(0==assignVal("Q",McdmN[k1]+McdmN[k2])) calcMainFunc();  // ?? Qstat?? 
     
     for(int i1=0;i1<2;i1++) for(int i2=0;i2<2;i2++) if(ww[i1][i2]) 
     { 
        double c=w*ww[i1][i2];
        int nAnCh0=nAnCh;       
        int key0=key&1;
        if(key&2 /*&& checkGam[i][j]*/ ) key0+=2;  // ???       
        
        Mcdm0=0.5*(pMass(p1[i1])+pMass(p2[i2]));
        if(key&2) plusAptr=&plusA; else plusAptr=NULL; 
        vcs+=c*calcSpectrum0(p1[i1],p2[i2], key0, Spectra_,plusAptr);
        for(l=0;l<6;l++) if(Spectra[l])
        {  for(int k=1;k<NZ;k++)Spectra_[l][k]*=c;
           addSpectrum(Spectra[l], Spectra_[l]);
        }
        for(int k=nAnCh0;k<nAnCh;k++) vSigmaCh[k].weight*=c;   

       if(plusA)
       {  nAnCh0=nAnCh;
          if(Spectra[0]) dSigmadE_x1=zInterp(log(X1),Spectra_[0])/(X1*Mcdm0);  
          if(Spectra[1]) dSigmadE_x1_e=zInterp(log(X1),Spectra_[1])/(X1*Mcdm0);
          vcs+=c*Spectra22A(name,aname,Spectra_,plusA);
          
         for(l=0;l<6;l++) if(Spectra[l])
         { 
           for(int k=1;k<NZ;k++)Spectra_[l][k]*=c;
           addSpectrum(Spectra[l], Spectra_[l]);
         }
        
         for(int k=nAnCh0;k<nAnCh;k++) vSigmaCh[k].weight*=c;
         cleanTxtList(plusA);
         plusA=NULL;
       }    
                  
     }
  }

//sorting
  vSigmaCh[nAnCh].weight=0;
  for(j=0;j<5;j++) vSigmaCh[nAnCh].prtcl[j]=NULL;
  for(i=0;i<nAnCh-1;)
  {  if(vSigmaCh[i].weight >= vSigmaCh[i+1].weight) i++; 
     else
     {  aChannel buff;
        buff=vSigmaCh[i+1];vSigmaCh[i+1]=vSigmaCh[i];vSigmaCh[i]=buff;
        if(i)i--;else i++;
     }
  }           
// normalization 
  if(vcs)
  {  for(l=0;l<6;l++) if(Spectra[l])for(i=1;i<NZ;i++)Spectra[l][i]/=vcs;
     for(i=0;i<nAnCh;i++) vSigmaCh[i].weight/= vcs;
  } 
  
  
  vcs*=2.9979E-26;
  

  if(PrintOn )
  { int i=0;
    char txt[100];
    printf("==================================\n annihilation cross section %.2E cm^3/s\n",vcs  );
    printf(" contribution of processes\n");
//double sss=0;         
    for(i=0;vSigmaCh[i].weight>0*1.E-4;i++)
    {
    sprintf(txt,"%s,%s -> %s %s ", vSigmaCh[i].prtcl[0],
                                   vSigmaCh[i].prtcl[1],
                                   vSigmaCh[i].prtcl[2],
                                   vSigmaCh[i].prtcl[3]);
    if(vSigmaCh[i].prtcl[4]) strcat(txt,vSigmaCh[i].prtcl[4]);                               
    printf("  %-30.30s  %.2E\n",txt,vSigmaCh[i].weight);

//sss+=vSigmaCh[i].weight;    
     
    }
//    printf("chanStat: 11 -> %.2E  12-> %.2E 22->%.2E\n",    2.9979E-26*chStat[0],2.9979E-26*chStat[1],2.9979E-26*chStat[2]);
  
//printf("sss=%E\n", sss);    
    
  }     

                                                                     
  if(SpectraFlag==0 &&  Mcdm0 < 2) printf("WARNING! Spectra obtained at Mcdm=2GeV are used !\n");
  if(SpectraFlag!=0 &&  Mcdm0 < 5) printf("WARNING! Spectra obtained at Mcdm=5GeV are used !\n");	 

  if(errMess) { printf(" %s\n", errMess); free(errMess); errMess=NULL;} 

  return vcs; 
}

static double ME(double m0,double m1,double m2)
{ 
  return  m2*m2*m1*m1 +(m0*m0 -m2*m2 -m1*m1)*(m0*m0 -m2*m2 -m1*m1)/8;
}


double calcSpectrumPlus(char*proc22,int outP,double *spect, int*Err)
{  
   int in,out,err;
   char p[4][P_NAME_SIZE];
   for(int i=0;i<NZ;i++) spect[i]=0;
      
   err=proc2names(proc22, &in, &out, p[0],p[1],p[2],p[3],NULL);
   if(err) { if(Err) *Err=1; return 0;}

   int d1,d2;
   char*ap1=antiParticle(p[0]);
   for(d1=1;d1<=Ncdm;d1++) if(strcmp(CDM[d1],p[0])==0 || strcmp(CDM[d1],ap1)==0) break;
   char*ap2=antiParticle(p[2]);
   for(d2=1;d2<=Ncdm;d2++) if(strcmp(CDM[d2],p[0])==0 || strcmp(CDM[d2],ap2)==0) break; 
   if(d1>Ncdm || d2>Ncdm) {if(Err) *Err=2; return 0;}
   double C; // combinatoric coefficient
   { double n=0;
     for(int i=1;i<=Ncdm;i++) n+=fracCDM[i]/McdmN[i];
//printf("n=%e fracCDM[1]=%E \n",n,fracCDM[1]);     
     C=2*fracCDM[d1]*fracCDM[d2]/(McdmN[d1]*McdmN[d2]*n*n);   
     if(strcmp(p[0],ap1)) C/=2;
     if(strcmp(p[1],ap2)) C/=2;
     if(strcmp(p[0],p[1])==0) C/=2;
     if(strcmp(p[0],ap2) && (strcmp(p[0],ap1) || strcmp(p[1],ap2))) C*=2; 
   }  
     
   numout*cc22=newProcess(proc22);
   if(!cc22) { if(Err) *Err=3; return 0;}
  
   double m1=pMass(p[2]), m2=pMass(p[3]), w1=pWidth(p[2],NULL), w2=pWidth(p[3],NULL); 
   if(w1==0 &&  w2==0) { if(Err) *Err=4; return 0;}   

   double M1=pMass(p[0]), M2=pMass(p[1]);
   #define V0 (vRot/299792.*1.732)  /*  1.732=sqrt(3)  */
   double Pcm=V0*M1*M2/(M1+M2);
   double m0=sqrt(M1*M1+Pcm*Pcm) + sqrt(M2*M2+Pcm*Pcm);   
   spect[0]=m0/2;
   
   double const X7[7]={2.544604E-02,1.292344E-01 ,2.970774E-01 ,5.000000E-01 ,7.029226E-01 ,8.707656E-01 ,9.745540E-01 };
   double const F7[7]={6.474248E-02,1.398527E-01 ,1.909150E-01 ,2.089796E-01 ,1.909150E-01 ,1.398527E-01 ,6.474248E-02 };

   char *q1=p[2], *q2=p[3];
   if(strcmp(q1,q2)==0 || strcmp(q1,antiParticle(q2))==0) q2=q1;
   int n1=2;
   double sum=0;
   for(int k=0;k<2;k++)
   {  if(k) 
      { if(q1==q2){  sum*=2; break;}
        else { double m=m1; m1=m2; m2=m;  
                      m=w1; w1=w2; w2=m;
               char*q=q1;q1=q2;q2=q1;       
             }
      }                

      double y_11=atan(-m1/w1), y_12=atan( (pow(m0/(1+m2/m1),2)-m1*m1)/(m1*w1));        //  m1_^2 =  m1^2 + w1*m1*tan(y1)
      double y_1m=atan((pow(m0-m2-w2,2)-m1*m1)/(m1*w1));
      if(y_1m<y_11) y_1m=y_11; 
//double X[1000];   
      for(int l=0;l<2;l++) 
      { double ya,yb;
        if(l==0)
        { if(y_1m==y_11) continue; else { ya=y_11; yb=y_1m;} }
        else { ya=y_1m; yb=y_12;}
      
        for(int i=0; i<7;i++) 
        {  
//X[i]=0;   
          double x1=X7[i];
//       double x1=(i+0.5)/100.;
          double y1=ya+  pow(x1,n1)*(yb-ya);
          double m1_=m1*m1 + w1*m1*tan(y1);
          if(m1_<0) m1_=0;
          m1_=sqrt(m1_);
          double J1=F7[i] * n1*pow(x1,n1-1)*(yb-ya);
          double y_21= atan( ( pow(m1_*m2/m1,2)-m2*m2)/(m2*w2)) , y_22=atan( ( pow(m0-m1_,2)-m2*m2)/(m2*w2));   // m2_^2 = m2^2 + w2*m2*tan(y2)
          int n2=2;

          for(int j=0;j<7;j++)
          {  
            double outMass[]= {0.,0.511E-3,0.939,0.,0.,0.};
//          double x2=(j+0.5)/1000.;
            double x2=X7[j]; 
            double y2=y_21+ pow(x2,n2)*(y_22-y_21);
            double J2=F7[j]*n2*pow(x2,n2-1)*(y_22-y_21);
            double m2_=m2*m2 + w2*m2*tan(y2);
            if(m2_<=0) continue;
            m2_=sqrt(m2_);
            double pcm=decayPcm(m0,m1_,m2_);
            double tab1[NZ], tab2[NZ]; 
            decaySpectrum(p[2],m1_,outP, tab1);
            decaySpectrum(p[3],m2_,outP, tab2);
            double Y=acosh(sqrt(m1_*m1_ + pcm*pcm)/m1_);
            boost(Y,m0/2,outMass[outP],tab1);
                   Y=acosh(sqrt(m2_*m2_ + pcm*pcm)/m2_);
            boost(Y,m0/2,outMass[outP],tab2);
            double f=J1*J2*ME(m0,m1_,m2_)*pcm;
            sum+=f;
            for(int i=1;i<NZ;i++) spect[i]+=f*(tab1[i]+tab2[i]);
//          X[i]+=J1*J2*ME(m0,m1_,m2_)*decayPcm(m0,m1_,m2_);
          } 
        }   
//   displayPlot("small m1 ","x",0,1,0,1,"X",100,X,NULL);
      }
   }
   for(int i=1;i<NZ;i++) spect[i]/=sum;
   
//   calculation vSigma 

   txtList L3, L4; 
   double  w3=pWidth(p[2],&L3), w4=pWidth(p[3],&L4),b3=0,b4=0; 
   if(w3==0 || w4==0) return 0;
   
   double  M3,M4;

   char proc24[100];      
   sprintf(proc24,"%s,%s->", p[0],p[1]);

   char  q3[4][P_NAME_SIZE],q4[4][P_NAME_SIZE];
   
   char*ch=proc24+strlen(proc24);   
   for(txtList L=L3;L;L=L->next)
   {  double b,ms;
      int n=sscanf(L->txt,"%lf %s -> %[^,], %[^,], %[^,]",&b,q3[0],q3[1],q3[2],q3[3]);
      ms=pMass(q3[1])+pMass(q3[2]);
      if(b>b3 && n==4 && ms < m0 )
      {
         sprintf(ch,"%s,%s", q3[1],q3[2]); 
         b3=b; 
         M3=ms;     
      } 
   }
   if(b3==0) { if(Err) *Err=5; return NAN;} 
   
   ch=proc24+strlen(proc24);
   for(txtList L=L4;L;L=L->next)
   {  double b,ms;
      int n=sscanf(L->txt,"%lf %s -> %[^,], %[^,], %[^,]",&b,q4[0],q4[1],q4[2],q4[3]);
      ms=pMass(q4[1])+pMass(q4[2]); 
      if(b>b4 && n==4 && ms<m0)
      {
        if(strcmp(q3[1],q4[1]) && strcmp(q3[2],q4[2]) && strcmp(q3[1],q4[2])&& strcmp(q3[2],q4[1]) ) 
        { sprintf(ch,",%s,%s",q4[1],q4[2]); 
          b4=b; 
          M4=m0;
        }   
      } 
   }
   if(b4==0) { if(Err) *Err=6; return NAN;} 
  
//printf("         proc=%s\n", proc24);
   
   char lib[100];
   err=process2Lib(proc24,lib);
 
   char excludeV[2*P_NAME_SIZE+20];
   if( abs(pNum(p[2]))==abs(pNum(p[3]))) 
   {  sprintf(excludeV,"%s<2", p[2]); 
      sprintf(lib+strlen(lib),"El2%d", abs(pTabPos(p[2])));
   }  else 
   {  sprintf(excludeV,"%s<1,%s<1", p[2],p[3]); 
      sprintf(lib+strlen(lib),"El1%dEl1d",abs(pTabPos(p[2])),abs(pTabPos(p[3])));
   }  
   
   numout*cc24=getMEcode(0,ForceUG, proc24, excludeV, NULL,lib);
//   printf("cc24=%p proc23=|%s| excludeV=|%s| lib=|%s|\n",cc24, proc24, excludeV, lib);

   if(cc24==NULL) return 0;
   passParameters(cc24);
 
   int i3=2,i4=3;
   for(i3=2;     strcmp(q3[1], cc24->interface->pinf(1,i3+1,NULL,NULL)); i3++) /*printf("%d %s %s\n", i3, q3[1], cc24->interface->pinf(1,i3+1,NULL,NULL))*/  ;
   for(i4=i3+1;  strcmp(q3[2], cc24->interface->pinf(1,i4+1,NULL,NULL)); i4++);

//printf("i3=%d i4=%d q3[1]=%s q3[2]=%s  \n", i3,i4,q3[1],q3[2]);
   
//   printf(" %s %s\n", cc24->interface->pinf(1,i3+1,NULL,NULL), cc24->interface->pinf(1,i4+1,NULL,NULL) );
   double vcs, dcs,chi2;
   int N;            
   vcs= V0*cs24Vegas(cc24, 1, Pcm, i3, i4, 
       2 , 100000,  2, 100000, 
      &dcs, &chi2, &N)/b3/b4;
      printf(" Pcm=%.2E   vcs=%.2E dcs=%.2E chi2=%.2E err=%d  C=%e\n", Pcm, vcs, dcs, chi2, err,C);   
   if(strcmp(p[2],p[3])==0)   vcs/=2;  

/*
   passParameters(cc22);  
   double Ms=pMass(p[2])+pMass(p[3]);
*/   
   return vcs*C*2.9979E-26;   // cm^3/s ??     
}




double SpectdNdE(double E, double *tab){
  double z;
  if(E>tab[0]) return 0;
  z=log(E/tab[0]); 
  return 1/E*zInterp(z,tab); 
}

double eSpectdNdE(double E, double *tab){
  double z;
  if(E>tab[0]) return 0;
  z=log(E/tab[0]); 
  return  zInterp(z,tab); 
}


double FSRdNdE(double E, double p,double m, double q, int spin2)
{
    if(E>p || E<=0) return 0;
    double x=E/p;
    
    double kappa=sqrt(1/(1+m*m/p/p));
    double f=q*q/M_PI/137.*log((1+kappa)/(1-kappa))/E;
    if(spin2&1) return f*(1-x+x*x/2); else return f*(1-x);     
}

// ======================================  CDM spectrum limit by Planck  ============================

static double  fe(double E) //  Fig.3   1506.03811
{
  double ln10E[40]= { 3.6990, 3.9315, 4.1640, 4.3965, 4.6291, 4.8616, 5.0941, 5.3267, 5.5592, 5.7917, 6.0242, 6.2568, 6.4893, 6.7218, 6.9543, 7.1869, 7.4194, 7.6519, 7.8844, 8.1170, 8.3495, 8.5820, 8.8145, 9.0471, 9.2796, 9.5121, 9.7446, 9.9772, 10.2097, 10.4422, 10.6747, 10.9073, 11.1398, 11.3723, 11.6048, 11.8374, 12.0699, 12.3024, 12.5349, 12.7675};
  double  e[40]=    { 0.3084, 0.3072, 0.3044, 0.2983, 0.2944, 0.2819, 0.2593, 0.2283, 0.1904, 0.1527, 0.1770, 0.3680, 0.6173, 0.8228, 0.9298, 0.9756, 0.9932, 0.9889, 0.9288, 0.8256, 0.7032, 0.5589, 0.4296, 0.3952, 0.4660, 0.5714, 0.6388, 0.6243,  0.5593,  0.4868,  0.4502,  0.4357,  0.4083,  0.3994,  0.4049,  0.4064,  0.4070,  0.4016,  0.4036,  0.4056};
  return polint3(log10(E), 40,ln10E,e);
}

static double fA(double E)  //  Fig.3  1506.03811
{
  double ln10E[40]= { 3.6990, 3.9315, 4.1640, 4.3965, 4.6291, 4.8616, 5.0941, 5.3267, 5.5592, 5.7917, 6.0242, 6.2568, 6.4893, 6.7218, 6.9543, 7.1869, 7.4194, 7.6519, 7.8844, 8.1170, 8.3495, 8.5820, 8.8145, 9.0471, 9.2796, 9.5121, 9.7446, 9.9772, 10.2097, 10.4422, 10.6747, 10.9073, 11.1398, 11.3723, 11.6048, 11.8374, 12.0699, 12.3024, 12.5349, 12.7675};
  double A[40]=     { 0.8767, 0.7697, 0.7021, 0.6842, 0.6963, 0.6748, 0.5851, 0.4591, 0.3256, 0.2161, 0.1448, 0.1631, 0.2924, 0.4681, 0.5957, 0.6781, 0.7180, 0.7260, 0.7195, 0.6858, 0.6288, 0.5462, 0.4469, 0.3591, 0.3184, 0.3378, 0.3905, 0.4228,  0.4212,  0.3902,  0.3307,  0.3901,  0.4332,  0.4148,  0.4029,  0.4047,  0.4069,  0.4017,  0.4037,  0.4056};
  return polint3(log10(E), 40,ln10E,A);
}

static double intEe(double E, double *tab) { return eSpectdNdE(E,tab)*fe(E*1E9);}
static double intEA(double E, double *tab) { return eSpectdNdE(E,tab)*fA(E*1E9);} 


double   PlanckCMB(double vSigma, double *Sg,double*Se)  // 1506.03811
{ int err;
  double Egmin=5E-6, Eemin=5E-6; 
  if(Eemin<1E-7*Se[0]) Eemin=1E-7*Se[0];
  double r=simpson_arg(intEA,Sg, Egmin,Sg[0],1E-3,NULL)/2 + simpson_arg(intEe,Se,Eemin,Se[0], 1E-3,NULL);
  double nEff=0;
  for(int i=1;i<=Ncdm;i++) nEff+=fracCDM[i]/McdmN[i];  
  r*=vSigma*nEff*nEff;
  return r/3.2E-28;  
}
