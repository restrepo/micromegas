//#define NOSTATISTICS
//#define ONSHELL
//#define NoWeff

#define e0 1.E-6

#define nGAUSS 4

//#define ERROR_PLOT

#include "micromegas.h"
#include "micromegas_aux.h"
#include "micromegas_f.h"
#include "../CalcHEP_src/c_source/ntools/include/vegas.h"

static int tGridDim=0;
static double *ltGrid=NULL,*yfi=NULL;
double YFi(double T)
{  
   if(tGridDim==0) return NAN;
   double lnT=log(T);      
   
   if(lnT<ltGrid[0] || lnT>ltGrid[tGridDim-1]) return NAN;

   if(lnT>ltGrid[tGridDim-2]) return (exp(ltGrid[tGridDim-1])-T)/(exp(ltGrid[tGridDim-1]) - exp(ltGrid[tGridDim-2]))*yfi[tGridDim-2];  

   return polint3(lnT,tGridDim-1,ltGrid,yfi);  
} 

void printYFi(void)
{ for(int i=0;i<tGridDim;i++) printf("T=%E Y=%E\n",exp(ltGrid[i]),yfi[i]); 
}

static int makeLnTgrid(double TR)
{
//printf("make grid\n");
  if(Tend>=TR) return 1;
  if(Tend==0) return 1;
  tGridDim=5*log(TR/Tend)/log(10)+1;
  if(tGridDim<10) tGridDim=10;
  ltGrid=realloc(ltGrid, tGridDim*sizeof(double));
  yfi=realloc(yfi, tGridDim*sizeof(double));
  for(int i=0;i<tGridDim;i++) { ltGrid[i]= log( Tend*pow(TR/Tend,(double)i/(tGridDim-1.))); yfi[i]=0;}
  return 0;
}  
   
char*NAME=NULL,*aNAME=NULL;

int nFeeble=0;
char ** FeebleP=NULL;
double * FeebleY=NULL;

REAL * Taddress=NULL;

static double tTcut=0;
double mTcut=0;


//=====================  SECTION:  Function for Statistics.    Stat2,  K1. K2, gammaLor

typedef struct { double m1,m2,m3, s1,s2,s3, pcm, e2,e3,e2_,e3_,gamma;}  K1to2_struct;


static double K1to2_int(double x, void * par_)  
{
//  NTOT++;
  double eps=1E-8;
  if(x==0 || x==1) return 0;
  K1to2_struct * par=par_;
  double pcm=par->pcm,s1=par->s1,s2=par->s2,s3=par->s3,e2=par->e2,
         e3=par->e3,e2_=par->e2_,e3_=par->e3_,m1=par->m1;

  double n=par->gamma;

  double J=exp(-m1)*(n+1)*(n+2)*pow(1-x,n)*x;
  double z=pow(1-x,n+1)*(1+x*(n+1));

  double eE= exp(-m1)*(1- z);
  double Es=m1;
  if(z<eps) Es+=z; else Es-=log(1-z); 
  
  double chY,shY,dchY;
  
  if(z<eps)  dchY=z*(1+z/2)/m1; else dchY= -log(1-z)/m1;
  chY=1+dchY; 
  shY=sqrt(dchY*(1+chY));

  if(shY==0) return 0;
  double D1,D11;

  if(Es > eps)
  {   D1= (1+s1*eE);
      D11=(1-s2*s3*eE);
  } else 
  {  D1= (1+s1)-s1*Es*(1-Es);     
     D11=(1-s2*s3) + s2*s3*Es*(1-Es);
  }   
 
  if(s2==0 && s3==0) return J*2*pcm*shY /D1/D11;
  
  double intS;
  long double D2p,D3p,D2m,D3m,E2p,E3p,E2m,E3m;

  E2m=chY*e2_+pcm/(chY+shY);
  E3m=chY*e3_+pcm/(chY+shY);

  E2p=E2m+2*shY*pcm;
  E3p=E3m+2*shY*pcm;

  if(E2m> eps) D2m=1+s2*expl(-E2m); else  D2m=(1+s2)-s2*E2m*(1-E2m); 
  if(E3m> eps) D3m=1+s3*expl(-E3m); else  D3m=(1+s3)-s3*E3m*(1-E3m);
  if(E2p>eps)  D2p=1+s2*expl(-E2p); else  D2p=(1+s2)-s2*E2p*(1-E2p);
  if(E3p>eps)  D3p=1+s3*expl(-E3p); else  D3p=(1+s3)-s3*E3p*(1-E3p);

  intS= ( 2*pcm*shY  +logl(D2p*D3p/D2m/D3m) )/D1/D11;

  return  J*intS;
}


double K1to2(double m1,double m2,double m3, double s1,double s2,double s3)
{
#ifdef NOSTATISTICS
    return bessK1(m1);
#endif
// !! for photons s=1 
    s1*=-1; s2*=-1;s3*=-1;    
// !! for photons s=-1     
if(s1<-exp(m1)) s1=-exp(m1)*0.9999;  
if(s2<-exp(m2)) s2=-exp(m2)*0.9999;
if(s3<-exp(m3)) s3=-exp(m3)*0.9999;

    if(m1<=m2+m3) return 0; 
    double pcm=decayPcm(m1,m2,m3); 
    if(m1>25 ) return bessK1(m1);
    double delta1,delta23;
    K1to2_struct  par;
    par.m1=m1;par.m2=m2;par.m3=m3; par.s1=s1;par.s2=s2;par.s3=s3;
    par.pcm=pcm;
    par.e2=sqrt(m2*m2+pcm*pcm);
    par.e3=sqrt(m3*m3+pcm*pcm);
    par.e2_=m2*m2/(par.e2+pcm);
    par.e3_=m3*m3/(par.e3+pcm);

    if(s1>=0) delta1=1; 
    else if(m1 > 1E-6) delta1 = -1/s1-exp(-m1); 
    else delta1=-(1+s1)/s1+m1;
    
    if(s2*s3<=0 || s2>0) delta23=1;
    else if(m1 > 1E-6) delta23 = 1/s2/s3-exp(-m1);
    else delta23= (1-s2*s3)/s2/s3+m1;

    double gamma23; 
    if(delta23>0.1) gamma23=0;
    else
    { gamma23=1;
      for(int i=0;i<5;i++) gamma23=log(gamma23*delta23)/log(0.05);
    }      
    
    double gamma123;
    if(delta1>0.1 || delta23>0.1) gamma123=0;
    else
    {   
       double delta123=delta1>delta23? delta1:delta23;
       gamma123=log(delta123)/log(0.5);
    } 
    par.gamma= gamma123>gamma23 ? gamma123 : gamma23;
    if(par.gamma<1) par.gamma=1;   

    int err;
    double res=  0.5*simpson_arg(K1to2_int,&par,0,1,1E-3,&err)/pcm;
    if(err) 
    {  printf("Problem in K1to2(%e, %e, %e, %E, %E, %E) err=%d File  freezein.c, line 134)\n", m1,m2,m3,-s1,-s2,-s3,err);
#ifdef ERROR_PLOT
    displayPlot("K1to2_int", "x", 0,1, 0,1,"K1to2_int", 0,K1to2_int,&par);  
    exit(0);
#endif   
    }    
    return res;
}


double Stat2(double P, double M, double m1,double m2,double eta1, double eta2)
{

#ifdef NOSTATISTICS
   return 1;
#else   
   if(eta1==0 && eta2==0) return 1;
   double E=sqrt(P*P+M*M);
   if(E>25) return 1;
    
   double pcm=decayPcm(M,m1,m2);
   double e1=sqrt(m1*m1+pcm*pcm),e2=sqrt(m2*m2+pcm*pcm);
   double sht=P/M, cht=E/M;
   double ex_1=eta1*exp(-cht*e1), ex_2=eta2*exp(-cht*e2), ex_p=exp(-sht*pcm);
   
   double f;
   if(pcm*sht < 1E-5*(m1+m2)*cht) f=1+ex_1/(1-ex_1)+ ex_2/(1-ex_2); else
   f= 1+log( (1-ex_1*ex_p)*(1-ex_2*ex_p)/(1-ex_1/ex_p)/(1-ex_2/ex_p))/(2*pcm*sht);
            
   return f/(1-eta1*eta2*exp(-E));     
#endif   
}


static double K2_int(double y, void *arg)
{  double *arg_=arg;
   double mu=arg_[0];
   double s=arg_[1];
   double n=arg_[2];
   if(y==0) return 0; 
   double x=log(y);
   double J=1/y;
   return J*x*x/(exp(sqrt(x*x+mu*mu)) - s)*pow(mu/sqrt(x*x+mu*mu),n) ;
}

static double K2prime(double mu, double s)
{ double arg[3];
#ifdef NOSTATISTICS
  return bessK2(mu);  
#endif
  if(mu+1<s) s=mu+0.9999;
  if(-s>exp(mu)*0.999) s=-exp(mu)*0.999;   
  arg[0]=mu,arg[1]=s; arg[2]=0;
  int err;
  double res=simpson_arg(K2_int,arg, 0, 1,1E-5,&err)/mu/mu;
  if(err) 
  { printf("Warning from Simpson (freezein.c line 193) code=%d. Precison is not reached\n",err);
#ifdef ERROR_PLOT
    displayPlot("Error in simpson", "x", 0,1,0,1,"K2_int",0,K2_int,arg);
    exit(0);      
#endif
  }
  return res;  
}  

static double gammaLor(double mu, double s)
{ 
#ifdef NOSTATISTICS
   s=0;
#endif     
  if(fabs(s)< exp(mu)*1E-3) s=0;

  if( (s <= 1 && s >= -1)|| mu>20 ) 
  { double t=1/mu;
    if(s>1)   s=1;    
    if(s<=-1) s=-1;  
    double fb[21]={1 ,1.1927E+00 ,1.7945E+00 ,2.7716E+00 ,4.1010E+00 ,5.7833E+00 ,7.8218E+00 ,1.0219E+01 ,1.2977E+01 ,1.6097E+01 ,1.9579E+01 ,2.3425E+01 ,2.7635E+01 ,3.2209E+01 ,3.7147E+01 ,4.2450E+01 ,4.8118E+01 ,5.4150E+01 ,6.0548E+01 ,6.7310E+01 ,7.4438E+01};
    double fm[21]={1 ,1.1927E+00 ,1.8143E+00 ,2.9269E+00 ,4.5580E+00 ,6.7126E+00 ,9.3876E+00 ,1.2579E+01 ,1.6282E+01 ,2.0494E+01 ,2.5214E+01 ,3.0439E+01 ,3.6168E+01 ,4.2400E+01 ,4.9136E+01 ,5.6373E+01 ,6.4112E+01 ,7.2352E+01 ,8.1094E+01 ,9.0337E+01 ,1.0008E+02};
    double ff[21]={1 ,1.1927E+00 ,1.8322E+00 ,3.0401E+00 ,4.8384E+00 ,7.2169E+00 ,1.0166E+01 ,1.3679E+01 ,1.7751E+01 ,2.2380E+01 ,2.7563E+01 ,3.3298E+01 ,3.9585E+01 ,4.6422E+01 ,5.3810E+01 ,6.1748E+01 ,7.0235E+01 ,7.9271E+01 ,8.8856E+01 ,9.8990E+01 ,1.0967E+02};
    double x0[21]={0 ,1.2500E-01 ,5.0000E-01 ,1.1250E+00 ,2.0000E+00 ,3.1250E+00 ,4.5000E+00 ,6.1250E+00 ,8.0000E+00 ,1.0125E+01 ,1.2500E+01 ,1.5125E+01 ,1.8000E+01 ,2.1125E+01 ,2.4500E+01 ,2.8125E+01 ,3.2000E+01 ,3.6125E+01 ,4.0500E+01 ,4.5125E+01 ,5.0000E+01};

    double be[3];
    if(t>50) { be[0]=fb[20]+1.461824*(t-50); be[1]=fm[20]+2*(t-50);     be[2]=ff[20]+2.191793*(t-50);}
    else     { be[0]=polint3(t,21,x0,fb);    be[1]=polint3(t,21,x0,fm); be[2]=polint3(t,21,x0,ff);}

    return  0.5*(1+s)*s*be[0]+(1-s)*(1+s)*be[1]+ 0.5*s*(s-1)*be[2]; 
  }    
  return  K2prime(mu,s)/K1to2(mu,0,0,s,0,0);
}


//==============================   SECTION DECAY =============================

static double decayIntegrand(double lnT, void *arg_)
{   
    double s1=((double*)arg_)[0];
    double M= ((double*)arg_)[1];
    double T=exp(lnT);
    double J=T;
    double E=2*M_PI*M_PI/45*hEff(T)*T*T*T;
    double H=sqrt(8*M_PI/3.*M_PI*M_PI/30.*gEff(T))*T*T/MPlanck;  //Hubble 
  
    double dtdT=(1+hEffLnDiff(T)/3)/(H*T);              // dt => dT
    return  K1to2(M/T,0,0,s1,0,0)*T/E*dtdT*J;  
}

double decayAbundance(double TR, double M, double w, int Ndf,double eta, int plot)
{ 
    double Tmin=M/100;
    double arg[2]={eta,M};     
    int err;
    
    if(plot)
    {  double arr[100];
       for(int i=0;i<100;i++) 
       { double T=Tmin*pow(TR/Tmin,(i+0.5)/100); 
         arr[i]=decayIntegrand(log(T),&arg)*w*Ndf*fabs(eta)*M*M/(M_PI*M_PI)*log(10);
       }
       char txt[50];
       sprintf(txt," decay Abundance. Mass=%.2E Width=%.2E\n",M,w); 
       displayPlot(txt,"T",Tmin,TR,1,1,"dY/dlog10(T)",100,arr,NULL); 
    }
    double res= w*Ndf*fabs(eta)*M*M/(M_PI*M_PI)*simpson_arg( decayIntegrand,&arg,log(Tmin),log(TR),0.001,&err);
    if(err)
    { printf("Warning from Simpson. Precison is not reached (freezein.c line 287 code=%d)\n",err);
#ifdef ERROR_PLOT
      displayPlot("Error in simpson", "x", log(Tmin), log(TR),0,1,"decayIntegrand",0,decayIntegrand,&arg);
      exit(0);      
#endif
    }
    return res;
}


static double Yeq_1(double T, double M ,double eta)  // M=0 is not supported 
{ 
  double mu=M/T;   
  if(eta && log(fabs(eta)*1E4)< mu) eta=0;
  if(mu<25 && fabs(eta)>0 ) return mu*mu/((2*M_PI*M_PI)*(2*M_PI*M_PI/45)*hEff(T))*K2prime(mu,eta);
  else return mu*mu/((2*M_PI*M_PI)*(2*M_PI*M_PI/45)*hEff(T))*K2pol(1/mu)*exp(-mu)*sqrt(M_PI/2/mu);
}


static double getEta(double Y, double T, double M, int s)
{ 
  if(Y==0) return 0;
  double mu=M/T;
  double eta=Y/Yeq_1(T,M,0); 
  if(eta < exp(mu)*1E-3) return s*eta;
  for(int i=0;i<3;i++) { if(s*eta>exp(mu)*0.999 ) eta=exp(mu)*0.999; eta=Y/Yeq_1(T,M,s*eta);}
  if(s*eta>exp(mu)) 
  { 
//  printf("getEta:  unphysical  eta=%E : Y=%E  T=%E  exp(mu)=%E \n",eta,Y, T,exp(mu)); 
    eta=exp(mu)*0.999;}    
  return s*eta;
}    


static double  medMass, medWidth;
static txtList medDecayList;  
static int     medEta, medNDF;
static int     etaOn;            
static double _p_,T0;
aChannel*omegaFiCh=NULL;

static int getDecBrKE(double T, txtList L, double eta, double * brTot, double * brBath,  double  *brSig1,int * nFiCh)
{ 
  *brBath=0;
  *brTot=0;
  *brSig1=0;
  int nFiCh_;
  if(nFiCh) nFiCh_=*nFiCh;
  int s[5];
  double mu[5]; 


  for(txtList l=L; l;l=l->next)
  {
    char name[5][20];
    double br;
    int n=sscanf(l->txt,"%lf %s -> %[^,], %[^,], %[^,], %[^,]",&br,name[0],name[1],name[2],name[3],name[4]);
//printf("l->txt=%s\n", l->txt);    
    n-=2; // number of outgoing particles
    for(int i=0;i<=n;i++)mu[i]=pMass(name[i])/T;
    for(int i=1;i<=n;i++)if(isFeeble(name[i])) s[i]=0; else
    {  
       qNumbers(name[i], s+i ,NULL,NULL);
       if(s[i]&1) s[i]=-1; else s[i]=1;
    }

//     if(mu[0]<25)  printf("mu=%E s[1]=%d,s[2]=%d eta=%e K=%E\n",mu[0],s[1],s[2], eta, K1to2(mu[0],mu[1],mu[2],eta,s[1],s[2])/K1to2(mu[0],0,0,eta,0,0));    
    if(mu[0]<25)
    {  double mu_[2]={0,0};
       int s_[2]={0,0};
       int ii=0;
       for(int i=1;i<=n;i++) if(s[i]) {s_[ii]=s[i]; mu_[ii]=mu[i];ii++; if(ii>1) break;} 
       if(ii) br*=K1to2(mu[0],mu_[0],mu_[1],eta,s_[0],s_[1])/K1to2(mu[0],0,0,eta,0,0);
    } 
     *brTot+=br;
     {
       int ba=1;
       for(int i=1;i<=n;i++) if(s[i]==0)ba=0;
       if(ba) *brBath+=br;
     }  
     {  double brS1=0;
        for(int i=1; i<=n;i++) if(s[i]==0 && (strcmp(name[i],NAME)==0||strcmp(name[i],aNAME)==0)) brS1+=br; 
        *brSig1+=brS1;

        if(brS1>0 && nFiCh)
        { omegaFiCh[nFiCh_].weight=brS1*pMass(NAME);
          omegaFiCh[nFiCh_].err=0;
          omegaFiCh[nFiCh_].nin=1;
          for(int i=0;i<=n;i++) 
          {  int k=pTabPos(name[i]);
             if(k>0)  omegaFiCh[nFiCh_].prtcl[i]= ModelPrtcls[k-1].name; 
              else    omegaFiCh[nFiCh_].prtcl[i]= ModelPrtcls[-k-1].aname;
          }
          for(int i=n+1;i<5;i++) omegaFiCh[nFiCh_].prtcl[i]=NULL;
          nFiCh_++; 
          omegaFiCh=realloc(omegaFiCh,(nFiCh_+1)*sizeof(aChannel));
        }      
     }          
  }
  
  if(nFiCh) 
  { double s=0;
    for(int i=*nFiCh;i<nFiCh_;i++) s+=omegaFiCh[i].weight;
    for(int i=*nFiCh;i<nFiCh_;i++) omegaFiCh[i].weight/=s;
    *nFiCh=nFiCh_;
  }
  return 0;
}


static  void  mediatorDerivs( double lnT, double *Y,double * dY)
{ 
   double T=exp(lnT);
   double eta;
  if(etaOn) eta=getEta(*Y,T,medMass,medEta); else eta=medEta;

   double brBath,brTot,brSig1;
   getDecBrKE(T, medDecayList, eta,  &brTot, &brBath, &brSig1,NULL);     

   double H=sqrt(8*M_PI/3.*M_PI*M_PI/30.*gEff(T))*T*T/MPlanck;  //Huble 
   double dtdLnT=(1+hEffLnDiff(T)/3)/H;                             // dt => dlnT

   double alpha=medWidth*brTot/gammaLor(medMass/T,eta)*dtdLnT;
   double beta=brBath/brTot;
   if(Y[0]<0) Y[0]=0;
   dY[0]= (alpha*(Y[0]-beta*Yeq_1(T,medMass,eta)));
   double c=medWidth*Y[0]/gammaLor(medMass/T,eta)*dtdLnT;
   dY[1]= -c*brSig1;
}    

static double derEq0(double lnT)
{
   double T=exp(lnT);
   double brBath,brTot,brSig1;
   getDecBrKE(T,medDecayList,medEta,&brTot,&brBath,&brSig1,NULL);     

   double H=sqrt(8*M_PI/3.*M_PI*M_PI/30.*gEff(T))*T*T/MPlanck;  //Huble 
   double dtdLnT=(1+hEffLnDiff(T)/3)/H;                       // dt => dlnT

   double alpha=medWidth*brTot/gammaLor(medMass/T,medEta)*dtdLnT;
//   double beta=brBath/brTot;
  
   double  Y=Yeq_1(T,medMass,medEta);
   
   return medWidth*brSig1*Y/gammaLor(medMass/T,medEta)*dtdLnT; 
}

static int saveMediators=0;

double darkOmegaFiDecay(double TR, char * pname, char* fiPrtcl)
{ 
  int plot=0; 
  int KE=1;
  int p=pTabPos(fiPrtcl);
  if(p==0) { printf("darkOmegaFiDecay: unknown particle %s\n", fiPrtcl);  return 0; }
  if(! isFeeble(fiPrtcl)) { printf("darkOmegaFiDecay: particle %s is not feeble\n", fiPrtcl);  return 0; } 
  p=abs(p)-1;
  NAME=ModelPrtcls[p].name; aNAME=ModelPrtcls[p].aname;

  if(TR<=Tend) { printf("Reheating temperature is  too small\n"); return 0;} 
  int tabPos=pTabPos(pname);
  if(tabPos==0) { printf(" '%s' unknown particle\n",pname); return 0; }
  tabPos=abs(tabPos)-1;
  medNDF=(ModelPrtcls[tabPos].cdim)*(ModelPrtcls[tabPos].spin2+1);
  if(strcmp(ModelPrtcls[tabPos].name,ModelPrtcls[tabPos].aname)) medNDF*=2; // particle +antiparticle 
  
  medEta=  ModelPrtcls[tabPos].spin2&1? -1:1;

  int VZdecayMem=VZdecay, VWdecayMem=VWdecay;  VZdecay=0;VWdecay=0;
  if(VZdecayMem|| VWdecayMem) cleanDecayTable();
   Taddress=NULL;  
   for(int i=0;i<nModelVars;i++) if(strcmp(varNames[i],"T")==0) { Taddress=varValues+i; break;}
  

  if(Taddress)     
  { *Taddress=0;
    calcMainFunc();
    cleanDecayTable();
  }

  medMass=pMass(pname);
  medWidth=pWidth(pname,&medDecayList);     


  double  *Y1=NULL,*Y2=NULL, *Ym=NULL;
   
  int nStep=tGridDim-1;   
  if(!saveMediators) makeLnTgrid(TR);
  
  int mSize=sizeof(double)*tGridDim;
  if(plot)
  {  Y1=malloc(mSize);
     Y2=malloc(mSize);
     Ym=malloc(mSize);
  }    

  double omega;    
  {
     double eta;
     double T=TR;
     int s0,nT;
     double y[3]={0,0,0};
     
     if(isFeeble(pname)) { s0=0; eta=0;} else { s0=medEta; eta=medEta; y[0]=Yeq_1(TR,medMass,medEta);}
     if(Taddress)
     { *Taddress=TR;
       calcMainFunc();
       cleanDecayTable();
     }

     medWidth=pWidth(pname,&medDecayList);
     medMass=pMass(pname);

     double brTot,brBath,brSig1;

     for(nT=tGridDim-2,T=TR;nT>=0;nT--)
     {  double Tnext=exp(ltGrid[nT]);

        double H=sqrt(8*M_PI/3.*M_PI*M_PI/30.*gEff(Tnext))*Tnext*Tnext/MPlanck;  //Hubble
        double dtdLnT=(1+hEffLnDiff(Tnext)/3)/(H);              // dt => dT

        if(Taddress)
        {  *Taddress=Tnext;
           calcMainFunc();
           cleanDecayTable();
           medWidth=pWidth(pname,&medDecayList);
           medMass=pMass(pname);    
        }

        getDecBrKE(Tnext,medDecayList, eta,  &brTot, &brBath, &brSig1,NULL);
        if(s0)                                                                          
        {  int n;
           eta=s0;                              
                                                   
           y[0]= Yeq_1(Tnext,medMass,s0);                                                  
           y[1]+=gauss(derEq0,log(Tnext),log(T),3);
        }                                                                               
        else                                                                            
        {                                                                               
           double a=(T-Tnext)*medWidth*brTot*dtdLnT/gammaLor(medMass/Tnext,eta);        
           if(a>100) // Y=Yeq up to beta factor;                                        
           {  int n;                                                                          
              double beta=brBath/brTot;                                                 
              eta=beta*medEta;   
              y[0]= beta*Yeq_1(Tnext,medMass,eta);                                         
              y[1]+=gauss(derEq0, log(Tnext),log(T),3);
           }else                                                                        
           {                                                                            
              if(medMass/T<100) etaOn=1; else etaOn=0;                                                                              
              odeint(y,2, log(T), log(Tnext), 1E-3, log(T/Tnext), mediatorDerivs);                                              
              if(etaOn)eta=getEta(y[0],Tnext,medMass,1);                                   
           }                                                                            
        }                                                                               
        yfi[nT]+=y[1]*medNDF;
        T=Tnext;  
      }
      omega=medNDF*(y[1]*pMass(NAME) )*EntropyNow/RhoCrit100;
   }                                              
   if(VZdecayMem|| VWdecayMem)
   {  cleanDecayTable();
      VZdecay=VZdecayMem, VWdecay=VWdecayMem;
   }
   Tstart=TR;
   return omega;
}


// ============ Get Mediator   ===========
static double T;


//static int tGridDim;
//static double*ltGrid=NULL;

static int nMed=0,nMedMax=0;
typedef struct { char *name;
                 int my,massPos,widthPos,spin2,cdim; 
                 double  wDD, wIn, Twidth; 
                 double * wBath;
                 double * wSig;
                 double * wEff;
                 double * mass;
                 double * wIwD;
                 double * eta;
                 double * fKK;
               } s_mediator; 
static s_mediator* mediatorArr=NULL;


static void getGauss(int* nG,double **xG,double **fG)
{
 static  double X1[1]={0.5};
 static  double F1[1]={1};
 static  double X2[2]={2.113249E-01,7.886751E-01 };
 static  double F2[2]={5.000000E-01,5.000000E-01 };
 static  double X3[3]={1.127017E-01,5.000000E-01 ,8.872983E-01 };
 static  double F3[3]={2.777778E-01,4.444444E-01 ,2.777778E-01 };
 static  double X4[4]={6.943185E-02,3.300095E-01 ,6.699905E-01 ,9.305682E-01 };
 static  double F4[4]={1.739274E-01,3.260726E-01 ,3.260726E-01 ,1.739274E-01 };

   if(*nG<1)*nG=1; 
   if(*nG>4)*nG=4;
   switch(*nG)
   { case 1: *xG=X1; *fG=F1;break;
     case 2: *xG=X2; *fG=F2;break;
     case 3: *xG=X3; *fG=F3;break;
     case 4: *xG=X4; *fG=F4;break;
   }  
}            

static void cleanMediatorArr(void)
{ 
   for(int i=0;i<nMed;i++) 
   { 
     if(mediatorArr[i].wBath) free(mediatorArr[i].wBath) ;
     if(mediatorArr[i].wSig ) free(mediatorArr[i].wSig ) ;
     if(mediatorArr[i].wEff ) free(mediatorArr[i].wEff ) ;
     if(mediatorArr[i].mass ) free(mediatorArr[i].mass ) ;
     if(mediatorArr[i].wIwD ) free(mediatorArr[i].wIwD ) ;
     if(mediatorArr[i].eta  ) free(mediatorArr[i].eta  ) ;
     if(mediatorArr[i].fKK  ) free(mediatorArr[i].fKK  ) ;
   }
   free(mediatorArr);
   mediatorArr=NULL;
   nMed=0;  
}



static void fillMediatorArr(double Tr,par_22 *arg)
{ 
  int mSize=tGridDim*sizeof(double);
  int nT;
  double T;

  int *pdg=arg->pdg;
  int i;

  if(nMed)for(i=0;i<nMed;i++) mediatorArr[i].my=0;
  double  step= exp(ltGrid[tGridDim-1]-ltGrid[tGridDim-2]);
//  printf("step=%E\n",step);
  for(int n=1;;n++)
  {  
     int m,w,pnum;
     char*s=arg->cc->interface->den_info(1,n,&m,&w,&pnum);
     if(!s) break;
     if(s[0]==1 && s[1]==2 && m && w)
     {  char*name=ModelPrtcls[pnum].name;

        for(i=0;i<nMed;i++) { if(mediatorArr[i].name==name)  break;}
        if(i<nMed) 
        { mediatorArr[i].widthPos=w;
          mediatorArr[i].my=1;     
          continue; 
        }    
        mediatorArr=realloc(mediatorArr,(nMed+1)*sizeof(s_mediator));
        mediatorArr[nMed].wBath=mediatorArr[nMed].wEff=mediatorArr[nMed].eta=mediatorArr[nMed].mass=mediatorArr[nMed].wIwD=NULL;
        mediatorArr[nMed].fKK=NULL; mediatorArr[nMed].wSig=NULL;
      
        mediatorArr[nMed].widthPos=w;   // position in local 2->2 code
        for(int i=0;i<nModelVars+nModelFunc;i++) if(strcmp(varNames[i],arg->cc->interface->varName[m])==0)  { m=i;  break;} 
        mediatorArr[nMed].massPos=m;    // position in global parameter list
        mediatorArr[nMed].name=ModelPrtcls[pnum].name;
        mediatorArr[nMed].spin2=ModelPrtcls[pnum].spin2;
        mediatorArr[nMed].cdim=ModelPrtcls[pnum].cdim;
        mediatorArr[nMed].wBath =malloc(mSize);
        mediatorArr[nMed].wSig  =malloc(mSize);
        mediatorArr[nMed].wEff  =malloc(mSize);
        mediatorArr[nMed].eta   =malloc(mSize);
        mediatorArr[nMed].mass  =malloc(mSize);
        mediatorArr[nMed].wIwD  =malloc(mSize);
        mediatorArr[nMed].fKK   =malloc(mSize);
        mediatorArr[nMed].my=-1;
        nMed++;  
     }
  } 

  for(nT=tGridDim-1;nT>=0;nT--)   
  {
     T=exp(ltGrid[nT]);
     
     if(Taddress) 
     { *Taddress=T;
        calcMainFunc();
        cleanDecayTable();
     }
       
     for(int n=0;n<nMed;n++) if(mediatorArr[n].my==-1)
     {  
        txtList L;
        double Mm=fabs(varValues[mediatorArr[n].massPos]);
        mediatorArr[n].mass[nT]=Mm;
        double wI=0,wD=0,eta,eta_;
        if(nT==tGridDim-1)
        { if(isFeeble(mediatorArr[n].name)) eta_=0; else
          if(mediatorArr[n].spin2&1)eta_=-1; else eta_=1;
        }
        else eta_= mediatorArr[n].eta[nT+1];
        
        eta=eta_;  
        double w=pWidth(mediatorArr[n].name,&L);
//printf("   T=%E w=%E\n", T,w);        
        double wBath=0,wBath_=0;
        double wSig;
//        for(int k= (nT==tGridDim-1)? 3:1;k;k--)
        {
          if(w>0)for(txtList l=L; l;l=l->next)
          {
             char name[4][11];
             int s[3],pdgD[3],q3[3],cdim[3];
             double mu[3];
             double br;
             if(4 < sscanf(l->txt,"%lf %s -> %[^,], %[^,], %[^,]",&br,name[0],name[1],name[2],name[3])) continue;          
             for(int i=0;i<3;i++) mu[i]=pMass(name[i])/T;
             for(int i=1;i<3;i++)
             { pdgD[i]=qNumbers(name[i], s+i ,q3+i,cdim+i);
              if(isFeeble(name[i])) s[i]=0; else { if(s[i]&1) s[i]=-1; else s[i]=1;}
             }
             int sig=0;
             if((pdgD[1]==pdg[0] && pdgD[2]==pdg[1])|| (pdgD[1]==pdg[1] && pdgD[2]==pdg[0]))  wI=w*br;
             else if((pdgD[1]==pdg[2] && pdgD[2]==pdg[3])|| (pdgD[1]==pdg[3] && pdgD[2]==pdg[2])) { wD=w*br; sig=1;}
             if(mu[0]>mu[1]+mu[2] && ( (s[1]|| s[2]) && mu[0]<25 )) br*=K1to2(mu[0],mu[1],mu[2],eta,s[1],s[2])/K1to2(mu[0],0,0,eta,0,0);
             if(sig) wSig=w*br;
             if(s[1] && s[2]) wBath_+=w*br; 
             wBath+=w*br;
          }
          if(wBath==0) wBath=w;
          if(wBath>0)
          { double mu=Mm/T;
            double alpha=wBath/Hubble(T)/gammaLor(mu ,eta)*log(2),eAlpha;
            double beta=wBath_/wBath;
  //printf("T=%.2E wBath=%.2E H=%.2E massPos=%d    Mm=%.2E    mu=%.2E eta=%E  gammaLor=%.2E  alpha=%.2E \n", T, wBath, Hubble(T), mediatorArr[n].massPos,  Mm, mu, eta,gammaLor(mu ,eta),alpha);        
            if(alpha>20) eAlpha=0; else  eAlpha=exp(-alpha);
                 if(mu<0.1) eta=beta*(1-eAlpha)+fabs(eta_)*eAlpha;
            else if(mu<10)  eta=beta*(1-eAlpha)+fabs(eta_)*eAlpha*bessK2(mu/step)/bessK2(mu)/step/step;
            else         {  double delta=T/Mm -alpha;
                            if(delta<-20) eta=beta*(1-eAlpha);
                                          eta=beta*(1-eAlpha)+fabs(eta_)*exp(delta)*M_SQRT2*K2pol(step/mu)/K2pol(1/mu)/step/step;
                         }
           if(mediatorArr[n].spin2&1)  eta*=-1;
          }
        }
        mediatorArr[n].wSig[nT]=wSig;    
        mediatorArr[n].wBath[nT]=wBath;
        mediatorArr[n].eta[nT]=eta;
        mediatorArr[n].wEff[nT]=wBath;
        
if(!isfinite(eta))
{        
printf("mediatorArr[%d].eta[%d]=%E\n",n,nT,eta);
exit(0);
}        
//printf("T=%E eta=%E wSig=%E wBath=%E \n", T, mediatorArr[n].eta[nT], mediatorArr[n].wSig[nT],mediatorArr[n].wBath[nT]);        
//        mediatorArr[n].my=1;   
     }
   }
   
//return;
   int ng=nGAUSS; 
   double *P,*cG;
   getGauss(&ng,&P,&cG);
          
   T=Tr;
    
  for(nT=tGridDim-1;nT>=0;nT--)   
  {
    T=exp(ltGrid[nT]);
   for(int n=0;n<nMed;n++) if(mediatorArr[n].my)
   { 
      double Mm=mediatorArr[n].mass[nT];

      if(T<Mm/25) mediatorArr[n].fKK[nT]=1; else
      { double mu[3];
         mu[0]=Mm/T;
         mu[1]=pMass(arg->cc->interface->pinf(1,1,NULL,NULL))/T;
         mu[2]=pMass(arg->cc->interface->pinf(1,2,NULL,NULL))/T;
         mediatorArr[n].fKK[nT]= K1to2(mu[0],mu[1],mu[2],0,arg->eta[0],arg->eta[1])/K1to2(mu[0],mu[1],mu[2],mediatorArr[n].eta[nT],arg->eta[0],arg->eta[1]);
      }  
      
      {  double wI=0,wD=0;
         txtList L;
         double w=pWidth(mediatorArr[n].name,&L);
         if(w>0)for(txtList l=L; l;l=l->next)
         {
            char name[4][11];
            int s[3],pdgD[3],q3[3],cdim[3];
            double mu[3];
            double br;
            if(4 < sscanf(l->txt,"%lf %s -> %[^,], %[^,], %[^,]",&br,name[0],name[1],name[2],name[3])) continue;

            for(int i=1;i<3;i++) pdgD[i]=qNumbers(name[i], s+i ,q3+i,cdim+i);
    
            if((pdgD[1]==pdg[0] && pdgD[2]==pdg[1])|| (pdgD[1]==pdg[1] && pdgD[2]==pdg[0]))  wI=w*br;
            else if((pdgD[1]==pdg[2] && pdgD[2]==pdg[3])|| (pdgD[1]==pdg[3] && pdgD[2]==pdg[2])) { wD=w*br;}
         }
         mediatorArr[n].wIwD[nT]=wI*wD; 
     }
     if(mediatorArr[n].wIwD[nT]==0) { mediatorArr[n].wEff[nT]=mediatorArr[n].wBath[nT]; continue;} 
#ifndef NoWeff
     double L[10], z[10];
     int i;
     double zc=1/T/T; 
       
     for(i=0;i<ng;i++) L[i]=-log(P[i]);
     
     double wEff= mediatorArr[n].wBath[nT];
     double eta = mediatorArr[n].eta[nT];
     double I0=0.5*wEff*T*T/Hubble(T)/gammaLor(mediatorArr[n].mass[nT]/T,eta);


     for(i=0;i<ng;i++) z[i]=zc+L[i]/I0;
     for(;;)
     {  double bD=0; 
        for(i=0;i<ng;i++)
        {  double lTi=log(1/sqrt(z[i]));
           if(lTi<ltGrid[0]) bD+=cG[i]*mediatorArr[n].wSig[0]/mediatorArr[n].wBath[0];
           else bD += cG[i]*polint3(lTi,tGridDim,ltGrid,mediatorArr[n].wSig)
                           /polint3(lTi,tGridDim,ltGrid,mediatorArr[n].wBath);
        }                              
        double wD= mediatorArr[n].wSig[nT]/bD;
        
//printf("wD=%E wEff=%E  bD=%E\n", wD, wEff,bD);        
        if(fabs(wD-wEff)<1.E-3*wD) { wEff=wD; break;}
        wEff=wD;
        for(i=0;i<ng;i++)
        {
          int j;
          for(j=0,I0=0;j<ng;j++) 
          { double Tj=1/sqrt(zc+P[j]*(z[i]-zc)), lTj=log(Tj);
            double wBath,Mm;
            if(lTj<=ltGrid[0])
            { wBath=mediatorArr[n].wBath[0];
              Mm=mediatorArr[n].mass[0];
              eta=mediatorArr[n].eta[0];
            }  
            else
            {  
               wBath=polint3(lTj,tGridDim,ltGrid,mediatorArr[n].wBath);
               Mm=polint3(lTj,tGridDim,ltGrid,mediatorArr[n].mass);
               eta=polint3(lTj,tGridDim,ltGrid,mediatorArr[n].eta);
            }  
//printf("i=%d j=%d zc=%E z[i]=%E P[j]=%E Tij=%E  eta=%E\n",i,j,zc,z[i], P[j], Tj,eta);            
            I0+=0.5*cG[j]*wBath*Tj*Tj/Hubble(Tj)/gammaLor(Mm/Tj,eta);            
          }
//printf("I0=%e\n", I0);          
          z[i]=zc+L[i]/I0; 
        }
     }
     mediatorArr[n].wEff[nT]=wEff;
#else      
     mediatorArr[n].wEff[nT]=mediatorArr[n].wBath[nT];
#endif
     mediatorArr[n].wEff[nT]*=mediatorArr[n].fKK[nT]; 
//printf(">T=%E eta=%E  wBat=%E  \n",T,mediatorArr[n].eta[nT],mediatorArr[n].wBath[nT]);        
     
   }    
 }
}

//======================= Intervals ===================

static double * intervals=NULL;
static int nIntervals=0; 


static void printIntervals(void) 
{ printf("Intervals: ");
  for(int i=0;i<nIntervals;i++) printf("[%e,%e] ",intervals[2*i],intervals[2*i+1]);
  printf("\n");
} 




// ========   SECTION  COLLISIONS  ======== 

static double sqrtSminT(double m1,double m2, double m3, double m4, double T)
{

  T=0;
  double smin=m1+m2;
  if(smin< m3+m4) smin=m3+m4;
  if(T<=0) return  smin;
//return T*sqrt(2);
  
  double T2=T*T;
  double s1=smin,s2=smin+T,sx,pp;
  for(;;)
  {  
    pp=2*decayPcm(s2,m1,m2)*decayPcm(s2,m3,m4);
    if(pp>T2) break;
    s1=s2; s2+=T/2;
  }
  sx=(s1+s2)/2;
  for(;;)
  {
    pp=2*decayPcm(sx,m1,m2)*decayPcm(sx,m3,m4);
    if(pp<T2) s1=sx;  else  s2=sx;  
    if(fabs(s2-s1)<1E-5*s2) { return s2; }
    sx=0.5*(s1+s2);
  }
}

//========= Multiple integrtion

static double eps=1E-2;


static double sIntegrand_y(double y, void*arg_)
{
  if(y==0) return 0;
  
  par_22* arg=arg_;
//  if(arg->err) return 0;
  
  double sqrtS,J;

  long double z=y*(2-y);
  
  sqrtS=-3*T*logl(z);
  
//  y=  1-sqrt(1-exp(-sqrtS/(3*T)));
  
  J=6*T*(1-y)/z;
//printf("1-y=%E Tx=%E sqrtS=%E\n", 1-y,Tx,sqrtS);  
  double sqme_Int;
  int err=0;


  if(Qaddress)
  {  *Qaddress =sqrtS;
     calcMainFunc();  
     passParameters(arg->cc);
     for(int n=0;n<nMed;n++) if(mediatorArr[n].my) arg->cc->interface->va[mediatorArr[n].widthPos]=mediatorArr[n].Twidth;
     mass22_par(arg,T);
  }


  if(kin22_par(arg,sqrtS, sqrt(4*M_PI*parton_alpha(sqrtS)))){  printf("kin22 problem sqrtS=%e \n",sqrtS); return 0;}
  err=0;
  arg->T=T;
  sqme_Int=sqmeIntDel(arg,eps/3);
  
  double res=  sqme_Int/(32*M_PI*sqrtS)*T/(8*M_PI*M_PI*M_PI*M_PI)*arg->PcmIn*arg->PcmOut
  *K1to2(sqrtS/T,arg->pmass[0]/T,arg->pmass[1]/T,0,arg->eta[0],arg->eta[1])
//*bessK1(sqrtS/T)
*2*sqrtS*J;

  

  if((arg->eta[2] || arg->eta[3]) && sqrtS/T<25  ) 
      res*=K1to2(sqrtS/T,arg->pmass[2]/T,arg->pmass[3]/T,0,arg->eta[2],arg->eta[3])/
           K1to2(sqrtS/T,arg->pmass[2]/T,arg->pmass[3]/T,0,0,0);

  return res;
 
}

static double lnTIntegrand(double lnT,void*arg_)
{ 
  T=exp(lnT);
  par_22 *arg=arg_;
  double E=2*M_PI*M_PI/45*hEff(T)*T*T*T;                      // entropy       
  double H=sqrt(8*M_PI/3.*M_PI*M_PI/30.*gEff(T))*T*T/MPlanck;  //Huble 
  double J=(1+hEffLnDiff(T)/3)/H;                             // dt => dlnT
 

  if(Taddress) *Taddress=T;
  if(Qaddress) *Qaddress=1;
  if(Taddress || Qaddress)
  {  calcMainFunc(); 
     cleanDecayTable();
     passParameters(arg->cc);
  }   
  mass22_parDel(arg,T);
  double sqrtSmin;
/*  
  sqrtSmin=arg->pmass[0]+arg->pmass[1];  
//printf("arg->pmass[0]+arg->pmass[1] = %E %E\n", (double)(arg->pmass[0]), (double)(arg->pmass[1]));  
  double dM=0;
  for(int i=0;i<2;i++) if(arg->pmass[i]==0)
  {  if(arg->pdg[i]==21) dM+= T*sqrt(4*M_PI*alphaQCD(T));
     if(arg->pdg[i]==22) dM+= T*0.31/sqrt(6);
  }   
  sqrtSmin+=dM;   
  
  double mout=arg->pmass[2]+arg->pmass[3];
  if(mout>sqrtSmin) sqrtSmin=mout;
*/  
  sqrtSmin=sqrtSminT( arg->pmass[0], arg->pmass[1],arg->pmass[2],arg->pmass[3],T*tTcut);

//printf("Masses=%e %e %e %e\n", (double)arg->pmass[0], (double)arg->pmass[1],(double)arg->pmass[2],(double)arg->pmass[3]);  
  
  sqrtSmin*=1.000001;
  nIntervals=1;
  intervals=realloc(intervals,2*sizeof(double)); 
  intervals[0]=0;  intervals[1]=1-sqrt(1-exp(-sqrtSmin/(3*T)));  intervals[1]=exp(-sqrtSmin/(3*T))/(1+sqrt(1-exp(-sqrtSmin/(3*T))));
//printf("sqrtSmin=%E  T=%E  intervals[0]=%E  intervals[1]=%E\n", sqrtSmin,T,intervals[0], intervals[1]);
  
  double sum=0;     
  for(int n=0;n<nMed;n++) if(mediatorArr[n].my)
  { double Mm=varValues[mediatorArr[n].massPos];
    
//    if(sqrtSmin<=Mm)
    { 
      double mu1=arg->pmass[0]/T, mu2=arg->pmass[1]/T;
      double wEff,wIwD, eta;
      if(lnT<ltGrid[0]) 
      {  wEff=mediatorArr[n].wEff[0];
         wIwD=mediatorArr[n].wIwD[0];
         eta=mediatorArr[n].eta[0];
      } else   
      {
         wEff=polint3(lnT,tGridDim,ltGrid,mediatorArr[n].wEff);
         wIwD=polint3(lnT,tGridDim,ltGrid,mediatorArr[n].wIwD); 
         eta=polint3(lnT,tGridDim,ltGrid,mediatorArr[n].eta);
      }
//if(Mm/T <25)   printf("eta=%E K=%e\n",eta,K1to2(Mm/T,mu1,mu2,0,arg->eta[0],arg->eta[1])/K1to2(Mm/T,mu1,mu2,eta,arg->eta[0],arg->eta[1]));      
//      if(Mm/T <25)  wEff*=K1to2(Mm/T,mu1,mu2,0,arg->eta[0],arg->eta[1])/K1to2(Mm/T,mu1,mu2,eta,arg->eta[0],arg->eta[1]);   
       arg->cc->interface->va[mediatorArr[n].widthPos]=wEff;
//       printf("T=%E  %s=%E\n", T,   arg->cc->interface->varName[mediatorArr[n].widthPos],wEff);
       mediatorArr[n].Twidth=wEff;
//printf("T+%E M=%E width=%E eta=%E\n",T,Mm,wEff,eta);
       
/*   
      {      
         double d0 = e0>100*wEff/Tx? e0: 100*wEff/Tx; 
//printf("T=%E d0*Tx=%e wEff=%E\n", T,d0*Tx,wEff);           
         delInterval(exp(-Mm/Tx -d0)  ,exp(-Mm/Tx+d0),&intervals,&nIntervals);
         double c=1;   
         if((arg->eta[2] || arg->eta[3]) && Mm/T<25  ) 
            c=K1to2(Mm/T,arg->pmass[2]/T,arg->pmass[3]/T,0,arg->eta[2],arg->eta[3])/
              K1to2(Mm/T,arg->pmass[2]/T,arg->pmass[3]/T,0,0,0);
         sum+= c*(mediatorArr[n].spin2+1)*mediatorArr[n].cdim*T*Mm*Mm*wIwD*K1to2(Mm/T,mu1,mu2,0,arg->eta[0],arg->eta[1])/( (2*M_PI*M_PI)*wEff) ;
      }
*/      
    
    } 
//    else 
//    { mediatorArr[n].Twidth=polint3(lnT,tGridDim,ltGrid,mediatorArr[n].wEff);     
//    }
  }
  double NdfIn=arg->ndf[0]*arg->ndf[1];
  if(arg->pdg[0]==arg->pdg[1]) NdfIn/=2;
  sum/=NdfIn; 
//printf("nIntervals=%d\n", nIntervals);  
#ifndef ONSHELL 
  for(int i=0;i<nIntervals;i++)
  { int err;
    sum+= simpson_arg(sIntegrand_y,arg,intervals[2*i], intervals[2*i+1],1E-3,&err);
    if(err)
    {  arg->err=arg->err|(2*8*err);
#ifdef ERROR_PLOT
      { printf("Warning from Simpson. Precision is not reached (freezein.c line 1222 code=%d  arg->err=%d ErrIS=%d )\n",err,arg->err,ErrIS);
       printf("i=%d interval: %e %e  \n",i, intervals[2*i],intervals[2*i+1]);
       displayPlot("s_Integrand_y", "y",intervals[2*i], intervals[2*i+1], 0,1, "dI/dy", 0, sIntegrand_y,arg);
       exit(0);
      } 
#endif
    } 
  }
#endif

  return sum*J/E; 
}  



static int printFi22Error(FILE*f, int err)
{
 switch(err)
     { 
       case  1: fprintf(f,"process  is absent\n"); break;
       case  2: fprintf(f,"2->2 type process is expected\n");break;
       case  3: fprintf(f,"can not calculate local parameters\n");break;
       case  4: fprintf(f,"reheating temperature is smaler than Tend  too or Tend=0\n");break;       
       case  5: fprintf(f,"there are incoming feeble particles\n"); break;
       case  6: fprintf(f,"there is no odd  feeble particles among outgoing ones\n");break;

       case   7: fprintf(f,"Lost of precision in  temperature integrand\n");break;
       case   8: fprintf(f,"Pole in temperature integrand\n");break;
       case   9: fprintf(f,"NaN in temperature intergrand\n");break;
       case  10: fprintf(f,"Loss of precision in  sqrt(s) integrand\n");break;
       case  11: fprintf(f,"Pole in sqrt(s) integrand\n");break;
       case  12: fprintf(f,"NaN  in sqrt(s) intergrand\n");break;
       case  13: fprintf(f,"Lost of precision in angle integration\n");break;
       case  14: fprintf(f,"Pole in angle  integrand\n");break;
       case  15: fprintf(f,"NaN  in angle  intergrand\n");break;
       case  16: fprintf(f,"loss of precision caused by diagramm cancelation\n");
     }
}

static double totOmegaCoeff;

double  darkOmegaFi22(double Tr, char *Proc, char * fiPrtcl,  int *err_)
{  
   int plot=0;
   int err=0, intErr=0;
   double omega=0; 
   int  NDM=0; // number of "fiPrtcl"-type particles anong outgoing ones
   int NConj=1;
   double Cin;
   par_22 arg;
   double Tmin=Tend;

   if(err_) *err_=0;   
   numout*cc=newProcess(Proc);
   if(cc->interface->nvar==0) return 0;
// printf("Proc=%s\n",Proc);      
   int i;
   for(Qaddress=NULL,i=0;i<nModelVars;i++) if(strcmp(varNames[i],"Q")==0) { Qaddress=varValues+i; break;}
   for(Taddress=NULL,i=0;i<nModelVars;i++) if(strcmp(varNames[i],"T")==0) { Taddress=varValues+i; break;}

   err=init22_par(&arg,cc,1);  //process   is absent || 2->2 type process is expected 
   if(!err)
   { 
      if(Taddress)     
      { *Taddress=0;
         calcMainFunc();
         cleanDecayTable();
      }
      if(passParameters(cc)) err=3; // Can not calculate local parameters for process 
   }
     
   if(!err)   
   { 
      for(int i=0;i<4;i++) if(isFeeble(cc->interface->pinf(1,1+i,arg.pmass+i,NULL))) arg.eta[i]=0; else arg.eta[i]=1-2*(arg.spin2[i]&1);
      if(arg.eta[0]==0 || arg.eta[1]==0)  err=5; //  There are incoming feeble particles
      for(int i=2;i<4;i++) 
      { char*p=cc->interface->pinf(1,1+i,NULL,NULL);
        if(strcmp(p,fiPrtcl)==0) NDM++;
        else 
        {
           char*ap=antiParticle(p);
           if(strcmp(ap,fiPrtcl)==0) NDM++;
        }        
      } 
      if(NDM==0) err=6;  // There is no odd  feeble particles among outgoing ones   
   }

   if(!err)  if(!saveMediators && makeLnTgrid(Tr)) err=4;

   if(!err)
   { *(cc->interface->BWrange)=1000;  //!!!         
     fillMediatorArr(Tr,&arg);
     int neutral[4];
     for(i=0;i<4;i++) cc->interface->pinfAux(1,1+i,NULL,NULL,neutral+i,NULL);
     if( arg.pdg[0]+arg.pdg[1]!=0 && (!neutral[0] || !neutral[1]))NConj=2;
     if( arg.pdg[2]+arg.pdg[3]!=0 && (!neutral[2] || !neutral[3]))NConj=2;  //??
     
     Cin=1;
     if(arg.pdg[0]==arg.pdg[1]) Cin=0.5;
     totOmegaCoeff=Cin*arg.ndf[0]*arg.ndf[1]*NConj*NDM;
   }
   if(!err)
   {   
//displayPlot("lnTIntegrand", "log(T)",1.,7.,0,1,"integrand", 0, lnTIntegrand,&arg);        
       double Y=0;
       for(int nT=tGridDim-2;nT>=0;nT--)  
       {  double lnT1= ltGrid[nT+1], lnT2=ltGrid[nT];
          if(nT==tGridDim-2) Y+=simpson_arg(lnTIntegrand,&arg,lnT2,lnT1,1E-3,NULL);  
          else              Y+=gauss_arg(lnTIntegrand,&arg,lnT2,lnT1,4);
          yfi[nT]+=Y*totOmegaCoeff;                                                                                                                                                                            
       }      
       arg.err=0; 
       omega= pMass(fiPrtcl)*EntropyNow/RhoCrit100*totOmegaCoeff*Y;
       if(err) arg.err=arg.err|(2*8*8*err);
             if(arg.err & 2*8*8*4)  err=7; 
        else if(arg.err & 2*8*8*2)  err=8;
        else if(arg.err & 2*8*8*1)  err=9;
        else if(arg.err & 2*8*  4)  err=10;
        else if(arg.err & 2*8*  2)  err=11;
        else if(arg.err & 2*8*  1)  err=12;
        else if(arg.err & 2*    4)  err=13;
        else if(arg.err & 2*    2)  err=14;
        else if(arg.err & 2*    1)  err=15;
        else if(arg.err &       1)  err=16;                                                                                                                                                                                                                                                                                    
   }
   if(!saveMediators) cleanMediatorArr();
   if(err_)  *err_=err; else  printFi22Error(stdout,err);
   Tstart=Tr;
   return omega;
}
//===================== Summation =====================

static void addPrtcl(char **all, int n)
{ 
      *all=realloc(*all,strlen(*all)+2+strlen(ModelPrtcls[n].name));
      sprintf(*all+strlen(*all),",%s",ModelPrtcls[n].name);
      if(strcmp(ModelPrtcls[n].name,ModelPrtcls[n].aname))
      { *all=realloc(*all,strlen(*all)+2+strlen(ModelPrtcls[n].aname));
        sprintf(*all+strlen(*all),",%s",ModelPrtcls[n].aname);
      }    
}       



double darkOmegaFi(double TR,char* fiPrtcl,   int *err)
{  

  int p=pTabPos(fiPrtcl);
  if(p==0) { printf("darkOmegaFi: unknown particle %s\n", fiPrtcl); if(err) *err=20; return 0; }
  if(! isFeeble(fiPrtcl)) { printf("darkOmegaFi: particle %s is not feeble\n", fiPrtcl); if(err) *err=21; return 0; } 
  p=abs(p)-1;
  NAME=ModelPrtcls[p].name; aNAME=ModelPrtcls[p].aname;

  char * feebleDm=malloc(1);
  feebleDm[0]=0;
  addPrtcl(&feebleDm,p);
  double minFeebleMass=pMass(fiPrtcl);
  
  char * bath=malloc(1);
  bath[0]=0;
  double omega=0,dOmega;
  
  int nFiCh=0;
  omegaFiCh=realloc(omegaFiCh,sizeof(aChannel));
  
  int VZdecayMem=VZdecay, VWdecayMem=VWdecay;  VZdecay=0;VWdecay=0;
  if(VZdecayMem|| VWdecayMem) cleanDecayTable();

  
 
  if(minFeebleMass<0) 
  { printf("Feeble particles are not detected\n"); 
    free(feebleDm);
    free(bath);
    if(err)*err=1;
    return 0;
  }
  *err=0;
  
  makeLnTgrid(TR); 
  
  for(int n=0;n<nModelParticles;n++) if(!isFeeble(ModelPrtcls[n].name))
  { double m=pMass(ModelPrtcls[n].name);
    if(m>minFeebleMass) 
    {  double  brTot, brBath, brSig1;
       txtList L;
       pWidth(ModelPrtcls[n].name,&L);
       int nFiCh_=nFiCh;
       getDecBrKE(m/3, L, 1, &brTot, &brBath,&brSig1,&nFiCh);
       if(brBath==0 && brSig1>0)
       {
         dOmega=darkOmegaFiDecay(TR,ModelPrtcls[n].name,NAME);
         omega+=dOmega;         
         if(dOmega) for(int i=nFiCh_;i<nFiCh;i++) omegaFiCh[i].weight*=dOmega; else nFiCh=nFiCh_;
       }  else { addPrtcl(&bath,n); nFiCh=nFiCh_;}
    } else  addPrtcl(&bath,n);
  }
//printf("FiCh=%d\n", nFiCh);
//printf("feebleDm=%s\n", feebleDm);
//printf("bath=%s\n", bath);

  if(!strlen(feebleDm) || !strlen(bath)) return 0;
  
  char * command=malloc(strlen(compDir) + strlen(calchepDir) + strlen(libDir) 
                       +strlen(feebleDm)+  strlen(bath)+ 200);
  int delWorkDir=prepareWorkPlace();
      
  sprintf(command,"cd %s;"
           " %s/bin/s_calchep -blind \"{{allBath,allBath->feeblDm,1*x{%s{%s{{{[[{0\" >/dev/null;"
           " if(test $? -eq 0) then mv results/list_prc.txt %s; fi",
            compDir,calchepDir,bath+1,feebleDm+1,libDir);
  system(command); 
  free(command);
//printf("allBath=%s\n", bath+1);  
  if(delWorkDir) cleanWorkPlace(); 

  char*fname0=malloc(strlen(libDir)+50);   
  sprintf(fname0,"%s/list_prc.txt",libDir);
  FILE*f0=fopen(fname0,"r");
  free(fname0);
  saveMediators=1;
  
  if(f0)
  { 
    char p[4][P_NAME_SIZE+1];
    while(4==fscanf(f0," %[^ ,] , %s -> %[^ ,] , %s",p[0],p[1],p[2],p[3]))
    {  int i,err; 
//printf(" %s %s -> %s %s\n", p[0],p[1],p[2],p[3]);
       if(strcmp(p[0],p[1])>0) continue;   
       int n[4],n_[4];
       for(i=0;i<4;i++)
       {  n[i]=pTabPos(p[i]);
          if(strcmp(ModelPrtcls[abs(n[i])-1].name, ModelPrtcls[abs(n[i])-1].aname)==0) n_[i]=n[i];else n_[i]=-n[i];
       }
       if(n[0] >n[1] ) { int nn=n[0];  n[0]=n[1];   n[1]=nn;} 
       if(n[2] >n[3] ) { int nn=n[2];  n[2]=n[3];   n[2]=nn;}
       if(n_[0]>n_[1]) { int nn=n_[0]; n_[0]=n_[1]; n_[1]=nn;} 
       if(n_[2]>n_[3]) { int nn=n_[2]; n_[2]=n_[3]; n_[2]=nn;}
       
       for(i=0;i<4;i++) if(n[i]!=n_[i]) break;
       if(i<4 && n[i]<n_[i]) continue;
       char proc[100];
       sprintf(proc,"%s,%s->%s,%s", p[0],p[1],p[2],p[3]);
       
       dOmega=darkOmegaFi22(TR,proc,NAME,&err);
       
       if(dOmega||err)
       {         
         for(int i=0;i<4;i++) if(n[i]>0) omegaFiCh[nFiCh].prtcl[i]=ModelPrtcls[n[i]-1].name;
         else  omegaFiCh[nFiCh].prtcl[i]=ModelPrtcls[-n[i]-1].aname;
         omegaFiCh[nFiCh].prtcl[4]=NULL;
         omegaFiCh[nFiCh].weight=dOmega;
         omegaFiCh[nFiCh].err=err;
         omegaFiCh[nFiCh].nin=2;

         nFiCh++;
         omegaFiCh=realloc(omegaFiCh,(nFiCh+1)*sizeof(aChannel));
                  omega+=dOmega;

       } 
    }
    fclose(f0);  
  }
  
//  if(omega) for(int i=0;i<nFiCh;i++) omegaFiCh[i].weight/=omega;
  omegaFiCh[nFiCh].weight=0; omegaFiCh[nFiCh].nin=0;
  
  
  for(int i=0;i<5;i++) omegaFiCh[nFiCh].prtcl[i]=NULL;
  free(feebleDm);
  free(bath);
  cleanMediatorArr();

   if(VZdecayMem|| VWdecayMem)
   {  cleanDecayTable();
      VZdecay=VZdecayMem, VWdecay=VWdecayMem;
   }
     
  if(err)*err=0;      
  saveMediators=0;
  Tstart=TR;
  return omega;
}

void printChannelsFi(double cut, int prcn, FILE * f)
{  
  if(omegaFiCh && omegaFiCh[0].prtcl[0])
  {  
    fprintf(f,"# Channels which contribute to omega h^2 via freeze-in\n");

// calc omega 
    double omega=0; 
    for(int i=0; omegaFiCh[i].prtcl[0];i++) omega+=omegaFiCh[i].weight;

// sorting 

    for(int i=0; omegaFiCh[i+1].prtcl[0];)
    { 
       if(omegaFiCh[i].weight < omegaFiCh[i+1].weight)
       { aChannel tmp=omegaFiCh[i];
         omegaFiCh[i]=omegaFiCh[i+1];
         omegaFiCh[i+1]=tmp;
         if(i) i--; else i++;
       } else i++;
    } 

    for(int i=0; omegaFiCh[i].prtcl[0];i++) if(omegaFiCh[i].weight/omega>cut || omegaFiCh[i].err )
    {  if(prcn) fprintf(f," %6.2f%%  ", 100*omegaFiCh[i].weight/omega);
        else    fprintf(f," %.3E  ", omegaFiCh[i].weight/omega);
        char proc[40];
        sprintf(proc,"%s",omegaFiCh[i].prtcl[0]);
        if(omegaFiCh[i].nin==1)
        {  sprintf(proc+strlen(proc)," -> %s, %s", omegaFiCh[i].prtcl[1],omegaFiCh[i].prtcl[2]);
           for(int k=3; omegaFiCh[i].prtcl[k];k++) sprintf(proc+strlen(proc),", %s", omegaFiCh[i].prtcl[k]);
        }
        else  sprintf(proc+strlen(proc),", %s ->  %s,%s ", omegaFiCh[i].prtcl[1],omegaFiCh[i].prtcl[2],omegaFiCh[i].prtcl[3]);
        fprintf(f,"%-30.30s ",proc); 
        if(omegaFiCh[i].err) { fprintf(f,"Warning %d:",omegaFiCh[i].err);  printFi22Error(f,omegaFiCh[i].err);}
        else fprintf(f,"\n"); 
    }
  }  
}


//================================= testing 22 ==================

typedef struct{ numout*cc; REAL m[4]; double eta[4], sqrtSmin, C, nWidth, T;}  frin22Par;


static double NabDmDmInt(double x, frin22Par* arg)
{  int err;

   if(x==0 || x==1) return 0;

   double T=arg->T;
   
   double z=x*(2-x);
   double sqrtS=arg->sqrtSmin-3*T*log(z); 
   
   double J=6*T*(1-x)/z;
   double pcm=decayPcm(sqrtS,arg->m[0],arg->m[1]);
   double M[4];
   for(int i=0;i<4;i++) M[i]=arg->m[i];
   double Kout=1;
   if(mTcut>0)
   { double dm2=pow(mTcut*T,2);
     for(int i=0;i<4;i++) if(arg->eta[i]) M[i]=sqrt(M[i]*M[i]+dm2);
     Kout=decayPcm(sqrtS,M[2],M[3])/decayPcm(sqrtS,arg->m[2],arg->m[3]);
   }  
   double Pcm=decayPcm(sqrtS,M[0],M[1]); 
   double cos1=-1,cos2=1;
   double pp;
   if(tTcut>0)
   { pp=2*pcm*decayPcm(sqrtS,arg->m[2],arg->m[3]);
     
     for(int n=1;;n++)
     {  
        int m,w,pnum;
        char*s=arg->cc->interface->den_info(1,n,&m,&w,&pnum);
        if(!s) break;
//        if(!m) continue;
        double kappa=Tkappa[pnum]*T;
        if(kappa==0) continue;
         double cos0=1-kappa*kappa/pp; 
        if(s[0]==1 && s[1]==3)
        {
          if(cos2>cos0) cos2=cos0;       
        } else  if(s[0]==1 && s[1]==4)    
        {
          if(-cos1>cos0) cos1=-cos0;
        }   
     }

//   cos1=-cos0;
//   cos2=cos0;   
   } 
//printf("T=%E pcm/T=%E cos= %E %E\n", T,pp/T,cos1,cos2);   
   if(cos1>=cos2) return 0;
//   if(cos1>0.999 || cos2<-0.999) return 0;  
   GGscale=sqrtS;
   double K1=K1to2(sqrtS/T,M[0]/T,M[1]/T,0,arg->eta[0],arg->eta[1]);
//double K1=bessK1(sqrtS/T);
   double  res= J*pcm*Pcm*sqrtS*sqrtS*K1*cs22(arg->cc,1,pcm, cos1, cos2, NULL)*Kout/3.8937966E8/(cos2-cos1)*2;

#ifndef NOSTISTICS
   if((arg->eta[2] || arg->eta[3]) && sqrtS/T<5  ) 
      res*=K1to2(sqrtS/T,M[2]/T,M[3]/T,0,arg->eta[2],arg->eta[3])/
           K1to2(sqrtS/T,M[2]/T,M[3]/T,0,0,0);
#endif 

   return res;
}


static double sqrtS2x(double sqrtS, frin22Par*arg)
{ double z=exp( (arg->sqrtSmin-sqrtS)/(3*arg->T));
//printf("z=%E\n",z);
  return 1-sqrt(1-z);  
}


static double TdYdTfreezeIn(double T, frin22Par*arg)      // integration over s 
{ 
   double s=2*M_PI*M_PI*T*T*T*hEff(T)/45;
   double H=Hubble(T);
   int err;

   if(Taddress) 
   {  *Taddress=T; 
      calcMainFunc(); cleanDecayTable();
      passParameters(arg->cc); 
      for(int i=0;i<4;i++) { arg->cc->interface->pinf(1,i+1,arg->m+i,NULL); arg->m[i]=Fabs(arg->m[i]);}
   }  
   arg->sqrtSmin=sqrtSminT( arg->m[0], arg->m[1], arg->m[2], arg->m[3],T*tTcut);  
   if(mTcut>0)
   {  double dm2=pow(mTcut*T,2);
      double sqrtS0=sqrt( pow(arg->m[0],2)+dm2) + sqrt(pow(arg->m[1],2)+dm2);
      if(sqrtS0>arg->sqrtSmin) arg->sqrtSmin=sqrtS0;
   } 
   
//printf("T=%E sqrtSmin/T=%E\n",T, arg-> sqrtSmin/T);   
   arg->T=T;
   
   int nIntervals=1;
   double *intervals=malloc(2*sizeof(double));
   intervals[0]=0;  intervals[1]=1;

   if(arg->nWidth>0)for(int n=1;;n++)
   {  
      int m,w,pnum;
      char*s=arg->cc->interface->den_info(1,n,&m,&w,&pnum);
      if(!s) break;
      if(!m) continue;
      if(s[0]!=1 || s[1]!=2) continue;
           
      double mm= fabs(arg->cc->interface->va[m]);
      double ww=0; 
      if(w) ww=arg->cc->interface->va[w];  
      if(ww==0) 
      { if( mm>arg->sqrtSmin ) { printf(" mm=%E sqrtSmin=%E     zero width in s-channel propagator \n",mm, (double)arg->sqrtSmin, ModelPrtcls[abs(pnum)].name ) ; exit(0); }
         else continue;
      }  
      
      if(sqrt(mm*mm+ arg->nWidth*ww*mm)<=arg->sqrtSmin) continue;
      
      double x1=sqrtS2x( sqrt(mm*mm+ arg->nWidth*ww*mm)  , arg);
      double x2;
      if( mm*mm- arg->nWidth*ww*mm < arg->sqrtSmin) x2=1; else x2=sqrtS2x( sqrt(mm*mm- arg->nWidth*ww*mm), arg);
      delInterval(x1,x2, &intervals, &nIntervals);      
   }

   double res=0;
   for(int k=0;k<nIntervals;k++) res+=T*arg->C/(4*pow(M_PI,4))*simpson_arg( NabDmDmInt,arg,intervals[2*k],intervals[2*k+1],1E-3,&err)*(1+hEffLnDiff(T)/3)/(s*H);
//   printf("nIntervals=%d intervals[0]=%E, intervals[1]=%E\n",nIntervals, intervals[0], intervals[1]);  
//   res=T*arg->C/(4*pow(M_PI,4))*simpson_arg( NabDmDmInt,arg,0,1,1E-3,&err)*(1+hEffLnDiff(T)/3)/(s*H);
   free(intervals);
   return res;
}

static double TdYdTfreezeInLn(double lnT, frin22Par*arg){   double T=exp(lnT);   return  TdYdTfreezeIn(T, arg);}
static double dYfreezeIn10(double T, frin22Par*arg) { return  log(10)*TdYdTfreezeIn(T, arg);}


double  YfreezeIn22(char*proc, double T0, double TR, int plot, int * err)
{
   numout*cc=newProcess(proc); 
   if(cc==NULL)                                          { printf("YfreezeIn22: no such process '%s'\n",proc);          if(err)*err=1;  return 0;}
   if(cc->interface->nprc!=1)                            { printf("YfreezeIn22 does not work with multiple process\n"); if(err)*err=2;  return 0;} 
   if(cc->interface->nin!=2 || cc->interface->nout!=2)   { printf("YfreezeIn22  works only with 2->2 processes\n");     if(err)*err=3;  return 0;}  
   if(isFeeble(cc->interface->pinf(1,1,NULL,NULL)) 
    ||isFeeble(cc->interface->pinf(1,2,NULL,NULL)))      { printf("freezeIn22: incoming particles can not be feeble\n");if(err)*err=4;  return 0;}
      
   frin22Par par;
   par.cc=cc;
   double Tmem;
   Taddress=varAddress("T");
   if(!Taddress)  printf(" 'T' variable is absent\n");  
   else Tmem=*Taddress; 
   passParameters(par.cc);
   par.nWidth=0;       
   int pdg[4];
   char*name[4];
   for(int i=0;i<4;i++) 
   {  name[i]=par.cc->interface->pinf(1,i+1,par.m+i,pdg+i); 
      par.m[i]=Fabs(par.m[i]);
      int spin2;
      qNumbers(name[i], &spin2,NULL, NULL);
      if(isFeeble(name[i])) par.eta[i]=0; else par.eta[i]=1-2*(spin2&1);
   }
   
   int ndf[2];
   for(int i=0;i<2;i++)  par.cc->interface->pinfAux(1, i+1,NULL,NULL,NULL,ndf+i);
   par.C=ndf[0]*ndf[1];
   if(pdg[0]==pdg[1]) par.C/=2;
   par.C*=2; // because we have assume 2 DM particles
   double res;
   res=simpson_arg(TdYdTfreezeInLn, &par, log(T0),log(TR),1E-3,NULL);    // temperature integration

//   if(plot) displayPlot("dY/dT","T",T0,TR,1,1,"dY/dT",0,dYfreezeIn,&par); 
   if(plot) displayPlot("dY/dlog(T)","T",T0,TR,1,1,"dY/dlog10(T)",0,dYfreezeIn10,&par); 
   if(Taddress) { *Taddress=Tmem;  calcMainFunc(); cleanDecayTable(); } 
   return res;
}
 