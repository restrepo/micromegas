#include "micromegas.h"
#include "micromegas_aux.h"

//====================   g1eff  & h1eff ================

typedef   struct{ double mu,eta;}    geff_struct;


static double g1eff_int(double x,void *par_)
{  if(x==0 || x==1) return 0;
   geff_struct*par=par_;
   double p=2*log(x);
   double e=sqrt(p*p+par->mu*par->mu);
   double den;

   if(par->eta==-1 && e<1E-8) den=e*(1+e/2); else den=exp(e)+par->eta;

   return   e*p*p/den/x/M_PI/M_PI;
}


static double h1eff_int(double x,void *par_)
{  if(x==0 || x==1) return 0;
   geff_struct*par=par_;
   double p=2*log(x);
   double e=sqrt(p*p+par->mu*par->mu);
   double den;

   if(par->eta==-1 && e<1E-8) den=e*(1+e/2); else den=exp(e)+par->eta;
   return  (par->mu*par->mu+4./3.*p*p)*p*p/den/x/M_PI/M_PI/e;
}


double g1eff(double mu,int eta)
{ if(!isfinite(mu)) return 0;
  if(mu>80) return 0.2*pow(mu,2.5)*exp(-mu);
  geff_struct par;
  par.eta=-eta;  //!!!
  par.mu=mu;
  return simpson_arg(g1eff_int,&par, 0, 1,1E-4,NULL)/(M_PI*M_PI/30);
}


double h1eff(double mu,int eta)
{ if(!isfinite(mu)) return 0;
  if(mu>80) return 0.15*pow(mu,2.5)*exp(-mu);
  geff_struct par;
  par.eta=-eta;  //!!!
  par.mu=mu;
  return simpson_arg(h1eff_int,&par, 0, 1,1E-4,NULL)/(2*M_PI*M_PI/45);
}


double Hubble(double T)
{
   double muDE=2.24E-12, muDM=5.19E-10;    //GeV
//muDM=5.18E-10;
   return sqrt(8*M_PI/3.*(
      M_PI*M_PI/30.*gEff(T)*pow(T,4)         // radiation
     + pow(muDE,4)                           // darkEnergy
     + muDM* 2*M_PI*M_PI/45*hEff(T)*pow(T,3) // darkMatter
    ))/MPlanck;
}

static double HubbleIntegrand( double lnT) {  double T=exp(lnT); return (1+ hEffLnDiff(T)/3)/Hubble(T); }
double HubbleTime(double T1, double T2) { return   simpson(HubbleIntegrand, log(T2), log(T1), 1E-4,NULL)*6.582E-25; }

static double freeSteamInt( double lnT, void*par_)
{ double *par=par_;
  double  pi0=par[0], T1=par[1], T2=par[2];
  double T=exp(lnT);

  double a=T2/T*pow(hEff(T2)/hEff(T), 1./3.);
  double pi=pi0*T/T1*pow(hEff(T)/hEff(T1), 1./3.);
  double v=pi/sqrt(1+pi*pi);

  return v/a/Hubble(T)*(1+hEffLnDiff(T)/3);
}

double freeStreaming(double pi, double T1, double T2)
{
   double par[3];
   par[0]=pi, par[1]=T1, par[2]=T2;
   double res= simpson_arg(freeSteamInt, par, log(T2), log(T1), 1E-3,NULL);
//1/GeV = 6.395 E-39 Mpc
   return res*T2/T2_73K* pow(hEff(T2)/hEff(T2_73K),1./3.)*6.395E-39;
}

static int Tdim=0;
static double *t_=NULL, *heff_=NULL, *geff_=NULL, *geff2_=NULL, *s3_=NULL,
 *heff_ln_diff_=NULL;

int loadHeffGeff(char*fname)
{  double *tabs[5];
   int nRec,nCol,i;

   if(fname==NULL)
   {
     if(heff_==NULL )
     {
        char*fName=malloc(strlen(micrO)+50);
        sprintf(fName,"%s/Data/hgEff/DHS.thg",micrO); // DHS
//        sprintf(fName,"%s/Data/hgEff/GG.thg",micrO); 
        nRec=readTable(fName,&nCol,tabs);
        free(fName);
        if(nRec<=0) return nRec;
     }
     else return Tdim;
   }
   else nRec=readTable(fname,&nCol,tabs);
   if(nRec==0)
   {
      char*fName=malloc(strlen(micrO)+strlen(fname)+20);
      sprintf(fName,"%s/Data/hgEff/%s",micrO,fname);
      nRec=readTable(fName,&nCol,tabs);
      free(fName);
   }

   if(nRec<=0) return nRec;
   if(nCol!=3)
   { for(i=0;i<nCol;i++) free(tabs[i]);
     return 0;
   }
   if(t_)free(t_);  if(heff_)free(heff_);   if(geff_)free(geff_);  if(s3_)free(s3_);
   if(heff_ln_diff_) free(heff_ln_diff_);

   if(tabs[0][0]> tabs[0][nRec-1])
   for(int i=0;i<nRec/2;i++)
   { int im=nRec-1-i;
     double b[3];
     for(int j=0;j<3;j++) b[j]=tabs[j][i];
     for(int j=0;j<3;j++) tabs[j][i]=tabs[j][im];
     for(int j=0;j<3;j++) tabs[j][im]=b[j];
   }

   t_=tabs[0];
   heff_=tabs[1];
   geff_=tabs[2];

   s3_=malloc(nRec*sizeof(double));
   for(i=0;i<nRec;i++) s3_[i]=t_[i]*pow(2*M_PI*M_PI/45.*heff_[i],0.3333333333333333);
   heff_ln_diff_=malloc(nRec*sizeof(double));
   polintDiff(nRec,t_, heff_,  heff_ln_diff_);
   for(i=0;i<nRec;i++) heff_ln_diff_[i]*=t_[i]/heff_[i];
   Tdim=nRec;
   free(geff2_); geff2_=NULL;
   return nRec;
}

double gEff(double T)
{
  if(Tdim==0) loadHeffGeff(NULL);
  if(t_[0] < t_[Tdim-1])
  { if(T< t_[0]) T=t_[0];
    if(T> t_[Tdim-1]) T=t_[Tdim-1];
  } else
  { if(T> t_[0]) T=t_[0];
    if(T< t_[Tdim-1]) T=t_[Tdim-1];
  }

  return polint3(T,Tdim,t_,geff_);
}

double hEffLnDiff(double T)
{
  if(Tdim==0) loadHeffGeff(NULL);
  if(t_[0] < t_[Tdim-1])
  { if(T< t_[0]) T=t_[0];
    if(T> t_[Tdim-1]) T=t_[Tdim-1];
  } else
  { if(T> t_[0]) T=t_[0];
    if(T< t_[Tdim-1]) T=t_[Tdim-1];
  }
  return polint3(T,Tdim,t_,heff_ln_diff_);
}


double hEff(double T)
{
  if(Tdim==0) loadHeffGeff(NULL);
  if(t_[0] < t_[Tdim-1])
  { if(T< t_[0]) T=t_[0];
    if(T> t_[Tdim-1]) T=t_[Tdim-1];
  } else
  { if(T> t_[0]) T=t_[0];
    if(T< t_[Tdim-1]) T=t_[Tdim-1];
  }
  return polint3(T,Tdim,t_,heff_);
}

double T_s3(double s3) {if(s3>s3_[Tdim-1]) return t_[Tdim-1]*(s3/s3_[Tdim-1]); else  return polint3(s3,Tdim,s3_,t_);}
double s3_T(double T)  { if(T> t_[Tdim-1]) return s3_[Tdim-1]*T/t_[Tdim-1];    else  return polint3(T,Tdim,t_,s3_);}


static void  geffDeriv( double T, double *g, double *dg)
{   double eps=T*0.01;
    *dg= 4*(hEff(T)- (*g))/T + 4./3.* (hEff(T+eps) - hEff(T-eps))/(2*eps);
}

double gEff2(double T)
{
  if(geff2_==NULL)
  {  geff2_=malloc(Tdim*sizeof(double));
     for(int n=0; n<Tdim; n++)
     {
        int m,m1;
        if(t_[0] < t_[Tdim-1]) { m=n; m1=m-1;} else { m=Tdim-n-1; m1=m+1;}
        if( t_[m]<1E-3) geff2_[m]=geff_[m];
        else
        { double g=geff2_[m1];
          odeint( &g, 1, t_[m1],t_[m], 1E-4, (t_[m]-t_[m1])/2, geffDeriv);
          geff2_[m]=g;
        }
     }
/*
printf("==================\n");
for(int i=0;i<Tdim;i++) printf(" %E\n",geff2_[i]);
printf("=====================\n");
*/
  }
  return polint3(T,Tdim,t_,geff2_);
}
