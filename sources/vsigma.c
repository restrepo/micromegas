#include <sys/utsname.h>
#include "micromegas.h"
#include "micromegas_aux.h"
#include "micromegas_f.h"

//======================  vSigmaCC =======================
#define LATER     // it is about 2->4. I do not remember why it was commented out and restore it for version 5.4.1 

static double T_;
static double M1,M2,sqrtSmin,sqrtSmax;
static CalcHEP_interface * CI;
static int i3=2,i4=3,i5=4,i6=5;
static REAL pmass[6];
static int pdg[6];

static double sing2(char *s, int nout, double m, double w)
{  int i;
   int nin=2;
   
   double sum0=0,sum1=0;
   if(s[0]<=nin) return 0;
   if(strlen(s)<2) return 0;

   for(i=0;s[i];i++) sum1+=pmass[s[i]-1];
   if(m<=sum1) return 0;
   for(i=nin;i<nin+nout;i++) sum0+=pmass[i];      
   if(m>= sqrtSmax -sum0 +sum1) return 0;   
   return 1/(m*w);     
}

static double s_integrandT_(double  sqrtS )
{  double sv_tot,t,bess, x1,x2,y;
   REAL ms,md,PcmIn;
   double Rm;
   
   ms = M1 + M2; 
   if(ms>=sqrtS)  return 0;
   x1=M1/T_; x2=M2/T_; y=sqrtS/T_;

   if(y-x1-x2>50) return 0;   
      
   md = M1 - M2;
   PcmIn = Sqrt((sqrtS-ms)*(sqrtS+ms)*(sqrtS-md)*(sqrtS+md))/(2*sqrtS);
   kin22(PcmIn,pmass);
   
   sv_tot=simpson(dSigma_dCos,-1.,1.,1E-3,NULL); 
   improveCrossSection(pdg[0],pdg[1],pdg[2],pdg[3],PcmIn,&sv_tot);
   bess=    sqrt(2*x1*x2/y/M_PI)*exp(-(y-x1-x2))*K1pol(1/y)/(K2pol(1/x1)*K2pol(1/x2));
   Rm=PcmIn*sqrtS/M1/M2;   
   return  bess*Rm*Rm*sv_tot/T_;    
}   

/*
bessK2(x) = exp(-x)*sqrt(M_PI/2/x)*K2pol(1/x)
bessK1(x) = exp(-x)*sqrt(M_PI/2/x)*K1pol(1/x) 
*/

static double u_integrand_( double u)
{  double z,y,sv_tot,w;
   double Xf_1;
   double ms,md,sqrtS,PcmIn,res0;
   
   if(u==0. || u==1.) return 0.;
   z=1-u*u;
   sqrtS=M1+M2-3*T_*log(z);
   if(sqrtS<=M1+M2 || sqrtS<=pmass[2]+pmass[3]) return 0;
   return s_integrandT_(sqrtS )*6*T_*u/z;
}

static double vsigma23integrandT(double *x, double w)
{
   double pcmIn,sqrtS,M45;
   double M45_min=pmass[i4]+pmass[i5],M45_max=sqrtSmax-pmass[i3];
   int err;
   double r, x1,x2,y,z, bess,Rm;
   double GG;
   REAL pvect[20];

   z=1-x[0]*x[0];
   sqrtS=M1+M2-3*T_*log(z);
   
   if(sqrtS<=sqrtSmin || sqrtS>=sqrtSmax) return 0;    
   pcmIn=decayPcm(sqrtS,pmass[0], pmass[1]);

   M45=M45_min+x[1]*(M45_max-M45_min);
   
   x1=M1/T_; x2=M2/T_; y=sqrtS/T_;
   if(y-x1-x2>50) return 0;

   r=kinematic_23(pcmIn,i3, M45, 2*x[2]-1 ,2*x[3]-1,M_PI*x[4],pmass,pvect)*8*M_PI*(M45_max-M45_min); //  /pcmIn
   if(r==0) return 0;
   GG=sqrt(4*M_PI*alphaQCD(sqrtS));
   r*= CI->sqme(1,GG, pvect,NULL,&err);
   bess= sqrt(2*x1*x2/y/M_PI)*exp(-(y-x1-x2))*K1pol(1/y)/(K2pol(1/x1)*K2pol(1/x2));
//   bess=bessK1(sqrtS/T_)/bessK2(M1/T_)/bessK2(M2/T_);

   Rm=sqrtS/M1/M2;   
   
   r*= pcmIn*bess*Rm*Rm*6*x[0]/z;
   return r; 
}

static double vsigma23integrand0(double *x, double w)
{
   double r,sqrtS=M1+M2, M45_min=pmass[i4]+pmass[i5],M45_max=M1+M2-pmass[i3];
   int err;
   double GG;
   REAL pvect[20];
   
   r=kinematic_23(0.,i3,M45_min+x[0]*(M45_max-M45_min),0.5,2*x[1]-1,M_PI*x[2],pmass,pvect)*8*M_PI*(M45_max-M45_min);
   if(r==0) return 0;
   GG=sqrt(4*M_PI*alphaQCD(sqrtS));
   r*= CI->sqme(1,GG, pvect,NULL,&err);
   r*= (M1+M2)/M1/M2;
//printf("r=%e\n",r);   
   return r; 
}

#ifdef LATER
static double vsigma24integrandT(double *x, double w)
{
   double pcmIn,sqrtS,M34,M56;
   double M34_min=pmass[i3]+pmass[i4],M34_max=sqrtSmax-pmass[i5]-pmass[i6];
   double M56_min=pmass[i5]+pmass[i6],M56_max=sqrtSmax-pmass[i3]-pmass[i4];
   int err;
   double r, x1,x2,y,z, bess,Rm;
   double GG;
   REAL pvect[24];

   z=1-x[0]*x[0];
   sqrtS=M1+M2-3*T_*log(z);
   
   if(sqrtS<=sqrtSmin || sqrtS>=sqrtSmax) return 0;    
   pcmIn=decayPcm(sqrtS,pmass[0], pmass[1]);

   M34=M34_min+x[1]*(M34_max-M34_min);
   M56=M56_min+x[2]*(M56_max-M56_min);
   
   x1=M1/T_; x2=M2/T_; y=sqrtS/T_;
   if(y-x1-x2>50) return 0;

   r=kinematic_24(pcmIn,i3,i4, M34, M56,  2*x[3]-1 ,2*x[4]-1,2*x[5]-1, 2*M_PI*x[6], 2*M_PI*x[7],pmass,pvect)
                  *(M34_max-M34_min)*(M56_max-M56_min)*2*4*M_PI*4*M_PI; //  /pcmIn
   GG=sqrt(4*M_PI*alphaQCD(sqrtS));
   r*= CI->sqme(1,GG, pvect,NULL,&err);
   bess=    sqrt(2*x1*x2/y/M_PI)*exp(-(y-x1-x2))*K1pol(1/y)/(K2pol(1/x1)*K2pol(1/x2));
//   bess=bessK1(sqrtS/T_)/bessK2(M1/T_)/bessK2(M2/T_);

   Rm=sqrtS/M1/M2;   
   
   r*= pcmIn*bess*Rm*Rm*6*x[0]/z;
   return r; 
}


static double vsigma24integrand0(double *x, double w)
{
   double r,sqrtS=M1+M2,M34,M56;
   double  M34_min=pmass[i3]+pmass[i4],M34_max=M1+M2-pmass[i5]-pmass[i6],
           M56_min=pmass[i5]+pmass[i6],M56_max=M1+M2-pmass[i3]-pmass[i4];
   int err;
   double GG;
   REAL pvect[24];
   
   M34=M34_min+x[0]*(M34_max-M34_min);
   M56=M56_min+x[1]*(M56_max-M56_min);
   
   r=kinematic_24(0., i3,i4,M34, M56, 0.5 ,2*x[2]-1,2*x[3]-1, 2*M_PI*x[4], 2*M_PI*x[5],pmass,pvect)
                     *(M34_max-M34_min)*(M56_max-M56_min)*2*4*M_PI*4*M_PI; 
   
   GG=sqrt(4*M_PI*alphaQCD(sqrtS));
   r*= CI->sqme(1,GG, pvect,NULL,&err);
   r*= (M1+M2)/M1/M2;
//printf("r=%e\n",r);   
   return r; 
}
#endif

static double chFraction(double T, char*p0,char*p1)
{
  double M1=pMass(p0),M2=pMass(p1);
  char *c0=NULL,*c1=NULL;
  double s=0,k=M1*M1*M2*M2/sqrt(M1*M2)*K2pol(T/M1)*K2pol(T/M2);            
  for(int i=0;i<nModelParticles;i++) if(ModelPrtcls[i].name[0]=='~')
  {  double m=pMass(ModelPrtcls[i].name);
     int dim=ModelPrtcls[i].cdim*(ModelPrtcls[i].spin2+1);
         
          if(strcmp(p0,ModelPrtcls[i].name) ==0){ c0=ModelPrtcls[i].aname;k*=dim;} 
     else if(strcmp(p0,ModelPrtcls[i].aname)==0){ c0=ModelPrtcls[i].name; k*=dim;}
     
          if(strcmp(p1,ModelPrtcls[i].name) ==0){ c1=ModelPrtcls[i].aname;k*=dim;} 
     else if(strcmp(p1,ModelPrtcls[i].aname)==0){ c1=ModelPrtcls[i].name; k*=dim;}

         if(strcmp(ModelPrtcls[i].name,ModelPrtcls[i].aname)!=0) dim*=2;
         if(0.5*(M1+M2)-m >30*T){ k=0;s=1;break;}
         s+=dim*m*m/sqrt(m)*K2pol(T/m)*exp(-(m-0.5*(M1+M2))/T); 
      }
      if(k)      
      {  if(strcmp(p0,p1)) k*=2;
         if(! (  (strcmp(p0,c0)==0 && strcmp(p1,c1)==0)
               ||(strcmp(p0,c1)==0 && strcmp(p1,c0)==0)
              ) ) k*=2;   
      }     
//      printf("dmFactor %s(%s) %s(%s) =%e\n",p0,c0,p1,c1,k/s/s); 
      return k/s/s;
}



double vSigmaCC(double T,numout* cc, int mode)
{
  int i,err,n,n0,m,w;
  char*s, *pname[6];
  double msum;
  double a=0,factor,dMax=0;
  int spin2,cdim,neutral1,neutral2;
  double oldQ;
  
  double bEps=1.E-4;
  double dI;
  int dmOut;
      
  CI=cc->interface; 
  T_=T;
  


  if(passParameters(cc)) return -1;
  if(Qaddress && CI->nout==2) 
  {  oldQ=*Qaddress;
     for(i=0;i<2;i++) pname[i]=CI->pinf(1,i+1,pmass+i,pdg+i);
     *Qaddress=pmass[0]+pmass[1];
     calcMainFunc();
     if(passParameters(cc)) return -1;
  }
  
  for(i=0;i<2+CI->nout;i++) pname[i]=CI->pinf(1,i+1,pmass+i,pdg+i);  

  M1=pmass[0];
  M2=pmass[1];
  

  if(mode) 
  { if(pname[0][0]!='~' || pname[1][0]!='~') return 0;
    if(T==0 && (M1< Mcdm || M2<Mcdm)) return 0;
    dmOut=0; 
    for(i=2;i<2+CI->nout;i++) if(pname[i][0]=='~') dmOut++;
    if(dmOut==2) return 0; 
  }    
 
  for(i=2,msum=0;i<CI->nout;i++) msum+=pmass[i];
  
  sqrtSmin=M1+M2;
  
  if(msum > sqrtSmin)
  { if(T==0) return 0; else sqrtSmin=msum; }
  sqrtSmax=sqrtSmin-T*log(bEps); 

  n0=0; 
  if(CI->nout>2) for(n=1;(s=CI->den_info(1,n,&m,&w,NULL));n++)
  { double mm=0,ww=0;
    if(m) mm=fabs(CI->va[m]); if(w) ww=CI->va[w];
    double d=sing2(s,CI->nout,mm,ww); 
    if(!isfinite(d)) { printf("non-integrable pole\n"); return 0;}
    if(d>dMax){ dMax=d; n0=n;} 
  }

  switch(CI->nout)
  { 
     case 2:
       if(T==0) a=vcs22(cc,1,&err); else
       {  double eps=1.E-3;
          sqme22=CI->sqme;
          nsub22=1;
          a=simpson(u_integrand_,0.,1.,eps,NULL)*3.8937966E8;
       }   
       break;
     case 3:
     {  
        if(n0)
        {  s=CI->den_info(1,n0,&m,&w,NULL);
           for(i3=2;i3<5;i3++) if(i3!=s[0]-1 && i3!=s[1]-1) break;
           for(i4=2;i4<4;i4++) if(i4!=i3) break;   
           for(i5=i4+1;i5<=4;i5++) if(i5!=i3) break;    
        } else {i3=2;i4=3;i5=4;}
   
        
        if(T==0) a=vegas_chain(3, vsigma23integrand0 ,2000,1., 0.03,&dI,NULL);
        else     a=vegas_chain(5, vsigma23integrandT ,2000,1., 0.03,&dI,NULL);
        break;
     }
     case 4:
#ifdef LATER      
         if(n0) 
         {  s=CI->den_info(1,n0,&m,&w,NULL);
            i3=s[0]-1;
            i4=s[1]-1;
            for(i5=2;i5<5;i5++)  if(i5!=i3 && i5!=i4) break;
            for(i6=i5+1;i6<=5;i6++) if(i6!=i3 && i6!=i4) break;
         }else { i3=2;i4=3;i5=4;i6=5;} 
//         printf("i3,i4,i5,i6= %d %d %d %d\n", i3,i4,i5,i6); 
         if(T==0) a=vegas_chain(6, vsigma24integrand0 ,4000,1., 0.03,&dI,NULL);
         else     a=vegas_chain(8, vsigma24integrandT ,20000,1., 0.001,&dI,NULL);
#else 
    printf("2->4 is not implemented yet\n");
#endif         
                                
     break;
     default:
        printf("Too many outgoing particles\n");
        a=0;  
   }  

//   WIDTH_FOR_OMEGA=0;
  
   if(mode)
   {  a*= 1-0.5*dmOut;
      char *p0=pname[0],*p1=pname[1],*c0=NULL,*c1=NULL;
      double s=0,k=M1*M1*M2*M2/sqrt(M1*M2)*K2pol(T/M1)*K2pol(T/M2);            
      for(i=0;i<nModelParticles;i++) if(ModelPrtcls[i].name[0]=='~')
      {  double m=pMass(ModelPrtcls[i].name);
         int dim=ModelPrtcls[i].cdim*(ModelPrtcls[i].spin2+1);
         
              if(strcmp(p0,ModelPrtcls[i].name) ==0){ c0=ModelPrtcls[i].aname;k*=dim;} 
         else if(strcmp(p0,ModelPrtcls[i].aname)==0){ c0=ModelPrtcls[i].name; k*=dim;}
              if(strcmp(p1,ModelPrtcls[i].name) ==0){ c1=ModelPrtcls[i].aname;k*=dim;} 
         else if(strcmp(p1,ModelPrtcls[i].aname)==0){ c1=ModelPrtcls[i].name; k*=dim;}
         if(strcmp(ModelPrtcls[i].name,ModelPrtcls[i].aname)!=0) dim*=2;
         if(0.5*(M1+M2)-m >30*T){ k=0;s=1;break;}
         s+=dim*m*m/sqrt(m)*K2pol(T/m)*exp(-(m-0.5*(M1+M2))/T); 
      }
      if(k)      
      {  if(strcmp(p0,p1)) k*=2;
         if(! (  (strcmp(p0,c0)==0 && strcmp(p1,c1)==0)
               ||(strcmp(p0,c1)==0 && strcmp(p1,c0)==0)
              ) ) k*=2;   
      } 
      a*=k/s/s;
   }  

   if(Qaddress && CI->nout==2) { *Qaddress=oldQ; calcMainFunc();}
   
   return a;
}

double vsigmacc_(double *T, int *ccf,int*mode)
{
  numout*cc;
  memcpy(&cc,ccf,sizeof(cc));
  return vSigmaCC(*T,cc,*mode); 
}

//#define P_NAME_SIZE 11


int proc2names(char*process, int * in, int * out, ...)
{ 
   if(strlen(process)>100) return 1;
   char proc[100];
   strcpy(proc, process);
   char *ch1=strstr(proc,"->");
   if(ch1==NULL) return 2;
   ch1[0]=','; ch1[1]=' ';

   *in=0; *out=0;
   va_list ap;
   va_start(ap, out);
   int err=0;

   for(  char*ch=proc;;)
   { char buff[100];
     int r=sscanf(ch,"%[^ , ]",buff);
     if(r==0) break;
     ch+=strlen(buff);
     if(ch<=ch1) (*in)++; else (*out)++;
     while(ch[0]==' '||ch[0]==',') ch++;
     trim(buff);
     char*p=va_arg(ap, char*);
     if(p==NULL) {err=3; break;}
     if(strlen(buff)>=P_NAME_SIZE) { err=4; break;} 
     strcpy(p,buff);
     if(ch[0]==0) break;
   }      
   va_end(ap);
   return err;        
} 

// ================  vSigmaPlus ==================

static numout*cc23,*cc24;
static int *l_23;
static double cs23_sum(double Pcm)
{ //if(Pcm==0) return 0;
  double cs=0;
  int err;
  for(int k=1;k<=cc23->interface->nprc;k++) cs+=cs23Pcm(cc23, k, Pcm, l_23[k-1] ,&err);
  return cs;
}  


typedef struct{ char proc22[100]; numout*cc22;  double T,M1,M2;  int Npcm; double*pcmArr; double *csArr; double sqrtSmin,sqrtSmin22,sqrtSmax;} tabulated2N; 


double tabVsInt(double sqrtS, void * par_)
{
  tabulated2N*par=par_;
  double M1=par->M1, M2=par->M2, T=par->T;
  
  double pcm=decayPcm(sqrtS,M1,M2);
  if(pcm==0) return 0;
  double cs=polint3(pcm, par->Npcm, par->pcmArr,par->csArr);
  
  double bess,x1,x2,y,Rm;
   
  x1=M1/T; x2=M2/T; y=sqrtS/T;

  if(y-x1-x2>100) return 0;   
 
  bess=    sqrt(2*x1*x2/y/M_PI)*exp(-(y-x1-x2))*K1pol(1/y)/(K2pol(1/x1)*K2pol(1/x2));
  Rm=sqrtS/M1/M2;  
  
  return  pcm*bess*Rm*Rm*cs/T;
}

double cc22Int(double sqrtS, void * par_)
{ int err;
  tabulated2N*par=par_;
  double M1=par->M1, M2=par->M2, T=par->T;
  
  double pcm=decayPcm(sqrtS,M1,M2);
  if(pcm==0) return 0;
  double cs=cs22(par->cc22,1,pcm,-1,1,&err);
//printf("pcm=%E cs=%E\n",pcm,cs);  
  double bess,x1,x2,y,Rm;
   
  x1=M1/T; x2=M2/T; y=sqrtS/T;

//  if(y-x1-x2>100) return 0;   
 
  bess=    sqrt(2*x1*x2/y/M_PI)*exp(-(y-x1-x2))*K1pol(1/y)/(K2pol(1/x1)*K2pol(1/x2));
  Rm=sqrtS/M1/M2;  
     
  return  pcm*bess*Rm*Rm*cs*pcm/T;
}


double  vSigmaPlus23(char* proc22,  double T, int *Err)
{   
   static int myTime=0;
   static tabulated2N *all23=NULL;
   static int N23proc=0;
   if(myTime<sortOddTime) 
   {  myTime=sortOddTime;
      if(all23)
      { for(int k=0;k<N23proc;k++)  { free(all23[k].pcmArr); free(all23[k].csArr);}
        free(all23);
        all23=NULL;
        N23proc=0; 
      } 
   }

   int in,out,err;
   char p[4][P_NAME_SIZE],q[3][P_NAME_SIZE];   
   err=proc2names(proc22, &in, &out, p[0],p[1],p[2],p[3],NULL);
   if(err) { if(Err) *Err=err; return NAN;}
   
   if(p[0][0]!='~' || p[1][0]!='~') { if(Err) *Err=4; return NAN;}
   
   for(int k=0;k<N23proc;k++) if(strcmp(proc22,all23[k].proc22)==0) 
   {  all23[k].T=T;
      double Mmax=all23[k].sqrtSmax;
      if(Mmax>all23[k].sqrtSmin+20*T) Mmax=all23[k].sqrtSmin+20*T; 
      double r= simpson_arg(tabVsInt, all23+k, all23[k].sqrtSmin, Mmax,1E-3,NULL);
      if(all23[k].sqrtSmin22<Mmax) r-= simpson_arg( cc22Int, all23+k, all23[k].sqrtSmin22, Mmax,1E-3,NULL);
      return  r*chFraction(T,p[0],p[1]);
   }
   char lib[100];
   err=process2Lib(proc22,lib);
   if(err) { if(Err) *Err=5;  return NAN;}
   
   int l=2,l_=3;
   txtList lL,l_L;
   double  w=pWidth(p[l],&lL), w_=pWidth(p[l_],&l_L); 
   if(w==0 && w_==0) return 0;
//   if(w_>w) 
    if(w_<w)
   {  double m=w;w=w_;w_=m; l=5-l;l_=5-l_;
      void *mL=lL; lL=l_L; l_L=mL;
   }  
      
   char proc23[100];
   char A[100]="",B[100]="";
   
   for(txtList L=lL; L; L=L->next)
   {  
     sscanf(L->txt,"%*f %s -> %[^,], %[^,], %[^,]",q[0],q[1],q[2],q[3]);
//     printf("L->txt=%s  err=%d q =  %s %s %s\n",L->txt,q[0],q[1],q[2]); 
     strcat(q[1],",");  strcat(q[2],","); 
     { char*ch=strstr(A,q[1]);
       if(ch && (ch!=A && ch[-1]!=',')) ch=NULL; 
       if(ch==NULL)  strcat(A,q[1]); 
       
       ch=strstr(B,q[2]);
       if(ch && ch!=B && ch[-1]!=',') ch=NULL;    
       if(ch==NULL) strcat(B,q[2]); 
     }
   }  

   A[strlen(A)-1]=0;  
   B[strlen(B)-1]=0; 
   
   sprintf(proc23,"%s,%s->%s,_A_,_B_{%s{%s", p[0],p[1],p[l_],A,B);
   char excludeV[P_NAME_SIZE+20], excludeO[P_NAME_SIZE+20];
   sprintf(excludeV,"%s<1",p[l]);
   if(strcmp(p[l],p[l_])==0) sprintf(excludeO,"%s>1",p[l]); else sprintf(excludeO,"%s",p[l]); 
   
   int n=pTabPos(p[l]);
   if(n>0) sprintf(lib+strlen(lib),"_P%d",n); else sprintf(lib+strlen(lib),"_A%d",n);
   cc23=getMEcode(0,ForceUG, proc23, excludeV, excludeO,lib);

   if(cc23==NULL) return 0;
   passParameters(cc23);

   double  M1=pMass(p[0]),M2=pMass(p[1]),  M=pMass(p[l]), M_=pMass(p[l_]);
   double sqrtSmin22=M+M_; sqrtSmax=M+10*w+M_+10*w_, sqrtSmin=sqrtSmax;

      
   l_23=malloc(cc23->interface->nprc*sizeof(int));   
   for(int k=1;k<=cc23->interface->nprc;k++)
   {
      REAL m[5];
      for(int i=2;i<5;i++) { char* p23=cc23->interface->pinf(k,i+1,m+i,NULL); if(strcmp(p23,p[l_])==0) l_23[k-1]=i;}
      double Ms=m[2]+m[3]+m[4];
      if(Ms<sqrtSmin) sqrtSmin=Ms;
   }  

   if(sqrtSmin  <=M1+M2) sqrtSmin=  (M1+M2);  
   if(sqrtSmin22<=M1+M2) sqrtSmin22=(M1+M2); //sqrtSmin22*=1.00000001;

   if(sqrtSmax<=sqrtSmin) {free(l_23); return 0;}
    
   double pcmMin=decayPcm(sqrtSmin,M1,M2), pcmMax=decayPcm(sqrtSmax,M1,M2);
  
   int Npcm;
   double *pcmArr,*csArr;  
  
   err= buildInterpolation( cs23_sum, pcmMin , pcmMax, 
               -0.01, 0.01, &Npcm, &pcmArr, &csArr);
    
   polintStr Arg;
   Arg.dim=Npcm; 
   Arg.x=pcmArr; 
   Arg.y=csArr;
//   displayPlot("cs_23","pcm",pcmMin , 3, 0, 2, "cs23",0, cs23_sum,NULL, "interpalation",0,polint_arg,&Arg);  
   free(l_23);
       
   numout*cc22=newProcess(proc22);
   passParameters(cc22);  
   for(int i=0;i<Npcm;i++)
   {  double pcm=pcmArr[i]; 
      double sqrtS=sqrt(M1*M1+ pcm*pcm)+sqrt(M2*M2+ pcm*pcm);
      if(w_>0)
      { csArr[i]*=decayPcmW(sqrtS, M,M_,w,w_,0)/decayPcmW(sqrtS,M,M_,w,0,0);
        if(strcmp(p[l],p[l_])==0) csArr[i]/=2;
      }    
   }  
   all23=realloc(all23, (N23proc+1)*sizeof(tabulated2N));
   strcpy(all23[N23proc].proc22,proc22);
   all23[N23proc].cc22=cc22;
   all23[N23proc].sqrtSmin=sqrtSmin;
   all23[N23proc].sqrtSmax=sqrtSmax;
   all23[N23proc].sqrtSmin22=sqrtSmin22;
   all23[N23proc].pcmArr=pcmArr;
   all23[N23proc].csArr=csArr;
   all23[N23proc].Npcm=Npcm;
   all23[N23proc].M1=M1;
   all23[N23proc].M2=M2;
   all23[N23proc].T=T;
   
   int k=N23proc;
   N23proc++;
   
   double Mmax=all23[k].sqrtSmax;
   if(Mmax>sqrtSmin+20*T) Mmax=sqrtSmin+20*T; 
   double r=               simpson_arg( tabVsInt, all23+k, sqrtSmin, Mmax,1E-3,NULL);
   if(Mmax>sqrtSmin22) r-= simpson_arg( cc22Int, all23+k, sqrtSmin22,Mmax,1E-3,NULL);
   return  r*chFraction(T,p[0],p[1]);
}


double  vSigmaPlus24(char* proc22,  double T, int *Err)
{   

//printf("vSigmaPlus24 T=%e\n",T);
   static int myTime=0;
   static tabulated2N *all24=NULL;
   static int N24proc=0;
   if(myTime<sortOddTime) 
   {  myTime=sortOddTime;
      if(all24)
      { for(int k=0;k<N24proc;k++)  { free(all24[k].pcmArr); free(all24[k].csArr);}
        free(all24);
        all24=NULL;
        N24proc=0; 
      } 
   }
  
   int in,out,err;
   char p[4][P_NAME_SIZE],q3[4][P_NAME_SIZE],q4[4][P_NAME_SIZE];
      
   err=proc2names(proc22, &in, &out, p[0],p[1],p[2],p[3],NULL);
   if(err) { if(Err) *Err=err; return NAN;}
   
   if(p[0][0]!='~' || p[1][0]!='~') { if(Err) *Err=4; return NAN;}
   
   for(int k=0;k<N24proc;k++) if(strcmp(proc22,all24[k].proc22)==0) 
   {  all24[k].T=T;

      double Mmax=all24[k].sqrtSmax;
      if(Mmax>all24[k].sqrtSmin+20*T) Mmax=all24[k].sqrtSmin+20*T; 
      double r2=0,r1=                       simpson_arg(tabVsInt, all24+k, all24[k].sqrtSmin,  Mmax,1E-3,NULL);
      if(all24[k].sqrtSmin22<Mmax) r2=simpson_arg( cc22Int, all24+k, all24[k].sqrtSmin22,Mmax,1E-3,NULL);
double      r=r1-r2; 
/*
if(T<2E-3)
{ char TT[20];
  sprintf(TT,"T=%E r1=%E r=%E sqrtSmin=%E %E   ",T,r1,r, all24[k].sqrtSmin, all24[k].sqrtSmin22);
  displayPlot(TT,"sqrt(s)", sqrtSmin, Mmax, 0,2,"tabVsInt",0, tabVsInt, all24+k, "cc22Int",0, cc22Int, all24+k);   
}
*/
      return r*chFraction(T,p[0],p[1]);
   }
   txtList L3, L4; 
   double  w3=pWidth(p[2],&L3), w4=pWidth(p[3],&L4),b3=0,b4=0; 
   if(w3==0 || w4==0) return 0;
   
   double M1=pMass(p[0]),M2=pMass(p[1]), sqrtSmin22=pMass(p[2])+pMass(p[3]), M3,M4;
   char proc24[100];      
   sprintf(proc24,"%s,%s->", p[0],p[1]);

//printf("proc=%s\n", proc24);
   
   char*ch=proc24+strlen(proc24);   
   for(txtList L=L3;L;L=L->next)
   {  double b;
      sscanf(L->txt,"%lf",&b);
      if(b>b3)
      {
        int n=sscanf(L->txt,"%lf %s -> %[^,], %[^,], %[^,]",&b,q3[0],q3[1],q3[2],q3[3]);
//printf("n=%d\n",n);        
        if(n==4) 
        { sprintf(ch,"%s,%s", q3[1],q3[2]); 
          b3=b; 
          M3=pMass(q3[1])+pMass(q3[2]);
        }     
      } 
   }
   if(b3==0) { if(Err) *Err=5; return NAN;} 

//printf("    proc=%s\n", proc24);
   
   ch=proc24+strlen(proc24);
   for(txtList L=L4;L;L=L->next)
   { double b;
     sscanf(L->txt,"%lf",&b);
     if(b>b4)
     {
        int n=sscanf(L->txt,"%lf %s -> %[^,], %[^,], %[^,]",&b,q4[0],q4[1],q4[2],q4[3]);
        if(n==4 && strcmp(q3[1],q4[1]) && strcmp(q3[2],q4[2]) && strcmp(q3[1],q4[2])&& strcmp(q3[2],q4[1]) ) 
        { sprintf(ch,",%s,%s",q4[1],q4[2]); 
          b4=b; 
          M4=pMass(q4[1])+pMass(q4[2]);
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
   
   if(sqrtSmin22<M1+M2) sqrtSmin22=M1+M2; //sqrtSmin22*=1.00001;
   double  sqrtSmax=sqrtSmin22+10*(w3+w4), sqrtSmin=M1+M2;
   if(sqrtSmin<M3+M4) sqrtSmin=M3+M4;
   if(sqrtSmin<=M1+M2) sqrtSmin=(M1+M2);
//printf("sqrtSmax=%E sqrtSmin=%E\n", sqrtSmax, sqrtSmin);   
   if(sqrtSmax<=sqrtSmin) { return 0;}

   cc24=getMEcode(0,ForceUG, proc24, excludeV, NULL,lib);
//   printf("cc24=%p proc23=|%s| excludeV=|%s| lib=|%s|\n",cc24, proc24, excludeV, lib);

   if(cc24==NULL) return 0;
   passParameters(cc24);
      
    
   double pcmMin=decayPcm(sqrtSmin*1.0001,M1,M2), pcmMax=decayPcm(sqrtSmax,M1,M2);
   int i3=2,i4=3;
   for(i3=2;     strcmp(q3[1], cc24->interface->pinf(1,i3+1,NULL,NULL)); i3++) /*printf("%d %s %s\n", i3, q3[1], cc24->interface->pinf(1,i3+1,NULL,NULL))*/  ;
   for(i4=i3+1;  strcmp(q3[2], cc24->interface->pinf(1,i4+1,NULL,NULL)); i4++);

//printf("i3=%d i4=%d q3[1]=%s q3[2]=%s  \n", i3,i4,q3[1],q3[2]);
   
//   printf(" %s %s\n", cc24->interface->pinf(1,i3+1,NULL,NULL), cc24->interface->pinf(1,i4+1,NULL,NULL) );

   int Npcm=20;
   double *pcmArr=malloc(sizeof(double)*Npcm) ,*csArr=malloc(sizeof(double)*Npcm);  
      
   for(int i=0;i<Npcm;i++) 
   {  double dcs,chi2;
      int N;
      pcmArr[i]=pcmMin+ i*(pcmMax-pcmMin)/(Npcm-1.); 
      csArr[i]=pcmArr[i]*cs24Vegas(cc24, 1, pcmArr[i], i3, i4, 
       2 , 100000,  2, 100000, 
      &dcs, &chi2, &N)/b3/b4;
//      printf(" pcm=%.2E   cs=%.2E dcs=%.2E chi2=%.2E err=%d\n", pcmArr[i], csArr[i], dcs, chi2, err);   
   }
   
   if(strcmp(p[2],p[3])==0) for(int i=0;i< Npcm;i++) csArr[i]/=2;  

   numout*cc22=newProcess(proc22);
   passParameters(cc22);  
   double Ms=pMass(p[2])+pMass(p[3]);
   
//printf("M1=%E M2=%E Ms=%E\n", M1,M2,Ms);   
/*
   for(int i=0;i<Npcm;i++)
   {  double pcm=pcmArr[i]; 
      double sqrtS=sqrt(M1*M1+ pcm*pcm)+sqrt(M2*M2+ pcm*pcm);   
      if(sqrtS>Ms)   csArr[i]-=pcm*cs22(cc22,1, pcmArr[i],-1,1,NULL);
      printf("sqrtS=%E      pcm=%.2E   cs=%.2E \n", sqrtS, pcmArr[i], csArr[i]);   
   }
*/
   
//   printf("csArr[Npcm-1]=%e %e\n", csArr[Npcm-1], pcmArr[Npcm-1]*cs22(cc22,1, pcmArr[Npcm-1],-1,1,NULL));
//   printf("dmFactor %s %s =%e\n", p[0],p[1], chFraction(arrT[40],p[0],p[1]));
   
   all24=realloc(all24, (N24proc+1)*sizeof(tabulated2N));
   strcpy(all24[N24proc].proc22,proc22);
   all24[N24proc].cc22=cc22;
   all24[N24proc].sqrtSmin=sqrtSmin;
   all24[N24proc].sqrtSmax=sqrtSmax;
   all24[N24proc].sqrtSmin22=sqrtSmin22;
   all24[N24proc].pcmArr=pcmArr;
   all24[N24proc].csArr=csArr;
   all24[N24proc].Npcm=Npcm;
   all24[N24proc].M1=M1;
   all24[N24proc].M2=M2;
   all24[N24proc].T=T;
   
   int k=N24proc;
   N24proc++;
   
   double Mmax=sqrtSmax; 
   if(Mmax > sqrtSmin+20*T) Mmax=sqrtSmin+20*T; 
   double r= simpson_arg( tabVsInt, all24+k, sqrtSmin, Mmax,1E-3,NULL);
   if(sqrtSmin22<Mmax) r-=simpson_arg( cc22Int,  all24+k, sqrtSmin22, Mmax,1E-3,NULL);


   return  r*chFraction(T,p[0],p[1]);
}


