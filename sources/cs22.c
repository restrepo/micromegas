#include "micromegas.h"
#include "micromegas_aux.h"
#include "micromegas_f.h"
#include"../CalcHEP_src/c_source/ntools/include/vegas.h"

double (*sqme22)(int nsub, double GG, REAL *pvect, REAL*cb_coeff, int * err_code)=NULL; 
int  nsub22=0;

/*===========================================================*/
static double Q_ren,Q_fact;
static double GG=1.23;
static double PcmOut, totcoef;
static REAL pvect[20]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};


static double eps=0.001;

double GGscale=91.187;
/*
double  decayPcm(double am0,  double  am1,  double  am2)
{
  double  summ, diffm, pout;
  summ = am1 + am2;
  diffm = am1 - am2;
  if(am0<summ) return 0;
  return sqrt((am0-summ)*(am0+ summ)*(am0-diffm)*(am0+diffm))/(am0*2);
}
*/          

int  kin22(double PcmIn,REAL * pmass)
{  
   double sqrtS;
   int i;
   for(i=0;i<16;i++) pvect[i]=0;
   sqrtS=sqrt(pmass[0]*pmass[0]+PcmIn*PcmIn)+sqrt(pmass[1]*pmass[1]+PcmIn*PcmIn);
   PcmOut = decayPcm(sqrtS,pmass[2],pmass[3]);
//printf(" PcmOut =%E (%E %E %E) \n", PcmOut,sqrtS,pmass[2],pmass[3] );   
   if(PcmOut<sqrtS*1.E-4) return 1;
   totcoef = PcmOut /(32.0*M_PI*PcmIn*sqrtS*sqrtS);
   pvect[3] = PcmIn;
   pvect[7] =-PcmIn;
   pvect[0] = sqrt(PcmIn*PcmIn   + pmass[0]*pmass[0]);
   pvect[4] = sqrt(PcmIn*PcmIn   + pmass[1]*pmass[1]);
   pvect[8] = sqrt(PcmOut*PcmOut + pmass[2]*pmass[2]);
   pvect[12]= sqrt(PcmOut*PcmOut + pmass[3]*pmass[3]);

   return 0;
}

double  dSigma_dCos(double  cos_f)
{
   double  r;
   double sin_f=sqrt(fabs((1-cos_f)*(1+cos_f)));
   int err_code=0;
   
   pvect[11]=PcmOut*cos_f;
   pvect[15]=-pvect[11];
   pvect[10]=PcmOut*sin_f;
   pvect[14]=-pvect[10];
   
   
   r = (*sqme22)(nsub22,sqrt(4*M_PI*parton_alpha(GGscale)),pvect,NULL,&err_code);
   err_code=0;
   return r * totcoef;
}


double cs22(numout * cc, int nsub, double P, double cos1, double cos2 , int * err) 
{
  int i,k;
  REAL pmass[4];
  GG=sqrt(4*M_PI*parton_alpha(GGscale));
  passParameters(cc);
  
  *(cc->interface->gtwidth)=0;
  *(cc->interface->twidth)=0;
  *(cc->interface->gswidth)=0;
  
  for(i=0;i<4;i++) cc->interface->pinf(nsub,1+i,pmass+i,NULL);  
  *err=0;
  sqme22=cc->interface->sqme;
  nsub22=nsub; 
  if(kin22(P,pmass)) return 0; else return 3.8937966E8*simpson(dSigma_dCos,cos1,cos2,0.3*eps);
}

/*===================  Collider production ==========*/

static double sMin,sMax,pcmOut;
static REAL pmass[5];
static int pc1_,pc2_;
static int ppFlag;
static int i3,i4,i5;
static double pTmin_;

static numout * colliderProduction(char * name1,char *name2, int nf, int J)
{ 
  char libname[100], process[100], lName1[20], lName2[20];
  numout *cc;
  int i,first;

  
  if(name1==NULL && name2==NULL) return NULL;
  if(name1) pname2lib(name1,lName1);
  if(name2) pname2lib(name2,lName2);
  if(name1==NULL)
  {
    sprintf(libname,"PP_%s",lName2);
    sprintf(process,"proton,proton->%s",name2);
  } else if(name2==NULL)
  { 
    sprintf(libname,"PP_%s",lName1);
    sprintf(process,"proton,proton->%s",name1); 
  } else
  {               
     if(strcmp(lName1,lName2)>0)sprintf(libname,"PP_%s%s",lName1,lName2);
     else                       sprintf(libname,"PP_%s%s",lName2,lName1); 
     sprintf(process,"proton,proton->%s,%s",name1,name2);
  }
     
  sprintf(libname+strlen(libname),"_nf%d",nf);
 
  if(J)
  {  sprintf(process+strlen(process),",proton");
     strcat(libname,"_J");
  }
  
  sprintf(process+strlen(process),"{");
  
  for(i=0,first=1;i<nModelParticles;i++) 
  { int pdg=abs(ModelPrtcls[i].NPDG);
    if(pdg && (pdg==21 ||(nf>=2 &&  pdg<=nf)|| (nf==1 && pdg==2) )) 
    { 
       if(!first) strcat(process,","); else  first=0;
       sprintf(process+strlen(process),"%s",ModelPrtcls[i].name);
       if(strcmp(ModelPrtcls[i].name,ModelPrtcls[i].aname))
       sprintf(process+strlen(process),",%s",ModelPrtcls[i].aname);
    }                              
  }
  
  cc=getMEcode(0,ForceUG,process,NULL,"",libname);
  
  return cc;
}


static double  cos_integrand(double xcos)
{ int err;
  double xsin=sqrt(1-xcos*xcos);
  double q;
  pvect[9]=pcmOut*xcos;
  pvect[10]=pcmOut*xsin;
  pvect[13]=-pvect[9];
  pvect[14]=-pvect[10];
  q=Q_ren>0? Q_ren: pvect[0]+pvect[4];  
  return  sqme22(nsub22,sqrt(4*M_PI*parton_alpha(q)),pvect,NULL,&err);  
}


static double  s_integrand(double y)
{  double r,pcmIn,q;

   double pp=1.5;
   
   double  s=sMin*pow(1+ y*( pow(sMax/sMin,1-pp) -1) ,1/(1-pp));
   double x0=s/sMax;
   
   pcmIn=decayPcm(sqrt(s),pmass[0], pmass[1]);
   if(pcmIn==0) return 0;
   pvect[0]=sqrt(pmass[0]*pmass[0]+pcmIn*pcmIn);
   pvect[1]=pcmIn; pvect[2]=0; pvect[3]=0;
   pvect[4]=sqrt(pmass[1]*pmass[1]+pcmIn*pcmIn);
   pvect[5]=-pcmIn; pvect[6]=0; pvect[7]=0;
   pcmOut=decayPcm(sqrt(s),pmass[2], pmass[3]);
   pvect[8]=sqrt(pmass[2]*pmass[2]+pcmOut*pcmOut);
   pvect[11]=0;
   pvect[12]=sqrt(pmass[3]*pmass[3]+pcmOut*pcmOut);
   pvect[15]=0;

   q=Q_fact>0? Q_fact: sqrt(s);
   r=  3.8937966E8*pcmOut/(32*M_PI*pcmIn*s)*simpson(cos_integrand,-1.,1.,1.E-3);
   r*=convStrFun2(x0,q,pc1_,pc2_,ppFlag);
   r*=pow(s/sMax,pp)*(1- pow(sMin/sMax,1-pp))/(1-pp);
   return r; 
}


static double M45_min,M45_max,S34_min,S34_max,S35_min,S35_max;
static int npole34=0,npole35=0,npole45,npole12;
static double * pole34=NULL,*pole35=NULL,*pole45=NULL,*pole12=NULL;


static double veg_intergrand(double *x, double w)
{
   double r;
//   double pp=1.5;
   double M12,pcmIn;
   int err;
   double P,M45,m3q,sn_,cs_,J45,cs45,S34,S35;
   double M34_min=sqrt(S34_min),M34_max=sqrt(S34_max),
          M35_min=sqrt(S35_min),M35_max=sqrt(S35_max);
   double M34,M35;
          
   
   M12=sqrt(sMin)+x[0]*sqrt(sMax); 
  
   pcmIn=decayPcm(M12,pmass[0], pmass[1]);
   if(pcmIn==0) return 0;
 
   M45=M45_min+x[1]*(M45_max-M45_min);
   
   if(M12<=pmass[i3]+M45) return 0;
   P=decayPcm(M12,pmass[i3], M45); 
   if(P<=pTmin_) return 0; 
   sn_=pTmin_/P;

   cs_=sqrt(1-sn_*sn_);
   if(npole34==0 && npole35==0)  { cs45=(2*x[3]-1); J45=2;}  else 
   {  
      double E3= 0.5*(M12*M12-M45*M45-pmass[i3]*pmass[i3])/M45; 
      double p= sqrt(E3*E3-pmass[i3]*pmass[i3]);
      double pcm2=decayPcm(M45,pmass[i4],pmass[i5]);
      double E4=sqrt(pcm2*pcm2+pmass[i4]*pmass[i4]);
      double E5=sqrt(pcm2*pcm2+pmass[i4]*pmass[i5]);

      if(npole35==0)
      { 
         M34=M34_min+x[3]*(M34_max-M34_min);
         cs45=(M34*M34-pmass[i4]*pmass[i4]-pmass[i3]*pmass[i3]-2*E3*E4)/(2*p*pcm2);
         J45=M34*(M34_max-M34_min)/(p*pcm2);
      } else if(npole34==0)
      { 
         M35=M35_min+x[3]*(M35_max-M35_min);
         cs45=-(M35*M35 -pmass[i5]*pmass[i5]-pmass[i3]*pmass[i3]-2*E3*E4)/(2*p*pcm2);
         J45=M35*(M34_max-M34_min)/(p*pcm2);
      } else 
      {  int i;
         double pr34=0,pr35=0;
         if(x[3]<0.5)
         {   M34=M34_min+2*x[3]*(M34_max-M34_min);
             S34=M34*M34;  
             cs45=(S34 -pmass[i4]*pmass[i4]-pmass[i3]*pmass[i3]-2*E3*E4)/(2*p*pcm2);
             S35=pmass[i5]*pmass[i5]+pmass[i3]*pmass[i3]+2*E3*E5 -2*p*pcm2*cs45;
             J45=2*M34*(M34_max-M34_min)/(p*pcm2);
         } else 
         {   M35=M35_min+(2*x[3]-1)*(M35_max-M35_min);
             S35=M35*M35;  
             cs45=-(S35 -pmass[i5]*pmass[i5]-pmass[i3]*pmass[i3]-2*E3*E5)/(2*p*pcm2);
             S34=pmass[i4]*pmass[i4]+pmass[i3]*pmass[i3]+2*E3*E4 +2*p*pcm2*cs45; 
             J45=2*M35*(M35_max-M35_min)/(p*pcm2);
         } 
         for(i=0;i<npole34;i++)
         { double m=pole34[2*i],w=pole34[2*i+1];
           pr34+=  M35*m*w/( (m*m-S34)*(m*m-S34) +m*m*w*w);
         }
         pr34+=1/(M34_max-M34_min); 
           
         for(i=0;i<npole35;i++)
         { double m=pole35[2*i],w=pole35[2*i+1];
           pr35+=  M35*m*w/( (m*m-S35)*(m*m-S35) +m*m*w*w);
         } 
         pr35+=1/(M35_max-M35_min); 
         if(x[3]<0.5) J45*=pr34/(pr34+pr35); else  J45*=pr35/(pr34+pr35);  
      }    
   } 
   if(fabs(cs45)>1) return 0;
   r=  kinematic_23(pcmIn,i3,M45, cs_*(2*x[2]-1) ,cs45,M_PI*x[4],pmass, pvect)*4*M_PI*(M45_max-M45_min)*J45*cs_/pcmIn;
   {   double q= Q_ren>0? Q_ren : M12;
       double x0=M12*M12/sMax;
      
       r*= sqme22(nsub22,sqrt(4*M_PI*parton_alpha(q)),pvect,NULL,&err);
       q= Q_fact>0? Q_fact: M12;
       r*= convStrFun2(x0,q,pc1_,pc2_,ppFlag);
   }
   r*= 2*M12/sMax*(sqrt(sMax) - sqrt(sMin));

   return r; 
}

   
static  void impGrid(int N, double*x, double *y)
{ int i;
  double s,rc;
  double alph=1.5;
  double X[100],r[100];

  for(s=0, i=0; i<N; i++) s+=y[i];
  for(rc=0,i=0; i<N; i++)
  {
     r[i] = 0;
     if( y[i]  > 0)
     {  double xoln = log(s/y[i]);
           if(xoln<0.01)         r[i]=1;
           else if (xoln <= 70)  r[i] = pow( (1 - exp(-xoln))/xoln, alph);
           else                  r[i] = pow(  1/xoln,               alph);
     }
     rc += r[i];  
  }
  rc /= N;
  if(rc)
  {  double dr=0,xo=0,xn=0;
     int k=0;
     for(i=1;i<N;) 
     {  for(;dr<rc;k++) dr+=r[k];
        xo=x[k-1];
        xn=x[k];
        for(;dr>=rc; i++) {dr-=rc; X[i] = xn-(xn-xo)*(x[N]-x[0])*dr/r[k-1];}
     }
     for(i=1;i<N;i++) x[i]=X[i];
  }
}


static void setGrid(int N, double*x, double x0, double x1,double M0,double M1,int Npole,double *poles)
{ 
  int i,k,nStep=10; 
  double *y=malloc(N*sizeof(double));
  for(i=0;i<=N;i++) x[i]=x0+i*(x1-x0)/N;
  for(k=0;k<nStep;k++)
  {  int l;
     for(i=0;i<N;i++)
     {  double m1=M0+(x[i]  -x0)*(M1-M0)/(x1-x0);
        double m2=M0+(x[i+1]-x0)*(M1-M0)/(x1-x0); 
        y[i]= M0*M1/(M1-M0)*(1/m1-1/m2);    
        for(l=0;l<Npole;l++) 
        {
           double m=poles[2*l],w=poles[2*l+1]; 
           double x1=(m1*m1/m -m)/w, x2=(m2*m2/m -m)/w;
           y[i]+=atan(x2)-atan(x1);              
        }       
     }
     impGrid(N, x, y);           
   }
   free(y);
}

static void getPoles(numout*cc, int nsub, char * s0,double mMin, double mMax, int*npole,double**poles)
{ int i,n,m,w;
  char*s;
  double mass;
  CalcHEP_interface * CI=cc->interface;
  
//  printf("Get Pole (%d %d %d)\n",s0[0],s0[1],s0[2]);
  
  *npole=0;
  if(s0[0]==0 || s0[1]==0) return;
      
  for(i=0; s0[i+1]!=0;)
  {  if(s0[i]>s0[i+1]) 
     { int a=s0[i]; s0[i]=s0[i+1];s0[i+1]=a; 
       if(i) i--;else i++;
     } else i++;
  }
  
  for(n=1;(s=CI->den_info(nsub,n,&m,&w));n++) 
  if(strcmp(s0,s)==0 && fabs(CI->va[m])>mMin && fabs(CI->va[m]) < mMax)
  { 
//printf(" %s %E %E  x= %E \n",CI->varName[m],CI->va[m], CI->va[w],(CI->va[m]-mMin)/(mMax-mMin) );

     (*npole)++;

//printf("npole=%d size=%d\n", *npole, sizeof(double)*2*(*npole));     
     *poles=realloc(*poles,sizeof(double)*2*(*npole));
     (*poles)[2*((*npole)-1)]= CI->va[m]; 
     (*poles)[2*(*npole-1)+1]= CI->va[w];  
   }
//    else  printf("-- %s %E %E  x= %E \n",CI->varName[m],CI->va[m], CI->va[w],(CI->va[m]-mMin)/(mMax-mMin) );

}


static double vegas_cycle(vegasGrid *vegPtr, double eps, double aeps,int maxStep, int NN, double fact,double *dI)
{ int k,l;
  double *ti=malloc(maxStep*sizeof(double));
  double *dti=malloc(maxStep*sizeof(double));
    double ii,dii,chi2;
  
  for(k=0;k<maxStep;k++)
  { 
    double s0=0,s1=0,s2=0;    
    vegas_int(vegPtr, NN , 1.5, nPROCSS  , ti+k, dti+k);
//    printf("ti=%E dti=%E  NN=%d \n",ti[k], dti[k],NN);
    if(dti[k]==0){ dii=0; ii=ti[k];  break;}
    for(l=k;l>=k/2;l--)
    { s0+=1/(dti[l]*dti[l]);
      s1+=ti[l]/(dti[l]*dti[l]);
      s2+=ti[l]*ti[l]/(dti[l]*dti[l]);
      if(l!=k)
      { 
        ii=s1/s0;
        dii=1/sqrt(s0);
        chi2=(s2-s1*s1/s0)/(k-l+1);
        if(chi2> 1 )dii*=sqrt(chi2);

        if(dii<eps*fabs(ii)) break;
        if(dii<aeps) break;
      }  
    }
    if(k && (dii<eps*fabs(ii) || dii<aeps )) break;
    NN*=fact;    
  }  
  free(ti); free(dti);
  *dI=dii;
  return ii;
}

double hCollider(double Pcm, int pp, int nf, double Qren,double Qfact, char * name1,char *name2,double pTmin,int wrt)
{ 
  double  sigma_tot=0, Qstat;
  int i;
  numout *cc;
  int n1,n2;
  int nout;
  double dI,m1,m2;
  
  if(name1==NULL  && name2==NULL) return 0;
  
  if(nf>5)nf=5; 
  if(nf<0)nf=0;
  
  
 if(name1) 
 { n1=pTabPos(name1);  
    if(n1==0) { printf("%s - no such particle\n",name1); return 0;}
    m1=pMass(name1);
 }
 else { n1=0; m1=0;}
 if(name2) 
 {  n2=pTabPos(name2);  
    if(n2==0) { printf("%s - no such particle\n",name2); return 0;}
    m2=+pMass(name2);   
 }
 else { n2=0; m2=0;}
 
  
  sMax=4*Pcm*Pcm; 
  if(pTmin<=0) sMin=m1+m2;
  else 
  {
    double MM=m1+m2;
    sMin=sqrt(MM*MM+pTmin*pTmin)+pTmin;
    pTmin_=pTmin;
  }
  sMin*=sMin; sMin+=1; 
  
  ppFlag=pp;   
  cc=colliderProduction( name1,name2, nf, pTmin>0);
  
  if(!cc) return 0;
   
  Q_fact=Qfact;
  Q_ren=Qren;
   
  if(Qaddress)
  { Qstat=*Qaddress;
    *Qaddress=sqrt(sMin);
//    printf("Q=%E\n",findValW("Q"));
    calcMainFunc();
  }  
  
  if(passParameters(cc)) return 0;

  *(cc->interface->gtwidth)=0;
  *(cc->interface->twidth)=0;
  *(cc->interface->gswidth)=0;
  sqme22=cc->interface->sqme;
      
  sigma_tot=0;
  if(pTmin>0) nout=1;else nout=0;
  if(n1)nout++;
  if(n2)nout++;
  for(nsub22=1;nsub22<=cc->interface->nprc; nsub22++) 
  { int pc[5];
    char*n[5];
    double tmp,dI;
    for(i=0;i<2+nout;i++) n[i]=cc->interface->pinf(nsub22,i+1,pmass+i,pc+i); 
    
    if(pc[0]<=pc[1])
    { pc1_=pc[0];
      pc2_=pc[1];
      
     if(wrt)for(i=0;i<2+nout;i++) {printf("%s ",n[i]); if(i==1) printf(" -> ");}
     
     switch(nout)
     { case 1: 
       { 
          double pcmIn=decayPcm(pmass[2],pmass[0], pmass[1]);
          double q;
          int err;
          
          if(pcmIn==0) return 0;
           pvect[0]=sqrt(pmass[0]*pmass[0]+pcmIn*pcmIn);
           pvect[1]=pcmIn; pvect[2]=0; pvect[3]=0;
           pvect[4]=sqrt(pmass[1]*pmass[1]+pcmIn*pcmIn);
           pvect[5]=-pcmIn; pvect[6]=0; pvect[7]=0;
           pvect[8]=pmass[2];pvect[9]=pvect[10]=pvect[11]=0;
           q=Q_fact>0? Q_fact: pmass[2];
          tmp=convStrFun2(pmass[2]*pmass[2]/sMax,q,pc1_,pc2_,ppFlag);
           q=Q_ren>0? Q_ren: pmass[2];
          tmp*=cc->interface->sqme(nsub22,sqrt(4*M_PI*parton_alpha(q)),pvect,NULL,&err);
          tmp*=389379660.0*M_PI/(2*pcmIn*pmass[2]*sMax);
          if(wrt)printf("cs=%E \n",tmp); 
          break;
       }
       case 2:  tmp=simpson(s_integrand,0.,1.,1.E-2); if(wrt)printf("cs=%E \n", tmp);  break;
       case 3:
       { double m3q;
         vegasGrid *vegPtr=vegas_init(5,veg_intergrand,50);
         char s0[3]={0,0,0};
         double eps,aEps;
         
         for(i3=2;i3<5;i++)  if( (pc[i]<6 && pc[i]>-6) || pc[i]==21) break;
         for(i4=2;i4<5;i4++) if(i4!=i3) break;
         for(i5=2;i5<5;i5++) if(i5!=i3 && i5!=i4) break; 

         M45_min=pmass[i4]+pmass[i5];
         S35_min=pmass[i3]+pmass[i5]; S35_min*=S35_min;
         S34_min=pmass[i3]+pmass[i4]; S34_min*=S34_min;
         
         m3q=pmass[i3]*pmass[i3];
         M45_max= sMax-2*sqrt((pTmin_*pTmin_+m3q)*sMax)+m3q;
         if(M45_max<=M45_min*M45_min) return 0;
         M45_max=sqrt(M45_max);
         S34_max=sqrt(sMax-pmass[i5]); S34_max*=S34_max;
         S35_max=sqrt(sMax-pmass[i4]); S35_max*=S35_max;
         
         s0[0]=i3+1;s0[1]=i4+1; getPoles(cc,nsub22,s0,sqrt(S34_min),sqrt(S34_max),&npole34,&pole34);
         s0[0]=i3+1;s0[1]=i5+1; getPoles(cc,nsub22,s0,sqrt(S35_min),sqrt(S35_max),&npole35,&pole35);        
         s0[0]=i4+1;s0[1]=i5+1; getPoles(cc,nsub22,s0,M45_min,      M45_max,      &npole45,&pole45); 
         s0[0]=1;   s0[1]=2;    getPoles(cc,nsub22,s0,sqrt(sMin),   sqrt(sMax),   &npole12,&pole12);                 

         setGrid(50,  vegPtr->x_grid[0]   , 0, 1,sqrt(sMin),sqrt(sMax),npole12,pole12);
         setGrid(50,  vegPtr->x_grid[1]   , 0, 1,M45_min,M45_max,npole45,pole45);
         
         if(npole34 && npole35)
         { setGrid(25,  vegPtr->x_grid[3], 0, 0.5, sqrt(S34_min),sqrt(S34_max),npole34,pole34);
           setGrid(25,  vegPtr->x_grid[3]+25, 0.5, 1, sqrt(S35_min),sqrt(S35_max),npole35,pole35);
         }else if(npole34)
           setGrid(50,  vegPtr->x_grid[3], 0, 1, sqrt(S34_min),sqrt(S34_max),npole34,pole34);
          else if (npole35)
           setGrid(50,  vegPtr->x_grid[3], 0, 1, sqrt(S35_min),sqrt(S35_max),npole35,pole35);
         eps=0.01;
         aEps=1E-6;
         if(fabs(sigma_tot)*eps>aEps) aEps=sigma_tot*eps;
         tmp=vegas_cycle(vegPtr,0.05, aEps, 20,5000, 1.1,&dI); 
         vegas_finish(vegPtr);
         if(wrt)printf("cs=%E +/-%E \n", tmp,dI);
         free(pole12); free(pole34); free(pole35); free(pole45);
         pole12=pole34=pole35=pole45=NULL;
         break;
       }  
     }
     sigma_tot+=tmp;
    }
  }  

  if(Qaddress){ *Qaddress=Qstat; calcMainFunc();} 
            
  return sigma_tot;
}

#ifdef plazmaWidth
static numout* plazmaWidth_cc;
static double plazmaWidth_T;
static double plazmaWidth_m[4];
static double plazmaWidth_integrand(double Pcm)
{ int err;
  double E1,E2,sqrt_s; 
  if(Pcm==0) return 0;  
  E1=sqrt(Pcm*Pcm+plazmaWidth_m[0]*plazmaWidth_m[0]);
  E2=sqrt(Pcm*Pcm+plazmaWidth_m[1]*plazmaWidth_m[1]);
  sqrt_s=E1+E2;
  if(sqrt_s<=plazmaWidth_m[2]+plazmaWidth_m[3]) return 0;
  
  return 4*bessk1(sqrt_s/plazmaWidth_T)*cs22(plazmaWidth_cc,1,Pcm, -1., 1. , &err)*pow(sqrt_s*Pcm,3.)/E1/E2; 
}

double plazmaWidth(char *process,double T)
{  char libName[40];
   plazmaWidth_T=T;
   process2Lib(process,libName);
   process2Mass(process,plazmaWidth_m);
   plazmaWidth_cc=getMEcode(0,ForceUG,process,NULL,NULL,libName); 
   return simpson(plazmaWidth_integrand,0., 5*T,1.E-3)/(4*M_PI*M_PI*plazmaWidth_m[1]*plazmaWidth_m[1]*bessk2(plazmaWidth_m[1]/T))
   /3.8937966E8;    
}
#endif
/*============ Fortran ==========*/

double cs22_(int*ccf,int*nsub,double*P,double*cos1,double*cos2,int*err)
{ numout*cc;
  memcpy(&cc,ccf,sizeof(cc));
  return cs22(cc,*nsub,*P,*cos1,*cos2 ,err);
} 

void  sethelicities_(double *h1,double *h2) { Helicity[0]=*h1; Helicity[1]=*h2;}


double hcollider_(double*Pcm, int*pp, int* nf, double*Qren,double*Qfact, char * name1,char *name2,double*pTmin,int*wrt,int len1,int len2)
{ 
  char cname1[20], cname2[20];
  char *cname1_,*cname2_;
  fName2c(name1,cname1,len1);
  fName2c(name2,cname2,len2);
  if(strlen(cname1)==0) cname1_=NULL; else  cname1_=cname1;
  if(strlen(cname2)==0) cname2_=NULL; else  cname2_=cname2;
  return  hCollider(*Pcm, *pp, *nf, *Qren,*Qfact, cname1_,cname2_,*pTmin,*wrt);  
}