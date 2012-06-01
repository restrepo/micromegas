#include"../sources/micromegas.h"
#include"../sources/micromegas_aux.h"

static double amotry(double *p, double *y, int ndim,
	     double (*f)(double *), int ilo, double fac)
{
   int i,j;
   double  ytry, fac1=(1.0-fac)/ndim;
   double * p_buff=p+(ndim+1)*ndim;
   double * p_ilo =p+ilo*ndim;
   
   for(j=0;j<ndim;j++) p_buff[j]=p_ilo[j]*fac;
   for(i=0;i<=ndim;i++)  if(i!=ilo) 
     {double *p_i=p+i*ndim;  for(j=0;j<ndim;j++) p_buff[j] +=p_i[j]*fac1;} 
   ytry=(*f)(p_buff);
   
   if (ytry > y[ilo]) 
   {  for(j=0;j<ndim;j++) p_ilo[j]=p_buff[j];
      y[ilo]=ytry;
   }
//printf("amotry returns fac=%f %E\n",fac,  ytry);   
   return ytry;
}

static double amoeba(double *p, double * y, int ndim, double (*f)(double *), 
                                                    double eps, int *nCalls)
{
   int i,ilo,ihi,inlo,j;
   double ysave,ytry;

   for (;;) 
   {
      ihi=0;									     
      ilo = y[0]<y[1] ? (inlo=1,0) : (inlo=0,1);				     
      for (i=0;i<=ndim;i++)							     
      {										     
     	 if (y[i] >= y[ihi]) ihi=i;
     	 if (y[i] < y[ilo]) { inlo=ilo; ilo=i; } 
     	 else if (y[i] < y[inlo] && i != ilo) inlo=i;
      }										     

//printf("nCall=%d  ndim=%E\n",*nCalls,y[ilo]);   
     										     
      if((*nCalls)<=0||2*(y[ihi]-y[ilo])/(fabs(y[ilo])+fabs(y[ihi]))<eps)break;
     										     
      ytry=amotry(p,y,ndim,f,ilo,-1.0); (*nCalls)--;				     
      if (ytry >= y[ihi]) {ytry=amotry(p,y,ndim,f,ilo,2.); (*nCalls)--;}	     
      else if (ytry <= y[inlo])							     
      {										     
         ysave=y[ilo];								     
     	 ytry=amotry(p,y,ndim,f,ilo,0.5);  (*nCalls)--;
     	 if (ytry <= ysave)
     	 {  
     	    for (i=0;i<=ndim;i++)
     	    {  double * p_ihi=p+ihi*ndim;
               if (i != ihi)
     	       {  double * p_i=p+i*ndim;
     		  for(j=0;j<ndim;j++) p_i[j]=0.5*(p_i[j]+p_ihi[j]);
     		  y[i]=(*f)(p_i);
     	       }
            }
/*printf("srink\n");            */
     	    (*nCalls) -= ndim;
         }									     
      }										     
   }
   return y[ihi];
}
/*========================== end of amoeba ================*/



void DAMA_NAI( double vRotation, double vSun, double vEsc,
double Mwimp, double csSIp,double csSIn,double csSDp,double csSDn )
{ double dNdE[200];
  double pNa, pI ,Exposure,NaQF,IQF;
  double A4,A5,A6,A14, Y4,Y5,Y6,Y14, 
  Na4_s,Na5_s,Na6_s,Na14_s, Na4_w,Na5_w,Na6_w,Na14_w,
  I4_s,I5_s,I6_s,I14_s,     I4_w,I5_w,I6_w,I14_w;
/*
     arXiv:astro-ph/0307403v1
     arXiv:astro-ph/0501412
     arXiv:astro-ph/0804.2741
*/ 
   pNa=(23./150.);
   pI =(127./150.);  
   NaQF=0.3;    /* Quenching factor:  E_e/E_Na */     
   IQF=0.09;    /* Quenching factor:  E_e/E_I  */
   SetfMaxwell(vRotation,/*vSun+15,*/vEsc);  /* Summer  signal */
   Vearth=vSun+15;
   
    nucleusRecoilAux(Maxwell,23,Z_Na,J_Na23,S00Na23,S01Na23,S11Na23,
     csSIp,csSIn,csSDp, csSDn, dNdE); 
     Na4_s =cutRecoilResult(dNdE,2/NaQF, 4/NaQF);
     Na5_s =cutRecoilResult(dNdE,2/NaQF, 5/NaQF);
     Na6_s =cutRecoilResult(dNdE,2/NaQF, 6/NaQF);
     Na14_s=cutRecoilResult(dNdE,6/NaQF,14/NaQF);

    nucleusRecoilAux(Maxwell,127,Z_I,J_I127,S00I127,S01I127,S11I127,
    csSIp,csSIn,csSDp, csSDn, dNdE); 
     I4_s =cutRecoilResult(dNdE,2/IQF, 4/IQF);
     I5_s =cutRecoilResult(dNdE,2/IQF, 5/IQF);
     I6_s =cutRecoilResult(dNdE,2/IQF, 6/IQF);
     I14_s=cutRecoilResult(dNdE,6/IQF,14/IQF);
    
   SetfMaxwell(vRotation/*,vSun-15*/,vEsc);  /* winter signal */
   Vearth=vSun-15;
    nucleusRecoil0Aux(Maxwell,23,Z_Na,J_Na23,Sp_Na23,Sn_Na23,  
     csSIp,csSIn,csSDp, csSDn, dNdE); 
     Na4_w =cutRecoilResult(dNdE,2/NaQF, 4/NaQF);
     Na5_w =cutRecoilResult(dNdE,2/NaQF, 5/NaQF);
     Na6_w =cutRecoilResult(dNdE,2/NaQF, 6/NaQF);
     Na14_w=cutRecoilResult(dNdE,6/NaQF,14/NaQF);

    nucleusRecoilAux(Maxwell,127,Z_I,J_I127,S00I127,S01I127,S11I127,
    csSIp,csSIn,csSDp, csSDn, dNdE); 
     I4_w =cutRecoilResult(dNdE,2/IQF, 4/IQF);
     I5_w =cutRecoilResult(dNdE,2/IQF, 5/IQF);
     I6_w =cutRecoilResult(dNdE,2/IQF, 6/IQF);
     I14_w=cutRecoilResult(dNdE,6/IQF,14/IQF);
/* Modulation */
   A4= 0.5*( (Na4_s -Na4_w )*pNa+ (I4_s -I4_w)*pI )/2;
   A5= 0.5*( (Na5_s -Na5_w )*pNa+ (I5_s -I5_w)*pI )/3;
   A6= 0.5*( (Na6_s -Na6_w )*pNa+ (I6_s -I6_w)*pI )/4;
   A14=0.5*( (Na14_s-Na14_w)*pNa+ (I14_s-I14_w)*pI)/8;
/* Yield */
   Y4= 0.5*( (Na4_s +Na4_w )*pNa+ (I4_s +I4_w)*pI )/2;
   Y5= 0.5*( (Na5_s +Na5_w )*pNa+ (I5_s +I5_w)*pI )/3;
   Y6= 0.5*( (Na6_s +Na6_w )*pNa+ (I6_s +I6_w)*pI )/4;
   Y14=0.5*( (Na14_s+Na14_w)*pNa+ (I14_s+I14_w)*pI)/8;

   
   printf("DAMA_NaI signal\n");
   printf(" A4 =%.2E  A5 =%.2E A6 =%.2E A14 =%.2E\n",A4,A5,A6,A14);
   printf(" Y4 =%.2E  Y5 =%.2E Y6 =%.2E Y14 =%.2E\n",Y4,Y5,Y6,Y14);
   printf("xi^2(LIBRA)=%.2E\n",  
    pow((A4-0.0213)/0.0032,2)+pow((A5-0.0165)/0.0024,2)+pow((A6-0.0107)/0.0019,2));
   printf("xi^2(LIBRA+NaI)=%.2E\n",  
    pow((A4-0.0223)/0.0027,2)+pow((A5-0.0178)/0.0020,2)+pow((A6-0.0131)/0.0016,2));

}


void Xenon10(double Mwimp, double csSIp,double csSIn,double csSDp,double csSDn )
{ 
/*
   astro-ph/0706.0039v2   
*/   
  double dNdE[200]; 
  double N129,N131,N129_,N131_,NXX,NXX_,N,N_;
  double Emin=4.5;
  
   nucleusRecoilAux(Maxwell,129,Z_Xe,J_Xe129,S00Xe129,S01Xe129,S11Xe129,
   csSIp,csSIn,csSDp, csSDn, dNdE);  

    N129=cutRecoilResult(dNdE,Emin,26.9);
    N129_=cutRecoilResult(dNdE,Emin,15.*1.9/1.23);

   nucleusRecoilAux(Maxwell,131,Z_Xe,J_Xe131,S00Xe131,S01Xe131,S11Xe131,
   csSIp,csSIn,csSDp, csSDn, dNdE);  

    N131=cutRecoilResult(dNdE,Emin,26.9);
    N131_=cutRecoilResult(dNdE,Emin,15.);

    nucleusRecoil0Aux(Maxwell,132,Z_Xe,0,0.,0.,
    csSIp,csSIn,csSDp, csSDn, dNdE);  

    NXX=cutRecoilResult(dNdE,Emin,26.9);
    NXX_=cutRecoilResult(dNdE,Emin,15.);


    N= 136*(N129 *0.264 + N131 *0.212 + NXX *(1-0.264-0.212) );
    N_=136*(N129_*0.264 + N131_*0.212 + NXX_*(1-0.264-0.212) );

    printf("Xenon10 signal:\n");

    printf(" Events expected: [4.5-26.9]KeV=%.2E  [4.5-15.]KeV=%.2E \n",N,N_ );
    printf(" Excluded with probability %.2E\n", MaxGapLim(N_,N));  
}

void CDMS_123(double Mwimp, double csSIp,double csSIn,double csSDp,double csSDn )
{ 
/* 
   astro-ph/0808.3530v2
*/
  double dNdE[200]; 
  double NevSD,NevSI, Nev;
    
  nucleusRecoilAux(Maxwell,73,Z_Ge,J_Ge73,S00Ge73A,S01Ge73A,S11Ge73A,
    csSIp,csSIn,csSDp, csSDn, dNdE);
  NevSD=cutRecoilResult(dNdE,10.,100.)*0.0773;
 
  nucleusRecoil0Aux(Maxwell,73,Z_Ge,0,0,0,
   csSIp,csSIn,csSDp, csSDn, dNdE);
  NevSI=cutRecoilResult(dNdE,10.,100.)*(1-0.0773);   
  Nev=121.3*(NevSD+NevSI);
  printf("CDMS_123 events %.2E\n", Nev); 
}     

void CDMS_2006(double Mwimp, double csSIp,double csSIn,double csSDp,double csSDn )
{ 
/* astro-ph/0509269v2 */
  double dNdE[200]; 
  double NevGe,NevSi;

    nucleusRecoilAux(Maxwell,73,Z_Ge,J_Ge73,S00Ge73,S01Ge73,S11Ge73,
       csSIp,csSIn,csSDp, csSDn, dNdE);        
    NevGe=11.5*cutRecoilResult(dNdE,10.,20);

    nucleusRecoilAux(Maxwell,73,Z_Ge,0,S00Ge73,S01Ge73,S11Ge73,
       csSIp,csSIn,csSDp, csSDn, dNdE);        
    NevGe+=11.5*(1/0.0773-1)*cutRecoilResult(dNdE,10.,20);
    
    
    nucleusRecoilAux(Maxwell,29,Z_Si,J_Si29,S00Si29,S01Si29,S11Si29,
     csSIp,csSIn,csSDp, csSDn, dNdE);  
    
    NevSi=1.7*cutRecoilResult(dNdE,10.,20);

    nucleusRecoilAux(Maxwell,29,Z_Si,0.,S00Si29,S01Si29,S11Si29,
     csSIp,csSIn,csSDp, csSDn, dNdE);  
    
    NevSi+=1.7*(1/0.0468-1)*cutRecoilResult(dNdE,10.,20);
    
    printf("CDMS_2006 events:= %.2E (%.2E Ge, %.2E Si)\n",NevGe+NevSi, NevGe,NevSi);
}     

void He3_100g(double Emin, int nMonths, double Mwimp, double csSIp,double csSIn,double csSDp,double csSDn )
{ 
  double dNdE[200]; 
  double Yield;

  nucleusRecoil0Aux(Maxwell,3,Z_He,J_He3,Sp_He3,Sn_He3,
   csSIp,csSIn,csSDp, csSDn, dNdE);        
  Yield=cutRecoilResult(dNdE,Emin,6.);
    
  printf("He3 signal  events:= %.2E (100g x %d months)\n",0.1*30*nMonths*Yield,
  nMonths);
}     

static CoGeNT_plot(double quenching)
{
/*
1004.0697
*/
  int data[28]={49,36,25,21,19,20,18,40,44,36,
                10, 5, 5,10, 4, 4,10, 3, 2, 5,
                3,1,6,3,6,3,3,1};
  double Emin=0.4,Emax=3.2,dE=0.1,exposure=0.33*56;

  
  int i;

  FILE*f=fopen("_plot_","w");
  
exposure=1;  
  
  fprintf(f,"#type 1 %%histogram\n");
  fprintf(f,"#title %s %.2f\n","CoGeNT results for quenching",quenching);
  fprintf(f,"#yName dN/dE\n");
  fprintf(f,"#xMin %f\n", Emin/quenching);
  fprintf(f,"#xMax %f\n", Emax/quenching);
  fprintf(f,"#xDim %d\n",28);
  fprintf(f,"#xName  Recoil Energy[KeV]\n");
  for(i=0;i<=28;i++) fprintf(f,"%E  %E\n",data[i]/(dE/quenching)/exposure, 
                                     sqrt(data[i])/(dE/quenching)/exposure );
  fclose(f);
  system("trap \"rm -f _plot_\" 0 1 2 3 9; ../CalcHEP_src/bin/plot_view < _plot_ &");
}



static double GeandZn(double E)
{ double  q1=(E-1.298)/0.0975  ,q2=(E-1.1)/0.0975;
  return exp(-q1*q1)+0.4*exp(-q2*q2);
}


static double C[4]={2.216547E+01, 4.136623E+02, 1.155849E+03, 2.279175E+00};

static double dNdE_stat[200];

static double FFS(double E)
{double   Q=0.2;
     return  dNdERecoil(dNdE_stat,E/Q)/Q;
}     

static double FFB(double E)
{ return  C[0]+ C[1]*GeandZn(E) + C[2]*exp(-C[3]*(E-0.45)); }     

double CoGeNTchi2(double *c)
{
  int data[28]={49,36,25,21,19,20,18,40,44,36,
                10, 5, 5,10, 4, 4,10, 3, 2, 5,
                3,1,6,3,6,3,3,1};
  double Emin=0.4,Emax=3.2,dE=0.1,exposure=0.33*56;
  double Q=0.199935;
  int i;
  double chi2=0;

  for(i=0;i<3;i++)C[i]=fabs(c[i]);
  C[3]=fabs( c[3]-0.1) +0.1;
   
  for(i=0;i<28;i++)
  {  double dchi=data[i];
     double E1=pow( (Emin+    i*dE)/Q,1./1.1204);
     double E2=pow( (Emin+(i+1)*dE)/Q,1./1.1204);

/*
     double E1=(Emin+    i*dE)/Q;
     double E2=(Emin+(i+1)*dE)/Q;
*/
   
       
     dchi -= cutRecoilResult(dNdE_stat,E1,E2); 
     dchi -= C[0]*dE;
     dchi -= C[1]*simpson(GeandZn,Emin+i*dE,Emin+(i+1)*dE,1.E-3);
     dchi -= C[2]/C[3]* (exp(-C[3]*(Emin+i*dE-0.45))- exp(-C[3]*(Emin+(i+1)*dE-0.45)));
     chi2+= dchi*dchi/data[i]; 
  }
  return -chi2;
}



void CoGeNT(double Mwimp, double csSIp,double csSIn,double csSDp,double csSDn )
{ 
/*
  1002.4703
  1003.0014
  1004.0697   
*/
  double dNdE[200]; 
  double Yield;
  double quenching=0.2;
  int i;
  
  setRecoilEnergyGrid(0.1,200);
  
  nucleusRecoilAux(Maxwell,73,Z_Ge,J_Ge73,S00Ge73A,S01Ge73A,S11Ge73A,
      csSIp,csSIn,csSDp, csSDn, dNdE);

  for(i=0;i<200;i++) dNdE_stat[i]= 0.33*56*dNdE[i];       
  
  CoGeNT_plot(1.);  
  displayFunc(FFS, 0.4,3.2, "signal");  
  { /*  Backgraund  fitting */
     double X[24], Y[5],chi2,yield1;
     int i,j, N; 
     yield1=0.33*56*cutRecoilResult(dNdE,0.4/quenching,0.5/quenching);
     if(yield1>100) { printf("Excluded by first bin\n"); return;}   
  
     C[0]=20;
     C[1]=400;
     C[3]=1;
     if(yield1>49) C[2]=0; else C[2]=(49-Yield)/0.1;
     for(i=0;i<5;i++) for(j=0;j<4;j++)  X[i*4+j]=C[j];
     for(i=0;i<4;i++) X[i*4+i]+=1;

     printf("Chi^2 for starting points: "); 
     for(i=0;i<5;i++) { Y[i]=CoGeNTchi2(X+i*4); printf(" %.1E ", -Y[i]);}
     printf("\n");
  
     N=1000; chi2=-amoeba(X, Y, 4, CoGeNTchi2,1.e-4, &N);
                                                      
     printf("Chi^2=%e\n",chi2);
     if(N<0) printf(" No convergence   for 1000 points!\n"); 
     for(i=0;i<3;i++) C[i]=fabs(X[i]);
     C[3]= 0.1+fabs(X[3]-0.1);
     {  char txt[50];
        sprintf(txt,"BG: %.0f + %.0f*exp(-%.1f*(E-0.45))+%.0f*Pole(E); ch^2=%.1f",C[0],C[2],C[3],C[1], chi2); 
        displayFunc(FFB, 0.4,3.2, txt);
     }     
  }  
  setRecoilEnergyGrid(-1,-1);
}     



void COUPP(double Mwimp, double csSIp,double csSIn, double csSDp,double csSDn)
{ /*
    astro-ph/0804.2886v1
  */  
  double dNdE[200]; 
  double pC,pF,pI,NC,NF,NI;
  double Emin=10, Exposure=52; /* GeV*Day */ 
 
  pC=12 /(12.+3*19.+127.);
  pF=19*3/(12.+3*19.+127.);
  pI=127/(12.+3*19.+127.);

  nucleusRecoil0Aux(Maxwell,12,6,0,0,0,
    csSIp,csSIn,csSDp, csSDn, dNdE);
  NC=cutRecoilResult(dNdE,Emin,30);

  nucleusRecoilAux(Maxwell,19,Z_F,J_F19,S00F19,S01F19,S11F19,
     csSIp,csSIn,csSDp, csSDn, dNdE);  
  NF=cutRecoilResult(dNdE,Emin,30);


  nucleusRecoil0Aux(Maxwell,19,Z_F,J_F19,Sp_F19,Sn_F19,
     csSIp,csSIn,csSDp, csSDn, dNdE);  
  NF=cutRecoilResult(dNdE,Emin,30);
 
  nucleusRecoilAux(Maxwell,127,Z_I,J_I127,S00I127,S01I127,S11I127,
     csSIp,csSIn,csSDp, csSDn, dNdE);  
  NI=cutRecoilResult(dNdE,Emin,30);
  
  printf("COUP events:=%.2E\n", (NC*pC+NF*pF+NI*pI)*Exposure);
   
}


int main(int narg, char** args)
{ double  csSIp,csSIn,csSDp,csSDn;
  double vRotation=230, vSun=225.2,vEsc=600.;

  if(narg < 6)
  {
    printf("Parameters expected: Mwinp[GeV], csSIp, csSIn, csSDp,csSDn[pb]\n"); 
    return 1;
  }   

  if(sscanf(args[1],"%lf",&Mcdm)==1 &&
     sscanf(args[2],"%lf",&csSIp)==1 &&
     sscanf(args[3],"%lf",&csSIn)==1 &&
     sscanf(args[4],"%lf",&csSDp)==1 &&
     sscanf(args[5],"%lf",&csSDn)==1) 
    printf("Model parameters:\n Mcdm=%.2E csSIp=%.2E, csSIn=%.2E,csSDp=%.2E,csSDn=%.2E\n",
                          Mcdm, csSIp, csSIn,csSDp,csSDn);
  else { printf("Wrong parameters\n"); return 2;}  

  DAMA_NAI(vRotation, vSun, vEsc,Mcdm, csSIp,csSIn,csSDp,csSDn);

  SetfMaxwell(vRotation/*,vSun*/,vEsc);
  Vearth=vSun; 
  Xenon10(Mcdm, csSIp,csSIn,csSDp,csSDn );
  CDMS_123(Mcdm, csSIp,csSIn,csSDp,csSDn ); 
  CDMS_2006(Mcdm, csSIp,csSIn,csSDp,csSDn);
  COUPP(Mcdm, csSIp,csSIn,csSDp,csSDn);   

  CoGeNT(Mcdm, csSIp,csSIn,csSDp,csSDn);   
  
/* He3_100g(1., 2, Mcdm, csSIp,csSIn,csSDp,csSDn); */

}
