/*
 Copyright (C) 2002, Alexander Pukhov pukhov@theory.sinp.msu.ru
*/
#include "syst.h"
#include "interface.h"
#include "alphas2.h"
#include "crt_util.h"
#include "plot.h"
#include "kinaux.h"
#include "const.h"
#include "strfun.h"
#include "sf_pdt.h"
#include "read_func.h"
#include "phys_val.h"
#include "parser.h"
#include<math.h>

double (*sf_alpha)(double)=NULL;

static double alphaMZ=0.1172;
static int alphaOrder=2;
static int alphaNF=5;
static double MbMb=4.2;
static double Mtp=175;
static int alphaPDF=1;
static char scale_str[61]="91.187";

static double alpha(int nf, int odr, double lambda,  double dscale)
{
    double b0 = 11. -  (2./3.)*nf;
    double b1 = 51. - (19./3.)*nf;
    double b2 = 2857. - (5033./9.)*nf + (325./27.)*nf*nf;
    double rl = 2*log(dscale / lambda);
    double alpha0= 4*M_PI/(b0*rl);
    double d__4 = log(rl) - .5;
    double d__2 = 2*b1/(b0*b0*rl);


    if(odr==1) return alpha0;
    else if(odr==2) return  alpha0*(1 - 2*b1*log(rl)/(b0*b0*rl));
    else if(odr==3) return  alpha0*(1 - 2*b1*log(rl)/(b0*b0*rl)
         + d__2*d__2 *(d__4*d__4 + b2*b0 /(8*b1*b1) - 1.25)  );
    else { fprintf(stderr,"Can not evaluate alpha in so large oder (%d).\n",odr);
           exit(1);
         }
}

static double findLambda(int nf,int odr, double alpha0, double M)
{ double l1=0.1, l2=0.3;
  double l,a,a1,a2;

  while((a1=alpha(nf,odr,l1,M)-alpha0) > 0)  l1*=0.7;

  while((a2=alpha(nf,odr,l2,M)-alpha0) < 0) l2*=1.3;

  do{ l=(l1*a2-l2*a1)/(a2-a1);
      a=alpha(nf,odr, l,M)-alpha0;
      if(fabs(a1)>fabs(a2)){ a1=a;l1=l;} else {a2=a;l2=l;}
    } while (fabs(a) > 0.00001*alpha0);
  return l;
}

#define  MCHARM   1.3

static double L3, L4, L5, L6;
static int   nf3,nf4,nf5,nf6;

static void init_alpha(void)
{
  double lambda,al; 
  int nf=alphaNF;

  if(nf==6) nf--;
  
  lambda=findLambda(nf,alphaOrder,alphaMZ,91.187);

  if(alphaNF>6) alphaNF=6; if(alphaNF<3) alphaNF=3;
  if(alphaOrder>3)alphaOrder=3; if(alphaOrder<1)alphaOrder=1;

  nf5=nf; L5=lambda;
  if(nf5==5) 
  {  if(MbMb<1.5*L5) {L4=L5;nf4=nf5;} else 
     { al=alpha(5,alphaOrder,L5,MbMb  );
       if(alphaOrder==3) al*=(1.+(11./72.)*pow(al/M_PI,2.));
       L4=findLambda(4,alphaOrder,al,MbMb);nf4=4; ;
     }
  } else {nf4=nf;L4=lambda;}
  if(nf4==4) 
  { if(MCHARM<1.5*L4) {L3=L4;nf3=nf4;} else
    { al=alpha(4,alphaOrder,L4,MCHARM);
      if(alphaOrder==3) al*=(1.+(11./72.)*pow(al/M_PI,2.));
      L3=findLambda(3,alphaOrder,al,MCHARM);nf3=3;
    }
  } else {nf3=nf; L3=lambda;}

  if(alphaNF==6)
  {  nf6=6;
     al=alpha(5,alphaOrder,L5,Mtp);
     if(alphaOrder==3) al*=(1.-(11./72.)*pow(al/M_PI,2.));
     L6=findLambda(6,alphaOrder,al,Mtp);
  }
  else {nf6=nf5; L6=L5;}
}

double alpha_2(double Q)
{
  if(alphaPDF && sf_alpha)  return (*sf_alpha)(Q); 

  if(Q<L3*1.5) Q=L3*1.5;
  
        if(Q<MCHARM)  return alpha(nf3, alphaOrder, L3, Q);
  else  if(Q<MbMb)    return alpha(nf4, alphaOrder, L4, Q);
  else  if(Q<Mtp)     return alpha(nf5, alphaOrder, L5, Q);
  else                return alpha(nf6, alphaOrder, L6, Q);
}

static int qq0(char *s, double *v)
{
  char key;
  char plist[20];
  int i;
  
  if(isdigit(*s)) { return 1==sscanf(s,"%lf%s",v,plist);}
  *v=100;  
  
  if(checkPhysVal(s,&key,plist) && strchr("STZM",key)) return 1;
   
  for(i=1;i<=nvar_int+nfunc_int;i++)
   if(strcmp(varName_int[i],s)==0 && strcmp(s,"GG")!=0) return 1;
   
  return 0;
}

int qcdmen_(void)
{
   void * pscr=NULL;
   int mode;
   int returnCode=0;

   initStrFun(0);
L10:{  char strmen[]="\030"
      " parton dist. alpha OFF "
      " alpha(MZ)=  ZZZZ       "
      " nf =        NF         "
      " order=      NNLO       "
      " mb(mb)=     MbMb       "
      " Mtop(pole)= Mtp        "
      " Q[Gev] = YYY           "
      " Alpha(Q) plot          ";

      if(alphaPDF)
      { if(sf_alpha)improveStr(strmen,"OFF","%-.3s","ON");
        else         improveStr(strmen,"OFF","%-.3s","!ON");
      }
      improveStr(strmen,"ZZZZ","%.4f", alphaMZ);
      improveStr(strmen,"NF","%d",alphaNF);

           if(alphaOrder==1) improveStr(strmen,"NNLO","%-.4s","LO");
      else if(alphaOrder==2) improveStr(strmen,"NNLO","%-.4s","NLO");
      else alphaOrder=3;
      
      improveStr(strmen,"MbMb","%.3f", MbMb);
      improveStr(strmen,"Mtp","%.2f", Mtp);
      
      improveStr(strmen,"YYY","%-.14s", scale_str);

      menu1(54,8,"QCD alpha",strmen,"n_alpha",&pscr,&mode);
    }
    switch (mode)
    { case 0: if(returnCode) init_alpha(); return returnCode;
      case 1: alphaPDF=!alphaPDF;
              if(alphaPDF && !sf_alpha)
               messanykey(20,20,"WARNING! Parton distribution does not\n"
                                "define alphaQCD ");
              returnCode=1;
              break;
      case 2:
         { double alphaMZ_old=alphaMZ; 
           if(correctDouble(3,15,"Enter new value ",&alphaMZ,1)) returnCode=1;
           if(alphaMZ>0 && alphaMZ<0.3) returnCode=1; else  
           { alphaMZ=alphaMZ_old;
              messanykey(5,15,"Your input is out of alphaMZ range");
           }
         }
         break;
      case 3:
         { int NF_old=alphaNF;
           if(correctInt(3,15,"Enter new value ",&alphaNF,1))
           {
              if(alphaNF<=6 && alphaNF>=3) returnCode=1;
              else { messanykey(5,15,"NF out of range"); alphaNF=NF_old;}
           }   
         }
         break;
      case 4:
         {  char lomen[]="\010"
               " LO     "
               " NLO    "
               " NNLO   ";
            void *pscrlo=NULL;
            int k=0;
            menu1(52,12,"",lomen,"",&pscrlo,&k);
            if(k) { alphaOrder=k; returnCode=1; put_text(&pscrlo);}
         }
         break;
      case 5: correctDouble(3,15,"Enter new value ",&MbMb,1);  break;
      case 6: correctDouble(3,15,"Enter new value ",&Mtp,1);   break;   
      case 7:
         { int npos=1,rc;
           do
           {  double  dscale;
              goto_xy(2,12); print("QCD scale: ");
              if(str_redact(scale_str,npos,60)==KB_ENTER) returnCode=1;
              goto_xy(2,12); clr_eol();
              rc=calcExpression(scale_str, qq0, &dscale);
              if(rc && rc != cannotevaluate)
              { npos=rderrpos;
                if(rc==unknownidentifier) messanykey(10,10," Unknown parameter");
                else messanykey(10,10," Syntax error");
              }
           }  while(rc && rc != cannotevaluate);
         }
	 break;
      case 8:
	{ void * screen;
	  int i;
	  double f[150];

	  static double qMin=1, qMax=1000;
	  static int nPoints=100;
	  if(returnCode) init_alpha();
	  get_text(1,1,maxCol(),maxRow(),&screen);

          if(correctDouble(40 ,15 ,"Q_min=",&qMin,0)&& qMin>=0.5
          && correctDouble(40 ,16 ,"Q_max=",&qMax,0)&& qMax>qMin
          && correctInt(33,17,"number of points=" ,&nPoints,0)
          && nPoints>3&& nPoints<=150)
	  {
	    for(i=0;i<nPoints;i++)
	    { double z=(double)i/(double)(nPoints-1);
	      double q=pow(qMin,1-z)*pow(qMax,z);
	      f[i]=alpha_2(q);
	    }
	    plot_1(log10(qMin),log10(qMax),nPoints,f,NULL, " ", "log10(Q/GeV)", "Alpha(Q)");
	  } else  messanykey(40,18,
	          " Correct input is \n"
	          " 0.5<= Q_min <Q_max\n"
	          " number of points <=150 and >=4");
	  put_text(&screen);
	}	
    }
    goto L10;
}

int w_qcd__(FILE * mode)
{
  fprintf(mode,"alphaPDF=%d alpha(MZ)=%E NF=%d Order=%d MbMb=%E Mtp=%E Scale= %s",
              alphaPDF, alphaMZ, alphaNF, alphaOrder, MbMb, Mtp, scale_str);
  return 0;
}

void i_qcd(void)
{
   if(nin_int==1) strcpy(scale_str,"M1"); else strcpy(scale_str,"M12");
   init_alpha();
}

int r_qcd__(FILE *mode)
{ int n;
  n=fscanf(mode, "alphaPDF=%d alpha(MZ)=%lf NF=%d Order=%d MbMb=%lf Mtp=%lf Scale= %[^\n]",
            &alphaPDF, &alphaMZ,&alphaNF,&alphaOrder,&MbMb,&Mtp, scale_str);
  trim(scale_str);
  init_alpha();
  if(n!=7) return 1; else return 0;
}

void alf_(double q)
{
  static int k=-1;
  if(!k) return;
  if(k==-1)
  {  for(k = 1; k <= nvar_int; ++k) 
      if (strcmp("GG", varName_int[k]) == 0) break; 
     if(k>nvar_int) k=0;
  }

  if(k>0) va_int[k]=sqrt(4*M_PI*alpha_2(q));    
}


static int qq1(char *s,double *v)
{
  char key, plist[20];
  int i;

  if(isdigit(*s)) { sscanf(s,"%lf",v); return 1;}
  if(checkPhysVal(s,&key,plist)) {*v=calcPhysVal(key,plist); return 1;}
 
  for(i=1;i<=nvar_int+nfunc_int;i++)
    if(strcmp(varName_int[i],s)==0) { *v=va_int[i]; return 1;}
 
  return 0;
}

double Scale(void)
{
   double dscale;
   if(calcExpression(scale_str ,qq1, &dscale))
   {  fprintf(stderr, " ERROR in evaluation of  QCD scale\n");
      sortie(50);
   }
   if (dscale < 1.5*L3 )  dscale = 1.5*L3;
   return dscale;
}
