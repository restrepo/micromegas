/*
EDM   in MSSM with CP violation
Notation of  hep-ph/
*/

#ifdef __cplusplus

#include <cmath>
#include <complex>
using std::complex;

typedef complex<double> Complex;
#define creal real 
#define cimag imag
#define carg  arg
#define clog  log
#define cexp  exp
#define cpow  pow
const Complex I=Complex(0,1);

#else 

#include <complex.h>
#define Complex double complex

#endif

#include"../../sources/micromegas_aux.h"
#include"../../sources/micromegas.h"

#include"pmodel.h"
/*
#include"cpsh.h"
*/

static double findValW2(char * name)
{ int err; 
  double m;
  err=findVal(name,&m);
  if(err==0) return m; else 
  { printf(" Parameter %s is not found\n", name);
    return 0;
  }
}

static double f1(double z) 
{ return 3.0*z/2./(z-1.)/(z-1.)/(z-1.)*(z*z-4.*z+3.+2*log(z)); }
   				  
static double f2(double z) 
{  return  3./(z-1.)/(z-1.)/(z-1.)*(z*z-1.-2.*z*log(z)); }

static double _zz;             /* to pass argument */ 
static double pEps=1.E-10;      /* precision of integral calculation */
                                  
static double _F(double u)
{ double x=u*u, xx=x*(1-x)/_zz;
  double eps;
  if(x==0.) return 0; 
  eps=1-xx;

  if(eps<0.01 && eps>-0.01) return -u*xx*(1+eps*(0.5+eps*(1./3.+eps*(0.25+eps*0.2))));

  return  u*xx/eps*log(xx);
}   



static double F(double z)                      
{ _zz=z;
  return 4*simpson(_F,0.,M_SQRT1_2 ,pEps); 
}
                                  
static double _f(double u)
{ double x=u*u, xx=x*(1-x)/_zz;
  double eps;
  if(x==0.) return 0; 
  eps=1-xx;
  
  if(eps<0.01 && eps>-0.01) return u*(1-2*_zz*xx)*(1+eps*(0.5+eps*(1./3.+eps*(0.25+eps*0.2))));

  return  -u*(1-2*_zz*xx)/eps*log(xx);
}   

static double f(double z)                      /* function 3.5 */
{ _zz=z;
  return 2*simpson(_f,0.,M_SQRT1_2 ,pEps); 
}

static double _g(double u)
{ double x=u*u, xx=x*(1-x)/_zz;
  double eps;
  if(x==0.) return 0; 
  eps=1-xx;
  
  if(eps<0.01 && eps>-0.01) return u*(1+eps*(0.5+eps*(1./3.+eps*(0.25+eps*0.2))));

  return  -u/eps*log(xx);
}   

static double g(double z)                      /* function 3.6 */
{ _zz=z;
  return 2*simpson(_g,0.,M_SQRT1_2 ,pEps); 
}

static double I2(double x1, double x2, double x3)
{  return (x1*x2*log(x1/x2)+x2*x3*log(x2/x3)+x1*x3*log(x3/x1))/
   (x1-x2)/(x2-x3)/(x1-x3);
}

static Complex Zafi(double A, double fi)
{  fi*=M_PI/180.;
   return A*(cos(fi)+I*sin(fi));
}
        				
static	double alphem =  1.0/137.036;
static	double alphew; 
static  double alpsmz  ; /* alpha_s(Mz) */
static	double  b, sb, cb;

/*static double alps(double mu, double mz, double alpsmz)
{ double v=1.-23.*alpsmz*log(mz/mu)/6./M_PI;
  return alpsmz*(1-116.*alpsmz*log(v)/v/4./M_PI/23.)/v;
}
*/


void edm_(double* el,double * thallium )
{
   char namer[20], namei[20];
   double nlr[5][3],y[5][3],mn[5],Zh[4][4],mh[4];
   double xit[4],xib[4],zetat[4],zetab[4],xt1h[4],xt2h[4],xb1h[4],xb2h[4];		
   double me=0.511e-3;
   double mw,mz,msn,mch1,mch2,sw,tw,ee,x1,x2,c1,c2,edm_ch,edm_nt,mse1,mse2,fneut,ye,edm1loop;
   double mselsq,msersq,ml1,mr1,ylep,xlep,delta,delta2,th2sel,c2b,tb;		
   double ht,hb,mst1,mst2,msb1,msb2,edm2loop,edmh,mb,mt;
   double phil1;
   double csg,Cs,xth[4],xbh[4],gStt[4],gPtt[4],gSbb[4],gPbb[4],edmhq;
   int i,j,a;
   Complex mc11,mc12,mc21,mc22;

   Complex ze;
   Complex nl[5][3],nr[5][3];
   Complex
   Zu11,Zu21,Zu12,Zu22,Zv11,Zv12,Zv21,Zv22,N[5][5],E[3][3],Ae,mu,At,Xt,Ab,Xb,hbc,htc;
   Complex deltahb,Deltahb,deltaht,Deltaht,mgl; 
   Complex md11,md12,md21,md22;
   double csgt,csgb,csgst,csgsb,edmtq,edmbq;
   double ac1[4],ac2[4],bc1[4],bc2[4],edmhch,xc1h[4],xc2h[4];
   double edmhch1,edmhch2,edmhch3;
   double alpss;
   double mbsusy,mtsusy,hbs,hts;
   b=atan(findValW("tb"));
   mw=findValW2("MW");
   mz=findValW2("MZ");
   ee=findValW2("EE");
   alpsmz=findValW2("alfSMZ");
   alphew = ee*ee/4./M_PI;
   sb=sin(b);
   cb=cos(b);
   tb=tan(b);
   c2b=cb*cb-sb*sb;
   
   ye=me/(M_SQRT2*mw*cb);
   msn=findValW2("MSnl");
   mch1=findValW2("MC1");
   mch2=findValW2("MC2");
   sw=findValW2("SW");
   ee=findValW2("EE");
   mt=findValW2("Mtp");
   mb=findValW2("MbMt");
   hb=ee*mb/(mw*sw*cb);
   ht=ee*mt/(mw*sw*sb);

    
   tw=sqrt(sw*sw/(1.-sw*sw));
   
   Zu11=findValW2("Zu11r")+I*findValW2("Zu11i");
   Zu12=findValW2("Zu12r")+I*findValW2("Zu12i");
   Zu21=findValW2("Zu21r")+I*findValW2("Zu21i");
   Zu22=findValW2("Zu22r")+I*findValW2("Zu22i");
   Zv11=findValW2("Zv11r")+I*findValW2("Zv11i");
   Zv12=findValW2("Zv12r")+I*findValW2("Zv12i");
   Zv21=findValW2("Zv21r")+I*findValW2("Zv21i");
   Zv22=findValW2("Zv22r")+I*findValW2("Zv22i");
   x1=mch1/msn; x1*=x1;
   x2=mch2/msn; x2*=x2;
   
   

   c1=cimag(conj(-Zv11)*Zu12);
   c2=cimag(conj(-Zv21)*Zu22);
   
   edm_ch=alphew/12./M_PI/sw/sw*ye*(f1(x1)/mch1*c1+f1(x2)/mch2*c2);


/* define the selectron mixing matrix*/
/* for the moment take A universal for the sleptons -need to introduce new
parameter*/
   ml1=findValW2("Ml1");
   mr1=findValW2("Mr1");
   mselsq=ml1*ml1+me*me+mz*mz*c2b*(-0.5+sw*sw);
   msersq=mr1*mr1-me*me+mz*mz*c2b*(sw*sw);
   Ae=Zafi(findValW2("aAe"),findValW2("fiAe")); 
   mu=Zafi(findValW2("aMu"),findValW2("fiMu"));
   ze=me*(Ae-tb*conj(mu));
   
   phil1=carg(ze);
   delta2=(mselsq-msersq)*(mselsq-msersq)+4*creal(ze*conj(ze));
   
   if(delta2>0) delta=sqrt(delta2);
   else printf("wrong parameter in selectron matrix");
   
    th2sel=acos(-(mselsq-msersq)/delta);
   

   /*printf("sel2th', %.3E  phil=%.3e\n",th2sel,phil1);	
*/
   E[1][1]= E[2][2]= cos(th2sel/2.);  
   E[1][2]=-sin(th2sel/2.)*(cos(phil1)-I*sin(phil1));
   E[2][1]=sin(th2sel/2.)*(cos(phil1)+I*sin(phil1));
   
   
  /* check the diagonalisation of the matrix*/
   md11= conj(E[1][1])*E[1][1]*mselsq+ conj(E[2][1])*E[2][1]*msersq+
   conj(E[1][1])*E[2][1]*conj(ze)+ conj(E[2][1])*E[1][1]*ze;
   md12= conj(E[1][1])*E[1][2]*mselsq+ conj(E[2][1])*E[2][2]*msersq+
   conj(E[1][1])*E[2][2]*conj(ze)+ conj(E[2][1])*E[1][2]*ze;
   md21= conj(E[1][2])*E[1][1]*mselsq+ conj(E[2][2])*E[2][1]*msersq+
   conj(E[1][2])*E[2][1]*conj(ze)+ conj(E[2][2])*E[1][1]*ze;
   md22= conj(E[1][2])*E[1][2]*mselsq+ conj(E[2][2])*E[2][2]*msersq+
   conj(E[1][2])*E[2][2]*conj(ze)+ conj(E[2][2])*E[1][2]*ze;
  

/* Neutralino contribution */
   mse1=findValW2("MSeL");
   mse2=findValW2("MSeR");
   for(i=1;i<=4;i++)
   { sprintf(namer,"MNE%d",i);
     mn[i]=findValW2(namer);
   }  

   for(i=1;i<=4;i++) for(j=1;j<=4;j++)
   { sprintf(namer,"Zn%d%dr",i,j);
     sprintf(namei,"Zn%d%di",i,j);
     N[i][j]=findValW2(namer)+I*findValW2(namei);
   }  
      
   for(i=1;i<=4;i++)
   {  y[i][1]= mse1*mse1/(mn[i]*mn[i]);
      y[i][2]= mse2*mse2/(mn[i]*mn[i]);
   
      for(a=1;a<=2;a++)
      {  nl[i][a]= (N[i][2]+tw*N[i][1])*E[1][a] - M_SQRT2*ye*N[i][3]*E[2][a];
         nr[i][a]=2*tw*N[i][1]*conj(E[2][a])+M_SQRT2*ye*(N[i][3]*conj(E[1][a]));  	  
         nlr[i][a]=cimag(nl[i][a]*nr[i][a]); 
      }
   } 
  /* take the formula of Choi et al*/   
   for(i=1,fneut=0;i<=4;i++)for(j=1;j<=2;j++)
   { fneut+=f2(y[i][j])/(mn[i])*nlr[i][j];
     /*printf("fneut, %E %E  \n",fneut,nlr[i][j]);*/
   }

   edm_nt=alphew/48./M_PI/sw/sw*fneut;
/* Conversion Gev-1->cm 1.97e-14 Gevcm*/	
   edm1loop=1.97e-14*(edm_ch+edm_nt);
   
   
/* Compute the two-loop Higgs contribution  including mixing
Reference Pilaftsis hep-ph/0207277 */

   for(i=1;i<=3;i++) 
   {  sprintf(namer,"Mh%d",i);
      mh[i]=findValW2(namer);
      for(j=1;j<=3;j++)
      { sprintf(namer,"Zh%d%d",i,j);
        Zh[i][j]=findValW2(namer);
      }  
   }
    

   mst1=findValW2("MSt1");
   mst2=findValW2("MSt2");
   At=Zafi(findValW2("aAt"),findValW2("fiAt"));
   Xt=At-conj(mu)/tb;
   msb1=findValW2("MSb1");
   msb2=findValW2("MSb2");
   Ab=Zafi(findValW2("aAb"),findValW2("fiAb"));
   Xb=Ab-tb*conj(mu);
   alpss=alphaQCD(sqrt(mst1*mst2));
   
   mbsusy=MbRun(sqrt(msb1*msb2));
   mtsusy=MtRun(sqrt(mst1*mst2));
   hbs=ee*mbsusy/(mw*sw*cb);
   hts=ee*mtsusy/(mw*sw*sb);

   
   for(edmh=0.,i=1;i<=3;i++)
   { xit[i]=ht*ht*sb*sb/2/(mst2*mst2-mst1*mst1)*(
         cimag(mu*At)*Zh[3][i]/sb
        -creal(mu*Xt)*Zh[1][i]
        +creal(conj(At)*Xt)*Zh[2][i]               )/sb;
     xib[i]=hb*hb*cb*cb/2/(msb2*msb2-msb1*msb1)*(
          cimag(mu*Ab)*Zh[3][i]/cb   
         -creal(mu*Xb)*Zh[2][i] 
         +creal(conj(Ab)*Xb)*Zh[1][i]             )/cb;
     zetat[i]=-ht*ht*sb*Zh[2][i]/2.;
     zetab[i]=-hb*hb*cb*Zh[1][i]/2.;
     xt1h[i]=mst1*mst1/mh[i]/mh[i];
     xt2h[i]=mst2*mst2/mh[i]/mh[i];
     xb1h[i]=msb1*msb1/mh[i]/mh[i];
     xb2h[i]=msb2*msb2/mh[i]/mh[i];
     edmh+= -(
     4/9.*(xit[i]*(F(xt1h[i])-F(xt2h[i])) +zetat[i]*(F(xt1h[i])+F(xt2h[i])))      
    +1/9.*(xib[i]*(F(xb1h[i])-F(xb2h[i])) +zetab[i]*(F(xb1h[i])+F(xb2h[i])))
             )/mh[i]/mh[i]*Zh[3][i]*tb;
   
  }


/* Top and bottom two-loop contributions*/

/* First calculate coefficients of effective Lagrangian for Hqq couplings*/
	
mgl=Zafi(findValW2("aM3"),findValW2("fiM3"));
    
deltahb=-2./3./M_PI*alpss*conj(mgl)*Ab*I2(msb1*msb1,msb2*msb2,creal(mgl*conj(mgl)))
-hts*hts/2./16./M_PI/M_PI*mu*conj(mu)*I2(mst1*mst1,mst2*mst2,creal(mu*conj(mu))); 
Deltahb=2./3./M_PI*alpss*conj(mgl*mu)*I2(msb1*msb1,msb2*msb2,creal(mgl*conj(mgl)))
+hts*hts/16./2./M_PI/M_PI*conj(At*mu)*I2(mst1*mst1,mst2*mst2,creal(mu*conj(mu))); 	 

deltaht=-2./3./M_PI*alpss*conj(mgl)*At*I2(mst1*mst1,mst2*mst2,creal(mgl*conj(mgl)))
-hbs*hbs/2./16./M_PI/M_PI*mu*conj(mu)*I2(msb1*msb1,msb2*msb2,creal(mu*conj(mu))); 
Deltaht=2./3./M_PI*alpss*conj(mgl*mu)*I2(mst1*mst1,mst2*mst2,creal(mgl*conj(mgl)))
+hbs*hbs/16./2./M_PI/M_PI*conj(Ab*mu)*I2(msb1*msb1,msb2*msb2,creal(mu*conj(mu))); 

/*
printf("deltahb %.3e %.3e \n",creal(deltahb),cimag(deltahb));
printf("Deltahb %.3e %.3e \n",creal(Deltahb),cimag(Deltahb));	
printf("deltaht %.3e %.3e \n",creal(deltaht),cimag(deltaht));
printf("Deltaht %.3e %.3e \n",creal(Deltaht),cimag(Deltaht));	
*/

   edmhq=0;  
   edmtq=0;
   edmbq=0;
   for(i=1;i<=3;i++)
   {
      xth[i]=mt*mt/mh[i]/mh[i];	
      xbh[i]=mb*mb/mh[i]/mh[i];

      gSbb[i]= Zh[1][i]/cb*creal( (1.+deltahb)/(1.+deltahb+tb*Deltahb))
              +Zh[2][i]/cb*creal( (Deltahb)/(1.+deltahb+tb*Deltahb))
              +Zh[3][i]*cimag(Deltahb*(1+tb*tb)/(1.+deltahb+Deltahb*tb));
      gPbb[i]=-Zh[3][i]*creal( ((1.+deltahb)*tb-Deltahb)/(1.+deltahb+tb*Deltahb))
	      +Zh[1][i]/cb*cimag( (Deltahb*tb)/(1.+deltahb+tb*Deltahb))
              -Zh[2][i]/cb*cimag(Deltahb/(1.+deltahb+Deltahb*tb));
	
/*	
      gStt[i]= Zh[2][i]/sb*creal( (1.+deltaht)/(1.+deltaht+Deltaht/tb))
	      +Zh[1][i]/sb*creal( (Deltaht)/(1.+deltaht+Deltaht/tb))
              +Zh[3][i]*cimag(Deltaht*(1+1./tb/tb)/(1.+deltaht+Deltaht/tb));
*/	
      double a= Zh[2][i]/sb*creal( (1.+deltaht)/(1.+deltaht+Deltaht/tb));
      double b= Zh[1][i]/sb*creal( (Deltaht)/(1.+deltaht+Deltaht/tb));
      double c= Zh[3][i]*cimag(Deltaht*(1+1./tb/tb)/(1.+deltaht+Deltaht/tb));
         
      gStt[i]=a+b+c;
      gPtt[i]=-Zh[3][i]*creal(((1.+deltaht)/tb-Deltaht)/(1.+deltaht+Deltaht/tb))
              +Zh[2][i]/sb*cimag( (Deltaht/tb)/(1.+deltaht+Deltaht/tb))
              -Zh[1][i]/sb*cimag(Deltaht/(1.+deltaht+Deltaht/tb));


	
	edmhq+=
           4/9.*(Zh[3][i]*tb*gStt[i]*f(xth[i])-Zh[1][i]/cb*gPtt[i]*g(xth[i]))
          +1/9.*(Zh[3][i]*tb*gSbb[i]*f(xbh[i])-Zh[1][i]/cb*gPbb[i]*g(xbh[i]));
	  
	  edmtq=
           4/9.*(Zh[3][i]*tb*gStt[i]*f(xth[i])-Zh[1][i]/cb*gPtt[i]*g(xth[i]));
	  
	  edmbq=
          +1/9.*(Zh[3][i]*tb*gSbb[i]*f(xbh[i])-Zh[1][i]/cb*gPbb[i]*g(xbh[i]));

/*printf("edmtq = %e edmbq= %e  \n",edmtq,edmbq);	
printf("edmtqS = %e edmtqP= %e  \n", 4/9.*(Zh[3][i]*tb*gStt[i]*f(xth[i])),
4/9.*(-Zh[1][i]/cb*gPtt[i]*g(xth[i]))   );	
printf("edmtqLO = %e edmbqLO= %e  \n",
4/9.*(Zh[3][i]*tb*Zh[2][i]/sb*f(xth[i])+Zh[1][i]/cb*Zh[3][i]/tb*g(xth[i])),
1/9.*(Zh[3][i]*tb*Zh[1][i]/cb*f(xbh[i])+Zh[3][i]*tb*Zh[1][i]/cb*g(xbh[i])));
	  
printf("gSbb[i] = %e gPbb[i]= %e  \n",gSbb[i],gPbb[i]);	
printf("gSbbLO[i] = %e gPbbLO[i]= %e  \n", Zh[1][i]/cb,-Zh[3][i]*tb);	
printf("gStt[i] = %e gPtt[i]= %e  \n",gStt[i],gPtt[i]);	
printf("gSttLO[i] = %e gPttLO[i]= %e  \n", Zh[2][i]/sb,-Zh[3][i]/tb  );	
*/
   
}
	
	

/* Chargino two-loop contribution*/
   edmhch=0;
  
   for(i=1;i<=3;i++)
   { xc1h[i]=mch1*mch1/mh[i]/mh[i];
     xc2h[i]=mch2*mch2/mh[i]/mh[i];
     /*printf("xc1%.3e  \n",xc1h[i]);
   */
      ac1[i]=Zh[1][i]*2.*creal(conj(Zu12)*Zv11)
           +Zh[2][i]*2.*creal(conj(Zu11)*Zv12)
	   +Zh[3][i]*2.*(sb*cimag(conj(Zu12)*Zv11)+cb*cimag(conj(Zu11)*Zv12));
	
      ac2[i]=Zh[1][i]*2.*creal(conj(Zu22)*Zv21)
           +Zh[2][i]*2.*creal(conj(Zu21)*Zv22)
	   +Zh[3][i]*2.*(sb*cimag(conj(Zu22)*Zv21)+cb*cimag(conj(Zu21)*Zv22));
	 
     bc1[i]=-Zh[1][i]*2.*cimag(conj(Zu12)*Zv11)
            -Zh[2][i]*2.*cimag(conj(Zu11)*Zv12)
	    +Zh[3][i]*2.*(sb*creal(conj(Zu12)*Zv11)+cb*creal(conj(Zu11)*Zv12));
	 
     bc2[i]=-Zh[1][i]*2.*cimag(conj(Zu22)*Zv21)
            -Zh[2][i]*2.*cimag(conj(Zu21)*Zv22)
	    +Zh[3][i]*2.*(sb*creal(conj(Zu22)*Zv21)+cb*creal(conj(Zu21)*Zv22));; 

    edmhch+=(-Zh[3][i]*tb*ac1[i]*f(xc1h[i])+Zh[1][i]/cb*bc1[i]*g(xc1h[i]))/mch1+
          (-Zh[3][i]*tb*ac2[i]*f(xc2h[i])+Zh[1][i]/cb*bc2[i]*g(xc2h[i]))/mch2;
	
	edm2loop=1.97e-14*(3.*alphem/32/M_PI/M_PI/M_PI*me*edmh+
	3.*alphem*alphem/8/M_PI/M_PI/sw/sw*me/mw/mw*edmhq
	-alphem*alphem/8/M_SQRT2/M_PI/M_PI/sw/sw*me/mw*edmhch);	

	
 /*printf("edm2loopsquark %.3e  \n", 1.97e-14*(3.*alphem/32/M_PI/M_PI/M_PI*me*edmh));

 printf("edmhq %.3e  \n", -alphem*alphem/8/M_SQRT2/M_PI/M_PI/sw/sw*me/mw*edmhch*1.97e-14);
 printf("edm2loopquark %.3e  \n", 1.97e-14*3.*alphem*alphem/8/M_PI/M_PI/sw/sw*me/mw/mw*edmhq);
*/
  



/* CS CP odd operator*/
   for(csg=0,i=1;i<=3;i++)
   {
     double ggg=2./3.*(gStt[i]+gSbb[i])+2.*mw*mw/ee/ee*sw*sw/3.*
  ( xit[i]*(1/(mst1*mst1)-1/(mst2*mst2))+zetat[i]*(1/(mst1*mst1)+1/(mst2*mst2))    
   +xib[i]*(1/(msb1*msb1)-1/(msb2*msb2))+zetab[i]*(1/(msb1*msb1)+1/(msb2*msb2))
  );
	csg+= ggg*Zh[3][i]/mh[i]/mh[i];	
   }
     	  
   }
	
   
   Cs=-0.1*tb*me*M_PI*alphem/sw/sw/mw/mw*csg;
   *el=edm1loop+edm2loop;
   *thallium=(Cs*8.5e-13)-585*(edm1loop+edm2loop); 

}
