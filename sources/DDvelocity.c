#include"micromegas_aux.h"
#include"micromegas_f.h"


/*===== Maxwell velosity distribution =====*/ 

double Maxwell(double v) 
{  double res,vsum,vdif,DV2;

   if(v<=0|| v>=Vesc+Vearth) return 0;

   double rv=Vesc/Vrot;
   double NR= pow(M_PI,1.5)*Vrot*Vrot*Vrot*(erf(rv)- 2*rv/sqrt(M_PI)*exp(-rv*rv));
             
   DV2=Vrot*Vrot;
   if(Vearth*v<0.001*DV2) res= 4*v*exp((-v*v-Vearth*Vearth)/(DV2))/(DV2);
   else 
   {   
     vsum=Vearth+v;
     if(vsum>Vesc) vsum=Vesc;
     vdif=Vearth-v;
     res=(exp(-vdif*vdif/DV2)-exp(-vsum*vsum/DV2))/Vearth;
   }
   return M_PI*Vrot*Vrot/NR*res;
}

static double erfi_int(double x){ return exp(x*x);} 
static double erfi( double x){ return 2/sqrt(M_PI)*simpson(erfi_int,0,x,1E-4,NULL);}

static double vR,vT,v_,sn; // integration parameters 
static double fiInt(double fi) { return exp( -pow(v_*sn*sin(fi)/vR,2) - pow(v_*sn*cos(fi)/vT,2));} // azimuth angle integrand 
                                          
static double cosInt(double cs) 
{  sn=sqrt(1-cs*cs);  return  simpson( fiInt,0,2*M_PI,1E-5,NULL)*exp(-pow((v_*cs-Vearth)/vT,2));} // polar anfle intergand
   

#define DIM 100

double SHMpp(double v) 
{ 
  static double vTab[DIM], SHMppTab[DIM] ;
  if(v<0 || v>Vearth+Vesc) return 0;
  double static beta=-1, eta=-1, vRot=-1,vE=-1,vEsc=-1;
  if( beta!=betaSHMpp || eta!=etaSHMpp ||   vEsc!=Vesc || vRot!=Vrot || vE!=Vearth)
  {  beta=betaSHMpp;
     eta=etaSHMpp;
     vEsc=Vesc;
     vRot=Vrot;
     vE=Vearth; 
     for(int i=0;i<DIM;i++) vTab[i]=i*(vE+vEsc)/(DIM-1);
     SHMppTab[0]=SHMppTab[DIM-1]=0;  
     double rv=vEsc/vRot;
     double NR= pow(M_PI,1.5)*vRot*vRot*vRot*(erf(rv)-     2*rv/sqrt(M_PI)*exp(-rv*rv));
     vR=vRot/sqrt(1-2./3.*beta), vT=vRot*sqrt((1-beta)/(1-2./3.*beta));
     double NS= pow(M_PI,1.5)*vR*vT*vT*(erf(vEsc/vR) - sqrt( (1/beta-1))*exp(-vEsc*vEsc/vT/vT)*erfi(vEsc/vR/sqrt(1/beta-1)));

     for(int i=1;i<DIM-1;i++)
     {  v_=vTab[i];  
       double vm= vEsc>v_+vE ? v_+vE : vEsc;
       double ediff; 
       if(vE*v_<0.001*vRot*vRot) ediff= 4*v*exp((-v_*v_-vE*vE)/vRot/vRot);
       else  ediff=exp(-pow((v_-vE)/vRot,2))-exp(-pow(vm/vRot,2));       
       SHMppTab[i]= (1-etaSHMpp)*M_PI*vRot*vRot/vE/NR*ediff;
              vm= vEsc<v_+vE ? v_+vE : vEsc;
       SHMppTab[i]+= etaSHMpp* v_*simpson(cosInt,-1+(vm*vm-vEsc*vEsc)/(2*vE*v_),1,1E-5,NULL)/NS;       
     }
  } 
  return polint3(v,DIM,vTab, SHMppTab);
}  
  