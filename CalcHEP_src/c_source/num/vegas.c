#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include "vegas.h"
#include "crt_util.h"
#include "drandXX.h"

#define MAX_DIM   15
#define MAX_NDMX  50

static void drand_arr(int dim, double * x)
{ int i;
  unsigned int umax= UINT_MAX;
  int rest=dim;
  unsigned int randpos=umax*drandXX();
  
  for(i=0;i<dim;i++) x[i]=-1; 
  
  for(i=0;i<dim;i++)
  {
     int pos=randpos%rest;
     int j;
     for(j=0;;j++) if(x[j]<0)
     {if(pos)pos--;else {x[j]=drandXX();break;}}
     randpos/=rest;
     umax/=rest;
     rest--;
     if(umax<rest){ umax= UINT_MAX;randpos=umax*drandXX();}
  }
}


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
/*printf("amotry returns fac=%f %E\n",fac,  ytry);   */
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



static int Kg[MAX_DIM], Ng[MAX_DIM];
static int Ndim;
static double (*f_)(double*, double);
static int Ndmx;
static double *Xgrid;
static double *Cgrid;

#define XG(j,i) Xgrid[(i)+(j)*(Ndmx)]


static int nroot(long N, int n)
{  int i,r;
   long N_;
   if(n==1) return N;
   r=pow(N,1./n);
   
   for(i=1,N_=r; i<n; i++) N_*=r; 
   
/*   printf("N_0=%d\n",N_);*/

   for(; N_<N; ) { r++; for(i=1,N_=r; i<n; i++) N_*=r;}
   for(; N_>N; ) { r--; for(i=1,N_=r; i<n; i++) N_*=r;} 
   return r;   
}
static void  generateVegasCubs(vegasGrid * vegPtr,long * nCubs) 
{  int i;
   double nCubs_=*nCubs;

   Ndim=vegPtr->ndim;
   Ndmx = vegPtr->ndmx+1;
   Xgrid= vegPtr->x_grid;
   Cgrid= vegPtr->c_grid;
   *nCubs=1;
   for(i=0;i<Ndim;i++)
   {
      Ng[i]   =nroot(nCubs_,Ndim -i);
      nCubs_ /= Ng[i];
      *nCubs *= Ng[i];
      Kg[i]=0;
   }
}



static void  Local2Global(double *XLOC,double * XGLOB,double * JACOB, 
int *GRID_LOC)
{ int j,n;
  double xlj,xn,xn_;                               
  *JACOB=1;                                         
  for (j = 0; j < Ndim; ++j)                       
  {  xlj = (Kg[j ] + XLOC[j])/Ng[j];               
     n=(int)(xlj*(Ndmx-1));                              
     if (n) xn_= XG(j,n);else  xn_=0;               
     xn = XG(j,n+1);                                 
     XGLOB[j] = xn_ +(xn-xn_)*(xlj*(Ndmx-1)-n);          
     *JACOB *= (xn-xn_)*(Ndmx-1);                         
     if(GRID_LOC) GRID_LOC[j] = n;                 
  }                                                
}


vegasGrid *  vegas_init(int dim,int nd)
{ 
   vegasGrid * vegPtr;

   if((dim>MAX_DIM)||(nd>MAX_NDMX)) return NULL;
   vegPtr=(vegasGrid * )malloc(sizeof(vegasGrid));
   if(vegPtr)
   { 
      vegPtr->ndmx=nd;
      vegPtr->ndim = dim;
      vegPtr->x_grid=malloc(dim*(nd+1)*sizeof(double));
      vegPtr->c_grid=malloc(dim*(nd+1)*sizeof(double));
      
      Xgrid=vegPtr->x_grid;
      Ndmx = vegPtr->ndmx+1;
      if(vegPtr->x_grid && vegPtr->c_grid)
      { int i, j;
        double * x_=vegPtr->x_grid;
        double * c_=vegPtr->c_grid;
        for(j=0;j<dim;j++) for(i=0;i<=nd;i++,x_++,c_++) 
        {  *x_=i/(double)nd;
           *c_=1./nd;
        }
      }else 
      { 
        if(vegPtr->x_grid) free(vegPtr->x_grid);
        if(vegPtr->c_grid) free(vegPtr->c_grid); 
        free(vegPtr); 
        return NULL;
      }
      vegPtr->nCubs=0;
      vegPtr->fMax=NULL;
   }
   return vegPtr;
}


void vegas_finish(vegasGrid * vegPtr) 
{   if(vegPtr)
    {  free(vegPtr->x_grid); 
       free(vegPtr->c_grid); 
       if(vegPtr->fMax) free(vegPtr->fMax);
       free(vegPtr); 
       vegPtr=NULL;
    }
}

/*     			*  VEGAS  *
      SUBROUTINE PERFORMS NDIM-DIMENSIONAL MONTE CARLO INTEG'N 
      - BY G.P. LEPAGE    SEPT 1976/(REV)AUG 1979 
      - ALGORITHM DESCRIBED IN J COMP PHYS 27,192(1978) 
*/

int vegas_int(vegasGrid * vegPtr, long ncall0, double alph, 
 double (*fxn)( double *,double), double *ti, double *tsi)
{
   int dim= vegPtr->ndim;

   double *d=malloc((vegPtr->ndmx+1)*dim*sizeof(double));
#define DD(j,i) d[(i)+(j)*Ndmx]

   double x[MAX_DIM];
   double xlocal[MAX_DIM];
   int    ia[MAX_DIM];
   
   int i,j;
   double  f2,fb,f2b;
   int  npg=2;
   long nCubs=ncall0/npg;
   long cCub;
   int ret_code=0;
   generateVegasCubs(vegPtr,&nCubs);

   npg=ncall0/nCubs;
    
   *ti  = 0;
   *tsi = 0;
   for (j = 0; j < dim; ++j) { for (i = 0; i < Ndmx-1; ++i)  DD(j,i) = 0;}

/*    - MAIN INTEGRATION LOOP */
   for(cCub=0; cCub<nCubs; cCub++) 
   {  
      if(informline(cCub,nCubs))  { ret_code=1; goto exi;}

      fb  = 0;
      f2b = 0;
      
      for(i = 0; i<npg; i++) 
      {   double f;
          drand_arr(dim,xlocal);
          Local2Global(xlocal,x, &f, ia);
          f *= (*fxn)(x,f);
          fb += f;
          f2= f*f;
          f2b += f2;
          for (j = 0;j<dim;++j) DD(j,ia[j]) += f2;
      }
       
      f2b = sqrt(f2b/npg);
      fb /=npg;
       
      f2b = (f2b - fb) * (f2b + fb)/(npg-1);             
       
      *ti  += fb/nCubs;
      *tsi += f2b/((double)nCubs * (double)nCubs);   
    
      for(i=dim-1; i>=0; i--){if(++Kg[i]<Ng[i]) break; else Kg[i]=0;}
    }
    
    *tsi = sqrt(fabs(*tsi));

    if(*tsi < 1.E-6*fabs(*ti)) *tsi=1.E-6*fabs(*ti);

    if(alph>0)  /* REFINE GRID */
    { 
        double r[MAX_NDMX]; 					       
        double dt[MAX_DIM];
        double xin[MAX_NDMX];
        
        double  xn, xo,dr;
        int k,ndm= Ndmx - 2;
        
        for (j = 0; j < dim ; ++j)
        {
            xo = DD(j,0);
            xn = DD(j,1);
            DD(j,0) = (xo + xn) / 2;
            dt[j] = DD(j,0);
            for (i = 1; i < ndm; ++i)
            {
                DD(j,i) = xo + xn;
                xo = xn;
                xn = DD(j,i+1);
                DD(j,i) = (DD(j,i) + xn) / 3;
                dt[j] += DD(j,i);
            }
            DD(j,ndm) = (xn + xo) / 2;
            dt[j] += DD(j,ndm);
        }
        for (j = 0; j < dim; ++j)
        {  double rc = 0;
    	   for (i = 0; i <= ndm; ++i)
    	   {
    	      r[i] = 0;
    	      if (DD(j,i) > 0)
    	      {  double xoln = log(dt[j]/DD(j,i));
    	         if (xoln <= 70.f)  r[i] = pow( (1 - exp(-xoln))/xoln, alph);
    	         else               r[i] = pow(  1/xoln,               alph);
    	      }
    	      rc += r[i];  
    	   }

    	   rc /= (Ndmx-1);
           if(rc)
	   { 
             for(i=0,k=0,xn=0,dr=0;i<ndm;) 
       	     {
                do
                {  dr += r[k];
    	           xo = xn;
    	           xn = XG(j,k+1);
    	           k++;
    	         } while (rc > dr);
    	         do
    	         {  dr -= rc;
    	            xin[i] = xn-(xn-xo)*dr/r[k-1];
                    i++;
                 } while (rc<=dr);
    	      }
    	      for (i=0;i<ndm;++i)  XG(j,i+1) = xin[i];
    	      XG(j,ndm+1) = 1;
    	      XG(j,0)=0;
           }
        }
        if(vegPtr->fMax){free(vegPtr->fMax);vegPtr->fMax=NULL;vegPtr->nCubs=0;}
     }
exi:  free(d); return ret_code;
#undef DD
} 

#define INCUB(x)   ((x)>0   ? (0.5+(x))/((x)+1.) : 0.5/(1.-(x)))
#define OUTCUB(x)  ((x)>0.5 ? (0.5-(x))/((x)-1.) : 1-0.5/x)

static double f_max(double* x)
{
  int i;
  double f;
  double xg[MAX_DIM];
  double xl[MAX_DIM];
  double tmp;

  for(i=0;i<Ndim;i++) xl[i]=INCUB(x[i]);
  Local2Global(xl,xg, &f,NULL);

  tmp=f_(xg,f);
  f*=tmp;
  
  if(f<0) return -f; else return f;
}


static double run_amoeba(int ndim, double *xx, double *y, double step, double eps, int nCalls)
{
   int i,j;
   for (i=1; i<=ndim;++i)
   { double * x_i=xx+i*ndim; 
     for(j=0;j<ndim;++j) x_i[j]=xx[j];
     if(x_i[i-1] >0.5)x_i[i-1] -= step; else x_i[i-1] += step; 
     for(j=0;j<ndim;j++) x_i[j]=OUTCUB(x_i[j]);
     y[i]=f_max(x_i);
   }

for(j=0;j<ndim;j++) xx[j]=OUTCUB(xx[j]);

/*printf(" test y[0]: %E %E \n",y[0], f_max(xx));*/

   {/* int nn=nCalls;*/
      double r=amoeba(xx,y,ndim,f_max,eps,&nCalls);
/*      printf("try: %d ",nn-nCalls); */
      return r; 
   }
}


int vegas_max(vegasGrid * vegPtr, long  nCubs, long nRandom, long nSimplex,
 double (*fxn)(double *,double), double * eff)
{
   int dim= vegPtr->ndim;
   
   double x[MAX_DIM]; 
   double xlocal[MAX_DIM];
   double *xx=malloc((dim+2)*dim*sizeof(double));
   double *y=malloc((dim+2)*sizeof(double));
   double average =0;
   double smax;
   float *fmax;

   int i;
   long  cCub;
   int ret_code=0;
   
   if(nRandom<=0) nRandom=2; 
   generateVegasCubs(vegPtr,&nCubs);
   if(vegPtr->fMax) free(vegPtr->fMax);
   vegPtr->nCubs=nCubs;    
   vegPtr->fMax=malloc(nCubs*sizeof(float));
   fmax=vegPtr->fMax;
   f_=fxn;

   for(cCub=0; cCub<nCubs; cCub++) 
   {  double cMax=0; 
      long k; 
      double f;

      if(informline(cCub,nCubs)) {ret_code=1; goto exi;}
/* if(cCub==36) nRandom*=1000;  */    
      for(k = 0; k<nRandom ; k++) 
      {  
         drand_arr(dim,xlocal);
         Local2Global(xlocal,x, &f, NULL);f *=(*fxn)(x,f); if(f<0) f=-f;
         average+=f;
         if(f>cMax){cMax=f;  for(i=0;i<dim;i++)xx[i]=xlocal[i]; y[0]=cMax; }
      }
/* if(cCub==36) nRandom/=1000; */
{
/* double m=cMax;*/ 
/*printf("L=%d max1=%E\n",cCub,cMax);*/
      if(cMax&&nSimplex) cMax=run_amoeba(dim, xx, y, 0.1 , 0.01, nSimplex); 
/*printf("L=%d max1=%E   max2=%E   r=%f %d %d %d \n",cCub,m, cMax, cMax/m,Kg[0],Kg[1],Kg[2]);*/
}
      fmax[cCub]=cMax;
      for(i=dim-1; i>=0; i--){if(++Kg[i]<Ng[i]) break; else Kg[i]=0;}
   }

   { double sum=0;
     float *p=fmax+nCubs;
     while(p!=fmax) sum+=*(--p);
     smax=sum;
/*     if(milk)
     { sum*=milk/nCubs;
       p=fmax+nCubs;
       while(p!=fmax) if(*(--p)<sum) *p=sum;
     }
*/
   }
   *eff=(average/nRandom)/smax;
exi:    
   free(xx); free(y); return ret_code;
}


int vegas_events(vegasGrid * vegPtr,  long  nEvents, double gmax, double milk,
int nSimplex,   double (*fxn)( double *,double), 
   void (*out)(long,int,char*),
   double * eff, double * rmax, double * mult, double * neg)
{
   int dim= vegPtr->ndim;

   double x[MAX_DIM];
   double xlocal[MAX_DIM];
   long  cEvent;
   long  cCub;
   double sum=0;
   long L0;     
   long nCubs=vegPtr->nCubs;
   float * smax,*fmax_,*fmax;
   double *xx,*y;
   int i;
   int ret_code=0;


   if(!(vegPtr->fMax)) return -1;
   fmax=vegPtr->fMax;
   smax=malloc(nCubs*sizeof(float));
   fmax_=malloc(nCubs*sizeof(float));

   xx=malloc((dim+2)*dim*sizeof(double));
   y=malloc((dim+2)*sizeof(double));

   f_=fxn;
   
   generateVegasCubs(vegPtr,&nCubs);
   
   *eff=0;
   *rmax=0;
   *mult=0; 
   *neg=0;
   
   for(cCub=0;cCub<nCubs;cCub++)fmax_[cCub]=fmax[cCub];
   if(milk)
   { double minFmax;
     for(sum=0,cCub=0; cCub<nCubs; cCub++) sum+=fmax_[cCub];
     minFmax=sum*milk/nCubs;
     for(cCub=0; cCub<nCubs; cCub++)
     if(fmax_[cCub]<minFmax) fmax_[cCub]=minFmax;
   }
   for(sum=0,cCub=0; cCub<nCubs; cCub++) 
   {  sum+=fmax_[cCub];
      smax[cCub]=sum;
   }
   for(cCub=0; cCub<nCubs; cCub++) smax[cCub]/=sum; 

   for(cEvent=0; cEvent<nEvents; ) 
   {  long L;
      double f;
      int n,sgn=0;
      int rebuild=0;
      char *drandXXstate ;      
      double newMax;

      { double rc=drandXX();
        long L0=0;
        long L1=nCubs-1;
        while(L0+1 < L1)     
        {  L=(L0+L1)/2;
           if(smax[L]<=rc) L0=L;  else L1=L;  
        }  
        if(smax[L0]>rc) L=L0; else L=L1;
      }

 
      if(informline(cEvent,nEvents)) { ret_code=1; goto exi;}

      L0=L; for (i =dim-1;i>=0; i--) {Kg[i]=L0%Ng[i]; L0=L0/Ng[i];}
/*
 printf("L=%d dim=%d Ng={%d %d %d}, Kg={%d %d %d}\n",
        L,dim, Ng[0],Ng[1],Ng[2],Kg[0],Kg[1],Kg[2]);      
*/
      drandXXstate=seedXX(NULL);                                               
      drand_arr(dim,xlocal);                              
      Local2Global(xlocal,x,&f,NULL);f*=(*fxn)(x,f);if(f<0){f=-f;sgn=1;}      
      (*eff)++;     
      if(f>fmax[L])                                                           
      {  newMax=f;                                                            
         fmax[L]=newMax;                                                      
         if(newMax>gmax*fmax_[L])                                             
         {  rebuild=1;                                                        
            if(nSimplex)                                                     
            {                                                                 
               for(i=0;i<dim;i++)xx[i]=xlocal[i];                             
               y[0]=newMax; 
/*printf("L=%d old=%.2E new1=%.2E ",L,fmax_[L],newMax); */                                                  
               newMax=run_amoeba(dim,xx,y,0.1,0.01,nSimplex);
/*printf(" new2=%.2E  r=%f %d %d %d\n",newMax, newMax/fmax_[L],Kg[0],Kg[1],Kg[2]);*/
               fmax[L]=newMax;                                                
            }                                                                 
         }                                                                    
      }                                                                       
      f/=fmax_[L];                                                            
      if(f>*rmax) *rmax=f;                                                    
      f/=gmax;                                                                
      n=(int)(f);                                                             
      f=f-n;                                                                  
      if(f>drandXX()) n++;                                                    
      if(n)                                                                   
      {                                                                       
         cEvent+=n;                                                           
         if(cEvent>nEvents) n-=(cEvent-nEvents);                              
         *mult+=n-1;                                                          
         if(sgn) {n*=-1; *neg +=n;}                                           
         (*out)(L,n,drandXXstate);                                            
      }                                                                       


      if(rebuild)
      {  sum=0;
         for(cCub=0; cCub<nCubs; cCub++)
         {  sum+=fmax[cCub];
            fmax_[cCub]=fmax[cCub];
            smax[cCub]=sum;
         }
         for(cCub=0; cCub<nCubs; cCub++) smax[cCub]/=sum;
         rebuild=0;                     
         be_be();
      }
   }
     
exi:    

   *neg=(*neg)/cEvent;
   *mult=(*mult)/cEvent;
   *eff=nEvents/(*eff);
   free(xx); free(y); free(smax);free(fmax_);return ret_code;
   
} /* vegas_ */


/*
int unpak_events(vegasGrid * vegPtr, long  nCubs, long  nEvents)
{
   int dim= vegPtr->ndim;

   double x[MAX_DIM];
   double xlocal[MAX_DIM];
   long  cEvent;

   int i;

   generateVegasCubs(vegPtr,&nCubs);

   for(cEvent=0; cEvent<nEvents;cEvent++ ) 
   {  long L;
      int w;
      char drandXXstate[100];      
      double f;
      
      if(EOF==fscanf(stdin, "   ",drandXXstate,&L,&w)) break;
       
      for (i =dim-1;i>=0; i--) {Kg[i]=L%Ng[i]; L=L/Ng[i];}

      drand_arr(dim,xlocal);                               
      Local2Global(xlocal,x,&f,NULL); 
      mkmom(x, &f);
      lorrot(rapidity,nin_+nout_);
             
   }
   
}
*/
