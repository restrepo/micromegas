/*
 Copyright (C) 1997 by Alexander Pukhov, email: pukhov@theory.npi.msu.su 
*/
#define real double
#include "chep_crt.h"
#include "procvar.h"
#include "n_compil.h"
#include "out_int.h"
#include "out_ext.h"
#include "pvars.h"

double DP[15]; /* not used */


static double  calcpolynom(varptr  q)
{ double        s=0, m;
  short *  i; 

   if (q == NULL) return 0.0;
   do
   {
      m = q->coef->rval;
      i=q->vars;
      while(*i) {m*=vararr[*i].tmpvalue; i++;} 
      if (q->sgn == '+')  s += m;  else  s -= m;
      q = q->next;
   }  while (q);
   return s;
}



static void  calcconstants(void)
{infoptr  i_ptr;

   i_ptr = info;
   while (i_ptr)
   {
      if (i_ptr->consttype == expr) i_ptr->rval = calcpolynom(i_ptr->const_);
		else                i_ptr->rval = i_ptr->ival;
      i_ptr = i_ptr->next;
   }
}
 

static double  calconediagr(prgcodeptr  prg_,   denlist denomi, 
                                           double * momenta, int * err)
{ int     i;
 double   r, rr;

   rr = calcpolynom(prg_->totd);
   if (rr == 0)
   {  *err = 2;   /* Division by zero */
      return 0;
   }
  
   r = calcpolynom(prg_->rnum)*calcpolynom(prg_->totn)/rr;

   for (i = 0; i < prg_->denorno; i++)
   { denlist den=denomi;
      while( prg_->order_num[i] !=den->order_num 
          || prg_->width[i]!=den->width ) den=den->next;    
      if (prg_->power[i] ==2) r*=den->val2;else r*=den->val1;
   }
   return r;
}

static double  matrixelement(prgcodeptr  prg, denlist denomi, 
                                            double * momenta, int *err)
{double        sum=0.0; 
 denlist den=denomi;

   while (den)
   { double rr, irr;

      rr = sqrMom(den->momStr,momenta);
      irr=vararr[den->mass].tmpvalue;
      rr=irr*irr - rr;
      
      irr=irr*vararr[den->width].tmpvalue;

      if (fabs(rr) + fabs(irr) == 0)
      {  *err = 2;    /* Division by zero */
         return 0;
      }
      if(den->width)
      {  den->val2=1.0/(rr*rr+irr*irr);
         den->val1=rr*den->val2;
         den->val0=rr*den->val1;            
      }else
      {
         den->val0=1;
         den->val1=1/rr;
         den->val2=den->val1/rr;
      }
      den=den->next;
   }
   while (prg) 
   {  
     sum += calconediagr(prg,denomi,momenta,err); 
     if (*err) return 0.0; 
     prg= prg->next; 
   } 
/*   if (escpressed()) *err = 10;  */
   return sum; 
} 

static void sprod_(double*);
static double smpl(int nsub,long double*momenta_,int * err)
{
 int i,j;
 double momenta[4*MAXINOUT];

 for(i=0;i<4*(nin+nout);i++) momenta[i]=momenta_[i]; 

 for(i=0;  i<nin+nout-1;i++)
 for(j=i+1;j<nin+nout;j++)  vararr[scalarProductPos(i+1,j+1)].tmpvalue=momenta[4*i]*momenta[4*j]-
  momenta[4*i+1]*momenta[4*j+1]-momenta[4*i+2]*momenta[4*j+2]-momenta[4*i+3]*momenta[4*j+3];
 
 if(calcCoef[1]) { calcconstants(); calcCoef[1]=0;}

 return  matrixelement(allcanal[nsub-1].codeptr,    
                       allcanal[nsub-1].denominators,  momenta, err);
}
#include"sqme0.c"
