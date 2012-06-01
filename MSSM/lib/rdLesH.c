#include<math.h>
#include"pmodel.h"
#include"pmodel_aux.h"
#include"pmodel_f.h"
#include"../../sources/micromegas.h"



void CheckNCsector(double *zero_, double*m1_, double*m2_,double*mu_, double*tb_,
double *mz_, double* sw_)
{
 double Z[4][4], M[4],MM[4][4];
 char name[20];
 int i,j,k,l;
 double mZero,mTot;
 double w,tw,cw,sw,b,tb,cb,sb,mz;

 
 for(i=0;i<4;i++) 
 {  sprintf(name,"MNE%d",i+1);  
    M[i]=findValW(name);
 }

 for(i=0;i<4;i++) for(j=0;j<4;j++)
 {  sprintf(name,"Zn%d%d",i+1,j+1); 
    Z[j][i]=findValW(name);
 }

 mTot=0; 

 for( i=0;i<4;i++)
 { double s;

   for(l=0;l<4;l++)
   { for(s=0,k=0;k<4;k++) s+=Z[i][k]*M[k]*Z[l][k];
     MM[i][l]=s; 
     if(fabs(s)>mTot) mTot=fabs(s);
   }
 }

 mZero=0; 
 if(mZero<fabs(MM[0][1])) mZero=fabs(MM[0][1]);
 if(mZero<fabs(MM[2][2])) mZero=fabs(MM[2][2]);
 if(mZero<fabs(MM[3][3])) mZero=fabs(MM[3][3]);

 if(zero_) *zero_=mZero;


 tw=-MM[0][2]/MM[1][2];
 w=atan(tw);
 sw=sin(w);
 cw=cos(w);
 if(sw_) *sw_=sw;
 tb=-MM[0][3]/MM[0][2];
 if(tb_) *tb_=tb;
 b=atan(tb);
 sb=sin(b);
 cb=cos(b);
 mz= -MM[2][0]/cb/sw;
 if(mz_) *mz_=mz;

 if(m1_) *m1_=MM[0][0];
 if(m2_) *m2_=MM[1][1];
 if(mu_) *mu_=-MM[3][2]; 
}


int readLesH(char *fname, int mode)
{ /* mode 0 - EWSB;  1 - SUGRA/AMBS; 2 - FILE */
  int err;  
  err=slhaRead(fname,0);
  if(err) return err;




#ifdef IMPROVED
{ 
  double zero,m1,m2,mu,tb,mz,sw;                                                                                   
  CheckNCsector(&zero,&m1,&m2,&mu,&tb,&mz,&sw);
  printf("IMPROVED!\n");   
  printf("zero=%e  m1=%e ,m2=%E,mu=%E ,tb=%E,mz=%E,sw=%E\n",
    zero,m1,m2,mu,tb,mz,sw);

  assignValW("mu",mu);
  assignValW("tb",tb); 
  assignValW("MZ",mz);
  assignValW("SW",sw);
}
#endif

  return 0;
}

int lesHinput(char * fname)
{ int err= readLesH(fname, 2);
  if(err==0) FillVal(2);
  return err;
}
