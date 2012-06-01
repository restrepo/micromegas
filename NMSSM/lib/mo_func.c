#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include"../../sources/micromegas.h"
#include"../../sources/micromegas_aux.h"
#include"pmodel.h"
#include"pmodel_aux.h"

int readVarNMSSM(char * fname)
{
  char*vlist[36]={"alfEMZ","alfSMZ","McMc","MbMb","Mtp","tb","MG1","MG2","MG3",
  "Ml2","Ml3","Mr2","Mr3","Mq2","Mq3","Mu2","Mu3","Md2","Md3","At","Ab","Al",
  "Au","Ad","Am","mu","LambdQ","KappaQ","aLmbdQ","aKappQ","MZ","Mm","Ml","wt",
  "wZ","wW"};
  return readVarSpecial(fname,36,vlist);
} 

void o1Contents(FILE * f)
{ double val; 
  int err; 
 
  fprintf(f,"\n~o1 = ");
  err=findVal("Zn11",&val);
  if(err==0) fprintf(f,"%.3f*",val); else  fprintf(f,"???*");  
  fprintf(f,"bino");

  err=findVal("Zn12",&val);
  if(err==0) fprintf(f,"%+.3f*",val); else  fprintf(f,"???*");  
  fprintf(f,"wino");
  
  err=findVal("Zn13",&val);
  if(err==0) fprintf(f,"%+.3f*",val); else  fprintf(f,"???*");  
  fprintf(f,"higgsino1");
  
  err=findVal("Zn14",&val);
  if(err==0) fprintf(f,"%+.3f*",val); else  fprintf(f,"???*");  
  fprintf(f,"higgsino2");
  
  err=findVal("Zn15",&val);
  if(err==0) fprintf(f,"%+.3f*",val); else  fprintf(f,"???*");  
  fprintf(f,"singlino\n");
}


static void FillVal(int mode)
{ char name[10];
  int i,j,k;
  char* Zf[3]={"Zb","Zt","Zl"};
  char* Qmix[3]={"SBOTMIX","STOPMIX","STAUMIX"};
  double Pa[2][3]={{0,0,0},{0,0,0}};
  double Q;
  
  char* massName[35]={"Mha","Mhb","Mh1","Mh2","Mh3","MHc", "MNE1", "MNE2", "MNE3", "MNE4", "MNE5",    "MC1",  "MC2",  "MSG", "MSne", "MSnm", "MSnl", "MSeL", "MSeR", "MSmL", "MSmR", "MSl1", "MSl2", "MSdL", "MSdR", "MSuL", "MSuR", "MSsL", "MSsR", "MScL", "MScR", "MSb1", "MSb2", "MSt1","MSt2"};
  int   massId[35]  ={ 36  ,  46,  25,   35,   45,     37,1000022,1000023,1000025,1000035, 1000045, 1000024,1000037,1000021,1000012,1000014,1000016,1000011,2000011,1000013,2000013,1000015,2000015,1000001,2000001,1000002,2000002,1000003,2000003,1000004,2000004,1000005,2000005,1000006,2000006};   

  char * softName[13]={"MG1","MG2","MG3","Ml2","Ml3","Mr2","Mr3","Mq2","Mq3","Mu2","Mu3","Md2","Md3"};
  int      softId[13]={   1 ,   2 ,   3 ,  32 ,  33 ,  35 ,  36 ,  42 ,  43 ,  45 ,  46 ,  48 ,  49};

  char * mohName[11]={"Lambda","Kappa","aLmbd0","dMb","vev","Mh1R","Mh2R","Mh3R","MhaR","MhbR","MHcR"};
  int    mohId[11]  ={      1 ,     2 ,      3 ,   4 ,   5 ,   25 ,   35 ,   45 ,   36 ,   46 , 37 };

  for(i=0;i<35;i++) assignValW(massName[i],slhaVal("MASS",0.,1,massId[i]));
  Q=sqrt(fabs(findValW("MSt1")*findValW("MSt2"))); 
  for(i=0;i<11;i++) assignValW(mohName[i],slhaVal("MO_HIGGS",Q,1,mohId[i]));
  
  for(i=1;i<=5;i++) for(j=1;j<=5;j++) 
  { sprintf(name,"Zn%d%d",i,j); assignValW(name,slhaVal("NMNMIX",Q,2,i,j));}
  for(i=1;i<=3;i++) for(j=1;j<=3;j++) 
  { sprintf(name,"Zh%d%d",i,j); assignValW(name,slhaVal("NMHMIX",Q,2,i,j));}

  for(i=1;i<=2;i++) for(j=1;j<=3;j++) Pa[i-1][j-1]=slhaVal("NMAMIX",Q,2,i,j);
  { double tb;
    assignValW("Pa12",Pa[0][2]);
    assignValW("Pa22",Pa[1][2]);
    if(Pa[0][1]==0)  assignValW("Pa11",0); else
    {
      tb=Pa[0][0]/Pa[0][1];
      assignValW("Pa11",Pa[0][1]*sqrt(1+tb*tb));
    }
    if(Pa[1][1]==0)  assignValW("Pa21",0); else
    {
      tb=Pa[1][0]/Pa[1][1];
      assignValW("Pa21",Pa[1][1]*sqrt(1+tb*tb));
    }
  }
  
  for(i=1;i<=2;i++) for(j=1;j<=2;j++)
  { sprintf(name,"Zu%d%d",i,j);assignValW(name,slhaVal("UMIX",Q,2,i,j));
    sprintf(name,"Zv%d%d",i,j);assignValW(name,slhaVal("VMIX",Q,2,i,j));
  }
   
  for(k=0;k<3;k++)
  { double M[3];
    for(i=1;i<=2;i++) 
    {M[i]=slhaVal(Qmix[k],Q,2,1,i);
     sprintf(name,"%s1%d",Zf[k],i); assignValW(name,M[i]);
    }  
    M[2]*=-1;
    for(i=1;i<=2;i++) 
    { sprintf(name,"%s2%d",Zf[k],i); assignValW(name,M[3-i]); }  
  }   
  
  if(mode>0)
  {  
    for(i=0;i<13;i++) assignValW(softName[i],slhaVal("MSOFT",Q,1,softId[i])); 
    assignValW("mu", slhaVal("HMIX",Q,1,1));    
    assignValW("Al", slhaVal("Ae",Q,2,3,3));
    assignValW("Ab", slhaVal("Ad",Q,2,3,3));
    assignValW("At", slhaVal("Au",Q,2,3,3));
    assignValW("Am", slhaValExists("Ae",2,2,2)>0 ? slhaVal("Ae",Q,2,2,2):slhaVal("Ae",Q,2,3,3));
    assignValW("Au", slhaValExists("Au",2,2,2)>0 ? slhaVal("Au",0.,2,2,2):slhaVal("Au",0.,2,3,3));
    assignValW("Ad", slhaValExists("Ad",2,2,2)>0 ? slhaVal("Ad",0.,2,2,2):slhaVal("Ad",0.,2,3,3));
    assignValW("tb",slhaVal("MINPAR",Q,1,3) );
  }
  
  if(mode==2)
  {    assignValW("alfSMZ",slhaVal("SMINPUTS",Q,1,3) );
    assignValW("MbMb",  slhaVal("SMINPUTS",Q,1,5) );
    assignValW("Mtp",   slhaVal("SMINPUTS",Q,1,6) );
    assignValW("Ml",    slhaVal("SMINPUTS",Q,1,7) );
  }   

}


static void CharginoZM(void)
{
  double M2,mu,g2v1,g2v2,offQ,TrX2,detX,D,tU,tV,CU,SU,CV,SV;  
  double Zn_[5][5],NMassM[5][5],M[5];
  double mc[2],Zv_[2][2],Zu_[2][2];
  char name[10];
  int i,j,k;

  for(i=1;i<=5;i++) for(j=1;j<=5;j++) 
  { sprintf(name,"Zn%d%d",i,j); Zn_[i-1][j-1]=findValW(name);}  
  for(i=1;i<=5;i++) { sprintf(name,"MNE%d",i); M[i-1]=findValW(name);} 
  for(i=0;i<5;i++) for(j=0;j<5;j++) for(k=0,NMassM[i][j]=0;k<5;k++)
  NMassM[i][j]+=Zn_[k][i]*M[k]*Zn_[k][j];

   M2=NMassM[1][1]; 
   mu=-NMassM[2][3];
   g2v1= -NMassM[1][3]*sqrt(2.);
   g2v2=  NMassM[1][2]*sqrt(2.);
        
   offQ =g2v1*g2v1+g2v2*g2v2;
   TrX2 =offQ +M2*M2+mu*mu;
   detX =mu*M2 - g2v1*g2v2;
   D=TrX2*TrX2 - 4.*detX*detX;

   tU=(g2v2*g2v2-g2v1*g2v1-M2*M2+mu*mu-sqrt(D))/2./(M2*g2v2+mu*g2v1);
   tV=(g2v1*g2v1-g2v2*g2v2-M2*M2+mu*mu-sqrt(D))/2./(M2*g2v1+mu*g2v2);

   CU=cos(atan(tU));
   SU=sin(atan(tU));
   CV=cos(atan(tV));
   SV=sin(atan(tV));
  
  Zu_[0][0]=CU;
  Zu_[0][1]=SU;
  Zu_[1][0]=-SU;
  Zu_[1][1]=CU;
  Zv_[0][0]=CV;
  Zv_[0][1]=SV;
  Zv_[1][0]=-SV;
  Zv_[1][1]=CV;

  for(i=0;i<2;i++) mc[i]=g2v1*Zu_[i][0]*Zv_[i][1]
                        +g2v2*Zu_[i][1]*Zv_[i][0]
                        +  M2*Zu_[i][0]*Zv_[i][0]
                        +  mu*Zu_[i][1]*Zv_[i][1];
 
 for(i=1;i<=2;i++)
 { sprintf(name,"MC%d",i);
   assignValW(name,mc[i-1]);
   for(j=1;j<=2;j++)
   { 
      sprintf(name,"Zu%d%d",i,j);
      assignValW(name,Zu_[i-1][j-1]);
      sprintf(name,"Zv%d%d",i,j);
      assignValW(name,Zv_[i-1][j-1]);
   }    
 }

}


#define V(N) findValW(#N)

int nmssmEWSB(void)
{  int err;

   err=ewsbNMSSM(V(tb),V(MG1),V(MG2),V(MG3),V(Ml2),V(Ml3),V(Mr2),V(Mr3),
     V(Mq2),V(Mq3),V(Mu2),V(Mu3),V(Md2),V(Md3),V(At),V(Ab),V(Al),V(mu),
     V(LambdQ),V(KappaQ),V(aLmbdQ),V(aKappQ));

printf("ewsbNMSSM finishes\n");

   if(delFiles) system("rm -f  inp.dat spectr.dat omega.dat decay.dat out.dat");
   if(err) return err;
   FillVal(0);
   CharginoZM();   
   return err; 
}

#undef V

int nmssmSUGRA(double  m0,double mhf, double a0,double tb, double sgn,
double  Lambda, double aLambda, double aKappa)
{  int err;

   err= sugraNMSSM(m0, mhf, a0, tb, sgn, Lambda, aLambda,aKappa);
   if(delFiles) system("rm -f  inp.dat spectr.dat omega.dat decay.dat out.dat");
   if(err==0)
   { 
     FillVal(1);
     CharginoZM();
   }
   return err;
}

int readSLHA(char * fname)
{  double maxl;
   int err=slhaRead(fname, 0);   
   if(err) return err;

printf("readSLHA calls higgspotent\n");
   
   maxl=higgspotent(0);

   if(maxl<0) return -1;
   if(maxl>=1)  printf("Warning! max of Higgs coupling restored by micromegas"
                       " is %.1E\n",maxl);
   
   FillVal(2);
   CharginoZM();
   return 0;
}


double bsgnlo_(double *M, double *P)
{ *M= slhaVal("LOWEN",0.,1,12);
  *P= slhaVal("LOWEN",0.,1,11);
  return slhaVal("LOWEN",0.,1,1);
}

double deltamd_(double *M, double *P)
{ *M= slhaVal("LOWEN",0.,1,22);
  *P= slhaVal("LOWEN",0.,1,21);
  return slhaVal("LOWEN",0.,1,2);
}

double deltams_(double *M, double *P)
{ *M= slhaVal("LOWEN",0.,1,32);
  *P= slhaVal("LOWEN",0.,1,31);
  return slhaVal("LOWEN",0.,1,3);
}

double bsmumu_(double *M, double *P)
{ *M= slhaVal("LOWEN",0.,1,42);
  *P= slhaVal("LOWEN",0.,1,41);
  return slhaVal("LOWEN",0.,1,4);
}

double btaunu_(double *M, double *P)
{ *M= slhaVal("LOWEN",0.,1,52);
  *P= slhaVal("LOWEN",0.,1,51);
  return slhaVal("LOWEN",0.,1,5);
}

double gmuon_(double *M, double *P)
{ *M= slhaVal("LOWEN",0.,1,62);
  *P= slhaVal("LOWEN",0.,1,61);
  return slhaVal("LOWEN",0.,1,6);
}
    
