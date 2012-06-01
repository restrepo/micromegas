#include<stdio.h>
#include<math.h>
#include <stdlib.h>
#include<string.h>
/*#include"../../CalcHEP_src/c_source/model_aux/include/alpha_s.h" */
#include"../../sources/micromegas.h"
#include "pmodel.h"
#include "localpath.h"

extern int FError;


#define  NVAR 126

typedef struct{char * name; double val;} varSTR;


static varSTR modelVar[NVAR]=
{
{"la1",   0 },{"la2",   0 },{"la3",   0 },{"la4",   0 },{"la5r",  0 },{"la5i",  0 },
{"la6r",  0 },{"la6i",  0 },{"la7r",  0 },{"la7i",  0 },{"MHc",   0 },{"Mh1",   0 },
{"Mh2",   0 },{"Mh3",   0 },{"Zh11",  0 },{"Zh12",  0 },{"Zh13",  0 },{"Zh21",  0 },
{"Zh22",  0 },{"Zh23",  0 },{"Zh31",  0 },{"Zh32",  0 },{"Zh33",  0 },{"MSt1",  0 },
{"MSt2",  0 },{"Zt11r", 0 },{"Zt12r", 0 },{"Zt21r", 0 },{"Zt22r", 0 },{"Zt11i", 0 },
{"Zt12i", 0 },{"Zt21i", 0 },{"Zt22i", 0 },{"MSb1",  0 },{"MSb2",  0 },{"Zb11r", 0 },
{"Zb12r", 0 },{"Zb21r", 0 },{"Zb22r", 0 },{"Zb11i", 0 },{"Zb12i", 0 },{"Zb21i", 0 },
{"Zb22i", 0 },{"MSnl",  0 },{"MSl1",  0 },{"MSl2",  0 },{"Zl11r", 0 },{"Zl12r", 0 },
{"Zl21r", 0 },{"Zl22r", 0 },{"Zl11i", 0 },{"Zl12i", 0 },{"Zl21i", 0 },{"Zl22i", 0 },
{"MNE1",  0 },{"MNE2",  0 },{"MNE3",  0 },{"MNE4",  0 },{"Zn11r", 0 },{"Zn11i", 0 },
{"Zn12r", 0 },{"Zn12i", 0 },{"Zn13r", 0 },{"Zn13i", 0 },{"Zn14r", 0 },{"Zn14i", 0 },
{"Zn21r", 0 },{"Zn21i", 0 },{"Zn22r", 0 },{"Zn22i", 0 },{"Zn23r", 0 },{"Zn23i", 0 },
{"Zn24r", 0 },{"Zn24i", 0 },{"Zn31r", 0 },{"Zn31i", 0 },{"Zn32r", 0 },{"Zn32i", 0 },
{"Zn33r", 0 },{"Zn33i", 0 },{"Zn34r", 0 },{"Zn34i", 0 },{"Zn41r", 0 },{"Zn41i", 0 },
{"Zn42r", 0 },{"Zn42i", 0 },{"Zn43r", 0 },{"Zn43i", 0 },{"Zn44r", 0 },{"Zn44i", 0 },
{"MC1",   0 },{"MC2",   0 },{"Zu11r", 0 },{"Zu11i", 0 },{"Zu12r", 0 },{"Zu12i", 0 },
{"Zu21r", 0 },{"Zu21i", 0 },{"Zu22r", 0 },{"Zu22i", 0 },{"Zv11r", 0 },{"Zv11i", 0 },
{"Zv12r", 0 },{"Zv12i", 0 },{"Zv21r", 0 },{"Zv21i", 0 },{"Zv22r", 0 },{"Zv22i", 0 },
{"Tu3r",  0 },{"Tu3i",  0 },{"Td3r",  0 },{"Td3i",  0 },{"Bsmumu",0 },{"Bsg",   0 },
{"Bulnu", 0 },{"Bdll",  0 },{"ABsg",  0 },{"deltMd",0 },{"deltMs",0 },{"EDMth", 0 },
{"EDMHg", 0 },{"EDMel", 0 },{"EDMmu", 0 },{"LEPex", 0 },{"LEPrat",0 },{"gmuon", 0 }
};


/*
void printVar(FILE *f, int N)
{
  int i;
  fprintf(f,"\n# MSSM parameters:\n");
  if(N>NVAR) N=NVAR;
  for(i=0;i<N;i++) fprintf(f,"%-6.6s   %f\n", modelVar[i].name, modelVar[i].val);
}
*/

static double * varAddress_cp(char *name)
{int i;
 for(i=0;i<NVAR;i++)if(!strcmp(name,modelVar[i].name))return &(modelVar[i].val);
 return NULL;
}

static int assignVal_cp(char * name, double val)
{
  double * a=varAddress_cp(name);
  if(a)  {*a=val; return 0;} else return 1;
}


static int findVal_cp(char * name, double * val)
{
  double * va=varAddress_cp(name); 
  if(va){*val=*va; return 0;} else  return 1; 
}



static double  findValW_cp(char*name)
{
  double val;

  if(findVal_cp(name,&val))
  {  printf(" name %s not found\n",name);
     return 0;
  } else return val;
}




static int readVar_cp(char *fname)
{
  double val;
  char name[80];
  int n;
  FILE * f=fopen(fname,"r");
  if(f==NULL) return -1;
  
  for(n=1;;n++)
  { if(fscanf(f,"%s",name)!=1) { if(n<NVAR) n=-2; else n=0;  break;}
    if(name[0]=='#') { fscanf(f,"%*[^\n]"); continue;}
    if(fscanf(f,"%lf",&val)!=1) break;
    fscanf(f,"%*[^\n]");
    { int err=assignVal_cp(name,val);
      if(err)printf("unknown name in micromegas.in '%s'\n",name);
    }
  }
  fclose(f);
  return n;                                         
}

#define NTOFc(X) double X(double ok){return findValW_cp(#X);}  

NTOFc(la1)   NTOFc(la2)   NTOFc(la3)   NTOFc(la4)   NTOFc(la5r)  NTOFc(la5i)  
NTOFc(la6r)  NTOFc(la6i)  NTOFc(la7r)  NTOFc(la7i)  NTOFc(MHc)   NTOFc(Mh1)   
NTOFc(Mh2)   NTOFc(Mh3)   NTOFc(Zh11)  NTOFc(Zh12)  NTOFc(Zh13)  NTOFc(Zh21)  
NTOFc(Zh22)  NTOFc(Zh23)  NTOFc(Zh31)  NTOFc(Zh32)  NTOFc(Zh33)  NTOFc(MSt1)  
NTOFc(MSt2)  NTOFc(Zt11r) NTOFc(Zt12r) NTOFc(Zt21r) NTOFc(Zt22r) NTOFc(Zt11i) 
NTOFc(Zt12i) NTOFc(Zt21i) NTOFc(Zt22i) NTOFc(MSb1)  NTOFc(MSb2)  NTOFc(Zb11r) 
NTOFc(Zb12r) NTOFc(Zb21r) NTOFc(Zb22r) NTOFc(Zb11i) NTOFc(Zb12i) NTOFc(Zb21i) 
NTOFc(Zb22i) NTOFc(MSnl)  NTOFc(MSl1)  NTOFc(MSl2)  NTOFc(Zl11r) NTOFc(Zl12r) 
NTOFc(Zl21r) NTOFc(Zl22r) NTOFc(Zl11i) NTOFc(Zl12i) NTOFc(Zl21i) NTOFc(Zl22i) 
NTOFc(MNE1)  NTOFc(MNE2)  NTOFc(MNE3)  NTOFc(MNE4)  NTOFc(Zn11r) NTOFc(Zn11i) 
NTOFc(Zn12r) NTOFc(Zn12i) NTOFc(Zn13r) NTOFc(Zn13i) NTOFc(Zn14r) NTOFc(Zn14i) 
NTOFc(Zn21r) NTOFc(Zn21i) NTOFc(Zn22r) NTOFc(Zn22i) NTOFc(Zn23r) NTOFc(Zn23i) 
NTOFc(Zn24r) NTOFc(Zn24i) NTOFc(Zn31r) NTOFc(Zn31i) NTOFc(Zn32r) NTOFc(Zn32i) 
NTOFc(Zn33r) NTOFc(Zn33i) NTOFc(Zn34r) NTOFc(Zn34i) NTOFc(Zn41r) NTOFc(Zn41i)
NTOFc(Zn42r) NTOFc(Zn42i) NTOFc(Zn43r) NTOFc(Zn43i) NTOFc(Zn44r) NTOFc(Zn44i) 
NTOFc(MC1)   NTOFc(MC2)   NTOFc(Zu11r) NTOFc(Zu11i) NTOFc(Zu12r) NTOFc(Zu12i) 
NTOFc(Zu21r) NTOFc(Zu21i) NTOFc(Zu22r) NTOFc(Zu22i) NTOFc(Zv11r) NTOFc(Zv11i) 
NTOFc(Zv12r) NTOFc(Zv12i) NTOFc(Zv21r) NTOFc(Zv21i) NTOFc(Zv22r) NTOFc(Zv22i) 
NTOFc(Tu3r)  NTOFc(Tu3i)  NTOFc(Td3r)  NTOFc(Td3i)
#undef NTOF


double cpHiggs(double EE,double alfaSMZ,double Mtp, double McMt,double MbMt,double Ml, double SW, double tb,  double MHc,  
               double aMu, double fiMu, double aM1, double fiM1, double aM2, double fiM2,
               double  aM3, double  fiM3, double Mq3, double Mu3,  double Md3,  double Ml3,  double  Mr3, 
               double aAt, double fiAt, double aAb,  double fiAb, double  aAl, double fiAl,
               double aAe, double fiAe, double Au, double Ad)
{ FILE *f;
  char * command;
  int err;  
 f=fopen("CPsuperH.in","w");  
 
 fprintf(f,"%E              ! SMPARA( 1) = 1/AEM(MZ)  \n",4*M_PI/EE/EE );
 fprintf(f,"%E             ! SMPARA( 2) = AS(MZ)  \n",alfaSMZ);
 fprintf(f,"91.187D0            ! SMPARA( 3) = MZ in GeV \n");
 fprintf(f,"%E             ! SMPARA( 4) = sin^2 Theta_W \n",SW*SW);
 fprintf(f,"0.5D-3              ! SMPARA( 5) = m_e in GeV \n");
 fprintf(f,"0.1065D0            ! SMPARA( 6) = m_mu in GeV \n");
 fprintf(f,"%E             ! SMPARA( 7) = m_tau in GeV \n",Ml);
 fprintf(f,"6.D-3               ! SMPARA( 8) = m_d (m_t) in GeV \n");
 fprintf(f,"0.115D0             ! SMPARA( 9) = m_s (m_t) in GeV \n");
 fprintf(f,"%E             ! SMPARA(10) = m_b (m_t) in GeV \n",MbMt);
 fprintf(f,"3.D-3               ! SMPARA(11) = m_u (m_t) in GeV \n");
 fprintf(f,"%E              ! SMPARA(12) = m_c (m_t) in GeV \n",McMt);
 fprintf(f,"%E             ! SMPARA(13) = m_t^POLE  in GeV \n",Mtp);
 fprintf(f,"2.118D0             ! SMPARA(14) = Gam_W  in GeV \n" );
 fprintf(f,"2.4952D0            ! SMPARA(15) = Gam_Z  in GeV \n");
/*new*/ 
 fprintf(f,"0.2272D0            ! SMPARA(16) = lambda_CKM\n");
 fprintf(f,"0.8180D0            ! SMPARA(17) = A_CKM\n");
 fprintf(f,"0.2210D0            ! SMPARA(18) = rho^bar_CKM\n");
 fprintf(f,"0.3400D0            ! SMPARA(19) = eta^bar_CKM\n");
/* end of new*/

 fprintf(f,"%E      ! SSPARA( 1) = tan beta \n",             tb      );
 fprintf(f,"%E      ! SSPARA( 2) = m_H^pm^POLE in GeV \n",   MHc     );
 fprintf(f,"%E      ! SSPARA( 3) = |mu| in GeV \n",          aMu     );
 fprintf(f,"%E      ! SSPARA( 4) = Phi_mu in Degree \n",     fiMu    );
 fprintf(f,"%E      ! SSPARA( 5) = |M_1| in GeV \n",         aM1     );
 fprintf(f,"%E      ! SSPARA( 6) = Phi_1 in Degree \n",      fiM1    );
 fprintf(f,"%E      ! SSPARA( 7) = |M_2| in GeV \n",         aM2     );
 fprintf(f,"%E      ! SSPARA( 8) = Phi_2 in Degree \n",      fiM2    );
 fprintf(f,"%E      ! SSPARA( 9) = |M_3| in GeV \n",         aM3     );
 fprintf(f,"%E      ! SSPARA(10) = Phi_3 in Degree \n",      fiM3    );
 fprintf(f,"%E      ! SSPARA(11) = m_Q3 in GeV \n",          Mq3     );
 fprintf(f,"%E      ! SSPARA(12) = m_U3 in GeV \n",          Mu3     );
 fprintf(f,"%E      ! SSPARA(13) = m_D3 in GeV \n",          Md3     );
 fprintf(f,"%E      ! SSPARA(14) = m_L3 in GeV \n",          Ml3     );
 fprintf(f,"%E      ! SSPARA(15) = m_E3 in GeV \n",          Mr3     );
 fprintf(f,"%E      ! SSPARA(16) = |A_t| in GeV \n",         aAt     );
 fprintf(f,"%E      ! SSPARA(17) = Phi_{A_t} in Degree \n",  fiAt    );
 fprintf(f,"%E      ! SSPARA(18) = |A_b| in GeV \n",         aAb     );
 fprintf(f,"%E      ! SSPARA(19) = Phi_{A_b} in Degree \n",  fiAb    );
 fprintf(f,"%E      ! SSPARA(20) = |A_tau| in GeV \n",       aAl     );
 fprintf(f,"%E      ! SSPARA(21) = Phi_{A_tau} in Degree \n",fiAl    );
/* new */
fprintf(f,"%E              ! SSPARA(22) = Hierarchy factor between first 2 and third generations M_Q\n",findValW("Mq2")/Mq3);
fprintf(f,"%E              ! SSPARA(23) = Hierarchy factor between first 2 and third generations M_U\n",findValW("Mu2")/Mu3);
fprintf(f,"%E              ! SSPARA(24) = Hierarchy factor between first 2 and third generations M_D\n",findValW("Md2")/Md3);
fprintf(f,"%E              ! SSPARA(25) = Hierarchy factor between first 2 and third generations M_L\n",findValW("Ml2")/Ml3);
fprintf(f,"%E              ! SSPARA(26) = Hierarchy factor between first 2 and third generations M_E\n",findValW("Mr2")/Mr3);
/*===*/

fprintf(f,"%E               ! SSPARA(27) = |A_e| in GeV\n",aAe);
fprintf(f,"%E               ! SSPARA(28) = Phi_{A_e} in Degree\n",fiAe);
fprintf(f,"%E               ! SSPARA(29) = |A_mu| in GeV\n",aAe);
fprintf(f,"%E               ! SSPARA(30) = Phi_{A_mu} in Degree\n",fiAe);
fprintf(f,"%E               ! SSPARA(31) = |A_u| in GeV\n",Au);
fprintf(f,"0.               ! SSPARA(32) = Phi_{A_u} in Degree\n");
fprintf(f,"%E               ! SSPARA(33) = |A_c| in GeV\n",Au);
fprintf(f,"0.               ! SSPARA(34) = Phi_{A_c} in Degree\n");
fprintf(f,"%E               ! SSPARA(35) = |A_d| in GeV\n",Ad);
fprintf(f,"0.               ! SSPARA(36) = Phi_{A_d} in Degree\n");
fprintf(f,"%E               ! SSPARA(37) = |A_s| in GeV\n",Ad);
fprintf(f,"0.               ! SSPARA(38) = Phi_{A_s} in Degree\n");


 fprintf(f,"0                   ! IFLAG_H(1)  if 1, print input parameters \n");
 fprintf(f,"1                   ! IFLAG_H(2)  if 1, print Higgs sector \n");
 fprintf(f,"1                   ! IFLAG_H(3)  if 1, print sfermion sector \n");
 fprintf(f,"1                   ! IFLAG_H(4)  if 1, print -ino sector \n");
 fprintf(f,"0                   ! IFLAG_H(5)  print couplings : 1=H_1 2=H_2 3=H_3 4=H^pm 5=Higgs self 6=ALL \n");
 fprintf(f,"5                   ! IFLAG_H(6)  print decay widths and brs : 1=H_1 2=H_2 3=H_3 4=H^pm 5=ALL \n");
 fprintf(f,"0                   ! IFLAG_H(10) if 0, include the threshold corrections  \n");
 fprintf(f,"0                   ! IFLAG_H(11) use pole mass (0) or eff. pot. mass (1) \n");
/* new */


 fprintf(f,"5                   ! IFLAG_H(12) 5 or 0 for full improvement\n");
 fprintf(f,"0                   ! IFLAG_H(13) 1 Not to include the off-diagonal absorptive parts\n");
 fprintf(f,"1                   ! IFLAG_H(14) 1 to print FILLDHPG results\n");
/* fprintf(f,"1                   ! IFLAG_H(15) 1 to print HiggsEDM results\n");*/
 fprintf(f,"1                   ! IFLAG_H(16) 1 to print FILLBOBS results\n");
 fprintf(f,"1                   ! IFLAG_H(17) 1 to print b -> s gamma details\n");
 fprintf(f,"1                   ! IFLAG_H(18) 1 to print EDM results\n");
 fprintf(f,"1                   ! IFLAG_H(19) 1 to print fllmuon results\n");

/*===*/
 fclose(f);

 command=malloc(strlen(lPath)+200);
 sprintf(command,"lPath=\"%s\"; HBpath=\"$lPath/HiggsBound/\"; export HBpath;  rm -f micromegas.in; $lPath/CPsuperH2/cpsuperh.exe < CPsuperH.in > CPsuperH.out; rm Key.dat",lPath);
 err=System(command);
 free(command);
 if(err>=0) err= readVar_cp("micromegas.in");
 
 if(delFiles) system("rm ./micromegas.in  ./CPsuperH.in ./CPsuperH.out");     
 return err;
}

static double   fiRe(double tau)
{
  double x;
  if(tau<1)
  { x=asin(sqrt(tau));
    return x*x;
  }else if(tau==1) return 0;
  else
  {
    x=sqrt(1-1/tau);
    x=log((1+x)/(1-x));
    return -0.25*(x*x - M_PI*M_PI);
  }   
}

static double fiIm(double tau)
{
  double x;
  if(tau<=1) return 0;
  else
  {
     x=sqrt(1-1/tau);
     x=log((1+x)/(1-x));
     return 0.5*x*M_PI;
  }
}

double FevenFr(double z)
{
  double tau= z*z/4;
  return (tau+(tau-1)*fiRe(tau))/(tau*tau);
}

double FevenFi(double z)
{
  double tau= z*z/4;
  return (tau-1)*fiIm(tau)/(tau*tau);
}

double Foddr(double z) { double tau=z*z/4; return fiRe(tau)/tau;}
double Foddi(double z) { double tau=z*z/4; return fiIm(tau)/tau;} 

double F0r(double z)   { double tau=z*z/4; return -(tau-fiRe(tau))/(tau*tau);}
double F0i(double z)   { double tau=z*z/4; return fiIm(tau)/(tau*tau);}

double mbPole(double ok){ return MbPole;}


 double bsmumu(void) { return findValW_cp("Bsmumu"); }
 double bsgnlo(void) { return findValW_cp("Bsg");    }
 double Bulnu(void)  { return findValW_cp("Bulnu"); }
 double Bdll(void)   { return findValW_cp("Bdll");  }
 double ABsg(void)   { return findValW_cp("ABsg");  }
 double deltaMd(void){ return findValW_cp("deltMd");}
 double deltaMs(void){ return findValW_cp("deltMs");}
 double EDMth(void)  { return findValW_cp("EDMth"); }
 double EDMel(void)  { return findValW_cp("EDMel"); }
 double EDMmu(void)  { return findValW_cp("EDMmu"); }
 double EDMHg(void)  { return findValW_cp("EDMHg"); }
 double gmuon(void)  { return findValW_cp("gmuon"); }
 int LepBound(double * ratio)
 {
    if(ratio) *ratio= findValW_cp("LEPrat");
    return findValW_cp("LEPex");   
 }
 

 double bsmumu_(void) { return findValW_cp("Bsmumu"); }
 double bsglno_(void) { return findValW_cp("Bsg");    }
 double bulnu_(void)  { return findValW_cp("Bulnu"); }
 double bdll_(void)   { return findValW_cp("Bdll");  }
 double absg_(void)   { return findValW_cp("ABsg");  }
 double deltamd_(void){ return findValW_cp("deltMd");}
 double deltams_(void){ return findValW_cp("deltMs");}
 double edmth_(void)  { return findValW_cp("EDMth"); }
 double edmel_(void)  { return findValW_cp("EDMel"); }
 double edmmu_(void)  { return findValW_cp("EDMmu"); }
 double edmhg_(void)  { return findValW_cp("EDMHg"); }
 double gmuon_(void)   { return findValW_cp("gmuon"); }

int lepbound_(double * ratio) { return LepBound(ratio);}
