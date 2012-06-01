/*
 Copyright (C) 1997, Alexander Pukhov 
*/

#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>

#include "files.h"
#include "interface.h"
#include "subproc.h"
#include "chep_crt.h"
#include "strfun.h"
#include "n_calchep_.h"
#include "pdf.h"
#include "alphas2.h"
#include "sf_pdf.h"


/*static char * pdfName[2]={NULL,NULL}; */

static int nGroup[2]={3,3}, nSet[2]={41,41}, sgn[2]={1,1};
static int parton[15]={83,81,2,1,4,3,5,21,-5,-3,-4,-1,-2,-81,-83};
static int pnum[2]={0,0};


static char *param=NULL;

static double sc2=0.22*0.22;

static double alpha_pdf(double q ) { return alphas2_(&q); }

int mc_pdf(int i) { return 2212*sgn[i-1];}

static int initPDF(void)
{ int i;
  double v[20];
  
  static int dummy=0;

  if(dummy) return 0;
  if(param) return 1;

  param=malloc(400);
  for(i=0;i<400;i++) param[i]=' ';
  memcpy(param,"Init0",strlen("Init0"));
  v[0]=0;
  pdfset_(param,v,20);

  for(i=0;i<NGRMAX;i++) if(w505120_.NPGSMX[i][0]) break;
  if(i==NGRMAX) 
  { free(param);
    param=NULL;
    dummy=1;
    return 0;
  }
  
  if(w505120_.NPGSMX[2][0]<41) 
  { nGroup[0]=i+1,nGroup[1]=i+1;
    nSet[0]=1,    nSet[1]=1;
  }  
    
  memcpy(param,   "NPTYPE",strlen("NPTYPE")); 
  memcpy(param+20, "NGROUP",strlen("NGROUP"));
  memcpy(param+40,"NSET",strlen("NSET"));     
  return 1;
}

int p_pdf(long pNum) 
{  
  int i;
  if(!initPDF())return 0; 
  for(i=0;i<15;i++) if(pNum==parton[i]) return 1;
  return 0;
}

void n_pdf(int i, char *name) 
{  
  char txt[10];
  i--;
  memcpy(txt,w505110_.SFNAME[nSet[i]-1][nGroup[i]-1][0],8); 
  txt[8]=0; trim(txt); 
  sprintf(name,"PDFLIB: %8.8s Sign=%d",txt,sgn[i]);
}

int r_pdf(int i, char *name)
{ int k,l; 
  char txt[10],txt2[10];
  i--;

  if(!initPDF())return 0; 
  
  if(2!=sscanf(name,"PDFLIB:%*c%8c Sign=%d",txt,sgn+i)) return 0;
  txt[8]=0;
  trim(txt);

  for(k=0;k<NGRMAX;k++)
  for(l=w505120_.NPGSMX[k][0]-1;l+1;l--)
  {  
     memcpy(txt2,w505110_.SFNAME[l][k][0],8);
     txt2[8]=0;
     trim(txt2);
     if(strcmp(txt,txt2)==0)
     { nGroup[i]=k+1;
       nSet[i]=l+1;
       return 1;
     }
  }
  return 0;
}


int be_pdf(int i,double * be, double * mass) 
{  int k;
   long N,N1,N2;
   pinf_int( Nsub,1,NULL,&N1);
   pinf_int( Nsub,2,NULL,&N2);

   if(i==1) N=N1; else N=N2;
   if(abs(N1)>80 &&  abs(N2)>80) {if(N>0) N-=80; else N+=80;}
    
   *mass=0.9383;
   *be=1;
   sf_alpha=&(alpha_pdf);   
   for(k=0;k<15;k++) if(parton[k]==N) {pnum[i-1]=N; return 1;}
   pnum[i-1]=0; 
   return 0;
}
 

int m_pdf(int i)
{ 
  int size=2+8*500;
  int size_=1;
  void *pscr=NULL;
  void *pscr0=NULL;

  char * strmen=malloc(size);

  int n1=1,n0=1,k,l;
  i--;
  initPDF();

  strcpy(strmen,"\10");

  for(k=0;k<NGRMAX;k++)
  for(l=w505120_.NPGSMX[k][0]-1;l+1;l--)
  { 
     if(size<size_+10) { size+=8*50; strmen=realloc(strmen,size);}
     memcpy(strmen+size_,w505110_.SFNAME[l][k][0],8);
     size_+=8;
     if(k+1<nGroup[i] ||(k+1==nGroup[i] &&  l+1>nSet[i])) n1++;
  }

  strmen[size_]=0;

  for(;n0!=0 && n0!=3;)
  {  
     char strmen0[]="\020"
                    " SET =          "
                    " Charge  = NNN  "
                    " OK             ";
                    
     memcpy(strmen0+8,w505110_.SFNAME[nSet[i]-1][nGroup[i]-1][0],8);
     improveStr(strmen0,"NNN","%d",sgn[i]);
     menu1(54,14,"",strmen0,"",&pscr0,&n0);
     switch(n0) 
     { case 0: free(strmen);return 0;
       case 1: 
       { int n=n1;
          
          menu1(5,5,"",strmen,"",&pscr,&n);
          if(n)
          { n1=n;
            for(k=0;k<NGRMAX;k++)
            { for(l=w505120_.NPGSMX[k][0];l;l--){n--; if(!n) break;}
              if(l) break;
            } 
             nSet[i]=l;
             nGroup[i]=k+1;
             put_text(&pscr);
          } 
       } break;
       case 2: sgn[i]=-sgn[i]; break;
       case 3: put_text(&pscr0); break;
     }
  }  
  free(strmen);
  return 1;
}

double c_pdf(int i, double x, double q)
{
  static int nGroup_=-1, nSet_=-1;
  double upv, dnv,usea, dsea, str, chm, bot, top, gl;
  int p;
  i--;
  p=pnum[i];

  if(nGroup_!=nGroup[i] || nSet_ != nSet[i])
  { 
     double v[20];
     nGroup_=nGroup[i]; nSet_=nSet[i];
     v[0]=1;
     v[1]=nGroup_;
     v[2]=nSet_;
     pdfset_(param,v,20);
  }

  structm_(&x,&q,&upv,&dnv,&usea,&dsea,&str,&chm,&bot,&top,&gl); 
  if(sgn[i]<0) p=-p;

  switch(p)
  { case 81:           return ((dnv+dsea)*(1-sc2) +str*sc2)/x;
    case 83:           return ((dnv+dsea)*sc2 +str*(1-sc2))/x;
    case 2 :           return (upv+usea)/x;
    case 1 :           return (dnv+dsea)/x;
    case 3 : case -3 : return str/x;
    case 4 : case -4 : return chm/x;
    case 5 : case -5 : return bot/x; 
    case 21: case -21: return gl/x;
    case -1:           return dsea/x;
    case -2:           return usea/x;
    case -81:          return (dsea*(1-sc2) +str*sc2)/x;
    case -83:          return (dsea*sc2 +str*(1-sc2))/x;
  } 
}

