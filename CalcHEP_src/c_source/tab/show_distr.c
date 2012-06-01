#include <dirent.h>
#include "histogram.h"
#include "chep_crt.h"
#include "files.h" 
#include "rootDir.h"    
#include "subproc.h"
#include "viewdir.h"
#include "interface.h"
#include "V_and_P.h"

/*
#define const
#include"num_out.h"
#undef const

int  nvar_ext=0;
char * varName_ext[1];
double  va_ext[1];
int nfunc_ext=0; 

int nin_ext, nout_ext;
*/
static char* pinf_ext(int nsub,int num , double * mass, long *N ){ return NULL;}

int nModelVars=0;
int nModelFunc=0;
double aWidth(char * txt){ return 0;}

int nModelParticles=0;

ModelPrtclsStr ModelPrtcls[1];


double pWidth(char*pName,void * LL,int *dim){ return 0;}

char*varNames[1];
double varValues[1];
int calcMainFunc(void){return 0;}

static void f6_key_prog (int x){  viewDir("."); }

static void f10_key_prog (int x)
{ if( mess_y_n(15,15," Quit session? ")) { finish(); exit(0);} }    
    

 
int main(int argc, char ** argv)
{ 
   FILE*f;
   int nf=1;

   f3_key[3]=f6_key_prog;   f3_mess[3]="Results";
   f3_key[7]=f10_key_prog;  f3_mess[7]="Quit";
   blind=0;

   if(argc<2) 
   { printf(" This  routine is intended to display distributions produced\n" 
            " in  calcHEP numerical sessions (files dist_#).\n"
            " The name of file should submited to this routine as parameter\n");
     return 1;
   }

   if(strcmp(argv[1],"-blind")==0)
   { blind=1;     
     nf=3; 
     inkeyString=argv[2];
   } else  if (strcmp(argv[1],"+blind")==0 )
   { blind=2;
     nf=2;
   }
   nvar_int=0;
   nfunc_int=0; 
   pinf_int=&pinf_ext;
  if(nf!=argc-1) 
  { printf(" This  routine is intended to display distributions produced\n" 
           " in  calcHEP numerical sessions (files dist_#).\n"
           " The name of file should submited to this routine as parameter\n");
    return 1;
  }
  { char *p=getenv("CALCHEP");
    if(!p) p=rootDir;  
    sprintf(pathtocomphep,"%s%c",p,f_slash);
    sprintf(pathtohelp,"%s%chelp%c",p,f_slash,f_slash);
  }                          
              
  f=fopen(argv[nf],"r");
   
  if(!f)
  {  printf("Can not open the '%s' file\n",argv[nf]);
     exit(1);
  }   

  if(rdr_hist2(f,Process))
  {  printf("Wrong  format for %s \n",argv[nf]);
         exit(1);
  }else     
  {  char*c1,*c2;
     c1=strstr(Process,",");
     c2=strstr(Process,"->");
     if(c1>c2) nin_int=1; else nin_int=2;
     nout_int=2;
  }
  fclose(f);

  { char txt[80];
    sprintf(txt," %s Viewer for  distributions ",version);   
    start1(txt,pathtocomphep,"calchep.ini;../calchep.ini",NULL);
  }    
    
  goto_xy(20,2);  
  scrcolor(Blue,LightGray);     print(" File:"); 
  scrcolor(Red,LightGray);      print(" %s     ",argv[nf]);     
  scrcolor(Blue,BGmain);        print(" Process: ");
  scrcolor(Red,LightGray);      print(" %s",Process);

  showHist(54,5);
  
  finish();

  return 0;
}
