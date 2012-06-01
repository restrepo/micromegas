/*
 Copyright (C) 1997, Alexander Pukhov, e-mail pukhov@theory.npi.msu.su 
*/

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <signal.h>

#include "chep_crt.h"
#include "files.h" 
#include "rw_sess.h"
#include "interface.h"
/*#include "num_in.h" */
#include "n_calchep_.h"
#include "viewdir.h"
#include"rootDir.h"    
#include"../../include/num_out.h"
#include "read_func.h"
#include "rd_num.h"
#include "parser.h"
                                                                    
static void f6_key_prog (int x){  viewDir("."); }
                                                                        

static void f10_key_prog (int x)
{
    if( mess_y_n(15,15," Quit session? ")) 
    {
        w_sess__(NULL);
        finish();
        sortie(0);
    }
}

static void f8_key_prog(int x)
{  
  static char FUNC[75]="2*2  % Press ESC to finish, F1 for help, ENTER to calculate";  
  int npos=1;
  void * pscr;
  get_text(1,20,80,21,&pscr);
  scrcolor(Red,White);   
  goto_xy(1,20); print("CALC : ");
  goto_xy(1,21); print("result=");
  scrcolor(Black,White);
  for(;;)
  { double res;
    int err; 
    int key;
    goto_xy(9,20); key=str_redact(FUNC,npos,70); 
    if(key==KB_ESC) break;
    else if(key==KB_F1)  show_help("n_calc");
    else if(key==KB_ENTER)
    {  goto_xy(8,21); 
       err=calcExpression(FUNC,rd_num,&res);
       goto_xy(9,21);
       if(err) {print("Erorr: %s",errmesstxt(err)); npos=rderrpos;}
       else  print("%E                     ",res); 
    }
  }
  put_text(&pscr);
}

static void f9_key_prog(int x)
{
  FILE*f;
  char fname[200];
  sprintf(fname,"%s%cCITE",pathtocomphep,f_slash);
  f=fopen(fname,"r");
  showtext (1, 1, 80,1,"",f);
  fclose(f);
}

static void xw_error(void) {sortie(80);}
  
int main(int argc,char** argv)
{
  int n;
  char icon[200];
/*
  char cwd[200];
  char version_[100];
*/  
  blind=0;
  for( n=1;n<argc;n++) if(strcmp(argv[n],"-blind")==0&& n<argc-1 )
  { blind=1;     
    inkeyString=argv[++n];
  } else  if (strcmp(argv[n],"+blind")==0 ) blind=2;
  

  if(!writeLockFile(".lock"))
  { fprintf(stderr,"locked by other n_calchep. See .lock\n");
    exit(100);
  }
                 
 
  sprintf(pathtocomphep,"%s%c",rootDir,f_slash);
  sprintf(pathtohelp,"%shelp%c",pathtocomphep,f_slash);
  sprintf(icon,"%sicon",pathtocomphep);

  f3_key[3]=f6_key_prog;   f3_mess[3]="Results";
  f3_key[5]=f8_key_prog;   f3_mess[5]="Calc";
  f3_key[6]=f9_key_prog;   f3_mess[6]="Ref";
  f3_key[7]=f10_key_prog;  f3_mess[7]="Quit";

/* **  initialization of the session */
  link_process(PtrInterface_ext);  

/*
  if(!getcwd(cwd,200)) strcpy(version_,version); else
  {
     for(n=strlen(cwd)-1; cwd[n]!=f_slash; n--);
     sprintf(version_,"%s:%s",version,cwd+n+1);
  }
*/  
  start1("CalcHEP/num",icon,"calchep.ini;../calchep.ini",&xw_error); 
  r_sess__(NULL); 

  goto_xy(10,10); print("Calculation of constaints.  Please, be patient.");
  escpressed();
  n_comphep();
  finish();  
  sortie(0);
}
