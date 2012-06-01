/*
 Copyright (C) 2002, Alexander Pukhov
*/

#include <unistd.h>
#include <stdarg.h>
#include <math.h>

#include "interface.h"
#include "cut.h"
#include "4_vector.h"
#include "q_kin.h"
#include "regul.h"
#include "runVegas.h"
#include "rw_sess.h"
#include "subproc.h"
#include "vegas.h"
#include "strfun.h"
#include "alphas2.h"
#include "crt_util.h"
#include "histogram.h"
#include "drandXX.h"
#include "n_calchep_.h"
#include "files.h"
#include "events.h"

int nSess=1;

double inP1=7000, inP2=7000;

static vegasGrid * veg_Ptr=NULL;
static int hFill=0;
vegas_integral integral={{5,0},{10000,10000},0,0.,0.,0.,0.,0.,0.,0,0,0}; 


static void clearStatistics(void)
{
  integral.I=0;
  integral.dI=0;
  integral.khi2=0;
  integral.s0=0; 
  integral.s1=0; 
  integral.s2=0; 
  integral.n_it=0; 
  integral.nCallTot=0; 

  clearHists();  
  { char fname[20];
    sprintf(fname,"distr_%d",nSess);
    unlink(fname);
  }
}                         

void clearGrid(void){ vegas_finish(veg_Ptr); veg_Ptr=NULL;}
void clearEventMax(void)
 { if(veg_Ptr && veg_Ptr->fMax) {free(veg_Ptr->fMax); veg_Ptr->fMax=NULL;}}

void newSession(void)
{
   if(integral.old)
   { char fname[20];

     messanykey(15,15,
     "Some parameters where changed.\nSo integral and statictics for\n"
     "distribushions is forgotten!\nSession number is increased.");
     integral.old=0;
     nSess++;
     clearStatistics();
     sprintf(fname,"prt_%d",nSess);
     unlink(fname);
   }
}



int saveVegasGrid( FILE * f)
{
  if(veg_Ptr)
  {  int i,j;
     double * x=veg_Ptr->x_grid;
     fprintf(f," Vegas_grid: dim=%d  size=%d\n", veg_Ptr->ndim, veg_Ptr->ndmx);
     for(i=0;i<veg_Ptr->ndim;i++)
     { for(j=0;j<=veg_Ptr->ndmx;j++) fprintf(f," %.15E",*(x++));
       fprintf(f,"\n");
     }
     if(veg_Ptr->fMax)
     { long l;
       fprintf(f,"Max(%d):\n",veg_Ptr->nCubs);
       for(l=0;l<veg_Ptr->nCubs;l++) fprintf(f,"%.1E\n",veg_Ptr->fMax[l]);
     } else fprintf(f,"Max(0):\n");
  }else  fprintf(f," Vegas_grid: dim=%d  size=%d\n", 0, 0);
  return 0;
}

int readVegasGrid(FILE * f)
{
  int i,j,ndim,ndmx;
  double * x;
  
  if(veg_Ptr) {vegas_finish(veg_Ptr);veg_Ptr=NULL;}  
  fscanf(f," Vegas_grid: dim=%d  size=%d\n", &ndim, &ndmx);
  if(ndim && ndmx)
  { 
    veg_Ptr=vegas_init(ndim,ndmx);
    x=veg_Ptr->x_grid;
    for(i=0;i<ndim;i++)for(j=0;j<=ndmx;j++) fscanf(f," %lf",(x++));
    fscanf(f," Max(%ld):\n",&(veg_Ptr->nCubs));       
    if(veg_Ptr->nCubs) 
    { long l;
      veg_Ptr->fMax=malloc(sizeof(float)*veg_Ptr->nCubs);
      for(l=0;l<veg_Ptr->nCubs;l++) fscanf(f,"%f",veg_Ptr->fMax+l);
    } else veg_Ptr->fMax=NULL;
  }
  return 0;
}

static int nCall;
static double badPoints;
static double negPoints;

static void printLn(FILE * iprt,int *line,char * format, ...)
{  
   va_list args;
   char dump[STRSIZ];
   va_start(args, format);
   vsprintf(dump,format,args);
   va_end(args);

   goto_xy(1,*line); print("%53s","");
   goto_xy(1,*line);
   print("%s\n",dump); 
   (*line)++;
   if (*line >= maxRow()-2 ) *line=8; else 
   {
     scrcolor(Blue, BGmain);
     print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
   }   
   if(iprt) {fprintf(iprt,"%s\n",dump); fflush(iprt);}
}

static double func_(double *x, double wgt)
{
    double ret_val=0.;
    int err=0;
    double factor_0;
    double q;
    double xx0=x[0],xx1=x[1],x1,x2;    
    nCall++; 
/* ** call kinematics preparation of scalar products */
    mkmom(x, &factor_0);

    if(sf_num[0]){x1=x[0]; x[0]=xx0;}
    if(sf_num[1]){x2=x[1]; x[1]=xx1;}
    if (!factor_0) goto exi;
    
    factor_0 *= calcCutFactor(); 
    if (!factor_0)   goto exi;

    q=Scale();
/* **  structure function  multiplication */
    if (nin_int == 2) 
    {
	if(sf_num[0]) { factor_0 *= strfun_(1, x1,q);  if(factor_0==0.) {/*printf("|x1=%.2f|",x1);*/ goto exi;}}
	if(sf_num[1]) { factor_0 *= strfun_(2, x2,q);  if(factor_0==0.) {/*printf("|x2=%.2f|",x2);*/ goto exi;}} 
    }   
    if (!factor_0)  { printf("strf");  goto exi;}
/* ** call for 'running strong coupling constant' */
    alf_(q);
    ret_val = factor_0 * sqme_int(Nsub,pvect,&err);
    if(err)       badPoints+=  (ret_val>0 ? ret_val*wgt : - ret_val*wgt); 
    if(ret_val<0) negPoints+=ret_val*wgt;
exi:

    if(hFill) fillHists(ret_val*wgt);

    return ret_val;
} /* func_ */


int runVegas(void)
{
    int i;
    double sd;
    double avgi;
    char mess[25];
    FILE * iprt = NULL;
    int mode=1;
    void * pscr=NULL;
    static int n_Line=7;

    i=imkmom(inP1,inP2);
    if(veg_Ptr&&veg_Ptr->ndim!=i)clearGrid();
    if(!veg_Ptr) veg_Ptr=vegas_init(i,50);     

    if(nin_int == 2) strcpy(mess, " Cross section [pb] ");
      else           strcpy(mess, "  Width     [Gev]   ");
    
/* ** save current session parameters */
     w_sess__(NULL);
/* ** open protocol and resulting files */
       
    {  char fname[50];
       sprintf(fname,"%sprt_%d",outputDir,nSess);
       iprt=fopen(fname,"a");
       if(ftell(iprt)==0) 
       { fprintf(iprt,"    CalcHEP kinematics module \n The session parameters:\n");
         w_sess__(iprt);
         fprintf(iprt,"===================================\n");   
       }
    }

/* **  initkinematics */

    correctHistList();
    
/* *** Main cycle */
    if(!integral.old || n_Line==7)
    { n_Line=7;
      scrcolor(Blue, BGmain);
      printLn(iprt,&n_Line," #IT  %20s Error %%    nCall   chi**2",mess); 
    }

    for(;;)
    {   int k;
        char strmen[]="\030"
         " nSess_1  = N2_1        "
         " nCalls_1 = N1_1        "
         " nSess_2  = N2_2        "
         " nCalls_2 = N1_2        "                            
         " Set  Distributions     "
         "*Start integration      "
         " Display Distributions  "
         " Clear statistic        "
         " Freeze grid        OFF " 
	 " Clear  grid            ";

        improveStr(strmen,"N1_1","%d",integral.ncall[0]);
        improveStr(strmen,"N2_1","%d",integral.itmx[0]);
        improveStr(strmen,"N1_2","%d",integral.ncall[1]);
        improveStr(strmen,"N2_2","%d",integral.itmx[1]);
        
        if(integral.freeze) improveStr(strmen,"OFF","ON");

        menu1(54,7,"",strmen,"n_veg_*",&pscr,&mode);
        switch(mode)
        {     
        case 0:
          if(iprt) fclose(iprt);
          return 0;           
        case 1:  
          correctInt(50,12,"Enter new value ",&integral.itmx[0],1); break;
        case 2: 
          correctLong(50,12,"Enter new value ",&integral.ncall[0],1); break;
        case 3:  
          correctInt(50,12,"Enter new value ",&integral.itmx[1],1); break;
        case 4: 
          correctLong(50,12,"Enter new value ",&integral.ncall[1],1); break;

        case 5:  editHist(); break;
        case 6:
          if(veg_Ptr->fMax && !integral.freeze)
          {  if(!mess_y_n(15,15,"You have event generator prepared.\n"
             " The  answer 'Y'  will start Vegas session \nwhich destroys it."
             " To save the event generator answer 'N' \nand set "
             " ' Freeze grid' ON")) break;
             else { free(veg_Ptr->fMax); veg_Ptr->fMax=NULL; }  
          }
          if(!blind && integral.old && integral.ncall[1]*integral.itmx[1] &&
           !mess_y_n(15,15," The obtained results will be"
             " lost because nSess_2!=0.\n Continue?")) break;     
          for(k=0;k<2;k++)for (i = 1; i <= integral.itmx[k]; ++i)                                       
          { char  errtxt[100]="";                                                     
            if(integral.ncall[k]==0) break;
            if(k==1 && i==1)  
            { clearStatistics();
              sprintf(errtxt,"clear statistics.");
            }                                                                           
            nCall=0;                                                                  
            negPoints=0;                                                              
            badPoints=0; 
            hFill=1;  
            if(vegas_int(veg_Ptr, integral.ncall[k],1.5*(!integral.freeze), 
                 func_, &avgi, &sd)        
              ) break;
            integral.old=1;                                              
            negPoints/=nCall;                                                         
            badPoints/=nCall;                                                         
            integral.nCallTot+=nCall;                                                          
            scrcolor(FGmain,BGmain);                                                  
            printLn(iprt,&n_Line,"%3d    %12.4E     %10.2E %8d ",                     
                 ++integral.n_it, avgi,100*sd/fabs(avgi),nCall);                                 
            if(negPoints<0) sprintf(errtxt+strlen(errtxt)," Negative points %.1G%%;",                
                                      -100*negPoints/(avgi-2*negPoints));             
            if(badPoints)  sprintf(errtxt+strlen(errtxt),                             
                 "Bad Precision %.1G%%;",100*badPoints/(avgi-2*negPoints));           
                                                                                      
            if(errtxt[0])                                                             
            {                                                                         
               scrcolor(Red,BGmain);                                                  
               printLn(iprt,&n_Line,"%s",errtxt);                                     
            }                                                                         
            sd=1/(sd*sd);                                                             
            integral.s0+=sd;                                                                  
            integral.s1+=avgi*sd;                                                             
            integral.s2+=avgi*avgi*sd;                                                        
          } 
          
          integral.I=integral.s1/integral.s0; 
          integral.dI=1/sqrt(integral.s0);
          if(integral.n_it<=1) integral.khi2=0; else 
          integral.khi2=(integral.s2-integral.s1*integral.s1/integral.s0)/(integral.n_it-1); 

          scrcolor(FGmain,BGmain);                                                 
          printLn(iprt,&n_Line," < >   %12.4E     %10.2E %8d %8.1G" ,         
          integral.I, 100*integral.dI/fabs(integral.I), integral.nCallTot, integral.khi2);

          if(histTab.strings)
          { char  fname[20];
            FILE * d;
            sprintf(fname,"distr_%d",nSess);
            d=fopen(fname,"w");  
            wrt_hist2(d,Process);
            fclose(d);
          }
          messanykey(54,11,"Integration is over");
/*          integral.freeze=0; */
          break;

        case 7: showHist(54,10); break;
        case 8: clearStatistics();
                messanykey(54,13,"Old results for integral\n"
                "and distributions\nare deleted.");
                break;
        case 9: integral.freeze=!integral.freeze; break; 
        case 10: { int ndim=veg_Ptr->ndim;
                  vegas_finish(veg_Ptr);
                  veg_Ptr=vegas_init(ndim,50);
                }   
                messanykey(57,11,"OK");
                break;
       }
    }    
}


int runEvents(void)
{
    FILE * iprt = NULL;
    int i;

    i=imkmom(inP1,inP2);
    if(veg_Ptr&&veg_Ptr->ndim!=i)clearGrid();
    if(!veg_Ptr) veg_Ptr=vegas_init(i,50);    

    w_sess__(NULL);
/* ** open protocol and resulting files */
       
    {  char fname[50];
       sprintf(fname,"%sprt_%d",outputDir,nSess);
       iprt=fopen(fname,"a");
       if(ftell(iprt)==0) 
       { fprintf(iprt,"    CalcHEP kinematics module \n The session parameters:\n");
         w_sess__(iprt);
         fprintf(iprt,"===================================\n");   
       }
    }

/* **  initkinematics */


    { char fname[50];

      sprintf(fname,"%sevents_%d.txt",outputDir,nSess);

      hFill=0;
      generateEvents(veg_Ptr,func_,fname, iprt);
    }
    fclose(iprt);
    return 0;
}
