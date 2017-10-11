
#include"micromegas.h"
#include"micromegas_aux.h"
#include"../CalcHEP_src/c_source/chep_crt/include/crt.h"

#include <sys/types.h>
#include <unistd.h>
#include <signal.h>              
#include <sys/wait.h>
#include <stdarg.h>
 

 
static int First=1;
  
static void disconnect(int N) { setsid();}
  

extern void  plot_1(double xMin, double xMax, int dim,
                       double *f, double *ff,char* upstr, 
                       char* xstr, char* ystr);
extern int blind;

  
static int newPID=0;
static int pidList[100];

extern char pathtocalchep[], pathtohelp[];
  
void displayPlotN(char * title, double xMin, double xMax,  char*xName,  int dim, int N, double**f,double**ff,char**Y)
{ int pid;
  
  if(First) { First=0;   signal(SIGUSR1, disconnect);}  
  
  pid=fork();
  if(pid==0) 
  {  int err,i;
     va_list ap;   
     
     blind=0; 
     err=start1("micrOMEGAs Plot",NULL ,"calchep.ini",NULL);
     if(err) 
     { printf("Can not display plot because micromegas is compiled without X11\n");
       exit(0);
     }
     sprintf(pathtocalchep,"%s/",calchepDir);
     sprintf(pathtohelp,"%s/help/",pathtocalchep);
     clearTypeAhead();     
     plot_Nar(title,xMin,xMax,xName, dim, N, f,ff,Y);
     finish();
     exit(0);
  } else pidList[newPID++]=pid;
}

void displayPlot(char * title, double xMin, double xMax,  char*xName,  int dim, int N, ...)
{
  int i;
  double **f; double**ff; char**Y;  
  va_list ap;   
        
  f =malloc(N*sizeof(double*));
  ff=malloc(N*sizeof(double*));
  Y = malloc(N*sizeof(char*));
  
  va_start(ap,N);
  for(i=0;i<N;i++) 
  { f[i]=va_arg(ap,double*);
    ff[i]=va_arg(ap,double*);
    Y[i]=va_arg(ap,char*);
  }   
  va_end(ap);
  
  displayPlotN(title,xMin,xMax, xName, dim,N, f,ff,Y);

  free(f); free(ff);free(Y);
}


void  killPlots(void)
{
  int  C,i;
  
  for(i=0;i<newPID;i++) if(waitpid(pidList[i],NULL,WNOHANG)) break;

  if(newPID && i==newPID)
  {
    printf("Kill all plots (Y/N)? "); 
    C=getchar();
    if(C=='y'|| C=='Y')  kill(0,SIGKILL); else kill(0,SIGUSR1);
  }
  newPID=0;
}



void  killplots_(void) { killPlots();}

void displayFunc(double (*F)(double), double x1  ,double x2, char * mess)
{
  int i;
  double f[100];
  
  for(i=0;i<100;i++) f[i]=F(x1+(i+0.5)*(x2-x1)/100.);
  displayPlot(mess,x1,x2,"x",100,1,f,NULL,"F(x)");
}  

void displayFunc10(double (*F)(double), double x1  ,double x2, char * mess)
{
  int i;
  double f[100];
  
  for(i=0;i<100;i++) f[i]=F(pow(10,x1+(i+0.5)*(x2-x1)/99.));
  displayPlot(mess,x1,x2,"x",100,1,f,NULL,"F(10^x)");
}  
