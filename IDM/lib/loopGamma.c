#include"../../include/micromegas.h"
#include"../../include/micromegas_aux.h"
#include"pmodel.h"


#include<sys/wait.h>
#include<unistd.h>

#define FIN  "nngg.in"
#define FOUT "nngg.out"

int loopGamma(double * csAZ, double *csAA)
{
  double sigmav;
  char buff[2000];
  int err,i;
  double Qstat;
  FILE*f;

  *csAA=0,*csAZ=0; 
  
  sprintf(buff,"%s/../lib/nnggidm/lGamma.exe",WORK);
  if(access( buff,X_OK))
  { char buf[2000]; 
    sprintf(buf, "dir=%s/../lib/nnggidm;  which  gmake; if(test $? -eq 0) then  gmake -C $dir; else make -C $dir; fi",WORK);
    system(buf);
  } 
  if(access( buff,X_OK)) 
  {  
    printf("Can not find/compile executable %s\n",buff);
    return 10;
  }  
   
  if(strcmp(CDM1,"~X")) 
  {
     printf("Error: loopGamma routine expects that DM is presented by  ~X, but for the current point DM=%s\n",CDM1); 
     return -2;  
  }
  if(Qaddress){ Qstat=*Qaddress; *Qaddress=2*Mcdm1; calcMainFunc();}
  
  f=fopen(FIN,"w");
  for(i=0;i<nModelVars;i++) fprintf(f,"%-6.6s   %f\n", varNames[i], varValues[i]);
  fprintf(f,"%-6.6s   %f\n", "GG", sqrt(4*M_PI*parton_alpha(2*Mcdm1/3.)));
  fprintf(f,"%-6.6s   %f\n", "wZ", pWidth("Z",NULL));
  fprintf(f,"%-6.6s   %f\n", "wW", pWidth("W+",NULL));
  fprintf(f,"%-6.6s   %f\n", "wh", pWidth("h",NULL)); 
  fclose(f); 
   
  if(Qaddress){ *Qaddress=Qstat; calcMainFunc();}  
     
  if(!access(FOUT,R_OK)) unlink(FOUT);
  
  sprintf(buff+strlen(buff)," %s %s",FIN,FOUT);
  err=System(buff);   
  double ee=findValW("EE");
  
  double aFac=4*M_PI/ee/ee/137.036;
  
//  printf("aFac=%E\n",aFac);
  
  if(err>=0) 
  {  err=slhaRead(FOUT,1);
     *csAZ=slhaVal("Lgamma",0.,1,1)*aFac*2.9979E-26;
     *csAA=slhaVal("Lgamma",0.,1,2)*aFac*aFac*2.9979E-26;
  }  

  return err;
}  

extern int  loopgamma_(double * cs1, double *cs2); /* fortran */
int         loopgamma_(double * cs1, double *cs2){ return loopGamma(cs1,cs2); }
 