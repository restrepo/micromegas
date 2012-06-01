#include"../../sources/micromegas.h"
#include"../../sources/micromegas_aux.h"
#include"pmodel.h"


#include<sys/wait.h>
#include<unistd.h>

#define FIN  "slha_for_LG.txt"
#define FOUT "LG_out.txt"

int loopgamma_(double * cs1, double *cs2)
{
  double sigmav;
  char buff[2000];
  int err;

  *cs1=0,*cs2=0; 

  if(!access(FOUT,R_OK)) unlink(FOUT);
  
  sprintf(buff,"%s/../nngg/lGamma.exe",WORK);
  if(access( buff,X_OK))
  { char buf[2000]; 
    sprintf(buf, "make -C %s/../nngg",WORK);
    system(buf);
  } 
  if(access( buff,X_OK)) 
  {  
    printf("Can not found/compile executable %s\n",buff);
    return 10;
  }  

  err=slhaWrite(FIN);
  if(err) return err; 

  if(!access(FOUT,R_OK)) unlink(FOUT);
  
  sprintf(buff+strlen(buff)," %s %s",FIN,FOUT);
  err=System(buff);   
  
  if(err>=0) 
  {  err=slhaRead(FOUT,1);
     *cs1=slhaVal("Lgamma",0.,1,1)*2.9979E-26;
     *cs2=slhaVal("Lgamma",0.,1,2)*2.9979E-26;
  }  

  if(!access(FOUT,R_OK)) unlink(FOUT);
  if(!access(FIN,R_OK)) unlink(FIN);
  return err;
}  
