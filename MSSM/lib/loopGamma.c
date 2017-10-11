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
  int err;
  FILE*f;
   
  *csAA=0,*csAZ=0; 

  if(!access(FOUT,R_OK)) unlink(FOUT);
  
  sprintf(buff,"%s/../lib/nngg12/lGamma.exe",WORK);
  if(access( buff,X_OK))
  { char buf[2000]; 
    sprintf(buf, "dir=%s/../lib/nngg12;  which  gmake; if(test $? -eq 0) then  gmake -C $dir; else make -C $dir; fi",WORK);
    system(buf);
  } 
  if(access( buff,X_OK)) 
  {  
    printf("Can not found/compile executable %s\n",buff);
    return 10;
  }  

  err=slhaWrite(FIN);
  if(err) return err; 
  f=fopen(FIN,"a");
  if(slhaDecayExists(36)<0) slhaDecayPrint("H3",0,f);
  if(slhaDecayExists(25)<0) slhaDecayPrint("h", 0,f);  
  if(slhaDecayExists(35)<0) slhaDecayPrint("H", 0,f);  
  fclose(f);    
  if(!access(FOUT,R_OK)) unlink(FOUT);
  
  sprintf(buff+strlen(buff)," %s %s",FIN,FOUT);
  err=System(buff);   
  
  if(err>=0) 
  {  err=slhaRead(FOUT,1);
     *csAZ=slhaVal("Lgamma",0.,1,1)*2.9979E-26;
     *csAA=slhaVal("Lgamma",0.,1,2)*2.9979E-26;
  }  

//  if(!access(FOUT,R_OK)) unlink(FOUT);
//  if(!access(FIN,R_OK)) unlink(FIN);
  return err;
}  

extern int  loopgamma_(double * cs1, double *cs2); /* fortran */
int         loopgamma_(double * cs1, double *cs2){ return loopGamma(cs1,cs2); }
 