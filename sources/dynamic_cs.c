#include "micromegas.h"
#include "micromegas_aux.h"
#include "micromegas_f.h"

numout*newProcess(char*Process,char*lib)
{  int i;
   char * reserved[3]={"omg", "pwidth", "dir"};
   
   for(i=0;  Process[i]&& Process[i]!=' '  ;i++);   

   if(Process[i])  for(i=0;i<3;i++)
   {
     char * c=strstr(lib,reserved[i]);
     if(c==lib)
     {  printf("Library names started from '%s' are reserved for internal use\n",
                   reserved[i]);
        if(i==0) printf("The libraries of WIMP annihilation  become available automatically\n"
                      "after calculation of relic density.\n");              
        return NULL;
      }
   } 
   return getMEcode(0,ForceUG,Process,NULL,"",lib);
}


void newprocess_(char*Process,char*lib, int * address, int len1,int len2)
{ char cProcess[100], clib[100];
  numout*cc;
  fName2c(Process,cProcess,len1);
  fName2c(lib,clib,len2);
  cc=newProcess(cProcess,clib);
  if(!cc) *address=0;else memcpy(address,&cc,sizeof(cc)); 
}
  
     
void  procinfo1_(int*ccf, int *ntot, int * nin, int *nout)
{  numout*cc;
   memcpy(&cc,ccf,sizeof(cc));
   procInfo1(cc, ntot, nin, nout);}

void procinfo2_(int*ccf,int*nsub,char*name,double*mass,int len)
{ int ntot, nin, nout,i;
  numout*cc;
  char ** cname;  
  memcpy(&cc,ccf,sizeof(cc));
  procInfo1(cc, &ntot, &nin, &nout);
  cname=malloc((nin+nout)*sizeof(char*));
  
  procInfo2(cc,*nsub,cname, mass);
  for(i=0;i<nin+nout;i++) cName2f(cname[i],name+i*len,len);
}

