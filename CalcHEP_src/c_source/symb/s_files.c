/*
 Copyright (C) 1997, Alexander Pukhov 
*/

#include <unistd.h>
#include <dirent.h>
#include <pthread.h>

#include "n_proc.h"

#include"syst.h"
#include"syst2.h" 
#include"s_files.h"

 FILE * menup;
 FILE * menuq;
 FILE * diagrp;   /* file of Adiagram; */
 FILE * diagrq;   /* file of CSdiagram; */

 FILE * catalog;
 FILE * archiv;


char*mdFls[5] = {"vars","func","prtcls","lgrng","extlib"};
shortstr  pathtouser;

void wrt_menu(int menutype,int k,char*txt,int ndel,int ncalc,int nrest,long recpos)
{
   if (menutype == 1)
   {
      fseek(menup,(k - 1)*71 + 2,SEEK_SET);
      fprintf(menup,"%4d| %-44.44s|%5d|%5d|%-8ld",k,txt,ndel,nrest,recpos);
   }
   else
   {
      fseek(menuq,(k - 1)*77 + 2,SEEK_SET);
      fprintf(menuq,"%4d| %-44.44s|%5d|%5d|%5d|%-8ld",k,txt,ndel,ncalc,nrest,recpos);      
   }
}


int rd_menu(int menutype,int k,char*txt,int*ndel,int*ncalc,int*nrest,long*recpos)
{
   if (menutype == 1)
   {
      fseek(menup,(k - 1)*71 + 2,SEEK_SET);
      if(5!=fscanf(menup,"%d| %[^|]%*c%d|%d|%ld",&k,txt,ndel,nrest,recpos)) return 0;
      *ncalc=0;
   }
   else
   {
      fseek(menuq,(k - 1)*77 + 2,SEEK_SET);
      if(6!=fscanf(menuq,"%4d| %[^|]%*c%d|%d|%d|%ld",&k,txt,ndel,ncalc,nrest,recpos))
           return 0;
   }
   return 1;
}

int whichArchive(int nFile,int rw)
{ static int ArchNum=0, rw_=0;
  char archivName[40];  
  if(nFile==0)
  { if(ArchNum!=0) {ArchNum=0; fclose(archiv);}
    return 0;
  }
  
  if(ArchNum==nFile  && rw==rw_ )
  {  if(rw=='w' && ftell(archiv) >= MAXARCHIVE) 
     { fclose(archiv);
       ArchNum++;
       sprintf(archivName,ARCHIV_NAME,ArchNum);
       archiv=fopen(archivName,"w");
     } 
     return ArchNum;
  } 
  
  if(ArchNum) fclose(archiv); 
  ArchNum=nFile;  rw_=rw;
  sprintf(archivName,ARCHIV_NAME,ArchNum);
  if(rw=='w')
  { if(access(archivName,F_OK)==0) archiv=fopen(archivName,"a");
    else                           archiv=fopen(archivName,"w");
  }
  else  archiv=fopen(archivName,"r");
   
  return ArchNum;
} 

  typedef char nameRec[12];
  static nameRec * list=NULL;
  static int nFiles=0,cFile;


static pthread_mutex_t keyN=PTHREAD_MUTEX_INITIALIZER;

static int dN,fN,N;

static void * pCompile_cycle(int * input)
{  
  int i,ibeg,iend;
  char csrc[400];
  char osrc[400];
  char command[800];
  
    
  for(;;)
  { char *c;
    pthread_mutex_lock(&keyN);
    if(N>=fN) {  pthread_mutex_unlock(&keyN); return NULL;}
    ibeg=N;
    N+=20;
    iend=N;
    pthread_mutex_unlock(&keyN);
    if(iend>fN) iend=fN;
    csrc[0]=0; osrc[0]=0;
    for(i=ibeg+1;i<=iend;i++) sprintf(csrc+strlen(csrc)," f%d.c",i); 
    strcpy(osrc,csrc);
    for(c=strstr(osrc,".c"); c; c=strstr(c,".c")) c[1]='o';
     
    sprintf(command," cSrc=\"%s\"\n oSrc=\"%s\"\n"
                    " $CC $CFLAGS -c -I$CALCHEP/include $cSrc \n rm $cSrc\n ar r proclib_%d.a  $oSrc\n rm $oSrc\n $RANLIB proclib_%d.a",
        csrc, osrc, input[0], input[0]);
    system(command);
  } 
  return NULL;
}

static void* comp_control(void * par_)
{  
  infoLine(0.);
  for(;;)
  { double r;
    sleep(1);
    pthread_mutex_lock(&keyN);
    if(N>=fN) { pthread_mutex_unlock(&keyN); infoLine(2.);  return NULL; }
    r=(double)N/(double)fN;
    pthread_mutex_unlock(&keyN);
    infoLine(r);
  }   
}


int pCompile(void)
{
  DIR *dirPtr=opendir("./");
  struct dirent * dp;
  int i,k,n;
  char command[200], csrc[50], osrc[40];
  char *c;
    
  if(!dirPtr) return 1;
  
  fN=dN=0;
  while((dp=readdir(dirPtr)))
  { char *c=strstr(dp->d_name,".c");
    if(c&&c[2]==0) 
    { if(sscanf(dp->d_name,"f%d",&n)==1 && n>fN) fN=n;
      if(sscanf(dp->d_name,"d%d",&n)==1 && n>dN) dN=n;
    }  
  }
  closedir(dirPtr);
  N=0;
  if(fN)
  {  pthread_t *threads = malloc(nPROCSS*sizeof(pthread_t));
     int*num=malloc(nPROCSS*sizeof(int));
     pthread_t control;  
     for (k=0;k<nPROCSS;k++) num[k]=k+1;  
     for (k=0;k<nPROCSS;k++) pthread_create(threads+k,NULL,pCompile_cycle, num+k);
      comp_control(NULL);     
     for (k=0;k<nPROCSS;k++)   pthread_join(threads[k],NULL);
     free(threads);
  }
  
  if(dN)strcpy(csrc," d*.c"); else csrc[0]=0;
  if(access("service.c",R_OK)==0) strcat(csrc," service.c");
  if(access("sqme.c",R_OK)==0)    strcat(csrc," sqme.c"); 
  if(access("VandP.c",R_OK)==0)   strcat(csrc," VandP.c");
  strcpy(osrc,csrc);
  for(c=strstr(osrc,".c"); c; c=strstr(c,".c")) c[1]='o';

  if(strlen(csrc))
  {        
    sprintf(command,"csrc=\"%s\"\n"
                    "osrc=\"%s\"\n"
                    "$CC $CFLAGS -c -I$CALCHEP/include $csrc\n" 
                    "ar r proclib_0.a $osrc\n"
                    "$RANLIB proclib_0.a \n" 
                    "rm $csrc $osrc ",csrc,osrc);
    system(command);
  }                 
  system(". ./EXTLIBsh\n"
         "$CALCHEP/sbin/ld_n");
  
  return 0;
}
