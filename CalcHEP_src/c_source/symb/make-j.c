#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <pthread.h>


static FILE*f=NULL;
static pthread_mutex_t keyN=PTHREAD_MUTEX_INITIALIZER;

static void * pCompile_cycle(int * input)
{  
  char command[100];
  char fName[12];    
  for(;;)
  { char fName[12];
    pthread_mutex_lock(&keyN);
      if(f && 1!=fscanf(f,"%s",fName)) { fclose(f); f=NULL; }    
    pthread_mutex_unlock(&keyN);
    if(!f) return NULL;
    sprintf(command," $CC -c $CFLAGS %s; if(test $? -eq 0) then  rm %s;fi\n",
    fName,fName);
    system(command);
  } 
  return NULL;
}

int main(int argc, char** argv)
{  
   DIR *dirPtr=NULL;
   struct dirent * dp=NULL;
   int Ntot=0;
   static int nParProc;

   if(argc!=2) return 1;
   if(1!=sscanf(argv[1],"%d",&nParProc)) return 2;
   if(nParProc<=0) nParProc=1;

   f=fopen("make-j_list.txt","w");
   
   dirPtr=opendir("./");   
   while((dp=readdir(dirPtr)))
   { char *c=strstr(dp->d_name,".c");
     if(c&&c[2]==0) { fprintf(f,"%s\n",dp->d_name); Ntot++;}
   }  
   closedir(dirPtr);
   fclose(f);
   if(!Ntot) return 0;
   
   f=fopen("make-j_list.txt","r");
   { pthread_t *threads = malloc(nParProc*sizeof(pthread_t));
     int k;  
     for (k=0;k<nParProc;k++) pthread_create(threads+k,NULL,pCompile_cycle,NULL);  
     for (k=0;k<nParProc;k++) pthread_join(threads[k],NULL);
     free(threads);
  }
   
  unlink("make-j_list.txt"); // closed by threads
  return 0;   
}
