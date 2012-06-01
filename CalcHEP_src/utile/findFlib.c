#include<stdio.h>
#include<stdlib.h>
#include<string.h>

typedef struct  {void * next; char * word;} list;

static list * getLink(FILE*f)
{
  char word[200];
  list * flags=NULL,*new;
  while(fscanf(f,"%s",word)==1)
  { if(word[0]=='-' && strlen(word)>2 && toupper(word[1])=='L')
    {
      new=(list*)malloc(sizeof(list));
      new->next=flags;
      new->word=(char*)malloc(1+strlen(word));
      strcpy(new->word,word);
      flags=new;
    } 
  } 
  return flags;  
}

static void printlist(list*l) {for(;l;l=l->next) printf("%s\n",l->word);} 

int main(int argc, char**argv ) 
{ FILE *f;
  char command[100];
  int err;
  list*flagsC,*flagsF,*flagsD,*i;
  if(argc!=5)
  { printf("This progran needs 4 arguments like: cc -v f77 -v\n");
    return 1;
  }
  system("echo \"int main(void){return 1;}\">test.c");
  sprintf(command,"%s  %s test.c -o test.exe 2>Err 1>Err",argv[1],argv[2]);
  err=system(command);
  if(err) 
  {  printf("Error in execution \n %s\nsee file Err for details\n",command);
     return 1;
  } 
  
  f=fopen("Err","r");
  flagsC=getLink(f);
  fclose(f);
  if(flagsC==NULL) 
  { printf("Warning: No -L -l flags were founded,\n    it looks like %s "
     "is not the appropriate flag\n", argv[2]); 
  }
  system("echo \"       write(*,*) 'Ok'\n      end\n\">test.f");
  sprintf(command,"%s  %s test.f -o test.exe 2>Err 1>Err",argv[3],argv[4]);
  err=system(command);
  if(err) 
  {  printf("Error in execution \n %s\nsee file Err for details\n",command);
     return 1;
  } 
  
  f=fopen("Err","r");
  flagsF=getLink(f);
  fclose(f);
  if(flagsF==NULL) 
  { printf("Warrning: No -L -l flags were founded,\n     it looks like %s "
     "is  not the appropriate flag\n", argv[4]); 
  }
 
  flagsD=NULL;
  for(i=flagsF;i;i=i->next)
  { list *j,*k;
     for(j=flagsC;j;j=j->next) if(strcmp(j->word,i->word)==0) break;
     if(j==NULL) 
     { for(k=flagsD;k;k=k->next) if(strcmp(k->word,i->word)==0) break;
       if(k==NULL) 
       { list*new=(list*)malloc(sizeof(list));
         new->next=flagsD;
         new->word=(char*)malloc(1+strlen(i->word));
         strcpy(new->word,i->word);
         flagsD=new;
       }  
     }
  }    
  for(i=flagsD;i;i=i->next)
  printf("%s\n",i->word);
  system("rm test.c test.f test.exe Err"); 
  return 0;
}
