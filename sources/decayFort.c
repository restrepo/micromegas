
#include"micromegas.h"
#include"micromegas_aux.h"



double pwidth_(char * name, txtList *L , int *dim, int len)
{
   char cName[10];
   fName2c(name,cName, len);
   return pWidth(cName, L , dim);
}

void printtxtlist_(txtList * list, int * Nch)
{
  char fname[20];
  FILE*f;
  sprintf(fname,"%d.tmptxt",getpid());
  f=fopen(fname,"w");

  printTxtList(*list,f);
  fclose(f);
  fortreread_(Nch,fname,strlen(fname));
  unlink(fname); 
}

double  findbr_(txtList*list, char * name,int len)
{
    char cName[30];
    fName2c(name,cName, len); 
    return  findBr(*list,cName);
}
