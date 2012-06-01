/*
 Copyright (C) 1997, Alexander Pukhov 
*/

#include <unistd.h>

#include "chep_crt.h"
#include "files.h"


 char * outputDir = "";

 char  pathtocomphep[STRSIZ];

 char  version[36]="CalcHEP version 2.5b    ";

void copyfile(char* namefrom,char* nameto)
{ 

  FILE * filefrom;
  FILE * fileto;
  char  s [STRSIZ];
  filefrom= fopen(namefrom,"r");
  fileto  = fopen(nameto,"w");
  if ( (filefrom==NULL) || (fileto==NULL) ) return;

  for(;;)
  if (fgets(s,STRSIZ,filefrom) != NULL) f_printf(fileto,s) ; else
  {
     fclose(fileto);
     fclose(filefrom);
     return;
  }
}

void nextFileName(char* f_name,char * firstname,char * ext)
{  int tabnum=0;
   FILE * f;

   for(;;)
   {  tabnum++;
      sprintf(f_name,"%s%s%d%s",outputDir,firstname,tabnum,ext);
      f = fopen(f_name,"r");
      if(f) fclose(f); else return;
   } 
}

