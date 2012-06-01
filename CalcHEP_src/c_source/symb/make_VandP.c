#include<stdio.h>
#include <string.h>
#include <math.h>
#include "model.h"
#include "read_mdl.h"
#include "procvar.h"
#include "reader_c.h"
#include "parser.h"

extern int depQ1;

static void  readEXTFunc(FILE*f)
{ char buff[200];
  for(;fgets(buff,199,f);)
  { 
    trim(buff);
    if(strstr(buff,"extern ")==buff);
    { char *c;
      c=strchr(buff,'(');
      if(c)
      {
         c[0]=0;c--;
         while(c[0]==' ') c--;
         while(c[0]!=' ' && c[0]!='*' ) c--;
         c[0]=' ';
         EXTFunc=realloc(EXTFunc, strlen(EXTFunc)+3+strlen(c));
         sprintf(EXTFunc+strlen(EXTFunc),"%s ",c);
      } 
    }   
  }     
}   

static void  readModelFunc(FILE*f)
{ linelist  ln;
  fprintf(f,"/*  Special model functions  */\n"); 
  for(ln=modelTab[4].strings ; ln; ln=ln->next)
  { char buff[100];
    if(sscanf(ln->line,"%[^%\n|]", buff)!=1) continue;
    trim(buff);
    if(strstr(buff,"extern ")==buff);
    { char *c;
      c=strchr(buff,'(');
      if(c)
      {
         fprintf(f,"%s\n",buff);
         c[0]=0;c--;
         while(c[0]==' ') c--;
         while(c[0]!=' ' && c[0]!='*' ) c--;
         c[0]=' ';
         EXTFunc=realloc(EXTFunc, strlen(EXTFunc)+3+strlen(c));
         sprintf(EXTFunc+strlen(EXTFunc),"%s ",c);
      } 
    }   
  }
  fprintf(f,"\n");     
}   



int main(int argv, char**argc)
{
  int i,i10,nv,nLn,L;
  FILE*f,*fExt;  
  int nVar=0,nFunc=0,first; 
  char path[200];
  char * CalcHEP=NULL;
  
  if(argv!=3) { printf("Arguments expected: 1)path to model files; 2) model number.\n"); return 1;}

  L=strlen(argc[0]);
  CalcHEP=malloc(L+10);
  strcpy(CalcHEP,argc[0]);
  CalcHEP[L-15]=0;
  sscanf(argc[1],"%s",path);
  if(sscanf(argc[2],"%d",&L)!=1) { printf("Second argument should be a number\n"); return 1;}

  readModelFiles(path,L);
  if(!loadModel(0,1)) {printf("Error in model\n"); return 3;}


  f=fopen("VandP.c","w");
  fprintf(f,"#include <stdio.h>\n");
  fprintf(f,"#include <stdlib.h>\n");
  fprintf(f,"#include <string.h>\n");
  fprintf(f,"#include \"%s/include/extern.h\"\n", CalcHEP);
  fprintf(f,"#include \"%s/include/V_and_P.h\"\n",CalcHEP);
  fprintf(f,"#include \"autoprot.h\"\n");
  fprintf(f,"extern int  FError;\n");  

  if(EXTFunc) free(EXTFunc);
  EXTFunc=malloc(2); EXTFunc[0]=' '; EXTFunc[1]=0;
  sprintf(CalcHEP+strlen(CalcHEP),"/include/extern.h");  
  fExt=fopen(CalcHEP,"r");
  readEXTFunc(fExt);
  fclose(fExt);
  readModelFunc(f);
  ext_h=fopen("autoprot.h","w");

  for(nLn=0,i=0;i<nparticles;i++)
     if(!strchr("*fcCtT",prtclbase[i].hlp)&&i+1<=prtclbase[i].anti)nLn++;

  fprintf(f,"int nModelParticles=%d;\n",nLn);
  fprintf(f,"ModelPrtclsStr ModelPrtcls[%d]=\n{\n",nLn);

  for(nLn=0,i=0;i<nparticles;i++) if(!strchr("*fcCtT",prtclbase[i].hlp))
  { int anti=prtclbase[i].anti; 
    if(i+1>anti)  continue;
    if(nLn)fprintf(f,","); else fprintf(f," "); nLn++; 
      fprintf(f," {\"%s\",",prtclbase[i].name);
      if(i+1==anti)   fprintf(f,"\"%s\", ",prtclbase[i].name);
           else       fprintf(f,"\"%s\", ",prtclbase[anti-1].name);
     fprintf(f,"%ld, \"%s\",\"%s\",%d,%d,%d}\n",
       prtclbase[i].N,  prtclbase[i].massidnt, prtclbase[i].imassidnt,  
       prtclbase[i].spin, prtclbase[i].cdim,prtclbase[i].q3);
    
  }
  fprintf(f,"};\n");

  if (vararr) free(vararr);
  vararr = (singlevardescription*)m_alloc((nCommonVars+1)
                                            * sizeof(singlevardescription));

  sprintf(vararr[0].alias,"XXX");
  vararr[0].tmpvalue=vararr[0].num=vararr[0].used = 0;

  for(i=1;i<=nCommonVars;i++)
  {  if(strcmp(modelvars[i].varname,"i")==0) 
     { sprintf(vararr[i].alias,"I");
       vararr[i].tmpvalue=0;
       vararr[i].num=0;
       vararr[i].used = 0;
     }else
     { sprintf(vararr[i].alias,"V[%d]",nVar+nFunc);
       vararr[i].tmpvalue=modelvars[i].varvalue;
       vararr[i].num=nVar+nFunc;
       vararr[i].used = 1;
       if(modelvars[i].func) nFunc++; else nVar++;
     }
  }

  fprintf(f,"int nModelVars=%d;\n",nVar);
  fprintf(f,"int nModelFunc=%d;\n",nFunc);
  fprintf(f,"int LastVar=%d;\n",nVar);
  fprintf(f,"char*varNames[%d]={\n ", nVar+nFunc);
 
  for(first=1,i10=0, i=0;i<=nCommonVars;i++) if(vararr[i].used) 
  { if(first)  first=0; else fprintf(f,",");
    fprintf(f,"\"%s\"",modelvars[i].varname);
    if(++i10==10) {fprintf(f,"\n"); i10=0;}
  }
  fprintf(f,"};\n");
  fprintf(f,"double varValues[%d]={\n ", nVar+nFunc);
  for(first=1,nv=0,i10=0,i=0 ;i<=nCommonVars;i++) if(vararr[i].used) 
  { if(first)  first=0;else fprintf(f,",");
    fprintf(f,"%14.6E",modelvars[i].varvalue);
    if(++i10==10) { fprintf(f,"\n"); i10=0; }
    if(++nv==nVar) break;
  }
  fprintf(f,"};\n");
  

  fprintf(f,"int calcMainFunc(void)\n{\n");
  fprintf(f,"   int i;\n");
  fprintf(f,"   static double * VV=NULL;\n");
  fprintf(f,"   static int iQ=-1;\n");
  fprintf(f,"   static int cErr=0;\n");
  fprintf(f,"   double *V=varValues;\n");
  fprintf(f,"   FError=0;\n");
  fprintf(f,"   if(VV&&!cErr)\n"); 
  fprintf(f,"   { for(i=0;i<nModelVars;i++) if(i!=iQ && VV[i]!=V[i]) break;\n");
  fprintf(f,"     if(i==nModelVars) ");

  if(depQ1<=nCommonVars)
  { 
    fprintf(f,"     {if(iQ>=0 && VV[iQ]!=V[iQ]) goto FirstQ; else return 0;} \n");
  } else fprintf(f,"    return 0;\n");
  
  fprintf(f,"   }\n");
  
/*  fprintf(f," printf(\" callMainFunc\\n\");\n");  */

  fprintf(f,"  cErr=1;\n");
  for(i=1;i<=nCommonVars;i++)
  {
     if (vararr[i].used &&  modelvars[i].func)
     { checkNaN=0;
if(i==depQ1) fprintf(f," FirstQ:\n  cErr=1;\n");     
        {  char * ss=(char *)readExpression(modelvars[i].func,rd_c,act_c,free);
           fprintf(f,"   %s=%s;\n",vararr[i].alias,ss+3);
           free(ss);
           fprintf(f," LastVar=%d; ",vararr[i].num);
        }
        if(checkNaN)
        fprintf(f,"   if(!finite(%s) || FError) return %d;\n",vararr[i].alias,vararr[i].num);
        else fprintf(f,"\n");
     }
  }

  fprintf(f,"   if(VV==NULL) \n");
  fprintf(f,"   {  VV=malloc(sizeof(double)*nModelVars);\n");
  fprintf(f,"      for(i=0;i<nModelVars;i++) if(strcmp(varNames[i],\"Q\")==0) iQ=i;\n");
  fprintf(f,"   }\n");
  fprintf(f,"   for(i=0;i<nModelVars;i++) VV[i]=V[i];\n");
  fprintf(f,"   cErr=0;\n");
  fprintf(f,"return 0;\n}\n");


  fclose(f);
  fclose(ext_h);
  return 0;
}
