#include<string.h>
#include"micromegas.h"
#include"../CalcHEP_src/c_source/dynamicME/include/dynamic_cs.h"

int slhaDecayPrint(char * name,FILE*f)
{
   double w;
   txtList all;
   int dim;
   int i;
   int PDG;        
   PDG=qNumbers(name,NULL,NULL,NULL);
   if(!PDG) return 0;
   w=pWidth(name,&all,&dim);
   fprintf(f,"DECAY %d  %E  # %s\n",PDG,w,name);
   for(;all;all=all->next)
   {  
      char pn[20], buff[100], *chB,*chE;
      strcpy(buff,all->txt);
      sscanf(buff,"%s", pn);
      fprintf(f," %s %d  ",pn,dim);
      chB=strstr(buff,"->");
      chB+=2;
      for(i=0;i<dim;i++)
      { 
         chE=strchr(chB,',');
         if(chE)chE[0]=0;
         sscanf(chB,"%s",pn);
         fprintf(f," %d", qNumbers(pn,NULL,NULL,NULL));
         if(chE)chB=chE+1;else break;           
      }
      chB=strstr(all->txt,"->");
      fprintf(f,"  # %s \n",chB+2);
   } 
   fprintf(f,"\n");
   return PDG;
} 
