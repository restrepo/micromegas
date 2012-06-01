/*
 Copyright (C) 1997, Alexander Pukhov 
*/

#include "chep_crt.h"
#include "files.h"
#include "plot.h"
#include <math.h>

int main(int argc,char** argv)
{  
 char   procName[100],xName[60],yName[60],buff[200];
 double xMin,xMax;
 int xDim=1;
 double yMin,yMax;
 int yDim=1;
 int ptype=-1; 
 double *f=NULL, *df=NULL;

 int i; 
 int n;
 char icon_name[STRSIZ];
 blind=0;

 strcpy(pathtocomphep,argv[0]);
 n = strlen(pathtocomphep)-1;
 while (n>=0 &&  pathtocomphep[n] != f_slash) n--;
 pathtocomphep[n-3]=0;
 sprintf(pathtohelp,"%shelp%c",pathtocomphep,f_slash);

 while(1==fscanf(stdin,"%[^\n]%*c",  buff))
 { char word[40];
   if(buff[0]!='#') break;
   if(1==sscanf(buff+1,"%s",word))
   {      if(strcmp(word,"type")==0)   sscanf(buff+1,"%*s %d",&ptype);  
     else if(strcmp(word,"type")==0)   sscanf(buff+1,"%*s %d",&ptype);  
     else if(strcmp(word,"xMin")==0)   sscanf(buff+1,"%*s %lf",&xMin);  
     else if(strcmp(word,"xMax")==0)   sscanf(buff+1,"%*s %lf",&xMax);  
     else if(strcmp(word,"yMin")==0)   sscanf(buff+1,"%*s %lf",&yMin);  
     else if(strcmp(word,"yMax")==0)   sscanf(buff+1,"%*s %lf",&yMax);  
     else if(strcmp(word,"xDim")==0)   sscanf(buff+1,"%*s %d",&xDim);  
     else if(strcmp(word,"yDim")==0)   sscanf(buff+1,"%*s %d",&yDim);    
     else if(strcmp(word,"title")==0)  
     { for(i=strlen(buff);buff[i-1]==' '; i--);
       sscanf(buff+1,"%*s %[^\n]",procName);  
     }
     else if(strcmp(word,"xName")==0) 
     { for(i=strlen(buff);buff[i-1]==' '; i--);
       sscanf(buff+1,"%*s %[^\n]",xName);  
     }   
     else if(strcmp(word,"yName")==0) 
     { for(i=strlen(buff);buff[i-1]==' '; i--);
       sscanf(buff+1,"%*s %[^\n]",yName);  
     }
   } else break;
   
 }  
  
 f=(double*)malloc(xDim*yDim*sizeof(double));   
 if(ptype)
 { 
   df=(double*)malloc(xDim*yDim*sizeof(double)); 
   sscanf(buff,"%lf  %lf",f,df);  
   for(i=1;i<xDim*yDim;i++) fscanf(stdin,"%lf  %lf",f+i,df+i);
   for(i=1;i<xDim*yDim;i++) if( !finite(f[i])|| !finite(df[i]) )
      { printf(" NAN in table %s\n",procName); return 0;}
 }
 else 
 { 
   sscanf(buff,"%lf",f);
   for(i=1;i<xDim*yDim;i++) fscanf(stdin,"%lf",f+i);
   for(i=1;i<xDim*yDim;i++) if( !finite(f[i])){ printf(" NAN in table %s\n",procName); return 0;} 
 }             
 sprintf(icon_name,"%sicon",pathtocomphep);
 start1(version,icon_name,"calchep.ini;../calchep.ini",NULL);  
 clearTypeAhead();
 if(ptype==2) plot_2(xMin,xMax,xDim,yMin,yMax,yDim,f,df,procName,xName,yName);
 else plot_1(xMin,xMax,xDim,f,df,procName, xName, yName);
 finish();
 return 0;
}
