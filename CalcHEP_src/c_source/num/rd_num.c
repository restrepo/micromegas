/*
 Copyright (C) 1997, Alexander Pukhov 
*/
#include <stdio.h>
#include <ctype.h>
#include "interface.h"
#include "rd_num.h"

#include "V_and_P.h"
//#include "dynamic_cs.h"

int  rd_num(char* s, double *p)
{
  int n;   

  if (isdigit(*s)) { sscanf(s,"%lf",p); return 1;}

  for(n=0;n<nModelVars+nModelFunc;n++) if(strcmp(varNames[n],s)==0) 
  {*p=varValues[n];return 1;}
  
  
  if(strcmp("PI",s)==0) { *p=M_PI; return 1;}  
    
  return 0;   
}
