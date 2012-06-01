#ifndef __PROCVAR_
#define __PROCVAR_

#include "physics.h"

typedef struct
   {
      char    alias[15];
      double  tmpvalue;
      int     num; 
      int     used;
   }  singlevardescription;

extern singlevardescription *vararr;
extern int nProcessVar;

extern int  initvararray(int nsub, char key,int width);

#endif
