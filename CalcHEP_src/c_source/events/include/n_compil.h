#ifndef __N_COMPIL_
#define __N_COMPIL_

#include"syst2.h"
#include "sets.h"
#include "optimise.h"
#include "physics.h"
#include "denominators.h"

typedef struct prgcoderec
   {
      struct prgcoderec *   next;
      varptr  totn, totd, rnum;
      int     denorno;
      int     power[2*MAXINOUT-6];
      int     order_num[2*MAXINOUT-6];
      int     width[2*MAXINOUT-6];
   }  prgcoderec;

typedef struct prgcoderec *prgcodeptr;

typedef struct canal
   {
      char         prtclnames[MAXINOUT][6];
      double *      prtclmasses[MAXINOUT];
      prgcodeptr   codeptr;
      denlist      denominators;
   }  canal;
   
extern int  compileall(void);
extern  struct canal *  allcanal;
#endif


