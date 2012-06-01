/*
 Copyright (C) 1997 by Alexander Pukhov, email: pukhov@theory.npi.msu.su  
*/
#define real double
#include"chep_crt.h"
#include"interpret.h"
#include"syst2.h"
#include"procvar.h"
#include"n_compil.h"
#include"n_comphep_.h"
#include"out_int.h"
#include "getmem.h"

void interpretator(void)
{
 marktp   heapbg;
 void * pscr;
 double sqrts_mem=sqrts;
   get_text(1,1,maxCol(),maxRow(),&pscr);   
   clr_scr(FGmain,BGmain);
   
   allcanal=(canal*)m_alloc(subproc_sq*sizeof(*allcanal));
   mark_(&heapbg);

   if(compileall()==0) { n_comphep(); free(calcCoef); calcCoef=NULL;} 


   release_(&heapbg);
   free(allcanal);

   clr_scr(FGmain,BGmain);
   put_text(&pscr);
   sqrts=sqrts_mem;
}
