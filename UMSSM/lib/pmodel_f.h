#ifndef  __MSSM_F_
#define  __MSSM_F_

/*===============================================
      UMSSMTools
 =================================================*/

extern int umssmtools_(int *PDG_LSP);
/* 
     integer function umssmtools(PDG_LSP)
     int PDG_LSP
*/ 

/*===============================================
      Les Houches accord interface
 =================================================*/
extern int readlesh_(char *fname , int*SM,  int len );
/* 
     integer function readLesH(fname,SM)
     character*(*) fname
     int SM
*/  

extern int leshinput_(char *fname, int len);
/* 
     integer function leshinput(fname)
     character*(*) fname
*/ 

/*===============================================
      Higgs bounds
 =================================================*/
extern int  hbblocks_(char * fname,int len);
/* 
     integer function hbblocks(fname)
     character*(*) fname
*/ 

extern int lilithf_(char*fname,int len);
/* 
     integer function lilithf(fname)
     character*(*) fname
*/ 

void o1contents_(int *N);
/*
  subroutine  o1Contents(int *N )
*/
    
#endif

