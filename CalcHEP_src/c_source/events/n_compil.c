/*
 Copyright (C) 1997 by Alexander Pukhov, email: pukhov@theory.npi.msu.su 
*/

#include "os.h"
#include "syst2.h"
#include "optimise.h" 
#include "crt_util.h"
#include "procvar.h"
#include "getmem.h"
#include "out_serv.h"
#include "saveres.h"
#include "pvars.h"
#define INTERPRET
#include "out_ext.h"
#undef  INTERPRET
#include "n_compil.h"


int   nin_; 
int   nout_;
int   nprc_;

int   nvar_;
int   nfunc_;

struct canal *  allcanal=NULL;

static void  louddiagrexpr(prgcodeptr*  prg, deninforec * dendescript)
{prgcodeptr  pntr;
 int i;
 catrec cr; 
 polyvars  allVars={0,NULL};
    vardef=&allVars;
  
   fseek(catalog,dendescript->cr_pos,SEEK_SET);
   FREAD1(cr,catalog);
             
       
   pntr = *prg; 
   *prg = (prgcodeptr)getmem_(sizeof(prgcoderec)); 
   (*prg)->next = pntr; 

   fseek(archiv,cr.factpos,SEEK_SET); 
   readvardef();
   readpolynom(&(*prg)->totn);  
   readpolynom(&(*prg)->totd);
   clearvardef(); 

   fseek(archiv,cr.rnumpos,SEEK_SET);
   readvardef(); 
   readpolynom(&(*prg)->rnum); 
   clearvardef();

   (*prg)->denorno=dendescript->tot_den;

   for (i = 0; i < (*prg)->denorno ; i++) 
   { 
     (*prg)->power[i] = dendescript->denarr[i].power;
     (*prg)->order_num[i] = dendescript->denarr[i].order_num;  
     (*prg)->width[i] = dendescript->denarr[i].width;      
   } 
}


static void  compileprocess(int  nsub, prgcodeptr*  prg, denlist * denomi )
{

   FILE * fd;                /* file of (deninforec)  */
   char fd_name[STRSIZ];
   int nden_w,nden_0;        
   deninforec   dendescript;


   *prg = NULL;


   sprintf(fd_name,"%stmp%cden.inf",pathtouser,f_slash);
   fd=fopen(fd_name,"wb");
   denominatorStatistic(nsub, &nden_w, &nden_0, denomi, fd);
   fclose(fd);


   fd=fopen(fd_name,"rb"); 
   
   while(FREAD1(dendescript,fd) == 1) louddiagrexpr(prg,&dendescript);  
   
   fclose(fd);
   
  unlink(fd_name); 


}



int  compileall(void)
{int  n_sub, ndel, ncalc, nrest; long recpos; 
 int       i,k, pNum[MAXINOUT+1]; 
 char      procname[STRSIZ]; 
 infoptr      inf;
 
   outputLanguage='I';
   catalog=fopen(CATALOG_NAME,"rb");
   archiv =fopen(ARCHIV_NAME,"rb");

   initvararray(0,'R');

   firstVar=nmodelvar;
   if( !strcmp( modelvars[firstVar].varname,strongconst))  firstVar--;
      
   initinfo(); 

   menuq=fopen(MENUQ_NAME,"rb"); 
   nprc_=0; 
   for(n_sub=1;n_sub<=subproc_sq;n_sub++) 
   { 
      rd_menu(2,n_sub,procname,&ndel,&ncalc,&nrest,&recpos); 
      if (nrest != 0) 
      { 
          messanykey(10,10,
"All diagrams in the subprocess must be\ncalculated or marked as deleted"); 
          return -1;
      } 
      if (ncalc != 0) 
      {
         findPrtclNum(procname,pNum);
         for(i=0;i<nin+nout;i++)
         { strcpy(allcanal[nprc_].prtclnames[i], prtclbase[pNum[i+1]-1].name);
            allcanal[nprc_].prtclmasses[i]= 
              &(vararr[modelVarPos(prtclbase[pNum[i+1]-1].massidnt)].tmpvalue);
         }
         compileprocess(n_sub,&allcanal[nprc_].codeptr,
                              &allcanal[nprc_].denominators); 
         nprc_++;
      } 
   } 
   fclose(menuq); 
   fclose(catalog);
   fclose(archiv);


   if(nprc_ == 0) 
   { 
      messanykey(10,10,
      " Set of calculated diagrams for this processis is empty");  
      return -1;
   } 

   inf = info;

   while (inf != NULL)
   {
      if (inf->consttype == numb)
      {
         inf->rval = inf->ival;
         inf->consttype =  rnumb ;
      }
      inf = inf->next;
   }

   
   nin_ = nin; 
   nout_= nout;

   nvar_=0;
   nfunc_=0;
     
   for (k=1;k<=nmodelvar; k++)  
   {  if(vararr[k].used)
      { if( modelvars[k].func) nfunc_++;else  nvar_++;}
   } 
      
   return 0;
} 
