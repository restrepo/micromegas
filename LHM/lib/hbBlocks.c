#include"../../include/micromegas.h"
#include"../../include/micromegas_aux.h"
#include "pmodel.h"


#define SQR(x) (x)*(x)
int  hbBlocksMDL(char*fname,int *nHch)
{ FILE * f=fopen(fname,"w");
  double Q;
  if(!f) return 0;
  Q=findValW("Q");

 fprintf(f,"Block Mass\n 25  %E # Higgs Mass\n\n",findValW("MH"));

  slhaDecayPrint("H",0,f);
  slhaDecayPrint("t",0,f);

  fprintf(f,"Block HiggsBoundsInputHiggsCouplingsBosons\n");
  fprintf(f,"# Effective coupling normalised to SM one and squared\n");
  fprintf(f,"# For (*) normalized on Sin(2*W)\n");
  fprintf(f," %12.4E  3    25    24    24 # higgs-W-W \n",       SQR(( 1.+findValW("del")*findValW("del")/12.)*
( 1.-findValW("del")*findValW("del")*findValW("vh")*findValW("vh")/findValW("v")/findValW("v")/3.)) );
  fprintf(f," %12.4E  3    25    23    23 # higgs-Z-Z \n",       SQR(( 1.+findValW("del")*findValW("del")/12.)*
( 1.-findValW("del")*findValW("del")*findValW("vh")*findValW("vh")/findValW("v")/findValW("v")/3.)) );
  fprintf(f," %12.4E  3    25    25    23 # higgs-higgs-Z \n",   0. );

  { double vev = 2*findValW("MW")*findValW("SW")/findValW("EE"),
    MH = findValW("MH"),
    aQCD=alphaQCD(MH)/M_PI,
    LGGSM=lGGhSM(MH,aQCD, findValW("Mcp"),findValW("Mbp"),findValW("Mtp"),vev),
    LAASM=lAAhSM(MH,aQCD, findValW("Mcp"),findValW("Mbp"),findValW("Mtp"),vev);


    fprintf(f," %12.4E  3    25    21    21 # higgs-gluon-gluon\n",  SQR(findValW("LGGH")/LGGSM) );
    fprintf(f," %12.4E  3    25    22    22 # higgs-gamma-gamma\n",  SQR(findValW("LAAH")/LAASM) );
  }


  fprintf(f,"Block HiggsBoundsInputHiggsCouplingsFermions\n");
  fprintf(f,"# Effective coupling normalised to SM one and squared\n");
  fprintf(f," %12.4E   %12.4E   3    25     5    5 # higgs-b-b \n"    ,SQR(1-findValW("vh")*findValW("vh")/findValW("f")/findValW("f")/4.) ,0.);
  fprintf(f," %12.4E   %12.4E   3    25     6    6 # higgs-top-top \n",SQR(1+findValW("mtcorr")/2/findValW("f")/findValW("f")*findValW("B00014")),0.);
  fprintf(f," %12.4E   %12.4E   3    25    15   15 # higgs-tau-tau \n",1.,0.);

  fclose(f);

  assignValW("Q",Q);
  calcMainFunc();
  if(nHch) *nHch=0;
  return 1;
}

int  hbblocksmdl_(char *fname, int * nHch,int len)
{ 
  char * cname=malloc(len+2);
  int nHiggs;  
  fName2c(fname,cname,len);    
  nHiggs= hbBlocksMDL(cname,nHch);      
  free(cname);        
  return nHiggs;          
}   