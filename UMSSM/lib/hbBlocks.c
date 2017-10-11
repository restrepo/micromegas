#include"../../sources/micromegas.h"
#include"../../sources/micromegas_aux.h"
#include "pmodel.h"

#define SQR(x) (x)*(x)
int  HBblocks(char * fname)
{ FILE * f=fopen(fname,"a");
  int pdg,i;
  char *h[3]={"h1","h2","h3"};

  if(!f) return 1;

  fprintf(f,"Block HiggsBoundsInputHiggsCouplingsBosons\n");
  fprintf(f,"# Effective coupling normalised to SM one and squared\n");
  fprintf(f,"# For (*) normalized on Sin(2*W)\n"); 
  for(i=0;i<3;i++)
  { pdg=pNum(h[i]);
    fprintf(f," %12.4E  3  %d    24    24 # %s-W-W \n",       SQR(slhaVal("REDCOUP",0.,2,1+i,4)),pdg,h[i]);
    fprintf(f," %12.4E  3  %d    23    23 # %s-Z-Z \n",       SQR(slhaVal("REDCOUP",0.,2,1+i,5)),pdg,h[i]);
    fprintf(f," %12.4E  3  %d    21    21 # %s-gluon-gluon\n",SQR(slhaVal("REDCOUP",0.,2,1+i,6)),pdg,h[i]);
    fprintf(f," %12.4E  3  %d    22    22 # %s-gamma-gamma\n",SQR(slhaVal("REDCOUP",0.,2,1+i,7)),pdg,h[i]);
  }
  { pdg=pNum("ha");
    fprintf(f," %12.4E  3  %d    24    24 # %s-W-W \n",       SQR(slhaVal("REDCOUP",0.,2,4,4)),pdg,"ha");
    fprintf(f," %12.4E  3  %d    23    23 # %s-Z-Z \n",       SQR(slhaVal("REDCOUP",0.,2,4,5)),pdg,"ha");
    fprintf(f," %12.4E  3  %d    21    21 # %s-gluon-gluon\n",SQR(slhaVal("REDCOUP",0.,2,4,6)),pdg,"ha");
    fprintf(f," %12.4E  3  %d    22    22 # %s-gamma-gamma\n",SQR(slhaVal("REDCOUP",0.,2,4,7)),pdg,"ha");
  }

  fprintf(f,"Block HiggsBoundsInputHiggsCouplingsFermions\n");
  fprintf(f,"# Effective coupling normalised to SM one and squared\n");
  for(i=0;i<3;i++)
  { pdg=pNum(h[i]);
    fprintf(f," %12.4E  0.  3  %d   5    5  # %s-b-b    \n",SQR(slhaVal("REDCOUP",0.,2,1+i,3)),pdg,h[i]);
    fprintf(f," %12.4E  0.  3  %d   6    6  # %s-top-top\n",SQR(slhaVal("REDCOUP",0.,2,1+i,1)),pdg,h[i]);
    fprintf(f," %12.4E  0.  3  %d  15   15  # %s-tau-tau\n",SQR(slhaVal("REDCOUP",0.,2,1+i,2)),pdg,h[i]);
  }
  { pdg=pNum("ha");
    fprintf(f," 0.  %12.4E  3  %d   5    5  # %s-b-b    \n",SQR(slhaVal("REDCOUP",0.,2,4,3)),pdg,"ha");
    fprintf(f," 0.  %12.4E  3  %d   6    6  # %s-top-top\n",SQR(slhaVal("REDCOUP",0.,2,4,1)),pdg,"ha");
    fprintf(f," 0.  %12.4E  3  %d  15   15  # %s-tau-tau\n",SQR(slhaVal("REDCOUP",0.,2,4,2)),pdg,"ha");
  }

  fclose(f);
  return 0;
}
