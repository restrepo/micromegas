#include"../../sources/micromegas_aux.h"
#include"pmodel_f.h"
#include"pmodel.h"


void o1contents_(int *file)
{
  char fname[20];
  FILE*f;

  sprintf(fname,"%d.tmptxt",getpid());
  f=fopen(fname,"w");
  o1Contents(f);
  fclose(f);
  fortreread_(file,fname,strlen(fname));
  unlink(fname);
}
