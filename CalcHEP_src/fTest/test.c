#include<stdlib.h>
double *Q=NULL;
void xxx(void)
{ int i;
  Q=(double*) realloc(Q,10*sizeof(double));
  for(i=0;i<10;i++) Q[i]=0;
}
