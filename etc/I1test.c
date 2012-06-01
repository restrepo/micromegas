#include<math.h>
#include<stdlib.h>

#include"../sources/micromegas.h"
#include"../sources/micromegas_aux.h"

int main(int argc,char** argv)
{ 
  double mq,msq,mne;
  if(argc<4) 
  { printf("The program needs 3 arguments: mq,msq,mne\n");
    return 1;
  }
  sscanf(argv[1],"%lf",&mq);
  sscanf(argv[2],"%lf",&msq);
  sscanf(argv[3],"%lf",&mne);
  printf("BOX/TREE= %E, (1-0.5*mq/(msq-mne))=%E \n",LintIk(1,msq,mq,mne)*3*mq*mq*(msq*msq-mne*mne)/2,
  1-0.5*mq/(msq-mne));
  return 0;
}
