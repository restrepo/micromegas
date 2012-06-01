#include"../sources/micromegas.h"
#include"lib/pmodel.h"

static double SpA[NZ];
static int pdg=24;
static int out=0;


int main(int argc,char** argv)
{ int i;
  Mcdm=1000;
  char txt[40];
printf("argc=%d\n",argc);  
  if(argc!=4) { printf("program needs 3 arguments:\n"
                       " Mcdm  pdg  out\n");
                exit(1);
              }
  sscanf(argv[1],"%lf",&Mcdm);           
  sscanf(argv[2],"%d",&pdg);
  sscanf(argv[3],"%d",&out);
  
  basicSpectra(pdg,out,SpA);
  sprintf(txt,"Mcdm=%f, pdg=%d out=%d",Mcdm,pdg,out);
  displaySpectrum(SpA,txt,Mcdm*1.E-1,Mcdm, 1);
  if(pdg==24 || pdg==23)
  { double SpA_t[NZ],SpA_l[NZ];
    int i;
    basicSpectra(pdg+'T',out,SpA_t);
    sprintf(txt,"Mcdm=%f, pdg=%d+'T' out=%d",Mcdm,pdg,out);
    displaySpectrum(SpA_t,txt,Mcdm*1.E-1,Mcdm, 1);
    basicSpectra(pdg+'L',out,SpA_l);
    sprintf(txt,"Mcdm=%f, pdg=%d+'L' out=%d",Mcdm,pdg,out);
    displaySpectrum(SpA_l,txt,Mcdm*1.E-1,Mcdm, 1);
    for(i=0;i<NZ;i++) SpA[i]=(SpA_l[i]+2*SpA_t[i])/3;
    sprintf(txt,"Mcdm=%f, pdg=%d+<av> out=%d",Mcdm,pdg,out);
    displaySpectrum(SpA,txt,Mcdm*1.E-1,Mcdm, 1);    
  }        
  return 0; 
}
