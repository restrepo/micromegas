#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"

int main(int argc,char** argv)
{ int err,n,i;
  SpectraFlag=0;      
  double pi0[NZ],pich[NZ],kl[NZ],ks[NZ],kch[NZ];
  double M=1;
  basicSpectra(M,111,0,pi0);
  basicSpectra(M,211,0,pich);
  basicSpectra(M,130,0,kl);
  basicSpectra(M,310,0,ks);
  basicSpectra(M,321,0,kch);
  displayPlot(" Edgamma/dE","E",0.001,M ,1,5,"p0",0, eSpectdNdE,pi0
                                            ,"pich",0, eSpectdNdE,pich
                                            ,"KL",0, eSpectdNdE,kl
                                            ,"Ks",0, eSpectdNdE,ks
                                            ,"K+/-",0, eSpectdNdE,kch
                                          );  
}
