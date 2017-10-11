#include <math.h>
#include "pmodel.h"
#include "../../include/micromegas.h"


void FillVal(int mode)
{ 
  double Q=0.;

  if(mode) 
  {  assignValW("MG1",slhaVal("EXTPAR", Q, 1, 1));
     assignValW("MG2",slhaVal("EXTPAR", Q, 1, 2));
     assignValW("MG3",slhaVal("EXTPAR", Q, 1, 3));
     assignValW("At", slhaVal("EXTPAR", Q, 1, 11));
     assignValW("Ab", slhaVal("EXTPAR", Q, 1, 12));
     assignValW("Al", slhaVal("EXTPAR", Q, 1, 13));
     assignValW("Ml2",slhaVal("EXTPAR", Q, 1, 32));
     assignValW("Ml3",slhaVal("EXTPAR", Q, 1, 33));
     assignValW("Mr2",slhaVal("EXTPAR", Q, 1, 35));
     assignValW("Mr3",slhaVal("EXTPAR", Q, 1, 36));
     assignValW("Mq2",slhaVal("EXTPAR", Q, 1, 42));
     assignValW("Mq3",slhaVal("EXTPAR", Q, 1, 43));
     assignValW("Mu2",slhaVal("EXTPAR", Q, 1, 45));
     assignValW("Mu3",slhaVal("EXTPAR", Q, 1, 46));
     assignValW("Md2",slhaVal("EXTPAR", Q, 1, 48));
     assignValW("Md3",slhaVal("EXTPAR", Q, 1, 49));
     assignValW("Mn2",slhaVal("EXTPAR", Q, 1, 58));
     assignValW("Mnlr",slhaVal("EXTPAR", Q, 1,59));
     assignValW("Alda",slhaVal("EXTPAR", Q, 1,63));
     assignValW("mu", slhaVal("EXTPAR", Q, 1, 65));
     assignValW("MZ2",slhaVal("EXTPAR", Q, 1,101));
     assignValW("aZZ",slhaVal("EXTPAR", Q, 1,102));
     assignValW("MK", slhaVal("EXTPAR", Q, 1,103));
     assignValW("M1p",slhaVal("EXTPAR", Q, 1,104));
     assignValW("tE6",slhaVal("EXTPAR", Q, 1,105));
  }  
  if(mode>1)     
  {
     assignValW("alfSMZ",slhaVal("SMINPUTS",Q,1,3) );
     assignValW("MbMb",  slhaVal("SMINPUTS",Q,1,5) );
     assignValW("Mtp",   slhaVal("SMINPUTS",Q,1,6) );
     assignValW("Ml",    slhaVal("SMINPUTS",Q,1,7) );
  }
}