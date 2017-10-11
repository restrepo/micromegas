#include<math.h>
#include"pmodel.h"
#include"pmodel_aux.h"
#include"pmodel_f.h"
#include"../../include/micromegas.h"


void FillVal(int mode)
{ 
  double Q=91.;
  
  if(mode) 
  {  assignValW("MH3",slhaVal("MASS", Q,  1, 36));
     assignValW("mu", slhaVal("HMIX", Q,  1, 1) );
     assignValW("MG1",slhaVal("MSOFT", Q, 1, 1));
     assignValW("MG2",slhaVal("MSOFT", Q, 1, 2));
     assignValW("MG3",slhaVal("MSOFT", Q, 1, 3));
     assignValW("Ml1",slhaVal("MSOFT", Q, 1, 31));
     assignValW("Ml2",slhaVal("MSOFT", Q, 1, 32));
     assignValW("Ml3",slhaVal("MSOFT", Q, 1, 33));
     assignValW("Mr1",slhaVal("MSOFT", Q, 1, 34));
     assignValW("Mr2",slhaVal("MSOFT", Q, 1, 35));
     assignValW("Mr3",slhaVal("MSOFT", Q, 1, 36));
     assignValW("Mq1",slhaVal("MSOFT", Q, 1, 41));
     assignValW("Mq2",slhaVal("MSOFT", Q, 1, 42));
     assignValW("Mq3",slhaVal("MSOFT", Q, 1, 43));
     assignValW("Mu1",slhaVal("MSOFT", Q, 1, 44));
     assignValW("Mu2",slhaVal("MSOFT", Q, 1, 45));
     assignValW("Mu3",slhaVal("MSOFT", Q, 1, 46));
     assignValW("Md1",slhaVal("MSOFT", Q, 1, 47));
     assignValW("Md2",slhaVal("MSOFT", Q, 1, 48));
     assignValW("Md3",slhaVal("MSOFT", Q, 1, 49));
     assignValW("At", slhaVal("Au", Q, 2, 3, 3) );
     assignValW("Ab", slhaVal("Ad", Q, 2, 3, 3) );
     assignValW("Al", slhaVal("Ae", Q, 2, 3, 3) );
     assignValW("Am", slhaValExists("Ae",2,2,2)>0 ? slhaVal("Ae",Q,2,2,2):slhaVal("Al",Q,2,3,3));
     assignValW("Ad", slhaValExists("Ad",2,2,2)>0 ? slhaVal("Ad",Q,2,2,2):slhaVal("Ad",Q,2,3,3));
     assignValW("At", slhaVal("Au",Q,2,3,3));
     assignValW("Au", slhaValExists("Au",2,2,2)>0 ? slhaVal("Au",Q,2,2,2):slhaVal("Au",Q,2,3,3));
  }  
  if(mode>1)     
  {
     assignValW("alfSMZ",slhaVal("SMINPUTS",Q,1,3) );
     assignValW("MbMb",  slhaVal("SMINPUTS",Q,1,5) );
     assignValW("Mtp",   slhaVal("SMINPUTS",Q,1,6) );
     assignValW("Ml",    slhaVal("SMINPUTS",Q,1,7) );
  }
}