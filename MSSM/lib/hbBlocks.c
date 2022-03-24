#include"../../include/micromegas.h"
#include"../../include/micromegas_aux.h"
#include "pmodel.h"

#define SQR(x) (x)*(x)
int  hbBlocksMDL(char*fname,int *nHch)
{ if(nHch) *nHch=1;
  if(useSLHAwidth && 
       (  blockExists("HiggsBoundsInputHiggsCouplingsBosons")
        ||blockExists("HiggsBoundsHiggsCouplingsBosons")
       )
    )
  {   
       slhaWrite(fname);
       return 3;
  }        
  FILE * f=fopen(fname,"w");
   
  fprintf(f," Block MODSEL  # Model selection\n");
  fprintf(f,"  1    0       # general MSSM\n");
  fprintf(f,"  3    0       # MSSM particles\n");
      
  double tb,sb,cb,alpha,sa,ca,ta,samb,camb,dMb,MbHl,MbSM,MbH,MbH3;
  double vev= 2*findValW("MW")*findValW("SW")/findValW("EE"),
  Mcp=findValW("Mcp"),Mbp=findValW("Mbp"),Mtp=findValW("Mtp"),
  Mh=findValW("Mh"),MH=findValW("MH"),MH3=findValW("MH3");
  double LGGSM,LAASM; 
    
  fprintf(f,"BLOCK MASS\n");
  fprintf(f,"  25 %E  # h\n",pMass("h"));
  fprintf(f,"  35 %E  # h\n",pMass("H"));
  fprintf(f,"  36 %E  # h\n",pMass("H3"));
  fprintf(f,"  37 %E  # h\n",pMass("H+"));
  

  slhaDecayPrint("h", 0, f);
  slhaDecayPrint("H", 0, f);
  slhaDecayPrint("H3",0, f);
  slhaDecayPrint("t", 0, f);
  slhaDecayPrint("H+",0, f);

  tb=findValW("tB");  
    sb=tb/sqrt(1+tb*tb);
    cb=1/sqrt(1+tb*tb);
  alpha=findValW("alpha");
    sa=sin(alpha);
    ca=cos(alpha);
    ta=sa/ca;
    samb=sa*cb-ca*sb;
    camb=ca*cb+sa*sb;
  dMb=findValW("dMb");

  MbSM=findValW("Mb");
  MbH= MbSM/(1+dMb)*(1+dMb*ta/tb);  
  MbH3=MbSM/(1+dMb)*(1-dMb/tb/tb);
  MbHl=MbSM/(1+dMb)*(1-dMb/ta/tb);
 
  fprintf(f,"Block HiggsBoundsInputHiggsCouplingsBosons\n");
  fprintf(f,"# Effective coupling normalised to SM one and squared\n");
  fprintf(f,"# For (*) normalized on MZ/vev \n"); 
  fprintf(f," %12.4E  3    25    24    24 # higgs-W-W \n",        SQR(samb)  );
  fprintf(f," %12.4E  3    25    23    23 # higgs-Z-Z \n",        SQR(samb)  );
  fprintf(f," %12.4E  3    25    25    23 # higgs-higgs-Z \n",    0.   );

  LGGSM=lGGhSM(Mh,alphaQCD(Mh)/M_PI, Mcp,Mbp,Mtp,vev);
  LAASM=lAAhSM(Mh,alphaQCD(Mh)/M_PI, Mcp,Mbp,Mtp,vev);

    fprintf(f," %12.4E  3    25    21    21 # higgs-gluon-gluon\n",  SQR(findValW("LGGh")/LGGSM) );           
    fprintf(f," %12.4E  3    25    22    22 # higgs-gamma-gamma\n",  SQR(findValW("LAAh")/LAASM) );
 
  fprintf(f," %12.4E  3    35    24    24 # higgs-W-W \n",        SQR(camb)  );
  fprintf(f," %12.4E  3    35    23    23 # higgs-Z-Z \n",        SQR(camb)  );
  fprintf(f," %12.4E  3    35    25    23 # higgs-higgs-Z \n",    0.  );
  fprintf(f," %12.4E  3    35    35    23 # higgs-higgs-Z \n",    0.  );
  
  LGGSM=lGGhSM(MH,alphaQCD(MH)/M_PI, Mcp,Mbp,Mtp,vev);
  LAASM=lAAhSM(MH,alphaQCD(MH)/M_PI, Mcp,Mbp,Mtp,vev);
  
  fprintf(f," %12.4E  3    35    21    21 # higgs-gluon-gluon\n",SQR(findValW("LGGH")/LGGSM)  );   
  fprintf(f," %12.4E  3    35    22    22 # higgs-gamma-gamma\n",SQR(findValW("LAAH")/LAASM)  );  

  fprintf(f," %12.4E  3    36    24    24 # higgs-W-W \n",        0.  );
  fprintf(f," %12.4E  3    36    23    23 # higgs-Z-Z \n",        0.  );

  LGGSM=lGGhSM(MH3,alphaQCD(MH3)/M_PI, Mcp,Mbp,Mtp,vev);
  LAASM=lAAhSM(MH3,alphaQCD(MH3)/M_PI, Mcp,Mbp,Mtp,vev);
  
  fprintf(f," %12.4E  3    36    21    21 # higgs-gluon-gluon\n",SQR(2*findValW("LGGH3")/LGGSM) );
  fprintf(f," %12.4E  3    36    22    22 # higgs-gamma-gamma\n",SQR(2*findValW("LAAH3")/LAASM) );             
  
  fprintf(f," %12.4E  3    36    25    23 #*higgs-higgs-Z \n",    SQR(camb)  );
  fprintf(f," %12.4E  3    36    35    23 #*higgs-higgs-Z \n",    SQR(samb)  );
  fprintf(f," %12.4E  3    36    36    23 #* higgs-higgs-Z \n",   0.  );

  fprintf(f,"Block HiggsBoundsInputHiggsCouplingsFermions\n");
  fprintf(f,"# Effective coupling normalised to SM one and squared\n");
  fprintf(f," %12.4E   %12.4E   3    25     5    5 # higgs-b-b \n"    ,SQR((sa/cb)*(MbHl/MbSM)),0.);
  fprintf(f," %12.4E   %12.4E   3    25     6    6 # higgs-top-top \n",SQR(ca/sb)              ,0.);
  fprintf(f," %12.4E   %12.4E   3    25    15   15 # higgs-tau-tau \n",SQR(sa/cb)              ,0.);

  fprintf(f," %12.4E   %12.4E   3    35     5    5 # higgs-b-b \n"    ,SQR((ca/cb)*(MbH/MbSM))  ,0.);
  fprintf(f," %12.4E   %12.4E   3    35     6    6 # higgs-top-top \n",SQR(sa/sb)              ,0.);  
  fprintf(f," %12.4E   %12.4E   3    35    15   15 # higgs-tau-tau \n",SQR(ca/cb)  ,0.);

  fprintf(f," %12.4E   %12.4E   3    36     5    5 # higgs-b-b \n"    ,0.,SQR(tb*(MbH3/MbSM)));
  fprintf(f," %12.4E   %12.4E   3    36     6    6 # higgs-top-top \n",0.,SQR(1/tb)          );
  fprintf(f," %12.4E   %12.4E   3    36    15   15 # higgs-tau-tau \n",0.,SQR(tb)            );
     
  fclose(f);
  return 3;
}
