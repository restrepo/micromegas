#include<math.h>
#include"pmodel.h"
#include"pmodel_aux.h"
#include"pmodel_f.h"
#include"../../sources/micromegas.h"

void FillVal(int mode) 
{ char name[10];
  int i,j,k;
  char* Zf[3]={"Zb","Zt","Zl"};
  char* Qmix[3]={"SBOTMIX","STOPMIX","STAUMIX"};
  char* massName[32]={"MH3","Mh","MHH","MHc", "MNE1", "MNE2", "MNE3", "MNE4",  "MC1",  "MC2",  "MSG", "MSne", "MSnm", "MSnl", "MSeL", "MSeR", "MSmL", "MSmR", "MSl1", "MSl2", "MSdL", "MSdR", "MSuL", "MSuR", "MSsL", "MSsR", "MScL", "MScR", "MSb1", "MSb2", "MSt1","MSt2"};
  int   massId[32]  ={ 36  ,  25,   35,   37,1000022,1000023,1000025,1000035,1000024,1000037,1000021,1000012,1000014,1000016,1000011,2000011,1000013,2000013,1000015,2000015,1000001,2000001,1000002,2000002,1000003,2000003,1000004,2000004,1000005,2000005,1000006,2000006};   
  double Q;
  if(mode==0)i=1; else i=0; 
  for(;i<32;i++) assignValW(massName[i],slhaVal("MASS",0.,1,massId[i]));
  Q=sqrt(fabs(findValW("MSt1")*findValW("MSt2")));
    
  for(i=1;i<=4;i++) for(j=1;j<=4;j++) 
  { sprintf(name,"Zn%d%d",i,j); assignValW(name,slhaVal("NMIX",Q,2,i,j));}
  
  for(i=1;i<=2;i++) for(j=1;j<=2;j++)
  { sprintf(name,"Zu%d%d",i,j);assignValW(name,slhaVal("UMIX",Q,2,i,j));
    sprintf(name,"Zv%d%d",i,j);assignValW(name,slhaVal("VMIX",Q,2,i,j));
  }
   
  for(k=0;k<3;k++)
  { double M[3];
    for(i=1;i<=2;i++) 
    {M[i]=slhaVal(Qmix[k],Q,2,1,i);
     sprintf(name,"%s1%d",Zf[k],i); assignValW(name,M[i]);
    }  
    M[2]*=-1;
    for(i=1;i<=2;i++) 
    {
     sprintf(name,"%s2%d",Zf[k],i); assignValW(name,M[3-i]);
    }  
  }   
  assignValW("alpha", slhaVal("ALPHA",Q,0));
  
  if(mode>0)
  { int MGok[3];
    
    assignValW("tb",   
    slhaValExists("MINPAR",1,3)>0 ? slhaVal("MINPAR",Q,1,3):
    slhaVal("HMIX",Q,1,2)   );
     
    assignValW("mu", slhaVal("HMIX",Q,1,1));
    assignValW("Am", slhaValExists("Ae",2,2,2)>0 ? slhaVal("Ae",Q,2,2,2):slhaVal("Ae",Q,2,3,3));
    assignValW("Al", slhaVal("Ae",Q,2,3,3));
    assignValW("Ab", slhaVal("Ad",Q,2,3,3));
    assignValW("Ad", slhaValExists("Ad",2,2,2)>0 ? slhaVal("Ad",Q,2,2,2):slhaVal("Ad",Q,2,3,3));
    assignValW("At", slhaVal("Au",Q,2,3,3));
    assignValW("Au", slhaValExists("Au",2,2,2)>0 ? slhaVal("Au",Q,2,2,2):slhaVal("Au",Q,2,3,3));

    for(i=1;i<=2;i++) 
    { MGok[i]=(slhaValExists("MSOFT",1,i)>0);
      if(MGok[i]){ sprintf(name,"MG%d",i);assignValW(name,slhaVal("MSOFT",Q,1,i));}
    } 
    if(!(MGok[1] && MGok[2])) 
    { double mg1,mg2;
      CheckNCsector(NULL, &mg1, &mg2,NULL,NULL,NULL,NULL);
      assignValW("MG1",mg1);
      assignValW("MG2",mg2);
    }
  }
  
  if(mode==2)
  {
    assignValW("alfSMZ",slhaVal("SMINPUTS",Q,1,3) );
    assignValW("MbMb",  slhaVal("SMINPUTS",Q,1,5) );
    assignValW("Mtp",   slhaVal("SMINPUTS",Q,1,6) );
    assignValW("Ml",    slhaVal("SMINPUTS",Q,1,7) );
  }
   
  assignValW("dMb",deltaMb());
  assignValW("dMs",deltaMs());
  assignValW("dMd",deltaMd());    
}

