#include<math.h>
#include"../sources/micromegas.h"

#define a (0.52)

struct 
{ double A, Z, J; double (*S00)(double);double (*S01)(double);double (*S11)(double);}

sdInfo[20]={
{  19,Z_F  ,J_F19   ,S00F19    ,S01F19    ,S11F19   },
{  23,Z_Na ,J_Na23  ,S00Na23   ,S01Na23   ,S11Na23  },
{  29,Z_Si ,J_Si29  ,S00Si29A   ,S01Si29A   ,S11Si29A  },
{  27,Z_Al ,J_Al27  ,S00Al27   ,S01Al27   ,S11Al27  },
{  39,Z_K  ,J_K39   ,S00K39    ,S01K39    ,S11K39   },
{  73,Z_Ge ,J_Ge73  ,S00Ge73   ,S01Ge73   ,S11Ge73  },
{  93,Z_Nb ,J_Nb93  ,S00Nb93   ,S01Nb93   ,S11Nb93  },
{  125,Z_Te ,J_Te125 ,S00Te125  ,S01Te125  ,S11Te125 },
{  127,Z_I  ,J_I127  ,S00I127   ,S01I127   ,S11I127  },
{  129,Z_Xe, J_Xe129, S00Xe129  ,S01Xe129  ,S11Xe129 },
{  131,Z_Xe ,J_Xe131 ,S00Xe131  ,S01Xe131  ,S11Xe131 },
{  207,Z_Pb ,J_Pb207 ,S00Pb207  ,S01Pb207  ,S11Pb207 },
     
{  23,Z_Na ,J_Na23  ,S00Na23A   ,S01Na23A   ,S11Na23A  },
{  29,Z_Si ,J_Si29  ,S00Si29    ,S01Si29    ,S11Si29  },
{  73,Z_Ge ,J_Ge73  ,S00Ge73A   ,S01Ge73A   ,S11Ge73A  },
{  125,Z_Te ,J_Te125 ,S00Te125A  ,S01Te125A  ,S11Te125A },
{  127,Z_I  ,J_I127  ,S00I127A   ,S01I127A   ,S11I127A  },
{  129,Z_Xe, J_Xe129, S00Xe129A  ,S01Xe129A  ,S11Xe129A },
{  131,Z_Xe ,J_Xe131 ,S00Xe131A  ,S01Xe131A  ,S11Xe131A },
{  131,Z_Xe ,J_Xe131 ,S00Xe131B  ,S01Xe131B  ,S11Xe131B }
};


int main(void)
{ int i;
  double E0=15E-6;
  double Sp,Sn;
  
  double Ap=-1.14, An=1;  /* Z-induced interaction */
  
  double Aplus=Ap+An, Aminus=Ap-An;
  
  
 printf("plot [2.:6.0] \"-\", 1.7*x-0.28 -0.82* ( x-3.7 + sqrt((x-3.7)**2+0.2)) \n");
    
  for(i=0;i<11;i++)
  { double p=sqrt(2*sdInfo[i].A*0.939*E0)/0.197327;
    double J=sdInfo[i].J;
    double coeff;
    double R0,R;
    double FFth0,FFth1,FF;
 
    FFth0=sdInfo[i].S00(0.)*Aplus*Aplus+ Aminus*Aminus*sdInfo[i].S11(0.)+
    Aplus*Aminus*sdInfo[i].S01(0.);

    FFth1=sdInfo[i].S00(p)*Aplus*Aplus+ Aminus*Aminus*sdInfo[i].S11(p)+
    Aplus*Aminus*sdInfo[i].S01(p);
    
    R=-log(FFth1/FFth0)/(p*p);
    R0=pow(sdInfo[i].A,1./3.);
    printf("  %E %E \n", pow(sdInfo[i].A,1./3.),  2*sqrt(R)) ; 
  }
}

