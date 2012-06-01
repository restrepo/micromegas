#include"../../sources/micromegas.h"
#include"../../sources/micromegas_aux.h"
#include"pmodel.h"
#include"pmodel_f.h"

int readVarlHiggs(char * fname)
{
  char*vlist[24]={"EE","SW","alfSMZ","McMc","MbMb","Mtpole",
  "s12","s23","s13","MZ","f","sa","zero0","MH","kappa","kappal",
  "Mu","Md","Ms","wW","wZ","Mm","Ml","wtop"};

  return readVarSpecial(fname,24,vlist);
} 

int  readvarlhiggs_(char * f_name,int len)
{
  char c_name[100];
  fName2c(f_name,c_name,len);
  
  return readVarlHiggs(c_name);
}
