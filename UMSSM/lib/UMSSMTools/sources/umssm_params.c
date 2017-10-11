/* UMSSM parameters for UMSSMTools */
#include "../../pmodel.h"
#include "../../../../include/micromegas.h"
#include "../../../../CalcHEP_src/include/VandP.h"


int assignValFunc(char * name, double val)
{
  double * a=varAddress(name);
  if(a && a<=varValues+nModelVars+nModelFunc )  {*a=val; return 0;} else {printf(" %s not found\n", name);return 1;}
}

char UPARFnames[304][10];
//if additional higgses->sfermions functions are considered : char UPARFnames[364][10];
int PDG_LSP;


int read_prmU(int r)
{
  double val;
  char name[20];
  int n,err=0;
  FILE * f=fopen("UMSSMTools.par","r");
  if(f==NULL) {puts("Error : UMSSMTools.par not found.\n"); return -1;}
  for(n=1;;n++)
  { if(fscanf(f,"%s",name)!=1) { n=0; break;}
    if(name[0]=='#') { fscanf(f,"%*[^\n]"); continue;}
    if(fscanf(f,"%lf",&val)!=1) break;
    fscanf(f,"%*[^\n]");
    if(n==304) PDG_LSP=val;
    else
    { err=assignValFunc(name,val); /*if (n%2!=0) {printf("\n");} printf("***** %s = %f   ",name,val);*/
      strcpy(UPARFnames[n-1],name);
//      printf("Name of UPARF(%d) : %s = %f\n",n-1,UPARFnames[n-1],val);
      if(err==1) break;
    }
  }
  fclose(f);
  return 0;
}
/******************************************************************************************	
*******************************************************************************************
******************************************************************************************/	

double UparC(int r)
{
if(findValW("SW")<=0)
{
if(read_prmU(r))
{
puts("Error : Can not read one (or several) parameter(s) from UMSSMTools.par.");
return 0;
}
}
if (r==303) return PDG_LSP;
else if(r>=0 && r<=302)return findValW(UPARFnames[r]);
//if additional higgses->sfermions functions are considered : else if(r>=0 && r<=363)return findValW(UPARFnames[r]);
else {printf("Error : UPARF(%d) does not exist",r);return -1.;}
}

double uparctof_(int *r){return UparC(*r);}
