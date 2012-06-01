#include"micromegas_aux.h"
#include"micromegas.h"

void printHiggs(FILE* f)
{
   int i;
   fprintf(f,"Higgs masses and widths\n");
   for(i=0;i<nModelParticles;i++)
   { if (ModelPrtcls[i].name[0]!='~' && ModelPrtcls[i].spin2==0 && ModelPrtcls[i].cdim==1)
     { int j;
       double width;
       txtList LL;
       int dim;
       for(j=0;j<nModelVars;j++)if(strcmp(ModelPrtcls[i].width,varNames[j])==0)
       { width=varValues[j]; break;} 
       if(j==nModelVars) width=pWidth(ModelPrtcls[i].name, &LL,&dim);
        fprintf(f,"%5.5s  %7.2f %.2E\n",ModelPrtcls[i].name, 
            findValW(ModelPrtcls[i].mass), width);
     }
   }  
}
