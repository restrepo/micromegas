#include"../../sources/micromegas.h"
#include"../../sources/micromegas_aux.h"
#include"pmodel.h"

void smodels(double Pcm, int nf,double csMinFb, char*fileName,int wrt) 
{ 
  int SMP[16]={1,2,3,4,5,6, 11,12,13,14,15,16, 21,22,23,24};
  int i,j;
  FILE*f=fopen(fileName,"w");
  int np=0;
  char**plist=NULL;
    
  fprintf(f,"BLOCK MASS\n");
  for(i=0;i<nModelParticles;i++) if(pMass(ModelPrtcls[i].name) <Pcm)
  { 
    for(j=0;j<16;j++) if(abs(ModelPrtcls[i].NPDG)==SMP[j]) break; 
    if(j==16 )
    { 
       np++; 
       plist=realloc(plist,np*sizeof(char*));
       plist[np-1]=ModelPrtcls[i].name;
       if(strcmp(ModelPrtcls[i].name,ModelPrtcls[i].aname))
       { np++;
         plist=realloc(plist,np*sizeof(char*));
         plist[np-1]=ModelPrtcls[i].aname;
       }    
      fprintf(f,"  %d  %E  # %s  \n",ModelPrtcls[i].NPDG,findValW(ModelPrtcls[i].mass),ModelPrtcls[i].name);   
    }
  }
  fprintf(f,"\n");

  for(i=0;i<nModelParticles;i++) 
  {  for(j=0;j<16;j++)
    { 
       if(ModelPrtcls[i].NPDG==SMP[j]) break;
    }
    if(j==16) slhaDecayPrint(ModelPrtcls[i].name,1,f); 
  }

    for(i=0;i<np;i++) for(j=i;j<np;j++) if(pMass(plist[i])+pMass(plist[j])<Pcm)
    if(plist[i][0]=='~' && plist[j][0]=='~')
    {  int q31,q32,q3,c1,c2;

       qNumbers(plist[i], NULL, &q31,&c1);
       qNumbers(plist[j], NULL, &q32,&c2);
       q3=q31+q32;
       if(q3<0) { q3*=-1; if(abs(c1)==3) c1*=-1; if(abs(c2)==3)  c2*=-1;}
       if(c1>c2){ int c=c1; c1=c2;c2=c;}
       
       if (  (c2==1 || (c1==1 && c2==8) || (c1==-3 && c2==3) || (c1==8 && c2==8) ) 
        
                                       && (q3!=0 && q3 !=3) ) continue;
                                       
       if ( ((c1==-3 && c2== 3)||(c1== 1 && c2== 1)||
             (c1== 8 && c2== 8)||(c1== 1 && c2== 8))  && (q3!=0 && q3!=3) ) continue;                            
       if ( ((c1== 3 && c2== 8)||(c1== 1 && c2== 3))  && (q3!=2)          ) continue;
       if ( ((c1==-3 && c2== 8)||(c1==-3 && c2== 1))  && (q3!=1)          ) continue;
       if (  (c1== 3 && c2== 3)                       && (q3!=4 && q3!=1) ) continue;
       if (  (c1==-3 && c2==-3)                       && (q3!=2)          ) continue;
        
       {  double dcs;
          double Qf=0.5*(pMass(plist[i])+pMass(plist[j]));
          dcs=hCollider(Pcm,1,nf,Qf,Qf,plist[i],plist[j],0,wrt);
          if(dcs>csMinFb*0.001)
          {
            fprintf(f,"XSECTION  %E   2212  2212  2  %d  %d\n",2*Pcm, pNum(plist[i]),pNum(plist[j])); 
/*pb*/      fprintf(f,"0  0  0  0  0  0 %E micrOMEGAs 3.6\n\n", dcs);
          }
       }
    }    

  fclose(f);
  free(plist);
  
  f=fopen("particles.py","w");
  fprintf(f,"#!/usr/bin/env python\n");
  
  fprintf(f,"rOdd ={\n");
  for(np=0,i=0;i<nModelParticles;i++) if(ModelPrtcls[i].name[0]=='~'  && pMass(ModelPrtcls[i].name) <Pcm )
  {  
     { if(np) fprintf(f,",\n");
       fprintf(f, " %d : \"%s\"", ModelPrtcls[i].NPDG,ModelPrtcls[i].name);
       if(strcmp(ModelPrtcls[i].name,ModelPrtcls[i].aname))
       { fprintf(f,",\n");
         fprintf(f, " %d : \"%s\"", -ModelPrtcls[i].NPDG,ModelPrtcls[i].aname);
       }
       np++;
     }           
  }
  fprintf(f,"\n}\n");

  fprintf(f,"rEven ={\n");
  for(np=0,i=0;i<nModelParticles;i++) if(ModelPrtcls[i].name[0]!='~' && pMass(ModelPrtcls[i].name) <Pcm  )
  {  
     for(j=0;j<16;j++) if(abs(ModelPrtcls[i].NPDG)==SMP[j]) break;
     if(j==16 )
     { if(np) fprintf(f,",\n");
       fprintf(f, " %d : \"%s\"", ModelPrtcls[i].NPDG,ModelPrtcls[i].name);
       if(strcmp(ModelPrtcls[i].name,ModelPrtcls[i].aname))
       { fprintf(f,",\n");
         fprintf(f, " %d : \"%s\"", -ModelPrtcls[i].NPDG,ModelPrtcls[i].aname);
       }
       np++;
     }
  }

     fprintf(f,",\n"
"  23 : \"Z\",\n"
"  22 : \"photon\",\n"
"  24 : \"W+\",\n"
" -24 : \"W-\",\n"
"  16 : \"nu\",\n"
" -16 : \"nu\",\n"
"  15 : \"ta-\",\n"
" -15 : \"ta+\",\n"
"  14 : \"nu\",\n"
" -14 : \"nu\",\n"
"  13 : \"mu-\",\n"
" -13 : \"mu+\",\n"
"  12 : \"nu\",\n"
" -12 : \"nu\",\n"
"  11 : \"e-\",\n"
" -11 : \"e+\",\n"
"  5  : \"b\",\n"
" -5  : \"b\",\n"
"  6  : \"t+\",\n"
" -6  : \"t-\",\n"
"  1  : \"jet\",\n"
"  2  : \"jet\",\n"
"  3  : \"jet\",\n"
"  4  : \"jet\",\n"
"  21 : \"jet\",\n"
" -1  : \"jet\",\n"
" -2  : \"jet\",\n"
" -3  : \"jet\",\n"
" -4  : \"jet\""  );
  
  fprintf(f,"\n}\n");

fprintf(f,  
" ptcDic = {\"e\"  : [\"e+\",  \"e-\"],\n"
"          \"mu\" : [\"mu+\", \"mu-\"],\n"
"          \"ta\" : [\"ta+\", \"ta-\"],\n"
"          \"l+\" : [\"e+\",  \"mu+\"],\n"
"          \"l-\" : [\"e-\",  \"mu-\"],\n"
"          \"l\"  : [\"e-\",  \"mu-\", \"e+\", \"mu+\"],\n"
"          \"W\"  : [\"W+\",  \"W-\"],\n"
"          \"t\"  : [\"t+\",  \"t-\"],\n"
"          \"L+\" : [\"e+\",  \"mu+\", \"ta+\"],\n"
"          \"L-\" : [\"e-\",  \"mu-\", \"ta-\"],\n"
"          \"L\"  : [\"e+\",  \"mu+\", \"ta+\", \"e-\", \"mu-\", \"ta-\"]}\n"
);  
  
  
  

  fprintf(f,"qNumbers ={\n");
  for(np=0,i=0;i<nModelParticles;i++) if(pMass(ModelPrtcls[i].name) <Pcm  )
  {  
     for(j=0;j<16;j++) if(abs(ModelPrtcls[i].NPDG)==SMP[j]) break;
     if(j==16 )
     { if(np) fprintf(f,",\n");
       fprintf(f, " %d : [%d,%d,%d]", ModelPrtcls[i].NPDG, ModelPrtcls[i].spin2, ModelPrtcls[i].q3, ModelPrtcls[i].cdim);
       np++;
     }           
  }
  fprintf(f,"\n}\n");

  
  fclose(f); 
}
