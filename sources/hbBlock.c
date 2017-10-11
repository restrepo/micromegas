#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"

static int txt2plist(char*txt,double *br, int * ind)
{  
   if(0==sscanf(txt,"%lf",br)) return 0;
   char *chE, *chB=strstr(txt,"->");    
   if(!chB) return 0;
   chB+=2;
   int n;
   
   for(n=0, chE=chB; chE; n++, chE=strchr(chE+1,','))
   {  char buff[20];
         sscanf(chE+1,"%[^,]",buff);
         trim(buff); 
         ind[n]=pNum(buff);
   }
   return n;
}



int hbBlocksMO(char *fname, int *nHiggsCh)  
{ int i,j,k;

  FILE*f=fopen(fname,"w");

  double Mcp=1.478662E+00, Mbp=bPoleMass(), Mtp=tPoleMass();
  int VZdecay_save=VZdecay, VWdecay_save=VWdecay, changeVirtual=0;
  
  double coeff[10];

  char  *G =pdg2name(21);  if(!G) { printf("Gluon  is absent in the model");    return -1;} 
  char  *A =pdg2name(22);  if(!A) { printf("Photon is absent in the model");    return -1;}
  char  *Z = pdg2name(23); if(!Z) { printf("Z is absent in the model");         return -1;}
  char  *Wp= pdg2name(24); if(!Wp){ printf("W+ is absent in the model");        return -1;}
  char  *Wm=  antiParticle(Wp);
  double MW=pMass(Wp);
   
  lVert* WWA=getLagrVertex(Wp,Wm,A,NULL); if(!WWA) { printf(" %s %s %s is absent\n", Wp,Wm,A); return -2;}
  getNumCoeff(WWA,coeff);
//for(k=0;k<WWA->nTerms;k++) printf("%s ",WWA->SymbVert[k]);
//printf("\n");  
  for(k=0;k<WWA->nTerms;k++) if( strcmp(WWA->SymbVert[k],"p1.m1*m3.m2")==0) break;
  double EE=coeff[k];
//for(k=0;k<WWA->nTerms;k++) printf( "%e ", coeff[k]);
//printf("\n"); 
  
  lVert* WWZ=getLagrVertex(Wp,Wm,Z,NULL); if(!WWZ) { printf(" %s %s %s is absent\n", Wp,Wm,Z); return -2;}
  getNumCoeff(WWZ,coeff);
  for(k=0;k<WWZ->nTerms;k++) if( strcmp(WWZ->SymbVert[k],"p1.m1*m3.m2")==0) break;
  double SW=sin(atan(EE/coeff[k]));
  double vev= 2*MW*SW/EE;

  fprintf(f,"# SM basic parameters: MW=%.3E alphaEM^{-1}=%.3E SW^2=%.3E  VEV=%.3E \n", MW,4*M_PI/EE/EE,SW*SW,vev);

// Find Higgses 

   char ** Higgs=NULL;
   int nHiggs=0;
   char ** chHiggs=NULL;
   int nHch=0;
   
   for(i=0;i<nModelParticles;i++) if(ModelPrtcls[i].name[0]!='~' && ModelPrtcls[i].spin2==0 
                                   && ModelPrtcls[i].cdim==1)
   { 
       if(strcmp(ModelPrtcls[i].name,ModelPrtcls[i].aname)==0)
       {  Higgs=(char**)realloc(Higgs, (sizeof(char*))*(nHiggs+1));
          Higgs[nHiggs++]=ModelPrtcls[i].name;
       } else
       {  chHiggs=(char**)realloc(chHiggs, (sizeof(char*))*(nHch+1));
          chHiggs[nHch++]=ModelPrtcls[i].name;
       }   
   }
   fprintf(f," BLOCK MASS\n");
   for(i=0;i<nHiggs;i++) fprintf(f," %d  %E # %s\n", pNum(Higgs[i]), pMass(Higgs[i]),Higgs[i]);    
   for(i=0;i<nHch;i++)   fprintf(f," %d  %E # %s\n", pNum(chHiggs[i]), pMass(chHiggs[i]),chHiggs[i]);
 
   fprintf(f,"Block HiggsBoundsInputHiggsCouplingsFermions\n");
   int ferm[3]={5,6,15};
   for(i=0;i<nHiggs;i++) for(j=0;j<3;j++)
   {  char *fe=pdg2name(ferm[j]); if(!fe) continue;
      char *Fe=antiParticle(fe);  if(!Fe) continue; 
      double mf=pMass(fe); if(mf==0) continue;
      lVert *hff=getLagrVertex(Fe,fe,Higgs[i],NULL); 
      double c[2]={0,0};
      if(hff)
      {  getNumCoeff(hff,coeff); 
         for(k=0;k<hff->nTerms;k++)
         { if(strcmp(hff->SymbVert[k],"1")==0)    c[0]=coeff[k]*vev/mf; else 
           if(strcmp(hff->SymbVert[k],"G5*i")==0) c[1]=coeff[k]*vev/mf;
         }
      }       
      fprintf(f," %E  %E  3  %d  %d  %d  # %s %s %s\n", c[0]*c[0], c[1]*c[1], pNum(Higgs[i]), ferm[j], ferm[j], Higgs[i],Fe,fe); 
   }   
   
   fprintf(f,"Block HiggsBoundsInputHiggsCouplingsBosons\n");
   int vbm[2]={23,24};
   for(i=0;i<nHiggs;i++) for(j=0;j<2;j++)
   {  char *ve=pdg2name(vbm[j]); if(!ve) continue;
      char *Ve=antiParticle(ve); if(!Ve) continue; 
      lVert *hvv=getLagrVertex(Ve,ve,Higgs[i],NULL);
//printf("bosons: %s %s\n", ve,Ve);             
      double c=0;
      if(hvv)
      {  getNumCoeff(hvv,coeff); 
         for(k=0;k<hvv->nTerms;k++) if(strcmp(hvv->SymbVert[k],"m2.m1")==0) c=0.5*coeff[k]*vev/pMass(ve)/pMass(ve);  
      }    
      fprintf(f," %E  3  %d  %d  %d  # %s %s %s\n", c*c, pNum(Higgs[i]), vbm[j], vbm[j], Higgs[i],Ve,ve); 
   }

   for(i=0;i<nHiggs;i++) for(j=i+1;j<nHiggs;j++)
   {  char *ve=pdg2name(23); if(!ve) continue;
      lVert *hhv=getLagrVertex(Higgs[i],Higgs[j],ve,NULL);
      double c=0;
      if(hhv)
      {  getNumCoeff(hhv,coeff); 
         for(k=0;k<hhv->nTerms;k++) if(strcmp(hhv->SymbVert[k],"p1.m3*i")==0) c=coeff[k]*vev/pMass(ve); 
      }    
      fprintf(f," %E  3  %d  %d  %d  # %s %s %s\n", c*c, pNum(Higgs[i]),pNum(Higgs[j]) , 23, Higgs[i],Higgs[j] ,ve); 
   }

   for(i=0;i<nHiggs;i++) for(j=0;j<nHch;j++)
   {  int charge3;
      qNumbers(chHiggs[j],NULL, &charge3, NULL); 
      char *ve=NULL;
      if(charge3==3) ve=pdg2name(-24);
      else if(charge3==-3) ve=pdg2name(24);
      if(!ve) continue;
      lVert *hhv=getLagrVertex(Higgs[i],chHiggs[j],ve,NULL);
      double c[2]={0,0};
      if(hhv)
      {  getNumCoeff(hhv,coeff); 
         for(k=0;k<hhv->nTerms;k++) if(strcmp(hhv->SymbVert[k],"p1.m3*i")==0) c[0]=coeff[k]*vev/pMass(ve);
                               else if(strcmp(hhv->SymbVert[k],"p1.m3")==0) c[1]=coeff[k]*vev/pMass(ve); 
      }    
      fprintf(f," %E  3  %d  %d  %d  # %s %s %s\n", c[0]*c[0]+c[1]*c[1], pNum(Higgs[i]),pNum(chHiggs[j]) ,pNum(ve), Higgs[i],chHiggs[j] ,ve); 
   }
   
   
// Gluon  and photon decays  
   double complex * ffE=(double complex *)malloc(nHiggs*sizeof(double complex));  // FmunuFmunu couplings      for photons
   double complex * faE=(double complex *)malloc(nHiggs*sizeof(double complex));  // FmunuEps()Fmunu couplings for photons
   double complex * ffC=(double complex *)malloc(nHiggs*sizeof(double complex));  // FmunuFmunu couplings      for gluons
   double complex * faC=(double complex *)malloc(nHiggs*sizeof(double complex));  // FmunuEps()Fmunu couplings for gluons
   double *lamQGG      =(double *)malloc(nHiggs*sizeof(double));
   double *lamQAA      =(double *)malloc(nHiggs*sizeof(double));          
   for(i=0;i< nHiggs;i++)
   {  ffE[i]=faE[i]=ffC[i]=faC[i]=0; 
      txtList L = makeDecayList(Higgs[i],2), l=L;
      double mH=pMass(Higgs[i]);
      double a=alphaQCD(mH)/M_PI;
      for(l=L;l;l=l->next)
      { char Xp[10], Xm[10];
        sscanf(strstr(l->txt,"->")+2, "%[^,],%s",Xp,Xm);
        trim(Xp); trim(Xm); 
        if(strcmp(Xp,antiParticle(Xm))==0)
        {  
          int pdg,spin2,charge3,cdim;
          double mX=pMass(Xp);
          double dffE=0,dfaE=0,dffC=0,dfaC=0;
          pdg=qNumbers(Xp, &spin2, &charge3, &cdim);
          if((charge3 !=0 || cdim!=1) && mX>0.5)
          {  double mXp; // pole mass
             switch(abs(pdg))
             { case 4: mXp=Mcp;    break;
               case 5: mXp=Mbp;      break;
               case 6: mXp=Mtp; break;
               default:mXp=mX;
             }             
             double mN= (spin2&1)?  mX : mX*mX;   
             lVert *xxh=getLagrVertex(Xm,Xp,Higgs[i],NULL);
             if(!xxh) continue;
             getNumCoeff(xxh,coeff);
             for(k=0;k<xxh->nTerms;k++)  
             { int addFF=0,addFA=0;
               switch(spin2)
               {
                 case 0:  if(strcmp(xxh->SymbVert[k],"1")==0) addFF=1;    break;
                 case 1:  if(strcmp(xxh->SymbVert[k],"1")==0) addFF=1; else
                          if(strcmp(xxh->SymbVert[k],"G5*i")==0) addFA=1; break;
                 case 2:  if(strcmp(xxh->SymbVert[k],"m2.m1")==0) addFF=1;break; 
               }
               if(addFF)
               { if(cdim!=1) ffC[i]+=hGGeven(mH,a,1,spin2,cdim,mXp,coeff[k]/mN); 
                 if(charge3) ffE[i]+=hAAeven(mH,a,1,spin2,cdim,mXp,coeff[k]/mN)*charge3*charge3/9.;
               }
               if(addFA) 
               { if(cdim!=1) faC[i]+=0.5*hGGodd(mH,a,1,spin2,cdim,mXp,coeff[k]/mN); 
                 if(charge3) faE[i]+=0.5*hAAodd(mH,a,1,spin2,cdim,mXp,coeff[k]/mN)*charge3*charge3/9.; 
               } 
             }   
          }                     
        }    
      }
      
      double  LGGSM=lGGhSM(mH,alphaQCD(mH)/M_PI, Mcp,Mbp,Mtp,vev);
      double  LAASM=lAAhSM(mH,alphaQCD(mH)/M_PI, Mcp,Mbp,Mtp,vev);
              lamQGG[i]=(ffC[i]*conj(ffC[i]) +4*faC[i]*conj(faC[i]));
              lamQAA[i]=(ffE[i]*conj(ffE[i]) +4*faE[i]*conj(faE[i]));      
      fprintf(f," %E  3  %d  %d  %d  # %s %s %s\n", lamQGG[i]/LGGSM/LGGSM, pNum(Higgs[i]) ,pNum(G), pNum(G), Higgs[i] ,G ,G);
      fprintf(f," %E  3  %d  %d  %d  # %s %s %s\n", lamQAA[i]/LAASM/LAASM, pNum(Higgs[i]) ,pNum(A), pNum(A), Higgs[i] ,A ,A);
      cleanTxtList(L); 
   }
   free(ffE); free(faE), free(ffC); free(faC);    
  
   if(VZdecay==0 || VWdecay==0) 
   { VZdecay=1; VWdecay=1; cleanDecayTable(); changeVirtual=1;}

   char  *t =pdg2name(6);
   if(t) slhaDecayPrint(t, 0, f);

   for(i=0;i< nHiggs;i++)
   {  txtList L,l;
      double  BrAA=0, BrGG=0;
      double width,widthP;
      width=pWidth(Higgs[i],&L);
      for(l=L;l;l=l->next)
      {  double br;
         int id[10],nd;
         nd=txt2plist(l->txt,&br,id); 
         if(nd==2)
         { if(id[0]==21 && id[1]==21) BrGG=br;
           if(id[0]==22 && id[1]==22) BrAA=br;
         }   
      }
//BrGG=0;  BrAA=0;
      if(BrGG>0 && BrAA>0) slhaDecayPrint(Higgs[i], 0, f); else
      { double wGG,wAA,brGG_=0,brAA_=0;
        double mH=pMass(Higgs[i]);
        double a=alphaQCD(mH)/M_PI;
        double Rqcd=(1+a*(149./12.+a*(68.6482-a*212.447))); 
        wGG=2*Rqcd*pow(mH,3)*lamQGG[i]/M_PI;
        wAA=0.25*pow(mH,3)*lamQAA[i]/M_PI;        
        double K=1;
        if(BrGG<=0 && wGG>0) K+=wGG/width;
        if(BrAA<=0 && wAA>0) K+=wAA/width; 
        fprintf(f,"DECAY %d %E # %s\n", pNum(Higgs[i]), K*width, Higgs[i]);
        if(BrGG<=0 && wGG>0) fprintf(f," %E 2 %4d %4d # %s %s\n", wGG/(K*width),pNum(G),pNum(G),G,G);
        if(BrAA<=0 && wAA>0) fprintf(f," %E 2 %4d %4d # %s %s\n", wAA/(K*width),pNum(A),pNum(A),A,A); 
         
        for(l=L;l;l=l->next)
        { double br;
          int id[10],nd;
          nd=txt2plist(l->txt,&br,id);
          fprintf(f," %e %d ", br/K, nd);
          for(int j=0;j<nd;j++) fprintf(f," %4d", id[j]);
          fprintf(f," # %s\n", strstr(l->txt,"->")+2);
        }                                                                      
      }
   }
   
   for(i=0;i<nHch;i++) slhaDecayPrint(chHiggs[i], 0, f);
      
   free(Higgs); free(chHiggs); free(lamQGG); free(lamQAA);
   fclose(f);   
   if(nHiggsCh) *nHiggsCh=nHch;
   if(changeVirtual)
   { VZdecay=VZdecay_save; VWdecay=VWdecay_save; cleanDecayTable();}

   return nHiggs;
}  

int  hbblocksmo_(char *fname, int * nHch,int len)
{ 
  char * cname=malloc(len+2);
  int nHiggs;  
  fName2c(fname,cname,len);    
  nHiggs= hbBlocksMO(cname,nHch);      
  free(cname);        
  return nHiggs;          
}
  