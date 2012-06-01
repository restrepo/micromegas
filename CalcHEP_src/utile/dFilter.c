#include<stdio.h>
#include<math.h>
               
typedef struct { char version[100];
                      int nIn,nOut;
                      double cs;
                      int  csOk;
                      int nEvents;
                      struct {double P;
                              long num;
                              char strfun[40];
                             } inInf[2];
                            
                      struct { double mass; 
                               long   num;
                               char   name[10];
                             } pInf[100];  
               }eventHeader;
               
             
int readEventHeader(FILE *f,eventHeader * H)
{  char word[100];
   int i;
   for(;1==fscanf(f,"%s",word);)
   {
     if(strcmp(word,"#CalcHEP")==0){fscanf(f," %[^\n]",H->version);}
     if(strcmp(word,"#Type")==0){fscanf(f,"%d -> %d",&(H->nIn),
                                                     &(H->nOut));}
     if(strcmp(word,"#Initial_state")==0)
     { if(H->nIn==2)
       {  
         fscanf(f," P1_3=%lf P2_3=%lf",&(H->inInf[0].P), &(H->inInf[1].P));
         for(i=0;i<2;i++) 
         fscanf(f," StrFun%*c=\"%[^\"]\" %ld",H->inInf[i].strfun,
                                             &(H->inInf[i].num));
       }else  fscanf(f," P1_3=%lf", &(H->inInf[0].P));                                         
     }
     if(strcmp(word,"#PROCESS")==0)
     { for(i=0;i<H->nIn;i++) fscanf(f," %d(%[^)])",&(H->pInf[i].num),
                                                   H->pInf[i].name);
       fscanf(f," ->");
       for(i=H->nIn;i<H->nIn+H->nOut;i++) fscanf(f," %ld(%[^)])",&(H->pInf[i].num),
                                                              H->pInf[i].name);             
     }
     if(strcmp(word,"#MASSES")==0)
        for(i=0;i<H->nIn+H->nOut;i++) fscanf(f," %lf",&(H->pInf[i].mass));
                                                        
     if(strcmp(word,"#Cross_section(Width)")==0)
        {if(fscanf(f," %lf",&(H->cs))!=1) H->csOk=0; else H->csOk=1;}
     if(strcmp(word,"#Number_of_events")==0){fscanf(f," %d",&(H->nEvents));}
     if(strcmp(word,"#Events")==0) { for(; fgetc(f)!='\n';); return 1;}
   }
   return 0;
}


int writeEventHeader(FILE *f,eventHeader * H)
{  char word[100];
   int i;
 
   fprintf(f,"#CalcHEP %s\n",H->version);
   fprintf(f,"#Type %d -> %d\n",H->nIn,H->nOut);
   
   fprintf(f,"#Initial_state\n");   
   for(i=0;i<H->nIn;i++) fprintf(f,"  P%d_3=%.6E",i+1,H->inInf[i].P);
   fprintf(f,"\n");
   if(H->nIn==2)
     for(i=0;i<2;i++) fprintf(f,"StrFun%d=\"%s\" %d\n",i+1, H->inInf[i].strfun,
                                                            H->inInf[i].num); 

   fprintf(f,"#PROCESS  ");
   for(i=0;i<H->nIn;i++) fprintf(f," %d(%s)",H->pInf[i].num,H->pInf[i].name);
   fprintf(f," -> ");
   for(i=H->nIn;i<H->nIn+H->nOut;i++) fprintf(f," %d(%s)",H->pInf[i].num,
                                                           H->pInf[i].name);
   fprintf(f,"\n");

   fprintf(f,"#MASSES");
   for(i=0;i<H->nIn+H->nOut;i++) fprintf(f," %.10E",H->pInf[i].mass);
   fprintf(f,"\n");

   fprintf(f,"#Cross_section(Width) ");
   if(H->csOk) fprintf(f,"%.6E\n",H->cs); else fprintf(f,"Unknown\n"); 
   fprintf(f,"#Number_of_events   %d\n",H->nEvents);
   fprintf(f,"#Events\n");

   return 1;
}

int readEvent(FILE*f, int *Nmom, double * mom, int * clr1, int * clr2, double *Q, int * w)
{ int n=0,c1,c2;
  if(fscanf(f,"%d",w)!=1) return 0;
  for(;1==fscanf(f," %lf",mom+n);n++);
  *Nmom=n;
  fscanf(f,"| %lf",Q);
  for(n=0;2==fscanf(f," (%d %d)",clr1+n,clr2+n);n++);
  clr1[n]=0; clr2[n]=0;
  return 1;
}

void writeEvent(FILE*f, int Nmom, double * mom, int * clr1, int * clr2, double Q,int w)
{ int n=0,c1,c2;
  fprintf(f,"%8d ",w);
  for(n=0;n<Nmom;n++) fprintf(f," %17.10E", mom[n]);
  fprintf(f,"| %10.3E   ",Q);
  for(n=0;clr1[n];n++)  fprintf(f,"(%d %d)",clr1[n],clr2[n]);
  fprintf(f,"\n");
}


static void  boost(double * n, double *p)
{

   double M=p[0]; /* particle mass */
   double shY=n[0]; /* h-sine of rapidity */
   double chY=sqrt(1+ shY*shY);
   double f= sqrt(M*M+p[1]*p[1]+p[2]*p[2]+p[3]*p[3])*shY 
                   + (p[1]*n[1]+p[2]*n[2]+p[3]*n[3])*(chY-1);

   int i;
   for(i=1;i<=3;i++) p[i]+=n[i]*f;
} 
  

static void findBoost(double * p, double *n)
{
  double p2=p[1]*p[1]+p[2]*p[2]+p[3]*p[3];
  double M=p[0];  
  int i;

  if(p2 == 0) {for(i=0;i<=3;i++) n[i]=0; return;}  
  p2=sqrt(p2);
  n[0]=p2/M;
  for(i=1;i<=3;i++) n[i]=p[i]/p2;
}


int main(int argc,char** argv)
{ FILE*f;
  eventHeader inProc,Decay,outProc;
  double brRat;
  int err,mult,Shift;
    
  int n,nMom;
  double mom[100], Q;
  int colr1[100],colr2[100],w;
  
  if(argc != 3) 
  { fprintf(stderr,"this program needs 2 arguments:\n"
    " a) the branching ration;\n b) file with decay events\n");
    return 1;
  }
  
  if(1!=sscanf(argv[1],"%lf",&brRat))
  { fprintf(stderr,"The first argument is not a number\n"); return 2; }
        
  f=fopen(argv[2],"r");
  if(!f) { fprintf(stderr,"Can not open file %e\n",argv[2]); return 3;} 

  if(!readEventHeader(f,&Decay))
  { fprintf(stderr,"Can not read header of Decay file\n");   return 4;} 
  
  if(Decay.nIn!=1) 
  { fprintf(stderr,"%s is not a file of decays\n", argv[2]);  return 5;}  
   
  if(!readEventHeader(stdin,&inProc)) 
  { fprintf(stderr," no event file in stdin \n");  return 6;}
  
  for(mult=0,n=inProc.nIn;n<inProc.nIn+inProc.nOut;n++)
  { if(strcmp(Decay.pInf[0].name, inProc.pInf[n].name)==0) mult++;}
  
  if(mult==0) 
  { fprintf(stderr,"Warning: no paricle to disintegrate in the incoming list!\n");
    writeEventHeader(stdout,&inProc);
    while(readEvent(stdin, &nMom,mom, colr1, colr2, &Q, &w))
          writeEvent(stdout, nMom,mom, colr1, colr2, Q, w);
    return 0;      
  }

  if(inProc.nEvents*mult>Decay.nEvents) 
  { fprintf(stderr,"Not enough events in Decay file \n");  return 7;} 
  
  outProc=inProc;
  if(outProc.csOk) outProc.cs*=pow(brRat,mult);  
  outProc.nOut+=mult*(Decay.nOut-1);     
  for(Shift=0,n=inProc.nIn;n<inProc.nIn+inProc.nOut;n++)
    if(strcmp(Decay.pInf[0].name, inProc.pInf[n].name)==0)
    { int k; 
      for(k=1;k<=Decay.nOut;k++,Shift++)outProc.pInf[n+Shift]=Decay.pInf[k];
      Shift--;    
    } else outProc.pInf[n+Shift]=inProc.pInf[n];

  writeEventHeader(stdout,&outProc);

  while(readEvent(stdin, &nMom,mom, colr1, colr2, &Q, &w))
  { double momOut[300];
    double nMomOut=nMom+(Decay.nOut -1)*3*mult;
    int m0=0,i,j;
  
    if(inProc.nIn==2){ m0=2; momOut[0]=mom[0]; momOut[1]=mom[1];} 
    
    for(Shift=0,n=0;n<inProc.nOut;n++)
    if(strcmp(Decay.pInf[0].name, inProc.pInf[n+inProc.nIn].name)==0)
    {  int k;
       int nMom_, colr1_[100], colr2_[100], w_;
       double mom_[300],Q_, b[4],p[4]; 
       readEvent(f, &nMom_,mom_, colr1_, colr2_, &Q_, &w_);
       w*=w_;

       for(i=0;colr1[i];i++)
       { if(colr1[i]>n+Shift+inProc.nIn+1)colr1[i]+=Decay.nOut-1;
         if(colr2[i]>n+Shift+inProc.nIn+1)colr2[i]+=Decay.nOut-1;
       }  
       
       for(i=0;colr1_[i];i++) 
       { if(colr1_[i]==1) 
         { for(j=0;colr2[j];j++) if(colr2[j] ==n+Shift+inProc.nIn+1)
           { colr2[j]=colr2_[i]+n+Shift+inProc.nIn-1;break;} 
         } else if(colr2_[i]==1)
         { for(j=0;colr1[j];j++) if(colr1[j] ==n+Shift+inProc.nIn+1)
           { colr1[j]=colr1_[i]+n+Shift+inProc.nIn-1;break;} 
         } else
         { for(j=0;colr1[j];j++);
           colr1[j]=colr1_[i]+n+Shift+inProc.nIn-1;
           colr2[j]=colr2_[i]+n+Shift+inProc.nIn-1;
           j++;
           colr1[j]=0;
           colr2[j]=0;
         }
       }
       
       p[0]=inProc.pInf[n+inProc.nIn].mass;
       for(i=0;i<3;i++) p[i+1]=mom[m0+3*n+i];  
       findBoost(p,b);		     
       
       for(k=1;k<=Decay.nOut;k++,Shift++)
       {  
          p[0]=Decay.pInf[k].mass;
          for(i=0;i<3;i++) p[i+1]=mom_[3*(k-1) +i]; 
          boost(b,p); 		  
 
	  for(i=0;i<3;i++)momOut[m0 +(n+Shift)*3 +i]=p[1+i]; 
	  
       }      
       Shift--;
    }else for(i=0;i<3;i++)momOut[m0 +(n+Shift)*3 +i]=mom[m0 +n*3 +i];
    writeEvent(stdout, nMomOut,momOut, colr1, colr2, Q, w);   
  }
  fclose(f);           
  return 0;  
}
