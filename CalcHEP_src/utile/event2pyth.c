#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include <dirent.h>

extern double id2mass_(int*);

static double totCS=0, maxLum=0;

typedef struct eventfile_info
{ 
  struct eventfile_info* next;
  FILE *F;  
  char * fileName;
  int Nin,Nout;
  double cs;
  double inMom[2];
  int inPID[2];
  int PIDs[20];  
  double pmass[20];
  long nEvents;
  long cEvent;
  long FirstEventPos;  
} eventfile_info;

eventfile_info * All=NULL;

typedef struct decay_info
{ 
  struct decay_info* next;
  int ID;
  double  width;
  eventfile_info * List;
 } decay_info;

decay_info * Decays=NULL;


#define MAXNUP 500
extern struct  hepeup_
{    int     NUP,IDPRUP;
     double  XWGTUP,SCALUP,AQEDUP,AQCDUP;
     int     IDUP[MAXNUP],ISTUP[MAXNUP],MOTHUP[MAXNUP][2],ICOLUP[MAXNUP][2];
     double  PUP[MAXNUP][5], VTIMUP[MAXNUP],SPINUP[MAXNUP];
} hepeup_;

#define MAXPUP 100
extern struct heprup_
{  int IDBMUP[2];
   double  EBMUP[2];
   int PDFGUP[2],PDFSUP[2];
   int  IDWTUP,NPRUP;
   double XSECUP[MAXPUP],XERRUP[MAXPUP],XMAXUP[MAXPUP];
   int LPRUP[MAXPUP];
} heprup_;


void fName2c(char*f_name,char*c_name,int len)
{ int i; for(i=len-1;i>=0 &&f_name[i]==' ';i--);
  c_name[i+1]=0;
  for(;i>=0;i--) c_name[i]=f_name[i];
}


static eventfile_info * initEventFile(char* fname)
{
   char buff[200];
   eventfile_info * Finfo;
   FILE*F;
   int i,n;
   int ntot;
   
   F=fopen(fname,"r");
   if(F==NULL) return NULL;
   
   Finfo=(eventfile_info *)malloc(sizeof(eventfile_info)); 

   Finfo->F=F;
        
   for(fscanf(F,"%s",buff);!feof(F); )
   {
      if(strcmp(buff,"#Type")==0) 
      { fscanf(F,"%d -> %d", &Finfo->Nin, &Finfo->Nout);
        ntot=Finfo->Nin+Finfo->Nout;        
      } else if(strcmp(buff,"#Initial_state")==0)
      {
        for(i=0;i<2;i++) { Finfo->inPID[i]=0;Finfo->inMom[i]=0;}
        if(Finfo->Nin==2)
        { 
          for(i=0;i<2;i++){fscanf(F," P%d_3=",&n); fscanf(F,"%lf",Finfo->inMom+n-1);}
          for(i=0;i<2;i++) 
          { fscanf(F," StrFun%d=",&n); 
            if(fscanf(F,"\"%*[^\"]\" %d",Finfo->inPID+n-1)==0)Finfo->inPID[n-1]=0;
          }
        }  
      }
      else if(strcmp(buff,"#PROCESS")==0)
      {  int np;
         for(np=0;np<ntot;np++)
         { fscanf(F,"%d",&Finfo->PIDs[np]);
           fscanf(F,"%*s");
           if(np==Finfo->Nin-1)fscanf(F," -> ");
         }  
      }  else if(strcmp(buff,"#MASSES")==0)
      { int i;
         for(i=0;i<ntot;i++) fscanf(F," %lf",&Finfo->pmass[i]);
      } else if(strcmp(buff,"#Cross_section(Width)")==0)
      {  if( fscanf(F," %lf",&Finfo->cs)==0) Finfo->cs=0;
      } else if(strcmp(buff,"#Number_of_events")==0)
      {  fscanf(F," %ld",&Finfo->nEvents);      
      } else if(strcmp(buff,"#Events")==0)
      {  fscanf(F,"%*[^\n]%*c");
         for(i=0;i<Finfo->Nin;i++) if(Finfo->inPID[i]==0) Finfo->inPID[i]=Finfo->PIDs[i];

         if(Finfo->cs == 0) 
         { 
            fprintf(stderr,"Error: zero cross section/width in %s\n",
            Finfo->fileName);
            exit(5);
         } 
         if(Finfo->cs<0) Finfo->cs*=-1; 
         Finfo->cEvent=1;
         Finfo->FirstEventPos=ftell(F);
         Finfo->fileName=malloc(strlen(fname)+1);
         strcpy(Finfo->fileName,fname);
         return Finfo;
      }
      fscanf(F,"%s",buff);
   }
   return NULL; 
}

static eventfile_info* secondCopy(eventfile_info* new, eventfile_info* old)
{ int i;
  for(;old; old=old->next) if(new->Nout==old->Nout)
  { 
     for(i=0;i<new->Nin+new->Nout;i++) if(new->PIDs[i]!=old->PIDs[i]) break;
     if(i>=new->Nin+new->Nout) return  old;
  }
  return NULL;
}

int scandir_( char * dirname, int len)
{
  DIR *dirPtr;
  struct dirent * dp;
  char dirName[200];
  eventfile_info *Finfo, *Copy;
  int count=0;
  int i;
   
  
  fName2c(dirname,dirName,len);
      
  dirPtr=opendir(dirName);
  if(!dirPtr)
  { fprintf(stderr,"ERROR:%s - no such directory\n",dirName);
    exit(1);
  }  
  while((dp=readdir(dirPtr)))
  { char fullName[300];
    int N;
     
    if( sscanf(dp->d_name,"events_%d",&N)==1)
    { sprintf(fullName,"%s/%s", dirName,dp->d_name);  
      Finfo=initEventFile(fullName);
      if(!Finfo) continue;
      if(Finfo->Nin==2) 
      { double lum=Finfo->nEvents/Finfo->cs;

        if(All)
        { for(i=0;i<2;i++) 
          if( All->inMom[i]!=Finfo->inMom[i] ||
              All->inPID[i]!=Finfo->inPID[i]   ) 
          { fprintf(stderr, "Error: different incoming states in\n"
             "%s \n and  %s\n", All->fileName, Finfo->fileName);
            exit(1);
          }
          Copy= secondCopy(Finfo,All); 
          if(Copy) 
          { fprintf(stderr,"Error: Indentical processes in files\n%s\n%s\n",
            Copy->fileName, Finfo->fileName);
            exit(2);
          }
          if(maxLum > lum) maxLum=lum;          
        }else  maxLum=lum; 
          
        Finfo->next=All;
        All=Finfo;
        totCS+=All->cs;      
        count ++;
      }
      if(Finfo->Nin==1) 
      { 
        decay_info * D; 
        for(D=Decays;D;D=D->next)if(D->ID==Finfo->PIDs[0])
        { 
          Copy=secondCopy(Finfo,D->List);
          if(Copy) 
          { fprintf(stderr,"Error: Indentical processes in files\n%s\n%s\n",
             Copy->fileName, Finfo->fileName);
             exit(2);
          }  
          Finfo->next=D->List;
          D->List=Finfo;
          D->width+=Finfo->cs;
          break;
        }  
        if(!D)  
        { D=malloc(sizeof(decay_info));
          D->next=Decays;
          D->ID=Finfo->PIDs[0]; 
          D->width=Finfo->cs;
          D->List=Finfo;
          Decays=D;
          count++;
        }   
      }
    } 
  }  
  closedir(dirPtr);
  if(All)
  {  
     heprup_.XSECUP[0] = totCS;
     for(i=0;i<2;i++)
     { double m,p;
       heprup_.IDBMUP[i]=All->inPID[i];
       if(All->inPID[i]!=All->PIDs[i]) m=id2mass_(heprup_.IDBMUP+i);
                   else                m=All->pmass[i];
       p=All->inMom[i];          
       heprup_.EBMUP[i]=sqrt(m*m+p*p);           
     }  
  }
  return count;
}

void eventstat_(double * cs, int * nevents)
{ *cs=totCS;  *nevents=totCS*maxLum; }

static void readEvent(eventfile_info *Finfo, int *Nmom, double * mom, int * clr1, int * clr2, double *Q, int * w)
{ int n=0;
  if(Finfo->cEvent >= Finfo->nEvents)
  { Finfo->cEvent=1;
    fseek(Finfo->F, Finfo->FirstEventPos, SEEK_SET);  
    fprintf(stderr,"Warning: File %s finishes. Return to the beginning.\n",
    Finfo->fileName);
  } else Finfo->cEvent++; 

  fscanf(Finfo->F,"%d",w);
  for(;1==fscanf(Finfo->F," %lf",mom+n);n++);
  *Nmom=n;
  fscanf(Finfo->F,"| %lf",Q);
  for(n=0;2==fscanf(Finfo->F," (%d %d)",clr1+n,clr2+n);n++);
  clr1[n]=0; clr2[n]=0;
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


void upevnt_(void)
{
   int i,I,j,N, CC,Nmom,w, clr1[10],clr2[10];
   double cs, mom[32], Q;
   eventfile_info * Finfo;
   
   if(!All) 
   { fprintf(stderr,"Error:No event file is found\n");
     exit(2);
   }

   for(Finfo=All,cs=totCS*drand48(); cs>Finfo->cs; )
   { cs-=Finfo->cs;
     if(Finfo->next) Finfo=Finfo->next;
   }
   hepeup_.IDPRUP=1;
   hepeup_.XWGTUP=1.;
   N=Finfo->Nin+Finfo->Nout;
   hepeup_.NUP=N;
   
   for(I=0;I<N;I++)
   {  hepeup_.PUP[I][4]=Finfo->pmass[I];
      hepeup_.IDUP[I] =Finfo->PIDs[I];
      hepeup_.SPINUP[I]=9;
      hepeup_.VTIMUP[I]=0;
      for(j=0;j<2;j++) hepeup_.ICOLUP[I][j]=0;
   }
   readEvent(Finfo,&Nmom, mom, clr1, clr2, &Q, &w);
   hepeup_.SCALUP=Q;
   
   for(I=0;I<2;I++)
   { for(j=0;j<2;j++) {hepeup_.PUP[I][j]=0;hepeup_.MOTHUP[I][j]=0;}
     hepeup_.PUP[I][2]=mom[I];
     hepeup_.ISTUP[I]=-1;
   } 

   for(I=2;I<N;I++) 
   {  for(j=0;j<2;j++) hepeup_.MOTHUP[I][j]=j+1;
      for(j=0;j<3;j++) hepeup_.PUP[I][j]=mom[2+3*(I-2)+j];
      hepeup_.ISTUP[I]=1;    
   } 

   for(i=0,CC=500;clr1[i];i++,CC++)
   { int k1=clr1[i]-1,k2=clr2[i]-1; 
     if(k1<2) hepeup_.ICOLUP[k1][0]=CC; else hepeup_.ICOLUP[k1][1]=CC;
     if(k2<2) hepeup_.ICOLUP[k2][1]=CC; else hepeup_.ICOLUP[k2][0]=CC;
   }
   
   for(I=2;I<hepeup_.NUP;I++) if(hepeup_.ISTUP[I]==1)
   { decay_info *D;
     for(D=Decays;D;D=D->next)
     { printf("D=%p\n",D); if(D->ID==hepeup_.IDUP[I]) 
     { double width; 
       eventfile_info *Finfo=D->List;
       double p[4],b[4];
       int k;
       int shift=hepeup_.NUP-1;

       for(width=D->width*drand48();  Finfo->cs <= width; )
       { printf("width=%e,  Finfo->cs= %e\n",width,Finfo->cs );
         width-=Finfo->cs;   
         if(Finfo->next) Finfo=Finfo->next;
       }
       
       readEvent(Finfo,&Nmom, mom, clr1, clr2, &Q, &w);
       hepeup_.ISTUP[I]=2;
       p[0]=hepeup_.PUP[I][4]; 
       for(k=0;k<3;k++) p[k+1]=hepeup_.PUP[I][k]; 
       findBoost(p,b);
       for(k=1;k<=Finfo->Nout;k++)
       {  
         p[0]=Finfo->pmass[k];
         for(j=0;j<3;j++) p[j+1]=mom[3*(k-1)+j];
         boost(b,p);

         for(j=0;j<3;j++) hepeup_.PUP[shift+k][j]=p[j+1];
          
         hepeup_.MOTHUP[shift+k][0]=I+1;
         hepeup_.MOTHUP[shift+k][1]=I+1;
         hepeup_.SPINUP[shift+k]=9;
         hepeup_.VTIMUP[shift+k]=0;
         hepeup_.ICOLUP[shift+k][0]=0;
         hepeup_.ICOLUP[shift+k][1]=0;
         hepeup_.PUP[shift+k][4]=Finfo->pmass[k]; 
         hepeup_.ISTUP[shift+k]=1;
         hepeup_.IDUP[shift+k] =Finfo->PIDs[k];
       }
      
       for(i=0;clr1[i];i++)
       { int k1=clr1[i]-1, k2=clr2[i]-1;
         if(k1 && k2) 
         { hepeup_.ICOLUP[k1+shift][1]=++CC;
           hepeup_.ICOLUP[k2+shift][0]=CC;
         } else if(k1)  hepeup_.ICOLUP[k1+shift][1]=hepeup_.ICOLUP[I][1];
         else   if(k2)  hepeup_.ICOLUP[k2+shift][0]=hepeup_.ICOLUP[I][0];                      
                      
       }
       
       hepeup_.NUP+=Finfo->Nout;
       break;
     }
     }
           
   }
   
   for(I=0;I<hepeup_.NUP;I++)
   {
     hepeup_.PUP[I][3]=hepeup_.PUP[I][4]*hepeup_.PUP[I][4];
     for(j=0;j<3;j++)  hepeup_.PUP[I][3]+=hepeup_.PUP[I][j]*hepeup_.PUP[I][j];
     hepeup_.PUP[I][3]=sqrt(hepeup_.PUP[I][3]);
   }   
}


static void cleanList(eventfile_info * list)
{  eventfile_info * Finfo;
   for(;list;)
   { fclose(list->F); 
     free(list->fileName);
     Finfo=list;
     list=list->next;
     free(Finfo);
   }
}

void closeevents_(void)
{ decay_info * D;

  totCS=0; maxLum=0;
  cleanList(All); All=NULL;
  for(;Decays;)
  { cleanList(Decays->List);
    D=Decays;
    Decays=Decays->next;
    free(D);
  }   
}
