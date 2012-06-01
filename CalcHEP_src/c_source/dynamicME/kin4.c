#include<math.h>
#include<stdlib.h>
#include<stdio.h>

#include "../../include/V_and_P.h"
#include"dynamic_cs.h"

extern double simpson(double (*func)(double), double a, double b, double eps);
extern double alpha_2(double Q);
/*===========================================================*/

double (*sqme)(int nsub,double *pvect, int * err_code)=NULL;

static int  nsub_stat;

static double*Q=NULL;

double  decayPcm(double am0,  double  am1,  double  am2)
{
  double  summ, diffm;
  summ = am1 + am2;
  diffm = am1 - am2;
  if(am0<summ) return 0;
  return sqrt((am0-summ)*(am0+ summ)*(am0-diffm)*(am0+diffm))/(am0*2);
}
 



int pTabPos(char * name)
{
  int i;
  for(i=0;i<nModelParticles;i++)
  { 
    if(!strcmp(name,ModelPrtcls[i].name )) return   i+1;
    if(!strcmp(name,ModelPrtcls[i].aname)) return -(i+1);
  }
  return 0;
}

char * pdg2name(int pdg)
{
  int i;
  if(pdg==0) return NULL;

  for(i=0;i<nModelParticles;i++)
  {          if(ModelPrtcls[i].NPDG==pdg) return ModelPrtcls[i].name;
     else  { if(ModelPrtcls[i].NPDG==-pdg) return ModelPrtcls[i].aname;}
  }   
  return NULL;
} 

static void setQforParticle(char*pname)
{
  char *nm; 
  double *ma;  
  int n,i,pdg,cdim;
  if(!Q) return;

  n=pTabPos(pname);
  if(!n){printf("Wrong particle name '%s'\n",pname); return ;}
  nm=ModelPrtcls[abs(n)-1].mass;
  if(nm[0]=='0') return ; else ma=varAddress(nm);

  cdim=abs(ModelPrtcls[abs(n)-1].cdim);
  pdg=abs(ModelPrtcls[abs(n)-1].NPDG);

  if(cdim==1){ calcMainFunc(); *Q=*(varAddress(nm)); calcMainFunc(); return;}
  switch(pdg)
  { case 1:case 2:case 3: *Q=1; return;
    case 4: *Q=1.5; break;
    case 5: *Q=5;   break;
    case 6: *Q=175; break;
  }
  calcMainFunc();  
  for(i=0;i<10;i++)
  { 
    if( fabs(*Q-(*ma)) < 1E-2*(fabs(*Q))) break;
    *Q=*ma;
    calcMainFunc();
  }
}


double pMass(char * name)
{
  char *nm;
  int n=pTabPos(name);
  if(!n){printf("Wrong particle name '%s'\n",name); return 0;}
  nm=ModelPrtcls[abs(n)-1].mass;
  if(nm[0]=='0') return 0; else 
  { double *ma=varAddress(nm);
    return fabs(*ma);
  }
}


int ForceUG=0;



int  procInfo1(numout*cc, int *nsub, int * nin, int *nout)
{
  if(nin) *nin=cc->interface->nin;
  if(nout)*nout=cc->interface->nout;
  if(nsub)*nsub=cc->interface->nprc;
  return 0;
}

int procInfo2(numout*cc,int nsub,char**name,double*mass)
{
  int i;
  int ntot=cc->interface->nin+cc->interface->nout;
    
  if(nsub<1 || nsub> cc->interface->nprc) return 2;

  if(name)for(i=0;i<ntot ;i++) 
  name[i]=(cc->interface->pinf)(nsub,i+1,NULL,NULL);

  if(mass)
  {  
    for(i=1;i<=cc->interface->nvar;i++) if(cc->link[i]) cc->interface->va[i]=*(cc->link[i]);  
    if(cc->interface->calcFunc()>0) { printf("cannot calculate constr\n"); return 4;}       
    for(i=0;i<ntot ;i++) cc->interface->pinf(nsub,i+1,mass+i,NULL);     
  }
  return 0;
}




long pNum(char * name)
{
  int n=pTabPos(name);
  if(!n){printf("Wrong particle name ''%s''\n",name); return 0;}
  if(n>0)  return  ModelPrtcls[abs(n)-1].NPDG;
  else     return -ModelPrtcls[abs(n)-1].NPDG;
}

                    
void pname2lib(char*pname, char * libname)
{
  int n;
  n=pTabPos(pname);
  if(!n) {printf("Wrong particle name '''%s'''\n",pname); libname[0]=0; return;}
  if(n>0) sprintf(libname,"p%d",n); else sprintf(libname,"a%d",-n);
}


static int decodeProcess(char *txt,int*inList,int*outList)
{ char name[20];
   char *ch_,*ch;
   int i,p;   
   ch_=strstr(txt,"->");
   if(!ch_) { inList[0]=0; ch=txt;}
   else
   { 
     for(p=0,ch=txt;; )
     { sscanf(ch," %[^,]",name);
       ch_=strstr(name,"->");
       if(ch_) *ch_=0;
       for(i=strlen(name)-1; i>=0 && name[i]==' '; i--) name[i]=0;
       inList[p]=pTabPos(name);
       if(!inList[p]) return -(p+1);
       p++;
       if(ch_) break;     
       ch=strchr(ch,',');
       if(!ch) break; else ch++;
     }  
     inList[p]=0; 
     ch=strstr(txt,"->")+2; 
   }  

   for(p=0;ch; )
   { sscanf(ch," %[^,]",name);
     for(i=strlen(name)-1;i>=0 && name[i]==' '; i--) name[i]=0;
     outList[p]=pTabPos(name);
     if(!outList[p]) return p+1;
     p++;    
     ch=strchr(ch,',');
     if(!ch) break;
     ch++;
   }
   outList[p]=0;
   return 0;
}


 
void massFilter(double M, txtList * List)
{
  txtList lold=*List, lnew=NULL,lnext;

  while(lold)
  {  double Msum=0;
     char *ch;
     ch=strstr(lold->txt,"->")+2;
     for( ; ch; ch=strchr(ch,','))
     { char buff[10];
       int n;   
       char *nm;

       ch++;
       sscanf(ch,"%[^,]",buff);
       n=pTabPos(buff);
       nm=ModelPrtcls[abs(n)-1].mass;
    
       if(nm[0]=='0')
       { 
          switch(abs(ModelPrtcls[abs(n)-1].NPDG))
          {
            case 1: case 2: Msum+=0.07; break;
            case 3: Msum+=0.3; break;
            case 4: Msum+=1.5; break;
            case 5: Msum+=5. ; break;
          }
       }   
       else
       { double*ma =varAddress(nm);
         Msum+= fabs( *ma);
       }  

       Msum+=pMass(buff);
     } 
     lnext=lold->next;
     if(M>Msum) {lold->next=lnew; lnew=lold;}
     else {free(lold->txt); free(lold);}
     lold=lnext;
 } 
 *List=lnew;
}

void gammaGluFilter(txtList * List)
{
  txtList lold=*List, lnew=NULL,lnext;

  while(lold)
  {  int del=0,code;
     char *ch;
     ch=strstr(lold->txt,"->")+2;
     for( ; !del && ch; ch=strchr(ch,','))
     { char buff[10];
       ch++;
       sscanf(ch,"%[^,]",buff); 
       code=pNum(buff);
       if(code==22 || code ==21) { del=1;}
     } 
     lnext=lold->next;
     if(del) {free(lold->txt); free(lold);}
     else    {lold->next=lnew; lnew=lold;}
  
     lold=lnext;
 } 
 *List=lnew;
}


void process2Lib(char * process,char * lib)
{ 
  char * ch, *ch1;
  char bufflib[20];
  char *process_;
  
  process_=malloc(strlen(process)+1);
  strcpy(process_,process);
    
  ch= strstr(process_,"->");
  ch[0]=' '; ch[1]=','; ch=ch+2; for(;*ch==' ';ch++);

  lib[0]=0;
  ch1=strtok(process_," ,"); 
  for(;ch1;)
  { if(ch1==ch) strcat(lib,"_");
    pname2lib(ch1,bufflib);
    strcat(lib,bufflib);
    ch1=strtok(NULL," ,");
  }
  free(process_);
}

void process2Mass(char * process,double * mass)
{ 
  char * ch, *ch1;
  char *process_;
  int i;  

  process_=malloc(strlen(process)+1);
  strcpy(process_,process);
    
  ch= strstr(process_,"->");
  ch[0]=' '; ch[1]=','; ch=ch+2; for(;*ch==' ';ch++);

  ch1=strtok(process_," ,"); 
  for(i=0;ch1;i++,ch1=strtok(NULL," ,"))
  { /*if(ch1==ch) strcat(lib,"_");*/
    mass[i]=pMass(ch1);
  }
  free(process_);
}


/*===========================================================*/
static double pmass[4];

/* ======================decayTable ===================*/

 decayTableStr* decayTable=NULL;
 
 void cleanDecayTable(void)
 { int i,j;
   if(decayTable)
   {  
     for(i=0;i<nModelParticles;i++)
     {  
        decayTable[i].width=0;
        decayTable[i].dim=0;       
        for(j=0;j<2;j++)  
        if(decayTable[i].pdList[j]) cleanTxtList(decayTable[i].pdList[j]);
     }
   }  
   else  decayTable=malloc(nModelParticles*sizeof(decayTableStr));
   for(i=0;i<nModelParticles;i++)
   { for(j=0;j<2;j++) decayTable[i].pdList[j]=NULL;
     decayTable[i].width=0;
     decayTable[i].dim=0; 
   }
 }


/*======================  1->2 decay ==================*/

double pWidth2(numout * cc, int nsub)
{
  double pvect[12];
  double width=0.;
  double m1,m2,m3; 
  int i,ntot,nin,nout;

  procInfo1(cc,&ntot,&nin,&nout);
  if(nsub<1 ||  nsub>ntot|| nin!=1||nout !=2)  return 0;
     
  for(i=1;i<=cc->interface->nvar;i++) if(cc->link[i]) cc->interface->va[i]=*(cc->link[i]);

  if( cc->interface->calcFunc()>0 ) { return -1;}
  cc->interface->pinf(nsub,1,&m1,NULL);

  if(cc->GG) if(cc->GG) *(cc->GG)=sqrt(4*M_PI*alpha_2(m1));
  if(cc->Q) 
  {
    *(cc->Q)=m1;
    if( cc->interface->calcFunc()>0 ) { return -1;}
    if(cc->SC && cc->GG)  *(cc->GG)=*(cc->SC);           
  }  
  
  cc->interface->pinf(nsub,1,&m1,NULL);
  cc->interface->pinf(nsub,2,&m2,NULL); 
  cc->interface->pinf(nsub,3,&m3,NULL);
  if(m1 >m2 + m3)
  {   int i,err_code; 
      double md=m2-m3;
      double ms=m2+m3;
      double pRestOut=sqrt((m1*m1 - ms*ms)*(m1*m1-md*md))/(2*m1);
      double totcoef= pRestOut/(8. * M_PI * m1*m1);
           
      for(i=1;i<12;i++) pvect[i]=0;
      pvect[0]=m1;
      pvect[7]=pRestOut;
      pvect[4]=sqrt(pRestOut*pRestOut+m2*m2);
      pvect[11]=-pRestOut;
      pvect[8]=sqrt(pRestOut*pRestOut+m3*m3);
      width = totcoef * (cc->interface->sqme)(nsub,pvect,&err_code);
  }
  return width;
}

 
double decay2Info(char * pname, FILE* f)
{ int i,j,ntot;
  numout * cc;
  double wtot;
  char pname2[20],process[20],plib[20];
  char * dname[8];

  for(i=0,j=0;pname[i];i++)
  if(pname[i]!=' ') pname2[j++]=pname[i];
  pname2[j]=0;
  strcpy(plib,"2width_");
  pname2lib(pname2,plib+7);
  sprintf(process,"%s->2*x",pname2);
  cc=getMEcode(0,ForceUG,process,NULL,"",plib);
  if(!cc) return -1; 
  procInfo1(cc,&ntot,NULL,NULL);
  if(f) fprintf(f,"\n Partial width for %s->2x decays in GeV\n",pname2); 
  for(wtot=0,i=1;i<=ntot;i++)
  { double w;
    procInfo2(cc,i,dname,NULL);
    w=pWidth2(cc,i);
    if(w!=0)
    { wtot+=w;
      if(f) fprintf(f,"%3.3s %3.3s  %.2E\n",dname[1],dname[2],w); 
    }
  }
  if(f) fprintf(f," Total width %.2E GeV\n",wtot);
  return  wtot;
}

static double decay22List(char * pname, txtList *LL)
{ int i,j,ntot;
  numout * cc;
  double wtot,w;
  char pname2[20],process[20],plib[20];
  char * dname[8];
  txtList L=NULL,l;
  char buff[100];

  for(i=0,j=0;pname[i];i++)
  if(pname[i]!=' ') pname2[j++]=pname[i];
  pname2[j]=0;
  strcpy(plib,"2width_");
  pname2lib(pname2,plib+7);
  sprintf(process,"%s->2*x",pname2);
  cc=getMEcode(0,ForceUG,process,NULL,"",plib);
  if(!cc) { if(LL) *LL=NULL; return -1;} 
  procInfo1(cc,&ntot,NULL,NULL);
  for(wtot=0,i=1;i<=ntot;i++)
  { 
    procInfo2(cc,i,dname,NULL);
    w=pWidth2(cc,i);
    if(w!=0)
    {  if(LL)
       { l=malloc(sizeof(txtListStr));
         l->next=L;
         L=l;
         sprintf(buff,"%E  %s -> %s,%s",w,pname2,dname[1],dname[2]);
         l->txt=malloc(20+strlen(buff));
         strcpy(l->txt,buff);
       } 
       wtot+=w; 
    }
  }

  if(LL)
  {  for(l=L;l;l=l->next)
     { 
       sscanf(l->txt,"%lf %[^\n]",&w,buff);
       sprintf(l->txt,"%E %s",w/wtot,buff);
     }   
    *LL=L;
  }
  return  wtot;
}


/*==================  1->3 decay ===================*/

static double _x_;

static double kimematic_1_3(double *pmass, int kin, double xm2, double xcos, double * P)
{ 
  double factor,pout,mQmin,mQmax,mQ,m12,chY,shY,xsin;
  double p4,p8,p5,p9;
  double m0=pmass[0],m1=pmass[1],m2=pmass[2],m3=pmass[3];
  int i;
  
  P[0]=m0; P[1]=P[2]=0;
  for(i=0;i<4;i++) P[3+i*4]=0;
  
  mQmin=(m1+m2)*(m1+m2);
  mQmax=(m0-m3)*(m0-m3);
  mQ=mQmin*(1-xm2)+mQmax*xm2;
  m12=sqrt(mQ);
   
  factor=(mQmax-mQmin)/(128*M_PI*M_PI*M_PI*m0*m0*m12);
  
  pout=decayPcm(m0,m12,m3);
  P[12]=sqrt(pout*pout+m3*m3); P[13]=-pout; P[14]=0;

  factor*=pout;  
  
  chY=sqrt(1+pout*pout/mQ);
  shY=sqrt(pout*pout/mQ);
  
  pout=decayPcm(m12,m1,m2);
  factor*=pout;
  xsin=sqrt(1-xcos*xcos);
  p4=sqrt(m1*m1+pout*pout);    p8=sqrt(m2*m2+pout*pout);
  p5=xcos*pout;                p9=-p5;
  P[6]=xsin*pout;              P[10]=-P[6];
  
  P[4]=chY*p4 + shY*p5;    P[8]=chY*p8 + shY*p9;
  P[5]=shY*p4 + chY*p5;    P[9]=shY*p8 + chY*p9;
    
  return factor;
}


static double dWidthdCos(double xcos)
{
  double factor;
  double P[16];
  int err_code;
  factor=kimematic_1_3(pmass,1,_x_,xcos, P); 
  return factor*(*sqme)(nsub_stat,P,&err_code);
}

static double dWidthdM(double x)
{ _x_=x; return simpson(dWidthdCos,-1.,1.,1.E-4); }



static double width13(numout * cc, int nsub, int * err) 
{
  int i;
  for(i=1;i<=cc->interface->nvar;i++) if(cc->link[i]) cc->interface->va[i]=*(cc->link[i]);

  if( cc->interface->calcFunc()>0 ) {*err=4; return 0;}

  for(i=0;i<4;i++) cc->interface->pinf(nsub,1+i,pmass+i,NULL); 


  if(cc->GG)
  {  if(cc->SC) *(cc->GG)=*(cc->SC); 
     else       *(cc->GG)=sqrt(4*M_PI*alpha_2(pmass[0])) ;
  }
  
  *(cc->interface->gtwidth)=0;
  *(cc->interface->twidth)=0;
  *(cc->interface->gswidth)=0;
    
  
  *err=0;
  sqme=cc->interface->sqme;
  nsub_stat=nsub; 
  return simpson(dWidthdM,0.000001,0.99999,1.E-2);
}

static txtList conBrList(txtList BrList)
{ txtList out=NULL;
  char buff[100];
  double br;
  int inCode[10], outCode[10],i;
  for(;BrList;BrList=BrList->next)
  { txtList new=malloc(sizeof(txtListStr));
    new->next=out;out=new;
    sscanf(BrList->txt,"%lf %[^\n]",&br,buff); 
    decodeProcess(buff,inCode,outCode);
    if(inCode[0]>0) sprintf(buff,"%E  %s -> ",br,ModelPrtcls[inCode[0]-1].aname);   
    else            sprintf(buff,"%E  %s -> ",br,ModelPrtcls[-inCode[0]-1].name);
    for(i=0;outCode[i];i++)
    { if(i) strcat(buff,",");
      if(outCode[i]>0) strcat(buff,ModelPrtcls[outCode[i]-1].aname);   
      else             strcat(buff,ModelPrtcls[-outCode[i]-1].name);  
    }
    new->txt=malloc(strlen(buff)+1);
    strcpy(new->txt,buff);
  }
  return out;  
}


double pWidth(char *name, txtList * LL,int *dim)
{
  txtList L,l,Lout;
  char libName[100];
  double sum=0,width;
  int i,i0,j,j0;
  double Qstat;
  
  for(i=0;i<nModelParticles;i++)
  { char *pnames[2]={ModelPrtcls[i].name,ModelPrtcls[i].aname};
    for(j=0;j<2;j++)
    if(strcmp(name,pnames[j])==0) 
    { if(decayTable[i].dim)
      { if(dim)*dim=decayTable[i].dim;
        if(LL) *LL=decayTable[i].pdList[j];
        return decayTable[i].width;
      } else break;
    } if(j!=2) break;    
  }    
  i0=i,j0=j;
  if(i0==nModelParticles)
  { printf("%s out of model particles\n",name);
    if(LL) *LL=NULL;
    if(dim)*dim=0;
    return 0;
  }  
 
  if(Q==NULL)  for(i=0;i<nModelVars;i++)
                   if(strcmp(varNames[i],"Q")==0){ Q= varValues+i; break;}
  if(Q) { Qstat=*Q; setQforParticle(name);}
 
  width=decay22List(name,&L);
  
  if(L) 
  { if(dim) *dim=2;
    if(LL) *LL=L;
    decayTable[i0].pdList[j0]=L;
    if(strcmp(ModelPrtcls[i0].name,ModelPrtcls[i0].aname)) 
                 decayTable[i0].pdList[1-j0]=conBrList(L);
    
    decayTable[i0].dim=2;
    decayTable[i0].width=width;  
    if(Q) { *Q=Qstat; calcMainFunc();}
    return width;
  }

  L= makeDecayList(name,3);
  massFilter(pMass(name),&L);
  gammaGluFilter(&L);
  Lout=NULL;

  decayTable[i0].dim=3;
  for(sum=0,l=L;l;l=l->next)  
  { numout* cc;
    int err;
    txtList newr;
    
    process2Lib(l->txt ,libName);
    cc=getMEcode(0,ForceUG,l->txt,NULL,"",libName);
    width=width13(cc, 1, &err);
    sum+=width;
    newr=malloc(sizeof(txtListStr));
    newr->next=Lout;
    Lout=newr;
    newr->txt=malloc(strlen(l->txt)+20);
    sprintf(newr->txt,"%E  %s",width,l->txt);
  }
  cleanTxtList(L);
  if(Lout)
  for(L=Lout;L;L=L->next)
  { char buff[100];
    sscanf(L->txt,"%lf %[^\n]",&width,buff);
    sprintf(L->txt,"%E %s",width/sum,buff);  
  }   
  if(LL) *LL=Lout;
  if(dim) *dim=3;
  decayTable[i0].pdList[j0]=Lout;
  if(strcmp(ModelPrtcls[i0].name,ModelPrtcls[i0].aname)) 
               decayTable[i0].pdList[1-j0]=conBrList(Lout);
  decayTable[i0].dim=3;
  decayTable[i0].width=sum;  
  if(Q) { *Q=Qstat; calcMainFunc();}
  return sum;
}

static int pListEq(char * txt1, char * txt2)  
{  char buff[100];
   char rd1[10][10];
   char rd2[10][10];
   int n1,n2,i1,i2;
   char *ch;
    
   strcpy(buff,txt1); while((ch=strchr(buff,','))) ch[0]=' ';
   
   n1=sscanf(buff,"%s %s %s %s %s %s %s %s %s %s",
   rd1[0],rd1[1],rd1[2],rd1[3],rd1[4],rd1[5],rd1[6],rd1[7],rd1[8],rd1[9]); 
   
   strcpy(buff,txt2); while((ch=strchr(buff,','))) ch[0]=' ';
   
   n2=sscanf(buff,"%s %s %s %s %s %s %s %s %s %s",
   rd2[0],rd2[1],rd2[2],rd2[3],rd2[4],rd2[5],rd2[6],rd2[7],rd2[8],rd2[9]); 
   
   if(n1!=n2) return 0;
   for(i1=0;i1<n1;i1++)
   { for(i2=0;i2<n2;i2++) if(strcmp(rd1[i1],rd2[i2])==0){rd2[i2][0]=0; break;}
     if(i2==n2) return 0;
   } 
   return 1;
}      

double findBr(txtList L, char * pattern)
{ char buff[100];
  char *ch;
  double width;
  
  for(;L;L=L->next)
  { 
     sscanf(L->txt,"%lf %[^\n]",&width,buff);
     ch=strstr(buff,"->");
     ch+=2;
     if( pListEq(ch,pattern)) return width;
  }
  return 0;   
}


