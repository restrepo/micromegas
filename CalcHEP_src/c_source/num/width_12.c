/*
 Copyright (C) 1997,2006  Alexander Pukhov 
*/

#include <math.h>
#include "chep_crt.h"
#include "err_code.h"
#include "plot.h"
#include "num_serv.h"
#include "width_12.h"
#include "const.h"
#include "interface.h"
#include "param.h"
#include "4_vector.h"
#include "files.h"
#include "alphas2.h"

static char inParticle[10];
static double * widths=NULL;
static int EffQmass=1;
static  double *Q=NULL,*GG=NULL,*SC=NULL;
static int nsubSel=0;


static void  decay12information(double totwidth,int Branchings)
{  
   clrbox(1,1,53,16);  
   clrbox(1,16, maxRow(), maxCol());
   
   goto_xy(5,3);scrcolor(Red,BGmain); print(" Decay ");
   scrcolor(Blue,BGmain);print("%s -> 2*x      ",inParticle);
   goto_xy(1,4);
   scrcolor(Red,BGmain); 
   print(" Total width : "); 
   scrcolor(FGmain,BGmain); 
   if(err_code) {print(" incorrect       "); return;}

   print("%.3E GeV     ",totwidth);
   if (totwidth > 0. ) 
   {  int i,xcount = 31, ycount = 15;
      int * sort= (int*) malloc(sizeof(int)*nprc_int);

      for(i=0;i<nprc_int;i++) sort[i]=i;  
      for(i=0;i<nprc_int-1;)
      {   
         if(widths[sort[i]] < widths[sort[i+1]])
         {  int buff=sort[i];
            sort[i]=sort[i+1];
            sort[i+1]=buff;
            if (i!=0) i--;
         } else i++;
      }   

      goto_xy(1,6);
      scrcolor(Red,BGmain);
      if(Branchings) print(" Modes and fractions :");
      else           print(" Partial widths [Gev] :"); 
      for(i=0,xcount=5,ycount=7;i<nprc_int;i++,ycount++) 
      if(widths[sort[i]]>0)
      {  if(ycount >= maxRow()) 
         { ycount =16; xcount += 25;
           if(xcount > maxCol()-20) break;
           if(xcount >30) ycount=16; else ycount=7;
         }
	 goto_xy(xcount,ycount);
         scrcolor(Blue,BGmain);
         print("%3s ",pinf_int(sort[i]+1,2,NULL,NULL));
         print("%3s - ", pinf_int(sort[i]+1,3,NULL,NULL));
         scrcolor(FGmain,BGmain);
         if(Branchings) print("%8.2E%%", 100*widths[sort[i]]/totwidth); 
         else           print("%9.3E", widths[sort[i]]); 
      }  
      free(sort);
   } 
   scrcolor(FGmain,BGmain); 
} 


static void  writeLesHdecays(void)
{ 
 int i,nsub;
 int first=1;
 
 long N1,N2,N3;
 double S; 
 double *widths=(double *) malloc(sizeof(double)*nprc_int);
 double *sum=(double *) malloc(sizeof(double)*nprc_int);

 double *masses=(double *) malloc(3*sizeof(double)*nprc_int);
 long  *codes=(long *) malloc(3*sizeof(long)*nprc_int);
 char ** names=malloc(3*sizeof(char*)*nprc_int); 
 int npTot=0;
 char f_name[50];
 FILE *f;
 char*inP=" ";

 nextFileName(f_name,"decaySLHA",".txt");
 f=fopen(f_name,"w");
fprintf(f,"#Masses, widths and branchings in SLHA format:\n");
fprintf(f,"BLOCK MODELPARAMETERS:\n");
for(i=1;i<=nvar_int;i++) 
fprintf(f,"# %8s  %E\n", varName_int[i], va_int[i]);
fprintf(f,"BLOCK MASS  # Mass Spectrum\n");
for(nsub=1;nsub<=nprc_int;nsub++)
{ for(i=1;i<4;i++) 
  {  double m;
     long N; 
     char * nm=pinf_int(nsub,i,&m,&N);
     int j;
     for(j=0;j<npTot;j++) if(strcmp(nm,names[j])==0 || N==codes[j] ||
     N==-codes[j]) break;
     if(j==npTot) 
     { masses[j]=m;
       codes[j]=N;
       names[j]=nm;
       npTot++;
     }
  }  
} 
for(i=0;i<npTot;i++) fprintf(f,"%10d  %E  # %s\n",codes[i],masses[i],names[i]);  
free(masses); free(codes); free(names);
  
 err_code = 0; 

 for(nsub=0;nsub<nprc_int;nsub++) sum[nsub]=0;   

 for(nsub=1;nsub<=nprc_int;nsub++) 
 {  double m1, m2, m3;
    char *inP_new=pinf_int(nsub,1,&m1,&N1);
    int Nsum;
        
    first=strcmp(inP_new,inP); 
      
    if(first)
    { Nsum=nsub-1; 
      inP=inP_new;
      if(EffQmass&&Q) *Q=m1;
      if(calcFunc_int()>0)
      {  messanykey(15,15,"Can not calculate constraints");
           return;
      }
      if(GG)
      { if(SC) *GG=*SC; 
        else {if(Q) *GG=sqrt(4*M_PI*alpha_2(*Q)); else *GG=sqrt(4*M_PI*alpha_2(m1));}
      }
    }  
    
    pinf_int(nsub,1,&m1,&N1); 
    pinf_int(nsub,2,&m2,NULL);  
    pinf_int(nsub,3,&m3,NULL);

          
    if (m1 <=m2 + m3) widths[nsub-1] = 0.0; 
    else 
    { 
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
            
        widths[nsub-1] = totcoef * sqme_int(nsub,pvect,&err_code);
        if(err_code != 0) {  errormessage(); widths[nsub-1]=0; err_code=0;}
        sum[Nsum] += widths[nsub-1];      
    }
 }

 inP=" "; 
 for(nsub=1;nsub<=nprc_int;nsub++) 
 {     
    char *inP_new=pinf_int(nsub,1,NULL,&N1);
    
    first=strcmp(inP_new,inP); 
      
    if(first)
    {  inP=inP_new;
       S=sum[nsub-1];
       if(S>0) fprintf(f,"DECAY   %d   %E # %s width\n",  N1, S,  inP);
    }
    if(S>0 && widths[nsub-1]!=0)
    {  char * inP2=pinf_int(nsub,2,NULL,&N2);
       char * inP3=pinf_int(nsub,3,NULL,&N3);
       fprintf(f," %12.4E 2  %8d  %8d # Br(%-3s -> %-3s %-3s) %11.3E[Gev]\n",
        widths[nsub-1]/S,N2,N3,inP,inP2,inP3,widths[nsub-1]);
    }
 }
 fclose(f);
{ char mess[100];
  sprintf(mess,"Output is saved  in %s",f_name);   
  messanykey( 12, 12,mess);  
}
 free(widths); free(sum);

} 

static double calcwidth12(void)
{ 
 int i,nsub;
 double width12 = 0.;
 double selChan=0;
 int first=1;
 
 long N1;

 err_code = 0; 

 for(nsub=1;nsub<=nprc_int;nsub++) widths[nsub-1]=0;
 for(nsub=1;nsub<=nprc_int;nsub++) 
 {  double m1, m2, m3;

    if(strcmp(pinf_int(nsub,1,&m1,&N1),inParticle)==0) 
    { 
      if(first)
      { 
        if(EffQmass&&Q) *Q=m1;
        if(calcFunc_int()>0)
        {  messanykey(15,15,"Can not  calculate constraints");
           return 0;
        } 
        if(GG)
        { if(SC) *GG=*SC; 
          else {if(Q) *GG=sqrt(4*M_PI*alpha_2(*Q)); else *GG=sqrt(4*M_PI*alpha_2(m1));}
        }

        first=0;
      }  
      pinf_int(nsub,1,&m1,NULL);pinf_int(nsub,2,&m2,NULL);pinf_int(nsub,3,&m3,NULL);
          
      if (m1 <=m2 + m3) widths[nsub-1] = 0.0; 
      else 
      { 
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
            
        widths[nsub-1] = totcoef * sqme_int(nsub,pvect,&err_code);
        if(err_code != 0) {  errormessage(); widths[nsub-1]=0; err_code=0;}
        width12 += widths[nsub-1];
        if(nsubSel==nsub) selChan= widths[nsub-1];
      }
    }  
 }
 if(nsubSel) { if(width12) return selChan/width12; else return 0;}  
 return width12; 
} 

static void inmenutxt(char ** menutxt)
{ 
  char * n1,*n2;

  int i,pos=11;

  *menutxt=(char *)malloc(2+10*nprc_int);
  *menutxt[0]=10;

  n1=pinf_int(1,1,NULL,NULL);
  sprintf(*menutxt+1," %-9.9s",n1);
  for(i=2;i<=nprc_int;i++) 
  { n2=pinf_int(i,1,NULL,NULL);
    if(strcmp(n1,n2)!=0)
    { n1=n2;
       sprintf(*menutxt+pos," %-9.9s",n1);
       pos+=10;  
    }
  }
}


static void * selectChan(void)
{  
  int nsub, ntot;
  char * menutxt;
  void * pscr=NULL;
   
  for(nsub=1,ntot=1;nsub<=nprc_int;nsub++) 
      if(strcmp(pinf_int(nsub,1,NULL,NULL),inParticle)==0) ntot++;
   
  menutxt=malloc(2 +15*ntot);
  menutxt[0]=15; menutxt[1]=0;
  sprintf(menutxt+1," total width   ");

  for(nsub=1,ntot=1;nsub<=nprc_int;nsub++)
         if(strcmp(pinf_int(nsub,1,NULL,NULL),inParticle)==0)
  { sprintf(menutxt+1+15*ntot++," BR(%-3.3s , %-3.3s) ", pinf_int(nsub,2,NULL,NULL),
     pinf_int(nsub,3,NULL,NULL));
  } 

  for(ntot=0;ntot==0;) menu1(56,7,"Select",menutxt,"",&pscr,&ntot); 
  if(ntot==0) return NULL;  
  ntot--;   
  if(ntot==0) nsubSel=0; else 
  {
    for(nsub=1;nsub<=nprc_int&& ntot;nsub++)
         if(strcmp(pinf_int(nsub,1,NULL,NULL),inParticle)==0) ntot--;
    nsubSel=nsub-1;
  }
      
  return pscr; 
}

void  decay12(void)
{ 
   int  i, k,L;
   void * pscr=NULL; 
   char * mlist;
   static int Branch=1;

   widths=(double*)malloc(sizeof(double)*nprc_int);

   for(i=1;i<=nvar_int;i++)
   {   if(!strcmp(varName_int[i],"Q"))  Q=va_int+i;
       else if(!strcmp(varName_int[i],"GG")) GG=va_int+i;
   }
 
   if(GG)for(i=1;i<=nvar_int+nfunc_int;i++) 
      if(!strcmp(varName_int[i],"SC")){ SC=va_int+i; break;}
   
   inmenutxt(&mlist);
   L=mlist[0];
   sscanf(mlist+1,"%s",inParticle);   

   for(k=1;k;) 
   {  
      char strmen[]="\030"        
         " Incoming particle      "
         " Show Branchings        "
         " QCD Scale Q= Free      "
         " Model parameters       "
         " Constraints            "
         " Parameter dependence   "
         " Les Houches output     ";

      clrbox(1,13, maxCol(), maxRow());
      nsubSel=0;
      decay12information(calcwidth12(),Branch);

      if(EffQmass) improveStr(strmen,"Free ","M1");
      if(!Branch)  improveStr(strmen,"Branchings","Partial widths"); 
      menu1(54,4,"",strmen,"n_12_*",&pscr,&k);

      switch (k)
      { 
        case 1:
           {
             if(strlen(mlist)>L+2)
             { void * pscr2=NULL;
               int k=1;
                menu1(56,5,"",mlist,"",&pscr2,&k);
               if(k)  sscanf(mlist+(k-1)*L+1,"%s",inParticle);
               put_text(&pscr2);
             }
           }
           break;
        case 2: Branch=!Branch;     break;
        case 3: EffQmass=!EffQmass; break;
        case 4: change_parameter(54,8,0); break;
        case 5: show_depend(54,8); break;
        case 6:
           { char proc[20];
             char dimInfo[20]="Width  [GeV]";
             void * pscr=selectChan();
             if(!pscr) break;
             if(nsubSel==0) sprintf(proc,"%s -> 2*x",inParticle); else
             { sprintf(proc," BR(%s ->  %s %s)",inParticle, 
                pinf_int(nsubSel,2,NULL,NULL),  pinf_int(nsubSel,3,NULL,NULL));
                dimInfo[0]=0;
             }  
	     paramdependence( calcwidth12,proc,dimInfo);
	     put_text(&pscr);
	   } break;
        case 7: writeLesHdecays(); break; 
      } 
   }
   free(widths);
   free(mlist);   
   clrbox(1,1,53,16);
   clrbox(1,16,maxCol(),maxRow());
}
