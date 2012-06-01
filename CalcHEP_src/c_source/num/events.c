/*
 Copyright (C) 2001, Alexander Pukhov, e-mail pukhov@theory.sinp.msu.ru
*/
 
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"vegas.h"
#include"drandXX.h"
#include"chep_crt.h"
#include"interface.h"
#include"4_vector.h"
#include"subproc.h"
#include"alphas2.h"
#include<unistd.h>
#include"strfun.h"
#include"files.h"
#include"runVegas.h"
#include"events.h"

static FILE * events_;

static void writeEvent(long cCub,  int n, char * rand_state)
{ 
   int i;
   int icc;
#ifdef PARKED 
 fprintf(events_,"%05X|%s|%d\n",cCub,rand_state,n);
#else

   if(*cb_pow_int)
   {  double sum=0;
      for(i=0;i<*cb_pow_int;i++) sum+=fabs((*cb_coeff_int)[i]);
      sum*=drandXX();
      for(i=0;i<*cb_pow_int;i++)
      { sum-=fabs((*cb_coeff_int)[i]);
        if(sum<=0) break;
      } 
      if(i==*cb_pow_int) i--;
      if((*cb_coeff_int)[i]<0) n*=-1;
      icc=i;
   }
   

   fprintf(events_,"%8d ",n);  
   if(nin_int==2) fprintf(events_," %17.10E %17.10E",pvect[3],pvect[7]);

   for(i=0;i<nout_int;i++) fprintf(events_," %17.10E %17.10E %17.10E",
   pvect[4*(i+nin_int)+1],pvect[4*(i+nin_int)+2],pvect[4*(i+nin_int)+3]);

   fprintf(events_,"| %10.3E ", Scale());

   
   if(*cb_pow_int)
   {  int j;
      fprintf(events_,"  ");
      for(j=0;j<*cb_nc_int;j++)
      fprintf(events_,"(%d %d)",(*cb_chains_int)[2*(*cb_nc_int)*icc+2*j], 
                                (*cb_chains_int)[2*(*cb_nc_int)*icc+2*j+1]);
   }


   fprintf(events_,"\n");

#endif   
}

static long n_cube=1000;
   

static long nRandom=100;
static long nSimplex1=50;
static long nSimplex2=50;
static long nEvents=10000;
static double max =2;
static double milk=0.; 

int saveEventSettings(FILE * f)
{
  fprintf(f, "%ld %ld %ld %ld %ld %f",n_cube, nRandom, nSimplex1, nSimplex2,  nEvents, max);
  return 0;

}

int  readEventSettings(FILE * f)
{
  fscanf(f, "%ld %d %d %ld %ld %lf", &n_cube, &nRandom, &nSimplex1, &nSimplex2, &nEvents, &max);
  return 0;
}


static void write_event_cap(void) 
{
  int i,j;

  fprintf(events_,"#%s\n",version);
  fprintf(events_,"#Type %d -> %d\n",nin_int,nout_int); 
  fprintf(events_,"#Initial_state\n");
  fprintf(events_,"  P1_3=%E" , inP1); 
  if(nin_int>1) fprintf(events_,"  P2_3=%E\n",-inP2);else fprintf(events_,"\n");
  if(nin_int>1) wrt_sf__(events_);

  fprintf(events_,"#PROCESS  ");
  for(i=1;i<=nin_int+nout_int; i++)
  { long pcode;
    char * pname=pinf_int(Nsub,i,NULL,&pcode);
    switch(pcode)
    { case  81: pcode= 1; break;
      case -81: pcode=-1; break;
      case  83: pcode= 3; break;
      case -83: pcode=-3; break;
    }
    fprintf(events_," %d(%s)", pcode, pname);
    if(i==nin_int)  fprintf(events_," ->");
  } 
  fprintf(events_,"\n");    
  fprintf(events_,"#MASSES ");
  for(i=0;i<nin_int+nout_int;i++)
  {  double m;
     pinf_int(Nsub,i+1,&m,NULL); 
     fprintf(events_," %.10E", m);
  }   
  fprintf(events_,"\n");
  
  if(integral.n_it) fprintf(events_,"#Cross_section(Width) %E\n",
                                                     integral.s1/integral.s0);
  else   fprintf(events_,"#Cross_section(Width) Unknown\n"); 
  fprintf(events_,"#Number_of_events %10d\n",0);

  fprintf(events_,"#Events  "); 
  if(nin_int==2) fprintf(events_,"     P1_3 [Gev]        P2_3 [Gev]   ");
  for(i=1;i<=nout_int; i++) for(j=1;j<=3;j++) 
                          fprintf(events_,"     P%d_%d [Gev]   ",i+nin_int,j);

  fprintf(events_,"  QCD SCALE    Color chains\n");

}



   
void  generateEvents( vegasGrid * vegPtr,  
                   double (*func)(double *,double),  char *fname,  FILE * iprt)
{                                                   
   int mode=1;
   void * pscr=NULL;
   void * pscr_=NULL;
   double max=2;
   double eff0;

ret:  
   if(!vegPtr->fMax)for(mode=1;mode!=4; )
   {   
      char strmen[]="\030"
                    " sub-cubes = N1         "
                    " random search = N2     "
                    " simplex search= N3     "
                    " Start search of maxima "; 
      improveStr(strmen,"N1","%d",n_cube);
      improveStr(strmen,"N2","%d",nRandom);
      improveStr(strmen,"N3","%d",nSimplex1);
      menu1(54,11,"Preparing of generator",strmen,"n_prep_gen_*",&pscr_,&mode);
      switch(mode)
      {
        case 0: return;
        case 1: correctLong(50,15,"Number of sub-cubes:",&n_cube,1); break;
        case 2: correctLong(50,15,"Random search:",&nRandom,1);break;
        case 3: correctLong(50,15,"Simplex steps :",&nSimplex1,1); break; 

/*                    " milk      = N3         "  
      improveStr(strmen,"N3","%.1f",milk);
        case 4: correctDouble(50,15,"Content of milk:",&milk,1); break;
*/
        case 4: if(n_cube < 1 ) n_cube=1; 
        { 
          float * t1=malloc(n_cube*sizeof(float));
          float * t2=malloc(n_cube*sizeof(float));
          float * t3=malloc(n_cube*sizeof(float));
          int mCheck=(t3!=NULL);
          if(t1)free(t1);if(t2)free(t2);if(t3)free(t3); 

          if(mCheck)                                                         
          {                                                   
             if(vegas_max(vegPtr,n_cube, nRandom,nSimplex1, func,&eff0)==0)
             { char mess[50];
               sprintf(mess,"Expected efficiency %f",eff0);
               messanykey(25,15,mess);
               put_text(&pscr_);
             } else mode=1;
          } else                                                                       
          { 
             messanykey(25,15,"Not enough memory.\n"                                    
                              "Decrease the number of sub-cubes");                     
             mode=1;                                                                  
          }                                                                            
        }                                                                             
      }
   }
    
   for(mode=1;;)
   {
     char strmen[]="\030"
                   "Number of events=N1     "
                   " MAX*N2                 "
                   " find new MAX = N3      "
                   " Launch generator       "
		   " New search of maxima   ";
     improveStr(strmen,"N1","%d",nEvents);
     improveStr(strmen,"N2","%.2G",max);
     improveStr(strmen,"N3","%d",nSimplex2);
     menu1(54,10,"",strmen,"n_gen_*",&pscr,&mode);
     switch(mode)
     { case 0: return; 
       case 1: correctLong(50,15,"",&nEvents,1); break;
       case 2: correctDouble(50,15,"Increase max in ",&max,1); break;
       case 3: correctLong(50,15,"Number of steps:",&nSimplex2,1); break;
       case 4: 
       { long  nGenerated=0; 
         double eff, rmax, mult, neg;   
         char mess[200];
         int broken;
         long fileEnd;
         events_= fopen(fname,"a");       
         if(ftell(events_)==0) write_event_cap();
         fileEnd=ftell(events_); 
         build_cb_int(Nsub);
         broken= vegas_events(vegPtr,nEvents,max,milk,nSimplex2,func,writeEvent,
         &eff,&rmax,&mult,&neg);
         fclose(events_);

         sprintf(mess,"Statistic\n efficiency: %.1E\nReached max: %.1E\n"
                      "Mult. events: %.1E \nNeg.events: %.1E \n",
                      eff,rmax,mult, neg);

         if(broken)  
         {  strcat(mess,"---------------\n Events are not saved"); 
            messanykey(25,15,mess);
            truncate(fname,fileEnd);
         }
         else   
         {  int l=strlen(mess);
            strcat(mess,"---------------\n Accept events? ");    
            if(mess_y_n(25,15,mess)) 
            {  
               long  nEvPos=0;
               integral.old=1;
               mess[l]=0;
               events_=fopen(fname,"r+");
               while(nEvPos==0)
               { char ch;
                 char word[100];
                 do fscanf(events_,"%c",&ch); while(ch !='#');
                 fscanf(events_,"%s",word);
                 if(strcmp(word,"Number_of_events")==0) nEvPos=ftell(events_);
               }
               fscanf(events_,"%ld",&nGenerated);
               nGenerated+=nEvents;
               fseek(events_,nEvPos,SEEK_SET);
               fprintf(events_," %10d",nGenerated);
               fclose(events_);

               fprintf(iprt," %d events are stored in '%s'\n",nEvents,fname);
               fprintf(iprt,"%s\n",mess);
               fflush(iprt);
            } else  truncate(fname,fileEnd);
         }
         destroy_cb_int();
       } 
       break;

       case 5:
       free(vegPtr->fMax); 
       vegPtr->fMax=NULL;
       vegPtr->nCubs=0;
       put_text(&pscr);
       goto ret;

     }
   }
}                                                                     
