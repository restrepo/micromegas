/*
 Copyright (C) 1997, Alexander Pukhov 
*/

#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>

#include "pdt.h"
#include "files.h"
#include "interface.h"
#include "subproc.h"
#include "chep_crt.h"
#include "strfun.h"
#include "n_calchep_.h"
#include "alphas2.h"
#include "sf_pdt.h"


static int alphaMode=0;

static char * pdtName[2]={NULL,NULL}; 
static pdtStr pdtData[2];
static pdtStr pdtDataAux[2];
static int  statAux[2]={0,0};

static double alpha_pdt(double q ){ return interAlpha(q, pdtData+alphaMode-1);}

static void comphepPdtList(long pNum,  pdtList **list)
{  int i;
   char * dNames[3];
   char fname[STRSIZ],rootPdt[STRSIZ];
   long pNum1,pNum2;

   sprintf(rootPdt,"%s/pdTables", pathtocomphep);


   switch(pNum)
   { case  81: pNum1= 1; pNum2= 3; break;
     case -81: pNum1=-1; pNum2= 3; break;
     case  83: pNum1= 3; pNum2= 1; break;
     case -83: pNum1=-3; pNum2=-1; break;
     default:  pNum1=pNum; pNum2=0;
   }

   dNames[0]="..";
   dNames[1]=".";
   
   dNames[2]=rootPdt;

   for(i=0 ;i<3;i++)
   {  DIR *dirPtr=opendir(dNames[i]);
      struct dirent * dp;

      if(!dirPtr) continue;
 
      while(dp=readdir(dirPtr))
      { char *c=dp->d_name;
        int l=strlen(c);
        if(l>=4 && strcmp(c+l-4,".pdt")==0) 
        { sprintf(fname,"%s%c%s",dNames[i],f_slash,c); 
          if(pNum2==0) makePdtList(fname, pNum1, list);
          else
          { pdtList * list1=NULL,*l1,*l;
            pdtList * list2=NULL,*l2;
            makePdtList(fname, pNum1, &list1);
            makePdtList(fname, pNum2, &list2);
            for(l1=list1;l1;)
            {
               for(l2=list2;l2;l2=l2->next) if(strcmp(l1->name,l2->name)==0)
               { l1->posaux=l2->position;
                 l=l1; l1=l1->next;
                 l->next=*list;
                 *list=l;
                 break;
               }
               if(l2==NULL) 
               { l=l1;
                 l1=l1->next;
                 free(l->name); free(l->file); free(l); 
               }
            }
            delPdtList(l2); 
          }
        }
      } 
      closedir(dirPtr);
   }
} 

int mc_pdt(int i) {return pdtData[i-1].beamP;}

int p_pdt(long pNum) 
{  
   pdtList *list=NULL; 
   comphepPdtList(pNum, &list);
   if(list)  {delPdtList(list); return 1;} else return  0; 
}

void n_pdt(int i, char *name) 
{i--;  if(pdtName[i]) strcpy(name,pdtName[i]); else strcpy(name,"PDT:"); }


int r_pdt(int i, char *name)
{   
   i--;
   if(name==strstr(name,"PDT") ) 
   {
      if(pdtName[i]) free(pdtName[i]);
      pdtName[i]=malloc(strlen(name)+1);
      strcpy(pdtName[i],name);
      return 1;
   }else return 0; 
}


int be_pdt(int i,double * be, double * mass) 
{  
   pdtList *list=NULL;
   pdtList *list_;
   int ret_code=0; 
   static int  stat[2]={0,0};
   long pNum, N1, N2;
   i--;

   if(stat[i])    freePdtData(pdtData+i); 
   if(statAux[i]) { freePdtData(pdtDataAux+i); statAux[i]=0;}  

   pinf_int(Nsub,1,NULL,&N1);
   pinf_int(Nsub,2,NULL,&N2);

   if(i==0) pNum=N1; else pNum=N2;
   if( (abs(N1)==81 ||abs(N1)==83) && (abs(N2)==81 ||abs(N2)==83))
   { if(pNum>0) pNum-=80; else pNum+=80;}
   
   comphepPdtList(pNum, &list);   
   list_=list;
   while(list_ && ret_code==0)  
   { if(strcmp(list_->name,pdtName[i]+4)==0) 
     { /*printf("pNum=%d  name=%s N1=%d ",pNum,list_->name,list_->position);*/
       ret_code= !getPdtData(list_->file,list_->position,pdtData+i);
       if(ret_code) 
       { stat[i]=1;
         pdtData[i].beamP=list_->beamP;
         if(list_->posaux)
         { /*  printf(" N2=%d\n",list_->posaux); */
            ret_code= !getPdtData(list_->file,list_->posaux,pdtDataAux+i);
            if(ret_code)statAux[i]=1;
         }   
       }
       /* printf("\n");*/
       break;  
     } 
     list_=list_->next;
   }    
   if(list) delPdtList(list);
   if(ret_code) 
   {
     if(pdtData[i].alpha) { alphaMode=i+1; sf_alpha=&(alpha_pdt);}
     if(pdtData[i].pow1<0) *be=pdtData[i].pow1+1; else *be=1.;
     *mass=pdtData[i].mass;
   } 
   return ret_code; 
}

#define WIDTH 30
int m_pdt(int i)
{
   pdtList *list=NULL;
   pdtList *list_; 
   int k=0;
   char * strmen; 
   long pNum;

   i--;
   pinf_int(Nsub,i+1,NULL,&pNum);

   comphepPdtList(pNum, &list);
   if(!list) return 0;
   for(list_=list;list_;list_=list_->next,k++);
   strmen=malloc(2+WIDTH*k);
   strmen[0]=WIDTH;strmen[1]=0; 
   for(list_=list; list_; list_=list_->next) sprintf(strmen+strlen(strmen),
       " PDT:%-*.*s",WIDTH-5,WIDTH-5,list_->name);
   menu1(80-4-WIDTH,7,"",strmen,"",NULL,&k);
   if(k==0) { free(strmen);delPdtList(list);return 0;}  
   k--;   
   for(list_=list; k ; list_=list_->next, k--);
   if(pdtName[i]) free(pdtName[i]);
   pdtName[i]=malloc(5+strlen(list_->name));   
   sprintf(pdtName[i],"PDT:%s",list_->name); 

   free(strmen); delPdtList(list); return 1;
}

double c_pdt(int i, double x, double q)
{  double r;
   i--;
   r=interFunc( x, q , pdtData+i);
   if(statAux[i])
   {  double  sc2=0.221*0.221;
      double rAux=interFunc( x, q , pdtDataAux+i);
      r=r*(1-sc2)+rAux*sc2;
   }

   if(r<0)r=0;
   return r;
}
