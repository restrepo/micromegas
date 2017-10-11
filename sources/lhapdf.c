#include <dlfcn.h>
#include"../sources/micromegas.h"
#include"../sources/micromegas_aux.h"

extern void   initpdfsetbynamem_(int*P,char*name,int len);
extern void   numberpdfm_(int *P,int*max);
extern void   initpdfm_(int*P,int*nSet);
extern void   evolvepdfm_(int*P,double* x, double*q, double*fxq);
extern double alphaspdfm_(int*,double*);       

static int P=1;
static double lhaPartonPdf(int p, double x,double q)
{
   double z,f[14];
   evolvepdfm_(&P,&x,&q,f);
   switch(p)
   {
     case 2 :           z=f[8]/x;  break;
     case 1 :           z=f[7]/x;  break;
     case 3 : case -3 : z=f[9]/x;  break;
     case 4 : case -4 : z=f[10]/x; break;
     case 5 : case -5 : z=f[11]/x; break;
     case 21: case -21: z=f[6]/x;  break;
     case -1:           z=f[5]/x;  break;
     case -2:           z=f[4]/x;  break;
     default:           z=0;
   }
   if(z<0) return 0;
   return z;                                 
}

static double lhaAlpha(double q) { return  alphaspdfm_(&P,&q);}


void setLHAPDF(int nset, char *name)
{   int nmax;
    initpdfsetbynamem_(&P,name,strlen(name));
    numberpdfm_(&P,&nmax);
    printf("\n %s   sets. members   0 - %d  are available\n",name,nmax);
    initpdfm_(&P,&nset); 
    parton_distr=lhaPartonPdf;
    parton_alpha=lhaAlpha;
    sprintf(pdfName,"LHA:%s:%d",name,nset);
}

void LHAPDFList(void)
{ int i;
  void (*getpdfsetlist)(char* s, size_t len)=NULL;
  void* h=dlopen(NULL,RTLD_NOW); 
  
  getpdfsetlist=dlsym(h,"lhapdf_getpdfsetlist_"); 
  if(!getpdfsetlist) getpdfsetlist=dlsym(h,"_lhapdf_getpdfsetlist_");
  if(getpdfsetlist)
  {  char buff[10000];
     char *ch1;
     getpdfsetlist(buff,10000);
     for(i=9999; i>=0 && buff[i]==' '; i--) ;
     buff[i+1]=0;
     i=strlen(buff);

     printf("List of LHAPDF sets\n"); 
     for(ch1=buff;ch1; ch1=strchr(ch1,' '))
     {  char name[50];
        while(ch1[0]==' ') ch1++;
        sscanf(ch1,"%s",name);
        printf("%s\n",name);
    }
  }

}

//  FORTRAN 

extern void setlhapdf_(int *nset, char *fname,int len);

void  setlhapdf_(int *nset, char *fname,int len)
{ 
  char cname[50];
  fName2c(fname,cname,len);
  setLHAPDF(*nset, cname);
}  

extern void lhapdflist_(void);

void lhapdflist_(void) { LHAPDFList(); }
