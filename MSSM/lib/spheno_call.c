
/*===================================*/
#include"pmodel.h"
#include"pmodel_aux.h"
#include"../../include/micromegas_aux.h"
#include<sys/wait.h>
#include<unistd.h>

#define FIN  "LesHouches.in"
#define FOUT "SPheno.spc"

#define VERSION "4.0.4" 

static int SystemCall(int mode)
{ 
  char buff[2000];
  int err;
  FILE*f;

  if(!access(FOUT,R_OK)) unlink(FOUT);

  sprintf(buff,"%s/Packages/SPheno-%s/bin/SPheno",micrO,VERSION);
  if(access( buff,X_OK))
  { sprintf(buff,"cd %s/Packages;  make -f SPHENO.makef  VERSION=%s",micrO,VERSION);
    system(buff);
    sprintf(buff,"%s/Packages/SPheno-%s/bin/SPheno",micrO,VERSION);
    if(access( buff,X_OK))
    {
      printf("Executable \n %s\n is not found. Program stops.\n",buff);
      exit(13);
    }  
  }  
/*
   f=fopen("Control.in","w");
   fprintf(f,"0       | ErrorLevel\n"
             ".false. ! Calculation of branching ratios\n"
             ".false. ! Calculation of cross sections\n");
   fclose(f);
*/
  err=system(buff);   
  if(err>=0)   err=slhaRead(FOUT,0); else cleanSLHAdata();
  return err;
}


double  sphenoSUGRAc(double tb, double gMG1,double gMG2,double gMG3,
             double gAl, double gAt, double gAb, double sgn, double gMHu, double gMHd,
             double gMl2,double gMl3,double gMr2,double gMr3,
             double gMq2,double gMq3,double gMu2,double gMu3,double gMd2,double gMd3)
{ 
   int err;

   err=sugraLesH(FIN, tb, gMG1,gMG2,gMG3,
             gAl, gAt, gAb, sgn, gMHu, gMHd,
             gMl2,gMl3,gMr2,gMr3,
             gMq2,gMq3,gMu2,gMu3,gMd2,gMd3);
              
   if(err) {printf("can not write down LesHouches.in file\n"); exit(10);}

   return SystemCall(1);
}

double  sphenoSUGRAnuhc(double tb, double gMG1,double gMG2,double gMG3,
             double gAl, double gAt, double gAb,double gMl2,double gMl3,double gMr2,double gMr3,
             double gMq2,double gMq3,double gMu2,double gMu3,double gMd2,double gMd3,double mu,double MA)
{ 
   int err=sugraHiggsLesH(FIN, tb, gMG1,gMG2,gMG3, gAl, gAt, gAb, gMl2,gMl3,gMr2,gMr3,gMq2,gMq3,gMu2,gMu3,gMd2,gMd3,mu,MA);
   FILE*f=fopen(FIN,"a");
   
   fprintf(f," 0  -1  # EWSB\n");
   fclose(f);  
  
              
   if(err) {printf("can not write down LesHouches.in file\n"); exit(10);}

   return SystemCall(1);
}


double  sphenoAMSBc(double m0,double m32, double tb, double sgn)
{ 
   int err;
   err=amsbLesH(FIN, m0,m32, tb, (int)sgn);
   if(err) { printf("can not write down LesHouches.in file\n"); exit(10);}

   return SystemCall(1);
}



double  sphenoEwsbMSSMc(double tb, double MG1, double MG2, double MG3, double Al, double At, double Ab,
  double mu, double MH3, double Ml1, double Ml2, double Ml3, double Mr1, double Mr2, double Mr3,
    double Mq1, double Mq2, double Mq3, double Mu1, double Mu2, double Mu3,
                                        double Md1, double Md2, double Md3)
{  int err;
   if(EWSBLesH(FIN,tb,MG1,MG2,MG3,Al,At,Ab,mu,MH3,Ml1,Ml2,Ml3,Mr1,Mr2,Mr3, Mq1,Mq2,Mq3,Mu1,Mu2,Mu3,
                                        Md1,Md2,Md3))
   { printf("can not write down LesHouches.in file\n"); exit(10);}
     err=SystemCall(0);
   return err;
}

