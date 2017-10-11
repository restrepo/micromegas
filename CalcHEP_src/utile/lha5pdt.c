
#include<stdio.h>
#include<math.h>
#include<string.h>

extern void   initpdfsetbynamem_(int*P,char*name,int len);
extern void   getlhapdfversion_(char * version, int len);
extern void   numberpdfm_(int *P,int*max);
extern void   initpdfm_(int*P,int*nSet);
extern void   evolvepdfm_(int*P,double* x, double*q, double*fxq);
extern double alphaspdfm_(int*,double*);       
extern void   getdatapath_(char* dirpath, int len);


#define NQ  50
#define NX 100
#define NP 9

double q_grid[NQ]={
 1.0000E+00,1.2068E+00,1.4563E+00,1.7575E+00,2.1210E+00,2.5595E+00,3.0888E+00,3.7276E+00,4.4984E+00,5.4287E+00
,6.5513E+00,7.9060E+00,9.5410E+00,1.1514E+01,1.3895E+01,1.6768E+01,2.0236E+01,2.4421E+01,2.9471E+01,3.5565E+01
,4.2919E+01,5.1795E+01,6.2506E+01,7.5431E+01,9.1030E+01,1.0985E+02,1.3257E+02,1.5999E+02,1.9307E+02,2.3300E+02
,2.8118E+02,3.3932E+02,4.0949E+02,4.9417E+02,5.9636E+02,7.1969E+02,8.6851E+02,1.0481E+03,1.2649E+03,1.5264E+03
,1.8421E+03,2.2230E+03,2.6827E+03,3.2375E+03,3.9069E+03,4.7149E+03,5.6899E+03,6.8665E+03,8.2864E+03,1.0000E+04
                  };

double x_grid[NX]={
 1.0000E-09,1.4508E-09,2.1049E-09,3.0539E-09,4.4306E-09,6.4281E-09,9.3260E-09,1.3530E-08,1.9630E-08,2.8480E-08
,4.1320E-08,5.9948E-08,8.6975E-08,1.2619E-07,1.8307E-07,2.6561E-07,3.8535E-07,5.5908E-07,8.1113E-07,1.1768E-06
,1.7074E-06,2.4771E-06,3.5938E-06,5.2140E-06,7.5646E-06,1.0975E-05,1.5923E-05,2.3101E-05,3.3516E-05,4.8626E-05
,7.0548E-05,1.0235E-04,1.4850E-04,2.1544E-04,3.1257E-04,4.5349E-04,6.5793E-04,9.5455E-04,1.3849E-03,2.0092E-03
,2.9151E-03,4.2292E-03,6.1359E-03,8.9022E-03,1.2915E-02,1.8738E-02,2.7186E-02,3.9442E-02,5.7224E-02,8.3022E-02
,1.0000E-01,1.1837E-01,1.3673E-01,1.5510E-01,1.7347E-01,1.9184E-01,2.1020E-01,2.2857E-01,2.4694E-01,2.6531E-01
,2.8367E-01,3.0204E-01,3.2041E-01,3.3878E-01,3.5714E-01,3.7551E-01,3.9388E-01,4.1224E-01,4.3061E-01,4.4898E-01
,4.6735E-01,4.8571E-01,5.0408E-01,5.2245E-01,5.4082E-01,5.5918E-01,5.7755E-01,5.9592E-01,6.1429E-01,6.3265E-01
,6.5102E-01,6.6939E-01,6.8776E-01,7.0612E-01,7.2449E-01,7.4286E-01,7.6122E-01,7.7959E-01,7.9796E-01,8.1633E-01
,8.3469E-01,8.5306E-01,8.7143E-01,8.8980E-01,9.0816E-01,9.2653E-01,9.4490E-01,9.6327E-01,9.8163E-01,1.0000E+00
                  };  
static int pdg[NP]={-2,-1,1,2,3,4,5,6,21};

static   char*pdf_name="NNPDF23_lo_as_0130.LHgrid";
static   int mem=0; 
static   int Index=247000;
static   char*reference="arXiv:1207.1303";
 
int main(int argc,char** argv)
{  
   static int P=1;        
   char pdf_name_[100],fname[100];
   if(mem==0) strcpy(pdf_name_,pdf_name); else  sprintf(pdf_name_,"%s:%d",pdf_name,mem);
   sprintf(fname,"%s.pdt",pdf_name_);
    
   initpdfsetbynamem_(&P,pdf_name,strlen(pdf_name));
   initpdfm_(&P,&mem);
   
int i,ix,iq;

  FILE*f=fopen(fname,"w");

  fprintf(f,"#distribution \"%s(proton)\"        2212 => ",pdf_name_);
  for(i=0;i<NP;i++) if( abs(pdg[i])>=3 && abs(pdg[i])<=6) fprintf(f," (%d %d) ", abs(pdg[i]),-abs(pdg[i])); else fprintf(f," %d ",pdg[i]);
  fprintf(f,"\n");
  fprintf(f,"#distribution \"%s(anti-proton)\"  -2212 => ",pdf_name_);
  for(i=0;i<NP;i++) if( abs(pdg[i])>=3 && abs(pdg[i])<=6) fprintf(f," (%d %d) ", abs(pdg[i]),-abs(pdg[i])); else 
                               {if(pdg[i]>20)   fprintf(f," %d ",pdg[i]); else fprintf(f," %d ",-pdg[i]);}
  fprintf(f,"\n");
  fprintf(f,"#Index %d\n",Index);
  fprintf(f,"#Set %d\n",mem);         
  fprintf(f,"#Interpolation biCubicLogXQ\n");
  fprintf(f,"#Reference %s\n",reference);
  fprintf(f,"#Source lhapdf5\n");
  
  fprintf(f,"#Q_grid\n");
  for(i=0;i<NQ;i++) { fprintf(f," %.4E",q_grid[i]); if((i+1)%10==0) fprintf(f,"\n");}
  fprintf(f,"\n#X_grid\n");
  for(i=0;i<NX;i++) { fprintf(f," %.4E",x_grid[i]); if((i+1)%10==0) fprintf(f,"\n");}
  
  fprintf(f,"\n#Alpha\n");
  for(i=0;i<NQ;i++) { fprintf(f," %.4E",alphaspdfm_(&P,q_grid+i)); if((i+1)%10==0) fprintf(f,"\n");}

  
for(i=0;i<9;i++)
{ fprintf(f,"\n#%d-parton\n",i+1);


  for(iq=0;iq<NQ;iq++)
  {  for(ix=0;ix<NX;ix++)
     { double ff[14],z;
       evolvepdfm_(&P,x_grid+ix,q_grid+iq,ff);
       switch (pdg[i])
       {
         case 2 :           z=ff[8];  break;
         case 1 :           z=ff[7];  break;
         case 3 : case -3 : z=ff[9];  break;
         case 4 : case -4 : z=ff[10]; break;
         case 5 : case -5 : z=ff[11]; break;
         case 21: case -21: z=ff[6];  break;
         case -1:           z=ff[5];  break;
         case -2:           z=ff[4];  break;
         default:           z=0;
       }
       if(z<0) z=0; 
       fprintf(f," %.4E",z);
     } 
     fprintf(f,"\n");
  }     
}
  fclose(f);   
  return 0;
}
