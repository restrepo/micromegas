#include "../include/micromegas.h"
#include "../include/micromegas_aux.h"

static double  fe(double E)  //  Fig.3   1506.03811
{
  double ln10E[40]= { 3.6990, 3.9315, 4.1640, 4.3965, 4.6291, 4.8616, 5.0941, 5.3267, 5.5592, 5.7917, 6.0242, 6.2568, 6.4893, 6.7218, 6.9543, 7.1869, 7.4194, 7.6519, 7.8844, 8.1170, 8.3495, 8.5820, 8.8145, 9.0471, 9.2796, 9.5121, 9.7446, 9.9772, 10.2097, 10.4422, 10.6747, 10.9073, 11.1398, 11.3723, 11.6048, 11.8374, 12.0699, 12.3024, 12.5349, 12.7675};
  double  e[40]=    { 0.3084, 0.3072, 0.3044, 0.2983, 0.2944, 0.2819, 0.2593, 0.2283, 0.1904, 0.1527, 0.1770, 0.3680, 0.6173, 0.8228, 0.9298, 0.9756, 0.9932, 0.9889, 0.9288, 0.8256, 0.7032, 0.5589, 0.4296, 0.3952, 0.4660, 0.5714, 0.6388, 0.6243,  0.5593,  0.4868,  0.4502,  0.4357,  0.4083,  0.3994,  0.4049,  0.4064,  0.4070,  0.4016,  0.4036,  0.4056};
  return polint3(log10(E), 40,ln10E,e);
}

static double fA(double E)   //  Fig.3  1506.03811
{
  double ln10E[40]= { 3.6990, 3.9315, 4.1640, 4.3965, 4.6291, 4.8616, 5.0941, 5.3267, 5.5592, 5.7917, 6.0242, 6.2568, 6.4893, 6.7218, 6.9543, 7.1869, 7.4194, 7.6519, 7.8844, 8.1170, 8.3495, 8.5820, 8.8145, 9.0471, 9.2796, 9.5121, 9.7446, 9.9772, 10.2097, 10.4422, 10.6747, 10.9073, 11.1398, 11.3723, 11.6048, 11.8374, 12.0699, 12.3024, 12.5349, 12.7675};
  double A[40]=     { 0.8767, 0.7697, 0.7021, 0.6842, 0.6963, 0.6748, 0.5851, 0.4591, 0.3256, 0.2161, 0.1448, 0.1631, 0.2924, 0.4681, 0.5957, 0.6781, 0.7180, 0.7260, 0.7195, 0.6858, 0.6288, 0.5462, 0.4469, 0.3591, 0.3184, 0.3378, 0.3905, 0.4228,  0.4212,  0.3902,  0.3307,  0.3901,  0.4332,  0.4148,  0.4029,  0.4047,  0.4069,  0.4017,  0.4037,  0.4056};
  return polint3(log10(E), 40,ln10E,A);
}


static double intEe(double E, double *tab) { return eSpectdNdE(E,tab)*fe(E*1E9);}
static double intEA(double E, double *tab) { return eSpectdNdE(E,tab)*fA(E*1E9);} 


double  pann(double Mdm,int pdg)
{ int err;
  double tab[NZ],r=0;
  basicSpectra(Mdm, pdg, 0, tab);
  r+=simpson_arg(intEA,tab, 5E-6,Mdm, 1E-3,&err)/(2*Mdm);
  basicSpectra(Mdm, pdg, 1, tab); 
  r+=simpson_arg(intEe,tab,1E-7*Mdm,Mdm, 1E-3,&err)/Mdm;
  return r;            
} 
     
double pannTxt(double Mdm,char*pdgTxt)
{ 
  int  pdg; 
  sscanf(pdgTxt,"%d",&pdg);
  pann(Mdm,pdg);
}


int main(void)
{
  displayPlot("f_eff, FIG.3","E",5E3,1E13,1,2,"electron",0,fe,NULL    
                                             ,"photon",  0,fA,NULL);

  double tabA[NZ],tabE[NZ];
  int pdg=21;
  double MM=10000;
  char txt[20];  
  basicSpectra(MM,pdg,0,tabA);
  basicSpectra(MM,pdg,1, tabE);
  sprintf(txt, "pdg=%d ->A/e",pdg);
  
  displayPlot(txt,"E",5E-6,MM,1,2,"A",0,intEA,tabA,"e",0,intEe,tabE);

  char WT_txt[10];
  sprintf(WT_txt,"%d",24+'T'); 
  SpectraFlag=2;
  displayPlot("Fig.4","Mdm", 5,1000,1,5
                                         , "e",0,pannTxt,"11"
                                         , "g",0,pannTxt,"21"
                                         , "A",0,pannTxt,"22"
                                         , "W",0,pannTxt,"24"
                                         , "W_T",0,pannTxt,WT_txt
                                         );

                                       
                                         
    
}