#ifndef __CPSH__
#define __CPSH__

#ifdef __cplusplus
extern "C" {
#endif

#include<stdio.h>

#include"../../CalcHEP_src/c_source/SLHAplus/include/SLHAplus.h"

extern int delFiles;

extern void HiggsMasses(FILE *f);
extern void o1Contents(FILE * f);
 
extern double
la1(double)  ,la2(double)  ,la3(double)  ,la4(double)  ,la5r(double) ,la5i(double),
la6r(double) ,la6i(double) ,la7r(double) ,la7i(double) ,MHc(double)  ,Mh1(double),
Mh2(double)  ,Mh3(double)  ,Zh11(double) ,Zh12(double) ,Zh13(double) ,Zh21(double),
Zh22(double) ,Zh23(double) ,Zh31(double) ,Zh32(double) ,Zh33(double) ,MSt1(double),
MSt2(double) ,Zt11r(double),Zt12r(double),Zt21r(double),Zt22r(double),Zt11i(double),
Zt12i(double),Zt21i(double),Zt22i(double),MSb1(double) ,MSb2(double) ,Zb11r(double),
Zb12r(double),Zb21r(double),Zb22r(double),Zb11i(double),Zb12i(double),Zb21i(double),
Zb22i(double),MSnl(double) ,MSl1(double) ,MSl2(double) ,Zl11r(double),Zl12r(double),
Zl21r(double),Zl22r(double),Zl11i(double),Zl12i(double),Zl21i(double),Zl22i(double),
MNE1(double) ,MNE2(double) ,MNE3(double) ,MNE4(double) ,Zn11r(double),Zn11i(double),
Zn12r(double),Zn12i(double),Zn13r(double),Zn13i(double),Zn14r(double),Zn14i(double),
Zn21r(double),Zn21i(double),Zn22r(double),Zn22i(double),Zn23r(double),Zn23i(double),
Zn24r(double),Zn24i(double),Zn31r(double),Zn31i(double),Zn32r(double),Zn32i(double),
Zn33r(double),Zn33i(double),Zn34r(double),Zn34i(double),Zn41r(double),Zn41i(double),
Zn42r(double),Zn42i(double),Zn43r(double),Zn43i(double),Zn44r(double),Zn44i(double),
MC1(double)  ,MC2(double)  ,Zu11r(double),Zu11i(double),Zu12r(double),Zu12i(double),
Zu21r(double),Zu21i(double),Zu22r(double),Zu22i(double),Zv11r(double),Zv11i(double),
Zv12r(double),Zv12i(double),Zv21r(double),Zv21i(double),Zv22r(double),Zv22i(double),
Tu3r(double) ,Tu3r(double), Td3r(double) ,Td3r(double); 


extern double cpHiggs(double,double,double,double,double,double,double,double,
 double,double,double,double,double,double,double,double,double,double,double,
 double,double,double,double,double,double,double,double,double,double,double,
 double,double);


extern double  FevenFr(double),FevenFi(double),Foddr(double),Foddi(double),
F0r(double),F0i(double),mbPole(double);

extern double MbPole;
extern void edm_(double* , double *);
extern int readVarCPVMSSM(char*fname);

extern double bsmumu(void);
extern double bsgnlo(void);   
extern double Bulnu(void); 
extern double Bdll(void);  
extern double ABsg(void);  
extern double deltaMd(void);
extern double deltaMs(void);
extern double EDMth(void) ;
extern double EDMel(void) ;
extern double EDMmu(void) ;
extern double EDMHg(void);
extern double gmuon(void);
extern int LepBound(double * ratio);

#ifdef __cplusplus
}
#endif 

#endif
