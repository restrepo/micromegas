
#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"

static double IC22nuArm2(double E) { return 1E6*IC22nuAr(E);}
static double IC22nuBarArm2(double E) { return 1E6*IC22nuBarAr(E);}


int main(void)
{
  displayPlot("Effective Area[m^2]", "Enu", 10, 1000, 1,4
                                            ,"icDCnu",   0, icDCnuAr,NULL
                                             ,"icDCnuBar",0,icDCnuBarAr,NULL
                                             ,"ic22nu",   0, IC22nuArm2,NULL 
                                             ,"ic22nuBar",0, IC22nuBarArm2,NULL 
                                            );
 

}