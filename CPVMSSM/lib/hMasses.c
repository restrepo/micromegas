#include"../../sources/micromegas.h"
#include"pmodel.h"


void o1Contents(FILE * f)
{ double ok=0; 
  fprintf(f,"\n~o1 = ");

  fprintf(f,"(%.3f%+.3f*i)*bino+",Zn11r(ok),Zn11i(ok)); 
  fprintf(f,"(%.3f%+.3f*i)*wino+",Zn12r(ok),Zn12i(ok)); 
  fprintf(f,"(%.3f%+.3f*i)*higgsino1+",Zn13r(ok),Zn13i(ok)); 
  fprintf(f,"(%.3f%+.3f*i)*higgsino2",Zn14r(ok),Zn14i(ok)); 
  fprintf(f,"\n");
}
