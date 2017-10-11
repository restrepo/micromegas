#include"../../sources/micromegas.h"
#include"../../sources/micromegas_aux.h"
#include "pmodel.h"

int LiLithF(char*fname)
{
  int i, npart=0;

  double CU, Cb, Ctau, CV, Cgamma, Cg;
  double Mcp=findValW("Mcp"), Mbp=findValW("Mbp"), Mtp=findValW("Mtp");
  double LGGSM, LAASM;
  double vev = 2*findValW("MW")*findValW("SW")/findValW("EE");
  FILE*f;
  double w;
  txtList L;
 

  double mass = pMass("h");

  if(mass < 123. || mass > 128.) return 0;

  f=fopen(fname,"w");
  
  fprintf(f, "<?xml version=\"1.0\"?>\n");
  fprintf(f, "<lilithinput>\n");
    // compute invisible and undetected branching ratios
    double invBR = 0., undBR = 0.;
    w=pWidth("h",&L);

    if(Mcdm1 < 0.5*mass) {
      char invdecay[50];
      char cdmName[50];
//      sortOddParticles(cdmName);
      strcpy(invdecay, CDM1);
      strcat(invdecay, ",");
      strcat(invdecay, CDM1);
      invBR = findBr(L, invdecay);
    }
    undBR = 1 - invBR - findBr(L, "b B") - findBr(L, "c C") - findBr(L, "l L") -
            findBr(L, "W+ W-") - findBr(L, "A A") - findBr(L, "Z Z") -
            findBr(L, "G G") - findBr(L, "m M") - findBr(L, "A Z") -
            findBr(L, "u U") - findBr(L, "d D") - findBr(L, "s S");

    LAASM=lAAhSM(mass, alphaQCD(mass)/M_PI, Mcp, Mbp, Mtp, vev);

    CU = 1;
    Cb = 1;
    Ctau=1;
    CV = 1;
    Cgamma = findValW("LAAH")/LAASM;
    Cg =1;
      
    fprintf(f, "  <reducedcouplings part=\"%s\">\n", "H");
    fprintf(f, "    <mass>%f</mass>\n", mass);
    fprintf(f, "    <C to=\"uu\">%f</C>\n", CU);
    fprintf(f, "    <C to=\"bb\">%f</C>\n", Cb);
    fprintf(f, "    <C to=\"mumu\">%f</C>\n", Ctau);
    fprintf(f, "    <C to=\"tautau\">%f</C>\n", Ctau);
    fprintf(f, "    <C to=\"VV\">%f</C>\n", CV);
    fprintf(f, "    <C to=\"gammagamma\">%f</C>\n", Cgamma);
    fprintf(f, "    <C to=\"gg\">%f</C>\n", Cg);
//    fprintf(f, "    <C to=\"Zgamma\">%f</C>\n", 1.);
    fprintf(f, "    <precision>%s</precision>\n", "BEST-QCD");
    fprintf(f, "    <extraBR>\n");
    fprintf(f, "      <BR to=\"invisible\">%f</BR>\n", invBR);
    fprintf(f, "      <BR to=\"undetected\">%f</BR>\n", undBR);
    fprintf(f, "    </extraBR>\n");
    fprintf(f, "  </reducedcouplings>\n");


  fprintf(f, "</lilithinput>\n");
  fclose(f);
  return 1;
}

