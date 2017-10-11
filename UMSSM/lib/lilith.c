#include"../../sources/micromegas.h"
#include"../../sources/micromegas_aux.h"
#include "pmodel.h"

int LiLithF(char*fname)
{
  unsigned int i, npart=0;

  double CU, Cb, Ctau, CV, Cgamma, Cg;

  char *parts[4]={"h1","h2","h3","ha"};
  FILE*f; 
  f=fopen(fname,"w");
  
  fprintf(f, "<?xml version=\"1.0\"?>\n");
  fprintf(f, "<lilithinput>\n");

  for(i=0; i<4; i++) 
  {
    double mass = pMass(parts[i]);
    if(mass < 123. || mass > 128.) {
      continue;
    }
    ++npart;

    // compute invisible and undetected branching ratios
    double invBR = 0., undBR = 0.;
    double w;
    txtList L;
    w=pWidth((char*)parts[i], &L);

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
            findBr(L, "W+ W-") - findBr(L, "A A") - findBr(L, "Z1 Z1") -
            findBr(L, "G G") - findBr(L, "m M") - findBr(L, "A Z1") -
            findBr(L, "u U") - findBr(L, "d D") - findBr(L, "s S");

    CU =   slhaVal("REDCOUP",0.,2,1+i,1);
    Cb =   slhaVal("REDCOUP",0.,2,1+i,3);
    Ctau = slhaVal("REDCOUP",0.,2,1+i,2);
    CV =   slhaVal("REDCOUP",0.,2,1+i,4);
    Cg=    slhaVal("REDCOUP",0.,2,1+i,5);
    Cgamma=slhaVal("REDCOUP",0.,2,1+i,6);

    fprintf(f, "  <reducedcouplings part=\"%s\">\n", parts[i]);
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
  }

  fprintf(f, "</lilithinput>\n");
  fclose(f);
  return npart;
}

