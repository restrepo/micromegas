#include"../../sources/micromegas.h"
#include"../../sources/micromegas_aux.h"
#include "pmodel.h"

int LiLithF(char*fname)
{
  unsigned int i, npart=0;

  double tb, sb, cb, alpha, sa, ca, ta, samb, camb, dMb, MbHl, MbH, MbH3, MbSM;
  double CU, Cb, Ctau, CV, Cgamma, Cg;
  double Mcp=findValW("Mcp"), Mbp=findValW("Mbp"), Mtp=findValW("Mtp");
  double LGGSM, LAASM;
  double vev = 2*findValW("MW")*findValW("SW")/findValW("EE");
  FILE*f;
  

  tb=findValW("tB");
  sb=tb/sqrt(1+tb*tb);
  cb=1/sqrt(1+tb*tb);
  alpha=findValW("alpha");
  sa=sin(alpha);
  ca=cos(alpha);
  ta=sa/ca;
  samb=sa*cb-ca*sb;
  camb=ca*cb+sa*sb;
  dMb=findValW("dMb");
  MbSM=findValW("Mb");
  MbHl = MbSM/(1+dMb)*(1-dMb/ta/tb);
  MbH = MbSM/(1+dMb)*(1+dMb*ta/tb);  
  MbH3 = MbSM/(1+dMb)*(1-dMb/tb/tb);

  // define Higgs states possibly contributing to the signal
  char *parts[3]={"h","H","H3"};

  f=fopen(fname,"w");
  
  fprintf(f, "<?xml version=\"1.0\"?>\n");
  fprintf(f, "<lilithinput>\n");

  for(i=0; i<3; i++) {
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
            findBr(L, "W+ W-") - findBr(L, "A A") - findBr(L, "Z Z") -
            findBr(L, "G G") - findBr(L, "m M") - findBr(L, "A Z") -
            findBr(L, "u U") - findBr(L, "d D") - findBr(L, "s S");

    LGGSM=lGGhSM(mass, alphaQCD(mass)/M_PI, Mcp, Mbp, Mtp, vev);
    LAASM=lAAhSM(mass, alphaQCD(mass)/M_PI, Mcp, Mbp, Mtp, vev);

    if(strcmp(parts[i], "h") == 0) {
      CU = ca/sb;
      Cb = -(sa/cb)*(MbHl/MbSM);
      Ctau = -(sa/cb);
      CV = -samb;
      Cgamma = findValW("LAAh")/LAASM;
      Cg = findValW("LGGh")/LGGSM;
    } else if(strcmp(parts[i], "H") == 0) {
      CU = sa/sb;
      Cb = (ca/cb)*(MbH/MbSM);
      Ctau = ca/cb;
      CV = camb;
      Cgamma = findValW("LAAH")/LAASM;
      Cg = findValW("LGGH")/LGGSM;
    } else { // for H3 (i.e. A)
      CU = 1./tb;
      Cb = tb*(MbH3/MbSM);
      Ctau = tb;
      CV = 0.;
      Cgamma = 0.5*findValW("LAAH3")/LAASM;
      Cg = 0.5*findValW("LGGH3")/LGGSM;
    }


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

