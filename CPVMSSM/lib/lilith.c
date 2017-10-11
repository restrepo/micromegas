#include"../../sources/micromegas.h"
#include"../../sources/micromegas_aux.h"
#include "pmodel.h"

int LiLithF(char*fname)
{
  unsigned int i, npart=0;

  double CU, Cb, Cl,  CU5, Cb5, Cl5,  CV,Cgamma, Cg;
  double Mcp=findValW("Mcp"), Mbp=findValW("Mbp"), Mtp=findValW("Mtp");
  double vev = 2*findValW("MW")*findValW("SW")/findValW("EE");

  FILE*f;
  
  char *parts[3]={"h1","h2","h3"};
  int hPdgn[3]={25,35,36};
  

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
    char format[100];
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

      sprintf(format,"%%lf %%*f  3 %d 6 6", hPdgn[i]);
      CU=   slhaValFormat("LiLithInputHiggsCouplingsFermions",0.,format);
      sprintf(format,"%%*f %%lf  3 %d 6 6", hPdgn[i]);
      CU5=  slhaValFormat("LiLithInputHiggsCouplingsFermions",0.,format);

      sprintf(format,"%%lf %%*f  3 %d 5 5", hPdgn[i]);
      Cb=   slhaValFormat("LiLithInputHiggsCouplingsFermions",0.,format);
      sprintf(format,"%%*f %%lf  3 %d 5 5", hPdgn[i]);
      Cb5=  slhaValFormat("LiLithInputHiggsCouplingsFermions",0.,format);
    
      sprintf(format,"%%lf %%*f  3 %d 15 15", hPdgn[i]);
      Cl=   slhaValFormat("LiLithInputHiggsCouplingsFermions",0.,format);
      sprintf(format,"%%*f %%lf  3 %d 15 15", hPdgn[i]);
      Cl5=  slhaValFormat("LiLithInputHiggsCouplingsFermions",0.,format);
      
      sprintf(format,"%%lf  3 %d  24 24", hPdgn[i]);
      CV=   slhaValFormat("LiLithInputHiggsCouplingsBosons",0.,format);
      sprintf(format,"%%lf  3 %d  21 21", hPdgn[i]);
      Cg=   slhaValFormat("LiLithInputHiggsCouplingsBosons",0.,format);

//      sprintf(format,"%%lf  3 %d  22 22", hPdgn[i]);
//      Cgamma= slhaValFormat("LiLithInputHiggsCouplingsBosons",0.,format);

      char LaTxt[20];
      double LaV,LaV5,LaSM;
      sprintf(LaTxt,"LAA%s",parts[i]);
      LaV=findValW(LaTxt);
      sprintf(LaTxt,"imLAA%s",parts[i]);
      LaV5=findValW(LaTxt); 
             
      Cgamma=sqrt(LaV*LaV+LaV5*LaV5/4)/lAAhSM(mass,alphaQCD(mass)/M_PI, Mcp,Mbp,Mtp,vev);

    
      fprintf(f, "  <reducedcouplings part=\"%s\">\n", parts[i]);
      fprintf(f, "    <mass>%f</mass>\n", mass);
      fprintf(f, "    <C to=\"uu\" part=\"re\">  %f</C>\n", CU);
      fprintf(f, "    <C to=\"uu\" part=\"im\">  %f</C>\n", CU5);
      fprintf(f, "    <C to=\"bb\" part=\"re\">%f</C>\n", Cb);
      fprintf(f, "    <C to=\"bb\" part=\"im\">%f</C>\n", Cb5); 
      fprintf(f, "    <C to=\"mumu\" part=\"re\">%f</C>\n", Cl);
      fprintf(f, "    <C to=\"mumu\" part=\"im\">%f</C>\n", Cl5); 
      fprintf(f, "    <C to=\"tautau\" part=\"re\">%f</C>\n", Cl);
      fprintf(f, "    <C to=\"tautau\" part=\"im\">%f</C>\n", Cl5); 
    
      fprintf(f, "    <C to=\"VV\">%f</C>\n", CV);
      fprintf(f, "    <C to=\"gammagamma\">%f</C>\n", Cgamma);
      fprintf(f, "    <C to=\"gg\">%f</C>\n", Cg);
//    fprintf(f, "    <C to=\"Zgamma\">%f</C>\n", 1.);
      fprintf(f, "    <precision>%s</precision>\n", "LO");
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

