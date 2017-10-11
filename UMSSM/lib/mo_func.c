#include "pmodel.h"
#include "../../include/micromegas.h"
#include "../../include/micromegas_aux.h"
/*#include"../../CalcHEP_src/c_source/model_aux/include/SLHAreader.h"*/

#include "lpath.h"
#include <math.h>

void o1Contents(FILE * f)
{
  fprintf(f,"\n~o1 = %.3f*bino %+.3f*wino %+.3f*higgsino1 %+.3f*higgsino2 %+.3f*singlino %+.3f*binoprime\n",
      findValW("Zn11"), findValW("Zn12"), findValW("Zn13"), findValW("Zn14"), findValW("Zn15"), findValW("Zn16"));
}

double SignParam (double x)
{
  if (x>=0) return 1;
  else return -1;
}

//Random number generator : first don't forget to initialize the pseudo random number generator that rand() uses, for instance with srand(time(NULL));
//randpar("x=name of a given variable of the model (or not)",x_start,x_max,precision,sign,logarithmic scan,experimental contraint);
//examples : M1p=randpar("M1p",0,20000,2,1,0,0); Mtp=randpar("Mtp",173.34,1,3,1,0,1); msd=randpar("",27.5,0.3,3,1);

//caution when scanning on mu and Alda with negative values allowed (flag sign=1) :
//start using randpar for mu and then use it for Alda since they must get the same sign.

double randpar(char * xx, double start, double amp, int pwr, int sign, int logmod, int expval)
{
double XX;double sgn=0;int FloorAmp;

if (sign==0) sgn=1;
else if (sign==1) { while(sgn==0){sgn=(rand() % 3) - 1;} }
else {printf("**** ERROR **** Fifth must be either 0 or 1.\n");exit(1);}
if (!strcmp("Alda",xx)) {if (findValW("mu")>=0) {sgn=1;} else {sgn=-1;} }


if (logmod==0) {
FloorAmp=floor(amp*pow(10,pwr));
if (expval==0) {XX= sgn*( start + (rand() % (abs(FloorAmp)+1))*pow(10,-pwr) );}
else if (expval==1) {XX= start + sgn*( (rand() % (abs(FloorAmp)+1))*pow(10,-pwr) );}
else {printf("**** ERROR **** Seventh argument must be either 0 or 1.\n");exit(1);}
	       }
else if (logmod==1) {
FloorAmp=(start-amp)*pow(10,pwr);
XX= sgn*pow(10,-(amp + (rand() % (abs(FloorAmp)+1))*pow(10,-pwr)) );
                    }
else {printf("**** ERROR **** Sixth argument must be either 0 or 1.\n");exit(1);}

if(strcmp("",xx)) {assignValW(xx,XX);printf("%s = %+5e\n",xx,XX);}

return XX;
}

static int runTools(char * cmd, char * fout,int mode)
{  char * command;
   int err;
   double maxl;
   
   if(access(fout,F_OK)==0) unlink(fout);
   command=malloc(100+strlen(LPATH));

   sprintf(command,
    "lPath=%s;EXPCON_PATH=$lPath/UMSSMTools/EXPCON;export EXPCON_PATH;$lPath/%s"
    ,LPATH,cmd);printf("\n*****************************************\n");
    printf("*********** UMSSMTools_1.1.2 ************\n");
    printf("***** Last update : 4 november 2016 *****\n");
   err=System(command);printf("*****************************************\n"); 
   free(command);
   if(err>=0) err=slhaRead(fout,0);

   return err; 
}



/* =====  end of header part ========= */


#define V(N) findValW(#N)
  int umssmtools(int PDG_LSP)
{

   int err,nw;
   FILE*  f=fopen("UMSSM_inp.dat","w");
   if(f==NULL) return -1;

   fprintf(f,"Block MODSEL    # Select model\n"
             "  1    0           # EWSB input\n"
             "  3    1           # UMSSM PARTICLE CONTENT\n"
           "Block SMINPUTS    # Standard Model inputs\n");
   fprintf(f," 1   %.8E       # alpha_em^(-1)(MZ) SM MSbar\n",1/V(alfEMZ));
   fprintf(f," 2   %.8E       # G_Fermi \n",1.16637E-5); 
   fprintf(f," 3   %.8E       # alpha_s(MZ) SM MSbar\n",V(alfSMZ));
   fprintf(f," 4   %.8E       # MZ1\n",V(MZ1));
   fprintf(f," 5   %.8E       # mb(mb) SM MSbar\n",V(MbMb));
   fprintf(f," 6   %.8E       # mtop (pole mass)\n",V(Mtp));
   fprintf(f," 7   %.8E       # Mtau     \n",V(Ml));

   fprintf(f,"Block EXTPAR\n");
   fprintf(f," 0   -1        # EWSB\n"); 
   fprintf(f," 1   %.8E      # MG1\n",V(MG1));
   fprintf(f," 2   %.8E      # MG2\n",V(MG2));
   fprintf(f," 3   %.8E      # MG3\n",V(MG3));

   fprintf(f," 11  %.8E      # At \n",V(At));
   fprintf(f," 12  %.8E      # Ab \n",V(Ab));
   fprintf(f," 13  %.8E      # Atau\n",V(Al));
   fprintf(f," 31  %.8E      # Ml1\n",V(Ml2));
   fprintf(f," 32  %.8E      # Ml2\n",V(Ml2));
   fprintf(f," 33  %.8E      # Ml3\n",V(Ml3));
   fprintf(f," 34  %.8E      # MR2\n",V(Mr2));
   fprintf(f," 35  %.8E      # MR2\n",V(Mr2));
   fprintf(f," 36  %.8E      # MR3\n",V(Mr3));

   fprintf(f," 41  %.8E      # Mq1\n",V(Mq2));
   fprintf(f," 42  %.8E      # Mq2\n",V(Mq2));
   fprintf(f," 43  %.8E      # Mq3\n",V(Mq3));
   fprintf(f," 44  %.8E      # Mu1\n",V(Mu2));
   fprintf(f," 45  %.8E      # Mu2\n",V(Mu2));
   fprintf(f," 46  %.8E      # Mu3\n",V(Mu3));
   fprintf(f," 47  %.8E      # Md1\n",V(Md2));
   fprintf(f," 48  %.8E      # Md2\n",V(Md2));
   fprintf(f," 49  %.8E      # Md3\n",V(Md3));
   fprintf(f," 57  %.8E      # Mn1\n",V(Mn2));
   fprintf(f," 58  %.8E      # Mn2\n",V(Mn2));
   fprintf(f," 59  %.8E      # Mnlr\n",V(Mnlr));
   fprintf(f," 63  %.8E       # Alda\n",V(Alda));
   fprintf(f," 65  %.8E      # mu\n", V(mu));
   fprintf(f," 101 %.8E      # MZ2\n",V(MZ2));
   fprintf(f," 102 %.8E      # aZZ\n",V(aZZ));
   fprintf(f," 103 %.8E      # MK\n",V(MK));
   fprintf(f," 104 %.8E      # M1p\n",V(M1p));
   fprintf(f," 105 %.8E      # tE6\n",V(tE6));

   
   fclose(f);


   FILE*  ff=fopen("UMSSMTools.par","w");
   if(ff==NULL) return -1;
   fprintf(ff,"MG1     %+.9e\n",V(MG1));
   fprintf(ff,"MG2     %+.9e\n",V(MG2));
   fprintf(ff,"MG3     %+.9e\n",V(MG3));
   fprintf(ff,"M1p     %+.9e\n",V(M1p));
   fprintf(ff,"MK      %+.9e\n",V(MK));
   fprintf(ff,"Ml2     %+.9e\n",V(Ml2));
   fprintf(ff,"Ml3     %+.9e\n",V(Ml3));
   fprintf(ff,"Mr2     %+.9e\n",V(Mr2));
   fprintf(ff,"Mr3     %+.9e\n",V(Mr3));
   fprintf(ff,"Mq2     %+.9e\n",V(Mq2));
   fprintf(ff,"Mq3     %+.9e\n",V(Mq3));
   fprintf(ff,"Mu2     %+.9e\n",V(Mu2));
   fprintf(ff,"Mu3     %+.9e\n",V(Mu3));
   fprintf(ff,"Md2     %+.9e\n",V(Md2));
   fprintf(ff,"Md3     %+.9e\n",V(Md3));
   fprintf(ff,"Mn2     %+.9e\n",V(Mn2));
   fprintf(ff,"Mn32    %+.9e\n",V(Mn32));
   fprintf(ff,"At      %+.9e\n",V(At));
   fprintf(ff,"Ab      %+.9e\n",V(Ab));
   fprintf(ff,"Al      %+.9e\n",V(Al));
   fprintf(ff,"Am      %+.9e\n",V(Am));
   fprintf(ff,"SW      %+.9e\n",V(SW));
   fprintf(ff,"g2      %+.9e\n",V(g2));
   fprintf(ff,"gp      %+.9e\n",V(gp));
   fprintf(ff,"g1p     %+.9e\n",V(g1p));
   fprintf(ff,"v       %+.9e\n",V(v));
   fprintf(ff,"vs      %+.9e\n",V(vs));
   fprintf(ff,"NCp     %+.9e\n",V(NCp));
   fprintf(ff,"Q1      %+.9e\n",V(Q1));
   fprintf(ff,"Q2      %+.9e\n",V(Q2));
   fprintf(ff,"QQ      %+.9e\n",V(QQ));
   fprintf(ff,"QL      %+.9e\n",V(QL));
   fprintf(ff,"MW      %+.9e\n",V(MW));
   fprintf(ff,"MZ1     %+.9e\n",V(MZ1));
   fprintf(ff,"MZ2     %+.9e\n",V(MZ2));
   fprintf(ff,"aZZ     %+.9e\n",V(aZZ));
   fprintf(ff,"saZZ    %+.9e\n",V(saZZ));
   fprintf(ff,"caZZ    %+.9e\n",V(caZZ));
   fprintf(ff,"mu      %+.9e\n",V(mu));
   fprintf(ff,"Alda    %+.9e\n",V(Alda));
   fprintf(ff,"tE6     %+.9e\n",V(tE6));
   fprintf(ff,"Mtp     %+.9e\n",V(Mtp));
   fprintf(ff,"MbMb    %+.9e\n",V(MbMb));
   fprintf(ff,"Ms2GeV  %+.9e\n",V(Ms2GeV));
   fprintf(ff,"tb      %+.9e\n",V(tb));
   fprintf(ff,"lda     %+.9e\n",V(lda));
   fprintf(ff,"ZA11    %+.9e\n",V(ZA11));
   fprintf(ff,"ZA12    %+.9e\n",V(ZA12));
   fprintf(ff,"ZA13    %+.9e\n",V(ZA13));
   fprintf(ff,"ZA21    %+.9e\n",V(ZA21));
   fprintf(ff,"ZA22    %+.9e\n",V(ZA22));
   fprintf(ff,"ZA23    %+.9e\n",V(ZA23));
   fprintf(ff,"ZA31    %+.9e\n",V(ZA31));
   fprintf(ff,"ZA32    %+.9e\n",V(ZA32));
   fprintf(ff,"ZA33    %+.9e\n",V(ZA33));
   fprintf(ff,"Mha     %+.9e\n",V(Mha));
   fprintf(ff,"Zh11    %+.9e\n",V(Zh11));
   fprintf(ff,"Zh12    %+.9e\n",V(Zh12));
   fprintf(ff,"Zh13    %+.9e\n",V(Zh13));
   fprintf(ff,"Zh21    %+.9e\n",V(Zh21));
   fprintf(ff,"Zh22    %+.9e\n",V(Zh22));
   fprintf(ff,"Zh23    %+.9e\n",V(Zh23));
   fprintf(ff,"Zh31    %+.9e\n",V(Zh31));
   fprintf(ff,"Zh32    %+.9e\n",V(Zh32));
   fprintf(ff,"Zh33    %+.9e\n",V(Zh33));
   fprintf(ff,"Mh1     %+.9e\n",V(Mh1));
   fprintf(ff,"Mh2     %+.9e\n",V(Mh2));
   fprintf(ff,"Mh3     %+.9e\n",V(Mh3));
   fprintf(ff,"MHc     %+.9e\n",V(MHc));
   fprintf(ff,"Zu11    %+.9e\n",V(Zu11));
   fprintf(ff,"Zu12    %+.9e\n",V(Zu12));
   fprintf(ff,"Zu21    %+.9e\n",V(Zu21));
   fprintf(ff,"Zu22    %+.9e\n",V(Zu22));
   fprintf(ff,"Zv11    %+.9e\n",V(Zv11));
   fprintf(ff,"Zv12    %+.9e\n",V(Zv12));
   fprintf(ff,"Zv21    %+.9e\n",V(Zv21));
   fprintf(ff,"Zv22    %+.9e\n",V(Zv22));
   fprintf(ff,"MC1     %+.9e\n",V(MC1));
   fprintf(ff,"MC2     %+.9e\n",V(MC2));
   fprintf(ff,"MSuL    %+.9e\n",V(MSuL));
   fprintf(ff,"MSuR    %+.9e\n",V(MSuR));
   fprintf(ff,"MSdL    %+.9e\n",V(MSdL));
   fprintf(ff,"MSdR    %+.9e\n",V(MSdR));
   fprintf(ff,"Zt11    %+.9e\n",V(Zt11));
   fprintf(ff,"Zt12    %+.9e\n",V(Zt12));
   fprintf(ff,"MSt1    %+.9e\n",V(MSt1));
   fprintf(ff,"MSt2    %+.9e\n",V(MSt2));
   fprintf(ff,"Zb11    %+.9e\n",V(Zb11));
   fprintf(ff,"Zb12    %+.9e\n",V(Zb12));
   fprintf(ff,"MSb1    %+.9e\n",V(MSb1));
   fprintf(ff,"MSb2    %+.9e\n",V(MSb2));
   fprintf(ff,"MSeL    %+.9e\n",V(MSeL));
   fprintf(ff,"MSeR    %+.9e\n",V(MSeR));
   fprintf(ff,"MSmL    %+.9e\n",V(MSmL));
   fprintf(ff,"MSmR    %+.9e\n",V(MSmR));
   fprintf(ff,"Zl11    %+.9e\n",V(Zl11));
   fprintf(ff,"Zl12    %+.9e\n",V(Zl12));
   fprintf(ff,"MSl1    %+.9e\n",V(MSl1));
   fprintf(ff,"MSl2    %+.9e\n",V(MSl2));
   fprintf(ff,"MSne    %+.9e\n",V(MSne));
   fprintf(ff,"MSnm    %+.9e\n",V(MSnm));
   fprintf(ff,"MSnl    %+.9e\n",V(MSnl));
   fprintf(ff,"Mner    %+.9e\n",V(Mner));
   fprintf(ff,"Mnmr    %+.9e\n",V(Mnmr));
   fprintf(ff,"Mnlr    %+.9e\n",V(Mnlr));
   fprintf(ff,"Zn11    %+.9e\n",V(Zn11));
   fprintf(ff,"Zn12    %+.9e\n",V(Zn12));
   fprintf(ff,"Zn13    %+.9e\n",V(Zn13));
   fprintf(ff,"Zn14    %+.9e\n",V(Zn14));
   fprintf(ff,"Zn15    %+.9e\n",V(Zn15));
   fprintf(ff,"Zn16    %+.9e\n",V(Zn16));
   fprintf(ff,"Zn21    %+.9e\n",V(Zn21));
   fprintf(ff,"Zn22    %+.9e\n",V(Zn22));
   fprintf(ff,"Zn23    %+.9e\n",V(Zn23));
   fprintf(ff,"Zn24    %+.9e\n",V(Zn24));
   fprintf(ff,"Zn25    %+.9e\n",V(Zn25));
   fprintf(ff,"Zn26    %+.9e\n",V(Zn26));
   fprintf(ff,"Zn31    %+.9e\n",V(Zn31));
   fprintf(ff,"Zn32    %+.9e\n",V(Zn32));
   fprintf(ff,"Zn33    %+.9e\n",V(Zn33));
   fprintf(ff,"Zn34    %+.9e\n",V(Zn34));
   fprintf(ff,"Zn35    %+.9e\n",V(Zn35));
   fprintf(ff,"Zn36    %+.9e\n",V(Zn36));
   fprintf(ff,"Zn41    %+.9e\n",V(Zn41));
   fprintf(ff,"Zn42    %+.9e\n",V(Zn42));
   fprintf(ff,"Zn43    %+.9e\n",V(Zn43));
   fprintf(ff,"Zn44    %+.9e\n",V(Zn44));
   fprintf(ff,"Zn45    %+.9e\n",V(Zn45));
   fprintf(ff,"Zn46    %+.9e\n",V(Zn46));
   fprintf(ff,"Zn51    %+.9e\n",V(Zn51));
   fprintf(ff,"Zn52    %+.9e\n",V(Zn52));
   fprintf(ff,"Zn53    %+.9e\n",V(Zn53));
   fprintf(ff,"Zn54    %+.9e\n",V(Zn54));
   fprintf(ff,"Zn55    %+.9e\n",V(Zn55));
   fprintf(ff,"Zn56    %+.9e\n",V(Zn56));
   fprintf(ff,"Zn61    %+.9e\n",V(Zn61));
   fprintf(ff,"Zn62    %+.9e\n",V(Zn62));
   fprintf(ff,"Zn63    %+.9e\n",V(Zn63));
   fprintf(ff,"Zn64    %+.9e\n",V(Zn64));
   fprintf(ff,"Zn65    %+.9e\n",V(Zn65));
   fprintf(ff,"Zn66    %+.9e\n",V(Zn66));
   fprintf(ff,"MNE1    %+.9e\n",V(MNE1));
   fprintf(ff,"MNE2    %+.9e\n",V(MNE2));
   fprintf(ff,"MNE3    %+.9e\n",V(MNE3));
   fprintf(ff,"MNE4    %+.9e\n",V(MNE4));
   fprintf(ff,"MNE5    %+.9e\n",V(MNE5));
   fprintf(ff,"MNE6    %+.9e\n",V(MNE6));
   fprintf(ff,"MSG     %+.9e\n",V(MSG));
/* for h -> hh : */
   fprintf(ff,"h2toh11 %+.9e\n",V(h2toh11));
   fprintf(ff,"h3toh11 %+.9e\n",V(h3toh11));
   fprintf(ff,"h3toh12 %+.9e\n",V(h3toh12));
   fprintf(ff,"h3toh22 %+.9e\n",V(h3toh22));
/* for h -> aa : */
   fprintf(ff,"h1tohaa %+.9e\n",V(h1tohaa));
   fprintf(ff,"h2tohaa %+.9e\n",V(h2tohaa));
   fprintf(ff,"h3tohaa %+.9e\n",V(h3tohaa));
/* for h -> h+h-: */
   fprintf(ff,"h1toHpm %+.9e\n",V(h1toHpm));
   fprintf(ff,"h2toHpm %+.9e\n",V(h2toHpm));
   fprintf(ff,"h3toHpm %+.9e\n",V(h3toHpm));
/* for h -> aZ : */
   fprintf(ff,"h1toaZ1 %+.9e\n",V(h1toaZ1));
   fprintf(ff,"h2toaZ1 %+.9e\n",V(h2toaZ1));
   fprintf(ff,"h3toaZ1 %+.9e\n",V(h3toaZ1));
/* for h -> chi1chi1 : */
   fprintf(ff,"h1toN11 %+.9e\n",V(h1toN11));
   fprintf(ff,"h2toN11 %+.9e\n",V(h2toN11));
   fprintf(ff,"h3toN11 %+.9e\n",V(h3toN11));
/* for h -> chi1chi2 : */
   fprintf(ff,"h1toN12 %+.9e\n",V(h1toN12));
   fprintf(ff,"h2toN12 %+.9e\n",V(h2toN12));
   fprintf(ff,"h3toN12 %+.9e\n",V(h3toN12));
/* for h -> chi1chi3 : */
   fprintf(ff,"h1toN13 %+.9e\n",V(h1toN13));
   fprintf(ff,"h2toN13 %+.9e\n",V(h2toN13));
   fprintf(ff,"h3toN13 %+.9e\n",V(h3toN13));
/* for h -> chi1chi4 : */
   fprintf(ff,"h1toN14 %+.9e\n",V(h1toN14));
   fprintf(ff,"h2toN14 %+.9e\n",V(h2toN14));
   fprintf(ff,"h3toN14 %+.9e\n",V(h3toN14));
/* for h -> chi1chi5 : */
   fprintf(ff,"h1toN15 %+.9e\n",V(h1toN15));
   fprintf(ff,"h2toN15 %+.9e\n",V(h2toN15));
   fprintf(ff,"h3toN15 %+.9e\n",V(h3toN15));
/* for h -> chi1chi6 : */
   fprintf(ff,"h1toN16 %+.9e\n",V(h1toN16));
   fprintf(ff,"h2toN16 %+.9e\n",V(h2toN16));
   fprintf(ff,"h3toN16 %+.9e\n",V(h3toN16));
/* for h -> chi2chi2 : */
   fprintf(ff,"h1toN22 %+.9e\n",V(h1toN22));
   fprintf(ff,"h2toN22 %+.9e\n",V(h2toN22));
   fprintf(ff,"h3toN22 %+.9e\n",V(h3toN22));
/* for h -> chi2chi3 : */
   fprintf(ff,"h1toN23 %+.9e\n",V(h1toN23));
   fprintf(ff,"h2toN23 %+.9e\n",V(h2toN23));
   fprintf(ff,"h3toN23 %+.9e\n",V(h3toN23));
/* for h -> chi2chi4 : */
   fprintf(ff,"h1toN24 %+.9e\n",V(h1toN24));
   fprintf(ff,"h2toN24 %+.9e\n",V(h2toN24));
   fprintf(ff,"h3toN24 %+.9e\n",V(h3toN24));
/* for h -> chi2chi5 : */
   fprintf(ff,"h1toN25 %+.9e\n",V(h1toN25));
   fprintf(ff,"h2toN25 %+.9e\n",V(h2toN25));
   fprintf(ff,"h3toN25 %+.9e\n",V(h3toN25));
/* for h -> chi2chi6 : */
   fprintf(ff,"h1toN26 %+.9e\n",V(h1toN26));
   fprintf(ff,"h2toN26 %+.9e\n",V(h2toN26));
   fprintf(ff,"h3toN26 %+.9e\n",V(h3toN26));
/* for h -> chi3chi3 : */
   fprintf(ff,"h1toN33 %+.9e\n",V(h1toN33));
   fprintf(ff,"h2toN33 %+.9e\n",V(h2toN33));
   fprintf(ff,"h3toN33 %+.9e\n",V(h3toN33));
/* for h -> chi3chi4 : */
   fprintf(ff,"h1toN34 %+.9e\n",V(h1toN34));
   fprintf(ff,"h2toN34 %+.9e\n",V(h2toN34));
   fprintf(ff,"h3toN34 %+.9e\n",V(h3toN34));
/* for h -> chi3chi5 : */
   fprintf(ff,"h1toN35 %+.9e\n",V(h1toN35));
   fprintf(ff,"h2toN35 %+.9e\n",V(h2toN35));
   fprintf(ff,"h3toN35 %+.9e\n",V(h3toN35));
/* for h -> chi3chi6 : */
   fprintf(ff,"h1toN36 %+.9e\n",V(h1toN36));
   fprintf(ff,"h2toN36 %+.9e\n",V(h2toN36));
   fprintf(ff,"h3toN36 %+.9e\n",V(h3toN36));
/* for h -> chi4chi4 : */
   fprintf(ff,"h1toN44 %+.9e\n",V(h1toN44));
   fprintf(ff,"h2toN44 %+.9e\n",V(h2toN44));
   fprintf(ff,"h3toN44 %+.9e\n",V(h3toN44));
/* for h -> chi4chi5 : */
   fprintf(ff,"h1toN45 %+.9e\n",V(h1toN45));
   fprintf(ff,"h2toN45 %+.9e\n",V(h2toN45));
   fprintf(ff,"h3toN45 %+.9e\n",V(h3toN45));
/* for h -> chi4chi6 : */
   fprintf(ff,"h1toN46 %+.9e\n",V(h1toN46));
   fprintf(ff,"h2toN46 %+.9e\n",V(h2toN46));
   fprintf(ff,"h3toN46 %+.9e\n",V(h3toN46));
/* for h -> chi5chi5 : */
   fprintf(ff,"h1toN55 %+.9e\n",V(h1toN55));
   fprintf(ff,"h2toN55 %+.9e\n",V(h2toN55));
   fprintf(ff,"h3toN55 %+.9e\n",V(h3toN55));
/* for h -> chi5chi6 : */
   fprintf(ff,"h1toN56 %+.9e\n",V(h1toN56));
   fprintf(ff,"h2toN56 %+.9e\n",V(h2toN56));
   fprintf(ff,"h3toN56 %+.9e\n",V(h3toN56));
/* for h -> chi6chi6 : */
   fprintf(ff,"h1toN66 %+.9e\n",V(h1toN66));
   fprintf(ff,"h2toN66 %+.9e\n",V(h2toN66));
   fprintf(ff,"h3toN66 %+.9e\n",V(h3toN66));
/* for h -> chi1+chi1- : */
   fprintf(ff,"h1toC11 %+.9e\n",V(h1toC11));
   fprintf(ff,"h2toC11 %+.9e\n",V(h2toC11));
   fprintf(ff,"h3toC11 %+.9e\n",V(h3toC11));
/* for h -> chi1+chi2- : */
   fprintf(ff,"h1toCLH %+.9e\n",V(h1toCLH));
   fprintf(ff,"h2toCLH %+.9e\n",V(h2toCLH));
   fprintf(ff,"h3toCLH %+.9e\n",V(h3toCLH));
   fprintf(ff,"h1toCRH %+.9e\n",V(h1toCRH));
   fprintf(ff,"h2toCRH %+.9e\n",V(h2toCRH));
   fprintf(ff,"h3toCRH %+.9e\n",V(h3toCRH));
/* for h -> chi2+chi2- : */      
   fprintf(ff,"h1toC22 %+.9e\n",V(h1toC22));
   fprintf(ff,"h2toC22 %+.9e\n",V(h2toC22));
   fprintf(ff,"h3toC22 %+.9e\n",V(h3toC22));
/* for a -> chi1chiN : */        
   fprintf(ff,"hatoN11 %+.9e\n",V(hatoN11));
   fprintf(ff,"hatoN12 %+.9e\n",V(hatoN12));
   fprintf(ff,"hatoN13 %+.9e\n",V(hatoN13));
   fprintf(ff,"hatoN14 %+.9e\n",V(hatoN14));
   fprintf(ff,"hatoN15 %+.9e\n",V(hatoN15));
   fprintf(ff,"hatoN16 %+.9e\n",V(hatoN16));
/* for a -> chi2chiN : */        
   fprintf(ff,"hatoN22 %+.9e\n",V(hatoN22));
   fprintf(ff,"hatoN23 %+.9e\n",V(hatoN23));
   fprintf(ff,"hatoN24 %+.9e\n",V(hatoN24));
   fprintf(ff,"hatoN25 %+.9e\n",V(hatoN25));
   fprintf(ff,"hatoN26 %+.9e\n",V(hatoN26));
/* for a -> chi3chiN : */        
   fprintf(ff,"hatoN33 %+.9e\n",V(hatoN33));
   fprintf(ff,"hatoN34 %+.9e\n",V(hatoN34));
   fprintf(ff,"hatoN35 %+.9e\n",V(hatoN35));
   fprintf(ff,"hatoN36 %+.9e\n",V(hatoN36));
/* for a -> chi4chiN : */        
   fprintf(ff,"hatoN44 %+.9e\n",V(hatoN44));
   fprintf(ff,"hatoN45 %+.9e\n",V(hatoN45));
   fprintf(ff,"hatoN46 %+.9e\n",V(hatoN46));
/* for a -> chi5chiN : */        
   fprintf(ff,"hatoN55 %+.9e\n",V(hatoN55));
   fprintf(ff,"hatoN56 %+.9e\n",V(hatoN56));
/* for a -> chi6chiN : */        
   fprintf(ff,"hatoN66 %+.9e\n",V(hatoN66));
/* for a -> chi1+chi1- : */      
   fprintf(ff,"hatoC11 %+.9e\n",V(hatoC11));
/* for a -> chi1+chi2- : */      
   fprintf(ff,"hatoCLH %+.9e\n",V(hatoCLH));
   fprintf(ff,"hatoCRH %+.9e\n",V(hatoCRH));
/* for a -> chi2+chi2- : */      
   fprintf(ff,"hatoC22 %+.9e\n",V(hatoC22));
/* for h+ -> chi1L+chiN : */     
   fprintf(ff,"LtoN1C1 %+.9e\n",V(LtoN1C1));
   fprintf(ff,"LtoN2C1 %+.9e\n",V(LtoN2C1));
   fprintf(ff,"LtoN3C1 %+.9e\n",V(LtoN3C1));
   fprintf(ff,"LtoN4C1 %+.9e\n",V(LtoN4C1));
   fprintf(ff,"LtoN5C1 %+.9e\n",V(LtoN5C1));
   fprintf(ff,"LtoN6C1 %+.9e\n",V(LtoN6C1));
/* for h+ -> chi2L+chiN : */     
   fprintf(ff,"LtoN1C2 %+.9e\n",V(LtoN1C2));
   fprintf(ff,"LtoN2C2 %+.9e\n",V(LtoN2C2));
   fprintf(ff,"LtoN3C2 %+.9e\n",V(LtoN3C2));
   fprintf(ff,"LtoN4C2 %+.9e\n",V(LtoN4C2));
   fprintf(ff,"LtoN5C2 %+.9e\n",V(LtoN5C2));
   fprintf(ff,"LtoN6C2 %+.9e\n",V(LtoN6C2));
/* for h+ -> chi1R+chiN : */     
   fprintf(ff,"RtoN1C1 %+.9e\n",V(RtoN1C1));
   fprintf(ff,"RtoN2C1 %+.9e\n",V(RtoN2C1));
   fprintf(ff,"RtoN3C1 %+.9e\n",V(RtoN3C1));
   fprintf(ff,"RtoN4C1 %+.9e\n",V(RtoN4C1));
   fprintf(ff,"RtoN5C1 %+.9e\n",V(RtoN5C1));
   fprintf(ff,"RtoN6C1 %+.9e\n",V(RtoN6C1));
/* for h+ -> chi2R+chiN : */     
   fprintf(ff,"RtoN1C2 %+.9e\n",V(RtoN1C2));
   fprintf(ff,"RtoN2C2 %+.9e\n",V(RtoN2C2));
   fprintf(ff,"RtoN3C2 %+.9e\n",V(RtoN3C2));
   fprintf(ff,"RtoN4C2 %+.9e\n",V(RtoN4C2));
   fprintf(ff,"RtoN5C2 %+.9e\n",V(RtoN5C2));
   fprintf(ff,"RtoN6C2 %+.9e\n",V(RtoN6C2));
/* for checkmin : */ 
   fprintf(ff,"M2Hd    %+.9e\n",V(M2Hd));
   fprintf(ff,"M2Hu    %+.9e\n",V(M2Hu));
   fprintf(ff,"M2S     %+.9e\n",V(M2S));
   fprintf(ff,"la1     %+.9e\n",V(la1));
   fprintf(ff,"la2     %+.9e\n",V(la2));
   fprintf(ff,"la3     %+.9e\n",V(la3));
   fprintf(ff,"la4     %+.9e\n",V(la4));
   fprintf(ff,"la5     %+.9e\n",V(la5));
   fprintf(ff,"la6     %+.9e\n",V(la6));
   fprintf(ff,"la7     %+.9e\n",V(la7));
   fprintf(ff,"aa5     %+.9e\n",V(aa5));
   fprintf(ff,"la1s    %+.9e\n",V(la1s));
   fprintf(ff,"la2s    %+.9e\n",V(la2s));
/* terms for UMSSMTools/sources/runpar.f : */ 
   fprintf(ff,"M2Q3    %+.9e\n",V(M2Q3));
   fprintf(ff,"M2U3    %+.9e\n",V(M2U3));
   fprintf(ff,"M2D3    %+.9e\n",V(M2D3));
   fprintf(ff,"M2L3    %+.9e\n",V(M2L3));
   fprintf(ff,"M2R3    %+.9e\n",V(M2R3));
/* PDG of the LSP to know how to compute invisible branching of higgses: */ 
   fprintf(ff,"PDG_LSP %7d\n",PDG_LSP);
// if no running of vevs is considered for higgses->sfermions : */
// First change the value of "keys haHpmtosfer" into the LanHEP definition of the UMSSM and then reproduce the modified model files
// Then :
/* for h -> ~uL ~UL : */         
//   fprintf(ff,"h1toSuL %+.9e\n",V(h1toSuL));
//   fprintf(ff,"h2toSuL %+.9e\n",V(h2toSuL));
//   fprintf(ff,"h3toSuL %+.9e\n",V(h3toSuL));
/* for h -> ~uR ~UR : */         
//   fprintf(ff,"h1toSuR %+.9e\n",V(h1toSuR));
//   fprintf(ff,"h2toSuR %+.9e\n",V(h2toSuR));
//   fprintf(ff,"h3toSuR %+.9e\n",V(h3toSuR));
/* for h -> ~dL ~DL : */         
//   fprintf(ff,"h1toSdL %+.9e\n",V(h1toSdL));
//   fprintf(ff,"h2toSdL %+.9e\n",V(h2toSdL));
//   fprintf(ff,"h3toSdL %+.9e\n",V(h3toSdL));
/* for h -> ~dR ~DR : */         
//   fprintf(ff,"h1toSdR %+.9e\n",V(h1toSdR));
//   fprintf(ff,"h2toSdR %+.9e\n",V(h2toSdR));
//   fprintf(ff,"h3toSdR %+.9e\n",V(h3toSdR));
/* for h -> ~t1 ~T1 : */         
//   fprintf(ff,"h1toSt1 %+.9e\n",V(h1toSt1));
//   fprintf(ff,"h2toSt1 %+.9e\n",V(h2toSt1));
//   fprintf(ff,"h3toSt1 %+.9e\n",V(h3toSt1));
/* for h -> ~t2 ~T2 : */         
//   fprintf(ff,"h1toSt2 %+.9e\n",V(h1toSt2));
//   fprintf(ff,"h2toSt2 %+.9e\n",V(h2toSt2));
//   fprintf(ff,"h3toSt2 %+.9e\n",V(h3toSt2));
/* for h -> ~t1 ~T2 : */         
//   fprintf(ff,"h1toStm %+.9e\n",V(h1toStm));
//   fprintf(ff,"h2toStm %+.9e\n",V(h2toStm));
//   fprintf(ff,"h3toStm %+.9e\n",V(h3toStm));
/* for h -> ~b1 ~B1 : */         
//   fprintf(ff,"h1toSb1 %+.9e\n",V(h1toSb1));
//   fprintf(ff,"h2toSb1 %+.9e\n",V(h2toSb1));
//   fprintf(ff,"h3toSb1 %+.9e\n",V(h3toSb1));
/* for h -> ~b2 ~B2 : */         
//   fprintf(ff,"h1toSb2 %+.9e\n",V(h1toSb2));
//   fprintf(ff,"h2toSb2 %+.9e\n",V(h2toSb2));
//   fprintf(ff,"h3toSb2 %+.9e\n",V(h3toSb2));
/* for h -> ~b1 ~B2 : */         
//   fprintf(ff,"h1toSbm %+.9e\n",V(h1toSbm));
//   fprintf(ff,"h2toSbm %+.9e\n",V(h2toSbm));
//   fprintf(ff,"h3toSbm %+.9e\n",V(h3toSbm));
/* for h -> ~eL ~EL : */         
//   fprintf(ff,"h1toSeL %+.9e\n",V(h1toSeL));
//   fprintf(ff,"h2toSeL %+.9e\n",V(h2toSeL));
//   fprintf(ff,"h3toSeL %+.9e\n",V(h3toSeL));
/* for h -> ~eR ~ER : */         
//   fprintf(ff,"h1toSeR %+.9e\n",V(h1toSeR));
//   fprintf(ff,"h2toSeR %+.9e\n",V(h2toSeR));
//   fprintf(ff,"h3toSeR %+.9e\n",V(h3toSeR));
/* for h -> ~nl ~Nl : */         
//   fprintf(ff,"h1toSnL %+.9e\n",V(h1toSnL));
//   fprintf(ff,"h2toSnL %+.9e\n",V(h2toSnL));
//   fprintf(ff,"h3toSnL %+.9e\n",V(h3toSnL));
/* for h -> ~l1 ~L1 : */         
//   fprintf(ff,"h1toSl1 %+.9e\n",V(h1toSl1));
//   fprintf(ff,"h2toSl1 %+.9e\n",V(h2toSl1));
//   fprintf(ff,"h3toSl1 %+.9e\n",V(h3toSl1));
/* for h -> ~l2 ~L2 : */         
//   fprintf(ff,"h1toSl2 %+.9e\n",V(h1toSl2));
//   fprintf(ff,"h2toSl2 %+.9e\n",V(h2toSl2));
//   fprintf(ff,"h3toSl2 %+.9e\n",V(h3toSl2));
/* for h -> ~l1 ~L2 : */         
//   fprintf(ff,"h1toSlm %+.9e\n",V(h1toSlm));
//   fprintf(ff,"h2toSlm %+.9e\n",V(h2toSlm));
//   fprintf(ff,"h3toSlm %+.9e\n",V(h3toSlm));
/* for h -> ~nr ~Nr : */         
//   fprintf(ff,"h1toSnR %+.9e\n",V(h1toSnR));
//   fprintf(ff,"h2toSnR %+.9e\n",V(h2toSnR));
//   fprintf(ff,"h3toSnR %+.9e\n",V(h3toSnR));
/* for a -> ~f1 ~F2 : */         
//   fprintf(ff,"hatoStm %+.9e\n",V(hatoStm));
//   fprintf(ff,"hatoSbm %+.9e\n",V(hatoSbm));
//   fprintf(ff,"hatoSlm %+.9e\n",V(hatoSlm));
/* for h+ -> ~fu1 ~fd2 : */      
//   fprintf(ff,"HcSq311 %+.9e\n",V(HcSq311));
//   fprintf(ff,"HcSq312 %+.9e\n",V(HcSq312));
//   fprintf(ff,"HcSq321 %+.9e\n",V(HcSq321));
//   fprintf(ff,"HcSq322 %+.9e\n",V(HcSq322));
//   fprintf(ff,"HcSnl1  %+.9e\n",V(HcSnl1));
//   fprintf(ff,"HcSnl2  %+.9e\n",V(HcSnl2));


   fclose(ff);

   err= runTools("umssmtoolslib","UMSSM_spectr.dat",0);

   if(err) {FError=1;}
//   nw= slhaWarnings(NULL);
   return err;
}

double bsg_(double *M, double *P)
{ *M= slhaVal("LOWEN",0.,1,12);
  *P= slhaVal("LOWEN",0.,1,11);
  return slhaVal("LOWEN",0.,1,1);
}

double deltamd_(double *M, double *P)
{ *M= slhaVal("LOWEN",0.,1,22);
  *P= slhaVal("LOWEN",0.,1,21);
  return slhaVal("LOWEN",0.,1,2);
}

double deltams_(double *M, double *P)
{ *M= slhaVal("LOWEN",0.,1,32);
  *P= slhaVal("LOWEN",0.,1,31);
  return slhaVal("LOWEN",0.,1,3);
}

double bsmumu_(double *M, double *P)
{ *M= slhaVal("LOWEN",0.,1,42);
  *P= slhaVal("LOWEN",0.,1,41);
  return slhaVal("LOWEN",0.,1,4);
}

double btaunu_(double *M, double *P)
{ *M= slhaVal("LOWEN",0.,1,52);
  *P= slhaVal("LOWEN",0.,1,51);
  return slhaVal("LOWEN",0.,1,5);
}

double gmuon_(double *M, double *P)
{ *M= slhaVal("LOWEN",0.,1,62);
  *P= slhaVal("LOWEN",0.,1,61);
  return slhaVal("LOWEN",0.,1,6);
}

double bxislllow_(double *M, double *P)
{ *M= slhaVal("LOWEN",0.,1,72);
  *P= slhaVal("LOWEN",0.,1,71);
  return slhaVal("LOWEN",0.,1,7);
}

double bxisllhigh_(double *M, double *P)
{ *M= slhaVal("LOWEN",0.,1,82);
  *P= slhaVal("LOWEN",0.,1,81);
  return slhaVal("LOWEN",0.,1,8);
}

double bdg_(double *M, double *P)
{ *M= slhaVal("LOWEN",0.,1,92);
  *P= slhaVal("LOWEN",0.,1,91);
  return slhaVal("LOWEN",0.,1,9);
}

double bdmumu_(double *M, double *P)
{ *M= slhaVal("LOWEN",0.,1,102);
  *P= slhaVal("LOWEN",0.,1,101);
  return slhaVal("LOWEN",0.,1,10);
}

double bxisnunu_(double *M, double *P)
{ *M= slhaVal("LOWEN",0.,1,112);
  *P= slhaVal("LOWEN",0.,1,111);
  return slhaVal("LOWEN",0.,1,110);
}

double bpkpnunu_(double *M, double *P)
{ *M= slhaVal("LOWEN",0.,1,122);
  *P= slhaVal("LOWEN",0.,1,121);
  return slhaVal("LOWEN",0.,1,120);
}

double bksnunu_(double *M, double *P)
{ *M= slhaVal("LOWEN",0.,1,132);
  *P= slhaVal("LOWEN",0.,1,131);
  return slhaVal("LOWEN",0.,1,13);
}

double rdtaul_(double *M, double *P)
{ *M= slhaVal("LOWEN",0.,1,142);
  *P= slhaVal("LOWEN",0.,1,141);
  return slhaVal("LOWEN",0.,1,14);
}

double rdstaul_(double *M, double *P)
{ *M= slhaVal("LOWEN",0.,1,152);
  *P= slhaVal("LOWEN",0.,1,151);
  return slhaVal("LOWEN",0.,1,15);
}

double kppipnunu_(double *M, double *P)
{ *M= slhaVal("LOWEN",0.,1,162);
  *P= slhaVal("LOWEN",0.,1,161);
  return slhaVal("LOWEN",0.,1,16);
}

double klpi0nunu_(double *M, double *P)
{ *M= slhaVal("LOWEN",0.,1,172);
  *P= slhaVal("LOWEN",0.,1,171);
  return slhaVal("LOWEN",0.,1,17);
}

double deltamk_(double *M, double *P)
{ *M= slhaVal("LOWEN",0.,1,182);
  *P= slhaVal("LOWEN",0.,1,181);
  return slhaVal("LOWEN",0.,1,18);
}

double epsk_(double *M, double *P)
{ *M= slhaVal("LOWEN",0.,1,192);
  *P= slhaVal("LOWEN",0.,1,191);
  return slhaVal("LOWEN",0.,1,19);
}
