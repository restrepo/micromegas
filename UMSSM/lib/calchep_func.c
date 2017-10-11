#include "pmodel.h"
#include "../../sources/micromegas.h"
#include "../../sources/micromegas_aux.h"
/*#include"../../CalcHEP_src/c_source/model_aux/include/SLHAreader.h"*/

#include "lpath.h"
#include <math.h>

void o1Contents(FILE * f)
{
  fprintf(f,"\n~o1 = %.3f*bino %+.3f*wino %+.3f*higgsino1 %+.3f*higgsino2 %+.3f*singlino %+.3f*binoprime\n",
      findValW("Zn11"), findValW("Zn12"), findValW("Zn13"), findValW("Zn14"), findValW("Zn15"), findValW("Zn16"));
}

//Constraint on the invisible Z-width = 0.5 MeV (see remaining uncertainties on Gamma_Z in A. Freitas, 1401.2447)
int Zinv(void)
{
txtList L;
double width;
char pname[3]="", LightestOdd[4]="", LightestaOdd[4]="", Inv[9]="";
strcat(pname,pdg2name(23));
strcat(LightestOdd,CDM1);
strcat(LightestaOdd,pdg2name(-pNum(CDM1)));
if(!strcmp(LightestaOdd,""))strcat(LightestaOdd,LightestOdd);
strcat(Inv,LightestOdd);strcat(Inv,",");strcat(Inv,LightestaOdd);
width=pWidth(pname,&L);
//printf("%s-width into %s : %f MeV\n",pdg2name(23),Inv,findBr(L,Inv)*width*1000.);
if(findBr(L,Inv)*width*1000. > 0.5){printf("WARNING: invisible %s-width above 0.5 MeV \n",pdg2name(23));return 1;}
else return 0;
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
    printf("*********** UMSSMTools_1.0.1 ************\n");
    printf("**** Last update : 18 september 2015 ****\n");
   err=System(command);printf("*****************************************\n"); 
   free(command);
   if(err>=0) err=slhaRead(fout,0);

   return err; 
}



/* =====  end of header part ========= */


#define V(N) findValW(#N)
  int UMSSMTools(void)
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
   fprintf(f," 47  %.8E      # Mn1\n",V(Mn2));
   fprintf(f," 48  %.8E      # Mn2\n",V(Mn2));
   fprintf(f," 49  %.8E      # Mnlr\n",V(Mnlr));
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
   fprintf(ff,"MG1    %+.9e\n",V(MG1));
   fprintf(ff,"MG2    %+.9e\n",V(MG2));
   fprintf(ff,"MG3    %+.9e\n",V(MG3));
   fprintf(ff,"M1p    %+.9e\n",V(M1p));
   fprintf(ff,"MK     %+.9e\n",V(MK));
   fprintf(ff,"Ml2    %+.9e\n",V(Ml2));
   fprintf(ff,"Ml3    %+.9e\n",V(Ml3));
   fprintf(ff,"Mr2    %+.9e\n",V(Mr2));
   fprintf(ff,"Mr3    %+.9e\n",V(Mr3));
   fprintf(ff,"Mq2    %+.9e\n",V(Mq2));
   fprintf(ff,"Mq3    %+.9e\n",V(Mq3));
   fprintf(ff,"Mu2    %+.9e\n",V(Mu2));
   fprintf(ff,"Mu3    %+.9e\n",V(Mu3));
   fprintf(ff,"Md2    %+.9e\n",V(Md2));
   fprintf(ff,"Md3    %+.9e\n",V(Md3));
   fprintf(ff,"Mn2    %+.9e\n",V(Mn2));
   fprintf(ff,"Mn32   %+.9e\n",V(Mn32));
   fprintf(ff,"At     %+.9e\n",V(At));
   fprintf(ff,"Ab     %+.9e\n",V(Ab));
   fprintf(ff,"Al     %+.9e\n",V(Al));
   fprintf(ff,"Am     %+.9e\n",V(Am));
   fprintf(ff,"SW     %+.9e\n",V(SW));
   fprintf(ff,"g2     %+.9e\n",V(g2));
   fprintf(ff,"gp     %+.9e\n",V(gp));
   fprintf(ff,"g1p    %+.9e\n",V(g1p));
   fprintf(ff,"v      %+.9e\n",V(v));
   fprintf(ff,"vs     %+.9e\n",V(vs));
   fprintf(ff,"NCp    %+.9e\n",V(NCp));
   fprintf(ff,"Q1     %+.9e\n",V(Q1));
   fprintf(ff,"Q2     %+.9e\n",V(Q2));
   fprintf(ff,"QQ     %+.9e\n",V(QQ));
   fprintf(ff,"QL     %+.9e\n",V(QL));
   fprintf(ff,"MW     %+.9e\n",V(MW));
   fprintf(ff,"MZ1    %+.9e\n",V(MZ1));
   fprintf(ff,"MZ2    %+.9e\n",V(MZ2));
   fprintf(ff,"aZZ    %+.9e\n",V(aZZ));
   fprintf(ff,"saZZ   %+.9e\n",V(saZZ));
   fprintf(ff,"caZZ   %+.9e\n",V(caZZ));
   fprintf(ff,"mu     %+.9e\n",V(mu));
   fprintf(ff,"Alda   %+.9e\n",V(Alda));
   fprintf(ff,"tE6    %+.9e\n",V(tE6));
   fprintf(ff,"Mtp    %+.9e\n",V(Mtp));
   fprintf(ff,"MbMb   %+.9e\n",V(MbMb));
   fprintf(ff,"Ms2GeV %+.9e\n",V(Ms2GeV));
   fprintf(ff,"tb     %+.9e\n",V(tb));
   fprintf(ff,"lda    %+.9e\n",V(lda));
   fprintf(ff,"ZA11   %+.9e\n",V(ZA11));
   fprintf(ff,"ZA12   %+.9e\n",V(ZA12));
   fprintf(ff,"ZA13   %+.9e\n",V(ZA13));
   fprintf(ff,"ZA21   %+.9e\n",V(ZA21));
   fprintf(ff,"ZA22   %+.9e\n",V(ZA22));
   fprintf(ff,"ZA23   %+.9e\n",V(ZA23));
   fprintf(ff,"ZA31   %+.9e\n",V(ZA31));
   fprintf(ff,"ZA32   %+.9e\n",V(ZA32));
   fprintf(ff,"ZA33   %+.9e\n",V(ZA33));
   fprintf(ff,"Mha    %+.9e\n",V(Mha));
   fprintf(ff,"Zh11   %+.9e\n",V(Zh11));
   fprintf(ff,"Zh12   %+.9e\n",V(Zh12));
   fprintf(ff,"Zh13   %+.9e\n",V(Zh13));
   fprintf(ff,"Zh21   %+.9e\n",V(Zh21));
   fprintf(ff,"Zh22   %+.9e\n",V(Zh22));
   fprintf(ff,"Zh23   %+.9e\n",V(Zh23));
   fprintf(ff,"Zh31   %+.9e\n",V(Zh31));
   fprintf(ff,"Zh32   %+.9e\n",V(Zh32));
   fprintf(ff,"Zh33   %+.9e\n",V(Zh33));
   fprintf(ff,"Mh1    %+.9e\n",V(Mh1));
   fprintf(ff,"Mh2    %+.9e\n",V(Mh2));
   fprintf(ff,"Mh3    %+.9e\n",V(Mh3));
   fprintf(ff,"MHc    %+.9e\n",V(MHc));
   fprintf(ff,"Zu11   %+.9e\n",V(Zu11));
   fprintf(ff,"Zu12   %+.9e\n",V(Zu12));
   fprintf(ff,"Zu21   %+.9e\n",V(Zu21));
   fprintf(ff,"Zu22   %+.9e\n",V(Zu22));
   fprintf(ff,"Zv11   %+.9e\n",V(Zv11));
   fprintf(ff,"Zv12   %+.9e\n",V(Zv12));
   fprintf(ff,"Zv21   %+.9e\n",V(Zv21));
   fprintf(ff,"Zv22   %+.9e\n",V(Zv22));
   fprintf(ff,"MC1    %+.9e\n",V(MC1));
   fprintf(ff,"MC2    %+.9e\n",V(MC2));
   fprintf(ff,"MSuL   %+.9e\n",V(MSuL));
   fprintf(ff,"MSuR   %+.9e\n",V(MSuR));
   fprintf(ff,"MSdL   %+.9e\n",V(MSdL));
   fprintf(ff,"MSdR   %+.9e\n",V(MSdR));
   fprintf(ff,"Zt11   %+.9e\n",V(Zt11));
   fprintf(ff,"Zt12   %+.9e\n",V(Zt12));
   fprintf(ff,"MSt1   %+.9e\n",V(MSt1));
   fprintf(ff,"MSt2   %+.9e\n",V(MSt2));
   fprintf(ff,"Zb11   %+.9e\n",V(Zb11));
   fprintf(ff,"Zb12   %+.9e\n",V(Zb12));
   fprintf(ff,"MSb1   %+.9e\n",V(MSb1));
   fprintf(ff,"MSb2   %+.9e\n",V(MSb2));
   fprintf(ff,"MSeL   %+.9e\n",V(MSeL));
   fprintf(ff,"MSeR   %+.9e\n",V(MSeR));
   fprintf(ff,"MSmL   %+.9e\n",V(MSmL));
   fprintf(ff,"MSmR   %+.9e\n",V(MSmR));
   fprintf(ff,"Zl11   %+.9e\n",V(Zl11));
   fprintf(ff,"Zl12   %+.9e\n",V(Zl12));
   fprintf(ff,"MSl1   %+.9e\n",V(MSl1));
   fprintf(ff,"MSl2   %+.9e\n",V(MSl2));
   fprintf(ff,"MSne   %+.9e\n",V(MSne));
   fprintf(ff,"MSnm   %+.9e\n",V(MSnm));
   fprintf(ff,"MSnl   %+.9e\n",V(MSnl));
   fprintf(ff,"Mner   %+.9e\n",V(Mner));
   fprintf(ff,"Mnmr   %+.9e\n",V(Mnmr));
   fprintf(ff,"Mnlr   %+.9e\n",V(Mnlr));
   fprintf(ff,"Zn11   %+.9e\n",V(Zn11));
   fprintf(ff,"Zn12   %+.9e\n",V(Zn12));
   fprintf(ff,"Zn13   %+.9e\n",V(Zn13));
   fprintf(ff,"Zn14   %+.9e\n",V(Zn14));
   fprintf(ff,"Zn15   %+.9e\n",V(Zn15));
   fprintf(ff,"Zn16   %+.9e\n",V(Zn16));
   fprintf(ff,"Zn21   %+.9e\n",V(Zn21));
   fprintf(ff,"Zn22   %+.9e\n",V(Zn22));
   fprintf(ff,"Zn23   %+.9e\n",V(Zn23));
   fprintf(ff,"Zn24   %+.9e\n",V(Zn24));
   fprintf(ff,"Zn25   %+.9e\n",V(Zn25));
   fprintf(ff,"Zn26   %+.9e\n",V(Zn26));
   fprintf(ff,"Zn31   %+.9e\n",V(Zn31));
   fprintf(ff,"Zn32   %+.9e\n",V(Zn32));
   fprintf(ff,"Zn33   %+.9e\n",V(Zn33));
   fprintf(ff,"Zn34   %+.9e\n",V(Zn34));
   fprintf(ff,"Zn35   %+.9e\n",V(Zn35));
   fprintf(ff,"Zn36   %+.9e\n",V(Zn36));
   fprintf(ff,"Zn41   %+.9e\n",V(Zn41));
   fprintf(ff,"Zn42   %+.9e\n",V(Zn42));
   fprintf(ff,"Zn43   %+.9e\n",V(Zn43));
   fprintf(ff,"Zn44   %+.9e\n",V(Zn44));
   fprintf(ff,"Zn45   %+.9e\n",V(Zn45));
   fprintf(ff,"Zn46   %+.9e\n",V(Zn46));
   fprintf(ff,"Zn51   %+.9e\n",V(Zn51));
   fprintf(ff,"Zn52   %+.9e\n",V(Zn52));
   fprintf(ff,"Zn53   %+.9e\n",V(Zn53));
   fprintf(ff,"Zn54   %+.9e\n",V(Zn54));
   fprintf(ff,"Zn55   %+.9e\n",V(Zn55));
   fprintf(ff,"Zn56   %+.9e\n",V(Zn56));
   fprintf(ff,"Zn61   %+.9e\n",V(Zn61));
   fprintf(ff,"Zn62   %+.9e\n",V(Zn62));
   fprintf(ff,"Zn63   %+.9e\n",V(Zn63));
   fprintf(ff,"Zn64   %+.9e\n",V(Zn64));
   fprintf(ff,"Zn65   %+.9e\n",V(Zn65));
   fprintf(ff,"Zn66   %+.9e\n",V(Zn66));
   fprintf(ff,"MNE1   %+.9e\n",V(MNE1));
   fprintf(ff,"MNE2   %+.9e\n",V(MNE2));
   fprintf(ff,"MNE3   %+.9e\n",V(MNE3));
   fprintf(ff,"MNE4   %+.9e\n",V(MNE4));
   fprintf(ff,"MNE5   %+.9e\n",V(MNE5));
   fprintf(ff,"MNE6   %+.9e\n",V(MNE6));
   fprintf(ff,"MSG    %+.9e\n",V(MSG));
/* for h -> hh : */
   fprintf(ff,"B00874 %+.9e\n",V(B00874));
   fprintf(ff,"B00888 %+.9e\n",V(B00888));
   fprintf(ff,"B00922 %+.9e\n",V(B00922));
   fprintf(ff,"B01023 %+.9e\n",V(B01023));
/* for h -> aa : */
   fprintf(ff,"B00945 %+.9e\n",V(B00945));
   fprintf(ff,"B01046 %+.9e\n",V(B01046));
   fprintf(ff,"B01119 %+.9e\n",V(B01119));
/* for h -> h+h-: */
   fprintf(ff,"B00086 %+.9e\n",V(B00086));
   fprintf(ff,"B00093 %+.9e\n",V(B00093));
   fprintf(ff,"B00100 %+.9e\n",V(B00100));
/* for h -> aZ : */
   fprintf(ff,"B00301 %+.9e\n",V(B00301));
   fprintf(ff,"B00303 %+.9e\n",V(B00303));
   fprintf(ff,"B00305 %+.9e\n",V(B00305));
/* for h -> chi1chi1 : */
   fprintf(ff,"B01343 %+.9e\n",V(B01343));
   fprintf(ff,"B01345 %+.9e\n",V(B01345));
   fprintf(ff,"B01347 %+.9e\n",V(B01347));
/* for h -> chi1chi2 : */
   fprintf(ff,"B01371 %+.9e\n",V(B01371));
   fprintf(ff,"B01375 %+.9e\n",V(B01375));
   fprintf(ff,"B01379 %+.9e\n",V(B01379));
/* for h -> chi1chi3 : */
   fprintf(ff,"B01405 %+.9e\n",V(B01405));
   fprintf(ff,"B01409 %+.9e\n",V(B01409));
   fprintf(ff,"B01413 %+.9e\n",V(B01413));
/* for h -> chi1chi4 : */
   fprintf(ff,"B01439 %+.9e\n",V(B01439));
   fprintf(ff,"B01443 %+.9e\n",V(B01443));
   fprintf(ff,"B01447 %+.9e\n",V(B01447));
/* for h -> chi1chi5 : */
   fprintf(ff,"B01473 %+.9e\n",V(B01473));
   fprintf(ff,"B01477 %+.9e\n",V(B01477));
   fprintf(ff,"B01481 %+.9e\n",V(B01481));
/* for h -> chi1chi6 : */
   fprintf(ff,"B01507 %+.9e\n",V(B01507));
   fprintf(ff,"B01511 %+.9e\n",V(B01511));
   fprintf(ff,"B01515 %+.9e\n",V(B01515));
/* for h -> chi2chi2 : */
   fprintf(ff,"B01531 %+.9e\n",V(B01531));
   fprintf(ff,"B01533 %+.9e\n",V(B01533));
   fprintf(ff,"B01535 %+.9e\n",V(B01535));
/* for h -> chi2chi3 : */
   fprintf(ff,"B01559 %+.9e\n",V(B01559));
   fprintf(ff,"B01563 %+.9e\n",V(B01563));
   fprintf(ff,"B01567 %+.9e\n",V(B01567));
/* for h -> chi2chi4 : */
   fprintf(ff,"B01593 %+.9e\n",V(B01593));
   fprintf(ff,"B01597 %+.9e\n",V(B01597));
   fprintf(ff,"B01601 %+.9e\n",V(B01601));
/* for h -> chi2chi5 : */
   fprintf(ff,"B01627 %+.9e\n",V(B01627));
   fprintf(ff,"B01631 %+.9e\n",V(B01631));
   fprintf(ff,"B01635 %+.9e\n",V(B01635));
/* for h -> chi2chi6 : */
   fprintf(ff,"B01661 %+.9e\n",V(B01661));
   fprintf(ff,"B01665 %+.9e\n",V(B01665));
   fprintf(ff,"B01669 %+.9e\n",V(B01669));
/* for h -> chi3chi3 : */
   fprintf(ff,"B01685 %+.9e\n",V(B01685));
   fprintf(ff,"B01687 %+.9e\n",V(B01687));
   fprintf(ff,"B01689 %+.9e\n",V(B01689));
/* for h -> chi3chi4 : */
   fprintf(ff,"B01713 %+.9e\n",V(B01713));
   fprintf(ff,"B01717 %+.9e\n",V(B01717));
   fprintf(ff,"B01721 %+.9e\n",V(B01721));
/* for h -> chi3chi5 : */
   fprintf(ff,"B01747 %+.9e\n",V(B01747));
   fprintf(ff,"B01751 %+.9e\n",V(B01751));
   fprintf(ff,"B01755 %+.9e\n",V(B01755));
/* for h -> chi3chi6 : */
   fprintf(ff,"B01781 %+.9e\n",V(B01781));
   fprintf(ff,"B01785 %+.9e\n",V(B01785));
   fprintf(ff,"B01789 %+.9e\n",V(B01789));
/* for h -> chi4chi4 : */
   fprintf(ff,"B01805 %+.9e\n",V(B01805));
   fprintf(ff,"B01807 %+.9e\n",V(B01807));
   fprintf(ff,"B01809 %+.9e\n",V(B01809));
/* for h -> chi4chi5 : */
   fprintf(ff,"B01833 %+.9e\n",V(B01833));
   fprintf(ff,"B01837 %+.9e\n",V(B01837));
   fprintf(ff,"B01841 %+.9e\n",V(B01841));
/* for h -> chi4chi6 : */
   fprintf(ff,"B01867 %+.9e\n",V(B01867));
   fprintf(ff,"B01871 %+.9e\n",V(B01871));
   fprintf(ff,"B01875 %+.9e\n",V(B01875));
/* for h -> chi5chi5 : */
   fprintf(ff,"B01891 %+.9e\n",V(B01891));
   fprintf(ff,"B01893 %+.9e\n",V(B01893));
   fprintf(ff,"B01895 %+.9e\n",V(B01895));
/* for h -> chi5chi6 : */
   fprintf(ff,"B01919 %+.9e\n",V(B01919));
   fprintf(ff,"B01923 %+.9e\n",V(B01923));
   fprintf(ff,"B01927 %+.9e\n",V(B01927));
/* for h -> chi6chi6 : */
   fprintf(ff,"B01943 %+.9e\n",V(B01943));
   fprintf(ff,"B01945 %+.9e\n",V(B01945));
   fprintf(ff,"B01947 %+.9e\n",V(B01947));
/* for h -> chi1+chi1- : */
   fprintf(ff,"B01250 %+.9e\n",V(B01250));
   fprintf(ff,"B01251 %+.9e\n",V(B01251));
   fprintf(ff,"B01252 %+.9e\n",V(B01252));
/* for h -> chi1+chi2- : */
   fprintf(ff,"B01266 %+.9e\n",V(B01266));
   fprintf(ff,"B01268 %+.9e\n",V(B01268));
   fprintf(ff,"B01270 %+.9e\n",V(B01270));
   fprintf(ff,"B01267 %+.9e\n",V(B01267));
   fprintf(ff,"B01269 %+.9e\n",V(B01269));
   fprintf(ff,"B01271 %+.9e\n",V(B01271));
/* for h -> chi2+chi2- : */      
   fprintf(ff,"B01328 %+.9e\n",V(B01328));
   fprintf(ff,"B01329 %+.9e\n",V(B01329));
   fprintf(ff,"B01330 %+.9e\n",V(B01330));
/* for h -> ~uL ~UL : */         
   fprintf(ff,"B00963 %+.9e\n",V(B00963));
   fprintf(ff,"B01064 %+.9e\n",V(B01064));
   fprintf(ff,"B01137 %+.9e\n",V(B01137));
/* for h -> ~uR ~UR : */         
   fprintf(ff,"B00964 %+.9e\n",V(B00964));
   fprintf(ff,"B01065 %+.9e\n",V(B01065));
   fprintf(ff,"B01138 %+.9e\n",V(B01138));
/* for h -> ~dL ~DL : */         
   fprintf(ff,"B00966 %+.9e\n",V(B00966));
   fprintf(ff,"B01067 %+.9e\n",V(B01067));
   fprintf(ff,"B01140 %+.9e\n",V(B01140));
/* for h -> ~dR ~DR : */         
   fprintf(ff,"B00967 %+.9e\n",V(B00967));
   fprintf(ff,"B01068 %+.9e\n",V(B01068));
   fprintf(ff,"B01141 %+.9e\n",V(B01141));
/* for h -> ~t1 ~T1 : */         
   fprintf(ff,"B00989 %+.9e\n",V(B00989));
   fprintf(ff,"B01090 %+.9e\n",V(B01090));
   fprintf(ff,"B01163 %+.9e\n",V(B01163));
/* for h -> ~t2 ~T2 : */         
   fprintf(ff,"B01001 %+.9e\n",V(B01001));
   fprintf(ff,"B01102 %+.9e\n",V(B01102));
   fprintf(ff,"B01175 %+.9e\n",V(B01175));
/* for h -> ~t1 ~T2 : */         
   fprintf(ff,"B00997 %+.9e\n",V(B00997));
   fprintf(ff,"B01098 %+.9e\n",V(B01098));
   fprintf(ff,"B01171 %+.9e\n",V(B01171));
/* for h -> ~b1 ~B1 : */         
   fprintf(ff,"B00949 %+.9e\n",V(B00949));
   fprintf(ff,"B01050 %+.9e\n",V(B01050));
   fprintf(ff,"B01123 %+.9e\n",V(B01123));
/* for h -> ~b2 ~B2 : */         
   fprintf(ff,"B00961 %+.9e\n",V(B00961));
   fprintf(ff,"B01062 %+.9e\n",V(B01062));
   fprintf(ff,"B01135 %+.9e\n",V(B01135));
/* for h -> ~b1 ~B2 : */         
   fprintf(ff,"B00957 %+.9e\n",V(B00957));
   fprintf(ff,"B01058 %+.9e\n",V(B01058));
   fprintf(ff,"B01131 %+.9e\n",V(B01131));
/* for h -> ~eL ~EL : */         
   fprintf(ff,"B00968 %+.9e\n",V(B00968));
   fprintf(ff,"B01069 %+.9e\n",V(B01069));
   fprintf(ff,"B01142 %+.9e\n",V(B01142));
/* for h -> ~eR ~ER : */         
   fprintf(ff,"B00969 %+.9e\n",V(B00969));
   fprintf(ff,"B01070 %+.9e\n",V(B01070));
   fprintf(ff,"B01143 %+.9e\n",V(B01143));
/* for h -> ~nl ~Nl : */         
   fprintf(ff,"B00985 %+.9e\n",V(B00985));
   fprintf(ff,"B01086 %+.9e\n",V(B01086));
   fprintf(ff,"B01159 %+.9e\n",V(B01159));
/* for h -> ~l1 ~L1 : */         
   fprintf(ff,"B00972 %+.9e\n",V(B00972));
   fprintf(ff,"B01073 %+.9e\n",V(B01073));
   fprintf(ff,"B01146 %+.9e\n",V(B01146));
/* for h -> ~l2 ~L2 : */         
   fprintf(ff,"B00983 %+.9e\n",V(B00983));
   fprintf(ff,"B01084 %+.9e\n",V(B01084));
   fprintf(ff,"B01157 %+.9e\n",V(B01157));
/* for h -> ~l1 ~L2 : */         
   fprintf(ff,"B00980 %+.9e\n",V(B00980));
   fprintf(ff,"B01081 %+.9e\n",V(B01081));
   fprintf(ff,"B01154 %+.9e\n",V(B01154));
/* for h -> ~nr ~Nr : */         
   fprintf(ff,"B00984 %+.9e\n",V(B00984));
   fprintf(ff,"B01085 %+.9e\n",V(B01085));
   fprintf(ff,"B01158 %+.9e\n",V(B01158));
/* for a -> chi1chiN : */        
   fprintf(ff,"B01349 %+.9e\n",V(B01349));
   fprintf(ff,"B01383 %+.9e\n",V(B01383));
   fprintf(ff,"B01417 %+.9e\n",V(B01417));
   fprintf(ff,"B01451 %+.9e\n",V(B01451));
   fprintf(ff,"B01485 %+.9e\n",V(B01485));
   fprintf(ff,"B01519 %+.9e\n",V(B01519));
/* for a -> chi2chiN : */        
   fprintf(ff,"B01537 %+.9e\n",V(B01537));
   fprintf(ff,"B01571 %+.9e\n",V(B01571));
   fprintf(ff,"B01605 %+.9e\n",V(B01605));
   fprintf(ff,"B01639 %+.9e\n",V(B01639));
   fprintf(ff,"B01673 %+.9e\n",V(B01673));
/* for a -> chi3chiN : */        
   fprintf(ff,"B01691 %+.9e\n",V(B01691));
   fprintf(ff,"B01725 %+.9e\n",V(B01725));
   fprintf(ff,"B01759 %+.9e\n",V(B01759));
   fprintf(ff,"B01793 %+.9e\n",V(B01793));
/* for a -> chi4chiN : */        
   fprintf(ff,"B01811 %+.9e\n",V(B01811));
   fprintf(ff,"B01845 %+.9e\n",V(B01845));
   fprintf(ff,"B01879 %+.9e\n",V(B01879));
/* for a -> chi5chiN : */        
   fprintf(ff,"B01897 %+.9e\n",V(B01897));
   fprintf(ff,"B01931 %+.9e\n",V(B01931));
/* for a -> chi6chiN : */        
   fprintf(ff,"B01949 %+.9e\n",V(B01949));
/* for a -> chi1+chi1- : */      
   fprintf(ff,"B01253 %+.9e\n",V(B01253));
/* for a -> chi1+chi2- : */      
   fprintf(ff,"B01272 %+.9e\n",V(B01272));
   fprintf(ff,"B01273 %+.9e\n",V(B01273));
/* for a -> chi2+chi2- : */      
   fprintf(ff,"B01331 %+.9e\n",V(B01331));
/* for a -> ~f1 ~F2 : */         
   fprintf(ff,"B01187 %+.9e\n",V(B01187));
   fprintf(ff,"B01179 %+.9e\n",V(B01179));
   fprintf(ff,"B01183 %+.9e\n",V(B01183));
/* for h+ -> chi1L+chiN : */     
   fprintf(ff,"B01206 %+.9e\n",V(B01206));
   fprintf(ff,"B01212 %+.9e\n",V(B01212));
   fprintf(ff,"B01218 %+.9e\n",V(B01218));
   fprintf(ff,"B01224 %+.9e\n",V(B01224));
   fprintf(ff,"B01230 %+.9e\n",V(B01230));
   fprintf(ff,"B01236 %+.9e\n",V(B01236));
/* for h+ -> chi2L+chiN : */     
   fprintf(ff,"B01274 %+.9e\n",V(B01274));
   fprintf(ff,"B01280 %+.9e\n",V(B01280));
   fprintf(ff,"B01286 %+.9e\n",V(B01286));
   fprintf(ff,"B01292 %+.9e\n",V(B01292));
   fprintf(ff,"B01298 %+.9e\n",V(B01298));
   fprintf(ff,"B01304 %+.9e\n",V(B01304));
/* for h+ -> chi1R+chiN : */     
   fprintf(ff,"B01207 %+.9e\n",V(B01207));
   fprintf(ff,"B01213 %+.9e\n",V(B01213));
   fprintf(ff,"B01219 %+.9e\n",V(B01219));
   fprintf(ff,"B01225 %+.9e\n",V(B01225));
   fprintf(ff,"B01231 %+.9e\n",V(B01231));
   fprintf(ff,"B01237 %+.9e\n",V(B01237));
/* for h+ -> chi2R+chiN : */     
   fprintf(ff,"B01275 %+.9e\n",V(B01275));
   fprintf(ff,"B01281 %+.9e\n",V(B01281));
   fprintf(ff,"B01287 %+.9e\n",V(B01287));
   fprintf(ff,"B01293 %+.9e\n",V(B01293));
   fprintf(ff,"B01299 %+.9e\n",V(B01299));
   fprintf(ff,"B01305 %+.9e\n",V(B01305));
/* for h+ -> ~fu1 ~fd2 : */      
   fprintf(ff,"B00146 %+.9e\n",V(B00146));
   fprintf(ff,"B00148 %+.9e\n",V(B00148));
   fprintf(ff,"B00150 %+.9e\n",V(B00150));
   fprintf(ff,"B00152 %+.9e\n",V(B00152));
   fprintf(ff,"B00143 %+.9e\n",V(B00143));
   fprintf(ff,"B00144 %+.9e\n",V(B00144));
/* for checkmin : */ 
   fprintf(ff,"M2Hd   %+.9e\n",V(M2Hd));
   fprintf(ff,"M2Hu   %+.9e\n",V(M2Hu));
   fprintf(ff,"M2S    %+.9e\n",V(M2S));
   fprintf(ff,"la1    %+.9e\n",V(la1));
   fprintf(ff,"la2    %+.9e\n",V(la2));
   fprintf(ff,"la3    %+.9e\n",V(la3));
   fprintf(ff,"la4    %+.9e\n",V(la4));
   fprintf(ff,"la5    %+.9e\n",V(la5));
   fprintf(ff,"la6    %+.9e\n",V(la6));
   fprintf(ff,"la7    %+.9e\n",V(la7));
   fprintf(ff,"aa5    %+.9e\n",V(aa5));
   fprintf(ff,"la1s   %+.9e\n",V(la1s));
   fprintf(ff,"la2s   %+.9e\n",V(la2s));
/* last terms for runpar : */ 
   fprintf(ff,"M2Q3   %+.9e\n",V(M2Q3));
   fprintf(ff,"M2U3   %+.9e\n",V(M2U3));
   fprintf(ff,"M2D3   %+.9e\n",V(M2D3));
   fprintf(ff,"M2L3   %+.9e\n",V(M2L3));
   fprintf(ff,"M2R3   %+.9e\n",V(M2R3));


   fclose(ff);

   err= runTools("umhdecay","UMSSM_spectr.dat",0);

   if(err) {FError=1;}
//   nw= slhaWarnings(NULL);
   return err;
}




