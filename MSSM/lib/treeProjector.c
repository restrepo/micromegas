#include"../../include/micromegas.h"

#include"pmodel.h"



int treeMSSM(void )
{
  double MZ= findValW("MZ");
  double MW= findValW("MW");  
  double EE=findValW("EE");
  double McMc=1.4;
  double Sqrt2=sqrt(2);
  
  double  alfSMZ=slhaVal("SMINPUTS",0., 1,3);     
  double  MbMb  =slhaVal("SMINPUTS",0., 1,5); 
  double  Mtp   =slhaVal("SMINPUTS",0., 1,6); 
  double  Ml    =slhaVal("SMINPUTS",0., 1,7); 

  double QSUSY=sqrt(slhaVal("MASS",91,1,1000006)*slhaVal("MASS",91,1,2000006));
   
  double       MG1=slhaVal("MSOFT",QSUSY,1, 1);
  double       MG2=slhaVal("MSOFT",QSUSY,1, 2);
  double       MG3=slhaVal("MSOFT",QSUSY,1, 3);
  double       Ml1=slhaVal("MSOFT",QSUSY,1,31 );
  double       Ml2=slhaVal("MSOFT",QSUSY,1,32 );
  double       Ml3=slhaVal("MSOFT",QSUSY,1,33 );
  double       Mr1=slhaVal("MSOFT",QSUSY,1,34 );
  double       Mr2=slhaVal("MSOFT",QSUSY,1,35 );
  double       Mr3=slhaVal("MSOFT",QSUSY,1,36 );
  double       Mq1=slhaVal("MSOFT",QSUSY,1,41 );
  double       Mq2=slhaVal("MSOFT",QSUSY,1,42 );
  double       Mq3=slhaVal("MSOFT",QSUSY,1,43 );
  double       Mu1=slhaVal("MSOFT",QSUSY,1,44 );
  double       Mu2=slhaVal("MSOFT",QSUSY,1,45 );
  double       Mu3=slhaVal("MSOFT",QSUSY,1,46 );
  double       Md1=slhaVal("MSOFT",QSUSY,1,47 );
  double       Md2=slhaVal("MSOFT",QSUSY,1,48 );
  double       Md3=slhaVal("MSOFT",QSUSY,1,49 );

  double       At= slhaVal("Au",QSUSY,2,3,3 );
  double       Ab= slhaVal("Ad",QSUSY,2,3,3 );
  double       Al= slhaVal("Ae",QSUSY,2,3,3 );

  double       mu= slhaVal("HMIX",QSUSY,1,1 );
  double       MH3=slhaVal("MASS",QSUSY,1,36);
  double       tb= slhaVal("HMIX",QSUSY,1,2 );


  double CW     =MW/MZ ;// cos of the Weinberg angle;
  double SW     =sqrt(1-CW*CW) ;// sin of the Weinberg angle;
  double S2W    =2*SW*CW;
  double C2W    =CW*CW-SW*SW;
  double VEV    =2*MW*SW/EE;
  double GF     =(EE*EE)/pow(2*SW*MW,2)/Sqrt2 ;//  Fermi constant/ 1.166379E-5;
  double Lqcd   =initQCD5(alfSMZ, McMc, MbMb, Mtp);
  double pi     =acos(-1);
  double Mbp    =MbMb*(1+4./3.*alphaQCD(MbMb)/pi);
  double Mcp    =McMc*(1+4./3.*alphaQCD(McMc)/pi);
  double sb     =tb/sqrt(1+tb*tb) ;// Sine beta;
  double cb     =sqrt(1-sb*sb) ;// Cosine beta;
  double c2b    =cb*cb-sb*sb ;// cos(2b);
  double MSne   =sqrt(c2b*(MW*MW)/2+(CW*CW)*Ml1*Ml1)/CW ;// e-sNeutrino mass;
  double MSnm   =sqrt(c2b*(MW*MW)/2+(CW*CW)*Ml2*Ml2)/CW ;// mu-sneutrino mass;
  double MSeL   =sqrt(-c2b*((CW*CW)-(SW*SW))*(MW*MW)/2+(CW*CW)*Ml1*Ml1)/CW ;// 1st selectron mass;
  double MSeR   =sqrt(-(SW*SW)*c2b*(MW*MW)+(CW*CW)*(Mr1*Mr1))/CW ;// 2nd selectron mass;
  double MSmL   =sqrt(-c2b*((CW*CW)-(SW*SW))*(MW*MW)/2+(CW*CW)*Ml2*Ml2)/CW ;// 1st smuon mass;
  double MSmR   =sqrt(-(SW*SW)*c2b*(MW*MW)+(CW*CW)*Mr2*Mr2)/CW ;// 2nd smuon mass;
  double MSuL   =sqrt((Mq1*Mq1)+(MW*MW)/(CW*CW)*(1./2.-2./3.*(SW*SW))*c2b);
  double MScL   =sqrt((Mq2*Mq2)+(MW*MW)/(CW*CW)*(1./2.-2./3.*(SW*SW))*c2b);
  double MSuR   =sqrt((Mu1*Mu1)+(MW*MW)/(CW*CW)*(2./3.*(SW*SW))*c2b);
  double MScR   =sqrt((Mu2*Mu2)+(MW*MW)/(CW*CW)*(2./3.*(SW*SW))*c2b);
  double MSdL   =sqrt((Mq1*Mq1)-(MW*MW)/(CW*CW)*(1./2.-1./3.*(SW*SW))*c2b);
  double MSsL   =sqrt((Mq2*Mq2)-(MW*MW)/(CW*CW)*(1./2.-1./3.*(SW*SW))*c2b);
  double MSdR   =sqrt((Mu1*Mu1)+(MW*MW)/(CW*CW)*(-1./3.*(SW*SW))*c2b);
  double MSsR   =sqrt((Mu2*Mu2)+(MW*MW)/(CW*CW)*(-1./3.*(SW*SW))*c2b);
  double MSnl   =sqrt(c2b*(MW*MW)/2+(CW*CW)*(Ml3*Ml3))/CW ;// tau-sneutrino mass;
  double lDiag  =rDiagonal(2, Ml3*Ml3+Ml*Ml-(1./2.-(SW*SW))*MZ*MZ*c2b, Ml*(Al-mu*tb), Mr3*Mr3+Ml*Ml-(SW*SW)*MZ*MZ*c2b);

  double MSl1   =sqrt(MassArray(lDiag, 1));
  double MSl2   =sqrt(MassArray(lDiag, 2));
  double Zl11   =MixMatrix(lDiag, 1, 1);
  double Zl12   =MixMatrix(lDiag, 1, 2);
  double Zl21   =MixMatrix(lDiag, 2, 1);
  double Zl22   =MixMatrix(lDiag, 2, 2);
  double Mt11_  =(Mq3*Mq3)+(1./2.-2./3.*(SW*SW))*MZ*MZ*c2b;
  double Mt12_  =At-mu/tb;
  double Mt22_  =(Mu3*Mu3)+2./3.*(SW*SW)*(MZ*MZ)*c2b;
  double QSUSY_  =pow((Mt11_+(Mtp*Mtp))*(Mt22_+(Mtp*Mtp))-Mt12_*Mt12_*(Mtp*Mtp),0.25);
  double MtMM   =MtRun(QSUSY_);

  double tDiag  =rDiagonal(2, Mt11_+(MtMM*MtMM), MtMM*Mt12_, Mt22_+MtMM*MtMM);
  double MSt1   =sqrt(MassArray(tDiag, 1));
  double MSt2   =sqrt(MassArray(tDiag, 2));
  double Zt11   =MixMatrix(tDiag, 1, 1);
  double Zt12   =MixMatrix(tDiag, 1, 2);
  double Zt21   =MixMatrix(tDiag, 2, 1);
  double Zt22   =MixMatrix(tDiag, 2, 2);
  double MbMM   =MbRun(QSUSY_);
  double bDiag  =rDiagonal(2, (Mq3*Mq3)+MbMM*MbMM-(1./2.-1./3.*(SW*SW))*(MZ*MZ)*c2b, MbMM*(Ab-mu*tb), Md3*Md3+MbMM*MbMM-1/3*(SW*SW)*(MZ*MZ)*c2b);
  double MSb1   =sqrt(MassArray(bDiag, 1));
  double MSb2   =sqrt(MassArray(bDiag, 2));
  double Zb11   =MixMatrix(bDiag, 1, 1);
  double Zb12   =MixMatrix(bDiag, 1, 2);
  double Zb21   =MixMatrix(bDiag, 2, 1);
  double Zb22   =MixMatrix(bDiag, 2, 2);
  double MSG    =MG3;
  double zero   =initDiagonal();
         zero=0;
  double MNE[5];
  MNE[1]=slhaVal("MASS",QSUSY,1,1000022);
  MNE[2]=slhaVal("MASS",QSUSY,1,1000023);
  MNE[3]=slhaVal("MASS",QSUSY,1,1000025);
  MNE[4]=slhaVal("MASS",QSUSY,1,1000035);
  double Zn[5][5];
  for(int i=1;i<=4;i++) for(int j=1;j<=4;j++) Zn[i][j]= slhaVal("Nmix",QSUSY,2,i,j);
    
         
  double MG1_= Zn[1][1]*MNE[1]*Zn[1][1]+Zn[2][1]*MNE[2]*Zn[2][1]+Zn[3][1]*MNE[3]*Zn[3][1]+Zn[4][1]*MNE[4]*Zn[4][1];
  double MG2_= Zn[1][2]*MNE[1]*Zn[1][2]+Zn[2][2]*MNE[2]*Zn[2][2]+Zn[3][2]*MNE[3]*Zn[3][2]+Zn[4][2]*MNE[4]*Zn[4][2];        
// printf("MG1=%E MG1_=%E MG2=%E MG2_=%E\n", MG1,MG1_,MG2,MG2_);         
  MG1=MG1_;
  MG2=MG2_;

  double NeDiag =rDiagonal(4, MG1, zero, -MZ*SW*cb, MZ*SW*sb, MG2, MZ*CW*cb, -MZ*CW*sb, zero, -mu, zero);
  double MNE1   =MassArray(NeDiag, 1);
  double MNE2   =MassArray(NeDiag, 2);
  double MNE3   =MassArray(NeDiag, 3);
  double MNE4   =MassArray(NeDiag, 4);
  double Zn11   =MixMatrix(NeDiag, 1, 1);
  double Zn12   =MixMatrix(NeDiag, 1, 2);
  double Zn13   =MixMatrix(NeDiag, 1, 3);
  double Zn14   =MixMatrix(NeDiag, 1, 4);
  double Zn21   =MixMatrix(NeDiag, 2, 1);
  double Zn22   =MixMatrix(NeDiag, 2, 2);
  double Zn23   =MixMatrix(NeDiag, 2, 3);
  double Zn24   =MixMatrix(NeDiag, 2, 4);
  double Zn31   =MixMatrix(NeDiag, 3, 1);
  double Zn32   =MixMatrix(NeDiag, 3, 2);
  double Zn33   =MixMatrix(NeDiag, 3, 3);
  double Zn34   =MixMatrix(NeDiag, 3, 4);
  double Zn41   =MixMatrix(NeDiag, 4, 1);
  double Zn42   =MixMatrix(NeDiag, 4, 2);
  double Zn43   =MixMatrix(NeDiag, 4, 3);
  double Zn44   =MixMatrix(NeDiag, 4, 4);
  double chDiag =rDiagonalA(2, MG2, Sqrt2*MW*sb, Sqrt2*MW*cb, mu);
  double MC1    =MassArray(chDiag, 1);
  double MC2    =MassArray(chDiag, 2);
  double Zu11   =MixMatrixU(chDiag, 1, 1);
  double Zu12   =MixMatrixU(chDiag, 1, 2);
  double Zu21   =MixMatrixU(chDiag, 2, 1);
  double Zu22   =MixMatrixU(chDiag, 2, 2);
  double Zv11   =MixMatrix(chDiag, 1, 1);
  double Zv12   =MixMatrix(chDiag, 1, 2);
  double Zv21   =MixMatrix(chDiag, 2, 1);
  double Zv22   =MixMatrix(chDiag, 2, 2);

  MH3+=2*(fabs(MNE1)-fabs(MNE[1]));
  
  
  double calcL  =calcLambdas(tb, MH3, mu, QSUSY, Mtp, MSG, MC2, Mq3, Mu3, Md3, Ab, At);
  double la1    =Lambda1();
  double la2    =Lambda2();
  double la3    =Lambda3();
  double la4    =Lambda4();
  double la5    =Lambda5();
  double la6    =-Lambda6();
  double la7    =-Lambda7();
  double MHc    =sqrt((MH3*MH3)+2/(EE*EE)*(MW*MW)*(SW*SW)*(la5-la4));
  double Mhh11  =(VEV*VEV)*((cb*cb)*la1+(sb*sb)*la5-2*cb*la6*sb)+(MH3*MH3)*(sb*sb);
  double Mhh12  =(VEV*VEV)*(sb*cb*(la3+la4)-la6*(cb*cb)-la7*(sb*sb))-(MH3*MH3)*cb*sb;
  double Mhh22  =(VEV*VEV)*((sb*sb)*la2+(cb*cb)*la5-2*sb*cb*la7)+(MH3*MH3)*(cb*cb);
  double hDiag  =rDiagonal(2, Mhh11, Mhh12, Mhh22);
  double Mh     =sqrt(MassArray(hDiag, 1));
  double MH     =sqrt(MassArray(hDiag, 2));
  double Zh11   =MixMatrix(hDiag, 2, 1);
  double Zh12   =MixMatrix(hDiag, 2, 2);
  double Zh21   =MixMatrix(hDiag, 1, 1);
  double Zh22   =MixMatrix(hDiag, 1, 2);

{
  FILE *f;
  double Q;

  f=fopen("tree.slha","w");

/*
   fName2c(CODE,codeTxt,len1);

  fprintf(f,"Block SPINFO\n");
  sscanf(codeTxt,"%s %[^\n]",codeN,codeV);
  fprintf(f,"   1  %s\n",codeN);
  fprintf(f,"   2  %s\n",codeV);

  if(*err<0) fprintf(f,"   3  %d\n",-(*err)  );
  if(*err>0) {fprintf(f,"   4  %d\n", (*err)  ); fclose(f); return 1;}
*/
  Q=findValW("QSUSY");
  fprintf(f,"Block MASS   # Mass spectrum\n");
  fprintf(f,"#PDG code      mass           particle\n");
  fprintf(f,"       24   %16.8E    # MW\n",              80.423);
  fprintf(f,"       25   %16.8E    # h0\n",              Mh); 
  fprintf(f,"       35   %16.8E    # H0\n",              MH); 
  fprintf(f,"       36   %16.8E    # A0\n",              MH3); 
  fprintf(f,"       37   %16.8E    # H+\n",              MHc); 
  fprintf(f,"  1000001   %16.8E    # ~d_L\n",            MSdL); 
  fprintf(f,"  1000002   %16.8E    # ~u_L\n",            MSuL); 
  fprintf(f,"  1000003   %16.8E    # ~s_L\n",            MSsL); 
  fprintf(f,"  1000004   %16.8E    # ~c_L\n",            MScL); 
  fprintf(f,"  1000005   %16.8E    # ~b_1\n",            MSb1); 
  fprintf(f,"  1000006   %16.8E    # ~t_1\n",            MSt1); 
  fprintf(f,"  1000011   %16.8E    # ~e_L\n",            MSeL); 
  fprintf(f,"  1000012   %16.8E    # ~nue_L\n",          MSne); 
  fprintf(f,"  1000013   %16.8E    # ~mu_L\n",           MSmL); 
  fprintf(f,"  1000014   %16.8E    # ~numu_L\n",         MSnm); 
  fprintf(f,"  1000015   %16.8E    # ~stau_1\n",         MSl1); 
  fprintf(f,"  1000016   %16.8E    # ~nu_tau_L\n",       MSnl); 
  fprintf(f,"  1000021   %16.8E    # ~g\n",              MSG); 
  fprintf(f,"  1000022   %16.8E    # ~neutralino(1)\n",  MNE1); 
  fprintf(f,"  1000023   %16.8E    # ~neutralino(2)\n",  MNE2); 
  fprintf(f,"  1000024   %16.8E    # ~chargino(1)\n",    MC1); 
  fprintf(f,"  1000025   %16.8E    # ~neutralino(3)\n",  MNE3); 
  fprintf(f,"  1000035   %16.8E    # ~neutralino(4)\n",  MNE4); 
  fprintf(f,"  1000037   %16.8E    # ~chargino(2)\n",    MC2); 
  fprintf(f,"  2000001   %16.8E    # ~d_R\n",            MSdR); 
  fprintf(f,"  2000002   %16.8E    # ~u_R\n",            MSuR); 
  fprintf(f,"  2000003   %16.8E    # ~s_R\n",            MSsR); 
  fprintf(f,"  2000004   %16.8E    # ~c_R\n",            MScR); 
  fprintf(f,"  2000005   %16.8E    # ~b_2\n",            MSb2); 
  fprintf(f,"  2000006   %16.8E    # ~t_2\n",            MSt2); 
  fprintf(f,"  2000011   %16.8E    # ~e_R\n",            MSeR); 
  fprintf(f,"  2000013   %16.8E    # ~mu_R\n",           MSmR); 
  fprintf(f,"  2000015   %16.8E    # ~stau_2\n",         MSl2); 

  fprintf(f,"# Higgs mixing\n");
  fprintf(f,"Block ALPHA\n");
  fprintf(f,"     %16.8E\n",asin(Zh12));

  fprintf(f,"Block HMIX Q= %16.8E\n",Q);
  fprintf(f,"  1  %16.8E  # mu(Q)MSSM DRbar\n",mu);
  fprintf(f,"  2  %16.8E  # tan beta(Q)MSSM DRbar\n",tb);
  fprintf(f,"  3  %16.8E  # vev\n",  2.41971868E+02);
  
  fprintf(f,"Block MSOFT Q= %16.8E\n",Q);
    fprintf(f,"    1   %16.8E #MG1 \n", MG1  );  
    fprintf(f,"    2   %16.8E #MG2 \n", MG2  );  
    fprintf(f,"    3   %16.8E #MG3 \n", MG3  );  
    fprintf(f,"   31   %16.8E #Ml1 \n", Ml1  );  
    fprintf(f,"   32   %16.8E #Ml2 \n", Ml2  );  
    fprintf(f,"   33   %16.8E #Ml3 \n", Ml3  );  
    fprintf(f,"   34   %16.8E #Mr1 \n", Mr1  );  
    fprintf(f,"   35   %16.8E #Mr2 \n", Mr2  );  
    fprintf(f,"   36   %16.8E #Mr3 \n", Mr3  );  
    fprintf(f,"   41   %16.8E #Mq1 \n", Mq1  );  
    fprintf(f,"   42   %16.8E #Mq2 \n", Mq2  );  
    fprintf(f,"   43   %16.8E #Mq3 \n", Mq3  );  
    fprintf(f,"   44   %16.8E #Mu1 \n", Mu1  );  
    fprintf(f,"   45   %16.8E #Mu2 \n", Mu2  );  
    fprintf(f,"   46   %16.8E #Mu3 \n", Mu3  );  
    fprintf(f,"   47   %16.8E #Md1 \n", Md1  );  
    fprintf(f,"   48   %16.8E #Md2 \n", Md2  );  
    fprintf(f,"   49   %16.8E #Md3 \n", Md3  );  
    
  fprintf(f,"Block STOPMIX\n");
    fprintf(f," 1 1   %14E   # Zt11\n",Zt11);
    fprintf(f," 1 2   %14E   # Zt12\n",Zt12);
    fprintf(f," 2 1   %14E   # Zt21\n",Zt21);
    fprintf(f," 2 2   %14E   # Zt22\n",Zt22);
  
  fprintf(f,"Block SBOTMIX\n"); 
    fprintf(f," 1 1    %14E   # Zb11\n",Zb11);
    fprintf(f," 1 2    %14E   # Zb12\n",Zb12);
    fprintf(f," 2 1    %14E   # Zb21\n",Zb21);
    fprintf(f," 2 2    %14E   # Zb22\n",Zb22);


  fprintf(f,"Block STAUMIX\n");
    fprintf(f," 1 1    %14E   # Zl11\n",Zl11);
    fprintf(f," 1 2    %14E   # Zl12\n",Zl12);
    fprintf(f," 2 1    %14E   # Zl21\n",Zl21);
    fprintf(f," 2 2    %14E   # Zl22\n",Zl22);
    
  fprintf(f,"Block NMIX\n");  
  fprintf(f," 1 1    %14E   # Zn11\n",Zn11);
  fprintf(f," 1 2    %14E   # Zn12\n",Zn12);
  fprintf(f," 1 3    %14E   # Zn13\n",Zn13);
  fprintf(f," 1 4    %14E   # Zn14\n",Zn14);
  fprintf(f," 2 1    %14E   # Zn21\n",Zn21);
  fprintf(f," 2 2    %14E   # Zn22\n",Zn22);
  fprintf(f," 2 3    %14E   # Zn23\n",Zn23);
  fprintf(f," 2 4    %14E   # Zn24\n",Zn24);
  fprintf(f," 3 1    %14E   # Zn31\n",Zn31);
  fprintf(f," 3 2    %14E   # Zn32\n",Zn32);
  fprintf(f," 3 3    %14E   # Zn33\n",Zn33);
  fprintf(f," 3 4    %14E   # Zn34\n",Zn34);
  fprintf(f," 4 1    %14E   # Zn31\n",Zn41);
  fprintf(f," 4 2    %14E   # Zn32\n",Zn42);
  fprintf(f," 4 3    %14E   # Zn33\n",Zn43);
  fprintf(f," 4 4    %14E   # Zn34\n",Zn44);

  
  fprintf(f,"Block UMIX\n");
  fprintf(f," 1 1    %14E   # Zu11\n",Zu11);
  fprintf(f," 1 2    %14E   # Zu12\n",Zu12);
  fprintf(f," 2 1    %14E   # Zu21\n",Zu21);
  fprintf(f," 2 2    %14E   # Zu22\n",Zu22);

  fprintf(f,"Block VMIX\n");
  fprintf(f," 1 1    %14E   # Zv11\n",Zv11);
  fprintf(f," 1 2    %14E   # Zv12\n",Zv12);
  fprintf(f," 2 1    %14E   # Zv21\n",Zv21);
  fprintf(f," 2 2    %14E   # Zv22\n",Zv22);
  

  fprintf(f,"Block AU \n");
  fprintf(f,"  3 3 %16.8E # At\n",At);
  fprintf(f,"Block AD \n"); 
  fprintf(f,"  3 3 %16.8E # Ab\n",Ab);
  fprintf(f,"Block AE \n");
  fprintf(f,"  3 3 %16.8E # Al\n",Al);
  
/*  fprintf(f,"  2 2 %16.8E # Am\n",Am); */
  fprintf(f,"Block SMINPUTS\n"); 
  fprintf(f,"  1  %16.8E # 1/alfEMZ\n",1.27932432E+02); 
  fprintf(f,"  2  %16.8E # G_Fermi\n",1.16637E-5);
  fprintf(f,"  3  %16.8E # alfSMZ\n", alfSMZ); 
  fprintf(f,"  5  %16.8E # MbMb\n", MbMb);
  fprintf(f,"  6  %16.8E # Mtp\n", Mtp);
  fprintf(f,"  7  %16.8E # Mtau\n",Ml);
  fprintf(f,"Block MODSEL\n");
  fprintf(f,"  1    0    # General MSSM\n");
  fprintf(f,"Block MINPAR\n"); 
  fprintf(f,"  3  %16.8E # tb\n",tb); 
  fclose(f);
}
 int  err=lesHinput("tree.slha");
 if(err) return err;
 char cdmName[20];
 err=sortOddParticles(cdmName);
 if(err)  printf("Can't calculate %s\n",cdmName);
 return err;
    
}

int treeEwsbMSSM(void )
{
  double MZ= findValW("MZ");
  double MW= findValW("MW");  
  double EE=findValW("EE");
  double McMc=1.4;
  double Sqrt2=sqrt(2);

  
  double  alfSMZ=findValW("alfSMZ");     
  double  MbMb  =findValW("MbMb"); 
  double  Mtp   =findValW("Mtp"); 
  double  Ml    =findValW("Ml"); 

//  double QSUSY=sqrt(findValW("MASS",91,1,1000006)*findValW("MASS",91,1,2000006));
   
  double       MG1=findValW("MG1");
  double       MG2=findValW("MG2");
  double       MG3=findValW("MG3");
  double       Ml1=findValW("Ml1");
  double       Ml2=findValW("Ml2");
  double       Ml3=findValW("Ml3");
  double       Mr1=findValW("Mr1");
  double       Mr2=findValW("Mr2");
  double       Mr3=findValW("Mr3");
  double       Mq1=findValW("Mq1");
  double       Mq2=findValW("Mq2");
  double       Mq3=findValW("Mq3");
  double       Mu1=findValW("Mu1");
  double       Mu2=findValW("Mu2");
  double       Mu3=findValW("Mu3");
  double       Md1=findValW("Md1");
  double       Md2=findValW("Md2");
  double       Md3=findValW("Md3");

  double       At= findValW("At"  );
  double       Ab= findValW("Ab"  );
  double       Al= findValW("Al"  );

  double       mu= findValW("mu");
  double       MH3=findValW("MH3");
  double       tb= findValW("tb");


  double CW     =MW/MZ ;// cos of the Weinberg angle;
  double SW     =sqrt(1-CW*CW) ;// sin of the Weinberg angle;
  double S2W    =2*SW*CW;
  double C2W    =CW*CW-SW*SW;
  double VEV    =2*MW*SW/EE;
  double GF     =(EE*EE)/pow(2*SW*MW,2)/Sqrt2 ;//  Fermi constant/ 1.166379E-5;
  double Lqcd   =initQCD5(alfSMZ, McMc, MbMb, Mtp);
  double pi     =acos(-1);
  double Mbp    =MbMb*(1+4./3.*alphaQCD(MbMb)/pi);
  double Mcp    =McMc*(1+4./3.*alphaQCD(McMc)/pi);
  double sb     =tb/sqrt(1+tb*tb) ;// Sine beta;
  double cb     =sqrt(1-sb*sb) ;// Cosine beta;
  double c2b    =cb*cb-sb*sb ;// cos(2b);
  double MSne   =sqrt(c2b*(MW*MW)/2+(CW*CW)*Ml1*Ml1)/CW ;// e-sNeutrino mass;
  double MSnm   =sqrt(c2b*(MW*MW)/2+(CW*CW)*Ml2*Ml2)/CW ;// mu-sneutrino mass;
  double MSeL   =sqrt(-c2b*((CW*CW)-(SW*SW))*(MW*MW)/2+(CW*CW)*Ml1*Ml1)/CW ;// 1st selectron mass;
  double MSeR   =sqrt(-(SW*SW)*c2b*(MW*MW)+(CW*CW)*(Mr1*Mr1))/CW ;// 2nd selectron mass;
  double MSmL   =sqrt(-c2b*((CW*CW)-(SW*SW))*(MW*MW)/2+(CW*CW)*Ml2*Ml2)/CW ;// 1st smuon mass;
  double MSmR   =sqrt(-(SW*SW)*c2b*(MW*MW)+(CW*CW)*Mr2*Mr2)/CW ;// 2nd smuon mass;
  double MSuL   =sqrt((Mq1*Mq1)+(MW*MW)/(CW*CW)*(1./2.-2./3.*(SW*SW))*c2b);
  double MScL   =sqrt((Mq2*Mq2)+(MW*MW)/(CW*CW)*(1./2.-2./3.*(SW*SW))*c2b);
  double MSuR   =sqrt((Mu1*Mu1)+(MW*MW)/(CW*CW)*(2./3.*(SW*SW))*c2b);
  double MScR   =sqrt((Mu2*Mu2)+(MW*MW)/(CW*CW)*(2./3.*(SW*SW))*c2b);
  double MSdL   =sqrt((Mq1*Mq1)-(MW*MW)/(CW*CW)*(1./2.-1./3.*(SW*SW))*c2b);
  double MSsL   =sqrt((Mq2*Mq2)-(MW*MW)/(CW*CW)*(1./2.-1./3.*(SW*SW))*c2b);
  double MSdR   =sqrt((Mu1*Mu1)+(MW*MW)/(CW*CW)*(-1./3.*(SW*SW))*c2b);
  double MSsR   =sqrt((Mu2*Mu2)+(MW*MW)/(CW*CW)*(-1./3.*(SW*SW))*c2b);
  double MSnl   =sqrt(c2b*(MW*MW)/2+(CW*CW)*(Ml3*Ml3))/CW ;// tau-sneutrino mass;
  double lDiag  =rDiagonal(2, Ml3*Ml3+Ml*Ml-(1./2.-(SW*SW))*MZ*MZ*c2b, Ml*(Al-mu*tb), Mr3*Mr3+Ml*Ml-(SW*SW)*MZ*MZ*c2b);

  double MSl1   =sqrt(MassArray(lDiag, 1));
  double MSl2   =sqrt(MassArray(lDiag, 2));
  double Zl11   =MixMatrix(lDiag, 1, 1);
  double Zl12   =MixMatrix(lDiag, 1, 2);
  double Zl21   =MixMatrix(lDiag, 2, 1);
  double Zl22   =MixMatrix(lDiag, 2, 2);
  double Mt11_  =(Mq3*Mq3)+(1./2.-2./3.*(SW*SW))*MZ*MZ*c2b;
  double Mt12_  =At-mu/tb;
  double Mt22_  =(Mu3*Mu3)+2./3.*(SW*SW)*(MZ*MZ)*c2b;
  double QSUSY_  =pow((Mt11_+(Mtp*Mtp))*(Mt22_+(Mtp*Mtp))-Mt12_*Mt12_*(Mtp*Mtp),0.25);
  double MtMM   =MtRun(QSUSY_);

  double tDiag  =rDiagonal(2, Mt11_+(MtMM*MtMM), MtMM*Mt12_, Mt22_+MtMM*MtMM);
  double MSt1   =sqrt(MassArray(tDiag, 1));
  double MSt2   =sqrt(MassArray(tDiag, 2));
  double Zt11   =MixMatrix(tDiag, 1, 1);
  double Zt12   =MixMatrix(tDiag, 1, 2);
  double Zt21   =MixMatrix(tDiag, 2, 1);
  double Zt22   =MixMatrix(tDiag, 2, 2);
  double MbMM   =MbRun(QSUSY_);
  double bDiag  =rDiagonal(2, (Mq3*Mq3)+MbMM*MbMM-(1./2.-1./3.*(SW*SW))*(MZ*MZ)*c2b, MbMM*(Ab-mu*tb), Md3*Md3+MbMM*MbMM-1/3*(SW*SW)*(MZ*MZ)*c2b);
  double MSb1   =sqrt(MassArray(bDiag, 1));
  double MSb2   =sqrt(MassArray(bDiag, 2));
  double Zb11   =MixMatrix(bDiag, 1, 1);
  double Zb12   =MixMatrix(bDiag, 1, 2);
  double Zb21   =MixMatrix(bDiag, 2, 1);
  double Zb22   =MixMatrix(bDiag, 2, 2);
  double MSG    =MG3;
  double zero   =initDiagonal();
         zero=0;

  double NeDiag =rDiagonal(4, MG1, zero, -MZ*SW*cb, MZ*SW*sb, MG2, MZ*CW*cb, -MZ*CW*sb, zero, -mu, zero);
  double MNE1   =MassArray(NeDiag, 1);
  double MNE2   =MassArray(NeDiag, 2);
  double MNE3   =MassArray(NeDiag, 3);
  double MNE4   =MassArray(NeDiag, 4);
  double Zn11   =MixMatrix(NeDiag, 1, 1);
  double Zn12   =MixMatrix(NeDiag, 1, 2);
  double Zn13   =MixMatrix(NeDiag, 1, 3);
  double Zn14   =MixMatrix(NeDiag, 1, 4);
  double Zn21   =MixMatrix(NeDiag, 2, 1);
  double Zn22   =MixMatrix(NeDiag, 2, 2);
  double Zn23   =MixMatrix(NeDiag, 2, 3);
  double Zn24   =MixMatrix(NeDiag, 2, 4);
  double Zn31   =MixMatrix(NeDiag, 3, 1);
  double Zn32   =MixMatrix(NeDiag, 3, 2);
  double Zn33   =MixMatrix(NeDiag, 3, 3);
  double Zn34   =MixMatrix(NeDiag, 3, 4);
  double Zn41   =MixMatrix(NeDiag, 4, 1);
  double Zn42   =MixMatrix(NeDiag, 4, 2);
  double Zn43   =MixMatrix(NeDiag, 4, 3);
  double Zn44   =MixMatrix(NeDiag, 4, 4);
  double chDiag =rDiagonalA(2, MG2, Sqrt2*MW*sb, Sqrt2*MW*cb, mu);
  double MC1    =MassArray(chDiag, 1);
  double MC2    =MassArray(chDiag, 2);
  double Zu11   =MixMatrixU(chDiag, 1, 1);
  double Zu12   =MixMatrixU(chDiag, 1, 2);
  double Zu21   =MixMatrixU(chDiag, 2, 1);
  double Zu22   =MixMatrixU(chDiag, 2, 2);
  double Zv11   =MixMatrix(chDiag, 1, 1);
  double Zv12   =MixMatrix(chDiag, 1, 2);
  double Zv21   =MixMatrix(chDiag, 2, 1);
  double Zv22   =MixMatrix(chDiag, 2, 2);
  
  
  double calcL  =calcLambdas(tb, MH3, mu, sqrt(MSt1*MSt2), Mtp, MSG, MC2, Mq3, Mu3, Md3, Ab, At);
  double la1    =Lambda1();
  double la2    =Lambda2();
  double la3    =Lambda3();
  double la4    =Lambda4();
  double la5    =Lambda5();
  double la6    =-Lambda6();
  double la7    =-Lambda7();
  double MHc    =sqrt((MH3*MH3)+2/(EE*EE)*(MW*MW)*(SW*SW)*(la5-la4));
  double Mhh11  =(VEV*VEV)*((cb*cb)*la1+(sb*sb)*la5-2*cb*la6*sb)+(MH3*MH3)*(sb*sb);
  double Mhh12  =(VEV*VEV)*(sb*cb*(la3+la4)-la6*(cb*cb)-la7*(sb*sb))-(MH3*MH3)*cb*sb;
  double Mhh22  =(VEV*VEV)*((sb*sb)*la2+(cb*cb)*la5-2*sb*cb*la7)+(MH3*MH3)*(cb*cb);
  double hDiag  =rDiagonal(2, Mhh11, Mhh12, Mhh22);
  double Mh     =sqrt(MassArray(hDiag, 1));
  double MH     =sqrt(MassArray(hDiag, 2));
  double Zh11   =MixMatrix(hDiag, 2, 1);
  double Zh12   =MixMatrix(hDiag, 2, 2);
  double Zh21   =MixMatrix(hDiag, 1, 1);
  double Zh22   =MixMatrix(hDiag, 1, 2);

{
  FILE *f;
  double Q;

  f=fopen("tree.slha","w");

/*
   fName2c(CODE,codeTxt,len1);

  fprintf(f,"Block SPINFO\n");
  sscanf(codeTxt,"%s %[^\n]",codeN,codeV);
  fprintf(f,"   1  %s\n",codeN);
  fprintf(f,"   2  %s\n",codeV);

  if(*err<0) fprintf(f,"   3  %d\n",-(*err)  );
  if(*err>0) {fprintf(f,"   4  %d\n", (*err)  ); fclose(f); return 1;}
*/
  Q=sqrt(MSt1*MSt2);
  fprintf(f,"Block MASS   # Mass spectrum\n");
  fprintf(f,"#PDG code      mass           particle\n");
  fprintf(f,"       24   %16.8E    # MW\n",              80.423);
  fprintf(f,"       25   %16.8E    # h0\n",              Mh); 
  fprintf(f,"       35   %16.8E    # H0\n",              MH); 
  fprintf(f,"       36   %16.8E    # A0\n",              MH3); 
  fprintf(f,"       37   %16.8E    # H+\n",              MHc); 
  fprintf(f,"  1000001   %16.8E    # ~d_L\n",            MSdL); 
  fprintf(f,"  1000002   %16.8E    # ~u_L\n",            MSuL); 
  fprintf(f,"  1000003   %16.8E    # ~s_L\n",            MSsL); 
  fprintf(f,"  1000004   %16.8E    # ~c_L\n",            MScL); 
  fprintf(f,"  1000005   %16.8E    # ~b_1\n",            MSb1); 
  fprintf(f,"  1000006   %16.8E    # ~t_1\n",            MSt1); 
  fprintf(f,"  1000011   %16.8E    # ~e_L\n",            MSeL); 
  fprintf(f,"  1000012   %16.8E    # ~nue_L\n",          MSne); 
  fprintf(f,"  1000013   %16.8E    # ~mu_L\n",           MSmL); 
  fprintf(f,"  1000014   %16.8E    # ~numu_L\n",         MSnm); 
  fprintf(f,"  1000015   %16.8E    # ~stau_1\n",         MSl1); 
  fprintf(f,"  1000016   %16.8E    # ~nu_tau_L\n",       MSnl); 
  fprintf(f,"  1000021   %16.8E    # ~g\n",              MSG); 
  fprintf(f,"  1000022   %16.8E    # ~neutralino(1)\n",  MNE1); 
  fprintf(f,"  1000023   %16.8E    # ~neutralino(2)\n",  MNE2); 
  fprintf(f,"  1000024   %16.8E    # ~chargino(1)\n",    MC1); 
  fprintf(f,"  1000025   %16.8E    # ~neutralino(3)\n",  MNE3); 
  fprintf(f,"  1000035   %16.8E    # ~neutralino(4)\n",  MNE4); 
  fprintf(f,"  1000037   %16.8E    # ~chargino(2)\n",    MC2); 
  fprintf(f,"  2000001   %16.8E    # ~d_R\n",            MSdR); 
  fprintf(f,"  2000002   %16.8E    # ~u_R\n",            MSuR); 
  fprintf(f,"  2000003   %16.8E    # ~s_R\n",            MSsR); 
  fprintf(f,"  2000004   %16.8E    # ~c_R\n",            MScR); 
  fprintf(f,"  2000005   %16.8E    # ~b_2\n",            MSb2); 
  fprintf(f,"  2000006   %16.8E    # ~t_2\n",            MSt2); 
  fprintf(f,"  2000011   %16.8E    # ~e_R\n",            MSeR); 
  fprintf(f,"  2000013   %16.8E    # ~mu_R\n",           MSmR); 
  fprintf(f,"  2000015   %16.8E    # ~stau_2\n",         MSl2); 

  fprintf(f,"# Higgs mixing\n");
  fprintf(f,"Block ALPHA\n");
  fprintf(f,"     %16.8E\n",asin(Zh12));

  fprintf(f,"Block HMIX Q= %16.8E\n",Q);
  fprintf(f,"  1  %16.8E  # mu(Q)MSSM DRbar\n",mu);
  fprintf(f,"  2  %16.8E  # tan beta(Q)MSSM DRbar\n",tb);
  fprintf(f,"  3  %16.8E  # vev\n",  2.41971868E+02);
  
  fprintf(f,"Block MSOFT Q= %16.8E\n",Q);
    fprintf(f,"    1   %16.8E #MG1 \n", MG1  );  
    fprintf(f,"    2   %16.8E #MG2 \n", MG2  );  
    fprintf(f,"    3   %16.8E #MG3 \n", MG3  );  
    fprintf(f,"   31   %16.8E #Ml1 \n", Ml1  );  
    fprintf(f,"   32   %16.8E #Ml2 \n", Ml2  );  
    fprintf(f,"   33   %16.8E #Ml3 \n", Ml3  );  
    fprintf(f,"   34   %16.8E #Mr1 \n", Mr1  );  
    fprintf(f,"   35   %16.8E #Mr2 \n", Mr2  );  
    fprintf(f,"   36   %16.8E #Mr3 \n", Mr3  );  
    fprintf(f,"   41   %16.8E #Mq1 \n", Mq1  );  
    fprintf(f,"   42   %16.8E #Mq2 \n", Mq2  );  
    fprintf(f,"   43   %16.8E #Mq3 \n", Mq3  );  
    fprintf(f,"   44   %16.8E #Mu1 \n", Mu1  );  
    fprintf(f,"   45   %16.8E #Mu2 \n", Mu2  );  
    fprintf(f,"   46   %16.8E #Mu3 \n", Mu3  );  
    fprintf(f,"   47   %16.8E #Md1 \n", Md1  );  
    fprintf(f,"   48   %16.8E #Md2 \n", Md2  );  
    fprintf(f,"   49   %16.8E #Md3 \n", Md3  );  
    
  fprintf(f,"Block STOPMIX\n");
    fprintf(f," 1 1   %14E   # Zt11\n",Zt11);
    fprintf(f," 1 2   %14E   # Zt12\n",Zt12);
    fprintf(f," 2 1   %14E   # Zt21\n",Zt21);
    fprintf(f," 2 2   %14E   # Zt22\n",Zt22);
  
  fprintf(f,"Block SBOTMIX\n"); 
    fprintf(f," 1 1    %14E   # Zb11\n",Zb11);
    fprintf(f," 1 2    %14E   # Zb12\n",Zb12);
    fprintf(f," 2 1    %14E   # Zb21\n",Zb21);
    fprintf(f," 2 2    %14E   # Zb22\n",Zb22);


  fprintf(f,"Block STAUMIX\n");
    fprintf(f," 1 1    %14E   # Zl11\n",Zl11);
    fprintf(f," 1 2    %14E   # Zl12\n",Zl12);
    fprintf(f," 2 1    %14E   # Zl21\n",Zl21);
    fprintf(f," 2 2    %14E   # Zl22\n",Zl22);
    
  fprintf(f,"Block NMIX\n");  
  fprintf(f," 1 1    %14E   # Zn11\n",Zn11);
  fprintf(f," 1 2    %14E   # Zn12\n",Zn12);
  fprintf(f," 1 3    %14E   # Zn13\n",Zn13);
  fprintf(f," 1 4    %14E   # Zn14\n",Zn14);
  fprintf(f," 2 1    %14E   # Zn21\n",Zn21);
  fprintf(f," 2 2    %14E   # Zn22\n",Zn22);
  fprintf(f," 2 3    %14E   # Zn23\n",Zn23);
  fprintf(f," 2 4    %14E   # Zn24\n",Zn24);
  fprintf(f," 3 1    %14E   # Zn31\n",Zn31);
  fprintf(f," 3 2    %14E   # Zn32\n",Zn32);
  fprintf(f," 3 3    %14E   # Zn33\n",Zn33);
  fprintf(f," 3 4    %14E   # Zn34\n",Zn34);
  fprintf(f," 4 1    %14E   # Zn31\n",Zn41);
  fprintf(f," 4 2    %14E   # Zn32\n",Zn42);
  fprintf(f," 4 3    %14E   # Zn33\n",Zn43);
  fprintf(f," 4 4    %14E   # Zn34\n",Zn44);

  
  fprintf(f,"Block UMIX\n");
  fprintf(f," 1 1    %14E   # Zu11\n",Zu11);
  fprintf(f," 1 2    %14E   # Zu12\n",Zu12);
  fprintf(f," 2 1    %14E   # Zu21\n",Zu21);
  fprintf(f," 2 2    %14E   # Zu22\n",Zu22);

  fprintf(f,"Block VMIX\n");
  fprintf(f," 1 1    %14E   # Zv11\n",Zv11);
  fprintf(f," 1 2    %14E   # Zv12\n",Zv12);
  fprintf(f," 2 1    %14E   # Zv21\n",Zv21);
  fprintf(f," 2 2    %14E   # Zv22\n",Zv22);
  

  fprintf(f,"Block AU \n");
  fprintf(f,"  3 3 %16.8E # At\n",At);
  fprintf(f,"Block AD \n"); 
  fprintf(f,"  3 3 %16.8E # Ab\n",Ab);
  fprintf(f,"Block AE \n");
  fprintf(f,"  3 3 %16.8E # Al\n",Al);
  
/*  fprintf(f,"  2 2 %16.8E # Am\n",Am); */
  fprintf(f,"Block SMINPUTS\n"); 
  fprintf(f,"  1  %16.8E # 1/alfEMZ\n",1.27932432E+02); 
  fprintf(f,"  2  %16.8E # G_Fermi\n",1.16637E-5);
  fprintf(f,"  3  %16.8E # alfSMZ\n", alfSMZ); 
  fprintf(f,"  5  %16.8E # MbMb\n", MbMb);
  fprintf(f,"  6  %16.8E # Mtp\n", Mtp);
  fprintf(f,"  7  %16.8E # Mtau\n",Ml);
  fprintf(f,"Block MODSEL\n");
  fprintf(f,"  1    0    # General MSSM\n");
  fprintf(f,"Block MINPAR\n"); 
  fprintf(f,"  3  %16.8E # tb\n",tb); 
  fclose(f);
}

 int err=slhaRead("tree.slha",0);
printf("err=%d\n",err); 
 if(err)  cleanSLHAdata();
 return err;
    
}
