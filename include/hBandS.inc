   if(access( "../Packages/HiggsSignals/HiggsSignals",X_OK )) system("cd ../Packages;  make -f hBandS.make");
   system("rm -f HB.out HS.out");
   system("cp HB.in HS.in");
   char buff[200];
#ifdef HIGGSBOUNDS
   sprintf(buff,"../Packages/HiggsBounds/HiggsBounds  LandH SLHA %d  %d HB.in HB.out > hb.stdout",NH0,NHch);
   system(buff);


   slhaRead("HB.out",1+4);
   slhaSTRFormat("HiggsBoundsResults","0 %s %*s",HB_version);

   HB_result=slhaValFormat("HiggsBoundsResults",0.,"1 2 %lf");
   HB_obsratio=slhaValFormat("HiggsBoundsResults",0.,"1 3 %lf" );
   slhaSTRFormat("HiggsBoundsResults","1 5 ||%[^|]||",HB_chan);
#endif 

#ifdef HIGGSSIGNALS
#define DataSet " latestresults "
//#define Method  " peak " 
//#define  Method " mass "
#define  Method " both "
#define PDF  " 2 "  // Gaussian
//#define PDF " 1 "  // box 
//#define PDF " 3 "  // box+Gaussia

  sprintf(buff,"../Packages/HiggsSignals/HiggsSignals "  DataSet Method   PDF  " SLHA %d %d HS.in > hs.stdout",NH0,NHch); 
//#../Packages/HiggsSignals-1.3.0/HiggsSignals   latestresults peak  2 SLHA 3 1 HS.in  > hs.stdout
  system(buff);

  system("grep -A 10000  HiggsSignalsResults HS.in > HS.out");
     slhaRead("HS.out",1+4);
     slhaSTRFormat("HiggsSignalsResults"," 0 ||%[^|]||",HS_version);
     HS_observ=slhaVal("HiggsSignalsResults",0.,1,7);
     HS_chi2=slhaVal("HiggsSignalsResults",0.,1,12);
     HS_pval=slhaVal("HiggsSignalsResults",0.,1,13);
#undef PDF
#undef Method
#undef DataSet
#endif
