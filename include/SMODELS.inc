{  int n1,n2;
   char smodels_database[100];
   char smodels_version[100];
   pythonversion_(&n1,&n2);
   system("rm -f smodels.*");
   if(n1==2)
   {
     printf("\nSMODELS:\n");
     smodels(5/*nf*/, 0.1/*fb*/, "smodels.in",0);
     system("make -C ../Packages/smodels-v1.1.0patch1 1>/dev/null"); 
     system("../Packages/smodels-v1.1.0patch1/smodelsTools.py -particles ./ xseccomputer -p -s 8 -N -v warning  -O -f smodels.in 2>smodels.out");
     system("../Packages/smodels-v1.1.0patch1/runSModelS.py -f smodels.in -o ./ -particles ./ -v error 2>> smodels.out");
     slhaRead("smodels.in.smodelsslha", 1);
     slhaSTRFormat("SModelS_Settings"," 0 %s",smodels_version);
     slhaSTRFormat("SModelS_Settings"," 1 %s",smodels_database);
     printf("%s with database %s \n",smodels_version, smodels_database);
     result=slhaVal("SModelS_Exclusion",0.,2,0,0); 
     if(result>=0)
     { 
       slhaSTRFormat("SModelS_Exclusion"," 1 0  %s",topology);
       Rvalue=slhaVal("SModelS_Exclusion",0.,2,1,1);
       slhaSTRFormat("SModelS_Exclusion"," 1 4 %s",analysis);
       printf("highest R=%.2E from %s, topology %s\n", Rvalue,analysis,topology);
     }
     switch(result)
     { case -1: printf("no channels for testing.\n");break;
       case  0: printf("not excluded.\n");break; 
       case  1: printf("excluded.\n ");break;
     }
   } else 
   { printf("SMODELS can not work with Python version > 2\n"); 
     result=0;
     Rvalue=0;
   }  
}
