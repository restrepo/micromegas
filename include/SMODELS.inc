{
   char smodels_database[100];
   char smodels_version[100];
   system("rm -f smodels.*");
   char *VERSION="1.2.2";
 
   printf("\nSMODELS:\n");
   smodels(LHCrun,5/*nf*/, 0.1/*fb*/,  "smodels.in",0);
   
   char*command=(char*)malloc(strlen(micrO)+200 );
    
   sprintf(command,"%s/Packages/smodels-%s/runSModelS.py",micrO,VERSION);
   if(access(command,X_OK ))
   {
      sprintf(command,"cd  %s/Packages;  make -f  SMODELS.makef VERSION=%s", micrO,VERSION); 
      system(command);
   }   

//   sprintf(command,"micrO=%s; python3  $micrO/Packages/smodels-%s/smodelsTools.py xseccomputer -p -N -v warning  -O -f smodels.in 2>smodels.out",micrO,VERSION),
//   system(command);
//printf("%s\n", command);
   sprintf(command,"micrO=%s; python3  $micrO/Packages/smodels-%s/runSModelS.py -f smodels.in -p $micrO/include/micromegasInterfaceParams.ini -o ./ -v error 2>>smodels.out",micrO,VERSION);
   system(command);
//printf("%s\n", command);   
   free(command);   
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
}
