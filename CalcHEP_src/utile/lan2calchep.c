#include<stdlib.h>
#include<stdio.h>
#include<unistd.h>
#include<string.h>

static void wtrim(char * s)
{ int k,l;
  for(l=0; s[l]==' '      ; l++);
  for(k=l; s[k]&&s[k]!=' '; k++) s[k-l]=s[k];
  s[k-l]=0;      
}

static int mcnum(char* symb,char* name)
{
struct  {char *symb; int num;} mc[104]=
{
{"A",  22},
{"Z",  23},
{"W+", 24},
{"W-",-24},
{"G",  21},
{"e",  11},
{"e1", 11},
{"E", -11},
{"E1",-11},
{"m",  13},
{"e2", 13},
{"M", -13},
{"E2",-13},
{"ne", 12},
{"n1", 12},
{"Ne",-12},
{"N1",-12},
{"nm", 14},
{"n2", 14},
{"Nm",-14},
{"N2",-14},
{"nl", 16},
{"n3", 16},
{"Nl",-16},
{"N3",-16},
{"l",  15},
{"e3", 15},
{"L", -15},
{"E3",-15},
{"s",   3},
{"S",  -3},
{"c",   4},
{"C",  -4},
{"u",   2},
{"U",  -2},
{"d",   1},
{"D",  -1},
{"t",   6},
{"T",  -6},
{"b",   5},
{"B",  -5},
{"h",  25},
{"h1", 25},
{"H",  35},
{"h2", 35},
{"h3", 45},
{"H3", 36},
{"ha", 36},
{"hb", 46},
{"H+", 37},
{"H-",-37},
{"~1+", 1000024},
{"~1-",-1000024},
{"~2+", 1000037},
{"~2-",-1000037},
{"~o1", 1000022},
{"~o2", 1000023},
{"~o3", 1000025},
{"~o4", 1000035},
{"~o5", 1000045},
{"~G",  1000021},
{"~eL", 1000011},
{"~EL",-1000011},
{"~eR", 2000011},
{"~ER",-2000011},
{"~mL", 1000013},
{"~ML",-1000013},
{"~mR", 2000013},
{"~MR",-2000013},
{"~l1", 1000015},
{"~L1",-1000015},
{"~l2", 2000015},
{"~L2",-2000015},
{"~ne", 1000012},
{"~Ne",-1000012},
{"~nm", 1000014},
{"~Nm",-1000014},
{"~nl", 1000016},
{"~Nl",-1000016},
{"~uL", 1000002},
{"~UL",-1000002},
{"~uR", 2000002},
{"~UR",-2000002},
{"~cL", 1000004},
{"~CL",-1000004},
{"~cR", 2000004},
{"~CR",-2000004},
{"~t1", 1000006},
{"~T1",-1000006},
{"~t2", 2000006},
{"~T2",-2000006},
{"~dR", 2000001},
{"~DR",-2000001},
{"~dL", 1000001},
{"~DL",-1000001},
{"~sL", 1000003},
{"~SL",-1000003},
{"~sR", 2000003},
{"~SR",-2000003},
{"~b1", 1000005},
{"~B1",-1000005},
{"~b2", 2000005},
{"~B2",-2000005},
{"~g",  1000021}   };
 
int i;

for(i=0;i<104;i++) if (strcmp(mc[i].symb,symb)==0) return mc[i].num;


return 0;
}


#define BSIZE 5000
int main(int argc, char ** argv)
{  int from, to, i,j;
   FILE *fFrom, *fTo;
   char buff[BSIZE];
   char * part[4]={"vars","func","prtcls","lgrng"};
   int ok=1;

   if(argc!=4) ok=0; 
   else if(1!=(i=sscanf(argv[2],"%d",&from))) { ok=0;printf("i=%d\n",i);}
   else if(1!=(i=sscanf(argv[3],"%d",&to)  )) { ok=0;printf("i_=%d\n",i);}

   if(!ok)
   { printf("This programs needs 3 arguments\n"
     "1)path to destination of models in LanHEP format\n"
     "2)Model number N (like varsN.mdl)\n"
     "3)Target model number. The model will be created in the current directory\n");
     return 1;
   }        

   for(i=0;i<4;i++)
   {
     sprintf(buff,"%s/%s%d.mdl",argv[1],part[i],from);
     if(access(buff,R_OK))
     {  printf("Source file '%s' not found\n",buff);
        return 2;
     }
   }
   
   sprintf(buff,"cp %s/%s%d.mdl %s%d.mdl",argv[1],part[3],from,part[3],to);
   system(buff);
   
   for(i=0;i<3;i++)
   { int k=0;
     sprintf(buff,"%s/%s%d.mdl",argv[1],part[i],from);
     fFrom=fopen(buff,"r");
     sprintf(buff,"%s%d.mdl",part[i],to);
     fTo=fopen(buff,"w");
     fgets(buff,BSIZE,fFrom);  fputs(buff,fTo);
     fgets(buff,BSIZE,fFrom);  fputs(buff,fTo);

     switch(i)
     { 
       case 0:   /* vars*.md */
       for(;fgets(buff,BSIZE,fFrom);)
       { if( buff[0] != '=' ) fputs(" ",fTo);
         fputs(buff,fTo);
       }
       break;
       
       case 1:   /* func*.md */
       for(k=0;1==fscanf(fFrom,"%[^|]%*c",buff) ;)
       { if(buff[0]=='=') {fprintf(fTo,"%s\n",buff); break;}
         fprintf(fTo," %s|",buff);
         fscanf(fFrom,"%[^|]%*c",buff); 
         if(k==0) {fprintf(fTo,"%s&",buff); k=1;} else fprintf(fTo,"%s%%",buff);
         fgets(buff,BSIZE,fFrom); fputs(buff,fTo);
       }
       break;
       
       case 2:  /* Particles */
       for(j=0;j<3;j++){ fscanf(fFrom,"%[^|]%*c",buff); fprintf(fTo,"%s|",buff);}
       fprintf(fTo," number |");
       fgets(buff,BSIZE,fFrom);  fputs(buff,fTo);       
  
       for(;;)
       { char name[40],symb[40],aux[40];
         int nn;
         if( fscanf(fFrom,"%[^|%]%*c",buff)==EOF ) break;
         fprintf(fTo,"%s",buff);
         if(buff[0]=='=')break; else fprintf(fTo,"|");
         strcpy(name,buff); 
         
         fscanf(fFrom,"%[^|]%*c",symb); fprintf(fTo,"%s|",symb); 
         fscanf(fFrom,"%[^|]%*c",buff); fprintf(fTo,"%s|",buff);
         
         wtrim(symb);
         nn=mcnum(symb,name);
         fgets(buff,BSIZE,fFrom);
         
         if(nn==0)
         { sscanf(buff,"%*[^|]%*c%*[^|]%*c%*[^|]%*c%*[^|]%*c%[^|]",aux);
           wtrim(aux);
           if(strcmp(aux,"*")) printf("Warning! Monte Caro code for %s (%s) is unknown.\n"
                               " Replaced by zero. Improve it!\n",symb, name);
         }   
           
         fprintf(fTo,"%8d|",nn);
         fputs(buff,fTo);         
       }
       break; 
     }   
     fclose(fFrom); fclose(fTo);
   }

   sprintf(buff,"extlib%d.mdl",to);
   fTo=fopen(buff,"w");
   sprintf(buff,"%s/%s%d.mdl",argv[1],part[0],from);
   fFrom=fopen(buff,"r");
   fscanf(fFrom,"%[^\n]",buff);
   fprintf(fTo,"%s\n",buff);
   fprintf(fTo,"Libraries\n"
 "External libraries  and citation                                     <|\n"
 "=======================================================================\n");
   fclose(fTo),fclose(fFrom);          
   return 0;
}
