#include<stdlib.h>
#include<string.h>
#include<stdio.h> 

#include "polint.h"

int main(int argc, char** argv)
{ 
   if(argc<2) { printf("One argument is expected,  name  of directory with  LHAPDF distribution.\n");
                printf("Second  argument if presented, corresponds to set number.\n");
                exit(1);
              }
              
   char *dirname=malloc(strlen(argv[1])+1);
   strcpy(dirname,argv[1]);
   char *ch=dirname+strlen(dirname)-1;
   if(ch[0]=='/') { ch[0]=0; ch--;}
   while(ch>dirname && ch[0]!='/') ch--;
   if(ch[0]=='/') ch++;
   char* finfo=malloc(strlen(dirname)+1+strlen(ch)+6);
   sprintf(finfo,"%s/%s.info",dirname,ch);
   char* fdat=malloc(strlen(dirname)+1+strlen(ch)+10);
   sprintf(fdat,"%s/%s_",dirname,ch);
   int memb=0;
   
   if(argc==2) strcat(fdat,"0000.dat");
   else  
   {  int i;
      for(i=0;i< 4-strlen(argv[2]);i++) strcat(fdat,"0");
       strcat(fdat,argv[2]);
       strcat(fdat,".dat");
       sscanf(argv[2],"%d", &memb); 
   }                
   FILE*Finf=fopen(finfo,"r");
   FILE*Fdat=fopen(fdat,"r");
   if(!Finf || !Fdat){ printf("Can not open files\n"); exit(2);}

   char key[50],val[20];
   

// Parameters which we read form .info file 
   int Index=0;
   char Ref[100]={""};
   char Format[100]={""};
   int  Particle=0;
   int  Adim=0;
   double *Al=NULL;
   double *Qs=NULL;
//==============
   int i;

// Read info file    
   while(1==fscanf(Finf,"%[^:]:",key))
   {
//     printf("key='%s'\n", key);
     if(0==strcmp(key,"SetIndex"))        fscanf(Finf,"%d",&Index); 
     else if(0==strcmp(key,"Reference"))  fscanf(Finf,"%[^\n]",Ref);        
     else if(0==strcmp(key,"Format"))     
     {
       fscanf(Finf,"%s",Format); 
       if(0!=strcmp(Format,"lhagrid1")) { printf("Unknown format '%s'\n",Format); exit(3);}
     }
     else if(0==strcmp(key,"AlphaS_Qs")) 
     {  fscanf(Finf,"%*[^[][");
         Adim=0;
        double al;       
         while(1==fscanf(Finf,"%lf",&al))
         {  Qs=(double*) realloc(Qs, (Adim+1)*sizeof(double));
            Qs[Adim++]=al;
            fscanf(Finf," ,");
         } 
     }
     else if(0==strcmp(key,"AlphaS_Vals")) 
     {  fscanf(Finf,"%*[^[][");
        Adim=0;
        double al;       
         while(1==fscanf(Finf,"%lf",&al))
         {  Al=(double*) realloc(Al, (Adim+1)*sizeof(double));
            Al[Adim++]=al;
            fscanf(Finf," ,");
         } 
     }
     else if(0==strcmp(key,"Particle")) fscanf(Finf," %d",&Particle);
     char c;
     do  fscanf(Finf,"%c",&c); while(c!='\n');
   }
   fclose(Finf);
   
   printf(" Index=%d\n", Index);
   printf(" Set=%d\n",memb);
   printf(" Particle=%d\n", Particle);
   printf(" Format=%s\n", Format);
   printf(" Reference=%s\n", Ref);
   if(Particle!=2212){printf("Proton is expected\n"); exit(4);} 
   

//  Read data file        
   for(;;)
   { 
     fscanf(Fdat,"%s ",key);
     if(key[0]=='-') break;
     if(strcmp(key,"Format")==0)
     {  fscanf(Fdat,"%s ",val);
        if(strcmp(val,"lhagrid1")!=0)
        { printf("Unknown LHAPDF fprmat\n");
          exit(3);
        }
     }     
   }

   long posB=ftell(Fdat);
   fscanf(Fdat,"%*[^\n]");
   long posE=ftell(Fdat);
   
   fseek(Fdat,posB,SEEK_SET);
   int Xdim=0;
   double *X=NULL;

   for(;ftell(Fdat)<posE;)
   {
     X=(double*)realloc(X,(Xdim+1)*sizeof(double));
     fscanf(Fdat,"%lf ", X+Xdim);
     Xdim++;
   }  

   posB=ftell(Fdat);
   fscanf(Fdat,"%*[^\n]");
   posE=ftell(Fdat);   

   
   int Qdim=0;
   double *Q=NULL;
   for(fseek(Fdat,posB,SEEK_SET);ftell(Fdat)<posE;)
   { 
     Q=(double*)realloc(Q,(Qdim+1)*sizeof(double));
     fscanf(Fdat,"%lf ", Q+Qdim);
     Qdim++;
   }
   
   
   int Pdim=0;
   int  *P=NULL;
   
   posB=ftell(Fdat);
   fscanf(Fdat,"%*[^\n]");
   posE=ftell(Fdat);
   
   for(fseek(Fdat,posB,SEEK_SET);ftell(Fdat)<posE;)
   { 
       P=(int*)realloc(P,(Pdim+1)*sizeof(int));
       fscanf(Fdat,"%d ", P+Pdim);
       Pdim++;
   }             
 
   double *buff=(double*)malloc(Xdim*Qdim*Pdim*sizeof(double));

   for(i=0;i<Xdim*Qdim*Pdim;i++) if(1!=fscanf(Fdat," %lf ",buff+i))
   {  
      printf(" Unexpected end of file\n"); exit(5);     
   }
   
   fscanf(Fdat," %s ",key);
   printf("last key=%s\n",key);
   fclose(Fdat);


//  Output 
   int ip,iq,ix;

   if(argc>2)
   { char*ch_=malloc(strlen(ch)+strlen(argv[2])+2);
     sprintf(ch_,"%s:%s",ch,argv[2]);
     ch=ch_; 
   }
   
   char*fout=malloc(strlen(ch)+5);
   sprintf(fout,"%s.pdt",ch);   
   FILE*Fout=fopen(fout,"w");

   fprintf(Fout,"\n#distribution \"%s(proton)\"       2212 => ", ch);
   posB=ftell(Fout); for(i=0;i<Pdim;i++)  fprintf(Fout,"    ");    
   fprintf(Fout,"\n#distribution \"%s(anti-proton)\" -2212 => ", ch);
   posE=ftell(Fout); for(i=0;i<Pdim;i++)  fprintf(Fout,"    ");

   fprintf(Fout,"\n#Index %d",Index);
   fprintf(Fout,"\n#Memb %d",memb);
   fprintf(Fout,"\n#Source  LHAPDF6");
   fprintf(Fout,"\n#Reference %s",Ref); 
   fprintf(Fout,"\n#Interpolation biCubicLogXQ");
   
   fprintf(Fout,"\n#Q_grid\n");
   for(iq=0;iq<Qdim;iq++)
   { fprintf(Fout," %.4E",Q[iq]);
     if( (iq+1)%10==0) fprintf(Fout,"\n");
   }

   fprintf(Fout,"\n#Alpha\n");
   for(iq=0;iq<Qdim;iq++)
   { fprintf(Fout," %.4E",  polint3(Q[iq] , Adim, Qs, Al));
     if( (iq+1)%10==0) fprintf(Fout,"\n");
   }

   fprintf(Fout,"\n#X_grid\n");
   for(ix=0;ix<Xdim;ix++)
   { fprintf(Fout," %.4E",X[ix]);
     if( (ix+1)%10==0) fprintf(Fout,"\n");
   }
   
   char * Pstr=(char*)malloc(Pdim*4); Pstr[0]=0;
   char *aPstr=(char*)malloc(Pdim*4);aPstr[0]=0;
   int pk;
   for(ip=0,pk=0;ip<Pdim;ip++)
   {  int wrt=1;
      switch(P[ip])
      { case 21: case 22: case 23: 
            sprintf(Pstr+strlen(Pstr)," %d",P[ip]);
            sprintf(aPstr+strlen(aPstr)," %d",P[ip]);
            break;
        case 3: case 4: case 5: case 6: case 24:
            sprintf(Pstr+strlen(Pstr)," (%d %d)",P[ip],-P[ip]);
            sprintf(aPstr+strlen(aPstr)," (%d %d)",P[ip],-P[ip]);
            break;
        case -2: case -1: case 1: case 2: 
            sprintf(Pstr+strlen(Pstr)," %d",P[ip]);
            sprintf(aPstr+strlen(aPstr)," %d",-P[ip]);
            break;
        default : wrt=0;    
      }
      if(!wrt) continue;
      pk++; 
      fprintf(Fout,"\n#%d-parton\n",pk);
      for(iq=0;iq<Qdim;iq++)
      {  double d;
         for(ix=0;ix<Xdim;ix++) { fprintf(Fout," %.4E",buff[ip+ iq*Pdim +ix*Pdim*Qdim] );} 
         fprintf(Fout,"\n");
      }    
   } 
   fseek(Fout,posB,SEEK_SET); fprintf(Fout,"%s",Pstr); 
   fseek(Fout,posE,SEEK_SET); fprintf(Fout,"%s",aPstr); 
        
   fclose(Fout);
   return 0;   
}
