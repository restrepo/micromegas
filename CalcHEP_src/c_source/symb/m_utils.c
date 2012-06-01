/*
 Copyright (C) 1997, Alexander Pukhov 
*/

#include <unistd.h>
#include <stdlib.h>
#include <dirent.h>
#include <pwd.h>   

#include"chep_crt.h"
#include "syst2.h"
#include "physics.h"
#include "screen.h"
#include "s_files.h"
#include "file_scr.h"
#include "read_mdl.h"

#include "m_utils.h"


#define nhardmdl 4
#define newmodeltxt "   IMPORT OF MODELS   "

#define MAX_MODEL_NUM  ((STRSIZ-2)/22 -2)
  
void fillModelMenu(void)
{ int i,j;
  FILE * txt;
  char name[80];
  strcpy(modelmenu,"\026");
  maxmodel=0;
  for (i=1;i<=MAX_MODEL_NUM; i++)
  {  char fname[100];
     sprintf(fname,"./models%c%s%d.mdl",f_slash,mdFls[0],i);     
     txt=fopen(fname,"r");		
     if(txt==NULL) break;
     fgets(name,60,txt);
     trim(name);
     for (j=strlen(name); j<21;j++) name[j]=' ';
     name[21]=0;
     strcat(modelmenu," ");
     strcat(modelmenu,name);
     fclose(txt);
     maxmodel++;
  }  
  if ( maxmodel < MAX_MODEL_NUM ) strcat(modelmenu,newmodeltxt) ;
}



int  deletemodel(int n)
{
  int   i;
  char from[100],to[100];

  sprintf(from,"%smodels%c%s%d.mdl",pathtocomphep,f_slash,mdFls[0],n);
 
  if( mess_y_n(56,10,"Delete model?"))
  {
     for(i=0;i<5;i++)
     {
        sprintf(to,"%smodels%c%s%d.mdl",pathtouser,f_slash,mdFls[i],n);
	unlink(to);
     }
     for(; n < maxmodel;n++)
     {
        for(i=0;i<5;i++)
        {
          sprintf(from,"%smodels%c%s%d.mdl",pathtouser,f_slash,mdFls[i],n+1);
          sprintf(to,  "%smodels%c%s%d.mdl",pathtouser,f_slash,mdFls[i],n);
          rename(from,to);
        }
     }
     return 1;

     return 1;
  } 
  return 0;
}


static char * UnixNameConv(char * sIn)
{ 
  static struct passwd *passwdPrt;
  if(sIn[0]=='$') return getenv(sIn+1);
  if(strcmp(sIn,"~")==0)   return getenv("HOME");
  if(sIn[0]!='~') return NULL;
  passwdPrt=getpwnam(sIn+1);
  if(passwdPrt) return passwdPrt->pw_dir; else return 0;   
}


int  makenewmodel(void)
{ static shortstr  dirName="";
  int       key, nmdl=0;
  int  n, i;
  DIR *dirPtr;
  FILE*f,*g;
  struct dirent * dp;  
  char new[4000];
  void *pscr = NULL;

  char fname[STRSIZ];
  char buff[STRSIZ];

  for(;;)
  {  goto_xy(1,7); scrcolor(Red,BGmain);
     print("Enter name of directory with models, or press Esc to Exit, or F1 for Help\n");
     goto_xy(1,9); scrcolor(Yellow,Blue);
     print("Dir= "); 
     key = str_redact(dirName,1,72);
     scrcolor(FGmain,BGmain);
     if(key==KB_ESC) { clrbox(1,7,80,12); return 1;}
     if(key==KB_F1) show_help("s_addmodel");
     if(key==KB_ENTER)
     {  char Path[STRSIZ];
        trim(dirName);
        if(dirName[0]=='~' ||dirName[0]=='$')
        {  char * upath,* cslash;
           cslash=strchr(dirName,'/');
           if(cslash) cslash[0]=0;
           upath=UnixNameConv(dirName);
           if(cslash){cslash[0]='/';cslash++;}else cslash="";
           if(upath)  sprintf(Path,"%s/%s",upath,cslash);else 
           {  messanykey(10,10,"Can not interpret the beginning of Path");
              continue;
           }
        } else strcpy(Path,dirName);
        strcpy(new,"\040");
        trim(dirName);
        dirPtr=opendir(Path);
        if(!dirPtr)
        {  messanykey(10,10,"Such directory is absent");
           continue;
        }
        while(dp=readdir(dirPtr))
        { char tail[100];
          if(sscanf(dp->d_name,"vars%d.%s",&n,tail)==2 && n<100 && strcmp(tail,"mdl")==0)
          { 
             for(i=3;i>=0;i--)
             {   
                sprintf(fname,"%s/%s%d.mdl",Path,mdFls[i],n);
                f=fopen(fname,"r");
                if(!f) break;
                if(i==0)
                {  fscanf(f,"%[^\n]",buff);
                   trim(buff);
                   buff[22]=0;
                   sprintf(new+strlen(new),"%22.22s |%3d.mdl ",buff,n);
                }
                fclose(f);
             }
          }
        }
     
        if(strlen(new)==1) 
        { messanykey(12,12," This directory does not contain\n"
                           " model files");
          continue;
        }
        menu1(5,12,"Choose a model",new,NULL,&pscr,&nmdl);
        if(nmdl==0) continue;

        sscanf(new+1+new[0]*(nmdl-1), "%22c |%d",buff,&n);
        buff[22]=0; trim(buff);
        correctStr(5,16,"Correct name ",buff,22,1);
        trim(buff);
        for(i=0;i<5;i++)
        {  int s; char ch[1000];
           sprintf(fname,"%s/%s%d.mdl",Path,mdFls[i],n);
           f=fopen(fname,"r");
           sprintf(fname,"models/%s%d.mdl",mdFls[i],maxmodel+1);
           g=fopen(fname,"w"); fprintf(g,"%s\n",buff);
           if(f==NULL && i==4)
           { fprintf(g,"Libraries\n");
             fprintf(g,"External libraries  and citation                                      <|\n");
           }else 
           {
             fscanf(f,"%*[^\n]%*c");
             while( s=fread(ch,1,1000,f)) fwrite(ch,1,s,g);
             fclose(f); 
           }  
           fclose(g);  
        }
        fillModelMenu();
        messanykey(10,16,"The model is added"); 
        put_text(&pscr);
     }     
  }
}



int  continuetest(void)
{shortstr  txt;
 int  k, ndel, ncalc, nrest; 
 long recpos;

   menuq=fopen(MENUQ_NAME,"rb");
  
   for(k=1;k<=subproc_sq;k++)
   {
      rd_menu(2,k,txt,&ndel,&ncalc,&nrest,&recpos);
      if (nrest != 0)
      {
         fclose(menuq);
         return 0;
      }
   }
   fclose(menuq);
   return 1;
}


void  clear_tmp(void)
{  int i;
   char name[40];
   for(i=1;;i++) { sprintf(name,ARCHIV_NAME,i);    if(unlink(name)) break;}
   for(i=1;;i++) { sprintf(name,ARCHIV_NAME "2",i);if(unlink(name)) break;}  
   unlink(CATALOG_NAME);
   unlink(CATALOG_NAME "2");
}
