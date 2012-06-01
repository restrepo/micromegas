/*
 Copyright (C) 1997, Alexander Pukhov 
*/
#include"4_vector.h"
#include<math.h>
#include"phys_val.h"
#include "interface.h"
#include"subproc.h"
#include<ctype.h>
#include"kinaux.h"
#include"usrfun.h"
#include "const.h"
#include "syst.h"

double calcPhysVal(char key,char * lv)
{
   double p1,p2,p3, q1,q2,q3, mp,mq ,cs,dl; 
   int j;


   int np1=4*(lv[0]-1), np2=4*(lv[1]-1);
   double s=0;
   int i=0;
   double pp[4]={0,0,0,0};
 
   switch(key)
   {  case 'A':
      case 'C': 
            p1 = pvect[np1+1];
            p2 = pvect[np1+2];
            p3 = pvect[np1+3];
            if(np2==-4) np2=0;
            q1 = pvect[np2+1];
            q2 = pvect[np2+2];
            q3 = pvect[np2+3];

            cs = (p1*q1+p2*q2+p3*q3)/
                sqrt( (p1*p1+p2*p2+p3*p3) * (q1*q1+q2*q2+q3*q3) );
            if (key == 'A') return  acos(cs)*180/M_PI; 
            return cs;
        case 'J':
            p1 = pvect[np1+1];
            p2 = pvect[np1+2];
            p3 = pvect[np1+3];
            mp = sqrt(p1*p1 + p2*p2 + p3*p3);

            q1 = pvect[np2+1];
            q2 = pvect[np2+2];
            q3 = pvect[np2+3];
            mq = sqrt(q1*q1 + q2*q2 + q3*q3);
            
            cs = (p1*q1 + p2*q2)/sqrt( (p1*p1 + p2*p2)*(q1*q1 + q2*q2) );
            cs = acos(cs);
                   
            dl = (mp + p3) * (mq - q3)  / (mp - p3) / (mq +q3);
            dl = 0.5*log(dl);
            
            return sqrt(dl*dl + cs*cs);

        case 'P':   
	{  
	   double mtot, mtot2, ms,md, pcm,p2;
	 
	   pinf_int(Nsub, lv[0], &mp, NULL);
           for(i=1;lv[i];i++)
           for(j=0;j<4;j++) pp[j]+=pvect[4*(lv[i]-1)+j];

	   p2 = pp[1]*pp[1] + pp[2]*pp[2] + pp[3]*pp[3];
	   mq=sqrt(fabs(pp[0]*pp[0]-p2));

	   for(j=0;j<4;j++) pp[j]+=pvect[np1+j];
	   p2 = pp[1]*pp[1] + pp[2]*pp[2] + pp[3]*pp[3];
	   
           mtot2 = fabs(pp[0]*pp[0] -p2);
           mtot=sqrt(mtot2);
           
	   ms = mp + mq;
	   md = mp - mq;

	   pcm = sqrt( (mtot2-ms*ms) * (mtot2-md*md) )/(2*mtot);
	   s =pp[1]*pvect[np1+1]+pp[2]*pvect[np1+2]+pp[3]*pvect[np1+3];
              
           return  (s*pp[0] - pvect[np1]*p2 ) /(sqrt(p2) * mtot * pcm);

        }
        case 'E':
            while(lv[i]!=0)   s += pvect[  (lv[i++]<<2)  -4];
            return s;

        case 'T':
            {  i=0; 
               do for(j=1;j<3;j++) pp[j] += pvect[4*(lv[i]-1)+j]; while(lv[++i]);
               return sqrt( pp[1]*pp[1]+pp[2]*pp[2]);        
            } 
        case 'Z':
            {  double E=0, Pl=0; 
               i=0; 
               do 
               { E += pvect[4*(lv[i]-1)];
                 Pl+= pvect[4*(lv[i]-1)+3]; 
               }
               while(lv[++i]);
               if(nin_int==1) return E; else return sqrt(E*E-Pl*Pl);
            }
        case 'S':
        case 'M': 
            {   
               do 
               {
                  if(lv[i]>nin_int) for(j=0;j<4;j++) pp[j] += pvect[4*(lv[i]-1)+j];
                  else for(j=0;j<4;j++)           pp[j] -= pvect[4*(lv[i]-1)+j];
               } while(lv[++i]); 
               s=pp[0]*pp[0]; for(j=1;j<4;j++) s -=pp[j]*pp[j];
               if(key=='M') return sqrt(s);
               return s ;
            }
         case 'Y': 
               do for(j=0;j<4;j += 3) pp[j] += pvect[4*(lv[i]-1)+j]; while(lv[++i]); 
               return  log(( pp[0]+pp[3])/(pp[0]-pp[3]))/2; 
         case 'N':
               do for(j=0;j<4;j++) pp[j] += pvect[4*(lv[i]-1)+j]; while(lv[++i]);
               mp = sqrt(pp[1]*pp[1]+pp[2]*pp[2]+pp[3]*pp[3]);
               return  log(( mp+pp[3])/(mp-pp[3]))/2;
                                                                  
        case 'U': return usrfun(lv);
    }
    return 0;
}

int  checkPhysVal(char * name, char * key, char *plist)
{ int i=0,j=0;
  int n,k;

  
  while(name[i]==' '&&name[i]!=0) i++;
  *key=name[i++];
  
  if(*key==0) return 0;
  
   *key= toupper(*key);
  if(strchr("ACEJMPSTUYNZ",*key)==NULL) return 0; 
    
  if(*key == 'U')
  {  for( ;name[i] && name[i] != ' '  && i<6; i++) plist[j++]=name[i];
     plist[j]=0;
     for( ;name[i];i++) if(name[i] != ' ') return 0;
     return 1;
  }
    

  for( ; name[i]&&name[i]!=' '; i++)
  {  
    n=name[i]-'0';
    if(n<=0 || n>nin_int+nout_int) return 0;
    for(k=0; k<j;k++) if(plist[k] == n ) return 0;     
    plist[j++]=n;         
  }
  plist[j]=0;                                            
  for( ;name[i];i++) if(name[i] != ' ') return 0; 

  if(strchr("CAJ",*key)!=NULL && strlen(plist)!=2 )  return 0;
  if(strchr("P",*key)!=NULL && strlen(plist)<2 )  return 0;
  
  if(strchr("MS",*key)!=NULL && strlen(plist)<1)  return 0;
  
  if(strchr("JPTZ",*key)!=NULL)
           for(i=0;i<strlen(plist);i++) {if(plist[i]<=nin_int) return 0;} 
   
  if(strchr("MNYZ",*key) && !spole_(plist)) return 0;

  if(nin_int==1)
  { if( strchr("TYN",*key)) return 0;  
    if( strchr("ACP",*key) && (plist[0]==1 || plist[1]==1)) return 0;
  }    

  return 1;
}

#define MAXOUT 10

void cleanPVlist(physValRec * p)
{ physValRec *p1=p;
  while(p){ p1=p; p=p->next; free(p1);} 
}

int  checkPhysValN(char * name, char * key, physValRec **pLists)
{ int i=0,k,ln;
  char *chB;
  char pnum[10][10];
  int kk[10];
  *pLists=NULL;
  for(chB=name;*chB==' ';chB++) continue;
  
  *key=toupper(*chB); key[1]=key[2]=0;
  if(*key==0) return 0;


  if(!strchr("ACEJMPTUYNZ",*key)) return 0;
  chB++;  

  if(*key=='E' && nin_int==2 && (strcmp(chB,"1")==0 ||strcmp(chB,"2")==0)) 
  { *pLists=malloc(sizeof(physValRec)); 
    (*pLists)->next=NULL;
    (*pLists)->pstr[0]=chB[0]-'0'; (*pLists)->pstr[1]=0;
    return 1; 
  }
      
  if(*key=='M' && nin_int==2 && (strcmp(chB,"12")==0 ||strcmp(chB,"21")==0)) 
  { *pLists=malloc(sizeof(physValRec)); 
    (*pLists)->next=NULL;
    strcpy((*pLists)->pstr,"\1\2");
    printf("   => 12\n");
    return 1; 
  }    
   
  if(*key=='U')
  {
    if(strlen(chB)>9) return 0;
    *pLists=malloc(sizeof(physValRec));
    (*pLists)->next=NULL;
    strcpy((*pLists)->pstr,chB);
    return 1;
  }

  if(*chB=='^'||*chB=='_') { key[1]=*chB; chB++; } 
   
  if(*chB!='(') return 0;
  if(!strchr(chB,')')) return 0;

  for(ln=0; chB ;chB=strchr(chB,','),ln++)
  { char pname[100];
    if(ln==9) return 0;
    chB++;
    sscanf(chB,"%[^,)]",pname);
    trim(pname);
    pnum[ln][0]=0;
    k=0;
    if(strcmp(pname,"Jet"))
    {
      for(i=nin_int+1;i<=nin_int+nout_int;i++)
      if(strcmp(pname,pinf_int(Nsub,i,NULL,NULL))==0)
      { pnum[ln][k]=i;k++; pnum[ln][k]=0;}
      if(!k)
      { int n;
        for(n=1;n<=nprc_int;n++) 
        { for(i=nin_int+1;i<=nin_int+nout_int;i++)
             if(strcmp(pname,pinf_int(n,i,NULL,NULL))==0) break;
          if(i<=nin_int+nout_int) break;
        }
        if(n>nprc_int) return 0;
      }
    }else
    { long pNum;
      for(i=nin_int+1;i<=nin_int+nout_int;i++)
      { pinf_int(Nsub,i,NULL,&pNum);
        switch(labs(pNum))
        { case 1: case 2: case 3: case 4: case 5: case 21: case 81:case 83:
          pnum[ln][k]=i;k++; pnum[ln][k]=0;
        } 
      }  
      if(!k)
      { int n;
        for(n=1;n<=nprc_int;n++) 
        { for(i=nin_int+1;i<=nin_int+nout_int;i++)
          {  pinf_int(n,i,NULL,&pNum);
             pNum=labs(pNum);
             if((pNum>0 && pNum<6) ||pNum==81||pNum==83 || pNum==21 ) break;
          }   
          if(i<=nin_int+nout_int) break;
        } 
        if(n>nprc_int) return 0;
      }
    }  
  }
  if(ln<1) return 0;
  
  if(*key=='J' && ln!=2) return 0;
  if(*key=='P' && ln<2 ) return 0;
  if(strchr("CA",*key) &&(( ln==1 && nin_int==1)||ln>2)) return 0;   
  
  if(nin_int==1)
  { if( strchr("TYNZ",*key)) return 0;  
  }    
  
  for(i=0;i<ln;i++) if(pnum[i][0]==0) return 1;
  for(i=0;i<ln;i++)kk[i]=0;  
  for(;;)
  { char tmp[10];
    int i0;
    i=0;
    for(i=0;i<ln;i++) tmp[i]=pnum[i][kk[i]];
    tmp[i]=0;
    i0= (*key=='P')? 1:0;
    for(i=i0;i<ln-1;) if(tmp[i]<=tmp[i+1]) i++; else
    { int a=tmp[i+1]; tmp[i+1]=tmp[i]; tmp[i]=a;
      if(i==i0) i++;  else i--;
    } 
    for(i=i0;i<ln-1;i++) if(tmp[i]==tmp[i+1]) break;
    if(i==ln-1)
    { 
      if(*key=='P' && strchr(tmp+1,tmp[0]))continue;
      { physValRec * cc= *pLists;
        for(;cc;cc=cc->next) if(strcmp(tmp,cc->pstr)==0)break; 
        if(!cc)
        { physValRec * new=malloc(sizeof(physValRec));
          new->next=*pLists;
          strcpy(new->pstr,tmp);
          *pLists=new; 
        }
      }
    }
    for(i=ln-1;i>=0;i--) 
    {
     kk[i]++;
     if(pnum[i][kk[i]]==0) kk[i]=0; else break;
    } 
    if(i<0) break;
  }
  return 1;
}



