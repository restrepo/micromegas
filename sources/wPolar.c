#include "micromegas.h"
#include "micromegas_aux.h"

//xVtoxll(int Nin,int Nout,char**name,int *pdg, int lV, double *wV,  double *br)

int vPolar( char**N,  double*lng)
{ 
  char process[50];
  sprintf(process,"%s,%s->",N[0],N[1]);


//printf("Wpolar %s %s -> %s %s\n", N[0],N[1],N[2],N[3]);


  int v;   // position of vector particle W,or Z
  int pdg[4]; // of vector
    
  for(v=0;v<4;v++)  pdg[v]=pNum(N[v]); 
  for(v=2;v<4;v++)  if(pdg[v]==23 || abs(pdg[v])==24) break;
  
  if(v==4) return 1;
  double Mv=pMass(N[v]);

  int v_=5-v; // other outgoing particle
  
  double wV,br;
  
  numout*cc=xVtoxll(2,2,N,pdg, v, &wV,  &br);
  
  if(!cc) return 2;
  if(passParameters(cc)) return 3;

  REAL m[5],pvect[20],pcm1,pcm2,chY,shY;
  char *p[5];
  for(int i=0;i<5;i++) { p[i]=cc->interface->pinf(1,i+1,m+i,NULL); m[i]=Fabs(m[i]);}

  int i3,ie,in;
  for(int i=2;i<5;i++)    if(strcmp(p[i],N[v_])==0)  { i3=i; break;}
  for(int i=2;i<5;i++)    if(i!=i3)   {ie=i; break;} 
  for(int i=ie+1;i<5;i++) if(i!=i3)   { in=i; break;} 
//printf("proc=%s %s -> %s %s %s ; i3=%d ie=%d in=%d\n",p[0],p[1],p[2],p[3],p[4],i3,ie,in);   
  
  if( m[0]+m[1] < m[i3]+ Mv) return 4;

  pcm1=decayPcm(m[0]+m[1],m[i3],Mv);
  pcm2=decayPcm(Mv,m[ie],m[in]);

  for(int i=0;i<20;i++) pvect[i]=0;
  pvect[0]=m[0];
  pvect[4]=m[1];
  pvect[4*i3]= Sqrt(m[i3]*m[i3]+pcm1*pcm1);
  pvect[4*i3+3] = -pcm1;  

  double r[3];
  for(int i=0;i<3;i++)
  {  REAL csfi=i-1;  

     pvect[4*ie]=Sqrt(m[ie]*m[ie]+pcm2*pcm2);
     pvect[4*ie+3]=pcm2*csfi;
     pvect[4*ie+2]=pcm2*Sqrt(1-csfi*csfi);
  
     pvect[4*in]=Sqrt(m[in]*m[in]+pcm2*pcm2);
     pvect[4*in+3]=-pcm2*csfi;
     pvect[4*in+2]=-pcm2*Sqrt(1-csfi*csfi);  

     chY=Sqrt(1+pcm1*pcm1/m[i3]/m[i3]);
     shY=Sqrt(pcm1*pcm1/m[i3]/m[i3]);  
  
     { REAL p0=pvect[4*ie], p3=pvect[4*ie+3];
       pvect[4*ie]=  chY*p0 + shY*p3;
       pvect[4*ie+3]=shY*p0 + chY*p3;

       p0=pvect[4*in]; p3=pvect[4*in+3];
       pvect[4*in]=  chY*p0 + shY*p3;
       pvect[4*in+3]=shY*p0 + chY*p3;
     }
     int err_code=0;
     r[i]=(cc->interface->sqme)(1,1./*GG*/,pvect,NULL,&err_code);
     if(err_code) return 5;
  }
//printf("r[0]=%e r[1]=%e r[3]=%e\n", r[0],r[1],r[2]);  
  r[2]/=4;
  r[0]/=4;
  r[1]=(r[1]-r[0]-r[2])/2;
  double  s=r[0]+r[1]+r[2];                 
  *lng=r[1]/s;  //    *left=r[2]/s;  *right=r[0]/s;
  if(*lng<0) *lng=0;
// printf("       lng=%E\n",*lng); 
  return 0;
}  
