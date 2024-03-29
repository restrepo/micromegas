/*
  Author Alexander Pukhov.
*/

#include "include/pdt.h"
#include<math.h>
/* STEQ INTERPOLATION CODES {*/

/*extern double qql[],xxl[],f3[],f7[],*cc3,*ccb,qqlb[],fb[];*/

static double  polint(double x, int n,  double *xa, double *ya)
{  double z[10];
   int i,m;

   for(i=0;i<n;i++) z[i]=ya[i];

   for(m=1;m<n;m++) for(i=0;i<n-m;i++)
   z[i]=(z[i]*(xa[i+m]-x) - z[i+1]*(xa[i]-x))/(xa[i+m]-xa[i]);
   return z[0];
}    


static int  leftX(int dim, double * xa, double x)
{  int k1,k2,k3;

   if(x<xa[0]) return 0;
   if(x>xa[dim-1]) return dim-3;

   k1=0; 
   k2=dim-1;
      
   while(k2-k1>1)
   { k3=(k1+k2)/2;
     if(xa[k3]>x)k2=k3; else k1=k3;
   } 

   k3=k1;
   if(k3<0) k3=0;
   if(k3>dim-2) k3=dim-2; 
   return k3;
}

static  double cubicInt1D(double x,  double * X, double *F)
{

   if(x<X[0]&& F[0]>0 && F[1]>0)  return pow(F[0], (x-X[1])/(X[0]-X[1]))*pow(F[1], (X[0]-x)/(X[0]-X[1]));
   if(x>X[3]&& F[2]>0 && F[3]>0)  return pow(F[2], (x-X[3])/(X[2]-X[3]))*pow(F[3], (X[2]-x)/(X[2]-X[3]));
//   return polint(x,4,X,F);
   
   if(x<=X[1]) return polint(x,3,X,F);
   if(x>=X[2]) return polint(x,3,X+1,F+1);

   double d0=(F[1]-F[0])/(X[1]-X[0]),d1=(F[2]-F[1])/(X[2]-X[1]),d2=(F[3]-F[2])/(X[3]-X[2]);

   double dL=d0*(X[2]-X[1])/(X[2]-X[0])+d1*(X[1]-X[0])/(X[2]-X[0]);
   double dR=d1*(X[3]-X[2])/(X[3]-X[1])+d2*(X[2]-X[1])/(X[3]-X[1]);
   double dx=X[2]-X[1];

   double t= (x-X[1])/dx, t2=t*t,t3=t2*t;
   dL*=dx; dR*=dx;

//printf(" F1=%e F2=%e dL=%E dR=%E\n", F[1],F[2],dL, dR); 
   return F[1]*(1-3*t2+2*t3) + F[2]*t2*(3-2*t) +dL*t*(1-t)*(1-t)+ dR*t2*(t-1);
}


static  double biCubicLogXQ(double x,double q,pdtStr * W)
{  
   double logQ=log(q);
   double logX=log(x);
   double f[4];
   int i; 
   int px = leftX(W->nx, W->x_grid, x);   if(px>0) px--; if(px>W->nx-4) px=W->nx-4; 
   int pq = leftX(W->nq, W->q_grid, q);   if(pq>0) pq--; if(pq>W->nq-4) pq=W->nq-4; 
  
   for(i=0;i<4;i++) f[i]= cubicInt1D(logX, W->lx_grid +px, W->strfun +W->nx*(pq+i)+px );
   if(q>1) return cubicInt1D(logQ,W->lq_grid+pq,f); 
   else  
   {  for(i=0;i<4;i++) f[i]/=W->q_grid[pq+i]*W->q_grid[pq+i];   
      return cubicInt1D(logQ,W->lq_grid+pq,f)*q*q;
   }   
}


static double qSplineCteq6(double x,double *xa, double *ya)
{
  double s12=xa[0]-xa[1], s13=xa[0]-xa[2], s23=xa[1]-xa[2], s24=xa[1]-xa[3], 
         s34=xa[2]-xa[3];
  double sy2=x-xa[1], sy3=x-xa[2];

  double c1=s13/s23, c2=s12/s23, c3=s34/s23, c4=s24/s23;

  double cxx=sy2*sy3/( s12*s34 - (s12+s13)*(s24+s34));
  double c5=(s34*sy2-(s24+s34)*sy3)*cxx/s12, c6=((s12+s13)*sy2-s12*sy3)*cxx/s34;
  double res =

        ( c5*(ya[0] -  ya[1]*c1 +ya[2]*c2 )
         +c6*(ya[3]+ya[1]*c3 - ya[2]*c4) 
         +ya[1]*sy3 - ya[2]*sy2
        )/s23;
  return res;
}


double int_cteq6(double x, double q, pdtStr * W)
{
  double x3=pow(x,0.3);
  double loglogQ=log(log(q/0.22));
  double * q_grid=W->q_grid_aux;
  double * x_grid=W->x_grid_aux;
  int pq = leftX(W->nq, q_grid, loglogQ)-1;
  int px = leftX(W->nx, W->x_grid_aux, x3)-1;
  int i;
  int qExt=0;
  int xExt=0;
  double tmp[4];
  
  if(pq<0)  { pq=0; qExt=-1;} else if(pq > W->nq-4) { pq=W->nq-4; qExt=1;}
  if(px<=0) { px=0; xExt=-1;} else if(px > W->nx-4) { px=W->nx-4; xExt=1;}

       if(xExt==0) for(i=0;i<4;i++) 
            tmp[i]=qSplineCteq6(x3, x_grid+px, W->strfun+W->nx*(pq+i)+px);
  else if(xExt>0) for(i=0;i<4;i++) 
            tmp[i]=polint(x3,4, x_grid+px, W->strfun+W->nx*(pq+i)+px); 
  else  
  {   double ftmp[4];
      int j;
      ftmp[0]=0;
      for(i=0;i<4;i++)
      {  for(j=1;j<4;j++) ftmp[j]= pow(x_grid[j],2/0.3)*W->strfun[W->nx*(pq+i)+j];
         tmp[i]= polint(x3,4,x_grid,ftmp)/(x*x);
      }         
  }
             
  if(qExt) return polint(loglogQ,4,q_grid+pq,tmp); 
      else return qSplineCteq6(loglogQ,q_grid+pq,tmp);
}


double interFunc(double x, double q, pdtStr * W)
{ 

  if(W->q_grid)
  { if(q<=W->q_threshold) return 0; 
    if(q< W->q_min)
    { double q0=W->q_min;
      double q1=q0*1.005;
//             q1=W->q_grid[1];
      double f0=(W->interpolation)(x,q0,W); 
      double f1=(W->interpolation)(x,q1,W); 
      double k1=q1/q0, k2=k1*k1;
      double C=1+log(f1/f0/k2)/k2/log(k2); 
      double k=q/q0;
      if(W->index){ f0/=x;f1/=x;}
//printf("C=%E\n",C);
      if(C<-2.5) C=-2.5;            
      if(f0<=0) return 0;
//      double C=(f1-f0)/f0/(0.1*q0);  
//printf("f0=%E f1=%E C=%e f1'=%E  f(%e)=%e\n",f0,f1,C,f0*pow(k1,2+k1*(C-2)),k, f0*pow(k,2+ k*k*(C-2)));          
      return f0*pow(k*k,1+ k*k*(C-1)); 
    }
    if(q>W->q_max) 
    { W->nLargeQ++; 
      if(W->nLargeQ==1) printf("Too large scale Q=%.2E < Qmax=%.2E for  structure function  \"PDT:%s\" \n",q,W->q_max, W->source); 
       else if(W->nLargeQ==101) printf(" More then 100 calls of structure function \"PDT:%s\" with large Q \n",W->source);  
       q=W->q_max;  
    } else if(q<W->q_min) W->nSmallQ++;  
  }
  if(x<W->x_min)  W->nSmallX++;

  double res=(W->interpolation)(x,q,W);
//  printf("res=%e\n",res);
  if(W->index) res/=x;  
  return res; 
}


double interAlpha(double q, pdtStr * W )
{ 
  double logQ;
  int pq;
  if(!W->alpha) return -1.;
  logQ=log(q);
  pq=leftX(W->nq,W->lq_grid,logQ); 
  if(pq>W->nq-3) pq = W->nq-3;
  return polint(logQ, 3, W->lq_grid+pq, W->alpha+pq);
}


void freePdtData( pdtStr * data)
{ if(!data)return;
  if(data->source)     {free(data->source);      data->source=NULL;     } 
  if(data->x_grid)     {free(data->x_grid);      data->x_grid=NULL;     }
  if(data->q_grid)     {free(data->q_grid);      data->q_grid=NULL;     }
 
  if(data->lx_grid)     {free(data->lx_grid);      data->lx_grid=NULL;     }
  if(data->lq_grid)     {free(data->lq_grid);      data->lq_grid=NULL;     }
       
        
  
  if(data->alpha)      {free(data->alpha);       data->alpha =NULL;     }
  if(data->strfun)     {free(data->strfun);      data->strfun=NULL;     }
  if(data->q_grid_aux) {free(data->q_grid_aux);  data->q_grid_aux=NULL; }
  if(data->x_grid_aux) {free(data->x_grid_aux);  data->x_grid_aux=NULL; }
//  if(data->aux)        {free(data->aux);         data->aux=NULL;        }
}


int getPdtData(char * file, int n_parton, pdtStr * data )
{ char pattern[20];
  char buff[100]; 
  char c;
  int  nx=0,nq=1;  
  FILE *f=fopen(file,"r"); 
  int errNo=0;
  
  if(!f) return -1;
//printf("getPdtData (%d)   %s\n",n_parton,file);
  data->index=0;
  data->set=0;
  data->nq=0;
  data->nx=0;
  data->source=NULL;
  data->x_grid=NULL;
  data->q_grid=NULL;
  data->lx_grid=NULL;
  data->lq_grid=NULL; 
  data->x_grid_aux=NULL;
  data->q_grid_aux=NULL;
  data->alpha =NULL;
  data->strfun=NULL;
  data->interpolation=&biCubicLogXQ;
  data->q_threshold=0;
  data->mass=1;
  data->qt0=0;
 
  data->pow0=0;
  data->pow1=0;
  data->x_min=0;

  data->nSmallX=0;
  data->nSmallQ=0;
  data->nLargeQ=0;
  data->nLargeX=0;
  
  sprintf(pattern,"%d-parton",n_parton);

  while(1==fscanf(f,"%c",&c))   if(c=='#')
  { double qq; int i;
    
    fscanf(f,"%s",buff);
    if(!strcmp(buff,"distribution") && !data->source )
    {
      fscanf(f," \"%[^\"]",buff);
      char* ch=strchr(buff,'(');
      if(ch) ch[0]=0;
      data->source=malloc(1+strlen(buff));
      strcpy(data->source, buff);      
    }
    else if(!strcmp(buff,"Index"))
    {  if(1!=fscanf(f,"%d",&data->index)) goto errexit;}
    else if(!strcmp(buff,"Set"))
    { if(1!=fscanf(f,"%d",&data->set)) goto errexit;} 
    else if(!strcmp(buff,"Mass"))
    {  if(1!=fscanf(f,"%lf",&data->mass)) goto errexit;}  
    else if(!strcmp(buff,"Q_grid"))
    {  long fpos=ftell(f);
       if(data->q_grid || data->strfun) goto errexit;
       nq=0;
       while(fscanf(f,"%lf",&qq)) nq++;
       if(nq<3)  goto errexit; 
       data->nq=nq;
       data->q_grid=malloc(nq*sizeof(double));
       fseek(f,fpos,SEEK_SET);
       for(i=0;i<nq;i++) 
       { if(fscanf(f,"%lf", data->q_grid+i)!=1) goto errexit;
         if(i){ if(data->q_grid[i-1]>=data->q_grid[i]) goto errexit;}
         else if(data->q_grid[0]<=0)  goto errexit;
       }
    } 
    else if(!strcmp(buff,"Interpolation"))
    {                  
       fscanf(f,"%s",buff);
       if(0==strcmp(buff,"biCubicLogXQ")) data->interpolation=&biCubicLogXQ;
       else if(0==strcmp(buff,"CTEQ6")) {  data->interpolation= &int_cteq6;}      
       else printf("Inknown interpolation '%s', biCubicLogXQ  is used instead\n",buff); 
    }     
    else if(!strcmp(buff,"X_grid"))
    {  long fpos=ftell(f);
       if(data->x_grid)  goto errexit; 
       while(fscanf(f,"%lf",&qq)) nx++; 
       data->nx=nx;
       if(nx<3)  goto errexit; 
       data->x_grid=malloc(nx*sizeof(double));
       fseek(f,fpos,SEEK_SET);
       for(i=0;i<nx;i++) fscanf(f,"%lf", data->x_grid+i);
       for(i=1;i<nx;i++) if(data->x_grid[i-1]>=data->x_grid[i]) goto errexit;   
    }
    else if(!strcmp(buff,"Alpha"))    
    {  if(!data->q_grid)  goto errexit; 
       data->alpha=malloc(nq*sizeof(double));
       for(i=0;i<nq;i++)  
       if(fscanf(f,"%lf", data->alpha+i)!=1)  goto errexit; 
       if(fscanf(f,"%lf", &qq)==1)  goto errexit; 
    } else if(!strcmp(buff,pattern))
    {  int nn=nq*nx;
       data->strfun=malloc(nn*sizeof(double));
       for(i=0;i<nn;i++) 
       if(fscanf(f,"%lf", data->strfun+i)!=1) goto errexit; 
       if(fscanf(f,"%lf", &qq)==1)  goto errexit;
//{ int i; for(i=0;i<10;i++) printf(" %E ", data->strfun[i]);
//  printf("\n");
// }       
       break;
    } else if(!strcmp(buff,"q_threshold"))
         { if(fscanf(f,"%lf", &(data->q_threshold))!=1)  goto errexit; } 
      else if(!strcmp(buff,"x_min"))
         { if(fscanf(f,"%lf", &(data->x_min))!=1)  goto errexit; }
  }
  
//printf("q_grid=%p x_grid-%p\n", data->q_grid, data->x_grid);
  
  if(data->strfun==NULL) {errNo=-2; goto errexit;}  
//printf("ok11\n"); 

  if(data->x_min==0.) data->x_min=data->x_grid[0];
      
  if(data->q_grid)
  { int i;
    data->q_min=data->q_grid[0];
    data->q_max=data->q_grid[data->nq-1];

    for(i=0;i<nq;i++) if(data->q_grid[i]<data->q_threshold) data->qt0=i+1;
 
    data->lq_grid=(double*)malloc(sizeof(double)*nq);

    for(i=0;i<nq;i++)data->lq_grid[i]=log(data->q_grid[i]);
  }
  
  if(data->x_grid)
  { int i;
    data->lx_grid=(double*)malloc(sizeof(double)*nx);
    for(i=0;i<nx;i++)data->lx_grid[i]=log(data->x_grid[i]);
  }
   
  if(data->interpolation == &int_cteq6) 
  { int i; 
    if(!data->q_grid) { errNo=-3; goto errexit;}

    data->x_grid_aux=(double*)malloc(sizeof(double)*nx);

    if(data->x_grid[0]==0) data->x_grid_aux[0]=0; else data->x_grid_aux[0]=pow( data->x_grid[0], 0.3);
    for(i=1;i<nx;i++) data->x_grid_aux[i]=pow( data->x_grid[i], 0.3);

    data->q_grid_aux=(double*)malloc(sizeof(double)*nq);
    for(i=0;i<nq;i++) data->q_grid_aux[i]=log( data->lq_grid[i]-log(0.22)); 
  }

  fclose(f); 
  return 0;
  errexit:
  { if(errNo==0) errNo=ftell(f); 
    fclose(f);   
    freePdtData(data);
    printf("error Exit %d\n",errNo);
    return errNo;
  } 
}


void delPdtList(pdtList * list)
{ 
  while(list)
  { pdtList * next=list->next;;
    free(list->file);
    free(list->name);
    free(list->partons);
    free(list->items);
    free(list);
    list=next; 
  }
}


long  makePdtList(char * file,  pdtList ** list)
{ char s[100];
  char dName[100];
  FILE * f=fopen(file,"r");  
  int partons_[100], positions_[100];
  int N,L;

  if(!f) return 0;
  while(fscanf(f,"%s",s) ==1 && strcmp(s,"#distribution")!=0) ;
  for(L=0;!feof(f);)
  { long mother;
    int bcount=0,pos;
    char ch;
    L=0;
    if(2!=fscanf(f," \"%[^\"]%*c %ld => ",dName,&mother)) break; 
    for(pos=1,ch=0;ch!='#';)
    {  if(1==fscanf(f," %d ",&N))
       {
         partons_[L]=N;   positions_[L]=pos; L++;       
         if(bcount==0) pos++; 
       }
       else 
       {  if(1!=fscanf(f,"%c",&ch)||!strchr("(#)",ch)
          ||(bcount &&strchr("(#",ch)) || (!bcount && ch==')') )goto exi;
          if(ch=='(') bcount=1; 
          else if(ch==')') { bcount=0; pos++;}
       }      
    } 
    {  pdtList * new=malloc(sizeof(pdtList));
       new->name=malloc(strlen(dName)+1);
       strcpy(new->name,dName);
       new->file=malloc(strlen(file)+1);
       strcpy(new->file,file);
       new->partons=malloc(sizeof(int)*(L+1));
       memcpy(new->partons,partons_,sizeof(int)*L);
       new->partons[L]=0;

       new->items=malloc(sizeof(int)*(L+1));
       memcpy(new->items,positions_,sizeof(int)*L);
       new->items[L]=0;
       new->beamP=mother;
       new->next=*list;
       *list=new;   
    }
    if(1!=fscanf(f,"%s",s)) break;
    if(strcmp(s,"distribution")==0) continue; else { fclose(f); return 0;} 
  }
exi:
  N=ftell(f);
  fclose(f); return N; 
}

int checkPartons( int * pNum, pdtList * L)
{ int *p;
  for(;*pNum; pNum++)
  { 
    if(*pNum==81 || *pNum==83) 
    {  for(p=L->partons ;*p;p++) if(*p==3) break; if(*p==0) return 0;
       for(p=L->partons ;*p;p++) if(*p==1) break; if(*p==0) return 0;
    } else if(*pNum==-81 || *pNum==-83)     
    {  for(p=L->partons ;*p;p++) if(*p==-3) break; if(*p==0) return 0;
       for(p=L->partons ;*p;p++) if(*p==-1) break; if(*p==0) return 0;
    }  else  
    {
      for(p=L->partons ;*p;p++) if(*pNum==*p) break;
      if(*p==0) return 0;
    }  
  }
  return 1;
}

