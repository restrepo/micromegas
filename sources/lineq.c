#include "../include/micromegas_aux.h"
int solveLinEq( int N,double *A, double *c)
{
  int *index=(int*) malloc(N*sizeof(int));
  int i,j,k;
  for(i=0;i<N;i++) index[i]=i;
  
  for(k=0;k<N;k++)
  {  double max=fabs(A[index[k]*N+k]);
     int i0=k;

     for(i=k+1;i<N;i++) 
     {  double  m=fabs(A[index[i]*N+k]);
        if(m>max) { i0=i; max=m;}
     } 
     if(max==0){ free(index); return k+1;}
     if(i0!=k) { int mem=index[k]; index[k]=index[i0];index[i0]=mem;}
     
     for(i=k+1;i<N;i++) 
     { double b=A[index[i]*N+k]/A[index[k]*N+k];
       for(j=k;j<N;j++) A[index[i]*N+j]-=b*A[index[k]*N+j];
       c[index[i]]-=b*c[index[k]];
     }
  }

  for(k=N-1;k>=0;k--)
  {  
    for(j=N-1;j>k;j--) c[index[k]]-=A[index[k]*N+j]*c[index[j]];
    c[index[k]]/=A[index[k]*N+k];
  }
                               
  for(k=0;k<N-1;k++) if(k!=index[k])
  {  
    i=index[k];
    double ck=c[i];
    c[i]=c[k];
    c[k]=ck;  
    for(j=k+1;index[j]!=k;j++) continue;
    index[j]=i;
    index[k]=k;
  }
                             
  free(index);
  return 0; 
}
