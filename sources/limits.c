
#include"micromegas.h"
#include"micromegas_aux.h"

int Zinvisible(void)
{
  char* zname;
  txtList L;
  int i;
  double width,br;
  
  zname=pdg2name(23);
  if(!zname) { printf("Z boson  is absent\n"); return 0;}   
  width=pWidth(zname,&L);
  for(i=0,br=0;i<Ncdm;i++) if(CDM[i] && 2*pMass(CDM[i])< pMass(zname))
  { char channel[40];
    sprintf(channel,"%s %s",CDM[i],antiParticle(CDM[i]));
    br+=findBr(L,channel);   
  }
  
//  printf(" partial width Z->DM=%.1EMeV\n", br*1000*width);  
  if(br*1000*width>0.5)
  { printf(" partial width Z->DM,DM is %.2EMeV,  more than 0.5 MeV. See 1401.2447\n",br*1000*width);
    return 1;
  } else return 0;
}

extern int zinvisible_(void); //Fortran 
int zinvisible_(void) { return Zinvisible();}

int LspNlsp_LEP(double *cs_out)
{  
  char *e=pdg2name(11),*e_=pdg2name(-11),*z_=pdg2name(23);
  if(!e || !e_ || z_) { if(cs_out) *cs_out=0;  return 0;} 
  int VZ_tmp=VZdecay;
  if(!VZ_tmp) { VZdecay=1; cleanDecayTable();}

  double res=0;
  double P=104;

  for(int  n=1;n<=Ncdm;n++)  if(CDM[n] && pMass(CDM[n])<P/2.)   // cycle over  DM sectors
  {   
     char process[100];
     sprintf(process,"%s,%s->%s,1*x",e,e_,CDM[n]);
     numout*cc= newProcess(process);
     if(!cc) continue;
     for(int k=1;k<=cc->interface->nprc;k++)   // cycle over second second odd particles
     {  
         int l,pdg_;
         char*cdm_;
         REAL m_;
         for(l=3;l<=4;l++)
         { cdm_=cc->interface->pinf(k,l,&m_,&pdg_);
           if(strcmp(cdm_, CDM[k])) break;
         }  
         if(l>4) continue;
         if(m_+pMass(CDM[n])>2*P) continue;
         if(abs(pdg_) == abs(pNum(CDM[n]))) continue;
         double dM=pMass(cdm_)-pMass(CDM[n]);
         if(dM< 2 ) continue;  // for smaller dM we have to treat Z->p0 by some special way
    
         txtList L;
         double width=pWidth(cdm_,&L);
         double brZ=findBr(L,z_);
         double brZu=brZ*0.119, brZd=brZ*0.152;
         brZ=brZu+2*brZd; // d,u,s
         if(dM< 3)  brZ +=brZu; // c
         if(dM< 10) brZ +=brZd; // b
         int err;
         double cs=cs22(cc,k,P,-1,1 ,&err);
         if(m_>100) res+=cs*brZ/0.1; else res+=cs*brZ/0.5;
      }   
  }      
  if(VZ_tmp!=VZdecay) { VZdecay=VZ_tmp; cleanDecayTable();}
  if(cs_out) *cs_out=res;
  if( res>1) return 1; else return 0;    
}

int  lspnlsp_lep_(double *cs_out) { return  LspNlsp_LEP(cs_out);}