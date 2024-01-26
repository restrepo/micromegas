#ifdef PAPERS
2303.17021
1503.03513  DHS
1503.04935  ML
1407.6387   Lattice calculations
1309.5258
1103.2528
1007.2580   Lattice for Nf=3 and 4
0501232     HP
9510408     free energy

Thermodynamics

 e(T)= pi^2/30 T^4 gEff(T)  energy density    gEff(T)= Sum n_i*g1(M_i/T,eta_i)
 s(T)=2pi^2/45 T^3 hEff(T)  entropy density   hEff(T)= Sum n_i*h1(M_i/T,eta_i)

 de/dT = T ds/dT                              gEff(T)-hEff(T) =  dhEff(T)/dT /3 - dgEff(T)/dT/4
 p=s*T-e - pressure = free energy          =    M_PI*M_PI/90* (4*heff(T)-3*gEff(T) )
 I=e-3p  - trace of energy-momentum tensor =  2*M_PI*M_PI/15* (  gEff(T)-  hEff(T) )
 I/T^4 = T d(p(T)/T^4)/dT

#endif


#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"


static double Mc=0.85;

static  double hadrons(double T, double (* F)(double,int), int mode);  // performs summation  ndf(h)*F(M/T,eta)   over hadrons constructed from u,d,s  with M<2.5GeV
                                                                       // mode 0 - mesons, mode=1 mesons and baryons

static double mesons( double T, double (* F)(double,int)) { return hadrons(T, F, 0);}
static double allHadr( double T, double (* F)(double,int)){ return hadrons(T, F, 1);}

static double press1divT4(double t, int eta) { return M_PI*M_PI/90*(4*h1eff(t,eta)-3*g1eff(t,eta));}    // pressure*T^4  for 1 degree of freedom,  t=M/T
static double I1T4(double t,double eta){ return 2*M_PI*M_PI/15*(g1eff(t,eta)-h1eff(t,eta));}            // e-m trace *T^4 for 1 degree of freedom,  t=M/T


double static pdivT4(double T)  // Eq.16 1407.6387   - pressure obtained in  for Nf=3 in  lattice caculation.
                                //  Formula has a pole at T=45MeV  and below T=85MeV is replaced on summation over hadrons
{
  if(T>0.085)
  {
     double pid=95*M_PI*M_PI/180, Tc=0.154, t=T/Tc;
     double ct=3.8706, an=-8.7704, bn=3.9200, cn=0, dn=0.3419,
            t0=0.9761, ad=-1.2600, bd=0.8425, cd=0, dd=-0.0474;
     return  0.5*(1+tanh(ct*(t-t0))) * (pid +(an+(bn+(cn+ dn/t)/t)/t)/t ) / (1 +(ad+(bd+(cd+ dd/t)/t)/t)/t );
  } else return allHadr(T,press1divT4);

}


double static IdivT4(double T)     // (e-m trace)*T^4 obtained by  T* d pdivT4(T)/dT
{
  if(T>0.085)
  {
     double pid=95*M_PI*M_PI/180, Tc=0.154, t=T/Tc;
     double ct=3.8706, an=-8.7704, bn=3.9200, cn=0, dn=0.3419,
            t0=0.9761, ad=-1.2600, bd=0.8425, cd=0, dd=-0.0474;

     double A=(1+tanh(ct*(t-t0))),        B=(pid +(an+(bn+(cn+ dn/t)/t)/t)/t ),      C=(1 +(ad+(bd+(cd+ dd/t)/t)/t)/t );
     double A_=ct/pow(cosh(ct*(t-t0)),2), B_=  -(an+(2*bn+(3*cn+ 4*dn/t)/t)/t)/t/t,  C_=-(ad+(2*bd+(3*cd+ 4*dd/t)/t)/t)/t/t;
     return 0.5*(A_*B/C + A*B_/C - A*B*C_/C/C)*(T/Tc);
  } else return allHadr(T,I1T4);
}


static double hLat(double T)
{  return 45/(2*M_PI*M_PI)*(IdivT4(T)+4*pdivT4(T))+     // Nf=3 lattice calculation
  (
     2 + 6*h1eff(80/T,1) + 3*h1eff(91/T,1) + h1eff(125/T,1)        //  bosons
     + 7./8.*2*3                                                   //  neutrinos
     + 4*( h1eff(5E-4/T,-1) + h1eff(0.1/T,-1) + h1eff(1.7/T,-1))   //  charged leptons
     + 12*(h1eff(Mc/T,-1) +h1eff(5/T,-1) + h1eff(175/T,-1) )       //  heavy quarks
  );

}

static double cQuarkPressureDivT4(double T)   // Lattice calculation   1007.2580 Table 8 Q=mc/ms=11.85
{
  double  t[22]={ 147, 152, 158, 162, 166, 170, 175, 185, 200, 215, 228, 250, 275, 299, 330, 366, 400, 450, 500, 600, 800,1000};
  double  p[22]={0.01,0.02,0.04,0.05,0.06,0.07,0.08,0.10,0.14,0.18,0.21,0.26,0.32,0.37,0.43,0.49,0.55,0.63,0.70,0.81,0.96,1.05};

  if(T<0.147) return 0;
  if(T>1)
  {  double r=p[21]+ (p[21]-p[20])* (T-1)/( t[21]-t[20])*1000;  // linear interpolation abuve 1GeV
     double IG=7./8. * 12 * M_PI*M_PI/90;                       // Ideal Gas massles limit
     if(r> IG) return IG; else return r;
  }

  return polint3(T*1000,22,t,p);
}

static double cQuarkIGpressure(double T){ return 12*press1divT4(Mc/T,-1);}   // c-quark pressure


static double hEffI(double T) // ideal gas approximation
{

  double sum=  2+16+6*h1eff(80/T,1) + 3*h1eff(90/T,1) + h1eff(125/T,1)         // bosons
              +7./8.*2*3                                                       // neutrinos
              +7./8.*2*2*2 + 4*h1eff(1.7/T,-1)                                 // charged leptons
              +12*( 3*7./8.+ h1eff(1.4/T,-1) +h1eff(5/T,-1) + h1eff(175/T,-1)) // quarks
            ;
  return sum;
}


int main(int argc,char** argv)
{
  sortOddParticles(NULL);
  displayPlot("c-quark pressure","T", 0.1, 2, 0,2,"Latt",0,cQuarkPressureDivT4,NULL, "Ideal Gas",0,cQuarkIGpressure,NULL); // there is a good agreement for Mc=0.85.

  double Tmin=0.1,Tmax=1000;
  int scale=1;

  char * gh[3] = {"GG", "DHS", "LM"};
  double heff[3][200];
  for(int i=0;i<3;i++)
  {  char fname[100];
     sprintf(fname,"%s.thg",gh[i]);
     int err=loadHeffGeff(fname);
     if(err>0) printf("Table is   %s loaded. Number of records=%d \n", fname,err);
     else if(err==0) printf("Wrong path or file name\n");
     else  printf("wrong file format line %d\n", -err);
     for(int j=0;j<200;j++)
     { double T;
         if(scale) T=Tmin*pow(Tmax/Tmin,(j+0.5)/200);
         else      T=Tmin+  (Tmax-Tmin)*(j+0.5)/200;
         heff[i][j]=hEff(T);
     }
  }

  displayPlot("hEff", "T", Tmin,Tmax, scale,4, gh[0],200,heff[0],NULL
                                             , gh[1],200,heff[1],NULL
                                             , gh[2],200,heff[2],NULL
                                             , "hLat",0,hLat,NULL
                                                            );

}

static int C(int N )
{   int d[4];

    for(int i=0;i<4;i++)
    { int NN=N;
      N/=10;
      d[i]=NN-10*N;
    }
    if(d[3]!=0 || d[1]!=d[2] ) d[0]*=2;
    return d[0];
}

static  double hadrons(double T, double (* F)(double,int), int mode )
{  double sum=
//                        PDG            name
//LIGHT I= 1MESONS
  C(    111 )*F(140*1E-3/T,1)    //         π0
 +C(    211 )*F(140*1E-3/T,1)    //         π+
 +C(9000111 )*F(980*1E-3/T,1)    //0        a0
 +C(9000211 )*F(980*1E-3/T,1)    //+        a0
 +C(0100111 )*F(1300*1E-3/T,1)   //          π
 +C( 100211 )*F(1300*1E-3/T,1)   //+         π
 +C(  10111 )*F(1450*1E-3/T,1)   //0        a0
 +C(  10211 )*F(1450*1E-3/T,1)   //+        a0
 +C(9010111 )*F(1800*1E-3/T,1)   //0         π
 +C(9010211 )*F(1800*1E-3/T,1)   //+         π
 +C(    113 )*F(770*1E-3/T,1)    //0         ρ
 +C(    213 )*F(770*1E-3/T,1)    //+         ρ
 +C(  10113 )*F(1235*1E-3/T,1)   //0        b1
 +C(  10213 )*F(1235*1E-3/T,1)   //+        b1
 +C(  20113 )*F(1260*1E-3/T,1)   //0        a1
 +C(  20213 )*F(1260*1E-3/T,1)   //+        a1
 +C(9000113 )*F(1400*1E-3/T,1)   //0        π1
 +C(9000213 )*F(1400*1E-3/T,1)   //+        π1
 +C( 100113 )*F(1450*1E-3/T,1)   //0         ρ
 +C( 100213 )*F(1450*1E-3/T,1)   //+         ρ
 +C(9010113 )*F(1600*1E-3/T,1)   //0        π1
 +C(9010213 )*F(1600*1E-3/T,1)   //+        π1
 +C(9020113 )*F(1640*1E-3/T,1)   //0        a1
 +C(9020213 )*F(1640*1E-3/T,1)   //+        a1
 +C(  30113 )*F(1700*1E-3/T,1)   //0         ρ
 +C(  30213 )*F(1700*1E-3/T,1)   //+         ρ
 +C(9030113 )*F(1900*1E-3/T,1)   //0         ρ
 +C(9030213 )*F(1900*1E-3/T,1)   //+         ρ
 +C(9040113 )*F(2150*1E-3/T,1)   //0         ρ
 +C(9040213 )*F(2150*1E-3/T,1)   //+         ρ
 +C(    115 )*F(1320*1E-3/T,1)   //0        a2
 +C(    215 )*F(1320*1E-3/T,1)   //+        a2
 +C(  10115 )*F(1670*1E-3/T,1)   //0        π2
 +C(  10215 )*F(1670*1E-3/T,1)   //+        π2
 +C(9000115 )*F(1700*1E-3/T,1)   //0        a2
 +C(9000215 )*F(1700*1E-3/T,1)   //+        a2
 +C(9010115 )*F(2100*1E-3/T,1)   //0        π2
 +C(9010215 )*F(2100*1E-3/T,1)   //+        π2
 +C(    117 )*F(1690*1E-3/T,1)   //0        ρ3
 +C(    217 )*F(1690*1E-3/T,1)   //+        ρ3
 +C(9000117 )*F(1990*1E-3/T,1)   //0        ρ3
 +C(9000217 )*F(1990*1E-3/T,1)   //+        ρ3
 +C(9010117 )*F(2250*1E-3/T,1)   //0        ρ3
 +C(9010217 )*F(2250*1E-3/T,1)   //+        ρ3
 +C(    119 )*F(2040*1E-3/T,1)   //0        a4
 +C(    219 )*F(2040*1E-3/T,1)   //+        a4
//        LIGHT I= 0MESONS
 +C(    221 )*F(548*1E-3/T,1)    //          η
 +C(    331 )*F(958*1E-3/T,1)    //          η′
 +C(9000221 )*F(500*1E-3/T,1)    //         f0
 +C(9010221 )*F(980*1E-3/T,1)    //         f0
 +C( 100221 )*F(1295*1E-3/T,1)   //          η
 +C(  10221 )*F(1370*1E-3/T,1)   //         f0
 +C(9020221 )*F(1405*1E-3/T,1)   //          η
 +C( 100331 )*F(1475*1E-3/T,1)   //          η
 +C(9030221 )*F(1500*1E-3/T,1)   //         f0
 +C(  10331 )*F(1710*1E-3/T,1)   //         f0
 +C(9040221 )*F(1760*1E-3/T,1)   //          η
 +C(9050221 )*F(2020*1E-3/T,1)   //         f0
 +C(9060221 )*F(2100*1E-3/T,1)   //         f0
 +C(9070221 )*F(2200*1E-3/T,1)   //         f0
 +C(9080221 )*F(2225*1E-3/T,1)   //          η
 +C(    223 )*F(782*1E-3/T,1)    //          ω
 +C(    333 )*F(1020*1E-3/T,1)   //          φ
 +C(  10223 )*F(1170*1E-3/T,1)   //         h1
 +C(  20223 )*F(1285*1E-3/T,1)   //         f1
 +C(  10333 )*F(1380*1E-3/T,1)   //         h1
 +C(  20333 )*F(1420*1E-3/T,1)   //         f1
 +C( 100223 )*F(1420*1E-3/T,1)   //          ω
 +C(9000223 )*F(1510*1E-3/T,1)   //         f1
 +C(9010223 )*F(1595*1E-3/T,1)   //         h1
 +C(  30223 )*F(1650*1E-3/T,1)   //          ω
 +C( 100333 )*F(1680*1E-3/T,1)   //          φ
 +C(    225 )*F(1270*1E-3/T,1)   //         f2
 +C(9000225 )*F(1430*1E-3/T,1)   //         f2
 +C(    335 )*F(1525*1E-3/T,1)   //        f′2
 +C(9010225 )*F(1565*1E-3/T,1)   //         f2
 +C(9020225 )*F(1640*1E-3/T,1)   //         f2
 +C(  10225 )*F(1645*1E-3/T,1)   //         η2
 +C(9030225 )*F(1810*1E-3/T,1)   //         f2
 +C(  10335 )*F(1870*1E-3/T,1)   //         η2
 +C(9040225 )*F(1910*1E-3/T,1)   //         f2
 +C(9050225 )*F(1950*1E-3/T,1)   //         f2
 +C(9060225 )*F(2010*1E-3/T,1)   //         f2
 +C(9070225 )*F(2150*1E-3/T,1)   //         f2
 +C(9080225 )*F(2300*1E-3/T,1)   //         f2
 +C(9090225 )*F(2340*1E-3/T,1)   //         f2
 +C(    227 )*F(1670*1E-3/T,1)   //         ω3
 +C(    337 )*F(1850*1E-3/T,1)   //         φ3
 +C(    229 )*F(2050*1E-3/T,1)   //         f4
 +C(9000229 )*F(2220*1E-3/T,1)   //         fJ
 +C(9010229 )*F(2300*1E-3/T,1)   //         f4

//        STRANGE MESONS
//     13  +F( 498 *1E-3/T,1)   // 0     K0_L
//     31  +F( 498 *1E-3/T,1)   // 0     K0_S
 +C(    311 )*F(498*1E-3/T,1)    //        K0
 +C(    321 )*F( 494*1E-3/T,1)   //        K+
 +C(9000311 )*F(700*1E-3/T,1)    //0      K∗0
 +C(9000321 )*F(700*1E-3/T,1)    //+      K∗0
 +C(  10311 )*F(1430*1E-3/T,1)   //0      K∗0
 +C(  10321 )*F(1430*1E-3/T,1)   //+      K∗0
 +C( 100311 )*F(1460*1E-3/T,1)   //0        K
 +C( 100321 )*F(1460*1E-3/T,1)   //+        K
 +C(9010311 )*F(1830*1E-3/T,1)   //0        K
 +C(9010321 )*F(1830*1E-3/T,1)   //+        K
 +C(9020311 )*F(1950*1E-3/T,1)   //0      K∗0
 +C(9020321 )*F(1950*1E-3/T,1)   //+      K∗0
 +C(    313 )*F(892*1E-3/T,1)    //0       K∗
 +C(    323 )*F(892*1E-3/T,1)    //+       K∗
 +C(  10313 )*F(1270*1E-3/T,1)   //0       K1
 +C(  10323 )*F(1270*1E-3/T,1)   //+       K1
 +C(  20313 )*F(1400*1E-3/T,1)   //0       K1
 +C(  20323 )*F(1400*1E-3/T,1)   //+       K1
 +C( 100313 )*F(1410*1E-3/T,1)   //0       K∗
 +C( 100323 )*F(1410*1E-3/T,1)   //+       K∗
 +C(9000313 )*F(1650*1E-3/T,1)   //0       K1
 +C(9000323 )*F(1650*1E-3/T,1)   //+       K1
 +C(  30313 )*F(1680*1E-3/T,1)   //0       K∗
 +C(  30323 )*F(1680*1E-3/T,1)   //+       K∗
 +C(    315 )*F(1430*1E-3/T,1)   //0      K∗2
 +C(    325 )*F(1430*1E-3/T,1)   //+      K∗2
 +C(9000315 )*F(1580*1E-3/T,1)   //0       K2
 +C(9000325 )*F(1580*1E-3/T,1)   //+       K2
 +C(  10315 )*F(1770*1E-3/T,1)   //0       K2
 +C(  10325 )*F(1770*1E-3/T,1)   //+       K2
 +C(  20315 )*F(1820*1E-3/T,1)   //0       K2
 +C(  20325 )*F(1820*1E-3/T,1)   //+       K2
 +C(9010315 )*F(1980*1E-3/T,1)   //0      K∗2
 +C(9010325 )*F(1980*1E-3/T,1)   //+      K∗2
 +C(9020315 )*F(2250*1E-3/T,1)   //0       K2
 +C(9020325 )*F(2250*1E-3/T,1)   //+       K2
 +C(    317 )*F(1780*1E-3/T,1)   //0      K∗3
 +C(    327 )*F(1780*1E-3/T,1)   //+      K∗3
 +C(9010317 )*F(2320*1E-3/T,1)   //0       K3
 +C(9010327 )*F(2320*1E-3/T,1)   //+       K3
 +C(    319 )*F(2045*1E-3/T,1)   //0      K∗4
 +C(    329 )*F(2045*1E-3/T,1)   //+      K∗4
// +C(9000319 )*F(2500*1E-3/T,1)   //0       K4
// +C(9000329 )*F(2500*1E-3/T,1)   //+       K4
;

if(mode)
{ sum +=
// N-barions
 +C(2212 )*F(938 *1E-3/T,-1)   //p
 +C(2112 )*F(940 *1E-3/T,-1)   //n
 +2*(1+1)*F(1440 *1E-3/T,-1)  //N(1440) 1/2+
 +2*(3+1)*F(1520 *1E-3/T,-1)  //N(1520) 3/2−
 +2*(1+1)*F(1535 *1E-3/T,-1)  //N(1535) 1/2−
 +2*(1+1)*F(1650 *1E-3/T,-1)  //N(1650) 1/2−
 +2*(5+1)*F(1675 *1E-3/T,-1)  //N(1675) 5/2−
 +2*(5+1)*F(1680 *1E-3/T,-1)  //N(1680) 5/2+
 +2*(3+1)*F(1700 *1E-3/T,-1)  //N(1700) 3/2−
 +2*(1+1)*F(1710 *1E-3/T,-1)  //N(1710) 1/2+
 +2*(3+1)*F(1720 *1E-3/T,-1)  //N(1720) 3/2+
 +2*(5+1)*F(1860 *1E-3/T,-1)  //N(1860) 5/2+
 +2*(3+1)*F(1875 *1E-3/T,-1)  //N(1875) 3/2−
 +2*(1+1)*F(1880 *1E-3/T,-1)  //N(1880) 1/2+
 +2*(1+1)*F(1895 *1E-3/T,-1)  //N(1895) 1/2−
 +2*(3+1)*F(1900 *1E-3/T,-1)  //N(1900) 3/2+
 +2*(7+1)*F(1990 *1E-3/T,-1)  //N(1990) 7/2+
 +2*(5+1)*F(2000 *1E-3/T,-1)  //N(2000) 5/2+
 +2*(3+1)*F(2040 *1E-3/T,-1)  //N(2040) 3/2+
 +2*(5+1)*F(2060 *1E-3/T,-1)  //N(2060) 5/2−
 +2*(1+1)*F(2100 *1E-3/T,-1)  //N(2100) 1/2+
 +2*(3+1)*F(2120 *1E-3/T,-1)  //N(2120) 3/2−
 +2*(7+1)*F(2190 *1E-3/T,-1)  //N(2190) 7/2−
 +2*(9+1)*F(2220 *1E-3/T,-1)  //N(2220) 9/2+
 +2*(9+1)*F(2250 *1E-3/T,-1)  //N(2250) 9/2−
 +2*(1+1)*F(2300 *1E-3/T,-1)  //N(2300) 1/2+
// +2*(5+1)*F(2570 *1E-3/T,-1)  //N(2570) 5/2−
//+2*(11+1)*F(2600 *1E-3/T,-1) // N(2600) 11/2−
//+2*(13+1)*F(2700 *1E-3/T,-1) // N(2700) 13/2+

//  Delta-baryons
  +2*(3+1)*F(1232  *1E-3/T,-1) // Delta(1232) 3/2+
  +2*(3+1)*F(1600  *1E-3/T,-1) // Delta(1600) 3/2+
  +2*(1+1)*F(1620  *1E-3/T,-1) // Delta(1620) 1/2−
  +2*(3+1)*F(1700  *1E-3/T,-1) // Delta(1700) 3/2−
  +2*(1+1)*F(1750  *1E-3/T,-1) // Delta(1750) 1/2+
  +2*(1+1)*F(1900  *1E-3/T,-1) // Delta(1900) 1/2−
  +2*(5+1)*F(1905  *1E-3/T,-1) // Delta(1905) 5/2+
  +2*(1+1)*F(1910  *1E-3/T,-1) // Delta(1910) 1/2+
  +2*(3+1)*F(1920  *1E-3/T,-1) // Delta(1920) 3/2+
  +2*(5+1)*F(1930  *1E-3/T,-1) // Delta(1930) 5/2−
  +2*(3+1)*F(1940  *1E-3/T,-1) // Delta(1940) 3/2−
  +2*(7+1)*F(1950  *1E-3/T,-1) // Delta(1950) 7/2+
  +2*(5+1)*F(2000  *1E-3/T,-1) // Delta(2000) 5/2+
  +2*(1+1)*F(2150  *1E-3/T,-1) // Delta(2150) 1/2−
  +2*(7+1)*F(2200  *1E-3/T,-1) // Delta(2200) 7/2−
  +2*(9+1)*F(2300  *1E-3/T,-1) // Delta(2300) 9/2+
  +2*(5+1)*F(2350  *1E-3/T,-1) // Delta(2350) 5/2−
  +2*(7+1)*F(2390  *1E-3/T,-1) // Delta(2390) 7/2+
  +2*(9+1)*F(2400  *1E-3/T,-1) // Delta(2400) 9/2−
 +2*(11+1)*F(2420  *1E-3/T,-1) // Delta(2420) 11/2+
// +2*(13+1)*F(2750  *1E-3/T,-1) // Delta(2750) 13/2−
// +2*(15+1)*F(2950  *1E-3/T,-1) // Delta(2950) 15/2+

 +2*(1+1)*F(1116 *1E-3/T,-1) // Lambda(1116) 1/2+
 +2*(1+1)*F(1380 *1E-3/T,-1) // Lambda(1380) 1/2−
 +2*(1+1)*F(1405 *1E-3/T,-1) // Lambda(1405) 1/2−
 +2*(3+1)*F(1520 *1E-3/T,-1) // Lambda(1520) 3/2−
 +2*(1+1)*F(1600 *1E-3/T,-1) // Lambda(1600) 1/2+
 +2*(1+1)*F(1670 *1E-3/T,-1) // Lambda(1670) 1/2−
 +2*(3+1)*F(1690 *1E-3/T,-1) // Lambda(1690) 3/2−
 +2*(1+1)*F(1710 *1E-3/T,-1) // Lambda(1710) 1/2+
 +2*(1+1)*F(1800 *1E-3/T,-1) // Lambda(1800) 1/2−
 +2*(1+1)*F(1810 *1E-3/T,-1) // Lambda(1810) 1/2+
 +2*(5+1)*F(1820 *1E-3/T,-1) // Lambda(1820) 5/2+
 +2*(5+1)*F(1830 *1E-3/T,-1) // Lambda(1830) 5/2−
 +2*(3+1)*F(1890 *1E-3/T,-1) // Lambda(1890) 3/2+
 +2*(1+1)*F(2000 *1E-3/T,-1) // Lambda(2000) 1/2-
 +2*(3+1)*F(2050 *1E-3/T,-1) // Lambda(2050) 3/2−
 +2*(3+1)*F(2070 *1E-3/T,-1) // Lambda(2070) 3/2+
 +2*(5+1)*F(2080 *1E-3/T,-1) // Lambda(2080) 5/2−
 +2*(7+1)*F(2085 *1E-3/T,-1) // Lambda(2085) 7/2+
 +2*(7+1)*F(2100 *1E-3/T,-1) // Lambda(2100) 7/2−
 +2*(5+1)*F(2110 *1E-3/T,-1) // Lambda(2110) 5/2+
 +2*(3+1)*F(2325 *1E-3/T,-1) // Lambda(2325) 3/2−
 +2*(9+1)*F(2350 *1E-3/T,-1) // Lambda(2350) 9/2+
 +2*(1+1)*F(1189 *1E-3/T,-1) // Sigma+(1189) 1/2
 +2*(1+1)*F(1192 *1E-3/T,-1) // Sigma0(1192) 1/2
 +2*(1+1)*F(1197 *1E-3/T,-1) // Sigma−(1197) 1/2
 +2*(3+1)*F(1385 *1E-3/T,-1) // Sigma (1385) 3/2+
 +2*(3+1)*F(1580 *1E-3/T,-1) // Sigma (1580) 3/2−
 +2*(1+1)*F(1620 *1E-3/T,-1) // Sigma (1620) 1/2−
 +2*(1+1)*F(1660 *1E-3/T,-1) // Sigma (1660) 1/2+
 +2*(3+1)*F(1670 *1E-3/T,-1) // Sigma (1670) 3/2−
 +2*(1+1)*F(1750 *1E-3/T,-1) // Sigma (1750) 1/2−
 +2*(5+1)*F(1775 *1E-3/T,-1) // Sigma (1775) 5/2−
 +2*(3+1)*F(1780 *1E-3/T,-1) // Sigma (1780) 3/2+
 +2*(1+1)*F(1880 *1E-3/T,-1) // Sigma (1880) 1/2+
 +2*(1+1)*F(1900 *1E-3/T,-1) // Sigma (1900) 1/2−
 +2*(3+1)*F(1910 *1E-3/T,-1) // Sigma (1910) 3/2−
 +2*(5+1)*F(1915 *1E-3/T,-1) // Sigma (1915) 5/2+
 +2*(3+1)*F(1940 *1E-3/T,-1) // Sigma (1940) 3/2+
 +2*(3+1)*F(2010 *1E-3/T,-1) // Sigma (2010) 3/2−
 +2*(7+1)*F(2030 *1E-3/T,-1) // Sigma (2030) 7/2+
 +2*(5+1)*F(2070 *1E-3/T,-1) // Sigma (2070) 5/2+
 +2*(3+1)*F(2080 *1E-3/T,-1) // Sigma (2080) 3/2+
 +2*(7+1)*F(2100 *1E-3/T,-1) // Sigma (2100) 7/2−
 +2*(1+1)*F(2110 *1E-3/T,-1) // Sigma (2110) 1/2−
 +2*(3+1)*F(2230 *1E-3/T,-1) // Sigma (2230) 3/2+
 +2*(1+1)*F(1315 *1E-3/T,-1) // Xi0   (1315) 1/2
 +2*(1+1)*F(1323 *1E-3/T,-1) // Xi-   (1323) 1/2
 +2*(3+1)*F(1530 *1E-3/T,-1) // Xi    (1530) 3/2
 +2*(3+1)*F(1672 *1E-3/T,-1) // Omega-(1672) 3/2
;
}

  return sum;

}
