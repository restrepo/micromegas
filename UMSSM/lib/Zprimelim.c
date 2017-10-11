#include<math.h>
#include"../../sources/micromegas.h"
#include "pmodel.h"

// Routines to derive ATLAS Z' limits with 20 fb-1 at sqrt(s) = 8 TeV for the UMSSM and extrapolation for CMS Z' limits

// For a fixed number of tE6 (36 exactly), get the cubic interpolation of sigma*Br(Z2 -> l+l-)_th (in fb) for a Z2 mass between 2 and 2.7 TeV.
// Case of some sigma*Br(Z2 -> l+l-)_th obtained for aZZ different from 0 commented.

double sBthtE6Fixed(double xx, int jj)
{
if (xx>2700 || xx<2000) {printf("**** ERROR **** First argument outside range [2000,2700] GeV. Please try a value in this range.\n");exit(1);}
double Zpmass[5]={2000,2200,2400,2600,2700};

if (jj==1)       {double sigB[5]={0.71815,0.36425,0.1823,0.09406,0.08}; return polint4(xx,5,Zpmass,sigB);} //te6 == 1.5708
//else if (jj==2)  {double sigB[5]={0.71155,0.35745,0.1837,0.094955,0.08};return polint4(xx,5,Zpmass,sigB);} //azz == -0.001
//else if (jj==3)  {double sigB[5]={0.72375,0.36195,0.1824,0.095945,0.08};return polint4(xx,5,Zpmass,sigB);} //azz ==  0.001
else if (jj==2)  {double sigB[5]={0.7401,0.3650,0.1827,0.09556,0.08};   return polint4(xx,5,Zpmass,sigB);} //te6 == 1.4708
//else if (jj==5)  {double sigB[5]={0.7439,0.3737,0.1851,0.09715,0.08};   return polint4(xx,5,Zpmass,sigB);} //azz == -0.001
//else if (jj==6)  {double sigB[5]={0.7350,0.3659,0.1874,0.09713,0.08};   return polint4(xx,5,Zpmass,sigB);} //azz ==  0.001
else if (jj==3)  {double sigB[5]={0.7566,0.3652,0.1858,0.09609,0.08};   return polint4(xx,5,Zpmass,sigB);} //te6 == 1.3708
//else if (jj==8)  {double sigB[5]={0.7406,0.3686,0.1861,0.09676,0.08};   return polint4(xx,5,Zpmass,sigB);} //azz == -0.001
//else if (jj==9)  {double sigB[5]={0.7585,0.3733,0.1885,0.09781,0.08};   return polint4(xx,5,Zpmass,sigB);} //azz ==  0.001
else if (jj==4) {double sigB[5]={0.7611,0.3731,0.1866,0.09486,0.08};   return polint4(xx,5,Zpmass,sigB);} //te6 == 1.3181
//else if (jj==11) {double sigB[5]={0.7570,0.3697,0.1858,0.09501,0.08};   return polint4(xx,5,Zpmass,sigB);} //azz == -0.001
//else if (jj==12) {double sigB[5]={0.7673,0.3768,0.1870,0.0970,0.08};    return polint4(xx,5,Zpmass,sigB);} //azz ==  0.001
else if (jj==5) {double sigB[5]={0.7613,0.3749,0.1862,0.09578,0.08};   return polint4(xx,5,Zpmass,sigB);} //te6 == 1.2708
//else if (jj==14) {double sigB[5]={0.7674,0.3717,0.1876,0.09499,0.08};   return polint4(xx,5,Zpmass,sigB);} //azz == -0.001
//else if (jj==15) {double sigB[5]={0.7641,0.3682,0.1835,0.09556,0.08};   return polint4(xx,5,Zpmass,sigB);} //azz ==  0.001
else if (jj==6) {double sigB[5]={0.7640,0.3707,0.1867,0.09403,0.075};  return polint4(xx,5,Zpmass,sigB);} //te6 == 1.1708
//else if (jj==17) {double sigB[5]={0.7625,0.3738,0.1874,0.09367,0.08};   return polint4(xx,5,Zpmass,sigB);} //azz == -0.001
//else if (jj==18) {double sigB[5]={0.7725,0.3718,0.1840,0.09531,0.08};   return polint4(xx,5,Zpmass,sigB);} //azz ==  0.001
else if (jj==7) {double sigB[5]={0.7809,0.3720,0.1854,0.09417,0.075};  return polint4(xx,5,Zpmass,sigB);} //te6 == 1.0708
//else if (jj==20) {double sigB[5]={0.7708,0.3734,0.1853,0.09464,0.08};   return polint4(xx,5,Zpmass,sigB);} //azz == -0.001
//else if (jj==21) {double sigB[5]={0.7801,0.3747,0.1814,0.09418,0.08};   return polint4(xx,5,Zpmass,sigB);} //azz ==  0.001
else if (jj==8) {double sigB[5]={0.7940,0.3755,0.1849,0.09308,0.073};  return polint4(xx,5,Zpmass,sigB);} //te6 == 0.9708
//else if (jj==23) {double sigB[5]={0.8015,0.3719,0.1861,0.09409,0.08};   return polint4(xx,5,Zpmass,sigB);} //azz == -0.001
//else if (jj==24) {double sigB[5]={0.7906,0.3785,0.1860,0.09345,0.08};   return polint4(xx,5,Zpmass,sigB);} //azz ==  0.001
else if (jj==9) {double sigB[5]={0.8095,0.3887,0.1893,0.09483,0.08};   return polint4(xx,5,Zpmass,sigB);} //te6 == 0.8708
//else if (jj==26) {double sigB[5]={0.8105,0.3836,0.1919,0.09602,0.08};   return polint4(xx,5,Zpmass,sigB);} //azz == -0.001
//else if (jj==27) {double sigB[5]={0.8175,0.3896,0.1933,0.09417,0.08};   return polint4(xx,5,Zpmass,sigB);} //azz ==  0.001
else if (jj==10) {double sigB[5]={0.8553,0.4013,0.1951,0.09893,0.08};   return polint4(xx,5,Zpmass,sigB);} //te6 == 0.7708
//else if (jj==29) {double sigB[5]={0.8478,0.4047,0.1965,0.1007,0.08};    return polint4(xx,5,Zpmass,sigB);} //azz == -0.001
//else if (jj==30) {double sigB[5]={0.8629,0.4063,0.1964,0.1004,0.08};    return polint4(xx,5,Zpmass,sigB);} //azz ==  0.001
else if (jj==11) {double sigB[5]={0.8904,0.4225,0.2084,0.1063,0.08};    return polint4(xx,5,Zpmass,sigB);} //te6 == 0.6708
//else if (jj==32) {double sigB[5]={0.8902,0.4248,0.2074,0.1060,0.085};   return polint4(xx,5,Zpmass,sigB);} //azz == -0.001
//else if (jj==33) {double sigB[5]={0.8921,0.4224,0.2095,0.1068,0.08};    return polint4(xx,5,Zpmass,sigB);} //azz ==  0.001
else if (jj==12) {double sigB[5]={0.9401,0.4546,0.2197,0.1142,0.085};   return polint4(xx,5,Zpmass,sigB);} //te6 == 0.5708
//else if (jj==35) {double sigB[5]={0.9488,0.4480,0.2219,0.1138,0.09};    return polint4(xx,5,Zpmass,sigB);} //azz == -0.001
//else if (jj==36) {double sigB[5]={0.9524,0.4525,0.2244,0.1152,0.09};    return polint4(xx,5,Zpmass,sigB);} //azz ==  0.001
else if (jj==13) {double sigB[5]={1.020,0.4835,0.2396,0.1233,0.09};     return polint4(xx,5,Zpmass,sigB);} //te6 == 0.4708
//else if (jj==38) {double sigB[5]={1.024,0.4832,0.2387,0.1237,0.095};    return polint4(xx,5,Zpmass,sigB);} //azz == -0.001
//else if (jj==39) {double sigB[5]={1.005,0.4864,0.2385,0.1227,0.09};     return polint4(xx,5,Zpmass,sigB);} //azz ==  0.001
else if (jj==14) {double sigB[5]={1.034,0.5043,0.2555,0.1280,0.095};    return polint4(xx,5,Zpmass,sigB);} //te6 == 0.4064
//else if (jj==41) {double sigB[5]={1.047,0.5076,0.2516,0.1298,0.10};     return polint4(xx,5,Zpmass,sigB);} //azz == -0.001
//else if (jj==42) {double sigB[5]={1.053,0.5040,0.2533,0.1293,0.105};    return polint4(xx,5,Zpmass,sigB);} //azz ==  0.001
else if (jj==15) {double sigB[5]={1.071,0.5216,0.2567,0.1332,0.10};     return polint4(xx,5,Zpmass,sigB);} //te6 == 0.3708
//else if (jj==44) {double sigB[5]={1.077,0.5166,0.2533,0.1342,0.105};    return polint4(xx,5,Zpmass,sigB);} //azz == -0.001
//else if (jj==45) {double sigB[5]={1.079,0.5228,0.2604,0.1338,0.115};    return polint4(xx,5,Zpmass,sigB);} //azz ==  0.001
else if (jj==16) {double sigB[5]={1.137,0.5560,0.2747,0.1442,0.115};    return polint4(xx,5,Zpmass,sigB);} //te6 == 0.2708
//else if (jj==47) {double sigB[5]={1.125,0.5528,0.2785,0.1447,0.11};     return polint4(xx,5,Zpmass,sigB);} //azz == -0.001
//else if (jj==48) {double sigB[5]={1.143,0.5490,0.2802,0.1457,0.115};    return polint4(xx,5,Zpmass,sigB);}  //azz ==  0.001
else if (jj==17) {double sigB[5]={1.188,0.5829,0.2951,0.1529,0.1122};   return polint4(xx,5,Zpmass,sigB);} //te6 == 0.1708
//else if (jj==50) {double sigB[5]={1.205,0.5927,0.2949,0.1542,0.1130};   return polint4(xx,5,Zpmass,sigB);} //azz == -0.001
//else if (jj==51) {double sigB[5]={1.207,0.5898,0.2955,0.1546,0.1148};   return polint4(xx,5,Zpmass,sigB);} //azz ==  0.001
else if (jj==18) {double sigB[5]={1.246,0.6056,0.3059,0.1613,0.1198};   return polint4(xx,5,Zpmass,sigB);} //te6 == 0.0708
//else if (jj==53) {double sigB[5]={1.243,0.6153,0.3101,0.1630,0.1192};   return polint4(xx,5,Zpmass,sigB);} //azz == -0.001
//else if (jj==54) {double sigB[5]={1.252,0.6252,0.3122,0.1634,0.1204};   return polint4(xx,5,Zpmass,sigB);} //azz ==  0.001
else if (jj==19) {double sigB[5]={1.276,0.6310,0.3175,0.1694,0.1245};   return polint4(xx,5,Zpmass,sigB);} //te6 == 0
//else if (jj==56) {double sigB[5]={1.280,0.6254,0.3196,0.1677,0.1243};   return polint4(xx,5,Zpmass,sigB);} //azz == -0.001
//else if (jj==57) {double sigB[5]={1.278,0.6271,0.3205,0.1685,0.1243};   return polint4(xx,5,Zpmass,sigB);}  //azz ==  0.001// NOW ONLY FOR azz == 0
else if (jj==20) {double sigB[5]={1.301,0.6413,0.3259,0.1726,0.1256};   return polint4(xx,5,Zpmass,sigB);} //te6 == -0.0708
else if (jj==21) {double sigB[5]={1.307,0.6506,0.3337,0.1754,0.1283};   return polint4(xx,5,Zpmass,sigB);} //te6 == -0.1708
else if (jj==22) {double sigB[5]={1.308,0.6551,0.3328,0.1759,0.1290};   return polint4(xx,5,Zpmass,sigB);} //te6 == -0.2708
else if (jj==23) {double sigB[5]={1.267,0.6452,0.3260,0.1721,0.1272};   return polint4(xx,5,Zpmass,sigB);} //te6 == -0.3708
else if (jj==24) {double sigB[5]={1.240,0.6192,0.3201,0.1676,0.1228};   return polint4(xx,5,Zpmass,sigB);} //te6 == -0.4708
else if (jj==25) {double sigB[5]={1.163,0.5955,0.3022,0.1596,0.1162};   return polint4(xx,5,Zpmass,sigB);} //te6 == -0.5708
else if (jj==26) {double sigB[5]={1.076,0.5461,0.2809,0.1482,0.1081};   return polint4(xx,5,Zpmass,sigB);} //te6 == -0.6708
else if (jj==27) {double sigB[5]={1.012,0.5022,0.2597,0.1341,0.0996};   return polint4(xx,5,Zpmass,sigB);} //te6 == -0.7708
else if (jj==28) {double sigB[5]={0.8948,0.4588,0.2353,0.1229,0.09079}; return polint4(xx,5,Zpmass,sigB);} //te6 == -0.8708
else if (jj==29) {double sigB[5]={0.8791,0.4387,0.2279,0.1181,0.08666}; return polint4(xx,5,Zpmass,sigB);} //te6 == -0.9113
else if (jj==30) {double sigB[5]={0.8145,0.4138,0.2104,0.1114,0.08052}; return polint4(xx,5,Zpmass,sigB);} //te6 == -0.9708
else if (jj==31) {double sigB[5]={0.7446,0.3809,0.1962,0.1024,0.07422}; return polint4(xx,5,Zpmass,sigB);} //te6 == -1.0708
else if (jj==32) {double sigB[5]={0.6964,0.3518,0.1825,0.09542,0.06946};return polint4(xx,5,Zpmass,sigB);} //te6 == -1.1708
else if (jj==33) {double sigB[5]={0.6766,0.3459,0.1746,0.09178,0.06622};return polint4(xx,5,Zpmass,sigB);} //te6 == -1.2708
else if (jj==34) {double sigB[5]={0.6715,0.3426,0.1776,0.09175,0.067};  return polint4(xx,5,Zpmass,sigB);} //te6 == -1.3708
else if (jj==35) {double sigB[5]={0.698,0.3521,0.1789,0.09364,0.06739}; return polint4(xx,5,Zpmass,sigB);} //te6 == -1.4708
else if (jj==36) {double sigB[5]={0.7166,0.3641,0.1809,0.09376,0.06856};return polint4(xx,5,Zpmass,sigB);} //te6 == -1.5708
else {printf("**** ERROR **** Second argument can take an integer value from 1 to 36. Please try one of these values.\n");exit(1);}
}


// Get the cubic interpolation of either sigma*Br(Z2 -> l+l-)_exp (in fb) or MZ2_lim for any tE6 in the range [-pi/2,pi/2].
// *********************************************
// ***************** CAUTION : *****************
// *********************************************
// This derivation of sigma*Br(Z2 -> l+l-)_exp as a function of tE6 assume that sigma*Br(Z2 -> l+l-)_exp = alpha*sigma*Br(Z2 -> l+l-)_th
// where alpha is a real.

double tE6vssBorMZ2lim(double xx, int jj)
{
if (xx>1.5708 || xx<-1.5708) {printf("**** ERROR **** First argument outside range [-pi/2,pi/2]. Please try a value in this range.\n");exit(1);}
double tE6[36]={1.5708,1.4708,1.3708,1.3181,1.2708,1.1708,1.0708,0.9708,0.8708,0.7708,0.6708,0.5708,0.4708,0.4064,0.3708,0.2708,0.1708,0.0708,0,-0.0708,-0.1708,-0.2708,-0.3708,-0.4708,-0.5708,-0.6708,-0.7708,-0.8708,-0.9113,-0.9708,-1.0708,-1.1708,-1.2708,-1.3708,-1.4708,-1.5708};
if (jj==0)       {
double sigB[36]={0.145,0.145697,0.145943,0.145372,0.145799,0.144986,0.145051,0.144545,0.145358,0.147262,0.150686,0.154356,0.158584,0.160767,0.163183,0.168293,0.172335,0.176237,0.18,0.181480,0.182776,0.183006,0.181246,0.179166,0.175464,0.170189,0.163665,0.158483,0.156263,0.153163,0.148998,0.145768,0.144084,0.144070,0.144945,0.145};return polint4(xx,36,tE6,sigB);
}
else if (jj==1)  {double Mlim[36]={2463.32,2463.89,2468.15,2468.25,2467.66,2470.06,2468.25,2467.72,2470.15,2476.82,2492.34,2502.79,2521.22,2531.83,2533.93,2546.02,2561.96,2570.58,2579.57,2583.58,2586.68,2586.97,2582.80,2578.64,2569.11,2555.35,2537.75,2519.31,2513.29,2497.51,2483.12,2468.06,2457.68,2462.36,2463.26,2464.55};return polint4(xx,36,tE6,Mlim);
}
else {printf("**** ERROR **** Second argument is either 0 or 1. Please try only one of these values.\n");exit(1);}
}


// Get the cubic interpolation of sigma*Br(Z2 -> l+l-)_th (in fb) for any tE6 in the range [-pi/2,pi/2] and any Z2 mass between 2 and 2.7 TeV.

double sBthForAnytE6MZ2(double xx, double tE6)
{
if (xx>2700 || xx<2000) {printf("**** ERROR **** First argument outside range [2000,2700] GeV. Please try a value in this range.\n");exit(1);}
else if (tE6>1.5708 || tE6<-1.5708) {printf("**** ERROR **** First argument outside range [-pi/2,pi/2]. Please try a value in this range.\n");exit(1);}
double Zpmass[5]={2000,2200,2400,2600,2700};
double sigBth[3];
double tE6th[3];
if      (tE6<=1.5708 && tE6>1.4708)
{
double sigBA[5]={0.71815,0.36425,0.1823,0.09406,0.08};sigBth[0] = polint4(xx,5,Zpmass,sigBA);tE6th[0]=1.5708;
double sigBB[5]={0.7401,0.3650,0.1827,0.09556,0.08};sigBth[1] = polint4(xx,5,Zpmass,sigBB);tE6th[1]=1.4708;
double sigBC[5]={0.7566,0.3652,0.1858,0.09609,0.08};sigBth[2] = polint4(xx,5,Zpmass,sigBC);tE6th[2]=1.3708;
return polint4(tE6,3,tE6th,sigBth);
}
else if      (tE6<=1.4708 && tE6>1.3708)
{
double sigBA[5]={0.71815,0.36425,0.1823,0.09406,0.08};sigBth[0] = polint4(xx,5,Zpmass,sigBA);tE6th[0]=1.5708;
double sigBB[5]={0.7401,0.3650,0.1827,0.09556,0.08};sigBth[1] = polint4(xx,5,Zpmass,sigBB);tE6th[1]=1.4708;
double sigBC[5]={0.7566,0.3652,0.1858,0.09609,0.08};sigBth[2] = polint4(xx,5,Zpmass,sigBC);tE6th[2]=1.3708;
return polint4(tE6,3,tE6th,sigBth);
}
else if      (tE6<=1.3708 && tE6>1.3181)
{
double sigBA[5]={0.7401,0.3650,0.1827,0.09556,0.08};sigBth[0] = polint4(xx,5,Zpmass,sigBA);tE6th[0]=1.4708;
double sigBB[5]={0.7566,0.3652,0.1858,0.09609,0.08};sigBth[1] = polint4(xx,5,Zpmass,sigBB);tE6th[1]=1.3708;
double sigBC[5]={0.7611,0.3731,0.1866,0.09486,0.08};sigBth[2] = polint4(xx,5,Zpmass,sigBC);tE6th[2]=1.3181;
return polint4(tE6,3,tE6th,sigBth);
}
else if      (tE6<=1.3181 && tE6>1.2708)
{
double sigBA[5]={0.7566,0.3652,0.1858,0.09609,0.08};sigBth[0] = polint4(xx,5,Zpmass,sigBA);tE6th[0]=1.3708;
double sigBB[5]={0.7611,0.3731,0.1866,0.09486,0.08};sigBth[1] = polint4(xx,5,Zpmass,sigBB);tE6th[1]=1.3181;
double sigBC[5]={0.7613,0.3749,0.1862,0.09578,0.08};sigBth[2] = polint4(xx,5,Zpmass,sigBC);tE6th[2]=1.2708;
return polint4(tE6,3,tE6th,sigBth);
}
else if      (tE6<=1.2708 && tE6>1.1708)
{
double sigBA[5]={0.7611,0.3731,0.1866,0.09486,0.08};sigBth[1] = polint4(xx,5,Zpmass,sigBA);tE6th[1]=1.3181;
double sigBB[5]={0.7613,0.3749,0.1862,0.09578,0.08};sigBth[1] = polint4(xx,5,Zpmass,sigBB);tE6th[1]=1.2708;
double sigBC[5]={0.7640,0.3707,0.1867,0.09403,0.075};sigBth[2] = polint4(xx,5,Zpmass,sigBC);tE6th[2]=1.1708;
return polint4(tE6,3,tE6th,sigBth);
}
else if      (tE6<=1.1708 && tE6>1.0708)
{
double sigBA[5]={0.7613,0.3749,0.1862,0.09578,0.08};sigBth[0] = polint4(xx,5,Zpmass,sigBA);tE6th[0]=1.2708;
double sigBB[5]={0.7640,0.3707,0.1867,0.09403,0.075};sigBth[1] = polint4(xx,5,Zpmass,sigBB);tE6th[1]=1.1708;
double sigBC[5]={0.7809,0.3720,0.1854,0.09417,0.075};sigBth[2] = polint4(xx,5,Zpmass,sigBC);tE6th[2]=1.0708;
return polint4(tE6,3,tE6th,sigBth);
}
else if      (tE6<=1.0708 && tE6>0.9708)
{
double sigBA[5]={0.7640,0.3707,0.1867,0.09403,0.075};sigBth[0] = polint4(xx,5,Zpmass,sigBA);tE6th[0]=1.1708;
double sigBB[5]={0.7809,0.3720,0.1854,0.09417,0.075};sigBth[1] = polint4(xx,5,Zpmass,sigBB);tE6th[1]=1.0708;
double sigBC[5]={0.7940,0.3755,0.1849,0.09308,0.073};sigBth[2] = polint4(xx,5,Zpmass,sigBC);tE6th[2]=0.9708;
return polint4(tE6,3,tE6th,sigBth);
}
else if      (tE6<=0.9708 && tE6>0.8708)
{
double sigBA[5]={0.7809,0.3720,0.1854,0.09417,0.075};sigBth[0] = polint4(xx,5,Zpmass,sigBA);tE6th[0]=1.0708;
double sigBB[5]={0.7940,0.3755,0.1849,0.09308,0.073};sigBth[1] = polint4(xx,5,Zpmass,sigBB);tE6th[1]=0.9708;
double sigBC[5]={0.8095,0.3887,0.1893,0.09483,0.08};sigBth[2] = polint4(xx,5,Zpmass,sigBC);tE6th[2]=0.8708;
return polint4(tE6,3,tE6th,sigBth);
}
else if      (tE6<=0.8708 && tE6>0.7708)
{
double sigBA[5]={0.7940,0.3755,0.1849,0.09308,0.073};sigBth[0] = polint4(xx,5,Zpmass,sigBA);tE6th[0]=0.9708;
double sigBB[5]={0.8095,0.3887,0.1893,0.09483,0.08};sigBth[1] = polint4(xx,5,Zpmass,sigBB);tE6th[1]=0.8708;
double sigBC[5]={0.8553,0.4013,0.1951,0.09893,0.08};sigBth[2] = polint4(xx,5,Zpmass,sigBC);tE6th[2]=0.7708;
return polint4(tE6,3,tE6th,sigBth);
}
else if      (tE6<=0.7708 && tE6>0.6708)
{
double sigBA[5]={0.8095,0.3887,0.1893,0.09483,0.08};sigBth[0] = polint4(xx,5,Zpmass,sigBA);tE6th[0]=0.8708;
double sigBB[5]={0.8553,0.4013,0.1951,0.09893,0.08};sigBth[1] = polint4(xx,5,Zpmass,sigBB);tE6th[1]=0.7708;
double sigBC[5]={0.8904,0.4225,0.2084,0.1063,0.08};sigBth[2] = polint4(xx,5,Zpmass,sigBC);tE6th[2]=0.6708;
return polint4(tE6,3,tE6th,sigBth);
}
else if      (tE6<=0.6708 && tE6>0.5708)
{
double sigBA[5]={0.8553,0.4013,0.1951,0.09893,0.08};sigBth[0] = polint4(xx,5,Zpmass,sigBA);tE6th[0]=0.7708;
double sigBB[5]={0.8904,0.4225,0.2084,0.1063,0.08};sigBth[1] = polint4(xx,5,Zpmass,sigBB);tE6th[1]=0.6708;
double sigBC[5]={0.9401,0.4546,0.2197,0.1142,0.085};sigBth[2] = polint4(xx,5,Zpmass,sigBC);tE6th[2]=0.5708;
return polint4(tE6,3,tE6th,sigBth);
}
else if      (tE6<=0.5708 && tE6>0.4708)
{
double sigBA[5]={0.8904,0.4225,0.2084,0.1063,0.08};sigBth[0] = polint4(xx,5,Zpmass,sigBA);tE6th[0]=0.6708;
double sigBB[5]={0.9401,0.4546,0.2197,0.1142,0.085};sigBth[1] = polint4(xx,5,Zpmass,sigBB);tE6th[1]=0.5708;
double sigBC[5]={1.020,0.4835,0.2396,0.1233,0.09};sigBth[2] = polint4(xx,5,Zpmass,sigBC);tE6th[2]=0.4708;
return polint4(tE6,3,tE6th,sigBth);
}
else if      (tE6<=0.4708 && tE6>0.4064)
{
double sigBA[5]={0.9401,0.4546,0.2197,0.1142,0.085};sigBth[0] = polint4(xx,5,Zpmass,sigBA);tE6th[0]=0.5708;
double sigBB[5]={1.020,0.4835,0.2396,0.1233,0.09};sigBth[1] = polint4(xx,5,Zpmass,sigBB);tE6th[1]=0.4708;
double sigBC[5]={1.034,0.5043,0.2555,0.1280,0.095};sigBth[2] = polint4(xx,5,Zpmass,sigBC);tE6th[2]=0.4064;
return polint4(tE6,3,tE6th,sigBth);
}
else if      (tE6<=0.4064 && tE6>0.3708)
{
double sigBA[5]={1.020,0.4835,0.2396,0.1233,0.09};sigBth[0] = polint4(xx,5,Zpmass,sigBA);tE6th[0]=0.4708;
double sigBB[5]={1.034,0.5043,0.2555,0.1280,0.095};sigBth[1] = polint4(xx,5,Zpmass,sigBB);tE6th[1]=0.4064;
double sigBC[5]={1.071,0.5216,0.2567,0.1332,0.10};sigBth[2] = polint4(xx,5,Zpmass,sigBC);tE6th[2]=0.3708;
return polint4(tE6,3,tE6th,sigBth);
}
else if      (tE6<=0.3708 && tE6>0.2708)
{
double sigBA[5]={1.034,0.5043,0.2555,0.1280,0.095};sigBth[0] = polint4(xx,5,Zpmass,sigBA);tE6th[0]=0.4064;
double sigBB[5]={1.071,0.5216,0.2567,0.1332,0.10};sigBth[1] = polint4(xx,5,Zpmass,sigBB);tE6th[1]=0.3708;
double sigBC[5]={1.137,0.5560,0.2747,0.1442,0.115};sigBth[2] = polint4(xx,5,Zpmass,sigBC);tE6th[2]=0.2708;
return polint4(tE6,3,tE6th,sigBth);
}
else if      (tE6<=0.2708 && tE6>0.1708)
{
double sigBA[5]={1.071,0.5216,0.2567,0.1332,0.10};sigBth[0] = polint4(xx,5,Zpmass,sigBA);tE6th[0]=0.3708;
double sigBB[5]={1.137,0.5560,0.2747,0.1442,0.115};sigBth[1] = polint4(xx,5,Zpmass,sigBB);tE6th[1]=0.2708;
double sigBC[5]={1.188,0.5829,0.2951,0.1529,0.1122};sigBth[2] = polint4(xx,5,Zpmass,sigBC);tE6th[2]=0.1708;
return polint4(tE6,3,tE6th,sigBth);
}
else if      (tE6<=0.1708 && tE6>0.0708)
{
double sigBA[5]={1.137,0.5560,0.2747,0.1442,0.115};sigBth[0] = polint4(xx,5,Zpmass,sigBA);tE6th[0]=0.2708;
double sigBB[5]={1.188,0.5829,0.2951,0.1529,0.1122};sigBth[1] = polint4(xx,5,Zpmass,sigBB);tE6th[1]=0.1708;
double sigBC[5]={1.246,0.6056,0.3059,0.1613,0.1198};sigBth[2] = polint4(xx,5,Zpmass,sigBC);tE6th[2]=0.0708;
return polint4(tE6,3,tE6th,sigBth);
}
else if      (tE6<=0.0708 && tE6>0)
{
double sigBA[5]={1.188,0.5829,0.2951,0.1529,0.1122};sigBth[0] = polint4(xx,5,Zpmass,sigBA);tE6th[0]=0.1708;
double sigBB[5]={1.246,0.6056,0.3059,0.1613,0.1198};sigBth[1] = polint4(xx,5,Zpmass,sigBB);tE6th[1]=0.0708;
double sigBC[5]={1.276,0.6310,0.3175,0.1694,0.1245};sigBth[2] = polint4(xx,5,Zpmass,sigBC);tE6th[2]=0;
return polint4(tE6,3,tE6th,sigBth);
}
else if      (tE6<=0 && tE6>-0.0708)
{
double sigBA[5]={1.246,0.6056,0.3059,0.1613,0.1198};sigBth[0] = polint4(xx,5,Zpmass,sigBA);tE6th[0]=0.0708;
double sigBB[5]={1.276,0.6310,0.3175,0.1694,0.1245};sigBth[1] = polint4(xx,5,Zpmass,sigBB);tE6th[1]=0;
double sigBC[5]={1.301,0.6413,0.3259,0.1726,0.1256};sigBth[2] = polint4(xx,5,Zpmass,sigBC);tE6th[2]=-0.0708;
return polint4(tE6,3,tE6th,sigBth);
}
else if      (tE6<=-0.0708 && tE6>-0.1708)
{
double sigBA[5]={1.276,0.6310,0.3175,0.1694,0.1245};sigBth[0] = polint4(xx,5,Zpmass,sigBA);tE6th[0]=0;
double sigBB[5]={1.301,0.6413,0.3259,0.1726,0.1256};sigBth[1] = polint4(xx,5,Zpmass,sigBB);tE6th[1]=-0.0708;
double sigBC[5]={1.307,0.6506,0.3337,0.1754,0.1283};sigBth[2] = polint4(xx,5,Zpmass,sigBC);tE6th[2]=-0.1708;
return polint4(tE6,3,tE6th,sigBth);
}
else if      (tE6<=-0.1708 && tE6>-0.2708)
{
double sigBA[5]={1.301,0.6413,0.3259,0.1726,0.1256};sigBth[0] = polint4(xx,5,Zpmass,sigBA);tE6th[0]=-0.0708;
double sigBB[5]={1.307,0.6506,0.3337,0.1754,0.1283};sigBth[1] = polint4(xx,5,Zpmass,sigBB);tE6th[1]=-0.1708;
double sigBC[5]={1.308,0.6551,0.3328,0.1759,0.1290};sigBth[2] = polint4(xx,5,Zpmass,sigBC);tE6th[2]=-0.2708;
return polint4(tE6,3,tE6th,sigBth);
}
else if      (tE6<=-0.2708 && tE6>-0.3708)
{
double sigBA[5]={1.307,0.6506,0.3337,0.1754,0.1283};sigBth[0] = polint4(xx,5,Zpmass,sigBA);tE6th[0]=-0.1708;
double sigBB[5]={1.308,0.6551,0.3328,0.1759,0.1290};sigBth[1] = polint4(xx,5,Zpmass,sigBB);tE6th[1]=-0.2708;
double sigBC[5]={1.267,0.6452,0.3260,0.1721,0.1272};sigBth[2] = polint4(xx,5,Zpmass,sigBC);tE6th[2]=-0.3708;
return polint4(tE6,3,tE6th,sigBth);
}
else if      (tE6<=-0.3708 && tE6>-0.4708)
{
double sigBA[5]={1.308,0.6551,0.3328,0.1759,0.1290};sigBth[0] = polint4(xx,5,Zpmass,sigBA);tE6th[0]=-0.2708;
double sigBB[5]={1.267,0.6452,0.3260,0.1721,0.1272};sigBth[1] = polint4(xx,5,Zpmass,sigBB);tE6th[1]=-0.3708;
double sigBC[5]={1.240,0.6192,0.3201,0.1676,0.1228};sigBth[2] = polint4(xx,5,Zpmass,sigBC);tE6th[2]=-0.4708;
return polint4(tE6,3,tE6th,sigBth);
}
else if      (tE6<=-0.4708 && tE6>-0.5708)
{
double sigBA[5]={1.267,0.6452,0.3260,0.1721,0.1272};sigBth[0] = polint4(xx,5,Zpmass,sigBA);tE6th[0]=-0.3708;
double sigBB[5]={1.240,0.6192,0.3201,0.1676,0.1228};sigBth[1] = polint4(xx,5,Zpmass,sigBB);tE6th[1]=-0.4708;
double sigBC[5]={1.163,0.5955,0.3022,0.1596,0.1162};sigBth[2] = polint4(xx,5,Zpmass,sigBC);tE6th[2]=-0.5708;
return polint4(tE6,3,tE6th,sigBth);
}
else if      (tE6<=-0.5708 && tE6>-0.6708)
{
double sigBA[5]={1.240,0.6192,0.3201,0.1676,0.1228};sigBth[0] = polint4(xx,5,Zpmass,sigBA);tE6th[0]=-0.4708;
double sigBB[5]={1.163,0.5955,0.3022,0.1596,0.1162};sigBth[1] = polint4(xx,5,Zpmass,sigBB);tE6th[1]=-0.5708;
double sigBC[5]={1.076,0.5461,0.2809,0.1482,0.1081};sigBth[2] = polint4(xx,5,Zpmass,sigBC);tE6th[2]=-0.6708;
return polint4(tE6,3,tE6th,sigBth);
}
else if      (tE6<=-0.6708 && tE6>-0.7708)
{
double sigBA[5]={1.163,0.5955,0.3022,0.1596,0.1162};sigBth[0] = polint4(xx,5,Zpmass,sigBA);tE6th[0]=-0.5708;
double sigBB[5]={1.076,0.5461,0.2809,0.1482,0.1081};sigBth[1] = polint4(xx,5,Zpmass,sigBB);tE6th[1]=-0.6708;
double sigBC[5]={1.012,0.5022,0.2597,0.1341,0.0996};sigBth[2] = polint4(xx,5,Zpmass,sigBC);tE6th[2]=-0.7708;
return polint4(tE6,3,tE6th,sigBth);
}
else if      (tE6<=-0.7708 && tE6>-0.8708)
{
double sigBA[5]={1.076,0.5461,0.2809,0.1482,0.1081};sigBth[0] = polint4(xx,5,Zpmass,sigBA);tE6th[0]=-0.6708;
double sigBB[5]={1.012,0.5022,0.2597,0.1341,0.0996};sigBth[1] = polint4(xx,5,Zpmass,sigBB);tE6th[1]=-0.7708;
double sigBC[5]={0.8948,0.4588,0.2353,0.1229,0.09079};sigBth[2] = polint4(xx,5,Zpmass,sigBC);tE6th[2]=-0.8708;
return polint4(tE6,3,tE6th,sigBth);
}
else if      (tE6<=-0.8708 && tE6>-0.9113)
{
double sigBA[5]={1.012,0.5022,0.2597,0.1341,0.0996};sigBth[0] = polint4(xx,5,Zpmass,sigBA);tE6th[0]=-0.7708;
double sigBB[5]={0.8948,0.4588,0.2353,0.1229,0.09079};sigBth[1] = polint4(xx,5,Zpmass,sigBB);tE6th[1]=-0.8708;
double sigBC[5]={0.8791,0.4387,0.2279,0.1181,0.08666};sigBth[2] = polint4(xx,5,Zpmass,sigBC);tE6th[2]=-0.9113;
return polint4(tE6,3,tE6th,sigBth);
}
else if      (tE6<=-0.9113 && tE6>-0.9708)
{
double sigBA[5]={0.8948,0.4588,0.2353,0.1229,0.09079};sigBth[0] = polint4(xx,5,Zpmass,sigBA);tE6th[0]=-0.8708;
double sigBB[5]={0.8791,0.4387,0.2279,0.1181,0.08666};sigBth[1] = polint4(xx,5,Zpmass,sigBB);tE6th[1]=-0.9113;
double sigBC[5]={0.8145,0.4138,0.2104,0.1114,0.08052};sigBth[2] = polint4(xx,5,Zpmass,sigBC);tE6th[2]=-0.9708;
return polint4(tE6,3,tE6th,sigBth);
}
else if      (tE6<=-0.9708 && tE6>-1.0708)
{
double sigBA[5]={0.8791,0.4387,0.2279,0.1181,0.08666};sigBth[0] = polint4(xx,5,Zpmass,sigBA);tE6th[0]=-0.9113;
double sigBB[5]={0.8145,0.4138,0.2104,0.1114,0.08052};sigBth[1] = polint4(xx,5,Zpmass,sigBB);tE6th[1]=-0.9708;
double sigBC[5]={0.7446,0.3809,0.1962,0.1024,0.07422};sigBth[2] = polint4(xx,5,Zpmass,sigBC);tE6th[2]=-1.0708;
return polint4(tE6,3,tE6th,sigBth);
}
else if      (tE6<=-1.0708 && tE6>-1.1708)
{
double sigBA[5]={0.8145,0.4138,0.2104,0.1114,0.08052};sigBth[0] = polint4(xx,5,Zpmass,sigBA);tE6th[0]=-0.9708;
double sigBB[5]={0.7446,0.3809,0.1962,0.1024,0.07422};sigBth[1] = polint4(xx,5,Zpmass,sigBB);tE6th[1]=-1.0708;
double sigBC[5]={0.6964,0.3518,0.1825,0.09542,0.06946};sigBth[2] = polint4(xx,5,Zpmass,sigBC);tE6th[2]=-1.1708;
return polint4(tE6,3,tE6th,sigBth);
}
else if      (tE6<=-1.1708 && tE6>-1.2708)
{
double sigBA[5]={0.7446,0.3809,0.1962,0.1024,0.07422};sigBth[0] = polint4(xx,5,Zpmass,sigBA);tE6th[0]=-1.0708;
double sigBB[5]={0.6964,0.3518,0.1825,0.09542,0.06946};sigBth[1] = polint4(xx,5,Zpmass,sigBB);tE6th[1]=-1.1708;
double sigBC[5]={0.6766,0.3459,0.1746,0.09178,0.06622};sigBth[2] = polint4(xx,5,Zpmass,sigBC);tE6th[2]=-1.2708;
return polint4(tE6,3,tE6th,sigBth);
}
else if      (tE6<=-1.2708 && tE6>-1.3708)
{
double sigBA[5]={0.6964,0.3518,0.1825,0.09542,0.06946};sigBth[0] = polint4(xx,5,Zpmass,sigBA);tE6th[0]=-1.1708;
double sigBB[5]={0.6766,0.3459,0.1746,0.09178,0.06622};sigBth[1] = polint4(xx,5,Zpmass,sigBB);tE6th[1]=-1.2708;
double sigBC[5]={0.6715,0.3426,0.1776,0.09175,0.067};sigBth[2] = polint4(xx,5,Zpmass,sigBC);tE6th[2]=-1.3708;
return polint4(tE6,3,tE6th,sigBth);
}
else if      (tE6<=-1.3708 && tE6>-1.4708)
{
double sigBA[5]={0.6766,0.3459,0.1746,0.09178,0.06622};sigBth[0] = polint4(xx,5,Zpmass,sigBA);tE6th[0]=-1.2708;
double sigBB[5]={0.6715,0.3426,0.1776,0.09175,0.067};sigBth[1] = polint4(xx,5,Zpmass,sigBB);tE6th[1]=-1.3708;
double sigBC[5]={0.698,0.3521,0.1789,0.09364,0.06739};sigBth[2] = polint4(xx,5,Zpmass,sigBC);tE6th[2]=-1.4708;
return polint4(tE6,3,tE6th,sigBth);
}
else
{
double sigBA[5]={0.6715,0.3426,0.1776,0.09175,0.067};sigBth[0] = polint4(xx,5,Zpmass,sigBA);tE6th[0]=-1.3708;
double sigBB[5]={0.698,0.3521,0.1789,0.09364,0.06739};sigBth[1] = polint4(xx,5,Zpmass,sigBB);tE6th[1]=-1.4708;
double sigBC[5]={0.7166,0.3641,0.1809,0.09376,0.06856};sigBth[2] = polint4(xx,5,Zpmass,sigBC);tE6th[2]=-1.5708;
return polint4(tE6,3,tE6th,sigBth);
}
}


int Zprimelimits(void)
{
double MZ2=findValW("MZ2");
printf("\nZ2 mass : %5.3f GeV \n",MZ2);    
printf("==================================\n");
printf("==== Limits on the Z2 boson : ====\n");
printf("==================================\n");
printf("Using Z' limits with 20 fb-1 at sqrt(s) = 8 TeV from\n");
printf("- ATLAS, Phys. Rev. D90 no. 5, (2014) 052005, http://arxiv.org/abs/1405.4123;\n");
printf("- CMS, JHEP 1504 (2015) 025, http://arxiv.org/abs/1412.6302:\n");
printf("==================================\n");
double tE6=findValW("tE6");
double Mlim=tE6vssBorMZ2lim(tE6,1)+104.;//shift of 44 GeV to get Mlim(psi)\sim2.51 TeV and Mlim(chi)\sim2.62 TeV (ATLAS_8TeV:1405.4123);
                                        //then shift of 60 GeV to get Mlim(psi) \sim 2.57 TeV (CMS_8TeV:1412.6302)
printf("Interpolated MZ2 limit for tE6=%+5.4f rad : %5.3f GeV\n",tE6,Mlim);



if(MZ2 >= Mlim)
{
printf("-> below the MZ2 chosen\n");
printf("==================================\n");
}

else
{
printf("-> above the MZ2 chosen\n");
printf("Check whether this point is safe after including BSM branching fractions\n");
char * pname;
txtList L;
double width,BrSM=0.,Brn=0.,BrSUSY=0.,Brh=0.,BrT=0.,sBexp=0.,sBth=0.;
pname = "Z2";
width=pWidth(pname,&L);
printf("%s-> :   total width=%E\n",pname,width);

BrSM=findBr(L,"e,E")+findBr(L,"m,M")+findBr(L,"l,L")+findBr(L,"ne,Ne")+findBr(L,"nm,Nm")+findBr(L,"nl,Nl")+findBr(L,"u,U")+findBr(L,"d,D")+findBr(L,"c,C")+findBr(L,"s,S")+findBr(L,"t,T")+findBr(L,"b,B")+findBr(L,"W+,W-");
printf("SM branching of Z2 : %4.2f %% (%4.2f %% leptonic)\n",100*BrSM,100*(findBr(L,"e,E") + findBr(L,"m,M")));

Brn=findBr(L,"ner,Ner")+findBr(L,"nmr,Nmr")+findBr(L,"nlr,Nlr");
printf("RH neutrino branching of Z2 : %4.2f %%\n",100*Brn);

BrSUSY=findBr(L,"~1+,~1-")+findBr(L,"~1+,~2-")+findBr(L,"~2+,~1-")+findBr(L,"~2+,~2-")+findBr(L,"~o1,~o1")+findBr(L,"~o1,~o2")+findBr(L,"~o1,~o3")+findBr(L,"~o1,~o4")+findBr(L,"~o1,~o5")+findBr(L,"~o1,~o6")+findBr(L,"~o2,~o2")+findBr(L,"~o2,~o3")+findBr(L,"~o2,~o4")+findBr(L,"~o2,~o5")+findBr(L,"~o2,~o6")+findBr(L,"~o3,~o3")+findBr(L,"~o3,~o4")+findBr(L,"~o3,~o5")+findBr(L,"~o3,~o6")+findBr(L,"~o4,~o4")+findBr(L,"~o4,~o5")+findBr(L,"~o4,~o6")+findBr(L,"~o5,~o5")+findBr(L,"~o5,~o6")+findBr(L,"~o6,~o6")+findBr(L,"~nE,~NE")+findBr(L,"~nM,~NM")+findBr(L,"~nL,~NL")+findBr(L,"~ne,~Ne")+findBr(L,"~nm,~Nm")+findBr(L,"~nl,~Nl")+findBr(L,"~t1,~T1")+findBr(L,"~t1,~T2")+findBr(L,"~t2,~T1")+findBr(L,"~t2,~T2")+findBr(L,"~b1,~B1")+findBr(L,"~b1,~B2")+findBr(L,"~b2,~B1")+findBr(L,"~b2,~B2")+findBr(L,"~uL,~UL")+findBr(L,"~uR,~UR")+findBr(L,"~dL,~DL")+findBr(L,"~dR,~DR")+findBr(L,"~cL,~CL")+findBr(L,"~cR,~CR")+findBr(L,"~sL,~SL")+findBr(L,"~sR,~SR")+findBr(L,"~l1,~L1")+findBr(L,"~l1,~L2")+findBr(L,"~l2,~L1")+findBr(L,"~l2,~L2")+findBr(L,"~mL,~ML")+findBr(L,"~mR,~MR")+findBr(L,"~eL,~EL")+findBr(L,"~eR,~ER");
printf("SUSY branching of Z2 : %4.2f %%\n",100*BrSUSY);

Brh=findBr(L,"Z1,h1")+findBr(L,"Z1,h2")+findBr(L,"h1,ha")+findBr(L,"h2,ha")+findBr(L,"W-,H+")+findBr(L,"W+,H-")+findBr(L,"H+,H-");
printf("Higgs or Higgs-gauge boson branching of Z2 : %4.2f %%\n",100*Brh);

BrT=BrSM+Brn+BrSUSY+Brh;
printf("Ratio B(Z2->UMSSM)/B(Z2->SM) = %5.3f \n\n",BrT/BrSM);


sBexp=tE6vssBorMZ2lim(tE6,0)*BrT/BrSM;
sBth=sBthForAnytE6MZ2(MZ2,tE6);

printf("sigma*Br(Z2 -> l+l-)_exp for UMSSM is %.3e fb,\n",sBexp);
printf("sigma*Br(Z2 -> l+l-)_th for (tE6,MZ2)=(%+5.4f rad, %5.3f GeV) is %.3e fb\n",tE6,MZ2,sBth);


if(sBexp >= sBth)
{
printf("==================================\n");
}
else
{
printf("Z2 too light, excluded by LHC limits on Z'\n");
printf("==================================\n");
return 1;
}
}
return 0;
}
