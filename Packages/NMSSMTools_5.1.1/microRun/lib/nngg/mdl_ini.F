*     LanHEP output produced at Thu Oct 23 23:48:20 2014
*     Model named 'NMSSM(../spect.slha)+hgg'

#include "DgetNEWv0.inc"
#include "CgetNEWv0.inc"


      double precision function myif(a,b,c)
      double precision a,b,c
      if(a.lt.0D0) then
      myif=b
      else
      myif=c
      endif
      return 
      end

      subroutine ModelDefaults
      implicit none

#include "model.h"

      EE = 0.31344D0
      MZ = 91.187D0
      MW = 80.385D0
      MW2=MW*MW
      GG = 1.117D0
      alfSMZ = 0.1172D0
      MbMb = 4.214D0
      McMc = 1.3D0
      Mtp = 173.1D0
      Ml = 1.777D0
      Mq = 0.05D0
      Q = 100D0
      Au = 0D0
      Ad = 0D0
      Mele = 0.000511D0
      Mmu = 0.1D0

      end


      subroutine ModelConstIni(mode,Wh,err )
      implicit none
      integer mode
      double precision  Wh(5)
      integer err
#include "model.h"

      double precision 
     > aPrintF,openAppend,System,slhaVal,
     > rDiagonal,MtRun,slhaValFormat,MassArray,MixMatrix,
     > MixMatrixU,slhaVal1,slhaVal2,slhaWidth,initQCD5,
     > MtEff,MbEff,McEff,alphaQCD,MbRun,myif,findValW

      call ModelDefaults

      hisW(1)=Wh(1) 
      hisW(2)=Wh(2)
      hisW(3)=Wh(3)
      hiaW(1)=Wh(4)
      hiaW(2)=Wh(5)


      CW = MW/MZ
      SW = sqrt(1D0-CW**2)
      C2W = CW**2-SW**2
      LamQCD = initQCD5(alfSMZ, McMc, MbMb, Mtp)

      if(mode.lt.2) then
        Mtp= slhaVal1('SMINPUTS',MZ,6)
        mu = slhaVal1('HMIX', MZ, 1)
        Lambda = slhaVal1('NMSSMRUN', MZ, 1)
        Kappa  = slhaVal1('NMSSMRUN', MZ, 2)
        aLambda= slhaVal1('NMSSMRUN', MZ, 3)
        aKappa = slhaVal1('NMSSMRUN', MZ, 4)
        muP = slhaVal1('NMSSMRUN', MZ, 8)
        At  = slhaVal2('Au', MZ, 3, 3)
        Ab  = slhaVal2('Ad', MZ, 3, 3)
        Al  = slhaVal2('Ae', MZ, 3, 3)
        Ml2 = slhaVal1('MSOFT', MZ, 32)
        Ml3 = slhaVal1('MSOFT', MZ, 33)
        Mr2 = slhaVal1('MSOFT', MZ, 35)
        Mr3 = slhaVal1('MSOFT', MZ, 36)
        Mq2 = slhaVal1('MSOFT', MZ, 42)
        Mq3 = slhaVal1('MSOFT', MZ, 43)
        Mu2 = slhaVal1('MSOFT', MZ, 45)
        Mu3 = slhaVal1('MSOFT', MZ, 46)
        Md2 = slhaVal1('MSOFT', MZ, 48)
        Md3 = slhaVal1('MSOFT', MZ, 49)
        tB  = slhaVal1('HMIX', MZ, 2)
        Mh1 = slhaVal1('MASS', MZ, 25) 
        Mh2 = slhaVal1('MASS', MZ, 35)
        Mh3 = slhaVal1('MASS', MZ, 45) 
        Mha = slhaVal1('MASS', MZ, 36)
        Mhb = slhaVal1('MASS', MZ, 46) 
        MHc = slhaVal1('MASS', MZ, 37) 
        MNE1= slhaVal1('MASS', MZ, 1000022)
        MNE2= slhaVal1('MASS', MZ, 1000023)
        MNE3= slhaVal1('MASS', MZ, 1000025)
        MNE4= slhaVal1('MASS', MZ, 1000035)
        MNE5= slhaVal1('MASS', MZ, 1000045)
        MC1 = slhaVal1('MASS', MZ, 1000024)
        MC2 = slhaVal1('MASS', MZ, 1000037)
        MSG = slhaVal1('MASS', MZ, 1000021)
        MSuL= slhaVal1('MASS', MZ, 1000002)
        MSuR= slhaVal1('MASS', MZ, 2000002)
        MSdL= slhaVal1('MASS', MZ, 1000001)
        MSdR= slhaVal1('MASS', MZ, 2000001)
        MScL= slhaVal1('MASS', MZ, 1000004)
        MScR = slhaVal1('MASS', MZ, 2000004)
        MSsL = slhaVal1('MASS', MZ, 1000003)
        MSsR = slhaVal1('MASS', MZ, 2000003)
        MSt1 = slhaVal1('MASS', MZ, 1000006)
        MSt2 = slhaVal1('MASS', MZ, 2000006)
        MSb1 = slhaVal1('MASS', MZ, 1000005)
        MSb2 = slhaVal1('MASS', MZ, 2000005)
        QSUSY = sqrt(MSt1*MSt2)
        MSeL = slhaVal1('MASS', MZ, 1000011)
        MSeR = slhaVal1('MASS', MZ, 2000011)
        MSmL = slhaVal1('MASS', MZ, 1000013)
        MSmR = slhaVal1('MASS', MZ, 2000013)
        MSl1 = slhaVal1('MASS', MZ, 1000015)
        MSl2 = slhaVal1('MASS', MZ, 2000015)
        MSne = slhaVal1('MASS', MZ, 1000012)
        MSnm = slhaVal1('MASS', MZ, 1000014)
        MSnl = slhaVal1('MASS', MZ, 1000016)
        Zl11 = slhaVal2('STAUMIX', QSUSY, 1, 1)
        Zl12 = slhaVal2('STAUMIX', QSUSY, 1, 2)
        Zl21 = slhaVal2('STAUMIX', QSUSY, 2, 1)
        Zl22 = slhaVal2('STAUMIX', QSUSY, 2, 2)
        Zh11 = slhaVal2('NMHMIX', QSUSY, 1, 1)
        Zh12 = slhaVal2('NMHMIX', QSUSY, 1, 2)
        Zh13 = slhaVal2('NMHMIX', QSUSY, 1, 3)
        Zh21 = slhaVal2('NMHMIX', QSUSY, 2, 1)
        Zh22 = slhaVal2('NMHMIX', QSUSY, 2, 2)
        Zh23 = slhaVal2('NMHMIX', QSUSY, 2, 3)
        Zh31 = slhaVal2('NMHMIX', QSUSY, 3, 1)
        Zh32 = slhaVal2('NMHMIX', QSUSY, 3, 2)
        Zh33 = slhaVal2('NMHMIX', QSUSY, 3, 3)
        Za11 = slhaVal2('NMAMIX', QSUSY, 1, 1)
        Za12 = slhaVal2('NMAMIX', QSUSY, 1, 2)
        Za13 = slhaVal2('NMAMIX', QSUSY, 1, 3)
        Za21 = slhaVal2('NMAMIX', QSUSY, 2, 1)
        Za22 = slhaVal2('NMAMIX', QSUSY, 2, 2)
        Za23 = slhaVal2('NMAMIX', QSUSY, 2, 3)
        Zn11 = slhaVal2('NMNMIX', QSUSY, 1, 1)
        Zn12 = slhaVal2('NMNMIX', QSUSY, 1, 2)
        Zn13 = slhaVal2('NMNMIX', QSUSY, 1, 3)
        Zn14 = slhaVal2('NMNMIX', QSUSY, 1, 4)
        Zn15 = slhaVal2('NMNMIX', QSUSY, 1, 5)
        Zn21 = slhaVal2('NMNMIX', QSUSY, 2, 1)
        Zn22 = slhaVal2('NMNMIX', QSUSY, 2, 2)
        Zn23 = slhaVal2('NMNMIX', QSUSY, 2, 3)
        Zn24 = slhaVal2('NMNMIX', QSUSY, 2, 4)
        Zn25 = slhaVal2('NMNMIX', QSUSY, 2, 5)
        Zn31 = slhaVal2('NMNMIX', QSUSY, 3, 1)
        Zn32 = slhaVal2('NMNMIX', QSUSY, 3, 2)
        Zn33 = slhaVal2('NMNMIX', QSUSY, 3, 3)
        Zn34 = slhaVal2('NMNMIX', QSUSY, 3, 4)
        Zn35 = slhaVal2('NMNMIX', QSUSY, 3, 5)
        Zn41 = slhaVal2('NMNMIX', QSUSY, 4, 1)
        Zn42 = slhaVal2('NMNMIX', QSUSY, 4, 2)
        Zn43 = slhaVal2('NMNMIX', QSUSY, 4, 3)
        Zn44 = slhaVal2('NMNMIX', QSUSY, 4, 4)
        Zn45 = slhaVal2('NMNMIX', QSUSY, 4, 5)
        Zn51 = slhaVal2('NMNMIX', QSUSY, 5, 1)
        Zn52 = slhaVal2('NMNMIX', QSUSY, 5, 2)
        Zn53 = slhaVal2('NMNMIX', QSUSY, 5, 3)
        Zn54 = slhaVal2('NMNMIX', QSUSY, 5, 4)
        Zn55 = slhaVal2('NMNMIX', QSUSY, 5, 5)
        Zu11 = slhaVal2('UMIX', QSUSY, 1, 1)
        Zu12 = slhaVal2('UMIX', QSUSY, 1, 2)
        Zu21 = slhaVal2('UMIX', QSUSY, 2, 1)
        Zu22 = slhaVal2('UMIX', QSUSY, 2, 2)
        Zv11 = slhaVal2('VMIX', QSUSY, 1, 1)
        Zv12 = slhaVal2('VMIX', QSUSY, 1, 2)
        Zv21 = slhaVal2('VMIX', QSUSY, 2, 1)
        Zv22 = slhaVal2('VMIX', QSUSY, 2, 2)
        Zt11 = slhaVal2('STOPMIX', QSUSY, 1, 1)
        Zt12 = slhaVal2('STOPMIX', QSUSY, 1, 2)
        Zt21 = slhaVal2('STOPMIX', QSUSY, 2, 1)
        Zt22 = slhaVal2('STOPMIX', QSUSY, 2, 2)
        Zb11 = slhaVal2('SBOTMIX', QSUSY, 1, 1)
        Zb12 = slhaVal2('SBOTMIX', QSUSY, 1, 2)
        Zb21 = slhaVal2('SBOTMIX', QSUSY, 2, 1)
        Zb22 = slhaVal2('SBOTMIX', QSUSY, 2, 2)
        la1 = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'L1 %lf')
        la2 = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'L2 %lf')
        la3 = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'L3 %lf')
        la4 = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'L4 %lf')
        la5 = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'L5 %lf')
        la6 = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'L6 %lf')
        la7 = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'L7 %lf')
        la1s = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'K1 %lf')
        la2s = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'K2 %lf')
        la3s = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'K3 %lf')
        la4s = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'K4 %lf')
        la5s = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'K5 %lf')
        la6s = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'K6 %lf')
        la7s = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'K7 %lf')
        la8s = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'K8 %lf')
        aa1 = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'A1 %lf')
        aa2 = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'A2 %lf')
        aa3 = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'A3  %lf')
        aa4 = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'A4 %lf')
        aa5 = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'A5 %lf')
        aa6 = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'A6 %lf')
        mB1 = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'B1 %lf')
        mB2 = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'B2 %lf')
        X = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'X %lf')
        dMb = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'DELMB %lf')
        xvev= slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'XVEV %lf')
        xif = slhaVal1('NMSSMRUN', QSUSY, 6)
        xis = slhaVal1('NMSSMRUN', QSUSY, 7)
        MM3 = slhaVal1('NMSSMRUN', QSUSY, 9)
        MSP = slhaVal1('NMSSMRUN', QSUSY, 12)
      else 
        Mtp=findValW('Mtp')
        mu=findValW('mu')
        Lambda  = findValW('Lambda')
        Kappa   = findValW('Kappa')
        aLambda = findValW('aLambda')
        aKappa  = findValW('aKappa') 
        muP = findValW('muP')
        At= findValW('At')
        Ab= findValW('Ab')
        Al= findValW('Al')
        Ml2=findValW('Ml2')
        Ml3=findValW('Ml3')  
        Mr2=findValW('Mr2')
        Mr3=findValW('Mr3')
        Mq2=findValW('Mq2') 
        Mq3=findValW('Mq2')
        Mu2=findValW('Mu2')  
        Mu3=findValW('Mu3')      
        Md2=findValW('Md2')
        Md3=findValW('Md3')
        tB=findValW('tB')
        Mh1  = findValW('Mh1')
        Mh2  = findValW('Mh2')
        Mh3  = findValW('Mh3')
        Mha  = findValW('Mha')
        Mhb  = findValW('Mhb')
        MHc  = findValW('MHc') 
        MNE1 =findValW('MNE1')
        MNE2 =findValW('MNE2')
        MNE3 =findValW('MNE3')
        MNE4 =findValW('MNE4')
        MNE5 =findValW('MNE5')
        MC1   =findValW('MC1')
        MC2   =findValW('MC2')
        MSG = findValW('MSG')
        MSuL= findValW('MSuL')
        MSuR= findValW('MSuR')
        MSdL= findValW('MSdL')
        MSdR= findValW('MSdR')
        MScL= findValW('MScL')
        MScR= findValW('MScR')
        MSsL= findValW('MSsL')
        MSsR= findValW('MSsR')
        MSt1= findValW('MSt1')
        MSt2= findValW('MSt2')      
        MSb1= findValW('MSb1')
        MSb2= findValW('MSb2')
        MSeL= findValW('MSeL')
        MSeR= findValW('MSeR')
        MSmL= findValW('MSmL')
        MSmR= findValW('MSmR')
        MSl1= findValW('MSl1')
        MSl2= findValW('MSl2')
        MSne= findValW('MSne')
        MSnm= findValW('MSnm')
        MSnl= findValW('MSnl')
        Zl11 = findValW('Zl11')
        Zl12 = findValW('Zl12')
        Zl21 = findValW('Zl21')
        Zl22 = findValW('Zl22')
        Zh11 = findValW('Zh11')
        Zh12 = findValW('Zh21')
        Zh13 = findValW('Zh31')
        Zh21 = findValW('Zh12')
        Zh22 = findValW('Zh22')
        Zh23 = findValW('Zh32')
        Zh31 = findValW('Zh13')
        Zh32 = findValW('Zh23')
        Zh33 = findValW('Zh33')
        Zh11 = findValW('Zh11')
        Zh12 = findValW('Zh12')
        Zh13 = findValW('Zh13')
        Zh21 = findValW('Zh21')
        Zh22 = findValW('Zh22')
        Zh23 = findValW('Zh23')
        Zh31 = findValW('Zh31')
        Zh32 = findValW('Zh32')
        Zh33 = findValW('Zh33')
        Za11 = findValW('Za11')
        Za12 = findValW('Za12')
        Za13 = findValW('Za13')
        Za21 = findValW('Za21')
        Za22 = findValW('Za22')
        Za23 = findValW('Za23')  
        Zn11 = findValW('Zn11')
        Zn12 = findValW('Zn12')
        Zn13 = findValW('Zn13')
        Zn14 = findValW('Zn14')
        Zn15 = findValW('Zn15')
        Zn21 = findValW('Zn21')
        Zn22 = findValW('Zn22')
        Zn23 = findValW('Zn23')
        Zn24 = findValW('Zn24')
        Zn25 = findValW('Zn25')
        Zn31 = findValW('Zn31')
        Zn32 = findValW('Zn32')
        Zn33 = findValW('Zn33')
        Zn34 = findValW('Zn34')
        Zn35 = findValW('Zn35')
        Zn41 = findValW('Zn41')
        Zn42 = findValW('Zn42')
        Zn43 = findValW('Zn43')
        Zn44 = findValW('Zn44')
        Zn45 = findValW('Zn45')
        Zn51 = findValW('Zn51')
        Zn52 = findValW('Zn52')
        Zn53 = findValW('Zn53')
        Zn54 = findValW('Zn54')
        Zn55 = findValW('Zn55')
        Zu11 = findValW('Zu11')
        Zu12 = findValW('Zu12')
        Zu21 = findValW('Zu21')
        Zu22 = findValW('Zu22')
        Zv11 = findValW('Zv11')
        Zv12 = findValW('Zv12')
        Zv21 = findValW('Zv21')
        Zv22 = findValW('Zv22')
        Zt11 = findValW('Zt11')
        Zt12 = findValW('Zt12')
        Zt21 = findValW('Zt21')
        Zt22 = findValW('Zt22')
        Zb11 = findValW('Zb11')
        Zb12 = findValW('Zb12')
        Zb21 = findValW('Zb21')
        Zb22 = findValW('Zb22')
        la1  = findValW('la1')
        la2  = findValW('la2')
        la3  = findValW('la3')
        la4  = findValW('la4')
        la5  = findValW('la5')
        la6  = findValW('la6')
        la7  = findValW('la7')
        la1s = findValW('la1s')
        la2s = findValW('la2s')
        la3s = findValW('la3s')
        la4s = findValW('la4s')
        la5s = findValW('la5s')
        la6s = findValW('la6s')
        la7s = findValW('la7s')
        la8s = findValW('la8s')
        aa1  = findValW('aa1')
        aa2  = findValW('aa2')
        aa3  = findValW('aa3')
        aa4  = findValW('aa4')
        aa5  = findValW('aa5')
        aa6  = findValW('aa6')
        mB1  = findValW('B1')
        mB2  = findValW('B2')
        X    = findValW('X')
        dMb  = findValW('dMb')
        xif  = findValW('xif')
        xis  = findValW('xis')
        MM3  = findValW('MM3')
        MSP  = findValW('MSP')
        xvev = findValW('xvev')
      endif
C      call  ModelDigest
      sb = tB/sqrt(1D0+tB**2)
      cb = sqrt(1D0-sb**2)
      t2b = 2D0*tB/(1D0-tB**2)
      Pa12 = Za13
      Pa22 = Za23
      Pa11 = myif(Za11*Pa22, Pa22, -Pa22)
      Pa21 = myif(Za21*Pa12, Pa12, -Pa12)
      Td3 = dMb/(1D0+dMb)
      Q=2*MNE1
      Mt = MtEff(Q)
      Mb = MbEff(Q)
      Mc = McEff(Q)
      xH2 = MSsL**2/MSG**2
      yH2 = MSsR**2/MSG**2
      dMd = -3D0/2D0*alphaQCD(MSG)/3.1415D0*tB*mu/MSG*(xH2
     &      *log(xH2)/(1D0-xH2)-yH2*log(yH2)/(1D0-yH2))/(xH2
     &      -yH2)
      Td2 = dMd/(1D0+dMd)
      fiuu = atan(-2D0*Mq*(Au-mu/tB)/(MSuR**2-MSuL**2))/2D0
      fidd = atan(-2D0*Mq*(1D0-Td2)*(Ad-mu*tB)/(MSdR**2-MSdL
     &      **2))/2D0
      ficc = atan(-2D0*Mc*(Au-mu/tB)/(MScR**2-MScL**2))/2D0
      fiss = atan(-2D0*Mq*(1D0-Td2)*(Ad-mu*tB)/(MSsR**2-MSsL
     &      **2))/2D0
      Zuu11 = cos(fiuu)
      Zuu12 = sin(fiuu)
      Zuu21 = -sin(fiuu)
      Zuu22 = cos(fiuu)
      Zdd11 = cos(fidd)
      Zdd12 = sin(fidd)
      Zdd21 = -sin(fidd)
      Zdd22 = cos(fidd)
      Zcc11 = cos(ficc)
      Zcc12 = sin(ficc)
      Zcc21 = -sin(ficc)
      Zcc22 = cos(ficc)
      Zss11 = cos(fiss)
      Zss12 = sin(fiss)
      Zss21 = -sin(fiss)
      Zss22 = cos(fiss)
      MbMM = MbRun(sqrt(MSb1*MSb2))
      MtMM = MtRun(sqrt(MSt1*MSt2))
      NMM11 = Zn11*MNE1*Zn11+Zn21*MNE2*Zn21+Zn31*MNE3*Zn31
     &      +Zn41*MNE4*Zn41+Zn51*MNE5*Zn51
      NMM12 = Zn11*MNE1*Zn12+Zn21*MNE2*Zn22+Zn31*MNE3*Zn32
     &      +Zn41*MNE4*Zn42+Zn51*MNE5*Zn52
      NMM13 = Zn11*MNE1*Zn13+Zn21*MNE2*Zn23+Zn31*MNE3*Zn33
     &      +Zn41*MNE4*Zn43+Zn51*MNE5*Zn53
      NMM14 = Zn11*MNE1*Zn14+Zn21*MNE2*Zn24+Zn31*MNE3*Zn34
     &      +Zn41*MNE4*Zn44+Zn51*MNE5*Zn54
      NMM15 = Zn11*MNE1*Zn15+Zn21*MNE2*Zn25+Zn31*MNE3*Zn35
     &      +Zn41*MNE4*Zn45+Zn51*MNE5*Zn55
      NMM22 = Zn12*MNE1*Zn12+Zn22*MNE2*Zn22+Zn32*MNE3*Zn32
     &      +Zn42*MNE4*Zn42+Zn52*MNE5*Zn52
      NMM23 = Zn12*MNE1*Zn13+Zn22*MNE2*Zn23+Zn32*MNE3*Zn33
     &      +Zn42*MNE4*Zn43+Zn52*MNE5*Zn53
      NMM24 = Zn12*MNE1*Zn14+Zn22*MNE2*Zn24+Zn32*MNE3*Zn34
     &      +Zn42*MNE4*Zn44+Zn52*MNE5*Zn54
      NMM25 = Zn12*MNE1*Zn15+Zn22*MNE2*Zn25+Zn32*MNE3*Zn35
     &      +Zn42*MNE4*Zn45+Zn52*MNE5*Zn55
      NMM33 = Zn13*MNE1*Zn13+Zn23*MNE2*Zn23+Zn33*MNE3*Zn33
     &      +Zn43*MNE4*Zn43+Zn53*MNE5*Zn53
      NMM34 = Zn13*MNE1*Zn14+Zn23*MNE2*Zn24+Zn33*MNE3*Zn34
     &      +Zn43*MNE4*Zn44+Zn53*MNE5*Zn54
      NMM35 = Zn13*MNE1*Zn15+Zn23*MNE2*Zn25+Zn33*MNE3*Zn35
     &      +Zn43*MNE4*Zn45+Zn53*MNE5*Zn55
      NMM44 = Zn14*MNE1*Zn14+Zn24*MNE2*Zn24+Zn34*MNE3*Zn34
     &      +Zn44*MNE4*Zn44+Zn54*MNE5*Zn54
      NMM45 = Zn14*MNE1*Zn15+Zn24*MNE2*Zn25+Zn34*MNE3*Zn35
     &      +Zn44*MNE4*Zn45+Zn54*MNE5*Zn55
      NMM55 = Zn15*MNE1*Zn15+Zn25*MNE2*Zn25+Zn35*MNE3*Zn35
     &      +Zn45*MNE4*Zn45+Zn55*MNE5*Zn55
      MG1I = NMM11
      MG2I = NMM22

      call aaini01
      call aaini02
      call aaini03
      call mtrini
      end

      subroutine aaini01
      implicit none
#include "model.h"

      AAABR(1) = MW*(1D0/CW-CW-SW**2/CW)
      AAABR(2) = CW**2+SW**2
      AAABR(3) = 1D0-CW**2-SW**2
      AAABR(4) = MW**2*(CW**2+2D0*SW**2+SW**4/CW**2)
      AAABR(5) = EE*MW
      AAABR(6) = EE*MW/SW
      AAABR(7) = EE*MW/SW*(Zh21*cb+Zh22*sb)
      AAABR(8) = EE*MW/SW*(Zh31*cb+Zh32*sb)
      AAABR(9) = EE*MW/SW*(Zh11*cb+Zh12*sb)
      AAABR(10) = EE*MW*(SW/CW-CW/SW)
      AAABR(11) = EE/CW*MW/SW
      AAABR(12) = EE/CW**2*MW/SW*(Zh11*cb+Zh12*sb)
      AAABR(13) = EE/CW**2*MW/SW*(Zh31*cb+Zh32*sb)
      AAABR(14) = EE/CW**2*MW/SW*(Zh21*cb+Zh22*sb)
      AAABR(15) = CW*EE/SW
      AAABR(16) = Sqrt2*Zh23*aa2+Sqrt2*Zh23*aa1*sb**2-Sqrt2
     &      *Zh23*aa2*sb**2+Sqrt2*Zh23*la2s*xvev+Sqrt2*Zh23
     &      *aa5*cb*sb+Sqrt2*Zh23*aa6*cb*sb+Sqrt2*Zh23*la1s
     &      *sb**2*xvev-Sqrt2*Zh23*la2s*sb**2*xvev+2D0*MW
     &      /EE*SW*Zh21*cb*la3+4D0*MW/EE*SW*Zh21*la6*sb-6D0
     &      *MW/EE*SW*Zh21*la6*sb**3-2D0*MW/EE*SW*Zh21*la7
     &      *sb+2D0*MW/EE*SW*Zh21*la7*sb**3-2D0*MW/EE*SW
     &      *Zh22*cb*la7+2D0*MW/EE*SW*Zh22*la2*sb-2D0*MW
     &      /EE*SW*Zh22*la2*sb**3+2D0*MW/EE*SW*Zh22*la3*sb
     &      **3-2D0*MW/EE*SW*Zh22*la4*sb+2D0*MW/EE*SW*Zh22
     &      *la4*sb**3-2D0*MW/EE*SW*Zh22*la5*sb+2D0*MW/EE
     &      *SW*Zh22*la5*sb**3+2D0*Sqrt2*Zh23*cb*la4s*sb
     &      *xvev+2D0*Sqrt2*Zh23*cb*la5s*sb*xvev+2D0*Sqrt2
     &      *Zh23*cb*la6s*sb*xvev+2D0*MW/EE*SW*Zh21*cb*la1
     &      *sb**2-2D0*MW/EE*SW*Zh21*cb*la3*sb**2-2D0*MW
     &      /EE*SW*Zh21*cb*la4*sb**2-2D0*MW/EE*SW*Zh21*cb
     &      *la5*sb**2-2D0*MW/EE*SW*Zh22*cb*la6*sb**2+6D0
     &      *MW/EE*SW*Zh22*cb*la7*sb**2
      AAABR(17) = Sqrt2*Zh33*aa2+Sqrt2*Zh33*aa1*sb**2-Sqrt2
     &      *Zh33*aa2*sb**2+Sqrt2*Zh33*la2s*xvev+Sqrt2*Zh33
     &      *aa5*cb*sb+Sqrt2*Zh33*aa6*cb*sb+Sqrt2*Zh33*la1s
     &      *sb**2*xvev-Sqrt2*Zh33*la2s*sb**2*xvev+2D0*MW
     &      /EE*SW*Zh31*cb*la3+4D0*MW/EE*SW*Zh31*la6*sb-6D0
     &      *MW/EE*SW*Zh31*la6*sb**3-2D0*MW/EE*SW*Zh31*la7
     &      *sb+2D0*MW/EE*SW*Zh31*la7*sb**3-2D0*MW/EE*SW
     &      *Zh32*cb*la7+2D0*MW/EE*SW*Zh32*la2*sb-2D0*MW
     &      /EE*SW*Zh32*la2*sb**3+2D0*MW/EE*SW*Zh32*la3*sb
     &      **3-2D0*MW/EE*SW*Zh32*la4*sb+2D0*MW/EE*SW*Zh32
     &      *la4*sb**3-2D0*MW/EE*SW*Zh32*la5*sb+2D0*MW/EE
     &      *SW*Zh32*la5*sb**3+2D0*Sqrt2*Zh33*cb*la4s*sb
     &      *xvev+2D0*Sqrt2*Zh33*cb*la5s*sb*xvev+2D0*Sqrt2
     &      *Zh33*cb*la6s*sb*xvev+2D0*MW/EE*SW*Zh31*cb*la1
     &      *sb**2-2D0*MW/EE*SW*Zh31*cb*la3*sb**2-2D0*MW
     &      /EE*SW*Zh31*cb*la4*sb**2-2D0*MW/EE*SW*Zh31*cb
     &      *la5*sb**2-2D0*MW/EE*SW*Zh32*cb*la6*sb**2+6D0
     &      *MW/EE*SW*Zh32*cb*la7*sb**2
      AAABR(18) = Sqrt2*Zh13*aa2+Sqrt2*Zh13*aa1*sb**2-Sqrt2
     &      *Zh13*aa2*sb**2+Sqrt2*Zh13*la2s*xvev+Sqrt2*Zh13
     &      *aa5*cb*sb+Sqrt2*Zh13*aa6*cb*sb+Sqrt2*Zh13*la1s
     &      *sb**2*xvev-Sqrt2*Zh13*la2s*sb**2*xvev+2D0*MW
     &      /EE*SW*Zh11*cb*la3+4D0*MW/EE*SW*Zh11*la6*sb-6D0
     &      *MW/EE*SW*Zh11*la6*sb**3-2D0*MW/EE*SW*Zh11*la7
     &      *sb+2D0*MW/EE*SW*Zh11*la7*sb**3-2D0*MW/EE*SW
     &      *Zh12*cb*la7+2D0*MW/EE*SW*Zh12*la2*sb-2D0*MW
     &      /EE*SW*Zh12*la2*sb**3+2D0*MW/EE*SW*Zh12*la3*sb
     &      **3-2D0*MW/EE*SW*Zh12*la4*sb+2D0*MW/EE*SW*Zh12
     &      *la4*sb**3-2D0*MW/EE*SW*Zh12*la5*sb+2D0*MW/EE
     &      *SW*Zh12*la5*sb**3+2D0*Sqrt2*Zh13*cb*la4s*sb
     &      *xvev+2D0*Sqrt2*Zh13*cb*la5s*sb*xvev+2D0*Sqrt2
     &      *Zh13*cb*la6s*sb*xvev+2D0*MW/EE*SW*Zh11*cb*la1
     &      *sb**2-2D0*MW/EE*SW*Zh11*cb*la3*sb**2-2D0*MW
     &      /EE*SW*Zh11*cb*la4*sb**2-2D0*MW/EE*SW*Zh11*cb
     &      *la5*sb**2-2D0*MW/EE*SW*Zh12*cb*la6*sb**2+6D0
     &      *MW/EE*SW*Zh12*cb*la7*sb**2
      AAABR(19) = EE*MW/SW*Sqrt2*cb*sb
      AAABR(20) = EE/SW*Sqrt2*(2D0*MW*Zl11/cb*sb-2D0*MW*Zl11
     &      /cb*sb**3-Lambda/MW*Ml*Zl12*xvev-Ml**2/MW*Zl11
     &      /cb*sb-Al/MW*Ml*Zl12/cb*sb)
      AAABR(21) = EE/SW*Sqrt2*(2D0*MW*Zl21/cb*sb-2D0*MW*Zl21
     &      /cb*sb**3-Lambda/MW*Ml*Zl22*xvev-Ml**2/MW*Zl21
     &      /cb*sb-Al/MW*Ml*Zl22/cb*sb)
      AAABR(22) = EE/SW*Sqrt2*(2D0*MW*Zdd11*Zuu11/cb*sb-2D0
     &      *MW*Zdd11*Zuu11/cb*sb**3-Lambda/MW*Mq*Zdd11*Zuu12
     &      *xvev-Lambda/MW*Mq*Zdd12*Zuu11*xvev-Mq**2/MW
     &      *Zdd11*Zuu11/cb/sb-Mq**2/MW*Zdd12*Zuu12/cb/sb
     &      -Ad/MW*Mq*Zdd12*Zuu11/cb*sb-Au/MW*Mq*Zdd11*Zuu12
     &      /cb/sb+Au/MW*Mq*Zdd11*Zuu12/cb*sb)
      AAABR(23) = EE/SW*Sqrt2*(2D0*MW*Zb11*Zt11/cb*sb-2D0
     &      *MW*Zb11*Zt11/cb*sb**3-Lambda/MW*MbMM*Zb12*Zt11
     &      *xvev-Lambda/MW*MtMM*Zb11*Zt12*xvev-MbMM**2/MW
     &      *Zb11*Zt11/cb*sb-MtMM**2/MW*Zb11*Zt11/cb/sb+MtMM
     &      **2/MW*Zb11*Zt11/cb*sb-Ab/MW*MbMM*Zb12*Zt11/cb
     &      *sb-At/MW*MtMM*Zb11*Zt12/cb/sb+At/MW*MtMM*Zb11
     &      *Zt12/cb*sb-MbMM/MW*MtMM*Zb12*Zt12/cb/sb)
      AAABR(24) = EE/SW*Sqrt2*(2D0*MW*Zcc11*Zss11/cb*sb-2D0
     &      *MW*Zcc11*Zss11/cb*sb**3-Lambda/MW*Mc*Zcc12*Zss11
     &      *xvev-Lambda/MW*Mq*Zcc11*Zss12*xvev-Mc**2/MW
     &      *Zcc11*Zss11/cb/sb+Mc**2/MW*Zcc11*Zss11/cb*sb
     &      -Mq**2/MW*Zcc11*Zss11/cb*sb-Ad/MW*Mq*Zcc11*Zss12
     &      /cb*sb-Au/MW*Mc*Zcc12*Zss11/cb/sb+Au/MW*Mc*Zcc12
     &      *Zss11/cb*sb-Mc/MW*Mq*Zcc12*Zss12/cb/sb)
      AAABR(25) = EE/SW*Sqrt2*(2D0*MW*Zcc11*Zss21/cb*sb-2D0
     &      *MW*Zcc11*Zss21/cb*sb**3-Lambda/MW*Mc*Zcc12*Zss21
     &      *xvev-Lambda/MW*Mq*Zcc11*Zss22*xvev-Mc**2/MW
     &      *Zcc11*Zss21/cb/sb+Mc**2/MW*Zcc11*Zss21/cb*sb
     &      -Mq**2/MW*Zcc11*Zss21/cb*sb-Ad/MW*Mq*Zcc11*Zss22
     &      /cb*sb-Au/MW*Mc*Zcc12*Zss21/cb/sb+Au/MW*Mc*Zcc12
     &      *Zss21/cb*sb-Mc/MW*Mq*Zcc12*Zss22/cb/sb)
      AAABR(26) = EE/SW*Sqrt2*(2D0*MW*Zb21*Zt11/cb*sb-2D0
     &      *MW*Zb21*Zt11/cb*sb**3-Lambda/MW*MbMM*Zb22*Zt11
     &      *xvev-Lambda/MW*MtMM*Zb21*Zt12*xvev-MbMM**2/MW
     &      *Zb21*Zt11/cb*sb-MtMM**2/MW*Zb21*Zt11/cb/sb+MtMM
     &      **2/MW*Zb21*Zt11/cb*sb-Ab/MW*MbMM*Zb22*Zt11/cb
     &      *sb-At/MW*MtMM*Zb21*Zt12/cb/sb+At/MW*MtMM*Zb21
     &      *Zt12/cb*sb-MbMM/MW*MtMM*Zb22*Zt12/cb/sb)
      AAABR(27) = EE/SW*Sqrt2*(2D0*MW*Zdd21*Zuu11/cb*sb-2D0
     &      *MW*Zdd21*Zuu11/cb*sb**3-Lambda/MW*Mq*Zdd21*Zuu12
     &      *xvev-Lambda/MW*Mq*Zdd22*Zuu11*xvev-Mq**2/MW
     &      *Zdd21*Zuu11/cb/sb-Mq**2/MW*Zdd22*Zuu12/cb/sb
     &      -Ad/MW*Mq*Zdd22*Zuu11/cb*sb-Au/MW*Mq*Zdd21*Zuu12
     &      /cb/sb+Au/MW*Mq*Zdd21*Zuu12/cb*sb)
      AAABR(28) = EE/SW*Sqrt2*(2D0*MW*Zb11*Zt21/cb*sb-2D0
     &      *MW*Zb11*Zt21/cb*sb**3-Lambda/MW*MbMM*Zb12*Zt21
     &      *xvev-Lambda/MW*MtMM*Zb11*Zt22*xvev-MbMM**2/MW
     &      *Zb11*Zt21/cb*sb-MtMM**2/MW*Zb11*Zt21/cb/sb+MtMM
     &      **2/MW*Zb11*Zt21/cb*sb-Ab/MW*MbMM*Zb12*Zt21/cb
     &      *sb-At/MW*MtMM*Zb11*Zt22/cb/sb+At/MW*MtMM*Zb11
     &      *Zt22/cb*sb-MbMM/MW*MtMM*Zb12*Zt22/cb/sb)
      AAABR(29) = EE/SW*Sqrt2*(2D0*MW*Zcc21*Zss11/cb*sb-2D0
     &      *MW*Zcc21*Zss11/cb*sb**3-Lambda/MW*Mc*Zcc22*Zss11
     &      *xvev-Lambda/MW*Mq*Zcc21*Zss12*xvev-Mc**2/MW
     &      *Zcc21*Zss11/cb/sb+Mc**2/MW*Zcc21*Zss11/cb*sb
     &      -Mq**2/MW*Zcc21*Zss11/cb*sb-Ad/MW*Mq*Zcc21*Zss12
     &      /cb*sb-Au/MW*Mc*Zcc22*Zss11/cb/sb+Au/MW*Mc*Zcc22
     &      *Zss11/cb*sb-Mc/MW*Mq*Zcc22*Zss12/cb/sb)
      AAABR(30) = EE/SW*Sqrt2*(2D0*MW*Zdd11*Zuu21/cb*sb-2D0
     &      *MW*Zdd11*Zuu21/cb*sb**3-Lambda/MW*Mq*Zdd11*Zuu22
     &      *xvev-Lambda/MW*Mq*Zdd12*Zuu21*xvev-Mq**2/MW
     &      *Zdd11*Zuu21/cb/sb-Mq**2/MW*Zdd12*Zuu22/cb/sb
     &      -Ad/MW*Mq*Zdd12*Zuu21/cb*sb-Au/MW*Mq*Zdd11*Zuu22
     &      /cb/sb+Au/MW*Mq*Zdd11*Zuu22/cb*sb)
      AAABR(31) = EE/SW*Sqrt2*(2D0*MW*Zb21*Zt21/cb*sb-2D0
     &      *MW*Zb21*Zt21/cb*sb**3-Lambda/MW*MbMM*Zb22*Zt21
     &      *xvev-Lambda/MW*MtMM*Zb21*Zt22*xvev-MbMM**2/MW
     &      *Zb21*Zt21/cb*sb-MtMM**2/MW*Zb21*Zt21/cb/sb+MtMM
     &      **2/MW*Zb21*Zt21/cb*sb-Ab/MW*MbMM*Zb22*Zt21/cb
     &      *sb-At/MW*MtMM*Zb21*Zt22/cb/sb+At/MW*MtMM*Zb21
     &      *Zt22/cb*sb-MbMM/MW*MtMM*Zb22*Zt22/cb/sb)
      AAABR(32) = EE/SW*Sqrt2*(2D0*MW*Zcc21*Zss21/cb*sb-2D0
     &      *MW*Zcc21*Zss21/cb*sb**3-Lambda/MW*Mc*Zcc22*Zss21
     &      *xvev-Lambda/MW*Mq*Zcc21*Zss22*xvev-Mc**2/MW
     &      *Zcc21*Zss21/cb/sb+Mc**2/MW*Zcc21*Zss21/cb*sb
     &      -Mq**2/MW*Zcc21*Zss21/cb*sb-Ad/MW*Mq*Zcc21*Zss22
     &      /cb*sb-Au/MW*Mc*Zcc22*Zss21/cb/sb+Au/MW*Mc*Zcc22
     &      *Zss21/cb*sb-Mc/MW*Mq*Zcc22*Zss22/cb/sb)
      AAABR(33) = EE/SW*Sqrt2*(2D0*MW*Zdd21*Zuu21/cb*sb-2D0
     &      *MW*Zdd21*Zuu21/cb*sb**3-Lambda/MW*Mq*Zdd21*Zuu22
     &      *xvev-Lambda/MW*Mq*Zdd22*Zuu21*xvev-Mq**2/MW
     &      *Zdd21*Zuu21/cb/sb-Mq**2/MW*Zdd22*Zuu22/cb/sb
     &      -Ad/MW*Mq*Zdd22*Zuu21/cb*sb-Au/MW*Mq*Zdd21*Zuu22
     &      /cb/sb+Au/MW*Mq*Zdd21*Zuu22/cb*sb)
      AAABR(34) = Pa12*Sqrt2*aa5-Pa12*Sqrt2*aa6+2D0*Pa12
     &      *Sqrt2*la5s*xvev-2D0*Pa12*Sqrt2*la6s*xvev+2D0
     &      *MW/EE*Pa11*SW*la4-2D0*MW/EE*Pa11*SW*la5
      AAABR(35) = Pa22*Sqrt2*aa5-Pa22*Sqrt2*aa6+2D0*Pa22
     &      *Sqrt2*la5s*xvev-2D0*Pa22*Sqrt2*la6s*xvev+2D0
     &      *MW/EE*Pa21*SW*la4-2D0*MW/EE*Pa21*SW*la5
      AAABR(36) = Sqrt2*Zh13*aa5+Sqrt2*Zh13*aa6-2D0*Sqrt2
     &      *Zh13*aa5*sb**2-2D0*Sqrt2*Zh13*aa6*sb**2+2D0
     &      *Sqrt2*Zh13*la4s*xvev+2D0*Sqrt2*Zh13*la5s*xvev
     &      +2D0*Sqrt2*Zh13*la6s*xvev+2D0*Sqrt2*Zh13*aa1
     &      *cb*sb-2D0*Sqrt2*Zh13*aa2*cb*sb-4D0*Sqrt2*Zh13
     &      *la4s*sb**2*xvev-4D0*Sqrt2*Zh13*la5s*sb**2*xvev
     &      -4D0*Sqrt2*Zh13*la6s*sb**2*xvev+4D0*MW/EE*SW
     &      *Zh11*cb*la6+4D0*MW/EE*SW*Zh11*la1*sb-4D0*MW
     &      /EE*SW*Zh11*la1*sb**3-4D0*MW/EE*SW*Zh11*la3*sb
     &      +4D0*MW/EE*SW*Zh11*la3*sb**3-2D0*MW/EE*SW*Zh11
     &      *la4*sb+4D0*MW/EE*SW*Zh11*la4*sb**3-2D0*MW/EE
     &      *SW*Zh11*la5*sb+4D0*MW/EE*SW*Zh11*la5*sb**3-2D0
     &      *MW/EE*SW*Zh12*cb*la4-2D0*MW/EE*SW*Zh12*cb*la5
     &      -4D0*MW/EE*SW*Zh12*la6*sb+4D0*MW/EE*SW*Zh12*la6
     &      *sb**3+8D0*MW/EE*SW*Zh12*la7*sb-12D0*MW/EE*SW
     &      *Zh12*la7*sb**3+2D0*Sqrt2*Zh13*cb*la1s*sb*xvev
     &      -2D0*Sqrt2*Zh13*cb*la2s*sb*xvev-12D0*MW/EE*SW
     &      *Zh11*cb*la6*sb**2+4D0*MW/EE*SW*Zh11*cb*la7*sb
     &      **2-4D0*MW/EE*SW*Zh12*cb*la2*sb**2+4D0*MW/EE
     &      *SW*Zh12*cb*la3*sb**2+4D0*MW/EE*SW*Zh12*cb*la4
     &      *sb**2+4D0*MW/EE*SW*Zh12*cb*la5*sb**2
      AAABR(37) = Sqrt2*Zh23*aa5+Sqrt2*Zh23*aa6-2D0*Sqrt2
     &      *Zh23*aa5*sb**2-2D0*Sqrt2*Zh23*aa6*sb**2+2D0
     &      *Sqrt2*Zh23*la4s*xvev+2D0*Sqrt2*Zh23*la5s*xvev
     &      +2D0*Sqrt2*Zh23*la6s*xvev+2D0*Sqrt2*Zh23*aa1
     &      *cb*sb-2D0*Sqrt2*Zh23*aa2*cb*sb-4D0*Sqrt2*Zh23
     &      *la4s*sb**2*xvev-4D0*Sqrt2*Zh23*la5s*sb**2*xvev
     &      -4D0*Sqrt2*Zh23*la6s*sb**2*xvev+4D0*MW/EE*SW
     &      *Zh21*cb*la6+4D0*MW/EE*SW*Zh21*la1*sb-4D0*MW
     &      /EE*SW*Zh21*la1*sb**3-4D0*MW/EE*SW*Zh21*la3*sb
     &      +4D0*MW/EE*SW*Zh21*la3*sb**3-2D0*MW/EE*SW*Zh21
     &      *la4*sb+4D0*MW/EE*SW*Zh21*la4*sb**3-2D0*MW/EE
     &      *SW*Zh21*la5*sb+4D0*MW/EE*SW*Zh21*la5*sb**3-2D0
     &      *MW/EE*SW*Zh22*cb*la4-2D0*MW/EE*SW*Zh22*cb*la5
     &      -4D0*MW/EE*SW*Zh22*la6*sb+4D0*MW/EE*SW*Zh22*la6
     &      *sb**3+8D0*MW/EE*SW*Zh22*la7*sb-12D0*MW/EE*SW
     &      *Zh22*la7*sb**3+2D0*Sqrt2*Zh23*cb*la1s*sb*xvev
     &      -2D0*Sqrt2*Zh23*cb*la2s*sb*xvev-12D0*MW/EE*SW
     &      *Zh21*cb*la6*sb**2+4D0*MW/EE*SW*Zh21*cb*la7*sb
     &      **2-4D0*MW/EE*SW*Zh22*cb*la2*sb**2+4D0*MW/EE
     &      *SW*Zh22*cb*la3*sb**2+4D0*MW/EE*SW*Zh22*cb*la4
     &      *sb**2+4D0*MW/EE*SW*Zh22*cb*la5*sb**2
      AAABR(38) = Sqrt2*Zh33*aa5+Sqrt2*Zh33*aa6-2D0*Sqrt2
     &      *Zh33*aa5*sb**2-2D0*Sqrt2*Zh33*aa6*sb**2+2D0
     &      *Sqrt2*Zh33*la4s*xvev+2D0*Sqrt2*Zh33*la5s*xvev
     &      +2D0*Sqrt2*Zh33*la6s*xvev+2D0*Sqrt2*Zh33*aa1
     &      *cb*sb-2D0*Sqrt2*Zh33*aa2*cb*sb-4D0*Sqrt2*Zh33
     &      *la4s*sb**2*xvev-4D0*Sqrt2*Zh33*la5s*sb**2*xvev
     &      -4D0*Sqrt2*Zh33*la6s*sb**2*xvev+4D0*MW/EE*SW
     &      *Zh31*cb*la6+4D0*MW/EE*SW*Zh31*la1*sb-4D0*MW
     &      /EE*SW*Zh31*la1*sb**3-4D0*MW/EE*SW*Zh31*la3*sb
     &      +4D0*MW/EE*SW*Zh31*la3*sb**3-2D0*MW/EE*SW*Zh31
     &      *la4*sb+4D0*MW/EE*SW*Zh31*la4*sb**3-2D0*MW/EE
     &      *SW*Zh31*la5*sb+4D0*MW/EE*SW*Zh31*la5*sb**3-2D0
     &      *MW/EE*SW*Zh32*cb*la4-2D0*MW/EE*SW*Zh32*cb*la5
     &      -4D0*MW/EE*SW*Zh32*la6*sb+4D0*MW/EE*SW*Zh32*la6
     &      *sb**3+8D0*MW/EE*SW*Zh32*la7*sb-12D0*MW/EE*SW
     &      *Zh32*la7*sb**3+2D0*Sqrt2*Zh33*cb*la1s*sb*xvev
     &      -2D0*Sqrt2*Zh33*cb*la2s*sb*xvev-12D0*MW/EE*SW
     &      *Zh31*cb*la6*sb**2+4D0*MW/EE*SW*Zh31*cb*la7*sb
     &      **2-4D0*MW/EE*SW*Zh32*cb*la2*sb**2+4D0*MW/EE
     &      *SW*Zh32*cb*la3*sb**2+4D0*MW/EE*SW*Zh32*cb*la4
     &      *sb**2+4D0*MW/EE*SW*Zh32*cb*la5*sb**2
      AAABR(39) = EE*MW/SW*Sqrt2*(1D0-2D0*sb**2)
      AAABR(40) = EE/SW*Sqrt2*(MW*Zl11-Ml**2/MW*Zl11-2D0
     &      *MW*Zl11*sb**2-Al/MW*Ml*Zl12+Lambda/MW*Ml*Zl12
     &      /cb*sb*xvev)
      AAABR(41) = EE/MW*Ml/SW*(Al*Zl11*Zl22-Al*Zl12*Zl21
     &      -Lambda*Zl11*Zl22/cb*sb*xvev+Lambda*Zl12*Zl21
     &      /cb*sb*xvev)
      AAABR(42) = Ml*(Lambda*Pa12*Sqrt2*Zl11*Zl22/cb*sb-Lambda
     &      *Pa12*Sqrt2*Zl12*Zl21/cb*sb+EE*Lambda/MW*Pa11
     &      /SW*Zl11*Zl22*xvev-EE*Lambda/MW*Pa11/SW*Zl12
     &      *Zl21*xvev+Al*EE/MW*Pa11/SW*Zl11*Zl22/cb*sb-Al
     &      *EE/MW*Pa11/SW*Zl12*Zl21/cb*sb)
      AAABR(43) = Ml*(Lambda*Pa22*Sqrt2*Zl11*Zl22/cb*sb-Lambda
     &      *Pa22*Sqrt2*Zl12*Zl21/cb*sb+EE*Lambda/MW*Pa21
     &      /SW*Zl11*Zl22*xvev-EE*Lambda/MW*Pa21/SW*Zl12
     &      *Zl21*xvev+Al*EE/MW*Pa21/SW*Zl11*Zl22/cb*sb-Al
     &      *EE/MW*Pa21/SW*Zl12*Zl21/cb*sb)
      AAABR(44) = EE*MW*(Zh21/SW*cb-Zh22/SW*sb-SW/CW**2*Zh21
     &      *cb+SW/CW**2*Zh22*sb)
      AAABR(45) = 2D0*EE/MW*Ml**2/SW*Zh11/cb-EE*MW/SW*Zh11
     &      *Zl11**2/cb+EE*MW/SW*Zh12*Zl11**2*sb+EE/CW**2
     &      *MW*SW*Zh11*Zl11**2/cb-2D0*EE/CW**2*MW*SW*Zh11
     &      *Zl12**2/cb-EE/CW**2*MW*SW*Zh12*Zl11**2*sb+2D0
     &      *EE/CW**2*MW*SW*Zh12*Zl12**2*sb+EE*MW/SW*Zh11
     &      *Zl11**2/cb*sb**2-EE/CW**2*MW*SW*Zh11*Zl11**2
     &      /cb*sb**2+2D0*EE/CW**2*MW*SW*Zh11*Zl12**2/cb
     &      *sb**2-2D0*Lambda*Ml*Sqrt2*Zh13*Zl11*Zl12/cb
     &      *sb+2D0*Al*EE/MW*Ml/SW*Zh11*Zl11*Zl12/cb-2D0
     &      *EE*Lambda/MW*Ml/SW*Zh12*Zl11*Zl12/cb*xvev
      AAABR(46) = EE*MW*(Zh11/SW*cb-Zh12/SW*sb-SW/CW**2*Zh11
     &      *cb+SW/CW**2*Zh12*sb)
      AAABR(47) = 2D0*EE/MW*Ml**2/SW*Zh21/cb-EE*MW/SW*Zh21
     &      *Zl11**2/cb+EE*MW/SW*Zh22*Zl11**2*sb+EE/CW**2
     &      *MW*SW*Zh21*Zl11**2/cb-2D0*EE/CW**2*MW*SW*Zh21
     &      *Zl12**2/cb-EE/CW**2*MW*SW*Zh22*Zl11**2*sb+2D0
     &      *EE/CW**2*MW*SW*Zh22*Zl12**2*sb+EE*MW/SW*Zh21
     &      *Zl11**2/cb*sb**2-EE/CW**2*MW*SW*Zh21*Zl11**2
     &      /cb*sb**2+2D0*EE/CW**2*MW*SW*Zh21*Zl12**2/cb
     &      *sb**2-2D0*Lambda*Ml*Sqrt2*Zh23*Zl11*Zl12/cb
     &      *sb+2D0*Al*EE/MW*Ml/SW*Zh21*Zl11*Zl12/cb-2D0
     &      *EE*Lambda/MW*Ml/SW*Zh22*Zl11*Zl12/cb*xvev
      AAABR(48) = EE*MW*(Zh31/SW*cb-Zh32/SW*sb-SW/CW**2*Zh31
     &      *cb+SW/CW**2*Zh32*sb)
      AAABR(49) = 2D0*EE/MW*Ml**2/SW*Zh31/cb-EE*MW/SW*Zh31
     &      *Zl11**2/cb+EE*MW/SW*Zh32*Zl11**2*sb+EE/CW**2
     &      *MW*SW*Zh31*Zl11**2/cb-2D0*EE/CW**2*MW*SW*Zh31
     &      *Zl12**2/cb-EE/CW**2*MW*SW*Zh32*Zl11**2*sb+2D0
     &      *EE/CW**2*MW*SW*Zh32*Zl12**2*sb+EE*MW/SW*Zh31
     &      *Zl11**2/cb*sb**2-EE/CW**2*MW*SW*Zh31*Zl11**2
     &      /cb*sb**2+2D0*EE/CW**2*MW*SW*Zh31*Zl12**2/cb
     &      *sb**2-2D0*Lambda*Ml*Sqrt2*Zh33*Zl11*Zl12/cb
     &      *sb+2D0*Al*EE/MW*Ml/SW*Zh31*Zl11*Zl12/cb-2D0
     &      *EE*Lambda/MW*Ml/SW*Zh32*Zl11*Zl12/cb*xvev
      AAABR(50) = EE*MW/SW*Zh11*Zl11*Zl21/cb-EE*MW/SW*Zh12
     &      *Zl11*Zl21*sb-EE/CW**2*MW*SW*Zh11*Zl11*Zl21/cb
     &      +2D0*EE/CW**2*MW*SW*Zh11*Zl12*Zl22/cb+EE/CW**2
     &      *MW*SW*Zh12*Zl11*Zl21*sb-2D0*EE/CW**2*MW*SW*Zh12
     &      *Zl12*Zl22*sb-EE*MW/SW*Zh11*Zl11*Zl21/cb*sb**2
     &      +Lambda*Ml*Sqrt2*Zh13*Zl11*Zl22/cb*sb+Lambda
     &      *Ml*Sqrt2*Zh13*Zl12*Zl21/cb*sb-Al*EE/MW*Ml/SW
     &      *Zh11*Zl11*Zl22/cb-Al*EE/MW*Ml/SW*Zh11*Zl12*Zl21
     &      /cb+EE/CW**2*MW*SW*Zh11*Zl11*Zl21/cb*sb**2-2D0
     &      *EE/CW**2*MW*SW*Zh11*Zl12*Zl22/cb*sb**2+EE*Lambda
     &      /MW*Ml/SW*Zh12*Zl11*Zl22/cb*xvev+EE*Lambda/MW
     &      *Ml/SW*Zh12*Zl12*Zl21/cb*xvev
      AAABR(51) = EE*MW/SW*Zh31*Zl11*Zl21/cb-EE*MW/SW*Zh32
     &      *Zl11*Zl21*sb-EE/CW**2*MW*SW*Zh31*Zl11*Zl21/cb
     &      +2D0*EE/CW**2*MW*SW*Zh31*Zl12*Zl22/cb+EE/CW**2
     &      *MW*SW*Zh32*Zl11*Zl21*sb-2D0*EE/CW**2*MW*SW*Zh32
     &      *Zl12*Zl22*sb-EE*MW/SW*Zh31*Zl11*Zl21/cb*sb**2
     &      +Lambda*Ml*Sqrt2*Zh33*Zl11*Zl22/cb*sb+Lambda
     &      *Ml*Sqrt2*Zh33*Zl12*Zl21/cb*sb-Al*EE/MW*Ml/SW
     &      *Zh31*Zl11*Zl22/cb-Al*EE/MW*Ml/SW*Zh31*Zl12*Zl21
     &      /cb+EE/CW**2*MW*SW*Zh31*Zl11*Zl21/cb*sb**2-2D0
     &      *EE/CW**2*MW*SW*Zh31*Zl12*Zl22/cb*sb**2+EE*Lambda
     &      /MW*Ml/SW*Zh32*Zl11*Zl22/cb*xvev+EE*Lambda/MW
     &      *Ml/SW*Zh32*Zl12*Zl21/cb*xvev
      AAABR(52) = EE*MW/SW*Zh21*Zl11*Zl21/cb-EE*MW/SW*Zh22
     &      *Zl11*Zl21*sb-EE/CW**2*MW*SW*Zh21*Zl11*Zl21/cb
     &      +2D0*EE/CW**2*MW*SW*Zh21*Zl12*Zl22/cb+EE/CW**2
     &      *MW*SW*Zh22*Zl11*Zl21*sb-2D0*EE/CW**2*MW*SW*Zh22
     &      *Zl12*Zl22*sb-EE*MW/SW*Zh21*Zl11*Zl21/cb*sb**2
     &      +Lambda*Ml*Sqrt2*Zh23*Zl11*Zl22/cb*sb+Lambda
     &      *Ml*Sqrt2*Zh23*Zl12*Zl21/cb*sb-Al*EE/MW*Ml/SW
     &      *Zh21*Zl11*Zl22/cb-Al*EE/MW*Ml/SW*Zh21*Zl12*Zl21
     &      /cb+EE/CW**2*MW*SW*Zh21*Zl11*Zl21/cb*sb**2-2D0
     &      *EE/CW**2*MW*SW*Zh21*Zl12*Zl22/cb*sb**2+EE*Lambda
     &      /MW*Ml/SW*Zh22*Zl11*Zl22/cb*xvev+EE*Lambda/MW
     &      *Ml/SW*Zh22*Zl12*Zl21/cb*xvev
      AAABR(53) = EE/SW*Sqrt2*(MW*Zl21-Ml**2/MW*Zl21-2D0
     &      *MW*Zl21*sb**2-Al/MW*Ml*Zl22+Lambda/MW*Ml*Zl22
     &      /cb*sb*xvev)
      AAABR(54) = EE/CW**2*MW*SW*(Zh31*cb-Zh32*sb)
      AAABR(55) = 2D0*EE/MW*Ml**2/SW*Zh21/cb-EE*MW/SW*Zh21
     &      *Zl21**2/cb+EE*MW/SW*Zh22*Zl21**2*sb+EE/CW**2
     &      *MW*SW*Zh21*Zl21**2/cb-2D0*EE/CW**2*MW*SW*Zh21
     &      *Zl22**2/cb-EE/CW**2*MW*SW*Zh22*Zl21**2*sb+2D0
     &      *EE/CW**2*MW*SW*Zh22*Zl22**2*sb+EE*MW/SW*Zh21
     &      *Zl21**2/cb*sb**2-EE/CW**2*MW*SW*Zh21*Zl21**2
     &      /cb*sb**2+2D0*EE/CW**2*MW*SW*Zh21*Zl22**2/cb
     &      *sb**2-2D0*Lambda*Ml*Sqrt2*Zh23*Zl21*Zl22/cb
     &      *sb+2D0*Al*EE/MW*Ml/SW*Zh21*Zl21*Zl22/cb-2D0
     &      *EE*Lambda/MW*Ml/SW*Zh22*Zl21*Zl22/cb*xvev
      AAABR(56) = EE/CW**2*MW*SW*(Zh11*cb-Zh12*sb)
      AAABR(57) = 2D0*EE/MW*Ml**2/SW*Zh31/cb-EE*MW/SW*Zh31
     &      *Zl21**2/cb+EE*MW/SW*Zh32*Zl21**2*sb+EE/CW**2
     &      *MW*SW*Zh31*Zl21**2/cb-2D0*EE/CW**2*MW*SW*Zh31
     &      *Zl22**2/cb-EE/CW**2*MW*SW*Zh32*Zl21**2*sb+2D0
     &      *EE/CW**2*MW*SW*Zh32*Zl22**2*sb+EE*MW/SW*Zh31
     &      *Zl21**2/cb*sb**2-EE/CW**2*MW*SW*Zh31*Zl21**2
     &      /cb*sb**2+2D0*EE/CW**2*MW*SW*Zh31*Zl22**2/cb
     &      *sb**2-2D0*Lambda*Ml*Sqrt2*Zh33*Zl21*Zl22/cb
     &      *sb+2D0*Al*EE/MW*Ml/SW*Zh31*Zl21*Zl22/cb-2D0
     &      *EE*Lambda/MW*Ml/SW*Zh32*Zl21*Zl22/cb*xvev
      AAABR(58) = EE/CW**2*MW*SW*(Zh21*cb-Zh22*sb)
      AAABR(59) = 2D0*EE/MW*Ml**2/SW*Zh11/cb-EE*MW/SW*Zh11
     &      *Zl21**2/cb+EE*MW/SW*Zh12*Zl21**2*sb+EE/CW**2
     &      *MW*SW*Zh11*Zl21**2/cb-2D0*EE/CW**2*MW*SW*Zh11
     &      *Zl22**2/cb-EE/CW**2*MW*SW*Zh12*Zl21**2*sb+2D0
     &      *EE/CW**2*MW*SW*Zh12*Zl22**2*sb+EE*MW/SW*Zh11
     &      *Zl21**2/cb*sb**2-EE/CW**2*MW*SW*Zh11*Zl21**2
     &      /cb*sb**2+2D0*EE/CW**2*MW*SW*Zh11*Zl22**2/cb
     &      *sb**2-2D0*Lambda*Ml*Sqrt2*Zh13*Zl21*Zl22/cb
     &      *sb+2D0*Al*EE/MW*Ml/SW*Zh11*Zl21*Zl22/cb-2D0
     &      *EE*Lambda/MW*Ml/SW*Zh12*Zl21*Zl22/cb*xvev
      AAABR(60) = EE*MW*(Zh31/SW*cb-Zh32/SW*sb+SW/CW**2*Zh31
     &      *cb-SW/CW**2*Zh32*sb)
      AAABR(61) = EE*MW*(Zh11/SW*cb-Zh12/SW*sb+SW/CW**2*Zh11
     &      *cb-SW/CW**2*Zh12*sb)
      AAABR(62) = EE*MW*(Zh21/SW*cb-Zh22/SW*sb+SW/CW**2*Zh21
     &      *cb-SW/CW**2*Zh22*sb)
      AAABR(63) = EE/SW*Sqrt2*(MW*Zcc11*Zss11+Mc**2/MW*Zcc11
     &      *Zss11-Mq**2/MW*Zcc11*Zss11-2D0*MW*Zcc11*Zss11
     &      *sb**2-Ad/MW*Mq*Zcc11*Zss12+Au/MW*Mc*Zcc12*Zss11
     &      -Lambda/MW*Mc*Zcc12*Zss11/cb/sb*xvev+Lambda/MW
     &      *Mc*Zcc12*Zss11/cb*sb*xvev+Lambda/MW*Mq*Zcc11
     &      *Zss12/cb*sb*xvev)
      AAABR(64) = EE/SW*Sqrt2*(MW*Zdd11*Zuu11-2D0*MW*Zdd11
     &      *Zuu11*sb**2-Ad/MW*Mq*Zdd12*Zuu11+Au/MW*Mq*Zdd11
     &      *Zuu12-Lambda/MW*Mq*Zdd11*Zuu12/cb/sb*xvev+Lambda
     &      /MW*Mq*Zdd11*Zuu12/cb*sb*xvev+Lambda/MW*Mq*Zdd12
     &      *Zuu11/cb*sb*xvev)
      AAABR(65) = EE/SW*Sqrt2*(MW*Zb11*Zt11-MbMM**2/MW*Zb11
     &      *Zt11+MtMM**2/MW*Zb11*Zt11-2D0*MW*Zb11*Zt11*sb
     &      **2-Ab/MW*MbMM*Zb12*Zt11+At/MW*MtMM*Zb11*Zt12
     &      +Lambda/MW*MbMM*Zb12*Zt11/cb*sb*xvev-Lambda/MW
     &      *MtMM*Zb11*Zt12/cb/sb*xvev+Lambda/MW*MtMM*Zb11
     &      *Zt12/cb*sb*xvev)
      AAABR(66) = EE/SW*Sqrt2*(MW*Zb11*Zt21-MbMM**2/MW*Zb11
     &      *Zt21+MtMM**2/MW*Zb11*Zt21-2D0*MW*Zb11*Zt21*sb
     &      **2-Ab/MW*MbMM*Zb12*Zt21+At/MW*MtMM*Zb11*Zt22
     &      +Lambda/MW*MbMM*Zb12*Zt21/cb*sb*xvev-Lambda/MW
     &      *MtMM*Zb11*Zt22/cb/sb*xvev+Lambda/MW*MtMM*Zb11
     &      *Zt22/cb*sb*xvev)
      AAABR(67) = EE/SW*Sqrt2*(MW*Zcc21*Zss11+Mc**2/MW*Zcc21
     &      *Zss11-Mq**2/MW*Zcc21*Zss11-2D0*MW*Zcc21*Zss11
     &      *sb**2-Ad/MW*Mq*Zcc21*Zss12+Au/MW*Mc*Zcc22*Zss11
     &      -Lambda/MW*Mc*Zcc22*Zss11/cb/sb*xvev+Lambda/MW
     &      *Mc*Zcc22*Zss11/cb*sb*xvev+Lambda/MW*Mq*Zcc21
     &      *Zss12/cb*sb*xvev)
      AAABR(68) = EE/SW*Sqrt2*(MW*Zdd11*Zuu21-2D0*MW*Zdd11
     &      *Zuu21*sb**2-Ad/MW*Mq*Zdd12*Zuu21+Au/MW*Mq*Zdd11
     &      *Zuu22-Lambda/MW*Mq*Zdd11*Zuu22/cb/sb*xvev+Lambda
     &      /MW*Mq*Zdd11*Zuu22/cb*sb*xvev+Lambda/MW*Mq*Zdd12
     &      *Zuu21/cb*sb*xvev)
      AAABR(69) = EE/MW*Mq/SW*(Ad*Zss11*Zss22-Ad*Zss12*Zss21
     &      -Lambda*Zss11*Zss22/cb*sb*xvev+Lambda*Zss12*Zss21
     &      /cb*sb*xvev)
      AAABR(70) = EE/MW*MbMM/SW*(Ab*Zb11*Zb22-Ab*Zb12*Zb21
     &      -Lambda*Zb11*Zb22/cb*sb*xvev+Lambda*Zb12*Zb21
     &      /cb*sb*xvev)
      AAABR(71) = EE/MW*Mq/SW*(Ad*Zdd11*Zdd22-Ad*Zdd12*Zdd21
     &      -Lambda*Zdd11*Zdd22/cb*sb*xvev+Lambda*Zdd12*Zdd21
     &      /cb*sb*xvev)
      AAABR(72) = Mq*(Lambda*Pa22*Sqrt2*Zss11*Zss22/cb*sb
     &      -Lambda*Pa22*Sqrt2*Zss12*Zss21/cb*sb+EE*Lambda
     &      /MW*Pa21/SW*Zss11*Zss22*xvev-EE*Lambda/MW*Pa21
     &      /SW*Zss12*Zss21*xvev+Ad*EE/MW*Pa21/SW*Zss11*Zss22
     &      /cb*sb-Ad*EE/MW*Pa21/SW*Zss12*Zss21/cb*sb)
      AAABR(73) = Mq*(Lambda*Pa12*Sqrt2*Zdd11*Zdd22/cb*sb
     &      -Lambda*Pa12*Sqrt2*Zdd12*Zdd21/cb*sb+EE*Lambda
     &      /MW*Pa11/SW*Zdd11*Zdd22*xvev-EE*Lambda/MW*Pa11
     &      /SW*Zdd12*Zdd21*xvev+Ad*EE/MW*Pa11/SW*Zdd11*Zdd22
     &      /cb*sb-Ad*EE/MW*Pa11/SW*Zdd12*Zdd21/cb*sb)
      AAABR(74) = MbMM*(Lambda*Pa22*Sqrt2*Zb11*Zb22/cb*sb
     &      -Lambda*Pa22*Sqrt2*Zb12*Zb21/cb*sb+EE*Lambda
     &      /MW*Pa21/SW*Zb11*Zb22*xvev-EE*Lambda/MW*Pa21
     &      /SW*Zb12*Zb21*xvev+Ab*EE/MW*Pa21/SW*Zb11*Zb22
     &      /cb*sb-Ab*EE/MW*Pa21/SW*Zb12*Zb21/cb*sb)
      AAABR(75) = Mq*(Lambda*Pa12*Sqrt2*Zss11*Zss22/cb*sb
     &      -Lambda*Pa12*Sqrt2*Zss12*Zss21/cb*sb+EE*Lambda
     &      /MW*Pa11/SW*Zss11*Zss22*xvev-EE*Lambda/MW*Pa11
     &      /SW*Zss12*Zss21*xvev+Ad*EE/MW*Pa11/SW*Zss11*Zss22
     &      /cb*sb-Ad*EE/MW*Pa11/SW*Zss12*Zss21/cb*sb)
      AAABR(76) = MbMM*(Lambda*Pa12*Sqrt2*Zb11*Zb22/cb*sb
     &      -Lambda*Pa12*Sqrt2*Zb12*Zb21/cb*sb+EE*Lambda
     &      /MW*Pa11/SW*Zb11*Zb22*xvev-EE*Lambda/MW*Pa11
     &      /SW*Zb12*Zb21*xvev+Ab*EE/MW*Pa11/SW*Zb11*Zb22
     &      /cb*sb-Ab*EE/MW*Pa11/SW*Zb12*Zb21/cb*sb)
      AAABR(77) = Mq*(Lambda*Pa22*Sqrt2*Zdd11*Zdd22/cb*sb
     &      -Lambda*Pa22*Sqrt2*Zdd12*Zdd21/cb*sb+EE*Lambda
     &      /MW*Pa21/SW*Zdd11*Zdd22*xvev-EE*Lambda/MW*Pa21
     &      /SW*Zdd12*Zdd21*xvev+Ad*EE/MW*Pa21/SW*Zdd11*Zdd22
     &      /cb*sb-Ad*EE/MW*Pa21/SW*Zdd12*Zdd21/cb*sb)
      AAABR(78) = 6D0*EE/MW*Mq**2/SW*Zh31/cb-3D0*EE*MW/SW
     &      *Zh31*Zss11**2/cb+3D0*EE*MW/SW*Zh32*Zss11**2
     &      *sb-EE/CW**2*MW*SW*Zh31*Zss11**2/cb-2D0*EE/CW
     &      **2*MW*SW*Zh31*Zss12**2/cb+EE/CW**2*MW*SW*Zh32
     &      *Zss11**2*sb+2D0*EE/CW**2*MW*SW*Zh32*Zss12**2
     &      *sb+3D0*EE*MW/SW*Zh31*Zss11**2/cb*sb**2+EE/CW
     &      **2*MW*SW*Zh31*Zss11**2/cb*sb**2+2D0*EE/CW**2
     &      *MW*SW*Zh31*Zss12**2/cb*sb**2-6D0*Lambda*Mq*Sqrt2
     &      *Zh33*Zss11*Zss12/cb*sb+6D0*Ad*EE/MW*Mq/SW*Zh31
     &      *Zss11*Zss12/cb-6D0*EE*Lambda/MW*Mq/SW*Zh32*Zss11
     &      *Zss12/cb*xvev
      AAABR(79) = 6D0*EE/MW*Mq**2/SW*Zh11/cb-3D0*EE*MW/SW
     &      *Zh11*Zss11**2/cb+3D0*EE*MW/SW*Zh12*Zss11**2
     &      *sb-EE/CW**2*MW*SW*Zh11*Zss11**2/cb-2D0*EE/CW
     &      **2*MW*SW*Zh11*Zss12**2/cb+EE/CW**2*MW*SW*Zh12
     &      *Zss11**2*sb+2D0*EE/CW**2*MW*SW*Zh12*Zss12**2
     &      *sb+3D0*EE*MW/SW*Zh11*Zss11**2/cb*sb**2+EE/CW
     &      **2*MW*SW*Zh11*Zss11**2/cb*sb**2+2D0*EE/CW**2
     &      *MW*SW*Zh11*Zss12**2/cb*sb**2-6D0*Lambda*Mq*Sqrt2
     &      *Zh13*Zss11*Zss12/cb*sb+6D0*Ad*EE/MW*Mq/SW*Zh11
     &      *Zss11*Zss12/cb-6D0*EE*Lambda/MW*Mq/SW*Zh12*Zss11
     &      *Zss12/cb*xvev
      AAABR(80) = 6D0*EE/MW*Mq**2/SW*Zh11/cb-3D0*EE*MW/SW
     &      *Zdd11**2*Zh11/cb+3D0*EE*MW/SW*Zdd11**2*Zh12
     &      *sb-EE/CW**2*MW*SW*Zdd11**2*Zh11/cb+EE/CW**2
     &      *MW*SW*Zdd11**2*Zh12*sb-2D0*EE/CW**2*MW*SW*Zdd12
     &      **2*Zh11/cb+2D0*EE/CW**2*MW*SW*Zdd12**2*Zh12
     &      *sb+3D0*EE*MW/SW*Zdd11**2*Zh11/cb*sb**2+EE/CW
     &      **2*MW*SW*Zdd11**2*Zh11/cb*sb**2+2D0*EE/CW**2
     &      *MW*SW*Zdd12**2*Zh11/cb*sb**2-6D0*Lambda*Mq*Sqrt2
     &      *Zdd11*Zdd12*Zh13/cb*sb+6D0*Ad*EE/MW*Mq/SW*Zdd11
     &      *Zdd12*Zh11/cb-6D0*EE*Lambda/MW*Mq/SW*Zdd11*Zdd12
     &      *Zh12/cb*xvev
      AAABR(81) = 6D0*EE/MW*Mq**2/SW*Zh21/cb-3D0*EE*MW/SW
     &      *Zdd11**2*Zh21/cb+3D0*EE*MW/SW*Zdd11**2*Zh22
     &      *sb-EE/CW**2*MW*SW*Zdd11**2*Zh21/cb+EE/CW**2
     &      *MW*SW*Zdd11**2*Zh22*sb-2D0*EE/CW**2*MW*SW*Zdd12
     &      **2*Zh21/cb+2D0*EE/CW**2*MW*SW*Zdd12**2*Zh22
     &      *sb+3D0*EE*MW/SW*Zdd11**2*Zh21/cb*sb**2+EE/CW
     &      **2*MW*SW*Zdd11**2*Zh21/cb*sb**2+2D0*EE/CW**2
     &      *MW*SW*Zdd12**2*Zh21/cb*sb**2-6D0*Lambda*Mq*Sqrt2
     &      *Zdd11*Zdd12*Zh23/cb*sb+6D0*Ad*EE/MW*Mq/SW*Zdd11
     &      *Zdd12*Zh21/cb-6D0*EE*Lambda/MW*Mq/SW*Zdd11*Zdd12
     &      *Zh22/cb*xvev
      AAABR(82) = 6D0*EE/MW*Mq**2/SW*Zh31/cb-3D0*EE*MW/SW
     &      *Zdd11**2*Zh31/cb+3D0*EE*MW/SW*Zdd11**2*Zh32
     &      *sb-EE/CW**2*MW*SW*Zdd11**2*Zh31/cb+EE/CW**2
     &      *MW*SW*Zdd11**2*Zh32*sb-2D0*EE/CW**2*MW*SW*Zdd12
     &      **2*Zh31/cb+2D0*EE/CW**2*MW*SW*Zdd12**2*Zh32
     &      *sb+3D0*EE*MW/SW*Zdd11**2*Zh31/cb*sb**2+EE/CW
     &      **2*MW*SW*Zdd11**2*Zh31/cb*sb**2+2D0*EE/CW**2
     &      *MW*SW*Zdd12**2*Zh31/cb*sb**2-6D0*Lambda*Mq*Sqrt2
     &      *Zdd11*Zdd12*Zh33/cb*sb+6D0*Ad*EE/MW*Mq/SW*Zdd11
     &      *Zdd12*Zh31/cb-6D0*EE*Lambda/MW*Mq/SW*Zdd11*Zdd12
     &      *Zh32/cb*xvev
      AAABR(83) = 6D0*EE/MW*MbMM**2/SW*Zh11/cb-3D0*EE*MW
     &      /SW*Zb11**2*Zh11/cb+3D0*EE*MW/SW*Zb11**2*Zh12
     &      *sb-EE/CW**2*MW*SW*Zb11**2*Zh11/cb+EE/CW**2*MW
     &      *SW*Zb11**2*Zh12*sb-2D0*EE/CW**2*MW*SW*Zb12**2
     &      *Zh11/cb+2D0*EE/CW**2*MW*SW*Zb12**2*Zh12*sb+3D0
     &      *EE*MW/SW*Zb11**2*Zh11/cb*sb**2+EE/CW**2*MW*SW
     &      *Zb11**2*Zh11/cb*sb**2+2D0*EE/CW**2*MW*SW*Zb12
     &      **2*Zh11/cb*sb**2-6D0*Lambda*MbMM*Sqrt2*Zb11
     &      *Zb12*Zh13/cb*sb+6D0*Ab*EE/MW*MbMM/SW*Zb11*Zb12
     &      *Zh11/cb-6D0*EE*Lambda/MW*MbMM/SW*Zb11*Zb12*Zh12
     &      /cb*xvev
      AAABR(84) = 6D0*EE/MW*Mq**2/SW*Zh21/cb-3D0*EE*MW/SW
     &      *Zh21*Zss11**2/cb+3D0*EE*MW/SW*Zh22*Zss11**2
     &      *sb-EE/CW**2*MW*SW*Zh21*Zss11**2/cb-2D0*EE/CW
     &      **2*MW*SW*Zh21*Zss12**2/cb+EE/CW**2*MW*SW*Zh22
     &      *Zss11**2*sb+2D0*EE/CW**2*MW*SW*Zh22*Zss12**2
     &      *sb+3D0*EE*MW/SW*Zh21*Zss11**2/cb*sb**2+EE/CW
     &      **2*MW*SW*Zh21*Zss11**2/cb*sb**2+2D0*EE/CW**2
     &      *MW*SW*Zh21*Zss12**2/cb*sb**2-6D0*Lambda*Mq*Sqrt2
     &      *Zh23*Zss11*Zss12/cb*sb+6D0*Ad*EE/MW*Mq/SW*Zh21
     &      *Zss11*Zss12/cb-6D0*EE*Lambda/MW*Mq/SW*Zh22*Zss11
     &      *Zss12/cb*xvev
      AAABR(85) = 6D0*EE/MW*MbMM**2/SW*Zh31/cb-3D0*EE*MW
     &      /SW*Zb11**2*Zh31/cb+3D0*EE*MW/SW*Zb11**2*Zh32
     &      *sb-EE/CW**2*MW*SW*Zb11**2*Zh31/cb+EE/CW**2*MW
     &      *SW*Zb11**2*Zh32*sb-2D0*EE/CW**2*MW*SW*Zb12**2
     &      *Zh31/cb+2D0*EE/CW**2*MW*SW*Zb12**2*Zh32*sb+3D0
     &      *EE*MW/SW*Zb11**2*Zh31/cb*sb**2+EE/CW**2*MW*SW
     &      *Zb11**2*Zh31/cb*sb**2+2D0*EE/CW**2*MW*SW*Zb12
     &      **2*Zh31/cb*sb**2-6D0*Lambda*MbMM*Sqrt2*Zb11
     &      *Zb12*Zh33/cb*sb+6D0*Ab*EE/MW*MbMM/SW*Zb11*Zb12
     &      *Zh31/cb-6D0*EE*Lambda/MW*MbMM/SW*Zb11*Zb12*Zh32
     &      /cb*xvev
      AAABR(86) = 6D0*EE/MW*MbMM**2/SW*Zh21/cb-3D0*EE*MW
     &      /SW*Zb11**2*Zh21/cb+3D0*EE*MW/SW*Zb11**2*Zh22
     &      *sb-EE/CW**2*MW*SW*Zb11**2*Zh21/cb+EE/CW**2*MW
     &      *SW*Zb11**2*Zh22*sb-2D0*EE/CW**2*MW*SW*Zb12**2
     &      *Zh21/cb+2D0*EE/CW**2*MW*SW*Zb12**2*Zh22*sb+3D0
     &      *EE*MW/SW*Zb11**2*Zh21/cb*sb**2+EE/CW**2*MW*SW
     &      *Zb11**2*Zh21/cb*sb**2+2D0*EE/CW**2*MW*SW*Zb12
     &      **2*Zh21/cb*sb**2-6D0*Lambda*MbMM*Sqrt2*Zb11
     &      *Zb12*Zh23/cb*sb+6D0*Ab*EE/MW*MbMM/SW*Zb11*Zb12
     &      *Zh21/cb-6D0*EE*Lambda/MW*MbMM/SW*Zb11*Zb12*Zh22
     &      /cb*xvev
      AAABR(87) = 3D0*EE*MW/SW*Zdd11*Zdd21*Zh21/cb-3D0*EE
     &      *MW/SW*Zdd11*Zdd21*Zh22*sb+EE/CW**2*MW*SW*Zdd11
     &      *Zdd21*Zh21/cb-EE/CW**2*MW*SW*Zdd11*Zdd21*Zh22
     &      *sb+2D0*EE/CW**2*MW*SW*Zdd12*Zdd22*Zh21/cb-2D0
     &      *EE/CW**2*MW*SW*Zdd12*Zdd22*Zh22*sb-3D0*EE*MW
     &      /SW*Zdd11*Zdd21*Zh21/cb*sb**2+3D0*Lambda*Mq*Sqrt2
     &      *Zdd11*Zdd22*Zh23/cb*sb+3D0*Lambda*Mq*Sqrt2*Zdd12
     &      *Zdd21*Zh23/cb*sb-3D0*Ad*EE/MW*Mq/SW*Zdd11*Zdd22
     &      *Zh21/cb-3D0*Ad*EE/MW*Mq/SW*Zdd12*Zdd21*Zh21
     &      /cb-EE/CW**2*MW*SW*Zdd11*Zdd21*Zh21/cb*sb**2
     &      -2D0*EE/CW**2*MW*SW*Zdd12*Zdd22*Zh21/cb*sb**2
     &      +3D0*EE*Lambda/MW*Mq/SW*Zdd11*Zdd22*Zh22/cb*xvev
     &      +3D0*EE*Lambda/MW*Mq/SW*Zdd12*Zdd21*Zh22/cb*xvev
      AAABR(88) = 3D0*EE*MW/SW*Zb11*Zb21*Zh11/cb-3D0*EE*MW
     &      /SW*Zb11*Zb21*Zh12*sb+EE/CW**2*MW*SW*Zb11*Zb21
     &      *Zh11/cb-EE/CW**2*MW*SW*Zb11*Zb21*Zh12*sb+2D0
     &      *EE/CW**2*MW*SW*Zb12*Zb22*Zh11/cb-2D0*EE/CW**2
     &      *MW*SW*Zb12*Zb22*Zh12*sb-3D0*EE*MW/SW*Zb11*Zb21
     &      *Zh11/cb*sb**2+3D0*Lambda*MbMM*Sqrt2*Zb11*Zb22
     &      *Zh13/cb*sb+3D0*Lambda*MbMM*Sqrt2*Zb12*Zb21*Zh13
     &      /cb*sb-3D0*Ab*EE/MW*MbMM/SW*Zb11*Zb22*Zh11/cb
     &      -3D0*Ab*EE/MW*MbMM/SW*Zb12*Zb21*Zh11/cb-EE/CW
     &      **2*MW*SW*Zb11*Zb21*Zh11/cb*sb**2-2D0*EE/CW**2
     &      *MW*SW*Zb12*Zb22*Zh11/cb*sb**2+3D0*EE*Lambda
     &      /MW*MbMM/SW*Zb11*Zb22*Zh12/cb*xvev+3D0*EE*Lambda
     &      /MW*MbMM/SW*Zb12*Zb21*Zh12/cb*xvev
      AAABR(89) = 3D0*EE*MW/SW*Zb11*Zb21*Zh31/cb-3D0*EE*MW
     &      /SW*Zb11*Zb21*Zh32*sb+EE/CW**2*MW*SW*Zb11*Zb21
     &      *Zh31/cb-EE/CW**2*MW*SW*Zb11*Zb21*Zh32*sb+2D0
     &      *EE/CW**2*MW*SW*Zb12*Zb22*Zh31/cb-2D0*EE/CW**2
     &      *MW*SW*Zb12*Zb22*Zh32*sb-3D0*EE*MW/SW*Zb11*Zb21
     &      *Zh31/cb*sb**2+3D0*Lambda*MbMM*Sqrt2*Zb11*Zb22
     &      *Zh33/cb*sb+3D0*Lambda*MbMM*Sqrt2*Zb12*Zb21*Zh33
     &      /cb*sb-3D0*Ab*EE/MW*MbMM/SW*Zb11*Zb22*Zh31/cb
     &      -3D0*Ab*EE/MW*MbMM/SW*Zb12*Zb21*Zh31/cb-EE/CW
     &      **2*MW*SW*Zb11*Zb21*Zh31/cb*sb**2-2D0*EE/CW**2
     &      *MW*SW*Zb12*Zb22*Zh31/cb*sb**2+3D0*EE*Lambda
     &      /MW*MbMM/SW*Zb11*Zb22*Zh32/cb*xvev+3D0*EE*Lambda
     &      /MW*MbMM/SW*Zb12*Zb21*Zh32/cb*xvev
      AAABR(90) = 3D0*EE*MW/SW*Zdd11*Zdd21*Zh11/cb-3D0*EE
     &      *MW/SW*Zdd11*Zdd21*Zh12*sb+EE/CW**2*MW*SW*Zdd11
     &      *Zdd21*Zh11/cb-EE/CW**2*MW*SW*Zdd11*Zdd21*Zh12
     &      *sb+2D0*EE/CW**2*MW*SW*Zdd12*Zdd22*Zh11/cb-2D0
     &      *EE/CW**2*MW*SW*Zdd12*Zdd22*Zh12*sb-3D0*EE*MW
     &      /SW*Zdd11*Zdd21*Zh11/cb*sb**2+3D0*Lambda*Mq*Sqrt2
     &      *Zdd11*Zdd22*Zh13/cb*sb+3D0*Lambda*Mq*Sqrt2*Zdd12
     &      *Zdd21*Zh13/cb*sb-3D0*Ad*EE/MW*Mq/SW*Zdd11*Zdd22
     &      *Zh11/cb-3D0*Ad*EE/MW*Mq/SW*Zdd12*Zdd21*Zh11
     &      /cb-EE/CW**2*MW*SW*Zdd11*Zdd21*Zh11/cb*sb**2
     &      -2D0*EE/CW**2*MW*SW*Zdd12*Zdd22*Zh11/cb*sb**2
     &      +3D0*EE*Lambda/MW*Mq/SW*Zdd11*Zdd22*Zh12/cb*xvev
     &      +3D0*EE*Lambda/MW*Mq/SW*Zdd12*Zdd21*Zh12/cb*xvev
      AAABR(91) = 3D0*EE*MW/SW*Zdd11*Zdd21*Zh31/cb-3D0*EE
     &      *MW/SW*Zdd11*Zdd21*Zh32*sb+EE/CW**2*MW*SW*Zdd11
     &      *Zdd21*Zh31/cb-EE/CW**2*MW*SW*Zdd11*Zdd21*Zh32
     &      *sb+2D0*EE/CW**2*MW*SW*Zdd12*Zdd22*Zh31/cb-2D0
     &      *EE/CW**2*MW*SW*Zdd12*Zdd22*Zh32*sb-3D0*EE*MW
     &      /SW*Zdd11*Zdd21*Zh31/cb*sb**2+3D0*Lambda*Mq*Sqrt2
     &      *Zdd11*Zdd22*Zh33/cb*sb+3D0*Lambda*Mq*Sqrt2*Zdd12
     &      *Zdd21*Zh33/cb*sb-3D0*Ad*EE/MW*Mq/SW*Zdd11*Zdd22
     &      *Zh31/cb-3D0*Ad*EE/MW*Mq/SW*Zdd12*Zdd21*Zh31
     &      /cb-EE/CW**2*MW*SW*Zdd11*Zdd21*Zh31/cb*sb**2
     &      -2D0*EE/CW**2*MW*SW*Zdd12*Zdd22*Zh31/cb*sb**2
     &      +3D0*EE*Lambda/MW*Mq/SW*Zdd11*Zdd22*Zh32/cb*xvev
     &      +3D0*EE*Lambda/MW*Mq/SW*Zdd12*Zdd21*Zh32/cb*xvev
      AAABR(92) = 3D0*EE*MW/SW*Zb11*Zb21*Zh21/cb-3D0*EE*MW
     &      /SW*Zb11*Zb21*Zh22*sb+EE/CW**2*MW*SW*Zb11*Zb21
     &      *Zh21/cb-EE/CW**2*MW*SW*Zb11*Zb21*Zh22*sb+2D0
     &      *EE/CW**2*MW*SW*Zb12*Zb22*Zh21/cb-2D0*EE/CW**2
     &      *MW*SW*Zb12*Zb22*Zh22*sb-3D0*EE*MW/SW*Zb11*Zb21
     &      *Zh21/cb*sb**2+3D0*Lambda*MbMM*Sqrt2*Zb11*Zb22
     &      *Zh23/cb*sb+3D0*Lambda*MbMM*Sqrt2*Zb12*Zb21*Zh23
     &      /cb*sb-3D0*Ab*EE/MW*MbMM/SW*Zb11*Zb22*Zh21/cb
     &      -3D0*Ab*EE/MW*MbMM/SW*Zb12*Zb21*Zh21/cb-EE/CW
     &      **2*MW*SW*Zb11*Zb21*Zh21/cb*sb**2-2D0*EE/CW**2
     &      *MW*SW*Zb12*Zb22*Zh21/cb*sb**2+3D0*EE*Lambda
     &      /MW*MbMM/SW*Zb11*Zb22*Zh22/cb*xvev+3D0*EE*Lambda
     &      /MW*MbMM/SW*Zb12*Zb21*Zh22/cb*xvev
      AAABR(93) = 3D0*EE*MW/SW*Zh31*Zss11*Zss21/cb-3D0*EE
     &      *MW/SW*Zh32*Zss11*Zss21*sb+EE/CW**2*MW*SW*Zh31
     &      *Zss11*Zss21/cb+2D0*EE/CW**2*MW*SW*Zh31*Zss12
     &      *Zss22/cb-EE/CW**2*MW*SW*Zh32*Zss11*Zss21*sb
     &      -2D0*EE/CW**2*MW*SW*Zh32*Zss12*Zss22*sb-3D0*EE
     &      *MW/SW*Zh31*Zss11*Zss21/cb*sb**2+3D0*Lambda*Mq
     &      *Sqrt2*Zh33*Zss11*Zss22/cb*sb+3D0*Lambda*Mq*Sqrt2
     &      *Zh33*Zss12*Zss21/cb*sb-3D0*Ad*EE/MW*Mq/SW*Zh31
     &      *Zss11*Zss22/cb-3D0*Ad*EE/MW*Mq/SW*Zh31*Zss12
     &      *Zss21/cb-EE/CW**2*MW*SW*Zh31*Zss11*Zss21/cb
     &      *sb**2-2D0*EE/CW**2*MW*SW*Zh31*Zss12*Zss22/cb
     &      *sb**2+3D0*EE*Lambda/MW*Mq/SW*Zh32*Zss11*Zss22
     &      /cb*xvev+3D0*EE*Lambda/MW*Mq/SW*Zh32*Zss12*Zss21
     &      /cb*xvev
      AAABR(94) = 3D0*EE*MW/SW*Zh11*Zss11*Zss21/cb-3D0*EE
     &      *MW/SW*Zh12*Zss11*Zss21*sb+EE/CW**2*MW*SW*Zh11
     &      *Zss11*Zss21/cb+2D0*EE/CW**2*MW*SW*Zh11*Zss12
     &      *Zss22/cb-EE/CW**2*MW*SW*Zh12*Zss11*Zss21*sb
     &      -2D0*EE/CW**2*MW*SW*Zh12*Zss12*Zss22*sb-3D0*EE
     &      *MW/SW*Zh11*Zss11*Zss21/cb*sb**2+3D0*Lambda*Mq
     &      *Sqrt2*Zh13*Zss11*Zss22/cb*sb+3D0*Lambda*Mq*Sqrt2
     &      *Zh13*Zss12*Zss21/cb*sb-3D0*Ad*EE/MW*Mq/SW*Zh11
     &      *Zss11*Zss22/cb-3D0*Ad*EE/MW*Mq/SW*Zh11*Zss12
     &      *Zss21/cb-EE/CW**2*MW*SW*Zh11*Zss11*Zss21/cb
     &      *sb**2-2D0*EE/CW**2*MW*SW*Zh11*Zss12*Zss22/cb
     &      *sb**2+3D0*EE*Lambda/MW*Mq/SW*Zh12*Zss11*Zss22
     &      /cb*xvev+3D0*EE*Lambda/MW*Mq/SW*Zh12*Zss12*Zss21
     &      /cb*xvev
      AAABR(95) = 3D0*EE*MW/SW*Zh21*Zss11*Zss21/cb-3D0*EE
     &      *MW/SW*Zh22*Zss11*Zss21*sb+EE/CW**2*MW*SW*Zh21
     &      *Zss11*Zss21/cb+2D0*EE/CW**2*MW*SW*Zh21*Zss12
     &      *Zss22/cb-EE/CW**2*MW*SW*Zh22*Zss11*Zss21*sb
     &      -2D0*EE/CW**2*MW*SW*Zh22*Zss12*Zss22*sb-3D0*EE
     &      *MW/SW*Zh21*Zss11*Zss21/cb*sb**2+3D0*Lambda*Mq
     &      *Sqrt2*Zh23*Zss11*Zss22/cb*sb+3D0*Lambda*Mq*Sqrt2
     &      *Zh23*Zss12*Zss21/cb*sb-3D0*Ad*EE/MW*Mq/SW*Zh21
     &      *Zss11*Zss22/cb-3D0*Ad*EE/MW*Mq/SW*Zh21*Zss12
     &      *Zss21/cb-EE/CW**2*MW*SW*Zh21*Zss11*Zss21/cb
     &      *sb**2-2D0*EE/CW**2*MW*SW*Zh21*Zss12*Zss22/cb
     &      *sb**2+3D0*EE*Lambda/MW*Mq/SW*Zh22*Zss11*Zss22
     &      /cb*xvev+3D0*EE*Lambda/MW*Mq/SW*Zh22*Zss12*Zss21
     &      /cb*xvev
      AAABR(96) = EE/SW*Sqrt2*(MW*Zcc11*Zss21+Mc**2/MW*Zcc11
     &      *Zss21-Mq**2/MW*Zcc11*Zss21-2D0*MW*Zcc11*Zss21
     &      *sb**2-Ad/MW*Mq*Zcc11*Zss22+Au/MW*Mc*Zcc12*Zss21
     &      -Lambda/MW*Mc*Zcc12*Zss21/cb/sb*xvev+Lambda/MW
     &      *Mc*Zcc12*Zss21/cb*sb*xvev+Lambda/MW*Mq*Zcc11
     &      *Zss22/cb*sb*xvev)
      AAABR(97) = EE/SW*Sqrt2*(MW*Zdd21*Zuu11-2D0*MW*Zdd21
     &      *Zuu11*sb**2-Ad/MW*Mq*Zdd22*Zuu11+Au/MW*Mq*Zdd21
     &      *Zuu12-Lambda/MW*Mq*Zdd21*Zuu12/cb/sb*xvev+Lambda
     &      /MW*Mq*Zdd21*Zuu12/cb*sb*xvev+Lambda/MW*Mq*Zdd22
     &      *Zuu11/cb*sb*xvev)
      AAABR(98) = EE/SW*Sqrt2*(MW*Zb21*Zt11-MbMM**2/MW*Zb21
     &      *Zt11+MtMM**2/MW*Zb21*Zt11-2D0*MW*Zb21*Zt11*sb
     &      **2-Ab/MW*MbMM*Zb22*Zt11+At/MW*MtMM*Zb21*Zt12
     &      +Lambda/MW*MbMM*Zb22*Zt11/cb*sb*xvev-Lambda/MW
     &      *MtMM*Zb21*Zt12/cb/sb*xvev+Lambda/MW*MtMM*Zb21
     &      *Zt12/cb*sb*xvev)
      AAABR(99) = EE/SW*Sqrt2*(MW*Zb21*Zt21-MbMM**2/MW*Zb21
     &      *Zt21+MtMM**2/MW*Zb21*Zt21-2D0*MW*Zb21*Zt21*sb
     &      **2-Ab/MW*MbMM*Zb22*Zt21+At/MW*MtMM*Zb21*Zt22
     &      +Lambda/MW*MbMM*Zb22*Zt21/cb*sb*xvev-Lambda/MW
     &      *MtMM*Zb21*Zt22/cb/sb*xvev+Lambda/MW*MtMM*Zb21
     &      *Zt22/cb*sb*xvev)
      AAABR(100) = EE/SW*Sqrt2*(MW*Zdd21*Zuu21-2D0*MW*Zdd21
     &      *Zuu21*sb**2-Ad/MW*Mq*Zdd22*Zuu21+Au/MW*Mq*Zdd21
     &      *Zuu22-Lambda/MW*Mq*Zdd21*Zuu22/cb/sb*xvev+Lambda
     &      /MW*Mq*Zdd21*Zuu22/cb*sb*xvev+Lambda/MW*Mq*Zdd22
     &      *Zuu21/cb*sb*xvev)
      AAABR(101) = EE/SW*Sqrt2*(MW*Zcc21*Zss21+Mc**2/MW*Zcc21
     &      *Zss21-Mq**2/MW*Zcc21*Zss21-2D0*MW*Zcc21*Zss21
     &      *sb**2-Ad/MW*Mq*Zcc21*Zss22+Au/MW*Mc*Zcc22*Zss21
     &      -Lambda/MW*Mc*Zcc22*Zss21/cb/sb*xvev+Lambda/MW
     &      *Mc*Zcc22*Zss21/cb*sb*xvev+Lambda/MW*Mq*Zcc21
     &      *Zss22/cb*sb*xvev)
      AAABR(102) = 6D0*EE/MW*Mq**2/SW*Zh11/cb-3D0*EE*MW/SW
     &      *Zh11*Zss21**2/cb+3D0*EE*MW/SW*Zh12*Zss21**2
     &      *sb-EE/CW**2*MW*SW*Zh11*Zss21**2/cb-2D0*EE/CW
     &      **2*MW*SW*Zh11*Zss22**2/cb+EE/CW**2*MW*SW*Zh12
     &      *Zss21**2*sb+2D0*EE/CW**2*MW*SW*Zh12*Zss22**2
     &      *sb+3D0*EE*MW/SW*Zh11*Zss21**2/cb*sb**2+EE/CW
     &      **2*MW*SW*Zh11*Zss21**2/cb*sb**2+2D0*EE/CW**2
     &      *MW*SW*Zh11*Zss22**2/cb*sb**2-6D0*Lambda*Mq*Sqrt2
     &      *Zh13*Zss21*Zss22/cb*sb+6D0*Ad*EE/MW*Mq/SW*Zh11
     &      *Zss21*Zss22/cb-6D0*EE*Lambda/MW*Mq/SW*Zh12*Zss21
     &      *Zss22/cb*xvev
      AAABR(103) = 6D0*EE/MW*MbMM**2/SW*Zh31/cb-3D0*EE*MW
     &      /SW*Zb21**2*Zh31/cb+3D0*EE*MW/SW*Zb21**2*Zh32
     &      *sb-EE/CW**2*MW*SW*Zb21**2*Zh31/cb+EE/CW**2*MW
     &      *SW*Zb21**2*Zh32*sb-2D0*EE/CW**2*MW*SW*Zb22**2
     &      *Zh31/cb+2D0*EE/CW**2*MW*SW*Zb22**2*Zh32*sb+3D0
     &      *EE*MW/SW*Zb21**2*Zh31/cb*sb**2+EE/CW**2*MW*SW
     &      *Zb21**2*Zh31/cb*sb**2+2D0*EE/CW**2*MW*SW*Zb22
     &      **2*Zh31/cb*sb**2-6D0*Lambda*MbMM*Sqrt2*Zb21
     &      *Zb22*Zh33/cb*sb+6D0*Ab*EE/MW*MbMM/SW*Zb21*Zb22
     &      *Zh31/cb-6D0*EE*Lambda/MW*MbMM/SW*Zb21*Zb22*Zh32
     &      /cb*xvev
      AAABR(104) = 6D0*EE/MW*Mq**2/SW*Zh31/cb-3D0*EE*MW/SW
     &      *Zdd21**2*Zh31/cb+3D0*EE*MW/SW*Zdd21**2*Zh32
     &      *sb-EE/CW**2*MW*SW*Zdd21**2*Zh31/cb+EE/CW**2
     &      *MW*SW*Zdd21**2*Zh32*sb-2D0*EE/CW**2*MW*SW*Zdd22
     &      **2*Zh31/cb+2D0*EE/CW**2*MW*SW*Zdd22**2*Zh32
     &      *sb+3D0*EE*MW/SW*Zdd21**2*Zh31/cb*sb**2+EE/CW
     &      **2*MW*SW*Zdd21**2*Zh31/cb*sb**2+2D0*EE/CW**2
     &      *MW*SW*Zdd22**2*Zh31/cb*sb**2-6D0*Lambda*Mq*Sqrt2
     &      *Zdd21*Zdd22*Zh33/cb*sb+6D0*Ad*EE/MW*Mq/SW*Zdd21
     &      *Zdd22*Zh31/cb-6D0*EE*Lambda/MW*Mq/SW*Zdd21*Zdd22
     &      *Zh32/cb*xvev
      AAABR(105) = 6D0*EE/MW*Mq**2/SW*Zh31/cb-3D0*EE*MW/SW
     &      *Zh31*Zss21**2/cb+3D0*EE*MW/SW*Zh32*Zss21**2
     &      *sb-EE/CW**2*MW*SW*Zh31*Zss21**2/cb-2D0*EE/CW
     &      **2*MW*SW*Zh31*Zss22**2/cb+EE/CW**2*MW*SW*Zh32
     &      *Zss21**2*sb+2D0*EE/CW**2*MW*SW*Zh32*Zss22**2
     &      *sb+3D0*EE*MW/SW*Zh31*Zss21**2/cb*sb**2+EE/CW
     &      **2*MW*SW*Zh31*Zss21**2/cb*sb**2+2D0*EE/CW**2
     &      *MW*SW*Zh31*Zss22**2/cb*sb**2-6D0*Lambda*Mq*Sqrt2
     &      *Zh33*Zss21*Zss22/cb*sb+6D0*Ad*EE/MW*Mq/SW*Zh31
     &      *Zss21*Zss22/cb-6D0*EE*Lambda/MW*Mq/SW*Zh32*Zss21
     &      *Zss22/cb*xvev
      AAABR(106) = 6D0*EE/MW*MbMM**2/SW*Zh21/cb-3D0*EE*MW
     &      /SW*Zb21**2*Zh21/cb+3D0*EE*MW/SW*Zb21**2*Zh22
     &      *sb-EE/CW**2*MW*SW*Zb21**2*Zh21/cb+EE/CW**2*MW
     &      *SW*Zb21**2*Zh22*sb-2D0*EE/CW**2*MW*SW*Zb22**2
     &      *Zh21/cb+2D0*EE/CW**2*MW*SW*Zb22**2*Zh22*sb+3D0
     &      *EE*MW/SW*Zb21**2*Zh21/cb*sb**2+EE/CW**2*MW*SW
     &      *Zb21**2*Zh21/cb*sb**2+2D0*EE/CW**2*MW*SW*Zb22
     &      **2*Zh21/cb*sb**2-6D0*Lambda*MbMM*Sqrt2*Zb21
     &      *Zb22*Zh23/cb*sb+6D0*Ab*EE/MW*MbMM/SW*Zb21*Zb22
     &      *Zh21/cb-6D0*EE*Lambda/MW*MbMM/SW*Zb21*Zb22*Zh22
     &      /cb*xvev
      AAABR(107) = 6D0*EE/MW*MbMM**2/SW*Zh11/cb-3D0*EE*MW
     &      /SW*Zb21**2*Zh11/cb+3D0*EE*MW/SW*Zb21**2*Zh12
     &      *sb-EE/CW**2*MW*SW*Zb21**2*Zh11/cb+EE/CW**2*MW
     &      *SW*Zb21**2*Zh12*sb-2D0*EE/CW**2*MW*SW*Zb22**2
     &      *Zh11/cb+2D0*EE/CW**2*MW*SW*Zb22**2*Zh12*sb+3D0
     &      *EE*MW/SW*Zb21**2*Zh11/cb*sb**2+EE/CW**2*MW*SW
     &      *Zb21**2*Zh11/cb*sb**2+2D0*EE/CW**2*MW*SW*Zb22
     &      **2*Zh11/cb*sb**2-6D0*Lambda*MbMM*Sqrt2*Zb21
     &      *Zb22*Zh13/cb*sb+6D0*Ab*EE/MW*MbMM/SW*Zb21*Zb22
     &      *Zh11/cb-6D0*EE*Lambda/MW*MbMM/SW*Zb21*Zb22*Zh12
     &      /cb*xvev
      AAABR(108) = 6D0*EE/MW*Mq**2/SW*Zh21/cb-3D0*EE*MW/SW
     &      *Zdd21**2*Zh21/cb+3D0*EE*MW/SW*Zdd21**2*Zh22
     &      *sb-EE/CW**2*MW*SW*Zdd21**2*Zh21/cb+EE/CW**2
     &      *MW*SW*Zdd21**2*Zh22*sb-2D0*EE/CW**2*MW*SW*Zdd22
     &      **2*Zh21/cb+2D0*EE/CW**2*MW*SW*Zdd22**2*Zh22
     &      *sb+3D0*EE*MW/SW*Zdd21**2*Zh21/cb*sb**2+EE/CW
     &      **2*MW*SW*Zdd21**2*Zh21/cb*sb**2+2D0*EE/CW**2
     &      *MW*SW*Zdd22**2*Zh21/cb*sb**2-6D0*Lambda*Mq*Sqrt2
     &      *Zdd21*Zdd22*Zh23/cb*sb+6D0*Ad*EE/MW*Mq/SW*Zdd21
     &      *Zdd22*Zh21/cb-6D0*EE*Lambda/MW*Mq/SW*Zdd21*Zdd22
     &      *Zh22/cb*xvev
      AAABR(109) = 6D0*EE/MW*Mq**2/SW*Zh21/cb-3D0*EE*MW/SW
     &      *Zh21*Zss21**2/cb+3D0*EE*MW/SW*Zh22*Zss21**2
     &      *sb-EE/CW**2*MW*SW*Zh21*Zss21**2/cb-2D0*EE/CW
     &      **2*MW*SW*Zh21*Zss22**2/cb+EE/CW**2*MW*SW*Zh22
     &      *Zss21**2*sb+2D0*EE/CW**2*MW*SW*Zh22*Zss22**2
     &      *sb+3D0*EE*MW/SW*Zh21*Zss21**2/cb*sb**2+EE/CW
     &      **2*MW*SW*Zh21*Zss21**2/cb*sb**2+2D0*EE/CW**2
     &      *MW*SW*Zh21*Zss22**2/cb*sb**2-6D0*Lambda*Mq*Sqrt2
     &      *Zh23*Zss21*Zss22/cb*sb+6D0*Ad*EE/MW*Mq/SW*Zh21
     &      *Zss21*Zss22/cb-6D0*EE*Lambda/MW*Mq/SW*Zh22*Zss21
     &      *Zss22/cb*xvev
      AAABR(110) = 6D0*EE/MW*Mq**2/SW*Zh11/cb-3D0*EE*MW/SW
     &      *Zdd21**2*Zh11/cb+3D0*EE*MW/SW*Zdd21**2*Zh12
     &      *sb-EE/CW**2*MW*SW*Zdd21**2*Zh11/cb+EE/CW**2
     &      *MW*SW*Zdd21**2*Zh12*sb-2D0*EE/CW**2*MW*SW*Zdd22
     &      **2*Zh11/cb+2D0*EE/CW**2*MW*SW*Zdd22**2*Zh12
     &      *sb+3D0*EE*MW/SW*Zdd21**2*Zh11/cb*sb**2+EE/CW
     &      **2*MW*SW*Zdd21**2*Zh11/cb*sb**2+2D0*EE/CW**2
     &      *MW*SW*Zdd22**2*Zh11/cb*sb**2-6D0*Lambda*Mq*Sqrt2
     &      *Zdd21*Zdd22*Zh13/cb*sb+6D0*Ad*EE/MW*Mq/SW*Zdd21
     &      *Zdd22*Zh11/cb-6D0*EE*Lambda/MW*Mq/SW*Zdd21*Zdd22
     &      *Zh12/cb*xvev
      AAABR(111) = EE/MW*Mc/SW*(Au*Zcc11*Zcc22-Au*Zcc12*Zcc21
     &      -Lambda*Zcc11*Zcc22*cb/sb*xvev+Lambda*Zcc12*Zcc21
     &      *cb/sb*xvev)
      AAABR(112) = EE/MW*MtMM/SW*(At*Zt11*Zt22-At*Zt12*Zt21
     &      -Lambda*Zt11*Zt22*cb/sb*xvev+Lambda*Zt12*Zt21
     &      *cb/sb*xvev)
      AAABR(113) = EE/MW*Mq/SW*(Au*Zuu11*Zuu22-Au*Zuu12*Zuu21
     &      -Lambda*Zuu11*Zuu22*cb/sb*xvev+Lambda*Zuu12*Zuu21
     &      *cb/sb*xvev)
      AAABR(114) = Mq*(Lambda*Pa22*Sqrt2*Zuu11*Zuu22*cb/sb
     &      -Lambda*Pa22*Sqrt2*Zuu12*Zuu21*cb/sb+EE*Lambda
     &      /MW*Pa21/SW*Zuu11*Zuu22*xvev-EE*Lambda/MW*Pa21
     &      /SW*Zuu12*Zuu21*xvev+Au*EE/MW*Pa21/SW*Zuu11*Zuu22
     &      *cb/sb-Au*EE/MW*Pa21/SW*Zuu12*Zuu21*cb/sb)
      AAABR(115) = Mc*(Lambda*Pa22*Sqrt2*Zcc11*Zcc22*cb/sb
     &      -Lambda*Pa22*Sqrt2*Zcc12*Zcc21*cb/sb+EE*Lambda
     &      /MW*Pa21/SW*Zcc11*Zcc22*xvev-EE*Lambda/MW*Pa21
     &      /SW*Zcc12*Zcc21*xvev+Au*EE/MW*Pa21/SW*Zcc11*Zcc22
     &      *cb/sb-Au*EE/MW*Pa21/SW*Zcc12*Zcc21*cb/sb)
      AAABR(116) = Mq*(Lambda*Pa12*Sqrt2*Zuu11*Zuu22*cb/sb
     &      -Lambda*Pa12*Sqrt2*Zuu12*Zuu21*cb/sb+EE*Lambda
     &      /MW*Pa11/SW*Zuu11*Zuu22*xvev-EE*Lambda/MW*Pa11
     &      /SW*Zuu12*Zuu21*xvev+Au*EE/MW*Pa11/SW*Zuu11*Zuu22
     &      *cb/sb-Au*EE/MW*Pa11/SW*Zuu12*Zuu21*cb/sb)
      AAABR(117) = MtMM*(Lambda*Pa22*Sqrt2*Zt11*Zt22*cb/sb
     &      -Lambda*Pa22*Sqrt2*Zt12*Zt21*cb/sb+EE*Lambda
     &      /MW*Pa21/SW*Zt11*Zt22*xvev-EE*Lambda/MW*Pa21
     &      /SW*Zt12*Zt21*xvev+At*EE/MW*Pa21/SW*Zt11*Zt22
     &      *cb/sb-At*EE/MW*Pa21/SW*Zt12*Zt21*cb/sb)
      AAABR(118) = Mc*(Lambda*Pa12*Sqrt2*Zcc11*Zcc22*cb/sb
     &      -Lambda*Pa12*Sqrt2*Zcc12*Zcc21*cb/sb+EE*Lambda
     &      /MW*Pa11/SW*Zcc11*Zcc22*xvev-EE*Lambda/MW*Pa11
     &      /SW*Zcc12*Zcc21*xvev+Au*EE/MW*Pa11/SW*Zcc11*Zcc22
     &      *cb/sb-Au*EE/MW*Pa11/SW*Zcc12*Zcc21*cb/sb)
      AAABR(119) = MtMM*(Lambda*Pa12*Sqrt2*Zt11*Zt22*cb/sb
     &      -Lambda*Pa12*Sqrt2*Zt12*Zt21*cb/sb+EE*Lambda
     &      /MW*Pa11/SW*Zt11*Zt22*xvev-EE*Lambda/MW*Pa11
     &      /SW*Zt12*Zt21*xvev+At*EE/MW*Pa11/SW*Zt11*Zt22
     &      *cb/sb-At*EE/MW*Pa11/SW*Zt12*Zt21*cb/sb)
      AAABR(120) = 6D0*EE/MW*Mc**2/SW*Zh12/sb+3D0*EE*MW/SW
     &      *Zcc11**2*Zh11*cb-3D0*EE*MW/SW*Zcc11**2*Zh12
     &      *sb-EE/CW**2*MW*SW*Zcc11**2*Zh11*cb+EE/CW**2
     &      *MW*SW*Zcc11**2*Zh12*sb+4D0*EE/CW**2*MW*SW*Zcc12
     &      **2*Zh11*cb-4D0*EE/CW**2*MW*SW*Zcc12**2*Zh12
     &      *sb-6D0*Lambda*Mc*Sqrt2*Zcc11*Zcc12*Zh13*cb/sb
     &      +6D0*Au*EE/MW*Mc/SW*Zcc11*Zcc12*Zh12/sb-6D0*EE
     &      *Lambda/MW*Mc/SW*Zcc11*Zcc12*Zh11/sb*xvev
      AAABR(121) = 6D0*EE/MW*Mq**2/SW*Zh12/sb+3D0*EE*MW/SW
     &      *Zh11*Zuu11**2*cb-3D0*EE*MW/SW*Zh12*Zuu11**2
     &      *sb-EE/CW**2*MW*SW*Zh11*Zuu11**2*cb+4D0*EE/CW
     &      **2*MW*SW*Zh11*Zuu12**2*cb+EE/CW**2*MW*SW*Zh12
     &      *Zuu11**2*sb-4D0*EE/CW**2*MW*SW*Zh12*Zuu12**2
     &      *sb-6D0*Lambda*Mq*Sqrt2*Zh13*Zuu11*Zuu12*cb/sb
     &      +6D0*Au*EE/MW*Mq/SW*Zh12*Zuu11*Zuu12/sb-6D0*EE
     &      *Lambda/MW*Mq/SW*Zh11*Zuu11*Zuu12/sb*xvev
      AAABR(122) = 6D0*EE/MW*MtMM**2/SW*Zh12/sb+3D0*EE*MW
     &      /SW*Zh11*Zt11**2*cb-3D0*EE*MW/SW*Zh12*Zt11**2
     &      *sb-EE/CW**2*MW*SW*Zh11*Zt11**2*cb+4D0*EE/CW
     &      **2*MW*SW*Zh11*Zt12**2*cb+EE/CW**2*MW*SW*Zh12
     &      *Zt11**2*sb-4D0*EE/CW**2*MW*SW*Zh12*Zt12**2*sb
     &      -6D0*Lambda*MtMM*Sqrt2*Zh13*Zt11*Zt12*cb/sb+6D0
     &      *At*EE/MW*MtMM/SW*Zh12*Zt11*Zt12/sb-6D0*EE*Lambda
     &      /MW*MtMM/SW*Zh11*Zt11*Zt12/sb*xvev
      AAABR(123) = 6D0*EE/MW*MtMM**2/SW*Zh32/sb+3D0*EE*MW
     &      /SW*Zh31*Zt11**2*cb-3D0*EE*MW/SW*Zh32*Zt11**2
     &      *sb-EE/CW**2*MW*SW*Zh31*Zt11**2*cb+4D0*EE/CW
     &      **2*MW*SW*Zh31*Zt12**2*cb+EE/CW**2*MW*SW*Zh32
     &      *Zt11**2*sb-4D0*EE/CW**2*MW*SW*Zh32*Zt12**2*sb
     &      -6D0*Lambda*MtMM*Sqrt2*Zh33*Zt11*Zt12*cb/sb+6D0
     &      *At*EE/MW*MtMM/SW*Zh32*Zt11*Zt12/sb-6D0*EE*Lambda
     &      /MW*MtMM/SW*Zh31*Zt11*Zt12/sb*xvev
      AAABR(124) = 6D0*EE/MW*Mc**2/SW*Zh22/sb+3D0*EE*MW/SW
     &      *Zcc11**2*Zh21*cb-3D0*EE*MW/SW*Zcc11**2*Zh22
     &      *sb-EE/CW**2*MW*SW*Zcc11**2*Zh21*cb+EE/CW**2
     &      *MW*SW*Zcc11**2*Zh22*sb+4D0*EE/CW**2*MW*SW*Zcc12
     &      **2*Zh21*cb-4D0*EE/CW**2*MW*SW*Zcc12**2*Zh22
     &      *sb-6D0*Lambda*Mc*Sqrt2*Zcc11*Zcc12*Zh23*cb/sb
     &      +6D0*Au*EE/MW*Mc/SW*Zcc11*Zcc12*Zh22/sb-6D0*EE
     &      *Lambda/MW*Mc/SW*Zcc11*Zcc12*Zh21/sb*xvev
      AAABR(125) = 6D0*EE/MW*MtMM**2/SW*Zh22/sb+3D0*EE*MW
     &      /SW*Zh21*Zt11**2*cb-3D0*EE*MW/SW*Zh22*Zt11**2
     &      *sb-EE/CW**2*MW*SW*Zh21*Zt11**2*cb+4D0*EE/CW
     &      **2*MW*SW*Zh21*Zt12**2*cb+EE/CW**2*MW*SW*Zh22
     &      *Zt11**2*sb-4D0*EE/CW**2*MW*SW*Zh22*Zt12**2*sb
     &      -6D0*Lambda*MtMM*Sqrt2*Zh23*Zt11*Zt12*cb/sb+6D0
     &      *At*EE/MW*MtMM/SW*Zh22*Zt11*Zt12/sb-6D0*EE*Lambda
     &      /MW*MtMM/SW*Zh21*Zt11*Zt12/sb*xvev
      AAABR(126) = 6D0*EE/MW*Mq**2/SW*Zh32/sb+3D0*EE*MW/SW
     &      *Zh31*Zuu11**2*cb-3D0*EE*MW/SW*Zh32*Zuu11**2
     &      *sb-EE/CW**2*MW*SW*Zh31*Zuu11**2*cb+4D0*EE/CW
     &      **2*MW*SW*Zh31*Zuu12**2*cb+EE/CW**2*MW*SW*Zh32
     &      *Zuu11**2*sb-4D0*EE/CW**2*MW*SW*Zh32*Zuu12**2
     &      *sb-6D0*Lambda*Mq*Sqrt2*Zh33*Zuu11*Zuu12*cb/sb
     &      +6D0*Au*EE/MW*Mq/SW*Zh32*Zuu11*Zuu12/sb-6D0*EE
     &      *Lambda/MW*Mq/SW*Zh31*Zuu11*Zuu12/sb*xvev
      AAABR(127) = 6D0*EE/MW*Mq**2/SW*Zh22/sb+3D0*EE*MW/SW
     &      *Zh21*Zuu11**2*cb-3D0*EE*MW/SW*Zh22*Zuu11**2
     &      *sb-EE/CW**2*MW*SW*Zh21*Zuu11**2*cb+4D0*EE/CW
     &      **2*MW*SW*Zh21*Zuu12**2*cb+EE/CW**2*MW*SW*Zh22
     &      *Zuu11**2*sb-4D0*EE/CW**2*MW*SW*Zh22*Zuu12**2
     &      *sb-6D0*Lambda*Mq*Sqrt2*Zh23*Zuu11*Zuu12*cb/sb
     &      +6D0*Au*EE/MW*Mq/SW*Zh22*Zuu11*Zuu12/sb-6D0*EE
     &      *Lambda/MW*Mq/SW*Zh21*Zuu11*Zuu12/sb*xvev
      AAABR(128) = 6D0*EE/MW*Mc**2/SW*Zh32/sb+3D0*EE*MW/SW
     &      *Zcc11**2*Zh31*cb-3D0*EE*MW/SW*Zcc11**2*Zh32
     &      *sb-EE/CW**2*MW*SW*Zcc11**2*Zh31*cb+EE/CW**2
     &      *MW*SW*Zcc11**2*Zh32*sb+4D0*EE/CW**2*MW*SW*Zcc12
     &      **2*Zh31*cb-4D0*EE/CW**2*MW*SW*Zcc12**2*Zh32
     &      *sb-6D0*Lambda*Mc*Sqrt2*Zcc11*Zcc12*Zh33*cb/sb
     &      +6D0*Au*EE/MW*Mc/SW*Zcc11*Zcc12*Zh32/sb-6D0*EE
     &      *Lambda/MW*Mc/SW*Zcc11*Zcc12*Zh31/sb*xvev
      AAABR(129) = 3D0*EE*MW/SW*Zh31*Zt11*Zt21*cb-3D0*EE
     &      *MW/SW*Zh32*Zt11*Zt21*sb-EE/CW**2*MW*SW*Zh31
     &      *Zt11*Zt21*cb+4D0*EE/CW**2*MW*SW*Zh31*Zt12*Zt22
     &      *cb+EE/CW**2*MW*SW*Zh32*Zt11*Zt21*sb-4D0*EE/CW
     &      **2*MW*SW*Zh32*Zt12*Zt22*sb-3D0*Lambda*MtMM*Sqrt2
     &      *Zh33*Zt11*Zt22*cb/sb-3D0*Lambda*MtMM*Sqrt2*Zh33
     &      *Zt12*Zt21*cb/sb+3D0*At*EE/MW*MtMM/SW*Zh32*Zt11
     &      *Zt22/sb+3D0*At*EE/MW*MtMM/SW*Zh32*Zt12*Zt21
     &      /sb-3D0*EE*Lambda/MW*MtMM/SW*Zh31*Zt11*Zt22/sb
     &      *xvev-3D0*EE*Lambda/MW*MtMM/SW*Zh31*Zt12*Zt21
     &      /sb*xvev
      AAABR(130) = 3D0*EE*MW/SW*Zcc11*Zcc21*Zh11*cb-3D0*EE
     &      *MW/SW*Zcc11*Zcc21*Zh12*sb-EE/CW**2*MW*SW*Zcc11
     &      *Zcc21*Zh11*cb+EE/CW**2*MW*SW*Zcc11*Zcc21*Zh12
     &      *sb+4D0*EE/CW**2*MW*SW*Zcc12*Zcc22*Zh11*cb-4D0
     &      *EE/CW**2*MW*SW*Zcc12*Zcc22*Zh12*sb-3D0*Lambda
     &      *Mc*Sqrt2*Zcc11*Zcc22*Zh13*cb/sb-3D0*Lambda*Mc
     &      *Sqrt2*Zcc12*Zcc21*Zh13*cb/sb+3D0*Au*EE/MW*Mc
     &      /SW*Zcc11*Zcc22*Zh12/sb+3D0*Au*EE/MW*Mc/SW*Zcc12
     &      *Zcc21*Zh12/sb-3D0*EE*Lambda/MW*Mc/SW*Zcc11*Zcc22
     &      *Zh11/sb*xvev-3D0*EE*Lambda/MW*Mc/SW*Zcc12*Zcc21
     &      *Zh11/sb*xvev
      AAABR(131) = 3D0*EE*MW/SW*Zh21*Zuu11*Zuu21*cb-3D0*EE
     &      *MW/SW*Zh22*Zuu11*Zuu21*sb-EE/CW**2*MW*SW*Zh21
     &      *Zuu11*Zuu21*cb+4D0*EE/CW**2*MW*SW*Zh21*Zuu12
     &      *Zuu22*cb+EE/CW**2*MW*SW*Zh22*Zuu11*Zuu21*sb
     &      -4D0*EE/CW**2*MW*SW*Zh22*Zuu12*Zuu22*sb-3D0*Lambda
     &      *Mq*Sqrt2*Zh23*Zuu11*Zuu22*cb/sb-3D0*Lambda*Mq
     &      *Sqrt2*Zh23*Zuu12*Zuu21*cb/sb+3D0*Au*EE/MW*Mq
     &      /SW*Zh22*Zuu11*Zuu22/sb+3D0*Au*EE/MW*Mq/SW*Zh22
     &      *Zuu12*Zuu21/sb-3D0*EE*Lambda/MW*Mq/SW*Zh21*Zuu11
     &      *Zuu22/sb*xvev-3D0*EE*Lambda/MW*Mq/SW*Zh21*Zuu12
     &      *Zuu21/sb*xvev
      AAABR(132) = 3D0*EE*MW/SW*Zcc11*Zcc21*Zh21*cb-3D0*EE
     &      *MW/SW*Zcc11*Zcc21*Zh22*sb-EE/CW**2*MW*SW*Zcc11
     &      *Zcc21*Zh21*cb+EE/CW**2*MW*SW*Zcc11*Zcc21*Zh22
     &      *sb+4D0*EE/CW**2*MW*SW*Zcc12*Zcc22*Zh21*cb-4D0
     &      *EE/CW**2*MW*SW*Zcc12*Zcc22*Zh22*sb-3D0*Lambda
     &      *Mc*Sqrt2*Zcc11*Zcc22*Zh23*cb/sb-3D0*Lambda*Mc
     &      *Sqrt2*Zcc12*Zcc21*Zh23*cb/sb+3D0*Au*EE/MW*Mc
     &      /SW*Zcc11*Zcc22*Zh22/sb+3D0*Au*EE/MW*Mc/SW*Zcc12
     &      *Zcc21*Zh22/sb-3D0*EE*Lambda/MW*Mc/SW*Zcc11*Zcc22
     &      *Zh21/sb*xvev-3D0*EE*Lambda/MW*Mc/SW*Zcc12*Zcc21
     &      *Zh21/sb*xvev
      AAABR(133) = 3D0*EE*MW/SW*Zh11*Zt11*Zt21*cb-3D0*EE
     &      *MW/SW*Zh12*Zt11*Zt21*sb-EE/CW**2*MW*SW*Zh11
     &      *Zt11*Zt21*cb+4D0*EE/CW**2*MW*SW*Zh11*Zt12*Zt22
     &      *cb+EE/CW**2*MW*SW*Zh12*Zt11*Zt21*sb-4D0*EE/CW
     &      **2*MW*SW*Zh12*Zt12*Zt22*sb-3D0*Lambda*MtMM*Sqrt2
     &      *Zh13*Zt11*Zt22*cb/sb-3D0*Lambda*MtMM*Sqrt2*Zh13
     &      *Zt12*Zt21*cb/sb+3D0*At*EE/MW*MtMM/SW*Zh12*Zt11
     &      *Zt22/sb+3D0*At*EE/MW*MtMM/SW*Zh12*Zt12*Zt21
     &      /sb-3D0*EE*Lambda/MW*MtMM/SW*Zh11*Zt11*Zt22/sb
     &      *xvev-3D0*EE*Lambda/MW*MtMM/SW*Zh11*Zt12*Zt21
     &      /sb*xvev
      AAABR(134) = 3D0*EE*MW/SW*Zh21*Zt11*Zt21*cb-3D0*EE
     &      *MW/SW*Zh22*Zt11*Zt21*sb-EE/CW**2*MW*SW*Zh21
     &      *Zt11*Zt21*cb+4D0*EE/CW**2*MW*SW*Zh21*Zt12*Zt22
     &      *cb+EE/CW**2*MW*SW*Zh22*Zt11*Zt21*sb-4D0*EE/CW
     &      **2*MW*SW*Zh22*Zt12*Zt22*sb-3D0*Lambda*MtMM*Sqrt2
     &      *Zh23*Zt11*Zt22*cb/sb-3D0*Lambda*MtMM*Sqrt2*Zh23
     &      *Zt12*Zt21*cb/sb+3D0*At*EE/MW*MtMM/SW*Zh22*Zt11
     &      *Zt22/sb+3D0*At*EE/MW*MtMM/SW*Zh22*Zt12*Zt21
     &      /sb-3D0*EE*Lambda/MW*MtMM/SW*Zh21*Zt11*Zt22/sb
     &      *xvev-3D0*EE*Lambda/MW*MtMM/SW*Zh21*Zt12*Zt21
     &      /sb*xvev
      AAABR(135) = 3D0*EE*MW/SW*Zh11*Zuu11*Zuu21*cb-3D0*EE
     &      *MW/SW*Zh12*Zuu11*Zuu21*sb-EE/CW**2*MW*SW*Zh11
     &      *Zuu11*Zuu21*cb+4D0*EE/CW**2*MW*SW*Zh11*Zuu12
     &      *Zuu22*cb+EE/CW**2*MW*SW*Zh12*Zuu11*Zuu21*sb
     &      -4D0*EE/CW**2*MW*SW*Zh12*Zuu12*Zuu22*sb-3D0*Lambda
     &      *Mq*Sqrt2*Zh13*Zuu11*Zuu22*cb/sb-3D0*Lambda*Mq
     &      *Sqrt2*Zh13*Zuu12*Zuu21*cb/sb+3D0*Au*EE/MW*Mq
     &      /SW*Zh12*Zuu11*Zuu22/sb+3D0*Au*EE/MW*Mq/SW*Zh12
     &      *Zuu12*Zuu21/sb-3D0*EE*Lambda/MW*Mq/SW*Zh11*Zuu11
     &      *Zuu22/sb*xvev-3D0*EE*Lambda/MW*Mq/SW*Zh11*Zuu12
     &      *Zuu21/sb*xvev
      AAABR(136) = 3D0*EE*MW/SW*Zcc11*Zcc21*Zh31*cb-3D0*EE
     &      *MW/SW*Zcc11*Zcc21*Zh32*sb-EE/CW**2*MW*SW*Zcc11
     &      *Zcc21*Zh31*cb+EE/CW**2*MW*SW*Zcc11*Zcc21*Zh32
     &      *sb+4D0*EE/CW**2*MW*SW*Zcc12*Zcc22*Zh31*cb-4D0
     &      *EE/CW**2*MW*SW*Zcc12*Zcc22*Zh32*sb-3D0*Lambda
     &      *Mc*Sqrt2*Zcc11*Zcc22*Zh33*cb/sb-3D0*Lambda*Mc
     &      *Sqrt2*Zcc12*Zcc21*Zh33*cb/sb+3D0*Au*EE/MW*Mc
     &      /SW*Zcc11*Zcc22*Zh32/sb+3D0*Au*EE/MW*Mc/SW*Zcc12
     &      *Zcc21*Zh32/sb-3D0*EE*Lambda/MW*Mc/SW*Zcc11*Zcc22
     &      *Zh31/sb*xvev-3D0*EE*Lambda/MW*Mc/SW*Zcc12*Zcc21
     &      *Zh31/sb*xvev
      AAABR(137) = 3D0*EE*MW/SW*Zh31*Zuu11*Zuu21*cb-3D0*EE
     &      *MW/SW*Zh32*Zuu11*Zuu21*sb-EE/CW**2*MW*SW*Zh31
     &      *Zuu11*Zuu21*cb+4D0*EE/CW**2*MW*SW*Zh31*Zuu12
     &      *Zuu22*cb+EE/CW**2*MW*SW*Zh32*Zuu11*Zuu21*sb
     &      -4D0*EE/CW**2*MW*SW*Zh32*Zuu12*Zuu22*sb-3D0*Lambda
     &      *Mq*Sqrt2*Zh33*Zuu11*Zuu22*cb/sb-3D0*Lambda*Mq
     &      *Sqrt2*Zh33*Zuu12*Zuu21*cb/sb+3D0*Au*EE/MW*Mq
     &      /SW*Zh32*Zuu11*Zuu22/sb+3D0*Au*EE/MW*Mq/SW*Zh32
     &      *Zuu12*Zuu21/sb-3D0*EE*Lambda/MW*Mq/SW*Zh31*Zuu11
     &      *Zuu22/sb*xvev-3D0*EE*Lambda/MW*Mq/SW*Zh31*Zuu12
     &      *Zuu21/sb*xvev
      AAABR(138) = 6D0*EE/MW*Mq**2/SW*Zh12/sb+3D0*EE*MW/SW
     &      *Zh11*Zuu21**2*cb-3D0*EE*MW/SW*Zh12*Zuu21**2
     &      *sb-EE/CW**2*MW*SW*Zh11*Zuu21**2*cb+4D0*EE/CW
     &      **2*MW*SW*Zh11*Zuu22**2*cb+EE/CW**2*MW*SW*Zh12
     &      *Zuu21**2*sb-4D0*EE/CW**2*MW*SW*Zh12*Zuu22**2
     &      *sb-6D0*Lambda*Mq*Sqrt2*Zh13*Zuu21*Zuu22*cb/sb
     &      +6D0*Au*EE/MW*Mq/SW*Zh12*Zuu21*Zuu22/sb-6D0*EE
     &      *Lambda/MW*Mq/SW*Zh11*Zuu21*Zuu22/sb*xvev
      AAABR(139) = 6D0*EE/MW*Mq**2/SW*Zh32/sb+3D0*EE*MW/SW
     &      *Zh31*Zuu21**2*cb-3D0*EE*MW/SW*Zh32*Zuu21**2
     &      *sb-EE/CW**2*MW*SW*Zh31*Zuu21**2*cb+4D0*EE/CW
     &      **2*MW*SW*Zh31*Zuu22**2*cb+EE/CW**2*MW*SW*Zh32
     &      *Zuu21**2*sb-4D0*EE/CW**2*MW*SW*Zh32*Zuu22**2
     &      *sb-6D0*Lambda*Mq*Sqrt2*Zh33*Zuu21*Zuu22*cb/sb
     &      +6D0*Au*EE/MW*Mq/SW*Zh32*Zuu21*Zuu22/sb-6D0*EE
     &      *Lambda/MW*Mq/SW*Zh31*Zuu21*Zuu22/sb*xvev
      AAABR(140) = 6D0*EE/MW*Mc**2/SW*Zh12/sb+3D0*EE*MW/SW
     &      *Zcc21**2*Zh11*cb-3D0*EE*MW/SW*Zcc21**2*Zh12
     &      *sb-EE/CW**2*MW*SW*Zcc21**2*Zh11*cb+EE/CW**2
     &      *MW*SW*Zcc21**2*Zh12*sb+4D0*EE/CW**2*MW*SW*Zcc22
     &      **2*Zh11*cb-4D0*EE/CW**2*MW*SW*Zcc22**2*Zh12
     &      *sb-6D0*Lambda*Mc*Sqrt2*Zcc21*Zcc22*Zh13*cb/sb
     &      +6D0*Au*EE/MW*Mc/SW*Zcc21*Zcc22*Zh12/sb-6D0*EE
     &      *Lambda/MW*Mc/SW*Zcc21*Zcc22*Zh11/sb*xvev
      AAABR(141) = 6D0*EE/MW*MtMM**2/SW*Zh22/sb+3D0*EE*MW
     &      /SW*Zh21*Zt21**2*cb-3D0*EE*MW/SW*Zh22*Zt21**2
     &      *sb-EE/CW**2*MW*SW*Zh21*Zt21**2*cb+4D0*EE/CW
     &      **2*MW*SW*Zh21*Zt22**2*cb+EE/CW**2*MW*SW*Zh22
     &      *Zt21**2*sb-4D0*EE/CW**2*MW*SW*Zh22*Zt22**2*sb
     &      -6D0*Lambda*MtMM*Sqrt2*Zh23*Zt21*Zt22*cb/sb+6D0
     &      *At*EE/MW*MtMM/SW*Zh22*Zt21*Zt22/sb-6D0*EE*Lambda
     &      /MW*MtMM/SW*Zh21*Zt21*Zt22/sb*xvev
      AAABR(142) = 6D0*EE/MW*MtMM**2/SW*Zh12/sb+3D0*EE*MW
     &      /SW*Zh11*Zt21**2*cb-3D0*EE*MW/SW*Zh12*Zt21**2
     &      *sb-EE/CW**2*MW*SW*Zh11*Zt21**2*cb+4D0*EE/CW
     &      **2*MW*SW*Zh11*Zt22**2*cb+EE/CW**2*MW*SW*Zh12
     &      *Zt21**2*sb-4D0*EE/CW**2*MW*SW*Zh12*Zt22**2*sb
     &      -6D0*Lambda*MtMM*Sqrt2*Zh13*Zt21*Zt22*cb/sb+6D0
     &      *At*EE/MW*MtMM/SW*Zh12*Zt21*Zt22/sb-6D0*EE*Lambda
     &      /MW*MtMM/SW*Zh11*Zt21*Zt22/sb*xvev
      AAABR(143) = 6D0*EE/MW*Mq**2/SW*Zh22/sb+3D0*EE*MW/SW
     &      *Zh21*Zuu21**2*cb-3D0*EE*MW/SW*Zh22*Zuu21**2
     &      *sb-EE/CW**2*MW*SW*Zh21*Zuu21**2*cb+4D0*EE/CW
     &      **2*MW*SW*Zh21*Zuu22**2*cb+EE/CW**2*MW*SW*Zh22
     &      *Zuu21**2*sb-4D0*EE/CW**2*MW*SW*Zh22*Zuu22**2
     &      *sb-6D0*Lambda*Mq*Sqrt2*Zh23*Zuu21*Zuu22*cb/sb
     &      +6D0*Au*EE/MW*Mq/SW*Zh22*Zuu21*Zuu22/sb-6D0*EE
     &      *Lambda/MW*Mq/SW*Zh21*Zuu21*Zuu22/sb*xvev
      AAABR(144) = 6D0*EE/MW*Mc**2/SW*Zh32/sb+3D0*EE*MW/SW
     &      *Zcc21**2*Zh31*cb-3D0*EE*MW/SW*Zcc21**2*Zh32
     &      *sb-EE/CW**2*MW*SW*Zcc21**2*Zh31*cb+EE/CW**2
     &      *MW*SW*Zcc21**2*Zh32*sb+4D0*EE/CW**2*MW*SW*Zcc22
     &      **2*Zh31*cb-4D0*EE/CW**2*MW*SW*Zcc22**2*Zh32
     &      *sb-6D0*Lambda*Mc*Sqrt2*Zcc21*Zcc22*Zh33*cb/sb
     &      +6D0*Au*EE/MW*Mc/SW*Zcc21*Zcc22*Zh32/sb-6D0*EE
     &      *Lambda/MW*Mc/SW*Zcc21*Zcc22*Zh31/sb*xvev
      AAABR(145) = 6D0*EE/MW*Mc**2/SW*Zh22/sb+3D0*EE*MW/SW
     &      *Zcc21**2*Zh21*cb-3D0*EE*MW/SW*Zcc21**2*Zh22
     &      *sb-EE/CW**2*MW*SW*Zcc21**2*Zh21*cb+EE/CW**2
     &      *MW*SW*Zcc21**2*Zh22*sb+4D0*EE/CW**2*MW*SW*Zcc22
     &      **2*Zh21*cb-4D0*EE/CW**2*MW*SW*Zcc22**2*Zh22
     &      *sb-6D0*Lambda*Mc*Sqrt2*Zcc21*Zcc22*Zh23*cb/sb
     &      +6D0*Au*EE/MW*Mc/SW*Zcc21*Zcc22*Zh22/sb-6D0*EE
     &      *Lambda/MW*Mc/SW*Zcc21*Zcc22*Zh21/sb*xvev
      AAABR(146) = 6D0*EE/MW*MtMM**2/SW*Zh32/sb+3D0*EE*MW
     &      /SW*Zh31*Zt21**2*cb-3D0*EE*MW/SW*Zh32*Zt21**2
     &      *sb-EE/CW**2*MW*SW*Zh31*Zt21**2*cb+4D0*EE/CW
     &      **2*MW*SW*Zh31*Zt22**2*cb+EE/CW**2*MW*SW*Zh32
     &      *Zt21**2*sb-4D0*EE/CW**2*MW*SW*Zh32*Zt22**2*sb
     &      -6D0*Lambda*MtMM*Sqrt2*Zh33*Zt21*Zt22*cb/sb+6D0
     &      *At*EE/MW*MtMM/SW*Zh32*Zt21*Zt22/sb-6D0*EE*Lambda
     &      /MW*MtMM/SW*Zh31*Zt21*Zt22/sb*xvev
      AAABR(147) = Sqrt2*Zh23*aa1-Sqrt2*Zh23*aa1*sb**2+Sqrt2
     &      *Zh23*aa2*sb**2+Sqrt2*Zh23*la1s*xvev-Sqrt2*Zh23
     &      *aa5*cb*sb-Sqrt2*Zh23*aa6*cb*sb-Sqrt2*Zh23*la1s
     &      *sb**2*xvev+Sqrt2*Zh23*la2s*sb**2*xvev+2D0*MW
     &      /EE*SW*Zh21*cb*la1-6D0*MW/EE*SW*Zh21*la6*sb+6D0
     &      *MW/EE*SW*Zh21*la6*sb**3-2D0*MW/EE*SW*Zh21*la7
     &      *sb**3-2D0*MW/EE*SW*Zh22*cb*la6+2D0*MW/EE*SW
     &      *Zh22*la2*sb**3+2D0*MW/EE*SW*Zh22*la3*sb-2D0
     &      *MW/EE*SW*Zh22*la3*sb**3+2D0*MW/EE*SW*Zh22*la4
     &      *sb-2D0*MW/EE*SW*Zh22*la4*sb**3+2D0*MW/EE*SW
     &      *Zh22*la5*sb-2D0*MW/EE*SW*Zh22*la5*sb**3-2D0
     &      *Sqrt2*Zh23*cb*la4s*sb*xvev-2D0*Sqrt2*Zh23*cb
     &      *la5s*sb*xvev-2D0*Sqrt2*Zh23*cb*la6s*sb*xvev
     &      -2D0*MW/EE*SW*Zh21*cb*la1*sb**2+2D0*MW/EE*SW
     &      *Zh21*cb*la3*sb**2+2D0*MW/EE*SW*Zh21*cb*la4*sb
     &      **2+2D0*MW/EE*SW*Zh21*cb*la5*sb**2+2D0*MW/EE
     &      *SW*Zh22*cb*la6*sb**2-6D0*MW/EE*SW*Zh22*cb*la7
     &      *sb**2
      AAABR(148) = Sqrt2*Zh33*aa1-Sqrt2*Zh33*aa1*sb**2+Sqrt2
     &      *Zh33*aa2*sb**2+Sqrt2*Zh33*la1s*xvev-Sqrt2*Zh33
     &      *aa5*cb*sb-Sqrt2*Zh33*aa6*cb*sb-Sqrt2*Zh33*la1s
     &      *sb**2*xvev+Sqrt2*Zh33*la2s*sb**2*xvev+2D0*MW
     &      /EE*SW*Zh31*cb*la1-6D0*MW/EE*SW*Zh31*la6*sb+6D0
     &      *MW/EE*SW*Zh31*la6*sb**3-2D0*MW/EE*SW*Zh31*la7
     &      *sb**3-2D0*MW/EE*SW*Zh32*cb*la6+2D0*MW/EE*SW
     &      *Zh32*la2*sb**3+2D0*MW/EE*SW*Zh32*la3*sb-2D0
     &      *MW/EE*SW*Zh32*la3*sb**3+2D0*MW/EE*SW*Zh32*la4
     &      *sb-2D0*MW/EE*SW*Zh32*la4*sb**3+2D0*MW/EE*SW
     &      *Zh32*la5*sb-2D0*MW/EE*SW*Zh32*la5*sb**3-2D0
     &      *Sqrt2*Zh33*cb*la4s*sb*xvev-2D0*Sqrt2*Zh33*cb
     &      *la5s*sb*xvev-2D0*Sqrt2*Zh33*cb*la6s*sb*xvev
     &      -2D0*MW/EE*SW*Zh31*cb*la1*sb**2+2D0*MW/EE*SW
     &      *Zh31*cb*la3*sb**2+2D0*MW/EE*SW*Zh31*cb*la4*sb
     &      **2+2D0*MW/EE*SW*Zh31*cb*la5*sb**2+2D0*MW/EE
     &      *SW*Zh32*cb*la6*sb**2-6D0*MW/EE*SW*Zh32*cb*la7
     &      *sb**2
      AAABR(149) = Sqrt2*Zh13*aa1-Sqrt2*Zh13*aa1*sb**2+Sqrt2
     &      *Zh13*aa2*sb**2+Sqrt2*Zh13*la1s*xvev-Sqrt2*Zh13
     &      *aa5*cb*sb-Sqrt2*Zh13*aa6*cb*sb-Sqrt2*Zh13*la1s
     &      *sb**2*xvev+Sqrt2*Zh13*la2s*sb**2*xvev+2D0*MW
     &      /EE*SW*Zh11*cb*la1-6D0*MW/EE*SW*Zh11*la6*sb+6D0
     &      *MW/EE*SW*Zh11*la6*sb**3-2D0*MW/EE*SW*Zh11*la7
     &      *sb**3-2D0*MW/EE*SW*Zh12*cb*la6+2D0*MW/EE*SW
     &      *Zh12*la2*sb**3+2D0*MW/EE*SW*Zh12*la3*sb-2D0
     &      *MW/EE*SW*Zh12*la3*sb**3+2D0*MW/EE*SW*Zh12*la4
     &      *sb-2D0*MW/EE*SW*Zh12*la4*sb**3+2D0*MW/EE*SW
     &      *Zh12*la5*sb-2D0*MW/EE*SW*Zh12*la5*sb**3-2D0
     &      *Sqrt2*Zh13*cb*la4s*sb*xvev-2D0*Sqrt2*Zh13*cb
     &      *la5s*sb*xvev-2D0*Sqrt2*Zh13*cb*la6s*sb*xvev
     &      -2D0*MW/EE*SW*Zh11*cb*la1*sb**2+2D0*MW/EE*SW
     &      *Zh11*cb*la3*sb**2+2D0*MW/EE*SW*Zh11*cb*la4*sb
     &      **2+2D0*MW/EE*SW*Zh11*cb*la5*sb**2+2D0*MW/EE
     &      *SW*Zh12*cb*la6*sb**2-6D0*MW/EE*SW*Zh12*cb*la7
     &      *sb**2
      AAABR(150) = Pa21*Sqrt2*Zh23*aa5+Pa21*Sqrt2*Zh23*aa6
     &      -2D0*Pa21*Sqrt2*Zh23*aa5*sb**2-2D0*Pa21*Sqrt2
     &      *Zh23*aa6*sb**2+2D0*Pa21*Sqrt2*Zh23*la4s*xvev
     &      +2D0*Pa21*Sqrt2*Zh23*la5s*xvev+2D0*Pa21*Sqrt2
     &      *Zh23*la6s*xvev-Pa22*Sqrt2*Zh21*aa5*sb+Pa22*Sqrt2
     &      *Zh21*aa6*sb+Pa22*Sqrt2*Zh22*aa5*cb-Pa22*Sqrt2
     &      *Zh22*aa6*cb+2D0*Pa21*Sqrt2*Zh23*aa1*cb*sb-2D0
     &      *Pa21*Sqrt2*Zh23*aa2*cb*sb-4D0*Pa21*Sqrt2*Zh23
     &      *la4s*sb**2*xvev-4D0*Pa21*Sqrt2*Zh23*la5s*sb
     &      **2*xvev-4D0*Pa21*Sqrt2*Zh23*la6s*sb**2*xvev
     &      -2D0*Pa22*Sqrt2*Zh21*la5s*sb*xvev+2D0*Pa22*Sqrt2
     &      *Zh21*la6s*sb*xvev+2D0*Pa22*Sqrt2*Zh22*cb*la5s
     &      *xvev-2D0*Pa22*Sqrt2*Zh22*cb*la6s*xvev+4D0*MW
     &      /EE*Pa21*SW*Zh21*cb*la6+4D0*MW/EE*Pa21*SW*Zh21
     &      *la1*sb-4D0*MW/EE*Pa21*SW*Zh21*la1*sb**3-4D0
     &      *MW/EE*Pa21*SW*Zh21*la3*sb+4D0*MW/EE*Pa21*SW
     &      *Zh21*la3*sb**3-4D0*MW/EE*Pa21*SW*Zh21*la4*sb
     &      +4D0*MW/EE*Pa21*SW*Zh21*la4*sb**3+4D0*MW/EE*Pa21
     &      *SW*Zh21*la5*sb**3-4D0*MW/EE*Pa21*SW*Zh22*cb
     &      *la5-4D0*MW/EE*Pa21*SW*Zh22*la6*sb+4D0*MW/EE
     &      *Pa21*SW*Zh22*la6*sb**3+8D0*MW/EE*Pa21*SW*Zh22
     &      *la7*sb-12D0*MW/EE*Pa21*SW*Zh22*la7*sb**3+2D0
     &      *Pa21*Sqrt2*Zh23*cb*la1s*sb*xvev-2D0*Pa21*Sqrt2
     &      *Zh23*cb*la2s*sb*xvev-12D0*MW/EE*Pa21*SW*Zh21
     &      *cb*la6*sb**2+4D0*MW/EE*Pa21*SW*Zh21*cb*la7*sb
     &      **2-4D0*MW/EE*Pa21*SW*Zh22*cb*la2*sb**2+4D0*MW
     &      /EE*Pa21*SW*Zh22*cb*la3*sb**2+4D0*MW/EE*Pa21
     &      *SW*Zh22*cb*la4*sb**2+4D0*MW/EE*Pa21*SW*Zh22
     &      *cb*la5*sb**2
      AAABR(151) = Pa11*Sqrt2*Zh13*aa5+Pa11*Sqrt2*Zh13*aa6
     &      -2D0*Pa11*Sqrt2*Zh13*aa5*sb**2-2D0*Pa11*Sqrt2
     &      *Zh13*aa6*sb**2+2D0*Pa11*Sqrt2*Zh13*la4s*xvev
     &      +2D0*Pa11*Sqrt2*Zh13*la5s*xvev+2D0*Pa11*Sqrt2
     &      *Zh13*la6s*xvev-Pa12*Sqrt2*Zh11*aa5*sb+Pa12*Sqrt2
     &      *Zh11*aa6*sb+Pa12*Sqrt2*Zh12*aa5*cb-Pa12*Sqrt2
     &      *Zh12*aa6*cb+2D0*Pa11*Sqrt2*Zh13*aa1*cb*sb-2D0
     &      *Pa11*Sqrt2*Zh13*aa2*cb*sb-4D0*Pa11*Sqrt2*Zh13
     &      *la4s*sb**2*xvev-4D0*Pa11*Sqrt2*Zh13*la5s*sb
     &      **2*xvev-4D0*Pa11*Sqrt2*Zh13*la6s*sb**2*xvev
     &      -2D0*Pa12*Sqrt2*Zh11*la5s*sb*xvev+2D0*Pa12*Sqrt2
     &      *Zh11*la6s*sb*xvev+2D0*Pa12*Sqrt2*Zh12*cb*la5s
     &      *xvev-2D0*Pa12*Sqrt2*Zh12*cb*la6s*xvev+4D0*MW
     &      /EE*Pa11*SW*Zh11*cb*la6+4D0*MW/EE*Pa11*SW*Zh11
     &      *la1*sb-4D0*MW/EE*Pa11*SW*Zh11*la1*sb**3-4D0
     &      *MW/EE*Pa11*SW*Zh11*la3*sb+4D0*MW/EE*Pa11*SW
     &      *Zh11*la3*sb**3-4D0*MW/EE*Pa11*SW*Zh11*la4*sb
     &      +4D0*MW/EE*Pa11*SW*Zh11*la4*sb**3+4D0*MW/EE*Pa11
     &      *SW*Zh11*la5*sb**3-4D0*MW/EE*Pa11*SW*Zh12*cb
     &      *la5-4D0*MW/EE*Pa11*SW*Zh12*la6*sb+4D0*MW/EE
     &      *Pa11*SW*Zh12*la6*sb**3+8D0*MW/EE*Pa11*SW*Zh12
     &      *la7*sb-12D0*MW/EE*Pa11*SW*Zh12*la7*sb**3+2D0
     &      *Pa11*Sqrt2*Zh13*cb*la1s*sb*xvev-2D0*Pa11*Sqrt2
     &      *Zh13*cb*la2s*sb*xvev-12D0*MW/EE*Pa11*SW*Zh11
     &      *cb*la6*sb**2+4D0*MW/EE*Pa11*SW*Zh11*cb*la7*sb
     &      **2-4D0*MW/EE*Pa11*SW*Zh12*cb*la2*sb**2+4D0*MW
     &      /EE*Pa11*SW*Zh12*cb*la3*sb**2+4D0*MW/EE*Pa11
     &      *SW*Zh12*cb*la4*sb**2+4D0*MW/EE*Pa11*SW*Zh12
     &      *cb*la5*sb**2
      AAABR(152) = Pa21*Sqrt2*Zh33*aa5+Pa21*Sqrt2*Zh33*aa6
     &      -2D0*Pa21*Sqrt2*Zh33*aa5*sb**2-2D0*Pa21*Sqrt2
     &      *Zh33*aa6*sb**2+2D0*Pa21*Sqrt2*Zh33*la4s*xvev
     &      +2D0*Pa21*Sqrt2*Zh33*la5s*xvev+2D0*Pa21*Sqrt2
     &      *Zh33*la6s*xvev-Pa22*Sqrt2*Zh31*aa5*sb+Pa22*Sqrt2
     &      *Zh31*aa6*sb+Pa22*Sqrt2*Zh32*aa5*cb-Pa22*Sqrt2
     &      *Zh32*aa6*cb+2D0*Pa21*Sqrt2*Zh33*aa1*cb*sb-2D0
     &      *Pa21*Sqrt2*Zh33*aa2*cb*sb-4D0*Pa21*Sqrt2*Zh33
     &      *la4s*sb**2*xvev-4D0*Pa21*Sqrt2*Zh33*la5s*sb
     &      **2*xvev-4D0*Pa21*Sqrt2*Zh33*la6s*sb**2*xvev
     &      -2D0*Pa22*Sqrt2*Zh31*la5s*sb*xvev+2D0*Pa22*Sqrt2
     &      *Zh31*la6s*sb*xvev+2D0*Pa22*Sqrt2*Zh32*cb*la5s
     &      *xvev-2D0*Pa22*Sqrt2*Zh32*cb*la6s*xvev+4D0*MW
     &      /EE*Pa21*SW*Zh31*cb*la6+4D0*MW/EE*Pa21*SW*Zh31
     &      *la1*sb-4D0*MW/EE*Pa21*SW*Zh31*la1*sb**3-4D0
     &      *MW/EE*Pa21*SW*Zh31*la3*sb+4D0*MW/EE*Pa21*SW
     &      *Zh31*la3*sb**3-4D0*MW/EE*Pa21*SW*Zh31*la4*sb
     &      +4D0*MW/EE*Pa21*SW*Zh31*la4*sb**3+4D0*MW/EE*Pa21
     &      *SW*Zh31*la5*sb**3-4D0*MW/EE*Pa21*SW*Zh32*cb
     &      *la5-4D0*MW/EE*Pa21*SW*Zh32*la6*sb+4D0*MW/EE
     &      *Pa21*SW*Zh32*la6*sb**3+8D0*MW/EE*Pa21*SW*Zh32
     &      *la7*sb-12D0*MW/EE*Pa21*SW*Zh32*la7*sb**3+2D0
     &      *Pa21*Sqrt2*Zh33*cb*la1s*sb*xvev-2D0*Pa21*Sqrt2
     &      *Zh33*cb*la2s*sb*xvev-12D0*MW/EE*Pa21*SW*Zh31
     &      *cb*la6*sb**2+4D0*MW/EE*Pa21*SW*Zh31*cb*la7*sb
     &      **2-4D0*MW/EE*Pa21*SW*Zh32*cb*la2*sb**2+4D0*MW
     &      /EE*Pa21*SW*Zh32*cb*la3*sb**2+4D0*MW/EE*Pa21
     &      *SW*Zh32*cb*la4*sb**2+4D0*MW/EE*Pa21*SW*Zh32
     &      *cb*la5*sb**2
      AAABR(153) = Pa11*Sqrt2*Zh23*aa5+Pa11*Sqrt2*Zh23*aa6
     &      -2D0*Pa11*Sqrt2*Zh23*aa5*sb**2-2D0*Pa11*Sqrt2
     &      *Zh23*aa6*sb**2+2D0*Pa11*Sqrt2*Zh23*la4s*xvev
     &      +2D0*Pa11*Sqrt2*Zh23*la5s*xvev+2D0*Pa11*Sqrt2
     &      *Zh23*la6s*xvev-Pa12*Sqrt2*Zh21*aa5*sb+Pa12*Sqrt2
     &      *Zh21*aa6*sb+Pa12*Sqrt2*Zh22*aa5*cb-Pa12*Sqrt2
     &      *Zh22*aa6*cb+2D0*Pa11*Sqrt2*Zh23*aa1*cb*sb-2D0
     &      *Pa11*Sqrt2*Zh23*aa2*cb*sb-4D0*Pa11*Sqrt2*Zh23
     &      *la4s*sb**2*xvev-4D0*Pa11*Sqrt2*Zh23*la5s*sb
     &      **2*xvev-4D0*Pa11*Sqrt2*Zh23*la6s*sb**2*xvev
     &      -2D0*Pa12*Sqrt2*Zh21*la5s*sb*xvev+2D0*Pa12*Sqrt2
     &      *Zh21*la6s*sb*xvev+2D0*Pa12*Sqrt2*Zh22*cb*la5s
     &      *xvev-2D0*Pa12*Sqrt2*Zh22*cb*la6s*xvev+4D0*MW
     &      /EE*Pa11*SW*Zh21*cb*la6+4D0*MW/EE*Pa11*SW*Zh21
     &      *la1*sb-4D0*MW/EE*Pa11*SW*Zh21*la1*sb**3-4D0
     &      *MW/EE*Pa11*SW*Zh21*la3*sb+4D0*MW/EE*Pa11*SW
     &      *Zh21*la3*sb**3-4D0*MW/EE*Pa11*SW*Zh21*la4*sb
     &      +4D0*MW/EE*Pa11*SW*Zh21*la4*sb**3+4D0*MW/EE*Pa11
     &      *SW*Zh21*la5*sb**3-4D0*MW/EE*Pa11*SW*Zh22*cb
     &      *la5-4D0*MW/EE*Pa11*SW*Zh22*la6*sb+4D0*MW/EE
     &      *Pa11*SW*Zh22*la6*sb**3+8D0*MW/EE*Pa11*SW*Zh22
     &      *la7*sb-12D0*MW/EE*Pa11*SW*Zh22*la7*sb**3+2D0
     &      *Pa11*Sqrt2*Zh23*cb*la1s*sb*xvev-2D0*Pa11*Sqrt2
     &      *Zh23*cb*la2s*sb*xvev-12D0*MW/EE*Pa11*SW*Zh21
     &      *cb*la6*sb**2+4D0*MW/EE*Pa11*SW*Zh21*cb*la7*sb
     &      **2-4D0*MW/EE*Pa11*SW*Zh22*cb*la2*sb**2+4D0*MW
     &      /EE*Pa11*SW*Zh22*cb*la3*sb**2+4D0*MW/EE*Pa11
     &      *SW*Zh22*cb*la4*sb**2+4D0*MW/EE*Pa11*SW*Zh22
     &      *cb*la5*sb**2
      AAABR(154) = Pa21*Sqrt2*Zh13*aa5+Pa21*Sqrt2*Zh13*aa6
     &      -2D0*Pa21*Sqrt2*Zh13*aa5*sb**2-2D0*Pa21*Sqrt2
     &      *Zh13*aa6*sb**2+2D0*Pa21*Sqrt2*Zh13*la4s*xvev
     &      +2D0*Pa21*Sqrt2*Zh13*la5s*xvev+2D0*Pa21*Sqrt2
     &      *Zh13*la6s*xvev-Pa22*Sqrt2*Zh11*aa5*sb+Pa22*Sqrt2
     &      *Zh11*aa6*sb+Pa22*Sqrt2*Zh12*aa5*cb-Pa22*Sqrt2
     &      *Zh12*aa6*cb+2D0*Pa21*Sqrt2*Zh13*aa1*cb*sb-2D0
     &      *Pa21*Sqrt2*Zh13*aa2*cb*sb-4D0*Pa21*Sqrt2*Zh13
     &      *la4s*sb**2*xvev-4D0*Pa21*Sqrt2*Zh13*la5s*sb
     &      **2*xvev-4D0*Pa21*Sqrt2*Zh13*la6s*sb**2*xvev
     &      -2D0*Pa22*Sqrt2*Zh11*la5s*sb*xvev+2D0*Pa22*Sqrt2
     &      *Zh11*la6s*sb*xvev+2D0*Pa22*Sqrt2*Zh12*cb*la5s
     &      *xvev-2D0*Pa22*Sqrt2*Zh12*cb*la6s*xvev+4D0*MW
     &      /EE*Pa21*SW*Zh11*cb*la6+4D0*MW/EE*Pa21*SW*Zh11
     &      *la1*sb-4D0*MW/EE*Pa21*SW*Zh11*la1*sb**3-4D0
     &      *MW/EE*Pa21*SW*Zh11*la3*sb+4D0*MW/EE*Pa21*SW
     &      *Zh11*la3*sb**3-4D0*MW/EE*Pa21*SW*Zh11*la4*sb
     &      +4D0*MW/EE*Pa21*SW*Zh11*la4*sb**3+4D0*MW/EE*Pa21
     &      *SW*Zh11*la5*sb**3-4D0*MW/EE*Pa21*SW*Zh12*cb
     &      *la5-4D0*MW/EE*Pa21*SW*Zh12*la6*sb+4D0*MW/EE
     &      *Pa21*SW*Zh12*la6*sb**3+8D0*MW/EE*Pa21*SW*Zh12
     &      *la7*sb-12D0*MW/EE*Pa21*SW*Zh12*la7*sb**3+2D0
     &      *Pa21*Sqrt2*Zh13*cb*la1s*sb*xvev-2D0*Pa21*Sqrt2
     &      *Zh13*cb*la2s*sb*xvev-12D0*MW/EE*Pa21*SW*Zh11
     &      *cb*la6*sb**2+4D0*MW/EE*Pa21*SW*Zh11*cb*la7*sb
     &      **2-4D0*MW/EE*Pa21*SW*Zh12*cb*la2*sb**2+4D0*MW
     &      /EE*Pa21*SW*Zh12*cb*la3*sb**2+4D0*MW/EE*Pa21
     &      *SW*Zh12*cb*la4*sb**2+4D0*MW/EE*Pa21*SW*Zh12
     &      *cb*la5*sb**2
      AAABR(155) = Pa11*Sqrt2*Zh33*aa5+Pa11*Sqrt2*Zh33*aa6
     &      -2D0*Pa11*Sqrt2*Zh33*aa5*sb**2-2D0*Pa11*Sqrt2
     &      *Zh33*aa6*sb**2+2D0*Pa11*Sqrt2*Zh33*la4s*xvev
     &      +2D0*Pa11*Sqrt2*Zh33*la5s*xvev+2D0*Pa11*Sqrt2
     &      *Zh33*la6s*xvev-Pa12*Sqrt2*Zh31*aa5*sb+Pa12*Sqrt2
     &      *Zh31*aa6*sb+Pa12*Sqrt2*Zh32*aa5*cb-Pa12*Sqrt2
     &      *Zh32*aa6*cb+2D0*Pa11*Sqrt2*Zh33*aa1*cb*sb-2D0
     &      *Pa11*Sqrt2*Zh33*aa2*cb*sb-4D0*Pa11*Sqrt2*Zh33
     &      *la4s*sb**2*xvev-4D0*Pa11*Sqrt2*Zh33*la5s*sb
     &      **2*xvev-4D0*Pa11*Sqrt2*Zh33*la6s*sb**2*xvev
     &      -2D0*Pa12*Sqrt2*Zh31*la5s*sb*xvev+2D0*Pa12*Sqrt2
     &      *Zh31*la6s*sb*xvev+2D0*Pa12*Sqrt2*Zh32*cb*la5s
     &      *xvev-2D0*Pa12*Sqrt2*Zh32*cb*la6s*xvev+4D0*MW
     &      /EE*Pa11*SW*Zh31*cb*la6+4D0*MW/EE*Pa11*SW*Zh31
     &      *la1*sb-4D0*MW/EE*Pa11*SW*Zh31*la1*sb**3-4D0
     &      *MW/EE*Pa11*SW*Zh31*la3*sb+4D0*MW/EE*Pa11*SW
     &      *Zh31*la3*sb**3-4D0*MW/EE*Pa11*SW*Zh31*la4*sb
     &      +4D0*MW/EE*Pa11*SW*Zh31*la4*sb**3+4D0*MW/EE*Pa11
     &      *SW*Zh31*la5*sb**3-4D0*MW/EE*Pa11*SW*Zh32*cb
     &      *la5-4D0*MW/EE*Pa11*SW*Zh32*la6*sb+4D0*MW/EE
     &      *Pa11*SW*Zh32*la6*sb**3+8D0*MW/EE*Pa11*SW*Zh32
     &      *la7*sb-12D0*MW/EE*Pa11*SW*Zh32*la7*sb**3+2D0
     &      *Pa11*Sqrt2*Zh33*cb*la1s*sb*xvev-2D0*Pa11*Sqrt2
     &      *Zh33*cb*la2s*sb*xvev-12D0*MW/EE*Pa11*SW*Zh31
     &      *cb*la6*sb**2+4D0*MW/EE*Pa11*SW*Zh31*cb*la7*sb
     &      **2-4D0*MW/EE*Pa11*SW*Zh32*cb*la2*sb**2+4D0*MW
     &      /EE*Pa11*SW*Zh32*cb*la3*sb**2+4D0*MW/EE*Pa11
     &      *SW*Zh32*cb*la4*sb**2+4D0*MW/EE*Pa11*SW*Zh32
     &      *cb*la5*sb**2
      AAABR(156) = Pa11**2*Sqrt2*Zh23*aa2+Pa12**2*Sqrt2*Zh23
     &      *aa3-3D0*Pa12**2*Sqrt2*Zh23*aa4+Pa11**2*Sqrt2
     &      *Zh23*aa1*sb**2-Pa11**2*Sqrt2*Zh23*aa2*sb**2
     &      +Pa11**2*Sqrt2*Zh23*la2s*xvev+2D0*Pa12**2*Sqrt2
     &      *Zh23*la3s*xvev-12D0*Pa12**2*Sqrt2*Zh23*la7s
     &      *xvev+Pa11*Pa12*Sqrt2*Zh21*aa5*cb-Pa11*Pa12*Sqrt2
     &      *Zh21*aa6*cb+Pa11*Pa12*Sqrt2*Zh22*aa5*sb-Pa11
     &      *Pa12*Sqrt2*Zh22*aa6*sb+Pa11**2*Sqrt2*Zh23*aa5
     &      *cb*sb+Pa11**2*Sqrt2*Zh23*aa6*cb*sb+Pa11**2*Sqrt2
     &      *Zh23*la1s*sb**2*xvev-Pa11**2*Sqrt2*Zh23*la2s
     &      *sb**2*xvev+4D0*MW/EE*Pa11*Pa12*SW*Zh23*la5s
     &      -4D0*MW/EE*Pa11*Pa12*SW*Zh23*la6s+2D0*MW/EE*Pa11
     &      **2*SW*Zh21*cb*la3+2D0*MW/EE*Pa11**2*SW*Zh21
     &      *cb*la4-2D0*MW/EE*Pa11**2*SW*Zh21*cb*la5+4D0
     &      *MW/EE*Pa11**2*SW*Zh21*la6*sb-6D0*MW/EE*Pa11
     &      **2*SW*Zh21*la6*sb**3-2D0*MW/EE*Pa11**2*SW*Zh21
     &      *la7*sb+2D0*MW/EE*Pa11**2*SW*Zh21*la7*sb**3-2D0
     &      *MW/EE*Pa11**2*SW*Zh22*cb*la7+2D0*MW/EE*Pa11
     &      **2*SW*Zh22*la2*sb-2D0*MW/EE*Pa11**2*SW*Zh22
     &      *la2*sb**3+2D0*MW/EE*Pa11**2*SW*Zh22*la3*sb**3
     &      +2D0*MW/EE*Pa11**2*SW*Zh22*la4*sb**3-4D0*MW/EE
     &      *Pa11**2*SW*Zh22*la5*sb+2D0*MW/EE*Pa11**2*SW
     &      *Zh22*la5*sb**3+2D0*MW/EE*Pa12**2*SW*Zh21*cb
     &      *la1s-2D0*MW/EE*Pa12**2*SW*Zh21*la4s*sb+2D0*MW
     &      /EE*Pa12**2*SW*Zh21*la5s*sb+2D0*MW/EE*Pa12**2
     &      *SW*Zh21*la6s*sb-2D0*MW/EE*Pa12**2*SW*Zh22*cb
     &      *la4s+2D0*MW/EE*Pa12**2*SW*Zh22*cb*la5s+2D0*MW
     &      /EE*Pa12**2*SW*Zh22*cb*la6s+2D0*MW/EE*Pa12**2
     &      *SW*Zh22*la2s*sb+2D0*Pa11*Pa12*Sqrt2*Zh21*cb
     &      *la5s*xvev-2D0*Pa11*Pa12*Sqrt2*Zh21*cb*la6s*xvev
     &      +2D0*Pa11*Pa12*Sqrt2*Zh22*la5s*sb*xvev-2D0*Pa11
     &      *Pa12*Sqrt2*Zh22*la6s*sb*xvev+2D0*Pa11**2*Sqrt2
     &      *Zh23*cb*la4s*sb*xvev+2D0*Pa11**2*Sqrt2*Zh23
     &      *cb*la5s*sb*xvev+2D0*Pa11**2*Sqrt2*Zh23*cb*la6s
     &      *sb*xvev+2D0*MW/EE*Pa11**2*SW*Zh21*cb*la1*sb
     &      **2-2D0*MW/EE*Pa11**2*SW*Zh21*cb*la3*sb**2-2D0
     &      *MW/EE*Pa11**2*SW*Zh21*cb*la4*sb**2-2D0*MW/EE
     &      *Pa11**2*SW*Zh21*cb*la5*sb**2-2D0*MW/EE*Pa11
     &      **2*SW*Zh22*cb*la6*sb**2+6D0*MW/EE*Pa11**2*SW
     &      *Zh22*cb*la7*sb**2
      AAABR(157) = Pa11**2*Sqrt2*Zh33*aa2+Pa12**2*Sqrt2*Zh33
     &      *aa3-3D0*Pa12**2*Sqrt2*Zh33*aa4+Pa11**2*Sqrt2
     &      *Zh33*aa1*sb**2-Pa11**2*Sqrt2*Zh33*aa2*sb**2
     &      +Pa11**2*Sqrt2*Zh33*la2s*xvev+2D0*Pa12**2*Sqrt2
     &      *Zh33*la3s*xvev-12D0*Pa12**2*Sqrt2*Zh33*la7s
     &      *xvev+Pa11*Pa12*Sqrt2*Zh31*aa5*cb-Pa11*Pa12*Sqrt2
     &      *Zh31*aa6*cb+Pa11*Pa12*Sqrt2*Zh32*aa5*sb-Pa11
     &      *Pa12*Sqrt2*Zh32*aa6*sb+Pa11**2*Sqrt2*Zh33*aa5
     &      *cb*sb+Pa11**2*Sqrt2*Zh33*aa6*cb*sb+Pa11**2*Sqrt2
     &      *Zh33*la1s*sb**2*xvev-Pa11**2*Sqrt2*Zh33*la2s
     &      *sb**2*xvev+4D0*MW/EE*Pa11*Pa12*SW*Zh33*la5s
     &      -4D0*MW/EE*Pa11*Pa12*SW*Zh33*la6s+2D0*MW/EE*Pa11
     &      **2*SW*Zh31*cb*la3+2D0*MW/EE*Pa11**2*SW*Zh31
     &      *cb*la4-2D0*MW/EE*Pa11**2*SW*Zh31*cb*la5+4D0
     &      *MW/EE*Pa11**2*SW*Zh31*la6*sb-6D0*MW/EE*Pa11
     &      **2*SW*Zh31*la6*sb**3-2D0*MW/EE*Pa11**2*SW*Zh31
     &      *la7*sb+2D0*MW/EE*Pa11**2*SW*Zh31*la7*sb**3-2D0
     &      *MW/EE*Pa11**2*SW*Zh32*cb*la7+2D0*MW/EE*Pa11
     &      **2*SW*Zh32*la2*sb-2D0*MW/EE*Pa11**2*SW*Zh32
     &      *la2*sb**3+2D0*MW/EE*Pa11**2*SW*Zh32*la3*sb**3
     &      +2D0*MW/EE*Pa11**2*SW*Zh32*la4*sb**3-4D0*MW/EE
     &      *Pa11**2*SW*Zh32*la5*sb+2D0*MW/EE*Pa11**2*SW
     &      *Zh32*la5*sb**3+2D0*MW/EE*Pa12**2*SW*Zh31*cb
     &      *la1s-2D0*MW/EE*Pa12**2*SW*Zh31*la4s*sb+2D0*MW
     &      /EE*Pa12**2*SW*Zh31*la5s*sb+2D0*MW/EE*Pa12**2
     &      *SW*Zh31*la6s*sb-2D0*MW/EE*Pa12**2*SW*Zh32*cb
     &      *la4s+2D0*MW/EE*Pa12**2*SW*Zh32*cb*la5s+2D0*MW
     &      /EE*Pa12**2*SW*Zh32*cb*la6s+2D0*MW/EE*Pa12**2
     &      *SW*Zh32*la2s*sb+2D0*Pa11*Pa12*Sqrt2*Zh31*cb
     &      *la5s*xvev-2D0*Pa11*Pa12*Sqrt2*Zh31*cb*la6s*xvev
     &      +2D0*Pa11*Pa12*Sqrt2*Zh32*la5s*sb*xvev-2D0*Pa11
     &      *Pa12*Sqrt2*Zh32*la6s*sb*xvev+2D0*Pa11**2*Sqrt2
     &      *Zh33*cb*la4s*sb*xvev+2D0*Pa11**2*Sqrt2*Zh33
     &      *cb*la5s*sb*xvev+2D0*Pa11**2*Sqrt2*Zh33*cb*la6s
     &      *sb*xvev+2D0*MW/EE*Pa11**2*SW*Zh31*cb*la1*sb
     &      **2-2D0*MW/EE*Pa11**2*SW*Zh31*cb*la3*sb**2-2D0
     &      *MW/EE*Pa11**2*SW*Zh31*cb*la4*sb**2-2D0*MW/EE
     &      *Pa11**2*SW*Zh31*cb*la5*sb**2-2D0*MW/EE*Pa11
     &      **2*SW*Zh32*cb*la6*sb**2+6D0*MW/EE*Pa11**2*SW
     &      *Zh32*cb*la7*sb**2
      AAABR(158) = Pa21**2*Sqrt2*Zh13*aa2+Pa22**2*Sqrt2*Zh13
     &      *aa3-3D0*Pa22**2*Sqrt2*Zh13*aa4+Pa21**2*Sqrt2
     &      *Zh13*aa1*sb**2-Pa21**2*Sqrt2*Zh13*aa2*sb**2
     &      +Pa21**2*Sqrt2*Zh13*la2s*xvev+2D0*Pa22**2*Sqrt2
     &      *Zh13*la3s*xvev-12D0*Pa22**2*Sqrt2*Zh13*la7s
     &      *xvev+Pa21*Pa22*Sqrt2*Zh11*aa5*cb-Pa21*Pa22*Sqrt2
     &      *Zh11*aa6*cb+Pa21*Pa22*Sqrt2*Zh12*aa5*sb-Pa21
     &      *Pa22*Sqrt2*Zh12*aa6*sb+Pa21**2*Sqrt2*Zh13*aa5
     &      *cb*sb+Pa21**2*Sqrt2*Zh13*aa6*cb*sb+Pa21**2*Sqrt2
     &      *Zh13*la1s*sb**2*xvev-Pa21**2*Sqrt2*Zh13*la2s
     &      *sb**2*xvev+4D0*MW/EE*Pa21*Pa22*SW*Zh13*la5s
     &      -4D0*MW/EE*Pa21*Pa22*SW*Zh13*la6s+2D0*MW/EE*Pa21
     &      **2*SW*Zh11*cb*la3+2D0*MW/EE*Pa21**2*SW*Zh11
     &      *cb*la4-2D0*MW/EE*Pa21**2*SW*Zh11*cb*la5+4D0
     &      *MW/EE*Pa21**2*SW*Zh11*la6*sb-6D0*MW/EE*Pa21
     &      **2*SW*Zh11*la6*sb**3-2D0*MW/EE*Pa21**2*SW*Zh11
     &      *la7*sb+2D0*MW/EE*Pa21**2*SW*Zh11*la7*sb**3-2D0
     &      *MW/EE*Pa21**2*SW*Zh12*cb*la7+2D0*MW/EE*Pa21
     &      **2*SW*Zh12*la2*sb-2D0*MW/EE*Pa21**2*SW*Zh12
     &      *la2*sb**3+2D0*MW/EE*Pa21**2*SW*Zh12*la3*sb**3
     &      +2D0*MW/EE*Pa21**2*SW*Zh12*la4*sb**3-4D0*MW/EE
     &      *Pa21**2*SW*Zh12*la5*sb+2D0*MW/EE*Pa21**2*SW
     &      *Zh12*la5*sb**3+2D0*MW/EE*Pa22**2*SW*Zh11*cb
     &      *la1s-2D0*MW/EE*Pa22**2*SW*Zh11*la4s*sb+2D0*MW
     &      /EE*Pa22**2*SW*Zh11*la5s*sb+2D0*MW/EE*Pa22**2
     &      *SW*Zh11*la6s*sb-2D0*MW/EE*Pa22**2*SW*Zh12*cb
     &      *la4s+2D0*MW/EE*Pa22**2*SW*Zh12*cb*la5s+2D0*MW
     &      /EE*Pa22**2*SW*Zh12*cb*la6s+2D0*MW/EE*Pa22**2
     &      *SW*Zh12*la2s*sb+2D0*Pa21*Pa22*Sqrt2*Zh11*cb
     &      *la5s*xvev-2D0*Pa21*Pa22*Sqrt2*Zh11*cb*la6s*xvev
     &      +2D0*Pa21*Pa22*Sqrt2*Zh12*la5s*sb*xvev-2D0*Pa21
     &      *Pa22*Sqrt2*Zh12*la6s*sb*xvev+2D0*Pa21**2*Sqrt2
     &      *Zh13*cb*la4s*sb*xvev+2D0*Pa21**2*Sqrt2*Zh13
     &      *cb*la5s*sb*xvev+2D0*Pa21**2*Sqrt2*Zh13*cb*la6s
     &      *sb*xvev+2D0*MW/EE*Pa21**2*SW*Zh11*cb*la1*sb
     &      **2-2D0*MW/EE*Pa21**2*SW*Zh11*cb*la3*sb**2-2D0
     &      *MW/EE*Pa21**2*SW*Zh11*cb*la4*sb**2-2D0*MW/EE
     &      *Pa21**2*SW*Zh11*cb*la5*sb**2-2D0*MW/EE*Pa21
     &      **2*SW*Zh12*cb*la6*sb**2+6D0*MW/EE*Pa21**2*SW
     &      *Zh12*cb*la7*sb**2
      AAABR(159) = 2D0*Pa11*Pa21*Sqrt2*Zh23*aa2+2D0*Pa12
     &      *Pa22*Sqrt2*Zh23*aa3-6D0*Pa12*Pa22*Sqrt2*Zh23
     &      *aa4+2D0*Pa11*Pa21*Sqrt2*Zh23*aa1*sb**2-2D0*Pa11
     &      *Pa21*Sqrt2*Zh23*aa2*sb**2+2D0*Pa11*Pa21*Sqrt2
     &      *Zh23*la2s*xvev+Pa11*Pa22*Sqrt2*Zh21*aa5*cb-Pa11
     &      *Pa22*Sqrt2*Zh21*aa6*cb+Pa11*Pa22*Sqrt2*Zh22
     &      *aa5*sb-Pa11*Pa22*Sqrt2*Zh22*aa6*sb+Pa12*Pa21
     &      *Sqrt2*Zh21*aa5*cb-Pa12*Pa21*Sqrt2*Zh21*aa6*cb
     &      +Pa12*Pa21*Sqrt2*Zh22*aa5*sb-Pa12*Pa21*Sqrt2
     &      *Zh22*aa6*sb+4D0*Pa12*Pa22*Sqrt2*Zh23*la3s*xvev
     &      -24D0*Pa12*Pa22*Sqrt2*Zh23*la7s*xvev+4D0*MW/EE
     &      *Pa11*Pa22*SW*Zh23*la5s-4D0*MW/EE*Pa11*Pa22*SW
     &      *Zh23*la6s+4D0*MW/EE*Pa12*Pa21*SW*Zh23*la5s-4D0
     &      *MW/EE*Pa12*Pa21*SW*Zh23*la6s+2D0*Pa11*Pa21*Sqrt2
     &      *Zh23*aa5*cb*sb+2D0*Pa11*Pa21*Sqrt2*Zh23*aa6
     &      *cb*sb+2D0*Pa11*Pa21*Sqrt2*Zh23*la1s*sb**2*xvev
     &      -2D0*Pa11*Pa21*Sqrt2*Zh23*la2s*sb**2*xvev+2D0
     &      *Pa11*Pa22*Sqrt2*Zh21*cb*la5s*xvev-2D0*Pa11*Pa22
     &      *Sqrt2*Zh21*cb*la6s*xvev+2D0*Pa11*Pa22*Sqrt2
     &      *Zh22*la5s*sb*xvev-2D0*Pa11*Pa22*Sqrt2*Zh22*la6s
     &      *sb*xvev+2D0*Pa12*Pa21*Sqrt2*Zh21*cb*la5s*xvev
     &      -2D0*Pa12*Pa21*Sqrt2*Zh21*cb*la6s*xvev+2D0*Pa12
     &      *Pa21*Sqrt2*Zh22*la5s*sb*xvev-2D0*Pa12*Pa21*Sqrt2
     &      *Zh22*la6s*sb*xvev+4D0*MW/EE*Pa11*Pa21*SW*Zh21
     &      *cb*la3+4D0*MW/EE*Pa11*Pa21*SW*Zh21*cb*la4-4D0
     &      *MW/EE*Pa11*Pa21*SW*Zh21*cb*la5+8D0*MW/EE*Pa11
     &      *Pa21*SW*Zh21*la6*sb-12D0*MW/EE*Pa11*Pa21*SW
     &      *Zh21*la6*sb**3-4D0*MW/EE*Pa11*Pa21*SW*Zh21*la7
     &      *sb+4D0*MW/EE*Pa11*Pa21*SW*Zh21*la7*sb**3-4D0
     &      *MW/EE*Pa11*Pa21*SW*Zh22*cb*la7+4D0*MW/EE*Pa11
     &      *Pa21*SW*Zh22*la2*sb-4D0*MW/EE*Pa11*Pa21*SW*Zh22
     &      *la2*sb**3+4D0*MW/EE*Pa11*Pa21*SW*Zh22*la3*sb
     &      **3+4D0*MW/EE*Pa11*Pa21*SW*Zh22*la4*sb**3-8D0
     &      *MW/EE*Pa11*Pa21*SW*Zh22*la5*sb+4D0*MW/EE*Pa11
     &      *Pa21*SW*Zh22*la5*sb**3+4D0*MW/EE*Pa12*Pa22*SW
     &      *Zh21*cb*la1s-4D0*MW/EE*Pa12*Pa22*SW*Zh21*la4s
     &      *sb+4D0*MW/EE*Pa12*Pa22*SW*Zh21*la5s*sb+4D0*MW
     &      /EE*Pa12*Pa22*SW*Zh21*la6s*sb-4D0*MW/EE*Pa12
     &      *Pa22*SW*Zh22*cb*la4s+4D0*MW/EE*Pa12*Pa22*SW
     &      *Zh22*cb*la5s+4D0*MW/EE*Pa12*Pa22*SW*Zh22*cb
     &      *la6s+4D0*MW/EE*Pa12*Pa22*SW*Zh22*la2s*sb+4D0
     &      *Pa11*Pa21*Sqrt2*Zh23*cb*la4s*sb*xvev+4D0*Pa11
     &      *Pa21*Sqrt2*Zh23*cb*la5s*sb*xvev+4D0*Pa11*Pa21
     &      *Sqrt2*Zh23*cb*la6s*sb*xvev+4D0*MW/EE*Pa11*Pa21
     &      *SW*Zh21*cb*la1*sb**2-4D0*MW/EE*Pa11*Pa21*SW
     &      *Zh21*cb*la3*sb**2-4D0*MW/EE*Pa11*Pa21*SW*Zh21
     &      *cb*la4*sb**2-4D0*MW/EE*Pa11*Pa21*SW*Zh21*cb
     &      *la5*sb**2-4D0*MW/EE*Pa11*Pa21*SW*Zh22*cb*la6
     &      *sb**2+12D0*MW/EE*Pa11*Pa21*SW*Zh22*cb*la7*sb**2
      AAABR(160) = Pa21**2*Sqrt2*Zh23*aa2+Pa22**2*Sqrt2*Zh23
     &      *aa3-3D0*Pa22**2*Sqrt2*Zh23*aa4+Pa21**2*Sqrt2
     &      *Zh23*aa1*sb**2-Pa21**2*Sqrt2*Zh23*aa2*sb**2
     &      +Pa21**2*Sqrt2*Zh23*la2s*xvev+2D0*Pa22**2*Sqrt2
     &      *Zh23*la3s*xvev-12D0*Pa22**2*Sqrt2*Zh23*la7s
     &      *xvev+Pa21*Pa22*Sqrt2*Zh21*aa5*cb-Pa21*Pa22*Sqrt2
     &      *Zh21*aa6*cb+Pa21*Pa22*Sqrt2*Zh22*aa5*sb-Pa21
     &      *Pa22*Sqrt2*Zh22*aa6*sb+Pa21**2*Sqrt2*Zh23*aa5
     &      *cb*sb+Pa21**2*Sqrt2*Zh23*aa6*cb*sb+Pa21**2*Sqrt2
     &      *Zh23*la1s*sb**2*xvev-Pa21**2*Sqrt2*Zh23*la2s
     &      *sb**2*xvev+4D0*MW/EE*Pa21*Pa22*SW*Zh23*la5s
     &      -4D0*MW/EE*Pa21*Pa22*SW*Zh23*la6s+2D0*MW/EE*Pa21
     &      **2*SW*Zh21*cb*la3+2D0*MW/EE*Pa21**2*SW*Zh21
     &      *cb*la4-2D0*MW/EE*Pa21**2*SW*Zh21*cb*la5+4D0
     &      *MW/EE*Pa21**2*SW*Zh21*la6*sb-6D0*MW/EE*Pa21
     &      **2*SW*Zh21*la6*sb**3-2D0*MW/EE*Pa21**2*SW*Zh21
     &      *la7*sb+2D0*MW/EE*Pa21**2*SW*Zh21*la7*sb**3-2D0
     &      *MW/EE*Pa21**2*SW*Zh22*cb*la7+2D0*MW/EE*Pa21
     &      **2*SW*Zh22*la2*sb-2D0*MW/EE*Pa21**2*SW*Zh22
     &      *la2*sb**3+2D0*MW/EE*Pa21**2*SW*Zh22*la3*sb**3
     &      +2D0*MW/EE*Pa21**2*SW*Zh22*la4*sb**3-4D0*MW/EE
     &      *Pa21**2*SW*Zh22*la5*sb+2D0*MW/EE*Pa21**2*SW
     &      *Zh22*la5*sb**3+2D0*MW/EE*Pa22**2*SW*Zh21*cb
     &      *la1s-2D0*MW/EE*Pa22**2*SW*Zh21*la4s*sb+2D0*MW
     &      /EE*Pa22**2*SW*Zh21*la5s*sb+2D0*MW/EE*Pa22**2
     &      *SW*Zh21*la6s*sb-2D0*MW/EE*Pa22**2*SW*Zh22*cb
     &      *la4s+2D0*MW/EE*Pa22**2*SW*Zh22*cb*la5s+2D0*MW
     &      /EE*Pa22**2*SW*Zh22*cb*la6s+2D0*MW/EE*Pa22**2
     &      *SW*Zh22*la2s*sb+2D0*Pa21*Pa22*Sqrt2*Zh21*cb
     &      *la5s*xvev-2D0*Pa21*Pa22*Sqrt2*Zh21*cb*la6s*xvev
     &      +2D0*Pa21*Pa22*Sqrt2*Zh22*la5s*sb*xvev-2D0*Pa21
     &      *Pa22*Sqrt2*Zh22*la6s*sb*xvev+2D0*Pa21**2*Sqrt2
     &      *Zh23*cb*la4s*sb*xvev+2D0*Pa21**2*Sqrt2*Zh23
     &      *cb*la5s*sb*xvev+2D0*Pa21**2*Sqrt2*Zh23*cb*la6s
     &      *sb*xvev+2D0*MW/EE*Pa21**2*SW*Zh21*cb*la1*sb
     &      **2-2D0*MW/EE*Pa21**2*SW*Zh21*cb*la3*sb**2-2D0
     &      *MW/EE*Pa21**2*SW*Zh21*cb*la4*sb**2-2D0*MW/EE
     &      *Pa21**2*SW*Zh21*cb*la5*sb**2-2D0*MW/EE*Pa21
     &      **2*SW*Zh22*cb*la6*sb**2+6D0*MW/EE*Pa21**2*SW
     &      *Zh22*cb*la7*sb**2
      AAABR(161) = 2D0*Pa11*Pa21*Sqrt2*Zh33*aa2+2D0*Pa12
     &      *Pa22*Sqrt2*Zh33*aa3-6D0*Pa12*Pa22*Sqrt2*Zh33
     &      *aa4+2D0*Pa11*Pa21*Sqrt2*Zh33*aa1*sb**2-2D0*Pa11
     &      *Pa21*Sqrt2*Zh33*aa2*sb**2+2D0*Pa11*Pa21*Sqrt2
     &      *Zh33*la2s*xvev+Pa11*Pa22*Sqrt2*Zh31*aa5*cb-Pa11
     &      *Pa22*Sqrt2*Zh31*aa6*cb+Pa11*Pa22*Sqrt2*Zh32
     &      *aa5*sb-Pa11*Pa22*Sqrt2*Zh32*aa6*sb+Pa12*Pa21
     &      *Sqrt2*Zh31*aa5*cb-Pa12*Pa21*Sqrt2*Zh31*aa6*cb
     &      +Pa12*Pa21*Sqrt2*Zh32*aa5*sb-Pa12*Pa21*Sqrt2
     &      *Zh32*aa6*sb+4D0*Pa12*Pa22*Sqrt2*Zh33*la3s*xvev
     &      -24D0*Pa12*Pa22*Sqrt2*Zh33*la7s*xvev+4D0*MW/EE
     &      *Pa11*Pa22*SW*Zh33*la5s-4D0*MW/EE*Pa11*Pa22*SW
     &      *Zh33*la6s+4D0*MW/EE*Pa12*Pa21*SW*Zh33*la5s-4D0
     &      *MW/EE*Pa12*Pa21*SW*Zh33*la6s+2D0*Pa11*Pa21*Sqrt2
     &      *Zh33*aa5*cb*sb+2D0*Pa11*Pa21*Sqrt2*Zh33*aa6
     &      *cb*sb+2D0*Pa11*Pa21*Sqrt2*Zh33*la1s*sb**2*xvev
     &      -2D0*Pa11*Pa21*Sqrt2*Zh33*la2s*sb**2*xvev+2D0
     &      *Pa11*Pa22*Sqrt2*Zh31*cb*la5s*xvev-2D0*Pa11*Pa22
     &      *Sqrt2*Zh31*cb*la6s*xvev+2D0*Pa11*Pa22*Sqrt2
     &      *Zh32*la5s*sb*xvev-2D0*Pa11*Pa22*Sqrt2*Zh32*la6s
     &      *sb*xvev+2D0*Pa12*Pa21*Sqrt2*Zh31*cb*la5s*xvev
     &      -2D0*Pa12*Pa21*Sqrt2*Zh31*cb*la6s*xvev+2D0*Pa12
     &      *Pa21*Sqrt2*Zh32*la5s*sb*xvev-2D0*Pa12*Pa21*Sqrt2
     &      *Zh32*la6s*sb*xvev+4D0*MW/EE*Pa11*Pa21*SW*Zh31
     &      *cb*la3+4D0*MW/EE*Pa11*Pa21*SW*Zh31*cb*la4-4D0
     &      *MW/EE*Pa11*Pa21*SW*Zh31*cb*la5+8D0*MW/EE*Pa11
     &      *Pa21*SW*Zh31*la6*sb-12D0*MW/EE*Pa11*Pa21*SW
     &      *Zh31*la6*sb**3-4D0*MW/EE*Pa11*Pa21*SW*Zh31*la7
     &      *sb+4D0*MW/EE*Pa11*Pa21*SW*Zh31*la7*sb**3-4D0
     &      *MW/EE*Pa11*Pa21*SW*Zh32*cb*la7+4D0*MW/EE*Pa11
     &      *Pa21*SW*Zh32*la2*sb-4D0*MW/EE*Pa11*Pa21*SW*Zh32
     &      *la2*sb**3+4D0*MW/EE*Pa11*Pa21*SW*Zh32*la3*sb
     &      **3+4D0*MW/EE*Pa11*Pa21*SW*Zh32*la4*sb**3-8D0
     &      *MW/EE*Pa11*Pa21*SW*Zh32*la5*sb+4D0*MW/EE*Pa11
     &      *Pa21*SW*Zh32*la5*sb**3+4D0*MW/EE*Pa12*Pa22*SW
     &      *Zh31*cb*la1s-4D0*MW/EE*Pa12*Pa22*SW*Zh31*la4s
     &      *sb+4D0*MW/EE*Pa12*Pa22*SW*Zh31*la5s*sb+4D0*MW
     &      /EE*Pa12*Pa22*SW*Zh31*la6s*sb-4D0*MW/EE*Pa12
     &      *Pa22*SW*Zh32*cb*la4s+4D0*MW/EE*Pa12*Pa22*SW
     &      *Zh32*cb*la5s+4D0*MW/EE*Pa12*Pa22*SW*Zh32*cb
     &      *la6s+4D0*MW/EE*Pa12*Pa22*SW*Zh32*la2s*sb+4D0
     &      *Pa11*Pa21*Sqrt2*Zh33*cb*la4s*sb*xvev+4D0*Pa11
     &      *Pa21*Sqrt2*Zh33*cb*la5s*sb*xvev+4D0*Pa11*Pa21
     &      *Sqrt2*Zh33*cb*la6s*sb*xvev+4D0*MW/EE*Pa11*Pa21
     &      *SW*Zh31*cb*la1*sb**2-4D0*MW/EE*Pa11*Pa21*SW
     &      *Zh31*cb*la3*sb**2-4D0*MW/EE*Pa11*Pa21*SW*Zh31
     &      *cb*la4*sb**2-4D0*MW/EE*Pa11*Pa21*SW*Zh31*cb
     &      *la5*sb**2-4D0*MW/EE*Pa11*Pa21*SW*Zh32*cb*la6
     &      *sb**2+12D0*MW/EE*Pa11*Pa21*SW*Zh32*cb*la7*sb**2
      AAABR(162) = Pa21**2*Sqrt2*Zh33*aa2+Pa22**2*Sqrt2*Zh33
     &      *aa3-3D0*Pa22**2*Sqrt2*Zh33*aa4+Pa21**2*Sqrt2
     &      *Zh33*aa1*sb**2-Pa21**2*Sqrt2*Zh33*aa2*sb**2
     &      +Pa21**2*Sqrt2*Zh33*la2s*xvev+2D0*Pa22**2*Sqrt2
     &      *Zh33*la3s*xvev-12D0*Pa22**2*Sqrt2*Zh33*la7s
     &      *xvev+Pa21*Pa22*Sqrt2*Zh31*aa5*cb-Pa21*Pa22*Sqrt2
     &      *Zh31*aa6*cb+Pa21*Pa22*Sqrt2*Zh32*aa5*sb-Pa21
     &      *Pa22*Sqrt2*Zh32*aa6*sb+Pa21**2*Sqrt2*Zh33*aa5
     &      *cb*sb+Pa21**2*Sqrt2*Zh33*aa6*cb*sb+Pa21**2*Sqrt2
     &      *Zh33*la1s*sb**2*xvev-Pa21**2*Sqrt2*Zh33*la2s
     &      *sb**2*xvev+4D0*MW/EE*Pa21*Pa22*SW*Zh33*la5s
     &      -4D0*MW/EE*Pa21*Pa22*SW*Zh33*la6s+2D0*MW/EE*Pa21
     &      **2*SW*Zh31*cb*la3+2D0*MW/EE*Pa21**2*SW*Zh31
     &      *cb*la4-2D0*MW/EE*Pa21**2*SW*Zh31*cb*la5+4D0
     &      *MW/EE*Pa21**2*SW*Zh31*la6*sb-6D0*MW/EE*Pa21
     &      **2*SW*Zh31*la6*sb**3-2D0*MW/EE*Pa21**2*SW*Zh31
     &      *la7*sb+2D0*MW/EE*Pa21**2*SW*Zh31*la7*sb**3-2D0
     &      *MW/EE*Pa21**2*SW*Zh32*cb*la7+2D0*MW/EE*Pa21
     &      **2*SW*Zh32*la2*sb-2D0*MW/EE*Pa21**2*SW*Zh32
     &      *la2*sb**3+2D0*MW/EE*Pa21**2*SW*Zh32*la3*sb**3
     &      +2D0*MW/EE*Pa21**2*SW*Zh32*la4*sb**3-4D0*MW/EE
     &      *Pa21**2*SW*Zh32*la5*sb+2D0*MW/EE*Pa21**2*SW
     &      *Zh32*la5*sb**3+2D0*MW/EE*Pa22**2*SW*Zh31*cb
     &      *la1s-2D0*MW/EE*Pa22**2*SW*Zh31*la4s*sb+2D0*MW
     &      /EE*Pa22**2*SW*Zh31*la5s*sb+2D0*MW/EE*Pa22**2
     &      *SW*Zh31*la6s*sb-2D0*MW/EE*Pa22**2*SW*Zh32*cb
     &      *la4s+2D0*MW/EE*Pa22**2*SW*Zh32*cb*la5s+2D0*MW
     &      /EE*Pa22**2*SW*Zh32*cb*la6s+2D0*MW/EE*Pa22**2
     &      *SW*Zh32*la2s*sb+2D0*Pa21*Pa22*Sqrt2*Zh31*cb
     &      *la5s*xvev-2D0*Pa21*Pa22*Sqrt2*Zh31*cb*la6s*xvev
     &      +2D0*Pa21*Pa22*Sqrt2*Zh32*la5s*sb*xvev-2D0*Pa21
     &      *Pa22*Sqrt2*Zh32*la6s*sb*xvev+2D0*Pa21**2*Sqrt2
     &      *Zh33*cb*la4s*sb*xvev+2D0*Pa21**2*Sqrt2*Zh33
     &      *cb*la5s*sb*xvev+2D0*Pa21**2*Sqrt2*Zh33*cb*la6s
     &      *sb*xvev+2D0*MW/EE*Pa21**2*SW*Zh31*cb*la1*sb
     &      **2-2D0*MW/EE*Pa21**2*SW*Zh31*cb*la3*sb**2-2D0
     &      *MW/EE*Pa21**2*SW*Zh31*cb*la4*sb**2-2D0*MW/EE
     &      *Pa21**2*SW*Zh31*cb*la5*sb**2-2D0*MW/EE*Pa21
     &      **2*SW*Zh32*cb*la6*sb**2+6D0*MW/EE*Pa21**2*SW
     &      *Zh32*cb*la7*sb**2
      AAABR(163) = 2D0*Pa11*Pa21*Sqrt2*Zh13*aa2+2D0*Pa12
     &      *Pa22*Sqrt2*Zh13*aa3-6D0*Pa12*Pa22*Sqrt2*Zh13
     &      *aa4+2D0*Pa11*Pa21*Sqrt2*Zh13*aa1*sb**2-2D0*Pa11
     &      *Pa21*Sqrt2*Zh13*aa2*sb**2+2D0*Pa11*Pa21*Sqrt2
     &      *Zh13*la2s*xvev+Pa11*Pa22*Sqrt2*Zh11*aa5*cb-Pa11
     &      *Pa22*Sqrt2*Zh11*aa6*cb+Pa11*Pa22*Sqrt2*Zh12
     &      *aa5*sb-Pa11*Pa22*Sqrt2*Zh12*aa6*sb+Pa12*Pa21
     &      *Sqrt2*Zh11*aa5*cb-Pa12*Pa21*Sqrt2*Zh11*aa6*cb
     &      +Pa12*Pa21*Sqrt2*Zh12*aa5*sb-Pa12*Pa21*Sqrt2
     &      *Zh12*aa6*sb+4D0*Pa12*Pa22*Sqrt2*Zh13*la3s*xvev
     &      -24D0*Pa12*Pa22*Sqrt2*Zh13*la7s*xvev+4D0*MW/EE
     &      *Pa11*Pa22*SW*Zh13*la5s-4D0*MW/EE*Pa11*Pa22*SW
     &      *Zh13*la6s+4D0*MW/EE*Pa12*Pa21*SW*Zh13*la5s-4D0
     &      *MW/EE*Pa12*Pa21*SW*Zh13*la6s+2D0*Pa11*Pa21*Sqrt2
     &      *Zh13*aa5*cb*sb+2D0*Pa11*Pa21*Sqrt2*Zh13*aa6
     &      *cb*sb+2D0*Pa11*Pa21*Sqrt2*Zh13*la1s*sb**2*xvev
     &      -2D0*Pa11*Pa21*Sqrt2*Zh13*la2s*sb**2*xvev+2D0
     &      *Pa11*Pa22*Sqrt2*Zh11*cb*la5s*xvev-2D0*Pa11*Pa22
     &      *Sqrt2*Zh11*cb*la6s*xvev+2D0*Pa11*Pa22*Sqrt2
     &      *Zh12*la5s*sb*xvev-2D0*Pa11*Pa22*Sqrt2*Zh12*la6s
     &      *sb*xvev+2D0*Pa12*Pa21*Sqrt2*Zh11*cb*la5s*xvev
     &      -2D0*Pa12*Pa21*Sqrt2*Zh11*cb*la6s*xvev+2D0*Pa12
     &      *Pa21*Sqrt2*Zh12*la5s*sb*xvev-2D0*Pa12*Pa21*Sqrt2
     &      *Zh12*la6s*sb*xvev+4D0*MW/EE*Pa11*Pa21*SW*Zh11
     &      *cb*la3+4D0*MW/EE*Pa11*Pa21*SW*Zh11*cb*la4-4D0
     &      *MW/EE*Pa11*Pa21*SW*Zh11*cb*la5+8D0*MW/EE*Pa11
     &      *Pa21*SW*Zh11*la6*sb-12D0*MW/EE*Pa11*Pa21*SW
     &      *Zh11*la6*sb**3-4D0*MW/EE*Pa11*Pa21*SW*Zh11*la7
     &      *sb+4D0*MW/EE*Pa11*Pa21*SW*Zh11*la7*sb**3-4D0
     &      *MW/EE*Pa11*Pa21*SW*Zh12*cb*la7+4D0*MW/EE*Pa11
     &      *Pa21*SW*Zh12*la2*sb-4D0*MW/EE*Pa11*Pa21*SW*Zh12
     &      *la2*sb**3+4D0*MW/EE*Pa11*Pa21*SW*Zh12*la3*sb
     &      **3+4D0*MW/EE*Pa11*Pa21*SW*Zh12*la4*sb**3-8D0
     &      *MW/EE*Pa11*Pa21*SW*Zh12*la5*sb+4D0*MW/EE*Pa11
     &      *Pa21*SW*Zh12*la5*sb**3+4D0*MW/EE*Pa12*Pa22*SW
     &      *Zh11*cb*la1s-4D0*MW/EE*Pa12*Pa22*SW*Zh11*la4s
     &      *sb+4D0*MW/EE*Pa12*Pa22*SW*Zh11*la5s*sb+4D0*MW
     &      /EE*Pa12*Pa22*SW*Zh11*la6s*sb-4D0*MW/EE*Pa12
     &      *Pa22*SW*Zh12*cb*la4s+4D0*MW/EE*Pa12*Pa22*SW
     &      *Zh12*cb*la5s+4D0*MW/EE*Pa12*Pa22*SW*Zh12*cb
     &      *la6s+4D0*MW/EE*Pa12*Pa22*SW*Zh12*la2s*sb+4D0
     &      *Pa11*Pa21*Sqrt2*Zh13*cb*la4s*sb*xvev+4D0*Pa11
     &      *Pa21*Sqrt2*Zh13*cb*la5s*sb*xvev+4D0*Pa11*Pa21
     &      *Sqrt2*Zh13*cb*la6s*sb*xvev+4D0*MW/EE*Pa11*Pa21
     &      *SW*Zh11*cb*la1*sb**2-4D0*MW/EE*Pa11*Pa21*SW
     &      *Zh11*cb*la3*sb**2-4D0*MW/EE*Pa11*Pa21*SW*Zh11
     &      *cb*la4*sb**2-4D0*MW/EE*Pa11*Pa21*SW*Zh11*cb
     &      *la5*sb**2-4D0*MW/EE*Pa11*Pa21*SW*Zh12*cb*la6
     &      *sb**2+12D0*MW/EE*Pa11*Pa21*SW*Zh12*cb*la7*sb**2
      AAABR(164) = Pa11**2*Sqrt2*Zh13*aa2+Pa12**2*Sqrt2*Zh13
     &      *aa3-3D0*Pa12**2*Sqrt2*Zh13*aa4+Pa11**2*Sqrt2
     &      *Zh13*aa1*sb**2-Pa11**2*Sqrt2*Zh13*aa2*sb**2
     &      +Pa11**2*Sqrt2*Zh13*la2s*xvev+2D0*Pa12**2*Sqrt2
     &      *Zh13*la3s*xvev-12D0*Pa12**2*Sqrt2*Zh13*la7s
     &      *xvev+Pa11*Pa12*Sqrt2*Zh11*aa5*cb-Pa11*Pa12*Sqrt2
     &      *Zh11*aa6*cb+Pa11*Pa12*Sqrt2*Zh12*aa5*sb-Pa11
     &      *Pa12*Sqrt2*Zh12*aa6*sb+Pa11**2*Sqrt2*Zh13*aa5
     &      *cb*sb+Pa11**2*Sqrt2*Zh13*aa6*cb*sb+Pa11**2*Sqrt2
     &      *Zh13*la1s*sb**2*xvev-Pa11**2*Sqrt2*Zh13*la2s
     &      *sb**2*xvev+4D0*MW/EE*Pa11*Pa12*SW*Zh13*la5s
     &      -4D0*MW/EE*Pa11*Pa12*SW*Zh13*la6s+2D0*MW/EE*Pa11
     &      **2*SW*Zh11*cb*la3+2D0*MW/EE*Pa11**2*SW*Zh11
     &      *cb*la4-2D0*MW/EE*Pa11**2*SW*Zh11*cb*la5+4D0
     &      *MW/EE*Pa11**2*SW*Zh11*la6*sb-6D0*MW/EE*Pa11
     &      **2*SW*Zh11*la6*sb**3-2D0*MW/EE*Pa11**2*SW*Zh11
     &      *la7*sb+2D0*MW/EE*Pa11**2*SW*Zh11*la7*sb**3-2D0
     &      *MW/EE*Pa11**2*SW*Zh12*cb*la7+2D0*MW/EE*Pa11
     &      **2*SW*Zh12*la2*sb-2D0*MW/EE*Pa11**2*SW*Zh12
     &      *la2*sb**3+2D0*MW/EE*Pa11**2*SW*Zh12*la3*sb**3
     &      +2D0*MW/EE*Pa11**2*SW*Zh12*la4*sb**3-4D0*MW/EE
     &      *Pa11**2*SW*Zh12*la5*sb+2D0*MW/EE*Pa11**2*SW
     &      *Zh12*la5*sb**3+2D0*MW/EE*Pa12**2*SW*Zh11*cb
     &      *la1s-2D0*MW/EE*Pa12**2*SW*Zh11*la4s*sb+2D0*MW
     &      /EE*Pa12**2*SW*Zh11*la5s*sb+2D0*MW/EE*Pa12**2
     &      *SW*Zh11*la6s*sb-2D0*MW/EE*Pa12**2*SW*Zh12*cb
     &      *la4s+2D0*MW/EE*Pa12**2*SW*Zh12*cb*la5s+2D0*MW
     &      /EE*Pa12**2*SW*Zh12*cb*la6s+2D0*MW/EE*Pa12**2
     &      *SW*Zh12*la2s*sb+2D0*Pa11*Pa12*Sqrt2*Zh11*cb
     &      *la5s*xvev-2D0*Pa11*Pa12*Sqrt2*Zh11*cb*la6s*xvev
     &      +2D0*Pa11*Pa12*Sqrt2*Zh12*la5s*sb*xvev-2D0*Pa11
     &      *Pa12*Sqrt2*Zh12*la6s*sb*xvev+2D0*Pa11**2*Sqrt2
     &      *Zh13*cb*la4s*sb*xvev+2D0*Pa11**2*Sqrt2*Zh13
     &      *cb*la5s*sb*xvev+2D0*Pa11**2*Sqrt2*Zh13*cb*la6s
     &      *sb*xvev+2D0*MW/EE*Pa11**2*SW*Zh11*cb*la1*sb
     &      **2-2D0*MW/EE*Pa11**2*SW*Zh11*cb*la3*sb**2-2D0
     &      *MW/EE*Pa11**2*SW*Zh11*cb*la4*sb**2-2D0*MW/EE
     &      *Pa11**2*SW*Zh11*cb*la5*sb**2-2D0*MW/EE*Pa11
     &      **2*SW*Zh12*cb*la6*sb**2+6D0*MW/EE*Pa11**2*SW
     &      *Zh12*cb*la7*sb**2
      AAABR(165) = Sqrt2*Zh11**2*Zh23*aa1+Sqrt2*Zh12**2*Zh23
     &      *aa2+3D0*Sqrt2*Zh13**2*Zh23*aa3+3D0*Sqrt2*Zh13
     &      **2*Zh23*aa4-Sqrt2*Zh11*Zh12*Zh23*aa5-Sqrt2*Zh11
     &      *Zh12*Zh23*aa6+2D0*Sqrt2*Zh11*Zh13*Zh21*aa1-Sqrt2
     &      *Zh11*Zh13*Zh22*aa5-Sqrt2*Zh11*Zh13*Zh22*aa6
     &      +Sqrt2*Zh11**2*Zh23*la1s*xvev-Sqrt2*Zh12*Zh13
     &      *Zh21*aa5-Sqrt2*Zh12*Zh13*Zh21*aa6+2D0*Sqrt2
     &      *Zh12*Zh13*Zh22*aa2+Sqrt2*Zh12**2*Zh23*la2s*xvev
     &      +6D0*Sqrt2*Zh13**2*Zh23*la3s*xvev+12D0*Sqrt2
     &      *Zh13**2*Zh23*la7s*xvev+12D0*Sqrt2*Zh13**2*Zh23
     &      *la8s*xvev-2D0*Sqrt2*Zh11*Zh12*Zh23*la4s*xvev
     &      -2D0*Sqrt2*Zh11*Zh12*Zh23*la5s*xvev-2D0*Sqrt2
     &      *Zh11*Zh12*Zh23*la6s*xvev+2D0*Sqrt2*Zh11*Zh13
     &      *Zh21*la1s*xvev-2D0*Sqrt2*Zh11*Zh13*Zh22*la4s
     &      *xvev-2D0*Sqrt2*Zh11*Zh13*Zh22*la5s*xvev-2D0
     &      *Sqrt2*Zh11*Zh13*Zh22*la6s*xvev-2D0*Sqrt2*Zh12
     &      *Zh13*Zh21*la4s*xvev-2D0*Sqrt2*Zh12*Zh13*Zh21
     &      *la5s*xvev-2D0*Sqrt2*Zh12*Zh13*Zh21*la6s*xvev
     &      +2D0*Sqrt2*Zh12*Zh13*Zh22*la2s*xvev+6D0*MW/EE
     &      *SW*Zh11**2*Zh21*cb*la1-6D0*MW/EE*SW*Zh11**2
     &      *Zh21*la6*sb-6D0*MW/EE*SW*Zh11**2*Zh22*cb*la6
     &      +2D0*MW/EE*SW*Zh11**2*Zh22*la3*sb+2D0*MW/EE*SW
     &      *Zh11**2*Zh22*la4*sb+2D0*MW/EE*SW*Zh11**2*Zh22
     &      *la5*sb+2D0*MW/EE*SW*Zh12**2*Zh21*cb*la3+2D0
     &      *MW/EE*SW*Zh12**2*Zh21*cb*la4+2D0*MW/EE*SW*Zh12
     &      **2*Zh21*cb*la5-6D0*MW/EE*SW*Zh12**2*Zh21*la7
     &      *sb-6D0*MW/EE*SW*Zh12**2*Zh22*cb*la7+6D0*MW/EE
     &      *SW*Zh12**2*Zh22*la2*sb+2D0*MW/EE*SW*Zh13**2
     &      *Zh21*cb*la1s-2D0*MW/EE*SW*Zh13**2*Zh21*la4s
     &      *sb-2D0*MW/EE*SW*Zh13**2*Zh21*la5s*sb-2D0*MW
     &      /EE*SW*Zh13**2*Zh21*la6s*sb-2D0*MW/EE*SW*Zh13
     &      **2*Zh22*cb*la4s-2D0*MW/EE*SW*Zh13**2*Zh22*cb
     &      *la5s-2D0*MW/EE*SW*Zh13**2*Zh22*cb*la6s+2D0*MW
     &      /EE*SW*Zh13**2*Zh22*la2s*sb-12D0*MW/EE*SW*Zh11
     &      *Zh12*Zh21*cb*la6+4D0*MW/EE*SW*Zh11*Zh12*Zh21
     &      *la3*sb+4D0*MW/EE*SW*Zh11*Zh12*Zh21*la4*sb+4D0
     &      *MW/EE*SW*Zh11*Zh12*Zh21*la5*sb+4D0*MW/EE*SW
     &      *Zh11*Zh12*Zh22*cb*la3+4D0*MW/EE*SW*Zh11*Zh12
     &      *Zh22*cb*la4+4D0*MW/EE*SW*Zh11*Zh12*Zh22*cb*la5
     &      -12D0*MW/EE*SW*Zh11*Zh12*Zh22*la7*sb+4D0*MW/EE
     &      *SW*Zh11*Zh13*Zh23*cb*la1s-4D0*MW/EE*SW*Zh11
     &      *Zh13*Zh23*la4s*sb-4D0*MW/EE*SW*Zh11*Zh13*Zh23
     &      *la5s*sb-4D0*MW/EE*SW*Zh11*Zh13*Zh23*la6s*sb
     &      -4D0*MW/EE*SW*Zh12*Zh13*Zh23*cb*la4s-4D0*MW/EE
     &      *SW*Zh12*Zh13*Zh23*cb*la5s-4D0*MW/EE*SW*Zh12
     &      *Zh13*Zh23*cb*la6s+4D0*MW/EE*SW*Zh12*Zh13*Zh23
     &      *la2s*sb
      AAABR(166) = Sqrt2*Zh13*Zh21**2*aa1+Sqrt2*Zh13*Zh22
     &      **2*aa2+3D0*Sqrt2*Zh13*Zh23**2*aa3+3D0*Sqrt2
     &      *Zh13*Zh23**2*aa4+2D0*Sqrt2*Zh11*Zh21*Zh23*aa1
     &      -Sqrt2*Zh11*Zh22*Zh23*aa5-Sqrt2*Zh11*Zh22*Zh23
     &      *aa6-Sqrt2*Zh12*Zh21*Zh23*aa5-Sqrt2*Zh12*Zh21
     &      *Zh23*aa6+2D0*Sqrt2*Zh12*Zh22*Zh23*aa2-Sqrt2
     &      *Zh13*Zh21*Zh22*aa5-Sqrt2*Zh13*Zh21*Zh22*aa6
     &      +Sqrt2*Zh13*Zh21**2*la1s*xvev+Sqrt2*Zh13*Zh22
     &      **2*la2s*xvev+6D0*Sqrt2*Zh13*Zh23**2*la3s*xvev
     &      +12D0*Sqrt2*Zh13*Zh23**2*la7s*xvev+12D0*Sqrt2
     &      *Zh13*Zh23**2*la8s*xvev+2D0*Sqrt2*Zh11*Zh21*Zh23
     &      *la1s*xvev-2D0*Sqrt2*Zh11*Zh22*Zh23*la4s*xvev
     &      -2D0*Sqrt2*Zh11*Zh22*Zh23*la5s*xvev-2D0*Sqrt2
     &      *Zh11*Zh22*Zh23*la6s*xvev-2D0*Sqrt2*Zh12*Zh21
     &      *Zh23*la4s*xvev-2D0*Sqrt2*Zh12*Zh21*Zh23*la5s
     &      *xvev-2D0*Sqrt2*Zh12*Zh21*Zh23*la6s*xvev+2D0
     &      *Sqrt2*Zh12*Zh22*Zh23*la2s*xvev-2D0*Sqrt2*Zh13
     &      *Zh21*Zh22*la4s*xvev-2D0*Sqrt2*Zh13*Zh21*Zh22
     &      *la5s*xvev-2D0*Sqrt2*Zh13*Zh21*Zh22*la6s*xvev
     &      +6D0*MW/EE*SW*Zh11*Zh21**2*cb*la1-6D0*MW/EE*SW
     &      *Zh11*Zh21**2*la6*sb+2D0*MW/EE*SW*Zh11*Zh22**2
     &      *cb*la3+2D0*MW/EE*SW*Zh11*Zh22**2*cb*la4+2D0
     &      *MW/EE*SW*Zh11*Zh22**2*cb*la5-6D0*MW/EE*SW*Zh11
     &      *Zh22**2*la7*sb+2D0*MW/EE*SW*Zh11*Zh23**2*cb
     &      *la1s-2D0*MW/EE*SW*Zh11*Zh23**2*la4s*sb-2D0*MW
     &      /EE*SW*Zh11*Zh23**2*la5s*sb-2D0*MW/EE*SW*Zh11
     &      *Zh23**2*la6s*sb-6D0*MW/EE*SW*Zh12*Zh21**2*cb
     &      *la6+2D0*MW/EE*SW*Zh12*Zh21**2*la3*sb+2D0*MW
     &      /EE*SW*Zh12*Zh21**2*la4*sb+2D0*MW/EE*SW*Zh12
     &      *Zh21**2*la5*sb-6D0*MW/EE*SW*Zh12*Zh22**2*cb
     &      *la7+6D0*MW/EE*SW*Zh12*Zh22**2*la2*sb-2D0*MW
     &      /EE*SW*Zh12*Zh23**2*cb*la4s-2D0*MW/EE*SW*Zh12
     &      *Zh23**2*cb*la5s-2D0*MW/EE*SW*Zh12*Zh23**2*cb
     &      *la6s+2D0*MW/EE*SW*Zh12*Zh23**2*la2s*sb-12D0
     &      *MW/EE*SW*Zh11*Zh21*Zh22*cb*la6+4D0*MW/EE*SW
     &      *Zh11*Zh21*Zh22*la3*sb+4D0*MW/EE*SW*Zh11*Zh21
     &      *Zh22*la4*sb+4D0*MW/EE*SW*Zh11*Zh21*Zh22*la5
     &      *sb+4D0*MW/EE*SW*Zh12*Zh21*Zh22*cb*la3+4D0*MW
     &      /EE*SW*Zh12*Zh21*Zh22*cb*la4+4D0*MW/EE*SW*Zh12
     &      *Zh21*Zh22*cb*la5-12D0*MW/EE*SW*Zh12*Zh21*Zh22
     &      *la7*sb+4D0*MW/EE*SW*Zh13*Zh21*Zh23*cb*la1s-4D0
     &      *MW/EE*SW*Zh13*Zh21*Zh23*la4s*sb-4D0*MW/EE*SW
     &      *Zh13*Zh21*Zh23*la5s*sb-4D0*MW/EE*SW*Zh13*Zh21
     &      *Zh23*la6s*sb-4D0*MW/EE*SW*Zh13*Zh22*Zh23*cb
     &      *la4s-4D0*MW/EE*SW*Zh13*Zh22*Zh23*cb*la5s-4D0
     &      *MW/EE*SW*Zh13*Zh22*Zh23*cb*la6s+4D0*MW/EE*SW
     &      *Zh13*Zh22*Zh23*la2s*sb
      AAABR(167) = 2D0*Sqrt2*Zh11*Zh21*Zh33*aa1-Sqrt2*Zh11
     &      *Zh22*Zh33*aa5-Sqrt2*Zh11*Zh22*Zh33*aa6+2D0*Sqrt2
     &      *Zh11*Zh23*Zh31*aa1-Sqrt2*Zh11*Zh23*Zh32*aa5
     &      -Sqrt2*Zh11*Zh23*Zh32*aa6-Sqrt2*Zh12*Zh21*Zh33
     &      *aa5-Sqrt2*Zh12*Zh21*Zh33*aa6+2D0*Sqrt2*Zh12
     &      *Zh22*Zh33*aa2-Sqrt2*Zh12*Zh23*Zh31*aa5-Sqrt2
     &      *Zh12*Zh23*Zh31*aa6+2D0*Sqrt2*Zh12*Zh23*Zh32
     &      *aa2+2D0*Sqrt2*Zh13*Zh21*Zh31*aa1-Sqrt2*Zh13
     &      *Zh21*Zh32*aa5-Sqrt2*Zh13*Zh21*Zh32*aa6-Sqrt2
     &      *Zh13*Zh22*Zh31*aa5-Sqrt2*Zh13*Zh22*Zh31*aa6
     &      +2D0*Sqrt2*Zh13*Zh22*Zh32*aa2+6D0*Sqrt2*Zh13
     &      *Zh23*Zh33*aa3+6D0*Sqrt2*Zh13*Zh23*Zh33*aa4+2D0
     &      *Sqrt2*Zh11*Zh21*Zh33*la1s*xvev-2D0*Sqrt2*Zh11
     &      *Zh22*Zh33*la4s*xvev-2D0*Sqrt2*Zh11*Zh22*Zh33
     &      *la5s*xvev-2D0*Sqrt2*Zh11*Zh22*Zh33*la6s*xvev
     &      +2D0*Sqrt2*Zh11*Zh23*Zh31*la1s*xvev-2D0*Sqrt2
     &      *Zh11*Zh23*Zh32*la4s*xvev-2D0*Sqrt2*Zh11*Zh23
     &      *Zh32*la5s*xvev-2D0*Sqrt2*Zh11*Zh23*Zh32*la6s
     &      *xvev-2D0*Sqrt2*Zh12*Zh21*Zh33*la4s*xvev-2D0
     &      *Sqrt2*Zh12*Zh21*Zh33*la5s*xvev-2D0*Sqrt2*Zh12
     &      *Zh21*Zh33*la6s*xvev+2D0*Sqrt2*Zh12*Zh22*Zh33
     &      *la2s*xvev-2D0*Sqrt2*Zh12*Zh23*Zh31*la4s*xvev
     &      -2D0*Sqrt2*Zh12*Zh23*Zh31*la5s*xvev-2D0*Sqrt2
     &      *Zh12*Zh23*Zh31*la6s*xvev+2D0*Sqrt2*Zh12*Zh23
     &      *Zh32*la2s*xvev+2D0*Sqrt2*Zh13*Zh21*Zh31*la1s
     &      *xvev-2D0*Sqrt2*Zh13*Zh21*Zh32*la4s*xvev-2D0
     &      *Sqrt2*Zh13*Zh21*Zh32*la5s*xvev-2D0*Sqrt2*Zh13
     &      *Zh21*Zh32*la6s*xvev-2D0*Sqrt2*Zh13*Zh22*Zh31
     &      *la4s*xvev-2D0*Sqrt2*Zh13*Zh22*Zh31*la5s*xvev
     &      -2D0*Sqrt2*Zh13*Zh22*Zh31*la6s*xvev+2D0*Sqrt2
     &      *Zh13*Zh22*Zh32*la2s*xvev+12D0*Sqrt2*Zh13*Zh23
     &      *Zh33*la3s*xvev+24D0*Sqrt2*Zh13*Zh23*Zh33*la7s
     &      *xvev+24D0*Sqrt2*Zh13*Zh23*Zh33*la8s*xvev+12D0
     &      *MW/EE*SW*Zh11*Zh21*Zh31*cb*la1-12D0*MW/EE*SW
     &      *Zh11*Zh21*Zh31*la6*sb-12D0*MW/EE*SW*Zh11*Zh21
     &      *Zh32*cb*la6+4D0*MW/EE*SW*Zh11*Zh21*Zh32*la3
     &      *sb+4D0*MW/EE*SW*Zh11*Zh21*Zh32*la4*sb+4D0*MW
     &      /EE*SW*Zh11*Zh21*Zh32*la5*sb-12D0*MW/EE*SW*Zh11
     &      *Zh22*Zh31*cb*la6+4D0*MW/EE*SW*Zh11*Zh22*Zh31
     &      *la3*sb+4D0*MW/EE*SW*Zh11*Zh22*Zh31*la4*sb+4D0
     &      *MW/EE*SW*Zh11*Zh22*Zh31*la5*sb+4D0*MW/EE*SW
     &      *Zh11*Zh22*Zh32*cb*la3+4D0*MW/EE*SW*Zh11*Zh22
     &      *Zh32*cb*la4+4D0*MW/EE*SW*Zh11*Zh22*Zh32*cb*la5
     &      -12D0*MW/EE*SW*Zh11*Zh22*Zh32*la7*sb+4D0*MW/EE
     &      *SW*Zh11*Zh23*Zh33*cb*la1s-4D0*MW/EE*SW*Zh11
     &      *Zh23*Zh33*la4s*sb-4D0*MW/EE*SW*Zh11*Zh23*Zh33
     &      *la5s*sb-4D0*MW/EE*SW*Zh11*Zh23*Zh33*la6s*sb
     &      -12D0*MW/EE*SW*Zh12*Zh21*Zh31*cb*la6+4D0*MW/EE
     &      *SW*Zh12*Zh21*Zh31*la3*sb+4D0*MW/EE*SW*Zh12*Zh21
     &      *Zh31*la4*sb+4D0*MW/EE*SW*Zh12*Zh21*Zh31*la5
     &      *sb+4D0*MW/EE*SW*Zh12*Zh21*Zh32*cb*la3+4D0*MW
     &      /EE*SW*Zh12*Zh21*Zh32*cb*la4+4D0*MW/EE*SW*Zh12
     &      *Zh21*Zh32*cb*la5-12D0*MW/EE*SW*Zh12*Zh21*Zh32
     &      *la7*sb+4D0*MW/EE*SW*Zh12*Zh22*Zh31*cb*la3+4D0
     &      *MW/EE*SW*Zh12*Zh22*Zh31*cb*la4+4D0*MW/EE*SW
     &      *Zh12*Zh22*Zh31*cb*la5-12D0*MW/EE*SW*Zh12*Zh22
     &      *Zh31*la7*sb-12D0*MW/EE*SW*Zh12*Zh22*Zh32*cb
     &      *la7+12D0*MW/EE*SW*Zh12*Zh22*Zh32*la2*sb-4D0
     &      *MW/EE*SW*Zh12*Zh23*Zh33*cb*la4s-4D0*MW/EE*SW
     &      *Zh12*Zh23*Zh33*cb*la5s-4D0*MW/EE*SW*Zh12*Zh23
     &      *Zh33*cb*la6s+4D0*MW/EE*SW*Zh12*Zh23*Zh33*la2s
     &      *sb+4D0*MW/EE*SW*Zh13*Zh21*Zh33*cb*la1s-4D0*MW
     &      /EE*SW*Zh13*Zh21*Zh33*la4s*sb-4D0*MW/EE*SW*Zh13
     &      *Zh21*Zh33*la5s*sb-4D0*MW/EE*SW*Zh13*Zh21*Zh33
     &      *la6s*sb-4D0*MW/EE*SW*Zh13*Zh22*Zh33*cb*la4s
     &      -4D0*MW/EE*SW*Zh13*Zh22*Zh33*cb*la5s-4D0*MW/EE
     &      *SW*Zh13*Zh22*Zh33*cb*la6s+4D0*MW/EE*SW*Zh13
     &      *Zh22*Zh33*la2s*sb+4D0*MW/EE*SW*Zh13*Zh23*Zh31
     &      *cb*la1s-4D0*MW/EE*SW*Zh13*Zh23*Zh31*la4s*sb
     &      -4D0*MW/EE*SW*Zh13*Zh23*Zh31*la5s*sb-4D0*MW/EE
     &      *SW*Zh13*Zh23*Zh31*la6s*sb-4D0*MW/EE*SW*Zh13
     &      *Zh23*Zh32*cb*la4s-4D0*MW/EE*SW*Zh13*Zh23*Zh32
     &      *cb*la5s-4D0*MW/EE*SW*Zh13*Zh23*Zh32*cb*la6s
     &      +4D0*MW/EE*SW*Zh13*Zh23*Zh32*la2s*sb
      AAABR(168) = Sqrt2*Zh13*Zh31**2*aa1+Sqrt2*Zh13*Zh32
     &      **2*aa2+3D0*Sqrt2*Zh13*Zh33**2*aa3+3D0*Sqrt2
     &      *Zh13*Zh33**2*aa4+2D0*Sqrt2*Zh11*Zh31*Zh33*aa1
     &      -Sqrt2*Zh11*Zh32*Zh33*aa5-Sqrt2*Zh11*Zh32*Zh33
     &      *aa6-Sqrt2*Zh12*Zh31*Zh33*aa5-Sqrt2*Zh12*Zh31
     &      *Zh33*aa6+2D0*Sqrt2*Zh12*Zh32*Zh33*aa2-Sqrt2
     &      *Zh13*Zh31*Zh32*aa5-Sqrt2*Zh13*Zh31*Zh32*aa6
     &      +Sqrt2*Zh13*Zh31**2*la1s*xvev+Sqrt2*Zh13*Zh32
     &      **2*la2s*xvev+6D0*Sqrt2*Zh13*Zh33**2*la3s*xvev
     &      +12D0*Sqrt2*Zh13*Zh33**2*la7s*xvev+12D0*Sqrt2
     &      *Zh13*Zh33**2*la8s*xvev+2D0*Sqrt2*Zh11*Zh31*Zh33
     &      *la1s*xvev-2D0*Sqrt2*Zh11*Zh32*Zh33*la4s*xvev
     &      -2D0*Sqrt2*Zh11*Zh32*Zh33*la5s*xvev-2D0*Sqrt2
     &      *Zh11*Zh32*Zh33*la6s*xvev-2D0*Sqrt2*Zh12*Zh31
     &      *Zh33*la4s*xvev-2D0*Sqrt2*Zh12*Zh31*Zh33*la5s
     &      *xvev-2D0*Sqrt2*Zh12*Zh31*Zh33*la6s*xvev+2D0
     &      *Sqrt2*Zh12*Zh32*Zh33*la2s*xvev-2D0*Sqrt2*Zh13
     &      *Zh31*Zh32*la4s*xvev-2D0*Sqrt2*Zh13*Zh31*Zh32
     &      *la5s*xvev-2D0*Sqrt2*Zh13*Zh31*Zh32*la6s*xvev
     &      +6D0*MW/EE*SW*Zh11*Zh31**2*cb*la1-6D0*MW/EE*SW
     &      *Zh11*Zh31**2*la6*sb+2D0*MW/EE*SW*Zh11*Zh32**2
     &      *cb*la3+2D0*MW/EE*SW*Zh11*Zh32**2*cb*la4+2D0
     &      *MW/EE*SW*Zh11*Zh32**2*cb*la5-6D0*MW/EE*SW*Zh11
     &      *Zh32**2*la7*sb+2D0*MW/EE*SW*Zh11*Zh33**2*cb
     &      *la1s-2D0*MW/EE*SW*Zh11*Zh33**2*la4s*sb-2D0*MW
     &      /EE*SW*Zh11*Zh33**2*la5s*sb-2D0*MW/EE*SW*Zh11
     &      *Zh33**2*la6s*sb-6D0*MW/EE*SW*Zh12*Zh31**2*cb
     &      *la6+2D0*MW/EE*SW*Zh12*Zh31**2*la3*sb+2D0*MW
     &      /EE*SW*Zh12*Zh31**2*la4*sb+2D0*MW/EE*SW*Zh12
     &      *Zh31**2*la5*sb-6D0*MW/EE*SW*Zh12*Zh32**2*cb
     &      *la7+6D0*MW/EE*SW*Zh12*Zh32**2*la2*sb-2D0*MW
     &      /EE*SW*Zh12*Zh33**2*cb*la4s-2D0*MW/EE*SW*Zh12
     &      *Zh33**2*cb*la5s-2D0*MW/EE*SW*Zh12*Zh33**2*cb
     &      *la6s+2D0*MW/EE*SW*Zh12*Zh33**2*la2s*sb-12D0
     &      *MW/EE*SW*Zh11*Zh31*Zh32*cb*la6+4D0*MW/EE*SW
     &      *Zh11*Zh31*Zh32*la3*sb+4D0*MW/EE*SW*Zh11*Zh31
     &      *Zh32*la4*sb+4D0*MW/EE*SW*Zh11*Zh31*Zh32*la5
     &      *sb+4D0*MW/EE*SW*Zh12*Zh31*Zh32*cb*la3+4D0*MW
     &      /EE*SW*Zh12*Zh31*Zh32*cb*la4+4D0*MW/EE*SW*Zh12
     &      *Zh31*Zh32*cb*la5-12D0*MW/EE*SW*Zh12*Zh31*Zh32
     &      *la7*sb+4D0*MW/EE*SW*Zh13*Zh31*Zh33*cb*la1s-4D0
     &      *MW/EE*SW*Zh13*Zh31*Zh33*la4s*sb-4D0*MW/EE*SW
     &      *Zh13*Zh31*Zh33*la5s*sb-4D0*MW/EE*SW*Zh13*Zh31
     &      *Zh33*la6s*sb-4D0*MW/EE*SW*Zh13*Zh32*Zh33*cb
     &      *la4s-4D0*MW/EE*SW*Zh13*Zh32*Zh33*cb*la5s-4D0
     &      *MW/EE*SW*Zh13*Zh32*Zh33*cb*la6s+4D0*MW/EE*SW
     &      *Zh13*Zh32*Zh33*la2s*sb
      AAABR(169) = Sqrt2*Zh23*Zh31**2*aa1+Sqrt2*Zh23*Zh32
     &      **2*aa2+3D0*Sqrt2*Zh23*Zh33**2*aa3+3D0*Sqrt2
     &      *Zh23*Zh33**2*aa4+2D0*Sqrt2*Zh21*Zh31*Zh33*aa1
     &      -Sqrt2*Zh21*Zh32*Zh33*aa5-Sqrt2*Zh21*Zh32*Zh33
     &      *aa6-Sqrt2*Zh22*Zh31*Zh33*aa5-Sqrt2*Zh22*Zh31
     &      *Zh33*aa6+2D0*Sqrt2*Zh22*Zh32*Zh33*aa2-Sqrt2
     &      *Zh23*Zh31*Zh32*aa5-Sqrt2*Zh23*Zh31*Zh32*aa6
     &      +Sqrt2*Zh23*Zh31**2*la1s*xvev+Sqrt2*Zh23*Zh32
     &      **2*la2s*xvev+6D0*Sqrt2*Zh23*Zh33**2*la3s*xvev
     &      +12D0*Sqrt2*Zh23*Zh33**2*la7s*xvev+12D0*Sqrt2
     &      *Zh23*Zh33**2*la8s*xvev+2D0*Sqrt2*Zh21*Zh31*Zh33
     &      *la1s*xvev-2D0*Sqrt2*Zh21*Zh32*Zh33*la4s*xvev
     &      -2D0*Sqrt2*Zh21*Zh32*Zh33*la5s*xvev-2D0*Sqrt2
     &      *Zh21*Zh32*Zh33*la6s*xvev-2D0*Sqrt2*Zh22*Zh31
     &      *Zh33*la4s*xvev-2D0*Sqrt2*Zh22*Zh31*Zh33*la5s
     &      *xvev-2D0*Sqrt2*Zh22*Zh31*Zh33*la6s*xvev+2D0
     &      *Sqrt2*Zh22*Zh32*Zh33*la2s*xvev-2D0*Sqrt2*Zh23
     &      *Zh31*Zh32*la4s*xvev-2D0*Sqrt2*Zh23*Zh31*Zh32
     &      *la5s*xvev-2D0*Sqrt2*Zh23*Zh31*Zh32*la6s*xvev
     &      +6D0*MW/EE*SW*Zh21*Zh31**2*cb*la1-6D0*MW/EE*SW
     &      *Zh21*Zh31**2*la6*sb+2D0*MW/EE*SW*Zh21*Zh32**2
     &      *cb*la3+2D0*MW/EE*SW*Zh21*Zh32**2*cb*la4+2D0
     &      *MW/EE*SW*Zh21*Zh32**2*cb*la5-6D0*MW/EE*SW*Zh21
     &      *Zh32**2*la7*sb+2D0*MW/EE*SW*Zh21*Zh33**2*cb
     &      *la1s-2D0*MW/EE*SW*Zh21*Zh33**2*la4s*sb-2D0*MW
     &      /EE*SW*Zh21*Zh33**2*la5s*sb-2D0*MW/EE*SW*Zh21
     &      *Zh33**2*la6s*sb-6D0*MW/EE*SW*Zh22*Zh31**2*cb
     &      *la6+2D0*MW/EE*SW*Zh22*Zh31**2*la3*sb+2D0*MW
     &      /EE*SW*Zh22*Zh31**2*la4*sb+2D0*MW/EE*SW*Zh22
     &      *Zh31**2*la5*sb-6D0*MW/EE*SW*Zh22*Zh32**2*cb
     &      *la7+6D0*MW/EE*SW*Zh22*Zh32**2*la2*sb-2D0*MW
     &      /EE*SW*Zh22*Zh33**2*cb*la4s-2D0*MW/EE*SW*Zh22
     &      *Zh33**2*cb*la5s-2D0*MW/EE*SW*Zh22*Zh33**2*cb
     &      *la6s+2D0*MW/EE*SW*Zh22*Zh33**2*la2s*sb-12D0
     &      *MW/EE*SW*Zh21*Zh31*Zh32*cb*la6+4D0*MW/EE*SW
     &      *Zh21*Zh31*Zh32*la3*sb+4D0*MW/EE*SW*Zh21*Zh31
     &      *Zh32*la4*sb+4D0*MW/EE*SW*Zh21*Zh31*Zh32*la5
     &      *sb+4D0*MW/EE*SW*Zh22*Zh31*Zh32*cb*la3+4D0*MW
     &      /EE*SW*Zh22*Zh31*Zh32*cb*la4+4D0*MW/EE*SW*Zh22
     &      *Zh31*Zh32*cb*la5-12D0*MW/EE*SW*Zh22*Zh31*Zh32
     &      *la7*sb+4D0*MW/EE*SW*Zh23*Zh31*Zh33*cb*la1s-4D0
     &      *MW/EE*SW*Zh23*Zh31*Zh33*la4s*sb-4D0*MW/EE*SW
     &      *Zh23*Zh31*Zh33*la5s*sb-4D0*MW/EE*SW*Zh23*Zh31
     &      *Zh33*la6s*sb-4D0*MW/EE*SW*Zh23*Zh32*Zh33*cb
     &      *la4s-4D0*MW/EE*SW*Zh23*Zh32*Zh33*cb*la5s-4D0
     &      *MW/EE*SW*Zh23*Zh32*Zh33*cb*la6s+4D0*MW/EE*SW
     &      *Zh23*Zh32*Zh33*la2s*sb
      AAABR(170) = Sqrt2*Zh33**3*aa3+Sqrt2*Zh33**3*aa4+Sqrt2
     &      *Zh31**2*Zh33*aa1+Sqrt2*Zh32**2*Zh33*aa2+2D0
     &      *Sqrt2*Zh33**3*la3s*xvev+4D0*Sqrt2*Zh33**3*la7s
     &      *xvev+4D0*Sqrt2*Zh33**3*la8s*xvev-Sqrt2*Zh31
     &      *Zh32*Zh33*aa5-Sqrt2*Zh31*Zh32*Zh33*aa6+Sqrt2
     &      *Zh31**2*Zh33*la1s*xvev+Sqrt2*Zh32**2*Zh33*la2s
     &      *xvev+2D0*MW/EE*SW*Zh31**3*cb*la1-2D0*MW/EE*SW
     &      *Zh31**3*la6*sb-2D0*MW/EE*SW*Zh32**3*cb*la7+2D0
     &      *MW/EE*SW*Zh32**3*la2*sb-2D0*Sqrt2*Zh31*Zh32
     &      *Zh33*la4s*xvev-2D0*Sqrt2*Zh31*Zh32*Zh33*la5s
     &      *xvev-2D0*Sqrt2*Zh31*Zh32*Zh33*la6s*xvev+2D0
     &      *MW/EE*SW*Zh31*Zh32**2*cb*la3+2D0*MW/EE*SW*Zh31
     &      *Zh32**2*cb*la4+2D0*MW/EE*SW*Zh31*Zh32**2*cb
     &      *la5-6D0*MW/EE*SW*Zh31*Zh32**2*la7*sb+2D0*MW
     &      /EE*SW*Zh31*Zh33**2*cb*la1s-2D0*MW/EE*SW*Zh31
     &      *Zh33**2*la4s*sb-2D0*MW/EE*SW*Zh31*Zh33**2*la5s
     &      *sb-2D0*MW/EE*SW*Zh31*Zh33**2*la6s*sb-6D0*MW
     &      /EE*SW*Zh31**2*Zh32*cb*la6+2D0*MW/EE*SW*Zh31
     &      **2*Zh32*la3*sb+2D0*MW/EE*SW*Zh31**2*Zh32*la4
     &      *sb+2D0*MW/EE*SW*Zh31**2*Zh32*la5*sb-2D0*MW/EE
     &      *SW*Zh32*Zh33**2*cb*la4s-2D0*MW/EE*SW*Zh32*Zh33
     &      **2*cb*la5s-2D0*MW/EE*SW*Zh32*Zh33**2*cb*la6s
     &      +2D0*MW/EE*SW*Zh32*Zh33**2*la2s*sb
      AAABR(171) = Sqrt2*Zh13**3*aa3+Sqrt2*Zh13**3*aa4+Sqrt2
     &      *Zh11**2*Zh13*aa1+Sqrt2*Zh12**2*Zh13*aa2+2D0
     &      *Sqrt2*Zh13**3*la3s*xvev+4D0*Sqrt2*Zh13**3*la7s
     &      *xvev+4D0*Sqrt2*Zh13**3*la8s*xvev-Sqrt2*Zh11
     &      *Zh12*Zh13*aa5-Sqrt2*Zh11*Zh12*Zh13*aa6+Sqrt2
     &      *Zh11**2*Zh13*la1s*xvev+Sqrt2*Zh12**2*Zh13*la2s
     &      *xvev+2D0*MW/EE*SW*Zh11**3*cb*la1-2D0*MW/EE*SW
     &      *Zh11**3*la6*sb-2D0*MW/EE*SW*Zh12**3*cb*la7+2D0
     &      *MW/EE*SW*Zh12**3*la2*sb-2D0*Sqrt2*Zh11*Zh12
     &      *Zh13*la4s*xvev-2D0*Sqrt2*Zh11*Zh12*Zh13*la5s
     &      *xvev-2D0*Sqrt2*Zh11*Zh12*Zh13*la6s*xvev+2D0
     &      *MW/EE*SW*Zh11*Zh12**2*cb*la3+2D0*MW/EE*SW*Zh11
     &      *Zh12**2*cb*la4+2D0*MW/EE*SW*Zh11*Zh12**2*cb
     &      *la5-6D0*MW/EE*SW*Zh11*Zh12**2*la7*sb+2D0*MW
     &      /EE*SW*Zh11*Zh13**2*cb*la1s-2D0*MW/EE*SW*Zh11
     &      *Zh13**2*la4s*sb-2D0*MW/EE*SW*Zh11*Zh13**2*la5s
     &      *sb-2D0*MW/EE*SW*Zh11*Zh13**2*la6s*sb-6D0*MW
     &      /EE*SW*Zh11**2*Zh12*cb*la6+2D0*MW/EE*SW*Zh11
     &      **2*Zh12*la3*sb+2D0*MW/EE*SW*Zh11**2*Zh12*la4
     &      *sb+2D0*MW/EE*SW*Zh11**2*Zh12*la5*sb-2D0*MW/EE
     &      *SW*Zh12*Zh13**2*cb*la4s-2D0*MW/EE*SW*Zh12*Zh13
     &      **2*cb*la5s-2D0*MW/EE*SW*Zh12*Zh13**2*cb*la6s
     &      +2D0*MW/EE*SW*Zh12*Zh13**2*la2s*sb
      AAABR(172) = Sqrt2*Zh21**2*Zh33*aa1+Sqrt2*Zh22**2*Zh33
     &      *aa2+3D0*Sqrt2*Zh23**2*Zh33*aa3+3D0*Sqrt2*Zh23
     &      **2*Zh33*aa4-Sqrt2*Zh21*Zh22*Zh33*aa5-Sqrt2*Zh21
     &      *Zh22*Zh33*aa6+2D0*Sqrt2*Zh21*Zh23*Zh31*aa1-Sqrt2
     &      *Zh21*Zh23*Zh32*aa5-Sqrt2*Zh21*Zh23*Zh32*aa6
     &      +Sqrt2*Zh21**2*Zh33*la1s*xvev-Sqrt2*Zh22*Zh23
     &      *Zh31*aa5-Sqrt2*Zh22*Zh23*Zh31*aa6+2D0*Sqrt2
     &      *Zh22*Zh23*Zh32*aa2+Sqrt2*Zh22**2*Zh33*la2s*xvev
     &      +6D0*Sqrt2*Zh23**2*Zh33*la3s*xvev+12D0*Sqrt2
     &      *Zh23**2*Zh33*la7s*xvev+12D0*Sqrt2*Zh23**2*Zh33
     &      *la8s*xvev-2D0*Sqrt2*Zh21*Zh22*Zh33*la4s*xvev
     &      -2D0*Sqrt2*Zh21*Zh22*Zh33*la5s*xvev-2D0*Sqrt2
     &      *Zh21*Zh22*Zh33*la6s*xvev+2D0*Sqrt2*Zh21*Zh23
     &      *Zh31*la1s*xvev-2D0*Sqrt2*Zh21*Zh23*Zh32*la4s
     &      *xvev-2D0*Sqrt2*Zh21*Zh23*Zh32*la5s*xvev-2D0
     &      *Sqrt2*Zh21*Zh23*Zh32*la6s*xvev-2D0*Sqrt2*Zh22
     &      *Zh23*Zh31*la4s*xvev-2D0*Sqrt2*Zh22*Zh23*Zh31
     &      *la5s*xvev-2D0*Sqrt2*Zh22*Zh23*Zh31*la6s*xvev
     &      +2D0*Sqrt2*Zh22*Zh23*Zh32*la2s*xvev+6D0*MW/EE
     &      *SW*Zh21**2*Zh31*cb*la1-6D0*MW/EE*SW*Zh21**2
     &      *Zh31*la6*sb-6D0*MW/EE*SW*Zh21**2*Zh32*cb*la6
     &      +2D0*MW/EE*SW*Zh21**2*Zh32*la3*sb+2D0*MW/EE*SW
     &      *Zh21**2*Zh32*la4*sb+2D0*MW/EE*SW*Zh21**2*Zh32
     &      *la5*sb+2D0*MW/EE*SW*Zh22**2*Zh31*cb*la3+2D0
     &      *MW/EE*SW*Zh22**2*Zh31*cb*la4+2D0*MW/EE*SW*Zh22
     &      **2*Zh31*cb*la5-6D0*MW/EE*SW*Zh22**2*Zh31*la7
     &      *sb-6D0*MW/EE*SW*Zh22**2*Zh32*cb*la7+6D0*MW/EE
     &      *SW*Zh22**2*Zh32*la2*sb+2D0*MW/EE*SW*Zh23**2
     &      *Zh31*cb*la1s-2D0*MW/EE*SW*Zh23**2*Zh31*la4s
     &      *sb-2D0*MW/EE*SW*Zh23**2*Zh31*la5s*sb-2D0*MW
     &      /EE*SW*Zh23**2*Zh31*la6s*sb-2D0*MW/EE*SW*Zh23
     &      **2*Zh32*cb*la4s-2D0*MW/EE*SW*Zh23**2*Zh32*cb
     &      *la5s-2D0*MW/EE*SW*Zh23**2*Zh32*cb*la6s+2D0*MW
     &      /EE*SW*Zh23**2*Zh32*la2s*sb-12D0*MW/EE*SW*Zh21
     &      *Zh22*Zh31*cb*la6+4D0*MW/EE*SW*Zh21*Zh22*Zh31
     &      *la3*sb+4D0*MW/EE*SW*Zh21*Zh22*Zh31*la4*sb+4D0
     &      *MW/EE*SW*Zh21*Zh22*Zh31*la5*sb+4D0*MW/EE*SW
     &      *Zh21*Zh22*Zh32*cb*la3+4D0*MW/EE*SW*Zh21*Zh22
     &      *Zh32*cb*la4+4D0*MW/EE*SW*Zh21*Zh22*Zh32*cb*la5
     &      -12D0*MW/EE*SW*Zh21*Zh22*Zh32*la7*sb+4D0*MW/EE
     &      *SW*Zh21*Zh23*Zh33*cb*la1s-4D0*MW/EE*SW*Zh21
     &      *Zh23*Zh33*la4s*sb-4D0*MW/EE*SW*Zh21*Zh23*Zh33
     &      *la5s*sb-4D0*MW/EE*SW*Zh21*Zh23*Zh33*la6s*sb
     &      -4D0*MW/EE*SW*Zh22*Zh23*Zh33*cb*la4s-4D0*MW/EE
     &      *SW*Zh22*Zh23*Zh33*cb*la5s-4D0*MW/EE*SW*Zh22
     &      *Zh23*Zh33*cb*la6s+4D0*MW/EE*SW*Zh22*Zh23*Zh33
     &      *la2s*sb
      AAABR(173) = Sqrt2*Zh11**2*Zh33*aa1+Sqrt2*Zh12**2*Zh33
     &      *aa2+3D0*Sqrt2*Zh13**2*Zh33*aa3+3D0*Sqrt2*Zh13
     &      **2*Zh33*aa4-Sqrt2*Zh11*Zh12*Zh33*aa5-Sqrt2*Zh11
     &      *Zh12*Zh33*aa6+2D0*Sqrt2*Zh11*Zh13*Zh31*aa1-Sqrt2
     &      *Zh11*Zh13*Zh32*aa5-Sqrt2*Zh11*Zh13*Zh32*aa6
     &      +Sqrt2*Zh11**2*Zh33*la1s*xvev-Sqrt2*Zh12*Zh13
     &      *Zh31*aa5-Sqrt2*Zh12*Zh13*Zh31*aa6+2D0*Sqrt2
     &      *Zh12*Zh13*Zh32*aa2+Sqrt2*Zh12**2*Zh33*la2s*xvev
     &      +6D0*Sqrt2*Zh13**2*Zh33*la3s*xvev+12D0*Sqrt2
     &      *Zh13**2*Zh33*la7s*xvev+12D0*Sqrt2*Zh13**2*Zh33
     &      *la8s*xvev-2D0*Sqrt2*Zh11*Zh12*Zh33*la4s*xvev
     &      -2D0*Sqrt2*Zh11*Zh12*Zh33*la5s*xvev-2D0*Sqrt2
     &      *Zh11*Zh12*Zh33*la6s*xvev+2D0*Sqrt2*Zh11*Zh13
     &      *Zh31*la1s*xvev-2D0*Sqrt2*Zh11*Zh13*Zh32*la4s
     &      *xvev-2D0*Sqrt2*Zh11*Zh13*Zh32*la5s*xvev-2D0
     &      *Sqrt2*Zh11*Zh13*Zh32*la6s*xvev-2D0*Sqrt2*Zh12
     &      *Zh13*Zh31*la4s*xvev-2D0*Sqrt2*Zh12*Zh13*Zh31
     &      *la5s*xvev-2D0*Sqrt2*Zh12*Zh13*Zh31*la6s*xvev
     &      +2D0*Sqrt2*Zh12*Zh13*Zh32*la2s*xvev+6D0*MW/EE
     &      *SW*Zh11**2*Zh31*cb*la1-6D0*MW/EE*SW*Zh11**2
     &      *Zh31*la6*sb-6D0*MW/EE*SW*Zh11**2*Zh32*cb*la6
     &      +2D0*MW/EE*SW*Zh11**2*Zh32*la3*sb+2D0*MW/EE*SW
     &      *Zh11**2*Zh32*la4*sb+2D0*MW/EE*SW*Zh11**2*Zh32
     &      *la5*sb+2D0*MW/EE*SW*Zh12**2*Zh31*cb*la3+2D0
     &      *MW/EE*SW*Zh12**2*Zh31*cb*la4+2D0*MW/EE*SW*Zh12
     &      **2*Zh31*cb*la5-6D0*MW/EE*SW*Zh12**2*Zh31*la7
     &      *sb-6D0*MW/EE*SW*Zh12**2*Zh32*cb*la7+6D0*MW/EE
     &      *SW*Zh12**2*Zh32*la2*sb+2D0*MW/EE*SW*Zh13**2
     &      *Zh31*cb*la1s-2D0*MW/EE*SW*Zh13**2*Zh31*la4s
     &      *sb-2D0*MW/EE*SW*Zh13**2*Zh31*la5s*sb-2D0*MW
     &      /EE*SW*Zh13**2*Zh31*la6s*sb-2D0*MW/EE*SW*Zh13
     &      **2*Zh32*cb*la4s-2D0*MW/EE*SW*Zh13**2*Zh32*cb
     &      *la5s-2D0*MW/EE*SW*Zh13**2*Zh32*cb*la6s+2D0*MW
     &      /EE*SW*Zh13**2*Zh32*la2s*sb-12D0*MW/EE*SW*Zh11
     &      *Zh12*Zh31*cb*la6+4D0*MW/EE*SW*Zh11*Zh12*Zh31
     &      *la3*sb+4D0*MW/EE*SW*Zh11*Zh12*Zh31*la4*sb+4D0
     &      *MW/EE*SW*Zh11*Zh12*Zh31*la5*sb+4D0*MW/EE*SW
     &      *Zh11*Zh12*Zh32*cb*la3+4D0*MW/EE*SW*Zh11*Zh12
     &      *Zh32*cb*la4+4D0*MW/EE*SW*Zh11*Zh12*Zh32*cb*la5
     &      -12D0*MW/EE*SW*Zh11*Zh12*Zh32*la7*sb+4D0*MW/EE
     &      *SW*Zh11*Zh13*Zh33*cb*la1s-4D0*MW/EE*SW*Zh11
     &      *Zh13*Zh33*la4s*sb-4D0*MW/EE*SW*Zh11*Zh13*Zh33
     &      *la5s*sb-4D0*MW/EE*SW*Zh11*Zh13*Zh33*la6s*sb
     &      -4D0*MW/EE*SW*Zh12*Zh13*Zh33*cb*la4s-4D0*MW/EE
     &      *SW*Zh12*Zh13*Zh33*cb*la5s-4D0*MW/EE*SW*Zh12
     &      *Zh13*Zh33*cb*la6s+4D0*MW/EE*SW*Zh12*Zh13*Zh33
     &      *la2s*sb
      AAABR(174) = Sqrt2*Zh23**3*aa3+Sqrt2*Zh23**3*aa4+Sqrt2
     &      *Zh21**2*Zh23*aa1+Sqrt2*Zh22**2*Zh23*aa2+2D0
     &      *Sqrt2*Zh23**3*la3s*xvev+4D0*Sqrt2*Zh23**3*la7s
     &      *xvev+4D0*Sqrt2*Zh23**3*la8s*xvev-Sqrt2*Zh21
     &      *Zh22*Zh23*aa5-Sqrt2*Zh21*Zh22*Zh23*aa6+Sqrt2
     &      *Zh21**2*Zh23*la1s*xvev+Sqrt2*Zh22**2*Zh23*la2s
     &      *xvev+2D0*MW/EE*SW*Zh21**3*cb*la1-2D0*MW/EE*SW
     &      *Zh21**3*la6*sb-2D0*MW/EE*SW*Zh22**3*cb*la7+2D0
     &      *MW/EE*SW*Zh22**3*la2*sb-2D0*Sqrt2*Zh21*Zh22
     &      *Zh23*la4s*xvev-2D0*Sqrt2*Zh21*Zh22*Zh23*la5s
     &      *xvev-2D0*Sqrt2*Zh21*Zh22*Zh23*la6s*xvev+2D0
     &      *MW/EE*SW*Zh21*Zh22**2*cb*la3+2D0*MW/EE*SW*Zh21
     &      *Zh22**2*cb*la4+2D0*MW/EE*SW*Zh21*Zh22**2*cb
     &      *la5-6D0*MW/EE*SW*Zh21*Zh22**2*la7*sb+2D0*MW
     &      /EE*SW*Zh21*Zh23**2*cb*la1s-2D0*MW/EE*SW*Zh21
     &      *Zh23**2*la4s*sb-2D0*MW/EE*SW*Zh21*Zh23**2*la5s
     &      *sb-2D0*MW/EE*SW*Zh21*Zh23**2*la6s*sb-6D0*MW
     &      /EE*SW*Zh21**2*Zh22*cb*la6+2D0*MW/EE*SW*Zh21
     &      **2*Zh22*la3*sb+2D0*MW/EE*SW*Zh21**2*Zh22*la4
     &      *sb+2D0*MW/EE*SW*Zh21**2*Zh22*la5*sb-2D0*MW/EE
     &      *SW*Zh22*Zh23**2*cb*la4s-2D0*MW/EE*SW*Zh22*Zh23
     &      **2*cb*la5s-2D0*MW/EE*SW*Zh22*Zh23**2*cb*la6s
     &      +2D0*MW/EE*SW*Zh22*Zh23**2*la2s*sb
      AAABR(175) = EE*(SW/CW-CW/SW)
      AAABR(176) = EE*Pa21/SW
      AAABR(177) = EE*Pa11/SW
      AAABR(178) = EE/SW*(Zh31*sb-Zh32*cb)
      AAABR(179) = EE/SW*(Zh21*sb-Zh22*cb)
      AAABR(180) = EE/SW*(Zh11*sb-Zh12*cb)
      AAABR(181) = EE*(SW/CW*Zl11**2+2D0*SW/CW*Zl12**2-CW
     &      /SW*Zl11**2)
      AAABR(182) = EE*(SW/CW*Zl11*Zl21+2D0*SW/CW*Zl12*Zl22
     &      -CW/SW*Zl11*Zl21)
      AAABR(183) = EE/SW*Sqrt2*Zl11
      AAABR(184) = EE/SW*Sqrt2
      AAABR(185) = EE/CW*SW
      AAABR(186) = EE*(SW/CW*Zl21**2+2D0*SW/CW*Zl22**2-CW
     &      /SW*Zl21**2)
      AAABR(187) = EE/SW*Sqrt2*Zl21
      AAABR(188) = EE*(SW/CW+CW/SW)
      AAABR(189) = EE*(SW/CW*Zdd11**2-2D0*SW/CW*Zdd12**2
     &      +3D0*CW/SW*Zdd11**2)
      AAABR(190) = EE*(SW/CW*Zb11**2-2D0*SW/CW*Zb12**2+3D0
     &      *CW/SW*Zb11**2)
      AAABR(191) = EE*(SW/CW*Zss11**2-2D0*SW/CW*Zss12**2
     &      +3D0*CW/SW*Zss11**2)
      AAABR(192) = EE*(SW/CW*Zdd11*Zdd21-2D0*SW/CW*Zdd12
     &      *Zdd22+3D0*CW/SW*Zdd11*Zdd21)
      AAABR(193) = EE*(SW/CW*Zss11*Zss21-2D0*SW/CW*Zss12
     &      *Zss22+3D0*CW/SW*Zss11*Zss21)
      AAABR(194) = EE*(SW/CW*Zb11*Zb21-2D0*SW/CW*Zb12*Zb22
     &      +3D0*CW/SW*Zb11*Zb21)
      AAABR(195) = EE/SW*Sqrt2*Zb11*Zt11
      AAABR(196) = EE/SW*Sqrt2*Zdd11*Zuu11
      AAABR(197) = EE/SW*Sqrt2*Zcc11*Zss11
      AAABR(198) = EE/SW*Sqrt2*Zdd11*Zuu21
      AAABR(199) = EE/SW*Sqrt2*Zcc21*Zss11
      AAABR(200) = EE/SW*Sqrt2*Zb11*Zt21
      AAABR(201) = EE*(SW/CW*Zss21**2-2D0*SW/CW*Zss22**2
     &      +3D0*CW/SW*Zss21**2)
      AAABR(202) = EE*(SW/CW*Zb21**2-2D0*SW/CW*Zb22**2+3D0
     &      *CW/SW*Zb21**2)
      AAABR(203) = EE*(SW/CW*Zdd21**2-2D0*SW/CW*Zdd22**2
     &      +3D0*CW/SW*Zdd21**2)
      AAABR(204) = EE/SW*Sqrt2*Zcc11*Zss21
      AAABR(205) = EE/SW*Sqrt2*Zdd21*Zuu11
      AAABR(206) = EE/SW*Sqrt2*Zb21*Zt11
      AAABR(207) = EE/SW*Sqrt2*Zcc21*Zss21
      AAABR(208) = EE/SW*Sqrt2*Zb21*Zt21
      AAABR(209) = EE/SW*Sqrt2*Zdd21*Zuu21
      AAABR(210) = EE*(SW/CW*Zcc11**2+4D0*SW/CW*Zcc12**2
     &      -3D0*CW/SW*Zcc11**2)
      AAABR(211) = EE*(SW/CW*Zuu11**2+4D0*SW/CW*Zuu12**2
     &      -3D0*CW/SW*Zuu11**2)
      AAABR(212) = EE*(SW/CW*Zt11**2+4D0*SW/CW*Zt12**2-3D0
     &      *CW/SW*Zt11**2)
      AAABR(213) = EE*(SW/CW*Zt11*Zt21+4D0*SW/CW*Zt12*Zt22
     &      -3D0*CW/SW*Zt11*Zt21)
      AAABR(214) = EE*(SW/CW*Zcc11*Zcc21+4D0*SW/CW*Zcc12
     &      *Zcc22-3D0*CW/SW*Zcc11*Zcc21)
      AAABR(215) = EE*(SW/CW*Zuu11*Zuu21+4D0*SW/CW*Zuu12
     &      *Zuu22-3D0*CW/SW*Zuu11*Zuu21)
      AAABR(216) = EE*(SW/CW*Zt21**2+4D0*SW/CW*Zt22**2-3D0
     &      *CW/SW*Zt21**2)
      AAABR(217) = EE*(SW/CW*Zcc21**2+4D0*SW/CW*Zcc22**2
     &      -3D0*CW/SW*Zcc21**2)
      AAABR(218) = EE*(SW/CW*Zuu21**2+4D0*SW/CW*Zuu22**2
     &      -3D0*CW/SW*Zuu21**2)
      AAABR(219) = EE/SW
      AAABR(220) = EE/SW*(Zh11*cb+Zh12*sb)
      AAABR(221) = EE/SW*(Zh21*cb+Zh22*sb)
      AAABR(222) = EE/SW*(Zh31*cb+Zh32*sb)
      AAABR(223) = EE*(SW/CW*Zh21*cb+SW/CW*Zh22*sb+CW/SW
     &      *Zh21*cb+CW/SW*Zh22*sb)
      AAABR(224) = EE*(SW/CW*Zh31*cb+SW/CW*Zh32*sb+CW/SW
     &      *Zh31*cb+CW/SW*Zh32*sb)
      AAABR(225) = EE*(SW/CW*Zh11*cb+SW/CW*Zh12*sb+CW/SW
     &      *Zh11*cb+CW/SW*Zh12*sb)
      AAABR(226) = EE*Pa11*(SW/CW*Zh11*sb-SW/CW*Zh12*cb+CW
     &      /SW*Zh11*sb-CW/SW*Zh12*cb)
      AAABR(227) = EE*Pa11*(SW/CW*Zh31*sb-SW/CW*Zh32*cb+CW
     &      /SW*Zh31*sb-CW/SW*Zh32*cb)
      AAABR(228) = EE*Pa21*(SW/CW*Zh21*sb-SW/CW*Zh22*cb+CW
     &      /SW*Zh21*sb-CW/SW*Zh22*cb)
      AAABR(229) = EE*Pa21*(SW/CW*Zh31*sb-SW/CW*Zh32*cb+CW
     &      /SW*Zh31*sb-CW/SW*Zh32*cb)
      AAABR(230) = EE*Pa21*(SW/CW*Zh11*sb-SW/CW*Zh12*cb+CW
     &      /SW*Zh11*sb-CW/SW*Zh12*cb)
      AAABR(231) = EE*Pa11*(SW/CW*Zh21*sb-SW/CW*Zh22*cb+CW
     &      /SW*Zh21*sb-CW/SW*Zh22*cb)
      AAABR(232) = EE/CW*MW*SW
      AAABR(233) = EE*MW*(2D0*SW*Zh11*cb+2D0*SW*Zh12*sb+SW
     &      **3/CW**2*Zh11*cb+SW**3/CW**2*Zh12*sb+CW**2/SW
     &      *Zh11*cb+CW**2/SW*Zh12*sb)
      AAABR(234) = EE*MW*(2D0*SW*Zh21*cb+2D0*SW*Zh22*sb+SW
     &      **3/CW**2*Zh21*cb+SW**3/CW**2*Zh22*sb+CW**2/SW
     &      *Zh21*cb+CW**2/SW*Zh22*sb)
      AAABR(235) = EE*MW*(2D0*SW*Zh31*cb+2D0*SW*Zh32*sb+SW
     &      **3/CW**2*Zh31*cb+SW**3/CW**2*Zh32*sb+CW**2/SW
     &      *Zh31*cb+CW**2/SW*Zh32*sb)
      AAABR(236) = EE/MW/SW*Sqrt2*(NMM23*Zu22*Zv21+NMM24
     &      *Zu21*Zv22)
      AAABR(237) = EE/MW/SW*Sqrt2*(NMM23*Zu12*Zv21+NMM24
     &      *Zu11*Zv22)
      AAABR(238) = EE/MW/SW*Sqrt2*(NMM23*Zu22*Zv11+NMM24
     &      *Zu21*Zv12)
      AAABR(239) = EE/MW/SW*Sqrt2*(NMM23*Zu12*Zv11+NMM24
     &      *Zu11*Zv12)
      AAABR(240) = Sqrt2*(NMM34*Pa22*Zu12*Zv12/xvev+EE/MW
     &      *NMM23*Pa21/SW*Zu12*Zv11/cb*sb-EE/MW*NMM24*Pa21
     &      /SW*Zu11*Zv12/cb/sb+EE/MW*NMM24*Pa21/SW*Zu11
     &      *Zv12/cb*sb)
      AAABR(241) = Sqrt2*(NMM34*Pa22*Zu22*Zv12/xvev+EE/MW
     &      *NMM23*Pa21/SW*Zu22*Zv11/cb*sb-EE/MW*NMM24*Pa21
     &      /SW*Zu21*Zv12/cb/sb+EE/MW*NMM24*Pa21/SW*Zu21
     &      *Zv12/cb*sb)
      AAABR(242) = Sqrt2*(NMM34*Pa22*Zu12*Zv22/xvev+EE/MW
     &      *NMM23*Pa21/SW*Zu12*Zv21/cb*sb-EE/MW*NMM24*Pa21
     &      /SW*Zu11*Zv22/cb/sb+EE/MW*NMM24*Pa21/SW*Zu11
     &      *Zv22/cb*sb)
      AAABR(243) = Sqrt2*(NMM34*Pa12*Zu22*Zv22/xvev+EE/MW
     &      *NMM23*Pa11/SW*Zu22*Zv21/cb*sb-EE/MW*NMM24*Pa11
     &      /SW*Zu21*Zv22/cb/sb+EE/MW*NMM24*Pa11/SW*Zu21
     &      *Zv22/cb*sb)
      AAABR(244) = Sqrt2*(NMM34*Pa22*Zu22*Zv22/xvev+EE/MW
     &      *NMM23*Pa21/SW*Zu22*Zv21/cb*sb-EE/MW*NMM24*Pa21
     &      /SW*Zu21*Zv22/cb/sb+EE/MW*NMM24*Pa21/SW*Zu21
     &      *Zv22/cb*sb)
      AAABR(245) = Sqrt2*(NMM34*Pa12*Zu12*Zv22/xvev+EE/MW
     &      *NMM23*Pa11/SW*Zu12*Zv21/cb*sb-EE/MW*NMM24*Pa11
     &      /SW*Zu11*Zv22/cb/sb+EE/MW*NMM24*Pa11/SW*Zu11
     &      *Zv22/cb*sb)
      AAABR(246) = Sqrt2*(NMM34*Pa12*Zu22*Zv12/xvev+EE/MW
     &      *NMM23*Pa11/SW*Zu22*Zv11/cb*sb-EE/MW*NMM24*Pa11
     &      /SW*Zu21*Zv12/cb/sb+EE/MW*NMM24*Pa11/SW*Zu21
     &      *Zv12/cb*sb)
      AAABR(247) = Sqrt2*(NMM34*Pa12*Zu12*Zv12/xvev+EE/MW
     &      *NMM23*Pa11/SW*Zu12*Zv11/cb*sb-EE/MW*NMM24*Pa11
     &      /SW*Zu11*Zv12/cb/sb+EE/MW*NMM24*Pa11/SW*Zu11
     &      *Zv12/cb*sb)
      AAABR(248) = Sqrt2*(NMM34*Zh13*Zu22*Zv22/xvev-EE/MW
     &      *NMM23/SW*Zh11*Zu22*Zv21/cb+EE/MW*NMM24/SW*Zh12
     &      *Zu21*Zv22/sb)
      AAABR(249) = Sqrt2*(NMM34*Zh23*Zu22*Zv22/xvev-EE/MW
     &      *NMM23/SW*Zh21*Zu22*Zv21/cb+EE/MW*NMM24/SW*Zh22
     &      *Zu21*Zv22/sb)
      AAABR(250) = Sqrt2*(NMM34*Zh13*Zu12*Zv12/xvev-EE/MW
     &      *NMM23/SW*Zh11*Zu12*Zv11/cb+EE/MW*NMM24/SW*Zh12
     &      *Zu11*Zv12/sb)
      AAABR(251) = Sqrt2*(NMM34*Zh23*Zu22*Zv12/xvev-EE/MW
     &      *NMM23/SW*Zh21*Zu22*Zv11/cb+EE/MW*NMM24/SW*Zh22
     &      *Zu21*Zv12/sb)
      AAABR(252) = Sqrt2*(NMM34*Zh23*Zu12*Zv22/xvev-EE/MW
     &      *NMM23/SW*Zh21*Zu12*Zv21/cb+EE/MW*NMM24/SW*Zh22
     &      *Zu11*Zv22/sb)
      AAABR(253) = Sqrt2*(NMM34*Zh33*Zu22*Zv22/xvev-EE/MW
     &      *NMM23/SW*Zh31*Zu22*Zv21/cb+EE/MW*NMM24/SW*Zh32
     &      *Zu21*Zv22/sb)
      AAABR(254) = Sqrt2*(NMM34*Zh23*Zu12*Zv12/xvev-EE/MW
     &      *NMM23/SW*Zh21*Zu12*Zv11/cb+EE/MW*NMM24/SW*Zh22
     &      *Zu11*Zv12/sb)
      AAABR(255) = Sqrt2*(NMM34*Zh13*Zu22*Zv12/xvev-EE/MW
     &      *NMM23/SW*Zh11*Zu22*Zv11/cb+EE/MW*NMM24/SW*Zh12
     &      *Zu21*Zv12/sb)
      AAABR(256) = Sqrt2*(NMM34*Zh13*Zu12*Zv22/xvev-EE/MW
     &      *NMM23/SW*Zh11*Zu12*Zv21/cb+EE/MW*NMM24/SW*Zh12
     &      *Zu11*Zv22/sb)
      AAABR(257) = Sqrt2*(NMM34*Zh33*Zu22*Zv12/xvev-EE/MW
     &      *NMM23/SW*Zh31*Zu22*Zv11/cb+EE/MW*NMM24/SW*Zh32
     &      *Zu21*Zv12/sb)
      AAABR(258) = Sqrt2*(NMM34*Zh33*Zu12*Zv22/xvev-EE/MW
     &      *NMM23/SW*Zh31*Zu12*Zv21/cb+EE/MW*NMM24/SW*Zh32
     &      *Zu11*Zv22/sb)
      AAABR(259) = Sqrt2*(NMM34*Zh33*Zu12*Zv12/xvev-EE/MW
     &      *NMM23/SW*Zh31*Zu12*Zv11/cb+EE/MW*NMM24/SW*Zh32
     &      *Zu11*Zv12/sb)
      AAABR(260) = EE/SW*Zu11
      AAABR(261) = EE/SW*Zu21
      AAABR(262) = EE/SW*(2D0*Zl11*Zu11-Ml/MW*Sqrt2*Zl12
     &      *Zu12/cb)
      AAABR(263) = EE/SW*(2D0*Zl11*Zu21-Ml/MW*Sqrt2*Zl12
     &      *Zu22/cb)
      AAABR(264) = EE/SW*(2D0*Zl21*Zu11-Ml/MW*Sqrt2*Zl22
     &      *Zu12/cb)
      AAABR(265) = EE/SW*(2D0*Zl21*Zu21-Ml/MW*Sqrt2*Zl22
     &      *Zu22/cb)
      AAABR(266) = EE/MW/SW/cb*(2D0*NMM24*Zn24*Zv11/sb-2D0
     &      *NMM24*Zn24*Zv11*sb-NMM14*Sqrt2*Zn21*Zv12/sb
     &      +NMM14*Sqrt2*Zn21*Zv12*sb+NMM24*Sqrt2*Zn22*Zv12
     &      /sb-NMM24*Sqrt2*Zn22*Zv12*sb+NMM45*Sqrt2*Zn25
     &      *Zv12*sb)
      AAABR(267) = EE/MW/SW/cb*(2D0*NMM23*Zn23*Zu11*sb+NMM13
     &      *Sqrt2*Zn21*Zu12*sb-NMM23*Sqrt2*Zn22*Zu12*sb
     &      -NMM35*Sqrt2*Zn25*Zu12/sb+NMM35*Sqrt2*Zn25*Zu12
     &      *sb)
      AAABR(268) = EE/MW/SW/cb*(2D0*NMM24*Zn54*Zv11/sb-2D0
     &      *NMM24*Zn54*Zv11*sb-NMM14*Sqrt2*Zn51*Zv12/sb
     &      +NMM14*Sqrt2*Zn51*Zv12*sb+NMM24*Sqrt2*Zn52*Zv12
     &      /sb-NMM24*Sqrt2*Zn52*Zv12*sb+NMM45*Sqrt2*Zn55
     &      *Zv12*sb)
      AAABR(269) = EE/MW/SW/cb*(2D0*NMM23*Zn53*Zu11*sb+NMM13
     &      *Sqrt2*Zn51*Zu12*sb-NMM23*Sqrt2*Zn52*Zu12*sb
     &      -NMM35*Sqrt2*Zn55*Zu12/sb+NMM35*Sqrt2*Zn55*Zu12
     &      *sb)
      AAABR(270) = EE/MW/SW/cb*(2D0*NMM24*Zn44*Zv11/sb-2D0
     &      *NMM24*Zn44*Zv11*sb-NMM14*Sqrt2*Zn41*Zv12/sb
     &      +NMM14*Sqrt2*Zn41*Zv12*sb+NMM24*Sqrt2*Zn42*Zv12
     &      /sb-NMM24*Sqrt2*Zn42*Zv12*sb+NMM45*Sqrt2*Zn45
     &      *Zv12*sb)
      AAABR(271) = EE/MW/SW/cb*(2D0*NMM23*Zn43*Zu11*sb+NMM13
     &      *Sqrt2*Zn41*Zu12*sb-NMM23*Sqrt2*Zn42*Zu12*sb
     &      -NMM35*Sqrt2*Zn45*Zu12/sb+NMM35*Sqrt2*Zn45*Zu12
     &      *sb)
      AAABR(272) = EE/MW/SW/cb*(2D0*NMM24*Zn34*Zv11/sb-2D0
     &      *NMM24*Zn34*Zv11*sb-NMM14*Sqrt2*Zn31*Zv12/sb
     &      +NMM14*Sqrt2*Zn31*Zv12*sb+NMM24*Sqrt2*Zn32*Zv12
     &      /sb-NMM24*Sqrt2*Zn32*Zv12*sb+NMM45*Sqrt2*Zn35
     &      *Zv12*sb)
      AAABR(273) = EE/MW/SW/cb*(2D0*NMM23*Zn33*Zu11*sb+NMM13
     &      *Sqrt2*Zn31*Zu12*sb-NMM23*Sqrt2*Zn32*Zu12*sb
     &      -NMM35*Sqrt2*Zn35*Zu12/sb+NMM35*Sqrt2*Zn35*Zu12
     &      *sb)
      AAABR(274) = EE/MW/SW/cb*(2D0*NMM24*Zn14*Zv11/sb-2D0
     &      *NMM24*Zn14*Zv11*sb-NMM14*Sqrt2*Zn11*Zv12/sb
     &      +NMM14*Sqrt2*Zn11*Zv12*sb+NMM24*Sqrt2*Zn12*Zv12
     &      /sb-NMM24*Sqrt2*Zn12*Zv12*sb+NMM45*Sqrt2*Zn15
     &      *Zv12*sb)
      AAABR(275) = EE/MW/SW/cb*(2D0*NMM23*Zn13*Zu11*sb+NMM13
     &      *Sqrt2*Zn11*Zu12*sb-NMM23*Sqrt2*Zn12*Zu12*sb
     &      -NMM35*Sqrt2*Zn15*Zu12/sb+NMM35*Sqrt2*Zn15*Zu12
     &      *sb)
      AAABR(276) = EE/MW/SW/cb*(2D0*NMM24*Zn34*Zv21/sb-2D0
     &      *NMM24*Zn34*Zv21*sb-NMM14*Sqrt2*Zn31*Zv22/sb
     &      +NMM14*Sqrt2*Zn31*Zv22*sb+NMM24*Sqrt2*Zn32*Zv22
     &      /sb-NMM24*Sqrt2*Zn32*Zv22*sb+NMM45*Sqrt2*Zn35
     &      *Zv22*sb)
      AAABR(277) = EE/MW/SW/cb*(2D0*NMM23*Zn33*Zu21*sb+NMM13
     &      *Sqrt2*Zn31*Zu22*sb-NMM23*Sqrt2*Zn32*Zu22*sb
     &      -NMM35*Sqrt2*Zn35*Zu22/sb+NMM35*Sqrt2*Zn35*Zu22
     &      *sb)
      AAABR(278) = EE/MW/SW/cb*(2D0*NMM24*Zn24*Zv21/sb-2D0
     &      *NMM24*Zn24*Zv21*sb-NMM14*Sqrt2*Zn21*Zv22/sb
     &      +NMM14*Sqrt2*Zn21*Zv22*sb+NMM24*Sqrt2*Zn22*Zv22
     &      /sb-NMM24*Sqrt2*Zn22*Zv22*sb+NMM45*Sqrt2*Zn25
     &      *Zv22*sb)
      AAABR(279) = EE/MW/SW/cb*(2D0*NMM23*Zn23*Zu21*sb+NMM13
     &      *Sqrt2*Zn21*Zu22*sb-NMM23*Sqrt2*Zn22*Zu22*sb
     &      -NMM35*Sqrt2*Zn25*Zu22/sb+NMM35*Sqrt2*Zn25*Zu22
     &      *sb)
      AAABR(280) = EE/MW/SW/cb*(2D0*NMM24*Zn14*Zv21/sb-2D0
     &      *NMM24*Zn14*Zv21*sb-NMM14*Sqrt2*Zn11*Zv22/sb
     &      +NMM14*Sqrt2*Zn11*Zv22*sb+NMM24*Sqrt2*Zn12*Zv22
     &      /sb-NMM24*Sqrt2*Zn12*Zv22*sb+NMM45*Sqrt2*Zn15
     &      *Zv22*sb)
      AAABR(281) = EE/MW/SW/cb*(2D0*NMM23*Zn13*Zu21*sb+NMM13
     &      *Sqrt2*Zn11*Zu22*sb-NMM23*Sqrt2*Zn12*Zu22*sb
     &      -NMM35*Sqrt2*Zn15*Zu22/sb+NMM35*Sqrt2*Zn15*Zu22
     &      *sb)
      AAABR(282) = EE/MW/SW/cb*(2D0*NMM24*Zn44*Zv21/sb-2D0
     &      *NMM24*Zn44*Zv21*sb-NMM14*Sqrt2*Zn41*Zv22/sb
     &      +NMM14*Sqrt2*Zn41*Zv22*sb+NMM24*Sqrt2*Zn42*Zv22
     &      /sb-NMM24*Sqrt2*Zn42*Zv22*sb+NMM45*Sqrt2*Zn45
     &      *Zv22*sb)
      AAABR(283) = EE/MW/SW/cb*(2D0*NMM23*Zn43*Zu21*sb+NMM13
     &      *Sqrt2*Zn41*Zu22*sb-NMM23*Sqrt2*Zn42*Zu22*sb
     &      -NMM35*Sqrt2*Zn45*Zu22/sb+NMM35*Sqrt2*Zn45*Zu22
     &      *sb)
      AAABR(284) = EE/MW/SW/cb*(2D0*NMM24*Zn54*Zv21/sb-2D0
     &      *NMM24*Zn54*Zv21*sb-NMM14*Sqrt2*Zn51*Zv22/sb
     &      +NMM14*Sqrt2*Zn51*Zv22*sb+NMM24*Sqrt2*Zn52*Zv22
     &      /sb-NMM24*Sqrt2*Zn52*Zv22*sb+NMM45*Sqrt2*Zn55
     &      *Zv22*sb)
      AAABR(285) = EE/MW/SW/cb*(2D0*NMM23*Zn53*Zu21*sb+NMM13
     &      *Sqrt2*Zn51*Zu22*sb-NMM23*Sqrt2*Zn52*Zu22*sb
     &      -NMM35*Sqrt2*Zn55*Zu22/sb+NMM35*Sqrt2*Zn55*Zu22
     &      *sb)
      AAABR(286) = EE/MW/SW*(2D0*NMM24*Zn54*Zv21-NMM14*Sqrt2
     &      *Zn51*Zv22+NMM24*Sqrt2*Zn52*Zv22-NMM45*Sqrt2
     &      *Zn55*Zv22)
      AAABR(287) = EE/MW/SW*(2D0*NMM23*Zn53*Zu21+NMM13*Sqrt2
     &      *Zn51*Zu22-NMM23*Sqrt2*Zn52*Zu22+NMM35*Sqrt2
     &      *Zn55*Zu22)
      AAABR(288) = EE/MW/SW*(2D0*NMM24*Zn34*Zv21-NMM14*Sqrt2
     &      *Zn31*Zv22+NMM24*Sqrt2*Zn32*Zv22-NMM45*Sqrt2
     &      *Zn35*Zv22)
      AAABR(289) = EE/MW/SW*(2D0*NMM23*Zn33*Zu21+NMM13*Sqrt2
     &      *Zn31*Zu22-NMM23*Sqrt2*Zn32*Zu22+NMM35*Sqrt2
     &      *Zn35*Zu22)
      AAABR(290) = EE/MW/SW*(2D0*NMM24*Zn24*Zv21-NMM14*Sqrt2
     &      *Zn21*Zv22+NMM24*Sqrt2*Zn22*Zv22-NMM45*Sqrt2
     &      *Zn25*Zv22)
      AAABR(291) = EE/MW/SW*(2D0*NMM23*Zn23*Zu21+NMM13*Sqrt2
     &      *Zn21*Zu22-NMM23*Sqrt2*Zn22*Zu22+NMM35*Sqrt2
     &      *Zn25*Zu22)
      AAABR(292) = EE/MW/SW*(2D0*NMM24*Zn44*Zv11-NMM14*Sqrt2
     &      *Zn41*Zv12+NMM24*Sqrt2*Zn42*Zv12-NMM45*Sqrt2
     &      *Zn45*Zv12)
      AAABR(293) = EE/MW/SW*(2D0*NMM23*Zn43*Zu11+NMM13*Sqrt2
     &      *Zn41*Zu12-NMM23*Sqrt2*Zn42*Zu12+NMM35*Sqrt2
     &      *Zn45*Zu12)
      AAABR(294) = EE/MW/SW*(2D0*NMM24*Zn34*Zv11-NMM14*Sqrt2
     &      *Zn31*Zv12+NMM24*Sqrt2*Zn32*Zv12-NMM45*Sqrt2
     &      *Zn35*Zv12)
      AAABR(295) = EE/MW/SW*(2D0*NMM23*Zn33*Zu11+NMM13*Sqrt2
     &      *Zn31*Zu12-NMM23*Sqrt2*Zn32*Zu12+NMM35*Sqrt2
     &      *Zn35*Zu12)
      AAABR(296) = EE/MW/SW*(2D0*NMM24*Zn14*Zv21-NMM14*Sqrt2
     &      *Zn11*Zv22+NMM24*Sqrt2*Zn12*Zv22-NMM45*Sqrt2
     &      *Zn15*Zv22)
      AAABR(297) = EE/MW/SW*(2D0*NMM23*Zn13*Zu21+NMM13*Sqrt2
     &      *Zn11*Zu22-NMM23*Sqrt2*Zn12*Zu22+NMM35*Sqrt2
     &      *Zn15*Zu22)
      AAABR(298) = EE/MW/SW*(2D0*NMM24*Zn14*Zv11-NMM14*Sqrt2
     &      *Zn11*Zv12+NMM24*Sqrt2*Zn12*Zv12-NMM45*Sqrt2
     &      *Zn15*Zv12)
      AAABR(299) = EE/MW/SW*(2D0*NMM23*Zn13*Zu11+NMM13*Sqrt2
     &      *Zn11*Zu12-NMM23*Sqrt2*Zn12*Zu12+NMM35*Sqrt2
     &      *Zn15*Zu12)
      AAABR(300) = EE/MW/SW*(2D0*NMM24*Zn54*Zv11-NMM14*Sqrt2
     &      *Zn51*Zv12+NMM24*Sqrt2*Zn52*Zv12-NMM45*Sqrt2
     &      *Zn55*Zv12)
      AAABR(301) = EE/MW/SW*(2D0*NMM23*Zn53*Zu11+NMM13*Sqrt2
     &      *Zn51*Zu12-NMM23*Sqrt2*Zn52*Zu12+NMM35*Sqrt2
     &      *Zn55*Zu12)
      AAABR(302) = EE/MW/SW*(2D0*NMM24*Zn24*Zv11-NMM14*Sqrt2
     &      *Zn21*Zv12+NMM24*Sqrt2*Zn22*Zv12-NMM45*Sqrt2
     &      *Zn25*Zv12)
      AAABR(303) = EE/MW/SW*(2D0*NMM23*Zn23*Zu11+NMM13*Sqrt2
     &      *Zn21*Zu12-NMM23*Sqrt2*Zn22*Zu12+NMM35*Sqrt2
     &      *Zn25*Zu12)
      AAABR(304) = EE/MW/SW*(2D0*NMM24*Zn44*Zv21-NMM14*Sqrt2
     &      *Zn41*Zv22+NMM24*Sqrt2*Zn42*Zv22-NMM45*Sqrt2
     &      *Zn45*Zv22)
      AAABR(305) = EE/MW/SW*(2D0*NMM23*Zn43*Zu21+NMM13*Sqrt2
     &      *Zn41*Zu22-NMM23*Sqrt2*Zn42*Zu22+NMM35*Sqrt2
     &      *Zn45*Zu22)
      AAABR(306) = EE/SW*(2D0*Zb11*Zu11-MbMM/MW*Sqrt2*Zb12
     &      *Zu12/cb)
      AAABR(307) = EE/MW*MtMM/SW*Sqrt2*Zb11*Zv12/sb
      AAABR(308) = EE/SW*(2D0*Zss11*Zu11-Mq/MW*Sqrt2*Zss12
     &      *Zu12/cb)
      AAABR(309) = EE/MW*Mc/SW*Sqrt2*Zss11*Zv12/sb
      AAABR(310) = EE/SW*(2D0*Zss11*Zu21-Mq/MW*Sqrt2*Zss12
     &      *Zu22/cb)
      AAABR(311) = EE/MW*Mc/SW*Sqrt2*Zss11*Zv22/sb
      AAABR(312) = EE/SW*(2D0*Zb11*Zu21-MbMM/MW*Sqrt2*Zb12
     &      *Zu22/cb)
      AAABR(313) = EE/MW*MtMM/SW*Sqrt2*Zb11*Zv22/sb
      AAABR(314) = EE/SW*(2D0*Zdd11*Zu21-Mq/MW*Sqrt2*Zdd12
     &      *Zu22/cb)
      AAABR(315) = EE/MW*Mq/SW*Sqrt2*Zdd11*Zv22/sb
      AAABR(316) = EE/SW*(2D0*Zdd11*Zu11-Mq/MW*Sqrt2*Zdd12
     &      *Zu12/cb)
      AAABR(317) = EE/MW*Mq/SW*Sqrt2*Zdd11*Zv12/sb
      AAABR(318) = EE/SW*(2D0*Zb21*Zu21-MbMM/MW*Sqrt2*Zb22
     &      *Zu22/cb)
      AAABR(319) = EE/MW*MtMM/SW*Sqrt2*Zb21*Zv22/sb
      AAABR(320) = EE/SW*(2D0*Zss21*Zu11-Mq/MW*Sqrt2*Zss22
     &      *Zu12/cb)
      AAABR(321) = EE/MW*Mc/SW*Sqrt2*Zss21*Zv12/sb
      AAABR(322) = EE/SW*(2D0*Zb21*Zu11-MbMM/MW*Sqrt2*Zb22
     &      *Zu12/cb)
      AAABR(323) = EE/MW*MtMM/SW*Sqrt2*Zb21*Zv12/sb
      AAABR(324) = EE/SW*(2D0*Zdd21*Zu21-Mq/MW*Sqrt2*Zdd22
     &      *Zu22/cb)
      AAABR(325) = EE/MW*Mq/SW*Sqrt2*Zdd21*Zv22/sb
      AAABR(326) = EE/SW*(2D0*Zdd21*Zu11-Mq/MW*Sqrt2*Zdd22
     &      *Zu12/cb)
      AAABR(327) = EE/MW*Mq/SW*Sqrt2*Zdd21*Zv12/sb
      AAABR(328) = EE/SW*(2D0*Zss21*Zu21-Mq/MW*Sqrt2*Zss22
     &      *Zu22/cb)
      AAABR(329) = EE/MW*Mc/SW*Sqrt2*Zss21*Zv22/sb
      AAABR(330) = EE/MW*Ml/SW
      AAABR(331) = EE/MW*Ml*Pa11/SW*tB
      AAABR(332) = EE/MW*Ml*Pa21/SW*tB
      AAABR(333) = EE/MW*Ml/SW*Zh31/cb
      AAABR(334) = EE/MW*Ml/SW*Zh11/cb
      AAABR(335) = EE/MW*Ml/SW*Zh21/cb
      AAABR(336) = EE/MW*Ml/SW*Sqrt2*tB
      AAABR(337) = EE/MW*Ml/SW*Sqrt2
      AAABR(338) = EE*Sqrt2*(Zn21/CW+Zn22/SW)
      AAABR(339) = EE*Sqrt2*(Zn11/CW+Zn12/SW)
      AAABR(340) = EE*Sqrt2*(Zn41/CW+Zn42/SW)
      AAABR(341) = EE*Sqrt2*(Zl11/CW*Zn41+Zl11/SW*Zn42-Ml
     &      /MW/SW*Zl12*Zn43/cb)
      AAABR(342) = EE*Sqrt2*(2D0*Zl12/CW*Zn41+Ml/MW/SW*Zl11
     &      *Zn43/cb)
      AAABR(343) = EE*Sqrt2*(Zl11/CW*Zn11+Zl11/SW*Zn12-Ml
     &      /MW/SW*Zl12*Zn13/cb)
      AAABR(344) = EE*Sqrt2*(2D0*Zl12/CW*Zn11+Ml/MW/SW*Zl11
     &      *Zn13/cb)
      AAABR(345) = EE*Sqrt2*(Zn51/CW+Zn52/SW)
      AAABR(346) = EE*Sqrt2*(Zn31/CW+Zn32/SW)
      AAABR(347) = EE*Sqrt2*(Zl11/CW*Zn21+Zl11/SW*Zn22-Ml
     &      /MW/SW*Zl12*Zn23/cb)
      AAABR(348) = EE*Sqrt2*(2D0*Zl12/CW*Zn21+Ml/MW/SW*Zl11
     &      *Zn23/cb)
      AAABR(349) = EE*Sqrt2*(Zl11/CW*Zn31+Zl11/SW*Zn32-Ml
     &      /MW/SW*Zl12*Zn33/cb)
      AAABR(350) = EE*Sqrt2*(2D0*Zl12/CW*Zn31+Ml/MW/SW*Zl11
     &      *Zn33/cb)
      AAABR(351) = EE*Sqrt2*(Zl11/CW*Zn51+Zl11/SW*Zn52-Ml
     &      /MW/SW*Zl12*Zn53/cb)
      AAABR(352) = EE*Sqrt2*(2D0*Zl12/CW*Zn51+Ml/MW/SW*Zl11
     &      *Zn53/cb)
      AAABR(353) = EE/CW*Sqrt2*Zn11
      AAABR(354) = EE*Sqrt2*(Zl21/CW*Zn41+Zl21/SW*Zn42-Ml
     &      /MW/SW*Zl22*Zn43/cb)
      AAABR(355) = EE*Sqrt2*(2D0*Zl22/CW*Zn41+Ml/MW/SW*Zl21
     &      *Zn43/cb)
      AAABR(356) = EE/CW*Sqrt2*Zn31
      AAABR(357) = EE/CW*Sqrt2*Zn41
      AAABR(358) = EE*Sqrt2*(Zl21/CW*Zn31+Zl21/SW*Zn32-Ml
     &      /MW/SW*Zl22*Zn33/cb)
      AAABR(359) = EE*Sqrt2*(2D0*Zl22/CW*Zn31+Ml/MW/SW*Zl21
     &      *Zn33/cb)
      AAABR(360) = EE/CW*Sqrt2*Zn51
      AAABR(361) = EE*Sqrt2*(Zl21/CW*Zn21+Zl21/SW*Zn22-Ml
     &      /MW/SW*Zl22*Zn23/cb)
      AAABR(362) = EE*Sqrt2*(2D0*Zl22/CW*Zn21+Ml/MW/SW*Zl21
     &      *Zn23/cb)
      AAABR(363) = EE/CW*Sqrt2*Zn21
      AAABR(364) = EE*Sqrt2*(Zl21/CW*Zn51+Zl21/SW*Zn52-Ml
     &      /MW/SW*Zl22*Zn53/cb)
      AAABR(365) = EE*Sqrt2*(2D0*Zl22/CW*Zn51+Ml/MW/SW*Zl21
     &      *Zn53/cb)
      AAABR(366) = EE*Sqrt2*(Zl21/CW*Zn11+Zl21/SW*Zn12-Ml
     &      /MW/SW*Zl22*Zn13/cb)
      AAABR(367) = EE*Sqrt2*(2D0*Zl22/CW*Zn11+Ml/MW/SW*Zl21
     &      *Zn13/cb)
      AAABR(368) = EE/SW*Zv11
      AAABR(369) = EE/MW*Ml/SW*Sqrt2*Zu12/cb
      AAABR(370) = EE/SW*Zv21
      AAABR(371) = EE/MW*Ml/SW*Sqrt2*Zu22/cb
      AAABR(372) = EE*Sqrt2*(Zn51/CW-Zn52/SW)
      AAABR(373) = EE*Sqrt2*(Zn41/CW-Zn42/SW)
      AAABR(374) = EE*Sqrt2*(Zn31/CW-Zn32/SW)
      AAABR(375) = EE*Sqrt2*(Zn21/CW-Zn22/SW)
      AAABR(376) = EE*Sqrt2*(Zn11/CW-Zn12/SW)
      AAABR(377) = EE*Sqrt2*(Zn21/CW*Zss11-3D0*Zn22/SW*Zss11
     &      +3D0*Mq/MW/SW*Zn23*Zss12/cb)
      AAABR(378) = EE*Sqrt2*(2D0*Zn21/CW*Zss12+3D0*Mq/MW
     &      /SW*Zn23*Zss11/cb)
      AAABR(379) = EE*Sqrt2*(Zdd11/CW*Zn41-3D0*Zdd11/SW*Zn42
     &      +3D0*Mq/MW/SW*Zdd12*Zn43/cb)
      AAABR(380) = EE*Sqrt2*(2D0*Zdd12/CW*Zn41+3D0*Mq/MW
     &      /SW*Zdd11*Zn43/cb)
      AAABR(381) = EE*Sqrt2*(Zb11/CW*Zn11-3D0*Zb11/SW*Zn12
     &      +3D0*MbMM/MW/SW*Zb12*Zn13/cb)
      AAABR(382) = EE*Sqrt2*(2D0*Zb12/CW*Zn11+3D0*MbMM/MW
     &      /SW*Zb11*Zn13/cb)
      AAABR(383) = EE*Sqrt2*(Zb11/CW*Zn41-3D0*Zb11/SW*Zn42
     &      +3D0*MbMM/MW/SW*Zb12*Zn43/cb)
      AAABR(384) = EE*Sqrt2*(2D0*Zb12/CW*Zn41+3D0*MbMM/MW
     &      /SW*Zb11*Zn43/cb)
      AAABR(385) = EE*Sqrt2*(Zn11/CW*Zss11-3D0*Zn12/SW*Zss11
     &      +3D0*Mq/MW/SW*Zn13*Zss12/cb)
      AAABR(386) = EE*Sqrt2*(2D0*Zn11/CW*Zss12+3D0*Mq/MW
     &      /SW*Zn13*Zss11/cb)
      AAABR(387) = EE*Sqrt2*(Zdd11/CW*Zn21-3D0*Zdd11/SW*Zn22
     &      +3D0*Mq/MW/SW*Zdd12*Zn23/cb)
      AAABR(388) = EE*Sqrt2*(2D0*Zdd12/CW*Zn21+3D0*Mq/MW
     &      /SW*Zdd11*Zn23/cb)
      AAABR(389) = EE*Sqrt2*(Zb11/CW*Zn31-3D0*Zb11/SW*Zn32
     &      +3D0*MbMM/MW/SW*Zb12*Zn33/cb)
      AAABR(390) = EE*Sqrt2*(2D0*Zb12/CW*Zn31+3D0*MbMM/MW
     &      /SW*Zb11*Zn33/cb)
      AAABR(391) = EE*Sqrt2*(Zdd11/CW*Zn51-3D0*Zdd11/SW*Zn52
     &      +3D0*Mq/MW/SW*Zdd12*Zn53/cb)
      AAABR(392) = EE*Sqrt2*(2D0*Zdd12/CW*Zn51+3D0*Mq/MW
     &      /SW*Zdd11*Zn53/cb)
      AAABR(393) = EE*Sqrt2*(Zn41/CW*Zss11-3D0*Zn42/SW*Zss11
     &      +3D0*Mq/MW/SW*Zn43*Zss12/cb)
      AAABR(394) = EE*Sqrt2*(2D0*Zn41/CW*Zss12+3D0*Mq/MW
     &      /SW*Zn43*Zss11/cb)
      AAABR(395) = EE*Sqrt2*(Zn31/CW*Zss11-3D0*Zn32/SW*Zss11
     &      +3D0*Mq/MW/SW*Zn33*Zss12/cb)
      AAABR(396) = EE*Sqrt2*(2D0*Zn31/CW*Zss12+3D0*Mq/MW
     &      /SW*Zn33*Zss11/cb)
      AAABR(397) = EE*Sqrt2*(Zdd11/CW*Zn11-3D0*Zdd11/SW*Zn12
     &      +3D0*Mq/MW/SW*Zdd12*Zn13/cb)
      AAABR(398) = EE*Sqrt2*(2D0*Zdd12/CW*Zn11+3D0*Mq/MW
     &      /SW*Zdd11*Zn13/cb)
      AAABR(399) = EE*Sqrt2*(Zdd11/CW*Zn31-3D0*Zdd11/SW*Zn32
     &      +3D0*Mq/MW/SW*Zdd12*Zn33/cb)
      AAABR(400) = EE*Sqrt2*(2D0*Zdd12/CW*Zn31+3D0*Mq/MW
     &      /SW*Zdd11*Zn33/cb)
      AAABR(401) = EE*Sqrt2*(Zb11/CW*Zn21-3D0*Zb11/SW*Zn22
     &      +3D0*MbMM/MW/SW*Zb12*Zn23/cb)
      AAABR(402) = EE*Sqrt2*(2D0*Zb12/CW*Zn21+3D0*MbMM/MW
     &      /SW*Zb11*Zn23/cb)
      AAABR(403) = EE*Sqrt2*(Zb11/CW*Zn51-3D0*Zb11/SW*Zn52
     &      +3D0*MbMM/MW/SW*Zb12*Zn53/cb)
      AAABR(404) = EE*Sqrt2*(2D0*Zb12/CW*Zn51+3D0*MbMM/MW
     &      /SW*Zb11*Zn53/cb)
      AAABR(405) = EE*Sqrt2*(Zn51/CW*Zss11-3D0*Zn52/SW*Zss11
     &      +3D0*Mq/MW/SW*Zn53*Zss12/cb)
      AAABR(406) = EE*Sqrt2*(2D0*Zn51/CW*Zss12+3D0*Mq/MW
     &      /SW*Zn53*Zss11/cb)
      AAABR(407) = EE*Sqrt2*(Zb21/CW*Zn51-3D0*Zb21/SW*Zn52
     &      +3D0*MbMM/MW/SW*Zb22*Zn53/cb)
      AAABR(408) = EE*Sqrt2*(2D0*Zb22/CW*Zn51+3D0*MbMM/MW
     &      /SW*Zb21*Zn53/cb)
      AAABR(409) = EE*Sqrt2*(Zn41/CW*Zss21-3D0*Zn42/SW*Zss21
     &      +3D0*Mq/MW/SW*Zn43*Zss22/cb)
      AAABR(410) = EE*Sqrt2*(2D0*Zn41/CW*Zss22+3D0*Mq/MW
     &      /SW*Zn43*Zss21/cb)
      AAABR(411) = EE*Sqrt2*(Zb21/CW*Zn41-3D0*Zb21/SW*Zn42
     &      +3D0*MbMM/MW/SW*Zb22*Zn43/cb)
      AAABR(412) = EE*Sqrt2*(2D0*Zb22/CW*Zn41+3D0*MbMM/MW
     &      /SW*Zb21*Zn43/cb)
      AAABR(413) = EE*Sqrt2*(Zdd21/CW*Zn31-3D0*Zdd21/SW*Zn32
     &      +3D0*Mq/MW/SW*Zdd22*Zn33/cb)
      AAABR(414) = EE*Sqrt2*(2D0*Zdd22/CW*Zn31+3D0*Mq/MW
     &      /SW*Zdd21*Zn33/cb)
      AAABR(415) = EE*Sqrt2*(Zdd21/CW*Zn21-3D0*Zdd21/SW*Zn22
     &      +3D0*Mq/MW/SW*Zdd22*Zn23/cb)
      AAABR(416) = EE*Sqrt2*(2D0*Zdd22/CW*Zn21+3D0*Mq/MW
     &      /SW*Zdd21*Zn23/cb)
      AAABR(417) = EE*Sqrt2*(Zn51/CW*Zss21-3D0*Zn52/SW*Zss21
     &      +3D0*Mq/MW/SW*Zn53*Zss22/cb)
      AAABR(418) = EE*Sqrt2*(2D0*Zn51/CW*Zss22+3D0*Mq/MW
     &      /SW*Zn53*Zss21/cb)
      AAABR(419) = EE*Sqrt2*(Zb21/CW*Zn11-3D0*Zb21/SW*Zn12
     &      +3D0*MbMM/MW/SW*Zb22*Zn13/cb)
      AAABR(420) = EE*Sqrt2*(2D0*Zb22/CW*Zn11+3D0*MbMM/MW
     &      /SW*Zb21*Zn13/cb)
      AAABR(421) = EE*Sqrt2*(Zb21/CW*Zn21-3D0*Zb21/SW*Zn22
     &      +3D0*MbMM/MW/SW*Zb22*Zn23/cb)
      AAABR(422) = EE*Sqrt2*(2D0*Zb22/CW*Zn21+3D0*MbMM/MW
     &      /SW*Zb21*Zn23/cb)
      AAABR(423) = EE*Sqrt2*(Zn31/CW*Zss21-3D0*Zn32/SW*Zss21
     &      +3D0*Mq/MW/SW*Zn33*Zss22/cb)
      AAABR(424) = EE*Sqrt2*(2D0*Zn31/CW*Zss22+3D0*Mq/MW
     &      /SW*Zn33*Zss21/cb)
      AAABR(425) = EE*Sqrt2*(Zb21/CW*Zn31-3D0*Zb21/SW*Zn32
     &      +3D0*MbMM/MW/SW*Zb22*Zn33/cb)
      AAABR(426) = EE*Sqrt2*(2D0*Zb22/CW*Zn31+3D0*MbMM/MW
     &      /SW*Zb21*Zn33/cb)
      AAABR(427) = EE*Sqrt2*(Zn21/CW*Zss21-3D0*Zn22/SW*Zss21
     &      +3D0*Mq/MW/SW*Zn23*Zss22/cb)
      AAABR(428) = EE*Sqrt2*(2D0*Zn21/CW*Zss22+3D0*Mq/MW
     &      /SW*Zn23*Zss21/cb)
      AAABR(429) = EE*Sqrt2*(Zdd21/CW*Zn11-3D0*Zdd21/SW*Zn12
     &      +3D0*Mq/MW/SW*Zdd22*Zn13/cb)
      AAABR(430) = EE*Sqrt2*(2D0*Zdd22/CW*Zn11+3D0*Mq/MW
     &      /SW*Zdd21*Zn13/cb)
      AAABR(431) = EE*Sqrt2*(Zn11/CW*Zss21-3D0*Zn12/SW*Zss21
     &      +3D0*Mq/MW/SW*Zn13*Zss22/cb)
      AAABR(432) = EE*Sqrt2*(2D0*Zn11/CW*Zss22+3D0*Mq/MW
     &      /SW*Zn13*Zss21/cb)
      AAABR(433) = EE*Sqrt2*(Zdd21/CW*Zn41-3D0*Zdd21/SW*Zn42
     &      +3D0*Mq/MW/SW*Zdd22*Zn43/cb)
      AAABR(434) = EE*Sqrt2*(2D0*Zdd22/CW*Zn41+3D0*Mq/MW
     &      /SW*Zdd21*Zn43/cb)
      AAABR(435) = EE*Sqrt2*(Zdd21/CW*Zn51-3D0*Zdd21/SW*Zn52
     &      +3D0*Mq/MW/SW*Zdd22*Zn53/cb)
      AAABR(436) = EE*Sqrt2*(2D0*Zdd22/CW*Zn51+3D0*Mq/MW
     &      /SW*Zdd21*Zn53/cb)
      AAABR(437) = EE/MW*Mq/SW
      AAABR(438) = EE/MW*Mb/SW
      AAABR(439) = EE/MW*Mq*Pa21/SW/cb*(sb-Td2/sb)
      AAABR(440) = EE/MW*Mq*Pa11/SW/cb*(sb-Td2/sb)
      AAABR(441) = EE/MW*Mb*Pa21/SW/cb*(sb-Td3/sb)
      AAABR(442) = EE/MW*Mb*Pa11/SW/cb*(sb-Td3/sb)
      AAABR(443) = EE/MW*Mq/SW*(Zh21/cb-Td2*Zh21/cb+Td2*Zh22
     &      /sb)
      AAABR(444) = EE/MW*Mq/SW*(Zh11/cb-Td2*Zh11/cb+Td2*Zh12
     &      /sb)
      AAABR(445) = EE/MW*Mq/SW*(Zh31/cb-Td2*Zh31/cb+Td2*Zh32
     &      /sb)
      AAABR(446) = EE/MW*Mb/SW*(Zh11/cb-Td3*Zh11/cb+Td3*Zh12
     &      /sb)
      AAABR(447) = EE/MW*Mb/SW*(Zh31/cb-Td3*Zh31/cb+Td3*Zh32
     &      /sb)
      AAABR(448) = EE/MW*Mb/SW*(Zh21/cb-Td3*Zh21/cb+Td3*Zh22
     &      /sb)
      AAABR(449) = EE/MW*Mq/SW*Sqrt2/cb*(1D0/sb-sb)
      AAABR(450) = EE/MW*Mq/SW*Sqrt2/cb*(sb-Td2/sb)
      AAABR(451) = EE/MW*Mc/SW*Sqrt2/cb*(1D0/sb-sb)
      AAABR(452) = EE/MW*Mb/SW*Sqrt2/cb*(sb-Td3/sb)
      AAABR(453) = EE/MW*Mt/SW*Sqrt2/cb*(1D0/sb-sb)
      AAABR(454) = EE/MW*Mc/SW*Sqrt2
      AAABR(455) = EE/MW*Mq/SW*Sqrt2
      AAABR(456) = EE/MW*Mb/SW*Sqrt2
      AAABR(457) = EE/MW*Mt/SW*Sqrt2
      AAABR(458) = EE/SW*(2D0*Zt11*Zv11-MtMM/MW*Sqrt2*Zt12
     &      *Zv12/sb)
      AAABR(459) = EE/MW*MbMM/SW*Sqrt2*Zt11*Zu12/cb
      AAABR(460) = EE/SW*(2D0*Zuu11*Zv11-Mq/MW*Sqrt2*Zuu12
     &      *Zv12/sb)
      AAABR(461) = EE/MW*Mq/SW*Sqrt2*Zu12*Zuu11/cb
      AAABR(462) = EE/SW*(2D0*Zcc11*Zv11-Mc/MW*Sqrt2*Zcc12
     &      *Zv12/sb)
      AAABR(463) = EE/MW*Mq/SW*Sqrt2*Zcc11*Zu12/cb
      AAABR(464) = EE/SW*(2D0*Zcc21*Zv11-Mc/MW*Sqrt2*Zcc22
     &      *Zv12/sb)
      AAABR(465) = EE/MW*Mq/SW*Sqrt2*Zcc21*Zu12/cb
      AAABR(466) = EE/SW*(2D0*Zuu21*Zv11-Mq/MW*Sqrt2*Zuu22
     &      *Zv12/sb)
      AAABR(467) = EE/MW*Mq/SW*Sqrt2*Zu12*Zuu21/cb
      AAABR(468) = EE/SW*(2D0*Zt21*Zv11-MtMM/MW*Sqrt2*Zt22
     &      *Zv12/sb)
      AAABR(469) = EE/MW*MbMM/SW*Sqrt2*Zt21*Zu12/cb
      AAABR(470) = EE/SW*(2D0*Zt11*Zv21-MtMM/MW*Sqrt2*Zt12
     &      *Zv22/sb)
      AAABR(471) = EE/MW*MbMM/SW*Sqrt2*Zt11*Zu22/cb
      AAABR(472) = EE/SW*(2D0*Zcc11*Zv21-Mc/MW*Sqrt2*Zcc12
     &      *Zv22/sb)
      AAABR(473) = EE/MW*Mq/SW*Sqrt2*Zcc11*Zu22/cb
      AAABR(474) = EE/SW*(2D0*Zuu11*Zv21-Mq/MW*Sqrt2*Zuu12
     &      *Zv22/sb)
      AAABR(475) = EE/MW*Mq/SW*Sqrt2*Zu22*Zuu11/cb
      AAABR(476) = EE/SW*(2D0*Zt21*Zv21-MtMM/MW*Sqrt2*Zt22
     &      *Zv22/sb)
      AAABR(477) = EE/MW*MbMM/SW*Sqrt2*Zt21*Zu22/cb
      AAABR(478) = EE/SW*(2D0*Zuu21*Zv21-Mq/MW*Sqrt2*Zuu22
     &      *Zv22/sb)
      AAABR(479) = EE/MW*Mq/SW*Sqrt2*Zu22*Zuu21/cb
      AAABR(480) = EE/SW*(2D0*Zcc21*Zv21-Mc/MW*Sqrt2*Zcc22
     &      *Zv22/sb)
      AAABR(481) = EE/MW*Mq/SW*Sqrt2*Zcc21*Zu22/cb
      AAABR(482) = GG*Sqrt2*Zdd11
      AAABR(483) = GG*Sqrt2*Zdd12
      AAABR(484) = GG*Sqrt2*Zss11
      AAABR(485) = GG*Sqrt2*Zss12
      AAABR(486) = GG*Sqrt2*Zb11
      AAABR(487) = GG*Sqrt2*Zb12
      AAABR(488) = GG*Sqrt2*Zss21
      AAABR(489) = GG*Sqrt2*Zss22
      AAABR(490) = GG*Sqrt2*Zb21
      AAABR(491) = GG*Sqrt2*Zb22
      AAABR(492) = GG*Sqrt2*Zdd21
      AAABR(493) = GG*Sqrt2*Zdd22
      AAABR(494) = EE*Sqrt2*(Zn11/CW*Zt11+3D0*Zn12/SW*Zt11
     &      +3D0*MtMM/MW/SW*Zn14*Zt12/sb)
      AAABR(495) = EE*Sqrt2*(4D0*Zn11/CW*Zt12-3D0*MtMM/MW
     &      /SW*Zn14*Zt11/sb)
      AAABR(496) = EE*Sqrt2*(Zn11/CW*Zuu11+3D0*Zn12/SW*Zuu11
     &      +3D0*Mq/MW/SW*Zn14*Zuu12/sb)
      AAABR(497) = EE*Sqrt2*(4D0*Zn11/CW*Zuu12-3D0*Mq/MW
     &      /SW*Zn14*Zuu11/sb)
      AAABR(498) = EE*Sqrt2*(Zcc11/CW*Zn21+3D0*Zcc11/SW*Zn22
     &      +3D0*Mc/MW/SW*Zcc12*Zn24/sb)
      AAABR(499) = EE*Sqrt2*(4D0*Zcc12/CW*Zn21-3D0*Mc/MW
     &      /SW*Zcc11*Zn24/sb)
      AAABR(500) = EE*Sqrt2*(Zn51/CW*Zt11+3D0*Zn52/SW*Zt11
     &      +3D0*MtMM/MW/SW*Zn54*Zt12/sb)
      AAABR(501) = EE*Sqrt2*(4D0*Zn51/CW*Zt12-3D0*MtMM/MW
     &      /SW*Zn54*Zt11/sb)
      AAABR(502) = EE*Sqrt2*(Zn21/CW*Zt11+3D0*Zn22/SW*Zt11
     &      +3D0*MtMM/MW/SW*Zn24*Zt12/sb)
      AAABR(503) = EE*Sqrt2*(4D0*Zn21/CW*Zt12-3D0*MtMM/MW
     &      /SW*Zn24*Zt11/sb)
      AAABR(504) = EE*Sqrt2*(Zcc11/CW*Zn41+3D0*Zcc11/SW*Zn42
     &      +3D0*Mc/MW/SW*Zcc12*Zn44/sb)
      AAABR(505) = EE*Sqrt2*(4D0*Zcc12/CW*Zn41-3D0*Mc/MW
     &      /SW*Zcc11*Zn44/sb)
      AAABR(506) = EE*Sqrt2*(Zcc11/CW*Zn51+3D0*Zcc11/SW*Zn52
     &      +3D0*Mc/MW/SW*Zcc12*Zn54/sb)
      AAABR(507) = EE*Sqrt2*(4D0*Zcc12/CW*Zn51-3D0*Mc/MW
     &      /SW*Zcc11*Zn54/sb)
      AAABR(508) = EE*Sqrt2*(Zcc11/CW*Zn11+3D0*Zcc11/SW*Zn12
     &      +3D0*Mc/MW/SW*Zcc12*Zn14/sb)
      AAABR(509) = EE*Sqrt2*(4D0*Zcc12/CW*Zn11-3D0*Mc/MW
     &      /SW*Zcc11*Zn14/sb)
      AAABR(510) = EE*Sqrt2*(Zn31/CW*Zt11+3D0*Zn32/SW*Zt11
     &      +3D0*MtMM/MW/SW*Zn34*Zt12/sb)
      AAABR(511) = EE*Sqrt2*(4D0*Zn31/CW*Zt12-3D0*MtMM/MW
     &      /SW*Zn34*Zt11/sb)
      AAABR(512) = EE*Sqrt2*(Zn41/CW*Zt11+3D0*Zn42/SW*Zt11
     &      +3D0*MtMM/MW/SW*Zn44*Zt12/sb)
      AAABR(513) = EE*Sqrt2*(4D0*Zn41/CW*Zt12-3D0*MtMM/MW
     &      /SW*Zn44*Zt11/sb)
      AAABR(514) = EE*Sqrt2*(Zn51/CW*Zuu11+3D0*Zn52/SW*Zuu11
     &      +3D0*Mq/MW/SW*Zn54*Zuu12/sb)
      AAABR(515) = EE*Sqrt2*(4D0*Zn51/CW*Zuu12-3D0*Mq/MW
     &      /SW*Zn54*Zuu11/sb)
      AAABR(516) = EE*Sqrt2*(Zcc11/CW*Zn31+3D0*Zcc11/SW*Zn32
     &      +3D0*Mc/MW/SW*Zcc12*Zn34/sb)
      AAABR(517) = EE*Sqrt2*(4D0*Zcc12/CW*Zn31-3D0*Mc/MW
     &      /SW*Zcc11*Zn34/sb)
      AAABR(518) = EE*Sqrt2*(Zn41/CW*Zuu11+3D0*Zn42/SW*Zuu11
     &      +3D0*Mq/MW/SW*Zn44*Zuu12/sb)
      AAABR(519) = EE*Sqrt2*(4D0*Zn41/CW*Zuu12-3D0*Mq/MW
     &      /SW*Zn44*Zuu11/sb)
      AAABR(520) = EE*Sqrt2*(Zn31/CW*Zuu11+3D0*Zn32/SW*Zuu11
     &      +3D0*Mq/MW/SW*Zn34*Zuu12/sb)
      AAABR(521) = EE*Sqrt2*(4D0*Zn31/CW*Zuu12-3D0*Mq/MW
     &      /SW*Zn34*Zuu11/sb)
      AAABR(522) = EE*Sqrt2*(Zn21/CW*Zuu11+3D0*Zn22/SW*Zuu11
     &      +3D0*Mq/MW/SW*Zn24*Zuu12/sb)
      AAABR(523) = EE*Sqrt2*(4D0*Zn21/CW*Zuu12-3D0*Mq/MW
     &      /SW*Zn24*Zuu11/sb)
      AAABR(524) = EE*Sqrt2*(Zn11/CW*Zuu21+3D0*Zn12/SW*Zuu21
     &      +3D0*Mq/MW/SW*Zn14*Zuu22/sb)
      AAABR(525) = EE*Sqrt2*(4D0*Zn11/CW*Zuu22-3D0*Mq/MW
     &      /SW*Zn14*Zuu21/sb)
      AAABR(526) = EE*Sqrt2*(Zcc21/CW*Zn21+3D0*Zcc21/SW*Zn22
     &      +3D0*Mc/MW/SW*Zcc22*Zn24/sb)
      AAABR(527) = EE*Sqrt2*(4D0*Zcc22/CW*Zn21-3D0*Mc/MW
     &      /SW*Zcc21*Zn24/sb)
      AAABR(528) = EE*Sqrt2*(Zcc21/CW*Zn51+3D0*Zcc21/SW*Zn52
     &      +3D0*Mc/MW/SW*Zcc22*Zn54/sb)
      AAABR(529) = EE*Sqrt2*(4D0*Zcc22/CW*Zn51-3D0*Mc/MW
     &      /SW*Zcc21*Zn54/sb)
      AAABR(530) = EE*Sqrt2*(Zn51/CW*Zuu21+3D0*Zn52/SW*Zuu21
     &      +3D0*Mq/MW/SW*Zn54*Zuu22/sb)
      AAABR(531) = EE*Sqrt2*(4D0*Zn51/CW*Zuu22-3D0*Mq/MW
     &      /SW*Zn54*Zuu21/sb)
      AAABR(532) = EE*Sqrt2*(Zn21/CW*Zuu21+3D0*Zn22/SW*Zuu21
     &      +3D0*Mq/MW/SW*Zn24*Zuu22/sb)
      AAABR(533) = EE*Sqrt2*(4D0*Zn21/CW*Zuu22-3D0*Mq/MW
     &      /SW*Zn24*Zuu21/sb)
      AAABR(534) = EE*Sqrt2*(Zcc21/CW*Zn41+3D0*Zcc21/SW*Zn42
     &      +3D0*Mc/MW/SW*Zcc22*Zn44/sb)
      AAABR(535) = EE*Sqrt2*(4D0*Zcc22/CW*Zn41-3D0*Mc/MW
     &      /SW*Zcc21*Zn44/sb)
      AAABR(536) = EE*Sqrt2*(Zcc21/CW*Zn31+3D0*Zcc21/SW*Zn32
     &      +3D0*Mc/MW/SW*Zcc22*Zn34/sb)
      AAABR(537) = EE*Sqrt2*(4D0*Zcc22/CW*Zn31-3D0*Mc/MW
     &      /SW*Zcc21*Zn34/sb)
      AAABR(538) = EE*Sqrt2*(Zn31/CW*Zt21+3D0*Zn32/SW*Zt21
     &      +3D0*MtMM/MW/SW*Zn34*Zt22/sb)
      AAABR(539) = EE*Sqrt2*(4D0*Zn31/CW*Zt22-3D0*MtMM/MW
     &      /SW*Zn34*Zt21/sb)
      AAABR(540) = EE*Sqrt2*(Zn51/CW*Zt21+3D0*Zn52/SW*Zt21
     &      +3D0*MtMM/MW/SW*Zn54*Zt22/sb)
      AAABR(541) = EE*Sqrt2*(4D0*Zn51/CW*Zt22-3D0*MtMM/MW
     &      /SW*Zn54*Zt21/sb)
      AAABR(542) = EE*Sqrt2*(Zcc21/CW*Zn11+3D0*Zcc21/SW*Zn12
     &      +3D0*Mc/MW/SW*Zcc22*Zn14/sb)
      AAABR(543) = EE*Sqrt2*(4D0*Zcc22/CW*Zn11-3D0*Mc/MW
     &      /SW*Zcc21*Zn14/sb)
      AAABR(544) = EE*Sqrt2*(Zn11/CW*Zt21+3D0*Zn12/SW*Zt21
     &      +3D0*MtMM/MW/SW*Zn14*Zt22/sb)
      AAABR(545) = EE*Sqrt2*(4D0*Zn11/CW*Zt22-3D0*MtMM/MW
     &      /SW*Zn14*Zt21/sb)
      AAABR(546) = EE*Sqrt2*(Zn41/CW*Zuu21+3D0*Zn42/SW*Zuu21
     &      +3D0*Mq/MW/SW*Zn44*Zuu22/sb)
      AAABR(547) = EE*Sqrt2*(4D0*Zn41/CW*Zuu22-3D0*Mq/MW
     &      /SW*Zn44*Zuu21/sb)
      AAABR(548) = EE*Sqrt2*(Zn21/CW*Zt21+3D0*Zn22/SW*Zt21
     &      +3D0*MtMM/MW/SW*Zn24*Zt22/sb)
      AAABR(549) = EE*Sqrt2*(4D0*Zn21/CW*Zt22-3D0*MtMM/MW
     &      /SW*Zn24*Zt21/sb)
      AAABR(550) = EE*Sqrt2*(Zn31/CW*Zuu21+3D0*Zn32/SW*Zuu21
     &      +3D0*Mq/MW/SW*Zn34*Zuu22/sb)
      AAABR(551) = EE*Sqrt2*(4D0*Zn31/CW*Zuu22-3D0*Mq/MW
     &      /SW*Zn34*Zuu21/sb)
      AAABR(552) = EE*Sqrt2*(Zn41/CW*Zt21+3D0*Zn42/SW*Zt21
     &      +3D0*MtMM/MW/SW*Zn44*Zt22/sb)
      AAABR(553) = EE*Sqrt2*(4D0*Zn41/CW*Zt22-3D0*MtMM/MW
     &      /SW*Zn44*Zt21/sb)
      AAABR(554) = EE/MW*Mt/SW
      AAABR(555) = EE/MW*Mc/SW
      AAABR(556) = EE/MW*Mt*Pa11/SW/tB
      AAABR(557) = EE/MW*Mt*Pa21/SW/tB
      AAABR(558) = EE/MW*Mq*Pa21/SW/tB
      AAABR(559) = EE/MW*Mc*Pa11/SW/tB
      AAABR(560) = EE/MW*Mc*Pa21/SW/tB
      AAABR(561) = EE/MW*Mq*Pa11/SW/tB
      AAABR(562) = EE/MW*Mt/SW*Zh32/sb
      AAABR(563) = EE/MW*Mq/SW*Zh32/sb
      AAABR(564) = EE/MW*Mq/SW*Zh12/sb
      AAABR(565) = EE/MW*Mc/SW*Zh12/sb
      AAABR(566) = EE/MW*Mc/SW*Zh22/sb
      AAABR(567) = EE/MW*Mt/SW*Zh22/sb
      AAABR(568) = EE/MW*Mc/SW*Zh32/sb
      AAABR(569) = EE/MW*Mq/SW*Zh22/sb
      AAABR(570) = EE/MW*Mt/SW*Zh12/sb
      AAABR(571) = GG*Sqrt2*Zt11
      AAABR(572) = GG*Sqrt2*Zt12
      AAABR(573) = GG*Sqrt2*Zuu11
      AAABR(574) = GG*Sqrt2*Zuu12
      AAABR(575) = GG*Sqrt2*Zcc11
      AAABR(576) = GG*Sqrt2*Zcc12
      AAABR(577) = GG*Sqrt2*Zt21
      AAABR(578) = GG*Sqrt2*Zt22
      AAABR(579) = GG*Sqrt2*Zuu21
      AAABR(580) = GG*Sqrt2*Zuu22
      AAABR(581) = GG*Sqrt2*Zcc21
      AAABR(582) = GG*Sqrt2*Zcc22
      AAABR(583) = EE/MW/SW*(NMM13*Zn11*Zn13-NMM14*Zn11*Zn14
     &      +NMM23*Zn12*Zn13-NMM24*Zn12*Zn14+NMM35*Zn13*Zn15
     &      -NMM45*Zn14*Zn15)
      AAABR(584) = EE/MW/SW*(NMM13*Zn11*Zn43+NMM13*Zn13*Zn41
     &      -NMM14*Zn11*Zn44-NMM14*Zn14*Zn41+NMM23*Zn12*Zn43
     &      +NMM23*Zn13*Zn42-NMM24*Zn12*Zn44-NMM24*Zn14*Zn42
     &      +NMM35*Zn13*Zn45+NMM35*Zn15*Zn43-NMM45*Zn14*Zn45
     &      -NMM45*Zn15*Zn44)
      AAABR(585) = EE/MW/SW*(NMM13*Zn11*Zn53+NMM13*Zn13*Zn51
     &      -NMM14*Zn11*Zn54-NMM14*Zn14*Zn51+NMM23*Zn12*Zn53
     &      +NMM23*Zn13*Zn52-NMM24*Zn12*Zn54-NMM24*Zn14*Zn52
     &      +NMM35*Zn13*Zn55+NMM35*Zn15*Zn53-NMM45*Zn14*Zn55
     &      -NMM45*Zn15*Zn54)
      AAABR(586) = EE/MW/SW*(NMM13*Zn11*Zn23+NMM13*Zn13*Zn21
     &      -NMM14*Zn11*Zn24-NMM14*Zn14*Zn21+NMM23*Zn12*Zn23
     &      +NMM23*Zn13*Zn22-NMM24*Zn12*Zn24-NMM24*Zn14*Zn22
     &      +NMM35*Zn13*Zn25+NMM35*Zn15*Zn23-NMM45*Zn14*Zn25
     &      -NMM45*Zn15*Zn24)
      AAABR(587) = EE/MW/SW*(NMM13*Zn11*Zn33+NMM13*Zn13*Zn31
     &      -NMM14*Zn11*Zn34-NMM14*Zn14*Zn31+NMM23*Zn12*Zn33
     &      +NMM23*Zn13*Zn32-NMM24*Zn12*Zn34-NMM24*Zn14*Zn32
     &      +NMM35*Zn13*Zn35+NMM35*Zn15*Zn33-NMM45*Zn14*Zn35
     &      -NMM45*Zn15*Zn34)
      AAABR(588) = NMM34*Pa12*Sqrt2*Zn13*Zn54/xvev+NMM34
     &      *Pa12*Sqrt2*Zn14*Zn53/xvev+NMM55*Pa12*Sqrt2*Zn15
     &      *Zn55/xvev-Pa12*Sqrt2*Zn15*Zn55*muP/xvev-EE/MW
     &      *NMM13*Pa11/SW*Zn11*Zn53/cb*sb-EE/MW*NMM13*Pa11
     &      /SW*Zn13*Zn51/cb*sb-EE/MW*NMM14*Pa11/SW*Zn11
     &      *Zn54/cb/sb+EE/MW*NMM14*Pa11/SW*Zn11*Zn54/cb
     &      *sb-EE/MW*NMM14*Pa11/SW*Zn14*Zn51/cb/sb+EE/MW
     &      *NMM14*Pa11/SW*Zn14*Zn51/cb*sb-EE/MW*NMM23*Pa11
     &      /SW*Zn12*Zn53/cb*sb-EE/MW*NMM23*Pa11/SW*Zn13
     &      *Zn52/cb*sb-EE/MW*NMM24*Pa11/SW*Zn12*Zn54/cb
     &      /sb+EE/MW*NMM24*Pa11/SW*Zn12*Zn54/cb*sb-EE/MW
     &      *NMM24*Pa11/SW*Zn14*Zn52/cb/sb+EE/MW*NMM24*Pa11
     &      /SW*Zn14*Zn52/cb*sb+EE/MW*NMM35*Pa11/SW*Zn13
     &      *Zn55/cb/sb-EE/MW*NMM35*Pa11/SW*Zn13*Zn55/cb
     &      *sb+EE/MW*NMM35*Pa11/SW*Zn15*Zn53/cb/sb-EE/MW
     &      *NMM35*Pa11/SW*Zn15*Zn53/cb*sb+EE/MW*NMM45*Pa11
     &      /SW*Zn14*Zn55/cb*sb+EE/MW*NMM45*Pa11/SW*Zn15
     &      *Zn54/cb*sb
      AAABR(589) = NMM34*Pa22*Sqrt2*Zn13*Zn24/xvev+NMM34
     &      *Pa22*Sqrt2*Zn14*Zn23/xvev+NMM55*Pa22*Sqrt2*Zn15
     &      *Zn25/xvev-Pa22*Sqrt2*Zn15*Zn25*muP/xvev-EE/MW
     &      *NMM13*Pa21/SW*Zn11*Zn23/cb*sb-EE/MW*NMM13*Pa21
     &      /SW*Zn13*Zn21/cb*sb-EE/MW*NMM14*Pa21/SW*Zn11
     &      *Zn24/cb/sb+EE/MW*NMM14*Pa21/SW*Zn11*Zn24/cb
     &      *sb-EE/MW*NMM14*Pa21/SW*Zn14*Zn21/cb/sb+EE/MW
     &      *NMM14*Pa21/SW*Zn14*Zn21/cb*sb-EE/MW*NMM23*Pa21
     &      /SW*Zn12*Zn23/cb*sb-EE/MW*NMM23*Pa21/SW*Zn13
     &      *Zn22/cb*sb-EE/MW*NMM24*Pa21/SW*Zn12*Zn24/cb
     &      /sb+EE/MW*NMM24*Pa21/SW*Zn12*Zn24/cb*sb-EE/MW
     &      *NMM24*Pa21/SW*Zn14*Zn22/cb/sb+EE/MW*NMM24*Pa21
     &      /SW*Zn14*Zn22/cb*sb+EE/MW*NMM35*Pa21/SW*Zn13
     &      *Zn25/cb/sb-EE/MW*NMM35*Pa21/SW*Zn13*Zn25/cb
     &      *sb+EE/MW*NMM35*Pa21/SW*Zn15*Zn23/cb/sb-EE/MW
     &      *NMM35*Pa21/SW*Zn15*Zn23/cb*sb+EE/MW*NMM45*Pa21
     &      /SW*Zn14*Zn25/cb*sb+EE/MW*NMM45*Pa21/SW*Zn15
     &      *Zn24/cb*sb
      AAABR(590) = NMM34*Pa22*Sqrt2*Zn13*Zn54/xvev+NMM34
     &      *Pa22*Sqrt2*Zn14*Zn53/xvev+NMM55*Pa22*Sqrt2*Zn15
     &      *Zn55/xvev-Pa22*Sqrt2*Zn15*Zn55*muP/xvev-EE/MW
     &      *NMM13*Pa21/SW*Zn11*Zn53/cb*sb-EE/MW*NMM13*Pa21
     &      /SW*Zn13*Zn51/cb*sb-EE/MW*NMM14*Pa21/SW*Zn11
     &      *Zn54/cb/sb+EE/MW*NMM14*Pa21/SW*Zn11*Zn54/cb
     &      *sb-EE/MW*NMM14*Pa21/SW*Zn14*Zn51/cb/sb+EE/MW
     &      *NMM14*Pa21/SW*Zn14*Zn51/cb*sb-EE/MW*NMM23*Pa21
     &      /SW*Zn12*Zn53/cb*sb-EE/MW*NMM23*Pa21/SW*Zn13
     &      *Zn52/cb*sb-EE/MW*NMM24*Pa21/SW*Zn12*Zn54/cb
     &      /sb+EE/MW*NMM24*Pa21/SW*Zn12*Zn54/cb*sb-EE/MW
     &      *NMM24*Pa21/SW*Zn14*Zn52/cb/sb+EE/MW*NMM24*Pa21
     &      /SW*Zn14*Zn52/cb*sb+EE/MW*NMM35*Pa21/SW*Zn13
     &      *Zn55/cb/sb-EE/MW*NMM35*Pa21/SW*Zn13*Zn55/cb
     &      *sb+EE/MW*NMM35*Pa21/SW*Zn15*Zn53/cb/sb-EE/MW
     &      *NMM35*Pa21/SW*Zn15*Zn53/cb*sb+EE/MW*NMM45*Pa21
     &      /SW*Zn14*Zn55/cb*sb+EE/MW*NMM45*Pa21/SW*Zn15
     &      *Zn54/cb*sb
      AAABR(591) = NMM34*Pa12*Sqrt2*Zn13*Zn34/xvev+NMM34
     &      *Pa12*Sqrt2*Zn14*Zn33/xvev+NMM55*Pa12*Sqrt2*Zn15
     &      *Zn35/xvev-Pa12*Sqrt2*Zn15*Zn35*muP/xvev-EE/MW
     &      *NMM13*Pa11/SW*Zn11*Zn33/cb*sb-EE/MW*NMM13*Pa11
     &      /SW*Zn13*Zn31/cb*sb-EE/MW*NMM14*Pa11/SW*Zn11
     &      *Zn34/cb/sb+EE/MW*NMM14*Pa11/SW*Zn11*Zn34/cb
     &      *sb-EE/MW*NMM14*Pa11/SW*Zn14*Zn31/cb/sb+EE/MW
     &      *NMM14*Pa11/SW*Zn14*Zn31/cb*sb-EE/MW*NMM23*Pa11
     &      /SW*Zn12*Zn33/cb*sb-EE/MW*NMM23*Pa11/SW*Zn13
     &      *Zn32/cb*sb-EE/MW*NMM24*Pa11/SW*Zn12*Zn34/cb
     &      /sb+EE/MW*NMM24*Pa11/SW*Zn12*Zn34/cb*sb-EE/MW
     &      *NMM24*Pa11/SW*Zn14*Zn32/cb/sb+EE/MW*NMM24*Pa11
     &      /SW*Zn14*Zn32/cb*sb+EE/MW*NMM35*Pa11/SW*Zn13
     &      *Zn35/cb/sb-EE/MW*NMM35*Pa11/SW*Zn13*Zn35/cb
     &      *sb+EE/MW*NMM35*Pa11/SW*Zn15*Zn33/cb/sb-EE/MW
     &      *NMM35*Pa11/SW*Zn15*Zn33/cb*sb+EE/MW*NMM45*Pa11
     &      /SW*Zn14*Zn35/cb*sb+EE/MW*NMM45*Pa11/SW*Zn15
     &      *Zn34/cb*sb
      AAABR(592) = NMM55*Pa12*Sqrt2*Zn15**2/xvev-Pa12*Sqrt2
     &      *Zn15**2*muP/xvev+2D0*NMM34*Pa12*Sqrt2*Zn13*Zn14
     &      /xvev-2D0*EE/MW*NMM13*Pa11/SW*Zn11*Zn13/cb*sb
     &      -2D0*EE/MW*NMM14*Pa11/SW*Zn11*Zn14/cb/sb+2D0
     &      *EE/MW*NMM14*Pa11/SW*Zn11*Zn14/cb*sb-2D0*EE/MW
     &      *NMM23*Pa11/SW*Zn12*Zn13/cb*sb-2D0*EE/MW*NMM24
     &      *Pa11/SW*Zn12*Zn14/cb/sb+2D0*EE/MW*NMM24*Pa11
     &      /SW*Zn12*Zn14/cb*sb+2D0*EE/MW*NMM35*Pa11/SW*Zn13
     &      *Zn15/cb/sb-2D0*EE/MW*NMM35*Pa11/SW*Zn13*Zn15
     &      /cb*sb+2D0*EE/MW*NMM45*Pa11/SW*Zn14*Zn15/cb*sb
      AAABR(593) = NMM34*Pa12*Sqrt2*Zn13*Zn44/xvev+NMM34
     &      *Pa12*Sqrt2*Zn14*Zn43/xvev+NMM55*Pa12*Sqrt2*Zn15
     &      *Zn45/xvev-Pa12*Sqrt2*Zn15*Zn45*muP/xvev-EE/MW
     &      *NMM13*Pa11/SW*Zn11*Zn43/cb*sb-EE/MW*NMM13*Pa11
     &      /SW*Zn13*Zn41/cb*sb-EE/MW*NMM14*Pa11/SW*Zn11
     &      *Zn44/cb/sb+EE/MW*NMM14*Pa11/SW*Zn11*Zn44/cb
     &      *sb-EE/MW*NMM14*Pa11/SW*Zn14*Zn41/cb/sb+EE/MW
     &      *NMM14*Pa11/SW*Zn14*Zn41/cb*sb-EE/MW*NMM23*Pa11
     &      /SW*Zn12*Zn43/cb*sb-EE/MW*NMM23*Pa11/SW*Zn13
     &      *Zn42/cb*sb-EE/MW*NMM24*Pa11/SW*Zn12*Zn44/cb
     &      /sb+EE/MW*NMM24*Pa11/SW*Zn12*Zn44/cb*sb-EE/MW
     &      *NMM24*Pa11/SW*Zn14*Zn42/cb/sb+EE/MW*NMM24*Pa11
     &      /SW*Zn14*Zn42/cb*sb+EE/MW*NMM35*Pa11/SW*Zn13
     &      *Zn45/cb/sb-EE/MW*NMM35*Pa11/SW*Zn13*Zn45/cb
     &      *sb+EE/MW*NMM35*Pa11/SW*Zn15*Zn43/cb/sb-EE/MW
     &      *NMM35*Pa11/SW*Zn15*Zn43/cb*sb+EE/MW*NMM45*Pa11
     &      /SW*Zn14*Zn45/cb*sb+EE/MW*NMM45*Pa11/SW*Zn15
     &      *Zn44/cb*sb
      AAABR(594) = NMM55*Pa22*Sqrt2*Zn15**2/xvev-Pa22*Sqrt2
     &      *Zn15**2*muP/xvev+2D0*NMM34*Pa22*Sqrt2*Zn13*Zn14
     &      /xvev-2D0*EE/MW*NMM13*Pa21/SW*Zn11*Zn13/cb*sb
     &      -2D0*EE/MW*NMM14*Pa21/SW*Zn11*Zn14/cb/sb+2D0
     &      *EE/MW*NMM14*Pa21/SW*Zn11*Zn14/cb*sb-2D0*EE/MW
     &      *NMM23*Pa21/SW*Zn12*Zn13/cb*sb-2D0*EE/MW*NMM24
     &      *Pa21/SW*Zn12*Zn14/cb/sb+2D0*EE/MW*NMM24*Pa21
     &      /SW*Zn12*Zn14/cb*sb+2D0*EE/MW*NMM35*Pa21/SW*Zn13
     &      *Zn15/cb/sb-2D0*EE/MW*NMM35*Pa21/SW*Zn13*Zn15
     &      /cb*sb+2D0*EE/MW*NMM45*Pa21/SW*Zn14*Zn15/cb*sb
      AAABR(595) = NMM34*Pa22*Sqrt2*Zn13*Zn34/xvev+NMM34
     &      *Pa22*Sqrt2*Zn14*Zn33/xvev+NMM55*Pa22*Sqrt2*Zn15
     &      *Zn35/xvev-Pa22*Sqrt2*Zn15*Zn35*muP/xvev-EE/MW
     &      *NMM13*Pa21/SW*Zn11*Zn33/cb*sb-EE/MW*NMM13*Pa21
     &      /SW*Zn13*Zn31/cb*sb-EE/MW*NMM14*Pa21/SW*Zn11
     &      *Zn34/cb/sb+EE/MW*NMM14*Pa21/SW*Zn11*Zn34/cb
     &      *sb-EE/MW*NMM14*Pa21/SW*Zn14*Zn31/cb/sb+EE/MW
     &      *NMM14*Pa21/SW*Zn14*Zn31/cb*sb-EE/MW*NMM23*Pa21
     &      /SW*Zn12*Zn33/cb*sb-EE/MW*NMM23*Pa21/SW*Zn13
     &      *Zn32/cb*sb-EE/MW*NMM24*Pa21/SW*Zn12*Zn34/cb
     &      /sb+EE/MW*NMM24*Pa21/SW*Zn12*Zn34/cb*sb-EE/MW
     &      *NMM24*Pa21/SW*Zn14*Zn32/cb/sb+EE/MW*NMM24*Pa21
     &      /SW*Zn14*Zn32/cb*sb+EE/MW*NMM35*Pa21/SW*Zn13
     &      *Zn35/cb/sb-EE/MW*NMM35*Pa21/SW*Zn13*Zn35/cb
     &      *sb+EE/MW*NMM35*Pa21/SW*Zn15*Zn33/cb/sb-EE/MW
     &      *NMM35*Pa21/SW*Zn15*Zn33/cb*sb+EE/MW*NMM45*Pa21
     &      /SW*Zn14*Zn35/cb*sb+EE/MW*NMM45*Pa21/SW*Zn15
     &      *Zn34/cb*sb
      AAABR(596) = NMM34*Pa22*Sqrt2*Zn13*Zn44/xvev+NMM34
     &      *Pa22*Sqrt2*Zn14*Zn43/xvev+NMM55*Pa22*Sqrt2*Zn15
     &      *Zn45/xvev-Pa22*Sqrt2*Zn15*Zn45*muP/xvev-EE/MW
     &      *NMM13*Pa21/SW*Zn11*Zn43/cb*sb-EE/MW*NMM13*Pa21
     &      /SW*Zn13*Zn41/cb*sb-EE/MW*NMM14*Pa21/SW*Zn11
     &      *Zn44/cb/sb+EE/MW*NMM14*Pa21/SW*Zn11*Zn44/cb
     &      *sb-EE/MW*NMM14*Pa21/SW*Zn14*Zn41/cb/sb+EE/MW
     &      *NMM14*Pa21/SW*Zn14*Zn41/cb*sb-EE/MW*NMM23*Pa21
     &      /SW*Zn12*Zn43/cb*sb-EE/MW*NMM23*Pa21/SW*Zn13
     &      *Zn42/cb*sb-EE/MW*NMM24*Pa21/SW*Zn12*Zn44/cb
     &      /sb+EE/MW*NMM24*Pa21/SW*Zn12*Zn44/cb*sb-EE/MW
     &      *NMM24*Pa21/SW*Zn14*Zn42/cb/sb+EE/MW*NMM24*Pa21
     &      /SW*Zn14*Zn42/cb*sb+EE/MW*NMM35*Pa21/SW*Zn13
     &      *Zn45/cb/sb-EE/MW*NMM35*Pa21/SW*Zn13*Zn45/cb
     &      *sb+EE/MW*NMM35*Pa21/SW*Zn15*Zn43/cb/sb-EE/MW
     &      *NMM35*Pa21/SW*Zn15*Zn43/cb*sb+EE/MW*NMM45*Pa21
     &      /SW*Zn14*Zn45/cb*sb+EE/MW*NMM45*Pa21/SW*Zn15
     &      *Zn44/cb*sb
      AAABR(597) = NMM34*Pa12*Sqrt2*Zn13*Zn24/xvev+NMM34
     &      *Pa12*Sqrt2*Zn14*Zn23/xvev+NMM55*Pa12*Sqrt2*Zn15
     &      *Zn25/xvev-Pa12*Sqrt2*Zn15*Zn25*muP/xvev-EE/MW
     &      *NMM13*Pa11/SW*Zn11*Zn23/cb*sb-EE/MW*NMM13*Pa11
     &      /SW*Zn13*Zn21/cb*sb-EE/MW*NMM14*Pa11/SW*Zn11
     &      *Zn24/cb/sb+EE/MW*NMM14*Pa11/SW*Zn11*Zn24/cb
     &      *sb-EE/MW*NMM14*Pa11/SW*Zn14*Zn21/cb/sb+EE/MW
     &      *NMM14*Pa11/SW*Zn14*Zn21/cb*sb-EE/MW*NMM23*Pa11
     &      /SW*Zn12*Zn23/cb*sb-EE/MW*NMM23*Pa11/SW*Zn13
     &      *Zn22/cb*sb-EE/MW*NMM24*Pa11/SW*Zn12*Zn24/cb
     &      /sb+EE/MW*NMM24*Pa11/SW*Zn12*Zn24/cb*sb-EE/MW
     &      *NMM24*Pa11/SW*Zn14*Zn22/cb/sb+EE/MW*NMM24*Pa11
     &      /SW*Zn14*Zn22/cb*sb+EE/MW*NMM35*Pa11/SW*Zn13
     &      *Zn25/cb/sb-EE/MW*NMM35*Pa11/SW*Zn13*Zn25/cb
     &      *sb+EE/MW*NMM35*Pa11/SW*Zn15*Zn23/cb/sb-EE/MW
     &      *NMM35*Pa11/SW*Zn15*Zn23/cb*sb+EE/MW*NMM45*Pa11
     &      /SW*Zn14*Zn25/cb*sb+EE/MW*NMM45*Pa11/SW*Zn15
     &      *Zn24/cb*sb
      AAABR(598) = NMM34*Sqrt2*Zh33*Zn13*Zn54/xvev+NMM34
     &      *Sqrt2*Zh33*Zn14*Zn53/xvev+NMM55*Sqrt2*Zh33*Zn15
     &      *Zn55/xvev-Sqrt2*Zh33*Zn15*Zn55*muP/xvev+EE/MW
     &      *NMM13/SW*Zh31*Zn11*Zn53/cb+EE/MW*NMM13/SW*Zh31
     &      *Zn13*Zn51/cb+EE/MW*NMM14/SW*Zh32*Zn11*Zn54/sb
     &      +EE/MW*NMM14/SW*Zh32*Zn14*Zn51/sb+EE/MW*NMM23
     &      /SW*Zh31*Zn12*Zn53/cb+EE/MW*NMM23/SW*Zh31*Zn13
     &      *Zn52/cb+EE/MW*NMM24/SW*Zh32*Zn12*Zn54/sb+EE
     &      /MW*NMM24/SW*Zh32*Zn14*Zn52/sb+EE/MW*NMM35/SW
     &      *Zh32*Zn13*Zn55/sb+EE/MW*NMM35/SW*Zh32*Zn15*Zn53
     &      /sb+EE/MW*NMM45/SW*Zh31*Zn14*Zn55/cb+EE/MW*NMM45
     &      /SW*Zh31*Zn15*Zn54/cb
      AAABR(599) = NMM34*Sqrt2*Zh33*Zn13*Zn34/xvev+NMM34
     &      *Sqrt2*Zh33*Zn14*Zn33/xvev+NMM55*Sqrt2*Zh33*Zn15
     &      *Zn35/xvev-Sqrt2*Zh33*Zn15*Zn35*muP/xvev+EE/MW
     &      *NMM13/SW*Zh31*Zn11*Zn33/cb+EE/MW*NMM13/SW*Zh31
     &      *Zn13*Zn31/cb+EE/MW*NMM14/SW*Zh32*Zn11*Zn34/sb
     &      +EE/MW*NMM14/SW*Zh32*Zn14*Zn31/sb+EE/MW*NMM23
     &      /SW*Zh31*Zn12*Zn33/cb+EE/MW*NMM23/SW*Zh31*Zn13
     &      *Zn32/cb+EE/MW*NMM24/SW*Zh32*Zn12*Zn34/sb+EE
     &      /MW*NMM24/SW*Zh32*Zn14*Zn32/sb+EE/MW*NMM35/SW
     &      *Zh32*Zn13*Zn35/sb+EE/MW*NMM35/SW*Zh32*Zn15*Zn33
     &      /sb+EE/MW*NMM45/SW*Zh31*Zn14*Zn35/cb+EE/MW*NMM45
     &      /SW*Zh31*Zn15*Zn34/cb
      AAABR(600) = NMM34*Sqrt2*Zh13*Zn13*Zn34/xvev+NMM34
     &      *Sqrt2*Zh13*Zn14*Zn33/xvev+NMM55*Sqrt2*Zh13*Zn15
     &      *Zn35/xvev-Sqrt2*Zh13*Zn15*Zn35*muP/xvev+EE/MW
     &      *NMM13/SW*Zh11*Zn11*Zn33/cb+EE/MW*NMM13/SW*Zh11
     &      *Zn13*Zn31/cb+EE/MW*NMM14/SW*Zh12*Zn11*Zn34/sb
     &      +EE/MW*NMM14/SW*Zh12*Zn14*Zn31/sb+EE/MW*NMM23
     &      /SW*Zh11*Zn12*Zn33/cb+EE/MW*NMM23/SW*Zh11*Zn13
     &      *Zn32/cb+EE/MW*NMM24/SW*Zh12*Zn12*Zn34/sb+EE
     &      /MW*NMM24/SW*Zh12*Zn14*Zn32/sb+EE/MW*NMM35/SW
     &      *Zh12*Zn13*Zn35/sb+EE/MW*NMM35/SW*Zh12*Zn15*Zn33
     &      /sb+EE/MW*NMM45/SW*Zh11*Zn14*Zn35/cb+EE/MW*NMM45
     &      /SW*Zh11*Zn15*Zn34/cb
      AAABR(601) = NMM34*Sqrt2*Zh23*Zn13*Zn54/xvev+NMM34
     &      *Sqrt2*Zh23*Zn14*Zn53/xvev+NMM55*Sqrt2*Zh23*Zn15
     &      *Zn55/xvev-Sqrt2*Zh23*Zn15*Zn55*muP/xvev+EE/MW
     &      *NMM13/SW*Zh21*Zn11*Zn53/cb+EE/MW*NMM13/SW*Zh21
     &      *Zn13*Zn51/cb+EE/MW*NMM14/SW*Zh22*Zn11*Zn54/sb
     &      +EE/MW*NMM14/SW*Zh22*Zn14*Zn51/sb+EE/MW*NMM23
     &      /SW*Zh21*Zn12*Zn53/cb+EE/MW*NMM23/SW*Zh21*Zn13
     &      *Zn52/cb+EE/MW*NMM24/SW*Zh22*Zn12*Zn54/sb+EE
     &      /MW*NMM24/SW*Zh22*Zn14*Zn52/sb+EE/MW*NMM35/SW
     &      *Zh22*Zn13*Zn55/sb+EE/MW*NMM35/SW*Zh22*Zn15*Zn53
     &      /sb+EE/MW*NMM45/SW*Zh21*Zn14*Zn55/cb+EE/MW*NMM45
     &      /SW*Zh21*Zn15*Zn54/cb
      AAABR(602) = NMM34*Sqrt2*Zh23*Zn13*Zn44/xvev+NMM34
     &      *Sqrt2*Zh23*Zn14*Zn43/xvev+NMM55*Sqrt2*Zh23*Zn15
     &      *Zn45/xvev-Sqrt2*Zh23*Zn15*Zn45*muP/xvev+EE/MW
     &      *NMM13/SW*Zh21*Zn11*Zn43/cb+EE/MW*NMM13/SW*Zh21
     &      *Zn13*Zn41/cb+EE/MW*NMM14/SW*Zh22*Zn11*Zn44/sb
     &      +EE/MW*NMM14/SW*Zh22*Zn14*Zn41/sb+EE/MW*NMM23
     &      /SW*Zh21*Zn12*Zn43/cb+EE/MW*NMM23/SW*Zh21*Zn13
     &      *Zn42/cb+EE/MW*NMM24/SW*Zh22*Zn12*Zn44/sb+EE
     &      /MW*NMM24/SW*Zh22*Zn14*Zn42/sb+EE/MW*NMM35/SW
     &      *Zh22*Zn13*Zn45/sb+EE/MW*NMM35/SW*Zh22*Zn15*Zn43
     &      /sb+EE/MW*NMM45/SW*Zh21*Zn14*Zn45/cb+EE/MW*NMM45
     &      /SW*Zh21*Zn15*Zn44/cb
      AAABR(603) = NMM34*Sqrt2*Zh23*Zn13*Zn24/xvev+NMM34
     &      *Sqrt2*Zh23*Zn14*Zn23/xvev+NMM55*Sqrt2*Zh23*Zn15
     &      *Zn25/xvev-Sqrt2*Zh23*Zn15*Zn25*muP/xvev+EE/MW
     &      *NMM13/SW*Zh21*Zn11*Zn23/cb+EE/MW*NMM13/SW*Zh21
     &      *Zn13*Zn21/cb+EE/MW*NMM14/SW*Zh22*Zn11*Zn24/sb
     &      +EE/MW*NMM14/SW*Zh22*Zn14*Zn21/sb+EE/MW*NMM23
     &      /SW*Zh21*Zn12*Zn23/cb+EE/MW*NMM23/SW*Zh21*Zn13
     &      *Zn22/cb+EE/MW*NMM24/SW*Zh22*Zn12*Zn24/sb+EE
     &      /MW*NMM24/SW*Zh22*Zn14*Zn22/sb+EE/MW*NMM35/SW
     &      *Zh22*Zn13*Zn25/sb+EE/MW*NMM35/SW*Zh22*Zn15*Zn23
     &      /sb+EE/MW*NMM45/SW*Zh21*Zn14*Zn25/cb+EE/MW*NMM45
     &      /SW*Zh21*Zn15*Zn24/cb
      AAABR(604) = NMM34*Sqrt2*Zh33*Zn13*Zn44/xvev+NMM34
     &      *Sqrt2*Zh33*Zn14*Zn43/xvev+NMM55*Sqrt2*Zh33*Zn15
     &      *Zn45/xvev-Sqrt2*Zh33*Zn15*Zn45*muP/xvev+EE/MW
     &      *NMM13/SW*Zh31*Zn11*Zn43/cb+EE/MW*NMM13/SW*Zh31
     &      *Zn13*Zn41/cb+EE/MW*NMM14/SW*Zh32*Zn11*Zn44/sb
     &      +EE/MW*NMM14/SW*Zh32*Zn14*Zn41/sb+EE/MW*NMM23
     &      /SW*Zh31*Zn12*Zn43/cb+EE/MW*NMM23/SW*Zh31*Zn13
     &      *Zn42/cb+EE/MW*NMM24/SW*Zh32*Zn12*Zn44/sb+EE
     &      /MW*NMM24/SW*Zh32*Zn14*Zn42/sb+EE/MW*NMM35/SW
     &      *Zh32*Zn13*Zn45/sb+EE/MW*NMM35/SW*Zh32*Zn15*Zn43
     &      /sb+EE/MW*NMM45/SW*Zh31*Zn14*Zn45/cb+EE/MW*NMM45
     &      /SW*Zh31*Zn15*Zn44/cb
      AAABR(605) = NMM34*Sqrt2*Zh13*Zn13*Zn24/xvev+NMM34
     &      *Sqrt2*Zh13*Zn14*Zn23/xvev+NMM55*Sqrt2*Zh13*Zn15
     &      *Zn25/xvev-Sqrt2*Zh13*Zn15*Zn25*muP/xvev+EE/MW
     &      *NMM13/SW*Zh11*Zn11*Zn23/cb+EE/MW*NMM13/SW*Zh11
     &      *Zn13*Zn21/cb+EE/MW*NMM14/SW*Zh12*Zn11*Zn24/sb
     &      +EE/MW*NMM14/SW*Zh12*Zn14*Zn21/sb+EE/MW*NMM23
     &      /SW*Zh11*Zn12*Zn23/cb+EE/MW*NMM23/SW*Zh11*Zn13
     &      *Zn22/cb+EE/MW*NMM24/SW*Zh12*Zn12*Zn24/sb+EE
     &      /MW*NMM24/SW*Zh12*Zn14*Zn22/sb+EE/MW*NMM35/SW
     &      *Zh12*Zn13*Zn25/sb+EE/MW*NMM35/SW*Zh12*Zn15*Zn23
     &      /sb+EE/MW*NMM45/SW*Zh11*Zn14*Zn25/cb+EE/MW*NMM45
     &      /SW*Zh11*Zn15*Zn24/cb
      AAABR(606) = NMM34*Sqrt2*Zh13*Zn13*Zn54/xvev+NMM34
     &      *Sqrt2*Zh13*Zn14*Zn53/xvev+NMM55*Sqrt2*Zh13*Zn15
     &      *Zn55/xvev-Sqrt2*Zh13*Zn15*Zn55*muP/xvev+EE/MW
     &      *NMM13/SW*Zh11*Zn11*Zn53/cb+EE/MW*NMM13/SW*Zh11
     &      *Zn13*Zn51/cb+EE/MW*NMM14/SW*Zh12*Zn11*Zn54/sb
     &      +EE/MW*NMM14/SW*Zh12*Zn14*Zn51/sb+EE/MW*NMM23
     &      /SW*Zh11*Zn12*Zn53/cb+EE/MW*NMM23/SW*Zh11*Zn13
     &      *Zn52/cb+EE/MW*NMM24/SW*Zh12*Zn12*Zn54/sb+EE
     &      /MW*NMM24/SW*Zh12*Zn14*Zn52/sb+EE/MW*NMM35/SW
     &      *Zh12*Zn13*Zn55/sb+EE/MW*NMM35/SW*Zh12*Zn15*Zn53
     &      /sb+EE/MW*NMM45/SW*Zh11*Zn14*Zn55/cb+EE/MW*NMM45
     &      /SW*Zh11*Zn15*Zn54/cb
      AAABR(607) = NMM55*Sqrt2*Zh33*Zn15**2/xvev-Sqrt2*Zh33
     &      *Zn15**2*muP/xvev+2D0*NMM34*Sqrt2*Zh33*Zn13*Zn14
     &      /xvev+2D0*EE/MW*NMM13/SW*Zh31*Zn11*Zn13/cb+2D0
     &      *EE/MW*NMM14/SW*Zh32*Zn11*Zn14/sb+2D0*EE/MW*NMM23
     &      /SW*Zh31*Zn12*Zn13/cb+2D0*EE/MW*NMM24/SW*Zh32
     &      *Zn12*Zn14/sb+2D0*EE/MW*NMM35/SW*Zh32*Zn13*Zn15
     &      /sb+2D0*EE/MW*NMM45/SW*Zh31*Zn14*Zn15/cb
      AAABR(608) = NMM34*Sqrt2*Zh33*Zn13*Zn24/xvev+NMM34
     &      *Sqrt2*Zh33*Zn14*Zn23/xvev+NMM55*Sqrt2*Zh33*Zn15
     &      *Zn25/xvev-Sqrt2*Zh33*Zn15*Zn25*muP/xvev+EE/MW
     &      *NMM13/SW*Zh31*Zn11*Zn23/cb+EE/MW*NMM13/SW*Zh31
     &      *Zn13*Zn21/cb+EE/MW*NMM14/SW*Zh32*Zn11*Zn24/sb
     &      +EE/MW*NMM14/SW*Zh32*Zn14*Zn21/sb+EE/MW*NMM23
     &      /SW*Zh31*Zn12*Zn23/cb+EE/MW*NMM23/SW*Zh31*Zn13
     &      *Zn22/cb+EE/MW*NMM24/SW*Zh32*Zn12*Zn24/sb+EE
     &      /MW*NMM24/SW*Zh32*Zn14*Zn22/sb+EE/MW*NMM35/SW
     &      *Zh32*Zn13*Zn25/sb+EE/MW*NMM35/SW*Zh32*Zn15*Zn23
     &      /sb+EE/MW*NMM45/SW*Zh31*Zn14*Zn25/cb+EE/MW*NMM45
     &      /SW*Zh31*Zn15*Zn24/cb
      AAABR(609) = NMM34*Sqrt2*Zh23*Zn13*Zn34/xvev+NMM34
     &      *Sqrt2*Zh23*Zn14*Zn33/xvev+NMM55*Sqrt2*Zh23*Zn15
     &      *Zn35/xvev-Sqrt2*Zh23*Zn15*Zn35*muP/xvev+EE/MW
     &      *NMM13/SW*Zh21*Zn11*Zn33/cb+EE/MW*NMM13/SW*Zh21
     &      *Zn13*Zn31/cb+EE/MW*NMM14/SW*Zh22*Zn11*Zn34/sb
     &      +EE/MW*NMM14/SW*Zh22*Zn14*Zn31/sb+EE/MW*NMM23
     &      /SW*Zh21*Zn12*Zn33/cb+EE/MW*NMM23/SW*Zh21*Zn13
     &      *Zn32/cb+EE/MW*NMM24/SW*Zh22*Zn12*Zn34/sb+EE
     &      /MW*NMM24/SW*Zh22*Zn14*Zn32/sb+EE/MW*NMM35/SW
     &      *Zh22*Zn13*Zn35/sb+EE/MW*NMM35/SW*Zh22*Zn15*Zn33
     &      /sb+EE/MW*NMM45/SW*Zh21*Zn14*Zn35/cb+EE/MW*NMM45
     &      /SW*Zh21*Zn15*Zn34/cb
      AAABR(610) = NMM55*Sqrt2*Zh13*Zn15**2/xvev-Sqrt2*Zh13
     &      *Zn15**2*muP/xvev+2D0*NMM34*Sqrt2*Zh13*Zn13*Zn14
     &      /xvev+2D0*EE/MW*NMM13/SW*Zh11*Zn11*Zn13/cb+2D0
     &      *EE/MW*NMM14/SW*Zh12*Zn11*Zn14/sb+2D0*EE/MW*NMM23
     &      /SW*Zh11*Zn12*Zn13/cb+2D0*EE/MW*NMM24/SW*Zh12
     &      *Zn12*Zn14/sb+2D0*EE/MW*NMM35/SW*Zh12*Zn13*Zn15
     &      /sb+2D0*EE/MW*NMM45/SW*Zh11*Zn14*Zn15/cb
      AAABR(611) = NMM34*Sqrt2*Zh13*Zn13*Zn44/xvev+NMM34
     &      *Sqrt2*Zh13*Zn14*Zn43/xvev+NMM55*Sqrt2*Zh13*Zn15
     &      *Zn45/xvev-Sqrt2*Zh13*Zn15*Zn45*muP/xvev+EE/MW
     &      *NMM13/SW*Zh11*Zn11*Zn43/cb+EE/MW*NMM13/SW*Zh11
     &      *Zn13*Zn41/cb+EE/MW*NMM14/SW*Zh12*Zn11*Zn44/sb
     &      +EE/MW*NMM14/SW*Zh12*Zn14*Zn41/sb+EE/MW*NMM23
     &      /SW*Zh11*Zn12*Zn43/cb+EE/MW*NMM23/SW*Zh11*Zn13
     &      *Zn42/cb+EE/MW*NMM24/SW*Zh12*Zn12*Zn44/sb+EE
     &      /MW*NMM24/SW*Zh12*Zn14*Zn42/sb+EE/MW*NMM35/SW
     &      *Zh12*Zn13*Zn45/sb+EE/MW*NMM35/SW*Zh12*Zn15*Zn43
     &      /sb+EE/MW*NMM45/SW*Zh11*Zn14*Zn45/cb+EE/MW*NMM45
     &      /SW*Zh11*Zn15*Zn44/cb
      AAABR(612) = NMM55*Sqrt2*Zh23*Zn15**2/xvev-Sqrt2*Zh23
     &      *Zn15**2*muP/xvev+2D0*NMM34*Sqrt2*Zh23*Zn13*Zn14
     &      /xvev+2D0*EE/MW*NMM13/SW*Zh21*Zn11*Zn13/cb+2D0
     &      *EE/MW*NMM14/SW*Zh22*Zn11*Zn14/sb+2D0*EE/MW*NMM23
     &      /SW*Zh21*Zn12*Zn13/cb+2D0*EE/MW*NMM24/SW*Zh22
     &      *Zn12*Zn14/sb+2D0*EE/MW*NMM35/SW*Zh22*Zn13*Zn15
     &      /sb+2D0*EE/MW*NMM45/SW*Zh21*Zn14*Zn15/cb
      AAABR(613) = EE/MW/SW*(NMM13*Zn21*Zn33+NMM13*Zn23*Zn31
     &      -NMM14*Zn21*Zn34-NMM14*Zn24*Zn31+NMM23*Zn22*Zn33
     &      +NMM23*Zn23*Zn32-NMM24*Zn22*Zn34-NMM24*Zn24*Zn32
     &      +NMM35*Zn23*Zn35+NMM35*Zn25*Zn33-NMM45*Zn24*Zn35
     &      -NMM45*Zn25*Zn34)
      AAABR(614) = EE/MW/SW*(NMM13*Zn21*Zn43+NMM13*Zn23*Zn41
     &      -NMM14*Zn21*Zn44-NMM14*Zn24*Zn41+NMM23*Zn22*Zn43
     &      +NMM23*Zn23*Zn42-NMM24*Zn22*Zn44-NMM24*Zn24*Zn42
     &      +NMM35*Zn23*Zn45+NMM35*Zn25*Zn43-NMM45*Zn24*Zn45
     &      -NMM45*Zn25*Zn44)
      AAABR(615) = EE/MW/SW*(NMM13*Zn21*Zn23-NMM14*Zn21*Zn24
     &      +NMM23*Zn22*Zn23-NMM24*Zn22*Zn24+NMM35*Zn23*Zn25
     &      -NMM45*Zn24*Zn25)
      AAABR(616) = EE/MW/SW*(NMM13*Zn21*Zn53+NMM13*Zn23*Zn51
     &      -NMM14*Zn21*Zn54-NMM14*Zn24*Zn51+NMM23*Zn22*Zn53
     &      +NMM23*Zn23*Zn52-NMM24*Zn22*Zn54-NMM24*Zn24*Zn52
     &      +NMM35*Zn23*Zn55+NMM35*Zn25*Zn53-NMM45*Zn24*Zn55
     &      -NMM45*Zn25*Zn54)
      AAABR(617) = NMM34*Pa12*Sqrt2*Zn23*Zn34/xvev+NMM34
     &      *Pa12*Sqrt2*Zn24*Zn33/xvev+NMM55*Pa12*Sqrt2*Zn25
     &      *Zn35/xvev-Pa12*Sqrt2*Zn25*Zn35*muP/xvev-EE/MW
     &      *NMM13*Pa11/SW*Zn21*Zn33/cb*sb-EE/MW*NMM13*Pa11
     &      /SW*Zn23*Zn31/cb*sb-EE/MW*NMM14*Pa11/SW*Zn21
     &      *Zn34/cb/sb+EE/MW*NMM14*Pa11/SW*Zn21*Zn34/cb
     &      *sb-EE/MW*NMM14*Pa11/SW*Zn24*Zn31/cb/sb+EE/MW
     &      *NMM14*Pa11/SW*Zn24*Zn31/cb*sb-EE/MW*NMM23*Pa11
     &      /SW*Zn22*Zn33/cb*sb-EE/MW*NMM23*Pa11/SW*Zn23
     &      *Zn32/cb*sb-EE/MW*NMM24*Pa11/SW*Zn22*Zn34/cb
     &      /sb+EE/MW*NMM24*Pa11/SW*Zn22*Zn34/cb*sb-EE/MW
     &      *NMM24*Pa11/SW*Zn24*Zn32/cb/sb+EE/MW*NMM24*Pa11
     &      /SW*Zn24*Zn32/cb*sb+EE/MW*NMM35*Pa11/SW*Zn23
     &      *Zn35/cb/sb-EE/MW*NMM35*Pa11/SW*Zn23*Zn35/cb
     &      *sb+EE/MW*NMM35*Pa11/SW*Zn25*Zn33/cb/sb-EE/MW
     &      *NMM35*Pa11/SW*Zn25*Zn33/cb*sb+EE/MW*NMM45*Pa11
     &      /SW*Zn24*Zn35/cb*sb+EE/MW*NMM45*Pa11/SW*Zn25
     &      *Zn34/cb*sb
      AAABR(618) = NMM34*Pa22*Sqrt2*Zn23*Zn44/xvev+NMM34
     &      *Pa22*Sqrt2*Zn24*Zn43/xvev+NMM55*Pa22*Sqrt2*Zn25
     &      *Zn45/xvev-Pa22*Sqrt2*Zn25*Zn45*muP/xvev-EE/MW
     &      *NMM13*Pa21/SW*Zn21*Zn43/cb*sb-EE/MW*NMM13*Pa21
     &      /SW*Zn23*Zn41/cb*sb-EE/MW*NMM14*Pa21/SW*Zn21
     &      *Zn44/cb/sb+EE/MW*NMM14*Pa21/SW*Zn21*Zn44/cb
     &      *sb-EE/MW*NMM14*Pa21/SW*Zn24*Zn41/cb/sb+EE/MW
     &      *NMM14*Pa21/SW*Zn24*Zn41/cb*sb-EE/MW*NMM23*Pa21
     &      /SW*Zn22*Zn43/cb*sb-EE/MW*NMM23*Pa21/SW*Zn23
     &      *Zn42/cb*sb-EE/MW*NMM24*Pa21/SW*Zn22*Zn44/cb
     &      /sb+EE/MW*NMM24*Pa21/SW*Zn22*Zn44/cb*sb-EE/MW
     &      *NMM24*Pa21/SW*Zn24*Zn42/cb/sb+EE/MW*NMM24*Pa21
     &      /SW*Zn24*Zn42/cb*sb+EE/MW*NMM35*Pa21/SW*Zn23
     &      *Zn45/cb/sb-EE/MW*NMM35*Pa21/SW*Zn23*Zn45/cb
     &      *sb+EE/MW*NMM35*Pa21/SW*Zn25*Zn43/cb/sb-EE/MW
     &      *NMM35*Pa21/SW*Zn25*Zn43/cb*sb+EE/MW*NMM45*Pa21
     &      /SW*Zn24*Zn45/cb*sb+EE/MW*NMM45*Pa21/SW*Zn25
     &      *Zn44/cb*sb
      AAABR(619) = NMM34*Pa12*Sqrt2*Zn23*Zn54/xvev+NMM34
     &      *Pa12*Sqrt2*Zn24*Zn53/xvev+NMM55*Pa12*Sqrt2*Zn25
     &      *Zn55/xvev-Pa12*Sqrt2*Zn25*Zn55*muP/xvev-EE/MW
     &      *NMM13*Pa11/SW*Zn21*Zn53/cb*sb-EE/MW*NMM13*Pa11
     &      /SW*Zn23*Zn51/cb*sb-EE/MW*NMM14*Pa11/SW*Zn21
     &      *Zn54/cb/sb+EE/MW*NMM14*Pa11/SW*Zn21*Zn54/cb
     &      *sb-EE/MW*NMM14*Pa11/SW*Zn24*Zn51/cb/sb+EE/MW
     &      *NMM14*Pa11/SW*Zn24*Zn51/cb*sb-EE/MW*NMM23*Pa11
     &      /SW*Zn22*Zn53/cb*sb-EE/MW*NMM23*Pa11/SW*Zn23
     &      *Zn52/cb*sb-EE/MW*NMM24*Pa11/SW*Zn22*Zn54/cb
     &      /sb+EE/MW*NMM24*Pa11/SW*Zn22*Zn54/cb*sb-EE/MW
     &      *NMM24*Pa11/SW*Zn24*Zn52/cb/sb+EE/MW*NMM24*Pa11
     &      /SW*Zn24*Zn52/cb*sb+EE/MW*NMM35*Pa11/SW*Zn23
     &      *Zn55/cb/sb-EE/MW*NMM35*Pa11/SW*Zn23*Zn55/cb
     &      *sb+EE/MW*NMM35*Pa11/SW*Zn25*Zn53/cb/sb-EE/MW
     &      *NMM35*Pa11/SW*Zn25*Zn53/cb*sb+EE/MW*NMM45*Pa11
     &      /SW*Zn24*Zn55/cb*sb+EE/MW*NMM45*Pa11/SW*Zn25
     &      *Zn54/cb*sb
      AAABR(620) = NMM55*Pa22*Sqrt2*Zn25**2/xvev-Pa22*Sqrt2
     &      *Zn25**2*muP/xvev+2D0*NMM34*Pa22*Sqrt2*Zn23*Zn24
     &      /xvev-2D0*EE/MW*NMM13*Pa21/SW*Zn21*Zn23/cb*sb
     &      -2D0*EE/MW*NMM14*Pa21/SW*Zn21*Zn24/cb/sb+2D0
     &      *EE/MW*NMM14*Pa21/SW*Zn21*Zn24/cb*sb-2D0*EE/MW
     &      *NMM23*Pa21/SW*Zn22*Zn23/cb*sb-2D0*EE/MW*NMM24
     &      *Pa21/SW*Zn22*Zn24/cb/sb+2D0*EE/MW*NMM24*Pa21
     &      /SW*Zn22*Zn24/cb*sb+2D0*EE/MW*NMM35*Pa21/SW*Zn23
     &      *Zn25/cb/sb-2D0*EE/MW*NMM35*Pa21/SW*Zn23*Zn25
     &      /cb*sb+2D0*EE/MW*NMM45*Pa21/SW*Zn24*Zn25/cb*sb
      AAABR(621) = NMM55*Pa12*Sqrt2*Zn25**2/xvev-Pa12*Sqrt2
     &      *Zn25**2*muP/xvev+2D0*NMM34*Pa12*Sqrt2*Zn23*Zn24
     &      /xvev-2D0*EE/MW*NMM13*Pa11/SW*Zn21*Zn23/cb*sb
     &      -2D0*EE/MW*NMM14*Pa11/SW*Zn21*Zn24/cb/sb+2D0
     &      *EE/MW*NMM14*Pa11/SW*Zn21*Zn24/cb*sb-2D0*EE/MW
     &      *NMM23*Pa11/SW*Zn22*Zn23/cb*sb-2D0*EE/MW*NMM24
     &      *Pa11/SW*Zn22*Zn24/cb/sb+2D0*EE/MW*NMM24*Pa11
     &      /SW*Zn22*Zn24/cb*sb+2D0*EE/MW*NMM35*Pa11/SW*Zn23
     &      *Zn25/cb/sb-2D0*EE/MW*NMM35*Pa11/SW*Zn23*Zn25
     &      /cb*sb+2D0*EE/MW*NMM45*Pa11/SW*Zn24*Zn25/cb*sb
      AAABR(622) = NMM34*Pa22*Sqrt2*Zn23*Zn34/xvev+NMM34
     &      *Pa22*Sqrt2*Zn24*Zn33/xvev+NMM55*Pa22*Sqrt2*Zn25
     &      *Zn35/xvev-Pa22*Sqrt2*Zn25*Zn35*muP/xvev-EE/MW
     &      *NMM13*Pa21/SW*Zn21*Zn33/cb*sb-EE/MW*NMM13*Pa21
     &      /SW*Zn23*Zn31/cb*sb-EE/MW*NMM14*Pa21/SW*Zn21
     &      *Zn34/cb/sb+EE/MW*NMM14*Pa21/SW*Zn21*Zn34/cb
     &      *sb-EE/MW*NMM14*Pa21/SW*Zn24*Zn31/cb/sb+EE/MW
     &      *NMM14*Pa21/SW*Zn24*Zn31/cb*sb-EE/MW*NMM23*Pa21
     &      /SW*Zn22*Zn33/cb*sb-EE/MW*NMM23*Pa21/SW*Zn23
     &      *Zn32/cb*sb-EE/MW*NMM24*Pa21/SW*Zn22*Zn34/cb
     &      /sb+EE/MW*NMM24*Pa21/SW*Zn22*Zn34/cb*sb-EE/MW
     &      *NMM24*Pa21/SW*Zn24*Zn32/cb/sb+EE/MW*NMM24*Pa21
     &      /SW*Zn24*Zn32/cb*sb+EE/MW*NMM35*Pa21/SW*Zn23
     &      *Zn35/cb/sb-EE/MW*NMM35*Pa21/SW*Zn23*Zn35/cb
     &      *sb+EE/MW*NMM35*Pa21/SW*Zn25*Zn33/cb/sb-EE/MW
     &      *NMM35*Pa21/SW*Zn25*Zn33/cb*sb+EE/MW*NMM45*Pa21
     &      /SW*Zn24*Zn35/cb*sb+EE/MW*NMM45*Pa21/SW*Zn25
     &      *Zn34/cb*sb
      AAABR(623) = NMM34*Pa22*Sqrt2*Zn23*Zn54/xvev+NMM34
     &      *Pa22*Sqrt2*Zn24*Zn53/xvev+NMM55*Pa22*Sqrt2*Zn25
     &      *Zn55/xvev-Pa22*Sqrt2*Zn25*Zn55*muP/xvev-EE/MW
     &      *NMM13*Pa21/SW*Zn21*Zn53/cb*sb-EE/MW*NMM13*Pa21
     &      /SW*Zn23*Zn51/cb*sb-EE/MW*NMM14*Pa21/SW*Zn21
     &      *Zn54/cb/sb+EE/MW*NMM14*Pa21/SW*Zn21*Zn54/cb
     &      *sb-EE/MW*NMM14*Pa21/SW*Zn24*Zn51/cb/sb+EE/MW
     &      *NMM14*Pa21/SW*Zn24*Zn51/cb*sb-EE/MW*NMM23*Pa21
     &      /SW*Zn22*Zn53/cb*sb-EE/MW*NMM23*Pa21/SW*Zn23
     &      *Zn52/cb*sb-EE/MW*NMM24*Pa21/SW*Zn22*Zn54/cb
     &      /sb+EE/MW*NMM24*Pa21/SW*Zn22*Zn54/cb*sb-EE/MW
     &      *NMM24*Pa21/SW*Zn24*Zn52/cb/sb+EE/MW*NMM24*Pa21
     &      /SW*Zn24*Zn52/cb*sb+EE/MW*NMM35*Pa21/SW*Zn23
     &      *Zn55/cb/sb-EE/MW*NMM35*Pa21/SW*Zn23*Zn55/cb
     &      *sb+EE/MW*NMM35*Pa21/SW*Zn25*Zn53/cb/sb-EE/MW
     &      *NMM35*Pa21/SW*Zn25*Zn53/cb*sb+EE/MW*NMM45*Pa21
     &      /SW*Zn24*Zn55/cb*sb+EE/MW*NMM45*Pa21/SW*Zn25
     &      *Zn54/cb*sb
      AAABR(624) = NMM34*Pa12*Sqrt2*Zn23*Zn44/xvev+NMM34
     &      *Pa12*Sqrt2*Zn24*Zn43/xvev+NMM55*Pa12*Sqrt2*Zn25
     &      *Zn45/xvev-Pa12*Sqrt2*Zn25*Zn45*muP/xvev-EE/MW
     &      *NMM13*Pa11/SW*Zn21*Zn43/cb*sb-EE/MW*NMM13*Pa11
     &      /SW*Zn23*Zn41/cb*sb-EE/MW*NMM14*Pa11/SW*Zn21
     &      *Zn44/cb/sb+EE/MW*NMM14*Pa11/SW*Zn21*Zn44/cb
     &      *sb-EE/MW*NMM14*Pa11/SW*Zn24*Zn41/cb/sb+EE/MW
     &      *NMM14*Pa11/SW*Zn24*Zn41/cb*sb-EE/MW*NMM23*Pa11
     &      /SW*Zn22*Zn43/cb*sb-EE/MW*NMM23*Pa11/SW*Zn23
     &      *Zn42/cb*sb-EE/MW*NMM24*Pa11/SW*Zn22*Zn44/cb
     &      /sb+EE/MW*NMM24*Pa11/SW*Zn22*Zn44/cb*sb-EE/MW
     &      *NMM24*Pa11/SW*Zn24*Zn42/cb/sb+EE/MW*NMM24*Pa11
     &      /SW*Zn24*Zn42/cb*sb+EE/MW*NMM35*Pa11/SW*Zn23
     &      *Zn45/cb/sb-EE/MW*NMM35*Pa11/SW*Zn23*Zn45/cb
     &      *sb+EE/MW*NMM35*Pa11/SW*Zn25*Zn43/cb/sb-EE/MW
     &      *NMM35*Pa11/SW*Zn25*Zn43/cb*sb+EE/MW*NMM45*Pa11
     &      /SW*Zn24*Zn45/cb*sb+EE/MW*NMM45*Pa11/SW*Zn25
     &      *Zn44/cb*sb
      AAABR(625) = NMM34*Sqrt2*Zh33*Zn23*Zn44/xvev+NMM34
     &      *Sqrt2*Zh33*Zn24*Zn43/xvev+NMM55*Sqrt2*Zh33*Zn25
     &      *Zn45/xvev-Sqrt2*Zh33*Zn25*Zn45*muP/xvev+EE/MW
     &      *NMM13/SW*Zh31*Zn21*Zn43/cb+EE/MW*NMM13/SW*Zh31
     &      *Zn23*Zn41/cb+EE/MW*NMM14/SW*Zh32*Zn21*Zn44/sb
     &      +EE/MW*NMM14/SW*Zh32*Zn24*Zn41/sb+EE/MW*NMM23
     &      /SW*Zh31*Zn22*Zn43/cb+EE/MW*NMM23/SW*Zh31*Zn23
     &      *Zn42/cb+EE/MW*NMM24/SW*Zh32*Zn22*Zn44/sb+EE
     &      /MW*NMM24/SW*Zh32*Zn24*Zn42/sb+EE/MW*NMM35/SW
     &      *Zh32*Zn23*Zn45/sb+EE/MW*NMM35/SW*Zh32*Zn25*Zn43
     &      /sb+EE/MW*NMM45/SW*Zh31*Zn24*Zn45/cb+EE/MW*NMM45
     &      /SW*Zh31*Zn25*Zn44/cb
      AAABR(626) = NMM34*Sqrt2*Zh23*Zn23*Zn34/xvev+NMM34
     &      *Sqrt2*Zh23*Zn24*Zn33/xvev+NMM55*Sqrt2*Zh23*Zn25
     &      *Zn35/xvev-Sqrt2*Zh23*Zn25*Zn35*muP/xvev+EE/MW
     &      *NMM13/SW*Zh21*Zn21*Zn33/cb+EE/MW*NMM13/SW*Zh21
     &      *Zn23*Zn31/cb+EE/MW*NMM14/SW*Zh22*Zn21*Zn34/sb
     &      +EE/MW*NMM14/SW*Zh22*Zn24*Zn31/sb+EE/MW*NMM23
     &      /SW*Zh21*Zn22*Zn33/cb+EE/MW*NMM23/SW*Zh21*Zn23
     &      *Zn32/cb+EE/MW*NMM24/SW*Zh22*Zn22*Zn34/sb+EE
     &      /MW*NMM24/SW*Zh22*Zn24*Zn32/sb+EE/MW*NMM35/SW
     &      *Zh22*Zn23*Zn35/sb+EE/MW*NMM35/SW*Zh22*Zn25*Zn33
     &      /sb+EE/MW*NMM45/SW*Zh21*Zn24*Zn35/cb+EE/MW*NMM45
     &      /SW*Zh21*Zn25*Zn34/cb
      AAABR(627) = NMM34*Sqrt2*Zh33*Zn23*Zn54/xvev+NMM34
     &      *Sqrt2*Zh33*Zn24*Zn53/xvev+NMM55*Sqrt2*Zh33*Zn25
     &      *Zn55/xvev-Sqrt2*Zh33*Zn25*Zn55*muP/xvev+EE/MW
     &      *NMM13/SW*Zh31*Zn21*Zn53/cb+EE/MW*NMM13/SW*Zh31
     &      *Zn23*Zn51/cb+EE/MW*NMM14/SW*Zh32*Zn21*Zn54/sb
     &      +EE/MW*NMM14/SW*Zh32*Zn24*Zn51/sb+EE/MW*NMM23
     &      /SW*Zh31*Zn22*Zn53/cb+EE/MW*NMM23/SW*Zh31*Zn23
     &      *Zn52/cb+EE/MW*NMM24/SW*Zh32*Zn22*Zn54/sb+EE
     &      /MW*NMM24/SW*Zh32*Zn24*Zn52/sb+EE/MW*NMM35/SW
     &      *Zh32*Zn23*Zn55/sb+EE/MW*NMM35/SW*Zh32*Zn25*Zn53
     &      /sb+EE/MW*NMM45/SW*Zh31*Zn24*Zn55/cb+EE/MW*NMM45
     &      /SW*Zh31*Zn25*Zn54/cb
      AAABR(628) = NMM34*Sqrt2*Zh23*Zn23*Zn54/xvev+NMM34
     &      *Sqrt2*Zh23*Zn24*Zn53/xvev+NMM55*Sqrt2*Zh23*Zn25
     &      *Zn55/xvev-Sqrt2*Zh23*Zn25*Zn55*muP/xvev+EE/MW
     &      *NMM13/SW*Zh21*Zn21*Zn53/cb+EE/MW*NMM13/SW*Zh21
     &      *Zn23*Zn51/cb+EE/MW*NMM14/SW*Zh22*Zn21*Zn54/sb
     &      +EE/MW*NMM14/SW*Zh22*Zn24*Zn51/sb+EE/MW*NMM23
     &      /SW*Zh21*Zn22*Zn53/cb+EE/MW*NMM23/SW*Zh21*Zn23
     &      *Zn52/cb+EE/MW*NMM24/SW*Zh22*Zn22*Zn54/sb+EE
     &      /MW*NMM24/SW*Zh22*Zn24*Zn52/sb+EE/MW*NMM35/SW
     &      *Zh22*Zn23*Zn55/sb+EE/MW*NMM35/SW*Zh22*Zn25*Zn53
     &      /sb+EE/MW*NMM45/SW*Zh21*Zn24*Zn55/cb+EE/MW*NMM45
     &      /SW*Zh21*Zn25*Zn54/cb
      AAABR(629) = NMM34*Sqrt2*Zh13*Zn23*Zn44/xvev+NMM34
     &      *Sqrt2*Zh13*Zn24*Zn43/xvev+NMM55*Sqrt2*Zh13*Zn25
     &      *Zn45/xvev-Sqrt2*Zh13*Zn25*Zn45*muP/xvev+EE/MW
     &      *NMM13/SW*Zh11*Zn21*Zn43/cb+EE/MW*NMM13/SW*Zh11
     &      *Zn23*Zn41/cb+EE/MW*NMM14/SW*Zh12*Zn21*Zn44/sb
     &      +EE/MW*NMM14/SW*Zh12*Zn24*Zn41/sb+EE/MW*NMM23
     &      /SW*Zh11*Zn22*Zn43/cb+EE/MW*NMM23/SW*Zh11*Zn23
     &      *Zn42/cb+EE/MW*NMM24/SW*Zh12*Zn22*Zn44/sb+EE
     &      /MW*NMM24/SW*Zh12*Zn24*Zn42/sb+EE/MW*NMM35/SW
     &      *Zh12*Zn23*Zn45/sb+EE/MW*NMM35/SW*Zh12*Zn25*Zn43
     &      /sb+EE/MW*NMM45/SW*Zh11*Zn24*Zn45/cb+EE/MW*NMM45
     &      /SW*Zh11*Zn25*Zn44/cb
      AAABR(630) = NMM34*Sqrt2*Zh23*Zn23*Zn44/xvev+NMM34
     &      *Sqrt2*Zh23*Zn24*Zn43/xvev+NMM55*Sqrt2*Zh23*Zn25
     &      *Zn45/xvev-Sqrt2*Zh23*Zn25*Zn45*muP/xvev+EE/MW
     &      *NMM13/SW*Zh21*Zn21*Zn43/cb+EE/MW*NMM13/SW*Zh21
     &      *Zn23*Zn41/cb+EE/MW*NMM14/SW*Zh22*Zn21*Zn44/sb
     &      +EE/MW*NMM14/SW*Zh22*Zn24*Zn41/sb+EE/MW*NMM23
     &      /SW*Zh21*Zn22*Zn43/cb+EE/MW*NMM23/SW*Zh21*Zn23
     &      *Zn42/cb+EE/MW*NMM24/SW*Zh22*Zn22*Zn44/sb+EE
     &      /MW*NMM24/SW*Zh22*Zn24*Zn42/sb+EE/MW*NMM35/SW
     &      *Zh22*Zn23*Zn45/sb+EE/MW*NMM35/SW*Zh22*Zn25*Zn43
     &      /sb+EE/MW*NMM45/SW*Zh21*Zn24*Zn45/cb+EE/MW*NMM45
     &      /SW*Zh21*Zn25*Zn44/cb
      AAABR(631) = NMM34*Sqrt2*Zh13*Zn23*Zn34/xvev+NMM34
     &      *Sqrt2*Zh13*Zn24*Zn33/xvev+NMM55*Sqrt2*Zh13*Zn25
     &      *Zn35/xvev-Sqrt2*Zh13*Zn25*Zn35*muP/xvev+EE/MW
     &      *NMM13/SW*Zh11*Zn21*Zn33/cb+EE/MW*NMM13/SW*Zh11
     &      *Zn23*Zn31/cb+EE/MW*NMM14/SW*Zh12*Zn21*Zn34/sb
     &      +EE/MW*NMM14/SW*Zh12*Zn24*Zn31/sb+EE/MW*NMM23
     &      /SW*Zh11*Zn22*Zn33/cb+EE/MW*NMM23/SW*Zh11*Zn23
     &      *Zn32/cb+EE/MW*NMM24/SW*Zh12*Zn22*Zn34/sb+EE
     &      /MW*NMM24/SW*Zh12*Zn24*Zn32/sb+EE/MW*NMM35/SW
     &      *Zh12*Zn23*Zn35/sb+EE/MW*NMM35/SW*Zh12*Zn25*Zn33
     &      /sb+EE/MW*NMM45/SW*Zh11*Zn24*Zn35/cb+EE/MW*NMM45
     &      /SW*Zh11*Zn25*Zn34/cb
      AAABR(632) = NMM55*Sqrt2*Zh33*Zn25**2/xvev-Sqrt2*Zh33
     &      *Zn25**2*muP/xvev+2D0*NMM34*Sqrt2*Zh33*Zn23*Zn24
     &      /xvev+2D0*EE/MW*NMM13/SW*Zh31*Zn21*Zn23/cb+2D0
     &      *EE/MW*NMM14/SW*Zh32*Zn21*Zn24/sb+2D0*EE/MW*NMM23
     &      /SW*Zh31*Zn22*Zn23/cb+2D0*EE/MW*NMM24/SW*Zh32
     &      *Zn22*Zn24/sb+2D0*EE/MW*NMM35/SW*Zh32*Zn23*Zn25
     &      /sb+2D0*EE/MW*NMM45/SW*Zh31*Zn24*Zn25/cb
      AAABR(633) = NMM55*Sqrt2*Zh23*Zn25**2/xvev-Sqrt2*Zh23
     &      *Zn25**2*muP/xvev+2D0*NMM34*Sqrt2*Zh23*Zn23*Zn24
     &      /xvev+2D0*EE/MW*NMM13/SW*Zh21*Zn21*Zn23/cb+2D0
     &      *EE/MW*NMM14/SW*Zh22*Zn21*Zn24/sb+2D0*EE/MW*NMM23
     &      /SW*Zh21*Zn22*Zn23/cb+2D0*EE/MW*NMM24/SW*Zh22
     &      *Zn22*Zn24/sb+2D0*EE/MW*NMM35/SW*Zh22*Zn23*Zn25
     &      /sb+2D0*EE/MW*NMM45/SW*Zh21*Zn24*Zn25/cb
      AAABR(634) = NMM55*Sqrt2*Zh13*Zn25**2/xvev-Sqrt2*Zh13
     &      *Zn25**2*muP/xvev+2D0*NMM34*Sqrt2*Zh13*Zn23*Zn24
     &      /xvev+2D0*EE/MW*NMM13/SW*Zh11*Zn21*Zn23/cb+2D0
     &      *EE/MW*NMM14/SW*Zh12*Zn21*Zn24/sb+2D0*EE/MW*NMM23
     &      /SW*Zh11*Zn22*Zn23/cb+2D0*EE/MW*NMM24/SW*Zh12
     &      *Zn22*Zn24/sb+2D0*EE/MW*NMM35/SW*Zh12*Zn23*Zn25
     &      /sb+2D0*EE/MW*NMM45/SW*Zh11*Zn24*Zn25/cb
      AAABR(635) = NMM34*Sqrt2*Zh33*Zn23*Zn34/xvev+NMM34
     &      *Sqrt2*Zh33*Zn24*Zn33/xvev+NMM55*Sqrt2*Zh33*Zn25
     &      *Zn35/xvev-Sqrt2*Zh33*Zn25*Zn35*muP/xvev+EE/MW
     &      *NMM13/SW*Zh31*Zn21*Zn33/cb+EE/MW*NMM13/SW*Zh31
     &      *Zn23*Zn31/cb+EE/MW*NMM14/SW*Zh32*Zn21*Zn34/sb
     &      +EE/MW*NMM14/SW*Zh32*Zn24*Zn31/sb+EE/MW*NMM23
     &      /SW*Zh31*Zn22*Zn33/cb+EE/MW*NMM23/SW*Zh31*Zn23
     &      *Zn32/cb+EE/MW*NMM24/SW*Zh32*Zn22*Zn34/sb+EE
     &      /MW*NMM24/SW*Zh32*Zn24*Zn32/sb+EE/MW*NMM35/SW
     &      *Zh32*Zn23*Zn35/sb+EE/MW*NMM35/SW*Zh32*Zn25*Zn33
     &      /sb+EE/MW*NMM45/SW*Zh31*Zn24*Zn35/cb+EE/MW*NMM45
     &      /SW*Zh31*Zn25*Zn34/cb
      AAABR(636) = NMM34*Sqrt2*Zh13*Zn23*Zn54/xvev+NMM34
     &      *Sqrt2*Zh13*Zn24*Zn53/xvev+NMM55*Sqrt2*Zh13*Zn25
     &      *Zn55/xvev-Sqrt2*Zh13*Zn25*Zn55*muP/xvev+EE/MW
     &      *NMM13/SW*Zh11*Zn21*Zn53/cb+EE/MW*NMM13/SW*Zh11
     &      *Zn23*Zn51/cb+EE/MW*NMM14/SW*Zh12*Zn21*Zn54/sb
     &      +EE/MW*NMM14/SW*Zh12*Zn24*Zn51/sb+EE/MW*NMM23
     &      /SW*Zh11*Zn22*Zn53/cb+EE/MW*NMM23/SW*Zh11*Zn23
     &      *Zn52/cb+EE/MW*NMM24/SW*Zh12*Zn22*Zn54/sb+EE
     &      /MW*NMM24/SW*Zh12*Zn24*Zn52/sb+EE/MW*NMM35/SW
     &      *Zh12*Zn23*Zn55/sb+EE/MW*NMM35/SW*Zh12*Zn25*Zn53
     &      /sb+EE/MW*NMM45/SW*Zh11*Zn24*Zn55/cb+EE/MW*NMM45
     &      /SW*Zh11*Zn25*Zn54/cb
      AAABR(637) = EE/MW/SW*(NMM13*Zn31*Zn33-NMM14*Zn31*Zn34
     &      +NMM23*Zn32*Zn33-NMM24*Zn32*Zn34+NMM35*Zn33*Zn35
     &      -NMM45*Zn34*Zn35)
      AAABR(638) = EE/MW/SW*(NMM13*Zn31*Zn43+NMM13*Zn33*Zn41
     &      -NMM14*Zn31*Zn44-NMM14*Zn34*Zn41+NMM23*Zn32*Zn43
     &      +NMM23*Zn33*Zn42-NMM24*Zn32*Zn44-NMM24*Zn34*Zn42
     &      +NMM35*Zn33*Zn45+NMM35*Zn35*Zn43-NMM45*Zn34*Zn45
     &      -NMM45*Zn35*Zn44)
      AAABR(639) = EE/MW/SW*(NMM13*Zn31*Zn53+NMM13*Zn33*Zn51
     &      -NMM14*Zn31*Zn54-NMM14*Zn34*Zn51+NMM23*Zn32*Zn53
     &      +NMM23*Zn33*Zn52-NMM24*Zn32*Zn54-NMM24*Zn34*Zn52
     &      +NMM35*Zn33*Zn55+NMM35*Zn35*Zn53-NMM45*Zn34*Zn55
     &      -NMM45*Zn35*Zn54)
      AAABR(640) = NMM34*Pa12*Sqrt2*Zn33*Zn54/xvev+NMM34
     &      *Pa12*Sqrt2*Zn34*Zn53/xvev+NMM55*Pa12*Sqrt2*Zn35
     &      *Zn55/xvev-Pa12*Sqrt2*Zn35*Zn55*muP/xvev-EE/MW
     &      *NMM13*Pa11/SW*Zn31*Zn53/cb*sb-EE/MW*NMM13*Pa11
     &      /SW*Zn33*Zn51/cb*sb-EE/MW*NMM14*Pa11/SW*Zn31
     &      *Zn54/cb/sb+EE/MW*NMM14*Pa11/SW*Zn31*Zn54/cb
     &      *sb-EE/MW*NMM14*Pa11/SW*Zn34*Zn51/cb/sb+EE/MW
     &      *NMM14*Pa11/SW*Zn34*Zn51/cb*sb-EE/MW*NMM23*Pa11
     &      /SW*Zn32*Zn53/cb*sb-EE/MW*NMM23*Pa11/SW*Zn33
     &      *Zn52/cb*sb-EE/MW*NMM24*Pa11/SW*Zn32*Zn54/cb
     &      /sb+EE/MW*NMM24*Pa11/SW*Zn32*Zn54/cb*sb-EE/MW
     &      *NMM24*Pa11/SW*Zn34*Zn52/cb/sb+EE/MW*NMM24*Pa11
     &      /SW*Zn34*Zn52/cb*sb+EE/MW*NMM35*Pa11/SW*Zn33
     &      *Zn55/cb/sb-EE/MW*NMM35*Pa11/SW*Zn33*Zn55/cb
     &      *sb+EE/MW*NMM35*Pa11/SW*Zn35*Zn53/cb/sb-EE/MW
     &      *NMM35*Pa11/SW*Zn35*Zn53/cb*sb+EE/MW*NMM45*Pa11
     &      /SW*Zn34*Zn55/cb*sb+EE/MW*NMM45*Pa11/SW*Zn35
     &      *Zn54/cb*sb
      AAABR(641) = NMM55*Pa12*Sqrt2*Zn35**2/xvev-Pa12*Sqrt2
     &      *Zn35**2*muP/xvev+2D0*NMM34*Pa12*Sqrt2*Zn33*Zn34
     &      /xvev-2D0*EE/MW*NMM13*Pa11/SW*Zn31*Zn33/cb*sb
     &      -2D0*EE/MW*NMM14*Pa11/SW*Zn31*Zn34/cb/sb+2D0
     &      *EE/MW*NMM14*Pa11/SW*Zn31*Zn34/cb*sb-2D0*EE/MW
     &      *NMM23*Pa11/SW*Zn32*Zn33/cb*sb-2D0*EE/MW*NMM24
     &      *Pa11/SW*Zn32*Zn34/cb/sb+2D0*EE/MW*NMM24*Pa11
     &      /SW*Zn32*Zn34/cb*sb+2D0*EE/MW*NMM35*Pa11/SW*Zn33
     &      *Zn35/cb/sb-2D0*EE/MW*NMM35*Pa11/SW*Zn33*Zn35
     &      /cb*sb+2D0*EE/MW*NMM45*Pa11/SW*Zn34*Zn35/cb*sb
      AAABR(642) = NMM34*Pa12*Sqrt2*Zn33*Zn44/xvev+NMM34
     &      *Pa12*Sqrt2*Zn34*Zn43/xvev+NMM55*Pa12*Sqrt2*Zn35
     &      *Zn45/xvev-Pa12*Sqrt2*Zn35*Zn45*muP/xvev-EE/MW
     &      *NMM13*Pa11/SW*Zn31*Zn43/cb*sb-EE/MW*NMM13*Pa11
     &      /SW*Zn33*Zn41/cb*sb-EE/MW*NMM14*Pa11/SW*Zn31
     &      *Zn44/cb/sb+EE/MW*NMM14*Pa11/SW*Zn31*Zn44/cb
     &      *sb-EE/MW*NMM14*Pa11/SW*Zn34*Zn41/cb/sb+EE/MW
     &      *NMM14*Pa11/SW*Zn34*Zn41/cb*sb-EE/MW*NMM23*Pa11
     &      /SW*Zn32*Zn43/cb*sb-EE/MW*NMM23*Pa11/SW*Zn33
     &      *Zn42/cb*sb-EE/MW*NMM24*Pa11/SW*Zn32*Zn44/cb
     &      /sb+EE/MW*NMM24*Pa11/SW*Zn32*Zn44/cb*sb-EE/MW
     &      *NMM24*Pa11/SW*Zn34*Zn42/cb/sb+EE/MW*NMM24*Pa11
     &      /SW*Zn34*Zn42/cb*sb+EE/MW*NMM35*Pa11/SW*Zn33
     &      *Zn45/cb/sb-EE/MW*NMM35*Pa11/SW*Zn33*Zn45/cb
     &      *sb+EE/MW*NMM35*Pa11/SW*Zn35*Zn43/cb/sb-EE/MW
     &      *NMM35*Pa11/SW*Zn35*Zn43/cb*sb+EE/MW*NMM45*Pa11
     &      /SW*Zn34*Zn45/cb*sb+EE/MW*NMM45*Pa11/SW*Zn35
     &      *Zn44/cb*sb
      AAABR(643) = NMM55*Pa22*Sqrt2*Zn35**2/xvev-Pa22*Sqrt2
     &      *Zn35**2*muP/xvev+2D0*NMM34*Pa22*Sqrt2*Zn33*Zn34
     &      /xvev-2D0*EE/MW*NMM13*Pa21/SW*Zn31*Zn33/cb*sb
     &      -2D0*EE/MW*NMM14*Pa21/SW*Zn31*Zn34/cb/sb+2D0
     &      *EE/MW*NMM14*Pa21/SW*Zn31*Zn34/cb*sb-2D0*EE/MW
     &      *NMM23*Pa21/SW*Zn32*Zn33/cb*sb-2D0*EE/MW*NMM24
     &      *Pa21/SW*Zn32*Zn34/cb/sb+2D0*EE/MW*NMM24*Pa21
     &      /SW*Zn32*Zn34/cb*sb+2D0*EE/MW*NMM35*Pa21/SW*Zn33
     &      *Zn35/cb/sb-2D0*EE/MW*NMM35*Pa21/SW*Zn33*Zn35
     &      /cb*sb+2D0*EE/MW*NMM45*Pa21/SW*Zn34*Zn35/cb*sb
      AAABR(644) = NMM34*Pa22*Sqrt2*Zn33*Zn54/xvev+NMM34
     &      *Pa22*Sqrt2*Zn34*Zn53/xvev+NMM55*Pa22*Sqrt2*Zn35
     &      *Zn55/xvev-Pa22*Sqrt2*Zn35*Zn55*muP/xvev-EE/MW
     &      *NMM13*Pa21/SW*Zn31*Zn53/cb*sb-EE/MW*NMM13*Pa21
     &      /SW*Zn33*Zn51/cb*sb-EE/MW*NMM14*Pa21/SW*Zn31
     &      *Zn54/cb/sb+EE/MW*NMM14*Pa21/SW*Zn31*Zn54/cb
     &      *sb-EE/MW*NMM14*Pa21/SW*Zn34*Zn51/cb/sb+EE/MW
     &      *NMM14*Pa21/SW*Zn34*Zn51/cb*sb-EE/MW*NMM23*Pa21
     &      /SW*Zn32*Zn53/cb*sb-EE/MW*NMM23*Pa21/SW*Zn33
     &      *Zn52/cb*sb-EE/MW*NMM24*Pa21/SW*Zn32*Zn54/cb
     &      /sb+EE/MW*NMM24*Pa21/SW*Zn32*Zn54/cb*sb-EE/MW
     &      *NMM24*Pa21/SW*Zn34*Zn52/cb/sb+EE/MW*NMM24*Pa21
     &      /SW*Zn34*Zn52/cb*sb+EE/MW*NMM35*Pa21/SW*Zn33
     &      *Zn55/cb/sb-EE/MW*NMM35*Pa21/SW*Zn33*Zn55/cb
     &      *sb+EE/MW*NMM35*Pa21/SW*Zn35*Zn53/cb/sb-EE/MW
     &      *NMM35*Pa21/SW*Zn35*Zn53/cb*sb+EE/MW*NMM45*Pa21
     &      /SW*Zn34*Zn55/cb*sb+EE/MW*NMM45*Pa21/SW*Zn35
     &      *Zn54/cb*sb
      AAABR(645) = NMM34*Pa22*Sqrt2*Zn33*Zn44/xvev+NMM34
     &      *Pa22*Sqrt2*Zn34*Zn43/xvev+NMM55*Pa22*Sqrt2*Zn35
     &      *Zn45/xvev-Pa22*Sqrt2*Zn35*Zn45*muP/xvev-EE/MW
     &      *NMM13*Pa21/SW*Zn31*Zn43/cb*sb-EE/MW*NMM13*Pa21
     &      /SW*Zn33*Zn41/cb*sb-EE/MW*NMM14*Pa21/SW*Zn31
     &      *Zn44/cb/sb+EE/MW*NMM14*Pa21/SW*Zn31*Zn44/cb
     &      *sb-EE/MW*NMM14*Pa21/SW*Zn34*Zn41/cb/sb+EE/MW
     &      *NMM14*Pa21/SW*Zn34*Zn41/cb*sb-EE/MW*NMM23*Pa21
     &      /SW*Zn32*Zn43/cb*sb-EE/MW*NMM23*Pa21/SW*Zn33
     &      *Zn42/cb*sb-EE/MW*NMM24*Pa21/SW*Zn32*Zn44/cb
     &      /sb+EE/MW*NMM24*Pa21/SW*Zn32*Zn44/cb*sb-EE/MW
     &      *NMM24*Pa21/SW*Zn34*Zn42/cb/sb+EE/MW*NMM24*Pa21
     &      /SW*Zn34*Zn42/cb*sb+EE/MW*NMM35*Pa21/SW*Zn33
     &      *Zn45/cb/sb-EE/MW*NMM35*Pa21/SW*Zn33*Zn45/cb
     &      *sb+EE/MW*NMM35*Pa21/SW*Zn35*Zn43/cb/sb-EE/MW
     &      *NMM35*Pa21/SW*Zn35*Zn43/cb*sb+EE/MW*NMM45*Pa21
     &      /SW*Zn34*Zn45/cb*sb+EE/MW*NMM45*Pa21/SW*Zn35
     &      *Zn44/cb*sb
      AAABR(646) = NMM55*Sqrt2*Zh23*Zn35**2/xvev-Sqrt2*Zh23
     &      *Zn35**2*muP/xvev+2D0*NMM34*Sqrt2*Zh23*Zn33*Zn34
     &      /xvev+2D0*EE/MW*NMM13/SW*Zh21*Zn31*Zn33/cb+2D0
     &      *EE/MW*NMM14/SW*Zh22*Zn31*Zn34/sb+2D0*EE/MW*NMM23
     &      /SW*Zh21*Zn32*Zn33/cb+2D0*EE/MW*NMM24/SW*Zh22
     &      *Zn32*Zn34/sb+2D0*EE/MW*NMM35/SW*Zh22*Zn33*Zn35
     &      /sb+2D0*EE/MW*NMM45/SW*Zh21*Zn34*Zn35/cb
      AAABR(647) = NMM34*Sqrt2*Zh23*Zn33*Zn54/xvev+NMM34
     &      *Sqrt2*Zh23*Zn34*Zn53/xvev+NMM55*Sqrt2*Zh23*Zn35
     &      *Zn55/xvev-Sqrt2*Zh23*Zn35*Zn55*muP/xvev+EE/MW
     &      *NMM13/SW*Zh21*Zn31*Zn53/cb+EE/MW*NMM13/SW*Zh21
     &      *Zn33*Zn51/cb+EE/MW*NMM14/SW*Zh22*Zn31*Zn54/sb
     &      +EE/MW*NMM14/SW*Zh22*Zn34*Zn51/sb+EE/MW*NMM23
     &      /SW*Zh21*Zn32*Zn53/cb+EE/MW*NMM23/SW*Zh21*Zn33
     &      *Zn52/cb+EE/MW*NMM24/SW*Zh22*Zn32*Zn54/sb+EE
     &      /MW*NMM24/SW*Zh22*Zn34*Zn52/sb+EE/MW*NMM35/SW
     &      *Zh22*Zn33*Zn55/sb+EE/MW*NMM35/SW*Zh22*Zn35*Zn53
     &      /sb+EE/MW*NMM45/SW*Zh21*Zn34*Zn55/cb+EE/MW*NMM45
     &      /SW*Zh21*Zn35*Zn54/cb
      AAABR(648) = NMM55*Sqrt2*Zh13*Zn35**2/xvev-Sqrt2*Zh13
     &      *Zn35**2*muP/xvev+2D0*NMM34*Sqrt2*Zh13*Zn33*Zn34
     &      /xvev+2D0*EE/MW*NMM13/SW*Zh11*Zn31*Zn33/cb+2D0
     &      *EE/MW*NMM14/SW*Zh12*Zn31*Zn34/sb+2D0*EE/MW*NMM23
     &      /SW*Zh11*Zn32*Zn33/cb+2D0*EE/MW*NMM24/SW*Zh12
     &      *Zn32*Zn34/sb+2D0*EE/MW*NMM35/SW*Zh12*Zn33*Zn35
     &      /sb+2D0*EE/MW*NMM45/SW*Zh11*Zn34*Zn35/cb
      AAABR(649) = NMM34*Sqrt2*Zh13*Zn33*Zn54/xvev+NMM34
     &      *Sqrt2*Zh13*Zn34*Zn53/xvev+NMM55*Sqrt2*Zh13*Zn35
     &      *Zn55/xvev-Sqrt2*Zh13*Zn35*Zn55*muP/xvev+EE/MW
     &      *NMM13/SW*Zh11*Zn31*Zn53/cb+EE/MW*NMM13/SW*Zh11
     &      *Zn33*Zn51/cb+EE/MW*NMM14/SW*Zh12*Zn31*Zn54/sb
     &      +EE/MW*NMM14/SW*Zh12*Zn34*Zn51/sb+EE/MW*NMM23
     &      /SW*Zh11*Zn32*Zn53/cb+EE/MW*NMM23/SW*Zh11*Zn33
     &      *Zn52/cb+EE/MW*NMM24/SW*Zh12*Zn32*Zn54/sb+EE
     &      /MW*NMM24/SW*Zh12*Zn34*Zn52/sb+EE/MW*NMM35/SW
     &      *Zh12*Zn33*Zn55/sb+EE/MW*NMM35/SW*Zh12*Zn35*Zn53
     &      /sb+EE/MW*NMM45/SW*Zh11*Zn34*Zn55/cb+EE/MW*NMM45
     &      /SW*Zh11*Zn35*Zn54/cb
      AAABR(650) = NMM34*Sqrt2*Zh33*Zn33*Zn44/xvev+NMM34
     &      *Sqrt2*Zh33*Zn34*Zn43/xvev+NMM55*Sqrt2*Zh33*Zn35
     &      *Zn45/xvev-Sqrt2*Zh33*Zn35*Zn45*muP/xvev+EE/MW
     &      *NMM13/SW*Zh31*Zn31*Zn43/cb+EE/MW*NMM13/SW*Zh31
     &      *Zn33*Zn41/cb+EE/MW*NMM14/SW*Zh32*Zn31*Zn44/sb
     &      +EE/MW*NMM14/SW*Zh32*Zn34*Zn41/sb+EE/MW*NMM23
     &      /SW*Zh31*Zn32*Zn43/cb+EE/MW*NMM23/SW*Zh31*Zn33
     &      *Zn42/cb+EE/MW*NMM24/SW*Zh32*Zn32*Zn44/sb+EE
     &      /MW*NMM24/SW*Zh32*Zn34*Zn42/sb+EE/MW*NMM35/SW
     &      *Zh32*Zn33*Zn45/sb+EE/MW*NMM35/SW*Zh32*Zn35*Zn43
     &      /sb+EE/MW*NMM45/SW*Zh31*Zn34*Zn45/cb+EE/MW*NMM45
     &      /SW*Zh31*Zn35*Zn44/cb
      AAABR(651) = NMM34*Sqrt2*Zh33*Zn33*Zn54/xvev+NMM34
     &      *Sqrt2*Zh33*Zn34*Zn53/xvev+NMM55*Sqrt2*Zh33*Zn35
     &      *Zn55/xvev-Sqrt2*Zh33*Zn35*Zn55*muP/xvev+EE/MW
     &      *NMM13/SW*Zh31*Zn31*Zn53/cb+EE/MW*NMM13/SW*Zh31
     &      *Zn33*Zn51/cb+EE/MW*NMM14/SW*Zh32*Zn31*Zn54/sb
     &      +EE/MW*NMM14/SW*Zh32*Zn34*Zn51/sb+EE/MW*NMM23
     &      /SW*Zh31*Zn32*Zn53/cb+EE/MW*NMM23/SW*Zh31*Zn33
     &      *Zn52/cb+EE/MW*NMM24/SW*Zh32*Zn32*Zn54/sb+EE
     &      /MW*NMM24/SW*Zh32*Zn34*Zn52/sb+EE/MW*NMM35/SW
     &      *Zh32*Zn33*Zn55/sb+EE/MW*NMM35/SW*Zh32*Zn35*Zn53
     &      /sb+EE/MW*NMM45/SW*Zh31*Zn34*Zn55/cb+EE/MW*NMM45
     &      /SW*Zh31*Zn35*Zn54/cb
      AAABR(652) = NMM34*Sqrt2*Zh13*Zn33*Zn44/xvev+NMM34
     &      *Sqrt2*Zh13*Zn34*Zn43/xvev+NMM55*Sqrt2*Zh13*Zn35
     &      *Zn45/xvev-Sqrt2*Zh13*Zn35*Zn45*muP/xvev+EE/MW
     &      *NMM13/SW*Zh11*Zn31*Zn43/cb+EE/MW*NMM13/SW*Zh11
     &      *Zn33*Zn41/cb+EE/MW*NMM14/SW*Zh12*Zn31*Zn44/sb
     &      +EE/MW*NMM14/SW*Zh12*Zn34*Zn41/sb+EE/MW*NMM23
     &      /SW*Zh11*Zn32*Zn43/cb+EE/MW*NMM23/SW*Zh11*Zn33
     &      *Zn42/cb+EE/MW*NMM24/SW*Zh12*Zn32*Zn44/sb+EE
     &      /MW*NMM24/SW*Zh12*Zn34*Zn42/sb+EE/MW*NMM35/SW
     &      *Zh12*Zn33*Zn45/sb+EE/MW*NMM35/SW*Zh12*Zn35*Zn43
     &      /sb+EE/MW*NMM45/SW*Zh11*Zn34*Zn45/cb+EE/MW*NMM45
     &      /SW*Zh11*Zn35*Zn44/cb
      AAABR(653) = NMM55*Sqrt2*Zh33*Zn35**2/xvev-Sqrt2*Zh33
     &      *Zn35**2*muP/xvev+2D0*NMM34*Sqrt2*Zh33*Zn33*Zn34
     &      /xvev+2D0*EE/MW*NMM13/SW*Zh31*Zn31*Zn33/cb+2D0
     &      *EE/MW*NMM14/SW*Zh32*Zn31*Zn34/sb+2D0*EE/MW*NMM23
     &      /SW*Zh31*Zn32*Zn33/cb+2D0*EE/MW*NMM24/SW*Zh32
     &      *Zn32*Zn34/sb+2D0*EE/MW*NMM35/SW*Zh32*Zn33*Zn35
     &      /sb+2D0*EE/MW*NMM45/SW*Zh31*Zn34*Zn35/cb
      AAABR(654) = NMM34*Sqrt2*Zh23*Zn33*Zn44/xvev+NMM34
     &      *Sqrt2*Zh23*Zn34*Zn43/xvev+NMM55*Sqrt2*Zh23*Zn35
     &      *Zn45/xvev-Sqrt2*Zh23*Zn35*Zn45*muP/xvev+EE/MW
     &      *NMM13/SW*Zh21*Zn31*Zn43/cb+EE/MW*NMM13/SW*Zh21
     &      *Zn33*Zn41/cb+EE/MW*NMM14/SW*Zh22*Zn31*Zn44/sb
     &      +EE/MW*NMM14/SW*Zh22*Zn34*Zn41/sb+EE/MW*NMM23
     &      /SW*Zh21*Zn32*Zn43/cb+EE/MW*NMM23/SW*Zh21*Zn33
     &      *Zn42/cb+EE/MW*NMM24/SW*Zh22*Zn32*Zn44/sb+EE
     &      /MW*NMM24/SW*Zh22*Zn34*Zn42/sb+EE/MW*NMM35/SW
     &      *Zh22*Zn33*Zn45/sb+EE/MW*NMM35/SW*Zh22*Zn35*Zn43
     &      /sb+EE/MW*NMM45/SW*Zh21*Zn34*Zn45/cb+EE/MW*NMM45
     &      /SW*Zh21*Zn35*Zn44/cb
      AAABR(655) = EE/MW/SW*(NMM13*Zn41*Zn53+NMM13*Zn43*Zn51
     &      -NMM14*Zn41*Zn54-NMM14*Zn44*Zn51+NMM23*Zn42*Zn53
     &      +NMM23*Zn43*Zn52-NMM24*Zn42*Zn54-NMM24*Zn44*Zn52
     &      +NMM35*Zn43*Zn55+NMM35*Zn45*Zn53-NMM45*Zn44*Zn55
     &      -NMM45*Zn45*Zn54)
      AAABR(656) = EE/MW/SW*(NMM13*Zn41*Zn43-NMM14*Zn41*Zn44
     &      +NMM23*Zn42*Zn43-NMM24*Zn42*Zn44+NMM35*Zn43*Zn45
     &      -NMM45*Zn44*Zn45)
      AAABR(657) = NMM34*Pa22*Sqrt2*Zn43*Zn54/xvev+NMM34
     &      *Pa22*Sqrt2*Zn44*Zn53/xvev+NMM55*Pa22*Sqrt2*Zn45
     &      *Zn55/xvev-Pa22*Sqrt2*Zn45*Zn55*muP/xvev-EE/MW
     &      *NMM13*Pa21/SW*Zn41*Zn53/cb*sb-EE/MW*NMM13*Pa21
     &      /SW*Zn43*Zn51/cb*sb-EE/MW*NMM14*Pa21/SW*Zn41
     &      *Zn54/cb/sb+EE/MW*NMM14*Pa21/SW*Zn41*Zn54/cb
     &      *sb-EE/MW*NMM14*Pa21/SW*Zn44*Zn51/cb/sb+EE/MW
     &      *NMM14*Pa21/SW*Zn44*Zn51/cb*sb-EE/MW*NMM23*Pa21
     &      /SW*Zn42*Zn53/cb*sb-EE/MW*NMM23*Pa21/SW*Zn43
     &      *Zn52/cb*sb-EE/MW*NMM24*Pa21/SW*Zn42*Zn54/cb
     &      /sb+EE/MW*NMM24*Pa21/SW*Zn42*Zn54/cb*sb-EE/MW
     &      *NMM24*Pa21/SW*Zn44*Zn52/cb/sb+EE/MW*NMM24*Pa21
     &      /SW*Zn44*Zn52/cb*sb+EE/MW*NMM35*Pa21/SW*Zn43
     &      *Zn55/cb/sb-EE/MW*NMM35*Pa21/SW*Zn43*Zn55/cb
     &      *sb+EE/MW*NMM35*Pa21/SW*Zn45*Zn53/cb/sb-EE/MW
     &      *NMM35*Pa21/SW*Zn45*Zn53/cb*sb+EE/MW*NMM45*Pa21
     &      /SW*Zn44*Zn55/cb*sb+EE/MW*NMM45*Pa21/SW*Zn45
     &      *Zn54/cb*sb
      AAABR(658) = NMM34*Pa12*Sqrt2*Zn43*Zn54/xvev+NMM34
     &      *Pa12*Sqrt2*Zn44*Zn53/xvev+NMM55*Pa12*Sqrt2*Zn45
     &      *Zn55/xvev-Pa12*Sqrt2*Zn45*Zn55*muP/xvev-EE/MW
     &      *NMM13*Pa11/SW*Zn41*Zn53/cb*sb-EE/MW*NMM13*Pa11
     &      /SW*Zn43*Zn51/cb*sb-EE/MW*NMM14*Pa11/SW*Zn41
     &      *Zn54/cb/sb+EE/MW*NMM14*Pa11/SW*Zn41*Zn54/cb
     &      *sb-EE/MW*NMM14*Pa11/SW*Zn44*Zn51/cb/sb+EE/MW
     &      *NMM14*Pa11/SW*Zn44*Zn51/cb*sb-EE/MW*NMM23*Pa11
     &      /SW*Zn42*Zn53/cb*sb-EE/MW*NMM23*Pa11/SW*Zn43
     &      *Zn52/cb*sb-EE/MW*NMM24*Pa11/SW*Zn42*Zn54/cb
     &      /sb+EE/MW*NMM24*Pa11/SW*Zn42*Zn54/cb*sb-EE/MW
     &      *NMM24*Pa11/SW*Zn44*Zn52/cb/sb+EE/MW*NMM24*Pa11
     &      /SW*Zn44*Zn52/cb*sb+EE/MW*NMM35*Pa11/SW*Zn43
     &      *Zn55/cb/sb-EE/MW*NMM35*Pa11/SW*Zn43*Zn55/cb
     &      *sb+EE/MW*NMM35*Pa11/SW*Zn45*Zn53/cb/sb-EE/MW
     &      *NMM35*Pa11/SW*Zn45*Zn53/cb*sb+EE/MW*NMM45*Pa11
     &      /SW*Zn44*Zn55/cb*sb+EE/MW*NMM45*Pa11/SW*Zn45
     &      *Zn54/cb*sb
      AAABR(659) = NMM55*Pa22*Sqrt2*Zn45**2/xvev-Pa22*Sqrt2
     &      *Zn45**2*muP/xvev+2D0*NMM34*Pa22*Sqrt2*Zn43*Zn44
     &      /xvev-2D0*EE/MW*NMM13*Pa21/SW*Zn41*Zn43/cb*sb
     &      -2D0*EE/MW*NMM14*Pa21/SW*Zn41*Zn44/cb/sb+2D0
     &      *EE/MW*NMM14*Pa21/SW*Zn41*Zn44/cb*sb-2D0*EE/MW
     &      *NMM23*Pa21/SW*Zn42*Zn43/cb*sb-2D0*EE/MW*NMM24
     &      *Pa21/SW*Zn42*Zn44/cb/sb+2D0*EE/MW*NMM24*Pa21
     &      /SW*Zn42*Zn44/cb*sb+2D0*EE/MW*NMM35*Pa21/SW*Zn43
     &      *Zn45/cb/sb-2D0*EE/MW*NMM35*Pa21/SW*Zn43*Zn45
     &      /cb*sb+2D0*EE/MW*NMM45*Pa21/SW*Zn44*Zn45/cb*sb
      AAABR(660) = NMM55*Pa12*Sqrt2*Zn45**2/xvev-Pa12*Sqrt2
     &      *Zn45**2*muP/xvev+2D0*NMM34*Pa12*Sqrt2*Zn43*Zn44
     &      /xvev-2D0*EE/MW*NMM13*Pa11/SW*Zn41*Zn43/cb*sb
     &      -2D0*EE/MW*NMM14*Pa11/SW*Zn41*Zn44/cb/sb+2D0
     &      *EE/MW*NMM14*Pa11/SW*Zn41*Zn44/cb*sb-2D0*EE/MW
     &      *NMM23*Pa11/SW*Zn42*Zn43/cb*sb-2D0*EE/MW*NMM24
     &      *Pa11/SW*Zn42*Zn44/cb/sb+2D0*EE/MW*NMM24*Pa11
     &      /SW*Zn42*Zn44/cb*sb+2D0*EE/MW*NMM35*Pa11/SW*Zn43
     &      *Zn45/cb/sb-2D0*EE/MW*NMM35*Pa11/SW*Zn43*Zn45
     &      /cb*sb+2D0*EE/MW*NMM45*Pa11/SW*Zn44*Zn45/cb*sb
      AAABR(661) = NMM55*Sqrt2*Zh23*Zn45**2/xvev-Sqrt2*Zh23
     &      *Zn45**2*muP/xvev+2D0*NMM34*Sqrt2*Zh23*Zn43*Zn44
     &      /xvev+2D0*EE/MW*NMM13/SW*Zh21*Zn41*Zn43/cb+2D0
     &      *EE/MW*NMM14/SW*Zh22*Zn41*Zn44/sb+2D0*EE/MW*NMM23
     &      /SW*Zh21*Zn42*Zn43/cb+2D0*EE/MW*NMM24/SW*Zh22
     &      *Zn42*Zn44/sb+2D0*EE/MW*NMM35/SW*Zh22*Zn43*Zn45
     &      /sb+2D0*EE/MW*NMM45/SW*Zh21*Zn44*Zn45/cb
      AAABR(662) = NMM34*Sqrt2*Zh33*Zn43*Zn54/xvev+NMM34
     &      *Sqrt2*Zh33*Zn44*Zn53/xvev+NMM55*Sqrt2*Zh33*Zn45
     &      *Zn55/xvev-Sqrt2*Zh33*Zn45*Zn55*muP/xvev+EE/MW
     &      *NMM13/SW*Zh31*Zn41*Zn53/cb+EE/MW*NMM13/SW*Zh31
     &      *Zn43*Zn51/cb+EE/MW*NMM14/SW*Zh32*Zn41*Zn54/sb
     &      +EE/MW*NMM14/SW*Zh32*Zn44*Zn51/sb+EE/MW*NMM23
     &      /SW*Zh31*Zn42*Zn53/cb+EE/MW*NMM23/SW*Zh31*Zn43
     &      *Zn52/cb+EE/MW*NMM24/SW*Zh32*Zn42*Zn54/sb+EE
     &      /MW*NMM24/SW*Zh32*Zn44*Zn52/sb+EE/MW*NMM35/SW
     &      *Zh32*Zn43*Zn55/sb+EE/MW*NMM35/SW*Zh32*Zn45*Zn53
     &      /sb+EE/MW*NMM45/SW*Zh31*Zn44*Zn55/cb+EE/MW*NMM45
     &      /SW*Zh31*Zn45*Zn54/cb
      AAABR(663) = NMM34*Sqrt2*Zh23*Zn43*Zn54/xvev+NMM34
     &      *Sqrt2*Zh23*Zn44*Zn53/xvev+NMM55*Sqrt2*Zh23*Zn45
     &      *Zn55/xvev-Sqrt2*Zh23*Zn45*Zn55*muP/xvev+EE/MW
     &      *NMM13/SW*Zh21*Zn41*Zn53/cb+EE/MW*NMM13/SW*Zh21
     &      *Zn43*Zn51/cb+EE/MW*NMM14/SW*Zh22*Zn41*Zn54/sb
     &      +EE/MW*NMM14/SW*Zh22*Zn44*Zn51/sb+EE/MW*NMM23
     &      /SW*Zh21*Zn42*Zn53/cb+EE/MW*NMM23/SW*Zh21*Zn43
     &      *Zn52/cb+EE/MW*NMM24/SW*Zh22*Zn42*Zn54/sb+EE
     &      /MW*NMM24/SW*Zh22*Zn44*Zn52/sb+EE/MW*NMM35/SW
     &      *Zh22*Zn43*Zn55/sb+EE/MW*NMM35/SW*Zh22*Zn45*Zn53
     &      /sb+EE/MW*NMM45/SW*Zh21*Zn44*Zn55/cb+EE/MW*NMM45
     &      /SW*Zh21*Zn45*Zn54/cb
      AAABR(664) = NMM34*Sqrt2*Zh13*Zn43*Zn54/xvev+NMM34
     &      *Sqrt2*Zh13*Zn44*Zn53/xvev+NMM55*Sqrt2*Zh13*Zn45
     &      *Zn55/xvev-Sqrt2*Zh13*Zn45*Zn55*muP/xvev+EE/MW
     &      *NMM13/SW*Zh11*Zn41*Zn53/cb+EE/MW*NMM13/SW*Zh11
     &      *Zn43*Zn51/cb+EE/MW*NMM14/SW*Zh12*Zn41*Zn54/sb
     &      +EE/MW*NMM14/SW*Zh12*Zn44*Zn51/sb+EE/MW*NMM23
     &      /SW*Zh11*Zn42*Zn53/cb+EE/MW*NMM23/SW*Zh11*Zn43
     &      *Zn52/cb+EE/MW*NMM24/SW*Zh12*Zn42*Zn54/sb+EE
     &      /MW*NMM24/SW*Zh12*Zn44*Zn52/sb+EE/MW*NMM35/SW
     &      *Zh12*Zn43*Zn55/sb+EE/MW*NMM35/SW*Zh12*Zn45*Zn53
     &      /sb+EE/MW*NMM45/SW*Zh11*Zn44*Zn55/cb+EE/MW*NMM45
     &      /SW*Zh11*Zn45*Zn54/cb
      AAABR(665) = NMM55*Sqrt2*Zh13*Zn45**2/xvev-Sqrt2*Zh13
     &      *Zn45**2*muP/xvev+2D0*NMM34*Sqrt2*Zh13*Zn43*Zn44
     &      /xvev+2D0*EE/MW*NMM13/SW*Zh11*Zn41*Zn43/cb+2D0
     &      *EE/MW*NMM14/SW*Zh12*Zn41*Zn44/sb+2D0*EE/MW*NMM23
     &      /SW*Zh11*Zn42*Zn43/cb+2D0*EE/MW*NMM24/SW*Zh12
     &      *Zn42*Zn44/sb+2D0*EE/MW*NMM35/SW*Zh12*Zn43*Zn45
     &      /sb+2D0*EE/MW*NMM45/SW*Zh11*Zn44*Zn45/cb
      AAABR(666) = NMM55*Sqrt2*Zh33*Zn45**2/xvev-Sqrt2*Zh33
     &      *Zn45**2*muP/xvev+2D0*NMM34*Sqrt2*Zh33*Zn43*Zn44
     &      /xvev+2D0*EE/MW*NMM13/SW*Zh31*Zn41*Zn43/cb+2D0
     &      *EE/MW*NMM14/SW*Zh32*Zn41*Zn44/sb+2D0*EE/MW*NMM23
     &      /SW*Zh31*Zn42*Zn43/cb+2D0*EE/MW*NMM24/SW*Zh32
     &      *Zn42*Zn44/sb+2D0*EE/MW*NMM35/SW*Zh32*Zn43*Zn45
     &      /sb+2D0*EE/MW*NMM45/SW*Zh31*Zn44*Zn45/cb
      AAABR(667) = EE/MW/SW*(NMM13*Zn51*Zn53-NMM14*Zn51*Zn54
     &      +NMM23*Zn52*Zn53-NMM24*Zn52*Zn54+NMM35*Zn53*Zn55
     &      -NMM45*Zn54*Zn55)
      AAABR(668) = NMM55*Pa12*Sqrt2*Zn55**2/xvev-Pa12*Sqrt2
     &      *Zn55**2*muP/xvev+2D0*NMM34*Pa12*Sqrt2*Zn53*Zn54
     &      /xvev-2D0*EE/MW*NMM13*Pa11/SW*Zn51*Zn53/cb*sb
     &      -2D0*EE/MW*NMM14*Pa11/SW*Zn51*Zn54/cb/sb+2D0
     &      *EE/MW*NMM14*Pa11/SW*Zn51*Zn54/cb*sb-2D0*EE/MW
     &      *NMM23*Pa11/SW*Zn52*Zn53/cb*sb-2D0*EE/MW*NMM24
     &      *Pa11/SW*Zn52*Zn54/cb/sb+2D0*EE/MW*NMM24*Pa11
     &      /SW*Zn52*Zn54/cb*sb+2D0*EE/MW*NMM35*Pa11/SW*Zn53
     &      *Zn55/cb/sb-2D0*EE/MW*NMM35*Pa11/SW*Zn53*Zn55
     &      /cb*sb+2D0*EE/MW*NMM45*Pa11/SW*Zn54*Zn55/cb*sb
      AAABR(669) = NMM55*Pa22*Sqrt2*Zn55**2/xvev-Pa22*Sqrt2
     &      *Zn55**2*muP/xvev+2D0*NMM34*Pa22*Sqrt2*Zn53*Zn54
     &      /xvev-2D0*EE/MW*NMM13*Pa21/SW*Zn51*Zn53/cb*sb
     &      -2D0*EE/MW*NMM14*Pa21/SW*Zn51*Zn54/cb/sb+2D0
     &      *EE/MW*NMM14*Pa21/SW*Zn51*Zn54/cb*sb-2D0*EE/MW
     &      *NMM23*Pa21/SW*Zn52*Zn53/cb*sb-2D0*EE/MW*NMM24
     &      *Pa21/SW*Zn52*Zn54/cb/sb+2D0*EE/MW*NMM24*Pa21
     &      /SW*Zn52*Zn54/cb*sb+2D0*EE/MW*NMM35*Pa21/SW*Zn53
     &      *Zn55/cb/sb-2D0*EE/MW*NMM35*Pa21/SW*Zn53*Zn55
     &      /cb*sb+2D0*EE/MW*NMM45*Pa21/SW*Zn54*Zn55/cb*sb
      AAABR(670) = NMM55*Sqrt2*Zh33*Zn55**2/xvev-Sqrt2*Zh33
     &      *Zn55**2*muP/xvev+2D0*NMM34*Sqrt2*Zh33*Zn53*Zn54
     &      /xvev+2D0*EE/MW*NMM13/SW*Zh31*Zn51*Zn53/cb+2D0
     &      *EE/MW*NMM14/SW*Zh32*Zn51*Zn54/sb+2D0*EE/MW*NMM23
     &      /SW*Zh31*Zn52*Zn53/cb+2D0*EE/MW*NMM24/SW*Zh32
     &      *Zn52*Zn54/sb+2D0*EE/MW*NMM35/SW*Zh32*Zn53*Zn55
     &      /sb+2D0*EE/MW*NMM45/SW*Zh31*Zn54*Zn55/cb
      AAABR(671) = NMM55*Sqrt2*Zh13*Zn55**2/xvev-Sqrt2*Zh13
     &      *Zn55**2*muP/xvev+2D0*NMM34*Sqrt2*Zh13*Zn53*Zn54
     &      /xvev+2D0*EE/MW*NMM13/SW*Zh11*Zn51*Zn53/cb+2D0
     &      *EE/MW*NMM14/SW*Zh12*Zn51*Zn54/sb+2D0*EE/MW*NMM23
     &      /SW*Zh11*Zn52*Zn53/cb+2D0*EE/MW*NMM24/SW*Zh12
     &      *Zn52*Zn54/sb+2D0*EE/MW*NMM35/SW*Zh12*Zn53*Zn55
     &      /sb+2D0*EE/MW*NMM45/SW*Zh11*Zn54*Zn55/cb
      AAABR(672) = NMM55*Sqrt2*Zh23*Zn55**2/xvev-Sqrt2*Zh23
     &      *Zn55**2*muP/xvev+2D0*NMM34*Sqrt2*Zh23*Zn53*Zn54
     &      /xvev+2D0*EE/MW*NMM13/SW*Zh21*Zn51*Zn53/cb+2D0
     &      *EE/MW*NMM14/SW*Zh22*Zn51*Zn54/sb+2D0*EE/MW*NMM23
     &      /SW*Zh21*Zn52*Zn53/cb+2D0*EE/MW*NMM24/SW*Zh22
     &      *Zn52*Zn54/sb+2D0*EE/MW*NMM35/SW*Zh22*Zn53*Zn55
     &      /sb+2D0*EE/MW*NMM45/SW*Zh21*Zn54*Zn55/cb
      AAABR(673) = EE*(SW/CW*Zu12*Zu22-2D0*CW/SW*Zu11*Zu21
     &      -CW/SW*Zu12*Zu22)
      AAABR(674) = EE*(SW/CW*Zv12*Zv22-2D0*CW/SW*Zv11*Zv21
     &      -CW/SW*Zv12*Zv22)
      AAABR(675) = EE*(SW/CW*Zu12**2-2D0*CW/SW*Zu11**2-CW
     &      /SW*Zu12**2)
      AAABR(676) = EE*(SW/CW*Zv12**2-2D0*CW/SW*Zv11**2-CW
     &      /SW*Zv12**2)
      AAABR(677) = EE*(SW/CW*Zu22**2-2D0*CW/SW*Zu21**2-CW
     &      /SW*Zu22**2)
      AAABR(678) = EE*(SW/CW*Zv22**2-2D0*CW/SW*Zv21**2-CW
     &      /SW*Zv22**2)
      AAABR(679) = EE/SW*(2D0*Zn42*Zu21+Sqrt2*Zn43*Zu22)
      AAABR(680) = EE/SW*(2D0*Zn42*Zv21-Sqrt2*Zn44*Zv22)
      AAABR(681) = EE/SW*(2D0*Zn52*Zu11+Sqrt2*Zn53*Zu12)
      AAABR(682) = EE/SW*(2D0*Zn52*Zv11-Sqrt2*Zn54*Zv12)
      AAABR(683) = EE/SW*(2D0*Zn32*Zu21+Sqrt2*Zn33*Zu22)
      AAABR(684) = EE/SW*(2D0*Zn32*Zv21-Sqrt2*Zn34*Zv22)
      AAABR(685) = EE/SW*(2D0*Zn12*Zu21+Sqrt2*Zn13*Zu22)
      AAABR(686) = EE/SW*(2D0*Zn12*Zv21-Sqrt2*Zn14*Zv22)
      AAABR(687) = EE/SW*(2D0*Zn22*Zu11+Sqrt2*Zn23*Zu12)
      AAABR(688) = EE/SW*(2D0*Zn22*Zv11-Sqrt2*Zn24*Zv12)
      AAABR(689) = EE/SW*(2D0*Zn22*Zu21+Sqrt2*Zn23*Zu22)
      AAABR(690) = EE/SW*(2D0*Zn22*Zv21-Sqrt2*Zn24*Zv22)
      AAABR(691) = EE/SW*(2D0*Zn42*Zu11+Sqrt2*Zn43*Zu12)
      AAABR(692) = EE/SW*(2D0*Zn42*Zv11-Sqrt2*Zn44*Zv12)
      AAABR(693) = EE/SW*(2D0*Zn52*Zu21+Sqrt2*Zn53*Zu22)
      AAABR(694) = EE/SW*(2D0*Zn52*Zv21-Sqrt2*Zn54*Zv22)
      AAABR(695) = EE/SW*(2D0*Zn32*Zu11+Sqrt2*Zn33*Zu12)
      AAABR(696) = EE/SW*(2D0*Zn32*Zv11-Sqrt2*Zn34*Zv12)
      AAABR(697) = EE/SW*(2D0*Zn12*Zu11+Sqrt2*Zn13*Zu12)
      AAABR(698) = EE/SW*(2D0*Zn12*Zv11-Sqrt2*Zn14*Zv12)
      AAABR(699) = EE*(SW/CW+3D0*CW/SW)
      AAABR(700) = EE*(SW/CW-3D0*CW/SW)
      AAABR(701) = EE*(SW/CW*Zn13*Zn23-SW/CW*Zn14*Zn24+CW
     &      /SW*Zn13*Zn23-CW/SW*Zn14*Zn24)
      AAABR(702) = EE*(SW/CW*Zn13*Zn33-SW/CW*Zn14*Zn34+CW
     &      /SW*Zn13*Zn33-CW/SW*Zn14*Zn34)
      AAABR(703) = EE*(SW/CW*Zn13*Zn53-SW/CW*Zn14*Zn54+CW
     &      /SW*Zn13*Zn53-CW/SW*Zn14*Zn54)
      AAABR(704) = EE*(SW/CW*Zn13*Zn43-SW/CW*Zn14*Zn44+CW
     &      /SW*Zn13*Zn43-CW/SW*Zn14*Zn44)
      AAABR(705) = EE*(SW/CW*Zn13**2-SW/CW*Zn14**2+CW/SW
     &      *Zn13**2-CW/SW*Zn14**2)
      AAABR(706) = EE*(SW/CW*Zn23*Zn53-SW/CW*Zn24*Zn54+CW
     &      /SW*Zn23*Zn53-CW/SW*Zn24*Zn54)
      AAABR(707) = EE*(SW/CW*Zn23**2-SW/CW*Zn24**2+CW/SW
     &      *Zn23**2-CW/SW*Zn24**2)
      AAABR(708) = EE*(SW/CW*Zn23*Zn43-SW/CW*Zn24*Zn44+CW
     &      /SW*Zn23*Zn43-CW/SW*Zn24*Zn44)
      AAABR(709) = EE*(SW/CW*Zn23*Zn33-SW/CW*Zn24*Zn34+CW
     &      /SW*Zn23*Zn33-CW/SW*Zn24*Zn34)
      AAABR(710) = EE*(SW/CW*Zn33*Zn43-SW/CW*Zn34*Zn44+CW
     &      /SW*Zn33*Zn43-CW/SW*Zn34*Zn44)
      AAABR(711) = EE*(SW/CW*Zn33*Zn53-SW/CW*Zn34*Zn54+CW
     &      /SW*Zn33*Zn53-CW/SW*Zn34*Zn54)
      AAABR(712) = EE*(SW/CW*Zn33**2-SW/CW*Zn34**2+CW/SW
     &      *Zn33**2-CW/SW*Zn34**2)
      AAABR(713) = EE*(SW/CW*Zn43**2-SW/CW*Zn44**2+CW/SW
     &      *Zn43**2-CW/SW*Zn44**2)
      AAABR(714) = EE*(SW/CW*Zn43*Zn53-SW/CW*Zn44*Zn54+CW
     &      /SW*Zn43*Zn53-CW/SW*Zn44*Zn54)
      AAABR(715) = EE*(SW/CW*Zn53**2-SW/CW*Zn54**2+CW/SW
     &      *Zn53**2-CW/SW*Zn54**2)
      AAABR(716) = la2+la1*sb**4-2D0*la2*sb**2+la2*sb**4
     &      +2D0*la3*sb**2-2D0*la3*sb**4+2D0*la4*sb**2-2D0
     &      *la4*sb**4+2D0*la5*sb**2-2D0*la5*sb**4+4D0*cb
     &      *la6*sb**3+4D0*cb*la7*sb-4D0*cb*la7*sb**3
      AAABR(717) = la7+3D0*la6*sb**2-4D0*la6*sb**4-5D0*la7
     &      *sb**2+4D0*la7*sb**4+cb*la1*sb**3-cb*la2*sb+cb
     &      *la2*sb**3+cb*la3*sb-2D0*cb*la3*sb**3+cb*la4
     &      *sb-2D0*cb*la4*sb**3+cb*la5*sb-2D0*cb*la5*sb**3
      AAABR(718) = la5+la1*sb**2-la1*sb**4+la2*sb**2-la2
     &      *sb**4-2D0*la3*sb**2+2D0*la3*sb**4-2D0*la4*sb
     &      **2+2D0*la4*sb**4-2D0*la5*sb**2+2D0*la5*sb**4
     &      +2D0*cb*la6*sb-4D0*cb*la6*sb**3-2D0*cb*la7*sb
     &      +4D0*cb*la7*sb**3
      AAABR(719) = EE**2*(1D0/CW**2+1D0/SW**2-2D0*sb**2/CW
     &      **2-2D0*sb**2/SW**2)
      AAABR(720) = EE**2/cb**2*(Zl11**2/CW**2-2D0*Zl12**2
     &      /CW**2+Zl11**2/SW**2-3D0*Zl11**2/CW**2*sb**2
     &      +2D0*Zl11**2/CW**2*sb**4+6D0*Zl12**2/CW**2*sb
     &      **2-4D0*Zl12**2/CW**2*sb**4-3D0*Zl11**2/SW**2
     &      *sb**2+2D0*Zl11**2/SW**2*sb**4-2D0*Ml**2/MW**2
     &      /SW**2*Zl12**2*sb**2)
      AAABR(721) = EE**2/cb**2*(Zl11/CW**2*Zl21-2D0*Zl12
     &      /CW**2*Zl22+Zl11/SW**2*Zl21-3D0*Zl11/CW**2*Zl21
     &      *sb**2+2D0*Zl11/CW**2*Zl21*sb**4+6D0*Zl12/CW
     &      **2*Zl22*sb**2-4D0*Zl12/CW**2*Zl22*sb**4-3D0
     &      *Zl11/SW**2*Zl21*sb**2+2D0*Zl11/SW**2*Zl21*sb
     &      **4-2D0*Ml**2/MW**2/SW**2*Zl12*Zl22*sb**2)
      AAABR(722) = EE**2/CW**2*(1D0-2D0*sb**2)
      AAABR(723) = EE**2/cb**2*(Zl21**2/CW**2-2D0*Zl22**2
     &      /CW**2+Zl21**2/SW**2-3D0*Zl21**2/CW**2*sb**2
     &      +2D0*Zl21**2/CW**2*sb**4+6D0*Zl22**2/CW**2*sb
     &      **2-4D0*Zl22**2/CW**2*sb**4-3D0*Zl21**2/SW**2
     &      *sb**2+2D0*Zl21**2/SW**2*sb**4-2D0*Ml**2/MW**2
     &      /SW**2*Zl22**2*sb**2)
      AAABR(724) = EE**2/cb**2*(1D0/CW**2-1D0/SW**2-3D0*sb
     &      **2/CW**2+2D0*sb**4/CW**2+3D0*sb**2/SW**2-2D0
     &      *sb**4/SW**2-2D0*Ml**2/MW**2/SW**2*sb**2)
      AAABR(725) = EE**2*(1D0/CW**2-1D0/SW**2-2D0*sb**2/CW
     &      **2+2D0*sb**2/SW**2)
      AAABR(726) = EE**2/cb**2*(Zss11**2/CW**2+2D0*Zss12
     &      **2/CW**2-3D0*Zss11**2/SW**2-3D0*Zss11**2/CW
     &      **2*sb**2+2D0*Zss11**2/CW**2*sb**4-6D0*Zss12
     &      **2/CW**2*sb**2+4D0*Zss12**2/CW**2*sb**4+9D0
     &      *Zss11**2/SW**2*sb**2-6D0*Zss11**2/SW**2*sb**4
     &      -12D0*Mc**2/MW**2/SW**2*Zss11**2+6D0*Mc**2/MW
     &      **2/SW**2*Zss11**2/sb**2+6D0*Mc**2/MW**2/SW**2
     &      *Zss11**2*sb**2+6D0*Mq**2/MW**2/SW**2*Zss12**2
     &      *sb**2)
      AAABR(727) = EE**2/cb**2*(Zdd11**2/CW**2+2D0*Zdd12
     &      **2/CW**2-3D0*Zdd11**2/SW**2-3D0*Zdd11**2/CW
     &      **2*sb**2+2D0*Zdd11**2/CW**2*sb**4-6D0*Zdd12
     &      **2/CW**2*sb**2+4D0*Zdd12**2/CW**2*sb**4+9D0
     &      *Zdd11**2/SW**2*sb**2-6D0*Zdd11**2/SW**2*sb**4
     &      -12D0*Mq**2/MW**2/SW**2*Zdd11**2+6D0*Mq**2/MW
     &      **2/SW**2*sb**2+6D0*Mq**2/MW**2/SW**2*Zdd11**2
     &      /sb**2)
      AAABR(728) = EE**2/cb**2*(Zb11**2/CW**2+2D0*Zb12**2
     &      /CW**2-3D0*Zb11**2/SW**2-3D0*Zb11**2/CW**2*sb
     &      **2+2D0*Zb11**2/CW**2*sb**4-6D0*Zb12**2/CW**2
     &      *sb**2+4D0*Zb12**2/CW**2*sb**4+9D0*Zb11**2/SW
     &      **2*sb**2-6D0*Zb11**2/SW**2*sb**4-12D0*MtMM**2
     &      /MW**2/SW**2*Zb11**2+6D0*MbMM**2/MW**2/SW**2
     &      *Zb12**2*sb**2+6D0*MtMM**2/MW**2/SW**2*Zb11**2
     &      /sb**2+6D0*MtMM**2/MW**2/SW**2*Zb11**2*sb**2)
      AAABR(729) = EE**2/cb**2*(Zss11/CW**2*Zss21+2D0*Zss12
     &      /CW**2*Zss22-3D0*Zss11/SW**2*Zss21-3D0*Zss11
     &      /CW**2*Zss21*sb**2+2D0*Zss11/CW**2*Zss21*sb**4
     &      -6D0*Zss12/CW**2*Zss22*sb**2+4D0*Zss12/CW**2
     &      *Zss22*sb**4+9D0*Zss11/SW**2*Zss21*sb**2-6D0
     &      *Zss11/SW**2*Zss21*sb**4-12D0*Mc**2/MW**2/SW
     &      **2*Zss11*Zss21+6D0*Mc**2/MW**2/SW**2*Zss11*Zss21
     &      /sb**2+6D0*Mc**2/MW**2/SW**2*Zss11*Zss21*sb**2
     &      +6D0*Mq**2/MW**2/SW**2*Zss12*Zss22*sb**2)
      AAABR(730) = EE**2/cb**2*(Zdd11/CW**2*Zdd21+2D0*Zdd12
     &      /CW**2*Zdd22-3D0*Zdd11/SW**2*Zdd21-3D0*Zdd11
     &      /CW**2*Zdd21*sb**2+2D0*Zdd11/CW**2*Zdd21*sb**4
     &      -6D0*Zdd12/CW**2*Zdd22*sb**2+4D0*Zdd12/CW**2
     &      *Zdd22*sb**4+9D0*Zdd11/SW**2*Zdd21*sb**2-6D0
     &      *Zdd11/SW**2*Zdd21*sb**4-12D0*Mq**2/MW**2/SW
     &      **2*Zdd11*Zdd21+6D0*Mq**2/MW**2/SW**2*Zdd11*Zdd21
     &      /sb**2)
      AAABR(731) = EE**2/cb**2*(Zb11/CW**2*Zb21+2D0*Zb12
     &      /CW**2*Zb22-3D0*Zb11/SW**2*Zb21-3D0*Zb11/CW**2
     &      *Zb21*sb**2+2D0*Zb11/CW**2*Zb21*sb**4-6D0*Zb12
     &      /CW**2*Zb22*sb**2+4D0*Zb12/CW**2*Zb22*sb**4+9D0
     &      *Zb11/SW**2*Zb21*sb**2-6D0*Zb11/SW**2*Zb21*sb
     &      **4-12D0*MtMM**2/MW**2/SW**2*Zb11*Zb21+6D0*MbMM
     &      **2/MW**2/SW**2*Zb12*Zb22*sb**2+6D0*MtMM**2/MW
     &      **2/SW**2*Zb11*Zb21/sb**2+6D0*MtMM**2/MW**2/SW
     &      **2*Zb11*Zb21*sb**2)
      AAABR(732) = EE**2/cb**2*(Zb21**2/CW**2+2D0*Zb22**2
     &      /CW**2-3D0*Zb21**2/SW**2-3D0*Zb21**2/CW**2*sb
     &      **2+2D0*Zb21**2/CW**2*sb**4-6D0*Zb22**2/CW**2
     &      *sb**2+4D0*Zb22**2/CW**2*sb**4+9D0*Zb21**2/SW
     &      **2*sb**2-6D0*Zb21**2/SW**2*sb**4-12D0*MtMM**2
     &      /MW**2/SW**2*Zb21**2+6D0*MbMM**2/MW**2/SW**2
     &      *Zb22**2*sb**2+6D0*MtMM**2/MW**2/SW**2*Zb21**2
     &      /sb**2+6D0*MtMM**2/MW**2/SW**2*Zb21**2*sb**2)
      AAABR(733) = EE**2/cb**2*(Zdd21**2/CW**2+2D0*Zdd22
     &      **2/CW**2-3D0*Zdd21**2/SW**2-3D0*Zdd21**2/CW
     &      **2*sb**2+2D0*Zdd21**2/CW**2*sb**4-6D0*Zdd22
     &      **2/CW**2*sb**2+4D0*Zdd22**2/CW**2*sb**4+9D0
     &      *Zdd21**2/SW**2*sb**2-6D0*Zdd21**2/SW**2*sb**4
     &      -12D0*Mq**2/MW**2/SW**2*Zdd21**2+6D0*Mq**2/MW
     &      **2/SW**2*sb**2+6D0*Mq**2/MW**2/SW**2*Zdd21**2
     &      /sb**2)
      AAABR(734) = EE**2/cb**2*(Zss21**2/CW**2+2D0*Zss22
     &      **2/CW**2-3D0*Zss21**2/SW**2-3D0*Zss21**2/CW
     &      **2*sb**2+2D0*Zss21**2/CW**2*sb**4-6D0*Zss22
     &      **2/CW**2*sb**2+4D0*Zss22**2/CW**2*sb**4+9D0
     &      *Zss21**2/SW**2*sb**2-6D0*Zss21**2/SW**2*sb**4
     &      -12D0*Mc**2/MW**2/SW**2*Zss21**2+6D0*Mc**2/MW
     &      **2/SW**2*Zss21**2/sb**2+6D0*Mc**2/MW**2/SW**2
     &      *Zss21**2*sb**2+6D0*Mq**2/MW**2/SW**2*Zss22**2
     &      *sb**2)
      AAABR(735) = EE**2/cb**2*(Zcc11**2/CW**2-4D0*Zcc12
     &      **2/CW**2+3D0*Zcc11**2/SW**2-3D0*Zcc11**2/CW
     &      **2*sb**2+2D0*Zcc11**2/CW**2*sb**4+12D0*Zcc12
     &      **2/CW**2*sb**2-8D0*Zcc12**2/CW**2*sb**4-9D0
     &      *Zcc11**2/SW**2*sb**2+6D0*Zcc11**2/SW**2*sb**4
     &      -12D0*Mc**2/MW**2/SW**2*Zcc12**2+6D0*Mc**2/MW
     &      **2/SW**2*Zcc12**2/sb**2+6D0*Mc**2/MW**2/SW**2
     &      *Zcc12**2*sb**2+6D0*Mq**2/MW**2/SW**2*Zcc11**2
     &      *sb**2)
      AAABR(736) = EE**2/cb**2*(Zt11**2/CW**2-4D0*Zt12**2
     &      /CW**2+3D0*Zt11**2/SW**2-3D0*Zt11**2/CW**2*sb
     &      **2+2D0*Zt11**2/CW**2*sb**4+12D0*Zt12**2/CW**2
     &      *sb**2-8D0*Zt12**2/CW**2*sb**4-9D0*Zt11**2/SW
     &      **2*sb**2+6D0*Zt11**2/SW**2*sb**4-12D0*MtMM**2
     &      /MW**2/SW**2*Zt12**2+6D0*MbMM**2/MW**2/SW**2
     &      *Zt11**2*sb**2+6D0*MtMM**2/MW**2/SW**2*Zt12**2
     &      /sb**2+6D0*MtMM**2/MW**2/SW**2*Zt12**2*sb**2)
      AAABR(737) = EE**2/cb**2*(Zuu11**2/CW**2-4D0*Zuu12
     &      **2/CW**2+3D0*Zuu11**2/SW**2-3D0*Zuu11**2/CW
     &      **2*sb**2+2D0*Zuu11**2/CW**2*sb**4+12D0*Zuu12
     &      **2/CW**2*sb**2-8D0*Zuu12**2/CW**2*sb**4-9D0
     &      *Zuu11**2/SW**2*sb**2+6D0*Zuu11**2/SW**2*sb**4
     &      -12D0*Mq**2/MW**2/SW**2*Zuu12**2+6D0*Mq**2/MW
     &      **2/SW**2*sb**2+6D0*Mq**2/MW**2/SW**2*Zuu12**2
     &      /sb**2)
      AAABR(738) = EE**2/cb**2*(Zuu11/CW**2*Zuu21-4D0*Zuu12
     &      /CW**2*Zuu22+3D0*Zuu11/SW**2*Zuu21-3D0*Zuu11
     &      /CW**2*Zuu21*sb**2+2D0*Zuu11/CW**2*Zuu21*sb**4
     &      +12D0*Zuu12/CW**2*Zuu22*sb**2-8D0*Zuu12/CW**2
     &      *Zuu22*sb**4-9D0*Zuu11/SW**2*Zuu21*sb**2+6D0
     &      *Zuu11/SW**2*Zuu21*sb**4-12D0*Mq**2/MW**2/SW
     &      **2*Zuu12*Zuu22+6D0*Mq**2/MW**2/SW**2*Zuu12*Zuu22
     &      /sb**2)
      AAABR(739) = EE**2/cb**2*(Zt11/CW**2*Zt21-4D0*Zt12
     &      /CW**2*Zt22+3D0*Zt11/SW**2*Zt21-3D0*Zt11/CW**2
     &      *Zt21*sb**2+2D0*Zt11/CW**2*Zt21*sb**4+12D0*Zt12
     &      /CW**2*Zt22*sb**2-8D0*Zt12/CW**2*Zt22*sb**4-9D0
     &      *Zt11/SW**2*Zt21*sb**2+6D0*Zt11/SW**2*Zt21*sb
     &      **4-12D0*MtMM**2/MW**2/SW**2*Zt12*Zt22+6D0*MbMM
     &      **2/MW**2/SW**2*Zt11*Zt21*sb**2+6D0*MtMM**2/MW
     &      **2/SW**2*Zt12*Zt22/sb**2+6D0*MtMM**2/MW**2/SW
     &      **2*Zt12*Zt22*sb**2)
      AAABR(740) = EE**2/cb**2*(Zcc11/CW**2*Zcc21-4D0*Zcc12
     &      /CW**2*Zcc22+3D0*Zcc11/SW**2*Zcc21-3D0*Zcc11
     &      /CW**2*Zcc21*sb**2+2D0*Zcc11/CW**2*Zcc21*sb**4
     &      +12D0*Zcc12/CW**2*Zcc22*sb**2-8D0*Zcc12/CW**2
     &      *Zcc22*sb**4-9D0*Zcc11/SW**2*Zcc21*sb**2+6D0
     &      *Zcc11/SW**2*Zcc21*sb**4-12D0*Mc**2/MW**2/SW
     &      **2*Zcc12*Zcc22+6D0*Mc**2/MW**2/SW**2*Zcc12*Zcc22
     &      /sb**2+6D0*Mc**2/MW**2/SW**2*Zcc12*Zcc22*sb**2
     &      +6D0*Mq**2/MW**2/SW**2*Zcc11*Zcc21*sb**2)
      AAABR(741) = EE**2/cb**2*(Zcc21**2/CW**2-4D0*Zcc22
     &      **2/CW**2+3D0*Zcc21**2/SW**2-3D0*Zcc21**2/CW
     &      **2*sb**2+2D0*Zcc21**2/CW**2*sb**4+12D0*Zcc22
     &      **2/CW**2*sb**2-8D0*Zcc22**2/CW**2*sb**4-9D0
     &      *Zcc21**2/SW**2*sb**2+6D0*Zcc21**2/SW**2*sb**4
     &      -12D0*Mc**2/MW**2/SW**2*Zcc22**2+6D0*Mc**2/MW
     &      **2/SW**2*Zcc22**2/sb**2+6D0*Mc**2/MW**2/SW**2
     &      *Zcc22**2*sb**2+6D0*Mq**2/MW**2/SW**2*Zcc21**2
     &      *sb**2)
      AAABR(742) = EE**2/cb**2*(Zt21**2/CW**2-4D0*Zt22**2
     &      /CW**2+3D0*Zt21**2/SW**2-3D0*Zt21**2/CW**2*sb
     &      **2+2D0*Zt21**2/CW**2*sb**4+12D0*Zt22**2/CW**2
     &      *sb**2-8D0*Zt22**2/CW**2*sb**4-9D0*Zt21**2/SW
     &      **2*sb**2+6D0*Zt21**2/SW**2*sb**4-12D0*MtMM**2
     &      /MW**2/SW**2*Zt22**2+6D0*MbMM**2/MW**2/SW**2
     &      *Zt21**2*sb**2+6D0*MtMM**2/MW**2/SW**2*Zt22**2
     &      /sb**2+6D0*MtMM**2/MW**2/SW**2*Zt22**2*sb**2)
      AAABR(743) = EE**2/cb**2*(Zuu21**2/CW**2-4D0*Zuu22
     &      **2/CW**2+3D0*Zuu21**2/SW**2-3D0*Zuu21**2/CW
     &      **2*sb**2+2D0*Zuu21**2/CW**2*sb**4+12D0*Zuu22
     &      **2/CW**2*sb**2-8D0*Zuu22**2/CW**2*sb**4-9D0
     &      *Zuu21**2/SW**2*sb**2+6D0*Zuu21**2/SW**2*sb**4
     &      -12D0*Mq**2/MW**2/SW**2*Zuu22**2+6D0*Mq**2/MW
     &      **2/SW**2*sb**2+6D0*Mq**2/MW**2/SW**2*Zuu22**2
     &      /sb**2)
      AAABR(744) = la3+la4+2D0*la1*sb**2-2D0*la1*sb**4+2D0
     &      *la2*sb**2-2D0*la2*sb**4-4D0*la3*sb**2+4D0*la3
     &      *sb**4-4D0*la4*sb**2+4D0*la4*sb**4-4D0*la5*sb
     &      **2+4D0*la5*sb**4+4D0*cb*la6*sb-8D0*cb*la6*sb
     &      **3-4D0*cb*la7*sb+8D0*cb*la7*sb**3
      AAABR(745) = la3+la1*sb**2-la1*sb**4+la2*sb**2-la2
     &      *sb**4-2D0*la3*sb**2+2D0*la3*sb**4-2D0*la4*sb
     &      **2+2D0*la4*sb**4-2D0*la5*sb**2+2D0*la5*sb**4
     &      +2D0*cb*la6*sb-4D0*cb*la6*sb**3-2D0*cb*la7*sb
     &      +4D0*cb*la7*sb**3
      AAABR(746) = Pa21*(la7+3D0*la6*sb**2-4D0*la6*sb**4
     &      -5D0*la7*sb**2+4D0*la7*sb**4+cb*la1*sb**3-cb
     &      *la2*sb+cb*la2*sb**3+cb*la3*sb-2D0*cb*la3*sb
     &      **3+cb*la4*sb-2D0*cb*la4*sb**3+cb*la5*sb-2D0
     &      *cb*la5*sb**3)
      AAABR(747) = Pa11*(la7+3D0*la6*sb**2-4D0*la6*sb**4
     &      -5D0*la7*sb**2+4D0*la7*sb**4+cb*la1*sb**3-cb
     &      *la2*sb+cb*la2*sb**3+cb*la3*sb-2D0*cb*la3*sb
     &      **3+cb*la4*sb-2D0*cb*la4*sb**3+cb*la5*sb-2D0
     &      *cb*la5*sb**3)
      AAABR(748) = Pa11*Pa21*la2+Pa12*Pa22*la2s+Pa11*Pa21
     &      *la1*sb**4-2D0*Pa11*Pa21*la2*sb**2+Pa11*Pa21
     &      *la2*sb**4+2D0*Pa11*Pa21*la3*sb**2-2D0*Pa11*Pa21
     &      *la3*sb**4+2D0*Pa11*Pa21*la4*sb**2-2D0*Pa11*Pa21
     &      *la4*sb**4+2D0*Pa11*Pa21*la5*sb**2-2D0*Pa11*Pa21
     &      *la5*sb**4+Pa12*Pa22*la1s*sb**2-Pa12*Pa22*la2s
     &      *sb**2+4D0*Pa11*Pa21*cb*la6*sb**3+4D0*Pa11*Pa21
     &      *cb*la7*sb-4D0*Pa11*Pa21*cb*la7*sb**3+2D0*Pa12
     &      *Pa22*cb*la4s*sb-2D0*Pa12*Pa22*cb*la5s*sb-2D0
     &      *Pa12*Pa22*cb*la6s*sb
      AAABR(749) = Pa21**2*la2+Pa22**2*la2s+Pa21**2*la1*sb
     &      **4-2D0*Pa21**2*la2*sb**2+Pa21**2*la2*sb**4+2D0
     &      *Pa21**2*la3*sb**2-2D0*Pa21**2*la3*sb**4+2D0
     &      *Pa21**2*la4*sb**2-2D0*Pa21**2*la4*sb**4+2D0
     &      *Pa21**2*la5*sb**2-2D0*Pa21**2*la5*sb**4+Pa22
     &      **2*la1s*sb**2-Pa22**2*la2s*sb**2+4D0*Pa21**2
     &      *cb*la6*sb**3+4D0*Pa21**2*cb*la7*sb-4D0*Pa21
     &      **2*cb*la7*sb**3+2D0*Pa22**2*cb*la4s*sb-2D0*Pa22
     &      **2*cb*la5s*sb-2D0*Pa22**2*cb*la6s*sb
      AAABR(750) = Pa11**2*la2+Pa12**2*la2s+Pa11**2*la1*sb
     &      **4-2D0*Pa11**2*la2*sb**2+Pa11**2*la2*sb**4+2D0
     &      *Pa11**2*la3*sb**2-2D0*Pa11**2*la3*sb**4+2D0
     &      *Pa11**2*la4*sb**2-2D0*Pa11**2*la4*sb**4+2D0
     &      *Pa11**2*la5*sb**2-2D0*Pa11**2*la5*sb**4+Pa12
     &      **2*la1s*sb**2-Pa12**2*la2s*sb**2+4D0*Pa11**2
     &      *cb*la6*sb**3+4D0*Pa11**2*cb*la7*sb-4D0*Pa11
     &      **2*cb*la7*sb**3+2D0*Pa12**2*cb*la4s*sb-2D0*Pa12
     &      **2*cb*la5s*sb-2D0*Pa12**2*cb*la6s*sb
      AAABR(751) = Zh21*Zh31*la3-Zh21*Zh32*la7-Zh22*Zh31
     &      *la7+Zh22*Zh32*la2+Zh23*Zh33*la2s+Zh21*Zh31*la1
     &      *sb**2-Zh21*Zh31*la3*sb**2-Zh21*Zh32*la6*sb**2
     &      +Zh21*Zh32*la7*sb**2-Zh22*Zh31*la6*sb**2+Zh22
     &      *Zh31*la7*sb**2-Zh22*Zh32*la2*sb**2+Zh22*Zh32
     &      *la3*sb**2+Zh23*Zh33*la1s*sb**2-Zh23*Zh33*la2s
     &      *sb**2+2D0*Zh21*Zh31*cb*la6*sb-Zh21*Zh32*cb*la4
     &      *sb-Zh21*Zh32*cb*la5*sb-Zh22*Zh31*cb*la4*sb-Zh22
     &      *Zh31*cb*la5*sb+2D0*Zh22*Zh32*cb*la7*sb+2D0*Zh23
     &      *Zh33*cb*la4s*sb+2D0*Zh23*Zh33*cb*la5s*sb+2D0
     &      *Zh23*Zh33*cb*la6s*sb
      AAABR(752) = Zh11*Zh31*la3-Zh11*Zh32*la7-Zh12*Zh31
     &      *la7+Zh12*Zh32*la2+Zh13*Zh33*la2s+Zh11*Zh31*la1
     &      *sb**2-Zh11*Zh31*la3*sb**2-Zh11*Zh32*la6*sb**2
     &      +Zh11*Zh32*la7*sb**2-Zh12*Zh31*la6*sb**2+Zh12
     &      *Zh31*la7*sb**2-Zh12*Zh32*la2*sb**2+Zh12*Zh32
     &      *la3*sb**2+Zh13*Zh33*la1s*sb**2-Zh13*Zh33*la2s
     &      *sb**2+2D0*Zh11*Zh31*cb*la6*sb-Zh11*Zh32*cb*la4
     &      *sb-Zh11*Zh32*cb*la5*sb-Zh12*Zh31*cb*la4*sb-Zh12
     &      *Zh31*cb*la5*sb+2D0*Zh12*Zh32*cb*la7*sb+2D0*Zh13
     &      *Zh33*cb*la4s*sb+2D0*Zh13*Zh33*cb*la5s*sb+2D0
     &      *Zh13*Zh33*cb*la6s*sb
      AAABR(753) = Zh11*Zh21*la3-Zh11*Zh22*la7-Zh12*Zh21
     &      *la7+Zh12*Zh22*la2+Zh13*Zh23*la2s+Zh11*Zh21*la1
     &      *sb**2-Zh11*Zh21*la3*sb**2-Zh11*Zh22*la6*sb**2
     &      +Zh11*Zh22*la7*sb**2-Zh12*Zh21*la6*sb**2+Zh12
     &      *Zh21*la7*sb**2-Zh12*Zh22*la2*sb**2+Zh12*Zh22
     &      *la3*sb**2+Zh13*Zh23*la1s*sb**2-Zh13*Zh23*la2s
     &      *sb**2+2D0*Zh11*Zh21*cb*la6*sb-Zh11*Zh22*cb*la4
     &      *sb-Zh11*Zh22*cb*la5*sb-Zh12*Zh21*cb*la4*sb-Zh12
     &      *Zh21*cb*la5*sb+2D0*Zh12*Zh22*cb*la7*sb+2D0*Zh13
     &      *Zh23*cb*la4s*sb+2D0*Zh13*Zh23*cb*la5s*sb+2D0
     &      *Zh13*Zh23*cb*la6s*sb
      AAABR(754) = Zh11**2*la3+Zh12**2*la2+Zh13**2*la2s-2D0
     &      *Zh11*Zh12*la7+Zh11**2*la1*sb**2-Zh11**2*la3
     &      *sb**2-Zh12**2*la2*sb**2+Zh12**2*la3*sb**2+Zh13
     &      **2*la1s*sb**2-Zh13**2*la2s*sb**2-2D0*Zh11*Zh12
     &      *la6*sb**2+2D0*Zh11*Zh12*la7*sb**2+2D0*Zh11**2
     &      *cb*la6*sb+2D0*Zh12**2*cb*la7*sb+2D0*Zh13**2
     &      *cb*la4s*sb+2D0*Zh13**2*cb*la5s*sb+2D0*Zh13**2
     &      *cb*la6s*sb-2D0*Zh11*Zh12*cb*la4*sb-2D0*Zh11
     &      *Zh12*cb*la5*sb
      AAABR(755) = Zh31**2*la3+Zh32**2*la2+Zh33**2*la2s-2D0
     &      *Zh31*Zh32*la7+Zh31**2*la1*sb**2-Zh31**2*la3
     &      *sb**2-Zh32**2*la2*sb**2+Zh32**2*la3*sb**2+Zh33
     &      **2*la1s*sb**2-Zh33**2*la2s*sb**2-2D0*Zh31*Zh32
     &      *la6*sb**2+2D0*Zh31*Zh32*la7*sb**2+2D0*Zh31**2
     &      *cb*la6*sb+2D0*Zh32**2*cb*la7*sb+2D0*Zh33**2
     &      *cb*la4s*sb+2D0*Zh33**2*cb*la5s*sb+2D0*Zh33**2
     &      *cb*la6s*sb-2D0*Zh31*Zh32*cb*la4*sb-2D0*Zh31
     &      *Zh32*cb*la5*sb
      AAABR(756) = Zh21**2*la3+Zh22**2*la2+Zh23**2*la2s-2D0
     &      *Zh21*Zh22*la7+Zh21**2*la1*sb**2-Zh21**2*la3
     &      *sb**2-Zh22**2*la2*sb**2+Zh22**2*la3*sb**2+Zh23
     &      **2*la1s*sb**2-Zh23**2*la2s*sb**2-2D0*Zh21*Zh22
     &      *la6*sb**2+2D0*Zh21*Zh22*la7*sb**2+2D0*Zh21**2
     &      *cb*la6*sb+2D0*Zh22**2*cb*la7*sb+2D0*Zh23**2
     &      *cb*la4s*sb+2D0*Zh23**2*cb*la5s*sb+2D0*Zh23**2
     &      *cb*la6s*sb-2D0*Zh21*Zh22*cb*la4*sb-2D0*Zh21
     &      *Zh22*cb*la5*sb
      AAABR(757) = EE**2*tB*(Zl11**2/CW**2-2D0*Zl12**2/CW
     &      **2+Zl11**2/SW**2-Zl11**2/CW**2*sb**2+2D0*Zl12
     &      **2/CW**2*sb**2-Zl11**2/SW**2*sb**2+Ml**2/MW
     &      **2/SW**2*Zl12**2)
      AAABR(758) = EE**2*cb*sb*(1D0/CW**2+1D0/SW**2)
      AAABR(759) = EE**2*tB*(Zl11/CW**2*Zl21-2D0*Zl12/CW
     &      **2*Zl22+Zl11/SW**2*Zl21-Zl11/CW**2*Zl21*sb**2
     &      +2D0*Zl12/CW**2*Zl22*sb**2-Zl11/SW**2*Zl21*sb
     &      **2+Ml**2/MW**2/SW**2*Zl12*Zl22)
      AAABR(760) = EE**2/CW**2*cb*sb
      AAABR(761) = EE**2*tB*(Zl21**2/CW**2-2D0*Zl22**2/CW
     &      **2+Zl21**2/SW**2-Zl21**2/CW**2*sb**2+2D0*Zl22
     &      **2/CW**2*sb**2-Zl21**2/SW**2*sb**2+Ml**2/MW
     &      **2/SW**2*Zl22**2)
      AAABR(762) = EE**2*tB*(1D0/CW**2-1D0/SW**2-sb**2/CW
     &      **2+sb**2/SW**2+Ml**2/MW**2/SW**2)
      AAABR(763) = EE**2*cb*sb*(1D0/CW**2-1D0/SW**2)
      AAABR(764) = EE**2/SW**2*Sqrt2*Zl11*tB*(2D0-2D0*sb
     &      **2-Ml**2/MW**2)
      AAABR(765) = EE**2/SW**2*Sqrt2*cb*sb
      AAABR(766) = EE**2/SW**2*Sqrt2*Zl21*tB*(2D0-2D0*sb
     &      **2-Ml**2/MW**2)
      AAABR(767) = EE**2*Pa21/SW**2*Sqrt2*(1D0-2D0*sb**2)
      AAABR(768) = EE**2*Pa11/SW**2*Sqrt2*(1D0-2D0*sb**2)
      AAABR(769) = 1D0/cb**2*(EE**2*Pa11/SW**2*Sqrt2*Zl11
     &      -3D0*EE**2*Pa11/SW**2*Sqrt2*Zl11*sb**2+2D0*EE
     &      **2*Pa11/SW**2*Sqrt2*Zl11*sb**4+2D0*EE*Lambda
     &      /MW*Ml*Pa12/SW*Zl12-2D0*EE*Lambda/MW*Ml*Pa12
     &      /SW*Zl12*sb**2+EE**2/MW**2*Ml**2*Pa11/SW**2*Sqrt2
     &      *Zl11*sb**2)
      AAABR(770) = 1D0/cb**2*(EE**2*Pa21/SW**2*Sqrt2*Zl11
     &      -3D0*EE**2*Pa21/SW**2*Sqrt2*Zl11*sb**2+2D0*EE
     &      **2*Pa21/SW**2*Sqrt2*Zl11*sb**4+2D0*EE*Lambda
     &      /MW*Ml*Pa22/SW*Zl12-2D0*EE*Lambda/MW*Ml*Pa22
     &      /SW*Zl12*sb**2+EE**2/MW**2*Ml**2*Pa21/SW**2*Sqrt2
     &      *Zl11*sb**2)
      AAABR(771) = 1D0/cb**2*(EE**2*Pa11/SW**2*Sqrt2*Zl21
     &      -3D0*EE**2*Pa11/SW**2*Sqrt2*Zl21*sb**2+2D0*EE
     &      **2*Pa11/SW**2*Sqrt2*Zl21*sb**4+2D0*EE*Lambda
     &      /MW*Ml*Pa12/SW*Zl22-2D0*EE*Lambda/MW*Ml*Pa12
     &      /SW*Zl22*sb**2+EE**2/MW**2*Ml**2*Pa11/SW**2*Sqrt2
     &      *Zl21*sb**2)
      AAABR(772) = 1D0/cb**2*(EE**2*Pa21/SW**2*Sqrt2*Zl21
     &      -3D0*EE**2*Pa21/SW**2*Sqrt2*Zl21*sb**2+2D0*EE
     &      **2*Pa21/SW**2*Sqrt2*Zl21*sb**4+2D0*EE*Lambda
     &      /MW*Ml*Pa22/SW*Zl22-2D0*EE*Lambda/MW*Ml*Pa22
     &      /SW*Zl22*sb**2+EE**2/MW**2*Ml**2*Pa21/SW**2*Sqrt2
     &      *Zl21*sb**2)
      AAABR(773) = EE**2/SW**2*Sqrt2*Zh12*Zl11/cb+EE**2/SW
     &      **2*Sqrt2*Zh11*Zl11/cb**2*sb-EE**2/SW**2*Sqrt2
     &      *Zh11*Zl11/cb**2*sb**3-EE**2/SW**2*Sqrt2*Zh12
     &      *Zl11/cb*sb**2-2D0*EE*Lambda/MW*Ml/SW*Zh13*Zl12
     &      /cb**2+2D0*EE*Lambda/MW*Ml/SW*Zh13*Zl12/cb**2
     &      *sb**2-EE**2/MW**2*Ml**2/SW**2*Sqrt2*Zh11*Zl11
     &      /cb**2*sb
      AAABR(774) = EE**2/SW**2*Sqrt2*(Zh31*sb+Zh32*cb)
      AAABR(775) = EE**2/SW**2*Sqrt2*(Zh11*sb+Zh12*cb)
      AAABR(776) = EE**2/SW**2*Sqrt2*Zh22*Zl11/cb+EE**2/SW
     &      **2*Sqrt2*Zh21*Zl11/cb**2*sb-EE**2/SW**2*Sqrt2
     &      *Zh21*Zl11/cb**2*sb**3-EE**2/SW**2*Sqrt2*Zh22
     &      *Zl11/cb*sb**2-2D0*EE*Lambda/MW*Ml/SW*Zh23*Zl12
     &      /cb**2+2D0*EE*Lambda/MW*Ml/SW*Zh23*Zl12/cb**2
     &      *sb**2-EE**2/MW**2*Ml**2/SW**2*Sqrt2*Zh21*Zl11
     &      /cb**2*sb
      AAABR(777) = EE**2/SW**2*Sqrt2*(Zh21*sb+Zh22*cb)
      AAABR(778) = EE**2/SW**2*Sqrt2*Zh32*Zl11/cb+EE**2/SW
     &      **2*Sqrt2*Zh31*Zl11/cb**2*sb-EE**2/SW**2*Sqrt2
     &      *Zh31*Zl11/cb**2*sb**3-EE**2/SW**2*Sqrt2*Zh32
     &      *Zl11/cb*sb**2-2D0*EE*Lambda/MW*Ml/SW*Zh33*Zl12
     &      /cb**2+2D0*EE*Lambda/MW*Ml/SW*Zh33*Zl12/cb**2
     &      *sb**2-EE**2/MW**2*Ml**2/SW**2*Sqrt2*Zh31*Zl11
     &      /cb**2*sb
      AAABR(779) = EE**2/SW**2*Sqrt2*Zh22*Zl21/cb+EE**2/SW
     &      **2*Sqrt2*Zh21*Zl21/cb**2*sb-EE**2/SW**2*Sqrt2
     &      *Zh21*Zl21/cb**2*sb**3-EE**2/SW**2*Sqrt2*Zh22
     &      *Zl21/cb*sb**2-2D0*EE*Lambda/MW*Ml/SW*Zh23*Zl22
     &      /cb**2+2D0*EE*Lambda/MW*Ml/SW*Zh23*Zl22/cb**2
     &      *sb**2-EE**2/MW**2*Ml**2/SW**2*Sqrt2*Zh21*Zl21
     &      /cb**2*sb
      AAABR(780) = EE**2/SW**2*Sqrt2*Zh12*Zl21/cb+EE**2/SW
     &      **2*Sqrt2*Zh11*Zl21/cb**2*sb-EE**2/SW**2*Sqrt2
     &      *Zh11*Zl21/cb**2*sb**3-EE**2/SW**2*Sqrt2*Zh12
     &      *Zl21/cb*sb**2-2D0*EE*Lambda/MW*Ml/SW*Zh13*Zl22
     &      /cb**2+2D0*EE*Lambda/MW*Ml/SW*Zh13*Zl22/cb**2
     &      *sb**2-EE**2/MW**2*Ml**2/SW**2*Sqrt2*Zh11*Zl21
     &      /cb**2*sb
      AAABR(781) = EE**2/SW**2*Sqrt2*Zh32*Zl21/cb+EE**2/SW
     &      **2*Sqrt2*Zh31*Zl21/cb**2*sb-EE**2/SW**2*Sqrt2
     &      *Zh31*Zl21/cb**2*sb**3-EE**2/SW**2*Sqrt2*Zh32
     &      *Zl21/cb*sb**2-2D0*EE*Lambda/MW*Ml/SW*Zh33*Zl22
     &      /cb**2+2D0*EE*Lambda/MW*Ml/SW*Zh33*Zl22/cb**2
     &      *sb**2-EE**2/MW**2*Ml**2/SW**2*Sqrt2*Zh31*Zl21
     &      /cb**2*sb
      AAABR(782) = EE**2/cb*(Zb11**2/CW**2*sb-Zb11**2/CW
     &      **2*sb**3+2D0*Zb12**2/CW**2*sb-2D0*Zb12**2/CW
     &      **2*sb**3-3D0*Zb11**2/SW**2*sb+3D0*Zb11**2/SW
     &      **2*sb**3-3D0*MbMM**2/MW**2/SW**2*Zb12**2*sb
     &      +3D0*MtMM**2/MW**2/SW**2*Zb11**2/sb-3D0*MtMM
     &      **2/MW**2/SW**2*Zb11**2*sb)
      AAABR(783) = EE**2/cb*(Zdd11**2/CW**2*sb-Zdd11**2/CW
     &      **2*sb**3+2D0*Zdd12**2/CW**2*sb-2D0*Zdd12**2
     &      /CW**2*sb**3-3D0*Zdd11**2/SW**2*sb+3D0*Zdd11
     &      **2/SW**2*sb**3-3D0*Mq**2/MW**2/SW**2*sb+3D0
     &      *Mq**2/MW**2/SW**2*Zdd11**2/sb)
      AAABR(784) = EE**2/cb*(Zss11**2/CW**2*sb-Zss11**2/CW
     &      **2*sb**3+2D0*Zss12**2/CW**2*sb-2D0*Zss12**2
     &      /CW**2*sb**3-3D0*Zss11**2/SW**2*sb+3D0*Zss11
     &      **2/SW**2*sb**3+3D0*Mc**2/MW**2/SW**2*Zss11**2
     &      /sb-3D0*Mc**2/MW**2/SW**2*Zss11**2*sb-3D0*Mq
     &      **2/MW**2/SW**2*Zss12**2*sb)
      AAABR(785) = EE**2/cb*(Zdd11/CW**2*Zdd21*sb-Zdd11/CW
     &      **2*Zdd21*sb**3+2D0*Zdd12/CW**2*Zdd22*sb-2D0
     &      *Zdd12/CW**2*Zdd22*sb**3-3D0*Zdd11/SW**2*Zdd21
     &      *sb+3D0*Zdd11/SW**2*Zdd21*sb**3+3D0*Mq**2/MW
     &      **2/SW**2*Zdd11*Zdd21/sb)
      AAABR(786) = EE**2/cb*(Zb11/CW**2*Zb21*sb-Zb11/CW**2
     &      *Zb21*sb**3+2D0*Zb12/CW**2*Zb22*sb-2D0*Zb12/CW
     &      **2*Zb22*sb**3-3D0*Zb11/SW**2*Zb21*sb+3D0*Zb11
     &      /SW**2*Zb21*sb**3-3D0*MbMM**2/MW**2/SW**2*Zb12
     &      *Zb22*sb+3D0*MtMM**2/MW**2/SW**2*Zb11*Zb21/sb
     &      -3D0*MtMM**2/MW**2/SW**2*Zb11*Zb21*sb)
      AAABR(787) = EE**2/cb*(Zss11/CW**2*Zss21*sb-Zss11/CW
     &      **2*Zss21*sb**3+2D0*Zss12/CW**2*Zss22*sb-2D0
     &      *Zss12/CW**2*Zss22*sb**3-3D0*Zss11/SW**2*Zss21
     &      *sb+3D0*Zss11/SW**2*Zss21*sb**3+3D0*Mc**2/MW
     &      **2/SW**2*Zss11*Zss21/sb-3D0*Mc**2/MW**2/SW**2
     &      *Zss11*Zss21*sb-3D0*Mq**2/MW**2/SW**2*Zss12*Zss22
     &      *sb)
      AAABR(788) = EE**2/cb*(Zdd21**2/CW**2*sb-Zdd21**2/CW
     &      **2*sb**3+2D0*Zdd22**2/CW**2*sb-2D0*Zdd22**2
     &      /CW**2*sb**3-3D0*Zdd21**2/SW**2*sb+3D0*Zdd21
     &      **2/SW**2*sb**3-3D0*Mq**2/MW**2/SW**2*sb+3D0
     &      *Mq**2/MW**2/SW**2*Zdd21**2/sb)
      AAABR(789) = EE**2/cb*(Zss21**2/CW**2*sb-Zss21**2/CW
     &      **2*sb**3+2D0*Zss22**2/CW**2*sb-2D0*Zss22**2
     &      /CW**2*sb**3-3D0*Zss21**2/SW**2*sb+3D0*Zss21
     &      **2/SW**2*sb**3+3D0*Mc**2/MW**2/SW**2*Zss21**2
     &      /sb-3D0*Mc**2/MW**2/SW**2*Zss21**2*sb-3D0*Mq
     &      **2/MW**2/SW**2*Zss22**2*sb)
      AAABR(790) = EE**2/cb*(Zb21**2/CW**2*sb-Zb21**2/CW
     &      **2*sb**3+2D0*Zb22**2/CW**2*sb-2D0*Zb22**2/CW
     &      **2*sb**3-3D0*Zb21**2/SW**2*sb+3D0*Zb21**2/SW
     &      **2*sb**3-3D0*MbMM**2/MW**2/SW**2*Zb22**2*sb
     &      +3D0*MtMM**2/MW**2/SW**2*Zb21**2/sb-3D0*MtMM
     &      **2/MW**2/SW**2*Zb21**2*sb)
      AAABR(791) = EE**2/cb*(Zcc11**2/CW**2*sb-Zcc11**2/CW
     &      **2*sb**3-4D0*Zcc12**2/CW**2*sb+4D0*Zcc12**2
     &      /CW**2*sb**3+3D0*Zcc11**2/SW**2*sb-3D0*Zcc11
     &      **2/SW**2*sb**3+3D0*Mc**2/MW**2/SW**2*Zcc12**2
     &      /sb-3D0*Mc**2/MW**2/SW**2*Zcc12**2*sb-3D0*Mq
     &      **2/MW**2/SW**2*Zcc11**2*sb)
      AAABR(792) = EE**2/cb*(Zuu11**2/CW**2*sb-Zuu11**2/CW
     &      **2*sb**3-4D0*Zuu12**2/CW**2*sb+4D0*Zuu12**2
     &      /CW**2*sb**3+3D0*Zuu11**2/SW**2*sb-3D0*Zuu11
     &      **2/SW**2*sb**3-3D0*Mq**2/MW**2/SW**2*sb+3D0
     &      *Mq**2/MW**2/SW**2*Zuu12**2/sb)
      AAABR(793) = EE**2/cb*(Zt11**2/CW**2*sb-Zt11**2/CW
     &      **2*sb**3-4D0*Zt12**2/CW**2*sb+4D0*Zt12**2/CW
     &      **2*sb**3+3D0*Zt11**2/SW**2*sb-3D0*Zt11**2/SW
     &      **2*sb**3-3D0*MbMM**2/MW**2/SW**2*Zt11**2*sb
     &      +3D0*MtMM**2/MW**2/SW**2*Zt12**2/sb-3D0*MtMM
     &      **2/MW**2/SW**2*Zt12**2*sb)
      AAABR(794) = EE**2/cb*(Zcc11/CW**2*Zcc21*sb-Zcc11/CW
     &      **2*Zcc21*sb**3-4D0*Zcc12/CW**2*Zcc22*sb+4D0
     &      *Zcc12/CW**2*Zcc22*sb**3+3D0*Zcc11/SW**2*Zcc21
     &      *sb-3D0*Zcc11/SW**2*Zcc21*sb**3+3D0*Mc**2/MW
     &      **2/SW**2*Zcc12*Zcc22/sb-3D0*Mc**2/MW**2/SW**2
     &      *Zcc12*Zcc22*sb-3D0*Mq**2/MW**2/SW**2*Zcc11*Zcc21
     &      *sb)
      AAABR(795) = EE**2/cb*(Zuu11/CW**2*Zuu21*sb-Zuu11/CW
     &      **2*Zuu21*sb**3-4D0*Zuu12/CW**2*Zuu22*sb+4D0
     &      *Zuu12/CW**2*Zuu22*sb**3+3D0*Zuu11/SW**2*Zuu21
     &      *sb-3D0*Zuu11/SW**2*Zuu21*sb**3+3D0*Mq**2/MW
     &      **2/SW**2*Zuu12*Zuu22/sb)
      AAABR(796) = EE**2/cb*(Zt11/CW**2*Zt21*sb-Zt11/CW**2
     &      *Zt21*sb**3-4D0*Zt12/CW**2*Zt22*sb+4D0*Zt12/CW
     &      **2*Zt22*sb**3+3D0*Zt11/SW**2*Zt21*sb-3D0*Zt11
     &      /SW**2*Zt21*sb**3-3D0*MbMM**2/MW**2/SW**2*Zt11
     &      *Zt21*sb+3D0*MtMM**2/MW**2/SW**2*Zt12*Zt22/sb
     &      -3D0*MtMM**2/MW**2/SW**2*Zt12*Zt22*sb)
      AAABR(797) = EE**2/SW**2*Sqrt2/cb*(2D0*Zcc11*Zss11
     &      *sb-2D0*Zcc11*Zss11*sb**3-Mc**2/MW**2*Zcc11*Zss11
     &      /sb+Mc**2/MW**2*Zcc11*Zss11*sb-Mq**2/MW**2*Zcc11
     &      *Zss11*sb+Mc/MW**2*Mq*Zcc12*Zss12/sb)
      AAABR(798) = EE**2/SW**2*Sqrt2/cb*(2D0*Zdd11*Zuu11
     &      *sb-2D0*Zdd11*Zuu11*sb**3-Mq**2/MW**2*Zdd11*Zuu11
     &      /sb+Mq**2/MW**2*Zdd12*Zuu12/sb)
      AAABR(799) = EE**2/SW**2*Sqrt2/cb*(2D0*Zb11*Zt11*sb
     &      -2D0*Zb11*Zt11*sb**3-MbMM**2/MW**2*Zb11*Zt11
     &      *sb-MtMM**2/MW**2*Zb11*Zt11/sb+MtMM**2/MW**2
     &      *Zb11*Zt11*sb+MbMM/MW**2*MtMM*Zb12*Zt12/sb)
      AAABR(800) = EE**2/SW**2*Sqrt2/cb*(2D0*Zcc11*Zss21
     &      *sb-2D0*Zcc11*Zss21*sb**3-Mc**2/MW**2*Zcc11*Zss21
     &      /sb+Mc**2/MW**2*Zcc11*Zss21*sb-Mq**2/MW**2*Zcc11
     &      *Zss21*sb+Mc/MW**2*Mq*Zcc12*Zss22/sb)
      AAABR(801) = EE**2/SW**2*Sqrt2/cb*(2D0*Zb21*Zt11*sb
     &      -2D0*Zb21*Zt11*sb**3-MbMM**2/MW**2*Zb21*Zt11
     &      *sb-MtMM**2/MW**2*Zb21*Zt11/sb+MtMM**2/MW**2
     &      *Zb21*Zt11*sb+MbMM/MW**2*MtMM*Zb22*Zt12/sb)
      AAABR(802) = EE**2/SW**2*Sqrt2/cb*(2D0*Zdd21*Zuu11
     &      *sb-2D0*Zdd21*Zuu11*sb**3-Mq**2/MW**2*Zdd21*Zuu11
     &      /sb+Mq**2/MW**2*Zdd22*Zuu12/sb)
      AAABR(803) = 1D0/cb**2*(EE**2*Pa11/SW**2*Sqrt2*Zdd11
     &      *Zuu11-3D0*EE**2*Pa11/SW**2*Sqrt2*Zdd11*Zuu11
     &      *sb**2+2D0*EE**2*Pa11/SW**2*Sqrt2*Zdd11*Zuu11
     &      *sb**4-2D0*EE*Lambda/MW*Mq*Pa12/SW*Zdd11*Zuu12
     &      +2D0*EE*Lambda/MW*Mq*Pa12/SW*Zdd12*Zuu11+2D0
     &      *EE**2/MW**2*Mq**2*Pa11/SW**2*Sqrt2*Zdd11*Zuu11
     &      +2D0*EE*Lambda/MW*Mq*Pa12/SW*Zdd11*Zuu12*sb**2
     &      -2D0*EE*Lambda/MW*Mq*Pa12/SW*Zdd12*Zuu11*sb**2
     &      -EE**2/MW**2*Mq**2*Pa11/SW**2*Sqrt2*Zdd11*Zuu11
     &      /sb**2)
      AAABR(804) = 1D0/cb**2*(EE**2*Pa11/SW**2*Sqrt2*Zb11
     &      *Zt11-3D0*EE**2*Pa11/SW**2*Sqrt2*Zb11*Zt11*sb
     &      **2+2D0*EE**2*Pa11/SW**2*Sqrt2*Zb11*Zt11*sb**4
     &      +2D0*EE*Lambda/MW*MbMM*Pa12/SW*Zb12*Zt11-2D0
     &      *EE*Lambda/MW*MtMM*Pa12/SW*Zb11*Zt12+2D0*EE**2
     &      /MW**2*MtMM**2*Pa11/SW**2*Sqrt2*Zb11*Zt11-2D0
     &      *EE*Lambda/MW*MbMM*Pa12/SW*Zb12*Zt11*sb**2+2D0
     &      *EE*Lambda/MW*MtMM*Pa12/SW*Zb11*Zt12*sb**2+EE
     &      **2/MW**2*MbMM**2*Pa11/SW**2*Sqrt2*Zb11*Zt11
     &      *sb**2-EE**2/MW**2*MtMM**2*Pa11/SW**2*Sqrt2*Zb11
     &      *Zt11/sb**2-EE**2/MW**2*MtMM**2*Pa11/SW**2*Sqrt2
     &      *Zb11*Zt11*sb**2)
      AAABR(805) = 1D0/cb**2*(EE**2*Pa21/SW**2*Sqrt2*Zdd11
     &      *Zuu11-3D0*EE**2*Pa21/SW**2*Sqrt2*Zdd11*Zuu11
     &      *sb**2+2D0*EE**2*Pa21/SW**2*Sqrt2*Zdd11*Zuu11
     &      *sb**4-2D0*EE*Lambda/MW*Mq*Pa22/SW*Zdd11*Zuu12
     &      +2D0*EE*Lambda/MW*Mq*Pa22/SW*Zdd12*Zuu11+2D0
     &      *EE**2/MW**2*Mq**2*Pa21/SW**2*Sqrt2*Zdd11*Zuu11
     &      +2D0*EE*Lambda/MW*Mq*Pa22/SW*Zdd11*Zuu12*sb**2
     &      -2D0*EE*Lambda/MW*Mq*Pa22/SW*Zdd12*Zuu11*sb**2
     &      -EE**2/MW**2*Mq**2*Pa21/SW**2*Sqrt2*Zdd11*Zuu11
     &      /sb**2)
      AAABR(806) = 1D0/cb**2*(EE**2*Pa21/SW**2*Sqrt2*Zb11
     &      *Zt11-3D0*EE**2*Pa21/SW**2*Sqrt2*Zb11*Zt11*sb
     &      **2+2D0*EE**2*Pa21/SW**2*Sqrt2*Zb11*Zt11*sb**4
     &      +2D0*EE*Lambda/MW*MbMM*Pa22/SW*Zb12*Zt11-2D0
     &      *EE*Lambda/MW*MtMM*Pa22/SW*Zb11*Zt12+2D0*EE**2
     &      /MW**2*MtMM**2*Pa21/SW**2*Sqrt2*Zb11*Zt11-2D0
     &      *EE*Lambda/MW*MbMM*Pa22/SW*Zb12*Zt11*sb**2+2D0
     &      *EE*Lambda/MW*MtMM*Pa22/SW*Zb11*Zt12*sb**2+EE
     &      **2/MW**2*MbMM**2*Pa21/SW**2*Sqrt2*Zb11*Zt11
     &      *sb**2-EE**2/MW**2*MtMM**2*Pa21/SW**2*Sqrt2*Zb11
     &      *Zt11/sb**2-EE**2/MW**2*MtMM**2*Pa21/SW**2*Sqrt2
     &      *Zb11*Zt11*sb**2)
      AAABR(807) = 1D0/cb**2*(EE**2*Pa11/SW**2*Sqrt2*Zcc11
     &      *Zss11-3D0*EE**2*Pa11/SW**2*Sqrt2*Zcc11*Zss11
     &      *sb**2+2D0*EE**2*Pa11/SW**2*Sqrt2*Zcc11*Zss11
     &      *sb**4-2D0*EE*Lambda/MW*Mc*Pa12/SW*Zcc12*Zss11
     &      +2D0*EE*Lambda/MW*Mq*Pa12/SW*Zcc11*Zss12+2D0
     &      *EE**2/MW**2*Mc**2*Pa11/SW**2*Sqrt2*Zcc11*Zss11
     &      +2D0*EE*Lambda/MW*Mc*Pa12/SW*Zcc12*Zss11*sb**2
     &      -2D0*EE*Lambda/MW*Mq*Pa12/SW*Zcc11*Zss12*sb**2
     &      -EE**2/MW**2*Mc**2*Pa11/SW**2*Sqrt2*Zcc11*Zss11
     &      /sb**2-EE**2/MW**2*Mc**2*Pa11/SW**2*Sqrt2*Zcc11
     &      *Zss11*sb**2+EE**2/MW**2*Mq**2*Pa11/SW**2*Sqrt2
     &      *Zcc11*Zss11*sb**2)
      AAABR(808) = 1D0/cb**2*(EE**2*Pa21/SW**2*Sqrt2*Zcc11
     &      *Zss11-3D0*EE**2*Pa21/SW**2*Sqrt2*Zcc11*Zss11
     &      *sb**2+2D0*EE**2*Pa21/SW**2*Sqrt2*Zcc11*Zss11
     &      *sb**4-2D0*EE*Lambda/MW*Mc*Pa22/SW*Zcc12*Zss11
     &      +2D0*EE*Lambda/MW*Mq*Pa22/SW*Zcc11*Zss12+2D0
     &      *EE**2/MW**2*Mc**2*Pa21/SW**2*Sqrt2*Zcc11*Zss11
     &      +2D0*EE*Lambda/MW*Mc*Pa22/SW*Zcc12*Zss11*sb**2
     &      -2D0*EE*Lambda/MW*Mq*Pa22/SW*Zcc11*Zss12*sb**2
     &      -EE**2/MW**2*Mc**2*Pa21/SW**2*Sqrt2*Zcc11*Zss11
     &      /sb**2-EE**2/MW**2*Mc**2*Pa21/SW**2*Sqrt2*Zcc11
     &      *Zss11*sb**2+EE**2/MW**2*Mq**2*Pa21/SW**2*Sqrt2
     &      *Zcc11*Zss11*sb**2)
      AAABR(809) = 1D0/cb**2*(EE**2*Pa11/SW**2*Sqrt2*Zcc11
     &      *Zss21-3D0*EE**2*Pa11/SW**2*Sqrt2*Zcc11*Zss21
     &      *sb**2+2D0*EE**2*Pa11/SW**2*Sqrt2*Zcc11*Zss21
     &      *sb**4-2D0*EE*Lambda/MW*Mc*Pa12/SW*Zcc12*Zss21
     &      +2D0*EE*Lambda/MW*Mq*Pa12/SW*Zcc11*Zss22+2D0
     &      *EE**2/MW**2*Mc**2*Pa11/SW**2*Sqrt2*Zcc11*Zss21
     &      +2D0*EE*Lambda/MW*Mc*Pa12/SW*Zcc12*Zss21*sb**2
     &      -2D0*EE*Lambda/MW*Mq*Pa12/SW*Zcc11*Zss22*sb**2
     &      -EE**2/MW**2*Mc**2*Pa11/SW**2*Sqrt2*Zcc11*Zss21
     &      /sb**2-EE**2/MW**2*Mc**2*Pa11/SW**2*Sqrt2*Zcc11
     &      *Zss21*sb**2+EE**2/MW**2*Mq**2*Pa11/SW**2*Sqrt2
     &      *Zcc11*Zss21*sb**2)
      AAABR(810) = 1D0/cb**2*(EE**2*Pa21/SW**2*Sqrt2*Zcc11
     &      *Zss21-3D0*EE**2*Pa21/SW**2*Sqrt2*Zcc11*Zss21
     &      *sb**2+2D0*EE**2*Pa21/SW**2*Sqrt2*Zcc11*Zss21
     &      *sb**4-2D0*EE*Lambda/MW*Mc*Pa22/SW*Zcc12*Zss21
     &      +2D0*EE*Lambda/MW*Mq*Pa22/SW*Zcc11*Zss22+2D0
     &      *EE**2/MW**2*Mc**2*Pa21/SW**2*Sqrt2*Zcc11*Zss21
     &      +2D0*EE*Lambda/MW*Mc*Pa22/SW*Zcc12*Zss21*sb**2
     &      -2D0*EE*Lambda/MW*Mq*Pa22/SW*Zcc11*Zss22*sb**2
     &      -EE**2/MW**2*Mc**2*Pa21/SW**2*Sqrt2*Zcc11*Zss21
     &      /sb**2-EE**2/MW**2*Mc**2*Pa21/SW**2*Sqrt2*Zcc11
     &      *Zss21*sb**2+EE**2/MW**2*Mq**2*Pa21/SW**2*Sqrt2
     &      *Zcc11*Zss21*sb**2)
      AAABR(811) = 1D0/cb**2*(EE**2*Pa11/SW**2*Sqrt2*Zdd21
     &      *Zuu11-3D0*EE**2*Pa11/SW**2*Sqrt2*Zdd21*Zuu11
     &      *sb**2+2D0*EE**2*Pa11/SW**2*Sqrt2*Zdd21*Zuu11
     &      *sb**4-2D0*EE*Lambda/MW*Mq*Pa12/SW*Zdd21*Zuu12
     &      +2D0*EE*Lambda/MW*Mq*Pa12/SW*Zdd22*Zuu11+2D0
     &      *EE**2/MW**2*Mq**2*Pa11/SW**2*Sqrt2*Zdd21*Zuu11
     &      +2D0*EE*Lambda/MW*Mq*Pa12/SW*Zdd21*Zuu12*sb**2
     &      -2D0*EE*Lambda/MW*Mq*Pa12/SW*Zdd22*Zuu11*sb**2
     &      -EE**2/MW**2*Mq**2*Pa11/SW**2*Sqrt2*Zdd21*Zuu11
     &      /sb**2)
      AAABR(812) = 1D0/cb**2*(EE**2*Pa21/SW**2*Sqrt2*Zb21
     &      *Zt11-3D0*EE**2*Pa21/SW**2*Sqrt2*Zb21*Zt11*sb
     &      **2+2D0*EE**2*Pa21/SW**2*Sqrt2*Zb21*Zt11*sb**4
     &      +2D0*EE*Lambda/MW*MbMM*Pa22/SW*Zb22*Zt11-2D0
     &      *EE*Lambda/MW*MtMM*Pa22/SW*Zb21*Zt12+2D0*EE**2
     &      /MW**2*MtMM**2*Pa21/SW**2*Sqrt2*Zb21*Zt11-2D0
     &      *EE*Lambda/MW*MbMM*Pa22/SW*Zb22*Zt11*sb**2+2D0
     &      *EE*Lambda/MW*MtMM*Pa22/SW*Zb21*Zt12*sb**2+EE
     &      **2/MW**2*MbMM**2*Pa21/SW**2*Sqrt2*Zb21*Zt11
     &      *sb**2-EE**2/MW**2*MtMM**2*Pa21/SW**2*Sqrt2*Zb21
     &      *Zt11/sb**2-EE**2/MW**2*MtMM**2*Pa21/SW**2*Sqrt2
     &      *Zb21*Zt11*sb**2)
      AAABR(813) = 1D0/cb**2*(EE**2*Pa11/SW**2*Sqrt2*Zb21
     &      *Zt11-3D0*EE**2*Pa11/SW**2*Sqrt2*Zb21*Zt11*sb
     &      **2+2D0*EE**2*Pa11/SW**2*Sqrt2*Zb21*Zt11*sb**4
     &      +2D0*EE*Lambda/MW*MbMM*Pa12/SW*Zb22*Zt11-2D0
     &      *EE*Lambda/MW*MtMM*Pa12/SW*Zb21*Zt12+2D0*EE**2
     &      /MW**2*MtMM**2*Pa11/SW**2*Sqrt2*Zb21*Zt11-2D0
     &      *EE*Lambda/MW*MbMM*Pa12/SW*Zb22*Zt11*sb**2+2D0
     &      *EE*Lambda/MW*MtMM*Pa12/SW*Zb21*Zt12*sb**2+EE
     &      **2/MW**2*MbMM**2*Pa11/SW**2*Sqrt2*Zb21*Zt11
     &      *sb**2-EE**2/MW**2*MtMM**2*Pa11/SW**2*Sqrt2*Zb21
     &      *Zt11/sb**2-EE**2/MW**2*MtMM**2*Pa11/SW**2*Sqrt2
     &      *Zb21*Zt11*sb**2)
      AAABR(814) = 1D0/cb**2*(EE**2*Pa21/SW**2*Sqrt2*Zdd21
     &      *Zuu11-3D0*EE**2*Pa21/SW**2*Sqrt2*Zdd21*Zuu11
     &      *sb**2+2D0*EE**2*Pa21/SW**2*Sqrt2*Zdd21*Zuu11
     &      *sb**4-2D0*EE*Lambda/MW*Mq*Pa22/SW*Zdd21*Zuu12
     &      +2D0*EE*Lambda/MW*Mq*Pa22/SW*Zdd22*Zuu11+2D0
     &      *EE**2/MW**2*Mq**2*Pa21/SW**2*Sqrt2*Zdd21*Zuu11
     &      +2D0*EE*Lambda/MW*Mq*Pa22/SW*Zdd21*Zuu12*sb**2
     &      -2D0*EE*Lambda/MW*Mq*Pa22/SW*Zdd22*Zuu11*sb**2
     &      -EE**2/MW**2*Mq**2*Pa21/SW**2*Sqrt2*Zdd21*Zuu11
     &      /sb**2)
      AAABR(815) = EE**2/SW**2*Sqrt2*Zb11*Zh12*Zt11/cb+EE
     &      **2/SW**2*Sqrt2*Zb11*Zh11*Zt11/cb**2*sb-EE**2
     &      /SW**2*Sqrt2*Zb11*Zh11*Zt11/cb**2*sb**3-EE**2
     &      /SW**2*Sqrt2*Zb11*Zh12*Zt11/cb*sb**2-2D0*EE*Lambda
     &      /MW*MbMM/SW*Zb12*Zh13*Zt11/cb**2-2D0*EE*Lambda
     &      /MW*MtMM/SW*Zb11*Zh13*Zt12/cb**2+EE**2/MW**2
     &      *MtMM**2/SW**2*Sqrt2*Zb11*Zh12*Zt11/cb+2D0*EE
     &      *Lambda/MW*MbMM/SW*Zb12*Zh13*Zt11/cb**2*sb**2
     &      +2D0*EE*Lambda/MW*MtMM/SW*Zb11*Zh13*Zt12/cb**2
     &      *sb**2-EE**2/MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb12
     &      *Zh12*Zt12/cb-EE**2/MW**2*MbMM**2/SW**2*Sqrt2
     &      *Zb11*Zh11*Zt11/cb**2*sb-EE**2/MW**2*MtMM**2
     &      /SW**2*Sqrt2*Zb11*Zh12*Zt11/cb/sb**2-EE**2/MW
     &      **2*MbMM*MtMM/SW**2*Sqrt2*Zb12*Zh11*Zt12/cb**2
     &      /sb+EE**2/MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb12*Zh11
     &      *Zt12/cb**2*sb
      AAABR(816) = EE**2/SW**2*Sqrt2*Zdd11*Zh12*Zuu11/cb
     &      +EE**2/SW**2*Sqrt2*Zdd11*Zh11*Zuu11/cb**2*sb
     &      -EE**2/SW**2*Sqrt2*Zdd11*Zh11*Zuu11/cb**2*sb
     &      **3-EE**2/SW**2*Sqrt2*Zdd11*Zh12*Zuu11/cb*sb
     &      **2-2D0*EE*Lambda/MW*Mq/SW*Zdd11*Zh13*Zuu12/cb
     &      **2-2D0*EE*Lambda/MW*Mq/SW*Zdd12*Zh13*Zuu11/cb
     &      **2+EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd11*Zh12
     &      *Zuu11/cb-EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd12
     &      *Zh12*Zuu12/cb+2D0*EE*Lambda/MW*Mq/SW*Zdd11*Zh13
     &      *Zuu12/cb**2*sb**2+2D0*EE*Lambda/MW*Mq/SW*Zdd12
     &      *Zh13*Zuu11/cb**2*sb**2-EE**2/MW**2*Mq**2/SW
     &      **2*Sqrt2*Zdd11*Zh11*Zuu11/cb**2*sb-EE**2/MW
     &      **2*Mq**2/SW**2*Sqrt2*Zdd11*Zh12*Zuu11/cb/sb
     &      **2-EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd12*Zh11
     &      *Zuu12/cb**2/sb+EE**2/MW**2*Mq**2/SW**2*Sqrt2
     &      *Zdd12*Zh11*Zuu12/cb**2*sb
      AAABR(817) = EE**2/SW**2*Sqrt2*Zcc11*Zh32*Zss11/cb
     &      +EE**2/SW**2*Sqrt2*Zcc11*Zh31*Zss11/cb**2*sb
     &      -EE**2/SW**2*Sqrt2*Zcc11*Zh31*Zss11/cb**2*sb
     &      **3-EE**2/SW**2*Sqrt2*Zcc11*Zh32*Zss11/cb*sb
     &      **2-2D0*EE*Lambda/MW*Mc/SW*Zcc12*Zh33*Zss11/cb
     &      **2-2D0*EE*Lambda/MW*Mq/SW*Zcc11*Zh33*Zss12/cb
     &      **2+EE**2/MW**2*Mc**2/SW**2*Sqrt2*Zcc11*Zh32
     &      *Zss11/cb+2D0*EE*Lambda/MW*Mc/SW*Zcc12*Zh33*Zss11
     &      /cb**2*sb**2+2D0*EE*Lambda/MW*Mq/SW*Zcc11*Zh33
     &      *Zss12/cb**2*sb**2-EE**2/MW**2*Mc*Mq/SW**2*Sqrt2
     &      *Zcc12*Zh32*Zss12/cb-EE**2/MW**2*Mc**2/SW**2
     &      *Sqrt2*Zcc11*Zh32*Zss11/cb/sb**2-EE**2/MW**2
     &      *Mq**2/SW**2*Sqrt2*Zcc11*Zh31*Zss11/cb**2*sb
     &      -EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc12*Zh31*Zss12
     &      /cb**2/sb+EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc12
     &      *Zh31*Zss12/cb**2*sb
      AAABR(818) = EE**2/SW**2*Sqrt2*Zdd11*Zh32*Zuu11/cb
     &      +EE**2/SW**2*Sqrt2*Zdd11*Zh31*Zuu11/cb**2*sb
     &      -EE**2/SW**2*Sqrt2*Zdd11*Zh31*Zuu11/cb**2*sb
     &      **3-EE**2/SW**2*Sqrt2*Zdd11*Zh32*Zuu11/cb*sb
     &      **2-2D0*EE*Lambda/MW*Mq/SW*Zdd11*Zh33*Zuu12/cb
     &      **2-2D0*EE*Lambda/MW*Mq/SW*Zdd12*Zh33*Zuu11/cb
     &      **2+EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd11*Zh32
     &      *Zuu11/cb-EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd12
     &      *Zh32*Zuu12/cb+2D0*EE*Lambda/MW*Mq/SW*Zdd11*Zh33
     &      *Zuu12/cb**2*sb**2+2D0*EE*Lambda/MW*Mq/SW*Zdd12
     &      *Zh33*Zuu11/cb**2*sb**2-EE**2/MW**2*Mq**2/SW
     &      **2*Sqrt2*Zdd11*Zh31*Zuu11/cb**2*sb-EE**2/MW
     &      **2*Mq**2/SW**2*Sqrt2*Zdd11*Zh32*Zuu11/cb/sb
     &      **2-EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd12*Zh31
     &      *Zuu12/cb**2/sb+EE**2/MW**2*Mq**2/SW**2*Sqrt2
     &      *Zdd12*Zh31*Zuu12/cb**2*sb
      AAABR(819) = EE**2/SW**2*Sqrt2*Zcc11*Zh12*Zss11/cb
     &      +EE**2/SW**2*Sqrt2*Zcc11*Zh11*Zss11/cb**2*sb
     &      -EE**2/SW**2*Sqrt2*Zcc11*Zh11*Zss11/cb**2*sb
     &      **3-EE**2/SW**2*Sqrt2*Zcc11*Zh12*Zss11/cb*sb
     &      **2-2D0*EE*Lambda/MW*Mc/SW*Zcc12*Zh13*Zss11/cb
     &      **2-2D0*EE*Lambda/MW*Mq/SW*Zcc11*Zh13*Zss12/cb
     &      **2+EE**2/MW**2*Mc**2/SW**2*Sqrt2*Zcc11*Zh12
     &      *Zss11/cb+2D0*EE*Lambda/MW*Mc/SW*Zcc12*Zh13*Zss11
     &      /cb**2*sb**2+2D0*EE*Lambda/MW*Mq/SW*Zcc11*Zh13
     &      *Zss12/cb**2*sb**2-EE**2/MW**2*Mc*Mq/SW**2*Sqrt2
     &      *Zcc12*Zh12*Zss12/cb-EE**2/MW**2*Mc**2/SW**2
     &      *Sqrt2*Zcc11*Zh12*Zss11/cb/sb**2-EE**2/MW**2
     &      *Mq**2/SW**2*Sqrt2*Zcc11*Zh11*Zss11/cb**2*sb
     &      -EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc12*Zh11*Zss12
     &      /cb**2/sb+EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc12
     &      *Zh11*Zss12/cb**2*sb
      AAABR(820) = EE**2/SW**2*Sqrt2*Zdd11*Zh22*Zuu11/cb
     &      +EE**2/SW**2*Sqrt2*Zdd11*Zh21*Zuu11/cb**2*sb
     &      -EE**2/SW**2*Sqrt2*Zdd11*Zh21*Zuu11/cb**2*sb
     &      **3-EE**2/SW**2*Sqrt2*Zdd11*Zh22*Zuu11/cb*sb
     &      **2-2D0*EE*Lambda/MW*Mq/SW*Zdd11*Zh23*Zuu12/cb
     &      **2-2D0*EE*Lambda/MW*Mq/SW*Zdd12*Zh23*Zuu11/cb
     &      **2+EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd11*Zh22
     &      *Zuu11/cb-EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd12
     &      *Zh22*Zuu12/cb+2D0*EE*Lambda/MW*Mq/SW*Zdd11*Zh23
     &      *Zuu12/cb**2*sb**2+2D0*EE*Lambda/MW*Mq/SW*Zdd12
     &      *Zh23*Zuu11/cb**2*sb**2-EE**2/MW**2*Mq**2/SW
     &      **2*Sqrt2*Zdd11*Zh21*Zuu11/cb**2*sb-EE**2/MW
     &      **2*Mq**2/SW**2*Sqrt2*Zdd11*Zh22*Zuu11/cb/sb
     &      **2-EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd12*Zh21
     &      *Zuu12/cb**2/sb+EE**2/MW**2*Mq**2/SW**2*Sqrt2
     &      *Zdd12*Zh21*Zuu12/cb**2*sb
      AAABR(821) = EE**2/SW**2*Sqrt2*Zb11*Zh32*Zt11/cb+EE
     &      **2/SW**2*Sqrt2*Zb11*Zh31*Zt11/cb**2*sb-EE**2
     &      /SW**2*Sqrt2*Zb11*Zh31*Zt11/cb**2*sb**3-EE**2
     &      /SW**2*Sqrt2*Zb11*Zh32*Zt11/cb*sb**2-2D0*EE*Lambda
     &      /MW*MbMM/SW*Zb12*Zh33*Zt11/cb**2-2D0*EE*Lambda
     &      /MW*MtMM/SW*Zb11*Zh33*Zt12/cb**2+EE**2/MW**2
     &      *MtMM**2/SW**2*Sqrt2*Zb11*Zh32*Zt11/cb+2D0*EE
     &      *Lambda/MW*MbMM/SW*Zb12*Zh33*Zt11/cb**2*sb**2
     &      +2D0*EE*Lambda/MW*MtMM/SW*Zb11*Zh33*Zt12/cb**2
     &      *sb**2-EE**2/MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb12
     &      *Zh32*Zt12/cb-EE**2/MW**2*MbMM**2/SW**2*Sqrt2
     &      *Zb11*Zh31*Zt11/cb**2*sb-EE**2/MW**2*MtMM**2
     &      /SW**2*Sqrt2*Zb11*Zh32*Zt11/cb/sb**2-EE**2/MW
     &      **2*MbMM*MtMM/SW**2*Sqrt2*Zb12*Zh31*Zt12/cb**2
     &      /sb+EE**2/MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb12*Zh31
     &      *Zt12/cb**2*sb
      AAABR(822) = EE**2/SW**2*Sqrt2*Zb11*Zh22*Zt11/cb+EE
     &      **2/SW**2*Sqrt2*Zb11*Zh21*Zt11/cb**2*sb-EE**2
     &      /SW**2*Sqrt2*Zb11*Zh21*Zt11/cb**2*sb**3-EE**2
     &      /SW**2*Sqrt2*Zb11*Zh22*Zt11/cb*sb**2-2D0*EE*Lambda
     &      /MW*MbMM/SW*Zb12*Zh23*Zt11/cb**2-2D0*EE*Lambda
     &      /MW*MtMM/SW*Zb11*Zh23*Zt12/cb**2+EE**2/MW**2
     &      *MtMM**2/SW**2*Sqrt2*Zb11*Zh22*Zt11/cb+2D0*EE
     &      *Lambda/MW*MbMM/SW*Zb12*Zh23*Zt11/cb**2*sb**2
     &      +2D0*EE*Lambda/MW*MtMM/SW*Zb11*Zh23*Zt12/cb**2
     &      *sb**2-EE**2/MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb12
     &      *Zh22*Zt12/cb-EE**2/MW**2*MbMM**2/SW**2*Sqrt2
     &      *Zb11*Zh21*Zt11/cb**2*sb-EE**2/MW**2*MtMM**2
     &      /SW**2*Sqrt2*Zb11*Zh22*Zt11/cb/sb**2-EE**2/MW
     &      **2*MbMM*MtMM/SW**2*Sqrt2*Zb12*Zh21*Zt12/cb**2
     &      /sb+EE**2/MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb12*Zh21
     &      *Zt12/cb**2*sb
      AAABR(823) = EE**2/SW**2*Sqrt2*Zcc11*Zh22*Zss11/cb
     &      +EE**2/SW**2*Sqrt2*Zcc11*Zh21*Zss11/cb**2*sb
     &      -EE**2/SW**2*Sqrt2*Zcc11*Zh21*Zss11/cb**2*sb
     &      **3-EE**2/SW**2*Sqrt2*Zcc11*Zh22*Zss11/cb*sb
     &      **2-2D0*EE*Lambda/MW*Mc/SW*Zcc12*Zh23*Zss11/cb
     &      **2-2D0*EE*Lambda/MW*Mq/SW*Zcc11*Zh23*Zss12/cb
     &      **2+EE**2/MW**2*Mc**2/SW**2*Sqrt2*Zcc11*Zh22
     &      *Zss11/cb+2D0*EE*Lambda/MW*Mc/SW*Zcc12*Zh23*Zss11
     &      /cb**2*sb**2+2D0*EE*Lambda/MW*Mq/SW*Zcc11*Zh23
     &      *Zss12/cb**2*sb**2-EE**2/MW**2*Mc*Mq/SW**2*Sqrt2
     &      *Zcc12*Zh22*Zss12/cb-EE**2/MW**2*Mc**2/SW**2
     &      *Sqrt2*Zcc11*Zh22*Zss11/cb/sb**2-EE**2/MW**2
     &      *Mq**2/SW**2*Sqrt2*Zcc11*Zh21*Zss11/cb**2*sb
     &      -EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc12*Zh21*Zss12
     &      /cb**2/sb+EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc12
     &      *Zh21*Zss12/cb**2*sb
      AAABR(824) = EE**2/SW**2*Sqrt2*Zcc11*Zh22*Zss21/cb
     &      +EE**2/SW**2*Sqrt2*Zcc11*Zh21*Zss21/cb**2*sb
     &      -EE**2/SW**2*Sqrt2*Zcc11*Zh21*Zss21/cb**2*sb
     &      **3-EE**2/SW**2*Sqrt2*Zcc11*Zh22*Zss21/cb*sb
     &      **2-2D0*EE*Lambda/MW*Mc/SW*Zcc12*Zh23*Zss21/cb
     &      **2-2D0*EE*Lambda/MW*Mq/SW*Zcc11*Zh23*Zss22/cb
     &      **2+EE**2/MW**2*Mc**2/SW**2*Sqrt2*Zcc11*Zh22
     &      *Zss21/cb+2D0*EE*Lambda/MW*Mc/SW*Zcc12*Zh23*Zss21
     &      /cb**2*sb**2+2D0*EE*Lambda/MW*Mq/SW*Zcc11*Zh23
     &      *Zss22/cb**2*sb**2-EE**2/MW**2*Mc*Mq/SW**2*Sqrt2
     &      *Zcc12*Zh22*Zss22/cb-EE**2/MW**2*Mc**2/SW**2
     &      *Sqrt2*Zcc11*Zh22*Zss21/cb/sb**2-EE**2/MW**2
     &      *Mq**2/SW**2*Sqrt2*Zcc11*Zh21*Zss21/cb**2*sb
     &      -EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc12*Zh21*Zss22
     &      /cb**2/sb+EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc12
     &      *Zh21*Zss22/cb**2*sb
      AAABR(825) = EE**2/SW**2*Sqrt2*Zdd21*Zh12*Zuu11/cb
     &      +EE**2/SW**2*Sqrt2*Zdd21*Zh11*Zuu11/cb**2*sb
     &      -EE**2/SW**2*Sqrt2*Zdd21*Zh11*Zuu11/cb**2*sb
     &      **3-EE**2/SW**2*Sqrt2*Zdd21*Zh12*Zuu11/cb*sb
     &      **2-2D0*EE*Lambda/MW*Mq/SW*Zdd21*Zh13*Zuu12/cb
     &      **2-2D0*EE*Lambda/MW*Mq/SW*Zdd22*Zh13*Zuu11/cb
     &      **2+EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd21*Zh12
     &      *Zuu11/cb-EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd22
     &      *Zh12*Zuu12/cb+2D0*EE*Lambda/MW*Mq/SW*Zdd21*Zh13
     &      *Zuu12/cb**2*sb**2+2D0*EE*Lambda/MW*Mq/SW*Zdd22
     &      *Zh13*Zuu11/cb**2*sb**2-EE**2/MW**2*Mq**2/SW
     &      **2*Sqrt2*Zdd21*Zh11*Zuu11/cb**2*sb-EE**2/MW
     &      **2*Mq**2/SW**2*Sqrt2*Zdd21*Zh12*Zuu11/cb/sb
     &      **2-EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd22*Zh11
     &      *Zuu12/cb**2/sb+EE**2/MW**2*Mq**2/SW**2*Sqrt2
     &      *Zdd22*Zh11*Zuu12/cb**2*sb
      AAABR(826) = EE**2/SW**2*Sqrt2*Zcc11*Zh12*Zss21/cb
     &      +EE**2/SW**2*Sqrt2*Zcc11*Zh11*Zss21/cb**2*sb
     &      -EE**2/SW**2*Sqrt2*Zcc11*Zh11*Zss21/cb**2*sb
     &      **3-EE**2/SW**2*Sqrt2*Zcc11*Zh12*Zss21/cb*sb
     &      **2-2D0*EE*Lambda/MW*Mc/SW*Zcc12*Zh13*Zss21/cb
     &      **2-2D0*EE*Lambda/MW*Mq/SW*Zcc11*Zh13*Zss22/cb
     &      **2+EE**2/MW**2*Mc**2/SW**2*Sqrt2*Zcc11*Zh12
     &      *Zss21/cb+2D0*EE*Lambda/MW*Mc/SW*Zcc12*Zh13*Zss21
     &      /cb**2*sb**2+2D0*EE*Lambda/MW*Mq/SW*Zcc11*Zh13
     &      *Zss22/cb**2*sb**2-EE**2/MW**2*Mc*Mq/SW**2*Sqrt2
     &      *Zcc12*Zh12*Zss22/cb-EE**2/MW**2*Mc**2/SW**2
     &      *Sqrt2*Zcc11*Zh12*Zss21/cb/sb**2-EE**2/MW**2
     &      *Mq**2/SW**2*Sqrt2*Zcc11*Zh11*Zss21/cb**2*sb
     &      -EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc12*Zh11*Zss22
     &      /cb**2/sb+EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc12
     &      *Zh11*Zss22/cb**2*sb
      AAABR(827) = EE**2/SW**2*Sqrt2*Zdd21*Zh22*Zuu11/cb
     &      +EE**2/SW**2*Sqrt2*Zdd21*Zh21*Zuu11/cb**2*sb
     &      -EE**2/SW**2*Sqrt2*Zdd21*Zh21*Zuu11/cb**2*sb
     &      **3-EE**2/SW**2*Sqrt2*Zdd21*Zh22*Zuu11/cb*sb
     &      **2-2D0*EE*Lambda/MW*Mq/SW*Zdd21*Zh23*Zuu12/cb
     &      **2-2D0*EE*Lambda/MW*Mq/SW*Zdd22*Zh23*Zuu11/cb
     &      **2+EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd21*Zh22
     &      *Zuu11/cb-EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd22
     &      *Zh22*Zuu12/cb+2D0*EE*Lambda/MW*Mq/SW*Zdd21*Zh23
     &      *Zuu12/cb**2*sb**2+2D0*EE*Lambda/MW*Mq/SW*Zdd22
     &      *Zh23*Zuu11/cb**2*sb**2-EE**2/MW**2*Mq**2/SW
     &      **2*Sqrt2*Zdd21*Zh21*Zuu11/cb**2*sb-EE**2/MW
     &      **2*Mq**2/SW**2*Sqrt2*Zdd21*Zh22*Zuu11/cb/sb
     &      **2-EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd22*Zh21
     &      *Zuu12/cb**2/sb+EE**2/MW**2*Mq**2/SW**2*Sqrt2
     &      *Zdd22*Zh21*Zuu12/cb**2*sb
      AAABR(828) = EE**2/SW**2*Sqrt2*Zb21*Zh22*Zt11/cb+EE
     &      **2/SW**2*Sqrt2*Zb21*Zh21*Zt11/cb**2*sb-EE**2
     &      /SW**2*Sqrt2*Zb21*Zh21*Zt11/cb**2*sb**3-EE**2
     &      /SW**2*Sqrt2*Zb21*Zh22*Zt11/cb*sb**2-2D0*EE*Lambda
     &      /MW*MbMM/SW*Zb22*Zh23*Zt11/cb**2-2D0*EE*Lambda
     &      /MW*MtMM/SW*Zb21*Zh23*Zt12/cb**2+EE**2/MW**2
     &      *MtMM**2/SW**2*Sqrt2*Zb21*Zh22*Zt11/cb+2D0*EE
     &      *Lambda/MW*MbMM/SW*Zb22*Zh23*Zt11/cb**2*sb**2
     &      +2D0*EE*Lambda/MW*MtMM/SW*Zb21*Zh23*Zt12/cb**2
     &      *sb**2-EE**2/MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb22
     &      *Zh22*Zt12/cb-EE**2/MW**2*MbMM**2/SW**2*Sqrt2
     &      *Zb21*Zh21*Zt11/cb**2*sb-EE**2/MW**2*MtMM**2
     &      /SW**2*Sqrt2*Zb21*Zh22*Zt11/cb/sb**2-EE**2/MW
     &      **2*MbMM*MtMM/SW**2*Sqrt2*Zb22*Zh21*Zt12/cb**2
     &      /sb+EE**2/MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb22*Zh21
     &      *Zt12/cb**2*sb
      AAABR(829) = EE**2/SW**2*Sqrt2*Zb21*Zh32*Zt11/cb+EE
     &      **2/SW**2*Sqrt2*Zb21*Zh31*Zt11/cb**2*sb-EE**2
     &      /SW**2*Sqrt2*Zb21*Zh31*Zt11/cb**2*sb**3-EE**2
     &      /SW**2*Sqrt2*Zb21*Zh32*Zt11/cb*sb**2-2D0*EE*Lambda
     &      /MW*MbMM/SW*Zb22*Zh33*Zt11/cb**2-2D0*EE*Lambda
     &      /MW*MtMM/SW*Zb21*Zh33*Zt12/cb**2+EE**2/MW**2
     &      *MtMM**2/SW**2*Sqrt2*Zb21*Zh32*Zt11/cb+2D0*EE
     &      *Lambda/MW*MbMM/SW*Zb22*Zh33*Zt11/cb**2*sb**2
     &      +2D0*EE*Lambda/MW*MtMM/SW*Zb21*Zh33*Zt12/cb**2
     &      *sb**2-EE**2/MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb22
     &      *Zh32*Zt12/cb-EE**2/MW**2*MbMM**2/SW**2*Sqrt2
     &      *Zb21*Zh31*Zt11/cb**2*sb-EE**2/MW**2*MtMM**2
     &      /SW**2*Sqrt2*Zb21*Zh32*Zt11/cb/sb**2-EE**2/MW
     &      **2*MbMM*MtMM/SW**2*Sqrt2*Zb22*Zh31*Zt12/cb**2
     &      /sb+EE**2/MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb22*Zh31
     &      *Zt12/cb**2*sb
      AAABR(830) = EE**2/SW**2*Sqrt2*Zcc11*Zh32*Zss21/cb
     &      +EE**2/SW**2*Sqrt2*Zcc11*Zh31*Zss21/cb**2*sb
     &      -EE**2/SW**2*Sqrt2*Zcc11*Zh31*Zss21/cb**2*sb
     &      **3-EE**2/SW**2*Sqrt2*Zcc11*Zh32*Zss21/cb*sb
     &      **2-2D0*EE*Lambda/MW*Mc/SW*Zcc12*Zh33*Zss21/cb
     &      **2-2D0*EE*Lambda/MW*Mq/SW*Zcc11*Zh33*Zss22/cb
     &      **2+EE**2/MW**2*Mc**2/SW**2*Sqrt2*Zcc11*Zh32
     &      *Zss21/cb+2D0*EE*Lambda/MW*Mc/SW*Zcc12*Zh33*Zss21
     &      /cb**2*sb**2+2D0*EE*Lambda/MW*Mq/SW*Zcc11*Zh33
     &      *Zss22/cb**2*sb**2-EE**2/MW**2*Mc*Mq/SW**2*Sqrt2
     &      *Zcc12*Zh32*Zss22/cb-EE**2/MW**2*Mc**2/SW**2
     &      *Sqrt2*Zcc11*Zh32*Zss21/cb/sb**2-EE**2/MW**2
     &      *Mq**2/SW**2*Sqrt2*Zcc11*Zh31*Zss21/cb**2*sb
     &      -EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc12*Zh31*Zss22
     &      /cb**2/sb+EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc12
     &      *Zh31*Zss22/cb**2*sb
      AAABR(831) = EE**2/SW**2*Sqrt2*Zb21*Zh12*Zt11/cb+EE
     &      **2/SW**2*Sqrt2*Zb21*Zh11*Zt11/cb**2*sb-EE**2
     &      /SW**2*Sqrt2*Zb21*Zh11*Zt11/cb**2*sb**3-EE**2
     &      /SW**2*Sqrt2*Zb21*Zh12*Zt11/cb*sb**2-2D0*EE*Lambda
     &      /MW*MbMM/SW*Zb22*Zh13*Zt11/cb**2-2D0*EE*Lambda
     &      /MW*MtMM/SW*Zb21*Zh13*Zt12/cb**2+EE**2/MW**2
     &      *MtMM**2/SW**2*Sqrt2*Zb21*Zh12*Zt11/cb+2D0*EE
     &      *Lambda/MW*MbMM/SW*Zb22*Zh13*Zt11/cb**2*sb**2
     &      +2D0*EE*Lambda/MW*MtMM/SW*Zb21*Zh13*Zt12/cb**2
     &      *sb**2-EE**2/MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb22
     &      *Zh12*Zt12/cb-EE**2/MW**2*MbMM**2/SW**2*Sqrt2
     &      *Zb21*Zh11*Zt11/cb**2*sb-EE**2/MW**2*MtMM**2
     &      /SW**2*Sqrt2*Zb21*Zh12*Zt11/cb/sb**2-EE**2/MW
     &      **2*MbMM*MtMM/SW**2*Sqrt2*Zb22*Zh11*Zt12/cb**2
     &      /sb+EE**2/MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb22*Zh11
     &      *Zt12/cb**2*sb
      AAABR(832) = EE**2/SW**2*Sqrt2*Zdd21*Zh32*Zuu11/cb
     &      +EE**2/SW**2*Sqrt2*Zdd21*Zh31*Zuu11/cb**2*sb
     &      -EE**2/SW**2*Sqrt2*Zdd21*Zh31*Zuu11/cb**2*sb
     &      **3-EE**2/SW**2*Sqrt2*Zdd21*Zh32*Zuu11/cb*sb
     &      **2-2D0*EE*Lambda/MW*Mq/SW*Zdd21*Zh33*Zuu12/cb
     &      **2-2D0*EE*Lambda/MW*Mq/SW*Zdd22*Zh33*Zuu11/cb
     &      **2+EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd21*Zh32
     &      *Zuu11/cb-EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd22
     &      *Zh32*Zuu12/cb+2D0*EE*Lambda/MW*Mq/SW*Zdd21*Zh33
     &      *Zuu12/cb**2*sb**2+2D0*EE*Lambda/MW*Mq/SW*Zdd22
     &      *Zh33*Zuu11/cb**2*sb**2-EE**2/MW**2*Mq**2/SW
     &      **2*Sqrt2*Zdd21*Zh31*Zuu11/cb**2*sb-EE**2/MW
     &      **2*Mq**2/SW**2*Sqrt2*Zdd21*Zh32*Zuu11/cb/sb
     &      **2-EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd22*Zh31
     &      *Zuu12/cb**2/sb+EE**2/MW**2*Mq**2/SW**2*Sqrt2
     &      *Zdd22*Zh31*Zuu12/cb**2*sb
      AAABR(833) = EE**2/cb*(Zcc21**2/CW**2*sb-Zcc21**2/CW
     &      **2*sb**3-4D0*Zcc22**2/CW**2*sb+4D0*Zcc22**2
     &      /CW**2*sb**3+3D0*Zcc21**2/SW**2*sb-3D0*Zcc21
     &      **2/SW**2*sb**3+3D0*Mc**2/MW**2/SW**2*Zcc22**2
     &      /sb-3D0*Mc**2/MW**2/SW**2*Zcc22**2*sb-3D0*Mq
     &      **2/MW**2/SW**2*Zcc21**2*sb)
      AAABR(834) = EE**2/cb*(Zuu21**2/CW**2*sb-Zuu21**2/CW
     &      **2*sb**3-4D0*Zuu22**2/CW**2*sb+4D0*Zuu22**2
     &      /CW**2*sb**3+3D0*Zuu21**2/SW**2*sb-3D0*Zuu21
     &      **2/SW**2*sb**3-3D0*Mq**2/MW**2/SW**2*sb+3D0
     &      *Mq**2/MW**2/SW**2*Zuu22**2/sb)
      AAABR(835) = EE**2/cb*(Zt21**2/CW**2*sb-Zt21**2/CW
     &      **2*sb**3-4D0*Zt22**2/CW**2*sb+4D0*Zt22**2/CW
     &      **2*sb**3+3D0*Zt21**2/SW**2*sb-3D0*Zt21**2/SW
     &      **2*sb**3-3D0*MbMM**2/MW**2/SW**2*Zt21**2*sb
     &      +3D0*MtMM**2/MW**2/SW**2*Zt22**2/sb-3D0*MtMM
     &      **2/MW**2/SW**2*Zt22**2*sb)
      AAABR(836) = EE**2/SW**2*Sqrt2/cb*(2D0*Zcc21*Zss11
     &      *sb-2D0*Zcc21*Zss11*sb**3-Mc**2/MW**2*Zcc21*Zss11
     &      /sb+Mc**2/MW**2*Zcc21*Zss11*sb-Mq**2/MW**2*Zcc21
     &      *Zss11*sb+Mc/MW**2*Mq*Zcc22*Zss12/sb)
      AAABR(837) = EE**2/SW**2*Sqrt2/cb*(2D0*Zdd11*Zuu21
     &      *sb-2D0*Zdd11*Zuu21*sb**3-Mq**2/MW**2*Zdd11*Zuu21
     &      /sb+Mq**2/MW**2*Zdd12*Zuu22/sb)
      AAABR(838) = EE**2/SW**2*Sqrt2/cb*(2D0*Zb11*Zt21*sb
     &      -2D0*Zb11*Zt21*sb**3-MbMM**2/MW**2*Zb11*Zt21
     &      *sb-MtMM**2/MW**2*Zb11*Zt21/sb+MtMM**2/MW**2
     &      *Zb11*Zt21*sb+MbMM/MW**2*MtMM*Zb12*Zt22/sb)
      AAABR(839) = EE**2/SW**2*Sqrt2/cb*(2D0*Zdd21*Zuu21
     &      *sb-2D0*Zdd21*Zuu21*sb**3-Mq**2/MW**2*Zdd21*Zuu21
     &      /sb+Mq**2/MW**2*Zdd22*Zuu22/sb)
      AAABR(840) = EE**2/SW**2*Sqrt2/cb*(2D0*Zb21*Zt21*sb
     &      -2D0*Zb21*Zt21*sb**3-MbMM**2/MW**2*Zb21*Zt21
     &      *sb-MtMM**2/MW**2*Zb21*Zt21/sb+MtMM**2/MW**2
     &      *Zb21*Zt21*sb+MbMM/MW**2*MtMM*Zb22*Zt22/sb)
      AAABR(841) = EE**2/SW**2*Sqrt2/cb*(2D0*Zcc21*Zss21
     &      *sb-2D0*Zcc21*Zss21*sb**3-Mc**2/MW**2*Zcc21*Zss21
     &      /sb+Mc**2/MW**2*Zcc21*Zss21*sb-Mq**2/MW**2*Zcc21
     &      *Zss21*sb+Mc/MW**2*Mq*Zcc22*Zss22/sb)
      AAABR(842) = 1D0/cb**2*(EE**2*Pa11/SW**2*Sqrt2*Zb11
     &      *Zt21-3D0*EE**2*Pa11/SW**2*Sqrt2*Zb11*Zt21*sb
     &      **2+2D0*EE**2*Pa11/SW**2*Sqrt2*Zb11*Zt21*sb**4
     &      +2D0*EE*Lambda/MW*MbMM*Pa12/SW*Zb12*Zt21-2D0
     &      *EE*Lambda/MW*MtMM*Pa12/SW*Zb11*Zt22+2D0*EE**2
     &      /MW**2*MtMM**2*Pa11/SW**2*Sqrt2*Zb11*Zt21-2D0
     &      *EE*Lambda/MW*MbMM*Pa12/SW*Zb12*Zt21*sb**2+2D0
     &      *EE*Lambda/MW*MtMM*Pa12/SW*Zb11*Zt22*sb**2+EE
     &      **2/MW**2*MbMM**2*Pa11/SW**2*Sqrt2*Zb11*Zt21
     &      *sb**2-EE**2/MW**2*MtMM**2*Pa11/SW**2*Sqrt2*Zb11
     &      *Zt21/sb**2-EE**2/MW**2*MtMM**2*Pa11/SW**2*Sqrt2
     &      *Zb11*Zt21*sb**2)
      AAABR(843) = 1D0/cb**2*(EE**2*Pa21/SW**2*Sqrt2*Zb11
     &      *Zt21-3D0*EE**2*Pa21/SW**2*Sqrt2*Zb11*Zt21*sb
     &      **2+2D0*EE**2*Pa21/SW**2*Sqrt2*Zb11*Zt21*sb**4
     &      +2D0*EE*Lambda/MW*MbMM*Pa22/SW*Zb12*Zt21-2D0
     &      *EE*Lambda/MW*MtMM*Pa22/SW*Zb11*Zt22+2D0*EE**2
     &      /MW**2*MtMM**2*Pa21/SW**2*Sqrt2*Zb11*Zt21-2D0
     &      *EE*Lambda/MW*MbMM*Pa22/SW*Zb12*Zt21*sb**2+2D0
     &      *EE*Lambda/MW*MtMM*Pa22/SW*Zb11*Zt22*sb**2+EE
     &      **2/MW**2*MbMM**2*Pa21/SW**2*Sqrt2*Zb11*Zt21
     &      *sb**2-EE**2/MW**2*MtMM**2*Pa21/SW**2*Sqrt2*Zb11
     &      *Zt21/sb**2-EE**2/MW**2*MtMM**2*Pa21/SW**2*Sqrt2
     &      *Zb11*Zt21*sb**2)
      AAABR(844) = 1D0/cb**2*(EE**2*Pa11/SW**2*Sqrt2*Zcc21
     &      *Zss11-3D0*EE**2*Pa11/SW**2*Sqrt2*Zcc21*Zss11
     &      *sb**2+2D0*EE**2*Pa11/SW**2*Sqrt2*Zcc21*Zss11
     &      *sb**4-2D0*EE*Lambda/MW*Mc*Pa12/SW*Zcc22*Zss11
     &      +2D0*EE*Lambda/MW*Mq*Pa12/SW*Zcc21*Zss12+2D0
     &      *EE**2/MW**2*Mc**2*Pa11/SW**2*Sqrt2*Zcc21*Zss11
     &      +2D0*EE*Lambda/MW*Mc*Pa12/SW*Zcc22*Zss11*sb**2
     &      -2D0*EE*Lambda/MW*Mq*Pa12/SW*Zcc21*Zss12*sb**2
     &      -EE**2/MW**2*Mc**2*Pa11/SW**2*Sqrt2*Zcc21*Zss11
     &      /sb**2-EE**2/MW**2*Mc**2*Pa11/SW**2*Sqrt2*Zcc21
     &      *Zss11*sb**2+EE**2/MW**2*Mq**2*Pa11/SW**2*Sqrt2
     &      *Zcc21*Zss11*sb**2)
      AAABR(845) = 1D0/cb**2*(EE**2*Pa11/SW**2*Sqrt2*Zdd11
     &      *Zuu21-3D0*EE**2*Pa11/SW**2*Sqrt2*Zdd11*Zuu21
     &      *sb**2+2D0*EE**2*Pa11/SW**2*Sqrt2*Zdd11*Zuu21
     &      *sb**4-2D0*EE*Lambda/MW*Mq*Pa12/SW*Zdd11*Zuu22
     &      +2D0*EE*Lambda/MW*Mq*Pa12/SW*Zdd12*Zuu21+2D0
     &      *EE**2/MW**2*Mq**2*Pa11/SW**2*Sqrt2*Zdd11*Zuu21
     &      +2D0*EE*Lambda/MW*Mq*Pa12/SW*Zdd11*Zuu22*sb**2
     &      -2D0*EE*Lambda/MW*Mq*Pa12/SW*Zdd12*Zuu21*sb**2
     &      -EE**2/MW**2*Mq**2*Pa11/SW**2*Sqrt2*Zdd11*Zuu21
     &      /sb**2)
      AAABR(846) = 1D0/cb**2*(EE**2*Pa21/SW**2*Sqrt2*Zdd11
     &      *Zuu21-3D0*EE**2*Pa21/SW**2*Sqrt2*Zdd11*Zuu21
     &      *sb**2+2D0*EE**2*Pa21/SW**2*Sqrt2*Zdd11*Zuu21
     &      *sb**4-2D0*EE*Lambda/MW*Mq*Pa22/SW*Zdd11*Zuu22
     &      +2D0*EE*Lambda/MW*Mq*Pa22/SW*Zdd12*Zuu21+2D0
     &      *EE**2/MW**2*Mq**2*Pa21/SW**2*Sqrt2*Zdd11*Zuu21
     &      +2D0*EE*Lambda/MW*Mq*Pa22/SW*Zdd11*Zuu22*sb**2
     &      -2D0*EE*Lambda/MW*Mq*Pa22/SW*Zdd12*Zuu21*sb**2
     &      -EE**2/MW**2*Mq**2*Pa21/SW**2*Sqrt2*Zdd11*Zuu21
     &      /sb**2)
      AAABR(847) = 1D0/cb**2*(EE**2*Pa21/SW**2*Sqrt2*Zcc21
     &      *Zss11-3D0*EE**2*Pa21/SW**2*Sqrt2*Zcc21*Zss11
     &      *sb**2+2D0*EE**2*Pa21/SW**2*Sqrt2*Zcc21*Zss11
     &      *sb**4-2D0*EE*Lambda/MW*Mc*Pa22/SW*Zcc22*Zss11
     &      +2D0*EE*Lambda/MW*Mq*Pa22/SW*Zcc21*Zss12+2D0
     &      *EE**2/MW**2*Mc**2*Pa21/SW**2*Sqrt2*Zcc21*Zss11
     &      +2D0*EE*Lambda/MW*Mc*Pa22/SW*Zcc22*Zss11*sb**2
     &      -2D0*EE*Lambda/MW*Mq*Pa22/SW*Zcc21*Zss12*sb**2
     &      -EE**2/MW**2*Mc**2*Pa21/SW**2*Sqrt2*Zcc21*Zss11
     &      /sb**2-EE**2/MW**2*Mc**2*Pa21/SW**2*Sqrt2*Zcc21
     &      *Zss11*sb**2+EE**2/MW**2*Mq**2*Pa21/SW**2*Sqrt2
     &      *Zcc21*Zss11*sb**2)
      AAABR(848) = 1D0/cb**2*(EE**2*Pa11/SW**2*Sqrt2*Zdd21
     &      *Zuu21-3D0*EE**2*Pa11/SW**2*Sqrt2*Zdd21*Zuu21
     &      *sb**2+2D0*EE**2*Pa11/SW**2*Sqrt2*Zdd21*Zuu21
     &      *sb**4-2D0*EE*Lambda/MW*Mq*Pa12/SW*Zdd21*Zuu22
     &      +2D0*EE*Lambda/MW*Mq*Pa12/SW*Zdd22*Zuu21+2D0
     &      *EE**2/MW**2*Mq**2*Pa11/SW**2*Sqrt2*Zdd21*Zuu21
     &      +2D0*EE*Lambda/MW*Mq*Pa12/SW*Zdd21*Zuu22*sb**2
     &      -2D0*EE*Lambda/MW*Mq*Pa12/SW*Zdd22*Zuu21*sb**2
     &      -EE**2/MW**2*Mq**2*Pa11/SW**2*Sqrt2*Zdd21*Zuu21
     &      /sb**2)
      AAABR(849) = 1D0/cb**2*(EE**2*Pa21/SW**2*Sqrt2*Zdd21
     &      *Zuu21-3D0*EE**2*Pa21/SW**2*Sqrt2*Zdd21*Zuu21
     &      *sb**2+2D0*EE**2*Pa21/SW**2*Sqrt2*Zdd21*Zuu21
     &      *sb**4-2D0*EE*Lambda/MW*Mq*Pa22/SW*Zdd21*Zuu22
     &      +2D0*EE*Lambda/MW*Mq*Pa22/SW*Zdd22*Zuu21+2D0
     &      *EE**2/MW**2*Mq**2*Pa21/SW**2*Sqrt2*Zdd21*Zuu21
     &      +2D0*EE*Lambda/MW*Mq*Pa22/SW*Zdd21*Zuu22*sb**2
     &      -2D0*EE*Lambda/MW*Mq*Pa22/SW*Zdd22*Zuu21*sb**2
     &      -EE**2/MW**2*Mq**2*Pa21/SW**2*Sqrt2*Zdd21*Zuu21
     &      /sb**2)
      AAABR(850) = 1D0/cb**2*(EE**2*Pa21/SW**2*Sqrt2*Zb21
     &      *Zt21-3D0*EE**2*Pa21/SW**2*Sqrt2*Zb21*Zt21*sb
     &      **2+2D0*EE**2*Pa21/SW**2*Sqrt2*Zb21*Zt21*sb**4
     &      +2D0*EE*Lambda/MW*MbMM*Pa22/SW*Zb22*Zt21-2D0
     &      *EE*Lambda/MW*MtMM*Pa22/SW*Zb21*Zt22+2D0*EE**2
     &      /MW**2*MtMM**2*Pa21/SW**2*Sqrt2*Zb21*Zt21-2D0
     &      *EE*Lambda/MW*MbMM*Pa22/SW*Zb22*Zt21*sb**2+2D0
     &      *EE*Lambda/MW*MtMM*Pa22/SW*Zb21*Zt22*sb**2+EE
     &      **2/MW**2*MbMM**2*Pa21/SW**2*Sqrt2*Zb21*Zt21
     &      *sb**2-EE**2/MW**2*MtMM**2*Pa21/SW**2*Sqrt2*Zb21
     &      *Zt21/sb**2-EE**2/MW**2*MtMM**2*Pa21/SW**2*Sqrt2
     &      *Zb21*Zt21*sb**2)
      AAABR(851) = 1D0/cb**2*(EE**2*Pa11/SW**2*Sqrt2*Zb21
     &      *Zt21-3D0*EE**2*Pa11/SW**2*Sqrt2*Zb21*Zt21*sb
     &      **2+2D0*EE**2*Pa11/SW**2*Sqrt2*Zb21*Zt21*sb**4
     &      +2D0*EE*Lambda/MW*MbMM*Pa12/SW*Zb22*Zt21-2D0
     &      *EE*Lambda/MW*MtMM*Pa12/SW*Zb21*Zt22+2D0*EE**2
     &      /MW**2*MtMM**2*Pa11/SW**2*Sqrt2*Zb21*Zt21-2D0
     &      *EE*Lambda/MW*MbMM*Pa12/SW*Zb22*Zt21*sb**2+2D0
     &      *EE*Lambda/MW*MtMM*Pa12/SW*Zb21*Zt22*sb**2+EE
     &      **2/MW**2*MbMM**2*Pa11/SW**2*Sqrt2*Zb21*Zt21
     &      *sb**2-EE**2/MW**2*MtMM**2*Pa11/SW**2*Sqrt2*Zb21
     &      *Zt21/sb**2-EE**2/MW**2*MtMM**2*Pa11/SW**2*Sqrt2
     &      *Zb21*Zt21*sb**2)
      AAABR(852) = 1D0/cb**2*(EE**2*Pa11/SW**2*Sqrt2*Zcc21
     &      *Zss21-3D0*EE**2*Pa11/SW**2*Sqrt2*Zcc21*Zss21
     &      *sb**2+2D0*EE**2*Pa11/SW**2*Sqrt2*Zcc21*Zss21
     &      *sb**4-2D0*EE*Lambda/MW*Mc*Pa12/SW*Zcc22*Zss21
     &      +2D0*EE*Lambda/MW*Mq*Pa12/SW*Zcc21*Zss22+2D0
     &      *EE**2/MW**2*Mc**2*Pa11/SW**2*Sqrt2*Zcc21*Zss21
     &      +2D0*EE*Lambda/MW*Mc*Pa12/SW*Zcc22*Zss21*sb**2
     &      -2D0*EE*Lambda/MW*Mq*Pa12/SW*Zcc21*Zss22*sb**2
     &      -EE**2/MW**2*Mc**2*Pa11/SW**2*Sqrt2*Zcc21*Zss21
     &      /sb**2-EE**2/MW**2*Mc**2*Pa11/SW**2*Sqrt2*Zcc21
     &      *Zss21*sb**2+EE**2/MW**2*Mq**2*Pa11/SW**2*Sqrt2
     &      *Zcc21*Zss21*sb**2)
      AAABR(853) = 1D0/cb**2*(EE**2*Pa21/SW**2*Sqrt2*Zcc21
     &      *Zss21-3D0*EE**2*Pa21/SW**2*Sqrt2*Zcc21*Zss21
     &      *sb**2+2D0*EE**2*Pa21/SW**2*Sqrt2*Zcc21*Zss21
     &      *sb**4-2D0*EE*Lambda/MW*Mc*Pa22/SW*Zcc22*Zss21
     &      +2D0*EE*Lambda/MW*Mq*Pa22/SW*Zcc21*Zss22+2D0
     &      *EE**2/MW**2*Mc**2*Pa21/SW**2*Sqrt2*Zcc21*Zss21
     &      +2D0*EE*Lambda/MW*Mc*Pa22/SW*Zcc22*Zss21*sb**2
     &      -2D0*EE*Lambda/MW*Mq*Pa22/SW*Zcc21*Zss22*sb**2
     &      -EE**2/MW**2*Mc**2*Pa21/SW**2*Sqrt2*Zcc21*Zss21
     &      /sb**2-EE**2/MW**2*Mc**2*Pa21/SW**2*Sqrt2*Zcc21
     &      *Zss21*sb**2+EE**2/MW**2*Mq**2*Pa21/SW**2*Sqrt2
     &      *Zcc21*Zss21*sb**2)
      AAABR(854) = EE**2/SW**2*Sqrt2*Zcc21*Zh32*Zss11/cb
     &      +EE**2/SW**2*Sqrt2*Zcc21*Zh31*Zss11/cb**2*sb
     &      -EE**2/SW**2*Sqrt2*Zcc21*Zh31*Zss11/cb**2*sb
     &      **3-EE**2/SW**2*Sqrt2*Zcc21*Zh32*Zss11/cb*sb
     &      **2-2D0*EE*Lambda/MW*Mc/SW*Zcc22*Zh33*Zss11/cb
     &      **2-2D0*EE*Lambda/MW*Mq/SW*Zcc21*Zh33*Zss12/cb
     &      **2+EE**2/MW**2*Mc**2/SW**2*Sqrt2*Zcc21*Zh32
     &      *Zss11/cb+2D0*EE*Lambda/MW*Mc/SW*Zcc22*Zh33*Zss11
     &      /cb**2*sb**2+2D0*EE*Lambda/MW*Mq/SW*Zcc21*Zh33
     &      *Zss12/cb**2*sb**2-EE**2/MW**2*Mc*Mq/SW**2*Sqrt2
     &      *Zcc22*Zh32*Zss12/cb-EE**2/MW**2*Mc**2/SW**2
     &      *Sqrt2*Zcc21*Zh32*Zss11/cb/sb**2-EE**2/MW**2
     &      *Mq**2/SW**2*Sqrt2*Zcc21*Zh31*Zss11/cb**2*sb
     &      -EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc22*Zh31*Zss12
     &      /cb**2/sb+EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc22
     &      *Zh31*Zss12/cb**2*sb
      AAABR(855) = EE**2/SW**2*Sqrt2*Zdd11*Zh22*Zuu21/cb
     &      +EE**2/SW**2*Sqrt2*Zdd11*Zh21*Zuu21/cb**2*sb
     &      -EE**2/SW**2*Sqrt2*Zdd11*Zh21*Zuu21/cb**2*sb
     &      **3-EE**2/SW**2*Sqrt2*Zdd11*Zh22*Zuu21/cb*sb
     &      **2-2D0*EE*Lambda/MW*Mq/SW*Zdd11*Zh23*Zuu22/cb
     &      **2-2D0*EE*Lambda/MW*Mq/SW*Zdd12*Zh23*Zuu21/cb
     &      **2+EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd11*Zh22
     &      *Zuu21/cb-EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd12
     &      *Zh22*Zuu22/cb+2D0*EE*Lambda/MW*Mq/SW*Zdd11*Zh23
     &      *Zuu22/cb**2*sb**2+2D0*EE*Lambda/MW*Mq/SW*Zdd12
     &      *Zh23*Zuu21/cb**2*sb**2-EE**2/MW**2*Mq**2/SW
     &      **2*Sqrt2*Zdd11*Zh21*Zuu21/cb**2*sb-EE**2/MW
     &      **2*Mq**2/SW**2*Sqrt2*Zdd11*Zh22*Zuu21/cb/sb
     &      **2-EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd12*Zh21
     &      *Zuu22/cb**2/sb+EE**2/MW**2*Mq**2/SW**2*Sqrt2
     &      *Zdd12*Zh21*Zuu22/cb**2*sb
      AAABR(856) = EE**2/SW**2*Sqrt2*Zb11*Zh32*Zt21/cb+EE
     &      **2/SW**2*Sqrt2*Zb11*Zh31*Zt21/cb**2*sb-EE**2
     &      /SW**2*Sqrt2*Zb11*Zh31*Zt21/cb**2*sb**3-EE**2
     &      /SW**2*Sqrt2*Zb11*Zh32*Zt21/cb*sb**2-2D0*EE*Lambda
     &      /MW*MbMM/SW*Zb12*Zh33*Zt21/cb**2-2D0*EE*Lambda
     &      /MW*MtMM/SW*Zb11*Zh33*Zt22/cb**2+EE**2/MW**2
     &      *MtMM**2/SW**2*Sqrt2*Zb11*Zh32*Zt21/cb+2D0*EE
     &      *Lambda/MW*MbMM/SW*Zb12*Zh33*Zt21/cb**2*sb**2
     &      +2D0*EE*Lambda/MW*MtMM/SW*Zb11*Zh33*Zt22/cb**2
     &      *sb**2-EE**2/MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb12
     &      *Zh32*Zt22/cb-EE**2/MW**2*MbMM**2/SW**2*Sqrt2
     &      *Zb11*Zh31*Zt21/cb**2*sb-EE**2/MW**2*MtMM**2
     &      /SW**2*Sqrt2*Zb11*Zh32*Zt21/cb/sb**2-EE**2/MW
     &      **2*MbMM*MtMM/SW**2*Sqrt2*Zb12*Zh31*Zt22/cb**2
     &      /sb+EE**2/MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb12*Zh31
     &      *Zt22/cb**2*sb
      AAABR(857) = EE**2/SW**2*Sqrt2*Zcc21*Zh22*Zss11/cb
     &      +EE**2/SW**2*Sqrt2*Zcc21*Zh21*Zss11/cb**2*sb
     &      -EE**2/SW**2*Sqrt2*Zcc21*Zh21*Zss11/cb**2*sb
     &      **3-EE**2/SW**2*Sqrt2*Zcc21*Zh22*Zss11/cb*sb
     &      **2-2D0*EE*Lambda/MW*Mc/SW*Zcc22*Zh23*Zss11/cb
     &      **2-2D0*EE*Lambda/MW*Mq/SW*Zcc21*Zh23*Zss12/cb
     &      **2+EE**2/MW**2*Mc**2/SW**2*Sqrt2*Zcc21*Zh22
     &      *Zss11/cb+2D0*EE*Lambda/MW*Mc/SW*Zcc22*Zh23*Zss11
     &      /cb**2*sb**2+2D0*EE*Lambda/MW*Mq/SW*Zcc21*Zh23
     &      *Zss12/cb**2*sb**2-EE**2/MW**2*Mc*Mq/SW**2*Sqrt2
     &      *Zcc22*Zh22*Zss12/cb-EE**2/MW**2*Mc**2/SW**2
     &      *Sqrt2*Zcc21*Zh22*Zss11/cb/sb**2-EE**2/MW**2
     &      *Mq**2/SW**2*Sqrt2*Zcc21*Zh21*Zss11/cb**2*sb
     &      -EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc22*Zh21*Zss12
     &      /cb**2/sb+EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc22
     &      *Zh21*Zss12/cb**2*sb
      AAABR(858) = EE**2/SW**2*Sqrt2*Zdd11*Zh32*Zuu21/cb
     &      +EE**2/SW**2*Sqrt2*Zdd11*Zh31*Zuu21/cb**2*sb
     &      -EE**2/SW**2*Sqrt2*Zdd11*Zh31*Zuu21/cb**2*sb
     &      **3-EE**2/SW**2*Sqrt2*Zdd11*Zh32*Zuu21/cb*sb
     &      **2-2D0*EE*Lambda/MW*Mq/SW*Zdd11*Zh33*Zuu22/cb
     &      **2-2D0*EE*Lambda/MW*Mq/SW*Zdd12*Zh33*Zuu21/cb
     &      **2+EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd11*Zh32
     &      *Zuu21/cb-EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd12
     &      *Zh32*Zuu22/cb+2D0*EE*Lambda/MW*Mq/SW*Zdd11*Zh33
     &      *Zuu22/cb**2*sb**2+2D0*EE*Lambda/MW*Mq/SW*Zdd12
     &      *Zh33*Zuu21/cb**2*sb**2-EE**2/MW**2*Mq**2/SW
     &      **2*Sqrt2*Zdd11*Zh31*Zuu21/cb**2*sb-EE**2/MW
     &      **2*Mq**2/SW**2*Sqrt2*Zdd11*Zh32*Zuu21/cb/sb
     &      **2-EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd12*Zh31
     &      *Zuu22/cb**2/sb+EE**2/MW**2*Mq**2/SW**2*Sqrt2
     &      *Zdd12*Zh31*Zuu22/cb**2*sb
      AAABR(859) = EE**2/SW**2*Sqrt2*Zb11*Zh22*Zt21/cb+EE
     &      **2/SW**2*Sqrt2*Zb11*Zh21*Zt21/cb**2*sb-EE**2
     &      /SW**2*Sqrt2*Zb11*Zh21*Zt21/cb**2*sb**3-EE**2
     &      /SW**2*Sqrt2*Zb11*Zh22*Zt21/cb*sb**2-2D0*EE*Lambda
     &      /MW*MbMM/SW*Zb12*Zh23*Zt21/cb**2-2D0*EE*Lambda
     &      /MW*MtMM/SW*Zb11*Zh23*Zt22/cb**2+EE**2/MW**2
     &      *MtMM**2/SW**2*Sqrt2*Zb11*Zh22*Zt21/cb+2D0*EE
     &      *Lambda/MW*MbMM/SW*Zb12*Zh23*Zt21/cb**2*sb**2
     &      +2D0*EE*Lambda/MW*MtMM/SW*Zb11*Zh23*Zt22/cb**2
     &      *sb**2-EE**2/MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb12
     &      *Zh22*Zt22/cb-EE**2/MW**2*MbMM**2/SW**2*Sqrt2
     &      *Zb11*Zh21*Zt21/cb**2*sb-EE**2/MW**2*MtMM**2
     &      /SW**2*Sqrt2*Zb11*Zh22*Zt21/cb/sb**2-EE**2/MW
     &      **2*MbMM*MtMM/SW**2*Sqrt2*Zb12*Zh21*Zt22/cb**2
     &      /sb+EE**2/MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb12*Zh21
     &      *Zt22/cb**2*sb
      AAABR(860) = EE**2/SW**2*Sqrt2*Zb11*Zh12*Zt21/cb+EE
     &      **2/SW**2*Sqrt2*Zb11*Zh11*Zt21/cb**2*sb-EE**2
     &      /SW**2*Sqrt2*Zb11*Zh11*Zt21/cb**2*sb**3-EE**2
     &      /SW**2*Sqrt2*Zb11*Zh12*Zt21/cb*sb**2-2D0*EE*Lambda
     &      /MW*MbMM/SW*Zb12*Zh13*Zt21/cb**2-2D0*EE*Lambda
     &      /MW*MtMM/SW*Zb11*Zh13*Zt22/cb**2+EE**2/MW**2
     &      *MtMM**2/SW**2*Sqrt2*Zb11*Zh12*Zt21/cb+2D0*EE
     &      *Lambda/MW*MbMM/SW*Zb12*Zh13*Zt21/cb**2*sb**2
     &      +2D0*EE*Lambda/MW*MtMM/SW*Zb11*Zh13*Zt22/cb**2
     &      *sb**2-EE**2/MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb12
     &      *Zh12*Zt22/cb-EE**2/MW**2*MbMM**2/SW**2*Sqrt2
     &      *Zb11*Zh11*Zt21/cb**2*sb-EE**2/MW**2*MtMM**2
     &      /SW**2*Sqrt2*Zb11*Zh12*Zt21/cb/sb**2-EE**2/MW
     &      **2*MbMM*MtMM/SW**2*Sqrt2*Zb12*Zh11*Zt22/cb**2
     &      /sb+EE**2/MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb12*Zh11
     &      *Zt22/cb**2*sb
      AAABR(861) = EE**2/SW**2*Sqrt2*Zdd11*Zh12*Zuu21/cb
     &      +EE**2/SW**2*Sqrt2*Zdd11*Zh11*Zuu21/cb**2*sb
     &      -EE**2/SW**2*Sqrt2*Zdd11*Zh11*Zuu21/cb**2*sb
     &      **3-EE**2/SW**2*Sqrt2*Zdd11*Zh12*Zuu21/cb*sb
     &      **2-2D0*EE*Lambda/MW*Mq/SW*Zdd11*Zh13*Zuu22/cb
     &      **2-2D0*EE*Lambda/MW*Mq/SW*Zdd12*Zh13*Zuu21/cb
     &      **2+EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd11*Zh12
     &      *Zuu21/cb-EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd12
     &      *Zh12*Zuu22/cb+2D0*EE*Lambda/MW*Mq/SW*Zdd11*Zh13
     &      *Zuu22/cb**2*sb**2+2D0*EE*Lambda/MW*Mq/SW*Zdd12
     &      *Zh13*Zuu21/cb**2*sb**2-EE**2/MW**2*Mq**2/SW
     &      **2*Sqrt2*Zdd11*Zh11*Zuu21/cb**2*sb-EE**2/MW
     &      **2*Mq**2/SW**2*Sqrt2*Zdd11*Zh12*Zuu21/cb/sb
     &      **2-EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd12*Zh11
     &      *Zuu22/cb**2/sb+EE**2/MW**2*Mq**2/SW**2*Sqrt2
     &      *Zdd12*Zh11*Zuu22/cb**2*sb
      AAABR(862) = EE**2/SW**2*Sqrt2*Zcc21*Zh12*Zss11/cb
     &      +EE**2/SW**2*Sqrt2*Zcc21*Zh11*Zss11/cb**2*sb
     &      -EE**2/SW**2*Sqrt2*Zcc21*Zh11*Zss11/cb**2*sb
     &      **3-EE**2/SW**2*Sqrt2*Zcc21*Zh12*Zss11/cb*sb
     &      **2-2D0*EE*Lambda/MW*Mc/SW*Zcc22*Zh13*Zss11/cb
     &      **2-2D0*EE*Lambda/MW*Mq/SW*Zcc21*Zh13*Zss12/cb
     &      **2+EE**2/MW**2*Mc**2/SW**2*Sqrt2*Zcc21*Zh12
     &      *Zss11/cb+2D0*EE*Lambda/MW*Mc/SW*Zcc22*Zh13*Zss11
     &      /cb**2*sb**2+2D0*EE*Lambda/MW*Mq/SW*Zcc21*Zh13
     &      *Zss12/cb**2*sb**2-EE**2/MW**2*Mc*Mq/SW**2*Sqrt2
     &      *Zcc22*Zh12*Zss12/cb-EE**2/MW**2*Mc**2/SW**2
     &      *Sqrt2*Zcc21*Zh12*Zss11/cb/sb**2-EE**2/MW**2
     &      *Mq**2/SW**2*Sqrt2*Zcc21*Zh11*Zss11/cb**2*sb
     &      -EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc22*Zh11*Zss12
     &      /cb**2/sb+EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc22
     &      *Zh11*Zss12/cb**2*sb
      AAABR(863) = EE**2/SW**2*Sqrt2*Zdd21*Zh12*Zuu21/cb
     &      +EE**2/SW**2*Sqrt2*Zdd21*Zh11*Zuu21/cb**2*sb
     &      -EE**2/SW**2*Sqrt2*Zdd21*Zh11*Zuu21/cb**2*sb
     &      **3-EE**2/SW**2*Sqrt2*Zdd21*Zh12*Zuu21/cb*sb
     &      **2-2D0*EE*Lambda/MW*Mq/SW*Zdd21*Zh13*Zuu22/cb
     &      **2-2D0*EE*Lambda/MW*Mq/SW*Zdd22*Zh13*Zuu21/cb
     &      **2+EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd21*Zh12
     &      *Zuu21/cb-EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd22
     &      *Zh12*Zuu22/cb+2D0*EE*Lambda/MW*Mq/SW*Zdd21*Zh13
     &      *Zuu22/cb**2*sb**2+2D0*EE*Lambda/MW*Mq/SW*Zdd22
     &      *Zh13*Zuu21/cb**2*sb**2-EE**2/MW**2*Mq**2/SW
     &      **2*Sqrt2*Zdd21*Zh11*Zuu21/cb**2*sb-EE**2/MW
     &      **2*Mq**2/SW**2*Sqrt2*Zdd21*Zh12*Zuu21/cb/sb
     &      **2-EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd22*Zh11
     &      *Zuu22/cb**2/sb+EE**2/MW**2*Mq**2/SW**2*Sqrt2
     &      *Zdd22*Zh11*Zuu22/cb**2*sb
      AAABR(864) = EE**2/SW**2*Sqrt2*Zb21*Zh32*Zt21/cb+EE
     &      **2/SW**2*Sqrt2*Zb21*Zh31*Zt21/cb**2*sb-EE**2
     &      /SW**2*Sqrt2*Zb21*Zh31*Zt21/cb**2*sb**3-EE**2
     &      /SW**2*Sqrt2*Zb21*Zh32*Zt21/cb*sb**2-2D0*EE*Lambda
     &      /MW*MbMM/SW*Zb22*Zh33*Zt21/cb**2-2D0*EE*Lambda
     &      /MW*MtMM/SW*Zb21*Zh33*Zt22/cb**2+EE**2/MW**2
     &      *MtMM**2/SW**2*Sqrt2*Zb21*Zh32*Zt21/cb+2D0*EE
     &      *Lambda/MW*MbMM/SW*Zb22*Zh33*Zt21/cb**2*sb**2
     &      +2D0*EE*Lambda/MW*MtMM/SW*Zb21*Zh33*Zt22/cb**2
     &      *sb**2-EE**2/MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb22
     &      *Zh32*Zt22/cb-EE**2/MW**2*MbMM**2/SW**2*Sqrt2
     &      *Zb21*Zh31*Zt21/cb**2*sb-EE**2/MW**2*MtMM**2
     &      /SW**2*Sqrt2*Zb21*Zh32*Zt21/cb/sb**2-EE**2/MW
     &      **2*MbMM*MtMM/SW**2*Sqrt2*Zb22*Zh31*Zt22/cb**2
     &      /sb+EE**2/MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb22*Zh31
     &      *Zt22/cb**2*sb
      AAABR(865) = EE**2/SW**2*Sqrt2*Zcc21*Zh22*Zss21/cb
     &      +EE**2/SW**2*Sqrt2*Zcc21*Zh21*Zss21/cb**2*sb
     &      -EE**2/SW**2*Sqrt2*Zcc21*Zh21*Zss21/cb**2*sb
     &      **3-EE**2/SW**2*Sqrt2*Zcc21*Zh22*Zss21/cb*sb
     &      **2-2D0*EE*Lambda/MW*Mc/SW*Zcc22*Zh23*Zss21/cb
     &      **2-2D0*EE*Lambda/MW*Mq/SW*Zcc21*Zh23*Zss22/cb
     &      **2+EE**2/MW**2*Mc**2/SW**2*Sqrt2*Zcc21*Zh22
     &      *Zss21/cb+2D0*EE*Lambda/MW*Mc/SW*Zcc22*Zh23*Zss21
     &      /cb**2*sb**2+2D0*EE*Lambda/MW*Mq/SW*Zcc21*Zh23
     &      *Zss22/cb**2*sb**2-EE**2/MW**2*Mc*Mq/SW**2*Sqrt2
     &      *Zcc22*Zh22*Zss22/cb-EE**2/MW**2*Mc**2/SW**2
     &      *Sqrt2*Zcc21*Zh22*Zss21/cb/sb**2-EE**2/MW**2
     &      *Mq**2/SW**2*Sqrt2*Zcc21*Zh21*Zss21/cb**2*sb
     &      -EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc22*Zh21*Zss22
     &      /cb**2/sb+EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc22
     &      *Zh21*Zss22/cb**2*sb
      AAABR(866) = EE**2/SW**2*Sqrt2*Zb21*Zh12*Zt21/cb+EE
     &      **2/SW**2*Sqrt2*Zb21*Zh11*Zt21/cb**2*sb-EE**2
     &      /SW**2*Sqrt2*Zb21*Zh11*Zt21/cb**2*sb**3-EE**2
     &      /SW**2*Sqrt2*Zb21*Zh12*Zt21/cb*sb**2-2D0*EE*Lambda
     &      /MW*MbMM/SW*Zb22*Zh13*Zt21/cb**2-2D0*EE*Lambda
     &      /MW*MtMM/SW*Zb21*Zh13*Zt22/cb**2+EE**2/MW**2
     &      *MtMM**2/SW**2*Sqrt2*Zb21*Zh12*Zt21/cb+2D0*EE
     &      *Lambda/MW*MbMM/SW*Zb22*Zh13*Zt21/cb**2*sb**2
     &      +2D0*EE*Lambda/MW*MtMM/SW*Zb21*Zh13*Zt22/cb**2
     &      *sb**2-EE**2/MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb22
     &      *Zh12*Zt22/cb-EE**2/MW**2*MbMM**2/SW**2*Sqrt2
     &      *Zb21*Zh11*Zt21/cb**2*sb-EE**2/MW**2*MtMM**2
     &      /SW**2*Sqrt2*Zb21*Zh12*Zt21/cb/sb**2-EE**2/MW
     &      **2*MbMM*MtMM/SW**2*Sqrt2*Zb22*Zh11*Zt22/cb**2
     &      /sb+EE**2/MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb22*Zh11
     &      *Zt22/cb**2*sb
      AAABR(867) = EE**2/SW**2*Sqrt2*Zb21*Zh22*Zt21/cb+EE
     &      **2/SW**2*Sqrt2*Zb21*Zh21*Zt21/cb**2*sb-EE**2
     &      /SW**2*Sqrt2*Zb21*Zh21*Zt21/cb**2*sb**3-EE**2
     &      /SW**2*Sqrt2*Zb21*Zh22*Zt21/cb*sb**2-2D0*EE*Lambda
     &      /MW*MbMM/SW*Zb22*Zh23*Zt21/cb**2-2D0*EE*Lambda
     &      /MW*MtMM/SW*Zb21*Zh23*Zt22/cb**2+EE**2/MW**2
     &      *MtMM**2/SW**2*Sqrt2*Zb21*Zh22*Zt21/cb+2D0*EE
     &      *Lambda/MW*MbMM/SW*Zb22*Zh23*Zt21/cb**2*sb**2
     &      +2D0*EE*Lambda/MW*MtMM/SW*Zb21*Zh23*Zt22/cb**2
     &      *sb**2-EE**2/MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb22
     &      *Zh22*Zt22/cb-EE**2/MW**2*MbMM**2/SW**2*Sqrt2
     &      *Zb21*Zh21*Zt21/cb**2*sb-EE**2/MW**2*MtMM**2
     &      /SW**2*Sqrt2*Zb21*Zh22*Zt21/cb/sb**2-EE**2/MW
     &      **2*MbMM*MtMM/SW**2*Sqrt2*Zb22*Zh21*Zt22/cb**2
     &      /sb+EE**2/MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb22*Zh21
     &      *Zt22/cb**2*sb
      AAABR(868) = EE**2/SW**2*Sqrt2*Zdd21*Zh32*Zuu21/cb
     &      +EE**2/SW**2*Sqrt2*Zdd21*Zh31*Zuu21/cb**2*sb
     &      -EE**2/SW**2*Sqrt2*Zdd21*Zh31*Zuu21/cb**2*sb
     &      **3-EE**2/SW**2*Sqrt2*Zdd21*Zh32*Zuu21/cb*sb
     &      **2-2D0*EE*Lambda/MW*Mq/SW*Zdd21*Zh33*Zuu22/cb
     &      **2-2D0*EE*Lambda/MW*Mq/SW*Zdd22*Zh33*Zuu21/cb
     &      **2+EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd21*Zh32
     &      *Zuu21/cb-EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd22
     &      *Zh32*Zuu22/cb+2D0*EE*Lambda/MW*Mq/SW*Zdd21*Zh33
     &      *Zuu22/cb**2*sb**2+2D0*EE*Lambda/MW*Mq/SW*Zdd22
     &      *Zh33*Zuu21/cb**2*sb**2-EE**2/MW**2*Mq**2/SW
     &      **2*Sqrt2*Zdd21*Zh31*Zuu21/cb**2*sb-EE**2/MW
     &      **2*Mq**2/SW**2*Sqrt2*Zdd21*Zh32*Zuu21/cb/sb
     &      **2-EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd22*Zh31
     &      *Zuu22/cb**2/sb+EE**2/MW**2*Mq**2/SW**2*Sqrt2
     &      *Zdd22*Zh31*Zuu22/cb**2*sb
      AAABR(869) = EE**2/SW**2*Sqrt2*Zdd21*Zh22*Zuu21/cb
     &      +EE**2/SW**2*Sqrt2*Zdd21*Zh21*Zuu21/cb**2*sb
     &      -EE**2/SW**2*Sqrt2*Zdd21*Zh21*Zuu21/cb**2*sb
     &      **3-EE**2/SW**2*Sqrt2*Zdd21*Zh22*Zuu21/cb*sb
     &      **2-2D0*EE*Lambda/MW*Mq/SW*Zdd21*Zh23*Zuu22/cb
     &      **2-2D0*EE*Lambda/MW*Mq/SW*Zdd22*Zh23*Zuu21/cb
     &      **2+EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd21*Zh22
     &      *Zuu21/cb-EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd22
     &      *Zh22*Zuu22/cb+2D0*EE*Lambda/MW*Mq/SW*Zdd21*Zh23
     &      *Zuu22/cb**2*sb**2+2D0*EE*Lambda/MW*Mq/SW*Zdd22
     &      *Zh23*Zuu21/cb**2*sb**2-EE**2/MW**2*Mq**2/SW
     &      **2*Sqrt2*Zdd21*Zh21*Zuu21/cb**2*sb-EE**2/MW
     &      **2*Mq**2/SW**2*Sqrt2*Zdd21*Zh22*Zuu21/cb/sb
     &      **2-EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd22*Zh21
     &      *Zuu22/cb**2/sb+EE**2/MW**2*Mq**2/SW**2*Sqrt2
     &      *Zdd22*Zh21*Zuu22/cb**2*sb
      AAABR(870) = EE**2/SW**2*Sqrt2*Zcc21*Zh32*Zss21/cb
     &      +EE**2/SW**2*Sqrt2*Zcc21*Zh31*Zss21/cb**2*sb
     &      -EE**2/SW**2*Sqrt2*Zcc21*Zh31*Zss21/cb**2*sb
     &      **3-EE**2/SW**2*Sqrt2*Zcc21*Zh32*Zss21/cb*sb
     &      **2-2D0*EE*Lambda/MW*Mc/SW*Zcc22*Zh33*Zss21/cb
     &      **2-2D0*EE*Lambda/MW*Mq/SW*Zcc21*Zh33*Zss22/cb
     &      **2+EE**2/MW**2*Mc**2/SW**2*Sqrt2*Zcc21*Zh32
     &      *Zss21/cb+2D0*EE*Lambda/MW*Mc/SW*Zcc22*Zh33*Zss21
     &      /cb**2*sb**2+2D0*EE*Lambda/MW*Mq/SW*Zcc21*Zh33
     &      *Zss22/cb**2*sb**2-EE**2/MW**2*Mc*Mq/SW**2*Sqrt2
     &      *Zcc22*Zh32*Zss22/cb-EE**2/MW**2*Mc**2/SW**2
     &      *Sqrt2*Zcc21*Zh32*Zss21/cb/sb**2-EE**2/MW**2
     &      *Mq**2/SW**2*Sqrt2*Zcc21*Zh31*Zss21/cb**2*sb
     &      -EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc22*Zh31*Zss22
     &      /cb**2/sb+EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc22
     &      *Zh31*Zss22/cb**2*sb
      AAABR(871) = EE**2/SW**2*Sqrt2*Zcc21*Zh12*Zss21/cb
     &      +EE**2/SW**2*Sqrt2*Zcc21*Zh11*Zss21/cb**2*sb
     &      -EE**2/SW**2*Sqrt2*Zcc21*Zh11*Zss21/cb**2*sb
     &      **3-EE**2/SW**2*Sqrt2*Zcc21*Zh12*Zss21/cb*sb
     &      **2-2D0*EE*Lambda/MW*Mc/SW*Zcc22*Zh13*Zss21/cb
     &      **2-2D0*EE*Lambda/MW*Mq/SW*Zcc21*Zh13*Zss22/cb
     &      **2+EE**2/MW**2*Mc**2/SW**2*Sqrt2*Zcc21*Zh12
     &      *Zss21/cb+2D0*EE*Lambda/MW*Mc/SW*Zcc22*Zh13*Zss21
     &      /cb**2*sb**2+2D0*EE*Lambda/MW*Mq/SW*Zcc21*Zh13
     &      *Zss22/cb**2*sb**2-EE**2/MW**2*Mc*Mq/SW**2*Sqrt2
     &      *Zcc22*Zh12*Zss22/cb-EE**2/MW**2*Mc**2/SW**2
     &      *Sqrt2*Zcc21*Zh12*Zss21/cb/sb**2-EE**2/MW**2
     &      *Mq**2/SW**2*Sqrt2*Zcc21*Zh11*Zss21/cb**2*sb
     &      -EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc22*Zh11*Zss22
     &      /cb**2/sb+EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc22
     &      *Zh11*Zss22/cb**2*sb
      AAABR(872) = la6-5D0*la6*sb**2+4D0*la6*sb**4+3D0*la7
     &      *sb**2-4D0*la7*sb**4+cb*la1*sb-cb*la1*sb**3-cb
     &      *la2*sb**3-cb*la3*sb+2D0*cb*la3*sb**3-cb*la4
     &      *sb+2D0*cb*la4*sb**3-cb*la5*sb+2D0*cb*la5*sb**3
      AAABR(873) = Pa21*(la4+la5+2D0*la1*sb**2-2D0*la1*sb
     &      **4+2D0*la2*sb**2-2D0*la2*sb**4-4D0*la3*sb**2
     &      +4D0*la3*sb**4-4D0*la4*sb**2+4D0*la4*sb**4-4D0
     &      *la5*sb**2+4D0*la5*sb**4+4D0*cb*la6*sb-8D0*cb
     &      *la6*sb**3-4D0*cb*la7*sb+8D0*cb*la7*sb**3)
      AAABR(874) = Pa11*(la4+la5+2D0*la1*sb**2-2D0*la1*sb
     &      **4+2D0*la2*sb**2-2D0*la2*sb**4-4D0*la3*sb**2
     &      +4D0*la3*sb**4-4D0*la4*sb**2+4D0*la4*sb**4-4D0
     &      *la5*sb**2+4D0*la5*sb**4+4D0*cb*la6*sb-8D0*cb
     &      *la6*sb**3-4D0*cb*la7*sb+8D0*cb*la7*sb**3)
      AAABR(875) = Zh11*la4*sb-Zh11*la5*sb-Zh12*cb*la4+Zh12
     &      *cb*la5
      AAABR(876) = Zh31*la4*sb-Zh31*la5*sb-Zh32*cb*la4+Zh32
     &      *cb*la5
      AAABR(877) = Zh21*la4*sb-Zh21*la5*sb-Zh22*cb*la4+Zh22
     &      *cb*la5
      AAABR(878) = Pa11*Pa21*la7+Pa12*Pa22*la4s-Pa12*Pa22
     &      *la5s-Pa12*Pa22*la6s+3D0*Pa11*Pa21*la6*sb**2
     &      -4D0*Pa11*Pa21*la6*sb**4-5D0*Pa11*Pa21*la7*sb
     &      **2+4D0*Pa11*Pa21*la7*sb**4-2D0*Pa12*Pa22*la4s
     &      *sb**2+2D0*Pa12*Pa22*la5s*sb**2+2D0*Pa12*Pa22
     &      *la6s*sb**2+Pa11*Pa21*cb*la1*sb**3-Pa11*Pa21
     &      *cb*la2*sb+Pa11*Pa21*cb*la2*sb**3+Pa11*Pa21*cb
     &      *la3*sb-2D0*Pa11*Pa21*cb*la3*sb**3+Pa11*Pa21
     &      *cb*la4*sb-2D0*Pa11*Pa21*cb*la4*sb**3+Pa11*Pa21
     &      *cb*la5*sb-2D0*Pa11*Pa21*cb*la5*sb**3+Pa12*Pa22
     &      *cb*la1s*sb-Pa12*Pa22*cb*la2s*sb
      AAABR(879) = Pa21**2*la7+Pa22**2*la4s-Pa22**2*la5s
     &      -Pa22**2*la6s+3D0*Pa21**2*la6*sb**2-4D0*Pa21
     &      **2*la6*sb**4-5D0*Pa21**2*la7*sb**2+4D0*Pa21
     &      **2*la7*sb**4-2D0*Pa22**2*la4s*sb**2+2D0*Pa22
     &      **2*la5s*sb**2+2D0*Pa22**2*la6s*sb**2+Pa21**2
     &      *cb*la1*sb**3-Pa21**2*cb*la2*sb+Pa21**2*cb*la2
     &      *sb**3+Pa21**2*cb*la3*sb-2D0*Pa21**2*cb*la3*sb
     &      **3+Pa21**2*cb*la4*sb-2D0*Pa21**2*cb*la4*sb**3
     &      +Pa21**2*cb*la5*sb-2D0*Pa21**2*cb*la5*sb**3+Pa22
     &      **2*cb*la1s*sb-Pa22**2*cb*la2s*sb
      AAABR(880) = Pa11**2*la7+Pa12**2*la4s-Pa12**2*la5s
     &      -Pa12**2*la6s+3D0*Pa11**2*la6*sb**2-4D0*Pa11
     &      **2*la6*sb**4-5D0*Pa11**2*la7*sb**2+4D0*Pa11
     &      **2*la7*sb**4-2D0*Pa12**2*la4s*sb**2+2D0*Pa12
     &      **2*la5s*sb**2+2D0*Pa12**2*la6s*sb**2+Pa11**2
     &      *cb*la1*sb**3-Pa11**2*cb*la2*sb+Pa11**2*cb*la2
     &      *sb**3+Pa11**2*cb*la3*sb-2D0*Pa11**2*cb*la3*sb
     &      **3+Pa11**2*cb*la4*sb-2D0*Pa11**2*cb*la4*sb**3
     &      +Pa11**2*cb*la5*sb-2D0*Pa11**2*cb*la5*sb**3+Pa12
     &      **2*cb*la1s*sb-Pa12**2*cb*la2s*sb
      AAABR(881) = 2D0*Pa22*Zh13*la5s-2D0*Pa22*Zh13*la6s
     &      +Pa21*Zh11*cb*la4-Pa21*Zh11*cb*la5+Pa21*Zh12
     &      *la4*sb-Pa21*Zh12*la5*sb
      AAABR(882) = 2D0*Pa22*Zh33*la5s-2D0*Pa22*Zh33*la6s
     &      +Pa21*Zh31*cb*la4-Pa21*Zh31*cb*la5+Pa21*Zh32
     &      *la4*sb-Pa21*Zh32*la5*sb
      AAABR(883) = 2D0*Pa22*Zh23*la5s-2D0*Pa22*Zh23*la6s
     &      +Pa21*Zh21*cb*la4-Pa21*Zh21*cb*la5+Pa21*Zh22
     &      *la4*sb-Pa21*Zh22*la5*sb
      AAABR(884) = 2D0*Pa12*Zh33*la5s-2D0*Pa12*Zh33*la6s
     &      +Pa11*Zh31*cb*la4-Pa11*Zh31*cb*la5+Pa11*Zh32
     &      *la4*sb-Pa11*Zh32*la5*sb
      AAABR(885) = 2D0*Pa12*Zh13*la5s-2D0*Pa12*Zh13*la6s
     &      +Pa11*Zh11*cb*la4-Pa11*Zh11*cb*la5+Pa11*Zh12
     &      *la4*sb-Pa11*Zh12*la5*sb
      AAABR(886) = 2D0*Pa12*Zh23*la5s-2D0*Pa12*Zh23*la6s
     &      +Pa11*Zh21*cb*la4-Pa11*Zh21*cb*la5+Pa11*Zh22
     &      *la4*sb-Pa11*Zh22*la5*sb
      AAABR(887) = 2D0*Zh11*Zh21*la6-Zh11*Zh22*la4-Zh11*Zh22
     &      *la5-Zh12*Zh21*la4-Zh12*Zh21*la5+2D0*Zh12*Zh22
     &      *la7+2D0*Zh13*Zh23*la4s+2D0*Zh13*Zh23*la5s+2D0
     &      *Zh13*Zh23*la6s-4D0*Zh11*Zh21*la6*sb**2+2D0*Zh11
     &      *Zh22*la4*sb**2+2D0*Zh11*Zh22*la5*sb**2+2D0*Zh12
     &      *Zh21*la4*sb**2+2D0*Zh12*Zh21*la5*sb**2-4D0*Zh12
     &      *Zh22*la7*sb**2-4D0*Zh13*Zh23*la4s*sb**2-4D0
     &      *Zh13*Zh23*la5s*sb**2-4D0*Zh13*Zh23*la6s*sb**2
     &      +2D0*Zh11*Zh21*cb*la1*sb-2D0*Zh11*Zh21*cb*la3
     &      *sb-2D0*Zh11*Zh22*cb*la6*sb+2D0*Zh11*Zh22*cb
     &      *la7*sb-2D0*Zh12*Zh21*cb*la6*sb+2D0*Zh12*Zh21
     &      *cb*la7*sb-2D0*Zh12*Zh22*cb*la2*sb+2D0*Zh12*Zh22
     &      *cb*la3*sb+2D0*Zh13*Zh23*cb*la1s*sb-2D0*Zh13
     &      *Zh23*cb*la2s*sb
      AAABR(888) = Zh31**2*la6+Zh32**2*la7+Zh33**2*la4s+Zh33
     &      **2*la5s+Zh33**2*la6s-Zh31*Zh32*la4-Zh31*Zh32
     &      *la5-2D0*Zh31**2*la6*sb**2-2D0*Zh32**2*la7*sb
     &      **2-2D0*Zh33**2*la4s*sb**2-2D0*Zh33**2*la5s*sb
     &      **2-2D0*Zh33**2*la6s*sb**2+2D0*Zh31*Zh32*la4
     &      *sb**2+2D0*Zh31*Zh32*la5*sb**2+Zh31**2*cb*la1
     &      *sb-Zh31**2*cb*la3*sb-Zh32**2*cb*la2*sb+Zh32
     &      **2*cb*la3*sb+Zh33**2*cb*la1s*sb-Zh33**2*cb*la2s
     &      *sb-2D0*Zh31*Zh32*cb*la6*sb+2D0*Zh31*Zh32*cb
     &      *la7*sb
      AAABR(889) = 2D0*Zh11*Zh31*la6-Zh11*Zh32*la4-Zh11*Zh32
     &      *la5-Zh12*Zh31*la4-Zh12*Zh31*la5+2D0*Zh12*Zh32
     &      *la7+2D0*Zh13*Zh33*la4s+2D0*Zh13*Zh33*la5s+2D0
     &      *Zh13*Zh33*la6s-4D0*Zh11*Zh31*la6*sb**2+2D0*Zh11
     &      *Zh32*la4*sb**2+2D0*Zh11*Zh32*la5*sb**2+2D0*Zh12
     &      *Zh31*la4*sb**2+2D0*Zh12*Zh31*la5*sb**2-4D0*Zh12
     &      *Zh32*la7*sb**2-4D0*Zh13*Zh33*la4s*sb**2-4D0
     &      *Zh13*Zh33*la5s*sb**2-4D0*Zh13*Zh33*la6s*sb**2
     &      +2D0*Zh11*Zh31*cb*la1*sb-2D0*Zh11*Zh31*cb*la3
     &      *sb-2D0*Zh11*Zh32*cb*la6*sb+2D0*Zh11*Zh32*cb
     &      *la7*sb-2D0*Zh12*Zh31*cb*la6*sb+2D0*Zh12*Zh31
     &      *cb*la7*sb-2D0*Zh12*Zh32*cb*la2*sb+2D0*Zh12*Zh32
     &      *cb*la3*sb+2D0*Zh13*Zh33*cb*la1s*sb-2D0*Zh13
     &      *Zh33*cb*la2s*sb
      AAABR(890) = Zh21**2*la6+Zh22**2*la7+Zh23**2*la4s+Zh23
     &      **2*la5s+Zh23**2*la6s-Zh21*Zh22*la4-Zh21*Zh22
     &      *la5-2D0*Zh21**2*la6*sb**2-2D0*Zh22**2*la7*sb
     &      **2-2D0*Zh23**2*la4s*sb**2-2D0*Zh23**2*la5s*sb
     &      **2-2D0*Zh23**2*la6s*sb**2+2D0*Zh21*Zh22*la4
     &      *sb**2+2D0*Zh21*Zh22*la5*sb**2+Zh21**2*cb*la1
     &      *sb-Zh21**2*cb*la3*sb-Zh22**2*cb*la2*sb+Zh22
     &      **2*cb*la3*sb+Zh23**2*cb*la1s*sb-Zh23**2*cb*la2s
     &      *sb-2D0*Zh21*Zh22*cb*la6*sb+2D0*Zh21*Zh22*cb
     &      *la7*sb
      AAABR(891) = Zh11**2*la6+Zh12**2*la7+Zh13**2*la4s+Zh13
     &      **2*la5s+Zh13**2*la6s-Zh11*Zh12*la4-Zh11*Zh12
     &      *la5-2D0*Zh11**2*la6*sb**2-2D0*Zh12**2*la7*sb
     &      **2-2D0*Zh13**2*la4s*sb**2-2D0*Zh13**2*la5s*sb
     &      **2-2D0*Zh13**2*la6s*sb**2+2D0*Zh11*Zh12*la4
     &      *sb**2+2D0*Zh11*Zh12*la5*sb**2+Zh11**2*cb*la1
     &      *sb-Zh11**2*cb*la3*sb-Zh12**2*cb*la2*sb+Zh12
     &      **2*cb*la3*sb+Zh13**2*cb*la1s*sb-Zh13**2*cb*la2s
     &      *sb-2D0*Zh11*Zh12*cb*la6*sb+2D0*Zh11*Zh12*cb
     &      *la7*sb
      AAABR(892) = 2D0*Zh21*Zh31*la6-Zh21*Zh32*la4-Zh21*Zh32
     &      *la5-Zh22*Zh31*la4-Zh22*Zh31*la5+2D0*Zh22*Zh32
     &      *la7+2D0*Zh23*Zh33*la4s+2D0*Zh23*Zh33*la5s+2D0
     &      *Zh23*Zh33*la6s-4D0*Zh21*Zh31*la6*sb**2+2D0*Zh21
     &      *Zh32*la4*sb**2+2D0*Zh21*Zh32*la5*sb**2+2D0*Zh22
     &      *Zh31*la4*sb**2+2D0*Zh22*Zh31*la5*sb**2-4D0*Zh22
     &      *Zh32*la7*sb**2-4D0*Zh23*Zh33*la4s*sb**2-4D0
     &      *Zh23*Zh33*la5s*sb**2-4D0*Zh23*Zh33*la6s*sb**2
     &      +2D0*Zh21*Zh31*cb*la1*sb-2D0*Zh21*Zh31*cb*la3
     &      *sb-2D0*Zh21*Zh32*cb*la6*sb+2D0*Zh21*Zh32*cb
     &      *la7*sb-2D0*Zh22*Zh31*cb*la6*sb+2D0*Zh22*Zh31
     &      *cb*la7*sb-2D0*Zh22*Zh32*cb*la2*sb+2D0*Zh22*Zh32
     &      *cb*la3*sb+2D0*Zh23*Zh33*cb*la1s*sb-2D0*Zh23
     &      *Zh33*cb*la2s*sb
      AAABR(893) = EE**2*(Zl11**2/CW**2-2D0*Zl12**2/CW**2
     &      +Zl11**2/SW**2)
      AAABR(894) = EE**2*(1D0/CW**2+1D0/SW**2)
      AAABR(895) = EE**2/cb**2*(Zl11**4/CW**2+4D0*Zl12**4
     &      /CW**2+Zl11**4/SW**2-4D0*Zl11**2/CW**2*Zl12**2
     &      -Zl11**4/CW**2*sb**2-4D0*Zl12**4/CW**2*sb**2
     &      -Zl11**4/SW**2*sb**2+4D0*Zl11**2/CW**2*Zl12**2
     &      *sb**2+4D0*Ml**2/MW**2/SW**2*Zl11**2*Zl12**2)
      AAABR(896) = EE**2/cb**2*(Zl11**3/CW**2*Zl21+4D0*Zl12
     &      **3/CW**2*Zl22+Zl11**3/SW**2*Zl21-2D0*Zl11/CW
     &      **2*Zl12**2*Zl21-2D0*Zl11**2/CW**2*Zl12*Zl22
     &      -Zl11**3/CW**2*Zl21*sb**2-4D0*Zl12**3/CW**2*Zl22
     &      *sb**2-Zl11**3/SW**2*Zl21*sb**2+2D0*Zl11/CW**2
     &      *Zl12**2*Zl21*sb**2+2D0*Zl11**2/CW**2*Zl12*Zl22
     &      *sb**2+2D0*Ml**2/MW**2/SW**2*Zl11*Zl12**2*Zl21
     &      +2D0*Ml**2/MW**2/SW**2*Zl11**2*Zl12*Zl22)
      AAABR(897) = EE**2*(Zl11/CW**2*Zl21-2D0*Zl12/CW**2
     &      *Zl22+Zl11/SW**2*Zl21)
      AAABR(898) = EE**2/cb**2*(Zl11**2/CW**2*Zl21**2+4D0
     &      *Zl12**2/CW**2*Zl22**2+Zl11**2/SW**2*Zl21**2
     &      -Zl11**2/CW**2*Zl21**2*sb**2-4D0*Zl12**2/CW**2
     &      *Zl22**2*sb**2-Zl11**2/SW**2*Zl21**2*sb**2-4D0
     &      *Zl11/CW**2*Zl12*Zl21*Zl22+4D0*Zl11/CW**2*Zl12
     &      *Zl21*Zl22*sb**2+4D0*Ml**2/MW**2/SW**2*Zl11*Zl12
     &      *Zl21*Zl22)
      AAABR(899) = EE**2/CW**2*(Zl11**2-2D0*Zl12**2)
      AAABR(900) = EE**2/CW**2
      AAABR(901) = EE**2/cb**2*(Zl11**2/CW**2*Zl21**2-Zl11
     &      **2/CW**2*Zl22**2-Zl12**2/CW**2*Zl21**2+4D0*Zl12
     &      **2/CW**2*Zl22**2+Zl11**2/SW**2*Zl21**2-Zl11
     &      **2/CW**2*Zl21**2*sb**2+Zl11**2/CW**2*Zl22**2
     &      *sb**2+Zl12**2/CW**2*Zl21**2*sb**2-4D0*Zl12**2
     &      /CW**2*Zl22**2*sb**2-Zl11**2/SW**2*Zl21**2*sb
     &      **2-2D0*Zl11/CW**2*Zl12*Zl21*Zl22+Ml**2/MW**2
     &      /SW**2*Zl11**2*Zl22**2+Ml**2/MW**2/SW**2*Zl12
     &      **2*Zl21**2+2D0*Zl11/CW**2*Zl12*Zl21*Zl22*sb
     &      **2+2D0*Ml**2/MW**2/SW**2*Zl11*Zl12*Zl21*Zl22)
      AAABR(902) = EE**2*(Zl21**2/CW**2-2D0*Zl22**2/CW**2
     &      +Zl21**2/SW**2)
      AAABR(903) = EE**2/CW**2*(Zl11*Zl21-2D0*Zl12*Zl22)
      AAABR(904) = EE**2/cb**2*(Zl11/CW**2*Zl21**3+4D0*Zl12
     &      /CW**2*Zl22**3+Zl11/SW**2*Zl21**3-2D0*Zl11/CW
     &      **2*Zl21*Zl22**2-Zl11/CW**2*Zl21**3*sb**2-2D0
     &      *Zl12/CW**2*Zl21**2*Zl22-4D0*Zl12/CW**2*Zl22
     &      **3*sb**2-Zl11/SW**2*Zl21**3*sb**2+2D0*Zl11/CW
     &      **2*Zl21*Zl22**2*sb**2+2D0*Zl12/CW**2*Zl21**2
     &      *Zl22*sb**2+2D0*Ml**2/MW**2/SW**2*Zl11*Zl21*Zl22
     &      **2+2D0*Ml**2/MW**2/SW**2*Zl12*Zl21**2*Zl22)
      AAABR(905) = EE**2/SW**2*Zl11
      AAABR(906) = EE**2*(1D0/CW**2-1D0/SW**2)
      AAABR(907) = EE**2/SW**2
      AAABR(908) = EE**2*(Zl11**2/CW**2-2D0*Zl12**2/CW**2
     &      -Zl11**2/SW**2)
      AAABR(909) = EE**2/cb**2*(Zl11**2/CW**2-2D0*Zl12**2
     &      /CW**2+Zl11**2/SW**2-Zl11**2/CW**2*sb**2+2D0
     &      *Zl12**2/CW**2*sb**2-Zl11**2/SW**2*sb**2+2D0
     &      *Ml**2/MW**2/SW**2*Zl12**2)
      AAABR(910) = EE**2*(Zl11/CW**2*Zl21-2D0*Zl12/CW**2
     &      *Zl22-Zl11/SW**2*Zl21)
      AAABR(911) = EE**2/SW**2*Zl21
      AAABR(912) = EE**2/cb**2*(Zl11/CW**2*Zl21-2D0*Zl12
     &      /CW**2*Zl22+Zl11/SW**2*Zl21-Zl11/CW**2*Zl21*sb
     &      **2+2D0*Zl12/CW**2*Zl22*sb**2-Zl11/SW**2*Zl21
     &      *sb**2+2D0*Ml**2/MW**2/SW**2*Zl12*Zl22)
      AAABR(913) = EE**2*(Zdd11**2/CW**2+2D0*Zdd12**2/CW
     &      **2-3D0*Zdd11**2/SW**2)
      AAABR(914) = EE**2*(Zss11**2/CW**2+2D0*Zss12**2/CW
     &      **2-3D0*Zss11**2/SW**2)
      AAABR(915) = EE**2*(Zb11**2/CW**2+2D0*Zb12**2/CW**2
     &      -3D0*Zb11**2/SW**2)
      AAABR(916) = EE**2/cb**2*(Zl11**2/CW**2*Zss11**2+2D0
     &      *Zl11**2/CW**2*Zss12**2-2D0*Zl12**2/CW**2*Zss11
     &      **2-4D0*Zl12**2/CW**2*Zss12**2-3D0*Zl11**2/SW
     &      **2*Zss11**2-Zl11**2/CW**2*Zss11**2*sb**2-2D0
     &      *Zl11**2/CW**2*Zss12**2*sb**2+2D0*Zl12**2/CW
     &      **2*Zss11**2*sb**2+4D0*Zl12**2/CW**2*Zss12**2
     &      *sb**2+3D0*Zl11**2/SW**2*Zss11**2*sb**2-12D0
     &      *Ml/MW**2*Mq/SW**2*Zl11*Zl12*Zss11*Zss12)
      AAABR(917) = EE**2/cb**2*(Zdd11**2/CW**2*Zl11**2-2D0
     &      *Zdd11**2/CW**2*Zl12**2+2D0*Zdd12**2/CW**2*Zl11
     &      **2-4D0*Zdd12**2/CW**2*Zl12**2-3D0*Zdd11**2/SW
     &      **2*Zl11**2-Zdd11**2/CW**2*Zl11**2*sb**2+2D0
     &      *Zdd11**2/CW**2*Zl12**2*sb**2-2D0*Zdd12**2/CW
     &      **2*Zl11**2*sb**2+4D0*Zdd12**2/CW**2*Zl12**2
     &      *sb**2+3D0*Zdd11**2/SW**2*Zl11**2*sb**2-12D0
     &      *Ml/MW**2*Mq/SW**2*Zdd11*Zdd12*Zl11*Zl12)
      AAABR(918) = EE**2/cb**2*(Zb11**2/CW**2*Zl11**2-2D0
     &      *Zb11**2/CW**2*Zl12**2+2D0*Zb12**2/CW**2*Zl11
     &      **2-4D0*Zb12**2/CW**2*Zl12**2-3D0*Zb11**2/SW
     &      **2*Zl11**2-Zb11**2/CW**2*Zl11**2*sb**2+2D0*Zb11
     &      **2/CW**2*Zl12**2*sb**2-2D0*Zb12**2/CW**2*Zl11
     &      **2*sb**2+4D0*Zb12**2/CW**2*Zl12**2*sb**2+3D0
     &      *Zb11**2/SW**2*Zl11**2*sb**2-12D0*MbMM/MW**2
     &      *Ml/SW**2*Zb11*Zb12*Zl11*Zl12)
      AAABR(919) = EE**2/cb**2*(Zl11**2/CW**2*Zss11*Zss21
     &      +2D0*Zl11**2/CW**2*Zss12*Zss22-2D0*Zl12**2/CW
     &      **2*Zss11*Zss21-4D0*Zl12**2/CW**2*Zss12*Zss22
     &      -3D0*Zl11**2/SW**2*Zss11*Zss21-Zl11**2/CW**2
     &      *Zss11*Zss21*sb**2-2D0*Zl11**2/CW**2*Zss12*Zss22
     &      *sb**2+2D0*Zl12**2/CW**2*Zss11*Zss21*sb**2+4D0
     &      *Zl12**2/CW**2*Zss12*Zss22*sb**2+3D0*Zl11**2
     &      /SW**2*Zss11*Zss21*sb**2-6D0*Ml/MW**2*Mq/SW**2
     &      *Zl11*Zl12*Zss11*Zss22-6D0*Ml/MW**2*Mq/SW**2
     &      *Zl11*Zl12*Zss12*Zss21)
      AAABR(920) = EE**2*(Zdd11/CW**2*Zdd21+2D0*Zdd12/CW
     &      **2*Zdd22-3D0*Zdd11/SW**2*Zdd21)
      AAABR(921) = EE**2/cb**2*(Zb11/CW**2*Zb21*Zl11**2-2D0
     &      *Zb11/CW**2*Zb21*Zl12**2+2D0*Zb12/CW**2*Zb22
     &      *Zl11**2-4D0*Zb12/CW**2*Zb22*Zl12**2-3D0*Zb11
     &      /SW**2*Zb21*Zl11**2-Zb11/CW**2*Zb21*Zl11**2*sb
     &      **2+2D0*Zb11/CW**2*Zb21*Zl12**2*sb**2-2D0*Zb12
     &      /CW**2*Zb22*Zl11**2*sb**2+4D0*Zb12/CW**2*Zb22
     &      *Zl12**2*sb**2+3D0*Zb11/SW**2*Zb21*Zl11**2*sb
     &      **2-6D0*MbMM/MW**2*Ml/SW**2*Zb11*Zb22*Zl11*Zl12
     &      -6D0*MbMM/MW**2*Ml/SW**2*Zb12*Zb21*Zl11*Zl12)
      AAABR(922) = EE**2*(Zb11/CW**2*Zb21+2D0*Zb12/CW**2
     &      *Zb22-3D0*Zb11/SW**2*Zb21)
      AAABR(923) = EE**2/cb**2*(Zdd11/CW**2*Zdd21*Zl11**2
     &      -2D0*Zdd11/CW**2*Zdd21*Zl12**2+2D0*Zdd12/CW**2
     &      *Zdd22*Zl11**2-4D0*Zdd12/CW**2*Zdd22*Zl12**2
     &      -3D0*Zdd11/SW**2*Zdd21*Zl11**2-Zdd11/CW**2*Zdd21
     &      *Zl11**2*sb**2+2D0*Zdd11/CW**2*Zdd21*Zl12**2
     &      *sb**2-2D0*Zdd12/CW**2*Zdd22*Zl11**2*sb**2+4D0
     &      *Zdd12/CW**2*Zdd22*Zl12**2*sb**2+3D0*Zdd11/SW
     &      **2*Zdd21*Zl11**2*sb**2-6D0*Ml/MW**2*Mq/SW**2
     &      *Zdd11*Zdd22*Zl11*Zl12-6D0*Ml/MW**2*Mq/SW**2
     &      *Zdd12*Zdd21*Zl11*Zl12)
      AAABR(924) = EE**2*(Zss11/CW**2*Zss21+2D0*Zss12/CW
     &      **2*Zss22-3D0*Zss11/SW**2*Zss21)
      AAABR(925) = EE**2/cb**2*(Zl11/CW**2*Zl21*Zss11**2
     &      +2D0*Zl11/CW**2*Zl21*Zss12**2-2D0*Zl12/CW**2
     &      *Zl22*Zss11**2-4D0*Zl12/CW**2*Zl22*Zss12**2-3D0
     &      *Zl11/SW**2*Zl21*Zss11**2-Zl11/CW**2*Zl21*Zss11
     &      **2*sb**2-2D0*Zl11/CW**2*Zl21*Zss12**2*sb**2
     &      +2D0*Zl12/CW**2*Zl22*Zss11**2*sb**2+4D0*Zl12
     &      /CW**2*Zl22*Zss12**2*sb**2+3D0*Zl11/SW**2*Zl21
     &      *Zss11**2*sb**2-6D0*Ml/MW**2*Mq/SW**2*Zl11*Zl22
     &      *Zss11*Zss12-6D0*Ml/MW**2*Mq/SW**2*Zl12*Zl21
     &      *Zss11*Zss12)
      AAABR(926) = EE**2/cb**2*(Zb11**2/CW**2*Zl11*Zl21-2D0
     &      *Zb11**2/CW**2*Zl12*Zl22+2D0*Zb12**2/CW**2*Zl11
     &      *Zl21-4D0*Zb12**2/CW**2*Zl12*Zl22-3D0*Zb11**2
     &      /SW**2*Zl11*Zl21-Zb11**2/CW**2*Zl11*Zl21*sb**2
     &      +2D0*Zb11**2/CW**2*Zl12*Zl22*sb**2-2D0*Zb12**2
     &      /CW**2*Zl11*Zl21*sb**2+4D0*Zb12**2/CW**2*Zl12
     &      *Zl22*sb**2+3D0*Zb11**2/SW**2*Zl11*Zl21*sb**2
     &      -6D0*MbMM/MW**2*Ml/SW**2*Zb11*Zb12*Zl11*Zl22
     &      -6D0*MbMM/MW**2*Ml/SW**2*Zb11*Zb12*Zl12*Zl21)
      AAABR(927) = EE**2/cb**2*(Zdd11**2/CW**2*Zl11*Zl21
     &      -2D0*Zdd11**2/CW**2*Zl12*Zl22+2D0*Zdd12**2/CW
     &      **2*Zl11*Zl21-4D0*Zdd12**2/CW**2*Zl12*Zl22-3D0
     &      *Zdd11**2/SW**2*Zl11*Zl21-Zdd11**2/CW**2*Zl11
     &      *Zl21*sb**2+2D0*Zdd11**2/CW**2*Zl12*Zl22*sb**2
     &      -2D0*Zdd12**2/CW**2*Zl11*Zl21*sb**2+4D0*Zdd12
     &      **2/CW**2*Zl12*Zl22*sb**2+3D0*Zdd11**2/SW**2
     &      *Zl11*Zl21*sb**2-6D0*Ml/MW**2*Mq/SW**2*Zdd11
     &      *Zdd12*Zl11*Zl22-6D0*Ml/MW**2*Mq/SW**2*Zdd11
     &      *Zdd12*Zl12*Zl21)
      AAABR(928) = EE**2/cb**2*(Zl11/CW**2*Zl21*Zss11*Zss21
     &      +2D0*Zl11/CW**2*Zl21*Zss12*Zss22-2D0*Zl12/CW
     &      **2*Zl22*Zss11*Zss21-4D0*Zl12/CW**2*Zl22*Zss12
     &      *Zss22-3D0*Zl11/SW**2*Zl21*Zss11*Zss21-Zl11/CW
     &      **2*Zl21*Zss11*Zss21*sb**2-2D0*Zl11/CW**2*Zl21
     &      *Zss12*Zss22*sb**2+2D0*Zl12/CW**2*Zl22*Zss11
     &      *Zss21*sb**2+4D0*Zl12/CW**2*Zl22*Zss12*Zss22
     &      *sb**2+3D0*Zl11/SW**2*Zl21*Zss11*Zss21*sb**2
     &      -6D0*Ml/MW**2*Mq/SW**2*Zl11*Zl22*Zss12*Zss21
     &      -6D0*Ml/MW**2*Mq/SW**2*Zl12*Zl21*Zss11*Zss22)
      AAABR(929) = EE**2/cb**2*(Zb11/CW**2*Zb21*Zl11*Zl21
     &      -2D0*Zb11/CW**2*Zb21*Zl12*Zl22+2D0*Zb12/CW**2
     &      *Zb22*Zl11*Zl21-4D0*Zb12/CW**2*Zb22*Zl12*Zl22
     &      -3D0*Zb11/SW**2*Zb21*Zl11*Zl21-Zb11/CW**2*Zb21
     &      *Zl11*Zl21*sb**2+2D0*Zb11/CW**2*Zb21*Zl12*Zl22
     &      *sb**2-2D0*Zb12/CW**2*Zb22*Zl11*Zl21*sb**2+4D0
     &      *Zb12/CW**2*Zb22*Zl12*Zl22*sb**2+3D0*Zb11/SW
     &      **2*Zb21*Zl11*Zl21*sb**2-6D0*MbMM/MW**2*Ml/SW
     &      **2*Zb11*Zb22*Zl12*Zl21-6D0*MbMM/MW**2*Ml/SW
     &      **2*Zb12*Zb21*Zl11*Zl22)
      AAABR(930) = EE**2/cb**2*(Zdd11/CW**2*Zdd21*Zl11*Zl21
     &      -2D0*Zdd11/CW**2*Zdd21*Zl12*Zl22+2D0*Zdd12/CW
     &      **2*Zdd22*Zl11*Zl21-4D0*Zdd12/CW**2*Zdd22*Zl12
     &      *Zl22-3D0*Zdd11/SW**2*Zdd21*Zl11*Zl21-Zdd11/CW
     &      **2*Zdd21*Zl11*Zl21*sb**2+2D0*Zdd11/CW**2*Zdd21
     &      *Zl12*Zl22*sb**2-2D0*Zdd12/CW**2*Zdd22*Zl11*Zl21
     &      *sb**2+4D0*Zdd12/CW**2*Zdd22*Zl12*Zl22*sb**2
     &      +3D0*Zdd11/SW**2*Zdd21*Zl11*Zl21*sb**2-6D0*Ml
     &      /MW**2*Mq/SW**2*Zdd11*Zdd22*Zl12*Zl21-6D0*Ml
     &      /MW**2*Mq/SW**2*Zdd12*Zdd21*Zl11*Zl22)
      AAABR(931) = EE**2*(Zss21**2/CW**2+2D0*Zss22**2/CW
     &      **2-3D0*Zss21**2/SW**2)
      AAABR(932) = EE**2*(Zdd21**2/CW**2+2D0*Zdd22**2/CW
     &      **2-3D0*Zdd21**2/SW**2)
      AAABR(933) = EE**2/cb**2*(Zb21**2/CW**2*Zl11**2-2D0
     &      *Zb21**2/CW**2*Zl12**2+2D0*Zb22**2/CW**2*Zl11
     &      **2-4D0*Zb22**2/CW**2*Zl12**2-3D0*Zb21**2/SW
     &      **2*Zl11**2-Zb21**2/CW**2*Zl11**2*sb**2+2D0*Zb21
     &      **2/CW**2*Zl12**2*sb**2-2D0*Zb22**2/CW**2*Zl11
     &      **2*sb**2+4D0*Zb22**2/CW**2*Zl12**2*sb**2+3D0
     &      *Zb21**2/SW**2*Zl11**2*sb**2-12D0*MbMM/MW**2
     &      *Ml/SW**2*Zb21*Zb22*Zl11*Zl12)
      AAABR(934) = EE**2/cb**2*(Zl11**2/CW**2*Zss21**2+2D0
     &      *Zl11**2/CW**2*Zss22**2-2D0*Zl12**2/CW**2*Zss21
     &      **2-4D0*Zl12**2/CW**2*Zss22**2-3D0*Zl11**2/SW
     &      **2*Zss21**2-Zl11**2/CW**2*Zss21**2*sb**2-2D0
     &      *Zl11**2/CW**2*Zss22**2*sb**2+2D0*Zl12**2/CW
     &      **2*Zss21**2*sb**2+4D0*Zl12**2/CW**2*Zss22**2
     &      *sb**2+3D0*Zl11**2/SW**2*Zss21**2*sb**2-12D0
     &      *Ml/MW**2*Mq/SW**2*Zl11*Zl12*Zss21*Zss22)
      AAABR(935) = EE**2*(Zb21**2/CW**2+2D0*Zb22**2/CW**2
     &      -3D0*Zb21**2/SW**2)
      AAABR(936) = EE**2/cb**2*(Zdd21**2/CW**2*Zl11**2-2D0
     &      *Zdd21**2/CW**2*Zl12**2+2D0*Zdd22**2/CW**2*Zl11
     &      **2-4D0*Zdd22**2/CW**2*Zl12**2-3D0*Zdd21**2/SW
     &      **2*Zl11**2-Zdd21**2/CW**2*Zl11**2*sb**2+2D0
     &      *Zdd21**2/CW**2*Zl12**2*sb**2-2D0*Zdd22**2/CW
     &      **2*Zl11**2*sb**2+4D0*Zdd22**2/CW**2*Zl12**2
     &      *sb**2+3D0*Zdd21**2/SW**2*Zl11**2*sb**2-12D0
     &      *Ml/MW**2*Mq/SW**2*Zdd21*Zdd22*Zl11*Zl12)
      AAABR(937) = EE**2/cb**2*(Zb11/CW**2*Zb21*Zl11*Zl21
     &      -2D0*Zb11/CW**2*Zb21*Zl12*Zl22+2D0*Zb12/CW**2
     &      *Zb22*Zl11*Zl21-4D0*Zb12/CW**2*Zb22*Zl12*Zl22
     &      -3D0*Zb11/SW**2*Zb21*Zl11*Zl21-Zb11/CW**2*Zb21
     &      *Zl11*Zl21*sb**2+2D0*Zb11/CW**2*Zb21*Zl12*Zl22
     &      *sb**2-2D0*Zb12/CW**2*Zb22*Zl11*Zl21*sb**2+4D0
     &      *Zb12/CW**2*Zb22*Zl12*Zl22*sb**2+3D0*Zb11/SW
     &      **2*Zb21*Zl11*Zl21*sb**2-6D0*MbMM/MW**2*Ml/SW
     &      **2*Zb11*Zb22*Zl11*Zl22-6D0*MbMM/MW**2*Ml/SW
     &      **2*Zb12*Zb21*Zl12*Zl21)
      AAABR(938) = EE**2/cb**2*(Zl11/CW**2*Zl21*Zss11*Zss21
     &      +2D0*Zl11/CW**2*Zl21*Zss12*Zss22-2D0*Zl12/CW
     &      **2*Zl22*Zss11*Zss21-4D0*Zl12/CW**2*Zl22*Zss12
     &      *Zss22-3D0*Zl11/SW**2*Zl21*Zss11*Zss21-Zl11/CW
     &      **2*Zl21*Zss11*Zss21*sb**2-2D0*Zl11/CW**2*Zl21
     &      *Zss12*Zss22*sb**2+2D0*Zl12/CW**2*Zl22*Zss11
     &      *Zss21*sb**2+4D0*Zl12/CW**2*Zl22*Zss12*Zss22
     &      *sb**2+3D0*Zl11/SW**2*Zl21*Zss11*Zss21*sb**2
     &      -6D0*Ml/MW**2*Mq/SW**2*Zl11*Zl22*Zss11*Zss22
     &      -6D0*Ml/MW**2*Mq/SW**2*Zl12*Zl21*Zss12*Zss21)
      AAABR(939) = EE**2/cb**2*(Zdd11/CW**2*Zdd21*Zl11*Zl21
     &      -2D0*Zdd11/CW**2*Zdd21*Zl12*Zl22+2D0*Zdd12/CW
     &      **2*Zdd22*Zl11*Zl21-4D0*Zdd12/CW**2*Zdd22*Zl12
     &      *Zl22-3D0*Zdd11/SW**2*Zdd21*Zl11*Zl21-Zdd11/CW
     &      **2*Zdd21*Zl11*Zl21*sb**2+2D0*Zdd11/CW**2*Zdd21
     &      *Zl12*Zl22*sb**2-2D0*Zdd12/CW**2*Zdd22*Zl11*Zl21
     &      *sb**2+4D0*Zdd12/CW**2*Zdd22*Zl12*Zl22*sb**2
     &      +3D0*Zdd11/SW**2*Zdd21*Zl11*Zl21*sb**2-6D0*Ml
     &      /MW**2*Mq/SW**2*Zdd11*Zdd22*Zl11*Zl22-6D0*Ml
     &      /MW**2*Mq/SW**2*Zdd12*Zdd21*Zl12*Zl21)
      AAABR(940) = EE**2/cb**2*(Zl11/CW**2*Zl21*Zss21**2
     &      +2D0*Zl11/CW**2*Zl21*Zss22**2-2D0*Zl12/CW**2
     &      *Zl22*Zss21**2-4D0*Zl12/CW**2*Zl22*Zss22**2-3D0
     &      *Zl11/SW**2*Zl21*Zss21**2-Zl11/CW**2*Zl21*Zss21
     &      **2*sb**2-2D0*Zl11/CW**2*Zl21*Zss22**2*sb**2
     &      +2D0*Zl12/CW**2*Zl22*Zss21**2*sb**2+4D0*Zl12
     &      /CW**2*Zl22*Zss22**2*sb**2+3D0*Zl11/SW**2*Zl21
     &      *Zss21**2*sb**2-6D0*Ml/MW**2*Mq/SW**2*Zl11*Zl22
     &      *Zss21*Zss22-6D0*Ml/MW**2*Mq/SW**2*Zl12*Zl21
     &      *Zss21*Zss22)
      AAABR(941) = EE**2/cb**2*(Zb21**2/CW**2*Zl11*Zl21-2D0
     &      *Zb21**2/CW**2*Zl12*Zl22+2D0*Zb22**2/CW**2*Zl11
     &      *Zl21-4D0*Zb22**2/CW**2*Zl12*Zl22-3D0*Zb21**2
     &      /SW**2*Zl11*Zl21-Zb21**2/CW**2*Zl11*Zl21*sb**2
     &      +2D0*Zb21**2/CW**2*Zl12*Zl22*sb**2-2D0*Zb22**2
     &      /CW**2*Zl11*Zl21*sb**2+4D0*Zb22**2/CW**2*Zl12
     &      *Zl22*sb**2+3D0*Zb21**2/SW**2*Zl11*Zl21*sb**2
     &      -6D0*MbMM/MW**2*Ml/SW**2*Zb21*Zb22*Zl11*Zl22
     &      -6D0*MbMM/MW**2*Ml/SW**2*Zb21*Zb22*Zl12*Zl21)
      AAABR(942) = EE**2/cb**2*(Zdd21**2/CW**2*Zl11*Zl21
     &      -2D0*Zdd21**2/CW**2*Zl12*Zl22+2D0*Zdd22**2/CW
     &      **2*Zl11*Zl21-4D0*Zdd22**2/CW**2*Zl12*Zl22-3D0
     &      *Zdd21**2/SW**2*Zl11*Zl21-Zdd21**2/CW**2*Zl11
     &      *Zl21*sb**2+2D0*Zdd21**2/CW**2*Zl12*Zl22*sb**2
     &      -2D0*Zdd22**2/CW**2*Zl11*Zl21*sb**2+4D0*Zdd22
     &      **2/CW**2*Zl12*Zl22*sb**2+3D0*Zdd21**2/SW**2
     &      *Zl11*Zl21*sb**2-6D0*Ml/MW**2*Mq/SW**2*Zdd21
     &      *Zdd22*Zl11*Zl22-6D0*Ml/MW**2*Mq/SW**2*Zdd21
     &      *Zdd22*Zl12*Zl21)
      AAABR(943) = EE**2*(Zcc11**2/CW**2-4D0*Zcc12**2/CW
     &      **2+3D0*Zcc11**2/SW**2)
      AAABR(944) = EE**2*(Zuu11**2/CW**2-4D0*Zuu12**2/CW
     &      **2+3D0*Zuu11**2/SW**2)
      AAABR(945) = EE**2*(Zt11**2/CW**2-4D0*Zt12**2/CW**2
     &      +3D0*Zt11**2/SW**2)
      AAABR(946) = EE**2*(Zl11**2/CW**2*Zuu11**2-4D0*Zl11
     &      **2/CW**2*Zuu12**2-2D0*Zl12**2/CW**2*Zuu11**2
     &      +8D0*Zl12**2/CW**2*Zuu12**2+3D0*Zl11**2/SW**2
     &      *Zuu11**2)
      AAABR(947) = EE**2*(Zcc11**2/CW**2*Zl11**2-2D0*Zcc11
     &      **2/CW**2*Zl12**2-4D0*Zcc12**2/CW**2*Zl11**2
     &      +8D0*Zcc12**2/CW**2*Zl12**2+3D0*Zcc11**2/SW**2
     &      *Zl11**2)
      AAABR(948) = EE**2*(Zl11**2/CW**2*Zt11**2-4D0*Zl11
     &      **2/CW**2*Zt12**2-2D0*Zl12**2/CW**2*Zt11**2+8D0
     &      *Zl12**2/CW**2*Zt12**2+3D0*Zl11**2/SW**2*Zt11
     &      **2)
      AAABR(949) = EE**2*(Zl11**2/CW**2*Zt11*Zt21-4D0*Zl11
     &      **2/CW**2*Zt12*Zt22-2D0*Zl12**2/CW**2*Zt11*Zt21
     &      +8D0*Zl12**2/CW**2*Zt12*Zt22+3D0*Zl11**2/SW**2
     &      *Zt11*Zt21)
      AAABR(950) = EE**2*(Zuu11/CW**2*Zuu21-4D0*Zuu12/CW
     &      **2*Zuu22+3D0*Zuu11/SW**2*Zuu21)
      AAABR(951) = EE**2*(Zt11/CW**2*Zt21-4D0*Zt12/CW**2
     &      *Zt22+3D0*Zt11/SW**2*Zt21)
      AAABR(952) = EE**2*(Zl11**2/CW**2*Zuu11*Zuu21-4D0*Zl11
     &      **2/CW**2*Zuu12*Zuu22-2D0*Zl12**2/CW**2*Zuu11
     &      *Zuu21+8D0*Zl12**2/CW**2*Zuu12*Zuu22+3D0*Zl11
     &      **2/SW**2*Zuu11*Zuu21)
      AAABR(953) = EE**2*(Zcc11/CW**2*Zcc21-4D0*Zcc12/CW
     &      **2*Zcc22+3D0*Zcc11/SW**2*Zcc21)
      AAABR(954) = EE**2*(Zcc11/CW**2*Zcc21*Zl11**2-2D0*Zcc11
     &      /CW**2*Zcc21*Zl12**2-4D0*Zcc12/CW**2*Zcc22*Zl11
     &      **2+8D0*Zcc12/CW**2*Zcc22*Zl12**2+3D0*Zcc11/SW
     &      **2*Zcc21*Zl11**2)
      AAABR(955) = EE**2*(Zl11/CW**2*Zl21*Zuu11**2-4D0*Zl11
     &      /CW**2*Zl21*Zuu12**2-2D0*Zl12/CW**2*Zl22*Zuu11
     &      **2+8D0*Zl12/CW**2*Zl22*Zuu12**2+3D0*Zl11/SW
     &      **2*Zl21*Zuu11**2)
      AAABR(956) = EE**2*(Zl11/CW**2*Zl21*Zt11**2-4D0*Zl11
     &      /CW**2*Zl21*Zt12**2-2D0*Zl12/CW**2*Zl22*Zt11
     &      **2+8D0*Zl12/CW**2*Zl22*Zt12**2+3D0*Zl11/SW**2
     &      *Zl21*Zt11**2)
      AAABR(957) = EE**2*(Zcc11**2/CW**2*Zl11*Zl21-2D0*Zcc11
     &      **2/CW**2*Zl12*Zl22-4D0*Zcc12**2/CW**2*Zl11*Zl21
     &      +8D0*Zcc12**2/CW**2*Zl12*Zl22+3D0*Zcc11**2/SW
     &      **2*Zl11*Zl21)
      AAABR(958) = EE**2*(Zcc11/CW**2*Zcc21*Zl11*Zl21-2D0
     &      *Zcc11/CW**2*Zcc21*Zl12*Zl22-4D0*Zcc12/CW**2
     &      *Zcc22*Zl11*Zl21+8D0*Zcc12/CW**2*Zcc22*Zl12*Zl22
     &      +3D0*Zcc11/SW**2*Zcc21*Zl11*Zl21)
      AAABR(959) = EE**2*(Zl11/CW**2*Zl21*Zt11*Zt21-4D0*Zl11
     &      /CW**2*Zl21*Zt12*Zt22-2D0*Zl12/CW**2*Zl22*Zt11
     &      *Zt21+8D0*Zl12/CW**2*Zl22*Zt12*Zt22+3D0*Zl11
     &      /SW**2*Zl21*Zt11*Zt21)
      AAABR(960) = EE**2*(Zl11/CW**2*Zl21*Zuu11*Zuu21-4D0
     &      *Zl11/CW**2*Zl21*Zuu12*Zuu22-2D0*Zl12/CW**2*Zl22
     &      *Zuu11*Zuu21+8D0*Zl12/CW**2*Zl22*Zuu12*Zuu22
     &      +3D0*Zl11/SW**2*Zl21*Zuu11*Zuu21)
      AAABR(961) = EE**2/SW**2*Zt11/cb**2*(Zb11*Zl11-Zb11
     &      *Zl11*sb**2+MbMM/MW**2*Ml*Zb12*Zl12)
      AAABR(962) = EE**2/SW**2*Zcc11*Zss11
      AAABR(963) = EE**2/SW**2*Zb11*Zt11
      AAABR(964) = EE**2/SW**2*Zuu11/cb**2*(Zdd11*Zl11-Zdd11
     &      *Zl11*sb**2+Ml/MW**2*Mq*Zdd12*Zl12)
      AAABR(965) = EE**2/SW**2*Zcc11/cb**2*(Zl11*Zss11-Zl11
     &      *Zss11*sb**2+Ml/MW**2*Mq*Zl12*Zss12)
      AAABR(966) = EE**2/SW**2*Zdd11*Zuu11
      AAABR(967) = EE**2/SW**2*Zdd21*Zuu11
      AAABR(968) = EE**2/SW**2*Zcc11*Zss21
      AAABR(969) = EE**2/SW**2*Zb21*Zt11
      AAABR(970) = EE**2/SW**2*Zt11/cb**2*(Zb21*Zl11-Zb21
     &      *Zl11*sb**2+MbMM/MW**2*Ml*Zb22*Zl12)
      AAABR(971) = EE**2/SW**2*Zuu11/cb**2*(Zdd21*Zl11-Zdd21
     &      *Zl11*sb**2+Ml/MW**2*Mq*Zdd22*Zl12)
      AAABR(972) = EE**2/SW**2*Zcc11/cb**2*(Zl11*Zss21-Zl11
     &      *Zss21*sb**2+Ml/MW**2*Mq*Zl12*Zss22)
      AAABR(973) = EE**2*(Zl11**2/CW**2*Zt21**2-4D0*Zl11
     &      **2/CW**2*Zt22**2-2D0*Zl12**2/CW**2*Zt21**2+8D0
     &      *Zl12**2/CW**2*Zt22**2+3D0*Zl11**2/SW**2*Zt21
     &      **2)
      AAABR(974) = EE**2*(Zcc21**2/CW**2-4D0*Zcc22**2/CW
     &      **2+3D0*Zcc21**2/SW**2)
      AAABR(975) = EE**2*(Zt21**2/CW**2-4D0*Zt22**2/CW**2
     &      +3D0*Zt21**2/SW**2)
      AAABR(976) = EE**2*(Zuu21**2/CW**2-4D0*Zuu22**2/CW
     &      **2+3D0*Zuu21**2/SW**2)
      AAABR(977) = EE**2*(Zcc21**2/CW**2*Zl11**2-2D0*Zcc21
     &      **2/CW**2*Zl12**2-4D0*Zcc22**2/CW**2*Zl11**2
     &      +8D0*Zcc22**2/CW**2*Zl12**2+3D0*Zcc21**2/SW**2
     &      *Zl11**2)
      AAABR(978) = EE**2*(Zl11**2/CW**2*Zuu21**2-4D0*Zl11
     &      **2/CW**2*Zuu22**2-2D0*Zl12**2/CW**2*Zuu21**2
     &      +8D0*Zl12**2/CW**2*Zuu22**2+3D0*Zl11**2/SW**2
     &      *Zuu21**2)
      AAABR(979) = EE**2*(Zcc21**2/CW**2*Zl11*Zl21-2D0*Zcc21
     &      **2/CW**2*Zl12*Zl22-4D0*Zcc22**2/CW**2*Zl11*Zl21
     &      +8D0*Zcc22**2/CW**2*Zl12*Zl22+3D0*Zcc21**2/SW
     &      **2*Zl11*Zl21)
      AAABR(980) = EE**2*(Zl11/CW**2*Zl21*Zt21**2-4D0*Zl11
     &      /CW**2*Zl21*Zt22**2-2D0*Zl12/CW**2*Zl22*Zt21
     &      **2+8D0*Zl12/CW**2*Zl22*Zt22**2+3D0*Zl11/SW**2
     &      *Zl21*Zt21**2)
      AAABR(981) = EE**2*(Zl11/CW**2*Zl21*Zuu21**2-4D0*Zl11
     &      /CW**2*Zl21*Zuu22**2-2D0*Zl12/CW**2*Zl22*Zuu21
     &      **2+8D0*Zl12/CW**2*Zl22*Zuu22**2+3D0*Zl11/SW
     &      **2*Zl21*Zuu21**2)
      AAABR(982) = EE**2/SW**2*Zcc21/cb**2*(Zl11*Zss11-Zl11
     &      *Zss11*sb**2+Ml/MW**2*Mq*Zl12*Zss12)
      AAABR(983) = EE**2/SW**2*Zdd11*Zuu21
      AAABR(984) = EE**2/SW**2*Zb11*Zt21
      AAABR(985) = EE**2/SW**2*Zcc21*Zss11
      AAABR(986) = EE**2/SW**2*Zuu21/cb**2*(Zdd11*Zl11-Zdd11
     &      *Zl11*sb**2+Ml/MW**2*Mq*Zdd12*Zl12)
      AAABR(987) = EE**2/SW**2*Zt21/cb**2*(Zb11*Zl11-Zb11
     &      *Zl11*sb**2+MbMM/MW**2*Ml*Zb12*Zl12)
      AAABR(988) = EE**2/SW**2*Zdd21*Zuu21
      AAABR(989) = EE**2/SW**2*Zuu21/cb**2*(Zdd21*Zl11-Zdd21
     &      *Zl11*sb**2+Ml/MW**2*Mq*Zdd22*Zl12)
      AAABR(990) = EE**2/SW**2*Zt21/cb**2*(Zb21*Zl11-Zb21
     &      *Zl11*sb**2+MbMM/MW**2*Ml*Zb22*Zl12)
      AAABR(991) = EE**2/SW**2*Zb21*Zt21
      AAABR(992) = EE**2/SW**2*Zcc21*Zss21
      AAABR(993) = EE**2/SW**2*Zcc21/cb**2*(Zl11*Zss21-Zl11
     &      *Zss21*sb**2+Ml/MW**2*Mq*Zl12*Zss22)
      AAABR(994) = EE**2*(Zl11**2/CW**2-2D0*Zl12**2/CW**2
     &      +Zl11**2/SW**2-2D0*Zl11**2/CW**2*sb**2+4D0*Zl12
     &      **2/CW**2*sb**2-2D0*Zl11**2/SW**2*sb**2+2D0*Ml
     &      **2/MW**2/SW**2*Zl12**2)
      AAABR(995) = EE**2*(Zl11/CW**2*Zl21-2D0*Zl12/CW**2
     &      *Zl22+Zl11/SW**2*Zl21-2D0*Zl11/CW**2*Zl21*sb
     &      **2+4D0*Zl12/CW**2*Zl22*sb**2-2D0*Zl11/SW**2
     &      *Zl21*sb**2+2D0*Ml**2/MW**2/SW**2*Zl12*Zl22)
      AAABR(996) = EE**2/SW**2*Sqrt2*(1D0-2D0*sb**2)
      AAABR(997) = EE**2/SW**2*Sqrt2*Zl11*(1D0-2D0*sb**2
     &      -Ml**2/MW**2)
      AAABR(998) = EE**2*Pa11/SW**2*Sqrt2*cb*sb
      AAABR(999) = EE**2*Pa21/SW**2*Sqrt2*cb*sb
      AAABR(1000) = tB*(2D0*EE**2*Pa11/SW**2*Sqrt2*Zl11-2D0
     &      *EE**2*Pa11/SW**2*Sqrt2*Zl11*sb**2+2D0*EE*Lambda
     &      /MW*Ml*Pa12/SW*Zl12-EE**2/MW**2*Ml**2*Pa11/SW
     &      **2*Sqrt2*Zl11)
      end

      subroutine aaini02
      implicit none
#include "model.h"

      AAABR(1001) = tB*(2D0*EE**2*Pa21/SW**2*Sqrt2*Zl11-2D0
     &      *EE**2*Pa21/SW**2*Sqrt2*Zl11*sb**2+2D0*EE*Lambda
     &      /MW*Ml*Pa22/SW*Zl12-EE**2/MW**2*Ml**2*Pa21/SW
     &      **2*Sqrt2*Zl11)
      AAABR(1002) = EE**2/SW**2*Sqrt2*(Zh11*cb-Zh12*sb)
      AAABR(1003) = EE**2/SW**2*Sqrt2*(Zh21*cb-Zh22*sb)
      AAABR(1004) = EE**2/SW**2*Sqrt2*(Zh31*cb-Zh32*sb)
      AAABR(1005) = EE**2/SW**2*Sqrt2*Zh21*Zl11/cb-EE**2
     &      /SW**2*Sqrt2*Zh22*Zl11*sb-EE**2/SW**2*Sqrt2*Zh21
     &      *Zl11/cb*sb**2-EE**2/MW**2*Ml**2/SW**2*Sqrt2
     &      *Zh21*Zl11/cb+2D0*EE*Lambda/MW*Ml/SW*Zh23*Zl12
     &      /cb*sb
      AAABR(1006) = EE**2/SW**2*Sqrt2*Zh31*Zl11/cb-EE**2
     &      /SW**2*Sqrt2*Zh32*Zl11*sb-EE**2/SW**2*Sqrt2*Zh31
     &      *Zl11/cb*sb**2-EE**2/MW**2*Ml**2/SW**2*Sqrt2
     &      *Zh31*Zl11/cb+2D0*EE*Lambda/MW*Ml/SW*Zh33*Zl12
     &      /cb*sb
      AAABR(1007) = EE**2/SW**2*Sqrt2*Zh11*Zl11/cb-EE**2
     &      /SW**2*Sqrt2*Zh12*Zl11*sb-EE**2/SW**2*Sqrt2*Zh11
     &      *Zl11/cb*sb**2-EE**2/MW**2*Ml**2/SW**2*Sqrt2
     &      *Zh11*Zl11/cb+2D0*EE*Lambda/MW*Ml/SW*Zh13*Zl12
     &      /cb*sb
      AAABR(1008) = EE**2*(Zl11**2/CW**2-2D0*Zl12**2/CW**2
     &      -Zl11**2/SW**2-2D0*Zl11**2/CW**2*sb**2+4D0*Zl12
     &      **2/CW**2*sb**2+2D0*Ml**2/MW**2/SW**2+2D0*Zl11
     &      **2/SW**2*sb**2)
      AAABR(1009) = EE**2*(Zl11/CW**2*Zl21-2D0*Zl12/CW**2
     &      *Zl22-Zl11/SW**2*Zl21-2D0*Zl11/CW**2*Zl21*sb
     &      **2+4D0*Zl12/CW**2*Zl22*sb**2+2D0*Zl11/SW**2
     &      *Zl21*sb**2)
      AAABR(1010) = EE**2*Pa11*cb*sb*(1D0/CW**2-1D0/SW**2)
      AAABR(1011) = EE**2*Pa21*cb*sb*(1D0/CW**2-1D0/SW**2)
      AAABR(1012) = tB*(EE**2/CW**2*Pa21*Zl11**2-2D0*EE**2
     &      /CW**2*Pa21*Zl12**2-EE**2*Pa21/SW**2*Zl11**2
     &      -EE**2/CW**2*Pa21*Zl11**2*sb**2+2D0*EE**2/CW
     &      **2*Pa21*Zl12**2*sb**2+EE**2/MW**2*Ml**2*Pa21
     &      /SW**2+EE**2*Pa21/SW**2*Zl11**2*sb**2-EE*Lambda
     &      /MW*Ml*Pa22/SW*Sqrt2*Zl11*Zl12)
      AAABR(1013) = tB*(EE**2/CW**2*Pa11*Zl11**2-2D0*EE**2
     &      /CW**2*Pa11*Zl12**2-EE**2*Pa11/SW**2*Zl11**2
     &      -EE**2/CW**2*Pa11*Zl11**2*sb**2+2D0*EE**2/CW
     &      **2*Pa11*Zl12**2*sb**2+EE**2/MW**2*Ml**2*Pa11
     &      /SW**2+EE**2*Pa11/SW**2*Zl11**2*sb**2-EE*Lambda
     &      /MW*Ml*Pa12/SW*Sqrt2*Zl11*Zl12)
      AAABR(1014) = tB*(2D0*EE**2/CW**2*Pa11*Zl11*Zl21-4D0
     &      *EE**2/CW**2*Pa11*Zl12*Zl22-2D0*EE**2*Pa11/SW
     &      **2*Zl11*Zl21-2D0*EE**2/CW**2*Pa11*Zl11*Zl21
     &      *sb**2+4D0*EE**2/CW**2*Pa11*Zl12*Zl22*sb**2+2D0
     &      *EE**2*Pa11/SW**2*Zl11*Zl21*sb**2-EE*Lambda/MW
     &      *Ml*Pa12/SW*Sqrt2*Zl11*Zl22-EE*Lambda/MW*Ml*Pa12
     &      /SW*Sqrt2*Zl12*Zl21)
      AAABR(1015) = tB*(2D0*EE**2/CW**2*Pa21*Zl11*Zl21-4D0
     &      *EE**2/CW**2*Pa21*Zl12*Zl22-2D0*EE**2*Pa21/SW
     &      **2*Zl11*Zl21-2D0*EE**2/CW**2*Pa21*Zl11*Zl21
     &      *sb**2+4D0*EE**2/CW**2*Pa21*Zl12*Zl22*sb**2+2D0
     &      *EE**2*Pa21/SW**2*Zl11*Zl21*sb**2-EE*Lambda/MW
     &      *Ml*Pa22/SW*Sqrt2*Zl11*Zl22-EE*Lambda/MW*Ml*Pa22
     &      /SW*Sqrt2*Zl12*Zl21)
      AAABR(1016) = EE*Lambda/MW*Ml/SW*Sqrt2*Zh33*tB*(Zl11
     &      *Zl22-Zl12*Zl21)
      AAABR(1017) = EE*Lambda/MW*Ml/SW*Sqrt2*Zh23*tB*(Zl11
     &      *Zl22-Zl12*Zl21)
      AAABR(1018) = EE*Lambda/MW*Ml/SW*Sqrt2*Zh13*tB*(Zl11
     &      *Zl22-Zl12*Zl21)
      AAABR(1019) = EE**2*Pa21**2*(1D0/CW**2-1D0/SW**2-2D0
     &      *sb**2/CW**2+2D0*sb**2/SW**2)
      AAABR(1020) = EE**2*Pa11**2*(1D0/CW**2-1D0/SW**2-2D0
     &      *sb**2/CW**2+2D0*sb**2/SW**2)
      AAABR(1021) = EE**2*Pa11*Pa21*(1D0/CW**2-1D0/SW**2
     &      -2D0*sb**2/CW**2+2D0*sb**2/SW**2)
      AAABR(1022) = 1D0/cb**2*(EE**2/CW**2*Pa11*Pa21*Zl11
     &      **2-2D0*EE**2/CW**2*Pa11*Pa21*Zl12**2-EE**2*Pa11
     &      *Pa21/SW**2*Zl11**2-3D0*EE**2/CW**2*Pa11*Pa21
     &      *Zl11**2*sb**2+2D0*EE**2/CW**2*Pa11*Pa21*Zl11
     &      **2*sb**4+6D0*EE**2/CW**2*Pa11*Pa21*Zl12**2*sb
     &      **2-4D0*EE**2/CW**2*Pa11*Pa21*Zl12**2*sb**4+3D0
     &      *EE**2*Pa11*Pa21/SW**2*Zl11**2*sb**2-2D0*EE**2
     &      *Pa11*Pa21/SW**2*Zl11**2*sb**4-2D0*EE**2/MW**2
     &      *Ml**2*Pa11*Pa21/SW**2*sb**2-2D0*EE*Lambda/MW
     &      *Ml*Pa11*Pa22/SW*Sqrt2*Zl11*Zl12-2D0*EE*Lambda
     &      /MW*Ml*Pa12*Pa21/SW*Sqrt2*Zl11*Zl12+2D0*EE*Lambda
     &      /MW*Ml*Pa11*Pa22/SW*Sqrt2*Zl11*Zl12*sb**2+2D0
     &      *EE*Lambda/MW*Ml*Pa12*Pa21/SW*Sqrt2*Zl11*Zl12
     &      *sb**2)
      AAABR(1023) = 1D0/cb**2*(EE**2/CW**2*Pa11**2*Zl11**2
     &      -2D0*EE**2/CW**2*Pa11**2*Zl12**2-EE**2*Pa11**2
     &      /SW**2*Zl11**2-3D0*EE**2/CW**2*Pa11**2*Zl11**2
     &      *sb**2+2D0*EE**2/CW**2*Pa11**2*Zl11**2*sb**4
     &      +6D0*EE**2/CW**2*Pa11**2*Zl12**2*sb**2-4D0*EE
     &      **2/CW**2*Pa11**2*Zl12**2*sb**4+3D0*EE**2*Pa11
     &      **2/SW**2*Zl11**2*sb**2-2D0*EE**2*Pa11**2/SW
     &      **2*Zl11**2*sb**4-2D0*EE**2/MW**2*Ml**2*Pa11
     &      **2/SW**2*sb**2-4D0*EE*Lambda/MW*Ml*Pa11*Pa12
     &      /SW*Sqrt2*Zl11*Zl12+4D0*EE*Lambda/MW*Ml*Pa11
     &      *Pa12/SW*Sqrt2*Zl11*Zl12*sb**2)
      AAABR(1024) = 1D0/cb**2*(EE**2/CW**2*Pa21**2*Zl11**2
     &      -2D0*EE**2/CW**2*Pa21**2*Zl12**2-EE**2*Pa21**2
     &      /SW**2*Zl11**2-3D0*EE**2/CW**2*Pa21**2*Zl11**2
     &      *sb**2+2D0*EE**2/CW**2*Pa21**2*Zl11**2*sb**4
     &      +6D0*EE**2/CW**2*Pa21**2*Zl12**2*sb**2-4D0*EE
     &      **2/CW**2*Pa21**2*Zl12**2*sb**4+3D0*EE**2*Pa21
     &      **2/SW**2*Zl11**2*sb**2-2D0*EE**2*Pa21**2/SW
     &      **2*Zl11**2*sb**4-2D0*EE**2/MW**2*Ml**2*Pa21
     &      **2/SW**2*sb**2-4D0*EE*Lambda/MW*Ml*Pa21*Pa22
     &      /SW*Sqrt2*Zl11*Zl12+4D0*EE*Lambda/MW*Ml*Pa21
     &      *Pa22/SW*Sqrt2*Zl11*Zl12*sb**2)
      AAABR(1025) = 1D0/cb**2*(EE**2/CW**2*Pa11*Pa21*Zl11
     &      *Zl21-2D0*EE**2/CW**2*Pa11*Pa21*Zl12*Zl22-EE
     &      **2*Pa11*Pa21/SW**2*Zl11*Zl21-3D0*EE**2/CW**2
     &      *Pa11*Pa21*Zl11*Zl21*sb**2+2D0*EE**2/CW**2*Pa11
     &      *Pa21*Zl11*Zl21*sb**4+6D0*EE**2/CW**2*Pa11*Pa21
     &      *Zl12*Zl22*sb**2-4D0*EE**2/CW**2*Pa11*Pa21*Zl12
     &      *Zl22*sb**4+3D0*EE**2*Pa11*Pa21/SW**2*Zl11*Zl21
     &      *sb**2-2D0*EE**2*Pa11*Pa21/SW**2*Zl11*Zl21*sb
     &      **4-EE*Lambda/MW*Ml*Pa11*Pa22/SW*Sqrt2*Zl11*Zl22
     &      -EE*Lambda/MW*Ml*Pa11*Pa22/SW*Sqrt2*Zl12*Zl21
     &      -EE*Lambda/MW*Ml*Pa12*Pa21/SW*Sqrt2*Zl11*Zl22
     &      -EE*Lambda/MW*Ml*Pa12*Pa21/SW*Sqrt2*Zl12*Zl21
     &      +EE*Lambda/MW*Ml*Pa11*Pa22/SW*Sqrt2*Zl11*Zl22
     &      *sb**2+EE*Lambda/MW*Ml*Pa11*Pa22/SW*Sqrt2*Zl12
     &      *Zl21*sb**2+EE*Lambda/MW*Ml*Pa12*Pa21/SW*Sqrt2
     &      *Zl11*Zl22*sb**2+EE*Lambda/MW*Ml*Pa12*Pa21/SW
     &      *Sqrt2*Zl12*Zl21*sb**2)
      AAABR(1026) = 1D0/cb**2*(EE**2/CW**2*Pa11**2*Zl11*Zl21
     &      -2D0*EE**2/CW**2*Pa11**2*Zl12*Zl22-EE**2*Pa11
     &      **2/SW**2*Zl11*Zl21-3D0*EE**2/CW**2*Pa11**2*Zl11
     &      *Zl21*sb**2+2D0*EE**2/CW**2*Pa11**2*Zl11*Zl21
     &      *sb**4+6D0*EE**2/CW**2*Pa11**2*Zl12*Zl22*sb**2
     &      -4D0*EE**2/CW**2*Pa11**2*Zl12*Zl22*sb**4+3D0
     &      *EE**2*Pa11**2/SW**2*Zl11*Zl21*sb**2-2D0*EE**2
     &      *Pa11**2/SW**2*Zl11*Zl21*sb**4-2D0*EE*Lambda
     &      /MW*Ml*Pa11*Pa12/SW*Sqrt2*Zl11*Zl22-2D0*EE*Lambda
     &      /MW*Ml*Pa11*Pa12/SW*Sqrt2*Zl12*Zl21+2D0*EE*Lambda
     &      /MW*Ml*Pa11*Pa12/SW*Sqrt2*Zl11*Zl22*sb**2+2D0
     &      *EE*Lambda/MW*Ml*Pa11*Pa12/SW*Sqrt2*Zl12*Zl21
     &      *sb**2)
      AAABR(1027) = 1D0/cb**2*(EE**2/CW**2*Pa21**2*Zl11*Zl21
     &      -2D0*EE**2/CW**2*Pa21**2*Zl12*Zl22-EE**2*Pa21
     &      **2/SW**2*Zl11*Zl21-3D0*EE**2/CW**2*Pa21**2*Zl11
     &      *Zl21*sb**2+2D0*EE**2/CW**2*Pa21**2*Zl11*Zl21
     &      *sb**4+6D0*EE**2/CW**2*Pa21**2*Zl12*Zl22*sb**2
     &      -4D0*EE**2/CW**2*Pa21**2*Zl12*Zl22*sb**4+3D0
     &      *EE**2*Pa21**2/SW**2*Zl11*Zl21*sb**2-2D0*EE**2
     &      *Pa21**2/SW**2*Zl11*Zl21*sb**4-2D0*EE*Lambda
     &      /MW*Ml*Pa21*Pa22/SW*Sqrt2*Zl11*Zl22-2D0*EE*Lambda
     &      /MW*Ml*Pa21*Pa22/SW*Sqrt2*Zl12*Zl21+2D0*EE*Lambda
     &      /MW*Ml*Pa21*Pa22/SW*Sqrt2*Zl11*Zl22*sb**2+2D0
     &      *EE*Lambda/MW*Ml*Pa21*Pa22/SW*Sqrt2*Zl12*Zl21
     &      *sb**2)
      AAABR(1028) = EE*Lambda/MW*Ml/SW*Sqrt2*(Pa11*Zh13*Zl11
     &      *Zl22-Pa11*Zh13*Zl12*Zl21+Pa12*Zh12*Zl11*Zl22
     &      /cb-Pa12*Zh12*Zl12*Zl21/cb)
      AAABR(1029) = EE*Lambda/MW*Ml/SW*Sqrt2*(Pa11*Zh33*Zl11
     &      *Zl22-Pa11*Zh33*Zl12*Zl21+Pa12*Zh32*Zl11*Zl22
     &      /cb-Pa12*Zh32*Zl12*Zl21/cb)
      AAABR(1030) = EE*Lambda/MW*Ml/SW*Sqrt2*(Pa21*Zh23*Zl11
     &      *Zl22-Pa21*Zh23*Zl12*Zl21+Pa22*Zh22*Zl11*Zl22
     &      /cb-Pa22*Zh22*Zl12*Zl21/cb)
      AAABR(1031) = EE*Lambda/MW*Ml/SW*Sqrt2*(Pa21*Zh13*Zl11
     &      *Zl22-Pa21*Zh13*Zl12*Zl21+Pa22*Zh12*Zl11*Zl22
     &      /cb-Pa22*Zh12*Zl12*Zl21/cb)
      AAABR(1032) = EE*Lambda/MW*Ml/SW*Sqrt2*(Pa11*Zh23*Zl11
     &      *Zl22-Pa11*Zh23*Zl12*Zl21+Pa12*Zh22*Zl11*Zl22
     &      /cb-Pa12*Zh22*Zl12*Zl21/cb)
      AAABR(1033) = EE*Lambda/MW*Ml/SW*Sqrt2*(Pa21*Zh33*Zl11
     &      *Zl22-Pa21*Zh33*Zl12*Zl21+Pa22*Zh32*Zl11*Zl22
     &      /cb-Pa22*Zh32*Zl12*Zl21/cb)
      AAABR(1034) = EE**2*(Zh21/CW**2*Zh31-Zh22/CW**2*Zh32
     &      -Zh21/SW**2*Zh31+Zh22/SW**2*Zh32)
      AAABR(1035) = EE**2/CW**2*Zh21*Zh31*Zl11**2/cb**2-2D0
     &      *EE**2/CW**2*Zh21*Zh31*Zl12**2/cb**2-EE**2/CW
     &      **2*Zh22*Zh32*Zl11**2/cb**2+2D0*EE**2/CW**2*Zh22
     &      *Zh32*Zl12**2/cb**2-EE**2/SW**2*Zh21*Zh31*Zl11
     &      **2/cb**2+EE**2/SW**2*Zh22*Zh32*Zl11**2/cb**2
     &      -EE**2/CW**2*Zh21*Zh31*Zl11**2/cb**2*sb**2+2D0
     &      *EE**2/CW**2*Zh21*Zh31*Zl12**2/cb**2*sb**2+EE
     &      **2/CW**2*Zh22*Zh32*Zl11**2/cb**2*sb**2-2D0*EE
     &      **2/CW**2*Zh22*Zh32*Zl12**2/cb**2*sb**2+2D0*EE
     &      **2/MW**2*Ml**2/SW**2*Zh21*Zh31/cb**2+EE**2/SW
     &      **2*Zh21*Zh31*Zl11**2/cb**2*sb**2-EE**2/SW**2
     &      *Zh22*Zh32*Zl11**2/cb**2*sb**2-2D0*EE*Lambda
     &      /MW*Ml/SW*Sqrt2*Zh22*Zh33*Zl11*Zl12/cb-2D0*EE
     &      *Lambda/MW*Ml/SW*Sqrt2*Zh23*Zh32*Zl11*Zl12/cb
      AAABR(1036) = EE**2*(Zh31**2/CW**2-Zh32**2/CW**2-Zh31
     &      **2/SW**2+Zh32**2/SW**2)
      AAABR(1037) = EE**2/CW**2*Zh21**2*Zl11**2/cb**2-2D0
     &      *EE**2/CW**2*Zh21**2*Zl12**2/cb**2-EE**2/CW**2
     &      *Zh22**2*Zl11**2/cb**2+2D0*EE**2/CW**2*Zh22**2
     &      *Zl12**2/cb**2-EE**2/SW**2*Zh21**2*Zl11**2/cb
     &      **2+EE**2/SW**2*Zh22**2*Zl11**2/cb**2-EE**2/CW
     &      **2*Zh21**2*Zl11**2/cb**2*sb**2+2D0*EE**2/CW
     &      **2*Zh21**2*Zl12**2/cb**2*sb**2+EE**2/CW**2*Zh22
     &      **2*Zl11**2/cb**2*sb**2-2D0*EE**2/CW**2*Zh22
     &      **2*Zl12**2/cb**2*sb**2+2D0*EE**2/MW**2*Ml**2
     &      /SW**2*Zh21**2/cb**2+EE**2/SW**2*Zh21**2*Zl11
     &      **2/cb**2*sb**2-EE**2/SW**2*Zh22**2*Zl11**2/cb
     &      **2*sb**2-4D0*EE*Lambda/MW*Ml/SW*Sqrt2*Zh22*Zh23
     &      *Zl11*Zl12/cb
      AAABR(1038) = EE**2*(Zh21**2/CW**2-Zh22**2/CW**2-Zh21
     &      **2/SW**2+Zh22**2/SW**2)
      AAABR(1039) = EE**2/CW**2*Zh31**2*Zl11**2/cb**2-2D0
     &      *EE**2/CW**2*Zh31**2*Zl12**2/cb**2-EE**2/CW**2
     &      *Zh32**2*Zl11**2/cb**2+2D0*EE**2/CW**2*Zh32**2
     &      *Zl12**2/cb**2-EE**2/SW**2*Zh31**2*Zl11**2/cb
     &      **2+EE**2/SW**2*Zh32**2*Zl11**2/cb**2-EE**2/CW
     &      **2*Zh31**2*Zl11**2/cb**2*sb**2+2D0*EE**2/CW
     &      **2*Zh31**2*Zl12**2/cb**2*sb**2+EE**2/CW**2*Zh32
     &      **2*Zl11**2/cb**2*sb**2-2D0*EE**2/CW**2*Zh32
     &      **2*Zl12**2/cb**2*sb**2+2D0*EE**2/MW**2*Ml**2
     &      /SW**2*Zh31**2/cb**2+EE**2/SW**2*Zh31**2*Zl11
     &      **2/cb**2*sb**2-EE**2/SW**2*Zh32**2*Zl11**2/cb
     &      **2*sb**2-4D0*EE*Lambda/MW*Ml/SW*Sqrt2*Zh32*Zh33
     &      *Zl11*Zl12/cb
      AAABR(1040) = EE**2*(Zh11/CW**2*Zh21-Zh12/CW**2*Zh22
     &      -Zh11/SW**2*Zh21+Zh12/SW**2*Zh22)
      AAABR(1041) = EE**2/CW**2*Zh11*Zh21*Zl11**2/cb**2-2D0
     &      *EE**2/CW**2*Zh11*Zh21*Zl12**2/cb**2-EE**2/CW
     &      **2*Zh12*Zh22*Zl11**2/cb**2+2D0*EE**2/CW**2*Zh12
     &      *Zh22*Zl12**2/cb**2-EE**2/SW**2*Zh11*Zh21*Zl11
     &      **2/cb**2+EE**2/SW**2*Zh12*Zh22*Zl11**2/cb**2
     &      -EE**2/CW**2*Zh11*Zh21*Zl11**2/cb**2*sb**2+2D0
     &      *EE**2/CW**2*Zh11*Zh21*Zl12**2/cb**2*sb**2+EE
     &      **2/CW**2*Zh12*Zh22*Zl11**2/cb**2*sb**2-2D0*EE
     &      **2/CW**2*Zh12*Zh22*Zl12**2/cb**2*sb**2+2D0*EE
     &      **2/MW**2*Ml**2/SW**2*Zh11*Zh21/cb**2+EE**2/SW
     &      **2*Zh11*Zh21*Zl11**2/cb**2*sb**2-EE**2/SW**2
     &      *Zh12*Zh22*Zl11**2/cb**2*sb**2-2D0*EE*Lambda
     &      /MW*Ml/SW*Sqrt2*Zh12*Zh23*Zl11*Zl12/cb-2D0*EE
     &      *Lambda/MW*Ml/SW*Sqrt2*Zh13*Zh22*Zl11*Zl12/cb
      AAABR(1042) = EE**2*(Zh11/CW**2*Zh31-Zh12/CW**2*Zh32
     &      -Zh11/SW**2*Zh31+Zh12/SW**2*Zh32)
      AAABR(1043) = EE**2*(Zh11**2/CW**2-Zh12**2/CW**2-Zh11
     &      **2/SW**2+Zh12**2/SW**2)
      AAABR(1044) = EE**2/CW**2*Zh11*Zh31*Zl11**2/cb**2-2D0
     &      *EE**2/CW**2*Zh11*Zh31*Zl12**2/cb**2-EE**2/CW
     &      **2*Zh12*Zh32*Zl11**2/cb**2+2D0*EE**2/CW**2*Zh12
     &      *Zh32*Zl12**2/cb**2-EE**2/SW**2*Zh11*Zh31*Zl11
     &      **2/cb**2+EE**2/SW**2*Zh12*Zh32*Zl11**2/cb**2
     &      -EE**2/CW**2*Zh11*Zh31*Zl11**2/cb**2*sb**2+2D0
     &      *EE**2/CW**2*Zh11*Zh31*Zl12**2/cb**2*sb**2+EE
     &      **2/CW**2*Zh12*Zh32*Zl11**2/cb**2*sb**2-2D0*EE
     &      **2/CW**2*Zh12*Zh32*Zl12**2/cb**2*sb**2+2D0*EE
     &      **2/MW**2*Ml**2/SW**2*Zh11*Zh31/cb**2+EE**2/SW
     &      **2*Zh11*Zh31*Zl11**2/cb**2*sb**2-EE**2/SW**2
     &      *Zh12*Zh32*Zl11**2/cb**2*sb**2-2D0*EE*Lambda
     &      /MW*Ml/SW*Sqrt2*Zh12*Zh33*Zl11*Zl12/cb-2D0*EE
     &      *Lambda/MW*Ml/SW*Sqrt2*Zh13*Zh32*Zl11*Zl12/cb
      AAABR(1045) = EE**2/CW**2*Zh11**2*Zl11**2/cb**2-2D0
     &      *EE**2/CW**2*Zh11**2*Zl12**2/cb**2-EE**2/CW**2
     &      *Zh12**2*Zl11**2/cb**2+2D0*EE**2/CW**2*Zh12**2
     &      *Zl12**2/cb**2-EE**2/SW**2*Zh11**2*Zl11**2/cb
     &      **2+EE**2/SW**2*Zh12**2*Zl11**2/cb**2-EE**2/CW
     &      **2*Zh11**2*Zl11**2/cb**2*sb**2+2D0*EE**2/CW
     &      **2*Zh11**2*Zl12**2/cb**2*sb**2+EE**2/CW**2*Zh12
     &      **2*Zl11**2/cb**2*sb**2-2D0*EE**2/CW**2*Zh12
     &      **2*Zl12**2/cb**2*sb**2+2D0*EE**2/MW**2*Ml**2
     &      /SW**2*Zh11**2/cb**2+EE**2/SW**2*Zh11**2*Zl11
     &      **2/cb**2*sb**2-EE**2/SW**2*Zh12**2*Zl11**2/cb
     &      **2*sb**2-4D0*EE*Lambda/MW*Ml/SW*Sqrt2*Zh12*Zh13
     &      *Zl11*Zl12/cb
      AAABR(1046) = EE**2/CW**2*Zh11*Zh21*Zl11*Zl21/cb**2
     &      -2D0*EE**2/CW**2*Zh11*Zh21*Zl12*Zl22/cb**2-EE
     &      **2/CW**2*Zh12*Zh22*Zl11*Zl21/cb**2+2D0*EE**2
     &      /CW**2*Zh12*Zh22*Zl12*Zl22/cb**2-EE**2/SW**2
     &      *Zh11*Zh21*Zl11*Zl21/cb**2+EE**2/SW**2*Zh12*Zh22
     &      *Zl11*Zl21/cb**2-EE**2/CW**2*Zh11*Zh21*Zl11*Zl21
     &      /cb**2*sb**2+2D0*EE**2/CW**2*Zh11*Zh21*Zl12*Zl22
     &      /cb**2*sb**2+EE**2/CW**2*Zh12*Zh22*Zl11*Zl21
     &      /cb**2*sb**2-2D0*EE**2/CW**2*Zh12*Zh22*Zl12*Zl22
     &      /cb**2*sb**2+EE**2/SW**2*Zh11*Zh21*Zl11*Zl21
     &      /cb**2*sb**2-EE**2/SW**2*Zh12*Zh22*Zl11*Zl21
     &      /cb**2*sb**2-EE*Lambda/MW*Ml/SW*Sqrt2*Zh12*Zh23
     &      *Zl11*Zl22/cb-EE*Lambda/MW*Ml/SW*Sqrt2*Zh12*Zh23
     &      *Zl12*Zl21/cb-EE*Lambda/MW*Ml/SW*Sqrt2*Zh13*Zh22
     &      *Zl11*Zl22/cb-EE*Lambda/MW*Ml/SW*Sqrt2*Zh13*Zh22
     &      *Zl12*Zl21/cb
      AAABR(1047) = EE**2/CW**2*Zh11*Zh31*Zl11*Zl21/cb**2
     &      -2D0*EE**2/CW**2*Zh11*Zh31*Zl12*Zl22/cb**2-EE
     &      **2/CW**2*Zh12*Zh32*Zl11*Zl21/cb**2+2D0*EE**2
     &      /CW**2*Zh12*Zh32*Zl12*Zl22/cb**2-EE**2/SW**2
     &      *Zh11*Zh31*Zl11*Zl21/cb**2+EE**2/SW**2*Zh12*Zh32
     &      *Zl11*Zl21/cb**2-EE**2/CW**2*Zh11*Zh31*Zl11*Zl21
     &      /cb**2*sb**2+2D0*EE**2/CW**2*Zh11*Zh31*Zl12*Zl22
     &      /cb**2*sb**2+EE**2/CW**2*Zh12*Zh32*Zl11*Zl21
     &      /cb**2*sb**2-2D0*EE**2/CW**2*Zh12*Zh32*Zl12*Zl22
     &      /cb**2*sb**2+EE**2/SW**2*Zh11*Zh31*Zl11*Zl21
     &      /cb**2*sb**2-EE**2/SW**2*Zh12*Zh32*Zl11*Zl21
     &      /cb**2*sb**2-EE*Lambda/MW*Ml/SW*Sqrt2*Zh12*Zh33
     &      *Zl11*Zl22/cb-EE*Lambda/MW*Ml/SW*Sqrt2*Zh12*Zh33
     &      *Zl12*Zl21/cb-EE*Lambda/MW*Ml/SW*Sqrt2*Zh13*Zh32
     &      *Zl11*Zl22/cb-EE*Lambda/MW*Ml/SW*Sqrt2*Zh13*Zh32
     &      *Zl12*Zl21/cb
      AAABR(1048) = EE**2/CW**2*Zh31**2*Zl11*Zl21/cb**2-2D0
     &      *EE**2/CW**2*Zh31**2*Zl12*Zl22/cb**2-EE**2/CW
     &      **2*Zh32**2*Zl11*Zl21/cb**2+2D0*EE**2/CW**2*Zh32
     &      **2*Zl12*Zl22/cb**2-EE**2/SW**2*Zh31**2*Zl11
     &      *Zl21/cb**2+EE**2/SW**2*Zh32**2*Zl11*Zl21/cb
     &      **2-EE**2/CW**2*Zh31**2*Zl11*Zl21/cb**2*sb**2
     &      +2D0*EE**2/CW**2*Zh31**2*Zl12*Zl22/cb**2*sb**2
     &      +EE**2/CW**2*Zh32**2*Zl11*Zl21/cb**2*sb**2-2D0
     &      *EE**2/CW**2*Zh32**2*Zl12*Zl22/cb**2*sb**2+EE
     &      **2/SW**2*Zh31**2*Zl11*Zl21/cb**2*sb**2-EE**2
     &      /SW**2*Zh32**2*Zl11*Zl21/cb**2*sb**2-2D0*EE*Lambda
     &      /MW*Ml/SW*Sqrt2*Zh32*Zh33*Zl11*Zl22/cb-2D0*EE
     &      *Lambda/MW*Ml/SW*Sqrt2*Zh32*Zh33*Zl12*Zl21/cb
      AAABR(1049) = EE**2/CW**2*Zh11**2*Zl11*Zl21/cb**2-2D0
     &      *EE**2/CW**2*Zh11**2*Zl12*Zl22/cb**2-EE**2/CW
     &      **2*Zh12**2*Zl11*Zl21/cb**2+2D0*EE**2/CW**2*Zh12
     &      **2*Zl12*Zl22/cb**2-EE**2/SW**2*Zh11**2*Zl11
     &      *Zl21/cb**2+EE**2/SW**2*Zh12**2*Zl11*Zl21/cb
     &      **2-EE**2/CW**2*Zh11**2*Zl11*Zl21/cb**2*sb**2
     &      +2D0*EE**2/CW**2*Zh11**2*Zl12*Zl22/cb**2*sb**2
     &      +EE**2/CW**2*Zh12**2*Zl11*Zl21/cb**2*sb**2-2D0
     &      *EE**2/CW**2*Zh12**2*Zl12*Zl22/cb**2*sb**2+EE
     &      **2/SW**2*Zh11**2*Zl11*Zl21/cb**2*sb**2-EE**2
     &      /SW**2*Zh12**2*Zl11*Zl21/cb**2*sb**2-2D0*EE*Lambda
     &      /MW*Ml/SW*Sqrt2*Zh12*Zh13*Zl11*Zl22/cb-2D0*EE
     &      *Lambda/MW*Ml/SW*Sqrt2*Zh12*Zh13*Zl12*Zl21/cb
      AAABR(1050) = EE**2/CW**2*Zh21**2*Zl11*Zl21/cb**2-2D0
     &      *EE**2/CW**2*Zh21**2*Zl12*Zl22/cb**2-EE**2/CW
     &      **2*Zh22**2*Zl11*Zl21/cb**2+2D0*EE**2/CW**2*Zh22
     &      **2*Zl12*Zl22/cb**2-EE**2/SW**2*Zh21**2*Zl11
     &      *Zl21/cb**2+EE**2/SW**2*Zh22**2*Zl11*Zl21/cb
     &      **2-EE**2/CW**2*Zh21**2*Zl11*Zl21/cb**2*sb**2
     &      +2D0*EE**2/CW**2*Zh21**2*Zl12*Zl22/cb**2*sb**2
     &      +EE**2/CW**2*Zh22**2*Zl11*Zl21/cb**2*sb**2-2D0
     &      *EE**2/CW**2*Zh22**2*Zl12*Zl22/cb**2*sb**2+EE
     &      **2/SW**2*Zh21**2*Zl11*Zl21/cb**2*sb**2-EE**2
     &      /SW**2*Zh22**2*Zl11*Zl21/cb**2*sb**2-2D0*EE*Lambda
     &      /MW*Ml/SW*Sqrt2*Zh22*Zh23*Zl11*Zl22/cb-2D0*EE
     &      *Lambda/MW*Ml/SW*Sqrt2*Zh22*Zh23*Zl12*Zl21/cb
      AAABR(1051) = EE**2/CW**2*Zh21*Zh31*Zl11*Zl21/cb**2
     &      -2D0*EE**2/CW**2*Zh21*Zh31*Zl12*Zl22/cb**2-EE
     &      **2/CW**2*Zh22*Zh32*Zl11*Zl21/cb**2+2D0*EE**2
     &      /CW**2*Zh22*Zh32*Zl12*Zl22/cb**2-EE**2/SW**2
     &      *Zh21*Zh31*Zl11*Zl21/cb**2+EE**2/SW**2*Zh22*Zh32
     &      *Zl11*Zl21/cb**2-EE**2/CW**2*Zh21*Zh31*Zl11*Zl21
     &      /cb**2*sb**2+2D0*EE**2/CW**2*Zh21*Zh31*Zl12*Zl22
     &      /cb**2*sb**2+EE**2/CW**2*Zh22*Zh32*Zl11*Zl21
     &      /cb**2*sb**2-2D0*EE**2/CW**2*Zh22*Zh32*Zl12*Zl22
     &      /cb**2*sb**2+EE**2/SW**2*Zh21*Zh31*Zl11*Zl21
     &      /cb**2*sb**2-EE**2/SW**2*Zh22*Zh32*Zl11*Zl21
     &      /cb**2*sb**2-EE*Lambda/MW*Ml/SW*Sqrt2*Zh22*Zh33
     &      *Zl11*Zl22/cb-EE*Lambda/MW*Ml/SW*Sqrt2*Zh22*Zh33
     &      *Zl12*Zl21/cb-EE*Lambda/MW*Ml/SW*Sqrt2*Zh23*Zh32
     &      *Zl11*Zl22/cb-EE*Lambda/MW*Ml/SW*Sqrt2*Zh23*Zh32
     &      *Zl12*Zl21/cb
      AAABR(1052) = EE**2/cb**2*(Zl21**4/CW**2+4D0*Zl22**4
     &      /CW**2+Zl21**4/SW**2-4D0*Zl21**2/CW**2*Zl22**2
     &      -Zl21**4/CW**2*sb**2-4D0*Zl22**4/CW**2*sb**2
     &      -Zl21**4/SW**2*sb**2+4D0*Zl21**2/CW**2*Zl22**2
     &      *sb**2+4D0*Ml**2/MW**2/SW**2*Zl21**2*Zl22**2)
      AAABR(1053) = EE**2/CW**2*(Zl21**2-2D0*Zl22**2)
      AAABR(1054) = EE**2/cb**2*(Zl21**2/CW**2-2D0*Zl22**2
     &      /CW**2+Zl21**2/SW**2-Zl21**2/CW**2*sb**2+2D0
     &      *Zl22**2/CW**2*sb**2-Zl21**2/SW**2*sb**2+2D0
     &      *Ml**2/MW**2/SW**2*Zl22**2)
      AAABR(1055) = EE**2*(Zl21**2/CW**2-2D0*Zl22**2/CW**2
     &      -Zl21**2/SW**2)
      AAABR(1056) = EE**2/CW**2*(Zdd11**2+2D0*Zdd12**2)
      AAABR(1057) = EE**2/CW**2*(Zb11**2+2D0*Zb12**2)
      AAABR(1058) = EE**2/cb**2*(Zdd11**2/CW**2*Zl21**2-2D0
     &      *Zdd11**2/CW**2*Zl22**2+2D0*Zdd12**2/CW**2*Zl21
     &      **2-4D0*Zdd12**2/CW**2*Zl22**2-3D0*Zdd11**2/SW
     &      **2*Zl21**2-Zdd11**2/CW**2*Zl21**2*sb**2+2D0
     &      *Zdd11**2/CW**2*Zl22**2*sb**2-2D0*Zdd12**2/CW
     &      **2*Zl21**2*sb**2+4D0*Zdd12**2/CW**2*Zl22**2
     &      *sb**2+3D0*Zdd11**2/SW**2*Zl21**2*sb**2-12D0
     &      *Ml/MW**2*Mq/SW**2*Zdd11*Zdd12*Zl21*Zl22)
      AAABR(1059) = EE**2/cb**2*(Zb11**2/CW**2*Zl21**2-2D0
     &      *Zb11**2/CW**2*Zl22**2+2D0*Zb12**2/CW**2*Zl21
     &      **2-4D0*Zb12**2/CW**2*Zl22**2-3D0*Zb11**2/SW
     &      **2*Zl21**2-Zb11**2/CW**2*Zl21**2*sb**2+2D0*Zb11
     &      **2/CW**2*Zl22**2*sb**2-2D0*Zb12**2/CW**2*Zl21
     &      **2*sb**2+4D0*Zb12**2/CW**2*Zl22**2*sb**2+3D0
     &      *Zb11**2/SW**2*Zl21**2*sb**2-12D0*MbMM/MW**2
     &      *Ml/SW**2*Zb11*Zb12*Zl21*Zl22)
      AAABR(1060) = EE**2/CW**2*(Zss11**2+2D0*Zss12**2)
      AAABR(1061) = EE**2/cb**2*(Zl21**2/CW**2*Zss11**2+2D0
     &      *Zl21**2/CW**2*Zss12**2-2D0*Zl22**2/CW**2*Zss11
     &      **2-4D0*Zl22**2/CW**2*Zss12**2-3D0*Zl21**2/SW
     &      **2*Zss11**2-Zl21**2/CW**2*Zss11**2*sb**2-2D0
     &      *Zl21**2/CW**2*Zss12**2*sb**2+2D0*Zl22**2/CW
     &      **2*Zss11**2*sb**2+4D0*Zl22**2/CW**2*Zss12**2
     &      *sb**2+3D0*Zl21**2/SW**2*Zss11**2*sb**2-12D0
     &      *Ml/MW**2*Mq/SW**2*Zl21*Zl22*Zss11*Zss12)
      AAABR(1062) = EE**2/CW**2*(Zdd11*Zdd21+2D0*Zdd12*Zdd22
     &      )
      AAABR(1063) = EE**2/CW**2*(Zss11*Zss21+2D0*Zss12*Zss22
     &      )
      AAABR(1064) = EE**2/CW**2*(Zb11*Zb21+2D0*Zb12*Zb22)
      AAABR(1065) = EE**2/cb**2*(Zl21**2/CW**2*Zss11*Zss21
     &      +2D0*Zl21**2/CW**2*Zss12*Zss22-2D0*Zl22**2/CW
     &      **2*Zss11*Zss21-4D0*Zl22**2/CW**2*Zss12*Zss22
     &      -3D0*Zl21**2/SW**2*Zss11*Zss21-Zl21**2/CW**2
     &      *Zss11*Zss21*sb**2-2D0*Zl21**2/CW**2*Zss12*Zss22
     &      *sb**2+2D0*Zl22**2/CW**2*Zss11*Zss21*sb**2+4D0
     &      *Zl22**2/CW**2*Zss12*Zss22*sb**2+3D0*Zl21**2
     &      /SW**2*Zss11*Zss21*sb**2-6D0*Ml/MW**2*Mq/SW**2
     &      *Zl21*Zl22*Zss11*Zss22-6D0*Ml/MW**2*Mq/SW**2
     &      *Zl21*Zl22*Zss12*Zss21)
      AAABR(1066) = EE**2/cb**2*(Zdd11/CW**2*Zdd21*Zl21**2
     &      -2D0*Zdd11/CW**2*Zdd21*Zl22**2+2D0*Zdd12/CW**2
     &      *Zdd22*Zl21**2-4D0*Zdd12/CW**2*Zdd22*Zl22**2
     &      -3D0*Zdd11/SW**2*Zdd21*Zl21**2-Zdd11/CW**2*Zdd21
     &      *Zl21**2*sb**2+2D0*Zdd11/CW**2*Zdd21*Zl22**2
     &      *sb**2-2D0*Zdd12/CW**2*Zdd22*Zl21**2*sb**2+4D0
     &      *Zdd12/CW**2*Zdd22*Zl22**2*sb**2+3D0*Zdd11/SW
     &      **2*Zdd21*Zl21**2*sb**2-6D0*Ml/MW**2*Mq/SW**2
     &      *Zdd11*Zdd22*Zl21*Zl22-6D0*Ml/MW**2*Mq/SW**2
     &      *Zdd12*Zdd21*Zl21*Zl22)
      AAABR(1067) = EE**2/cb**2*(Zb11/CW**2*Zb21*Zl21**2
     &      -2D0*Zb11/CW**2*Zb21*Zl22**2+2D0*Zb12/CW**2*Zb22
     &      *Zl21**2-4D0*Zb12/CW**2*Zb22*Zl22**2-3D0*Zb11
     &      /SW**2*Zb21*Zl21**2-Zb11/CW**2*Zb21*Zl21**2*sb
     &      **2+2D0*Zb11/CW**2*Zb21*Zl22**2*sb**2-2D0*Zb12
     &      /CW**2*Zb22*Zl21**2*sb**2+4D0*Zb12/CW**2*Zb22
     &      *Zl22**2*sb**2+3D0*Zb11/SW**2*Zb21*Zl21**2*sb
     &      **2-6D0*MbMM/MW**2*Ml/SW**2*Zb11*Zb22*Zl21*Zl22
     &      -6D0*MbMM/MW**2*Ml/SW**2*Zb12*Zb21*Zl21*Zl22)
      AAABR(1068) = EE**2/cb**2*(Zl21**2/CW**2*Zss21**2+2D0
     &      *Zl21**2/CW**2*Zss22**2-2D0*Zl22**2/CW**2*Zss21
     &      **2-4D0*Zl22**2/CW**2*Zss22**2-3D0*Zl21**2/SW
     &      **2*Zss21**2-Zl21**2/CW**2*Zss21**2*sb**2-2D0
     &      *Zl21**2/CW**2*Zss22**2*sb**2+2D0*Zl22**2/CW
     &      **2*Zss21**2*sb**2+4D0*Zl22**2/CW**2*Zss22**2
     &      *sb**2+3D0*Zl21**2/SW**2*Zss21**2*sb**2-12D0
     &      *Ml/MW**2*Mq/SW**2*Zl21*Zl22*Zss21*Zss22)
      AAABR(1069) = EE**2/CW**2*(Zb21**2+2D0*Zb22**2)
      AAABR(1070) = EE**2/CW**2*(Zdd21**2+2D0*Zdd22**2)
      AAABR(1071) = EE**2/CW**2*(Zss21**2+2D0*Zss22**2)
      AAABR(1072) = EE**2/cb**2*(Zb21**2/CW**2*Zl21**2-2D0
     &      *Zb21**2/CW**2*Zl22**2+2D0*Zb22**2/CW**2*Zl21
     &      **2-4D0*Zb22**2/CW**2*Zl22**2-3D0*Zb21**2/SW
     &      **2*Zl21**2-Zb21**2/CW**2*Zl21**2*sb**2+2D0*Zb21
     &      **2/CW**2*Zl22**2*sb**2-2D0*Zb22**2/CW**2*Zl21
     &      **2*sb**2+4D0*Zb22**2/CW**2*Zl22**2*sb**2+3D0
     &      *Zb21**2/SW**2*Zl21**2*sb**2-12D0*MbMM/MW**2
     &      *Ml/SW**2*Zb21*Zb22*Zl21*Zl22)
      AAABR(1073) = EE**2/cb**2*(Zdd21**2/CW**2*Zl21**2-2D0
     &      *Zdd21**2/CW**2*Zl22**2+2D0*Zdd22**2/CW**2*Zl21
     &      **2-4D0*Zdd22**2/CW**2*Zl22**2-3D0*Zdd21**2/SW
     &      **2*Zl21**2-Zdd21**2/CW**2*Zl21**2*sb**2+2D0
     &      *Zdd21**2/CW**2*Zl22**2*sb**2-2D0*Zdd22**2/CW
     &      **2*Zl21**2*sb**2+4D0*Zdd22**2/CW**2*Zl22**2
     &      *sb**2+3D0*Zdd21**2/SW**2*Zl21**2*sb**2-12D0
     &      *Ml/MW**2*Mq/SW**2*Zdd21*Zdd22*Zl21*Zl22)
      AAABR(1074) = EE**2/CW**2*(Zt11**2-4D0*Zt12**2)
      AAABR(1075) = EE**2/CW**2*(Zcc11**2-4D0*Zcc12**2)
      AAABR(1076) = EE**2*(Zl21**2/CW**2*Zuu11**2-4D0*Zl21
     &      **2/CW**2*Zuu12**2-2D0*Zl22**2/CW**2*Zuu11**2
     &      +8D0*Zl22**2/CW**2*Zuu12**2+3D0*Zl21**2/SW**2
     &      *Zuu11**2)
      AAABR(1077) = EE**2/CW**2*(Zuu11**2-4D0*Zuu12**2)
      AAABR(1078) = EE**2*(Zcc11**2/CW**2*Zl21**2-2D0*Zcc11
     &      **2/CW**2*Zl22**2-4D0*Zcc12**2/CW**2*Zl21**2
     &      +8D0*Zcc12**2/CW**2*Zl22**2+3D0*Zcc11**2/SW**2
     &      *Zl21**2)
      AAABR(1079) = EE**2*(Zl21**2/CW**2*Zt11**2-4D0*Zl21
     &      **2/CW**2*Zt12**2-2D0*Zl22**2/CW**2*Zt11**2+8D0
     &      *Zl22**2/CW**2*Zt12**2+3D0*Zl21**2/SW**2*Zt11
     &      **2)
      AAABR(1080) = EE**2/CW**2*(Zt11*Zt21-4D0*Zt12*Zt22)
      AAABR(1081) = EE**2/CW**2*(Zcc11*Zcc21-4D0*Zcc12*Zcc22
     &      )
      AAABR(1082) = EE**2/CW**2*(Zuu11*Zuu21-4D0*Zuu12*Zuu22
     &      )
      AAABR(1083) = EE**2*(Zcc11/CW**2*Zcc21*Zl21**2-2D0
     &      *Zcc11/CW**2*Zcc21*Zl22**2-4D0*Zcc12/CW**2*Zcc22
     &      *Zl21**2+8D0*Zcc12/CW**2*Zcc22*Zl22**2+3D0*Zcc11
     &      /SW**2*Zcc21*Zl21**2)
      AAABR(1084) = EE**2*(Zl21**2/CW**2*Zuu11*Zuu21-4D0
     &      *Zl21**2/CW**2*Zuu12*Zuu22-2D0*Zl22**2/CW**2
     &      *Zuu11*Zuu21+8D0*Zl22**2/CW**2*Zuu12*Zuu22+3D0
     &      *Zl21**2/SW**2*Zuu11*Zuu21)
      AAABR(1085) = EE**2*(Zl21**2/CW**2*Zt11*Zt21-4D0*Zl21
     &      **2/CW**2*Zt12*Zt22-2D0*Zl22**2/CW**2*Zt11*Zt21
     &      +8D0*Zl22**2/CW**2*Zt12*Zt22+3D0*Zl21**2/SW**2
     &      *Zt11*Zt21)
      AAABR(1086) = EE**2/SW**2*Zcc11/cb**2*(Zl21*Zss11-Zl21
     &      *Zss11*sb**2+Ml/MW**2*Mq*Zl22*Zss12)
      AAABR(1087) = EE**2/SW**2*Zt11/cb**2*(Zb11*Zl21-Zb11
     &      *Zl21*sb**2+MbMM/MW**2*Ml*Zb12*Zl22)
      AAABR(1088) = EE**2/SW**2*Zuu11/cb**2*(Zdd11*Zl21-Zdd11
     &      *Zl21*sb**2+Ml/MW**2*Mq*Zdd12*Zl22)
      AAABR(1089) = EE**2/SW**2*Zuu11/cb**2*(Zdd21*Zl21-Zdd21
     &      *Zl21*sb**2+Ml/MW**2*Mq*Zdd22*Zl22)
      AAABR(1090) = EE**2/SW**2*Zcc11/cb**2*(Zl21*Zss21-Zl21
     &      *Zss21*sb**2+Ml/MW**2*Mq*Zl22*Zss22)
      AAABR(1091) = EE**2/SW**2*Zt11/cb**2*(Zb21*Zl21-Zb21
     &      *Zl21*sb**2+MbMM/MW**2*Ml*Zb22*Zl22)
      AAABR(1092) = EE**2/CW**2*(Zt21**2-4D0*Zt22**2)
      AAABR(1093) = EE**2/CW**2*(Zcc21**2-4D0*Zcc22**2)
      AAABR(1094) = EE**2*(Zl21**2/CW**2*Zt21**2-4D0*Zl21
     &      **2/CW**2*Zt22**2-2D0*Zl22**2/CW**2*Zt21**2+8D0
     &      *Zl22**2/CW**2*Zt22**2+3D0*Zl21**2/SW**2*Zt21
     &      **2)
      AAABR(1095) = EE**2/CW**2*(Zuu21**2-4D0*Zuu22**2)
      AAABR(1096) = EE**2*(Zcc21**2/CW**2*Zl21**2-2D0*Zcc21
     &      **2/CW**2*Zl22**2-4D0*Zcc22**2/CW**2*Zl21**2
     &      +8D0*Zcc22**2/CW**2*Zl22**2+3D0*Zcc21**2/SW**2
     &      *Zl21**2)
      AAABR(1097) = EE**2*(Zl21**2/CW**2*Zuu21**2-4D0*Zl21
     &      **2/CW**2*Zuu22**2-2D0*Zl22**2/CW**2*Zuu21**2
     &      +8D0*Zl22**2/CW**2*Zuu22**2+3D0*Zl21**2/SW**2
     &      *Zuu21**2)
      AAABR(1098) = EE**2/SW**2*Zcc21/cb**2*(Zl21*Zss11-Zl21
     &      *Zss11*sb**2+Ml/MW**2*Mq*Zl22*Zss12)
      AAABR(1099) = EE**2/SW**2*Zt21/cb**2*(Zb11*Zl21-Zb11
     &      *Zl21*sb**2+MbMM/MW**2*Ml*Zb12*Zl22)
      AAABR(1100) = EE**2/SW**2*Zuu21/cb**2*(Zdd11*Zl21-Zdd11
     &      *Zl21*sb**2+Ml/MW**2*Mq*Zdd12*Zl22)
      AAABR(1101) = EE**2/SW**2*Zt21/cb**2*(Zb21*Zl21-Zb21
     &      *Zl21*sb**2+MbMM/MW**2*Ml*Zb22*Zl22)
      AAABR(1102) = EE**2/SW**2*Zuu21/cb**2*(Zdd21*Zl21-Zdd21
     &      *Zl21*sb**2+Ml/MW**2*Mq*Zdd22*Zl22)
      AAABR(1103) = EE**2/SW**2*Zcc21/cb**2*(Zl21*Zss21-Zl21
     &      *Zss21*sb**2+Ml/MW**2*Mq*Zl22*Zss22)
      AAABR(1104) = EE**2*(Zl21**2/CW**2-2D0*Zl22**2/CW**2
     &      +Zl21**2/SW**2-2D0*Zl21**2/CW**2*sb**2+4D0*Zl22
     &      **2/CW**2*sb**2-2D0*Zl21**2/SW**2*sb**2+2D0*Ml
     &      **2/MW**2/SW**2*Zl22**2)
      AAABR(1105) = EE**2/SW**2*Sqrt2*Zl21*(1D0-2D0*sb**2
     &      -Ml**2/MW**2)
      AAABR(1106) = tB*(2D0*EE**2*Pa11/SW**2*Sqrt2*Zl21-2D0
     &      *EE**2*Pa11/SW**2*Sqrt2*Zl21*sb**2+2D0*EE*Lambda
     &      /MW*Ml*Pa12/SW*Zl22-EE**2/MW**2*Ml**2*Pa11/SW
     &      **2*Sqrt2*Zl21)
      AAABR(1107) = tB*(2D0*EE**2*Pa21/SW**2*Sqrt2*Zl21-2D0
     &      *EE**2*Pa21/SW**2*Sqrt2*Zl21*sb**2+2D0*EE*Lambda
     &      /MW*Ml*Pa22/SW*Zl22-EE**2/MW**2*Ml**2*Pa21/SW
     &      **2*Sqrt2*Zl21)
      AAABR(1108) = EE**2/SW**2*Sqrt2*Zh21*Zl21/cb-EE**2
     &      /SW**2*Sqrt2*Zh22*Zl21*sb-EE**2/SW**2*Sqrt2*Zh21
     &      *Zl21/cb*sb**2-EE**2/MW**2*Ml**2/SW**2*Sqrt2
     &      *Zh21*Zl21/cb+2D0*EE*Lambda/MW*Ml/SW*Zh23*Zl22
     &      /cb*sb
      AAABR(1109) = EE**2/SW**2*Sqrt2*Zh11*Zl21/cb-EE**2
     &      /SW**2*Sqrt2*Zh12*Zl21*sb-EE**2/SW**2*Sqrt2*Zh11
     &      *Zl21/cb*sb**2-EE**2/MW**2*Ml**2/SW**2*Sqrt2
     &      *Zh11*Zl21/cb+2D0*EE*Lambda/MW*Ml/SW*Zh13*Zl22
     &      /cb*sb
      AAABR(1110) = EE**2/SW**2*Sqrt2*Zh31*Zl21/cb-EE**2
     &      /SW**2*Sqrt2*Zh32*Zl21*sb-EE**2/SW**2*Sqrt2*Zh31
     &      *Zl21/cb*sb**2-EE**2/MW**2*Ml**2/SW**2*Sqrt2
     &      *Zh31*Zl21/cb+2D0*EE*Lambda/MW*Ml/SW*Zh33*Zl22
     &      /cb*sb
      AAABR(1111) = EE**2*(Zl21**2/CW**2-2D0*Zl22**2/CW**2
     &      -Zl21**2/SW**2-2D0*Zl21**2/CW**2*sb**2+4D0*Zl22
     &      **2/CW**2*sb**2+2D0*Ml**2/MW**2/SW**2+2D0*Zl21
     &      **2/SW**2*sb**2)
      AAABR(1112) = EE**2/CW**2*Pa11*cb*sb
      AAABR(1113) = tB*(EE**2/CW**2*Pa11*Zl21**2-2D0*EE**2
     &      /CW**2*Pa11*Zl22**2-EE**2*Pa11/SW**2*Zl21**2
     &      -EE**2/CW**2*Pa11*Zl21**2*sb**2+2D0*EE**2/CW
     &      **2*Pa11*Zl22**2*sb**2+EE**2/MW**2*Ml**2*Pa11
     &      /SW**2+EE**2*Pa11/SW**2*Zl21**2*sb**2-EE*Lambda
     &      /MW*Ml*Pa12/SW*Sqrt2*Zl21*Zl22)
      AAABR(1114) = EE**2/CW**2*Pa21*cb*sb
      AAABR(1115) = tB*(EE**2/CW**2*Pa21*Zl21**2-2D0*EE**2
     &      /CW**2*Pa21*Zl22**2-EE**2*Pa21/SW**2*Zl21**2
     &      -EE**2/CW**2*Pa21*Zl21**2*sb**2+2D0*EE**2/CW
     &      **2*Pa21*Zl22**2*sb**2+EE**2/MW**2*Ml**2*Pa21
     &      /SW**2+EE**2*Pa21/SW**2*Zl21**2*sb**2-EE*Lambda
     &      /MW*Ml*Pa22/SW*Sqrt2*Zl21*Zl22)
      AAABR(1116) = EE**2/CW**2*Pa21**2*(1D0-2D0*sb**2)
      AAABR(1117) = 1D0/cb**2*(EE**2/CW**2*Pa11*Pa21*Zl21
     &      **2-2D0*EE**2/CW**2*Pa11*Pa21*Zl22**2-EE**2*Pa11
     &      *Pa21/SW**2*Zl21**2-3D0*EE**2/CW**2*Pa11*Pa21
     &      *Zl21**2*sb**2+2D0*EE**2/CW**2*Pa11*Pa21*Zl21
     &      **2*sb**4+6D0*EE**2/CW**2*Pa11*Pa21*Zl22**2*sb
     &      **2-4D0*EE**2/CW**2*Pa11*Pa21*Zl22**2*sb**4+3D0
     &      *EE**2*Pa11*Pa21/SW**2*Zl21**2*sb**2-2D0*EE**2
     &      *Pa11*Pa21/SW**2*Zl21**2*sb**4-2D0*EE**2/MW**2
     &      *Ml**2*Pa11*Pa21/SW**2*sb**2-2D0*EE*Lambda/MW
     &      *Ml*Pa11*Pa22/SW*Sqrt2*Zl21*Zl22-2D0*EE*Lambda
     &      /MW*Ml*Pa12*Pa21/SW*Sqrt2*Zl21*Zl22+2D0*EE*Lambda
     &      /MW*Ml*Pa11*Pa22/SW*Sqrt2*Zl21*Zl22*sb**2+2D0
     &      *EE*Lambda/MW*Ml*Pa12*Pa21/SW*Sqrt2*Zl21*Zl22
     &      *sb**2)
      AAABR(1118) = 1D0/cb**2*(EE**2/CW**2*Pa21**2*Zl21**2
     &      -2D0*EE**2/CW**2*Pa21**2*Zl22**2-EE**2*Pa21**2
     &      /SW**2*Zl21**2-3D0*EE**2/CW**2*Pa21**2*Zl21**2
     &      *sb**2+2D0*EE**2/CW**2*Pa21**2*Zl21**2*sb**4
     &      +6D0*EE**2/CW**2*Pa21**2*Zl22**2*sb**2-4D0*EE
     &      **2/CW**2*Pa21**2*Zl22**2*sb**4+3D0*EE**2*Pa21
     &      **2/SW**2*Zl21**2*sb**2-2D0*EE**2*Pa21**2/SW
     &      **2*Zl21**2*sb**4-2D0*EE**2/MW**2*Ml**2*Pa21
     &      **2/SW**2*sb**2-4D0*EE*Lambda/MW*Ml*Pa21*Pa22
     &      /SW*Sqrt2*Zl21*Zl22+4D0*EE*Lambda/MW*Ml*Pa21
     &      *Pa22/SW*Sqrt2*Zl21*Zl22*sb**2)
      AAABR(1119) = EE**2/CW**2*Pa11**2*(1D0-2D0*sb**2)
      AAABR(1120) = 1D0/cb**2*(EE**2/CW**2*Pa11**2*Zl21**2
     &      -2D0*EE**2/CW**2*Pa11**2*Zl22**2-EE**2*Pa11**2
     &      /SW**2*Zl21**2-3D0*EE**2/CW**2*Pa11**2*Zl21**2
     &      *sb**2+2D0*EE**2/CW**2*Pa11**2*Zl21**2*sb**4
     &      +6D0*EE**2/CW**2*Pa11**2*Zl22**2*sb**2-4D0*EE
     &      **2/CW**2*Pa11**2*Zl22**2*sb**4+3D0*EE**2*Pa11
     &      **2/SW**2*Zl21**2*sb**2-2D0*EE**2*Pa11**2/SW
     &      **2*Zl21**2*sb**4-2D0*EE**2/MW**2*Ml**2*Pa11
     &      **2/SW**2*sb**2-4D0*EE*Lambda/MW*Ml*Pa11*Pa12
     &      /SW*Sqrt2*Zl21*Zl22+4D0*EE*Lambda/MW*Ml*Pa11
     &      *Pa12/SW*Sqrt2*Zl21*Zl22*sb**2)
      AAABR(1121) = EE**2/CW**2*Pa11*Pa21*(1D0-2D0*sb**2)
      AAABR(1122) = EE**2/CW**2*Zh21*Zh31*Zl21**2/cb**2-2D0
     &      *EE**2/CW**2*Zh21*Zh31*Zl22**2/cb**2-EE**2/CW
     &      **2*Zh22*Zh32*Zl21**2/cb**2+2D0*EE**2/CW**2*Zh22
     &      *Zh32*Zl22**2/cb**2-EE**2/SW**2*Zh21*Zh31*Zl21
     &      **2/cb**2+EE**2/SW**2*Zh22*Zh32*Zl21**2/cb**2
     &      -EE**2/CW**2*Zh21*Zh31*Zl21**2/cb**2*sb**2+2D0
     &      *EE**2/CW**2*Zh21*Zh31*Zl22**2/cb**2*sb**2+EE
     &      **2/CW**2*Zh22*Zh32*Zl21**2/cb**2*sb**2-2D0*EE
     &      **2/CW**2*Zh22*Zh32*Zl22**2/cb**2*sb**2+2D0*EE
     &      **2/MW**2*Ml**2/SW**2*Zh21*Zh31/cb**2+EE**2/SW
     &      **2*Zh21*Zh31*Zl21**2/cb**2*sb**2-EE**2/SW**2
     &      *Zh22*Zh32*Zl21**2/cb**2*sb**2-2D0*EE*Lambda
     &      /MW*Ml/SW*Sqrt2*Zh22*Zh33*Zl21*Zl22/cb-2D0*EE
     &      *Lambda/MW*Ml/SW*Sqrt2*Zh23*Zh32*Zl21*Zl22/cb
      AAABR(1123) = EE**2/CW**2*Zh31**2*Zl21**2/cb**2-2D0
     &      *EE**2/CW**2*Zh31**2*Zl22**2/cb**2-EE**2/CW**2
     &      *Zh32**2*Zl21**2/cb**2+2D0*EE**2/CW**2*Zh32**2
     &      *Zl22**2/cb**2-EE**2/SW**2*Zh31**2*Zl21**2/cb
     &      **2+EE**2/SW**2*Zh32**2*Zl21**2/cb**2-EE**2/CW
     &      **2*Zh31**2*Zl21**2/cb**2*sb**2+2D0*EE**2/CW
     &      **2*Zh31**2*Zl22**2/cb**2*sb**2+EE**2/CW**2*Zh32
     &      **2*Zl21**2/cb**2*sb**2-2D0*EE**2/CW**2*Zh32
     &      **2*Zl22**2/cb**2*sb**2+2D0*EE**2/MW**2*Ml**2
     &      /SW**2*Zh31**2/cb**2+EE**2/SW**2*Zh31**2*Zl21
     &      **2/cb**2*sb**2-EE**2/SW**2*Zh32**2*Zl21**2/cb
     &      **2*sb**2-4D0*EE*Lambda/MW*Ml/SW*Sqrt2*Zh32*Zh33
     &      *Zl21*Zl22/cb
      AAABR(1124) = EE**2/CW**2*(Zh11*Zh21-Zh12*Zh22)
      AAABR(1125) = EE**2/CW**2*(Zh11*Zh31-Zh12*Zh32)
      AAABR(1126) = EE**2/CW**2*(Zh21**2-Zh22**2)
      AAABR(1127) = EE**2/CW**2*(Zh11**2-Zh12**2)
      AAABR(1128) = EE**2/CW**2*(Zh21*Zh31-Zh22*Zh32)
      AAABR(1129) = EE**2/CW**2*(Zh31**2-Zh32**2)
      AAABR(1130) = EE**2/CW**2*Zh11**2*Zl21**2/cb**2-2D0
     &      *EE**2/CW**2*Zh11**2*Zl22**2/cb**2-EE**2/CW**2
     &      *Zh12**2*Zl21**2/cb**2+2D0*EE**2/CW**2*Zh12**2
     &      *Zl22**2/cb**2-EE**2/SW**2*Zh11**2*Zl21**2/cb
     &      **2+EE**2/SW**2*Zh12**2*Zl21**2/cb**2-EE**2/CW
     &      **2*Zh11**2*Zl21**2/cb**2*sb**2+2D0*EE**2/CW
     &      **2*Zh11**2*Zl22**2/cb**2*sb**2+EE**2/CW**2*Zh12
     &      **2*Zl21**2/cb**2*sb**2-2D0*EE**2/CW**2*Zh12
     &      **2*Zl22**2/cb**2*sb**2+2D0*EE**2/MW**2*Ml**2
     &      /SW**2*Zh11**2/cb**2+EE**2/SW**2*Zh11**2*Zl21
     &      **2/cb**2*sb**2-EE**2/SW**2*Zh12**2*Zl21**2/cb
     &      **2*sb**2-4D0*EE*Lambda/MW*Ml/SW*Sqrt2*Zh12*Zh13
     &      *Zl21*Zl22/cb
      AAABR(1131) = EE**2/CW**2*Zh11*Zh21*Zl21**2/cb**2-2D0
     &      *EE**2/CW**2*Zh11*Zh21*Zl22**2/cb**2-EE**2/CW
     &      **2*Zh12*Zh22*Zl21**2/cb**2+2D0*EE**2/CW**2*Zh12
     &      *Zh22*Zl22**2/cb**2-EE**2/SW**2*Zh11*Zh21*Zl21
     &      **2/cb**2+EE**2/SW**2*Zh12*Zh22*Zl21**2/cb**2
     &      -EE**2/CW**2*Zh11*Zh21*Zl21**2/cb**2*sb**2+2D0
     &      *EE**2/CW**2*Zh11*Zh21*Zl22**2/cb**2*sb**2+EE
     &      **2/CW**2*Zh12*Zh22*Zl21**2/cb**2*sb**2-2D0*EE
     &      **2/CW**2*Zh12*Zh22*Zl22**2/cb**2*sb**2+2D0*EE
     &      **2/MW**2*Ml**2/SW**2*Zh11*Zh21/cb**2+EE**2/SW
     &      **2*Zh11*Zh21*Zl21**2/cb**2*sb**2-EE**2/SW**2
     &      *Zh12*Zh22*Zl21**2/cb**2*sb**2-2D0*EE*Lambda
     &      /MW*Ml/SW*Sqrt2*Zh12*Zh23*Zl21*Zl22/cb-2D0*EE
     &      *Lambda/MW*Ml/SW*Sqrt2*Zh13*Zh22*Zl21*Zl22/cb
      AAABR(1132) = EE**2/CW**2*Zh21**2*Zl21**2/cb**2-2D0
     &      *EE**2/CW**2*Zh21**2*Zl22**2/cb**2-EE**2/CW**2
     &      *Zh22**2*Zl21**2/cb**2+2D0*EE**2/CW**2*Zh22**2
     &      *Zl22**2/cb**2-EE**2/SW**2*Zh21**2*Zl21**2/cb
     &      **2+EE**2/SW**2*Zh22**2*Zl21**2/cb**2-EE**2/CW
     &      **2*Zh21**2*Zl21**2/cb**2*sb**2+2D0*EE**2/CW
     &      **2*Zh21**2*Zl22**2/cb**2*sb**2+EE**2/CW**2*Zh22
     &      **2*Zl21**2/cb**2*sb**2-2D0*EE**2/CW**2*Zh22
     &      **2*Zl22**2/cb**2*sb**2+2D0*EE**2/MW**2*Ml**2
     &      /SW**2*Zh21**2/cb**2+EE**2/SW**2*Zh21**2*Zl21
     &      **2/cb**2*sb**2-EE**2/SW**2*Zh22**2*Zl21**2/cb
     &      **2*sb**2-4D0*EE*Lambda/MW*Ml/SW*Sqrt2*Zh22*Zh23
     &      *Zl21*Zl22/cb
      AAABR(1133) = EE**2/CW**2*Zh11*Zh31*Zl21**2/cb**2-2D0
     &      *EE**2/CW**2*Zh11*Zh31*Zl22**2/cb**2-EE**2/CW
     &      **2*Zh12*Zh32*Zl21**2/cb**2+2D0*EE**2/CW**2*Zh12
     &      *Zh32*Zl22**2/cb**2-EE**2/SW**2*Zh11*Zh31*Zl21
     &      **2/cb**2+EE**2/SW**2*Zh12*Zh32*Zl21**2/cb**2
     &      -EE**2/CW**2*Zh11*Zh31*Zl21**2/cb**2*sb**2+2D0
     &      *EE**2/CW**2*Zh11*Zh31*Zl22**2/cb**2*sb**2+EE
     &      **2/CW**2*Zh12*Zh32*Zl21**2/cb**2*sb**2-2D0*EE
     &      **2/CW**2*Zh12*Zh32*Zl22**2/cb**2*sb**2+2D0*EE
     &      **2/MW**2*Ml**2/SW**2*Zh11*Zh31/cb**2+EE**2/SW
     &      **2*Zh11*Zh31*Zl21**2/cb**2*sb**2-EE**2/SW**2
     &      *Zh12*Zh32*Zl21**2/cb**2*sb**2-2D0*EE*Lambda
     &      /MW*Ml/SW*Sqrt2*Zh12*Zh33*Zl21*Zl22/cb-2D0*EE
     &      *Lambda/MW*Ml/SW*Sqrt2*Zh13*Zh32*Zl21*Zl22/cb
      AAABR(1134) = EE**2*(Zss11**2/CW**2+2D0*Zss12**2/CW
     &      **2+3D0*Zss11**2/SW**2)
      AAABR(1135) = EE**2*(Zdd11**2/CW**2+2D0*Zdd12**2/CW
     &      **2+3D0*Zdd11**2/SW**2)
      AAABR(1136) = EE**2*(Zb11**2/CW**2+2D0*Zb12**2/CW**2
     &      +3D0*Zb11**2/SW**2)
      AAABR(1137) = EE**2*(Zb11/CW**2*Zb21+2D0*Zb12/CW**2
     &      *Zb22+3D0*Zb11/SW**2*Zb21)
      AAABR(1138) = EE**2*(Zdd11/CW**2*Zdd21+2D0*Zdd12/CW
     &      **2*Zdd22+3D0*Zdd11/SW**2*Zdd21)
      AAABR(1139) = EE**2*(Zss11/CW**2*Zss21+2D0*Zss12/CW
     &      **2*Zss22+3D0*Zss11/SW**2*Zss21)
      AAABR(1140) = EE**2*(Zss21**2/CW**2+2D0*Zss22**2/CW
     &      **2+3D0*Zss21**2/SW**2)
      AAABR(1141) = EE**2*(Zb21**2/CW**2+2D0*Zb22**2/CW**2
     &      +3D0*Zb21**2/SW**2)
      AAABR(1142) = EE**2*(Zdd21**2/CW**2+2D0*Zdd22**2/CW
     &      **2+3D0*Zdd21**2/SW**2)
      AAABR(1143) = EE**2*(Zt11**2/CW**2-4D0*Zt12**2/CW**2
     &      -3D0*Zt11**2/SW**2)
      AAABR(1144) = EE**2*(Zcc11**2/CW**2-4D0*Zcc12**2/CW
     &      **2-3D0*Zcc11**2/SW**2)
      AAABR(1145) = EE**2*(Zuu11**2/CW**2-4D0*Zuu12**2/CW
     &      **2-3D0*Zuu11**2/SW**2)
      AAABR(1146) = EE**2*(Zt11/CW**2*Zt21-4D0*Zt12/CW**2
     &      *Zt22-3D0*Zt11/SW**2*Zt21)
      AAABR(1147) = EE**2*(Zcc11/CW**2*Zcc21-4D0*Zcc12/CW
     &      **2*Zcc22-3D0*Zcc11/SW**2*Zcc21)
      AAABR(1148) = EE**2*(Zuu11/CW**2*Zuu21-4D0*Zuu12/CW
     &      **2*Zuu22-3D0*Zuu11/SW**2*Zuu21)
      AAABR(1149) = EE**2*(Zt21**2/CW**2-4D0*Zt22**2/CW**2
     &      -3D0*Zt21**2/SW**2)
      AAABR(1150) = EE**2*(Zcc21**2/CW**2-4D0*Zcc22**2/CW
     &      **2-3D0*Zcc21**2/SW**2)
      AAABR(1151) = EE**2*(Zuu21**2/CW**2-4D0*Zuu22**2/CW
     &      **2-3D0*Zuu21**2/SW**2)
      AAABR(1152) = EE**2*(1D0/CW**2-1D0/SW**2-2D0*sb**2
     &      /CW**2+2D0*sb**2/SW**2+2D0*Ml**2/MW**2/SW**2)
      AAABR(1153) = EE**2*Pa11*cb*sb*(1D0/CW**2+1D0/SW**2)
      AAABR(1154) = EE**2*Pa21*cb*sb*(1D0/CW**2+1D0/SW**2)
      AAABR(1155) = EE**2*Pa11*Pa21*(1D0/CW**2+1D0/SW**2
     &      -2D0*sb**2/CW**2-2D0*sb**2/SW**2)
      AAABR(1156) = EE**2*Pa21**2*(1D0/CW**2+1D0/SW**2-2D0
     &      *sb**2/CW**2-2D0*sb**2/SW**2)
      AAABR(1157) = EE**2*Pa11**2*(1D0/CW**2+1D0/SW**2-2D0
     &      *sb**2/CW**2-2D0*sb**2/SW**2)
      AAABR(1158) = EE**2*(Zh21**2/CW**2-Zh22**2/CW**2+Zh21
     &      **2/SW**2-Zh22**2/SW**2)
      AAABR(1159) = EE**2*(Zh11/CW**2*Zh21-Zh12/CW**2*Zh22
     &      +Zh11/SW**2*Zh21-Zh12/SW**2*Zh22)
      AAABR(1160) = EE**2*(Zh11**2/CW**2-Zh12**2/CW**2+Zh11
     &      **2/SW**2-Zh12**2/SW**2)
      AAABR(1161) = EE**2*(Zh21/CW**2*Zh31-Zh22/CW**2*Zh32
     &      +Zh21/SW**2*Zh31-Zh22/SW**2*Zh32)
      AAABR(1162) = EE**2*(Zh31**2/CW**2-Zh32**2/CW**2+Zh31
     &      **2/SW**2-Zh32**2/SW**2)
      AAABR(1163) = EE**2*(Zh11/CW**2*Zh31-Zh12/CW**2*Zh32
     &      +Zh11/SW**2*Zh31-Zh12/SW**2*Zh32)
      AAABR(1164) = EE**2/cb**2*(Zb11**2/CW**2*Zdd11**2+2D0
     &      *Zb11**2/CW**2*Zdd12**2+2D0*Zb12**2/CW**2*Zdd11
     &      **2+4D0*Zb12**2/CW**2*Zdd12**2+9D0*Zb11**2/SW
     &      **2*Zdd11**2-Zb11**2/CW**2*Zdd11**2*sb**2-2D0
     &      *Zb11**2/CW**2*Zdd12**2*sb**2-2D0*Zb12**2/CW
     &      **2*Zdd11**2*sb**2-4D0*Zb12**2/CW**2*Zdd12**2
     &      *sb**2-9D0*Zb11**2/SW**2*Zdd11**2*sb**2+36D0
     &      *MbMM/MW**2*Mq/SW**2*Zb11*Zb12*Zdd11*Zdd12)
      AAABR(1165) = GG**2/cb**2*(1D0-sb**2)
      AAABR(1166) = GG**2/cb**2*(Zb11**4+Zb12**4+2D0*Zb11
     &      **2*Zb12**2-Zb11**4*sb**2-Zb12**4*sb**2-2D0*Zb11
     &      **2*Zb12**2*sb**2)
      AAABR(1167) = EE**2/cb**2*(Zb11**4/CW**2+4D0*Zb12**4
     &      /CW**2+9D0*Zb11**4/SW**2+4D0*Zb11**2/CW**2*Zb12
     &      **2-Zb11**4/CW**2*sb**2-4D0*Zb12**4/CW**2*sb
     &      **2-9D0*Zb11**4/SW**2*sb**2-4D0*Zb11**2/CW**2
     &      *Zb12**2*sb**2+36D0*MbMM**2/MW**2/SW**2*Zb11
     &      **2*Zb12**2)
      AAABR(1168) = GG**2/cb**2*(Zss11**4+Zss12**4+2D0*Zss11
     &      **2*Zss12**2-Zss11**4*sb**2-Zss12**4*sb**2-2D0
     &      *Zss11**2*Zss12**2*sb**2)
      AAABR(1169) = EE**2/cb**2*(Zss11**4/CW**2+4D0*Zss12
     &      **4/CW**2+9D0*Zss11**4/SW**2+4D0*Zss11**2/CW
     &      **2*Zss12**2-Zss11**4/CW**2*sb**2-4D0*Zss12**4
     &      /CW**2*sb**2-9D0*Zss11**4/SW**2*sb**2-4D0*Zss11
     &      **2/CW**2*Zss12**2*sb**2+36D0*Mq**2/MW**2/SW
     &      **2*Zss11**2*Zss12**2)
      AAABR(1170) = EE**2/cb**2*(Zb11**2/CW**2*Zss11**2+2D0
     &      *Zb11**2/CW**2*Zss12**2+2D0*Zb12**2/CW**2*Zss11
     &      **2+4D0*Zb12**2/CW**2*Zss12**2+9D0*Zb11**2/SW
     &      **2*Zss11**2-Zb11**2/CW**2*Zss11**2*sb**2-2D0
     &      *Zb11**2/CW**2*Zss12**2*sb**2-2D0*Zb12**2/CW
     &      **2*Zss11**2*sb**2-4D0*Zb12**2/CW**2*Zss12**2
     &      *sb**2-9D0*Zb11**2/SW**2*Zss11**2*sb**2+36D0
     &      *MbMM/MW**2*Mq/SW**2*Zb11*Zb12*Zss11*Zss12)
      AAABR(1171) = EE**2/cb**2*(Zdd11**2/CW**2*Zss11**2
     &      +2D0*Zdd11**2/CW**2*Zss12**2+2D0*Zdd12**2/CW
     &      **2*Zss11**2+4D0*Zdd12**2/CW**2*Zss12**2+9D0
     &      *Zdd11**2/SW**2*Zss11**2-Zdd11**2/CW**2*Zss11
     &      **2*sb**2-2D0*Zdd11**2/CW**2*Zss12**2*sb**2-2D0
     &      *Zdd12**2/CW**2*Zss11**2*sb**2-4D0*Zdd12**2/CW
     &      **2*Zss12**2*sb**2-9D0*Zdd11**2/SW**2*Zss11**2
     &      *sb**2+36D0*Mq**2/MW**2/SW**2*Zdd11*Zdd12*Zss11
     &      *Zss12)
      AAABR(1172) = GG**2/cb**2*(Zdd11**4+Zdd12**4+2D0*Zdd11
     &      **2*Zdd12**2-Zdd11**4*sb**2-Zdd12**4*sb**2-2D0
     &      *Zdd11**2*Zdd12**2*sb**2)
      AAABR(1173) = EE**2/cb**2*(Zdd11**4/CW**2+4D0*Zdd12
     &      **4/CW**2+9D0*Zdd11**4/SW**2+4D0*Zdd11**2/CW
     &      **2*Zdd12**2-Zdd11**4/CW**2*sb**2-4D0*Zdd12**4
     &      /CW**2*sb**2-9D0*Zdd11**4/SW**2*sb**2-4D0*Zdd11
     &      **2/CW**2*Zdd12**2*sb**2+36D0*Mq**2/MW**2/SW
     &      **2*Zdd11**2*Zdd12**2)
      AAABR(1174) = EE**2/cb**2*(Zb11/CW**2*Zb21*Zdd11**2
     &      +2D0*Zb11/CW**2*Zb21*Zdd12**2+2D0*Zb12/CW**2
     &      *Zb22*Zdd11**2+4D0*Zb12/CW**2*Zb22*Zdd12**2+9D0
     &      *Zb11/SW**2*Zb21*Zdd11**2-Zb11/CW**2*Zb21*Zdd11
     &      **2*sb**2-2D0*Zb11/CW**2*Zb21*Zdd12**2*sb**2
     &      -2D0*Zb12/CW**2*Zb22*Zdd11**2*sb**2-4D0*Zb12
     &      /CW**2*Zb22*Zdd12**2*sb**2-9D0*Zb11/SW**2*Zb21
     &      *Zdd11**2*sb**2+18D0*MbMM/MW**2*Mq/SW**2*Zb11
     &      *Zb22*Zdd11*Zdd12+18D0*MbMM/MW**2*Mq/SW**2*Zb12
     &      *Zb21*Zdd11*Zdd12)
      AAABR(1175) = GG**2/cb**2*(Zb11**3*Zb21+Zb12**3*Zb22
     &      +Zb11*Zb12**2*Zb21+Zb11**2*Zb12*Zb22-Zb11**3
     &      *Zb21*sb**2-Zb12**3*Zb22*sb**2-Zb11*Zb12**2*Zb21
     &      *sb**2-Zb11**2*Zb12*Zb22*sb**2)
      AAABR(1176) = EE**2/cb**2*(Zb11**3/CW**2*Zb21+4D0*Zb12
     &      **3/CW**2*Zb22+9D0*Zb11**3/SW**2*Zb21+2D0*Zb11
     &      /CW**2*Zb12**2*Zb21+2D0*Zb11**2/CW**2*Zb12*Zb22
     &      -Zb11**3/CW**2*Zb21*sb**2-4D0*Zb12**3/CW**2*Zb22
     &      *sb**2-9D0*Zb11**3/SW**2*Zb21*sb**2-2D0*Zb11
     &      /CW**2*Zb12**2*Zb21*sb**2-2D0*Zb11**2/CW**2*Zb12
     &      *Zb22*sb**2+18D0*MbMM**2/MW**2/SW**2*Zb11*Zb12
     &      **2*Zb21+18D0*MbMM**2/MW**2/SW**2*Zb11**2*Zb12
     &      *Zb22)
      AAABR(1177) = EE**2/cb**2*(Zdd11/CW**2*Zdd21*Zss11
     &      **2+2D0*Zdd11/CW**2*Zdd21*Zss12**2+2D0*Zdd12
     &      /CW**2*Zdd22*Zss11**2+4D0*Zdd12/CW**2*Zdd22*Zss12
     &      **2+9D0*Zdd11/SW**2*Zdd21*Zss11**2-Zdd11/CW**2
     &      *Zdd21*Zss11**2*sb**2-2D0*Zdd11/CW**2*Zdd21*Zss12
     &      **2*sb**2-2D0*Zdd12/CW**2*Zdd22*Zss11**2*sb**2
     &      -4D0*Zdd12/CW**2*Zdd22*Zss12**2*sb**2-9D0*Zdd11
     &      /SW**2*Zdd21*Zss11**2*sb**2+18D0*Mq**2/MW**2
     &      /SW**2*Zdd11*Zdd22*Zss11*Zss12+18D0*Mq**2/MW
     &      **2/SW**2*Zdd12*Zdd21*Zss11*Zss12)
      AAABR(1178) = GG**2/cb**2*(Zss11**3*Zss21+Zss12**3
     &      *Zss22+Zss11*Zss12**2*Zss21+Zss11**2*Zss12*Zss22
     &      -Zss11**3*Zss21*sb**2-Zss12**3*Zss22*sb**2-Zss11
     &      *Zss12**2*Zss21*sb**2-Zss11**2*Zss12*Zss22*sb
     &      **2)
      AAABR(1179) = EE**2/cb**2*(Zss11**3/CW**2*Zss21+4D0
     &      *Zss12**3/CW**2*Zss22+9D0*Zss11**3/SW**2*Zss21
     &      +2D0*Zss11/CW**2*Zss12**2*Zss21+2D0*Zss11**2
     &      /CW**2*Zss12*Zss22-Zss11**3/CW**2*Zss21*sb**2
     &      -4D0*Zss12**3/CW**2*Zss22*sb**2-9D0*Zss11**3
     &      /SW**2*Zss21*sb**2-2D0*Zss11/CW**2*Zss12**2*Zss21
     &      *sb**2-2D0*Zss11**2/CW**2*Zss12*Zss22*sb**2+18D0
     &      *Mq**2/MW**2/SW**2*Zss11*Zss12**2*Zss21+18D0
     &      *Mq**2/MW**2/SW**2*Zss11**2*Zss12*Zss22)
      AAABR(1180) = EE**2/cb**2*(Zb11/CW**2*Zb21*Zss11**2
     &      +2D0*Zb11/CW**2*Zb21*Zss12**2+2D0*Zb12/CW**2
     &      *Zb22*Zss11**2+4D0*Zb12/CW**2*Zb22*Zss12**2+9D0
     &      *Zb11/SW**2*Zb21*Zss11**2-Zb11/CW**2*Zb21*Zss11
     &      **2*sb**2-2D0*Zb11/CW**2*Zb21*Zss12**2*sb**2
     &      -2D0*Zb12/CW**2*Zb22*Zss11**2*sb**2-4D0*Zb12
     &      /CW**2*Zb22*Zss12**2*sb**2-9D0*Zb11/SW**2*Zb21
     &      *Zss11**2*sb**2+18D0*MbMM/MW**2*Mq/SW**2*Zb11
     &      *Zb22*Zss11*Zss12+18D0*MbMM/MW**2*Mq/SW**2*Zb12
     &      *Zb21*Zss11*Zss12)
      AAABR(1181) = EE**2/cb**2*(Zb11**2/CW**2*Zdd11*Zdd21
     &      +2D0*Zb11**2/CW**2*Zdd12*Zdd22+2D0*Zb12**2/CW
     &      **2*Zdd11*Zdd21+4D0*Zb12**2/CW**2*Zdd12*Zdd22
     &      +9D0*Zb11**2/SW**2*Zdd11*Zdd21-Zb11**2/CW**2
     &      *Zdd11*Zdd21*sb**2-2D0*Zb11**2/CW**2*Zdd12*Zdd22
     &      *sb**2-2D0*Zb12**2/CW**2*Zdd11*Zdd21*sb**2-4D0
     &      *Zb12**2/CW**2*Zdd12*Zdd22*sb**2-9D0*Zb11**2
     &      /SW**2*Zdd11*Zdd21*sb**2+18D0*MbMM/MW**2*Mq/SW
     &      **2*Zb11*Zb12*Zdd11*Zdd22+18D0*MbMM/MW**2*Mq
     &      /SW**2*Zb11*Zb12*Zdd12*Zdd21)
      AAABR(1182) = EE**2/cb**2*(Zdd11**2/CW**2*Zss11*Zss21
     &      +2D0*Zdd11**2/CW**2*Zss12*Zss22+2D0*Zdd12**2
     &      /CW**2*Zss11*Zss21+4D0*Zdd12**2/CW**2*Zss12*Zss22
     &      +9D0*Zdd11**2/SW**2*Zss11*Zss21-Zdd11**2/CW**2
     &      *Zss11*Zss21*sb**2-2D0*Zdd11**2/CW**2*Zss12*Zss22
     &      *sb**2-2D0*Zdd12**2/CW**2*Zss11*Zss21*sb**2-4D0
     &      *Zdd12**2/CW**2*Zss12*Zss22*sb**2-9D0*Zdd11**2
     &      /SW**2*Zss11*Zss21*sb**2+18D0*Mq**2/MW**2/SW
     &      **2*Zdd11*Zdd12*Zss11*Zss22+18D0*Mq**2/MW**2
     &      /SW**2*Zdd11*Zdd12*Zss12*Zss21)
      AAABR(1183) = EE**2/cb**2*(Zb11**2/CW**2*Zss11*Zss21
     &      +2D0*Zb11**2/CW**2*Zss12*Zss22+2D0*Zb12**2/CW
     &      **2*Zss11*Zss21+4D0*Zb12**2/CW**2*Zss12*Zss22
     &      +9D0*Zb11**2/SW**2*Zss11*Zss21-Zb11**2/CW**2
     &      *Zss11*Zss21*sb**2-2D0*Zb11**2/CW**2*Zss12*Zss22
     &      *sb**2-2D0*Zb12**2/CW**2*Zss11*Zss21*sb**2-4D0
     &      *Zb12**2/CW**2*Zss12*Zss22*sb**2-9D0*Zb11**2
     &      /SW**2*Zss11*Zss21*sb**2+18D0*MbMM/MW**2*Mq/SW
     &      **2*Zb11*Zb12*Zss11*Zss22+18D0*MbMM/MW**2*Mq
     &      /SW**2*Zb11*Zb12*Zss12*Zss21)
      AAABR(1184) = GG**2/cb**2*(Zdd11**3*Zdd21+Zdd12**3
     &      *Zdd22+Zdd11*Zdd12**2*Zdd21+Zdd11**2*Zdd12*Zdd22
     &      -Zdd11**3*Zdd21*sb**2-Zdd12**3*Zdd22*sb**2-Zdd11
     &      *Zdd12**2*Zdd21*sb**2-Zdd11**2*Zdd12*Zdd22*sb
     &      **2)
      AAABR(1185) = EE**2/cb**2*(Zdd11**3/CW**2*Zdd21+4D0
     &      *Zdd12**3/CW**2*Zdd22+9D0*Zdd11**3/SW**2*Zdd21
     &      +2D0*Zdd11/CW**2*Zdd12**2*Zdd21+2D0*Zdd11**2
     &      /CW**2*Zdd12*Zdd22-Zdd11**3/CW**2*Zdd21*sb**2
     &      -4D0*Zdd12**3/CW**2*Zdd22*sb**2-9D0*Zdd11**3
     &      /SW**2*Zdd21*sb**2-2D0*Zdd11/CW**2*Zdd12**2*Zdd21
     &      *sb**2-2D0*Zdd11**2/CW**2*Zdd12*Zdd22*sb**2+18D0
     &      *Mq**2/MW**2/SW**2*Zdd11*Zdd12**2*Zdd21+18D0
     &      *Mq**2/MW**2/SW**2*Zdd11**2*Zdd12*Zdd22)
      AAABR(1186) = GG**2/cb**2*(Zb11**2*Zb21**2+Zb12**2
     &      *Zb22**2-Zb11**2*Zb21**2*sb**2-Zb12**2*Zb22**2
     &      *sb**2+2D0*Zb11*Zb12*Zb21*Zb22-2D0*Zb11*Zb12
     &      *Zb21*Zb22*sb**2)
      AAABR(1187) = EE**2/cb**2*(Zb11**2/CW**2*Zb21**2+4D0
     &      *Zb12**2/CW**2*Zb22**2+9D0*Zb11**2/SW**2*Zb21
     &      **2-Zb11**2/CW**2*Zb21**2*sb**2-4D0*Zb12**2/CW
     &      **2*Zb22**2*sb**2-9D0*Zb11**2/SW**2*Zb21**2*sb
     &      **2+4D0*Zb11/CW**2*Zb12*Zb21*Zb22-4D0*Zb11/CW
     &      **2*Zb12*Zb21*Zb22*sb**2+36D0*MbMM**2/MW**2/SW
     &      **2*Zb11*Zb12*Zb21*Zb22)
      AAABR(1188) = EE**2/cb**2*(Zb11/CW**2*Zb21*Zss11*Zss21
     &      +2D0*Zb11/CW**2*Zb21*Zss12*Zss22+2D0*Zb12/CW
     &      **2*Zb22*Zss11*Zss21+4D0*Zb12/CW**2*Zb22*Zss12
     &      *Zss22+9D0*Zb11/SW**2*Zb21*Zss11*Zss21-Zb11/CW
     &      **2*Zb21*Zss11*Zss21*sb**2-2D0*Zb11/CW**2*Zb21
     &      *Zss12*Zss22*sb**2-2D0*Zb12/CW**2*Zb22*Zss11
     &      *Zss21*sb**2-4D0*Zb12/CW**2*Zb22*Zss12*Zss22
     &      *sb**2-9D0*Zb11/SW**2*Zb21*Zss11*Zss21*sb**2
     &      +18D0*MbMM/MW**2*Mq/SW**2*Zb11*Zb22*Zss12*Zss21
     &      +18D0*MbMM/MW**2*Mq/SW**2*Zb12*Zb21*Zss11*Zss22)
      AAABR(1189) = GG**2/cb**2*(Zdd11**2*Zdd21**2+Zdd12
     &      **2*Zdd22**2-Zdd11**2*Zdd21**2*sb**2-Zdd12**2
     &      *Zdd22**2*sb**2+2D0*Zdd11*Zdd12*Zdd21*Zdd22-2D0
     &      *Zdd11*Zdd12*Zdd21*Zdd22*sb**2)
      AAABR(1190) = EE**2/cb**2*(Zdd11**2/CW**2*Zdd21**2
     &      +4D0*Zdd12**2/CW**2*Zdd22**2+9D0*Zdd11**2/SW
     &      **2*Zdd21**2-Zdd11**2/CW**2*Zdd21**2*sb**2-4D0
     &      *Zdd12**2/CW**2*Zdd22**2*sb**2-9D0*Zdd11**2/SW
     &      **2*Zdd21**2*sb**2+4D0*Zdd11/CW**2*Zdd12*Zdd21
     &      *Zdd22-4D0*Zdd11/CW**2*Zdd12*Zdd21*Zdd22*sb**2
     &      +36D0*Mq**2/MW**2/SW**2*Zdd11*Zdd12*Zdd21*Zdd22)
      AAABR(1191) = EE**2/cb**2*(Zdd11/CW**2*Zdd21*Zss11
     &      *Zss21+2D0*Zdd11/CW**2*Zdd21*Zss12*Zss22+2D0
     &      *Zdd12/CW**2*Zdd22*Zss11*Zss21+4D0*Zdd12/CW**2
     &      *Zdd22*Zss12*Zss22+9D0*Zdd11/SW**2*Zdd21*Zss11
     &      *Zss21-Zdd11/CW**2*Zdd21*Zss11*Zss21*sb**2-2D0
     &      *Zdd11/CW**2*Zdd21*Zss12*Zss22*sb**2-2D0*Zdd12
     &      /CW**2*Zdd22*Zss11*Zss21*sb**2-4D0*Zdd12/CW**2
     &      *Zdd22*Zss12*Zss22*sb**2-9D0*Zdd11/SW**2*Zdd21
     &      *Zss11*Zss21*sb**2+18D0*Mq**2/MW**2/SW**2*Zdd11
     &      *Zdd22*Zss12*Zss21+18D0*Mq**2/MW**2/SW**2*Zdd12
     &      *Zdd21*Zss11*Zss22)
      AAABR(1192) = EE**2/cb**2*(Zb11/CW**2*Zb21*Zdd11*Zdd21
     &      +2D0*Zb11/CW**2*Zb21*Zdd12*Zdd22+2D0*Zb12/CW
     &      **2*Zb22*Zdd11*Zdd21+4D0*Zb12/CW**2*Zb22*Zdd12
     &      *Zdd22+9D0*Zb11/SW**2*Zb21*Zdd11*Zdd21-Zb11/CW
     &      **2*Zb21*Zdd11*Zdd21*sb**2-2D0*Zb11/CW**2*Zb21
     &      *Zdd12*Zdd22*sb**2-2D0*Zb12/CW**2*Zb22*Zdd11
     &      *Zdd21*sb**2-4D0*Zb12/CW**2*Zb22*Zdd12*Zdd22
     &      *sb**2-9D0*Zb11/SW**2*Zb21*Zdd11*Zdd21*sb**2
     &      +18D0*MbMM/MW**2*Mq/SW**2*Zb11*Zb22*Zdd12*Zdd21
     &      +18D0*MbMM/MW**2*Mq/SW**2*Zb12*Zb21*Zdd11*Zdd22)
      AAABR(1193) = GG**2/cb**2*(Zss11**2*Zss21**2+Zss12
     &      **2*Zss22**2-Zss11**2*Zss21**2*sb**2-Zss12**2
     &      *Zss22**2*sb**2+2D0*Zss11*Zss12*Zss21*Zss22-2D0
     &      *Zss11*Zss12*Zss21*Zss22*sb**2)
      AAABR(1194) = EE**2/cb**2*(Zss11**2/CW**2*Zss21**2
     &      +4D0*Zss12**2/CW**2*Zss22**2+9D0*Zss11**2/SW
     &      **2*Zss21**2-Zss11**2/CW**2*Zss21**2*sb**2-4D0
     &      *Zss12**2/CW**2*Zss22**2*sb**2-9D0*Zss11**2/SW
     &      **2*Zss21**2*sb**2+4D0*Zss11/CW**2*Zss12*Zss21
     &      *Zss22-4D0*Zss11/CW**2*Zss12*Zss21*Zss22*sb**2
     &      +36D0*Mq**2/MW**2/SW**2*Zss11*Zss12*Zss21*Zss22)
      AAABR(1195) = EE**2/cb**2*(Zdd11/CW**2*Zdd21*Zss11
     &      *Zss21+2D0*Zdd11/CW**2*Zdd21*Zss12*Zss22+2D0
     &      *Zdd12/CW**2*Zdd22*Zss11*Zss21+4D0*Zdd12/CW**2
     &      *Zdd22*Zss12*Zss22+9D0*Zdd11/SW**2*Zdd21*Zss11
     &      *Zss21-Zdd11/CW**2*Zdd21*Zss11*Zss21*sb**2-2D0
     &      *Zdd11/CW**2*Zdd21*Zss12*Zss22*sb**2-2D0*Zdd12
     &      /CW**2*Zdd22*Zss11*Zss21*sb**2-4D0*Zdd12/CW**2
     &      *Zdd22*Zss12*Zss22*sb**2-9D0*Zdd11/SW**2*Zdd21
     &      *Zss11*Zss21*sb**2+18D0*Mq**2/MW**2/SW**2*Zdd11
     &      *Zdd22*Zss11*Zss22+18D0*Mq**2/MW**2/SW**2*Zdd12
     &      *Zdd21*Zss12*Zss21)
      AAABR(1196) = GG**2/cb**2*(Zdd11**2*Zdd21**2+Zdd11
     &      **2*Zdd22**2+Zdd12**2*Zdd21**2+Zdd12**2*Zdd22
     &      **2-Zdd11**2*Zdd21**2*sb**2-Zdd11**2*Zdd22**2
     &      *sb**2-Zdd12**2*Zdd21**2*sb**2-Zdd12**2*Zdd22
     &      **2*sb**2)
      AAABR(1197) = EE**2/cb**2*(Zdd11**2/CW**2*Zdd21**2
     &      +2D0*Zdd11**2/CW**2*Zdd22**2+2D0*Zdd12**2/CW
     &      **2*Zdd21**2+4D0*Zdd12**2/CW**2*Zdd22**2+9D0
     &      *Zdd11**2/SW**2*Zdd21**2-Zdd11**2/CW**2*Zdd21
     &      **2*sb**2-2D0*Zdd11**2/CW**2*Zdd22**2*sb**2-2D0
     &      *Zdd12**2/CW**2*Zdd21**2*sb**2-4D0*Zdd12**2/CW
     &      **2*Zdd22**2*sb**2-9D0*Zdd11**2/SW**2*Zdd21**2
     &      *sb**2+36D0*Mq**2/MW**2/SW**2*Zdd11*Zdd12*Zdd21
     &      *Zdd22)
      AAABR(1198) = EE**2/cb**2*(Zdd11**2/CW**2*Zdd21**2
     &      +4D0*Zdd12**2/CW**2*Zdd22**2+9D0*Zdd11**2/SW
     &      **2*Zdd21**2-Zdd11**2/CW**2*Zdd21**2*sb**2-4D0
     &      *Zdd12**2/CW**2*Zdd22**2*sb**2-9D0*Zdd11**2/SW
     &      **2*Zdd21**2*sb**2+4D0*Zdd11/CW**2*Zdd12*Zdd21
     &      *Zdd22+18D0*Mq**2/MW**2/SW**2*Zdd11**2*Zdd22
     &      **2+18D0*Mq**2/MW**2/SW**2*Zdd12**2*Zdd21**2
     &      -4D0*Zdd11/CW**2*Zdd12*Zdd21*Zdd22*sb**2)
      AAABR(1199) = EE**2/cb**2*(Zb21**2/CW**2*Zss11**2+2D0
     &      *Zb21**2/CW**2*Zss12**2+2D0*Zb22**2/CW**2*Zss11
     &      **2+4D0*Zb22**2/CW**2*Zss12**2+9D0*Zb21**2/SW
     &      **2*Zss11**2-Zb21**2/CW**2*Zss11**2*sb**2-2D0
     &      *Zb21**2/CW**2*Zss12**2*sb**2-2D0*Zb22**2/CW
     &      **2*Zss11**2*sb**2-4D0*Zb22**2/CW**2*Zss12**2
     &      *sb**2-9D0*Zb21**2/SW**2*Zss11**2*sb**2+36D0
     &      *MbMM/MW**2*Mq/SW**2*Zb21*Zb22*Zss11*Zss12)
      AAABR(1200) = EE**2/cb**2*(Zb11/CW**2*Zb21*Zdd11*Zdd21
     &      +2D0*Zb11/CW**2*Zb21*Zdd12*Zdd22+2D0*Zb12/CW
     &      **2*Zb22*Zdd11*Zdd21+4D0*Zb12/CW**2*Zb22*Zdd12
     &      *Zdd22+9D0*Zb11/SW**2*Zb21*Zdd11*Zdd21-Zb11/CW
     &      **2*Zb21*Zdd11*Zdd21*sb**2-2D0*Zb11/CW**2*Zb21
     &      *Zdd12*Zdd22*sb**2-2D0*Zb12/CW**2*Zb22*Zdd11
     &      *Zdd21*sb**2-4D0*Zb12/CW**2*Zb22*Zdd12*Zdd22
     &      *sb**2-9D0*Zb11/SW**2*Zb21*Zdd11*Zdd21*sb**2
     &      +18D0*MbMM/MW**2*Mq/SW**2*Zb11*Zb22*Zdd11*Zdd22
     &      +18D0*MbMM/MW**2*Mq/SW**2*Zb12*Zb21*Zdd12*Zdd21)
      AAABR(1201) = EE**2/cb**2*(Zdd11**2/CW**2*Zss21**2
     &      +2D0*Zdd11**2/CW**2*Zss22**2+2D0*Zdd12**2/CW
     &      **2*Zss21**2+4D0*Zdd12**2/CW**2*Zss22**2+9D0
     &      *Zdd11**2/SW**2*Zss21**2-Zdd11**2/CW**2*Zss21
     &      **2*sb**2-2D0*Zdd11**2/CW**2*Zss22**2*sb**2-2D0
     &      *Zdd12**2/CW**2*Zss21**2*sb**2-4D0*Zdd12**2/CW
     &      **2*Zss22**2*sb**2-9D0*Zdd11**2/SW**2*Zss21**2
     &      *sb**2+36D0*Mq**2/MW**2/SW**2*Zdd11*Zdd12*Zss21
     &      *Zss22)
      AAABR(1202) = EE**2/cb**2*(Zb11/CW**2*Zb21*Zss11*Zss21
     &      +2D0*Zb11/CW**2*Zb21*Zss12*Zss22+2D0*Zb12/CW
     &      **2*Zb22*Zss11*Zss21+4D0*Zb12/CW**2*Zb22*Zss12
     &      *Zss22+9D0*Zb11/SW**2*Zb21*Zss11*Zss21-Zb11/CW
     &      **2*Zb21*Zss11*Zss21*sb**2-2D0*Zb11/CW**2*Zb21
     &      *Zss12*Zss22*sb**2-2D0*Zb12/CW**2*Zb22*Zss11
     &      *Zss21*sb**2-4D0*Zb12/CW**2*Zb22*Zss12*Zss22
     &      *sb**2-9D0*Zb11/SW**2*Zb21*Zss11*Zss21*sb**2
     &      +18D0*MbMM/MW**2*Mq/SW**2*Zb11*Zb22*Zss11*Zss22
     &      +18D0*MbMM/MW**2*Mq/SW**2*Zb12*Zb21*Zss12*Zss21)
      AAABR(1203) = EE**2/cb**2*(Zdd21**2/CW**2*Zss11**2
     &      +2D0*Zdd21**2/CW**2*Zss12**2+2D0*Zdd22**2/CW
     &      **2*Zss11**2+4D0*Zdd22**2/CW**2*Zss12**2+9D0
     &      *Zdd21**2/SW**2*Zss11**2-Zdd21**2/CW**2*Zss11
     &      **2*sb**2-2D0*Zdd21**2/CW**2*Zss12**2*sb**2-2D0
     &      *Zdd22**2/CW**2*Zss11**2*sb**2-4D0*Zdd22**2/CW
     &      **2*Zss12**2*sb**2-9D0*Zdd21**2/SW**2*Zss11**2
     &      *sb**2+36D0*Mq**2/MW**2/SW**2*Zdd21*Zdd22*Zss11
     &      *Zss12)
      AAABR(1204) = EE**2/cb**2*(Zb21**2/CW**2*Zdd11**2+2D0
     &      *Zb21**2/CW**2*Zdd12**2+2D0*Zb22**2/CW**2*Zdd11
     &      **2+4D0*Zb22**2/CW**2*Zdd12**2+9D0*Zb21**2/SW
     &      **2*Zdd11**2-Zb21**2/CW**2*Zdd11**2*sb**2-2D0
     &      *Zb21**2/CW**2*Zdd12**2*sb**2-2D0*Zb22**2/CW
     &      **2*Zdd11**2*sb**2-4D0*Zb22**2/CW**2*Zdd12**2
     &      *sb**2-9D0*Zb21**2/SW**2*Zdd11**2*sb**2+36D0
     &      *MbMM/MW**2*Mq/SW**2*Zb21*Zb22*Zdd11*Zdd12)
      AAABR(1205) = GG**2/cb**2*(Zb11**2*Zb21**2+Zb11**2
     &      *Zb22**2+Zb12**2*Zb21**2+Zb12**2*Zb22**2-Zb11
     &      **2*Zb21**2*sb**2-Zb11**2*Zb22**2*sb**2-Zb12
     &      **2*Zb21**2*sb**2-Zb12**2*Zb22**2*sb**2)
      AAABR(1206) = EE**2/cb**2*(Zb11**2/CW**2*Zb21**2+2D0
     &      *Zb11**2/CW**2*Zb22**2+2D0*Zb12**2/CW**2*Zb21
     &      **2+4D0*Zb12**2/CW**2*Zb22**2+9D0*Zb11**2/SW
     &      **2*Zb21**2-Zb11**2/CW**2*Zb21**2*sb**2-2D0*Zb11
     &      **2/CW**2*Zb22**2*sb**2-2D0*Zb12**2/CW**2*Zb21
     &      **2*sb**2-4D0*Zb12**2/CW**2*Zb22**2*sb**2-9D0
     &      *Zb11**2/SW**2*Zb21**2*sb**2+36D0*MbMM**2/MW
     &      **2/SW**2*Zb11*Zb12*Zb21*Zb22)
      AAABR(1207) = EE**2/cb**2*(Zb11**2/CW**2*Zb21**2+4D0
     &      *Zb12**2/CW**2*Zb22**2+9D0*Zb11**2/SW**2*Zb21
     &      **2-Zb11**2/CW**2*Zb21**2*sb**2-4D0*Zb12**2/CW
     &      **2*Zb22**2*sb**2-9D0*Zb11**2/SW**2*Zb21**2*sb
     &      **2+4D0*Zb11/CW**2*Zb12*Zb21*Zb22+18D0*MbMM**2
     &      /MW**2/SW**2*Zb11**2*Zb22**2+18D0*MbMM**2/MW
     &      **2/SW**2*Zb12**2*Zb21**2-4D0*Zb11/CW**2*Zb12
     &      *Zb21*Zb22*sb**2)
      AAABR(1208) = EE**2/cb**2*(Zb11**2/CW**2*Zdd21**2+2D0
     &      *Zb11**2/CW**2*Zdd22**2+2D0*Zb12**2/CW**2*Zdd21
     &      **2+4D0*Zb12**2/CW**2*Zdd22**2+9D0*Zb11**2/SW
     &      **2*Zdd21**2-Zb11**2/CW**2*Zdd21**2*sb**2-2D0
     &      *Zb11**2/CW**2*Zdd22**2*sb**2-2D0*Zb12**2/CW
     &      **2*Zdd21**2*sb**2-4D0*Zb12**2/CW**2*Zdd22**2
     &      *sb**2-9D0*Zb11**2/SW**2*Zdd21**2*sb**2+36D0
     &      *MbMM/MW**2*Mq/SW**2*Zb11*Zb12*Zdd21*Zdd22)
      AAABR(1209) = GG**2/cb**2*(Zss11**2*Zss21**2+Zss11
     &      **2*Zss22**2+Zss12**2*Zss21**2+Zss12**2*Zss22
     &      **2-Zss11**2*Zss21**2*sb**2-Zss11**2*Zss22**2
     &      *sb**2-Zss12**2*Zss21**2*sb**2-Zss12**2*Zss22
     &      **2*sb**2)
      AAABR(1210) = EE**2/cb**2*(Zss11**2/CW**2*Zss21**2
     &      +2D0*Zss11**2/CW**2*Zss22**2+2D0*Zss12**2/CW
     &      **2*Zss21**2+4D0*Zss12**2/CW**2*Zss22**2+9D0
     &      *Zss11**2/SW**2*Zss21**2-Zss11**2/CW**2*Zss21
     &      **2*sb**2-2D0*Zss11**2/CW**2*Zss22**2*sb**2-2D0
     &      *Zss12**2/CW**2*Zss21**2*sb**2-4D0*Zss12**2/CW
     &      **2*Zss22**2*sb**2-9D0*Zss11**2/SW**2*Zss21**2
     &      *sb**2+36D0*Mq**2/MW**2/SW**2*Zss11*Zss12*Zss21
     &      *Zss22)
      AAABR(1211) = EE**2/cb**2*(Zss11**2/CW**2*Zss21**2
     &      +4D0*Zss12**2/CW**2*Zss22**2+9D0*Zss11**2/SW
     &      **2*Zss21**2-Zss11**2/CW**2*Zss21**2*sb**2-4D0
     &      *Zss12**2/CW**2*Zss22**2*sb**2-9D0*Zss11**2/SW
     &      **2*Zss21**2*sb**2+4D0*Zss11/CW**2*Zss12*Zss21
     &      *Zss22+18D0*Mq**2/MW**2/SW**2*Zss11**2*Zss22
     &      **2+18D0*Mq**2/MW**2/SW**2*Zss12**2*Zss21**2
     &      -4D0*Zss11/CW**2*Zss12*Zss21*Zss22*sb**2)
      AAABR(1212) = EE**2/cb**2*(Zb11**2/CW**2*Zss21**2+2D0
     &      *Zb11**2/CW**2*Zss22**2+2D0*Zb12**2/CW**2*Zss21
     &      **2+4D0*Zb12**2/CW**2*Zss22**2+9D0*Zb11**2/SW
     &      **2*Zss21**2-Zb11**2/CW**2*Zss21**2*sb**2-2D0
     &      *Zb11**2/CW**2*Zss22**2*sb**2-2D0*Zb12**2/CW
     &      **2*Zss21**2*sb**2-4D0*Zb12**2/CW**2*Zss22**2
     &      *sb**2-9D0*Zb11**2/SW**2*Zss21**2*sb**2+36D0
     &      *MbMM/MW**2*Mq/SW**2*Zb11*Zb12*Zss21*Zss22)
      AAABR(1213) = GG**2/cb**2*(Zss11*Zss21**3+Zss12*Zss22
     &      **3+Zss11*Zss21*Zss22**2-Zss11*Zss21**3*sb**2
     &      +Zss12*Zss21**2*Zss22-Zss12*Zss22**3*sb**2-Zss11
     &      *Zss21*Zss22**2*sb**2-Zss12*Zss21**2*Zss22*sb
     &      **2)
      AAABR(1214) = EE**2/cb**2*(Zss11/CW**2*Zss21**3+4D0
     &      *Zss12/CW**2*Zss22**3+9D0*Zss11/SW**2*Zss21**3
     &      +2D0*Zss11/CW**2*Zss21*Zss22**2-Zss11/CW**2*Zss21
     &      **3*sb**2+2D0*Zss12/CW**2*Zss21**2*Zss22-4D0
     &      *Zss12/CW**2*Zss22**3*sb**2-9D0*Zss11/SW**2*Zss21
     &      **3*sb**2-2D0*Zss11/CW**2*Zss21*Zss22**2*sb**2
     &      -2D0*Zss12/CW**2*Zss21**2*Zss22*sb**2+18D0*Mq
     &      **2/MW**2/SW**2*Zss11*Zss21*Zss22**2+18D0*Mq
     &      **2/MW**2/SW**2*Zss12*Zss21**2*Zss22)
      AAABR(1215) = EE**2/cb**2*(Zb11/CW**2*Zb21*Zdd21**2
     &      +2D0*Zb11/CW**2*Zb21*Zdd22**2+2D0*Zb12/CW**2
     &      *Zb22*Zdd21**2+4D0*Zb12/CW**2*Zb22*Zdd22**2+9D0
     &      *Zb11/SW**2*Zb21*Zdd21**2-Zb11/CW**2*Zb21*Zdd21
     &      **2*sb**2-2D0*Zb11/CW**2*Zb21*Zdd22**2*sb**2
     &      -2D0*Zb12/CW**2*Zb22*Zdd21**2*sb**2-4D0*Zb12
     &      /CW**2*Zb22*Zdd22**2*sb**2-9D0*Zb11/SW**2*Zb21
     &      *Zdd21**2*sb**2+18D0*MbMM/MW**2*Mq/SW**2*Zb11
     &      *Zb22*Zdd21*Zdd22+18D0*MbMM/MW**2*Mq/SW**2*Zb12
     &      *Zb21*Zdd21*Zdd22)
      AAABR(1216) = GG**2/cb**2*(Zb11*Zb21**3+Zb12*Zb22**3
     &      +Zb11*Zb21*Zb22**2-Zb11*Zb21**3*sb**2+Zb12*Zb21
     &      **2*Zb22-Zb12*Zb22**3*sb**2-Zb11*Zb21*Zb22**2
     &      *sb**2-Zb12*Zb21**2*Zb22*sb**2)
      AAABR(1217) = EE**2/cb**2*(Zb11/CW**2*Zb21**3+4D0*Zb12
     &      /CW**2*Zb22**3+9D0*Zb11/SW**2*Zb21**3+2D0*Zb11
     &      /CW**2*Zb21*Zb22**2-Zb11/CW**2*Zb21**3*sb**2
     &      +2D0*Zb12/CW**2*Zb21**2*Zb22-4D0*Zb12/CW**2*Zb22
     &      **3*sb**2-9D0*Zb11/SW**2*Zb21**3*sb**2-2D0*Zb11
     &      /CW**2*Zb21*Zb22**2*sb**2-2D0*Zb12/CW**2*Zb21
     &      **2*Zb22*sb**2+18D0*MbMM**2/MW**2/SW**2*Zb11
     &      *Zb21*Zb22**2+18D0*MbMM**2/MW**2/SW**2*Zb12*Zb21
     &      **2*Zb22)
      AAABR(1218) = EE**2/cb**2*(Zb21**2/CW**2*Zss11*Zss21
     &      +2D0*Zb21**2/CW**2*Zss12*Zss22+2D0*Zb22**2/CW
     &      **2*Zss11*Zss21+4D0*Zb22**2/CW**2*Zss12*Zss22
     &      +9D0*Zb21**2/SW**2*Zss11*Zss21-Zb21**2/CW**2
     &      *Zss11*Zss21*sb**2-2D0*Zb21**2/CW**2*Zss12*Zss22
     &      *sb**2-2D0*Zb22**2/CW**2*Zss11*Zss21*sb**2-4D0
     &      *Zb22**2/CW**2*Zss12*Zss22*sb**2-9D0*Zb21**2
     &      /SW**2*Zss11*Zss21*sb**2+18D0*MbMM/MW**2*Mq/SW
     &      **2*Zb21*Zb22*Zss11*Zss22+18D0*MbMM/MW**2*Mq
     &      /SW**2*Zb21*Zb22*Zss12*Zss21)
      AAABR(1219) = EE**2/cb**2*(Zdd11/CW**2*Zdd21*Zss21
     &      **2+2D0*Zdd11/CW**2*Zdd21*Zss22**2+2D0*Zdd12
     &      /CW**2*Zdd22*Zss21**2+4D0*Zdd12/CW**2*Zdd22*Zss22
     &      **2+9D0*Zdd11/SW**2*Zdd21*Zss21**2-Zdd11/CW**2
     &      *Zdd21*Zss21**2*sb**2-2D0*Zdd11/CW**2*Zdd21*Zss22
     &      **2*sb**2-2D0*Zdd12/CW**2*Zdd22*Zss21**2*sb**2
     &      -4D0*Zdd12/CW**2*Zdd22*Zss22**2*sb**2-9D0*Zdd11
     &      /SW**2*Zdd21*Zss21**2*sb**2+18D0*Mq**2/MW**2
     &      /SW**2*Zdd11*Zdd22*Zss21*Zss22+18D0*Mq**2/MW
     &      **2/SW**2*Zdd12*Zdd21*Zss21*Zss22)
      AAABR(1220) = EE**2/cb**2*(Zb11/CW**2*Zb21*Zss21**2
     &      +2D0*Zb11/CW**2*Zb21*Zss22**2+2D0*Zb12/CW**2
     &      *Zb22*Zss21**2+4D0*Zb12/CW**2*Zb22*Zss22**2+9D0
     &      *Zb11/SW**2*Zb21*Zss21**2-Zb11/CW**2*Zb21*Zss21
     &      **2*sb**2-2D0*Zb11/CW**2*Zb21*Zss22**2*sb**2
     &      -2D0*Zb12/CW**2*Zb22*Zss21**2*sb**2-4D0*Zb12
     &      /CW**2*Zb22*Zss22**2*sb**2-9D0*Zb11/SW**2*Zb21
     &      *Zss21**2*sb**2+18D0*MbMM/MW**2*Mq/SW**2*Zb11
     &      *Zb22*Zss21*Zss22+18D0*MbMM/MW**2*Mq/SW**2*Zb12
     &      *Zb21*Zss21*Zss22)
      AAABR(1221) = EE**2/cb**2*(Zdd21**2/CW**2*Zss11*Zss21
     &      +2D0*Zdd21**2/CW**2*Zss12*Zss22+2D0*Zdd22**2
     &      /CW**2*Zss11*Zss21+4D0*Zdd22**2/CW**2*Zss12*Zss22
     &      +9D0*Zdd21**2/SW**2*Zss11*Zss21-Zdd21**2/CW**2
     &      *Zss11*Zss21*sb**2-2D0*Zdd21**2/CW**2*Zss12*Zss22
     &      *sb**2-2D0*Zdd22**2/CW**2*Zss11*Zss21*sb**2-4D0
     &      *Zdd22**2/CW**2*Zss12*Zss22*sb**2-9D0*Zdd21**2
     &      /SW**2*Zss11*Zss21*sb**2+18D0*Mq**2/MW**2/SW
     &      **2*Zdd21*Zdd22*Zss11*Zss22+18D0*Mq**2/MW**2
     &      /SW**2*Zdd21*Zdd22*Zss12*Zss21)
      AAABR(1222) = GG**2/cb**2*(Zdd11*Zdd21**3+Zdd12*Zdd22
     &      **3+Zdd11*Zdd21*Zdd22**2-Zdd11*Zdd21**3*sb**2
     &      +Zdd12*Zdd21**2*Zdd22-Zdd12*Zdd22**3*sb**2-Zdd11
     &      *Zdd21*Zdd22**2*sb**2-Zdd12*Zdd21**2*Zdd22*sb
     &      **2)
      AAABR(1223) = EE**2/cb**2*(Zdd11/CW**2*Zdd21**3+4D0
     &      *Zdd12/CW**2*Zdd22**3+9D0*Zdd11/SW**2*Zdd21**3
     &      +2D0*Zdd11/CW**2*Zdd21*Zdd22**2-Zdd11/CW**2*Zdd21
     &      **3*sb**2+2D0*Zdd12/CW**2*Zdd21**2*Zdd22-4D0
     &      *Zdd12/CW**2*Zdd22**3*sb**2-9D0*Zdd11/SW**2*Zdd21
     &      **3*sb**2-2D0*Zdd11/CW**2*Zdd21*Zdd22**2*sb**2
     &      -2D0*Zdd12/CW**2*Zdd21**2*Zdd22*sb**2+18D0*Mq
     &      **2/MW**2/SW**2*Zdd11*Zdd21*Zdd22**2+18D0*Mq
     &      **2/MW**2/SW**2*Zdd12*Zdd21**2*Zdd22)
      AAABR(1224) = EE**2/cb**2*(Zb21**2/CW**2*Zdd11*Zdd21
     &      +2D0*Zb21**2/CW**2*Zdd12*Zdd22+2D0*Zb22**2/CW
     &      **2*Zdd11*Zdd21+4D0*Zb22**2/CW**2*Zdd12*Zdd22
     &      +9D0*Zb21**2/SW**2*Zdd11*Zdd21-Zb21**2/CW**2
     &      *Zdd11*Zdd21*sb**2-2D0*Zb21**2/CW**2*Zdd12*Zdd22
     &      *sb**2-2D0*Zb22**2/CW**2*Zdd11*Zdd21*sb**2-4D0
     &      *Zb22**2/CW**2*Zdd12*Zdd22*sb**2-9D0*Zb21**2
     &      /SW**2*Zdd11*Zdd21*sb**2+18D0*MbMM/MW**2*Mq/SW
     &      **2*Zb21*Zb22*Zdd11*Zdd22+18D0*MbMM/MW**2*Mq
     &      /SW**2*Zb21*Zb22*Zdd12*Zdd21)
      AAABR(1225) = EE**2/SW**2/cb**2*(Zcc11*Zdd11*Zss11
     &      *Zuu11-Zcc11*Zdd11*Zss11*Zuu11*sb**2+Mq**2/MW
     &      **2*Zcc11*Zdd12*Zss12*Zuu11-Mc/MW**2*Mq*Zcc12
     &      *Zdd11*Zss11*Zuu12+Mc/MW**2*Mq*Zcc12*Zdd11*Zss11
     &      *Zuu12/sb**2)
      AAABR(1226) = EE**2*(Zb11**2/CW**2*Zcc11**2-4D0*Zb11
     &      **2/CW**2*Zcc12**2+2D0*Zb12**2/CW**2*Zcc11**2
     &      -8D0*Zb12**2/CW**2*Zcc12**2-9D0*Zb11**2/SW**2
     &      *Zcc11**2)
      AAABR(1227) = GG**2
      AAABR(1228) = EE**2/SW**2/cb**2*(Zb11*Zdd11*Zt11*Zuu11
     &      -Zb11*Zdd11*Zt11*Zuu11*sb**2+MbMM/MW**2*Mq*Zb12
     &      *Zdd12*Zt11*Zuu11-Mq/MW**2*MtMM*Zb11*Zdd11*Zt12
     &      *Zuu12+Mq/MW**2*MtMM*Zb11*Zdd11*Zt12*Zuu12/sb
     &      **2)
      AAABR(1229) = EE**2*(Zcc11**2/CW**2*Zdd11**2+2D0*Zcc11
     &      **2/CW**2*Zdd12**2-4D0*Zcc12**2/CW**2*Zdd11**2
     &      -8D0*Zcc12**2/CW**2*Zdd12**2-9D0*Zcc11**2/SW
     &      **2*Zdd11**2)
      AAABR(1230) = EE**2*(Zb11**2/CW**2*Zuu11**2-4D0*Zb11
     &      **2/CW**2*Zuu12**2+2D0*Zb12**2/CW**2*Zuu11**2
     &      -8D0*Zb12**2/CW**2*Zuu12**2-9D0*Zb11**2/SW**2
     &      *Zuu11**2)
      AAABR(1231) = EE**2/cb**2*(Zcc11**2/CW**2*Zss11**2
     &      +2D0*Zcc11**2/CW**2*Zss12**2-4D0*Zcc12**2/CW
     &      **2*Zss11**2-8D0*Zcc12**2/CW**2*Zss12**2-9D0
     &      *Zcc11**2/SW**2*Zss11**2-Zcc11**2/CW**2*Zss11
     &      **2*sb**2-2D0*Zcc11**2/CW**2*Zss12**2*sb**2+4D0
     &      *Zcc12**2/CW**2*Zss11**2*sb**2+8D0*Zcc12**2/CW
     &      **2*Zss12**2*sb**2+9D0*Zcc11**2/SW**2*Zss11**2
     &      *sb**2)
      AAABR(1232) = EE**2/SW**2/cb**2*(Zcc11**2*Zss11**2
     &      -Zcc11**2*Zss11**2*sb**2-Mc**2/MW**2*Zcc12**2
     &      *Zss11**2+Mq**2/MW**2*Zcc11**2*Zss12**2+Mc**2
     &      /MW**2*Zcc12**2*Zss11**2/sb**2)
      AAABR(1233) = EE**2*(Zdd11**2/CW**2*Zt11**2-4D0*Zdd11
     &      **2/CW**2*Zt12**2+2D0*Zdd12**2/CW**2*Zt11**2
     &      -8D0*Zdd12**2/CW**2*Zt12**2-9D0*Zdd11**2/SW**2
     &      *Zt11**2)
      AAABR(1234) = EE**2/SW**2/cb**2*(Zb11*Zcc11*Zss11*Zt11
     &      -Zb11*Zcc11*Zss11*Zt11*sb**2+MbMM/MW**2*Mq*Zb12
     &      *Zcc11*Zss12*Zt11-Mc/MW**2*MtMM*Zb11*Zcc12*Zss11
     &      *Zt12+Mc/MW**2*MtMM*Zb11*Zcc12*Zss11*Zt12/sb**2)
      AAABR(1235) = EE**2*(Zss11**2/CW**2*Zuu11**2-4D0*Zss11
     &      **2/CW**2*Zuu12**2+2D0*Zss12**2/CW**2*Zuu11**2
     &      -8D0*Zss12**2/CW**2*Zuu12**2-9D0*Zss11**2/SW
     &      **2*Zuu11**2)
      AAABR(1236) = EE**2/cb**2*(Zdd11**2/CW**2*Zuu11**2
     &      -4D0*Zdd11**2/CW**2*Zuu12**2+2D0*Zdd12**2/CW
     &      **2*Zuu11**2-8D0*Zdd12**2/CW**2*Zuu12**2-9D0
     &      *Zdd11**2/SW**2*Zuu11**2-Zdd11**2/CW**2*Zuu11
     &      **2*sb**2+4D0*Zdd11**2/CW**2*Zuu12**2*sb**2-2D0
     &      *Zdd12**2/CW**2*Zuu11**2*sb**2+8D0*Zdd12**2/CW
     &      **2*Zuu12**2*sb**2+9D0*Zdd11**2/SW**2*Zuu11**2
     &      *sb**2)
      AAABR(1237) = EE**2/SW**2/cb**2*(Zdd11**2*Zuu11**2
     &      -Zdd11**2*Zuu11**2*sb**2-Mq**2/MW**2*Zdd11**2
     &      *Zuu12**2+Mq**2/MW**2*Zdd12**2*Zuu11**2+Mq**2
     &      /MW**2*Zdd11**2*Zuu12**2/sb**2)
      AAABR(1238) = EE**2/cb**2*(Zb11**2/CW**2*Zt11**2-4D0
     &      *Zb11**2/CW**2*Zt12**2+2D0*Zb12**2/CW**2*Zt11
     &      **2-8D0*Zb12**2/CW**2*Zt12**2-9D0*Zb11**2/SW
     &      **2*Zt11**2-Zb11**2/CW**2*Zt11**2*sb**2+4D0*Zb11
     &      **2/CW**2*Zt12**2*sb**2-2D0*Zb12**2/CW**2*Zt11
     &      **2*sb**2+8D0*Zb12**2/CW**2*Zt12**2*sb**2+9D0
     &      *Zb11**2/SW**2*Zt11**2*sb**2)
      AAABR(1239) = EE**2/SW**2/cb**2*(Zb11**2*Zt11**2-Zb11
     &      **2*Zt11**2*sb**2+MbMM**2/MW**2*Zb12**2*Zt11
     &      **2-MtMM**2/MW**2*Zb11**2*Zt12**2+MtMM**2/MW
     &      **2*Zb11**2*Zt12**2/sb**2)
      AAABR(1240) = EE**2*(Zss11**2/CW**2*Zt11**2-4D0*Zss11
     &      **2/CW**2*Zt12**2+2D0*Zss12**2/CW**2*Zt11**2
     &      -8D0*Zss12**2/CW**2*Zt12**2-9D0*Zss11**2/SW**2
     &      *Zt11**2)
      AAABR(1241) = EE**2/SW**2/cb**2*(Zb11*Zdd11*Zt11*Zuu21
     &      -Zb11*Zdd11*Zt11*Zuu21*sb**2+MbMM/MW**2*Mq*Zb12
     &      *Zdd12*Zt11*Zuu21-Mq/MW**2*MtMM*Zb11*Zdd11*Zt12
     &      *Zuu22+Mq/MW**2*MtMM*Zb11*Zdd11*Zt12*Zuu22/sb
     &      **2)
      AAABR(1242) = EE**2*(Zb11**2/CW**2*Zcc11*Zcc21-4D0
     &      *Zb11**2/CW**2*Zcc12*Zcc22+2D0*Zb12**2/CW**2
     &      *Zcc11*Zcc21-8D0*Zb12**2/CW**2*Zcc12*Zcc22-9D0
     &      *Zb11**2/SW**2*Zcc11*Zcc21)
      AAABR(1243) = EE**2/SW**2/cb**2*(Zcc21*Zdd11*Zss11
     &      *Zuu11-Zcc21*Zdd11*Zss11*Zuu11*sb**2+Mq**2/MW
     &      **2*Zcc21*Zdd12*Zss12*Zuu11-Mc/MW**2*Mq*Zcc22
     &      *Zdd11*Zss11*Zuu12+Mc/MW**2*Mq*Zcc22*Zdd11*Zss11
     &      *Zuu12/sb**2)
      AAABR(1244) = EE**2/cb**2*(Zdd11**2/CW**2*Zuu11*Zuu21
     &      -4D0*Zdd11**2/CW**2*Zuu12*Zuu22+2D0*Zdd12**2
     &      /CW**2*Zuu11*Zuu21-8D0*Zdd12**2/CW**2*Zuu12*Zuu22
     &      -9D0*Zdd11**2/SW**2*Zuu11*Zuu21-Zdd11**2/CW**2
     &      *Zuu11*Zuu21*sb**2+4D0*Zdd11**2/CW**2*Zuu12*Zuu22
     &      *sb**2-2D0*Zdd12**2/CW**2*Zuu11*Zuu21*sb**2+8D0
     &      *Zdd12**2/CW**2*Zuu12*Zuu22*sb**2+9D0*Zdd11**2
     &      /SW**2*Zuu11*Zuu21*sb**2)
      AAABR(1245) = EE**2/SW**2/cb**2*(Zdd11**2*Zuu11*Zuu21
     &      -Zdd11**2*Zuu11*Zuu21*sb**2-Mq**2/MW**2*Zdd11
     &      **2*Zuu12*Zuu22+Mq**2/MW**2*Zdd12**2*Zuu11*Zuu21
     &      +Mq**2/MW**2*Zdd11**2*Zuu12*Zuu22/sb**2)
      AAABR(1246) = EE**2*(Zb11**2/CW**2*Zuu11*Zuu21-4D0
     &      *Zb11**2/CW**2*Zuu12*Zuu22+2D0*Zb12**2/CW**2
     &      *Zuu11*Zuu21-8D0*Zb12**2/CW**2*Zuu12*Zuu22-9D0
     &      *Zb11**2/SW**2*Zuu11*Zuu21)
      AAABR(1247) = EE**2*(Zdd11**2/CW**2*Zt11*Zt21-4D0*Zdd11
     &      **2/CW**2*Zt12*Zt22+2D0*Zdd12**2/CW**2*Zt11*Zt21
     &      -8D0*Zdd12**2/CW**2*Zt12*Zt22-9D0*Zdd11**2/SW
     &      **2*Zt11*Zt21)
      AAABR(1248) = EE**2/cb**2*(Zcc11/CW**2*Zcc21*Zss11
     &      **2+2D0*Zcc11/CW**2*Zcc21*Zss12**2-4D0*Zcc12
     &      /CW**2*Zcc22*Zss11**2-8D0*Zcc12/CW**2*Zcc22*Zss12
     &      **2-9D0*Zcc11/SW**2*Zcc21*Zss11**2-Zcc11/CW**2
     &      *Zcc21*Zss11**2*sb**2-2D0*Zcc11/CW**2*Zcc21*Zss12
     &      **2*sb**2+4D0*Zcc12/CW**2*Zcc22*Zss11**2*sb**2
     &      +8D0*Zcc12/CW**2*Zcc22*Zss12**2*sb**2+9D0*Zcc11
     &      /SW**2*Zcc21*Zss11**2*sb**2)
      AAABR(1249) = EE**2/SW**2/cb**2*(Zcc11*Zcc21*Zss11
     &      **2-Zcc11*Zcc21*Zss11**2*sb**2-Mc**2/MW**2*Zcc12
     &      *Zcc22*Zss11**2+Mq**2/MW**2*Zcc11*Zcc21*Zss12
     &      **2+Mc**2/MW**2*Zcc12*Zcc22*Zss11**2/sb**2)
      AAABR(1250) = EE**2*(Zss11**2/CW**2*Zuu11*Zuu21-4D0
     &      *Zss11**2/CW**2*Zuu12*Zuu22+2D0*Zss12**2/CW**2
     &      *Zuu11*Zuu21-8D0*Zss12**2/CW**2*Zuu12*Zuu22-9D0
     &      *Zss11**2/SW**2*Zuu11*Zuu21)
      AAABR(1251) = EE**2/SW**2/cb**2*(Zb11*Zcc11*Zss11*Zt21
     &      -Zb11*Zcc11*Zss11*Zt21*sb**2+MbMM/MW**2*Mq*Zb12
     &      *Zcc11*Zss12*Zt21-Mc/MW**2*MtMM*Zb11*Zcc12*Zss11
     &      *Zt22+Mc/MW**2*MtMM*Zb11*Zcc12*Zss11*Zt22/sb**2)
      AAABR(1252) = EE**2/SW**2/cb**2*(Zcc11*Zdd11*Zss11
     &      *Zuu21-Zcc11*Zdd11*Zss11*Zuu21*sb**2+Mq**2/MW
     &      **2*Zcc11*Zdd12*Zss12*Zuu21-Mc/MW**2*Mq*Zcc12
     &      *Zdd11*Zss11*Zuu22+Mc/MW**2*Mq*Zcc12*Zdd11*Zss11
     &      *Zuu22/sb**2)
      AAABR(1253) = EE**2*(Zcc11/CW**2*Zcc21*Zdd11**2+2D0
     &      *Zcc11/CW**2*Zcc21*Zdd12**2-4D0*Zcc12/CW**2*Zcc22
     &      *Zdd11**2-8D0*Zcc12/CW**2*Zcc22*Zdd12**2-9D0
     &      *Zcc11/SW**2*Zcc21*Zdd11**2)
      AAABR(1254) = EE**2/SW**2/cb**2*(Zb11*Zdd11*Zt21*Zuu11
     &      -Zb11*Zdd11*Zt21*Zuu11*sb**2+MbMM/MW**2*Mq*Zb12
     &      *Zdd12*Zt21*Zuu11-Mq/MW**2*MtMM*Zb11*Zdd11*Zt22
     &      *Zuu12+Mq/MW**2*MtMM*Zb11*Zdd11*Zt22*Zuu12/sb
     &      **2)
      AAABR(1255) = EE**2/cb**2*(Zb11**2/CW**2*Zt11*Zt21
     &      -4D0*Zb11**2/CW**2*Zt12*Zt22+2D0*Zb12**2/CW**2
     &      *Zt11*Zt21-8D0*Zb12**2/CW**2*Zt12*Zt22-9D0*Zb11
     &      **2/SW**2*Zt11*Zt21-Zb11**2/CW**2*Zt11*Zt21*sb
     &      **2+4D0*Zb11**2/CW**2*Zt12*Zt22*sb**2-2D0*Zb12
     &      **2/CW**2*Zt11*Zt21*sb**2+8D0*Zb12**2/CW**2*Zt12
     &      *Zt22*sb**2+9D0*Zb11**2/SW**2*Zt11*Zt21*sb**2)
      AAABR(1256) = EE**2/SW**2/cb**2*(Zb11**2*Zt11*Zt21
     &      -Zb11**2*Zt11*Zt21*sb**2+MbMM**2/MW**2*Zb12**2
     &      *Zt11*Zt21-MtMM**2/MW**2*Zb11**2*Zt12*Zt22+MtMM
     &      **2/MW**2*Zb11**2*Zt12*Zt22/sb**2)
      AAABR(1257) = EE**2*(Zss11**2/CW**2*Zt11*Zt21-4D0*Zss11
     &      **2/CW**2*Zt12*Zt22+2D0*Zss12**2/CW**2*Zt11*Zt21
     &      -8D0*Zss12**2/CW**2*Zt12*Zt22-9D0*Zss11**2/SW
     &      **2*Zt11*Zt21)
      AAABR(1258) = EE**2/SW**2/cb**2*(Zb11*Zcc21*Zss11*Zt11
     &      -Zb11*Zcc21*Zss11*Zt11*sb**2+MbMM/MW**2*Mq*Zb12
     &      *Zcc21*Zss12*Zt11-Mc/MW**2*MtMM*Zb11*Zcc22*Zss11
     &      *Zt12+Mc/MW**2*MtMM*Zb11*Zcc22*Zss11*Zt12/sb**2)
      AAABR(1259) = EE**2*(Zdd11/CW**2*Zdd21*Zt11**2-4D0
     &      *Zdd11/CW**2*Zdd21*Zt12**2+2D0*Zdd12/CW**2*Zdd22
     &      *Zt11**2-8D0*Zdd12/CW**2*Zdd22*Zt12**2-9D0*Zdd11
     &      /SW**2*Zdd21*Zt11**2)
      AAABR(1260) = EE**2*(Zcc11**2/CW**2*Zdd11*Zdd21+2D0
     &      *Zcc11**2/CW**2*Zdd12*Zdd22-4D0*Zcc12**2/CW**2
     &      *Zdd11*Zdd21-8D0*Zcc12**2/CW**2*Zdd12*Zdd22-9D0
     &      *Zcc11**2/SW**2*Zdd11*Zdd21)
      AAABR(1261) = EE**2/SW**2/cb**2*(Zb21*Zdd11*Zt11*Zuu11
     &      -Zb21*Zdd11*Zt11*Zuu11*sb**2+MbMM/MW**2*Mq*Zb22
     &      *Zdd12*Zt11*Zuu11-Mq/MW**2*MtMM*Zb21*Zdd11*Zt12
     &      *Zuu12+Mq/MW**2*MtMM*Zb21*Zdd11*Zt12*Zuu12/sb
     &      **2)
      AAABR(1262) = EE**2/SW**2/cb**2*(Zcc11*Zdd11*Zss21
     &      *Zuu11-Zcc11*Zdd11*Zss21*Zuu11*sb**2+Mq**2/MW
     &      **2*Zcc11*Zdd12*Zss22*Zuu11-Mc/MW**2*Mq*Zcc12
     &      *Zdd11*Zss21*Zuu12+Mc/MW**2*Mq*Zcc12*Zdd11*Zss21
     &      *Zuu12/sb**2)
      AAABR(1263) = EE**2/cb**2*(Zb11/CW**2*Zb21*Zt11**2
     &      -4D0*Zb11/CW**2*Zb21*Zt12**2+2D0*Zb12/CW**2*Zb22
     &      *Zt11**2-8D0*Zb12/CW**2*Zb22*Zt12**2-9D0*Zb11
     &      /SW**2*Zb21*Zt11**2-Zb11/CW**2*Zb21*Zt11**2*sb
     &      **2+4D0*Zb11/CW**2*Zb21*Zt12**2*sb**2-2D0*Zb12
     &      /CW**2*Zb22*Zt11**2*sb**2+8D0*Zb12/CW**2*Zb22
     &      *Zt12**2*sb**2+9D0*Zb11/SW**2*Zb21*Zt11**2*sb
     &      **2)
      AAABR(1264) = EE**2/SW**2/cb**2*(Zb11*Zb21*Zt11**2
     &      -Zb11*Zb21*Zt11**2*sb**2+MbMM**2/MW**2*Zb12*Zb22
     &      *Zt11**2-MtMM**2/MW**2*Zb11*Zb21*Zt12**2+MtMM
     &      **2/MW**2*Zb11*Zb21*Zt12**2/sb**2)
      AAABR(1265) = EE**2*(Zb11/CW**2*Zb21*Zuu11**2-4D0*Zb11
     &      /CW**2*Zb21*Zuu12**2+2D0*Zb12/CW**2*Zb22*Zuu11
     &      **2-8D0*Zb12/CW**2*Zb22*Zuu12**2-9D0*Zb11/SW
     &      **2*Zb21*Zuu11**2)
      AAABR(1266) = EE**2/SW**2/cb**2*(Zcc11*Zdd21*Zss11
     &      *Zuu11-Zcc11*Zdd21*Zss11*Zuu11*sb**2+Mq**2/MW
     &      **2*Zcc11*Zdd22*Zss12*Zuu11-Mc/MW**2*Mq*Zcc12
     &      *Zdd21*Zss11*Zuu12+Mc/MW**2*Mq*Zcc12*Zdd21*Zss11
     &      *Zuu12/sb**2)
      AAABR(1267) = EE**2*(Zss11/CW**2*Zss21*Zuu11**2-4D0
     &      *Zss11/CW**2*Zss21*Zuu12**2+2D0*Zss12/CW**2*Zss22
     &      *Zuu11**2-8D0*Zss12/CW**2*Zss22*Zuu12**2-9D0
     &      *Zss11/SW**2*Zss21*Zuu11**2)
      AAABR(1268) = EE**2*(Zb11/CW**2*Zb21*Zcc11**2-4D0*Zb11
     &      /CW**2*Zb21*Zcc12**2+2D0*Zb12/CW**2*Zb22*Zcc11
     &      **2-8D0*Zb12/CW**2*Zb22*Zcc12**2-9D0*Zb11/SW
     &      **2*Zb21*Zcc11**2)
      AAABR(1269) = EE**2/cb**2*(Zcc11**2/CW**2*Zss11*Zss21
     &      +2D0*Zcc11**2/CW**2*Zss12*Zss22-4D0*Zcc12**2
     &      /CW**2*Zss11*Zss21-8D0*Zcc12**2/CW**2*Zss12*Zss22
     &      -9D0*Zcc11**2/SW**2*Zss11*Zss21-Zcc11**2/CW**2
     &      *Zss11*Zss21*sb**2-2D0*Zcc11**2/CW**2*Zss12*Zss22
     &      *sb**2+4D0*Zcc12**2/CW**2*Zss11*Zss21*sb**2+8D0
     &      *Zcc12**2/CW**2*Zss12*Zss22*sb**2+9D0*Zcc11**2
     &      /SW**2*Zss11*Zss21*sb**2)
      AAABR(1270) = EE**2/SW**2/cb**2*(Zcc11**2*Zss11*Zss21
     &      -Zcc11**2*Zss11*Zss21*sb**2-Mc**2/MW**2*Zcc12
     &      **2*Zss11*Zss21+Mq**2/MW**2*Zcc11**2*Zss12*Zss22
     &      +Mc**2/MW**2*Zcc12**2*Zss11*Zss21/sb**2)
      AAABR(1271) = EE**2/SW**2/cb**2*(Zb11*Zdd21*Zt11*Zuu11
     &      -Zb11*Zdd21*Zt11*Zuu11*sb**2+MbMM/MW**2*Mq*Zb12
     &      *Zdd22*Zt11*Zuu11-Mq/MW**2*MtMM*Zb11*Zdd21*Zt12
     &      *Zuu12+Mq/MW**2*MtMM*Zb11*Zdd21*Zt12*Zuu12/sb
     &      **2)
      AAABR(1272) = EE**2*(Zss11/CW**2*Zss21*Zt11**2-4D0
     &      *Zss11/CW**2*Zss21*Zt12**2+2D0*Zss12/CW**2*Zss22
     &      *Zt11**2-8D0*Zss12/CW**2*Zss22*Zt12**2-9D0*Zss11
     &      /SW**2*Zss21*Zt11**2)
      AAABR(1273) = EE**2/SW**2/cb**2*(Zb21*Zcc11*Zss11*Zt11
     &      -Zb21*Zcc11*Zss11*Zt11*sb**2+MbMM/MW**2*Mq*Zb22
     &      *Zcc11*Zss12*Zt11-Mc/MW**2*MtMM*Zb21*Zcc12*Zss11
     &      *Zt12+Mc/MW**2*MtMM*Zb21*Zcc12*Zss11*Zt12/sb**2)
      AAABR(1274) = EE**2/SW**2/cb**2*(Zb11*Zcc11*Zss21*Zt11
     &      -Zb11*Zcc11*Zss21*Zt11*sb**2+MbMM/MW**2*Mq*Zb12
     &      *Zcc11*Zss22*Zt11-Mc/MW**2*MtMM*Zb11*Zcc12*Zss21
     &      *Zt12+Mc/MW**2*MtMM*Zb11*Zcc12*Zss21*Zt12/sb**2)
      AAABR(1275) = EE**2/cb**2*(Zdd11/CW**2*Zdd21*Zuu11
     &      **2-4D0*Zdd11/CW**2*Zdd21*Zuu12**2+2D0*Zdd12
     &      /CW**2*Zdd22*Zuu11**2-8D0*Zdd12/CW**2*Zdd22*Zuu12
     &      **2-9D0*Zdd11/SW**2*Zdd21*Zuu11**2-Zdd11/CW**2
     &      *Zdd21*Zuu11**2*sb**2+4D0*Zdd11/CW**2*Zdd21*Zuu12
     &      **2*sb**2-2D0*Zdd12/CW**2*Zdd22*Zuu11**2*sb**2
     &      +8D0*Zdd12/CW**2*Zdd22*Zuu12**2*sb**2+9D0*Zdd11
     &      /SW**2*Zdd21*Zuu11**2*sb**2)
      AAABR(1276) = EE**2/SW**2/cb**2*(Zdd11*Zdd21*Zuu11
     &      **2-Zdd11*Zdd21*Zuu11**2*sb**2-Mq**2/MW**2*Zdd11
     &      *Zdd21*Zuu12**2+Mq**2/MW**2*Zdd12*Zdd22*Zuu11
     &      **2+Mq**2/MW**2*Zdd11*Zdd21*Zuu12**2/sb**2)
      AAABR(1277) = EE**2/cb**2*(Zcc11/CW**2*Zcc21*Zss11
     &      *Zss21+2D0*Zcc11/CW**2*Zcc21*Zss12*Zss22-4D0
     &      *Zcc12/CW**2*Zcc22*Zss11*Zss21-8D0*Zcc12/CW**2
     &      *Zcc22*Zss12*Zss22-9D0*Zcc11/SW**2*Zcc21*Zss11
     &      *Zss21-Zcc11/CW**2*Zcc21*Zss11*Zss21*sb**2-2D0
     &      *Zcc11/CW**2*Zcc21*Zss12*Zss22*sb**2+4D0*Zcc12
     &      /CW**2*Zcc22*Zss11*Zss21*sb**2+8D0*Zcc12/CW**2
     &      *Zcc22*Zss12*Zss22*sb**2+9D0*Zcc11/SW**2*Zcc21
     &      *Zss11*Zss21*sb**2)
      AAABR(1278) = EE**2/SW**2/cb**2*(Zcc11*Zcc21*Zss11
     &      *Zss21-Zcc11*Zcc21*Zss11*Zss21*sb**2-Mc**2/MW
     &      **2*Zcc12*Zcc22*Zss11*Zss21+Mq**2/MW**2*Zcc11
     &      *Zcc21*Zss12*Zss22+Mc**2/MW**2*Zcc12*Zcc22*Zss11
     &      *Zss21/sb**2)
      AAABR(1279) = EE**2*(Zdd11/CW**2*Zdd21*Zt11*Zt21-4D0
     &      *Zdd11/CW**2*Zdd21*Zt12*Zt22+2D0*Zdd12/CW**2
     &      *Zdd22*Zt11*Zt21-8D0*Zdd12/CW**2*Zdd22*Zt12*Zt22
     &      -9D0*Zdd11/SW**2*Zdd21*Zt11*Zt21)
      AAABR(1280) = EE**2*(Zb11/CW**2*Zb21*Zcc11*Zcc21-4D0
     &      *Zb11/CW**2*Zb21*Zcc12*Zcc22+2D0*Zb12/CW**2*Zb22
     &      *Zcc11*Zcc21-8D0*Zb12/CW**2*Zb22*Zcc12*Zcc22
     &      -9D0*Zb11/SW**2*Zb21*Zcc11*Zcc21)
      AAABR(1281) = EE**2/cb**2*(Zdd11/CW**2*Zdd21*Zuu11
     &      *Zuu21-4D0*Zdd11/CW**2*Zdd21*Zuu12*Zuu22+2D0
     &      *Zdd12/CW**2*Zdd22*Zuu11*Zuu21-8D0*Zdd12/CW**2
     &      *Zdd22*Zuu12*Zuu22-9D0*Zdd11/SW**2*Zdd21*Zuu11
     &      *Zuu21-Zdd11/CW**2*Zdd21*Zuu11*Zuu21*sb**2+4D0
     &      *Zdd11/CW**2*Zdd21*Zuu12*Zuu22*sb**2-2D0*Zdd12
     &      /CW**2*Zdd22*Zuu11*Zuu21*sb**2+8D0*Zdd12/CW**2
     &      *Zdd22*Zuu12*Zuu22*sb**2+9D0*Zdd11/SW**2*Zdd21
     &      *Zuu11*Zuu21*sb**2)
      AAABR(1282) = EE**2/SW**2/cb**2*(Zdd11*Zdd21*Zuu11
     &      *Zuu21-Zdd11*Zdd21*Zuu11*Zuu21*sb**2-Mq**2/MW
     &      **2*Zdd11*Zdd21*Zuu12*Zuu22+Mq**2/MW**2*Zdd12
     &      *Zdd22*Zuu11*Zuu21+Mq**2/MW**2*Zdd11*Zdd21*Zuu12
     &      *Zuu22/sb**2)
      AAABR(1283) = EE**2/SW**2/cb**2*(Zcc11*Zdd11*Zss21
     &      *Zuu21-Zcc11*Zdd11*Zss21*Zuu21*sb**2+Mq**2/MW
     &      **2*Zcc11*Zdd12*Zss22*Zuu21-Mc/MW**2*Mq*Zcc12
     &      *Zdd11*Zss21*Zuu22+Mc/MW**2*Mq*Zcc12*Zdd11*Zss21
     &      *Zuu22/sb**2)
      AAABR(1284) = EE**2/SW**2/cb**2*(Zb21*Zcc21*Zss11*Zt11
     &      -Zb21*Zcc21*Zss11*Zt11*sb**2+MbMM/MW**2*Mq*Zb22
     &      *Zcc21*Zss12*Zt11-Mc/MW**2*MtMM*Zb21*Zcc22*Zss11
     &      *Zt12+Mc/MW**2*MtMM*Zb21*Zcc22*Zss11*Zt12/sb**2)
      AAABR(1285) = EE**2*(Zb11/CW**2*Zb21*Zuu11*Zuu21-4D0
     &      *Zb11/CW**2*Zb21*Zuu12*Zuu22+2D0*Zb12/CW**2*Zb22
     &      *Zuu11*Zuu21-8D0*Zb12/CW**2*Zb22*Zuu12*Zuu22
     &      -9D0*Zb11/SW**2*Zb21*Zuu11*Zuu21)
      AAABR(1286) = EE**2/SW**2/cb**2*(Zb11*Zcc11*Zss21*Zt21
     &      -Zb11*Zcc11*Zss21*Zt21*sb**2+MbMM/MW**2*Mq*Zb12
     &      *Zcc11*Zss22*Zt21-Mc/MW**2*MtMM*Zb11*Zcc12*Zss21
     &      *Zt22+Mc/MW**2*MtMM*Zb11*Zcc12*Zss21*Zt22/sb**2)
      AAABR(1287) = EE**2/cb**2*(Zb11/CW**2*Zb21*Zt11*Zt21
     &      -4D0*Zb11/CW**2*Zb21*Zt12*Zt22+2D0*Zb12/CW**2
     &      *Zb22*Zt11*Zt21-8D0*Zb12/CW**2*Zb22*Zt12*Zt22
     &      -9D0*Zb11/SW**2*Zb21*Zt11*Zt21-Zb11/CW**2*Zb21
     &      *Zt11*Zt21*sb**2+4D0*Zb11/CW**2*Zb21*Zt12*Zt22
     &      *sb**2-2D0*Zb12/CW**2*Zb22*Zt11*Zt21*sb**2+8D0
     &      *Zb12/CW**2*Zb22*Zt12*Zt22*sb**2+9D0*Zb11/SW
     &      **2*Zb21*Zt11*Zt21*sb**2)
      AAABR(1288) = EE**2/SW**2/cb**2*(Zb11*Zb21*Zt11*Zt21
     &      -Zb11*Zb21*Zt11*Zt21*sb**2+MbMM**2/MW**2*Zb12
     &      *Zb22*Zt11*Zt21-MtMM**2/MW**2*Zb11*Zb21*Zt12
     &      *Zt22+MtMM**2/MW**2*Zb11*Zb21*Zt12*Zt22/sb**2)
      AAABR(1289) = EE**2/SW**2/cb**2*(Zb11*Zdd21*Zt21*Zuu11
     &      -Zb11*Zdd21*Zt21*Zuu11*sb**2+MbMM/MW**2*Mq*Zb12
     &      *Zdd22*Zt21*Zuu11-Mq/MW**2*MtMM*Zb11*Zdd21*Zt22
     &      *Zuu12+Mq/MW**2*MtMM*Zb11*Zdd21*Zt22*Zuu12/sb
     &      **2)
      AAABR(1290) = EE**2/SW**2/cb**2*(Zcc21*Zdd21*Zss11
     &      *Zuu11-Zcc21*Zdd21*Zss11*Zuu11*sb**2+Mq**2/MW
     &      **2*Zcc21*Zdd22*Zss12*Zuu11-Mc/MW**2*Mq*Zcc22
     &      *Zdd21*Zss11*Zuu12+Mc/MW**2*Mq*Zcc22*Zdd21*Zss11
     &      *Zuu12/sb**2)
      AAABR(1291) = EE**2*(Zcc11/CW**2*Zcc21*Zdd11*Zdd21
     &      +2D0*Zcc11/CW**2*Zcc21*Zdd12*Zdd22-4D0*Zcc12
     &      /CW**2*Zcc22*Zdd11*Zdd21-8D0*Zcc12/CW**2*Zcc22
     &      *Zdd12*Zdd22-9D0*Zcc11/SW**2*Zcc21*Zdd11*Zdd21)
      AAABR(1292) = EE**2/SW**2/cb**2*(Zb21*Zdd11*Zt11*Zuu21
     &      -Zb21*Zdd11*Zt11*Zuu21*sb**2+MbMM/MW**2*Mq*Zb22
     &      *Zdd12*Zt11*Zuu21-Mq/MW**2*MtMM*Zb21*Zdd11*Zt12
     &      *Zuu22+Mq/MW**2*MtMM*Zb21*Zdd11*Zt12*Zuu22/sb
     &      **2)
      AAABR(1293) = EE**2*(Zss11/CW**2*Zss21*Zt11*Zt21-4D0
     &      *Zss11/CW**2*Zss21*Zt12*Zt22+2D0*Zss12/CW**2
     &      *Zss22*Zt11*Zt21-8D0*Zss12/CW**2*Zss22*Zt12*Zt22
     &      -9D0*Zss11/SW**2*Zss21*Zt11*Zt21)
      AAABR(1294) = EE**2*(Zss11/CW**2*Zss21*Zuu11*Zuu21
     &      -4D0*Zss11/CW**2*Zss21*Zuu12*Zuu22+2D0*Zss12
     &      /CW**2*Zss22*Zuu11*Zuu21-8D0*Zss12/CW**2*Zss22
     &      *Zuu12*Zuu22-9D0*Zss11/SW**2*Zss21*Zuu11*Zuu21)
      AAABR(1295) = EE**2/SW**2/cb**2*(Zb11*Zdd11*Zt21*Zuu21
     &      -Zb11*Zdd11*Zt21*Zuu21*sb**2+MbMM/MW**2*Mq*Zb12
     &      *Zdd12*Zt21*Zuu21-Mq/MW**2*MtMM*Zb11*Zdd11*Zt22
     &      *Zuu22+Mq/MW**2*MtMM*Zb11*Zdd11*Zt22*Zuu22/sb
     &      **2)
      AAABR(1296) = EE**2/cb**2*(Zcc21**2/CW**2*Zss11**2
     &      +2D0*Zcc21**2/CW**2*Zss12**2-4D0*Zcc22**2/CW
     &      **2*Zss11**2-8D0*Zcc22**2/CW**2*Zss12**2-9D0
     &      *Zcc21**2/SW**2*Zss11**2-Zcc21**2/CW**2*Zss11
     &      **2*sb**2-2D0*Zcc21**2/CW**2*Zss12**2*sb**2+4D0
     &      *Zcc22**2/CW**2*Zss11**2*sb**2+8D0*Zcc22**2/CW
     &      **2*Zss12**2*sb**2+9D0*Zcc21**2/SW**2*Zss11**2
     &      *sb**2)
      AAABR(1297) = EE**2/SW**2/cb**2*(Zcc21**2*Zss11**2
     &      -Zcc21**2*Zss11**2*sb**2-Mc**2/MW**2*Zcc22**2
     &      *Zss11**2+Mq**2/MW**2*Zcc21**2*Zss12**2+Mc**2
     &      /MW**2*Zcc22**2*Zss11**2/sb**2)
      AAABR(1298) = EE**2*(Zcc21**2/CW**2*Zdd11**2+2D0*Zcc21
     &      **2/CW**2*Zdd12**2-4D0*Zcc22**2/CW**2*Zdd11**2
     &      -8D0*Zcc22**2/CW**2*Zdd12**2-9D0*Zcc21**2/SW
     &      **2*Zdd11**2)
      AAABR(1299) = EE**2*(Zss11**2/CW**2*Zt21**2-4D0*Zss11
     &      **2/CW**2*Zt22**2+2D0*Zss12**2/CW**2*Zt21**2
     &      -8D0*Zss12**2/CW**2*Zt22**2-9D0*Zss11**2/SW**2
     &      *Zt21**2)
      AAABR(1300) = EE**2*(Zb11**2/CW**2*Zuu21**2-4D0*Zb11
     &      **2/CW**2*Zuu22**2+2D0*Zb12**2/CW**2*Zuu21**2
     &      -8D0*Zb12**2/CW**2*Zuu22**2-9D0*Zb11**2/SW**2
     &      *Zuu21**2)
      AAABR(1301) = EE**2/SW**2/cb**2*(Zcc21*Zdd11*Zss11
     &      *Zuu21-Zcc21*Zdd11*Zss11*Zuu21*sb**2+Mq**2/MW
     &      **2*Zcc21*Zdd12*Zss12*Zuu21-Mc/MW**2*Mq*Zcc22
     &      *Zdd11*Zss11*Zuu22+Mc/MW**2*Mq*Zcc22*Zdd11*Zss11
     &      *Zuu22/sb**2)
      AAABR(1302) = EE**2/SW**2/cb**2*(Zb11*Zcc21*Zss11*Zt21
     &      -Zb11*Zcc21*Zss11*Zt21*sb**2+MbMM/MW**2*Mq*Zb12
     &      *Zcc21*Zss12*Zt21-Mc/MW**2*MtMM*Zb11*Zcc22*Zss11
     &      *Zt22+Mc/MW**2*MtMM*Zb11*Zcc22*Zss11*Zt22/sb**2)
      AAABR(1303) = EE**2*(Zss11**2/CW**2*Zuu21**2-4D0*Zss11
     &      **2/CW**2*Zuu22**2+2D0*Zss12**2/CW**2*Zuu21**2
     &      -8D0*Zss12**2/CW**2*Zuu22**2-9D0*Zss11**2/SW
     &      **2*Zuu21**2)
      AAABR(1304) = EE**2*(Zb11**2/CW**2*Zcc21**2-4D0*Zb11
     &      **2/CW**2*Zcc22**2+2D0*Zb12**2/CW**2*Zcc21**2
     &      -8D0*Zb12**2/CW**2*Zcc22**2-9D0*Zb11**2/SW**2
     &      *Zcc21**2)
      AAABR(1305) = EE**2/cb**2*(Zb11**2/CW**2*Zt21**2-4D0
     &      *Zb11**2/CW**2*Zt22**2+2D0*Zb12**2/CW**2*Zt21
     &      **2-8D0*Zb12**2/CW**2*Zt22**2-9D0*Zb11**2/SW
     &      **2*Zt21**2-Zb11**2/CW**2*Zt21**2*sb**2+4D0*Zb11
     &      **2/CW**2*Zt22**2*sb**2-2D0*Zb12**2/CW**2*Zt21
     &      **2*sb**2+8D0*Zb12**2/CW**2*Zt22**2*sb**2+9D0
     &      *Zb11**2/SW**2*Zt21**2*sb**2)
      AAABR(1306) = EE**2/SW**2/cb**2*(Zb11**2*Zt21**2-Zb11
     &      **2*Zt21**2*sb**2+MbMM**2/MW**2*Zb12**2*Zt21
     &      **2-MtMM**2/MW**2*Zb11**2*Zt22**2+MtMM**2/MW
     &      **2*Zb11**2*Zt22**2/sb**2)
      AAABR(1307) = EE**2/cb**2*(Zdd11**2/CW**2*Zuu21**2
     &      -4D0*Zdd11**2/CW**2*Zuu22**2+2D0*Zdd12**2/CW
     &      **2*Zuu21**2-8D0*Zdd12**2/CW**2*Zuu22**2-9D0
     &      *Zdd11**2/SW**2*Zuu21**2-Zdd11**2/CW**2*Zuu21
     &      **2*sb**2+4D0*Zdd11**2/CW**2*Zuu22**2*sb**2-2D0
     &      *Zdd12**2/CW**2*Zuu21**2*sb**2+8D0*Zdd12**2/CW
     &      **2*Zuu22**2*sb**2+9D0*Zdd11**2/SW**2*Zuu21**2
     &      *sb**2)
      AAABR(1308) = EE**2/SW**2/cb**2*(Zdd11**2*Zuu21**2
     &      -Zdd11**2*Zuu21**2*sb**2-Mq**2/MW**2*Zdd11**2
     &      *Zuu22**2+Mq**2/MW**2*Zdd12**2*Zuu21**2+Mq**2
     &      /MW**2*Zdd11**2*Zuu22**2/sb**2)
      AAABR(1309) = EE**2*(Zdd11**2/CW**2*Zt21**2-4D0*Zdd11
     &      **2/CW**2*Zt22**2+2D0*Zdd12**2/CW**2*Zt21**2
     &      -8D0*Zdd12**2/CW**2*Zt22**2-9D0*Zdd11**2/SW**2
     &      *Zt21**2)
      AAABR(1310) = EE**2/SW**2/cb**2*(Zb21*Zcc11*Zss11*Zt21
     &      -Zb21*Zcc11*Zss11*Zt21*sb**2+MbMM/MW**2*Mq*Zb22
     &      *Zcc11*Zss12*Zt21-Mc/MW**2*MtMM*Zb21*Zcc12*Zss11
     &      *Zt22+Mc/MW**2*MtMM*Zb21*Zcc12*Zss11*Zt22/sb**2)
      AAABR(1311) = EE**2/SW**2/cb**2*(Zb11*Zcc21*Zss21*Zt11
     &      -Zb11*Zcc21*Zss21*Zt11*sb**2+MbMM/MW**2*Mq*Zb12
     &      *Zcc21*Zss22*Zt11-Mc/MW**2*MtMM*Zb11*Zcc22*Zss21
     &      *Zt12+Mc/MW**2*MtMM*Zb11*Zcc22*Zss21*Zt12/sb**2)
      AAABR(1312) = EE**2/SW**2/cb**2*(Zcc21*Zdd11*Zss21
     &      *Zuu11-Zcc21*Zdd11*Zss21*Zuu11*sb**2+Mq**2/MW
     &      **2*Zcc21*Zdd12*Zss22*Zuu11-Mc/MW**2*Mq*Zcc22
     &      *Zdd11*Zss21*Zuu12+Mc/MW**2*Mq*Zcc22*Zdd11*Zss21
     &      *Zuu12/sb**2)
      AAABR(1313) = EE**2/SW**2/cb**2*(Zb11*Zdd21*Zt11*Zuu21
     &      -Zb11*Zdd21*Zt11*Zuu21*sb**2+MbMM/MW**2*Mq*Zb12
     &      *Zdd22*Zt11*Zuu21-Mq/MW**2*MtMM*Zb11*Zdd21*Zt12
     &      *Zuu22+Mq/MW**2*MtMM*Zb11*Zdd21*Zt12*Zuu22/sb
     &      **2)
      AAABR(1314) = EE**2/SW**2/cb**2*(Zcc11*Zdd21*Zss11
     &      *Zuu21-Zcc11*Zdd21*Zss11*Zuu21*sb**2+Mq**2/MW
     &      **2*Zcc11*Zdd22*Zss12*Zuu21-Mc/MW**2*Mq*Zcc12
     &      *Zdd21*Zss11*Zuu22+Mc/MW**2*Mq*Zcc12*Zdd21*Zss11
     &      *Zuu22/sb**2)
      AAABR(1315) = EE**2/SW**2/cb**2*(Zb21*Zdd11*Zt21*Zuu11
     &      -Zb21*Zdd11*Zt21*Zuu11*sb**2+MbMM/MW**2*Mq*Zb22
     &      *Zdd12*Zt21*Zuu11-Mq/MW**2*MtMM*Zb21*Zdd11*Zt22
     &      *Zuu12+Mq/MW**2*MtMM*Zb21*Zdd11*Zt22*Zuu12/sb
     &      **2)
      AAABR(1316) = EE**2/cb**2*(Zcc21**2/CW**2*Zss11*Zss21
     &      +2D0*Zcc21**2/CW**2*Zss12*Zss22-4D0*Zcc22**2
     &      /CW**2*Zss11*Zss21-8D0*Zcc22**2/CW**2*Zss12*Zss22
     &      -9D0*Zcc21**2/SW**2*Zss11*Zss21-Zcc21**2/CW**2
     &      *Zss11*Zss21*sb**2-2D0*Zcc21**2/CW**2*Zss12*Zss22
     &      *sb**2+4D0*Zcc22**2/CW**2*Zss11*Zss21*sb**2+8D0
     &      *Zcc22**2/CW**2*Zss12*Zss22*sb**2+9D0*Zcc21**2
     &      /SW**2*Zss11*Zss21*sb**2)
      AAABR(1317) = EE**2/SW**2/cb**2*(Zcc21**2*Zss11*Zss21
     &      -Zcc21**2*Zss11*Zss21*sb**2-Mc**2/MW**2*Zcc22
     &      **2*Zss11*Zss21+Mq**2/MW**2*Zcc21**2*Zss12*Zss22
     &      +Mc**2/MW**2*Zcc22**2*Zss11*Zss21/sb**2)
      AAABR(1318) = EE**2/cb**2*(Zb11/CW**2*Zb21*Zt21**2
     &      -4D0*Zb11/CW**2*Zb21*Zt22**2+2D0*Zb12/CW**2*Zb22
     &      *Zt21**2-8D0*Zb12/CW**2*Zb22*Zt22**2-9D0*Zb11
     &      /SW**2*Zb21*Zt21**2-Zb11/CW**2*Zb21*Zt21**2*sb
     &      **2+4D0*Zb11/CW**2*Zb21*Zt22**2*sb**2-2D0*Zb12
     &      /CW**2*Zb22*Zt21**2*sb**2+8D0*Zb12/CW**2*Zb22
     &      *Zt22**2*sb**2+9D0*Zb11/SW**2*Zb21*Zt21**2*sb
     &      **2)
      AAABR(1319) = EE**2/SW**2/cb**2*(Zb11*Zb21*Zt21**2
     &      -Zb11*Zb21*Zt21**2*sb**2+MbMM**2/MW**2*Zb12*Zb22
     &      *Zt21**2-MtMM**2/MW**2*Zb11*Zb21*Zt22**2+MtMM
     &      **2/MW**2*Zb11*Zb21*Zt22**2/sb**2)
      AAABR(1320) = EE**2/SW**2/cb**2*(Zcc21*Zdd11*Zss21
     &      *Zuu21-Zcc21*Zdd11*Zss21*Zuu21*sb**2+Mq**2/MW
     &      **2*Zcc21*Zdd12*Zss22*Zuu21-Mc/MW**2*Mq*Zcc22
     &      *Zdd11*Zss21*Zuu22+Mc/MW**2*Mq*Zcc22*Zdd11*Zss21
     &      *Zuu22/sb**2)
      AAABR(1321) = EE**2/SW**2/cb**2*(Zcc21*Zdd21*Zss11
     &      *Zuu21-Zcc21*Zdd21*Zss11*Zuu21*sb**2+Mq**2/MW
     &      **2*Zcc21*Zdd22*Zss12*Zuu21-Mc/MW**2*Mq*Zcc22
     &      *Zdd21*Zss11*Zuu22+Mc/MW**2*Mq*Zcc22*Zdd21*Zss11
     &      *Zuu22/sb**2)
      AAABR(1322) = EE**2/SW**2/cb**2*(Zb11*Zdd21*Zt21*Zuu21
     &      -Zb11*Zdd21*Zt21*Zuu21*sb**2+MbMM/MW**2*Mq*Zb12
     &      *Zdd22*Zt21*Zuu21-Mq/MW**2*MtMM*Zb11*Zdd21*Zt22
     &      *Zuu22+Mq/MW**2*MtMM*Zb11*Zdd21*Zt22*Zuu22/sb
     &      **2)
      AAABR(1323) = EE**2/SW**2/cb**2*(Zb21*Zcc21*Zss11*Zt21
     &      -Zb21*Zcc21*Zss11*Zt21*sb**2+MbMM/MW**2*Mq*Zb22
     &      *Zcc21*Zss12*Zt21-Mc/MW**2*MtMM*Zb21*Zcc22*Zss11
     &      *Zt22+Mc/MW**2*MtMM*Zb21*Zcc22*Zss11*Zt22/sb**2)
      AAABR(1324) = EE**2*(Zb11/CW**2*Zb21*Zuu21**2-4D0*Zb11
     &      /CW**2*Zb21*Zuu22**2+2D0*Zb12/CW**2*Zb22*Zuu21
     &      **2-8D0*Zb12/CW**2*Zb22*Zuu22**2-9D0*Zb11/SW
     &      **2*Zb21*Zuu21**2)
      AAABR(1325) = EE**2*(Zb11/CW**2*Zb21*Zcc21**2-4D0*Zb11
     &      /CW**2*Zb21*Zcc22**2+2D0*Zb12/CW**2*Zb22*Zcc21
     &      **2-8D0*Zb12/CW**2*Zb22*Zcc22**2-9D0*Zb11/SW
     &      **2*Zb21*Zcc21**2)
      AAABR(1326) = EE**2/SW**2/cb**2*(Zb11*Zcc21*Zss21*Zt21
     &      -Zb11*Zcc21*Zss21*Zt21*sb**2+MbMM/MW**2*Mq*Zb12
     &      *Zcc21*Zss22*Zt21-Mc/MW**2*MtMM*Zb11*Zcc22*Zss21
     &      *Zt22+Mc/MW**2*MtMM*Zb11*Zcc22*Zss21*Zt22/sb**2)
      AAABR(1327) = EE**2*(Zss11/CW**2*Zss21*Zt21**2-4D0
     &      *Zss11/CW**2*Zss21*Zt22**2+2D0*Zss12/CW**2*Zss22
     &      *Zt21**2-8D0*Zss12/CW**2*Zss22*Zt22**2-9D0*Zss11
     &      /SW**2*Zss21*Zt21**2)
      AAABR(1328) = EE**2*(Zcc21**2/CW**2*Zdd11*Zdd21+2D0
     &      *Zcc21**2/CW**2*Zdd12*Zdd22-4D0*Zcc22**2/CW**2
     &      *Zdd11*Zdd21-8D0*Zcc22**2/CW**2*Zdd12*Zdd22-9D0
     &      *Zcc21**2/SW**2*Zdd11*Zdd21)
      AAABR(1329) = EE**2/SW**2/cb**2*(Zb21*Zdd11*Zt21*Zuu21
     &      -Zb21*Zdd11*Zt21*Zuu21*sb**2+MbMM/MW**2*Mq*Zb22
     &      *Zdd12*Zt21*Zuu21-Mq/MW**2*MtMM*Zb21*Zdd11*Zt22
     &      *Zuu22+Mq/MW**2*MtMM*Zb21*Zdd11*Zt22*Zuu22/sb
     &      **2)
      AAABR(1330) = EE**2*(Zdd11/CW**2*Zdd21*Zt21**2-4D0
     &      *Zdd11/CW**2*Zdd21*Zt22**2+2D0*Zdd12/CW**2*Zdd22
     &      *Zt21**2-8D0*Zdd12/CW**2*Zdd22*Zt22**2-9D0*Zdd11
     &      /SW**2*Zdd21*Zt21**2)
      AAABR(1331) = EE**2*(Zss11/CW**2*Zss21*Zuu21**2-4D0
     &      *Zss11/CW**2*Zss21*Zuu22**2+2D0*Zss12/CW**2*Zss22
     &      *Zuu21**2-8D0*Zss12/CW**2*Zss22*Zuu22**2-9D0
     &      *Zss11/SW**2*Zss21*Zuu21**2)
      AAABR(1332) = EE**2/cb**2*(Zdd11/CW**2*Zdd21*Zuu21
     &      **2-4D0*Zdd11/CW**2*Zdd21*Zuu22**2+2D0*Zdd12
     &      /CW**2*Zdd22*Zuu21**2-8D0*Zdd12/CW**2*Zdd22*Zuu22
     &      **2-9D0*Zdd11/SW**2*Zdd21*Zuu21**2-Zdd11/CW**2
     &      *Zdd21*Zuu21**2*sb**2+4D0*Zdd11/CW**2*Zdd21*Zuu22
     &      **2*sb**2-2D0*Zdd12/CW**2*Zdd22*Zuu21**2*sb**2
     &      +8D0*Zdd12/CW**2*Zdd22*Zuu22**2*sb**2+9D0*Zdd11
     &      /SW**2*Zdd21*Zuu21**2*sb**2)
      AAABR(1333) = EE**2/SW**2/cb**2*(Zdd11*Zdd21*Zuu21
     &      **2-Zdd11*Zdd21*Zuu21**2*sb**2-Mq**2/MW**2*Zdd11
     &      *Zdd21*Zuu22**2+Mq**2/MW**2*Zdd12*Zdd22*Zuu21
     &      **2+Mq**2/MW**2*Zdd11*Zdd21*Zuu22**2/sb**2)
      AAABR(1334) = EE**2*(Zdd11**2/CW**2+2D0*Zdd12**2/CW
     &      **2-3D0*Zdd11**2/SW**2-2D0*Zdd11**2/CW**2*sb
     &      **2-4D0*Zdd12**2/CW**2*sb**2-6D0*Mq**2/MW**2
     &      /SW**2+6D0*Zdd11**2/SW**2*sb**2)
      AAABR(1335) = EE**2*(Zss11**2/CW**2+2D0*Zss12**2/CW
     &      **2-3D0*Zss11**2/SW**2-2D0*Zss11**2/CW**2*sb
     &      **2-4D0*Zss12**2/CW**2*sb**2+6D0*Zss11**2/SW
     &      **2*sb**2-6D0*Mc**2/MW**2/SW**2*Zss11**2-6D0
     &      *Mq**2/MW**2/SW**2*Zss12**2)
      AAABR(1336) = EE**2*(Zb11**2/CW**2+2D0*Zb12**2/CW**2
     &      -3D0*Zb11**2/SW**2-2D0*Zb11**2/CW**2*sb**2-4D0
     &      *Zb12**2/CW**2*sb**2+6D0*Zb11**2/SW**2*sb**2
     &      -6D0*MbMM**2/MW**2/SW**2*Zb12**2-6D0*MtMM**2
     &      /MW**2/SW**2*Zb11**2)
      AAABR(1337) = EE**2*(Zss11/CW**2*Zss21+2D0*Zss12/CW
     &      **2*Zss22-3D0*Zss11/SW**2*Zss21-2D0*Zss11/CW
     &      **2*Zss21*sb**2-4D0*Zss12/CW**2*Zss22*sb**2+6D0
     &      *Zss11/SW**2*Zss21*sb**2-6D0*Mc**2/MW**2/SW**2
     &      *Zss11*Zss21-6D0*Mq**2/MW**2/SW**2*Zss12*Zss22)
      AAABR(1338) = EE**2*(Zb11/CW**2*Zb21+2D0*Zb12/CW**2
     &      *Zb22-3D0*Zb11/SW**2*Zb21-2D0*Zb11/CW**2*Zb21
     &      *sb**2-4D0*Zb12/CW**2*Zb22*sb**2+6D0*Zb11/SW
     &      **2*Zb21*sb**2-6D0*MbMM**2/MW**2/SW**2*Zb12*Zb22
     &      -6D0*MtMM**2/MW**2/SW**2*Zb11*Zb21)
      AAABR(1339) = EE**2*(Zdd11/CW**2*Zdd21+2D0*Zdd12/CW
     &      **2*Zdd22-3D0*Zdd11/SW**2*Zdd21-2D0*Zdd11/CW
     &      **2*Zdd21*sb**2-4D0*Zdd12/CW**2*Zdd22*sb**2+6D0
     &      *Zdd11/SW**2*Zdd21*sb**2)
      AAABR(1340) = EE**2/SW**2*Sqrt2*Zb11*Zt11*(1D0-2D0
     &      *sb**2-MbMM**2/MW**2+MtMM**2/MW**2)
      AAABR(1341) = EE**2/SW**2*Sqrt2*Zdd11*Zuu11*(1D0-2D0
     &      *sb**2)
      AAABR(1342) = EE**2/SW**2*Sqrt2*Zcc11*Zss11*(1D0-2D0
     &      *sb**2+Mc**2/MW**2-Mq**2/MW**2)
      AAABR(1343) = EE**2/SW**2*Sqrt2*Zdd11*Zuu21*(1D0-2D0
     &      *sb**2)
      AAABR(1344) = EE**2/SW**2*Sqrt2*Zcc21*Zss11*(1D0-2D0
     &      *sb**2+Mc**2/MW**2-Mq**2/MW**2)
      AAABR(1345) = EE**2/SW**2*Sqrt2*Zb11*Zt21*(1D0-2D0
     &      *sb**2-MbMM**2/MW**2+MtMM**2/MW**2)
      AAABR(1346) = 1D0/cb*(2D0*EE**2*Pa21/SW**2*Sqrt2*Zcc11
     &      *Zss11*sb-2D0*EE**2*Pa21/SW**2*Sqrt2*Zcc11*Zss11
     &      *sb**3+2D0*EE*Lambda/MW*Mc*Pa22/SW*Zcc12*Zss11
     &      /sb-2D0*EE*Lambda/MW*Mc*Pa22/SW*Zcc12*Zss11*sb
     &      +2D0*EE*Lambda/MW*Mq*Pa22/SW*Zcc11*Zss12*sb-EE
     &      **2/MW**2*Mc**2*Pa21/SW**2*Sqrt2*Zcc11*Zss11
     &      /sb+EE**2/MW**2*Mc**2*Pa21/SW**2*Sqrt2*Zcc11
     &      *Zss11*sb-EE**2/MW**2*Mq**2*Pa21/SW**2*Sqrt2
     &      *Zcc11*Zss11*sb-EE**2/MW**2*Mc*Mq*Pa21/SW**2
     &      *Sqrt2*Zcc12*Zss12/sb)
      AAABR(1347) = 1D0/cb*(2D0*EE**2*Pa11/SW**2*Sqrt2*Zb11
     &      *Zt11*sb-2D0*EE**2*Pa11/SW**2*Sqrt2*Zb11*Zt11
     &      *sb**3+2D0*EE*Lambda/MW*MbMM*Pa12/SW*Zb12*Zt11
     &      *sb+2D0*EE*Lambda/MW*MtMM*Pa12/SW*Zb11*Zt12/sb
     &      -2D0*EE*Lambda/MW*MtMM*Pa12/SW*Zb11*Zt12*sb-EE
     &      **2/MW**2*MbMM**2*Pa11/SW**2*Sqrt2*Zb11*Zt11
     &      *sb-EE**2/MW**2*MtMM**2*Pa11/SW**2*Sqrt2*Zb11
     &      *Zt11/sb+EE**2/MW**2*MtMM**2*Pa11/SW**2*Sqrt2
     &      *Zb11*Zt11*sb-EE**2/MW**2*MbMM*MtMM*Pa11/SW**2
     &      *Sqrt2*Zb12*Zt12/sb)
      AAABR(1348) = 1D0/cb*(2D0*EE**2*Pa21/SW**2*Sqrt2*Zdd11
     &      *Zuu11*sb-2D0*EE**2*Pa21/SW**2*Sqrt2*Zdd11*Zuu11
     &      *sb**3+2D0*EE*Lambda/MW*Mq*Pa22/SW*Zdd11*Zuu12
     &      /sb-2D0*EE*Lambda/MW*Mq*Pa22/SW*Zdd11*Zuu12*sb
     &      +2D0*EE*Lambda/MW*Mq*Pa22/SW*Zdd12*Zuu11*sb-EE
     &      **2/MW**2*Mq**2*Pa21/SW**2*Sqrt2*Zdd11*Zuu11
     &      /sb-EE**2/MW**2*Mq**2*Pa21/SW**2*Sqrt2*Zdd12
     &      *Zuu12/sb)
      AAABR(1349) = 1D0/cb*(2D0*EE**2*Pa21/SW**2*Sqrt2*Zb11
     &      *Zt11*sb-2D0*EE**2*Pa21/SW**2*Sqrt2*Zb11*Zt11
     &      *sb**3+2D0*EE*Lambda/MW*MbMM*Pa22/SW*Zb12*Zt11
     &      *sb+2D0*EE*Lambda/MW*MtMM*Pa22/SW*Zb11*Zt12/sb
     &      -2D0*EE*Lambda/MW*MtMM*Pa22/SW*Zb11*Zt12*sb-EE
     &      **2/MW**2*MbMM**2*Pa21/SW**2*Sqrt2*Zb11*Zt11
     &      *sb-EE**2/MW**2*MtMM**2*Pa21/SW**2*Sqrt2*Zb11
     &      *Zt11/sb+EE**2/MW**2*MtMM**2*Pa21/SW**2*Sqrt2
     &      *Zb11*Zt11*sb-EE**2/MW**2*MbMM*MtMM*Pa21/SW**2
     &      *Sqrt2*Zb12*Zt12/sb)
      AAABR(1350) = 1D0/cb*(2D0*EE**2*Pa11/SW**2*Sqrt2*Zdd11
     &      *Zuu11*sb-2D0*EE**2*Pa11/SW**2*Sqrt2*Zdd11*Zuu11
     &      *sb**3+2D0*EE*Lambda/MW*Mq*Pa12/SW*Zdd11*Zuu12
     &      /sb-2D0*EE*Lambda/MW*Mq*Pa12/SW*Zdd11*Zuu12*sb
     &      +2D0*EE*Lambda/MW*Mq*Pa12/SW*Zdd12*Zuu11*sb-EE
     &      **2/MW**2*Mq**2*Pa11/SW**2*Sqrt2*Zdd11*Zuu11
     &      /sb-EE**2/MW**2*Mq**2*Pa11/SW**2*Sqrt2*Zdd12
     &      *Zuu12/sb)
      AAABR(1351) = 1D0/cb*(2D0*EE**2*Pa11/SW**2*Sqrt2*Zcc11
     &      *Zss11*sb-2D0*EE**2*Pa11/SW**2*Sqrt2*Zcc11*Zss11
     &      *sb**3+2D0*EE*Lambda/MW*Mc*Pa12/SW*Zcc12*Zss11
     &      /sb-2D0*EE*Lambda/MW*Mc*Pa12/SW*Zcc12*Zss11*sb
     &      +2D0*EE*Lambda/MW*Mq*Pa12/SW*Zcc11*Zss12*sb-EE
     &      **2/MW**2*Mc**2*Pa11/SW**2*Sqrt2*Zcc11*Zss11
     &      /sb+EE**2/MW**2*Mc**2*Pa11/SW**2*Sqrt2*Zcc11
     &      *Zss11*sb-EE**2/MW**2*Mq**2*Pa11/SW**2*Sqrt2
     &      *Zcc11*Zss11*sb-EE**2/MW**2*Mc*Mq*Pa11/SW**2
     &      *Sqrt2*Zcc12*Zss12/sb)
      AAABR(1352) = 1D0/cb*(2D0*EE**2*Pa11/SW**2*Sqrt2*Zcc21
     &      *Zss11*sb-2D0*EE**2*Pa11/SW**2*Sqrt2*Zcc21*Zss11
     &      *sb**3+2D0*EE*Lambda/MW*Mc*Pa12/SW*Zcc22*Zss11
     &      /sb-2D0*EE*Lambda/MW*Mc*Pa12/SW*Zcc22*Zss11*sb
     &      +2D0*EE*Lambda/MW*Mq*Pa12/SW*Zcc21*Zss12*sb-EE
     &      **2/MW**2*Mc**2*Pa11/SW**2*Sqrt2*Zcc21*Zss11
     &      /sb+EE**2/MW**2*Mc**2*Pa11/SW**2*Sqrt2*Zcc21
     &      *Zss11*sb-EE**2/MW**2*Mq**2*Pa11/SW**2*Sqrt2
     &      *Zcc21*Zss11*sb-EE**2/MW**2*Mc*Mq*Pa11/SW**2
     &      *Sqrt2*Zcc22*Zss12/sb)
      AAABR(1353) = 1D0/cb*(2D0*EE**2*Pa21/SW**2*Sqrt2*Zb11
     &      *Zt21*sb-2D0*EE**2*Pa21/SW**2*Sqrt2*Zb11*Zt21
     &      *sb**3+2D0*EE*Lambda/MW*MbMM*Pa22/SW*Zb12*Zt21
     &      *sb+2D0*EE*Lambda/MW*MtMM*Pa22/SW*Zb11*Zt22/sb
     &      -2D0*EE*Lambda/MW*MtMM*Pa22/SW*Zb11*Zt22*sb-EE
     &      **2/MW**2*MbMM**2*Pa21/SW**2*Sqrt2*Zb11*Zt21
     &      *sb-EE**2/MW**2*MtMM**2*Pa21/SW**2*Sqrt2*Zb11
     &      *Zt21/sb+EE**2/MW**2*MtMM**2*Pa21/SW**2*Sqrt2
     &      *Zb11*Zt21*sb-EE**2/MW**2*MbMM*MtMM*Pa21/SW**2
     &      *Sqrt2*Zb12*Zt22/sb)
      AAABR(1354) = 1D0/cb*(2D0*EE**2*Pa21/SW**2*Sqrt2*Zdd11
     &      *Zuu21*sb-2D0*EE**2*Pa21/SW**2*Sqrt2*Zdd11*Zuu21
     &      *sb**3+2D0*EE*Lambda/MW*Mq*Pa22/SW*Zdd11*Zuu22
     &      /sb-2D0*EE*Lambda/MW*Mq*Pa22/SW*Zdd11*Zuu22*sb
     &      +2D0*EE*Lambda/MW*Mq*Pa22/SW*Zdd12*Zuu21*sb-EE
     &      **2/MW**2*Mq**2*Pa21/SW**2*Sqrt2*Zdd11*Zuu21
     &      /sb-EE**2/MW**2*Mq**2*Pa21/SW**2*Sqrt2*Zdd12
     &      *Zuu22/sb)
      AAABR(1355) = 1D0/cb*(2D0*EE**2*Pa11/SW**2*Sqrt2*Zb11
     &      *Zt21*sb-2D0*EE**2*Pa11/SW**2*Sqrt2*Zb11*Zt21
     &      *sb**3+2D0*EE*Lambda/MW*MbMM*Pa12/SW*Zb12*Zt21
     &      *sb+2D0*EE*Lambda/MW*MtMM*Pa12/SW*Zb11*Zt22/sb
     &      -2D0*EE*Lambda/MW*MtMM*Pa12/SW*Zb11*Zt22*sb-EE
     &      **2/MW**2*MbMM**2*Pa11/SW**2*Sqrt2*Zb11*Zt21
     &      *sb-EE**2/MW**2*MtMM**2*Pa11/SW**2*Sqrt2*Zb11
     &      *Zt21/sb+EE**2/MW**2*MtMM**2*Pa11/SW**2*Sqrt2
     &      *Zb11*Zt21*sb-EE**2/MW**2*MbMM*MtMM*Pa11/SW**2
     &      *Sqrt2*Zb12*Zt22/sb)
      AAABR(1356) = 1D0/cb*(2D0*EE**2*Pa11/SW**2*Sqrt2*Zdd11
     &      *Zuu21*sb-2D0*EE**2*Pa11/SW**2*Sqrt2*Zdd11*Zuu21
     &      *sb**3+2D0*EE*Lambda/MW*Mq*Pa12/SW*Zdd11*Zuu22
     &      /sb-2D0*EE*Lambda/MW*Mq*Pa12/SW*Zdd11*Zuu22*sb
     &      +2D0*EE*Lambda/MW*Mq*Pa12/SW*Zdd12*Zuu21*sb-EE
     &      **2/MW**2*Mq**2*Pa11/SW**2*Sqrt2*Zdd11*Zuu21
     &      /sb-EE**2/MW**2*Mq**2*Pa11/SW**2*Sqrt2*Zdd12
     &      *Zuu22/sb)
      AAABR(1357) = 1D0/cb*(2D0*EE**2*Pa21/SW**2*Sqrt2*Zcc21
     &      *Zss11*sb-2D0*EE**2*Pa21/SW**2*Sqrt2*Zcc21*Zss11
     &      *sb**3+2D0*EE*Lambda/MW*Mc*Pa22/SW*Zcc22*Zss11
     &      /sb-2D0*EE*Lambda/MW*Mc*Pa22/SW*Zcc22*Zss11*sb
     &      +2D0*EE*Lambda/MW*Mq*Pa22/SW*Zcc21*Zss12*sb-EE
     &      **2/MW**2*Mc**2*Pa21/SW**2*Sqrt2*Zcc21*Zss11
     &      /sb+EE**2/MW**2*Mc**2*Pa21/SW**2*Sqrt2*Zcc21
     &      *Zss11*sb-EE**2/MW**2*Mq**2*Pa21/SW**2*Sqrt2
     &      *Zcc21*Zss11*sb-EE**2/MW**2*Mc*Mq*Pa21/SW**2
     &      *Sqrt2*Zcc22*Zss12/sb)
      AAABR(1358) = EE**2/SW**2*Sqrt2*Zdd11*Zh11*Zuu11/cb
     &      -EE**2/SW**2*Sqrt2*Zdd11*Zh12*Zuu11*sb-EE**2
     &      /SW**2*Sqrt2*Zdd11*Zh11*Zuu11/cb*sb**2-EE**2
     &      /MW**2*Mq**2/SW**2*Sqrt2*Zdd11*Zh11*Zuu11/cb
     &      +EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd11*Zh12*Zuu11
     &      /sb+EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd12*Zh11
     &      *Zuu12/cb-EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd12
     &      *Zh12*Zuu12/sb-2D0*EE*Lambda/MW*Mq/SW*Zdd11*Zh13
     &      *Zuu12/cb/sb+2D0*EE*Lambda/MW*Mq/SW*Zdd11*Zh13
     &      *Zuu12/cb*sb+2D0*EE*Lambda/MW*Mq/SW*Zdd12*Zh13
     &      *Zuu11/cb*sb
      AAABR(1359) = EE**2/SW**2*Sqrt2*Zcc11*Zh31*Zss11/cb
     &      -EE**2/SW**2*Sqrt2*Zcc11*Zh32*Zss11*sb-EE**2
     &      /SW**2*Sqrt2*Zcc11*Zh31*Zss11/cb*sb**2+EE**2
     &      /MW**2*Mc**2/SW**2*Sqrt2*Zcc11*Zh32*Zss11/sb
     &      -EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zcc11*Zh31*Zss11
     &      /cb-2D0*EE*Lambda/MW*Mc/SW*Zcc12*Zh33*Zss11/cb
     &      /sb+2D0*EE*Lambda/MW*Mc/SW*Zcc12*Zh33*Zss11/cb
     &      *sb+2D0*EE*Lambda/MW*Mq/SW*Zcc11*Zh33*Zss12/cb
     &      *sb+EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc12*Zh31
     &      *Zss12/cb-EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc12
     &      *Zh32*Zss12/sb
      AAABR(1360) = EE**2/SW**2*Sqrt2*Zdd11*Zh21*Zuu11/cb
     &      -EE**2/SW**2*Sqrt2*Zdd11*Zh22*Zuu11*sb-EE**2
     &      /SW**2*Sqrt2*Zdd11*Zh21*Zuu11/cb*sb**2-EE**2
     &      /MW**2*Mq**2/SW**2*Sqrt2*Zdd11*Zh21*Zuu11/cb
     &      +EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd11*Zh22*Zuu11
     &      /sb+EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd12*Zh21
     &      *Zuu12/cb-EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd12
     &      *Zh22*Zuu12/sb-2D0*EE*Lambda/MW*Mq/SW*Zdd11*Zh23
     &      *Zuu12/cb/sb+2D0*EE*Lambda/MW*Mq/SW*Zdd11*Zh23
     &      *Zuu12/cb*sb+2D0*EE*Lambda/MW*Mq/SW*Zdd12*Zh23
     &      *Zuu11/cb*sb
      AAABR(1361) = EE**2/SW**2*Sqrt2*Zb11*Zh21*Zt11/cb-EE
     &      **2/SW**2*Sqrt2*Zb11*Zh22*Zt11*sb-EE**2/SW**2
     &      *Sqrt2*Zb11*Zh21*Zt11/cb*sb**2-EE**2/MW**2*MbMM
     &      **2/SW**2*Sqrt2*Zb11*Zh21*Zt11/cb+EE**2/MW**2
     &      *MtMM**2/SW**2*Sqrt2*Zb11*Zh22*Zt11/sb+2D0*EE
     &      *Lambda/MW*MbMM/SW*Zb12*Zh23*Zt11/cb*sb-2D0*EE
     &      *Lambda/MW*MtMM/SW*Zb11*Zh23*Zt12/cb/sb+2D0*EE
     &      *Lambda/MW*MtMM/SW*Zb11*Zh23*Zt12/cb*sb+EE**2
     &      /MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb12*Zh21*Zt12/cb
     &      -EE**2/MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb12*Zh22
     &      *Zt12/sb
      AAABR(1362) = EE**2/SW**2*Sqrt2*Zcc11*Zh21*Zss11/cb
     &      -EE**2/SW**2*Sqrt2*Zcc11*Zh22*Zss11*sb-EE**2
     &      /SW**2*Sqrt2*Zcc11*Zh21*Zss11/cb*sb**2+EE**2
     &      /MW**2*Mc**2/SW**2*Sqrt2*Zcc11*Zh22*Zss11/sb
     &      -EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zcc11*Zh21*Zss11
     &      /cb-2D0*EE*Lambda/MW*Mc/SW*Zcc12*Zh23*Zss11/cb
     &      /sb+2D0*EE*Lambda/MW*Mc/SW*Zcc12*Zh23*Zss11/cb
     &      *sb+2D0*EE*Lambda/MW*Mq/SW*Zcc11*Zh23*Zss12/cb
     &      *sb+EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc12*Zh21
     &      *Zss12/cb-EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc12
     &      *Zh22*Zss12/sb
      AAABR(1363) = EE**2/SW**2*Sqrt2*Zdd11*Zh31*Zuu11/cb
     &      -EE**2/SW**2*Sqrt2*Zdd11*Zh32*Zuu11*sb-EE**2
     &      /SW**2*Sqrt2*Zdd11*Zh31*Zuu11/cb*sb**2-EE**2
     &      /MW**2*Mq**2/SW**2*Sqrt2*Zdd11*Zh31*Zuu11/cb
     &      +EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd11*Zh32*Zuu11
     &      /sb+EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd12*Zh31
     &      *Zuu12/cb-EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd12
     &      *Zh32*Zuu12/sb-2D0*EE*Lambda/MW*Mq/SW*Zdd11*Zh33
     &      *Zuu12/cb/sb+2D0*EE*Lambda/MW*Mq/SW*Zdd11*Zh33
     &      *Zuu12/cb*sb+2D0*EE*Lambda/MW*Mq/SW*Zdd12*Zh33
     &      *Zuu11/cb*sb
      AAABR(1364) = EE**2/SW**2*Sqrt2*Zcc11*Zh11*Zss11/cb
     &      -EE**2/SW**2*Sqrt2*Zcc11*Zh12*Zss11*sb-EE**2
     &      /SW**2*Sqrt2*Zcc11*Zh11*Zss11/cb*sb**2+EE**2
     &      /MW**2*Mc**2/SW**2*Sqrt2*Zcc11*Zh12*Zss11/sb
     &      -EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zcc11*Zh11*Zss11
     &      /cb-2D0*EE*Lambda/MW*Mc/SW*Zcc12*Zh13*Zss11/cb
     &      /sb+2D0*EE*Lambda/MW*Mc/SW*Zcc12*Zh13*Zss11/cb
     &      *sb+2D0*EE*Lambda/MW*Mq/SW*Zcc11*Zh13*Zss12/cb
     &      *sb+EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc12*Zh11
     &      *Zss12/cb-EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc12
     &      *Zh12*Zss12/sb
      AAABR(1365) = EE**2/SW**2*Sqrt2*Zb11*Zh11*Zt11/cb-EE
     &      **2/SW**2*Sqrt2*Zb11*Zh12*Zt11*sb-EE**2/SW**2
     &      *Sqrt2*Zb11*Zh11*Zt11/cb*sb**2-EE**2/MW**2*MbMM
     &      **2/SW**2*Sqrt2*Zb11*Zh11*Zt11/cb+EE**2/MW**2
     &      *MtMM**2/SW**2*Sqrt2*Zb11*Zh12*Zt11/sb+2D0*EE
     &      *Lambda/MW*MbMM/SW*Zb12*Zh13*Zt11/cb*sb-2D0*EE
     &      *Lambda/MW*MtMM/SW*Zb11*Zh13*Zt12/cb/sb+2D0*EE
     &      *Lambda/MW*MtMM/SW*Zb11*Zh13*Zt12/cb*sb+EE**2
     &      /MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb12*Zh11*Zt12/cb
     &      -EE**2/MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb12*Zh12
     &      *Zt12/sb
      AAABR(1366) = EE**2/SW**2*Sqrt2*Zb11*Zh31*Zt11/cb-EE
     &      **2/SW**2*Sqrt2*Zb11*Zh32*Zt11*sb-EE**2/SW**2
     &      *Sqrt2*Zb11*Zh31*Zt11/cb*sb**2-EE**2/MW**2*MbMM
     &      **2/SW**2*Sqrt2*Zb11*Zh31*Zt11/cb+EE**2/MW**2
     &      *MtMM**2/SW**2*Sqrt2*Zb11*Zh32*Zt11/sb+2D0*EE
     &      *Lambda/MW*MbMM/SW*Zb12*Zh33*Zt11/cb*sb-2D0*EE
     &      *Lambda/MW*MtMM/SW*Zb11*Zh33*Zt12/cb/sb+2D0*EE
     &      *Lambda/MW*MtMM/SW*Zb11*Zh33*Zt12/cb*sb+EE**2
     &      /MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb12*Zh31*Zt12/cb
     &      -EE**2/MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb12*Zh32
     &      *Zt12/sb
      AAABR(1367) = EE**2/SW**2*Sqrt2*Zcc21*Zh11*Zss11/cb
     &      -EE**2/SW**2*Sqrt2*Zcc21*Zh12*Zss11*sb-EE**2
     &      /SW**2*Sqrt2*Zcc21*Zh11*Zss11/cb*sb**2+EE**2
     &      /MW**2*Mc**2/SW**2*Sqrt2*Zcc21*Zh12*Zss11/sb
     &      -EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zcc21*Zh11*Zss11
     &      /cb-2D0*EE*Lambda/MW*Mc/SW*Zcc22*Zh13*Zss11/cb
     &      /sb+2D0*EE*Lambda/MW*Mc/SW*Zcc22*Zh13*Zss11/cb
     &      *sb+2D0*EE*Lambda/MW*Mq/SW*Zcc21*Zh13*Zss12/cb
     &      *sb+EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc22*Zh11
     &      *Zss12/cb-EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc22
     &      *Zh12*Zss12/sb
      AAABR(1368) = EE**2/SW**2*Sqrt2*Zcc21*Zh21*Zss11/cb
     &      -EE**2/SW**2*Sqrt2*Zcc21*Zh22*Zss11*sb-EE**2
     &      /SW**2*Sqrt2*Zcc21*Zh21*Zss11/cb*sb**2+EE**2
     &      /MW**2*Mc**2/SW**2*Sqrt2*Zcc21*Zh22*Zss11/sb
     &      -EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zcc21*Zh21*Zss11
     &      /cb-2D0*EE*Lambda/MW*Mc/SW*Zcc22*Zh23*Zss11/cb
     &      /sb+2D0*EE*Lambda/MW*Mc/SW*Zcc22*Zh23*Zss11/cb
     &      *sb+2D0*EE*Lambda/MW*Mq/SW*Zcc21*Zh23*Zss12/cb
     &      *sb+EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc22*Zh21
     &      *Zss12/cb-EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc22
     &      *Zh22*Zss12/sb
      AAABR(1369) = EE**2/SW**2*Sqrt2*Zdd11*Zh21*Zuu21/cb
     &      -EE**2/SW**2*Sqrt2*Zdd11*Zh22*Zuu21*sb-EE**2
     &      /SW**2*Sqrt2*Zdd11*Zh21*Zuu21/cb*sb**2-EE**2
     &      /MW**2*Mq**2/SW**2*Sqrt2*Zdd11*Zh21*Zuu21/cb
     &      +EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd11*Zh22*Zuu21
     &      /sb+EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd12*Zh21
     &      *Zuu22/cb-EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd12
     &      *Zh22*Zuu22/sb-2D0*EE*Lambda/MW*Mq/SW*Zdd11*Zh23
     &      *Zuu22/cb/sb+2D0*EE*Lambda/MW*Mq/SW*Zdd11*Zh23
     &      *Zuu22/cb*sb+2D0*EE*Lambda/MW*Mq/SW*Zdd12*Zh23
     &      *Zuu21/cb*sb
      AAABR(1370) = EE**2/SW**2*Sqrt2*Zcc21*Zh31*Zss11/cb
     &      -EE**2/SW**2*Sqrt2*Zcc21*Zh32*Zss11*sb-EE**2
     &      /SW**2*Sqrt2*Zcc21*Zh31*Zss11/cb*sb**2+EE**2
     &      /MW**2*Mc**2/SW**2*Sqrt2*Zcc21*Zh32*Zss11/sb
     &      -EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zcc21*Zh31*Zss11
     &      /cb-2D0*EE*Lambda/MW*Mc/SW*Zcc22*Zh33*Zss11/cb
     &      /sb+2D0*EE*Lambda/MW*Mc/SW*Zcc22*Zh33*Zss11/cb
     &      *sb+2D0*EE*Lambda/MW*Mq/SW*Zcc21*Zh33*Zss12/cb
     &      *sb+EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc22*Zh31
     &      *Zss12/cb-EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc22
     &      *Zh32*Zss12/sb
      AAABR(1371) = EE**2/SW**2*Sqrt2*Zb11*Zh11*Zt21/cb-EE
     &      **2/SW**2*Sqrt2*Zb11*Zh12*Zt21*sb-EE**2/SW**2
     &      *Sqrt2*Zb11*Zh11*Zt21/cb*sb**2-EE**2/MW**2*MbMM
     &      **2/SW**2*Sqrt2*Zb11*Zh11*Zt21/cb+EE**2/MW**2
     &      *MtMM**2/SW**2*Sqrt2*Zb11*Zh12*Zt21/sb+2D0*EE
     &      *Lambda/MW*MbMM/SW*Zb12*Zh13*Zt21/cb*sb-2D0*EE
     &      *Lambda/MW*MtMM/SW*Zb11*Zh13*Zt22/cb/sb+2D0*EE
     &      *Lambda/MW*MtMM/SW*Zb11*Zh13*Zt22/cb*sb+EE**2
     &      /MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb12*Zh11*Zt22/cb
     &      -EE**2/MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb12*Zh12
     &      *Zt22/sb
      AAABR(1372) = EE**2/SW**2*Sqrt2*Zdd11*Zh11*Zuu21/cb
     &      -EE**2/SW**2*Sqrt2*Zdd11*Zh12*Zuu21*sb-EE**2
     &      /SW**2*Sqrt2*Zdd11*Zh11*Zuu21/cb*sb**2-EE**2
     &      /MW**2*Mq**2/SW**2*Sqrt2*Zdd11*Zh11*Zuu21/cb
     &      +EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd11*Zh12*Zuu21
     &      /sb+EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd12*Zh11
     &      *Zuu22/cb-EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd12
     &      *Zh12*Zuu22/sb-2D0*EE*Lambda/MW*Mq/SW*Zdd11*Zh13
     &      *Zuu22/cb/sb+2D0*EE*Lambda/MW*Mq/SW*Zdd11*Zh13
     &      *Zuu22/cb*sb+2D0*EE*Lambda/MW*Mq/SW*Zdd12*Zh13
     &      *Zuu21/cb*sb
      AAABR(1373) = EE**2/SW**2*Sqrt2*Zb11*Zh21*Zt21/cb-EE
     &      **2/SW**2*Sqrt2*Zb11*Zh22*Zt21*sb-EE**2/SW**2
     &      *Sqrt2*Zb11*Zh21*Zt21/cb*sb**2-EE**2/MW**2*MbMM
     &      **2/SW**2*Sqrt2*Zb11*Zh21*Zt21/cb+EE**2/MW**2
     &      *MtMM**2/SW**2*Sqrt2*Zb11*Zh22*Zt21/sb+2D0*EE
     &      *Lambda/MW*MbMM/SW*Zb12*Zh23*Zt21/cb*sb-2D0*EE
     &      *Lambda/MW*MtMM/SW*Zb11*Zh23*Zt22/cb/sb+2D0*EE
     &      *Lambda/MW*MtMM/SW*Zb11*Zh23*Zt22/cb*sb+EE**2
     &      /MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb12*Zh21*Zt22/cb
     &      -EE**2/MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb12*Zh22
     &      *Zt22/sb
      AAABR(1374) = EE**2/SW**2*Sqrt2*Zdd11*Zh31*Zuu21/cb
     &      -EE**2/SW**2*Sqrt2*Zdd11*Zh32*Zuu21*sb-EE**2
     &      /SW**2*Sqrt2*Zdd11*Zh31*Zuu21/cb*sb**2-EE**2
     &      /MW**2*Mq**2/SW**2*Sqrt2*Zdd11*Zh31*Zuu21/cb
     &      +EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd11*Zh32*Zuu21
     &      /sb+EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd12*Zh31
     &      *Zuu22/cb-EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd12
     &      *Zh32*Zuu22/sb-2D0*EE*Lambda/MW*Mq/SW*Zdd11*Zh33
     &      *Zuu22/cb/sb+2D0*EE*Lambda/MW*Mq/SW*Zdd11*Zh33
     &      *Zuu22/cb*sb+2D0*EE*Lambda/MW*Mq/SW*Zdd12*Zh33
     &      *Zuu21/cb*sb
      AAABR(1375) = EE**2/SW**2*Sqrt2*Zb11*Zh31*Zt21/cb-EE
     &      **2/SW**2*Sqrt2*Zb11*Zh32*Zt21*sb-EE**2/SW**2
     &      *Sqrt2*Zb11*Zh31*Zt21/cb*sb**2-EE**2/MW**2*MbMM
     &      **2/SW**2*Sqrt2*Zb11*Zh31*Zt21/cb+EE**2/MW**2
     &      *MtMM**2/SW**2*Sqrt2*Zb11*Zh32*Zt21/sb+2D0*EE
     &      *Lambda/MW*MbMM/SW*Zb12*Zh33*Zt21/cb*sb-2D0*EE
     &      *Lambda/MW*MtMM/SW*Zb11*Zh33*Zt22/cb/sb+2D0*EE
     &      *Lambda/MW*MtMM/SW*Zb11*Zh33*Zt22/cb*sb+EE**2
     &      /MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb12*Zh31*Zt22/cb
     &      -EE**2/MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb12*Zh32
     &      *Zt22/sb
      AAABR(1376) = EE**2*(Zb11**2/CW**2+2D0*Zb12**2/CW**2
     &      +3D0*Zb11**2/SW**2-2D0*Zb11**2/CW**2*sb**2-4D0
     &      *Zb12**2/CW**2*sb**2-6D0*MbMM**2/MW**2/SW**2
     &      -6D0*Zb11**2/SW**2*sb**2)
      AAABR(1377) = EE**2*(Zdd11**2/CW**2+2D0*Zdd12**2/CW
     &      **2+3D0*Zdd11**2/SW**2-2D0*Zdd11**2/CW**2*sb
     &      **2-4D0*Zdd12**2/CW**2*sb**2-6D0*Mq**2/MW**2
     &      /SW**2-6D0*Zdd11**2/SW**2*sb**2)
      AAABR(1378) = EE**2*(Zss11**2/CW**2+2D0*Zss12**2/CW
     &      **2+3D0*Zss11**2/SW**2-2D0*Zss11**2/CW**2*sb
     &      **2-4D0*Zss12**2/CW**2*sb**2-6D0*Mq**2/MW**2
     &      /SW**2-6D0*Zss11**2/SW**2*sb**2)
      AAABR(1379) = EE**2*(Zdd11/CW**2*Zdd21+2D0*Zdd12/CW
     &      **2*Zdd22+3D0*Zdd11/SW**2*Zdd21-2D0*Zdd11/CW
     &      **2*Zdd21*sb**2-4D0*Zdd12/CW**2*Zdd22*sb**2-6D0
     &      *Zdd11/SW**2*Zdd21*sb**2)
      AAABR(1380) = EE**2*(Zss11/CW**2*Zss21+2D0*Zss12/CW
     &      **2*Zss22+3D0*Zss11/SW**2*Zss21-2D0*Zss11/CW
     &      **2*Zss21*sb**2-4D0*Zss12/CW**2*Zss22*sb**2-6D0
     &      *Zss11/SW**2*Zss21*sb**2)
      AAABR(1381) = EE**2*(Zb11/CW**2*Zb21+2D0*Zb12/CW**2
     &      *Zb22+3D0*Zb11/SW**2*Zb21-2D0*Zb11/CW**2*Zb21
     &      *sb**2-4D0*Zb12/CW**2*Zb22*sb**2-6D0*Zb11/SW
     &      **2*Zb21*sb**2)
      AAABR(1382) = tB*(EE**2/CW**2*Pa21*Zss11**2+2D0*EE
     &      **2/CW**2*Pa21*Zss12**2+3D0*EE**2*Pa21/SW**2
     &      *Zss11**2-EE**2/CW**2*Pa21*Zss11**2*sb**2-2D0
     &      *EE**2/CW**2*Pa21*Zss12**2*sb**2-3D0*EE**2/MW
     &      **2*Mq**2*Pa21/SW**2-3D0*EE**2*Pa21/SW**2*Zss11
     &      **2*sb**2+3D0*EE*Lambda/MW*Mq*Pa22/SW*Sqrt2*Zss11
     &      *Zss12)
      AAABR(1383) = tB*(EE**2/CW**2*Pa11*Zss11**2+2D0*EE
     &      **2/CW**2*Pa11*Zss12**2+3D0*EE**2*Pa11/SW**2
     &      *Zss11**2-EE**2/CW**2*Pa11*Zss11**2*sb**2-2D0
     &      *EE**2/CW**2*Pa11*Zss12**2*sb**2-3D0*EE**2/MW
     &      **2*Mq**2*Pa11/SW**2-3D0*EE**2*Pa11/SW**2*Zss11
     &      **2*sb**2+3D0*EE*Lambda/MW*Mq*Pa12/SW*Sqrt2*Zss11
     &      *Zss12)
      AAABR(1384) = tB*(EE**2/CW**2*Pa11*Zdd11**2+2D0*EE
     &      **2/CW**2*Pa11*Zdd12**2+3D0*EE**2*Pa11/SW**2
     &      *Zdd11**2-EE**2/CW**2*Pa11*Zdd11**2*sb**2-2D0
     &      *EE**2/CW**2*Pa11*Zdd12**2*sb**2-3D0*EE**2/MW
     &      **2*Mq**2*Pa11/SW**2-3D0*EE**2*Pa11/SW**2*Zdd11
     &      **2*sb**2+3D0*EE*Lambda/MW*Mq*Pa12/SW*Sqrt2*Zdd11
     &      *Zdd12)
      AAABR(1385) = tB*(EE**2/CW**2*Pa21*Zdd11**2+2D0*EE
     &      **2/CW**2*Pa21*Zdd12**2+3D0*EE**2*Pa21/SW**2
     &      *Zdd11**2-EE**2/CW**2*Pa21*Zdd11**2*sb**2-2D0
     &      *EE**2/CW**2*Pa21*Zdd12**2*sb**2-3D0*EE**2/MW
     &      **2*Mq**2*Pa21/SW**2-3D0*EE**2*Pa21/SW**2*Zdd11
     &      **2*sb**2+3D0*EE*Lambda/MW*Mq*Pa22/SW*Sqrt2*Zdd11
     &      *Zdd12)
      AAABR(1386) = tB*(EE**2/CW**2*Pa21*Zb11**2+2D0*EE**2
     &      /CW**2*Pa21*Zb12**2+3D0*EE**2*Pa21/SW**2*Zb11
     &      **2-EE**2/CW**2*Pa21*Zb11**2*sb**2-2D0*EE**2
     &      /CW**2*Pa21*Zb12**2*sb**2-3D0*EE**2/MW**2*MbMM
     &      **2*Pa21/SW**2-3D0*EE**2*Pa21/SW**2*Zb11**2*sb
     &      **2+3D0*EE*Lambda/MW*MbMM*Pa22/SW*Sqrt2*Zb11
     &      *Zb12)
      AAABR(1387) = tB*(EE**2/CW**2*Pa11*Zb11**2+2D0*EE**2
     &      /CW**2*Pa11*Zb12**2+3D0*EE**2*Pa11/SW**2*Zb11
     &      **2-EE**2/CW**2*Pa11*Zb11**2*sb**2-2D0*EE**2
     &      /CW**2*Pa11*Zb12**2*sb**2-3D0*EE**2/MW**2*MbMM
     &      **2*Pa11/SW**2-3D0*EE**2*Pa11/SW**2*Zb11**2*sb
     &      **2+3D0*EE*Lambda/MW*MbMM*Pa12/SW*Sqrt2*Zb11
     &      *Zb12)
      AAABR(1388) = tB*(2D0*EE**2/CW**2*Pa11*Zss11*Zss21
     &      +4D0*EE**2/CW**2*Pa11*Zss12*Zss22+6D0*EE**2*Pa11
     &      /SW**2*Zss11*Zss21-2D0*EE**2/CW**2*Pa11*Zss11
     &      *Zss21*sb**2-4D0*EE**2/CW**2*Pa11*Zss12*Zss22
     &      *sb**2-6D0*EE**2*Pa11/SW**2*Zss11*Zss21*sb**2
     &      +3D0*EE*Lambda/MW*Mq*Pa12/SW*Sqrt2*Zss11*Zss22
     &      +3D0*EE*Lambda/MW*Mq*Pa12/SW*Sqrt2*Zss12*Zss21)
      AAABR(1389) = tB*(2D0*EE**2/CW**2*Pa11*Zdd11*Zdd21
     &      +4D0*EE**2/CW**2*Pa11*Zdd12*Zdd22+6D0*EE**2*Pa11
     &      /SW**2*Zdd11*Zdd21-2D0*EE**2/CW**2*Pa11*Zdd11
     &      *Zdd21*sb**2-4D0*EE**2/CW**2*Pa11*Zdd12*Zdd22
     &      *sb**2-6D0*EE**2*Pa11/SW**2*Zdd11*Zdd21*sb**2
     &      +3D0*EE*Lambda/MW*Mq*Pa12/SW*Sqrt2*Zdd11*Zdd22
     &      +3D0*EE*Lambda/MW*Mq*Pa12/SW*Sqrt2*Zdd12*Zdd21)
      AAABR(1390) = tB*(2D0*EE**2/CW**2*Pa21*Zb11*Zb21+4D0
     &      *EE**2/CW**2*Pa21*Zb12*Zb22+6D0*EE**2*Pa21/SW
     &      **2*Zb11*Zb21-2D0*EE**2/CW**2*Pa21*Zb11*Zb21
     &      *sb**2-4D0*EE**2/CW**2*Pa21*Zb12*Zb22*sb**2-6D0
     &      *EE**2*Pa21/SW**2*Zb11*Zb21*sb**2+3D0*EE*Lambda
     &      /MW*MbMM*Pa22/SW*Sqrt2*Zb11*Zb22+3D0*EE*Lambda
     &      /MW*MbMM*Pa22/SW*Sqrt2*Zb12*Zb21)
      AAABR(1391) = tB*(2D0*EE**2/CW**2*Pa21*Zss11*Zss21
     &      +4D0*EE**2/CW**2*Pa21*Zss12*Zss22+6D0*EE**2*Pa21
     &      /SW**2*Zss11*Zss21-2D0*EE**2/CW**2*Pa21*Zss11
     &      *Zss21*sb**2-4D0*EE**2/CW**2*Pa21*Zss12*Zss22
     &      *sb**2-6D0*EE**2*Pa21/SW**2*Zss11*Zss21*sb**2
     &      +3D0*EE*Lambda/MW*Mq*Pa22/SW*Sqrt2*Zss11*Zss22
     &      +3D0*EE*Lambda/MW*Mq*Pa22/SW*Sqrt2*Zss12*Zss21)
      AAABR(1392) = tB*(2D0*EE**2/CW**2*Pa21*Zdd11*Zdd21
     &      +4D0*EE**2/CW**2*Pa21*Zdd12*Zdd22+6D0*EE**2*Pa21
     &      /SW**2*Zdd11*Zdd21-2D0*EE**2/CW**2*Pa21*Zdd11
     &      *Zdd21*sb**2-4D0*EE**2/CW**2*Pa21*Zdd12*Zdd22
     &      *sb**2-6D0*EE**2*Pa21/SW**2*Zdd11*Zdd21*sb**2
     &      +3D0*EE*Lambda/MW*Mq*Pa22/SW*Sqrt2*Zdd11*Zdd22
     &      +3D0*EE*Lambda/MW*Mq*Pa22/SW*Sqrt2*Zdd12*Zdd21)
      AAABR(1393) = tB*(2D0*EE**2/CW**2*Pa11*Zb11*Zb21+4D0
     &      *EE**2/CW**2*Pa11*Zb12*Zb22+6D0*EE**2*Pa11/SW
     &      **2*Zb11*Zb21-2D0*EE**2/CW**2*Pa11*Zb11*Zb21
     &      *sb**2-4D0*EE**2/CW**2*Pa11*Zb12*Zb22*sb**2-6D0
     &      *EE**2*Pa11/SW**2*Zb11*Zb21*sb**2+3D0*EE*Lambda
     &      /MW*MbMM*Pa12/SW*Sqrt2*Zb11*Zb22+3D0*EE*Lambda
     &      /MW*MbMM*Pa12/SW*Sqrt2*Zb12*Zb21)
      AAABR(1394) = EE*Lambda/MW*Mq/SW*Sqrt2*Zh33*tB*(Zss11
     &      *Zss22-Zss12*Zss21)
      AAABR(1395) = EE*Lambda/MW*Mq/SW*Sqrt2*Zh33*tB*(Zdd11
     &      *Zdd22-Zdd12*Zdd21)
      AAABR(1396) = EE*Lambda/MW*Mq/SW*Sqrt2*Zh13*tB*(Zss11
     &      *Zss22-Zss12*Zss21)
      AAABR(1397) = EE*Lambda/MW*Mq/SW*Sqrt2*Zh23*tB*(Zss11
     &      *Zss22-Zss12*Zss21)
      AAABR(1398) = EE*Lambda/MW*MbMM/SW*Sqrt2*Zh33*tB*(Zb11
     &      *Zb22-Zb12*Zb21)
      AAABR(1399) = EE*Lambda/MW*MbMM/SW*Sqrt2*Zh23*tB*(Zb11
     &      *Zb22-Zb12*Zb21)
      AAABR(1400) = EE*Lambda/MW*Mq/SW*Sqrt2*Zh13*tB*(Zdd11
     &      *Zdd22-Zdd12*Zdd21)
      AAABR(1401) = EE*Lambda/MW*MbMM/SW*Sqrt2*Zh13*tB*(Zb11
     &      *Zb22-Zb12*Zb21)
      AAABR(1402) = EE*Lambda/MW*Mq/SW*Sqrt2*Zh23*tB*(Zdd11
     &      *Zdd22-Zdd12*Zdd21)
      AAABR(1403) = 1D0/cb**2*(EE**2/CW**2*Pa21**2*Zb11**2
     &      +2D0*EE**2/CW**2*Pa21**2*Zb12**2+3D0*EE**2*Pa21
     &      **2/SW**2*Zb11**2-3D0*EE**2/CW**2*Pa21**2*Zb11
     &      **2*sb**2+2D0*EE**2/CW**2*Pa21**2*Zb11**2*sb
     &      **4-6D0*EE**2/CW**2*Pa21**2*Zb12**2*sb**2+4D0
     &      *EE**2/CW**2*Pa21**2*Zb12**2*sb**4-9D0*EE**2
     &      *Pa21**2/SW**2*Zb11**2*sb**2+6D0*EE**2*Pa21**2
     &      /SW**2*Zb11**2*sb**4+6D0*EE**2/MW**2*MbMM**2
     &      *Pa21**2/SW**2*sb**2+12D0*EE*Lambda/MW*MbMM*Pa21
     &      *Pa22/SW*Sqrt2*Zb11*Zb12-12D0*EE*Lambda/MW*MbMM
     &      *Pa21*Pa22/SW*Sqrt2*Zb11*Zb12*sb**2)
      AAABR(1404) = 1D0/cb**2*(EE**2/CW**2*Pa11**2*Zss11
     &      **2+2D0*EE**2/CW**2*Pa11**2*Zss12**2+3D0*EE**2
     &      *Pa11**2/SW**2*Zss11**2-3D0*EE**2/CW**2*Pa11
     &      **2*Zss11**2*sb**2+2D0*EE**2/CW**2*Pa11**2*Zss11
     &      **2*sb**4-6D0*EE**2/CW**2*Pa11**2*Zss12**2*sb
     &      **2+4D0*EE**2/CW**2*Pa11**2*Zss12**2*sb**4-9D0
     &      *EE**2*Pa11**2/SW**2*Zss11**2*sb**2+6D0*EE**2
     &      *Pa11**2/SW**2*Zss11**2*sb**4+6D0*EE**2/MW**2
     &      *Mq**2*Pa11**2/SW**2*sb**2+12D0*EE*Lambda/MW
     &      *Mq*Pa11*Pa12/SW*Sqrt2*Zss11*Zss12-12D0*EE*Lambda
     &      /MW*Mq*Pa11*Pa12/SW*Sqrt2*Zss11*Zss12*sb**2)
      AAABR(1405) = 1D0/cb**2*(EE**2/CW**2*Pa11*Pa21*Zdd11
     &      **2+2D0*EE**2/CW**2*Pa11*Pa21*Zdd12**2+3D0*EE
     &      **2*Pa11*Pa21/SW**2*Zdd11**2-3D0*EE**2/CW**2
     &      *Pa11*Pa21*Zdd11**2*sb**2+2D0*EE**2/CW**2*Pa11
     &      *Pa21*Zdd11**2*sb**4-6D0*EE**2/CW**2*Pa11*Pa21
     &      *Zdd12**2*sb**2+4D0*EE**2/CW**2*Pa11*Pa21*Zdd12
     &      **2*sb**4-9D0*EE**2*Pa11*Pa21/SW**2*Zdd11**2
     &      *sb**2+6D0*EE**2*Pa11*Pa21/SW**2*Zdd11**2*sb
     &      **4+6D0*EE**2/MW**2*Mq**2*Pa11*Pa21/SW**2*sb
     &      **2+6D0*EE*Lambda/MW*Mq*Pa11*Pa22/SW*Sqrt2*Zdd11
     &      *Zdd12+6D0*EE*Lambda/MW*Mq*Pa12*Pa21/SW*Sqrt2
     &      *Zdd11*Zdd12-6D0*EE*Lambda/MW*Mq*Pa11*Pa22/SW
     &      *Sqrt2*Zdd11*Zdd12*sb**2-6D0*EE*Lambda/MW*Mq
     &      *Pa12*Pa21/SW*Sqrt2*Zdd11*Zdd12*sb**2)
      AAABR(1406) = 1D0/cb**2*(EE**2/CW**2*Pa21**2*Zss11
     &      **2+2D0*EE**2/CW**2*Pa21**2*Zss12**2+3D0*EE**2
     &      *Pa21**2/SW**2*Zss11**2-3D0*EE**2/CW**2*Pa21
     &      **2*Zss11**2*sb**2+2D0*EE**2/CW**2*Pa21**2*Zss11
     &      **2*sb**4-6D0*EE**2/CW**2*Pa21**2*Zss12**2*sb
     &      **2+4D0*EE**2/CW**2*Pa21**2*Zss12**2*sb**4-9D0
     &      *EE**2*Pa21**2/SW**2*Zss11**2*sb**2+6D0*EE**2
     &      *Pa21**2/SW**2*Zss11**2*sb**4+6D0*EE**2/MW**2
     &      *Mq**2*Pa21**2/SW**2*sb**2+12D0*EE*Lambda/MW
     &      *Mq*Pa21*Pa22/SW*Sqrt2*Zss11*Zss12-12D0*EE*Lambda
     &      /MW*Mq*Pa21*Pa22/SW*Sqrt2*Zss11*Zss12*sb**2)
      AAABR(1407) = 1D0/cb**2*(EE**2/CW**2*Pa11*Pa21*Zss11
     &      **2+2D0*EE**2/CW**2*Pa11*Pa21*Zss12**2+3D0*EE
     &      **2*Pa11*Pa21/SW**2*Zss11**2-3D0*EE**2/CW**2
     &      *Pa11*Pa21*Zss11**2*sb**2+2D0*EE**2/CW**2*Pa11
     &      *Pa21*Zss11**2*sb**4-6D0*EE**2/CW**2*Pa11*Pa21
     &      *Zss12**2*sb**2+4D0*EE**2/CW**2*Pa11*Pa21*Zss12
     &      **2*sb**4-9D0*EE**2*Pa11*Pa21/SW**2*Zss11**2
     &      *sb**2+6D0*EE**2*Pa11*Pa21/SW**2*Zss11**2*sb
     &      **4+6D0*EE**2/MW**2*Mq**2*Pa11*Pa21/SW**2*sb
     &      **2+6D0*EE*Lambda/MW*Mq*Pa11*Pa22/SW*Sqrt2*Zss11
     &      *Zss12+6D0*EE*Lambda/MW*Mq*Pa12*Pa21/SW*Sqrt2
     &      *Zss11*Zss12-6D0*EE*Lambda/MW*Mq*Pa11*Pa22/SW
     &      *Sqrt2*Zss11*Zss12*sb**2-6D0*EE*Lambda/MW*Mq
     &      *Pa12*Pa21/SW*Sqrt2*Zss11*Zss12*sb**2)
      AAABR(1408) = 1D0/cb**2*(EE**2/CW**2*Pa11**2*Zb11**2
     &      +2D0*EE**2/CW**2*Pa11**2*Zb12**2+3D0*EE**2*Pa11
     &      **2/SW**2*Zb11**2-3D0*EE**2/CW**2*Pa11**2*Zb11
     &      **2*sb**2+2D0*EE**2/CW**2*Pa11**2*Zb11**2*sb
     &      **4-6D0*EE**2/CW**2*Pa11**2*Zb12**2*sb**2+4D0
     &      *EE**2/CW**2*Pa11**2*Zb12**2*sb**4-9D0*EE**2
     &      *Pa11**2/SW**2*Zb11**2*sb**2+6D0*EE**2*Pa11**2
     &      /SW**2*Zb11**2*sb**4+6D0*EE**2/MW**2*MbMM**2
     &      *Pa11**2/SW**2*sb**2+12D0*EE*Lambda/MW*MbMM*Pa11
     &      *Pa12/SW*Sqrt2*Zb11*Zb12-12D0*EE*Lambda/MW*MbMM
     &      *Pa11*Pa12/SW*Sqrt2*Zb11*Zb12*sb**2)
      AAABR(1409) = 1D0/cb**2*(EE**2/CW**2*Pa11*Pa21*Zb11
     &      **2+2D0*EE**2/CW**2*Pa11*Pa21*Zb12**2+3D0*EE
     &      **2*Pa11*Pa21/SW**2*Zb11**2-3D0*EE**2/CW**2*Pa11
     &      *Pa21*Zb11**2*sb**2+2D0*EE**2/CW**2*Pa11*Pa21
     &      *Zb11**2*sb**4-6D0*EE**2/CW**2*Pa11*Pa21*Zb12
     &      **2*sb**2+4D0*EE**2/CW**2*Pa11*Pa21*Zb12**2*sb
     &      **4-9D0*EE**2*Pa11*Pa21/SW**2*Zb11**2*sb**2+6D0
     &      *EE**2*Pa11*Pa21/SW**2*Zb11**2*sb**4+6D0*EE**2
     &      /MW**2*MbMM**2*Pa11*Pa21/SW**2*sb**2+6D0*EE*Lambda
     &      /MW*MbMM*Pa11*Pa22/SW*Sqrt2*Zb11*Zb12+6D0*EE
     &      *Lambda/MW*MbMM*Pa12*Pa21/SW*Sqrt2*Zb11*Zb12
     &      -6D0*EE*Lambda/MW*MbMM*Pa11*Pa22/SW*Sqrt2*Zb11
     &      *Zb12*sb**2-6D0*EE*Lambda/MW*MbMM*Pa12*Pa21/SW
     &      *Sqrt2*Zb11*Zb12*sb**2)
      AAABR(1410) = 1D0/cb**2*(EE**2/CW**2*Pa11**2*Zdd11
     &      **2+2D0*EE**2/CW**2*Pa11**2*Zdd12**2+3D0*EE**2
     &      *Pa11**2/SW**2*Zdd11**2-3D0*EE**2/CW**2*Pa11
     &      **2*Zdd11**2*sb**2+2D0*EE**2/CW**2*Pa11**2*Zdd11
     &      **2*sb**4-6D0*EE**2/CW**2*Pa11**2*Zdd12**2*sb
     &      **2+4D0*EE**2/CW**2*Pa11**2*Zdd12**2*sb**4-9D0
     &      *EE**2*Pa11**2/SW**2*Zdd11**2*sb**2+6D0*EE**2
     &      *Pa11**2/SW**2*Zdd11**2*sb**4+6D0*EE**2/MW**2
     &      *Mq**2*Pa11**2/SW**2*sb**2+12D0*EE*Lambda/MW
     &      *Mq*Pa11*Pa12/SW*Sqrt2*Zdd11*Zdd12-12D0*EE*Lambda
     &      /MW*Mq*Pa11*Pa12/SW*Sqrt2*Zdd11*Zdd12*sb**2)
      AAABR(1411) = 1D0/cb**2*(EE**2/CW**2*Pa21**2*Zdd11
     &      **2+2D0*EE**2/CW**2*Pa21**2*Zdd12**2+3D0*EE**2
     &      *Pa21**2/SW**2*Zdd11**2-3D0*EE**2/CW**2*Pa21
     &      **2*Zdd11**2*sb**2+2D0*EE**2/CW**2*Pa21**2*Zdd11
     &      **2*sb**4-6D0*EE**2/CW**2*Pa21**2*Zdd12**2*sb
     &      **2+4D0*EE**2/CW**2*Pa21**2*Zdd12**2*sb**4-9D0
     &      *EE**2*Pa21**2/SW**2*Zdd11**2*sb**2+6D0*EE**2
     &      *Pa21**2/SW**2*Zdd11**2*sb**4+6D0*EE**2/MW**2
     &      *Mq**2*Pa21**2/SW**2*sb**2+12D0*EE*Lambda/MW
     &      *Mq*Pa21*Pa22/SW*Sqrt2*Zdd11*Zdd12-12D0*EE*Lambda
     &      /MW*Mq*Pa21*Pa22/SW*Sqrt2*Zdd11*Zdd12*sb**2)
      AAABR(1412) = 1D0/cb**2*(EE**2/CW**2*Pa11**2*Zss11
     &      *Zss21+2D0*EE**2/CW**2*Pa11**2*Zss12*Zss22+3D0
     &      *EE**2*Pa11**2/SW**2*Zss11*Zss21-3D0*EE**2/CW
     &      **2*Pa11**2*Zss11*Zss21*sb**2+2D0*EE**2/CW**2
     &      *Pa11**2*Zss11*Zss21*sb**4-6D0*EE**2/CW**2*Pa11
     &      **2*Zss12*Zss22*sb**2+4D0*EE**2/CW**2*Pa11**2
     &      *Zss12*Zss22*sb**4-9D0*EE**2*Pa11**2/SW**2*Zss11
     &      *Zss21*sb**2+6D0*EE**2*Pa11**2/SW**2*Zss11*Zss21
     &      *sb**4+6D0*EE*Lambda/MW*Mq*Pa11*Pa12/SW*Sqrt2
     &      *Zss11*Zss22+6D0*EE*Lambda/MW*Mq*Pa11*Pa12/SW
     &      *Sqrt2*Zss12*Zss21-6D0*EE*Lambda/MW*Mq*Pa11*Pa12
     &      /SW*Sqrt2*Zss11*Zss22*sb**2-6D0*EE*Lambda/MW
     &      *Mq*Pa11*Pa12/SW*Sqrt2*Zss12*Zss21*sb**2)
      AAABR(1413) = 1D0/cb**2*(EE**2/CW**2*Pa11**2*Zb11*Zb21
     &      +2D0*EE**2/CW**2*Pa11**2*Zb12*Zb22+3D0*EE**2
     &      *Pa11**2/SW**2*Zb11*Zb21-3D0*EE**2/CW**2*Pa11
     &      **2*Zb11*Zb21*sb**2+2D0*EE**2/CW**2*Pa11**2*Zb11
     &      *Zb21*sb**4-6D0*EE**2/CW**2*Pa11**2*Zb12*Zb22
     &      *sb**2+4D0*EE**2/CW**2*Pa11**2*Zb12*Zb22*sb**4
     &      -9D0*EE**2*Pa11**2/SW**2*Zb11*Zb21*sb**2+6D0
     &      *EE**2*Pa11**2/SW**2*Zb11*Zb21*sb**4+6D0*EE*Lambda
     &      /MW*MbMM*Pa11*Pa12/SW*Sqrt2*Zb11*Zb22+6D0*EE
     &      *Lambda/MW*MbMM*Pa11*Pa12/SW*Sqrt2*Zb12*Zb21
     &      -6D0*EE*Lambda/MW*MbMM*Pa11*Pa12/SW*Sqrt2*Zb11
     &      *Zb22*sb**2-6D0*EE*Lambda/MW*MbMM*Pa11*Pa12/SW
     &      *Sqrt2*Zb12*Zb21*sb**2)
      AAABR(1414) = 1D0/cb**2*(EE**2/CW**2*Pa21**2*Zdd11
     &      *Zdd21+2D0*EE**2/CW**2*Pa21**2*Zdd12*Zdd22+3D0
     &      *EE**2*Pa21**2/SW**2*Zdd11*Zdd21-3D0*EE**2/CW
     &      **2*Pa21**2*Zdd11*Zdd21*sb**2+2D0*EE**2/CW**2
     &      *Pa21**2*Zdd11*Zdd21*sb**4-6D0*EE**2/CW**2*Pa21
     &      **2*Zdd12*Zdd22*sb**2+4D0*EE**2/CW**2*Pa21**2
     &      *Zdd12*Zdd22*sb**4-9D0*EE**2*Pa21**2/SW**2*Zdd11
     &      *Zdd21*sb**2+6D0*EE**2*Pa21**2/SW**2*Zdd11*Zdd21
     &      *sb**4+6D0*EE*Lambda/MW*Mq*Pa21*Pa22/SW*Sqrt2
     &      *Zdd11*Zdd22+6D0*EE*Lambda/MW*Mq*Pa21*Pa22/SW
     &      *Sqrt2*Zdd12*Zdd21-6D0*EE*Lambda/MW*Mq*Pa21*Pa22
     &      /SW*Sqrt2*Zdd11*Zdd22*sb**2-6D0*EE*Lambda/MW
     &      *Mq*Pa21*Pa22/SW*Sqrt2*Zdd12*Zdd21*sb**2)
      AAABR(1415) = 1D0/cb**2*(EE**2/CW**2*Pa11**2*Zdd11
     &      *Zdd21+2D0*EE**2/CW**2*Pa11**2*Zdd12*Zdd22+3D0
     &      *EE**2*Pa11**2/SW**2*Zdd11*Zdd21-3D0*EE**2/CW
     &      **2*Pa11**2*Zdd11*Zdd21*sb**2+2D0*EE**2/CW**2
     &      *Pa11**2*Zdd11*Zdd21*sb**4-6D0*EE**2/CW**2*Pa11
     &      **2*Zdd12*Zdd22*sb**2+4D0*EE**2/CW**2*Pa11**2
     &      *Zdd12*Zdd22*sb**4-9D0*EE**2*Pa11**2/SW**2*Zdd11
     &      *Zdd21*sb**2+6D0*EE**2*Pa11**2/SW**2*Zdd11*Zdd21
     &      *sb**4+6D0*EE*Lambda/MW*Mq*Pa11*Pa12/SW*Sqrt2
     &      *Zdd11*Zdd22+6D0*EE*Lambda/MW*Mq*Pa11*Pa12/SW
     &      *Sqrt2*Zdd12*Zdd21-6D0*EE*Lambda/MW*Mq*Pa11*Pa12
     &      /SW*Sqrt2*Zdd11*Zdd22*sb**2-6D0*EE*Lambda/MW
     &      *Mq*Pa11*Pa12/SW*Sqrt2*Zdd12*Zdd21*sb**2)
      AAABR(1416) = 1D0/cb**2*(EE**2/CW**2*Pa11*Pa21*Zss11
     &      *Zss21+2D0*EE**2/CW**2*Pa11*Pa21*Zss12*Zss22
     &      +3D0*EE**2*Pa11*Pa21/SW**2*Zss11*Zss21-3D0*EE
     &      **2/CW**2*Pa11*Pa21*Zss11*Zss21*sb**2+2D0*EE
     &      **2/CW**2*Pa11*Pa21*Zss11*Zss21*sb**4-6D0*EE
     &      **2/CW**2*Pa11*Pa21*Zss12*Zss22*sb**2+4D0*EE
     &      **2/CW**2*Pa11*Pa21*Zss12*Zss22*sb**4-9D0*EE
     &      **2*Pa11*Pa21/SW**2*Zss11*Zss21*sb**2+6D0*EE
     &      **2*Pa11*Pa21/SW**2*Zss11*Zss21*sb**4+3D0*EE
     &      *Lambda/MW*Mq*Pa11*Pa22/SW*Sqrt2*Zss11*Zss22
     &      +3D0*EE*Lambda/MW*Mq*Pa11*Pa22/SW*Sqrt2*Zss12
     &      *Zss21+3D0*EE*Lambda/MW*Mq*Pa12*Pa21/SW*Sqrt2
     &      *Zss11*Zss22+3D0*EE*Lambda/MW*Mq*Pa12*Pa21/SW
     &      *Sqrt2*Zss12*Zss21-3D0*EE*Lambda/MW*Mq*Pa11*Pa22
     &      /SW*Sqrt2*Zss11*Zss22*sb**2-3D0*EE*Lambda/MW
     &      *Mq*Pa11*Pa22/SW*Sqrt2*Zss12*Zss21*sb**2-3D0
     &      *EE*Lambda/MW*Mq*Pa12*Pa21/SW*Sqrt2*Zss11*Zss22
     &      *sb**2-3D0*EE*Lambda/MW*Mq*Pa12*Pa21/SW*Sqrt2
     &      *Zss12*Zss21*sb**2)
      AAABR(1417) = 1D0/cb**2*(EE**2/CW**2*Pa11*Pa21*Zdd11
     &      *Zdd21+2D0*EE**2/CW**2*Pa11*Pa21*Zdd12*Zdd22
     &      +3D0*EE**2*Pa11*Pa21/SW**2*Zdd11*Zdd21-3D0*EE
     &      **2/CW**2*Pa11*Pa21*Zdd11*Zdd21*sb**2+2D0*EE
     &      **2/CW**2*Pa11*Pa21*Zdd11*Zdd21*sb**4-6D0*EE
     &      **2/CW**2*Pa11*Pa21*Zdd12*Zdd22*sb**2+4D0*EE
     &      **2/CW**2*Pa11*Pa21*Zdd12*Zdd22*sb**4-9D0*EE
     &      **2*Pa11*Pa21/SW**2*Zdd11*Zdd21*sb**2+6D0*EE
     &      **2*Pa11*Pa21/SW**2*Zdd11*Zdd21*sb**4+3D0*EE
     &      *Lambda/MW*Mq*Pa11*Pa22/SW*Sqrt2*Zdd11*Zdd22
     &      +3D0*EE*Lambda/MW*Mq*Pa11*Pa22/SW*Sqrt2*Zdd12
     &      *Zdd21+3D0*EE*Lambda/MW*Mq*Pa12*Pa21/SW*Sqrt2
     &      *Zdd11*Zdd22+3D0*EE*Lambda/MW*Mq*Pa12*Pa21/SW
     &      *Sqrt2*Zdd12*Zdd21-3D0*EE*Lambda/MW*Mq*Pa11*Pa22
     &      /SW*Sqrt2*Zdd11*Zdd22*sb**2-3D0*EE*Lambda/MW
     &      *Mq*Pa11*Pa22/SW*Sqrt2*Zdd12*Zdd21*sb**2-3D0
     &      *EE*Lambda/MW*Mq*Pa12*Pa21/SW*Sqrt2*Zdd11*Zdd22
     &      *sb**2-3D0*EE*Lambda/MW*Mq*Pa12*Pa21/SW*Sqrt2
     &      *Zdd12*Zdd21*sb**2)
      AAABR(1418) = 1D0/cb**2*(EE**2/CW**2*Pa11*Pa21*Zb11
     &      *Zb21+2D0*EE**2/CW**2*Pa11*Pa21*Zb12*Zb22+3D0
     &      *EE**2*Pa11*Pa21/SW**2*Zb11*Zb21-3D0*EE**2/CW
     &      **2*Pa11*Pa21*Zb11*Zb21*sb**2+2D0*EE**2/CW**2
     &      *Pa11*Pa21*Zb11*Zb21*sb**4-6D0*EE**2/CW**2*Pa11
     &      *Pa21*Zb12*Zb22*sb**2+4D0*EE**2/CW**2*Pa11*Pa21
     &      *Zb12*Zb22*sb**4-9D0*EE**2*Pa11*Pa21/SW**2*Zb11
     &      *Zb21*sb**2+6D0*EE**2*Pa11*Pa21/SW**2*Zb11*Zb21
     &      *sb**4+3D0*EE*Lambda/MW*MbMM*Pa11*Pa22/SW*Sqrt2
     &      *Zb11*Zb22+3D0*EE*Lambda/MW*MbMM*Pa11*Pa22/SW
     &      *Sqrt2*Zb12*Zb21+3D0*EE*Lambda/MW*MbMM*Pa12*Pa21
     &      /SW*Sqrt2*Zb11*Zb22+3D0*EE*Lambda/MW*MbMM*Pa12
     &      *Pa21/SW*Sqrt2*Zb12*Zb21-3D0*EE*Lambda/MW*MbMM
     &      *Pa11*Pa22/SW*Sqrt2*Zb11*Zb22*sb**2-3D0*EE*Lambda
     &      /MW*MbMM*Pa11*Pa22/SW*Sqrt2*Zb12*Zb21*sb**2-3D0
     &      *EE*Lambda/MW*MbMM*Pa12*Pa21/SW*Sqrt2*Zb11*Zb22
     &      *sb**2-3D0*EE*Lambda/MW*MbMM*Pa12*Pa21/SW*Sqrt2
     &      *Zb12*Zb21*sb**2)
      AAABR(1419) = 1D0/cb**2*(EE**2/CW**2*Pa21**2*Zb11*Zb21
     &      +2D0*EE**2/CW**2*Pa21**2*Zb12*Zb22+3D0*EE**2
     &      *Pa21**2/SW**2*Zb11*Zb21-3D0*EE**2/CW**2*Pa21
     &      **2*Zb11*Zb21*sb**2+2D0*EE**2/CW**2*Pa21**2*Zb11
     &      *Zb21*sb**4-6D0*EE**2/CW**2*Pa21**2*Zb12*Zb22
     &      *sb**2+4D0*EE**2/CW**2*Pa21**2*Zb12*Zb22*sb**4
     &      -9D0*EE**2*Pa21**2/SW**2*Zb11*Zb21*sb**2+6D0
     &      *EE**2*Pa21**2/SW**2*Zb11*Zb21*sb**4+6D0*EE*Lambda
     &      /MW*MbMM*Pa21*Pa22/SW*Sqrt2*Zb11*Zb22+6D0*EE
     &      *Lambda/MW*MbMM*Pa21*Pa22/SW*Sqrt2*Zb12*Zb21
     &      -6D0*EE*Lambda/MW*MbMM*Pa21*Pa22/SW*Sqrt2*Zb11
     &      *Zb22*sb**2-6D0*EE*Lambda/MW*MbMM*Pa21*Pa22/SW
     &      *Sqrt2*Zb12*Zb21*sb**2)
      AAABR(1420) = 1D0/cb**2*(EE**2/CW**2*Pa21**2*Zss11
     &      *Zss21+2D0*EE**2/CW**2*Pa21**2*Zss12*Zss22+3D0
     &      *EE**2*Pa21**2/SW**2*Zss11*Zss21-3D0*EE**2/CW
     &      **2*Pa21**2*Zss11*Zss21*sb**2+2D0*EE**2/CW**2
     &      *Pa21**2*Zss11*Zss21*sb**4-6D0*EE**2/CW**2*Pa21
     &      **2*Zss12*Zss22*sb**2+4D0*EE**2/CW**2*Pa21**2
     &      *Zss12*Zss22*sb**4-9D0*EE**2*Pa21**2/SW**2*Zss11
     &      *Zss21*sb**2+6D0*EE**2*Pa21**2/SW**2*Zss11*Zss21
     &      *sb**4+6D0*EE*Lambda/MW*Mq*Pa21*Pa22/SW*Sqrt2
     &      *Zss11*Zss22+6D0*EE*Lambda/MW*Mq*Pa21*Pa22/SW
     &      *Sqrt2*Zss12*Zss21-6D0*EE*Lambda/MW*Mq*Pa21*Pa22
     &      /SW*Sqrt2*Zss11*Zss22*sb**2-6D0*EE*Lambda/MW
     &      *Mq*Pa21*Pa22/SW*Sqrt2*Zss12*Zss21*sb**2)
      AAABR(1421) = EE*Lambda/MW*Mq/SW*Sqrt2*(Pa11*Zh13*Zss11
     &      *Zss22-Pa11*Zh13*Zss12*Zss21+Pa12*Zh12*Zss11
     &      *Zss22/cb-Pa12*Zh12*Zss12*Zss21/cb)
      AAABR(1422) = EE*Lambda/MW*Mq/SW*Sqrt2*(Pa11*Zh33*Zss11
     &      *Zss22-Pa11*Zh33*Zss12*Zss21+Pa12*Zh32*Zss11
     &      *Zss22/cb-Pa12*Zh32*Zss12*Zss21/cb)
      AAABR(1423) = EE*Lambda/MW*Mq/SW*Sqrt2*(Pa11*Zdd11
     &      *Zdd22*Zh13-Pa11*Zdd12*Zdd21*Zh13+Pa12*Zdd11
     &      *Zdd22*Zh12/cb-Pa12*Zdd12*Zdd21*Zh12/cb)
      AAABR(1424) = EE*Lambda/MW*MbMM/SW*Sqrt2*(Pa21*Zb11
     &      *Zb22*Zh13-Pa21*Zb12*Zb21*Zh13+Pa22*Zb11*Zb22
     &      *Zh12/cb-Pa22*Zb12*Zb21*Zh12/cb)
      AAABR(1425) = EE*Lambda/MW*Mq/SW*Sqrt2*(Pa11*Zdd11
     &      *Zdd22*Zh23-Pa11*Zdd12*Zdd21*Zh23+Pa12*Zdd11
     &      *Zdd22*Zh22/cb-Pa12*Zdd12*Zdd21*Zh22/cb)
      AAABR(1426) = EE*Lambda/MW*MbMM/SW*Sqrt2*(Pa21*Zb11
     &      *Zb22*Zh23-Pa21*Zb12*Zb21*Zh23+Pa22*Zb11*Zb22
     &      *Zh22/cb-Pa22*Zb12*Zb21*Zh22/cb)
      AAABR(1427) = EE*Lambda/MW*Mq/SW*Sqrt2*(Pa21*Zdd11
     &      *Zdd22*Zh33-Pa21*Zdd12*Zdd21*Zh33+Pa22*Zdd11
     &      *Zdd22*Zh32/cb-Pa22*Zdd12*Zdd21*Zh32/cb)
      AAABR(1428) = EE*Lambda/MW*MbMM/SW*Sqrt2*(Pa11*Zb11
     &      *Zb22*Zh23-Pa11*Zb12*Zb21*Zh23+Pa12*Zb11*Zb22
     &      *Zh22/cb-Pa12*Zb12*Zb21*Zh22/cb)
      AAABR(1429) = EE*Lambda/MW*Mq/SW*Sqrt2*(Pa21*Zh23*Zss11
     &      *Zss22-Pa21*Zh23*Zss12*Zss21+Pa22*Zh22*Zss11
     &      *Zss22/cb-Pa22*Zh22*Zss12*Zss21/cb)
      AAABR(1430) = EE*Lambda/MW*Mq/SW*Sqrt2*(Pa11*Zh23*Zss11
     &      *Zss22-Pa11*Zh23*Zss12*Zss21+Pa12*Zh22*Zss11
     &      *Zss22/cb-Pa12*Zh22*Zss12*Zss21/cb)
      AAABR(1431) = EE*Lambda/MW*MbMM/SW*Sqrt2*(Pa21*Zb11
     &      *Zb22*Zh33-Pa21*Zb12*Zb21*Zh33+Pa22*Zb11*Zb22
     &      *Zh32/cb-Pa22*Zb12*Zb21*Zh32/cb)
      AAABR(1432) = EE*Lambda/MW*Mq/SW*Sqrt2*(Pa11*Zdd11
     &      *Zdd22*Zh33-Pa11*Zdd12*Zdd21*Zh33+Pa12*Zdd11
     &      *Zdd22*Zh32/cb-Pa12*Zdd12*Zdd21*Zh32/cb)
      AAABR(1433) = EE*Lambda/MW*MbMM/SW*Sqrt2*(Pa11*Zb11
     &      *Zb22*Zh13-Pa11*Zb12*Zb21*Zh13+Pa12*Zb11*Zb22
     &      *Zh12/cb-Pa12*Zb12*Zb21*Zh12/cb)
      AAABR(1434) = EE*Lambda/MW*Mq/SW*Sqrt2*(Pa21*Zdd11
     &      *Zdd22*Zh13-Pa21*Zdd12*Zdd21*Zh13+Pa22*Zdd11
     &      *Zdd22*Zh12/cb-Pa22*Zdd12*Zdd21*Zh12/cb)
      AAABR(1435) = EE*Lambda/MW*Mq/SW*Sqrt2*(Pa21*Zh33*Zss11
     &      *Zss22-Pa21*Zh33*Zss12*Zss21+Pa22*Zh32*Zss11
     &      *Zss22/cb-Pa22*Zh32*Zss12*Zss21/cb)
      AAABR(1436) = EE*Lambda/MW*Mq/SW*Sqrt2*(Pa21*Zdd11
     &      *Zdd22*Zh23-Pa21*Zdd12*Zdd21*Zh23+Pa22*Zdd11
     &      *Zdd22*Zh22/cb-Pa22*Zdd12*Zdd21*Zh22/cb)
      AAABR(1437) = EE*Lambda/MW*MbMM/SW*Sqrt2*(Pa11*Zb11
     &      *Zb22*Zh33-Pa11*Zb12*Zb21*Zh33+Pa12*Zb11*Zb22
     &      *Zh32/cb-Pa12*Zb12*Zb21*Zh32/cb)
      AAABR(1438) = EE*Lambda/MW*Mq/SW*Sqrt2*(Pa21*Zh13*Zss11
     &      *Zss22-Pa21*Zh13*Zss12*Zss21+Pa22*Zh12*Zss11
     &      *Zss22/cb-Pa22*Zh12*Zss12*Zss21/cb)
      AAABR(1439) = EE**2/CW**2*Zh11**2*Zss11**2/cb**2+2D0
     &      *EE**2/CW**2*Zh11**2*Zss12**2/cb**2-EE**2/CW
     &      **2*Zh12**2*Zss11**2/cb**2-2D0*EE**2/CW**2*Zh12
     &      **2*Zss12**2/cb**2+3D0*EE**2/SW**2*Zh11**2*Zss11
     &      **2/cb**2-3D0*EE**2/SW**2*Zh12**2*Zss11**2/cb
     &      **2-EE**2/CW**2*Zh11**2*Zss11**2/cb**2*sb**2
     &      -2D0*EE**2/CW**2*Zh11**2*Zss12**2/cb**2*sb**2
     &      +EE**2/CW**2*Zh12**2*Zss11**2/cb**2*sb**2+2D0
     &      *EE**2/CW**2*Zh12**2*Zss12**2/cb**2*sb**2-6D0
     &      *EE**2/MW**2*Mq**2/SW**2*Zh11**2/cb**2-3D0*EE
     &      **2/SW**2*Zh11**2*Zss11**2/cb**2*sb**2+3D0*EE
     &      **2/SW**2*Zh12**2*Zss11**2/cb**2*sb**2+12D0*EE
     &      *Lambda/MW*Mq/SW*Sqrt2*Zh12*Zh13*Zss11*Zss12/cb
      AAABR(1440) = EE**2/CW**2*Zh11*Zh31*Zss11**2/cb**2
     &      +2D0*EE**2/CW**2*Zh11*Zh31*Zss12**2/cb**2-EE
     &      **2/CW**2*Zh12*Zh32*Zss11**2/cb**2-2D0*EE**2
     &      /CW**2*Zh12*Zh32*Zss12**2/cb**2+3D0*EE**2/SW
     &      **2*Zh11*Zh31*Zss11**2/cb**2-3D0*EE**2/SW**2
     &      *Zh12*Zh32*Zss11**2/cb**2-EE**2/CW**2*Zh11*Zh31
     &      *Zss11**2/cb**2*sb**2-2D0*EE**2/CW**2*Zh11*Zh31
     &      *Zss12**2/cb**2*sb**2+EE**2/CW**2*Zh12*Zh32*Zss11
     &      **2/cb**2*sb**2+2D0*EE**2/CW**2*Zh12*Zh32*Zss12
     &      **2/cb**2*sb**2-6D0*EE**2/MW**2*Mq**2/SW**2*Zh11
     &      *Zh31/cb**2-3D0*EE**2/SW**2*Zh11*Zh31*Zss11**2
     &      /cb**2*sb**2+3D0*EE**2/SW**2*Zh12*Zh32*Zss11
     &      **2/cb**2*sb**2+6D0*EE*Lambda/MW*Mq/SW*Sqrt2
     &      *Zh12*Zh33*Zss11*Zss12/cb+6D0*EE*Lambda/MW*Mq
     &      /SW*Sqrt2*Zh13*Zh32*Zss11*Zss12/cb
      AAABR(1441) = EE**2/CW**2*Zb11**2*Zh21**2/cb**2-EE
     &      **2/CW**2*Zb11**2*Zh22**2/cb**2+2D0*EE**2/CW
     &      **2*Zb12**2*Zh21**2/cb**2-2D0*EE**2/CW**2*Zb12
     &      **2*Zh22**2/cb**2+3D0*EE**2/SW**2*Zb11**2*Zh21
     &      **2/cb**2-3D0*EE**2/SW**2*Zb11**2*Zh22**2/cb
     &      **2-EE**2/CW**2*Zb11**2*Zh21**2/cb**2*sb**2+EE
     &      **2/CW**2*Zb11**2*Zh22**2/cb**2*sb**2-2D0*EE
     &      **2/CW**2*Zb12**2*Zh21**2/cb**2*sb**2+2D0*EE
     &      **2/CW**2*Zb12**2*Zh22**2/cb**2*sb**2-6D0*EE
     &      **2/MW**2*MbMM**2/SW**2*Zh21**2/cb**2-3D0*EE
     &      **2/SW**2*Zb11**2*Zh21**2/cb**2*sb**2+3D0*EE
     &      **2/SW**2*Zb11**2*Zh22**2/cb**2*sb**2+12D0*EE
     &      *Lambda/MW*MbMM/SW*Sqrt2*Zb11*Zb12*Zh22*Zh23/cb
      AAABR(1442) = EE**2/CW**2*Zb11**2*Zh11*Zh31/cb**2-EE
     &      **2/CW**2*Zb11**2*Zh12*Zh32/cb**2+2D0*EE**2/CW
     &      **2*Zb12**2*Zh11*Zh31/cb**2-2D0*EE**2/CW**2*Zb12
     &      **2*Zh12*Zh32/cb**2+3D0*EE**2/SW**2*Zb11**2*Zh11
     &      *Zh31/cb**2-3D0*EE**2/SW**2*Zb11**2*Zh12*Zh32
     &      /cb**2-EE**2/CW**2*Zb11**2*Zh11*Zh31/cb**2*sb
     &      **2+EE**2/CW**2*Zb11**2*Zh12*Zh32/cb**2*sb**2
     &      -2D0*EE**2/CW**2*Zb12**2*Zh11*Zh31/cb**2*sb**2
     &      +2D0*EE**2/CW**2*Zb12**2*Zh12*Zh32/cb**2*sb**2
     &      -6D0*EE**2/MW**2*MbMM**2/SW**2*Zh11*Zh31/cb**2
     &      -3D0*EE**2/SW**2*Zb11**2*Zh11*Zh31/cb**2*sb**2
     &      +3D0*EE**2/SW**2*Zb11**2*Zh12*Zh32/cb**2*sb**2
     &      +6D0*EE*Lambda/MW*MbMM/SW*Sqrt2*Zb11*Zb12*Zh12
     &      *Zh33/cb+6D0*EE*Lambda/MW*MbMM/SW*Sqrt2*Zb11
     &      *Zb12*Zh13*Zh32/cb
      AAABR(1443) = EE**2/CW**2*Zh21**2*Zss11**2/cb**2+2D0
     &      *EE**2/CW**2*Zh21**2*Zss12**2/cb**2-EE**2/CW
     &      **2*Zh22**2*Zss11**2/cb**2-2D0*EE**2/CW**2*Zh22
     &      **2*Zss12**2/cb**2+3D0*EE**2/SW**2*Zh21**2*Zss11
     &      **2/cb**2-3D0*EE**2/SW**2*Zh22**2*Zss11**2/cb
     &      **2-EE**2/CW**2*Zh21**2*Zss11**2/cb**2*sb**2
     &      -2D0*EE**2/CW**2*Zh21**2*Zss12**2/cb**2*sb**2
     &      +EE**2/CW**2*Zh22**2*Zss11**2/cb**2*sb**2+2D0
     &      *EE**2/CW**2*Zh22**2*Zss12**2/cb**2*sb**2-6D0
     &      *EE**2/MW**2*Mq**2/SW**2*Zh21**2/cb**2-3D0*EE
     &      **2/SW**2*Zh21**2*Zss11**2/cb**2*sb**2+3D0*EE
     &      **2/SW**2*Zh22**2*Zss11**2/cb**2*sb**2+12D0*EE
     &      *Lambda/MW*Mq/SW*Sqrt2*Zh22*Zh23*Zss11*Zss12/cb
      AAABR(1444) = EE**2/CW**2*Zh21*Zh31*Zss11**2/cb**2
     &      +2D0*EE**2/CW**2*Zh21*Zh31*Zss12**2/cb**2-EE
     &      **2/CW**2*Zh22*Zh32*Zss11**2/cb**2-2D0*EE**2
     &      /CW**2*Zh22*Zh32*Zss12**2/cb**2+3D0*EE**2/SW
     &      **2*Zh21*Zh31*Zss11**2/cb**2-3D0*EE**2/SW**2
     &      *Zh22*Zh32*Zss11**2/cb**2-EE**2/CW**2*Zh21*Zh31
     &      *Zss11**2/cb**2*sb**2-2D0*EE**2/CW**2*Zh21*Zh31
     &      *Zss12**2/cb**2*sb**2+EE**2/CW**2*Zh22*Zh32*Zss11
     &      **2/cb**2*sb**2+2D0*EE**2/CW**2*Zh22*Zh32*Zss12
     &      **2/cb**2*sb**2-6D0*EE**2/MW**2*Mq**2/SW**2*Zh21
     &      *Zh31/cb**2-3D0*EE**2/SW**2*Zh21*Zh31*Zss11**2
     &      /cb**2*sb**2+3D0*EE**2/SW**2*Zh22*Zh32*Zss11
     &      **2/cb**2*sb**2+6D0*EE*Lambda/MW*Mq/SW*Sqrt2
     &      *Zh22*Zh33*Zss11*Zss12/cb+6D0*EE*Lambda/MW*Mq
     &      /SW*Sqrt2*Zh23*Zh32*Zss11*Zss12/cb
      AAABR(1445) = EE**2/CW**2*Zb11**2*Zh31**2/cb**2-EE
     &      **2/CW**2*Zb11**2*Zh32**2/cb**2+2D0*EE**2/CW
     &      **2*Zb12**2*Zh31**2/cb**2-2D0*EE**2/CW**2*Zb12
     &      **2*Zh32**2/cb**2+3D0*EE**2/SW**2*Zb11**2*Zh31
     &      **2/cb**2-3D0*EE**2/SW**2*Zb11**2*Zh32**2/cb
     &      **2-EE**2/CW**2*Zb11**2*Zh31**2/cb**2*sb**2+EE
     &      **2/CW**2*Zb11**2*Zh32**2/cb**2*sb**2-2D0*EE
     &      **2/CW**2*Zb12**2*Zh31**2/cb**2*sb**2+2D0*EE
     &      **2/CW**2*Zb12**2*Zh32**2/cb**2*sb**2-6D0*EE
     &      **2/MW**2*MbMM**2/SW**2*Zh31**2/cb**2-3D0*EE
     &      **2/SW**2*Zb11**2*Zh31**2/cb**2*sb**2+3D0*EE
     &      **2/SW**2*Zb11**2*Zh32**2/cb**2*sb**2+12D0*EE
     &      *Lambda/MW*MbMM/SW*Sqrt2*Zb11*Zb12*Zh32*Zh33/cb
      AAABR(1446) = EE**2/CW**2*Zdd11**2*Zh11*Zh31/cb**2
     &      -EE**2/CW**2*Zdd11**2*Zh12*Zh32/cb**2+2D0*EE
     &      **2/CW**2*Zdd12**2*Zh11*Zh31/cb**2-2D0*EE**2
     &      /CW**2*Zdd12**2*Zh12*Zh32/cb**2+3D0*EE**2/SW
     &      **2*Zdd11**2*Zh11*Zh31/cb**2-3D0*EE**2/SW**2
     &      *Zdd11**2*Zh12*Zh32/cb**2-EE**2/CW**2*Zdd11**2
     &      *Zh11*Zh31/cb**2*sb**2+EE**2/CW**2*Zdd11**2*Zh12
     &      *Zh32/cb**2*sb**2-2D0*EE**2/CW**2*Zdd12**2*Zh11
     &      *Zh31/cb**2*sb**2+2D0*EE**2/CW**2*Zdd12**2*Zh12
     &      *Zh32/cb**2*sb**2-6D0*EE**2/MW**2*Mq**2/SW**2
     &      *Zh11*Zh31/cb**2-3D0*EE**2/SW**2*Zdd11**2*Zh11
     &      *Zh31/cb**2*sb**2+3D0*EE**2/SW**2*Zdd11**2*Zh12
     &      *Zh32/cb**2*sb**2+6D0*EE*Lambda/MW*Mq/SW*Sqrt2
     &      *Zdd11*Zdd12*Zh12*Zh33/cb+6D0*EE*Lambda/MW*Mq
     &      /SW*Sqrt2*Zdd11*Zdd12*Zh13*Zh32/cb
      AAABR(1447) = EE**2/CW**2*Zdd11**2*Zh11**2/cb**2-EE
     &      **2/CW**2*Zdd11**2*Zh12**2/cb**2+2D0*EE**2/CW
     &      **2*Zdd12**2*Zh11**2/cb**2-2D0*EE**2/CW**2*Zdd12
     &      **2*Zh12**2/cb**2+3D0*EE**2/SW**2*Zdd11**2*Zh11
     &      **2/cb**2-3D0*EE**2/SW**2*Zdd11**2*Zh12**2/cb
     &      **2-EE**2/CW**2*Zdd11**2*Zh11**2/cb**2*sb**2
     &      +EE**2/CW**2*Zdd11**2*Zh12**2/cb**2*sb**2-2D0
     &      *EE**2/CW**2*Zdd12**2*Zh11**2/cb**2*sb**2+2D0
     &      *EE**2/CW**2*Zdd12**2*Zh12**2/cb**2*sb**2-6D0
     &      *EE**2/MW**2*Mq**2/SW**2*Zh11**2/cb**2-3D0*EE
     &      **2/SW**2*Zdd11**2*Zh11**2/cb**2*sb**2+3D0*EE
     &      **2/SW**2*Zdd11**2*Zh12**2/cb**2*sb**2+12D0*EE
     &      *Lambda/MW*Mq/SW*Sqrt2*Zdd11*Zdd12*Zh12*Zh13/cb
      AAABR(1448) = EE**2/CW**2*Zdd11**2*Zh21**2/cb**2-EE
     &      **2/CW**2*Zdd11**2*Zh22**2/cb**2+2D0*EE**2/CW
     &      **2*Zdd12**2*Zh21**2/cb**2-2D0*EE**2/CW**2*Zdd12
     &      **2*Zh22**2/cb**2+3D0*EE**2/SW**2*Zdd11**2*Zh21
     &      **2/cb**2-3D0*EE**2/SW**2*Zdd11**2*Zh22**2/cb
     &      **2-EE**2/CW**2*Zdd11**2*Zh21**2/cb**2*sb**2
     &      +EE**2/CW**2*Zdd11**2*Zh22**2/cb**2*sb**2-2D0
     &      *EE**2/CW**2*Zdd12**2*Zh21**2/cb**2*sb**2+2D0
     &      *EE**2/CW**2*Zdd12**2*Zh22**2/cb**2*sb**2-6D0
     &      *EE**2/MW**2*Mq**2/SW**2*Zh21**2/cb**2-3D0*EE
     &      **2/SW**2*Zdd11**2*Zh21**2/cb**2*sb**2+3D0*EE
     &      **2/SW**2*Zdd11**2*Zh22**2/cb**2*sb**2+12D0*EE
     &      *Lambda/MW*Mq/SW*Sqrt2*Zdd11*Zdd12*Zh22*Zh23/cb
      AAABR(1449) = EE**2/CW**2*Zdd11**2*Zh21*Zh31/cb**2
     &      -EE**2/CW**2*Zdd11**2*Zh22*Zh32/cb**2+2D0*EE
     &      **2/CW**2*Zdd12**2*Zh21*Zh31/cb**2-2D0*EE**2
     &      /CW**2*Zdd12**2*Zh22*Zh32/cb**2+3D0*EE**2/SW
     &      **2*Zdd11**2*Zh21*Zh31/cb**2-3D0*EE**2/SW**2
     &      *Zdd11**2*Zh22*Zh32/cb**2-EE**2/CW**2*Zdd11**2
     &      *Zh21*Zh31/cb**2*sb**2+EE**2/CW**2*Zdd11**2*Zh22
     &      *Zh32/cb**2*sb**2-2D0*EE**2/CW**2*Zdd12**2*Zh21
     &      *Zh31/cb**2*sb**2+2D0*EE**2/CW**2*Zdd12**2*Zh22
     &      *Zh32/cb**2*sb**2-6D0*EE**2/MW**2*Mq**2/SW**2
     &      *Zh21*Zh31/cb**2-3D0*EE**2/SW**2*Zdd11**2*Zh21
     &      *Zh31/cb**2*sb**2+3D0*EE**2/SW**2*Zdd11**2*Zh22
     &      *Zh32/cb**2*sb**2+6D0*EE*Lambda/MW*Mq/SW*Sqrt2
     &      *Zdd11*Zdd12*Zh22*Zh33/cb+6D0*EE*Lambda/MW*Mq
     &      /SW*Sqrt2*Zdd11*Zdd12*Zh23*Zh32/cb
      AAABR(1450) = EE**2/CW**2*Zb11**2*Zh11*Zh21/cb**2-EE
     &      **2/CW**2*Zb11**2*Zh12*Zh22/cb**2+2D0*EE**2/CW
     &      **2*Zb12**2*Zh11*Zh21/cb**2-2D0*EE**2/CW**2*Zb12
     &      **2*Zh12*Zh22/cb**2+3D0*EE**2/SW**2*Zb11**2*Zh11
     &      *Zh21/cb**2-3D0*EE**2/SW**2*Zb11**2*Zh12*Zh22
     &      /cb**2-EE**2/CW**2*Zb11**2*Zh11*Zh21/cb**2*sb
     &      **2+EE**2/CW**2*Zb11**2*Zh12*Zh22/cb**2*sb**2
     &      -2D0*EE**2/CW**2*Zb12**2*Zh11*Zh21/cb**2*sb**2
     &      +2D0*EE**2/CW**2*Zb12**2*Zh12*Zh22/cb**2*sb**2
     &      -6D0*EE**2/MW**2*MbMM**2/SW**2*Zh11*Zh21/cb**2
     &      -3D0*EE**2/SW**2*Zb11**2*Zh11*Zh21/cb**2*sb**2
     &      +3D0*EE**2/SW**2*Zb11**2*Zh12*Zh22/cb**2*sb**2
     &      +6D0*EE*Lambda/MW*MbMM/SW*Sqrt2*Zb11*Zb12*Zh12
     &      *Zh23/cb+6D0*EE*Lambda/MW*MbMM/SW*Sqrt2*Zb11
     &      *Zb12*Zh13*Zh22/cb
      AAABR(1451) = EE**2/CW**2*Zdd11**2*Zh31**2/cb**2-EE
     &      **2/CW**2*Zdd11**2*Zh32**2/cb**2+2D0*EE**2/CW
     &      **2*Zdd12**2*Zh31**2/cb**2-2D0*EE**2/CW**2*Zdd12
     &      **2*Zh32**2/cb**2+3D0*EE**2/SW**2*Zdd11**2*Zh31
     &      **2/cb**2-3D0*EE**2/SW**2*Zdd11**2*Zh32**2/cb
     &      **2-EE**2/CW**2*Zdd11**2*Zh31**2/cb**2*sb**2
     &      +EE**2/CW**2*Zdd11**2*Zh32**2/cb**2*sb**2-2D0
     &      *EE**2/CW**2*Zdd12**2*Zh31**2/cb**2*sb**2+2D0
     &      *EE**2/CW**2*Zdd12**2*Zh32**2/cb**2*sb**2-6D0
     &      *EE**2/MW**2*Mq**2/SW**2*Zh31**2/cb**2-3D0*EE
     &      **2/SW**2*Zdd11**2*Zh31**2/cb**2*sb**2+3D0*EE
     &      **2/SW**2*Zdd11**2*Zh32**2/cb**2*sb**2+12D0*EE
     &      *Lambda/MW*Mq/SW*Sqrt2*Zdd11*Zdd12*Zh32*Zh33/cb
      AAABR(1452) = EE**2/CW**2*Zb11**2*Zh21*Zh31/cb**2-EE
     &      **2/CW**2*Zb11**2*Zh22*Zh32/cb**2+2D0*EE**2/CW
     &      **2*Zb12**2*Zh21*Zh31/cb**2-2D0*EE**2/CW**2*Zb12
     &      **2*Zh22*Zh32/cb**2+3D0*EE**2/SW**2*Zb11**2*Zh21
     &      *Zh31/cb**2-3D0*EE**2/SW**2*Zb11**2*Zh22*Zh32
     &      /cb**2-EE**2/CW**2*Zb11**2*Zh21*Zh31/cb**2*sb
     &      **2+EE**2/CW**2*Zb11**2*Zh22*Zh32/cb**2*sb**2
     &      -2D0*EE**2/CW**2*Zb12**2*Zh21*Zh31/cb**2*sb**2
     &      +2D0*EE**2/CW**2*Zb12**2*Zh22*Zh32/cb**2*sb**2
     &      -6D0*EE**2/MW**2*MbMM**2/SW**2*Zh21*Zh31/cb**2
     &      -3D0*EE**2/SW**2*Zb11**2*Zh21*Zh31/cb**2*sb**2
     &      +3D0*EE**2/SW**2*Zb11**2*Zh22*Zh32/cb**2*sb**2
     &      +6D0*EE*Lambda/MW*MbMM/SW*Sqrt2*Zb11*Zb12*Zh22
     &      *Zh33/cb+6D0*EE*Lambda/MW*MbMM/SW*Sqrt2*Zb11
     &      *Zb12*Zh23*Zh32/cb
      AAABR(1453) = EE**2/CW**2*Zh31**2*Zss11**2/cb**2+2D0
     &      *EE**2/CW**2*Zh31**2*Zss12**2/cb**2-EE**2/CW
     &      **2*Zh32**2*Zss11**2/cb**2-2D0*EE**2/CW**2*Zh32
     &      **2*Zss12**2/cb**2+3D0*EE**2/SW**2*Zh31**2*Zss11
     &      **2/cb**2-3D0*EE**2/SW**2*Zh32**2*Zss11**2/cb
     &      **2-EE**2/CW**2*Zh31**2*Zss11**2/cb**2*sb**2
     &      -2D0*EE**2/CW**2*Zh31**2*Zss12**2/cb**2*sb**2
     &      +EE**2/CW**2*Zh32**2*Zss11**2/cb**2*sb**2+2D0
     &      *EE**2/CW**2*Zh32**2*Zss12**2/cb**2*sb**2-6D0
     &      *EE**2/MW**2*Mq**2/SW**2*Zh31**2/cb**2-3D0*EE
     &      **2/SW**2*Zh31**2*Zss11**2/cb**2*sb**2+3D0*EE
     &      **2/SW**2*Zh32**2*Zss11**2/cb**2*sb**2+12D0*EE
     &      *Lambda/MW*Mq/SW*Sqrt2*Zh32*Zh33*Zss11*Zss12/cb
      AAABR(1454) = EE**2/CW**2*Zh11*Zh21*Zss11**2/cb**2
     &      +2D0*EE**2/CW**2*Zh11*Zh21*Zss12**2/cb**2-EE
     &      **2/CW**2*Zh12*Zh22*Zss11**2/cb**2-2D0*EE**2
     &      /CW**2*Zh12*Zh22*Zss12**2/cb**2+3D0*EE**2/SW
     &      **2*Zh11*Zh21*Zss11**2/cb**2-3D0*EE**2/SW**2
     &      *Zh12*Zh22*Zss11**2/cb**2-EE**2/CW**2*Zh11*Zh21
     &      *Zss11**2/cb**2*sb**2-2D0*EE**2/CW**2*Zh11*Zh21
     &      *Zss12**2/cb**2*sb**2+EE**2/CW**2*Zh12*Zh22*Zss11
     &      **2/cb**2*sb**2+2D0*EE**2/CW**2*Zh12*Zh22*Zss12
     &      **2/cb**2*sb**2-6D0*EE**2/MW**2*Mq**2/SW**2*Zh11
     &      *Zh21/cb**2-3D0*EE**2/SW**2*Zh11*Zh21*Zss11**2
     &      /cb**2*sb**2+3D0*EE**2/SW**2*Zh12*Zh22*Zss11
     &      **2/cb**2*sb**2+6D0*EE*Lambda/MW*Mq/SW*Sqrt2
     &      *Zh12*Zh23*Zss11*Zss12/cb+6D0*EE*Lambda/MW*Mq
     &      /SW*Sqrt2*Zh13*Zh22*Zss11*Zss12/cb
      AAABR(1455) = EE**2/CW**2*Zb11**2*Zh11**2/cb**2-EE
     &      **2/CW**2*Zb11**2*Zh12**2/cb**2+2D0*EE**2/CW
     &      **2*Zb12**2*Zh11**2/cb**2-2D0*EE**2/CW**2*Zb12
     &      **2*Zh12**2/cb**2+3D0*EE**2/SW**2*Zb11**2*Zh11
     &      **2/cb**2-3D0*EE**2/SW**2*Zb11**2*Zh12**2/cb
     &      **2-EE**2/CW**2*Zb11**2*Zh11**2/cb**2*sb**2+EE
     &      **2/CW**2*Zb11**2*Zh12**2/cb**2*sb**2-2D0*EE
     &      **2/CW**2*Zb12**2*Zh11**2/cb**2*sb**2+2D0*EE
     &      **2/CW**2*Zb12**2*Zh12**2/cb**2*sb**2-6D0*EE
     &      **2/MW**2*MbMM**2/SW**2*Zh11**2/cb**2-3D0*EE
     &      **2/SW**2*Zb11**2*Zh11**2/cb**2*sb**2+3D0*EE
     &      **2/SW**2*Zb11**2*Zh12**2/cb**2*sb**2+12D0*EE
     &      *Lambda/MW*MbMM/SW*Sqrt2*Zb11*Zb12*Zh12*Zh13/cb
      AAABR(1456) = EE**2/CW**2*Zdd11**2*Zh11*Zh21/cb**2
     &      -EE**2/CW**2*Zdd11**2*Zh12*Zh22/cb**2+2D0*EE
     &      **2/CW**2*Zdd12**2*Zh11*Zh21/cb**2-2D0*EE**2
     &      /CW**2*Zdd12**2*Zh12*Zh22/cb**2+3D0*EE**2/SW
     &      **2*Zdd11**2*Zh11*Zh21/cb**2-3D0*EE**2/SW**2
     &      *Zdd11**2*Zh12*Zh22/cb**2-EE**2/CW**2*Zdd11**2
     &      *Zh11*Zh21/cb**2*sb**2+EE**2/CW**2*Zdd11**2*Zh12
     &      *Zh22/cb**2*sb**2-2D0*EE**2/CW**2*Zdd12**2*Zh11
     &      *Zh21/cb**2*sb**2+2D0*EE**2/CW**2*Zdd12**2*Zh12
     &      *Zh22/cb**2*sb**2-6D0*EE**2/MW**2*Mq**2/SW**2
     &      *Zh11*Zh21/cb**2-3D0*EE**2/SW**2*Zdd11**2*Zh11
     &      *Zh21/cb**2*sb**2+3D0*EE**2/SW**2*Zdd11**2*Zh12
     &      *Zh22/cb**2*sb**2+6D0*EE*Lambda/MW*Mq/SW*Sqrt2
     &      *Zdd11*Zdd12*Zh12*Zh23/cb+6D0*EE*Lambda/MW*Mq
     &      /SW*Sqrt2*Zdd11*Zdd12*Zh13*Zh22/cb
      AAABR(1457) = EE**2/CW**2*Zb11*Zb21*Zh11**2/cb**2-EE
     &      **2/CW**2*Zb11*Zb21*Zh12**2/cb**2+2D0*EE**2/CW
     &      **2*Zb12*Zb22*Zh11**2/cb**2-2D0*EE**2/CW**2*Zb12
     &      *Zb22*Zh12**2/cb**2+3D0*EE**2/SW**2*Zb11*Zb21
     &      *Zh11**2/cb**2-3D0*EE**2/SW**2*Zb11*Zb21*Zh12
     &      **2/cb**2-EE**2/CW**2*Zb11*Zb21*Zh11**2/cb**2
     &      *sb**2+EE**2/CW**2*Zb11*Zb21*Zh12**2/cb**2*sb
     &      **2-2D0*EE**2/CW**2*Zb12*Zb22*Zh11**2/cb**2*sb
     &      **2+2D0*EE**2/CW**2*Zb12*Zb22*Zh12**2/cb**2*sb
     &      **2-3D0*EE**2/SW**2*Zb11*Zb21*Zh11**2/cb**2*sb
     &      **2+3D0*EE**2/SW**2*Zb11*Zb21*Zh12**2/cb**2*sb
     &      **2+6D0*EE*Lambda/MW*MbMM/SW*Sqrt2*Zb11*Zb22
     &      *Zh12*Zh13/cb+6D0*EE*Lambda/MW*MbMM/SW*Sqrt2
     &      *Zb12*Zb21*Zh12*Zh13/cb
      AAABR(1458) = EE**2/CW**2*Zh11**2*Zss11*Zss21/cb**2
     &      +2D0*EE**2/CW**2*Zh11**2*Zss12*Zss22/cb**2-EE
     &      **2/CW**2*Zh12**2*Zss11*Zss21/cb**2-2D0*EE**2
     &      /CW**2*Zh12**2*Zss12*Zss22/cb**2+3D0*EE**2/SW
     &      **2*Zh11**2*Zss11*Zss21/cb**2-3D0*EE**2/SW**2
     &      *Zh12**2*Zss11*Zss21/cb**2-EE**2/CW**2*Zh11**2
     &      *Zss11*Zss21/cb**2*sb**2-2D0*EE**2/CW**2*Zh11
     &      **2*Zss12*Zss22/cb**2*sb**2+EE**2/CW**2*Zh12
     &      **2*Zss11*Zss21/cb**2*sb**2+2D0*EE**2/CW**2*Zh12
     &      **2*Zss12*Zss22/cb**2*sb**2-3D0*EE**2/SW**2*Zh11
     &      **2*Zss11*Zss21/cb**2*sb**2+3D0*EE**2/SW**2*Zh12
     &      **2*Zss11*Zss21/cb**2*sb**2+6D0*EE*Lambda/MW
     &      *Mq/SW*Sqrt2*Zh12*Zh13*Zss11*Zss22/cb+6D0*EE
     &      *Lambda/MW*Mq/SW*Sqrt2*Zh12*Zh13*Zss12*Zss21/cb
      AAABR(1459) = EE**2/CW**2*Zdd11*Zdd21*Zh21*Zh31/cb
     &      **2-EE**2/CW**2*Zdd11*Zdd21*Zh22*Zh32/cb**2+2D0
     &      *EE**2/CW**2*Zdd12*Zdd22*Zh21*Zh31/cb**2-2D0
     &      *EE**2/CW**2*Zdd12*Zdd22*Zh22*Zh32/cb**2+3D0
     &      *EE**2/SW**2*Zdd11*Zdd21*Zh21*Zh31/cb**2-3D0
     &      *EE**2/SW**2*Zdd11*Zdd21*Zh22*Zh32/cb**2-EE**2
     &      /CW**2*Zdd11*Zdd21*Zh21*Zh31/cb**2*sb**2+EE**2
     &      /CW**2*Zdd11*Zdd21*Zh22*Zh32/cb**2*sb**2-2D0
     &      *EE**2/CW**2*Zdd12*Zdd22*Zh21*Zh31/cb**2*sb**2
     &      +2D0*EE**2/CW**2*Zdd12*Zdd22*Zh22*Zh32/cb**2
     &      *sb**2-3D0*EE**2/SW**2*Zdd11*Zdd21*Zh21*Zh31
     &      /cb**2*sb**2+3D0*EE**2/SW**2*Zdd11*Zdd21*Zh22
     &      *Zh32/cb**2*sb**2+3D0*EE*Lambda/MW*Mq/SW*Sqrt2
     &      *Zdd11*Zdd22*Zh22*Zh33/cb+3D0*EE*Lambda/MW*Mq
     &      /SW*Sqrt2*Zdd11*Zdd22*Zh23*Zh32/cb+3D0*EE*Lambda
     &      /MW*Mq/SW*Sqrt2*Zdd12*Zdd21*Zh22*Zh33/cb+3D0
     &      *EE*Lambda/MW*Mq/SW*Sqrt2*Zdd12*Zdd21*Zh23*Zh32
     &      /cb
      AAABR(1460) = EE**2/CW**2*Zdd11*Zdd21*Zh31**2/cb**2
     &      -EE**2/CW**2*Zdd11*Zdd21*Zh32**2/cb**2+2D0*EE
     &      **2/CW**2*Zdd12*Zdd22*Zh31**2/cb**2-2D0*EE**2
     &      /CW**2*Zdd12*Zdd22*Zh32**2/cb**2+3D0*EE**2/SW
     &      **2*Zdd11*Zdd21*Zh31**2/cb**2-3D0*EE**2/SW**2
     &      *Zdd11*Zdd21*Zh32**2/cb**2-EE**2/CW**2*Zdd11
     &      *Zdd21*Zh31**2/cb**2*sb**2+EE**2/CW**2*Zdd11
     &      *Zdd21*Zh32**2/cb**2*sb**2-2D0*EE**2/CW**2*Zdd12
     &      *Zdd22*Zh31**2/cb**2*sb**2+2D0*EE**2/CW**2*Zdd12
     &      *Zdd22*Zh32**2/cb**2*sb**2-3D0*EE**2/SW**2*Zdd11
     &      *Zdd21*Zh31**2/cb**2*sb**2+3D0*EE**2/SW**2*Zdd11
     &      *Zdd21*Zh32**2/cb**2*sb**2+6D0*EE*Lambda/MW*Mq
     &      /SW*Sqrt2*Zdd11*Zdd22*Zh32*Zh33/cb+6D0*EE*Lambda
     &      /MW*Mq/SW*Sqrt2*Zdd12*Zdd21*Zh32*Zh33/cb
      AAABR(1461) = EE**2/CW**2*Zh21*Zh31*Zss11*Zss21/cb
     &      **2+2D0*EE**2/CW**2*Zh21*Zh31*Zss12*Zss22/cb
     &      **2-EE**2/CW**2*Zh22*Zh32*Zss11*Zss21/cb**2-2D0
     &      *EE**2/CW**2*Zh22*Zh32*Zss12*Zss22/cb**2+3D0
     &      *EE**2/SW**2*Zh21*Zh31*Zss11*Zss21/cb**2-3D0
     &      *EE**2/SW**2*Zh22*Zh32*Zss11*Zss21/cb**2-EE**2
     &      /CW**2*Zh21*Zh31*Zss11*Zss21/cb**2*sb**2-2D0
     &      *EE**2/CW**2*Zh21*Zh31*Zss12*Zss22/cb**2*sb**2
     &      +EE**2/CW**2*Zh22*Zh32*Zss11*Zss21/cb**2*sb**2
     &      +2D0*EE**2/CW**2*Zh22*Zh32*Zss12*Zss22/cb**2
     &      *sb**2-3D0*EE**2/SW**2*Zh21*Zh31*Zss11*Zss21
     &      /cb**2*sb**2+3D0*EE**2/SW**2*Zh22*Zh32*Zss11
     &      *Zss21/cb**2*sb**2+3D0*EE*Lambda/MW*Mq/SW*Sqrt2
     &      *Zh22*Zh33*Zss11*Zss22/cb+3D0*EE*Lambda/MW*Mq
     &      /SW*Sqrt2*Zh22*Zh33*Zss12*Zss21/cb+3D0*EE*Lambda
     &      /MW*Mq/SW*Sqrt2*Zh23*Zh32*Zss11*Zss22/cb+3D0
     &      *EE*Lambda/MW*Mq/SW*Sqrt2*Zh23*Zh32*Zss12*Zss21
     &      /cb
      AAABR(1462) = EE**2/CW**2*Zh21**2*Zss11*Zss21/cb**2
     &      +2D0*EE**2/CW**2*Zh21**2*Zss12*Zss22/cb**2-EE
     &      **2/CW**2*Zh22**2*Zss11*Zss21/cb**2-2D0*EE**2
     &      /CW**2*Zh22**2*Zss12*Zss22/cb**2+3D0*EE**2/SW
     &      **2*Zh21**2*Zss11*Zss21/cb**2-3D0*EE**2/SW**2
     &      *Zh22**2*Zss11*Zss21/cb**2-EE**2/CW**2*Zh21**2
     &      *Zss11*Zss21/cb**2*sb**2-2D0*EE**2/CW**2*Zh21
     &      **2*Zss12*Zss22/cb**2*sb**2+EE**2/CW**2*Zh22
     &      **2*Zss11*Zss21/cb**2*sb**2+2D0*EE**2/CW**2*Zh22
     &      **2*Zss12*Zss22/cb**2*sb**2-3D0*EE**2/SW**2*Zh21
     &      **2*Zss11*Zss21/cb**2*sb**2+3D0*EE**2/SW**2*Zh22
     &      **2*Zss11*Zss21/cb**2*sb**2+6D0*EE*Lambda/MW
     &      *Mq/SW*Sqrt2*Zh22*Zh23*Zss11*Zss22/cb+6D0*EE
     &      *Lambda/MW*Mq/SW*Sqrt2*Zh22*Zh23*Zss12*Zss21/cb
      AAABR(1463) = EE**2/CW**2*Zb11*Zb21*Zh21*Zh31/cb**2
     &      -EE**2/CW**2*Zb11*Zb21*Zh22*Zh32/cb**2+2D0*EE
     &      **2/CW**2*Zb12*Zb22*Zh21*Zh31/cb**2-2D0*EE**2
     &      /CW**2*Zb12*Zb22*Zh22*Zh32/cb**2+3D0*EE**2/SW
     &      **2*Zb11*Zb21*Zh21*Zh31/cb**2-3D0*EE**2/SW**2
     &      *Zb11*Zb21*Zh22*Zh32/cb**2-EE**2/CW**2*Zb11*Zb21
     &      *Zh21*Zh31/cb**2*sb**2+EE**2/CW**2*Zb11*Zb21
     &      *Zh22*Zh32/cb**2*sb**2-2D0*EE**2/CW**2*Zb12*Zb22
     &      *Zh21*Zh31/cb**2*sb**2+2D0*EE**2/CW**2*Zb12*Zb22
     &      *Zh22*Zh32/cb**2*sb**2-3D0*EE**2/SW**2*Zb11*Zb21
     &      *Zh21*Zh31/cb**2*sb**2+3D0*EE**2/SW**2*Zb11*Zb21
     &      *Zh22*Zh32/cb**2*sb**2+3D0*EE*Lambda/MW*MbMM
     &      /SW*Sqrt2*Zb11*Zb22*Zh22*Zh33/cb+3D0*EE*Lambda
     &      /MW*MbMM/SW*Sqrt2*Zb11*Zb22*Zh23*Zh32/cb+3D0
     &      *EE*Lambda/MW*MbMM/SW*Sqrt2*Zb12*Zb21*Zh22*Zh33
     &      /cb+3D0*EE*Lambda/MW*MbMM/SW*Sqrt2*Zb12*Zb21
     &      *Zh23*Zh32/cb
      AAABR(1464) = EE**2/CW**2*Zb11*Zb21*Zh21**2/cb**2-EE
     &      **2/CW**2*Zb11*Zb21*Zh22**2/cb**2+2D0*EE**2/CW
     &      **2*Zb12*Zb22*Zh21**2/cb**2-2D0*EE**2/CW**2*Zb12
     &      *Zb22*Zh22**2/cb**2+3D0*EE**2/SW**2*Zb11*Zb21
     &      *Zh21**2/cb**2-3D0*EE**2/SW**2*Zb11*Zb21*Zh22
     &      **2/cb**2-EE**2/CW**2*Zb11*Zb21*Zh21**2/cb**2
     &      *sb**2+EE**2/CW**2*Zb11*Zb21*Zh22**2/cb**2*sb
     &      **2-2D0*EE**2/CW**2*Zb12*Zb22*Zh21**2/cb**2*sb
     &      **2+2D0*EE**2/CW**2*Zb12*Zb22*Zh22**2/cb**2*sb
     &      **2-3D0*EE**2/SW**2*Zb11*Zb21*Zh21**2/cb**2*sb
     &      **2+3D0*EE**2/SW**2*Zb11*Zb21*Zh22**2/cb**2*sb
     &      **2+6D0*EE*Lambda/MW*MbMM/SW*Sqrt2*Zb11*Zb22
     &      *Zh22*Zh23/cb+6D0*EE*Lambda/MW*MbMM/SW*Sqrt2
     &      *Zb12*Zb21*Zh22*Zh23/cb
      AAABR(1465) = EE**2/CW**2*Zdd11*Zdd21*Zh11*Zh21/cb
     &      **2-EE**2/CW**2*Zdd11*Zdd21*Zh12*Zh22/cb**2+2D0
     &      *EE**2/CW**2*Zdd12*Zdd22*Zh11*Zh21/cb**2-2D0
     &      *EE**2/CW**2*Zdd12*Zdd22*Zh12*Zh22/cb**2+3D0
     &      *EE**2/SW**2*Zdd11*Zdd21*Zh11*Zh21/cb**2-3D0
     &      *EE**2/SW**2*Zdd11*Zdd21*Zh12*Zh22/cb**2-EE**2
     &      /CW**2*Zdd11*Zdd21*Zh11*Zh21/cb**2*sb**2+EE**2
     &      /CW**2*Zdd11*Zdd21*Zh12*Zh22/cb**2*sb**2-2D0
     &      *EE**2/CW**2*Zdd12*Zdd22*Zh11*Zh21/cb**2*sb**2
     &      +2D0*EE**2/CW**2*Zdd12*Zdd22*Zh12*Zh22/cb**2
     &      *sb**2-3D0*EE**2/SW**2*Zdd11*Zdd21*Zh11*Zh21
     &      /cb**2*sb**2+3D0*EE**2/SW**2*Zdd11*Zdd21*Zh12
     &      *Zh22/cb**2*sb**2+3D0*EE*Lambda/MW*Mq/SW*Sqrt2
     &      *Zdd11*Zdd22*Zh12*Zh23/cb+3D0*EE*Lambda/MW*Mq
     &      /SW*Sqrt2*Zdd11*Zdd22*Zh13*Zh22/cb+3D0*EE*Lambda
     &      /MW*Mq/SW*Sqrt2*Zdd12*Zdd21*Zh12*Zh23/cb+3D0
     &      *EE*Lambda/MW*Mq/SW*Sqrt2*Zdd12*Zdd21*Zh13*Zh22
     &      /cb
      AAABR(1466) = EE**2/CW**2*Zb11*Zb21*Zh11*Zh21/cb**2
     &      -EE**2/CW**2*Zb11*Zb21*Zh12*Zh22/cb**2+2D0*EE
     &      **2/CW**2*Zb12*Zb22*Zh11*Zh21/cb**2-2D0*EE**2
     &      /CW**2*Zb12*Zb22*Zh12*Zh22/cb**2+3D0*EE**2/SW
     &      **2*Zb11*Zb21*Zh11*Zh21/cb**2-3D0*EE**2/SW**2
     &      *Zb11*Zb21*Zh12*Zh22/cb**2-EE**2/CW**2*Zb11*Zb21
     &      *Zh11*Zh21/cb**2*sb**2+EE**2/CW**2*Zb11*Zb21
     &      *Zh12*Zh22/cb**2*sb**2-2D0*EE**2/CW**2*Zb12*Zb22
     &      *Zh11*Zh21/cb**2*sb**2+2D0*EE**2/CW**2*Zb12*Zb22
     &      *Zh12*Zh22/cb**2*sb**2-3D0*EE**2/SW**2*Zb11*Zb21
     &      *Zh11*Zh21/cb**2*sb**2+3D0*EE**2/SW**2*Zb11*Zb21
     &      *Zh12*Zh22/cb**2*sb**2+3D0*EE*Lambda/MW*MbMM
     &      /SW*Sqrt2*Zb11*Zb22*Zh12*Zh23/cb+3D0*EE*Lambda
     &      /MW*MbMM/SW*Sqrt2*Zb11*Zb22*Zh13*Zh22/cb+3D0
     &      *EE*Lambda/MW*MbMM/SW*Sqrt2*Zb12*Zb21*Zh12*Zh23
     &      /cb+3D0*EE*Lambda/MW*MbMM/SW*Sqrt2*Zb12*Zb21
     &      *Zh13*Zh22/cb
      AAABR(1467) = EE**2/CW**2*Zdd11*Zdd21*Zh21**2/cb**2
     &      -EE**2/CW**2*Zdd11*Zdd21*Zh22**2/cb**2+2D0*EE
     &      **2/CW**2*Zdd12*Zdd22*Zh21**2/cb**2-2D0*EE**2
     &      /CW**2*Zdd12*Zdd22*Zh22**2/cb**2+3D0*EE**2/SW
     &      **2*Zdd11*Zdd21*Zh21**2/cb**2-3D0*EE**2/SW**2
     &      *Zdd11*Zdd21*Zh22**2/cb**2-EE**2/CW**2*Zdd11
     &      *Zdd21*Zh21**2/cb**2*sb**2+EE**2/CW**2*Zdd11
     &      *Zdd21*Zh22**2/cb**2*sb**2-2D0*EE**2/CW**2*Zdd12
     &      *Zdd22*Zh21**2/cb**2*sb**2+2D0*EE**2/CW**2*Zdd12
     &      *Zdd22*Zh22**2/cb**2*sb**2-3D0*EE**2/SW**2*Zdd11
     &      *Zdd21*Zh21**2/cb**2*sb**2+3D0*EE**2/SW**2*Zdd11
     &      *Zdd21*Zh22**2/cb**2*sb**2+6D0*EE*Lambda/MW*Mq
     &      /SW*Sqrt2*Zdd11*Zdd22*Zh22*Zh23/cb+6D0*EE*Lambda
     &      /MW*Mq/SW*Sqrt2*Zdd12*Zdd21*Zh22*Zh23/cb
      AAABR(1468) = EE**2/CW**2*Zh11*Zh31*Zss11*Zss21/cb
     &      **2+2D0*EE**2/CW**2*Zh11*Zh31*Zss12*Zss22/cb
     &      **2-EE**2/CW**2*Zh12*Zh32*Zss11*Zss21/cb**2-2D0
     &      *EE**2/CW**2*Zh12*Zh32*Zss12*Zss22/cb**2+3D0
     &      *EE**2/SW**2*Zh11*Zh31*Zss11*Zss21/cb**2-3D0
     &      *EE**2/SW**2*Zh12*Zh32*Zss11*Zss21/cb**2-EE**2
     &      /CW**2*Zh11*Zh31*Zss11*Zss21/cb**2*sb**2-2D0
     &      *EE**2/CW**2*Zh11*Zh31*Zss12*Zss22/cb**2*sb**2
     &      +EE**2/CW**2*Zh12*Zh32*Zss11*Zss21/cb**2*sb**2
     &      +2D0*EE**2/CW**2*Zh12*Zh32*Zss12*Zss22/cb**2
     &      *sb**2-3D0*EE**2/SW**2*Zh11*Zh31*Zss11*Zss21
     &      /cb**2*sb**2+3D0*EE**2/SW**2*Zh12*Zh32*Zss11
     &      *Zss21/cb**2*sb**2+3D0*EE*Lambda/MW*Mq/SW*Sqrt2
     &      *Zh12*Zh33*Zss11*Zss22/cb+3D0*EE*Lambda/MW*Mq
     &      /SW*Sqrt2*Zh12*Zh33*Zss12*Zss21/cb+3D0*EE*Lambda
     &      /MW*Mq/SW*Sqrt2*Zh13*Zh32*Zss11*Zss22/cb+3D0
     &      *EE*Lambda/MW*Mq/SW*Sqrt2*Zh13*Zh32*Zss12*Zss21
     &      /cb
      AAABR(1469) = EE**2/CW**2*Zh11*Zh21*Zss11*Zss21/cb
     &      **2+2D0*EE**2/CW**2*Zh11*Zh21*Zss12*Zss22/cb
     &      **2-EE**2/CW**2*Zh12*Zh22*Zss11*Zss21/cb**2-2D0
     &      *EE**2/CW**2*Zh12*Zh22*Zss12*Zss22/cb**2+3D0
     &      *EE**2/SW**2*Zh11*Zh21*Zss11*Zss21/cb**2-3D0
     &      *EE**2/SW**2*Zh12*Zh22*Zss11*Zss21/cb**2-EE**2
     &      /CW**2*Zh11*Zh21*Zss11*Zss21/cb**2*sb**2-2D0
     &      *EE**2/CW**2*Zh11*Zh21*Zss12*Zss22/cb**2*sb**2
     &      +EE**2/CW**2*Zh12*Zh22*Zss11*Zss21/cb**2*sb**2
     &      +2D0*EE**2/CW**2*Zh12*Zh22*Zss12*Zss22/cb**2
     &      *sb**2-3D0*EE**2/SW**2*Zh11*Zh21*Zss11*Zss21
     &      /cb**2*sb**2+3D0*EE**2/SW**2*Zh12*Zh22*Zss11
     &      *Zss21/cb**2*sb**2+3D0*EE*Lambda/MW*Mq/SW*Sqrt2
     &      *Zh12*Zh23*Zss11*Zss22/cb+3D0*EE*Lambda/MW*Mq
     &      /SW*Sqrt2*Zh12*Zh23*Zss12*Zss21/cb+3D0*EE*Lambda
     &      /MW*Mq/SW*Sqrt2*Zh13*Zh22*Zss11*Zss22/cb+3D0
     &      *EE*Lambda/MW*Mq/SW*Sqrt2*Zh13*Zh22*Zss12*Zss21
     &      /cb
      AAABR(1470) = EE**2/CW**2*Zb11*Zb21*Zh31**2/cb**2-EE
     &      **2/CW**2*Zb11*Zb21*Zh32**2/cb**2+2D0*EE**2/CW
     &      **2*Zb12*Zb22*Zh31**2/cb**2-2D0*EE**2/CW**2*Zb12
     &      *Zb22*Zh32**2/cb**2+3D0*EE**2/SW**2*Zb11*Zb21
     &      *Zh31**2/cb**2-3D0*EE**2/SW**2*Zb11*Zb21*Zh32
     &      **2/cb**2-EE**2/CW**2*Zb11*Zb21*Zh31**2/cb**2
     &      *sb**2+EE**2/CW**2*Zb11*Zb21*Zh32**2/cb**2*sb
     &      **2-2D0*EE**2/CW**2*Zb12*Zb22*Zh31**2/cb**2*sb
     &      **2+2D0*EE**2/CW**2*Zb12*Zb22*Zh32**2/cb**2*sb
     &      **2-3D0*EE**2/SW**2*Zb11*Zb21*Zh31**2/cb**2*sb
     &      **2+3D0*EE**2/SW**2*Zb11*Zb21*Zh32**2/cb**2*sb
     &      **2+6D0*EE*Lambda/MW*MbMM/SW*Sqrt2*Zb11*Zb22
     &      *Zh32*Zh33/cb+6D0*EE*Lambda/MW*MbMM/SW*Sqrt2
     &      *Zb12*Zb21*Zh32*Zh33/cb
      AAABR(1471) = EE**2/CW**2*Zdd11*Zdd21*Zh11*Zh31/cb
     &      **2-EE**2/CW**2*Zdd11*Zdd21*Zh12*Zh32/cb**2+2D0
     &      *EE**2/CW**2*Zdd12*Zdd22*Zh11*Zh31/cb**2-2D0
     &      *EE**2/CW**2*Zdd12*Zdd22*Zh12*Zh32/cb**2+3D0
     &      *EE**2/SW**2*Zdd11*Zdd21*Zh11*Zh31/cb**2-3D0
     &      *EE**2/SW**2*Zdd11*Zdd21*Zh12*Zh32/cb**2-EE**2
     &      /CW**2*Zdd11*Zdd21*Zh11*Zh31/cb**2*sb**2+EE**2
     &      /CW**2*Zdd11*Zdd21*Zh12*Zh32/cb**2*sb**2-2D0
     &      *EE**2/CW**2*Zdd12*Zdd22*Zh11*Zh31/cb**2*sb**2
     &      +2D0*EE**2/CW**2*Zdd12*Zdd22*Zh12*Zh32/cb**2
     &      *sb**2-3D0*EE**2/SW**2*Zdd11*Zdd21*Zh11*Zh31
     &      /cb**2*sb**2+3D0*EE**2/SW**2*Zdd11*Zdd21*Zh12
     &      *Zh32/cb**2*sb**2+3D0*EE*Lambda/MW*Mq/SW*Sqrt2
     &      *Zdd11*Zdd22*Zh12*Zh33/cb+3D0*EE*Lambda/MW*Mq
     &      /SW*Sqrt2*Zdd11*Zdd22*Zh13*Zh32/cb+3D0*EE*Lambda
     &      /MW*Mq/SW*Sqrt2*Zdd12*Zdd21*Zh12*Zh33/cb+3D0
     &      *EE*Lambda/MW*Mq/SW*Sqrt2*Zdd12*Zdd21*Zh13*Zh32
     &      /cb
      AAABR(1472) = EE**2/CW**2*Zb11*Zb21*Zh11*Zh31/cb**2
     &      -EE**2/CW**2*Zb11*Zb21*Zh12*Zh32/cb**2+2D0*EE
     &      **2/CW**2*Zb12*Zb22*Zh11*Zh31/cb**2-2D0*EE**2
     &      /CW**2*Zb12*Zb22*Zh12*Zh32/cb**2+3D0*EE**2/SW
     &      **2*Zb11*Zb21*Zh11*Zh31/cb**2-3D0*EE**2/SW**2
     &      *Zb11*Zb21*Zh12*Zh32/cb**2-EE**2/CW**2*Zb11*Zb21
     &      *Zh11*Zh31/cb**2*sb**2+EE**2/CW**2*Zb11*Zb21
     &      *Zh12*Zh32/cb**2*sb**2-2D0*EE**2/CW**2*Zb12*Zb22
     &      *Zh11*Zh31/cb**2*sb**2+2D0*EE**2/CW**2*Zb12*Zb22
     &      *Zh12*Zh32/cb**2*sb**2-3D0*EE**2/SW**2*Zb11*Zb21
     &      *Zh11*Zh31/cb**2*sb**2+3D0*EE**2/SW**2*Zb11*Zb21
     &      *Zh12*Zh32/cb**2*sb**2+3D0*EE*Lambda/MW*MbMM
     &      /SW*Sqrt2*Zb11*Zb22*Zh12*Zh33/cb+3D0*EE*Lambda
     &      /MW*MbMM/SW*Sqrt2*Zb11*Zb22*Zh13*Zh32/cb+3D0
     &      *EE*Lambda/MW*MbMM/SW*Sqrt2*Zb12*Zb21*Zh12*Zh33
     &      /cb+3D0*EE*Lambda/MW*MbMM/SW*Sqrt2*Zb12*Zb21
     &      *Zh13*Zh32/cb
      AAABR(1473) = EE**2/CW**2*Zh31**2*Zss11*Zss21/cb**2
     &      +2D0*EE**2/CW**2*Zh31**2*Zss12*Zss22/cb**2-EE
     &      **2/CW**2*Zh32**2*Zss11*Zss21/cb**2-2D0*EE**2
     &      /CW**2*Zh32**2*Zss12*Zss22/cb**2+3D0*EE**2/SW
     &      **2*Zh31**2*Zss11*Zss21/cb**2-3D0*EE**2/SW**2
     &      *Zh32**2*Zss11*Zss21/cb**2-EE**2/CW**2*Zh31**2
     &      *Zss11*Zss21/cb**2*sb**2-2D0*EE**2/CW**2*Zh31
     &      **2*Zss12*Zss22/cb**2*sb**2+EE**2/CW**2*Zh32
     &      **2*Zss11*Zss21/cb**2*sb**2+2D0*EE**2/CW**2*Zh32
     &      **2*Zss12*Zss22/cb**2*sb**2-3D0*EE**2/SW**2*Zh31
     &      **2*Zss11*Zss21/cb**2*sb**2+3D0*EE**2/SW**2*Zh32
     &      **2*Zss11*Zss21/cb**2*sb**2+6D0*EE*Lambda/MW
     &      *Mq/SW*Sqrt2*Zh32*Zh33*Zss11*Zss22/cb+6D0*EE
     &      *Lambda/MW*Mq/SW*Sqrt2*Zh32*Zh33*Zss12*Zss21/cb
      AAABR(1474) = EE**2/CW**2*Zdd11*Zdd21*Zh11**2/cb**2
     &      -EE**2/CW**2*Zdd11*Zdd21*Zh12**2/cb**2+2D0*EE
     &      **2/CW**2*Zdd12*Zdd22*Zh11**2/cb**2-2D0*EE**2
     &      /CW**2*Zdd12*Zdd22*Zh12**2/cb**2+3D0*EE**2/SW
     &      **2*Zdd11*Zdd21*Zh11**2/cb**2-3D0*EE**2/SW**2
     &      *Zdd11*Zdd21*Zh12**2/cb**2-EE**2/CW**2*Zdd11
     &      *Zdd21*Zh11**2/cb**2*sb**2+EE**2/CW**2*Zdd11
     &      *Zdd21*Zh12**2/cb**2*sb**2-2D0*EE**2/CW**2*Zdd12
     &      *Zdd22*Zh11**2/cb**2*sb**2+2D0*EE**2/CW**2*Zdd12
     &      *Zdd22*Zh12**2/cb**2*sb**2-3D0*EE**2/SW**2*Zdd11
     &      *Zdd21*Zh11**2/cb**2*sb**2+3D0*EE**2/SW**2*Zdd11
     &      *Zdd21*Zh12**2/cb**2*sb**2+6D0*EE*Lambda/MW*Mq
     &      /SW*Sqrt2*Zdd11*Zdd22*Zh12*Zh13/cb+6D0*EE*Lambda
     &      /MW*Mq/SW*Sqrt2*Zdd12*Zdd21*Zh12*Zh13/cb
      AAABR(1475) = GG**2/cb**2*(Zss21**4+Zss22**4+2D0*Zss21
     &      **2*Zss22**2-Zss21**4*sb**2-Zss22**4*sb**2-2D0
     &      *Zss21**2*Zss22**2*sb**2)
      AAABR(1476) = EE**2/cb**2*(Zss21**4/CW**2+4D0*Zss22
     &      **4/CW**2+9D0*Zss21**4/SW**2+4D0*Zss21**2/CW
     &      **2*Zss22**2-Zss21**4/CW**2*sb**2-4D0*Zss22**4
     &      /CW**2*sb**2-9D0*Zss21**4/SW**2*sb**2-4D0*Zss21
     &      **2/CW**2*Zss22**2*sb**2+36D0*Mq**2/MW**2/SW
     &      **2*Zss21**2*Zss22**2)
      AAABR(1477) = EE**2/cb**2*(Zb21**2/CW**2*Zdd21**2+2D0
     &      *Zb21**2/CW**2*Zdd22**2+2D0*Zb22**2/CW**2*Zdd21
     &      **2+4D0*Zb22**2/CW**2*Zdd22**2+9D0*Zb21**2/SW
     &      **2*Zdd21**2-Zb21**2/CW**2*Zdd21**2*sb**2-2D0
     &      *Zb21**2/CW**2*Zdd22**2*sb**2-2D0*Zb22**2/CW
     &      **2*Zdd21**2*sb**2-4D0*Zb22**2/CW**2*Zdd22**2
     &      *sb**2-9D0*Zb21**2/SW**2*Zdd21**2*sb**2+36D0
     &      *MbMM/MW**2*Mq/SW**2*Zb21*Zb22*Zdd21*Zdd22)
      AAABR(1478) = GG**2/cb**2*(Zb21**4+Zb22**4+2D0*Zb21
     &      **2*Zb22**2-Zb21**4*sb**2-Zb22**4*sb**2-2D0*Zb21
     &      **2*Zb22**2*sb**2)
      AAABR(1479) = EE**2/cb**2*(Zb21**4/CW**2+4D0*Zb22**4
     &      /CW**2+9D0*Zb21**4/SW**2+4D0*Zb21**2/CW**2*Zb22
     &      **2-Zb21**4/CW**2*sb**2-4D0*Zb22**4/CW**2*sb
     &      **2-9D0*Zb21**4/SW**2*sb**2-4D0*Zb21**2/CW**2
     &      *Zb22**2*sb**2+36D0*MbMM**2/MW**2/SW**2*Zb21
     &      **2*Zb22**2)
      AAABR(1480) = EE**2/cb**2*(Zdd21**2/CW**2*Zss21**2
     &      +2D0*Zdd21**2/CW**2*Zss22**2+2D0*Zdd22**2/CW
     &      **2*Zss21**2+4D0*Zdd22**2/CW**2*Zss22**2+9D0
     &      *Zdd21**2/SW**2*Zss21**2-Zdd21**2/CW**2*Zss21
     &      **2*sb**2-2D0*Zdd21**2/CW**2*Zss22**2*sb**2-2D0
     &      *Zdd22**2/CW**2*Zss21**2*sb**2-4D0*Zdd22**2/CW
     &      **2*Zss22**2*sb**2-9D0*Zdd21**2/SW**2*Zss21**2
     &      *sb**2+36D0*Mq**2/MW**2/SW**2*Zdd21*Zdd22*Zss21
     &      *Zss22)
      AAABR(1481) = EE**2/cb**2*(Zb21**2/CW**2*Zss21**2+2D0
     &      *Zb21**2/CW**2*Zss22**2+2D0*Zb22**2/CW**2*Zss21
     &      **2+4D0*Zb22**2/CW**2*Zss22**2+9D0*Zb21**2/SW
     &      **2*Zss21**2-Zb21**2/CW**2*Zss21**2*sb**2-2D0
     &      *Zb21**2/CW**2*Zss22**2*sb**2-2D0*Zb22**2/CW
     &      **2*Zss21**2*sb**2-4D0*Zb22**2/CW**2*Zss22**2
     &      *sb**2-9D0*Zb21**2/SW**2*Zss21**2*sb**2+36D0
     &      *MbMM/MW**2*Mq/SW**2*Zb21*Zb22*Zss21*Zss22)
      AAABR(1482) = GG**2/cb**2*(Zdd21**4+Zdd22**4+2D0*Zdd21
     &      **2*Zdd22**2-Zdd21**4*sb**2-Zdd22**4*sb**2-2D0
     &      *Zdd21**2*Zdd22**2*sb**2)
      AAABR(1483) = EE**2/cb**2*(Zdd21**4/CW**2+4D0*Zdd22
     &      **4/CW**2+9D0*Zdd21**4/SW**2+4D0*Zdd21**2/CW
     &      **2*Zdd22**2-Zdd21**4/CW**2*sb**2-4D0*Zdd22**4
     &      /CW**2*sb**2-9D0*Zdd21**4/SW**2*sb**2-4D0*Zdd21
     &      **2/CW**2*Zdd22**2*sb**2+36D0*Mq**2/MW**2/SW
     &      **2*Zdd21**2*Zdd22**2)
      AAABR(1484) = EE**2*(Zcc11**2/CW**2*Zdd21**2+2D0*Zcc11
     &      **2/CW**2*Zdd22**2-4D0*Zcc12**2/CW**2*Zdd21**2
     &      -8D0*Zcc12**2/CW**2*Zdd22**2-9D0*Zcc11**2/SW
     &      **2*Zdd21**2)
      AAABR(1485) = EE**2/cb**2*(Zdd21**2/CW**2*Zuu11**2
     &      -4D0*Zdd21**2/CW**2*Zuu12**2+2D0*Zdd22**2/CW
     &      **2*Zuu11**2-8D0*Zdd22**2/CW**2*Zuu12**2-9D0
     &      *Zdd21**2/SW**2*Zuu11**2-Zdd21**2/CW**2*Zuu11
     &      **2*sb**2+4D0*Zdd21**2/CW**2*Zuu12**2*sb**2-2D0
     &      *Zdd22**2/CW**2*Zuu11**2*sb**2+8D0*Zdd22**2/CW
     &      **2*Zuu12**2*sb**2+9D0*Zdd21**2/SW**2*Zuu11**2
     &      *sb**2)
      AAABR(1486) = EE**2/SW**2/cb**2*(Zdd21**2*Zuu11**2
     &      -Zdd21**2*Zuu11**2*sb**2-Mq**2/MW**2*Zdd21**2
     &      *Zuu12**2+Mq**2/MW**2*Zdd22**2*Zuu11**2+Mq**2
     &      /MW**2*Zdd21**2*Zuu12**2/sb**2)
      AAABR(1487) = EE**2/cb**2*(Zb21**2/CW**2*Zt11**2-4D0
     &      *Zb21**2/CW**2*Zt12**2+2D0*Zb22**2/CW**2*Zt11
     &      **2-8D0*Zb22**2/CW**2*Zt12**2-9D0*Zb21**2/SW
     &      **2*Zt11**2-Zb21**2/CW**2*Zt11**2*sb**2+4D0*Zb21
     &      **2/CW**2*Zt12**2*sb**2-2D0*Zb22**2/CW**2*Zt11
     &      **2*sb**2+8D0*Zb22**2/CW**2*Zt12**2*sb**2+9D0
     &      *Zb21**2/SW**2*Zt11**2*sb**2)
      AAABR(1488) = EE**2/SW**2/cb**2*(Zb21**2*Zt11**2-Zb21
     &      **2*Zt11**2*sb**2+MbMM**2/MW**2*Zb22**2*Zt11
     &      **2-MtMM**2/MW**2*Zb21**2*Zt12**2+MtMM**2/MW
     &      **2*Zb21**2*Zt12**2/sb**2)
      AAABR(1489) = EE**2*(Zb21**2/CW**2*Zuu11**2-4D0*Zb21
     &      **2/CW**2*Zuu12**2+2D0*Zb22**2/CW**2*Zuu11**2
     &      -8D0*Zb22**2/CW**2*Zuu12**2-9D0*Zb21**2/SW**2
     &      *Zuu11**2)
      AAABR(1490) = EE**2*(Zb21**2/CW**2*Zcc11**2-4D0*Zb21
     &      **2/CW**2*Zcc12**2+2D0*Zb22**2/CW**2*Zcc11**2
     &      -8D0*Zb22**2/CW**2*Zcc12**2-9D0*Zb21**2/SW**2
     &      *Zcc11**2)
      AAABR(1491) = EE**2*(Zdd21**2/CW**2*Zt11**2-4D0*Zdd21
     &      **2/CW**2*Zt12**2+2D0*Zdd22**2/CW**2*Zt11**2
     &      -8D0*Zdd22**2/CW**2*Zt12**2-9D0*Zdd21**2/SW**2
     &      *Zt11**2)
      AAABR(1492) = EE**2*(Zss21**2/CW**2*Zuu11**2-4D0*Zss21
     &      **2/CW**2*Zuu12**2+2D0*Zss22**2/CW**2*Zuu11**2
     &      -8D0*Zss22**2/CW**2*Zuu12**2-9D0*Zss21**2/SW
     &      **2*Zuu11**2)
      AAABR(1493) = EE**2/SW**2/cb**2*(Zcc11*Zdd21*Zss21
     &      *Zuu11-Zcc11*Zdd21*Zss21*Zuu11*sb**2+Mq**2/MW
     &      **2*Zcc11*Zdd22*Zss22*Zuu11-Mc/MW**2*Mq*Zcc12
     &      *Zdd21*Zss21*Zuu12+Mc/MW**2*Mq*Zcc12*Zdd21*Zss21
     &      *Zuu12/sb**2)
      AAABR(1494) = EE**2/SW**2/cb**2*(Zb21*Zdd21*Zt11*Zuu11
     &      -Zb21*Zdd21*Zt11*Zuu11*sb**2+MbMM/MW**2*Mq*Zb22
     &      *Zdd22*Zt11*Zuu11-Mq/MW**2*MtMM*Zb21*Zdd21*Zt12
     &      *Zuu12+Mq/MW**2*MtMM*Zb21*Zdd21*Zt12*Zuu12/sb
     &      **2)
      AAABR(1495) = EE**2/SW**2/cb**2*(Zb21*Zcc11*Zss21*Zt11
     &      -Zb21*Zcc11*Zss21*Zt11*sb**2+MbMM/MW**2*Mq*Zb22
     &      *Zcc11*Zss22*Zt11-Mc/MW**2*MtMM*Zb21*Zcc12*Zss21
     &      *Zt12+Mc/MW**2*MtMM*Zb21*Zcc12*Zss21*Zt12/sb**2)
      AAABR(1496) = EE**2*(Zss21**2/CW**2*Zt11**2-4D0*Zss21
     &      **2/CW**2*Zt12**2+2D0*Zss22**2/CW**2*Zt11**2
     &      -8D0*Zss22**2/CW**2*Zt12**2-9D0*Zss21**2/SW**2
     &      *Zt11**2)
      AAABR(1497) = EE**2/cb**2*(Zcc11**2/CW**2*Zss21**2
     &      +2D0*Zcc11**2/CW**2*Zss22**2-4D0*Zcc12**2/CW
     &      **2*Zss21**2-8D0*Zcc12**2/CW**2*Zss22**2-9D0
     &      *Zcc11**2/SW**2*Zss21**2-Zcc11**2/CW**2*Zss21
     &      **2*sb**2-2D0*Zcc11**2/CW**2*Zss22**2*sb**2+4D0
     &      *Zcc12**2/CW**2*Zss21**2*sb**2+8D0*Zcc12**2/CW
     &      **2*Zss22**2*sb**2+9D0*Zcc11**2/SW**2*Zss21**2
     &      *sb**2)
      AAABR(1498) = EE**2/SW**2/cb**2*(Zcc11**2*Zss21**2
     &      -Zcc11**2*Zss21**2*sb**2-Mc**2/MW**2*Zcc12**2
     &      *Zss21**2+Mq**2/MW**2*Zcc11**2*Zss22**2+Mc**2
     &      /MW**2*Zcc12**2*Zss21**2/sb**2)
      AAABR(1499) = EE**2/cb**2*(Zb21**2/CW**2*Zt11*Zt21
     &      -4D0*Zb21**2/CW**2*Zt12*Zt22+2D0*Zb22**2/CW**2
     &      *Zt11*Zt21-8D0*Zb22**2/CW**2*Zt12*Zt22-9D0*Zb21
     &      **2/SW**2*Zt11*Zt21-Zb21**2/CW**2*Zt11*Zt21*sb
     &      **2+4D0*Zb21**2/CW**2*Zt12*Zt22*sb**2-2D0*Zb22
     &      **2/CW**2*Zt11*Zt21*sb**2+8D0*Zb22**2/CW**2*Zt12
     &      *Zt22*sb**2+9D0*Zb21**2/SW**2*Zt11*Zt21*sb**2)
      AAABR(1500) = EE**2/SW**2/cb**2*(Zb21**2*Zt11*Zt21
     &      -Zb21**2*Zt11*Zt21*sb**2+MbMM**2/MW**2*Zb22**2
     &      *Zt11*Zt21-MtMM**2/MW**2*Zb21**2*Zt12*Zt22+MtMM
     &      **2/MW**2*Zb21**2*Zt12*Zt22/sb**2)
      AAABR(1501) = EE**2/SW**2/cb**2*(Zb21*Zcc21*Zss21*Zt11
     &      -Zb21*Zcc21*Zss21*Zt11*sb**2+MbMM/MW**2*Mq*Zb22
     &      *Zcc21*Zss22*Zt11-Mc/MW**2*MtMM*Zb21*Zcc22*Zss21
     &      *Zt12+Mc/MW**2*MtMM*Zb21*Zcc22*Zss21*Zt12/sb**2)
      AAABR(1502) = EE**2/SW**2/cb**2*(Zb21*Zdd21*Zt11*Zuu21
     &      -Zb21*Zdd21*Zt11*Zuu21*sb**2+MbMM/MW**2*Mq*Zb22
     &      *Zdd22*Zt11*Zuu21-Mq/MW**2*MtMM*Zb21*Zdd21*Zt12
     &      *Zuu22+Mq/MW**2*MtMM*Zb21*Zdd21*Zt12*Zuu22/sb
     &      **2)
      AAABR(1503) = EE**2/cb**2*(Zdd21**2/CW**2*Zuu11*Zuu21
     &      -4D0*Zdd21**2/CW**2*Zuu12*Zuu22+2D0*Zdd22**2
     &      /CW**2*Zuu11*Zuu21-8D0*Zdd22**2/CW**2*Zuu12*Zuu22
     &      -9D0*Zdd21**2/SW**2*Zuu11*Zuu21-Zdd21**2/CW**2
     &      *Zuu11*Zuu21*sb**2+4D0*Zdd21**2/CW**2*Zuu12*Zuu22
     &      *sb**2-2D0*Zdd22**2/CW**2*Zuu11*Zuu21*sb**2+8D0
     &      *Zdd22**2/CW**2*Zuu12*Zuu22*sb**2+9D0*Zdd21**2
     &      /SW**2*Zuu11*Zuu21*sb**2)
      AAABR(1504) = EE**2/SW**2/cb**2*(Zdd21**2*Zuu11*Zuu21
     &      -Zdd21**2*Zuu11*Zuu21*sb**2-Mq**2/MW**2*Zdd21
     &      **2*Zuu12*Zuu22+Mq**2/MW**2*Zdd22**2*Zuu11*Zuu21
     &      +Mq**2/MW**2*Zdd21**2*Zuu12*Zuu22/sb**2)
      AAABR(1505) = EE**2*(Zss21**2/CW**2*Zt11*Zt21-4D0*Zss21
     &      **2/CW**2*Zt12*Zt22+2D0*Zss22**2/CW**2*Zt11*Zt21
     &      -8D0*Zss22**2/CW**2*Zt12*Zt22-9D0*Zss21**2/SW
     &      **2*Zt11*Zt21)
      AAABR(1506) = EE**2*(Zss21**2/CW**2*Zuu11*Zuu21-4D0
     &      *Zss21**2/CW**2*Zuu12*Zuu22+2D0*Zss22**2/CW**2
     &      *Zuu11*Zuu21-8D0*Zss22**2/CW**2*Zuu12*Zuu22-9D0
     &      *Zss21**2/SW**2*Zuu11*Zuu21)
      AAABR(1507) = EE**2/cb**2*(Zcc11/CW**2*Zcc21*Zss21
     &      **2+2D0*Zcc11/CW**2*Zcc21*Zss22**2-4D0*Zcc12
     &      /CW**2*Zcc22*Zss21**2-8D0*Zcc12/CW**2*Zcc22*Zss22
     &      **2-9D0*Zcc11/SW**2*Zcc21*Zss21**2-Zcc11/CW**2
     &      *Zcc21*Zss21**2*sb**2-2D0*Zcc11/CW**2*Zcc21*Zss22
     &      **2*sb**2+4D0*Zcc12/CW**2*Zcc22*Zss21**2*sb**2
     &      +8D0*Zcc12/CW**2*Zcc22*Zss22**2*sb**2+9D0*Zcc11
     &      /SW**2*Zcc21*Zss21**2*sb**2)
      AAABR(1508) = EE**2/SW**2/cb**2*(Zcc11*Zcc21*Zss21
     &      **2-Zcc11*Zcc21*Zss21**2*sb**2-Mc**2/MW**2*Zcc12
     &      *Zcc22*Zss21**2+Mq**2/MW**2*Zcc11*Zcc21*Zss22
     &      **2+Mc**2/MW**2*Zcc12*Zcc22*Zss21**2/sb**2)
      AAABR(1509) = EE**2*(Zdd21**2/CW**2*Zt11*Zt21-4D0*Zdd21
     &      **2/CW**2*Zt12*Zt22+2D0*Zdd22**2/CW**2*Zt11*Zt21
     &      -8D0*Zdd22**2/CW**2*Zt12*Zt22-9D0*Zdd21**2/SW
     &      **2*Zt11*Zt21)
      AAABR(1510) = EE**2/SW**2/cb**2*(Zcc21*Zdd21*Zss21
     &      *Zuu11-Zcc21*Zdd21*Zss21*Zuu11*sb**2+Mq**2/MW
     &      **2*Zcc21*Zdd22*Zss22*Zuu11-Mc/MW**2*Mq*Zcc22
     &      *Zdd21*Zss21*Zuu12+Mc/MW**2*Mq*Zcc22*Zdd21*Zss21
     &      *Zuu12/sb**2)
      AAABR(1511) = EE**2*(Zb21**2/CW**2*Zuu11*Zuu21-4D0
     &      *Zb21**2/CW**2*Zuu12*Zuu22+2D0*Zb22**2/CW**2
     &      *Zuu11*Zuu21-8D0*Zb22**2/CW**2*Zuu12*Zuu22-9D0
     &      *Zb21**2/SW**2*Zuu11*Zuu21)
      AAABR(1512) = EE**2*(Zb21**2/CW**2*Zcc11*Zcc21-4D0
     &      *Zb21**2/CW**2*Zcc12*Zcc22+2D0*Zb22**2/CW**2
     &      *Zcc11*Zcc21-8D0*Zb22**2/CW**2*Zcc12*Zcc22-9D0
     &      *Zb21**2/SW**2*Zcc11*Zcc21)
      AAABR(1513) = EE**2/SW**2/cb**2*(Zcc11*Zdd21*Zss21
     &      *Zuu21-Zcc11*Zdd21*Zss21*Zuu21*sb**2+Mq**2/MW
     &      **2*Zcc11*Zdd22*Zss22*Zuu21-Mc/MW**2*Mq*Zcc12
     &      *Zdd21*Zss21*Zuu22+Mc/MW**2*Mq*Zcc12*Zdd21*Zss21
     &      *Zuu22/sb**2)
      AAABR(1514) = EE**2/SW**2/cb**2*(Zb21*Zdd21*Zt21*Zuu11
     &      -Zb21*Zdd21*Zt21*Zuu11*sb**2+MbMM/MW**2*Mq*Zb22
     &      *Zdd22*Zt21*Zuu11-Mq/MW**2*MtMM*Zb21*Zdd21*Zt22
     &      *Zuu12+Mq/MW**2*MtMM*Zb21*Zdd21*Zt22*Zuu12/sb
     &      **2)
      AAABR(1515) = EE**2/SW**2/cb**2*(Zb21*Zcc11*Zss21*Zt21
     &      -Zb21*Zcc11*Zss21*Zt21*sb**2+MbMM/MW**2*Mq*Zb22
     &      *Zcc11*Zss22*Zt21-Mc/MW**2*MtMM*Zb21*Zcc12*Zss21
     &      *Zt22+Mc/MW**2*MtMM*Zb21*Zcc12*Zss21*Zt22/sb**2)
      AAABR(1516) = EE**2*(Zcc11/CW**2*Zcc21*Zdd21**2+2D0
     &      *Zcc11/CW**2*Zcc21*Zdd22**2-4D0*Zcc12/CW**2*Zcc22
     &      *Zdd21**2-8D0*Zcc12/CW**2*Zcc22*Zdd22**2-9D0
     &      *Zcc11/SW**2*Zcc21*Zdd21**2)
      AAABR(1517) = EE**2/cb**2*(Zb21**2/CW**2*Zt21**2-4D0
     &      *Zb21**2/CW**2*Zt22**2+2D0*Zb22**2/CW**2*Zt21
     &      **2-8D0*Zb22**2/CW**2*Zt22**2-9D0*Zb21**2/SW
     &      **2*Zt21**2-Zb21**2/CW**2*Zt21**2*sb**2+4D0*Zb21
     &      **2/CW**2*Zt22**2*sb**2-2D0*Zb22**2/CW**2*Zt21
     &      **2*sb**2+8D0*Zb22**2/CW**2*Zt22**2*sb**2+9D0
     &      *Zb21**2/SW**2*Zt21**2*sb**2)
      AAABR(1518) = EE**2/SW**2/cb**2*(Zb21**2*Zt21**2-Zb21
     &      **2*Zt21**2*sb**2+MbMM**2/MW**2*Zb22**2*Zt21
     &      **2-MtMM**2/MW**2*Zb21**2*Zt22**2+MtMM**2/MW
     &      **2*Zb21**2*Zt22**2/sb**2)
      AAABR(1519) = EE**2*(Zb21**2/CW**2*Zuu21**2-4D0*Zb21
     &      **2/CW**2*Zuu22**2+2D0*Zb22**2/CW**2*Zuu21**2
     &      -8D0*Zb22**2/CW**2*Zuu22**2-9D0*Zb21**2/SW**2
     &      *Zuu21**2)
      AAABR(1520) = EE**2/SW**2/cb**2*(Zb21*Zcc21*Zss21*Zt21
     &      -Zb21*Zcc21*Zss21*Zt21*sb**2+MbMM/MW**2*Mq*Zb22
     &      *Zcc21*Zss22*Zt21-Mc/MW**2*MtMM*Zb21*Zcc22*Zss21
     &      *Zt22+Mc/MW**2*MtMM*Zb21*Zcc22*Zss21*Zt22/sb**2)
      AAABR(1521) = EE**2/SW**2/cb**2*(Zcc21*Zdd21*Zss21
     &      *Zuu21-Zcc21*Zdd21*Zss21*Zuu21*sb**2+Mq**2/MW
     &      **2*Zcc21*Zdd22*Zss22*Zuu21-Mc/MW**2*Mq*Zcc22
     &      *Zdd21*Zss21*Zuu22+Mc/MW**2*Mq*Zcc22*Zdd21*Zss21
     &      *Zuu22/sb**2)
      AAABR(1522) = EE**2/cb**2*(Zcc21**2/CW**2*Zss21**2
     &      +2D0*Zcc21**2/CW**2*Zss22**2-4D0*Zcc22**2/CW
     &      **2*Zss21**2-8D0*Zcc22**2/CW**2*Zss22**2-9D0
     &      *Zcc21**2/SW**2*Zss21**2-Zcc21**2/CW**2*Zss21
     &      **2*sb**2-2D0*Zcc21**2/CW**2*Zss22**2*sb**2+4D0
     &      *Zcc22**2/CW**2*Zss21**2*sb**2+8D0*Zcc22**2/CW
     &      **2*Zss22**2*sb**2+9D0*Zcc21**2/SW**2*Zss21**2
     &      *sb**2)
      AAABR(1523) = EE**2/SW**2/cb**2*(Zcc21**2*Zss21**2
     &      -Zcc21**2*Zss21**2*sb**2-Mc**2/MW**2*Zcc22**2
     &      *Zss21**2+Mq**2/MW**2*Zcc21**2*Zss22**2+Mc**2
     &      /MW**2*Zcc22**2*Zss21**2/sb**2)
      AAABR(1524) = EE**2*(Zss21**2/CW**2*Zt21**2-4D0*Zss21
     &      **2/CW**2*Zt22**2+2D0*Zss22**2/CW**2*Zt21**2
     &      -8D0*Zss22**2/CW**2*Zt22**2-9D0*Zss21**2/SW**2
     &      *Zt21**2)
      AAABR(1525) = EE**2/SW**2/cb**2*(Zb21*Zdd21*Zt21*Zuu21
     &      -Zb21*Zdd21*Zt21*Zuu21*sb**2+MbMM/MW**2*Mq*Zb22
     &      *Zdd22*Zt21*Zuu21-Mq/MW**2*MtMM*Zb21*Zdd21*Zt22
     &      *Zuu22+Mq/MW**2*MtMM*Zb21*Zdd21*Zt22*Zuu22/sb
     &      **2)
      AAABR(1526) = EE**2*(Zb21**2/CW**2*Zcc21**2-4D0*Zb21
     &      **2/CW**2*Zcc22**2+2D0*Zb22**2/CW**2*Zcc21**2
     &      -8D0*Zb22**2/CW**2*Zcc22**2-9D0*Zb21**2/SW**2
     &      *Zcc21**2)
      AAABR(1527) = EE**2*(Zcc21**2/CW**2*Zdd21**2+2D0*Zcc21
     &      **2/CW**2*Zdd22**2-4D0*Zcc22**2/CW**2*Zdd21**2
     &      -8D0*Zcc22**2/CW**2*Zdd22**2-9D0*Zcc21**2/SW
     &      **2*Zdd21**2)
      AAABR(1528) = EE**2/cb**2*(Zdd21**2/CW**2*Zuu21**2
     &      -4D0*Zdd21**2/CW**2*Zuu22**2+2D0*Zdd22**2/CW
     &      **2*Zuu21**2-8D0*Zdd22**2/CW**2*Zuu22**2-9D0
     &      *Zdd21**2/SW**2*Zuu21**2-Zdd21**2/CW**2*Zuu21
     &      **2*sb**2+4D0*Zdd21**2/CW**2*Zuu22**2*sb**2-2D0
     &      *Zdd22**2/CW**2*Zuu21**2*sb**2+8D0*Zdd22**2/CW
     &      **2*Zuu22**2*sb**2+9D0*Zdd21**2/SW**2*Zuu21**2
     &      *sb**2)
      AAABR(1529) = EE**2/SW**2/cb**2*(Zdd21**2*Zuu21**2
     &      -Zdd21**2*Zuu21**2*sb**2-Mq**2/MW**2*Zdd21**2
     &      *Zuu22**2+Mq**2/MW**2*Zdd22**2*Zuu21**2+Mq**2
     &      /MW**2*Zdd21**2*Zuu22**2/sb**2)
      AAABR(1530) = EE**2*(Zss21**2/CW**2*Zuu21**2-4D0*Zss21
     &      **2/CW**2*Zuu22**2+2D0*Zss22**2/CW**2*Zuu21**2
     &      -8D0*Zss22**2/CW**2*Zuu22**2-9D0*Zss21**2/SW
     &      **2*Zuu21**2)
      AAABR(1531) = EE**2*(Zdd21**2/CW**2*Zt21**2-4D0*Zdd21
     &      **2/CW**2*Zt22**2+2D0*Zdd22**2/CW**2*Zt21**2
     &      -8D0*Zdd22**2/CW**2*Zt22**2-9D0*Zdd21**2/SW**2
     &      *Zt21**2)
      AAABR(1532) = EE**2*(Zdd21**2/CW**2+2D0*Zdd22**2/CW
     &      **2-3D0*Zdd21**2/SW**2-2D0*Zdd21**2/CW**2*sb
     &      **2-4D0*Zdd22**2/CW**2*sb**2-6D0*Mq**2/MW**2
     &      /SW**2+6D0*Zdd21**2/SW**2*sb**2)
      AAABR(1533) = EE**2*(Zb21**2/CW**2+2D0*Zb22**2/CW**2
     &      -3D0*Zb21**2/SW**2-2D0*Zb21**2/CW**2*sb**2-4D0
     &      *Zb22**2/CW**2*sb**2+6D0*Zb21**2/SW**2*sb**2
     &      -6D0*MbMM**2/MW**2/SW**2*Zb22**2-6D0*MtMM**2
     &      /MW**2/SW**2*Zb21**2)
      AAABR(1534) = EE**2*(Zss21**2/CW**2+2D0*Zss22**2/CW
     &      **2-3D0*Zss21**2/SW**2-2D0*Zss21**2/CW**2*sb
     &      **2-4D0*Zss22**2/CW**2*sb**2+6D0*Zss21**2/SW
     &      **2*sb**2-6D0*Mc**2/MW**2/SW**2*Zss21**2-6D0
     &      *Mq**2/MW**2/SW**2*Zss22**2)
      AAABR(1535) = EE**2/SW**2*Sqrt2*Zdd21*Zuu11*(1D0-2D0
     &      *sb**2)
      AAABR(1536) = EE**2/SW**2*Sqrt2*Zcc11*Zss21*(1D0-2D0
     &      *sb**2+Mc**2/MW**2-Mq**2/MW**2)
      AAABR(1537) = EE**2/SW**2*Sqrt2*Zb21*Zt11*(1D0-2D0
     &      *sb**2-MbMM**2/MW**2+MtMM**2/MW**2)
      AAABR(1538) = EE**2/SW**2*Sqrt2*Zdd21*Zuu21*(1D0-2D0
     &      *sb**2)
      AAABR(1539) = EE**2/SW**2*Sqrt2*Zcc21*Zss21*(1D0-2D0
     &      *sb**2+Mc**2/MW**2-Mq**2/MW**2)
      AAABR(1540) = EE**2/SW**2*Sqrt2*Zb21*Zt21*(1D0-2D0
     &      *sb**2-MbMM**2/MW**2+MtMM**2/MW**2)
      AAABR(1541) = 1D0/cb*(2D0*EE**2*Pa11/SW**2*Sqrt2*Zb21
     &      *Zt11*sb-2D0*EE**2*Pa11/SW**2*Sqrt2*Zb21*Zt11
     &      *sb**3+2D0*EE*Lambda/MW*MbMM*Pa12/SW*Zb22*Zt11
     &      *sb+2D0*EE*Lambda/MW*MtMM*Pa12/SW*Zb21*Zt12/sb
     &      -2D0*EE*Lambda/MW*MtMM*Pa12/SW*Zb21*Zt12*sb-EE
     &      **2/MW**2*MbMM**2*Pa11/SW**2*Sqrt2*Zb21*Zt11
     &      *sb-EE**2/MW**2*MtMM**2*Pa11/SW**2*Sqrt2*Zb21
     &      *Zt11/sb+EE**2/MW**2*MtMM**2*Pa11/SW**2*Sqrt2
     &      *Zb21*Zt11*sb-EE**2/MW**2*MbMM*MtMM*Pa11/SW**2
     &      *Sqrt2*Zb22*Zt12/sb)
      AAABR(1542) = 1D0/cb*(2D0*EE**2*Pa21/SW**2*Sqrt2*Zb21
     &      *Zt11*sb-2D0*EE**2*Pa21/SW**2*Sqrt2*Zb21*Zt11
     &      *sb**3+2D0*EE*Lambda/MW*MbMM*Pa22/SW*Zb22*Zt11
     &      *sb+2D0*EE*Lambda/MW*MtMM*Pa22/SW*Zb21*Zt12/sb
     &      -2D0*EE*Lambda/MW*MtMM*Pa22/SW*Zb21*Zt12*sb-EE
     &      **2/MW**2*MbMM**2*Pa21/SW**2*Sqrt2*Zb21*Zt11
     &      *sb-EE**2/MW**2*MtMM**2*Pa21/SW**2*Sqrt2*Zb21
     &      *Zt11/sb+EE**2/MW**2*MtMM**2*Pa21/SW**2*Sqrt2
     &      *Zb21*Zt11*sb-EE**2/MW**2*MbMM*MtMM*Pa21/SW**2
     &      *Sqrt2*Zb22*Zt12/sb)
      AAABR(1543) = 1D0/cb*(2D0*EE**2*Pa21/SW**2*Sqrt2*Zdd21
     &      *Zuu11*sb-2D0*EE**2*Pa21/SW**2*Sqrt2*Zdd21*Zuu11
     &      *sb**3+2D0*EE*Lambda/MW*Mq*Pa22/SW*Zdd21*Zuu12
     &      /sb-2D0*EE*Lambda/MW*Mq*Pa22/SW*Zdd21*Zuu12*sb
     &      +2D0*EE*Lambda/MW*Mq*Pa22/SW*Zdd22*Zuu11*sb-EE
     &      **2/MW**2*Mq**2*Pa21/SW**2*Sqrt2*Zdd21*Zuu11
     &      /sb-EE**2/MW**2*Mq**2*Pa21/SW**2*Sqrt2*Zdd22
     &      *Zuu12/sb)
      AAABR(1544) = 1D0/cb*(2D0*EE**2*Pa11/SW**2*Sqrt2*Zcc11
     &      *Zss21*sb-2D0*EE**2*Pa11/SW**2*Sqrt2*Zcc11*Zss21
     &      *sb**3+2D0*EE*Lambda/MW*Mc*Pa12/SW*Zcc12*Zss21
     &      /sb-2D0*EE*Lambda/MW*Mc*Pa12/SW*Zcc12*Zss21*sb
     &      +2D0*EE*Lambda/MW*Mq*Pa12/SW*Zcc11*Zss22*sb-EE
     &      **2/MW**2*Mc**2*Pa11/SW**2*Sqrt2*Zcc11*Zss21
     &      /sb+EE**2/MW**2*Mc**2*Pa11/SW**2*Sqrt2*Zcc11
     &      *Zss21*sb-EE**2/MW**2*Mq**2*Pa11/SW**2*Sqrt2
     &      *Zcc11*Zss21*sb-EE**2/MW**2*Mc*Mq*Pa11/SW**2
     &      *Sqrt2*Zcc12*Zss22/sb)
      AAABR(1545) = 1D0/cb*(2D0*EE**2*Pa11/SW**2*Sqrt2*Zdd21
     &      *Zuu11*sb-2D0*EE**2*Pa11/SW**2*Sqrt2*Zdd21*Zuu11
     &      *sb**3+2D0*EE*Lambda/MW*Mq*Pa12/SW*Zdd21*Zuu12
     &      /sb-2D0*EE*Lambda/MW*Mq*Pa12/SW*Zdd21*Zuu12*sb
     &      +2D0*EE*Lambda/MW*Mq*Pa12/SW*Zdd22*Zuu11*sb-EE
     &      **2/MW**2*Mq**2*Pa11/SW**2*Sqrt2*Zdd21*Zuu11
     &      /sb-EE**2/MW**2*Mq**2*Pa11/SW**2*Sqrt2*Zdd22
     &      *Zuu12/sb)
      AAABR(1546) = 1D0/cb*(2D0*EE**2*Pa21/SW**2*Sqrt2*Zcc11
     &      *Zss21*sb-2D0*EE**2*Pa21/SW**2*Sqrt2*Zcc11*Zss21
     &      *sb**3+2D0*EE*Lambda/MW*Mc*Pa22/SW*Zcc12*Zss21
     &      /sb-2D0*EE*Lambda/MW*Mc*Pa22/SW*Zcc12*Zss21*sb
     &      +2D0*EE*Lambda/MW*Mq*Pa22/SW*Zcc11*Zss22*sb-EE
     &      **2/MW**2*Mc**2*Pa21/SW**2*Sqrt2*Zcc11*Zss21
     &      /sb+EE**2/MW**2*Mc**2*Pa21/SW**2*Sqrt2*Zcc11
     &      *Zss21*sb-EE**2/MW**2*Mq**2*Pa21/SW**2*Sqrt2
     &      *Zcc11*Zss21*sb-EE**2/MW**2*Mc*Mq*Pa21/SW**2
     &      *Sqrt2*Zcc12*Zss22/sb)
      AAABR(1547) = 1D0/cb*(2D0*EE**2*Pa11/SW**2*Sqrt2*Zcc21
     &      *Zss21*sb-2D0*EE**2*Pa11/SW**2*Sqrt2*Zcc21*Zss21
     &      *sb**3+2D0*EE*Lambda/MW*Mc*Pa12/SW*Zcc22*Zss21
     &      /sb-2D0*EE*Lambda/MW*Mc*Pa12/SW*Zcc22*Zss21*sb
     &      +2D0*EE*Lambda/MW*Mq*Pa12/SW*Zcc21*Zss22*sb-EE
     &      **2/MW**2*Mc**2*Pa11/SW**2*Sqrt2*Zcc21*Zss21
     &      /sb+EE**2/MW**2*Mc**2*Pa11/SW**2*Sqrt2*Zcc21
     &      *Zss21*sb-EE**2/MW**2*Mq**2*Pa11/SW**2*Sqrt2
     &      *Zcc21*Zss21*sb-EE**2/MW**2*Mc*Mq*Pa11/SW**2
     &      *Sqrt2*Zcc22*Zss22/sb)
      AAABR(1548) = 1D0/cb*(2D0*EE**2*Pa11/SW**2*Sqrt2*Zb21
     &      *Zt21*sb-2D0*EE**2*Pa11/SW**2*Sqrt2*Zb21*Zt21
     &      *sb**3+2D0*EE*Lambda/MW*MbMM*Pa12/SW*Zb22*Zt21
     &      *sb+2D0*EE*Lambda/MW*MtMM*Pa12/SW*Zb21*Zt22/sb
     &      -2D0*EE*Lambda/MW*MtMM*Pa12/SW*Zb21*Zt22*sb-EE
     &      **2/MW**2*MbMM**2*Pa11/SW**2*Sqrt2*Zb21*Zt21
     &      *sb-EE**2/MW**2*MtMM**2*Pa11/SW**2*Sqrt2*Zb21
     &      *Zt21/sb+EE**2/MW**2*MtMM**2*Pa11/SW**2*Sqrt2
     &      *Zb21*Zt21*sb-EE**2/MW**2*MbMM*MtMM*Pa11/SW**2
     &      *Sqrt2*Zb22*Zt22/sb)
      AAABR(1549) = 1D0/cb*(2D0*EE**2*Pa21/SW**2*Sqrt2*Zdd21
     &      *Zuu21*sb-2D0*EE**2*Pa21/SW**2*Sqrt2*Zdd21*Zuu21
     &      *sb**3+2D0*EE*Lambda/MW*Mq*Pa22/SW*Zdd21*Zuu22
     &      /sb-2D0*EE*Lambda/MW*Mq*Pa22/SW*Zdd21*Zuu22*sb
     &      +2D0*EE*Lambda/MW*Mq*Pa22/SW*Zdd22*Zuu21*sb-EE
     &      **2/MW**2*Mq**2*Pa21/SW**2*Sqrt2*Zdd21*Zuu21
     &      /sb-EE**2/MW**2*Mq**2*Pa21/SW**2*Sqrt2*Zdd22
     &      *Zuu22/sb)
      AAABR(1550) = 1D0/cb*(2D0*EE**2*Pa11/SW**2*Sqrt2*Zdd21
     &      *Zuu21*sb-2D0*EE**2*Pa11/SW**2*Sqrt2*Zdd21*Zuu21
     &      *sb**3+2D0*EE*Lambda/MW*Mq*Pa12/SW*Zdd21*Zuu22
     &      /sb-2D0*EE*Lambda/MW*Mq*Pa12/SW*Zdd21*Zuu22*sb
     &      +2D0*EE*Lambda/MW*Mq*Pa12/SW*Zdd22*Zuu21*sb-EE
     &      **2/MW**2*Mq**2*Pa11/SW**2*Sqrt2*Zdd21*Zuu21
     &      /sb-EE**2/MW**2*Mq**2*Pa11/SW**2*Sqrt2*Zdd22
     &      *Zuu22/sb)
      AAABR(1551) = 1D0/cb*(2D0*EE**2*Pa21/SW**2*Sqrt2*Zcc21
     &      *Zss21*sb-2D0*EE**2*Pa21/SW**2*Sqrt2*Zcc21*Zss21
     &      *sb**3+2D0*EE*Lambda/MW*Mc*Pa22/SW*Zcc22*Zss21
     &      /sb-2D0*EE*Lambda/MW*Mc*Pa22/SW*Zcc22*Zss21*sb
     &      +2D0*EE*Lambda/MW*Mq*Pa22/SW*Zcc21*Zss22*sb-EE
     &      **2/MW**2*Mc**2*Pa21/SW**2*Sqrt2*Zcc21*Zss21
     &      /sb+EE**2/MW**2*Mc**2*Pa21/SW**2*Sqrt2*Zcc21
     &      *Zss21*sb-EE**2/MW**2*Mq**2*Pa21/SW**2*Sqrt2
     &      *Zcc21*Zss21*sb-EE**2/MW**2*Mc*Mq*Pa21/SW**2
     &      *Sqrt2*Zcc22*Zss22/sb)
      AAABR(1552) = 1D0/cb*(2D0*EE**2*Pa21/SW**2*Sqrt2*Zb21
     &      *Zt21*sb-2D0*EE**2*Pa21/SW**2*Sqrt2*Zb21*Zt21
     &      *sb**3+2D0*EE*Lambda/MW*MbMM*Pa22/SW*Zb22*Zt21
     &      *sb+2D0*EE*Lambda/MW*MtMM*Pa22/SW*Zb21*Zt22/sb
     &      -2D0*EE*Lambda/MW*MtMM*Pa22/SW*Zb21*Zt22*sb-EE
     &      **2/MW**2*MbMM**2*Pa21/SW**2*Sqrt2*Zb21*Zt21
     &      *sb-EE**2/MW**2*MtMM**2*Pa21/SW**2*Sqrt2*Zb21
     &      *Zt21/sb+EE**2/MW**2*MtMM**2*Pa21/SW**2*Sqrt2
     &      *Zb21*Zt21*sb-EE**2/MW**2*MbMM*MtMM*Pa21/SW**2
     &      *Sqrt2*Zb22*Zt22/sb)
      AAABR(1553) = EE**2/SW**2*Sqrt2*Zcc11*Zh31*Zss21/cb
     &      -EE**2/SW**2*Sqrt2*Zcc11*Zh32*Zss21*sb-EE**2
     &      /SW**2*Sqrt2*Zcc11*Zh31*Zss21/cb*sb**2+EE**2
     &      /MW**2*Mc**2/SW**2*Sqrt2*Zcc11*Zh32*Zss21/sb
     &      -EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zcc11*Zh31*Zss21
     &      /cb-2D0*EE*Lambda/MW*Mc/SW*Zcc12*Zh33*Zss21/cb
     &      /sb+2D0*EE*Lambda/MW*Mc/SW*Zcc12*Zh33*Zss21/cb
     &      *sb+2D0*EE*Lambda/MW*Mq/SW*Zcc11*Zh33*Zss22/cb
     &      *sb+EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc12*Zh31
     &      *Zss22/cb-EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc12
     &      *Zh32*Zss22/sb
      AAABR(1554) = EE**2/SW**2*Sqrt2*Zb21*Zh31*Zt11/cb-EE
     &      **2/SW**2*Sqrt2*Zb21*Zh32*Zt11*sb-EE**2/SW**2
     &      *Sqrt2*Zb21*Zh31*Zt11/cb*sb**2-EE**2/MW**2*MbMM
     &      **2/SW**2*Sqrt2*Zb21*Zh31*Zt11/cb+EE**2/MW**2
     &      *MtMM**2/SW**2*Sqrt2*Zb21*Zh32*Zt11/sb+2D0*EE
     &      *Lambda/MW*MbMM/SW*Zb22*Zh33*Zt11/cb*sb-2D0*EE
     &      *Lambda/MW*MtMM/SW*Zb21*Zh33*Zt12/cb/sb+2D0*EE
     &      *Lambda/MW*MtMM/SW*Zb21*Zh33*Zt12/cb*sb+EE**2
     &      /MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb22*Zh31*Zt12/cb
     &      -EE**2/MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb22*Zh32
     &      *Zt12/sb
      AAABR(1555) = EE**2/SW**2*Sqrt2*Zb21*Zh11*Zt11/cb-EE
     &      **2/SW**2*Sqrt2*Zb21*Zh12*Zt11*sb-EE**2/SW**2
     &      *Sqrt2*Zb21*Zh11*Zt11/cb*sb**2-EE**2/MW**2*MbMM
     &      **2/SW**2*Sqrt2*Zb21*Zh11*Zt11/cb+EE**2/MW**2
     &      *MtMM**2/SW**2*Sqrt2*Zb21*Zh12*Zt11/sb+2D0*EE
     &      *Lambda/MW*MbMM/SW*Zb22*Zh13*Zt11/cb*sb-2D0*EE
     &      *Lambda/MW*MtMM/SW*Zb21*Zh13*Zt12/cb/sb+2D0*EE
     &      *Lambda/MW*MtMM/SW*Zb21*Zh13*Zt12/cb*sb+EE**2
     &      /MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb22*Zh11*Zt12/cb
     &      -EE**2/MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb22*Zh12
     &      *Zt12/sb
      AAABR(1556) = EE**2/SW**2*Sqrt2*Zcc11*Zh11*Zss21/cb
     &      -EE**2/SW**2*Sqrt2*Zcc11*Zh12*Zss21*sb-EE**2
     &      /SW**2*Sqrt2*Zcc11*Zh11*Zss21/cb*sb**2+EE**2
     &      /MW**2*Mc**2/SW**2*Sqrt2*Zcc11*Zh12*Zss21/sb
     &      -EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zcc11*Zh11*Zss21
     &      /cb-2D0*EE*Lambda/MW*Mc/SW*Zcc12*Zh13*Zss21/cb
     &      /sb+2D0*EE*Lambda/MW*Mc/SW*Zcc12*Zh13*Zss21/cb
     &      *sb+2D0*EE*Lambda/MW*Mq/SW*Zcc11*Zh13*Zss22/cb
     &      *sb+EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc12*Zh11
     &      *Zss22/cb-EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc12
     &      *Zh12*Zss22/sb
      AAABR(1557) = EE**2/SW**2*Sqrt2*Zdd21*Zh11*Zuu11/cb
     &      -EE**2/SW**2*Sqrt2*Zdd21*Zh12*Zuu11*sb-EE**2
     &      /SW**2*Sqrt2*Zdd21*Zh11*Zuu11/cb*sb**2-EE**2
     &      /MW**2*Mq**2/SW**2*Sqrt2*Zdd21*Zh11*Zuu11/cb
     &      +EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd21*Zh12*Zuu11
     &      /sb+EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd22*Zh11
     &      *Zuu12/cb-EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd22
     &      *Zh12*Zuu12/sb-2D0*EE*Lambda/MW*Mq/SW*Zdd21*Zh13
     &      *Zuu12/cb/sb+2D0*EE*Lambda/MW*Mq/SW*Zdd21*Zh13
     &      *Zuu12/cb*sb+2D0*EE*Lambda/MW*Mq/SW*Zdd22*Zh13
     &      *Zuu11/cb*sb
      AAABR(1558) = EE**2/SW**2*Sqrt2*Zb21*Zh21*Zt11/cb-EE
     &      **2/SW**2*Sqrt2*Zb21*Zh22*Zt11*sb-EE**2/SW**2
     &      *Sqrt2*Zb21*Zh21*Zt11/cb*sb**2-EE**2/MW**2*MbMM
     &      **2/SW**2*Sqrt2*Zb21*Zh21*Zt11/cb+EE**2/MW**2
     &      *MtMM**2/SW**2*Sqrt2*Zb21*Zh22*Zt11/sb+2D0*EE
     &      *Lambda/MW*MbMM/SW*Zb22*Zh23*Zt11/cb*sb-2D0*EE
     &      *Lambda/MW*MtMM/SW*Zb21*Zh23*Zt12/cb/sb+2D0*EE
     &      *Lambda/MW*MtMM/SW*Zb21*Zh23*Zt12/cb*sb+EE**2
     &      /MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb22*Zh21*Zt12/cb
     &      -EE**2/MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb22*Zh22
     &      *Zt12/sb
      AAABR(1559) = EE**2/SW**2*Sqrt2*Zdd21*Zh21*Zuu11/cb
     &      -EE**2/SW**2*Sqrt2*Zdd21*Zh22*Zuu11*sb-EE**2
     &      /SW**2*Sqrt2*Zdd21*Zh21*Zuu11/cb*sb**2-EE**2
     &      /MW**2*Mq**2/SW**2*Sqrt2*Zdd21*Zh21*Zuu11/cb
     &      +EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd21*Zh22*Zuu11
     &      /sb+EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd22*Zh21
     &      *Zuu12/cb-EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd22
     &      *Zh22*Zuu12/sb-2D0*EE*Lambda/MW*Mq/SW*Zdd21*Zh23
     &      *Zuu12/cb/sb+2D0*EE*Lambda/MW*Mq/SW*Zdd21*Zh23
     &      *Zuu12/cb*sb+2D0*EE*Lambda/MW*Mq/SW*Zdd22*Zh23
     &      *Zuu11/cb*sb
      AAABR(1560) = EE**2/SW**2*Sqrt2*Zdd21*Zh31*Zuu11/cb
     &      -EE**2/SW**2*Sqrt2*Zdd21*Zh32*Zuu11*sb-EE**2
     &      /SW**2*Sqrt2*Zdd21*Zh31*Zuu11/cb*sb**2-EE**2
     &      /MW**2*Mq**2/SW**2*Sqrt2*Zdd21*Zh31*Zuu11/cb
     &      +EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd21*Zh32*Zuu11
     &      /sb+EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd22*Zh31
     &      *Zuu12/cb-EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd22
     &      *Zh32*Zuu12/sb-2D0*EE*Lambda/MW*Mq/SW*Zdd21*Zh33
     &      *Zuu12/cb/sb+2D0*EE*Lambda/MW*Mq/SW*Zdd21*Zh33
     &      *Zuu12/cb*sb+2D0*EE*Lambda/MW*Mq/SW*Zdd22*Zh33
     &      *Zuu11/cb*sb
      AAABR(1561) = EE**2/SW**2*Sqrt2*Zcc11*Zh21*Zss21/cb
     &      -EE**2/SW**2*Sqrt2*Zcc11*Zh22*Zss21*sb-EE**2
     &      /SW**2*Sqrt2*Zcc11*Zh21*Zss21/cb*sb**2+EE**2
     &      /MW**2*Mc**2/SW**2*Sqrt2*Zcc11*Zh22*Zss21/sb
     &      -EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zcc11*Zh21*Zss21
     &      /cb-2D0*EE*Lambda/MW*Mc/SW*Zcc12*Zh23*Zss21/cb
     &      /sb+2D0*EE*Lambda/MW*Mc/SW*Zcc12*Zh23*Zss21/cb
     &      *sb+2D0*EE*Lambda/MW*Mq/SW*Zcc11*Zh23*Zss22/cb
     &      *sb+EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc12*Zh21
     &      *Zss22/cb-EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc12
     &      *Zh22*Zss22/sb
      AAABR(1562) = EE**2/SW**2*Sqrt2*Zcc21*Zh31*Zss21/cb
     &      -EE**2/SW**2*Sqrt2*Zcc21*Zh32*Zss21*sb-EE**2
     &      /SW**2*Sqrt2*Zcc21*Zh31*Zss21/cb*sb**2+EE**2
     &      /MW**2*Mc**2/SW**2*Sqrt2*Zcc21*Zh32*Zss21/sb
     &      -EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zcc21*Zh31*Zss21
     &      /cb-2D0*EE*Lambda/MW*Mc/SW*Zcc22*Zh33*Zss21/cb
     &      /sb+2D0*EE*Lambda/MW*Mc/SW*Zcc22*Zh33*Zss21/cb
     &      *sb+2D0*EE*Lambda/MW*Mq/SW*Zcc21*Zh33*Zss22/cb
     &      *sb+EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc22*Zh31
     &      *Zss22/cb-EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc22
     &      *Zh32*Zss22/sb
      AAABR(1563) = EE**2/SW**2*Sqrt2*Zb21*Zh11*Zt21/cb-EE
     &      **2/SW**2*Sqrt2*Zb21*Zh12*Zt21*sb-EE**2/SW**2
     &      *Sqrt2*Zb21*Zh11*Zt21/cb*sb**2-EE**2/MW**2*MbMM
     &      **2/SW**2*Sqrt2*Zb21*Zh11*Zt21/cb+EE**2/MW**2
     &      *MtMM**2/SW**2*Sqrt2*Zb21*Zh12*Zt21/sb+2D0*EE
     &      *Lambda/MW*MbMM/SW*Zb22*Zh13*Zt21/cb*sb-2D0*EE
     &      *Lambda/MW*MtMM/SW*Zb21*Zh13*Zt22/cb/sb+2D0*EE
     &      *Lambda/MW*MtMM/SW*Zb21*Zh13*Zt22/cb*sb+EE**2
     &      /MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb22*Zh11*Zt22/cb
     &      -EE**2/MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb22*Zh12
     &      *Zt22/sb
      AAABR(1564) = EE**2/SW**2*Sqrt2*Zb21*Zh31*Zt21/cb-EE
     &      **2/SW**2*Sqrt2*Zb21*Zh32*Zt21*sb-EE**2/SW**2
     &      *Sqrt2*Zb21*Zh31*Zt21/cb*sb**2-EE**2/MW**2*MbMM
     &      **2/SW**2*Sqrt2*Zb21*Zh31*Zt21/cb+EE**2/MW**2
     &      *MtMM**2/SW**2*Sqrt2*Zb21*Zh32*Zt21/sb+2D0*EE
     &      *Lambda/MW*MbMM/SW*Zb22*Zh33*Zt21/cb*sb-2D0*EE
     &      *Lambda/MW*MtMM/SW*Zb21*Zh33*Zt22/cb/sb+2D0*EE
     &      *Lambda/MW*MtMM/SW*Zb21*Zh33*Zt22/cb*sb+EE**2
     &      /MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb22*Zh31*Zt22/cb
     &      -EE**2/MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb22*Zh32
     &      *Zt22/sb
      AAABR(1565) = EE**2/SW**2*Sqrt2*Zdd21*Zh21*Zuu21/cb
     &      -EE**2/SW**2*Sqrt2*Zdd21*Zh22*Zuu21*sb-EE**2
     &      /SW**2*Sqrt2*Zdd21*Zh21*Zuu21/cb*sb**2-EE**2
     &      /MW**2*Mq**2/SW**2*Sqrt2*Zdd21*Zh21*Zuu21/cb
     &      +EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd21*Zh22*Zuu21
     &      /sb+EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd22*Zh21
     &      *Zuu22/cb-EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd22
     &      *Zh22*Zuu22/sb-2D0*EE*Lambda/MW*Mq/SW*Zdd21*Zh23
     &      *Zuu22/cb/sb+2D0*EE*Lambda/MW*Mq/SW*Zdd21*Zh23
     &      *Zuu22/cb*sb+2D0*EE*Lambda/MW*Mq/SW*Zdd22*Zh23
     &      *Zuu21/cb*sb
      AAABR(1566) = EE**2/SW**2*Sqrt2*Zcc21*Zh21*Zss21/cb
     &      -EE**2/SW**2*Sqrt2*Zcc21*Zh22*Zss21*sb-EE**2
     &      /SW**2*Sqrt2*Zcc21*Zh21*Zss21/cb*sb**2+EE**2
     &      /MW**2*Mc**2/SW**2*Sqrt2*Zcc21*Zh22*Zss21/sb
     &      -EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zcc21*Zh21*Zss21
     &      /cb-2D0*EE*Lambda/MW*Mc/SW*Zcc22*Zh23*Zss21/cb
     &      /sb+2D0*EE*Lambda/MW*Mc/SW*Zcc22*Zh23*Zss21/cb
     &      *sb+2D0*EE*Lambda/MW*Mq/SW*Zcc21*Zh23*Zss22/cb
     &      *sb+EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc22*Zh21
     &      *Zss22/cb-EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc22
     &      *Zh22*Zss22/sb
      AAABR(1567) = EE**2/SW**2*Sqrt2*Zdd21*Zh31*Zuu21/cb
     &      -EE**2/SW**2*Sqrt2*Zdd21*Zh32*Zuu21*sb-EE**2
     &      /SW**2*Sqrt2*Zdd21*Zh31*Zuu21/cb*sb**2-EE**2
     &      /MW**2*Mq**2/SW**2*Sqrt2*Zdd21*Zh31*Zuu21/cb
     &      +EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd21*Zh32*Zuu21
     &      /sb+EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd22*Zh31
     &      *Zuu22/cb-EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd22
     &      *Zh32*Zuu22/sb-2D0*EE*Lambda/MW*Mq/SW*Zdd21*Zh33
     &      *Zuu22/cb/sb+2D0*EE*Lambda/MW*Mq/SW*Zdd21*Zh33
     &      *Zuu22/cb*sb+2D0*EE*Lambda/MW*Mq/SW*Zdd22*Zh33
     &      *Zuu21/cb*sb
      AAABR(1568) = EE**2/SW**2*Sqrt2*Zb21*Zh21*Zt21/cb-EE
     &      **2/SW**2*Sqrt2*Zb21*Zh22*Zt21*sb-EE**2/SW**2
     &      *Sqrt2*Zb21*Zh21*Zt21/cb*sb**2-EE**2/MW**2*MbMM
     &      **2/SW**2*Sqrt2*Zb21*Zh21*Zt21/cb+EE**2/MW**2
     &      *MtMM**2/SW**2*Sqrt2*Zb21*Zh22*Zt21/sb+2D0*EE
     &      *Lambda/MW*MbMM/SW*Zb22*Zh23*Zt21/cb*sb-2D0*EE
     &      *Lambda/MW*MtMM/SW*Zb21*Zh23*Zt22/cb/sb+2D0*EE
     &      *Lambda/MW*MtMM/SW*Zb21*Zh23*Zt22/cb*sb+EE**2
     &      /MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb22*Zh21*Zt22/cb
     &      -EE**2/MW**2*MbMM*MtMM/SW**2*Sqrt2*Zb22*Zh22
     &      *Zt22/sb
      AAABR(1569) = EE**2/SW**2*Sqrt2*Zcc21*Zh11*Zss21/cb
     &      -EE**2/SW**2*Sqrt2*Zcc21*Zh12*Zss21*sb-EE**2
     &      /SW**2*Sqrt2*Zcc21*Zh11*Zss21/cb*sb**2+EE**2
     &      /MW**2*Mc**2/SW**2*Sqrt2*Zcc21*Zh12*Zss21/sb
     &      -EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zcc21*Zh11*Zss21
     &      /cb-2D0*EE*Lambda/MW*Mc/SW*Zcc22*Zh13*Zss21/cb
     &      /sb+2D0*EE*Lambda/MW*Mc/SW*Zcc22*Zh13*Zss21/cb
     &      *sb+2D0*EE*Lambda/MW*Mq/SW*Zcc21*Zh13*Zss22/cb
     &      *sb+EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc22*Zh11
     &      *Zss22/cb-EE**2/MW**2*Mc*Mq/SW**2*Sqrt2*Zcc22
     &      *Zh12*Zss22/sb
      AAABR(1570) = EE**2/SW**2*Sqrt2*Zdd21*Zh11*Zuu21/cb
     &      -EE**2/SW**2*Sqrt2*Zdd21*Zh12*Zuu21*sb-EE**2
     &      /SW**2*Sqrt2*Zdd21*Zh11*Zuu21/cb*sb**2-EE**2
     &      /MW**2*Mq**2/SW**2*Sqrt2*Zdd21*Zh11*Zuu21/cb
     &      +EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd21*Zh12*Zuu21
     &      /sb+EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd22*Zh11
     &      *Zuu22/cb-EE**2/MW**2*Mq**2/SW**2*Sqrt2*Zdd22
     &      *Zh12*Zuu22/sb-2D0*EE*Lambda/MW*Mq/SW*Zdd21*Zh13
     &      *Zuu22/cb/sb+2D0*EE*Lambda/MW*Mq/SW*Zdd21*Zh13
     &      *Zuu22/cb*sb+2D0*EE*Lambda/MW*Mq/SW*Zdd22*Zh13
     &      *Zuu21/cb*sb
      AAABR(1571) = EE**2*(Zss21**2/CW**2+2D0*Zss22**2/CW
     &      **2+3D0*Zss21**2/SW**2-2D0*Zss21**2/CW**2*sb
     &      **2-4D0*Zss22**2/CW**2*sb**2-6D0*Mq**2/MW**2
     &      /SW**2-6D0*Zss21**2/SW**2*sb**2)
      AAABR(1572) = EE**2*(Zb21**2/CW**2+2D0*Zb22**2/CW**2
     &      +3D0*Zb21**2/SW**2-2D0*Zb21**2/CW**2*sb**2-4D0
     &      *Zb22**2/CW**2*sb**2-6D0*MbMM**2/MW**2/SW**2
     &      -6D0*Zb21**2/SW**2*sb**2)
      AAABR(1573) = EE**2*(Zdd21**2/CW**2+2D0*Zdd22**2/CW
     &      **2+3D0*Zdd21**2/SW**2-2D0*Zdd21**2/CW**2*sb
     &      **2-4D0*Zdd22**2/CW**2*sb**2-6D0*Mq**2/MW**2
     &      /SW**2-6D0*Zdd21**2/SW**2*sb**2)
      AAABR(1574) = tB*(EE**2/CW**2*Pa21*Zb21**2+2D0*EE**2
     &      /CW**2*Pa21*Zb22**2+3D0*EE**2*Pa21/SW**2*Zb21
     &      **2-EE**2/CW**2*Pa21*Zb21**2*sb**2-2D0*EE**2
     &      /CW**2*Pa21*Zb22**2*sb**2-3D0*EE**2/MW**2*MbMM
     &      **2*Pa21/SW**2-3D0*EE**2*Pa21/SW**2*Zb21**2*sb
     &      **2+3D0*EE*Lambda/MW*MbMM*Pa22/SW*Sqrt2*Zb21
     &      *Zb22)
      AAABR(1575) = tB*(EE**2/CW**2*Pa11*Zdd21**2+2D0*EE
     &      **2/CW**2*Pa11*Zdd22**2+3D0*EE**2*Pa11/SW**2
     &      *Zdd21**2-EE**2/CW**2*Pa11*Zdd21**2*sb**2-2D0
     &      *EE**2/CW**2*Pa11*Zdd22**2*sb**2-3D0*EE**2/MW
     &      **2*Mq**2*Pa11/SW**2-3D0*EE**2*Pa11/SW**2*Zdd21
     &      **2*sb**2+3D0*EE*Lambda/MW*Mq*Pa12/SW*Sqrt2*Zdd21
     &      *Zdd22)
      AAABR(1576) = tB*(EE**2/CW**2*Pa21*Zss21**2+2D0*EE
     &      **2/CW**2*Pa21*Zss22**2+3D0*EE**2*Pa21/SW**2
     &      *Zss21**2-EE**2/CW**2*Pa21*Zss21**2*sb**2-2D0
     &      *EE**2/CW**2*Pa21*Zss22**2*sb**2-3D0*EE**2/MW
     &      **2*Mq**2*Pa21/SW**2-3D0*EE**2*Pa21/SW**2*Zss21
     &      **2*sb**2+3D0*EE*Lambda/MW*Mq*Pa22/SW*Sqrt2*Zss21
     &      *Zss22)
      AAABR(1577) = tB*(EE**2/CW**2*Pa11*Zss21**2+2D0*EE
     &      **2/CW**2*Pa11*Zss22**2+3D0*EE**2*Pa11/SW**2
     &      *Zss21**2-EE**2/CW**2*Pa11*Zss21**2*sb**2-2D0
     &      *EE**2/CW**2*Pa11*Zss22**2*sb**2-3D0*EE**2/MW
     &      **2*Mq**2*Pa11/SW**2-3D0*EE**2*Pa11/SW**2*Zss21
     &      **2*sb**2+3D0*EE*Lambda/MW*Mq*Pa12/SW*Sqrt2*Zss21
     &      *Zss22)
      AAABR(1578) = tB*(EE**2/CW**2*Pa11*Zb21**2+2D0*EE**2
     &      /CW**2*Pa11*Zb22**2+3D0*EE**2*Pa11/SW**2*Zb21
     &      **2-EE**2/CW**2*Pa11*Zb21**2*sb**2-2D0*EE**2
     &      /CW**2*Pa11*Zb22**2*sb**2-3D0*EE**2/MW**2*MbMM
     &      **2*Pa11/SW**2-3D0*EE**2*Pa11/SW**2*Zb21**2*sb
     &      **2+3D0*EE*Lambda/MW*MbMM*Pa12/SW*Sqrt2*Zb21
     &      *Zb22)
      AAABR(1579) = tB*(EE**2/CW**2*Pa21*Zdd21**2+2D0*EE
     &      **2/CW**2*Pa21*Zdd22**2+3D0*EE**2*Pa21/SW**2
     &      *Zdd21**2-EE**2/CW**2*Pa21*Zdd21**2*sb**2-2D0
     &      *EE**2/CW**2*Pa21*Zdd22**2*sb**2-3D0*EE**2/MW
     &      **2*Mq**2*Pa21/SW**2-3D0*EE**2*Pa21/SW**2*Zdd21
     &      **2*sb**2+3D0*EE*Lambda/MW*Mq*Pa22/SW*Sqrt2*Zdd21
     &      *Zdd22)
      AAABR(1580) = 1D0/cb**2*(EE**2/CW**2*Pa11*Pa21*Zss21
     &      **2+2D0*EE**2/CW**2*Pa11*Pa21*Zss22**2+3D0*EE
     &      **2*Pa11*Pa21/SW**2*Zss21**2-3D0*EE**2/CW**2
     &      *Pa11*Pa21*Zss21**2*sb**2+2D0*EE**2/CW**2*Pa11
     &      *Pa21*Zss21**2*sb**4-6D0*EE**2/CW**2*Pa11*Pa21
     &      *Zss22**2*sb**2+4D0*EE**2/CW**2*Pa11*Pa21*Zss22
     &      **2*sb**4-9D0*EE**2*Pa11*Pa21/SW**2*Zss21**2
     &      *sb**2+6D0*EE**2*Pa11*Pa21/SW**2*Zss21**2*sb
     &      **4+6D0*EE**2/MW**2*Mq**2*Pa11*Pa21/SW**2*sb
     &      **2+6D0*EE*Lambda/MW*Mq*Pa11*Pa22/SW*Sqrt2*Zss21
     &      *Zss22+6D0*EE*Lambda/MW*Mq*Pa12*Pa21/SW*Sqrt2
     &      *Zss21*Zss22-6D0*EE*Lambda/MW*Mq*Pa11*Pa22/SW
     &      *Sqrt2*Zss21*Zss22*sb**2-6D0*EE*Lambda/MW*Mq
     &      *Pa12*Pa21/SW*Sqrt2*Zss21*Zss22*sb**2)
      AAABR(1581) = 1D0/cb**2*(EE**2/CW**2*Pa21**2*Zdd21
     &      **2+2D0*EE**2/CW**2*Pa21**2*Zdd22**2+3D0*EE**2
     &      *Pa21**2/SW**2*Zdd21**2-3D0*EE**2/CW**2*Pa21
     &      **2*Zdd21**2*sb**2+2D0*EE**2/CW**2*Pa21**2*Zdd21
     &      **2*sb**4-6D0*EE**2/CW**2*Pa21**2*Zdd22**2*sb
     &      **2+4D0*EE**2/CW**2*Pa21**2*Zdd22**2*sb**4-9D0
     &      *EE**2*Pa21**2/SW**2*Zdd21**2*sb**2+6D0*EE**2
     &      *Pa21**2/SW**2*Zdd21**2*sb**4+6D0*EE**2/MW**2
     &      *Mq**2*Pa21**2/SW**2*sb**2+12D0*EE*Lambda/MW
     &      *Mq*Pa21*Pa22/SW*Sqrt2*Zdd21*Zdd22-12D0*EE*Lambda
     &      /MW*Mq*Pa21*Pa22/SW*Sqrt2*Zdd21*Zdd22*sb**2)
      AAABR(1582) = 1D0/cb**2*(EE**2/CW**2*Pa11*Pa21*Zdd21
     &      **2+2D0*EE**2/CW**2*Pa11*Pa21*Zdd22**2+3D0*EE
     &      **2*Pa11*Pa21/SW**2*Zdd21**2-3D0*EE**2/CW**2
     &      *Pa11*Pa21*Zdd21**2*sb**2+2D0*EE**2/CW**2*Pa11
     &      *Pa21*Zdd21**2*sb**4-6D0*EE**2/CW**2*Pa11*Pa21
     &      *Zdd22**2*sb**2+4D0*EE**2/CW**2*Pa11*Pa21*Zdd22
     &      **2*sb**4-9D0*EE**2*Pa11*Pa21/SW**2*Zdd21**2
     &      *sb**2+6D0*EE**2*Pa11*Pa21/SW**2*Zdd21**2*sb
     &      **4+6D0*EE**2/MW**2*Mq**2*Pa11*Pa21/SW**2*sb
     &      **2+6D0*EE*Lambda/MW*Mq*Pa11*Pa22/SW*Sqrt2*Zdd21
     &      *Zdd22+6D0*EE*Lambda/MW*Mq*Pa12*Pa21/SW*Sqrt2
     &      *Zdd21*Zdd22-6D0*EE*Lambda/MW*Mq*Pa11*Pa22/SW
     &      *Sqrt2*Zdd21*Zdd22*sb**2-6D0*EE*Lambda/MW*Mq
     &      *Pa12*Pa21/SW*Sqrt2*Zdd21*Zdd22*sb**2)
      AAABR(1583) = 1D0/cb**2*(EE**2/CW**2*Pa11**2*Zdd21
     &      **2+2D0*EE**2/CW**2*Pa11**2*Zdd22**2+3D0*EE**2
     &      *Pa11**2/SW**2*Zdd21**2-3D0*EE**2/CW**2*Pa11
     &      **2*Zdd21**2*sb**2+2D0*EE**2/CW**2*Pa11**2*Zdd21
     &      **2*sb**4-6D0*EE**2/CW**2*Pa11**2*Zdd22**2*sb
     &      **2+4D0*EE**2/CW**2*Pa11**2*Zdd22**2*sb**4-9D0
     &      *EE**2*Pa11**2/SW**2*Zdd21**2*sb**2+6D0*EE**2
     &      *Pa11**2/SW**2*Zdd21**2*sb**4+6D0*EE**2/MW**2
     &      *Mq**2*Pa11**2/SW**2*sb**2+12D0*EE*Lambda/MW
     &      *Mq*Pa11*Pa12/SW*Sqrt2*Zdd21*Zdd22-12D0*EE*Lambda
     &      /MW*Mq*Pa11*Pa12/SW*Sqrt2*Zdd21*Zdd22*sb**2)
      AAABR(1584) = 1D0/cb**2*(EE**2/CW**2*Pa21**2*Zb21**2
     &      +2D0*EE**2/CW**2*Pa21**2*Zb22**2+3D0*EE**2*Pa21
     &      **2/SW**2*Zb21**2-3D0*EE**2/CW**2*Pa21**2*Zb21
     &      **2*sb**2+2D0*EE**2/CW**2*Pa21**2*Zb21**2*sb
     &      **4-6D0*EE**2/CW**2*Pa21**2*Zb22**2*sb**2+4D0
     &      *EE**2/CW**2*Pa21**2*Zb22**2*sb**4-9D0*EE**2
     &      *Pa21**2/SW**2*Zb21**2*sb**2+6D0*EE**2*Pa21**2
     &      /SW**2*Zb21**2*sb**4+6D0*EE**2/MW**2*MbMM**2
     &      *Pa21**2/SW**2*sb**2+12D0*EE*Lambda/MW*MbMM*Pa21
     &      *Pa22/SW*Sqrt2*Zb21*Zb22-12D0*EE*Lambda/MW*MbMM
     &      *Pa21*Pa22/SW*Sqrt2*Zb21*Zb22*sb**2)
      AAABR(1585) = 1D0/cb**2*(EE**2/CW**2*Pa11*Pa21*Zb21
     &      **2+2D0*EE**2/CW**2*Pa11*Pa21*Zb22**2+3D0*EE
     &      **2*Pa11*Pa21/SW**2*Zb21**2-3D0*EE**2/CW**2*Pa11
     &      *Pa21*Zb21**2*sb**2+2D0*EE**2/CW**2*Pa11*Pa21
     &      *Zb21**2*sb**4-6D0*EE**2/CW**2*Pa11*Pa21*Zb22
     &      **2*sb**2+4D0*EE**2/CW**2*Pa11*Pa21*Zb22**2*sb
     &      **4-9D0*EE**2*Pa11*Pa21/SW**2*Zb21**2*sb**2+6D0
     &      *EE**2*Pa11*Pa21/SW**2*Zb21**2*sb**4+6D0*EE**2
     &      /MW**2*MbMM**2*Pa11*Pa21/SW**2*sb**2+6D0*EE*Lambda
     &      /MW*MbMM*Pa11*Pa22/SW*Sqrt2*Zb21*Zb22+6D0*EE
     &      *Lambda/MW*MbMM*Pa12*Pa21/SW*Sqrt2*Zb21*Zb22
     &      -6D0*EE*Lambda/MW*MbMM*Pa11*Pa22/SW*Sqrt2*Zb21
     &      *Zb22*sb**2-6D0*EE*Lambda/MW*MbMM*Pa12*Pa21/SW
     &      *Sqrt2*Zb21*Zb22*sb**2)
      AAABR(1586) = 1D0/cb**2*(EE**2/CW**2*Pa11**2*Zb21**2
     &      +2D0*EE**2/CW**2*Pa11**2*Zb22**2+3D0*EE**2*Pa11
     &      **2/SW**2*Zb21**2-3D0*EE**2/CW**2*Pa11**2*Zb21
     &      **2*sb**2+2D0*EE**2/CW**2*Pa11**2*Zb21**2*sb
     &      **4-6D0*EE**2/CW**2*Pa11**2*Zb22**2*sb**2+4D0
     &      *EE**2/CW**2*Pa11**2*Zb22**2*sb**4-9D0*EE**2
     &      *Pa11**2/SW**2*Zb21**2*sb**2+6D0*EE**2*Pa11**2
     &      /SW**2*Zb21**2*sb**4+6D0*EE**2/MW**2*MbMM**2
     &      *Pa11**2/SW**2*sb**2+12D0*EE*Lambda/MW*MbMM*Pa11
     &      *Pa12/SW*Sqrt2*Zb21*Zb22-12D0*EE*Lambda/MW*MbMM
     &      *Pa11*Pa12/SW*Sqrt2*Zb21*Zb22*sb**2)
      AAABR(1587) = 1D0/cb**2*(EE**2/CW**2*Pa21**2*Zss21
     &      **2+2D0*EE**2/CW**2*Pa21**2*Zss22**2+3D0*EE**2
     &      *Pa21**2/SW**2*Zss21**2-3D0*EE**2/CW**2*Pa21
     &      **2*Zss21**2*sb**2+2D0*EE**2/CW**2*Pa21**2*Zss21
     &      **2*sb**4-6D0*EE**2/CW**2*Pa21**2*Zss22**2*sb
     &      **2+4D0*EE**2/CW**2*Pa21**2*Zss22**2*sb**4-9D0
     &      *EE**2*Pa21**2/SW**2*Zss21**2*sb**2+6D0*EE**2
     &      *Pa21**2/SW**2*Zss21**2*sb**4+6D0*EE**2/MW**2
     &      *Mq**2*Pa21**2/SW**2*sb**2+12D0*EE*Lambda/MW
     &      *Mq*Pa21*Pa22/SW*Sqrt2*Zss21*Zss22-12D0*EE*Lambda
     &      /MW*Mq*Pa21*Pa22/SW*Sqrt2*Zss21*Zss22*sb**2)
      AAABR(1588) = 1D0/cb**2*(EE**2/CW**2*Pa11**2*Zss21
     &      **2+2D0*EE**2/CW**2*Pa11**2*Zss22**2+3D0*EE**2
     &      *Pa11**2/SW**2*Zss21**2-3D0*EE**2/CW**2*Pa11
     &      **2*Zss21**2*sb**2+2D0*EE**2/CW**2*Pa11**2*Zss21
     &      **2*sb**4-6D0*EE**2/CW**2*Pa11**2*Zss22**2*sb
     &      **2+4D0*EE**2/CW**2*Pa11**2*Zss22**2*sb**4-9D0
     &      *EE**2*Pa11**2/SW**2*Zss21**2*sb**2+6D0*EE**2
     &      *Pa11**2/SW**2*Zss21**2*sb**4+6D0*EE**2/MW**2
     &      *Mq**2*Pa11**2/SW**2*sb**2+12D0*EE*Lambda/MW
     &      *Mq*Pa11*Pa12/SW*Sqrt2*Zss21*Zss22-12D0*EE*Lambda
     &      /MW*Mq*Pa11*Pa12/SW*Sqrt2*Zss21*Zss22*sb**2)
      AAABR(1589) = EE**2/CW**2*Zdd21**2*Zh11*Zh31/cb**2
     &      -EE**2/CW**2*Zdd21**2*Zh12*Zh32/cb**2+2D0*EE
     &      **2/CW**2*Zdd22**2*Zh11*Zh31/cb**2-2D0*EE**2
     &      /CW**2*Zdd22**2*Zh12*Zh32/cb**2+3D0*EE**2/SW
     &      **2*Zdd21**2*Zh11*Zh31/cb**2-3D0*EE**2/SW**2
     &      *Zdd21**2*Zh12*Zh32/cb**2-EE**2/CW**2*Zdd21**2
     &      *Zh11*Zh31/cb**2*sb**2+EE**2/CW**2*Zdd21**2*Zh12
     &      *Zh32/cb**2*sb**2-2D0*EE**2/CW**2*Zdd22**2*Zh11
     &      *Zh31/cb**2*sb**2+2D0*EE**2/CW**2*Zdd22**2*Zh12
     &      *Zh32/cb**2*sb**2-6D0*EE**2/MW**2*Mq**2/SW**2
     &      *Zh11*Zh31/cb**2-3D0*EE**2/SW**2*Zdd21**2*Zh11
     &      *Zh31/cb**2*sb**2+3D0*EE**2/SW**2*Zdd21**2*Zh12
     &      *Zh32/cb**2*sb**2+6D0*EE*Lambda/MW*Mq/SW*Sqrt2
     &      *Zdd21*Zdd22*Zh12*Zh33/cb+6D0*EE*Lambda/MW*Mq
     &      /SW*Sqrt2*Zdd21*Zdd22*Zh13*Zh32/cb
      AAABR(1590) = EE**2/CW**2*Zb21**2*Zh11*Zh31/cb**2-EE
     &      **2/CW**2*Zb21**2*Zh12*Zh32/cb**2+2D0*EE**2/CW
     &      **2*Zb22**2*Zh11*Zh31/cb**2-2D0*EE**2/CW**2*Zb22
     &      **2*Zh12*Zh32/cb**2+3D0*EE**2/SW**2*Zb21**2*Zh11
     &      *Zh31/cb**2-3D0*EE**2/SW**2*Zb21**2*Zh12*Zh32
     &      /cb**2-EE**2/CW**2*Zb21**2*Zh11*Zh31/cb**2*sb
     &      **2+EE**2/CW**2*Zb21**2*Zh12*Zh32/cb**2*sb**2
     &      -2D0*EE**2/CW**2*Zb22**2*Zh11*Zh31/cb**2*sb**2
     &      +2D0*EE**2/CW**2*Zb22**2*Zh12*Zh32/cb**2*sb**2
     &      -6D0*EE**2/MW**2*MbMM**2/SW**2*Zh11*Zh31/cb**2
     &      -3D0*EE**2/SW**2*Zb21**2*Zh11*Zh31/cb**2*sb**2
     &      +3D0*EE**2/SW**2*Zb21**2*Zh12*Zh32/cb**2*sb**2
     &      +6D0*EE*Lambda/MW*MbMM/SW*Sqrt2*Zb21*Zb22*Zh12
     &      *Zh33/cb+6D0*EE*Lambda/MW*MbMM/SW*Sqrt2*Zb21
     &      *Zb22*Zh13*Zh32/cb
      AAABR(1591) = EE**2/CW**2*Zdd21**2*Zh11*Zh21/cb**2
     &      -EE**2/CW**2*Zdd21**2*Zh12*Zh22/cb**2+2D0*EE
     &      **2/CW**2*Zdd22**2*Zh11*Zh21/cb**2-2D0*EE**2
     &      /CW**2*Zdd22**2*Zh12*Zh22/cb**2+3D0*EE**2/SW
     &      **2*Zdd21**2*Zh11*Zh21/cb**2-3D0*EE**2/SW**2
     &      *Zdd21**2*Zh12*Zh22/cb**2-EE**2/CW**2*Zdd21**2
     &      *Zh11*Zh21/cb**2*sb**2+EE**2/CW**2*Zdd21**2*Zh12
     &      *Zh22/cb**2*sb**2-2D0*EE**2/CW**2*Zdd22**2*Zh11
     &      *Zh21/cb**2*sb**2+2D0*EE**2/CW**2*Zdd22**2*Zh12
     &      *Zh22/cb**2*sb**2-6D0*EE**2/MW**2*Mq**2/SW**2
     &      *Zh11*Zh21/cb**2-3D0*EE**2/SW**2*Zdd21**2*Zh11
     &      *Zh21/cb**2*sb**2+3D0*EE**2/SW**2*Zdd21**2*Zh12
     &      *Zh22/cb**2*sb**2+6D0*EE*Lambda/MW*Mq/SW*Sqrt2
     &      *Zdd21*Zdd22*Zh12*Zh23/cb+6D0*EE*Lambda/MW*Mq
     &      /SW*Sqrt2*Zdd21*Zdd22*Zh13*Zh22/cb
      AAABR(1592) = EE**2/CW**2*Zb21**2*Zh11*Zh21/cb**2-EE
     &      **2/CW**2*Zb21**2*Zh12*Zh22/cb**2+2D0*EE**2/CW
     &      **2*Zb22**2*Zh11*Zh21/cb**2-2D0*EE**2/CW**2*Zb22
     &      **2*Zh12*Zh22/cb**2+3D0*EE**2/SW**2*Zb21**2*Zh11
     &      *Zh21/cb**2-3D0*EE**2/SW**2*Zb21**2*Zh12*Zh22
     &      /cb**2-EE**2/CW**2*Zb21**2*Zh11*Zh21/cb**2*sb
     &      **2+EE**2/CW**2*Zb21**2*Zh12*Zh22/cb**2*sb**2
     &      -2D0*EE**2/CW**2*Zb22**2*Zh11*Zh21/cb**2*sb**2
     &      +2D0*EE**2/CW**2*Zb22**2*Zh12*Zh22/cb**2*sb**2
     &      -6D0*EE**2/MW**2*MbMM**2/SW**2*Zh11*Zh21/cb**2
     &      -3D0*EE**2/SW**2*Zb21**2*Zh11*Zh21/cb**2*sb**2
     &      +3D0*EE**2/SW**2*Zb21**2*Zh12*Zh22/cb**2*sb**2
     &      +6D0*EE*Lambda/MW*MbMM/SW*Sqrt2*Zb21*Zb22*Zh12
     &      *Zh23/cb+6D0*EE*Lambda/MW*MbMM/SW*Sqrt2*Zb21
     &      *Zb22*Zh13*Zh22/cb
      AAABR(1593) = EE**2/CW**2*Zdd21**2*Zh21*Zh31/cb**2
     &      -EE**2/CW**2*Zdd21**2*Zh22*Zh32/cb**2+2D0*EE
     &      **2/CW**2*Zdd22**2*Zh21*Zh31/cb**2-2D0*EE**2
     &      /CW**2*Zdd22**2*Zh22*Zh32/cb**2+3D0*EE**2/SW
     &      **2*Zdd21**2*Zh21*Zh31/cb**2-3D0*EE**2/SW**2
     &      *Zdd21**2*Zh22*Zh32/cb**2-EE**2/CW**2*Zdd21**2
     &      *Zh21*Zh31/cb**2*sb**2+EE**2/CW**2*Zdd21**2*Zh22
     &      *Zh32/cb**2*sb**2-2D0*EE**2/CW**2*Zdd22**2*Zh21
     &      *Zh31/cb**2*sb**2+2D0*EE**2/CW**2*Zdd22**2*Zh22
     &      *Zh32/cb**2*sb**2-6D0*EE**2/MW**2*Mq**2/SW**2
     &      *Zh21*Zh31/cb**2-3D0*EE**2/SW**2*Zdd21**2*Zh21
     &      *Zh31/cb**2*sb**2+3D0*EE**2/SW**2*Zdd21**2*Zh22
     &      *Zh32/cb**2*sb**2+6D0*EE*Lambda/MW*Mq/SW*Sqrt2
     &      *Zdd21*Zdd22*Zh22*Zh33/cb+6D0*EE*Lambda/MW*Mq
     &      /SW*Sqrt2*Zdd21*Zdd22*Zh23*Zh32/cb
      AAABR(1594) = EE**2/CW**2*Zdd21**2*Zh11**2/cb**2-EE
     &      **2/CW**2*Zdd21**2*Zh12**2/cb**2+2D0*EE**2/CW
     &      **2*Zdd22**2*Zh11**2/cb**2-2D0*EE**2/CW**2*Zdd22
     &      **2*Zh12**2/cb**2+3D0*EE**2/SW**2*Zdd21**2*Zh11
     &      **2/cb**2-3D0*EE**2/SW**2*Zdd21**2*Zh12**2/cb
     &      **2-EE**2/CW**2*Zdd21**2*Zh11**2/cb**2*sb**2
     &      +EE**2/CW**2*Zdd21**2*Zh12**2/cb**2*sb**2-2D0
     &      *EE**2/CW**2*Zdd22**2*Zh11**2/cb**2*sb**2+2D0
     &      *EE**2/CW**2*Zdd22**2*Zh12**2/cb**2*sb**2-6D0
     &      *EE**2/MW**2*Mq**2/SW**2*Zh11**2/cb**2-3D0*EE
     &      **2/SW**2*Zdd21**2*Zh11**2/cb**2*sb**2+3D0*EE
     &      **2/SW**2*Zdd21**2*Zh12**2/cb**2*sb**2+12D0*EE
     &      *Lambda/MW*Mq/SW*Sqrt2*Zdd21*Zdd22*Zh12*Zh13/cb
      AAABR(1595) = EE**2/CW**2*Zb21**2*Zh11**2/cb**2-EE
     &      **2/CW**2*Zb21**2*Zh12**2/cb**2+2D0*EE**2/CW
     &      **2*Zb22**2*Zh11**2/cb**2-2D0*EE**2/CW**2*Zb22
     &      **2*Zh12**2/cb**2+3D0*EE**2/SW**2*Zb21**2*Zh11
     &      **2/cb**2-3D0*EE**2/SW**2*Zb21**2*Zh12**2/cb
     &      **2-EE**2/CW**2*Zb21**2*Zh11**2/cb**2*sb**2+EE
     &      **2/CW**2*Zb21**2*Zh12**2/cb**2*sb**2-2D0*EE
     &      **2/CW**2*Zb22**2*Zh11**2/cb**2*sb**2+2D0*EE
     &      **2/CW**2*Zb22**2*Zh12**2/cb**2*sb**2-6D0*EE
     &      **2/MW**2*MbMM**2/SW**2*Zh11**2/cb**2-3D0*EE
     &      **2/SW**2*Zb21**2*Zh11**2/cb**2*sb**2+3D0*EE
     &      **2/SW**2*Zb21**2*Zh12**2/cb**2*sb**2+12D0*EE
     &      *Lambda/MW*MbMM/SW*Sqrt2*Zb21*Zb22*Zh12*Zh13/cb
      AAABR(1596) = EE**2/CW**2*Zb21**2*Zh21*Zh31/cb**2-EE
     &      **2/CW**2*Zb21**2*Zh22*Zh32/cb**2+2D0*EE**2/CW
     &      **2*Zb22**2*Zh21*Zh31/cb**2-2D0*EE**2/CW**2*Zb22
     &      **2*Zh22*Zh32/cb**2+3D0*EE**2/SW**2*Zb21**2*Zh21
     &      *Zh31/cb**2-3D0*EE**2/SW**2*Zb21**2*Zh22*Zh32
     &      /cb**2-EE**2/CW**2*Zb21**2*Zh21*Zh31/cb**2*sb
     &      **2+EE**2/CW**2*Zb21**2*Zh22*Zh32/cb**2*sb**2
     &      -2D0*EE**2/CW**2*Zb22**2*Zh21*Zh31/cb**2*sb**2
     &      +2D0*EE**2/CW**2*Zb22**2*Zh22*Zh32/cb**2*sb**2
     &      -6D0*EE**2/MW**2*MbMM**2/SW**2*Zh21*Zh31/cb**2
     &      -3D0*EE**2/SW**2*Zb21**2*Zh21*Zh31/cb**2*sb**2
     &      +3D0*EE**2/SW**2*Zb21**2*Zh22*Zh32/cb**2*sb**2
     &      +6D0*EE*Lambda/MW*MbMM/SW*Sqrt2*Zb21*Zb22*Zh22
     &      *Zh33/cb+6D0*EE*Lambda/MW*MbMM/SW*Sqrt2*Zb21
     &      *Zb22*Zh23*Zh32/cb
      AAABR(1597) = EE**2/CW**2*Zh11**2*Zss21**2/cb**2+2D0
     &      *EE**2/CW**2*Zh11**2*Zss22**2/cb**2-EE**2/CW
     &      **2*Zh12**2*Zss21**2/cb**2-2D0*EE**2/CW**2*Zh12
     &      **2*Zss22**2/cb**2+3D0*EE**2/SW**2*Zh11**2*Zss21
     &      **2/cb**2-3D0*EE**2/SW**2*Zh12**2*Zss21**2/cb
     &      **2-EE**2/CW**2*Zh11**2*Zss21**2/cb**2*sb**2
     &      -2D0*EE**2/CW**2*Zh11**2*Zss22**2/cb**2*sb**2
     &      +EE**2/CW**2*Zh12**2*Zss21**2/cb**2*sb**2+2D0
     &      *EE**2/CW**2*Zh12**2*Zss22**2/cb**2*sb**2-6D0
     &      *EE**2/MW**2*Mq**2/SW**2*Zh11**2/cb**2-3D0*EE
     &      **2/SW**2*Zh11**2*Zss21**2/cb**2*sb**2+3D0*EE
     &      **2/SW**2*Zh12**2*Zss21**2/cb**2*sb**2+12D0*EE
     &      *Lambda/MW*Mq/SW*Sqrt2*Zh12*Zh13*Zss21*Zss22/cb
      AAABR(1598) = EE**2/CW**2*Zdd21**2*Zh21**2/cb**2-EE
     &      **2/CW**2*Zdd21**2*Zh22**2/cb**2+2D0*EE**2/CW
     &      **2*Zdd22**2*Zh21**2/cb**2-2D0*EE**2/CW**2*Zdd22
     &      **2*Zh22**2/cb**2+3D0*EE**2/SW**2*Zdd21**2*Zh21
     &      **2/cb**2-3D0*EE**2/SW**2*Zdd21**2*Zh22**2/cb
     &      **2-EE**2/CW**2*Zdd21**2*Zh21**2/cb**2*sb**2
     &      +EE**2/CW**2*Zdd21**2*Zh22**2/cb**2*sb**2-2D0
     &      *EE**2/CW**2*Zdd22**2*Zh21**2/cb**2*sb**2+2D0
     &      *EE**2/CW**2*Zdd22**2*Zh22**2/cb**2*sb**2-6D0
     &      *EE**2/MW**2*Mq**2/SW**2*Zh21**2/cb**2-3D0*EE
     &      **2/SW**2*Zdd21**2*Zh21**2/cb**2*sb**2+3D0*EE
     &      **2/SW**2*Zdd21**2*Zh22**2/cb**2*sb**2+12D0*EE
     &      *Lambda/MW*Mq/SW*Sqrt2*Zdd21*Zdd22*Zh22*Zh23/cb
      AAABR(1599) = EE**2/CW**2*Zh11*Zh31*Zss21**2/cb**2
     &      +2D0*EE**2/CW**2*Zh11*Zh31*Zss22**2/cb**2-EE
     &      **2/CW**2*Zh12*Zh32*Zss21**2/cb**2-2D0*EE**2
     &      /CW**2*Zh12*Zh32*Zss22**2/cb**2+3D0*EE**2/SW
     &      **2*Zh11*Zh31*Zss21**2/cb**2-3D0*EE**2/SW**2
     &      *Zh12*Zh32*Zss21**2/cb**2-EE**2/CW**2*Zh11*Zh31
     &      *Zss21**2/cb**2*sb**2-2D0*EE**2/CW**2*Zh11*Zh31
     &      *Zss22**2/cb**2*sb**2+EE**2/CW**2*Zh12*Zh32*Zss21
     &      **2/cb**2*sb**2+2D0*EE**2/CW**2*Zh12*Zh32*Zss22
     &      **2/cb**2*sb**2-6D0*EE**2/MW**2*Mq**2/SW**2*Zh11
     &      *Zh31/cb**2-3D0*EE**2/SW**2*Zh11*Zh31*Zss21**2
     &      /cb**2*sb**2+3D0*EE**2/SW**2*Zh12*Zh32*Zss21
     &      **2/cb**2*sb**2+6D0*EE*Lambda/MW*Mq/SW*Sqrt2
     &      *Zh12*Zh33*Zss21*Zss22/cb+6D0*EE*Lambda/MW*Mq
     &      /SW*Sqrt2*Zh13*Zh32*Zss21*Zss22/cb
      AAABR(1600) = EE**2/CW**2*Zdd21**2*Zh31**2/cb**2-EE
     &      **2/CW**2*Zdd21**2*Zh32**2/cb**2+2D0*EE**2/CW
     &      **2*Zdd22**2*Zh31**2/cb**2-2D0*EE**2/CW**2*Zdd22
     &      **2*Zh32**2/cb**2+3D0*EE**2/SW**2*Zdd21**2*Zh31
     &      **2/cb**2-3D0*EE**2/SW**2*Zdd21**2*Zh32**2/cb
     &      **2-EE**2/CW**2*Zdd21**2*Zh31**2/cb**2*sb**2
     &      +EE**2/CW**2*Zdd21**2*Zh32**2/cb**2*sb**2-2D0
     &      *EE**2/CW**2*Zdd22**2*Zh31**2/cb**2*sb**2+2D0
     &      *EE**2/CW**2*Zdd22**2*Zh32**2/cb**2*sb**2-6D0
     &      *EE**2/MW**2*Mq**2/SW**2*Zh31**2/cb**2-3D0*EE
     &      **2/SW**2*Zdd21**2*Zh31**2/cb**2*sb**2+3D0*EE
     &      **2/SW**2*Zdd21**2*Zh32**2/cb**2*sb**2+12D0*EE
     &      *Lambda/MW*Mq/SW*Sqrt2*Zdd21*Zdd22*Zh32*Zh33/cb
      AAABR(1601) = EE**2/CW**2*Zh21*Zh31*Zss21**2/cb**2
     &      +2D0*EE**2/CW**2*Zh21*Zh31*Zss22**2/cb**2-EE
     &      **2/CW**2*Zh22*Zh32*Zss21**2/cb**2-2D0*EE**2
     &      /CW**2*Zh22*Zh32*Zss22**2/cb**2+3D0*EE**2/SW
     &      **2*Zh21*Zh31*Zss21**2/cb**2-3D0*EE**2/SW**2
     &      *Zh22*Zh32*Zss21**2/cb**2-EE**2/CW**2*Zh21*Zh31
     &      *Zss21**2/cb**2*sb**2-2D0*EE**2/CW**2*Zh21*Zh31
     &      *Zss22**2/cb**2*sb**2+EE**2/CW**2*Zh22*Zh32*Zss21
     &      **2/cb**2*sb**2+2D0*EE**2/CW**2*Zh22*Zh32*Zss22
     &      **2/cb**2*sb**2-6D0*EE**2/MW**2*Mq**2/SW**2*Zh21
     &      *Zh31/cb**2-3D0*EE**2/SW**2*Zh21*Zh31*Zss21**2
     &      /cb**2*sb**2+3D0*EE**2/SW**2*Zh22*Zh32*Zss21
     &      **2/cb**2*sb**2+6D0*EE*Lambda/MW*Mq/SW*Sqrt2
     &      *Zh22*Zh33*Zss21*Zss22/cb+6D0*EE*Lambda/MW*Mq
     &      /SW*Sqrt2*Zh23*Zh32*Zss21*Zss22/cb
      AAABR(1602) = EE**2/CW**2*Zb21**2*Zh21**2/cb**2-EE
     &      **2/CW**2*Zb21**2*Zh22**2/cb**2+2D0*EE**2/CW
     &      **2*Zb22**2*Zh21**2/cb**2-2D0*EE**2/CW**2*Zb22
     &      **2*Zh22**2/cb**2+3D0*EE**2/SW**2*Zb21**2*Zh21
     &      **2/cb**2-3D0*EE**2/SW**2*Zb21**2*Zh22**2/cb
     &      **2-EE**2/CW**2*Zb21**2*Zh21**2/cb**2*sb**2+EE
     &      **2/CW**2*Zb21**2*Zh22**2/cb**2*sb**2-2D0*EE
     &      **2/CW**2*Zb22**2*Zh21**2/cb**2*sb**2+2D0*EE
     &      **2/CW**2*Zb22**2*Zh22**2/cb**2*sb**2-6D0*EE
     &      **2/MW**2*MbMM**2/SW**2*Zh21**2/cb**2-3D0*EE
     &      **2/SW**2*Zb21**2*Zh21**2/cb**2*sb**2+3D0*EE
     &      **2/SW**2*Zb21**2*Zh22**2/cb**2*sb**2+12D0*EE
     &      *Lambda/MW*MbMM/SW*Sqrt2*Zb21*Zb22*Zh22*Zh23/cb
      AAABR(1603) = EE**2/CW**2*Zh31**2*Zss21**2/cb**2+2D0
     &      *EE**2/CW**2*Zh31**2*Zss22**2/cb**2-EE**2/CW
     &      **2*Zh32**2*Zss21**2/cb**2-2D0*EE**2/CW**2*Zh32
     &      **2*Zss22**2/cb**2+3D0*EE**2/SW**2*Zh31**2*Zss21
     &      **2/cb**2-3D0*EE**2/SW**2*Zh32**2*Zss21**2/cb
     &      **2-EE**2/CW**2*Zh31**2*Zss21**2/cb**2*sb**2
     &      -2D0*EE**2/CW**2*Zh31**2*Zss22**2/cb**2*sb**2
     &      +EE**2/CW**2*Zh32**2*Zss21**2/cb**2*sb**2+2D0
     &      *EE**2/CW**2*Zh32**2*Zss22**2/cb**2*sb**2-6D0
     &      *EE**2/MW**2*Mq**2/SW**2*Zh31**2/cb**2-3D0*EE
     &      **2/SW**2*Zh31**2*Zss21**2/cb**2*sb**2+3D0*EE
     &      **2/SW**2*Zh32**2*Zss21**2/cb**2*sb**2+12D0*EE
     &      *Lambda/MW*Mq/SW*Sqrt2*Zh32*Zh33*Zss21*Zss22/cb
      AAABR(1604) = EE**2/CW**2*Zh21**2*Zss21**2/cb**2+2D0
     &      *EE**2/CW**2*Zh21**2*Zss22**2/cb**2-EE**2/CW
     &      **2*Zh22**2*Zss21**2/cb**2-2D0*EE**2/CW**2*Zh22
     &      **2*Zss22**2/cb**2+3D0*EE**2/SW**2*Zh21**2*Zss21
     &      **2/cb**2-3D0*EE**2/SW**2*Zh22**2*Zss21**2/cb
     &      **2-EE**2/CW**2*Zh21**2*Zss21**2/cb**2*sb**2
     &      -2D0*EE**2/CW**2*Zh21**2*Zss22**2/cb**2*sb**2
     &      +EE**2/CW**2*Zh22**2*Zss21**2/cb**2*sb**2+2D0
     &      *EE**2/CW**2*Zh22**2*Zss22**2/cb**2*sb**2-6D0
     &      *EE**2/MW**2*Mq**2/SW**2*Zh21**2/cb**2-3D0*EE
     &      **2/SW**2*Zh21**2*Zss21**2/cb**2*sb**2+3D0*EE
     &      **2/SW**2*Zh22**2*Zss21**2/cb**2*sb**2+12D0*EE
     &      *Lambda/MW*Mq/SW*Sqrt2*Zh22*Zh23*Zss21*Zss22/cb
      AAABR(1605) = EE**2/CW**2*Zb21**2*Zh31**2/cb**2-EE
     &      **2/CW**2*Zb21**2*Zh32**2/cb**2+2D0*EE**2/CW
     &      **2*Zb22**2*Zh31**2/cb**2-2D0*EE**2/CW**2*Zb22
     &      **2*Zh32**2/cb**2+3D0*EE**2/SW**2*Zb21**2*Zh31
     &      **2/cb**2-3D0*EE**2/SW**2*Zb21**2*Zh32**2/cb
     &      **2-EE**2/CW**2*Zb21**2*Zh31**2/cb**2*sb**2+EE
     &      **2/CW**2*Zb21**2*Zh32**2/cb**2*sb**2-2D0*EE
     &      **2/CW**2*Zb22**2*Zh31**2/cb**2*sb**2+2D0*EE
     &      **2/CW**2*Zb22**2*Zh32**2/cb**2*sb**2-6D0*EE
     &      **2/MW**2*MbMM**2/SW**2*Zh31**2/cb**2-3D0*EE
     &      **2/SW**2*Zb21**2*Zh31**2/cb**2*sb**2+3D0*EE
     &      **2/SW**2*Zb21**2*Zh32**2/cb**2*sb**2+12D0*EE
     &      *Lambda/MW*MbMM/SW*Sqrt2*Zb21*Zb22*Zh32*Zh33/cb
      AAABR(1606) = EE**2/CW**2*Zh11*Zh21*Zss21**2/cb**2
     &      +2D0*EE**2/CW**2*Zh11*Zh21*Zss22**2/cb**2-EE
     &      **2/CW**2*Zh12*Zh22*Zss21**2/cb**2-2D0*EE**2
     &      /CW**2*Zh12*Zh22*Zss22**2/cb**2+3D0*EE**2/SW
     &      **2*Zh11*Zh21*Zss21**2/cb**2-3D0*EE**2/SW**2
     &      *Zh12*Zh22*Zss21**2/cb**2-EE**2/CW**2*Zh11*Zh21
     &      *Zss21**2/cb**2*sb**2-2D0*EE**2/CW**2*Zh11*Zh21
     &      *Zss22**2/cb**2*sb**2+EE**2/CW**2*Zh12*Zh22*Zss21
     &      **2/cb**2*sb**2+2D0*EE**2/CW**2*Zh12*Zh22*Zss22
     &      **2/cb**2*sb**2-6D0*EE**2/MW**2*Mq**2/SW**2*Zh11
     &      *Zh21/cb**2-3D0*EE**2/SW**2*Zh11*Zh21*Zss21**2
     &      /cb**2*sb**2+3D0*EE**2/SW**2*Zh12*Zh22*Zss21
     &      **2/cb**2*sb**2+6D0*EE*Lambda/MW*Mq/SW*Sqrt2
     &      *Zh12*Zh23*Zss21*Zss22/cb+6D0*EE*Lambda/MW*Mq
     &      /SW*Sqrt2*Zh13*Zh22*Zss21*Zss22/cb
      AAABR(1607) = EE**2*(Zcc11**2/CW**2*Zuu11**2-4D0*Zcc11
     &      **2/CW**2*Zuu12**2-4D0*Zcc12**2/CW**2*Zuu11**2
     &      +16D0*Zcc12**2/CW**2*Zuu12**2+9D0*Zcc11**2/SW
     &      **2*Zuu11**2+36D0*Mc/MW**2*Mq/SW**2*Zcc11*Zcc12
     &      *Zuu11*Zuu12/sb**2)
      AAABR(1608) = GG**2*(Zcc11**4+Zcc12**4+2D0*Zcc11**2
     &      *Zcc12**2)
      AAABR(1609) = EE**2*(Zcc11**4/CW**2+16D0*Zcc12**4/CW
     &      **2+9D0*Zcc11**4/SW**2-8D0*Zcc11**2/CW**2*Zcc12
     &      **2+36D0*Mc**2/MW**2/SW**2*Zcc11**2*Zcc12**2
     &      /sb**2)
      AAABR(1610) = EE**2*(Zt11**2/CW**2*Zuu11**2-4D0*Zt11
     &      **2/CW**2*Zuu12**2-4D0*Zt12**2/CW**2*Zuu11**2
     &      +16D0*Zt12**2/CW**2*Zuu12**2+9D0*Zt11**2/SW**2
     &      *Zuu11**2+36D0*Mq/MW**2*MtMM/SW**2*Zt11*Zt12
     &      *Zuu11*Zuu12/sb**2)
      AAABR(1611) = EE**2*(Zcc11**2/CW**2*Zt11**2-4D0*Zcc11
     &      **2/CW**2*Zt12**2-4D0*Zcc12**2/CW**2*Zt11**2
     &      +16D0*Zcc12**2/CW**2*Zt12**2+9D0*Zcc11**2/SW
     &      **2*Zt11**2+36D0*Mc/MW**2*MtMM/SW**2*Zcc11*Zcc12
     &      *Zt11*Zt12/sb**2)
      AAABR(1612) = GG**2*(Zuu11**4+Zuu12**4+2D0*Zuu11**2
     &      *Zuu12**2)
      AAABR(1613) = EE**2*(Zuu11**4/CW**2+16D0*Zuu12**4/CW
     &      **2+9D0*Zuu11**4/SW**2-8D0*Zuu11**2/CW**2*Zuu12
     &      **2+36D0*Mq**2/MW**2/SW**2*Zuu11**2*Zuu12**2
     &      /sb**2)
      AAABR(1614) = GG**2*(Zt11**4+Zt12**4+2D0*Zt11**2*Zt12
     &      **2)
      AAABR(1615) = EE**2*(Zt11**4/CW**2+16D0*Zt12**4/CW
     &      **2+9D0*Zt11**4/SW**2-8D0*Zt11**2/CW**2*Zt12
     &      **2+36D0*MtMM**2/MW**2/SW**2*Zt11**2*Zt12**2
     &      /sb**2)
      AAABR(1616) = EE**2*(Zt11/CW**2*Zt21*Zuu11**2-4D0*Zt11
     &      /CW**2*Zt21*Zuu12**2-4D0*Zt12/CW**2*Zt22*Zuu11
     &      **2+16D0*Zt12/CW**2*Zt22*Zuu12**2+9D0*Zt11/SW
     &      **2*Zt21*Zuu11**2+18D0*Mq/MW**2*MtMM/SW**2*Zt11
     &      *Zt22*Zuu11*Zuu12/sb**2+18D0*Mq/MW**2*MtMM/SW
     &      **2*Zt12*Zt21*Zuu11*Zuu12/sb**2)
      AAABR(1617) = GG**2*(Zcc11**3*Zcc21+Zcc12**3*Zcc22
     &      +Zcc11*Zcc12**2*Zcc21+Zcc11**2*Zcc12*Zcc22)
      AAABR(1618) = EE**2*(Zcc11**3/CW**2*Zcc21+16D0*Zcc12
     &      **3/CW**2*Zcc22+9D0*Zcc11**3/SW**2*Zcc21-4D0
     &      *Zcc11/CW**2*Zcc12**2*Zcc21-4D0*Zcc11**2/CW**2
     &      *Zcc12*Zcc22+18D0*Mc**2/MW**2/SW**2*Zcc11*Zcc12
     &      **2*Zcc21/sb**2+18D0*Mc**2/MW**2/SW**2*Zcc11
     &      **2*Zcc12*Zcc22/sb**2)
      AAABR(1619) = GG**2*(Zt11**3*Zt21+Zt12**3*Zt22+Zt11
     &      *Zt12**2*Zt21+Zt11**2*Zt12*Zt22)
      AAABR(1620) = EE**2*(Zt11**3/CW**2*Zt21+16D0*Zt12**3
     &      /CW**2*Zt22+9D0*Zt11**3/SW**2*Zt21-4D0*Zt11/CW
     &      **2*Zt12**2*Zt21-4D0*Zt11**2/CW**2*Zt12*Zt22
     &      +18D0*MtMM**2/MW**2/SW**2*Zt11*Zt12**2*Zt21/sb
     &      **2+18D0*MtMM**2/MW**2/SW**2*Zt11**2*Zt12*Zt22
     &      /sb**2)
      AAABR(1621) = EE**2*(Zt11**2/CW**2*Zuu11*Zuu21-4D0
     &      *Zt11**2/CW**2*Zuu12*Zuu22-4D0*Zt12**2/CW**2
     &      *Zuu11*Zuu21+16D0*Zt12**2/CW**2*Zuu12*Zuu22+9D0
     &      *Zt11**2/SW**2*Zuu11*Zuu21+18D0*Mq/MW**2*MtMM
     &      /SW**2*Zt11*Zt12*Zuu11*Zuu22/sb**2+18D0*Mq/MW
     &      **2*MtMM/SW**2*Zt11*Zt12*Zuu12*Zuu21/sb**2)
      AAABR(1622) = EE**2*(Zcc11**2/CW**2*Zt11*Zt21-4D0*Zcc11
     &      **2/CW**2*Zt12*Zt22-4D0*Zcc12**2/CW**2*Zt11*Zt21
     &      +16D0*Zcc12**2/CW**2*Zt12*Zt22+9D0*Zcc11**2/SW
     &      **2*Zt11*Zt21+18D0*Mc/MW**2*MtMM/SW**2*Zcc11
     &      *Zcc12*Zt11*Zt22/sb**2+18D0*Mc/MW**2*MtMM/SW
     &      **2*Zcc11*Zcc12*Zt12*Zt21/sb**2)
      AAABR(1623) = EE**2*(Zcc11**2/CW**2*Zuu11*Zuu21-4D0
     &      *Zcc11**2/CW**2*Zuu12*Zuu22-4D0*Zcc12**2/CW**2
     &      *Zuu11*Zuu21+16D0*Zcc12**2/CW**2*Zuu12*Zuu22
     &      +9D0*Zcc11**2/SW**2*Zuu11*Zuu21+18D0*Mc/MW**2
     &      *Mq/SW**2*Zcc11*Zcc12*Zuu11*Zuu22/sb**2+18D0
     &      *Mc/MW**2*Mq/SW**2*Zcc11*Zcc12*Zuu12*Zuu21/sb
     &      **2)
      AAABR(1624) = EE**2*(Zcc11/CW**2*Zcc21*Zuu11**2-4D0
     &      *Zcc11/CW**2*Zcc21*Zuu12**2-4D0*Zcc12/CW**2*Zcc22
     &      *Zuu11**2+16D0*Zcc12/CW**2*Zcc22*Zuu12**2+9D0
     &      *Zcc11/SW**2*Zcc21*Zuu11**2+18D0*Mc/MW**2*Mq
     &      /SW**2*Zcc11*Zcc22*Zuu11*Zuu12/sb**2+18D0*Mc
     &      /MW**2*Mq/SW**2*Zcc12*Zcc21*Zuu11*Zuu12/sb**2)
      AAABR(1625) = EE**2*(Zcc11/CW**2*Zcc21*Zt11**2-4D0
     &      *Zcc11/CW**2*Zcc21*Zt12**2-4D0*Zcc12/CW**2*Zcc22
     &      *Zt11**2+16D0*Zcc12/CW**2*Zcc22*Zt12**2+9D0*Zcc11
     &      /SW**2*Zcc21*Zt11**2+18D0*Mc/MW**2*MtMM/SW**2
     &      *Zcc11*Zcc22*Zt11*Zt12/sb**2+18D0*Mc/MW**2*MtMM
     &      /SW**2*Zcc12*Zcc21*Zt11*Zt12/sb**2)
      AAABR(1626) = GG**2*(Zuu11**3*Zuu21+Zuu12**3*Zuu22
     &      +Zuu11*Zuu12**2*Zuu21+Zuu11**2*Zuu12*Zuu22)
      AAABR(1627) = EE**2*(Zuu11**3/CW**2*Zuu21+16D0*Zuu12
     &      **3/CW**2*Zuu22+9D0*Zuu11**3/SW**2*Zuu21-4D0
     &      *Zuu11/CW**2*Zuu12**2*Zuu21-4D0*Zuu11**2/CW**2
     &      *Zuu12*Zuu22+18D0*Mq**2/MW**2/SW**2*Zuu11*Zuu12
     &      **2*Zuu21/sb**2+18D0*Mq**2/MW**2/SW**2*Zuu11
     &      **2*Zuu12*Zuu22/sb**2)
      AAABR(1628) = GG**2*(Zuu11**2*Zuu21**2+Zuu12**2*Zuu22
     &      **2+2D0*Zuu11*Zuu12*Zuu21*Zuu22)
      AAABR(1629) = EE**2*(Zuu11**2/CW**2*Zuu21**2+16D0*Zuu12
     &      **2/CW**2*Zuu22**2+9D0*Zuu11**2/SW**2*Zuu21**2
     &      -8D0*Zuu11/CW**2*Zuu12*Zuu21*Zuu22+36D0*Mq**2
     &      /MW**2/SW**2*Zuu11*Zuu12*Zuu21*Zuu22/sb**2)
      AAABR(1630) = EE**2*(Zcc11/CW**2*Zcc21*Zuu11*Zuu21
     &      -4D0*Zcc11/CW**2*Zcc21*Zuu12*Zuu22-4D0*Zcc12
     &      /CW**2*Zcc22*Zuu11*Zuu21+16D0*Zcc12/CW**2*Zcc22
     &      *Zuu12*Zuu22+9D0*Zcc11/SW**2*Zcc21*Zuu11*Zuu21
     &      +18D0*Mc/MW**2*Mq/SW**2*Zcc11*Zcc22*Zuu12*Zuu21
     &      /sb**2+18D0*Mc/MW**2*Mq/SW**2*Zcc12*Zcc21*Zuu11
     &      *Zuu22/sb**2)
      AAABR(1631) = GG**2*(Zcc11**2*Zcc21**2+Zcc12**2*Zcc22
     &      **2+2D0*Zcc11*Zcc12*Zcc21*Zcc22)
      AAABR(1632) = EE**2*(Zcc11**2/CW**2*Zcc21**2+16D0*Zcc12
     &      **2/CW**2*Zcc22**2+9D0*Zcc11**2/SW**2*Zcc21**2
     &      -8D0*Zcc11/CW**2*Zcc12*Zcc21*Zcc22+36D0*Mc**2
     &      /MW**2/SW**2*Zcc11*Zcc12*Zcc21*Zcc22/sb**2)
      AAABR(1633) = GG**2*(Zt11**2*Zt21**2+Zt12**2*Zt22**2
     &      +2D0*Zt11*Zt12*Zt21*Zt22)
      AAABR(1634) = EE**2*(Zt11**2/CW**2*Zt21**2+16D0*Zt12
     &      **2/CW**2*Zt22**2+9D0*Zt11**2/SW**2*Zt21**2-8D0
     &      *Zt11/CW**2*Zt12*Zt21*Zt22+36D0*MtMM**2/MW**2
     &      /SW**2*Zt11*Zt12*Zt21*Zt22/sb**2)
      AAABR(1635) = EE**2*(Zt11/CW**2*Zt21*Zuu11*Zuu21-4D0
     &      *Zt11/CW**2*Zt21*Zuu12*Zuu22-4D0*Zt12/CW**2*Zt22
     &      *Zuu11*Zuu21+16D0*Zt12/CW**2*Zt22*Zuu12*Zuu22
     &      +9D0*Zt11/SW**2*Zt21*Zuu11*Zuu21+18D0*Mq/MW**2
     &      *MtMM/SW**2*Zt11*Zt22*Zuu12*Zuu21/sb**2+18D0
     &      *Mq/MW**2*MtMM/SW**2*Zt12*Zt21*Zuu11*Zuu22/sb
     &      **2)
      AAABR(1636) = EE**2*(Zcc11/CW**2*Zcc21*Zt11*Zt21-4D0
     &      *Zcc11/CW**2*Zcc21*Zt12*Zt22-4D0*Zcc12/CW**2
     &      *Zcc22*Zt11*Zt21+16D0*Zcc12/CW**2*Zcc22*Zt12
     &      *Zt22+9D0*Zcc11/SW**2*Zcc21*Zt11*Zt21+18D0*Mc
     &      /MW**2*MtMM/SW**2*Zcc11*Zcc22*Zt12*Zt21/sb**2
     &      +18D0*Mc/MW**2*MtMM/SW**2*Zcc12*Zcc21*Zt11*Zt22
     &      /sb**2)
      AAABR(1637) = EE**2*(Zcc11/CW**2*Zcc21*Zt11*Zt21-4D0
     &      *Zcc11/CW**2*Zcc21*Zt12*Zt22-4D0*Zcc12/CW**2
     &      *Zcc22*Zt11*Zt21+16D0*Zcc12/CW**2*Zcc22*Zt12
     &      *Zt22+9D0*Zcc11/SW**2*Zcc21*Zt11*Zt21+18D0*Mc
     &      /MW**2*MtMM/SW**2*Zcc11*Zcc22*Zt11*Zt22/sb**2
     &      +18D0*Mc/MW**2*MtMM/SW**2*Zcc12*Zcc21*Zt12*Zt21
     &      /sb**2)
      AAABR(1638) = GG**2*(Zt11**2*Zt21**2+Zt11**2*Zt22**2
     &      +Zt12**2*Zt21**2+Zt12**2*Zt22**2)
      AAABR(1639) = EE**2*(Zt11**2/CW**2*Zt21**2-4D0*Zt11
     &      **2/CW**2*Zt22**2-4D0*Zt12**2/CW**2*Zt21**2+16D0
     &      *Zt12**2/CW**2*Zt22**2+9D0*Zt11**2/SW**2*Zt21
     &      **2+36D0*MtMM**2/MW**2/SW**2*Zt11*Zt12*Zt21*Zt22
     &      /sb**2)
      AAABR(1640) = EE**2*(Zt11**2/CW**2*Zt21**2+16D0*Zt12
     &      **2/CW**2*Zt22**2+9D0*Zt11**2/SW**2*Zt21**2-8D0
     &      *Zt11/CW**2*Zt12*Zt21*Zt22+18D0*MtMM**2/MW**2
     &      /SW**2*Zt11**2*Zt22**2/sb**2+18D0*MtMM**2/MW
     &      **2/SW**2*Zt12**2*Zt21**2/sb**2)
      AAABR(1641) = EE**2*(Zcc11/CW**2*Zcc21*Zuu11*Zuu21
     &      -4D0*Zcc11/CW**2*Zcc21*Zuu12*Zuu22-4D0*Zcc12
     &      /CW**2*Zcc22*Zuu11*Zuu21+16D0*Zcc12/CW**2*Zcc22
     &      *Zuu12*Zuu22+9D0*Zcc11/SW**2*Zcc21*Zuu11*Zuu21
     &      +18D0*Mc/MW**2*Mq/SW**2*Zcc11*Zcc22*Zuu11*Zuu22
     &      /sb**2+18D0*Mc/MW**2*Mq/SW**2*Zcc12*Zcc21*Zuu12
     &      *Zuu21/sb**2)
      AAABR(1642) = EE**2*(Zcc21**2/CW**2*Zuu11**2-4D0*Zcc21
     &      **2/CW**2*Zuu12**2-4D0*Zcc22**2/CW**2*Zuu11**2
     &      +16D0*Zcc22**2/CW**2*Zuu12**2+9D0*Zcc21**2/SW
     &      **2*Zuu11**2+36D0*Mc/MW**2*Mq/SW**2*Zcc21*Zcc22
     &      *Zuu11*Zuu12/sb**2)
      AAABR(1643) = EE**2*(Zcc21**2/CW**2*Zt11**2-4D0*Zcc21
     &      **2/CW**2*Zt12**2-4D0*Zcc22**2/CW**2*Zt11**2
     &      +16D0*Zcc22**2/CW**2*Zt12**2+9D0*Zcc21**2/SW
     &      **2*Zt11**2+36D0*Mc/MW**2*MtMM/SW**2*Zcc21*Zcc22
     &      *Zt11*Zt12/sb**2)
      AAABR(1644) = EE**2*(Zcc11**2/CW**2*Zuu21**2-4D0*Zcc11
     &      **2/CW**2*Zuu22**2-4D0*Zcc12**2/CW**2*Zuu21**2
     &      +16D0*Zcc12**2/CW**2*Zuu22**2+9D0*Zcc11**2/SW
     &      **2*Zuu21**2+36D0*Mc/MW**2*Mq/SW**2*Zcc11*Zcc12
     &      *Zuu21*Zuu22/sb**2)
      AAABR(1645) = EE**2*(Zt11/CW**2*Zt21*Zuu11*Zuu21-4D0
     &      *Zt11/CW**2*Zt21*Zuu12*Zuu22-4D0*Zt12/CW**2*Zt22
     &      *Zuu11*Zuu21+16D0*Zt12/CW**2*Zt22*Zuu12*Zuu22
     &      +9D0*Zt11/SW**2*Zt21*Zuu11*Zuu21+18D0*Mq/MW**2
     &      *MtMM/SW**2*Zt11*Zt22*Zuu11*Zuu22/sb**2+18D0
     &      *Mq/MW**2*MtMM/SW**2*Zt12*Zt21*Zuu12*Zuu21/sb
     &      **2)
      AAABR(1646) = EE**2*(Zcc11**2/CW**2*Zt21**2-4D0*Zcc11
     &      **2/CW**2*Zt22**2-4D0*Zcc12**2/CW**2*Zt21**2
     &      +16D0*Zcc12**2/CW**2*Zt22**2+9D0*Zcc11**2/SW
     &      **2*Zt21**2+36D0*Mc/MW**2*MtMM/SW**2*Zcc11*Zcc12
     &      *Zt21*Zt22/sb**2)
      AAABR(1647) = EE**2*(Zt21**2/CW**2*Zuu11**2-4D0*Zt21
     &      **2/CW**2*Zuu12**2-4D0*Zt22**2/CW**2*Zuu11**2
     &      +16D0*Zt22**2/CW**2*Zuu12**2+9D0*Zt21**2/SW**2
     &      *Zuu11**2+36D0*Mq/MW**2*MtMM/SW**2*Zt21*Zt22
     &      *Zuu11*Zuu12/sb**2)
      AAABR(1648) = GG**2*(Zcc11**2*Zcc21**2+Zcc11**2*Zcc22
     &      **2+Zcc12**2*Zcc21**2+Zcc12**2*Zcc22**2)
      AAABR(1649) = EE**2*(Zcc11**2/CW**2*Zcc21**2-4D0*Zcc11
     &      **2/CW**2*Zcc22**2-4D0*Zcc12**2/CW**2*Zcc21**2
     &      +16D0*Zcc12**2/CW**2*Zcc22**2+9D0*Zcc11**2/SW
     &      **2*Zcc21**2+36D0*Mc**2/MW**2/SW**2*Zcc11*Zcc12
     &      *Zcc21*Zcc22/sb**2)
      AAABR(1650) = EE**2*(Zcc11**2/CW**2*Zcc21**2+16D0*Zcc12
     &      **2/CW**2*Zcc22**2+9D0*Zcc11**2/SW**2*Zcc21**2
     &      -8D0*Zcc11/CW**2*Zcc12*Zcc21*Zcc22+18D0*Mc**2
     &      /MW**2/SW**2*Zcc11**2*Zcc22**2/sb**2+18D0*Mc
     &      **2/MW**2/SW**2*Zcc12**2*Zcc21**2/sb**2)
      AAABR(1651) = EE**2*(Zt11**2/CW**2*Zuu21**2-4D0*Zt11
     &      **2/CW**2*Zuu22**2-4D0*Zt12**2/CW**2*Zuu21**2
     &      +16D0*Zt12**2/CW**2*Zuu22**2+9D0*Zt11**2/SW**2
     &      *Zuu21**2+36D0*Mq/MW**2*MtMM/SW**2*Zt11*Zt12
     &      *Zuu21*Zuu22/sb**2)
      AAABR(1652) = GG**2*(Zuu11**2*Zuu21**2+Zuu11**2*Zuu22
     &      **2+Zuu12**2*Zuu21**2+Zuu12**2*Zuu22**2)
      AAABR(1653) = EE**2*(Zuu11**2/CW**2*Zuu21**2-4D0*Zuu11
     &      **2/CW**2*Zuu22**2-4D0*Zuu12**2/CW**2*Zuu21**2
     &      +16D0*Zuu12**2/CW**2*Zuu22**2+9D0*Zuu11**2/SW
     &      **2*Zuu21**2+36D0*Mq**2/MW**2/SW**2*Zuu11*Zuu12
     &      *Zuu21*Zuu22/sb**2)
      AAABR(1654) = EE**2*(Zuu11**2/CW**2*Zuu21**2+16D0*Zuu12
     &      **2/CW**2*Zuu22**2+9D0*Zuu11**2/SW**2*Zuu21**2
     &      -8D0*Zuu11/CW**2*Zuu12*Zuu21*Zuu22+18D0*Mq**2
     &      /MW**2/SW**2*Zuu11**2*Zuu22**2/sb**2+18D0*Mq
     &      **2/MW**2/SW**2*Zuu12**2*Zuu21**2/sb**2)
      AAABR(1655) = EE**2*(Zcc21**2/CW**2*Zuu11*Zuu21-4D0
     &      *Zcc21**2/CW**2*Zuu12*Zuu22-4D0*Zcc22**2/CW**2
     &      *Zuu11*Zuu21+16D0*Zcc22**2/CW**2*Zuu12*Zuu22
     &      +9D0*Zcc21**2/SW**2*Zuu11*Zuu21+18D0*Mc/MW**2
     &      *Mq/SW**2*Zcc21*Zcc22*Zuu11*Zuu22/sb**2+18D0
     &      *Mc/MW**2*Mq/SW**2*Zcc21*Zcc22*Zuu12*Zuu21/sb
     &      **2)
      AAABR(1656) = GG**2*(Zcc11*Zcc21**3+Zcc12*Zcc22**3
     &      +Zcc11*Zcc21*Zcc22**2+Zcc12*Zcc21**2*Zcc22)
      AAABR(1657) = EE**2*(Zcc11/CW**2*Zcc21**3+16D0*Zcc12
     &      /CW**2*Zcc22**3+9D0*Zcc11/SW**2*Zcc21**3-4D0
     &      *Zcc11/CW**2*Zcc21*Zcc22**2-4D0*Zcc12/CW**2*Zcc21
     &      **2*Zcc22+18D0*Mc**2/MW**2/SW**2*Zcc11*Zcc21
     &      *Zcc22**2/sb**2+18D0*Mc**2/MW**2/SW**2*Zcc12
     &      *Zcc21**2*Zcc22/sb**2)
      AAABR(1658) = EE**2*(Zt21**2/CW**2*Zuu11*Zuu21-4D0
     &      *Zt21**2/CW**2*Zuu12*Zuu22-4D0*Zt22**2/CW**2
     &      *Zuu11*Zuu21+16D0*Zt22**2/CW**2*Zuu12*Zuu22+9D0
     &      *Zt21**2/SW**2*Zuu11*Zuu21+18D0*Mq/MW**2*MtMM
     &      /SW**2*Zt21*Zt22*Zuu11*Zuu22/sb**2+18D0*Mq/MW
     &      **2*MtMM/SW**2*Zt21*Zt22*Zuu12*Zuu21/sb**2)
      AAABR(1659) = EE**2*(Zt11/CW**2*Zt21*Zuu21**2-4D0*Zt11
     &      /CW**2*Zt21*Zuu22**2-4D0*Zt12/CW**2*Zt22*Zuu21
     &      **2+16D0*Zt12/CW**2*Zt22*Zuu22**2+9D0*Zt11/SW
     &      **2*Zt21*Zuu21**2+18D0*Mq/MW**2*MtMM/SW**2*Zt11
     &      *Zt22*Zuu21*Zuu22/sb**2+18D0*Mq/MW**2*MtMM/SW
     &      **2*Zt12*Zt21*Zuu21*Zuu22/sb**2)
      AAABR(1660) = GG**2*(Zt11*Zt21**3+Zt12*Zt22**3+Zt11
     &      *Zt21*Zt22**2+Zt12*Zt21**2*Zt22)
      AAABR(1661) = EE**2*(Zt11/CW**2*Zt21**3+16D0*Zt12/CW
     &      **2*Zt22**3+9D0*Zt11/SW**2*Zt21**3-4D0*Zt11/CW
     &      **2*Zt21*Zt22**2-4D0*Zt12/CW**2*Zt21**2*Zt22
     &      +18D0*MtMM**2/MW**2/SW**2*Zt11*Zt21*Zt22**2/sb
     &      **2+18D0*MtMM**2/MW**2/SW**2*Zt12*Zt21**2*Zt22
     &      /sb**2)
      AAABR(1662) = EE**2*(Zcc11/CW**2*Zcc21*Zt21**2-4D0
     &      *Zcc11/CW**2*Zcc21*Zt22**2-4D0*Zcc12/CW**2*Zcc22
     &      *Zt21**2+16D0*Zcc12/CW**2*Zcc22*Zt22**2+9D0*Zcc11
     &      /SW**2*Zcc21*Zt21**2+18D0*Mc/MW**2*MtMM/SW**2
     &      *Zcc11*Zcc22*Zt21*Zt22/sb**2+18D0*Mc/MW**2*MtMM
     &      /SW**2*Zcc12*Zcc21*Zt21*Zt22/sb**2)
      AAABR(1663) = GG**2*(Zuu11*Zuu21**3+Zuu12*Zuu22**3
     &      +Zuu11*Zuu21*Zuu22**2+Zuu12*Zuu21**2*Zuu22)
      AAABR(1664) = EE**2*(Zuu11/CW**2*Zuu21**3+16D0*Zuu12
     &      /CW**2*Zuu22**3+9D0*Zuu11/SW**2*Zuu21**3-4D0
     &      *Zuu11/CW**2*Zuu21*Zuu22**2-4D0*Zuu12/CW**2*Zuu21
     &      **2*Zuu22+18D0*Mq**2/MW**2/SW**2*Zuu11*Zuu21
     &      *Zuu22**2/sb**2+18D0*Mq**2/MW**2/SW**2*Zuu12
     &      *Zuu21**2*Zuu22/sb**2)
      AAABR(1665) = EE**2*(Zcc11/CW**2*Zcc21*Zuu21**2-4D0
     &      *Zcc11/CW**2*Zcc21*Zuu22**2-4D0*Zcc12/CW**2*Zcc22
     &      *Zuu21**2+16D0*Zcc12/CW**2*Zcc22*Zuu22**2+9D0
     &      *Zcc11/SW**2*Zcc21*Zuu21**2+18D0*Mc/MW**2*Mq
     &      /SW**2*Zcc11*Zcc22*Zuu21*Zuu22/sb**2+18D0*Mc
     &      /MW**2*Mq/SW**2*Zcc12*Zcc21*Zuu21*Zuu22/sb**2)
      AAABR(1666) = EE**2*(Zcc21**2/CW**2*Zt11*Zt21-4D0*Zcc21
     &      **2/CW**2*Zt12*Zt22-4D0*Zcc22**2/CW**2*Zt11*Zt21
     &      +16D0*Zcc22**2/CW**2*Zt12*Zt22+9D0*Zcc21**2/SW
     &      **2*Zt11*Zt21+18D0*Mc/MW**2*MtMM/SW**2*Zcc21
     &      *Zcc22*Zt11*Zt22/sb**2+18D0*Mc/MW**2*MtMM/SW
     &      **2*Zcc21*Zcc22*Zt12*Zt21/sb**2)
      AAABR(1667) = EE**2*(Zuu11**2/CW**2-4D0*Zuu12**2/CW
     &      **2+3D0*Zuu11**2/SW**2-2D0*Zuu11**2/CW**2*sb
     &      **2+8D0*Zuu12**2/CW**2*sb**2-6D0*Mq**2/MW**2
     &      /SW**2-6D0*Zuu11**2/SW**2*sb**2)
      AAABR(1668) = EE**2*(Zcc11**2/CW**2-4D0*Zcc12**2/CW
     &      **2+3D0*Zcc11**2/SW**2-2D0*Zcc11**2/CW**2*sb
     &      **2+8D0*Zcc12**2/CW**2*sb**2-6D0*Zcc11**2/SW
     &      **2*sb**2-6D0*Mc**2/MW**2/SW**2*Zcc12**2-6D0
     &      *Mq**2/MW**2/SW**2*Zcc11**2)
      AAABR(1669) = EE**2*(Zt11**2/CW**2-4D0*Zt12**2/CW**2
     &      +3D0*Zt11**2/SW**2-2D0*Zt11**2/CW**2*sb**2+8D0
     &      *Zt12**2/CW**2*sb**2-6D0*Zt11**2/SW**2*sb**2
     &      -6D0*MbMM**2/MW**2/SW**2*Zt11**2-6D0*MtMM**2
     &      /MW**2/SW**2*Zt12**2)
      AAABR(1670) = EE**2*(Zcc11/CW**2*Zcc21-4D0*Zcc12/CW
     &      **2*Zcc22+3D0*Zcc11/SW**2*Zcc21-2D0*Zcc11/CW
     &      **2*Zcc21*sb**2+8D0*Zcc12/CW**2*Zcc22*sb**2-6D0
     &      *Zcc11/SW**2*Zcc21*sb**2-6D0*Mc**2/MW**2/SW**2
     &      *Zcc12*Zcc22-6D0*Mq**2/MW**2/SW**2*Zcc11*Zcc21)
      AAABR(1671) = EE**2*(Zt11/CW**2*Zt21-4D0*Zt12/CW**2
     &      *Zt22+3D0*Zt11/SW**2*Zt21-2D0*Zt11/CW**2*Zt21
     &      *sb**2+8D0*Zt12/CW**2*Zt22*sb**2-6D0*Zt11/SW
     &      **2*Zt21*sb**2-6D0*MbMM**2/MW**2/SW**2*Zt11*Zt21
     &      -6D0*MtMM**2/MW**2/SW**2*Zt12*Zt22)
      AAABR(1672) = EE**2*(Zuu11/CW**2*Zuu21-4D0*Zuu12/CW
     &      **2*Zuu22+3D0*Zuu11/SW**2*Zuu21-2D0*Zuu11/CW
     &      **2*Zuu21*sb**2+8D0*Zuu12/CW**2*Zuu22*sb**2-6D0
     &      *Zuu11/SW**2*Zuu21*sb**2)
      AAABR(1673) = EE**2*(Zcc11**2/CW**2-4D0*Zcc12**2/CW
     &      **2-3D0*Zcc11**2/SW**2-2D0*Zcc11**2/CW**2*sb
     &      **2+8D0*Zcc12**2/CW**2*sb**2-6D0*Mc**2/MW**2
     &      /SW**2+6D0*Zcc11**2/SW**2*sb**2)
      AAABR(1674) = EE**2*(Zt11**2/CW**2-4D0*Zt12**2/CW**2
     &      -3D0*Zt11**2/SW**2-2D0*Zt11**2/CW**2*sb**2+8D0
     &      *Zt12**2/CW**2*sb**2-6D0*MtMM**2/MW**2/SW**2
     &      +6D0*Zt11**2/SW**2*sb**2)
      AAABR(1675) = EE**2*(Zuu11**2/CW**2-4D0*Zuu12**2/CW
     &      **2-3D0*Zuu11**2/SW**2-2D0*Zuu11**2/CW**2*sb
     &      **2+8D0*Zuu12**2/CW**2*sb**2-6D0*Mq**2/MW**2
     &      /SW**2+6D0*Zuu11**2/SW**2*sb**2)
      AAABR(1676) = EE**2*(Zt11/CW**2*Zt21-4D0*Zt12/CW**2
     &      *Zt22-3D0*Zt11/SW**2*Zt21-2D0*Zt11/CW**2*Zt21
     &      *sb**2+8D0*Zt12/CW**2*Zt22*sb**2+6D0*Zt11/SW
     &      **2*Zt21*sb**2)
      AAABR(1677) = EE**2*(Zuu11/CW**2*Zuu21-4D0*Zuu12/CW
     &      **2*Zuu22-3D0*Zuu11/SW**2*Zuu21-2D0*Zuu11/CW
     &      **2*Zuu21*sb**2+8D0*Zuu12/CW**2*Zuu22*sb**2+6D0
     &      *Zuu11/SW**2*Zuu21*sb**2)
      AAABR(1678) = EE**2*(Zcc11/CW**2*Zcc21-4D0*Zcc12/CW
     &      **2*Zcc22-3D0*Zcc11/SW**2*Zcc21-2D0*Zcc11/CW
     &      **2*Zcc21*sb**2+8D0*Zcc12/CW**2*Zcc22*sb**2+6D0
     &      *Zcc11/SW**2*Zcc21*sb**2)
      AAABR(1679) = 1D0/tB*(EE**2/CW**2*Pa11*Zcc11**2*sb
     &      **2-4D0*EE**2/CW**2*Pa11*Zcc12**2*sb**2+3D0*EE
     &      **2/MW**2*Mc**2*Pa11/SW**2-3D0*EE**2*Pa11/SW
     &      **2*Zcc11**2*sb**2-3D0*EE*Lambda/MW*Mc*Pa12/SW
     &      *Sqrt2*Zcc11*Zcc12)
      AAABR(1680) = 1D0/tB*(EE**2/CW**2*Pa11*Zuu11**2*sb
     &      **2-4D0*EE**2/CW**2*Pa11*Zuu12**2*sb**2+3D0*EE
     &      **2/MW**2*Mq**2*Pa11/SW**2-3D0*EE**2*Pa11/SW
     &      **2*Zuu11**2*sb**2-3D0*EE*Lambda/MW*Mq*Pa12/SW
     &      *Sqrt2*Zuu11*Zuu12)
      AAABR(1681) = 1D0/tB*(EE**2/CW**2*Pa21*Zt11**2*sb**2
     &      -4D0*EE**2/CW**2*Pa21*Zt12**2*sb**2+3D0*EE**2
     &      /MW**2*MtMM**2*Pa21/SW**2-3D0*EE**2*Pa21/SW**2
     &      *Zt11**2*sb**2-3D0*EE*Lambda/MW*MtMM*Pa22/SW
     &      *Sqrt2*Zt11*Zt12)
      AAABR(1682) = 1D0/tB*(EE**2/CW**2*Pa11*Zt11**2*sb**2
     &      -4D0*EE**2/CW**2*Pa11*Zt12**2*sb**2+3D0*EE**2
     &      /MW**2*MtMM**2*Pa11/SW**2-3D0*EE**2*Pa11/SW**2
     &      *Zt11**2*sb**2-3D0*EE*Lambda/MW*MtMM*Pa12/SW
     &      *Sqrt2*Zt11*Zt12)
      AAABR(1683) = 1D0/tB*(EE**2/CW**2*Pa21*Zcc11**2*sb
     &      **2-4D0*EE**2/CW**2*Pa21*Zcc12**2*sb**2+3D0*EE
     &      **2/MW**2*Mc**2*Pa21/SW**2-3D0*EE**2*Pa21/SW
     &      **2*Zcc11**2*sb**2-3D0*EE*Lambda/MW*Mc*Pa22/SW
     &      *Sqrt2*Zcc11*Zcc12)
      AAABR(1684) = 1D0/tB*(EE**2/CW**2*Pa21*Zuu11**2*sb
     &      **2-4D0*EE**2/CW**2*Pa21*Zuu12**2*sb**2+3D0*EE
     &      **2/MW**2*Mq**2*Pa21/SW**2-3D0*EE**2*Pa21/SW
     &      **2*Zuu11**2*sb**2-3D0*EE*Lambda/MW*Mq*Pa22/SW
     &      *Sqrt2*Zuu11*Zuu12)
      AAABR(1685) = 1D0/tB*(2D0*EE**2/CW**2*Pa21*Zcc11*Zcc21
     &      *sb**2-8D0*EE**2/CW**2*Pa21*Zcc12*Zcc22*sb**2
     &      -6D0*EE**2*Pa21/SW**2*Zcc11*Zcc21*sb**2-3D0*EE
     &      *Lambda/MW*Mc*Pa22/SW*Sqrt2*Zcc11*Zcc22-3D0*EE
     &      *Lambda/MW*Mc*Pa22/SW*Sqrt2*Zcc12*Zcc21)
      AAABR(1686) = 1D0/tB*(2D0*EE**2/CW**2*Pa21*Zuu11*Zuu21
     &      *sb**2-8D0*EE**2/CW**2*Pa21*Zuu12*Zuu22*sb**2
     &      -6D0*EE**2*Pa21/SW**2*Zuu11*Zuu21*sb**2-3D0*EE
     &      *Lambda/MW*Mq*Pa22/SW*Sqrt2*Zuu11*Zuu22-3D0*EE
     &      *Lambda/MW*Mq*Pa22/SW*Sqrt2*Zuu12*Zuu21)
      AAABR(1687) = 1D0/tB*(2D0*EE**2/CW**2*Pa11*Zcc11*Zcc21
     &      *sb**2-8D0*EE**2/CW**2*Pa11*Zcc12*Zcc22*sb**2
     &      -6D0*EE**2*Pa11/SW**2*Zcc11*Zcc21*sb**2-3D0*EE
     &      *Lambda/MW*Mc*Pa12/SW*Sqrt2*Zcc11*Zcc22-3D0*EE
     &      *Lambda/MW*Mc*Pa12/SW*Sqrt2*Zcc12*Zcc21)
      AAABR(1688) = 1D0/tB*(2D0*EE**2/CW**2*Pa11*Zt11*Zt21
     &      *sb**2-8D0*EE**2/CW**2*Pa11*Zt12*Zt22*sb**2-6D0
     &      *EE**2*Pa11/SW**2*Zt11*Zt21*sb**2-3D0*EE*Lambda
     &      /MW*MtMM*Pa12/SW*Sqrt2*Zt11*Zt22-3D0*EE*Lambda
     &      /MW*MtMM*Pa12/SW*Sqrt2*Zt12*Zt21)
      AAABR(1689) = 1D0/tB*(2D0*EE**2/CW**2*Pa11*Zuu11*Zuu21
     &      *sb**2-8D0*EE**2/CW**2*Pa11*Zuu12*Zuu22*sb**2
     &      -6D0*EE**2*Pa11/SW**2*Zuu11*Zuu21*sb**2-3D0*EE
     &      *Lambda/MW*Mq*Pa12/SW*Sqrt2*Zuu11*Zuu22-3D0*EE
     &      *Lambda/MW*Mq*Pa12/SW*Sqrt2*Zuu12*Zuu21)
      AAABR(1690) = 1D0/tB*(2D0*EE**2/CW**2*Pa21*Zt11*Zt21
     &      *sb**2-8D0*EE**2/CW**2*Pa21*Zt12*Zt22*sb**2-6D0
     &      *EE**2*Pa21/SW**2*Zt11*Zt21*sb**2-3D0*EE*Lambda
     &      /MW*MtMM*Pa22/SW*Sqrt2*Zt11*Zt22-3D0*EE*Lambda
     &      /MW*MtMM*Pa22/SW*Sqrt2*Zt12*Zt21)
      AAABR(1691) = EE*Lambda/MW*Mq/SW*Sqrt2*Zh13/tB*(Zuu11
     &      *Zuu22-Zuu12*Zuu21)
      AAABR(1692) = EE*Lambda/MW*Mc/SW*Sqrt2*Zh23/tB*(Zcc11
     &      *Zcc22-Zcc12*Zcc21)
      AAABR(1693) = EE*Lambda/MW*MtMM/SW*Sqrt2*Zh23/tB*(Zt11
     &      *Zt22-Zt12*Zt21)
      AAABR(1694) = EE*Lambda/MW*Mc/SW*Sqrt2*Zh13/tB*(Zcc11
     &      *Zcc22-Zcc12*Zcc21)
      AAABR(1695) = EE*Lambda/MW*Mq/SW*Sqrt2*Zh23/tB*(Zuu11
     &      *Zuu22-Zuu12*Zuu21)
      AAABR(1696) = EE*Lambda/MW*Mq/SW*Sqrt2*Zh33/tB*(Zuu11
     &      *Zuu22-Zuu12*Zuu21)
      AAABR(1697) = EE*Lambda/MW*MtMM/SW*Sqrt2*Zh13/tB*(Zt11
     &      *Zt22-Zt12*Zt21)
      AAABR(1698) = EE*Lambda/MW*MtMM/SW*Sqrt2*Zh33/tB*(Zt11
     &      *Zt22-Zt12*Zt21)
      AAABR(1699) = EE*Lambda/MW*Mc/SW*Sqrt2*Zh33/tB*(Zcc11
     &      *Zcc22-Zcc12*Zcc21)
      AAABR(1700) = EE**2/CW**2*Pa11*Pa21*Zuu11**2-4D0*EE
     &      **2/CW**2*Pa11*Pa21*Zuu12**2-3D0*EE**2*Pa11*Pa21
     &      /SW**2*Zuu11**2-2D0*EE**2/CW**2*Pa11*Pa21*Zuu11
     &      **2*sb**2+8D0*EE**2/CW**2*Pa11*Pa21*Zuu12**2
     &      *sb**2-6D0*EE**2/MW**2*Mq**2*Pa11*Pa21/SW**2
     &      +6D0*EE**2*Pa11*Pa21/SW**2*Zuu11**2*sb**2+6D0
     &      *EE**2/MW**2*Mq**2*Pa11*Pa21/SW**2/sb**2+6D0
     &      *EE*Lambda/MW*Mq*Pa11*Pa22/SW*Sqrt2*Zuu11*Zuu12
     &      +6D0*EE*Lambda/MW*Mq*Pa12*Pa21/SW*Sqrt2*Zuu11
     &      *Zuu12
      AAABR(1701) = EE**2/CW**2*Pa11*Pa21*Zt11**2-4D0*EE
     &      **2/CW**2*Pa11*Pa21*Zt12**2-3D0*EE**2*Pa11*Pa21
     &      /SW**2*Zt11**2-2D0*EE**2/CW**2*Pa11*Pa21*Zt11
     &      **2*sb**2+8D0*EE**2/CW**2*Pa11*Pa21*Zt12**2*sb
     &      **2-6D0*EE**2/MW**2*MtMM**2*Pa11*Pa21/SW**2+6D0
     &      *EE**2*Pa11*Pa21/SW**2*Zt11**2*sb**2+6D0*EE**2
     &      /MW**2*MtMM**2*Pa11*Pa21/SW**2/sb**2+6D0*EE*Lambda
     &      /MW*MtMM*Pa11*Pa22/SW*Sqrt2*Zt11*Zt12+6D0*EE
     &      *Lambda/MW*MtMM*Pa12*Pa21/SW*Sqrt2*Zt11*Zt12
      AAABR(1702) = EE**2/CW**2*Pa21**2*Zuu11**2-4D0*EE**2
     &      /CW**2*Pa21**2*Zuu12**2-3D0*EE**2*Pa21**2/SW
     &      **2*Zuu11**2-2D0*EE**2/CW**2*Pa21**2*Zuu11**2
     &      *sb**2+8D0*EE**2/CW**2*Pa21**2*Zuu12**2*sb**2
     &      -6D0*EE**2/MW**2*Mq**2*Pa21**2/SW**2+6D0*EE**2
     &      *Pa21**2/SW**2*Zuu11**2*sb**2+6D0*EE**2/MW**2
     &      *Mq**2*Pa21**2/SW**2/sb**2+12D0*EE*Lambda/MW
     &      *Mq*Pa21*Pa22/SW*Sqrt2*Zuu11*Zuu12
      AAABR(1703) = EE**2/CW**2*Pa11*Pa21*Zcc11**2-4D0*EE
     &      **2/CW**2*Pa11*Pa21*Zcc12**2-3D0*EE**2*Pa11*Pa21
     &      /SW**2*Zcc11**2-2D0*EE**2/CW**2*Pa11*Pa21*Zcc11
     &      **2*sb**2+8D0*EE**2/CW**2*Pa11*Pa21*Zcc12**2
     &      *sb**2-6D0*EE**2/MW**2*Mc**2*Pa11*Pa21/SW**2
     &      +6D0*EE**2*Pa11*Pa21/SW**2*Zcc11**2*sb**2+6D0
     &      *EE**2/MW**2*Mc**2*Pa11*Pa21/SW**2/sb**2+6D0
     &      *EE*Lambda/MW*Mc*Pa11*Pa22/SW*Sqrt2*Zcc11*Zcc12
     &      +6D0*EE*Lambda/MW*Mc*Pa12*Pa21/SW*Sqrt2*Zcc11
     &      *Zcc12
      AAABR(1704) = EE**2/CW**2*Pa21**2*Zt11**2-4D0*EE**2
     &      /CW**2*Pa21**2*Zt12**2-3D0*EE**2*Pa21**2/SW**2
     &      *Zt11**2-2D0*EE**2/CW**2*Pa21**2*Zt11**2*sb**2
     &      +8D0*EE**2/CW**2*Pa21**2*Zt12**2*sb**2-6D0*EE
     &      **2/MW**2*MtMM**2*Pa21**2/SW**2+6D0*EE**2*Pa21
     &      **2/SW**2*Zt11**2*sb**2+6D0*EE**2/MW**2*MtMM
     &      **2*Pa21**2/SW**2/sb**2+12D0*EE*Lambda/MW*MtMM
     &      *Pa21*Pa22/SW*Sqrt2*Zt11*Zt12
      AAABR(1705) = EE**2/CW**2*Pa11**2*Zcc11**2-4D0*EE**2
     &      /CW**2*Pa11**2*Zcc12**2-3D0*EE**2*Pa11**2/SW
     &      **2*Zcc11**2-2D0*EE**2/CW**2*Pa11**2*Zcc11**2
     &      *sb**2+8D0*EE**2/CW**2*Pa11**2*Zcc12**2*sb**2
     &      -6D0*EE**2/MW**2*Mc**2*Pa11**2/SW**2+6D0*EE**2
     &      *Pa11**2/SW**2*Zcc11**2*sb**2+6D0*EE**2/MW**2
     &      *Mc**2*Pa11**2/SW**2/sb**2+12D0*EE*Lambda/MW
     &      *Mc*Pa11*Pa12/SW*Sqrt2*Zcc11*Zcc12
      AAABR(1706) = EE**2/CW**2*Pa11**2*Zt11**2-4D0*EE**2
     &      /CW**2*Pa11**2*Zt12**2-3D0*EE**2*Pa11**2/SW**2
     &      *Zt11**2-2D0*EE**2/CW**2*Pa11**2*Zt11**2*sb**2
     &      +8D0*EE**2/CW**2*Pa11**2*Zt12**2*sb**2-6D0*EE
     &      **2/MW**2*MtMM**2*Pa11**2/SW**2+6D0*EE**2*Pa11
     &      **2/SW**2*Zt11**2*sb**2+6D0*EE**2/MW**2*MtMM
     &      **2*Pa11**2/SW**2/sb**2+12D0*EE*Lambda/MW*MtMM
     &      *Pa11*Pa12/SW*Sqrt2*Zt11*Zt12
      AAABR(1707) = EE**2/CW**2*Pa11**2*Zuu11**2-4D0*EE**2
     &      /CW**2*Pa11**2*Zuu12**2-3D0*EE**2*Pa11**2/SW
     &      **2*Zuu11**2-2D0*EE**2/CW**2*Pa11**2*Zuu11**2
     &      *sb**2+8D0*EE**2/CW**2*Pa11**2*Zuu12**2*sb**2
     &      -6D0*EE**2/MW**2*Mq**2*Pa11**2/SW**2+6D0*EE**2
     &      *Pa11**2/SW**2*Zuu11**2*sb**2+6D0*EE**2/MW**2
     &      *Mq**2*Pa11**2/SW**2/sb**2+12D0*EE*Lambda/MW
     &      *Mq*Pa11*Pa12/SW*Sqrt2*Zuu11*Zuu12
      AAABR(1708) = EE**2/CW**2*Pa21**2*Zcc11**2-4D0*EE**2
     &      /CW**2*Pa21**2*Zcc12**2-3D0*EE**2*Pa21**2/SW
     &      **2*Zcc11**2-2D0*EE**2/CW**2*Pa21**2*Zcc11**2
     &      *sb**2+8D0*EE**2/CW**2*Pa21**2*Zcc12**2*sb**2
     &      -6D0*EE**2/MW**2*Mc**2*Pa21**2/SW**2+6D0*EE**2
     &      *Pa21**2/SW**2*Zcc11**2*sb**2+6D0*EE**2/MW**2
     &      *Mc**2*Pa21**2/SW**2/sb**2+12D0*EE*Lambda/MW
     &      *Mc*Pa21*Pa22/SW*Sqrt2*Zcc11*Zcc12
      AAABR(1709) = EE**2/CW**2*Pa11*Pa21*Zcc11*Zcc21-4D0
     &      *EE**2/CW**2*Pa11*Pa21*Zcc12*Zcc22-3D0*EE**2
     &      *Pa11*Pa21/SW**2*Zcc11*Zcc21-2D0*EE**2/CW**2
     &      *Pa11*Pa21*Zcc11*Zcc21*sb**2+8D0*EE**2/CW**2
     &      *Pa11*Pa21*Zcc12*Zcc22*sb**2+6D0*EE**2*Pa11*Pa21
     &      /SW**2*Zcc11*Zcc21*sb**2+3D0*EE*Lambda/MW*Mc
     &      *Pa11*Pa22/SW*Sqrt2*Zcc11*Zcc22+3D0*EE*Lambda
     &      /MW*Mc*Pa11*Pa22/SW*Sqrt2*Zcc12*Zcc21+3D0*EE
     &      *Lambda/MW*Mc*Pa12*Pa21/SW*Sqrt2*Zcc11*Zcc22
     &      +3D0*EE*Lambda/MW*Mc*Pa12*Pa21/SW*Sqrt2*Zcc12
     &      *Zcc21
      AAABR(1710) = EE**2/CW**2*Pa11**2*Zuu11*Zuu21-4D0*EE
     &      **2/CW**2*Pa11**2*Zuu12*Zuu22-3D0*EE**2*Pa11
     &      **2/SW**2*Zuu11*Zuu21-2D0*EE**2/CW**2*Pa11**2
     &      *Zuu11*Zuu21*sb**2+8D0*EE**2/CW**2*Pa11**2*Zuu12
     &      *Zuu22*sb**2+6D0*EE**2*Pa11**2/SW**2*Zuu11*Zuu21
     &      *sb**2+6D0*EE*Lambda/MW*Mq*Pa11*Pa12/SW*Sqrt2
     &      *Zuu11*Zuu22+6D0*EE*Lambda/MW*Mq*Pa11*Pa12/SW
     &      *Sqrt2*Zuu12*Zuu21
      AAABR(1711) = EE**2/CW**2*Pa21**2*Zuu11*Zuu21-4D0*EE
     &      **2/CW**2*Pa21**2*Zuu12*Zuu22-3D0*EE**2*Pa21
     &      **2/SW**2*Zuu11*Zuu21-2D0*EE**2/CW**2*Pa21**2
     &      *Zuu11*Zuu21*sb**2+8D0*EE**2/CW**2*Pa21**2*Zuu12
     &      *Zuu22*sb**2+6D0*EE**2*Pa21**2/SW**2*Zuu11*Zuu21
     &      *sb**2+6D0*EE*Lambda/MW*Mq*Pa21*Pa22/SW*Sqrt2
     &      *Zuu11*Zuu22+6D0*EE*Lambda/MW*Mq*Pa21*Pa22/SW
     &      *Sqrt2*Zuu12*Zuu21
      AAABR(1712) = EE**2/CW**2*Pa21**2*Zcc11*Zcc21-4D0*EE
     &      **2/CW**2*Pa21**2*Zcc12*Zcc22-3D0*EE**2*Pa21
     &      **2/SW**2*Zcc11*Zcc21-2D0*EE**2/CW**2*Pa21**2
     &      *Zcc11*Zcc21*sb**2+8D0*EE**2/CW**2*Pa21**2*Zcc12
     &      *Zcc22*sb**2+6D0*EE**2*Pa21**2/SW**2*Zcc11*Zcc21
     &      *sb**2+6D0*EE*Lambda/MW*Mc*Pa21*Pa22/SW*Sqrt2
     &      *Zcc11*Zcc22+6D0*EE*Lambda/MW*Mc*Pa21*Pa22/SW
     &      *Sqrt2*Zcc12*Zcc21
      AAABR(1713) = EE**2/CW**2*Pa11*Pa21*Zuu11*Zuu21-4D0
     &      *EE**2/CW**2*Pa11*Pa21*Zuu12*Zuu22-3D0*EE**2
     &      *Pa11*Pa21/SW**2*Zuu11*Zuu21-2D0*EE**2/CW**2
     &      *Pa11*Pa21*Zuu11*Zuu21*sb**2+8D0*EE**2/CW**2
     &      *Pa11*Pa21*Zuu12*Zuu22*sb**2+6D0*EE**2*Pa11*Pa21
     &      /SW**2*Zuu11*Zuu21*sb**2+3D0*EE*Lambda/MW*Mq
     &      *Pa11*Pa22/SW*Sqrt2*Zuu11*Zuu22+3D0*EE*Lambda
     &      /MW*Mq*Pa11*Pa22/SW*Sqrt2*Zuu12*Zuu21+3D0*EE
     &      *Lambda/MW*Mq*Pa12*Pa21/SW*Sqrt2*Zuu11*Zuu22
     &      +3D0*EE*Lambda/MW*Mq*Pa12*Pa21/SW*Sqrt2*Zuu12
     &      *Zuu21
      AAABR(1714) = EE**2/CW**2*Pa11*Pa21*Zt11*Zt21-4D0*EE
     &      **2/CW**2*Pa11*Pa21*Zt12*Zt22-3D0*EE**2*Pa11
     &      *Pa21/SW**2*Zt11*Zt21-2D0*EE**2/CW**2*Pa11*Pa21
     &      *Zt11*Zt21*sb**2+8D0*EE**2/CW**2*Pa11*Pa21*Zt12
     &      *Zt22*sb**2+6D0*EE**2*Pa11*Pa21/SW**2*Zt11*Zt21
     &      *sb**2+3D0*EE*Lambda/MW*MtMM*Pa11*Pa22/SW*Sqrt2
     &      *Zt11*Zt22+3D0*EE*Lambda/MW*MtMM*Pa11*Pa22/SW
     &      *Sqrt2*Zt12*Zt21+3D0*EE*Lambda/MW*MtMM*Pa12*Pa21
     &      /SW*Sqrt2*Zt11*Zt22+3D0*EE*Lambda/MW*MtMM*Pa12
     &      *Pa21/SW*Sqrt2*Zt12*Zt21
      AAABR(1715) = EE**2/CW**2*Pa11**2*Zt11*Zt21-4D0*EE
     &      **2/CW**2*Pa11**2*Zt12*Zt22-3D0*EE**2*Pa11**2
     &      /SW**2*Zt11*Zt21-2D0*EE**2/CW**2*Pa11**2*Zt11
     &      *Zt21*sb**2+8D0*EE**2/CW**2*Pa11**2*Zt12*Zt22
     &      *sb**2+6D0*EE**2*Pa11**2/SW**2*Zt11*Zt21*sb**2
     &      +6D0*EE*Lambda/MW*MtMM*Pa11*Pa12/SW*Sqrt2*Zt11
     &      *Zt22+6D0*EE*Lambda/MW*MtMM*Pa11*Pa12/SW*Sqrt2
     &      *Zt12*Zt21
      AAABR(1716) = EE**2/CW**2*Pa11**2*Zcc11*Zcc21-4D0*EE
     &      **2/CW**2*Pa11**2*Zcc12*Zcc22-3D0*EE**2*Pa11
     &      **2/SW**2*Zcc11*Zcc21-2D0*EE**2/CW**2*Pa11**2
     &      *Zcc11*Zcc21*sb**2+8D0*EE**2/CW**2*Pa11**2*Zcc12
     &      *Zcc22*sb**2+6D0*EE**2*Pa11**2/SW**2*Zcc11*Zcc21
     &      *sb**2+6D0*EE*Lambda/MW*Mc*Pa11*Pa12/SW*Sqrt2
     &      *Zcc11*Zcc22+6D0*EE*Lambda/MW*Mc*Pa11*Pa12/SW
     &      *Sqrt2*Zcc12*Zcc21
      AAABR(1717) = EE**2/CW**2*Pa21**2*Zt11*Zt21-4D0*EE
     &      **2/CW**2*Pa21**2*Zt12*Zt22-3D0*EE**2*Pa21**2
     &      /SW**2*Zt11*Zt21-2D0*EE**2/CW**2*Pa21**2*Zt11
     &      *Zt21*sb**2+8D0*EE**2/CW**2*Pa21**2*Zt12*Zt22
     &      *sb**2+6D0*EE**2*Pa21**2/SW**2*Zt11*Zt21*sb**2
     &      +6D0*EE*Lambda/MW*MtMM*Pa21*Pa22/SW*Sqrt2*Zt11
     &      *Zt22+6D0*EE*Lambda/MW*MtMM*Pa21*Pa22/SW*Sqrt2
     &      *Zt12*Zt21
      AAABR(1718) = EE*Lambda/MW*Mc/SW*Sqrt2*(Pa21*Zcc11
     &      *Zcc22*Zh33-Pa21*Zcc12*Zcc21*Zh33+Pa22*Zcc11
     &      *Zcc22*Zh31/sb-Pa22*Zcc12*Zcc21*Zh31/sb)
      AAABR(1719) = EE*Lambda/MW*Mc/SW*Sqrt2*(Pa11*Zcc11
     &      *Zcc22*Zh23-Pa11*Zcc12*Zcc21*Zh23+Pa12*Zcc11
     &      *Zcc22*Zh21/sb-Pa12*Zcc12*Zcc21*Zh21/sb)
      AAABR(1720) = EE*Lambda/MW*Mq/SW*Sqrt2*(Pa21*Zh23*Zuu11
     &      *Zuu22-Pa21*Zh23*Zuu12*Zuu21+Pa22*Zh21*Zuu11
     &      *Zuu22/sb-Pa22*Zh21*Zuu12*Zuu21/sb)
      AAABR(1721) = EE*Lambda/MW*Mq/SW*Sqrt2*(Pa11*Zh13*Zuu11
     &      *Zuu22-Pa11*Zh13*Zuu12*Zuu21+Pa12*Zh11*Zuu11
     &      *Zuu22/sb-Pa12*Zh11*Zuu12*Zuu21/sb)
      AAABR(1722) = EE*Lambda/MW*Mq/SW*Sqrt2*(Pa21*Zh13*Zuu11
     &      *Zuu22-Pa21*Zh13*Zuu12*Zuu21+Pa22*Zh11*Zuu11
     &      *Zuu22/sb-Pa22*Zh11*Zuu12*Zuu21/sb)
      AAABR(1723) = EE*Lambda/MW*Mq/SW*Sqrt2*(Pa11*Zh23*Zuu11
     &      *Zuu22-Pa11*Zh23*Zuu12*Zuu21+Pa12*Zh21*Zuu11
     &      *Zuu22/sb-Pa12*Zh21*Zuu12*Zuu21/sb)
      AAABR(1724) = EE*Lambda/MW*Mc/SW*Sqrt2*(Pa11*Zcc11
     &      *Zcc22*Zh13-Pa11*Zcc12*Zcc21*Zh13+Pa12*Zcc11
     &      *Zcc22*Zh11/sb-Pa12*Zcc12*Zcc21*Zh11/sb)
      AAABR(1725) = EE*Lambda/MW*Mc/SW*Sqrt2*(Pa11*Zcc11
     &      *Zcc22*Zh33-Pa11*Zcc12*Zcc21*Zh33+Pa12*Zcc11
     &      *Zcc22*Zh31/sb-Pa12*Zcc12*Zcc21*Zh31/sb)
      AAABR(1726) = EE*Lambda/MW*Mq/SW*Sqrt2*(Pa11*Zh33*Zuu11
     &      *Zuu22-Pa11*Zh33*Zuu12*Zuu21+Pa12*Zh31*Zuu11
     &      *Zuu22/sb-Pa12*Zh31*Zuu12*Zuu21/sb)
      AAABR(1727) = EE*Lambda/MW*Mc/SW*Sqrt2*(Pa21*Zcc11
     &      *Zcc22*Zh23-Pa21*Zcc12*Zcc21*Zh23+Pa22*Zcc11
     &      *Zcc22*Zh21/sb-Pa22*Zcc12*Zcc21*Zh21/sb)
      AAABR(1728) = EE*Lambda/MW*MtMM/SW*Sqrt2*(Pa21*Zh23
     &      *Zt11*Zt22-Pa21*Zh23*Zt12*Zt21+Pa22*Zh21*Zt11
     &      *Zt22/sb-Pa22*Zh21*Zt12*Zt21/sb)
      AAABR(1729) = EE*Lambda/MW*MtMM/SW*Sqrt2*(Pa21*Zh13
     &      *Zt11*Zt22-Pa21*Zh13*Zt12*Zt21+Pa22*Zh11*Zt11
     &      *Zt22/sb-Pa22*Zh11*Zt12*Zt21/sb)
      AAABR(1730) = EE*Lambda/MW*Mq/SW*Sqrt2*(Pa21*Zh33*Zuu11
     &      *Zuu22-Pa21*Zh33*Zuu12*Zuu21+Pa22*Zh31*Zuu11
     &      *Zuu22/sb-Pa22*Zh31*Zuu12*Zuu21/sb)
      AAABR(1731) = EE*Lambda/MW*Mc/SW*Sqrt2*(Pa21*Zcc11
     &      *Zcc22*Zh13-Pa21*Zcc12*Zcc21*Zh13+Pa22*Zcc11
     &      *Zcc22*Zh11/sb-Pa22*Zcc12*Zcc21*Zh11/sb)
      AAABR(1732) = EE*Lambda/MW*MtMM/SW*Sqrt2*(Pa21*Zh33
     &      *Zt11*Zt22-Pa21*Zh33*Zt12*Zt21+Pa22*Zh31*Zt11
     &      *Zt22/sb-Pa22*Zh31*Zt12*Zt21/sb)
      AAABR(1733) = EE*Lambda/MW*MtMM/SW*Sqrt2*(Pa11*Zh23
     &      *Zt11*Zt22-Pa11*Zh23*Zt12*Zt21+Pa12*Zh21*Zt11
     &      *Zt22/sb-Pa12*Zh21*Zt12*Zt21/sb)
      AAABR(1734) = EE*Lambda/MW*MtMM/SW*Sqrt2*(Pa11*Zh33
     &      *Zt11*Zt22-Pa11*Zh33*Zt12*Zt21+Pa12*Zh31*Zt11
     &      *Zt22/sb-Pa12*Zh31*Zt12*Zt21/sb)
      AAABR(1735) = EE*Lambda/MW*MtMM/SW*Sqrt2*(Pa11*Zh13
     &      *Zt11*Zt22-Pa11*Zh13*Zt12*Zt21+Pa12*Zh11*Zt11
     &      *Zt22/sb-Pa12*Zh11*Zt12*Zt21/sb)
      AAABR(1736) = EE**2/CW**2*Zh11*Zh21*Zuu11**2-4D0*EE
     &      **2/CW**2*Zh11*Zh21*Zuu12**2-EE**2/CW**2*Zh12
     &      *Zh22*Zuu11**2+4D0*EE**2/CW**2*Zh12*Zh22*Zuu12
     &      **2-3D0*EE**2/SW**2*Zh11*Zh21*Zuu11**2+3D0*EE
     &      **2/SW**2*Zh12*Zh22*Zuu11**2-6D0*EE**2/MW**2
     &      *Mq**2/SW**2*Zh12*Zh22/sb**2+6D0*EE*Lambda/MW
     &      *Mq/SW*Sqrt2*Zh11*Zh23*Zuu11*Zuu12/sb+6D0*EE
     &      *Lambda/MW*Mq/SW*Sqrt2*Zh13*Zh21*Zuu11*Zuu12/sb
      AAABR(1737) = EE**2/CW**2*Zh11*Zh21*Zt11**2-4D0*EE
     &      **2/CW**2*Zh11*Zh21*Zt12**2-EE**2/CW**2*Zh12
     &      *Zh22*Zt11**2+4D0*EE**2/CW**2*Zh12*Zh22*Zt12
     &      **2-3D0*EE**2/SW**2*Zh11*Zh21*Zt11**2+3D0*EE
     &      **2/SW**2*Zh12*Zh22*Zt11**2-6D0*EE**2/MW**2*MtMM
     &      **2/SW**2*Zh12*Zh22/sb**2+6D0*EE*Lambda/MW*MtMM
     &      /SW*Sqrt2*Zh11*Zh23*Zt11*Zt12/sb+6D0*EE*Lambda
     &      /MW*MtMM/SW*Sqrt2*Zh13*Zh21*Zt11*Zt12/sb
      AAABR(1738) = EE**2/CW**2*Zh11*Zh31*Zuu11**2-4D0*EE
     &      **2/CW**2*Zh11*Zh31*Zuu12**2-EE**2/CW**2*Zh12
     &      *Zh32*Zuu11**2+4D0*EE**2/CW**2*Zh12*Zh32*Zuu12
     &      **2-3D0*EE**2/SW**2*Zh11*Zh31*Zuu11**2+3D0*EE
     &      **2/SW**2*Zh12*Zh32*Zuu11**2-6D0*EE**2/MW**2
     &      *Mq**2/SW**2*Zh12*Zh32/sb**2+6D0*EE*Lambda/MW
     &      *Mq/SW*Sqrt2*Zh11*Zh33*Zuu11*Zuu12/sb+6D0*EE
     &      *Lambda/MW*Mq/SW*Sqrt2*Zh13*Zh31*Zuu11*Zuu12/sb
      AAABR(1739) = EE**2/CW**2*Zh21**2*Zuu11**2-4D0*EE**2
     &      /CW**2*Zh21**2*Zuu12**2-EE**2/CW**2*Zh22**2*Zuu11
     &      **2+4D0*EE**2/CW**2*Zh22**2*Zuu12**2-3D0*EE**2
     &      /SW**2*Zh21**2*Zuu11**2+3D0*EE**2/SW**2*Zh22
     &      **2*Zuu11**2-6D0*EE**2/MW**2*Mq**2/SW**2*Zh22
     &      **2/sb**2+12D0*EE*Lambda/MW*Mq/SW*Sqrt2*Zh21
     &      *Zh23*Zuu11*Zuu12/sb
      AAABR(1740) = EE**2/CW**2*Zcc11**2*Zh31**2-EE**2/CW
     &      **2*Zcc11**2*Zh32**2-4D0*EE**2/CW**2*Zcc12**2
     &      *Zh31**2+4D0*EE**2/CW**2*Zcc12**2*Zh32**2-3D0
     &      *EE**2/SW**2*Zcc11**2*Zh31**2+3D0*EE**2/SW**2
     &      *Zcc11**2*Zh32**2-6D0*EE**2/MW**2*Mc**2/SW**2
     &      *Zh32**2/sb**2+12D0*EE*Lambda/MW*Mc/SW*Sqrt2
     &      *Zcc11*Zcc12*Zh31*Zh33/sb
      AAABR(1741) = EE**2/CW**2*Zh31**2*Zt11**2-4D0*EE**2
     &      /CW**2*Zh31**2*Zt12**2-EE**2/CW**2*Zh32**2*Zt11
     &      **2+4D0*EE**2/CW**2*Zh32**2*Zt12**2-3D0*EE**2
     &      /SW**2*Zh31**2*Zt11**2+3D0*EE**2/SW**2*Zh32**2
     &      *Zt11**2-6D0*EE**2/MW**2*MtMM**2/SW**2*Zh32**2
     &      /sb**2+12D0*EE*Lambda/MW*MtMM/SW*Sqrt2*Zh31*Zh33
     &      *Zt11*Zt12/sb
      AAABR(1742) = EE**2/CW**2*Zcc11**2*Zh11*Zh21-EE**2
     &      /CW**2*Zcc11**2*Zh12*Zh22-4D0*EE**2/CW**2*Zcc12
     &      **2*Zh11*Zh21+4D0*EE**2/CW**2*Zcc12**2*Zh12*Zh22
     &      -3D0*EE**2/SW**2*Zcc11**2*Zh11*Zh21+3D0*EE**2
     &      /SW**2*Zcc11**2*Zh12*Zh22-6D0*EE**2/MW**2*Mc
     &      **2/SW**2*Zh12*Zh22/sb**2+6D0*EE*Lambda/MW*Mc
     &      /SW*Sqrt2*Zcc11*Zcc12*Zh11*Zh23/sb+6D0*EE*Lambda
     &      /MW*Mc/SW*Sqrt2*Zcc11*Zcc12*Zh13*Zh21/sb
      AAABR(1743) = EE**2/CW**2*Zcc11**2*Zh11**2-EE**2/CW
     &      **2*Zcc11**2*Zh12**2-4D0*EE**2/CW**2*Zcc12**2
     &      *Zh11**2+4D0*EE**2/CW**2*Zcc12**2*Zh12**2-3D0
     &      *EE**2/SW**2*Zcc11**2*Zh11**2+3D0*EE**2/SW**2
     &      *Zcc11**2*Zh12**2-6D0*EE**2/MW**2*Mc**2/SW**2
     &      *Zh12**2/sb**2+12D0*EE*Lambda/MW*Mc/SW*Sqrt2
     &      *Zcc11*Zcc12*Zh11*Zh13/sb
      AAABR(1744) = EE**2/CW**2*Zcc11**2*Zh21**2-EE**2/CW
     &      **2*Zcc11**2*Zh22**2-4D0*EE**2/CW**2*Zcc12**2
     &      *Zh21**2+4D0*EE**2/CW**2*Zcc12**2*Zh22**2-3D0
     &      *EE**2/SW**2*Zcc11**2*Zh21**2+3D0*EE**2/SW**2
     &      *Zcc11**2*Zh22**2-6D0*EE**2/MW**2*Mc**2/SW**2
     &      *Zh22**2/sb**2+12D0*EE*Lambda/MW*Mc/SW*Sqrt2
     &      *Zcc11*Zcc12*Zh21*Zh23/sb
      AAABR(1745) = EE**2/CW**2*Zcc11**2*Zh21*Zh31-EE**2
     &      /CW**2*Zcc11**2*Zh22*Zh32-4D0*EE**2/CW**2*Zcc12
     &      **2*Zh21*Zh31+4D0*EE**2/CW**2*Zcc12**2*Zh22*Zh32
     &      -3D0*EE**2/SW**2*Zcc11**2*Zh21*Zh31+3D0*EE**2
     &      /SW**2*Zcc11**2*Zh22*Zh32-6D0*EE**2/MW**2*Mc
     &      **2/SW**2*Zh22*Zh32/sb**2+6D0*EE*Lambda/MW*Mc
     &      /SW*Sqrt2*Zcc11*Zcc12*Zh21*Zh33/sb+6D0*EE*Lambda
     &      /MW*Mc/SW*Sqrt2*Zcc11*Zcc12*Zh23*Zh31/sb
      AAABR(1746) = EE**2/CW**2*Zh21**2*Zt11**2-4D0*EE**2
     &      /CW**2*Zh21**2*Zt12**2-EE**2/CW**2*Zh22**2*Zt11
     &      **2+4D0*EE**2/CW**2*Zh22**2*Zt12**2-3D0*EE**2
     &      /SW**2*Zh21**2*Zt11**2+3D0*EE**2/SW**2*Zh22**2
     &      *Zt11**2-6D0*EE**2/MW**2*MtMM**2/SW**2*Zh22**2
     &      /sb**2+12D0*EE*Lambda/MW*MtMM/SW*Sqrt2*Zh21*Zh23
     &      *Zt11*Zt12/sb
      AAABR(1747) = EE**2/CW**2*Zh11**2*Zt11**2-4D0*EE**2
     &      /CW**2*Zh11**2*Zt12**2-EE**2/CW**2*Zh12**2*Zt11
     &      **2+4D0*EE**2/CW**2*Zh12**2*Zt12**2-3D0*EE**2
     &      /SW**2*Zh11**2*Zt11**2+3D0*EE**2/SW**2*Zh12**2
     &      *Zt11**2-6D0*EE**2/MW**2*MtMM**2/SW**2*Zh12**2
     &      /sb**2+12D0*EE*Lambda/MW*MtMM/SW*Sqrt2*Zh11*Zh13
     &      *Zt11*Zt12/sb
      AAABR(1748) = EE**2/CW**2*Zh11*Zh31*Zt11**2-4D0*EE
     &      **2/CW**2*Zh11*Zh31*Zt12**2-EE**2/CW**2*Zh12
     &      *Zh32*Zt11**2+4D0*EE**2/CW**2*Zh12*Zh32*Zt12
     &      **2-3D0*EE**2/SW**2*Zh11*Zh31*Zt11**2+3D0*EE
     &      **2/SW**2*Zh12*Zh32*Zt11**2-6D0*EE**2/MW**2*MtMM
     &      **2/SW**2*Zh12*Zh32/sb**2+6D0*EE*Lambda/MW*MtMM
     &      /SW*Sqrt2*Zh11*Zh33*Zt11*Zt12/sb+6D0*EE*Lambda
     &      /MW*MtMM/SW*Sqrt2*Zh13*Zh31*Zt11*Zt12/sb
      AAABR(1749) = EE**2/CW**2*Zcc11**2*Zh11*Zh31-EE**2
     &      /CW**2*Zcc11**2*Zh12*Zh32-4D0*EE**2/CW**2*Zcc12
     &      **2*Zh11*Zh31+4D0*EE**2/CW**2*Zcc12**2*Zh12*Zh32
     &      -3D0*EE**2/SW**2*Zcc11**2*Zh11*Zh31+3D0*EE**2
     &      /SW**2*Zcc11**2*Zh12*Zh32-6D0*EE**2/MW**2*Mc
     &      **2/SW**2*Zh12*Zh32/sb**2+6D0*EE*Lambda/MW*Mc
     &      /SW*Sqrt2*Zcc11*Zcc12*Zh11*Zh33/sb+6D0*EE*Lambda
     &      /MW*Mc/SW*Sqrt2*Zcc11*Zcc12*Zh13*Zh31/sb
      AAABR(1750) = EE**2/CW**2*Zh21*Zh31*Zuu11**2-4D0*EE
     &      **2/CW**2*Zh21*Zh31*Zuu12**2-EE**2/CW**2*Zh22
     &      *Zh32*Zuu11**2+4D0*EE**2/CW**2*Zh22*Zh32*Zuu12
     &      **2-3D0*EE**2/SW**2*Zh21*Zh31*Zuu11**2+3D0*EE
     &      **2/SW**2*Zh22*Zh32*Zuu11**2-6D0*EE**2/MW**2
     &      *Mq**2/SW**2*Zh22*Zh32/sb**2+6D0*EE*Lambda/MW
     &      *Mq/SW*Sqrt2*Zh21*Zh33*Zuu11*Zuu12/sb+6D0*EE
     &      *Lambda/MW*Mq/SW*Sqrt2*Zh23*Zh31*Zuu11*Zuu12/sb
      AAABR(1751) = EE**2/CW**2*Zh31**2*Zuu11**2-4D0*EE**2
     &      /CW**2*Zh31**2*Zuu12**2-EE**2/CW**2*Zh32**2*Zuu11
     &      **2+4D0*EE**2/CW**2*Zh32**2*Zuu12**2-3D0*EE**2
     &      /SW**2*Zh31**2*Zuu11**2+3D0*EE**2/SW**2*Zh32
     &      **2*Zuu11**2-6D0*EE**2/MW**2*Mq**2/SW**2*Zh32
     &      **2/sb**2+12D0*EE*Lambda/MW*Mq/SW*Sqrt2*Zh31
     &      *Zh33*Zuu11*Zuu12/sb
      AAABR(1752) = EE**2/CW**2*Zh21*Zh31*Zt11**2-4D0*EE
     &      **2/CW**2*Zh21*Zh31*Zt12**2-EE**2/CW**2*Zh22
     &      *Zh32*Zt11**2+4D0*EE**2/CW**2*Zh22*Zh32*Zt12
     &      **2-3D0*EE**2/SW**2*Zh21*Zh31*Zt11**2+3D0*EE
     &      **2/SW**2*Zh22*Zh32*Zt11**2-6D0*EE**2/MW**2*MtMM
     &      **2/SW**2*Zh22*Zh32/sb**2+6D0*EE*Lambda/MW*MtMM
     &      /SW*Sqrt2*Zh21*Zh33*Zt11*Zt12/sb+6D0*EE*Lambda
     &      /MW*MtMM/SW*Sqrt2*Zh23*Zh31*Zt11*Zt12/sb
      AAABR(1753) = EE**2/CW**2*Zh11**2*Zuu11**2-4D0*EE**2
     &      /CW**2*Zh11**2*Zuu12**2-EE**2/CW**2*Zh12**2*Zuu11
     &      **2+4D0*EE**2/CW**2*Zh12**2*Zuu12**2-3D0*EE**2
     &      /SW**2*Zh11**2*Zuu11**2+3D0*EE**2/SW**2*Zh12
     &      **2*Zuu11**2-6D0*EE**2/MW**2*Mq**2/SW**2*Zh12
     &      **2/sb**2+12D0*EE*Lambda/MW*Mq/SW*Sqrt2*Zh11
     &      *Zh13*Zuu11*Zuu12/sb
      AAABR(1754) = EE**2/CW**2*Zh11*Zh21*Zt11*Zt21-4D0*EE
     &      **2/CW**2*Zh11*Zh21*Zt12*Zt22-EE**2/CW**2*Zh12
     &      *Zh22*Zt11*Zt21+4D0*EE**2/CW**2*Zh12*Zh22*Zt12
     &      *Zt22-3D0*EE**2/SW**2*Zh11*Zh21*Zt11*Zt21+3D0
     &      *EE**2/SW**2*Zh12*Zh22*Zt11*Zt21+3D0*EE*Lambda
     &      /MW*MtMM/SW*Sqrt2*Zh11*Zh23*Zt11*Zt22/sb+3D0
     &      *EE*Lambda/MW*MtMM/SW*Sqrt2*Zh11*Zh23*Zt12*Zt21
     &      /sb+3D0*EE*Lambda/MW*MtMM/SW*Sqrt2*Zh13*Zh21
     &      *Zt11*Zt22/sb+3D0*EE*Lambda/MW*MtMM/SW*Sqrt2
     &      *Zh13*Zh21*Zt12*Zt21/sb
      AAABR(1755) = EE**2/CW**2*Zh11*Zh31*Zuu11*Zuu21-4D0
     &      *EE**2/CW**2*Zh11*Zh31*Zuu12*Zuu22-EE**2/CW**2
     &      *Zh12*Zh32*Zuu11*Zuu21+4D0*EE**2/CW**2*Zh12*Zh32
     &      *Zuu12*Zuu22-3D0*EE**2/SW**2*Zh11*Zh31*Zuu11
     &      *Zuu21+3D0*EE**2/SW**2*Zh12*Zh32*Zuu11*Zuu21
     &      +3D0*EE*Lambda/MW*Mq/SW*Sqrt2*Zh11*Zh33*Zuu11
     &      *Zuu22/sb+3D0*EE*Lambda/MW*Mq/SW*Sqrt2*Zh11*Zh33
     &      *Zuu12*Zuu21/sb+3D0*EE*Lambda/MW*Mq/SW*Sqrt2
     &      *Zh13*Zh31*Zuu11*Zuu22/sb+3D0*EE*Lambda/MW*Mq
     &      /SW*Sqrt2*Zh13*Zh31*Zuu12*Zuu21/sb
      AAABR(1756) = EE**2/CW**2*Zcc11*Zcc21*Zh11*Zh31-EE
     &      **2/CW**2*Zcc11*Zcc21*Zh12*Zh32-4D0*EE**2/CW
     &      **2*Zcc12*Zcc22*Zh11*Zh31+4D0*EE**2/CW**2*Zcc12
     &      *Zcc22*Zh12*Zh32-3D0*EE**2/SW**2*Zcc11*Zcc21
     &      *Zh11*Zh31+3D0*EE**2/SW**2*Zcc11*Zcc21*Zh12*Zh32
     &      +3D0*EE*Lambda/MW*Mc/SW*Sqrt2*Zcc11*Zcc22*Zh11
     &      *Zh33/sb+3D0*EE*Lambda/MW*Mc/SW*Sqrt2*Zcc11*Zcc22
     &      *Zh13*Zh31/sb+3D0*EE*Lambda/MW*Mc/SW*Sqrt2*Zcc12
     &      *Zcc21*Zh11*Zh33/sb+3D0*EE*Lambda/MW*Mc/SW*Sqrt2
     &      *Zcc12*Zcc21*Zh13*Zh31/sb
      AAABR(1757) = EE**2/CW**2*Zcc11*Zcc21*Zh21**2-EE**2
     &      /CW**2*Zcc11*Zcc21*Zh22**2-4D0*EE**2/CW**2*Zcc12
     &      *Zcc22*Zh21**2+4D0*EE**2/CW**2*Zcc12*Zcc22*Zh22
     &      **2-3D0*EE**2/SW**2*Zcc11*Zcc21*Zh21**2+3D0*EE
     &      **2/SW**2*Zcc11*Zcc21*Zh22**2+6D0*EE*Lambda/MW
     &      *Mc/SW*Sqrt2*Zcc11*Zcc22*Zh21*Zh23/sb+6D0*EE
     &      *Lambda/MW*Mc/SW*Sqrt2*Zcc12*Zcc21*Zh21*Zh23/sb
      AAABR(1758) = EE**2/CW**2*Zcc11*Zcc21*Zh11**2-EE**2
     &      /CW**2*Zcc11*Zcc21*Zh12**2-4D0*EE**2/CW**2*Zcc12
     &      *Zcc22*Zh11**2+4D0*EE**2/CW**2*Zcc12*Zcc22*Zh12
     &      **2-3D0*EE**2/SW**2*Zcc11*Zcc21*Zh11**2+3D0*EE
     &      **2/SW**2*Zcc11*Zcc21*Zh12**2+6D0*EE*Lambda/MW
     &      *Mc/SW*Sqrt2*Zcc11*Zcc22*Zh11*Zh13/sb+6D0*EE
     &      *Lambda/MW*Mc/SW*Sqrt2*Zcc12*Zcc21*Zh11*Zh13/sb
      AAABR(1759) = EE**2/CW**2*Zh21*Zh31*Zt11*Zt21-4D0*EE
     &      **2/CW**2*Zh21*Zh31*Zt12*Zt22-EE**2/CW**2*Zh22
     &      *Zh32*Zt11*Zt21+4D0*EE**2/CW**2*Zh22*Zh32*Zt12
     &      *Zt22-3D0*EE**2/SW**2*Zh21*Zh31*Zt11*Zt21+3D0
     &      *EE**2/SW**2*Zh22*Zh32*Zt11*Zt21+3D0*EE*Lambda
     &      /MW*MtMM/SW*Sqrt2*Zh21*Zh33*Zt11*Zt22/sb+3D0
     &      *EE*Lambda/MW*MtMM/SW*Sqrt2*Zh21*Zh33*Zt12*Zt21
     &      /sb+3D0*EE*Lambda/MW*MtMM/SW*Sqrt2*Zh23*Zh31
     &      *Zt11*Zt22/sb+3D0*EE*Lambda/MW*MtMM/SW*Sqrt2
     &      *Zh23*Zh31*Zt12*Zt21/sb
      AAABR(1760) = EE**2/CW**2*Zh11*Zh21*Zuu11*Zuu21-4D0
     &      *EE**2/CW**2*Zh11*Zh21*Zuu12*Zuu22-EE**2/CW**2
     &      *Zh12*Zh22*Zuu11*Zuu21+4D0*EE**2/CW**2*Zh12*Zh22
     &      *Zuu12*Zuu22-3D0*EE**2/SW**2*Zh11*Zh21*Zuu11
     &      *Zuu21+3D0*EE**2/SW**2*Zh12*Zh22*Zuu11*Zuu21
     &      +3D0*EE*Lambda/MW*Mq/SW*Sqrt2*Zh11*Zh23*Zuu11
     &      *Zuu22/sb+3D0*EE*Lambda/MW*Mq/SW*Sqrt2*Zh11*Zh23
     &      *Zuu12*Zuu21/sb+3D0*EE*Lambda/MW*Mq/SW*Sqrt2
     &      *Zh13*Zh21*Zuu11*Zuu22/sb+3D0*EE*Lambda/MW*Mq
     &      /SW*Sqrt2*Zh13*Zh21*Zuu12*Zuu21/sb
      AAABR(1761) = EE**2/CW**2*Zh31**2*Zt11*Zt21-4D0*EE
     &      **2/CW**2*Zh31**2*Zt12*Zt22-EE**2/CW**2*Zh32
     &      **2*Zt11*Zt21+4D0*EE**2/CW**2*Zh32**2*Zt12*Zt22
     &      -3D0*EE**2/SW**2*Zh31**2*Zt11*Zt21+3D0*EE**2
     &      /SW**2*Zh32**2*Zt11*Zt21+6D0*EE*Lambda/MW*MtMM
     &      /SW*Sqrt2*Zh31*Zh33*Zt11*Zt22/sb+6D0*EE*Lambda
     &      /MW*MtMM/SW*Sqrt2*Zh31*Zh33*Zt12*Zt21/sb
      AAABR(1762) = EE**2/CW**2*Zh11*Zh31*Zt11*Zt21-4D0*EE
     &      **2/CW**2*Zh11*Zh31*Zt12*Zt22-EE**2/CW**2*Zh12
     &      *Zh32*Zt11*Zt21+4D0*EE**2/CW**2*Zh12*Zh32*Zt12
     &      *Zt22-3D0*EE**2/SW**2*Zh11*Zh31*Zt11*Zt21+3D0
     &      *EE**2/SW**2*Zh12*Zh32*Zt11*Zt21+3D0*EE*Lambda
     &      /MW*MtMM/SW*Sqrt2*Zh11*Zh33*Zt11*Zt22/sb+3D0
     &      *EE*Lambda/MW*MtMM/SW*Sqrt2*Zh11*Zh33*Zt12*Zt21
     &      /sb+3D0*EE*Lambda/MW*MtMM/SW*Sqrt2*Zh13*Zh31
     &      *Zt11*Zt22/sb+3D0*EE*Lambda/MW*MtMM/SW*Sqrt2
     &      *Zh13*Zh31*Zt12*Zt21/sb
      AAABR(1763) = EE**2/CW**2*Zcc11*Zcc21*Zh31**2-EE**2
     &      /CW**2*Zcc11*Zcc21*Zh32**2-4D0*EE**2/CW**2*Zcc12
     &      *Zcc22*Zh31**2+4D0*EE**2/CW**2*Zcc12*Zcc22*Zh32
     &      **2-3D0*EE**2/SW**2*Zcc11*Zcc21*Zh31**2+3D0*EE
     &      **2/SW**2*Zcc11*Zcc21*Zh32**2+6D0*EE*Lambda/MW
     &      *Mc/SW*Sqrt2*Zcc11*Zcc22*Zh31*Zh33/sb+6D0*EE
     &      *Lambda/MW*Mc/SW*Sqrt2*Zcc12*Zcc21*Zh31*Zh33/sb
      AAABR(1764) = EE**2/CW**2*Zh11**2*Zuu11*Zuu21-4D0*EE
     &      **2/CW**2*Zh11**2*Zuu12*Zuu22-EE**2/CW**2*Zh12
     &      **2*Zuu11*Zuu21+4D0*EE**2/CW**2*Zh12**2*Zuu12
     &      *Zuu22-3D0*EE**2/SW**2*Zh11**2*Zuu11*Zuu21+3D0
     &      *EE**2/SW**2*Zh12**2*Zuu11*Zuu21+6D0*EE*Lambda
     &      /MW*Mq/SW*Sqrt2*Zh11*Zh13*Zuu11*Zuu22/sb+6D0
     &      *EE*Lambda/MW*Mq/SW*Sqrt2*Zh11*Zh13*Zuu12*Zuu21
     &      /sb
      AAABR(1765) = EE**2/CW**2*Zh21**2*Zuu11*Zuu21-4D0*EE
     &      **2/CW**2*Zh21**2*Zuu12*Zuu22-EE**2/CW**2*Zh22
     &      **2*Zuu11*Zuu21+4D0*EE**2/CW**2*Zh22**2*Zuu12
     &      *Zuu22-3D0*EE**2/SW**2*Zh21**2*Zuu11*Zuu21+3D0
     &      *EE**2/SW**2*Zh22**2*Zuu11*Zuu21+6D0*EE*Lambda
     &      /MW*Mq/SW*Sqrt2*Zh21*Zh23*Zuu11*Zuu22/sb+6D0
     &      *EE*Lambda/MW*Mq/SW*Sqrt2*Zh21*Zh23*Zuu12*Zuu21
     &      /sb
      AAABR(1766) = EE**2/CW**2*Zh21**2*Zt11*Zt21-4D0*EE
     &      **2/CW**2*Zh21**2*Zt12*Zt22-EE**2/CW**2*Zh22
     &      **2*Zt11*Zt21+4D0*EE**2/CW**2*Zh22**2*Zt12*Zt22
     &      -3D0*EE**2/SW**2*Zh21**2*Zt11*Zt21+3D0*EE**2
     &      /SW**2*Zh22**2*Zt11*Zt21+6D0*EE*Lambda/MW*MtMM
     &      /SW*Sqrt2*Zh21*Zh23*Zt11*Zt22/sb+6D0*EE*Lambda
     &      /MW*MtMM/SW*Sqrt2*Zh21*Zh23*Zt12*Zt21/sb
      AAABR(1767) = EE**2/CW**2*Zcc11*Zcc21*Zh11*Zh21-EE
     &      **2/CW**2*Zcc11*Zcc21*Zh12*Zh22-4D0*EE**2/CW
     &      **2*Zcc12*Zcc22*Zh11*Zh21+4D0*EE**2/CW**2*Zcc12
     &      *Zcc22*Zh12*Zh22-3D0*EE**2/SW**2*Zcc11*Zcc21
     &      *Zh11*Zh21+3D0*EE**2/SW**2*Zcc11*Zcc21*Zh12*Zh22
     &      +3D0*EE*Lambda/MW*Mc/SW*Sqrt2*Zcc11*Zcc22*Zh11
     &      *Zh23/sb+3D0*EE*Lambda/MW*Mc/SW*Sqrt2*Zcc11*Zcc22
     &      *Zh13*Zh21/sb+3D0*EE*Lambda/MW*Mc/SW*Sqrt2*Zcc12
     &      *Zcc21*Zh11*Zh23/sb+3D0*EE*Lambda/MW*Mc/SW*Sqrt2
     &      *Zcc12*Zcc21*Zh13*Zh21/sb
      AAABR(1768) = EE**2/CW**2*Zh11**2*Zt11*Zt21-4D0*EE
     &      **2/CW**2*Zh11**2*Zt12*Zt22-EE**2/CW**2*Zh12
     &      **2*Zt11*Zt21+4D0*EE**2/CW**2*Zh12**2*Zt12*Zt22
     &      -3D0*EE**2/SW**2*Zh11**2*Zt11*Zt21+3D0*EE**2
     &      /SW**2*Zh12**2*Zt11*Zt21+6D0*EE*Lambda/MW*MtMM
     &      /SW*Sqrt2*Zh11*Zh13*Zt11*Zt22/sb+6D0*EE*Lambda
     &      /MW*MtMM/SW*Sqrt2*Zh11*Zh13*Zt12*Zt21/sb
      AAABR(1769) = EE**2/CW**2*Zh21*Zh31*Zuu11*Zuu21-4D0
     &      *EE**2/CW**2*Zh21*Zh31*Zuu12*Zuu22-EE**2/CW**2
     &      *Zh22*Zh32*Zuu11*Zuu21+4D0*EE**2/CW**2*Zh22*Zh32
     &      *Zuu12*Zuu22-3D0*EE**2/SW**2*Zh21*Zh31*Zuu11
     &      *Zuu21+3D0*EE**2/SW**2*Zh22*Zh32*Zuu11*Zuu21
     &      +3D0*EE*Lambda/MW*Mq/SW*Sqrt2*Zh21*Zh33*Zuu11
     &      *Zuu22/sb+3D0*EE*Lambda/MW*Mq/SW*Sqrt2*Zh21*Zh33
     &      *Zuu12*Zuu21/sb+3D0*EE*Lambda/MW*Mq/SW*Sqrt2
     &      *Zh23*Zh31*Zuu11*Zuu22/sb+3D0*EE*Lambda/MW*Mq
     &      /SW*Sqrt2*Zh23*Zh31*Zuu12*Zuu21/sb
      AAABR(1770) = EE**2/CW**2*Zh31**2*Zuu11*Zuu21-4D0*EE
     &      **2/CW**2*Zh31**2*Zuu12*Zuu22-EE**2/CW**2*Zh32
     &      **2*Zuu11*Zuu21+4D0*EE**2/CW**2*Zh32**2*Zuu12
     &      *Zuu22-3D0*EE**2/SW**2*Zh31**2*Zuu11*Zuu21+3D0
     &      *EE**2/SW**2*Zh32**2*Zuu11*Zuu21+6D0*EE*Lambda
     &      /MW*Mq/SW*Sqrt2*Zh31*Zh33*Zuu11*Zuu22/sb+6D0
     &      *EE*Lambda/MW*Mq/SW*Sqrt2*Zh31*Zh33*Zuu12*Zuu21
     &      /sb
      AAABR(1771) = EE**2/CW**2*Zcc11*Zcc21*Zh21*Zh31-EE
     &      **2/CW**2*Zcc11*Zcc21*Zh22*Zh32-4D0*EE**2/CW
     &      **2*Zcc12*Zcc22*Zh21*Zh31+4D0*EE**2/CW**2*Zcc12
     &      *Zcc22*Zh22*Zh32-3D0*EE**2/SW**2*Zcc11*Zcc21
     &      *Zh21*Zh31+3D0*EE**2/SW**2*Zcc11*Zcc21*Zh22*Zh32
     &      +3D0*EE*Lambda/MW*Mc/SW*Sqrt2*Zcc11*Zcc22*Zh21
     &      *Zh33/sb+3D0*EE*Lambda/MW*Mc/SW*Sqrt2*Zcc11*Zcc22
     &      *Zh23*Zh31/sb+3D0*EE*Lambda/MW*Mc/SW*Sqrt2*Zcc12
     &      *Zcc21*Zh21*Zh33/sb+3D0*EE*Lambda/MW*Mc/SW*Sqrt2
     &      *Zcc12*Zcc21*Zh23*Zh31/sb
      AAABR(1772) = EE**2*(Zt21**2/CW**2*Zuu21**2-4D0*Zt21
     &      **2/CW**2*Zuu22**2-4D0*Zt22**2/CW**2*Zuu21**2
     &      +16D0*Zt22**2/CW**2*Zuu22**2+9D0*Zt21**2/SW**2
     &      *Zuu21**2+36D0*Mq/MW**2*MtMM/SW**2*Zt21*Zt22
     &      *Zuu21*Zuu22/sb**2)
      AAABR(1773) = GG**2*(Zt21**4+Zt22**4+2D0*Zt21**2*Zt22
     &      **2)
      AAABR(1774) = EE**2*(Zt21**4/CW**2+16D0*Zt22**4/CW
     &      **2+9D0*Zt21**4/SW**2-8D0*Zt21**2/CW**2*Zt22
     &      **2+36D0*MtMM**2/MW**2/SW**2*Zt21**2*Zt22**2
     &      /sb**2)
      AAABR(1775) = EE**2*(Zcc21**2/CW**2*Zuu21**2-4D0*Zcc21
     &      **2/CW**2*Zuu22**2-4D0*Zcc22**2/CW**2*Zuu21**2
     &      +16D0*Zcc22**2/CW**2*Zuu22**2+9D0*Zcc21**2/SW
     &      **2*Zuu21**2+36D0*Mc/MW**2*Mq/SW**2*Zcc21*Zcc22
     &      *Zuu21*Zuu22/sb**2)
      AAABR(1776) = GG**2*(Zcc21**4+Zcc22**4+2D0*Zcc21**2
     &      *Zcc22**2)
      AAABR(1777) = EE**2*(Zcc21**4/CW**2+16D0*Zcc22**4/CW
     &      **2+9D0*Zcc21**4/SW**2-8D0*Zcc21**2/CW**2*Zcc22
     &      **2+36D0*Mc**2/MW**2/SW**2*Zcc21**2*Zcc22**2
     &      /sb**2)
      AAABR(1778) = GG**2*(Zuu21**4+Zuu22**4+2D0*Zuu21**2
     &      *Zuu22**2)
      AAABR(1779) = EE**2*(Zuu21**4/CW**2+16D0*Zuu22**4/CW
     &      **2+9D0*Zuu21**4/SW**2-8D0*Zuu21**2/CW**2*Zuu22
     &      **2+36D0*Mq**2/MW**2/SW**2*Zuu21**2*Zuu22**2
     &      /sb**2)
      AAABR(1780) = EE**2*(Zcc21**2/CW**2*Zt21**2-4D0*Zcc21
     &      **2/CW**2*Zt22**2-4D0*Zcc22**2/CW**2*Zt21**2
     &      +16D0*Zcc22**2/CW**2*Zt22**2+9D0*Zcc21**2/SW
     &      **2*Zt21**2+36D0*Mc/MW**2*MtMM/SW**2*Zcc21*Zcc22
     &      *Zt21*Zt22/sb**2)
      AAABR(1781) = EE**2*(Zcc21**2/CW**2-4D0*Zcc22**2/CW
     &      **2+3D0*Zcc21**2/SW**2-2D0*Zcc21**2/CW**2*sb
     &      **2+8D0*Zcc22**2/CW**2*sb**2-6D0*Zcc21**2/SW
     &      **2*sb**2-6D0*Mc**2/MW**2/SW**2*Zcc22**2-6D0
     &      *Mq**2/MW**2/SW**2*Zcc21**2)
      AAABR(1782) = EE**2*(Zuu21**2/CW**2-4D0*Zuu22**2/CW
     &      **2+3D0*Zuu21**2/SW**2-2D0*Zuu21**2/CW**2*sb
     &      **2+8D0*Zuu22**2/CW**2*sb**2-6D0*Mq**2/MW**2
     &      /SW**2-6D0*Zuu21**2/SW**2*sb**2)
      AAABR(1783) = EE**2*(Zt21**2/CW**2-4D0*Zt22**2/CW**2
     &      +3D0*Zt21**2/SW**2-2D0*Zt21**2/CW**2*sb**2+8D0
     &      *Zt22**2/CW**2*sb**2-6D0*Zt21**2/SW**2*sb**2
     &      -6D0*MbMM**2/MW**2/SW**2*Zt21**2-6D0*MtMM**2
     &      /MW**2/SW**2*Zt22**2)
      AAABR(1784) = EE**2*(Zcc21**2/CW**2-4D0*Zcc22**2/CW
     &      **2-3D0*Zcc21**2/SW**2-2D0*Zcc21**2/CW**2*sb
     &      **2+8D0*Zcc22**2/CW**2*sb**2-6D0*Mc**2/MW**2
     &      /SW**2+6D0*Zcc21**2/SW**2*sb**2)
      AAABR(1785) = EE**2*(Zt21**2/CW**2-4D0*Zt22**2/CW**2
     &      -3D0*Zt21**2/SW**2-2D0*Zt21**2/CW**2*sb**2+8D0
     &      *Zt22**2/CW**2*sb**2-6D0*MtMM**2/MW**2/SW**2
     &      +6D0*Zt21**2/SW**2*sb**2)
      AAABR(1786) = EE**2*(Zuu21**2/CW**2-4D0*Zuu22**2/CW
     &      **2-3D0*Zuu21**2/SW**2-2D0*Zuu21**2/CW**2*sb
     &      **2+8D0*Zuu22**2/CW**2*sb**2-6D0*Mq**2/MW**2
     &      /SW**2+6D0*Zuu21**2/SW**2*sb**2)
      AAABR(1787) = 1D0/tB*(EE**2/CW**2*Pa21*Zt21**2*sb**2
     &      -4D0*EE**2/CW**2*Pa21*Zt22**2*sb**2+3D0*EE**2
     &      /MW**2*MtMM**2*Pa21/SW**2-3D0*EE**2*Pa21/SW**2
     &      *Zt21**2*sb**2-3D0*EE*Lambda/MW*MtMM*Pa22/SW
     &      *Sqrt2*Zt21*Zt22)
      AAABR(1788) = 1D0/tB*(EE**2/CW**2*Pa11*Zuu21**2*sb
     &      **2-4D0*EE**2/CW**2*Pa11*Zuu22**2*sb**2+3D0*EE
     &      **2/MW**2*Mq**2*Pa11/SW**2-3D0*EE**2*Pa11/SW
     &      **2*Zuu21**2*sb**2-3D0*EE*Lambda/MW*Mq*Pa12/SW
     &      *Sqrt2*Zuu21*Zuu22)
      AAABR(1789) = 1D0/tB*(EE**2/CW**2*Pa21*Zcc21**2*sb
     &      **2-4D0*EE**2/CW**2*Pa21*Zcc22**2*sb**2+3D0*EE
     &      **2/MW**2*Mc**2*Pa21/SW**2-3D0*EE**2*Pa21/SW
     &      **2*Zcc21**2*sb**2-3D0*EE*Lambda/MW*Mc*Pa22/SW
     &      *Sqrt2*Zcc21*Zcc22)
      AAABR(1790) = 1D0/tB*(EE**2/CW**2*Pa21*Zuu21**2*sb
     &      **2-4D0*EE**2/CW**2*Pa21*Zuu22**2*sb**2+3D0*EE
     &      **2/MW**2*Mq**2*Pa21/SW**2-3D0*EE**2*Pa21/SW
     &      **2*Zuu21**2*sb**2-3D0*EE*Lambda/MW*Mq*Pa22/SW
     &      *Sqrt2*Zuu21*Zuu22)
      AAABR(1791) = 1D0/tB*(EE**2/CW**2*Pa11*Zt21**2*sb**2
     &      -4D0*EE**2/CW**2*Pa11*Zt22**2*sb**2+3D0*EE**2
     &      /MW**2*MtMM**2*Pa11/SW**2-3D0*EE**2*Pa11/SW**2
     &      *Zt21**2*sb**2-3D0*EE*Lambda/MW*MtMM*Pa12/SW
     &      *Sqrt2*Zt21*Zt22)
      AAABR(1792) = 1D0/tB*(EE**2/CW**2*Pa11*Zcc21**2*sb
     &      **2-4D0*EE**2/CW**2*Pa11*Zcc22**2*sb**2+3D0*EE
     &      **2/MW**2*Mc**2*Pa11/SW**2-3D0*EE**2*Pa11/SW
     &      **2*Zcc21**2*sb**2-3D0*EE*Lambda/MW*Mc*Pa12/SW
     &      *Sqrt2*Zcc21*Zcc22)
      AAABR(1793) = EE**2/CW**2*Pa11**2*Zcc21**2-4D0*EE**2
     &      /CW**2*Pa11**2*Zcc22**2-3D0*EE**2*Pa11**2/SW
     &      **2*Zcc21**2-2D0*EE**2/CW**2*Pa11**2*Zcc21**2
     &      *sb**2+8D0*EE**2/CW**2*Pa11**2*Zcc22**2*sb**2
     &      -6D0*EE**2/MW**2*Mc**2*Pa11**2/SW**2+6D0*EE**2
     &      *Pa11**2/SW**2*Zcc21**2*sb**2+6D0*EE**2/MW**2
     &      *Mc**2*Pa11**2/SW**2/sb**2+12D0*EE*Lambda/MW
     &      *Mc*Pa11*Pa12/SW*Sqrt2*Zcc21*Zcc22
      AAABR(1794) = EE**2/CW**2*Pa11**2*Zt21**2-4D0*EE**2
     &      /CW**2*Pa11**2*Zt22**2-3D0*EE**2*Pa11**2/SW**2
     &      *Zt21**2-2D0*EE**2/CW**2*Pa11**2*Zt21**2*sb**2
     &      +8D0*EE**2/CW**2*Pa11**2*Zt22**2*sb**2-6D0*EE
     &      **2/MW**2*MtMM**2*Pa11**2/SW**2+6D0*EE**2*Pa11
     &      **2/SW**2*Zt21**2*sb**2+6D0*EE**2/MW**2*MtMM
     &      **2*Pa11**2/SW**2/sb**2+12D0*EE*Lambda/MW*MtMM
     &      *Pa11*Pa12/SW*Sqrt2*Zt21*Zt22
      AAABR(1795) = EE**2/CW**2*Pa21**2*Zcc21**2-4D0*EE**2
     &      /CW**2*Pa21**2*Zcc22**2-3D0*EE**2*Pa21**2/SW
     &      **2*Zcc21**2-2D0*EE**2/CW**2*Pa21**2*Zcc21**2
     &      *sb**2+8D0*EE**2/CW**2*Pa21**2*Zcc22**2*sb**2
     &      -6D0*EE**2/MW**2*Mc**2*Pa21**2/SW**2+6D0*EE**2
     &      *Pa21**2/SW**2*Zcc21**2*sb**2+6D0*EE**2/MW**2
     &      *Mc**2*Pa21**2/SW**2/sb**2+12D0*EE*Lambda/MW
     &      *Mc*Pa21*Pa22/SW*Sqrt2*Zcc21*Zcc22
      AAABR(1796) = EE**2/CW**2*Pa11*Pa21*Zcc21**2-4D0*EE
     &      **2/CW**2*Pa11*Pa21*Zcc22**2-3D0*EE**2*Pa11*Pa21
     &      /SW**2*Zcc21**2-2D0*EE**2/CW**2*Pa11*Pa21*Zcc21
     &      **2*sb**2+8D0*EE**2/CW**2*Pa11*Pa21*Zcc22**2
     &      *sb**2-6D0*EE**2/MW**2*Mc**2*Pa11*Pa21/SW**2
     &      +6D0*EE**2*Pa11*Pa21/SW**2*Zcc21**2*sb**2+6D0
     &      *EE**2/MW**2*Mc**2*Pa11*Pa21/SW**2/sb**2+6D0
     &      *EE*Lambda/MW*Mc*Pa11*Pa22/SW*Sqrt2*Zcc21*Zcc22
     &      +6D0*EE*Lambda/MW*Mc*Pa12*Pa21/SW*Sqrt2*Zcc21
     &      *Zcc22
      AAABR(1797) = EE**2/CW**2*Pa11*Pa21*Zuu21**2-4D0*EE
     &      **2/CW**2*Pa11*Pa21*Zuu22**2-3D0*EE**2*Pa11*Pa21
     &      /SW**2*Zuu21**2-2D0*EE**2/CW**2*Pa11*Pa21*Zuu21
     &      **2*sb**2+8D0*EE**2/CW**2*Pa11*Pa21*Zuu22**2
     &      *sb**2-6D0*EE**2/MW**2*Mq**2*Pa11*Pa21/SW**2
     &      +6D0*EE**2*Pa11*Pa21/SW**2*Zuu21**2*sb**2+6D0
     &      *EE**2/MW**2*Mq**2*Pa11*Pa21/SW**2/sb**2+6D0
     &      *EE*Lambda/MW*Mq*Pa11*Pa22/SW*Sqrt2*Zuu21*Zuu22
     &      +6D0*EE*Lambda/MW*Mq*Pa12*Pa21/SW*Sqrt2*Zuu21
     &      *Zuu22
      AAABR(1798) = EE**2/CW**2*Pa11*Pa21*Zt21**2-4D0*EE
     &      **2/CW**2*Pa11*Pa21*Zt22**2-3D0*EE**2*Pa11*Pa21
     &      /SW**2*Zt21**2-2D0*EE**2/CW**2*Pa11*Pa21*Zt21
     &      **2*sb**2+8D0*EE**2/CW**2*Pa11*Pa21*Zt22**2*sb
     &      **2-6D0*EE**2/MW**2*MtMM**2*Pa11*Pa21/SW**2+6D0
     &      *EE**2*Pa11*Pa21/SW**2*Zt21**2*sb**2+6D0*EE**2
     &      /MW**2*MtMM**2*Pa11*Pa21/SW**2/sb**2+6D0*EE*Lambda
     &      /MW*MtMM*Pa11*Pa22/SW*Sqrt2*Zt21*Zt22+6D0*EE
     &      *Lambda/MW*MtMM*Pa12*Pa21/SW*Sqrt2*Zt21*Zt22
      AAABR(1799) = EE**2/CW**2*Pa11**2*Zuu21**2-4D0*EE**2
     &      /CW**2*Pa11**2*Zuu22**2-3D0*EE**2*Pa11**2/SW
     &      **2*Zuu21**2-2D0*EE**2/CW**2*Pa11**2*Zuu21**2
     &      *sb**2+8D0*EE**2/CW**2*Pa11**2*Zuu22**2*sb**2
     &      -6D0*EE**2/MW**2*Mq**2*Pa11**2/SW**2+6D0*EE**2
     &      *Pa11**2/SW**2*Zuu21**2*sb**2+6D0*EE**2/MW**2
     &      *Mq**2*Pa11**2/SW**2/sb**2+12D0*EE*Lambda/MW
     &      *Mq*Pa11*Pa12/SW*Sqrt2*Zuu21*Zuu22
      AAABR(1800) = EE**2/CW**2*Pa21**2*Zuu21**2-4D0*EE**2
     &      /CW**2*Pa21**2*Zuu22**2-3D0*EE**2*Pa21**2/SW
     &      **2*Zuu21**2-2D0*EE**2/CW**2*Pa21**2*Zuu21**2
     &      *sb**2+8D0*EE**2/CW**2*Pa21**2*Zuu22**2*sb**2
     &      -6D0*EE**2/MW**2*Mq**2*Pa21**2/SW**2+6D0*EE**2
     &      *Pa21**2/SW**2*Zuu21**2*sb**2+6D0*EE**2/MW**2
     &      *Mq**2*Pa21**2/SW**2/sb**2+12D0*EE*Lambda/MW
     &      *Mq*Pa21*Pa22/SW*Sqrt2*Zuu21*Zuu22
      AAABR(1801) = EE**2/CW**2*Pa21**2*Zt21**2-4D0*EE**2
     &      /CW**2*Pa21**2*Zt22**2-3D0*EE**2*Pa21**2/SW**2
     &      *Zt21**2-2D0*EE**2/CW**2*Pa21**2*Zt21**2*sb**2
     &      +8D0*EE**2/CW**2*Pa21**2*Zt22**2*sb**2-6D0*EE
     &      **2/MW**2*MtMM**2*Pa21**2/SW**2+6D0*EE**2*Pa21
     &      **2/SW**2*Zt21**2*sb**2+6D0*EE**2/MW**2*MtMM
     &      **2*Pa21**2/SW**2/sb**2+12D0*EE*Lambda/MW*MtMM
     &      *Pa21*Pa22/SW*Sqrt2*Zt21*Zt22
      AAABR(1802) = EE**2/CW**2*Zcc21**2*Zh31**2-EE**2/CW
     &      **2*Zcc21**2*Zh32**2-4D0*EE**2/CW**2*Zcc22**2
     &      *Zh31**2+4D0*EE**2/CW**2*Zcc22**2*Zh32**2-3D0
     &      *EE**2/SW**2*Zcc21**2*Zh31**2+3D0*EE**2/SW**2
     &      *Zcc21**2*Zh32**2-6D0*EE**2/MW**2*Mc**2/SW**2
     &      *Zh32**2/sb**2+12D0*EE*Lambda/MW*Mc/SW*Sqrt2
     &      *Zcc21*Zcc22*Zh31*Zh33/sb
      AAABR(1803) = EE**2/CW**2*Zh21*Zh31*Zt21**2-4D0*EE
     &      **2/CW**2*Zh21*Zh31*Zt22**2-EE**2/CW**2*Zh22
     &      *Zh32*Zt21**2+4D0*EE**2/CW**2*Zh22*Zh32*Zt22
     &      **2-3D0*EE**2/SW**2*Zh21*Zh31*Zt21**2+3D0*EE
     &      **2/SW**2*Zh22*Zh32*Zt21**2-6D0*EE**2/MW**2*MtMM
     &      **2/SW**2*Zh22*Zh32/sb**2+6D0*EE*Lambda/MW*MtMM
     &      /SW*Sqrt2*Zh21*Zh33*Zt21*Zt22/sb+6D0*EE*Lambda
     &      /MW*MtMM/SW*Sqrt2*Zh23*Zh31*Zt21*Zt22/sb
      AAABR(1804) = EE**2/CW**2*Zcc21**2*Zh21**2-EE**2/CW
     &      **2*Zcc21**2*Zh22**2-4D0*EE**2/CW**2*Zcc22**2
     &      *Zh21**2+4D0*EE**2/CW**2*Zcc22**2*Zh22**2-3D0
     &      *EE**2/SW**2*Zcc21**2*Zh21**2+3D0*EE**2/SW**2
     &      *Zcc21**2*Zh22**2-6D0*EE**2/MW**2*Mc**2/SW**2
     &      *Zh22**2/sb**2+12D0*EE*Lambda/MW*Mc/SW*Sqrt2
     &      *Zcc21*Zcc22*Zh21*Zh23/sb
      AAABR(1805) = EE**2/CW**2*Zh11*Zh31*Zuu21**2-4D0*EE
     &      **2/CW**2*Zh11*Zh31*Zuu22**2-EE**2/CW**2*Zh12
     &      *Zh32*Zuu21**2+4D0*EE**2/CW**2*Zh12*Zh32*Zuu22
     &      **2-3D0*EE**2/SW**2*Zh11*Zh31*Zuu21**2+3D0*EE
     &      **2/SW**2*Zh12*Zh32*Zuu21**2-6D0*EE**2/MW**2
     &      *Mq**2/SW**2*Zh12*Zh32/sb**2+6D0*EE*Lambda/MW
     &      *Mq/SW*Sqrt2*Zh11*Zh33*Zuu21*Zuu22/sb+6D0*EE
     &      *Lambda/MW*Mq/SW*Sqrt2*Zh13*Zh31*Zuu21*Zuu22/sb
      AAABR(1806) = EE**2/CW**2*Zcc21**2*Zh11*Zh21-EE**2
     &      /CW**2*Zcc21**2*Zh12*Zh22-4D0*EE**2/CW**2*Zcc22
     &      **2*Zh11*Zh21+4D0*EE**2/CW**2*Zcc22**2*Zh12*Zh22
     &      -3D0*EE**2/SW**2*Zcc21**2*Zh11*Zh21+3D0*EE**2
     &      /SW**2*Zcc21**2*Zh12*Zh22-6D0*EE**2/MW**2*Mc
     &      **2/SW**2*Zh12*Zh22/sb**2+6D0*EE*Lambda/MW*Mc
     &      /SW*Sqrt2*Zcc21*Zcc22*Zh11*Zh23/sb+6D0*EE*Lambda
     &      /MW*Mc/SW*Sqrt2*Zcc21*Zcc22*Zh13*Zh21/sb
      AAABR(1807) = EE**2/CW**2*Zh21*Zh31*Zuu21**2-4D0*EE
     &      **2/CW**2*Zh21*Zh31*Zuu22**2-EE**2/CW**2*Zh22
     &      *Zh32*Zuu21**2+4D0*EE**2/CW**2*Zh22*Zh32*Zuu22
     &      **2-3D0*EE**2/SW**2*Zh21*Zh31*Zuu21**2+3D0*EE
     &      **2/SW**2*Zh22*Zh32*Zuu21**2-6D0*EE**2/MW**2
     &      *Mq**2/SW**2*Zh22*Zh32/sb**2+6D0*EE*Lambda/MW
     &      *Mq/SW*Sqrt2*Zh21*Zh33*Zuu21*Zuu22/sb+6D0*EE
     &      *Lambda/MW*Mq/SW*Sqrt2*Zh23*Zh31*Zuu21*Zuu22/sb
      AAABR(1808) = EE**2/CW**2*Zh21**2*Zuu21**2-4D0*EE**2
     &      /CW**2*Zh21**2*Zuu22**2-EE**2/CW**2*Zh22**2*Zuu21
     &      **2+4D0*EE**2/CW**2*Zh22**2*Zuu22**2-3D0*EE**2
     &      /SW**2*Zh21**2*Zuu21**2+3D0*EE**2/SW**2*Zh22
     &      **2*Zuu21**2-6D0*EE**2/MW**2*Mq**2/SW**2*Zh22
     &      **2/sb**2+12D0*EE*Lambda/MW*Mq/SW*Sqrt2*Zh21
     &      *Zh23*Zuu21*Zuu22/sb
      AAABR(1809) = EE**2/CW**2*Zh31**2*Zt21**2-4D0*EE**2
     &      /CW**2*Zh31**2*Zt22**2-EE**2/CW**2*Zh32**2*Zt21
     &      **2+4D0*EE**2/CW**2*Zh32**2*Zt22**2-3D0*EE**2
     &      /SW**2*Zh31**2*Zt21**2+3D0*EE**2/SW**2*Zh32**2
     &      *Zt21**2-6D0*EE**2/MW**2*MtMM**2/SW**2*Zh32**2
     &      /sb**2+12D0*EE*Lambda/MW*MtMM/SW*Sqrt2*Zh31*Zh33
     &      *Zt21*Zt22/sb
      AAABR(1810) = EE**2/CW**2*Zh11*Zh31*Zt21**2-4D0*EE
     &      **2/CW**2*Zh11*Zh31*Zt22**2-EE**2/CW**2*Zh12
     &      *Zh32*Zt21**2+4D0*EE**2/CW**2*Zh12*Zh32*Zt22
     &      **2-3D0*EE**2/SW**2*Zh11*Zh31*Zt21**2+3D0*EE
     &      **2/SW**2*Zh12*Zh32*Zt21**2-6D0*EE**2/MW**2*MtMM
     &      **2/SW**2*Zh12*Zh32/sb**2+6D0*EE*Lambda/MW*MtMM
     &      /SW*Sqrt2*Zh11*Zh33*Zt21*Zt22/sb+6D0*EE*Lambda
     &      /MW*MtMM/SW*Sqrt2*Zh13*Zh31*Zt21*Zt22/sb
      AAABR(1811) = EE**2/CW**2*Zcc21**2*Zh21*Zh31-EE**2
     &      /CW**2*Zcc21**2*Zh22*Zh32-4D0*EE**2/CW**2*Zcc22
     &      **2*Zh21*Zh31+4D0*EE**2/CW**2*Zcc22**2*Zh22*Zh32
     &      -3D0*EE**2/SW**2*Zcc21**2*Zh21*Zh31+3D0*EE**2
     &      /SW**2*Zcc21**2*Zh22*Zh32-6D0*EE**2/MW**2*Mc
     &      **2/SW**2*Zh22*Zh32/sb**2+6D0*EE*Lambda/MW*Mc
     &      /SW*Sqrt2*Zcc21*Zcc22*Zh21*Zh33/sb+6D0*EE*Lambda
     &      /MW*Mc/SW*Sqrt2*Zcc21*Zcc22*Zh23*Zh31/sb
      AAABR(1812) = EE**2/CW**2*Zh11**2*Zt21**2-4D0*EE**2
     &      /CW**2*Zh11**2*Zt22**2-EE**2/CW**2*Zh12**2*Zt21
     &      **2+4D0*EE**2/CW**2*Zh12**2*Zt22**2-3D0*EE**2
     &      /SW**2*Zh11**2*Zt21**2+3D0*EE**2/SW**2*Zh12**2
     &      *Zt21**2-6D0*EE**2/MW**2*MtMM**2/SW**2*Zh12**2
     &      /sb**2+12D0*EE*Lambda/MW*MtMM/SW*Sqrt2*Zh11*Zh13
     &      *Zt21*Zt22/sb
      AAABR(1813) = EE**2/CW**2*Zh21**2*Zt21**2-4D0*EE**2
     &      /CW**2*Zh21**2*Zt22**2-EE**2/CW**2*Zh22**2*Zt21
     &      **2+4D0*EE**2/CW**2*Zh22**2*Zt22**2-3D0*EE**2
     &      /SW**2*Zh21**2*Zt21**2+3D0*EE**2/SW**2*Zh22**2
     &      *Zt21**2-6D0*EE**2/MW**2*MtMM**2/SW**2*Zh22**2
     &      /sb**2+12D0*EE*Lambda/MW*MtMM/SW*Sqrt2*Zh21*Zh23
     &      *Zt21*Zt22/sb
      AAABR(1814) = EE**2/CW**2*Zcc21**2*Zh11**2-EE**2/CW
     &      **2*Zcc21**2*Zh12**2-4D0*EE**2/CW**2*Zcc22**2
     &      *Zh11**2+4D0*EE**2/CW**2*Zcc22**2*Zh12**2-3D0
     &      *EE**2/SW**2*Zcc21**2*Zh11**2+3D0*EE**2/SW**2
     &      *Zcc21**2*Zh12**2-6D0*EE**2/MW**2*Mc**2/SW**2
     &      *Zh12**2/sb**2+12D0*EE*Lambda/MW*Mc/SW*Sqrt2
     &      *Zcc21*Zcc22*Zh11*Zh13/sb
      AAABR(1815) = EE**2/CW**2*Zh11**2*Zuu21**2-4D0*EE**2
     &      /CW**2*Zh11**2*Zuu22**2-EE**2/CW**2*Zh12**2*Zuu21
     &      **2+4D0*EE**2/CW**2*Zh12**2*Zuu22**2-3D0*EE**2
     &      /SW**2*Zh11**2*Zuu21**2+3D0*EE**2/SW**2*Zh12
     &      **2*Zuu21**2-6D0*EE**2/MW**2*Mq**2/SW**2*Zh12
     &      **2/sb**2+12D0*EE*Lambda/MW*Mq/SW*Sqrt2*Zh11
     &      *Zh13*Zuu21*Zuu22/sb
      AAABR(1816) = EE**2/CW**2*Zh31**2*Zuu21**2-4D0*EE**2
     &      /CW**2*Zh31**2*Zuu22**2-EE**2/CW**2*Zh32**2*Zuu21
     &      **2+4D0*EE**2/CW**2*Zh32**2*Zuu22**2-3D0*EE**2
     &      /SW**2*Zh31**2*Zuu21**2+3D0*EE**2/SW**2*Zh32
     &      **2*Zuu21**2-6D0*EE**2/MW**2*Mq**2/SW**2*Zh32
     &      **2/sb**2+12D0*EE*Lambda/MW*Mq/SW*Sqrt2*Zh31
     &      *Zh33*Zuu21*Zuu22/sb
      AAABR(1817) = EE**2/CW**2*Zh11*Zh21*Zuu21**2-4D0*EE
     &      **2/CW**2*Zh11*Zh21*Zuu22**2-EE**2/CW**2*Zh12
     &      *Zh22*Zuu21**2+4D0*EE**2/CW**2*Zh12*Zh22*Zuu22
     &      **2-3D0*EE**2/SW**2*Zh11*Zh21*Zuu21**2+3D0*EE
     &      **2/SW**2*Zh12*Zh22*Zuu21**2-6D0*EE**2/MW**2
     &      *Mq**2/SW**2*Zh12*Zh22/sb**2+6D0*EE*Lambda/MW
     &      *Mq/SW*Sqrt2*Zh11*Zh23*Zuu21*Zuu22/sb+6D0*EE
     &      *Lambda/MW*Mq/SW*Sqrt2*Zh13*Zh21*Zuu21*Zuu22/sb
      AAABR(1818) = EE**2/CW**2*Zcc21**2*Zh11*Zh31-EE**2
     &      /CW**2*Zcc21**2*Zh12*Zh32-4D0*EE**2/CW**2*Zcc22
     &      **2*Zh11*Zh31+4D0*EE**2/CW**2*Zcc22**2*Zh12*Zh32
     &      -3D0*EE**2/SW**2*Zcc21**2*Zh11*Zh31+3D0*EE**2
     &      /SW**2*Zcc21**2*Zh12*Zh32-6D0*EE**2/MW**2*Mc
     &      **2/SW**2*Zh12*Zh32/sb**2+6D0*EE*Lambda/MW*Mc
     &      /SW*Sqrt2*Zcc21*Zcc22*Zh11*Zh33/sb+6D0*EE*Lambda
     &      /MW*Mc/SW*Sqrt2*Zcc21*Zcc22*Zh13*Zh31/sb
      AAABR(1819) = EE**2/CW**2*Zh11*Zh21*Zt21**2-4D0*EE
     &      **2/CW**2*Zh11*Zh21*Zt22**2-EE**2/CW**2*Zh12
     &      *Zh22*Zt21**2+4D0*EE**2/CW**2*Zh12*Zh22*Zt22
     &      **2-3D0*EE**2/SW**2*Zh11*Zh21*Zt21**2+3D0*EE
     &      **2/SW**2*Zh12*Zh22*Zt21**2-6D0*EE**2/MW**2*MtMM
     &      **2/SW**2*Zh12*Zh22/sb**2+6D0*EE*Lambda/MW*MtMM
     &      /SW*Sqrt2*Zh11*Zh23*Zt21*Zt22/sb+6D0*EE*Lambda
     &      /MW*MtMM/SW*Sqrt2*Zh13*Zh21*Zt21*Zt22/sb
      AAABR(1820) = la1-2D0*la1*sb**2+la1*sb**4+la2*sb**4
     &      +2D0*la3*sb**2-2D0*la3*sb**4+2D0*la4*sb**2-2D0
     &      *la4*sb**4+2D0*la5*sb**2-2D0*la5*sb**4-4D0*cb
     &      *la6*sb+4D0*cb*la6*sb**3-4D0*cb*la7*sb**3
      AAABR(1821) = Pa21*(la6-5D0*la6*sb**2+4D0*la6*sb**4
     &      +3D0*la7*sb**2-4D0*la7*sb**4+cb*la1*sb-cb*la1
     &      *sb**3-cb*la2*sb**3-cb*la3*sb+2D0*cb*la3*sb**3
     &      -cb*la4*sb+2D0*cb*la4*sb**3-cb*la5*sb+2D0*cb
     &      *la5*sb**3)
      AAABR(1822) = Pa11*(la6-5D0*la6*sb**2+4D0*la6*sb**4
     &      +3D0*la7*sb**2-4D0*la7*sb**4+cb*la1*sb-cb*la1
     &      *sb**3-cb*la2*sb**3-cb*la3*sb+2D0*cb*la3*sb**3
     &      -cb*la4*sb+2D0*cb*la4*sb**3-cb*la5*sb+2D0*cb
     &      *la5*sb**3)
      AAABR(1823) = Pa21**2*la3+Pa22**2*la1s+Pa21**2*la1
     &      *sb**2-Pa21**2*la1*sb**4+Pa21**2*la2*sb**2-Pa21
     &      **2*la2*sb**4-2D0*Pa21**2*la3*sb**2+2D0*Pa21
     &      **2*la3*sb**4-2D0*Pa21**2*la4*sb**2+2D0*Pa21
     &      **2*la4*sb**4-2D0*Pa21**2*la5*sb**2+2D0*Pa21
     &      **2*la5*sb**4-Pa22**2*la1s*sb**2+Pa22**2*la2s
     &      *sb**2+2D0*Pa21**2*cb*la6*sb-4D0*Pa21**2*cb*la6
     &      *sb**3-2D0*Pa21**2*cb*la7*sb+4D0*Pa21**2*cb*la7
     &      *sb**3-2D0*Pa22**2*cb*la4s*sb+2D0*Pa22**2*cb
     &      *la5s*sb+2D0*Pa22**2*cb*la6s*sb
      AAABR(1824) = Pa11**2*la3+Pa12**2*la1s+Pa11**2*la1
     &      *sb**2-Pa11**2*la1*sb**4+Pa11**2*la2*sb**2-Pa11
     &      **2*la2*sb**4-2D0*Pa11**2*la3*sb**2+2D0*Pa11
     &      **2*la3*sb**4-2D0*Pa11**2*la4*sb**2+2D0*Pa11
     &      **2*la4*sb**4-2D0*Pa11**2*la5*sb**2+2D0*Pa11
     &      **2*la5*sb**4-Pa12**2*la1s*sb**2+Pa12**2*la2s
     &      *sb**2+2D0*Pa11**2*cb*la6*sb-4D0*Pa11**2*cb*la6
     &      *sb**3-2D0*Pa11**2*cb*la7*sb+4D0*Pa11**2*cb*la7
     &      *sb**3-2D0*Pa12**2*cb*la4s*sb+2D0*Pa12**2*cb
     &      *la5s*sb+2D0*Pa12**2*cb*la6s*sb
      AAABR(1825) = Pa11*Pa21*la3+Pa12*Pa22*la1s+Pa11*Pa21
     &      *la1*sb**2-Pa11*Pa21*la1*sb**4+Pa11*Pa21*la2
     &      *sb**2-Pa11*Pa21*la2*sb**4-2D0*Pa11*Pa21*la3
     &      *sb**2+2D0*Pa11*Pa21*la3*sb**4-2D0*Pa11*Pa21
     &      *la4*sb**2+2D0*Pa11*Pa21*la4*sb**4-2D0*Pa11*Pa21
     &      *la5*sb**2+2D0*Pa11*Pa21*la5*sb**4-Pa12*Pa22
     &      *la1s*sb**2+Pa12*Pa22*la2s*sb**2+2D0*Pa11*Pa21
     &      *cb*la6*sb-4D0*Pa11*Pa21*cb*la6*sb**3-2D0*Pa11
     &      *Pa21*cb*la7*sb+4D0*Pa11*Pa21*cb*la7*sb**3-2D0
     &      *Pa12*Pa22*cb*la4s*sb+2D0*Pa12*Pa22*cb*la5s*sb
     &      +2D0*Pa12*Pa22*cb*la6s*sb
      AAABR(1826) = Zh11*Zh31*la1-Zh11*Zh32*la6-Zh12*Zh31
     &      *la6+Zh12*Zh32*la3+Zh13*Zh33*la1s-Zh11*Zh31*la1
     &      *sb**2+Zh11*Zh31*la3*sb**2+Zh11*Zh32*la6*sb**2
     &      -Zh11*Zh32*la7*sb**2+Zh12*Zh31*la6*sb**2-Zh12
     &      *Zh31*la7*sb**2+Zh12*Zh32*la2*sb**2-Zh12*Zh32
     &      *la3*sb**2-Zh13*Zh33*la1s*sb**2+Zh13*Zh33*la2s
     &      *sb**2-2D0*Zh11*Zh31*cb*la6*sb+Zh11*Zh32*cb*la4
     &      *sb+Zh11*Zh32*cb*la5*sb+Zh12*Zh31*cb*la4*sb+Zh12
     &      *Zh31*cb*la5*sb-2D0*Zh12*Zh32*cb*la7*sb-2D0*Zh13
     &      *Zh33*cb*la4s*sb-2D0*Zh13*Zh33*cb*la5s*sb-2D0
     &      *Zh13*Zh33*cb*la6s*sb
      AAABR(1827) = Zh11*Zh21*la1-Zh11*Zh22*la6-Zh12*Zh21
     &      *la6+Zh12*Zh22*la3+Zh13*Zh23*la1s-Zh11*Zh21*la1
     &      *sb**2+Zh11*Zh21*la3*sb**2+Zh11*Zh22*la6*sb**2
     &      -Zh11*Zh22*la7*sb**2+Zh12*Zh21*la6*sb**2-Zh12
     &      *Zh21*la7*sb**2+Zh12*Zh22*la2*sb**2-Zh12*Zh22
     &      *la3*sb**2-Zh13*Zh23*la1s*sb**2+Zh13*Zh23*la2s
     &      *sb**2-2D0*Zh11*Zh21*cb*la6*sb+Zh11*Zh22*cb*la4
     &      *sb+Zh11*Zh22*cb*la5*sb+Zh12*Zh21*cb*la4*sb+Zh12
     &      *Zh21*cb*la5*sb-2D0*Zh12*Zh22*cb*la7*sb-2D0*Zh13
     &      *Zh23*cb*la4s*sb-2D0*Zh13*Zh23*cb*la5s*sb-2D0
     &      *Zh13*Zh23*cb*la6s*sb
      AAABR(1828) = Zh11**2*la1+Zh12**2*la3+Zh13**2*la1s
     &      -2D0*Zh11*Zh12*la6-Zh11**2*la1*sb**2+Zh11**2
     &      *la3*sb**2+Zh12**2*la2*sb**2-Zh12**2*la3*sb**2
     &      -Zh13**2*la1s*sb**2+Zh13**2*la2s*sb**2+2D0*Zh11
     &      *Zh12*la6*sb**2-2D0*Zh11*Zh12*la7*sb**2-2D0*Zh11
     &      **2*cb*la6*sb-2D0*Zh12**2*cb*la7*sb-2D0*Zh13
     &      **2*cb*la4s*sb-2D0*Zh13**2*cb*la5s*sb-2D0*Zh13
     &      **2*cb*la6s*sb+2D0*Zh11*Zh12*cb*la4*sb+2D0*Zh11
     &      *Zh12*cb*la5*sb
      AAABR(1829) = Zh21**2*la1+Zh22**2*la3+Zh23**2*la1s
     &      -2D0*Zh21*Zh22*la6-Zh21**2*la1*sb**2+Zh21**2
     &      *la3*sb**2+Zh22**2*la2*sb**2-Zh22**2*la3*sb**2
     &      -Zh23**2*la1s*sb**2+Zh23**2*la2s*sb**2+2D0*Zh21
     &      *Zh22*la6*sb**2-2D0*Zh21*Zh22*la7*sb**2-2D0*Zh21
     &      **2*cb*la6*sb-2D0*Zh22**2*cb*la7*sb-2D0*Zh23
     &      **2*cb*la4s*sb-2D0*Zh23**2*cb*la5s*sb-2D0*Zh23
     &      **2*cb*la6s*sb+2D0*Zh21*Zh22*cb*la4*sb+2D0*Zh21
     &      *Zh22*cb*la5*sb
      AAABR(1830) = Zh21*Zh31*la1-Zh21*Zh32*la6-Zh22*Zh31
     &      *la6+Zh22*Zh32*la3+Zh23*Zh33*la1s-Zh21*Zh31*la1
     &      *sb**2+Zh21*Zh31*la3*sb**2+Zh21*Zh32*la6*sb**2
     &      -Zh21*Zh32*la7*sb**2+Zh22*Zh31*la6*sb**2-Zh22
     &      *Zh31*la7*sb**2+Zh22*Zh32*la2*sb**2-Zh22*Zh32
     &      *la3*sb**2-Zh23*Zh33*la1s*sb**2+Zh23*Zh33*la2s
     &      *sb**2-2D0*Zh21*Zh31*cb*la6*sb+Zh21*Zh32*cb*la4
     &      *sb+Zh21*Zh32*cb*la5*sb+Zh22*Zh31*cb*la4*sb+Zh22
     &      *Zh31*cb*la5*sb-2D0*Zh22*Zh32*cb*la7*sb-2D0*Zh23
     &      *Zh33*cb*la4s*sb-2D0*Zh23*Zh33*cb*la5s*sb-2D0
     &      *Zh23*Zh33*cb*la6s*sb
      AAABR(1831) = Zh31**2*la1+Zh32**2*la3+Zh33**2*la1s
     &      -2D0*Zh31*Zh32*la6-Zh31**2*la1*sb**2+Zh31**2
     &      *la3*sb**2+Zh32**2*la2*sb**2-Zh32**2*la3*sb**2
     &      -Zh33**2*la1s*sb**2+Zh33**2*la2s*sb**2+2D0*Zh31
     &      *Zh32*la6*sb**2-2D0*Zh31*Zh32*la7*sb**2-2D0*Zh31
     &      **2*cb*la6*sb-2D0*Zh32**2*cb*la7*sb-2D0*Zh33
     &      **2*cb*la4s*sb-2D0*Zh33**2*cb*la5s*sb-2D0*Zh33
     &      **2*cb*la6s*sb+2D0*Zh31*Zh32*cb*la4*sb+2D0*Zh31
     &      *Zh32*cb*la5*sb
      AAABR(1832) = Pa21**2*la3+Pa21**2*la4+Pa21**2*la5+Pa22
     &      **2*la1s+3D0*Pa21**2*la1*sb**2-3D0*Pa21**2*la1
     &      *sb**4+3D0*Pa21**2*la2*sb**2-3D0*Pa21**2*la2
     &      *sb**4-6D0*Pa21**2*la3*sb**2+6D0*Pa21**2*la3
     &      *sb**4-6D0*Pa21**2*la4*sb**2+6D0*Pa21**2*la4
     &      *sb**4-6D0*Pa21**2*la5*sb**2+6D0*Pa21**2*la5
     &      *sb**4-Pa22**2*la1s*sb**2+Pa22**2*la2s*sb**2
     &      +6D0*Pa21**2*cb*la6*sb-12D0*Pa21**2*cb*la6*sb
     &      **3-6D0*Pa21**2*cb*la7*sb+12D0*Pa21**2*cb*la7
     &      *sb**3-2D0*Pa22**2*cb*la4s*sb+2D0*Pa22**2*cb
     &      *la5s*sb+2D0*Pa22**2*cb*la6s*sb
      AAABR(1833) = Pa11**2*la3+Pa11**2*la4+Pa11**2*la5+Pa12
     &      **2*la1s+3D0*Pa11**2*la1*sb**2-3D0*Pa11**2*la1
     &      *sb**4+3D0*Pa11**2*la2*sb**2-3D0*Pa11**2*la2
     &      *sb**4-6D0*Pa11**2*la3*sb**2+6D0*Pa11**2*la3
     &      *sb**4-6D0*Pa11**2*la4*sb**2+6D0*Pa11**2*la4
     &      *sb**4-6D0*Pa11**2*la5*sb**2+6D0*Pa11**2*la5
     &      *sb**4-Pa12**2*la1s*sb**2+Pa12**2*la2s*sb**2
     &      +6D0*Pa11**2*cb*la6*sb-12D0*Pa11**2*cb*la6*sb
     &      **3-6D0*Pa11**2*cb*la7*sb+12D0*Pa11**2*cb*la7
     &      *sb**3-2D0*Pa12**2*cb*la4s*sb+2D0*Pa12**2*cb
     &      *la5s*sb+2D0*Pa12**2*cb*la6s*sb
      AAABR(1834) = Pa11*Pa21*la3+Pa11*Pa21*la4+Pa11*Pa21
     &      *la5+Pa12*Pa22*la1s+3D0*Pa11*Pa21*la1*sb**2-3D0
     &      *Pa11*Pa21*la1*sb**4+3D0*Pa11*Pa21*la2*sb**2
     &      -3D0*Pa11*Pa21*la2*sb**4-6D0*Pa11*Pa21*la3*sb
     &      **2+6D0*Pa11*Pa21*la3*sb**4-6D0*Pa11*Pa21*la4
     &      *sb**2+6D0*Pa11*Pa21*la4*sb**4-6D0*Pa11*Pa21
     &      *la5*sb**2+6D0*Pa11*Pa21*la5*sb**4-Pa12*Pa22
     &      *la1s*sb**2+Pa12*Pa22*la2s*sb**2+6D0*Pa11*Pa21
     &      *cb*la6*sb-12D0*Pa11*Pa21*cb*la6*sb**3-6D0*Pa11
     &      *Pa21*cb*la7*sb+12D0*Pa11*Pa21*cb*la7*sb**3-2D0
     &      *Pa12*Pa22*cb*la4s*sb+2D0*Pa12*Pa22*cb*la5s*sb
     &      +2D0*Pa12*Pa22*cb*la6s*sb
      AAABR(1835) = Zh31**2*la1+Zh32**2*la3+Zh32**2*la4-Zh32
     &      **2*la5+Zh33**2*la1s-2D0*Zh31*Zh32*la6-Zh31**2
     &      *la1*sb**2+Zh31**2*la3*sb**2+Zh31**2*la4*sb**2
     &      -Zh31**2*la5*sb**2+Zh32**2*la2*sb**2-Zh32**2
     &      *la3*sb**2-Zh32**2*la4*sb**2+Zh32**2*la5*sb**2
     &      -Zh33**2*la1s*sb**2+Zh33**2*la2s*sb**2+2D0*Zh31
     &      *Zh32*la6*sb**2-2D0*Zh31*Zh32*la7*sb**2-2D0*Zh31
     &      **2*cb*la6*sb-2D0*Zh32**2*cb*la7*sb-2D0*Zh33
     &      **2*cb*la4s*sb-2D0*Zh33**2*cb*la5s*sb-2D0*Zh33
     &      **2*cb*la6s*sb+4D0*Zh31*Zh32*cb*la5*sb
      AAABR(1836) = Zh11**2*la1+Zh12**2*la3+Zh12**2*la4-Zh12
     &      **2*la5+Zh13**2*la1s-2D0*Zh11*Zh12*la6-Zh11**2
     &      *la1*sb**2+Zh11**2*la3*sb**2+Zh11**2*la4*sb**2
     &      -Zh11**2*la5*sb**2+Zh12**2*la2*sb**2-Zh12**2
     &      *la3*sb**2-Zh12**2*la4*sb**2+Zh12**2*la5*sb**2
     &      -Zh13**2*la1s*sb**2+Zh13**2*la2s*sb**2+2D0*Zh11
     &      *Zh12*la6*sb**2-2D0*Zh11*Zh12*la7*sb**2-2D0*Zh11
     &      **2*cb*la6*sb-2D0*Zh12**2*cb*la7*sb-2D0*Zh13
     &      **2*cb*la4s*sb-2D0*Zh13**2*cb*la5s*sb-2D0*Zh13
     &      **2*cb*la6s*sb+4D0*Zh11*Zh12*cb*la5*sb
      AAABR(1837) = Zh21**2*la1+Zh22**2*la3+Zh22**2*la4-Zh22
     &      **2*la5+Zh23**2*la1s-2D0*Zh21*Zh22*la6-Zh21**2
     &      *la1*sb**2+Zh21**2*la3*sb**2+Zh21**2*la4*sb**2
     &      -Zh21**2*la5*sb**2+Zh22**2*la2*sb**2-Zh22**2
     &      *la3*sb**2-Zh22**2*la4*sb**2+Zh22**2*la5*sb**2
     &      -Zh23**2*la1s*sb**2+Zh23**2*la2s*sb**2+2D0*Zh21
     &      *Zh22*la6*sb**2-2D0*Zh21*Zh22*la7*sb**2-2D0*Zh21
     &      **2*cb*la6*sb-2D0*Zh22**2*cb*la7*sb-2D0*Zh23
     &      **2*cb*la4s*sb-2D0*Zh23**2*cb*la5s*sb-2D0*Zh23
     &      **2*cb*la6s*sb+4D0*Zh21*Zh22*cb*la5*sb
      AAABR(1838) = Zh11*Zh21*la1-Zh11*Zh22*la6-Zh12*Zh21
     &      *la6+Zh12*Zh22*la3+Zh12*Zh22*la4-Zh12*Zh22*la5
     &      +Zh13*Zh23*la1s-Zh11*Zh21*la1*sb**2+Zh11*Zh21
     &      *la3*sb**2+Zh11*Zh21*la4*sb**2-Zh11*Zh21*la5
     &      *sb**2+Zh11*Zh22*la6*sb**2-Zh11*Zh22*la7*sb**2
     &      +Zh12*Zh21*la6*sb**2-Zh12*Zh21*la7*sb**2+Zh12
     &      *Zh22*la2*sb**2-Zh12*Zh22*la3*sb**2-Zh12*Zh22
     &      *la4*sb**2+Zh12*Zh22*la5*sb**2-Zh13*Zh23*la1s
     &      *sb**2+Zh13*Zh23*la2s*sb**2-2D0*Zh11*Zh21*cb
     &      *la6*sb+2D0*Zh11*Zh22*cb*la5*sb+2D0*Zh12*Zh21
     &      *cb*la5*sb-2D0*Zh12*Zh22*cb*la7*sb-2D0*Zh13*Zh23
     &      *cb*la4s*sb-2D0*Zh13*Zh23*cb*la5s*sb-2D0*Zh13
     &      *Zh23*cb*la6s*sb
      AAABR(1839) = Zh21*Zh31*la1-Zh21*Zh32*la6-Zh22*Zh31
     &      *la6+Zh22*Zh32*la3+Zh22*Zh32*la4-Zh22*Zh32*la5
     &      +Zh23*Zh33*la1s-Zh21*Zh31*la1*sb**2+Zh21*Zh31
     &      *la3*sb**2+Zh21*Zh31*la4*sb**2-Zh21*Zh31*la5
     &      *sb**2+Zh21*Zh32*la6*sb**2-Zh21*Zh32*la7*sb**2
     &      +Zh22*Zh31*la6*sb**2-Zh22*Zh31*la7*sb**2+Zh22
     &      *Zh32*la2*sb**2-Zh22*Zh32*la3*sb**2-Zh22*Zh32
     &      *la4*sb**2+Zh22*Zh32*la5*sb**2-Zh23*Zh33*la1s
     &      *sb**2+Zh23*Zh33*la2s*sb**2-2D0*Zh21*Zh31*cb
     &      *la6*sb+2D0*Zh21*Zh32*cb*la5*sb+2D0*Zh22*Zh31
     &      *cb*la5*sb-2D0*Zh22*Zh32*cb*la7*sb-2D0*Zh23*Zh33
     &      *cb*la4s*sb-2D0*Zh23*Zh33*cb*la5s*sb-2D0*Zh23
     &      *Zh33*cb*la6s*sb
      AAABR(1840) = Zh11*Zh31*la1-Zh11*Zh32*la6-Zh12*Zh31
     &      *la6+Zh12*Zh32*la3+Zh12*Zh32*la4-Zh12*Zh32*la5
     &      +Zh13*Zh33*la1s-Zh11*Zh31*la1*sb**2+Zh11*Zh31
     &      *la3*sb**2+Zh11*Zh31*la4*sb**2-Zh11*Zh31*la5
     &      *sb**2+Zh11*Zh32*la6*sb**2-Zh11*Zh32*la7*sb**2
     &      +Zh12*Zh31*la6*sb**2-Zh12*Zh31*la7*sb**2+Zh12
     &      *Zh32*la2*sb**2-Zh12*Zh32*la3*sb**2-Zh12*Zh32
     &      *la4*sb**2+Zh12*Zh32*la5*sb**2-Zh13*Zh33*la1s
     &      *sb**2+Zh13*Zh33*la2s*sb**2-2D0*Zh11*Zh31*cb
     &      *la6*sb+2D0*Zh11*Zh32*cb*la5*sb+2D0*Zh12*Zh31
     &      *cb*la5*sb-2D0*Zh12*Zh32*cb*la7*sb-2D0*Zh13*Zh33
     &      *cb*la4s*sb-2D0*Zh13*Zh33*cb*la5s*sb-2D0*Zh13
     &      *Zh33*cb*la6s*sb
      AAABR(1841) = Pa21**3*la7+Pa21*Pa22**2*la4s-Pa21*Pa22
     &      **2*la5s-Pa21*Pa22**2*la6s+3D0*Pa21**3*la6*sb
     &      **2-4D0*Pa21**3*la6*sb**4-5D0*Pa21**3*la7*sb
     &      **2+4D0*Pa21**3*la7*sb**4-2D0*Pa21*Pa22**2*la4s
     &      *sb**2+2D0*Pa21*Pa22**2*la5s*sb**2+2D0*Pa21*Pa22
     &      **2*la6s*sb**2+Pa21**3*cb*la1*sb**3-Pa21**3*cb
     &      *la2*sb+Pa21**3*cb*la2*sb**3+Pa21**3*cb*la3*sb
     &      -2D0*Pa21**3*cb*la3*sb**3+Pa21**3*cb*la4*sb-2D0
     &      *Pa21**3*cb*la4*sb**3+Pa21**3*cb*la5*sb-2D0*Pa21
     &      **3*cb*la5*sb**3+Pa21*Pa22**2*cb*la1s*sb-Pa21
     &      *Pa22**2*cb*la2s*sb
      AAABR(1842) = 3D0*Pa11*Pa21**2*la7+Pa11*Pa22**2*la4s
     &      -Pa11*Pa22**2*la5s-Pa11*Pa22**2*la6s+9D0*Pa11
     &      *Pa21**2*la6*sb**2-12D0*Pa11*Pa21**2*la6*sb**4
     &      -15D0*Pa11*Pa21**2*la7*sb**2+12D0*Pa11*Pa21**2
     &      *la7*sb**4-2D0*Pa11*Pa22**2*la4s*sb**2+2D0*Pa11
     &      *Pa22**2*la5s*sb**2+2D0*Pa11*Pa22**2*la6s*sb
     &      **2+2D0*Pa12*Pa21*Pa22*la4s-2D0*Pa12*Pa21*Pa22
     &      *la5s-2D0*Pa12*Pa21*Pa22*la6s+3D0*Pa11*Pa21**2
     &      *cb*la1*sb**3-3D0*Pa11*Pa21**2*cb*la2*sb+3D0
     &      *Pa11*Pa21**2*cb*la2*sb**3+3D0*Pa11*Pa21**2*cb
     &      *la3*sb-6D0*Pa11*Pa21**2*cb*la3*sb**3+3D0*Pa11
     &      *Pa21**2*cb*la4*sb-6D0*Pa11*Pa21**2*cb*la4*sb
     &      **3+3D0*Pa11*Pa21**2*cb*la5*sb-6D0*Pa11*Pa21
     &      **2*cb*la5*sb**3+Pa11*Pa22**2*cb*la1s*sb-Pa11
     &      *Pa22**2*cb*la2s*sb-4D0*Pa12*Pa21*Pa22*la4s*sb
     &      **2+4D0*Pa12*Pa21*Pa22*la5s*sb**2+4D0*Pa12*Pa21
     &      *Pa22*la6s*sb**2+2D0*Pa12*Pa21*Pa22*cb*la1s*sb
     &      -2D0*Pa12*Pa21*Pa22*cb*la2s*sb
      AAABR(1843) = Pa11**3*la7+Pa11*Pa12**2*la4s-Pa11*Pa12
     &      **2*la5s-Pa11*Pa12**2*la6s+3D0*Pa11**3*la6*sb
     &      **2-4D0*Pa11**3*la6*sb**4-5D0*Pa11**3*la7*sb
     &      **2+4D0*Pa11**3*la7*sb**4-2D0*Pa11*Pa12**2*la4s
     &      *sb**2+2D0*Pa11*Pa12**2*la5s*sb**2+2D0*Pa11*Pa12
     &      **2*la6s*sb**2+Pa11**3*cb*la1*sb**3-Pa11**3*cb
     &      *la2*sb+Pa11**3*cb*la2*sb**3+Pa11**3*cb*la3*sb
     &      -2D0*Pa11**3*cb*la3*sb**3+Pa11**3*cb*la4*sb-2D0
     &      *Pa11**3*cb*la4*sb**3+Pa11**3*cb*la5*sb-2D0*Pa11
     &      **3*cb*la5*sb**3+Pa11*Pa12**2*cb*la1s*sb-Pa11
     &      *Pa12**2*cb*la2s*sb
      AAABR(1844) = 3D0*Pa11**2*Pa21*la7+Pa12**2*Pa21*la4s
     &      -Pa12**2*Pa21*la5s-Pa12**2*Pa21*la6s+2D0*Pa11
     &      *Pa12*Pa22*la4s-2D0*Pa11*Pa12*Pa22*la5s-2D0*Pa11
     &      *Pa12*Pa22*la6s+9D0*Pa11**2*Pa21*la6*sb**2-12D0
     &      *Pa11**2*Pa21*la6*sb**4-15D0*Pa11**2*Pa21*la7
     &      *sb**2+12D0*Pa11**2*Pa21*la7*sb**4-2D0*Pa12**2
     &      *Pa21*la4s*sb**2+2D0*Pa12**2*Pa21*la5s*sb**2
     &      +2D0*Pa12**2*Pa21*la6s*sb**2-4D0*Pa11*Pa12*Pa22
     &      *la4s*sb**2+4D0*Pa11*Pa12*Pa22*la5s*sb**2+4D0
     &      *Pa11*Pa12*Pa22*la6s*sb**2+3D0*Pa11**2*Pa21*cb
     &      *la1*sb**3-3D0*Pa11**2*Pa21*cb*la2*sb+3D0*Pa11
     &      **2*Pa21*cb*la2*sb**3+3D0*Pa11**2*Pa21*cb*la3
     &      *sb-6D0*Pa11**2*Pa21*cb*la3*sb**3+3D0*Pa11**2
     &      *Pa21*cb*la4*sb-6D0*Pa11**2*Pa21*cb*la4*sb**3
     &      +3D0*Pa11**2*Pa21*cb*la5*sb-6D0*Pa11**2*Pa21
     &      *cb*la5*sb**3+Pa12**2*Pa21*cb*la1s*sb-Pa12**2
     &      *Pa21*cb*la2s*sb+2D0*Pa11*Pa12*Pa22*cb*la1s*sb
     &      -2D0*Pa11*Pa12*Pa22*cb*la2s*sb
      AAABR(1845) = Pa11*Zh31**2*la6+Pa11*Zh32**2*la7+Pa11
     &      *Zh33**2*la4s+Pa11*Zh33**2*la5s+Pa11*Zh33**2
     &      *la6s-2D0*Pa11*Zh31*Zh32*la5-2D0*Pa11*Zh31**2
     &      *la6*sb**2-2D0*Pa11*Zh32**2*la7*sb**2-2D0*Pa11
     &      *Zh33**2*la4s*sb**2-2D0*Pa11*Zh33**2*la5s*sb
     &      **2-2D0*Pa11*Zh33**2*la6s*sb**2+4D0*Pa11*Zh31
     &      *Zh32*la5*sb**2+Pa11*Zh31**2*cb*la1*sb-Pa11*Zh31
     &      **2*cb*la3*sb-Pa11*Zh31**2*cb*la4*sb+Pa11*Zh31
     &      **2*cb*la5*sb-Pa11*Zh32**2*cb*la2*sb+Pa11*Zh32
     &      **2*cb*la3*sb+Pa11*Zh32**2*cb*la4*sb-Pa11*Zh32
     &      **2*cb*la5*sb+Pa11*Zh33**2*cb*la1s*sb-Pa11*Zh33
     &      **2*cb*la2s*sb-2D0*Pa12*Zh31*Zh33*la5s*sb+2D0
     &      *Pa12*Zh31*Zh33*la6s*sb+2D0*Pa12*Zh32*Zh33*cb
     &      *la5s-2D0*Pa12*Zh32*Zh33*cb*la6s-2D0*Pa11*Zh31
     &      *Zh32*cb*la6*sb+2D0*Pa11*Zh31*Zh32*cb*la7*sb
      AAABR(1846) = Pa21*Zh11*Zh31*la6-Pa21*Zh11*Zh32*la5
     &      -Pa21*Zh12*Zh31*la5+Pa21*Zh12*Zh32*la7+Pa21*Zh13
     &      *Zh33*la4s+Pa21*Zh13*Zh33*la5s+Pa21*Zh13*Zh33
     &      *la6s-2D0*Pa21*Zh11*Zh31*la6*sb**2+2D0*Pa21*Zh11
     &      *Zh32*la5*sb**2+2D0*Pa21*Zh12*Zh31*la5*sb**2
     &      -2D0*Pa21*Zh12*Zh32*la7*sb**2-2D0*Pa21*Zh13*Zh33
     &      *la4s*sb**2-2D0*Pa21*Zh13*Zh33*la5s*sb**2-2D0
     &      *Pa21*Zh13*Zh33*la6s*sb**2-Pa22*Zh11*Zh33*la5s
     &      *sb+Pa22*Zh11*Zh33*la6s*sb+Pa22*Zh12*Zh33*cb
     &      *la5s-Pa22*Zh12*Zh33*cb*la6s-Pa22*Zh13*Zh31*la5s
     &      *sb+Pa22*Zh13*Zh31*la6s*sb+Pa22*Zh13*Zh32*cb
     &      *la5s-Pa22*Zh13*Zh32*cb*la6s+Pa21*Zh11*Zh31*cb
     &      *la1*sb-Pa21*Zh11*Zh31*cb*la3*sb-Pa21*Zh11*Zh31
     &      *cb*la4*sb+Pa21*Zh11*Zh31*cb*la5*sb-Pa21*Zh11
     &      *Zh32*cb*la6*sb+Pa21*Zh11*Zh32*cb*la7*sb-Pa21
     &      *Zh12*Zh31*cb*la6*sb+Pa21*Zh12*Zh31*cb*la7*sb
     &      -Pa21*Zh12*Zh32*cb*la2*sb+Pa21*Zh12*Zh32*cb*la3
     &      *sb+Pa21*Zh12*Zh32*cb*la4*sb-Pa21*Zh12*Zh32*cb
     &      *la5*sb+Pa21*Zh13*Zh33*cb*la1s*sb-Pa21*Zh13*Zh33
     &      *cb*la2s*sb
      AAABR(1847) = Pa11*Zh21*Zh31*la6-Pa11*Zh21*Zh32*la5
     &      -Pa11*Zh22*Zh31*la5+Pa11*Zh22*Zh32*la7+Pa11*Zh23
     &      *Zh33*la4s+Pa11*Zh23*Zh33*la5s+Pa11*Zh23*Zh33
     &      *la6s-2D0*Pa11*Zh21*Zh31*la6*sb**2+2D0*Pa11*Zh21
     &      *Zh32*la5*sb**2+2D0*Pa11*Zh22*Zh31*la5*sb**2
     &      -2D0*Pa11*Zh22*Zh32*la7*sb**2-2D0*Pa11*Zh23*Zh33
     &      *la4s*sb**2-2D0*Pa11*Zh23*Zh33*la5s*sb**2-2D0
     &      *Pa11*Zh23*Zh33*la6s*sb**2-Pa12*Zh21*Zh33*la5s
     &      *sb+Pa12*Zh21*Zh33*la6s*sb+Pa12*Zh22*Zh33*cb
     &      *la5s-Pa12*Zh22*Zh33*cb*la6s-Pa12*Zh23*Zh31*la5s
     &      *sb+Pa12*Zh23*Zh31*la6s*sb+Pa12*Zh23*Zh32*cb
     &      *la5s-Pa12*Zh23*Zh32*cb*la6s+Pa11*Zh21*Zh31*cb
     &      *la1*sb-Pa11*Zh21*Zh31*cb*la3*sb-Pa11*Zh21*Zh31
     &      *cb*la4*sb+Pa11*Zh21*Zh31*cb*la5*sb-Pa11*Zh21
     &      *Zh32*cb*la6*sb+Pa11*Zh21*Zh32*cb*la7*sb-Pa11
     &      *Zh22*Zh31*cb*la6*sb+Pa11*Zh22*Zh31*cb*la7*sb
     &      -Pa11*Zh22*Zh32*cb*la2*sb+Pa11*Zh22*Zh32*cb*la3
     &      *sb+Pa11*Zh22*Zh32*cb*la4*sb-Pa11*Zh22*Zh32*cb
     &      *la5*sb+Pa11*Zh23*Zh33*cb*la1s*sb-Pa11*Zh23*Zh33
     &      *cb*la2s*sb
      AAABR(1848) = Pa11*Zh11*Zh31*la6-Pa11*Zh11*Zh32*la5
     &      -Pa11*Zh12*Zh31*la5+Pa11*Zh12*Zh32*la7+Pa11*Zh13
     &      *Zh33*la4s+Pa11*Zh13*Zh33*la5s+Pa11*Zh13*Zh33
     &      *la6s-2D0*Pa11*Zh11*Zh31*la6*sb**2+2D0*Pa11*Zh11
     &      *Zh32*la5*sb**2+2D0*Pa11*Zh12*Zh31*la5*sb**2
     &      -2D0*Pa11*Zh12*Zh32*la7*sb**2-2D0*Pa11*Zh13*Zh33
     &      *la4s*sb**2-2D0*Pa11*Zh13*Zh33*la5s*sb**2-2D0
     &      *Pa11*Zh13*Zh33*la6s*sb**2-Pa12*Zh11*Zh33*la5s
     &      *sb+Pa12*Zh11*Zh33*la6s*sb+Pa12*Zh12*Zh33*cb
     &      *la5s-Pa12*Zh12*Zh33*cb*la6s-Pa12*Zh13*Zh31*la5s
     &      *sb+Pa12*Zh13*Zh31*la6s*sb+Pa12*Zh13*Zh32*cb
     &      *la5s-Pa12*Zh13*Zh32*cb*la6s+Pa11*Zh11*Zh31*cb
     &      *la1*sb-Pa11*Zh11*Zh31*cb*la3*sb-Pa11*Zh11*Zh31
     &      *cb*la4*sb+Pa11*Zh11*Zh31*cb*la5*sb-Pa11*Zh11
     &      *Zh32*cb*la6*sb+Pa11*Zh11*Zh32*cb*la7*sb-Pa11
     &      *Zh12*Zh31*cb*la6*sb+Pa11*Zh12*Zh31*cb*la7*sb
     &      -Pa11*Zh12*Zh32*cb*la2*sb+Pa11*Zh12*Zh32*cb*la3
     &      *sb+Pa11*Zh12*Zh32*cb*la4*sb-Pa11*Zh12*Zh32*cb
     &      *la5*sb+Pa11*Zh13*Zh33*cb*la1s*sb-Pa11*Zh13*Zh33
     &      *cb*la2s*sb
      AAABR(1849) = Pa11*Zh11*Zh21*la6-Pa11*Zh11*Zh22*la5
     &      -Pa11*Zh12*Zh21*la5+Pa11*Zh12*Zh22*la7+Pa11*Zh13
     &      *Zh23*la4s+Pa11*Zh13*Zh23*la5s+Pa11*Zh13*Zh23
     &      *la6s-2D0*Pa11*Zh11*Zh21*la6*sb**2+2D0*Pa11*Zh11
     &      *Zh22*la5*sb**2+2D0*Pa11*Zh12*Zh21*la5*sb**2
     &      -2D0*Pa11*Zh12*Zh22*la7*sb**2-2D0*Pa11*Zh13*Zh23
     &      *la4s*sb**2-2D0*Pa11*Zh13*Zh23*la5s*sb**2-2D0
     &      *Pa11*Zh13*Zh23*la6s*sb**2-Pa12*Zh11*Zh23*la5s
     &      *sb+Pa12*Zh11*Zh23*la6s*sb+Pa12*Zh12*Zh23*cb
     &      *la5s-Pa12*Zh12*Zh23*cb*la6s-Pa12*Zh13*Zh21*la5s
     &      *sb+Pa12*Zh13*Zh21*la6s*sb+Pa12*Zh13*Zh22*cb
     &      *la5s-Pa12*Zh13*Zh22*cb*la6s+Pa11*Zh11*Zh21*cb
     &      *la1*sb-Pa11*Zh11*Zh21*cb*la3*sb-Pa11*Zh11*Zh21
     &      *cb*la4*sb+Pa11*Zh11*Zh21*cb*la5*sb-Pa11*Zh11
     &      *Zh22*cb*la6*sb+Pa11*Zh11*Zh22*cb*la7*sb-Pa11
     &      *Zh12*Zh21*cb*la6*sb+Pa11*Zh12*Zh21*cb*la7*sb
     &      -Pa11*Zh12*Zh22*cb*la2*sb+Pa11*Zh12*Zh22*cb*la3
     &      *sb+Pa11*Zh12*Zh22*cb*la4*sb-Pa11*Zh12*Zh22*cb
     &      *la5*sb+Pa11*Zh13*Zh23*cb*la1s*sb-Pa11*Zh13*Zh23
     &      *cb*la2s*sb
      AAABR(1850) = Pa21*Zh11*Zh21*la6-Pa21*Zh11*Zh22*la5
     &      -Pa21*Zh12*Zh21*la5+Pa21*Zh12*Zh22*la7+Pa21*Zh13
     &      *Zh23*la4s+Pa21*Zh13*Zh23*la5s+Pa21*Zh13*Zh23
     &      *la6s-2D0*Pa21*Zh11*Zh21*la6*sb**2+2D0*Pa21*Zh11
     &      *Zh22*la5*sb**2+2D0*Pa21*Zh12*Zh21*la5*sb**2
     &      -2D0*Pa21*Zh12*Zh22*la7*sb**2-2D0*Pa21*Zh13*Zh23
     &      *la4s*sb**2-2D0*Pa21*Zh13*Zh23*la5s*sb**2-2D0
     &      *Pa21*Zh13*Zh23*la6s*sb**2-Pa22*Zh11*Zh23*la5s
     &      *sb+Pa22*Zh11*Zh23*la6s*sb+Pa22*Zh12*Zh23*cb
     &      *la5s-Pa22*Zh12*Zh23*cb*la6s-Pa22*Zh13*Zh21*la5s
     &      *sb+Pa22*Zh13*Zh21*la6s*sb+Pa22*Zh13*Zh22*cb
     &      *la5s-Pa22*Zh13*Zh22*cb*la6s+Pa21*Zh11*Zh21*cb
     &      *la1*sb-Pa21*Zh11*Zh21*cb*la3*sb-Pa21*Zh11*Zh21
     &      *cb*la4*sb+Pa21*Zh11*Zh21*cb*la5*sb-Pa21*Zh11
     &      *Zh22*cb*la6*sb+Pa21*Zh11*Zh22*cb*la7*sb-Pa21
     &      *Zh12*Zh21*cb*la6*sb+Pa21*Zh12*Zh21*cb*la7*sb
     &      -Pa21*Zh12*Zh22*cb*la2*sb+Pa21*Zh12*Zh22*cb*la3
     &      *sb+Pa21*Zh12*Zh22*cb*la4*sb-Pa21*Zh12*Zh22*cb
     &      *la5*sb+Pa21*Zh13*Zh23*cb*la1s*sb-Pa21*Zh13*Zh23
     &      *cb*la2s*sb
      AAABR(1851) = Pa11*Zh21**2*la6+Pa11*Zh22**2*la7+Pa11
     &      *Zh23**2*la4s+Pa11*Zh23**2*la5s+Pa11*Zh23**2
     &      *la6s-2D0*Pa11*Zh21*Zh22*la5-2D0*Pa11*Zh21**2
     &      *la6*sb**2-2D0*Pa11*Zh22**2*la7*sb**2-2D0*Pa11
     &      *Zh23**2*la4s*sb**2-2D0*Pa11*Zh23**2*la5s*sb
     &      **2-2D0*Pa11*Zh23**2*la6s*sb**2+4D0*Pa11*Zh21
     &      *Zh22*la5*sb**2+Pa11*Zh21**2*cb*la1*sb-Pa11*Zh21
     &      **2*cb*la3*sb-Pa11*Zh21**2*cb*la4*sb+Pa11*Zh21
     &      **2*cb*la5*sb-Pa11*Zh22**2*cb*la2*sb+Pa11*Zh22
     &      **2*cb*la3*sb+Pa11*Zh22**2*cb*la4*sb-Pa11*Zh22
     &      **2*cb*la5*sb+Pa11*Zh23**2*cb*la1s*sb-Pa11*Zh23
     &      **2*cb*la2s*sb-2D0*Pa12*Zh21*Zh23*la5s*sb+2D0
     &      *Pa12*Zh21*Zh23*la6s*sb+2D0*Pa12*Zh22*Zh23*cb
     &      *la5s-2D0*Pa12*Zh22*Zh23*cb*la6s-2D0*Pa11*Zh21
     &      *Zh22*cb*la6*sb+2D0*Pa11*Zh21*Zh22*cb*la7*sb
      AAABR(1852) = Pa21*Zh31**2*la6+Pa21*Zh32**2*la7+Pa21
     &      *Zh33**2*la4s+Pa21*Zh33**2*la5s+Pa21*Zh33**2
     &      *la6s-2D0*Pa21*Zh31*Zh32*la5-2D0*Pa21*Zh31**2
     &      *la6*sb**2-2D0*Pa21*Zh32**2*la7*sb**2-2D0*Pa21
     &      *Zh33**2*la4s*sb**2-2D0*Pa21*Zh33**2*la5s*sb
     &      **2-2D0*Pa21*Zh33**2*la6s*sb**2+4D0*Pa21*Zh31
     &      *Zh32*la5*sb**2+Pa21*Zh31**2*cb*la1*sb-Pa21*Zh31
     &      **2*cb*la3*sb-Pa21*Zh31**2*cb*la4*sb+Pa21*Zh31
     &      **2*cb*la5*sb-Pa21*Zh32**2*cb*la2*sb+Pa21*Zh32
     &      **2*cb*la3*sb+Pa21*Zh32**2*cb*la4*sb-Pa21*Zh32
     &      **2*cb*la5*sb+Pa21*Zh33**2*cb*la1s*sb-Pa21*Zh33
     &      **2*cb*la2s*sb-2D0*Pa22*Zh31*Zh33*la5s*sb+2D0
     &      *Pa22*Zh31*Zh33*la6s*sb+2D0*Pa22*Zh32*Zh33*cb
     &      *la5s-2D0*Pa22*Zh32*Zh33*cb*la6s-2D0*Pa21*Zh31
     &      *Zh32*cb*la6*sb+2D0*Pa21*Zh31*Zh32*cb*la7*sb
      AAABR(1853) = Pa21*Zh11**2*la6+Pa21*Zh12**2*la7+Pa21
     &      *Zh13**2*la4s+Pa21*Zh13**2*la5s+Pa21*Zh13**2
     &      *la6s-2D0*Pa21*Zh11*Zh12*la5-2D0*Pa21*Zh11**2
     &      *la6*sb**2-2D0*Pa21*Zh12**2*la7*sb**2-2D0*Pa21
     &      *Zh13**2*la4s*sb**2-2D0*Pa21*Zh13**2*la5s*sb
     &      **2-2D0*Pa21*Zh13**2*la6s*sb**2+4D0*Pa21*Zh11
     &      *Zh12*la5*sb**2+Pa21*Zh11**2*cb*la1*sb-Pa21*Zh11
     &      **2*cb*la3*sb-Pa21*Zh11**2*cb*la4*sb+Pa21*Zh11
     &      **2*cb*la5*sb-Pa21*Zh12**2*cb*la2*sb+Pa21*Zh12
     &      **2*cb*la3*sb+Pa21*Zh12**2*cb*la4*sb-Pa21*Zh12
     &      **2*cb*la5*sb+Pa21*Zh13**2*cb*la1s*sb-Pa21*Zh13
     &      **2*cb*la2s*sb-2D0*Pa22*Zh11*Zh13*la5s*sb+2D0
     &      *Pa22*Zh11*Zh13*la6s*sb+2D0*Pa22*Zh12*Zh13*cb
     &      *la5s-2D0*Pa22*Zh12*Zh13*cb*la6s-2D0*Pa21*Zh11
     &      *Zh12*cb*la6*sb+2D0*Pa21*Zh11*Zh12*cb*la7*sb
      AAABR(1854) = Pa21*Zh21**2*la6+Pa21*Zh22**2*la7+Pa21
     &      *Zh23**2*la4s+Pa21*Zh23**2*la5s+Pa21*Zh23**2
     &      *la6s-2D0*Pa21*Zh21*Zh22*la5-2D0*Pa21*Zh21**2
     &      *la6*sb**2-2D0*Pa21*Zh22**2*la7*sb**2-2D0*Pa21
     &      *Zh23**2*la4s*sb**2-2D0*Pa21*Zh23**2*la5s*sb
     &      **2-2D0*Pa21*Zh23**2*la6s*sb**2+4D0*Pa21*Zh21
     &      *Zh22*la5*sb**2+Pa21*Zh21**2*cb*la1*sb-Pa21*Zh21
     &      **2*cb*la3*sb-Pa21*Zh21**2*cb*la4*sb+Pa21*Zh21
     &      **2*cb*la5*sb-Pa21*Zh22**2*cb*la2*sb+Pa21*Zh22
     &      **2*cb*la3*sb+Pa21*Zh22**2*cb*la4*sb-Pa21*Zh22
     &      **2*cb*la5*sb+Pa21*Zh23**2*cb*la1s*sb-Pa21*Zh23
     &      **2*cb*la2s*sb-2D0*Pa22*Zh21*Zh23*la5s*sb+2D0
     &      *Pa22*Zh21*Zh23*la6s*sb+2D0*Pa22*Zh22*Zh23*cb
     &      *la5s-2D0*Pa22*Zh22*Zh23*cb*la6s-2D0*Pa21*Zh21
     &      *Zh22*cb*la6*sb+2D0*Pa21*Zh21*Zh22*cb*la7*sb
      AAABR(1855) = Pa21*Zh21*Zh31*la6-Pa21*Zh21*Zh32*la5
     &      -Pa21*Zh22*Zh31*la5+Pa21*Zh22*Zh32*la7+Pa21*Zh23
     &      *Zh33*la4s+Pa21*Zh23*Zh33*la5s+Pa21*Zh23*Zh33
     &      *la6s-2D0*Pa21*Zh21*Zh31*la6*sb**2+2D0*Pa21*Zh21
     &      *Zh32*la5*sb**2+2D0*Pa21*Zh22*Zh31*la5*sb**2
     &      -2D0*Pa21*Zh22*Zh32*la7*sb**2-2D0*Pa21*Zh23*Zh33
     &      *la4s*sb**2-2D0*Pa21*Zh23*Zh33*la5s*sb**2-2D0
     &      *Pa21*Zh23*Zh33*la6s*sb**2-Pa22*Zh21*Zh33*la5s
     &      *sb+Pa22*Zh21*Zh33*la6s*sb+Pa22*Zh22*Zh33*cb
     &      *la5s-Pa22*Zh22*Zh33*cb*la6s-Pa22*Zh23*Zh31*la5s
     &      *sb+Pa22*Zh23*Zh31*la6s*sb+Pa22*Zh23*Zh32*cb
     &      *la5s-Pa22*Zh23*Zh32*cb*la6s+Pa21*Zh21*Zh31*cb
     &      *la1*sb-Pa21*Zh21*Zh31*cb*la3*sb-Pa21*Zh21*Zh31
     &      *cb*la4*sb+Pa21*Zh21*Zh31*cb*la5*sb-Pa21*Zh21
     &      *Zh32*cb*la6*sb+Pa21*Zh21*Zh32*cb*la7*sb-Pa21
     &      *Zh22*Zh31*cb*la6*sb+Pa21*Zh22*Zh31*cb*la7*sb
     &      -Pa21*Zh22*Zh32*cb*la2*sb+Pa21*Zh22*Zh32*cb*la3
     &      *sb+Pa21*Zh22*Zh32*cb*la4*sb-Pa21*Zh22*Zh32*cb
     &      *la5*sb+Pa21*Zh23*Zh33*cb*la1s*sb-Pa21*Zh23*Zh33
     &      *cb*la2s*sb
      AAABR(1856) = Pa11*Zh11**2*la6+Pa11*Zh12**2*la7+Pa11
     &      *Zh13**2*la4s+Pa11*Zh13**2*la5s+Pa11*Zh13**2
     &      *la6s-2D0*Pa11*Zh11*Zh12*la5-2D0*Pa11*Zh11**2
     &      *la6*sb**2-2D0*Pa11*Zh12**2*la7*sb**2-2D0*Pa11
     &      *Zh13**2*la4s*sb**2-2D0*Pa11*Zh13**2*la5s*sb
     &      **2-2D0*Pa11*Zh13**2*la6s*sb**2+4D0*Pa11*Zh11
     &      *Zh12*la5*sb**2+Pa11*Zh11**2*cb*la1*sb-Pa11*Zh11
     &      **2*cb*la3*sb-Pa11*Zh11**2*cb*la4*sb+Pa11*Zh11
     &      **2*cb*la5*sb-Pa11*Zh12**2*cb*la2*sb+Pa11*Zh12
     &      **2*cb*la3*sb+Pa11*Zh12**2*cb*la4*sb-Pa11*Zh12
     &      **2*cb*la5*sb+Pa11*Zh13**2*cb*la1s*sb-Pa11*Zh13
     &      **2*cb*la2s*sb-2D0*Pa12*Zh11*Zh13*la5s*sb+2D0
     &      *Pa12*Zh11*Zh13*la6s*sb+2D0*Pa12*Zh12*Zh13*cb
     &      *la5s-2D0*Pa12*Zh12*Zh13*cb*la6s-2D0*Pa11*Zh11
     &      *Zh12*cb*la6*sb+2D0*Pa11*Zh11*Zh12*cb*la7*sb
      AAABR(1857) = Pa11**4*la2+2D0*Pa12**4*la3s+4D0*Pa12
     &      **4*la7s-4D0*Pa12**4*la8s+2D0*Pa11**2*Pa12**2
     &      *la2s+Pa11**4*la1*sb**4-2D0*Pa11**4*la2*sb**2
     &      +Pa11**4*la2*sb**4+2D0*Pa11**4*la3*sb**2-2D0
     &      *Pa11**4*la3*sb**4+2D0*Pa11**4*la4*sb**2-2D0
     &      *Pa11**4*la4*sb**4+2D0*Pa11**4*la5*sb**2-2D0
     &      *Pa11**4*la5*sb**4+2D0*Pa11**2*Pa12**2*la1s*sb
     &      **2-2D0*Pa11**2*Pa12**2*la2s*sb**2+4D0*Pa11**4
     &      *cb*la6*sb**3+4D0*Pa11**4*cb*la7*sb-4D0*Pa11
     &      **4*cb*la7*sb**3+4D0*Pa11**2*Pa12**2*cb*la4s
     &      *sb-4D0*Pa11**2*Pa12**2*cb*la5s*sb-4D0*Pa11**2
     &      *Pa12**2*cb*la6s*sb
      AAABR(1858) = Pa21**4*la2+2D0*Pa22**4*la3s+4D0*Pa22
     &      **4*la7s-4D0*Pa22**4*la8s+2D0*Pa21**2*Pa22**2
     &      *la2s+Pa21**4*la1*sb**4-2D0*Pa21**4*la2*sb**2
     &      +Pa21**4*la2*sb**4+2D0*Pa21**4*la3*sb**2-2D0
     &      *Pa21**4*la3*sb**4+2D0*Pa21**4*la4*sb**2-2D0
     &      *Pa21**4*la4*sb**4+2D0*Pa21**4*la5*sb**2-2D0
     &      *Pa21**4*la5*sb**4+2D0*Pa21**2*Pa22**2*la1s*sb
     &      **2-2D0*Pa21**2*Pa22**2*la2s*sb**2+4D0*Pa21**4
     &      *cb*la6*sb**3+4D0*Pa21**4*cb*la7*sb-4D0*Pa21
     &      **4*cb*la7*sb**3+4D0*Pa21**2*Pa22**2*cb*la4s
     &      *sb-4D0*Pa21**2*Pa22**2*cb*la5s*sb-4D0*Pa21**2
     &      *Pa22**2*cb*la6s*sb
      AAABR(1859) = Pa11**3*Pa21*la2+2D0*Pa12**3*Pa22*la3s
     &      +4D0*Pa12**3*Pa22*la7s-4D0*Pa12**3*Pa22*la8s
     &      +Pa11*Pa12**2*Pa21*la2s+Pa11**2*Pa12*Pa22*la2s
     &      +Pa11**3*Pa21*la1*sb**4-2D0*Pa11**3*Pa21*la2
     &      *sb**2+Pa11**3*Pa21*la2*sb**4+2D0*Pa11**3*Pa21
     &      *la3*sb**2-2D0*Pa11**3*Pa21*la3*sb**4+2D0*Pa11
     &      **3*Pa21*la4*sb**2-2D0*Pa11**3*Pa21*la4*sb**4
     &      +2D0*Pa11**3*Pa21*la5*sb**2-2D0*Pa11**3*Pa21
     &      *la5*sb**4+Pa11*Pa12**2*Pa21*la1s*sb**2-Pa11
     &      *Pa12**2*Pa21*la2s*sb**2+Pa11**2*Pa12*Pa22*la1s
     &      *sb**2-Pa11**2*Pa12*Pa22*la2s*sb**2+4D0*Pa11
     &      **3*Pa21*cb*la6*sb**3+4D0*Pa11**3*Pa21*cb*la7
     &      *sb-4D0*Pa11**3*Pa21*cb*la7*sb**3+2D0*Pa11*Pa12
     &      **2*Pa21*cb*la4s*sb-2D0*Pa11*Pa12**2*Pa21*cb
     &      *la5s*sb-2D0*Pa11*Pa12**2*Pa21*cb*la6s*sb+2D0
     &      *Pa11**2*Pa12*Pa22*cb*la4s*sb-2D0*Pa11**2*Pa12
     &      *Pa22*cb*la5s*sb-2D0*Pa11**2*Pa12*Pa22*cb*la6s
     &      *sb
      AAABR(1860) = Pa11*Pa21**3*la2+2D0*Pa12*Pa22**3*la3s
     &      +4D0*Pa12*Pa22**3*la7s-4D0*Pa12*Pa22**3*la8s
     &      +Pa11*Pa21*Pa22**2*la2s+Pa11*Pa21**3*la1*sb**4
     &      -2D0*Pa11*Pa21**3*la2*sb**2+Pa11*Pa21**3*la2
     &      *sb**4+2D0*Pa11*Pa21**3*la3*sb**2-2D0*Pa11*Pa21
     &      **3*la3*sb**4+2D0*Pa11*Pa21**3*la4*sb**2-2D0
     &      *Pa11*Pa21**3*la4*sb**4+2D0*Pa11*Pa21**3*la5
     &      *sb**2-2D0*Pa11*Pa21**3*la5*sb**4+Pa12*Pa21**2
     &      *Pa22*la2s+Pa11*Pa21*Pa22**2*la1s*sb**2-Pa11
     &      *Pa21*Pa22**2*la2s*sb**2+4D0*Pa11*Pa21**3*cb
     &      *la6*sb**3+4D0*Pa11*Pa21**3*cb*la7*sb-4D0*Pa11
     &      *Pa21**3*cb*la7*sb**3+Pa12*Pa21**2*Pa22*la1s
     &      *sb**2-Pa12*Pa21**2*Pa22*la2s*sb**2+2D0*Pa11
     &      *Pa21*Pa22**2*cb*la4s*sb-2D0*Pa11*Pa21*Pa22**2
     &      *cb*la5s*sb-2D0*Pa11*Pa21*Pa22**2*cb*la6s*sb
     &      +2D0*Pa12*Pa21**2*Pa22*cb*la4s*sb-2D0*Pa12*Pa21
     &      **2*Pa22*cb*la5s*sb-2D0*Pa12*Pa21**2*Pa22*cb
     &      *la6s*sb
      AAABR(1861) = 3D0*Pa11**2*Pa21**2*la2+Pa11**2*Pa22
     &      **2*la2s+Pa12**2*Pa21**2*la2s+6D0*Pa12**2*Pa22
     &      **2*la3s+12D0*Pa12**2*Pa22**2*la7s-12D0*Pa12
     &      **2*Pa22**2*la8s+3D0*Pa11**2*Pa21**2*la1*sb**4
     &      -6D0*Pa11**2*Pa21**2*la2*sb**2+3D0*Pa11**2*Pa21
     &      **2*la2*sb**4+6D0*Pa11**2*Pa21**2*la3*sb**2-6D0
     &      *Pa11**2*Pa21**2*la3*sb**4+6D0*Pa11**2*Pa21**2
     &      *la4*sb**2-6D0*Pa11**2*Pa21**2*la4*sb**4+6D0
     &      *Pa11**2*Pa21**2*la5*sb**2-6D0*Pa11**2*Pa21**2
     &      *la5*sb**4+Pa11**2*Pa22**2*la1s*sb**2-Pa11**2
     &      *Pa22**2*la2s*sb**2+Pa12**2*Pa21**2*la1s*sb**2
     &      -Pa12**2*Pa21**2*la2s*sb**2+4D0*Pa11*Pa12*Pa21
     &      *Pa22*la2s+12D0*Pa11**2*Pa21**2*cb*la6*sb**3
     &      +12D0*Pa11**2*Pa21**2*cb*la7*sb-12D0*Pa11**2
     &      *Pa21**2*cb*la7*sb**3+2D0*Pa11**2*Pa22**2*cb
     &      *la4s*sb-2D0*Pa11**2*Pa22**2*cb*la5s*sb-2D0*Pa11
     &      **2*Pa22**2*cb*la6s*sb+2D0*Pa12**2*Pa21**2*cb
     &      *la4s*sb-2D0*Pa12**2*Pa21**2*cb*la5s*sb-2D0*Pa12
     &      **2*Pa21**2*cb*la6s*sb+4D0*Pa11*Pa12*Pa21*Pa22
     &      *la1s*sb**2-4D0*Pa11*Pa12*Pa21*Pa22*la2s*sb**2
     &      +8D0*Pa11*Pa12*Pa21*Pa22*cb*la4s*sb-8D0*Pa11
     &      *Pa12*Pa21*Pa22*cb*la5s*sb-8D0*Pa11*Pa12*Pa21
     &      *Pa22*cb*la6s*sb
      AAABR(1862) = Pa11*Pa21*Zh31**2*la3+Pa11*Pa21*Zh31
     &      **2*la4-Pa11*Pa21*Zh31**2*la5+Pa11*Pa21*Zh32
     &      **2*la2+Pa11*Pa21*Zh33**2*la2s+Pa12*Pa22*Zh31
     &      **2*la1s+Pa12*Pa22*Zh32**2*la2s+2D0*Pa12*Pa22
     &      *Zh33**2*la3s-12D0*Pa12*Pa22*Zh33**2*la7s-2D0
     &      *Pa11*Pa21*Zh31*Zh32*la7+Pa11*Pa21*Zh31**2*la1
     &      *sb**2-Pa11*Pa21*Zh31**2*la3*sb**2-Pa11*Pa21
     &      *Zh31**2*la4*sb**2+Pa11*Pa21*Zh31**2*la5*sb**2
     &      -Pa11*Pa21*Zh32**2*la2*sb**2+Pa11*Pa21*Zh32**2
     &      *la3*sb**2+Pa11*Pa21*Zh32**2*la4*sb**2-Pa11*Pa21
     &      *Zh32**2*la5*sb**2+Pa11*Pa21*Zh33**2*la1s*sb
     &      **2-Pa11*Pa21*Zh33**2*la2s*sb**2-2D0*Pa12*Pa22
     &      *Zh31*Zh32*la4s+2D0*Pa12*Pa22*Zh31*Zh32*la5s
     &      +2D0*Pa12*Pa22*Zh31*Zh32*la6s-2D0*Pa11*Pa21*Zh31
     &      *Zh32*la6*sb**2+2D0*Pa11*Pa21*Zh31*Zh32*la7*sb
     &      **2+2D0*Pa11*Pa21*Zh31**2*cb*la6*sb+2D0*Pa11
     &      *Pa21*Zh32**2*cb*la7*sb+2D0*Pa11*Pa21*Zh33**2
     &      *cb*la4s*sb+2D0*Pa11*Pa21*Zh33**2*cb*la5s*sb
     &      +2D0*Pa11*Pa21*Zh33**2*cb*la6s*sb+2D0*Pa11*Pa22
     &      *Zh31*Zh33*cb*la5s-2D0*Pa11*Pa22*Zh31*Zh33*cb
     &      *la6s+2D0*Pa11*Pa22*Zh32*Zh33*la5s*sb-2D0*Pa11
     &      *Pa22*Zh32*Zh33*la6s*sb+2D0*Pa12*Pa21*Zh31*Zh33
     &      *cb*la5s-2D0*Pa12*Pa21*Zh31*Zh33*cb*la6s+2D0
     &      *Pa12*Pa21*Zh32*Zh33*la5s*sb-2D0*Pa12*Pa21*Zh32
     &      *Zh33*la6s*sb-4D0*Pa11*Pa21*Zh31*Zh32*cb*la5*sb
      AAABR(1863) = Pa11*Pa21*Zh11*Zh31*la3+Pa11*Pa21*Zh11
     &      *Zh31*la4-Pa11*Pa21*Zh11*Zh31*la5-Pa11*Pa21*Zh11
     &      *Zh32*la7-Pa11*Pa21*Zh12*Zh31*la7+Pa11*Pa21*Zh12
     &      *Zh32*la2+Pa11*Pa21*Zh13*Zh33*la2s+Pa12*Pa22
     &      *Zh11*Zh31*la1s-Pa12*Pa22*Zh11*Zh32*la4s+Pa12
     &      *Pa22*Zh11*Zh32*la5s+Pa12*Pa22*Zh11*Zh32*la6s
     &      -Pa12*Pa22*Zh12*Zh31*la4s+Pa12*Pa22*Zh12*Zh31
     &      *la5s+Pa12*Pa22*Zh12*Zh31*la6s+Pa12*Pa22*Zh12
     &      *Zh32*la2s+2D0*Pa12*Pa22*Zh13*Zh33*la3s-12D0
     &      *Pa12*Pa22*Zh13*Zh33*la7s+Pa11*Pa21*Zh11*Zh31
     &      *la1*sb**2-Pa11*Pa21*Zh11*Zh31*la3*sb**2-Pa11
     &      *Pa21*Zh11*Zh31*la4*sb**2+Pa11*Pa21*Zh11*Zh31
     &      *la5*sb**2-Pa11*Pa21*Zh11*Zh32*la6*sb**2+Pa11
     &      *Pa21*Zh11*Zh32*la7*sb**2-Pa11*Pa21*Zh12*Zh31
     &      *la6*sb**2+Pa11*Pa21*Zh12*Zh31*la7*sb**2-Pa11
     &      *Pa21*Zh12*Zh32*la2*sb**2+Pa11*Pa21*Zh12*Zh32
     &      *la3*sb**2+Pa11*Pa21*Zh12*Zh32*la4*sb**2-Pa11
     &      *Pa21*Zh12*Zh32*la5*sb**2+Pa11*Pa21*Zh13*Zh33
     &      *la1s*sb**2-Pa11*Pa21*Zh13*Zh33*la2s*sb**2+Pa11
     &      *Pa22*Zh11*Zh33*cb*la5s-Pa11*Pa22*Zh11*Zh33*cb
     &      *la6s+Pa11*Pa22*Zh12*Zh33*la5s*sb-Pa11*Pa22*Zh12
     &      *Zh33*la6s*sb+Pa11*Pa22*Zh13*Zh31*cb*la5s-Pa11
     &      *Pa22*Zh13*Zh31*cb*la6s+Pa11*Pa22*Zh13*Zh32*la5s
     &      *sb-Pa11*Pa22*Zh13*Zh32*la6s*sb+Pa12*Pa21*Zh11
     &      *Zh33*cb*la5s-Pa12*Pa21*Zh11*Zh33*cb*la6s+Pa12
     &      *Pa21*Zh12*Zh33*la5s*sb-Pa12*Pa21*Zh12*Zh33*la6s
     &      *sb+Pa12*Pa21*Zh13*Zh31*cb*la5s-Pa12*Pa21*Zh13
     &      *Zh31*cb*la6s+Pa12*Pa21*Zh13*Zh32*la5s*sb-Pa12
     &      *Pa21*Zh13*Zh32*la6s*sb+2D0*Pa11*Pa21*Zh11*Zh31
     &      *cb*la6*sb-2D0*Pa11*Pa21*Zh11*Zh32*cb*la5*sb
     &      -2D0*Pa11*Pa21*Zh12*Zh31*cb*la5*sb+2D0*Pa11*Pa21
     &      *Zh12*Zh32*cb*la7*sb+2D0*Pa11*Pa21*Zh13*Zh33
     &      *cb*la4s*sb+2D0*Pa11*Pa21*Zh13*Zh33*cb*la5s*sb
     &      +2D0*Pa11*Pa21*Zh13*Zh33*cb*la6s*sb
      AAABR(1864) = Pa21**2*Zh21*Zh31*la3+Pa21**2*Zh21*Zh31
     &      *la4-Pa21**2*Zh21*Zh31*la5-Pa21**2*Zh21*Zh32
     &      *la7-Pa21**2*Zh22*Zh31*la7+Pa21**2*Zh22*Zh32
     &      *la2+Pa21**2*Zh23*Zh33*la2s+Pa22**2*Zh21*Zh31
     &      *la1s-Pa22**2*Zh21*Zh32*la4s+Pa22**2*Zh21*Zh32
     &      *la5s+Pa22**2*Zh21*Zh32*la6s-Pa22**2*Zh22*Zh31
     &      *la4s+Pa22**2*Zh22*Zh31*la5s+Pa22**2*Zh22*Zh31
     &      *la6s+Pa22**2*Zh22*Zh32*la2s+2D0*Pa22**2*Zh23
     &      *Zh33*la3s-12D0*Pa22**2*Zh23*Zh33*la7s+Pa21**2
     &      *Zh21*Zh31*la1*sb**2-Pa21**2*Zh21*Zh31*la3*sb
     &      **2-Pa21**2*Zh21*Zh31*la4*sb**2+Pa21**2*Zh21
     &      *Zh31*la5*sb**2-Pa21**2*Zh21*Zh32*la6*sb**2+Pa21
     &      **2*Zh21*Zh32*la7*sb**2-Pa21**2*Zh22*Zh31*la6
     &      *sb**2+Pa21**2*Zh22*Zh31*la7*sb**2-Pa21**2*Zh22
     &      *Zh32*la2*sb**2+Pa21**2*Zh22*Zh32*la3*sb**2+Pa21
     &      **2*Zh22*Zh32*la4*sb**2-Pa21**2*Zh22*Zh32*la5
     &      *sb**2+Pa21**2*Zh23*Zh33*la1s*sb**2-Pa21**2*Zh23
     &      *Zh33*la2s*sb**2+2D0*Pa21*Pa22*Zh21*Zh33*cb*la5s
     &      -2D0*Pa21*Pa22*Zh21*Zh33*cb*la6s+2D0*Pa21*Pa22
     &      *Zh22*Zh33*la5s*sb-2D0*Pa21*Pa22*Zh22*Zh33*la6s
     &      *sb+2D0*Pa21*Pa22*Zh23*Zh31*cb*la5s-2D0*Pa21
     &      *Pa22*Zh23*Zh31*cb*la6s+2D0*Pa21*Pa22*Zh23*Zh32
     &      *la5s*sb-2D0*Pa21*Pa22*Zh23*Zh32*la6s*sb+2D0
     &      *Pa21**2*Zh21*Zh31*cb*la6*sb-2D0*Pa21**2*Zh21
     &      *Zh32*cb*la5*sb-2D0*Pa21**2*Zh22*Zh31*cb*la5
     &      *sb+2D0*Pa21**2*Zh22*Zh32*cb*la7*sb+2D0*Pa21
     &      **2*Zh23*Zh33*cb*la4s*sb+2D0*Pa21**2*Zh23*Zh33
     &      *cb*la5s*sb+2D0*Pa21**2*Zh23*Zh33*cb*la6s*sb
      AAABR(1865) = Pa11*Pa21*Zh11*Zh21*la3+Pa11*Pa21*Zh11
     &      *Zh21*la4-Pa11*Pa21*Zh11*Zh21*la5-Pa11*Pa21*Zh11
     &      *Zh22*la7-Pa11*Pa21*Zh12*Zh21*la7+Pa11*Pa21*Zh12
     &      *Zh22*la2+Pa11*Pa21*Zh13*Zh23*la2s+Pa12*Pa22
     &      *Zh11*Zh21*la1s-Pa12*Pa22*Zh11*Zh22*la4s+Pa12
     &      *Pa22*Zh11*Zh22*la5s+Pa12*Pa22*Zh11*Zh22*la6s
     &      -Pa12*Pa22*Zh12*Zh21*la4s+Pa12*Pa22*Zh12*Zh21
     &      *la5s+Pa12*Pa22*Zh12*Zh21*la6s+Pa12*Pa22*Zh12
     &      *Zh22*la2s+2D0*Pa12*Pa22*Zh13*Zh23*la3s-12D0
     &      *Pa12*Pa22*Zh13*Zh23*la7s+Pa11*Pa21*Zh11*Zh21
     &      *la1*sb**2-Pa11*Pa21*Zh11*Zh21*la3*sb**2-Pa11
     &      *Pa21*Zh11*Zh21*la4*sb**2+Pa11*Pa21*Zh11*Zh21
     &      *la5*sb**2-Pa11*Pa21*Zh11*Zh22*la6*sb**2+Pa11
     &      *Pa21*Zh11*Zh22*la7*sb**2-Pa11*Pa21*Zh12*Zh21
     &      *la6*sb**2+Pa11*Pa21*Zh12*Zh21*la7*sb**2-Pa11
     &      *Pa21*Zh12*Zh22*la2*sb**2+Pa11*Pa21*Zh12*Zh22
     &      *la3*sb**2+Pa11*Pa21*Zh12*Zh22*la4*sb**2-Pa11
     &      *Pa21*Zh12*Zh22*la5*sb**2+Pa11*Pa21*Zh13*Zh23
     &      *la1s*sb**2-Pa11*Pa21*Zh13*Zh23*la2s*sb**2+Pa11
     &      *Pa22*Zh11*Zh23*cb*la5s-Pa11*Pa22*Zh11*Zh23*cb
     &      *la6s+Pa11*Pa22*Zh12*Zh23*la5s*sb-Pa11*Pa22*Zh12
     &      *Zh23*la6s*sb+Pa11*Pa22*Zh13*Zh21*cb*la5s-Pa11
     &      *Pa22*Zh13*Zh21*cb*la6s+Pa11*Pa22*Zh13*Zh22*la5s
     &      *sb-Pa11*Pa22*Zh13*Zh22*la6s*sb+Pa12*Pa21*Zh11
     &      *Zh23*cb*la5s-Pa12*Pa21*Zh11*Zh23*cb*la6s+Pa12
     &      *Pa21*Zh12*Zh23*la5s*sb-Pa12*Pa21*Zh12*Zh23*la6s
     &      *sb+Pa12*Pa21*Zh13*Zh21*cb*la5s-Pa12*Pa21*Zh13
     &      *Zh21*cb*la6s+Pa12*Pa21*Zh13*Zh22*la5s*sb-Pa12
     &      *Pa21*Zh13*Zh22*la6s*sb+2D0*Pa11*Pa21*Zh11*Zh21
     &      *cb*la6*sb-2D0*Pa11*Pa21*Zh11*Zh22*cb*la5*sb
     &      -2D0*Pa11*Pa21*Zh12*Zh21*cb*la5*sb+2D0*Pa11*Pa21
     &      *Zh12*Zh22*cb*la7*sb+2D0*Pa11*Pa21*Zh13*Zh23
     &      *cb*la4s*sb+2D0*Pa11*Pa21*Zh13*Zh23*cb*la5s*sb
     &      +2D0*Pa11*Pa21*Zh13*Zh23*cb*la6s*sb
      AAABR(1866) = Pa11**2*Zh31**2*la3+Pa11**2*Zh31**2*la4
     &      -Pa11**2*Zh31**2*la5+Pa11**2*Zh32**2*la2+Pa11
     &      **2*Zh33**2*la2s+Pa12**2*Zh31**2*la1s+Pa12**2
     &      *Zh32**2*la2s+2D0*Pa12**2*Zh33**2*la3s-12D0*Pa12
     &      **2*Zh33**2*la7s-2D0*Pa11**2*Zh31*Zh32*la7+Pa11
     &      **2*Zh31**2*la1*sb**2-Pa11**2*Zh31**2*la3*sb
     &      **2-Pa11**2*Zh31**2*la4*sb**2+Pa11**2*Zh31**2
     &      *la5*sb**2-Pa11**2*Zh32**2*la2*sb**2+Pa11**2
     &      *Zh32**2*la3*sb**2+Pa11**2*Zh32**2*la4*sb**2
     &      -Pa11**2*Zh32**2*la5*sb**2+Pa11**2*Zh33**2*la1s
     &      *sb**2-Pa11**2*Zh33**2*la2s*sb**2-2D0*Pa12**2
     &      *Zh31*Zh32*la4s+2D0*Pa12**2*Zh31*Zh32*la5s+2D0
     &      *Pa12**2*Zh31*Zh32*la6s-2D0*Pa11**2*Zh31*Zh32
     &      *la6*sb**2+2D0*Pa11**2*Zh31*Zh32*la7*sb**2+2D0
     &      *Pa11**2*Zh31**2*cb*la6*sb+2D0*Pa11**2*Zh32**2
     &      *cb*la7*sb+2D0*Pa11**2*Zh33**2*cb*la4s*sb+2D0
     &      *Pa11**2*Zh33**2*cb*la5s*sb+2D0*Pa11**2*Zh33
     &      **2*cb*la6s*sb+4D0*Pa11*Pa12*Zh31*Zh33*cb*la5s
     &      -4D0*Pa11*Pa12*Zh31*Zh33*cb*la6s+4D0*Pa11*Pa12
     &      *Zh32*Zh33*la5s*sb-4D0*Pa11*Pa12*Zh32*Zh33*la6s
     &      *sb-4D0*Pa11**2*Zh31*Zh32*cb*la5*sb
      AAABR(1867) = Pa11**2*Zh11*Zh21*la3+Pa11**2*Zh11*Zh21
     &      *la4-Pa11**2*Zh11*Zh21*la5-Pa11**2*Zh11*Zh22
     &      *la7-Pa11**2*Zh12*Zh21*la7+Pa11**2*Zh12*Zh22
     &      *la2+Pa11**2*Zh13*Zh23*la2s+Pa12**2*Zh11*Zh21
     &      *la1s-Pa12**2*Zh11*Zh22*la4s+Pa12**2*Zh11*Zh22
     &      *la5s+Pa12**2*Zh11*Zh22*la6s-Pa12**2*Zh12*Zh21
     &      *la4s+Pa12**2*Zh12*Zh21*la5s+Pa12**2*Zh12*Zh21
     &      *la6s+Pa12**2*Zh12*Zh22*la2s+2D0*Pa12**2*Zh13
     &      *Zh23*la3s-12D0*Pa12**2*Zh13*Zh23*la7s+Pa11**2
     &      *Zh11*Zh21*la1*sb**2-Pa11**2*Zh11*Zh21*la3*sb
     &      **2-Pa11**2*Zh11*Zh21*la4*sb**2+Pa11**2*Zh11
     &      *Zh21*la5*sb**2-Pa11**2*Zh11*Zh22*la6*sb**2+Pa11
     &      **2*Zh11*Zh22*la7*sb**2-Pa11**2*Zh12*Zh21*la6
     &      *sb**2+Pa11**2*Zh12*Zh21*la7*sb**2-Pa11**2*Zh12
     &      *Zh22*la2*sb**2+Pa11**2*Zh12*Zh22*la3*sb**2+Pa11
     &      **2*Zh12*Zh22*la4*sb**2-Pa11**2*Zh12*Zh22*la5
     &      *sb**2+Pa11**2*Zh13*Zh23*la1s*sb**2-Pa11**2*Zh13
     &      *Zh23*la2s*sb**2+2D0*Pa11*Pa12*Zh11*Zh23*cb*la5s
     &      -2D0*Pa11*Pa12*Zh11*Zh23*cb*la6s+2D0*Pa11*Pa12
     &      *Zh12*Zh23*la5s*sb-2D0*Pa11*Pa12*Zh12*Zh23*la6s
     &      *sb+2D0*Pa11*Pa12*Zh13*Zh21*cb*la5s-2D0*Pa11
     &      *Pa12*Zh13*Zh21*cb*la6s+2D0*Pa11*Pa12*Zh13*Zh22
     &      *la5s*sb-2D0*Pa11*Pa12*Zh13*Zh22*la6s*sb+2D0
     &      *Pa11**2*Zh11*Zh21*cb*la6*sb-2D0*Pa11**2*Zh11
     &      *Zh22*cb*la5*sb-2D0*Pa11**2*Zh12*Zh21*cb*la5
     &      *sb+2D0*Pa11**2*Zh12*Zh22*cb*la7*sb+2D0*Pa11
     &      **2*Zh13*Zh23*cb*la4s*sb+2D0*Pa11**2*Zh13*Zh23
     &      *cb*la5s*sb+2D0*Pa11**2*Zh13*Zh23*cb*la6s*sb
      AAABR(1868) = Pa11*Pa21*Zh11**2*la3+Pa11*Pa21*Zh11
     &      **2*la4-Pa11*Pa21*Zh11**2*la5+Pa11*Pa21*Zh12
     &      **2*la2+Pa11*Pa21*Zh13**2*la2s+Pa12*Pa22*Zh11
     &      **2*la1s+Pa12*Pa22*Zh12**2*la2s+2D0*Pa12*Pa22
     &      *Zh13**2*la3s-12D0*Pa12*Pa22*Zh13**2*la7s-2D0
     &      *Pa11*Pa21*Zh11*Zh12*la7+Pa11*Pa21*Zh11**2*la1
     &      *sb**2-Pa11*Pa21*Zh11**2*la3*sb**2-Pa11*Pa21
     &      *Zh11**2*la4*sb**2+Pa11*Pa21*Zh11**2*la5*sb**2
     &      -Pa11*Pa21*Zh12**2*la2*sb**2+Pa11*Pa21*Zh12**2
     &      *la3*sb**2+Pa11*Pa21*Zh12**2*la4*sb**2-Pa11*Pa21
     &      *Zh12**2*la5*sb**2+Pa11*Pa21*Zh13**2*la1s*sb
     &      **2-Pa11*Pa21*Zh13**2*la2s*sb**2-2D0*Pa12*Pa22
     &      *Zh11*Zh12*la4s+2D0*Pa12*Pa22*Zh11*Zh12*la5s
     &      +2D0*Pa12*Pa22*Zh11*Zh12*la6s-2D0*Pa11*Pa21*Zh11
     &      *Zh12*la6*sb**2+2D0*Pa11*Pa21*Zh11*Zh12*la7*sb
     &      **2+2D0*Pa11*Pa21*Zh11**2*cb*la6*sb+2D0*Pa11
     &      *Pa21*Zh12**2*cb*la7*sb+2D0*Pa11*Pa21*Zh13**2
     &      *cb*la4s*sb+2D0*Pa11*Pa21*Zh13**2*cb*la5s*sb
     &      +2D0*Pa11*Pa21*Zh13**2*cb*la6s*sb+2D0*Pa11*Pa22
     &      *Zh11*Zh13*cb*la5s-2D0*Pa11*Pa22*Zh11*Zh13*cb
     &      *la6s+2D0*Pa11*Pa22*Zh12*Zh13*la5s*sb-2D0*Pa11
     &      *Pa22*Zh12*Zh13*la6s*sb+2D0*Pa12*Pa21*Zh11*Zh13
     &      *cb*la5s-2D0*Pa12*Pa21*Zh11*Zh13*cb*la6s+2D0
     &      *Pa12*Pa21*Zh12*Zh13*la5s*sb-2D0*Pa12*Pa21*Zh12
     &      *Zh13*la6s*sb-4D0*Pa11*Pa21*Zh11*Zh12*cb*la5*sb
      AAABR(1869) = Pa11*Pa21*Zh21**2*la3+Pa11*Pa21*Zh21
     &      **2*la4-Pa11*Pa21*Zh21**2*la5+Pa11*Pa21*Zh22
     &      **2*la2+Pa11*Pa21*Zh23**2*la2s+Pa12*Pa22*Zh21
     &      **2*la1s+Pa12*Pa22*Zh22**2*la2s+2D0*Pa12*Pa22
     &      *Zh23**2*la3s-12D0*Pa12*Pa22*Zh23**2*la7s-2D0
     &      *Pa11*Pa21*Zh21*Zh22*la7+Pa11*Pa21*Zh21**2*la1
     &      *sb**2-Pa11*Pa21*Zh21**2*la3*sb**2-Pa11*Pa21
     &      *Zh21**2*la4*sb**2+Pa11*Pa21*Zh21**2*la5*sb**2
     &      -Pa11*Pa21*Zh22**2*la2*sb**2+Pa11*Pa21*Zh22**2
     &      *la3*sb**2+Pa11*Pa21*Zh22**2*la4*sb**2-Pa11*Pa21
     &      *Zh22**2*la5*sb**2+Pa11*Pa21*Zh23**2*la1s*sb
     &      **2-Pa11*Pa21*Zh23**2*la2s*sb**2-2D0*Pa12*Pa22
     &      *Zh21*Zh22*la4s+2D0*Pa12*Pa22*Zh21*Zh22*la5s
     &      +2D0*Pa12*Pa22*Zh21*Zh22*la6s-2D0*Pa11*Pa21*Zh21
     &      *Zh22*la6*sb**2+2D0*Pa11*Pa21*Zh21*Zh22*la7*sb
     &      **2+2D0*Pa11*Pa21*Zh21**2*cb*la6*sb+2D0*Pa11
     &      *Pa21*Zh22**2*cb*la7*sb+2D0*Pa11*Pa21*Zh23**2
     &      *cb*la4s*sb+2D0*Pa11*Pa21*Zh23**2*cb*la5s*sb
     &      +2D0*Pa11*Pa21*Zh23**2*cb*la6s*sb+2D0*Pa11*Pa22
     &      *Zh21*Zh23*cb*la5s-2D0*Pa11*Pa22*Zh21*Zh23*cb
     &      *la6s+2D0*Pa11*Pa22*Zh22*Zh23*la5s*sb-2D0*Pa11
     &      *Pa22*Zh22*Zh23*la6s*sb+2D0*Pa12*Pa21*Zh21*Zh23
     &      *cb*la5s-2D0*Pa12*Pa21*Zh21*Zh23*cb*la6s+2D0
     &      *Pa12*Pa21*Zh22*Zh23*la5s*sb-2D0*Pa12*Pa21*Zh22
     &      *Zh23*la6s*sb-4D0*Pa11*Pa21*Zh21*Zh22*cb*la5*sb
      AAABR(1870) = Pa21**2*Zh11*Zh31*la3+Pa21**2*Zh11*Zh31
     &      *la4-Pa21**2*Zh11*Zh31*la5-Pa21**2*Zh11*Zh32
     &      *la7-Pa21**2*Zh12*Zh31*la7+Pa21**2*Zh12*Zh32
     &      *la2+Pa21**2*Zh13*Zh33*la2s+Pa22**2*Zh11*Zh31
     &      *la1s-Pa22**2*Zh11*Zh32*la4s+Pa22**2*Zh11*Zh32
     &      *la5s+Pa22**2*Zh11*Zh32*la6s-Pa22**2*Zh12*Zh31
     &      *la4s+Pa22**2*Zh12*Zh31*la5s+Pa22**2*Zh12*Zh31
     &      *la6s+Pa22**2*Zh12*Zh32*la2s+2D0*Pa22**2*Zh13
     &      *Zh33*la3s-12D0*Pa22**2*Zh13*Zh33*la7s+Pa21**2
     &      *Zh11*Zh31*la1*sb**2-Pa21**2*Zh11*Zh31*la3*sb
     &      **2-Pa21**2*Zh11*Zh31*la4*sb**2+Pa21**2*Zh11
     &      *Zh31*la5*sb**2-Pa21**2*Zh11*Zh32*la6*sb**2+Pa21
     &      **2*Zh11*Zh32*la7*sb**2-Pa21**2*Zh12*Zh31*la6
     &      *sb**2+Pa21**2*Zh12*Zh31*la7*sb**2-Pa21**2*Zh12
     &      *Zh32*la2*sb**2+Pa21**2*Zh12*Zh32*la3*sb**2+Pa21
     &      **2*Zh12*Zh32*la4*sb**2-Pa21**2*Zh12*Zh32*la5
     &      *sb**2+Pa21**2*Zh13*Zh33*la1s*sb**2-Pa21**2*Zh13
     &      *Zh33*la2s*sb**2+2D0*Pa21*Pa22*Zh11*Zh33*cb*la5s
     &      -2D0*Pa21*Pa22*Zh11*Zh33*cb*la6s+2D0*Pa21*Pa22
     &      *Zh12*Zh33*la5s*sb-2D0*Pa21*Pa22*Zh12*Zh33*la6s
     &      *sb+2D0*Pa21*Pa22*Zh13*Zh31*cb*la5s-2D0*Pa21
     &      *Pa22*Zh13*Zh31*cb*la6s+2D0*Pa21*Pa22*Zh13*Zh32
     &      *la5s*sb-2D0*Pa21*Pa22*Zh13*Zh32*la6s*sb+2D0
     &      *Pa21**2*Zh11*Zh31*cb*la6*sb-2D0*Pa21**2*Zh11
     &      *Zh32*cb*la5*sb-2D0*Pa21**2*Zh12*Zh31*cb*la5
     &      *sb+2D0*Pa21**2*Zh12*Zh32*cb*la7*sb+2D0*Pa21
     &      **2*Zh13*Zh33*cb*la4s*sb+2D0*Pa21**2*Zh13*Zh33
     &      *cb*la5s*sb+2D0*Pa21**2*Zh13*Zh33*cb*la6s*sb
      AAABR(1871) = Pa11**2*Zh11**2*la3+Pa11**2*Zh11**2*la4
     &      -Pa11**2*Zh11**2*la5+Pa11**2*Zh12**2*la2+Pa11
     &      **2*Zh13**2*la2s+Pa12**2*Zh11**2*la1s+Pa12**2
     &      *Zh12**2*la2s+2D0*Pa12**2*Zh13**2*la3s-12D0*Pa12
     &      **2*Zh13**2*la7s-2D0*Pa11**2*Zh11*Zh12*la7+Pa11
     &      **2*Zh11**2*la1*sb**2-Pa11**2*Zh11**2*la3*sb
     &      **2-Pa11**2*Zh11**2*la4*sb**2+Pa11**2*Zh11**2
     &      *la5*sb**2-Pa11**2*Zh12**2*la2*sb**2+Pa11**2
     &      *Zh12**2*la3*sb**2+Pa11**2*Zh12**2*la4*sb**2
     &      -Pa11**2*Zh12**2*la5*sb**2+Pa11**2*Zh13**2*la1s
     &      *sb**2-Pa11**2*Zh13**2*la2s*sb**2-2D0*Pa12**2
     &      *Zh11*Zh12*la4s+2D0*Pa12**2*Zh11*Zh12*la5s+2D0
     &      *Pa12**2*Zh11*Zh12*la6s-2D0*Pa11**2*Zh11*Zh12
     &      *la6*sb**2+2D0*Pa11**2*Zh11*Zh12*la7*sb**2+2D0
     &      *Pa11**2*Zh11**2*cb*la6*sb+2D0*Pa11**2*Zh12**2
     &      *cb*la7*sb+2D0*Pa11**2*Zh13**2*cb*la4s*sb+2D0
     &      *Pa11**2*Zh13**2*cb*la5s*sb+2D0*Pa11**2*Zh13
     &      **2*cb*la6s*sb+4D0*Pa11*Pa12*Zh11*Zh13*cb*la5s
     &      -4D0*Pa11*Pa12*Zh11*Zh13*cb*la6s+4D0*Pa11*Pa12
     &      *Zh12*Zh13*la5s*sb-4D0*Pa11*Pa12*Zh12*Zh13*la6s
     &      *sb-4D0*Pa11**2*Zh11*Zh12*cb*la5*sb
      AAABR(1872) = Pa21**2*Zh31**2*la3+Pa21**2*Zh31**2*la4
     &      -Pa21**2*Zh31**2*la5+Pa21**2*Zh32**2*la2+Pa21
     &      **2*Zh33**2*la2s+Pa22**2*Zh31**2*la1s+Pa22**2
     &      *Zh32**2*la2s+2D0*Pa22**2*Zh33**2*la3s-12D0*Pa22
     &      **2*Zh33**2*la7s-2D0*Pa21**2*Zh31*Zh32*la7+Pa21
     &      **2*Zh31**2*la1*sb**2-Pa21**2*Zh31**2*la3*sb
     &      **2-Pa21**2*Zh31**2*la4*sb**2+Pa21**2*Zh31**2
     &      *la5*sb**2-Pa21**2*Zh32**2*la2*sb**2+Pa21**2
     &      *Zh32**2*la3*sb**2+Pa21**2*Zh32**2*la4*sb**2
     &      -Pa21**2*Zh32**2*la5*sb**2+Pa21**2*Zh33**2*la1s
     &      *sb**2-Pa21**2*Zh33**2*la2s*sb**2-2D0*Pa22**2
     &      *Zh31*Zh32*la4s+2D0*Pa22**2*Zh31*Zh32*la5s+2D0
     &      *Pa22**2*Zh31*Zh32*la6s-2D0*Pa21**2*Zh31*Zh32
     &      *la6*sb**2+2D0*Pa21**2*Zh31*Zh32*la7*sb**2+2D0
     &      *Pa21**2*Zh31**2*cb*la6*sb+2D0*Pa21**2*Zh32**2
     &      *cb*la7*sb+2D0*Pa21**2*Zh33**2*cb*la4s*sb+2D0
     &      *Pa21**2*Zh33**2*cb*la5s*sb+2D0*Pa21**2*Zh33
     &      **2*cb*la6s*sb+4D0*Pa21*Pa22*Zh31*Zh33*cb*la5s
     &      -4D0*Pa21*Pa22*Zh31*Zh33*cb*la6s+4D0*Pa21*Pa22
     &      *Zh32*Zh33*la5s*sb-4D0*Pa21*Pa22*Zh32*Zh33*la6s
     &      *sb-4D0*Pa21**2*Zh31*Zh32*cb*la5*sb
      AAABR(1873) = Pa11*Pa21*Zh21*Zh31*la3+Pa11*Pa21*Zh21
     &      *Zh31*la4-Pa11*Pa21*Zh21*Zh31*la5-Pa11*Pa21*Zh21
     &      *Zh32*la7-Pa11*Pa21*Zh22*Zh31*la7+Pa11*Pa21*Zh22
     &      *Zh32*la2+Pa11*Pa21*Zh23*Zh33*la2s+Pa12*Pa22
     &      *Zh21*Zh31*la1s-Pa12*Pa22*Zh21*Zh32*la4s+Pa12
     &      *Pa22*Zh21*Zh32*la5s+Pa12*Pa22*Zh21*Zh32*la6s
     &      -Pa12*Pa22*Zh22*Zh31*la4s+Pa12*Pa22*Zh22*Zh31
     &      *la5s+Pa12*Pa22*Zh22*Zh31*la6s+Pa12*Pa22*Zh22
     &      *Zh32*la2s+2D0*Pa12*Pa22*Zh23*Zh33*la3s-12D0
     &      *Pa12*Pa22*Zh23*Zh33*la7s+Pa11*Pa21*Zh21*Zh31
     &      *la1*sb**2-Pa11*Pa21*Zh21*Zh31*la3*sb**2-Pa11
     &      *Pa21*Zh21*Zh31*la4*sb**2+Pa11*Pa21*Zh21*Zh31
     &      *la5*sb**2-Pa11*Pa21*Zh21*Zh32*la6*sb**2+Pa11
     &      *Pa21*Zh21*Zh32*la7*sb**2-Pa11*Pa21*Zh22*Zh31
     &      *la6*sb**2+Pa11*Pa21*Zh22*Zh31*la7*sb**2-Pa11
     &      *Pa21*Zh22*Zh32*la2*sb**2+Pa11*Pa21*Zh22*Zh32
     &      *la3*sb**2+Pa11*Pa21*Zh22*Zh32*la4*sb**2-Pa11
     &      *Pa21*Zh22*Zh32*la5*sb**2+Pa11*Pa21*Zh23*Zh33
     &      *la1s*sb**2-Pa11*Pa21*Zh23*Zh33*la2s*sb**2+Pa11
     &      *Pa22*Zh21*Zh33*cb*la5s-Pa11*Pa22*Zh21*Zh33*cb
     &      *la6s+Pa11*Pa22*Zh22*Zh33*la5s*sb-Pa11*Pa22*Zh22
     &      *Zh33*la6s*sb+Pa11*Pa22*Zh23*Zh31*cb*la5s-Pa11
     &      *Pa22*Zh23*Zh31*cb*la6s+Pa11*Pa22*Zh23*Zh32*la5s
     &      *sb-Pa11*Pa22*Zh23*Zh32*la6s*sb+Pa12*Pa21*Zh21
     &      *Zh33*cb*la5s-Pa12*Pa21*Zh21*Zh33*cb*la6s+Pa12
     &      *Pa21*Zh22*Zh33*la5s*sb-Pa12*Pa21*Zh22*Zh33*la6s
     &      *sb+Pa12*Pa21*Zh23*Zh31*cb*la5s-Pa12*Pa21*Zh23
     &      *Zh31*cb*la6s+Pa12*Pa21*Zh23*Zh32*la5s*sb-Pa12
     &      *Pa21*Zh23*Zh32*la6s*sb+2D0*Pa11*Pa21*Zh21*Zh31
     &      *cb*la6*sb-2D0*Pa11*Pa21*Zh21*Zh32*cb*la5*sb
     &      -2D0*Pa11*Pa21*Zh22*Zh31*cb*la5*sb+2D0*Pa11*Pa21
     &      *Zh22*Zh32*cb*la7*sb+2D0*Pa11*Pa21*Zh23*Zh33
     &      *cb*la4s*sb+2D0*Pa11*Pa21*Zh23*Zh33*cb*la5s*sb
     &      +2D0*Pa11*Pa21*Zh23*Zh33*cb*la6s*sb
      AAABR(1874) = Pa21**2*Zh11**2*la3+Pa21**2*Zh11**2*la4
     &      -Pa21**2*Zh11**2*la5+Pa21**2*Zh12**2*la2+Pa21
     &      **2*Zh13**2*la2s+Pa22**2*Zh11**2*la1s+Pa22**2
     &      *Zh12**2*la2s+2D0*Pa22**2*Zh13**2*la3s-12D0*Pa22
     &      **2*Zh13**2*la7s-2D0*Pa21**2*Zh11*Zh12*la7+Pa21
     &      **2*Zh11**2*la1*sb**2-Pa21**2*Zh11**2*la3*sb
     &      **2-Pa21**2*Zh11**2*la4*sb**2+Pa21**2*Zh11**2
     &      *la5*sb**2-Pa21**2*Zh12**2*la2*sb**2+Pa21**2
     &      *Zh12**2*la3*sb**2+Pa21**2*Zh12**2*la4*sb**2
     &      -Pa21**2*Zh12**2*la5*sb**2+Pa21**2*Zh13**2*la1s
     &      *sb**2-Pa21**2*Zh13**2*la2s*sb**2-2D0*Pa22**2
     &      *Zh11*Zh12*la4s+2D0*Pa22**2*Zh11*Zh12*la5s+2D0
     &      *Pa22**2*Zh11*Zh12*la6s-2D0*Pa21**2*Zh11*Zh12
     &      *la6*sb**2+2D0*Pa21**2*Zh11*Zh12*la7*sb**2+2D0
     &      *Pa21**2*Zh11**2*cb*la6*sb+2D0*Pa21**2*Zh12**2
     &      *cb*la7*sb+2D0*Pa21**2*Zh13**2*cb*la4s*sb+2D0
     &      *Pa21**2*Zh13**2*cb*la5s*sb+2D0*Pa21**2*Zh13
     &      **2*cb*la6s*sb+4D0*Pa21*Pa22*Zh11*Zh13*cb*la5s
     &      -4D0*Pa21*Pa22*Zh11*Zh13*cb*la6s+4D0*Pa21*Pa22
     &      *Zh12*Zh13*la5s*sb-4D0*Pa21*Pa22*Zh12*Zh13*la6s
     &      *sb-4D0*Pa21**2*Zh11*Zh12*cb*la5*sb
      AAABR(1875) = Pa21**2*Zh11*Zh21*la3+Pa21**2*Zh11*Zh21
     &      *la4-Pa21**2*Zh11*Zh21*la5-Pa21**2*Zh11*Zh22
     &      *la7-Pa21**2*Zh12*Zh21*la7+Pa21**2*Zh12*Zh22
     &      *la2+Pa21**2*Zh13*Zh23*la2s+Pa22**2*Zh11*Zh21
     &      *la1s-Pa22**2*Zh11*Zh22*la4s+Pa22**2*Zh11*Zh22
     &      *la5s+Pa22**2*Zh11*Zh22*la6s-Pa22**2*Zh12*Zh21
     &      *la4s+Pa22**2*Zh12*Zh21*la5s+Pa22**2*Zh12*Zh21
     &      *la6s+Pa22**2*Zh12*Zh22*la2s+2D0*Pa22**2*Zh13
     &      *Zh23*la3s-12D0*Pa22**2*Zh13*Zh23*la7s+Pa21**2
     &      *Zh11*Zh21*la1*sb**2-Pa21**2*Zh11*Zh21*la3*sb
     &      **2-Pa21**2*Zh11*Zh21*la4*sb**2+Pa21**2*Zh11
     &      *Zh21*la5*sb**2-Pa21**2*Zh11*Zh22*la6*sb**2+Pa21
     &      **2*Zh11*Zh22*la7*sb**2-Pa21**2*Zh12*Zh21*la6
     &      *sb**2+Pa21**2*Zh12*Zh21*la7*sb**2-Pa21**2*Zh12
     &      *Zh22*la2*sb**2+Pa21**2*Zh12*Zh22*la3*sb**2+Pa21
     &      **2*Zh12*Zh22*la4*sb**2-Pa21**2*Zh12*Zh22*la5
     &      *sb**2+Pa21**2*Zh13*Zh23*la1s*sb**2-Pa21**2*Zh13
     &      *Zh23*la2s*sb**2+2D0*Pa21*Pa22*Zh11*Zh23*cb*la5s
     &      -2D0*Pa21*Pa22*Zh11*Zh23*cb*la6s+2D0*Pa21*Pa22
     &      *Zh12*Zh23*la5s*sb-2D0*Pa21*Pa22*Zh12*Zh23*la6s
     &      *sb+2D0*Pa21*Pa22*Zh13*Zh21*cb*la5s-2D0*Pa21
     &      *Pa22*Zh13*Zh21*cb*la6s+2D0*Pa21*Pa22*Zh13*Zh22
     &      *la5s*sb-2D0*Pa21*Pa22*Zh13*Zh22*la6s*sb+2D0
     &      *Pa21**2*Zh11*Zh21*cb*la6*sb-2D0*Pa21**2*Zh11
     &      *Zh22*cb*la5*sb-2D0*Pa21**2*Zh12*Zh21*cb*la5
     &      *sb+2D0*Pa21**2*Zh12*Zh22*cb*la7*sb+2D0*Pa21
     &      **2*Zh13*Zh23*cb*la4s*sb+2D0*Pa21**2*Zh13*Zh23
     &      *cb*la5s*sb+2D0*Pa21**2*Zh13*Zh23*cb*la6s*sb
      AAABR(1876) = Pa21**2*Zh21**2*la3+Pa21**2*Zh21**2*la4
     &      -Pa21**2*Zh21**2*la5+Pa21**2*Zh22**2*la2+Pa21
     &      **2*Zh23**2*la2s+Pa22**2*Zh21**2*la1s+Pa22**2
     &      *Zh22**2*la2s+2D0*Pa22**2*Zh23**2*la3s-12D0*Pa22
     &      **2*Zh23**2*la7s-2D0*Pa21**2*Zh21*Zh22*la7+Pa21
     &      **2*Zh21**2*la1*sb**2-Pa21**2*Zh21**2*la3*sb
     &      **2-Pa21**2*Zh21**2*la4*sb**2+Pa21**2*Zh21**2
     &      *la5*sb**2-Pa21**2*Zh22**2*la2*sb**2+Pa21**2
     &      *Zh22**2*la3*sb**2+Pa21**2*Zh22**2*la4*sb**2
     &      -Pa21**2*Zh22**2*la5*sb**2+Pa21**2*Zh23**2*la1s
     &      *sb**2-Pa21**2*Zh23**2*la2s*sb**2-2D0*Pa22**2
     &      *Zh21*Zh22*la4s+2D0*Pa22**2*Zh21*Zh22*la5s+2D0
     &      *Pa22**2*Zh21*Zh22*la6s-2D0*Pa21**2*Zh21*Zh22
     &      *la6*sb**2+2D0*Pa21**2*Zh21*Zh22*la7*sb**2+2D0
     &      *Pa21**2*Zh21**2*cb*la6*sb+2D0*Pa21**2*Zh22**2
     &      *cb*la7*sb+2D0*Pa21**2*Zh23**2*cb*la4s*sb+2D0
     &      *Pa21**2*Zh23**2*cb*la5s*sb+2D0*Pa21**2*Zh23
     &      **2*cb*la6s*sb+4D0*Pa21*Pa22*Zh21*Zh23*cb*la5s
     &      -4D0*Pa21*Pa22*Zh21*Zh23*cb*la6s+4D0*Pa21*Pa22
     &      *Zh22*Zh23*la5s*sb-4D0*Pa21*Pa22*Zh22*Zh23*la6s
     &      *sb-4D0*Pa21**2*Zh21*Zh22*cb*la5*sb
      AAABR(1877) = Pa11**2*Zh21*Zh31*la3+Pa11**2*Zh21*Zh31
     &      *la4-Pa11**2*Zh21*Zh31*la5-Pa11**2*Zh21*Zh32
     &      *la7-Pa11**2*Zh22*Zh31*la7+Pa11**2*Zh22*Zh32
     &      *la2+Pa11**2*Zh23*Zh33*la2s+Pa12**2*Zh21*Zh31
     &      *la1s-Pa12**2*Zh21*Zh32*la4s+Pa12**2*Zh21*Zh32
     &      *la5s+Pa12**2*Zh21*Zh32*la6s-Pa12**2*Zh22*Zh31
     &      *la4s+Pa12**2*Zh22*Zh31*la5s+Pa12**2*Zh22*Zh31
     &      *la6s+Pa12**2*Zh22*Zh32*la2s+2D0*Pa12**2*Zh23
     &      *Zh33*la3s-12D0*Pa12**2*Zh23*Zh33*la7s+Pa11**2
     &      *Zh21*Zh31*la1*sb**2-Pa11**2*Zh21*Zh31*la3*sb
     &      **2-Pa11**2*Zh21*Zh31*la4*sb**2+Pa11**2*Zh21
     &      *Zh31*la5*sb**2-Pa11**2*Zh21*Zh32*la6*sb**2+Pa11
     &      **2*Zh21*Zh32*la7*sb**2-Pa11**2*Zh22*Zh31*la6
     &      *sb**2+Pa11**2*Zh22*Zh31*la7*sb**2-Pa11**2*Zh22
     &      *Zh32*la2*sb**2+Pa11**2*Zh22*Zh32*la3*sb**2+Pa11
     &      **2*Zh22*Zh32*la4*sb**2-Pa11**2*Zh22*Zh32*la5
     &      *sb**2+Pa11**2*Zh23*Zh33*la1s*sb**2-Pa11**2*Zh23
     &      *Zh33*la2s*sb**2+2D0*Pa11*Pa12*Zh21*Zh33*cb*la5s
     &      -2D0*Pa11*Pa12*Zh21*Zh33*cb*la6s+2D0*Pa11*Pa12
     &      *Zh22*Zh33*la5s*sb-2D0*Pa11*Pa12*Zh22*Zh33*la6s
     &      *sb+2D0*Pa11*Pa12*Zh23*Zh31*cb*la5s-2D0*Pa11
     &      *Pa12*Zh23*Zh31*cb*la6s+2D0*Pa11*Pa12*Zh23*Zh32
     &      *la5s*sb-2D0*Pa11*Pa12*Zh23*Zh32*la6s*sb+2D0
     &      *Pa11**2*Zh21*Zh31*cb*la6*sb-2D0*Pa11**2*Zh21
     &      *Zh32*cb*la5*sb-2D0*Pa11**2*Zh22*Zh31*cb*la5
     &      *sb+2D0*Pa11**2*Zh22*Zh32*cb*la7*sb+2D0*Pa11
     &      **2*Zh23*Zh33*cb*la4s*sb+2D0*Pa11**2*Zh23*Zh33
     &      *cb*la5s*sb+2D0*Pa11**2*Zh23*Zh33*cb*la6s*sb
      AAABR(1878) = Pa11**2*Zh21**2*la3+Pa11**2*Zh21**2*la4
     &      -Pa11**2*Zh21**2*la5+Pa11**2*Zh22**2*la2+Pa11
     &      **2*Zh23**2*la2s+Pa12**2*Zh21**2*la1s+Pa12**2
     &      *Zh22**2*la2s+2D0*Pa12**2*Zh23**2*la3s-12D0*Pa12
     &      **2*Zh23**2*la7s-2D0*Pa11**2*Zh21*Zh22*la7+Pa11
     &      **2*Zh21**2*la1*sb**2-Pa11**2*Zh21**2*la3*sb
     &      **2-Pa11**2*Zh21**2*la4*sb**2+Pa11**2*Zh21**2
     &      *la5*sb**2-Pa11**2*Zh22**2*la2*sb**2+Pa11**2
     &      *Zh22**2*la3*sb**2+Pa11**2*Zh22**2*la4*sb**2
     &      -Pa11**2*Zh22**2*la5*sb**2+Pa11**2*Zh23**2*la1s
     &      *sb**2-Pa11**2*Zh23**2*la2s*sb**2-2D0*Pa12**2
     &      *Zh21*Zh22*la4s+2D0*Pa12**2*Zh21*Zh22*la5s+2D0
     &      *Pa12**2*Zh21*Zh22*la6s-2D0*Pa11**2*Zh21*Zh22
     &      *la6*sb**2+2D0*Pa11**2*Zh21*Zh22*la7*sb**2+2D0
     &      *Pa11**2*Zh21**2*cb*la6*sb+2D0*Pa11**2*Zh22**2
     &      *cb*la7*sb+2D0*Pa11**2*Zh23**2*cb*la4s*sb+2D0
     &      *Pa11**2*Zh23**2*cb*la5s*sb+2D0*Pa11**2*Zh23
     &      **2*cb*la6s*sb+4D0*Pa11*Pa12*Zh21*Zh23*cb*la5s
     &      -4D0*Pa11*Pa12*Zh21*Zh23*cb*la6s+4D0*Pa11*Pa12
     &      *Zh22*Zh23*la5s*sb-4D0*Pa11*Pa12*Zh22*Zh23*la6s
     &      *sb-4D0*Pa11**2*Zh21*Zh22*cb*la5*sb
      AAABR(1879) = Pa11**2*Zh11*Zh31*la3+Pa11**2*Zh11*Zh31
     &      *la4-Pa11**2*Zh11*Zh31*la5-Pa11**2*Zh11*Zh32
     &      *la7-Pa11**2*Zh12*Zh31*la7+Pa11**2*Zh12*Zh32
     &      *la2+Pa11**2*Zh13*Zh33*la2s+Pa12**2*Zh11*Zh31
     &      *la1s-Pa12**2*Zh11*Zh32*la4s+Pa12**2*Zh11*Zh32
     &      *la5s+Pa12**2*Zh11*Zh32*la6s-Pa12**2*Zh12*Zh31
     &      *la4s+Pa12**2*Zh12*Zh31*la5s+Pa12**2*Zh12*Zh31
     &      *la6s+Pa12**2*Zh12*Zh32*la2s+2D0*Pa12**2*Zh13
     &      *Zh33*la3s-12D0*Pa12**2*Zh13*Zh33*la7s+Pa11**2
     &      *Zh11*Zh31*la1*sb**2-Pa11**2*Zh11*Zh31*la3*sb
     &      **2-Pa11**2*Zh11*Zh31*la4*sb**2+Pa11**2*Zh11
     &      *Zh31*la5*sb**2-Pa11**2*Zh11*Zh32*la6*sb**2+Pa11
     &      **2*Zh11*Zh32*la7*sb**2-Pa11**2*Zh12*Zh31*la6
     &      *sb**2+Pa11**2*Zh12*Zh31*la7*sb**2-Pa11**2*Zh12
     &      *Zh32*la2*sb**2+Pa11**2*Zh12*Zh32*la3*sb**2+Pa11
     &      **2*Zh12*Zh32*la4*sb**2-Pa11**2*Zh12*Zh32*la5
     &      *sb**2+Pa11**2*Zh13*Zh33*la1s*sb**2-Pa11**2*Zh13
     &      *Zh33*la2s*sb**2+2D0*Pa11*Pa12*Zh11*Zh33*cb*la5s
     &      -2D0*Pa11*Pa12*Zh11*Zh33*cb*la6s+2D0*Pa11*Pa12
     &      *Zh12*Zh33*la5s*sb-2D0*Pa11*Pa12*Zh12*Zh33*la6s
     &      *sb+2D0*Pa11*Pa12*Zh13*Zh31*cb*la5s-2D0*Pa11
     &      *Pa12*Zh13*Zh31*cb*la6s+2D0*Pa11*Pa12*Zh13*Zh32
     &      *la5s*sb-2D0*Pa11*Pa12*Zh13*Zh32*la6s*sb+2D0
     &      *Pa11**2*Zh11*Zh31*cb*la6*sb-2D0*Pa11**2*Zh11
     &      *Zh32*cb*la5*sb-2D0*Pa11**2*Zh12*Zh31*cb*la5
     &      *sb+2D0*Pa11**2*Zh12*Zh32*cb*la7*sb+2D0*Pa11
     &      **2*Zh13*Zh33*cb*la4s*sb+2D0*Pa11**2*Zh13*Zh33
     &      *cb*la5s*sb+2D0*Pa11**2*Zh13*Zh33*cb*la6s*sb
      AAABR(1880) = Zh11*Zh31**3*la1-Zh11*Zh32**3*la7-Zh12
     &      *Zh31**3*la6+Zh12*Zh32**3*la2+2D0*Zh13*Zh33**3
     &      *la3s+4D0*Zh13*Zh33**3*la7s+4D0*Zh13*Zh33**3
     &      *la8s+Zh11*Zh31*Zh32**2*la3+Zh11*Zh31*Zh32**2
     &      *la4+Zh11*Zh31*Zh32**2*la5+Zh11*Zh31*Zh33**2
     &      *la1s-3D0*Zh11*Zh31**2*Zh32*la6-Zh11*Zh32*Zh33
     &      **2*la4s-Zh11*Zh32*Zh33**2*la5s-Zh11*Zh32*Zh33
     &      **2*la6s-3D0*Zh12*Zh31*Zh32**2*la7-Zh12*Zh31
     &      *Zh33**2*la4s-Zh12*Zh31*Zh33**2*la5s-Zh12*Zh31
     &      *Zh33**2*la6s+Zh12*Zh31**2*Zh32*la3+Zh12*Zh31
     &      **2*Zh32*la4+Zh12*Zh31**2*Zh32*la5+Zh12*Zh32
     &      *Zh33**2*la2s+Zh13*Zh31**2*Zh33*la1s+Zh13*Zh32
     &      **2*Zh33*la2s-2D0*Zh13*Zh31*Zh32*Zh33*la4s-2D0
     &      *Zh13*Zh31*Zh32*Zh33*la5s-2D0*Zh13*Zh31*Zh32
     &      *Zh33*la6s
      AAABR(1881) = Zh11**3*Zh31*la1-Zh11**3*Zh32*la6-Zh12
     &      **3*Zh31*la7+Zh12**3*Zh32*la2+2D0*Zh13**3*Zh33
     &      *la3s+4D0*Zh13**3*Zh33*la7s+4D0*Zh13**3*Zh33
     &      *la8s+Zh11*Zh12**2*Zh31*la3+Zh11*Zh12**2*Zh31
     &      *la4+Zh11*Zh12**2*Zh31*la5-3D0*Zh11*Zh12**2*Zh32
     &      *la7+Zh11*Zh13**2*Zh31*la1s-Zh11*Zh13**2*Zh32
     &      *la4s-Zh11*Zh13**2*Zh32*la5s-Zh11*Zh13**2*Zh32
     &      *la6s-3D0*Zh11**2*Zh12*Zh31*la6+Zh11**2*Zh12
     &      *Zh32*la3+Zh11**2*Zh12*Zh32*la4+Zh11**2*Zh12
     &      *Zh32*la5+Zh11**2*Zh13*Zh33*la1s-Zh12*Zh13**2
     &      *Zh31*la4s-Zh12*Zh13**2*Zh31*la5s-Zh12*Zh13**2
     &      *Zh31*la6s+Zh12*Zh13**2*Zh32*la2s+Zh12**2*Zh13
     &      *Zh33*la2s-2D0*Zh11*Zh12*Zh13*Zh33*la4s-2D0*Zh11
     &      *Zh12*Zh13*Zh33*la5s-2D0*Zh11*Zh12*Zh13*Zh33
     &      *la6s
      AAABR(1882) = Zh21**4*la1+Zh22**4*la2+2D0*Zh23**4*la3s
     &      +4D0*Zh23**4*la7s+4D0*Zh23**4*la8s-4D0*Zh21*Zh22
     &      **3*la7+2D0*Zh21**2*Zh22**2*la3+2D0*Zh21**2*Zh22
     &      **2*la4+2D0*Zh21**2*Zh22**2*la5+2D0*Zh21**2*Zh23
     &      **2*la1s-4D0*Zh21**3*Zh22*la6+2D0*Zh22**2*Zh23
     &      **2*la2s-4D0*Zh21*Zh22*Zh23**2*la4s-4D0*Zh21
     &      *Zh22*Zh23**2*la5s-4D0*Zh21*Zh22*Zh23**2*la6s
      AAABR(1883) = Zh11**4*la1+Zh12**4*la2+2D0*Zh13**4*la3s
     &      +4D0*Zh13**4*la7s+4D0*Zh13**4*la8s-4D0*Zh11*Zh12
     &      **3*la7+2D0*Zh11**2*Zh12**2*la3+2D0*Zh11**2*Zh12
     &      **2*la4+2D0*Zh11**2*Zh12**2*la5+2D0*Zh11**2*Zh13
     &      **2*la1s-4D0*Zh11**3*Zh12*la6+2D0*Zh12**2*Zh13
     &      **2*la2s-4D0*Zh11*Zh12*Zh13**2*la4s-4D0*Zh11
     &      *Zh12*Zh13**2*la5s-4D0*Zh11*Zh12*Zh13**2*la6s
      AAABR(1884) = Zh31**4*la1+Zh32**4*la2+2D0*Zh33**4*la3s
     &      +4D0*Zh33**4*la7s+4D0*Zh33**4*la8s-4D0*Zh31*Zh32
     &      **3*la7+2D0*Zh31**2*Zh32**2*la3+2D0*Zh31**2*Zh32
     &      **2*la4+2D0*Zh31**2*Zh32**2*la5+2D0*Zh31**2*Zh33
     &      **2*la1s-4D0*Zh31**3*Zh32*la6+2D0*Zh32**2*Zh33
     &      **2*la2s-4D0*Zh31*Zh32*Zh33**2*la4s-4D0*Zh31
     &      *Zh32*Zh33**2*la5s-4D0*Zh31*Zh32*Zh33**2*la6s
      AAABR(1885) = Zh11**3*Zh21*la1-Zh11**3*Zh22*la6-Zh12
     &      **3*Zh21*la7+Zh12**3*Zh22*la2+2D0*Zh13**3*Zh23
     &      *la3s+4D0*Zh13**3*Zh23*la7s+4D0*Zh13**3*Zh23
     &      *la8s+Zh11*Zh12**2*Zh21*la3+Zh11*Zh12**2*Zh21
     &      *la4+Zh11*Zh12**2*Zh21*la5-3D0*Zh11*Zh12**2*Zh22
     &      *la7+Zh11*Zh13**2*Zh21*la1s-Zh11*Zh13**2*Zh22
     &      *la4s-Zh11*Zh13**2*Zh22*la5s-Zh11*Zh13**2*Zh22
     &      *la6s-3D0*Zh11**2*Zh12*Zh21*la6+Zh11**2*Zh12
     &      *Zh22*la3+Zh11**2*Zh12*Zh22*la4+Zh11**2*Zh12
     &      *Zh22*la5+Zh11**2*Zh13*Zh23*la1s-Zh12*Zh13**2
     &      *Zh21*la4s-Zh12*Zh13**2*Zh21*la5s-Zh12*Zh13**2
     &      *Zh21*la6s+Zh12*Zh13**2*Zh22*la2s+Zh12**2*Zh13
     &      *Zh23*la2s-2D0*Zh11*Zh12*Zh13*Zh23*la4s-2D0*Zh11
     &      *Zh12*Zh13*Zh23*la5s-2D0*Zh11*Zh12*Zh13*Zh23
     &      *la6s
      AAABR(1886) = Zh21*Zh31**3*la1-Zh21*Zh32**3*la7-Zh22
     &      *Zh31**3*la6+Zh22*Zh32**3*la2+2D0*Zh23*Zh33**3
     &      *la3s+4D0*Zh23*Zh33**3*la7s+4D0*Zh23*Zh33**3
     &      *la8s+Zh21*Zh31*Zh32**2*la3+Zh21*Zh31*Zh32**2
     &      *la4+Zh21*Zh31*Zh32**2*la5+Zh21*Zh31*Zh33**2
     &      *la1s-3D0*Zh21*Zh31**2*Zh32*la6-Zh21*Zh32*Zh33
     &      **2*la4s-Zh21*Zh32*Zh33**2*la5s-Zh21*Zh32*Zh33
     &      **2*la6s-3D0*Zh22*Zh31*Zh32**2*la7-Zh22*Zh31
     &      *Zh33**2*la4s-Zh22*Zh31*Zh33**2*la5s-Zh22*Zh31
     &      *Zh33**2*la6s+Zh22*Zh31**2*Zh32*la3+Zh22*Zh31
     &      **2*Zh32*la4+Zh22*Zh31**2*Zh32*la5+Zh22*Zh32
     &      *Zh33**2*la2s+Zh23*Zh31**2*Zh33*la1s+Zh23*Zh32
     &      **2*Zh33*la2s-2D0*Zh23*Zh31*Zh32*Zh33*la4s-2D0
     &      *Zh23*Zh31*Zh32*Zh33*la5s-2D0*Zh23*Zh31*Zh32
     &      *Zh33*la6s
      AAABR(1887) = 3D0*Zh11*Zh21**2*Zh31*la1-3D0*Zh11*Zh21
     &      **2*Zh32*la6+Zh11*Zh22**2*Zh31*la3+Zh11*Zh22
     &      **2*Zh31*la4+Zh11*Zh22**2*Zh31*la5-3D0*Zh11*Zh22
     &      **2*Zh32*la7+Zh11*Zh23**2*Zh31*la1s-Zh11*Zh23
     &      **2*Zh32*la4s-Zh11*Zh23**2*Zh32*la5s-Zh11*Zh23
     &      **2*Zh32*la6s-3D0*Zh12*Zh21**2*Zh31*la6+Zh12
     &      *Zh21**2*Zh32*la3+Zh12*Zh21**2*Zh32*la4+Zh12
     &      *Zh21**2*Zh32*la5-3D0*Zh12*Zh22**2*Zh31*la7+3D0
     &      *Zh12*Zh22**2*Zh32*la2-Zh12*Zh23**2*Zh31*la4s
     &      -Zh12*Zh23**2*Zh31*la5s-Zh12*Zh23**2*Zh31*la6s
     &      +Zh12*Zh23**2*Zh32*la2s+Zh13*Zh21**2*Zh33*la1s
     &      +Zh13*Zh22**2*Zh33*la2s+6D0*Zh13*Zh23**2*Zh33
     &      *la3s+12D0*Zh13*Zh23**2*Zh33*la7s+12D0*Zh13*Zh23
     &      **2*Zh33*la8s-6D0*Zh11*Zh21*Zh22*Zh31*la6+2D0
     &      *Zh11*Zh21*Zh22*Zh32*la3+2D0*Zh11*Zh21*Zh22*Zh32
     &      *la4+2D0*Zh11*Zh21*Zh22*Zh32*la5+2D0*Zh11*Zh21
     &      *Zh23*Zh33*la1s-2D0*Zh11*Zh22*Zh23*Zh33*la4s
     &      -2D0*Zh11*Zh22*Zh23*Zh33*la5s-2D0*Zh11*Zh22*Zh23
     &      *Zh33*la6s+2D0*Zh12*Zh21*Zh22*Zh31*la3+2D0*Zh12
     &      *Zh21*Zh22*Zh31*la4+2D0*Zh12*Zh21*Zh22*Zh31*la5
     &      -6D0*Zh12*Zh21*Zh22*Zh32*la7-2D0*Zh12*Zh21*Zh23
     &      *Zh33*la4s-2D0*Zh12*Zh21*Zh23*Zh33*la5s-2D0*Zh12
     &      *Zh21*Zh23*Zh33*la6s+2D0*Zh12*Zh22*Zh23*Zh33
     &      *la2s-2D0*Zh13*Zh21*Zh22*Zh33*la4s-2D0*Zh13*Zh21
     &      *Zh22*Zh33*la5s-2D0*Zh13*Zh21*Zh22*Zh33*la6s
     &      +2D0*Zh13*Zh21*Zh23*Zh31*la1s-2D0*Zh13*Zh21*Zh23
     &      *Zh32*la4s-2D0*Zh13*Zh21*Zh23*Zh32*la5s-2D0*Zh13
     &      *Zh21*Zh23*Zh32*la6s-2D0*Zh13*Zh22*Zh23*Zh31
     &      *la4s-2D0*Zh13*Zh22*Zh23*Zh31*la5s-2D0*Zh13*Zh22
     &      *Zh23*Zh31*la6s+2D0*Zh13*Zh22*Zh23*Zh32*la2s
      AAABR(1888) = 3D0*Zh21**2*Zh31**2*la1+Zh21**2*Zh32
     &      **2*la3+Zh21**2*Zh32**2*la4+Zh21**2*Zh32**2*la5
     &      +Zh21**2*Zh33**2*la1s+Zh22**2*Zh31**2*la3+Zh22
     &      **2*Zh31**2*la4+Zh22**2*Zh31**2*la5+3D0*Zh22
     &      **2*Zh32**2*la2+Zh22**2*Zh33**2*la2s+Zh23**2
     &      *Zh31**2*la1s+Zh23**2*Zh32**2*la2s+6D0*Zh23**2
     &      *Zh33**2*la3s+12D0*Zh23**2*Zh33**2*la7s+12D0
     &      *Zh23**2*Zh33**2*la8s-6D0*Zh21*Zh22*Zh31**2*la6
     &      -6D0*Zh21*Zh22*Zh32**2*la7-2D0*Zh21*Zh22*Zh33
     &      **2*la4s-2D0*Zh21*Zh22*Zh33**2*la5s-2D0*Zh21
     &      *Zh22*Zh33**2*la6s-6D0*Zh21**2*Zh31*Zh32*la6
     &      -6D0*Zh22**2*Zh31*Zh32*la7-2D0*Zh23**2*Zh31*Zh32
     &      *la4s-2D0*Zh23**2*Zh31*Zh32*la5s-2D0*Zh23**2
     &      *Zh31*Zh32*la6s+4D0*Zh21*Zh22*Zh31*Zh32*la3+4D0
     &      *Zh21*Zh22*Zh31*Zh32*la4+4D0*Zh21*Zh22*Zh31*Zh32
     &      *la5+4D0*Zh21*Zh23*Zh31*Zh33*la1s-4D0*Zh21*Zh23
     &      *Zh32*Zh33*la4s-4D0*Zh21*Zh23*Zh32*Zh33*la5s
     &      -4D0*Zh21*Zh23*Zh32*Zh33*la6s-4D0*Zh22*Zh23*Zh31
     &      *Zh33*la4s-4D0*Zh22*Zh23*Zh31*Zh33*la5s-4D0*Zh22
     &      *Zh23*Zh31*Zh33*la6s+4D0*Zh22*Zh23*Zh32*Zh33
     &      *la2s
      AAABR(1889) = Zh11*Zh21**3*la1-Zh11*Zh22**3*la7-Zh12
     &      *Zh21**3*la6+Zh12*Zh22**3*la2+2D0*Zh13*Zh23**3
     &      *la3s+4D0*Zh13*Zh23**3*la7s+4D0*Zh13*Zh23**3
     &      *la8s+Zh11*Zh21*Zh22**2*la3+Zh11*Zh21*Zh22**2
     &      *la4+Zh11*Zh21*Zh22**2*la5+Zh11*Zh21*Zh23**2
     &      *la1s-3D0*Zh11*Zh21**2*Zh22*la6-Zh11*Zh22*Zh23
     &      **2*la4s-Zh11*Zh22*Zh23**2*la5s-Zh11*Zh22*Zh23
     &      **2*la6s-3D0*Zh12*Zh21*Zh22**2*la7-Zh12*Zh21
     &      *Zh23**2*la4s-Zh12*Zh21*Zh23**2*la5s-Zh12*Zh21
     &      *Zh23**2*la6s+Zh12*Zh21**2*Zh22*la3+Zh12*Zh21
     &      **2*Zh22*la4+Zh12*Zh21**2*Zh22*la5+Zh12*Zh22
     &      *Zh23**2*la2s+Zh13*Zh21**2*Zh23*la1s+Zh13*Zh22
     &      **2*Zh23*la2s-2D0*Zh13*Zh21*Zh22*Zh23*la4s-2D0
     &      *Zh13*Zh21*Zh22*Zh23*la5s-2D0*Zh13*Zh21*Zh22
     &      *Zh23*la6s
      AAABR(1890) = 3D0*Zh11**2*Zh21**2*la1+Zh11**2*Zh22
     &      **2*la3+Zh11**2*Zh22**2*la4+Zh11**2*Zh22**2*la5
     &      +Zh11**2*Zh23**2*la1s+Zh12**2*Zh21**2*la3+Zh12
     &      **2*Zh21**2*la4+Zh12**2*Zh21**2*la5+3D0*Zh12
     &      **2*Zh22**2*la2+Zh12**2*Zh23**2*la2s+Zh13**2
     &      *Zh21**2*la1s+Zh13**2*Zh22**2*la2s+6D0*Zh13**2
     &      *Zh23**2*la3s+12D0*Zh13**2*Zh23**2*la7s+12D0
     &      *Zh13**2*Zh23**2*la8s-6D0*Zh11*Zh12*Zh21**2*la6
     &      -6D0*Zh11*Zh12*Zh22**2*la7-2D0*Zh11*Zh12*Zh23
     &      **2*la4s-2D0*Zh11*Zh12*Zh23**2*la5s-2D0*Zh11
     &      *Zh12*Zh23**2*la6s-6D0*Zh11**2*Zh21*Zh22*la6
     &      -6D0*Zh12**2*Zh21*Zh22*la7-2D0*Zh13**2*Zh21*Zh22
     &      *la4s-2D0*Zh13**2*Zh21*Zh22*la5s-2D0*Zh13**2
     &      *Zh21*Zh22*la6s+4D0*Zh11*Zh12*Zh21*Zh22*la3+4D0
     &      *Zh11*Zh12*Zh21*Zh22*la4+4D0*Zh11*Zh12*Zh21*Zh22
     &      *la5+4D0*Zh11*Zh13*Zh21*Zh23*la1s-4D0*Zh11*Zh13
     &      *Zh22*Zh23*la4s-4D0*Zh11*Zh13*Zh22*Zh23*la5s
     &      -4D0*Zh11*Zh13*Zh22*Zh23*la6s-4D0*Zh12*Zh13*Zh21
     &      *Zh23*la4s-4D0*Zh12*Zh13*Zh21*Zh23*la5s-4D0*Zh12
     &      *Zh13*Zh21*Zh23*la6s+4D0*Zh12*Zh13*Zh22*Zh23
     &      *la2s
      AAABR(1891) = Zh21**3*Zh31*la1-Zh21**3*Zh32*la6-Zh22
     &      **3*Zh31*la7+Zh22**3*Zh32*la2+2D0*Zh23**3*Zh33
     &      *la3s+4D0*Zh23**3*Zh33*la7s+4D0*Zh23**3*Zh33
     &      *la8s+Zh21*Zh22**2*Zh31*la3+Zh21*Zh22**2*Zh31
     &      *la4+Zh21*Zh22**2*Zh31*la5-3D0*Zh21*Zh22**2*Zh32
     &      *la7+Zh21*Zh23**2*Zh31*la1s-Zh21*Zh23**2*Zh32
     &      *la4s-Zh21*Zh23**2*Zh32*la5s-Zh21*Zh23**2*Zh32
     &      *la6s-3D0*Zh21**2*Zh22*Zh31*la6+Zh21**2*Zh22
     &      *Zh32*la3+Zh21**2*Zh22*Zh32*la4+Zh21**2*Zh22
     &      *Zh32*la5+Zh21**2*Zh23*Zh33*la1s-Zh22*Zh23**2
     &      *Zh31*la4s-Zh22*Zh23**2*Zh31*la5s-Zh22*Zh23**2
     &      *Zh31*la6s+Zh22*Zh23**2*Zh32*la2s+Zh22**2*Zh23
     &      *Zh33*la2s-2D0*Zh21*Zh22*Zh23*Zh33*la4s-2D0*Zh21
     &      *Zh22*Zh23*Zh33*la5s-2D0*Zh21*Zh22*Zh23*Zh33
     &      *la6s
      AAABR(1892) = 3D0*Zh11**2*Zh31**2*la1+Zh11**2*Zh32
     &      **2*la3+Zh11**2*Zh32**2*la4+Zh11**2*Zh32**2*la5
     &      +Zh11**2*Zh33**2*la1s+Zh12**2*Zh31**2*la3+Zh12
     &      **2*Zh31**2*la4+Zh12**2*Zh31**2*la5+3D0*Zh12
     &      **2*Zh32**2*la2+Zh12**2*Zh33**2*la2s+Zh13**2
     &      *Zh31**2*la1s+Zh13**2*Zh32**2*la2s+6D0*Zh13**2
     &      *Zh33**2*la3s+12D0*Zh13**2*Zh33**2*la7s+12D0
     &      *Zh13**2*Zh33**2*la8s-6D0*Zh11*Zh12*Zh31**2*la6
     &      -6D0*Zh11*Zh12*Zh32**2*la7-2D0*Zh11*Zh12*Zh33
     &      **2*la4s-2D0*Zh11*Zh12*Zh33**2*la5s-2D0*Zh11
     &      *Zh12*Zh33**2*la6s-6D0*Zh11**2*Zh31*Zh32*la6
     &      -6D0*Zh12**2*Zh31*Zh32*la7-2D0*Zh13**2*Zh31*Zh32
     &      *la4s-2D0*Zh13**2*Zh31*Zh32*la5s-2D0*Zh13**2
     &      *Zh31*Zh32*la6s+4D0*Zh11*Zh12*Zh31*Zh32*la3+4D0
     &      *Zh11*Zh12*Zh31*Zh32*la4+4D0*Zh11*Zh12*Zh31*Zh32
     &      *la5+4D0*Zh11*Zh13*Zh31*Zh33*la1s-4D0*Zh11*Zh13
     &      *Zh32*Zh33*la4s-4D0*Zh11*Zh13*Zh32*Zh33*la5s
     &      -4D0*Zh11*Zh13*Zh32*Zh33*la6s-4D0*Zh12*Zh13*Zh31
     &      *Zh33*la4s-4D0*Zh12*Zh13*Zh31*Zh33*la5s-4D0*Zh12
     &      *Zh13*Zh31*Zh33*la6s+4D0*Zh12*Zh13*Zh32*Zh33
     &      *la2s
      AAABR(1893) = 3D0*Zh11**2*Zh21*Zh31*la1-3D0*Zh11**2
     &      *Zh21*Zh32*la6-3D0*Zh11**2*Zh22*Zh31*la6+Zh11
     &      **2*Zh22*Zh32*la3+Zh11**2*Zh22*Zh32*la4+Zh11
     &      **2*Zh22*Zh32*la5+Zh11**2*Zh23*Zh33*la1s+Zh12
     &      **2*Zh21*Zh31*la3+Zh12**2*Zh21*Zh31*la4+Zh12
     &      **2*Zh21*Zh31*la5-3D0*Zh12**2*Zh21*Zh32*la7-3D0
     &      *Zh12**2*Zh22*Zh31*la7+3D0*Zh12**2*Zh22*Zh32
     &      *la2+Zh12**2*Zh23*Zh33*la2s+Zh13**2*Zh21*Zh31
     &      *la1s-Zh13**2*Zh21*Zh32*la4s-Zh13**2*Zh21*Zh32
     &      *la5s-Zh13**2*Zh21*Zh32*la6s-Zh13**2*Zh22*Zh31
     &      *la4s-Zh13**2*Zh22*Zh31*la5s-Zh13**2*Zh22*Zh31
     &      *la6s+Zh13**2*Zh22*Zh32*la2s+6D0*Zh13**2*Zh23
     &      *Zh33*la3s+12D0*Zh13**2*Zh23*Zh33*la7s+12D0*Zh13
     &      **2*Zh23*Zh33*la8s-6D0*Zh11*Zh12*Zh21*Zh31*la6
     &      +2D0*Zh11*Zh12*Zh21*Zh32*la3+2D0*Zh11*Zh12*Zh21
     &      *Zh32*la4+2D0*Zh11*Zh12*Zh21*Zh32*la5+2D0*Zh11
     &      *Zh12*Zh22*Zh31*la3+2D0*Zh11*Zh12*Zh22*Zh31*la4
     &      +2D0*Zh11*Zh12*Zh22*Zh31*la5-6D0*Zh11*Zh12*Zh22
     &      *Zh32*la7-2D0*Zh11*Zh12*Zh23*Zh33*la4s-2D0*Zh11
     &      *Zh12*Zh23*Zh33*la5s-2D0*Zh11*Zh12*Zh23*Zh33
     &      *la6s+2D0*Zh11*Zh13*Zh21*Zh33*la1s-2D0*Zh11*Zh13
     &      *Zh22*Zh33*la4s-2D0*Zh11*Zh13*Zh22*Zh33*la5s
     &      -2D0*Zh11*Zh13*Zh22*Zh33*la6s+2D0*Zh11*Zh13*Zh23
     &      *Zh31*la1s-2D0*Zh11*Zh13*Zh23*Zh32*la4s-2D0*Zh11
     &      *Zh13*Zh23*Zh32*la5s-2D0*Zh11*Zh13*Zh23*Zh32
     &      *la6s-2D0*Zh12*Zh13*Zh21*Zh33*la4s-2D0*Zh12*Zh13
     &      *Zh21*Zh33*la5s-2D0*Zh12*Zh13*Zh21*Zh33*la6s
     &      +2D0*Zh12*Zh13*Zh22*Zh33*la2s-2D0*Zh12*Zh13*Zh23
     &      *Zh31*la4s-2D0*Zh12*Zh13*Zh23*Zh31*la5s-2D0*Zh12
     &      *Zh13*Zh23*Zh31*la6s+2D0*Zh12*Zh13*Zh23*Zh32
     &      *la2s
      AAABR(1894) = 3D0*Zh11*Zh21*Zh31**2*la1+Zh11*Zh21*Zh32
     &      **2*la3+Zh11*Zh21*Zh32**2*la4+Zh11*Zh21*Zh32
     &      **2*la5+Zh11*Zh21*Zh33**2*la1s-3D0*Zh11*Zh22
     &      *Zh31**2*la6-3D0*Zh11*Zh22*Zh32**2*la7-Zh11*Zh22
     &      *Zh33**2*la4s-Zh11*Zh22*Zh33**2*la5s-Zh11*Zh22
     &      *Zh33**2*la6s-3D0*Zh12*Zh21*Zh31**2*la6-3D0*Zh12
     &      *Zh21*Zh32**2*la7-Zh12*Zh21*Zh33**2*la4s-Zh12
     &      *Zh21*Zh33**2*la5s-Zh12*Zh21*Zh33**2*la6s+Zh12
     &      *Zh22*Zh31**2*la3+Zh12*Zh22*Zh31**2*la4+Zh12
     &      *Zh22*Zh31**2*la5+3D0*Zh12*Zh22*Zh32**2*la2+Zh12
     &      *Zh22*Zh33**2*la2s+Zh13*Zh23*Zh31**2*la1s+Zh13
     &      *Zh23*Zh32**2*la2s+6D0*Zh13*Zh23*Zh33**2*la3s
     &      +12D0*Zh13*Zh23*Zh33**2*la7s+12D0*Zh13*Zh23*Zh33
     &      **2*la8s-6D0*Zh11*Zh21*Zh31*Zh32*la6+2D0*Zh11
     &      *Zh22*Zh31*Zh32*la3+2D0*Zh11*Zh22*Zh31*Zh32*la4
     &      +2D0*Zh11*Zh22*Zh31*Zh32*la5+2D0*Zh11*Zh23*Zh31
     &      *Zh33*la1s-2D0*Zh11*Zh23*Zh32*Zh33*la4s-2D0*Zh11
     &      *Zh23*Zh32*Zh33*la5s-2D0*Zh11*Zh23*Zh32*Zh33
     &      *la6s+2D0*Zh12*Zh21*Zh31*Zh32*la3+2D0*Zh12*Zh21
     &      *Zh31*Zh32*la4+2D0*Zh12*Zh21*Zh31*Zh32*la5-6D0
     &      *Zh12*Zh22*Zh31*Zh32*la7-2D0*Zh12*Zh23*Zh31*Zh33
     &      *la4s-2D0*Zh12*Zh23*Zh31*Zh33*la5s-2D0*Zh12*Zh23
     &      *Zh31*Zh33*la6s+2D0*Zh12*Zh23*Zh32*Zh33*la2s
     &      +2D0*Zh13*Zh21*Zh31*Zh33*la1s-2D0*Zh13*Zh21*Zh32
     &      *Zh33*la4s-2D0*Zh13*Zh21*Zh32*Zh33*la5s-2D0*Zh13
     &      *Zh21*Zh32*Zh33*la6s-2D0*Zh13*Zh22*Zh31*Zh33
     &      *la4s-2D0*Zh13*Zh22*Zh31*Zh33*la5s-2D0*Zh13*Zh22
     &      *Zh31*Zh33*la6s+2D0*Zh13*Zh22*Zh32*Zh33*la2s
     &      -2D0*Zh13*Zh23*Zh31*Zh32*la4s-2D0*Zh13*Zh23*Zh31
     &      *Zh32*la5s-2D0*Zh13*Zh23*Zh31*Zh32*la6s
      AAABR(1895) = EE**2*(SW/CW-CW/SW)
      AAABR(1896) = EE**2
      AAABR(1897) = EE**2*(2D0-SW**2/CW**2-CW**2/SW**2)
      AAABR(1898) = EE**2/CW*Pa21
      AAABR(1899) = EE**2*Pa21/SW
      AAABR(1900) = EE**2/CW*Pa11
      AAABR(1901) = EE**2*Pa11/SW
      AAABR(1902) = EE**2/CW*(Zh21*sb-Zh22*cb)
      AAABR(1903) = EE**2/SW*(Zh21*sb-Zh22*cb)
      AAABR(1904) = EE**2/CW*(Zh11*sb-Zh12*cb)
      AAABR(1905) = EE**2/SW*(Zh11*sb-Zh12*cb)
      AAABR(1906) = EE**2/CW*(Zh31*sb-Zh32*cb)
      AAABR(1907) = EE**2/SW*(Zh31*sb-Zh32*cb)
      AAABR(1908) = EE**2*(SW/CW*Zl11**2+2D0*SW/CW*Zl12**2
     &      -CW/SW*Zl11**2)
      AAABR(1909) = EE**2*(2D0*Zl11**2-SW**2/CW**2*Zl11**2
     &      -4D0*SW**2/CW**2*Zl12**2-CW**2/SW**2*Zl11**2)
      AAABR(1910) = EE**2/SW**2*Zl11**2
      AAABR(1911) = EE**2*(SW/CW*Zl11*Zl21+2D0*SW/CW*Zl12
     &      *Zl22-CW/SW*Zl11*Zl21)
      AAABR(1912) = EE**2*(2D0*Zl11*Zl21-SW**2/CW**2*Zl11
     &      *Zl21-4D0*SW**2/CW**2*Zl12*Zl22-CW**2/SW**2*Zl11
     &      *Zl21)
      AAABR(1913) = EE**2/SW**2*Zl11*Zl21
      AAABR(1914) = EE**2/CW*Sqrt2*Zl11
      AAABR(1915) = EE**2/SW*Sqrt2*Zl11
      AAABR(1916) = EE**2/CW*Sqrt2
      AAABR(1917) = EE**2/SW*Sqrt2
      AAABR(1918) = EE**2/CW*SW
      AAABR(1919) = EE**2*(SW/CW*Zl21**2+2D0*SW/CW*Zl22**2
     &      -CW/SW*Zl21**2)
      AAABR(1920) = EE**2/CW**2*SW**2
      AAABR(1921) = EE**2*(2D0*Zl21**2-SW**2/CW**2*Zl21**2
     &      -4D0*SW**2/CW**2*Zl22**2-CW**2/SW**2*Zl21**2)
      AAABR(1922) = EE**2/SW**2*Zl21**2
      AAABR(1923) = EE**2/CW*Sqrt2*Zl21
      AAABR(1924) = EE**2/SW*Sqrt2*Zl21
      AAABR(1925) = EE**2*(2D0+SW**2/CW**2+CW**2/SW**2)
      AAABR(1926) = EE**2*(SW/CW*Zss11**2-2D0*SW/CW*Zss12
     &      **2+3D0*CW/SW*Zss11**2)
      AAABR(1927) = EE**2*(SW/CW*Zb11**2-2D0*SW/CW*Zb12**2
     &      +3D0*CW/SW*Zb11**2)
      AAABR(1928) = EE**2*(SW/CW*Zdd11**2-2D0*SW/CW*Zdd12
     &      **2+3D0*CW/SW*Zdd11**2)
      AAABR(1929) = EE*GG*(SW/CW*Zb11**2-2D0*SW/CW*Zb12**2
     &      +3D0*CW/SW*Zb11**2)
      AAABR(1930) = EE*GG
      AAABR(1931) = EE*GG*(SW/CW*Zss11**2-2D0*SW/CW*Zss12
     &      **2+3D0*CW/SW*Zss11**2)
      AAABR(1932) = EE*GG*(SW/CW*Zdd11**2-2D0*SW/CW*Zdd12
     &      **2+3D0*CW/SW*Zdd11**2)
      AAABR(1933) = EE**2*(6D0*Zdd11**2+SW**2/CW**2*Zdd11
     &      **2+4D0*SW**2/CW**2*Zdd12**2+9D0*CW**2/SW**2
     &      *Zdd11**2)
      AAABR(1934) = EE**2*(6D0*Zss11**2+SW**2/CW**2*Zss11
     &      **2+4D0*SW**2/CW**2*Zss12**2+9D0*CW**2/SW**2
     &      *Zss11**2)
      AAABR(1935) = EE**2*(6D0*Zb11**2+SW**2/CW**2*Zb11**2
     &      +4D0*SW**2/CW**2*Zb12**2+9D0*CW**2/SW**2*Zb11
     &      **2)
      AAABR(1936) = EE**2/SW**2*Zdd11**2
      AAABR(1937) = EE**2/SW**2*Zb11**2
      AAABR(1938) = EE**2/SW**2*Zss11**2
      AAABR(1939) = EE**2*(SW/CW*Zss11*Zss21-2D0*SW/CW*Zss12
     &      *Zss22+3D0*CW/SW*Zss11*Zss21)
      AAABR(1940) = EE**2*(SW/CW*Zb11*Zb21-2D0*SW/CW*Zb12
     &      *Zb22+3D0*CW/SW*Zb11*Zb21)
      AAABR(1941) = EE**2*(SW/CW*Zdd11*Zdd21-2D0*SW/CW*Zdd12
     &      *Zdd22+3D0*CW/SW*Zdd11*Zdd21)
      AAABR(1942) = EE*GG*(SW/CW*Zdd11*Zdd21-2D0*SW/CW*Zdd12
     &      *Zdd22+3D0*CW/SW*Zdd11*Zdd21)
      AAABR(1943) = EE*GG*(SW/CW*Zb11*Zb21-2D0*SW/CW*Zb12
     &      *Zb22+3D0*CW/SW*Zb11*Zb21)
      AAABR(1944) = EE*GG*(SW/CW*Zss11*Zss21-2D0*SW/CW*Zss12
     &      *Zss22+3D0*CW/SW*Zss11*Zss21)
      AAABR(1945) = EE**2*(6D0*Zb11*Zb21+SW**2/CW**2*Zb11
     &      *Zb21+4D0*SW**2/CW**2*Zb12*Zb22+9D0*CW**2/SW
     &      **2*Zb11*Zb21)
      AAABR(1946) = EE**2*(6D0*Zdd11*Zdd21+SW**2/CW**2*Zdd11
     &      *Zdd21+4D0*SW**2/CW**2*Zdd12*Zdd22+9D0*CW**2
     &      /SW**2*Zdd11*Zdd21)
      AAABR(1947) = EE**2*(6D0*Zss11*Zss21+SW**2/CW**2*Zss11
     &      *Zss21+4D0*SW**2/CW**2*Zss12*Zss22+9D0*CW**2
     &      /SW**2*Zss11*Zss21)
      AAABR(1948) = EE**2/SW**2*Zb11*Zb21
      AAABR(1949) = EE**2/SW**2*Zdd11*Zdd21
      AAABR(1950) = EE**2/SW**2*Zss11*Zss21
      AAABR(1951) = EE**2/CW*Sqrt2*Zcc11*Zss11
      AAABR(1952) = EE**2/SW*Sqrt2*Zcc11*Zss11
      AAABR(1953) = EE**2/CW*Sqrt2*Zdd11*Zuu11
      AAABR(1954) = EE**2/SW*Sqrt2*Zdd11*Zuu11
      AAABR(1955) = EE**2/CW*Sqrt2*Zb11*Zt11
      AAABR(1956) = EE**2/SW*Sqrt2*Zb11*Zt11
      AAABR(1957) = EE*GG/SW*Sqrt2*Zcc11*Zss11
      AAABR(1958) = EE*GG/SW*Sqrt2*Zdd11*Zuu11
      AAABR(1959) = EE*GG/SW*Sqrt2*Zb11*Zt11
      AAABR(1960) = EE**2/CW*Sqrt2*Zdd11*Zuu21
      AAABR(1961) = EE**2/SW*Sqrt2*Zdd11*Zuu21
      AAABR(1962) = EE**2/CW*Sqrt2*Zcc21*Zss11
      AAABR(1963) = EE**2/SW*Sqrt2*Zcc21*Zss11
      AAABR(1964) = EE**2/CW*Sqrt2*Zb11*Zt21
      AAABR(1965) = EE**2/SW*Sqrt2*Zb11*Zt21
      AAABR(1966) = EE*GG/SW*Sqrt2*Zdd11*Zuu21
      AAABR(1967) = EE*GG/SW*Sqrt2*Zb11*Zt21
      AAABR(1968) = EE*GG/SW*Sqrt2*Zcc21*Zss11
      AAABR(1969) = EE**2*(SW/CW*Zdd21**2-2D0*SW/CW*Zdd22
     &      **2+3D0*CW/SW*Zdd21**2)
      AAABR(1970) = EE**2*(SW/CW*Zb21**2-2D0*SW/CW*Zb22**2
     &      +3D0*CW/SW*Zb21**2)
      AAABR(1971) = EE**2*(SW/CW*Zss21**2-2D0*SW/CW*Zss22
     &      **2+3D0*CW/SW*Zss21**2)
      AAABR(1972) = EE*GG*(SW/CW*Zss21**2-2D0*SW/CW*Zss22
     &      **2+3D0*CW/SW*Zss21**2)
      AAABR(1973) = EE*GG*(SW/CW*Zdd21**2-2D0*SW/CW*Zdd22
     &      **2+3D0*CW/SW*Zdd21**2)
      AAABR(1974) = EE*GG*(SW/CW*Zb21**2-2D0*SW/CW*Zb22**2
     &      +3D0*CW/SW*Zb21**2)
      AAABR(1975) = EE**2*(6D0*Zdd21**2+SW**2/CW**2*Zdd21
     &      **2+4D0*SW**2/CW**2*Zdd22**2+9D0*CW**2/SW**2
     &      *Zdd21**2)
      AAABR(1976) = EE**2*(6D0*Zb21**2+SW**2/CW**2*Zb21**2
     &      +4D0*SW**2/CW**2*Zb22**2+9D0*CW**2/SW**2*Zb21
     &      **2)
      AAABR(1977) = EE**2*(6D0*Zss21**2+SW**2/CW**2*Zss21
     &      **2+4D0*SW**2/CW**2*Zss22**2+9D0*CW**2/SW**2
     &      *Zss21**2)
      AAABR(1978) = EE**2/SW**2*Zdd21**2
      AAABR(1979) = EE**2/SW**2*Zb21**2
      AAABR(1980) = EE**2/SW**2*Zss21**2
      AAABR(1981) = EE**2/CW*Sqrt2*Zcc11*Zss21
      AAABR(1982) = EE**2/SW*Sqrt2*Zcc11*Zss21
      AAABR(1983) = EE**2/CW*Sqrt2*Zb21*Zt11
      AAABR(1984) = EE**2/SW*Sqrt2*Zb21*Zt11
      AAABR(1985) = EE**2/CW*Sqrt2*Zdd21*Zuu11
      AAABR(1986) = EE**2/SW*Sqrt2*Zdd21*Zuu11
      AAABR(1987) = EE*GG/SW*Sqrt2*Zdd21*Zuu11
      AAABR(1988) = EE*GG/SW*Sqrt2*Zb21*Zt11
      AAABR(1989) = EE*GG/SW*Sqrt2*Zcc11*Zss21
      AAABR(1990) = EE**2/CW*Sqrt2*Zdd21*Zuu21
      AAABR(1991) = EE**2/SW*Sqrt2*Zdd21*Zuu21
      AAABR(1992) = EE**2/CW*Sqrt2*Zcc21*Zss21
      AAABR(1993) = EE**2/SW*Sqrt2*Zcc21*Zss21
      AAABR(1994) = EE**2/CW*Sqrt2*Zb21*Zt21
      AAABR(1995) = EE**2/SW*Sqrt2*Zb21*Zt21
      AAABR(1996) = EE*GG/SW*Sqrt2*Zcc21*Zss21
      AAABR(1997) = EE*GG/SW*Sqrt2*Zb21*Zt21
      AAABR(1998) = EE*GG/SW*Sqrt2*Zdd21*Zuu21
      AAABR(1999) = EE**2*(SW/CW*Zt11**2+4D0*SW/CW*Zt12**2
     &      -3D0*CW/SW*Zt11**2)
      AAABR(2000) = EE**2*(SW/CW*Zcc11**2+4D0*SW/CW*Zcc12
     &      **2-3D0*CW/SW*Zcc11**2)
      end

      subroutine aaini03
      implicit none
#include "model.h"

      AAABR(2001) = EE**2*(SW/CW*Zuu11**2+4D0*SW/CW*Zuu12
     &      **2-3D0*CW/SW*Zuu11**2)
      AAABR(2002) = EE*GG*(SW/CW*Zuu11**2+4D0*SW/CW*Zuu12
     &      **2-3D0*CW/SW*Zuu11**2)
      AAABR(2003) = EE*GG*(SW/CW*Zt11**2+4D0*SW/CW*Zt12**2
     &      -3D0*CW/SW*Zt11**2)
      AAABR(2004) = EE*GG*(SW/CW*Zcc11**2+4D0*SW/CW*Zcc12
     &      **2-3D0*CW/SW*Zcc11**2)
      AAABR(2005) = EE**2*(6D0*Zcc11**2-SW**2/CW**2*Zcc11
     &      **2-16D0*SW**2/CW**2*Zcc12**2-9D0*CW**2/SW**2
     &      *Zcc11**2)
      AAABR(2006) = EE**2*(6D0*Zt11**2-SW**2/CW**2*Zt11**2
     &      -16D0*SW**2/CW**2*Zt12**2-9D0*CW**2/SW**2*Zt11
     &      **2)
      AAABR(2007) = EE**2*(6D0*Zuu11**2-SW**2/CW**2*Zuu11
     &      **2-16D0*SW**2/CW**2*Zuu12**2-9D0*CW**2/SW**2
     &      *Zuu11**2)
      AAABR(2008) = EE**2/SW**2*Zt11**2
      AAABR(2009) = EE**2/SW**2*Zuu11**2
      AAABR(2010) = EE**2/SW**2*Zcc11**2
      AAABR(2011) = EE**2*(SW/CW*Zt11*Zt21+4D0*SW/CW*Zt12
     &      *Zt22-3D0*CW/SW*Zt11*Zt21)
      AAABR(2012) = EE**2*(SW/CW*Zuu11*Zuu21+4D0*SW/CW*Zuu12
     &      *Zuu22-3D0*CW/SW*Zuu11*Zuu21)
      AAABR(2013) = EE**2*(SW/CW*Zcc11*Zcc21+4D0*SW/CW*Zcc12
     &      *Zcc22-3D0*CW/SW*Zcc11*Zcc21)
      AAABR(2014) = EE*GG*(SW/CW*Zcc11*Zcc21+4D0*SW/CW*Zcc12
     &      *Zcc22-3D0*CW/SW*Zcc11*Zcc21)
      AAABR(2015) = EE*GG*(SW/CW*Zt11*Zt21+4D0*SW/CW*Zt12
     &      *Zt22-3D0*CW/SW*Zt11*Zt21)
      AAABR(2016) = EE*GG*(SW/CW*Zuu11*Zuu21+4D0*SW/CW*Zuu12
     &      *Zuu22-3D0*CW/SW*Zuu11*Zuu21)
      AAABR(2017) = EE**2*(6D0*Zcc11*Zcc21-SW**2/CW**2*Zcc11
     &      *Zcc21-16D0*SW**2/CW**2*Zcc12*Zcc22-9D0*CW**2
     &      /SW**2*Zcc11*Zcc21)
      AAABR(2018) = EE**2*(6D0*Zuu11*Zuu21-SW**2/CW**2*Zuu11
     &      *Zuu21-16D0*SW**2/CW**2*Zuu12*Zuu22-9D0*CW**2
     &      /SW**2*Zuu11*Zuu21)
      AAABR(2019) = EE**2*(6D0*Zt11*Zt21-SW**2/CW**2*Zt11
     &      *Zt21-16D0*SW**2/CW**2*Zt12*Zt22-9D0*CW**2/SW
     &      **2*Zt11*Zt21)
      AAABR(2020) = EE**2/SW**2*Zcc11*Zcc21
      AAABR(2021) = EE**2/SW**2*Zt11*Zt21
      AAABR(2022) = EE**2/SW**2*Zuu11*Zuu21
      AAABR(2023) = EE**2*(SW/CW*Zcc21**2+4D0*SW/CW*Zcc22
     &      **2-3D0*CW/SW*Zcc21**2)
      AAABR(2024) = EE**2*(SW/CW*Zt21**2+4D0*SW/CW*Zt22**2
     &      -3D0*CW/SW*Zt21**2)
      AAABR(2025) = EE**2*(SW/CW*Zuu21**2+4D0*SW/CW*Zuu22
     &      **2-3D0*CW/SW*Zuu21**2)
      AAABR(2026) = EE*GG*(SW/CW*Zcc21**2+4D0*SW/CW*Zcc22
     &      **2-3D0*CW/SW*Zcc21**2)
      AAABR(2027) = EE*GG*(SW/CW*Zt21**2+4D0*SW/CW*Zt22**2
     &      -3D0*CW/SW*Zt21**2)
      AAABR(2028) = EE*GG*(SW/CW*Zuu21**2+4D0*SW/CW*Zuu22
     &      **2-3D0*CW/SW*Zuu21**2)
      AAABR(2029) = EE**2*(6D0*Zt21**2-SW**2/CW**2*Zt21**2
     &      -16D0*SW**2/CW**2*Zt22**2-9D0*CW**2/SW**2*Zt21
     &      **2)
      AAABR(2030) = EE**2*(6D0*Zuu21**2-SW**2/CW**2*Zuu21
     &      **2-16D0*SW**2/CW**2*Zuu22**2-9D0*CW**2/SW**2
     &      *Zuu21**2)
      AAABR(2031) = EE**2*(6D0*Zcc21**2-SW**2/CW**2*Zcc21
     &      **2-16D0*SW**2/CW**2*Zcc22**2-9D0*CW**2/SW**2
     &      *Zcc21**2)
      AAABR(2032) = EE**2/SW**2*Zt21**2
      AAABR(2033) = EE**2/SW**2*Zuu21**2
      AAABR(2034) = EE**2/SW**2*Zcc21**2
      AAABR(2035) = EE**2/CW
      AAABR(2036) = EE**2/SW
      AAABR(2037) = EE**2/CW*(Zh31*cb+Zh32*sb)
      AAABR(2038) = EE**2/SW*(Zh31*cb+Zh32*sb)
      AAABR(2039) = EE**2/CW*(Zh11*cb+Zh12*sb)
      AAABR(2040) = EE**2/SW*(Zh11*cb+Zh12*sb)
      AAABR(2041) = EE**2/CW*(Zh21*cb+Zh22*sb)
      AAABR(2042) = EE**2/SW*(Zh21*cb+Zh22*sb)
      AAABR(2043) = EE**2*Pa11**2*(2D0+SW**2/CW**2+CW**2
     &      /SW**2)
      AAABR(2044) = EE**2*Pa11*Pa21*(2D0+SW**2/CW**2+CW**2
     &      /SW**2)
      AAABR(2045) = EE**2*Pa21**2*(2D0+SW**2/CW**2+CW**2
     &      /SW**2)
      AAABR(2046) = EE**2*Pa11**2/SW**2
      AAABR(2047) = EE**2*Pa11*Pa21/SW**2
      AAABR(2048) = EE**2*Pa21**2/SW**2
      AAABR(2049) = EE**2*(2D0*Zh21*Zh31+2D0*Zh22*Zh32+SW
     &      **2/CW**2*Zh21*Zh31+SW**2/CW**2*Zh22*Zh32+CW
     &      **2/SW**2*Zh21*Zh31+CW**2/SW**2*Zh22*Zh32)
      AAABR(2050) = EE**2*(2D0*Zh11*Zh21+2D0*Zh12*Zh22+SW
     &      **2/CW**2*Zh11*Zh21+SW**2/CW**2*Zh12*Zh22+CW
     &      **2/SW**2*Zh11*Zh21+CW**2/SW**2*Zh12*Zh22)
      AAABR(2051) = EE**2*(2D0*Zh31**2+2D0*Zh32**2+SW**2
     &      /CW**2*Zh31**2+SW**2/CW**2*Zh32**2+CW**2/SW**2
     &      *Zh31**2+CW**2/SW**2*Zh32**2)
      AAABR(2052) = EE**2*(2D0*Zh11*Zh31+2D0*Zh12*Zh32+SW
     &      **2/CW**2*Zh11*Zh31+SW**2/CW**2*Zh12*Zh32+CW
     &      **2/SW**2*Zh11*Zh31+CW**2/SW**2*Zh12*Zh32)
      AAABR(2053) = EE**2*(2D0*Zh21**2+2D0*Zh22**2+SW**2
     &      /CW**2*Zh21**2+SW**2/CW**2*Zh22**2+CW**2/SW**2
     &      *Zh21**2+CW**2/SW**2*Zh22**2)
      AAABR(2054) = EE**2*(2D0*Zh11**2+2D0*Zh12**2+SW**2
     &      /CW**2*Zh11**2+SW**2/CW**2*Zh12**2+CW**2/SW**2
     &      *Zh11**2+CW**2/SW**2*Zh12**2)
      AAABR(2055) = EE**2/SW**2*(Zh21**2+Zh22**2)
      AAABR(2056) = EE**2/SW**2*(Zh11*Zh31+Zh12*Zh32)
      AAABR(2057) = EE**2/SW**2*(Zh11*Zh21+Zh12*Zh22)
      AAABR(2058) = EE**2/SW**2*(Zh11**2+Zh12**2)
      AAABR(2059) = EE**2/SW**2*(Zh31**2+Zh32**2)
      AAABR(2060) = EE**2/SW**2*(Zh21*Zh31+Zh22*Zh32)
      AAABR(2061) = CW*EE**2/SW
      AAABR(2062) = CW**2*EE**2/SW**2
      end

      subroutine mtrini
      implicit none
#include "model.h"

      integer m1,m2,m3,m4

      do m1 = 1,3
      MTR005(m1)=0D0
      MTR006(m1)=0D0
      MTR007(m1)=0D0
      MTR008(m1)=0D0
      MTR009(m1)=0D0
      MTR011(m1)=0D0
      MTR012(m1)=0D0
      MTR018(m1)=0D0
      MTR019(m1)=0D0
      MTR020(m1)=0D0
      MTR024(m1)=0D0
      MTR025(m1)=0D0
      MTR027(m1)=0D0
      MTR032(m1)=0D0
      MTR037(m1)=0D0
      MTR038(m1)=0D0
      MTR039(m1)=0D0
      MTR040(m1)=0D0
      MTR041(m1)=0D0
      MTR042(m1)=0D0
      MTR043(m1)=0D0
      MTR044(m1)=0D0
      MTR045(m1)=0D0
      MTR046(m1)=0D0
      MTR047(m1)=0D0
      MTR048(m1)=0D0
      MTR049(m1)=0D0
      MTR050(m1)=0D0
      MTR051(m1)=0D0
      MTR052(m1)=0D0
      MTR054(m1)=0D0
      MTR084(m1)=0D0
      MTR087(m1)=0D0
      MTR088(m1)=0D0
      MTR089(m1)=0D0
      MTR090(m1)=0D0
      MTR095(m1)=0D0
      MTR096(m1)=0D0
      MTR097(m1)=0D0
      MTR098(m1)=0D0
      MTR103(m1)=0D0
      MTR106(m1)=0D0
      MTR107(m1)=0D0
      MTR108(m1)=0D0
      MTR109(m1)=0D0
      MTR118(m1)=0D0
      MTR119(m1)=0D0
      MTR120(m1)=0D0
      MTR121(m1)=0D0
      MTR122(m1)=0D0
      MTR123(m1)=0D0
      MTR124(m1)=0D0
      MTR125(m1)=0D0
      MTR126(m1)=0D0
      MTR130(m1)=0D0
      MTR131(m1)=0D0
      MTR132(m1)=0D0
      MTR133(m1)=0D0
      MTR138(m1)=0D0
      MTR139(m1)=0D0
      MTR140(m1)=0D0
      MTR141(m1)=0D0
      MTR142(m1)=0D0
      MTR143(m1)=0D0
      MTR144(m1)=0D0
      MTR149(m1)=0D0
      MTR150(m1)=0D0
      MTR151(m1)=0D0
      MTR157(m1)=0D0
      MTR190(m1)=0D0
      MTR191(m1)=0D0
      MTR194(m1)=0D0
      MTR218(m1)=0D0
      MTR221(m1)=0D0
      MTR233(m1)=0D0
      MTR265(m1)=0D0
      MTR266(m1)=0D0
      MTR267(m1)=0D0
      MTR268(m1)=0D0
      MTR273(m1)=0D0
      MTR274(m1)=0D0
      MTR296(m1)=0D0
      MTR297(m1)=0D0
      MTR298(m1)=0D0
      MTR303(m1)=0D0
      MTR320(m1)=0D0
      MTR321(m1)=0D0
      MTR322(m1)=0D0
      MTR323(m1)=0D0
      MTR334(m1)=0D0
      MTR335(m1)=0D0
      MTR351(m1)=0D0
      MTR352(m1)=0D0
      MTR353(m1)=0D0
      MTR354(m1)=0D0
      MTR355(m1)=0D0
      MTR356(m1)=0D0
      MTR357(m1)=0D0
      MTR358(m1)=0D0
      MTR359(m1)=0D0
      MTR360(m1)=0D0
      MTR361(m1)=0D0
      MTR362(m1)=0D0
      MTR363(m1)=0D0
      MTR364(m1)=0D0
      MTR365(m1)=0D0
      MTR366(m1)=0D0
      MTR367(m1)=0D0
      MTR368(m1)=0D0
      MTR369(m1)=0D0
      MTR370(m1)=0D0
      MTR371(m1)=0D0
      MTR372(m1)=0D0
      MTR373(m1)=0D0
      MTR374(m1)=0D0
      MTR375(m1)=0D0
      MTR376(m1)=0D0
      MTR377(m1)=0D0
      MTR378(m1)=0D0
      MTR379(m1)=0D0
      MTR380(m1)=0D0
      MTR381(m1)=0D0
      MTR382(m1)=0D0
      MTR383(m1)=0D0
      MTR384(m1)=0D0
      MTR385(m1)=0D0
      MTR386(m1)=0D0
      MTR387(m1)=0D0
      MTR388(m1)=0D0
      MTR389(m1)=0D0
      MTR390(m1)=0D0
      MTR391(m1)=0D0
      MTR392(m1)=0D0
      MTR393(m1)=0D0
      MTR394(m1)=0D0
      MTR395(m1)=0D0
      MTR396(m1)=0D0
      MTR397(m1)=0D0
      enddo
      do m1 = 1,2
      do m2 = 1,3
      MTR013(m1,m2)=0D0
      MTR021(m1,m2)=0D0
      MTR028(m1,m2)=0D0
      MTR033(m1,m2)=0D0
      MTR053(m1,m2)=0D0
      MTR061(m1,m2)=0D0
      MTR062(m1,m2)=0D0
      MTR067(m1,m2)=0D0
      MTR068(m1,m2)=0D0
      MTR069(m1,m2)=0D0
      MTR070(m1,m2)=0D0
      MTR077(m1,m2)=0D0
      MTR078(m1,m2)=0D0
      MTR091(m1,m2)=0D0
      MTR092(m1,m2)=0D0
      MTR093(m1,m2)=0D0
      MTR094(m1,m2)=0D0
      MTR134(m1,m2)=0D0
      MTR135(m1,m2)=0D0
      MTR145(m1,m2)=0D0
      MTR146(m1,m2)=0D0
      MTR152(m1,m2)=0D0
      MTR153(m1,m2)=0D0
      MTR159(m1,m2)=0D0
      MTR192(m1,m2)=0D0
      MTR195(m1,m2)=0D0
      MTR196(m1,m2)=0D0
      MTR219(m1,m2)=0D0
      MTR222(m1,m2)=0D0
      MTR234(m1,m2)=0D0
      MTR269(m1,m2)=0D0
      MTR270(m1,m2)=0D0
      MTR275(m1,m2)=0D0
      MTR276(m1,m2)=0D0
      MTR299(m1,m2)=0D0
      MTR300(m1,m2)=0D0
      MTR304(m1,m2)=0D0
      MTR324(m1,m2)=0D0
      MTR325(m1,m2)=0D0
      MTR336(m1,m2)=0D0
      enddo
      enddo
      do m1 = 1,3
      do m2 = 1,3
      MTR015(m1,m2)=0D0
      MTR016(m1,m2)=0D0
      MTR017(m1,m2)=0D0
      MTR022(m1,m2)=0D0
      MTR023(m1,m2)=0D0
      MTR026(m1,m2)=0D0
      MTR029(m1,m2)=0D0
      MTR030(m1,m2)=0D0
      MTR031(m1,m2)=0D0
      MTR072(m1,m2)=0D0
      MTR086(m1,m2)=0D0
      MTR105(m1,m2)=0D0
      MTR129(m1,m2)=0D0
      MTR136(m1,m2)=0D0
      MTR137(m1,m2)=0D0
      MTR147(m1,m2)=0D0
      MTR148(m1,m2)=0D0
      MTR154(m1,m2)=0D0
      MTR155(m1,m2)=0D0
      MTR160(m1,m2)=0D0
      MTR161(m1,m2)=0D0
      MTR162(m1,m2)=0D0
      MTR163(m1,m2)=0D0
      MTR164(m1,m2)=0D0
      MTR165(m1,m2)=0D0
      MTR166(m1,m2)=0D0
      MTR167(m1,m2)=0D0
      MTR168(m1,m2)=0D0
      MTR169(m1,m2)=0D0
      MTR170(m1,m2)=0D0
      MTR171(m1,m2)=0D0
      MTR172(m1,m2)=0D0
      MTR173(m1,m2)=0D0
      MTR174(m1,m2)=0D0
      MTR175(m1,m2)=0D0
      MTR176(m1,m2)=0D0
      MTR177(m1,m2)=0D0
      MTR178(m1,m2)=0D0
      MTR179(m1,m2)=0D0
      MTR180(m1,m2)=0D0
      MTR181(m1,m2)=0D0
      MTR182(m1,m2)=0D0
      MTR183(m1,m2)=0D0
      MTR184(m1,m2)=0D0
      MTR185(m1,m2)=0D0
      MTR186(m1,m2)=0D0
      MTR187(m1,m2)=0D0
      MTR188(m1,m2)=0D0
      MTR189(m1,m2)=0D0
      MTR193(m1,m2)=0D0
      MTR197(m1,m2)=0D0
      MTR203(m1,m2)=0D0
      MTR204(m1,m2)=0D0
      MTR205(m1,m2)=0D0
      MTR206(m1,m2)=0D0
      MTR207(m1,m2)=0D0
      MTR208(m1,m2)=0D0
      MTR209(m1,m2)=0D0
      MTR210(m1,m2)=0D0
      MTR211(m1,m2)=0D0
      MTR212(m1,m2)=0D0
      MTR213(m1,m2)=0D0
      MTR214(m1,m2)=0D0
      MTR215(m1,m2)=0D0
      MTR216(m1,m2)=0D0
      MTR217(m1,m2)=0D0
      MTR220(m1,m2)=0D0
      MTR225(m1,m2)=0D0
      MTR226(m1,m2)=0D0
      MTR227(m1,m2)=0D0
      MTR228(m1,m2)=0D0
      MTR229(m1,m2)=0D0
      MTR230(m1,m2)=0D0
      MTR231(m1,m2)=0D0
      MTR232(m1,m2)=0D0
      MTR237(m1,m2)=0D0
      MTR238(m1,m2)=0D0
      MTR239(m1,m2)=0D0
      MTR240(m1,m2)=0D0
      MTR241(m1,m2)=0D0
      MTR242(m1,m2)=0D0
      MTR243(m1,m2)=0D0
      MTR244(m1,m2)=0D0
      MTR245(m1,m2)=0D0
      MTR246(m1,m2)=0D0
      MTR247(m1,m2)=0D0
      MTR248(m1,m2)=0D0
      MTR249(m1,m2)=0D0
      MTR250(m1,m2)=0D0
      MTR251(m1,m2)=0D0
      MTR252(m1,m2)=0D0
      MTR253(m1,m2)=0D0
      MTR254(m1,m2)=0D0
      MTR255(m1,m2)=0D0
      MTR256(m1,m2)=0D0
      MTR257(m1,m2)=0D0
      MTR258(m1,m2)=0D0
      MTR259(m1,m2)=0D0
      MTR260(m1,m2)=0D0
      MTR261(m1,m2)=0D0
      MTR262(m1,m2)=0D0
      MTR263(m1,m2)=0D0
      MTR264(m1,m2)=0D0
      MTR271(m1,m2)=0D0
      MTR272(m1,m2)=0D0
      MTR277(m1,m2)=0D0
      MTR283(m1,m2)=0D0
      MTR284(m1,m2)=0D0
      MTR285(m1,m2)=0D0
      MTR286(m1,m2)=0D0
      MTR287(m1,m2)=0D0
      MTR288(m1,m2)=0D0
      MTR289(m1,m2)=0D0
      MTR290(m1,m2)=0D0
      MTR291(m1,m2)=0D0
      MTR292(m1,m2)=0D0
      MTR293(m1,m2)=0D0
      MTR294(m1,m2)=0D0
      MTR295(m1,m2)=0D0
      MTR301(m1,m2)=0D0
      MTR302(m1,m2)=0D0
      MTR307(m1,m2)=0D0
      MTR308(m1,m2)=0D0
      MTR309(m1,m2)=0D0
      MTR310(m1,m2)=0D0
      MTR311(m1,m2)=0D0
      MTR312(m1,m2)=0D0
      MTR313(m1,m2)=0D0
      MTR314(m1,m2)=0D0
      MTR315(m1,m2)=0D0
      MTR316(m1,m2)=0D0
      MTR317(m1,m2)=0D0
      MTR318(m1,m2)=0D0
      MTR319(m1,m2)=0D0
      MTR326(m1,m2)=0D0
      MTR332(m1,m2)=0D0
      MTR333(m1,m2)=0D0
      MTR341(m1,m2)=0D0
      MTR343(m1,m2)=0D0
      MTR400(m1,m2)=0D0
      MTR401(m1,m2)=0D0
      enddo
      enddo
      do m1 = 1,3
      do m2 = 1,2
      MTR071(m1,m2)=0D0
      MTR085(m1,m2)=0D0
      MTR104(m1,m2)=0D0
      enddo
      enddo
      do m1 = 1,5
      do m2 = 1,3
      MTR074(m1,m2)=0D0
      MTR075(m1,m2)=0D0
      MTR076(m1,m2)=0D0
      MTR079(m1,m2)=0D0
      MTR080(m1,m2)=0D0
      MTR081(m1,m2)=0D0
      MTR082(m1,m2)=0D0
      MTR083(m1,m2)=0D0
      MTR099(m1,m2)=0D0
      MTR100(m1,m2)=0D0
      MTR101(m1,m2)=0D0
      MTR102(m1,m2)=0D0
      enddo
      enddo
      do m1 = 1,2
      do m2 = 1,2
      do m3 = 1,3
      MTR199(m1,m2,m3)=0D0
      MTR223(m1,m2,m3)=0D0
      MTR235(m1,m2,m3)=0D0
      MTR278(m1,m2,m3)=0D0
      MTR279(m1,m2,m3)=0D0
      MTR305(m1,m2,m3)=0D0
      MTR327(m1,m2,m3)=0D0
      MTR328(m1,m2,m3)=0D0
      MTR337(m1,m2,m3)=0D0
      enddo
      enddo
      enddo
      do m1 = 1,2
      do m2 = 1,3
      do m3 = 1,3
      MTR200(m1,m2,m3)=0D0
      MTR280(m1,m2,m3)=0D0
      MTR329(m1,m2,m3)=0D0
      MTR345(m1,m2,m3)=0D0
      enddo
      enddo
      enddo
      do m1 = 1,3
      do m2 = 1,3
      do m3 = 1,3
      MTR202(m1,m2,m3)=0D0
      MTR224(m1,m2,m3)=0D0
      MTR236(m1,m2,m3)=0D0
      MTR281(m1,m2,m3)=0D0
      MTR282(m1,m2,m3)=0D0
      MTR306(m1,m2,m3)=0D0
      MTR330(m1,m2,m3)=0D0
      MTR331(m1,m2,m3)=0D0
      MTR338(m1,m2,m3)=0D0
      enddo
      enddo
      enddo
      quuMass(1)=Mq
      quuMass(2)=Mc
      quuMass(3)=Mt
      qudMass(1)=Mq
      qudMass(2)=Mq
      qudMass(3)=Mb
      lpdMass(1)=Mele
      lpdMass(2)=Mmu
      lpdMass(3)=Ml
      neuMass(1)=MNE1
      neuMass(2)=MNE2
      neuMass(3)=MNE3
      neuMass(4)=MNE4
      neuMass(5)=MNE5
      chaMass(1)=MC1
      chaMass(2)=MC2
      sluMass(1)=MSne
      sluMass(2)=MSnm
      sluMass(3)=MSnl
      sldMass(1)=MSeL
      sldMass(2)=MSmL
      sldMass(3)=MSl1
      sleMass(1)=MSeR
      sleMass(2)=MSmR
      sleMass(3)=MSl2
      squMass(1)=MSuL
      squMass(2)=MScL
      squMass(3)=MSt1
      sqvMass(1)=MSuR
      sqvMass(2)=MScR
      sqvMass(3)=MSt2
      sqdMass(1)=MSdL
      sqdMass(2)=MSsL
      sqdMass(3)=MSb1
      sqeMass(1)=MSdR
      sqeMass(2)=MSsR
      sqeMass(3)=MSb2
      hisMass(1)=Mh1
      hisMass(2)=Mh2
      hisMass(3)=Mh3
      hiaMass(1)=Mha
      hiaMass(2)=Mhb
      MTR001(1)=AAABR(9)
      MTR001(2)=AAABR(7)
      MTR001(3)=AAABR(8)
      MTR002(1)=AAABR(12)
      MTR002(2)=AAABR(14)
      MTR002(3)=AAABR(13)
      MTR003(1)=AAABR(18)
      MTR003(2)=AAABR(16)
      MTR003(3)=AAABR(17)
      MTR004(1)=2D0*AAABR(19)
      MTR004(2)=2D0*AAABR(19)
      MTR004(3)=AAABR(20)
      MTR005(3)=1D0
      MTR006(1)=AAABR(22)
      MTR006(2)=AAABR(24)
      MTR006(3)=AAABR(23)
      MTR007(1)=AAABR(27)
      MTR007(2)=AAABR(25)
      MTR007(3)=AAABR(26)
      MTR008(1)=AAABR(30)
      MTR008(2)=AAABR(29)
      MTR008(3)=AAABR(28)
      MTR009(1)=AAABR(33)
      MTR009(2)=AAABR(32)
      MTR009(3)=AAABR(31)
      MTR010(1)=AAABR(34)
      MTR010(2)=AAABR(35)
      MTR011(1)=AAABR(36)
      MTR011(2)=AAABR(37)
      MTR011(3)=AAABR(38)
      MTR012(1)=AAABR(39)
      MTR012(2)=AAABR(39)
      MTR012(3)=AAABR(40)
      MTR013(1,3)=AAABR(42)
      MTR013(2,3)=AAABR(43)
      MTR014(1,1)=AAABR(46)
      MTR014(1,2)=AAABR(46)
      MTR014(1,3)=-AAABR(45)
      MTR014(2,1)=AAABR(44)
      MTR014(2,2)=AAABR(44)
      MTR014(2,3)=-AAABR(47)
      MTR014(3,1)=AAABR(48)
      MTR014(3,2)=AAABR(48)
      MTR014(3,3)=-AAABR(49)
      MTR015(1,3)=AAABR(50)
      MTR015(2,3)=AAABR(52)
      MTR015(3,3)=AAABR(51)
      MTR016(1,1)=2D0*AAABR(56)
      MTR016(1,2)=2D0*AAABR(56)
      MTR016(1,3)=-AAABR(59)
      MTR016(2,1)=2D0*AAABR(58)
      MTR016(2,2)=2D0*AAABR(58)
      MTR016(2,3)=-AAABR(55)
      MTR016(3,1)=2D0*AAABR(54)
      MTR016(3,2)=2D0*AAABR(54)
      MTR016(3,3)=-AAABR(57)
      MTR017(1,1)=AAABR(61)
      MTR017(1,2)=AAABR(61)
      MTR017(1,3)=AAABR(61)
      MTR017(2,1)=AAABR(62)
      MTR017(2,2)=AAABR(62)
      MTR017(2,3)=AAABR(62)
      MTR017(3,1)=AAABR(60)
      MTR017(3,2)=AAABR(60)
      MTR017(3,3)=AAABR(60)
      MTR018(1)=AAABR(64)
      MTR018(2)=AAABR(63)
      MTR018(3)=AAABR(65)
      MTR019(1)=AAABR(68)
      MTR019(2)=AAABR(67)
      MTR019(3)=AAABR(66)
      MTR020(1)=AAABR(71)
      MTR020(2)=AAABR(69)
      MTR020(3)=AAABR(70)
      MTR021(1,1)=AAABR(73)
      MTR021(1,2)=AAABR(75)
      MTR021(1,3)=AAABR(76)
      MTR021(2,1)=AAABR(77)
      MTR021(2,2)=AAABR(72)
      MTR021(2,3)=AAABR(74)
      MTR022(1,1)=AAABR(80)
      MTR022(1,2)=AAABR(79)
      MTR022(1,3)=AAABR(83)
      MTR022(2,1)=AAABR(81)
      MTR022(2,2)=AAABR(84)
      MTR022(2,3)=AAABR(86)
      MTR022(3,1)=AAABR(82)
      MTR022(3,2)=AAABR(78)
      MTR022(3,3)=AAABR(85)
      MTR023(1,1)=AAABR(90)
      MTR023(1,2)=AAABR(94)
      MTR023(1,3)=AAABR(88)
      MTR023(2,1)=AAABR(87)
      MTR023(2,2)=AAABR(95)
      MTR023(2,3)=AAABR(92)
      MTR023(3,1)=AAABR(91)
      MTR023(3,2)=AAABR(93)
      MTR023(3,3)=AAABR(89)
      MTR024(1)=AAABR(97)
      MTR024(2)=AAABR(96)
      MTR024(3)=AAABR(98)
      MTR025(1)=AAABR(100)
      MTR025(2)=AAABR(101)
      MTR025(3)=AAABR(99)
      MTR026(1,1)=AAABR(110)
      MTR026(1,2)=AAABR(102)
      MTR026(1,3)=AAABR(107)
      MTR026(2,1)=AAABR(108)
      MTR026(2,2)=AAABR(109)
      MTR026(2,3)=AAABR(106)
      MTR026(3,1)=AAABR(104)
      MTR026(3,2)=AAABR(105)
      MTR026(3,3)=AAABR(103)
      MTR027(1)=AAABR(113)
      MTR027(2)=AAABR(111)
      MTR027(3)=AAABR(112)
      MTR028(1,1)=AAABR(116)
      MTR028(1,2)=AAABR(118)
      MTR028(1,3)=AAABR(119)
      MTR028(2,1)=AAABR(114)
      MTR028(2,2)=AAABR(115)
      MTR028(2,3)=AAABR(117)
      MTR029(1,1)=AAABR(121)
      MTR029(1,2)=AAABR(120)
      MTR029(1,3)=AAABR(122)
      MTR029(2,1)=AAABR(127)
      MTR029(2,2)=AAABR(124)
      MTR029(2,3)=AAABR(125)
      MTR029(3,1)=AAABR(126)
      MTR029(3,2)=AAABR(128)
      MTR029(3,3)=AAABR(123)
      MTR030(1,1)=AAABR(135)
      MTR030(1,2)=AAABR(130)
      MTR030(1,3)=AAABR(133)
      MTR030(2,1)=AAABR(131)
      MTR030(2,2)=AAABR(132)
      MTR030(2,3)=AAABR(134)
      MTR030(3,1)=AAABR(137)
      MTR030(3,2)=AAABR(136)
      MTR030(3,3)=AAABR(129)
      MTR031(1,1)=AAABR(138)
      MTR031(1,2)=AAABR(140)
      MTR031(1,3)=AAABR(142)
      MTR031(2,1)=AAABR(143)
      MTR031(2,2)=AAABR(145)
      MTR031(2,3)=AAABR(141)
      MTR031(3,1)=AAABR(139)
      MTR031(3,2)=AAABR(144)
      MTR031(3,3)=AAABR(146)
      MTR032(1)=AAABR(149)
      MTR032(2)=AAABR(147)
      MTR032(3)=AAABR(148)
      MTR033(1,1)=AAABR(151)
      MTR033(1,2)=AAABR(153)
      MTR033(1,3)=AAABR(155)
      MTR033(2,1)=AAABR(154)
      MTR033(2,2)=AAABR(150)
      MTR033(2,3)=AAABR(152)
      MTR034(1,1,1)=2D0*AAABR(164)
      MTR034(1,1,2)=2D0*AAABR(156)
      MTR034(1,1,3)=2D0*AAABR(157)
      MTR034(1,2,1)=AAABR(163)
      MTR034(1,2,2)=AAABR(159)
      MTR034(1,2,3)=AAABR(161)
      MTR034(2,1,1)=MTR034(1,2,1)
      MTR034(2,1,2)=MTR034(1,2,2)
      MTR034(2,1,3)=MTR034(1,2,3)
      MTR034(2,2,1)=2D0*AAABR(158)
      MTR034(2,2,2)=2D0*AAABR(160)
      MTR034(2,2,3)=2D0*AAABR(162)
      MTR035(1,1,1)=6D0*AAABR(171)
      MTR035(1,1,2)=2D0*AAABR(165)
      MTR035(1,1,3)=2D0*AAABR(173)
      MTR035(1,2,1)=MTR035(1,1,2)
      MTR035(1,2,2)=2D0*AAABR(166)
      MTR035(1,2,3)=AAABR(167)
      MTR035(1,3,1)=MTR035(1,1,3)
      MTR035(1,3,2)=MTR035(1,2,3)
      MTR035(1,3,3)=2D0*AAABR(168)
      MTR035(2,1,1)=MTR035(1,1,2)
      MTR035(2,1,2)=MTR035(1,2,2)
      MTR035(2,1,3)=MTR035(1,2,3)
      MTR035(2,2,1)=MTR035(1,2,2)
      MTR035(2,2,2)=6D0*AAABR(174)
      MTR035(2,2,3)=2D0*AAABR(172)
      MTR035(2,3,1)=MTR035(1,2,3)
      MTR035(2,3,2)=MTR035(2,2,3)
      MTR035(2,3,3)=2D0*AAABR(169)
      MTR035(3,1,1)=MTR035(1,1,3)
      MTR035(3,1,2)=MTR035(1,2,3)
      MTR035(3,1,3)=MTR035(1,3,3)
      MTR035(3,2,1)=MTR035(1,2,3)
      MTR035(3,2,2)=MTR035(2,2,3)
      MTR035(3,2,3)=MTR035(2,3,3)
      MTR035(3,3,1)=MTR035(1,3,3)
      MTR035(3,3,2)=MTR035(2,3,3)
      MTR035(3,3,3)=6D0*AAABR(170)
      MTR036(1)=AAABR(177)
      MTR036(2)=AAABR(176)
      MTR037(1)=AAABR(180)
      MTR037(2)=AAABR(179)
      MTR037(3)=AAABR(178)
      MTR038(1)=AAABR(175)
      MTR038(2)=AAABR(175)
      MTR038(3)=AAABR(181)
      MTR039(1)=AAABR(184)
      MTR039(2)=AAABR(184)
      MTR039(3)=AAABR(183)
      MTR040(1)=2D0*AAABR(185)
      MTR040(2)=2D0*AAABR(185)
      MTR040(3)=AAABR(186)
      MTR041(1)=AAABR(189)
      MTR041(2)=AAABR(191)
      MTR041(3)=AAABR(190)
      MTR042(1)=AAABR(192)
      MTR042(2)=AAABR(193)
      MTR042(3)=AAABR(194)
      MTR043(1)=AAABR(196)
      MTR043(2)=AAABR(197)
      MTR043(3)=AAABR(195)
      MTR044(1)=AAABR(198)
      MTR044(2)=AAABR(199)
      MTR044(3)=AAABR(200)
      MTR045(1)=AAABR(203)
      MTR045(2)=AAABR(201)
      MTR045(3)=AAABR(202)
      MTR046(1)=AAABR(205)
      MTR046(2)=AAABR(204)
      MTR046(3)=AAABR(206)
      MTR047(1)=AAABR(209)
      MTR047(2)=AAABR(207)
      MTR047(3)=AAABR(208)
      MTR048(1)=AAABR(211)
      MTR048(2)=AAABR(210)
      MTR048(3)=AAABR(212)
      MTR049(1)=AAABR(215)
      MTR049(2)=AAABR(214)
      MTR049(3)=AAABR(213)
      MTR050(1)=AAABR(218)
      MTR050(2)=AAABR(217)
      MTR050(3)=AAABR(216)
      MTR051(1)=AAABR(220)
      MTR051(2)=AAABR(221)
      MTR051(3)=AAABR(222)
      MTR052(1)=AAABR(225)
      MTR052(2)=AAABR(223)
      MTR052(3)=AAABR(224)
      MTR053(1,1)=AAABR(226)
      MTR053(1,2)=AAABR(231)
      MTR053(1,3)=AAABR(227)
      MTR053(2,1)=AAABR(230)
      MTR053(2,2)=AAABR(228)
      MTR053(2,3)=AAABR(229)
      MTR054(1)=AAABR(233)
      MTR054(2)=AAABR(234)
      MTR054(3)=AAABR(235)
      MTR055(1,1)=AAABR(239)
      MTR055(1,2)=AAABR(238)
      MTR055(2,1)=AAABR(237)
      MTR055(2,2)=AAABR(236)
      MTR056(1,1)=AAABR(239)
      MTR056(1,2)=AAABR(237)
      MTR056(2,1)=AAABR(238)
      MTR056(2,2)=AAABR(236)
      MTR057(1,1,1)=AAABR(247)
      MTR057(1,1,2)=AAABR(240)
      MTR057(1,2,1)=AAABR(246)
      MTR057(1,2,2)=AAABR(241)
      MTR057(2,1,1)=AAABR(245)
      MTR057(2,1,2)=AAABR(242)
      MTR057(2,2,1)=AAABR(243)
      MTR057(2,2,2)=AAABR(244)
      MTR058(1,1,1)=AAABR(247)
      MTR058(1,1,2)=AAABR(240)
      MTR058(1,2,1)=AAABR(245)
      MTR058(1,2,2)=AAABR(242)
      MTR058(2,1,1)=AAABR(246)
      MTR058(2,1,2)=AAABR(241)
      MTR058(2,2,1)=AAABR(243)
      MTR058(2,2,2)=AAABR(244)
      MTR059(1,1,1)=AAABR(250)
      MTR059(1,1,2)=AAABR(254)
      MTR059(1,1,3)=AAABR(259)
      MTR059(1,2,1)=AAABR(255)
      MTR059(1,2,2)=AAABR(251)
      MTR059(1,2,3)=AAABR(257)
      MTR059(2,1,1)=AAABR(256)
      MTR059(2,1,2)=AAABR(252)
      MTR059(2,1,3)=AAABR(258)
      MTR059(2,2,1)=AAABR(248)
      MTR059(2,2,2)=AAABR(249)
      MTR059(2,2,3)=AAABR(253)
      MTR060(1,1,1)=AAABR(250)
      MTR060(1,1,2)=AAABR(254)
      MTR060(1,1,3)=AAABR(259)
      MTR060(1,2,1)=AAABR(256)
      MTR060(1,2,2)=AAABR(252)
      MTR060(1,2,3)=AAABR(258)
      MTR060(2,1,1)=AAABR(255)
      MTR060(2,1,2)=AAABR(251)
      MTR060(2,1,3)=AAABR(257)
      MTR060(2,2,1)=AAABR(248)
      MTR060(2,2,2)=AAABR(249)
      MTR060(2,2,3)=AAABR(253)
      MTR061(1,1)=2D0*AAABR(260)
      MTR061(1,2)=2D0*AAABR(260)
      MTR061(1,3)=AAABR(262)
      MTR061(2,1)=2D0*AAABR(261)
      MTR061(2,2)=2D0*AAABR(261)
      MTR061(2,3)=AAABR(263)
      MTR062(1,3)=AAABR(264)
      MTR062(2,3)=AAABR(265)
      MTR063(1,1)=AAABR(274)
      MTR063(1,2)=AAABR(266)
      MTR063(1,3)=AAABR(272)
      MTR063(1,4)=AAABR(270)
      MTR063(1,5)=AAABR(268)
      MTR063(2,1)=AAABR(280)
      MTR063(2,2)=AAABR(278)
      MTR063(2,3)=AAABR(276)
      MTR063(2,4)=AAABR(282)
      MTR063(2,5)=AAABR(284)
      MTR064(1,1)=AAABR(275)
      MTR064(1,2)=AAABR(267)
      MTR064(1,3)=AAABR(273)
      MTR064(1,4)=AAABR(271)
      MTR064(1,5)=AAABR(269)
      MTR064(2,1)=AAABR(281)
      MTR064(2,2)=AAABR(279)
      MTR064(2,3)=AAABR(277)
      MTR064(2,4)=AAABR(283)
      MTR064(2,5)=AAABR(285)
      MTR065(1,1)=AAABR(298)
      MTR065(1,2)=AAABR(302)
      MTR065(1,3)=AAABR(294)
      MTR065(1,4)=AAABR(292)
      MTR065(1,5)=AAABR(300)
      MTR065(2,1)=AAABR(296)
      MTR065(2,2)=AAABR(290)
      MTR065(2,3)=AAABR(288)
      MTR065(2,4)=AAABR(304)
      MTR065(2,5)=AAABR(286)
      MTR066(1,1)=AAABR(299)
      MTR066(1,2)=AAABR(303)
      MTR066(1,3)=AAABR(295)
      MTR066(1,4)=AAABR(293)
      MTR066(1,5)=AAABR(301)
      MTR066(2,1)=AAABR(297)
      MTR066(2,2)=AAABR(291)
      MTR066(2,3)=AAABR(289)
      MTR066(2,4)=AAABR(305)
      MTR066(2,5)=AAABR(287)
      MTR067(1,1)=AAABR(316)
      MTR067(1,2)=AAABR(308)
      MTR067(1,3)=AAABR(306)
      MTR067(2,1)=AAABR(314)
      MTR067(2,2)=AAABR(310)
      MTR067(2,3)=AAABR(312)
      MTR068(1,1)=AAABR(317)
      MTR068(1,2)=AAABR(309)
      MTR068(1,3)=AAABR(307)
      MTR068(2,1)=AAABR(315)
      MTR068(2,2)=AAABR(311)
      MTR068(2,3)=AAABR(313)
      MTR069(1,1)=AAABR(326)
      MTR069(1,2)=AAABR(320)
      MTR069(1,3)=AAABR(322)
      MTR069(2,1)=AAABR(324)
      MTR069(2,2)=AAABR(328)
      MTR069(2,3)=AAABR(318)
      MTR070(1,1)=AAABR(327)
      MTR070(1,2)=AAABR(321)
      MTR070(1,3)=AAABR(323)
      MTR070(2,1)=AAABR(325)
      MTR070(2,2)=AAABR(329)
      MTR070(2,3)=AAABR(319)
      MTR071(3,1)=AAABR(331)
      MTR071(3,2)=AAABR(332)
      MTR072(3,1)=AAABR(334)
      MTR072(3,2)=AAABR(335)
      MTR072(3,3)=AAABR(333)
      MTR073(1,1)=AAABR(339)
      MTR073(1,2)=AAABR(339)
      MTR073(1,3)=AAABR(343)
      MTR073(2,1)=AAABR(338)
      MTR073(2,2)=AAABR(338)
      MTR073(2,3)=AAABR(347)
      MTR073(3,1)=AAABR(346)
      MTR073(3,2)=AAABR(346)
      MTR073(3,3)=AAABR(349)
      MTR073(4,1)=AAABR(340)
      MTR073(4,2)=AAABR(340)
      MTR073(4,3)=AAABR(341)
      MTR073(5,1)=AAABR(345)
      MTR073(5,2)=AAABR(345)
      MTR073(5,3)=AAABR(351)
      MTR074(1,3)=AAABR(344)
      MTR074(2,3)=AAABR(348)
      MTR074(3,3)=AAABR(350)
      MTR074(4,3)=AAABR(342)
      MTR074(5,3)=AAABR(352)
      MTR075(1,1)=2D0*AAABR(353)
      MTR075(1,2)=2D0*AAABR(353)
      MTR075(1,3)=AAABR(367)
      MTR075(2,1)=2D0*AAABR(363)
      MTR075(2,2)=2D0*AAABR(363)
      MTR075(2,3)=AAABR(362)
      MTR075(3,1)=2D0*AAABR(356)
      MTR075(3,2)=2D0*AAABR(356)
      MTR075(3,3)=AAABR(359)
      MTR075(4,1)=2D0*AAABR(357)
      MTR075(4,2)=2D0*AAABR(357)
      MTR075(4,3)=AAABR(355)
      MTR075(5,1)=2D0*AAABR(360)
      MTR075(5,2)=2D0*AAABR(360)
      MTR075(5,3)=AAABR(365)
      MTR076(1,3)=AAABR(366)
      MTR076(2,3)=AAABR(361)
      MTR076(3,3)=AAABR(358)
      MTR076(4,3)=AAABR(354)
      MTR076(5,3)=AAABR(364)
      MTR077(1,1)=AAABR(368)
      MTR077(1,2)=AAABR(368)
      MTR077(1,3)=AAABR(368)
      MTR077(2,1)=AAABR(370)
      MTR077(2,2)=AAABR(370)
      MTR077(2,3)=AAABR(370)
      MTR078(1,3)=AAABR(369)
      MTR078(2,3)=AAABR(371)
      MTR079(1,1)=AAABR(376)
      MTR079(1,2)=AAABR(376)
      MTR079(1,3)=AAABR(376)
      MTR079(2,1)=AAABR(375)
      MTR079(2,2)=AAABR(375)
      MTR079(2,3)=AAABR(375)
      MTR079(3,1)=AAABR(374)
      MTR079(3,2)=AAABR(374)
      MTR079(3,3)=AAABR(374)
      MTR079(4,1)=AAABR(373)
      MTR079(4,2)=AAABR(373)
      MTR079(4,3)=AAABR(373)
      MTR079(5,1)=AAABR(372)
      MTR079(5,2)=AAABR(372)
      MTR079(5,3)=AAABR(372)
      MTR080(1,1)=AAABR(397)
      MTR080(1,2)=AAABR(385)
      MTR080(1,3)=AAABR(381)
      MTR080(2,1)=AAABR(387)
      MTR080(2,2)=AAABR(377)
      MTR080(2,3)=AAABR(401)
      MTR080(3,1)=AAABR(399)
      MTR080(3,2)=AAABR(395)
      MTR080(3,3)=AAABR(389)
      MTR080(4,1)=AAABR(379)
      MTR080(4,2)=AAABR(393)
      MTR080(4,3)=AAABR(383)
      MTR080(5,1)=AAABR(391)
      MTR080(5,2)=AAABR(405)
      MTR080(5,3)=AAABR(403)
      MTR081(1,1)=AAABR(398)
      MTR081(1,2)=AAABR(386)
      MTR081(1,3)=AAABR(382)
      MTR081(2,1)=AAABR(388)
      MTR081(2,2)=AAABR(378)
      MTR081(2,3)=AAABR(402)
      MTR081(3,1)=AAABR(400)
      MTR081(3,2)=AAABR(396)
      MTR081(3,3)=AAABR(390)
      MTR081(4,1)=AAABR(380)
      MTR081(4,2)=AAABR(394)
      MTR081(4,3)=AAABR(384)
      MTR081(5,1)=AAABR(392)
      MTR081(5,2)=AAABR(406)
      MTR081(5,3)=AAABR(404)
      MTR082(1,1)=AAABR(429)
      MTR082(1,2)=AAABR(431)
      MTR082(1,3)=AAABR(419)
      MTR082(2,1)=AAABR(415)
      MTR082(2,2)=AAABR(427)
      MTR082(2,3)=AAABR(421)
      MTR082(3,1)=AAABR(413)
      MTR082(3,2)=AAABR(423)
      MTR082(3,3)=AAABR(425)
      MTR082(4,1)=AAABR(433)
      MTR082(4,2)=AAABR(409)
      MTR082(4,3)=AAABR(411)
      MTR082(5,1)=AAABR(435)
      MTR082(5,2)=AAABR(417)
      MTR082(5,3)=AAABR(407)
      MTR083(1,1)=AAABR(430)
      MTR083(1,2)=AAABR(432)
      MTR083(1,3)=AAABR(420)
      MTR083(2,1)=AAABR(416)
      MTR083(2,2)=AAABR(428)
      MTR083(2,3)=AAABR(422)
      MTR083(3,1)=AAABR(414)
      MTR083(3,2)=AAABR(424)
      MTR083(3,3)=AAABR(426)
      MTR083(4,1)=AAABR(434)
      MTR083(4,2)=AAABR(410)
      MTR083(4,3)=AAABR(412)
      MTR083(5,1)=AAABR(436)
      MTR083(5,2)=AAABR(418)
      MTR083(5,3)=AAABR(408)
      MTR084(1)=AAABR(437)
      MTR084(2)=AAABR(437)
      MTR084(3)=AAABR(438)
      MTR085(1,1)=AAABR(440)
      MTR085(1,2)=AAABR(439)
      MTR085(2,1)=AAABR(440)
      MTR085(2,2)=AAABR(439)
      MTR085(3,1)=AAABR(442)
      MTR085(3,2)=AAABR(441)
      MTR086(1,1)=AAABR(444)
      MTR086(1,2)=AAABR(443)
      MTR086(1,3)=AAABR(445)
      MTR086(2,1)=AAABR(444)
      MTR086(2,2)=AAABR(443)
      MTR086(2,3)=AAABR(445)
      MTR086(3,1)=AAABR(446)
      MTR086(3,2)=AAABR(448)
      MTR086(3,3)=AAABR(447)
      MTR087(1)=AAABR(449)
      MTR087(2)=AAABR(451)
      MTR087(3)=AAABR(453)
      MTR088(1)=AAABR(450)
      MTR088(2)=AAABR(450)
      MTR088(3)=AAABR(452)
      MTR089(1)=AAABR(455)
      MTR089(2)=AAABR(454)
      MTR089(3)=AAABR(457)
      MTR090(1)=AAABR(455)
      MTR090(2)=AAABR(455)
      MTR090(3)=AAABR(456)
      MTR091(1,1)=AAABR(460)
      MTR091(1,2)=AAABR(462)
      MTR091(1,3)=AAABR(458)
      MTR091(2,1)=AAABR(474)
      MTR091(2,2)=AAABR(472)
      MTR091(2,3)=AAABR(470)
      MTR092(1,1)=AAABR(461)
      MTR092(1,2)=AAABR(463)
      MTR092(1,3)=AAABR(459)
      MTR092(2,1)=AAABR(475)
      MTR092(2,2)=AAABR(473)
      MTR092(2,3)=AAABR(471)
      MTR093(1,1)=AAABR(466)
      MTR093(1,2)=AAABR(464)
      MTR093(1,3)=AAABR(468)
      MTR093(2,1)=AAABR(478)
      MTR093(2,2)=AAABR(480)
      MTR093(2,3)=AAABR(476)
      MTR094(1,1)=AAABR(467)
      MTR094(1,2)=AAABR(465)
      MTR094(1,3)=AAABR(469)
      MTR094(2,1)=AAABR(479)
      MTR094(2,2)=AAABR(481)
      MTR094(2,3)=AAABR(477)
      MTR095(1)=AAABR(482)
      MTR095(2)=AAABR(484)
      MTR095(3)=AAABR(486)
      MTR096(1)=AAABR(483)
      MTR096(2)=AAABR(485)
      MTR096(3)=AAABR(487)
      MTR097(1)=AAABR(492)
      MTR097(2)=AAABR(488)
      MTR097(3)=AAABR(490)
      MTR098(1)=AAABR(493)
      MTR098(2)=AAABR(489)
      MTR098(3)=AAABR(491)
      MTR099(1,1)=AAABR(496)
      MTR099(1,2)=AAABR(508)
      MTR099(1,3)=AAABR(494)
      MTR099(2,1)=AAABR(522)
      MTR099(2,2)=AAABR(498)
      MTR099(2,3)=AAABR(502)
      MTR099(3,1)=AAABR(520)
      MTR099(3,2)=AAABR(516)
      MTR099(3,3)=AAABR(510)
      MTR099(4,1)=AAABR(518)
      MTR099(4,2)=AAABR(504)
      MTR099(4,3)=AAABR(512)
      MTR099(5,1)=AAABR(514)
      MTR099(5,2)=AAABR(506)
      MTR099(5,3)=AAABR(500)
      MTR100(1,1)=AAABR(497)
      MTR100(1,2)=AAABR(509)
      MTR100(1,3)=AAABR(495)
      MTR100(2,1)=AAABR(523)
      MTR100(2,2)=AAABR(499)
      MTR100(2,3)=AAABR(503)
      MTR100(3,1)=AAABR(521)
      MTR100(3,2)=AAABR(517)
      MTR100(3,3)=AAABR(511)
      MTR100(4,1)=AAABR(519)
      MTR100(4,2)=AAABR(505)
      MTR100(4,3)=AAABR(513)
      MTR100(5,1)=AAABR(515)
      MTR100(5,2)=AAABR(507)
      MTR100(5,3)=AAABR(501)
      MTR101(1,1)=AAABR(524)
      MTR101(1,2)=AAABR(542)
      MTR101(1,3)=AAABR(544)
      MTR101(2,1)=AAABR(532)
      MTR101(2,2)=AAABR(526)
      MTR101(2,3)=AAABR(548)
      MTR101(3,1)=AAABR(550)
      MTR101(3,2)=AAABR(536)
      MTR101(3,3)=AAABR(538)
      MTR101(4,1)=AAABR(546)
      MTR101(4,2)=AAABR(534)
      MTR101(4,3)=AAABR(552)
      MTR101(5,1)=AAABR(530)
      MTR101(5,2)=AAABR(528)
      MTR101(5,3)=AAABR(540)
      MTR102(1,1)=AAABR(525)
      MTR102(1,2)=AAABR(543)
      MTR102(1,3)=AAABR(545)
      MTR102(2,1)=AAABR(533)
      MTR102(2,2)=AAABR(527)
      MTR102(2,3)=AAABR(549)
      MTR102(3,1)=AAABR(551)
      MTR102(3,2)=AAABR(537)
      MTR102(3,3)=AAABR(539)
      MTR102(4,1)=AAABR(547)
      MTR102(4,2)=AAABR(535)
      MTR102(4,3)=AAABR(553)
      MTR102(5,1)=AAABR(531)
      MTR102(5,2)=AAABR(529)
      MTR102(5,3)=AAABR(541)
      MTR103(1)=AAABR(437)
      MTR103(2)=AAABR(555)
      MTR103(3)=AAABR(554)
      MTR104(1,1)=AAABR(561)
      MTR104(1,2)=AAABR(558)
      MTR104(2,1)=AAABR(559)
      MTR104(2,2)=AAABR(560)
      MTR104(3,1)=AAABR(556)
      MTR104(3,2)=AAABR(557)
      MTR105(1,1)=AAABR(564)
      MTR105(1,2)=AAABR(569)
      MTR105(1,3)=AAABR(563)
      MTR105(2,1)=AAABR(565)
      MTR105(2,2)=AAABR(566)
      MTR105(2,3)=AAABR(568)
      MTR105(3,1)=AAABR(570)
      MTR105(3,2)=AAABR(567)
      MTR105(3,3)=AAABR(562)
      MTR106(1)=AAABR(573)
      MTR106(2)=AAABR(575)
      MTR106(3)=AAABR(571)
      MTR107(1)=AAABR(574)
      MTR107(2)=AAABR(576)
      MTR107(3)=AAABR(572)
      MTR108(1)=AAABR(579)
      MTR108(2)=AAABR(581)
      MTR108(3)=AAABR(577)
      MTR109(1)=AAABR(580)
      MTR109(2)=AAABR(582)
      MTR109(3)=AAABR(578)
      MTR110(1,1)=2D0*AAABR(583)
      MTR110(1,2)=AAABR(586)
      MTR110(1,3)=AAABR(587)
      MTR110(1,4)=AAABR(584)
      MTR110(1,5)=AAABR(585)
      MTR110(2,1)=AAABR(586)
      MTR110(2,2)=2D0*AAABR(615)
      MTR110(2,3)=AAABR(613)
      MTR110(2,4)=AAABR(614)
      MTR110(2,5)=AAABR(616)
      MTR110(3,1)=AAABR(587)
      MTR110(3,2)=AAABR(613)
      MTR110(3,3)=2D0*AAABR(637)
      MTR110(3,4)=AAABR(638)
      MTR110(3,5)=AAABR(639)
      MTR110(4,1)=AAABR(584)
      MTR110(4,2)=AAABR(614)
      MTR110(4,3)=AAABR(638)
      MTR110(4,4)=2D0*AAABR(656)
      MTR110(4,5)=AAABR(655)
      MTR110(5,1)=AAABR(585)
      MTR110(5,2)=AAABR(616)
      MTR110(5,3)=AAABR(639)
      MTR110(5,4)=AAABR(655)
      MTR110(5,5)=2D0*AAABR(667)
      MTR111(1,1,1)=AAABR(592)
      MTR111(1,1,2)=AAABR(594)
      MTR111(1,2,1)=AAABR(597)
      MTR111(1,2,2)=AAABR(589)
      MTR111(1,3,1)=AAABR(591)
      MTR111(1,3,2)=AAABR(595)
      MTR111(1,4,1)=AAABR(593)
      MTR111(1,4,2)=AAABR(596)
      MTR111(1,5,1)=AAABR(588)
      MTR111(1,5,2)=AAABR(590)
      MTR111(2,1,1)=AAABR(597)
      MTR111(2,1,2)=AAABR(589)
      MTR111(2,2,1)=AAABR(621)
      MTR111(2,2,2)=AAABR(620)
      MTR111(2,3,1)=AAABR(617)
      MTR111(2,3,2)=AAABR(622)
      MTR111(2,4,1)=AAABR(624)
      MTR111(2,4,2)=AAABR(618)
      MTR111(2,5,1)=AAABR(619)
      MTR111(2,5,2)=AAABR(623)
      MTR111(3,1,1)=AAABR(591)
      MTR111(3,1,2)=AAABR(595)
      MTR111(3,2,1)=AAABR(617)
      MTR111(3,2,2)=AAABR(622)
      MTR111(3,3,1)=AAABR(641)
      MTR111(3,3,2)=AAABR(643)
      MTR111(3,4,1)=AAABR(642)
      MTR111(3,4,2)=AAABR(645)
      MTR111(3,5,1)=AAABR(640)
      MTR111(3,5,2)=AAABR(644)
      MTR111(4,1,1)=AAABR(593)
      MTR111(4,1,2)=AAABR(596)
      MTR111(4,2,1)=AAABR(624)
      MTR111(4,2,2)=AAABR(618)
      MTR111(4,3,1)=AAABR(642)
      MTR111(4,3,2)=AAABR(645)
      MTR111(4,4,1)=AAABR(660)
      MTR111(4,4,2)=AAABR(659)
      MTR111(4,5,1)=AAABR(658)
      MTR111(4,5,2)=AAABR(657)
      MTR111(5,1,1)=AAABR(588)
      MTR111(5,1,2)=AAABR(590)
      MTR111(5,2,1)=AAABR(619)
      MTR111(5,2,2)=AAABR(623)
      MTR111(5,3,1)=AAABR(640)
      MTR111(5,3,2)=AAABR(644)
      MTR111(5,4,1)=AAABR(658)
      MTR111(5,4,2)=AAABR(657)
      MTR111(5,5,1)=AAABR(668)
      MTR111(5,5,2)=AAABR(669)
      MTR112(1,1,1)=AAABR(610)
      MTR112(1,1,2)=AAABR(612)
      MTR112(1,1,3)=AAABR(607)
      MTR112(1,2,1)=AAABR(605)
      MTR112(1,2,2)=AAABR(603)
      MTR112(1,2,3)=AAABR(608)
      MTR112(1,3,1)=AAABR(600)
      MTR112(1,3,2)=AAABR(609)
      MTR112(1,3,3)=AAABR(599)
      MTR112(1,4,1)=AAABR(611)
      MTR112(1,4,2)=AAABR(602)
      MTR112(1,4,3)=AAABR(604)
      MTR112(1,5,1)=AAABR(606)
      MTR112(1,5,2)=AAABR(601)
      MTR112(1,5,3)=AAABR(598)
      MTR112(2,1,1)=AAABR(605)
      MTR112(2,1,2)=AAABR(603)
      MTR112(2,1,3)=AAABR(608)
      MTR112(2,2,1)=AAABR(634)
      MTR112(2,2,2)=AAABR(633)
      MTR112(2,2,3)=AAABR(632)
      MTR112(2,3,1)=AAABR(631)
      MTR112(2,3,2)=AAABR(626)
      MTR112(2,3,3)=AAABR(635)
      MTR112(2,4,1)=AAABR(629)
      MTR112(2,4,2)=AAABR(630)
      MTR112(2,4,3)=AAABR(625)
      MTR112(2,5,1)=AAABR(636)
      MTR112(2,5,2)=AAABR(628)
      MTR112(2,5,3)=AAABR(627)
      MTR112(3,1,1)=AAABR(600)
      MTR112(3,1,2)=AAABR(609)
      MTR112(3,1,3)=AAABR(599)
      MTR112(3,2,1)=AAABR(631)
      MTR112(3,2,2)=AAABR(626)
      MTR112(3,2,3)=AAABR(635)
      MTR112(3,3,1)=AAABR(648)
      MTR112(3,3,2)=AAABR(646)
      MTR112(3,3,3)=AAABR(653)
      MTR112(3,4,1)=AAABR(652)
      MTR112(3,4,2)=AAABR(654)
      MTR112(3,4,3)=AAABR(650)
      MTR112(3,5,1)=AAABR(649)
      MTR112(3,5,2)=AAABR(647)
      MTR112(3,5,3)=AAABR(651)
      MTR112(4,1,1)=AAABR(611)
      MTR112(4,1,2)=AAABR(602)
      MTR112(4,1,3)=AAABR(604)
      MTR112(4,2,1)=AAABR(629)
      MTR112(4,2,2)=AAABR(630)
      MTR112(4,2,3)=AAABR(625)
      MTR112(4,3,1)=AAABR(652)
      MTR112(4,3,2)=AAABR(654)
      MTR112(4,3,3)=AAABR(650)
      MTR112(4,4,1)=AAABR(665)
      MTR112(4,4,2)=AAABR(661)
      MTR112(4,4,3)=AAABR(666)
      MTR112(4,5,1)=AAABR(664)
      MTR112(4,5,2)=AAABR(663)
      MTR112(4,5,3)=AAABR(662)
      MTR112(5,1,1)=AAABR(606)
      MTR112(5,1,2)=AAABR(601)
      MTR112(5,1,3)=AAABR(598)
      MTR112(5,2,1)=AAABR(636)
      MTR112(5,2,2)=AAABR(628)
      MTR112(5,2,3)=AAABR(627)
      MTR112(5,3,1)=AAABR(649)
      MTR112(5,3,2)=AAABR(647)
      MTR112(5,3,3)=AAABR(651)
      MTR112(5,4,1)=AAABR(664)
      MTR112(5,4,2)=AAABR(663)
      MTR112(5,4,3)=AAABR(662)
      MTR112(5,5,1)=AAABR(671)
      MTR112(5,5,2)=AAABR(672)
      MTR112(5,5,3)=AAABR(670)
      MTR113(1,1)=AAABR(675)
      MTR113(1,2)=AAABR(673)
      MTR113(2,1)=AAABR(673)
      MTR113(2,2)=AAABR(677)
      MTR114(1,1)=AAABR(676)
      MTR114(1,2)=AAABR(674)
      MTR114(2,1)=AAABR(674)
      MTR114(2,2)=AAABR(678)
      MTR115(1,1)=AAABR(697)
      MTR115(1,2)=AAABR(687)
      MTR115(1,3)=AAABR(695)
      MTR115(1,4)=AAABR(691)
      MTR115(1,5)=AAABR(681)
      MTR115(2,1)=AAABR(685)
      MTR115(2,2)=AAABR(689)
      MTR115(2,3)=AAABR(683)
      MTR115(2,4)=AAABR(679)
      MTR115(2,5)=AAABR(693)
      MTR116(1,1)=AAABR(698)
      MTR116(1,2)=AAABR(688)
      MTR116(1,3)=AAABR(696)
      MTR116(1,4)=AAABR(692)
      MTR116(1,5)=AAABR(682)
      MTR116(2,1)=AAABR(686)
      MTR116(2,2)=AAABR(690)
      MTR116(2,3)=AAABR(684)
      MTR116(2,4)=AAABR(680)
      MTR116(2,5)=AAABR(694)
      MTR117(1,1)=AAABR(705)
      MTR117(1,2)=AAABR(701)
      MTR117(1,3)=AAABR(702)
      MTR117(1,4)=AAABR(704)
      MTR117(1,5)=AAABR(703)
      MTR117(2,1)=AAABR(701)
      MTR117(2,2)=AAABR(707)
      MTR117(2,3)=AAABR(709)
      MTR117(2,4)=AAABR(708)
      MTR117(2,5)=AAABR(706)
      MTR117(3,1)=AAABR(702)
      MTR117(3,2)=AAABR(709)
      MTR117(3,3)=AAABR(712)
      MTR117(3,4)=AAABR(710)
      MTR117(3,5)=AAABR(711)
      MTR117(4,1)=AAABR(704)
      MTR117(4,2)=AAABR(708)
      MTR117(4,3)=AAABR(710)
      MTR117(4,4)=AAABR(713)
      MTR117(4,5)=AAABR(714)
      MTR117(5,1)=AAABR(703)
      MTR117(5,2)=AAABR(706)
      MTR117(5,3)=AAABR(711)
      MTR117(5,4)=AAABR(714)
      MTR117(5,5)=AAABR(715)
      MTR118(1)=AAABR(719)
      MTR118(2)=AAABR(719)
      MTR118(3)=AAABR(720)
      MTR119(1)=2D0*AAABR(722)
      MTR119(2)=2D0*AAABR(722)
      MTR119(3)=-AAABR(723)
      MTR120(1)=AAABR(725)
      MTR120(2)=AAABR(725)
      MTR120(3)=AAABR(724)
      MTR121(1)=AAABR(727)
      MTR121(2)=AAABR(726)
      MTR121(3)=AAABR(728)
      MTR122(1)=AAABR(730)
      MTR122(2)=AAABR(729)
      MTR122(3)=AAABR(731)
      MTR123(1)=AAABR(733)
      MTR123(2)=AAABR(734)
      MTR123(3)=AAABR(732)
      MTR124(1)=AAABR(737)
      MTR124(2)=AAABR(735)
      MTR124(3)=AAABR(736)
      MTR125(1)=AAABR(738)
      MTR125(2)=AAABR(740)
      MTR125(3)=AAABR(739)
      MTR126(1)=AAABR(743)
      MTR126(2)=AAABR(741)
      MTR126(3)=AAABR(742)
      MTR127(1)=AAABR(747)
      MTR127(2)=AAABR(746)
      MTR128(1,1)=AAABR(750)
      MTR128(1,2)=AAABR(748)
      MTR128(2,1)=MTR128(1,2)
      MTR128(2,2)=AAABR(749)
      MTR129(1,1)=AAABR(754)
      MTR129(1,2)=AAABR(753)
      MTR129(1,3)=AAABR(752)
      MTR129(2,1)=MTR129(1,2)
      MTR129(2,2)=AAABR(756)
      MTR129(2,3)=AAABR(751)
      MTR129(3,1)=MTR129(1,3)
      MTR129(3,2)=MTR129(2,3)
      MTR129(3,3)=AAABR(755)
      MTR130(1)=AAABR(758)
      MTR130(2)=AAABR(758)
      MTR130(3)=AAABR(757)
      MTR131(1)=2D0*AAABR(760)
      MTR131(2)=2D0*AAABR(760)
      MTR131(3)=-AAABR(761)
      MTR132(1)=AAABR(763)
      MTR132(2)=AAABR(763)
      MTR132(3)=AAABR(762)
      MTR133(1)=2D0*AAABR(765)
      MTR133(2)=2D0*AAABR(765)
      MTR133(3)=AAABR(764)
      MTR134(1,1)=AAABR(768)
      MTR134(1,2)=AAABR(768)
      MTR134(1,3)=AAABR(769)
      MTR134(2,1)=AAABR(767)
      MTR134(2,2)=AAABR(767)
      MTR134(2,3)=AAABR(770)
      MTR135(1,3)=AAABR(771)
      MTR135(2,3)=AAABR(772)
      MTR136(1,1)=AAABR(775)
      MTR136(1,2)=AAABR(775)
      MTR136(1,3)=AAABR(773)
      MTR136(2,1)=AAABR(777)
      MTR136(2,2)=AAABR(777)
      MTR136(2,3)=AAABR(776)
      MTR136(3,1)=AAABR(774)
      MTR136(3,2)=AAABR(774)
      MTR136(3,3)=AAABR(778)
      MTR137(1,3)=AAABR(780)
      MTR137(2,3)=AAABR(779)
      MTR137(3,3)=AAABR(781)
      MTR138(1)=AAABR(783)
      MTR138(2)=AAABR(784)
      MTR138(3)=AAABR(782)
      MTR139(1)=AAABR(785)
      MTR139(2)=AAABR(787)
      MTR139(3)=AAABR(786)
      MTR140(1)=AAABR(788)
      MTR140(2)=AAABR(789)
      MTR140(3)=AAABR(790)
      MTR141(1)=AAABR(792)
      MTR141(2)=AAABR(791)
      MTR141(3)=AAABR(793)
      MTR142(1)=AAABR(795)
      MTR142(2)=AAABR(794)
      MTR142(3)=AAABR(796)
      MTR143(1)=AAABR(798)
      MTR143(2)=AAABR(797)
      MTR143(3)=AAABR(799)
      MTR144(1)=AAABR(802)
      MTR144(2)=AAABR(800)
      MTR144(3)=AAABR(801)
      MTR145(1,1)=AAABR(803)
      MTR145(1,2)=AAABR(807)
      MTR145(1,3)=AAABR(804)
      MTR145(2,1)=AAABR(805)
      MTR145(2,2)=AAABR(808)
      MTR145(2,3)=AAABR(806)
      MTR146(1,1)=AAABR(811)
      MTR146(1,2)=AAABR(809)
      MTR146(1,3)=AAABR(813)
      MTR146(2,1)=AAABR(814)
      MTR146(2,2)=AAABR(810)
      MTR146(2,3)=AAABR(812)
      MTR147(1,1)=AAABR(816)
      MTR147(1,2)=AAABR(819)
      MTR147(1,3)=AAABR(815)
      MTR147(2,1)=AAABR(820)
      MTR147(2,2)=AAABR(823)
      MTR147(2,3)=AAABR(822)
      MTR147(3,1)=AAABR(818)
      MTR147(3,2)=AAABR(817)
      MTR147(3,3)=AAABR(821)
      MTR148(1,1)=AAABR(825)
      MTR148(1,2)=AAABR(826)
      MTR148(1,3)=AAABR(831)
      MTR148(2,1)=AAABR(827)
      MTR148(2,2)=AAABR(824)
      MTR148(2,3)=AAABR(828)
      MTR148(3,1)=AAABR(832)
      MTR148(3,2)=AAABR(830)
      MTR148(3,3)=AAABR(829)
      MTR149(1)=AAABR(834)
      MTR149(2)=AAABR(833)
      MTR149(3)=AAABR(835)
      MTR150(1)=AAABR(837)
      MTR150(2)=AAABR(836)
      MTR150(3)=AAABR(838)
      MTR151(1)=AAABR(839)
      MTR151(2)=AAABR(841)
      MTR151(3)=AAABR(840)
      MTR152(1,1)=AAABR(845)
      MTR152(1,2)=AAABR(844)
      MTR152(1,3)=AAABR(842)
      MTR152(2,1)=AAABR(846)
      MTR152(2,2)=AAABR(847)
      MTR152(2,3)=AAABR(843)
      MTR153(1,1)=AAABR(848)
      MTR153(1,2)=AAABR(852)
      MTR153(1,3)=AAABR(851)
      MTR153(2,1)=AAABR(849)
      MTR153(2,2)=AAABR(853)
      MTR153(2,3)=AAABR(850)
      MTR154(1,1)=AAABR(861)
      MTR154(1,2)=AAABR(862)
      MTR154(1,3)=AAABR(860)
      MTR154(2,1)=AAABR(855)
      MTR154(2,2)=AAABR(857)
      MTR154(2,3)=AAABR(859)
      MTR154(3,1)=AAABR(858)
      MTR154(3,2)=AAABR(854)
      MTR154(3,3)=AAABR(856)
      MTR155(1,1)=AAABR(863)
      MTR155(1,2)=AAABR(871)
      MTR155(1,3)=AAABR(866)
      MTR155(2,1)=AAABR(869)
      MTR155(2,2)=AAABR(865)
      MTR155(2,3)=AAABR(867)
      MTR155(3,1)=AAABR(868)
      MTR155(3,2)=AAABR(870)
      MTR155(3,3)=AAABR(864)
      MTR156(1)=AAABR(874)
      MTR156(2)=AAABR(873)
      MTR157(1)=AAABR(875)
      MTR157(2)=AAABR(877)
      MTR157(3)=AAABR(876)
      MTR158(1,1)=AAABR(880)
      MTR158(1,2)=AAABR(878)
      MTR158(2,1)=MTR158(1,2)
      MTR158(2,2)=AAABR(879)
      MTR159(1,1)=AAABR(885)
      MTR159(1,2)=AAABR(886)
      MTR159(1,3)=AAABR(884)
      MTR159(2,1)=AAABR(881)
      MTR159(2,2)=AAABR(883)
      MTR159(2,3)=AAABR(882)
      MTR160(1,1)=2D0*AAABR(891)
      MTR160(1,2)=AAABR(887)
      MTR160(1,3)=AAABR(889)
      MTR160(2,1)=MTR160(1,2)
      MTR160(2,2)=2D0*AAABR(890)
      MTR160(2,3)=AAABR(892)
      MTR160(3,1)=MTR160(1,3)
      MTR160(3,2)=MTR160(2,3)
      MTR160(3,3)=2D0*AAABR(888)
      MTR161(1,2)=AAABR(894)
      MTR161(1,3)=AAABR(893)
      MTR161(2,1)=AAABR(894)
      MTR161(2,3)=AAABR(893)
      MTR161(3,1)=AAABR(893)
      MTR161(3,2)=AAABR(893)
      MTR162(1,1)=2D0*AAABR(894)
      MTR162(1,2)=AAABR(894)
      MTR162(1,3)=AAABR(893)
      MTR162(2,1)=AAABR(894)
      MTR162(2,2)=2D0*AAABR(894)
      MTR162(2,3)=AAABR(893)
      MTR162(3,1)=AAABR(893)
      MTR162(3,2)=AAABR(893)
      MTR162(3,3)=2D0*AAABR(895)
      MTR163(1,3)=1D0
      MTR163(2,3)=1D0
      MTR164(1,3)=AAABR(897)
      MTR164(2,3)=AAABR(897)
      MTR164(3,3)=2D0*AAABR(896)
      MTR165(3,3)=1D0
      MTR166(3,1)=AAABR(897)
      MTR166(3,2)=AAABR(897)
      MTR166(3,3)=2D0*AAABR(896)
      MTR167(1,1)=2D0*AAABR(900)
      MTR167(1,2)=2D0*AAABR(900)
      MTR167(1,3)=-AAABR(902)
      MTR167(2,1)=2D0*AAABR(900)
      MTR167(2,2)=2D0*AAABR(900)
      MTR167(2,3)=-AAABR(902)
      MTR167(3,1)=2D0*AAABR(899)
      MTR167(3,2)=2D0*AAABR(899)
      MTR167(3,3)=-2D0*AAABR(901)
      MTR168(3,1)=1D0
      MTR168(3,2)=1D0
      MTR169(1,3)=AAABR(903)
      MTR169(2,3)=AAABR(903)
      MTR169(3,3)=-AAABR(904)
      MTR170(1,2)=AAABR(906)
      MTR170(1,3)=AAABR(906)
      MTR170(2,1)=AAABR(906)
      MTR170(2,3)=AAABR(906)
      MTR170(3,1)=AAABR(908)
      MTR170(3,2)=AAABR(908)
      MTR171(1,1)=AAABR(894)
      MTR171(1,2)=2D0*AAABR(907)
      MTR171(1,3)=2D0*AAABR(905)
      MTR171(2,1)=2D0*AAABR(907)
      MTR171(2,2)=AAABR(894)
      MTR171(2,3)=2D0*AAABR(905)
      MTR171(3,1)=2D0*AAABR(905)
      MTR171(3,2)=2D0*AAABR(905)
      MTR171(3,3)=AAABR(909)
      MTR172(3,1)=2D0*AAABR(911)
      MTR172(3,2)=2D0*AAABR(911)
      MTR172(3,3)=AAABR(912)
      MTR173(1,1)=AAABR(913)
      MTR173(1,2)=AAABR(914)
      MTR173(1,3)=AAABR(915)
      MTR173(2,1)=AAABR(913)
      MTR173(2,2)=AAABR(914)
      MTR173(2,3)=AAABR(915)
      MTR173(3,1)=AAABR(917)
      MTR173(3,2)=AAABR(916)
      MTR173(3,3)=AAABR(918)
      MTR174(1,1)=AAABR(920)
      MTR174(1,2)=AAABR(924)
      MTR174(1,3)=AAABR(922)
      MTR174(2,1)=AAABR(920)
      MTR174(2,2)=AAABR(924)
      MTR174(2,3)=AAABR(922)
      MTR174(3,1)=AAABR(923)
      MTR174(3,2)=AAABR(919)
      MTR174(3,3)=AAABR(921)
      MTR175(3,1)=AAABR(927)
      MTR175(3,2)=AAABR(925)
      MTR175(3,3)=AAABR(926)
      MTR176(3,1)=AAABR(930)
      MTR176(3,2)=AAABR(928)
      MTR176(3,3)=AAABR(929)
      MTR177(1,1)=AAABR(932)
      MTR177(1,2)=AAABR(931)
      MTR177(1,3)=AAABR(935)
      MTR177(2,1)=AAABR(932)
      MTR177(2,2)=AAABR(931)
      MTR177(2,3)=AAABR(935)
      MTR177(3,1)=AAABR(936)
      MTR177(3,2)=AAABR(934)
      MTR177(3,3)=AAABR(933)
      MTR178(3,1)=AAABR(939)
      MTR178(3,2)=AAABR(938)
      MTR178(3,3)=AAABR(937)
      MTR179(3,1)=AAABR(942)
      MTR179(3,2)=AAABR(940)
      MTR179(3,3)=AAABR(941)
      MTR180(1,1)=AAABR(944)
      MTR180(1,2)=AAABR(943)
      MTR180(1,3)=AAABR(945)
      MTR180(2,1)=AAABR(944)
      MTR180(2,2)=AAABR(943)
      MTR180(2,3)=AAABR(945)
      MTR180(3,1)=AAABR(946)
      MTR180(3,2)=AAABR(947)
      MTR180(3,3)=AAABR(948)
      MTR181(1,1)=AAABR(950)
      MTR181(1,2)=AAABR(953)
      MTR181(1,3)=AAABR(951)
      MTR181(2,1)=AAABR(950)
      MTR181(2,2)=AAABR(953)
      MTR181(2,3)=AAABR(951)
      MTR181(3,1)=AAABR(952)
      MTR181(3,2)=AAABR(954)
      MTR181(3,3)=AAABR(949)
      MTR182(3,1)=AAABR(955)
      MTR182(3,2)=AAABR(957)
      MTR182(3,3)=AAABR(956)
      MTR183(3,1)=AAABR(960)
      MTR183(3,2)=AAABR(958)
      MTR183(3,3)=AAABR(959)
      MTR184(1,1)=AAABR(966)
      MTR184(1,2)=AAABR(962)
      MTR184(1,3)=AAABR(963)
      MTR184(2,1)=AAABR(966)
      MTR184(2,2)=AAABR(962)
      MTR184(2,3)=AAABR(963)
      MTR184(3,1)=AAABR(964)
      MTR184(3,2)=AAABR(965)
      MTR184(3,3)=AAABR(961)
      MTR185(1,1)=AAABR(967)
      MTR185(1,2)=AAABR(968)
      MTR185(1,3)=AAABR(969)
      MTR185(2,1)=AAABR(967)
      MTR185(2,2)=AAABR(968)
      MTR185(2,3)=AAABR(969)
      MTR185(3,1)=AAABR(971)
      MTR185(3,2)=AAABR(972)
      MTR185(3,3)=AAABR(970)
      MTR186(1,1)=AAABR(976)
      MTR186(1,2)=AAABR(974)
      MTR186(1,3)=AAABR(975)
      MTR186(2,1)=AAABR(976)
      MTR186(2,2)=AAABR(974)
      MTR186(2,3)=AAABR(975)
      MTR186(3,1)=AAABR(978)
      MTR186(3,2)=AAABR(977)
      MTR186(3,3)=AAABR(973)
      MTR187(3,1)=AAABR(981)
      MTR187(3,2)=AAABR(979)
      MTR187(3,3)=AAABR(980)
      MTR188(1,1)=AAABR(983)
      MTR188(1,2)=AAABR(985)
      MTR188(1,3)=AAABR(984)
      MTR188(2,1)=AAABR(983)
      MTR188(2,2)=AAABR(985)
      MTR188(2,3)=AAABR(984)
      MTR188(3,1)=AAABR(986)
      MTR188(3,2)=AAABR(982)
      MTR188(3,3)=AAABR(987)
      MTR189(1,1)=AAABR(988)
      MTR189(1,2)=AAABR(992)
      MTR189(1,3)=AAABR(991)
      MTR189(2,1)=AAABR(988)
      MTR189(2,2)=AAABR(992)
      MTR189(2,3)=AAABR(991)
      MTR189(3,1)=AAABR(989)
      MTR189(3,2)=AAABR(993)
      MTR189(3,3)=AAABR(990)
      MTR190(1)=AAABR(719)
      MTR190(2)=AAABR(719)
      MTR190(3)=AAABR(994)
      MTR191(1)=AAABR(996)
      MTR191(2)=AAABR(996)
      MTR191(3)=AAABR(997)
      MTR192(1,1)=2D0*AAABR(998)
      MTR192(1,2)=2D0*AAABR(998)
      MTR192(1,3)=AAABR(1000)
      MTR192(2,1)=2D0*AAABR(999)
      MTR192(2,2)=2D0*AAABR(999)
      MTR192(2,3)=AAABR(1001)
      MTR193(1,1)=AAABR(1002)
      MTR193(1,2)=AAABR(1002)
      MTR193(1,3)=AAABR(1007)
      MTR193(2,1)=AAABR(1003)
      MTR193(2,2)=AAABR(1003)
      MTR193(2,3)=AAABR(1005)
      MTR193(3,1)=AAABR(1004)
      MTR193(3,2)=AAABR(1004)
      MTR193(3,3)=AAABR(1006)
      MTR194(1)=AAABR(725)
      MTR194(2)=AAABR(725)
      MTR194(3)=AAABR(1008)
      MTR195(1,1)=AAABR(1010)
      MTR195(1,2)=AAABR(1010)
      MTR195(1,3)=AAABR(1013)
      MTR195(2,1)=AAABR(1011)
      MTR195(2,2)=AAABR(1011)
      MTR195(2,3)=AAABR(1012)
      MTR196(1,3)=AAABR(1014)
      MTR196(2,3)=AAABR(1015)
      MTR197(1,3)=AAABR(1018)
      MTR197(2,3)=AAABR(1017)
      MTR197(3,3)=AAABR(1016)
      MTR198(1,1,1)=AAABR(1020)
      MTR198(1,1,2)=AAABR(1020)
      MTR198(1,1,3)=AAABR(1023)
      MTR198(1,2,1)=AAABR(1021)
      MTR198(1,2,2)=AAABR(1021)
      MTR198(1,2,3)=AAABR(1022)
      MTR198(2,1,1)=MTR198(1,2,1)
      MTR198(2,1,2)=MTR198(1,2,2)
      MTR198(2,1,3)=MTR198(1,2,3)
      MTR198(2,2,1)=AAABR(1019)
      MTR198(2,2,2)=AAABR(1019)
      MTR198(2,2,3)=AAABR(1024)
      MTR199(1,1,3)=AAABR(1026)
      MTR199(1,2,3)=AAABR(1025)
      MTR199(2,1,3)=MTR199(1,2,3)
      MTR199(2,2,3)=AAABR(1027)
      MTR200(1,1,3)=AAABR(1028)
      MTR200(1,2,3)=AAABR(1032)
      MTR200(1,3,3)=AAABR(1029)
      MTR200(2,1,3)=AAABR(1031)
      MTR200(2,2,3)=AAABR(1030)
      MTR200(2,3,3)=AAABR(1033)
      MTR201(1,1,1)=AAABR(1043)
      MTR201(1,1,2)=AAABR(1043)
      MTR201(1,1,3)=AAABR(1045)
      MTR201(1,2,1)=AAABR(1040)
      MTR201(1,2,2)=AAABR(1040)
      MTR201(1,2,3)=AAABR(1041)
      MTR201(1,3,1)=AAABR(1042)
      MTR201(1,3,2)=AAABR(1042)
      MTR201(1,3,3)=AAABR(1044)
      MTR201(2,1,1)=MTR201(1,2,1)
      MTR201(2,1,2)=MTR201(1,2,2)
      MTR201(2,1,3)=MTR201(1,2,3)
      MTR201(2,2,1)=AAABR(1038)
      MTR201(2,2,2)=AAABR(1038)
      MTR201(2,2,3)=AAABR(1037)
      MTR201(2,3,1)=AAABR(1034)
      MTR201(2,3,2)=AAABR(1034)
      MTR201(2,3,3)=AAABR(1035)
      MTR201(3,1,1)=MTR201(1,3,1)
      MTR201(3,1,2)=MTR201(1,3,2)
      MTR201(3,1,3)=MTR201(1,3,3)
      MTR201(3,2,1)=MTR201(2,3,1)
      MTR201(3,2,2)=MTR201(2,3,2)
      MTR201(3,2,3)=MTR201(2,3,3)
      MTR201(3,3,1)=AAABR(1036)
      MTR201(3,3,2)=AAABR(1036)
      MTR201(3,3,3)=AAABR(1039)
      MTR202(1,1,3)=AAABR(1049)
      MTR202(1,2,3)=AAABR(1046)
      MTR202(1,3,3)=AAABR(1047)
      MTR202(2,1,3)=MTR202(1,2,3)
      MTR202(2,2,3)=AAABR(1050)
      MTR202(2,3,3)=AAABR(1051)
      MTR202(3,1,3)=MTR202(1,3,3)
      MTR202(3,2,3)=MTR202(2,3,3)
      MTR202(3,3,3)=AAABR(1048)
      MTR203(3,1)=AAABR(903)
      MTR203(3,2)=AAABR(903)
      MTR203(3,3)=-AAABR(904)
      MTR204(1,2)=2D0*AAABR(900)
      MTR204(1,3)=-AAABR(1053)
      MTR204(2,1)=2D0*AAABR(900)
      MTR204(2,3)=-AAABR(1053)
      MTR204(3,1)=-AAABR(1053)
      MTR204(3,2)=-AAABR(1053)
      MTR205(1,1)=4D0*AAABR(900)
      MTR205(1,2)=2D0*AAABR(900)
      MTR205(1,3)=-AAABR(1053)
      MTR205(2,1)=2D0*AAABR(900)
      MTR205(2,2)=4D0*AAABR(900)
      MTR205(2,3)=-AAABR(1053)
      MTR205(3,1)=-AAABR(1053)
      MTR205(3,2)=-AAABR(1053)
      MTR205(3,3)=AAABR(1052)
      MTR206(1,3)=2D0*AAABR(911)
      MTR206(2,3)=2D0*AAABR(911)
      MTR206(3,3)=AAABR(912)
      MTR207(1,1)=2D0*AAABR(900)
      MTR207(1,2)=2D0*AAABR(900)
      MTR207(1,3)=2D0*AAABR(900)
      MTR207(2,1)=2D0*AAABR(900)
      MTR207(2,2)=2D0*AAABR(900)
      MTR207(2,3)=2D0*AAABR(900)
      MTR207(3,1)=-AAABR(1055)
      MTR207(3,2)=-AAABR(1055)
      MTR207(3,3)=-AAABR(1054)
      MTR208(1,1)=2D0*AAABR(1056)
      MTR208(1,2)=2D0*AAABR(1060)
      MTR208(1,3)=2D0*AAABR(1057)
      MTR208(2,1)=2D0*AAABR(1056)
      MTR208(2,2)=2D0*AAABR(1060)
      MTR208(2,3)=2D0*AAABR(1057)
      MTR208(3,1)=-AAABR(1058)
      MTR208(3,2)=-AAABR(1061)
      MTR208(3,3)=-AAABR(1059)
      MTR209(1,1)=2D0*AAABR(1062)
      MTR209(1,2)=2D0*AAABR(1063)
      MTR209(1,3)=2D0*AAABR(1064)
      MTR209(2,1)=2D0*AAABR(1062)
      MTR209(2,2)=2D0*AAABR(1063)
      MTR209(2,3)=2D0*AAABR(1064)
      MTR209(3,1)=-AAABR(1066)
      MTR209(3,2)=-AAABR(1065)
      MTR209(3,3)=-AAABR(1067)
      MTR210(1,1)=2D0*AAABR(1070)
      MTR210(1,2)=2D0*AAABR(1071)
      MTR210(1,3)=2D0*AAABR(1069)
      MTR210(2,1)=2D0*AAABR(1070)
      MTR210(2,2)=2D0*AAABR(1071)
      MTR210(2,3)=2D0*AAABR(1069)
      MTR210(3,1)=-AAABR(1073)
      MTR210(3,2)=-AAABR(1068)
      MTR210(3,3)=-AAABR(1072)
      MTR211(1,1)=2D0*AAABR(1077)
      MTR211(1,2)=2D0*AAABR(1075)
      MTR211(1,3)=2D0*AAABR(1074)
      MTR211(2,1)=2D0*AAABR(1077)
      MTR211(2,2)=2D0*AAABR(1075)
      MTR211(2,3)=2D0*AAABR(1074)
      MTR211(3,1)=-AAABR(1076)
      MTR211(3,2)=-AAABR(1078)
      MTR211(3,3)=-AAABR(1079)
      MTR212(1,1)=2D0*AAABR(1082)
      MTR212(1,2)=2D0*AAABR(1081)
      MTR212(1,3)=2D0*AAABR(1080)
      MTR212(2,1)=2D0*AAABR(1082)
      MTR212(2,2)=2D0*AAABR(1081)
      MTR212(2,3)=2D0*AAABR(1080)
      MTR212(3,1)=-AAABR(1084)
      MTR212(3,2)=-AAABR(1083)
      MTR212(3,3)=-AAABR(1085)
      MTR213(3,1)=AAABR(1088)
      MTR213(3,2)=AAABR(1086)
      MTR213(3,3)=AAABR(1087)
      MTR214(3,1)=AAABR(1089)
      MTR214(3,2)=AAABR(1090)
      MTR214(3,3)=AAABR(1091)
      MTR215(1,1)=2D0*AAABR(1095)
      MTR215(1,2)=2D0*AAABR(1093)
      MTR215(1,3)=2D0*AAABR(1092)
      MTR215(2,1)=2D0*AAABR(1095)
      MTR215(2,2)=2D0*AAABR(1093)
      MTR215(2,3)=2D0*AAABR(1092)
      MTR215(3,1)=-AAABR(1097)
      MTR215(3,2)=-AAABR(1096)
      MTR215(3,3)=-AAABR(1094)
      MTR216(3,1)=AAABR(1100)
      MTR216(3,2)=AAABR(1098)
      MTR216(3,3)=AAABR(1099)
      MTR217(3,1)=AAABR(1102)
      MTR217(3,2)=AAABR(1103)
      MTR217(3,3)=AAABR(1101)
      MTR218(1)=2D0*AAABR(722)
      MTR218(2)=2D0*AAABR(722)
      MTR218(3)=-AAABR(1104)
      MTR219(1,3)=AAABR(1106)
      MTR219(2,3)=AAABR(1107)
      MTR220(1,3)=AAABR(1109)
      MTR220(2,3)=AAABR(1108)
      MTR220(3,3)=AAABR(1110)
      MTR221(1)=2D0*AAABR(722)
      MTR221(2)=2D0*AAABR(722)
      MTR221(3)=-AAABR(1111)
      MTR222(1,1)=2D0*AAABR(1112)
      MTR222(1,2)=2D0*AAABR(1112)
      MTR222(1,3)=-AAABR(1113)
      MTR222(2,1)=2D0*AAABR(1114)
      MTR222(2,2)=2D0*AAABR(1114)
      MTR222(2,3)=-AAABR(1115)
      MTR223(1,1,1)=2D0*AAABR(1119)
      MTR223(1,1,2)=2D0*AAABR(1119)
      MTR223(1,1,3)=-AAABR(1120)
      MTR223(1,2,1)=2D0*AAABR(1121)
      MTR223(1,2,2)=2D0*AAABR(1121)
      MTR223(1,2,3)=-AAABR(1117)
      MTR223(2,1,1)=MTR223(1,2,1)
      MTR223(2,1,2)=MTR223(1,2,2)
      MTR223(2,1,3)=MTR223(1,2,3)
      MTR223(2,2,1)=2D0*AAABR(1116)
      MTR223(2,2,2)=2D0*AAABR(1116)
      MTR223(2,2,3)=-AAABR(1118)
      MTR224(1,1,1)=2D0*AAABR(1127)
      MTR224(1,1,2)=2D0*AAABR(1127)
      MTR224(1,1,3)=-AAABR(1130)
      MTR224(1,2,1)=2D0*AAABR(1124)
      MTR224(1,2,2)=2D0*AAABR(1124)
      MTR224(1,2,3)=-AAABR(1131)
      MTR224(1,3,1)=2D0*AAABR(1125)
      MTR224(1,3,2)=2D0*AAABR(1125)
      MTR224(1,3,3)=-AAABR(1133)
      MTR224(2,1,1)=MTR224(1,2,1)
      MTR224(2,1,2)=MTR224(1,2,2)
      MTR224(2,1,3)=MTR224(1,2,3)
      MTR224(2,2,1)=2D0*AAABR(1126)
      MTR224(2,2,2)=2D0*AAABR(1126)
      MTR224(2,2,3)=-AAABR(1132)
      MTR224(2,3,1)=2D0*AAABR(1128)
      MTR224(2,3,2)=2D0*AAABR(1128)
      MTR224(2,3,3)=-AAABR(1122)
      MTR224(3,1,1)=MTR224(1,3,1)
      MTR224(3,1,2)=MTR224(1,3,2)
      MTR224(3,1,3)=MTR224(1,3,3)
      MTR224(3,2,1)=MTR224(2,3,1)
      MTR224(3,2,2)=MTR224(2,3,2)
      MTR224(3,2,3)=MTR224(2,3,3)
      MTR224(3,3,1)=2D0*AAABR(1129)
      MTR224(3,3,2)=2D0*AAABR(1129)
      MTR224(3,3,3)=-AAABR(1123)
      MTR225(1,2)=1D0
      MTR225(1,3)=1D0
      MTR225(2,1)=1D0
      MTR225(2,3)=1D0
      MTR225(3,1)=1D0
      MTR225(3,2)=1D0
      MTR226(1,1)=2D0
      MTR226(1,2)=1D0
      MTR226(1,3)=1D0
      MTR226(2,1)=1D0
      MTR226(2,2)=2D0
      MTR226(2,3)=1D0
      MTR226(3,1)=1D0
      MTR226(3,2)=1D0
      MTR226(3,3)=2D0
      MTR227(1,1)=AAABR(1135)
      MTR227(1,2)=AAABR(1134)
      MTR227(1,3)=AAABR(1136)
      MTR227(2,1)=AAABR(1135)
      MTR227(2,2)=AAABR(1134)
      MTR227(2,3)=AAABR(1136)
      MTR227(3,1)=AAABR(1135)
      MTR227(3,2)=AAABR(1134)
      MTR227(3,3)=AAABR(1136)
      MTR228(1,1)=AAABR(1138)
      MTR228(1,2)=AAABR(1139)
      MTR228(1,3)=AAABR(1137)
      MTR228(2,1)=AAABR(1138)
      MTR228(2,2)=AAABR(1139)
      MTR228(2,3)=AAABR(1137)
      MTR228(3,1)=AAABR(1138)
      MTR228(3,2)=AAABR(1139)
      MTR228(3,3)=AAABR(1137)
      MTR229(1,1)=AAABR(1142)
      MTR229(1,2)=AAABR(1140)
      MTR229(1,3)=AAABR(1141)
      MTR229(2,1)=AAABR(1142)
      MTR229(2,2)=AAABR(1140)
      MTR229(2,3)=AAABR(1141)
      MTR229(3,1)=AAABR(1142)
      MTR229(3,2)=AAABR(1140)
      MTR229(3,3)=AAABR(1141)
      MTR230(1,1)=AAABR(1145)
      MTR230(1,2)=AAABR(1144)
      MTR230(1,3)=AAABR(1143)
      MTR230(2,1)=AAABR(1145)
      MTR230(2,2)=AAABR(1144)
      MTR230(2,3)=AAABR(1143)
      MTR230(3,1)=AAABR(1145)
      MTR230(3,2)=AAABR(1144)
      MTR230(3,3)=AAABR(1143)
      MTR231(1,1)=AAABR(1148)
      MTR231(1,2)=AAABR(1147)
      MTR231(1,3)=AAABR(1146)
      MTR231(2,1)=AAABR(1148)
      MTR231(2,2)=AAABR(1147)
      MTR231(2,3)=AAABR(1146)
      MTR231(3,1)=AAABR(1148)
      MTR231(3,2)=AAABR(1147)
      MTR231(3,3)=AAABR(1146)
      MTR232(1,1)=AAABR(1151)
      MTR232(1,2)=AAABR(1150)
      MTR232(1,3)=AAABR(1149)
      MTR232(2,1)=AAABR(1151)
      MTR232(2,2)=AAABR(1150)
      MTR232(2,3)=AAABR(1149)
      MTR232(3,1)=AAABR(1151)
      MTR232(3,2)=AAABR(1150)
      MTR232(3,3)=AAABR(1149)
      MTR233(1)=AAABR(725)
      MTR233(2)=AAABR(725)
      MTR233(3)=AAABR(1152)
      MTR234(1,1)=AAABR(1153)
      MTR234(1,2)=AAABR(1153)
      MTR234(1,3)=AAABR(1153)
      MTR234(2,1)=AAABR(1154)
      MTR234(2,2)=AAABR(1154)
      MTR234(2,3)=AAABR(1154)
      MTR235(1,1,1)=AAABR(1157)
      MTR235(1,1,2)=AAABR(1157)
      MTR235(1,1,3)=AAABR(1157)
      MTR235(1,2,1)=AAABR(1155)
      MTR235(1,2,2)=AAABR(1155)
      MTR235(1,2,3)=AAABR(1155)
      MTR235(2,1,1)=MTR235(1,2,1)
      MTR235(2,1,2)=MTR235(1,2,2)
      MTR235(2,1,3)=MTR235(1,2,3)
      MTR235(2,2,1)=AAABR(1156)
      MTR235(2,2,2)=AAABR(1156)
      MTR235(2,2,3)=AAABR(1156)
      MTR236(1,1,1)=AAABR(1160)
      MTR236(1,1,2)=AAABR(1160)
      MTR236(1,1,3)=AAABR(1160)
      MTR236(1,2,1)=AAABR(1159)
      MTR236(1,2,2)=AAABR(1159)
      MTR236(1,2,3)=AAABR(1159)
      MTR236(1,3,1)=AAABR(1163)
      MTR236(1,3,2)=AAABR(1163)
      MTR236(1,3,3)=AAABR(1163)
      MTR236(2,1,1)=MTR236(1,2,1)
      MTR236(2,1,2)=MTR236(1,2,2)
      MTR236(2,1,3)=MTR236(1,2,3)
      MTR236(2,2,1)=AAABR(1158)
      MTR236(2,2,2)=AAABR(1158)
      MTR236(2,2,3)=AAABR(1158)
      MTR236(2,3,1)=AAABR(1161)
      MTR236(2,3,2)=AAABR(1161)
      MTR236(2,3,3)=AAABR(1161)
      MTR236(3,1,1)=MTR236(1,3,1)
      MTR236(3,1,2)=MTR236(1,3,2)
      MTR236(3,1,3)=MTR236(1,3,3)
      MTR236(3,2,1)=MTR236(2,3,1)
      MTR236(3,2,2)=MTR236(2,3,2)
      MTR236(3,2,3)=MTR236(2,3,3)
      MTR236(3,3,1)=AAABR(1162)
      MTR236(3,3,2)=AAABR(1162)
      MTR236(3,3,3)=AAABR(1162)
      MTR237(1,1)=AAABR(1173)
      MTR237(1,2)=AAABR(1171)
      MTR237(1,3)=AAABR(1164)
      MTR237(2,1)=AAABR(1171)
      MTR237(2,2)=AAABR(1169)
      MTR237(2,3)=AAABR(1170)
      MTR237(3,1)=AAABR(1164)
      MTR237(3,2)=AAABR(1170)
      MTR237(3,3)=AAABR(1167)
      MTR238(1,1)=AAABR(1172)
      MTR238(1,2)=AAABR(1165)
      MTR238(1,3)=AAABR(1165)
      MTR238(2,1)=AAABR(1165)
      MTR238(2,2)=AAABR(1168)
      MTR238(2,3)=AAABR(1165)
      MTR238(3,1)=AAABR(1165)
      MTR238(3,2)=AAABR(1165)
      MTR238(3,3)=AAABR(1166)
      MTR239(1,1)=AAABR(1185)
      MTR239(1,2)=AAABR(1182)
      MTR239(1,3)=AAABR(1174)
      MTR239(2,1)=AAABR(1177)
      MTR239(2,2)=AAABR(1179)
      MTR239(2,3)=AAABR(1180)
      MTR239(3,1)=AAABR(1181)
      MTR239(3,2)=AAABR(1183)
      MTR239(3,3)=AAABR(1176)
      MTR240(1,1)=AAABR(1184)
      MTR240(2,2)=AAABR(1178)
      MTR240(3,3)=AAABR(1175)
      MTR241(1,1)=AAABR(1189)
      MTR241(2,2)=AAABR(1193)
      MTR241(3,3)=AAABR(1186)
      MTR242(1,1)=AAABR(1190)
      MTR242(1,2)=AAABR(1191)
      MTR242(1,3)=AAABR(1192)
      MTR242(2,1)=AAABR(1191)
      MTR242(2,2)=AAABR(1194)
      MTR242(2,3)=AAABR(1188)
      MTR242(3,1)=AAABR(1192)
      MTR242(3,2)=AAABR(1188)
      MTR242(3,3)=AAABR(1187)
      MTR243(1,1)=AAABR(1185)
      MTR243(1,2)=AAABR(1177)
      MTR243(1,3)=AAABR(1181)
      MTR243(2,1)=AAABR(1182)
      MTR243(2,2)=AAABR(1179)
      MTR243(2,3)=AAABR(1183)
      MTR243(3,1)=AAABR(1174)
      MTR243(3,2)=AAABR(1180)
      MTR243(3,3)=AAABR(1176)
      MTR244(1,1)=AAABR(1198)
      MTR244(1,2)=AAABR(1195)
      MTR244(1,3)=AAABR(1200)
      MTR244(2,1)=AAABR(1195)
      MTR244(2,2)=AAABR(1211)
      MTR244(2,3)=AAABR(1202)
      MTR244(3,1)=AAABR(1200)
      MTR244(3,2)=AAABR(1202)
      MTR244(3,3)=AAABR(1207)
      MTR245(1,1)=AAABR(1196)
      MTR245(1,2)=AAABR(1165)
      MTR245(1,3)=AAABR(1165)
      MTR245(2,1)=AAABR(1165)
      MTR245(2,2)=AAABR(1209)
      MTR245(2,3)=AAABR(1165)
      MTR245(3,1)=AAABR(1165)
      MTR245(3,2)=AAABR(1165)
      MTR245(3,3)=AAABR(1205)
      MTR246(1,1)=AAABR(1197)
      MTR246(1,2)=AAABR(1201)
      MTR246(1,3)=AAABR(1204)
      MTR246(2,1)=AAABR(1203)
      MTR246(2,2)=AAABR(1210)
      MTR246(2,3)=AAABR(1199)
      MTR246(3,1)=AAABR(1208)
      MTR246(3,2)=AAABR(1212)
      MTR246(3,3)=AAABR(1206)
      MTR247(1,1)=AAABR(1222)
      MTR247(2,2)=AAABR(1213)
      MTR247(3,3)=AAABR(1216)
      MTR248(1,1)=AAABR(1223)
      MTR248(1,2)=AAABR(1221)
      MTR248(1,3)=AAABR(1215)
      MTR248(2,1)=AAABR(1219)
      MTR248(2,2)=AAABR(1214)
      MTR248(2,3)=AAABR(1220)
      MTR248(3,1)=AAABR(1224)
      MTR248(3,2)=AAABR(1218)
      MTR248(3,3)=AAABR(1217)
      MTR249(1,1)=AAABR(1223)
      MTR249(1,2)=AAABR(1219)
      MTR249(1,3)=AAABR(1224)
      MTR249(2,1)=AAABR(1221)
      MTR249(2,2)=AAABR(1214)
      MTR249(2,3)=AAABR(1218)
      MTR249(3,1)=AAABR(1215)
      MTR249(3,2)=AAABR(1220)
      MTR249(3,3)=AAABR(1217)
      MTR250(1,1)=AAABR(1237)
      MTR250(1,2)=AAABR(1225)
      MTR250(1,3)=AAABR(1228)
      MTR250(2,1)=AAABR(1225)
      MTR250(2,2)=AAABR(1232)
      MTR250(2,3)=AAABR(1234)
      MTR250(3,1)=AAABR(1228)
      MTR250(3,2)=AAABR(1234)
      MTR250(3,3)=AAABR(1239)
      MTR251(1,1)=AAABR(1236)
      MTR251(1,2)=AAABR(1229)
      MTR251(1,3)=AAABR(1233)
      MTR251(2,1)=AAABR(1235)
      MTR251(2,2)=AAABR(1231)
      MTR251(2,3)=AAABR(1240)
      MTR251(3,1)=AAABR(1230)
      MTR251(3,2)=AAABR(1226)
      MTR251(3,3)=AAABR(1238)
      MTR252(1,1)=AAABR(1165)
      MTR252(1,2)=AAABR(1227)
      MTR252(1,3)=AAABR(1227)
      MTR252(2,1)=AAABR(1227)
      MTR252(2,2)=AAABR(1165)
      MTR252(2,3)=AAABR(1227)
      MTR252(3,1)=AAABR(1227)
      MTR252(3,2)=AAABR(1227)
      MTR252(3,3)=AAABR(1165)
      MTR253(1,1)=AAABR(1245)
      MTR253(1,2)=AAABR(1243)
      MTR253(1,3)=AAABR(1254)
      MTR253(2,1)=AAABR(1252)
      MTR253(2,2)=AAABR(1249)
      MTR253(2,3)=AAABR(1251)
      MTR253(3,1)=AAABR(1241)
      MTR253(3,2)=AAABR(1258)
      MTR253(3,3)=AAABR(1256)
      MTR254(1,1)=AAABR(1244)
      MTR254(1,2)=AAABR(1253)
      MTR254(1,3)=AAABR(1247)
      MTR254(2,1)=AAABR(1250)
      MTR254(2,2)=AAABR(1248)
      MTR254(2,3)=AAABR(1257)
      MTR254(3,1)=AAABR(1246)
      MTR254(3,2)=AAABR(1242)
      MTR254(3,3)=AAABR(1255)
      MTR255(1,1)=AAABR(1275)
      MTR255(1,2)=AAABR(1260)
      MTR255(1,3)=AAABR(1259)
      MTR255(2,1)=AAABR(1267)
      MTR255(2,2)=AAABR(1269)
      MTR255(2,3)=AAABR(1272)
      MTR255(3,1)=AAABR(1265)
      MTR255(3,2)=AAABR(1268)
      MTR255(3,3)=AAABR(1263)
      MTR256(1,1)=AAABR(1276)
      MTR256(1,2)=AAABR(1266)
      MTR256(1,3)=AAABR(1271)
      MTR256(2,1)=AAABR(1262)
      MTR256(2,2)=AAABR(1270)
      MTR256(2,3)=AAABR(1274)
      MTR256(3,1)=AAABR(1261)
      MTR256(3,2)=AAABR(1273)
      MTR256(3,3)=AAABR(1264)
      MTR257(1,1)=AAABR(1281)
      MTR257(1,2)=AAABR(1291)
      MTR257(1,3)=AAABR(1279)
      MTR257(2,1)=AAABR(1294)
      MTR257(2,2)=AAABR(1277)
      MTR257(2,3)=AAABR(1293)
      MTR257(3,1)=AAABR(1285)
      MTR257(3,2)=AAABR(1280)
      MTR257(3,3)=AAABR(1287)
      MTR258(1,1)=AAABR(1282)
      MTR258(1,2)=AAABR(1290)
      MTR258(1,3)=AAABR(1289)
      MTR258(2,1)=AAABR(1283)
      MTR258(2,2)=AAABR(1278)
      MTR258(2,3)=AAABR(1286)
      MTR258(3,1)=AAABR(1292)
      MTR258(3,2)=AAABR(1284)
      MTR258(3,3)=AAABR(1288)
      MTR259(1,1)=AAABR(1245)
      MTR259(1,2)=AAABR(1252)
      MTR259(1,3)=AAABR(1241)
      MTR259(2,1)=AAABR(1243)
      MTR259(2,2)=AAABR(1249)
      MTR259(2,3)=AAABR(1258)
      MTR259(3,1)=AAABR(1254)
      MTR259(3,2)=AAABR(1251)
      MTR259(3,3)=AAABR(1256)
      MTR260(1,1)=AAABR(1308)
      MTR260(1,2)=AAABR(1301)
      MTR260(1,3)=AAABR(1295)
      MTR260(2,1)=AAABR(1301)
      MTR260(2,2)=AAABR(1297)
      MTR260(2,3)=AAABR(1302)
      MTR260(3,1)=AAABR(1295)
      MTR260(3,2)=AAABR(1302)
      MTR260(3,3)=AAABR(1306)
      MTR261(1,1)=AAABR(1307)
      MTR261(1,2)=AAABR(1298)
      MTR261(1,3)=AAABR(1309)
      MTR261(2,1)=AAABR(1303)
      MTR261(2,2)=AAABR(1296)
      MTR261(2,3)=AAABR(1299)
      MTR261(3,1)=AAABR(1300)
      MTR261(3,2)=AAABR(1304)
      MTR261(3,3)=AAABR(1305)
      MTR262(1,1)=AAABR(1282)
      MTR262(1,2)=AAABR(1314)
      MTR262(1,3)=AAABR(1313)
      MTR262(2,1)=AAABR(1312)
      MTR262(2,2)=AAABR(1278)
      MTR262(2,3)=AAABR(1311)
      MTR262(3,1)=AAABR(1315)
      MTR262(3,2)=AAABR(1310)
      MTR262(3,3)=AAABR(1288)
      MTR263(1,1)=AAABR(1332)
      MTR263(1,2)=AAABR(1328)
      MTR263(1,3)=AAABR(1330)
      MTR263(2,1)=AAABR(1331)
      MTR263(2,2)=AAABR(1316)
      MTR263(2,3)=AAABR(1327)
      MTR263(3,1)=AAABR(1324)
      MTR263(3,2)=AAABR(1325)
      MTR263(3,3)=AAABR(1318)
      MTR264(1,1)=AAABR(1333)
      MTR264(1,2)=AAABR(1321)
      MTR264(1,3)=AAABR(1322)
      MTR264(2,1)=AAABR(1320)
      MTR264(2,2)=AAABR(1317)
      MTR264(2,3)=AAABR(1326)
      MTR264(3,1)=AAABR(1329)
      MTR264(3,2)=AAABR(1323)
      MTR264(3,3)=AAABR(1319)
      MTR265(1)=AAABR(1334)
      MTR265(2)=AAABR(1335)
      MTR265(3)=AAABR(1336)
      MTR266(1)=AAABR(1339)
      MTR266(2)=AAABR(1337)
      MTR266(3)=AAABR(1338)
      MTR267(1)=AAABR(1341)
      MTR267(2)=AAABR(1342)
      MTR267(3)=AAABR(1340)
      MTR268(1)=AAABR(1343)
      MTR268(2)=AAABR(1344)
      MTR268(3)=AAABR(1345)
      MTR269(1,1)=AAABR(1350)
      MTR269(1,2)=AAABR(1351)
      MTR269(1,3)=AAABR(1347)
      MTR269(2,1)=AAABR(1348)
      MTR269(2,2)=AAABR(1346)
      MTR269(2,3)=AAABR(1349)
      MTR270(1,1)=AAABR(1356)
      MTR270(1,2)=AAABR(1352)
      MTR270(1,3)=AAABR(1355)
      MTR270(2,1)=AAABR(1354)
      MTR270(2,2)=AAABR(1357)
      MTR270(2,3)=AAABR(1353)
      MTR271(1,1)=AAABR(1358)
      MTR271(1,2)=AAABR(1364)
      MTR271(1,3)=AAABR(1365)
      MTR271(2,1)=AAABR(1360)
      MTR271(2,2)=AAABR(1362)
      MTR271(2,3)=AAABR(1361)
      MTR271(3,1)=AAABR(1363)
      MTR271(3,2)=AAABR(1359)
      MTR271(3,3)=AAABR(1366)
      MTR272(1,1)=AAABR(1372)
      MTR272(1,2)=AAABR(1367)
      MTR272(1,3)=AAABR(1371)
      MTR272(2,1)=AAABR(1369)
      MTR272(2,2)=AAABR(1368)
      MTR272(2,3)=AAABR(1373)
      MTR272(3,1)=AAABR(1374)
      MTR272(3,2)=AAABR(1370)
      MTR272(3,3)=AAABR(1375)
      MTR273(1)=AAABR(1377)
      MTR273(2)=AAABR(1378)
      MTR273(3)=AAABR(1376)
      MTR274(1)=AAABR(1379)
      MTR274(2)=AAABR(1380)
      MTR274(3)=AAABR(1381)
      MTR275(1,1)=AAABR(1384)
      MTR275(1,2)=AAABR(1383)
      MTR275(1,3)=AAABR(1387)
      MTR275(2,1)=AAABR(1385)
      MTR275(2,2)=AAABR(1382)
      MTR275(2,3)=AAABR(1386)
      MTR276(1,1)=AAABR(1389)
      MTR276(1,2)=AAABR(1388)
      MTR276(1,3)=AAABR(1393)
      MTR276(2,1)=AAABR(1392)
      MTR276(2,2)=AAABR(1391)
      MTR276(2,3)=AAABR(1390)
      MTR277(1,1)=AAABR(1400)
      MTR277(1,2)=AAABR(1396)
      MTR277(1,3)=AAABR(1401)
      MTR277(2,1)=AAABR(1402)
      MTR277(2,2)=AAABR(1397)
      MTR277(2,3)=AAABR(1399)
      MTR277(3,1)=AAABR(1395)
      MTR277(3,2)=AAABR(1394)
      MTR277(3,3)=AAABR(1398)
      MTR278(1,1,1)=AAABR(1410)
      MTR278(1,1,2)=AAABR(1404)
      MTR278(1,1,3)=AAABR(1408)
      MTR278(1,2,1)=AAABR(1405)
      MTR278(1,2,2)=AAABR(1407)
      MTR278(1,2,3)=AAABR(1409)
      MTR278(2,1,1)=MTR278(1,2,1)
      MTR278(2,1,2)=MTR278(1,2,2)
      MTR278(2,1,3)=MTR278(1,2,3)
      MTR278(2,2,1)=AAABR(1411)
      MTR278(2,2,2)=AAABR(1406)
      MTR278(2,2,3)=AAABR(1403)
      MTR279(1,1,1)=AAABR(1415)
      MTR279(1,1,2)=AAABR(1412)
      MTR279(1,1,3)=AAABR(1413)
      MTR279(1,2,1)=AAABR(1417)
      MTR279(1,2,2)=AAABR(1416)
      MTR279(1,2,3)=AAABR(1418)
      MTR279(2,1,1)=MTR279(1,2,1)
      MTR279(2,1,2)=MTR279(1,2,2)
      MTR279(2,1,3)=MTR279(1,2,3)
      MTR279(2,2,1)=AAABR(1414)
      MTR279(2,2,2)=AAABR(1420)
      MTR279(2,2,3)=AAABR(1419)
      MTR280(1,1,1)=AAABR(1423)
      MTR280(1,1,2)=AAABR(1421)
      MTR280(1,1,3)=AAABR(1433)
      MTR280(1,2,1)=AAABR(1425)
      MTR280(1,2,2)=AAABR(1430)
      MTR280(1,2,3)=AAABR(1428)
      MTR280(1,3,1)=AAABR(1432)
      MTR280(1,3,2)=AAABR(1422)
      MTR280(1,3,3)=AAABR(1437)
      MTR280(2,1,1)=AAABR(1434)
      MTR280(2,1,2)=AAABR(1438)
      MTR280(2,1,3)=AAABR(1424)
      MTR280(2,2,1)=AAABR(1436)
      MTR280(2,2,2)=AAABR(1429)
      MTR280(2,2,3)=AAABR(1426)
      MTR280(2,3,1)=AAABR(1427)
      MTR280(2,3,2)=AAABR(1435)
      MTR280(2,3,3)=AAABR(1431)
      MTR281(1,1,1)=AAABR(1447)
      MTR281(1,1,2)=AAABR(1439)
      MTR281(1,1,3)=AAABR(1455)
      MTR281(1,2,1)=AAABR(1456)
      MTR281(1,2,2)=AAABR(1454)
      MTR281(1,2,3)=AAABR(1450)
      MTR281(1,3,1)=AAABR(1446)
      MTR281(1,3,2)=AAABR(1440)
      MTR281(1,3,3)=AAABR(1442)
      MTR281(2,1,1)=MTR281(1,2,1)
      MTR281(2,1,2)=MTR281(1,2,2)
      MTR281(2,1,3)=MTR281(1,2,3)
      MTR281(2,2,1)=AAABR(1448)
      MTR281(2,2,2)=AAABR(1443)
      MTR281(2,2,3)=AAABR(1441)
      MTR281(2,3,1)=AAABR(1449)
      MTR281(2,3,2)=AAABR(1444)
      MTR281(2,3,3)=AAABR(1452)
      MTR281(3,1,1)=MTR281(1,3,1)
      MTR281(3,1,2)=MTR281(1,3,2)
      MTR281(3,1,3)=MTR281(1,3,3)
      MTR281(3,2,1)=MTR281(2,3,1)
      MTR281(3,2,2)=MTR281(2,3,2)
      MTR281(3,2,3)=MTR281(2,3,3)
      MTR281(3,3,1)=AAABR(1451)
      MTR281(3,3,2)=AAABR(1453)
      MTR281(3,3,3)=AAABR(1445)
      MTR282(1,1,1)=AAABR(1474)
      MTR282(1,1,2)=AAABR(1458)
      MTR282(1,1,3)=AAABR(1457)
      MTR282(1,2,1)=AAABR(1465)
      MTR282(1,2,2)=AAABR(1469)
      MTR282(1,2,3)=AAABR(1466)
      MTR282(1,3,1)=AAABR(1471)
      MTR282(1,3,2)=AAABR(1468)
      MTR282(1,3,3)=AAABR(1472)
      MTR282(2,1,1)=MTR282(1,2,1)
      MTR282(2,1,2)=MTR282(1,2,2)
      MTR282(2,1,3)=MTR282(1,2,3)
      MTR282(2,2,1)=AAABR(1467)
      MTR282(2,2,2)=AAABR(1462)
      MTR282(2,2,3)=AAABR(1464)
      MTR282(2,3,1)=AAABR(1459)
      MTR282(2,3,2)=AAABR(1461)
      MTR282(2,3,3)=AAABR(1463)
      MTR282(3,1,1)=MTR282(1,3,1)
      MTR282(3,1,2)=MTR282(1,3,2)
      MTR282(3,1,3)=MTR282(1,3,3)
      MTR282(3,2,1)=MTR282(2,3,1)
      MTR282(3,2,2)=MTR282(2,3,2)
      MTR282(3,2,3)=MTR282(2,3,3)
      MTR282(3,3,1)=AAABR(1460)
      MTR282(3,3,2)=AAABR(1473)
      MTR282(3,3,3)=AAABR(1470)
      MTR283(1,1)=AAABR(1482)
      MTR283(1,2)=AAABR(1165)
      MTR283(1,3)=AAABR(1165)
      MTR283(2,1)=AAABR(1165)
      MTR283(2,2)=AAABR(1475)
      MTR283(2,3)=AAABR(1165)
      MTR283(3,1)=AAABR(1165)
      MTR283(3,2)=AAABR(1165)
      MTR283(3,3)=AAABR(1478)
      MTR284(1,1)=AAABR(1483)
      MTR284(1,2)=AAABR(1480)
      MTR284(1,3)=AAABR(1477)
      MTR284(2,1)=AAABR(1480)
      MTR284(2,2)=AAABR(1476)
      MTR284(2,3)=AAABR(1481)
      MTR284(3,1)=AAABR(1477)
      MTR284(3,2)=AAABR(1481)
      MTR284(3,3)=AAABR(1479)
      MTR285(1,1)=AAABR(1276)
      MTR285(1,2)=AAABR(1262)
      MTR285(1,3)=AAABR(1261)
      MTR285(2,1)=AAABR(1266)
      MTR285(2,2)=AAABR(1270)
      MTR285(2,3)=AAABR(1273)
      MTR285(3,1)=AAABR(1271)
      MTR285(3,2)=AAABR(1274)
      MTR285(3,3)=AAABR(1264)
      MTR286(1,1)=AAABR(1282)
      MTR286(1,2)=AAABR(1312)
      MTR286(1,3)=AAABR(1315)
      MTR286(2,1)=AAABR(1314)
      MTR286(2,2)=AAABR(1278)
      MTR286(2,3)=AAABR(1310)
      MTR286(3,1)=AAABR(1313)
      MTR286(3,2)=AAABR(1311)
      MTR286(3,3)=AAABR(1288)
      MTR287(1,1)=AAABR(1485)
      MTR287(1,2)=AAABR(1484)
      MTR287(1,3)=AAABR(1491)
      MTR287(2,1)=AAABR(1492)
      MTR287(2,2)=AAABR(1497)
      MTR287(2,3)=AAABR(1496)
      MTR287(3,1)=AAABR(1489)
      MTR287(3,2)=AAABR(1490)
      MTR287(3,3)=AAABR(1487)
      MTR288(1,1)=AAABR(1486)
      MTR288(1,2)=AAABR(1493)
      MTR288(1,3)=AAABR(1494)
      MTR288(2,1)=AAABR(1493)
      MTR288(2,2)=AAABR(1498)
      MTR288(2,3)=AAABR(1495)
      MTR288(3,1)=AAABR(1494)
      MTR288(3,2)=AAABR(1495)
      MTR288(3,3)=AAABR(1488)
      MTR289(1,1)=AAABR(1503)
      MTR289(1,2)=AAABR(1516)
      MTR289(1,3)=AAABR(1509)
      MTR289(2,1)=AAABR(1506)
      MTR289(2,2)=AAABR(1507)
      MTR289(2,3)=AAABR(1505)
      MTR289(3,1)=AAABR(1511)
      MTR289(3,2)=AAABR(1512)
      MTR289(3,3)=AAABR(1499)
      MTR290(1,1)=AAABR(1504)
      MTR290(1,2)=AAABR(1510)
      MTR290(1,3)=AAABR(1514)
      MTR290(2,1)=AAABR(1513)
      MTR290(2,2)=AAABR(1508)
      MTR290(2,3)=AAABR(1515)
      MTR290(3,1)=AAABR(1502)
      MTR290(3,2)=AAABR(1501)
      MTR290(3,3)=AAABR(1500)
      MTR291(1,1)=AAABR(1282)
      MTR291(1,2)=AAABR(1283)
      MTR291(1,3)=AAABR(1292)
      MTR291(2,1)=AAABR(1290)
      MTR291(2,2)=AAABR(1278)
      MTR291(2,3)=AAABR(1284)
      MTR291(3,1)=AAABR(1289)
      MTR291(3,2)=AAABR(1286)
      MTR291(3,3)=AAABR(1288)
      MTR292(1,1)=AAABR(1333)
      MTR292(1,2)=AAABR(1320)
      MTR292(1,3)=AAABR(1329)
      MTR292(2,1)=AAABR(1321)
      MTR292(2,2)=AAABR(1317)
      MTR292(2,3)=AAABR(1323)
      MTR292(3,1)=AAABR(1322)
      MTR292(3,2)=AAABR(1326)
      MTR292(3,3)=AAABR(1319)
      MTR293(1,1)=AAABR(1504)
      MTR293(1,2)=AAABR(1513)
      MTR293(1,3)=AAABR(1502)
      MTR293(2,1)=AAABR(1510)
      MTR293(2,2)=AAABR(1508)
      MTR293(2,3)=AAABR(1501)
      MTR293(3,1)=AAABR(1514)
      MTR293(3,2)=AAABR(1515)
      MTR293(3,3)=AAABR(1500)
      MTR294(1,1)=AAABR(1528)
      MTR294(1,2)=AAABR(1527)
      MTR294(1,3)=AAABR(1531)
      MTR294(2,1)=AAABR(1530)
      MTR294(2,2)=AAABR(1522)
      MTR294(2,3)=AAABR(1524)
      MTR294(3,1)=AAABR(1519)
      MTR294(3,2)=AAABR(1526)
      MTR294(3,3)=AAABR(1517)
      MTR295(1,1)=AAABR(1529)
      MTR295(1,2)=AAABR(1521)
      MTR295(1,3)=AAABR(1525)
      MTR295(2,1)=AAABR(1521)
      MTR295(2,2)=AAABR(1523)
      MTR295(2,3)=AAABR(1520)
      MTR295(3,1)=AAABR(1525)
      MTR295(3,2)=AAABR(1520)
      MTR295(3,3)=AAABR(1518)
      MTR296(1)=AAABR(1532)
      MTR296(2)=AAABR(1534)
      MTR296(3)=AAABR(1533)
      MTR297(1)=AAABR(1535)
      MTR297(2)=AAABR(1536)
      MTR297(3)=AAABR(1537)
      MTR298(1)=AAABR(1538)
      MTR298(2)=AAABR(1539)
      MTR298(3)=AAABR(1540)
      MTR299(1,1)=AAABR(1545)
      MTR299(1,2)=AAABR(1544)
      MTR299(1,3)=AAABR(1541)
      MTR299(2,1)=AAABR(1543)
      MTR299(2,2)=AAABR(1546)
      MTR299(2,3)=AAABR(1542)
      MTR300(1,1)=AAABR(1550)
      MTR300(1,2)=AAABR(1547)
      MTR300(1,3)=AAABR(1548)
      MTR300(2,1)=AAABR(1549)
      MTR300(2,2)=AAABR(1551)
      MTR300(2,3)=AAABR(1552)
      MTR301(1,1)=AAABR(1557)
      MTR301(1,2)=AAABR(1556)
      MTR301(1,3)=AAABR(1555)
      MTR301(2,1)=AAABR(1559)
      MTR301(2,2)=AAABR(1561)
      MTR301(2,3)=AAABR(1558)
      MTR301(3,1)=AAABR(1560)
      MTR301(3,2)=AAABR(1553)
      MTR301(3,3)=AAABR(1554)
      MTR302(1,1)=AAABR(1570)
      MTR302(1,2)=AAABR(1569)
      MTR302(1,3)=AAABR(1563)
      MTR302(2,1)=AAABR(1565)
      MTR302(2,2)=AAABR(1566)
      MTR302(2,3)=AAABR(1568)
      MTR302(3,1)=AAABR(1567)
      MTR302(3,2)=AAABR(1562)
      MTR302(3,3)=AAABR(1564)
      MTR303(1)=AAABR(1573)
      MTR303(2)=AAABR(1571)
      MTR303(3)=AAABR(1572)
      MTR304(1,1)=AAABR(1575)
      MTR304(1,2)=AAABR(1577)
      MTR304(1,3)=AAABR(1578)
      MTR304(2,1)=AAABR(1579)
      MTR304(2,2)=AAABR(1576)
      MTR304(2,3)=AAABR(1574)
      MTR305(1,1,1)=AAABR(1583)
      MTR305(1,1,2)=AAABR(1588)
      MTR305(1,1,3)=AAABR(1586)
      MTR305(1,2,1)=AAABR(1582)
      MTR305(1,2,2)=AAABR(1580)
      MTR305(1,2,3)=AAABR(1585)
      MTR305(2,1,1)=MTR305(1,2,1)
      MTR305(2,1,2)=MTR305(1,2,2)
      MTR305(2,1,3)=MTR305(1,2,3)
      MTR305(2,2,1)=AAABR(1581)
      MTR305(2,2,2)=AAABR(1587)
      MTR305(2,2,3)=AAABR(1584)
      MTR306(1,1,1)=AAABR(1594)
      MTR306(1,1,2)=AAABR(1597)
      MTR306(1,1,3)=AAABR(1595)
      MTR306(1,2,1)=AAABR(1591)
      MTR306(1,2,2)=AAABR(1606)
      MTR306(1,2,3)=AAABR(1592)
      MTR306(1,3,1)=AAABR(1589)
      MTR306(1,3,2)=AAABR(1599)
      MTR306(1,3,3)=AAABR(1590)
      MTR306(2,1,1)=MTR306(1,2,1)
      MTR306(2,1,2)=MTR306(1,2,2)
      MTR306(2,1,3)=MTR306(1,2,3)
      MTR306(2,2,1)=AAABR(1598)
      MTR306(2,2,2)=AAABR(1604)
      MTR306(2,2,3)=AAABR(1602)
      MTR306(2,3,1)=AAABR(1593)
      MTR306(2,3,2)=AAABR(1601)
      MTR306(2,3,3)=AAABR(1596)
      MTR306(3,1,1)=MTR306(1,3,1)
      MTR306(3,1,2)=MTR306(1,3,2)
      MTR306(3,1,3)=MTR306(1,3,3)
      MTR306(3,2,1)=MTR306(2,3,1)
      MTR306(3,2,2)=MTR306(2,3,2)
      MTR306(3,2,3)=MTR306(2,3,3)
      MTR306(3,3,1)=AAABR(1600)
      MTR306(3,3,2)=AAABR(1603)
      MTR306(3,3,3)=AAABR(1605)
      MTR307(1,1)=AAABR(1613)
      MTR307(1,2)=AAABR(1607)
      MTR307(1,3)=AAABR(1610)
      MTR307(2,1)=AAABR(1607)
      MTR307(2,2)=AAABR(1609)
      MTR307(2,3)=AAABR(1611)
      MTR307(3,1)=AAABR(1610)
      MTR307(3,2)=AAABR(1611)
      MTR307(3,3)=AAABR(1615)
      MTR308(1,1)=AAABR(1612)
      MTR308(1,2)=AAABR(1227)
      MTR308(1,3)=AAABR(1227)
      MTR308(2,1)=AAABR(1227)
      MTR308(2,2)=AAABR(1608)
      MTR308(2,3)=AAABR(1227)
      MTR308(3,1)=AAABR(1227)
      MTR308(3,2)=AAABR(1227)
      MTR308(3,3)=AAABR(1614)
      MTR309(1,1)=AAABR(1627)
      MTR309(1,2)=AAABR(1624)
      MTR309(1,3)=AAABR(1616)
      MTR309(2,1)=AAABR(1623)
      MTR309(2,2)=AAABR(1618)
      MTR309(2,3)=AAABR(1622)
      MTR309(3,1)=AAABR(1621)
      MTR309(3,2)=AAABR(1625)
      MTR309(3,3)=AAABR(1620)
      MTR310(1,1)=AAABR(1626)
      MTR310(2,2)=AAABR(1617)
      MTR310(3,3)=AAABR(1619)
      MTR311(1,1)=AAABR(1628)
      MTR311(2,2)=AAABR(1631)
      MTR311(3,3)=AAABR(1633)
      MTR312(1,1)=AAABR(1629)
      MTR312(1,2)=AAABR(1630)
      MTR312(1,3)=AAABR(1635)
      MTR312(2,1)=AAABR(1630)
      MTR312(2,2)=AAABR(1632)
      MTR312(2,3)=AAABR(1636)
      MTR312(3,1)=AAABR(1635)
      MTR312(3,2)=AAABR(1636)
      MTR312(3,3)=AAABR(1634)
      MTR313(1,1)=AAABR(1627)
      MTR313(1,2)=AAABR(1623)
      MTR313(1,3)=AAABR(1621)
      MTR313(2,1)=AAABR(1624)
      MTR313(2,2)=AAABR(1618)
      MTR313(2,3)=AAABR(1625)
      MTR313(3,1)=AAABR(1616)
      MTR313(3,2)=AAABR(1622)
      MTR313(3,3)=AAABR(1620)
      MTR314(1,1)=AAABR(1654)
      MTR314(1,2)=AAABR(1641)
      MTR314(1,3)=AAABR(1645)
      MTR314(2,1)=AAABR(1641)
      MTR314(2,2)=AAABR(1650)
      MTR314(2,3)=AAABR(1637)
      MTR314(3,1)=AAABR(1645)
      MTR314(3,2)=AAABR(1637)
      MTR314(3,3)=AAABR(1640)
      MTR315(1,1)=AAABR(1652)
      MTR315(1,2)=AAABR(1227)
      MTR315(1,3)=AAABR(1227)
      MTR315(2,1)=AAABR(1227)
      MTR315(2,2)=AAABR(1648)
      MTR315(2,3)=AAABR(1227)
      MTR315(3,1)=AAABR(1227)
      MTR315(3,2)=AAABR(1227)
      MTR315(3,3)=AAABR(1638)
      MTR316(1,1)=AAABR(1653)
      MTR316(1,2)=AAABR(1642)
      MTR316(1,3)=AAABR(1647)
      MTR316(2,1)=AAABR(1644)
      MTR316(2,2)=AAABR(1649)
      MTR316(2,3)=AAABR(1646)
      MTR316(3,1)=AAABR(1651)
      MTR316(3,2)=AAABR(1643)
      MTR316(3,3)=AAABR(1639)
      MTR317(1,1)=AAABR(1664)
      MTR317(1,2)=AAABR(1655)
      MTR317(1,3)=AAABR(1658)
      MTR317(2,1)=AAABR(1665)
      MTR317(2,2)=AAABR(1657)
      MTR317(2,3)=AAABR(1662)
      MTR317(3,1)=AAABR(1659)
      MTR317(3,2)=AAABR(1666)
      MTR317(3,3)=AAABR(1661)
      MTR318(1,1)=AAABR(1664)
      MTR318(1,2)=AAABR(1665)
      MTR318(1,3)=AAABR(1659)
      MTR318(2,1)=AAABR(1655)
      MTR318(2,2)=AAABR(1657)
      MTR318(2,3)=AAABR(1666)
      MTR318(3,1)=AAABR(1658)
      MTR318(3,2)=AAABR(1662)
      MTR318(3,3)=AAABR(1661)
      MTR319(1,1)=AAABR(1663)
      MTR319(2,2)=AAABR(1656)
      MTR319(3,3)=AAABR(1660)
      MTR320(1)=AAABR(1667)
      MTR320(2)=AAABR(1668)
      MTR320(3)=AAABR(1669)
      MTR321(1)=AAABR(1672)
      MTR321(2)=AAABR(1670)
      MTR321(3)=AAABR(1671)
      MTR322(1)=AAABR(1675)
      MTR322(2)=AAABR(1673)
      MTR322(3)=AAABR(1674)
      MTR323(1)=AAABR(1677)
      MTR323(2)=AAABR(1678)
      MTR323(3)=AAABR(1676)
      MTR324(1,1)=AAABR(1680)
      MTR324(1,2)=AAABR(1679)
      MTR324(1,3)=AAABR(1682)
      MTR324(2,1)=AAABR(1684)
      MTR324(2,2)=AAABR(1683)
      MTR324(2,3)=AAABR(1681)
      MTR325(1,1)=AAABR(1689)
      MTR325(1,2)=AAABR(1687)
      MTR325(1,3)=AAABR(1688)
      MTR325(2,1)=AAABR(1686)
      MTR325(2,2)=AAABR(1685)
      MTR325(2,3)=AAABR(1690)
      MTR326(1,1)=AAABR(1691)
      MTR326(1,2)=AAABR(1694)
      MTR326(1,3)=AAABR(1697)
      MTR326(2,1)=AAABR(1695)
      MTR326(2,2)=AAABR(1692)
      MTR326(2,3)=AAABR(1693)
      MTR326(3,1)=AAABR(1696)
      MTR326(3,2)=AAABR(1699)
      MTR326(3,3)=AAABR(1698)
      MTR327(1,1,1)=AAABR(1707)
      MTR327(1,1,2)=AAABR(1705)
      MTR327(1,1,3)=AAABR(1706)
      MTR327(1,2,1)=AAABR(1700)
      MTR327(1,2,2)=AAABR(1703)
      MTR327(1,2,3)=AAABR(1701)
      MTR327(2,1,1)=MTR327(1,2,1)
      MTR327(2,1,2)=MTR327(1,2,2)
      MTR327(2,1,3)=MTR327(1,2,3)
      MTR327(2,2,1)=AAABR(1702)
      MTR327(2,2,2)=AAABR(1708)
      MTR327(2,2,3)=AAABR(1704)
      MTR328(1,1,1)=AAABR(1710)
      MTR328(1,1,2)=AAABR(1716)
      MTR328(1,1,3)=AAABR(1715)
      MTR328(1,2,1)=AAABR(1713)
      MTR328(1,2,2)=AAABR(1709)
      MTR328(1,2,3)=AAABR(1714)
      MTR328(2,1,1)=MTR328(1,2,1)
      MTR328(2,1,2)=MTR328(1,2,2)
      MTR328(2,1,3)=MTR328(1,2,3)
      MTR328(2,2,1)=AAABR(1711)
      MTR328(2,2,2)=AAABR(1712)
      MTR328(2,2,3)=AAABR(1717)
      MTR329(1,1,1)=AAABR(1721)
      MTR329(1,1,2)=AAABR(1724)
      MTR329(1,1,3)=AAABR(1735)
      MTR329(1,2,1)=AAABR(1723)
      MTR329(1,2,2)=AAABR(1719)
      MTR329(1,2,3)=AAABR(1733)
      MTR329(1,3,1)=AAABR(1726)
      MTR329(1,3,2)=AAABR(1725)
      MTR329(1,3,3)=AAABR(1734)
      MTR329(2,1,1)=AAABR(1722)
      MTR329(2,1,2)=AAABR(1731)
      MTR329(2,1,3)=AAABR(1729)
      MTR329(2,2,1)=AAABR(1720)
      MTR329(2,2,2)=AAABR(1727)
      MTR329(2,2,3)=AAABR(1728)
      MTR329(2,3,1)=AAABR(1730)
      MTR329(2,3,2)=AAABR(1718)
      MTR329(2,3,3)=AAABR(1732)
      MTR330(1,1,1)=AAABR(1753)
      MTR330(1,1,2)=AAABR(1743)
      MTR330(1,1,3)=AAABR(1747)
      MTR330(1,2,1)=AAABR(1736)
      MTR330(1,2,2)=AAABR(1742)
      MTR330(1,2,3)=AAABR(1737)
      MTR330(1,3,1)=AAABR(1738)
      MTR330(1,3,2)=AAABR(1749)
      MTR330(1,3,3)=AAABR(1748)
      MTR330(2,1,1)=MTR330(1,2,1)
      MTR330(2,1,2)=MTR330(1,2,2)
      MTR330(2,1,3)=MTR330(1,2,3)
      MTR330(2,2,1)=AAABR(1739)
      MTR330(2,2,2)=AAABR(1744)
      MTR330(2,2,3)=AAABR(1746)
      MTR330(2,3,1)=AAABR(1750)
      MTR330(2,3,2)=AAABR(1745)
      MTR330(2,3,3)=AAABR(1752)
      MTR330(3,1,1)=MTR330(1,3,1)
      MTR330(3,1,2)=MTR330(1,3,2)
      MTR330(3,1,3)=MTR330(1,3,3)
      MTR330(3,2,1)=MTR330(2,3,1)
      MTR330(3,2,2)=MTR330(2,3,2)
      MTR330(3,2,3)=MTR330(2,3,3)
      MTR330(3,3,1)=AAABR(1751)
      MTR330(3,3,2)=AAABR(1740)
      MTR330(3,3,3)=AAABR(1741)
      MTR331(1,1,1)=AAABR(1764)
      MTR331(1,1,2)=AAABR(1758)
      MTR331(1,1,3)=AAABR(1768)
      MTR331(1,2,1)=AAABR(1760)
      MTR331(1,2,2)=AAABR(1767)
      MTR331(1,2,3)=AAABR(1754)
      MTR331(1,3,1)=AAABR(1755)
      MTR331(1,3,2)=AAABR(1756)
      MTR331(1,3,3)=AAABR(1762)
      MTR331(2,1,1)=MTR331(1,2,1)
      MTR331(2,1,2)=MTR331(1,2,2)
      MTR331(2,1,3)=MTR331(1,2,3)
      MTR331(2,2,1)=AAABR(1765)
      MTR331(2,2,2)=AAABR(1757)
      MTR331(2,2,3)=AAABR(1766)
      MTR331(2,3,1)=AAABR(1769)
      MTR331(2,3,2)=AAABR(1771)
      MTR331(2,3,3)=AAABR(1759)
      MTR331(3,1,1)=MTR331(1,3,1)
      MTR331(3,1,2)=MTR331(1,3,2)
      MTR331(3,1,3)=MTR331(1,3,3)
      MTR331(3,2,1)=MTR331(2,3,1)
      MTR331(3,2,2)=MTR331(2,3,2)
      MTR331(3,2,3)=MTR331(2,3,3)
      MTR331(3,3,1)=AAABR(1770)
      MTR331(3,3,2)=AAABR(1763)
      MTR331(3,3,3)=AAABR(1761)
      MTR332(1,1)=AAABR(1779)
      MTR332(1,2)=AAABR(1775)
      MTR332(1,3)=AAABR(1772)
      MTR332(2,1)=AAABR(1775)
      MTR332(2,2)=AAABR(1777)
      MTR332(2,3)=AAABR(1780)
      MTR332(3,1)=AAABR(1772)
      MTR332(3,2)=AAABR(1780)
      MTR332(3,3)=AAABR(1774)
      MTR333(1,1)=AAABR(1778)
      MTR333(1,2)=AAABR(1227)
      MTR333(1,3)=AAABR(1227)
      MTR333(2,1)=AAABR(1227)
      MTR333(2,2)=AAABR(1776)
      MTR333(2,3)=AAABR(1227)
      MTR333(3,1)=AAABR(1227)
      MTR333(3,2)=AAABR(1227)
      MTR333(3,3)=AAABR(1773)
      MTR334(1)=AAABR(1782)
      MTR334(2)=AAABR(1781)
      MTR334(3)=AAABR(1783)
      MTR335(1)=AAABR(1786)
      MTR335(2)=AAABR(1784)
      MTR335(3)=AAABR(1785)
      MTR336(1,1)=AAABR(1788)
      MTR336(1,2)=AAABR(1792)
      MTR336(1,3)=AAABR(1791)
      MTR336(2,1)=AAABR(1790)
      MTR336(2,2)=AAABR(1789)
      MTR336(2,3)=AAABR(1787)
      MTR337(1,1,1)=AAABR(1799)
      MTR337(1,1,2)=AAABR(1793)
      MTR337(1,1,3)=AAABR(1794)
      MTR337(1,2,1)=AAABR(1797)
      MTR337(1,2,2)=AAABR(1796)
      MTR337(1,2,3)=AAABR(1798)
      MTR337(2,1,1)=MTR337(1,2,1)
      MTR337(2,1,2)=MTR337(1,2,2)
      MTR337(2,1,3)=MTR337(1,2,3)
      MTR337(2,2,1)=AAABR(1800)
      MTR337(2,2,2)=AAABR(1795)
      MTR337(2,2,3)=AAABR(1801)
      MTR338(1,1,1)=AAABR(1815)
      MTR338(1,1,2)=AAABR(1814)
      MTR338(1,1,3)=AAABR(1812)
      MTR338(1,2,1)=AAABR(1817)
      MTR338(1,2,2)=AAABR(1806)
      MTR338(1,2,3)=AAABR(1819)
      MTR338(1,3,1)=AAABR(1805)
      MTR338(1,3,2)=AAABR(1818)
      MTR338(1,3,3)=AAABR(1810)
      MTR338(2,1,1)=MTR338(1,2,1)
      MTR338(2,1,2)=MTR338(1,2,2)
      MTR338(2,1,3)=MTR338(1,2,3)
      MTR338(2,2,1)=AAABR(1808)
      MTR338(2,2,2)=AAABR(1804)
      MTR338(2,2,3)=AAABR(1813)
      MTR338(2,3,1)=AAABR(1807)
      MTR338(2,3,2)=AAABR(1811)
      MTR338(2,3,3)=AAABR(1803)
      MTR338(3,1,1)=MTR338(1,3,1)
      MTR338(3,1,2)=MTR338(1,3,2)
      MTR338(3,1,3)=MTR338(1,3,3)
      MTR338(3,2,1)=MTR338(2,3,1)
      MTR338(3,2,2)=MTR338(2,3,2)
      MTR338(3,2,3)=MTR338(2,3,3)
      MTR338(3,3,1)=AAABR(1816)
      MTR338(3,3,2)=AAABR(1802)
      MTR338(3,3,3)=AAABR(1809)
      MTR339(1)=AAABR(1822)
      MTR339(2)=AAABR(1821)
      MTR340(1,1)=AAABR(1824)
      MTR340(1,2)=AAABR(1825)
      MTR340(2,1)=MTR340(1,2)
      MTR340(2,2)=AAABR(1823)
      MTR341(1,1)=AAABR(1828)
      MTR341(1,2)=AAABR(1827)
      MTR341(1,3)=AAABR(1826)
      MTR341(2,1)=MTR341(1,2)
      MTR341(2,2)=AAABR(1829)
      MTR341(2,3)=AAABR(1830)
      MTR341(3,1)=MTR341(1,3)
      MTR341(3,2)=MTR341(2,3)
      MTR341(3,3)=AAABR(1831)
      MTR342(1,1)=AAABR(1833)
      MTR342(1,2)=AAABR(1834)
      MTR342(2,1)=MTR342(1,2)
      MTR342(2,2)=AAABR(1832)
      MTR343(1,1)=AAABR(1836)
      MTR343(1,2)=AAABR(1838)
      MTR343(1,3)=AAABR(1840)
      MTR343(2,1)=MTR343(1,2)
      MTR343(2,2)=AAABR(1837)
      MTR343(2,3)=AAABR(1839)
      MTR343(3,1)=MTR343(1,3)
      MTR343(3,2)=MTR343(2,3)
      MTR343(3,3)=AAABR(1835)
      MTR344(1,1,1)=3D0*AAABR(1843)
      MTR344(1,1,2)=AAABR(1844)
      MTR344(1,2,1)=MTR344(1,1,2)
      MTR344(1,2,2)=AAABR(1842)
      MTR344(2,1,1)=MTR344(1,1,2)
      MTR344(2,1,2)=MTR344(1,2,2)
      MTR344(2,2,1)=MTR344(1,2,2)
      MTR344(2,2,2)=3D0*AAABR(1841)
      MTR345(1,1,1)=AAABR(1856)
      MTR345(1,1,2)=AAABR(1849)
      MTR345(1,1,3)=AAABR(1848)
      MTR345(1,2,1)=MTR345(1,1,2)
      MTR345(1,2,2)=AAABR(1851)
      MTR345(1,2,3)=AAABR(1847)
      MTR345(1,3,1)=MTR345(1,1,3)
      MTR345(1,3,2)=MTR345(1,2,3)
      MTR345(1,3,3)=AAABR(1845)
      MTR345(2,1,1)=AAABR(1853)
      MTR345(2,1,2)=AAABR(1850)
      MTR345(2,1,3)=AAABR(1846)
      MTR345(2,2,1)=MTR345(2,1,2)
      MTR345(2,2,2)=AAABR(1854)
      MTR345(2,2,3)=AAABR(1855)
      MTR345(2,3,1)=MTR345(2,1,3)
      MTR345(2,3,2)=MTR345(2,2,3)
      MTR345(2,3,3)=AAABR(1852)
      MTR346(1,1,1,1)=3D0*AAABR(1857)
      MTR346(1,1,1,2)=3D0*AAABR(1859)
      MTR346(1,1,2,1)=MTR346(1,1,1,2)
      MTR346(1,1,2,2)=AAABR(1861)
      MTR346(1,2,1,1)=MTR346(1,1,1,2)
      MTR346(1,2,1,2)=MTR346(1,1,2,2)
      MTR346(1,2,2,1)=MTR346(1,1,2,2)
      MTR346(1,2,2,2)=3D0*AAABR(1860)
      MTR346(2,1,1,1)=MTR346(1,1,1,2)
      MTR346(2,1,1,2)=MTR346(1,1,2,2)
      MTR346(2,1,2,1)=MTR346(1,1,2,2)
      MTR346(2,1,2,2)=MTR346(1,2,2,2)
      MTR346(2,2,1,1)=MTR346(1,1,2,2)
      MTR346(2,2,1,2)=MTR346(1,2,2,2)
      MTR346(2,2,2,1)=MTR346(1,2,2,2)
      MTR346(2,2,2,2)=3D0*AAABR(1858)
      MTR347(1,1,1,1)=AAABR(1871)
      MTR347(1,1,1,2)=AAABR(1867)
      MTR347(1,1,1,3)=AAABR(1879)
      MTR347(1,1,2,1)=AAABR(1867)
      MTR347(1,1,2,2)=AAABR(1878)
      MTR347(1,1,2,3)=AAABR(1877)
      MTR347(1,1,3,1)=AAABR(1879)
      MTR347(1,1,3,2)=AAABR(1877)
      MTR347(1,1,3,3)=AAABR(1866)
      MTR347(1,2,1,1)=AAABR(1868)
      MTR347(1,2,1,2)=AAABR(1865)
      MTR347(1,2,1,3)=AAABR(1863)
      MTR347(1,2,2,1)=AAABR(1865)
      MTR347(1,2,2,2)=AAABR(1869)
      MTR347(1,2,2,3)=AAABR(1873)
      MTR347(1,2,3,1)=AAABR(1863)
      MTR347(1,2,3,2)=AAABR(1873)
      MTR347(1,2,3,3)=AAABR(1862)
      MTR347(2,1,1,1)=AAABR(1868)
      MTR347(2,1,1,2)=AAABR(1865)
      MTR347(2,1,1,3)=AAABR(1863)
      MTR347(2,1,2,1)=AAABR(1865)
      MTR347(2,1,2,2)=AAABR(1869)
      MTR347(2,1,2,3)=AAABR(1873)
      MTR347(2,1,3,1)=AAABR(1863)
      MTR347(2,1,3,2)=AAABR(1873)
      MTR347(2,1,3,3)=AAABR(1862)
      MTR347(2,2,1,1)=AAABR(1874)
      MTR347(2,2,1,2)=AAABR(1875)
      MTR347(2,2,1,3)=AAABR(1870)
      MTR347(2,2,2,1)=AAABR(1875)
      MTR347(2,2,2,2)=AAABR(1876)
      MTR347(2,2,2,3)=AAABR(1864)
      MTR347(2,2,3,1)=AAABR(1870)
      MTR347(2,2,3,2)=AAABR(1864)
      MTR347(2,2,3,3)=AAABR(1872)
      MTR348(1,1,1,1)=3D0*AAABR(1883)
      MTR348(1,1,1,2)=3D0*AAABR(1885)
      MTR348(1,1,1,3)=3D0*AAABR(1881)
      MTR348(1,1,2,1)=MTR348(1,1,1,2)
      MTR348(1,1,2,2)=AAABR(1890)
      MTR348(1,1,2,3)=AAABR(1893)
      MTR348(1,1,3,1)=MTR348(1,1,1,3)
      MTR348(1,1,3,2)=MTR348(1,1,2,3)
      MTR348(1,1,3,3)=AAABR(1892)
      MTR348(1,2,1,1)=MTR348(1,1,1,2)
      MTR348(1,2,1,2)=MTR348(1,1,2,2)
      MTR348(1,2,1,3)=MTR348(1,1,2,3)
      MTR348(1,2,2,1)=MTR348(1,1,2,2)
      MTR348(1,2,2,2)=3D0*AAABR(1889)
      MTR348(1,2,2,3)=AAABR(1887)
      MTR348(1,2,3,1)=MTR348(1,1,2,3)
      MTR348(1,2,3,2)=MTR348(1,2,2,3)
      MTR348(1,2,3,3)=AAABR(1894)
      MTR348(1,3,1,1)=MTR348(1,1,1,3)
      MTR348(1,3,1,2)=MTR348(1,1,2,3)
      MTR348(1,3,1,3)=MTR348(1,1,3,3)
      MTR348(1,3,2,1)=MTR348(1,1,2,3)
      MTR348(1,3,2,2)=MTR348(1,2,2,3)
      MTR348(1,3,2,3)=MTR348(1,2,3,3)
      MTR348(1,3,3,1)=MTR348(1,1,3,3)
      MTR348(1,3,3,2)=MTR348(1,2,3,3)
      MTR348(1,3,3,3)=3D0*AAABR(1880)
      MTR348(2,1,1,1)=MTR348(1,1,1,2)
      MTR348(2,1,1,2)=MTR348(1,1,2,2)
      MTR348(2,1,1,3)=MTR348(1,1,2,3)
      MTR348(2,1,2,1)=MTR348(1,1,2,2)
      MTR348(2,1,2,2)=MTR348(1,2,2,2)
      MTR348(2,1,2,3)=MTR348(1,2,2,3)
      MTR348(2,1,3,1)=MTR348(1,1,2,3)
      MTR348(2,1,3,2)=MTR348(1,2,2,3)
      MTR348(2,1,3,3)=MTR348(1,2,3,3)
      MTR348(2,2,1,1)=MTR348(1,1,2,2)
      MTR348(2,2,1,2)=MTR348(1,2,2,2)
      MTR348(2,2,1,3)=MTR348(1,2,2,3)
      MTR348(2,2,2,1)=MTR348(1,2,2,2)
      MTR348(2,2,2,2)=3D0*AAABR(1882)
      MTR348(2,2,2,3)=3D0*AAABR(1891)
      MTR348(2,2,3,1)=MTR348(1,2,2,3)
      MTR348(2,2,3,2)=MTR348(2,2,2,3)
      MTR348(2,2,3,3)=AAABR(1888)
      MTR348(2,3,1,1)=MTR348(1,1,2,3)
      MTR348(2,3,1,2)=MTR348(1,2,2,3)
      MTR348(2,3,1,3)=MTR348(1,2,3,3)
      MTR348(2,3,2,1)=MTR348(1,2,2,3)
      MTR348(2,3,2,2)=MTR348(2,2,2,3)
      MTR348(2,3,2,3)=MTR348(2,2,3,3)
      MTR348(2,3,3,1)=MTR348(1,2,3,3)
      MTR348(2,3,3,2)=MTR348(2,2,3,3)
      MTR348(2,3,3,3)=3D0*AAABR(1886)
      MTR348(3,1,1,1)=MTR348(1,1,1,3)
      MTR348(3,1,1,2)=MTR348(1,1,2,3)
      MTR348(3,1,1,3)=MTR348(1,1,3,3)
      MTR348(3,1,2,1)=MTR348(1,1,2,3)
      MTR348(3,1,2,2)=MTR348(1,2,2,3)
      MTR348(3,1,2,3)=MTR348(1,2,3,3)
      MTR348(3,1,3,1)=MTR348(1,1,3,3)
      MTR348(3,1,3,2)=MTR348(1,2,3,3)
      MTR348(3,1,3,3)=MTR348(1,3,3,3)
      MTR348(3,2,1,1)=MTR348(1,1,2,3)
      MTR348(3,2,1,2)=MTR348(1,2,2,3)
      MTR348(3,2,1,3)=MTR348(1,2,3,3)
      MTR348(3,2,2,1)=MTR348(1,2,2,3)
      MTR348(3,2,2,2)=MTR348(2,2,2,3)
      MTR348(3,2,2,3)=MTR348(2,2,3,3)
      MTR348(3,2,3,1)=MTR348(1,2,3,3)
      MTR348(3,2,3,2)=MTR348(2,2,3,3)
      MTR348(3,2,3,3)=MTR348(2,3,3,3)
      MTR348(3,3,1,1)=MTR348(1,1,3,3)
      MTR348(3,3,1,2)=MTR348(1,2,3,3)
      MTR348(3,3,1,3)=MTR348(1,3,3,3)
      MTR348(3,3,2,1)=MTR348(1,2,3,3)
      MTR348(3,3,2,2)=MTR348(2,2,3,3)
      MTR348(3,3,2,3)=MTR348(2,3,3,3)
      MTR348(3,3,3,1)=MTR348(1,3,3,3)
      MTR348(3,3,3,2)=MTR348(2,3,3,3)
      MTR348(3,3,3,3)=3D0*AAABR(1884)
      MTR349(1)=AAABR(1901)
      MTR349(2)=AAABR(1899)
      MTR350(1)=AAABR(1900)
      MTR350(2)=AAABR(1898)
      MTR351(1)=AAABR(1905)
      MTR351(2)=AAABR(1903)
      MTR351(3)=AAABR(1907)
      MTR352(1)=AAABR(1904)
      MTR352(2)=AAABR(1902)
      MTR352(3)=AAABR(1906)
      MTR353(1)=AAABR(1895)
      MTR353(2)=AAABR(1895)
      MTR353(3)=AAABR(1908)
      MTR354(1)=AAABR(907)
      MTR354(2)=AAABR(907)
      MTR354(3)=AAABR(1910)
      MTR355(1)=AAABR(1897)
      MTR355(2)=AAABR(1897)
      MTR355(3)=AAABR(1909)
      MTR356(1)=AAABR(1917)
      MTR356(2)=AAABR(1917)
      MTR356(3)=AAABR(1915)
      MTR357(1)=AAABR(1916)
      MTR357(2)=AAABR(1916)
      MTR357(3)=AAABR(1914)
      MTR358(1)=2D0*AAABR(1918)
      MTR358(2)=2D0*AAABR(1918)
      MTR358(3)=AAABR(1919)
      MTR359(1)=4D0*AAABR(1920)
      MTR359(2)=4D0*AAABR(1920)
      MTR359(3)=-AAABR(1921)
      MTR360(1)=AAABR(1928)
      MTR360(2)=AAABR(1926)
      MTR360(3)=AAABR(1927)
      MTR361(1)=AAABR(1932)
      MTR361(2)=AAABR(1931)
      MTR361(3)=AAABR(1929)
      MTR362(1)=AAABR(1936)
      MTR362(2)=AAABR(1938)
      MTR362(3)=AAABR(1937)
      MTR363(1)=AAABR(1933)
      MTR363(2)=AAABR(1934)
      MTR363(3)=AAABR(1935)
      MTR364(1)=AAABR(1941)
      MTR364(2)=AAABR(1939)
      MTR364(3)=AAABR(1940)
      MTR365(1)=AAABR(1942)
      MTR365(2)=AAABR(1944)
      MTR365(3)=AAABR(1943)
      MTR366(1)=AAABR(1949)
      MTR366(2)=AAABR(1950)
      MTR366(3)=AAABR(1948)
      MTR367(1)=AAABR(1946)
      MTR367(2)=AAABR(1947)
      MTR367(3)=AAABR(1945)
      MTR368(1)=AAABR(1954)
      MTR368(2)=AAABR(1952)
      MTR368(3)=AAABR(1956)
      MTR369(1)=AAABR(1958)
      MTR369(2)=AAABR(1957)
      MTR369(3)=AAABR(1959)
      MTR370(1)=AAABR(1953)
      MTR370(2)=AAABR(1951)
      MTR370(3)=AAABR(1955)
      MTR371(1)=AAABR(1961)
      MTR371(2)=AAABR(1963)
      MTR371(3)=AAABR(1965)
      MTR372(1)=AAABR(1966)
      MTR372(2)=AAABR(1968)
      MTR372(3)=AAABR(1967)
      MTR373(1)=AAABR(1960)
      MTR373(2)=AAABR(1962)
      MTR373(3)=AAABR(1964)
      MTR374(1)=AAABR(1969)
      MTR374(2)=AAABR(1971)
      MTR374(3)=AAABR(1970)
      MTR375(1)=AAABR(1973)
      MTR375(2)=AAABR(1972)
      MTR375(3)=AAABR(1974)
      MTR376(1)=AAABR(1978)
      MTR376(2)=AAABR(1980)
      MTR376(3)=AAABR(1979)
      MTR377(1)=AAABR(1975)
      MTR377(2)=AAABR(1977)
      MTR377(3)=AAABR(1976)
      MTR378(1)=AAABR(1986)
      MTR378(2)=AAABR(1982)
      MTR378(3)=AAABR(1984)
      MTR379(1)=AAABR(1987)
      MTR379(2)=AAABR(1989)
      MTR379(3)=AAABR(1988)
      MTR380(1)=AAABR(1985)
      MTR380(2)=AAABR(1981)
      MTR380(3)=AAABR(1983)
      MTR381(1)=AAABR(1991)
      MTR381(2)=AAABR(1993)
      MTR381(3)=AAABR(1995)
      MTR382(1)=AAABR(1998)
      MTR382(2)=AAABR(1996)
      MTR382(3)=AAABR(1997)
      MTR383(1)=AAABR(1990)
      MTR383(2)=AAABR(1992)
      MTR383(3)=AAABR(1994)
      MTR384(1)=AAABR(2001)
      MTR384(2)=AAABR(2000)
      MTR384(3)=AAABR(1999)
      MTR385(1)=AAABR(2002)
      MTR385(2)=AAABR(2004)
      MTR385(3)=AAABR(2003)
      MTR386(1)=AAABR(2009)
      MTR386(2)=AAABR(2010)
      MTR386(3)=AAABR(2008)
      MTR387(1)=AAABR(2007)
      MTR387(2)=AAABR(2005)
      MTR387(3)=AAABR(2006)
      MTR388(1)=AAABR(2012)
      MTR388(2)=AAABR(2013)
      MTR388(3)=AAABR(2011)
      MTR389(1)=AAABR(2016)
      MTR389(2)=AAABR(2014)
      MTR389(3)=AAABR(2015)
      MTR390(1)=AAABR(2022)
      MTR390(2)=AAABR(2020)
      MTR390(3)=AAABR(2021)
      MTR391(1)=AAABR(2018)
      MTR391(2)=AAABR(2017)
      MTR391(3)=AAABR(2019)
      MTR392(1)=AAABR(2025)
      MTR392(2)=AAABR(2023)
      MTR392(3)=AAABR(2024)
      MTR393(1)=AAABR(2028)
      MTR393(2)=AAABR(2026)
      MTR393(3)=AAABR(2027)
      MTR394(1)=AAABR(2033)
      MTR394(2)=AAABR(2034)
      MTR394(3)=AAABR(2032)
      MTR395(1)=AAABR(2030)
      MTR395(2)=AAABR(2031)
      MTR395(3)=AAABR(2029)
      MTR396(1)=AAABR(2040)
      MTR396(2)=AAABR(2042)
      MTR396(3)=AAABR(2038)
      MTR397(1)=AAABR(2039)
      MTR397(2)=AAABR(2041)
      MTR397(3)=AAABR(2037)
      MTR398(1,1)=AAABR(2046)
      MTR398(1,2)=AAABR(2047)
      MTR398(2,1)=MTR398(1,2)
      MTR398(2,2)=AAABR(2048)
      MTR399(1,1)=AAABR(2043)
      MTR399(1,2)=AAABR(2044)
      MTR399(2,1)=MTR399(1,2)
      MTR399(2,2)=AAABR(2045)
      MTR400(1,1)=AAABR(2058)
      MTR400(1,2)=AAABR(2057)
      MTR400(1,3)=AAABR(2056)
      MTR400(2,1)=MTR400(1,2)
      MTR400(2,2)=AAABR(2055)
      MTR400(2,3)=AAABR(2060)
      MTR400(3,1)=MTR400(1,3)
      MTR400(3,2)=MTR400(2,3)
      MTR400(3,3)=AAABR(2059)
      MTR401(1,1)=AAABR(2054)
      MTR401(1,2)=AAABR(2050)
      MTR401(1,3)=AAABR(2052)
      MTR401(2,1)=MTR401(1,2)
      MTR401(2,2)=AAABR(2053)
      MTR401(2,3)=AAABR(2049)
      MTR401(3,1)=MTR401(1,3)
      MTR401(3,2)=MTR401(2,3)
      MTR401(3,3)=AAABR(2051)

      end

***********************************************

      subroutine ModelVarIni(sqrtS, *)
      implicit none
      double precision sqrtS
      double precision Alfas

#include "model.h"

c      double precision ALPHAS2
c      external ALPHAS2

c      Alfas = ALPHAS2(sqrtS)
c      GG = sqrt(4*pi*Alfas)
      end

************************************************

      subroutine ModelDigest
      implicit none

#include "model.h"

      write(16,fmt='(A10,1PE12.4  )') 'EE=',EE
      write(16,fmt='(A10,1PE12.4  )') 'MZ=',MZ
      write(16,fmt='(A10,1PE12.4  )') 'MW=',MW
      write(16,fmt='(A10,1PE12.4  )') 'GG=',GG
      write(16,fmt='(A10,1PE12.4  )') 'alfSMZ=',alfSMZ
      write(16,fmt='(A10,1PE12.4  )') 'MbMb=',MbMb
      write(16,fmt='(A10,1PE12.4  )') 'McMc=',McMc
      write(16,fmt='(A10,1PE12.4  )') 'Mtp=',Mtp
      write(16,fmt='(A10,1PE12.4  )') 'Ml=',Ml
      write(16,fmt='(A10,1PE12.4  )') 'Mq=',Mq
      write(16,fmt='(A10,1PE12.4  )') 'CW=',CW
      write(16,fmt='(A10,1PE12.4  )') 'SW=',SW
      write(16,fmt='(A10,1PE12.4  )') 'C2W=',C2W
      write(16,fmt='(A10,1PE12.4  )') 'LamQCD=',LamQCD
      write(16,fmt='(A10,1PE12.4  )') 'rd=',rd
      write(16,fmt='(A10,1PE12.4  )') 'mu=',mu
      write(16,fmt='(A10,1PE12.4  )') 'Lambda=',Lambda
      write(16,fmt='(A10,1PE12.4  )') 'Kappa=',Kappa
      write(16,fmt='(A10,1PE12.4  )') 'aLambda=',aLambda
      write(16,fmt='(A10,1PE12.4  )') 'aKappa=',aKappa
      write(16,fmt='(A10,1PE12.4  )') 'muP=',muP
      write(16,fmt='(A10,1PE12.4  )') 'At=',At
      write(16,fmt='(A10,1PE12.4  )') 'Ab=',Ab
      write(16,fmt='(A10,1PE12.4  )') 'Al=',Al
      write(16,fmt='(A10,1PE12.4  )') 'Ml2=',Ml2
      write(16,fmt='(A10,1PE12.4  )') 'Ml3=',Ml3
      write(16,fmt='(A10,1PE12.4  )') 'Mr2=',Mr2
      write(16,fmt='(A10,1PE12.4  )') 'Mr3=',Mr3
      write(16,fmt='(A10,1PE12.4  )') 'Mq2=',Mq2
      write(16,fmt='(A10,1PE12.4  )') 'Mq3=',Mq3
      write(16,fmt='(A10,1PE12.4  )') 'Mu2=',Mu2
      write(16,fmt='(A10,1PE12.4  )') 'Mu3=',Mu3
      write(16,fmt='(A10,1PE12.4  )') 'Md2=',Md2
      write(16,fmt='(A10,1PE12.4  )') 'Md3=',Md3
      write(16,fmt='(A10,1PE12.4  )') 'tB=',tB
      write(16,fmt='(A10,1PE12.4  )') 'Mh1=',Mh1
      write(16,fmt='(A10,1PE12.4  )') 'Mh2=',Mh2
      write(16,fmt='(A10,1PE12.4  )') 'Mh3=',Mh3
      write(16,fmt='(A10,1PE12.4  )') 'Mha=',Mha
      write(16,fmt='(A10,1PE12.4  )') 'Mhb=',Mhb
      write(16,fmt='(A10,1PE12.4  )') 'MHc=',MHc
      write(16,fmt='(A10,1PE12.4  )') 'MNE1=',MNE1
      write(16,fmt='(A10,1PE12.4  )') 'MNE2=',MNE2
      write(16,fmt='(A10,1PE12.4  )') 'MNE3=',MNE3
      write(16,fmt='(A10,1PE12.4  )') 'MNE4=',MNE4
      write(16,fmt='(A10,1PE12.4  )') 'MNE5=',MNE5
      write(16,fmt='(A10,1PE12.4  )') 'MC1=',MC1
      write(16,fmt='(A10,1PE12.4  )') 'MC2=',MC2
      write(16,fmt='(A10,1PE12.4  )') 'MSG=',MSG
      write(16,fmt='(A10,1PE12.4  )') 'MSuL=',MSuL
      write(16,fmt='(A10,1PE12.4  )') 'MSuR=',MSuR
      write(16,fmt='(A10,1PE12.4  )') 'MSdL=',MSdL
      write(16,fmt='(A10,1PE12.4  )') 'MSdR=',MSdR
      write(16,fmt='(A10,1PE12.4  )') 'MScL=',MScL
      write(16,fmt='(A10,1PE12.4  )') 'MScR=',MScR
      write(16,fmt='(A10,1PE12.4  )') 'MSsL=',MSsL
      write(16,fmt='(A10,1PE12.4  )') 'MSsR=',MSsR
      write(16,fmt='(A10,1PE12.4  )') 'MSt1=',MSt1
      write(16,fmt='(A10,1PE12.4  )') 'MSt2=',MSt2
      write(16,fmt='(A10,1PE12.4  )') 'MSb1=',MSb1
      write(16,fmt='(A10,1PE12.4  )') 'MSb2=',MSb2
      write(16,fmt='(A10,1PE12.4  )') 'QSUSY=',QSUSY
      write(16,fmt='(A10,1PE12.4  )') 'MSeL=',MSeL
      write(16,fmt='(A10,1PE12.4  )') 'MSeR=',MSeR
      write(16,fmt='(A10,1PE12.4  )') 'MSmL=',MSmL
      write(16,fmt='(A10,1PE12.4  )') 'MSmR=',MSmR
      write(16,fmt='(A10,1PE12.4  )') 'MSl1=',MSl1
      write(16,fmt='(A10,1PE12.4  )') 'MSl2=',MSl2
      write(16,fmt='(A10,1PE12.4  )') 'MSne=',MSne
      write(16,fmt='(A10,1PE12.4  )') 'MSnm=',MSnm
      write(16,fmt='(A10,1PE12.4  )') 'MSnl=',MSnl
      write(16,fmt='(A10,1PE12.4  )') 'Zl11=',Zl11
      write(16,fmt='(A10,1PE12.4  )') 'Zl12=',Zl12
      write(16,fmt='(A10,1PE12.4  )') 'Zl21=',Zl21
      write(16,fmt='(A10,1PE12.4  )') 'Zl22=',Zl22
      write(16,fmt='(A10,1PE12.4  )') 'Zh11=',Zh11
      write(16,fmt='(A10,1PE12.4  )') 'Zh12=',Zh12
      write(16,fmt='(A10,1PE12.4  )') 'Zh13=',Zh13
      write(16,fmt='(A10,1PE12.4  )') 'Zh21=',Zh21
      write(16,fmt='(A10,1PE12.4  )') 'Zh22=',Zh22
      write(16,fmt='(A10,1PE12.4  )') 'Zh23=',Zh23
      write(16,fmt='(A10,1PE12.4  )') 'Zh31=',Zh31
      write(16,fmt='(A10,1PE12.4  )') 'Zh32=',Zh32
      write(16,fmt='(A10,1PE12.4  )') 'Zh33=',Zh33
      write(16,fmt='(A10,1PE12.4  )') 'Za11=',Za11
      write(16,fmt='(A10,1PE12.4  )') 'Za12=',Za12
      write(16,fmt='(A10,1PE12.4  )') 'Za13=',Za13
      write(16,fmt='(A10,1PE12.4  )') 'Za21=',Za21
      write(16,fmt='(A10,1PE12.4  )') 'Za22=',Za22
      write(16,fmt='(A10,1PE12.4  )') 'Za23=',Za23
      write(16,fmt='(A10,1PE12.4  )') 'Zn11=',Zn11
      write(16,fmt='(A10,1PE12.4  )') 'Zn12=',Zn12
      write(16,fmt='(A10,1PE12.4  )') 'Zn13=',Zn13
      write(16,fmt='(A10,1PE12.4  )') 'Zn14=',Zn14
      write(16,fmt='(A10,1PE12.4  )') 'Zn15=',Zn15
      write(16,fmt='(A10,1PE12.4  )') 'Zn21=',Zn21
      write(16,fmt='(A10,1PE12.4  )') 'Zn22=',Zn22
      write(16,fmt='(A10,1PE12.4  )') 'Zn23=',Zn23
      write(16,fmt='(A10,1PE12.4  )') 'Zn24=',Zn24
      write(16,fmt='(A10,1PE12.4  )') 'Zn25=',Zn25
      write(16,fmt='(A10,1PE12.4  )') 'Zn31=',Zn31
      write(16,fmt='(A10,1PE12.4  )') 'Zn32=',Zn32
      write(16,fmt='(A10,1PE12.4  )') 'Zn33=',Zn33
      write(16,fmt='(A10,1PE12.4  )') 'Zn34=',Zn34
      write(16,fmt='(A10,1PE12.4  )') 'Zn35=',Zn35
      write(16,fmt='(A10,1PE12.4  )') 'Zn41=',Zn41
      write(16,fmt='(A10,1PE12.4  )') 'Zn42=',Zn42
      write(16,fmt='(A10,1PE12.4  )') 'Zn43=',Zn43
      write(16,fmt='(A10,1PE12.4  )') 'Zn44=',Zn44
      write(16,fmt='(A10,1PE12.4  )') 'Zn45=',Zn45
      write(16,fmt='(A10,1PE12.4  )') 'Zn51=',Zn51
      write(16,fmt='(A10,1PE12.4  )') 'Zn52=',Zn52
      write(16,fmt='(A10,1PE12.4  )') 'Zn53=',Zn53
      write(16,fmt='(A10,1PE12.4  )') 'Zn54=',Zn54
      write(16,fmt='(A10,1PE12.4  )') 'Zn55=',Zn55
      write(16,fmt='(A10,1PE12.4  )') 'Zu11=',Zu11
      write(16,fmt='(A10,1PE12.4  )') 'Zu12=',Zu12
      write(16,fmt='(A10,1PE12.4  )') 'Zu21=',Zu21
      write(16,fmt='(A10,1PE12.4  )') 'Zu22=',Zu22
      write(16,fmt='(A10,1PE12.4  )') 'Zv11=',Zv11
      write(16,fmt='(A10,1PE12.4  )') 'Zv12=',Zv12
      write(16,fmt='(A10,1PE12.4  )') 'Zv21=',Zv21
      write(16,fmt='(A10,1PE12.4  )') 'Zv22=',Zv22
      write(16,fmt='(A10,1PE12.4  )') 'Zt11=',Zt11
      write(16,fmt='(A10,1PE12.4  )') 'Zt12=',Zt12
      write(16,fmt='(A10,1PE12.4  )') 'Zt21=',Zt21
      write(16,fmt='(A10,1PE12.4  )') 'Zt22=',Zt22
      write(16,fmt='(A10,1PE12.4  )') 'Zb11=',Zb11
      write(16,fmt='(A10,1PE12.4  )') 'Zb12=',Zb12
      write(16,fmt='(A10,1PE12.4  )') 'Zb21=',Zb21
      write(16,fmt='(A10,1PE12.4  )') 'Zb22=',Zb22
      write(16,fmt='(A10,1PE12.4  )') 'la1=',la1
      write(16,fmt='(A10,1PE12.4  )') 'la2=',la2
      write(16,fmt='(A10,1PE12.4  )') 'la3=',la3
      write(16,fmt='(A10,1PE12.4  )') 'la4=',la4
      write(16,fmt='(A10,1PE12.4  )') 'la5=',la5
      write(16,fmt='(A10,1PE12.4  )') 'la6=',la6
      write(16,fmt='(A10,1PE12.4  )') 'la7=',la7
      write(16,fmt='(A10,1PE12.4  )') 'la1s=',la1s
      write(16,fmt='(A10,1PE12.4  )') 'la2s=',la2s
      write(16,fmt='(A10,1PE12.4  )') 'la3s=',la3s
      write(16,fmt='(A10,1PE12.4  )') 'la4s=',la4s
      write(16,fmt='(A10,1PE12.4  )') 'la5s=',la5s
      write(16,fmt='(A10,1PE12.4  )') 'la6s=',la6s
      write(16,fmt='(A10,1PE12.4  )') 'la7s=',la7s
      write(16,fmt='(A10,1PE12.4  )') 'la8s=',la8s
      write(16,fmt='(A10,1PE12.4  )') 'aa1=',aa1
      write(16,fmt='(A10,1PE12.4  )') 'aa2=',aa2
      write(16,fmt='(A10,1PE12.4  )') 'aa3=',aa3
      write(16,fmt='(A10,1PE12.4  )') 'aa4=',aa4
      write(16,fmt='(A10,1PE12.4  )') 'aa5=',aa5
      write(16,fmt='(A10,1PE12.4  )') 'aa6=',aa6
      write(16,fmt='(A10,1PE12.4  )') 'B1=',mB1
      write(16,fmt='(A10,1PE12.4  )') 'B2=',mB2
      write(16,fmt='(A10,1PE12.4  )') 'X=',X
      write(16,fmt='(A10,1PE12.4  )') 'dMb=',dMb
      write(16,fmt='(A10,1PE12.4  )') 'xif=',xif
      write(16,fmt='(A10,1PE12.4  )') 'xis=',xis
      write(16,fmt='(A10,1PE12.4  )') 'MM3=',MM3
      write(16,fmt='(A10,1PE12.4  )') 'MSP=',MSP
C     write(16,fmt='(A10,1PE12.4  )') 'sb=',sb
C     write(16,fmt='(A10,1PE12.4  )') 'cb=',cb
C     write(16,fmt='(A10,1PE12.4  )') 't2b=',t2b
C     write(16,fmt='(A10,1PE12.4  )') 'xvev=',xvev
C     write(16,fmt='(A10,1PE12.4  )') 'Pa12=',Pa12
C     write(16,fmt='(A10,1PE12.4  )') 'Pa22=',Pa22
C     write(16,fmt='(A10,1PE12.4  )') 'Pa11=',Pa11
C     write(16,fmt='(A10,1PE12.4  )') 'Pa21=',Pa21
C     write(16,fmt='(A10,1PE12.4  )') 'Td3=',Td3
C     write(16,fmt='(A10,1PE12.4  )') 'Q=',Q
C     write(16,fmt='(A10,1PE12.4  )') 'Mt=',Mt
C     write(16,fmt='(A10,1PE12.4  )') 'Mb=',Mb
C     write(16,fmt='(A10,1PE12.4  )') 'Mc=',Mc
C     write(16,fmt='(A10,1PE12.4  )') 'xH2=',xH2
C     write(16,fmt='(A10,1PE12.4  )') 'yH2=',yH2
C     write(16,fmt='(A10,1PE12.4  )') 'dMd=',dMd
C     write(16,fmt='(A10,1PE12.4  )') 'Td2=',Td2
C     write(16,fmt='(A10,1PE12.4  )') 'Au=',Au
C     write(16,fmt='(A10,1PE12.4  )') 'Ad=',Ad
C     write(16,fmt='(A10,1PE12.4  )') 'fiuu=',fiuu
C     write(16,fmt='(A10,1PE12.4  )') 'fidd=',fidd
C     write(16,fmt='(A10,1PE12.4  )') 'ficc=',ficc
C     write(16,fmt='(A10,1PE12.4  )') 'fiss=',fiss
C     write(16,fmt='(A10,1PE12.4  )') 'Zuu11=',Zuu11
C     write(16,fmt='(A10,1PE12.4  )') 'Zuu12=',Zuu12
C     write(16,fmt='(A10,1PE12.4  )') 'Zuu21=',Zuu21
C     write(16,fmt='(A10,1PE12.4  )') 'Zuu22=',Zuu22
C     write(16,fmt='(A10,1PE12.4  )') 'Zdd11=',Zdd11
C     write(16,fmt='(A10,1PE12.4  )') 'Zdd12=',Zdd12
C     write(16,fmt='(A10,1PE12.4  )') 'Zdd21=',Zdd21
C     write(16,fmt='(A10,1PE12.4  )') 'Zdd22=',Zdd22
C     write(16,fmt='(A10,1PE12.4  )') 'Zcc11=',Zcc11
C     write(16,fmt='(A10,1PE12.4  )') 'Zcc12=',Zcc12
C     write(16,fmt='(A10,1PE12.4  )') 'Zcc21=',Zcc21
C     write(16,fmt='(A10,1PE12.4  )') 'Zcc22=',Zcc22
C     write(16,fmt='(A10,1PE12.4  )') 'Zss11=',Zss11
C     write(16,fmt='(A10,1PE12.4  )') 'Zss12=',Zss12
C     write(16,fmt='(A10,1PE12.4  )') 'Zss21=',Zss21
C     write(16,fmt='(A10,1PE12.4  )') 'Zss22=',Zss22
C     write(16,fmt='(A10,1PE12.4  )') 'Mele=',Mele
C     write(16,fmt='(A10,1PE12.4  )') 'Mmu=',Mmu
C     write(16,fmt='(A10,1PE12.4  )') 'MbMM=',MbMM
C     write(16,fmt='(A10,1PE12.4  )') 'MtMM=',MtMM
C     write(16,fmt='(A10,1PE12.4  )') 'NMM11=',NMM11
C     write(16,fmt='(A10,1PE12.4  )') 'NMM12=',NMM12
C     write(16,fmt='(A10,1PE12.4  )') 'NMM13=',NMM13
C     write(16,fmt='(A10,1PE12.4  )') 'NMM14=',NMM14
C     write(16,fmt='(A10,1PE12.4  )') 'NMM15=',NMM15
C     write(16,fmt='(A10,1PE12.4  )') 'NMM22=',NMM22
C     write(16,fmt='(A10,1PE12.4  )') 'NMM23=',NMM23
C     write(16,fmt='(A10,1PE12.4  )') 'NMM24=',NMM24
C     write(16,fmt='(A10,1PE12.4  )') 'NMM25=',NMM25
C     write(16,fmt='(A10,1PE12.4  )') 'NMM33=',NMM33
C     write(16,fmt='(A10,1PE12.4  )') 'NMM34=',NMM34
C     write(16,fmt='(A10,1PE12.4  )') 'NMM35=',NMM35
C     write(16,fmt='(A10,1PE12.4  )') 'NMM44=',NMM44
C     write(16,fmt='(A10,1PE12.4  )') 'NMM45=',NMM45
C     write(16,fmt='(A10,1PE12.4  )') 'NMM55=',NMM55
C     write(16,fmt='(A10,1PE12.4  )') 'MG1I=',MG1I
C     write(16,fmt='(A10,1PE12.4  )') 'MG2I=',MG2I

      end

