*     LanHEP output produced at Tue Apr 19 21:58:11 2016
*     Model named 'Inert Dublet'

      subroutine ModelDefaults

      implicit none


#include "model.h"

      EE = 0.31223D0
      GG = 1.117D0
      alfSMZ = 0.1184D0
      SW = 0.481D0
      MZ = 91.187D0
      Q = 100D0
      Mtp = 173.5D0
      MbMb = 4.23D0
      McMc = 1.27D0
      MHX = 111D0
      MH3 = 222D0
      MHC = 333D0
      laL = 0.01D0
      la2 = 0.01D0
      wZ = 2.502D0
      wW = 2.094D0
      Mm = 0.1057D0
      Ml = 1.777D0
      Mqu = 0.01D0
      Mqd = 0.01D0
      Ms = 0.2D0
      Mh = 200D0
      wh =0 
      end


      subroutine ModelConstIni(fail)
      implicit none
      integer fail

#include "model.h"

      double precision 
     > initQCD5,McRun,MbRun,MtRun,alphaQCD,
     > McEff,MbEff,MtEff

      CW = sqrt(1D0-SW**2)
      MW = MZ*CW
      Lqcd = initQCD5(alfSMZ, McMc, MbMb, Mtp)
      Mb = MbEff(Q)
      Mt = MtEff(Q)
      Mc = McEff(Q)
      Mbp = MbMb*(1D0+4D0/3D0*alphaQCD(MbMb)/PI)
      Mcp = Mc*(1D0+4D0/3D0*alphaQCD(Mc)/PI)
      mu2 = MHX**2-laL*(2D0*MW/EE*SW)**2
      la3 = 2D0*(MHC**2-mu2)/(2D0*MW/EE*SW)**2
      la5 = (MHX**2-MH3**2)/(2D0*MW/EE*SW)**2
      la4 = 2D0*laL-la3-la5
      Mqu2 = Mqu**2
      Mqd2 = Mqd**2
      Mc2 = Mc**2
      Ms2 = Ms**2
      Mt2 = Mt**2
      Mb2 = Mb**2
      Mm2 = Mm**2
      Ml2 = Ml**2
      Mh2 = Mh**2
      MH32 = MH3**2
      MHC2 = MHC**2
      MHX2 = MHX**2
      EE2 = EE**2
      MW2 = MW**2
      MZ2 = MZ**2

      call aaini01
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
      AAABR(7) = EE/CW*MW*(1D0/SW-2D0*SW)
      AAABR(8) = EE/CW*MW/SW
      AAABR(9) = EE/CW**2*MW/SW
      AAABR(10) = CW*EE/SW
      AAABR(11) = EE/MW*Mh**2/SW
      AAABR(12) = MW/EE*SW*(la4-la5)
      AAABR(13) = MW/EE*SW*(la4+la5)
      AAABR(14) = MW/EE*SW*la5
      AAABR(15) = MW/EE*SW*la3
      AAABR(16) = MW/EE*SW*(la3+la4-la5)
      AAABR(17) = MW/EE*SW*(la3+la4+la5)
      AAABR(18) = EE*(SW/CW-CW/SW)
      AAABR(19) = EE/SW
      AAABR(20) = EE*(SW/CW+CW/SW)
      AAABR(21) = EE/CW*MW*SW
      AAABR(22) = EE*MW*(2D0*SW+SW**3/CW**2+CW**2/SW)
      AAABR(23) = EE/MW*Mb/SW
      AAABR(24) = EE/MW*Mb/SW*Sqrt2
      AAABR(25) = EE/MW*Mt/SW*Sqrt2
      AAABR(26) = EE/MW*Mc/SW
      AAABR(27) = EE/MW*Ms/SW*Sqrt2
      AAABR(28) = EE/MW*Mc/SW*Sqrt2
      AAABR(29) = EE/MW*Mqd/SW
      AAABR(30) = EE/MW*Mqd/SW*Sqrt2
      AAABR(31) = EE/MW*Mqu/SW*Sqrt2
      AAABR(32) = EE/MW*Ml/SW
      AAABR(33) = EE/MW*Ml/SW*Sqrt2
      AAABR(34) = EE/MW*Mm/SW
      AAABR(35) = EE/MW*Mm/SW*Sqrt2
      AAABR(36) = EE/MW*Ms/SW
      AAABR(37) = EE/MW*Mt/SW
      AAABR(38) = EE/MW*Mqu/SW
      AAABR(39) = EE*(SW/CW+3D0*CW/SW)
      AAABR(40) = EE/CW*SW
      AAABR(41) = EE/SW*Sqrt2
      AAABR(42) = EE*(SW/CW-3D0*CW/SW)
      AAABR(43) = EE**2/MW**2*Mh**2/SW**2
      AAABR(44) = la3+la4
      AAABR(45) = la4+la5
      AAABR(46) = la4-la5
      AAABR(47) = la3+la4+la5
      AAABR(48) = la3+la4-la5
      AAABR(49) = EE**2*(SW/CW-CW/SW)
      AAABR(50) = EE**2
      AAABR(51) = EE**2*(2D0-SW**2/CW**2-CW**2/SW**2)
      AAABR(52) = EE**2/SW**2
      AAABR(53) = EE**2/CW
      AAABR(54) = EE**2/SW
      AAABR(55) = EE**2*(2D0+SW**2/CW**2+CW**2/SW**2)
      AAABR(56) = CW*EE**2/SW
      AAABR(57) = CW**2*EE**2/SW**2
      AAABR(58) = GG**2
      end

      subroutine mtrini
      implicit none
#include "model.h"

      integer m1,m2,m3,m4


      end

***********************************************

      subroutine ModelVarIni(nFile)
      implicit none
      double precision sqrtS
      integer nFile
      double precision Alfas

#include "model.h"
       character*10 name
       real*8  val

c      double precision ALPHAS2
c      external ALPHAS2

c      Alfas = ALPHAS2(sqrtS)
c      GG = sqrt(4*pi*Alfas)

      call ModelDefaults

123    continue         
       read(nFile,*,end=321) name, val
        if(name.eq."EE") then
            EE=val
        else if(name.eq."GG") then
            GG=val
        else if(name.eq."alfSMZ") then
            alfSMZ=val
        else if(name.eq."SW") then
            SW=val
        else if(name.eq."MZ") then
            MZ=val
        else if(name.eq."Q") then
            Q=val
        else if(name.eq."Mtp") then
            Mtp=val
        else if(name.eq."MbMb") then
            MbMb=val
        else if(name.eq."McMc") then
            McMc=val
        else if(name.eq."MHX") then
            MHX=val
        else if(name.eq."MH3") then
            MH3=val
        else if(name.eq."MHC") then
            MHC=val
        else if(name.eq."laL") then
            laL=val
        else if(name.eq."la2") then
            la2=val
        else if(name.eq."wZ") then
            wZ=val
        else if(name.eq."wW") then
            wW=val
        else if(name.eq."Mm") then
            Mm=val
        else if(name.eq."Ml") then
            Ml=val
        else if(name.eq."Mu") then
            Mqu=val
        else if(name.eq."Md") then
            Mqd=val
        else if(name.eq."Ms") then
            Ms=val
        else if(name.eq."Mh") then
            Mh=val
        else if(name.eq."wh") then
            wh=val
        endif        
       goto 123
321    continue
      end

************************************************


