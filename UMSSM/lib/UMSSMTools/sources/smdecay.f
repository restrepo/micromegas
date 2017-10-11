        SUBROUTINE SMDECAY()

        IMPLICIT NONE

        INTEGER N0,NF,NFEXT,NFGG,I,VFLAG,HFLAG

        DOUBLE PRECISION MH,S(3,3),SMASS(3)
        DOUBLE PRECISION PMASS,P(3,3),CMASS
        DOUBLE PRECISION CJ,CG,PI,HIGTOP,ASG,ASH,AS3,AS4,ASMT
        DOUBLE PRECISION SQR2,EPS,FQCD,XFAC,X,Y,RATCOUP,RAT
        DOUBLE PRECISION HJJ,HEE,HMM,HLL,HSS,HCC,HBB,HTT,HWW,HZZ,HGG,HZG
        DOUBLE PRECISION HS1,HS2,HC1,HC2,HB1,HB2,HT1,HT2,DCC,DBB
        DOUBLE PRECISION DLU,DLD,XM1,XM2,CWW,CZZ,XX(4),YY(4)
        DOUBLE PRECISION WIDTHSM(3),BRJJSM(3),BREESM(3),BRMMSM(3)
        DOUBLE PRECISION BRLLSM(3),HTWW,HTZZ,MPI,MEL
        DOUBLE PRECISION BRSSSM(3),BRCCSM(3),BRBBSM(3),BRTTSM(3)
        DOUBLE PRECISION BRWWSM(3),BRZZSM(3),BRGGSM(3),BRZGSM(3)
        DOUBLE PRECISION XLAMBDA,MC0,MB0,MT0,RMS,RMC
        DOUBLE PRECISION RMB,RMT,C2TW,T2TW,ALEM0,RMTTOP,FT,FB,RUNMB
        DOUBLE PRECISION HVV,HV,HFF,QCD0,HQCDM,HQCD,QCDH,TQCDH,HGGQCD
        DOUBLE PRECISION BETA,SP,ALPHAS,RUNM,QQINT,FINT,ACOUP
        DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
        DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW

        DOUBLE COMPLEX CTT,CTB,CTC,CTS,CTL,CTM,CTW
        DOUBLE COMPLEX CXT,CXB,CXC,CXS,CXL,CXM,CXW
        DOUBLE COMPLEX CLT,CLB,CLC,CLW,CXTZ,CXBZ,CXCZ,CXWZ
        DOUBLE COMPLEX CI1,CI2,CGZ,CF,CA,CB

        COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
        COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        COMMON/ALS/XLAMBDA,MC0,MB0,MT0,N0
        COMMON/BRSM/WIDTHSM,BRJJSM,BREESM,BRMMSM,BRLLSM,BRSSSM,BRCCSM,
     .      BRBBSM,BRTTSM,BRWWSM,BRZZSM,BRGGSM,BRZGSM 
        COMMON/HIGGSPEC/SMASS,S,PMASS,P,CMASS
        COMMON/FLAGS/VFLAG,HFLAG
        COMMON/SMEXT/MPI,MEL

        QQINT(RAT,X,Y)= RAT**2*X+(1.D0-RAT**2)*Y
        BETA(X)= DSQRT(1.D0-4.D0*X)
        CF(CA)= -CDLOG(-(1.D0+CDSQRT(1.D0-CA))
     .   / (1.D0-CDSQRT(1.D0-CA)))**2/4.D0
        CGZ(CA)= CDSQRT(1.D0-CA)/2.D0*CDLOG(-(1.D0+CDSQRT(1.D0-CA))
     .   / (1.D0-CDSQRT(1.D0-CA)))
        CI1(CA,CB)= CA*CB/2.D0/(CA-CB)
     .   + CA**2*CB**2/2/(CA-CB)**2*(CF(CA)-CF(CB))
     .   + CA**2*CB/(CA-CB)**2*(CGZ(CA)-CGZ(CB))
        CI2(CA,CB)= -CA*CB/2.D0/(CA-CB)*(CF(CA)-CF(CB))
        HV(X)= 3.D0*(1.D0-8.D0*X+20.D0*X**2)/DSQRT((4.D0*X-1.D0))
     .   * DACOS((3.D0*X-1.D0)/2.D0/DSQRT(X**3))
     .   - (1.D0-X)*(47.D0/2.D0*X-13.D0/2.D0+1.D0/X)
     .   - 3.D0/2.D0*(1.D0-6.D0*X+4.D0*X**2)*DLOG(X)
        HVV(X,Y)= GF/(4.D0*PI*SQR2)*X**3/2.D0*BETA(Y)
     .   * (1.D0-4.D0*Y+12.D0*Y**2)
        HFF(X,Y)= GF/(4.D0*PI*SQR2)*X**3*Y*(BETA(Y))**3
        QCD0(X)= (1.D0+X**2)*(4.D0*SP((1.D0-X)/(1.D0+X))
     .   + 2.D0*SP((X-1.D0)/(X+1.D0))
     .   - 3.D0*DLOG((1.D0+X)/(1.D0-X))*DLOG(2.D0/(1.D0+X))
     .   - 2.D0*DLOG((1.D0+X)/(1.D0-X))*DLOG(X))
     .   - 3.D0*X*DLOG(4.D0/(1.D0-X**2))-4.D0*X*DLOG(X)
        HQCDM(X)= QCD0(X)/X+(3.D0+34.D0*X**2-13.D0*X**4)/16.D0/X**3
     .   * DLOG((1.D0+X)/(1.D0-X))+3.D0/8.D0/X**2*(7.D0*X**2-1.D0)
        HQCD(X)= 5.67D0*ASH/PI
     .   + (29.14D0+RATCOUP*(1.57D0-2.D0*DLOG(HIGTOP)/3.D0
     .   + DLOG(X)**2/9.D0))*(ASH/PI)**2
     .   + (164.14D0-25.77D0*5.D0+0.259D0*5.D0**2)*(ASH/PI)**3
        QCDH(X)= 1.D0+HQCD(X)
        TQCDH(X)= 1.D0+4.D0/3.D0*HQCDM(BETA(X))*ASH/PI
        HGGQCD(ASG,NF)= 1.D0+ASG/PI*(95.D0/4.D0-NF*7.D0/6.D0)

        EPS= 1.D-8
        PI= 4.D0*DATAN(1.D0)
        SQR2= DSQRT(2.D0)

*   Number of light flavours included in the gluonic decays
*   Higgs -> gg* -> gqq (see hdecay): NFGG = 3
        NFGG= 3

*   Alpha_EM(0)
        ALEM0= 1.D0/137.04D0

*   Weak angle theta_W (S2TW = sin(theta_W)):
        C2TW= 1.D0-S2TW
        T2TW= S2TW/C2TW

*   Alpha_s at the top pole mass scales, used for the running
*   Yukawa coupling ht and running quark masses RMT below
*   NOTE: MT = top pole mass
        ASMT= ALPHAS(MT,2)

*   MT = Top pole mass; RMTTOP = running mass at Mtop (MS_bar):
        RMTTOP= MT/(1.D0+4.D0*ASMT/(3.D0*PI)+11.D0*(ASMT/PI)**2)

        DO I=1,3

         MH= SMASS(I)
         HIGTOP= (MAX(1.D0,MH)/MT)**2
         MT0= 3.D8
         ASH= ALPHAS(MAX(1.D0,MH),2)
         MC0= 1.D8
         MB0= 2.D8
         AS3= ALPHAS(MAX(1.D0,MH),2)
         MC0= MC
         AS4= ALPHAS(MAX(1.D0,MH),2)
         MB0= MBP
         MT0= MT

*  Running quark masses at MH

         RMS= RUNM(MAX(1.D0,MH),3)
 
         RMC= RUNM(MAX(1.D0,MH),4)

         RMB= RUNMB(MAX(1.D0,MH))

         IF(MH.GE.MT)THEN
          RMT= RMTTOP
     .      *(1.D0+7.D0/(4.D0*PI)*ASMT*DLOG(MH**2/MT**2))
     .      **(-4.D0/7.D0)
         ELSE
          RMT= RMTTOP
     .      *(1.D0+23.D0/(12.D0*PI)*ASMT*DLOG(MAX(1d0,MH)**2/MT**2))
     .      **(-12.D0/23.D0)
         ENDIF

*  Radiative couplings

         CTT= 4.D0*(MT/MH)**2*DCMPLX(1.D0,-EPS)
         CTB= 4.D0*(MBP/MH)**2*DCMPLX(1.D0,-EPS)
         CTC= 4.D0*(MC/MH)**2*DCMPLX(1.D0,-EPS)
         CTS= 4.d0*(MS/MH)**2*DCMPLX(1.D0,-EPS)
         CTL= 4.D0*(MTAU/MH)**2*DCMPLX(1.D0,-EPS)
         CTM= 4.D0*(MMUON/MH)**2*DCMPLX(1.D0,-EPS)
         CTW= 4.D0*(MW/MH)**2*DCMPLX(1.D0,-EPS)
         CXT= 2.D0*CTT*(1.D0+(1.D0-CTT)*CF(CTT))
         CXB= 2.D0*CTB*(1.D0+(1.D0-CTB)*CF(CTB))
         CXC= 2.D0*CTC*(1.D0+(1.D0-CTC)*CF(CTC))
         CXS= 2.D0*CTS*(1.D0+(1.D0-CTS)*CF(CTS))
         CXL= 2.D0*CTL*(1.D0+(1.D0-CTL)*CF(CTL))
         CXM= 2.D0*CTM*(1.D0+(1.D0-CTM)*CF(CTM))
         CXW= -(2.D0+3.D0*CTW+3.D0*CTW*(2.D0-CTW)*CF(CTW))
         CJ= CDABS(CXT+CXC+CXB)
         CG= CDABS(4.D0/3.D0*(CXT+CXC)+CXB/3.D0+CXL+CXW)

*  Partial widths

*   h -> gg

         NFEXT= 3
         ASG= AS3
         FQCD= HGGQCD(ASG,NFEXT)
         XFAC= CJ**2*FQCD
         IF(MH.LE.2.D0*MPI)THEN
          HJJ= 0.D0
         ELSE
           HJJ= GF/(64.D0*PI*SQR2)*MH**3*(ASG/PI)**2*XFAC
         ENDIF

*   h -> gg* -> gcc to be added to h -> cc

         NFEXT= 4
         ASG= AS4
         FQCD= HGGQCD(ASG,NFEXT)
         XFAC= CJ**2*FQCD
         DCC= GF/(64.D0*PI*SQR2)*MH**3*(ASG/PI)**2*XFAC-HJJ

*   h -> gg* -> gbb to be added to h -> bb

         NFEXT= 5
         ASG= ASH
         FQCD= HGGQCD(ASG,NFEXT)
         XFAC= CJ**2*FQCD
         DBB= GF/(64.D0*PI*SQR2)*MH**3*(ASG/PI)**2*XFAC-HJJ-DCC

         IF(NFGG.EQ.5)THEN
          HJJ= HJJ+DBB+DCC
          DBB= 0.D0
          DCC= 0.D0
         ELSEIF(NFGG.EQ.4)THEN
          HJJ= HJJ+DCC
          DCC= 0.D0
         ENDIF

*   h -> ee

         IF(MH.LE.2.D0*MEL)THEN
          HEE= 0.D0
         ELSE
          HEE= HFF(MH,(MEL/MH)**2)
         ENDIF

*   h -> mumu

         IF(MH.LE.2.D0*MMUON)THEN
          HMM= 0.D0
         ELSE
          HMM= HFF(MH,(MMUON/MH)**2)
         ENDIF

*   h -> tautau

         IF(MH.LE.2.D0*MTAU)THEN
          HLL= 0.D0
         ELSE
          HLL= HFF(MH,(MTAU/MH)**2)
         ENDIF

*   h -> ss

         IF(MH.LE.2.D0*MS)THEN
          HSS= 0.D0
         ELSE
          RATCOUP= 1.D0
          HS1= 3.D0*HFF(MH,(MS/MH)**2)
     .      * TQCDH((MS/MAX(1d0,MH))**2)
         HS2= 3.D0*HFF(MH,(RMS/MH)**2)
     .      * QCDH((RMS/MAX(1d0,MH))**2)
          IF(HS2.LT.0.D0) HS2=0.D0
          RAT= 2.D0*MS/MH
          HSS= QQINT(RAT,HS1,HS2)
         ENDIF

*   h -> cc

         IF(MH.LE.2.D0*MC)THEN
          HCC= 0.D0
         ELSE
          RATCOUP= 1.D0
          HC1= 3.D0*HFF(MH,(MC/MH)**2)
     .      * TQCDH((MC/MH)**2)
          HC2= 3.D0*HFF(MH,(RMC/MH)**2)
     .      * QCDH((RMC/MH)**2)
     .      + DCC
          IF(HC2.LT.0.D0) HC2=0.D0
          RAT= 2.D0*MC/MH
          HCC= QQINT(RAT,HC1,HC2)
         ENDIF

*   h -> bb

         IF(MH.LE.2.D0*MBP)THEN
          HBB= 0.D0
         ELSE
          RATCOUP= 1.D0
          HB1= 3.D0*HFF(MH,(MBP/MH)**2)
     .      * TQCDH((MBP/MH)**2)
          HB2= 3.D0*HFF(MH,(RMB/MH)**2)
     .      * QCDH((RMB/MH)**2)
     .      + DBB
          IF(HB2.LT.0.D0) HB2=0.D0
          RAT= 2.D0*MBP/MH
          HBB= QQINT(RAT,HB1,HB2)
         ENDIF

*   h -> tt

         IF (MH.LE.2.D0*MT)THEN
          HTT= 0.D0
         ELSE
          RATCOUP= 0.D0
          RMT= RUNM(MH,6)
          HT1= 3.D0*HFF(MH,(MT/MH)**2)
     .      * TQCDH((MT/MH)**2)
          HT2= 3.D0*HFF(MH,(RMT/MH)**2)
     .      * QCDH((RMT/MH)**2)
          IF(HT2.LT.0.D0) HT2=0.D0
          RAT= 2.D0*MT/MH
          HTT= QQINT(RAT,HT1,HT2)
         ENDIF

*   h -> WW

         IF(VFLAG.EQ.0)THEN
          DLD= 2.D0
          DLU= 2.D0
          XM1= 2.D0*MW-DLD
          XM2= 2.D0*MW+DLU
          IF(MH.LE.MW)THEN
           HWW= 0.D0
          ELSEIF(MH.LE.XM1)THEN
           CWW= 3.D0*GF**2*MW**4/16.D0/PI**3
           HWW= HV((MW/MH)**2)*CWW*MH
          ELSEIF(MH.LT.XM2)THEN
          CWW= 3.D0*GF**2*MW**4/16.D0/PI**3
           XX(1)= XM1-1.D0
           XX(2)= XM1
           XX(3)= XM2
           XX(4)= XM2+1.D0
           YY(1)= HV((MW/XX(1))**2)*CWW*XX(1)
           YY(2)= HV((MW/XX(2))**2)*CWW*XX(2)
           YY(3)= HVV(XX(3),(MW/XX(3))**2)
           YY(4)= HVV(XX(4),(MW/XX(4))**2)
           HWW= FINT(MH,XX,YY)
          ELSE
           HWW= HVV(MH,(MW/MH)**2)
          ENDIF
         ELSE
          CALL HTOVV(MW,2.08856D0,MH,HTWW)
          HWW = 3.D0/2.D0*GF*MW**4/DSQRT(2.D0)/PI/MH**3*HTWW
         ENDIF

*   h -> ZZ

         IF(VFLAG.EQ.0)THEN
          DLD= 2.D0
          DLU= 2.D0
          XM1= 2.D0*MZ-DLD
          XM2= 2.D0*MZ+DLU
          IF(MH.LE.MZ)THEN
           HZZ= 0.D0
          ELSEIF(MH.LE.XM1)THEN
           CZZ= 3.D0*GF**2*MZ**4/192.D0/PI**3
     .      * (7.D0-40.D0/3.D0*S2TW+160.D0/9.D0*S2TW**2)
           HZZ= HV((MZ/MH)**2)*CZZ*MH
          ELSEIF(MH.LT.XM2)THEN
           CZZ= 3.D0*GF**2*MZ**4/192.D0/PI**3
     .      * (7.D0-40.D0/3.D0*S2TW+160.D0/9.D0*S2TW**2)
           XX(1)= XM1-1.D0
           XX(2)= XM1
           XX(3)= XM2
           XX(4)= XM2+1.D0
           YY(1)= HV((MZ/XX(1))**2)*CZZ*XX(1)
           YY(2)= HV((MZ/XX(2))**2)*CZZ*XX(2)
           YY(3)= HVV(XX(3),(MZ/XX(3))**2)/2.D0
           YY(4)= HVV(XX(4),(MZ/XX(4))**2)/2.D0
           HZZ= FINT(MH,XX,YY)
          ELSE
           HZZ= HVV(MH,(MZ/MH)**2)/2.D0
          ENDIF
         ELSE
          CALL HTOVV(MZ,2.49581D0,MH,HTZZ)
          HZZ= 3.D0/4.D0*GF*MZ**4/DSQRT(2.D0)/PI/MH**3*HTZZ
         ENDIF
 
c   h -> gamma gamma

         XFAC= CG**2
         HGG= GF/(128.D0*PI*SQR2)*MH**3*(ALEM0/PI)**2*XFAC

*  h -> Z gamma

         IF(MH.LE.MZ)THEN
          HZG= 0.D0
         ELSE
          FT= -2.D0*(1.D0-8.D0/3.D0*S2TW)/DSQRT(S2TW*C2TW)
          FB= (-1.D0+4.D0/3.D0*S2TW)/DSQRT(S2TW*C2TW)
          CLT= 4.D0*(MT/MZ)**2*DCMPLX(1.D0,-EPS)
          CLB= 4.D0*(MBP/MZ)**2*DCMPLX(1.D0,-EPS)
          CLC= 4.D0*(MC/MZ)**2*DCMPLX(1.D0,-EPS)
          CLW= 4.D0*(MW/MZ)**2*DCMPLX(1.D0,-EPS)
          CXTZ= FT*(CI1(CTT,CLT) - CI2(CTT,CLT))
          CXBZ= FB*(CI1(CTB,CLB) - CI2(CTB,CLB))
          CXCZ= FT*(CI1(CTC,CLC) - CI2(CTC,CLC))
          CXWZ= -1.D0/DSQRT(T2TW)*(4.D0*(3.D0-T2TW)*CI2(CTW,CLW)
     .      + ((1.D0+2.D0/CTW)*T2TW - (5.D0+2.D0/CTW))*CI1(CTW,CLW))
          XFAC= CDABS(CXTZ+CXBZ+CXCZ+CXWZ)**2
          ACOUP= SQR2*GF*MZ**2*S2TW*C2TW/PI**2
          HZG= GF/(4.D0*PI*SQR2)*MH**3*(ALEM0/PI)*ACOUP/16.D0
     .      * XFAC*(1.D0-(MZ/MH)**2)**3
         ENDIF

*  Branching ratios

         WIDTHSM(I)=HJJ+HEE+HMM+HLL+HSS+HCC+HBB+HTT+HWW+HZZ+HGG+HZG
         BRJJSM(I)= HJJ/WIDTHSM(I)
         BREESM(I)= HEE/WIDTHSM(I)
         BRMMSM(I)= HMM/WIDTHSM(I)
         BRLLSM(I)= HLL/WIDTHSM(I)
         BRSSSM(I)= HSS/WIDTHSM(I)
         BRCCSM(I)= HCC/WIDTHSM(I)
         BRBBSM(I)= HBB/WIDTHSM(I)
         BRTTSM(I)= HTT/WIDTHSM(I)
         BRWWSM(I)= HWW/WIDTHSM(I)
         BRZZSM(I)= HZZ/WIDTHSM(I)
         BRGGSM(I)= HGG/WIDTHSM(I)
         BRZGSM(I)= HZG/WIDTHSM(I)

        ENDDO

        END
