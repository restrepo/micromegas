      SUBROUTINE RGESINVGM(PAR,IFAIL)

*   Subroutine to integrate the 2-loop RGEs for all soft terms
*   from MMESS down to Q2 = M_SUSY, through a CALL of the
*   subroutine ODEINTS that is part of the file integs.f
*
*   At MMESS:
*   MS^2 = MSMES (MAFLAG=-1,-3), MSINP (MAFLAG=-2,-4)
*   XIS = XISMES (MAFLAG=-2,-4), XISINP (MAFLAG=-1,-3)
*   XIF = XIFMES (MAFLAG=-3,-4), XIFINP (MAFLAG=-1,-2)
*
*   It uses COMMON/INPPAR and /MESEXT for the soft terms,
*   COMMON/MESCOUP for the gauge/Yukawa couplings at MMES.
*
*   The output (soft terms at Q2) is written into PAR(*) and
*   COMMON/SUSYMH and /SUSYEXT
*
***********************************************************************

      IMPLICIT NONE

      INTEGER IFAIL,OMGFLAG,MAFLAG,MOFLAG,GMFLAG,NN
      PARAMETER (NN=35)

      DOUBLE PRECISION PAR(*),EPS,X1,X2,Y(NN),PI,ALP1,ALP2,ALP3
      DOUBLE PRECISION Q2,COEF
      DOUBLE PRECISION G1MES,G2MES,G3MES,LMES,KMES,HTMES,HBMES,HLMES
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION MH1S,MH2S,MSS
      DOUBLE PRECISION MSUSYEFF,MMESS,N5
      DOUBLE PRECISION ALINP,XIFINP,XISINP,MSINP,MUPINP,MSPINP
      DOUBLE PRECISION XIFMES,XISMES,MSMES,MUPMES,MSPMES,M3HMES
      DOUBLE PRECISION LPPMES,LTTMES,LUMES,LDMES,LTMES,LBMES,LLMES,DHMES
      DOUBLE PRECISION XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY
      DOUBLE PRECISION M1INP,M2INP,M3INP,AKINP,ATINP,ABINP,
     .      ATAUINP,AMUINP,MH1INP,MH2INP,MQ3INP,MU3INP,MD3INP,
     .      MQINP,MUINP,MDINP,ML3INP,ME3INP,MLINP,MEINP
      DOUBLE PRECISION X,F1,F2,SP2,MSREF,D,DMIN

      COMMON/FLAGS/OMGFLAG,MAFLAG,MOFLAG
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/RENSCALE/Q2
      COMMON/INPPAR/ALINP,XIFINP,XISINP,MSINP,MUPINP,MSPINP
      COMMON/SOFTINP/M1INP,M2INP,M3INP,AKINP,ATINP,ABINP,
     .      ATAUINP,AMUINP,MH1INP,MH2INP,MQ3INP,MU3INP,MD3INP,
     .      MQINP,MUINP,MDINP,ML3INP,ME3INP,MLINP,MEINP
      COMMON/MESCAL/MSUSYEFF,MMESS,N5
      COMMON/MESCOUP/G1MES,G2MES,G3MES,LMES,KMES,HTMES,HBMES,HLMES      
      COMMON/MESEXT/XIFMES,XISMES,MSMES,MUPMES,MSPMES,M3HMES
      COMMON/MESGUT/LPPMES,LTTMES,LUMES,LDMES,LTMES,LBMES,LLMES,DHMES
      COMMON/SUSYEXT/XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY
      COMMON/SUSYMH/MH1S,MH2S,MSS
      COMMON/GMSCEN/MSREF,D,DMIN,GMFLAG
      
      EXTERNAL DERIVSS,RKQSS

      EPS=1d-8
      PI=4d0*DATAN(1d0)
      COEF=1d0/(16d0*PI**2)

* Definition of the couplings squared Y(I) at MMESS

      Y(1)=G1MES
      Y(2)=G2MES
      Y(3)=G3MES
      Y(4)=LMES
* NOTE: KMES = KAPPA, NOT KAPPA**2
      Y(5)=KMES
      Y(6)=HTMES
      Y(7)=HBMES
      Y(8)=HLMES

* Corrections higher order in X=MSUSYEFF/MMESS:

      X=MSUSYEFF/MMESS
      IF(X.GE.1.D-7) THEN
        IF(X.LT.1d0) THEN
          F1=((1d0+X)*DLOG(1d0+X)+(1d0-X)*DLOG(1d0-X))/X**2
          F2=(1d0+X)/X**2*(DLOG(1d0+X)-2d0*SP2(X/(1d0+X))
     .       +SP2(2d0*X/(1d0+X))/2d0)
     .       +(1d0-X)/X**2*(DLOG(1d0-X)-2d0*SP2(-X/(1d0-X))
     .       +SP2(-2d0*X/(1d0-X))/2d0)
        ELSE
          F1=0d0
          F2=0d0
        ENDIF
      ELSE
        F1=1d0
        F2=1d0
      ENDIF

* (SP2 = Li_2  is defined in the subroutine bsg.f)
      
* Input values for the soft terms at MMESS:

      ALP1=G1MES/(4d0*PI)
      ALP2=G2MES/(4d0*PI)
      ALP3=G3MES/(4d0*PI)
      
      Y(9)=5d0/3d0*ALP1*N5*MSUSYEFF*F1/(4d0*PI)
      Y(10)=ALP2*N5*MSUSYEFF*F1/(4d0*PI)
      Y(11)=ALP3*N5*MSUSYEFF*F1/(4d0*PI)
      IF(GMFLAG.EQ.1)THEN
       ALINP=-COEF*MSUSYEFF*(DHMES
     .       +2d0*LPPMES+3d0*LTTMES+3d0*LUMES+3d0*LDMES)
      ENDIF
      Y(12)=ALINP
      IF(PAR(2).NE.0d0)THEN
       IF(GMFLAG.EQ.1)THEN
        Y(13)=-3d0*COEF*MSUSYEFF*(DHMES
     .        +2d0*LPPMES+3d0*LTTMES+2d0*LUMES+2d0*LDMES)
       ELSE
        Y(13)=3d0*ALINP
       ENDIF
      ELSE
       Y(13)=0d0
      ENDIF
      Y(14)=-COEF*MSUSYEFF*(LUMES+3d0*LTMES+LBMES)
      Y(15)=-COEF*MSUSYEFF*(LDMES+LTMES+3d0*LBMES)
      Y(16)=-COEF*MSUSYEFF*(LDMES+3d0*LLMES)
      Y(17)=COEF*(5d0/6d0*ALP1**2+3d0/2d0*ALP2**2)*N5*F2*MSUSYEFF**2
     .     +(COEF*MSUSYEFF)**2*
     .      (-LMES*(DHMES+2d0*LPPMES+3d0*LTTMES+3d0*LDMES)
     .       +LUMES*(4d0*LUMES+2d0*LDMES+2d0*LMES+2d0*KMES
     .       +3d0*LBMES+LLMES-G1MES-3d0*G2MES)
     .       -3d0*HTMES*(3d0*LTMES+LBMES)
     .       +2d0*DSQRT(LMES*LPPMES*LUMES*LDMES)
     .       +6d0*DSQRT(LMES*LUMES*LBMES*HBMES)
     .       +2d0*DSQRT(LMES*LUMES*LLMES*HLMES))
      Y(18)=COEF*(5d0/6d0*ALP1**2+3d0/2d0*ALP2**2)*N5*F2*MSUSYEFF**2
     .     +(COEF*MSUSYEFF)**2*
     .      (-LMES*(DHMES+2d0*LPPMES+3d0*LTTMES+3d0*LUMES)
     .       +LDMES*(4d0*LDMES+2d0*LUMES+2d0*LMES+2d0*KMES
     .       +3d0*LTMES-G1MES-3d0*G2MES)
     .       -3d0*HBMES*(3d0*LBMES+LTMES)-3d0*HLMES*LLMES
     .       +2d0*DSQRT(LMES*LPPMES*LUMES*LDMES)
     .       +6d0*DSQRT(LMES*LDMES*LTMES*HTMES))
      IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3)THEN
       Y(19)=MSMES
      ELSE
       Y(19)=MSINP
      ENDIF
      Y(20)=COEF*(5d0/54d0*ALP1**2+3d0/2d0*ALP2**2
     .     +8d0/3d0*ALP3**2)*F2*N5*MSUSYEFF**2
     .     +(COEF*MSUSYEFF)**2*
     .      (-(13d0/9d0*G1MES+3d0*G2MES+16d0/3d0*G3MES)*LTMES
     .       -(7d0/9d0*G1MES+3d0*G2MES+16d0/3d0*G3MES)*LBMES
     .       +6d0*LTMES*(LTMES+HTMES)+LDMES*(LTMES-HBMES)
     .       +6d0*LBMES*(LBMES+HBMES)+LUMES*(LBMES-HTMES)
     .       +LBMES*LLMES+2d0*DSQRT(LMES*LDMES*LTMES*HTMES)
     .       +2d0*DSQRT(LMES*LUMES*LBMES*HBMES)
     .       +2d0*DSQRT(HBMES*HLMES*LBMES*LLMES))
      Y(21)=COEF*(40d0/27d0*ALP1**2
     .     +8d0/3d0*ALP3**2)*F2*N5*MSUSYEFF**2
     .     +(COEF*MSUSYEFF)**2*
     .      (-(26d0/9d0*G1MES+6d0*G2MES+32d0/3d0*G3MES)*LTMES
     .       +12d0*LTMES*(LTMES+HTMES)-2d0*HTMES*(LUMES+LBMES)
     .       +2d0*LTMES*(HBMES+LDMES)
     .       +4d0*DSQRT(LMES*LDMES*LTMES*HTMES))
      Y(22)=COEF*(10d0/27d0*ALP1**2
     .     +8d0/3d0*ALP3**2)*F2*N5*MSUSYEFF**2
     .     +(COEF*MSUSYEFF)**2*
     .      (-(14d0/9d0*G1MES+6d0*G2MES+32d0/3d0*G3MES)*LTMES
     .       +12d0*LBMES*(LBMES+HBMES)-2d0*HBMES*(LDMES+LTMES)
     .       +2d0*LBMES*(HTMES+LUMES+LLMES)
     .       +4d0*DSQRT(LMES*LUMES*LBMES*HBMES)
     .       +4d0*DSQRT(HBMES*HLMES*LBMES*LLMES))
      Y(23)=COEF*(5d0/54d0*ALP1**2+3d0/2d0*ALP2**2
     .     +8d0/3d0*ALP3**2)*F2*N5*MSUSYEFF**2
      Y(24)=COEF*(40d0/27d0*ALP1**2
     .     +8d0/3d0*ALP3**2)*F2*N5*MSUSYEFF**2
      Y(25)=COEF*(10d0/27d0*ALP1**2
     .     +8d0/3d0*ALP3**2)*F2*N5*MSUSYEFF**2
      Y(26)=COEF*(5d0/6d0*ALP1**2
     .     +3d0/2d0*ALP2**2)*F2*N5*MSUSYEFF**2
     .     +(COEF*MSUSYEFF)**2*
     .      (-(3d0*G1MES+3d0*G2MES)*LLMES
     .       +LLMES*(4d0*LLMES+3d0*LBMES+2d0*HLMES+LUMES)
     .       -LDMES*HLMES+2d0*DSQRT(LMES*LUMES*LLMES*HLMES)
     .       +6d0*DSQRT(HBMES*HLMES*LBMES*LLMES))
      Y(27)=COEF*(10d0/3d0*ALP1**2)*F2*N5*MSUSYEFF**2
     .     +(COEF*MSUSYEFF)**2*
     .      (-(6d0*G1MES+6d0*G2MES)*LLMES
     .       +LLMES*(8d0*LLMES+6d0*LBMES+4d0*HLMES+2d0*LUMES)
     .       -2d0*LDMES*HLMES+4d0*DSQRT(LMES*LUMES*LLMES*HLMES)
     .       +12d0*DSQRT(HBMES*HLMES*LBMES*LLMES))
      Y(28)=COEF*(5d0/6d0*ALP1**2
     .     +3d0/2d0*ALP2**2)*F2*N5*MSUSYEFF**2
      Y(29)=COEF*(10d0/3d0*ALP1**2)*F2*N5*MSUSYEFF**2
      IF(MAFLAG.EQ.-3 .OR. MAFLAG.EQ.-4)THEN
       Y(30)=XIFMES
      ELSE
       Y(30)=XIFINP
      ENDIF
      IF(MAFLAG.EQ.-2 .OR. MAFLAG.EQ.-4)THEN
       Y(31)=XISMES
      ELSE
       Y(31)=XISINP
      ENDIF
      Y(32)=MUPINP
      Y(33)=MSPINP
      Y(34)=0d0
      Y(35)=0d0

* Save the required input values for soft masses and couplings in
* COMMON/SOFTINP

      M1INP=Y(9)
      M2INP=Y(10)
      M3INP=Y(11)
      AKINP=Y(13)
      ATINP=Y(14)
      ABINP=Y(15)
      ATAUINP=Y(16)
      MH1INP=Y(17)
      MH2INP=Y(18)
      MQ3INP=Y(20)
      MU3INP=Y(21)
      MD3INP=Y(22)
      MQINP=Y(23)
      MUINP=Y(24)
      MDINP=Y(25)
      ML3INP=Y(26)
      ME3INP=Y(27)
      MLINP=Y(28)
      MEINP=Y(29)
      AMUINP=Y(35)

      X1=COEF*DLOG(MMESS**2/Q2)
      X2=0d0

      !WRITE(0,*)"CALL RGESINVGM"
      !WRITE(0,*)""
      !WRITE(0,*)"MMESS =",MMESS
      !WRITE(0,*)"G1MES =",Y(1)
      !WRITE(0,*)"G2MES =",Y(2)
      !WRITE(0,*)"G3MES =",Y(3)
      !WRITE(0,*)"LMES =",DSQRT(Y(4))
      !WRITE(0,*)"KMES =",Y(5)
      !WRITE(0,*)"HTMES =",Y(6)
      !WRITE(0,*)"HBMES =",Y(7)
      !WRITE(0,*)"HLMES =",Y(8)
      !WRITE(0,*)"M1MES =",Y(9)
      !WRITE(0,*)"M2MES =",Y(10)
      !WRITE(0,*)"M3MES =",Y(11)
      !WRITE(0,*)"ALMES =",Y(12)
      !WRITE(0,*)"AKMES =",Y(13)
      !WRITE(0,*)"ATOPMES =",Y(14)
      !WRITE(0,*)"ABOTMES =",Y(15)
      !WRITE(0,*)"ATAUMES =",Y(16)
      !WRITE(0,*)"AMUMES =",Y(35)
      !WRITE(0,*)"MH1MES =",Y(17)
      !WRITE(0,*)"MH2MES =",Y(18)
      !WRITE(0,*)"MSMES =",Y(19)
      !WRITE(0,*)"MQ3MES =",Y(20)
      !WRITE(0,*)"MU3MES =",Y(21)
      !WRITE(0,*)"MD3MES =",Y(22)
      !WRITE(0,*)"MQMES =",Y(23)
      !WRITE(0,*)"MUMES =",Y(24)
      !WRITE(0,*)"MDMES =",Y(25)
      !WRITE(0,*)"ML3MES =",Y(26)
      !WRITE(0,*)"ME3MES =",Y(27)
      !WRITE(0,*)"MLMES =",Y(28)
      !WRITE(0,*)"MEMES =",Y(29)
      !WRITE(0,*)"XIFMES =",Y(30)
      !WRITE(0,*)"XISMES =",Y(31)
      !WRITE(0,*)"MUPMES =",Y(32)
      !WRITE(0,*)"MSPMES =",Y(33)
      !WRITE(0,*)"M3HMES =",Y(34)
      !WRITE(0,*)""

      CALL ODEINTS(Y,NN,X1,X2,EPS,DERIVSS,RKQSS,IFAIL)      

      !WRITE(0,*)"MSUSY =",DSQRT(Q2)
      !WRITE(0,*)"G1 =",Y(1)
      !WRITE(0,*)"G2 =",Y(2)
      !WRITE(0,*)"G3 =",Y(3)
      !WRITE(0,*)"L =",DSQRT(Y(4))
      !WRITE(0,*)"K =",Y(5)
      !WRITE(0,*)"HT =",Y(6)
      !WRITE(0,*)"HB =",Y(7)
      !WRITE(0,*)"HL =",Y(8)
      !WRITE(0,*)"M1 =",Y(9)
      !WRITE(0,*)"M2 =",Y(10)
      !WRITE(0,*)"M3 =",Y(11)
      !WRITE(0,*)"AL =",Y(12)
      !WRITE(0,*)"AK =",Y(13)
      !WRITE(0,*)"ATOP =",Y(14)
      !WRITE(0,*)"ABOT =",Y(15)
      !WRITE(0,*)"ATAU =",Y(16)
      !WRITE(0,*)"AMUON =",Y(35)
      !WRITE(0,*)"MH1 =",Y(17)
      !WRITE(0,*)"MH2 =",Y(18)
      !WRITE(0,*)"MS =",Y(19)
      !WRITE(0,*)"MQ3 =",Y(20)
      !WRITE(0,*)"MU3 =",Y(21)
      !WRITE(0,*)"MD3 =",Y(22)
      !WRITE(0,*)"MQ =",Y(23)
      !WRITE(0,*)"MU =",Y(24)
      !WRITE(0,*)"MD =",Y(25)
      !WRITE(0,*)"ML3 =",Y(26)
      !WRITE(0,*)"ME3 =",Y(27)
      !WRITE(0,*)"ML =",Y(28)
      !WRITE(0,*)"ME =",Y(29)
      !WRITE(0,*)"XIF =",Y(30)
      !WRITE(0,*)"XIS =",Y(31)
      !WRITE(0,*)"MUP =",Y(32)
      !WRITE(0,*)"MSP =",Y(33)
      !WRITE(0,*)"M3H =",Y(34)
      !WRITE(0,*)""

      IF(IFAIL.NE.0)THEN
       IFAIL=13
       RETURN
      ENDIF

* Y(5) = KAPPA, NOT KAPPA**2

      IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-2)PAR(2)=Y(5)

* SOFT TERMS AT THE SUSY SCALE

      PAR(5)=Y(12)
      PAR(6)=Y(13)
      PAR(7)=Y(20)
      PAR(8)=Y(21)
      PAR(9)=Y(22)
      PAR(10)=Y(26)
      PAR(11)=Y(27)
      PAR(12)=Y(14)
      PAR(13)=Y(15)
      PAR(14)=Y(16)
      PAR(15)=Y(23)
      PAR(16)=Y(24)
      PAR(17)=Y(25)
      PAR(18)=Y(28)
      PAR(19)=Y(29)
      PAR(20)=Y(9)
      PAR(21)=Y(10)
      PAR(22)=Y(11)
      PAR(25)=Y(35)

* MH1S, MH2S AND MSS at Q2 are stored in COMMON/SUSYMH:
      
      MH1S=Y(17)
      MH2S=Y(18)
      MSS=Y(19)
      
* EXT parameters at Q2, stored in COMMON/SUSYEXT:

      XIFSUSY=Y(30)
      XISSUSY=Y(31)
      MUPSUSY=Y(32)
      MSPSUSY=Y(33)
      M3HSUSY=Y(34)
      
      END
