      SUBROUTINE MCMCSTEPGM(PAR,PROB,NPROB,IFAIL) 

      IMPLICIT NONE 
 
      INTEGER IFAIL,I,OMGFLAG,MAFLAG,MOFLAG,GMFLAG,NTOT,IDUM,NPROB
      INTEGER MCFLAG,NSUSY,NGUT,NMES
      PARAMETER (NSUSY=14,NGUT=21,NMES=22)
 
      DOUBLE PRECISION PAR(*),PROB(*),P,RAN2,XDUM
      DOUBLE PRECISION FTSUSY(NSUSY+2),FTGUT(NGUT+2),FTMES(NMES+2)
      DOUBLE PRECISION SMASS(3),PMASS(2),CMASS,SCOMP(3,3),PCOMP(2,2)
      DOUBLE PRECISION MGL,MCHA(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
      DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
      DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU,Q2
      DOUBLE PRECISION RMST1,RMST2,S2T,RMSB1,RMSB2,S2B,XT,XB
      DOUBLE PRECISION MSUSYEFF,MMESS,N5
      DOUBLE PRECISION ALINP,XIFINP,XISINP,MSINP,MUPINP,MSPINP
      DOUBLE PRECISION MHmin,MHmax,chi2max,chi2gam,chi2bb,chi2zz
      DOUBLE PRECISION MSUSYEFFCEN,MSUSYEFFDEV,MMESSCEN,MMESSDEV
      DOUBLE PRECISION TBCEN,TBDEV,LCEN,LDEV,KCEN,KDEV,ALCEN,ALDEV
      DOUBLE PRECISION XIFCEN,XIFDEV,XISCEN,XISDEV,MUPCEN,MUPDEV
      DOUBLE PRECISION MSPCEN,MSPDEV,MSCEN,MSDEV,DHCEN,DHDEV
      DOUBLE PRECISION LPPCEN,LPPDEV,LTTCEN,LTTDEV,LUCEN,LUDEV
      DOUBLE PRECISION LDCEN,LDDEV, LTCEN,LTDEV,LBCEN,LBDEV
      DOUBLE PRECISION LLCEN,LLDEV,XCEN,XDEV,X
      DOUBLE PRECISION MSUSYEFFMIN,MMESSMIN,TBMIN,LMIN,KMIN,ALMIN
      DOUBLE PRECISION XIFMIN,XISMIN,MUPMIN,MSPMIN,MSMIN,DHMIN
      DOUBLE PRECISION LPPMIN,LTTMIN,LUMIN,LDMIN,LTMIN,LBMIN,LLMIN
      DOUBLE PRECISION LPPMES,LTTMES,LUMES,LDMES,LTMES,LBMES,LLMES,DHMES
      DOUBLE PRECISION MSREF,D,DMIN

      COMMON/INPPAR/ALINP,XIFINP,XISINP,MSINP,MUPINP,MSPINP
      COMMON/MESGUT/LPPMES,LTTMES,LUMES,LDMES,LTMES,LBMES,LLMES,DHMES
      COMMON/MESCAL/MSUSYEFF,MMESS,N5
      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,NEU
      COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     . MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     . CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
      COMMON/RADCOR/RMST1,RMST2,S2T,RMSB1,RMSB2,S2B,XT,XB
      COMMON/HIGGSFIT/MHmin,MHmax,chi2max,chi2gam,chi2bb,chi2zz
      COMMON/MCMCPAR/MSUSYEFFCEN,MSUSYEFFDEV,MMESSCEN,MMESSDEV,
     . TBCEN,TBDEV,LCEN,LDEV,KCEN,KDEV,ALCEN,ALDEV,XIFCEN,
     . XIFDEV,XISCEN,XISDEV,MUPCEN,MUPDEV,MSPCEN,MSPDEV,
     . MSCEN,MSDEV,DHCEN,DHDEV,LPPCEN,LPPDEV,LTTCEN,LTTDEV,
     . LUCEN,LUDEV,LDCEN,LDDEV, LTCEN,LTDEV,LBCEN,LBDEV,
     . LLCEN,LLDEV,XCEN,XDEV,X,
     . MSUSYEFFMIN,MMESSMIN,TBMIN,LMIN,KMIN,ALMIN,
     . XIFMIN,XISMIN,MUPMIN,MSPMIN,MSMIN,DHMIN,
     . LPPMIN,LTTMIN,LUMIN,LDMIN,LTMIN,LBMIN,LLMIN
      COMMON/FINETUN/FTSUSY,FTGUT,FTMES
      COMMON/FLAGS/OMGFLAG,MAFLAG,MOFLAG
      COMMON/STEPS/NTOT,IDUM
      COMMON/GMSCEN/MSREF,D,DMIN,GMFLAG

      X=0d0

      IF(IFAIL.EQ.9.OR.(IFAIL.GE.11 .AND. IFAIL.LE.17))THEN
       X=1d40
      ENDIF

      IF(IFAIL.EQ.1.OR.IFAIL.EQ.3.OR.IFAIL.EQ.5.OR.IFAIL.EQ.7)THEN
       X=X+(2d0-SMASS(1))*1d20
      ENDIF
      IF(IFAIL.EQ.2.OR.IFAIL.EQ.3.OR.IFAIL.EQ.6.OR.IFAIL.EQ.7)THEN
       X=X+(2d0-PMASS(1))*1d20
      ENDIF
      IF(IFAIL.EQ.4.OR.IFAIL.EQ.5.OR.IFAIL.EQ.6.OR.IFAIL.EQ.7)THEN
       X=X+(2d0-CMASS)*1d20
      ENDIF

      IF(IFAIL.EQ.8)THEN
       X=MIN(X,RMST1)
       X=MIN(X,RMSB1)
       X=MIN(X,MST1)
       X=MIN(X,MST2)
       X=MIN(X,MSB1)
       X=MIN(X,MSB2)
       X=MIN(X,MUL)
       X=MIN(X,MUR)
       X=MIN(X,MDL)
       X=MIN(X,MDR)
       X=MIN(X,MSL1)
       X=MIN(X,MSNT)
       X=MIN(X,MSMU1)
       X=MIN(X,MSMUNT)
       X=MIN(X,MLR)
       X=MIN(X,MLL)
       X=MIN(X,MNL)
       X=(1d0-X)*1d20
      ENDIF

      IF(IFAIL.EQ.10 .OR. IFAIL.EQ.19)THEN
       DO I=1,NPROB
        X=X+DABS(PROB(I))
       ENDDO
       IF(X.NE.0d0)X=(1d0+X)*1d15
      ENDIF

      IF(IFAIL.EQ.18 .OR. IFAIL.EQ.19)THEN
       X=X+(1d0+D-DMIN)*1d10
      ENDIF

      IF(MCFLAG.EQ.1)THEN
       IF(FTMES(NMES+1).GT.1d3)THEN
        X=X+(FTMES(NMES+1)/1d3)*1d5
       ENDIF
       X=X+1d0
      ELSEIF(MCFLAG.EQ.2)THEN
       X=X+FTMES(NMES+1)
      ELSE
       X=X+1d0
      ENDIF

c      P=1d0/(1d0+DEXP(DLOG10(X/XCEN)/(.28d0*XDEV)))
      P=1d0/(1d0+DEXP((X-XCEN)/(.28d0*XDEV*MIN(X,XCEN))))
      XDUM=RAN2(IDUM)
!      WRITE(0,*)"X",XCEN,X
!      WRITE(0,*)"P",P,XDUM
      IF(P.GE.XDUM)THEN
!       WRITE(0,*)"OK",IFAIL
       XCEN=X
       MSUSYEFFCEN=MSUSYEFF
       MMESSCEN=MMESS
       TBCEN=PAR(3)
       LCEN=PAR(1)
       ALCEN=ALINP
       IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-2)THEN
        XIFCEN=XIFINP
       ELSE
        KCEN=PAR(2)
       ENDIF
       IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3)THEN
        XISCEN=XISINP
       ELSE
        MSCEN=MSINP
       ENDIF
       MUPCEN=MUPINP
       MSPCEN=MSPINP
       DHCEN=DHMES
       LPPCEN=LPPMES
       LTTCEN=LTTMES
       LUCEN=LUMES
       LDCEN=LDMES
       LTCEN=LTMES
       LBCEN=LBMES
       LLCEN=LLMES
      ELSE
!       WRITE(0,*)"NO",IFAIL
      ENDIF 

      END
