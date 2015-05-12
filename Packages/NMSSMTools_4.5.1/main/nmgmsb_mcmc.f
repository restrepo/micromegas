      PROGRAM MAIN

*  Program to compute the NMSSM Higgs masses, couplings and
*  Branching ratios, with experimental and theoretical constraints
*
*  Input data corresponds to GMSB model. In addition:
*
*      tan(beta) at the scale MZ, lambda at the scale Q2
*
*      The input file contains the starting point and number of MCMC steps
*
*  On output:
*
*      PAR(1) = lambda
*      PAR(2) = kappa
*      PAR(3) = tan(beta)
*      PAR(4) = mu (effective mu term = lambda*s)
*      PAR(5) = Alambda
*      PAR(6) = Akappa
*      PAR(7) = mQ3**2
*      PAR(8) = mU3**2
*      PAR(9) = mD3**2
*      PAR(10) = mL3**2
*      PAR(11) = mE3**2
*      PAR(12) = AU3
*      PAR(13) = AD3
*      PAR(14) = AE3
*      PAR(15) = mQ2**2
*      PAR(16) = mU2**2
*      PAR(17) = mD2**2
*      PAR(18) = mL2**2
*      PAR(19) = mE2**2
*      PAR(20) = M1
*      PAR(21) = M2
*      PAR(22) = M3
*      PAR(23) = MA (diagonal doublet CP-odd mass matrix element)
*      PAR(24) = MP (diagonal singlet CP-odd mass matrix element)
*      PAR(25) = AE2
*
*      All these parameters are assumed to be defined in DRbar at the scale
*      Q2 which is either user defined or computed as (2*mQ2+mU2+mD2)/4
*
*      SMASS(1-3): CP-even masses (ordered)
*
*      SCOMP(1-3,1-3): Mixing angles: if HB(I) are the bare states,
*        HB(I) = Re(H1), Re(H2), Re(S), and HM(I) are the mass eigenstates, 
*        the convention is HB(I) = SUM_(J=1,3) SCOMP(J,I)*HM(J)
*        which is equivalent to HM(I) = SUM_(J=1,3) SCOMP(I,J)*HB(J)
*
*      PMASS(1-2): CP-odd masses (ordered)
*
*      PCOMP(1-2,1-2): Mixing angles: if AB(I) are the bare states,
*        AB(I) = Im(H1), Im(H2), Im(S), and AM(I) are the mass eigenstates, 
*        the convention is 
*        AM(I) = PCOMP(I,1)*(COSBETA*AB(1)+SINBETA*AB(2))
*              + PCOMP(I,2)*AB(3)
*
*      CMASS: Charged Higgs mass
*
*      CU,CD,CV,CJ,CG(i) Reduced couplings of h1,h2,h3 (i=1,2,3) or
*                        a1,a2 (i=4,5) to up type fermions, down type
*                        fermions, gauge bosons, gluons and photons
*                        Note: CV(4)=CV(5)=0
*
*      WIDTH(i) Total decay width of h1,h2,h3,a1,a2 (i=1..5)
*               with the following branching ratios:
*      BRJJ(i) h1,h2,h3,a1,a2 -> gluon gluon
*      BRMM(i)        "       -> mu mu
*      BRLL(i)        "       -> tau tau
*      BRSS(i)        "       -> ss
*      BRCC(i)        "       -> cc
*      BRBB(i)        "       -> bb
*      BRTT(i)        "       -> tt
*      BRWW(i)        "       -> WW (BRWW(4)=BRWW(5)=0)
*      BRZZ(i)        "       -> ZZ (BRZZ(4)=BRZZ(5)=0)
*      BRGG(i)        "       -> gamma gamma
*      BRZG(i)        "       -> Z gamma
*      BRHIGGS(i)   (i=1..5)  -> other Higgses, including:
*        BRHAA(i,j)   hi      -> a1a1, a1a2, a2a2 (i=1..3, j=1..3)
*        BRHCHC(i)    hi      -> h+h- (i=1..3)
*        BRHAZ(i,j)   hi      -> Zaj  (i=1..3)
*        BRHCW(i)  h1,h2,h3   -> h+W- (i=1..3), a1,a2 -> h+W- (i=4,5)
*        BRHHH(i)     h2      -> h1h1, h3-> h1h1, h1h2, h2h2 (i=1..4)
*        BRAHA(i)     a2      -> a1hi (i=1..3)
*        BRAHZ(i,j)   ai      -> Zhj  (i=1,2, j=1..3)
*      BRSUSY(i)    (i=1..5)  -> susy particles, including:
*        BRNEU(i,j,k)         -> neutralinos j,k (i=1..5, j,k=1..5)
*        BRCHA(i,j)           -> charginos 11, 12, 22 (i=1..5, j=1..3)
*        BRHSQ(i,j)   hi      -> uLuL, uRuR, dLdL, dRdR, t1t1, t2t2,
*                                t1t2, b1b1, b2b2, b1b2 (i=1..3, j=1..10)
*        BRASQ(i,j)   ai      -> t1t2, b1b2 (i=1,2, j=1,2)
*        BRHSL(i,j)   hi      -> lLlL, lRlR, nLnL, l1l1, l2l2, l1l2,
*                                ntnt (i=1..3, j=1..7)
*        BRASL(i)     ai      -> l1l2 (i=1,2)
*
*      HCWIDTH  Total decay width of the charged Higgs
*               with the following branching ratios:
*      HCBRM         h+ -> mu nu_mu
*      HCBRL         "  -> tau nu_tau
*      HCBRSU        "  -> s u
*      HCBRBU        "  -> b u
*      HCBRSC        "  -> s c
*      HCBRBC        "  -> b c
*      HCBRBT        "  -> b t
*      HCBRWHT       "  -> neutral Higgs W+, including:
*        HCBRWH(i)   "  -> H1W+, H2W+, h3W+, a1W+, a2W+ (i=1..5)
*      HCBRSUSY      "  -> susy particles,including
*        HCBRNC(i,j) "  -> neutralino i chargino j (i=1..5, j=1,2)
*        HCBRSQ(i)   "  -> uLdL, t1b1, t1b2, t2b1, t2b2 (i=1..5)
*        HCBRSL(i)   "  -> lLnL, t1nt, t2nt (i=1..3)
*
*      MNEU(i)   Mass of neutralino chi_i (i=1,5, ordered in mass)
*      NEU(i,j)  chi_i components of bino, wino, higgsino u&d, singlino 
*                (i,j=1..5)
*
*      MCHA(i)       Chargino masses
*      U(i,j),V(i,j) Chargino mixing matrices
*
*  ERRORS: IFAIL = 0..19
*
*  IFAIL = 0         OK
*          1,3,5,7   m_h1^2 < 0
*          2,3,6,7   m_a1^2 < 0
*          4,5,6,7   m_h+^2 < 0
*          8         m_sfermion^2 < 0
*          9         l, tan(beta) or mu = 0
*          10        Violation of phenomenological constraint(s)
*          11,12,13  Problem in integration of RGEs
*          14,15,16  Convergence problem
*          17        No electroweak symmetry breaking
*          18        MSMES=/=MSREF (GMSB scenario at MMESS)
*          19        IFAIL = 10 & 18
*
*  Phenomenological constraints:
*
*      PROB(I)  = 0, I = 1..52: OK
*            
*      PROB(1) =/= 0   chargino too light
*      PROB(2) =/= 0   excluded by Z -> neutralinos
*      PROB(3) =/= 0   charged Higgs too light
*      PROB(4) =/= 0   excluded by ee -> hZ 
*      PROB(5) =/= 0   excluded by ee -> hZ, h -> bb
*      PROB(6) =/= 0   excluded by ee -> hZ, h -> tautau
*      PROB(7) =/= 0   excluded by ee -> hZ, h -> invisible 
*      PROB(8) =/= 0   excluded by ee -> hZ, h -> 2jets
*      PROB(9) =/= 0   excluded by ee -> hZ, h -> 2photons
*      PROB(10) =/= 0  excluded by ee -> hZ, h -> AA -> 4bs
*      PROB(11) =/= 0  excluded by ee -> hZ, h -> AA -> 4taus
*      PROB(12) =/= 0  excluded by ee -> hZ, h -> AA -> 2bs 2taus
*      PROB(13) =/= 0  excluded by Z -> hA (Z width)
*      PROB(14) =/= 0  excluded by ee -> hA -> 4bs
*      PROB(15) =/= 0  excluded by ee -> hA -> 4taus
*      PROB(16) =/= 0  excluded by ee -> hA -> 2bs 2taus
*      PROB(17) =/= 0  excluded by ee -> hA -> AAA -> 6bs
*      PROB(18) =/= 0  excluded by ee -> hA -> AAA -> 6taus
*      PROB(19) =/= 0  excluded by ee -> Zh -> ZAA -> Z + light pairs
*      PROB(20) =/= 0  excluded by stop -> b l sneutrino
*      PROB(21) =/= 0  excluded by stop -> neutralino c
*      PROB(22) =/= 0  excluded by sbottom -> neutralino b
*      PROB(23) =/= 0  squark/gluino too light
*      PROB(24) =/= 0  selectron/smuon too light
*      PROB(25) =/= 0  stau too light
*      PROB(27) =/= 0  Landau Pole in l, k, ht, hb below MGUT
*      PROB(28) =/= 0  unphysical global minimum
*      PROB(29) =/= 0  Higgs soft masses >> Msusy
*      PROB(32) =/= 0  b->s gamma more than 2 sigma away
*      PROB(33) =/= 0  Delta M_s more than 2 sigma away
*      PROB(34) =/= 0  Delta M_d more than 2 sigma away
*      PROB(35) =/= 0  B_s->mu+mu- more than 2 sigma away
*      PROB(36) =/= 0  B+-> tau+nu_tau more than 2 sigma away
*      PROB(37) =/= 0  (g-2)_muon more than 2 sigma away
*      PROB(38) =/= 0  excluded by Upsilon(1S) -> A gamma
*      PROB(39) =/= 0  excluded by eta_b(1S) mass measurement
*      PROB(40) =/= 0  BR(B-->X_s mu+ mu-) more than 2 sigma away
*      PROB(41) =/= 0  excluded by ee -> hZ, h -> AA -> 4taus (ALEPH analysis)
*      PROB(42) =/= 0  excluded by top -> b H+, H+ -> c s (CDF, D0)
*      PROB(43) =/= 0  excluded by top -> b H+, H+ -> tau nu_tau (D0)
*      PROB(44) =/= 0  excluded by top -> b H+, H+ -> W+ A1, A1 -> 2taus (CDF)
*      PROB(45) =/= 0  excluded by t -> bH+ (LHC)
*      PROB(46) =/= 0  No Higgs in the MHmin-MHmax GeV range
*      PROB(47) =/= 0  chi2gam > chi2max
*      PROB(48) =/= 0  chi2bb > chi2max
*      PROB(49) =/= 0  chi2zz > chi2max
*      PROB(50) =/= 0  excluded by sparticle searches at the LHC
*      PROB(51) =/= 0: excluded by H/A->tautau
*      PROB(52) =/= 0: Excluded H_125->AA->4mu (CMS)
*
************************************************************************

      IMPLICIT NONE

      INTEGER NFL,NPROB,NPAR
      PARAMETER (NFL=19,NPROB=52,NPAR=25)
      INTEGER NFAIL(NFL),IFAIL,I,TOT,ITOT,NTOT,IDUM,GMFLAG
      INTEGER ITER,ITERMU,Q2FIX,OMGFLAG,MAFLAG,MOFLAG,GMUFLAG,HFLAG
      INTEGER TOTMIN,TOTMAX,NMAX,IP

      DOUBLE PRECISION PAR(NPAR),PROB(NPROB),CHECK,MESTEST,DELMB
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
      DOUBLE PRECISION MSUSYEFFN,MSUSYEFFNN,MMESSN,MMESSNN,TBN,TBNN
      DOUBLE PRECISION LN,LNN,KN,KNN,ALN,ALNN,XIFN,XIFNN,XISN,XISNN
      DOUBLE PRECISION MUPN,MUPNN,MSPN,MSPNN,MSN,MSNN,DHN,DHNN
      DOUBLE PRECISION LPPN,LPPNN,LTTN,LTTNN,LUN,LUNN,LDN,LDNN
      DOUBLE PRECISION LTN,LTNN,LBN,LBNN,LLN,LLNN,MUN,MUNN
      DOUBLE PRECISION XIFMES,XISMES,MSMES,MUPMES,MSPMES,M3HMES
      DOUBLE PRECISION ALINP,XIFINP,XISINP,MSINP,MUPINP,MSPINP
      DOUBLE PRECISION MSUSYEFF,MMESS,N5,GAU
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION DETM,MUFAIL,SIGMU,Q2,Q2MIN
      DOUBLE PRECISION PI,ALP1,ALP2
      DOUBLE PRECISION G1MES,G2MES,G3MES,LMES,KMES,HTMES,HBMES,HLMES
      DOUBLE PRECISION LPPMES,LTTMES,LUMES,LDMES,LTMES,LBMES,LLMES,DHMES
      DOUBLE PRECISION MSREF,D,DMIN

      COMMON/FLAGS/OMGFLAG,MAFLAG,MOFLAG
      COMMON/GMUFLAG/GMUFLAG,HFLAG
      COMMON/GMSCEN/MSREF,D,DMIN,GMFLAG
      COMMON/SIGMU/SIGMU
      COMMON/RENSCALE/Q2
      COMMON/Q2FIX/Q2MIN,Q2FIX
      COMMON/DETM/DETM
      COMMON/MUFAIL/MUFAIL
      COMMON/DELMB/DELMB
      COMMON/MESCAL/MSUSYEFF,MMESS,N5
      COMMON/MESCOUP/G1MES,G2MES,G3MES,LMES,KMES,HTMES,HBMES,HLMES
      COMMON/MESEXT/XIFMES,XISMES,MSMES,MUPMES,MSPMES,M3HMES
      COMMON/MESGUT/LPPMES,LTTMES,LUMES,LDMES,LTMES,LBMES,LLMES,DHMES
      COMMON/INPPAR/ALINP,XIFINP,XISINP,MSINP,MUPINP,MSPINP
      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/MCMCPAR/MSUSYEFFCEN,MSUSYEFFDEV,MMESSCEN,MMESSDEV,
     . TBCEN,TBDEV,LCEN,LDEV,KCEN,KDEV,ALCEN,ALDEV,XIFCEN,
     . XIFDEV,XISCEN,XISDEV,MUPCEN,MUPDEV,MSPCEN,MSPDEV,
     . MSCEN,MSDEV,DHCEN,DHDEV,LPPCEN,LPPDEV,LTTCEN,LTTDEV,
     . LUCEN,LUDEV,LDCEN,LDDEV, LTCEN,LTDEV,LBCEN,LBDEV,
     . LLCEN,LLDEV,XCEN,XDEV,X,
     . MSUSYEFFMIN,MMESSMIN,TBMIN,LMIN,KMIN,ALMIN,
     . XIFMIN,XISMIN,MUPMIN,MSPMIN,MSMIN,DHMIN,
     . LPPMIN,LTTMIN,LUMIN,LDMIN,LTMIN,LBMIN,LLMIN
      COMMON/BOUNDS/MSUSYEFFN,MSUSYEFFNN,MMESSN,MMESSNN,TBN,TBNN,
     . LN,LNN,KN,KNN,ALN,ALNN,XIFN,XIFNN,XISN,XISNN,MUPN,MUPNN,
     . MSPN,MSPNN,MSN,MSNN,DHN,DHNN,LPPN,LPPNN,LTTN,LTTNN,LUN,LUNN,
     . LDN,LDNN,LTN,LTNN,LBN,LBNN,LLN,LLNN,MUN,MUNN
      COMMON/STEPS/NTOT,IDUM,TOTMIN,TOTMAX,NMAX

      PI=4d0*DATAN(1d0)

*   Initialization

      CALL INITIALIZE()
      DO I=1,NFL
       NFAIL(I)=0
      ENDDO
      TOT=0
      IFAIL=20
      IP=0

*   Reading of the input parameters

      CALL INPUT(PAR,NPAR)

*   Initialization of the range of parameters that has passed all tests

      MSUSYEFFN=1d99
      MSUSYEFFNN=-1d99
      MMESSN=1d99
      MMESSNN=-1d99
      TBN=1d99
      TBNN=-1d99
      LN=1d99
      LNN=-1d99
      ALN=1d99
      ALNN=-1d99
      MUN=1d99
      MUNN=-1d99
      XIFN=1d99
      XIFNN=-1d99
      KN=1d99
      KNN=-1d99
      MSN=1d99
      MSNN=-1d99
      XISN=1d99
      XISNN=-1d99
      MUPN=1d99
      MUPNN=-1d99
      MSPN=1d99
      MSPNN=-1d99
      DHN=1d99
      DHNN=-1d99
      LPPN=1d99
      LPPNN=-1d99
      LTTN=1d99
      LTTNN=-1d99
      LUN=1d99
      LUNN=-1d99
      LDN=1d99
      LDNN=-1d99
      LTN=1d99
      LTNN=-1d99
      LBN=1d99
      LBNN=-1d99
      LLN=1d99
      LLNN=-1d99

*   Beginning of the random scan

      DO ITOT=1,NTOT

 14   IF(IFAIL.EQ.20)THEN

       MSUSYEFF=MSUSYEFFCEN
       MMESS=MMESSCEN
       PAR(3)=TBCEN
       PAR(1)=LCEN
       IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-2)THEN
        XIFINP=XIFCEN
       ELSE
        PAR(2)=KCEN
       ENDIF
       IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3)THEN
        XISINP=XISCEN
       ELSE
        MSINP=MSCEN
       ENDIF
       MUPINP=MUPCEN
       MSPINP=MSPCEN
       DHMES=DHCEN
       IF(GMFLAG.EQ.1)THEN
        ALINP=-MSUSYEFF/(4d0*PI)**2*(DHMES
     .       +2d0*LPPMES+3d0*LTTMES+3d0*LUMES+3d0*LDMES)
       ELSE
        ALINP=ALCEN
       ENDIF
       LPPMES=LPPCEN
       LTTMES=LTTCEN
       LUMES=LUCEN
       LDMES=LDCEN
       LTMES=LTCEN
       LBMES=LBCEN
       LLMES=LLCEN

      ELSE

       IF(MSUSYEFFDEV.EQ.0d0)THEN
        MSUSYEFF=MSUSYEFFCEN
       ELSE
 101    MSUSYEFF=MSUSYEFFCEN
     .  +MAX(DABS(MSUSYEFFCEN),MSUSYEFFMIN)*MSUSYEFFDEV*GAU(IDUM)
        IF(MSUSYEFF.LE.0d0)GOTO 101
       ENDIF
       IF(MMESSDEV.EQ.0d0)THEN
        MMESS=MMESSCEN
       ELSE
 102    MMESS=MMESSCEN+MAX(DABS(MMESSCEN),MMESSMIN)*MMESSDEV*GAU(IDUM)
        IF(MMESS.LE.0d0)GOTO 102
       ENDIF
       IF(TBDEV.EQ.0d0)THEN
        PAR(3)=TBCEN
       ELSE
 103    PAR(3)=TBCEN+MAX(DABS(TBCEN),TBMIN)*TBDEV*GAU(IDUM)
       ENDIF
       IF(LDEV.EQ.0d0)THEN
        PAR(1)=LCEN
       ELSE
 104    PAR(1)=LCEN+MAX(DABS(LCEN),LMIN)*LDEV*GAU(IDUM)
        IF(PAR(1).LE.0d0)GOTO 104
       ENDIF
       IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-2)THEN
        IF(XIFDEV.EQ.0d0)THEN
         XIFINP=XIFCEN
        ELSE
 106     XIFINP=XIFCEN+MAX(DABS(XIFCEN),XIFMIN)*XIFDEV*GAU(IDUM)
        ENDIF
       ELSE
        IF(KDEV.EQ.0d0)THEN
         PAR(2)=KCEN
        ELSE
 107     PAR(2)=KCEN+MAX(DABS(KCEN),KMIN)*KDEV*GAU(IDUM)
        ENDIF
       ENDIF
       IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3)THEN
        IF(XISDEV.EQ.0d0)THEN
         XISINP=XISCEN
        ELSE
 108     XISINP=XISCEN+MAX(DABS(XISCEN),XISMIN)*XISDEV*GAU(IDUM)
        ENDIF
       ELSE
        IF(MSDEV.EQ.0d0)THEN
         MSINP=MSCEN
        ELSE
 109     MSINP=MSCEN+MAX(DABS(MSCEN),MSMIN)*MSDEV*GAU(IDUM)
        ENDIF
       ENDIF
       IF(MUPDEV.EQ.0d0)THEN
        MUPINP=MUPCEN
       ELSE
 110    MUPINP=MUPCEN+MAX(DABS(MUPCEN),MUPMIN)*MUPDEV*GAU(IDUM)
       ENDIF
       IF(MSPDEV.EQ.0d0)THEN
        MSPINP=MSPCEN
       ELSE
 111    MSPINP=MSPCEN+MAX(DABS(MSPCEN),MSPMIN)*MSPDEV*GAU(IDUM)
       ENDIF
       IF(DHDEV.EQ.0d0)THEN
        DHMES=DHCEN
       ELSE
 112    DHMES=DHCEN+MAX(DABS(DHCEN),DHMIN)*DHDEV*GAU(IDUM)
       ENDIF
       IF(LPPDEV.EQ.0d0)THEN
        LPPMES=LPPCEN
       ELSE
 113    LPPMES=LPPCEN+MAX(DABS(LPPCEN),LPPMIN)*LPPDEV*GAU(IDUM)
       ENDIF
       IF(LTTDEV.EQ.0d0)THEN
        LTTMES=LTTCEN
       ELSE
 114    LTTMES=LTTCEN+MAX(DABS(LTTCEN),LTTMIN)*LTTDEV*GAU(IDUM)
       ENDIF
       IF(LUDEV.EQ.0d0)THEN
        LUMES=LUCEN
       ELSE
 115    LUMES=LUCEN+MAX(DABS(LUCEN),LUMIN)*LUDEV*GAU(IDUM)
       ENDIF
       IF(LDDEV.EQ.0d0)THEN
        LDMES=LDCEN
       ELSE
 116    LDMES=LDCEN+MAX(DABS(LDCEN),LDMIN)*LDDEV*GAU(IDUM)
       ENDIF
       IF(LTDEV.EQ.0d0)THEN
        LTMES=LTCEN
       ELSE
 117    LTMES=LTCEN+MAX(DABS(LTCEN),LTMIN)*LTDEV*GAU(IDUM)
       ENDIF
       IF(LBDEV.EQ.0d0)THEN
        LBMES=LBCEN
       ELSE
 118    LBMES=LBCEN+MAX(DABS(LBCEN),LBMIN)*LBDEV*GAU(IDUM)
       ENDIF
       IF(LLDEV.EQ.0d0)THEN
        LLMES=LLCEN
       ELSE
 119    LLMES=LLCEN+MAX(DABS(LLCEN),LLMIN)*LLDEV*GAU(IDUM)
       ENDIF
       IF(GMFLAG.EQ.1)THEN
        ALINP=-MSUSYEFF/(4d0*PI)**2*(DHMES
     .       +2d0*LPPMES+3d0*LTTMES+3d0*LUMES+3d0*LDMES)
       ELSE
        IF(ALDEV.EQ.0d0)THEN
         ALINP=ALCEN
        ELSE
 105     ALINP=ALCEN+MAX(DABS(ALCEN),ALMIN)*ALDEV*GAU(IDUM)
        ENDIF
       ENDIF
      ENDIF

!      WRITE(0,*)"MAFLAG=",MAFLAG
!      WRITE(0,*)"MSUSYEFF =",MSUSYEFF
!      WRITE(0,*)"MMESS =",MMESS
!      WRITE(0,*)"TANB =",PAR(3)
!      WRITE(0,*)"LAMBDA =",PAR(1)
!      WRITE(0,*)"ALAMBDA=",ALINP
!      IF(MAFLAG.EQ.-3 .OR. MAFLAG.EQ.-4)WRITE(0,*)"KAPPA =",PAR(2)
!      IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-2)WRITE(0,*)"XIF =",XIFINP
!      IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3)WRITE(0,*)"XIS =",XISINP
!      IF(MAFLAG.EQ.-2 .OR. MAFLAG.EQ.-4)WRITE(0,*)"MS =",MSINP
!      WRITE(0,*)"MUP =",MUPINP
!      WRITE(0,*)"MSP =",MSPINP
!      WRITE(0,*)"DH =",DHMES
!      WRITE(0,*)"LPP =",LPPMES
!      WRITE(0,*)"LTT =",LTTMES
!      WRITE(0,*)"LU =",LUMES
!      WRITE(0,*)"LD =",LDMES
!      WRITE(0,*)"LT =",LTMES
!      WRITE(0,*)"LB =",LBMES
!      WRITE(0,*)"LL =",LLMES
!      WRITE(0,*)""
!      WRITE(0,*)""

*   Initialization of PROB

      DO I=1,NPROB
       PROB(I)=0d0
      ENDDO      

*   Guess parameters at Q2/MMESS

      MUFAIL=1d0
      IF(Q2FIX.EQ.0)THEN
       Q2=MAX(N5*MSUSYEFF**2/(4d0*PI)**4,Q2MIN)
      ENDIF
      IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-2)THEN
       PAR(2)=PAR(1)/5d0
      ELSE
       XIFMES=MSUSYEFF**2/(4d0*PI)**4
      ENDIF
      IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3)THEN
       MSMES=MSUSYEFF**2/(4d0*PI)**4
      ELSE
       XISMES=MSUSYEFF**3/(4d0*PI)**6
      ENDIF
      PAR(4)=SIGMU*N5*MSUSYEFF/(4d0*PI)**2
      IF(GMFLAG.EQ.1)THEN
       PAR(5)=-MSUSYEFF/(4d0*PI)**2
      ELSE
       PAR(5)=ALINP
      ENDIF
      IF(PAR(2).NE.0d0)THEN
       IF(GMFLAG.EQ.1)THEN
        PAR(6)=-3d0*MSUSYEFF/(4d0*PI)**2
       ELSE
        PAR(6)=3d0*ALINP
       ENDIF
      ELSE
       PAR(6)=0d0
      ENDIF
      DO I=7,11
       PAR(I)=N5*MSUSYEFF**2/(4d0*PI)**4
      ENDDO
      DO I=12,14
       PAR(I)=0d0
      ENDDO
      DO I=15,19
       PAR(I)=N5*MSUSYEFF**2/(4d0*PI)**4
      ENDDO
      DO I=20,22
       PAR(I)=N5*MSUSYEFF/(4d0*PI)**2
      ENDDO
      PAR(23)=MSUSYEFF/(4d0*PI)**2
      PAR(24)=MSUSYEFF/(4d0*PI)**2
      PAR(25)=0d0
      DELMB=.1d0
      IFAIL=0

!      WRITE(0,*)"Guesses"
!      WRITE(0,*)""
!      WRITE(0,*)"M1 =",PAR(20)
!      WRITE(0,*)"M2 =",PAR(21)
!      WRITE(0,*)"M3 =",PAR(22)
!      WRITE(0,*)"AL =",PAR(5)
!      WRITE(0,*)"AK =",PAR(6)
!      WRITE(0,*)"ATOP =",PAR(12)
!      WRITE(0,*)"ABOT =",PAR(13)
!      WRITE(0,*)"ATAU =",PAR(14)
!      WRITE(0,*)"AMUON =",PAR(25)
!      WRITE(0,*)"MQ3 =",PAR(7)
!      WRITE(0,*)"MU3 =",PAR(8)
!      WRITE(0,*)"MD3 =",PAR(9)
!      WRITE(0,*)"MQ =",PAR(15)
!      WRITE(0,*)"MU =",PAR(16)
!      WRITE(0,*)"MD =",PAR(17)
!      WRITE(0,*)"ML3 =",PAR(10)
!      WRITE(0,*)"ME3 =",PAR(11)
!      WRITE(0,*)"ML =",PAR(18)
!      WRITE(0,*)"ME =",PAR(19)
!      IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-2)THEN
!       WRITE(0,*)"K =",PAR(2)
!      ELSE
!      WRITE(0,*)"XIF =",XIFMES
!       ENDIF
!      IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3)THEN
!       WRITE(0,*)"MS =",MSMES
!      ELSE
!       WRITE(0,*)"XIS =",XISMES
!      ENDIF
!      WRITE(0,*)"MU =",PAR(4)
!      WRITE(0,*)""

*   Check for singular parameters l, tan(beta)

      IF(PAR(1)*PAR(3).EQ.0d0)THEN
       IFAIL=9
       GOTO 11
      ENDIF

*   Guess for couplings at MMESS

      CALL RGESGM(PAR,IFAIL)
      IF(IFAIL.NE.0)GOTO 11

*   External loop to compute the soft parameters at Q2

      ITER=0
 21   ITER=ITER+1
!      WRITE(0,*)"ITER =",ITER
!      WRITE(0,*)""
!      WRITE(0,*)""

      CALL RGESINVGM(PAR,IFAIL)
      IF(IFAIL.NE.0)GOTO 11

*   Internal loop to compute mu, k and ms

      ITERMU=0
 22   ITERMU=ITERMU+1
!      WRITE(0,*)"ITERMU =",ITERMU
!      WRITE(0,*)""
!      WRITE(0,*)""

      CALL RUNPAR(PAR)

      CALL MSFERM(PAR,IFAIL,0)
      IF(IFAIL.NE.0)GOTO 11

      CALL MINIMIZE(PAR,CHECK)
      IF(DETM.LT.0)IFAIL=17
      IF(IFAIL.NE.0)GOTO 11

      IF(CHECK.GT.1.D-12.AND.ITERMU.LE.100) GOTO 22
      IF(CHECK.GT.1.D-8) IFAIL=14
      IF(IFAIL.NE.0) GOTO 11
      
      CALL RGESGM(PAR,IFAIL)
      IF(PROB(27).NE.0.)IFAIL=11
      IF(IFAIL.NE.0) GOTO 11

      CALL RGESUNIGM(PAR,IFAIL,MESTEST)
      IF(IFAIL.NE.0)GOTO 11

      IF(MESTEST.GT.1.D-12.AND.ITER.LE.100) GOTO 21
      IF(MESTEST.GT.1.D-8) IFAIL=15
      IF(IFAIL.NE.0) GOTO 11

*   Computation of sfermion masses

      CALL MSFERM(PAR,IFAIL,1)
      IF(IFAIL.NE.0)GOTO 11
      
*   Computation of Higgs masses

      CALL MHIGGS(PAR,IFAIL)
      IF(IFAIL.EQ.-1)IFAIL=16
      IF(IFAIL.NE.0)GOTO 11

*   Computation of gluino mass

      CALL GLUINO(PAR)      
      
*   Computation of chargino masses

      CALL CHARGINO(PAR)

*   Computation of neutralino masses

      CALL NEUTRALINO(PAR)
      
*   Computation of Higgs + top branching ratios

      CALL DECAY(PAR)
      CALL TDECAY(PAR)

*   Exp. constraints on sparticles/Higgses

      CALL SUBEXP(PAR,PROB)
      CALL LHCHIG(PAR,PROB)
      CALL Higgs_CHI2(PAR,PROB)
      CALL CMS_TAUTAU(PAR,PROB)
      CALL CMS_AA_4MU(PAR,PROB)

*   b -> s gamma + B physics

      CALL BSG(PAR,PROB)

*   Anom. magn. moment of the Muon

      IF(GMUFLAG.EQ.1)CALL MAGNMU(PAR,PROB)

*   Global minimum?

      CALL CHECKMIN(PROB)

*  GUT scale

c      CALL RGESGMGUT(PROB,IFAIL)
c      IF(IFAIL.NE.0)GOTO 11

*   Computation of the fine-tuning

      CALL FTPAR(PAR,2)

*   Check for problems

      DO I=1,NPROB
       IF(PROB(I).NE.0d0)IFAIL=10
      ENDDO

*   Check MSMES=MSREF if GMFLAG=1  

      IF(GMFLAG.EQ.1)THEN
       MSREF=MSUSYEFF**2/(4d0*PI)**4*(7d0/5d0*DHMES**2-4d0*KMES**2*DHMES
     .      -(16d0*G3MES+6d0*G2MES+10d0/3d0*G1MES)*DHMES/5d0
     .      +LUMES*(8d0*LUMES+4d0*LMES-8d0*KMES+6d0*HTMES+6d0*LBMES
     .       +2d0*LLMES+12d0*LTTMES+16d0*LPPMES-2d0*G1MES-6d0*G2MES)
     .      +LDMES*(8d0*LDMES+4d0*LMES-8d0*KMES+6d0*LTMES+6d0*HBMES
     .       +2d0*HLMES+12d0*LTTMES+16d0*LPPMES-2d0*G1MES-6d0*G2MES)
     .      +LPPMES*(8d0*LPPMES-8d0*KMES-2d0*G1MES-6d0*G2MES)
     .      +LTTMES*(15d0*LTTMES-12d0*KMES-4d0/3d0*G1MES-16d0*G3MES)
     .      +8d0*LUMES*LDMES+12d0*LPPMES*LTTMES
     .      +16d0*DSQRT(LMES*LDMES*LUMES*LPPMES)
     .      +12d0*DSQRT(LMES*LDMES*LTMES*HTMES)
     .      +12d0*DSQRT(LMES*LUMES*LBMES*HBMES)
     .      +4d0*DSQRT(LMES*LUMES*LLMES*HLMES))
       D=2d0*DABS(MSMES-MSREF)/DABS(MSMES+MSREF)
       IF(D.GT.DMIN)THEN
        IF(IFAIL.EQ.0)IFAIL=18
        IF(IFAIL.EQ.10)IFAIL=19
       ENDIF
!       WRITE(0,*)MSREF,MSMES,D
      ENDIF

*   Recording of the results

11    CALL MCMCSTEPGM(PAR,PROB,NPROB,IFAIL)
      CALL OUTPUT(PAR,PROB,IFAIL)
      IF(IFAIL.EQ.0)THEN
       TOT=TOT+1
       MSUSYEFFN=MIN(MSUSYEFF,MSUSYEFFN)
       MSUSYEFFNN=MAX(MSUSYEFF,MSUSYEFFNN)
       MMESSN=MIN(MMESS,MMESSN)
       MMESSNN=MAX(MMESS,MMESSNN)
       TBN=MIN(PAR(3),TBN)
       TBNN=MAX(PAR(3),TBNN)
       LN=MIN(PAR(1),LN)
       LNN=MAX(PAR(1),LNN)
       KN=MIN(PAR(2),KN)
       KNN=MAX(PAR(2),KNN)
       ALN=MIN(ALINP,ALN)
       ALNN=MAX(ALINP,ALNN)
       MUN=MIN(PAR(4),MUN)
       MUNN=MAX(PAR(4),MUNN)
       XIFN=MIN(XIFMES,XIFN)
       XIFNN=MAX(XIFMES,XIFNN)
       XISN=MIN(XISMES,XISN)
       XISNN=MAX(XISMES,XISNN)
       MSN=MIN(MSMES,MSN)
       MSNN=MAX(MSMES,MSNN)
       MUPN=MIN(MUPMES,MUPN)
       MUPNN=MAX(MUPMES,MUPNN)
       MSPN=MIN(MSPMES,MSPN)
       MSPNN=MAX(MSPMES,MSPNN)
       DHN=MIN(DHMES,DHN)
       DHNN=MAX(DHMES,DHNN)
       LPPN=MIN(LPPMES,LPPN)
       LPPNN=MAX(LPPMES,LPPNN)
       LTTN=MIN(LTTMES,LTTN)
       LTTNN=MAX(LTTMES,LTTNN)
       LUN=MIN(LUMES,LUN)
       LUNN=MAX(LUMES,LUNN)
       LDN=MIN(LDMES,LDN)
       LDNN=MAX(LDMES,LDNN)
       LTN=MIN(LTMES,LTN)
       LTNN=MAX(LTMES,LTNN)
       LBN=MIN(LBMES,LBN)
       LBNN=MAX(LBMES,LBNN)
       LLN=MIN(LLMES,LLN)
       LLNN=MAX(LLMES,LLNN)
      ELSE
       NFAIL(IFAIL)=NFAIL(IFAIL)+1
      ENDIF

      IF(TOT.EQ.TOTMAX)THEN
       NTOT=ITOT
       GOTO 12
      ENDIF
      IF(IP.EQ.1)GOTO 13

      ENDDO

      IP=1
 13   IF(TOT.LT.TOTMIN .AND. NTOT.LT.NMAX)THEN
       NTOT=NTOT+1
       GOTO 14
      ENDIF

*   Summary of the scanning:
*   Number of points that passed/failed the tests
*   and range for scanned parameters

 12   CALL ERROR(TOT,NTOT,NFAIL)

      END


      SUBROUTINE INPUT(PAR,NPAR)
      
*******************************************************************
*   This subroutine reads SM and NMSSM parameters from input file   .
*******************************************************************

      IMPLICIT NONE

      CHARACTER CHINL*120,CHBLCK*60,CHDUM*120

      INTEGER I,NLINE,INL,ICH,IX,IVAL,Q2FIX,MCFLAG,TOTMIN,TOTMAX,NMAX
      INTEGER NTOT,ISEED,GMFLAG,N0,NLOOP,NBER,NPAR,ERR,VFLAG,Z3FLAG
      INTEGER OMGFLAG,MAFLAG,MOFLAG,PFLAG,GMUFLAG,HFLAG

      DOUBLE PRECISION PAR(*),VAL
      DOUBLE PRECISION ACC,XITLA,XLAMBDA,MC0,MB0,MT0
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION VUS,VCB,VUB,Q2,SIGMU,Q2MIN
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
      DOUBLE PRECISION MSUSYEFF,MMESS,N5,MSREF,D,DMIN

      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/CKM/VUS,VCB,VUB
      COMMON/ALS/XLAMBDA,MC0,MB0,MT0,N0
      COMMON/RENSCALE/Q2
      COMMON/Q2FIX/Q2MIN,Q2FIX
      COMMON/SIGMU/SIGMU
      COMMON/MCMCPAR/MSUSYEFFCEN,MSUSYEFFDEV,MMESSCEN,MMESSDEV,
     . TBCEN,TBDEV,LCEN,LDEV,KCEN,KDEV,ALCEN,ALDEV,XIFCEN,
     . XIFDEV,XISCEN,XISDEV,MUPCEN,MUPDEV,MSPCEN,MSPDEV,
     . MSCEN,MSDEV,DHCEN,DHDEV,LPPCEN,LPPDEV,LTTCEN,LTTDEV,
     . LUCEN,LUDEV,LDCEN,LDDEV, LTCEN,LTDEV,LBCEN,LBDEV,
     . LLCEN,LLDEV,XCEN,XDEV,X,
     . MSUSYEFFMIN,MMESSMIN,TBMIN,LMIN,KMIN,ALMIN,
     . XIFMIN,XISMIN,MUPMIN,MSPMIN,MSMIN,DHMIN,
     . LPPMIN,LTTMIN,LUMIN,LDMIN,LTMIN,LBMIN,LLMIN
      COMMON/STEPS/NTOT,ISEED,TOTMIN,TOTMAX,NMAX
      COMMON/MESCAL/MSUSYEFF,MMESS,N5
      COMMON/FLAGS/OMGFLAG,MAFLAG,MOFLAG
      COMMON/PFLAG/PFLAG
      COMMON/GMUFLAG/GMUFLAG,HFLAG
      COMMON/MCFLAG/MCFLAG
      COMMON/GMSCEN/MSREF,D,DMIN,GMFLAG
      COMMON/VFLAG/VFLAG

*   MAXIMUM DEV (%) FOR MSMES/MSREF (SCENARIO 2)
      DMIN=1d-1

*   INITIALIZATION OF THE SUSY PARAMETERS
      DO I=1,NPAR
       PAR(I)=0d0
      ENDDO      

*   INITIALIZATION OF THE SCANNING PARAMETERS
      N5=1d99
      SIGMU=1d99
      MSUSYEFFCEN=1d99
      MSUSYEFFDEV=1d99
      MMESSCEN=1d99
      MMESSDEV=1d99
      TBCEN=1d99
      TBDEV=1d99
      LCEN=1d99
      LDEV=1d99
      KCEN=1d99
      KDEV=1d99
      ALCEN=1d99
      ALDEV=1d99
      XIFCEN=1d99
      XIFDEV=1d99
      XISCEN=1d99
      XISDEV=1d99
      MUPCEN=0d0
      MUPDEV=1d99
      MSPCEN=0d0
      MSPDEV=1d99
      MSCEN=1d99
      MSDEV=1d99
      DHCEN=0d0
      DHDEV=1d99
      LPPCEN=0d0
      LPPDEV=1d99
      LTTCEN=0d0
      LTTDEV=1d99
      LUCEN=0d0
      LUDEV=1d99
      LDCEN=0d0
      LDDEV=1d99
      LTCEN=0d0
      LTDEV=1d99
      LBCEN=0d0
      LBDEV=1d99
      LLCEN=0d0
      LLDEV=1d99
      XCEN=1d99
      XDEV=1d-1
      LMIN=0d0
      KMIN=0d0
      ALMIN=0d0
      XIFMIN=0d0
      XISMIN=0d0
      MUPMIN=0d0
      MSPMIN=0d0
      MSMIN=0d0
      DHMIN=0d0
      LPPMIN=0d0
      LTTMIN=0d0
      LUMIN=0d0
      LDMIN=0d0
      LTMIN=0d0
      LBMIN=0d0
      LLMIN=0d0
      NTOT=0
      TOTMIN=0
      TOTMAX=1000000

*   DEFAULT VALUES FOR FLAGS
      OMGFLAG=0
      PFLAG=0
      GMUFLAG=1
      HFLAG=0
      GMFLAG=0
      MCFLAG=0
      VFLAG=0
      MOFLAG=0
      Z3FLAG=2

*   DEFAULT VALUE FOR THE RANDOM SEED
      ISEED=-1

*   DEFAULT VALUE FOR THE RENSCALE Q2
      Q2=0d0

*   INITIALIZE READ LOOP
      NLINE=0
      CHBLCK=' '

*   START TO READ NEW LINE INTO CHINL
 21   CHINL=' '

*   LINE NUMBER
      NLINE=NLINE+1

      READ(5,'(A120)',END=29,ERR=999) CHINL
      
*   CHECK FOR EMPTY OR COMMENT LINES
      IF(CHINL.EQ.' '.OR.CHINL(1:1).EQ.'#'
     .  .OR.CHINL(1:1).EQ.'*') GOTO 21

*   FORCE UPPER CASE LETTERS IN CHINL (AS REQUIRED BELOW)
      INL=0
 22   INL=INL+1
      IF(CHINL(INL:INL).NE.'#')THEN
       DO ICH=97,122
        IF(CHINL(INL:INL).EQ.CHAR(ICH)) CHINL(INL:INL)=CHAR(ICH-32)
       END DO
       IF(INL.LT.120) GOTO 22
      ENDIF

*   CHECK FOR BLOCK STATEMENT
      IF(CHINL(1:1).EQ.'B')THEN
       READ(CHINL,'(A6,A)',ERR=999) CHDUM,CHBLCK
       GOTO 21
      ENDIF

*   CHECK FOR NMSSM MODEL IN MODSEL
*   IF THE RELIC DENSITY SHOULD BE COMPUTED
*   THE BLOCK MODSEL MUST CONTAIN THE LINE "  9     1    "
      IF(CHBLCK(1:6).EQ.'MODSEL')THEN
       READ(CHINL,*,ERR=999) IX,IVAL
       IF(IX.EQ.1) Z3FLAG=IVAL
       IF(IX.EQ.8) PFLAG=IVAL
       IF(IX.EQ.11) GMUFLAG=IVAL
       IF(IX.EQ.12) GMFLAG=IVAL
       IF(IX.EQ.14) VFLAG=IVAL

*   READ SMINPUTS
      ELSEIF(CHBLCK(1:8).EQ.'SMINPUTS')THEN
       READ(CHINL,*,ERR=999) IX,VAL
       IF(IX.EQ.2) GF=VAL
       IF(IX.EQ.3) ALSMZ=VAL
       IF(IX.EQ.4) MZ=VAL
       IF(IX.EQ.5) MB=VAL
       IF(IX.EQ.6) MT=VAL
       IF(IX.EQ.7) MTAU=VAL
 
*   READ MESS PARAMETERS, SIGMU, Q2 AND TANBETA
      ELSEIF(CHBLCK(1:6).EQ.'MINPAR')THEN
       READ(CHINL,*,ERR=999) IX,VAL
       IF(IX.EQ.0) Q2=VAL**2
       IF(IX.EQ.1) MSUSYEFFCEN=VAL
       IF(IX.EQ.16) MSUSYEFFDEV=VAL
       IF(IX.EQ.17) MSUSYEFFMIN=VAL
       IF(IX.EQ.2) MMESSCEN=VAL
       IF(IX.EQ.26) MMESSDEV=VAL
       IF(IX.EQ.27) MMESSMIN=VAL
       IF(IX.EQ.3) TBCEN=VAL
       IF(IX.EQ.36) TBDEV=VAL
       IF(IX.EQ.37) TBMIN=VAL
       IF(IX.EQ.4)  SIGMU=VAL
       IF(IX.EQ.5)  N5=VAL
 
*   READ EXTPAR
      ELSEIF(CHBLCK(1:6).EQ.'EXTPAR')THEN
       READ(CHINL,*,ERR=999) IX,VAL
       IF(IX.EQ.-1) XDEV=VAL
       IF(IX.EQ.0) DMIN=VAL
       IF(IX.EQ.61) LCEN=VAL
       IF(IX.EQ.616) LDEV=VAL
       IF(IX.EQ.617) LMIN=VAL
       IF(IX.EQ.62) KCEN=VAL
       IF(IX.EQ.626) KDEV=VAL
       IF(IX.EQ.627) KMIN=VAL
       IF(IX.EQ.63) ALCEN=VAL
       IF(IX.EQ.636) ALDEV=VAL
       IF(IX.EQ.637) ALMIN=VAL
       IF(IX.EQ.66) XIFCEN=VAL
       IF(IX.EQ.666) XIFDEV=VAL
       IF(IX.EQ.667) XIFMIN=VAL
       IF(IX.EQ.67) XISCEN=VAL
       IF(IX.EQ.676) XISDEV=VAL
       IF(IX.EQ.677) XISMIN=VAL
       IF(IX.EQ.68) MUPCEN=VAL
       IF(IX.EQ.686) MUPDEV=VAL
       IF(IX.EQ.687) MUPMIN=VAL
       IF(IX.EQ.69) MSPCEN=VAL
       IF(IX.EQ.696) MSPDEV=VAL
       IF(IX.EQ.697) MSPMIN=VAL
       IF(IX.EQ.70) MSCEN=VAL
       IF(IX.EQ.706) MSDEV=VAL
       IF(IX.EQ.707) MSMIN=VAL
       IF(IX.EQ.71) DHCEN=VAL
       IF(IX.EQ.716) DHDEV=VAL
       IF(IX.EQ.717) DHMIN=VAL
       IF(IX.EQ.73) LPPCEN=VAL
       IF(IX.EQ.736) LPPDEV=VAL
       IF(IX.EQ.737) LPPMIN=VAL
       IF(IX.EQ.74) LTTCEN=VAL
       IF(IX.EQ.746) LTTDEV=VAL
       IF(IX.EQ.747) LTTMIN=VAL
       IF(IX.EQ.75) LUCEN=VAL
       IF(IX.EQ.756) LUDEV=VAL
       IF(IX.EQ.757) LUMIN=VAL
       IF(IX.EQ.76) LDCEN=VAL
       IF(IX.EQ.766) LDDEV=VAL
       IF(IX.EQ.767) LDMIN=VAL
       IF(IX.EQ.77) LTCEN=VAL
       IF(IX.EQ.776) LTDEV=VAL
       IF(IX.EQ.777) LTMIN=VAL
       IF(IX.EQ.78) LBCEN=VAL
       IF(IX.EQ.786) LBDEV=VAL
       IF(IX.EQ.787) LBMIN=VAL
       IF(IX.EQ.79) LLCEN=VAL
       IF(IX.EQ.796) LLDEV=VAL
       IF(IX.EQ.797) LLMIN=VAL

*   READ STEPS
       ELSEIF(CHBLCK(1:6).EQ.'STEPS')THEN
       READ(CHINL,*,ERR=999) IX,IVAL
       IF(IX.EQ.0) NTOT=IVAL
       IF(IX.EQ.1) ISEED=IVAL
       IF(IX.EQ.2) HFLAG=IVAL
       IF(IX.EQ.3) MCFLAG=IVAL
       IF(IX.EQ.4) TOTMIN=IVAL
       IF(IX.EQ.5) TOTMAX=IVAL
       IF(IX.EQ.6) NMAX=IVAL

      ENDIF

      GOTO 21

*   END OF READING FROM INPUT FILE

*   Check for errors

 29   ERR=0
      IF(GMFLAG.NE.0 .AND. GMFLAG.NE.1)THEN
       WRITE(0,1)"GMFLAG IS EITHER 0 OR 1"
       ERR=1
      ENDIF
      IF(XDEV.EQ.1d99)THEN
       WRITE(0,1)"XDEV MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(MSUSYEFFCEN.EQ.1d99)THEN
       WRITE(0,1)"MSUSYEFFCEN MUST BE GIVEN IN BLOCK MINPAR"
       ERR=1
      ENDIF
      IF(MMESSCEN.EQ.1d99)THEN
       WRITE(0,1)"MMESSCEN MUST BE GIVEN IN BLOCK MINPAR"
       ERR=1
      ENDIF
      IF(TBCEN.EQ.1d99)THEN
       WRITE(0,1)"TBCEN MUST BE GIVEN IN BLOCK MINPAR"
       ERR=1
      ENDIF
      IF(DABS(SIGMU).NE.1d0)THEN
       WRITE(0,1)"SIGMU IS EITHER 1 OR -1 IN BLOCK MINPAR"
       ERR=1
      ENDIF
      IF(N5.EQ.1d99)THEN
       WRITE(0,1)"N5 MUST BE GIVEN IN BLOCK MINPAR"
       ERR=1
      ENDIF
      IF(N5.LE.0d0.OR. ANINT(N5).NE.N5)THEN
       WRITE(0,1)"N5 MUST BE A POSITIVE INTEGER IN BLOCK MINPAR"
       ERR=1
      ENDIF
      IF(LCEN.EQ.1d99)THEN
       WRITE(0,1)"LCEN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(ALCEN.EQ.1d99 .AND. ALDEV.NE.1d99)THEN
       WRITE(0,1)"ALCEN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(KCEN.EQ.1d99 .AND. KDEV.NE.1d99)THEN
       WRITE(0,1)"KCEN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(XIFCEN.EQ.1d99 .AND. XIFDEV.NE.1d99)THEN
       WRITE(0,1)"XIFCEN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(XISCEN.EQ.1d99 .AND. XISDEV.NE.1d99)THEN
       WRITE(0,1)"XISCEN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(MSCEN.EQ.1d99 .AND. MSDEV.NE.1d99)THEN
       WRITE(0,1)"MSCEN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(XIFCEN.NE.1d99 .AND. KCEN.NE.1d99)THEN
       WRITE(0,1)"BOTH XIF AND KAPPA CANNOT BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(XISCEN.NE.1d99 .AND. MSCEN.NE.1d99)THEN
       WRITE(0,1)"BOTH XIS AND MS^2 CANNOT BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(GMFLAG.EQ.1 .AND. ALCEN.NE.1d99)THEN
       WRITE(0,1)"GMFLAG=1 => ALCEN CANNOT BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(GMFLAG.EQ.1 .AND. MSCEN.NE.1d99)THEN
       WRITE(0,1)"GMFLAG=1 => MSCEN CANNOT BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF

*   Set default values

      IF(KCEN.EQ.1d99 .AND. XIFCEN.EQ.1d99)XIFCEN=0d0
      IF(MSCEN.EQ.1d99 .AND. XISCEN.EQ.1d99)XISCEN=0d0
      IF(ALCEN.EQ.1d99)ALCEN=0d0

      IF(MSUSYEFFDEV.EQ.1d99)MSUSYEFFDEV=0d0
      IF(MMESSDEV.EQ.1d99)MMESSDEV=0d0
      IF(TBDEV.EQ.1d99)TBDEV=0d0
      IF(LDEV.EQ.1d99)LDEV=0d0
      IF(KDEV.EQ.1d99)KDEV=0d0
      IF(ALDEV.EQ.1d99)ALDEV=0d0
      IF(XIFDEV.EQ.1d99)XIFDEV=0d0
      IF(XISDEV.EQ.1d99)XISDEV=0d0
      IF(MUPDEV.EQ.1d99)MUPDEV=0d0
      IF(MSPDEV.EQ.1d99)MSPDEV=0d0
      IF(MSDEV.EQ.1d99)MSDEV=0d0
      IF(DHDEV.EQ.1d99)DHDEV=0d0
      IF(LPPDEV.EQ.1d99)LPPDEV=0d0
      IF(LTTDEV.EQ.1d99)LTTDEV=0d0
      IF(LUDEV.EQ.1d99)LUDEV=0d0
      IF(LDDEV.EQ.1d99)LDDEV=0d0
      IF(LTDEV.EQ.1d99)LTDEV=0d0
      IF(LBDEV.EQ.1d99)LBDEV=0d0
      IF(LLDEV.EQ.1d99)LLDEV=0d0

*   Set MAFLAG

      IF(KCEN.EQ.1d99 .AND. MSCEN.EQ.1d99)MAFLAG=-1
      IF(KCEN.EQ.1d99 .AND. XISCEN.EQ.1d99)MAFLAG=-2
      IF(XIFCEN.EQ.1d99 .AND. MSCEN.EQ.1d99)MAFLAG=-3
      IF(XIFCEN.EQ.1d99 .AND. XISCEN.EQ.1d99)MAFLAG=-4

*   Total number of points

      IF(NTOT.LE.0)THEN
       WRITE(0,1)"WRONG NUMBER OF POINTS IN BLOCK STEPS"
       ERR=1
      ENDIF
      IF(TOTMIN.GT.TOTMAX)THEN
       WRITE(0,1)"TOTMIN must be smaller than TOTMAX"
       ERR=1
      ENDIF

*   Check for Z3 breaking terms

      IF(MAFLAG.EQ.-2 .OR. MAFLAG.EQ.-3 .OR. MAFLAG.EQ.-4 .OR.
     . MUPCEN.NE.0d0 .OR. MUPDEV.NE.0d0 .OR. MSPCEN.NE.0d0 .OR. 
     . MSPDEV.NE.0d0 .OR. XIFCEN.NE.0d0 .OR. XIFDEV.NE.0d0 .OR.
     . XISCEN.NE.0d0 .OR. XISDEV.NE.0d0)THEN
       IF(PFLAG.NE.0)THEN
        WRITE(0,1)"HIGGS MASS PRECISION = 1 OR 2 ONLY FOR Z3-NMSSM"
        ERR=1
       ENDIF
       IF(Z3FLAG.GT.2)THEN
        WRITE(0,1)"PRESENCE OF Z3 BREAKING TERMS"
        ERR=1
       ENDIF
      ENDIF

*   Stop if error

      IF(ERR.EQ.1)THEN
       WRITE(0,1)"ERROR IN INPUT FILE"
       STOP 1
      ENDIF

*   Set Q2MIN, Q2FIX:
      Q2MIN=100d0**2
      Q2FIX=1
      IF(Q2.LE.Q2MIN)THEN
       Q2FIX=0
      ENDIF

*   Initialization for ALPHAS and RUNM (as in hdecay)
*   The bottom quark pole mass MBP is set in INIT and can be changed
*   only there (changing its running mass MB above has no effect
*   on MBP, since one would have to compute alpha_s(MB) first)

      MC0=MC
      MB0=MBP
      MT0=MT
      N0=5
      NLOOP=2
      NBER=18
      ACC=1d-8
      XLAMBDA=XITLA(NLOOP,ALSMZ,ACC)
      CALL ALSINI(ACC)
      CALL BERNINI(NBER)

*    g1,g2  and sin(theta)^2 in the on-shell scheme in terms of 
*    GF, MZ(pole) and MW(pole)

      g2=4d0*DSQRT(2d0)*GF*MW**2
      g1=4d0*DSQRT(2d0)*GF*(MZ**2-MW**2)
      S2TW=1d0-(MW/MZ)**2

      RETURN

 999  WRITE(0,1)"READ ERROR ON LINE:", NLINE
      WRITE(0,*)CHINL(1:80)
      STOP 1

 1    FORMAT(A)

      END


      SUBROUTINE OUTPUT(PAR,PROB,IFAIL) 

*********************************************************************      
*   Subroutine writing all the results in the the output file.
*********************************************************************      
 
      IMPLICIT NONE 
 
      INTEGER GMFLAG,I,J,NRES,IRES,NSUSY,NGUT,NMES,IMAX,IFAIL
      PARAMETER (NSUSY=14,NGUT=21,NMES=22,IMAX=200)
 
      DOUBLE PRECISION RES(IMAX),PAR(*),PROB(*),SIG(3,10)
      DOUBLE PRECISION SMASS(3),PMASS(2),CMASS,SCOMP(3,3),PCOMP(2,2)
      DOUBLE PRECISION MGL,MCHA(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5),BRSS(5)
      DOUBLE PRECISION BRCC(5),BRBB(5),BRTT(5),BRWW(3),BRZZ(3)
      DOUBLE PRECISION BRGG(5),BRZG(5),BRHHH(4),BRHAA(3,3)
      DOUBLE PRECISION BRHCHC(3),BRHAZ(3,2),BRAHA(3),BRAHZ(2,3)
      DOUBLE PRECISION BRHCW(5),BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5)
      DOUBLE PRECISION HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC
      DOUBLE PRECISION HCBRBT,HCBRWH(5),HCBRWHT,HCBRNC(5,2)
      DOUBLE PRECISION HCBRSQ(5),HCBRSL(3),HCBRSUSY,HCWIDTH
      DOUBLE PRECISION CU(5),CD(5),CV(3),CJ(5),CG(5)
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION VUS,VCB,VUB,MSREF,D,DMIN
      DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
      DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
      DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU,Q2
      DOUBLE PRECISION RMST1,RMST2,S2T,RMSB1,RMSB2,S2B,XT,XB
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION MHUS,MHDS,MSS
      DOUBLE PRECISION G1MES,G2MES,G3MES,LMES,KMES,HTMES,HBMES,HLMES
      DOUBLE PRECISION LPPMES,LTTMES,LUMES,LDMES,LTMES,LBMES,LLMES,DHMES
      DOUBLE PRECISION SIGMU,MGUT
      DOUBLE PRECISION MHUQ,MHDQ,MSQ,LQ,KQ,ALQ,AKQ,MUQ,NUQ
      DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ,QSTSB
      DOUBLE PRECISION BRSG,BRSGmax,BRSGmin,DMd,DMdmin,DMdmax,DMs,
     . DMsmax,DMsmin,BRBMUMU,BRBMUMUmax,BRBMUMUmin,BRBtaunu,
     . BRBtaunumax,BRBtaunumin
      DOUBLE PRECISION delmagmu,damumin,damumax,amuthmax,amuthmin
      DOUBLE PRECISION MSUSYEFF,MMESS,N5
      DOUBLE PRECISION M1MES,M2MES,M3MES,ALMES,AKMES,ATMES,ABMES,
     . ATAUMES,AMUMES,MHUMES,MHDMES,MQ3MES,MU3MES,MD3MES,
     . MQMES,MUMES,MDMES,ML3MES,ME3MES,MLMES,MEMES
      DOUBLE PRECISION XIFMES,XISMES,MSMES,MUPMES,MSPMES,M3HMES
      DOUBLE PRECISION XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY
      DOUBLE PRECISION ALINP,XIFINP,XISINP,MSINP,MUPINP,MSPINP
      DOUBLE PRECISION G1GUT,G2GUT,G3GUT,LGUT,KGUT,HTGUT,HBGUT,HLGUT
      DOUBLE PRECISION brtopbw,brtopbh,brtopneutrstop(5,2),toptot
      DOUBLE PRECISION FTSUSY(NSUSY+2),FTGUT(NGUT+2),FTMES(NMES+2)
      DOUBLE PRECISION DELMB,PX,PA(6),PB(2),PL(7),PK(8),MH(3),MMH(3)
      DOUBLE PRECISION DMH(3),MA(2),MMA(2),DMA(2),MHC,MMHC,DMHC
      DOUBLE PRECISION MHmin,MHmax,chi2max,chi2gam,chi2bb,chi2zz

      COMMON/SIGMU/SIGMU
      COMMON/INPPAR/ALINP,XIFINP,XISINP,MSINP,MUPINP,MSPINP
      COMMON/MESCAL/MSUSYEFF,MMESS,N5
      COMMON/MESEXT/XIFMES,XISMES,MSMES,MUPMES,MSPMES,M3HMES
      COMMON/SUSYEXT/XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY
      COMMON/SOFTMES/M1MES,M2MES,M3MES,ALMES,AKMES,ATMES,ABMES,
     . ATAUMES,AMUMES,MHUMES,MHDMES,MQ3MES,MU3MES,MD3MES,
     . MQMES,MUMES,MDMES,ML3MES,ME3MES,MLMES,MEMES
      COMMON/BRN/BRJJ,BREE,BRMM,BRLL,BRSS,BRCC,BRBB,BRTT,BRWW,
     . BRZZ,BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     . BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     . BRSUSY,WIDTH
      COMMON/BRC/HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC,
     . HCBRBT,HCBRWH,HCBRWHT,HCBRNC,HCBRSQ,HCBRSL,
     . HCBRSUSY,HCWIDTH
      COMMON/BRSG/BRSG,BRSGmax,BRSGmin,DMd,DMdmin,DMdmax,DMs,
     . DMsmax,DMsmin,BRBMUMU,BRBMUMUmax,BRBMUMUmin,BRBtaunu,
     . BRBtaunumax,BRBtaunumin
      COMMON/MAGMU/delmagmu,damumin,damumax,amuthmax,amuthmin
      COMMON/BR_top2body/brtopbw,brtopbh,brtopneutrstop
      COMMON/topwidth/toptot
      COMMON/REDCOUP/CU,CD,CV,CJ,CG
      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,NEU
      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/CKM/VUS,VCB,VUB
      COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     . MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     . CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
      COMMON/RADCOR/RMST1,RMST2,S2T,RMSB1,RMSB2,S2B,XT,XB
      COMMON/RENSCALE/Q2
      COMMON/STSBSCALE/QSTSB
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/SUSYMH/MHUS,MHDS,MSS
      COMMON/QMHIGGS/MHUQ,MHDQ,MSQ
      COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      COMMON/QPAR/LQ,KQ,ALQ,AKQ,MUQ,NUQ
      COMMON/MESCOUP/G1MES,G2MES,G3MES,LMES,KMES,HTMES,HBMES,HLMES
      COMMON/MESGUT/LPPMES,LTTMES,LUMES,LDMES,LTMES,LBMES,LLMES,DHMES
      COMMON/GUTCOUP/G1GUT,G2GUT,G3GUT,LGUT,KGUT,HTGUT,HBGUT,HLGUT
      COMMON/MGUT/MGUT
      COMMON/FINETUN/FTSUSY,FTGUT,FTMES
      COMMON/EFFHIGM/MH,MMH,DMH,MA,MMA,DMA,MHC,MMHC,DMHC
      COMMON/EFFCOUP/PX,PA,PB,PL,PK
      COMMON/DELMB/DELMB
      COMMON/LHCSIG/SIG
      COMMON/HIGGSFIT/MHmin,MHmax,chi2max,chi2gam,chi2bb,chi2zz
      COMMON/GMSCEN/MSREF,D,DMIN,GMFLAG

      IF(IFAIL.NE.0)RETURN

      IRES=11
      NRES=106+IRES

      RES(1)=PAR(1)
      RES(2)=PAR(2)
      RES(3)=PAR(3)
      RES(4)=PAR(4)
      RES(5)=MMESS
      RES(6)=MSUSYEFF
      RES(7)=ALMES
      RES(8)=MSMES
      RES(9)=XIFMES
      RES(10)=XISMES
      RES(11)=DHMES

      DO I=1,3
       RES(IRES-6+7*I)=SMASS(I)
       RES(IRES-5+7*I)=CV(I)
       RES(IRES-4+7*I)=CG(I)
       RES(IRES-3+7*I)=CJ(I)
       RES(IRES-2+7*I)=CU(I)
       RES(IRES-1+7*I)=CD(I)
       RES(IRES+7*I)=SCOMP(I,3)**2
      ENDDO
      DO I=1,2
       RES(IRES+16+6*I)=PMASS(I)
       RES(IRES+17+6*I)=CG(I+3)
       RES(IRES+18+6*I)=CJ(I+3)
       RES(IRES+19+6*I)=CU(I+3)
       RES(IRES+20+6*I)=CD(I+3)
       RES(IRES+21+6*I)=PCOMP(I,2)**2
      ENDDO
      DO I=1,3
       RES(IRES+33+I)=BRHAA(I,1)
      ENDDO
      DO I=1,4
       RES(IRES+36+I)=BRHHH(I)
      ENDDO
      DO I=1,5
       RES(IRES+40+I)=BRNEU(I,1,1)
      ENDDO
      DO I=1,3
       DO J=1,10
        RES(IRES+35+10*I+J)=SIG(I,J)
       ENDDO
      ENDDO
      RES(IRES+76)=chi2gam
      RES(IRES+77)=chi2bb
      RES(IRES+78)=chi2zz
      RES(IRES+79)=CMASS
      DO I=1,5
       RES(IRES+79+I)=DABS(MNEU(I))
      ENDDO
      RES(IRES+85)=NEU(1,1)**2
      RES(IRES+86)=NEU(1,3)**2+NEU(1,4)**2
      RES(IRES+87)=DABS(NEU(1,3)**2-NEU(1,4)**2)
      RES(IRES+88)=NEU(1,5)**2
      RES(IRES+89)=MCHA(1)
      RES(IRES+90)=MCHA(2)
      RES(IRES+91)=MGL
      RES(IRES+92)=MUR
      RES(IRES+93)=MUL
      RES(IRES+94)=MDR
      RES(IRES+95)=MDL
      RES(IRES+96)=MLR
      RES(IRES+97)=MLL
      RES(IRES+98)=MNL
      RES(IRES+99)=MST1
      RES(IRES+100)=MST2
      RES(IRES+101)=MSB1
      RES(IRES+102)=MSB2
      RES(IRES+103)=MSL1
      RES(IRES+104)=MSL2
      RES(IRES+105)=MSNT
      RES(IRES+106)=FTMES(NMES+1)

      WRITE(6,11)(RES(I),I=1,NRES)
 11   FORMAT(200E14.6)

      END


      SUBROUTINE ERROR(TOT,NTOT,NFAIL)

*********************************************************************      
*   Subroutine for the error file. It contains a summary of the scan:
*   Number of points that passed/failed the tests
*   and ranges for scanned parameters that passed the tests
*********************************************************************      
 
      IMPLICIT NONE

      INTEGER OMGFLAG,MAFLAG,MOFLAG,GMFLAG,I,S,TOT,NTOT,NFAIL(*)

      DOUBLE PRECISION MSUSYEFFN,MSUSYEFFNN,MMESSN,MMESSNN,TBN,TBNN
      DOUBLE PRECISION LN,LNN,KN,KNN,ALN,ALNN,XIFN,XIFNN,XISN,XISNN
      DOUBLE PRECISION MUPN,MUPNN,MSPN,MSPNN,MSN,MSNN,DHN,DHNN
      DOUBLE PRECISION LPPN,LPPNN,LTTN,LTTNN,LUN,LUNN,LDN,LDNN
      DOUBLE PRECISION LTN,LTNN,LBN,LBNN,LLN,LLNN,MUN,MUNN
      DOUBLE PRECISION MSREF,D,DMIN

      COMMON/FLAGS/OMGFLAG,MAFLAG,MOFLAG
      COMMON/BOUNDS/MSUSYEFFN,MSUSYEFFNN,MMESSN,MMESSNN,TBN,TBNN,
     . LN,LNN,KN,KNN,ALN,ALNN,XIFN,XIFNN,XISN,XISNN,MUPN,MUPNN,
     . MSPN,MSPNN,MSN,MSNN,DHN,DHNN,LPPN,LPPNN,LTTN,LTTNN,LUN,LUNN,
     . LDN,LDNN,LTN,LTNN,LBN,LBNN,LLN,LLNN,MUN,MUNN
      COMMON/GMSCEN/MSREF,D,DMIN,GMFLAG

      WRITE(0,20)"Number of points:                       "
      WRITE(0,*)
      WRITE(0,20)"  scanned                               ",NTOT
      WRITE(0,20)"  l, tan(beta) or mu=0                  ",NFAIL(9)
      WRITE(0,20)"  no electroweak symmetry breaking      ",NFAIL(17)
      S=0
      DO I=1,7
       S=S+NFAIL(I)
      ENDDO
      WRITE(0,20)"  with mh1^2 or ma1^2 or mhc^2 < 0      ",S
      WRITE(0,20)"  with m_sfermion^2 < 0                 ",NFAIL(8)
      WRITE(0,20)"  violating phenomenological constraints",NFAIL(10)
      IF(GMFLAG.EQ.1)THEN
       WRITE(0,20)"  MSMES=/=MSREF                         ",NFAIL(18)
       WRITE(0,20)"  MSMES=/=MSREF + violating constraints ",NFAIL(19)
      ENDIF
      S=NFAIL(11)+NFAIL(12)+NFAIL(13)
      WRITE(0,20)"  RGE integration problem               ",S
      S=NFAIL(14)+NFAIL(15)+NFAIL(16)
      WRITE(0,20)"  convergence problem                   ",S
      WRITE(0,*)
      WRITE(0,20)"Remaining good points                   ",TOT
      IF(TOT.GT.0)THEN
       WRITE(0,*)
       WRITE(0,20)"Parameter ranges for good points:       "
       WRITE(0,*)
       WRITE(0,30)" MSUSYEFF: ",MSUSYEFFN,MSUSYEFFNN
       WRITE(0,30)" MMESS: ",MMESSN,MMESSNN
       WRITE(0,30)" TANB: ",TBN,TBNN
       WRITE(0,30)" LAMBDA: ",LN,LNN
       WRITE(0,30)" ALAMBDA: ",ALN,ALNN
       IF(GMFLAG.EQ.1)
     .  WRITE(0,20)"(ALAMBDA is not an input parameter)"
       WRITE(0,30)"MUEFF: ",MUN,MUNN
       WRITE(0,20)"(MU is not an input parameter)"
       IF(MAFLAG.EQ.-3 .OR. MAFLAG.EQ.-4)THEN
        WRITE(0,30)" KAPPA: ",KN,KNN
        WRITE(0,30)" XIF: ",XIFN,XIFNN
        WRITE(0,20)"(XIF is not an input parameter)"
       ELSE
        WRITE(0,30)" KAPPA: ",KN,KNN
        WRITE(0,20)"(KAPPA is not an input parameter)"
        WRITE(0,30)" XIF: ",XIFN,XIFNN
       ENDIF
       IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3)THEN
        WRITE(0,30)" XIS: ",XISN,XISNN
        WRITE(0,30)" MS^2: ",MSN,MSNN
        WRITE(0,20)"(MS^2 is not an input parameter)"
       ELSE
        WRITE(0,30)" MS^2: ",MSN,MSNN
        WRITE(0,30)" XIS: ",XISN,XISNN
        WRITE(0,20)"(XIS is not an input parameter)"
       ENDIF
       IF(MUPN.NE.0d0 .OR. MUPNN.NE.0d0)
     . WRITE(0,30)" MU': ",MUPN,MUPNN
       IF(MSPN.NE.0d0 .OR. MSPNN.NE.0d0)
     . WRITE(0,30)" MS'^2: ",MSPN,MSPNN
       IF(DHN.NE.0d0 .OR. DHNN.NE.0d0)
     . WRITE(0,30)" DH: ",DHN,DHNN
       IF(LPPN.NE.0d0 .OR. LPPNN.NE.0d0)
     . WRITE(0,30)" LPP: ",LPPN,LPPNN
       IF(LTTN.NE.0d0 .OR. LTTNN.NE.0d0)
     . WRITE(0,30)" LTT: ",LTTN,LTTNN
       IF(LUN.NE.0d0 .OR. LUNN.NE.0d0)
     . WRITE(0,30)" LU: ",LUN,LUNN
       IF(LDN.NE.0d0 .OR. LDNN.NE.0d0)
     . WRITE(0,30)" LD: ",LDN,LDNN
       IF(LTN.NE.0d0 .OR. LTNN.NE.0d0)
     . WRITE(0,30)" LT: ",LTN,LBNN
       IF(LBN.NE.0d0 .OR. LBNN.NE.0d0)
     . WRITE(0,30)" LB: ",LBN,LBNN
       IF(LLN.NE.0d0 .OR. LLNN.NE.0d0)
     . WRITE(0,30)" LL: ",LLN,LLNN
       
      ENDIF

 20   FORMAT(A40,I10)
 30   FORMAT(A15,2E15.4)

      END
