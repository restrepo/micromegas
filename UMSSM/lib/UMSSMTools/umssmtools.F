        PROGRAM MAIN
************************************************************************
*
*  Phenomenological constraints:
*
*      PROB(I)  = 0    Not excluded by UMSSMTools
*            
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
*      PROB(28) =/= 0  unphysical global minimum
*      PROB(29) =/= 0  Higgs soft masses >> Msusy
*      PROB(32) =/= 0  b -> s gamma more than 2 sigma away
*      PROB(33) =/= 0  Delta M_s more than 2 sigma away
*      PROB(34) =/= 0  Delta M_d more than 2 sigma away
*      PROB(35) =/= 0  B_s -> mu+ mu- more than 2 sigma away
*      PROB(36) =/= 0  B+ -> tau nu_tau more than 2 sigma away
*      PROB(37) =/= 0  (g-2)_muon more than 2 sigma away
*      PROB(38) =/= 0  excluded by Upsilon(1S) -> A/H gamma
*      PROB(39) =/= 0  excluded by eta_b(1S) mass measurement
*      PROB(40) =/= 0  BR(B-->X_s mu+ mu-) more than 2 sigma away
*      PROB(41) =/= 0  excluded by ee -> hZ, h -> AA -> 4taus (new ALEPH analysis)
*      PROB(42) =/= 0  excluded by top -> b H+, H+ -> c s (CDF, D0)
*      PROB(43) =/= 0  excluded by top -> b H+, H+ -> tau nu_tau (D0)
*      PROB(44) =/= 0  excluded by top -> b H+, H+ -> W+ A, A -> 2taus (CDF)
*      PROB(45) =/= 0  excluded by t -> bH+ (LHC)"
*      PROB(46) =/= 0  no Higgs in the [MHMIN,MHMAX] GeV mass range
*      PROB(47) =/= 0  chi2(H->gg) > chi2MAX
*      PROB(48) =/= 0  chi2(H->bb) > chi2MAX
*      PROB(49) =/= 0  chi2(H->ZZ) > chi2MAX
*      PROB(51) =/= 0: excluded by H/A->tautau
*      PROB(52) =/= 0: excluded by H->AA->4leptons/2lept.+2b (LHC)
*      PROB(53) =/= 0: excluded by ggF->H/A->gamgam (65GeV < M < 122GeV, ATLAS)
*      PROB(55) =/= 0: b -> d gamma more than 2 sigma away
*      PROB(56) =/= 0: B_d -> mu+ mu- more than 2 sigma away
*      PROB(57) =/= 0: b -> s nu_L nubar_L more than 2 sigma away
*      PROB(58) =/= 0: b -> c tau nu_tau more than 2 sigma away (as SM)
*      PROB(59) =/= 0: K -> pi nu_L nubar_L more than 2 sigma away
*      PROB(60) =/= 0: Delta M_K / eps_K more than 2 sigma away
*
************************************************************************

        IMPLICIT NONE

        INTEGER NPROB
        PARAMETER (NPROB=60)
        INTEGER I

        DOUBLE PRECISION PROB(NPROB)

*   I/O files

        OPEN(15,FILE='UMSSM_inp.dat',STATUS='UNKNOWN')
        OPEN(16,FILE='SM_decay.dat', STATUS= 'UNKNOWN')
        OPEN(17,FILE='UMSSM_spectr.dat', STATUS= 'UNKNOWN')
        OPEN(18,FILE='UMSSM_decay.dat', STATUS= 'UNKNOWN')
	
*   Initialization, input parameters

        CALL INITIALIZE()
        CALL INPUT()

*   Initialization of PROB

        DO I=1,NPROB
         PROB(I)=0d0
        ENDDO      

*   Computation of parameters at QSTSB	
        CALL RUNPAR()

*   Computation of Higgs branching ratios in UMSSM and SM + t decays
        CALL UDECAY()
        CALL TDECAY()
        CALL SMDECAY()

*   Exp. constraints on sparticles/Higgses

        CALL OLDEXP(PROB)
        CALL LHCHIG(PROB)
        CALL Higgs_CHI2(PROB)
        CALL CMS_TAUTAU(PROB)
        CALL ATLAS_H_GAMGAM(PROB)
        CALL LHC_HSMAA_LEPTONS(PROB)

*   B + K physics

        CALL BPHYS(PROB)
        CALL KPHYS(PROB)

*   Anom. magn. moment of the Muon

        CALL MAGNMU(PROB)

*   Global minimum?

        CALL CHECKMIN(PROB)


*   Recording of the results
        CALL OUTPUT(PROB)

        END


        SUBROUTINE INPUT()

*******************************************************************
*   This subroutine reads SM and UMSSM parameters from input file   .
*******************************************************************

        IMPLICIT NONE

        INTEGER N0,NLOOP,NBER,Q2FIX,VFLAG,HFLAG

        DOUBLE PRECISION ACC,XITLA,XLAMBDA,MC0,MB0,MT0
        DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
        DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        DOUBLE PRECISION VUS,VCB,VUB,Q2,Q2MIN

        COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
        COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        COMMON/CKM/VUS,VCB,VUB
        COMMON/ALS/XLAMBDA,MC0,MB0,MT0,N0
        COMMON/RENSCALE/Q2
        COMMON/Q2FIX/Q2MIN,Q2FIX
        COMMON/FLAGS/VFLAG,HFLAG


*   DEFAULT VALUE FOR FLAGS
        HFLAG=0
        VFLAG=0

*   DEFAULT VALUE FOR THE RENSCALE Q2 (DEFAULT: 0.D0)
        Q2=0.D0


*   Set Q2MIN, Q2FIX:
        Q2MIN=100.D0**2
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
        ACC=1.D-8
        XLAMBDA=XITLA(NLOOP,ALSMZ,ACC)
        CALL ALSINI(ACC)
        CALL BERNINI(NBER)

        RETURN

        END

        SUBROUTINE OUTPUT(PROB)

*********************************************************************	
*   Subroutine writing all the results in the the output files.
*********************************************************************	
 
        IMPLICIT NONE

        INTEGER I,J,K,PDGNEU(6),PDGCHA(2),PDGST(2)
        DOUBLE PRECISION PROB(*)
        DOUBLE PRECISION SMASS(3),S(3,3),PMASS,P(3,3),CMASS
        DOUBLE PRECISION MGL,MCHA(2),U(2,2),V(2,2),MNEU(6),N(6,6)
        DOUBLE PRECISION BRJJ(4),BREE(4),BRMM(4),BRLL(4),BRSS(4),BRCC(4)
        DOUBLE PRECISION BRBB(4),BRTT(4),BRWW(3),BRZZ(3),BRGG(4)
        DOUBLE PRECISION BRZG(4),BRHHH(4),BRHAA(3),BRHCHC(3)
        DOUBLE PRECISION BRHAZ(3),BRAHZ(3),BRHCW(4)
        DOUBLE PRECISION BRHIGGS(4),BRNEU(4,6,6),BRCHA(4,3)
        DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,9),BRASQ(2),BRASL
        DOUBLE PRECISION BRSUSY(4),WIDTH(4)
        DOUBLE PRECISION HCBRE,HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC
        DOUBLE PRECISION HCBRBT,HCBRWH(4),HCBRWHT,HCBRNC(6,2)
        DOUBLE PRECISION HCBRSQ(5),HCBRSL(3),HCBRSUSY,HCWIDTH
        DOUBLE PRECISION CU(4),CD(4),CV(3),CVZ(3),CJ(4),CG(4),CB(4)
        DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
        DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        DOUBLE PRECISION VUS,VCB,VUB,TANB,SINB,COSB
        DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
        DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
        DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSNM,CSMU
        DOUBLE PRECISION SST,SSB,SSL,Q2
        DOUBLE PRECISION MH1X,MH2X,MSSX,LQ,KQ,ALQ,AKQ,MUQ,NUQ
        DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ,QSTSB,HSQ
        DOUBLE PRECISION delmagmu,damumin,damumax,amuthmax,amuthmin
        DOUBLE PRECISION UPARF,SAZZ,CAZZ,VEV,NCP,VEVS,G1P
        DOUBLE PRECISION QD,QU,QS,QQ,QUP,QDOW,QL,QE,QN
        DOUBLE PRECISION WIDTHSM(3),BRJJSM(3),BREESM(3),BRMMSM(3),
     .        BRLLSM(3)
        DOUBLE PRECISION BRSSSM(3),BRCCSM(3),BRBBSM(3),BRTTSM(3)
        DOUBLE PRECISION BRWWSM(3),BRZZSM(3),BRGGSM(3),BRZGSM(3)
        DOUBLE PRECISION brtopneutrstop(6,2),toptot,brtopbw,brtopbh
        DOUBLE PRECISION BRSG,BRSGmax,BRSGmin,DMd,DMdmin,DMdmax,DMs,
     .        DMsmax,DMsmin,BRBMUMU,BRBMUMUmax,BRBMUMUmin,BRBtaunu,
     .        BRBtaunumax,BRBtaunumin
        DOUBLE PRECISION BRBSll,BRBSllmin,BRBSllMax,
     .        BRBShll,BRBShllmin,BRBShllMax,
     .        BRDG,BRDGmin,BRDGmax,
     .        BRBdMUMU,BRBdMUMUmin,BRBdMUMUmax,
     .        BRBXsnunu,BRBXsnunumin,BRBXsnunumax,
     .        BRBpKpnunu,BRBpKpnunumin,BRBpKpnunuMax,
     .        BRBKsnunu,BRBKsnunumin,BRBKsnunuMax,
     .        RD_taul,RD_taulmin,RD_taulmax,
     .        RDs_taul,RDs_taulmin,RDs_taulmax
        DOUBLE PRECISION BRKpPipnunu,BRKpPipnunumin,BRKpPipnunuMax,
     .       BRKLPi0nunu,BRKLPi0nunumin,BRKLPi0nunuMax,
     .       DMK,DMKmin,DMKmax,epsK,epsKmin,epsKmax
        DOUBLE PRECISION eps0,epst0,epst1,epst2,epst3,epsts,epstb
        DOUBLE PRECISION epsY32,epsY31,epsY23,epsY13,epscs,epscb
        DOUBLE PRECISION SGNT,SGNB,SGNL,AU,AD,M2,LAMBDA,AL,ATAU
        DOUBLE PRECISION SIG(3,10),F2
        DOUBLE PRECISION MHmin,MHmax,chi2max,chi2gam,chi2bb,chi2zz

	
        COMMON/BRN/BRJJ,BREE,BRMM,BRLL,BRSS,BRCC,BRBB,BRTT,BRWW,BRZZ,
     .          BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHZ,
     .          BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     .          BRSUSY,WIDTH
        COMMON/BRC/HCBRE,HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC,
     .          HCBRBT,HCBRWH,HCBRWHT,HCBRNC,HCBRSQ,HCBRSL,
     .          HCBRSUSY,HCWIDTH
        COMMON/MAGMU/delmagmu,damumin,damumax,amuthmax,amuthmin
        COMMON/REDCOUP/CU,CD,CV,CVZ,CJ,CG,CB
        COMMON/HIGGSPEC/SMASS,S,PMASS,P,CMASS
        COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,N
        COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
        COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
        COMMON/CKM/VUS,VCB,VUB
        COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .          MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .          CST,CSB,CSL,MSMU1,MSMU2,MSNM,CSMU
        COMMON/RENSCALE/Q2
        COMMON/STSBSCALE/QSTSB
        COMMON/QMHIGGS/MH1X,MH2X,MSSX
        COMMON/QNMPAR/LQ,KQ,ALQ,AKQ,MUQ,NUQ
        COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ,HSQ
        COMMON/BR_top2body/brtopbw,brtopbh,brtopneutrstop
        COMMON/topwidth/toptot
        COMMON/UMSSM/SAZZ,CAZZ,VEV,NCP,QD,QU,QS,VEVS,G1P,QQ,
     .          QUP,QDOW,QL,QE,QN
        COMMON/BRSM/WIDTHSM,BRJJSM,BREESM,BRMMSM,BRLLSM,BRSSSM,BRCCSM,
     .          BRBBSM,BRTTSM,BRWWSM,BRZZSM,BRGGSM,BRZGSM
        COMMON/BOBS/BRSG,BRSGmax,BRSGmin,DMd,DMdmin,DMdmax,DMs,
     .          DMsmax,DMsmin,BRBMUMU,BRBMUMUmax,BRBMUMUmin,BRBtaunu,
     .          BRBtaunumax,BRBtaunumin,BRBSll,BRBSllmin,BRBSllMax,
     .          BRBShll,BRBShllmin,BRBShllMax,
     .          BRDG,BRDGmin,BRDGmax,
     .          BRBdMUMU,BRBdMUMUmin,BRBdMUMUmax,
     .          BRBXsnunu,BRBXsnunumin,BRBXsnunumax,
     .          BRBpKpnunu,BRBpKpnunumin,BRBpKpnunuMax,
     .          BRBKsnunu,BRBKsnunumin,BRBKsnunuMax,
     .          RD_taul,RD_taulmin,RD_taulmax,
     .          RDs_taul,RDs_taulmin,RDs_taulmax
        COMMON/KOBS/BRKpPipnunu,BRKpPipnunumin,BRKpPipnunuMax,
     .             BRKLPi0nunu,BRKLPi0nunumin,BRKLPi0nunuMax,
     .             DMK,DMKmin,DMKmax,epsK,epsKmin,epsKmax
        COMMON/EPSCOUP/eps0,epst0,epst1,epst2,epst3,epsts,epstb,
     .          epsY32,epsY31,epsY23,epsY13,epscs,epscb
        COMMON/LHCSIG/SIG
        COMMON/WARN/F2
        COMMON/HIGGSFIT/MHmin,MHmax,chi2max,chi2gam,chi2bb,chi2zz
        COMMON/NOBUG/TANB,AU,AD,M2,SST,SSB,LAMBDA,AL,ATAU

        COSB=1.D0/DSQRT(1.D0+TANB**2)
        SINB=TANB*COSB

        WRITE(17,899) "# UMSSMTools OUTPUT IN SLHA FORMAT"
        WRITE(17,899) "# Info about spectrum calculator"
        WRITE(17,899) "BLOCK SPINFO        # Program information"
        WRITE(17,900) 1,"UMSSMTools # spectrum calculator"
        WRITE(17,900) 2,"1.1.2      # version number"

        WRITE(17,900) 3,"WARNINGS FROM UMSSMTools :"
        WRITE(17,900)
     . 3,"Flags 2-23,32-40,42-49,51-53,55-60; else 0=OK"
        DO I=1,60
         IF(PROB(I).NE.0d0) GOTO 5
        ENDDO   
        WRITE(17,900) 3," 0  # Not excluded by UMSSMTools"
 5    IF(PROB(2).NE.0d0)
     . WRITE(17,900) 3," 2  # Neutralinos too light"
      IF(PROB(3).NE.0d0)
     . WRITE(17,900) 3," 3  # Charged Higgs too light"
      IF(PROB(4).NE.0d0)
     . WRITE(17,900) 3," 4  # Excluded by ee -> hZ, ind. of h decay"
      IF(PROB(5).NE.0d0)
     . WRITE(17,900) 3," 5  # Excluded by ee -> hZ, h -> bb"
      IF(PROB(6).NE.0d0)
     . WRITE(17,900) 3," 6  # Excluded by ee -> hZ, h -> tautau"
      IF(PROB(7).NE.0d0)
     . WRITE(17,900) 3," 7  # Excluded by ee -> hZ, h -> invisible"
      IF(PROB(8).NE.0d0)
     . WRITE(17,900) 3," 8  # Excluded by ee -> hZ, h -> 2jets"
      IF(PROB(9).NE.0d0)
     . WRITE(17,900) 3," 9  # Excluded by ee -> hZ, h -> 2photons"
      IF(PROB(10).NE.0d0)
     . WRITE(17,900) 3,"10  # Excluded by ee -> hZ, h -> AA -> 4bs"
      IF(PROB(11).NE.0d0 .OR. PROB(41).NE.0d0)
     . WRITE(17,900) 3,"11  # Excluded by ee -> hZ, h -> AA -> 4taus"
      IF(PROB(12).NE.0d0)
     . WRITE(17,900) 
     . 3,"12  # Excluded by ee -> hZ, h -> AA -> 2bs 2taus"
      IF(PROB(13).NE.0d0)
     . WRITE(17,900) 3,"13  # Excluded by ee -> Z -> hA (Z width)"
      IF(PROB(14).NE.0d0)
     . WRITE(17,900) 3,"14  # Excluded by ee -> hA -> 4bs"
      IF(PROB(15).NE.0d0)
     . WRITE(17,900) 3,"15  # Excluded by ee -> hA -> 4taus"
      IF(PROB(16).NE.0d0)
     . WRITE(17,900) 3,"16  # Excluded by ee -> hA -> 2bs 2taus"
      IF(PROB(17).NE.0d0)
     . WRITE(17,900) 3,"17  # Excluded by ee -> hA -> AAA -> 6bs"
      IF(PROB(18).NE.0d0)
     . WRITE(17,900) 3,"18  # Excluded by ee -> hA -> AAA -> 6taus"
      IF(PROB(19).NE.0d0)
     . WRITE(17,900)
     . 3,"19  # Excluded by ee -> hZ, h -> AA,A -> light pair"
      IF(PROB(20).NE.0d0)
     . WRITE(17,900) 3,"20  # Excluded by stop -> b l sneutrino"
      IF(PROB(21).NE.0d0)
     . WRITE(17,900) 3,"21  # Excluded by stop -> neutralino c"
      IF(PROB(22).NE.0d0)
     . WRITE(17,900) 3,"22  # Excluded by sbottom -> neutralino b"
      IF(PROB(23).NE.0d0)
     . WRITE(17,900) 3,"23  # Squark/gluino too light"
      IF(PROB(28).NE.0d0)
     . WRITE(17,900) 3,"28  # Unphysical global minimum"
      IF(PROB(29).NE.0d0)
     . WRITE(17,900) 3,"29  # Higgs soft masses >> Msusy"
      IF(PROB(32).NE.0d0)
     . WRITE(17,900) 3,"32  # b -> s gamma more than 2 sigma away"
      IF(PROB(33).NE.0d0)
     . WRITE(17,900) 3,"33  # Delta M_s more than 2 sigma away"
      IF(PROB(34).NE.0d0)
     . WRITE(17,900) 3,"34  # Delta M_d more than 2 sigma away"
      IF(PROB(35).NE.0d0)
     . WRITE(17,900) 3,"35  # B_s -> mu+ mu- more than 2 sigma away"
      IF(PROB(36).NE.0d0)
     . WRITE(17,900) 3,"36  # B+ -> tau nu_tau more than 2 sigma away"
      IF(PROB(37).NE.0d0)
     . WRITE(17,900) 3,"37  # Muon magn. mom. more than 2 sigma away"
      IF(PROB(38).NE.0d0)
     . WRITE(17,900) 
     . 3,"38  # Excluded by Upsilon(1S) -> A gamma (CLEO)"
      IF(PROB(38).LT.0d0)
     . WRITE(17,900) 3,"38  # (but A width> 10 MeV)"
      IF(PROB(39).NE.0d0)
     . WRITE(17,900)
     . 3,"39  # Excluded etab(1S) mass difference (BABAR-theory)"
       IF(PROB(40).NE.0d0)
     . WRITE(17,900) 3,"40  # Excluded by BR(B -> X_s mu +mu-)"
       IF(PROB(42).NE.0d0)
     . WRITE(17,900) 3,"42  # Excluded by top -> b H+, H+ -> c s"
       IF(PROB(43).NE.0d0)
     . WRITE(17,900) 
     . 3,"43  # Excluded by top -> b H+, H+ -> tau nu_tau"
       IF(PROB(44).NE.0d0)
     . WRITE(17,900)
     . 3,"44  # Excluded by top -> b H+, H+ -> W+ A, A -> 2taus"
       IF(PROB(45).NE.0d0)
     . WRITE(17,900) 3,"45  # Excluded by t -> bH+ (LHC)"
       IF(PROB(46).NE.0d0)
     . WRITE(17,918)
     . 3,"46  # No Higgs in the",MHMIN,MHMAX," GeV mass range"
       IF(PROB(47).NE.0d0)
     . WRITE(17,919) 3,"47  # chi2(H->gg) >",chi2MAX
       IF(PROB(48).NE.0d0)
     . WRITE(17,919) 3,"48  # chi2(H->bb) >",chi2MAX
       IF(PROB(49).NE.0d0)
     . WRITE(17,919) 3,"49  # chi2(H->ZZ) >",chi2MAX
       IF(PROB(51).NE.0d0)
     . WRITE(17,900)
     . 3,"51  # Excluded by ggF/bb->H/A->tautau at the LHC"
       IF(PROB(52).NE.0d0)
     . WRITE(17,900)
     . 3,"52  # Excluded by H->AA->4leptons/2lept.+2b (LHC)"
       IF(PROB(53).NE.0d0)
     . WRITE(17,900) 3,"53  # Excluded by ggF->H/A->gamgam (ATLAS)"
      IF(PROB(55).NE.0d0)
     . WRITE(17,900) 3,"55  # b -> d gamma more than 2 sigma away"
      IF(PROB(56).NE.0d0)
     . WRITE(17,900) 3,"56  # B_d -> mu+ mu- more than 2 sigma away"
      IF(PROB(57).NE.0d0)
     . WRITE(17,900)
     . 3,"57  # b -> s nu_L nubar_L more than 2 sigma away"
      IF(PROB(58).NE.0d0)
     . WRITE(17,900)
     . 3,"58  # b -> c tau nu_tau more than 2 sigma away (as SM)"
      IF(PROB(59).NE.0d0)
     . WRITE(17,900) 3,"59  # K -> pi nu nubar more than 2 sigma away"
      IF(PROB(60).NE.0d0)
     . WRITE(17,900) 3,"60  # DMK / epsK more than 2 sigma away"


        PDGNEU(1)=1000022
        PDGNEU(2)=1000023
        PDGNEU(3)=1000025
        PDGNEU(4)=1000035
        PDGNEU(5)=1000045
        PDGNEU(6)=1000032
        PDGCHA(1)=1000024
        PDGCHA(2)=1000037
        PDGST(1)=1000006
        PDGST(2)=2000006

 
        WRITE(17,899) "# Input parameters"
        WRITE(17,899) "BLOCK MODSEL"
        WRITE(17,899) "  3  3  # UMSSM PARTICLE CONTENT"
        WRITE(17,899) "BLOCK SMINPUTS"
        WRITE(17,901) 1,1.D0/ALEMMZ,"ALPHA_EM^-1(MZ)"
        WRITE(17,901) 2,GF,"GF"
        WRITE(17,901) 3,ALSMZ,"ALPHA_S(MZ)"
        WRITE(17,901) 4,MZ,"MZ1"
        WRITE(17,901) 5,MB,"MB"
        WRITE(17,901) 6,MT,"MTOP"
        WRITE(17,901) 7,MTAU,"MTAU"
        WRITE(17,899) "# SMINPUTS Beyond SLHA:"
        WRITE(17,906) "MW:",MW
        WRITE(17,906) "MS:",MS
        WRITE(17,906) "MC:",MC
        WRITE(17,906) "VUS:",VUS
        WRITE(17,906) "VCB:",VCB
        WRITE(17,906) "VUB:",VUB

        WRITE(17,899) "BLOCK MINPAR"
        WRITE(17,901) 3,TANB,"TANBETA"

        WRITE(17,899) "BLOCK EXTPAR"
        WRITE(17,901) 1,UPARF(0),"M1"
        WRITE(17,901) 2,M2,"M2"
        WRITE(17,901) 3,UPARF(2),"M3"
        WRITE(17,901) 11,AU,"ATOP"
        WRITE(17,901) 12,AD,"ABOTTOM"
        WRITE(17,901) 13,ATAU,"ATAU"
        WRITE(17,901) 16,UPARF(20),"AMUON"

        WRITE(17,901) 31,UPARF(5),"LEFT SELECTRON"
        WRITE(17,901) 32,UPARF(5),"LEFT SMUON"
        WRITE(17,901) 33,UPARF(6),"LEFT STAU"

        WRITE(17,901) 34,UPARF(7),"RIGHT SELECTRON"
        WRITE(17,901) 35,UPARF(7),"RIGHT SMUON"
        WRITE(17,901) 36,UPARF(8),"RIGHT STAU"

        WRITE(17,901) 41,UPARF(9),"LEFT 1ST GEN. SQUARKS"
        WRITE(17,901) 42,UPARF(9),"LEFT 2ND GEN. SQUARKS"
        WRITE(17,901) 43,UPARF(10),"LEFT 3RD GEN. SQUARKS"

        WRITE(17,901) 44,UPARF(11),"RIGHT U-SQUARKS"
        WRITE(17,901) 45,UPARF(11),"RIGHT C-SQUARKS"
        WRITE(17,901) 46,UPARF(12),"RIGHT T-SQUARKS"

        WRITE(17,901) 47,UPARF(13),"RIGHT D-SQUARKS"
        WRITE(17,901) 48,UPARF(13),"RIGHT S-SQUARKS"
        WRITE(17,901) 49,UPARF(14),"RIGHT B-SQUARKS"

        WRITE(17,901) 50,UPARF(15),"RIGHT SNEUTRINO E"
        WRITE(17,901) 51,UPARF(15),"RIGHT SNEUTRINO MU"
          IF(UPARF(16).LT.0.D0)THEN
            WRITE(17,901) 52,DSQRT(-UPARF(16)),"*I RIGHT SNEUTRINO TAU"
          ELSE
            WRITE(17,901) 52,DSQRT(UPARF(16)),"RIGHT SNEUTRINO TAU"     
          ENDIF   

        WRITE(17,901) 61,LAMBDA,"LAMBDA"
        WRITE(17,901) 63,AL,"A_LAMBDA"
        WRITE(17,901) 65,MUQ,"MU"

        WRITE(17,901) 102,UPARF(40),"THETA_E6"
        WRITE(17,901) 103,UPARF(4),"BINO-BINO' MIXING"
        WRITE(17,901) 104,UPARF(3),"BINO'"
        WRITE(17,901) 105,UPARF(35),"MIXING Z-Z'"

        WRITE(17,899) "# "
        WRITE(17,899) "BLOCK MASS   # Mass spectrum "
        WRITE(17,899) "#  PDG Code     mass             particle "
        WRITE(17,902) 25,SMASS(1),"lightest neutral scalar"
        WRITE(17,902) 35,SMASS(2),"second neutral scalar"
        WRITE(17,902) 45,SMASS(3),"third neutral scalar"
        WRITE(17,902) 36,PMASS,"pseudoscalar"
        WRITE(17,902) 37,CMASS,"charged Higgs"
        WRITE(17,902) 32,UPARF(34),"Z2"

        WRITE(17,902) 1000001,MDL," ~d_L"
        WRITE(17,902) 2000001,MDR," ~d_R"
        WRITE(17,902) 1000002,MUL," ~u_L"
        WRITE(17,902) 2000002,MUR," ~u_R"
        WRITE(17,902) 1000003,MDL," ~s_L"
        WRITE(17,902) 2000003,MDR," ~s_R"
        WRITE(17,902) 1000004,MUL," ~c_L"
        WRITE(17,902) 2000004,MUR," ~c_R"
        WRITE(17,902) 1000005,MSB1," ~b_1"
        WRITE(17,902) 2000005,MSB2," ~b_2"
        WRITE(17,902) 1000006,MST1," ~t_1"
        WRITE(17,902) 2000006,MST2," ~t_2"
        WRITE(17,902) 1000011,MLL," ~e_L"
        WRITE(17,902) 2000011,MLR," ~e_R"
        WRITE(17,902) 1000012,MNL," ~nue_L"
        WRITE(17,902) 1000013,UPARF(93)," ~mu_L"
        WRITE(17,902) 2000013,UPARF(94)," ~mu_R"
        WRITE(17,902) 1000014,MSNM," ~numu_L"
        WRITE(17,902) 1000015,MSL1," ~tau_1"
        WRITE(17,902) 2000015,MSL2," ~tau_2"
        WRITE(17,902) 1000016,MSNT," ~nutau_L"
        WRITE(17,902) 2000012,UPARF(102)," ~nue_R"
        WRITE(17,902) 2000014,UPARF(103)," ~numu_R"
        WRITE(17,902) 2000016,UPARF(104)," ~nutau_R"

        WRITE(17,902) 1000021,MGL," ~g"

      DO I=1,6
        WRITE(17,922) PDGNEU(I),MNEU(I),"neutralino(",I,")"
      ENDDO  
        WRITE(17,902) 1000024,MCHA(1),"chargino(1)"
        WRITE(17,902) 1000037,MCHA(2),"chargino(2)"

        WRITE(17,899) "# "
        WRITE(17,899) "# Low energy observables"
        WRITE(17,899) "BLOCK LOWEN"
        WRITE(17,899)
     .   "# Exp. 2 Sigma: 2.99E-4 < BR(b -> s gamma) < 3.87E-4:"
        WRITE(17,901) 1,BRSG,"BR(b -> s gamma)"
        WRITE(17,901) 11,BRSGMAX,"BR(b -> s gamma) + Theor.Err.)"
        WRITE(17,901) 12,BRSGMIN,"BR(b -> s gamma) - Theor.Err.)"
        WRITE(17,899) "# Exp. 2 Sigma: 5.04E-1 < Delta M_d < 5.16E-1:"
        WRITE(17,901) 2,DMD,"Delta M_d in ps^-1"
        WRITE(17,901) 21,DMdmax,"Delta M_d + Theor.Err."
        WRITE(17,901) 22,DMdmin,"Delta M_d - Theor.Err."
        WRITE(17,899) 
     .   "# Exp. 2 Sigma: 1.7715E+1 < Delta Ms < 1.7799E+1:"
        WRITE(17,901) 3,DMS,"Delta M_s in ps^-1"
        WRITE(17,901) 31,DMsmax,"Delta M_s + Theor.Err."
        WRITE(17,901) 32,DMsmin,"Delta M_s - Theor.Err."
        WRITE(17,899) 
     .   "# Exp. 2 Sigma: 1.7E-09 < BR(Bs -> mu+ mu-) < 4.5E-09:"
        WRITE(17,901) 4,BRBMUMU,"BR(Bs -> mu+ mu-)"
        WRITE(17,901) 41,BRBMUMUmax,"BR(Bs -> mu+ mu-) + Theor.Err."
        WRITE(17,901) 42,BRBMUMUmin,"BR(Bs -> mu+ mu-) - Theor.Err."
        WRITE(17,899) 
     .   "# Exp. 2 Sigma: 0.70E-4 < BR(B+ > tau+ + nu_tau) < 1.58E-4:"
        WRITE(17,901) 5,BRBtaunu,"BR(B+ -> tau+ nu_tau)"
        WRITE(17,901) 51,BRBtaunumax,
     .   "BR(B+ -> tau+ nu_tau) + Theor.Err."
        WRITE(17,901) 52,BRBtaunumin,
     .   "BR(B+ -> tau+ nu_tau) - Theor.Err."
        WRITE(17,899) "# BSM contr. to the muon anomalous magn. moment:"
        WRITE(17,901) 6,delmagmu,"Del_a_mu"
        WRITE(17,901) 61,amuthmax,"Del_a_mu + Theor.Err."
        WRITE(17,901) 62,amuthmin,"Del_a_mu - Theor.Err."
        WRITE(17,907) "# Minimal Exp. - SM (2 sigma):",damumin
        WRITE(17,907) "# Maximal Exp. - SM (2 sigma):",damumax
        WRITE(17,899) 
     .   "# Exp. 2 Sigma: 0.84E-6 < BR(B -> Xs l+ l-)low < 2.32E-6:"
        WRITE(17,899) "#(for 1 GeV^2 < M_{l+l-}^2 < 6 GeV^2)"
        WRITE(17,901) 7,BRBSll,"BR(B -> X_s l+ l-)low"
        WRITE(17,901)
     .   71,BRBSllmax,"BR(B -> X_s l+ l-)low + Theor.Err."
        WRITE(17,901)
     .   72,BRBSllmin,"BR(B -> X_s l+ l-)low - Theor.Err."
        WRITE(17,899) 
     .   "# Exp. 2 Sigma: 2.8E-7 < BR(B -> Xs l+ l-)high < 6.8E-7:"
        WRITE(17,899) "#(for 14.4 GeV^2 < M_{l+l-}^2 < mb^2)"
        WRITE(17,901) 8,BRBShll,"BR(B -> X_s l+ l-)high"
        WRITE(17,901)
     .   81,BRBShllmax,"BR(B -> X_s l+ l-)high + Theor.Err."
        WRITE(17,901)
     .   82,BRBShllmin,"BR(B -> X_s l+ l-)high - Theor.Err."
        WRITE(17,899)
     .   "# Exp. 2 Sigma: 2.7E-6 < BR(b -> d gamma) < 2.55E-5:"
        WRITE(17,901) 9,BRDG,"BR(b -> d gamma)"
        WRITE(17,901) 91,BRDGMAX,"BR(b -> d gamma) + Theor.Err.)"
        WRITE(17,901) 92,BRDGMIN,"BR(b -> d gamma) - Theor.Err.)"
        WRITE(17,899) 
     .   "# Exp. 2 Sigma: 1.1E-10 < BR(Bd -> mu+ mu-) < 7.1E-10:"
        WRITE(17,901) 10,BRBdMUMU,"BR(Bd -> mu+ mu-)"
        WRITE(17,901) 101,BRBdMUMUmax,"BR(Bd -> mu+ mu-) + Theor.Err."
        WRITE(17,901) 102,BRBdMUMUmin,"BR(Bd -> mu+ mu-) - Theor.Err."
        WRITE(17,899) 
     .   "# Exp. 2 Sigma: BR(B -> Xs nu_L nubar_L) < 6.4E-4:"
        WRITE(17,901) 110,BRBXsnunu,"BR(B-> Xs nu_L nubar_L)"
        WRITE(17,901)
     .   111,BRBXsnunumax,"BR(B -> Xs nu_L nubar_L) + Theor.Err."
        WRITE(17,901)
     .   112,BRBXsnunumin,"BR(B -> Xs nu_L nubar_L) - Theor.Err."
        WRITE(17,899) 
     .   "# Exp. 2 Sigma: BR(B+ -> K+ nu_L nubar_L) < 1.6E-5:"
        WRITE(17,901) 120,BRBpKpnunu,"   BR(B+-> K+ nu_L nubar_L)"
        WRITE(17,901)
     .   121,BRBpKpnunumax,"BR(B+ -> K+ nu_L nubar_L) + Theor.Err."
        WRITE(17,901)
     .   122,BRBpKpnunumin,"BR(B+ -> K+ nu_L nubar_L) - Theor.Err."
        WRITE(17,899) 
     .   "# Exp. 2 Sigma: BR(B -> Ks nu_L nubar_L) < 5.5E-5:"
        WRITE(17,901) 13,BRBKsnunu,"   BR(B-> Ks nu_L nubar_L)"
        WRITE(17,901)
     .   131,BRBKsnunumax,"BR(B -> Ks nu_L nubar_L) + Theor.Err."
        WRITE(17,901)
     .   132,BRBKsnunumin,"BR(B -> Ks nu_L nubar_L) - Theor.Err."
        WRITE(17,899)
     .   "# RD = BR[B+ -> D tau+ nu_tau]/BR[B+ -> D l+ nu_l]"
        WRITE(17,899) 
     .   "# Exp. 2 Sigma: 2.91E-1 < RD < 4.91E-1:"
        WRITE(17,901) 14,RD_taul,"RD"
        WRITE(17,901) 141,RD_taulmax,"RD + Theor.Err."
        WRITE(17,901) 142,RD_taulmin,"RD - Theor.Err."
        WRITE(17,899)
     .   "# RD* = BR[B+ -> D* tau+ nu_tau]/BR[B+ -> D* l+ nu_l]"
        WRITE(17,899) 
     .   "# Exp. 2 Sigma: 2.76E-1 < RD* < 3.68E-1:"
        WRITE(17,901) 15,RDs_taul,"RD*"
        WRITE(17,901) 151,RDs_taulmax,"RD* + Theor.Err."
        WRITE(17,901) 152,RDs_taulmin,"RD* - Theor.Err."
        WRITE(17,899) 
     .   "# Exp. 2 Sigma: BR(K+ -> Pi+ nu_L nubar_L) < 4.03E-10:"
        WRITE(17,901) 16,BRKpPipnunu,"BR(K+ -> Pi+ nu_L nubar_L)"
        WRITE(17,901)
     .   161,BRKpPipnunumax,"BR(K+ -> Pi+ nu_L nubar_L) + Theor.Err."
        WRITE(17,901)
     .   162,BRKpPipnunumin,"BR(K+ -> Pi+ nu_L nubar_L) - Theor.Err."
        WRITE(17,899) 
     .   "# Exp. 90%CL: BR(KL -> Pi0 nu_L nubar_L) < 2.6E-8:"
        WRITE(17,901) 17,BRKLPi0nunu,"BR(KL -> Pi0 nu_L nubar_L)"
        WRITE(17,901)
     .   171,BRKLPi0nunumax,"BR(KL -> Pi0 nu_L nubar_L) + Theor.Err."
        WRITE(17,901)
     .   172,BRKLPi0nunumin,"BR(KL -> Pi0 nu_L nubar_L) - Theor.Err."
        WRITE(17,899) 
     .   "# Exp. 2 Sigma: 5.275E-3 < Delta M_K < 5.311E-3:"
        WRITE(17,901) 18,DMK,"Delta M_K in ps^-1"
        WRITE(17,901)
     .   181,DMKmax,"Delta M_K + Theor.Err."
        WRITE(17,901)
     .   182,DMKmin,"Delta M_K - Theor.Err."
        WRITE(17,899) 
     .   "# Exp. 2 Sigma: 2.206E-3 < eps_K < 2.250E-3:"
        WRITE(17,901) 19,epsK,"eps_K"
        WRITE(17,901)
     .   191,epsKmax,"eps_K + Theor.Err."
        WRITE(17,901)
     .   192,epsKmin,"eps_K - Theor.Err."


        WRITE(17,899) "# "
        WRITE(17,899) "# 3*3 Higgs mixing"
        WRITE(17,899) "BLOCK NMHMIX"
      DO I=1,3
      DO J=1,3
        WRITE(17,921) I,J,S(I,J),"S_(",I,",",J,")"
      ENDDO  
      ENDDO

        WRITE(17,899) "# "
        WRITE(17,899) "# 3*3 Pseudoscalar Higgs mixing"
        WRITE(17,899) "BLOCK NMAMIX"
      DO I=1,3
      DO J=1,3
        WRITE(17,921) I,J,P(I,J),"P_(",I,",",J,")"
      ENDDO  
      ENDDO

        SGNT= 1.D0
         IF(CST.GT.SST)THEN
          SGNT= -1.D0   
         ENDIF
        SGNB= 1.D0
         IF(CSB.GT.SSB)THEN
          SGNB= -1.D0   
         ENDIF
        SSL= UPARF(96)
        SGNL= 1.D0
         IF(CSL.GT.SSL)THEN
          SGNL= -1.D0   
         ENDIF

        WRITE(17,899) "# "
        WRITE(17,899) "# 3rd generation sfermion mixing"
        WRITE(17,899) "BLOCK STOPMIX  # Stop mixing matrix"
        WRITE(17,903) 1,1,CST,"Rst_(1,1)"
        WRITE(17,903) 1,2,SST,"Rst_(1,2)"
        WRITE(17,903) 2,1,SST*SGNT,"Rst_(2,1)"
        WRITE(17,903) 2,2,-CST*SGNT,"Rst_(2,2)"
        WRITE(17,899) "BLOCK SBOTMIX  # Sbottom mixing matrix"
        WRITE(17,903) 1,1,CSB,"Rsb_(1,1)"
        WRITE(17,903) 1,2,SSB,"Rsb_(1,2)"
        WRITE(17,903) 2,1,SSB*SGNB,"Rsb_(2,1)"
        WRITE(17,903) 2,2,-CSB*SGNB,"Rsb_(2,2)"
        WRITE(17,899) "BLOCK STAUMIX  # Stau mixing matrix"
        WRITE(17,903) 1,1,CSL,"Rsl_(1,1)"
        WRITE(17,903) 1,2,SSL,"Rsl_(1,2)"
        WRITE(17,903) 2,1,SSL*SGNL,"Rsl_(2,1)"
        WRITE(17,903) 2,2,-CSL*SGNL,"Rsl_(2,2)"

        WRITE(17,899) "# "
        WRITE(17,899) "# Gaugino-Higgsino mixing"
        WRITE(17,899) "BLOCK NMNMIX  # 6*6 Neutralino Mixing Matrix"
      DO I=1,6
      DO J=1,6
        WRITE(17,921) I,J,N(I,J),"N_(",I,",",J,")"
      ENDDO  
      ENDDO

        WRITE(17,899) "# "
        WRITE(17,899) "BLOCK UMIX  # Chargino U Mixing Matrix"
        WRITE(17,903) 1,1,U(1,1),"U_(1,1)"
        WRITE(17,903) 1,2,U(1,2),"U_(1,2)"
        WRITE(17,903) 2,1,U(2,1),"U_(2,1)"
        WRITE(17,903) 2,2,U(2,2),"U_(2,2)"

        WRITE(17,899) "# "
        WRITE(17,899) "BLOCK VMIX  # Chargino V Mixing Matrix"
        WRITE(17,903) 1,1,V(1,1),"V_(1,1)"
        WRITE(17,903) 1,2,V(1,2),"V_(1,2)"
        WRITE(17,903) 2,1,V(2,1),"V_(2,1)"
        WRITE(17,903) 2,2,V(2,2),"V_(2,2)"

        WRITE(17,899) "# "
        WRITE(17,899) "# REDUCED CROSS SECTIONS AT LHC"
      DO I=1,3
        WRITE(17,922) 11+(I-1)*10,SIG(I,1),"VBF/VH -> H",I," -> tautau"
        WRITE(17,922) 12+(I-1)*10,SIG(I,2),"ggF -> H",I," -> tautau"
        WRITE(17,922) 13+(I-1)*10,SIG(I,3),"VBF/VH -> H",I," -> bb"
        WRITE(17,922) 14+(I-1)*10,SIG(I,4),"ttH -> H",I," -> bb"
        WRITE(17,922) 15+(I-1)*10,SIG(I,5),"VBF/VH -> H",I," -> ZZ/WW"
        WRITE(17,922) 16+(I-1)*10,SIG(I,6),"ggF -> H",I," -> ZZ/WW"
        WRITE(17,922)
     .   17+(I-1)*10,SIG(I,7),"VBF/VH -> H",I," -> gammagamma"
        WRITE(17,922)
     .   18+(I-1)*10,SIG(I,8),"ggF -> H",I," -> gammagamma"
        WRITE(17,922)
     .   19+(I-1)*10,SIG(I,9),"VBF/VH -> H",I," -> invisible"
        WRITE(17,922)
     .   110+(I-1)*100,SIG(I,10),"ggF -> H",I," -> invisible"
      ENDDO

        WRITE(18,899) "# HIGGS + TOP BRANCHING RATIOS IN SLHA FORMAT"
        WRITE(18,899) "# Info about decay package"
        WRITE(18,899) "BLOCK DCINFO        # Program information"
        WRITE(18,900) 1,"UMSSMTools # Decay package"
        WRITE(18,900) 2,"1.1.2      # Version number"

        WRITE(18,899) "#           PDG          Width"
      DO I=1,3
        IF(I.EQ.1)
     .  WRITE(18,904) 25,WIDTH(I),"Lightest neutral Higgs scalar"
        IF(I.EQ.2)
     .  WRITE(18,904) 35,WIDTH(I),"2nd neutral Higgs scalar"
        IF(I.EQ.3)
     .  WRITE(18,904) 45,WIDTH(I),"3rd neutral Higgs scalar"
        WRITE(18,923) BRJJ(I),2,21,21,"BR(H_",I," -> gluon gluon)"
        WRITE(18,923) BREE(I),2,11,-11,"BR(H_",I," -> e e)"
        WRITE(18,923) BRMM(I),2,13,-13,"BR(H_",I," -> muon muon)"
        WRITE(18,923) BRLL(I),2,15,-15,"BR(H_",I," -> tau tau)"
        WRITE(18,923) BRSS(I),2,3,-3,"BR(H_",I," -> s sbar)"
        WRITE(18,923) BRCC(I),2,4,-4,"BR(H_",I," -> c cbar)"
        WRITE(18,923) BRBB(I),2,5,-5,"BR(H_",I," -> b bbar)"
        WRITE(18,923) BRTT(I),2,6,-6,"BR(H_",I," -> t tbar)"
        WRITE(18,923) BRWW(I),2,24,-24,"BR(H_",I," -> W+ W-)"
        WRITE(18,923) BRZZ(I),2,23,23,"BR(H_",I," -> Z1 Z1)"
        WRITE(18,923) BRGG(I),2,22,22,"BR(H_",I," -> gamma gamma)"
        WRITE(18,923) BRZG(I),2,23,22,"BR(H_",I," -> Z1 gamma)"
        IF(I.GE.2)THEN
          IF(BRHHH(I-1).GT.0.D0)
     .    WRITE(18,923) BRHHH(I-1),2,25,25,"BR(H_",I," -> H_1 H_1)"
        ENDIF
        IF(I.EQ.3.AND.BRHHH(3).GT.0.D0)
     .  WRITE(18,923) BRHHH(3),2,25,35,"BR(H_",I," -> H_1 H_2)"
        IF(I.EQ.3.AND.BRHHH(4).GT.0.D0)
     .  WRITE(18,923) BRHHH(4),2,35,35,"BR(H_",I," -> H_2 H_2)"
        IF(BRHAA(I).GT.0.D0)
     .  WRITE(18,923) BRHAA(I),2,36,36,"BR(H_",I," -> H_A H_A)"
        IF(BRHAZ(I).GT.0.D0)
     .  WRITE(18,923) BRHAZ(I),2,36,23,"BR(H_",I," -> H_A Z1)"
        IF(BRHCHC(I).GT.0.D0)
     .  WRITE(18,923) BRHCHC(I),2,37,-37,"BR(H_",I," -> H+ H-)"
        IF(BRHCW(I).GT.0.D0)
     .  WRITE(18,923) BRHCW(I),2,24,-37,"BR(H_",I," -> W+ H-)"
      DO J=1,6
      DO K=J,6
        IF(BRNEU(I,J,K).GT.0.D0)
     .  WRITE(18,924) BRNEU(I,J,K),2,PDGNEU(J),PDGNEU(K),
     .  "BR(H_",I," -> neu_",J," neu_",K,")"
      ENDDO
      ENDDO
        IF(BRCHA(I,1).GT.0.D0)
     .  WRITE(18,923) BRCHA(I,1),2,1000024,-1000024,
     .  "BR(H_",I," -> cha_1 cha_1bar)"
        IF(BRCHA(I,2).GT.0.D0)
     .  WRITE(18,923) BRCHA(I,2),2,1000024,-1000037,
     .  "BR(H_",I," -> cha_1 cha_2bar)"
        IF(BRCHA(I,2).GT.0.D0)
     .  WRITE(18,923) BRCHA(I,2),2,1000037,-1000024,
     .  "BR(H_",I," -> cha_2 cha_1bar)"
        IF(BRCHA(I,3).GT.0.D0)
     .  WRITE(18,923) BRCHA(I,3),2,1000037,-1000037,
     .  "BR(H_",I," -> cha_2 cha_2bar)"
        IF(BRHSQ(I,1).GT.0.D0)
     .  WRITE(18,923) BRHSQ(I,1),2,1000002,-1000002,
     .  "BR(H_",I," -> ~u_L ~ubar_L)"
        IF(BRHSQ(I,1).GT.0.D0)
     .  WRITE(18,923) BRHSQ(I,1),2,1000004,-1000004,
     .  "BR(H_",I," -> ~c_L ~cbar_L)"
        IF(BRHSQ(I,2).GT.0.D0)
     .  WRITE(18,923) BRHSQ(I,2),2,2000002,-2000002,
     .  "BR(H_",I," -> ~u_R ~ubar_R)"
        IF(BRHSQ(I,2).GT.0.D0)
     .  WRITE(18,923) BRHSQ(I,2),2,2000004,-2000004,
     .  "BR(H_",I," -> ~c_R ~cbar_R)"
        IF(BRHSQ(I,3).GT.0.D0)
     .  WRITE(18,923) BRHSQ(I,3),2,1000001,-1000001,
     .  "BR(H_",I," -> ~d_L ~dbar_L)"
        IF(BRHSQ(I,3).GT.0.D0)
     .  WRITE(18,923) BRHSQ(I,3),2,1000003,-1000003,
     .  "BR(H_",I," -> ~s_L ~sbar_L)"
        IF(BRHSQ(I,4).GT.0.D0)
     .  WRITE(18,923) BRHSQ(I,4),2,2000001,-2000001,
     .  "BR(H_",I," -> ~d_R ~dbar_R)"
        IF(BRHSQ(I,4).GT.0.D0)
     .  WRITE(18,923) BRHSQ(I,4),2,2000003,-2000003,
     .  "BR(H_",I," -> ~s_R ~sbar_R)"
        IF(BRHSQ(I,5).GT.0.D0)
     .  WRITE(18,923) BRHSQ(I,5),2,1000006,-1000006,
     .  "BR(H_",I," -> ~t_1 ~tbar_1)"
        IF(BRHSQ(I,6).GT.0.D0)
     .  WRITE(18,923) BRHSQ(I,6),2,2000006,-2000006,
     .  "BR(H_",I," -> ~t_2 ~tbar_2)"
        IF(BRHSQ(I,7).GT.0.D0)
     .  WRITE(18,923) BRHSQ(I,7),2,1000006,-2000006,
     .  "BR(H_",I," -> ~t_1 ~tbar_2)"
        IF(BRHSQ(I,7).GT.0.D0)
     .  WRITE(18,923) BRHSQ(I,7),2,2000006,-1000006,
     .  "BR(H_",I," -> ~t_2 ~tbar_1)"
        IF(BRHSQ(I,8).GT.0.D0)
     .  WRITE(18,923) BRHSQ(I,8),2,1000005,-1000005,
     .  "BR(H_",I," -> ~b_1 ~bbar_1)"
        IF(BRHSQ(I,9).GT.0.D0)
     .  WRITE(18,923) BRHSQ(I,9),2,2000005,-2000005,
     .  "BR(H_",I," -> ~b_2 ~bbar_2)"
        IF(BRHSQ(I,10).GT.0.D0)
     .  WRITE(18,923) BRHSQ(I,10),2,1000005,-2000005,
     .  "BR(H_",I," -> ~b_1 ~bbar_2)"
        IF(BRHSQ(I,10).GT.0.D0)
     .  WRITE(18,923) BRHSQ(I,10),2,2000005,-1000005,
     .  "BR(H_",I," -> ~b_2 ~bbar_1)"
        IF(BRHSL(I,1).GT.0.D0)
     .  WRITE(18,923) BRHSL(I,1),2,1000011,-1000011,
     .  "BR(H_",I," -> ~e_L ~ebar_L)"
        IF(BRHSL(I,1).GT.0.D0)
     .  WRITE(18,923) BRHSL(I,1),2,1000013,-1000013,
     .  "BR(H_",I," -> ~mu_L ~mubar_L)"
        IF(BRHSL(I,2).GT.0.D0)
     .  WRITE(18,923) BRHSL(I,2),2,2000011,-2000011,
     .  "BR(H_",I," -> ~e_R ~ebar_R)"
        IF(BRHSL(I,2).GT.0.D0)
     .  WRITE(18,923) BRHSL(I,2),2,2000013,-2000013,
     .  "BR(H_",I," -> ~mu_R ~mubar_R)"
        IF(BRHSL(I,3).GT.0.D0)
     .  WRITE(18,923) BRHSL(I,3),2,1000012,-1000012,
     .  "BR(H_",I," -> ~nu_e_L ~nu_ebar_L)"
        IF(BRHSL(I,3).GT.0.D0)
     .  WRITE(18,923) BRHSL(I,3),2,1000014,-1000014,
     .  "BR(H_",I," -> ~nu_mu_L ~nu_mubar_L)"
        IF(BRHSL(I,4).GT.0.D0)
     .  WRITE(18,923) BRHSL(I,4),2,1000015,-1000015,
     .  "BR(H_",I," -> ~tau_1 ~taubar_1)"
        IF(BRHSL(I,5).GT.0.D0)
     .  WRITE(18,923) BRHSL(I,5),2,2000015,-2000015,
     .  "BR(H_",I," -> ~tau_2 ~taubar_2)"
        IF(BRHSL(I,6).GT.0.D0)
     .  WRITE(18,923) BRHSL(I,6),2,1000015,-2000015,
     .  "BR(H_",I," -> ~tau_1 ~taubar_2)"
        IF(BRHSL(I,6).GT.0.D0)
     .  WRITE(18,923) BRHSL(I,6),2,2000015,-1000015,
     .  "BR(H_",I," -> ~tau_2 ~taubar_1)"
        IF(BRHSL(I,7).GT.0.D0)
     .  WRITE(18,923) BRHSL(I,7),2,1000016,-1000016,
     .  "BR(H_",I," -> ~nu_tau_L ~nu_taubar_L)"
        IF(BRHSL(I,8).GT.0.D0)
     .  WRITE(18,923) BRHSL(I,8),2,2000012,-2000012,
     .  "BR(H_",I," -> ~nu_e_R ~nu_ebar_R)"
        IF(BRHSL(I,8).GT.0.D0)
     .  WRITE(18,923) BRHSL(I,8),2,2000014,-2000014,
     .  "BR(H_",I," -> ~nu_mu_R ~nu_mubar_R)"
        IF(BRHSL(I,9).GT.0.D0)
     .  WRITE(18,923) BRHSL(I,9),2,2000016,-2000016,
     .  "BR(H_",I," -> ~nu_tau_R ~nu_taubar_R)"
      ENDDO


        WRITE(18,904) 36,WIDTH(4),"Pseudoscalar"
        WRITE(18,905) BRJJ(4),2,21,21,"BR(H_A -> gluon gluon)"
        WRITE(18,905) BREE(4),2,11,-11,"BR(H_A -> e e)"
        WRITE(18,905) BRMM(4),2,13,-13,"BR(H_A -> muon muon)"
        WRITE(18,905) BRLL(4),2,15,-15,"BR(H_A -> tau tau)"
        WRITE(18,905) BRSS(4),2,3,-3,"BR(H_A -> s sbar)"
        WRITE(18,905) BRCC(4),2,4,-4,"BR(H_A -> c cbar)"
        WRITE(18,905) BRBB(4),2,5,-5,"BR(H_A -> b bbar)"
        WRITE(18,905) BRTT(4),2,6,-6,"BR(H_A -> t tbar)"
        WRITE(18,905) BRGG(4),2,22,22,"BR(H_A -> gamma gamma)"
        WRITE(18,905) BRZG(4),2,23,22,"BR(H_A -> Z1 gamma)"
        IF(BRAHZ(1).GT.0.D0)
     .  WRITE(18,905) BRAHZ(1),2,23,25,"BR(H_A -> Z1 H_1)"
        IF(BRAHZ(2).GT.0.D0)
     .  WRITE(18,905) BRAHZ(2),2,23,35,"BR(H_A -> Z1 H_2)"
      DO J=1,6
      DO K=J,6
        IF(BRNEU(4,J,K).GT.0.D0)
     .  WRITE(18,925) BRNEU(4,J,K),2,PDGNEU(J),PDGNEU(K),
     .  "BR(H_A -> neu_",J," neu_",K,")"
      ENDDO
      ENDDO
        IF(BRCHA(4,1).GT.0.D0)
     .  WRITE(18,905) BRCHA(4,1),2,1000024,-1000024,
     .  "BR(H_A -> cha_1 cha_1bar)"
        IF(BRCHA(4,2).GT.0.D0)
     .  WRITE(18,905) BRCHA(4,2),2,1000024,-1000037,
     .  "BR(H_A -> cha_1 cha_2bar)"
        IF(BRCHA(4,2).GT.0.D0)
     .  WRITE(18,905) BRCHA(4,2),2,1000037,-1000024,
     .  "BR(H_A -> cha_2 cha_1bar)"
        IF(BRCHA(4,3).GT.0.D0)
     .  WRITE(18,905) BRCHA(4,3),2,1000037,-1000037,
     .  "BR(H_A -> cha_2 cha_2bar)"
        IF(BRASQ(1).GT.0.D0)
     .  WRITE(18,905) BRASQ(1),2,1000006,-2000006,
     .  "BR(H_A -> ~t_1 ~tbar_2)"
        IF(BRASQ(1).GT.0.D0)
     .  WRITE(18,905) BRASQ(1),2,2000006,-1000006,
     .  "BR(H_A -> ~t_2 ~tbar_1)"
        IF(BRASQ(2).GT.0.D0)
     .  WRITE(18,905) BRASQ(2),2,1000005,-2000005,
     .  "BR(H_A -> ~b_1 ~bbar_2)"
        IF(BRASQ(2).GT.0.D0)
     .  WRITE(18,905) BRASQ(2),2,2000005,-1000005,
     .  "BR(H_A -> ~b_2 ~bbar_1)"
        IF(BRASL.GT.0.D0)
     .  WRITE(18,905) BRASL,2,1000015,-2000015,
     .  "BR(H_A -> ~tau_1 ~taubar_2)"
        IF(BRASL.GT.0.D0)
     .  WRITE(18,905) BRASL,2,2000015,-1000015,
     .  "BR(H_A -> ~tau_2 ~taubar_1)"


        WRITE(18,904) 37,HCWIDTH,"Charged Higgs"
        WRITE(18,905) HCBRE,2,-11,12,"BR(H+ -> e nu_e)"
        WRITE(18,905) HCBRM,2,-13,14,"BR(H+ -> muon nu_muon)"
        WRITE(18,905) HCBRL,2,-15,16,"BR(H+ -> tau nu_tau)"
        WRITE(18,905) HCBRSU,2,2,-3,"BR(H+ -> u sbar)"
        WRITE(18,905) HCBRSC,2,4,-3,"BR(H+ -> c sbar)"
        WRITE(18,905) HCBRBU,2,2,-5,"BR(H+ -> u bbar)"
        WRITE(18,905) HCBRBC,2,4,-5,"BR(H+ -> c bbar)"
        WRITE(18,905) HCBRBT,2,6,-5,"BR(H+ -> t bbar)"
        IF(HCBRWH(1).GT.0.D0)
     .  WRITE(18,905) HCBRWH(1),2,24,25,"BR(H+ -> W+ H_1)"
        IF(HCBRWH(2).GT.0.D0)
     .  WRITE(18,905) HCBRWH(2),2,24,35,"BR(H+ -> W+ H_2)"
        IF(HCBRWH(4).GT.0.D0)
     .  WRITE(18,905) HCBRWH(4),2,24,36,"BR(H+ -> W+ H_A)"
      DO J=1,2
      DO K=1,6
        IF(HCBRNC(K,J).GT.0.D0)
     .  WRITE(18,925) HCBRNC(K,J),2,PDGCHA(J),PDGNEU(K),
     .  "BR(H+ -> cha_",J," neu_",K,")"
      ENDDO
      ENDDO
        IF(HCBRSQ(1).GT.0.D0)
     .  WRITE(18,905) HCBRSQ(1),2,1000002,-1000001,
     .  "BR(H+ -> ~u_L ~dbar_L)"
        IF(HCBRSQ(1).GT.0.D0)
     .  WRITE(18,905) HCBRSQ(1),2,1000004,-1000003,
     .  "BR(H+ -> ~c_L ~sbar_L)"
        IF(HCBRSQ(2).GT.0.D0)
     .  WRITE(18,905) HCBRSQ(2),2,1000006,-1000005,
     .  "BR(H+ -> ~t_1 ~bbar_1)"
        IF(HCBRSQ(3).GT.0.D0)
     .  WRITE(18,905) HCBRSQ(3),2,1000006,-2000005,
     .  "BR(H+ -> ~t_1 ~bbar_2)"
        IF(HCBRSQ(4).GT.0.D0)
     .  WRITE(18,905) HCBRSQ(4),2,2000006,-1000005,
     .  "BR(H+ -> ~t_2 ~bbar_1)"
        IF(HCBRSQ(5).GT.0.D0)
     .  WRITE(18,905) HCBRSQ(5),2,2000006,-2000005,
     .  "BR(H+ -> ~t_2 ~bbar_2)"
        IF(HCBRSL(1).GT.0.D0)
     .  WRITE(18,905) HCBRSL(1),2,1000012,-1000011,
     .  "BR(H+ -> ~nu_e_L ~ebar_L)"
        IF(HCBRSL(1).GT.0.D0)
     .  WRITE(18,905) HCBRSL(1),2,1000014,-1000013,
     .  "BR(H+ -> ~nu_mu_L ~mubar_L)"
        IF(HCBRSL(2).GT.0.D0)
     .  WRITE(18,905) HCBRSL(2),2,1000016,-1000015,
     .  "BR(H+ -> ~nu_tau_L ~taubar_1)"
        IF(HCBRSL(3).GT.0.D0)
     .  WRITE(18,905) HCBRSL(3),2,1000016,-2000015,
     .  "BR(H+ -> ~nu_tau_L ~taubar_2)"




        WRITE(18,904) 6,toptot,"Top Quark"
        WRITE(18,905) brtopbw,2,5,24,"BR(t -> b W+)"
        IF(brtopbh.GT.0.D0)
     .  WRITE(18,905) brtopbh,2,5,37,"BR(t -> b H+)"
      DO J=1,6
      DO K=1,2
        IF(brtopneutrstop(J,K).GT.0.D0)
     .  WRITE(18,925) brtopneutrstop(J,K),2,PDGST(K),PDGNEU(J),
     .  "BR(t -> ~t_",K," neu_",J,")"
      ENDDO
      ENDDO



        WRITE(16,899) "# SM HIGGS BRANCHING RATIOS IN SLHA FORMAT"
        WRITE(16,899) "# Info about decay package"
        WRITE(16,899) "BLOCK DCINFO        # Program information"
        WRITE(16,900) 1,"UMSSMTools # Decay package"
        WRITE(16,900) 2,"1.1.2      # Version number"

        WRITE(16,899) "#           PDG          Width"
      DO I=1,3
        IF(I.EQ.1)
     .  WRITE(16,904) 25,WIDTHSM(I),"Lightest neutral Higgs scalar"
        IF(I.EQ.2)
     .  WRITE(16,904) 35,WIDTHSM(I),"2nd neutral Higgs scalar"
        IF(I.EQ.3)
     .  WRITE(16,904) 45,WIDTHSM(I),"3rd neutral Higgs scalar"
        WRITE(16,923) BRJJSM(I),2,21,21,"BR(H_",I," -> gluon gluon)"
        WRITE(16,923) BREESM(I),2,11,-11,"BR(H_",I," -> e e)"
        WRITE(16,923) BRMMSM(I),2,13,-13,"BR(H_",I," -> muon muon)"
        WRITE(16,923) BRLLSM(I),2,15,-15,"BR(H_",I," -> tau tau)"
        WRITE(16,923) BRSSSM(I),2,3,-3,"BR(H_",I," -> s sbar)"
        WRITE(16,923) BRCCSM(I),2,4,-4,"BR(H_",I," -> c cbar)"
        WRITE(16,923) BRBBSM(I),2,5,-5,"BR(H_",I," -> b bbar)"
        WRITE(16,923) BRTTSM(I),2,6,-6,"BR(H_",I," -> t tbar)"
        WRITE(16,923) BRWWSM(I),2,24,-24,"BR(H_",I," -> W+ W-)"
        WRITE(16,923) BRZZSM(I),2,23,23,"BR(H_",I," -> Z Z)"
        WRITE(16,923) BRGGSM(I),2,22,22,"BR(H_",I," -> gamma gamma)"
        WRITE(16,923) BRZGSM(I),2,23,22,"BR(H_",I," -> Z gamma)"
      ENDDO

        WRITE(17,899) "# "
        WRITE(17,899) "# Higgs reduced couplings"
        WRITE(17,899) "# (as compared to a SM Higgs with same mass)"
        WRITE(17,899) "BLOCK REDCOUP"
      DO I=1,3
        WRITE(17,926) "# H",I
        WRITE(17,903) I,1,CU(I),"U-type fermions"
        WRITE(17,903) I,2,CD(I),"D-type fermions"
        WRITE(17,903) I,3,CB(I),"b-quarks"
        WRITE(17,903) I,4,CV(I),"W boson"
        WRITE(17,903) I,5,CVZ(I),"Z1 boson"
        WRITE(17,903) I,6,CJ(I),"Gluons"
        WRITE(17,903) I,7,CG(I),"Photons"
      ENDDO
        WRITE(17,899) "# HA"
        WRITE(17,903) 4,1,CU(4),"U-type fermions"
        WRITE(17,903) 4,2,CD(4),"D-type fermions"
        WRITE(17,903) 4,3,CB(4),"b-quarks"
        WRITE(17,903) 4,4,0.,"W boson"
        WRITE(17,903) 4,5,0.,"Z1 boson"
        WRITE(17,903) 4,6,CJ(4),"Gluons"
        WRITE(17,903) 4,7,CG(4),"Photons"


 899  FORMAT(A)
 900  FORMAT(1X,I5,3X,A)
 901  FORMAT(1X,I5,3X,1P,E16.8,0P,3X,'#',1X,A)
 902  FORMAT(1X,I9,3X,1P,E16.8,0P,3X,'#',1X,A)
 903  FORMAT(1X,I2,1X,I2,3X,1P,E16.8,0P,3X,'#',1X,A)
 904  FORMAT('DECAY',1X,I9,3X,1P,E16.8,0P,3X,'#',1X,A)
 905  FORMAT(3X,1P,E16.8,0P,3X,I2,3X,I9,1X,I9,1X,2X,'#',1X,A)
 906  FORMAT('#',1X,A,3X,E16.8)
 907  FORMAT(A,1P,E16.8,A)
* 908	FORMAT(2E16.8)
* 909	FORMAT(8X,A,12X,A)
* 910   FORMAT(E16.8,3X,A)
* 911   FORMAT(A,F6.3,A,F6.3,A)
* 912   FORMAT(1X,1P,E12.4,0P,1X,I2,1X,I5,1X,I5,1X,I5,1X,'#',1X,A)
* 913   FORMAT(1X,1P,E12.4,3X,1P,E12.4,0P,1X,I3,1X,I5,1X,I5,1X,I4,
*     . 1X,'#',1X,A)
 918  FORMAT(1X,I5,3X,A,F6.1,'-',F5.1,A)
 919  FORMAT(1X,I5,3X,A,F6.2)
* 920  FORMAT('#',0P,I5,3X,1P,E16.8,0P,3X,'#',1X,A)
 921  FORMAT(1X,I2,1X,I2,3X,1P,E16.8,0P,3X,'#',1X,A,I1,A,I1,A)
 922  FORMAT(1X,I9,3X,1P,E16.8,0P,3X,'#',1X,A,I1,A)
 923  FORMAT(3X,1P,E16.8,0P,3X,I2,3X,I9,1X,I9,1X,2X,'#',1X,A,I1,A)
 924  FORMAT(3X,1P,E16.8,0P,3X,I2,3X,I9,1X,I9,1X,2X,'#',1X,A,I1,A,
     . I1,A,I1,A)
 925  FORMAT(3X,1P,E16.8,0P,3X,I2,3X,I9,1X,I9,1X,2X,'#',1X,A,I1,A,
     . I1,A)
 926  FORMAT(A,I1)

        END
