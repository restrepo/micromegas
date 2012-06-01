	SUBROUTINE RUNPAR(PAR)

*******************************************************************
* Subroutine to evolve the parameters from the scale Q2 to QSTSB
*
* The mass scale QSTSB of the 3rd generation squarks is
* computed here as QSTSB = MQ3*MU3 (QSTSB has the dimension mass^2)
*
* Several parameters at the scale QSTSB (relevant for the
* calculation of the Higgs masses) are computed here and stored in 
* COMMON/QGAUGE, /QHIGGS, /QQUARK, /QNMPAR:
* - The electroweak and strong gauge couplings (the electroweak 
*   couplings appear in the tree level Higgs mass matrix in MHIGGS,
*   alpha_s only in the two loop Higgs mass corrections)
* - the Higgs wave function normalization constants ZHU, ZHD, ZS 
*   and the Higgs vevs H1Q, H2Q and TANBQ
* - the top/bottom Yukawa couplings HTQ/HBQ and masses MTOPQ/MBOTQ
* - the NMSSM parameters LQ, KQ, ALQ, AKQ
*   (initially L, K, Mueff, AL and AK are defined at Q2)
* - if MAFLAG = 0: AL at Q2 is given in PAR(5),
*      then: MUQ=LQ*S, NUQ=KQ*S, MA (at QSTSB, stored in PAR(23))
*      and ALQ are computed here
* - if MAFLAG = 1: MA at QSTSB is given in PAR(23)
*      then: MUQ, NUQ, AL(at Q2, stored in PAR(5)) and ALQ 
*      are computed here
* - if MAFLAG < 0: this is the case NMSPEC and NMGMSB
*      then: ALQ, AKQ, LQ (and, provisionally, KQ) are still
*      computed here, but MUQ, NUQ, KQ and MA are subsequently
*      computed in the subroutine LOWMUK
*
* The SUSY scale Q2, where the soft terms are
* defined on input, is possibly much larger than QSTSB.
* Unless Q2 is defined by the user, Q2 is defined here in terms
* of the first generation squark masses as 
* MSUSY**2 == Q2 = MAX((2*MQ**2+MU**2+MD**2)/4,Q2MIN).
*
*******************************************************************

	IMPLICIT NONE

	INTEGER Q2FIX,OMGFLAG,MAFLAG,FLAG

	DOUBLE PRECISION PAR(*)
	DOUBLE PRECISION PI,COEF
	DOUBLE PRECISION tanb,SB2,CB2,h1,h2,M1,M2,HTAU
	DOUBLE PRECISION L,K,AL,AK,MU,NU,RUNMB,LQT,ALSMT,HT,HB	
	DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
	DOUBLE PRECISION Q2MIN,Q2,QSTSB,MA2
	DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
	DOUBLE PRECISION G1Q,G2Q,GQ,ALSQ
	DOUBLE PRECISION Lmu,Lnu,LM1mu,LM2MU,Lmunu,LQ2,LMAMT
	DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
	DOUBLE PRECISION HTQ,HBQ,MTOPQ,MBOTQ
	DOUBLE PRECISION LQ,KQ,ALQ,AKQ,MUQ,NUQ
		
	COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
	COMMON/RENSCALE/Q2
	COMMON/STSBSCALE/QSTSB
        COMMON/Q2FIX/Q2MIN,Q2FIX
	COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
	COMMON/QGAUGE/G1Q,G2Q,GQ,ALSQ
	COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
	COMMON/QQUARK/HTQ,HBQ,MTOPQ,MBOTQ
	COMMON/QNMPAR/LQ,KQ,ALQ,AKQ,MUQ,NUQ
	COMMON/FLAGS/OMGFLAG,MAFLAG

	PI=4.D0*DATAN(1.D0)
	COEF=1.D0/(16.D0*PI**2)
	FLAG=0

	!PRINT*,"CALL RUNPAR"
	!PRINT*,""

*   Definition of the SUSY scale Q2, unless defined initially
	IF(Q2FIX.EQ.0) THEN
	  Q2=MAX((2.D0*PAR(15)+PAR(16)+PAR(17))/4.D0,Q2MIN)
	ENDIF

	!PRINT*,"QSUSY =",DSQRT(Q2)

*   Definition of the scale QSTSB
	QSTSB=DSQRT(MAX(PAR(7)*PAR(8),Q2MIN**2))
c	QSTSB=MAX((2.D0*PAR(7)+PAR(8)+PAR(9))/4,Q2MIN)

	!PRINT*,"QSTSB =",DSQRT(QSTSB)

	LQ2=DLOG(QSTSB/Q2)
	LQT=DLOG(MAX(QSTSB,MT**2)/MT**2)

*   NMSSM parameters
	L=PAR(1)
	K=PAR(2)
	tanb=PAR(3)	
	MU=PAR(4)
	AK=PAR(6)
	NU=K/L*MU
	M1=PAR(20)
	M2=PAR(21)

	IF(MAFLAG.LT.0)THEN

	 AL=PAR(5)
	 MA2=PAR(23)**2

	ELSEIF(MAFLAG.EQ.0)THEN

	 AL=PAR(5)
	 MA2=MAX((AL+NU)*MU*(tanb+1.D0/tanb),0.D0)

	ELSEIF(MAFLAG.EQ.1)THEN

	 MA2=PAR(23)**2

	ENDIF

	!PRINT*,"MA =",DSQRT(MA2)

*   Trig. function of beta

	SB2=tanb**2/(1.D0+tanb**2)
	CB2=1.D0-SB2

*   Higgs vevs h1 = <hu>, h2 = <hd>:	
	h1=DSQRT(SB2/(2.D0*DSQRT(2.D0)*GF))
	h2=h1/tanb
	
	HTAU=MTAU/H2

*   Electroweak gauge couplings at QSTSB:

*   g_2**2 including the Higgs and sparticle thresholds:

	g2q=g2/(1.D0+g2*COEF*(DLOG(QSTSB/MZ**2)*19.D0/6.D0
     C	    -DLOG(QSTSB/MIN(QSTSB,MAX(MA2,MZ**2)))/6.D0
     C	    -DLOG(QSTSB/MIN(QSTSB,MAX(MU**2,MZ**2)))*2.D0/3.D0
     C	    -DLOG(QSTSB/MIN(QSTSB,MAX(PAR(7),MZ**2)))/2.D0
     C	    -DLOG(QSTSB/MIN(QSTSB,MAX(PAR(10),MZ**2)))/6.D0
     C	    -DLOG(QSTSB/MIN(QSTSB,MAX(PAR(15),MZ**2))) 
     C	    -DLOG(QSTSB/MIN(QSTSB,MAX(PAR(18),MZ**2)))/3.D0
     C	    -DLOG(QSTSB/MIN(QSTSB,MAX(M2**2,MZ**2)))*4.D0/3.D0))

*   g_1**2 including the top, Higgs and sparticle thresholds:

	g1q=g1/(1.D0-g1*COEF*(DLOG(QSTSB/MZ**2)*53.D0/9.D0
     C	    +DLOG(QSTSB/MT**2)*17.D0/18.D0
     C	    +DLOG(QSTSB/MIN(QSTSB,MAX(MA2,MZ**2)))/6.D0
     C	    +DLOG(QSTSB/MIN(QSTSB,MAX(MU**2,MZ**2)))*2.D0/3.D0
     C	    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(7),MZ**2)))/18.D0
     C	    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(8),MZ**2)))*4.D0/9.D0
     C	    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(9),MZ**2)))/9.D0
     C	    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(10),MZ**2)))/6.D0
     C	    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(11),MZ**2)))/3.D0
     C	    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(15),MZ**2)))/9.D0
     C	    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(16),MZ**2)))*8.D0/9.D0
     C	    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(17),MZ**2)))*2.D0/9.D0
     C	    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(18),MZ**2)))/3.D0
     C	    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(19),MZ**2)))*2.D0/3.D0))

	gq=(g1q+g2q)/2.D0

*   Alphas at MT and QSTSB 

	ALSMT=ALSMZ/(1.D0+23.D0/(12.D0*PI)*ALSMZ*DLOG(MT**2/MZ**2))
	ALSQ=ALSMT/(1.D0+ALSMT/(4.D0*PI)*(7.D0*LQT-2.D0*
     C    DLOG(MAX(QSTSB,PAR(22)**2)/MAX(PAR(22)**2,MT**2))))

*   Yukawas at MT (ht, hb: running MS_bar couplings)

	HT=MT/(1.D0+4.D0*ALSMT/(3.D0*PI)+11.D0*(ALSMT/PI)**2)/H1
	HB=RUNMB(MT)/H2

*   Logs for the Wave Function Renormalization Constants

	Lmu = DLOG(MIN(MAX(mu**2,MZ**2),QSTSB)/QSTSB)
	Lnu = DLOG(MIN(MAX(4.D0*nu**2,MZ**2),QSTSB)/QSTSB)
	LM1mu = DLOG(MIN(MAX(M1**2,mu**2,MZ**2),QSTSB)/QSTSB)
	LM2mu = DLOG(MIN(MAX(M2**2,mu**2,MZ**2),QSTSB)/QSTSB)
	Lmunu = DLOG(MIN(MAX(mu**2,4.D0*nu**2,MZ**2),QSTSB)/QSTSB)
	LMAMT = DLOG(MIN(MAX(MA2,MT**2),QSTSB)/MT**2)
	
	ZHU=1.D0+COEF*(3.D0*ht**2*LQT
     C      +CB2*(3*HB**2+(MTAU/H2)**2-3*HT**2)*LMAMT
     C      -G1Q/2.D0*LM1mu-3.D0*G2Q/2.D0*LM2mu
     C	    -3.D0*(G1Q+3.D0*G2Q)/4.D0*DLOG(QSTSB/MZ**2)
     C	    -L**2*LMUNU)
	
	ZHD=1.D0+COEF*(3.D0*hb**2*LQT+(MTAU/H2)**2*DLOG(QSTSB/MZ**2)
     C      +SB2*(3*HT**2-3*HB**2-(MTAU/H2)**2)*LMAMT
     C	    -G1Q/2.D0*LM1mu-3.D0*G2Q/2.D0*LM2mu
     C	    -3.D0*(G1Q+3.D0*G2Q)/4.D0*DLOG(QSTSB/MZ**2)
     C	    -L**2*LMUNU)
	    
	ZS=1.D0-2.D0*COEF*(L**2*Lmu+K**2*Lnu)

*   Higgs Vevs at QSTSB

	H1Q=H1/DSQRT(ZHU)
	H2Q=H2/DSQRT(ZHD)
	
	TANBQ=H1Q/H2Q

*   Top/Bottom Yukawas at QSTSB
*   (Note: RUNMB(Q) includes QCD corrections only)
*   including electroweak contributions
	HTQ=HT*(1.D0+7.D0/(4.D0*PI)*ALSMT*LQT)**(-4.D0/7.D0)
     C	 *(1.D0+COEF/4.D0*((-17.D0/6.D0*g1q-9.D0/2.D0*g2q
     C	     +9.D0*ht**2+hb**2)*LQT
     C   +(-9.D0*CB2*HT**2+(3.D0*CB2-1.D0)*HB**2+2.D0*HTAU**2*CB2)*LMAMT
     C    -2.D0*L**2*Lmunu-G1Q*LM1mu-3.D0*G2Q*LM2mu))

	HBQ=RUNMB(DSQRT(QSTSB))/H2Q
     C	 *(1.D0+COEF/4.D0*((-5.D0/6.D0*g1q-9.D0/2.D0*g2q
     C       +9.D0*hb**2+ht**2+2.D0*HTAU**2)*LQT
     C   +(-9.D0*SB2*HB**2+(3.D0*SB2-1.D0)*HT**2-2.D0*HTAU**2*SB2)*LMAMT
     C    -2.D0*L**2*Lmunu-G1Q*LM1mu-3.D0*G2Q*LM2mu))

*   Conversion to DR_bar:
*	HTQ=HTQ*(1.D0-ALSQ/(3.D0*PI)+g2q*COEF*3.D0/8.D0)
*	HBQ=HBQ*(1.D0-ALSQ/(3.D0*PI)+g2q*COEF*3.D0/8.D0)
	
*   Running Top and Bottom Quark Masses

     	MTOPQ=HTQ*H1Q
	MBOTQ=HBQ*H2Q

*   NMSSM Parameters at QSTSB, stored in COMMON/QNMPAR
	
	KQ=K*(1.D0+3.D0*COEF*(L**2+K**2)*LQ2)
	LQ=L*(1.D0+COEF/2.D0*(-G1Q-3.D0*G2Q+4.D0*L**2+2.D0*KQ**2
     C	  +3.D0*(HTQ**2+HBQ**2)+(MTAU/H2)**2)*LQ2)

	!PRINT*,"LQ =",LQ
	!PRINT*,"KQ =",KQ

	IF(MAFLAG.LT.0)THEN

	 ALQ=AL+COEF*(G1*M1+3.D0*G2*M2+4.D0*LQ**2*AL+2.D0*KQ**2*AK
     C	  +3.D0*(HTQ**2*PAR(12)+HBQ**2*PAR(13))
     C	  +(MTAU/H2)**2*PAR(14))*LQ2
	 AKQ=AK+6.D0*COEF*(LQ**2*ALQ+KQ**2*AK)*LQ2

	 !PRINT*,"ALQ =",ALQ
	 !PRINT*,"AKQ =",AKQ

	ELSEIF(MAFLAG.EQ.0)THEN

	 MUQ=MU*(1.D0+COEF/2.D0*(-G1Q-3.D0*G2Q+2.D0*LQ**2
     C	  +3.D0*(HTQ**2+HBQ**2)+(MTAU/H2)**2)*LQ2)
	 NUQ=MUQ*KQ/LQ
	 ALQ=AL+COEF*(G1*M1+3.D0*G2*M2+4.D0*LQ**2*AL+2.D0*KQ**2*AK
     C	  +3.D0*(HTQ**2*PAR(12)+HBQ**2*PAR(13))
     C	  +(MTAU/H2)**2*PAR(14))*LQ2
	 AKQ=AK+6.D0*COEF*(LQ**2*ALQ+KQ**2*AK)*LQ2
	 PAR(23)=DSQRT(MAX((ALQ+KQ/LQ*MUQ)*MUQ*(tanbQ+1.D0/tanbQ),0.D0))

	 !PRINT*,"MUQ =",MUQ
	 !PRINT*,"NUQ =",NUQ
	 !PRINT*,"ALQ =",ALQ
	 !PRINT*,"AKQ =",AKQ
	 !PRINT*,"MA =",PAR(23)

	ELSEIF(MAFLAG.EQ.1)THEN

	 MUQ=MU*(1.D0+COEF/2.D0*(-G1Q-3.D0*G2Q+2.D0*LQ**2
     C	    +3.D0*(HTQ**2+HBQ**2)+(MTAU/H2)**2)*LQ2)
	 NUQ=MUQ*KQ/LQ
	 ALQ=MA2*TANBQ/(MUQ*(1.D0+TANBQ**2))-MUQ*KQ/LQ
	 AKQ=AK+6.D0*COEF*(LQ**2*ALQ+KQ**2*AK)*LQ2
	 PAR(5)=ALQ-COEF*(G1*M1+3.D0*G2*M2+4.D0*LQ**2*ALQ+2.D0*KQ**2*AK
     C	  +3.D0*(HTQ**2*PAR(12)+HBQ**2*PAR(13))+(MTAU/H2)**2*PAR(14))*LQ2

	 !PRINT*,"MUQ =",MUQ
	 !PRINT*,"NUQ =",NUQ
	 !PRINT*,"ALQ =",ALQ
	 !PRINT*,"AKQ =",AKQ
	 !PRINT*,"ALS =",PAR(5)

	ENDIF

	!PRINT*,""
	!PRINT*,""

	END



