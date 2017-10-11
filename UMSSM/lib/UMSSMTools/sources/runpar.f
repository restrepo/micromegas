	SUBROUTINE RUNPAR()

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
*
* The SUSY scale Q2, where the soft terms are
* defined on input, is possibly much larger than QSTSB.
* Unless Q2 is defined by the user, Q2 is defined here in terms
* of the first generation squark masses as 
* MSUSY**2 == Q2 = MAX((2*MQ**2+MU**2+MD**2)/4,Q2MIN).
*
*******************************************************************

	IMPLICIT NONE

	INTEGER Q2FIX

	DOUBLE PRECISION PI,COEF
	DOUBLE PRECISION tanb,SB2,CB2,h1,h2,M1,M2,HTAU
	DOUBLE PRECISION L,K,AL,AK,MU,NU,RUNMB,LQT,ALSMT,HT,HB	
	DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
	DOUBLE PRECISION Q2MIN,Q2,QSTSB,MA2
	DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
	DOUBLE PRECISION G1Q,G2Q,GQ,ALSQ
	DOUBLE PRECISION Lmu,Lnu,LM1mu,LM2MU,Lmunu,LQ2,LMAMT
	DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ,HSQ
	DOUBLE PRECISION HTQ,HBQ,MTOPQ,MBOTQ
	DOUBLE PRECISION INTEG,DELMB,AU,M3
	DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
	DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
	DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU,SST,SSB
	DOUBLE PRECISION LQ,KQ,ALQ,AKQ,MUQ,NUQ
	DOUBLE PRECISION UPARF,SAZZ,CAZZ,VEV,NCP
	DOUBLE PRECISION QD,QU,QS,VEVS,G1P,QQ,QUP,AD,ATAU
	DOUBLE PRECISION QDOW,QL,QE,QN
	DOUBLE PRECISION ALSMA,DLA,DLQA,F1,F2,HTMA
	DOUBLE PRECISION M2Q3,M2U3,M2D3,M2L3,M2R3
		
	COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
	COMMON/RENSCALE/Q2
	COMMON/STSBSCALE/QSTSB
	COMMON/Q2FIX/Q2MIN,Q2FIX
	COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
	COMMON/QGAUGE/G1Q,G2Q,GQ,ALSQ
	COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     C		MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     C		CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
	COMMON/DELMB/DELMB
	COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ,HSQ
	COMMON/QQUARK/HTQ,HBQ,MTOPQ,MBOTQ
	COMMON/QNMPAR/LQ,KQ,ALQ,AKQ,MUQ,NUQ
	COMMON/UMSSM/SAZZ,CAZZ,VEV,NCP,QD,QU,QS,VEVS,G1P,QQ,
     C		QUP,QDOW,QL,QE,QN  
	COMMON/WARN/F2
	COMMON/NOBUG/TANB,AU,AD,M2,SST,SSB,L,AL,ATAU

	PI=4.D0*DATAN(1.D0)
	COEF=1.D0/(16.D0*PI**2)
	SST= UPARF(84)
	SSB= UPARF(88)
	AU= UPARF(17)
	AD= UPARF(18)
	ATAU= UPARF(19)

	!PRINT*,"CALL RUNPAR"
	!PRINT*,""

*   Definition of the SUSY scale Q2, unless defined initially
	IF(Q2FIX.EQ.0) THEN
	  Q2=MAX((2.D0*MUR**2+MDL**2+MDR**2)/4.D0,Q2MIN)
	ENDIF

	!PRINT*,"QSUSY =",DSQRT(Q2)

*   Definition of the scale QSTSB
	QSTSB=DSQRT(MAX((MST1*MST2)**2,Q2MIN**2))
c	QSTSB=MAX((2.D0*PAR(7)+PAR(8)+PAR(9))/4,Q2MIN)

	!PRINT*,"QSTSB =",DSQRT(QSTSB)

	LQ2=DLOG(QSTSB/Q2)
	LQT=DLOG(MAX(QSTSB,MT**2)/MT**2)

*   UMSSM parameters
	L=UPARF(45)
	K=0.D0
	tanb=UPARF(44)	
	MU=UPARF(38)
	AK=0.D0
	NU=K/L*MU
	M1=UPARF(0)
	M2=UPARF(1)
	AL=UPARF(39)

*   Trig. function of beta

	SB2=tanb**2/(1.D0+tanb**2)
	CB2=1.D0-SB2

	MA2=L*AL*VEVS/DSQRT(2d0)*
     C	    (TANB+1/TANB+DSQRT(SB2*CB2)*(VEV/VEVS)**2)

*   Higgs vevs h1 = <hd>, h2 = <hu>:	
	h2=DSQRT(SB2/(2.D0*DSQRT(2.D0)*GF))
	h1=h2/tanb
	
	HTAU=MTAU/H1

*   Electroweak gauge couplings at QSTSB:

	M2Q3=UPARF(358)
	M2U3=UPARF(359)
	M2D3=UPARF(360)
	M2L3=UPARF(361)
	M2R3=UPARF(362)


*   g_2**2 including the Higgs and sparticle thresholds:

	g2q=g2/(1.D0+g2*COEF*(DLOG(QSTSB/MZ**2)*19.D0/6.D0
     C	    -DLOG(QSTSB/MIN(QSTSB,MAX(MA2,MZ**2)))/6.D0
     C	    -DLOG(QSTSB/MIN(QSTSB,MAX(MU**2,MZ**2)))*2.D0/3.D0
     C	    -DLOG(QSTSB/MIN(QSTSB,MAX(M2Q3,MZ**2)))/2.D0
     C	    -DLOG(QSTSB/MIN(QSTSB,MAX(M2L3,MZ**2)))/6.D0
     C	    -DLOG(QSTSB/MIN(QSTSB,MAX(MDL**2,MZ**2))) 
     C	    -DLOG(QSTSB/MIN(QSTSB,MAX(UPARF(93)**2,MZ**2)))/3.D0
     C	    -DLOG(QSTSB/MIN(QSTSB,MAX(M2**2,MZ**2)))*4.D0/3.D0))

*   g_1**2 including the top, Higgs and sparticle thresholds:

	g1q=g1/(1.D0-g1*COEF*(DLOG(QSTSB/MZ**2)*53.D0/9.D0
     C	    +DLOG(QSTSB/MT**2)*17.D0/18.D0
     C	    +DLOG(QSTSB/MIN(QSTSB,MAX(MA2,MZ**2)))/6.D0
     C	    +DLOG(QSTSB/MIN(QSTSB,MAX(MU**2,MZ**2)))*2.D0/3.D0
     C	    +DLOG(QSTSB/MIN(QSTSB,MAX(M2Q3,MZ**2)))/18.D0
     C	    +DLOG(QSTSB/MIN(QSTSB,MAX(M2U3,MZ**2)))*4.D0/9.D0
     C	    +DLOG(QSTSB/MIN(QSTSB,MAX(M2D3,MZ**2)))/9.D0
     C	    +DLOG(QSTSB/MIN(QSTSB,MAX(M2L3,MZ**2)))/6.D0
     C	    +DLOG(QSTSB/MIN(QSTSB,MAX(M2R3,MZ**2)))/3.D0
     C	    +DLOG(QSTSB/MIN(QSTSB,MAX(MDL**2,MZ**2)))/9.D0
     C	    +DLOG(QSTSB/MIN(QSTSB,MAX(MUR**2,MZ**2)))*8.D0/9.D0
     C	    +DLOG(QSTSB/MIN(QSTSB,MAX(MDR**2,MZ**2)))*2.D0/9.D0
     C	    +DLOG(QSTSB/MIN(QSTSB,MAX(UPARF(93)**2,MZ**2)))/3.D0
     C	    +DLOG(QSTSB/MIN(QSTSB,MAX(UPARF(94)**2,MZ**2)))
     C	    *2.D0/3.D0))

	gq=(g1q+g2q)/2.D0

*   Alphas at MT and QSTSB 

	ALSMT=ALSMZ/(1.D0+23.D0/(12.D0*PI)*ALSMZ*DLOG(MT**2/MZ**2))
	ALSQ=ALSMT/(1.D0+ALSMT/(4.D0*PI)*(7.D0*LQT-2.D0*
     C    DLOG(MAX(QSTSB,UPARF(2)**2)/MAX(UPARF(2)**2,MT**2))))

*   Yukawas at MT (ht, hb: running MS_bar couplings)

	HT=MT/(1.D0+4.D0*ALSMT/(3.D0*PI)+11.D0*(ALSMT/PI)**2)/H2
	HB=RUNMB(MT)/H1

*   Logs for the Wave Function Renormalization Constants

	Lmu = DLOG(MIN(MAX(mu**2,MZ**2),QSTSB)/QSTSB)
	Lnu = DLOG(MIN(MAX(4.D0*nu**2,MZ**2),QSTSB)/QSTSB)
	LM1mu = DLOG(MIN(MAX(M1**2,mu**2,MZ**2),QSTSB)/QSTSB)
	LM2mu = DLOG(MIN(MAX(M2**2,mu**2,MZ**2),QSTSB)/QSTSB)
	Lmunu = DLOG(MIN(MAX(mu**2,4.D0*nu**2,MZ**2),QSTSB)/QSTSB)
	LMAMT = DLOG(MIN(MAX(MA2,MT**2),QSTSB)/MT**2)
    
* Aux. quantities for the resummation of logs ~ht^2*LQT
* and ~ht^2*LMAMT:

	ALSMA=ALSMT/(1.D0+ALSMT/(4.D0*PI)*(7.D0*LMAMT-2.D0*
     C      DLOG(MAX(MA2,UPARF(2)**2)/MAX(UPARF(2)**2,MT**2))))
	DLA=(ALSMA/ALSMT)**(1.D0/7.D0)
	DLQA=(ALSQ/ALSMA)**(1.D0/7.D0)
	F1=1d0-9.D0*SB2*HT**2*(1.D0-DLA)/(8.D0*PI*ALSMT)
	HTMA=HT*DLA**4/DSQRT(DABS(F1))
	F2=1d0-9.D0*SB2*HTMA**2*(1.D0-DLQA)/(8.D0*PI*ALSMA)	
	
	ZHU=(F1*F2)**(-2.D0/3.D0)*(1.D0+COEF*(
     C      +CB2*(3*HB**2+(MTAU/H1)**2)*LMAMT
     C      -G1Q/2.D0*LM1mu-3.D0*G2Q/2.D0*LM2mu
     C	    -3.D0*(G1Q+3.D0*G2Q)/4.D0*DLOG(QSTSB/MZ**2)
     C	    -L**2*LMUNU))
	
	ZHD=F1**(-2.D0/3.D0)*(1.D0+COEF*(3.D0*hb**2*LQT+(MTAU/H1)**2
     C      *DLOG(QSTSB/MZ**2)+SB2*(-3*HB**2-(MTAU/H1)**2)
     C      *LMAMT-G1Q/2.D0*LM1mu-3.D0*G2Q/2.D0*LM2mu
     C	    -3.D0*(G1Q+3.D0*G2Q)/4.D0*DLOG(QSTSB/MZ**2)
     C	    -L**2*LMUNU))

* not interesting for large MZ2, vs running (==1) :
	    
	ZS=1.D0-2.D0*COEF*(L**2*Lmu+K**2*Lnu) 

*   Higgs Vevs at QSTSB

	H2Q=H2/DSQRT(ZHU)
	H1Q=H1/DSQRT(ZHD)
	HSQ=VEVS/(DSQRT(ZS)*DSQRT(2.D0))
	
	TANBQ=H2Q/H1Q

*   Top/Bottom Yukawas at QSTSB
*   (Note: RUNMB(Q) includes QCD corrections only)
*   including electroweak contributions
	HTQ=HT*(1.D0+7.D0/(4.D0*PI)*ALSMT*LQT)**(-4.D0/7.D0)/DSQRT(F1*F2)
     C	 *(1.D0+COEF/4.D0*((-17.D0/6.D0*g1q-9.D0/2.D0*g2q+hb**2)*LQT
     C   +((3.D0*CB2-1.D0)*HB**2+2.D0*HTAU**2*CB2)*LMAMT
     C    -2.D0*L**2*Lmunu-G1Q*LM1mu-3.D0*G2Q*LM2mu))

	HBQ=RUNMB(DSQRT(QSTSB))/H1Q*F1**(-1d0/6d0)
     C	 *(1d0-3d0*HTMA**2*(1d0-DLQA)/(8d0*PI*ALSMA))**(-1d0/6d0)
     C	 *(1.D0+COEF/4.D0*((-5.D0/6.D0*g1q-9.D0/2.D0*g2q
     C   +9.D0*hb**2+2.D0*HTAU**2)*LQT
     C   +(-9.D0*SB2*HB**2-2.D0*HTAU**2*SB2)*LMAMT
     C    -2.D0*L**2*Lmunu-G1Q*LM1mu-3.D0*G2Q*LM2mu))

*   Conversion to DR_bar:
	HTQ=HTQ*(1.D0-ALSQ/(3.D0*PI)+g2q*COEF*3.D0/8.D0)
	HBQ=HBQ*(1.D0-ALSQ/(3.D0*PI)+g2q*COEF*3.D0/8.D0)
	
*   Running Top and Bottom Quark Masses

	MTOPQ=HTQ*H2Q
	MBOTQ=HBQ*H1Q

*   UMSSM Parameters at QSTSB, stored in COMMON/QNMPAR
	
	KQ=K*(1.D0+3.D0*COEF*(L**2+K**2)*LQ2)
	LQ=L*(1.D0+COEF/2.D0*(-G1Q-3.D0*G2Q+4.D0*L**2+2.D0*KQ**2
     C	  +3.D0*(HTQ**2+HBQ**2)+(MTAU/H1)**2)*LQ2)

	!PRINT*,"LQ =",LQ
	!PRINT*,"KQ =",KQ

	ALQ=AL+COEF*(G1*M1+3.D0*G2*M2+4.D0*LQ**2*AL+2.D0*KQ**2*AK
     C	  +3.D0*(HTQ**2*AU+HBQ**2*AD)
     C	  +(MTAU/H1)**2*ATAU)*LQ2
	AKQ=AK+6.D0*COEF*(LQ**2*ALQ+KQ**2*AK)*LQ2


	M3=UPARF(2)

* Calculation of the SUSY corrections to h_bot, DELMB, as in
* Carena et al., hep-ph/9912516

	DELMB=MUQ*TANBQ*(2.D0/(3.D0*PI)*ALSMZ*M3*INTEG(MSB1,MSB2,M3)
     C	  +COEF*HTQ**2*AU*INTEG(MST1,MST2,MUQ)
     C	  -COEF*G2Q*M2*(CST**2*INTEG(MST1,M2,MUQ)
     C	  +SST**2*INTEG(MST2,M2,MUQ)
     C	  +1.D0/2.D0*(CSB**2*INTEG(MSB1,M2,MUQ)
     C	  +SSB**2*INTEG(MSB2,M2,MUQ))))


	RETURN

	END



