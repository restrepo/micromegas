	SUBROUTINE UDECAY()

*   Subroutine to calculate the Higgs bosons BRs
*
*
*	CU,CD,CV,CJ,CG(i)   Reduced couplings of h1,h2,h3 (i=1,2,3) or
*			    ha (i=4) to up type fermions, down type
*			    fermions, gauge bosons, gluons and photons
*			    Note: CV(4) =0
*	CB(I)               Reduced couplings of h1,h2,h3 (i=1,2,3) or
*                           ha (i=4) to b-quarks including DELMB corrections
*
*
*	WIDTH(i) Total decay width of h1,h2,h3,ha (i=1..4)
*		 with the following branching ratios:
*	BRJJ(i)  h1,h2,h3,ha	-> gluon gluon
*	BREE(i)		  "	-> e e
*	BRMM(i)		  "	-> mu mu
*	BRLL(i)		  "	-> tau tau
*	BRSS(i)		  "	-> ss
*	BRCC(i)		  "	-> cc
*	BRBB(i)		  "	-> bb
*	BRTT(i)		  "	-> tt
*	BRWW(i)		  "	-> WW (BRWW(4)=0)
*	BRZZ(i)		  "	-> Z1Z1 (BRZZ(4)=0)
*	BRGG(i)		  "	-> gamma gamma
*	BRZG(i)		  "	-> Z1 gamma
*	BRHIGGS(i)    (i=1..4) 	-> other Higgses, including:
*	BRHAA(i)	hi      -> haha (i=1..3)
*	BRHCHC(i)	hi      -> h+h- (i=1..3)
*	BRHAZ(i)	hi      -> Z1ha  (i=1..3)
*	BRHCW(i)  h1,h2,h3      -> h+W- (i=1..3), ha -> h+W- (i=4)
*	BRHHH(i)	h2      -> h1h1, h3-> h1h1, h1h2, h2h2 (i=1..4)
*	BRAHZ(i)	ha      -> Z1hi  (i=1..3)
*	BRSUSY(i)     (i=1..4)	-> susy particles, including:
*	BRNEU(i,j,k)	        -> neutralinos j,k (i=1..4, j,k=1..6)
*	BRCHA(i,j)	        -> charginos 11, 12, 22 (i=1..4, j=1..3)
*	BRHSQ(i,j)	hi      -> uLuL, uRuR, dLdL, dRdR, t1t1, t2t2,
*			           t1t2, b1b1, b2b2, b1b2 (i=1..3, j=1..10)
*	BRASQ(i)	ha      -> t1t2, b1b2 (i=1,2)
*	BRHSL(i,j)	hi      -> lLlL, lRlR, nLnL, l1l1, l2l2, l1l2,
*			           nLtnLt, nRnR, nRtnRt (i=1..3, j=1..9)
*	BRASL      	ha      -> l1l2
*
*
*	HCWIDTH		Total decay width of the charged Higgs
*			with the following branching ratios:
*	HCBRE	        h+	-> e nu_e
*	HCBRM	             	-> mu nu_mu
*	HCBRL		        -> tau nu_tau
*	HCBRSU		        -> s u
*	HCBRBU		        -> b u
*	HCBRSC		        -> s c
*	HCBRBC		        -> b c
*	HCBRBT		        -> b t
*	HCBRWHT		        -> neutral Higgs W+, including:
*	HCBRWH(i)	        -> H1W+, H2W+, h3W+, haW+ (i=1..4)
*	HCBRSUSY	        -> susy particles,including
*	HCBRNC(i,j)             -> neutralino i chargino j (i=1..6, j=1,2)
*	HCBRSQ(i)	        -> uLdL, t1b1, t1b2, t2b1, t2b2 (i=1..5)
*	HCBRSL(i)	        -> lLnL, t1nt, t2nt (i=1..3)
*
*
*       GHHH(i)         4 couplings hi,hj,hk
*       GHAA(i)         coupling hi,ha,ha (i=1..3)
*       GHCC(i)         coupling hi,h+,h-(i=1..3)
*       GHNEUNEU(i,j,k) coupling hi,neutralino j,neutralino k (i=1..3, j,k=1..6)
*       GANEUNEU(j,k)   coupling ha,neutralino j,neutralino k (j,k=1..6)
*       GHCHACHA(i,j,k) coupling hi, chargino j, chargino k (i=1..3, j,k=1,2)
*       GACHACHA(j,k)   coupling ha, chargino j, chargino k (j,k=1,2)
*       GHCNEUCHAL(i,j) coupling h+, neutralino i, chargino_L j (i=1..6, j=1,2)
*       GHCNEUCHAR(i,j) coupling h+, neutralino i, chargino_R j (i=1..6, j=1,2)
*       HRULUL(i)       coupling hi,~u_L,~u_L (i=1..3) with running of vevs
*       HRDLDL(i)       coupling hi,~d_L,~d_L (i=1..3) with running of vevs
*       HRURUR(i)       coupling hi,~u_R,~u_R (i=1..3) with running of vevs
*       HRDRDR(i)       coupling hi,~d_R,~d_R (i=1..3) with running of vevs
*       HRLLLL(i)       coupling hi,~e_L,~e_L (i=1..3) with running of vevs
*       HRLRLR(i)       coupling hi,~e_R,~e_R (i=1..3) with running of vevs
*       HRNLNL(i)       coupling hi,~nu_L,~nu_L (i=1..3) with running of vevs
*       HRNRNR(i)       coupling hi,~nu_R,~nu_R (i=1..3) with running of vevs

	IMPLICIT NONE

	INTEGER A,B,D,I,J,K,N0,NF,NFEXT,NFGG,LOOP,VFLAG,HFLAG

	DOUBLE PRECISION MH,MA,S(3,3),SMASS(3)
	DOUBLE PRECISION PMASS,P(3,3),CMASS,C(2)
	DOUBLE PRECISION CU(4),CD(4),CV(3),CVZ(3),CJ(4),CG(4),CI(3),CB(4)
	DOUBLE PRECISION PI,HIGTOP,ASG,ASH,AS3,AS4,ASMT
	DOUBLE PRECISION EPS,FQCD,SQCD,XFAC,X,Y,RATCOUP,RAT
	DOUBLE PRECISION HJJ,HEE,HMM,HLL,HSS,HCC,HBB,HTT,HWW,HZZ,HGG,HZG
	DOUBLE PRECISION HS1,HS2,HC1,HC2,HB1,HB2,HT1,HT2,DCC,DBB
	DOUBLE PRECISION DLU,DLD,XM1,XM2,CWW,CZZ,XX(4),YY(4)
	DOUBLE PRECISION HHH(4),HAA,HHA(3),HHCHC
	DOUBLE PRECISION HAZ,HHCW,AHZ(3),HTOT,HNEU(6,6),HCHA(3)
	DOUBLE PRECISION HSQ(10),HSL(9),ASQ(2),ASL,STOT
	DOUBLE PRECISION BRJJ(4),BREE(4),BRMM(4),BRLL(4),BRSS(4),BRCC(4)
	DOUBLE PRECISION BRBB(4),BRTT(4),BRWW(3),BRZZ(3),BRGG(4),BRINV(3)
	DOUBLE PRECISION BRZG(4),BRHHH(4),BRHAA(3),BRHCHC(3)
	DOUBLE PRECISION BRHAZ(3),BRAHZ(3),BRHCW(4)
	DOUBLE PRECISION BRHIGGS(4),BRNEU(4,6,6),BRCHA(4,3)
	DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,9),BRASQ(2),BRASL
	DOUBLE PRECISION BRSUSY(4),WIDTH(4),GHHH,GHAA,GHCC,RH,CH
	DOUBLE PRECISION GHNEUNEU,GANEUNEU,GHCHACHA,GACHACHA
	DOUBLE PRECISION GHCNEUCHAL,GHCNEUCHAR
	DOUBLE PRECISION MGL,MCHA(2),U(2,2),V(2,2),MNEU(6),N(6,6)
	DOUBLE PRECISION TANBETA,SINBETA,COSBETA
	DOUBLE PRECISION XLAMBDA,MC0,MB0,MT0,RMS,RMC,RMB,RMT
	DOUBLE PRECISION SQR2,LAMBDA,ALAMBDA,H1,H2,ss
	DOUBLE PRECISION HCBRE,HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC
	DOUBLE PRECISION HCBRBT,HCBRWH(4),HCBRWHT,HCBRNC(6,2)
	DOUBLE PRECISION HCBRSQ(5),HCBRSL(3),HCBRSUSY,HCWIDTH
	DOUBLE PRECISION HEN,HMN,HLN,HSU,HSU1,HSU2,HSC,HSC1,HSC2
	DOUBLE PRECISION HBC,HBC1,HBC2,HBU,HBU1,HBU2,HBT,HBT1,HBT2
	DOUBLE PRECISION HCWH(4),HCNC(6,2),HCSQ(5),HCSL(3)
	DOUBLE PRECISION VUS,VCB,VUB,C2TW,T2TW,ALEM0
	DOUBLE PRECISION HVV,HV,HFF,QCD0,HQCDM,HQCD,QCDH,TQCDH,HGGQCD
	DOUBLE PRECISION AFF,AQCDM,AQCD,QCDA,TQCDA,AGGQCD,SGGQCD
	DOUBLE PRECISION QCDC,QCDCM,CQCD,CQCDM,QCDCI,QCDCMI
	DOUBLE PRECISION BETA,LAMB,SP,ALPHAS,RUNM,QQINT,FINT
	DOUBLE PRECISION T,Z,XI,BIJ,CFF
	DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
	DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
	DOUBLE PRECISION QSTSB,RMTTOP,FT,FB,ACOUP	
	DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
	DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
	DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU,SST,SSB,SSL
	DOUBLE PRECISION HF,HLQ,AU,AD,ATAU,LQ
	DOUBLE PRECISION HRULUL,HRDLDL,HRURUR,HRDRDR,HRT1T1,HRT2T2,HRT1T2
	DOUBLE PRECISION HRB1B1,HRB2B2,HRB1B2,HRLLLL,HRLRLR,HRNLNL,HRNRNR
	DOUBLE PRECISION HRL1L1,HRL2L2,HRL1L2    
	DOUBLE PRECISION HIT1T2,HIB1B2,HIL1L2,HPULDL,HPULDR,HPURDL,HPURDR
	DOUBLE PRECISION HPT1B1,HPT1B2,HPT2B1,HPT2B2,HPLLNL,HPL1NL,HPL2NL
	DOUBLE PRECISION HPLRNL
	DOUBLE PRECISION DELMB,RUNMB
	DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ,HSQQ
	DOUBLE PRECISION HUQ,HDQ,MTQ,MBQ
	DOUBLE PRECISION LLQ,KQ,ALQ,AKQ,MUEFFQ,NUQ
	DOUBLE PRECISION UPARF,SAZZ,CAZZ,VEV,NCP,VEVS,G1P,S2AZ
	DOUBLE PRECISION QD,QU,QS,QQ,QUP,QDOW,QL,QE,QN,M2,AL
	DOUBLE PRECISION SGNT,SGNB,SGNL
	DOUBLE PRECISION FCH1,FCH2,HTWW,HTZZ

	DOUBLE COMPLEX CTT,CTB,CTC,CTL,CTW,CTHC,CTCH1,CTCH2
	DOUBLE COMPLEX CXT,CXB,CXC,CXL,CXW,CXHC,CXCH1,CXCH2
	DOUBLE COMPLEX CXTSM,CXBSM,CXCSM,CXLSM
	DOUBLE COMPLEX CTUL,CTUR,CTDL,CTDR,CTST1,CTST2,CTSB1,CTSB2
	DOUBLE COMPLEX CXUL,CXUR,CXDL,CXDR,CXST1,CXST2,CXSB1,CXSB2
	DOUBLE COMPLEX CTLL,CTLR,CTSL1,CTSL2,CXLL,CXLR,CXSL1,CXSL2
	DOUBLE COMPLEX CLT,CLB,CLC,CLW,CLH,CXTZ,CXBZ,CXWZ,CXHZ,CXCZ
	DOUBLE COMPLEX CI1,CI2,CGZ,CF,CA,CBC
	DOUBLE COMPLEX CLCH1,CLCH2,CXCH1Z,CXCH2Z
		
	COMMON/ALEM0/ALEM0
	COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
	COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
	COMMON/CKM/VUS,VCB,VUB
	COMMON/ALS/XLAMBDA,MC0,MB0,MT0,N0
	COMMON/BRN/BRJJ,BREE,BRMM,BRLL,BRSS,BRCC,BRBB,BRTT,BRWW,BRZZ,
     C		BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHZ,
     C		BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     C		BRSUSY,WIDTH
	COMMON/BRC/HCBRE,HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC,
     C		HCBRBT,HCBRWH,HCBRWHT,HCBRNC,HCBRSQ,HCBRSL,
     C		HCBRSUSY,HCWIDTH
	COMMON/REDCOUP/CU,CD,CV,CVZ,CJ,CG,CB
	COMMON/HIGGSPEC/SMASS,S,PMASS,P,CMASS
	COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,N
	COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     C		MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     C		CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
	COMMON/DELMB/DELMB
	COMMON/STSBSCALE/QSTSB
	COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ,HSQQ
	COMMON/QQUARK/HUQ,HDQ,MTQ,MBQ
	COMMON/QNMPAR/LLQ,KQ,ALQ,AKQ,MUEFFQ,NUQ
	COMMON/UMSSM/SAZZ,CAZZ,VEV,NCP,QD,QU,QS,VEVS,G1P,QQ,
     C		QUP,QDOW,QL,QE,QN
	COMMON/NOBUG/TANBETA,AU,AD,M2,SST,SSB,LAMBDA,AL,ATAU
	COMMON/INV/BRINV
	COMMON/FLAGS/VFLAG,HFLAG

	QQINT(RAT,X,Y)= RAT**2*X+(1.D0-RAT**2)*Y
	BETA(X)= DSQRT(1.D0-4.D0*X)
	LAMB(X,Y)= DSQRT((1.D0-X-Y)**2-4.D0*X*Y)
	CF(CA)= -CDLOG(-(1.D0+CDSQRT(1.D0-CA))
     .	 / (1.D0-CDSQRT(1.D0-CA)))**2/4.D0
	CGZ(CA)= CDSQRT(1.D0-CA)/2.D0*CDLOG(-(1.D0+CDSQRT(1.D0-CA))
     .	 / (1.D0-CDSQRT(1.D0-CA)))
	CI1(CA,CBC)= CA*CBC/2.D0/(CA-CBC)
     .	 + CA**2*CBC**2/2/(CA-CBC)**2*(CF(CA)-CF(CBC))
     .	 + CA**2*CBC/(CA-CBC)**2*(CGZ(CA)-CGZ(CBC))
	CI2(CA,CBC)= -CA*CBC/2.D0/(CA-CBC)*(CF(CA)-CF(CBC))
	HV(X)= 3.D0*(1.D0-8.D0*X+20.D0*X**2)/DSQRT((4.D0*X-1.D0))
     .	 * DACOS((3.D0*X-1.D0)/2.D0/DSQRT(X**3))
     .	 - (1.D0-X)*(47.D0/2.D0*X-13.D0/2.D0+1.D0/X)
     .	 - 3.D0/2.D0*(1.D0-6.D0*X+4.D0*X**2)*DLOG(X)
	HVV(X,Y)= GF/(4.D0*PI*SQR2)*X**3/2.D0*BETA(Y)
     .	 * (1.D0-4.D0*Y+12.D0*Y**2)
	HFF(X,Y)= GF/(4.D0*PI*SQR2)*X**3*Y*(BETA(Y))**3
	AFF(X,Y)= GF/(4.D0*PI*SQR2)*X**3*Y*(BETA(Y))
	CFF(Z,T,X,Y)= GF/(4.D0*PI*SQR2)*Z**3*LAMB(X,Y)
     .	 * ((1.D0-X-Y)*(X*T**2+Y/T**2)-4.D0*X*Y)
	QCD0(X)= (1.D0+X**2)*(4.D0*SP((1.D0-X)/(1.D0+X))
     .	 +2.D0*SP((X-1.D0)/(X+1.D0))
     .	 - 3.D0*DLOG((1.D0+X)/(1.D0-X))*DLOG(2.D0/(1.D0+X))
     .	 - 2.D0*DLOG((1.D0+X)/(1.D0-X))*DLOG(X))
     .	 - 3.D0*X*DLOG(4.D0/(1.D0-X**2))-4.D0*X*DLOG(X)
	HQCDM(X)= QCD0(X)/X+(3.D0+34.D0*X**2-13.D0*X**4)/16.D0/X**3
     .	 * DLOG((1.D0+X)/(1.D0-X))+3.D0/8.D0/X**2*(7.D0*X**2-1.D0)
	AQCDM(X)= QCD0(X)/X+(19.D0+2.D0*X**2+3.D0*X**4)/16.D0/X
     .	 * DLOG((1.D0+X)/(1.D0-X))+3.D0/8.D0*(7.D0-X**2)
*	HQCD(X)= 5.67D0*ASH/PI
*     .	 + (29.14D0+RATCOUP*(1.57D0-2.D0*DLOG(HIGTOP)/3.D0
*     .	 + DLOG(X)**2/9.D0))*(ASH/PI)**2
*     .	 + (164.14D0-25.77D0*5.D0+0.259D0*5.D0**2)*(ASH/PI)**3
*	AQCD(X)= 5.67D0*ASH/PI
*     .	 + (29.14D0+RATCOUP*(3.83D0-DLOG(HIGTOP)
*     .	 + DLOG(X)**2/6.D0))*(ASH/PI)**2
*     .	 + (164.14D0-25.77D0*5.D0+0.259D0*5.D0**2)*(ASH/PI)**3
* New July 2010 (NMSSM):
      HQCD(X)=(4d0/3d0*HQCDM(BETA(X))
     .   +2d0*(4d0/3d0-DLOG(X))*(1d0-10d0*X)/(1d0-4d0*X))*ASH/PI
     .   + (29.14671d0 + RATCOUP*(1.570d0 - 2d0*DLOG(HIGTOP)/3d0
     .   + DLOG(X)**2/9d0))*(ASH/PI)**2
     .   +(164.14d0 - 25.77d0*5 + 0.259d0*5**2)*(ASH/PI)**3
      AQCD(X)=(4d0/3d0*AQCDM(BETA(X))
     .   +2d0*(4d0/3d0-DLOG(X))*(1d0-6d0*X)/(1d0-4d0*X))*ASH/PI
     .   + (29.14671d0 + RATCOUP*(23d0/6d0 - DLOG(HIGTOP)
     .   + DLOG(X)**2/6d0))*(ASH/PI)**2
     .   + (164.14d0 - 25.77d0*5 + 0.259d0*5**2)*(ASH/PI)**3
* End New
	QCDH(X)= 1.D0+HQCD(X)
	TQCDH(X)= 1.D0+4.D0/3.D0*HQCDM(BETA(X))*ASH/PI
	QCDA(X)= 1.D0+AQCD(X)
	TQCDA(X)= 1.D0+4.D0/3.D0*AQCDM(BETA(X))*ASH/PI
	HGGQCD(ASG,NF)= 1.D0+ASG/PI*(95.D0/4.D0-NF*7.D0/6.D0)
	AGGQCD(ASG,NF)= 1.D0+ASG/PI*(97.D0/4.D0-NF*7.D0/6.D0)
	SGGQCD(ASG)= ASG/PI*17.D0/6.D0
	XI(X,Y)= 2.D0*X/(1.D0-X-Y+LAMB(X,Y))
	BIJ(X,Y)= (1.D0-X-Y)/LAMB(X,Y)
     .	 * (4.D0*SP(XI(X,Y)*XI(Y,X))
     .	 - 2.D0*SP(-XI(X,Y))-2.D0*SP(-XI(Y,X))
     .	 + 2.D0*DLOG(XI(X,Y)*XI(Y,X))*DLOG(1.D0-XI(X,Y)*XI(Y,X))
     .	 - DLOG(XI(X,Y))*DLOG(1.D0+XI(X,Y))
     .	 - DLOG(XI(Y,X))*DLOG(1.D0+XI(Y,X)))
     .	 - 4.D0*(DLOG(1.D0-XI(X,Y)*XI(Y,X))
     .	 + XI(X,Y)*XI(Y,X)/(1.D0-XI(X,Y)*XI(Y,X))*DLOG(XI(X,Y)*XI(Y,X)))
     .	 + (LAMB(X,Y)+X-Y)/LAMB(X,Y)*(DLOG(1.D0+XI(X,Y))
     .	 - XI(X,Y)/(1.D0+XI(X,Y))*DLOG(XI(X,Y)))
     .	 + (LAMB(X,Y)-X+Y)/LAMB(X,Y)*(DLOG(1.D0+XI(Y,X))
     .	 - XI(Y,X)/(1.D0+XI(Y,X))*DLOG(XI(Y,X)))
	QCDC(X,Y)= 1.D0+4.D0/3.D0*ASH/PI*(9.D0/4.D0+(3.D0-2.D0*X+2.D0*Y)
     .	 / 4.D0*DLOG(X/Y)+((1.5D0-X-Y)*LAMB(X,Y)**2+5.D0*X*Y)/2.D0
     .	 / LAMB(X,Y)/(1.D0-X-Y)*DLOG(XI(X,Y)*XI(Y,X))+BIJ(X,Y))
     .	 + ASH/PI*(2.D0*(4.D0/3.D0-DLOG(X))
     .	 - (X*2.D0*(4.D0/3.D0-DLOG(X))+Y*2.D0*(4.D0/3.D0-DLOG(Y)))
     .	 / (1.D0-X-Y)-(X*2.D0*(4.D0/3.D0-DLOG(X))*(1.D0-X+Y)
     .	 + Y*2.D0*(4.D0/3.D0-DLOG(Y))*(1.D0+X-Y))/LAMB(X,Y)**2)
	QCDCI(X,Y)= 1.D0+4.D0/3.D0*ASH/PI*(3.D0+(Y-X)/2.D0*DLOG(X/Y)
     .	 + (2.D0*(1.D0-X-Y)+LAMB(X,Y)**2)/2.D0/LAMB(X,Y)
     .	 * DLOG(XI(X,Y)*XI(Y,X))+BIJ(X,Y))
     .	 + ASH/PI*(2.D0*(4.D0/3.D0-DLOG(X))+2.D0*(4.D0/3.D0-DLOG(Y))
     .	 - (X*2.D0*(4.D0/3.D0-DLOG(X))*(1.D0-X+Y)
     .	 + Y*2.D0*(4.D0/3.D0-DLOG(Y))*(1.D0+X-Y))/LAMB(X,Y)**2)
	QCDCM(X,Y)= 1.D0+4.D0/3.D0*ASH/PI*(9.D0/4.D0
     .	 + (3.D0-2.D0*X+2.D0*Y)/4.D0*DLOG(X/Y)+((1.5D0-X-Y)
     .	 * LAMB(X,Y)**2+5.D0*X*Y)/2.D0/LAMB(X,Y)/(1.D0-X-Y)
     .	 * DLOG(4.D0*X*Y/(1.D0-X-Y+LAMB(X,Y))**2)
     .	 + BIJ(X,Y))
	QCDCMI(X,Y)= 1.D0+4.D0/3.D0*ASH/PI*(3.D0+(Y-X)/2.D0*DLOG(X/Y)
     .	 + (2.D0*(1.D0-X-Y)*LAMB(X,Y)**2)/2.D0/LAMB(X,Y)
     .	 * DLOG(4.D0*X*Y/(1.D0-X-Y+LAMB(X,Y))**2)
     .	 + BIJ(X,Y))
	CQCD(Z,T,X,Y)= GF/(4.D0*PI*SQR2)*Z**3*LAMB(X,Y)
     .	 * ((1.D0-X-Y)*(X*T**2*QCDC(Y,X)
     .	 + Y/T**2*QCDC(X,Y))
     .	 - 4.D0*X*Y*QCDCI(X,Y))
	CQCDM(Z,T,X,Y)= GF/(4.D0*PI*SQR2)*Z**3*LAMB(X,Y)
     .	 * ((1.D0-X-Y)*(X*T**2*QCDCM(Y,X)
     .	 + Y/T**2*QCDCM(X,Y))
     .	 - 4.D0*X*Y*QCDCMI(X,Y))
	GHHH(A)= UPARF(147+A)/(4*S2TW*C2TW) 
	GHAA(A)= UPARF(151+A)/(4*S2TW*C2TW)
	GHCC(A)= -UPARF(154+A)/(4*S2TW*C2TW)     
	GHNEUNEU(A,B,D)= UPARF(160+A+3*(D-B)+LOOP)
     .	 /DSQRT(S2TW*C2TW)
	GANEUNEU(B,D)= UPARF(287+(D-B)+LOOP)
     .	 /DSQRT(S2TW*C2TW)
	GHCHACHA(A,B,D)= -UPARF(223+A+3*(2*B+D-3))/(SQR2*DSQRT(S2TW))
	GACHACHA(B,D)= -UPARF(308+(2*B+D-3))/(SQR2*DSQRT(S2TW))
	GHCNEUCHAL(A,B)= -UPARF(314+A+6*(B-1))/(2.D0*DSQRT(C2TW*S2TW))
	GHCNEUCHAR(A,B)= -UPARF(326+A+6*(B-1))/(2.D0*DSQRT(C2TW*S2TW))
* Relations commented : no running for the different parameters
	HRULUL(A)= -SQR2*g1/(12.D0*S2TW)*( (3.D0-4.D0*S2TW)
     .	 *(H1Q*S(A,1)-H2Q*S(A,2)) ) -SQR2*G1P**2*QQ*(QD*H1Q*S(A,1)
     .	 +QU*H2Q*S(A,2)+QS*HSQQ*S(A,3))
*	HRULUL(A)= UPARF(235+A)*g1/(12.D0*S2TW)
	HRDLDL(A)= SQR2*g1/(12.D0*S2TW)*( (3.D0-2.D0*S2TW)
     .	 *(H1Q*S(A,1)-H2Q*S(A,2)) ) -SQR2*G1P**2*QQ*(QD*H1Q*S(A,1)
     .	 +QU*H2Q*S(A,2)+QS*HSQQ*S(A,3))
*	HRDLDL(A)= UPARF(241+A)*g1/(12.D0*S2TW)
	HRURUR(A)= -SQR2*g1/3.D0*(H1Q*S(A,1)-H2Q*S(A,2))-SQR2*G1P**2
     .	 *QUP*(QD*H1Q*S(A,1)+QU*H2Q*S(A,2)+QS*HSQQ*S(A,3))
*	HRURUR(A)= -UPARF(238+A)*g1/3.D0
	HRDRDR(A)= SQR2*g1/6.D0*(H1Q*S(A,1)-H2Q*S(A,2))-SQR2*G1P**2
     .	 *QDOW*(QD*H1Q*S(A,1)+QU*H2Q*S(A,2)+QS*HSQQ*S(A,3))
*	HRDRDR(A)= UPARF(244+A)*g1/6.D0
	HRT1T1(A)= HRULUL(A)*CST**2 + HRURUR(A)*SST**2 +SQR2*HUQ
     .	 *( CST*SST*(MUEFFQ*S(A,1)-AU*S(A,2)+LLQ*H1Q*S(A,3))
     .	 - MTQ*S(A,2))
*	HRT1T1(A)= -UPARF(247+A)*DSQRT(g1)/(24.D0*SINBETA*S2TW*DSQRT(C2TW)*MW**2)
	HRT2T2(A)= HRULUL(A)*SST**2 + HRURUR(A)*CST**2 +SQR2*HUQ
     .	 *(-CST*SST*(MUEFFQ*S(A,1)-AU*S(A,2)+LLQ*H1Q*S(A,3))
     .	 - MTQ*S(A,2))
*	HRT2T2(A)= -UPARF(250+A)*DSQRT(g1)/(24.D0*SINBETA*S2TW*DSQRT(C2TW)*MW**2)
	HRT1T2(A)= SGNT*(HRULUL(A)*CST*SST - HRURUR(A)*SST*CST 
     .	 +HUQ/SQR2
     .	 *(SST**2-CST**2)*(MUEFFQ*S(A,1)-AU*S(A,2)+LLQ*H1Q*S(A,3)))
*	HRT1T2(A)= -UPARF(253+A)*DSQRT(g1)/(24.D0*SINBETA*S2TW*DSQRT(C2TW)*MW)
	HRB1B1(A)= HRDLDL(A)*CSB**2 + HRDRDR(A)*SSB**2 +SQR2*HDQ
     .	 *( CSB*SSB*(MUEFFQ*S(A,2)-AD*S(A,1)+LLQ*H2Q*S(A,3))
     .	 - MBQ*S(A,1))
*	HRB1B1(A)= -UPARF(256+A)*DSQRT(g1)/(24.D0*COSBETA*S2TW*DSQRT(C2TW)*MW**2)
	HRB2B2(A)= HRDLDL(A)*SSB**2 + HRDRDR(A)*CSB**2 +SQR2*HDQ
     .	 *(-CSB*SSB*(MUEFFQ*S(A,2)-AD*S(A,1)+LLQ*H2Q*S(A,3))
     .	 - MBQ*S(A,1))
*	HRB2B2(A)= -UPARF(259+A)*DSQRT(g1)/(24.D0*COSBETA*S2TW*DSQRT(C2TW)*MW**2)
	HRB1B2(A)= SGNB*(HRDLDL(A)*CSB*SSB - HRDRDR(A)*SSB*CSB
     .	  +HDQ/SQR2
     .	 *(SSB**2-CSB**2)*(MUEFFQ*S(A,2)-AD*S(A,1)+LLQ*H2Q*S(A,3)))
*	HRB1B2(A)= -UPARF(262+A)*DSQRT(g1)/(24.D0*COSBETA*S2TW*DSQRT(C2TW)*MW)
	HRLLLL(A)= SQR2*g1/(4.D0*S2TW)*( (1.D0-2.D0*S2TW)
     .	 *(H1Q*S(A,1)-H2Q*S(A,2)) ) -SQR2*G1P**2*QL*(QD*H1Q*S(A,1)
     .	 +QU*H2Q*S(A,2)+QS*HSQQ*S(A,3))
*	HRLLLL(A)= UPARF(265+A)*g1/(4.D0*S2TW)
	HRLRLR(A)= SQR2*g1/2.D0*(H1Q*S(A,1)-H2Q*S(A,2))-SQR2*G1P**2
     .	 *QE*(QD*H1Q*S(A,1)+QU*H2Q*S(A,2)+QS*HSQQ*S(A,3))
*	HRLRLR(A)= UPARF(268+A)*g1/2.D0
	HRL1L1(A)= HRLLLL(A)*CSL**2 + HRLRLR(A)*SSL**2 +SQR2*HLQ
     .	 *( CSL*SSL*(MUEFFQ*S(A,2)-ATAU*S(A,1)+LLQ*H2Q*S(A,3))
     .	 -MTAU*S(A,1))
*	HRL1L1(A)= UPARF(274+A)*DSQRT(g1)/(8.D0*COSBETA*S2TW*DSQRT(C2TW)*MW**2)
	HRL2L2(A)= HRLLLL(A)*SSL**2 + HRLRLR(A)*CSL**2 +SQR2*HLQ
     .	 *(-CSL*SSL*(MUEFFQ*S(A,2)-ATAU*S(A,1)+LLQ*H2Q*S(A,3))
     .	 -MTAU*S(A,1))
*	HRL2L2(A)= UPARF(277+A)*DSQRT(g1)/(8.D0*COSBETA*S2TW*DSQRT(C2TW)*MW**2)
	HRL1L2(A)= SGNL*(HRLLLL(A)*CSL*SSL - HRLRLR(A)*SSL*CSL
     .	  +HLQ/SQR2
     .	 *(SSL**2-CSL**2)*(MUEFFQ*S(A,2)-ATAU*S(A,1)+LLQ*H2Q*S(A,3)))
*	HRL1L2(A)= UPARF(280+A)*DSQRT(g1)/(8.D0*COSBETA*S2TW*DSQRT(C2TW)*MW)
	HRNLNL(A)= -SQR2*g1/(4.D0*S2TW)*(H1Q*S(A,1)-H2Q*S(A,2))
     .	 -SQR2*G1P**2*QL*(QD*H1Q*S(A,1)
     .	 +QU*H2Q*S(A,2)+QS*HSQQ*S(A,3))
*	HRNLNL(A)= -UPARF(271+A)*g1/(4.D0*S2TW)
	HRNRNR(A)= -SQR2*G1P**2*QN*(QD*H1Q*S(A,1)
     .	 +QU*H2Q*S(A,2)+QS*HSQQ*S(A,3))
*	HRNRNR(A)= -UPARF(283+A)*QN*G1P**2
	HIT1T2(A)= HUQ/SQR2*(MUEFFQ*P(A,1)+AU*P(A,2)+LLQ*H1Q*P(A,3))
*	HIT1T2(A)= -UPARF(312)*DSQRT(g2)*MT/(4.D0*SINBETA*MW)
	HIB1B2(A)= HDQ/SQR2*(MUEFFQ*P(A,2)+AD*P(A,1)+LLQ*H2Q*P(A,3))
*	HIB1B2(A)= -UPARF(313)*DSQRT(g2)*MBP/(4.D0*COSBETA*MW)
	HIL1L2(A)= HLQ/SQR2*(MUEFFQ*P(A,2)+ATAU*P(A,1)+LLQ*H2Q*P(A,3))
*	HIL1L2(A)= -UPARF(314)*DSQRT(g2)*MTAU/(4.D0*COSBETA*MW)

	EPS= 1.D-8
	PI= 4.D0*DATAN(1.D0)
	SQR2= DSQRT(2.D0)

*   Number of light flavours included in the gluonic decays
*   Higgs -> gg* -> gqq (see hdecay): NFGG = 3
	NFGG= 3


	C2TW= 1.D0-S2TW
	T2TW= S2TW/C2TW
    
	COSBETA= 1.D0/DSQRT(1.D0+TANBETA**2)
	SINBETA= TANBETA*COSBETA
        S2AZ=2*SAZZ*CAZZ

*   Higgs vevs H1 = <h_d>, H2 = <h_u>
	H2= SINBETA/DSQRT(2.D0*SQR2*GF)
	H1= H2/TANBETA
	
	ss= MUEFFQ/LAMBDA
	ALAMBDA= ALQ

*   Rotation angles in the stop/sbottom/stau sectors:
*	SST= DSQRT(1.D0-CST**2)
*	SSB= DSQRT(1.D0-CSB**2)
*	SSL= DSQRT(1.D0-CSL**2)
	SGNT= 1.D0 ! "ordinary" case here : SST>CST
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

*   CP odd/even mixing angles
	DO I=1,3
	 P(1,I)= UPARF(45+I)
	 P(2,I)= UPARF(48+I)
	 P(3,I)= UPARF(51+I)
	 S(1,I)= UPARF(55+I)
	 S(2,I)= UPARF(58+I)
	 S(3,I)= UPARF(61+I)
	 SMASS(I)= UPARF(64+I)   
	ENDDO
    
	PMASS= UPARF(55)
	CMASS= UPARF(68)
    
*   neutralino masses/mixing angles
	DO I=1,6
	 N(1,I)= UPARF(104+I)
	 N(2,I)= UPARF(110+I)
	 N(3,I)= UPARF(116+I)
	 N(4,I)= UPARF(122+I)
	 N(5,I)= UPARF(128+I)
	 N(6,I)= UPARF(134+I)
	 MNEU(I)=UPARF(140+I)
	ENDDO
    
*   chargino mixing angles
	DO I=1,2
	 U(1,I)= UPARF(68+I)
	 U(2,I)= UPARF(70+I)
	 V(1,I)= UPARF(72+I)
	 V(2,I)= UPARF(74+I)
	 MCHA(I)= UPARF(76+I)
	ENDDO

	MLL= UPARF(91)
	MLR= UPARF(92)
	MSL1= UPARF(97)
	MSL2= UPARF(98)
	MNL= UPARF(99)
	MSNT= UPARF(101)
        
	C(1)= COSBETA
	C(2)= SINBETA

*   Alpha_s at the top pole mass scales, used for the running
*   Yukawa coupling ht and running quark masses RMT below
*   NOTE: MT = top pole mass
	ASMT= ALPHAS(MT,2)
	
*   Tau Yukawa coupling for Higgs-stau coupling
	HLQ=MTAU/H1Q

*   MT = Top pole mass; RMTTOP = running mass at Mtop (MS_bar):
	RMTTOP= MT/(1.D0+4.D0*ASMT/(3.D0*PI)+11.D0*(ASMT/PI)**2)

* Loop over CP-even Higgs bosons
 
	DO I=1,3

	 MH= SMASS(I)
	 HIGTOP= (MH/MT)**2
	 MT0= 3.D8
	 ASH= ALPHAS(MH,2)
	 MC0= 1.D8
	 MB0= 2.D8
	 AS3= ALPHAS(MH,2)
	 MC0= MC
	 AS4= ALPHAS(MH,2)
	 MB0= MBP
	 MT0= MT

*  Running quark masses at MH

	 RMS= RUNM(MH,3)
         
	 RMC= RUNM(MH,4)

         RMB= RUNMB(MH)

	 IF(MH.GE.MT)THEN
	  RMT= RMTTOP
     C	   *(1.D0+7.D0/(4.D0*PI)*ASMT*DLOG(MH**2/MT**2))
     C	   **(-4.D0/7.D0)
	 ELSE
 	  RMT= RMTTOP
     C	   *(1.D0+23.D0/(12.D0*PI)*ASMT*DLOG(MH**2/MT**2))
     C	   **(-12.D0/23.D0)
	 ENDIF
 
*   Log for rad. corrs to Higgs self couplings
*   (QSTSB was initialized in the subroutine RUNPAR):
	LQ=DLOG(MAX(QSTSB,MH**2)/(MAX(MT,MH)**2))

*  Couplings relative to the standard model

	 CV(I)= S(I,1)*COSBETA+S(I,2)*SINBETA
	 CVZ(I)= CAZZ**2*(S(I,1)*COSBETA+S(I,2)*SINBETA)
     C	     + 2*DSQRT(S2TW)*NCP*S2AZ
     C	     *(COSBETA*S(I,1)*QD - SINBETA*S(I,2)*QU)
     C	     + 4*S2TW*(NCP*SAZZ)**2*(COSBETA*S(I,1)*QD**2
     C	     + SINBETA*S(I,2)*QU**2 + S(I,3)*QS**2*VEVS/VEV)
	 CU(I)= S(I,2)/SINBETA
	 CD(I)= S(I,1)/COSBETA
	 CB(I)=CD(I)*(1.D0+CU(I)/CD(I)*DELMB)/(1.D0+DELMB)
	 IF(CD(I).EQ.0.D0) CD(I)= EPS
	 IF(CB(I).EQ.0.D0) CB(I)= EPS

*  Effective coupling to 2 gluons/2 photons

	 CTT= 4.D0*(MT/MH)**2*DCMPLX(1.D0,-EPS)
	 CTB= 4.D0*(MBP/MH)**2*DCMPLX(1.D0,-EPS)
	 CTC= 4.D0*(MC/MH)**2*DCMPLX(1.D0,-EPS)
	 CTL= 4.D0*(MTAU/MH)**2*DCMPLX(1.D0,-EPS)
	 CTW= 4.D0*(MW/MH)**2*DCMPLX(1.D0,-EPS)
	 CTHC= 4.D0*(CMASS/MH)**2*DCMPLX(1.D0,-EPS)
	 CTCH1= 4.D0*(MCHA(1)/MH)**2*DCMPLX(1.D0,-EPS)
	 CTCH2= 4.D0*(MCHA(2)/MH)**2*DCMPLX(1.D0,-EPS)
	 CTUL= 4.D0*(MUL/MH)**2*DCMPLX(1.D0,-EPS)
	 CTUR= 4.D0*(MUR/MH)**2*DCMPLX(1.D0,-EPS)
	 CTDL= 4.D0*(MDL/MH)**2*DCMPLX(1.D0,-EPS)
	 CTDR= 4.D0*(MDR/MH)**2*DCMPLX(1.D0,-EPS)
	 CTST1= 4.D0*(MST1/MH)**2*DCMPLX(1.D0,-EPS)
	 CTST2= 4.D0*(MST2/MH)**2*DCMPLX(1.D0,-EPS)
	 CTSB1= 4.D0*(MSB1/MH)**2*DCMPLX(1.D0,-EPS)
	 CTSB2= 4.D0*(MSB2/MH)**2*DCMPLX(1.D0,-EPS)
	 CTLL= 4.D0*(MLL/MH)**2*DCMPLX(1.D0,-EPS)
	 CTLR= 4.D0*(MLR/MH)**2*DCMPLX(1.D0,-EPS)
	 CTSL1= 4.D0*(MSL1/MH)**2*DCMPLX(1.D0,-EPS)
	 CTSL2= 4.D0*(MSL2/MH)**2*DCMPLX(1.D0,-EPS)
	 CXT= 2.D0*CTT*(1.D0+(1.D0-CTT)*CF(CTT))
	 CXB= 2.D0*CTB*(1.D0+(1.D0-CTB)*CF(CTB))
	 CXC= 2.D0*CTC*(1.D0+(1.D0-CTC)*CF(CTC))
	 CXL= 2.D0*CTL*(1.D0+(1.D0-CTL)*CF(CTL))
	 CXW= -(2.D0+3.D0*CTW+3.D0*CTW*(2.D0-CTW)*CF(CTW))
	 CXHC= CTHC*(CTHC*CF(CTHC)-1.D0)
	 CXCH1= 2.D0*CTCH1*(1.D0+(1.D0-CTCH1)*CF(CTCH1))
	 CXCH2= 2.D0*CTCH2*(1.D0+(1.D0-CTCH2)*CF(CTCH2))
	 CXUL= CTUL*(CTUL*CF(CTUL)-1.D0)/(DSQRT(SQR2*GF)*MUL**2)
     C	       *HRULUL(I)
	 CXUR= CTUR*(CTUR*CF(CTUR)-1.D0)/(DSQRT(SQR2*GF)*MUR**2)
     C	       *HRURUR(I)
	 CXDL= CTDL*(CTDL*CF(CTDL)-1.D0)/(DSQRT(SQR2*GF)*MDL**2)
     C	       *HRDLDL(I)
	 CXDR= CTDR*(CTDR*CF(CTDR)-1.D0)/(DSQRT(SQR2*GF)*MDR**2)
     C	       *HRDRDR(I)
	 CXST1= CTST1*(CTST1*CF(CTST1)-1.D0)
     C	       /(2.D0*DSQRT(SQR2*GF)*MST1**2)*HRT1T1(I)
	 CXST2= CTST2*(CTST2*CF(CTST2)-1.D0)
     C	       /(2.D0*DSQRT(SQR2*GF)*MST2**2)*HRT2T2(I)
	 CXSB1= CTSB1*(CTSB1*CF(CTSB1)-1.D0)
     C	       /(2.D0*DSQRT(SQR2*GF)*MSB1**2)*HRB1B1(I)
	 CXSB2= CTSB2*(CTSB2*CF(CTSB2)-1.D0)
     C	       /(2.D0*DSQRT(SQR2*GF)*MSB2**2)*HRB2B2(I)
	 CXLL= CTLL*(CTLL*CF(CTLL)-1.D0)/(DSQRT(SQR2*GF)*MLL**2)
     C	      *HRLLLL(I)
	 CXLR= CTLR*(CTLR*CF(CTLR)-1.D0)/(DSQRT(SQR2*GF)*MLR**2)
     C	      *HRLRLR(I)
	 CXSL1= CTSL1*(CTSL1*CF(CTSL1)-1.D0)
     C	       /(2.D0*DSQRT(SQR2*GF)*MSL1**2)*HRL1L1(I)
	 CXSL2= CTSL2*(CTSL2*CF(CTSL2)-1.D0)
     C	       /(2.D0*DSQRT(SQR2*GF)*MSL2**2)*HRL2L2(I)

*   Here CJ and CG are the actual couplings. Later we divide
*   them by the SM couplings in order to obtain reduced couplings

	 CJ(I)= CDABS(CU(I)*(CXT+CXC)+CB(I)*CXB
     C	       +CXUL+CXUR+CXDL+CXDR+CXST1+CXST2+CXSB1+CXSB2)

	 CI(I)= DREAL(DCONJG(CU(I)*(CXT+CXC)+CB(I)*CXB
     C	       +CXUL+CXUR+CXDL+CXDR+CXST1+CXST2+CXSB1+CXSB2)
     C	      *(CXUL+CXUR+CXDL+CXDR+CXST1+CXST2+CXSB1+CXSB2))

	 CG(I)= CDABS(4.D0/3.D0*CU(I)*(CXT+CXC) + CB(I)*(CXB/3.D0+CXL)
     C	  + CV(I)*CXW + GHCC(I)/(2.D0*CMASS**2*DSQRT(SQR2*GF))*CXHC
     C	  + GHCHACHA(I,1,1)/(DSQRT(SQR2*GF)*MCHA(1))*CXCH1
     C	  + GHCHACHA(I,2,2)/(DSQRT(SQR2*GF)*MCHA(2))*CXCH2
     C	  + 4.D0/3.D0*(CXUL+CXUR+CXST1+CXST2)
     C	  + 1.D0/3.D0*(CXDL+CXDR+CXSB1+CXSB2)
     C	  + CXLL+CXLR+CXSL1+CXSL2)

*  Partial widths

*   h -> gg

	 NFEXT= 3
	 ASG= AS3
	 FQCD= HGGQCD(ASG,NFEXT)
	 SQCD= SGGQCD(ASG)
	 XFAC= MAX(0.D0,CJ(I)**2*FQCD+CI(I)*SQCD)
	 HJJ= HVV(MH,0.D0)*(ASG/PI)**2*XFAC/8.D0

*   h -> gg* -> gcc to be added to h -> cc

	 NFEXT= 4
	 ASG= AS4
	 FQCD= HGGQCD(ASG,NFEXT)
	 SQCD= SGGQCD(ASG)
	 XFAC= MAX(0.D0,CJ(I)**2*FQCD+CI(I)*SQCD)
	 DCC= HVV(MH,0.D0)*(ASG/PI)**2*XFAC/8.D0-HJJ

*   h -> gg* -> gbb to be added to h -> bb

	 NFEXT= 5
	 ASG= ASH
	 FQCD= HGGQCD(ASG,NFEXT)
	 SQCD= SGGQCD(ASG)
	 XFAC= MAX(0.D0,CJ(I)**2*FQCD+CI(I)*SQCD)
	 DBB= HVV(MH,0.D0)*(ASG/PI)**2*XFAC/8.D0-HJJ-DCC

	 IF(NFGG.EQ.5)THEN
	  HJJ= HJJ+DBB+DCC
	  DBB= 0.D0
	  DCC= 0.D0
	 ELSEIF(NFGG.EQ.4)THEN
	  HJJ= HJJ+DCC
	  DCC= 0.D0
	 ENDIF

*    Below CJ becomes the REDUCED coupling
*    to 2 gluons as compared to the SM coupling

	 NFEXT= 3
	 ASG= AS3
	 FQCD= HGGQCD(ASG,NFEXT)
	 SQCD= SGGQCD(ASG)
	 CJ(I)= DSQRT(MAX(0.D0,CJ(I)**2+CI(I)*SQCD/FQCD))/
     C	  CDABS(CXT+CXC+CXB)
     
*   h -> ee

	 IF(MH.LE.2.D0*0.000511)THEN
	  HEE= 0.D0
	 ELSE
	  HEE= CD(I)**2*HFF(MH,(0.000511/MH)**2)
	 ENDIF
     
*   h -> mumu

	 IF(MH.LE.2.D0*MMUON)THEN
	  HMM= 0.D0
	 ELSE
	  HMM= CD(I)**2*HFF(MH,(MMUON/MH)**2)
	 ENDIF

*   h -> tautau

	 IF(MH.LE.2.D0*MTAU)THEN
	  HLL= 0.D0
	 ELSE
	  HLL= CD(I)**2*HFF(MH,(MTAU/MH)**2)
	 ENDIF

*   h -> ss

	 RATCOUP= CU(I)/CD(I)
	 IF(MH.LE.2.D0*MS)THEN
	  HSS= 0.D0
	 ELSE
	  HS1= 3.D0*HFF(MH,(MS/MH)**2)*CD(I)**2
     .	     * TQCDH((MS/MH)**2)
	  HS2= 3.D0*HFF(MH,(RMS/MH)**2)*CD(I)**2
     .	     * QCDH((RMS/MH)**2)
	  IF(HS2.LT.0.D0) HS2=0.D0
	  RAT= 2.D0*MS/MH
	  HSS= QQINT(RAT,HS1,HS2)
	 ENDIF

*   h -> cc

	 RATCOUP= 1.D0
	 IF(MH.LE.2.D0*MC)THEN
	  HCC= 0.D0
	 ELSE
	  HC1= 3.D0*HFF(MH,(MC/MH)**2)*CU(I)**2
     .	     * TQCDH((MC/MH)**2)
	  HC2= 3.D0*HFF(MH,(RMC/MH)**2)*CU(I)**2
     .	     * QCDH((RMC/MH)**2)
     .     + DCC
	  IF(HC2.LT.0.D0) HC2=0.D0
	  RAT= 2.D0*MC/MH
	  HCC= QQINT(RAT,HC1,HC2)
	 ENDIF

*   h -> bb

	 RATCOUP= CU(I)/CD(I)
	 IF(MH.LE.2.D0*MBP)THEN
	  HBB= 0.D0
	 ELSE
	  HB1= 3.D0*HFF(MH,(MBP/MH)**2)*CB(I)**2
     .	     * TQCDH((MBP/MH)**2)
	  HB2= 3.D0*HFF(MH,(RMB/MH)**2)*CB(I)**2
     .	     * QCDH((RMB/MH)**2)
     .     + DBB
	  IF(HB2.LT.0.D0) HB2=0.D0
	  RAT= 2.D0*MBP/MH
	  HBB= QQINT(RAT,HB1,HB2)
	 ENDIF

*   h -> tt

	 RATCOUP= 0.D0
	 IF (MH.LE.2.D0*MT)THEN
	  HTT= 0.D0
	 ELSE
	  HT1= 3.D0*HFF(MH,(MT/MH)**2)*CU(I)**2
     .	     * TQCDH((MT/MH)**2)
	  IF (MH.LE.2.D0*RMT)THEN
	   HT2= 0.D0
	  ELSE
	   HT2= 3.D0*HFF(MH,(RMT/MH)**2)*CU(I)**2
     .	      * QCDH((RMT/MH)**2)
	  ENDIF
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
	   HWW= CV(I)**2*HV((MW/MH)**2)*CWW*MH
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
	   HWW= CV(I)**2*FINT(MH,XX,YY)
	  ELSE
	   HWW= CV(I)**2*HVV(MH,(MW/MH)**2)
	  ENDIF
	 ELSE
	  CALL HTOVV(MW,2.08856D0,MH,HTWW)
	  HWW= 3.D0/2.D0*GF*MW**4/DSQRT(2.D0)/PI/MH**3*HTWW*CV(I)**2
	 ENDIF

*   h -> Z1Z1

	 IF(VFLAG.EQ.0)THEN
	  DLD= 2.D0
	  DLU= 2.D0
	  XM1= 2.D0*MZ-DLD
	  XM2= 2.D0*MZ+DLU
	  IF(MH.LE.MZ)THEN
	   HZZ= 0.D0
	  ELSEIF(MH.LE.XM1)THEN
	   CZZ= 3.D0*GF**2*MZ**4/192.D0/PI**3
     C	  * (7.D0-40.D0/3.D0*S2TW+160.D0/9.D0*S2TW**2)
	   HZZ= CVZ(I)**2*HV((MZ/MH)**2)*CZZ*MH
	  ELSEIF(MH.LT.XM2)THEN
	   CZZ= 3.D0*GF**2*MZ**4/192.D0/PI**3
     C	  * (7.D0-40.D0/3.D0*S2TW+160.D0/9.D0*S2TW**2)
	   XX(1)= XM1-1.D0
	   XX(2)= XM1
	   XX(3)= XM2
	   XX(4)= XM2+1.D0
	   YY(1)= HV((MZ/XX(1))**2)*CZZ*XX(1)
	   YY(2)= HV((MZ/XX(2))**2)*CZZ*XX(2)
	   YY(3)= HVV(XX(3),(MZ/XX(3))**2)/2.D0
	   YY(4)= HVV(XX(4),(MZ/XX(4))**2)/2.D0
	   HZZ= CVZ(I)**2*FINT(MH,XX,YY)
	  ELSE
	   HZZ= CVZ(I)**2*HVV(MH,(MZ/MH)**2)/2.D0
	  ENDIF
	 ELSE
	  CALL HTOVV(MZ,2.49581D0,MH,HTZZ)
	  HZZ= 3.D0/4.D0*GF*MZ**4/DSQRT(2.D0)/PI/MH**3*HTZZ*CVZ(I)**2
	 ENDIF
 
*   h -> gamma gamma

	 XFAC= CG(I)**2
	 HGG= HVV(MH,0.D0)*(ALEM0/PI)**2/16.D0*XFAC

*    Below CG becomes the REDUCDED coupling
*    to 2 photons as compared to the SM coupling

	 CG(I)= CG(I)/CDABS(4.D0/3.D0*(CXT+CXC)+(CXB/3.D0+CXL)+CXW)
	 
*  h -> Z1 gamma

	IF(MH.LE.MZ)THEN
	 HZG= 0.D0
	ELSE
	 FT= -2.D0*(CAZZ*(1.D0-8.D0/3.D0*S2TW)
     .  +2.D0*SAZZ*(DSQRT(S2TW)*NCP*(QQ-QUP)))/DSQRT(S2TW*C2TW)*CU(I)
	 FB= (CAZZ*(-1.D0+4.D0/3.D0*S2TW)
     .  +2.D0*SAZZ*(DSQRT(S2TW)*NCP*(QQ-QDOW)))/DSQRT(S2TW*C2TW)*CB(I)
	 FCH1=4D0*MW/(MCHA(1)*DSQRT(G2*S2TW*C2TW))*GHCHACHA(I,1,1)
     .  *CAZZ*(-V(1,1)**2-0.5D0*V(1,2)**2
     .    -U(1,1)**2-0.5D0*U(1,2)**2+2.D0*S2TW)
     .  +SAZZ*DSQRT(S2TW)*NCP*(QD*U(1,2)**2-QU*V(1,2)**2)
	 FCH2=4D0*MW/(MCHA(2)*DSQRT(G2*S2TW*C2TW))*GHCHACHA(I,2,2)
     .  *CAZZ*(-V(2,1)**2-0.5D0*V(2,2)**2
     .    -U(2,1)**2-0.5D0*U(2,2)**2+2.D0*S2TW)
     .  +SAZZ*DSQRT(S2TW)*NCP*(QD*U(2,2)**2-QU*V(2,2)**2)
	 CLT= 4.D0*(MT/MZ)**2*DCMPLX(1.D0,-EPS)
	 CLB= 4.D0*(MBP/MZ)**2*DCMPLX(1.D0,-EPS)
	 CLC= 4.D0*(MC/MZ)**2*DCMPLX(1.D0,-EPS)
	 CLW= 4.D0*(MW/MZ)**2*DCMPLX(1.D0,-EPS)
	 CLH= 4.D0*(CMASS/MZ)**2*DCMPLX(1.D0,-EPS)
	 CLCH1= 4d0*(MCHA(1)/MZ)**2*DCMPLX(1d0,-EPS)
	 CLCH2= 4d0*(MCHA(2)/MZ)**2*DCMPLX(1d0,-EPS)
	 CXTZ= FT*(CI1(CTT,CLT) - CI2(CTT,CLT))
	 CXBZ= FB*(CI1(CTB,CLB) - CI2(CTB,CLB))
	 CXCZ= FT*(CI1(CTC,CLC) - CI2(CTC,CLC))
	 CXWZ= -1.D0/DSQRT(T2TW)*(4.D0*(3.D0-T2TW)*CI2(CTW,CLW)
     .	     + ((1.D0+2.D0/CTW)*T2TW
     .	     - (5.D0+2.D0/CTW))*CI1(CTW,CLW))*CV(I)
	 CXHZ= (1.D0-2.D0*S2TW)/(DSQRT(S2TW*C2TW)*2.D0*CMASS**2)
     .     * CI1(CTHC,CLH)*GHCC(I)/(DSQRT(SQR2*GF))
	 CXCH1Z=FCH1*(CI1(CTCH1,CLCH1) - CI2(CTCH1,CLCH1))
	 CXCH2Z=FCH2*(CI1(CTCH2,CLCH2) - CI2(CTCH2,CLCH2))
	 XFAC= CDABS(CXTZ+CXBZ+CXCZ+CXWZ+CXHZ+CXCH1Z+CXCH2Z)**2
	 ACOUP= SQR2*GF*MZ**2*S2TW*C2TW/PI**2
	 HZG= GF/(4.D0*PI*SQR2)*MH**3*(ALEM0/PI)*ACOUP/16.D0
     .	    * XFAC*(1.D0-(MZ/MH)**2)**3
	ENDIF

*   h -> hh

	 HTOT= 0.D0

	 IF(I.EQ.2)THEN
	  IF(MH.LE.2.D0*SMASS(1))THEN
	   HHH(1)= 0.D0
	  ELSE
	   RH= GHHH(1)
	   CH= BETA((SMASS(1)/MH)**2)
	   HHH(1)= CH*RH**2/(32.D0*PI*MH)
	   HTOT= HTOT+HHH(1)
	  ENDIF
	 ENDIF

	 IF(I.EQ.3)THEN
	  IF(MH.LE.2.D0*SMASS(1))THEN
	   HHH(2)= 0.D0
	  ELSE
	   RH= GHHH(2)
	   CH= BETA((SMASS(1)/MH)**2)
	   HHH(2)= CH*RH**2/(32.D0*PI*MH)
	   HTOT= HTOT+HHH(2)
	  ENDIF

	  IF(MH.LE.SMASS(1)+SMASS(2))THEN
	   HHH(3)= 0.D0
	  ELSE
	   RH= GHHH(3)
	   CH= LAMB((SMASS(1)/MH)**2,(SMASS(2)/MH)**2)
	   HHH(3)= CH*RH**2/(16.D0*PI*MH)
	   HTOT= HTOT+HHH(3)
	  ENDIF

	  IF(MH.LE.2.D0*SMASS(2))THEN
	   HHH(4)= 0.D0
	  ELSE
	   RH= GHHH(4)
	   CH= BETA((SMASS(2)/MH)**2)
	   HHH(4)= CH*RH**2/(32.D0*PI*MH)
	   HTOT= HTOT+HHH(4)
	  ENDIF
	 ENDIF

*   h -> aa

	 IF(MH.LE.2.D0*PMASS)THEN
	  HAA= 0.D0
	 ELSE
	  RH= GHAA(I)
	  CH= BETA((PMASS/MH)**2)
	  HAA= CH*RH**2/(32.D0*PI*MH)
	  HTOT= HTOT+HAA
	 ENDIF

*   h -> h+h-

	 IF(MH.LE.2.D0*CMASS)THEN
	  HHCHC= 0.D0
	 ELSE
	  RH= GHCC(I)
	  CH= BETA((CMASS/MH)**2)
	  HHCHC= CH*RH**2/(16.D0*PI*MH)
	  HTOT= HTOT+HHCHC
	 ENDIF

*   h -> aZ1

	  IF(MH.LE.PMASS+MZ)THEN
	   HAZ= 0.D0
	  ELSE
	   RH= UPARF(157+I)
	   CH= LAMB((PMASS/MH)**2,(MZ/MH)**2)
     C	     * LAMB((MH/MZ)**2,(PMASS/MZ)**2)**2
	   HAZ= GF/8.D0/SQR2/PI*MZ**4/MH*CH*RH**2
	   HTOT= HTOT+HAZ
	  ENDIF

*   h -> h+W-

	 IF(MH.LE.CMASS+MW)THEN
	  HHCW= 0.D0
	 ELSE
	  RH= S(I,1)*SINBETA-S(I,2)*COSBETA
	  CH= LAMB((CMASS/MH)**2,(MW/MH)**2)
     C	    * LAMB((MH/MW)**2,(CMASS/MW)**2)**2
	  HHCW= GF/8.D0/SQR2/PI*MW**4/MH*CH*RH**2
	  HTOT= HTOT+2.D0*HHCW
	 ENDIF

*   h -> neutralinos/charginos
	 STOT= 0.D0

	 LOOP= 0
	 DO J=1,6
	 DO K=J,6
	  IF(MH.LE.DABS(MNEU(J))+DABS(MNEU(K)))THEN
	   HNEU(J,K)= 0.D0
	  ELSE
	   HNEU(J,K)= 1.D0/(16.D0*PI)*(MH**2-(MNEU(J)+MNEU(K))**2)/MH
     C	        * LAMB((MNEU(J)/MH)**2,(MNEU(K)/MH)**2)
     C	        * GHNEUNEU(I,J,K)**2
	   IF(J.NE.K)HNEU(J,K)= HNEU(J,K)/2.D0
	   STOT= STOT+HNEU(J,K)
	  ENDIF
	  HNEU(K,J)= HNEU(J,K)
	 ENDDO
	 LOOP= LOOP+3*(7-J)     
	 ENDDO

	 IF(MH.LE.2.D0*DABS(MCHA(1)))THEN
	  HCHA(1)= 0.D0
	 ELSE
	  HCHA(1)= 1.D0/(8.D0*PI)*MH
     C	        * LAMB((MCHA(1)/MH)**2,(MCHA(1)/MH)**2)**3
     C	        * GHCHACHA(I,1,1)**2
	  STOT= STOT+HCHA(1)
	 ENDIF

	 IF(MH.LE.DABS(MCHA(1))+DABS(MCHA(2)))THEN
	  HCHA(2)= 0.D0
	 ELSE
	  HCHA(2)= 1.D0/(16.D0*PI)*MH
     C	        * LAMB((MCHA(1)/MH)**2,(MCHA(2)/MH)**2)
     C	        * ((GHCHACHA(I,1,2)**2+GHCHACHA(I,2,1)**2)
     C	        * (1.D0-(MCHA(1)/MH)**2-(MCHA(2)/MH)**2)
     C	        - 4.D0*GHCHACHA(I,1,2)*GHCHACHA(I,2,1)
     C	        * MCHA(1)*MCHA(2)/MH**2)
	  STOT= STOT+2.D0*HCHA(2)
	 ENDIF

	 IF(MH.LE.2.D0*DABS(MCHA(2)))THEN
	  HCHA(3)= 0.D0
	 ELSE
	  HCHA(3)= 1.D0/(8.D0*PI)*MH
     C	        * LAMB((MCHA(2)/MH)**2,(MCHA(2)/MH)**2)**3
     C	        * GHCHACHA(I,2,2)**2
	  STOT= STOT+HCHA(3)
	 ENDIF

*   h -> squarks
*   UL, UR, DL, DR are the first two families

	 IF(MH.LE.2.D0*MUL)THEN
	  HSQ(1)= 0.D0
	 ELSE
	  RH= HRULUL(I)
	  CH= BETA((MUL/MH)**2)
	  HSQ(1)= 3.D0*CH*RH**2/(16.D0*PI*MH)
	  STOT= STOT+2.D0*HSQ(1)
	 ENDIF

	 IF(MH.LE.2.D0*MUR)THEN
	  HSQ(2)= 0.D0
	 ELSE
	  RH= HRURUR(I)
	  CH= BETA((MUR/MH)**2)
	  HSQ(2)= 3.D0*CH*RH**2/(16.D0*PI*MH)
	  STOT= STOT+2.D0*HSQ(2)
	 ENDIF

	 IF(MH.LE.2.D0*MDL)THEN
	  HSQ(3)= 0.D0
	 ELSE
	  RH= HRDLDL(I)
	  CH= BETA((MDL/MH)**2)
	  HSQ(3)= 3.D0*CH*RH**2/(16.D0*PI*MH)
	  STOT= STOT+2.D0*HSQ(3)
	 ENDIF

	 IF(MH.LE.2.D0*MDR)THEN
	  HSQ(4)= 0.D0
	 ELSE
	  RH= HRDRDR(I)
	  CH= BETA((MDR/MH)**2)
	  HSQ(4)= 3.D0*CH*RH**2/(16.D0*PI*MH)
	  STOT= STOT+2.D0*HSQ(4)
	 ENDIF

	 IF(MH.LE.2.D0*MST1)THEN
	  HSQ(5)= 0.D0
	 ELSE
	  RH= HRT1T1(I)
	  CH= BETA((MST1/MH)**2)
	  HSQ(5)= 3.D0*CH*RH**2/(16.D0*PI*MH)
	  STOT= STOT+HSQ(5)
	 ENDIF

	 IF(MH.LE.2.D0*MST2)THEN
	  HSQ(6)= 0.D0
	 ELSE
	  RH= HRT2T2(I)
	  CH= BETA((MST2/MH)**2)
	  HSQ(6)= 3.D0*CH*RH**2/(16.D0*PI*MH)
	  STOT= STOT+HSQ(6)
	 ENDIF

	 IF(MH.LE.MST1+MST2)THEN
	  HSQ(7)= 0.D0
	 ELSE
	  RH= HRT1T2(I)
	  CH= LAMB((MST1/MH)**2,(MST2/MH)**2)
	  HSQ(7)= 3.D0*CH*RH**2/(16.D0*PI*MH)
	  STOT= STOT+2.D0*HSQ(7)
	 ENDIF

	 IF(MH.LE.2.D0*MSB1)THEN
	  HSQ(8)= 0.D0
	 ELSE
	  RH= HRB1B1(I)
	  CH= BETA((MSB1/MH)**2)
	  HSQ(8)= 3.D0*CH*RH**2/(16.D0*PI*MH)
	  STOT= STOT+HSQ(8)
	 ENDIF

	 IF(MH.LE.2.D0*MSB2)THEN
	  HSQ(9)= 0.D0
	 ELSE
	  RH= HRB2B2(I)
	  CH= BETA((MSB2/MH)**2)
	  HSQ(9)= 3.D0*CH*RH**2/(16.D0*PI*MH)
	  STOT= STOT+HSQ(9)
	 ENDIF

	 IF(MH.LE.MSB1+MSB2)THEN
	  HSQ(10)= 0.D0
	 ELSE
	  RH= HRB1B2(I)
	  CH= LAMB((MSB1/MH)**2,(MSB2/MH)**2)
	  HSQ(10)= 3.D0*CH*RH**2/(16.D0*PI*MH)
	  STOT= STOT+2.D0*HSQ(10)
	 ENDIF

*   h -> sleptons
*   LL, LR, NL are the first two families

	 IF(MH.LE.2.D0*MLL)THEN
	  HSL(1)= 0.D0
	 ELSE
	  RH= HRLLLL(I)
	  CH= BETA((MLL/MH)**2)
	  HSL(1)= CH*RH**2/(16.D0*PI*MH)
	  STOT= STOT+2.D0*HSL(1)
	 ENDIF

	 IF(MH.LE.2.D0*MLR)THEN
	  HSL(2)= 0.D0
	 ELSE
	  RH= HRLRLR(I)
	  CH= BETA((MLR/MH)**2)
	  HSL(2)= CH*RH**2/(16.D0*PI*MH)
	  STOT= STOT+2.D0*HSL(2)
	 ENDIF

	 IF(MH.LE.2.D0*MNL)THEN
	  HSL(3)= 0.D0
	 ELSE
	  RH= HRNLNL(I)
	  CH= BETA((MNL/MH)**2)
	  HSL(3)= CH*RH**2/(16.D0*PI*MH)
	  STOT= STOT+2.D0*HSL(3)
	 ENDIF

	 IF(MH.LE.2.D0*MSL1)THEN
	  HSL(4)= 0.D0
	 ELSE
	  RH= HRL1L1(I)
	  CH= BETA((MSL1/MH)**2)
	  HSL(4)= CH*RH**2/(16.D0*PI*MH)
	  STOT= STOT+HSL(4)
	 ENDIF

	 IF(MH.LE.2.D0*MSL2)THEN
	  HSL(5)= 0.D0
	 ELSE
	  RH= HRL2L2(I)
	  CH= BETA((MSL2/MH)**2)
	  HSL(5)= CH*RH**2/(16.D0*PI*MH)
	  STOT= STOT+HSL(5)
	 ENDIF

	 IF(MH.LE.MSL1+MSL2)THEN
	  HSL(6)= 0.D0
	 ELSE
	  RH= HRL1L2(I)
	  CH= LAMB((MSL1/MH)**2,(MSL2/MH)**2)
	  HSL(6)= CH*RH**2/(16.D0*PI*MH)
	  STOT= STOT+2.D0*HSL(6)
	 ENDIF

	 IF(MH.LE.2.D0*MSNT)THEN
	  HSL(7)= 0.D0
	 ELSE
	  RH= HRNLNL(I)
	  CH= BETA((MSNT/MH)**2)
	  HSL(7)= CH*RH**2/(16.D0*PI*MH)
	  STOT= STOT+HSL(7)
	 ENDIF

*   h -> RH sneutrinos

	 IF(MH.LE.2.D0*UPARF(102))THEN
	  HSL(8)= 0.D0
	 ELSE
	  RH= HRNRNR(I)
	  CH= BETA((UPARF(102)/MH)**2)
	  HSL(8)= CH*RH**2/(16.D0*PI*MH)
	  STOT= STOT+2.D0*HSL(8)
	 ENDIF

	 IF(MH.LE.2.D0*UPARF(104))THEN
	  HSL(9)= 0.D0
	 ELSE
	  RH= HRNRNR(I)
	  CH= BETA((UPARF(104)/MH)**2)
	  HSL(9)= CH*RH**2/(16.D0*PI*MH)
	  STOT= STOT+HSL(9)
	 ENDIF

     
*  Branching ratios

	 WIDTH(I)= HJJ+HEE+HMM+HLL+HSS+HCC+HBB
     .	     +HTT+HWW+HZZ+HGG+HZG+HTOT+STOT
	 BRJJ(I)= HJJ/WIDTH(I)
	 BREE(I)= HEE/WIDTH(I)
	 BRMM(I)= HMM/WIDTH(I)
	 BRLL(I)= HLL/WIDTH(I)
	 BRSS(I)= HSS/WIDTH(I)
	 BRCC(I)= HCC/WIDTH(I)
	 BRBB(I)= HBB/WIDTH(I)
	 BRTT(I)= HTT/WIDTH(I)
	 BRWW(I)= HWW/WIDTH(I)
	 BRZZ(I)= HZZ/WIDTH(I)
	 BRGG(I)= HGG/WIDTH(I)
	 BRZG(I)= HZG/WIDTH(I)
	 BRHIGGS(I)= HTOT/WIDTH(I)
	 BRSUSY(I)= STOT/WIDTH(I)
	 BRHAA(I)=HAA/WIDTH(I)
	 BRHCHC(I)= HHCHC/WIDTH(I)
	 BRHAZ(I)=HAZ/WIDTH(I)
	 BRHCW(I)= HHCW/WIDTH(I)
	 DO J=1,6
	 DO K=1,6
	  BRNEU(I,J,K)=HNEU(J,K)/WIDTH(I)
	 ENDDO
	 ENDDO
	 DO J=1,3
	  BRCHA(I,J)=HCHA(J)/WIDTH(I)
	 ENDDO
	 DO J=1,10
	  BRHSQ(I,J)=HSQ(J)/WIDTH(I)
	 ENDDO
	 DO J=1,9
	  BRHSL(I,J)=HSL(J)/WIDTH(I)
	 ENDDO

	 IF(MNEU(1).LE.UPARF(104))THEN
	  BRINV(I)= BRNEU(I,1,1)
	 ELSE
	  BRINV(I)= HSL(9)/WIDTH(I)
	 ENDIF

	ENDDO

	BRHHH(1)=HHH(1)/WIDTH(2)
	BRHHH(2)=HHH(2)/WIDTH(3)
	BRHHH(3)=HHH(3)/WIDTH(3)
	BRHHH(4)=HHH(4)/WIDTH(3)





* CP-odd Higgs boson

	 MA= PMASS
	 HIGTOP= (MA/MT)**2
	 MT0= 3.D8
	 ASH= ALPHAS(MA,2)
	 MC0= 1.D8
	 MB0= 2.D8
	 AS3= ALPHAS(MA,2)
	 MC0= MC
	 AS4= ALPHAS(MA,2)
	 MB0= MBP
	 MT0= MT

*  Running quark masses at MA

	 RMS= RUNM(MA,3)

	 RMC= RUNM(MA,4)

         RMB= RUNMB(MA)

	 IF(MA.GE.MT)THEN
	  RMT= RMTTOP
     C	   *(1.D0+7.D0/(4.D0*PI)*ASMT*DLOG(MA**2/MT**2))
     C	   **(-4.D0/7.D0)
	 ELSE
 	  RMT= RMTTOP
     C	   *(1.D0+23.D0/(12.D0*PI)*ASMT*DLOG(MA**2/MT**2))
     C	   **(-12.D0/23.D0)
	 ENDIF
 
*   Log for rad. corrs to Higgs self couplings
*   (QSTSB was initialized in the subroutine RUNPAR):
	LQ=DLOG(MAX(QSTSB,MA**2)/(MAX(MT,MA)**2))
	 
*  Relative couplings compared with equivalent SM

	 CU(4)= P(3,2)/SINBETA
	 CD(4)= P(3,1)/COSBETA
	 CB(4)= CD(4)/(1.D0+DELMB)
	 IF(CD(4).EQ.0.D0) CD(4)= EPS
	 IF(CB(4).EQ.0.D0) CB(4)= EPS

*  Effective coupling to 2 gluons/2 photons

	 CTT= 4.D0*(MT/MA)**2*DCMPLX(1.D0,-EPS)
	 CTB= 4.D0*(MBP/MA)**2*DCMPLX(1.D0,-EPS)
	 CTC= 4.D0*(MC/MA)**2*DCMPLX(1.D0,-EPS)
	 CTL= 4.D0*(MTAU/MA)**2*DCMPLX(1.D0,-EPS)
	 CTCH1= 4.D0*(MCHA(1)/MA)**2*DCMPLX(1.D0,-EPS)
	 CTCH2= 4.D0*(MCHA(2)/MA)**2*DCMPLX(1.D0,-EPS)
	 CXT= CTT*CF(CTT)
	 CXB= CTB*CF(CTB)
	 CXC= CTC*CF(CTC)
	 CXL= CTL*CF(CTL)
	 CXCH1= CTCH1*CF(CTCH1)
	 CXCH2= CTCH2*CF(CTCH2)

*   Here CJ and CG are the actual couplings. Later we divide
*   them by the SM couplings in order to obtain the reduced couplings 

	 CJ(4)= CDABS(CU(4)*(CXT+CXC) + CB(4)*CXB)

	 CG(4)= CDABS(4.D0/3.D0*CU(4)*(CXT+CXC)
     C	  + CB(4)*(CXB/3.D0+CXL)
     C	  + GACHACHA(1,1)/(DSQRT(SQR2*GF)*MCHA(1))*CXCH1
     C	  + GACHACHA(2,2)/(DSQRT(SQR2*GF)*MCHA(2))*CXCH2)

*  Partial widths

*   a -> gg    

	 NFEXT= 3
	 ASG= AS3
	 FQCD= AGGQCD(ASG,NFEXT)
	 XFAC= CJ(4)**2*FQCD
	 HJJ= GF/(16.D0*PI*SQR2)*MA**3*(ASG/PI)**2*XFAC

*   a -> gg* -> gcc to be added to a -> cc

	 NFEXT= 4
	 ASG= AS4
	 FQCD= AGGQCD(ASG,NFEXT)
	 XFAC= CJ(4)**2*FQCD
	 DCC= GF/(16.D0*PI*SQR2)*MA**3*(ASG/PI)**2*XFAC-HJJ

*   a -> gg* -> gbb to be added to a -> bb

	 NFEXT= 5
	 ASG= ASH
	 FQCD= AGGQCD(ASG,NFEXT)
	 XFAC= CJ(4)**2*FQCD
	 DBB= GF/(16.D0*PI*SQR2)*MA**3*(ASG/PI)**2*XFAC-HJJ-DCC

	 IF(NFGG.EQ.5)THEN
	  HJJ= HJJ+DBB+DCC
	  DBB= 0.D0
	  DCC= 0.D0
	 ELSEIF(NFGG.EQ.4)THEN
	  HJJ= HJJ+DCC
	  DCC= 0.D0
	 ENDIF

*    Below CJ becomes the REDUCED coupling
*    to 2 gluons as compared to the SM coupling

	 CXTSM= CTT*(1.D0+(1.D0-CTT)*CF(CTT))
	 CXBSM= CTB*(1.D0+(1.D0-CTB)*CF(CTB))
	 CXCSM= CTC*(1.D0+(1.D0-CTC)*CF(CTC))
	 CXLSM= 2.D0*CTL*(1.D0+(1.D0-CTL)*CF(CTL))
	 CJ(4)= CJ(4)/CDABS(CXTSM+CXCSM+CXBSM)

*   a -> ee

	 IF(MA.LE.2*0.000511)THEN
	 HEE= 0.D0
	 ELSE
	  HEE= CD(4)**2*AFF(MA,(0.000511/MA)**2)
	 ENDIF

*   a -> mumu

	 IF(MA.LE.2*MMUON)THEN
	 HMM= 0.D0
	 ELSE
	  HMM= CD(4)**2*AFF(MA,(MMUON/MA)**2)
	 ENDIF

*   a -> tautau

	 IF(MA.LE.2*MTAU)THEN
	  HLL= 0.D0
	 ELSE
	  HLL= CD(4)**2*AFF(MA,(MTAU/MA)**2)
	 ENDIF

*   a -> ss

	 RATCOUP= CU(4)/CD(4)
	 IF(MA.LE.2.D0*MS)THEN
	  HSS= 0.D0
	 ELSE
	  HS1= 3.D0*AFF(MA,(MS/MA)**2)*CD(4)**2
     .	     * TQCDA((MS/MA)**2)
	  HS2= 3.D0*AFF(MA,(RMS/MA)**2)*CD(4)**2
     .	     * QCDA((RMS/MA)**2)
	  IF(HS2.LT.0.D0) HS2=0.D0
	  RAT= 2.D0*MS/MA
	  HSS= QQINT(RAT,HS1,HS2)
	 ENDIF

*   a -> cc

	 RATCOUP= 1.D0
	 IF(MA.LE.2*MC)THEN
	  HCC= 0.D0
	 ELSE
	  HC1= 3.D0*AFF(MA,(MC/MA)**2)*CU(4)**2
     .	     * TQCDA((MC/MA)**2)
	  HC2= 3.D0*AFF(MA,(RMC/MA)**2)*CU(4)**2
     .	     * QCDA((RMC/MA)**2)
     .	     + DCC
	  IF(HC2.LT.0.D0) HC2=0.D0
	  RAT= 2.D0*MC/MA
	  HCC= QQINT(RAT,HC1,HC2)
	 ENDIF

*   a -> bb

	 RATCOUP= CU(4)/CD(4)
	 IF(MA.LE.2*MBP)THEN
	  HBB= 0.D0
	 ELSE
	  HB1= 3.D0*AFF(MA,(MBP/MA)**2)*CB(4)**2
     .	     * TQCDA((MBP/MA)**2)
	  HB2= 3.D0*AFF(MA,(RMB/MA)**2)*CB(4)**2
     .	     * QCDA((RMB/MA)**2)
     .	     + DBB
	  IF(HB2.LT.0.D0) HB2=0.D0
	  RAT= 2.D0*MBP/MA
	  HBB= QQINT(RAT,HB1,HB2)
	 ENDIF

*   a -> tt
 
	 RATCOUP= 0.D0
	 IF (MA.LE.2.D0*MT)THEN
	  HTT= 0.D0
	 ELSE
	  HT1= 3.D0*AFF(MA,(MT/MA)**2)*CU(4)**2
     .	     * TQCDA((MT/MA)**2)
	  IF(MA.LE.2.D0*RMT)THEN
	   HT2= 0.D0
	  ELSE
	   HT2= 3.D0*AFF(MA,(RMT/MA)**2)*CU(4)**2
     .	      * QCDA((RMT/MA)**2)
	  ENDIF
	  IF(HT2.LT.0.D0) HT2=0.D0
	  RAT= 2.D0*MT/MA
	  HTT= QQINT(RAT,HT1,HT2)
	 ENDIF
 
*   a -> gamma gamma

	 XFAC= CG(4)**2
	 HGG= GF/(32.D0*PI*SQR2)*MA**3*(ALEM0/PI)**2*XFAC

*    Below CG becomes the REDUCDED coupling
*    to 2 photons as compared to the SM coupling

	 CG(4)= CG(4)/CDABS(4.D0/3.D0*(CXTSM+CXCSM)+(CXBSM/3.D0+CXLSM))
	 
*   a -> Z1 gamma

	IF(MA.LE.MZ)THEN
	 HZG= 0.D0
	ELSE
	 FT= -2.D0*(CAZZ*(1.D0-8.D0/3.D0*S2TW)
     .  +2.D0*SAZZ*(DSQRT(S2TW)*NCP*(QQ-QUP)))/DSQRT(S2TW*C2TW)*CU(4)
	 FB= (CAZZ*(-1.D0+4.D0/3.D0*S2TW)
     .  +2.D0*SAZZ*(DSQRT(S2TW)*NCP*(QQ-QDOW)))/DSQRT(S2TW*C2TW)*CB(4)
	 FCH1=4D0*MW/(MCHA(1)*DSQRT(G2*S2TW*C2TW))*GACHACHA(1,1)
     .  *CAZZ*(-V(1,1)**2-0.5D0*V(1,2)**2
     .    -U(1,1)**2-0.5D0*U(1,2)**2+2.D0*S2TW)
     .  +SAZZ*DSQRT(S2TW)*NCP*(QD*U(1,2)**2-QU*V(1,2)**2)
	 FCH2=4D0*MW/(MCHA(2)*DSQRT(G2*S2TW*C2TW))*GACHACHA(2,2)
     .  *CAZZ*(-V(2,1)**2-0.5D0*V(2,2)**2
     .    -U(2,1)**2-0.5D0*U(2,2)**2+2.D0*S2TW)
     .  +SAZZ*DSQRT(S2TW)*NCP*(QD*U(2,2)**2-QU*V(2,2)**2)
	 CLT= 4.D0*(MT/MZ)**2*DCMPLX(1.D0,-EPS)
	 CLB= 4.D0*(MBP/MZ)**2*DCMPLX(1.D0,-EPS)
	 CLC= 4.D0*(MC/MZ)**2*DCMPLX(1.D0,-EPS)
	 CLCH1= 4d0*(MCHA(1)/MZ)**2*DCMPLX(1d0,-EPS)
	 CLCH2= 4d0*(MCHA(2)/MZ)**2*DCMPLX(1d0,-EPS)
	 CXTZ= FT*(-CI2(CTT,CLT))
	 CXBZ= FB*(-CI2(CTB,CLB))
	 CXCZ= FT*(-CI2(CTC,CLC))
	 CXCH1Z=FCH1*(-CI2(CTCH1,CLCH1))
	 CXCH2Z=FCH2*(-CI2(CTCH2,CLCH2))
	 XFAC= CDABS(CXTZ+CXBZ+CXCZ+CXCH1Z+CXCH2Z)**2
	 ACOUP= SQR2*GF*MZ**2*S2TW*C2TW/PI**2
	 HZG= GF/(4.D0*PI*SQR2)*MA**3*(ALEM0/PI)*ACOUP/16.D0
     .	    * XFAC*(1.D0-(MZ/MA)**2)**3
	ENDIF

*   a -> hZ1

	 DO J=1,3
	  IF(MA.LE.SMASS(J)+MZ)THEN
	   AHZ(J)= 0.D0
	  ELSE
	   RH= UPARF(157+J)
	   CH= LAMB((SMASS(J)/MA)**2,(MZ/MA)**2)
     C	     * LAMB((MA/MZ)**2,(SMASS(J)/MZ)**2)**2
	   AHZ(J)= GF/8.D0/SQR2/PI*MZ**4/MA*CH*RH**2
	   HTOT= HTOT+AHZ(J)
	  ENDIF
	 ENDDO

*   a -> h+W-

	 IF(MA.LE.CMASS+MW)THEN
	  HHCW= 0.D0
	 ELSE
	  RH= SINBETA*P(3,1)+COSBETA*P(3,2)
	  CH= LAMB((CMASS/MA)**2,(MW/MA)**2)
     C	    * LAMB((MA/MW)**2,(CMASS/MW)**2)**2
	  HHCW= GF/4.D0/SQR2/PI*MW**4/MA*CH*RH**2
	  HTOT= HTOT+HHCW
	 ENDIF

*   a -> neutralinos/charginos

	 STOT= 0.D0

	 LOOP= 0
	 DO J=1,6
	 DO K=J,6
	  IF(MA.LE.DABS(MNEU(J))+DABS(MNEU(K)))THEN
	   HNEU(J,K)= 0.D0
	  ELSE
	   HNEU(J,K)= 1.D0/(16.D0*PI)*(MA**2-(MNEU(J)-MNEU(K))**2)/MA
     C	        * LAMB((MNEU(J)/MA)**2,(MNEU(K)/MA)**2)
     C	        * GANEUNEU(J,K)**2
	   IF(J.NE.K)HNEU(J,K)= HNEU(J,K)/2.D0
	   STOT= STOT+HNEU(J,K)
	  ENDIF
	  HNEU(K,J)= HNEU(J,K)
	 ENDDO
	 LOOP= LOOP+(7-J) 
	 ENDDO

	 IF(MA.LE.2.D0*DABS(MCHA(1)))THEN
	  HCHA(1)= 0.D0
	 ELSE
	  HCHA(1)= 1.D0/(8.D0*PI)*MA
     C	       * LAMB((MCHA(1)/MA)**2,(MCHA(1)/MA)**2)
     C	       * GACHACHA(1,1)**2
	  STOT= STOT+HCHA(1)
	 ENDIF

	 IF(MA.LE.DABS(MCHA(1))+DABS(MCHA(2)))THEN
	  HCHA(2)= 0.D0
	 ELSE
	  HCHA(2)= 1.D0/(16.D0*PI)*MA
     C	       * LAMB((MCHA(1)/MA)**2,(MCHA(2)/MA)**2)
     C	       * ((GACHACHA(1,2)**2+GACHACHA(2,1)**2)
     C	       * (1.D0-(MCHA(1)/MA)**2-(MCHA(2)/MA)**2)
     C	       - 4.D0*GACHACHA(1,2)*GACHACHA(2,1)
     C	       * MCHA(1)*MCHA(2)/MA**2)
	  STOT= STOT+2.D0*HCHA(2)
	 ENDIF

	 IF(MA.LE.2.D0*DABS(MCHA(2)))THEN
	  HCHA(3)= 0.D0
	 ELSE
	  HCHA(3)= 1.D0/(8.D0*PI)*MA
     C	       * LAMB((MCHA(2)/MA)**2,(MCHA(2)/MA)**2)
     C	       * GACHACHA(2,2)**2
	  STOT= STOT+HCHA(3)
	 ENDIF

*   a -> squarks


	 IF(MA.LE.MST1+MST2)THEN
	  ASQ(1)= 0.D0
	 ELSE
	  RH= HIT1T2(3)
	  CH= LAMB((MST1/MA)**2,(MST2/MA)**2)
	  ASQ(1)= 3.D0*CH*RH**2/(16.D0*PI*MA)
	  STOT= STOT+2.D0*ASQ(1)
	 ENDIF

	 IF(MA.LE.MSB1+MSB2)THEN
	  ASQ(2)= 0.D0
	 ELSE
	  RH= HIB1B2(3)
	  CH= LAMB((MSB1/MA)**2,(MSB2/MA)**2)
	  ASQ(2)= 3.D0*CH*RH**2/(16.D0*PI*MA)
	  STOT= STOT+2.D0*ASQ(2)
	 ENDIF

*   a -> sleptons

	 IF(MA.LE.MSL1+MSL2)THEN
	  ASL= 0.D0
	 ELSE
	  RH= HIL1L2(3)
	  CH= LAMB((MSL1/MA)**2,(MSL2/MA)**2)
	  ASL= CH*RH**2/(16.D0*PI*MA)
	  STOT= STOT+2.D0*ASL
	 ENDIF

*  Branching ratios

	 WIDTH(4)= HJJ+HEE+HMM+HLL+HSS+HCC+HBB+HTT+HGG+HZG+HTOT+STOT
	 BRJJ(4)= HJJ/WIDTH(4)
	 BREE(4)= HEE/WIDTH(4)
	 BRMM(4)= HMM/WIDTH(4)
	 BRLL(4)= HLL/WIDTH(4)
	 BRSS(4)= HSS/WIDTH(4)
	 BRCC(4)= HCC/WIDTH(4)
	 BRBB(4)= HBB/WIDTH(4)
	 BRTT(4)= HTT/WIDTH(4)
	 BRGG(4)= HGG/WIDTH(4)
	 BRZG(4)= HZG/WIDTH(4)
	 BRHIGGS(4)= HTOT/WIDTH(4)
	 BRSUSY(4)= STOT/WIDTH(4)
	 DO J=1,3
	  BRAHZ(J)= AHZ(J)/WIDTH(4)
	 ENDDO
	 BRHCW(4)= HHCW/WIDTH(4)
	 DO J=1,6
	 DO K=1,6
	  BRNEU(4,J,K)=HNEU(J,K)/WIDTH(4)
	 ENDDO
	 ENDDO
	 DO J=1,3
	  BRCHA(4,J)=HCHA(J)/WIDTH(4)
	 ENDDO
	 DO J=1,2
	  BRASQ(J)=ASQ(J)/WIDTH(4)
	 ENDDO
	 BRASL=ASL/WIDTH(4)

* Charged Higgs boson

	ASH= ALPHAS(CMASS,2)

*  Running quark masses

	RMS= RUNM(CMASS,3)
	RMC= RUNM(CMASS,4)

*   Running bottom/top masses at CMASS:

	 IF(CMASS.GE.MT)THEN
	  RMT= RMTTOP
     C	   *(1.D0+7.D0/(4.D0*PI)*ASMT*DLOG(CMASS**2/MT**2))
     C	   **(-4.D0/7.D0)
	 ELSE
 	  RMT= RMTTOP
     C	   *(1.D0+23.D0/(12.D0*PI)*ASMT*DLOG(CMASS**2/MT**2))
     C	   **(-12.D0/23.D0)
     
	 ENDIF
     
         RMB=RUNMB(CMASS)

*  Partial widths

*   h+ -> enu

	IF(CMASS.LE.0.000511)THEN
	 HEN= 0.D0
	ELSE
	 HEN= CFF(CMASS,TANBETA,(0.000511/CMASS)**2,0.D0)
	ENDIF

*   h+ -> munu

	IF(CMASS.LE.MMUON)THEN
	 HMN= 0.D0
	ELSE
	 HMN= CFF(CMASS,TANBETA,(MMUON/CMASS)**2,0.D0)
	ENDIF

*   h+ -> taunu

	IF(CMASS.LE.MTAU)THEN
	 HLN= 0.D0
	ELSE
	 HLN= CFF(CMASS,TANBETA,(MTAU/CMASS)**2,0.D0)
	ENDIF

*   h+ -> su

	IF(CMASS.LE.MS+EPS)THEN
	 HSU= 0.D0
	ELSE
	 HSU1= 3.D0*VUS**2*CQCDM(CMASS,TANBETA,(MS/CMASS)**2,EPS)
	 HSU2= 3.D0*VUS**2*CQCD(CMASS,TANBETA,(RMS/CMASS)**2,EPS)
	 IF(HSU2.LT.0.D0) HSU2= 0.D0
	 RAT= MS/CMASS
	 HSU= QQINT(RAT,HSU1,HSU2)
	ENDIF

*   h+ -> cs

	IF(CMASS.LE.MS+MC)THEN
	 HSC= 0.D0
	ELSE
	 HSC1= 3.D0*CQCDM(CMASS,TANBETA,(MS/CMASS)**2,(MC/CMASS)**2)
	 HSC2= 3.D0*CQCD(CMASS,TANBETA,(RMS/CMASS)**2,(RMC/CMASS)**2)
	 IF(HSC2.LT.0.D0) HSC2= 0.D0
	 RAT= (MS+MC)/CMASS
	 HSC= QQINT(RAT,HSC1,HSC2)
	ENDIF

*   h+ -> cb

	IF(CMASS.LE.MBP+MC)THEN
	 HBC= 0.D0
	ELSE
	 HBC1= 3.D0*VCB**2
     C	  *CQCDM(CMASS,TANBETA,(MBP/CMASS)**2,(MC/CMASS)**2)
	 HBC2= 3.D0*VCB**2
     C	  *CQCD(CMASS,TANBETA,(RMB/CMASS)**2,(RMC/CMASS)**2)
	 IF(HBC2.LT.0.D0) HBC2= 0.D0
	 RAT= (MBP+MC)/CMASS
	 HBC= QQINT(RAT,HBC1,HBC2)
	ENDIF

*   Finite large tan(beta) SUSY correction:
	
	HBC=HBC/(1.D0+DELMB)**2


*   h+ -> bu

	IF(CMASS.LE.MBP+EPS)THEN
	 HBU= 0.D0
	ELSE
	 HBU1= 3.D0*VUB**2*CQCDM(CMASS,TANBETA,(MBP/CMASS)**2,EPS)
	 HBU2= 3.D0*VUB**2*CQCD(CMASS,TANBETA,(RMB/CMASS)**2,EPS)
	 IF(HBU2.LT.0.D0) HBU2= 0.D0
	 RAT= MBP/CMASS
	 HBU= QQINT(RAT,HBU1,HBU2)
	ENDIF

*   Finite large tan(beta) SUSY correction:

	HBU=HBU/(1.D0+DELMB)**2


*   h+ -> tb :

	IF(CMASS.LE.MT+MBP)THEN
	 HBT= 0.D0
	ELSE
	 HBT1= 3.D0*CQCDM(CMASS,TANBETA,(MBP/CMASS)**2,(MT/CMASS)**2)
	 IF(CMASS.LE.RMT+RMB)THEN
	  HBT2= 0.D0
	 ELSE
	  HBT2= 3.D0*CQCD(CMASS,TANBETA,(RMB/CMASS)**2,(RMT/CMASS)**2)
	 ENDIF
	 IF(HBT2.LT.0.D0) HBT2= 0.D0
	 RAT= (MBP+MT)/CMASS
	 HBT= QQINT(RAT,HBT1,HBT2)
	ENDIF

*   Finite large tan(beta) SUSY correction:

	HBT=HBT/(1.D0+DELMB)**2

*   h+ -> hW

	HTOT= 0.D0

	DO I=1,3
	 IF(CMASS.LE.MW+SMASS(I))THEN
	  HCWH(I)= 0.D0
	 ELSE
	  RH= S(I,1)*SINBETA-S(I,2)*COSBETA
	  CH= LAMB((SMASS(I)/CMASS)**2,(MW/CMASS)**2)
     C	   *LAMB((CMASS/MW)**2,(SMASS(I)/MW)**2)**2
	  HCWH(I)= GF/8.D0/SQR2/PI*MW**4/CMASS*CH*RH**2
	 HTOT= HTOT+HCWH(I)
	 ENDIF
	ENDDO

*   h+ -> aW

	 IF(CMASS.LE.MW+PMASS)THEN
	  HCWH(4)= 0.D0
	 ELSE
	  RH= SINBETA*P(3,1)+COSBETA*P(3,2)
	  CH= LAMB((PMASS/CMASS)**2,(MW/CMASS)**2)
     C	    * LAMB((CMASS/MW)**2,(PMASS/MW)**2)**2
	  HCWH(4)= GF/8.D0/SQR2/PI*MW**4/CMASS*CH*RH**2
	  HTOT= HTOT+HCWH(4)
	 ENDIF

*   h+ -> charginos+neutralinos

	STOT= 0.D0

	DO I=1,6
	DO J=1,2
	 IF (CMASS.LE.DABS(MNEU(I))+DABS(MCHA(J)))THEN
	  HCNC(I,J)= 0.D0
	 ELSE
	  HCNC(I,J)= 1.D0/(16.D0*PI)*CMASS
     C	       * LAMB((MNEU(I)/CMASS)**2,(MCHA(J)/CMASS)**2)
     C	       * ((GHCNEUCHAL(I,J)**2+GHCNEUCHAR(I,J)**2)
     C	       * (1.D0-(MNEU(I)/CMASS)**2-(MCHA(J)/CMASS)**2)
     C	       - 4.D0*GHCNEUCHAL(I,J)*GHCNEUCHAR(I,J)
     C	       * MNEU(I)*MCHA(J)/CMASS**2)

	  STOT= STOT+HCNC(I,J)
	 ENDIF
	ENDDO
	ENDDO

*   h+ -> squarks

*	 HPULDL=-g2*VEV*COSBETA*SINBETA/SQR2
	 HPULDL=-g2*H1Q*H2Q/DSQRT(H1Q**2+H2Q**2)
     

	 IF(CMASS.LE.MUL+MDL)THEN
	  HCSQ(1)= 0.D0
	 ELSE
	  RH= HPULDL
	  CH= LAMB((MUL/CMASS)**2,(MDL/CMASS)**2)
	  HCSQ(1)= 3.D0*CH*RH**2/(16.D0*PI*CMASS)
	  STOT= STOT+2.D0*HCSQ(1)
	 ENDIF

	 HPULDL=H1Q*((HDQ**2+HUQ**2)*SINBETA
     C	       -g2/DSQRT(1.D0+(H1Q/H2Q)**2))
	 HPULDR=HDQ*(MUEFFQ*COSBETA+AD*SINBETA)
	 HPURDL=HUQ*(MUEFFQ*SINBETA+AU*COSBETA)
	 HPURDR=HDQ*HUQ*DSQRT(H1Q**2+H2Q**2)
*	 HPT1B1=-UPARF(339)/(SINBETA*COSBETA*DSQRT(S2TW)*DSQRT(g2)*VEV**2)
*	 HPT1B2=-UPARF(340)/(SINBETA*COSBETA*DSQRT(S2TW)*DSQRT(g2)*VEV**2)
*	 HPT2B1=-UPARF(341)/(SINBETA*COSBETA*DSQRT(S2TW)*DSQRT(g2)*VEV**2)
*	 HPT2B2=-UPARF(342)/(SINBETA*COSBETA*DSQRT(S2TW)*DSQRT(g2)*VEV**2)

	 IF(CMASS.LE.MST1+MSB1)THEN
	  HCSQ(2)= 0.D0
	 ELSE
        RH= CST*CSB*HPULDL+SST*SSB*HPURDR
     .     +CST*SSB*HPULDR+SST*CSB*HPURDL
	  CH= LAMB((MST1/CMASS)**2,(MSB1/CMASS)**2)
	  HCSQ(2)= 3.D0*CH*RH**2/(16.D0*PI*CMASS)
	  STOT= STOT+HCSQ(2)
	 ENDIF

	 IF(CMASS.LE.MST1+MSB2)THEN
	  HCSQ(3)= 0.D0
	 ELSE
        RH= SGNB*(CST*SSB*HPULDL-SST*CSB*HPURDR
     .     -CST*CSB*HPULDR+SST*SSB*HPURDL)
	  CH= LAMB((MST1/CMASS)**2,(MSB2/CMASS)**2)
	  HCSQ(3)= 3.D0*CH*RH**2/(16.D0*PI*CMASS)
	  STOT= STOT+HCSQ(3)
	 ENDIF

	 IF(CMASS.LE.MST2+MSB1)THEN
	  HCSQ(4)= 0.D0
	 ELSE
        RH= SGNT*(SST*CSB*HPULDL-CST*SSB*HPURDR
     .     +SST*SSB*HPULDR-CST*CSB*HPURDL)
	  CH= LAMB((MST2/CMASS)**2,(MSB1/CMASS)**2)
	  HCSQ(4)= 3.D0*CH*RH**2/(16.D0*PI*CMASS)
	  STOT= STOT+HCSQ(4)
	 ENDIF

	 IF(CMASS.LE.MST2+MSB2)THEN
	  HCSQ(5)= 0.D0
	 ELSE
        RH= SGNT*SGNB*(SST*SSB*HPULDL+CST*CSB*HPURDR
     .     -SST*CSB*HPULDR-CST*SSB*HPURDL)
	  CH= LAMB((MST2/CMASS)**2,(MSB2/CMASS)**2)
	  HCSQ(5)= 3.D0*CH*RH**2/(16.D0*PI*CMASS)
	  STOT= STOT+HCSQ(5)
	 ENDIF

*   h+ -> sleptons

*	 HPLLNL=-g2*VEV*COSBETA*SINBETA/SQR2
	 HPLLNL=-g2*H1Q*H2Q/DSQRT(H1Q**2+H2Q**2)

	 IF(CMASS.LE.MLL+MNL)THEN
	  HCSL(1)= 0.D0
	 ELSE
	  RH= HPLLNL
	  CH= LAMB((MLL/CMASS)**2,(MNL/CMASS)**2)
	  HCSL(1)= CH*RH**2/(16.D0*PI*CMASS)
	  STOT= STOT+2.D0*HCSL(1)
	 ENDIF

*	 HPL1NL=-UPARF(343)*DSQRT(g2)/(4.D0*COSBETA*DSQRT(S2TW)*MW**2)
*	 HPL2NL=-UPARF(344)*DSQRT(g2)/(4.D0*COSBETA*DSQRT(S2TW)*MW**2)
	 HPLLNL=H1Q*H2Q/DSQRT(H1Q**2+H2Q**2)*(HLQ**2-g2)
	 HPLRNL=HLQ*(COSBETA*LAMBDA*ss+SINBETA*ATAU)     

	 IF(CMASS.LE.MSL1+MSNT)THEN
	  HCSL(2)= 0.D0
	 ELSE
	  RH= CSL*HPLLNL+SSL*HPLRNL
	  CH= LAMB((MSL1/CMASS)**2,(MSNT/CMASS)**2)
	  HCSL(2)= CH*RH**2/(16.D0*PI*CMASS)
	  STOT= STOT+HCSL(2)
	 ENDIF

	 IF(CMASS.LE.MSL2+MSNT)THEN
	  HCSL(3)= 0.D0
	 ELSE
	  RH= SGNL*(SSL*HPLLNL-CSL*HPLRNL)
	  CH= LAMB((MSL2/CMASS)**2,(MSNT/CMASS)**2)
	  HCSL(3)= CH*RH**2/(16.D0*PI*CMASS)
	  STOT= STOT+HCSL(3)
	 ENDIF

*  Branching ratios

	HCWIDTH= HLN+HMN+HEN+HSU+HBU+HSC+HBC+HBT+HTOT+STOT
      
	HCBRE= HEN/HCWIDTH
	HCBRM= HMN/HCWIDTH
	HCBRL= HLN/HCWIDTH
	HCBRSU= HSU/HCWIDTH
	HCBRBU= HBU/HCWIDTH
	HCBRSC= HSC/HCWIDTH
	HCBRBC= HBC/HCWIDTH
	HCBRBT= HBT/HCWIDTH
	HCBRWHT= HTOT/HCWIDTH
	HCBRSUSY= STOT/HCWIDTH
	DO I=1,4
	 HCBRWH(I)= HCWH(I)/HCWIDTH
	ENDDO
	DO I=1,6
	 DO J=1,2
	  HCBRNC(I,J)= HCNC(I,J)/HCWIDTH
	 ENDDO
	ENDDO
	DO I=1,5
	 HCBRSQ(I)= HCSQ(I)/HCWIDTH
	ENDDO
	DO I=1,3
	 HCBRSL(I)= HCSL(I)/HCWIDTH
	ENDDO

	RETURN

	END
    
    

	SUBROUTINE HTOVV(MV,GAMV,MH,HTVV)
	INTEGER I,IMAX
	PARAMETER (IMAX=20)
	DOUBLE PRECISION MV,GAMV,MH,HTVV,MV1,GAMV1,MH1
	DOUBLE PRECISION FTOVV1,DLT,SUM,UU,DD,RES
	COMMON/MV/MV1,GAMV1,MH1
	EXTERNAL FTOVV1
	MV1=MV
	GAMV1=GAMV
	MH1=MH
	DLT=1D0/DFLOAT(IMAX)
	SUM=0D0
	DO I=1,IMAX
	 UU=DLT*I
	 DD=UU-DLT
	 CALL QGAUS1(FTOVV1,DD,UU,RES)
	 SUM=SUM+RES
	ENDDO
	HTVV=SUM
 	END


	DOUBLE PRECISION FUNCTION FTOVV1(XX)
	INTEGER I,IMAX
	PARAMETER (IMAX=20)
	DOUBLE PRECISION DLT,SUM,UU,DD,RES
	DOUBLE PRECISION FTOVV,XX,X1
	COMMON/FIRST/X1
	EXTERNAL FTOVV2
	X1=XX
	DLT=1D0/DFLOAT(IMAX)
	SUM=0D0
	DO I=1,IMAX
	 UU=DLT*I
	 DD=UU-DLT
	 CALL QGAUS2(FTOVV2,DD,UU,RES)
	 SUM=SUM+RES
	ENDDO
	FTOVV1=SUM
	END


	DOUBLE PRECISION FUNCTION FTOVV2(XX)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION YY(2)
	COMMON/FIRST/X1
	YY(1)=X1
	YY(2)=XX
	FTOVV2=FTOVV(YY)
	END


	DOUBLE PRECISION FUNCTION FTOVV(XX)
	DOUBLE PRECISION LAMB,X,Y,XX,MV,GAMV,MH,PI
	DOUBLE PRECISION SP,SM,DJAC,PRO1,PRO2,Y1,Y2,Y3,T1,T2,AM2,GAM
	DIMENSION XX(2)
	COMMON/MV/MV,GAMV,MH
	LAMB(X,Y)=DSQRT((1.D0-X-Y)**2-4.D0*X*Y)
	PI=4D0*DATAN(1D0)
	IF(MH.LT.2*MV)THEN
	 SP = MH**2*XX(1)
	 SM = (MH-DSQRT(SP))**2*XX(2)
	 DJAC = MH**2*(MH-DSQRT(SP))**2
	 PRO1 = SP*GAMV/MV/((SP-MV**2)**2+MV**2*GAMV**2)
	 PRO2 = SM*GAMV/MV/((SM-MV**2)**2+MV**2*GAMV**2)
	ELSE
	 Y1 = DATAN((MH**2-MV**2)/MV/GAMV)
	 Y2 = -DATAN((MV**2)/MV/GAMV)
	 T1 = TAN(Y1*XX(1)+Y2*(1.D0-XX(1)))
	 SP = MV**2 + MV*GAMV*T1
	 Y3 = DATAN(((MH-DSQRT(SP))**2-MV**2)/MV/GAMV)
	 DJAC = (Y1-Y2)*(Y3-Y2)
	 T2 = TAN(Y3*XX(2)+Y2*(1.D0-XX(2)))
	 SM = MV**2 + MV*GAMV*T2
	 PRO1 = SP/MV**2
	 PRO2 = SM/MV**2
	ENDIF
	 AM2=MH**2
	 GAM = AM2*LAMB(SP/AM2,SM/AM2)*(1+LAMB(SP/AM2,SM/AM2)**2*MH**4
     C	  /SP/SM/12d0)
* In FTOVV: RADZZ=RADWW=1 (No EW rad. corrs. as in Hdecay)
	 FTOVV = PRO1*PRO2*GAM*DJAC/PI**2
	END


	SUBROUTINE QGAUS1(FUNC,A,B,SS)
c  Returns SS as integral of FUNC from A to B,
c  by 10-point Gauss-Legendre integration
	IMPLICIT REAL*8(A-Z)
	INTEGER J
	DIMENSION X(5),W(5)
	EXTERNAL FUNC
	DATA X/.1488743389D0,.4333953941D0,.6794095682D0
     C	  ,.8650633666D0,.9739065285D0/
	DATA W/.2955242247D0,.2692667193D0,.2190863625D0
     C	  ,.1494513491D0,.0666713443D0/
	XM=0.5D0*(B+A)
	XR=0.5D0*(B-A)
	SS=0.D0
	DO J=1,5
	  DX=XR*X(J)
	  SS=SS+W(J)*(FUNC(XM+DX)+FUNC(XM-DX))
	ENDDO
	SS=XR*SS
	END

	SUBROUTINE QGAUS2(FUNC,A,B,SS)
C     Returns SS as integral of FUNC from A to B, by 10-point Gauss-
C      Legendre integration
	IMPLICIT REAL*8(A-Z)
	INTEGER J
	DIMENSION X(5),W(5)
	EXTERNAL FUNC
	DATA X/.1488743389D0,.4333953941D0,.6794095682D0
     C	  ,.8650633666D0,.9739065285D0/
	DATA W/.2955242247D0,.2692667193D0,.2190863625D0
     C	  ,.1494513491D0,.0666713443D0/
	XM=0.5D0*(B+A)
	XR=0.5D0*(B-A)
	SS=0.D0
	DO J=1,5
	  DX=XR*X(J)
	  SS=SS+W(J)*(FUNC(XM+DX)+FUNC(XM-DX))
	ENDDO
	SS=XR*SS
	END
