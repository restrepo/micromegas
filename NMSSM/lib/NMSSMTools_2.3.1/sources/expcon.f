	SUBROUTINE SUBEXP(PAR,PROB)

***********************************************************************
*   Subroutine to check experimental constraints
*
*   The required data files and numbers (inv. Z width, lower bounds on
*   sparticle masses) are transferred via COMMON/LEP/...
*
*   On output:
*		PROB(I)  = 0, I = 1 - 26, 38: OK
*		
*		PROB(1) =/= 0	chargino too light
*		PROB(2) =/= 0	excluded by Z -> neutralinos
*		PROB(3) =/= 0	charged Higgs too light
*		PROB(4) =/= 0	excluded by ee -> hZ
*		PROB(5) =/= 0	excluded by ee -> hZ, h -> bb
*		PROB(6) =/= 0	excluded by ee -> hZ, h -> tautau
*		PROB(7) =/= 0	excluded by ee -> hZ, h -> invisible
*		PROB(8) =/= 0	excluded by ee -> hZ, h -> 2jets
*		PROB(9) =/= 0	excluded by ee -> hZ, h -> 2photons
*		PROB(10) =/= 0	excluded by ee -> hZ, h -> AA -> 4bs
*		PROB(11) =/= 0	excluded by ee -> hZ, h -> AA -> 4taus
*		PROB(12) =/= 0	excluded by ee -> hZ, h -> AA -> 2bs 2taus
*		PROB(13) =/= 0	excluded by Z -> hA (Z width)
*		PROB(14) =/= 0	excluded by ee -> hA -> 4bs
*		PROB(19) =/= 0	excluded by ee -> hA -> 4taus
*		PROB(16) =/= 0	excluded by ee -> hA -> 2bs 2taus
*		PROB(17) =/= 0	excluded by ee -> hA -> AAA -> 6bs
*		PROB(18) =/= 0	excluded by ee -> hA -> AAA -> 6taus
*		PROB(19) =/= 0	excluded by ee -> Zh -> ZAA -> Z+light pair
*		PROB(20) =/= 0	excluded by stop -> b l sneutrino
*		PROB(21) =/= 0	excluded by stop -> neutralino c
*		PROB(22) =/= 0	excluded by sbottom -> neutralino b
*		PROB(23) =/= 0	squark/gluino too light
*		PROB(24) =/= 0	selectron/smuon too light
*		PROB(25) =/= 0	stau too light
*		PROB(26) =/= 0	lightest neutralino is not LSP
*		PROB(38) =/= 0	excluded by Upsilon(1S) -> A gamma
*		PROB(39) =/= 0	excluded by etab(1S) mass difference
*
***********************************************************************

	IMPLICIT NONE

	INTEGER I,J,K
      	INTEGER NhZind,NhZbb,NhZll,NhZinv,NhZjj,NhZgg
	INTEGER NhA4b,NhA4tau,NhA2b2tau,NhA2tau2b
	INTEGER NAAA6b,NAAA6tau,NAAZ4b,NAAZ4tau,NAAZ2b2tau
	INTEGER Ncccc02,Ncccc04,Ncccc05,Ncccc06,Ncccc08,Ncccc1
	INTEGER Nccgg02,Nccgg04,Nccgg05,Nccgg06,Nccgg08,Nccgg1
	INTEGER Ncctt02,Ncctt04,Ncctt05,Ncctt06,Ncctt08,Ncctt1
	INTEGER Ngggg02,Ngggg04,Ngggg05,Ngggg06,Ngggg08,Ngggg1
	INTEGER Nttgg02,Nttgg04,Nttgg05,Nttgg06,Nttgg08,Nttgg1
	INTEGER Ntttt02,Ntttt04,Ntttt05,Ntttt06,Ntttt08,Ntttt1
	INTEGER Nstblsn,Nstnc,Nsbnb,Nglsq

	DOUBLE PRECISION PAR(*),PROB(*),PI,S,SQRS,SQR2,CONV
	DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
	DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2),CMASS
	DOUBLE PRECISION BRJJ(5),BRMM(5),BRLL(5),BRSS(5),BRCC(5)
	DOUBLE PRECISION BRBB(5),BRTT(5),BRWW(3),BRZZ(3),BRGG(5)
	DOUBLE PRECISION BRZG(5),BRHHH(4),BRHAA(3,3),BRHCHC(3)
	DOUBLE PRECISION BRHAZ(3,2),BRAHA(3),BRAHZ(2,3),BRHCW(5)
	DOUBLE PRECISION BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
	DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,6),BRASL(2,3)
	DOUBLE PRECISION BRSUSY(5),WIDTH(5),MBQM
	DOUBLE PRECISION MH,MA,TANB,h1,h2,R,Rmax,BRINV,M1,M2
	DOUBLE PRECISION GZ,GZMAX,ALSMZ,ALEMMZ,GF,g1,g2,S2TW,g
	DOUBLE PRECISION GZINV,GZINVMAX,MCHAMIN,MCMIN,SIGNEU1,SIGNEU
	DOUBLE PRECISION MSLMIN,MSTMIN,MSQMIN,MGLMIN,MMIN
	DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
	DOUBLE PRECISION LAMBDA,X,Y,Z,Q,O,DZ,E1,E2,SIG
	DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
	DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
	DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU,MAtest,ceff
	DOUBLE PRECISION hZind(1000,2),hZbb(1000,2),hZll(1000,2)
	DOUBLE PRECISION hZinv(1000,2),hZjj(1000,2),hZgg(1000,2)
	DOUBLE PRECISION hA4b(10000,3),hA4tau(10000,3)
	DOUBLE PRECISION hA2b2tau(10000,3),hA2tau2b(10000,3)
	DOUBLE PRECISION AAA6b(10000,3),AAA6tau(10000,3)
	DOUBLE PRECISION AAZ4b(10000,3),AAZ4tau(10000,3)
	DOUBLE PRECISION AAZ2b2tau(10000,3)
	DOUBLE PRECISION cccc02(100,2),cccc04(100,2),cccc05(100,2)
	DOUBLE PRECISION cccc06(100,2),cccc08(100,2),cccc1(100,2)
	DOUBLE PRECISION ccgg02(100,2),ccgg04(100,2),ccgg05(100,2)
	DOUBLE PRECISION ccgg06(100,2),ccgg08(100,2),ccgg1(100,2)
	DOUBLE PRECISION cctt02(100,2),cctt04(100,2),cctt05(100,2)
	DOUBLE PRECISION cctt06(100,2),cctt08(100,2),cctt1(100,2)
	DOUBLE PRECISION gggg02(100,2),gggg04(100,2),gggg05(100,2)
	DOUBLE PRECISION gggg06(100,2),gggg08(100,2),gggg1(100,2)
	DOUBLE PRECISION ttgg02(100,2),ttgg04(100,2),ttgg05(100,2)
	DOUBLE PRECISION ttgg06(100,2),ttgg08(100,2),ttgg1(100,2)
	DOUBLE PRECISION tttt02(100,2),tttt04(100,2),tttt05(100,2)
	DOUBLE PRECISION tttt06(100,2),tttt08(100,2),tttt1(100,2)
	DOUBLE PRECISION stblsn(100,2),stnc(100,2),sbnb(100,2)
	DOUBLE PRECISION glsq(100,2)
	DOUBLE PRECISION ALEM0,MY,BRYMUMU,C,ZZ,AP,GYEE
	DOUBLE PRECISION ALPHAS,ALSMY,DELTA,CLEOTAU,CLEOMU
	DOUBLE PRECISION M0,RETA,GAM2,XX,YY,F,YMAX,FMAX,D,MEMAX,UU,VV
	DOUBLE PRECISION CU(5),CD(5),CV(3),CJ(5),CG(5)

	COMMON/ALEM0/ALEM0
	COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
	COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
	COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
	COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
	COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     C		MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     C		CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
	COMMON/BRN/BRJJ,BRMM,BRLL,BRSS,BRCC,BRBB,BRTT,BRWW,BRZZ,
     C		BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     C		BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     C		BRSUSY,WIDTH
	COMMON/LEP/GZMAX,GZINVMAX,MCHAMIN,MCMIN,SIGNEU1,SIGNEU,
     C		MSLMIN,MSTMIN,MSQMIN,MGLMIN,
     C		hZind,hZbb,hZll,hZinv,hZjj,hZgg,
     C		hA4b,hA4tau,hA2b2tau,hA2tau2b,
     C		AAA6b,AAA6tau,AAZ4b,AAZ4tau,AAZ2b2tau,
     C          cccc02,cccc04,cccc05,cccc06,cccc08,cccc1,
     C          ccgg02,ccgg04,ccgg05,ccgg06,ccgg08,ccgg1,
     C          cctt02,cctt04,cctt05,cctt06,cctt08,cctt1,
     C          gggg02,gggg04,gggg05,gggg06,gggg08,gggg1,
     C          ttgg02,ttgg04,ttgg05,ttgg06,ttgg08,ttgg1,
     C          tttt02,tttt04,tttt05,tttt06,tttt08,tttt1,
     C		stblsn,stnc,sbnb,glsq,
     C		NhZind,NhZbb,NhZll,NhZinv,NhZjj,NhZgg,
     C		NhA4b,NhA4tau,NhA2b2tau,NhA2tau2b,
     C		NAAA6b,NAAA6tau,NAAZ4b,NAAZ4tau,NAAZ2b2tau,
     C          Ncccc02,Ncccc04,Ncccc05,Ncccc06,Ncccc08,Ncccc1,
     C          Nccgg02,Nccgg04,Nccgg05,Nccgg06,Nccgg08,Nccgg1,
     C          Ncctt02,Ncctt04,Ncctt05,Ncctt06,Ncctt08,Ncctt1,
     C          Ngggg02,Ngggg04,Ngggg05,Ngggg06,Ngggg08,Ngggg1,
     C          Nttgg02,Nttgg04,Nttgg05,Nttgg06,Nttgg08,Nttgg1,
     C          Ntttt02,Ntttt04,Ntttt05,Ntttt06,Ntttt08,Ntttt1,
     C		Nstblsn,Nstnc,Nsbnb,Nglsq
	COMMON/REDCOUP/CU,CD,CV,CJ,CG

	LAMBDA(X,Y,Z)= DSQRT(X**2+Y**2+Z**2-2.D0*X*Y-2.D0*X*Z-2.D0*Y*Z)

	PI=4.D0*DATAN(1D0)
	SQR2=DSQRT(2.D0)
	SQRS=209.D0
	S=SQRS**2
	CONV=.3894D9
	g=(g1+g2)/2.D0
	TANB=PAR(3)
	h2=1.D0/DSQRT(2.D0*DSQRT(2.D0)*(1.D0+TANB**2)*GF)
	h1=h2*TANB

* Test on LSP

	PROB(26)=DDIM(1.D0,MIN(DABS(MGL),DABS(MCH(1)),MUR,MUL,
     C	 MDR,MDL,MLR,MLL,MNL,MST1,MSB1,MSL1,MSNT)/DABS(MNEU(1)))

* Test on stop -> b l sneutrino

	I=1
	DOWHILE(stblsn(I,1).LE.MST1 .AND. I.LT.Nstblsn)
	 I=I+1
	ENDDO
	IF(I.GT.1 .AND. MST1.LT.stblsn(Nstblsn,1))THEN
	 MMIN=stblsn(I-1,2)+(MST1-stblsn(I-1,1))
     C	  /(stblsn(I,1)-stblsn(I-1,1))*(stblsn(I,2)-stblsn(I-1,2))
	 PROB(20)=DDIM(1.D0,MNL/MMIN)
	ENDIF

* Test on stop -> neutralino c

	I=1
	DOWHILE(stnc(I,1).LE.DABS(MNEU(1)) .AND. I.LT.Nstnc)
	 I=I+1
	ENDDO
	IF(I.GT.1 .AND. DABS(MNEU(1)).LT.stnc(Nstnc,1))THEN
	 MMIN=stnc(I-1,2)+(DABS(MNEU(1))-stnc(I-1,1))
     C	  /(stnc(I,1)-stnc(I-1,1))*(stnc(I,2)-stnc(I-1,2))
	 PROB(21)=DDIM(1.D0,MST1/MMIN)
	ENDIF

* Test on sbottom -> neutralino b

	I=1
	DOWHILE(sbnb(I,1).LE.MSB1 .AND. I.LT.Nsbnb)
	 I=I+1
	ENDDO
	IF(I.GT.1 .AND. MSB1.LT.sbnb(Nsbnb,1))THEN
	 MMIN=sbnb(I-1,2)+(MSB1-sbnb(I-1,1))/(sbnb(I,1)-sbnb(I-1,1))
     C	  *(sbnb(I,2)-sbnb(I-1,2))
	 PROB(22)=DDIM(1.D0,DABS(MNEU(1))/MMIN)
	ENDIF

* Test on gluino/squark masses

	I=1
	DOWHILE(glsq(I,1).LE.DABS(MGL) .AND. I.LT.Nglsq)
	 I=I+1
	ENDDO
	MMIN=MSQMIN
	IF(I.GT.1 .AND. DABS(MGL).LT.glsq(Nglsq,1))THEN
	 MMIN=glsq(I-1,2)+(DABS(MGL)-glsq(I-1,1))
     C	  /(glsq(I,1)-glsq(I-1,1))*(glsq(I,2)-glsq(I-1,2))
	ENDIF
	PROB(23)=DDIM(1.D0,MIN(MUR,MUL,MDR,MDL)/MMIN)
     C	 +DDIM(1.D0,DABS(MGL)/MGLMIN)

* Test on slepton masses

	PROB(24)=DDIM(1.D0,MIN(MLR,MLL)/MSLMIN)
	PROB(25)=DDIM(1.D0,MSL1/MSTMIN)

* Test on chargino mass

	PROB(1)=DDIM(1.D0,DABS(MCH(1))/MCHAMIN)

* Test on Z width into neutralinos

	IF(DABS(MNEU(1)).LT.MZ/2.D0)THEN
	 GZINV=MZ**3*GF/(12.D0*DSQRT(2.D0)*PI)
     C     *(NEU(1,3)**2-NEU(1,4)**2)**2*(1.D0-4*MNEU(1)**2/MZ**2)**1.5
 	 PROB(2)=DDIM(GZINV/GZINVMAX,1.D0)
	ENDIF

* Test on neutralinos

	DZ = 1.D0/(S-MZ**2)**2
	DO I=1,5
	 DO J=I,5
	  IF(DABS(MNEU(I))+DABS(MNEU(J)).LT.SQRS)THEN
	   O= (NEU(I,3)*NEU(J,3)-NEU(I,4)*NEU(J,4))
	   E1= (S+MNEU(I)**2-MNEU(J)**2)/(2.D0*SQRS)
	   E2= (S+MNEU(J)**2-MNEU(I)**2)/(2.D0*SQRS)
	   Q= LAMBDA(S,MNEU(I)**2,MNEU(J)**2)/(2.D0*SQRS)
	   SIG= 1.D0/(4.D0*PI)*Q/SQRS*DZ*(E1*E2+Q**2/3.D0
     C	       -MNEU(I)*MNEU(J))* g**2*O**2*(.25D0-S2TW+2.D0*S2TW**2)
	   SIG= CONV*SIG
	   IF(I.EQ.J)SIG= SIG/2.D0
	   IF(I.EQ.1)THEN
	    IF(J.GT.1)PROB(2)=PROB(2)+DDIM(SIG/SIGNEU1,1.D0)
	   ELSE
	    PROB(2)=PROB(2)+DDIM(SIG/SIGNEU,1.D0)
	   ENDIF
	  ENDIF
	 ENDDO
	ENDDO

* Test on charged Higgs mass

	 PROB(3)=DDIM(1.D0,CMASS/MCMIN)
	 
* Light A Physics

	MA=PMASS(1)

* Test on Upsilon(1S) -> A gamma (from CLEO)

	MY=9.46D0 ! Upsilon(1S) mass
	MBQM=4.9D0 ! b quark mass in quark models
	ALSMY=ALPHAS(MY,2) ! alpha_s at MY, 2 loop caclulation
	BRYMUMU=2.48D-2 ! BR(Upsilon(1S) -> mu mu)

	IF(MA.LT.MY)THEN

	ZZ=1.D0-MA**2/MY**2 ! energy fraction of the photon
	AP=6.D0*ZZ**.2D0 ! Nason function for QCD corrections
	C=1.D0+4.D0*ALSMY/(3.D0*PI)*(4.D0-AP) ! QCD corrections
	DELTA=1.2D0**2/MBQM**2 ! function for rel. corrections
	C=C* ! relativistic corrections (for MA<~8.8 GeV)
     C	  (MY**2-MA**2)**2/(4.D0*MBQM**2-MA**2)**2*(1.D0-
     C	  DELTA/3.D0*(36.D0*MBQM**2+MA**2)/(4.D0*MBQM**2-MA**2))

	C=MAX(C,1.D-6)

	RMAX=0.D0
	RMAX=SQR2*PI*ALEM0*CLEOTAU(MA)/(GF*MBQM**2*ZZ*C*BRYMUMU)
	IF(RMAX.NE.0.D0)PROB(38)=DDIM(CD(4)**2*BRLL(4)/RMAX,1.D0)

	RMAX=0.D0
	RMAX=SQR2*PI*ALEM0*CLEOMU(MA)/(GF*MBQM**2*ZZ*C*BRYMUMU)
	IF(RMAX.NE.0.D0)PROB(38)=DDIM(CD(4)**2*BRMM(4)/RMAX,1.D0)

	IF(WIDTH(4).GT.1.D-2)PROB(38)=-PROB(38)

	ENDIF

* Test on etab(1S) mass difference (BABAR - theory)

	M0=9.389D0 ! etab(1S) mass
	GYEE=1.34D-6 ! Gamma(Upsilon(1S) -> e+ e-)
	RETA=GYEE*9.D0*MY**2/(4.D0*ALEM0**2)*
     C	  (1.D0+16.D0*ALSMY/(3.D0*PI)) ! radial wave fun. at the origin
 	! Resolution of the 3rd degree eq. for the limit on Xd
	GAM2=(M0*2.D-2)**2
	XX=MA**2-M0**2

	IF(MA.LT.M0)THEN
	 MEMAX=M0-3.D-2
	ELSE
	 MEMAX=M0+4.D-2
	ENDIF
	YMAX=MEMAX**2-M0**2
	FMAX=XX*YMAX*(1.D0+GAM2/(XX+YMAX)**2)

	D=XX**2-GAM2/27.D0
	IF(D.LT.0.D0)THEN
	 UU=2.D0*DSQRT(GAM2)/DSQRT(3.D0)
	 VV=-3.D0*DSQRT(3.D0)*XX/DSQRT(GAM2)
	 YY=UU*DCOS(DACOS(VV)/3.D0 + 4.D0*PI/3.D0)
	 F=XX*YY*(1.D0+GAM2/(XX+YY)**2)
	 FMAX=MAX(FMAX,F)
	ENDIF

	RMAX=8.D0*PI*MZ**2/(3.D0*g*RETA*M0**3)*FMAX
	PROB(39)=DDIM(CD(4)**2/RMAX,1.D0)

* Higgs Strahlung

	DO I=1,3
	MH=SMASS(I)
	IF(MH+MZ.LT.SQRS)THEN

	 R=(SCOMP(I,1)*h1+SCOMP(I,2)*h2)**2/(h1**2+h2**2)

*  ee -> hZ flavor independent

	 Rmax=1.D0
	 J=1
	 DOWHILE(hZind(J,1).LE.MH .AND. J.LT.NhZind)
	  J=J+1
	 ENDDO
	 IF(J.GT.1 .AND. MH.LT.hZind(NhZind,1))
     C    Rmax=hZind(J-1,2)+(MH-hZind(J-1,1))/(hZind(J,1)
     C         -hZind(J-1,1))*(hZind(J,2)-hZind(J-1,2))
	 PROB(4)=PROB(4)+DDIM(R/Rmax,1.D0)

*  ee -> hZ, h -> bb

	 Rmax=1.D0
	 J=1
	 DOWHILE(hZbb(J,1).LE.MH .AND. J.LT.NhZbb)
	  J=J+1
	 ENDDO
	 IF(J.GT.1 .AND. MH.LT.hZbb(NhZbb,1))
     C    Rmax=hZbb(J-1,2)+(MH-hZbb(J-1,1))/(hZbb(J,1)-hZbb(J-1,1))
     C         *(hZbb(J,2)-hZbb(J-1,2))
	 PROB(5)=PROB(5)+DDIM(R*BRBB(I)/Rmax,1.D0)

*  ee -> hZ, h -> tautau

	 Rmax=1.D0
	 J=1
	 DOWHILE(hZll(J,1).LE.MH .AND. J.LT.NhZll)
	  J=J+1
	 ENDDO
	 IF(J.GT.1 .AND. MH.LT.hZll(NhZll,1))
     C    Rmax=hZll(J-1,2)+(MH-hZll(J-1,1))/(hZll(J,1)-hZll(J-1,1))
     C         *(hZll(J,2)-hZll(J-1,2))
	 PROB(6)=PROB(6)+DDIM(R*BRLL(I)/Rmax,1.D0)

*  ee -> hZ, h -> invisible

	 Rmax=1.D0
	 J=1
	 DOWHILE(hZinv(J,1).LE.MH .AND. J.LT.NhZinv)
	  J=J+1
	 ENDDO
	 IF(J.GT.1 .AND. MH.LT.hZinv(NhZinv,1))
     C    Rmax=hZinv(J-1,2)+(MH-hZinv(J-1,1))/(hZinv(J,1)-hZinv(J-1,1))
     C         *(hZinv(J,2)-hZinv(J-1,2))
	 BRINV=BRNEU(I,1,1)
     C	  +BRHAA(I,1)*BRNEU(4,1,1)**2
     C	  +BRHAA(I,2)*BRNEU(4,1,1)*BRNEU(5,1,1)
     C	  +BRHAA(I,3)*BRNEU(5,1,1)**2
	 IF(I.EQ.2)
     C	  BRINV=BRINV+BRHHH(1)*BRNEU(1,1,1)**2
	 IF(I.EQ.3)
     C	  BRINV=BRINV+BRHHH(2)*BRNEU(1,1,1)**2
     C	   +BRHHH(3)*BRNEU(1,1,1)*BRNEU(2,1,1)
     C	   +BRHHH(4)*BRNEU(2,1,1)**2
	 PROB(7)=PROB(7)+DDIM(R*BRINV/Rmax,1.D0)

*  ee -> hZ, h -> 2jets

	 Rmax=1.D0
	 J=1
	 DOWHILE(hZjj(J,1).LE.MH .AND. J.LT.NhZjj)
	  J=J+1
	 ENDDO
	 IF(J.GT.1 .AND. MH.LT.hZjj(NhZjj,1))
     C    Rmax=hZjj(J-1,2)+(MH-hZjj(J-1,1))/(hZjj(J,1)-hZjj(J-1,1))
     C         *(hZjj(J,2)-hZjj(J-1,2))
	 PROB(8)=PROB(8)
     C	         +DDIM(R*(BRJJ(I)+BRSS(I)+BRCC(I)+BRBB(I))/Rmax,1.D0)

*  ee -> hZ, h -> 2photons

	 Rmax=1.D0
	 J=1
	 DOWHILE(hZgg(J,1).LE.MH .AND. J.LT.NhZgg)
	  J=J+1
	 ENDDO
	 IF(J.GT.1 .AND. MH.LT.hZgg(NhZgg,1))
     C    Rmax=hZgg(J-1,2)+(MH-hZgg(J-1,1))/(hZgg(J,1)-hZgg(J-1,1))
     C         *(hZgg(J,2)-hZgg(J-1,2))
	 PROB(9)=PROB(9)+DDIM(R*BRGG(I)/Rmax,1.D0)

*  ee -> hZ, h -> AA -> 4bs

	 MA=PMASS(1)
	 IF(MH.GT.2.D0*MA)THEN
	  Rmax=1.D0
	  DO K=1,NAAZ4b
	   IF(AINT(MH).EQ.AAZ4b(K,1) .AND. AINT(MA).EQ.AAZ4b(K,2))THEN
	    Rmax=AAZ4b(K,3)
	    GOTO 1
	   ENDIF
	  ENDDO
 1	  PROB(10)=PROB(10)+DDIM(R*BRHAA(I,1)*BRBB(4)**2/Rmax,1.D0)
	 ENDIF

	 MA=SMASS(1)
	 IF(MH.GT.2.D0*MA)THEN
	  Rmax=1.D0
	  DO K=1,NAAZ4b
	   IF(AINT(MH).EQ.AAZ4b(K,1) .AND. AINT(MA).EQ.AAZ4b(K,2))THEN
	    Rmax=AAZ4b(K,3)
	    GOTO 2
	   ENDIF
	  ENDDO
 2	  PROB(10)=PROB(10)+DDIM(R*BRHHH(I-1)*BRBB(1)**2/Rmax,1.D0)
	 ENDIF

*  ee -> hZ, h -> AA -> 4taus

	 MA=PMASS(1)
	 IF(MH.GT.2.D0*MA)THEN
	  Rmax=1.D0
	  DO K=1,NAAZ4tau
	   IF(AINT(MH).EQ.AAZ4tau(K,1) .AND.
     C	     AINT(MA).EQ.AAZ4tau(K,2))THEN
	    Rmax=AAZ4tau(K,3)
	    GOTO 11
	   ENDIF
	  ENDDO
 11	  PROB(11)=PROB(11)+DDIM(R*BRHAA(I,1)*BRLL(4)**2/Rmax,1.D0)
	 ENDIF

	 MA=SMASS(1)
	 IF(MH.GT.2.D0*MA)THEN
	  Rmax=1.D0
	  DO K=1,NAAZ4tau
	   IF(AINT(MH).EQ.AAZ4tau(K,1) .AND.
     C	     AINT(MA).EQ.AAZ4tau(K,2))THEN
	    Rmax=AAZ4tau(K,3)
	    GOTO 12
	   ENDIF
	  ENDDO
 12	  PROB(11)=PROB(11)+DDIM(R*BRHHH(I-1)*BRLL(1)**2/Rmax,1.D0)
	 ENDIF

*  ee -> hZ, h -> AA -> 2bs 2taus

	 MA=PMASS(1)
	 IF(MH.GT.2.D0*MA)THEN
	  Rmax=1.D0
	  DO K=1,NAAZ2b2tau
	   IF(AINT(MH).EQ.AAZ2b2tau(K,1) .AND.
     C	     AINT(MA).EQ.AAZ2b2tau(K,2))THEN
	    Rmax=AAZ2b2tau(K,3)
	    GOTO 21
	   ENDIF
	  ENDDO
 21	  PROB(12)=PROB(12)
     C	   +DDIM(R*BRHAA(I,1)*2.D0*BRLL(4)*BRBB(4)/Rmax,1.D0)
	 ENDIF

	 MA=SMASS(1)
	 IF(MH.GT.2.D0*MA)THEN
	  Rmax=1.D0
	  DO K=1,NAAZ2b2tau
	   IF(AINT(MH).EQ.AAZ2b2tau(K,1) .AND.
     C	     AINT(MA).EQ.AAZ2b2tau(K,2))THEN
	    Rmax=AAZ2b2tau(K,3)
	    GOTO 22
	   ENDIF
	  ENDDO
 22	  PROB(12)=PROB(12)
     C	   +DDIM(R*BRHHH(I-1)*2.D0*BRBB(1)*BRLL(1)/Rmax,1.D0)
	 ENDIF

*  ee -> hZ -> AAZ -> Z + light pairs

	 MA=PMASS(1)
	 IF(MH.LT.2.D0*MA .OR. MH.LT.40.D0 .OR. MH.GT.90.D0
     C	   .OR. MA.LT.2.D0 .OR. MA.GT.12.D0)GOTO 752

*      AA -> cccc

	 ceff=R*BRHAA(I,1)*BRCC(4)**2
	 IF(ceff.GE.0.0.AND.ceff.LT.0.2)GOTO 102
	 IF(ceff.GE.0.2.AND.ceff.LT.0.4)GOTO 202
	 IF(ceff.GE.0.4.AND.ceff.LT.0.5)GOTO 302
	 IF(ceff.GE.0.5.AND.ceff.LT.0.6)GOTO 402
	 IF(ceff.GE.0.6.AND.ceff.LT.0.8)GOTO 502
	 IF(ceff.GE.0.8.AND.ceff.LT.1.0)GOTO 602
	 GOTO 702

 102	 Rmax=.2D0
	 J=1
	 DOWHILE(cccc02(J,1).LE.MH .AND. J.LT.Ncccc02)
	  J=J+1
	 ENDDO
	 MAtest=0.D0
	 IF(J.GT.1 .AND. MH.LT.cccc02(Ncccc02,1))
     C	  MAtest=cccc02(J-1,2)+(MH-cccc02(J-1,1))/
     C	   (cccc02(J,1)-cccc02(J-1,1))
     C	   *(cccc02(J,2)-cccc02(J-1,2))
	 PROB(19)=PROB(19)+5.D0*(Rmax-ceff)*DDIM(MAtest/MA,1.D0)
	 GOTO 702

 202	 Rmax=.4D0
	 J=1
	 DOWHILE(cccc04(J,1).LE.MH .AND. J.LT.Ncccc04)
	  J=J+1
	 ENDDO
	 MAtest=0.D0
	 IF(J.GT.1 .AND. MH.LT.cccc04(Ncccc04,1))
     C	  MAtest=cccc04(J-1,2)+(MH-cccc04(J-1,1))/
     C	   (cccc04(J,1)-cccc04(J-1,1))
     C	   *(cccc04(J,2)-cccc04(J-1,2))
	 PROB(19)=PROB(19)+5.D0*(Rmax-ceff)*DDIM(MAtest/MA,1.D0)
	 GOTO 702

 302	 Rmax=.5D0
	 J=1
	 DOWHILE(cccc05(J,1).LE.MH .AND. J.LT.Ncccc05)
	  J=J+1
	 ENDDO
	 MAtest=0.D0
	 IF(J.GT.1 .AND. MH.LT.cccc05(Ncccc05,1))
     C	  MAtest=cccc05(J-1,2)+(MH-cccc05(J-1,1))/
     C	   (cccc05(J,1)-cccc05(J-1,1))
     C	   *(cccc05(J,2)-cccc05(J-1,2))
	 PROB(19)=PROB(19)+10.D0*(Rmax-ceff)*DDIM(MAtest/MA,1.D0)
	 GOTO 702

 402	 Rmax=.6D0
	 J=1
	 DOWHILE(cccc06(J,1).LE.MH .AND. J.LT.Ncccc06)
	  J=J+1
	 ENDDO
	 MAtest=0.D0
	 IF(J.GT.1 .AND. MH.LT.cccc06(Ncccc06,1))
     C	  MAtest=cccc06(J-1,2)+(MH-cccc06(J-1,1))/
     C	   (cccc06(J,1)-cccc06(J-1,1))
     C	   *(cccc06(J,2)-cccc06(J-1,2))
	 PROB(19)=PROB(19)+10.D0*(Rmax-ceff)*DDIM(MAtest/MA,1.D0)
	 GOTO 702

 502	 Rmax=.8D0
	 J=1
	 DOWHILE(cccc08(J,1).LE.MH .AND. J.LT.Ncccc08)
	  J=J+1
	 ENDDO
	 MAtest=0.D0
	 IF(J.GT.1 .AND. MH.LT.cccc08(Ncccc08,1))
     C	  MAtest=cccc08(J-1,2)+(MH-cccc08(J-1,1))/
     C	   (cccc08(J,1)-cccc08(J-1,1))
     C	   *(cccc08(J,2)-cccc08(J-1,2))
	 PROB(19)=PROB(19)+5.D0*(Rmax-ceff)*DDIM(MAtest/MA,1.D0)
	 GOTO 702

 602	 Rmax=1.D0
	 J=1
	 DOWHILE(cccc1(J,1).LE.MH .AND. J.LT.Ncccc1)
	  J=J+1
	 ENDDO
	 MAtest=0.D0
	 IF(J.GT.1 .AND. MH.LT.cccc1(Ncccc1,1))
     C	  MAtest=cccc1(J-1,2)+(MH-cccc1(J-1,1))/
     C	   (cccc1(J,1)-cccc1(J-1,1))
     C	   *(cccc1(J,2)-cccc1(J-1,2))
	 PROB(19)=PROB(19)+5.D0*(Rmax-ceff)*DDIM(MAtest/MA,1.D0)

 702	 CONTINUE

*      AA -> ccjj

	 ceff=R*BRHAA(I,1)*BRCC(4)*(BRJJ(4)+BRSS(4))*2.D0
	 IF(ceff.GE.0.0.AND.ceff.LT.0.2)GOTO 112
	 IF(ceff.GE.0.2.AND.ceff.LT.0.4)GOTO 212
	 IF(ceff.GE.0.4.AND.ceff.LT.0.5)GOTO 312
	 IF(ceff.GE.0.5.AND.ceff.LT.0.6)GOTO 412
	 IF(ceff.GE.0.6.AND.ceff.LT.0.8)GOTO 512
	 IF(ceff.GE.0.8.AND.ceff.LT.1.0)GOTO 612
	 GOTO 712

 112	 Rmax=.2D0
	 J=1
	 DOWHILE(ccgg02(J,1).LE.MH .AND. J.LT.Nccgg02)
	  J=J+1
	 ENDDO
	 MAtest=0.D0
	 IF(J.GT.1 .AND. MH.LT.ccgg02(Nccgg02,1))
     C	  MAtest=ccgg02(J-1,2)+(MH-ccgg02(J-1,1))/
     C	   (ccgg02(J,1)-ccgg02(J-1,1))
     C	   *(ccgg02(J,2)-ccgg02(J-1,2))
	 PROB(19)=PROB(19)+5.D0*(Rmax-ceff)*DDIM(MAtest/MA,1.D0)
	 GOTO 712

 212	 Rmax=.4D0
	 J=1
	 DOWHILE(ccgg04(J,1).LE.MH .AND. J.LT.Nccgg04)
	  J=J+1
	 ENDDO
	 MAtest=0.D0
	 IF(J.GT.1 .AND. MH.LT.ccgg04(Nccgg04,1))
     C	  MAtest=ccgg04(J-1,2)+(MH-ccgg04(J-1,1))/
     C	   (ccgg04(J,1)-ccgg04(J-1,1))
     C	   *(ccgg04(J,2)-ccgg04(J-1,2))
	 PROB(19)=PROB(19)+5.D0*(Rmax-ceff)*DDIM(MAtest/MA,1.D0)
	 GOTO 712

 312	 Rmax=.5D0
	 J=1
	 DOWHILE(ccgg05(J,1).LE.MH .AND. J.LT.Nccgg05)
	  J=J+1
	 ENDDO
	 MAtest=0.D0
	 IF(J.GT.1 .AND. MH.LT.ccgg05(Nccgg05,1))
     C	  MAtest=ccgg05(J-1,2)+(MH-ccgg05(J-1,1))/
     C	   (ccgg05(J,1)-ccgg05(J-1,1))
     C	   *(ccgg05(J,2)-ccgg05(J-1,2))
	 PROB(19)=PROB(19)+10.D0*(Rmax-ceff)*DDIM(MAtest/MA,1.D0)
	 GOTO 712

 412	 Rmax=.6D0
	 J=1
	 DOWHILE(ccgg06(J,1).LE.MH .AND. J.LT.Nccgg06)
	  J=J+1
	 ENDDO
	 MAtest=0.D0
	 IF(J.GT.1 .AND. MH.LT.ccgg06(Nccgg06,1))
     C    MAtest=ccgg06(J-1,2)+(MH-ccgg06(J-1,1))/
     C	   (ccgg06(J,1)-ccgg06(J-1,1))
     C	   *(ccgg06(J,2)-ccgg06(J-1,2))
	 PROB(19)=PROB(19)+10.D0*(Rmax-ceff)*DDIM(MAtest/MA,1.D0)
	 GOTO 712

 512	 Rmax=.8D0
	 J=1
	 DOWHILE(ccgg08(J,1).LE.MH .AND. J.LT.Nccgg08)
	  J=J+1
	 ENDDO
	 MAtest=0.D0
	 IF(J.GT.1 .AND. MH.LT.ccgg08(Nccgg08,1))
     C	  MAtest=ccgg08(J-1,2)+(MH-ccgg08(J-1,1))/
     C	   (ccgg08(J,1)-ccgg08(J-1,1))
     C	   *(ccgg08(J,2)-ccgg08(J-1,2))
	 PROB(19)=PROB(19)+5.D0*(Rmax-ceff)*DDIM(MAtest/MA,1.D0)
	 GOTO 712

 612	 Rmax=1.D0
	 J=1
	 DOWHILE(ccgg1(J,1).LE.MH .AND. J.LT.Nccgg1)
	  J=J+1
	 ENDDO
	 MAtest=0.D0
	 IF(J.GT.1 .AND. MH.LT.ccgg1(Nccgg1,1))
     C	  MAtest=ccgg1(J-1,2)+(MH-ccgg1(J-1,1))/
     C	   (ccgg1(J,1)-ccgg1(J-1,1))
     C	   *(ccgg1(J,2)-ccgg1(J-1,2))
	 PROB(19)=PROB(19)+5.D0*(Rmax-ceff)*DDIM(MAtest/MA,1.D0)
	
 712	 CONTINUE

*      AA -> cctautau

	 ceff=R*BRHAA(I,1)*BRCC(4)*BRLL(4)*2.D0
	 IF(ceff.GE.0.0.AND.ceff.LT.0.2)GOTO 122
	 IF(ceff.GE.0.2.AND.ceff.LT.0.4)GOTO 222
	 IF(ceff.GE.0.4.AND.ceff.LT.0.5)GOTO 322
	 IF(ceff.GE.0.5.AND.ceff.LT.0.6)GOTO 422
	 IF(ceff.GE.0.6.AND.ceff.LT.0.8)GOTO 522
	 IF(ceff.GE.0.8.AND.ceff.LT.1.0)GOTO 622
	 GOTO 722

 122	 Rmax=.2D0
	 J=1
	 DOWHILE(cctt02(J,1).LE.MH .AND. J.LT.Ncctt02)
	  J=J+1
	 ENDDO
	 MAtest=0.D0
	 IF(J.GT.1 .AND. MH.LT.cctt02(Ncctt02,1))
     C	  MAtest=cctt02(J-1,2)+(MH-cctt02(J-1,1))/
     C	   (cctt02(J,1)-cctt02(J-1,1))
     C	   *(cctt02(J,2)-cctt02(J-1,2))
	 PROB(19)=PROB(19)+5.D0*(Rmax-ceff)*DDIM(MAtest/MA,1.D0)
	 GOTO 722

 222	 Rmax=.4D0
	 J=1
	 DOWHILE(cctt04(J,1).LE.MH .AND. J.LT.Ncctt04)
	  J=J+1
	 ENDDO
	 MAtest=0.D0
	 IF(J.GT.1 .AND. MH.LT.cctt04(Ncctt04,1))
     C	  MAtest=cctt04(J-1,2)+(MH-cctt04(J-1,1))/
     C	   (cctt04(J,1)-cctt04(J-1,1))
     C	   *(cctt04(J,2)-cctt04(J-1,2))
	 PROB(19)=PROB(19)+5.D0*(Rmax-ceff)*DDIM(MAtest/MA,1.D0)
	 GOTO 722

 322	 Rmax=.5D0
	 J=1
	 DOWHILE(cctt05(J,1).LE.MH .AND. J.LT.Ncctt05)
	  J=J+1
	 ENDDO
	 MAtest=0.D0
	 IF(J.GT.1 .AND. MH.LT.cctt05(Ncctt05,1))
     C	  MAtest=cctt05(J-1,2)+(MH-cctt05(J-1,1))/
     C	   (cctt05(J,1)-cctt05(J-1,1))
     C	   *(cctt05(J,2)-cctt05(J-1,2))
	 PROB(19)=PROB(19)+10.D0*(Rmax-ceff)*DDIM(MAtest/MA,1.D0)
	 GOTO 722

 422	 Rmax=.6D0
	 J=1
	 DOWHILE(cctt06(J,1).LE.MH .AND. J.LT.Ncctt06)
	  J=J+1
	 ENDDO
	 MAtest=0.D0
	 IF(J.GT.1 .AND. MH.LT.cctt06(Ncctt06,1))
     C	  MAtest=cctt06(J-1,2)+(MH-cctt06(J-1,1))/
     C	   (cctt06(J,1)-cctt06(J-1,1))
     C	   *(cctt06(J,2)-cctt06(J-1,2))
	 PROB(19)=PROB(19)+10.D0*(Rmax-ceff)*DDIM(MAtest/MA,1.D0)
	 GOTO 722

 522	 Rmax=.8D0
	 J=1
	 DOWHILE(cctt08(J,1).LE.MH .AND. J.LT.Ncctt08)
	  J=J+1
	 ENDDO
	 MAtest=0.D0
	 IF(J.GT.1 .AND. MH.LT.cctt08(Ncctt08,1))
     C	  MAtest=cctt08(J-1,2)+(MH-cctt08(J-1,1))/
     C	   (cctt08(J,1)-cctt08(J-1,1))
     C	   *(cctt08(J,2)-cctt08(J-1,2))
	 PROB(19)=PROB(19)+5.D0*(Rmax-ceff)*DDIM(MAtest/MA,1.D0)
	 GOTO 722

 622	 Rmax=1.D0
	 J=1
	 DOWHILE(cctt1(J,1).LE.MH .AND. J.LT.Ncctt1)
	  J=J+1
	 ENDDO
	 MAtest=0.D0
	 IF(J.GT.1 .AND. MH.LT.cctt1(Ncctt1,1))
     C	  MAtest=cctt1(J-1,2)+(MH-cctt1(J-1,1))/
     C	   (cctt1(J,1)-cctt1(J-1,1))
     C	   *(cctt1(J,2)-cctt1(J-1,2))
	 PROB(19)=PROB(19)+5.D0*(Rmax-ceff)*DDIM(MAtest/MA,1.D0)
	
 722	 CONTINUE

*      AA -> jjjj

	 ceff=R*BRHAA(I,1)*(BRJJ(4)+BRSS(4))**2
	 IF(ceff.GE.0.0.AND.ceff.LT.0.2)GOTO 132
	 IF(ceff.GE.0.2.AND.ceff.LT.0.4)GOTO 232
	 IF(ceff.GE.0.4.AND.ceff.LT.0.5)GOTO 332
	 IF(ceff.GE.0.5.AND.ceff.LT.0.6)GOTO 432
	 IF(ceff.GE.0.6.AND.ceff.LT.0.8)GOTO 532
	 IF(ceff.GE.0.8.AND.ceff.LT.1.0)GOTO 632
	 GOTO 732

 132	 Rmax=.2D0
	 J=1
	 DOWHILE(gggg02(J,1).LE.MH .AND. J.LT.Ngggg02)
	  J=J+1
	 ENDDO
	 MAtest=0.D0
	 IF(J.GT.1 .AND. MH.LT.gggg02(Ngggg02,1))
     C	  MAtest=gggg02(J-1,2)+(MH-gggg02(J-1,1))/
     C	   (gggg02(J,1)-gggg02(J-1,1))
     C	   *(gggg02(J,2)-gggg02(J-1,2))
	 PROB(19)=PROB(19)+5.D0*(Rmax-ceff)*DDIM(MAtest/MA,1.D0)
	 GOTO 732

 232	 Rmax=.4D0
	 J=1
	 DOWHILE(gggg04(J,1).LE.MH .AND. J.LT.Ngggg04)
	  J=J+1
	 ENDDO
	 MAtest=0.D0
	 IF(J.GT.1 .AND. MH.LT.gggg04(Ngggg04,1))
     C	  MAtest=gggg04(J-1,2)+(MH-gggg04(J-1,1))/
     C	   (gggg04(J,1)-gggg04(J-1,1))
     C	   *(gggg04(J,2)-gggg04(J-1,2))
	 PROB(19)=PROB(19)+5.D0*(Rmax-ceff)*DDIM(MAtest/MA,1.D0)
	 GOTO 732

 332	 Rmax=.5D0
	 J=1
	 DOWHILE(gggg05(J,1).LE.MH .AND. J.LT.Ngggg05)
	  J=J+1
	 ENDDO
	 MAtest=0.D0
	 IF(J.GT.1 .AND. MH.LT.gggg05(Ngggg05,1))
     C	  MAtest=gggg05(J-1,2)+(MH-gggg05(J-1,1))/
     C	   (gggg05(J,1)-gggg05(J-1,1))
     C	   *(gggg05(J,2)-gggg05(J-1,2))
	 PROB(19)=PROB(19)+10.D0*(Rmax-ceff)*DDIM(MAtest/MA,1.D0)
	 GOTO 732

 432	 Rmax=.6D0
	 J=1
	 DOWHILE(gggg06(J,1).LE.MH .AND. J.LT.Ngggg06)
	  J=J+1
	 ENDDO
	 MAtest=0.D0
	 IF(J.GT.1 .AND. MH.LT.gggg06(Ngggg06,1))
     C	  MAtest=gggg06(J-1,2)+(MH-gggg06(J-1,1))/
     C	   (gggg06(J,1)-gggg06(J-1,1))
     C	   *(gggg06(J,2)-gggg06(J-1,2))
	 PROB(19)=PROB(19)+10.D0*(Rmax-ceff)*DDIM(MAtest/MA,1.D0)
	 GOTO 732

 532	 Rmax=.8D0
	 J=1
	 DOWHILE(gggg08(J,1).LE.MH .AND. J.LT.Ngggg08)
	  J=J+1
	 ENDDO
	 MAtest=0.D0
	 IF(J.GT.1 .AND. MH.LT.gggg08(Ngggg08,1))
     C	  MAtest=gggg08(J-1,2)+(MH-gggg08(J-1,1))/
     C	   (gggg08(J,1)-gggg08(J-1,1))
     C	   *(gggg08(J,2)-gggg08(J-1,2))
	 PROB(19)=PROB(19)+5.D0*(Rmax-ceff)*DDIM(MAtest/MA,1.D0)
	 GOTO 732

 632	 Rmax=1.D0
	 J=1
	 DOWHILE(gggg1(J,1).LE.MH .AND. J.LT.Ngggg1)
	  J=J+1
	 ENDDO
	 MAtest=0.D0
	 IF(J.GT.1 .AND. MH.LT.gggg1(Ngggg1,1))
     C	  MAtest=gggg1(J-1,2)+(MH-gggg1(J-1,1))/
     C	   (gggg1(J,1)-gggg1(J-1,1))
     C	   *(gggg1(J,2)-gggg1(J-1,2))
	 PROB(19)=PROB(19)+5.D0*(Rmax-ceff)*DDIM(MAtest/MA,1.D0)
	
 732	 CONTINUE

*      AA -> tautaujj

	 ceff=R*BRHAA(I,1)*BRLL(4)*(BRJJ(4)+BRSS(4))*2.D0
	 IF(ceff.GE.0.0.AND.ceff.LT.0.2)GOTO 142
	 IF(ceff.GE.0.2.AND.ceff.LT.0.4)GOTO 242
	 IF(ceff.GE.0.4.AND.ceff.LT.0.5)GOTO 342
	 IF(ceff.GE.0.5.AND.ceff.LT.0.6)GOTO 442
	 IF(ceff.GE.0.6.AND.ceff.LT.0.8)GOTO 542
	 IF(ceff.GE.0.8.AND.ceff.LT.1.0)GOTO 642
	 GOTO 742

 142	 Rmax=.2D0
	 J=1
	 DOWHILE(ttgg02(J,1).LE.MH .AND. J.LT.Nttgg02)
	  J=J+1
	 ENDDO
	 MAtest=0.D0
	 IF(J.GT.1 .AND. MH.LT.ttgg02(Nttgg02,1))
     C	  MAtest=ttgg02(J-1,2)+(MH-ttgg02(J-1,1))/
     C	   (ttgg02(J,1)-ttgg02(J-1,1))
     C	   *(ttgg02(J,2)-ttgg02(J-1,2))
	 PROB(19)=PROB(19)+5.D0*(Rmax-ceff)*DDIM(MAtest/MA,1.D0)
	 GOTO 742

 242	 Rmax=.4D0
	 J=1
	 DOWHILE(ttgg04(J,1).LE.MH .AND. J.LT.Nttgg04)
	  J=J+1
	 ENDDO
	 MAtest=0.D0
	 IF(J.GT.1 .AND. MH.LT.ttgg04(Nttgg04,1))
     C	  MAtest=ttgg04(J-1,2)+(MH-ttgg04(J-1,1))/
     C	   (ttgg04(J,1)-ttgg04(J-1,1))
     C	   *(ttgg04(J,2)-ttgg04(J-1,2))
	 PROB(19)=PROB(19)+5.D0*(Rmax-ceff)*DDIM(MAtest/MA,1.D0)
	 GOTO 742

 342	 Rmax=.5D0
	 J=1
	 DOWHILE(ttgg05(J,1).LE.MH .AND. J.LT.Nttgg05)
	  J=J+1
	 ENDDO
	 MAtest=0.D0
	 IF(J.GT.1 .AND. MH.LT.ttgg05(Nttgg05,1))
     C	  MAtest=ttgg05(J-1,2)+(MH-ttgg05(J-1,1))/
     C	   (ttgg05(J,1)-ttgg05(J-1,1))
     C	   *(ttgg05(J,2)-ttgg05(J-1,2))
	 PROB(19)=PROB(19)+10.D0*(Rmax-ceff)*DDIM(MAtest/MA,1.D0)
	 GOTO 742

 442	 Rmax=.6D0
	 J=1
	 DOWHILE(ttgg06(J,1).LE.MH .AND. J.LT.Nttgg06)
	  J=J+1
	 ENDDO
	 MAtest=0.D0
	 IF(J.GT.1 .AND. MH.LT.ttgg06(Nttgg06,1))
     C	  MAtest=ttgg06(J-1,2)+(MH-ttgg06(J-1,1))/
     C	   (ttgg06(J,1)-ttgg06(J-1,1))
     C	   *(ttgg06(J,2)-ttgg06(J-1,2))
	 PROB(19)=PROB(19)+10.D0*(Rmax-ceff)*DDIM(MAtest/MA,1.D0)
	 GOTO 742

 542	 Rmax=.8D0
	 J=1
	 DOWHILE(ttgg08(J,1).LE.MH .AND. J.LT.Nttgg08)
	  J=J+1
	 ENDDO
	 MAtest=0.D0
	 IF(J.GT.1 .AND. MH.LT.ttgg08(Nttgg08,1))
     C	  MAtest=ttgg08(J-1,2)+(MH-ttgg08(J-1,1))/
     C	   (ttgg08(J,1)-ttgg08(J-1,1))
     C	   *(ttgg08(J,2)-ttgg08(J-1,2))
	 PROB(19)=PROB(19)+5.D0*(Rmax-ceff)*DDIM(MAtest/MA,1.D0)
	 GOTO 742

 642	 Rmax=1.D0
	 J=1
	 DOWHILE(ttgg1(J,1).LE.MH .AND. J.LT.Nttgg1)
	  J=J+1
	 ENDDO
	 MAtest=0.D0
	 IF(J.GT.1 .AND. MH.LT.ttgg1(Nttgg1,1))
     C	  MAtest=ttgg1(J-1,2)+(MH-ttgg1(J-1,1))/
     C	   (ttgg1(J,1)-ttgg1(J-1,1))
     C	   *(ttgg1(J,2)-ttgg1(J-1,2))
	 PROB(19)=PROB(19)+5.D0*(Rmax-ceff)*DDIM(MAtest/MA,1.D0)
	
 742	 CONTINUE

*      AA -> tautautautau

	 ceff=R*BRHAA(I,1)*BRLL(4)**2
	 IF(ceff.GE.0.0.AND.ceff.LT.0.2)GOTO 152
	 IF(ceff.GE.0.2.AND.ceff.LT.0.4)GOTO 252
	 IF(ceff.GE.0.4.AND.ceff.LT.0.5)GOTO 352
	 IF(ceff.GE.0.5.AND.ceff.LT.0.6)GOTO 452
	 IF(ceff.GE.0.6.AND.ceff.LT.0.8)GOTO 552
	 IF(ceff.GE.0.8.AND.ceff.LT.1.0)GOTO 652
	 GOTO 752

 152	 Rmax=.2D0
	 J=1
	 DOWHILE(tttt02(J,1).LE.MH .AND. J.LT.Ntttt02)
	  J=J+1
	 ENDDO
	 MAtest=0.D0
	 IF(J.GT.1 .AND. MH.LT.tttt02(Ntttt02,1))
     C	  MAtest=tttt02(J-1,2)+(MH-tttt02(J-1,1))/
     C	   (tttt02(J,1)-tttt02(J-1,1))
     C	   *(tttt02(J,2)-tttt02(J-1,2))
	 PROB(19)=PROB(19)+5.D0*(Rmax-ceff)*DDIM(MAtest/MA,1.D0)
	 GOTO 752

 252	 Rmax=.4D0
	 J=1
	 DOWHILE(tttt04(J,1).LE.MH .AND. J.LT.Ntttt04)
	  J=J+1
	 ENDDO
	 MAtest=0.D0
	 IF(J.GT.1 .AND. MH.LT.tttt04(Ntttt04,1))
     C	  MAtest=tttt04(J-1,2)+(MH-tttt04(J-1,1))/
     C	   (tttt04(J,1)-tttt04(J-1,1))
     C	   *(tttt04(J,2)-tttt04(J-1,2))
	 PROB(19)=PROB(19)+5.D0*(Rmax-ceff)*DDIM(MAtest/MA,1.D0)
	 GOTO 752

 352	 Rmax=.5D0
	 J=1
	 DOWHILE(tttt05(J,1).LE.MH .AND. J.LT.Ntttt05)
	  J=J+1
	 ENDDO
	 MAtest=0.D0
	 IF(J.GT.1 .AND. MH.LT.tttt05(Ntttt05,1))
     C	  MAtest=tttt05(J-1,2)+(MH-tttt05(J-1,1))/
     C	   (tttt05(J,1)-tttt05(J-1,1))
     C	   *(tttt05(J,2)-tttt05(J-1,2))
	 PROB(19)=PROB(19)+10.D0*(Rmax-ceff)*DDIM(MAtest/MA,1.D0)
	 GOTO 752

 452	 Rmax=.6D0
	 J=1
	 DOWHILE(tttt06(J,1).LE.MH .AND. J.LT.Ntttt06)
	  J=J+1
	 ENDDO
	 MAtest=0.D0
	 IF(J.GT.1 .AND. MH.LT.tttt06(Ntttt06,1))
     C	  MAtest=tttt06(J-1,2)+(MH-tttt06(J-1,1))/
     C	   (tttt06(J,1)-tttt06(J-1,1))
     C	   *(tttt06(J,2)-tttt06(J-1,2))
	 PROB(19)=PROB(19)+10.D0*(Rmax-ceff)*DDIM(MAtest/MA,1.D0)
	 GOTO 752

 552	 Rmax=.8D0
	 J=1
	 DOWHILE(tttt08(J,1).LE.MH .AND. J.LT.Ntttt08)
	  J=J+1
	 ENDDO
	 MAtest=0.D0
	 IF(J.GT.1 .AND. MH.LT.tttt08(Ntttt08,1))
     C	  MAtest=tttt08(J-1,2)+(MH-tttt08(J-1,1))/
     C	   (tttt08(J,1)-tttt08(J-1,1))
     C	   *(tttt08(J,2)-tttt08(J-1,2))
	 PROB(19)=PROB(19)+5.D0*(Rmax-ceff)*DDIM(MAtest/MA,1.D0)
	 GOTO 752

 652	 Rmax=1.D0
	 J=1
	 DOWHILE(tttt1(J,1).LE.MH .AND. J.LT.Ntttt1)
	  J=J+1
	 ENDDO
	 MAtest=0.D0
	 IF(J.GT.1 .AND. MH.LT.tttt1(Ntttt1,1))
     C	  MAtest=tttt1(J-1,2)+(MH-tttt1(J-1,1))/
     C	   (tttt1(J,1)-tttt1(J-1,1))
     C	   *(tttt1(J,2)-tttt1(J-1,2))
	 PROB(19)=PROB(19)+5.D0*(Rmax-ceff)*DDIM(MAtest/MA,1.D0)
	
 752	 CONTINUE

	ENDIF
	ENDDO

* Associated production

	DO I=1,3
	DO J=1,2
	MH=SMASS(I)
	MA=PMASS(J)
	R=PCOMP(J,1)**2*(SCOMP(I,1)*h2-SCOMP(I,2)*h1)**2
     C	  /(h1**2+h2**2)

*  Z width

	IF(MH+MA.LT.MZ)THEN
	 GZ=SQR2*GF*R*LAMBDA(MZ**2,MA**2,MH**2)**3/(48.D0*PI*MZ**3)
 	 PROB(13)=PROB(13)+DDIM(GZ/GZMAX,1.D0)
	ENDIF

*  ee -> hA -> 4bs

	IF(MH+MA.LT.SQRS)THEN

	 M1=MIN(MH,MA)
	 M2=MAX(MH,MA)

	 Rmax=1.D0
	 DO K=1,NhA4b
	  IF(AINT(M2).EQ.hA4b(K,1) .AND. AINT(M1).EQ.hA4b(K,2))THEN
	   Rmax=hA4b(K,3)
	   GOTO 3
	  ENDIF
	 ENDDO
 3	 PROB(14)=PROB(14)+DDIM(R*BRBB(I)*BRBB(J+3)/Rmax,1.D0)

*  ee -> hA -> 4taus

	 Rmax=1.D0
	 DO K=1,NhA4tau
	  IF(AINT(M2).EQ.hA4tau(K,1) .AND. AINT(M1).EQ.hA4tau(K,2))THEN
	   Rmax=hA4tau(K,3)
	   GOTO 4
	  ENDIF
	 ENDDO
 4	 PROB(15)=PROB(15)+DDIM(R*BRLL(I)*BRLL(J+3)/Rmax,1.D0)

*  ee -> hA -> 2b 2taus

	 Rmax=1.D0
	 IF(MA.GT.MH)THEN
	  DO K=1,NhA2b2tau
	   IF(AINT(M2).EQ.hA2b2tau(K,1) .AND.
     C	     AINT(M1).EQ.hA2b2tau(K,2))THEN
	    Rmax=hA2b2tau(K,3)
	    GOTO 5
	   ENDIF
	  ENDDO
 5	  PROB(16)=PROB(16)+DDIM(R*BRLL(I)*BRBB(J+3)/Rmax,1.D0)
	 ENDIF

	 Rmax=1.D0
	 IF(MH.GT.MA)THEN
	  DO K=1,NhA2b2tau
	   IF(AINT(M2).EQ.hA2b2tau(K,1) .AND.
     C	     AINT(M1).EQ.hA2b2tau(K,2))THEN
	    Rmax=hA2b2tau(K,3)
	    GOTO 6
	   ENDIF
	  ENDDO
 6	  PROB(16)=PROB(16)+DDIM(R*BRBB(I)*BRLL(J+3)/Rmax,1.D0)
	 ENDIF

	 Rmax=1.D0
	 IF(MA.GT.MH)THEN
	  DO K=1,NhA2tau2b
	   IF(AINT(M2).EQ.hA2tau2b(K,1) .AND.
     C	     AINT(M1).EQ.hA2tau2b(K,2))THEN
	    Rmax=hA2tau2b(K,3)
	    GOTO 7
	   ENDIF
	  ENDDO
 7	  PROB(16)=PROB(16)+DDIM(R*BRBB(I)*BRLL(J+3)/Rmax,1.D0)
	 ENDIF

	 Rmax=1.D0
	 IF(MH.GT.MA)THEN
	  DO K=1,NhA2tau2b
	   IF(AINT(M2).EQ.hA2tau2b(K,1) .AND.
     C	     AINT(M1).EQ.hA2tau2b(K,2))THEN
	    Rmax=hA2tau2b(K,3)
	    GOTO 8
	   ENDIF
	  ENDDO
 8	  PROB(16)=PROB(16)+DDIM(R*BRLL(I)*BRBB(J+3)/Rmax,1.D0)
	 ENDIF

*  ee -> hA -> AAA -> 6bs

	 IF(MH.GT.2.D0*MA)THEN
	  Rmax=1.D0
	  DO K=1,NAAA6b
	   IF(AINT(MH).EQ.AAA6b(K,1) .AND. AINT(MA).EQ.AAA6b(K,2))THEN
	    Rmax=AAA6b(K,3)
	    GOTO 9
	   ENDIF
	  ENDDO
 9	  PROB(17)=PROB(17)+DDIM(R*BRHAA(I,2*J-1)*BRBB(J+3)**3/Rmax,
     C	           1.D0)

*  ee -> hA -> AAA -> 6taus

	  Rmax=1.D0
	  DO K=1,NAAA6tau
	   IF(AINT(MH).EQ.AAA6tau(K,1) .AND.
     C	     AINT(MA).EQ.AAA6tau(K,2))THEN
	    Rmax=AAA6tau(K,3)
	    GOTO 10
	   ENDIF
	  ENDDO
 10	  PROB(18)=PROB(18)+DDIM(R*BRHAA(I,2*J-1)*BRLL(J+3)**3/Rmax,
     C	           1.D0)
	ENDIF

	ENDIF
	ENDDO
	ENDDO

	END


	DOUBLE PRECISION FUNCTION CLEOTAU(MX)

*  CLEO constraints on BR(Y -> A gamma)*BR(A -> tau tau)

	IMPLICIT NONE

	INTEGER I,N
	PARAMETER (N=17)
	DOUBLE PRECISION MX,X(N),M(N)

	DATA M/3.75D0,4.25D0,4.75D0,5.1D0,5.8D0,6.15D0,6.6D0,7.D0,
     C	      7.4D0,7.6D0,8.D0,8.25D0,8.6D0,9.D0,9.25D0,9.35D0,9.41D0/
	DATA X/2.9D-5,2.5D-5,2.D-5,2.3D-5,5.1D-5,2.5D-5,2.5D-5,2.7D-5,
     C	      4.5D-5,3.7D-5,2.7D-5,7.2D-5,6.8D-5,8.6D-5,2.1D-4,2.85D-4,
     C	      4.75D-4/

	CLEOTAU=0.D0

	IF(MX.LT.M(1).OR.MX.GT.M(N))RETURN

	DO I=2,N
	 IF(MX.LT.M(I))THEN
	  CLEOTAU=X(I-1)+(MX-M(I-1))/(M(I)-M(I-1))*(X(I)-X(I-1))
	  RETURN
	 ENDIF
	ENDDO

	END


	DOUBLE PRECISION FUNCTION CLEOMU(MX)

*  CLEO constraints on BR(Y -> A gamma)*BR(A -> mu mu)

	IMPLICIT NONE

	INTEGER I,N
	PARAMETER (N=2)
	DOUBLE PRECISION MX,X(N),M(N)

	DATA M/.25D0,3.75D0/
	DATA X/9.D-6,9.D-6/

	CLEOMU=0.D0

	IF(MX.LT.M(1).OR.MX.GT.M(N))RETURN

	DO I=2,N
	 IF(MX.LT.M(I))THEN
	  CLEOMU=X(I-1)+(MX-M(I-1))/(M(I)-M(I-1))*(X(I)-X(I-1))
	  RETURN
	 ENDIF
	ENDDO

	END
