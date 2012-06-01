	SUBROUTINE MHIGGS(PAR,IFAIL)

***********************************************************************
*
*	This subroutine computes the Higgs masses and couplings in
*	the NMSSM. The relevant parameters are read from
*       COMMON/QGAUGE, /QHIGGS, /QQUARK, /QNMPAR
*       (computed in RUNPAR) and  /GMSUSYPAR
*       The squark masses and mixing parameters are read from
*       COMMON/RADCOR (computed in MSFERM)
*
*	On output: 
*
*	SMASS(1-3): CP-even masses (ordered)
*
*	SCOMP(1-3,1-3): Mixing angles: if HB(I) are the bare states,
*	  HB(I) = Re(H1), Re(H2), Re(S), and HM(I) are the mass eigenstates, 
*	  the convention is HB(I) = SUM_(J=1,3) SCOMP(J,I)*HM(J)
*	  which is equivalent to HM(I) = SUM_(J=1,3) SCOMP(I,J)*HB(J)
*
*	PMASS(1-2): CP-odd masses (ordered)
*
*	PCOMP(1-2,1-2): Mixing angles: if AB(I) are the bare states,
*	  AB(I) = Im(H1), Im(H2), Im(S), and AM(I) are the mass eigenstates, 
*	  the convention is 
*	  AM(I) = PCOMP(I,1)*(COSBETA*AB(1)+SINBETA*AB(2))
*			+ PCOMP(I,2)*AB(3)
*
*	CMASS: Charged Higgs mass
*
*	IFAIL:  =   0	        OK
*		=   1,3,5,7	SMASS(1)**2 < 0
*		=   2,3,6,7	PMASS(1)**2 < 0
*		=   4,5,6,7	CMASS**2 < 0
*
*	The precision in the computation of the lightest Higgs mass is:
*
*	Terms ~ ht^4/hb^4 from (s)top/(s)bottom-loops are computed
*	exactly in the mixing parameters.
*
*	Terms ~ g^2*(ht^2/hb^2) (where g is an electro-weak gauge coupling)
*	are taken into account due to the wave function renormalizations
*	factors, finite self energies (pole masses) and corrections from
*       stop/sbottom D terms
*
*	Leading logs ~ (g, l, k)^4 are added explicitely.
*
*	Leading double logs from two loops ~ ht/hb^6 and ht/hb^4*alpha_s 
*	are taken into account.
*
*	For heavy higgses (with masses mhh > mtop) the leading log
*	contributions ~ (ht^2/hb^2)*log(mhh/mtop) to the pole masses
*	are taken into account, but not the corresponding effects on the
*	mixing angles.
*
*	All mixing angles are at the scale m_top.
*
*	The dominant errors come from one loop terms ~ (g,l,k)^4 without
*	large logs, and from two loop terms without large double logs.
*
************************************************************************

	IMPLICIT NONE

	INTEGER IFAIL,I,J

	DOUBLE PRECISION PAR(*),VEC3(3,3),VEC2(2,2),EPS
	DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2),CMASS
	DOUBLE PRECISION H(3,3),A(2,2),MH(3),MA(2),MHC
	DOUBLE PRECISION Alshift,B,X,BT,BB,MGAU
	DOUBLE PRECISION mst1,mst2,s2t,msb1,msb2,s2b,XT,XB,M1,M2,T
	DOUBLE PRECISION At,emt,fmt,gmt,Ab,emb,fmb,gmb
	DOUBLE PRECISION QSTSB,PI,C2TW,sferm,bos
	DOUBLE PRECISION sb,cb,s2,rt,rb,MS2,MP2,M12
	DOUBLE PRECISION LM2,Lmu,Lnu,LM2mu,Lmunu,LQZ,LA,LP,LS,LPP,P1,P2
	DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW,gb,gt,Db,Dt,D1,D2
	DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW,subdet
   	DOUBLE PRECISION G1Q,G2Q,GQ,ALSQ
	DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
	DOUBLE PRECISION HTQ,HBQ,MTOPQ,MBOTQ
	DOUBLE PRECISION LQ,KQ,ALQ,AKQ,MUQ,NUQ
	DOUBLE PRECISION MA2,COEF,GAUGE,htau
	DOUBLE PRECISION MHH,MAA,MSS,MHA,MHS,MAS,MPP,MPS,MPM
	DOUBLE PRECISION XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY
	DOUBLE PRECISION XIF,XIS,MUP,MSP,M3H
	DOUBLE PRECISION GMCOMB

	COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
	COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
	COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
	COMMON/RADCOR/mst1,mst2,s2t,msb1,msb2,s2b,XT,XB
 	COMMON/STSBSCALE/QSTSB	
	COMMON/QGAUGE/G1Q,G2Q,GQ,ALSQ
	COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
	COMMON/QQUARK/HTQ,HBQ,MTOPQ,MBOTQ
	COMMON/QNMPAR/LQ,KQ,ALQ,AKQ,MUQ,NUQ
	COMMON/ALSHIFT/ALSHIFT
 	COMMON/HOUT/MHH,MAA,MSS,MHA,MHS,MAS,MPP,MPS,MPM
	COMMON/HMO/MH,MA,MHC
	COMMON/GMSUSYPAR/XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY

	EPS=1.D-8
	PI=4.D0*DATAN(1.D0)
	COEF=1.D0/(16.D0*PI**2)
	
	IFAIL=0

*   Trig. functions of beta
	cb=1.D0/DSQRT(1.D0+tanbq**2)
	sb=tanbq*cb
	s2=2.D0*sb*cb

	M1=PAR(20)
	M2=PAR(21)
	MA2=PAR(23)**2
	At=PAR(12)
	Ab=PAR(13)
	htau=MTAU/H2Q

*   Identify the GMSB parameters at the Susy scale with QSTSB:
	XIF=XIFSUSY
	XIS=XISSUSY
	MUP=MUPSUSY
	MSP=MSPSUSY
	M3H=M3HSUSY

*   Weak angle theta_W (S2TW= sin(theta_W)**2):
	C2TW=1.D0-S2TW

*   Approximate value for the Singlet-like CP even Higgs mass squared:
	MS2=MAX(nuq*(4.D0*nuq+Akq),MZ**2)

*   Approximate value for the Singlet-like CP odd Higgs mass squared:
	MP2=-3.D0*nuq*Akq
     C   -XIF*(4.D0*KQ+LQ*MUP/MUQ)-2.D0*MSP-MUP*NUQ-LQ*XIS/MUQ

*   Approximate value for the off-diag. CP odd mass matrix element:
	M12=lq*DSQRT(h1q**2+h2q**2)*(Alq-2.D0*nuq)

*   One loop functions for stop/sbottom loop corrections

	rt= 3.D0/2.D0*COEF*htq**2
	IF(mst1.NE.mst2)THEN
	 fmt= (mst2*DLOG(mst2/QSTSB)-mst1*DLOG(mst1/QSTSB))/
     C     (mst2-mst1)-1.D0
	 gmt= s2t**2*((mst2+mst1)/(mst2-mst1)*DLOG(mst2/mst1)-2.D0)
	 emt= -mtopq*s2t*DLOG(mst2/mst1)
	ELSE
	 fmt= DLOG(mst1/QSTSB)
	 gmt= 0.D0
	 emt= 0.D0
	ENDIF

	rb= 3.D0/2.D0*COEF*hbq**2
	IF(msb1.NE.msb2)THEN
	 fmb= (msb2*DLOG(msb2/QSTSB)-msb1*DLOG(msb1/QSTSB))/
     C    (msb2-msb1)-1.D0
	 gmb= s2b**2*((msb2+msb1)/(msb2-msb1)*DLOG(msb2/msb1)-2.D0)
	 emb= -mbotq*s2b*DLOG(msb2/msb1)
	ELSE
	 fmb= DLOG(msb1/QSTSB)
	 gmb= 0.D0
	 emb= 0.D0
	ENDIF
	
*  The subsequent shifts in Alambda simplify the expressions for the 
*  one loop rad. corrs. below.
*  The parameter Alq is defined at the scale QSTSB

	IF(mst1.NE.mst2)THEN
         Alshift= Alq+2.D0*rt*At*
     C 	  ((mst2*DLOG(mst2/QSTSB)-mst1*DLOG(mst1/QSTSB))/
     C     (mst2-mst1)-1.D0)
 	ELSE
	 Alshift= Alq+2.D0*rt*At*DLOG(mst1/QSTSB)
	ENDIF

	IF(msb1.NE.msb2)THEN
         Alshift= Alshift+2.D0*rb*Ab*
     C 	  ((msb2*DLOG(msb2/QSTSB)-msb1*DLOG(msb1/QSTSB))/
     C     (msb2-msb1)-1.D0)
	ELSE
	 Alshift= Alshift+2.D0*rb*Ab*DLOG(msb1/QSTSB)
	ENDIF
	
	MGAU= (G1q*M1+3.D0*G2q*M2)
	LM2mu= DLOG(MAX(M2**2,muq**2,MZ**2)/QSTSB)

	Alshift= Alshift+MGAU*LM2MU*COEF
		
	B= Alshift+nuq

*   Tree level CP-even Higgs mass matrix

	GMCOMB=Lq*XIF+MUP*MUQ+M3H

	H(1,1) = gq*h1q**2 + (muq*B+GMCOMB)/tanbq
	H(2,2) = gq*h2q**2 + (muq*B+GMCOMB)*tanbq
	H(3,3) = Lq**2*(Alshift+MUP)*h1q*h2q/muq 
     C   + nuq*(AKQ+4.D0*NUQ+3.D0*MUP)-LQ/MUQ*(XIS+XIF*MUP)
	H(1,2) = (2.D0*Lq**2-gq)*h1q*h2q - muq*B-GMCOMB
	H(1,3) = Lq*(2.D0*muq*h1q - (B+nuq+MUP)*h2q)
	H(2,3) = Lq*(2.D0*muq*h2q - (B+nuq+MUP)*h1q)
	
*   1-loop radiative corrections

	H(1,1)= H(1,1) + rt*(4.D0*At*emt - At**2*gmt
     C		+ 4.D0*mtopq**2*DLOG(mst1*mst2/mtopq**4))
	H(2,2)= H(2,2) - rt*muq**2*gmt 
	H(3,3)= H(3,3) - rt*lq**2*h2q**2*gmt
	H(1,2)= H(1,2) + rt*muq*(At*gmt - 2.D0*emt)	
	H(1,3)= H(1,3) + rt*lq*h2q*(At*gmt - 2.D0*emt)
	H(2,3)= H(2,3) + rt*(4.D0*lq*h2q*muq*fmt - lq*h2q*muq*gmt)

	H(1,1)= H(1,1) - rb*muq**2*gmb
	H(2,2)= H(2,2) + rb*(4.D0*Ab*emb - Ab**2*gmb
     C		+ 4.D0*mbotq**2*DLOG(msb1*msb2/mbotq**4))
	H(3,3)= H(3,3) - rb*lq**2*h1q**2*gmb	
	H(1,2)= H(1,2) + rb*muq*(Ab*gmb - 2.D0*emb)	
	H(1,3)= H(1,3) + rb*(4.D0*lq*h1q*muq*fmb - lq*h1q*muq*gmb)
	H(2,3)= H(2,3) + rb*lq*h1q*(Ab*gmb - 2.D0*emb)
	
*   Corrections from higgs/stop/sbottom couplings from D-terms

	gt= gq/2.D0-2.D0*g1q/3.D0
	Dt= (PAR(7)-PAR(8)+gt*(h2q**2-h1q**2))/2.D0
	D1= 3.D0*COEF*(-gq/4.D0*emt + gt/2.D0*Dt/Xt*gmt)
	D2= 3.D0*COEF*(-gt*Dt/Xt*emt
     C	    -gq/2.D0*mtopq**2*DLOG(mst1*mst2/QSTSB**2))
	H(1,1)= H(1,1) + 2.D0*At*D1 + 2.D0*D2
	H(2,2)= H(2,2) + 2.D0*muq/tanbq*D1
	H(1,2)= H(1,2) - 1.D0/tanbq*D2 - (muq+At/tanbq)*D1
	H(1,3)= H(1,3) - lq*h2q*D1
	H(2,3)= H(2,3) + lq*h2q/tanbq*D1

	gb= gq/2.D0-g1q/3.D0
	Db= (PAR(7)-PAR(9)+gb*(h1q**2-h2q**2))/2.D0
	D1= 3.D0*COEF*(-gq/4.D0*emb + gb/2.D0*Db/Xb*gmb)
	D2= 3.D0*COEF*(-gb*Db/Xb*emb
     C	    -gq/2.D0*mbotq**2*DLOG(msb1*msb2/QSTSB**2))
	H(1,1)= H(1,1) + 2.D0*muq*tanbq*D1
	H(2,2)= H(2,2) + 2.D0*Ab*D1 + 2.D0*D2
	H(1,2)= H(1,2) - tanbq*D2 - (muq+Ab*tanbq)*D1
	H(1,3)= H(1,3) + lq*h1q*tanbq*D1
	H(2,3)= H(2,3) - lq*h1q*D1

*   2-loop terms

	T= DLOG(QSTSB/MT**2)
     
	H(1,1)= H(1,1) + 6.D0*COEF**2*htq**4*h1q**2*
     C	 (T**2*(64.D0*PI*ALSQ+4.D0/3.D0*g1q-3.D0*htq**2*sb**2
     C	 +3.D0*hbq**2*cb**2)+((dlog(MAX(QSTSB,MA2,MT**2)/MT**2))**2
     C	 -(dlog(MAX(MA2,MT**2)/MT**2))**2)*
     C	 (-3.D0*htq**2*cb**2-hbq**2*(3.D0*cb**2+1.D0)))
     
	H(2,2)= H(2,2) + 6*COEF**2*hbq**4*h2q**2*
     C   (T**2*(64.D0*PI*ALSQ-2.D0/3.D0*g1q-3.D0*hbq**2*cb**2
     C	 +3.D0*htq**2*sb**2)+((dlog(MAX(QSTSB,MA2,MT**2)/MT**2))**2
     C	 -(dlog(MAX(MA2,MT**2)/MT**2))**2)*
     C	 (-3.D0*hbq**2*sb**2-htq**2*(3.D0*sb**2+1.D0)))
   
*   Leading-log electroweak contributions (1 loop):
 
*   a) Sfermion contributions

	sferm= 2.D0*COEF*gq*MZ**2*
     C   ((S2TW**2/6.D0+C2TW**2*(3.D0/2.D0))
     C	 *DLOG(MAX(PAR(7),MZ**2)/QSTSB)
     C   +4.D0/3.D0*S2TW**2*DLOG(MAX(PAR(8),MZ**2)/QSTSB)
     C   +1.D0/3.D0*S2TW**2*DLOG(MAX(PAR(9),MZ**2)/QSTSB)
     C   +(S2TW**2/3.D0+C2TW**2*3.D0)*DLOG(MAX(PAR(15),MZ**2)/QSTSB)
     C   +8.D0/3.D0*S2TW**2*DLOG(MAX(PAR(16),MZ**2)/QSTSB)
     C   +2.D0/3.D0*S2TW**2*DLOG(MAX(PAR(17),MZ**2)/QSTSB)
     C   +(S2TW**2/2.D0+C2TW**2/2.D0)*DLOG(MAX(PAR(10),MZ**2)/QSTSB)
     C   +S2TW**2*DLOG(MAX(PAR(11),MZ**2)/QSTSB)
     C   +(S2TW**2+C2TW**2)*DLOG(MAX(PAR(18),MZ**2)/QSTSB)
     C   +2.D0*S2TW**2*DLOG(MAX(PAR(19),MZ**2)/QSTSB))

	H(1,1)= H(1,1) + sferm*sb**2
	H(2,2)= H(2,2) + sferm*cb**2
	H(1,2)= H(1,2) - sferm*sb*cb

*    b) Chargino/neutralino contributions
 
	LQZ= DLOG(QSTSB/MZ**2)
	LM2= DLOG(MAX(M2**2,MZ**2)/QSTSB)
	Lmu= DLOG(MAX(muq**2,MZ**2)/QSTSB)
	Lnu= DLOG(MAX(2.D0*(2.D0*nuq+MUP)**2,MZ**2)/QSTSB)
	Lmunu= DLOG(MAX(muq**2,2.D0*(2.D0*nuq+MUP)**2,MZ**2)/QSTSB)

	H(1,1)= H(1,1) + COEF*
     C	 ((MZ**2*SB**2*Gq*(-20.D0+32.D0*S2TW-16.D0*S2TW**2))*LM2MU
     C	 -4.D0*(MUq*(NUq+MUP/2.D0)*Lq**2/TANBQ+MZ**2*SB**2*Lq**4/Gq)
     C	 *LMUNU)

	H(2,2)= H(2,2) + COEF*
     C	 ((MZ**2*CB**2*Gq*(-20.D0+32.D0*S2TW-16.D0*S2TW**2))*LM2MU
     C	 -4.D0*(MUq*(NUq+MUP/2.D0)*Lq**2*TANBQ+MZ**2*CB**2*Lq**4/Gq)
     C	 *LMUNU)

	H(3,3)= H(3,3) + COEF*
     C	 (4.D0*KQ**2*(-8.D0*NUq**2-6.D0*NUq*MUP+MUP**3/NUq)*LNU
     C    -8.*Lq**2*MUq**2*LMU
     C    +Lq**3*MZ**2*MUP/(MUQ*GQ)*(4.D0*KQ-LQ*S2)*Lmunu)

	H(1,2)= H(1,2) + COEF*
     C	 (4.D0*(MUq*(NUq+MUP/2.D0)*Lq**2-MZ**2*SB*CB*Lq**4/Gq)*LMUNU
     C	 -(4.D0*Gq*MZ**2*SB*CB)*LM2MU)

	H(1,3)= H(1,3) + COEF*(MZ/DSQRT(Gq))*
     C	 (Lq*Gq*MUq*SB*(-12.D0+8.D0*S2TW)*LM2MU
     C	 +4.D0*Lq**2*(CB*LQ*(2.D0*NUQ+MUP/2.0D0)-SB*(LQ*MUQ
     C   +4.D0*KQ*(NUQ+MUP/2.0D0)))*LMUNU)
 
	H(2,3)= H(2,3) + COEF*(MZ/DSQRT(Gq))*
     C	 (Lq*Gq*MUq*CB*(-12.D0+8.D0*S2TW)*LM2MU
     C	 +4.D0*Lq**2*(SB*LQ*(2.D0*NUQ+MUP/2.0D0)-CB*(LQ*MUQ
     C   +4.D0*KQ*(NUQ+MUP/2.0D0)))*LMUNU)

*    c) Higgs loop contributions
*    (Only if all masses squared are positive, and only to the lighter
*    CP-even doublet-like state)

	subdet= H(3,3)*(sb**2*H(1,1)+cb**2*H(2,2)+s2*H(1,2))
     C	 -(sb*H(1,3)+cb*H(2,3))**2

	IF(subdet.GT.0.D0)THEN

	 P2= MAX(MA2+MP2,MZ**2)
	 P1= MAX((MA2*MP2-M12**2)/P2,MZ**2)
	 LA= DLOG(MAX(MA2,MZ**2)/MZ**2)
	 LS= DLOG(MS2/MZ**2)
	 LP= DLOG(P2/MZ**2)
	 LPP= DLOG(P2/P1)

	 bos= COEF*MZ**2/gq*((gq**2*(2.D0*S2TW**2
     C	 -2.D0*S2TW*(1.D0+s2**2)-11.D0/4.D0*s2**4+5.D0*s2**2+3.D0/4.D0)
     C	 +gq*lq**2*(2.D0*S2TW*s2**2+11.D0/2.D0*s2**4-15.D0/2.D0*s2**2
     C	 -1.D0)+lq**4*(-11.D0/4.D0*s2**4+5.D0/2.D0*s2**2+1.D0))*LA
     C	 +(lq**2*(lq-kq*s2)**2+3.D0*lq**2/MS2*(gq+(lq**2-gq)*s2**2)
     C   *(2.D0*muq-s2*(Alq+(2.D0*nuq+MUP)))**2
     C	 -lq**4/MS2**2*(2.D0*muq-s2*(Alq+(2.D0*nuq+MUP)))**4)*LS
     C	 +(gq**2/4.D0*(1.D0-s2**4)+gq*lq**2*(1.D0/2.D0*s2**4
     C	 +1.D0/2.D0*s2**2-1.D0)+lq**4*(-1.D0/4.D0*s2**4-1.D0/2.D0*s2**2
     C	 +1.D0)+lq**2*(lq+kq*s2)**2)*LP
     C	 -((gq-lq**2)**2/2.D0*MP2/P2*s2**2*(1.D0-s2**2)
     C	 +(lq*MA2*(lq+kq*s2)-lq**2*(Alq-(2.D0*nuq+MUP))**2
     C   -MP2/2.D0*(gq*(1.D0-s2**2)
     C	 -lq**2*(2.D0-s2**2)))**2/P2**2-lq*MA2*MP2*(lq+kq*s2)
     C	 *(gq*(1.D0-s2**2)-lq**2*(2.D0-s2**2))/P2**2)*LPP
     C   +(Gq**2*(-4.D0+S2**2+2.D0*S2TW*(1.D0+S2**2)-2.D0*S2TW**2)
     C   +Gq*Lq**2*(2.D0+S2**2-2.D0*S2**2*S2TW)
     C   -Lq**4*(4.D0+2.D0*S2**2)-2.D0*Lq**2*Kq**2*S2**2)*LQZ)

	 H(1,1)= H(1,1) + bos*sb**2
	 H(2,2)= H(2,2) + bos*cb**2
	 H(1,2)= H(1,2) + bos*sb*cb

	ENDIF

*   d) Gauge Loop Corrections

	GAUGE= MZ**2*Gq*COEF*(-9.D0+12.D0*S2TW-6.D0*S2TW**2)*LQZ
	
	H(1,1)= H(1,1) + GAUGE*sb**2
	H(2,2)= H(2,2) + GAUGE*cb**2
	H(1,2)= H(1,2) + GAUGE*sb*cb

*   Take care of the Z factors

	H(1,1)= H(1,1)/ZHU
	H(2,2)= H(2,2)/ZHD
	H(3,3)= H(3,3)/ZS
	H(1,2)= H(1,2)/DSQRT(ZHU*ZHD)
	H(1,3)= H(1,3)/DSQRT(ZHU*ZS)
	H(2,3)= H(2,3)/DSQRT(ZHD*ZS)

*   Diagonalization

	MHH= sb**2*H(1,1)+cb**2*H(2,2)+2.D0*cb*sb*H(1,2)
	MAA= cb**2*H(1,1)+sb**2*H(2,2)-2.D0*cb*sb*H(1,2)
	MSS= H(3,3)
	MHA= (H(1,1)-H(2,2))*cb*sb+(cb**2-sb**2)*H(1,2)
	MHS= sb*H(1,3)+cb*H(2,3)
	MAS= cb*H(1,3)-sb*H(2,3)
	MHH= MHH/DSQRT(DABS(MHH))
	MAA= MAA/DSQRT(DABS(MAA))
	MSS= MSS/DSQRT(DABS(MSS))
	MHA= MHA/DSQRT(DABS(MHA))
	MHS= MHS/DSQRT(DABS(MHS))
	MAS= MAS/DSQRT(DABS(MAS))

	CALL DIAGN(3,H,MH,VEC3,EPS)
	CALL SORTN(3,MH,VEC3)
	DO I= 1,3
	 DO J= 1,3
	  SCOMP(I,J)= VEC3(J,I)
	 ENDDO
	ENDDO

*   CP even pole masses

	DO I= 1,3
	IF(MH(I).GT.0.D0)THEN
	 IF(MH(I).GT.4.D0*MT**2)THEN
	  X= DSQRT(1.D0-4.D0*MT**2/MH(I))
	  BT= 2.D0-X*DLOG((1.D0+X)/(1.D0-X))
	 ELSE
	  X= DSQRT(MH(I)/(4.D0*MT**2-MH(I)))
	  BT= 2.D0*(1.D0-DATAN(X)/X)
	 ENDIF
	 IF(MH(I).GT.4.D0*MB**2)THEN
	  X= DSQRT(1.D0-4.D0*MB**2/MH(I))
	  BB= 2.D0-X*DLOG((1.D0+X)/(1.D0-X))
	 ELSE
	  X= DSQRT(MH(I)/(4.D0*MB**2-MH(I)))
	  BB= 2.D0*(1.D0-DATAN(X)/X)
	 ENDIF
	 SMASS(I)= MH(I)
     C	  - 2.D0*rt*SCOMP(I,1)**2*(MH(I)-4.D0*MT**2)*BT
     C	  - 2.D0*rb*SCOMP(I,2)**2*(MH(I)*DLOG(MT**2/MB**2)
     C	   +(MH(I)-4.D0*MB**2)*BB)
	 MH(I)= DSQRT(MH(I))
	ELSE
	 SMASS(I)= MH(I)
	ENDIF
	ENDDO

	IF(SMASS(1).LT.0.D0)THEN
	 IFAIL= IFAIL+1
	ELSE
	 DO I= 1,3
	  SMASS(I)= DSQRT(SMASS(I))
	 ENDDO
	ENDIF

*   CP-odd Higgs mass matrix including
*   1-loop top/bottom radiative corrections

	A(1,1)= (MUQ*B+GMCOMB)*(h1q/(ZHD*h2q)+h2q/(ZHU*h1q))
	A(2,2)= LQ**2*h1q*h2q/MUQ*(3.D0*NUQ+B)-3.D0*AKQ*NUQ
     C	 -XIF*(4.D0*KQ+LQ*MUP/MUQ)-2.D0*MSP-MUP*NUQ-LQ*XIS/MUQ
	A(1,2)= LQ*(Alshift-2.D0*NUQ-MUP)
     C	 *DSQRT(h1q**2/ZHD+h2q**2/ZHU)

*   Diagonalization

	MPP= A(1,1)/DSQRT(DABS(A(1,1)))
	MPS= A(2,2)/DSQRT(DABS(A(2,2)))
	MPM= A(1,2)/DSQRT(DABS(A(1,2)))

	CALL DIAGN(2,A,MA,VEC2,EPS)
	CALL SORTN(2,MA,VEC2)
	DO I= 1,2
	 DO J= 1,2
	  PCOMP(I,J)= VEC2(J,I)
	 ENDDO
	ENDDO

*   CP odd pole masses

	DO I= 1,2
	IF(MA(I).GT.0.D0)THEN
	 IF(MA(I).GT.4.D0*MT**2)THEN
	  X= DSQRT(1.D0-4.D0*MT**2/MA(I))
	  BT= 2.D0-X*DLOG((1.D0+X)/(1.D0-X))
	 ELSE
	  X= DSQRT(MA(I)/(4.D0*MT**2-MA(I)))
	  BT= 2.D0*(1.D0-DATAN(X)/X)
	 ENDIF
	 IF(MA(I).GT.4.D0*MB**2)THEN
	  X= DSQRT(1.D0-4.D0*MB**2/MA(I))
	  BB= 2.D0-X*DLOG((1.D0+X)/(1.D0-X))
	 ELSE
	  X= DSQRT(MA(I)/(4.D0*MB**2-MA(I)))
	  BB= 2.D0*(1.D0-DATAN(X)/X)
	 ENDIF
	 PMASS(I)= MA(I)
     C	  - 2.D0*rt*PCOMP(I,1)**2*cb**2*MA(I)*BT
     C	  - 2.D0*rb*PCOMP(I,1)**2*sb**2*(MA(I)*DLOG(MT**2/MB**2)
     C	   +MA(I)*BB)
	 MA(I)= DSQRT(MA(I))
	ELSE
	 PMASS(I)= MA(I)
	ENDIF
	ENDDO

	IF(PMASS(1).LT.0.D0)THEN
	 IFAIL= IFAIL+2
	ELSE
	 DO I= 1,2
	  PMASS(I)= DSQRT(PMASS(I))
	 ENDDO
	ENDIF

*   Charged Higgs mass including 1-loop radiative corrections

        MHC= ((G2Q/2.D0-LQ**2)*h1q*h2q+MUQ*B+GMCOMB)
     C	    *(h1q**2*ZHU+h2q**2*ZHD)/(h1q*h2q*ZHU*ZHD)
     C	    +3.D0*COEF*htq**2*hbq**2/(DSQRT(2.D0)*GF)*t
     C	    +COEF/3.D0*G2Q*MW**2*(12.D0*t+3.D0*(5*S2TW/C2TW-1.D0)*LM2mu
     C	    -4.D0*LM2-2.D0*Lmu)
	
*   Charged Higgs pole mass (in the LLA only)

	X=MHC/MT**2
	IF(X.EQ.0.D0)THEN
	 BT=1.D0
	ELSEIF(X.EQ.1.D0)THEN
	 BT=0.D0
	ELSE
	 BT=(1.D0-1.D0/X)*DLOG(DABS(X-1.D0))
	ENDIF
	CMASS= MHC+3.D0*COEF*((htq**2*cb**2+hbq**2*sb**2)
     C	      *(MHC*(BT-2.D0) - (MT**2+MB**2)*(BT-1.D0))
     C	      - 4.D0*HTQ*MT*HBQ*MB*cb*sb*(BT-1.D0))

	IF(MHC.GT.0)THEN
	 MHC= DSQRT(MHC)
	ENDIF

	IF(CMASS.LT.0.D0)THEN
	 IFAIL= IFAIL+4
	ELSE
	 CMASS= DSQRT(CMASS)
	ENDIF

	END


