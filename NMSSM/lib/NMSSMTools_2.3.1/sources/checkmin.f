	SUBROUTINE CHECKMIN(PAR,PROB)

**********************************************************************	
* Subroutine to check whether the physical minimum of the effective
* potential (<h1>, <h2> and <s> =/= 0) is deeper than minima with
* <h1>, <h2> or <s> = 0
*
* If not: PROB(28) =/= 0
*
* The effective potential includes 1 loop contributions (large logs 
* + finite) from (s)top and (s)bottom loops
*
* The soft masses squared mh1, mh2 and mss (at the scale QSTSB)
* are computed here and stored in COMMON/QMHIGGS (in NMHDECAY)
* or directly taken from COMMON/QMHIGGS (in NMSPEC and NMGMSB)
* If mh1, mh2 >> Q2 then PROB(29) =/= 0
*
**********************************************************************

	IMPLICIT NONE

	INTEGER OMGFLAG,MAFLAG

	DOUBLE PRECISION PAR(*),PROB(*)
	DOUBLE PRECISION QSTSB,pi,cc,mQ3,mU3,mD3,At,Ab
	DOUBLE PRECISION mst1,mst2,s2t,msb1,msb2,s2b,XT,XB
	DOUBLE PRECISION Mstop,dMst,Wt
	DOUBLE PRECISION Msbot,dMsb,Wb
	DOUBLE PRECISION ct,fmt1,fmt2,fmt,gmt
	DOUBLE PRECISION cb,fmb1,fmb2,fmb,gmb
	DOUBLE PRECISION V,V0,V1,V2,V3,D,mh1,mh2,mss
	DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
	DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
	DOUBLE PRECISION Bq,EPS
	DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
	DOUBLE PRECISION HTQ,HBQ,MTOPQ,MBOTQ
	DOUBLE PRECISION LQ,KQ,ALQ,AKQ,MUQ,NUQ
	DOUBLE PRECISION G1Q,G2Q,GQ,ALSQ,Q2
	DOUBLE PRECISION S1,S2,S3,S4,S5,S6,S7
	DOUBLE PRECISION XIF,XIS,MUP,MSP,M3H
	DOUBLE PRECISION A,B,C,P,Q,DET,AUX,XM
	DOUBLE PRECISION V31,V32,PHI,SIGQ,R

	COMMON/FLAGS/OMGFLAG,MAFLAG
	COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
	COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
	COMMON/STSBSCALE/QSTSB
	COMMON/QMHIGGS/MH1,MH2,MSS
	COMMON/RADCOR/mst1,mst2,s2t,msb1,msb2,s2b,XT,XB
	COMMON/RADCOR2/MQ3,MU3,MD3,AT,AB
	COMMON/QGAUGE/G1Q,G2Q,GQ,ALSQ
	COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
	COMMON/QQUARK/HTQ,HBQ,MTOPQ,MBOTQ
	COMMON/QNMPAR/LQ,KQ,ALQ,AKQ,MUQ,NUQ
	COMMON/RENSCALE/Q2
	COMMON/GMSUSYPAR/XIF,XIS,MUP,MSP,M3H

	pi=4.D0*DATAN(1.D0)
	cc=3.D0/(32.D0*pi**2)
	EPS=-1.D-2
	
* Store previous values:
	S1=MST1
	S2=MST2
	S3=MSB1
	S4=MSB2
	S5=MTOPQ
	S6=MBOTQ
	S7=NUQ
	
*   Parameters for (S)top/(S)bottom rad. corrs. in wrong minima:

	Mstop=(mQ3+mU3)/2.D0
	dMst=(mQ3-mU3)/2.D0
	Msbot=(mQ3+mD3)/2.D0
	dMsb=(mQ3-mD3)/2.D0
	
	IF(MIN(mst1,msb1).LE.0.D0)RETURN

	ct=3.D0*htq**2/(16.D0*pi**2)
	fmt1=mst1*(DLOG(mst1/QSTSB)-1.D0)
	fmt2=mst2*(DLOG(mst2/QSTSB)-1.D0)
	fmt=mtopq**2*(DLOG(mtopq**2/QSTSB)-1.D0)
	IF(mst1-mst2.NE.0.D0)THEN
	 gmt=(fmt2-fmt1)/(mst2-mst1)
	ELSE
	 gmt=DLOG(mst1/QSTSB)
	ENDIF

	cb=3.D0*hbq**2/(16.D0*pi**2)
	fmb1=msb1*(DLOG(msb1/QSTSB)-1.D0)
	fmb2=msb2*(DLOG(msb2/QSTSB)-1.D0)
	fmb=mbotq**2*(DLOG(mbotq**2/QSTSB)-1.D0)
	IF(msb1-msb2.NE.0.D0)THEN
	 gmb=(fmb2-fmb1)/(msb2-msb1)
	ELSE
	 gmb=DLOG(msb1/QSTSB)
	ENDIF

*   Soft masses

	IF(MAFLAG.GE.0)THEN
	 Bq=Alq+nuq
	 mh1= -lq**2*h2q**2 - muq**2 + muq*Bq/tanbq 
     C    + gq/2.D0*(h2q**2-h1q**2)
     C	  - ct*(fmt1+fmt2-2.D0*fmt+At*Xt*gmt)
     C	  + cb*muq/tanbq*Xb*gmb
	 mh2= -lq**2*h1q**2 - muq**2 + muq*Bq*tanbq 
     C    + gq/2.D0*(h1q**2-h2q**2)
     C	  + ct*muq*tanbq*Xt*gmt
     C	  - cb*(fmb1+fmb2-2.D0*fmb+Ab*Xb*gmb)
	 mss= -lq**2*(h1q**2+h2q**2) - 2.D0*nuq**2
     C	  + lq**2*h1q*h2q/muq*(Bq+nuq) - nuq*Akq
     C	  + ct*lq**2*h1q*h2q/muq*Xt*gmt
     C	  + cb*lq**2*h1q*h2q/muq*Xb*gmb
	ENDIF

*   Physical minimum

	V= muq**2*(h1q**2+h2q**2) + lq**2*h1q**2*h2q**2
     C	  + nuq**4/kq**2 + mup**2*muq**2/lq**2
     C	  - 2.D0*nuq*muq*h1q*h2q - 2.D0*lq*xif*h1q*h2q
     C	  - 2.D0*mup*muq*h1q*h2q + 2.D0*xif*nuq**2/kq
     C	  + 2.D0*mup*nuq**3/kq**2 + 2.D0*xif*mup*muq/lq   
     C	  - 2.D0*(Alq*muq+m3h)*h1q*h2q
     C    + 2.D0/3.D0*Akq*nuq**3/kq**2 
     C	  + 2.D0*xis*muq/lq + msp*muq**2/lq**2
     C    + gq/4.D0*(h1q**2-h2q**2)**2
     C	  + mh1*h1q**2 + mh2*h2q**2 + mss*muq**2/lq**2
     C	  + cc*
     C	  ( mst1**2*(DLOG(mst1/QSTSB)-1.5D0)
     C	  + mst2**2*(DLOG(mst2/QSTSB)-1.5D0)
     C	  - 2.D0*mtopq**4*(DLOG(mtopq**2/QSTSB)-1.5D0)
     C	  + msb1**2*(DLOG(msb1/QSTSB)-1.5D0)
     C	  + msb2**2*(DLOG(msb2/QSTSB)-1.5D0)
     C	  - 2.D0*mbotq**4*(DLOG(mbotq**2/QSTSB)-1.5D0))

*   Minimum with h1=h2=s=0

	mst1=mU3
	mst2=mQ3
	msb1=mD3
	msb2=mQ3
	V0= 0.D0
     C	  + cc*
     C	  ( mst1**2*(DLOG(mst1/QSTSB)-1.5D0)
     C	  + mst2**2*(DLOG(mst2/QSTSB)-1.5D0)
     C	  + msb1**2*(DLOG(msb1/QSTSB)-1.5D0)
     C	  + msb2**2*(DLOG(msb2/QSTSB)-1.5D0))
	IF(V.NE.0.D0)THEN
	 PROB(28)=DDIM(EPS,(V0-V)/DABS(V))
	ELSE
	 PROB(28)=DDIM(EPS,V0)
	ENDIF

*   Minimum with h2=s=0

	V1=0.D0
	IF(mh1.LT.0.D0)THEN
	 mtopq=htq*DSQRT(-2.D0*mh1/gq)
	 Wt=DSQRT(dMst**2+mtopq**2*At**2)
	 mst1=Mstop+mtopq**2-Wt
	 mst2=Mstop+mtopq**2+Wt
	 msb1=mD3
	 msb2=mQ3
	 IF(mst1.GT.0.D0)THEN
	  V1= -mh1**2/gq
     C	    + cc*
     C	    ( mst1**2*(DLOG(mst1/QSTSB)-1.5D0)
     C	    + mst2**2*(DLOG(mst2/QSTSB)-1.5D0)
     C	    - 2.D0*mtopq**4*(DLOG(mtopq**2/QSTSB)-1.5D0)
     C	    + msb1**2*(DLOG(msb1/QSTSB)-1.5D0)
     C	    + msb2**2*(DLOG(msb2/QSTSB)-1.5D0))
	  IF(V.NE.0.D0)THEN
	   PROB(28)=PROB(28)+DDIM(EPS,(V1-V)/DABS(V))
	  ELSE
	   PROB(28)=PROB(28)+DDIM(EPS,V1)
	  ENDIF
	 ENDIF
	ENDIF

*   Minimum with h1=s=0

	V2=0.D0
	IF(mh2.LT.0.D0)THEN
	 mst1=mU3
	 mst2=mQ3
	 mbotq=hbq*DSQRT(-2.D0*mh2/gq)
	 Wb=DSQRT(dMsb**2+mbotq**2*Ab**2)
	 msb1=Msbot+mbotq**2-Wb
	 msb2=Msbot+mbotq**2+Wb
	 IF(msb1.GT.0.D0)THEN
	  V2= -mh2**2/gq
     C	    + cc*
     C	    ( mst1**2*(DLOG(mst1/QSTSB)-1.5D0)
     C	    + mst2**2*(DLOG(mst2/QSTSB)-1.5D0)
     C	    + msb1**2*(DLOG(msb1/QSTSB)-1.5D0)
     C	    + msb2**2*(DLOG(msb2/QSTSB)-1.5D0)
     C	    - 2.D0*mbotq**4*(DLOG(mbotq**2/QSTSB)-1.5D0))
	  IF(V.NE.0.D0)THEN
	   PROB(28)=PROB(28)+DDIM(EPS,(V2-V)/DABS(V))
	  ELSE
	   PROB(28)=PROB(28)+DDIM(EPS,V2)
	  ENDIF
	 ENDIF
	ENDIF
 
 
*   Minimum with h1=h2=0
	A=4.D0*KQ**2
	B=KQ*(6.D0*MUP+2.*AKQ)
	C=2.D0*MUP**2+4.D0*KQ*XIF+2.D0*MSS+2.D0*MSP
	D=2.D0*MUP*XIF+2.D0*XIS
	P=C/(3.D0*A)-B**2/(9.D0*A**2)
	Q=B**3/(27.D0*A**3)-B*C/(6.D0*A**2)+D/(2.D0*A)
	DET=Q**2+P**3
	IF(Q.EQ.0.D0) Q=1.D0
	SIGQ=Q/DABS(Q)
	R=SIGQ*DSQRT(DABS(P))
	 mst1=mU3
	 mst2=mQ3
	 msb1=mD3
	 msb2=mQ3
	IF(DET.GE.0.D0) THEN
	  AUX=(DABS(Q)+DSQRT(DET))**(1.D0/3.D0)
	  XM=SIGQ*(P/AUX-AUX)-B/(3.D0*A)
	  V3=A*XM**4/4.D0+B*XM**3/3.D0+C*XM**2/2.D0+D*XM
     C	   + cc*
     C	   ( mst1**2*(DLOG(mst1/QSTSB)-1.5D0)
     C	   + mst2**2*(DLOG(mst2/QSTSB)-1.5D0)
     C	   + msb1**2*(DLOG(msb1/QSTSB)-1.5D0)
     C	   + msb2**2*(DLOG(msb2/QSTSB)-1.5D0))
     	 ELSE
	   PHI=DACOS(DABS(Q)/(DABS(P))**(3.D0/2.D0))
	   XM=-2.D0*R*DCOS(PHI/3.D0)-B/(3.D0*A)
	   V31=A*XM**4/4.D0+B*XM**3/3.D0+C*XM**2/2.D0+D*XM
     C	   + cc*
     C	   ( mst1**2*(DLOG(mst1/QSTSB)-1.5D0)
     C	   + mst2**2*(DLOG(mst2/QSTSB)-1.5D0)
     C	   + msb1**2*(DLOG(msb1/QSTSB)-1.5D0)
     C	   + msb2**2*(DLOG(msb2/QSTSB)-1.5D0))
	   XM=2.D0*R*DCOS((PHI-PI)/3.D0)-B/(3.D0*A)
	   V32=A*XM**4/4.D0+B*XM**3/3.D0+C*XM**2/2.D0+D*XM
     C	   + cc*
     C	   ( mst1**2*(DLOG(mst1/QSTSB)-1.5D0)
     C	   + mst2**2*(DLOG(mst2/QSTSB)-1.5D0)
     C	   + msb1**2*(DLOG(msb1/QSTSB)-1.5D0)
     C	   + msb2**2*(DLOG(msb2/QSTSB)-1.5D0))
     	   V3=MIN(V31,V32)
	 ENDIF
   	 IF(V.NE.0.D0)THEN
	  PROB(28)=PROB(28)+DDIM(EPS,(V3-V)/DABS(V))
	 ELSE
	  PROB(28)=PROB(28)+DDIM(EPS,V3)
	 ENDIF

	PROB(29)=DDIM(MAX(DABS(MH1),DABS(MH2))/Q2,10.D0)
	
* RETURN previous values:
	MST1=S1
	MST2=S2
	MSB1=S3
	MSB2=S4
	MTOPQ=S5
	MBOTQ=S6
	NUQ=S7

	!PRINT*,"CALL CHECKMIN"
	!PRINT*,""
	!PRINT*,"V =",V
	!PRINT*,"V0 =",V0
	!PRINT*,"V1 =",V1
	!PRINT*,"V2 =",V2
	!PRINT*,"V3 =",V3
	!PRINT*,""
	!PRINT*,""

	END


