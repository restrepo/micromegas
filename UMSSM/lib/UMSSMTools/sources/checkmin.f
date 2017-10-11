      SUBROUTINE CHECKMIN(PROB)

**********************************************************************      
* Subroutine to check whether the physical minimum of the effective
* potential (<hd>, <hu> and <s> =/= 0) is deeper than minima with
* <hd>, <hu> or <s> = 0
*
* If not: PROB(28) =/= 0
*
* For the soft masses squared mh1q, mh2q and msq
* (corrections at the scale QSTSB) :
* If mh1q, mh2q >> Q2 then PROB(29) =/= 0
*
**********************************************************************

      IMPLICIT NONE

      DOUBLE PRECISION PROB(*)
      DOUBLE PRECISION V,V1,V2,V3,mh1q,mh2q,msq
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION EPS,Q2,pi
      DOUBLE PRECISION LQ,KQ,ALQ,AKQ,MUEFFQ,NUQ 
      DOUBLE PRECISION UPARF,SAZZ,CAZZ,VEV,NCP,VEVS,G1P
      DOUBLE PRECISION QD,QU,QS,QQ,QUP,QDOW,QL,QE,QN
      DOUBLE PRECISION TANB,AU,AD,MU,M2
      DOUBLE PRECISION SST,SSB,LAMBDA,AL,ATAU
      DOUBLE PRECISION la1,la2,la3,la4,la5,la6,la7
      DOUBLE PRECISION aa5,la1s,la2s,cosb,sinb
      DOUBLE PRECISION vsq1,vsq2,vsq3

      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/QMHIGGS/MH1Q,MH2Q,MSQ
      COMMON/RENSCALE/Q2
      COMMON/QNMPAR/LQ,KQ,ALQ,AKQ,MUEFFQ,NUQ  
      COMMON/UMSSM/SAZZ,CAZZ,VEV,NCP,QD,QU,QS,VEVS,G1P,QQ,
     .		QUP,QDOW,QL,QE,QN
      COMMON/NOBUG/TANB,AU,AD,M2,SST,SSB,LAMBDA,AL,ATAU

      pi=4d0*DATAN(1d0)
      EPS=-1d-2
      sinb=tanb/dsqrt(1d0+tanb**2)
      cosb=sinb/tanb

*   Soft masses

      MH1Q=UPARF(345)
      MH2Q=UPARF(346)
      MSQ=UPARF(347)

*   Effective couplings

      la1=UPARF(348)
      la2=UPARF(349)
      la3=UPARF(350)
      la4=UPARF(351)
      la5=UPARF(352)
      la6=UPARF(353)
      la7=UPARF(354)
      aa5=UPARF(355)
      la1s=UPARF(356)
      la2s=UPARF(357)

*   Physical minimum

      V= (MUEFFQ*VEV)**2/2d0 + (LAMBDA*sinb*cosb*VEV**2/2d0)**2
     .  - AL*MUEFFQ*sinb*cosb*VEV**2
     .  + (g1+g2)/8d0*VEV**4*(sinb**2-cosb**2)**2/4d0
     .  + mh1q*(VEV*cosb)**2/2d0 + mh2q*(VEV*sinb)**2/2d0
     .  + msq*VEVS**2/2d0
     .  + G1P**2/2d0*(QD*(VEV*cosb)**2+QU*(VEV*sinb)**2
     .  + QS*VEVS**2)**2/4d0
     .  + la1*(VEV*cosb)**4/8d0 + la2*(VEV*sinb)**4/8d0
     .  + (la3+la4+la5)*VEV**4*(cosb*sinb)**2/4d0
     .  - (la6*(VEV*cosb)**2+la7*(VEV*sinb)**2)*VEV**2*cosb*sinb/2d0
     .  - aa5*VEVS*VEV**2*cosb*sinb/DSQRT(2d0)
     .  + (la1s*(VEV*cosb)**2+la2s*(VEV*sinb)**2)*VEVS**2/4d0

*   Minimum with h2=s=0
*   Min. V1 -> get vd**2 (vsq1) as a function of the other params. :

      vsq1=-mh1q/((g1+g2)/8d0 + (G1P*QD)**2/2d0 + la1/2d0)
      V1= mh1q*vsq1/2d0 + (g1+g2)/8d0*vsq1**2/4d0
     .  + (G1P*QD*vsq1)**2/8d0 + la1*vsq1**2/8d0 
      IF(V.NE.0d0)THEN
       PROB(28)=DDIM(EPS,(V1-V)/DABS(V))
      ELSE
       PROB(28)=DDIM(EPS,V1)
      ENDIF

*   Minimum with h1=s=0
*   Min. V2 -> get vu**2 (vsq2) as a function of the other params. :

      vsq2=-mh2q/((g1+g2)/8d0 + (G1P*QU)**2/2d0 + la2/2d0)
      V2= mh2q*vsq2/2d0 + (g1+g2)/8d0*vsq2**2/4d0
     .  + (G1P*QU*vsq2)**2/8d0 + la2*vsq2**2/8d0 
      IF(V.NE.0d0)THEN
       PROB(28)=PROB(28)+DDIM(EPS,(V2-V)/DABS(V))
      ELSE
       PROB(28)=PROB(28)+DDIM(EPS,V2)
      ENDIF

*   Minimum with h1=h2=0
*   Min. V3 -> get vs**2 (vsq3) as a function of the other params. :

      vsq3=-2d0*msq/(G1P*QS)**2
      V3= msq*vsq3/2d0 + (G1P*QS*vsq3)**2/8d0
      IF(V.NE.0d0)THEN
       PROB(28)=PROB(28)+DDIM(EPS,(V3-V)/DABS(V))
      ELSE
       PROB(28)=PROB(28)+DDIM(EPS,V3)
      ENDIF


*      WRITE(0,*)"CALL CHECKMIN"
*      WRITE(0,*)""
*      WRITE(*,121)"V  =",V
*      WRITE(*,121)"V1 =",V1
*      WRITE(*,121)"V2 =",V2
*      WRITE(*,121)"V3 =",V3
*      WRITE(*,121)"PROB(28) =",PROB(28)
*      WRITE(*,121)"mh1q =",mh1q
*      WRITE(*,121)"mh2q =",mh2q
*      WRITE(*,121)"msq =",msq
*      WRITE(0,*)""
*      WRITE(0,*)""

      PROB(29)=DDIM(MAX(DABS(MH1Q),DABS(MH2Q))/Q2,10d0)
      
 121  FORMAT(A,E20.3)

      END
