*   Function to handle UMSSM parameters/terms

        DOUBLE PRECISION FUNCTION UPARF(R)
        IMPLICIT NONE
        DOUBLE PRECISION UPARCTOF
        INTEGER R
        EXTERNAL UPARCTOF
        UPARF = UPARCTOF(R)
        END


*   Function for DELMB

        DOUBLE PRECISION FUNCTION INTEG(X,Y,Z)
        IMPLICIT NONE
        DOUBLE PRECISION X,Y,Z
        IF(DABS(X).EQ.DABS(Y) .AND. DABS(X).EQ.DABS(Z))THEN
         INTEG=.5d0/X
        ELSEIF(DABS(X).EQ.DABS(Y))THEN
         INTEG=(X**2-Z**2+Z**2*DLOG(Z**2/X**2))/(X**2-Z**2)**2
        ELSEIF(DABS(Y).EQ.DABS(Z))THEN
         INTEG=(Y**2-X**2+X**2*DLOG(X**2/Y**2))/(Y**2-X**2)**2
        ELSEIF(DABS(X).EQ.DABS(Z))THEN
         INTEG=(X**2-Y**2+Y**2*DLOG(Y**2/X**2))/(X**2-Y**2)**2
        ELSE
         INTEG=(X**2*Y**2*DLOG(X**2/Y**2)
     .  +Y**2*Z**2*DLOG(Y**2/Z**2)+Z**2*X**2*DLOG(Z**2/X**2))/
     .  ((X**2-Y**2)*(Y**2-Z**2)*(X**2-Z**2))
        ENDIF
        END




        DOUBLE PRECISION FUNCTION RUNM(Q,NF)

*   Subroutine to calculate the quark running masses

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IMPLICIT INTEGER (I-N)

        PARAMETER (NN=6)
        PARAMETER (ZETA3 = 1.202056903159594D0)

        DIMENSION AM(NN),YMSB(NN)

        COMMON/ALS/XLAMBDA,AMCA,AMBA,AMTA,N0A
        COMMON/GAUGE/ALSMZ,ALEMMZ,GF,gg1,gg2,S2TW
        COMMON/SMSPEC/AMS,AMC,AMB,AMBP,AMT,AMTAU,AMMUON,AMZ,AMW

        B0(NF)= (33.D0-2.D0*NF)/12.D0
        B1(NF)= (102.D0-38.D0/3.D0*NF)/16.D0
        B2(NF)= (2857.D0/2.D0-5033.D0/18.D0*NF+325.D0/54.D0*NF**2)/64.D0
        G0(NF)= 1.D0
        G1(NF)= (202.D0/3.D0-20.D0/9.D0*NF)/16.D0
        G2(NF)= (1249.D0-(2216.D0/27.D0+160.D0/3.D0*ZETA3)*NF
     .  - 140.D0/81.D0*NF**2)/64.D0
        C1(NF)= G1(NF)/B0(NF) - B1(NF)*G0(NF)/B0(NF)**2
        C2(NF)= ((G1(NF)/B0(NF) - B1(NF)*G0(NF)/B0(NF)**2)**2
     .  + G2(NF)/B0(NF) + B1(NF)**2*G0(NF)/B0(NF)**3
     .  - B1(NF)*G1(NF)/B0(NF)**2 - B2(NF)*G0(NF)/B0(NF)**2)/2.D0
        TRAN(X,XK)= 1.D0+4.D0/3.D0*ALPHAS(X,2)/PI+XK*(ALPHAS(X,2)/PI)**2
        CQ(X,NF)= (2.D0*B0(NF)*X)**(G0(NF)/B0(NF))
     .    * (1.D0+C1(NF)*X+C2(NF)*X**2)

        PI= 4.D0*DATAN(1.D0)
        ACC= 1.D-8
        AM(1)= 0
        AM(2)= 0
        AM(3)= AMS
        AM(4)= AMC
        AM(5)= AMBP
        AM(6)= AMT
        XK= 16.11D0
        DO 1 I=1,NF-1
         XK= XK - 1.04D0*(1.D0-AM(I)/AM(NF))
1       CONTINUE
        IF(NF.GE.4)THEN
         XMSB= AM(NF)/TRAN(AM(NF),0D0)
         XMHAT= XMSB/CQ(ALPHAS(AM(NF),2)/PI,NF)
        ELSE
         XMSB= 0
         XMHAT= 0
        ENDIF
        YMSB(3)= AMS
        IF(NF.EQ.3)THEN
         YMSB(4)= YMSB(3)*CQ(ALPHAS(AM(4),2)/PI,3)/
     .                  CQ(ALPHAS(1.D0,2)/PI,3)
         YMSB(5)= YMSB(4)*CQ(ALPHAS(AM(5),2)/PI,4)/
     .                  CQ(ALPHAS(AM(4),2)/PI,4)
         YMSB(6)= YMSB(5)*CQ(ALPHAS(AM(6),2)/PI,5)/
     .                  CQ(ALPHAS(AM(5),2)/PI,5)
        ELSEIF(NF.EQ.4)THEN
         YMSB(4)= XMSB
         YMSB(5)= YMSB(4)*CQ(ALPHAS(AM(5),2)/PI,4)/
     .                  CQ(ALPHAS(AM(4),2)/PI,4)
         YMSB(6)= YMSB(5)*CQ(ALPHAS(AM(6),2)/PI,5)/
     .                  CQ(ALPHAS(AM(5),2)/PI,5)
        ELSEIF(NF.EQ.5)THEN
         YMSB(5)= XMSB
         YMSB(4)= YMSB(5)*CQ(ALPHAS(AM(4),2)/PI,4)/
     .                  CQ(ALPHAS(AM(5),2)/PI,4)
         YMSB(6)= YMSB(5)*CQ(ALPHAS(AM(6),2)/PI,5)/
     .                  CQ(ALPHAS(AM(5),2)/PI,5)
        ELSEIF(NF.EQ.6)THEN
         YMSB(6)= XMSB
         YMSB(5)= YMSB(6)*CQ(ALPHAS(AM(5),2)/PI,5)/
     .                  CQ(ALPHAS(AM(6),2)/PI,5)
         YMSB(4)= YMSB(5)*CQ(ALPHAS(AM(4),2)/PI,4)/
     .                  CQ(ALPHAS(AM(5),2)/PI,4)
        ENDIF
        IF(Q.LT.AMC)THEN
         N0= 3
         Q0= 1.D0
        ELSEIF(Q.LE.AMBP)THEN
         N0= 4
         Q0= AMC
        ELSEIF(Q.LE.AMT)THEN
         N0= 5
         Q0= AMBP
        ELSE
         N0= 6
         Q0= AMT
        ENDIF
        IF(NF.GT.3)THEN
         XKFAC= TRAN(AM(NF),0D0)/TRAN(AM(NF),XK)
        ELSE
         XKFAC= 1D0
        ENDIF
        RUNM= YMSB(N0)*CQ(ALPHAS(Q,2)/PI,N0)/
     .             CQ(ALPHAS(Q0,2)/PI,N0)
     .   * XKFAC

        RETURN
        END


*  Running alpha_s and aux. subroutines/functions as in HDECAY

        DOUBLE PRECISION FUNCTION ALPHAS(Q,N)

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IMPLICIT INTEGER (I-N)

        DIMENSION XLB(6)

        COMMON/ALSLAM/XLB1(6),XLB2(6)
        COMMON/ALS/XLAMBDA,AMC,AMBP,AMT,N0

        B0(NF)=33.D0-2.D0*NF
        B1(NF)=6.D0*(153.D0-19.D0*NF)/B0(NF)**2
        ALS1(NF,X)=12.D0*PI/(B0(NF)*DLOG(X**2/XLB(NF)**2))
        ALS2(NF,X)=12.D0*PI/(B0(NF)*DLOG(X**2/XLB(NF)**2))
     .      *(1.D0-B1(NF)*DLOG(DLOG(X**2/XLB(NF)**2))
     .       /DLOG(X**2/XLB(NF)**2))

  
        PI=4.D0*DATAN(1.D0)

        IF(N.EQ.1)THEN
         DO 1 I=1,6
          XLB(I)=XLB1(I)
1       CONTINUE
        ELSE
         DO 2 I=1,6
          XLB(I)=XLB2(I)
2       CONTINUE
        ENDIF
        IF(Q.LT.AMC)THEN
         NF=3
        ELSEIF(Q.LE.AMBP)THEN
         NF=4
        ELSEIF(Q.LE.AMT)THEN
         NF=5
        ELSE
         NF=6
        ENDIF
        IF(N.EQ.1)THEN
          ALPHAS=ALS1(NF,Q)
        ELSE
          ALPHAS=ALS2(NF,Q)
        ENDIF

        RETURN
        END


        SUBROUTINE ALSINI(ACC)

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IMPLICIT INTEGER (I-N)

        DIMENSION XLB(6)

        COMMON/ALSLAM/XLB1(6),XLB2(6)
        COMMON/ALS/XLAMBDA,AMC,AMBP,AMT,N0
        
        PI=4.D0*DATAN(1.D0)
        XLB1(1)=0.D0
        XLB1(2)=0.D0
        XLB2(1)=0.D0
        XLB2(2)=0.D0
        IF(N0.EQ.3)THEN
         XLB(3)=XLAMBDA
         XLB(4)=XLB(3)*(XLB(3)/AMC)**(2.D0/25.D0)
         XLB(5)=XLB(4)*(XLB(4)/AMBP)**(2.D0/23.D0)
         XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
        ELSEIF(N0.EQ.4)THEN
         XLB(4)=XLAMBDA
         XLB(5)=XLB(4)*(XLB(4)/AMBP)**(2.D0/23.D0)
         XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
         XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
        ELSEIF(N0.EQ.5)THEN
         XLB(5)=XLAMBDA
         XLB(4)=XLB(5)*(XLB(5)/AMBP)**(-2.D0/25.D0)
         XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
         XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
        ELSEIF(N0.EQ.6)THEN
         XLB(6)=XLAMBDA
         XLB(5)=XLB(6)*(XLB(6)/AMT)**(-2.D0/23.D0)
         XLB(4)=XLB(5)*(XLB(5)/AMBP)**(-2.D0/25.D0)
         XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
        ENDIF
        DO 1 I=1,6
         XLB1(I)=XLB(I)
1       CONTINUE
        IF(N0.EQ.3)THEN
         XLB(3)=XLAMBDA
         XLB(4)=XLB(3)*(XLB(3)/AMC)**(2.D0/25.D0)
     .           *(2.D0*DLOG(AMC/XLB(3)))**(-107.D0/1875.D0)
         XLB(4)=XITER(AMC,XLB(3),3,XLB(4),4,ACC)
         XLB(5)=XLB(4)*(XLB(4)/AMBP)**(2.D0/23.D0)
     .           *(2.D0*DLOG(AMBP/XLB(4)))**(-963.D0/13225.D0)
         XLB(5)=XITER(AMBP,XLB(4),4,XLB(5),5,ACC)
         XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
     .           *(2.D0*DLOG(AMT/XLB(5)))**(-321.D0/3381.D0)
         XLB(6)=XITER(AMT,XLB(5),5,XLB(6),6,ACC)
        ELSEIF(N0.EQ.4)THEN
         XLB(4)=XLAMBDA
         XLB(5)=XLB(4)*(XLB(4)/AMBP)**(2.D0/23.D0)
     .           *(2.D0*DLOG(AMBP/XLB(4)))**(-963.D0/13225.D0)
         XLB(5)=XITER(AMBP,XLB(4),4,XLB(5),5,ACC)
         XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
     .           *(2.D0*DLOG(AMC/XLB(4)))**(107.D0/2025.D0)
         XLB(3)=XITER(AMC,XLB(4),4,XLB(3),3,ACC)
         XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
     .           *(2.D0*DLOG(AMT/XLB(5)))**(-321.D0/3381.D0)
         XLB(6)=XITER(AMT,XLB(5),5,XLB(6),6,ACC)
        ELSEIF(N0.EQ.5)THEN
         XLB(5)=XLAMBDA
         XLB(4)=XLB(5)*(XLB(5)/AMBP)**(-2.D0/25.D0)
     .           *(2.D0*DLOG(AMBP/XLB(5)))**(963.D0/14375.D0)
         XLB(4)=XITER(AMBP,XLB(5),5,XLB(4),4,ACC)
         XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
     .           *(2.D0*DLOG(AMC/XLB(4)))**(107.D0/2025.D0)
         XLB(3)=XITER(AMC,XLB(4),4,XLB(3),3,ACC)
         XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
     .           *(2.D0*DLOG(AMT/XLB(5)))**(-321.D0/3381.D0)
         XLB(6)=XITER(AMT,XLB(5),5,XLB(6),6,ACC)
        ELSEIF(N0.EQ.6)THEN
         XLB(6)=XLAMBDA
         XLB(5)=XLB(6)*(XLB(6)/AMT)**(-2.D0/23.D0)
     .           *(2.D0*DLOG(AMT/XLB(6)))**(321.D0/3703.D0)
         XLB(5)=XITER(AMT,XLB(6),6,XLB(5),5,ACC)
         XLB(4)=XLB(5)*(XLB(5)/AMBP)**(-2.D0/25.D0)
     .           *(2.D0*DLOG(AMBP/XLB(5)))**(963.D0/14375.D0)
         XLB(4)=XITER(AMBP,XLB(5),5,XLB(4),4,ACC)
         XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
     .           *(2.D0*DLOG(AMC/XLB(4)))**(107.D0/2025.D0)
         XLB(3)=XITER(AMC,XLB(4),4,XLB(3),3,ACC)
        ENDIF
        DO 2 I=1,6
         XLB2(I)=XLB(I)
2       CONTINUE

        RETURN
        END


        DOUBLE PRECISION FUNCTION XITER(Q,XLB1,NF1,XLB,NF2,ACC)

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IMPLICIT INTEGER (I-N)

        B0(NF)=33.D0-2.D0*NF
        B1(NF)=6.D0*(153.D0-19.D0*NF)/B0(NF)**2
        ALS2(NF,X,XLB)=12.D0*PI/(B0(NF)*DLOG(X**2/XLB**2))
     .            *(1.D0-B1(NF)*DLOG(DLOG(X**2/XLB**2))
     .            /DLOG(X**2/XLB**2))
        AA(NF)=12D0*PI/B0(NF)
        BB(NF)=B1(NF)/AA(NF)
        XIT(A,B,X)=A/2.D0*(1D0+DSQRT(1D0-4D0*B*DLOG(X)))
        PI=4.D0*DATAN(1.D0)
        XLB2=XLB
        II=0
1       II=II+1
        X=DLOG(Q**2/XLB2**2)
        ALP=ALS2(NF1,Q,XLB1)
        A=AA(NF2)/ALP
        B=BB(NF2)*ALP
        XX=XIT(A,B,X)
        XLB2=Q*DEXP(-XX/2.D0)
        Y1=ALS2(NF1,Q,XLB1)
        Y2=ALS2(NF2,Q,XLB2)
        DY=DABS(Y2-Y1)/Y1
        IF(DY.GE.ACC) GOTO 1
        XITER=XLB2

        RETURN
        END


        DOUBLE PRECISION FUNCTION XITLA(NO,ALP,ACC)

*  Iteration routine to determine improved Lambda's

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IMPLICIT INTEGER (I-N)

        COMMON/GAUGE/ALSMZ,ALEMMZ,GF,gg1,gg2,S2TW
        COMMON/SMSPEC/AMS,AMC,AMB,AMBP,AMT,AMTAU,AMMUON,AMZ,AMW

        B0(NF)=33.D0-2.D0*NF
        B1(NF)=6.D0*(153.D0-19.D0*NF)/B0(NF)**2
        ALS2(NF,X,XLB)=12.D0*PI/(B0(NF)*DLOG(X**2/XLB**2))
     .            *(1.D0-B1(NF)*DLOG(DLOG(X**2/XLB**2))
     .            /DLOG(X**2/XLB**2))
        AA(NF)=12D0*PI/B0(NF)
        BB(NF)=B1(NF)/AA(NF)
        XIT(A,B,X)=A/2.D0*(1D0+DSQRT(1D0-4D0*B*DLOG(X)))
        PI=4.D0*DATAN(1.D0)
        NF=5
        Q=AMZ
        XLB=Q*DEXP(-AA(NF)/ALP/2.D0)
        IF(NO.EQ.1)GOTO 111
        II=0
1       II=II+1
        X=DLOG(Q**2/XLB**2)
        A=AA(NF)/ALP
        B=BB(NF)*ALP
        XX=XIT(A,B,X)
        XLB=Q*DEXP(-XX/2.D0)
        Y1=ALP
        Y2=ALS2(NF,Q,XLB)
        DY=DABS(Y2-Y1)/Y1
        IF(DY.GE.ACC) GOTO 1
111     XITLA=XLB

        RETURN
        END


        DOUBLE PRECISION FUNCTION FINT(Z,XX,YY)

*  One-dimensional cubic interpolation
*  Z  = wanted point
*  XX = array of 4 discrete x-values around Z
*  YY = array of 4 discrete function-values around Z

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

        DIMENSION XX(4),YY(4)

        X = DLOG(Z)
        X0=DLOG(XX(1))
        X1=DLOG(XX(2))
        X2=DLOG(XX(3))
        X3=DLOG(XX(4))
        Y0=DLOG(YY(1))
        Y1=DLOG(YY(2))
        Y2=DLOG(YY(3))
        Y3=DLOG(YY(4))
        A0=(X-X1)*(X-X2)*(X-X3)/(X0-X1)/(X0-X2)/(X0-X3)
        A1=(X-X0)*(X-X2)*(X-X3)/(X1-X0)/(X1-X2)/(X1-X3)
        A2=(X-X0)*(X-X1)*(X-X3)/(X2-X0)/(X2-X1)/(X2-X3)
        A3=(X-X0)*(X-X1)*(X-X2)/(X3-X0)/(X3-X1)/(X3-X2)
        FINT=DEXP(A0*Y0+A1*Y1+A2*Y2+A3*Y3)

        RETURN
        END


*   Spence function and auxiliary functions as in HDECAY

        DOUBLE PRECISION FUNCTION SP(X)

*  REAL dilogarithm (Spence-function)

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DOUBLE COMPLEX CX,LI2

        CX = DCMPLX(X,0.D0)
        SP = DREAL(LI2(CX))

        RETURN
        END

        DOUBLE COMPLEX FUNCTION LI2(X)

*  COMPLEX dilogarithm (Spence-function)

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IMPLICIT INTEGER (I-N)
        DOUBLE COMPLEX X,Y,CLI2

        COMMON/CONST/ZETA2,ZETA3

        ZERO=1.D-16
        XR=DREAL(X)
        XI=DIMAG(X)
        R2=XR*XR+XI*XI
        LI2=0
        IF(R2.LE.ZERO)THEN
          LI2=X
          RETURN
        ENDIF
        RR=XR/R2
        IF(R2.EQ.1.D0.AND.XI.EQ.0.D0)THEN
          IF(XR.EQ.1.D0)THEN
            LI2=DCMPLX(ZETA2)
          ELSE
            LI2=-DCMPLX(ZETA2/2.D0)
          ENDIF
          RETURN
        ELSEIF(R2.GT.1.D0.AND.RR.GT.0.5D0)THEN
          Y=(X-1.D0)/X
          LI2=CLI2(Y)+ZETA2-CDLOG(X)*CDLOG(1.D0-X)+0.5D0*CDLOG(X)**2
          RETURN
        ELSEIF(R2.GT.1.D0.AND.RR.LE.0.5D0)THEN
          Y=1.D0/X
          LI2=-CLI2(Y)-ZETA2-0.5D0*CDLOG(-X)**2
          RETURN
        ELSEIF(R2.LE.1.D0.AND.XR.GT.0.5D0)THEN
          Y=1.D0-X
          LI2=-CLI2(Y)+ZETA2-CDLOG(X)*CDLOG(1.D0-X)
         RETURN
        ELSEIF(R2.LE.1.D0.AND.XR.LE.0.5D0)THEN
          Y=X
          LI2=CLI2(Y)
          RETURN

        ENDIF
        END
 

        DOUBLE COMPLEX FUNCTION CLI2(X)

*  Taylor-expansion for complex dilogarithm (Spence-function)

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IMPLICIT INTEGER (I-N)
        DOUBLE COMPLEX X,Z
        
        COMMON/BERNOULLI/B2(18),B12(18),B3(18)
        COMMON/POLY/NBER
        
        N=NBER-1
        Z=-CDLOG(1.D0-X)
        CLI2=B2(NBER)
        DO 111 I=N,1,-1
          CLI2=Z*CLI2+B2(I)
111     CONTINUE
        CLI2=Z**2*CLI2+Z

        RETURN
        END

 
        DOUBLE PRECISION FUNCTION FACULT(N)

*  DOUBLE PRECISION version of FACULTY

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IMPLICIT INTEGER (I-N)

        FACULT=1.D0
        IF(N.EQ.0)RETURN
        DO 999 I=1,N
          FACULT=FACULT*DFLOAT(I)
999     CONTINUE

        RETURN
        END

 
        SUBROUTINE BERNINI(N)

*  Initialization of coefficients for polylogarithms

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IMPLICIT INTEGER (I-N)
        DIMENSION B(18),PB(19)

        COMMON/BERNOULLI/B2(18),B12(18),B3(18)
        COMMON/CONST/ZETA2,ZETA3
        COMMON/POLY/NBER

        NBER=N
        PI=4.D0*DATAN(1.D0)

        B(1)=-1.D0/2.D0
        B(2)=1.D0/6.D0
        B(3)=0.D0
        B(4)=-1.D0/30.D0
        B(5)=0.D0
        B(6)=1.D0/42.D0
        B(7)=0.D0
        B(8)=-1.D0/30.D0
        B(9)=0.D0
        B(10)=5.D0/66.D0
        B(11)=0.D0
        B(12)=-691.D0/2730.D0
        B(13)=0.D0
        B(14)=7.D0/6.D0
        B(15)=0.D0
        B(16)=-3617.D0/510.D0
        B(17)=0.D0
        B(18)=43867.D0/798.D0
        ZETA2=PI**2/6.D0
        ZETA3=1.202056903159594D0

        DO 995 I=1,18
          B2(I)=B(I)/FACULT(I+1)
          B12(I)=DFLOAT(I+1)/FACULT(I+2)*B(I)/2.D0
          PB(I+1)=B(I)
          B3(I)=0.D0
995     CONTINUE
        PB(1)=1.D0
        DO 996 I=1,18
        DO 996 J=0,I
         B3(I)=B3(I)+PB(J+1)*PB(I-J+1)/FACULT(I-J)/FACULT(J+1)
     .         /DFLOAT(I+1)
996     CONTINUE

        RETURN
        END


c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c   Passarino-Veltman one- and two-points functions A0, B0 and B1  
c   orig from LoopTools, http://www.feynarts.de/looptools/
c   taken from Suspect2.3, modified by S. Kraml, 7 March 2005
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        double precision function NMA0(m2,q)

        implicit none
        double precision m2,q
        if(m2.ne.0.d0) then
         NMA0 = m2 * (1.d0-dlog( m2/q )) 
        else
         NMA0 = 0.d0
        endif
        end

        double precision function NMB0(p,m1,m2,q)

*     note: all input is quadratical: p=p^2, m1=m1^2, m2=m2^2, q=q^2

        implicit none
        double precision p, m1, m2
        double precision mudim2, divergence, lambda2, q
        double precision acc, eps, minacc
        double complex x1, x2, y1, y2, r, be0
        double complex Ieps, onePeps, oneMeps
        common/cutoff/mudim2, divergence, lambda2
        parameter (acc = 1.D-12)
        parameter (eps = 1.D-20)
        parameter (Ieps = (0.D0,1.D0)*eps)
        parameter (onePeps = 1.D0 + Ieps)
        parameter (oneMeps = 1.D0 - Ieps)

        double complex fpv, xlogx
        external fpv, xlogx

        divergence = 0.d0
        lambda2 = 0.d0
        mudim2 = q
        minacc = acc*(m1 + m2)

* general case
        if(abs(p) .gt. minacc) then
          call roots(p, m1, m2, x1, x2, y1, y2, r)
          if(abs(y1) .gt. .5D0 .and. abs(y2) .gt. .5D0) then
            be0 = -log(m2/mudim2) - 
     +        fpv(1, x1, y1) - fpv(1, x2, y2)
          else if(abs(x1) .lt. 10.D0 .and. abs(x2) .lt. 10.D0) then
            be0 = 2 - log(p*oneMeps/mudim2) +
     +        xlogx(-x1) + xlogx(-x2) - xlogx(y1) - xlogx(y2)
          else if(abs(x1) .gt. .5D0 .and. abs(x2) .gt. .5D0) then
            be0 = -log(m1/mudim2) -
     +        fpv(1, y1, x1) - fpv(1, y2, x2)
          else
            be0 = 1.D100
          endif

* zero momentum
        else if(abs(m1 - m2) .gt. minacc) then
          x2 = oneMeps*m1/(m1 - m2)
          y2 = oneMeps*m2/(m2 - m1)
          if(abs(y2) .gt. .5D0) then
            be0 = -log(m2/mudim2) - fpv(1, x2, y2)
          else
            be0 = -log(m1/mudim2) - fpv(1, y2, x2)
          endif
        else
          be0 = -log(m2/mudim2)
        endif

        NMB0 = dble(be0 + divergence)

        end


*---------------------------------------------------------------------
* auxiliary functions used by the B0,B1 two-point functions
* from Looptools http://www.feynarts.de/looptools/
*---------------------------------------------------------------------

        subroutine roots(p, m1, m2, x1, x2, y1, y2, r)

        implicit none
        double precision p, m1, m2
        double complex x1, x2, y1, y2, r
        double precision mudim2, divergence, lambda2
        common/cutoff/mudim2, divergence, lambda2
        double precision acc, eps
        double complex Ieps, onePeps, oneMeps
        parameter (acc = 1D-12)
        parameter (eps = 1D-20)
        parameter (Ieps = (0.D0,1.D0)*eps)
        parameter (onePeps = 1.D0 + Ieps)
        parameter (oneMeps = 1.D0 - Ieps)
        double precision q

        r = sqrt(dcmplx(p*(p - 2*(m1 + m2)) + (m1 - m2)**2))
        q = p + m1 - m2
        x1 = (q + r)/2.D0/p
        x2 = (q - r)/2.D0/p
        if(abs(x2) .gt. abs(x1)) then
          x1 = m1/p/x2
        else if(abs(x1) .gt. abs(x2)) then
          x2 = m1/p/x1
        endif
        x1 = x1 + abs(p*x1)/p*Ieps
        x2 = x2 - abs(p*x2)/p*Ieps
        q = p - m1 + m2
        y2 = (q + r)/2.D0/p
        y1 = (q - r)/2.D0/p
        if(abs(y2) .gt. abs(y1)) then
          y1 = m2/p/y2
        else if(abs(y1) .gt. abs(y2)) then
          y2 = m2/p/y1
        endif
        y1 = y1 - abs(p*y1)/p*Ieps
        y2 = y2 + abs(p*y2)/p*Ieps
        end


        double complex function fpv(n, x, y)

        implicit none
        integer n
        double complex x, y
        double precision mudim2, divergence, lambda2
        common/cutoff/mudim2, divergence, lambda2
        double precision acc, eps
        double complex Ieps, onePeps, oneMeps
        parameter (acc = 1D-12)
        parameter (eps = 1D-20)
        parameter (Ieps = (0,1)*eps)
        parameter (onePeps = 1 + Ieps)
        parameter (oneMeps = 1 - Ieps)
        integer m
        double complex xm
        if(abs(x) .lt. 10.D0) then
          if(n .eq. 0) then
            fpv = -log(-y/x)
          else if(abs(x) .lt. acc) then
            fpv = -1.D0/n
          else
            fpv = 0
            xm = 1
            do m = 0, n - 1
              fpv = fpv - xm/(n - m)
              xm = xm*x
            enddo
            fpv = fpv - xm*log(-y/x)
          endif
        else
          fpv = 0
          xm = 1
          do m = 1, 30
            xm = xm/x
            fpv = fpv + xm/(m + n)
            if(abs(xm/fpv) .lt. acc**2) return
          enddo
        endif
        end


        double complex function yfpv(n, x, y)

        implicit none
        integer n
        double complex x, y
        double complex fpv
        external fpv
        if(abs(y) .eq. 0.D0) then
          yfpv = 0
        else
          yfpv = y*fpv(n, x, y)
        endif
        end


        double complex function xlogx(x)

        implicit none
        double complex x
        if(abs(x) .eq. 0.D0) then
          xlogx = 0
        else
          xlogx = x*log(x)
        endif
        end


        DOUBLE PRECISION FUNCTION RUNMB(Q)

*   Subroutine to calculate the running b quark mass for Q > MB

        IMPLICIT NONE
        DOUBLE PRECISION Q
        DOUBLE PRECISION PI,ALPHAS,ALMB,ALMT,ALQ,U5MTMB,U6QMT
        DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW

        COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW

        PI=4.D0*DATAN(1.D0)

        ALQ=ALPHAS(Q,2)
        ALMB=ALPHAS(MB,2)

        IF(Q.LE.MT) THEN

         RUNMB=MB*(ALQ/ALMB)**(12.D0/23.D0)*(1.D0+7462.D0*(ALQ-ALMB)/
     .         (4.D0*PI*1587.D0))

        ELSE

         ALMT=ALPHAS(MT,2)
         U5MTMB=(ALMT/ALMB)**(12.D0/23.D0)*(1.D0+7462.D0*(ALMT-ALMB)/
     .         (4.D0*PI*1587.D0))
         U6QMT=(ALQ/ALMT)**(4.D0/7.D0)*(1.D0+7398.D0*(ALQ-ALMT)/
     .         (4.D0*PI*1323.D0))
         RUNMB=MB*U6QMT*U5MTMB

        ENDIF

        END
