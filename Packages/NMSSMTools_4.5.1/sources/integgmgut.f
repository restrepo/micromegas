*   Subroutines to integrate the RGEs for the gauge and Yukawa couplings
*   above the messenger scale: the 1-loop coefficients of the gauge beta
*   functions include the contributions from n5 pairs of messengers in
*   5- and 5_bar-representations of SU(5)

      SUBROUTINE
     .   ODEINTGMGUT(YSTART,NVAR,X1,X2,EPS,DERIVSGMGUT,RKQSGMGUT,IFAIL)

*   Driver subroutine to integrate Ordinary Differential Equations
*   using the Runge-Kutta 5th order method with adaptative step size
*   (from Numerical Recipes)
*   IFAIL=1 stepsize smaller than minimum in ODEINTGMGUT
*   IFAIL=2 too many steps in ODEINTGMGUT
*   IFAIL=3 stepsize underflow in RKQSGMGUT
*   IFAIL=4 max(|dy(i)/dx|*(1+x)/(1+|y(i)|)) > 1/eps^2

      IMPLICIT NONE

      INTEGER NVAR,NMAX,IFAIL,MAXSTP,NSTP,I
      PARAMETER(MAXSTP=10000,NMAX=500)

      DOUBLE PRECISION YSTART(NVAR),Y(NMAX),YSCAL(NMAX)
      DOUBLE PRECISION DYDX(NMAX),X,X1,X2,EPS,TINY
      DOUBLE PRECISION H1,HMIN,H,HDID,HNEXT

      EXTERNAL DERIVSGMGUT,RKQSGMGUT

      IFAIL=0
      TINY=EPS**4
      H1=DSQRT(EPS)*DABS(X2-X1)
      HMIN=EPS**2*DABS(X2-X1)

      X=X1
      H=DSIGN(H1,X2-X1)
      DO I=1,NVAR
       Y(I)=YSTART(I)
      ENDDO

      DO NSTP=1,MAXSTP

       CALL DERIVSGMGUT(NVAR,X,Y,DYDX)

       DO I=1,NVAR
        IF(DABS(DYDX(I))*(1d0+X)/(1d0+DABS(Y(I))).GT.EPS**(-2))IFAIL=4
       ENDDO
       IF(IFAIL.GT.0)RETURN

       DO I=1,NVAR
        YSCAL(I)=DABS(Y(I))+DABS(H*DYDX(I))+TINY
       ENDDO

       CALL RKQSGMGUT(Y,DYDX,NVAR,X,X2,H,EPS,YSCAL,
     .  HDID,HNEXT,DERIVSGMGUT,IFAIL)
       IF(IFAIL.GT.0)RETURN
      
       IF(DABS(Y(2)-5d0/3d0*Y(1)).LT.EPS)THEN
        X2=X
        DO I=1,NVAR
         YSTART(I)=Y(I)
        ENDDO
        RETURN
       ENDIF

       IF(ABS(HNEXT).LT.HMIN)THEN
        X2=X
        DO I=1,NVAR
         YSTART(I)=Y(I)
        ENDDO
        IFAIL=1
        RETURN
       ENDIF

       H=HNEXT

      ENDDO

      IFAIL=2

      RETURN
      END


      SUBROUTINE RKQSGMGUT(Y,DYDX,N,X,X2,HTRY,EPS,YSCAL,HDID,
     .    HNEXT,DERIVSGMGUT,IFAIL)

*   Stepper subroutine for ODEINTGMGUT

      IMPLICIT NONE

      INTEGER IFAIL,I,N,NMAX
      PARAMETER(NMAX=500)

      DOUBLE PRECISION EPS,HDID,HNEXT,HTRY,X,X2,DYDX(N),Y(N),YSCAL(N)
      DOUBLE PRECISION ERRMAX,H,HTEMP,XNEW,YERR(NMAX),YTEMP(NMAX)
      DOUBLE PRECISION SAFETY,PGROW,PSHRNK,ERRCON,G,GTEMP

      EXTERNAL DERIVSGMGUT

      SAFETY=.9d0
      PGROW=-.2d0
      PSHRNK=-.25d0
      ERRCON=(5d0/SAFETY)**(1d0/PGROW)
      H=HTRY
      G=Y(2)-5d0/3d0*Y(1)

1     CALL RKCKGMGUT(Y,DYDX,N,X,H,YTEMP,YERR,DERIVSGMGUT)

      ERRMAX=0d0
      DO I=1,N
       ERRMAX=MAX(ERRMAX,DABS(YERR(I)/YSCAL(I)))
      ENDDO
      ERRMAX=ERRMAX/EPS

      IF(ERRMAX.GT.1d0)THEN
       HTEMP=SAFETY*H*(ERRMAX**PSHRNK)
       H=SIGN(MAX(DABS(HTEMP),.1d0*DABS(H)),H)
       XNEW=X+H
       IF(XNEW.EQ.X)THEN
        IFAIL=3
        RETURN
       ENDIF
       GOTO 1
      ENDIF

      GTEMP=YTEMP(2)-5d0/3d0*YTEMP(1)
      IF(GTEMP.LT.-EPS)THEN
       IFAIL=-1
       X2=X+H
       H=G/(G-GTEMP)*H
       GOTO 1
      ENDIF

      HDID=H
      X=X+H
      DO I=1,N
       Y(I)=YTEMP(I)
      ENDDO

      IF(ERRMAX.GT.ERRCON)THEN
       HNEXT=SAFETY*H*(ERRMAX**PGROW)
      ELSE
        HNEXT=5d0*H
      ENDIF
      
      IF(X+HNEXT.GT.X2 .AND. IFAIL.EQ.-1) HNEXT=X2-X

      RETURN

      END


      SUBROUTINE RKCKGMGUT(Y,DYDX,N,X,H,YOUT,YERR,DERIVSGMGUT)

*   Algorithm subroutine for ODEINTGMGUT

      IMPLICIT NONE

      INTEGER I,N,NMAX
      PARAMETER(NMAX=500)

      DOUBLE PRECISION H,X,DYDX(N),Y(N),YERR(N),YOUT(N)
      DOUBLE PRECISION AK(NMAX),AK3(NMAX),AK4(NMAX)
      DOUBLE PRECISION AK5(NMAX),AK6(NMAX),YTEMP(NMAX)
      DOUBLE PRECISION A2,A3,A4,A5,A6,C1,C3,C4,C6
      DOUBLE PRECISION B21,B31,B32,B41,B42,B43,B51,B52,B53,B54
      DOUBLE PRECISION B61,B62,B63,B64,B65,DC1,DC3,DC4,DC5,DC6

      PARAMETER(A2=.2d0, A3=.3d0, A4=.6d0, A5=1d0, A6=.875d0,
     .       B21=.2d0, B31=3d0/40d0, B32=9d0/40d0, B41=.3d0, B42=-.9d0,
     .       B43=1.2d0, B51=-11d0/54d0, B52=2.5d0, B53=-70d0/27d0,
     .       B54=35d0/27d0, B61=1631d0/55296d0, B62=175d0/512d0,
     .       B63=575d0/13824d0, B64=44275d0/110592d0,
     .       B65=253d0/4096d0, C1=37d0/378d0, C3=250d0/621d0,
     .       C4=125d0/594d0, C6=512d0/1771d0, DC1=C1-2825d0/27648d0,
     .       DC3=C3-18575d0/48384d0, DC4=C4-13525d0/55296d0,
     .       DC5=-277d0/14336d0, DC6=C6-.25d0)

      EXTERNAL DERIVSGMGUT

      DO I=1,N
       YTEMP(I)=Y(I)+B21*H*DYDX(I)
      ENDDO
      CALL DERIVSGMGUT(N,X+A2*H,YTEMP,AK)
      DO I=1,N
       YTEMP(I)=Y(I)+H*(B31*DYDX(I)+B32*AK(I))
      ENDDO
      CALL DERIVSGMGUT(N,X+A3*H,YTEMP,AK3)
      DO I=1,N
       YTEMP(I)=Y(I)+H*(B41*DYDX(I)+B42*AK(I)+B43*AK3(I))
      ENDDO
      CALL DERIVSGMGUT(N,X+A4*H,YTEMP,AK4)
      DO I=1,N
       YTEMP(I)=Y(I)+H*(B51*DYDX(I)+B52*AK(I)+B53*AK3(I)+B54*AK4(I))
      ENDDO
      CALL DERIVSGMGUT(N,X+A5*H,YTEMP,AK5)
      DO I=1,N
       YTEMP(I)=Y(I)+H*(B61*DYDX(I)+B62*AK(I)+B63*AK3(I)+B64*AK4(I)+
     .  B65*AK5(I))
      ENDDO
      CALL DERIVSGMGUT(N,X+A6*H,YTEMP,AK6)
      DO I=1,N
       YOUT(I)=Y(I)+H*(C1*DYDX(I)+C3*AK3(I)+C4*AK4(I)+C6*AK6(I))
      ENDDO
      DO I=1,N
       YERR(I)=H*(DC1*DYDX(I)+DC3*AK3(I)+DC4*AK4(I)+DC5*AK5(I)+DC6*
     .  AK6(I))
      ENDDO

      RETURN
      END


      SUBROUTINE DERIVSGMGUT(N,X,Y,F)

*   2-loop Renormalization group equations for G1, G2, G3,
*   lambda, kappa, htop, hbot, htau to be integrated by ODEINT

      IMPLICIT NONE

      INTEGER N
      DOUBLE PRECISION X,Y(N),F(N),PI,c2
      DOUBLE PRECISION G1,G2,G3,L,K,HT,HB,HL
      DOUBLE PRECISION LPP,LTT,LU,LD,LT,LB,LL
      DOUBLE PRECISION MSUSYEFF,MMESS,N5

      COMMON/MESCAL/MSUSYEFF,MMESS,N5

      PI=4d0*DATAN(1d0)
      c2=0d0 ! 1d0/(16d0*PI**2)

      G1=Y(1)
      G2=Y(2)
      G3=Y(3)
      L=Y(4)
      K=Y(5)
      HT=Y(6)
      HB=Y(7)
      HL=Y(8)
      LPP=Y(9)
      LTT=Y(10)
      LU=Y(11)
      LD=Y(12)
      LT=Y(13)
      LB=Y(14)
      LL=Y(15)

      F(1)= (N5*5d0/3d0+11d0)*G1**2
     .    + c2*G1**2*(199d0/9d0*G1 + 9d0*G2 + 88d0/3d0*G3
     .    - 2d0*L - 26d0/3d0*HT - 14d0/3d0*HB - 6d0*HL)

      F(2)= (N5+1d0)*G2**2
     .    + c2*G2**2*(3d0*G1 + 25d0*G2 + 24d0*G3
     .    - 2d0*L - 6d0*HT - 6d0*HB - 2d0*HL)

      F(3)= (N5-3d0)*G3**2
     .    + c2*G3**2*(11d0/3d0*G1 + 9d0*G2 + 14d0*G3
     .    - 4 d0*HT - 4d0*HB)

      F(4)= L*(-G1 - 3d0*G2 + 4d0*L + 2d0*K + 3d0*HT
     .    + 3d0*HB + HL + 2d0*LPP + 3d0*LTT + 4d0*LU + 4d0*LD)
     .    + 3d0*DSQRT(L*LD*LT*HT) + 3d0*DSQRT(L*LU*LB*HB)
     .    + DSQRT(L*LU*LL*HL)
     .    + c2*L*(L*(2d0*G1 + 6d0*G2
     .    - 10d0*L - 12d0*K - 9d0*HT - 9d0*HB - 3d0*HL)
     .    + 23d0/2d0*G1**2 + 15d0/2d0*G2**2 + 3d0*G1*G2
     .    + 4d0/3d0*G1*HT - 2d0/3d0*G1*HB + 2d0*G1*HL
     .    + 16d0*G3*HT + 16d0*G3*HB - 3d0*HL**2
     .    - 9d0*HT**2 - 9d0*HB**2 - 6d0*HT*HB - 8d0*K**2)

      F(5)= 6d0*K*(L + K + LPP + 1.5d0*LTT + LU + LD)
     .    + c2*6d0*K*(-4d0*K**2 + L*(G1 + 3d0*G2
     .    - 2d0*L - 4d0*K - 3d0*HT - 3d0*HB - HL))

      F(6)= HT*(-13d0/9d0*G1 - 3d0*G2 - 16d0/3d0*G3
     .    + L + 6d0*HT + HB + LU + 6d0*LT + LB)
     .    + DSQRT(L*LD*LT*HT)
     .    + c2*HT*(HT*(2d0*G1 + 6d0*G2 + 16d0*G3
     .    - 3d0*L - 22d0*HT - 5d0*HB)
     .    + 2d0/3d0*G1*HB - 3d0*L**2 - 2d0*L*K
     .    - 4d0*L*HB - 5d0*HB**2 - HB*HL - L*HL
     .    + 2743d0/162d0*G1**2 + 15d0/2d0*G2**2- 16d0/9d0*G3**2
     .    + 5d0/3d0*G1*G2 + 136d0/27d0*G1*G3 + 8d0*G2*G3)

      F(7)= HB*(-7d0/9d0*G1 - 3d0*G2 - 16d0/3d0*G3
     .    + L + HT + 6d0*HB + HL + LD + 6d0*LB + LT)
     .    + DSQRT(L*LU*LB*HB) + DSQRT(LB*LL*HB*HL)
     .    + c2*HB*(HB*(2d0/3d0*G1 + 6d0*G2 + 16d0*G3
     .    - 3d0*L - 5d0*HT - 22d0*HB - 3d0*HL)
     .    + 4d0/3d0*G1*HT - 3d0*L**2 - 2d0*L*K
     .    - 4d0*L*HT - 5d0*HT**2 - 3d0*HL**2 + 2d0*G1*HL
     .    + 1435d0/162d0*G1**2 + 15d0/2d0*G2**2 - 16d0/9d0*G3**2
     .    + 5d0/3d0*G1*G2 + 40d0/27d0*G1*G3 + 8d0*G2*G3)

      F(8)= HL*(-3d0*G1 - 3d0*G2
     .    + L + 3d0*HB + 4d0*HL + LD + 4*LL)
     .    + DSQRT(L*LU*LL*HL) + 3d0*DSQRT(LB*LL*HB*HL)
     .    + c2*HL*(-10d0*HL**2 - 9d0*HL*HB - 9d0*HB**2
     .    - 3d0*HB*HT + HL*(2d0*G1 + 6d0*G2)
     .    + HB*(-2d0/3d0*G1 + 16*G3)
     .    - L*(3d0*HL + 3d0*HT + 3d0*L + 2d0*K)
     .    + 75d0/2d0*G1**2 + 3d0*G1*G2 + 15d0/2d0*G2**2)

      F(9)= LPP*(-G1 - 3d0*G2 + 2d0*L + 2d0*K
     .    + 4d0*LPP + 3d0*LTT + 2d0*LU + 2d0*LD)

      F(10)= LTT*(-4d0/9d0*G1 - 16d0/3d0*G3 + 2d0*L + 2d0*K
     .     + 2d0*LPP + 5d0*LTT + 2d0*LU + 2d0*LD)

      F(11)= LU*(-G1 - 3d0*G2 + 4d0*L + 2d0*K
     .    + 2d0*LPP + 3d0*LTT + 4d0*LU + 2d0*LD
     .    + 3d0*HT + 3d0*LB + LL
     .    + 3d0*DSQRT(L*LU*LB*HB) + DSQRT(L*LU*LL*HL))

      F(12)= LD*(-G1 - 3d0*G2 + 4d0*L + 2d0*K
     .    + 2d0*LPP + 3d0*LTT + 2d0*LU + 4d0*LD
     .    + 3d0*HB + HL + 3d0*LT
     .    + 3d0*DSQRT(L*LD*LT*HT))

      F(13)= LT*(-13d0/9d0*G1 - 3d0*G2 - 16d0/3d0*G3
     .     + 6d0*LT + 6d0*HT + HB + LD)
     .     + DSQRT(L*LD*LT*HT)

      F(14)= LB*(-7d0/9d0*G1 - 3d0*G2 - 16d0/3d0*G3
     .     + 6d0*LB + LL + 6d0*HB + HT + LU)
     .     + DSQRT(L*LU*LB*HB)
     .     + DSQRT(LB*HB*LL*HL)

      F(15)= LL*(-3d0*G1 - 3d0*G2
     .     + 4d0*LL + 3d0*LB + 4d0*HL + LU)
     .     + DSQRT(L*LU*LL*HL)
     .     + 3d0*DSQRT(LL*HL*LB*HB)

      RETURN
      END
