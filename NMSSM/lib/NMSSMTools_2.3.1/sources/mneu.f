	SUBROUTINE NEUTRALINO(PAR)

* Subroutine to compute the neutralino masses MNEU(i) (i=1..5, 
* ordered in mass) and the neutralino mixing matrix NEU(i,j)
* 1 loop rad. corrs. are taken into account (assuming degenerate 
* squarks/sleptons) as in BMPZ (Bagger et al., hep-ph/9606211).
* It is assumed that the input parameters M1, M2 and mu
* are the running DRbar masses at the scale Q2.
* The auxiliary functions B0 and B1 are defined in aux.f

	IMPLICIT NONE

	INTEGER I,J
	DOUBLE PRECISION PAR(*),PI
	DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
	DOUBLE PRECISION VALP(5),VECP(5,5),EPS
	DOUBLE PRECISION l,k,TanB,mu,M1,M2,v1,v2
	DOUBLE PRECISION ALSMZ,ALP0,GF,g1,g2,S2TW,Q2
	DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
	DOUBLE PRECISION MSQ,MSL,SIN2B,M22,MU2,MZ2,MW2,M12
	DOUBLE PRECISION QSTSB,B0,B1
	DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
	DOUBLE PRECISION HTQ,HBQ,MTQ,MBQ,MA2
	DOUBLE PRECISION XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY

	COMMON/GAUGE/ALSMZ,ALP0,GF,g1,g2,S2TW
	COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
	COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
	COMMON/RENSCALE/Q2
	COMMON/STSBSCALE/QSTSB
	COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
	COMMON/QQUARK/HTQ,HBQ,MTQ,MBQ
	COMMON/GMSUSYPAR/XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY

	EPS=1.D-8
	PI=4.D0*DATAN(1.D0)

	l=par(1)
	k=par(2)
	TanB=PAR(3)
	M1=PAR(20)
	M2=PAR(21)
	MA2=PAR(23)**2
	v2=1.D0/DSQRT(2.D0*DSQRT(2.D0)*(1.D0+TanB**2)*GF)
	v1=v2*TanB
	SIN2B=2.D0*TANB/(1.D0+TANB**2)
	M22=M2**2
	M12=M1**2
	MZ2=MZ**2
	MW2=MW**2

* For 1 loop rad. corrs.:
* Average squark mass squared:
	MSQ=(2.D0*PAR(7)+PAR(8)+PAR(9)
     C       +4.D0*PAR(15)+2.D0*PAR(16)+2.D0*PAR(17))/12.D0

* Average slepton mass squared:
	MSL=(2.D0*PAR(10)+PAR(11)+4.D0*PAR(18)
     C       +2.D0*PAR(19))/9.D0
 
* First: MU at the scale Q2:

	MU=PAR(4)
	MU2=MU**2

* MU incl. rad. corrs.:

	MU=MU*(1.D0-1.D0/(64.D0*PI**2)*(
     C         12.D0*(HTQ**2+HBQ**2)*B1(MU2,0.D0,QSTSB,Q2)
     C         +(3.D0*G2+G1)*(B1(MU2,M22,MA2,Q2)+B1(MU2,M22,MZ2,Q2)
     C              +2.D0*B1(MU2,MU2,MZ2,Q2)-4.D0*B0(MU2,MU2,MZ2,Q2))))

* M1 incl. rad. corrs.:

	M1=M1*(1.D0-G1/(16.D0*PI**2)*(11.D0*B1(M12,0.D0,MSQ,Q2) 
     C                            +9.D0* B1(M12,0.D0,MSL,Q2)
     C     +SIN2B*MU/M1*(B0(M12,MU2,MA2,Q2)-B0(M12,MU2,MZ2,Q2))
     C     +B1(M12,MU2,MA2,Q2)+B1(M12,MU2,MZ2,Q2)))
  
* M2 incl. rad. corrs.:

	M2=M2*(1.D0-G2/(16.D0*PI**2)*(9.D0*B1(M22,0.D0,MSQ,Q2) 
     C                            +3.D0* B1(M22,0.D0,MSL,Q2)
     C     +SIN2B*MU/M2*(B0(M22,MU2,MA2,Q2)-B0(M22,MU2,MZ2,Q2))
     C     +B1(M22,MU2,MA2,Q2)+B1(M22,MU2,MZ2,Q2)
     C     -8.D0*B0(M22,M22,MW2,Q2)+4.D0*B1(M22,M22,MW2,Q2)))

        NEU(1,1)=M1
        NEU(1,2)=0.D0
        NEU(1,3)=DSQRT(g1/2.D0)*v1
        NEU(1,4)=-DSQRT(g1/2.D0)*v2
        NEU(1,5)=0.D0
        NEU(2,2)=M2
        NEU(2,3)=-DSQRT(g2/2.D0)*v1
        NEU(2,4)=DSQRT(g2/2.D0)*v2
        NEU(2,5)=0.D0
        NEU(3,3)=0.D0
        NEU(3,4)=-mu
        NEU(3,5)=-l*v2
        NEU(4,4)=0.D0
        NEU(4,5)=-l*v1
        NEU(5,5)=2.D0*k/l*mu+MUPSUSY

	CALL DIAGN(5,NEU,VALP,VECP,EPS)
	CALL SORTNA(5,VALP,VECP)
	DO I=1,5
	 MNEU(I)=VALP(I)
	 DO J=1,5
	  NEU(I,J)=VECP(J,I)
	 ENDDO
	ENDDO

	END


