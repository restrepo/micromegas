C              Computation of the Muon Anomalous Moment
C
C	 - Literature:
C	[1]:  F. Jegerlehner, "Essentials of the Muon g-2.", 
C             Acta Phys.Polon.B38:3021,2007, hep-ph/0703125
C	[2]:  J. Bijnens and J. Prades, "The hadronic light-by-light
C	      contribution to the muon anomalous magnetic moment: Where
C	      do we stand?", Mod. Phys. Lett. A 22 (2007) 767,
C	      arXiv:hep-ph/0702170
C	[3]:  A. Czarnecki, W. J. Marciano, A. Vainshtein," Refinements in 
C             electroweak contributions to the muon anomalous magnetic moment."
C             Phys.Rev.D67:073006,2003, Erratum-ibid.D73:119901,2006, hep-ph/0212229
C	[4]:  S. P. Martin, J. D. Wells
C             "Muon anomalous magnetic dipole moment in supersymmetric theories."
C             Phys.Rev.D64:035003,2001, hep-ph/0103067
C	[5]:  J. P. Leveille, "The Second Order Weak Correction to (G-2) of the 
C             Muon in Arbitrary Gauge Models."
C             Nucl.Phys.B137:63,1978
C	[6]:  J. F. Gunion, D. Hooper, B. McElrath
C             "Light neutralino dark matter in the NMSSM"
C             Phys.Rev.D73:015011,2006, hep-ph/0509024
C	[7]:  G. Degrassi, G.F. Giudice, "QED logarithms in the electroweak
C             corrections to the muon anomalous magnetic moment."
C             Phys.Rev.D58:053007,1998, hep-ph/9803384
C	[8]:  S. Heinemeyer, D. Stockinger, G. Weiglein,
C             "Electroweak and supersymmetric two-loop corrections to (g-2)(mu)."
C             Nucl.Phys.B699:103-123,2004, hep-ph/0405255
C	[9]:  K. Cheung, C.-H. Chou, O.C.W. Kong, "Muon anomalous magnetic 
C             moment, two Higgs doublet model, and supersymmetry."
C             Phys.Rev.D64:111301,2001, hep-ph/0103183
C	[10]: A. Arhrib, S. Baek, "Two loop Barr-Zee type 
C             contributions to (g-2)(muon) in the MSSM."
C             Phys.Rev.D65:075002,2002, hep-ph/0104225
C	[11]: D. Stockinger, "The Muon Magnetic Moment and Supersymmetry."
C             J.Phys.G34:R45-R92,2007, hep-ph/0609168


        SUBROUTINE MAGNMU(PAR,PROB)
	IMPLICIT NONE

	INTEGER I,K
	DOUBLE PRECISION PAR(*),PROB(*)

	DOUBLE PRECISION MSMU(2),USMU(2,2),UST(2,2),USB(2,2),USL(2,2),
     C                   MST(2),MSB(2)
	DOUBLE PRECISION aux,Ymu,Pi,tanb,cosb,sinb,ALEM0,QNP
	DOUBLE PRECISION fc1,fc2,fn1,fn2,fhs,fhp,fhc,runmb,asf
	DOUBLE PRECISION At,Ab,Atau,lambda,Yl,fS,fPS,fSF,mtq,mbq
	DOUBLE PRECISION MSTAU(2),lambdHT(3,2),lambdHB(3,2),lambdHL(3,2)
	DOUBLE PRECISION amu2Lbos,amu2LHf,amu2LCh,amu2LSF,amu2L,amuerr,
     C                   amuChar,amuNeutr,amuHiggs
	DOUBLE PRECISION amuexp,amuqed,damuqed,amuhadlo,amuhadnlo
	DOUBLE PRECISION amulbl,amuweaklo,amufermnlo,amubosnlo,damusm
	DOUBLE PRECISION errexp,errqed,errhadlo,errhadnlo
	DOUBLE PRECISION errlbl,errfermnlo,errbosnlo,errsm
	DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
	DOUBLE PRECISION MS,MC,MBNP,MB,MT,MTAU,MMUON,MZ,MW
	DOUBLE PRECISION G1Q,G2Q,GQ,ALSQ
	DOUBLE PRECISION QSTSB
	DOUBLE PRECISION HTQ,HBQ,MTOPQ,MBOTQ
	DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
	DOUBLE PRECISION LQ,KQ,ALQ,AKQ,MUQ,NUQ
	DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2)
	DOUBLE PRECISION PCOMP(2,2),CMASS
	DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),N(5,5)
	DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
	DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
	DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
	DOUBLE PRECISION amu1L
	DOUBLE PRECISION delmagmu,damumin,damumax,amuthmax,amuthmin
	
	COMMON/ALEM0/ALEM0
	COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
	COMMON/SMSPEC/MS,MC,MBNP,MB,MT,MTAU,MMUON,MZ,MW
	COMMON/QGAUGE/G1Q,G2Q,GQ,ALSQ
	COMMON/STSBSCALE/QSTSB
	COMMON/QQUARK/HTQ,HBQ,MTOPQ,MBOTQ
	COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
	COMMON/QNMPAR/LQ,KQ,ALQ,AKQ,MUQ,NUQ
	COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
	COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,N
	COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     C		MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     C		CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
        COMMON/MAGMU/delmagmu,damumin,damumax,amuthmax,amuthmin
    
	pi=4.d0*datan(1.d0)

	TANB=PAR(3)
        sinb=tanb/dsqrt(1.d0+tanb**2)
        cosb=sinb/tanb

********************************************************************
C	   Evaluation of the Discrepancy between Experiment and 
C          SM Contributions
C          (except for 2-loop bosonic: (-2.2 +/- 0.2)d-10: see below)

C	       -> Experimental Measurement BNL:
	amuexp=11659208.0d-10
	errexp=(6.3d-10)**2
	
C	       -> QED Contributions (4 loops):
	amuqed=11658471.8113d-10
	errqed=(0.0162d-10)**2
	
C		-> Difference Exp.-QED:
	damuqed=amuexp-amuqed
	
C		-> LO Hadronic contr. from [1]:
	amuhadlo=692.1d-10
	errhadlo=(5.6d-10)**2
	
C		-> NLO Hadronic contr. from [1]:
	amuhadnlo=-10.03d-10
	errhadnlo=(.22d-10)**2
	
C		-> Light-by-Light contr. from [2]:
	amulbl=11.0D-10
	errlbl=(4.0D-10)**2
	
C		-> Weak SM contr., 1 Loop (error negl.)
	amuweaklo=19.482d-10
	
C		-> Weak SM fermionic 2 loop contr. from [3]
	amufermnlo=-1.512d-10
	errfermnlo=(.1d-10)**2
	
C	        -> Weak SM bosonic 2 loop contr., e.g. from [8]
C       (Serves just to estimate pure SM contribution.
C        It will be deduced from the 2-loop bosonic contr. below
C        in order to avoid double counting.)
	amubosnlo=GF*MMUON**2*ALEMMZ/(72.d0*dsqrt(2.d0)*pi**3)
     C        *(107.d0+23.d0*(1.d0-4.d0*s2TW)**2)*dlog(MMUON/MW)
        errbosnlo=(0.2D-10)**2
	
C		-> Difference Exp.-SM:
	damusm=damuqed-amuhadlo-amuhadnlo-amulbl-amuweaklo-amufermnlo
     C    -amubosnlo
     
C		-> 1-sigma error of this difference:
	errsm=dsqrt(errexp+errqed+errhadlo+errhadnlo+errlbl+errfermnlo
     C     +errbosnlo)
     
C		-> 2-sigma lower bound to the required SUSY contr.:
	damumin=damusm-2.D0*errsm
     
C		-> 2-sigma upper bound to the required SUSY contr.:
	damumax=damusm+2.D0*errsm
	
********************************************************************
    
C	 - SMuon masses and mixing matrix
        MSMU(1)=MSMU1
        MSMU(2)=MSMU2

	USMU(1,1)=CSMU
	USMU(1,2)=dsqrt(1.d0-CSMU**2)
	USMU(2,1)=-USMU(1,2)
	USMU(2,2)=USMU(1,1)

C	 - 1-Loop Chargino/SNeutrino Contribution (cf. [4]):
	Ymu=MMUON/H2Q
	aux=0.d0
	do k=1,2
	aux=aux+MMUON/(12.d0*MSMUNT**2)
     C     *(g2q*V(k,1)**2+(Ymu*U(k,2))**2)*Fc1((MCH(k)/MSMUNT)**2)
     C             +2.d0*MCH(k)/(3.d0*MSMUNT**2)
     C       *(-dsqrt(g2q)*V(k,1))*(Ymu*U(k,2))*Fc2((MCH(k)/MSMUNT)**2)
	enddo
	amuChar=MMUON/(16.d0*pi**2)*aux

C	 - 1-Loop Neutralino/SMuon Contribution (cf. [4]):
	aux=0.d0
	do i=1,5
	 do k=1,2
	 aux=aux-MMUON/(12.d0*MSMU(k)**2)*((1/dsqrt(2.d0)
     C          *(dsqrt(g1q)*N(i,1)+dsqrt(g2q)*N(i,2))*USMU(k,1)
     C             -Ymu*N(i,4)*USMU(k,2))**2
     C            +(dsqrt(2.d0*g1q)*N(i,1)*USMU(k,2)
     C             +Ymu*N(i,4)*USMU(k,1))**2)*fn1((MNEU(i)/MSMU(k))**2)
     C            +MNEU(i)/(3.d0*MSMU(k)**2)*(1/dsqrt(2.d0)
     C        *(dsqrt(g1q)*N(i,1)+dsqrt(g2q)*N(i,2))*USMU(k,1)
     C            -Ymu*N(i,4)*USMU(k,2))
     C        *(dsqrt(2.d0*g1q)*N(i,1)*USMU(k,2)
     C            +Ymu*N(i,4)*USMU(k,1))*fn2((MNEU(i)/MSMU(k))**2)
	 enddo
	enddo
	amuNeutr=MMUON/(16.d0*pi**2)*aux

C	 - 1-Loop Higgs/Muon Contribution (cf. [5,6]):
	aux=0.d0
	do i=1,3
	aux=aux+SCOMP(i,2)**2*fhs(SMASS(i)/MMUON)
	enddo
	do i=1,2
	aux=aux-PCOMP(i,1)**2*sinb**2*fhp(PMASS(i)/MMUON)
	enddo
	aux=aux+sinb**2*fhc(CMASS/MMUON)
	amuHiggs=(Ymu/(4.d0*pi))**2*aux
	
	amu1L=amuChar+amuNeutr+amuHiggs
	
C	 - Large Logarithms: 1-loop-like 2-loop contributions (cf. [7]):
	QNP=Max(MCH(2),MSMU(2))
	amu1L=amu1L*(1.d0-4.d0*ALEM0/pi*dlog(QNP/MMUON))
	
C	 - 2-loop bosonic contributions (Leading Log) (cf. [8]):
	aux=0.d0
	do i=1,3
	aux=aux+SCOMP(i,2)*(SCOMP(i,2)*cosb-SCOMP(i,1)*sinb)/SMASS(i)**2
	enddo
	aux=(cosb**2-sinb**2)*MZ**2/cosb*aux
	aux=(98.d0+9.d0*aux+23.d0*(1.d0-4.d0*s2TW)**2)/30.d0
	amu2Lbos=5.d0/(24.d0*dsqrt(2.d0)*pi**3)*GF*MMUON**2*ALEMMZ
     C           *(aux*dlog(MMUON**2/MW**2))
	
C	 - 2-loop Higgs/Fermion contribution (cf. [9]):
	mtq=MT/(1.d0+4.d0/(3.d0*pi)*asf(MT)+11.d0/pi**2*asf(MT)**2)
	mbq=runmb(MB)
	aux=0.d0
	do i=1,3
	aux=aux+4.d0/3.d0*SCOMP(i,1)*SCOMP(i,2)/sinb
     C                                      *FS((MTQ/SMASS(i))**2)
     C         +1.d0/3.d0*SCOMP(i,2)**2/cosb*FS((MBQ/SMASS(i))**2)
     C         +SCOMP(i,2)**2/cosb*FS((mtau/SMASS(i))**2)
	enddo
	aux=aux/cosb
	do i=1,2
	aux=aux+PCOMP(i,1)**2*(4.d0/3.d0*FPS((MTQ/PMASS(i))**2)
     C                        +1.d0/3.d0*FPS((MBQ/PMASS(i))**2)*tanb**2
     C                        +FPS((mtau/PMASS(i))**2)*tanb**2)
	enddo
	amu2LHf=GF*MMUON**2*ALEM0/(4.d0*dsqrt(2.d0)*pi**3)*aux
	
C	 - 2-loop Sfermion contribution (cf. [10]):
	AT=PAR(12)
	AB=PAR(13)
	Atau=PAR(14)
	lambda=PAR(1)
	Yl=MTAU/H2Q
	
	MST(1)=MST1
	MST(2)=MST2
	UST(1,1)=CST
	UST(1,2)=dsqrt(1.d0-CST**2)
	UST(2,2)=UST(1,1)
	UST(2,1)=-UST(1,2)
	
	MSB(1)=MSB1
	MSB(2)=MSB2
	USB(1,1)=CSB
	USB(1,2)=dsqrt(1.d0-CSB**2)
	USB(2,2)=USB(1,1)
	USB(2,1)=-USB(1,2)
	
	MSTAU(1)=MSL1
	MSTAU(2)=MSL2
	USL(1,1)=CSL
	USL(1,2)=dsqrt(1.d0-CSL**2)
	USL(2,2)=USL(1,1)
	USL(2,1)=-USL(1,2)
	
	aux=0.d0
	do i=1,3
	do k=1,2
	lambdHT(i,k)=HTQ*(At*SCOMP(i,1)-muq*SCOMP(i,2)
     C                       -lambda*H2Q*SCOMP(i,3))*UST(k,2)*UST(k,1)
     C               +(HTQ**2*H1Q*SCOMP(i,1)-g1q/3.d0
     C                  *(H1Q*SCOMP(i,1)-H2Q*SCOMP(i,2)))*UST(k,2)**2
     C               +(HTQ**2*H1Q*SCOMP(i,1)-(3.d0*g2q-g1q)/12.d0
     C                  *(H1Q*SCOMP(i,1)-H2Q*SCOMP(i,2)))*UST(k,1)**2
	lambdHB(i,k)=HBQ*(Ab*SCOMP(i,2)-muq*SCOMP(i,1)
     C                       -lambda*H1Q*SCOMP(i,3))*USB(k,2)*USB(k,1)
     C               +(HBQ**2*H2Q*SCOMP(i,2)+g1q/6.d0
     C                  *(H1Q*SCOMP(i,1)-H2Q*SCOMP(i,2)))*USB(k,2)**2
     C               +(HBQ**2*H2Q*SCOMP(i,2)+(3.d0*g2q+g1q)/12.d0
     C                  *(H1Q*SCOMP(i,1)-H2Q*SCOMP(i,2)))*USB(k,1)**2
	lambdHL(i,k)=Yl*(Atau*SCOMP(i,2)-muq*SCOMP(i,1)
     C                       -lambda*H1Q*SCOMP(i,3))*USL(k,1)*USL(k,2)
     C               +(Yl**2*H2Q*SCOMP(i,2)+g1q/2.d0
     C                  *(H1Q*SCOMP(i,1)-H2Q*SCOMP(i,2)))*USL(k,2)**2
     C               +(Yl**2*H2Q*SCOMP(i,2)+(g2q-g1q)/4.d0
     C                  *(H1Q*SCOMP(i,1)-H2Q*SCOMP(i,2)))*USL(k,1)**2
	
	lambdHT(i,k)=lambdHT(i,k)*dsqrt(2.d0)*2.d0*MW/dsqrt(g2)
	lambdHB(i,k)=lambdHB(i,k)*dsqrt(2.d0)*2.d0*MW/dsqrt(g2)
	lambdHL(i,k)=lambdHL(i,k)*dsqrt(2.d0)*2.d0*MW/dsqrt(g2)
	
	aux=aux+SCOMP(i,2)*(4.d0/3.d0*lambdHT(i,k)
     C                      *FSF((MST(k)/SMASS(i))**2)/MST(k)**2
     C                     +1.d0/3.d0*lambdHB(i,k)
     C                      *FSF((MSB(k)/SMASS(i))**2)/MSB(k)**2
     C                     +lambdHL(i,k)
     C                      *FSF((MSTAU(k)/SMASS(i))**2)/MSTAU(k)**2)
	enddo
	enddo
	amu2LSF=GF*MMUON**2*ALEMMZ/(4.d0*dsqrt(2.d0)*pi**3*cosb)*aux
	
C	 - 2-loop Chargino contribution (cf. [11]):
	aux=0.d0
	do k=1,2
	do i=1,3
	aux=aux+SCOMP(i,2)/cosb*dsqrt(2.d0)*MW/MCH(k)
     C          *((U(k,1)*V(k,2)*SCOMP(i,1)+U(k,2)*V(k,1)*SCOMP(i,2))
     C                    +lambda/dsqrt(g2)*U(k,2)*V(k,2)*SCOMP(i,3))
     C          *fS((MCH(k)/SMASS(i))**2)
	enddo
	do i=1,2
	aux=aux-PCOMP(i,1)*tanb*dsqrt(2.d0)*MW/MCH(k)
     C          *((U(k,1)*V(k,2)*cosb+U(k,2)*V(k,1)*sinb)*PCOMP(i,1)
     C                    -lambda/dsqrt(g2)*PCOMP(i,2)*U(k,2)*V(k,2))
     C          *fPS((MCH(k)/PMASS(i))**2)
	enddo
	enddo
	amu2LCh=GF*MMUON**2*dsqrt(2.d0)*ALEMMZ/(8.d0*pi**3)*aux
	
C	 - Conclusion and errors  (cf. [11])
C              -> Pure SUSY contribution to the discrepancy (G-2)mu:
C                 (Here the bosonic SM 2-loop contr. is deduced)
	amu2L=amu1L+amu2Lbos+amu2LHf+amu2LSF+amu2LCh-amubosnlo
C              -> Theoretical errors (added linearly):
	amuerr=2.8d-10+.02d0*dabs(amu1L)+.3d0*dabs(amu2L-amu1L)
C	       -> Theoretical Central Value:
	delmagmu=amu2L
C	       -> Theoretical Error Bounds:
	amuthmax=amu2L+amuerr
	amuthmin=amu2L-amuerr
	
C	 - Comparison with the required SUSY contr.:
	IF(damumin-amuthmax.ge.0.d0)THEN
	PROB(37)=damumin-amuthmax
	ENDIF
	IF(damuMax-amuthmin.le.0.d0)THEN
	PROB(37)=damuMax-amuthmin
	ENDIF

        return
        END
	
	
C*********************************************************************	
C*********************************************************************
C*********************************************************************

	double precision function fn1(x)

	implicit none
	double precision x
	if(dabs(x).gt.1.d-3)then
	if(dabs(x-1.d0).gt.1.d-3)then
	fn1=2.d0/(1.d0-x)**4
     C         *(1.d0-6.d0*x+3.d0*x**2+2.d0*x**3-6.d0*x**2*dlog(x))
	else
	fn1=1.d0
	endif
	else
	fn1=2.d0
	endif
	return
	end

C********************************************************************

	double precision function fn2(x)

	implicit none
	double precision x
	if(dabs(x).gt.1.d-3)then
	if(dabs(x-1.d0).gt.1.d-3)then
	fn2=3.d0/(1.d0-x)**3*(1.d0-x**2+2.d0*x*dlog(x))
	else
	fn2=1.d0
	endif
	else
	fn2=3.d0
	endif
	return
	end

C********************************************************************

	double precision function fc1(x)

	implicit none
	double precision x
	if(dabs(x).gt.1.d-3)then
	if(dabs(x-1.d0).gt.1.d-3)then
	fc1=2.d0/(1.d0-x)**4
     C         *(2.d0+3.d0*x-6.d0*x**2+x**3+6.d0*x*dlog(x))
	else
	fc1=1.d0
	endif
	else
	fc1=4.d0
	endif
	return
	end

C********************************************************************

	double precision function fc2(x)

	implicit none
	double precision x
	if(dabs(x-1.d0).gt.1.d-3)then
	fc2=3.d0/(1.d0-x)**3*(-3.d0+4.d0*x-x**2-2.d0*dlog(x))/2.d0
	else
	fc2=1.d0
	endif
	return
	end

C********************************************************************

	double precision function fhs(x)

	implicit none
	double precision x,aux
	if(abs(x).gt.2.001d0)then
	 if(abs(x).gt.10.d0)then
	 aux=(3.d0+x**2)*dlog(x**2)/x**4-7.d0/(6.d0*x**2)
	 else
	 aux=3.d0/2.d0-x**2+x**2/2.d0*(x**2-3.d0)*dlog(x**2)
     C       +(5.d0*x**4-x**6-4.d0*x**2)/(2.d0*dsqrt(x**2*(x**2-4.d0)))
     C       *dlog((2.d0-x**2-dsqrt(x**2*(x**2-4)))
     C              *(x**2-dsqrt(x**2*(x**2-4.d0)))
     C              /((2.d0-x**2+dsqrt(x**2*(x**2-4.d0)))
     C                *(x**2+dsqrt(x**2*(x**2-4.d0)))))
	 endif
	else
	aux=-5.d0/2.d0+4.d0*dlog(2.d0)
	endif
	fhs=aux
	return
	end

C********************************************************************

	double precision function fhp(x)

	implicit none
	double precision x,aux
	if(abs(x).gt.2.001d0)then
	 if(abs(x).gt.15.d0)then
	 aux=(5.d0+x**2)*dlog(x**2)/x**4-11.d0/(6.d0*x**2)
	 else
	 aux=1.d0/2.d0+x**2+x**2/2.d0*(1.d0-x**2)*dlog(x**2)
     C       +x**4*(x**2-3.d0)/(2.d0*dsqrt(x**2*(x**2-4.d0)))
     C       *dlog((2.d0-x**2-dsqrt(x**2*(x**2-4)))
     C              *(x**2-dsqrt(x**2*(x**2-4.d0)))
     C              /((2.d0-x**2+dsqrt(x**2*(x**2-4.d0)))
     C                *(x**2+dsqrt(x**2*(x**2-4.d0)))))
	 endif
	else
	aux=17.d0/2.d0-12.d0*dlog(2.d0)
	endif
	fhp=aux
	return
	end

C********************************************************************

	double precision function fhc(x)

	implicit none
	double precision x,aux
	if(abs(x).gt.2.001d0)then
	 if(abs(x).gt.10.d0)then
	 aux=-(1.d0+2.d0*x**2)/(12.d0*x**4)
	 else
	 aux=1.d0/2.d0-x**2+x**2*(x**2-1.d0)*dlog(x**2/(x**2-1.d0))
	 endif
	else
	aux=-7.d0/2.d0+24.d0*dlog(2.d0)-12.d0*dlog(3.d0)
	endif
	fhc=aux
	return
	end

C********************************************************************

	double precision function fPS(x)

	implicit none
	double precision x,aux
	IF(x.lt.0.002d0)THEN
	aux=0.4995501108d0*dsqrt(x)+1133.437039d0*x
     C      -371.1353082d0*x*dsqrt(x)-1087.028512d0*dlog(1.d0+x)
	ELSE
	 IF(x.lt.0.02d0)THEN
	 aux=-1.018546244d-3+7.096260514d-1*dsqrt(x)
     C     +38.60712022d0*x-331.8418712d0*x*dsqrt(x)
     C     +1522.657713d0*x**2-2874.048399d0*x**2*dsqrt(x)
	 ELSE
	  IF(x.lt.0.1d0)THEN
	  aux=4.306498122d-2+32.56225304d0*dlog(1.d0+x)
     C        +925.7514262d0*dlog(1.d0+x)**2
     C        -148.1323591d0*dsqrt(x)*dlog(1.d0+x)
     C        -648.4904538d0*x*dlog(1.d0+x)
	  ELSE
	   IF(x.lt.2.5d0)THEN
	   aux=-1.365013496d-1-1.858455623d0*dlog(1.d0+x)
     C         -5.996763746d-1*dlog(1.d0+x)**2
     C         +4.390843985d-1*dsqrt(x)*dlog(1.d0+x)
     C         -1.444359743d-1*x*dlog(1.d0+x)
     C         +3.852425143d0*dsqrt(x)
	   ELSE
	    IF(x.lt.100.d0)THEN
	    aux=4.304425955d-1+6.766323794d-2*dlog(1.d0+x)
     C          -1.584446296d-1*dlog(1.d0+x)**2
     C          -2.787080541d-1*dsqrt(x)*dlog(1.d0+x)
     C          +1.557845370d-3*x*dlog(1.d0+x)
     C          +2.139180566d0*dsqrt(x)
	    ELSE
	     IF(x.lt.10000.d0)THEN
	     aux=2.025445594d0+9.960866255d-1*dlog(x)
     C           +1.122896720d-4*dsqrt(x)
	     ELSE
	      aux=2.000835136d0+9.9992369d-1*dlog(x)
     C            +2.327105016d-7*dsqrt(x)
	     ENDIF
	    ENDIF
	   ENDIF
	  ENDIF
	 ENDIF
	ENDIF
	fPS=aux
	return
	end

C********************************************************************

	double precision function fS(x)

	implicit none
	double precision x,fPS
	fS=(2.d0*x-1.d0)*fPS(x)-2.d0*x*(2.d0+dlog(x))
	return
	end
		
C********************************************************************

	double precision function fSF(x)

	implicit none
	double precision x,fPS
	fSF=x/2.d0*(2.d0+dlog(x)-fPS(x))
	return
	end
	
