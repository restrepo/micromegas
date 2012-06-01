C           Subroutine for BR(b -> s gamma)


C Literature:
C   - THEORETICAL FORMULAE:
C [1] K.Chetyrkin, M.Misiak, M.Munz,
C     'Weak Radiative B-Meson Decay Beyond Leading Logarithms' 
C     Phys.Lett.B400:206-219,1997, 
C     Erratum-ibid.B425:414,1998, e-Print: hep-ph/9612313
C
C [2] M.Ciuchini, G.Degrassi, P.Gambino, G.Giudice
C     'Next to Leading QCD Corrections to B -> Xs gamma: Standard Model 
C     and Two-Higgs Doublet Model' Nucl.Phys.B534:3-20,1998, 
C     e-Print: hep-ph/9806308
C
C [3] G.Degrassi, P.Gambino, G.Giudice
C     'B -> Xs gamma in Supersymmetry: Large Contributions Beyond the 
C     Leading Order', JHEP 0012:009,2000, e-Print: hep-ph/0009337
C
C [4] P.Gambino, M.Misiak, 'Quark Mass Effects in B -> Xs gamma'
C     Nucl.Phys.B611:338-366,2001, e-Print: hep-ph/0104034
C
C [5] A.Buras, A.Czarnecki, M.Misiak, J.Urban,
C     'Completing the NLO QCD Calculation of B -> Xs gamma'
C     Nucl.Phys.B631:219-238,2002, e-Print: hep-ph/0203135
C
C [6] G.Belanger, F.Boudjema, A.Pukhov, A.Semenov,
C     'micrOMEGAs: Version 1.3', Appendix B,
C     Comput.Phys.Commun.174:577-604,2006, e-Print: hep-ph/0405253
C
C [7] T.Hurth, E.Lunghi, W.Porod,
C     'Untagged B -> Xs+d gamma CP asymmetry as a probe for New Physics'
C     Nucl.Phys.B704:56-74,2005, e-Print: hep-ph/0312260
C
C [8] M.Misiak et al.,
C     'Estimate of B(anti-B ---> X(s) gamma) at O(alpha(s)**2)'
C     Phys.Rev.Lett.98:022002,2007, e-Print: hep-ph/0609232
C
C [9] T.Becher and M.~Neubert,
C     Phys.\ Rev.\ Lett.\  {\bf 98} (2007) 022003 [arXiv:hep-ph/0610067].
C
C [10] A. J. Buras, P. H. Chankowski, J. Rosiek, L. Slawianowska,
C     ' Delta M(d, s), B0(d, s) ---> mu+ mu- and B ---> X(s) gamma 
C      in supersymmetry at large tan beta.'
C     Nucl.Phys.B659:3,2003, e-Print: hep-ph/0210145
C
C [11] C.Bobeth, A.J.Buras, F.Kruger and J.Urban,
C     'QCD corrections to anti-B --> X/d,s nu anti-nu, anti-B/d,s --> l+ l-,  
C      K--> pi nu anti-nu and K(L) --> mu+ mu- in the MSSM,'
C     Nucl.\ Phys.\  B {\bf 630} (2002) 87 [arXiv:hep-ph/0112305].
C
C [12] A.G. Akeroyd, S. Recksiegel,
C     ' The Effect of H+- on B+- ---> tau+- nu(tau) and 
C       B+- ---> mu+- muon neutrino'
C     J.Phys.G29:2311-2317,2003, e-Print: hep-ph/0306037
C
C
C   - SOURCES FOR EXPERIMENTAL/LATTICE QCD DATA:
C [13] E. Barberio et al. [Heavy Flavor Averaging Group (HFAG)],
C     ' Averages of b-hadron properties at the end of 2005'
C     arXiv:hep-ex/0603003, and\newline {\tt www.slac.stanford.edu/xorg/hfag}.
C
C [14] A. Abulencia et al. [CDF Collaboration],
C      Phys.\ Rev.\ Lett.\  {\bf 97} (2006) 242003
C      arXiv:hep-ex/0609040.
C
C [15] A. Gray et al. [HPQCD Collaboration],
C     ' The B meson decay constant from unquenched lattice QCD,'
C     Phys.\ Rev.\ Lett.\  {\bf 95} (2005) 212001 [arXiv:hep-lat/0507015].
C
C [16] E. Dalgic et al.,
C      Phys.\ Rev.\  D {\bf 76} (2007) 011501 [arXiv:hep-lat/0610104].
C
C [17] M. Okamoto,
C      'Full determination of the CKM matrix using recent results from lattice QCD,'
C      PoS {\bf LAT2005} (2006) 013 [arXiv:hep-lat/0510113].
C
C [18] P. Ball and R. Fleischer,
C      Eur.\ Phys.\ J.\  C {\bf 48} (2006) 413 [arXiv:hep-ph/0604249].
C
C [19] CDF Collaboration,
C      'Search for $B_s^0 \to \mu^+ \mu^-$ and $B_d^0 \to \mu^+ \mu^-$
C      Decays in 2 fb$^{-1}$ of $p\bar{p}$ Collisions with CDF II'', 
C      CDF Public Note 8956 (2007)


        SUBROUTINE BSG(PAR,PROB)
	IMPLICIT NONE

	INTEGER I,J,K
	DOUBLE PRECISION PAR(*),PROB(*)
	DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
	DOUBLE PRECISION MS,MC,MBNP,MB,MT,MTAU,MMUON,MZ,MW
	DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2)
	DOUBLE PRECISION PCOMP(2,2),CMASS
	DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),N(5,5)
	DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
	DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
	DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
	DOUBLE PRECISION TANB,COSB,SINB,au,ad,SST,SSB
	DOUBLE PRECISION ST(2),RST(2,2),SB(2),RSB(2,2),CCD(2),CCT(2,2)
	DOUBLE PRECISION FF1,FF2,FF3,FG1,FG2,FG3,ffh
	DOUBLE PRECISION fgh,esm,eh
	DOUBLE PRECISION asf,H2,C70SM,C70HIG
	DOUBLE PRECISION C80SM,C80HIG
	DOUBLE PRECISION C71HIG,C81HIG,dC7SM,dC8SM,dC7HIG
	DOUBLE PRECISION dC8HIG,C7CHARS,C8CHARS,C7CHAR,C8CHAR
	DOUBLE PRECISION XT,YT,XSQC(2),XSTC(2,2),PI,AAT,AAB,mu
	DOUBLE PRECISION xQg,xBg1,xBg2,xTg1,xTg2,xTneu(2,5),xBneu(2,5)
	DOUBLE PRECISION akk,etaS,epsb,epsbp,epst
	DOUBLE PRECISION MT0,MTH,etaH,fBssqBs,fBdsqBd,sqBBs
	DOUBLE PRECISION QSTSB
	DOUBLE PRECISION C70,C80,ALEM0,BRSL
	DOUBLE PRECISION eta,z,MB1S,etaB
	DOUBLE PRECISION delt,delt2,lndelt,lndeltp
	DOUBLE PRECISION C20b,C70b,C80b,C7EMb,C80S0,eta0
	DOUBLE PRECISION lambd2,gg1,gg2
	DOUBLE PRECISION ff11,ff12,ff17,ff18,ff22,ff27,ff28,ff77,ff78
	DOUBLE PRECISION ff88,Mch2
	DOUBLE PRECISION aa(8),bb(4),hh(8),h8(4),ee(8)
	DOUBLE PRECISION C7HIG,C8HIG,C41SM,C41HIG
	DOUBLE PRECISION dd(8),dt(8),dte(8),dim(8),da(8),db(8)
	DOUBLE PRECISION C10b,C70BSM,C80BSM
	DOUBLE PRECISION KC0,KT0,KBSM0,KC,KT,KCIM,KTIM,Ktot
	DOUBLE PRECISION KIM,KBSM,KBSMIM
	DOUBLE PRECISION sc0,VVtu,VVtuim,Vbsg,rmu
	DOUBLE PRECISION EPSew,HQET,BREMS,CCSL
	DOUBLE PRECISION af,bf,At1,A1,Ft1,F1,dAto,dA0,dFto,dF0
	DOUBLE PRECISION afim,bfim,sp2,aux
	DOUBLE PRECISION HTQ,HBQ,MTOPQ,MBOTQ
	DOUBLE PRECISION asmt,asmh,asmsusy,asc0,asmc,asmb
	DOUBLE PRECISION DC70BSM,DC80BSM,DKBSM
	DOUBLE PRECISION BRSG1,BRSG2
	DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
	DOUBLE PRECISION DSM,DH,Dchid,Dchis,DDP,ffp,ggp,gg0,S0
	DOUBLE PRECISION VVc,VVu,Vtb2,Vcs2,Vud2,VtdVtb2,VtsVtb2,runmb
	DOUBLE PRECISION MQU,MD,MS0,MC0,MB0,mmu,scR,scR2,runmass
	DOUBLE PRECISION epst0,epst1,epst2,epst3
	DOUBLE PRECISION epsY32,epsY31,epsY23,epsY13
	DOUBLE PRECISION DMdexpmin,DMdexpMax,DMsexpmin,DMsexpMax,
     C                   BRBMUMUexpMax,BRBTAUNUexpmin,BRBTAUNUexpMax,
     C                   BRSGexpmin,BRSGexpMax
	DOUBLE PRECISION sigRLbs,sigRLbd,sigLRbs,sigLRbd,BB0,BB1
	DOUBLE PRECISION C2LR,C1SLL,C1SRR,Ca,Cs,Cp,C70S0,sgn
	DOUBLE PRECISION Vub2,fB,mBu,tauB,rh
	DOUBLE PRECISION BRJJ(5),BRMM(5),BRLL(5),BRSS(5),BRCC(5)
	DOUBLE PRECISION BRBB(5),BRTT(5),BRWW(3),BRZZ(3),BRGG(5)
	DOUBLE PRECISION BRZG(5),BRHHH(4),BRHAA(3,3),BRHCHC(3)
	DOUBLE PRECISION BRHAZ(3,2),BRAHA(3),BRAHZ(2,3),BRHCW(5)
	DOUBLE PRECISION BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
	DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,6),BRASL(2,3)
	DOUBLE PRECISION BRSUSY(5),WIDTH(5)
	DOUBLE PRECISION DBRSGmax,DBRSGmin
	DOUBLE PRECISION dVub,dVtdVtb2,dVtsVtb2,dfB,dfBssqBs,dfBdsqBd
	DOUBLE PRECISION BRSG,BRSGmax,BRSGmin,DMd,DMdmin,DMdmax,DMs,
     C        DMsmax,DMsmin,BRBMUMU,BRBMUMUmax,BRBMUMUmin,BRBtaunu,
     C        BRBtaunumax,BRBtaunumin
	DOUBLE PRECISION BRBSll,BRBSllmin,BRBSllmax,CQ12,CQ22,Intpropa
	
	
	COMMON/ALEM0/ALEM0
	COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
	COMMON/SMSPEC/MS,MC,MBNP,MB,MT,MTAU,MMUON,MZ,MW
	COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
	COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,N
	COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     C		MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     C		CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
 	COMMON/STSBSCALE/QSTSB
	COMMON/BRSG/BRSG,BRSGmax,BRSGmin,DMd,DMdmin,DMdmax,DMs,
     C        DMsmax,DMsmin,BRBMUMU,BRBMUMUmax,BRBMUMUmin,BRBtaunu,
     C        BRBtaunumax,BRBtaunumin
	COMMON/QQUARK/HTQ,HBQ,MTOPQ,MBOTQ
	COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
	COMMON/BRN/BRJJ,BRMM,BRLL,BRSS,BRCC,BRBB,BRTT,BRWW,BRZZ,
     C		BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     C		BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     C		BRSUSY,WIDTH
     
        pi=4.d0*datan(1.d0)

	AAT=PAR(12)
	AAB=PAR(13)
	mu=PAR(4)
	Mch2=PAR(21)
	gg1=dsqrt(g1)
	gg2=dsqrt(g2)
	
C       Alpha_s at various scales:
C           Charm Quark Mass:
         asmc=asf(mc)
C           Bottom Quark Mass:
         asmb=asf(mb)	 
C           M_top:
         asmt=asf(MT)
C           Charged Higgs Mass:
         asmh=asf(CMASS)
C           Susy scale (squark masses):
         asmsusy=asf(dsqrt(QSTSB))

C**********************************************************************
C**********************************************************************
C      MASSES AND PARAMETERS
	
C       Matching Scale sc0 = m_top(m_top)(MSbar)
        MT0=MT/(1.d0+4.d0/(3.d0*pi)*asmt+11.d0/pi**2*asmt**2)
        sc0=MT0
	
C       Alphas at the matching scale:
	asc0=asf(sc0)

C       Trig. Functions of Beta
        TANB=PAR(3)
        sinb=tanb/dsqrt(1.d0+tanb**2)
        cosb=sinb/tanb
        au=1.d0/tanb
        ad=-tanb
	
C       M_top at the Charged Higgs Mass	Scale:
	mth=mt*(asmh/asmt)**(4.d0/7.d0)/(1.d0+4.d0/(3.d0*pi)*asmt)
	        
C       Stop Masses and Mixing Angles
	ST(1)=MST1
	ST(2)=MST2
        SST=DSQRT(1.d0-CST**2)
	RST(1,1)=CST
	RST(2,2)=CST
	RST(1,2)=SST
	RST(2,1)=-SST

C       Sbottom Masses and Mixing Angles
	SB(1)=MSB1
	SB(2)=MSB2
        SSB=DSQRT(1.d0-CSB**2)
	RSB(1,1)=CSB
	RSB(2,2)=CSB
	RSB(1,2)=SSB
	RSB(2,1)=-SSB

C       x Parameters (Squares of Mass Ratios)

        xt= (MT0/MW)**2                  !(m_t/m_W)^2
        yt=(MTH/CMASS)**2                !(m_t/m_H+)^2
		        
	do i=1,2
        xsqc(i)=(MUL/MCH(i))**2           !(m_Q/m_ch(i))^2
        do k=1,2
        xstc(k,i)=(ST(k)/MCH(i))**2       !(m_St(k)/m_ch(i))^2
        enddo
        enddo 
         
        xQg=(MUL/MGL)**2                  !(m_Q/m_gl)^2
	xBg1=(SB(1)/MGL)**2               !(m_Sb(1)/m_gl)^2
        xBg2=(SB(2)/MGL)**2               !(m_Sb(2)/m_gl)^2
        xTg1=(ST(1)/MGL)**2               !(m_St(1)/m_gl)^2
        xTg2=(ST(2)/MGL)**2               !(m_St(2)/m_gl)^2

        do i=1,5
        do k=1,2
        xTneu(k,i)=(ST(k)/MNEU(i))**2       !(m_St(k)/m_neu(i))^2
	xBneu(k,i)=(SB(k)/MNEU(i))**2       !(m_Sb(k)/m_neu(i))^2
	enddo
	enddo
	
C       CKM coefficients

	VVc=0.9879d0                     ! |V_cs.V_cb|/|V_ts.V_tb|
	VVu=0.4742d0                     ! |V_us.V_ub|/|V_ts.V_tb|
	Vtb2=(0.9991d0)**2                   ! (V_tb)^2
	Vcs2=(0.97296d0)**2                  ! (V_cs)^2
	Vud2=(0.97383d0)**2                  ! (V_ud)^2
	
C       Since (see, e.g., [18]) V_ub(incl) differs considerably from
C       V_ub(excl), we allow for the large range
C            3.3 10^-3 < V_ub < 4.7 10^-3:
C	(V_ub)^2
	Vub2=(4.0d-3)**2
C	uncertainty on V_ub
	dVub=0.7d-3

C       From [18]: V_tb*V_td=(8.6 +/- 2.8) 10^-3  (2sigma)
C	(V_tb*V_td)^2
	VtdVtb2=(8.6d-3)**2
C	uncertainty on (V_tb*V_td)^2
	dVtdVtb2=2.d0*dsqrt(VtdVtb2)*2.8d-3
C       From [18]: V_ts*V_tb=(41.3 +/- 1.4) 10^-3  (2sigma)
C	(Vtb.V_ts)^2
	VtsVtb2=(0.0413d0)**2
C	uncertainty on (Vtb.V_ts)^2
	dVtsVtb2=2.d0*dsqrt(VtsVtb2)*1.4d-3

C	(V_ts.V_tb/V_cb)^2
	Vbsg=VtsVtb2/(0.042d0)**2

C	V_us.V_ub/V_ts.V_tb
	VVtu=-0.011d0
	VVtuim=0.0180d0

C	Stop/Sbottom Scale
	scR=dsqrt(QSTSB)
	scR2=QSTSB
	
 	
C	Quark Masses at the Stop/Sbottom scale

	MS0=runmass(0.095d0,scR)
	MC0=runmass(1.25d0,scR)
	MB0=runmb(scR)
        MQU=runmass(0.002d0,scR)
	MD=runmass(0.005d0,scR)

C	Myon mass

	mmu=0.10566d0  !GeV
	
C	Hadronic parameters

C	From [16]:
C	f_Bs*sqrt(B_Bs)=(0.281 +/- 0.042) GeV  (2sigma)
	fBssqBs=0.281d0          ! GeV         f_Bs sqrt(B_Bs)
	dfBssqBs=0.042d0         ! GeV         uncertainty  (2sigma)

C	From [17]:
C	f_Bs*sqrt(B_Bs)/(f_Bd*sqrt(B_Bd))=1.216 +/- 0.0674 (2sigma)
C       so f_Bd*sqrt(B_Bd)=0.231*(1 +/- 0.1639) GeV  (2sigma)
	fBdsqBd=0.231d0            ! GeV         f_Bd*sqrt(B_Bd)
	dfBdsqBd=fBdsqBd*0.1639d0  ! GeV         uncertainty

C	From [15]:
C	f_B=(0.216 +/- 0.044) GeV (2sigma)
	fB=0.216d0               ! GeV         for B+ --> tau+ nu_tau
	dfB=0.044d0              ! GeV         uncertainty
	
C	Experimental data

C	Delta Md from [12], 2 sigma bounds:
C             0.499ps-1 < DMd=0.507ps-1 < 0.515ps-1
	DMdexpmin=0.499d0
	DMdexpMax=0.515d0

C	Delta Ms from [14], 2 sigma bounds: 
C              17.53ps-1 < DMs=17.77ps-1 < 18.01ps-1
	DMsexpmin=17.53d0
	DMsexpMax=18.01d0

C	From [19]: BR(Bs->mu+mu-) < 5.8 10^-8 (95% C.L.)
	BRBMUMUexpmax=5.8d-8
	
C       BR(B+ -> tau+ nu) from [12], 2 sigma bounds:
C	        0.34 10^-4 < BR(B+ -> tau+ nu)=1.32 10^-4 < 2.30 10^-4
	BRBTAUNUexpMax=2.3d-4
	BRBTAUNUexpmin=0.34d-4

C       BR(B -> Xs gamma) from [13], 2 sigma bounds:
C           3.03 10^-4 < BR(B -> Xs gamma)=3.52 10^-4 < 4.01 10^-4
	BRSGexpmin=3.03d-4
	BRSGexpMax=4.01d-4

C**********************************************************************
C	EFFECTIVE NEUTRAL HIGGS COUPLINGS
C       -> EPSILON COEFFICIENTS, notation following [10]:
C
C       epsilontilde_J as in eq. (5.1) in [10] (up to a factor tanb)
C         with Delta m_d as in eq.(2.5), and Sigma as in App. (A.2)
C
C       epsilon_Y^(JI) as in eq. (5.1) in [10] (up to a factor yt^2 
C         and a factor tanb), with (3.7) for lambda_0^(JI) and (3.53)
C         for the CKM matrix elements in terms of V^eff

C	 a) epst1 = (epsilontilde as in[10])*tanb, epst2 (* tanb)
	epst2=asf(scR)/(3.d0*pi)*(BB1(0.d0,MGL**2,MDL**2,scR2) 
     C      +BB1(0.d0,MGL**2,MDL**2,scR2))                  !gluinos

	aux=0.d0
	do i=1,5
	aux=aux+2.d0/H2Q*MNEU(i)                            !neutralinos
     C           *((gg1/3.d0*N(i,1)-gg2*N(i,2))*N(i,3)
     C                 *BB0(0.d0,MNEU(i)**2,MDL**2,scR2)
     C            +DSQRT(2.d0)/3.d0*gg1*N(i,1)*N(i,3)
     C                 *BB0(0.d0,MNEU(i)**2,MDL**2,scR2))
     C         +((MS0/H2Q)**2*N(i,3)**2
     C           +1.d0/2.d0*(gg1/3.d0*N(i,1)-gg2*N(i,2))**2)
     C                 *BB1(0.d0,MNEU(i)**2,MDL**2,scR2)
     C         +(2.d0/9.d0*gg1**2*N(i,1)**2
     C           +(MS0/H2Q)**2*N(i,3)**2)
     C                 *BB1(0.d0,MNEU(i)**2,MDL**2,scR2)
	enddo
	epst2=epst2+aux/(32.d0*pi**2)

	aux=0.d0
	do i=1,2
	aux=aux-2.d0*MCH(i)*gg2/H2Q                          !charginos
     C             *U(i,2)*V(i,1)*BB0(0.d0,MCH(i)**2,MUL**2,scR2)
     C         +(MC/H1Q)**2*V(i,2)**2
     C             *BB1(0.d0,MCH(i)**2,MUR**2,scR2)
     C         +((MS0/H2Q)**2*U(i,2)**2
     C           +gg2**2*V(i,1)**2)*BB1(0.d0,MCH(i)**2,MUL**2,scR2)
	enddo

	epst1=epst2+aux*Vud2/(32.d0*pi**2)
	epst2=epst2+aux*Vcs2/(32.d0*pi**2)
	
C	 b) epst3 (* tanb)
	epst3=asf(scR)/(3.d0*pi)                                     
     C *(BB1(0.d0,MGL**2,SB(1)**2,scR2)+BB1(0.d0,MGL**2,SB(2)**2,scR2)
     C -2.d0*MGL/MB0*SSB*CSB*(BB0(0.d0,MGL**2,SB(1)**2,scR2)
     C -BB0(0.d0,MGL**2,SB(2)**2,scR2)))                        !gluinos

	aux=0.d0
	do i=1,5
	aux=aux                                             !neutralinos
     C   +2.d0*((CSB*(gg1/3.d0*N(i,1)-gg2*N(i,2))
     C                     +HBQ*dsqrt(2.d0)*SSB*N(i,3))
     C        *(gg1/3.d0*SSB*N(i,1)+HBQ/dsqrt(2.d0)*CSB*N(i,3))
     C        *BB0(0.d0,MNEU(i)**2,SB(1)**2,scR2)
     C      -(SSB*(gg1/3.d0*N(i,1)-gg2*N(i,2))
     C                     -HBQ*dsqrt(2.d0)*CSB*N(i,3))
     C        *(gg1/3.d0*CSB*N(i,1)-HBQ/dsqrt(2.d0)*SSB*N(i,3))
     C        *BB0(0.d0,MNEU(i)**2,SB(2)**2,scR2))*MNEU(i)/MB0
     C   +((dsqrt(2.d0)*gg1/3.d0*SSB*N(i,1)+HBQ*CSB*N(i,3))**2
     C   +(CSB/dsqrt(2.d0)*(gg1/3.d0*N(i,1)-gg2*N(i,2))
     C                                  +HBQ*SSB*N(i,3))**2)
     C        *BB1(0.d0,MNEU(i)**2,SB(1)**2,scR2)
     C   +((dsqrt(2.d0)*gg1/3.d0*CSB*N(i,1)-HBQ*SSB*N(i,3))**2
     C   +(SSB/dsqrt(2.d0)*(gg1/3.d0*N(i,1)-gg2*N(i,2))
     C                                  -HBQ*CSB*N(i,3))**2)
     C        *BB1(0.d0,MNEU(i)**2,SB(2)**2,scR2)
	enddo
	epst3=epst3+1.d0/(32.d0*pi**2)*aux

	aux=0.d0
	do i=1,2                                            !charginos
	aux=aux-2.d0*MCH(i)/H2Q
     C    *(CST*U(i,2)*(gg2*CST*V(i,1)-HTQ*SST*V(i,2))
     C                       *BB0(0.d0,MCH(i)**2,ST(1)**2,scR2)
     C    +SST*U(i,2)*(gg2*SST*V(i,1)+HTQ*CST*V(i,2))
     C                       *BB0(0.d0,MCH(i)**2,ST(2)**2,scR2))
     C    +(HBQ**2*U(i,2)**2*CST**2+(gg2*CST*V(i,1)-HTQ*SST*V(i,2))**2)
     C                       *BB1(0.d0,MCH(i)**2,ST(1)**2,scR2)
     C    +(HBQ**2*U(i,2)**2*SST**2+(gg2*SST*V(i,1)+HTQ*CST*V(i,2))**2)
     C                       *BB1(0.d0,MCH(i)**2,ST(2)**2,scR2)
	enddo
	epst3=epst3+Vtb2/(32.d0*pi**2)*aux

C	 c) epsY32 (*Yt^2 tanb), epsY31 (*Yt^2 tanb)
	epsY32=0.d0
	do i=1,2
	epsY32=epsY32
     C        +(gg2*CST*V(i,1)-HTQ*SST*V(i,2))**2
     C                       *BB1(0.d0,MCH(i)**2,ST(1)**2,scR2)
     C        +(gg2*SST*V(i,1)+HTQ*CST*V(i,2))**2
     C                       *BB1(0.d0,MCH(i)**2,ST(2)**2,scR2)
     C        -2.d0*MCH(i)/H2Q*U(i,2)
     C                      *(CST*(gg2*CST*V(i,1)-HTQ*SST*V(i,2))
     C                           *BB0(0.d0,MCH(i)**2,ST(1)**2,scR2)
     C    		     +SST*(gg2*SST*V(i,1)+HTQ*CST*V(i,2))
     C                           *BB0(0.d0,MCH(i)**2,ST(2)**2,scR2))
     C        +(MS0/H2Q)**2*U(i,2)**2
     C         *(CST**2*BB1(0.d0,MCH(i)**2,ST(1)**2,scR2)
     C              +SST**2*BB1(0.d0,MCH(i)**2,ST(2)**2,scR2))
     C      +VVc*((gg2*V(i,1))**2*BB1(0.d0,MCH(i)**2,MUL**2,scR2)
     C            -2.d0*gg2*MCH(i)/H2Q*U(i,2)*V(i,1)
     C                              *BB0(0.d0,MCH(i)**2,MUL**2,scR2)
     C            +(MC/H1Q)**2*V(i,2)**2
     C                              *BB1(0.d0,MCH(i)**2,MUR**2,scR2)
     C            +(MS0/H2Q)**2*U(i,2)**2
     C                              *BB1(0.d0,MCH(i)**2,MUL**2,scR2))
	enddo
	epsY32=epsY32/(32.d0*pi**2)
	
	epsY31=0.d0
	do i=1,2
	epsY31=epsY31
     C        +(gg2*CST*V(i,1)-HTQ*SST*V(i,2))**2
     C                       *BB1(0.d0,MCH(i)**2,ST(1)**2,scR2)
     C        +(gg2*SST*V(i,1)+HTQ*CST*V(i,2))**2
     C                       *BB1(0.d0,MCH(i)**2,ST(2)**2,scR2)
     C        -2.d0*MCH(i)/H2Q*U(i,2)
     C                       *(CST*(gg2*CST*V(i,1)-HTQ*SST*V(i,2))
     C                          *BB0(0.d0,MCH(i)**2,ST(1)**2,scR2)
     C    		       +SST*(gg2*SST*V(i,1)+HTQ*CST*V(i,2))
     C                          *BB0(0.d0,MCH(i)**2,ST(2)**2,scR2))
     C        +(MD/H2Q)**2*U(i,2)**2
     C         *(CST**2*BB1(0.d0,MCH(i)**2,ST(1)**2,scR2)
     C              +SST**2*BB1(0.d0,MCH(i)**2,ST(2)**2,scR2))
     C      +VVu*((gg2*V(i,1))**2*BB1(0.d0,MCH(i)**2,MUL**2,scR2)
     C            -2.d0*gg2*MCH(i)/H2Q*U(i,2)*V(i,1)
     C                             *BB0(0.d0,MCH(i)**2,MUL**2,scR2)
     C            +(MQU/H1Q)**2*V(i,2)**2
     C                             *BB1(0.d0,MCH(i)**2,MUR**2,scR2)
     C            +(MD/H2Q)**2*U(i,2)**2
     C                            *BB1(0.d0,MCH(i)**2,MUL**2,scR2))
	enddo
	epsY31=epsY31/(32.d0*pi**2)
	
C	 d) epsY13 (*Yt^2 tanb), epsY23 (*Yt^2 tanb)
	epsY13=0.d0
	do i=1,2
	epsY13=epsY13-2.d0/H2Q*U(i,2)*MCH(i)
     C                 *(CST*(gg2*CST*V(i,1)-HTQ*SST*V(i,2))
     C                    *BB0(0.d0,MCH(i)**2,ST(1)**2,scR2)
     C		         +SST*(gg2*SST*V(i,1)+HTQ*CST*V(i,2))
     C                    *BB0(0.d0,MCH(i)**2,ST(2)**2,scR2))
     C          +(gg2*CST*V(i,1)-HTQ*SST*V(i,2))**2
     C                    *BB1(0.d0,MCH(i)**2,ST(1)**2,scR2)
     C          +(gg2*SST*V(i,1)+HTQ*CST*V(i,2))**2
     C                    *BB1(0.d0,MCH(i)**2,ST(2)**2,scR2)
     C          +HBQ/H2Q*MB0*U(i,2)**2
     C             *(CST**2*BB1(0.d0,MCH(i)**2,ST(1)**2,scR2)
     C                +SST**2*BB1(0.d0,MCH(i)**2,ST(2)**2,scR2))
     C      +VVu*(gg2**2*V(i,1)**2*BB1(0.d0,MCH(i)**2,MUL**2,scR2)
     C            +(MQU/H1Q)**2*V(i,2)**2
     C                           *BB1(0.d0,MCH(i)**2,MUR**2,scR2)
     C            -2.d0*gg2/H2Q*U(i,2)*V(i,1)*MCH(i)
     C                           *BB0(0.d0,MCH(i)**2,MUL**2,scR2)
     C            +HBQ/H2Q*MB0*U(i,2)**2
     C                           *BB1(0.d0,MCH(i)**2,MUL**2,scR2))
	enddo
	epsY13=epsY13/(32.d0*pi**2)

	epsY23=0.d0
	do i=1,2
	epsY23=epsY23-2.d0/H2Q*U(i,2)*MCH(i)
     C                 *(CST*(gg2*CST*V(i,1)-HTQ*SST*V(i,2))
     C                    *BB0(0.d0,MCH(i)**2,ST(1)**2,scR2)
     C		         +SST*(gg2*SST*V(i,1)+HTQ*CST*V(i,2))
     C                    *BB0(0.d0,MCH(i)**2,ST(2)**2,scR2))
     C          +(gg2*CST*V(i,1)-HTQ*SST*V(i,2))**2
     C                    *BB1(0.d0,MCH(i)**2,ST(1)**2,scR2)
     C          +(gg2*SST*V(i,1)+HTQ*CST*V(i,2))**2
     C                    *BB1(0.d0,MCH(i)**2,ST(2)**2,scR2)
     C          +HBQ/H2Q*MB0*U(i,2)**2
     C             *(CST**2*BB1(0.d0,MCH(i)**2,ST(1)**2,scR2)
     C                +SST**2*BB1(0.d0,MCH(i)**2,ST(2)**2,scR2))
     C      +VVc*(gg2**2*V(i,1)**2*BB1(0.d0,MCH(i)**2,MUL**2,scR2)
     C            +(MC0/H1Q)**2*V(i,2)**2
     C                           *BB1(0.d0,MCH(i)**2,MUR**2,scR2)
     C            -2.d0*gg2/H2Q*U(i,2)*V(i,1)*MCH(i)
     C                           *BB0(0.d0,MCH(i)**2,MUL**2,scR2)
     C            +HBQ/H2Q*MB0*U(i,2)**2
     C                           *BB1(0.d0,MCH(i)**2,MUL**2,scR2))
	enddo
	epsY23=epsY23/(32.d0*pi**2)
	
C	 e) epst0 (*tanb)
	epst0=epst3-Vtb2*epsY31


C        f) Couplings [X^s_RL]^JI as in eqs. (3.55) and (3.56), but
C         WITHOUT the S dependent mixing angles

	sigRLbs=MB0*epsY32/(H1Q*(1.d0+epst0)*(1.d0+epst3))
	
	sigRLbd=MB0*epsY31/(H1Q*(1.d0+epst0)*(1.d0+epst3))

	sigLRbs=MS0*epsY23/(H1Q*(1.d0+epst0)*(1.d0+epst3))
     C          *(1.d0+epst3+(epst2-epst3)*epsY32/epsY23)/(1.d0+epst2)

	sigLRbd=MD*epsY13/(H1Q*(1.d0+epst0)*(1.d0+epst3))
     C          *(1.d0+epst3+(epst1-epst3)*epsY31/epsY13)/(1.d0+epst1)


C**********************************************************************
C**********************************************************************

C     Towards DMs=m_Bs-m_Bbar_s and DMd=m_Bd-m_Bbar_d 
C      (still following [10])

C     Box contributions to both DMs and DMd
C          - Standard Model
	DSM=S0(xt)/xt
	
C	   - Charged Higgs
	DH=gg0(yt)/tanb**4+2.d0*xt*(ffp(xt,(CMASS/MW)**2)
     C                        +ggp(xt,xt,(CMASS/MW)**2)/4.d0)/tanb**2

C	   - Charginos
	do j=1,2
         do k=1,2
	 CCT(j,k)=V(j,1)*RST(k,1)
     C               -V(j,2)*RST(k,2)*MT0/(dsqrt(2.d0)*MW*sinb)
         enddo
        enddo
	
	aux=0.d0                       !3rd family contribution
	do i=1,2
	 do j=1,2
	  aux=aux+1.d0/MCH(j)**2
     C	     *(CCT(j,1)**2*CCT(i,1)**2
     C         *ggp(xstc(1,j),xstc(1,j),(MCH(i)/MCH(j))**2)
     C        +CCT(j,2)**2*CCT(i,2)**2
     C          *ggp(xstc(2,j),xstc(2,j),(MCH(i)/MCH(j))**2))
	  if(i.ne.j)aux=aux+1.d0/MCH(j)**2
     C        *CCT(j,1)*CCT(i,1)*CCT(i,1)*CCT(j,2)
     C         *(ggp(xstc(1,j),xstc(2,j),(MCH(i)/MCH(j))**2)
     C          +ggp(xstc(2,j),xstc(1,j),(MCH(i)/MCH(j))**2))
	 enddo
	enddo

	Dchid=MW**2/xt*aux
	Dchis=MW**2/xt*aux

	aux=0.d0                       !3rd family/1st-2nd interference
	do i=1,2
	 do j=1,2
	  if(i.ne.j)aux=aux+1.d0/MCH(j)**2
     C              *(V(j,1)*V(i,1)*CCT(i,1)*CCT(j,1)
     C                *ggp(xsqc(j),xstc(1,j),(MCH(i)/MCH(j))**2)
     C               +V(j,1)*V(i,1)*CCT(i,2)*CCT(j,2)
     C                *ggp(xsqc(j),xstc(2,j),(MCH(i)/MCH(j))**2)
     C               +V(i,1)*V(j,1)*CCT(j,1)*CCT(i,1)
     C                *ggp(xstc(1,j),xsqc(j),(MCH(i)/MCH(j))**2)
     C               +V(i,1)*V(j,1)*CCT(j,2)*CCT(i,2)
     C                *ggp(xstc(2,j),xsqc(j),(MCH(i)/MCH(j))**2))
	 enddo
	enddo

	Dchid=Dchid+VVu*MW**2/xt*aux
	Dchis=Dchis+VVc*MW**2/xt*aux

	aux=0.d0                       !1st-2nd family contribution
	do i=1,2
	 do j=1,2
	  aux=aux+1.d0/MCH(j)**2*V(j,1)**2*V(j,2)**2
     C              *ggp(xsqc(j),xsqc(j),(MCH(i)/MCH(j))**2)
	 enddo
	enddo

	Dchid=Dchid+VVu**2*MW**2/xt*aux
	Dchis=Dchis+VVc**2*MW**2/xt*aux

C**********************************************************************
	
C	Double Penguin Contributions to DMd
C	   - Wilson coefficients
	aux=0.d0
	do i=1,3
	aux=aux+(SCOMP(i,1)-SCOMP(i,2)*tanb)**2
     C *sgn(SMASS(i)**2-(5.2793d0)**2)/
     C dsqrt((SMASS(i)**2-(5.2793d0)**2)**2+(SMASS(i)*WIDTH(i))**2)
	enddo
	do i=1,2
	aux=aux+(PCOMP(i,1)*cosb+PCOMP(i,1)*sinb*tanb)**2
     C *sgn(PMASS(i)**2-(5.2793d0)**2)/
     C dsqrt((PMASS(i)**2-(5.2793d0)**2)**2+(PMASS(i)*WIDTH(3+i))**2)
	enddo

	C2LR=-(4.d0*pi/(GF*MW))**2*sigRLbd*sigLRbd*aux

	aux=0.d0
	do i=1,3
	aux=aux+(SCOMP(i,1)-SCOMP(i,2)*tanb)**2
     C *sgn(SMASS(i)**2-(5.2793d0)**2)/
     C dsqrt((SMASS(i)**2-(5.2793d0)**2)**2+(SMASS(i)*WIDTH(i))**2)
	enddo
	epst=aux
	do i=1,2
	aux=aux-(PCOMP(i,1)*cosb+PCOMP(i,1)*sinb*tanb)**2
     C *sgn(PMASS(i)**2-(5.2793d0)**2)/
     C dsqrt((PMASS(i)**2-(5.2793d0)**2)**2+(PMASS(i)*WIDTH(3+i))**2)
	enddo

	C1SLL=-(4.d0*pi/(GF*MW))**2*sigRLbd**2*aux/2.d0

	C1SRR=-(4.d0*pi/(GF*MW))**2*sigLRbd**2*aux/2.d0

	DDP=(0.90d0*C2LR-0.37d0*(C1SLL+C1SRR))/xt

C*********************************************************************	

C	Results for DMd
	etaB=0.551d0

	aux=(GF**2*MW**2/(6.d0*pi**2)*etaB*5.2794d0*(fBdsqBd)**2
     C         *xt*VtdVtb2/(6.58211915d-13))
     
	DMd=aux*abs(DSM+DH+Dchid+DDP)


C	Error Estimate: (2sigma on CKM and lattice QCD uncertainties)
C               * lattice QCD (sources [16,17]):
C                 1.134 < fBd sqrt(BBd) / fBs sqrt(BBs)=1.216 < 1.298
C                 0.239 GeV < fBs sqrt(BBs)=0.281 GeV < 0.323 GeV
C		* CKM factor (source [18]):
C                 6.2 10^-3 < VtbVtd=8.6 10^-3 < 11. 10^-3
C               * 30% on BSM (1st order QCD) contributions

C     First: 2 sigma "SM" (relative) error bars from CKM and 
C       lattice uncertainties: (2sigma, added quadratically)
	DMdMax=(1.d0+dsqrt((dVtdVtb2/VtdVtb2)**2
     C          +(2.d0*dfBdsqBd/fBdsqBd)**2))
	DMdmin=(1.d0-dsqrt((dVtdVtb2/VtdVtb2)**2
     C          +(2.d0*dfBdsqBd/fBdsqBd)**2))

C      Total error bars, allowing for 30% theory error 
C          on each BSM contribution:
	DMdmax=DMdMax*(DMd+aux*0.3d0*(dabs(DH)+dabs(Dchid)+dabs(DDP)))
	
	DMdmin=DMdmin*(DMd-aux*0.3d0*(dabs(DH)+dabs(Dchid)+dabs(DDP)))

C	Comparison with experimental data (source [12]):
C              (Recall: 2 sigma bounds: 0.499ps-1 < DMd < 0.515ps-1)

     	prob(34)=0.d0 

	IF(DMdmin.GE.DMdexpMax)
     C     PROB(34)=DMdmin-DMdexpMax     
     	IF(DMdmax.LE.DMdexpmin)
     C     PROB(34)=DMdmax-DMdexpmin

C*********************************************************************	
C*********************************************************************	

C	 Double Penguin Contributions to DMs
C	   - Wilson coefficients
	aux=0.d0
	do i=1,3
	aux=aux+(SCOMP(i,1)-SCOMP(i,2)*tanb)**2
     C *sgn(SMASS(i)**2-(5.3696d0)**2)/
     C dsqrt((SMASS(i)**2-(5.3696d0)**2)**2+(SMASS(i)*WIDTH(i))**2)
	enddo
	do i=1,2
	aux=aux+(PCOMP(i,1)*cosb+PCOMP(i,1)*sinb*tanb)**2
     C *sgn(PMASS(i)**2-(5.3696d0)**2)/
     C dsqrt((PMASS(i)**2-(5.3696d0)**2)**2+(PMASS(i)*WIDTH(3+i))**2)
	enddo

	C2LR=-(4.d0*pi/(GF*MW))**2*sigRLbs*sigLRbs*aux

	aux=0.d0
	do i=1,3
	aux=aux+(SCOMP(i,1)-SCOMP(i,2)*tanb)**2
     C *sgn(SMASS(i)**2-(5.3696d0)**2)/
     C dsqrt((SMASS(i)**2-(5.3696d0)**2)**2+(SMASS(i)*WIDTH(i))**2)
	enddo
	do i=1,2
	aux=aux-(PCOMP(i,1)*cosb+PCOMP(i,1)*sinb*tanb)**2
     C *sgn(PMASS(i)**2-(5.3696d0)**2)/
     C dsqrt((PMASS(i)**2-(5.3696d0)**2)**2+(PMASS(i)*WIDTH(3+i))**2)
	enddo

	C1SLL=-(4.d0*pi/(GF*MW))**2*sigRLbs**2*aux/2.d0

	C1SRR=-(4.d0*pi/(GF*MW))**2*sigLRbs**2*aux/2.d0	
	
	DDP=(0.90d0*C2LR-0.37d0*(C1SLL+C1SRR))/xt

C*********************************************************************	

C	Results for DMs	

	aux=GF**2*MW**2/(6.d0*pi**2)*0.55d0*5.3696d0*(fBssqBs)**2
     C         *VtsVtb2*xt/(6.58211915d-13)
	
	DMs=aux*dabs(DSM+DH+Dchis+DDP)

C	Error Estimate: (2sigma on CKM and lattice QCD uncertainties)
C	        * lattice QCD (source [16]):
C                 0.239 GeV < fBs sqrt(BBs)=0.281 GeV < 0.323 GeV
C		* CKM factor (source [18]):
C		  37.9 10^-3 < VtsVtb=41.3 10^-3 < 44.7 10^-3
C		* 30% on BSM (1st order QCD) contributions

C	First: 2 sigma "SM" (relative) error bars from CKM and 
C         lattice uncertainties: (2sigma, added quadratically)
	DMsMax=(1.d0+dsqrt((dVtsVtb2/VtsVtb2)**2
     C          +(2.d0*dfBssqBs/fBssqBs)**2))
	DMsmin=(1.d0-dsqrt((dVtsVtb2/VtsVtb2)**2
     C          +(2.d0*dfBssqBs/fBssqBs)**2))
	
C	Total error bars, allowing for 30% theory error 
C          on each BSM contribution:

	DMsmax=DMsMax*(DMs+aux*0.3d0*(dabs(DH)+dabs(Dchid)+dabs(DDP)))
	DMsmin=DMsmin*(DMs-aux*0.3d0*(dabs(DH)+dabs(Dchid)+dabs(DDP)))

C	Comparison with experimental data (source [14])
C               (Recall: 2 sigma bounds: 17.53ps-1 < DMs < 18.01ps-1)

	prob(33)=0.d0      

	IF(DMsmin.GE.DMsexpMax)
     C     PROB(33)=DMsmax-DMsexpMax
	IF(DMsmax.LE.DMsexpmin)
     C     PROB(33)=DMsmax-DMsexpmin

C*********************************************************************	
C*********************************************************************	

C	Contributions to BR(Bs -> mu+ mu-)    (see [10,11])
c	 a) Wilson coefficients - SM contribution
	Ca=1.d0/4.d0*(xt/(1.d0-xt)+xt/(xt-1.d0)**2*dlog(xt))
     C   -xt/8.d0*((xt-6.d0)/(xt-1.d0)
     C   +(3.d0*xt+2.d0)/(xt-1.d0)**2*dlog(xt))

C	 b) Wilson coefficients - DP contribution
	aux=0.d0
	do i=1,3
	aux=aux-(SCOMP(i,1)-SCOMP(i,2)*tanb)*SCOMP(i,2)
     C *sgn(SMASS(i)**2-(5.3696d0)**2)/(cosb
     C *dsqrt((SMASS(i)**2-(5.3696d0)**2)**2+(SMASS(i)*WIDTH(i))**2))
	enddo
	Cs=gg2/(2.d0*MW)*(pi/(GF*MW))**2*mmu/MB0*sigRLbs*aux

	aux=0.d0
	do i=1,2
	aux=aux-(-PCOMP(i,1)*cosb-PCOMP(i,1)*sinb*tanb)
     C          *PCOMP(i,1)*tanb*sgn(PMASS(i)**2-(5.3696d0)**2)/
     C dsqrt((PMASS(i)**2-(5.3696d0)**2)**2+(PMASS(i)*WIDTH(3+i))**2)
	enddo
	Cp=-gg2/(2.d0*MW)*(pi/(GF*MW))**2*mmu/MB0*sigRLbs*aux

C	 c) Branching ratio	
	sqBBs=dsqrt(1.3d0)       !([10])       sqrt(B_Bs)

	aux=GF**2*ALEM0**2/(s2tw**2*64.d0*pi**3)*(5.3696d0)**5
     C     *1.454d0/(6.58211915d-13)*(fBssqBs/sqBBs)**2
     C     *Vtsvtb2*dsqrt(1.d0-4.d0*mmu**2/5.3696d0**2)

	BRBMUMU=aux*
     C        ((1.d0-4.d0*mmu**2/(5.3696d0)**2)/(1.d0+MS/Mb)**2*Cs**2
     C         +(Cp/(1.d0+ms/mb)+2.d0*mmu/(5.3696d0)**2*Ca)**2)

C	Error Estimate: (2sigma on CKM and lattice QCD uncertainties)
C	        * lattice QCD (source [16]):
C                 0.239 GeV < fBs sqrt(BBs)=0.281 GeV < 0.323 GeV
C		* CKM factor (source [18]):
C		  37.9 10^-3 < VtsVtb=41.3 10^-3 < 44.7 10^-3
C		* 30% on BSM (1st order QCD) contributions

C       First: 2 sigma "SM" (relative) error bars from CKM and 
C       lattice uncertainties: (2sigma, added quadratically)
	BRBMUMUMAX=(1.d0+dsqrt((dVtsVtb2/VtsVtb2)**2
     C              +(2.d0*dfBssqBs/fBssqBs)**2))
	BRBMUMUmin=(1.d0-dsqrt((dVtsVtb2/VtsVtb2)**2
     C              +(2.d0*dfBssqBs/fBssqBs)**2))
	
C	Total error bars, allowing for 30% theory error 
C          on each BSM contribution:
	BRBMUMUmax=BRBMUMUMAX*(BRBMUMU+aux*0.3d0*
     C        ((1.d0-4.d0*mmu**2/(5.3696d0)**2)/(1.d0+MS/Mb)**2*Cs**2
     C         +(Cp/(1.d0+ms/mb)+2.d0*mmu/(5.3696d0)**2*Ca)**2))
	BRBMUMUmin=BRBMUMUmin*(BRBMUMU-aux*0.3d0*
     C        ((1.d0-4.d0*mmu**2/(5.3696d0)**2)/(1.d0+MS/Mb)**2*Cs**2
     C         +(Cp/(1.d0+ms/mb)+2.d0*mmu/(5.3696d0)**2*Ca)**2))

C	Comparison with experimental data (source [19])
C            (Recall 95% C.L.: BR(Bs->mu+mu-) < 5.8 10^-8)

	prob(35)=0.d0  

	IF(BRBMUMUmin.GE.BRBMUMUexpMax)
     C     PROB(35)=(BRBMUMUmin-BRBMUMUexpMax)*1.d8

C*********************************************************************	
C*********************************************************************	

C	Branching ratio BR(B+ -> tau+ nu_tau)
C                                                   following [12]
	tauB=1.638d-12/(6.58211915d-25)
	mBu=5.279d0

	rh=(1.d0-(mBu/CMASS)**2*tanb**2/(1.d0+epst0))**2
	aux=GF**2*mBu*MTAU**2/(8.d0*pi)*(1.d0-(MTAU/mBu)**2)**2
     C               *fB**2*Vub2*tauB
	BRBtaunu=aux*rh

C	Comparison with experimental data:
C          * 2 sigma bounds (source [12]): 
C                  0.34 10^-4 < BR(B+ -> tau+ nu) < 2.30 10^-4
C          * hadronic parameter (source [14]; 2sigma):
C                  0.172 GeV < fB=0.216 GeV <0.260 GeV
C          * CKM factor (2sigma): 3.3 10^-3 < Vub=4.0 10^-3 < 4.7 10^-3

C          2 sigma (absolute) error bars from CKM and 
C       lattice uncertainties:
     
	BRBtaunumax=(1.d0+dfB/fB)**2*(1.d0+dVub/dsqrt(Vub2))**2*BRBtaunu
	BRBtaunumin=(1.d0-dfB/fB)**2*(1.d0-dVub/dsqrt(Vub2))**2*BRBtaunu

	prob(36)=0.d0  

	IF(BRBtaunumin.GE.BRBTAUNUexpmax)
     C     PROB(36)=(BRBtaunumin-BRBTAUNUexpmax)/aux
   
	IF(BRBtaunumax.LE.BRBTAUNUexpmin)
     C     PROB(36)=(BRBtaunumax-BRBTAUNUexpmin)/aux
 
C*********************************************************************	
C*********************************************************************	
C  Contributions to BR(B->Xs gamma) at the Matching scale (sc0)
C**********************************************************************

C   Towards SUSY Contributions including large tan(beta) effects from 
C   charged Higgs couplings      (Used: [2], [3], [6])

C   First:  Effective Charged Higgs Yukawa Couplings
C    Epsilon_b, Epsilon'_b and Epsilon_t 

C      a) epsilon_b
         epsb=-asmsusy*2.d0/(3.d0*pi)*
     C         ((mu-AAB/tanb)/MGL*H2(xBg1,xBg2)
     C    +1.d0/(2.d0*tanb)*(1.d0-BB1(0.d0,MGL**2,SB(1)**2,QSTSB)
     C            -BB1(0.d0,MGL**2,SB(2)**2,QSTSB)))

         aux=0.d0
	 do j=1,2
	 aux=aux+U(j,2)*H2(xstc(1,j),xstc(2,j))*V(j,2)/MCH(j)
	 enddo
	 
         epsb=epsb-(HTQ/(4.d0*pi))**2*(AAT-mu/tanb)*aux      
	 
	 aux=(CST/ST(1))**2*H2(Mch2**2/ST(1)**2,mu**2/ST(1)**2)
     C       +(SST/ST(2))**2*H2(Mch2**2/ST(2)**2,mu**2/ST(2)**2)
     C       +(CSB/SB(1))**2*H2(Mch2**2/SB(1)**2,mu**2/SB(1)**2)/2.d0
     C       +(SSB/SB(2))**2*H2(Mch2**2/SB(2)**2,mu**2/SB(2)**2)/2.d0
	 
	 epsb=epsb+ALEMMZ/(4.d0*S2TW*pi)*mu*Mch2*aux	

C      b) epsilon'_b(t)
         epsbp=-asmsusy*2.d0/(3.d0*pi)*(mu-AAB/tanb)/MGL*
     C     (CST**2*(H2(xTg1,xBg2)*CSB**2+H2(xTg1,xBg1)*SSB**2)+
     C      SST**2*(H2(xTg2,xBg2)*CSB**2+H2(xTg2,xBg1)*SSB**2))      

         aux=0.d0
	 do j=1,5
	 aux=aux+1.d0/MNEU(j)*N(j,4)*N(j,3)*
     C     (CST**2*(H2(xTneu(2,j),xBneu(1,j))*CSB**2+
     C              H2(xTneu(2,j),xBneu(2,j))*SSB**2)+
     C      SST**2*(H2(xTneu(1,j),xBneu(1,j))*CSB**2+
     C              H2(xTneu(1,j),xBneu(2,j))*SSB**2))  
 	 enddo

	epsbp=epsbp+(HTQ/(4.d0*pi))**2*(AAT-mu/tanb)*aux  
	 
	aux=(CST/ST(1))**2*H2(Mch2**2/ST(1)**2,mu**2/ST(1)**2)/2.d0
     C       +(SST/ST(2))**2*H2(Mch2**2/ST(2)**2,mu**2/ST(2)**2)/2.d0
     C       +(CSB/SB(1))**2*H2(Mch2**2/SB(1)**2,mu**2/SB(1)**2)
     C       +(SSB/SB(2))**2*H2(Mch2**2/SB(2)**2,mu**2/SB(2)**2)
	 
	epsbp=epsbp+ALEMMZ/(4.d0*S2TW*pi)*mu*Mch2*aux

C      c) epsilon_t(s)
         epst= -asmsusy*2.d0/(3.d0*pi)*(mu+AAT/tanb)/MGL*
     C             (CST**2*H2(xTg2,xQg)+SST**2*H2(xTg1,xQg))
     
     	 aux=0.d0
	 do,i=1,5
	 aux=aux+N(i,4)*N(i,3)/MNEU(i)*(
     C    CST**2*CSB**2*H2(ST(1)**2/MNEU(i)**2,SB(2)**2/MNEU(i)**2)
     C   +CST**2*SSB**2*H2(ST(1)**2/MNEU(i)**2,SB(1)**2/MNEU(i)**2)
     C   +SST**2*CSB**2*H2(ST(2)**2/MNEU(i)**2,SB(2)**2/MNEU(i)**2)
     C   +SST**2*SSB**2*H2(ST(2)**2/MNEU(i)**2,SB(1)**2/MNEU(i)**2))
    
	 enddo

 	 epst=epst+HBQ**2/(16.d0*pi**2)*mu/tanb*aux
	 
C**********************************************************************
C	Chargino/Squark Contributions

C      1) Factors
         do j=1,2
	  CCD(J)=U(J,2)*MW/(dsqrt(2.d0)*COSB*MCH(J))
         do k=1,2
	  CCT(j,k)=V(j,1)*RST(k,1)
     C                    -V(j,2)*RST(k,2)*HTQ/dsqrt(g2)
         enddo
         enddo

C      2) Calculation of the Wilson Coefficients (SUSY Scale)
         c7CHARS=0.d0
         c8CHARS=0.d0
	 
	 akk=1.d0/(1.d0+epsb*tanb)
	 
         do j=1,2
	  C7CHARS=C7CHARS+2.d0/3.d0*V(J,1)**2*(MW/MUL)**2*FF1(XSQC(J))
     C         +akk*CCD(J)*V(J,1)*FF3(XSQC(J))

	  C8CHARS=C8CHARS+2.d0/3.d0*V(J,1)**2*(MW/MUL)**2*FG1(XSQC(J))
     C         +akk*CCD(J)*V(J,1)*FG3(XSQC(J))

         do k=1,2  ! k is the stop index
	   C7CHARS=C7CHARS-
     C       2.d0/3.d0*CCT(J,K)**2*(MW/ST(K))**2*FF1(XSTC(K,J))
     C       -akk*CCD(J)*CCT(J,K)*RST(K,1)*FF3(XSTC(K,J))
     
	   C8CHARS=C8CHARS-
     C       2.d0/3.d0*CCT(J,K)**2*(MW/ST(K))**2*FG1(XSTC(K,J))
     C       -akk*CCD(J)*CCT(J,K)*RST(K,1)*FG3(XSTC(K,J))
     
         enddo
         enddo
	 
C      3) Evolution from the SUSY scale to sc0

	 etaS= asmsusy/asc0

         C7CHAR=etaS**(16.d0/21.d0)*C7CHARS+
     C    8.d0/3.d0*(etaS**(14.d0/21.d0)-etaS**(16.d0/21.d0))*C8CHARS
         C8CHAR=etaS**(14.d0/21.d0)*C8CHARS
	 
C**********************************************************************
C         Charged Higgs Contributions
C      1) Lowest Order:

        c70HIG=au**2/3.d0*ff1(yt)-au*ad*ff2(yt)
        c80HIG=au**2/3.d0*fg1(yt)-au*ad*fg2(yt)

C      2) Order alpha_s:

         C41HIG=au**2*eh(yt)
         C71HIG=ffh(yt,tanb)-4.d0/9.d0*C41HIG
         C81HIG=fgh(yt,tanb)-1.d0/6.d0*C41HIG
	 
C      3) Large tan(beta) Corrections:

         dC7HIG=-tanb*(epsb+epst)*akk*ff2(yt)
         dC8HIG=-tanb*(epsb+epst)*akk*fg2(yt)
	 
C      4) Evolution from M_Higgs to sc0:

         etaH=asmh/asc0	

	 C7HIG=etaH**(16.d0/21.d0)*(C70HIG+dC7HIG)
     C       +8.d0/3.d0*(etaH**(2.d0/3.d0)-etaH**(16.d0/21.d0))
     C               *(C80HIG+dC8HIG)
	 C8HIG=etaH**(2.d0/3.d0)*(C80HIG+dC8HIG)

C**********************************************************************	
C       SM Contributions at the matching scale sc0 = m_top

C       1) Lowest Order: 
         c70SM=ff1(xt) 
         c80SM=fg1(xt) 

C       2) Order alpha_s correction to C41:

	 C41SM=esm(xt)-2.d0/3.d0+2.d0/3.d0*dlog(sc0**2/MW**2)
         
C	3) Large tan(beta) Corrections:

	 dC7SM=tanb*(epsb-epsbp)*akk*ff2(xt)
         dC8SM=tanb*(epsb-epsbp)*akk*fg2(xt)

C**********************************************************************	
C	Neutral Higgs Contribution to b -> s gamma 
C        (following [10], but including the evolution from the scale
C          PMASS(1) to the matching scale sc0)

	aux=0.d0
	do i=1,3
	aux=aux+(SCOMP(i,1)-SCOMP(i,2)*tanb)
     C   *(SCOMP(i,2)+SCOMP(i,1)*epst3/tanb)/SMASS(i)**2
	enddo
	do i=1,2
	aux=aux-(PCOMP(i,1)*cosb+PCOMP(i,1)*sinb*tanb)
     C       *(PCOMP(i,1)*sinb-PCOMP(i,1)*cosb*epst3/tanb)/PMASS(i)**2
	enddo
	
	C70S0=1.d0/18.d0*MW**2/g2*MB0/(H2Q*(1.d0+epst3))*sigRLbs*aux	

c	 - Running to mt
	IF(PMASS(1).le.sc0)THEN
	 IF(PMASS(1).ge.mb)THEN
	  eta0=asc0/asf(PMASS(1))
	 ELSE
	  eta0=asmb/asc0
	 ENDIF
	ELSE
 	 eta0=1.d0
	ENDIF
	C80S0=eta0**(14.d0/23.d0)*C70S0
	C70S0=eta0**(16.d0/23.d0)*C70S0+
     C     8.d0/3.d0*(eta0**(14.d0/23.d0)-eta0**(16.d0/23.d0))*C70S0

C**********************************************************************	
C       Sums at the matching scale sc0 = m_top
	     
	 C70BSM=dC7SM+C7HIG+C7CHAR+C70S0
	 C80BSM=dC8SM+C8HIG+C8CHAR+C80S0
	     	    
C     DC70BSM and DC80BSM are (conservative) error estimates:
C Dominant errors: HIG: Order (alphas**2) ~ 10%, 
C                  CHAR: Order(alphas) ~ 30%
	     	     
	 DC70BSM=.1d0*DABS(C7HIG)+.3d0*DABS(C7CHAR)+.3d0*dabs(C70S0)
	 DC80BSM=.1d0*DABS(C8HIG)+.3d0*DABS(C8CHAR)+.3d0*dabs(C80S0)

         C70=C70SM+C70BSM
	 C80=C80SM+C80BSM
	     
C**********************************************************************
C**********************************************************************
C      Coefficients for BR(B->Xs gamma) at the Scale m_b
C         Used: [1], [4], [5], [7], [8]
C**********************************************************************
C
C        Parameters:

	eta=asc0/asmb

C	A: LO COEFFICIENTS C10b, C20b, C70b, C80b
C          (Used in Bremsstrahlung and HQET Corrections)

C	    1) "Magic Numbers"
C	 a) Array: aa_i
	aa(1)=14.d0/23.d0
        aa(2)=16.d0/23.d0
        aa(3)=6.d0/23.d0
        aa(4)=-12.d0/23.d0
        aa(5)=0.4086d0
        aa(6)=-0.4230d0
        aa(7)=-0.8994d0
        aa(8)=0.1456d0

C	 b) Array: bb_i
	do i=1,4
	 j=4+i
         bb(i)=aa(j)
        enddo

C	 c) Array: hh_i
	hh(1)=626126.d0/272277.d0
        hh(2)=-56281.d0/51730.d0
        hh(3)=-3.d0/7.d0
        hh(4)=-1.d0/14.d0
        hh(5)=-0.6494d0
        hh(6)=-0.0380d0
        hh(7)=-0.0185d0
        hh(8)=-0.0057d0

C	 d) Array: h8_i
	h8(1)=-0.9135d0
	h8(2)=0.0873d0
	h8(3)=-0.0571d0
	h8(4)=0.0209d0
	
C	    2) Wilson Coefficients
	
	C10b=eta**(aa(3))-eta**(aa(4))

	C20b=1.d0/3.d0*(2.d0*eta**(aa(3))+eta**(aa(4)))
     
     	C70b=0.d0
	do i=1,8
	 C70b=C70b+hh(i)*eta**(aa(i))
	enddo
	C70b=C70b+eta**(16.d0/23.d0)*C70+
     C   8.d0/3.d0*(eta**(14.d0/23.d0)-eta**(16.d0/23.d0))*C80

	C80b=0.d0
	do i=1,4
	 C80b=C80b+h8(i)*eta**(bb(i))
	enddo
	C80b=C80b+eta**(14.d0/23.d0)*(C80+313063.d0/363036.d0)

C**********************************************************************
C	B: Lowest Order c-Quark (KC0) and t-Quark contributions (KT0)
C        (Only for Bremsstrahlung and HQET Corrections)

	KC0=0.d0
	
	do i=1,8
	 KC0=KC0+hh(i)*eta**(aa(i))
	enddo

	KC0=KC0-23.d0/36.d0*eta**(16.d0/23.d0)
     C     -8.d0/9.d0*(eta**(14.d0/23.d0)-eta**(16.d0/23.d0))
     
     	KT0=eta**(4.d0/23.d0)*(C70SM+23.d0/36.d0)
     C   +8.d0/3.d0*(eta**(2.d0/23.d0)-eta**(4.d0/23.d0))
     C                             *(C80SM+1.d0/3.d0)
     	
	KBSM0=C70BSM*ETA**(4.d0/23.d0)+
     C    C80BSM*(ETA**(2.d0/23.d0)-ETA**(4.d0/23.d0))*8.d0/3.d0
		
C**********************************************************************
C	C: Towards KC and KT including NLO in alpha_s

C	    1) "Magic Numbers" needed for NLO
C	 a) Array: e_i
	ee(1)=5.2620d0
        ee(2)=-8516.d0/2217.d0
        ee(3)=0.d0
        ee(4)=0.d0
        ee(5)=-1.9043d0
        ee(6)=-0.1008d0
        ee(7)=0.1216d0
        ee(8)=0.0183d0
	
C	 b) Array: dd_i
	do i=1,8
	 dd(i)=hh(i)
	enddo
	dd(1)=dd(1)-8.d0/9.d0
	dd(2)=dd(2)+8.d0/9.d0-23.d0/36.d0
	
C	 c) Array: dt_i
	dt(1)=-17.6507d0
	dt(2)=11.3460d0
	dt(3)=2.4692d0
	dt(4)=-0.8056d0
	dt(5)=4.8898d0
	dt(6)=-0.2308d0
	dt(7)=-0.5290d0
	dt(8)=0.1994d0

C	 d) Array: dte_i
	dte(1)=9.2746d0
	dte(2)=-6.9366d0
	dte(3)=-0.8740d0
	dte(4)=0.4218d0
	dte(5)=-2.7231d0
	dte(6)=0.4083d0
	dte(7)=0.1465d0
	dte(8)=0.0205d0
	
C	 e) Array: dim_i
	dim(1)=0.4702d0
	dim(2)=0.d0
	dim(3)=-0.4268d0
	dim(4)=-0.2222d0
	dim(5)=-0.9042d0
	dim(6)=0.1150d0
	dim(7)=-0.0975d0
	dim(8)=0.0115d0
	
C	 f) Array: da_i
	da(1)=0.d0
	da(2)=0.d0
	da(3)=0.8571d0
	da(4)=0.6667d0
	da(5)=0.1298d0
	da(6)=0.1951d0
	da(7)=0.1236d0
	da(8)=0.0276d0
	
C	 g) Array: db_i
	db(1)=0.d0
	db(2)=0.d0
	db(3)=0.8571d0
	db(4)=0.6667d0
	db(5)=0.2637d0
	db(6)=0.2906d0
	db(7)=-0.0611d0
	db(8)=-0.0171d0
	
C	    2) Charm Quark Contribution KC, LO + NLO (SM only)

C    For z=MC/MB, we choose a value that reproduces the NNLO result 
C    BRSG ~ 3.15 10^(-4) [8] in the SM:
	z=(0.307d0)**2  

	KC=0.d0
	do i=1,8
	KC=KC+eta**(aa(i))*(dd(i)*(1.d0+46.d0/3.d0*aa(i)
     C            *asmb/(4.d0*pi)*eta*dlog(sc0/MW))
     C            +asmb/(4.d0*pi)*(dt(i)+dte(i)*eta
     C            +(1.d0+VVtu)*(da(i)*af(z)+db(i)*bf(z))
     C            -VVtuim*(da(i)*afim(z)+db(i)*bfim(z))))
	enddo

C 	   3) Top Quark Contribution KT, LO + NLO (SM only)

	At1=A1(xt)
	Ft1=F1(xt)
	dAto=dA0(xt)
	dFto=dF0(xt)
	
	aux=0.d0
	do i=1,8
	 aux=aux+ee(i)*eta**(aa(i)+11.d0/23.d0)
	enddo
	
	KT=(1.d0-2.d0/9.d0*asmb**2)*
     C      (eta**(4.d0/23.d0)*(C70SM+23.d0/36.d0)
     C       +8.d0/3.d0*(eta**(2.d0/23.d0)-eta**(4.d0/23.d0))
     C                                     *(C80SM+1.d0/3.d0))
	KT=KT+asmb/(4.d0*pi)*(aux*(C41SM+7.d0/9.d0)
     C      +eta**(4.d0/23.d0)*(eta*(-1.d0/2.d0*At1)
     C      -2.d0*(12523.d0/3174.d0-7411.d0/4761.d0*eta
     C                   -2.d0/9.d0*pi**2)
     C      *(C70SM+23.d0/36.d0)-8.d0/3.d0*eta*(-1.d0/2.d0*Ft1)
     C      -2.d0*(-50092.d0/4761.d0+1110842.d0/357075.d0*eta
     C  +16.d0/27.d0*pi**2)*(C80SM+1.d0/3.d0))
     C      +eta**(2.d0/23.d0)*(8.d0/3.d0*eta*(-1.d0/2.d0*Ft1)
     C      -2.d0*(2745458.d0/357075.d0-38890.d0/14283.d0*eta
     C -4.d0/9.d0*pi*(pi))*(C80SM+1.d0/3.d0)))
	KT=KT+asc0/pi*dlog(sc0/MT)*4*xt*
     C          (-1.d0/2.d0*eta**(4.d0/23.d0)*dAto
     C           +4.d0/3.d0*(eta**(4.d0/23.d0)-eta**(2.d0/23.d0))*dFto)
     C       +asmb/(4.d0*pi)*((8.d0/3.d0*(C70SM+23.d0/36.d0)
     C             -64.d0/9.d0*(C80SM+1.d0/3.d0))*eta**(27.d0/23.d0)
     C             +32.d0/9.d0*eta**(25.d0/23.d0)*(C80SM+1.d0/3.d0))
     C                                       *dlog(sc0/MT)
     
C	 4) BSM Contribution to KT, LO + NLO 

	KBSM=(1.d0-2.d0/9.d0*asmb**2)*(eta**(4.d0/23.d0)*C70BSM
     C   +8.d0/3.d0*(eta**(2.d0/23.d0)-eta**(4.d0/23.d0))*C80BSM)
     
	KBSM=KBSM+asmb/(4.d0*pi)*(aux*C41HIG
     C      +eta**(4.d0/23.d0)*(eta*C71HIG
     C      -2.d0*(12523.d0/3174.d0-7411.d0/4761.d0*eta
     C             -2.d0/9.d0*pi**2)
     C      *C70BSM-8.d0/3.d0*eta*C81HIG
     C      -2.d0*(-50092.d0/4761.d0+1110842.d0/357075.d0*eta
     C    +16.d0/27.d0*pi**2)*C80BSM)
     C      +eta**(2.d0/23.d0)*(8.d0/3.d0*eta*C81HIG
     C      -2.d0*(2745458.d0/357075.d0-38890.d0/14283.d0*eta
     C    -4.d0/9.d0*pi*(pi))*C80BSM))
    
C           BSM error estimate from DC70BSM and DC80BSM:
	DKBSM=(1.d0-2.d0/9.d0*asmb**2)*(eta**(4.d0/23.d0)*DC70BSM
     C   +8.d0/3.d0*(eta**(2.d0/23.d0)-eta**(4.d0/23.d0))*DC80BSM)

C	 5) Imaginary Parts of KC, KT, KBSM

	KCIM=0.d0
	do i=1,8
	 KCIM=KCIM+eta**(aa(i))*(dim(i)*pi
     C                   +(1.d0+VVtu)*(da(i)*afim(z)+db(i)*bfim(z))
     C                   +VVtuim*(da(i)*af(z)+db(i)*bf(z)))
	enddo
	KCIM=asmb/(4.d0*pi)*KCIM

	KTIM=2.d0*asmb/9.d0*eta**(2.d0/23.d0)*(C80SM+1.d0/3.d0)

	KBSMIM=2.d0*asmb/9.d0*eta**(2.d0/23.d0)*C80BSM

**********************************************************************
C**********************************************************************
C    Towards BR(b -> s gamma) including Electroweal corrections,
C      Bremsstrahlung Corrections and Heavy Quark Effective Theory
C      Corrections  
C                                                   Used: [4], [7]
C**********************************************************************

C	   MB1S = "1S b-Quark Mass":
	MB1S=4.68d0   
	
C         Ratio rmu = m_b(sc0)/MB1S:

	rmu=0.578d0*0.1185d0/alsmz*(MB1S/4.69d0)**0.23d0
     C *(MC*(1.d0+asmc/pi*(-4.d0/3.d0))/1.25d0
     C  )**(-0.003d0)*(sc0/165.d0)**(-0.08d0)*(MB/4.69d0)**(0.006d0)

C**********************************************************************
C	   A: Summing the LO and NLO Contributions

	Ktot=KC+rmu*(KT+KBSM)
	KIM=KCIM+rmu*(KTIM+KBSMIM)

C**********************************************************************
C          B: ELECTROWEAK Corrections from
C   P.Gambino, U.Haisch, 
C   'Complete Electroweak Matching for Radiative B Decays'
C   JHEP 0110:020,2001, e-Print: hep-ph/0109058

	C7EMb=(88.d0/575.d0*eta**(16.d0/23.d0)
     C                      -40.d0/69.d0*eta**(-7.d0/23.d0)
     C                      +32.d0/75.d0*eta**(-9.d0/23.d0))
     C           *(C70BSM+asc0/(4.d0*pi)*C71HIG)
     C                     +(640.d0/1449.d0*eta**(14.d0/23.d0)
     C                     -704.d0/1725.d0*eta**(16.d0/23.d0)
     C                     +32.d0/1449.d0*eta**(-7.d0/23.d0)
     C                     -32.d0/575.d0*eta**(-9.d0/23.d0))
     C           *(C80BSM+asc0/(4.d0*pi)*C81HIG)
C	Summation of the Electroweak Corrections, the SM Contribution
C          is put to 0.0071

*	EPSew=0.0071d0
*     C     +ALEMMZ/(4.d0*pi)*(C7EMb-4.d0*rmu*KBSM0*dlog(MZ/MB))
 
	EPSew=0.0071d0
     C     +ALEMMZ*(C7EMb/asmb-rmu*KBSM0*dlog(MZ/MB)/pi)
    

C**********************************************************************
C	    C: Bremsstrahlung Contribution with delta=0.352
C For the cut on the photonic energy, we take E_ph > 1.6 GeV, that is
C delt=0.352. As for the f_ij coefficients, we will use numerical 
C values for ff_22 and ff_27.
C Formulae for delt=0.9 will be kept as comments.

	delt=0.352d0

C	 a) Coefficients: f_ij
	lndelt=dlog(delt)
	lndeltp=dlog(1.d0-delt)
	delt2=delt**2

C   If delt=0.9 we would have
C	ff22=0.107636d0-0.208484d0*dsqrt(z)-0.156146d0*z 
C   Now: delt=0.325:
	ff22=-0.01445908619d0+.9889254405d0*dsqrt(z)-6.618735072d0*z
     C       +17.52172828d0*z*dsqrt(z)-6.886559723d0*z**2
     C       -60.56748284d0*z**2*dsqrt(z)+89.45278858d0*z**3 
	ff22=ff22*((1.d0+VVtu)**2+VVtuim**2)+(VVtu**2+VVtuim**2)
     C   *(1.d0/3.d0*(1.d0-(1.d0-delt)**3)+1.d0/2.d0*(1.d0-delt)**2)
	ff11=1.d0/36.d0*ff22
	ff12=-1.d0/3.d0*ff22
	
C   If delt=0.9 we would have
C	ff27=-0.190805d0+0.948865d0*dsqrt(z)-0.787805d0*z 
C   Now: delt=0.325:
	ff27=-0.04047347658d0-.4349867931d0*dsqrt(z)+2.997393353d0*z
     C       -2.340077892d0*z*dsqrt(z)-4.025847473d0*z**2
     C       +7.412298665d0*z**2*dsqrt(z)-13.15396549d0*z**3 
	ff27=(1.d0+VVtu)*ff27
	ff17=-1.d0/6.d0*ff27
	ff28=-1.d0/3.d0*ff27
	ff18=-1.d0/6.d0*ff28
	ff77=1.d0/3.d0*(10.d0*delt+delt2-2.d0/3.d0*delt**3
     C             +delt*(delt-4.d0)*lndelt)
     C           -1.d0/3.d0*(2.d0*lndelt**2+7.d0*lndelt+31.d0/3.d0)
	ff78=8.d0/9.d0*(sp2(1.d0-delt)-pi**2/6.d0-delt*lndelt
     C        +9.d0/4.d0*delt-delt2/4.d0+delt**3/12.d0)
	ff88=1.d0/27.d0*(4.d0*sp2(1.d0-delt)-2.d0/3.d0*pi**2
     C       +8.d0*lndeltp-delt*(2.d0+delt)*lndelt
     C       +7.d0*delt+3.d0*delt2-2.d0/3.d0*delt**3
     C       -2.d0*(2.d0*delt+delt2+4.d0*lndeltp)*dlog(MB/0.95d0))
     
C	 b) Summing the Bremsstrahlung Contributions
	BREMS=asmb/pi*(ff11*C10b**2+ff12*C10b*C20b+ff17*C10b*C70b
     C                +ff18*C10b*C80b+ff22*C20b**2+ff27*C20b*C70b
     C                +ff28*C20b*C80b+ff77*C70b**2+ff78*C70b*C80b
     C                +ff88*C80b**2)

C**********************************************************************
C	    D: Heavy Quark Effective Theory Corrections

	lambd2=0.12d0  ! = 1/2(M_B*^2-M_B^2)
	aux=MC*(1.d0+asmc/pi*(-4.d0/3.d0))
	HQET=-lambd2/(9.d0*aux**2)*(C20b-1.d0/6.d0*C10b)
     C                                      *(KC0+rmu*(KT0+KBSM0))

C**********************************************************************
C	  E: Finally: The Branching Ratio B -> Xs gamma

C	     0) Parameters
	BRSL=0.1061d0

	CCSL=0.580d0
	
C	     1) Coefficient
	aux=6.d0*ALEM0/(pi*CCSL)*Vbsg*BRSL
	
C	     2) Result: BR(b -> Xs gamma)
	BRSG=aux*((Ktot+EPSew)**2+(KIM)**2+BREMS+HQET)
	
C        - Error estimate:
C     First: Variations of BRSG from BSM uncertainties DKBSM:        
	BRSG1=aux*((Ktot+DKBSM+EPSew)**2+(KIM)**2+BREMS+HQET)
	BRSG2=aux*((Ktot-DKBSM+EPSew)**2+(KIM)**2+BREMS+HQET)
C     Second: Add the SM uncertainties + 0.23 / - 0.43 (from [8,9]):
	DBRSGmax=(max(BRSG1,BRSG2)-BRSG)+0.23D-4
	DBRSGmin=(min(BRSG1,BRSG2)-BRSG)-0.43D-4

C     Comparison with experimental data:
C               * 2 sigma bounds (source [13])
C           3.07 10^-4 < BR(B -> Xs gamma)=3.55 10^-4 < 4.07 10^-4

	prob(32)=0.d0  

	IF(BRSG+DBRSGmin.GE.BRSGexpMax)
     C     PROB(32)=(BRSG+DBRSGmin-BRSGexpMax)*1.D4
	IF(BRSG+DBRSGmax.LE.BRSGexpmin)
     C     PROB(32)=(BRSG+DBRSGmax-BRSGexpmin)*1.D4
     
        BRSGmax= BRSG+DBRSGmax
	BRSGmin= BRSG+DBRSGmin

C**********************************************************************
C	  F: Constraints from BR(B-->X_s mu+ mu-)

C	     1) Wilson Coefficients
	C70=eta**(16/23)*(C70+58.d0/135.d0*(eta**(-10/23)-1.d0)
     C                    +29.d0/189.d0*(eta**(-28/23)-1))
	C70SM=eta**(16/23)*(C70SM
     C                      +58.d0/135.d0*(eta**(-10/23)-1.d0)
     C                      +29.d0/189.d0*(eta**(-28/23)-1.d0))

	aux=0.d0
	do i=1,3
	do j=1,3
	aux=aux+(SCOMP(i,1)-SCOMP(i,2)*tanb)*SCOMP(i,2)/cosb
     C          *(SCOMP(j,1)-SCOMP(j,2)*tanb)*SCOMP(j,2)/cosb
     C       *Intpropa(SMASS(i),WIDTH(i),SMASS(j),WIDTH(j),
     C                               1.d0/mb1s**2,6.d0/mb1s**2)
	enddo
	enddo
	CQ12=(gg2/(2.d0*MW)*(pi/(GF*MW))**2*mmu/MB0*sigRLbs)**2
	CQ12=(5.3696d0/s2tw)**2*CQ12*aux
	aux=0.d0
	do i=1,2
	do j=1,2
	aux=aux
     C      +(-PCOMP(i,1)*cosb-PCOMP(i,1)*sinb*tanb)*PCOMP(i,1)
     C      *(-PCOMP(j,1)*cosb-PCOMP(j,1)*sinb*tanb)*PCOMP(j,1)
     C       *tanb**2
     C      *Intpropa(PMASS(i),WIDTH(3+i),PMASS(j)
     C                    ,WIDTH(3+j),1.d0/mb1s**2,6.d0/mb1s**2)
	enddo
	enddo
	CQ22=(gg2/(2.d0*MW)*(pi/(GF*MW))**2*mmu/MB0*sigRLbs)**2
	CQ22=(5.3696d0/s2tw)**2*CQ22*aux

C	     2) Branching Ratio between 1GeV^2 and 6GeV^2
	BRBSll=1.59d-6+BRSL*4.d0/CCSL*Vbsg*(ALEM0/(4.d0*pi))**2
     C   *((8.d0*dlog(6.d0)-15.d0/mb1s**2+215.d0/(3.d0*mb1s**6))
     C                            *(C70**2-C70SM**2)
     C    +3.d0/2.d0*(cQ12+cQ22))

	BRBSllmax=(BRBSll-1.59d-6)*1.6d0+(1.59+.22d0)*1.d-6
	BRBSllmin=(BRBSll-1.59d-6)*.4d0+(1.59-.22d0)*1.d-6

C	     3) Comparison with experiment
C               * 2 sigma bounds
C           0.60 10^-6 < BR(B -> Xs l+l-)=1.60 10^-6 < 2.60 10^-6
	prob(40)=0.d0
	IF(BRBSllmin.GE.2.6d-6)
     C     PROB(40)=(BRBSllmin-2.6d-6)*1.D6
	IF(BRBSllmax.LE.0.6d-6)
     C     PROB(40)=(BRBSllmax-0.6d-6)*1.D6

C	     4) High M_{l+l-} region: 14.4 GeV^2 < M_{l+l-}^2 < mb^2
	aux=0.d0
	do i=1,3
	do j=1,3
	aux=aux+(SCOMP(i,1)-SCOMP(i,2)*tanb)*SCOMP(i,2)/cosb
     C          *(SCOMP(j,1)-SCOMP(j,2)*tanb)*SCOMP(j,2)/cosb
     C   *Intpropa(SMASS(i),WIDTH(i),SMASS(j),WIDTH(j),14.4d0/mb1s**2
     C                                    ,1.d0)
	enddo
	enddo
	CQ12=(gg2/(2.d0*MW)*(pi/(GF*MW))**2*mmu/MB0*sigRLbs)**2
	CQ12=(5.3696d0/s2tw)**2*CQ12*aux
	
	aux=0.d0
	do i=1,2
	do j=1,2
	aux=aux
     C      +(-PCOMP(i,1)*cosb-PCOMP(i,1)*sinb*tanb)*PCOMP(i,1)
     C      *(-PCOMP(j,1)*cosb-PCOMP(j,1)*sinb*tanb)*PCOMP(j,1)
     C       *tanb**2
     C*Intpropa(PMASS(i),WIDTH(3+i),PMASS(j),WIDTH(3+j),14.4d0/mb1s**2
     C                                         ,1.d0)
	enddo
	enddo
	CQ22=(gg2/(2.d0*MW)*(pi/(GF*MW))**2*mmu/MB0*sigRLbs)**2
	CQ22=(5.3696d0/s2tw)**2*CQ22*aux

	BRBSll=2.40d-7+BRSL*4.d0/CCSL*Vbsg*(ALEM0/(4.d0*pi))**2
     C   *((8.d0*dlog(mb1s**2/14.4d0)-4.d0*(1.d0-14.4d0/mb1s**2)
     C    +4.d0/3.d0*(1.d0-(14.4d0/mb1s**2)**3))*(C70**2-C70SM**2)
     C    +3.d0/2.d0*(cQ12+cQ22))

	BRBSllmax=(BRBSll-2.40d-7)*1.6d0+2.40d-7*(1.d0+2.d0*.29d0)
	BRBSllmin=(BRBSll-2.40d-7)*.4d0+2.40d-7*(1.d0-2.d0*.26d0)
c	prob(40)=0.d0
	IF(PROB(40).eq.0.d0)then
	IF(BRBSllmin.GE.6.8d-7)
     C     PROB(40)=(BRBSllmin-6.8d-7)*1.D7
	IF(BRBSllmax.LE.2.0d-7)
     C     PROB(40)=(BRBSllmax-2.0d-7)*1.D7
	ENDIF

C For comparison: the Pole contribution only (from Hiller):
	BRBSLL=BRSL/CCSL*Vbsg*3.d0*pi/(2.d0*dsqrt(2.d0)*GF)
     C      *mmu**2*PMASS(1)/(mb1s**8*WIDTH(4))
     C        *(mb1s**2-PMASS(1)**2)**2
     C   *sigRLbs**2*(-PCOMP(1,1)*cosb-PCOMP(1,1)*sinb*tanb)**2
     C   *(PCOMP(1,1)*tanb)**2*(5.3696d0/MB0)**2
	
        return
        END
	
	
C*********************************************************************	
C*********************************************************************
C*********************************************************************

C MSbar 5/6-flavour evolution for the Strong Coupling Constant
        double precision function asf(x)
	
C       ALPHA_S (x) - NB it uses 5 flavors for x<mt and 6 if x>=mt
	implicit none
        double precision x,pi,asc,fn,b0,b1,vvv,b0t,b1t,ast
	DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
	DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
	COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
	COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
	pi=4.d0*datan(1.d0)
	asc=ALSMZ
	
        fn=5.d0
        b0=11.d0-2.d0*fn/3.d0
        b1=102.d0-38.d0*fn/3.d0
        vvv=1.d0-b0*asc/(2.d0*pi)*dlog(MZ/x)
        asf=asc/vvv*(1.d0-b1/b0*asc/(4.d0*pi*vvv)*dlog(vvv))

        if(x.gt.MT)then
         vvv=1.d0-b0*asc/(2.d0*pi)*dlog(MZ/MT)
         ast=asc/vvv*(1.d0-b1/b0*asc/(4.d0*pi*vvv)*dlog(vvv))
         b0t=b0-2.d0/3.d0
         b1t=b1-38.d0/3.d0
         vvv=1.d0-b0t*ast/(2.d0*pi)*dlog(MT/x)
         asf=ast/vvv*(1.d0-b1t/b0t*ast/(4.d0*pi*vvv)*dlog(vvv))
        endif
        return
        end
		
C***************************************************************
C Dilogarithm
        double precision function sp2(x)
C       Li_2(x)
	implicit none
        double precision x,pi,f,bla
        external f
        pi=4.d0*datan(1.d0)
	bla=0.d0
        if(x.ge.-1.d0.and.x.le.0.5d0)then
         bla=f(x)
        else if(x.gt.0.5d0.and.x.lt.1.d0)then
         bla=-f(1.d0-x)+pi**2/6.d0-dlog(x)*dlog(1.d0-x)
        else if(x.lt.-1.d0)then
         bla=-f(1.d0/x)-pi**2/6.d0-.5d0*(dlog(-x))**2
C	else if(x.gt.-500.d0.and.x.lt.-10.d0)then
C	 bla=1.50616d0+0.153765d0*x-0.0000484249d0*x**2
C     C        -2.69934d-8*x**3
C     C        -1.97807d0*dlog(dabs(x))-0.0245271d0*x*dlog(x)
        else if(x.ge.1.d0)then
	bla=0.d0
*         write(6,*)'error in dilog',x
        endif
	sp2=bla
        return
        end

C***************************************************************

        double precision function sgn(x)
	
	implicit none
	double precision x
        if(x.ge.0.d0)then
	sgn=1.d0
	else
	sgn=-1.d0
	endif
        end
	
C***************************************************************

        double precision function f(x)
	
	implicit none
        integer i
	double precision b(12),x,z,cCc,sum
        z=-dlog(1.d0-x)
        b(1)=-.5d0
        b(2)=1.d0/6.d0
        b(3)=0.d0
        b(4)=-1.d0/30.d0
        b(5)=0.d0
        b(6)=1.d0/42.d0
        b(7)=0.d0
        b(8)=-1.d0/30.d0
        b(9)=0.d0
        b(10)=5.d0/66.d0
        b(11)=0.d0
        b(12)=-691.d0/2730.d0
        cCc=z
        sum=z   
        do i=1,12
        cCc=cCc*z/(i+1.d0)
        sum=sum+b(i)*cCc
        enddo
        f=sum
        end
	
C*******************************************************************

	double precision function ff1(x)

	implicit none
	double precision x,d,dg
	if(dabs(x-1.d0).gt.1.d-3)then
	d=1.d0/(x-1.d0)**3
	dg=dlog(x)/(x-1.d0)**4
	ff1=x*(7.d0-5.d0*x-8.d0*x**2)*d/24.d0
     C        +x**2*(3.d0*x-2.d0)*dg/4.d0
	else
	ff1=-5.d0/48.d0
	endif
	return
	end

C********************************************************************

	double precision function fg1(x)

	implicit none
	double precision x,d,dg
	if(dabs(x-1.d0).gt.1.d-2)then
	d=1.d0/(x-1.d0)**3
	dg=dlog(x)/(x-1.d0)**4
	fg1=x*(2.d0+5.d0*x-x**2)*d/8.d0-3.d0*x**2*dg/4.d0
	else
	fg1=-1.d0/16.d0
	endif
	return
	end

C********************************************************************

	double precision function ff2(x)

	implicit none
	double precision x,d,dg
	if(dabs(x-1.d0).gt.1.d-2)then
	d=1.d0/(x-1.d0)**2
	dg=dlog(x)/(x-1.d0)**3
	ff2=x*(3.d0-5.d0*x)*d/12.d0+x*(3.d0*x-2.d0)*dg/6.d0
	else
	ff2=-7.d0/36.d0
	endif
	return
	end

C********************************************************************

	double precision function fg2(x)

	implicit none
	double precision x,d,dg
	if(dabs(x-1.d0).gt.1.d-2)then
	d=1.d0/(x-1.d0)**2
	dg=dlog(x)/(x-1.d0)**3
	fg2=x*(3.d0-x)*d/4.d0-x*dg/2.d0
	else
	fg2=-1.d0/6.d0
	endif
	return
	end

C********************************************************************

	double precision function ff3(x)

	implicit none
	double precision x,ff1,ff2
	ff3=2.d0/3.d0*(1.d0-1.d0/x)*ff1(x)+ff2(x)+23.d0/36.d0
	return
	end
C********************************************************************

	double precision function fg3(x)

	implicit none
	double precision x,fg1,fg2
	fg3=2.d0/3.d0*(1.d0-1.d0/x)*fg1(x)+fg2(x)+1.d0/3.d0
	return
	end

C********************************************************************

        double precision function esm(x)
	
	implicit none
        double precision x,d,dg
        if(dabs(x-1.d0).gt.1.d-3)then
        d=1.d0/(x-1.d0)**3
        dg=dlog(x)/(x-1.d0)**4
        esm=x*(-18.d0+11.d0*x+x**2)*d/12.d0
     C        +x**2*(15.d0-16.d0*x+4.d0*x**2)*dg/6.d0
        esm=esm-2.d0/3.d0*dlog(x)
        else
        esm=43.d0/72.d0
        endif
        return
        end

C********************************************************************

        double precision function eh(x)
	
	implicit none
        double precision x,d,dg
        if(dabs(x-1.d0).gt.1.d-3)then
        d=1.d0/(x-1.d0)**3
        dg=dlog(x)/(x-1.d0)**4
        eh=x*(16.d0-29.d0*x+7.d0*x**2)*d/36.d0+x*(3.d0*x-2.d0)
     C       *dg/6.d0
        else
        eh=1.d0/8.d0
        endif
        return
        end
	
C********************************************************************

        double precision function af(x)
	
	implicit none
	double precision x,x2,x3,lnx,pi
        pi=4.d0*datan(1.d0)
	x2=x**2
	x3=x**3
	lnx=dlog(x)
        af=16.d0/9.d0*((5.d0/2.d0-1.d0/3.d0*pi**2-3.d0*1.2021d0
     C     +(5.d0/2.d0-3.d0/4.d0*pi**2)*lnx+lnx**2/4.d0
     C                 +lnx**3/12.d0)*x
     C    +(7.d0/4.d0+2.d0/3.d0*pi**2-1.d0/2.d0*pi**2*lnx-lnx**2/4.d0
     C     +lnx**3/12.d0)*x2
     C    +(-7.d0/6.d0-pi**2/4.d0+2.d0*lnx-3.d0/4.d0*lnx**2)*x3)    
        return
        end

C********************************************************************

        double precision function afim(x)
	
	implicit none
	double precision x,x2,x3,lnx,pi
        pi=4.d0*datan(1.d0)
	x2=x**2
	x3=x**3
	lnx=dlog(x)
        afim=16.d0/9.d0*pi*((2.d0-pi**2/6.d0+lnx/2.d0+lnx**2/2.d0)*x
     C          +(1.d0/2.d0-pi**2/6.d0-lnx+lnx**2/2.d0)*x2+x3)
        return
        end

C********************************************************************

        double precision function bf(x)
	
	implicit none
	double precision x,x2,x3,lnx,pi
        pi=4.d0*datan(1.d0)
	x2=x**2
	x3=x**3
	lnx=dlog(x)
        bf=-8.d0/9.d0*((-3.d0+pi**2/6.d0-lnx)*x
     C    +(1.d0/2.d0+pi**2-2.d0*lnx-lnx**2/2.d0)*x2
     C    +(-25.d0/12.d0-pi**2/9.d0-19.d0/18.d0*lnx+2.d0*lnx**2)*x3
     C    -2.d0/3.d0*pi**2*x**(3.d0/2.d0))
        return
        end

C********************************************************************

        double precision function bfim(x)
	
	implicit none
	double precision x,x2,x3,lnx,pi
        pi=4.d0*datan(1.d0)
	x2=x**2
	x3=x**3
	lnx=dlog(x)
        bfim=-8.d0/9.d0*pi*(-x+(1.d0-2.d0*lnx)*x2
     C       +(-10.d0/9.d0+4.d0/3.d0*lnx)*x3)
        return
        end

C********************************************************************
C NLO t-quark contribution to K
        double precision function A1(x)
	
	implicit none
	double precision x,x2,x3,x4,x5,dd4,dd5,lnx,sp2,spx
	x2=x**2
	x3=x**3
	x4=x**4
	x5=x**5
	dd4=(x-1.d0)**4
	dd5=(x-1.d0)**5
	lnx=dlog(x)
	spx=sp2(1.d0-1.d0/x)
        A1=(32.d0*x4+244.d0*x3-160.d0*x2+16.d0*x)/(9.d0*dd4)*spx
     C   -(-774.d0*x4-2826.d0*x3+1994.d0*x2-130.d0*x+8.d0)
     C                              /(81.d0*dd5)*lnx
     C   +(-94.d0*x4-18665.d0*x3+20682.d0*x2-9113.d0*x+2006.d0)
     C                              /(243.d0*dd4)
        return
        end

C********************************************************************
C NLO t-quark contribution to K
        double precision function F1(x)
	
	implicit none
	double precision x,x2,x3,x4,x5,dd4,dd5,lnx,sp2,spx
	x2=x**2
	x3=x**3
	x4=x**4
	x5=x**5
	dd4=(x-1.d0)**4
	dd5=(x-1.d0)**5
	lnx=dlog(x)
	spx=sp2(1.d0-1.d0/x)
        F1=(4.d0*x4-40.d0*x3-41.d0*x2-x)/(3.d0*dd4)*spx
     C   -(-144.d0*x4+3177.d0*x3+3661.d0*x2+250.d0*x-32.d0)
     C                   /(108.d0*dd5)*lnx
     C   +(-247.d0*x4+11890.d0*x3+31779.d0*x2-2966.d0*x+1016.d0)
     C                   /(648.d0*dd4)
        return
        end
	
C********************************************************************
C NLO t-quark contribution to K
        double precision function dA0(x)
	
	implicit none
	double precision x,x2,x3,x4,x5,dd4,dd5,lnx,sp2,spx
	x2=x**2
	x3=x**3
	x4=x**4
	x5=x**5
	dd4=(x-1.d0)**4
	dd5=(x-1.d0)**5
	lnx=dlog(x)
	spx=sp2(1.d0-1.d0/x)
        dA0=x*(3.d0*x2+5.d0*x-4)/(2.d0*dd5)*lnx
     C   +(-141.d0*x2+48.d0*x+21.d0)/(36.d0*dd4)
        return
        end
	
C*********************************************************************
C NLO t-quark contribution to K
        double precision function dF0(x)
	
	implicit none
	double precision x,x2,x3,x4,x5,dd4,dd5,lnx,sp2,spx
	x2=x**2
	x3=x**3
	x4=x**4
	x5=x**5
	dd4=(x-1.d0)**4
	dd5=(x-1.d0)**5
	lnx=dlog(x)
	spx=sp2(1.d0-1.d0/x)
        dF0=-3.d0*x*(x+1.d0)/dd5*lnx
     C   +(x2+10.d0*x+1.d0)/(2.d0*dd4)
        return
        end
		
C********************************************************************

	double precision function ffh(x,tanb)
	
	implicit none
	double precision x,x2,x3,x4,x5,dd3,dd4,dd5,lnx,ln2x,spx
	double precision sp2,sum1,sum2,tanb
	x2=x**2
	x3=x**3
	x4=x**4
	x5=x**5
	dd3=(x-1.d0)**3
	dd4=(x-1.d0)**4
	dd5=(x-1.d0)**5
	lnx=dlog(x)
	ln2x=lnx**2
	spx=sp2(1.d0-1.d0/x)
	sum1=4.d0*(-3.d0+7.d0*x-2.d0*x2)/(3.d0*dd3)*spx+
     C       (8.d0-14.d0*x-3.d0*x2)/(3.d0*dd4)*ln2x+
     C       2.d0*(-3.d0-x+12.d0*x2-2.d0*x3)/(3.d0*dd4)*lnx
     C       +(7.d0-13.d0*x+2.d0*x2)/dd3
	sum2=x*(18.d0-37.d0*x+8.d0*x2)/dd4*spx
     C        +x*(-14.d0+23.d0*x+3.d0*x2)/dd5*ln2x
     C        +(-50.d0+251.d0*x-174.d0*x2-192.d0*x3+21.d0*x4)
     C                        /(9.d0*dd5)*lnx+
     C        (797.d0-5436.d0*x+7569.d0*x2-1202.d0*x3)/(108.d0*dd4)
	ffh=-4.d0/3.d0*x*sum1+2.d0/9.d0*x/(tanb**2)*sum2
	return
	end

C**********************************************************************

	double precision function fgh(x,tanb)
	
	implicit none
	double precision x,x2,x3,x4,dd2,dd3,dd4,dd5,lnx,ln2x,sp2
	double precision spx,sum1,sum2,tanb
	x2=x**2
	x3=x**3
	x4=x**4	
	dd2=(x-1.d0)**2
	dd3=(x-1.d0)**3
	dd4=(x-1.d0)**4
	dd5=(x-1.d0)**5
	lnx=dlog(x)
	ln2x=lnx**2
	spx=sp2(1.d0-1.d0/x)
	sum1=(-36.d0+25.d0*x-17.d0*x**2)/(2.d0*dd3)*spx+(19.d0+17.d0*x)
     C     /(dd4)*ln2x+(-3.d0-187.d0*x+12.d0*x2-14.d0*x3)/(4.d0*dd4)*lnx
     C        +3.d0*(143.d0-44.d0*x+29.d0*x2)/(8.d0*dd3)
	sum2=x*(30.d0-17.d0*x+13.d0*x2)/dd4*spx
     C      -x*(31.d0+17.d0*x)/dd5*ln2x
     C      +(-226.d0+817.d0*x+1353.d0*x2+318.d0*x3+42.d0*x4)
     C       /(36.d0*dd5)*lnx+
     C      +(1130.d0-18153.d0*x+7650.d0*x2-4451.d0*x3)/(216.d0*dd4)
	fgh=-x/3.d0*sum1+x/(6.d0*tanb**2)*sum2
	return
	end

C*********************************************************************

C For effective Yukawa Couplings
	double precision function h2(x,y)
	
	implicit none
	double precision x,y
	if(dabs(x-y).gt.1.d-2)then
	  if(dabs(x-1.d0).lt.1.d-2)then
	         h2=1.d0/(-1.d0+y)-y*dlog(y)/(-1.d0+y)**2
	   else 
	      if(dabs(y-1.d0).lt.1.d-2)then
	         h2=1.d0/(-1.d0+x) - x*dlog(x)/(-1.d0+x)**2
	      else      
	         h2=x*dlog(x)/((1.d0-x)*(x-y))+ 
     C                y*dlog(y)/((1.d0-y)*(-x+y))
	      endif 
	   endif  
	else
	   if(dabs(x-1.d0).lt.1.d-2)then   
	      h2=-1.d0/2.d0                
	   else
	      h2=1.d0/(1.d0 - x) + dlog(x)/(-1.d0 + x)**2
	   endif
	endif      
	return
	end

C******************************************************************

	double precision function BB0(x,y,z,t)

	implicit none
	double precision x,y,z,t,B0
	BB0=1.d0-B0(x,y,z,t)
	return
	end

C******************************************************************

	double precision function BB1(x,y,z,t)

	implicit none
	double precision x,y,z,t,BB0
	if(dabs(y-z).lt.1.d-5)then
	BB1=1.d0/2.d0*BB0(x,y,z,t)
	else 
	BB1=1.d0/2.d0*BB0(x,y,z,t)-(y+z)/(4.d0*(y-z))
     C                      +y*z/(2.d0*(y-z)**2)*dlog(y/z)
	endif
	return
	end

C******************************************************************

     	double precision function S0(x)

	implicit none
	double precision x
	S0=x*(1.d0/4.d0+9.d0/(4.d0*(1.d0-x))-3.d0/(2.d0*(1.d0-x)**2)
     C       -3.d0*x**2*dlog(x)/(2.d0*(1.d0-x)**3))
	return
	end
	
C******************************************************************

	double precision function gg0(x)

	implicit none
	double precision x
	if(dabs(x-1.d0).lt.1.d-5)then
	gg0=1.d0/12.d0
	else
	gg0=x*(x**2-1.d0-2.d0*x*dlog(x))/(4.d0*(x-1.d0)**3)
	endif
	return
	end
	
C******************************************************************

	double precision function ffp(x,y)

	implicit none
	double precision x,y
	if(dabs(y-x).lt.1.d-5)then
	 if(dabs(y-1.d0).lt.1.d-5)then
	 ffp=1.d0/6.d0
	 else
	 ffp=(-1.d0+x**2-2.d0*x*dlog(x))
     C        /(-2.d0*x+6.d0*x**2+2.d0*x**4-6.d0*x**3)
	 endif
	else
	 if(dabs(y-1.d0).lt.1.d-5)then
	 ffp=(-2.d0*x+x*dlog(x)+2.d0+dlog(x))/(x-1.d0)**2
	 else
	  if(dabs(x-1.d0).lt.1.d-5)then
	  ffp=(-1.d0*y**2-2.d0*y*dlog(y))/(2.d0*(y-1.d0)**3)
	  else
	  ffp=(x**2-y)*dlog(x)/((x-y)**2*(x-1.d0)**2)
     C      -y*dlog(y)/((x-y)**2*(y-1.d0))-1.d0/((x-y)*(x-1.d0))
	  endif
	 endif
	endif
	return
	end
	
C******************************************************************

	double precision function ggp(x,y,z)

	implicit none
	double precision x,y,z
	if(dabs(z-1.d0).le.1.d-5)then
	  if(dabs(x-y).le.1.d-5)then
	   if(dabs(1.d0-x).le.1.d-5)then
	   ggp=1.d0/3.d0
	   else
	   ggp=(x**2-2.d0*dlog(x)*x-1.d0)/(3.d0*x-1.d0+x**3-3.d0*x**2)
	   endif
	  else
	   if(dabs(x-1.d0).le.1.d-5)then
	   ggp=(2.d0*y**2*dlog(x)+4.d0*y-1.d0-3.d0*y**2)
     C           /(-2.d0+6.d0*y+2.d0*y**3-6.d0*y**2)
	   else if(dabs(1.d0-y).le.1.d-5)then
	   ggp=(2.d0*x**2*dlog(x)+4.d0*x-1.d0-3.d0*x**2)
     C           /(-2.d0+6.d0*x+2.d0*x**3-6.d0*x**2)
	   else
	   ggp=(x**2*dlog(x)-y**2*dlog(y)-dlog(y)*y**2*x**2-x**2
     C         -2.d0*dlog(x)*y*x**2-x*y**2+2.d0*dlog(y)*y**2*x+x**2*y
     C         +x**2*dlog(x)*y**2+x+y**2-y)
	   endif
	  endif
	else
	  if(dabs(x-y).le.1.d-5)then
	   if(dabs(x-1).le.1.d-5)then
	   ggp=(-1.d0-3.d0*z**2+4.d0*z+2.d0*z**2*dlog(z))
     C          /(-2.d0-6.d0*z**2+2.d0*z**3+6.d0*z)
	   else
	   ggp=x/((x-z)*(x-1.d0))*(1.d0-(1.d0/(x-1.d0)
     C           +z/(x-z))*dlog(x))+z**2*dlog(z)/((x-z)**2*(z-1.d0))
	   endif
	  else if(abs(x-z).le.1.d-5)then
	   if(abs(y-1.d0).le.1.d-5)then
	   ggp=(x**2-2.d0*dlog(x)*x-1.d0)/(x**3-3.d0*x**2+3.d0*x-1.d0)
	   else
	   ggp=x/((x-y)*(y-1.d0))*(1.d0-(1.d0/(y-1.d0)+y/(x-y))*dlog(x))
     C           +y**2*dlog(y)/((x-y)**2*(y-1.d0))	
	   endif
	  else if(abs(y-z).le.1.d-5)then
	   if(abs(x-1.d0).le.1.d-5)then
	   ggp=(y**2-2.d0*dlog(y)*y-1.d0)/(y**3-3.d0*y**2+3.d0*y-1.d0)
	   else
	   ggp=y/((y-x)*(x-1.d0))*(1.d0-(1.d0/(x-1.d0)+x/(y-x))*dlog(y))
     C           +x**2*dlog(x)/((y-x)**2*(x-1.d0))
           endif
	  else
	   if(dabs(y-1.d0).le.1.d-5)then
	   ggp=(-x*z**2+x**2*dlog(x)+x**2*dlog(x)*z**2
     C          +2.d0*z**2*dlog(z)*x-z-z**2*dlog(z)-z**2*dlog(z)*x**2
     C          -x**2+x**2*z-2.d0*x**2*dlog(x)*z+z**2+x)
     C       /(2.d0*x*z**3+3.d0*x**2*z-3.d0*x*z**2-2.d0*x**2-z**3
     C         +2.d0*z**2+x-z+x**3*z**2-x**2*z**3-2.d0*x**3*z+x**3)
	   else if(dabs(x-1.d0).le.1.d-5)then
	   ggp=(-y*z**2+y**2*dlog(y)+y**2*dlog(y)*z**2
     C          +2.d0*z**2*dlog(z)*y-z-z**2*dlog(z)-z**2*dlog(z)*y**2
     C          -y**2+y**2*z-2.d0*y**2*dlog(y)*z+z**2+y)
     C       /(2.d0*y*z**3+3.d0*y**2*z-3.d0*y*z**2-2.d0*y**2-z**3
     C         +2.d0*z**2+y-z+y**3*z**2-y**2*z**3-2.d0*y**3*z+y**3)
	   else
	   ggp=1.d0/(x-y)*(1.d0/(x-z)*(x**2/(x-1.d0)*dlog(x)-3.d0/2.d0*x
     C                             -z**2/(z-1.d0)*dlog(z)+3.d0/2.d0*z)
     C                 -1.d0/(y-z)*(y**2/(y-1.d0)*dlog(y)-3.d0/2.d0*y
     C                             -z**2/(z-1.d0)*dlog(z)+3.d0/2.d0*z))
	   endif
	  endif
	endif
	return
	end
	
C******************************************************************

	double precision function RUNMass(y,x)

	implicit none
	double precision x,y
	DOUBLE PRECISION PI,nf,ge0,ge1,b0,b1,aux,aux1,asf
	DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
	
	COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW

	PI=4.d0*DATAN(1.d0)

C	 * x: renormalization scale
C	 * y: quark running mass at 2 GeV;
C	Running Masses at 2 GeV (PDG 2006):
C	- Mu= 1.5 to 3.0 MeV
C	- Md= 3. to 7. MeV
C	- Ms= 95 +/- 25 MeV
C	- Mc= 1.25 +/- 0.09 GeV
C	- Mb= 4.20 +/- 0.07 GeV
	
	
	aux=y
     
	 nf=4.d0
	 ge0=8.d0
	 ge1=404.d0/3.d0-38.d0/3.d0*nf
	 b0=11.d0-2.d0*nf/3.d0
	 b1=102.d0-38.d0/3.d0*nf
	 aux1=2.d0
	 
	 if(x.gt.MBP)then
	  aux=aux*(asf(MBP)/asf(MC))**(ge0/(2.d0*b0))
     C      *(1.d0+ge0/(8.d0*pi*b0)*(ge1/ge0-b1/b0)
     C                             *(asf(MBP)-asf(aux1)))
     
	  nf=5.d0
	  ge0=8.d0
	  ge1=404.d0/3.d0-38.d0/3.d0*nf
	  b0=11.d0-2.d0*nf/3.d0
	  b1=102.d0-38.d0/3.d0*nf
	  aux1=MBP
	  
	  if(x.gt.MT)then
	   aux=aux*(asf(MT)/asf(MBP))**(ge0/(2.d0*b0))
     C      *(1.d0+ge0/(8.d0*pi*b0)*(ge1/ge0-b1/b0)
     C                             *(asf(MT)-asf(aux1)))
     
	   nf=6.d0
	   ge0=8.d0
	   ge1=404.d0/3.d0-38.d0/3.d0*nf
	   b0=11.d0-2.d0*nf/3.d0
	   b1=102.d0-38.d0/3.d0*nf
	   aux1=MT
	  endif
	  
	 endif
	
	
	RUNMass=aux*(asf(x)/asf(aux1))**(ge0/(2.d0*b0))
     C      *(1.d0+ge0/(8.d0*pi*b0)*(ge1/ge0-b1/b0)
     C                             *(asf(x)-asf(aux1)))
	  
	  
	return
	end
	
C******************************************************************

C Additional function Intpropa for BR(B-->Xs l+l-):

     	double precision function Intpropa(mH1,Gam1,mH2,Gam2,bd1,bd2)

	implicit none
	double precision mH1,Gam1,mH2,Gam2,bd1,bd2
	double precision mHmin,mHmax,Gam,mb1s,aux,fac,fca
	
	MB1S=4.68d0 
	if(mH1.gt.mH2)then
	mHmax=mH1
	mHmin=mH2
	Gam=Gam2
	else
	mHmax=mH2
	mHmin=mH1
	Gam=Gam1
	endif
	
	fac=dsqrt(mHmin*mHmax)
	IF(fac.gt.150.d0)then
	 aux=(bd2-bd1)/2*(bd1+bd2-4.d0/3.d0*(bd2**2+bd1*bd2+bd1**2)
     C    +(bd2**3+bd2**2*bd1+bd2*bd1**2+bd1**3)/2.d0)*(mb1s/fac)**4
	ELSE
	 fac=(mb1s/mHmin)**2
	 fca=(mb1s/mHmax)**2
	 if(abs(mHmin-mHmax).gt.1.d-4)then
	  if(mHmin.gt.5.d0)then
	  aux=(2.d0*(bd2-bd1)
     C          +mHmin**2/(mHmin**2-mHmax**2)*(1-fac)**2
     C  *dlog((1.d0-fac*bd2)**2/(1.d0-fac*bd1)**2)/fac)/(2.d0*fac)
     C      +(2.d0*(bd2-bd1)
     C          +mHmax**2/(mHmax**2-mHmin**2)*(1-fca)**2
     C  *dlog((1.d0-fca*bd2)**2/(1.d0-fca*bd1)**2)/fca)/(2.d0*fca)
     C        +(bd2-bd1)*(bd2+bd1-4.d0)/2.d0
	  else
	  aux=(2.d0*(bd2-bd1)+mHmin**2*(mHmin**2-mHmax**2)
     C      /((mHmin**2-mHmax**2)**2+(mHmin*Gam)**2)
     C      *dlog(((1-fac*bd2)**2+(Gam/mHmin)**2)
     C                /((1-fac*bd1)**2+(Gam/mHmin)**2))/fac
     C        *((1-fac)**2+mHmax**2*Gam**2
     C                      /(mHmin**2*(mHmin**2-mHmax**2))
     C        *(3.d0+(1.d0-fac*Gam**2/mb1s**2)*fac*fca-2.d0*fac-2*fca)
     C                                                ))/(2.d0*fac)
	  aux=aux+(2.d0*(bd2-bd1)+mHmax**2*(mHmax**2-mHmin**2)
     C      /((mHmin**2-mHmax**2)**2+(mHmin*Gam)**2)
     C      *dlog((1-fca*bd2)**2/(1-fca*bd1)**2)/fca
     C        *(1-fca)**2)/(2.d0*fca)
	  aux=aux+Gam*mHmin*mHmax**2/fac
     C            /((mHmax**2-mHmin**2)**2+(mHmin*Gam)**2)
     C *(datan((mHmin**2-mb1s**2*bd2)/(mHmin*Gam))
     C        -datan((mHmin**2-mb1s**2*bd1)/(mHmin*Gam)))
     C   *(3.d0-4.d0*fac+fac**2+2.d0*fca-2*(mHmin/mHmax)**2
     C     -(Gam/mb1s)**2*(2.d0*fca+fac-2.d0*fac*fca))
	  aux=aux+(bd2-bd1)*(bd1+bd2-4.d0)/2.d0
	  endif
	 else
	  if(mHmin.gt.5.d0)then
	  aux=(2.d0+(1.d0-fac)**2/((1.d0-fac*bd2)*(1.d0-fac*bd1)))
     C       *(bd2-bd1)+(3.d0-4.d0*fac+fac**2)/2.d0
     C      *dlog((1.d0-fac*bd2)**2/((1.d0-fac*bd1)**2))/fac
	  aux=aux/fac+(bd2-bd1)/2.d0*(bd2+bd1-4.d0)
	 else
	  aux=-(bd2-bd1)*(bd2+bd1-4.d0)/2.d0-(2.d0*(bd2-bd1)
     C +(3.d0-4.d0*fac**2+(1.d0+mHmin**2*Gam**2/mb1s**4)*fac**2)
     C          *dlog(((1-fac*bd2)**2+(Gam/mHmin)**2)
     C           /((1-fac*bd1)**2+(Gam/mHmin)**2))/(2.d0*fac))/fac
     C      +(datan((mb1s**2*bd2-mHmin**2)/(mHmin*Gam))
     C       -datan((mb1s**2*bd1-mHmin**2)/(mHmin*Gam)))
     C       /(fac**3*Gam*mHmin/mb1s**2)
     C  *((1.d0-fac)**2
     C     +(mHmin*Gam/mb1s**2)**2*fac**2*(-3.d0+2.d0*fac))
	 endif
	endif
	ENDIF
	aux=aux/mb1s**4
	intpropa=aux
	return
	end
	
