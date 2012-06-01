	SUBROUTINE RELDEN(PAR,PROB)

**********************************************************************
*   Subroutine for the computation of the dark matter relic density
*   PROB(30) =/= 0 excluded by WMAP,
*   PROB(30) = -1  LSP is not the lightest neutralino,
*   PROB(31) =/= 0 Higgs eff. self-couplings in Micromegas > 1.
*
**********************************************************************

	IMPLICIT NONE

	CHARACTER name*10,mess*20

	INTEGER NORD(5),HORD(3),NBIN,OMGFLAG,MAFLAG
	INTEGER sortOddParticles,err,i,j
        INTEGER nucleonAmplitudes

	DOUBLE PRECISION PAR(*),PROB(*)
	DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW,PI
	DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
	DOUBLE PRECISION SMASS(3),PMASS(2),CMASS,SCOMP(3,3),PCOMP(2,2)
	DOUBLE PRECISION MH(3),MA(2),MHC
	DOUBLE PRECISION MGL,MCHA(2),UU(2,2),VV(2,2),MNEU(5),NEU(5,5)  
	DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
	DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
	DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
	DOUBLE PRECISION SST,SSB,SSL,COSB,SINB,TANB
	DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
	DOUBLE PRECISION LQ,KQ,ALQ,AKQ,MUQ,NUQ,ALS
	DOUBLE PRECISION tab(250),OMG,OMGMIN,OMGMAX,Xf
	DOUBLE PRECISION sigmaV,x(100),dNdx(100),EMIN,LAM
	DOUBLE PRECISION sigmaPiN,sigma0,csPsi,csNsi,csPsd,csNsd
	DOUBLE PRECISION higgsPotent,darkOmega,calcSpectrum,zInterp
	DOUBLE PRECISION FeScLoop,LOPmass,NOFF,Nmass,SCcoeff
	DOUBLE PRECISION pA0(2),pA5(2),nA0(2),nA5(2),ffS0P(3),ffS0N(3)

	COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW 
	COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW 
	COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS 
	COMMON/HMO/MH,MA,MHC
	COMMON/SUSYSPEC/MGL,MCHA,UU,VV,MNEU,NEU 
	COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL, 
     C		MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT, 
     C		CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
	COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
	COMMON/QNMPAR/LQ,KQ,ALQ,AKQ,MUQ,NUQ
	COMMON/ALSHIFT/ALS
	COMMON/MICROMG/OMG,OMGMIN,OMGMAX,Xf,sigmaV,x,dNdx,EMIN,NBIN
	COMMON/MICROMG2/sigmaPiN,sigma0,csPsi,csNsi,csPsd,csNsd
	COMMON/FLAGS/OMGFLAG,MAFLAG
	COMMON/LAM/LAM

	EXTERNAL FeScLoop,LOPmass

	DATA NOFF/-12345.d0/
	DATA NORD/1,2,4,3,5/
	DATA HORD/2,1,3/

	IF (OMGFLAG.EQ.0) RETURN

*   Input parameters:

        PI=4.D0*DATAN(1.D0)
	TANB=PAR(3)
	COSB=1.D0/DSQRT(1.D0+TANB**2)
	SINB=TANB*COSB

	SST=DSQRT(1-CST**2) 
	SSB=DSQRT(1-CSB**2) 
	SSL=DSQRT(1-CSL**2) 

	CALL assignValW('alfEMZ',ALEMMZ)
	CALL assignValW('alfSMZ',ALSMZ)
	CALL assignValW('MbMb',MB)
	CALL assignValW('Mtp',MT)

	CALL assignValW('At',PAR(12))
	CALL assignValW('Ab',PAR(13))
	CALL assignValW('Al',PAR(14))

	CALL assignValW('Lambda',LQ/DSQRT(ZHU*ZHD*ZS))
	CALL assignValW('Kappa',KQ/DSQRT(ZS**3))
	CALL assignValW('tb',PAR(3))
	CALL assignValW('aLmbd0',PAR(5))

	CALL assignValW('Mha',MA(1))
	CALL assignValW('Mhb',MA(2))
	CALL assignValW('MHc',MHC)
	DO i=1,3
	 WRITE(name,fmt='(A2,I1)') 'Mh',i
	 CALL assignValW(name,MH(i))
	DO j=1,3
	 WRITE(name,fmt='(A2,I1,I1)') 'Zh',i,j
	 CALL assignValW(name,SCOMP(i,HORD(j)))
	ENDDO
	ENDDO
        CALL assignValW('Pa11',PCOMP(1,1))
        CALL assignValW('Pa12',PCOMP(1,2))
        CALL assignValW('Pa21',PCOMP(2,1))
        CALL assignValW('Pa22',PCOMP(2,2))

	DO i=1,5
	 WRITE(name,fmt='(A3,I1)') 'MNE',i
	 CALL assignValW(name,MNEU(i))
	 DO j=1,5
	   WRITE(name,fmt='(A2,I1,I1)') 'Zn',i,j
	   CALL assignValW(name,NEU(i,NORD(j)))
	  ENDDO
	ENDDO

	CALL assignValW('MSl1',MSL1)
	CALL assignValW('MSl2',MSL2)
	CALL assignValW('Zl11',CSL)
	CALL assignValW('Zl12',SSL)
	CALL assignValW('Zl21',-SSL)
	CALL assignValW('Zl22',CSL)

	CALL assignValW('MSb1',MSB1)
	CALL assignValW('MSb2',MSB2)
	CALL assignValW('Zb11',CSB)
	CALL assignValW('Zb12',SSB)
	CALL assignValW('Zb21',-SSB)
	CALL assignValW('Zb22',CSB)

	CALL assignValW('MSt1',MST1)
	CALL assignValW('MSt2',MST2)
	CALL assignValW('Zt11',CST)
	CALL assignValW('Zt12',SST)
	CALL assignValW('Zt21',-SST)
	CALL assignValW('Zt22',CST)

	CALL assignValW('MSeL',MLL)
	CALL assignValW('MSeR',MLR)
	CALL assignValW('MSmL',MLL)
	CALL assignValW('MSmR',MLR)
	CALL assignValW('MSne',MNL)
	CALL assignValW('MSnm',MNL)
	CALL assignValW('MSnl',MSNT)
	CALL assignValW('MSuL',MUL)
	CALL assignValW('MSuR',MUR)
	CALL assignValW('MSdL',MDL)
	CALL assignValW('MSdR',MDR)
	CALL assignValW('MScL',MUL)
	CALL assignValW('MScR',MUR)
	CALL assignValW('MSsL',MDL)
	CALL assignValW('MSsR',MDR)
	CALL assignValW('MSG',MGL)

*   Improved Higgs potential

	LAM=higgsPotent()
	IF(LAM.EQ.-1.D0 .OR. LAM.GE.1.D0) THEN 
	  PROB(31)=LAM
	ENDIF

*   Sorting sparticles

	err=sortOddParticles(mess)
	IF(err.ne.0 .OR. mess.ne.'~o1') THEN 
	  OMG=-1.D0
	  PROB(30)=-1.D0
	  RETURN
	ENDIF
	IF(DABS(MNEU(1)).LT.1.D0)THEN
	  OMG=0.D0
	  RETURN
	ENDIF

*   Computing relic density

	OMG=darkOmega(Xf,1,1.D-6) 
	PROB(30)=DDIM(OMG,OMGMAX)-DDIM(OMGMIN,OMG)

	IF (OMGFLAG.EQ.1) RETURN
	IF (OMGFLAG.EQ.3) GOTO 1

*  Computing WIMP-Nucleon cross sections 
*  Muq/Mdq=0.553D0, Msq/Mdq=18.9D0

	CALL getScalarFF(0.553D0,18.9D0,sigmaPiN,sigma0,ffS0P,ffS0N)
	CALL setProtonFF(ffS0P,NOFF,NOFF)
	CALL setNeutronFF(ffS0N,NOFF,NOFF)
	err=nucleonAmplitudes(FeScLoop,pA0,pA5,nA0,nA5)
	Nmass=0.939D0

	SCcoeff=4/PI*3.8937966D8*(Nmass*lopmass()/(Nmass+ lopmass()))**2
	csPsi=SCcoeff*pA0(1)**2
	csNsi=SCcoeff*nA0(1)**2
	csPsd=3*SCcoeff*pA5(1)**2
	csNsd=3*SCcoeff*nA5(1)**2
	IF( pA0(1)*nA0(1) .lt. 0) csNsi=-csNsi
	IF( pA5(1)*nA5(1) .lt. 0) csNsd=-csNsd
       
	IF (OMGFLAG.EQ.2) RETURN

*  Computing indirect detection rate 

 1	sigmaV=calcSpectrum(1.D-3,0,tab,err)
	IF (err.EQ.0) sigmaV=0.D0
	IF (sigmaV.NE.0.D0) THEN
	 DO I=1,NBIN
	  dNdx(I)=zInterp(DLOG(10.D0)*x(I),tab)*DLOG(10.D0)
	 ENDDO
	ENDIF

	END


