	SUBROUTINE RGESUNI(PAR,IFAIL,GUTEST)

*   Subroutine to integrate the RGEs for all 21 soft terms
*   from the SUSY scale Q2 up to MGUT, through a call of the
*   subroutine ODEINTS that is part of the file integ.f
*   Q2 is either computed in terms of the first generation squarks, 
*   or put in by the user.
*
*   Present accuracy: 
*   Gauginos: two loops
*   Trilinear couplings and scalar masses: two loops, 
*   except for terms ~lambda**2 and kappa**2
*
*   MGUT and the gauge/Yukawa couplings at Q2 are read in from
*   COMMON/SUSYCOUP, initialized in RGES.
*
*   The soft terms at Q2 are read in from PAR(*).
*   MH1, MH2, MSS are read in from COMMON/SUSYMH.
*
*   PURPOSE: Once the gauge/Yukawa couplings at Q2 have been computed
*   in RGES (including SUSY thresholds), the agreement of the soft
*   terms at MGUT with the inputs is tested via GUTEST.
*
*   Also: MS**2 at MGUT is computed here and written into
*   COMMON/MSGUT/MSGUT
*********************************************************************** 

	IMPLICIT NONE

	INTEGER IFAIL,NN,I,IM
	PARAMETER (NN=30)

	DOUBLE PRECISION PAR(*),EPS,X1,X2,Y(NN),PI,COEF
	DOUBLE PRECISION MGUT,g1s,g2s,g3s,HTOPS,HBOTS,HTAUS,Q2
	DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
	DOUBLE PRECISION GUTEST,AKGUTIN,MSGUT,MHDGUTIN,MHUGUTIN,ALGUTIN
	DOUBLE PRECISION G1GUT,G2GUT,G3GUT,LGUT,KGUT,HTOPGUT
	DOUBLE PRECISION HBOTGUT,HTAUGUT,M1GUT,M2GUT,M3GUT,ALGUT,AKGUT
	DOUBLE PRECISION ATGUT,ABGUT,ATAUGUT,AMUGUT
	DOUBLE PRECISION MHUGUT,MHDGUT,MSSGUT,MQ3GUT
	DOUBLE PRECISION MU3GUT,MD3GUT,MQGUT,MUGUT,MDGUT,ML3GUT,ME3GUT
	DOUBLE PRECISION MLGUT,MEGUT,MSTOT
        DOUBLE PRECISION M0,M12,A0,MH1,MH2,MSS
	DOUBLE PRECISION M1GUTIN,M2GUTIN,M3GUTIN

	COMMON/MGUT/MGUT
	COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
	COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
	COMMON/RENSCALE/Q2
	COMMON/GUTCOUP/G1GUT,G2GUT,G3GUT,LGUT,KGUT,HTOPGUT,
     C        HBOTGUT,HTAUGUT
	COMMON/GUTPAR/M1GUT,M2GUT,M3GUT,ALGUT,AKGUT,ATGUT,ABGUT,
     C        ATAUGUT,AMUGUT,MHUGUT,MHDGUT,MSSGUT,MQ3GUT,MU3GUT,MD3GUT,
     C        MQGUT,MUGUT,MDGUT,ML3GUT,ME3GUT,MLGUT,MEGUT
        COMMON/MSGUT/MSGUT
	COMMON/SOFTGUT/M0,M12,A0
	COMMON/AKGUT/AKGUTIN
	COMMON/SUSYMH/MH1,MH2,MSS
	COMMON/MHGUT/MHDGUTIN,MHUGUTIN,ALGUTIN
	COMMON/MSAVE/MSTOT,IM
	COMMON/MGGUT/M1GUTIN,M2GUTIN,M3GUTIN

	EXTERNAL DERIVSS,RKQSS

	EPS=1.D-8
	PI=4.D0*DATAN(1.D0)
	COEF=1.D0/(16.D0*PI**2)
	IM=IM+1
	IF(IM.EQ.11)THEN
	 MSTOT=0.D0
	 IM=1
	ENDIF
	
* Definition of the couplings squared Y(I) at Q2=M_SUSY

	Y(1)=g1s
	Y(2)=g2s
	Y(3)=g3s
	Y(4)=PAR(1)**2
	Y(5)=PAR(2)**2
	Y(6)=HTOPS**2
	Y(7)=HBOTS**2
	Y(8)=HTAUS**2

* Definition of the soft terms Y(I) at Q2=M_SUSY

	Y(9)=PAR(20)
	Y(10)=PAR(21)
	Y(11)=PAR(22)
	Y(12)=PAR(5)
	Y(13)=PAR(6)
	Y(14)=PAR(12)
	Y(15)=PAR(13)
	Y(16)=PAR(14)
	Y(17)=MH1
	Y(18)=MH2
	Y(19)=MSS
	Y(20)=PAR(7)
	Y(21)=PAR(8)
	Y(22)=PAR(9)
	Y(23)=PAR(15)
	Y(24)=PAR(16)
	Y(25)=PAR(17)
	Y(26)=PAR(10)
	Y(27)=PAR(11)
	Y(28)=PAR(18)
	Y(29)=PAR(19)
	Y(30)=PAR(24)

	X1=0.D0
	X2=COEF*DLOG(MGUT**2/Q2)

	!PRINT*,"CALL RGESUNI"
	!PRINT*,""
	!PRINT*,"MSUSY =",DSQRT(Q2)
	!PRINT*,"G1 =",5.D0/3.D0*Y(1)
	!PRINT*,"G2 =",Y(2)
	!PRINT*,"G3 =",Y(3)
	!PRINT*,"L =",Y(4)
	!PRINT*,"K =",Y(5)
	!PRINT*,"HT =",Y(6)
	!PRINT*,"HB =",Y(7)
	!PRINT*,"HL =",Y(8)
	!PRINT*,"M1 =",Y(9)
	!PRINT*,"M2 =",Y(10)
	!PRINT*,"M3 =",Y(11)
	!PRINT*,"AL =",Y(12)
	!PRINT*,"AK =",Y(13)
	!PRINT*,"ATOP =",Y(14)
	!PRINT*,"ABOT =",Y(15)
	!PRINT*,"ATAU =",Y(16)
	!PRINT*,"AMUON =",Y(30)
	!PRINT*,"MHU =",Y(17)
	!PRINT*,"MHD =",Y(18)
	!PRINT*,"MS =",Y(19)
	!PRINT*,"MQ3 =",Y(20)
	!PRINT*,"MU3 =",Y(21)
	!PRINT*,"MD3 =",Y(22)
	!PRINT*,"MQ =",Y(23)
	!PRINT*,"MU =",Y(24)
	!PRINT*,"MD =",Y(25)
	!PRINT*,"ML3 =",Y(26)
	!PRINT*,"ME3 =",Y(27)
	!PRINT*,"ML =",Y(28)
	!PRINT*,"ME =",Y(29)
	!PRINT*,""

	CALL ODEINTS(Y,NN,X1,X2,EPS,DERIVSS,RKQSS,IFAIL)

	!PRINT*,"MGUT =",MGUT
	!PRINT*,"G1GUT =",5.D0/3.D0*Y(1)
	!PRINT*,"G2GUT =",Y(2)
	!PRINT*,"G3GUT =",Y(3)
	!PRINT*,"LGUT =",Y(4)
	!PRINT*,"KGUT =",Y(5)
	!PRINT*,"HTGUT =",Y(6)
	!PRINT*,"HBGUT =",Y(7)
	!PRINT*,"HLGUT =",Y(8)
	!PRINT*,"M1GUT =",Y(9)," M12 =",M12," ERR =",
!     C		(Y(9)-M12)**2/(1.D2+M12**2)
	!PRINT*,"M2GUT =",Y(10)," M12 =",M12," ERR =",
!     C		(Y(10)-M12)**2/(1.D2+M12**2)
	!PRINT*,"M3GUT =",Y(11)," M12 =",M12," ERR =",
!     C		(Y(11)-M12)**2/(1.D2+M12**2)
	!PRINT*,"ALGUT =",Y(12)," ALGUTIN =",ALGUTIN," ERR =",
!     C		(Y(12)-ALGUTIN)**2/(1.D2+ALGUTIN**2)
	!PRINT*,"AKGUT =",Y(13)," AKGUTIN =",AKGUTIN," ERR =",
!     C		(Y(13)-AKGUTIN)**2/(1.D2+AKGUTIN**2)
	!PRINT*,"ATOPGUT =",Y(14)," A0 =",A0," ERR =",
!     C		(Y(14)-A0)**2/(1.D2+A0**2)
	!PRINT*,"ABOTGUT =",Y(15)," A0 =",A0," ERR =",
!     C		(Y(15)-A0)**2/(1.D2+A0**2)
	!PRINT*,"ATAUGUT =",Y(16)," A0 =",A0," ERR =",
!     C		(Y(16)-A0)**2/(1.D2+A0**2)
	!PRINT*,"AMUGUT =",Y(30)," A0 =",A0," ERR =",
!     C		(Y(30)-A0)**2/(1.D2+A0**2)
	!PRINT*,"MHUGUT =",Y(17)," MHUGUTIN =",MHUGUTIN," ERR =",
!     C		(Y(17)-MHUGUTIN)**2/(1.D4+MHUGUTIN**2)
	!PRINT*,"MHDGUT =",Y(18)," MHDGUTIN =",MHDGUTIN," ERR =",
!     C		(Y(18)-MHDGUTIN)**2/(1.D4+MHDGUTIN**2)
	!PRINT*,"MSGUT =",Y(19)
	!PRINT*,"MQ3GUT =",Y(20)," M0 =",M0**2," ERR =",
!     C		(Y(20)-M0**2)**2/(1.D4+M0**4)
	!PRINT*,"MU3GUT =",Y(21)," M0 =",M0**2," ERR =",
!     C		(Y(21)-M0**2)**2/(1.D4+M0**4)
	!PRINT*,"MD3GUT =",Y(22)," M0 =",M0**2," ERR =",
!     C		(Y(22)-M0**2)**2/(1.D4+M0**4)
	!PRINT*,"MQGUT =",Y(23)," M0 =",M0**2," ERR =",
!     C		(Y(23)-M0**2)**2/(1.D4+M0**4)
	!PRINT*,"MUGUT =",Y(24)," M0 =",M0**2," ERR =",
!     C		(Y(24)-M0**2)**2/(1.D4+M0**4)
	!PRINT*,"MDGUT =",Y(25)," M0 =",M0**2," ERR =",
!     C		(Y(25)-M0**2)**2/(1.D4+M0**4)
	!PRINT*,"ML3GUT =",Y(26)," M0 =",M0**2," ERR =",
!     C		(Y(26)-M0**2)**2/(1.D4+M0**4)
	!PRINT*,"ME3GUT =",Y(27)," M0 =",M0**2," ERR =",
!     C		(Y(27)-M0**2)**2/(1.D4+M0**4)
	!PRINT*,"MLGUT =",Y(28)," M0 =",M0**2," ERR =",
!     C		(Y(28)-M0**2)**2/(1.D4+M0**4)
	!PRINT*,"MEGUT =",Y(29)," M0 =",M0**2," ERR =",
!     C		(Y(29)-M0**2)**2/(1.D4+M0**4)
	!PRINT*,""

	IF(IFAIL.NE.0)THEN
	 !PRINT*,"IFAIL =",IFAIL
	 !PRINT*,""
	 !PRINT*,""
	 IFAIL=12
	 RETURN
	ENDIF

* GUTEST:

	GUTEST=0.D0
	GUTEST=GUTEST+(Y(9)-M1GUTIN)**2/(1.D2+M1GUTIN**2)
	GUTEST=GUTEST+(Y(10)-M2GUTIN)**2/(1.D2+M2GUTIN**2)
	GUTEST=GUTEST+(Y(11)-M3GUTIN)**2/(1.D2+M3GUTIN**2)
	GUTEST=GUTEST+(Y(12)-ALGUTIN)**2/(1.D2+ALGUTIN**2)
	GUTEST=GUTEST+(Y(13)-AKGUTIN)**2/(1.D2+AKGUTIN**2)
	DO I=14,16
	 GUTEST=GUTEST+(Y(I)-A0)**2/(1.D2+A0**2)
	ENDDO
	GUTEST=GUTEST+(Y(30)-A0)**2/(1.D2+A0**2)
	GUTEST=GUTEST+(Y(17)-MHUGUTIN)**2/(1.D4+MHUGUTIN**2)
	GUTEST=GUTEST+(Y(18)-MHDGUTIN)**2/(1.D4+MHDGUTIN**2)
	DO I=20,29
	 GUTEST=GUTEST+(Y(I)-M0**2)**2/(1.D4+M0**4)
	ENDDO
	!PRINT*,"GUTEST =",GUTEST
	!PRINT*,""
	!PRINT*,""

* Couplings at the GUT scale

	G1GUT=Y(1)
	G2GUT=Y(2)
	G3GUT=Y(3)
	LGUT=Y(4)
	KGUT=Y(5)
	HTOPGUT=Y(6)
	HBOTGUT=Y(7)
	HTAUGUT=Y(8)
	
* Soft terms at the GUT scale

	M1GUT=Y(9)
	M2GUT=Y(10)
	M3GUT=Y(11)
	ALGUT=Y(12)
	AKGUT=Y(13)
	ATGUT=Y(14)
	ABGUT=Y(15)
	ATAUGUT=Y(16)
	MHUGUT=Y(17)
	MHDGUT=Y(18)
	MSSGUT=Y(19)
c	MSGUT=Y(19)
	MSTOT=MSTOT+Y(19)
	MSGUT=MSTOT/IM
	MQ3GUT=Y(20)
	MU3GUT=Y(21)
	MD3GUT=Y(22)
	MQGUT=Y(23)
	MUGUT=Y(24)
	MDGUT=Y(25)
	ML3GUT=Y(26)
	ME3GUT=Y(27)
	MLGUT=Y(28)
	MEGUT=Y(29)
	AMUGUT=Y(30)

	END
