      PROGRAM CPsuperH2
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
*-----------------------------------------------------------------------
*+CDE HC_ COMMON BLOCKS:
      COMMON /HC_SMPARA/ AEM_H,ASMZ_H,MZ_H,SW_H,ME_H,MMU_H,MTAU_H,MDMT_H
     .                  ,MSMT_H,MBMT_H,MUMT_H,MCMT_H,MTPOLE_H,GAMW_H
     .                  ,GAMZ_H,EEM_H,ASMT_H,CW_H,TW_H,MW_H,GW_H,GP_H
     .                  ,V_H,GF_H,MTMT_H
*
      COMMON /HC_RSUSYPARA/ TB_H,CB_H,SB_H,MQ3_H,MU3_H,MD3_H,ML3_H,ME3_H
*
      COMPLEX*16 MU_H,M1_H,M2_H,M3_H,AT_H,AB_H,ATAU_H
      COMMON /HC_CSUSYPARA/ MU_H,M1_H,M2_H,M3_H,AT_H,AB_H,ATAU_H
*
*NEW COMMON BLOCKS for V2
*
      REAL*8     RAUX_H(999)
      COMPLEX*16 CAUX_H(999)
      COMMON /HC_RAUX/ RAUX_H
      COMMON /HC_CAUX/ CAUX_H
      DATA NAUX/999/
*-----------------------------------------------------------------------
*ARRAYS:
      REAL*8 SMPARA_H(19),SSPARA_H(38)
      DATA NSMIN/19/
      DATA NSSIN/38/
*
      INTEGER*8 IFLAG_H(100)
      DATA NFLAG/100/
*
      REAL*8     HMASS_H(3),OMIX_H(3,3)
      REAL*8     STMASS_H(2),SBMASS_H(2),STAUMASS_H(2),SNU3MASS_H
      REAL*8     MC_H(2),MN_H(4)
      COMPLEX*16 STMIX_H(2,2),SBMIX_H(2,2),STAUMIX_H(2,2)
      COMPLEX*16 UL_H(2,2),UR_H(2,2),N_H(4,4)
*
      COMPLEX*16 NHC_H(100,3)  ! 100 = NCMAX
      REAL*8     SHC_H(100)
      COMPLEX*16 CHC_H(100)
      DATA NCMAX/100/
*
      REAL*8 GAMBRN(101,3,3)   ! 101 = IFLAG_H(20)+IFLAG_H(21)+1 = NMNH
*                                      ISMN       =ISUSYN        = 50
      REAL*8 GAMBRC(51,3)      !  51 = IFLAG_H(22)+IFLAG_H(23)+1 = NMCH
*                                      ISMC       =ISUSYC        = 25
      DATA NMNH/101/
      DATA NMCH/51/
      integer HBresult, chan,ncombined
      real*8 obsratio
*-----------------------------------------------------------------------
C! Micrmegas Print 
      open(321, FILE='micromegas.in',STATUS='UNKNOWN')
C! ===============
*read in input data from the file 'run'
*
      DO IP=1,NSMIN
       READ(*,*) SMPARA_H(IP)
      ENDDO
      DO IP=1,NSSIN
       READ(*,*) SSPARA_H(IP)
      ENDDO
*-----------------------------------------------------------------------
*initialize 
*
      DO IFLAG=1,NFLAG
       IFLAG_H(IFLAG)=0
      ENDDO
*
      DO INC=1,NCMAX
        DO IH=1,3
         NHC_H(INC,IH)=DCMPLX(0.D0,0.D0)
        ENDDO
         CHC_H(INC)=DCMPLX(0.D0,0.D0)
         SHC_H(INC)=0.D0
      ENDDO
*
      DO IM=1,NMNH
        DO IWB=1,3
          DO IH=1,3
            GAMBRN(IM,IWB,IH)=0.D0
          ENDDO
        ENDDO
      ENDDO
      DO IM=1,NMCH
        DO IWB=1,3
          GAMBRC(IM,IWB)=0.D0
        ENDDO
      ENDDO
*
      DO IAUX=1,NAUX
        RAUX_H(IAUX)=0.D0
        CAUX_H(IAUX)=DCMPLX(0.D0,0.D0)
      ENDDO
*-----------------------------------------------------------------------
*Set some flag's
*
      READ(*,*) IFLAG_H( 1) ! '1' will print input parameters
      READ(*,*) IFLAG_H( 2) ! '1' will print Higgs sector
      READ(*,*) IFLAG_H( 3) ! '1' will print masses and mixings of 
*                             stop and sbottom sectors
      READ(*,*) IFLAG_H( 4) ! '1' will print masses and mixings of chargino and 
*                             neutralino sectors
      READ(*,*) IFLAG_H( 5) ! '1-6' will print Higgs boson couplings
      READ(*,*) IFLAG_H( 6) ! '1-5' will print Higgs boson decays 
      READ(*,*) IFLAG_H(10) ! if 0, include radiative corrections to 
*                             H-top-top and H-bot-bot Yukawa couplings
      READ(*,*) IFLAG_H(11) ! if 0, use pole Higgs masses
*                             if 1, use effective potential Higgs mass
      READ(*,*) IFLAG_H(12) ! 5 For full improvemnt
       IF(IFLAG_H(12).EQ.0) IFLAG_H(12)=5
      READ(*,*) IFLAG_H(13) ! 1 Not to include the off-diagonal absorptive parts
      READ(*,*) IFLAG_H(14) ! 1 to print FILLDHPG results
      READ(*,*) IFLAG_H(16) ! 1 to print FILLBOBS results
      READ(*,*) IFLAG_H(17) ! 1 to print BTSGAM   results
      READ(*,*) IFLAG_H(18) ! 1 or 2 or 3 to print FILLEDMS results
      READ(*,*) IFLAG_H(19) ! 1 or 2 to print FILLMUON results
      
      write(*,*) 'after reading'
      write(*,*) 'IFLAG_H(13)=',IFLAG_H(13) 
      write(*,*) 'IFLAG_H(14)=',IFLAG_H(14) 
      write(*,*) 'IFLAG_H(15)=',IFLAG_H(15) 
      write(*,*) 'IFLAG_H(16)=',IFLAG_H(16) 
      write(*,*) 'IFLAG_H(18)=',IFLAG_H(18)
      write(*,*) 'IFLAG_H(19)=',IFLAG_H(19) 
*
      IFLAG_H(20) = (NMNH-1)/2      ! ISMN
      IFLAG_H(21) = (NMNH-1)/2      ! ISUSYN = ISMN
      IFLAG_H(22) = (NMCH-1)/2      ! ISMC
      IFLAG_H(23) = (NMCH-1)/2      ! ISUSYC = ISMC
*=======================================================================
      CALL FILLPARA2(NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H)
      MCH=SSPARA_H(2)
      IF(IFLAG_H(57).GT.0) THEN
       print*,'ERROR! IFLAG_H(57) = ',IFLAG_H(57)
       IFLAG_H(57)=0
       GOTO 99
      ENDIF
*=======================================================================
      CALL FILLHIGGS2(NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H,MCH
     .               ,HMASS_H,OMIX_H)
*
      MCH=RAUX_H(10) ! Charged Higgs pole mass or effetive-pot. mass
      IERR1=IFLAG_H(50)+IFLAG_H(51)+IFLAG_H(52)+IFLAG_H(54)
     .      +IFLAG_H(55)+IFLAG_H(60)
      IF(IERR1.GT.0) THEN
       IF(IFLAG_H(50).EQ.1) print*,'ERROR! IFLAG_H(50) = ',IFLAG_H(50)
       IF(IFLAG_H(51).EQ.1) print*,'ERROR! IFLAG_H(51) = ',IFLAG_H(51)
       IF(IFLAG_H(52).EQ.1) print*,'ERROR! IFLAG_H(52) = ',IFLAG_H(52)
       IF(IFLAG_H(54).EQ.1) print*,'ERROR! IFLAG_H(54) = ',IFLAG_H(54)
       IF(IFLAG_H(55).EQ.1) print*,'ERROR! IFLAG_H(55) = ',IFLAG_H(55)
       IF(IFLAG_H(60).EQ.1) print*,'ERROR! IFLAG_H(60) = ',IFLAG_H(60)
       IFLAG_H(50)=0
       IFLAG_H(51)=0
       IFLAG_H(52)=0
       IFLAG_H(54)=0
       IFLAG_H(55)=0
       IFLAG_H(60)=0
       GOTO 99
      ENDIF
      IF(IFLAG_H(53).EQ.1) THEN
       print*,'WARNING! IFLAG_H(53) = ',IFLAG_H(53)
       IFLAG_H(53)=0
      ENDIF
*=======================================================================
      CALL FILLCOUPL2(NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H,MCH
     . ,HMASS_H,OMIX_H,STMASS_H,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H
     . ,STAUMIX_H,SNU3MASS_H,MC_H,UL_H,UR_H,MN_H,N_H,NCMAX,NHC_H,SHC_H
     . ,CHC_H)
      IERR2=IFLAG_H(54)+IFLAG_H(56)
      IF(IERR2.GT.0) THEN
       IF(IFLAG_H(54).EQ.1) print*,'ERROR! IFLAG_H(54) = ',IFLAG_H(54)
       IF(IFLAG_H(56).EQ.1) print*,'ERROR! IFLAG_H(56) = ',IFLAG_H(56)
       IF(IFLAG_H(56).EQ.2) print*,'ERROR! IFLAG_H(56) = ',IFLAG_H(56)
       IF(IFLAG_H(56).EQ.3) print*,'ERROR! IFLAG_H(56) = ',IFLAG_H(56)
       IFLAG_H(54)=0
       IFLAG_H(56)=0
       GOTO 99
      ENDIF
*=======================================================================
      CALL FILLGAMBR2(NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H
     . ,MCH,HMASS_H,NCMAX,NHC_H,SHC_H,CHC_H,STMASS_H
     . ,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H,STAUMIX_H,SNU3MASS_H
     . ,MC_H,MN_H,NMNH,GAMBRN,NMCH,GAMBRC)
*=======================================================================
*For the subroutine FILLDHPG, we need \sqrt{s} value as an input. 
      SQRTS=HMASS_H(2)
*
      CALL FILLDHPG(SQRTS,NFLAG,IFLAG_H
     . ,MCH,HMASS_H,OMIX_H,NCMAX,NHC_H,SHC_H,CHC_H
     . ,STMASS_H,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H,STAUMIX_H
     . ,SNU3MASS_H,MC_H,UL_H,UR_H,MN_H,N_H,NMNH,GAMBRN,NMCH,GAMBRC)
*=======================================================================
      CALL FILLBOBS(NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H
     .             ,HMASS_H,OMIX_H,MCH,MC_H,UL_H,UR_H,STMASS_H,STMIX_H)
*=======================================================================
      CALL FILLEDMS(NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H
     . ,MCH,HMASS_H,OMIX_H
     . ,STMASS_H,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H
     . ,STAUMIX_H,SNU3MASS_H
     . ,MC_H,UL_H,UR_H,MN_H,N_H,NCMAX,NHC_H,SHC_H,CHC_H)
*=======================================================================
      CALL FILLMUON(NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H
     . ,MCH,HMASS_H,OMIX_H
     . ,STMASS_H,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H
     . ,STAUMIX_H,SNU3MASS_H
     . ,MC_H,UL_H,UR_H,MN_H,N_H,NCMAX,NHC_H,SHC_H,CHC_H)
*=======================================================================
 99   CONTINUE
*
       call applyHiggBounds(HMASS_H, GAMBRN,MHC_H,NHC_H, MBMT_H,  
     . HBresult, obsratio,chan,ncombined)
       write(321,fmt='(A5,2x,I2)') 'LEPex', HBresult-1
       write(321,fmt='(A6,2x,1PE16.8)') 'LEPrat', obsratio

C! Micromegas Print
* Higgs
       write(321,fmt='(A3,2x,1PE16.8)') 'MHc', SSPARA_H(2)

      do  I=1,3
        write(321,fmt='(A2,I1,2x,1PE16.8)') 'Mh', I,HMASS_H(I)  
      enddo  

      do I=1,3
      do J=1,3 
        write(321,fmt='(A2,I1,I1,2x,1PE16.8)') 'Zh', I,J,OMIX_H(I,J)
      enddo
      enddo
* sTop
      do  I=1,2
        write(321,fmt='(A3,I1,2x,1PE16.8)') 'MSt', I,STMASS_H(I)
      enddo

      do I=1,2
      do J=1,2
         write(321,fmt='(A2,I1,I1,A1,2x,1PE16.8)')
     >               'Zt',I,J,'r',DREAL(STMIX_H(J,I))
         write(321,fmt='(A2,I1,I1,A1,2x,1PE16.8)')
     >               'Zt',I,J,'i',-DIMAG(STMIX_H(J,I))
      enddo
      enddo

* sBot      
      do  I=1,2
        write(321,fmt='(A3,I1,2x,1PE16.8)') 'MSb', I,SBMASS_H(I)
      enddo

      do I=1,2
      do J=1,2
         write(321,fmt='(A2,I1,I1,A1,2x,1PE16.8)')
     >               'Zb',I,J,'r',DREAL(SBMIX_H(J,I))
         write(321,fmt='(A2,I1,I1,A1,2x,1PE16.8)')
     >               'Zb',I,J,'i',-DIMAG(SBMIX_H(J,I))
      enddo
      enddo

* sTau
      write(321,fmt='(A4,2x,1PE16.8)') 'MSnl', SNU3MASS_H
      do  I=1,2
        write(321,fmt='(A3,I1,2x,1PE16.8)') 'MSl', I,STAUMASS_H(I)
      enddo

      do I=1,2
      do J=1,2
         write(321,fmt='(A2,I1,I1,A1,2x,1PE16.8)')
     >               'Zl',I,J,'r',DREAL(STAUMIX_H(J,I))
         write(321,fmt='(A2,I1,I1,A1,2x,1PE16.8)')
     >               'Zl',I,J,'i',-DIMAG(STAUMIX_H(J,I))
      enddo
      enddo

* Neutralino
      do  I=1,4
        write(321,fmt='(A3,I1,2x,1PE16.8)') 'MNE', I,MN_H(I)
      enddo

      do I=1,4
      do J=1,4
        write(321,fmt='(A2,I1,I1,A1,2x,1PE16.8)') 
     >                 'Zn', I,J,'r',DREAL(N_H(I,J))
                write(321,fmt='(A2,I1,I1,A1,2x,1PE16.8)') 
     >                 'Zn', I,J,'i',DIMAG(N_H(I,J))
      enddo
      enddo
* Chargino 
      do  I=1,2
        write(321,fmt='(A2,I1,2x,1PE16.8)') 'MC', I,MC_H(I)
      enddo

      do I=1,2
      do J=1,2
        write(321,fmt='(A2,I1,I1,A1,2x,1PE16.8)')
     >                 'Zu', I,J,'r',DREAL(UL_H(I,J))
                write(321,fmt='(A2,I1,I1,A1,2x,1PE16.8)')
     >                 'Zu', I,J,'i',DIMAG(UL_H(I,J))
      enddo
      enddo

      do I=1,2
      do J=1,2
        write(321,fmt='(A2,I1,I1,A1,2x,1PE16.8)')
     >                 'Zv', I,J,'r',DREAL(UR_H(I,J))
                write(321,fmt='(A2,I1,I1,A1,2x,1PE16.8)')
     >                 'Zv', I,J,'i',DIMAG(UR_H(I,J))
      enddo
      enddo


      close(321)
C! End of Micromegas Print
      STOP
      END
