!*************************************************************************
! This file is part of
!
!       HiggsBounds 1.2.0
!
! by Philip Bechtle, Oliver Brein, Sven Heinemyer, Georg Weiglein
!    and Karina E. Williams.
!
!  Journal Reference: e-Print: arXiv:0811.4169 [hep-ph], submitted to CPC.
!  Web site: http://www.ippp.dur.ac.uk/HiggsBounds
!
!10/09/2009
!*************************************************************************


      PROGRAM HBwithCPsuperH
************************************************************************
* This is modified version of the cpsuperh2.f file which is supplied with 
* CPsuperH2.0 (downloaded from http://www.hep.man.ac.uk/u/jslee/CPsuperH.html)
* This file is part of the HiggsBounds distribution.
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
      REAL*8 SMPARA_H(19),SSPARA_H(26)
      DATA NSMIN/19/
      DATA NSSIN/26/
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

*-----------------------------------------------------------------------
*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *   *    
* used by initialize_HiggsBounds and run_HiggsBounds_part
* HB input:
        integer nHiggs        
        character*5 whichexperiment
        double precision Mh(3),                            
     &          CS_lep_hjZ_ratio(3),                              
     &          CS_lep_hjhi_ratio_nHbynH(3,3),                  
     &          CS_tev_gg_hj_ratio_dummy(3),
     &          CS_tev_bb_hj_ratio(3),   
     &          CS_tev_bg_hjb_ratio(3),                        
     &          CS_tev_ud_hjWp_ratio(3),CS_tev_cs_hjWp_ratio(3),
     &          CS_tev_ud_hjWm_ratio(3),CS_tev_cs_hjWm_ratio(3), 
     &          CS_tev_dd_hjZ_ratio(3),CS_tev_uu_hjZ_ratio(3),  
     &          CS_tev_ss_hjZ_ratio(3),CS_tev_cc_hjZ_ratio(3),  
     &          CS_tev_bb_hjZ_ratio(3),                       
     &          CS_tev_pp_vbf_ratio(3),                        
     &          BR_hjbb(3),BR_hjtautau(3),                      
     &          BR_hjWW(3),BR_hjgaga(3),                        
     &          BR_hjhihi_nHbynH(3,3)  


* HB output:
        integer HBresult,chan,ncombined  
        double precision obsratio
* misc:
        integer i,j
        double precision betasq
        double precision
     &          g2hjVV(3),g2hjbb(3),               
     &          g2hjhiZ_nHbynH(3,3)

c The number of neutral Higgs bosons in the MSSM is 3, therefore set
        nHiggs=3

c The string 'whichexperiment' determines which subset of experimental 
c results are used.
c In this example, we've used the option 'onlyL',
c which instructs HiggsBounds to use tables of results
c from LEP only (i.e. no Tevatron results).
        whichexperiment='onlyL'

c The subroutine initialize_HiggsBounds reads in all necessary
c tables etc.
c It must be called before calling the run_HiggsBounds_part subroutine.

        call initialize_HiggsBounds(nHiggs,whichexperiment) 

c If you would like to perform scans over variables, the subroutine
c initialize_HiggsBounds (and finish_HiggsBounds) should be called
c outside the do-loops in order to save time.
*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
*-----------------------------------------------------------------------
*read in input data from the file 'HBwithCPsuperH.input'
      write(*,*) 'call: HBwithCPsuperH < HBwithCPsuperH-before10June2009.input'
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
      READ(*,*) IFLAG_H(15) ! 1 to print HiggsEDM results
      READ(*,*) IFLAG_H(16) ! 1 to print FILLBOBS results
      READ(*,*) IFLAG_H(17) ! 1 to print BTSGAM   results
      READ(*,*) IFLAG_H(18) ! 1 or 2 or 3 to print FILLEDMS results
*
      IFLAG_H(20) = (NMNH-1)/2      ! ISMN
      IFLAG_H(21) = (NMNH-1)/2      ! ISUSYN = ISMN
      IFLAG_H(22) = (NMCH-1)/2      ! ISMC
      IFLAG_H(23) = (NMCH-1)/2      ! ISUSYC = ISMC
*=======================================================================
      CALL FILLPARA2(NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H)
      MCH=SSPARA_H(2)
      IF(IFLAG_H(57).GT.0) THEN
       write(*,*)'ERROR! IFLAG_H(57) = ',IFLAG_H(57)
       IFLAG_H(57)=0
       STOP
      ENDIF
*=======================================================================
      CALL FILLHIGGS2(NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H,MCH
     .               ,HMASS_H,OMIX_H)
*
      MCH=RAUX_H(10) ! Charged Higgs pole mass or effetive-pot. mass
      IERR1=IFLAG_H(50)+IFLAG_H(51)+IFLAG_H(52)+IFLAG_H(54)
     .      +IFLAG_H(55)+IFLAG_H(60)
      IF(IERR1.GT.0) THEN
       IF(IFLAG_H(50).EQ.1) write(*,*)'ERROR! IFLAG_H(50) = ',IFLAG_H(50)
       IF(IFLAG_H(51).EQ.1) write(*,*)'ERROR! IFLAG_H(51) = ',IFLAG_H(51)
       IF(IFLAG_H(52).EQ.1) write(*,*)'ERROR! IFLAG_H(52) = ',IFLAG_H(52)
       IF(IFLAG_H(54).EQ.1) write(*,*)'ERROR! IFLAG_H(54) = ',IFLAG_H(54)
       IF(IFLAG_H(55).EQ.1) write(*,*)'ERROR! IFLAG_H(55) = ',IFLAG_H(55)
       IF(IFLAG_H(60).EQ.1) write(*,*)'ERROR! IFLAG_H(60) = ',IFLAG_H(60)
       IFLAG_H(50)=0
       IFLAG_H(51)=0
       IFLAG_H(52)=0
       IFLAG_H(54)=0
       IFLAG_H(55)=0
       IFLAG_H(60)=0
       STOP
      ENDIF
      IF(IFLAG_H(53).EQ.1) THEN
       write(*,*)'WARNING! IFLAG_H(53) = ',IFLAG_H(53)
       IFLAG_H(53)=0
      ENDIF
*=======================================================================
      CALL FILLCOUPL2(NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H,MCH
     . ,HMASS_H,OMIX_H,STMASS_H,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H
     . ,STAUMIX_H,SNU3MASS_H,MC_H,UL_H,UR_H,MN_H,N_H,NCMAX,NHC_H,SHC_H
     . ,CHC_H)
      IERR2=IFLAG_H(54)+IFLAG_H(56)
      IF(IERR2.GT.0) THEN
       IF(IFLAG_H(54).EQ.1) write(*,*)'ERROR! IFLAG_H(54) = ',IFLAG_H(54)
       IF(IFLAG_H(56).EQ.1) write(*,*)'ERROR! IFLAG_H(56) = ',IFLAG_H(56)
       IFLAG_H(54)=0
       IFLAG_H(56)=0
       STOP
      ENDIF
*=======================================================================
      CALL FILLGAMBR2(NFLAG,IFLAG_H,MCH,HMASS_H,NCMAX,NHC_H,SHC_H
     . ,CHC_H,STMASS_H,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H,STAUMIX_H
     . ,SNU3MASS_H,MC_H,MN_H,NMNH,GAMBRN,NMCH,GAMBRC)
*=======================================================================
*For the subroutine FILLDHPG, we need \sqrt{s} value as an input. 
      SQRTS=HMASS_H(2)
*RAUX_H(101) is reserved for sqrt{s}
      RAUX_H(101)=SQRTS
*
      CALL FILLDHPG(SQRTS,NFLAG,IFLAG_H
     . ,MCH,HMASS_H,OMIX_H,NCMAX,NHC_H,SHC_H,CHC_H
     . ,STMASS_H,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H,STAUMIX_H
     . ,SNU3MASS_H,MC_H,UL_H,UR_H,MN_H,N_H,NMNH,GAMBRN,NMCH,GAMBRC)
*=======================================================================
      CALL HiggsEDM(NFLAG,IFLAG_H
     .             ,HMASS_H,STMASS_H,SBMASS_H,MC_H,OMIX_H,NCMAX,NHC_H)
*=======================================================================
      CALL FILLBOBS(NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H
     .             ,HMASS_H,OMIX_H,MCH,MC_H,UL_H,UR_H,STMASS_H,STMIX_H)
*=======================================================================
      CAUX_H(995)=(ATAU_H) ! A_E
      CAUX_H(996)=(AT_H)   ! A_U
      CAUX_H(997)=(AT_H)   ! A_C
      CAUX_H(998)=(AB_H)   ! A_D
      CAUX_H(999)=(AB_H)   ! A_S

      CALL FILLEDMS(NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H
     . ,MCH,HMASS_H,OMIX_H
     . ,STMASS_H,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H
     . ,STAUMIX_H,SNU3MASS_H
     . ,MC_H,UL_H,UR_H,MN_H,N_H,NCMAX,NHC_H,SHC_H,CHC_H
     . ,CAUX_H(995),CAUX_H(996),CAUX_H(997),CAUX_H(998),CAUX_H(999))
*=======================================================================
*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *      
c Set variables needed by HiggsBounds (using results from CPsuperH).
c See HiggsBounds documentation for definition of variables used
c as arguments to run_HiggsBounds_part and CPsuperH 
c documentation for all other variables.

c Note: It is slightly more accurate to use the subroutine run_HiggsBounds_part 
c rather than the subroutine run_HiggsBounds_effC because the SM branching ratios
c used internally in HiggsBounds (from HDecay) are not identical to the SM branching
c ratios used in CPsuperH

        do i=1,3

         Mh(i)=HMASS_H(i)         

         BR_hjbb(i)     = GAMBRN(6,2,i)
         BR_hjtautau(i) = GAMBRN(3,2,i) 
         BR_hjWW(i)     = GAMBRN(10,2,i) 
         BR_hjgaga(i)   = GAMBRN(17,2,i) 
         
         betasq=1.0D0-4.0D0*(MBMT_H/HMASS_H(i))**2.0D0
         g2hjbb(i)=     
     &      abs(NHC_H(17,i))**2.0D0   
     &    + abs(NHC_H(18,i))**2.0D0/betasq

         CS_tev_bg_hjb_ratio(i) = g2hjbb(i)
         CS_tev_bb_hj_ratio(i)  = g2hjbb(i)

         g2hjVV(i)= abs(NHC_H(70,i))**2.0D0

         CS_lep_hjZ_ratio(i)        = g2hjVV(i)
         CS_tev_dd_hjZ_ratio(i)     = g2hjVV(i)
         CS_tev_uu_hjZ_ratio(i)     = g2hjVV(i)
         CS_tev_ss_hjZ_ratio(i)     = g2hjVV(i)
         CS_tev_cc_hjZ_ratio(i)     = g2hjVV(i)
         CS_tev_bb_hjZ_ratio(i)     = g2hjVV(i)
         CS_tev_ud_hjWp_ratio(i)    = g2hjVV(i)
         CS_tev_cs_hjWp_ratio(i)    = g2hjVV(i)
         CS_tev_ud_hjWm_ratio(i)    = g2hjVV(i)
         CS_tev_cs_hjWm_ratio(i)    = g2hjVV(i)
         CS_tev_pp_vbf_ratio(i)     = g2hjVV(i)
c ------------------------------------------------------------------
c We set CS_tev_gg_hj_ratio to zero because it will not be required 
c by the option 'onlyL'
         CS_tev_gg_hj_ratio_dummy(i) = 0.0D0
c ------------------------------------------------------------------

	 BR_hjhihi_nHbynH(i,1)=GAMBRN(14,3,i)
	 BR_hjhihi_nHbynH(i,2)=GAMBRN(16,3,i)	 
	 BR_hjhihi_nHbynH(i,3)=0.0D0

        enddo
	
	do j=1,3    
	 do i=1,3 
	  if(i.lt.j)then	  
	   g2hjhiZ_nHbynH(j,i)=g2hjVV(6-j-i)
	   g2hjhiZ_nHbynH(i,j)=g2hjhiZ_nHbynH(j,i)
	  else
	   g2hjhiZ_nHbynH(j,i)=0.0D0	   
	  endif    
	 enddo
	enddo 

	do j=1,3    
	 do i=1,3        
          CS_lep_hjhi_ratio_nHbynH(j,i) = g2hjhiZ_nHbynH(j,i)
	 enddo
	enddo 

*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
c call to run_HiggsBounds_part subroutine:

        call run_HiggsBounds_part(nHiggs,Mh,                  
     &          CS_lep_hjZ_ratio,                             
     &          CS_lep_hjhi_ratio_nHbynH,                     
     &          CS_tev_gg_hj_ratio,CS_tev_bb_hj_ratio,        
     &          CS_tev_bg_hjb_ratio,                          
     &          CS_tev_ud_hjWp_ratio,CS_tev_cs_hjWp_ratio,     
     &          CS_tev_ud_hjWm_ratio,CS_tev_cs_hjWm_ratio,     
     &          CS_tev_dd_hjZ_ratio,CS_tev_uu_hjZ_ratio,      
     &          CS_tev_ss_hjZ_ratio,CS_tev_cc_hjZ_ratio,      
     &          CS_tev_bb_hjZ_ratio,                          
     &          CS_tev_pp_vbf_ratio,                          
     &          BR_hjbb,BR_hjtautau,                          
     &          BR_hjWW,BR_hjgaga,                             
     &          BR_hjhihi_nHbynH,                             
     &          HBresult,chan,                                
     &          obsratio, ncombined                          )

        write(*,*)        
        write(*,*)'*************    HiggsBounds Results  **************'
        write(*,*) 
        write(*,*)'Is this parameter point excluded by either LEP'
        write(*,*)'or Tevatron data?'         
        write(*,*) HBresult, ',  where'
        write(*,*)'               0 = yes, it is excluded'
        write(*,*)'               1 = no, it has not been excluded'
        write(*,*)'              -1 = invalid parameter set'    
        write(*,*)
        write(*,*)'The process with the highest statistical sensitivity'
        write(*,*)'is'
        write(*,*) chan,'(see Key.dat)'
        write(*,*)'This process has a theoretical rate vs. limit of'
        write(*,*) obsratio
        write(*,*)
        write(*,*)'The number of Higgs bosons which have contributed to'
        write(*,*)'the theoretical rate of this process was'
        write(*,*) ncombined
        write(*,*)
        write(*,*)'See HiggsBounds documentation for more information.'
        write(*,*)'****************************************************'
        write(*,*)
       
*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
c deallocates arrays used by HiggsBounds:

        call finish_HiggsBounds

      STOP
      END
