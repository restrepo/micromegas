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

c Definitions of NHMAX, TEV_NSCHN, TEV_NDCHN, LEP_NSCHN, LEP_NDCHN, etc. :
#include "tables-nhmax.F"

#include "input.F"
#include "LEP-bound.F"
#include "TEV-bound.F"



*************************************************************************
	program HiggsBounds
*************************************************************************
	implicit none
	logical useTEV,useLEP,useSingH,prefix
	character*100 fileprefix
	character*5 whichexpt
	integer id(20),i,j,inputmode,NHIGGS

	call write_HiggsBounds_info()

	call initialize_SM_decay_data()

	call parse_commandline(NHIGGS,whichexpt,
     &		inputmode,fileprefix,prefix)

	call calc_HiggsBounds(NHIGGS,whichexpt,
     &          inputmode,fileprefix,prefix)



	end

*************************************************************************
        subroutine write_HiggsBounds_info()
*************************************************************************
        implicit none

	write(*,*)
	write(*,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	write(*,*)'~                                                        ~'
	write(*,*)'~          HiggsBounds 1.1.0 (fortran77 version)         ~'
	write(*,*)'~                                                        ~'
	write(*,*)'~     Philip Bechtle, Oliver Brein, Sven Heinemeyer,     ~'
	write(*,*)'~            Georg Weiglein, Karina E. Williams          ~'
	write(*,*)'~                                                        ~'
	write(*,*)'~ arXiv:0811.4169, http://www.ippp.dur.ac.uk/HiggsBounds ~'
	write(*,*)'~                                                        ~'
	write(*,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	write(*,*)
	write(*,*)'HiggsBounds collects together results from '
	write(*,*)
	write(*,*)'    * the LEP collaborations and LEP Higgs Working Group'
	write(*,*)'    * the CDF and D0 Collaborations'
	write(*,*)'    * the program HDECAY (arXiv:hep-ph/9704448)'
	write(*,*)'    * TeV4LHC Higgs Working Group report'
	write(*,*)'      (see arXiv:hep-ph/0612172 and ref. therein)'
	write(*,*)

	end




*************************************************************************
	subroutine calc_HiggsBounds(NHIGGS,whichexpt,
     & 		inputmode,fileprefix,prefix)
*************************************************************************
	implicit none
	logical eof,useTEV,useLEP,useSingH,prefix,CMDLINE
	logical used_channel(0:LEP_NSCHN*NHMAX+LEP_NDCHN*NHMAX**2+TEV_NSCHN*NHMAX)
	character*100 fileprefix
	character*5 whichexpt
	integer id(20),i,j,k,inputmode,NHIGGS,maxid,minid,sanity_flag
	integer log2int,output_chan

	double precision Mh(NHMAX),Gamma_h(NHMAX)

	double precision BRhjWW(NHMAX),BRhjgamgam(NHMAX)
        double precision BRhjbbbar(NHMAX),BRhjtautau(NHMAX)

        double precision BRhjhihi(NHMAX,NHMAX)
        double precision BRhjhihi_compact(NHMAX**2-NHMAX)

	double precision R_hjZ_ee(NHMAX)	
        double precision R_hjhi_ee(NHMAX,NHMAX)
        double precision R_hjhi_ee_compact((NHMAX*(NHMAX+1))/2)

        double precision R_hj_gg(NHMAX),R_hj_bb(NHMAX)
        double precision R_hjb_bg(NHMAX)
        double precision R_hjWp_udb(NHMAX),R_hjWp_csb(NHMAX),
     &                   R_hjWm_dub(NHMAX),R_hjWm_scb(NHMAX)
        double precision R_hjZ_ddb(NHMAX),R_hjZ_uub(NHMAX),
     &                   R_hjZ_ssb(NHMAX),R_hjZ_ccb(NHMAX),R_hjZ_bbb(NHMAX)

        double precision R_sigma_H(NHMAX),R_sigma_Hb(NHMAX)
        double precision R_sigma_HW(NHMAX),R_sigma_HZ(NHMAX)
        double precision R_sigma_H_VBF(NHMAX)

	double precision g_hjbb_sq(NHMAX), g_hjtautau_sq(NHMAX), 
     &   	g_hjWW_sq(NHMAX), g_hjZZ_sq(NHMAX),
     &          g_hjgamgam_sq(NHMAX), g_hjgg_sq(NHMAX), 
     &		g_hjhiZ_sq(NHMAX,NHMAX)
        double precision g_hjhiZ_sq_compact((NHMAX*(NHMAX+1))/2)

	double precision add_dat(MAX_ADDITIONAL_COLUMNS)
	integer columns


        logical lep_allowed,lep_final_h_comb(NHMAX)
        integer lep_mlratio
        integer lep_result_flag,lep_num_H_combined,lep_sim_mass_sum(NHMAX)
        double precision lep_obsratio,lep_predratio,lep_model_R_sig_BR

        logical tev_allowed,tev_final_h_comb(NHMAX)
        integer tev_mlratio
        integer tev_result_flag,tev_num_H_combined
        double precision tev_obsratio,tev_predratio,tev_model_R_sig_BR

	logical allowed,H_combined(NHMAX)
        integer mlratio,experiment_id,best_result,num_H_combined
c experiment_id : 1 : LEP, 2 : Tevatron
        integer sim_mass_sum(NHMAX)
        double precision obsratio,predratio,model_R_sig_BR
#include "tables-tev_cb.h"
#include "tables-lep_cb.h"

	do i=0,LEP_NSCHN*NHIGGS+LEP_NDCHN*NHIGGS**2+TEV_NSCHN*NHIGGS
	  used_channel(i) = .FALSE.
	enddo

	call initialize_filestatus()

	call determine_topology_flags(whichexpt,
     &          useTEV,useLEP,useSingH)

	if(useTEV) call initialize_TEV_tables(NHIGGS)
	if(useLEP) call initialize_LEP_tables(NHIGGS)

	call open_output_file(30,fileprefix,prefix,NHIGGS,whichexpt)

	eof = .FALSE.

	do while(eof .eqv. .FALSE.)
	
	  call read_MH_GammaTot(NHIGGS,fileprefix,prefix,eof,
     &          id(1),Mh,Gamma_h)

	  maxid=id(1)
	  minid=id(1)

	  if(inputmode .eq. 1) then

	    call read_BR_1H(NHIGGS,fileprefix,prefix,eof,
     &            id(2),BRhjbbbar,BRhjtautau,BRhjWW,BRhjgamgam)
		maxid=max(maxid,id(2))
                minid=min(minid,id(2))
	    if(.not. useSingH) then
		call read_BR_2H(NHIGGS,fileprefix,prefix,eof,
     &		    id(3),BRhjhihi)
!         if NHIGGS =1 no file is read in, so id(3) has to get a value
	        if(NHIGGS .eq. 1) id(3)=id(1)
                maxid=max(maxid,id(3))
                minid=min(minid,id(3))
	    endif
	    if(useLEP) then 
		call read_LEP_HZ_CS_ratios(NHIGGS,fileprefix,prefix,eof,
     &	            id(4),R_hjZ_ee)
                maxid=max(maxid,id(4))
                minid=min(minid,id(4))
	    endif
	    if(.not. useSingH) then
	      if(useLEP) then
		call read_LEP_2H_CS_ratios(NHIGGS,fileprefix,prefix,eof,
     &	              id(5),R_hjhi_ee)
                maxid=max(maxid,id(5))
                minid=min(minid,id(5))
	      endif
	    endif
	    if(useTEV) then
		call read_TEV_H_0jet_partCS_ratios(NHIGGS,fileprefix,prefix,eof,
     &               id(6),R_hj_gg,R_hj_bb)
                maxid=max(maxid,id(6))
                minid=min(minid,id(6))
	    endif
	    if(useTEV) then
		call read_TEV_H_1jet_partCS_ratios(NHIGGS,fileprefix,prefix,eof,
     &              id(7),R_hjb_bg)
                maxid=max(maxid,id(7))
                minid=min(minid,id(7))
	    endif
	    if(useTEV) then 
		call read_TEV_H_vbf_hadCS_ratios(NHIGGS,fileprefix,prefix,eof,
     &              id(8),R_sigma_H_VBF)
                maxid=max(maxid,id(8))
                minid=min(minid,id(8))
            endif
	    if(useTEV) then 
		call read_TEV_HW_partCS_ratios(NHIGGS,fileprefix,prefix,eof,
     &              id(9),R_hjWp_udb,R_hjWp_csb,R_hjWm_dub,R_hjWm_scb)
                maxid=max(maxid,id(9))
                minid=min(minid,id(9))
	    endif
	    if(useTEV) then
		call read_TEV_HZ_partCS_ratios(NHIGGS,fileprefix,prefix,eof,
     &              id(10),R_hjZ_ddb,R_hjZ_uub,R_hjZ_ssb,R_hjZ_ccb,R_hjZ_bbb)
                maxid=max(maxid,id(10))
                minid=min(minid,id(10))
	    endif

	  elseif(inputmode .eq. 2) then
	    call read_couplings_sq_ratios(NHIGGS,fileprefix,prefix,eof,
     &          id(2), g_hjbb_sq, g_hjtautau_sq, g_hjWW_sq, g_hjZZ_sq,
     &          g_hjgamgam_sq, g_hjgg_sq, g_hjhiZ_sq)
            maxid=max(maxid,id(2))
            minid=min(minid,id(2))
	    if(.not. useSingH) then 
		call read_BR_2H(NHIGGS,fileprefix,prefix,eof,
     &	            id(3),BRhjhihi)
!         if NHIGGS =1 no file is read in, so id(3) has to get a value
	        if(NHIGGS .eq. 1) id(3)=id(1)
                maxid=max(maxid,id(3))
                minid=min(minid,id(3))
	    endif

!          elseif(inputmode .eq. 3) then
!c SLHA


	  elseif(inputmode .eq. 4) then

	    call read_BR_1H(NHIGGS,fileprefix,prefix,eof,
     &              id(2),BRhjbbbar,BRhjtautau,BRhjWW,BRhjgamgam)
                maxid=max(maxid,id(2))
                minid=min(minid,id(2))

	    if(.not. useSingH) then
		call read_BR_2H(NHIGGS,fileprefix,prefix,eof,
     &	            id(3),BRhjhihi)
!         if NHIGGS =1 no file is read in, so id(3) has to get a value
	        if(NHIGGS .eq. 1) id(3)=id(1)
                maxid=max(maxid,id(3))
                minid=min(minid,id(3))
	    endif

	    if(useLEP) then 
		call read_LEP_HZ_CS_ratios(NHIGGS,fileprefix,prefix,eof,
     &	            id(4),R_hjZ_ee)
                maxid=max(maxid,id(4))
                minid=min(minid,id(4))
	    endif

	    if(.not. useSingH) then
	      if(useLEP) then
		call read_LEP_2H_CS_ratios(NHIGGS,fileprefix,prefix,eof,
     &	            id(5),R_hjhi_ee)
                maxid=max(maxid,id(5))
                minid=min(minid,id(5))
	      endif
	    endif

	    if(useTEV) then
		call read_TEV_1H_hadCS_ratios(NHIGGS,fileprefix,prefix,eof,
     &              id(11),R_sigma_H,R_sigma_Hb,R_sigma_HW,R_sigma_HZ,
     &	            R_sigma_H_VBF)
                maxid=max(maxid,id(11))
                minid=min(minid,id(11))
	    endif
	  endif

	  call read_additional_data(NHIGGS,fileprefix,prefix,eof,
     &          id(12),add_dat,columns)
! if no columns in file additional.dat id(12) needs to get a value
	  if(columns .eq. 0) id(12)=id(1)
	  maxid=max(maxid,id(12))
          minid=min(minid,id(12))

	  if(maxid .ne. minid) then 
c	    write(*,*)'maxid=',maxid
c	    write(*,*)'minid=',minid
	    write(*,*)'missmatching reference numbers'
	    stop
	  endif
	

	  if(.not. eof) then

c	  write(*,*)
c	  write(*,345) id(1)
 345	  FORMAT('###### input line:',I5)

#ifdef VERBOSE
  	  write(*,99) id
 99	  FORMAT ('id=',20I2)

	  write(*,*)'Mh:'
	  call output_vector(NHIGGS,Mh)
	  write(*,*)'Gamma_h:'
	  call output_vector(NHIGGS,Gamma_h)

	  if((inputmode .eq. 1) .or. (inputmode .eq. 4)) then
	    write(*,*)'BRhjbbbar:'
	    call output_vector(NHIGGS,BRhjbbbar)
	    write(*,*)'BRhjtautau:'
	    call output_vector(NHIGGS,BRhjtautau)
	    write(*,*)'BRhjWW:'
	    call output_vector(NHIGGS,BRhjWW)
	    write(*,*)'BRhjgamgam:'
	    call output_vector(NHIGGS,BRhjgamgam)	
	      if(.not. useSingH) then 
		write(*,*)'BRhjhihi:'
	  	call output_matrix(NHIGGS,BRhjhihi)
	      endif
            if(useLEP) then
	      write(*,*)'R_hjZ_ee:'
              call output_vector(NHIGGS,R_hjZ_ee)
              if(.not. useSingH) then
		write(*,*)'R_hjhi_ee:'
		call output_matrix(NHIGGS,R_hjhi_ee)
	      endif
            endif
          endif

	  if(inputmode .eq. 2) then
	  	write(*,*)'g_hjbb_sq:'
	  	call output_vector(NHIGGS,g_hjbb_sq)
          	write(*,*)'g_hjtautau_sq:'
          	call output_vector(NHIGGS,g_hjtautau_sq)
          	write(*,*)'g_hjWW_sq:'
          	call output_vector(NHIGGS,g_hjWW_sq)
          	write(*,*)'g_hjZZ_sq:'
          	call output_vector(NHIGGS,g_hjZZ_sq)
          	write(*,*)'g_hjgamgam_sq:'
          	call output_vector(NHIGGS,g_hjgamgam_sq)
          	write(*,*)'g_hjgg_sq:'
          	call output_vector(NHIGGS,g_hjgg_sq)
	  	write(*,*)'g_hjhiZ_sq:'
	  	call output_matrix(NHIGGS,g_hjhiZ_sq)
	  endif

	  if(useTEV) then	  
            if(inputmode .eq. 1) then
	      write(*,*)'R_hj_gg:'
              call output_vector(NHIGGS,R_hj_gg)
	      write(*,*)'R_hj_bb:'
              call output_vector(NHIGGS,R_hj_bb)
	      write(*,*)'R_hjb_bg:'
              call output_vector(NHIGGS,R_hjb_bg)
	      write(*,*)'R_hjWp_udb:'
              call output_vector(NHIGGS,R_hjWp_udb)
	      write(*,*)'R_hjWp_csb:'
              call output_vector(NHIGGS,R_hjWp_csb)
	      write(*,*)'R_hjWm_dub:'
              call output_vector(NHIGGS,R_hjWm_dub)
	      write(*,*)'R_hjWm_scb:'
              call output_vector(NHIGGS,R_hjWm_scb)
	      write(*,*)'R_hjZ_ddb:'
              call output_vector(NHIGGS,R_hjZ_ddb)
	      write(*,*)'R_hjZ_uub:'
              call output_vector(NHIGGS,R_hjZ_uub)
	      write(*,*)'R_hjZ_ssb:'
              call output_vector(NHIGGS,R_hjZ_ssb)
	      write(*,*)'R_hjZ_ccb:'
              call output_vector(NHIGGS,R_hjZ_ccb)
	      write(*,*)'R_hj:Z_bbb'
              call output_vector(NHIGGS,R_hjZ_bbb)
	      write(*,*)'R_sigma_H_VBF:'
              call output_vector(NHIGGS,R_sigma_H_VBF)
            elseif(inputmode .eq. 4) then
              write(*,*)'R_sigma_H:' 
              call output_vector(NHIGGS,R_sigma_H)
              write(*,*)'R_sigma_Hb:'
              call output_vector(NHIGGS,R_sigma_Hb)
              write(*,*)'R_sigma_HW:'
              call output_vector(NHIGGS,R_sigma_HW)
              write(*,*)'R_sigma_HZ:'
              call output_vector(NHIGGS,R_sigma_HZ)
              write(*,*)'R_sigma_H_VBF:'
              call output_vector(NHIGGS,R_sigma_H_VBF)
	    endif
	  endif
#endif

	if(inputmode .eq. 2) then

! translate g_hjhiZ_sq to g_hjhiZ_sq_compact
        k=0
        do j=1,NHIGGS
          do i=1,NHIGGS
            if(i .le. j) then
              k=k+1
              g_hjhiZ_sq_compact(k)=g_hjhiZ_sq(j,i)
            endif
          enddo
        enddo

          call calc_CS_and_BR_from_couplings(NHIGGS,Mh,Gamma_h,
     &          g_hjbb_sq, g_hjtautau_sq, g_hjWW_sq, g_hjZZ_sq,
     &          g_hjgamgam_sq, g_hjgg_sq, g_hjhiZ_sq_compact,
     &          R_hjZ_ee,R_hjhi_ee_compact,
     &          R_hj_gg,R_hj_bb,R_hjb_bg,
     &          R_hjWp_udb,R_hjWp_csb,R_hjWm_dub,R_hjWm_scb,
     &          R_hjZ_ddb,R_hjZ_uub,R_hjZ_ssb,R_hjZ_ccb,R_hjZ_bbb,
     &          R_sigma_H_VBF,
     &          BRhjbbbar,BRhjtautau,BRhjWW,BRhjgamgam)

! translate R_hjhi_ee_compact to R_hjhi_ee
        k=0
        do j=1,NHIGGS
          do i=1,NHIGGS
            if(i .le. j) then
              k=k+1
              R_hjhi_ee(j,i)=R_hjhi_ee_compact(k)
              if(i .ne. j) R_hjhi_ee(i,j)=R_hjhi_ee(j,i)
            endif
          enddo
        enddo


	endif

	if((inputmode .eq. 1) .or. (inputmode .eq. 2)) then
	  call calc_hadCS_from_partCS(NHIGGS, Mh, 
     &  	R_hj_gg, R_hj_bb, R_hjb_bg,
     &          R_hjWp_udb, R_hjWp_csb, R_hjWm_dub, R_hjWm_scb,
     &          R_hjZ_ddb, R_hjZ_uub, R_hjZ_ssb, R_hjZ_ccb, R_hjZ_bbb,
     &		R_sigma_H, R_sigma_Hb, R_sigma_HW, R_sigma_HZ)
	endif

	CMDLINE=.TRUE.
	call check_input_sanity(NHIGGS,CMDLINE,id(1),useLEP,useTEV,useSingH,
     &		Mh,BRhjbbbar,BRhjtautau,BRhjWW,BRhjgamgam,BRhjhihi,sanity_flag)

        experiment_id=0
        if(sanity_flag .eq. 1) then
!tag1
	  if(useLEP) then
! translate BRhjhihi to BRhjhihi_compact
	k=0
        do j=1,NHIGGS
          do i=1,NHIGGS
            if(j .ne. i) then
              k=k+1
              BRhjhihi_compact(k)=BRhjhihi(j,i)
            endif
          enddo
        enddo

! translate R_hjhi_ee to R_hjhi_ee_compact
        k=0
        do j=1,NHIGGS
          do i=1,NHIGGS
            if(i .le. j) then
              k=k+1
              R_hjhi_ee_compact(k)=R_hjhi_ee(j,i)
            endif
          enddo
        enddo

c check LEP bounds
            call check_LEP_bound(NHIGGS,whichexpt,Mh,R_hjZ_ee,R_hjhi_ee_compact,
     &          BRhjbbbar,BRhjtautau,BRhjhihi_compact,BRhjgamgam,
     &          lep_allowed,lep_result_flag,lep_mlratio,lep_obsratio,lep_predratio,
     &		lep_num_H_combined,lep_final_h_comb,lep_model_R_sig_BR)
	  endif


	  if(useTEV) then
c check Tevatron bounds
	    call check_TEV_bound(NHIGGS, whichexpt, Mh, Gamma_h,
     &          R_sigma_H,R_sigma_Hb,R_sigma_HW,R_sigma_HZ,R_sigma_H_VBF,
     &          BRhjbbbar,BRhjtautau,BRhjWW,BRhjgamgam,
     &          tev_allowed,tev_result_flag,tev_mlratio,tev_obsratio,tev_predratio,
     &          tev_num_H_combined,tev_final_h_comb,tev_model_R_sig_BR)
	  endif
       	

	  if(useLEP) then
                experiment_id = 1
                allowed = lep_allowed
                mlratio = lep_mlratio
                obsratio = lep_obsratio
                predratio = lep_predratio
                best_result = lep_result_flag
	 	num_H_combined = lep_num_H_combined
		do i=1,NHIGGS
                        H_combined(i)=lep_final_h_comb(i)
                enddo
          endif
	  if(useTEV) then
                experiment_id = 2
                allowed = tev_allowed
                mlratio = tev_mlratio
                obsratio = tev_obsratio
		predratio = tev_predratio
                best_result = tev_result_flag
		num_H_combined = tev_num_H_combined
c		do i=1,NHIGGS
c			H_combined(i)=tev_final_h_comb(i)
c		enddo
	  endif

	  if(useLEP .and. useTEV) then
c  decide on which ratio to choose from LEP or TEV and set flag: experiment_id
#ifdef VERBOSE
	    write(*,*)'Overall channel with highest sensitivity:'
#endif
            if(lep_predratio .ge. tev_predratio) then
                experiment_id = 1
                allowed = lep_allowed
                mlratio = lep_mlratio
                obsratio = lep_obsratio
                predratio = lep_predratio
                best_result = lep_result_flag
		num_H_combined = lep_num_H_combined
c		do i=1,NHIGGS
c                        H_combined(i)=lep_final_h_comb(i)
c                enddo
            else
                experiment_id = 2
                allowed = tev_allowed
                mlratio = tev_mlratio
                obsratio = tev_obsratio
                predratio = tev_predratio
                best_result = tev_result_flag
		num_H_combined = tev_num_H_combined
c		do i=1,NHIGGS
c			H_combined(i)=tev_final_h_comb(i)
c		enddo
            endif
#ifdef VERBOSE
	        call determine_output_channel_number(NHIGGS,whichexpt,
     &          	experiment_id,mlratio,output_chan)
	        call write_output_channel(6,NHIGGS,whichexpt,output_chan)
        	write(*,*)'   R_sig_BR/S95_predicted =',predratio
        	write(*,*)'   R_sig_BR/S95_observed  =',obsratio
#endif
	  endif

	endif
!tag1
	if(sanity_flag .eq. 0) then
		best_result=-1
		experiment_id=0
		num_H_combined=0
	endif

! 16.1.2009
! If obsratio is zero, than output_chan is set to zero.
! Realization: set experiment_id=0
!	       then determine_output_channel_number will set 
!	       output_chan=0
	if(obsratio .eq. 0d0) then
	  experiment_id=0
	  num_H_combined=0
	endif

	call determine_output_channel_number(NHIGGS,whichexpt,
     &          experiment_id,mlratio,output_chan)

	used_channel(output_chan) = .TRUE.

#ifdef G77VERSION
	if(NHIGGS .eq. 1) write(30,401) id(1),(Mh(i),i=1,NHIGGS),best_result,output_chan,
     & 		obsratio,num_H_combined,(add_dat(j),j=1,columns)
	if(NHIGGS .eq. 2) write(30,402) id(1),(Mh(i),i=1,NHIGGS),best_result,output_chan,
     & 		obsratio,num_H_combined,(add_dat(j),j=1,columns)
	if(NHIGGS .eq. 3) write(30,403) id(1),(Mh(i),i=1,NHIGGS),best_result,output_chan,
     & 		obsratio,num_H_combined,(add_dat(j),j=1,columns)
	if(NHIGGS .eq. 4) write(30,404) id(1),(Mh(i),i=1,NHIGGS),best_result,output_chan,
     & 		obsratio,num_H_combined,(add_dat(j),j=1,columns)
	if(NHIGGS .eq. 5) write(30,405) id(1),(Mh(i),i=1,NHIGGS),best_result,output_chan,
     & 		obsratio,num_H_combined,(add_dat(j),j=1,columns)
	if(NHIGGS .eq. 6) write(30,406) id(1),(Mh(i),i=1,NHIGGS),best_result,output_chan,
     & 		obsratio,num_H_combined,(add_dat(j),j=1,columns)
	if(NHIGGS .eq. 7) write(30,407) id(1),(Mh(i),i=1,NHIGGS),best_result,output_chan,
     & 		obsratio,num_H_combined,(add_dat(j),j=1,columns)
	if(NHIGGS .eq. 8) write(30,408) id(1),(Mh(i),i=1,NHIGGS),best_result,output_chan,
     & 		obsratio,num_H_combined,(add_dat(j),j=1,columns)
	if(NHIGGS .eq. 9) write(30,409) id(1),(Mh(i),i=1,NHIGGS),best_result,output_chan,
     & 		obsratio,num_H_combined,(add_dat(j),j=1,columns)
	if(NHIGGS .eq. 10) write(30,410) id(1),(Mh(i),i=1,NHIGGS),best_result,output_chan,
     & 		obsratio,num_H_combined,(add_dat(j),j=1,columns)
 401	FORMAT (I7,G16.6,I5,I5,G16.6,I5,MAX_ADDITIONAL_COLUMNS G16.6)
 402	FORMAT (I7,2G16.6,I5,I5,G16.6,I5,MAX_ADDITIONAL_COLUMNS G16.6)
 403	FORMAT (I7,3G16.6,I5,I5,G16.6,I5,MAX_ADDITIONAL_COLUMNS G16.6)
 404	FORMAT (I7,4G16.6,I5,I5,G16.6,I5,MAX_ADDITIONAL_COLUMNS G16.6)
 405	FORMAT (I7,5G16.6,I5,I5,G16.6,I5,MAX_ADDITIONAL_COLUMNS G16.6)
 406	FORMAT (I7,6G16.6,I5,I5,G16.6,I5,MAX_ADDITIONAL_COLUMNS G16.6)
 407	FORMAT (I7,7G16.6,I5,I5,G16.6,I5,MAX_ADDITIONAL_COLUMNS G16.6)
 408	FORMAT (I7,8G16.6,I5,I5,G16.6,I5,MAX_ADDITIONAL_COLUMNS G16.6)
 409	FORMAT (I7,9G16.6,I5,I5,G16.6,I5,MAX_ADDITIONAL_COLUMNS G16.6)
 410	FORMAT (I7,10G16.6,I5,I5,G16.6,I5,MAX_ADDITIONAL_COLUMNS G16.6)

#else
	write(30,455) id(1),(Mh(i),i=1,NHIGGS),best_result,output_chan,
     & 		obsratio,num_H_combined,(add_dat(j),j=1,columns)
!,(log2int(H_combined(i)),i=1,NHIGGS)
 455	FORMAT (I7,<NHIGGS>G16.6,I5,I5,G16.6,I5,<columns>G16.6)
! 455	FORMAT (I7,<NHIGGS>G16.6,I5,I5,G16.6,I5,' ',<NHIGGS>I2)
#endif
	endif

	enddo

c write legend

	write(30,*)'# channel numbers used in this file'
	do i=0,LEP_NSCHN*NHIGGS+LEP_NDCHN*NHIGGS**2+TEV_NSCHN*NHIGGS
	  if(used_channel(i)) 
     &        call write_output_channel(30,NHIGGS,whichexpt,i)
	enddo
		
	close(30)

	end



*****************************************************************************
        subroutine check_input_sanity(NHIGGS,CMDLINE,id,useLEP,useTEV,useSingH,
     &		Mh,BRhjbbbar,BRhjtautau,BRhjWW,
     & 		BRhjgamgam,BRhjhihi,sanity_flag)
*****************************************************************************
* check_input_sanity ignores all BRs which are not relevant to the experiment.

	implicit none
	logical useLEP,useTEV,useSingH,CMDLINE
	integer NHIGGS,sanity_flag,i,j,id
	double precision sumBR(NHIGGS)
	double precision Mh(NHIGGS),BRhjbbbar(NHIGGS),BRhjtautau(NHIGGS)
        double precision BRhjWW(NHIGGS),BRhjgamgam(NHIGGS),BRhjhihi(NHIGGS,NHIGGS)
	double precision fuzziness

	fuzziness=0.05d0

	sanity_flag=1
	
	do j=1,NHIGGS
	  sumBR(j)=0d0
	  if(Mh(j) .le. 0d0) then
	    if(CMDLINE) write(*,7788) j,id
	    if(.not. CMDLINE) write(*,7789) j
 7788	    FORMAT('WARNING: zero or negative mass of Higgs #',I3,' in input line',I6)
 7789	    FORMAT('WARNING: zero or negative mass of Higgs #',I3)
            sanity_flag=0
	  endif
	enddo


	do j=1,NHIGGS

!>> check_input_sanity ignores all BRs which are not relevant 
!   to the chosen experiment.
!	  if(useLEP) then 
!		sumBR(j)=BRhjbbbar(j)+BRhjtautau(j)
!	    if(.not. useSingH) then
!	      do i=1,NHIGGS
!                sumBR(j)=sumBR(j)+BRhjhihi(j,i)
!	      enddo
!	    endif
!	  endif
!	  if(useTEV) then
!		sumBR(j)=BRhjbbbar(j)+BRhjtautau(j)
!     &			     +BRhjWW(j)+BRhjgamgam(j)
!	    if(useLEP .and. (.not. useSingH)) then
!	      do i=1,NHIGGS
!                sumBR(j)=sumBR(j)+BRhjhihi(j,i)
!              enddo
!	    endif
!    	  endif
!<<

!>> check_input_sanity takes all input into account
	  sumBR(j)=BRhjbbbar(j)+BRhjtautau(j)+BRhjWW(j)+BRhjgamgam(j)
	  do i=1,NHIGGS
            sumBR(j)=sumBR(j)+BRhjhihi(j,i)
          enddo
!<<

!>> 21.01.2009 fuzziness construction switched off
!          if(sumBR(j) .gt. 1d0+fuzziness) then
!	     if(CMDLINE) write(*,7790) j,id
!	     if(.not. CMDLINE) write(*,7791) j
! 7790	     FORMAT('WARNING: nonsensical sum of BRs for Higgs #',I3,' in input line',I6)
! 7791	     FORMAT('WARNING: nonsensical sum of BRs for Higgs #',I3)
!	     write(*,*)'BRhjbbbar(',j,')=',BRhjbbbar(j)
!            write(*,*)'BRhjtautau(',j,')=',BRhjtautau(j)
!	     write(*,*)'BRhjWW(',j,')=',BRhjWW(j)
!             write(*,*)'BRhjgamgam(',j,')=',BRhjgamgam(j)
!	     write(*,*)'sumBR(',j,')=',sumBR(j)
!	     sanity_flag=0
!	  endif
!<<
          if(sumBR(j) .gt. 1d0) then
c	     if(CMDLINE) write(*,7790) sumBR(j),j,id
c	     if(.not. CMDLINE) write(*,7791) sumBR(j),j
 7790	     FORMAT('WARNING: sum of BRs larger than 1 (',G10.4,') for Higgs #',I3,' in input line',I6)
 7791	     FORMAT('WARNING: sum of BRs larger than 1 (',G10.4,') for Higgs #',I3)
! 21.01.2009 currently only a warning is given
!	     sanity_flag=0
	  endif
	enddo

	end

*****************************************************************************
	subroutine parse_commandline(NHIGGS,whichexpt,
     &           inputmode,fileprefix,prefix)
*****************************************************************************
	implicit none 

	integer iargc, lnblnk
!        external iargc, lnblnk

	integer narg,i
	character*100 barearg,arg(4)
        character*5 whichexpt
	
	integer inputmode,NHIGGS
!	logical useTEV,useLEP,useSingHonly,prefix
	logical prefix
	character*100 fileprefix

	narg = iargc()
	
	if((narg .lt. 3) .or. (narg .gt. 4)) goto 100
	do i=1,narg,1
	  call getarg(i, barearg)
          arg(i)=barearg(:lnblnk(barearg))
	enddo

	if(arg(1). eq. 'onlyT') then
		whichexpt=arg(1)
	elseif(arg(1). eq. 'onlyL') then
		whichexpt=arg(1)
	elseif(arg(1). eq. 'LandT') then
		whichexpt=arg(1)
	elseif(arg(1). eq. 'singH') then
		whichexpt=arg(1)
	else
		goto 100
	endif
	
	if(arg(2) .eq. 'part') then
		inputmode=1
	elseif(arg(2) .eq. 'effC') then
		inputmode=2
!	elseif(arg(2) .eq. 'SLHA') then
!		inputmode=3
	elseif(arg(2) .eq. 'hadr') then
		inputmode=4
	else
		goto 100
	endif

#ifdef G77VERSION
	read(arg(3),'(I20)',err=100) NHIGGS
#else
	read(arg(3),'(I)',err=100) NHIGGS
#endif

	if((NHIGGS .lt. 1)	.or. (NHIGGS .gt. 10)) goto 100

	
	if(narg .eq. 3) then
		fileprefix='                    '
		prefix=.FALSE.
	else
		fileprefix=arg(4)
                prefix=.TRUE.
	endif



	return
	
 100 	continue

c>> check
c	do i=1,narg
c	  write(*,*)'arg(',i,')=',arg(i)
c	enddo
c        write(*,*)'NHIGGS=',NHIGGS
c	write(*,*)'inputmode=',inputmode
c	write(*,*)'fileprefix=',fileprefix
c        write(*,*)'prefix=',prefix
c<<

	write(*,*)'wrong input  !'
	write(*,*)'HiggsBounds usage :'
	write(*,*)'  HiggsBounds <analyses> <input mode> <nHiggs> [<fileprefix>]'
	write(*,*)'with'
	write(*,*)'     <analyses>   : LandT (LEP and Tevatron analyses)'
	write(*,*)'                  : onlyT (only Tevatron analyses)'
	write(*,*)'                  : onlyL (only LEP analyses)'
	write(*,*)'                  : singH (only single Higgs analyses)'
	write(*,*)'     <input mode> : part (BRs and CS ratios (partonic CS ratios for Tevatron))'
	write(*,*)'                  : hadr (BRs and CS ratios (hadronic CS ratios for Tevatron))'
	write(*,*)'                  : effC (effective couplings)'	
!	write(*,*)'                  : SLHA (SLHA format input)'
	write(*,*)'     <nHiggs>     : 1 to 10 (number of neutral Higgs bosons)'
	write(*,*)'     <fileprefix> : prefix for input files (optional,'
	write(*,*)'                           can also be a subdirectory)'

	stop

	end




*****************************************************************************
        subroutine open_output_file(file_id,fileprefix,prefix,NHIGGS,whichexpt)
*****************************************************************************
        implicit none

        integer iargc, lnblnk
!        external iargc, lnblnk

        integer narg,i,file_id,NHIGGS
        character*100 barearg,arg(4)
	character*5 whichexpt
        integer inputmode
        logical prefix
        character*100 fileprefix

	if(prefix) then
	  open(file_id,file=fileprefix(:lnblnk(fileprefix))//'HiggsBounds_results.dat')
	  open(file_id+1,file=fileprefix(:lnblnk(fileprefix))//'Key.dat')
	else
	  open(file_id,file='HiggsBounds_results.dat')
	  open(file_id+1,file='Key.dat')
	endif

	call write_table_key(file_id+1,NHIGGS,whichexpt)

	close(file_id+1)

        narg = iargc()
	do i=1,narg,1
          call getarg(i, barearg)
          arg(i)=barearg(:lnblnk(barearg))
        enddo
	
        write(file_id,*)'# File generated with HiggsBounds 1.1.0 (FORTRAN77 version)'
	write(file_id,*)'# HiggsBounds command used:'
	write(file_id,*)'#   HiggsBounds ',(arg(i)(:lnblnk(arg(i)))//' ',i=1,narg)
        write(file_id,*)'#'
	write(file_id,*)'# column description:'
	write(file_id,*)'#   0 : line id of input'
	write(file_id,459) NHIGGS
 459	FORMAT(' #   i : Higgs boson masses Mh(i) (i=1,...,',I3,')')
	write(file_id,460) NHIGGS+1
 460	FORMAT(' # ',I3,' : scenario allowed flag (1: allowed, 0: excluded, -1: unphysical)')
	write(file_id,461) NHIGGS+2
 461	FORMAT(' # ',I3,' : most sensitive channel (see files: Key_LEP-tables, Key_TEV-tables)')
!	write(file_id,462) NHIGGS+3
! 462	FORMAT(' # ',I3,' : experiment with best channel (1: LEP, 2: Tevatron)')
	write(file_id,463) NHIGGS+3
 463	FORMAT(' # ',I3,' : ratio [sig x BR]_model/[sig x BR]_limit (<1: allowed, >1: excluded)')
	write(file_id,464) NHIGGS+4
 464	FORMAT(' # ',I3,' : number of Higgs bosons combined in most sensitive channel')
	write(file_id,*)'# columns :'
#ifdef G77VERSION
	if(NHIGGS .eq. 1) write(file_id,401) (i,i=0,4+NHIGGS)
	if(NHIGGS .eq. 2) write(file_id,402) (i,i=0,4+NHIGGS)
	if(NHIGGS .eq. 3) write(file_id,403) (i,i=0,4+NHIGGS)
	if(NHIGGS .eq. 4) write(file_id,404) (i,i=0,4+NHIGGS)
	if(NHIGGS .eq. 5) write(file_id,405) (i,i=0,4+NHIGGS)
	if(NHIGGS .eq. 6) write(file_id,406) (i,i=0,4+NHIGGS)
	if(NHIGGS .eq. 7) write(file_id,407) (i,i=0,4+NHIGGS)
	if(NHIGGS .eq. 8) write(file_id,408) (i,i=0,4+NHIGGS)
	if(NHIGGS .eq. 9) write(file_id,409) (i,i=0,4+NHIGGS)
	if(NHIGGS .eq. 10) write(file_id,410) (i,i=0,4+NHIGGS)
 401	FORMAT (' #',I5,(I9,'       '),I5,I5,I9,'       ',I5)
 402	FORMAT (' #',I5,2(I9,'       '),I5,I5,I9,'       ',I5)
 403	FORMAT (' #',I5,3(I9,'       '),I5,I5,I9,'       ',I5)
 404	FORMAT (' #',I5,4(I9,'       '),I5,I5,I9,'       ',I5)
 405	FORMAT (' #',I5,5(I9,'       '),I5,I5,I9,'       ',I5)
 406	FORMAT (' #',I5,6(I9,'       '),I5,I5,I9,'       ',I5)
 407	FORMAT (' #',I5,7(I9,'       '),I5,I5,I9,'       ',I5)
 408	FORMAT (' #',I5,8(I9,'       '),I5,I5,I9,'       ',I5)
 409	FORMAT (' #',I5,9(I9,'       '),I5,I5,I9,'       ',I5)
 410	FORMAT (' #',I5,10(I9,'       '),I5,I5,I9,'       ',I5)
#else
	write(file_id,457) (i,i=0,4+NHIGGS)
 457	FORMAT (' #',I5,<NHIGGS>(I9,'       '),I5,I5,I9,'       ',I5)
#endif
        write(file_id,*)'#'


	end


*****************************************************************************
	subroutine determine_topology_flags(whichexpt,
     &		useTEV,useLEP,useSingH)
*****************************************************************************
	implicit none 

	character*5 whichexpt

	integer lnblnk
!        external lnblnk

	logical useTEV,useLEP,useSingH


!	whichexpt=whichexpt(::lnblnk(whichexpt))

	useSingH = .FALSE.

	if(whichexpt .eq. 'onlyT') then
		useTEV=.TRUE.
		useLEP=.FALSE.
	elseif(whichexpt .eq. 'onlyL') then
		useTEV=.FALSE.
                useLEP=.TRUE.
	elseif(whichexpt .eq. 'LandT') then
                useTEV=.TRUE.
                useLEP=.TRUE.
	elseif(whichexpt .eq. 'singH') then
               	useTEV=.TRUE.
               	useLEP=.TRUE.	
	       	useSingH = .TRUE.
	endif
	
	end



*****************************************************************************
        subroutine determine_output_channel_number(NHIGGS,whichexpt,
     &		experiment_id,exp_chan,output_chan)
*****************************************************************************
        implicit none
	character*5 whichexpt
	integer NHIGGS,experiment_id,exp_chan,output_chan,tev_offset

	if(experiment_id .eq. 0) then
		output_chan = 0
	else
	if(whichexpt .eq. 'onlyT') then
		tev_offset=0
		output_chan=exp_chan
c>> check
		if(experiment_id .eq. 1) then
		write(*,*)'error: LEP channel but setting is: onlyT'
		endif
c<< check
	elseif(whichexpt .eq. 'onlyL') then
		output_chan=exp_chan
c>> check
		if(experiment_id .eq. 2) then
		write(*,*)'error: TEV channel but setting is: onlyL'
		endif
c<< check		
	elseif(whichexpt .eq. 'LandT') then
		tev_offset=LEP_NSCHN*NHIGGS+LEP_NDCHN*NHIGGS**2
		if(experiment_id .eq. 1) then
			output_chan=exp_chan
		elseif(experiment_id .eq. 2) then
			output_chan=exp_chan+tev_offset
		endif
	elseif(whichexpt .eq. 'singH') then
		tev_offset=LEP_NSCHN*NHIGGS
		if(experiment_id .eq. 1) then
			output_chan=exp_chan
		elseif(experiment_id .eq. 2) then
			output_chan=exp_chan+tev_offset
		endif
	endif
	endif

	end




#include "analyses-descriptions.h"

