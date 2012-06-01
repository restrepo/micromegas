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

#include "tables-nhmax.F"
#include "LEP-bound.F"
#include "TEV-bound.F"


*************************************************************************
	subroutine run_HiggsBounds_effC(NHIGGS,Mh,GammaTotal, 
     & g2hjbb,g2hjtautau,g2hjWW,g2hjZZ, 
     & g2hjgaga,g2hjgg,g2hjhiZ, 
     & BR_hjhihi, 
     & HBresult,chan, 
     & obsratio, ncombined ) 
*************************************************************************
	implicit none
! input:
	integer NHIGGS,i,j,k
	double precision Mh(NHIGGS),GammaTotal(NHIGGS)
	double precision g2hjbb(NHIGGS),g2hjtautau(NHIGGS),
     &		g2hjWW(NHIGGS),g2hjZZ(NHIGGS),g2hjgaga(NHIGGS),g2hjgg(NHIGGS),
     &		g2hjhiZ(NHIGGS,NHIGGS),BR_hjhihi(NHIGGS,NHIGGS)
! output:
	integer HBresult,chan,ncombined
	double precision obsratio
! internal:
!  Tevatron partonic CS
	double precision CS_tev_gg_hj_ratio(NHIGGS),CS_tev_bb_hj_ratio(NHIGGS),
     & 		CS_tev_bg_hjb_ratio(NHIGGS),
     & 		CS_tev_ud_hjWp_ratio(NHIGGS),CS_tev_cs_hjWp_ratio(NHIGGS),
     & 		CS_tev_ud_hjWm_ratio(NHIGGS),CS_tev_cs_hjWm_ratio(NHIGGS),
     & 		CS_tev_dd_hjZ_ratio(NHIGGS),CS_tev_uu_hjZ_ratio(NHIGGS),
     & 		CS_tev_ss_hjZ_ratio(NHIGGS),CS_tev_cc_hjZ_ratio(NHIGGS),
     & 		CS_tev_bb_hjZ_ratio(NHIGGS)
!  LEP CS ratios & Tevatron hadronic CS
        double precision CS_lep_hjZ_ratio(NHIGGS),
     &		CS_lep_hjhi_ratio(NHIGGS,NHIGGS)

	double precision g2hjhiZ_compact((NHIGGS*(NHIGGS+1))/2)
	double precision CS_lep_hjhi_ratio_compact((NHIGGS*(NHIGGS+1))/2)

	double precision CS_tev_pp_hj_ratio(NHIGGS),
     &		CS_tev_pp_hjb_ratio(NHIGGS),CS_tev_pp_hjW_ratio(NHIGGS),
     &		CS_tev_pp_hjZ_ratio(NHIGGS),CS_tev_pp_vbf_ratio(NHIGGS)
!  BRs
       	double precision BR_hjbb(NHIGGS),BR_hjtautau(NHIGGS),BR_hjWW(NHIGGS),
     &		BR_hjgaga(NHIGGS)

! translate g2hjhiZ to g2hjhiZ_compact
        k=0
        do j=1,NHIGGS
          do i=1,NHIGGS
            if(i .le. j) then
              k=k+1
              g2hjhiZ_compact(k)=g2hjhiZ(j,i)
            endif
          enddo
        enddo


	  call calc_CS_and_BR_from_couplings(NHIGGS,Mh,GammaTotal,
     &		g2hjbb, g2hjtautau, g2hjWW, g2hjZZ,
     &          g2hjgaga, g2hjgg, g2hjhiZ_compact,
     &		CS_lep_hjZ_ratio,CS_lep_hjhi_ratio_compact,
     &		CS_tev_gg_hj_ratio,CS_tev_bb_hj_ratio,CS_tev_bg_hjb_ratio,
     &		CS_tev_ud_hjWp_ratio,CS_tev_cs_hjWp_ratio,
     &		CS_tev_ud_hjWm_ratio,CS_tev_cs_hjWm_ratio,
     &          CS_tev_dd_hjZ_ratio,CS_tev_uu_hjZ_ratio,CS_tev_ss_hjZ_ratio,
     &		CS_tev_cc_hjZ_ratio,CS_tev_bb_hjZ_ratio,
     &          CS_tev_pp_vbf_ratio,
     &		BR_hjbb,BR_hjtautau,BR_hjWW,BR_hjgaga)

! translate CS_lep_hjhi_ratio_compact to CS_lep_hjhi_ratio
        k=0
        do j=1,NHIGGS
          do i=1,NHIGGS
            if(i .le. j) then
              k=k+1
              CS_lep_hjhi_ratio(j,i)=CS_lep_hjhi_ratio_compact(k)
              if(i .ne. j) CS_lep_hjhi_ratio(i,j)=CS_lep_hjhi_ratio(j,i)
            endif
          enddo
        enddo


	call calc_hadCS_from_partCS(NHIGGS, Mh,
     &          CS_tev_gg_hj_ratio, CS_tev_bb_hj_ratio, CS_tev_bg_hjb_ratio,
     &		CS_tev_ud_hjWp_ratio, CS_tev_cs_hjWp_ratio, 
     &		CS_tev_ud_hjWm_ratio, CS_tev_cs_hjWm_ratio,
     &          CS_tev_dd_hjZ_ratio, CS_tev_uu_hjZ_ratio, CS_tev_ss_hjZ_ratio, 
     &		CS_tev_cc_hjZ_ratio, CS_tev_bb_hjZ_ratio,
     &          CS_tev_pp_hj_ratio, CS_tev_pp_hjb_ratio, 
     &		CS_tev_pp_hjW_ratio, CS_tev_pp_hjZ_ratio)

	call run_HiggsBounds_hadr(NHIGGS,Mh,
     & CS_lep_hjZ_ratio,CS_lep_hjhi_ratio,
     & CS_tev_pp_hj_ratio,CS_tev_pp_hjb_ratio,
     & CS_tev_pp_hjW_ratio,CS_tev_pp_hjZ_ratio,
     & CS_tev_pp_vbf_ratio,
     & BR_hjbb,BR_hjtautau,
     & BR_hjWW,BR_hjgaga,
     & BR_hjhihi,
     & HBresult,chan,obsratio, ncombined )

	end

*************************************************************************
	subroutine run_HiggsBounds_part(NHIGGS,Mh, 
     & CS_lep_hjZ_ratio, 
     & CS_lep_hjhi_ratio, 
     & CS_tev_gg_hj_ratio,CS_tev_bb_hj_ratio, 
     & CS_tev_bg_hjb_ratio, 
     & CS_tev_ud_hjWp_ratio,CS_tev_cs_hjWp_ratio, 
     & CS_tev_ud_hjWm_ratio,CS_tev_cs_hjWm_ratio, 
     & CS_tev_dd_hjZ_ratio,CS_tev_uu_hjZ_ratio, 
     & CS_tev_ss_hjZ_ratio,CS_tev_cc_hjZ_ratio, 
     & CS_tev_bb_hjZ_ratio, 
     & CS_tev_pp_vbf_ratio, 
     & BR_hjbb,BR_hjtautau, 
     & BR_hjWW,BR_hjgaga, 
     & BR_hjhihi, 
     & HBresult,chan, 
     & obsratio, ncombined ) 
*************************************************************************
	implicit none
! input:
	integer NHIGGS
	double precision Mh(NHIGGS),GammaTotal(NHIGGS)
        double precision CS_lep_hjZ_ratio(NHIGGS),
     &		CS_lep_hjhi_ratio(NHIGGS,NHIGGS)
	double precision CS_tev_gg_hj_ratio(NHIGGS),CS_tev_bb_hj_ratio(NHIGGS),
     & 		CS_tev_bg_hjb_ratio(NHIGGS),
     & 		CS_tev_ud_hjWp_ratio(NHIGGS),CS_tev_cs_hjWp_ratio(NHIGGS),
     & 		CS_tev_ud_hjWm_ratio(NHIGGS),CS_tev_cs_hjWm_ratio(NHIGGS),
     & 		CS_tev_dd_hjZ_ratio(NHIGGS),CS_tev_uu_hjZ_ratio(NHIGGS),
     & 		CS_tev_ss_hjZ_ratio(NHIGGS),CS_tev_cc_hjZ_ratio(NHIGGS),
     & 		CS_tev_bb_hjZ_ratio(NHIGGS)
       	double precision BR_hjbb(NHIGGS),BR_hjtautau(NHIGGS),BR_hjWW(NHIGGS),
     &		BR_hjgaga(NHIGGS),BR_hjhihi(NHIGGS,NHIGGS)
! output:
	integer HBresult,chan,ncombined
	double precision obsratio
! internal:
!  Tevatron hadronic CS
	double precision CS_tev_pp_hj_ratio(NHIGGS),
     &		CS_tev_pp_hjb_ratio(NHIGGS),CS_tev_pp_hjW_ratio(NHIGGS),
     &		CS_tev_pp_hjZ_ratio(NHIGGS),CS_tev_pp_vbf_ratio(NHIGGS)


	call calc_hadCS_from_partCS(NHIGGS, Mh,
     &          CS_tev_gg_hj_ratio, CS_tev_bb_hj_ratio, CS_tev_bg_hjb_ratio,
     &		CS_tev_ud_hjWp_ratio, CS_tev_cs_hjWp_ratio, 
     &		CS_tev_ud_hjWm_ratio, CS_tev_cs_hjWm_ratio,
     &          CS_tev_dd_hjZ_ratio, CS_tev_uu_hjZ_ratio, CS_tev_ss_hjZ_ratio, 
     &		CS_tev_cc_hjZ_ratio, CS_tev_bb_hjZ_ratio,
     &          CS_tev_pp_hj_ratio, CS_tev_pp_hjb_ratio, 
     &		CS_tev_pp_hjW_ratio, CS_tev_pp_hjZ_ratio)



	call run_HiggsBounds_hadr(NHIGGS,Mh,
     & 		CS_lep_hjZ_ratio,CS_lep_hjhi_ratio,
     & 		CS_tev_pp_hj_ratio,CS_tev_pp_hjb_ratio,
     & 		CS_tev_pp_hjW_ratio,CS_tev_pp_hjZ_ratio,
     & 		CS_tev_pp_vbf_ratio,
     & 		BR_hjbb,BR_hjtautau,
     & 		BR_hjWW,BR_hjgaga,
     & 		BR_hjhihi,
     & 		HBresult,chan,obsratio, ncombined )

	end

*************************************************************************
	subroutine run_HiggsBounds_hadr(NHIGGS,Mh, 
     & CS_lep_hjZ_ratio,CS_lep_hjhi_ratio, 
     & CS_tev_pp_hj_ratio,CS_tev_pp_hjb_ratio, 
     & CS_tev_pp_hjW_ratio,CS_tev_pp_hjZ_ratio, 
     & CS_tev_pp_vbf_ratio, 
     & BR_hjbb,BR_hjtautau, 
     & BR_hjWW,BR_hjgaga, 
     & BR_hjhihi, 
     & HBresult,chan,obsratio,ncombined ) 
*************************************************************************
	implicit none
! input:
	integer NHIGGS
	character*5 whichexpt
        common /chosen_analyses/ whichexpt

	double precision Mh(NHIGGS),GammaTotal(NHIGGS)
        double precision CS_lep_hjZ_ratio(NHIGGS),
     &		CS_lep_hjhi_ratio(NHIGGS,NHIGGS)
	double precision CS_tev_pp_hj_ratio(NHIGGS),
     &		CS_tev_pp_hjb_ratio(NHIGGS),CS_tev_pp_hjW_ratio(NHIGGS),
     &		CS_tev_pp_hjZ_ratio(NHIGGS),CS_tev_pp_vbf_ratio(NHIGGS)
       	double precision BR_hjbb(NHIGGS),BR_hjtautau(NHIGGS),BR_hjWW(NHIGGS),
     &		BR_hjgaga(NHIGGS),BR_hjhihi(NHIGGS,NHIGGS)
! output:
	integer HBresult,chan,ncombined
	double precision obsratio
! internal:
        double precision BR_hjhihi_compact(NHIGGS**2-NHIGGS)
        double precision CS_lep_hjhi_ratio_compact((NHIGGS*(NHIGGS+1))/2)

	logical useLEP,useTEV,useSingH,CMDLINE	
	integer iddummy,sanity_flag,i,j,k,output_chan

        logical lep_allowed,lep_final_h_comb(NHIGGS)
        integer lep_mlratio
        integer lep_result_flag,lep_num_H_combined,lep_sim_mass_sum(NHIGGS)
        double precision lep_obsratio,lep_predratio,lep_model_R_sig_BR

        logical tev_allowed,tev_final_h_comb(NHIGGS)
        integer tev_mlratio
        integer tev_result_flag,tev_num_H_combined
        double precision tev_obsratio,tev_predratio,tev_model_R_sig_BR

	logical allowed,H_combined(NHIGGS)
        integer mlratio,experiment_id
c experiment_id : 1 : LEP, 2 : Tevatron
        integer sim_mass_sum(NHIGGS)
        double precision predratio,model_R_sig_BR
#include "tables-tev_cb.h"
#include "tables-lep_cb.h"



	call determine_topology_flags(whichexpt,useTEV,useLEP,useSingH)



	iddummy=0
	CMDLINE=.FALSE.
	call check_input_sanity(NHIGGS,CMDLINE,iddummy,useLEP,useTEV,useSingH,
     &		Mh,BR_hjbb,BR_hjtautau,BR_hjWW,BR_hjgaga,BR_hjhihi,
     &		sanity_flag)



        if(sanity_flag .eq. 1) then
!tag1
	  if(useLEP) then
! translate BR_hjhihi to BR_hjhihi_compact
	k=0
        do j=1,NHIGGS
          do i=1,NHIGGS
            if(j .ne. i) then
              k=k+1
              BR_hjhihi_compact(k)=BR_hjhihi(j,i)
            endif
          enddo
        enddo

! translate CS_lep_hjhi_ratio to CS_lep_hjhi_ratio_compact
        k=0
        do j=1,NHIGGS
          do i=1,NHIGGS
            if(i .le. j) then
              k=k+1
              CS_lep_hjhi_ratio_compact(k)=CS_lep_hjhi_ratio(j,i)
            endif
          enddo
        enddo

c check LEP bounds
            call check_LEP_bound(NHIGGS, whichexpt, Mh,
     &		CS_lep_hjZ_ratio,CS_lep_hjhi_ratio_compact,
     &          BR_hjbb,BR_hjtautau,BR_hjhihi_compact,BR_hjgaga,
     &          lep_allowed,lep_result_flag,lep_mlratio,lep_obsratio,lep_predratio,
     &		lep_num_H_combined,lep_final_h_comb,lep_model_R_sig_BR)
	  endif



	  if(useTEV) then
c check Tevatron bounds
	    call check_TEV_bound(NHIGGS, whichexpt, Mh, GammaTotal,
     &          CS_tev_pp_hj_ratio,CS_tev_pp_hjb_ratio,CS_tev_pp_hjW_ratio,
     &		CS_tev_pp_hjZ_ratio,CS_tev_pp_vbf_ratio,
     &          BR_hjbb,BR_hjtautau,BR_hjWW,BR_hjgaga,
     &          tev_allowed,tev_result_flag,tev_mlratio,tev_obsratio,
     &		tev_predratio,tev_num_H_combined,tev_final_h_comb,
     &		tev_model_R_sig_BR)



	  endif

	  if(useLEP) then
                experiment_id = 1
                allowed = lep_allowed
                mlratio = lep_mlratio
                obsratio = lep_obsratio
                predratio = lep_predratio
                HBresult = lep_result_flag
	 	ncombined = lep_num_H_combined
          endif
	  if(useTEV) then
                experiment_id = 2
                allowed = tev_allowed
                mlratio = tev_mlratio
                obsratio = tev_obsratio
		predratio = tev_predratio
                HBresult = tev_result_flag
		ncombined = tev_num_H_combined
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
                HBresult = lep_result_flag
		ncombined = lep_num_H_combined
            else
                experiment_id = 2
                allowed = tev_allowed
                mlratio = tev_mlratio
                obsratio = tev_obsratio
                predratio = tev_predratio
                HBresult = tev_result_flag
		ncombined = tev_num_H_combined
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
		HBresult=-1
		experiment_id=0
		ncombined = 0
	endif

! If obsratio is zero, than output_chan is set to zero.
! Realization: set experiment_id=0
!	       then determine_output_channel_number will set 
!	       output_chan=0
	if(obsratio .eq. 0d0) then
	  experiment_id=0
	  ncombined=0
	endif

	call determine_output_channel_number(NHIGGS,whichexpt,
     &          experiment_id,mlratio,chan)

	end


***************************************************
        subroutine initialize_HiggsBounds(NHIGGS,whichexpt)
***************************************************
        implicit none
        integer NHIGGS
	character*5 whichexpt,whichexpt_val
	common /chosen_analyses/ whichexpt_val

	whichexpt_val=whichexpt

	call write_HiggsBounds_info()

        call initialize_SM_decay_data()
        call initialize_TEV_tables(NHIGGS)
        call initialize_LEP_tables(NHIGGS)

	open(87,file='Key.dat')
	call write_table_key(87,NHIGGS,whichexpt)
	close(87)

        end

***************************************************
        subroutine finish_HiggsBounds()
***************************************************
        implicit none


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




*****************************************************************************
	subroutine determine_topology_flags(whichexpt,
     &		useTEV,useLEP,useSingH)
*****************************************************************************
	implicit none 

	character*5 whichexpt

	integer lnblnk
        external lnblnk

	logical useTEV,useLEP,useSingH


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
	     if(CMDLINE) write(*,7790) sumBR(j),j,id
	     if(.not. CMDLINE) write(*,7791) sumBR(j),j
 7790	     FORMAT('WARNING: sum of BRs larger than 1 (',G10.4,') for Higgs #',I3,' in input line',I6)
 7791	     FORMAT('WARNING: sum of BRs larger than 1 (',G10.4,') for Higgs #',I3)
! 21.01.2009 currently only a warning is given
!	     sanity_flag=0
	  endif
	enddo

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
	elseif(whichexpt .eq. 'onlyL') then
		output_chan=exp_chan
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

