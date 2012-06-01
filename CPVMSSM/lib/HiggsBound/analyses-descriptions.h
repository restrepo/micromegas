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

*****************************************************************************
        subroutine det_chan_n_expid_from_output_chan(NHIGGS,whichexpt,
     &		output_chan,experiment_id,exp_chan)
*****************************************************************************
        implicit none
	character*5 whichexpt
	integer NHIGGS,experiment_id,exp_chan,output_chan,tev_offset


	if(output_chan .eq. 0) then
		experiment_id = 0
	else
	if(whichexpt .eq. 'onlyT') then
		tev_offset=0
		experiment_id=2
		exp_chan=output_chan
	elseif(whichexpt .eq. 'onlyL') then
		experiment_id=1
		exp_chan=output_chan
	elseif(whichexpt .eq. 'LandT') then
		tev_offset=LEP_NSCHN*NHIGGS+LEP_NDCHN*NHIGGS**2
		if(output_chan .le. tev_offset) then
			experiment_id=1
			exp_chan=output_chan
		elseif(output_chan .gt. tev_offset) then
			experiment_id=2
			exp_chan=output_chan-tev_offset
		endif
	elseif(whichexpt .eq. 'singH') then
		tev_offset=LEP_NSCHN*NHIGGS
		if(output_chan .le. tev_offset) then
			experiment_id=1
			exp_chan=output_chan
		elseif(output_chan .gt. tev_offset) then
			experiment_id=2
			exp_chan=output_chan-tev_offset
		endif
	endif
	endif

	end



!******************************************************************
	subroutine write_table_key(file_id,NHIGGS,whichexpt)
!******************************************************************
	implicit none
	character*5 whichexpt
	logical useTEV,useLEP,useSingH
        integer NHIGGS,file_id
        integer output_chan,channel,number_channels

#include "tables-lep_cb.h"
#include "tables-tev_cb.h"

	call determine_topology_flags(whichexpt,
     &          useTEV,useLEP,useSingH)

	number_channels = 0

	if(useLEP) then
	  if(useSingH) then
		number_channels = LEP_NSCHN*NHIGGS
	  else
		number_channels = LEP_NSCHN*NHIGGS+LEP_NDCHN*NHIGGS**2
	  endif
	endif

	if(useTEV) then
	  number_channels = number_channels 
     &		+ TEV_NSCHN*NHIGGS+TEV_NDCHN*NHIGGS**2
	endif

	do output_chan=0,number_channels
	  call write_output_channel(file_id,NHIGGS,whichexpt,
     &                                  output_chan)
	enddo

	end



!******************************************************************
	subroutine write_output_channel(file_id,NHIGGS,whichexpt,
     &					output_chan)
!******************************************************************
	implicit none
	character*5 whichexpt
	character*4 prefix
        integer NHIGGS,file_id,m,n,i,j,findi,findj,TEV_reaction_number
        integer LEP_reaction_number,number_LEP_channels
	integer output_chan,channel,experiment_id	

#include "tables-lep_cb.h"
#include "tables-tev_cb.h"

	if(file_id .eq. 6) then	
		prefix='    '	
	else
		prefix=' #  '
	endif

 	call det_chan_n_expid_from_output_chan(NHIGGS,whichexpt,
     &          output_chan,experiment_id,channel)

! 16.1.2009
	if(output_chan .eq. 0) then
	write(file_id,90) prefix
 90	FORMAT(A4,'  0: no analysis applies')
	endif

	if(experiment_id .eq. 1) then
	  n=LEP_reaction_number(NHIGGS,channel)

	  if(channel .le. LEP_NSCHN*NHIGGS) then
	    if(tlist_lep(channel) .eq. LEPhZ_bbZ) then
	      write(file_id,142) prefix,output_chan,n,n
 142          FORMAT (A4,I3,':  e+ e- -> H(',I1,') Z, H(',I1,') -> b b')
            elseif(tlist_lep(channel) .eq. LEPhZ_tautauZ)  then
	      write(file_id,143) prefix,output_chan,n,n
 143          FORMAT (A4,I3,':  e+ e- -> H(',I1,') Z, H(',I1,') -> tau tau')
            elseif(tlist_lep(channel) .eq. LEPhZ_anyZ)  then
	      write(file_id,1400) prefix,output_chan,n,n
 1400         FORMAT (A4,I3,':  e+ e- -> H(',I1,') Z, H(',I1,') -> anything (LEP, EPJC 27(2003)311)')
            elseif(tlist_lep(channel) .eq. LEPhZ_gamgamZ)  then
	      write(file_id,1401) prefix,output_chan,n,n
 1401         FORMAT (A4,I3,':  e+ e- -> H(',I1,') Z, H(',I1,') -> gamma gamma (LEP, LHWG note 2002-02)')
	    endif
	  else
	    i=findi(NHIGGS,n)
	    j=findj(NHIGGS,n)
	      if(tlist_lep(channel) .eq. LEPh2Z_h1h1Z_4bZ) then
		write(file_id,150) prefix,output_chan,j,j,i,i,i
 150 		FORMAT (A4,I3,':  e+ e- -> H(',I1,') Z, H(',I1,') -> H(',I1,') H(',I1,'), H(',I1,') -> b b')
              elseif(tlist_lep(channel) .eq. LEPh2Z_h1h1Z_4tauZ) then
		write(file_id,160) prefix,output_chan,j,j,i,i,i
 160 		FORMAT (A4,I3,':  e+ e- -> H(',I1,') Z, H(',I1,') -> H(',I1,') H(',I1,'), H(',I1,') -> tau tau')
              elseif(tlist_lep(channel) .eq. LEPh2h1_4b) then
		write(file_id,180) prefix,output_chan,j,i,j,i
 180 		FORMAT (A4,I3,':  e+ e- -> H(',I1,') H(',I1,'),  H(',I1,'), H(',I1,') -> b b')
              elseif(tlist_lep(channel) .eq. LEPh2h1_4tau) then 
		write(file_id,190) prefix,output_chan,j,i,j,i
 190 		FORMAT (A4,I3,':  e+ e- -> H(',I1,') H(',I1,'),  H(',I1,'), H(',I1,') -> tau tau')
              elseif(tlist_lep(channel) .eq. LEPh2h1_h1h1h1_6b) then 
		write(file_id,200) prefix,output_chan,j,i,j,i,i,i
 200 		FORMAT (A4,I3,':  e+ e- -> H(',I1,') H(',I1,'), H(',I1,') -> H(',I1,') H(',I1,'), H(',I1,') -> b b')
              elseif(tlist_lep(channel) .eq. LEPh2h1_h1h1h1_6tau) then
		write(file_id,210) prefix,output_chan,j,i,j,i,i,i
 210 		FORMAT (A4,I3,':  e+ e- -> H(',I1,') H(',I1,'), H(',I1,') -> H(',I1,') H(',I1,'), H(',I1,') -> tau tau')
              elseif(tlist_lep(channel) .eq. LEPh2Z_h1h1Z_2b2tau) then
		write(file_id,220) prefix,output_chan,j,j,i,i,i
 220 		FORMAT (A4,I3,':  e+ e- -> H(',I1,') Z, H(',I1,') -> H(',I1,') H(',I1,'), H(',I1,') -> b b, tau tau')
              elseif(tlist_lep(channel) .eq. LEPh2h1_2b2tau) then
		write(file_id,230) prefix,output_chan,j,i,j,i
 230 		FORMAT (A4,I3,':  e+ e- -> H(',I1,') H(',I1,'), H(',I1,') -> b b, H(',I1,') -> tau tau')
              elseif(tlist_lep(channel) .eq. LEPh2h1_2tau2b) then
		write(file_id,240) prefix,output_chan,j,i,j,i
 240 		FORMAT (A4,I3,':  e+ e- -> H(',I1,') H(',I1,'), H(',I1,') -> tau tau, H(',I1,') -> b b')
	      endif
	  endif
	endif

	if(experiment_id .eq. 2) then

	  n=TEV_reaction_number(NHIGGS,channel)

	  if(channel .le. TEV_NSCHN*NHIGGS) then
	  if(tlist_tev(channel) .eq. CDF0908_3534) write(file_id,1001) prefix,output_chan,n
	  if(tlist_tev(channel) .eq. CDF9889) write(file_id,1002) prefix,output_chan,n
          if(tlist_tev(channel) .eq. CDF9891) write(file_id,1003) prefix,output_chan,n
          if(tlist_tev(channel) .eq. D05876) write(file_id,1004) prefix,output_chan,n
          if(tlist_tev(channel) .eq. D05586) write(file_id,1005) prefix,output_chan,n
          if(tlist_tev(channel) .eq. CDF8961_D05536) write(file_id,1009) prefix,output_chan,n
          if(tlist_tev(channel) .eq. CDF9284) write(file_id,1010) prefix,output_chan,n
          if(tlist_tev(channel) .eq. D05972) write(file_id,1011) prefix,output_chan,n
          if(tlist_tev(channel) .eq. CDFwhlnubb_090814) write(file_id,1013) prefix,output_chan,n
          if(tlist_tev(channel) .eq. D05726) write(file_id,1014) prefix,output_chan,n
          if(tlist_tev(channel) .eq. D00805_3556) write(file_id,1015) prefix,output_chan,n
          if(tlist_tev(channel) .eq. D00805_2491) write(file_id,1017) prefix,output_chan,n
          if(tlist_tev(channel) .eq. CDF0906_1014) write(file_id,1018) prefix,output_chan,n
          if(tlist_tev(channel) .eq. D05873) write(file_id,1019) prefix,output_chan,n
          if(tlist_tev(channel) .eq. CDF7307v3) write(file_id,1020) prefix,output_chan,n
          if(tlist_tev(channel) .eq. D05858) write(file_id,1022) prefix,output_chan,n
          if(tlist_tev(channel) .eq. D00803_1514) write(file_id,1023) prefix,output_chan,n
          if(tlist_tev(channel) .eq. CDF9248) write(file_id,1024) prefix,output_chan,n
          if(tlist_tev(channel) .eq. CDF9290_D05645) write(file_id,1025) prefix,output_chan,n
          if(tlist_tev(channel) .eq. D05757) write(file_id,1026) prefix,output_chan,n
          if(tlist_tev(channel) .eq. CDF0809_3930) write(file_id,1027) prefix,output_chan,n
          if(tlist_tev(channel) .eq. CDF9465_D05754) write(file_id,1028) prefix,output_chan,n
          if(tlist_tev(channel) .eq. CDF9887) write(file_id,1029) prefix,output_chan,n
          if(tlist_tev(channel) .eq. D05871) write(file_id,1030) prefix,output_chan,n
          if(tlist_tev(channel) .eq. D00808_1970) write(file_id,1031) prefix,output_chan,n
          if(tlist_tev(channel) .eq. CDF0906_5613) write(file_id,1032) prefix,output_chan,n
          if(tlist_tev(channel) .eq. D00811_0024) write(file_id,1033) prefix,output_chan,n
          if(tlist_tev(channel) .eq. D05740) write(file_id,1034) prefix,output_chan,n
          if(tlist_tev(channel) .eq. D00808_1266) write(file_id,1035) prefix,output_chan,n
          if(tlist_tev(channel) .eq. CDF0802_0432) write(file_id,1036) prefix,output_chan,n
          if(tlist_tev(channel) .eq. D00901_1887) write(file_id,1037) prefix,output_chan,n
          if(tlist_tev(channel) .eq. CDF9674) write(file_id,1038) prefix,output_chan,n
          if(tlist_tev(channel) .eq. D00712_0598) write(file_id,1039) prefix,output_chan,n
          if(tlist_tev(channel) .eq. CDF9888_D05980) write(file_id,1040) prefix,output_chan,n
          if(tlist_tev(channel) .eq. CDF9713_D05889) write(file_id,1041) prefix,output_chan,n
          if(tlist_tev(channel) .eq. D05985) write(file_id,1042) prefix,output_chan,n
          if(tlist_tev(channel) .eq. CDF9897) write(file_id,1043) prefix,output_chan,n


 1001	FORMAT (A4,I3,': p p-bar -> Z H(',I1,') -> l l b b-bar, CDF with 2.7 fb^-1 (hep-ex/0908.3534)')
 1002	FORMAT (A4,I3,': p p-bar -> Z H(',I1,') -> l l b b-bar, CDF with 4.1 fb^-1 (CDF note 9889)')
 1003   FORMAT (A4,I3,': p p-bar -> V H(',I1,') -> b b-bar + miss. E_T (V=W,Z), CDF with 3.6 fb^-1 SM combined (CDF note 9891)')
 1004	FORMAT (A4,I3,': p p-bar -> Z H(',I1,') -> l l b b-bar, D0 with 4.2 fb^-1 (D0 note 5876)')
 1005	FORMAT (A4,I3,': p p-bar -> V H(',I1,') -> b b-bar + miss. E_T (V=W,Z), D0 with 2.1 fb^-1 (D0 note 5586)')
 1009	FORMAT (A4,I3,': p p-bar -> H(',I1,') + X, SM combined (CDF note 8961)')
! via HW -> l nu b b-bar, W W+ W- -> l+- nu l+- nu ..,
!                            HZ -> nu nu b b-bar, l l b b-bar,
!                            H -> W W -> l+- nu l-+ nu, SM combined
!                       (CDF note 8961
 1010	FORMAT (A4,I3,': p p-bar -> b H(',I1,') -> 3 b jets, CDF with 1.9 fb^-1 (CDF note 9284)')
 1011   FORMAT (A4,I3,': p p-bar -> W H(',I1,'), H -> b b-bar, D0 with 5.0 fb^-1 (D0 note 5972)')
 1013   FORMAT (A4,I3,': p p-bar -> W H(',I1,') -> l nu b b-bar, CDF with 4.3 fb^-1 (CDF http://www-cdf.fnal.gov/physics/new/hdg/results/whlnubb_090814/WH4.3fb.html)')
 1014   FORMAT (A4,I3,': p p-bar -> b H(',I1,') -> 3 b jets, D0 with 2.6 fb^-1 (D0 note 5726)')
 1015   FORMAT (A4,I3,': p p-bar -> b H(',I1,') -> 3 b jets, D0 (hep-ex/0805.3556)')
 1017   FORMAT (A4,I3,': p p-bar -> H(',I1,') -> tau tau, D0 with 1 fb^-1 absolute limits (D0 hep-ex/0805.2491)')
 1018   FORMAT (A4,I3,': p p-bar -> H(',I1,') -> tau tau, CDF with 1.8 fb^-1 absolute limits (CDF hep-ex/0906.1014)')
 1019   FORMAT (A4,I3,': p p-bar -> W H(',I1,') -> 3 Ws, D0 with 3.6 fb^-1 (D0 note 5873)') 
 1020   FORMAT (A4,I3,': p p-bar -> W H(',I1,') -> 3 Ws, CDF with 2.7 fb^-1 (CDF note 7307 version 3)')
 1022   FORMAT (A4,I3,': p p-bar -> H(',I1,') -> gamma gamma, D0 with 4.2 fb^-1 (D0 note 5858)')
 1023   FORMAT (A4,I3,': p p-bar -> H(',I1,') -> gamma gamma, D0 with 1.1 fb^-1 absolute limits (D0 hep-ex/0803.1514)')
 1024	FORMAT (A4,I3,': p p-bar -> H(',I1,') + X -> tau tau, CDF with 2.0 fb^-1 SM combined (CDF note 9248)')
!via HW with W -> j j, 
!                            HZ with Z -> j j,
!                            single H, VBF, SM combined
!                       (CDF note 9248,
 1025     FORMAT (A4,I3,': p p-bar -> H(',I1,') + X, CDF & D0 SM combined (CDF note 9290, D0 note 5645)')
! via HW -> l nu b b-bar, W W+ W- -> l+- nu l+- nu ..,
!                            HZ -> nu nu b b-bar, l l b b-bar,
!                            H -> tau tau via HW with W -> j j,
!                                             HZ with Z -> j j,
!                                             single H, VBF,
!			     H -> gamma gamma via ...
!                            SM combined (CDF note 9290, D0 note 5645)
 1026	FORMAT (A4,I3,': p p-bar -> H(',I1,') -> W W -> l l, D0 with 3.0 fb^-1 (D0 note 5757)')
 1027   FORMAT (A4,I3,': p p-bar -> H(',I1,') -> W W, CDF with 3.0 fb^-1 (hep-ex/0809.3930)')
 1028	FORMAT (A4,I3,': p p-bar -> H(',I1,') + X, CDF & D0 SM combined (CDF note 9465, D0 note 5754)')
 1029   FORMAT (A4,I3,': p p-bar -> H(',I1,') + X -> W W + X, CDF with 4.8 fb^-1 (CDF note 9887)')
 1030   FORMAT (A4,I3,': p p-bar -> H(',I1,') + X -> W W + X, D0 with 3.0-4.2 fb^-1 (D0 note 5871)')
 1031   FORMAT (A4,I3,': p p-bar -> W H(',I1,') -> l nu b b-bar, D0 with 1.1 fb^-1 (hep-ex/0808.1970)')
 1032   FORMAT (A4,I3,': p p-bar -> W H(',I1,') -> l nu b b-bar, CDF with 2.7 fb^-1 (hep-ex/0906.5613)')
 1033   FORMAT (A4,I3,': p p-bar -> b H(',I1,') -> b tau tau, D0 with 0.328 fb^-1 (hep-ex/0811.0024)')
 1034	FORMAT (A4,I3,': p p-bar -> H(',I1,') -> tau tau, D0 with 2.2 fb^-1 absolute limits (D0 note 5740)')
 1035	FORMAT (A4,I3,': p p-bar -> V H(',I1,') -> b b-bar + miss. E_T (V=W,Z), D0 with 0.93 fb^-1 (hep-ex/0808.1266)')
 1036	FORMAT (A4,I3,': p p-bar -> V H(',I1,') -> b b-bar + miss. E_T (V=W,Z), CDF with 1 fb^-1 (hep-ex/0802.0432)')
 1037   FORMAT (A4,I3,': p p-bar -> H(',I1,') -> gamma gamma, D0 with 2.7 fb^-1 (D0 hep-ex/0901.1887)')
 1038   FORMAT (A4,I3,': p p-bar -> H(',I1,') + X, SM combined (CDF note 9674)')
! via HW -> l nu b b-bar, W W+ W- -> l+- nu l+- nu ..,
!                            HZ -> nu nu b b-bar, l l b b-bar,
!                            H -> tau tau via HW with W -> j j,
!                                             HZ with Z -> j j,
!                                             single H, VBF,
!                            SM combined (CDF note 9674)
 1039	FORMAT (A4,I3,': p p-bar -> H(',I1,') + X, SM combined (D0 hep-ex/0712.0598)')
! via HW -> l nu b b-bar, W W+ W- -> l+- nu l+- nu .., 
!			     HZ -> nu nu b b-bar, l l b b-bar,
!			     H -> W W -> l+- nu l-+ nu, 
!			     SM combined (hep-ex/0712.0598)
 1040	FORMAT (A4,I3,': p p-bar -> H(',I1,') -> tau tau, CDF & D0 combined with 1.8 & 2.2 fb^-1 absolute limits (CDF note 9888, D0 note 5980)')
 1041   FORMAT (A4,I3,': p p-bar -> H(',I1,') + X, CDF & D0 SM combined for m_H >= 155 GeV (CDF note 9713, D0 note 5889)')
 1042   FORMAT (A4,I3,': p p-bar -> b H(',I1,') -> b tau tau, D0 with 2.7 fb^-1 (D0 note 5985)')
 1043   FORMAT (A4,I3,': p p-bar -> H(',I1,') + X, SM combined (CDF note 9897)')

	  else
!>> currently not needed. if needed: edit !
!	  i=findi(NHIGGS,n)
!	  j=findj(NHIGGS,n)
! ......
!<< currently not needed. if needed edit !
	  endif

	endif

	end
