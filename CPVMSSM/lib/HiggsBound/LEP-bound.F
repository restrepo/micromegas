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

!************************************************************
!      NOTES
!************************************************************
! LEP data files 
! * See hep-ex/0602042 'Seach for Neutral MSSM Higgs Bosons at LEP'
! * Program works out which channel is predicted to have the most power to
! exclude the parameter point, then checks the experimental result for that
! point only


!***********************************************************
	subroutine initialize_LEP_tables(NHIGGS)
!***********************************************************
	implicit none
	integer NHIGGS,i,j,m,n,imax,jmax,ntot
	real*8 Mhmin,Mhmax,sep,mh1min,mh2min,mh1max,mh2max,sep1,sep2
	real*8 mh1dummy,mh2dummy

#include "tables-lep_cb.h"

	write(*,*)'LEP-Constraints: initialization of data tables'

! tlist_lep(n) will tell the subroutines check_LEP_bound1 and check_LEP_bound2
! which LEP tables to look at
	do m=0,LEP_NSCHN-1,1
	do n=1,NHIGGS
	 tlist_lep(n+NHIGGS*m)=m+1
	enddo
	enddo
	do m=0,LEP_NDCHN-1,1
	do n=1,NHIGGS**2
	 tlist_lep(n+m*NHIGGS**2+LEP_NSCHN*NHIGGS)=m+1+LEP_NSCHN
	enddo				
	enddo

! initialize table 14b
	write(*,*)'reading table 14b (EPJC 46(2006)547)'

! information about how LEP data is stored in file:	 
	Mhmin14b=12.0D0
	Mhmax14b=120.0D0
	sep14b=0.5D0

	Mhmin=Mhmin14b
	Mhmax=Mhmax14b
        sep=sep14b
	ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 217. 
! If not change appropriate declarations in 'tables-lep_cb.h'

	open(10,file=INSTDIR//"/LEPtables/lep210_hbb.txt")

! first 5 lines are not needed	 
        do i=1,5
	  read(10,*)  	 
        end do	 	
		 
        do i=1,ntot
	  read(10,*)Mh1dummy,S95_14b_obs(i),S95_14b_pred(i)
        enddo

        close(10)	 	
	

! initialize table 14c
	write(*,*)'reading table 14c (EPJC 46(2006)547)'

! information about how LEP data is stored in file:	 
	Mhmin14c=4.0D0
	Mhmax14c=120.0D0
	sep14c=0.5D0

	Mhmin=Mhmin14c
	Mhmax=Mhmax14c
        sep=sep14c
	ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 233. 

	open(10,file=INSTDIR//"/LEPtables/lep210_htt_interpol.txt")

! first 5 lines are not needed	 
        do i=1,5
	  read(10,*)  	 
        end do	 	
		 
        do i=1,ntot
          read(10,*)Mh1dummy,S95_14c_obs(i),S95_14c_pred(i)
        enddo

        close(10)	 	
	

	write(*,*)'reading tables 15, 16 and 22 (EPJC 46(2006)547)'

! initialize tables 15,16,18,19,20,21,22,23,24 :

! fill S95 from files
! row 0 and column 0 in LEP file contain higgs masses 
! and (0,0) ie top left set to -100	

! tables 15, 16 and 22 :

	Mh1min15=1.0D0
        Mh1max15=60.0D0
        Mh2min15=2.0D0
        Mh2max15=120.0D0            

	sep1_15=1.0D0
	sep2_15=1.0D0
	
	Mh1min=Mh1min15
	Mh1max=Mh1max15
	Mh2min=Mh2min15
	Mh2max=Mh2max15
	sep1=sep1_15
	sep2=sep2_15

	imax=int((Mh2max-Mh2min))/sep2+1
	jmax=int((Mh1max-Mh1min))/sep1+1	
! (imax,jmax) is expected to be (119,60). 
	
	open(10,file=INSTDIR//"/LEPtables/table15full_pred.txt")
	open(11,file=INSTDIR//"/LEPtables/table15full_obs.txt")
        open(20,file=INSTDIR//"/LEPtables/table16full_pred.txt")
        open(21,file=INSTDIR//"/LEPtables/table16full_obs.txt")
        open(30,file=INSTDIR//"/LEPtables/table22full_pred.txt")
        open(31,file=INSTDIR//"/LEPtables/table22full_obs.txt")

	do i=0,imax  
         read(10,*) ( S95_15_pred(i,j), j=0,jmax )	   	 
         read(11,*) ( S95_15_obs(i,j), j=0,jmax )	   	 
         read(20,*) ( S95_16_pred(i,j), j=0,jmax )	   	 
         read(21,*) ( S95_16_obs(i,j), j=0,jmax )	   	 
         read(30,*) ( S95_22_pred(i,j), j=0,jmax )	   	 
         read(31,*) ( S95_22_obs(i,j), j=0,jmax )	   	 

        enddo
		
! set entries .ge. 100 to (-4): they will not be relevant
! Remark: If S95-value > 100 we don't apply the results anymore.
! In regions where LEP runs out of steam, S95 is sharply rising
! and there are usually a few more reasonable points with S95 > 100.
! Beyond that, S95 is abitrarily set to 10000.
! For S95 above about 100, interpolation is not very reliable anymore.
! Therefore, we stop using the tables there.
	do i=1,imax
	 do j=1,jmax
	  if(S95_15_pred(i,j) .ge. 100d0) S95_15_pred(i,j)= -4.0D0	
	  if(S95_15_obs(i,j)  .ge. 100d0) S95_15_obs(i,j) = -4.0D0	
	  if(S95_16_pred(i,j) .ge. 100d0) S95_16_pred(i,j)= -4.0D0	
	  if(S95_16_obs(i,j)  .ge. 100d0) S95_16_obs(i,j) = -4.0D0	
	  if(S95_22_pred(i,j) .ge. 100d0) S95_22_pred(i,j)= -4.0D0	
	  if(S95_22_obs(i,j)  .ge. 100d0) S95_22_obs(i,j) = -4.0D0	
	 enddo   	 
        enddo 	
	
	close(10)
        close(11)
        close(20)
        close(21)
        close(30)
        close(31)
	
	write(*,*)'reading tables 18, 19, 23, 24 (EPJC 46(2006)547)'

! tables 18, 19, 23, 24 :

	Mh1min18=1.0D0
        Mh1max18=180.0D0
        Mh2min18=1.0D0
        Mh2max18=180.0D0

	sep1_18=1.0D0
	sep2_18=1.0D0
	
	Mh1min=Mh1min18
	Mh1max=Mh1max18
	Mh2min=Mh2min18
	Mh2max=Mh2max18
	sep1=sep1_18
	sep2=sep2_18

	imax=int((Mh2max-Mh2min))/sep2+1
	jmax=int((Mh1max-Mh1min))/sep1+1	
! (imax,jmax) is expected to be (180,180). 

	open(10,file=INSTDIR//"/LEPtables/table18full_pred.txt")
	open(11,file=INSTDIR//"/LEPtables/table18full_obs.txt")
        open(20,file=INSTDIR//"/LEPtables/table19full_pred.txt")
        open(21,file=INSTDIR//"/LEPtables/table19full_obs.txt")
        open(30,file=INSTDIR//"/LEPtables/table23full_pred.txt")
        open(31,file=INSTDIR//"/LEPtables/table23full_obs.txt")
        open(40,file=INSTDIR//"/LEPtables/table24full_pred.txt")
        open(41,file=INSTDIR//"/LEPtables/table24full_obs.txt")

	do i=0,imax  
         read(10,*) ( S95_18_pred(i,j), j=0,jmax )	   	 
         read(11,*) ( S95_18_obs(i,j), j=0,jmax )	   	 
         read(20,*) ( S95_19_pred(i,j), j=0,jmax )	   	 
         read(21,*) ( S95_19_obs(i,j), j=0,jmax )	   	 
         read(30,*) ( S95_23_pred(i,j), j=0,jmax )	   	 
         read(31,*) ( S95_23_obs(i,j), j=0,jmax )	   	 
         read(40,*) ( S95_24_pred(i,j), j=0,jmax )	   	 
         read(41,*) ( S95_24_obs(i,j), j=0,jmax )	   	 
        enddo
		
! set entries .ge. 100 to (-4): they will not be relevant
! Remark: If S95-value > 100 we don't apply the results anymore.
! In regions where LEP runs out of steam, S95 is sharply rising
! and there are usually a few more reasonable points with S95 > 100.
! Beyond that, S95 is abitrarily set to 10000.
! For S95 above about 100, interpolation is not very reliable anymore.
! Therefore, we stop using the tables there.
	do i=1,imax
	 do j=1,jmax
	  if(S95_18_pred(i,j) .ge. 100d0) S95_18_pred(i,j)= -4.0D0	
	  if(S95_18_obs(i,j)  .ge. 100d0) S95_18_obs(i,j) = -4.0D0	
	  if(S95_19_pred(i,j) .ge. 100d0) S95_19_pred(i,j)= -4.0D0	
	  if(S95_19_obs(i,j)  .ge. 100d0) S95_19_obs(i,j) = -4.0D0	
	  if(S95_23_pred(i,j) .ge. 100d0) S95_23_pred(i,j)= -4.0D0	
	  if(S95_23_obs(i,j)  .ge. 100d0) S95_23_obs(i,j) = -4.0D0	
	  if(S95_24_pred(i,j) .ge. 100d0) S95_24_pred(i,j)= -4.0D0	
	  if(S95_24_obs(i,j)  .ge. 100d0) S95_24_obs(i,j) = -4.0D0	
	 enddo   	 
        enddo 	
	
	close(10)
        close(11)
        close(20)
        close(21)
        close(30)
        close(31)
        close(40)
        close(41)

	write(*,*)'reading tables 20, 21 (EPJC 46(2006)547)'

! tables 20, 21

	Mh1min20=1.0D0
        Mh1max20=90.0D0
        Mh2min20=2.0D0
        Mh2max20=180.0D0

	sep1_20=1.0D0
	sep2_20=1.0D0
	
	Mh1min=Mh1min20
	Mh1max=Mh1max20
	Mh2min=Mh2min20
	Mh2max=Mh2max20
	sep1=sep1_20
	sep2=sep2_20

	imax=int((Mh2max-Mh2min))/sep2+1
	jmax=int((Mh1max-Mh1min))/sep1+1	
! (imax,jmax) is expected to be (179,90). 

	open(10,file=INSTDIR//"/LEPtables/table20full_pred.txt")
	open(11,file=INSTDIR//"/LEPtables/table20full_obs.txt")
        open(20,file=INSTDIR//"/LEPtables/table21full_pred.txt")
        open(21,file=INSTDIR//"/LEPtables/table21full_obs.txt")

	do i=0,imax  
         read(10,*) ( S95_20_pred(i,j), j=0,jmax )	   	 
         read(11,*) ( S95_20_obs(i,j), j=0,jmax )	   	 
         read(20,*) ( S95_21_pred(i,j), j=0,jmax )	   	 
         read(21,*) ( S95_21_obs(i,j), j=0,jmax )	   	 
        enddo
		
! set entries .ge. 100 to (-4): they will not be relevant
! Remark: If S95-value > 100 we don't apply the results anymore.
! In regions where LEP runs out of steam, S95 is sharply rising
! and there are usually a few more reasonable points with S95 > 100.
! Beyond that, S95 is abitrarily set to 10000.
! For S95 above about 100, interpolation is not very reliable anymore.
! Therefore, we stop using the tables there.
	do i=1,imax
	 do j=1,jmax
	  if(S95_20_pred(i,j) .ge. 100d0) S95_20_pred(i,j)= -4.0D0	
	  if(S95_20_obs(i,j)  .ge. 100d0) S95_20_obs(i,j) = -4.0D0	
	  if(S95_21_pred(i,j) .ge. 100d0) S95_21_pred(i,j)= -4.0D0	
	  if(S95_21_obs(i,j)  .ge. 100d0) S95_21_obs(i,j) = -4.0D0	
	 enddo   	 
        enddo 	
	
	close(10)
        close(11)
        close(20)
        close(21)


!new additions
! new 1) channel e+ e- -> hZ, h -> anything
        write(*,*)'reading table for analysis: e+ e- -> hZ, h -> anything (EPJC 27(2003)311)'

! information about how LEP data is stored in file:	 
	Mhmin_hZ_anyZ=1.0D0
	Mhmax_hZ_anyZ=100.0D0
	sep_hZ_anyZ=1D0

	Mhmin=Mhmin_hZ_anyZ
	Mhmax=Mhmax_hZ_anyZ
        sep=sep_hZ_anyZ
	ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 100. 

	open(10,file=INSTDIR//"/LEPtables/lep_s95_ee_hZ_anythingZ.txt")

! first 5 lines are not needed	 
        do i=1,5
	  read(10,*)  	 
        end do	 	
		 
        do i=1,ntot
          read(10,*)Mh1dummy,S95_hZ_anyZ_pred(i),S95_hZ_anyZ_obs(i)
        enddo

        close(10)	 	
	
! new 2) channel e+ e- -> hZ, h -> gamma gamma
        write(*,*)'reading table for analysis: e+ e- -> hZ, h -> gamma gamma (LHWG note 2002-02)'

! information about how LEP data is stored in file:	 
	Mhmin_hZ_gamgamZ=20D0
	Mhmax_hZ_gamgamZ=116.0D0
	sep_hZ_gamgamZ=2D0

	Mhmin=Mhmin_hZ_gamgamZ
	Mhmax=Mhmax_hZ_gamgamZ
        sep=sep_hZ_gamgamZ
	ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 97. 

	open(10,file=INSTDIR//"/LEPtables/lep_s95_ee_hZ_gamgamZ.txt")

! first 5 lines are not needed	 
        do i=1,5
	  read(10,*)  	 
        end do	 	
		 
        do i=1,ntot
          read(10,*)Mh1dummy,S95_hZ_gamgamZ_obs(i),S95_hZ_gamgamZ_pred(i)
        enddo

        close(10)	 	
	

	write(*,*)'finished reading LEP tables'
	write(*,*)

	end




!***********************************************************
	subroutine check_LEP_bound(NHIGGS,whichexpt,
     &		Mh,R_hjZ_ee,R_hjhi_ee_compact, 
     &          BRhjbbbar,BRhjtautau,BRhjhihi_compact,BRhjgamgam,
     &		allowed,result_flag,
     &		mlratio,obsratio,most_sens_predratio,num_H_combined,
     &		final_h_comb,model_R_sig_BR)
!***********************************************************
! Input:
!	NHIGGS : # of neutral Higgs bosons 
!	Mh(NHMAX) : the neutral Higgs masses
! 	R_hjZ_ee(NHMAX) : ratio sigma(e+e-  -> H_j Z)/sigma(e+e-  -> H Z)_SM
!	R_hjhi_ee_compact(..) -> R_hjhi_ee(NHMAX,NHMAX) : 
!				ratio sigma(e+e-  -> H_j H_i)/sigma_ref
!	BRhjbbbar(NHMAX): branching ratio for H_j -> b b-bar
!	BRhjtautau(NHMAX): branching ratio for H_j -> tau+ tau-
!	BRhjhihi_compact(..) -> BRhjhihi(NHMAX,NHMAX) : 
!				branching ratio for H_j -> H_i H_i
!				only off-diagonal elements relevant
!	BRhjgamgam(NHMAX): branching ratio for H_j -> gamma gamma
! Output:
!	allowed : logical variable, TRUE if scenario is allowed at
!		95% confidence level (This is determined, using the procedure
!	 	described by Philipp Bechtle.), FALSE otherwise.
!	result_flag
!		1 if point is allowed, 
!		0 if point is not allowed
!	       -1 otherwise (theoretically inaccessible, etc.)
!	mlratio : gives information about which table was used for the 
!		experimental comparison.
!	obsratio, most_sens_predratio, num_H_combined, final_h_comb, model_R_sig_BR 
!--------------------------------
	implicit none 

	logical allowed,testKR,useTEV,useLEP,useSingH
	character*5 whichexpt

	integer NHIGGS,my_maxloc,output_chan
	integer result_flag,num_H_combined
	integer i,j,k,n,m,mlratio,nhelp
	integer mlratiotemp(1)	
	double precision R_sig_BR(LEP_NSCHN*NHMAX+LEP_NDCHN*NHMAX**2)
        double precision R_sig_BR_sing(LEP_NSCHN*NHMAX+LEP_NDCHN*NHMAX**2)
        double precision predratio(LEP_NSCHN*NHMAX+LEP_NDCHN*NHMAX**2)
	double precision obsratio,most_sens_predratio
	double precision BRhjbbbar(NHMAX),BRhjtautau(NHMAX),BRhjgamgam(NHMAX)
	double precision BRhjhihi(NHMAX,NHMAX)		
        double precision BRhjhihi_compact(NHMAX**2-NHMAX)

	double precision Mh(NHMAX)
	double precision R_hjZ_ee(NHMAX),R_hjhi_ee(NHMAX,NHMAX)
	double precision R_hjhi_ee_compact((NHMAX*(NHMAX+1))/2)

	logical higgs_combination(LEP_NSCHN*NHMAX+LEP_NDCHN*NHMAX**2,NHMAX),final_h_comb(NHMAX)
	integer sim_mass_sum(LEP_NSCHN*NHMAX+LEP_NDCHN*NHMAX**2)
        double precision Mh_avg(LEP_NSCHN*NHMAX),delta_Mh
	double precision log2dble,model_R_sig_BR

	integer tabsize,tabsizei,tabsizej
	double precision Mhmin,Mhmax,sep
	double precision Mh1min,Mh1max,sep1
	double precision Mh2min,Mh2max,sep2

	integer findi,findj,LEP_reaction_number,number_LEP_channels
#include "tables-lep_cb.h"

	delta_Mh = delta_Mh_LEP

        call determine_topology_flags(whichexpt,
     &          useTEV,useLEP,useSingH)

! translate BRhjhihi_compact to BRhjhihi
	k=0
        do j=1,NHIGGS
          do i=1,NHIGGS
            if(j .ne. i) then
              k=k+1
              BRhjhihi(j,i)=BRhjhihi_compact(k)
            else
              BRhjhihi(j,i)=0d0
            endif
          enddo
        enddo

! translate R_hjhi_ee_compact to  R_hjhi_ee
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

! tlist(n) will tell the subroutines check_LEP_bound1 and check_LEP_bound2
! which LEP tables to look at
	do m=0,LEP_NSCHN-1,1
	do n=1,NHIGGS
	 R_sig_BR_sing(n+NHIGGS*m)=0d0
         R_sig_BR(n+NHIGGS*m)=0d0
	 sim_mass_sum(n+NHIGGS*m)=0
	 do i=1,NHIGGS
          higgs_combination(n+NHIGGS*m,i)=.FALSE.
	 enddo
	enddo
	enddo
	do m=0,LEP_NDCHN-1,1
	do n=1,NHIGGS**2
	 R_sig_BR_sing(n+m*NHIGGS**2+LEP_NSCHN*NHIGGS)=0d0
	 R_sig_BR(n+m*NHIGGS**2+LEP_NSCHN*NHIGGS)=0d0
	 do i=1,NHIGGS
	  higgs_combination(n+m*NHIGGS**2+LEP_NSCHN*NHIGGS,i)=.FALSE.
         enddo
	enddo				
	enddo

! a0) calc single sigma x BR values (no combinations of near mass Higgs bosons)
!    single Higgs channels
	do n=1,NHIGGS
	  R_sig_BR_sing(n+NHIGGS*(LEPhZ_bbZ -1))=R_hjZ_ee(n)*BRhjbbbar(n)
          R_sig_BR_sing(n+NHIGGS*(LEPhZ_tautauZ -1))=R_hjZ_ee(n)*BRhjtautau(n)	
          R_sig_BR_sing(n+NHIGGS*(LEPhZ_anyZ -1))=R_hjZ_ee(n)
	  R_sig_BR_sing(n+NHIGGS*(LEPhZ_gamgamZ -1))=R_hjZ_ee(n)*BRhjgamgam(n)
	enddo
!    double Higgs channels: not needed because no adding of signals 
!    is (so far) done for double Higgs channels

! a) calc sum of sigma x BR contributions to R_sig_BR
!   single Higgs channels
!   R_sig_BR(n) will be compared to the S95 values in the LEP tables
	  do n=1,NHIGGS
	    do i=1,NHIGGS
	      if( (Mh(i) .ge. Mh(n)) .and. ( (Mh(i)-Mh(n)) .le. delta_Mh) ) then
		if(R_sig_BR_sing(i+NHIGGS*(LEPhZ_bbZ -1))*R_sig_BR_sing(n+NHIGGS*(LEPhZ_bbZ -1)) .ne. 0d0) then
	          R_sig_BR(n+NHIGGS*(LEPhZ_bbZ -1))=R_sig_BR(n+NHIGGS*(LEPhZ_bbZ -1))
     &			+R_sig_BR_sing(i+NHIGGS*(LEPhZ_bbZ -1))
		endif
		if(R_sig_BR_sing(i+NHIGGS*(LEPhZ_tautauZ -1))*R_sig_BR_sing(n+NHIGGS*(LEPhZ_tautauZ -1)) .ne. 0d0) then
	          R_sig_BR(n+NHIGGS*(LEPhZ_tautauZ -1))=R_sig_BR(n+NHIGGS*(LEPhZ_tautauZ -1))
     &			+R_sig_BR_sing(i+NHIGGS*(LEPhZ_tautauZ -1))
		endif
		if(R_sig_BR_sing(i+NHIGGS*(LEPhZ_anyZ -1))*R_sig_BR_sing(n+NHIGGS*(LEPhZ_anyZ -1)) .ne. 0d0) then
	          R_sig_BR(n+NHIGGS*(LEPhZ_anyZ -1))=R_sig_BR(n+NHIGGS*(LEPhZ_anyZ -1))
     &			+R_sig_BR_sing(i+NHIGGS*(LEPhZ_anyZ -1))
		endif
		if(R_sig_BR_sing(i+NHIGGS*(LEPhZ_gamgamZ -1))*R_sig_BR_sing(n+NHIGGS*(LEPhZ_gamgamZ -1)) .ne. 0d0) then
	          R_sig_BR(n+NHIGGS*(LEPhZ_gamgamZ -1))=R_sig_BR(n+NHIGGS*(LEPhZ_gamgamZ -1))
     &			+R_sig_BR_sing(i+NHIGGS*(LEPhZ_gamgamZ -1))
		endif
	      endif
	    enddo
	  enddo
!   double Higgs channels
	  if(.not. useSingH) then	  
	  do n=1,NHIGGS**2
	   i=findi(NHIGGS,n)
	   j=findj(NHIGGS,n)
           
	   R_sig_BR(n+(LEPh2Z_h1h1Z_4bZ -1-LEP_NSCHN)*NHIGGS**2+LEP_NSCHN*NHIGGS)=R_hjZ_ee(j)*BRhjhihi(j,i)*BRhjbbbar(i)**2.0D0	
	   R_sig_BR(n+(LEPh2Z_h1h1Z_4tauZ -1-LEP_NSCHN)*NHIGGS**2+LEP_NSCHN*NHIGGS)=R_hjZ_ee(j)*BRhjhihi(j,i)*BRhjtautau(i)**2.0D0
	   	 
	   R_sig_BR(n+(LEPh2h1_4b -1-LEP_NSCHN)*NHIGGS**2+LEP_NSCHN*NHIGGS)=R_hjhi_ee(j,i)*BRhjbbbar(j)*BRhjbbbar(i)	 	     
	   R_sig_BR(n+(LEPh2h1_4tau -1-LEP_NSCHN)*NHIGGS**2+LEP_NSCHN*NHIGGS)=R_hjhi_ee(j,i)*BRhjtautau(j)*BRhjtautau(i)	
	    
	   R_sig_BR(n+(LEPh2h1_h1h1h1_6b -1-LEP_NSCHN)*NHIGGS**2+LEP_NSCHN*NHIGGS)=R_hjhi_ee(j,i)*BRhjhihi(j,i)*BRhjbbbar(i)**3.0D0	 	   
	   R_sig_BR(n+(LEPh2h1_h1h1h1_6tau -1-LEP_NSCHN)*NHIGGS**2+LEP_NSCHN*NHIGGS)=R_hjhi_ee(j,i)*BRhjhihi(j,i)*BRhjtautau(i)**3.0D0	
	    	   
	   R_sig_BR(n+(LEPh2Z_h1h1Z_2b2tau -1-LEP_NSCHN)*NHIGGS**2+LEP_NSCHN*NHIGGS)=R_hjZ_ee(j)*BRhjhihi(j,i)*BRhjbbbar(i)*BRhjtautau(i)
	   	
	   R_sig_BR(n+(LEPh2h1_2b2tau -1-LEP_NSCHN)*NHIGGS**2+LEP_NSCHN*NHIGGS)=R_hjhi_ee(j,i)*BRhjbbbar(j)*BRhjtautau(i)	   	 
	   R_sig_BR(n+(LEPh2h1_2tau2b -1-LEP_NSCHN)*NHIGGS**2+LEP_NSCHN*NHIGGS)=R_hjhi_ee(j,i)*BRhjtautau(j)*BRhjbbbar(i)		   
	  enddo	               	  
	  endif
	
	  allowed = .TRUE.


!/!\: currently Mh_avg is not used and also not needed(!) 
!     for delta_Mh_LEP of the order of 1-2 GeV or below
! b) Determine the average Higgs mass for the Higgs combinations
!     (so far only done for the single Higgs channels)

	do m=0,LEP_NSCHN-1,1
	  do n=1,NHIGGS
	    sim_mass_sum(n+NHIGGS*m)=1
            higgs_combination(n+NHIGGS*m,n)=.TRUE.
	    Mh_avg(n+NHIGGS*m)=Mh(n)
	    do i=1,NHIGGS
              if(i .ne. n) then
	        if( (Mh(i) .ge. Mh(n)) .and. ( (Mh(i)-Mh(n)) .le. delta_Mh )
     &	 	  .and. ( R_sig_BR_sing(i+NHIGGS*m)*R_sig_BR_sing(n+NHIGGS*m) .ne. 0d0)  ) then
!not yet needed     &		    .and.  table_applicable(i+NHIGGS*m) ) then
	          Mh_avg(n+NHIGGS*m)=Mh_avg(n+NHIGGS*m)+Mh(i)
                  sim_mass_sum(n+NHIGGS*m)=sim_mass_sum(n+NHIGGS*m)+1
	          higgs_combination(n+NHIGGS*m,i)=.TRUE.
	        endif
	      endif
	    enddo
            Mh_avg(n+NHIGGS*m)=Mh_avg(n+NHIGGS*m)/dble(sim_mass_sum(n+NHIGGS*m))
          enddo
	enddo
! no Higgs mass combination done for the double Higgs channels 
	do m=0,LEP_NDCHN-1,1
	  do n=1,NHIGGS**2
	    sim_mass_sum(n+m*NHIGGS**2+LEP_NSCHN*NHIGGS)=1
	  enddo				
	enddo

	    
	   do n=1,NHIGGS
! comparing input to the predicted LEP results	   	  		      
	    do k=0,LEP_NSCHN-1,1 	   
	    call get_1d_LEPtable_data(tlist_lep(n+NHIGGS*k),tabsize,Mhmin,Mhmax,sep)

	     call check_LEP_bound1(Mh_avg(n+NHIGGS*k), 
     &	      R_sig_BR(n+NHIGGS*k),predratio(n+NHIGGS*k),tlist_lep(n+NHIGGS*k),
     &        tabsize,Mhmin,Mhmax,sep,.True.)
	    enddo
	   enddo



           if(.not. useSingH) then
           do n=1,NHIGGS**2	 	   
            i=findi(NHIGGS,n)
            j=findj(NHIGGS,n)

	    do k=0,LEP_NDCHN-1,1

	     call get_2d_LEPtable_data(tlist_lep(n+k*NHIGGS**2+LEP_NSCHN*NHIGGS),
     &      	tabsizei,tabsizej,Mh1min,Mh2min,Mh1max,Mh2max,sep1,sep2)

	     call check_LEP_bound2(Mh(i),Mh(j), 
     &	      R_sig_BR(n+k*NHIGGS**2+LEP_NSCHN*NHIGGS),predratio(n+k*NHIGGS**2+LEP_NSCHN*NHIGGS),
     &	      tlist_lep(n+k*NHIGGS**2+LEP_NSCHN*NHIGGS),
     &        tabsizei,tabsizej,Mh1min,Mh2min,Mh1max,Mh2max,sep1,sep2,
     &	      .True.)	    

	    enddo 	    
           enddo	  		  
	   endif

! find out which channel is best for comparing with experimental result	  

        if(useSingH) then
	  number_LEP_channels = LEP_NSCHN*NHIGGS 
	else
	  number_LEP_channels =	LEP_NSCHN*NHIGGS+LEP_NDCHN*NHIGGS**2
	endif

#ifdef DEBUG
	do nhelp=1,number_LEP_channels,1
	    call determine_output_channel_number(NHIGGS,whichexpt,
     &          1,nhelp,output_chan)
	    call write_output_channel(6,NHIGGS,whichexpt,output_chan)
	    write(*,*)'   R_sig_BR=sig*BR/sig*BR_ref =',R_sig_BR(nhelp)
            write(*,*)'   R_sig_BR/S95_pred.         =',predratio(nhelp) 
	enddo
#endif

	  mlratio=my_maxloc(predratio,number_LEP_channels)

! Determine the reaction number n connected to mlratio.
! For mlratio <= LEP_NSCHN*NHIGGS there are three reactions per table, 
! i.e. n is identical with the Higgs species.
! For mlratio >  LEP_NSCHN*NHIGGS there are NHIGGS**2-NHIGGS reactions per table,
! i.e. n lies between 1 and NHIGGS**2
	
	n=LEP_reaction_number(NHIGGS,mlratio)

! compare with R_sig_BR(mlratio) with experimental S95 value
	  if((mlratio .le. LEP_NSCHN*NHIGGS))then
	    call get_1d_LEPtable_data(tlist_lep(mlratio),tabsize,Mhmin,Mhmax,sep)
	    call check_LEP_bound1(Mh_avg(mlratio), 
     &	     R_sig_BR(mlratio),obsratio,tlist_lep(mlratio),
     &       tabsize,Mhmin,Mhmax,sep,.False.)        
	  else
	    i=findi(NHIGGS,n)	      
	    j=findj(NHIGGS,n)
	    call get_2d_LEPtable_data(tlist_lep(mlratio),
     &      	tabsizei,tabsizej,Mh1min,Mh2min,Mh1max,Mh2max,sep1,sep2)
	    call check_LEP_bound2(Mh(i),Mh(j), 
     &	     R_sig_BR(mlratio),obsratio,tlist_lep(mlratio),
     &        tabsizei,tabsizej,Mh1min,Mh2min,Mh1max,Mh2max,sep1,sep2,
     &	      .False.)
	  endif

! obsratio = -1 is considered "allowed". This is where the tables give no exclusion, e.g. because the 
!	Higgs mass is out of range.
          if(obsratio .lt. 1d0) then
	   result_flag=1
	   allowed=.True.
	  else
	   result_flag=0
	   allowed=.False.
	  endif

	model_R_sig_BR=R_sig_BR(mlratio)
	most_sens_predratio=predratio(mlratio)
	num_H_combined=sim_mass_sum(mlratio)
	do i=1,NHIGGS
		final_h_comb(i)=higgs_combination(mlratio,i)
	enddo
#ifdef VERBOSE
	write(*,*)'LEP channel with highest sensitivity:'
	    call determine_output_channel_number(NHIGGS,whichexpt,
     &          1,mlratio,output_chan)
	    call write_output_channel(6,NHIGGS,whichexpt,output_chan)

	write(*,*)'   R_sig_BR/S95_predicted =',most_sens_predratio
	write(*,*)'   R_sig_BR/S95_observed  =',obsratio
#endif

	end



!******************************************************************
	integer function findj(NHIGGS,n)
!******************************************************************
	implicit none
	integer n,NHIGGS

	findj=int((n-1)/NHIGGS)+1
	
	end


!******************************************************************
	integer function findi(NHIGGS,n)
!******************************************************************
	implicit none
	integer n,NHIGGS

	findi=n-NHIGGS*int((n-1)/NHIGGS)		

        end


!******************************************************************
	integer function LEP_reaction_number(NHIGGS,channel)
!******************************************************************
	implicit none
	integer channel,NHIGGS

	if(channel .le. LEP_NSCHN*NHIGGS) then
	  LEP_reaction_number=channel-int((channel-1)/NHIGGS)*NHIGGS
	else
	  LEP_reaction_number=channel-LEP_NSCHN*NHIGGS-int((channel-LEP_NSCHN*NHIGGS-1)/NHIGGS**2)*NHIGGS**2
	endif

	end



!******************************************************************
	subroutine get_1d_LEPtable_data(selec,tabsize,Mhmin,Mhmax,sep)
!******************************************************************
	implicit none
	integer tabsize,selec
        double precision sep
	double precision Mhmin,Mhmax

#include "tables-lep_cb.h"

	if(selec .eq. LEPhZ_bbZ) then
	  Mhmin=Mhmin14b
	  Mhmax=Mhmax14b
          sep=sep14b
	elseif(selec .eq. LEPhZ_tautauZ) then
       	  Mhmin=Mhmin14c
          Mhmax=Mhmax14c
          sep=sep14c
	elseif(selec .eq. LEPhZ_anyZ) then
       	  Mhmin=Mhmin_hZ_anyZ
          Mhmax=Mhmax_hZ_anyZ
          sep=sep_hZ_anyZ
	elseif(selec .eq. LEPhZ_gamgamZ) then
       	  Mhmin=Mhmin_hZ_gamgamZ
          Mhmax=Mhmax_hZ_gamgamZ
          sep=sep_hZ_gamgamZ
        endif

	tabsize=nint((Mhmax-Mhmin)/sep)+1

	end



!******************************************************************
	subroutine check_LEP_bound1(Mh1, 
     &	  factBR,ratio,selec,tabsize,Mhmin,Mhmax,sep,pred)
!******************************************************************
	implicit none
	logical pred
	integer tabsize
	integer i,ilow,ihigh,selec
	double precision Mh1,Mh1_test,Mh1_below,factBR
        double precision sep
	double precision Mhmin,Mhmax,ratio
	double precision Mh1dummy,Mh2dummy,S95obs,S95pred	
        double precision S95_above,S95_below,S95_interpol
        double precision S95(tabsize)

#include "tables-lep_cb.h"

        ratio=0.0D0
	
! check if mass is within mass range of table:
	if(Mh1 .gt. Mhmax) then
	  ratio=0D0
	elseif(Mh1 .lt. Mhmin) then	 
	  ratio=0D0
	else	

	if(selec .eq. LEPhZ_bbZ) then
	  if(pred) then
		call eval_1d_arr(1,tabsize,S95,S95_14b_pred)
	  else
		call eval_1d_arr(1,tabsize,S95,S95_14b_obs)
	  endif 
	elseif(selec .eq. LEPhZ_tautauZ) then
          if(pred) then 
		call eval_1d_arr(1,tabsize,S95,S95_14c_pred)
          else 
		call eval_1d_arr(1,tabsize,S95,S95_14c_obs)
          endif
	elseif(selec .eq. LEPhZ_anyZ) then
          if(pred) then 
                call eval_1d_arr(1,tabsize,S95,S95_hZ_anyZ_pred)
          else 
                call eval_1d_arr(1,tabsize,S95,S95_hZ_anyZ_obs)
          endif
	elseif(selec .eq. LEPhZ_gamgamZ) then
          if(pred) then 
                call eval_1d_arr(1,tabsize,S95,S95_hZ_gamgamZ_pred)
          else 
                call eval_1d_arr(1,tabsize,S95,S95_hZ_gamgamZ_obs)
          endif


	endif


! linear interpolate between points in table
        ilow=int((Mh1-Mhmin)/sep)+1

! check if interpolation is needed 
!	The expression in the following if-statement should be mathematically zero.
!	However, numerically, it is sometimes only nearly zero.
	if( ((Mh1-Mhmin)/sep+1d0-dble(ilow)) .lt. 1d-7) then
		S95_interpol=S95(ilow)
	else
	 i=ilow	 	 
	 
	 Mh1_below=dble(i-1)*sep+Mhmin
	 S95_below=S95(i)
	 S95_above=S95(i+1)
     	 S95_interpol=S95_below +(S95_above-S95_below)/sep*(Mh1-Mh1_below) 
	endif


        	 	 
! factBR is the ratio (sigma*BRs)_model/(sigma*BRs)_ref.
! factBR for the channel with the highest factBR/S95_interpol_predicted
! has to be smaller then S95_interpol_observed in order to yield an
! allowed scenario. 
! The latter is equivalent to ratio=factBR/S95_interpol_observed < 1.
	 ratio=factBR/S95_interpol		
	  
	endif 
	
	end


!******************************************************************
	subroutine get_2d_LEPtable_data(selec,
     &      tabsizei,tabsizej,Mh1min,Mh2min,Mh1max,Mh2max,sep1,sep2)
!******************************************************************
	implicit none
	integer tabsizei,tabsizej
	integer selec
	double precision Mh1min,Mh1max,Mh2min,Mh2max
        double precision sep1,sep2

#include "tables-lep_cb.h"


	if((selec .eq. LEPh2Z_h1h1Z_4bZ) .or. (selec .eq. LEPh2Z_h1h1Z_4tauZ) 
     &		.or. (selec .eq. LEPh2Z_h1h1Z_2b2tau)) then
	 Mh1min=Mh1min15
	 Mh1max=Mh1max15
         Mh2min=Mh2min15
         Mh2max=Mh2max15
	 sep1=sep1_15
	 sep2=sep2_15
	elseif((selec .eq. LEPh2h1_4b) .or. (selec .eq. LEPh2h1_4tau) .or. 
     &		(selec .eq. LEPh2h1_2b2tau) .or. (selec .eq. LEPh2h1_2tau2b)) then
         Mh1min=Mh1min18
         Mh1max=Mh1max18
         Mh2min=Mh2min18
         Mh2max=Mh2max18
         sep1=sep1_18
         sep2=sep2_18
	elseif((selec .eq. LEPh2h1_h1h1h1_6b) .or. (selec .eq. LEPh2h1_h1h1h1_6tau)) then
         Mh1min=Mh1min20
         Mh1max=Mh1max20
         Mh2min=Mh2min20
         Mh2max=Mh2max20
         sep1=sep1_20
         sep2=sep2_20
	else
         write(*,*)'ERROR B' 
	 write(*,*)'selec=',selec
         stop
        endif
	
	tabsizei=int((Mh2max-Mh2min))/sep2+1
	tabsizej=int((Mh1max-Mh1min))/sep1+1	

	end


!******************************************************************
	subroutine check_LEP_bound2(Mh1in,Mh2in, 
     &	  factBR,ratio,selec,
     &    imax,jmax,Mh1min,Mh2min,Mh1max,Mh2max,sep1,sep2,pred)
!******************************************************************	  
	implicit none
	logical pred

	integer imax,jmax,selec
	integer i,j
	integer ilow,jlow,ihigh,jhigh,ii,jj,k
	double precision Mh1in,Mh2in,Mh1,Mh2,ratio,Mhtemp
	double precision Mh1min,Mh1max,Mh2min,Mh2max
        double precision factBR,rvec(3),testtri
        double precision sep1,sep2,Mh1bit,Mh2bit,c1,c2,c3	
        double precision S95(0:imax,0:jmax)      

#include "tables-lep_cb.h"
		
	if(factBR .eq. 0d0) then
	ratio=-1d0
	return
	endif

!>> Ensure that Mh2 > Mh1 in order to make proper use of some 
!   of the tables.
!
!* tables 18 and 19:
!------------------
! table 18: process e+ e- -> Hi Hj with Hi, Hj -> b b
! table 19: process e+ e- -> Hi Hj with Hi, Hj -> tau tau
! Both cross sections are symmetric under the exchange m_Hi <--> m_Hj . 
! However, both tables are only given for Mh2 >= Mh1. 
! /!\: not anymore: >>Therefore, we swap the masses if Mh2in < Mh1in.<<
! Both permuatations of Higgs bosons will be checked anyway.
!  
!* tables 23 and 24:
!-------------------
! table 23: process e+ e- -> Hi Hj with, Hi -> b b, Hj -> tau tau 
!           and m_Hi <= m_Hj, i.e. Mh2 >= Mh1. 
!		Thus Mh2 -> m_Hj
!		     Mh1 -> m_Hi
! table 24: process e+ e- -> Hi Hj with, Hi -> tau tau, Hj -> b b
!           and m_Hi <= m_Hj, i.e. Mh2 >= Mh1.
!		Thus Mh2 -> m_Hj
!		     Mh1 -> m_Hi
! The external loop over n=1 ... NHIGSS**2 leads with 
!	i=findi(NHIGGS,n) and j=findj(NHIGGS,n) to a sum over all combinations
!	(i,j) in {(1,1), (1,2), (2,1), (2,2), (1,3), (3,1), (2,3), (3,2), ...}
! Thus in total one checks all possibilities
! for e+ e- -> Hi Hj with, Hi -> b b, Hj -> tau tau for all masses 
! m_Hi and m_Hj. Therefore, no swapping of masses is necessary.
!

	 Mh1=Mh1in
	 Mh2=Mh2in

	rvec(3)=0d0
	      	
! fill S95 from file
! row 0 and column 0 in LEP file contain higgs masses 
! and (0,0) ie top left set to -100	

       if(pred)then
	if(selec .eq. LEPh2Z_h1h1Z_4bZ) call eval_2d_arr(0,imax,0,jmax,S95,S95_15_pred)
        if(selec .eq. LEPh2Z_h1h1Z_4tauZ) call eval_2d_arr(0,imax,0,jmax,S95,S95_16_pred)
        if(selec .eq. LEPh2h1_4b) call eval_2d_arr(0,imax,0,jmax,S95,S95_18_pred)
        if(selec .eq. LEPh2h1_4tau) call eval_2d_arr(0,imax,0,jmax,S95,S95_19_pred)
        if(selec .eq. LEPh2h1_h1h1h1_6b) call eval_2d_arr(0,imax,0,jmax,S95,S95_20_pred)
        if(selec .eq. LEPh2h1_h1h1h1_6tau) call eval_2d_arr(0,imax,0,jmax,S95,S95_21_pred)
        if(selec .eq. LEPh2Z_h1h1Z_2b2tau) call eval_2d_arr(0,imax,0,jmax,S95,S95_22_pred)
        if(selec .eq. LEPh2h1_2b2tau) call eval_2d_arr(0,imax,0,jmax,S95,S95_23_pred)
        if(selec .eq. LEPh2h1_2tau2b) call eval_2d_arr(0,imax,0,jmax,S95,S95_24_pred)
       else
	if(selec .eq. LEPh2Z_h1h1Z_4bZ) call eval_2d_arr(0,imax,0,jmax,S95,S95_15_obs)
        if(selec .eq. LEPh2Z_h1h1Z_4tauZ) call eval_2d_arr(0,imax,0,jmax,S95,S95_16_obs)
        if(selec .eq. LEPh2h1_4b) call eval_2d_arr(0,imax,0,jmax,S95,S95_18_obs)
        if(selec .eq. LEPh2h1_4tau) call eval_2d_arr(0,imax,0,jmax,S95,S95_19_obs)
        if(selec .eq. LEPh2h1_h1h1h1_6b) call eval_2d_arr(0,imax,0,jmax,S95,S95_20_obs)
        if(selec .eq. LEPh2h1_h1h1h1_6tau) call eval_2d_arr(0,imax,0,jmax,S95,S95_21_obs)
        if(selec .eq. LEPh2Z_h1h1Z_2b2tau) call eval_2d_arr(0,imax,0,jmax,S95,S95_22_obs)
        if(selec .eq. LEPh2h1_2b2tau) call eval_2d_arr(0,imax,0,jmax,S95,S95_23_obs)
        if(selec .eq. LEPh2h1_2tau2b) call eval_2d_arr(0,imax,0,jmax,S95,S95_24_obs)
       endif

        ratio=0d0
	
! check if mass is within mass range of table:
	if(Mh1 .lt. Mh1min)ratio=-1.0D0
	if(Mh1 .gt. Mh1max)ratio=-1.0D0
	if(Mh2 .lt. Mh2min)ratio=-1.0D0
	if(Mh2 .gt. Mh2max)ratio=-1.0D0

! if hi->hj+hj is needed, check it's kinematically possible:         
	if((selec .eq. LEPh2Z_h1h1Z_4bZ) .or. (selec .eq. LEPh2Z_h1h1Z_4tauZ) 
     &		.or. (selec .eq. LEPh2Z_h1h1Z_2b2tau) 
     &		.or. (selec .eq. LEPh2h1_h1h1h1_6b) 
     &		.or. (selec .eq. LEPh2h1_h1h1h1_6tau)) then
	  if(Mh2 .lt. 2d0*Mh1) ratio=-1d0
	endif
	 		
        if(ratio .ge. 0d0)then	
! Interpolation: 
! points in table are treated as corners of flat triangular surfaces 	

	 
	 jlow=1+int((Mh1-Mh1min)/sep1)
	 Mh1bit=(Mh1-Mh1min)/sep1-dble(jlow-1)
  	 	 	 
	 ilow=1+int((Mh2-Mh2min)/sep2)
	 Mh2bit=(Mh2-Mh2min)/sep2-dble(ilow-1)	          			


! check if interpolation is needed, 
! (if  Mh1bit=Mh2bit=0 (or almost) no interpolation is needed)

	 if((Mh1bit**2+Mh2bit**2) .le. 1d-7) then !tag 1
		rvec(3)=S95(ilow,jlow)
	 else !tag 1
	  
!      *******    find out which triangle required point r is in
!      * B * *    i.e. triangle of orientation 'A' or 'B'
!      * * A *
!      *******	    	    
	   
	 testtri=Mh2bit/sqrt(Mh1bit**2.0D0+Mh2bit**2.0D0)
   	  	  
	  if(testtri.gt.0.7071D0)then !triangle orientation A
! check that all of the corners of the triangle are positive	 
	     if(S95(ilow+1,jlow).lt.0.0D0)then
	      rvec(3)=-6.0D0	  	   
             elseif(S95(ilow+1,jlow+1).lt.0.0D0)then
	      rvec(3)=-6.0D0	   
             elseif(S95(ilow,jlow).lt.0.0D0)then
	      rvec(3)=-6.0D0 
	     else
! find 'height' (i.e. S95 value) of each corner of triangle
	   
!             c3*       
!               * *          orth sides have length=1
!  rvec(2)^     * r * 
!         |     *     *      rvec(3) is S95_interpol
!         v   c1*********c2
!               <->
!             rvec(1)   
	   
	      c1=S95(ilow+1,jlow)
	      c2=S95(ilow+1,jlow+1)
	      c3=S95(ilow,jlow)
! find position of r relative to corner c1	    
	      rvec(1)=Mh1bit
	      rvec(2)=1.0D0-Mh2bit
! find 'height' of triangle at point r      
	      rvec(3)=rvec(1)*(c2-c1)+rvec(2)*(c3-c1)+c1
	    	    	    
	     endif
	  else !triangle orientation B
! check that all of the corners of the triangle are positive	 
	     if(S95(ilow,jlow+1).lt.0.0D0)then
	      rvec(3)=-6.0D0	  	   
             elseif(S95(ilow+1,jlow+1).lt.0.0D0)then
	      rvec(3)=-6.0D0	   
             elseif(S95(ilow,jlow).lt.0.0D0)then
	      rvec(3)=-6.0D0 
	     else
! find 'height' (i.e. S95 value) of each corner of triangle	   
	      c1=S95(ilow,jlow+1)
	      c2=S95(ilow,jlow)
	      c3=S95(ilow+1,jlow+1)
! find position of r relative to corner c1	    
	      rvec(1)=1.0D0-Mh1bit
	      rvec(2)=Mh2bit
! find 'height' of triangle at point r	      
	      rvec(3)=rvec(1)*(c2-c1)+rvec(2)*(c3-c1)+c1
	    	    	    		   
	  endif 	   	 	  
	 endif	!tag 1
	endif
	
	 if(rvec(3) .lt. 0d0)then
	   ratio = 0D0	  
	 else
	   ratio = factBR/rvec(3)
	 endif
	endif 

	end

