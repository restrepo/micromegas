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

!***********************************************************
	subroutine initialize_TEV_tables(NHIGGS)
!***********************************************************
	implicit none
	integer NHIGGS,i,j,m,n,imax,jmax,ntot
	real*8 Mhmin,Mhmax,sep,mh1min,mh2min,mh1max,mh2max,sep1,sep2
	real*8 mh1dummy,mh2dummy
#include "tables-tev_cb.h"

! fit function for BR_SM(H->bb) needed to convert tb-limits in sigma x BR-ratios
        double precision SMBR_Hbb

	write(*,*)'Tevatron-Constraints: initialization of data tables'

! tlist_tev(n) will tell the subroutines check_TEV_bound1
! which tables to look at
	do m=0,TEV_NSCHN-1,1
	do n=1,NHIGGS,1
	 tlist_tev(n+NHIGGS*m)=m+1
	enddo
	enddo
	do m=0,TEV_NDCHN-1,1
	do n=1,NHIGGS**2
	 tlist_tev(n+m*NHIGGS**2+TEV_NSCHN*NHIGGS)=m+1+TEV_NSCHN*NHIGGS
	enddo				
	enddo

! start values for Mhmin_overall/Mhmax_overall (min/max higher/lower than any realistic value)
	Mhmin_overall=10000d0
	Mhmax_overall=0d0


! initialize tables for the different channels

! 1) channel: ZH -> llbb
	write(*,*)'reading tables for channel ZH->llbb'
!>> obsolete
! c) D0 with 2.3fb^-1 (D0 note 5570)
! a) CDF with 2.4 fb^-1 (CDF note 9475)
! a) CDF with 2.7 fb^-1 (www-cdf.fnal.gov/physics/new/hdg/results/zhllbb_081024/index.htm) CDF_ZH_llbb_3_*
! b) CDF with 1 fb^-1 published (hep-ex/0807.4493)
!<<

! a) CDF with 4.1 fb^-1 (CDF note 9889)
	write(*,*)'  CDF with 4.1 fb^-1 (CDF note 9889)'
! information about how data is stored in file:	 
	Mhmin_CDF_ZH_llbb_4=100d0
	Mhmax_CDF_ZH_llbb_4=150d0
	sep_CDF_ZH_llbb_4=5d0

	Mhmin=Mhmin_CDF_ZH_llbb_4
	Mhmax=Mhmax_CDF_ZH_llbb_4
        sep=sep_CDF_ZH_llbb_4
	ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 11. 

	open(10,file=INSTDIR//"/TEVtables/CDF_ZH_llbb_4_exp.txt")
	open(11,file=INSTDIR//"/TEVtables/CDF_ZH_llbb_4_obs.txt")

        do i=1,ntot
	  read(10,*)Mh1dummy,CDF_ZH_llbb_4_pred(i)
          read(11,*)Mh1dummy,CDF_ZH_llbb_4_obs(i)
        enddo

        close(10)	 	
	close(11)

	Mhmin_overall=min(Mhmin_overall,Mhmin)
        Mhmax_overall=max(Mhmax_overall,Mhmax)


! b) CDF with 2.7 fb^-1 published (hep-ex/0908.3534)
	write(*,*)'  CDF table with 2.7 fb^-1 published (hep-ex/0908.3534)'
! information about how data is stored in file:	 
	Mhmin_CDF_ZH_llbb_pub_2=100d0
	Mhmax_CDF_ZH_llbb_pub_2=150d0
	sep_CDF_ZH_llbb_pub_2=5d0

	Mhmin=Mhmin_CDF_ZH_llbb_pub_2
	Mhmax=Mhmax_CDF_ZH_llbb_pub_2
        sep=sep_CDF_ZH_llbb_pub_2
	ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 11. 

	open(10,file=INSTDIR//"/TEVtables/CDF_ZH_llbb_pub_2_exp.txt")
	open(11,file=INSTDIR//"/TEVtables/CDF_ZH_llbb_pub_2_obs.txt")


        do i=1,ntot
	  read(10,*)Mh1dummy,CDF_ZH_llbb_pub_2_pred(i)
          read(11,*)Mh1dummy,CDF_ZH_llbb_pub_2_obs(i)
        enddo

        close(10)	 	
	close(11)

	Mhmin_overall=min(Mhmin_overall,Mhmin)
        Mhmax_overall=max(Mhmax_overall,Mhmax)


! c) D0 with 4.2 fb^-1 (D0 note 5876)
	write(*,*)'  D0 table with 4.2 fb^-1 (D0 note 5876)'
! information about how data is stored in file:	 
	Mhmin_D0_ZH_llbb_3=100d0
	Mhmax_D0_ZH_llbb_3=150d0
	sep_D0_ZH_llbb_3=5d0

	Mhmin=Mhmin_D0_ZH_llbb_3
	Mhmax=Mhmax_D0_ZH_llbb_3
        sep=sep_D0_ZH_llbb_3
	ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 5. 


	open(10,file=INSTDIR//"/TEVtables/D0_ZH_llbb_3_exp.txt")
	open(11,file=INSTDIR//"/TEVtables/D0_ZH_llbb_3_obs.txt")


        do i=1,ntot
	  read(10,*)Mh1dummy,D0_ZH_llbb_3_pred(i)
          read(11,*)Mh1dummy,D0_ZH_llbb_3_obs(i)
        enddo

        close(10)	 	
	close(11)

	Mhmin_overall=min(Mhmin_overall,Mhmin)
        Mhmax_overall=max(Mhmax_overall,Mhmax)




!2) channel: WH/ZH -> bb + missing E_T
	write(*,*)'reading tables for channel WH/ZH -> bb + missing E_T'
!>> obsolete
!a) CDF with 1 fb^-1 (CDF note 9166)
!a) CDF with 2.1 fb^-1 (CDF note 9483)
!a) CDF with 2.1 fb^-1 (CDF note 9642) CDF_VH_bbmET_3_*
!<<

!a) CDF with 3.6 fb^-1 (CDF note 9891)
	write(*,*)'  CDF table with 3.6 fb^-1 (CDF note 9891)'

! information about how data is stored in file:	 
	Mhmin_CDF_VH_bbmET_4=100d0
	Mhmax_CDF_VH_bbmET_4=150d0
	sep_CDF_VH_bbmET_4=5d0

	Mhmin=Mhmin_CDF_VH_bbmET_4
	Mhmax=Mhmax_CDF_VH_bbmET_4
        sep=sep_CDF_VH_bbmET_4
	ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 11. 


	open(10,file=INSTDIR//"/TEVtables/CDF_VH_bbmET_4_exp.txt")
	open(11,file=INSTDIR//"/TEVtables/CDF_VH_bbmET_4_obs.txt")


        do i=1,ntot
	  read(10,*)Mh1dummy,CDF_VH_bbmET_4_pred(i)
          read(11,*)Mh1dummy,CDF_VH_bbmET_4_obs(i)
        enddo

        close(10)	 	
	close(11)

	Mhmin_overall=min(Mhmin_overall,Mhmin)
        Mhmax_overall=max(Mhmax_overall,Mhmax)


!b) D0 with 2.1 fb^-1 (D0 note 5586)
	write(*,*)'  D0 with 2.1 fb^-1 (D0 note 5586)'
! information about how data is stored in file:	 
	Mhmin_D0_VH_bbmET=105d0
	Mhmax_D0_VH_bbmET=145d0
	sep_D0_VH_bbmET=10d0

	Mhmin=Mhmin_D0_VH_bbmET
	Mhmax=Mhmax_D0_VH_bbmET
        sep=sep_D0_VH_bbmET
	ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 5. 


	open(10,file=INSTDIR//"/TEVtables/D0_VH_bbmET_exp.txt")
	open(11,file=INSTDIR//"/TEVtables/D0_VH_bbmET_obs.txt")


        do i=1,ntot
	  read(10,*)Mh1dummy,D0_VH_bbmET_pred(i)
          read(11,*)Mh1dummy,D0_VH_bbmET_obs(i)
        enddo

        close(10)	 	
	close(11)

	Mhmin_overall=min(Mhmin_overall,Mhmin)
        Mhmax_overall=max(Mhmax_overall,Mhmax)


!c) D0 with 0.93 fb^-1 published (hep-ex/0808.1266)
	write(*,*)'  D0 with 0.93 fb^-1 published (hep-ex/0808.1266)'
! information about how data is stored in file:	 
	Mhmin_D0_VH_bbmET_pub=105d0
	Mhmax_D0_VH_bbmET_pub=135d0
	sep_D0_VH_bbmET_pub=10d0

	Mhmin=Mhmin_D0_VH_bbmET_pub
	Mhmax=Mhmax_D0_VH_bbmET_pub
        sep=sep_D0_VH_bbmET_pub
	ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 4. 


	open(10,file=INSTDIR//"/TEVtables/D0_VH_bbmET_pub_exp.txt")
	open(11,file=INSTDIR//"/TEVtables/D0_VH_bbmET_pub_obs.txt")


        do i=1,ntot
	  read(10,*)Mh1dummy,D0_VH_bbmET_pub_pred(i)
          read(11,*)Mh1dummy,D0_VH_bbmET_pub_obs(i)
        enddo

        close(10)	 	
	close(11)

	Mhmin_overall=min(Mhmin_overall,Mhmin)
        Mhmax_overall=max(Mhmax_overall,Mhmax)


!d) CDF with 1 fb^-1 published (hep-ex/0802.0432)
	write(*,*)'  CDF with 1 fb^-1 published (hep-ex/0802.0432)'

! information about how data is stored in file:	 
	Mhmin_CDF_VH_bbmET_pub=110d0
	Mhmax_CDF_VH_bbmET_pub=140d0
	sep_CDF_VH_bbmET_pub=5d0

	Mhmin=Mhmin_CDF_VH_bbmET_pub
	Mhmax=Mhmax_CDF_VH_bbmET_pub
        sep=sep_CDF_VH_bbmET_pub
	ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 7. 

	open(10,file=INSTDIR//"/TEVtables/CDF_VH_bbmET_pub_exp.txt")
	open(11,file=INSTDIR//"/TEVtables/CDF_VH_bbmET_pub_obs.txt")

        do i=1,ntot
	  read(10,*)Mh1dummy,CDF_VH_bbmET_pub_pred(i)
          read(11,*)Mh1dummy,CDF_VH_bbmET_pub_obs(i)
        enddo

        close(10)	 	
	close(11)

	Mhmin_overall=min(Mhmin_overall,Mhmin)
        Mhmax_overall=max(Mhmax_overall,Mhmax)




!3) channel: WH -> lnubb
	write(*,*)'reading tables for channel WH->lnubb'
!>> obsolete
!b) D0 with 1.7 fb^-1 (D0 note 5472)
!c) CDF with 2.7 fb^-1 (CDF note 9463)
!b) D0 with 2.7 fb^-1 (D0 note 5828) D0_WH_lnubb_2_*
!e) CDF with 1.0 fb^-1 published (hep-ex/0803.3493, PRD 78(2008)032008) 
!c) CDF with 2.7 fb^-1 (CDF note 9596) CDF_WH_lnubb_4_*
!<<

!b) D0 with 5.0 fb^-1 (D0 note 5972) 
	write(*,*)'  D0 table with 5.0 fb^-1 (D0 note 5972)'
! information about how data is stored in file:	 
	Mhmin_D0_WH_lnubb_3=100d0
	Mhmax_D0_WH_lnubb_3=150d0
	sep_D0_WH_lnubb_3=5d0

	Mhmin=Mhmin_D0_WH_lnubb_3
	Mhmax=Mhmax_D0_WH_lnubb_3
        sep=sep_D0_WH_lnubb_3
	ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 11

	open(10,file=INSTDIR//"/TEVtables/D0_WH_lnubb_3_exp.txt")
	open(11,file=INSTDIR//"/TEVtables/D0_WH_lnubb_3_obs.txt")

        do i=1,ntot
	  read(10,*)Mh1dummy,D0_WH_lnubb_3_pred(i)
          read(11,*)Mh1dummy,D0_WH_lnubb_3_obs(i)
        enddo

        close(10)	 	
	close(11)

	Mhmin_overall=min(Mhmin_overall,Mhmin)
        Mhmax_overall=max(Mhmax_overall,Mhmax)


!c) CDF with 4.3 fb^-1 (http://www-cdf.fnal.gov/physics/new/hdg/results/whlnubb_090814/WH4.3fb.html)
	write(*,*)'  CDF table with 4.3 fb^-1 (http://www-cdf.fnal.gov/physics/new/hdg/results/whlnubb_090814/WH4.3fb.html)'
! information about how data is stored in file:	 
	Mhmin_CDF_WH_lnubb_5=100d0
	Mhmax_CDF_WH_lnubb_5=150d0
	sep_CDF_WH_lnubb_5=5d0

	Mhmin=Mhmin_CDF_WH_lnubb_5
	Mhmax=Mhmax_CDF_WH_lnubb_5
        sep=sep_CDF_WH_lnubb_5
	ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 11. 

	open(10,file=INSTDIR//"/TEVtables/CDF_WH_lnubb_5_exp.txt")
	open(11,file=INSTDIR//"/TEVtables/CDF_WH_lnubb_5_obs.txt")

        do i=1,ntot
	  read(10,*)Mh1dummy,CDF_WH_lnubb_5_pred(i)
          read(11,*)Mh1dummy,CDF_WH_lnubb_5_obs(i)
        enddo

        close(10)	 	
	close(11)

	Mhmin_overall=min(Mhmin_overall,Mhmin)
        Mhmax_overall=max(Mhmax_overall,Mhmax)


!d) D0 with 1.1 fb^-1 published (hep-ex/0808.1970, PRL 102(2009)051803) 
	write(*,*)'  D0 table with 1.1 fb^-1 published (hep-ex/0808.1970)'
! information about how data is stored in file:	 
	Mhmin_D0_WH_lnubb_pub=100d0
	Mhmax_D0_WH_lnubb_pub=150d0
	sep_D0_WH_lnubb_pub=5d0

	Mhmin=Mhmin_D0_WH_lnubb_pub
	Mhmax=Mhmax_D0_WH_lnubb_pub
        sep=sep_D0_WH_lnubb_pub
	ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 11

	open(10,file=INSTDIR//"/TEVtables/D0_WH_lnubb_pub_exp.txt")
	open(11,file=INSTDIR//"/TEVtables/D0_WH_lnubb_pub_obs.txt")

        do i=1,ntot
	  read(10,*)Mh1dummy,D0_WH_lnubb_pub_pred(i)
          read(11,*)Mh1dummy,D0_WH_lnubb_pub_obs(i)
        enddo

        close(10)	 	
	close(11)

	Mhmin_overall=min(Mhmin_overall,Mhmin)
        Mhmax_overall=max(Mhmax_overall,Mhmax)


!e) CDF with 2.7 fb^-1 published (hep-ex/0906.5613)
	write(*,*)'  CDF table with 2.7 fb^-1 published (hep-ex/0906.5613)'
! information about how data is stored in file:	 
	Mhmin_CDF_WH_lnubb_pub_2=100d0
	Mhmax_CDF_WH_lnubb_pub_2=150d0
	sep_CDF_WH_lnubb_pub_2=5d0

	Mhmin=Mhmin_CDF_WH_lnubb_pub_2
	Mhmax=Mhmax_CDF_WH_lnubb_pub_2
        sep=sep_CDF_WH_lnubb_pub_2
	ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 11

	open(10,file=INSTDIR//"/TEVtables/CDF_WH_lnubb_pub_2_exp.txt")
	open(11,file=INSTDIR//"/TEVtables/CDF_WH_lnubb_pub_2_obs.txt")

        do i=1,ntot
	  read(10,*)Mh1dummy,CDF_WH_lnubb_pub_2_pred(i)
          read(11,*)Mh1dummy,CDF_WH_lnubb_pub_2_obs(i)
        enddo

        close(10)	 	
	close(11)

	Mhmin_overall=min(Mhmin_overall,Mhmin)
        Mhmax_overall=max(Mhmax_overall,Mhmax)





!4) channel: ppH -> WW
	write(*,*)'reading tables for channel ppH -> WW'
!>> obsolete
!a) CDF 1.9 fb^-1
!b) D0 note 5489
!c) D0  with 2.3 fb^-1 (D0 note 5624)
!d) CDF 2.4 fb^-1 (CDF note 9236)
!g) CDF 3.6 fb^-1 (www-cdf.fnal.gov/physics/new/hdg/results/hwwmenn_090227/)
!<< obsolete

!e) D0 3.0 fb^-1
	write(*,*)'  D0 table for ppH -> WW -> ll with 3.0 fb^-1 (D0 note 5757)'
! information about how data is stored in file:	 
	Mhmin_D0_ppH_WW_ll_2=115d0
	Mhmax_D0_ppH_WW_ll_2=200d0
	sep_D0_ppH_WW_ll_2=5d0

	Mhmin=Mhmin_D0_ppH_WW_ll_2
	Mhmax=Mhmax_D0_ppH_WW_ll_2
        sep=sep_D0_ppH_WW_ll_2
	ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 18. 

	open(10,file=INSTDIR//"/TEVtables/D0_ppH_WW_ll_2_exp.txt")
	open(11,file=INSTDIR//"/TEVtables/D0_ppH_WW_ll_2_obs.txt")

        do i=1,ntot
	  read(10,*)Mh1dummy,D0_ppH_WW_ll_2_pred(i)
          read(11,*)Mh1dummy,D0_ppH_WW_ll_2_obs(i)
        enddo

        close(10)	 	
	close(11)

	Mhmin_overall=min(Mhmin_overall,Mhmin)
        Mhmax_overall=max(Mhmax_overall,Mhmax)


!f) CDF 3.0 fb^-1
	write(*,*)'  CDF table for 3.0 fb^-1 (hep-ex/0809.3930)'
! information about how data is stored in file:	 
	Mhmin_CDF_ppH_WW_4=110d0
	Mhmax_CDF_ppH_WW_4=200d0
	sep_CDF_ppH_WW_4=10d0

	Mhmin=Mhmin_CDF_ppH_WW_4
	Mhmax=Mhmax_CDF_ppH_WW_4
        sep=sep_CDF_ppH_WW_4
	ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 10. 

	open(10,file=INSTDIR//"/TEVtables/CDF_ppH_WW_4_exp.txt")
	open(11,file=INSTDIR//"/TEVtables/CDF_ppH_WW_4_obs.txt")

        do i=1,ntot
	  read(10,*)Mh1dummy,CDF_ppH_WW_4_pred(i)
          read(11,*)Mh1dummy,CDF_ppH_WW_4_obs(i)
        enddo

        close(10)	 	
	close(11)

	Mhmin_overall=min(Mhmin_overall,Mhmin)
        Mhmax_overall=max(Mhmax_overall,Mhmax)


!g) CDF 4.8 fb^-1 (CDF note 9887, www-cdf.fnal.gov/physics/new/hdg/results/hwwmenn_090814/)
	write(*,*)'  CDF with 4.8 fb^-1 (CDF note 9887)'
! information about how data is stored in file:	 
	Mhmin_CDF_ppH_WW_6=110d0
	Mhmax_CDF_ppH_WW_6=200d0
	sep_CDF_ppH_WW_6=5d0

	Mhmin=Mhmin_CDF_ppH_WW_6
	Mhmax=Mhmax_CDF_ppH_WW_6
        sep=sep_CDF_ppH_WW_6
	ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 19. 

	open(10,file=INSTDIR//"/TEVtables/CDF_ppH_WW_6_exp.txt")
	open(11,file=INSTDIR//"/TEVtables/CDF_ppH_WW_6_obs.txt")

        do i=1,ntot
	  read(10,*)Mh1dummy,CDF_ppH_WW_6_pred(i)
          read(11,*)Mh1dummy,CDF_ppH_WW_6_obs(i)
        enddo

        close(10)	 	
	close(11)

	Mhmin_overall=min(Mhmin_overall,Mhmin)
        Mhmax_overall=max(Mhmax_overall,Mhmax)


!h) D0 3.0-4.2 fb^-1
	write(*,*)'  D0 table for ppH -> WW + X -> ll + X with 3.0-4.2 fb^-1 (D0 note 5871)'
! information about how data is stored in file:	 
	Mhmin_D0_ppH_WW_ll_3=115d0
	Mhmax_D0_ppH_WW_ll_3=200d0
	sep_D0_ppH_WW_ll_3=5d0

	Mhmin=Mhmin_D0_ppH_WW_ll_3
	Mhmax=Mhmax_D0_ppH_WW_ll_3
        sep=sep_D0_ppH_WW_ll_3
	ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 18. 

	open(10,file=INSTDIR//"/TEVtables/D0_ppH_WW_ll_3_exp.txt")
	open(11,file=INSTDIR//"/TEVtables/D0_ppH_WW_ll_3_obs.txt")

        do i=1,ntot
	  read(10,*)Mh1dummy,D0_ppH_WW_ll_3_pred(i)
          read(11,*)Mh1dummy,D0_ppH_WW_ll_3_obs(i)
        enddo

        close(10)	 	
	close(11)

	Mhmin_overall=min(Mhmin_overall,Mhmin)
        Mhmax_overall=max(Mhmax_overall,Mhmax)





!5) channel: SM combined limits

!a) CDF & D0 SM combined limits of Dec 2007
	write(*,*)'reading tables for the CDF & D0 SM combined limits (hep-ex/0712.2383)'
! information about how data is stored in file:	 
	Mhmin_CDF_D0_SM=110d0
	Mhmax_CDF_D0_SM=200d0
	sep_CDF_D0_SM=5d0

	Mhmin=Mhmin_CDF_D0_SM
	Mhmax=Mhmax_CDF_D0_SM
        sep=sep_CDF_D0_SM
	ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 19. 

	open(10,file=INSTDIR//"/TEVtables/CDF_D0_SM_exp.txt")
	open(11,file=INSTDIR//"/TEVtables/CDF_D0_SM_obs.txt")

        do i=1,ntot
	  read(10,*)Mh1dummy,CDF_D0_SM_pred(i)
          read(11,*)Mh1dummy,CDF_D0_SM_obs(i)
        enddo

        close(10)	 	
	close(11)

	Mhmin_overall=min(Mhmin_overall,Mhmin)
        Mhmax_overall=max(Mhmax_overall,Mhmax)


!b) CDF H -> tau tau SM combined limits 
	write(*,*)'reading tables for the CDF H -> tau tau SM combined limits (CDF note 9248)'
! information about how data is stored in file:	 
	Mhmin_CDF_SM_tautau=110d0
	Mhmax_CDF_SM_tautau=150d0
	sep_CDF_SM_tautau=5d0

	Mhmin=Mhmin_CDF_SM_tautau
	Mhmax=Mhmax_CDF_SM_tautau
        sep=sep_CDF_SM_tautau
	ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 9. 

	open(10,file=INSTDIR//"/TEVtables/CDF_SM_tautau_exp.txt")
	open(11,file=INSTDIR//"/TEVtables/CDF_SM_tautau_obs.txt")

        do i=1,ntot
	  read(10,*)Mh1dummy,CDF_SM_tautau_pred(i)
          read(11,*)Mh1dummy,CDF_SM_tautau_obs(i)
        enddo

        close(10)	 	
	close(11)

	Mhmin_overall=min(Mhmin_overall,Mhmin)
        Mhmax_overall=max(Mhmax_overall,Mhmax)


!c) CDF & D0 SM combined limits of Apr 2008
	write(*,*)'reading tables for the CDF & D0 SM combined limits (hep-ex/0804.3423)'
! information about how data is stored in file:	 
	Mhmin_CDF_D0_SM_2=110d0
	Mhmax_CDF_D0_SM_2=200d0
	sep_CDF_D0_SM_2=5d0

	Mhmin=Mhmin_CDF_D0_SM_2
	Mhmax=Mhmax_CDF_D0_SM_2
        sep=sep_CDF_D0_SM_2
	ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 19. 

	open(10,file=INSTDIR//"/TEVtables/CDF_D0_SM_2_exp.txt")
	open(11,file=INSTDIR//"/TEVtables/CDF_D0_SM_2_obs.txt")

        do i=1,ntot
	  read(10,*)Mh1dummy,CDF_D0_SM_2_pred(i)
          read(11,*)Mh1dummy,CDF_D0_SM_2_obs(i)
        enddo

        close(10)	 	
	close(11)

	Mhmin_overall=min(Mhmin_overall,Mhmin)
        Mhmax_overall=max(Mhmax_overall,Mhmax)


!d) CDF & D0 SM combined limits of August 2008
	write(*,*)'reading tables for the CDF & D0 SM combined limits (hep-ex/0808.0534)'
! information about how data is stored in file:	 
	Mhmin_CDF_D0_SM_3=155d0
	Mhmax_CDF_D0_SM_3=200d0
	sep_CDF_D0_SM_3=5d0

	Mhmin=Mhmin_CDF_D0_SM_3
	Mhmax=Mhmax_CDF_D0_SM_3
        sep=sep_CDF_D0_SM_3
	ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 10. 
! If not, change appropriate declarations in 'tables-tev_cb.h'

	open(10,file=INSTDIR//"/TEVtables/CDF_D0_SM_3_exp.txt")
	open(11,file=INSTDIR//"/TEVtables/CDF_D0_SM_3_obs.txt")

        do i=1,ntot
	  read(10,*)Mh1dummy,CDF_D0_SM_3_pred(i)
          read(11,*)Mh1dummy,CDF_D0_SM_3_obs(i)
        enddo

        close(10)	 	
	close(11)

	Mhmin_overall=min(Mhmin_overall,Mhmin)
        Mhmax_overall=max(Mhmax_overall,Mhmax)


!d) CDF SM combined limit with 2.0 - 3.0 fb^-1
	write(*,*)'reading tables for the CDF SM combined limit (CDF note 9674)'
! information about how data is stored in file:	 
	Mhmin_CDF_SM=100d0
	Mhmax_CDF_SM=200d0
	sep_CDF_SM=5d0

	Mhmin=Mhmin_CDF_SM
	Mhmax=Mhmax_CDF_SM
        sep=sep_CDF_SM
	ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 21. 
! If not, change appropriate declarations in 'tables-tev_cb.h'

	open(10,file=INSTDIR//"/TEVtables/CDF_SM_exp.txt")
	open(11,file=INSTDIR//"/TEVtables/CDF_SM_obs.txt")

        do i=1,ntot
	  read(10,*)Mh1dummy,CDF_SM_pred(i)
          read(11,*)Mh1dummy,CDF_SM_obs(i)
        enddo

        close(10)	 	
	close(11)

	Mhmin_overall=min(Mhmin_overall,Mhmin)
        Mhmax_overall=max(Mhmax_overall,Mhmax)


!e) D0 SM combined limit with 0.44 fb^-1, published
	write(*,*)'reading tables for the D0 SM combined limit (hep-ex/0712.0598)'
! information about how data is stored in file:	 
	Mhmin_D0_SM_pub=100d0
	Mhmax_D0_SM_pub=200d0
	sep_D0_SM_pub=5d0

	Mhmin=Mhmin_D0_SM_pub
	Mhmax=Mhmax_D0_SM_pub
        sep=sep_D0_SM_pub
	ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 21. 
! If not, change appropriate declarations in 'tables-tev_cb.h'

	open(10,file=INSTDIR//"/TEVtables/D0_SM_pub_exp.txt")
	open(11,file=INSTDIR//"/TEVtables/D0_SM_pub_obs.txt")

        do i=1,ntot
	  read(10,*)Mh1dummy,D0_SM_pub_pred(i)
          read(11,*)Mh1dummy,D0_SM_pub_obs(i)
        enddo

        close(10)	 	
	close(11)

	Mhmin_overall=min(Mhmin_overall,Mhmin)
        Mhmax_overall=max(Mhmax_overall,Mhmax)

!f) CDF & D0  combined limit with 0.9-4.2 fb^-1
	write(*,*)'reading tables for the CDF & D0 SM combined limits for m_H > 155 GeV (hep-ex/0903.4001)'
! information about how data is stored in file:	 
	Mhmin_CDF_D0_SM_4=155d0
	Mhmax_CDF_D0_SM_4=200d0
	sep_CDF_D0_SM_4=5d0

	Mhmin=Mhmin_CDF_D0_SM_4
	Mhmax=Mhmax_CDF_D0_SM_4
        sep=sep_CDF_D0_SM_4
	ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 10. 
! If not, change appropriate declarations in 'tables-tev_cb.h'

	open(10,file=INSTDIR//"/TEVtables/CDF_D0_SM_4_exp.txt")
	open(11,file=INSTDIR//"/TEVtables/CDF_D0_SM_4_obs.txt")

        do i=1,ntot
	  read(10,*)Mh1dummy,CDF_D0_SM_4_pred(i)
          read(11,*)Mh1dummy,CDF_D0_SM_4_obs(i)
        enddo

        close(10)	 	
	close(11)

	Mhmin_overall=min(Mhmin_overall,Mhmin)
        Mhmax_overall=max(Mhmax_overall,Mhmax)

!g) CDF SM combined limit with 2.0 - 4.8 fb^-1
	write(*,*)'reading tables for the CDF SM combined limit (CDF note 9897)'
! information about how data is stored in file:	 
	Mhmin_CDF_SM_2=100d0
	Mhmax_CDF_SM_2=200d0
	sep_CDF_SM_2=5d0

	Mhmin=Mhmin_CDF_SM_2
	Mhmax=Mhmax_CDF_SM_2
        sep=sep_CDF_SM_2
	ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 21. 
! If not, change appropriate declarations in 'tables-tev_cb.h'

	open(10,file=INSTDIR//"/TEVtables/CDF_SM_2_exp.txt")
	open(11,file=INSTDIR//"/TEVtables/CDF_SM_2_obs.txt")

        do i=1,ntot
	  read(10,*)Mh1dummy,CDF_SM_2_pred(i)
          read(11,*)Mh1dummy,CDF_SM_2_obs(i)
        enddo

        close(10)	 	
	close(11)

	Mhmin_overall=min(Mhmin_overall,Mhmin)
        Mhmax_overall=max(Mhmax_overall,Mhmax)




!6) channel: p p -> b H, H-> b b
	write(*,*)'reading tables for pp -> bH, H -> bb'
!a) CDF
	write(*,*)'  CDF tanb-table (CDF note 9284)'
! information about how data is stored in file:	 
	Mhmin_CDF_bH_bb=90d0
	Mhmax_CDF_bH_bb=210d0
	sep_CDF_bH_bb=10d0

	Mhmin=Mhmin_CDF_bH_bb
	Mhmax=Mhmax_CDF_bH_bb
        sep=sep_CDF_bH_bb
	ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 13. 

	open(10,file=INSTDIR//"/TEVtables/CDF_bH_bb_tb_limit_exp.txt")
	open(11,file=INSTDIR//"/TEVtables/CDF_bH_bb_tb_limit_obs.txt")

! reading tb-limits and converting to (sigma x BR)_limit/(sigma x BR)_SM
! see CDF note 9284
!   (sigma x BR)_limit                                           0.9
!  --------------------(m_H) = ( tanb_limit(m_H) )^2 x 2 x ------------------
!     (sigma x BR)_SM                                       BR_SM(H->bb,m_H)
!
! As this function is steeply rising with m_H it is not ideal for interpolating 
! between given points. We choose here (*)
!
!   (sigma x BR)_limit    
!  --------------------(m_H) = ( tanb_limit(m_H) )^2 x 2 x 0.9
!        sigma_SM                                       

        do i=1,ntot
	  read(10,*)Mh1dummy,CDF_bH_bb_pred(i)
	  CDF_bH_bb_pred(i)=CDF_bH_bb_pred(i)**2*2d0*0.9d0
!left out for choice (*) : /SMBR_Hbb(Mh1dummy)
          read(11,*)Mh1dummy,CDF_bH_bb_obs(i)
          CDF_bH_bb_obs(i)=CDF_bH_bb_obs(i)**2*2d0*0.9d0
!left out for choice (*) : /SMBR_Hbb(Mh1dummy)
        enddo

        close(10)	 	
	close(11)

	Mhmin_overall=min(Mhmin_overall,Mhmin)
        Mhmax_overall=max(Mhmax_overall,Mhmax)


!>> removed
!!b) CDF absolute sigma x BR table in pikobarn for 1.9fb^-1 (CDF note 9284)
!<< removed

!c) D0
	write(*,*)'  D0 tanb-table (D0 note 5726)'
! information about how data is stored in file:	 
	Mhmin_D0_bH_bb_tb_limit_2=90d0
	Mhmax_D0_bH_bb_tb_limit_2=220d0
	sep_D0_bH_bb_tb_limit_2=10d0

	Mhmin=Mhmin_D0_bH_bb_tb_limit_2
	Mhmax=Mhmax_D0_bH_bb_tb_limit_2
        sep=sep_D0_bH_bb_tb_limit_2
	ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 14. 

	open(10,file=INSTDIR//"/TEVtables/D0_bH_bb_tb_limit_2_exp.txt")
	open(11,file=INSTDIR//"/TEVtables/D0_bH_bb_tb_limit_2_obs.txt")

! reading tb-limits and converting to (sigma x BR)_limit/(sigma x BR)_SM
! see D0 note 5726 (currently unclear if factor of 2 below is right)
!   (sigma x BR)_limit                                           0.9
!  --------------------(m_H) = ( tanb_limit(m_H) )^2 x 2 x ------------------
!     (sigma x BR)_SM                                       BR_SM(H->bb,m_H)
!
! As this function is steeply rising with m_H it is not ideal for interpolating 
! between given points. We choose here (*):
!
!   (sigma x BR)_limit    
!  --------------------(m_H) = ( tanb_limit(m_H) )^2 x 2 x 0.9
!        sigma_SM                                       
        do i=1,ntot
	  read(10,*)Mh1dummy,D0_bH_bb_tb_limit_2_pred(i)
	  D0_bH_bb_tb_limit_2_pred(i)=D0_bH_bb_tb_limit_2_pred(i)**2*2d0*0.9d0
!left out for choice (*) : /SMBR_Hbb(Mh1dummy)
          read(11,*)Mh1dummy,D0_bH_bb_tb_limit_2_obs(i)
          D0_bH_bb_tb_limit_2_obs(i)=D0_bH_bb_tb_limit_2_obs(i)**2*2d0*0.9d0
!left out for choice (*) : /SMBR_Hbb(Mh1dummy)
        enddo

        close(10)	 	
	close(11)

	Mhmin_overall=min(Mhmin_overall,Mhmin)
        Mhmax_overall=max(Mhmax_overall,Mhmax)


!d) D0
	write(*,*)'  D0 absolute sigma x BR table in pikobarn for 1fb^-1 (hep-ex/0805.3556)'
! information about how data is stored in file:	 
	Mhmin_D0_bH_bb_2=90d0
	Mhmax_D0_bH_bb_2=220d0
	sep_D0_bH_bb_2=10d0

	Mhmin=Mhmin_D0_bH_bb_2
	Mhmax=Mhmax_D0_bH_bb_2
        sep=sep_D0_bH_bb_2
	ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 14. 

	open(10,file=INSTDIR//"/TEVtables/D0_bH_bb-abs-limit_exp.txt")
	open(11,file=INSTDIR//"/TEVtables/D0_bH_bb-abs-limit_obs.txt")

        do i=1,ntot
	  read(10,*)Mh1dummy,D0_bH_bb_2_pred(i)
          read(11,*)Mh1dummy,D0_bH_bb_2_obs(i)
        enddo

        close(10)	 	
	close(11)

	Mhmin_overall=min(Mhmin_overall,Mhmin)
        Mhmax_overall=max(Mhmax_overall,Mhmax)





!7) channel: ppH -> tau tau
!>> obsolete
!b) CDF absolute sigma x BR table in pikobarn for 1.8 fb^-1 (CDF note 9071)
!<<
	write(*,*)'reading tables for channel ppH -> tau tau'

!a) D0 1 fb^-1
	write(*,*)'  D0 absolute sigma x BR table in pikobarn for 1 fb^-1 (hep-ex/0805.2491)'
! information about how data is stored in file:	 
	Mhmin_D0_ppH_tautau_pub=90d0
	Mhmax_D0_ppH_tautau_pub=300d0
	sep_D0_ppH_tautau_pub=10d0

	Mhmin=Mhmin_D0_ppH_tautau_pub
	Mhmax=Mhmax_D0_ppH_tautau_pub
        sep=sep_D0_ppH_tautau_pub
	ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 22. 

	open(10,file=INSTDIR//"/TEVtables/D0_ppH_tautau_pub-abs-limit_exp.txt")
	open(11,file=INSTDIR//"/TEVtables/D0_ppH_tautau_pub-abs-limit_obs.txt")

        do i=1,ntot
	  read(10,*)Mh1dummy,D0_ppH_tautau_pub_pred(i)
          read(11,*)Mh1dummy,D0_ppH_tautau_pub_obs(i)
        enddo

        close(10)	 	
	close(11)

	Mhmin_overall=min(Mhmin_overall,Mhmin)
        Mhmax_overall=max(Mhmax_overall,Mhmax)


!b) CDF absolute sigma x BR table in pikobarn for 1.8 fb^-1 published (CDF0906_1014)
	write(*,*)'  CDF absolute sigma x BR table in pikobarn for 1.8 fb^-1 published (hep-ex/0906.1014)'
! information about how data is stored in file:	 
	Mhmin_CDF_ppH_tautau_pub=90d0
	Mhmax_CDF_ppH_tautau_pub=250d0
	sep_CDF_ppH_tautau_pub=10d0

	Mhmin=Mhmin_CDF_ppH_tautau_pub
	Mhmax=Mhmax_CDF_ppH_tautau_pub
        sep=sep_CDF_ppH_tautau_pub
	ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 17. 

	open(10,file=INSTDIR//"/TEVtables/CDF_ppH_tautau-abs-limit_pub_exp.txt")
	open(11,file=INSTDIR//"/TEVtables/CDF_ppH_tautau-abs-limit_pub_obs.txt")

        do i=1,ntot
	  read(10,*)Mh1dummy,CDF_ppH_tautau_pub_pred(i)
          read(11,*)Mh1dummy,CDF_ppH_tautau_pub_obs(i)
        enddo

        close(10)	 	
	close(11)

	Mhmin_overall=min(Mhmin_overall,Mhmin)
        Mhmax_overall=max(Mhmax_overall,Mhmax)

!c) D0 absolute sigma x BR table in pikobarn for 2.2 fb^-1 (D0 note 5740)
        write(*,*)'  D0 absolute sigma x BR table in pikobarn for 2.2 fb^-1 (D0 note 5740)'
! information about how data is stored in file:  
        Mhmin_D0_ppH_tautau=90d0
        Mhmax_D0_ppH_tautau=300d0
        sep_D0_ppH_tautau=10d0

        Mhmin=Mhmin_D0_ppH_tautau
        Mhmax=Mhmax_D0_ppH_tautau
        sep=sep_D0_ppH_tautau
        ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 22. 

        open(10,file=INSTDIR//"/TEVtables/D0_ppH_tautau-abs-limit_exp.txt")
        open(11,file=INSTDIR//"/TEVtables/D0_ppH_tautau-abs-limit_obs.txt")

        do i=1,ntot
          read(10,*)Mh1dummy,D0_ppH_tautau_pred(i)
          read(11,*)Mh1dummy,D0_ppH_tautau_obs(i)
        enddo

        close(10)               
        close(11)

        Mhmin_overall=min(Mhmin_overall,Mhmin)
        Mhmax_overall=max(Mhmax_overall,Mhmax)


!d) CDF & D0 combined limits for ppH -> tau tau, absolute sigma x BR table in pikobarn 
!   for 1.8 fb^-1 (CDF note 9071) & 2.2 fb^-1 (D0 note 5740), CDF note 9888, D0 note 5980
	write(*,*)'  CDF  & D0 combined limits (absolute sigma x BR in pikobarn) (CDF note 9888, D0 note 5980)'
! information about how data is stored in file:  
        Mhmin_CDF_D0_ppH_tautau=100d0
        Mhmax_CDF_D0_ppH_tautau=200d0
        sep_CDF_D0_ppH_tautau=10d0

        Mhmin=Mhmin_CDF_D0_ppH_tautau
        Mhmax=Mhmax_CDF_D0_ppH_tautau
        sep=sep_CDF_D0_ppH_tautau
        ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 11. 

        open(10,file=INSTDIR//"/TEVtables/CDF_D0_ppH_tautau-abs-limit_exp.txt")
        open(11,file=INSTDIR//"/TEVtables/CDF_D0_ppH_tautau-abs-limit_obs.txt")

        do i=1,ntot
          read(10,*)Mh1dummy,CDF_D0_ppH_tautau_pred(i)
          read(11,*)Mh1dummy,CDF_D0_ppH_tautau_obs(i)
        enddo

        close(10)               
        close(11)

        Mhmin_overall=min(Mhmin_overall,Mhmin)
        Mhmax_overall=max(Mhmax_overall,Mhmax)




!8) channel: WH -> WWW
	write(*,*)'reading tables for channel WH->WWW'
!>> obslete
!b) CDF with 1.9 fb^-1 (CDF note 7307 version 2)
!a) D0 absolute sigma x BR table in pikobarn with 1 fb^-1 (D0 note 5485)
!<<

!a) D0 with 3.6 fb^-1 (D0 note 5873)
	write(*,*)'  D0 with 3.6 fb^-1 (D0 note 5873)'
! information about how data is stored in file:	 
	Mhmin_D0_WH_WWW_2=120d0
	Mhmax_D0_WH_WWW_2=200d0
	sep_D0_WH_WWW_2=20d0

	Mhmin=Mhmin_D0_WH_WWW_2
	Mhmax=Mhmax_D0_WH_WWW_2
        sep=sep_D0_WH_WWW_2
	ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 5. 

	open(10,file=INSTDIR//"/TEVtables/D0_WH_WWW_2_exp.txt")
	open(11,file=INSTDIR//"/TEVtables/D0_WH_WWW_2_obs.txt")

        do i=1,ntot
	  read(10,*)Mh1dummy,D0_WH_WWW_2_pred(i)
          read(11,*)Mh1dummy,D0_WH_WWW_2_obs(i)
        enddo

        close(10)	 	
	close(11)

	Mhmin_overall=min(Mhmin_overall,Mhmin)
        Mhmax_overall=max(Mhmax_overall,Mhmax)


!b) CDF with 2.7 fb^-1 (CDF note 7307 version 3)
	write(*,*)'  CDF with 2.7 fb^-1 (CDF note 7307 version 3)'
! information about how data is stored in file:	 
	Mhmin_CDF_WH_WWW_2=110d0
	Mhmax_CDF_WH_WWW_2=200d0
	sep_CDF_WH_WWW_2=10d0

	Mhmin=Mhmin_CDF_WH_WWW_2
	Mhmax=Mhmax_CDF_WH_WWW_2
        sep=sep_CDF_WH_WWW_2
	ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 10. 

	open(10,file=INSTDIR//"/TEVtables/CDF_WH_WWW_2_exp.txt")
	open(11,file=INSTDIR//"/TEVtables/CDF_WH_WWW_2_obs.txt")

        do i=1,ntot
	  read(10,*)Mh1dummy,CDF_WH_WWW_2_pred(i)
          read(11,*)Mh1dummy,CDF_WH_WWW_2_obs(i)
        enddo

        close(10)	 	
	close(11)

	Mhmin_overall=min(Mhmin_overall,Mhmin)
        Mhmax_overall=max(Mhmax_overall,Mhmax)




!9) channel: ppH -> gamma gamma
	write(*,*)'reading tables for channel: ppH -> gamma gamma'
!>> obsolete
!a) D0 with 2.68 fb^-1 (D0 note 5737)
!<< obsolete

!a) D0 with 4.2 fb^-1
	write(*,*)'  D0 with 4.2 fb^-1 (D0 note 5858)'
! information about how data is stored in file:	 
	Mhmin_D0_ppH_gamgam_4=100d0
	Mhmax_D0_ppH_gamgam_4=150d0
	sep_D0_ppH_gamgam_4=2.5d0

	Mhmin=Mhmin_D0_ppH_gamgam_4
	Mhmax=Mhmax_D0_ppH_gamgam_4
        sep=sep_D0_ppH_gamgam_4
	ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 21. 

	open(10,file=INSTDIR//"/TEVtables/D0_ppH_gamgam_4_exp.txt")
	open(11,file=INSTDIR//"/TEVtables/D0_ppH_gamgam_4_obs.txt")

        do i=1,ntot
	  read(10,*)Mh1dummy,D0_ppH_gamgam_4_pred(i)
          read(11,*)Mh1dummy,D0_ppH_gamgam_4_obs(i)
        enddo

        close(10)	 	
	close(11)

	Mhmin_overall=min(Mhmin_overall,Mhmin)
        Mhmax_overall=max(Mhmax_overall,Mhmax)



!b) D0 with 1.1 fb^-1
	write(*,*)'  D0 absolute sigma x BR table in pikobarn for 1.1 fb^-1 published (D0 hep-ex/0803.1514)'
! information about how data is stored in file:	 
	Mhmin_D0_ppH_gamgam_pub=70d0
	Mhmax_D0_ppH_gamgam_pub=150d0
	sep_D0_ppH_gamgam_pub=10d0

	Mhmin=Mhmin_D0_ppH_gamgam_pub
	Mhmax=Mhmax_D0_ppH_gamgam_pub
        sep=sep_D0_ppH_gamgam_pub
	ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 9. 

	open(10,file=INSTDIR//"/TEVtables/D0_ppH_gamgam_pub-abs-limit_exp.txt")
	open(11,file=INSTDIR//"/TEVtables/D0_ppH_gamgam_pub-abs-limit_obs.txt")

        do i=1,ntot
	  read(10,*)Mh1dummy,D0_ppH_gamgam_pub_pred(i)
          read(11,*)Mh1dummy,D0_ppH_gamgam_pub_obs(i)
        enddo

        close(10)	 	
	close(11)

	Mhmin_overall=min(Mhmin_overall,Mhmin)
        Mhmax_overall=max(Mhmax_overall,Mhmax)


!b) D0 with 2.7 fb^-1
	write(*,*)'  D0 with 2.7 fb^-1 published (D0 hep-ex/0901.1887)'
! information about how data is stored in file:	 
	Mhmin_D0_ppH_gamgam_pub_2=100d0
	Mhmax_D0_ppH_gamgam_pub_2=150d0
	sep_D0_ppH_gamgam_pub_2=5d0

	Mhmin=Mhmin_D0_ppH_gamgam_pub_2
	Mhmax=Mhmax_D0_ppH_gamgam_pub_2
        sep=sep_D0_ppH_gamgam_pub_2
	ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 11. 

	open(10,file=INSTDIR//"/TEVtables/D0_ppH_gamgam_pub_2_exp.txt")
	open(11,file=INSTDIR//"/TEVtables/D0_ppH_gamgam_pub_2_obs.txt")

        do i=1,ntot
	  read(10,*)Mh1dummy,D0_ppH_gamgam_pub_2_pred(i)
          read(11,*)Mh1dummy,D0_ppH_gamgam_pub_2_obs(i)
        enddo

        close(10)	 	
	close(11)

	Mhmin_overall=min(Mhmin_overall,Mhmin)
        Mhmax_overall=max(Mhmax_overall,Mhmax)





!10) channel: p p -> b H, H-> tau tau
	write(*,*)'reading tables for pp -> bH, H -> tau tau'

!a) D0 with 0.328 fb^-1 published (hep-ex/0811.0024), absolute sigma x BR table in pikobarn
	write(*,*)'  D0 absolute sigma x BR table in pikobarn for 0.328 fb^-1, published (hep-ex/0811.0024)'
! information about how data is stored in file:	 
	Mhmin_D0_bH_tautau_pub=90d0
	Mhmax_D0_bH_tautau_pub=150d0
	sep_D0_bH_tautau_pub=10d0

	Mhmin=Mhmin_D0_bH_tautau_pub
	Mhmax=Mhmax_D0_bH_tautau_pub
        sep=sep_D0_bH_tautau_pub
	ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 7. 

	open(10,file=INSTDIR//"/TEVtables/D0_bH_tautau_pub-abs-limit_exp.txt")
	open(11,file=INSTDIR//"/TEVtables/D0_bH_tautau_pub-abs-limit_obs.txt")

        do i=1,ntot
	  read(10,*)Mh1dummy,D0_bH_tautau_pub_pred(i)
          read(11,*)Mh1dummy,D0_bH_tautau_pub_obs(i)
        enddo

        close(10)	 	
	close(11)

	Mhmin_overall=min(Mhmin_overall,Mhmin)
        Mhmax_overall=max(Mhmax_overall,Mhmax)


!b) D0 with 2.7 fb^-1 (D0 note 5985), absolute sigma x BR table in pikobarn
	write(*,*)'  D0 absolute sigma x BR table in pikobarn for 2.7 fb^-1, published (D0 note 5985)'
! information about how data is stored in file:	 
	Mhmin_D0_bH_tautau=90d0
	Mhmax_D0_bH_tautau=320d0
	sep_D0_bH_tautau=10d0

	Mhmin=Mhmin_D0_bH_tautau
	Mhmax=Mhmax_D0_bH_tautau
        sep=sep_D0_bH_tautau
	ntot=nint((Mhmax-Mhmin)/sep)+1
! ntot is expected to be 24. 

	open(10,file=INSTDIR//"/TEVtables/D0_bH_tautau-abs-limit_exp.txt")
	open(11,file=INSTDIR//"/TEVtables/D0_bH_tautau-abs-limit_obs.txt")

        do i=1,ntot
	  read(10,*)Mh1dummy,D0_bH_tautau_pred(i)
          read(11,*)Mh1dummy,D0_bH_tautau_obs(i)
        enddo

        close(10)	 	
	close(11)

	Mhmin_overall=min(Mhmin_overall,Mhmin)
        Mhmax_overall=max(Mhmax_overall,Mhmax)


	write(*,*)'finished reading Tevatron tables'

	write(*,123) Mhmin_overall,Mhmax_overall
 123	FORMAT(' Overall mass range of Tevatron tables: [',G10.4,',',G10.4,'] GeV.')
        write(*,*)

	end





! >> include fit functions for SM cross section ratios
#include "cs-ratios/Tevatron.R_H_bb.h"
#include "cs-ratios/Tevatron.R_H_gg.h"
#include "cs-ratios/Tevatron.R_HWp_csb.h"
#include "cs-ratios/Tevatron.R_HWm_dub.h"
#include "cs-ratios/Tevatron.R_HWm_scb.h"
#include "cs-ratios/Tevatron.R_HWp_udb.h"
#include "cs-ratios/Tevatron.R_HZ_bbb.h"
#include "cs-ratios/Tevatron.R_HZ_ccb.h"
#include "cs-ratios/Tevatron.R_HZ_ddb.h"
#include "cs-ratios/Tevatron.R_HZ_ssb.h"
#include "cs-ratios/Tevatron.R_HZ_uub.h"
! >> include fit functions for SM cross sections and branching ratios
#include "theory_BRfunctions.F"
#include "theory_XS_SM_functions.f"
#include "Tevatron.sigma_bg_Hb.h"
#include "Tevatron.sigma_bg_Hb.ptmin15.etamax2.h"
#include "Tevatron.sigma_bg_Hb.ptmin20.etamax2.5.h"
#include "Tevatron.sigma_bg_Hb.ptmin15.etamax2.5.h"
! <<


! TEV_NSCHN: # of single Higgs search channels
!   current channels:
!CDF0908_3534	: p p-bar -> ZH -> l l b b-bar, CDF with 2.7 fb^-1
!			(hep-ex/0908.3534, CDF_ZH_llbb_pub_2_exp.txt  CDF_ZH_llbb_pub_2_obs.txt)
!CDF9889        : p p-bar -> ZH -> l l b b-bar, CDF with 4.1 fb^-1
!			(CDF note 9889, CDF_ZH_llbb_4_exp.txt  CDF_ZH_llbb_4_obs.txt)
!CDF9891	: p p-bar -> H via HW -> b b-bar + miss. E_T,	
!			     HZ -> b b-bar + miss. E_T, SM combined, CDF with 3.6 fb^-1
!			(CDF note 9891, CDF_VH_bbmET_4_exp.txt  CDF_VH_bbmET_4_obs.txt)
!D05876         : p p-bar -> ZH -> l l b b-bar, D0 with 4.2 fb^-1
!			(D0 note 5876, D0_ZH_llbb_3_exp.txt  D0_ZH_llbb_3_obs.txt)
!D05586		: p p-bar -> H via HW -> b b-bar + miss. E_T,
!                            HZ -> b b-bar + miss. E_T, SM combined, D0 with 2.1 fb^-1
!			(D0 note 5586, D0_VH_bbmET_exp.txt  D0_VH_bbmET_obs.txt)
!CDF8961_D05536	: p p-bar -> H via HW -> l nu b b-bar, W W+ W- -> l+- nu l+- nu ..,
!			     HZ -> nu nu b b-bar, l l b b-bar,
!			     H -> W W -> l+- nu l-+ nu, SM combined 
!			(CDF note 8961, D0 note 5536, CDF_D0_SM_exp.txt  CDF_D0_SM_obs.txt)
!CDF9284	: p p-bar -> b H, H -> b b-bar, CDF with 1.9fb^-1
!			(CDF note 9284, CDF_bH_bb_tb_limit_exp.txt  CDF_bH_bb_tb_limit_obs.txt)
!D05972 	: p p-bar -> WH -> l nu b b-bar, D0 with 5.0 fb^-1
!			(D0 note 5828, D0_WH_lnubb_3_exp.txt  D0_WH_lnubb_3_obs.txt)
!CDFwhlnubb_090814	: p p-bar -> WH -> l nu b b-bar, CDF with 4.3 fb^-1
!			(CDF http://www-cdf.fnal.gov/physics/new/hdg/results/whlnubb_090814/WH4.3fb.html,
!			 CDF_WH_lnubb_5_exp.txt  CDF_WH_lnubb_5_obs.txt)
!D05726		: p p-bar -> b H, H -> b b-bar, D0
!			(D0 note 5726, D0_bH_bb_tb_limit_2_exp.txt  D0_bH_bb_tb_limit_2_obs.txt)
!D00805_3556	: p p-bar -> b H, H -> b b-bar, D0
!                       (hep-ex/0805.3556, D0_bH_bb-abs-limit_exp.txt  D0_bH_bb-abs-limit_obs.txt)
!	16: p p-bar -> b H, H -> b b-bar, CDF  with 1.9fb^-1 absolute limits
!                       (CDF note 9284, CDF_bH_bb-abs-limit_exp.txt   CDF_bH_bb-abs-limit_obs.txt)
!D00805_2491	: p p-bar -> H -> tau tau, D0 with 1fb^-1 absolute limits
!			(D0 hep-ex/0805.2491, D0_ppH_tautau_pub_exp.txt   D0_ppH_tautau_pub_obs.txt)
!CDF0906_1014	: p p-bar -> H -> tau tau, CDF with 1.8 fb^-1 absolute limits
!			(CDF note 9071, CDF_ppH_tautau-abs-limit_pub_exp.txt   CDF_ppH_tautau-abs-limit_pub_obs.txt)
!D05873		:  p p-bar -> WH, H -> WW, D0 with 3.6 fb^-1
!			(D0 note 5873, D0_WH_WWW_2_exp.txt  D0_WH_WWW_2_obs.txt)
!CDF7307v3	:  p p-bar -> WH, H -> WW, CDF with 2.7 fb^-1
!			(CDF note 7307, CDF_WH_WWW_2_exp.txt  CDF_WH_WWW_2_obs.txt)
!D05858		: p p-bar -> H -> gamma gamma via single H,
!			     HZ, HW, VBF, D0 with 4.2 fb^-1
!			(D0 note 5858, D0_ppH_gamgam_4_exp.txt  D0_ppH_gamgam_4_obs.txt)
!temporary removed !D00803_1514	: p p-bar -> H -> gamma gamma, D0 with 1.1 fb^-1 absolute limits
!temporary removed !			(D0 hep-ex/0803.1514, D0_ppH_gamgam_pub-abs-limit_exp.txt  D0_ppH_gamgam_pub-abs-limits_obs.txt)
!remark: D00803_1514 needs fermiophobic-ness check !
!CDF9248	: p p-bar -> H -> tau tau via HW with W -> j j, 
!			     HZ with Z -> j j,
!			     single H, VBF, SM combined 
!			(CDF note 9248, CDF_SM_tautau_exp.txt  CDF_SM_tautau_obs.txt)
!CDF9290_D05645	: p p-bar -> H via HW -> l nu b b-bar, W W+ W- -> l+- nu l+- nu .., 
!			     HZ -> nu nu b b-bar, l l b b-bar,
!			     H -> W W -> l+- nu l-+ nu, 
!			     H -> tau tau via HW with W -> j j,
!                                             HZ with Z -> j j,
!                                             single H, VBF,
!			     H -> gamma gamma via ...
!			     SM combined 
!			(CDF note 9290, D0 note 5645, CDF_D0_SM_2_exp.txt  CDF_D0_SM_2_obs.txt)
!D05757		: p p-bar -> H -> W W -> l l, D0 with 3.0 fb^-1
!			(D0 note 5757, D0_ppH_WW_ll_2_exp.txt  D0_ppH_WW_ll_2_obs.txt)
!CDF0809_3930	: p p-bar -> H -> W W, CDF with 3.0 fb^-1 (PRL draft)
!                       (hep-ex/0809.3930, CDF_ppH_WW_4_exp.txt  CDF_ppH_WW_4_obs.txt)
!CDF9465_D05754	: p p-bar -> H via HW -> W W+ W- -> l+- nu l+- nu .., 
!			            H -> W W -> l+- nu l-+ nu (D0); 
!		and p p-bar -> H via HW -> W W+ W- -> l+- nu l+- nu ..,
!			       H+0jet(single H),
!			       H+1jet(single H/WH/ZH/H via VBF)
!			       H+2jet(single H/WH/ZH/H via VBF) (CDF)
!			     SM combined 
!			(hep-ex/0808.0534, CDF note 9465, D0 note 5754, CDF_D0_SM_3_exp.txt  CDF_D0_SM_3_obs.txt)
!CDF9887  : p p-bar -> H + X -> W W + X via single H,
!				  HZ, HW, VBF, CDF with 4.8 fb^-1 
!                       (CDF note 9887, CDF_ppH_WW_6_exp.txt  CDF_ppH_WW_6_obs.txt)
!D05871		: p p-bar -> H + X -> W W + X-> l l + X via single H,
!				  HZ, HW, VBF, D0 with 3.0-4.2 fb^-1
!			(D0 note 5871, D0_ppH_WW_ll_3_exp.txt  D0_ppH_WW_ll_3_obs.txt)
!D00808_1970 	: p p-bar -> WH -> l nu b b-bar, D0 with 1.1 fb^-1
!			(hep-ex/0808.1970, D0_WH_lnubb_pub_exp.txt  D0_WH_lnubb_pub_obs.txt)
!CDF0906_5613 	: p p-bar -> WH -> l nu b b-bar, CDF with 2.7 fb^-1
!			(hep-ex/0906.5613, CDF_WH_lnubb_pub_2_exp.txt  CDF_WH_lnubb_pub_2_obs.txt)
!D00811_0024	: p p-bar -> b H, H -> tau tau, D0 with 0.328 fb^-1
!                       (hep-ex/0811.0024, D0_bH_tautau_pub-abs-limit_exp.txt  D0_bH_tautau_pub-abs-limit_obs.txt)
!D05985		: p p-bar -> b H, H -> tau tau, D0 with 2.7 fb^-1
!                       (D0 not 5985, D0_bH_tautau-abs-limit_exp.txt  D0_bH_tautau-abs-limit_obs.txt)
!D05740         : p p-bar -> H -> tau tau, D0 with 2.2 fb^-1 absolute limits
!                       (D0 note 5740, D0_ppH_tautau-abs-limit_exp.txt D0_ppH_tautau-abs-limit_obs.txt)
!CDF9888_D05980	: p p-bar -> H -> tau tau, CDF & D0 combined with 1.8 & 2.2 fb^-1 absolute limits
!                       (CDF note 9888, D0 note 5980, CDF_D0_ppH_tautau-abs-limit_exp.txt CDF_D0_ppH_tautau-abs-limit_obs.txt)
!D00808_1266	: p p-bar -> H via HW -> b b-bar + miss. E_T,
!                            HZ -> b b-bar + miss. E_T, SM combined, D0 with 0.93 fb^-1
!			(hep-ex/0808.1266, D0_VH_bbmET_pub_exp.txt  D0_VH_bbmET_pub_obs.txt)
!CDF0802_0432	: p p-bar -> H via HW -> b b-bar + miss. E_T,	
!			     HZ -> b b-bar + miss. E_T, SM combined, CDF with 1 fb^-1
!			(hep-ex/0802.0432, CDF_VH_bbmET_pub_exp.txt  CDF_VH_bbmET_pub_obs.txt)
!D00901_1887	: p p-bar -> H -> gamma gamma, D0 with 2.7 fb^-1
!			(hep-ex/0901_1887, D0_ppH_gamgam_pub_2_exp.txt  D0_ppH_gamgam_pub_2_obs.txt)
!CDF9674	: p p-bar -> H via HW -> l nu b b-bar, W W+ W- -> l+- nu l+- nu .., 
!			     HZ -> nu nu b b-bar, l l b b-bar,
!			     H -> W W -> l+- nu l-+ nu, 
!			     H -> tau tau via HW with W -> j j,
!                                             HZ with Z -> j j,
!                                             single H, VBF,
!			     SM combined 
!			(CDF note 9674, CDF_SM_exp.txt  CDF_SM_obs.txt)
!D00712_0598	: p p-bar -> H via HW -> l nu b b-bar, W W+ W- -> l+- nu l+- nu .., 
!			     HZ -> nu nu b b-bar, l l b b-bar,
!			     H -> W W -> l+- nu l-+ nu, 
!			     SM combined 
!			(hep-ex/0712.0598, D0_SM_pub_exp.txt  D0_SM_pub_obs.txt)
!CDF9713_D05889	: p p-bar -> HW -> b b l nu, b b-bar tau nu_tau, W+ W- W -> l+- nu l+- nu
!                         -> HV -> j j tau tau/tau tau j j (for m_H <= 145 GeV), b b + miss. E_T (V=W,Z)  
!                         -> HZ -> b b l l, 
!                         -> H via single H/HW/HZ/VBF + 0/1/2 jets, H -> W+ W- -> l+ nu l- nu,
!                         -> Htt -> b b t t (for m_H <= 155 GeV)
!                         -> H via single H/WH/ZH/VBF, H -> tau tau, gamma gamma
!                         -> HV -> b b j j (V=W,Z)
!			     SM combined 
!			(hep-ex/0808.0534, CDF note 9465, D0 note 5754, CDF_D0_SM_4_exp.txt  CDF_D0_SM_4_obs.txt)
! 		analysis currently only used for m_H > 150 GeV:
!			-> the HV -> j j tau tau/tau tau j j channel doesn't contribute
!			-> the Htt channel is ignored
!			-> the HW -> b b l nu, b b-bar tau nu_tau channels don't contribute
!                       -> the HZ -> b b l l channel doesn't contribute
!                       -> the H via single H/WH/ZH/VBF, H -> tau tau, gamma gamma channels don't contribute
!                       -> the HV -> b b j j (V=W,Z) channel doesn't contribute
!CDF9897	: p p-bar -> HW -> l nu b b-bar,
!			     HZ -> nu nu b b-bar, l l b b-bar,
!			     HV -> b b-bar j j,
!			     H via single H/HW/HZ/VBF, H -> W W
!			     H via single H/HW/HZ/VBF, H -> tau tau 
!			     SM combined 
!			(CDF note 9897, CDF_SM_2_exp.txt  CDF_SM_2_obs.txt)


!
! TEV_NDCHN: # of double Higgs search channels
!   currently: none

!***********************************************************
	subroutine check_TEV_bound(NHIGGS, whichexpt, Mh, Gamma_h,
     &		R_sigma_H,R_sigma_Hb,R_sigma_HW,R_sigma_HZ,R_sigma_H_VBF,
     &		BRhjbbbar,BRhjtautau,BRhjWW,BRhjgamgam,
     &		allowed,result_flag,mlratio,obsratio,most_sens_predratio,
     &		num_H_combined,final_h_comb,model_R_sig_BR)
!***********************************************************
! version for input option 2 (see manual)
! Input: 
!	Mh(NHMAX) : the NHMAX neutral Higgs masses
!    ratio of hadronic cross section for the processes:
!	R_sigma_H(NHMAX) 	: p p-bar -> H + X
!	R_sigma_Hb(NHMAX) 	: p p-bar -> H b/b-bar + X
!	R_sigma_HW(NHMAX) 	: p p-bar -> H W+/- + X
!	R_sigma_HZ(NHMAX) 	: p p-bar -> H Z + X
!	R_sigma_H_VBF(NHMAX) 	: p p-bar -> H + 2jet + X via VBF
!    model branching ratios for the decays:
!	BRhjbbbar(NHMAX) : H_j -> b b-bar
!	BRhjtautau(NHMAX): H_j -> tau+ tau-
!	BRhjWW(NHMAX)    : H_j -> W+ W-
!	BRhjgamgam(NHMAX): H_j -> gamma gamma
! Output:
!	allowed : logical variable, TRUE if scenario is allowed at
!		95% confidence level (This is determined, using the procedure
!	 	described by Philipp Bechtle.), FALSE otherwise.
!	result_flag
!		1 if point is allowed, 
!		0 if point is not allowed
!	       -1 otherwise (theoretically inaccessible, etc.)
!		TEV_NSCHN: # of single Higgs search channels
!		TEV_NDCHN: # of double Higgs search channels
!	mlratio : gives information about which table was used for the 
!		experimental comparison.
!	obsratio,most_sens_predratio,num_H_combined,final_h_comb,model_R_sig_BR
!--------------------------------
	implicit none 
	logical allowed
	integer NHIGGS,my_maxloc
	character*5 whichexpt
	
! This form only if ever double Higgs channels are included in the
! SM combined analysis
!	logical table_applicable(TEV_NSCHN*NHMAX+TEV_NDCHN*NHMAX**2)
! This form assumes that there only single Higgs channels are included in the
! SM combined analysis
	logical table_applicable(TEV_NSCHN*NHMAX)

	integer tabsize
	double precision Mhmin,Mhmax,sep

	integer result_flag,num_H_combined
	integer i,j,k,n,m,ndat,ios,mlratio,nhelp,chan,output_chan
	integer mlratiotemp(1)	
	double precision R_sig_BR(TEV_NSCHN*NHMAX+TEV_NDCHN*NHMAX**2)
        double precision R_sig_BR_sing(TEV_NSCHN*NHMAX+TEV_NDCHN*NHMAX**2)
        double precision predratio(TEV_NSCHN*NHMAX+TEV_NDCHN*NHMAX**2)
	double precision most_sens_predratio
	double precision obsratio
! ratio of BRs
	double precision R_BRhjbbbar(NHMAX),R_BRhjtautau(NHMAX)
	double precision R_BRhjWW(NHMAX),R_BRhjgamgam(NHMAX)
! BRs
	double precision BRhjbbbar(NHMAX),BRhjtautau(NHMAX)
	double precision BRhjWW(NHMAX),BRhjgamgam(NHMAX)
	double precision Mh(NHMAX),Gamma_h(NHMAX)
	double precision sum_sig_BR_SM, sum_sig_BR_SM_tautau, sum_sig_BR_SM_2
        double precision sum_sig_BR_SM_3,sum_sig_BR_SM_4,sum_sig_BR_CDF_SM,sum_sig_BR_D0_SM
	double precision sum_sig_BR_SM_5
! ratios of hadronic cross sections (new model/SM)
	double precision R_sigma_HW(NHMAX),R_sigma_HZ(NHMAX),R_sigma_H(NHMAX),
     &			 R_sigma_Hb(NHMAX),R_sigma_H_VBF(NHMAX)
! fit functions for SM cross sections and branching ratios
	double precision SMBR_HWW
	double precision SMBR_Hbb
	double precision SMBR_Htautau
        double precision SMBR_Hgamgam
	double precision SMCS_tev_pp_qq_HW
	double precision SMCS_tev_pp_qq_HZ
	double precision SMCS_tev_pp_gg_H
	double precision SMCS_tev_pp_bb_H
        double precision SMCS_tev_pp_bg_Hb
        double precision SMCS_tev_pp_bg_Hb_c1
        double precision SMCS_tev_pp_bg_Hb_c2
        double precision SMCS_tev_pp_bg_Hb_c3
        double precision SMCS_tev_pp_vbf_H
! Decay BRs of W+/- and Z
	double precision BR_Zll,BR_Znunu,BR_Zjj,BR_Wlnu,BR_Wjj,BR_Wtaunu
  
	integer findi,findj,TEV_reaction_number
	
	logical higgs_combination(TEV_NSCHN*NHMAX,NHMAX),final_h_comb(NHMAX)
	integer sim_mass_sum(TEV_NSCHN*NHMAX),similar_mass(TEV_NSCHN*NHMAX)
        double precision Mh_avg(TEV_NSCHN*NHMAX),delta_Mh
	double precision log2dble,model_R_sig_BR
	double precision Mhmin_SMcalc,Mhmax_SMcalc

#include "tables-tev_cb.h"

	Mhmin_SMcalc=65d0
	Mhmax_SMcalc=335d0

! decay BRs of Z and W bosons
!BR_Zll = sum(l=e,mu) BR(Z ->l+ l-) source: Yao et al. J Phys G 33 (2006) 1
	BR_Zll=(3.363d0+3.366d0)/100d0
!BR_Znunu = sum(l=e,mu,tau) BR(Z ->nu_l nu_l-bar)
	BR_Znunu=20d0/100d0
!BR_Wlnu = sum(l=e,mu) BR(W+ ->l + nu_l) = BR(W- ->l- + nu_l-bar)
	BR_Wlnu=(10.75d0+10.57d0)/100d0
!BR_Wtaunu =  BR(W+ ->tau+ nu_tau) = BR(W- -> tau- nu_tau-bar)
	BR_Wtaunu=(11.25d0)/100d0
! BRs W,Z->hadrons from Yao et al. J Phys G 33 (2006) 1
	BR_Wjj = 67.60d0/100d0
	BR_Zjj = 69.91d0/100d0

	delta_Mh = delta_Mh_TEV

	do m=0,TEV_NSCHN-1,1
	do n=1,NHIGGS,1
	 R_sig_BR_sing(n+NHIGGS*m)=0d0
	 R_sig_BR(n+NHIGGS*m)=0d0
	 sim_mass_sum(n+NHIGGS*m)=0
	 do i=1,NHIGGS
          higgs_combination(n+NHIGGS*m,i)=.FALSE.
	 enddo
	enddo
	enddo
	do m=0,TEV_NDCHN-1,1
	do n=1,NHIGGS**2
	 R_sig_BR_sing(n+m*NHIGGS**2+TEV_NSCHN*NHIGGS)=0d0
	 R_sig_BR(n+m*NHIGGS**2+TEV_NSCHN*NHIGGS)=0d0
	enddo				
	enddo
	

! R_sig_BR(n) will be compared to the S95 values in the Tevatron tables


	do n=1,NHIGGS
	    if( (Mh(n) .ge. Mhmin_SMcalc) .and. (Mh(n) .le. Mhmax_SMcalc) ) then 
!tag1
	      if(SMBR_Hbb(Mh(n)) .ne. 0d0) then
		R_BRhjbbbar(n)=BRhjbbbar(n)/SMBR_Hbb(Mh(n))
	      else
		R_BRhjbbbar(n)=0d0
	      endif

	      if(SMBR_HWW(Mh(n)) .ne. 0d0) then
		R_BRhjWW(n)=BRhjWW(n)/SMBR_HWW(Mh(n))
	      else
                R_BRhjWW(n)=0d0
              endif
	      if(SMBR_Htautau(Mh(n)) .ne. 0d0) then
		R_BRhjtautau(n)=BRhjtautau(n)/SMBR_Htautau(Mh(n))
	      else
		R_BRhjtautau(n)=0d0
	      endif
	      if(SMBR_Hgamgam(Mh(n)) .ne. 0d0) then
	    	R_BRhjgamgam(n)=BRhjgamgam(n)/SMBR_Hgamgam(Mh(n))
	      else
                R_BRhjgamgam(n)=0d0
              endif



	      do m=0,TEV_NSCHN-1,1
	       call check_SM_likeness(table_applicable(n+NHIGGS*m), m+1,Mh(n),
     &	 	R_sigma_HW(n),R_sigma_HZ(n),R_sigma_H(n),R_sigma_H_VBF(n), 
     &		R_BRhjbbbar(n),R_BRhjtautau(n),R_BRhjWW(n),R_BRhjgamgam(n))

	      enddo
	    endif 
!tag1
	enddo



! a0) calc single sigma x BR values (no combinations of near mass Higgs bosons)
!    single Higgs channels
	do i=1,NHIGGS
	   if( (Mh(i) .ge. Mhmin_SMcalc) .and. (Mh(i) .le. Mhmax_SMcalc) ) then 
! D05972, CDFwhlnubb_090814, D00808_1970, CDF0906_5613: p p-bar -> WH, H -> b b-bar, W -> l nu
	   	  R_sig_BR_sing(i+NHIGGS*(D05972 -1))=
     &			R_sigma_HW(i)*SMCS_tev_pp_qq_HW(Mh(i))*BRhjbbbar(i)
	   	  R_sig_BR_sing(i+NHIGGS*(CDFwhlnubb_090814 -1))=R_sig_BR_sing(i+NHIGGS*(D05972 -1))
	   	  R_sig_BR_sing(i+NHIGGS*(D00808_1970 -1))=R_sig_BR_sing(i+NHIGGS*(D05972 -1))
	   	  R_sig_BR_sing(i+NHIGGS*(CDF0906_5613 -1))=R_sig_BR_sing(i+NHIGGS*(D05972 -1))

! CDF0908_3534, CDF9889, D05876 : p p-bar -> ZH -> l l b b-bar
	          R_sig_BR_sing(i+NHIGGS*(CDF0908_3534 -1))=
     &                  R_sigma_HZ(i)*SMCS_tev_pp_qq_HZ(Mh(i))*BRhjbbbar(i)
	   	  R_sig_BR_sing(i+NHIGGS*(CDF9889 -1))= R_sig_BR_sing(i+NHIGGS*(CDF0908_3534 -1))
           	  R_sig_BR_sing(i+NHIGGS*(D05876 -1)) = R_sig_BR_sing(i+NHIGGS*(CDF0908_3534 -1))

! D05757, CDF0809_3930 : p p-bar -> H -> W W
           	  R_sig_BR_sing(i+NHIGGS*(D05757 -1))=
     &			R_sigma_H(i)*( SMCS_tev_pp_gg_H(Mh(i)) + SMCS_tev_pp_bb_H(Mh(i)) ) 
     &			 *BRhjWW(i)
	   	  R_sig_BR_sing(i+NHIGGS*(CDF0809_3930 -1))=R_sig_BR_sing(i+NHIGGS*(D05757 -1))

! CDF9887, D05871 : p p-bar -> H + X -> W W + X
           	  R_sig_BR_sing(i+NHIGGS*(CDF9887 -1))=
     &			log2dble(table_applicable(i+NHIGGS*(CDF9887 -1)))*( 
     &			  R_sigma_H(i)*( SMCS_tev_pp_gg_H(Mh(i)) + SMCS_tev_pp_bb_H(Mh(i)) ) 
     &			 +R_sigma_HW(i)*SMCS_tev_pp_qq_HW(Mh(i))
     &                   +R_sigma_HZ(i)*SMCS_tev_pp_qq_HZ(Mh(i))
     &                   +R_sigma_H_VBF(i)*SMCS_tev_pp_vbf_H(Mh(i))
     &			 ) *BRhjWW(i)
		   R_sig_BR_sing(i+NHIGGS*(D05871 -1))= R_sig_BR_sing(i+NHIGGS*(CDF9887 -1))

! CDF9284, D05726 : p p-bar -> b H, H -> b b-bar ([sigma x BR/sigma limit)
	   	  R_sig_BR_sing(i+NHIGGS*(CDF9284 -1))=
     &				R_sigma_Hb(i)*SMCS_tev_pp_bg_Hb(Mh(i))*BRhjbbbar(i) 
	   	  R_sig_BR_sing(i+NHIGGS*(D05726 -1))=R_sig_BR_sing(i+NHIGGS*(CDF9284 -1))

! D00805_3556 : p p-bar -> b H, H -> b b-bar (absolute sigma x BR limit)
	   	  R_sig_BR_sing(i+NHIGGS*(D00805_3556 -1))=
     &                                R_sigma_Hb(i)*SMCS_tev_pp_bg_Hb_c2(Mh(i))*BRhjbbbar(i) 

! 16 : p p-bar -> b H, H -> b b-bar (absolute sigma x BR limit)
!      the SM c.s. quoted in the CDF paper is a factor 0.7 lower than my c.s.
!      The CDF note 9284 result is also implemented as tab 10.
!      Table 16 was for checking whether the the absolute limit is consistent.
!      It is! So it should be switched off again.
!	   	  R_sig_BR_sing(i+NHIGGS*15)=
!     &               R_sigma_Hb(i)*0.7d0*SMCS_tev_pp_bg_Hb_c1(Mh(i))*BRhjbbbar(i) 

! D00805_2491, CDF0906_1014, D05740, CDF9888_D05980 : p p-bar -> H -> tau tau (absolute sigma x BR limit)
                  R_sig_BR_sing(i+NHIGGS*(D00805_2491 -1))=
     &             R_sigma_H(i)*( SMCS_tev_pp_gg_H(Mh(i)) + SMCS_tev_pp_bb_H(Mh(i)) ) 
     &                  *BRhjtautau(i)                          
                  R_sig_BR_sing(i+NHIGGS*(CDF0906_1014 -1))=R_sig_BR_sing(i+NHIGGS*(D00805_2491 -1))
                  R_sig_BR_sing(i+NHIGGS*(D05740 -1))=R_sig_BR_sing(i+NHIGGS*(D00805_2491 -1))
                  R_sig_BR_sing(i+NHIGGS*(CDF9888_D05980 -1))=R_sig_BR_sing(i+NHIGGS*(D00805_2491 -1))

! CDF7307v3, D05873 : p p-bar -> WH, H -> WW 
	   	  R_sig_BR_sing(i+NHIGGS*(CDF7307v3 -1))= 
     &			R_sigma_HW(i)*SMCS_tev_pp_qq_HW(Mh(i))*BRhjWW(i)
		  R_sig_BR_sing(i+NHIGGS*(D05873 -1))=R_sig_BR_sing(i+NHIGGS*(CDF7307v3 -1))

! CDF9891, D05586, D00808_1266, CDF0802_0432 : p p-bar -> WH/ZH, H -> b b-bar
	   	  R_sig_BR_sing(i+NHIGGS*(CDF9891 -1))=
     &			log2dble(table_applicable(i+NHIGGS*(CDF9891 -1)))*( 
     &			 R_sigma_HW(i)*SMCS_tev_pp_qq_HW(Mh(i))*BR_Wlnu
     &			+R_sigma_HZ(i)*SMCS_tev_pp_qq_HZ(Mh(i))*BR_Znunu
     &			 )*BRhjbbbar(i)
           	  R_sig_BR_sing(i+NHIGGS*(D05586 -1)) =R_sig_BR_sing(i+NHIGGS*(CDF9891 -1))
           	  R_sig_BR_sing(i+NHIGGS*(D00808_1266 -1)) =R_sig_BR_sing(i+NHIGGS*(CDF9891 -1))
           	  R_sig_BR_sing(i+NHIGGS*(CDF0802_0432 -1)) =R_sig_BR_sing(i+NHIGGS*(CDF9891 -1))

! D05858, D00901_1887 : p p-bar -> H -> gamma gamma
           	  R_sig_BR_sing(i+NHIGGS*(D05858 -1))=
     &			 log2dble(table_applicable(i+NHIGGS*(D05858 -1)))*( 
     &			 R_sigma_H(i)*( SMCS_tev_pp_gg_H(Mh(i))  + SMCS_tev_pp_bb_H(Mh(i)) )
     &		  	+R_sigma_HW(i)*SMCS_tev_pp_qq_HW(Mh(i)) 
     &	 	  	+R_sigma_HZ(i)*SMCS_tev_pp_qq_HZ(Mh(i)) 
     &		  	+R_sigma_H_VBF(i)*SMCS_tev_pp_vbf_H(Mh(i)) ) 
     &		 	*BRhjgamgam(i)
                  R_sig_BR_sing(i+NHIGGS*(D00901_1887 -1)) =R_sig_BR_sing(i+NHIGGS*(D05858 -1))


!temporary removed !D00803_1514 : p p-bar -> H -> gamma gamma (absolute sigma x BR limit)
!temporary removed            	  R_sig_BR_sing(i+NHIGGS*(D00803_1514 -1))=
!temporary removed      &			 ( R_sigma_H(i)*( SMCS_tev_pp_gg_H(Mh(i))  + SMCS_tev_pp_bb_H(Mh(i)) )
!temporary removed      &			  +R_sigma_HW(i)*SMCS_tev_pp_qq_HW(Mh(i)) 
!temporary removed      &	 	  	+R_sigma_HZ(i)*SMCS_tev_pp_qq_HZ(Mh(i)) 
!temporary removed      &		  	+R_sigma_H_VBF(i)*SMCS_tev_pp_vbf_H(Mh(i)) ) 
!temporary removed      &		 	*BRhjgamgam(i)

! D00811_0024, D05985 : p p-bar -> b H, H -> tau tau (absolute sigma x BR limit)
	   	  R_sig_BR_sing(i+NHIGGS*(D00811_0024 -1))=
     &                                R_sigma_Hb(i)*SMCS_tev_pp_bg_Hb_c3(Mh(i))*BRhjtautau(i) 
                  R_sig_BR_sing(i+NHIGGS*(D05985 -1)) =R_sig_BR_sing(i+NHIGGS*(D00811_0024 -1))



! SM combined analyses
! CDF8961_D05536, D00712_0598: : p p-bar -> (WH/ZH, H-> b b-bar)/(H,WH, H -> W W)

 	   	  R_sig_BR_sing(i+NHIGGS*(CDF8961_D05536 -1))=
     &			 log2dble(table_applicable(i+NHIGGS*(CDF8961_D05536 -1)))*( 
     &			 R_sigma_HW(i)*BRhjbbbar(i)*SMCS_tev_pp_qq_HW(Mh(i))*BR_Wlnu 
     &			+R_sigma_HW(i)*BRhjWW(i)*SMCS_tev_pp_qq_HW(Mh(i))*BR_Wlnu**2 
     &			+R_sigma_HZ(i)*BRhjbbbar(i)*SMCS_tev_pp_qq_HZ(Mh(i))*(BR_Znunu+BR_Zll) 
     &			+R_sigma_H(i)*BRhjWW(i)*( SMCS_tev_pp_gg_H(Mh(i)) + SMCS_tev_pp_bb_H(Mh(i)) ) 
     &			  *BR_Wlnu**2 )
                  R_sig_BR_sing(i+NHIGGS*(D00712_0598 -1)) =R_sig_BR_sing(i+NHIGGS*(CDF8961_D05536 -1))


! CDF9248 : p p-bar -> (single H, H via VBF ,WH/ZH), H -> tau tau

 	   	  R_sig_BR_sing(i+NHIGGS*(CDF9248 -1))=
     &			 log2dble(table_applicable(i+NHIGGS*(CDF9248 -1)))*( 
     &			 R_sigma_HW(i)*BR_Wjj*SMCS_tev_pp_qq_HW(Mh(i)) 
     &			+R_sigma_HZ(i)*BR_Zjj*SMCS_tev_pp_qq_HZ(Mh(i)) 
     &			+R_sigma_H(i)*( SMCS_tev_pp_gg_H(Mh(i)) + SMCS_tev_pp_bb_H(Mh(i)) ) 
     &			+R_sigma_H_VBF(i)*SMCS_tev_pp_vbf_H(Mh(i)) ) 
     &			*BRhjtautau(i)

! CDF9290_D05645 : p p-bar -> (WH/ZH, H-> b b-bar)/(H/WH, H -> W W)
!                  /(H/WH/ZH/H via VBF, H->tau tau)
!		   /(H/WH/ZH/H via VBF,H -> gamma gamma)

 	   	  R_sig_BR_sing(i+NHIGGS*(CDF9290_D05645 -1))=
     &			 log2dble(table_applicable(i+NHIGGS*(CDF9290_D05645 -1)))*( 
     &			 R_sigma_HW(i)*BRhjbbbar(i)*SMCS_tev_pp_qq_HW(Mh(i))*BR_Wlnu 
     &			+R_sigma_HW(i)*BRhjWW(i)*SMCS_tev_pp_qq_HW(Mh(i))*BR_Wlnu**2 
     &			+R_sigma_HZ(i)*BRhjbbbar(i)*SMCS_tev_pp_qq_HZ(Mh(i))*(BR_Znunu+BR_Zll) 
     &			+R_sigma_H(i)*BRhjWW(i)*( SMCS_tev_pp_gg_H(Mh(i)) + SMCS_tev_pp_bb_H(Mh(i)) ) 
     &			  *BR_Wlnu**2 
     &			+( R_sigma_HW(i)*BR_Wjj*SMCS_tev_pp_qq_HW(Mh(i)) 
     &			+R_sigma_HZ(i)*BR_Zjj*SMCS_tev_pp_qq_HZ(Mh(i)) 
     &			+R_sigma_H(i)*( SMCS_tev_pp_gg_H(Mh(i)) + SMCS_tev_pp_bb_H(Mh(i)) ) 
     &			+R_sigma_H_VBF(i)*SMCS_tev_pp_vbf_H(Mh(i)) ) 
     &			*BRhjtautau(i)
     &			+( R_sigma_H(i)*( SMCS_tev_pp_gg_H(Mh(i))  + SMCS_tev_pp_bb_H(Mh(i)) )
     &          	      +R_sigma_HW(i)*SMCS_tev_pp_qq_HW(Mh(i)) 
     &          	      +R_sigma_HZ(i)*SMCS_tev_pp_qq_HZ(Mh(i)) 
     &          	      +R_sigma_H_VBF(i)*SMCS_tev_pp_vbf_H(Mh(i)) ) 
     &          	       *BRhjgamgam(i)
     &			  	)

! CDF9465_D05754 : p p-bar -> H/WH/ZH/H via VBF, H -> W W
 	   	  R_sig_BR_sing(i+NHIGGS*(CDF9465_D05754 -1))=
     &			 log2dble(table_applicable(i+NHIGGS*(CDF9465_D05754 -1)))*( 
     &			 R_sigma_H(i)*( SMCS_tev_pp_gg_H(Mh(i)) + SMCS_tev_pp_bb_H(Mh(i)) ) 
     &			+R_sigma_HW(i)*SMCS_tev_pp_qq_HW(Mh(i))
     &                  +R_sigma_HZ(i)*SMCS_tev_pp_qq_HZ(Mh(i))
     &                  +R_sigma_H_VBF(i)*SMCS_tev_pp_vbf_H(Mh(i))
     &			  )*BRhjWW(i)

! CDF9674: p p-bar -> (WH/ZH, H-> b b-bar)/(H/WH, H -> W W)
!                  /(H/WH/ZH/H via VBF, H->tau tau)
 	   	  R_sig_BR_sing(i+NHIGGS*(CDF9674 -1))=
     &			 log2dble(table_applicable(i+NHIGGS*(CDF9674 -1)))*( 
     &			 R_sigma_HW(i)*BRhjbbbar(i)*SMC	S_tev_pp_qq_HW(Mh(i))*BR_Wlnu 
     &			+R_sigma_HW(i)*BRhjWW(i)*SMCS_tev_pp_qq_HW(Mh(i))*BR_Wlnu**2 
     &			+R_sigma_HZ(i)*BRhjbbbar(i)*SMCS_tev_pp_qq_HZ(Mh(i))*(BR_Znunu+BR_Zll) 
     &			+R_sigma_H(i)*BRhjWW(i)*( SMCS_tev_pp_gg_H(Mh(i)) + SMCS_tev_pp_bb_H(Mh(i)) ) 
     &			  *BR_Wlnu**2 
     &			+( R_sigma_HW(i)*BR_Wjj*SMCS_tev_pp_qq_HW(Mh(i)) 
     &			  +R_sigma_HZ(i)*BR_Zjj*SMCS_tev_pp_qq_HZ(Mh(i)) 
     &			  +R_sigma_H(i)*( SMCS_tev_pp_gg_H(Mh(i)) + SMCS_tev_pp_bb_H(Mh(i)) ) 
     &			  +R_sigma_H_VBF(i)*SMCS_tev_pp_vbf_H(Mh(i)) 
     &			  )*BRhjtautau(i)
     &			  	)

! CDF9713_D05889 : 
! p p-bar 
! -> HW -> 	b b l nu, 					not yet needed
!	   	b b-bar tau nu_tau, 				not yet needed
!		W+ W- W -> l+- nu l+- nu			x
! -> HV -> 	j j tau tau/tau tau j j (for m_H <= 145 GeV), 	not yet needed
!		b b + miss. E_T,				not yet needed
! 		b b j j 	(V=W,Z)				not yet needed
! -> HZ -> b b l l, 						not yet needed
! -> H via single H/HW/HZ/VBF + 0/1/2 jets, H -> W+ W- -> l+ nu l- nu,		x
! -> Htt -> b b t t (for m_H <= 155 GeV)			not yet needed
! -> H via single H/WH/ZH/VBF, 
!	H -> tau tau, 						not yet needed
!	H -> gamma gamma					not yet needed
 	   	  R_sig_BR_sing(i+NHIGGS*(CDF9713_D05889 -1))=
     &			 log2dble(table_applicable(i+NHIGGS*(CDF9713_D05889 -1)))*( 
     &			 ( R_sigma_H(i)*( SMCS_tev_pp_gg_H(Mh(i)) + SMCS_tev_pp_bb_H(Mh(i)) ) 
     &			  +R_sigma_HW(i)*SMCS_tev_pp_qq_HW(Mh(i))
     &                    +R_sigma_HZ(i)*SMCS_tev_pp_qq_HZ(Mh(i))
     &                    +R_sigma_H_VBF(i)*SMCS_tev_pp_vbf_H(Mh(i))
     &			 		) *BRhjWW(i)*BR_Wlnu**2
     &			 	)

! CDF9897: : p p-bar -> HW -> l nu b b-bar,
!			     HZ -> nu nu b b-bar, l l b b-bar,
!			     HV -> b b-bar j j,
!			     H via single H/HW/HZ/VBF, H -> W W
!			     H via single H/HW/HZ/VBF, H -> tau tau 
!			     SM combined 
 	   	  R_sig_BR_sing(i+NHIGGS*(CDF9897 -1))=
     &			 log2dble(table_applicable(i+NHIGGS*(CDF9897 -1)))*( 
     &			 ( R_sigma_HW(i)*SMCS_tev_pp_qq_HW(Mh(i))*(BR_Wlnu+BR_Wjj)
     &			  +R_sigma_HZ(i)*SMCS_tev_pp_qq_HZ(Mh(i))*(BR_Znunu+BR_Zll+BR_Zjj) 
     &			   )*BRhjbbbar(i)
     &			+( R_sigma_H(i)*( SMCS_tev_pp_gg_H(Mh(i)) + SMCS_tev_pp_bb_H(Mh(i)) ) 
     &			  +R_sigma_HW(i)*SMCS_tev_pp_qq_HW(Mh(i))
     &                    +R_sigma_HZ(i)*SMCS_tev_pp_qq_HZ(Mh(i))
     &                    +R_sigma_H_VBF(i)*SMCS_tev_pp_vbf_H(Mh(i))
     &			   )*BRhjWW(i)*BR_Wlnu**2
     &			+( R_sigma_HW(i)*SMCS_tev_pp_qq_HW(Mh(i))*BR_Wjj
     &			  +R_sigma_HZ(i)*SMCS_tev_pp_qq_HZ(Mh(i))*BR_Zjj
     &			  +R_sigma_H(i)*( SMCS_tev_pp_gg_H(Mh(i)) + SMCS_tev_pp_bb_H(Mh(i)) ) 
     &			  +R_sigma_H_VBF(i)*SMCS_tev_pp_vbf_H(Mh(i)) 
     &			   )*BRhjtautau(i)
     &			  	)

	   endif
        enddo


!    double Higgs channels: not needed because there are (so far) double Higgs channels



! a) calc sum of sigma x BR contributions to R_sig_BR
!   single Higgs channels
	do n=1,NHIGGS
!tag3
	    if( (Mh(n) .ge. Mhmin_SMcalc) .and. (Mh(n) .le. Mhmax_SMcalc) ) then 
!tag2
	      do i=1,NHIGGS
!tag5
	        if( (Mh(i) .ge. Mh(n)) .and. ( (Mh(i)-Mh(n)) .le. delta_Mh) )then
!tag4
! D05972, CDFwhlnubb_090814, D00808_1970, CDF0906_5613 : p p-bar -> WH, H -> b b-bar, W -> l nu
                if(R_sig_BR_sing(i+NHIGGS*(D05972 -1))*R_sig_BR_sing(n+NHIGGS*(D05972 -1)) .ne. 0d0) then
	   	  R_sig_BR(n+NHIGGS*(D05972 -1))= R_sig_BR(n+NHIGGS*(D05972 -1)) 
     &			+R_sig_BR_sing(i+NHIGGS*(D05972 -1))
	   	  R_sig_BR(n+NHIGGS*(CDFwhlnubb_090814 -1))=R_sig_BR(n+NHIGGS*(D05972 -1))
	   	  R_sig_BR(n+NHIGGS*(D00808_1970 -1))=R_sig_BR(n+NHIGGS*(D05972 -1))
	   	  R_sig_BR(n+NHIGGS*(CDF0906_5613 -1))=R_sig_BR(n+NHIGGS*(D05972 -1))
		endif

! CDF0908_3534, CDF9889, D05876 : p p-bar -> ZH -> l l b b-bar
		if(R_sig_BR_sing(i+NHIGGS*(CDF0908_3534 -1))*R_sig_BR_sing(n+NHIGGS*(CDF0908_3534 -1)) .ne. 0d0) then
	          R_sig_BR(n+NHIGGS*(CDF0908_3534 -1))=R_sig_BR(n+NHIGGS*(CDF0908_3534 -1))
     &			+R_sig_BR_sing(i+NHIGGS*(CDF0908_3534 -1))
	   	  R_sig_BR(n+NHIGGS*(CDF9889 -1))= R_sig_BR(n+NHIGGS*(CDF0908_3534 -1))
           	  R_sig_BR(n+NHIGGS*(D05876 -1))=R_sig_BR(n+NHIGGS*(CDF0908_3534 -1))
		endif

! D05757, CDF0809_3930 : p p-bar -> H -> W W
		if(R_sig_BR_sing(i+NHIGGS*(D05757 -1))*R_sig_BR_sing(n+NHIGGS*(D05757 -1)) .ne. 0d0) then
           	  R_sig_BR(n+NHIGGS*(D05757 -1))= R_sig_BR(n+NHIGGS*(D05757 -1)) 
     &			+R_sig_BR_sing(i+NHIGGS*(D05757 -1))
	   	  R_sig_BR(n+NHIGGS*(CDF0809_3930 -1))=R_sig_BR(n+NHIGGS*(D05757 -1))
		endif

! CDF9887, D05871 : p p-bar -> H + X -> W W + X
		if(R_sig_BR_sing(i+NHIGGS*(CDF9887 -1))*R_sig_BR_sing(n+NHIGGS*(CDF9887 -1)) .ne. 0d0) then
           	  R_sig_BR(n+NHIGGS*(CDF9887 -1))= R_sig_BR(n+NHIGGS*(CDF9887 -1)) 
     &			+R_sig_BR_sing(i+NHIGGS*(CDF9887 -1))
		  R_sig_BR(n+NHIGGS*(D05871 -1))=R_sig_BR(n+NHIGGS*(CDF9887 -1))
		endif


! CDF9284, D05726 : p p-bar -> b H, H -> b b-bar ([sigma x BR/sigma limit)
		if(R_sig_BR_sing(i+NHIGGS*(CDF9284 -1))*R_sig_BR_sing(n+NHIGGS*(CDF9284 -1)) .ne. 0d0) then
	   	  R_sig_BR(n+NHIGGS*(CDF9284 -1))=R_sig_BR(n+NHIGGS*(CDF9284 -1)) 
     &			+R_sig_BR_sing(i+NHIGGS*(CDF9284 -1))
	   	  R_sig_BR(n+NHIGGS*(D05726 -1))=R_sig_BR(n+NHIGGS*(CDF9284 -1))
		endif

! D00805_3556 : p p-bar -> b H, H -> b b-bar (absolute sigma x BR limit)
		if(R_sig_BR_sing(i+NHIGGS*(D00805_3556 -1))*R_sig_BR_sing(n+NHIGGS*(D00805_3556 -1)) .ne. 0d0) then
	   	  R_sig_BR(n+NHIGGS*(D00805_3556 -1))=R_sig_BR(n+NHIGGS*(D00805_3556 -1)) 
     &			+R_sig_BR_sing(i+NHIGGS*(D00805_3556 -1))
		endif

! 16 : p p-bar -> b H, H -> b b-bar (absolute sigma x BR limit)
!      the SM c.s. quoted in the CDF paper is a factor 0.7 lower than my c.s.
!      The CDF note 9284 result is also implemented as tab 10.
!      Table 16 was for checking whether the the absolute limit is consistent.
!      It is! So it should be switched off again.
!		if(R_sig_BR_sing(i+NHIGGS*15)*R_sig_BR_sing(n+NHIGGS*15) .ne. 0d0) then
!	   	  R_sig_BR(n+NHIGGS*15)=R_sig_BR(n+NHIGGS*15) 
!     &			+R_sig_BR_sing(i+NHIGGS*15)
!	   	  R_sig_BR(n+NHIGGS*15)=0d0
!		endif

! D00805_2491, CDF0906_1014, D05740, CDF9888_D05980 : p p-bar -> H -> tau tau (absolute sigma x BR limit)
                if(R_sig_BR_sing(i+NHIGGS*(D00805_2491 -1))*R_sig_BR_sing(n+NHIGGS*(D00805_2491 -1)) .ne. 0d0) then
                  R_sig_BR(n+NHIGGS*(D00805_2491 -1))=R_sig_BR(n+NHIGGS*(D00805_2491 -1)) 
     &                  +R_sig_BR_sing(i+NHIGGS*(D00805_2491 -1))
                  R_sig_BR(n+NHIGGS*(CDF0906_1014 -1))=R_sig_BR(n+NHIGGS*(D00805_2491 -1))
                  R_sig_BR(n+NHIGGS*(D05740 -1))=R_sig_BR(n+NHIGGS*(D00805_2491 -1))
                  R_sig_BR(n+NHIGGS*(CDF9888_D05980 -1))=R_sig_BR(n+NHIGGS*(D00805_2491 -1))
                endif

! CDF7307v3, D05873 : p p-bar -> WH, H -> WW 
		if(R_sig_BR_sing(i+NHIGGS*(CDF7307v3 -1))*R_sig_BR_sing(n+NHIGGS*(CDF7307v3 -1)) .ne. 0d0) then
	   	  R_sig_BR(n+NHIGGS*(CDF7307v3 -1))= R_sig_BR(n+NHIGGS*(CDF7307v3 -1)) 
     &			+R_sig_BR_sing(i+NHIGGS*(CDF7307v3 -1))
		  R_sig_BR(n+NHIGGS*(D05873 -1))=R_sig_BR(n+NHIGGS*(CDF7307v3 -1))
		endif

! CDF9891, D05586, D00808_1266, CDF0802_0432 : p p-bar -> WH/ZH, H -> b b-bar
		if(R_sig_BR_sing(i+NHIGGS*(CDF9891 -1))*R_sig_BR_sing(n+NHIGGS*(CDF9891 -1)) .ne. 0d0) then
	   	  R_sig_BR(n+NHIGGS*(CDF9891 -1))=R_sig_BR(n+NHIGGS*(CDF9891 -1)) 
     &			+R_sig_BR_sing(i+NHIGGS*(CDF9891 -1))
           	  R_sig_BR(n+NHIGGS*(D05586 -1)) =R_sig_BR(n+NHIGGS*(CDF9891 -1))
           	  R_sig_BR(n+NHIGGS*(D00808_1266 -1)) =R_sig_BR(n+NHIGGS*(CDF9891 -1))
           	  R_sig_BR(n+NHIGGS*(CDF0802_0432 -1)) =R_sig_BR(n+NHIGGS*(CDF9891 -1))
		endif

! D05858, D00901_1887 : p p-bar -> H -> gamma gamma
		if(R_sig_BR_sing(i+NHIGGS*(D05858 -1))*R_sig_BR_sing(n+NHIGGS*(D05858 -1)) .ne. 0d0) then
           	  R_sig_BR(n+NHIGGS*(D05858 -1))=R_sig_BR(n+NHIGGS*(D05858 -1)) 
     &			+R_sig_BR_sing(i+NHIGGS*(D05858 -1))
                  R_sig_BR(n+NHIGGS*(D00901_1887 -1)) =R_sig_BR(n+NHIGGS*(D05858 -1))
		endif

!temporary removed ! D00803_1514 : p p-bar -> H -> gamma gamma (absolute sigma x BR limit)
!temporary removed 		if(R_sig_BR_sing(i+NHIGGS*(D00803_1514 -1))*R_sig_BR_sing(n+NHIGGS*(D00803_1514 -1)) .ne. 0d0) then
!temporary removed            	  R_sig_BR(n+NHIGGS*(D00803_1514 -1))=R_sig_BR(n+NHIGGS*(D00803_1514 -1)) 
!temporary removed      &			+R_sig_BR_sing(i+NHIGGS*(D00803_1514 -1))
!temporary removed 		endif

! D00811_0024, D05985 : p p-bar -> b H, H -> tau tau (absolute sigma x BR limit)
		 if(R_sig_BR_sing(i+NHIGGS*(D00811_0024 -1))*R_sig_BR_sing(n+NHIGGS*(D00811_0024 -1)) .ne. 0d0) then
	   	  R_sig_BR(n+NHIGGS*(D00811_0024 -1))=R_sig_BR(n+NHIGGS*(D00811_0024 -1))
     &                                +R_sig_BR_sing(i+NHIGGS*(D00811_0024 -1))
                  R_sig_BR(n+NHIGGS*(D05985 -1)) =R_sig_BR(n+NHIGGS*(D00811_0024 -1))
		endif



! SM combined analyses
! CDF8961_D05536, D00712_0598 : p p-bar -> (WH/ZH, H-> b b-bar)/(H,WH, H -> W W)

		if(R_sig_BR_sing(i+NHIGGS*(CDF8961_D05536 -1))*R_sig_BR_sing(n+NHIGGS*(CDF8961_D05536 -1)) .ne. 0d0) then
 	   	  R_sig_BR(n+NHIGGS*(CDF8961_D05536 -1))= R_sig_BR(n+NHIGGS*(CDF8961_D05536 -1)) 
     &			+R_sig_BR_sing(i+NHIGGS*(CDF8961_D05536 -1))
		   R_sig_BR(n+NHIGGS*(D00712_0598 -1)) =R_sig_BR(n+NHIGGS*(CDF8961_D05536 -1))
		endif

! 24 : p p-bar -> (single H, H via VBF ,WH/ZH)

		if(R_sig_BR_sing(i+NHIGGS*(CDF9248 -1))*R_sig_BR_sing(n+NHIGGS*(CDF9248 -1)) .ne. 0d0) then
 	   	  R_sig_BR(n+NHIGGS*(CDF9248 -1))= R_sig_BR(n+NHIGGS*(CDF9248 -1)) 
     &			+R_sig_BR_sing(i+NHIGGS*(CDF9248 -1))
		endif

! 25 : p p-bar -> (WH/ZH, H-> b b-bar)/(H/WH, H -> W W)
!                  /(H/WH/ZH/H via VBF, H->tau tau)
!		   /(H/WH/ZH/H via VBF,H -> gamma gamma)

		if(R_sig_BR_sing(i+NHIGGS*(CDF9290_D05645 -1))*R_sig_BR_sing(n+NHIGGS*(CDF9290_D05645 -1)) .ne. 0d0) then
 	   	  R_sig_BR(n+NHIGGS*(CDF9290_D05645 -1))= R_sig_BR(n+NHIGGS*(CDF9290_D05645 -1)) 
     &			+R_sig_BR_sing(i+NHIGGS*(CDF9290_D05645 -1))
		endif

! 28 : p p-bar -> H/WH/ZH/H via VBF, H -> W W
		if(R_sig_BR_sing(i+NHIGGS*(CDF9465_D05754 -1))*R_sig_BR_sing(n+NHIGGS*(CDF9465_D05754 -1)) .ne. 0d0) then
 	   	  R_sig_BR(n+NHIGGS*(CDF9465_D05754 -1))= R_sig_BR(n+NHIGGS*(CDF9465_D05754 -1)) 
     &                 +R_sig_BR_sing(i+NHIGGS*(CDF9465_D05754 -1))		
		endif

! CDF9674 : p p-bar -> (WH/ZH, H-> b b-bar)/(H/WH, H -> W W)
!                  /(H/WH/ZH/H via VBF, H->tau tau)
		if(R_sig_BR_sing(i+NHIGGS*(CDF9674 -1))*R_sig_BR_sing(n+NHIGGS*(CDF9674 -1)) .ne. 0d0) then
 	   	  R_sig_BR(n+NHIGGS*(CDF9674 -1))= R_sig_BR(n+NHIGGS*(CDF9674 -1)) 
     &			+R_sig_BR_sing(i+NHIGGS*(CDF9674 -1))
		endif

! CDF9713_D05889 : 
! p p-bar 
! -> HW -> 	b b l nu, 					not yet needed
!	   	b b-bar tau nu_tau, 				not yet needed
!		W+ W- W -> l+- nu l+- nu			x
! -> HV -> 	j j tau tau/tau tau j j (for m_H <= 145 GeV), 	not yet needed
!		b b + miss. E_T,				not yet needed
! 		b b j j 	(V=W,Z)				not yet needed
! -> HZ -> b b l l, 						not yet needed
! -> H via single H/HW/HZ/VBF + 0/1/2 jets, H -> W+ W- -> l+ nu l- nu,		x
! -> Htt -> b b t t (for m_H <= 155 GeV)			not yet needed
! -> H via single H/WH/ZH/VBF, 
!	H -> tau tau, 						not yet needed
!	H -> gamma gamma					not yet needed
		if(R_sig_BR_sing(i+NHIGGS*(CDF9713_D05889 -1))*R_sig_BR_sing(n+NHIGGS*(CDF9713_D05889 -1)) .ne. 0d0) then
 	   	  R_sig_BR(n+NHIGGS*(CDF9713_D05889 -1))= R_sig_BR(n+NHIGGS*(CDF9713_D05889 -1)) 
     &			+R_sig_BR_sing(i+NHIGGS*(CDF9713_D05889 -1))
		endif

!CDF9897	: p p-bar -> HW -> l nu b b-bar,
!			     HZ -> nu nu b b-bar, l l b b-bar,
!			     HV -> b b-bar j j,
!			     H via single H/HW/HZ/VBF, H -> W W
!			     H via single H/HW/HZ/VBF, H -> tau tau 
!			     SM combined 
		if(R_sig_BR_sing(i+NHIGGS*(CDF9897 -1))*R_sig_BR_sing(n+NHIGGS*(CDF9897 -1)) .ne. 0d0) then
 	   	  R_sig_BR(n+NHIGGS*(CDF9897 -1))= R_sig_BR(n+NHIGGS*(CDF9897 -1)) 
     &			+R_sig_BR_sing(i+NHIGGS*(CDF9897 -1))
		endif



	        endif
!tag4
	      enddo
!tag5
	    endif	  
!tag2
        enddo
!tag3

! b) Determine the average Higgs mass for the Higgs combinations

	do m=0,TEV_NSCHN-1,1
	  do n=1,NHIGGS
	    sim_mass_sum(n+NHIGGS*m)=1
            higgs_combination(n+NHIGGS*m,n)=.TRUE.
	    Mh_avg(n+NHIGGS*m)=Mh(n)
	    do i=1,NHIGGS
              if(i .ne. n) then
	        if( (Mh(i) .ge. Mh(n)) .and. ( (Mh(i)-Mh(n)) .le. delta_Mh )
     &	 	    .and. ( R_sig_BR_sing(i+NHIGGS*m)*R_sig_BR_sing(n+NHIGGS*m) .ne. 0d0) 
     &		    .and.  table_applicable(i+NHIGGS*m) ) then
	          Mh_avg(n+NHIGGS*m)=Mh_avg(n+NHIGGS*m)+Mh(i)
                  sim_mass_sum(n+NHIGGS*m)=sim_mass_sum(n+NHIGGS*m)+1
	          higgs_combination(n+NHIGGS*m,i)=.TRUE.
	        endif
	      endif
	    enddo
            Mh_avg(n+NHIGGS*m)=Mh_avg(n+NHIGGS*m)/dble(sim_mass_sum(n+NHIGGS*m))
          enddo
	enddo


! c) multiply sigma x BR by normalization 

	do m=0,TEV_NSCHN-1,1
         do n=1,NHIGGS
	   if( (Mh_avg(n+NHIGGS*m) .ge. Mhmin_SMcalc) 
     &		.and. (Mh_avg(n+NHIGGS*m) .le. Mhmax_SMcalc) ) then 
	  if(m .eq. (D05972 -1)) then
! D05972, CDFwhlnubb_090814, D00808_1970, CDF0906_5613 : p p-bar -> WH, H -> b b-bar
	    R_sig_BR(n+NHIGGS*(D05972 -1))=R_sig_BR(n+NHIGGS*(D05972 -1))
     &		/SMCS_tev_pp_qq_HW(Mh_avg(n+NHIGGS*(D05972 -1)))/SMBR_Hbb(Mh_avg(n+NHIGGS*(D05972 -1)))

	  elseif(m .eq. (CDF0908_3534 -1)) then
! CDF0908_3534, CDF9889, D05876 : p p-bar -> ZH, H -> b b-bar
	    R_sig_BR(n+NHIGGS*(CDF0908_3534 -1))= R_sig_BR(n+NHIGGS*(CDF0908_3534 -1))
     &		/SMCS_tev_pp_qq_HZ(Mh_avg(n+NHIGGS*(CDF0908_3534 -1)))/SMBR_Hbb(Mh_avg(n+NHIGGS*(CDF0908_3534 -1)))

	  elseif(m .eq. (D05757 -1)) then
! D05757, CDF0809_3930 : p p-bar -> H -> W W
	    R_sig_BR(n+NHIGGS*(D05757 -1))= R_sig_BR(n+NHIGGS*(D05757 -1))
     &         	/( SMCS_tev_pp_gg_H(Mh_avg(n+NHIGGS*(D05757 -1))) 
     &		  +SMCS_tev_pp_bb_H(Mh_avg(n+NHIGGS*(D05757 -1))) )/SMBR_HWW(Mh_avg(n+NHIGGS*(D05757 -1)))

! CDF9887, D05871 : p p-bar -> H + X -> W W + X
	    R_sig_BR(n+NHIGGS*(CDF9887 -1))= R_sig_BR(n+NHIGGS*(CDF9887 -1))
     &         	/( SMCS_tev_pp_gg_H(Mh_avg(n+NHIGGS*(CDF9887 -1))) 
     &		  +SMCS_tev_pp_bb_H(Mh_avg(n+NHIGGS*(CDF9887 -1))) 
     &            +SMCS_tev_pp_qq_HZ(Mh_avg(n+NHIGGS*(CDF9887 -1)))
     &            +SMCS_tev_pp_qq_HW(Mh_avg(n+NHIGGS*(CDF9887 -1)))
     &            +SMCS_tev_pp_vbf_H(Mh_avg(n+NHIGGS*(CDF9887 -1)))
     &		  )/SMBR_HWW(Mh_avg(n+NHIGGS*(CDF9887 -1)))

	  elseif(m .eq. (CDF9284 -1)) then
! CDF9284, D05726 : p p-bar -> b H, H -> b b-bar ([sigma x BR/sigma limit)
            R_sig_BR(n+NHIGGS*(CDF9284 -1))=R_sig_BR(n+NHIGGS*(CDF9284 -1))
     &		/SMCS_tev_pp_bg_Hb(Mh_avg(n+NHIGGS*(CDF9284 -1)))

          elseif(m .eq. (CDF7307v3 -1)) then
! CDF7307v3, D05873 : p p-bar -> WH, H -> WW
            R_sig_BR(n+NHIGGS*(CDF7307v3 -1))= R_sig_BR(n+NHIGGS*(CDF7307v3 -1))
     &          /SMCS_tev_pp_qq_HW(Mh_avg(n+NHIGGS*(CDF7307v3 -1)))/SMBR_HWW(Mh_avg(n+NHIGGS*(CDF7307v3 -1)))

          elseif(m .eq. (CDF9891 -1)) then
! CDF9891, D05586, D00808_1266, CDF0802_0432 : p p-bar -> WH/ZH, H -> b b-bar
           R_sig_BR(n+NHIGGS*(CDF9891 -1))=R_sig_BR(n+NHIGGS*(CDF9891 -1))
     &		/( ( SMCS_tev_pp_qq_HW(Mh_avg(n+NHIGGS*(CDF9891 -1)))*BR_Wlnu
     &               +SMCS_tev_pp_qq_HZ(Mh_avg(n+NHIGGS*(CDF9891 -1)))*BR_Znunu 
     &		   )*SMBR_Hbb(Mh_avg(n+NHIGGS*(CDF9891 -1))) )

	  elseif(m .eq. (D05858 -1)) then
! D05858, D00901_1887 : p p-bar -> H -> gamma gamma
	  R_sig_BR(n+NHIGGS*(D05858 -1))=R_sig_BR(n+NHIGGS*(D05858 -1))
     &		/(
     &		  (  SMCS_tev_pp_gg_H(Mh_avg(n+NHIGGS*(D05858 -1)))  
     &	 	   + SMCS_tev_pp_bb_H(Mh_avg(n+NHIGGS*(D05858 -1)))
     &	 	   + SMCS_tev_pp_qq_HW(Mh_avg(n+NHIGGS*(D05858 -1)))
     &	 	   + SMCS_tev_pp_qq_HZ(Mh_avg(n+NHIGGS*(D05858 -1)))
     &	 	   + SMCS_tev_pp_vbf_H(Mh_avg(n+NHIGGS*(D05858 -1))) 
     &		  )*SMBR_Hgamgam(Mh_avg(n+NHIGGS*(D05858 -1)))
     &		 )
	
! SM combined analyses

          elseif(m .eq. (CDF8961_D05536 -1)) then
! CDF8961_D05536, D00712_0598 : p p-bar -> (WH/ZH, H-> b b-bar)/(H,WH, H -> W W)
           sum_sig_BR_SM=
     &           SMCS_tev_pp_qq_HW(Mh_avg(n+NHIGGS*(CDF8961_D05536 -1)))*( SMBR_HWW(Mh_avg(n+NHIGGS*(CDF8961_D05536 -1)))*BR_Wlnu
     &				+SMBR_Hbb(Mh_avg(n+NHIGGS*(CDF8961_D05536 -1))) )*BR_Wlnu
     &          +SMCS_tev_pp_qq_HZ(Mh_avg(n+NHIGGS*(CDF8961_D05536 -1)))*SMBR_Hbb(Mh_avg(n+NHIGGS*(CDF8961_D05536 -1)))*(BR_Znunu+BR_Zll)
     &          +( SMCS_tev_pp_gg_H(Mh_avg(n+NHIGGS*(CDF8961_D05536 -1))) 
     &			+ SMCS_tev_pp_bb_H(Mh_avg(n+NHIGGS*(CDF8961_D05536 -1))) )*SMBR_HWW(Mh_avg(n+NHIGGS*(CDF8961_D05536 -1)))*BR_Wlnu**2
           R_sig_BR(n+NHIGGS*(CDF8961_D05536 -1))= R_sig_BR(n+NHIGGS*(CDF8961_D05536 -1))/sum_sig_BR_SM

          elseif(m .eq. (CDF9248 -1)) then
! CDF9248 : p p-bar -> (single H, H via VBF ,WH/ZH)
           sum_sig_BR_SM_tautau=
     &           (  SMCS_tev_pp_qq_HW(Mh_avg(n+NHIGGS*(CDF9248 -1)))*BR_Wjj
     &             +SMCS_tev_pp_qq_HZ(Mh_avg(n+NHIGGS*(CDF9248 -1)))*BR_Zjj
     &             +SMCS_tev_pp_gg_H(Mh_avg(n+NHIGGS*(CDF9248 -1))) + SMCS_tev_pp_bb_H(Mh_avg(n+NHIGGS*(CDF9248 -1)))
     &             +SMCS_tev_pp_vbf_H(Mh_avg(n+NHIGGS*(CDF9248 -1))) )*SMBR_Htautau(Mh_avg(n+NHIGGS*(CDF9248 -1)))
           R_sig_BR(n+NHIGGS*(CDF9248 -1))= R_sig_BR(n+NHIGGS*(CDF9248 -1))/sum_sig_BR_SM_tautau

          elseif(m .eq. (CDF9290_D05645 -1)) then
! CDF9290_D05645 : p p-bar -> (WH/ZH, H-> b b-bar)/(H/WH, H -> W W)
!                  /(H/WH/ZH/H via VBF, H->tau tau)
!                  /(H/WH/ZH/H via VBF,H -> gamma gamma)
           sum_sig_BR_SM_2= 
     &           SMCS_tev_pp_qq_HW(Mh_avg(n+NHIGGS*(CDF9290_D05645 -1)))*( SMBR_HWW(Mh_avg(n+NHIGGS*(CDF9290_D05645 -1)))*BR_Wlnu
     &				+SMBR_Hbb(Mh_avg(n+NHIGGS*(CDF9290_D05645 -1))) )*BR_Wlnu
     &          +SMCS_tev_pp_qq_HZ(Mh_avg(n+NHIGGS*(CDF9290_D05645 -1)))*SMBR_Hbb(Mh_avg(n+NHIGGS*(CDF9290_D05645 -1)))*(BR_Znunu+BR_Zll)
     &          +( SMCS_tev_pp_gg_H(Mh_avg(n+NHIGGS*(CDF9290_D05645 -1))) 
     &			+ SMCS_tev_pp_bb_H(Mh_avg(n+NHIGGS*(CDF9290_D05645 -1))) )*SMBR_HWW(Mh_avg(n+NHIGGS*(CDF9290_D05645 -1)))*BR_Wlnu**2
     &           +(  
     &			 SMCS_tev_pp_qq_HW(Mh_avg(n+NHIGGS*(CDF9290_D05645 -1)))*BR_Wjj
     &             	+SMCS_tev_pp_qq_HZ(Mh_avg(n+NHIGGS*(CDF9290_D05645 -1)))*BR_Zjj
     &             	+SMCS_tev_pp_gg_H(Mh_avg(n+NHIGGS*(CDF9290_D05645 -1))) + SMCS_tev_pp_bb_H(Mh_avg(n+NHIGGS*(CDF9290_D05645 -1)))
     &             	+SMCS_tev_pp_vbf_H(Mh_avg(n+NHIGGS*(CDF9290_D05645 -1))) )*SMBR_Htautau(Mh_avg(n+NHIGGS*(CDF9290_D05645 -1)))
     &          +( 
     &			SMCS_tev_pp_gg_H(Mh_avg(n+NHIGGS*(CDF9290_D05645 -1)))  + SMCS_tev_pp_bb_H(Mh_avg(n+NHIGGS*(CDF9290_D05645 -1)))
     &                 +SMCS_tev_pp_qq_HW(Mh_avg(n+NHIGGS*(CDF9290_D05645 -1))) + SMCS_tev_pp_qq_HZ(Mh_avg(n+NHIGGS*(CDF9290_D05645 -1)))
     &                 +SMCS_tev_pp_vbf_H(Mh_avg(n+NHIGGS*(CDF9290_D05645 -1))) )*SMBR_Hgamgam(Mh_avg(n+NHIGGS*(CDF9290_D05645 -1)))
           R_sig_BR(n+NHIGGS*(CDF9290_D05645 -1))= R_sig_BR(n+NHIGGS*(CDF9290_D05645 -1))/sum_sig_BR_SM_2

          elseif(m .eq. (CDF9465_D05754 -1)) then
! CDF9465_D05754 : p p-bar -> H/WH/ZH/H via VBF, H -> W W
           sum_sig_BR_SM_3=
     &        (  SMCS_tev_pp_gg_H(Mh_avg(n+NHIGGS*(CDF9465_D05754 -1))) 
     &		+SMCS_tev_pp_bb_H(Mh_avg(n+NHIGGS*(CDF9465_D05754 -1))) 
     &		+SMCS_tev_pp_qq_HW(Mh_avg(n+NHIGGS*(CDF9465_D05754 -1)))
     &		+SMCS_tev_pp_qq_HZ(Mh_avg(n+NHIGGS*(CDF9465_D05754 -1)))
     &		+SMCS_tev_pp_vbf_H(Mh_avg(n+NHIGGS*(CDF9465_D05754 -1)))
     &	      )*SMBR_HWW(Mh_avg(n+NHIGGS*(CDF9465_D05754 -1)))
           R_sig_BR(n+NHIGGS*(CDF9465_D05754 -1))= R_sig_BR(n+NHIGGS*(CDF9465_D05754 -1))/sum_sig_BR_SM_3

          elseif(m .eq. (CDF9674 -1)) then
! CDF9674 : p p-bar -> (WH/ZH, H-> b b-bar)/(H/WH, H -> W W)
!                  /(H/WH/ZH/H via VBF, H->tau tau)
           sum_sig_BR_CDF_SM= 
     &           SMCS_tev_pp_qq_HW(Mh_avg(n+NHIGGS*(CDF9674 -1)))*( SMBR_HWW(Mh_avg(n+NHIGGS*(CDF9674 -1)))*BR_Wlnu
     &				+SMBR_Hbb(Mh_avg(n+NHIGGS*(CDF9674 -1))) )*BR_Wlnu
     &          +SMCS_tev_pp_qq_HZ(Mh_avg(n+NHIGGS*(CDF9674 -1)))*SMBR_Hbb(Mh_avg(n+NHIGGS*(CDF9674 -1)))*(BR_Znunu+BR_Zll)
     &          +( SMCS_tev_pp_gg_H(Mh_avg(n+NHIGGS*(CDF9674 -1))) 
     &			+ SMCS_tev_pp_bb_H(Mh_avg(n+NHIGGS*(CDF9674 -1))) )*SMBR_HWW(Mh_avg(n+NHIGGS*(CDF9674 -1)))*BR_Wlnu**2
     &           +(  
     &			 SMCS_tev_pp_qq_HW(Mh_avg(n+NHIGGS*(CDF9674 -1)))*BR_Wjj
     &             	+SMCS_tev_pp_qq_HZ(Mh_avg(n+NHIGGS*(CDF9674 -1)))*BR_Zjj
     &             	+SMCS_tev_pp_gg_H(Mh_avg(n+NHIGGS*(CDF9674 -1))) + SMCS_tev_pp_bb_H(Mh_avg(n+NHIGGS*(CDF9674 -1)))
     &             	+SMCS_tev_pp_vbf_H(Mh_avg(n+NHIGGS*(CDF9674 -1))) )*SMBR_Htautau(Mh_avg(n+NHIGGS*(CDF9674 -1)))
           R_sig_BR(n+NHIGGS*(CDF9674 -1))= R_sig_BR(n+NHIGGS*(CDF9674 -1))/sum_sig_BR_CDF_SM

          elseif(m .eq. (CDF9713_D05889 -1)) then
! CDF9713_D05889 : 
! p p-bar 
! -> HW -> 	b b l nu, 					not yet needed
!	   	b b-bar tau nu_tau, 				not yet needed
!		W+ W- W -> l+- nu l+- nu			x
! -> HV -> 	j j tau tau/tau tau j j (for m_H <= 145 GeV), 	not yet needed
!		b b + miss. E_T,				not yet needed
! 		b b j j 	(V=W,Z)				not yet needed
! -> HZ -> b b l l, 						not yet needed
! -> H via single H/HW/HZ/VBF + 0/1/2 jets, H -> W+ W- -> l+ nu l- nu,		x
! -> Htt -> b b t t (for m_H <= 155 GeV)			not yet needed
! -> H via single H/WH/ZH/VBF, 
!	H -> tau tau, 						not yet needed
!	H -> gamma gamma					not yet needed
           sum_sig_BR_SM_4= 
     &          (  SMCS_tev_pp_gg_H(Mh_avg(n+NHIGGS*(CDF9713_D05889 -1)))  + SMCS_tev_pp_bb_H(Mh_avg(n+NHIGGS*(CDF9713_D05889 -1)))
     &             +SMCS_tev_pp_qq_HW(Mh_avg(n+NHIGGS*(CDF9713_D05889 -1))) 
     &		   +SMCS_tev_pp_qq_HZ(Mh_avg(n+NHIGGS*(CDF9713_D05889 -1)))
     &             +SMCS_tev_pp_vbf_H(Mh_avg(n+NHIGGS*(CDF9713_D05889 -1))) )*SMBR_HWW(Mh_avg(n+NHIGGS*(CDF9713_D05889 -1)))*BR_Wlnu**2
           R_sig_BR(n+NHIGGS*(CDF9713_D05889 -1))= R_sig_BR(n+NHIGGS*(CDF9713_D05889 -1))/sum_sig_BR_SM_4

          elseif(m .eq. (CDF9897 -1)) then
! CDF9897: : p p-bar -> HW -> l nu b b-bar,
!			     HZ -> nu nu b b-bar, l l b b-bar,
!			     HV -> b b-bar j j,
!			     H via single H/HW/HZ/VBF, H -> W W
!			     H via single H/HW/HZ/VBF, H -> tau tau 
!			     SM combined 
           sum_sig_BR_SM_5=
     &                   ( SMCS_tev_pp_qq_HW(Mh_avg(n+NHIGGS*(CDF9897 -1)))*(BR_Wlnu+BR_Wjj)
     &                    +SMCS_tev_pp_qq_HZ(Mh_avg(n+NHIGGS*(CDF9897 -1)))*(BR_Znunu+BR_Zll+BR_Zjj) 
     &                      )*SMBR_Hbb(Mh_avg(n+NHIGGS*(CDF9897 -1)))
     &                  +( SMCS_tev_pp_gg_H(Mh_avg(n+NHIGGS*(CDF9897 -1))) 
     &                    +SMCS_tev_pp_bb_H(Mh_avg(n+NHIGGS*(CDF9897 -1)))
     &                    +SMCS_tev_pp_qq_HW(Mh_avg(n+NHIGGS*(CDF9897 -1)))
     &                    +SMCS_tev_pp_qq_HZ(Mh_avg(n+NHIGGS*(CDF9897 -1)))
     &                    +SMCS_tev_pp_vbf_H(Mh_avg(n+NHIGGS*(CDF9897 -1)))
     &                     )*SMBR_HWW(Mh_avg(n+NHIGGS*(CDF9897 -1)))*BR_Wlnu**2
     &                  +( SMCS_tev_pp_qq_HW(Mh_avg(n+NHIGGS*(CDF9897 -1)))*BR_Wjj 
     &                    +SMCS_tev_pp_qq_HZ(Mh_avg(n+NHIGGS*(CDF9897 -1)))*BR_Zjj
     &                    +SMCS_tev_pp_gg_H(Mh_avg(n+NHIGGS*(CDF9897 -1))) 
     &                    +SMCS_tev_pp_bb_H(Mh_avg(n+NHIGGS*(CDF9897 -1)))
     &                    +SMCS_tev_pp_vbf_H(Mh_avg(n+NHIGGS*(CDF9897 -1))) 
     &                     )*SMBR_Htautau(Mh_avg(n+NHIGGS*(CDF9897 -1)))

           R_sig_BR(n+NHIGGS*(CDF9897 -1))= R_sig_BR(n+NHIGGS*(CDF9897 -1))/sum_sig_BR_SM_5
	  endif	

	 else
	  R_sig_BR(n+NHIGGS*m)=0d0
	 endif

         enddo
        enddo
c signals with same values than others already evaluated above
	do n=1,NHIGGS
! D05972, CDFwhlnubb_090814, D00808_1970 : p p-bar -> WH, H -> b b-bar
           R_sig_BR(n+NHIGGS*(CDFwhlnubb_090814 -1))=R_sig_BR(n+NHIGGS*(D05972 -1))
           R_sig_BR(n+NHIGGS*(D00808_1970 -1))=R_sig_BR(n+NHIGGS*(D05972 -1))
           R_sig_BR(n+NHIGGS*(CDF0906_5613 -1))=R_sig_BR(n+NHIGGS*(D05972 -1))
! CDF0908_3534, CDF9889, D05876 : p p-bar -> ZH, H -> b b-bar
           R_sig_BR(n+NHIGGS*(CDF9889 -1))=R_sig_BR(n+NHIGGS*(CDF0908_3534 -1))
           R_sig_BR(n+NHIGGS*(D05876 -1))=R_sig_BR(n+NHIGGS*(CDF0908_3534 -1))
! D05757, CDF0809_3930 : p p-bar -> H -> W W
	   R_sig_BR(n+NHIGGS*(CDF0809_3930 -1))=R_sig_BR(n+NHIGGS*(D05757 -1))
! CDF9887, D05871 
           R_sig_BR(n+NHIGGS*(D05871 -1))=R_sig_BR(n+NHIGGS*(CDF9887 -1))
! CDF9284, D05726 : p p-bar -> b H, H -> b b-bar ([sigma x BR/sigma limit)
           R_sig_BR(n+NHIGGS*(D05726 -1))=R_sig_BR(n+NHIGGS*(CDF9284 -1))
! CDF7307v3, D05873 : p p-bar -> WH, H -> WW
	    R_sig_BR(n+NHIGGS*(D05873 -1))=R_sig_BR(n+NHIGGS*(CDF7307v3 -1))
! CDF9891, D05586, D00808_1266, CDF0802_0432 : p p-bar -> WH/ZH, H -> b b-bar
           R_sig_BR(n+NHIGGS*(D05586 -1)) =R_sig_BR(n+NHIGGS*(CDF9891 -1))
           R_sig_BR(n+NHIGGS*(D00808_1266 -1)) =R_sig_BR(n+NHIGGS*(CDF9891 -1))
           R_sig_BR(n+NHIGGS*(CDF0802_0432 -1)) =R_sig_BR(n+NHIGGS*(CDF9891 -1))
! D05858, D00901_1887 : p p-bar -> H -> gamma gamma
           R_sig_BR(n+NHIGGS*(D00901_1887 -1)) =R_sig_BR(n+NHIGGS*(D05858 -1))
! CDF8961_D05536, D00712_0598 : p p-bar -> (WH/ZH, H-> b b-bar)/(H,WH, H -> W W)
	   R_sig_BR(n+NHIGGS*(D00712_0598 -1)) =R_sig_BR(n+NHIGGS*(CDF8961_D05536 -1))
	enddo



! double Higgs channels
!>> currently none
! see LEP constraints for structure
!<< currently none

	  allowed = .TRUE.

! comparing input to the predicted Tevatron results	   	  		      
	   do n=1,NHIGGS
	    do k=0,TEV_NSCHN-1 	   
#ifndef ONE_SELECTED
	     if(table_applicable(n+NHIGGS*k)) then
#else
	     if(table_applicable(n+NHIGGS*k) .and. (tlist_tev(n+NHIGGS*k) .eq. ONE_SELECTED)) then
#endif
	      call get_1d_TEVtable_data(tlist_tev(n+NHIGGS*k),tabsize,Mhmin,Mhmax,sep)
	      call check_TEV_bound1(Mh_avg(n+NHIGGS*k), 
     &	       R_sig_BR(n+NHIGGS*k),predratio(n+NHIGGS*k),tlist_tev(n+NHIGGS*k),
     &	       tabsize,Mhmin,Mhmax,sep,.True.)
	     else
	       predratio(n+NHIGGS*k)=-1d0
	     endif
	    enddo
	   enddo

!>> (double Higgs channels)
!>> currently none
! see LEP constraints for structure
!<< currently none

#ifdef DEBUG
	do nhelp=1,TEV_NSCHN*NHIGGS+TEV_NDCHN*NHIGGS**2,1
	    call determine_output_channel_number(NHIGGS,whichexpt,
     &          2,nhelp,output_chan)
	    call write_output_channel(6,NHIGGS,whichexpt,output_chan)
	    write(*,*)'   R_sig_BR=sig*BR/sig*BR_ref =',R_sig_BR(nhelp)
            write(*,*)'   R_sig_BR/S95_pred.         =',predratio(nhelp) 
	enddo
#endif

	  mlratio=my_maxloc(predratio,TEV_NSCHN*NHIGGS+TEV_NDCHN*NHIGGS**2)


! Determine the reaction number n connected to mlratio.
! For mlratio <= TEV_NSCHN*NHIGGS there are three reactions per table, 
! i.e. n is identical with the Higgs species.
! For mlratio >  TEV_NSCHN*NHIGGS there are NHIGGS**2-NHIGGS reactions per table,
! i.e. n lies between 1 and NHIGGS**2
	
	n=TEV_reaction_number(NHIGGS,mlratio)


! compare with R_sig_BR(mlratio) with experimental S95 value
	  if((mlratio .le. TEV_NSCHN*NHIGGS))then
	    call get_1d_TEVtable_data(tlist_tev(mlratio),tabsize,Mhmin,Mhmax,sep)
	    call check_TEV_bound1(Mh_avg(mlratio), 
     &	     R_sig_BR(mlratio),obsratio,tlist_tev(mlratio),
     &       tabsize,Mhmin,Mhmax,sep,.False.)        
	  else
!>> not needed
!	    i=findi(n)	      
!	    j=findj(n)
!	    call check_TEV_bound2(Mh(i),Mh(j), 
!     &	     R_sig_BR(mlratio),obsratio,tlist_tev(mlratio),.False.)
!<< not needed
	  endif

! obsratio = -1 is considered "allowed". This is where the tables give no exclusion, e.g. because the 
!	Higgs mass is out of range.
          if(obsratio .lt. 1d0) then
	    allowed=.True.
	    result_flag=1
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
        write(*,*)'Tevatron channel with highest sensitivity:'
	call determine_output_channel_number(NHIGGS,whichexpt,
     &          2,mlratio,output_chan)
	call write_output_channel(6,NHIGGS,whichexpt,output_chan)
        write(*,*)'   R_sig_BR/S95_predicted =',most_sens_predratio
        write(*,*)'   R_sig_BR/S95_observed  =',obsratio
#endif

	end


!***********************************************************
	subroutine calc_CS_and_BR_from_couplings(NHIGGS,Mh,Gamma_h,
     &		g_hjbb_sq, g_hjtautau_sq, g_hjWW_sq, g_hjZZ_sq,
     &          g_hjgamgam_sq, g_hjgg_sq, g_hjhiZ_sq_compact,
     &          R_hjZ_ee,R_hjhi_ee_compact,
     &          R_hj_gg,R_hj_bb,R_hjb_bg,
     &          R_hjWp_udb,R_hjWp_csb,R_hjWm_dub,R_hjWm_scb,
     &          R_hjZ_ddb,R_hjZ_uub,R_hjZ_ssb,R_hjZ_ccb,R_hjZ_bbb,
     &		R_sigma_H_VBF,
     &          BRhjbbbar,BRhjtautau,BRhjWW,BRhjgamgam)
!***********************************************************
	implicit none
        integer NHIGGS,i,j,k,n

        double precision g_hjbb_sq(NHMAX), g_hjtautau_sq(NHMAX),
     &          g_hjWW_sq(NHMAX), g_hjZZ_sq(NHMAX),
     &          g_hjgamgam_sq(NHMAX), g_hjgg_sq(NHMAX),
     &          g_hjhiZ_sq(NHMAX,NHMAX)
	double precision g_hjhiZ_sq_compact((NHMAX*(NHMAX+1))/2)

        double precision Mh(NHMAX),Gamma_h(NHMAX),R_Gamma
        double precision BRhjWW(NHMAX),BRhjgamgam(NHMAX)
        double precision BRhjbbbar(NHMAX),BRhjtautau(NHMAX)

        double precision R_hjZ_ee(NHMAX)
        double precision R_hjhi_ee(NHMAX,NHMAX)
	double precision R_hjhi_ee_compact((NHMAX*(NHMAX+1))/2)

        double precision R_hj_gg(NHMAX),R_hj_bb(NHMAX)
        double precision R_hjb_bg(NHMAX)
        double precision R_hjWp_udb(NHMAX),R_hjWp_csb(NHMAX),
     &                   R_hjWm_dub(NHMAX),R_hjWm_scb(NHMAX)
        double precision R_hjZ_ddb(NHMAX),R_hjZ_uub(NHMAX),
     &                   R_hjZ_ssb(NHMAX),R_hjZ_ccb(NHMAX),R_hjZ_bbb(NHMAX)
        double precision R_sigma_H_VBF(NHMAX)
        double precision SM_ZZ_contrib_to_VBF,SM_WW_contrib_to_VBF
! fit functions for SM width and branching ratios
        double precision SMGamma_h
        double precision SMBR_HWW
        double precision SMBR_Hbb
        double precision SMBR_Htautau
        double precision SMBR_Hgamgam
        double precision Mhmin_SMcalc,Mhmax_SMcalc

	Mhmin_SMcalc=0.8d0
	Mhmax_SMcalc=500d0

	do i=1,NHIGGS
	BRhjbbbar(i)=0d0
	BRhjtautau(i)=0d0
	BRhjWW(i)=0d0
	BRhjgamgam(i)=0d0
	enddo

! translate g_hjhiZ_sq_compact to  g_hjhiZ_sq
        k=0
        do j=1,NHIGGS
          do i=1,NHIGGS
            if(i .le. j) then
              k=k+1
              g_hjhiZ_sq(j,i)=g_hjhiZ_sq_compact(k)
              if(i .ne. j) g_hjhiZ_sq(i,j)=g_hjhiZ_sq(j,i)
            endif
          enddo
        enddo

	do n=1,NHIGGS,1
         if( (Mh(n) .ge. Mhmin_SMcalc) .and. (Mh(n) .le. Mhmax_SMcalc) ) then 
! BRs
	  if(Gamma_h(n) .eq. 0d0) then
! if Gamma_h(n)=0 then no decay is possible.
	    R_Gamma=0d0
	  else
	    R_Gamma=SMGamma_h(Mh(n))/Gamma_h(n)
	  endif
	   BRhjbbbar(n)=g_hjbb_sq(n)*SMBR_Hbb(Mh(n))*R_Gamma
	   BRhjtautau(n)=g_hjtautau_sq(n)*SMBR_Htautau(Mh(n))*R_Gamma
	   BRhjWW(n)=g_hjWW_sq(n)*SMBR_HWW(Mh(n))*R_Gamma
	   BRhjgamgam(n)=g_hjgamgam_sq(n)*SMBR_Hgamgam(Mh(n))*R_Gamma
	 endif
! CS ratios
	
	  R_hjZ_ee(n)=g_hjZZ_sq(n)

	  do i=1,NHIGGS,1
	    R_hjhi_ee(n,i)=g_hjhiZ_sq(n,i)
	  enddo
	
	  R_hj_gg(n)=g_hjgg_sq(n)
          R_hj_bb(n)=g_hjbb_sq(n)

	  R_hjb_bg(n)=g_hjbb_sq(n)
	
          R_hjZ_ddb(n)=R_hjZ_ee(n)
	  R_hjZ_uub(n)=R_hjZ_ee(n)
          R_hjZ_ssb(n)=R_hjZ_ee(n)
          R_hjZ_ccb(n)=R_hjZ_ee(n)
          R_hjZ_bbb(n)=R_hjZ_ee(n)

	  R_hjWp_udb(n)=g_hjWW_sq(n)
	  R_hjWp_csb(n)=R_hjWp_udb(n)
	  R_hjWm_dub(n)=R_hjWp_udb(n)
	  R_hjWm_scb(n)=R_hjWp_udb(n)

c relative contributuion of WW- and ZZ-fusion to VBF (in LO) for 
c p p-bar collisions at SqrtS=1.96 TeV (calcuated by Terrance with VBF@NLO):
	  SM_ZZ_contrib_to_VBF=0.23d0
	  SM_WW_contrib_to_VBF=0.77d0
	  R_sigma_H_VBF(n)= g_hjWW_sq(n)*SM_WW_contrib_to_VBF
     &		      +g_hjZZ_sq(n)*SM_ZZ_contrib_to_VBF
	
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

	end

!***********************************************************
	subroutine calc_hadCS_from_partCS(NHIGGS, Mh, 
     &  	R_hj_gg, R_hj_bb, R_hjb_bg,
     &          R_hjWp_udb, R_hjWp_csb, R_hjWm_dub, R_hjWm_scb,
     &          R_hjZ_ddb, R_hjZ_uub, R_hjZ_ssb, R_hjZ_ccb, R_hjZ_bbb,
     &		R_sigma_H, R_sigma_Hb, R_sigma_HW, R_sigma_HZ)
!***********************************************************
! Input: 
!	Mh(NHMAX) : the NHMAX neutral Higgs masses
!    ratio of partonic cross sections (model/SM) for the process:
!	R_hj_gg	   : g g     -> H_j
!	R_hj_bb    : b b-bar -> H_j
!	R_hjb_bg   : g b/b-bar -> H_j b/b-bar
!	R_hjWp_udb : u d-bar -> H_j W+
!	R_hjWp_csb : c s-bar -> H_j W+
!	R_hjWm_dub : d u-bar -> H_j W-
!	R_hjWm_scb : s c-bar -> H_j W-
!	R_hjZ_ddb  : d d-bar -> H_j Z
!	R_hjZ_uub  : u u-bar -> H_j Z
!	R_hjZ_ssb  : s s-bar -> H_j Z
!	R_hjZ_ccb  : c c-bar -> H_j Z
!	R_hjZ_bbb  : b b-bar -> H_j Z
!--------------------------------
! Output:
!    ratio of hadronic cross section for the processes:
!	R_sigma_H(NHMAX) 	: p p-bar -> H + X
!	R_sigma_Hb(NHMAX) 	: p p-bar -> H b/b-bar + X
!	R_sigma_HW(NHMAX) 	: p p-bar -> H W+/- + X
!	R_sigma_HZ(NHMAX) 	: p p-bar -> H Z + X
!--------------------------------
	implicit none 
	integer NHIGGS,i,n
! ratios of partonic cross sections
	double precision R_hj_gg(NHMAX),R_hj_bb(NHMAX)
        double precision R_hjb_bg(NHMAX)
	double precision R_hjWp_udb(NHMAX),R_hjWp_csb(NHMAX),
     &			 R_hjWm_dub(NHMAX),R_hjWm_scb(NHMAX)
        double precision R_hjZ_ddb(NHMAX),R_hjZ_uub(NHMAX),
     &                   R_hjZ_ssb(NHMAX),R_hjZ_ccb(NHMAX),R_hjZ_bbb(NHMAX)
	double precision most_sens_predratio
	double precision obsratio
! ratios of hadronic cross sections (new model/SM)
	double precision R_sigma_HW(NHMAX),R_sigma_HZ(NHMAX),R_sigma_H(NHMAX),
     &			 R_sigma_Hb(NHMAX)
! fit functions for SM partial hadronic cross section ratios
	double precision tev_rH_bb
	double precision tev_rH_gg
	double precision tev_rHWpm_dub
	double precision tev_rHWpm_scb
	double precision tev_rHWpm_csb
	double precision tev_rHWpm_udb
	double precision tev_rHZ_bbb
	double precision tev_rHZ_ccb
	double precision tev_rHZ_ddb
	double precision tev_rHZ_ssb
	double precision tev_rHZ_uub

	integer sim_mass_sum(NHMAX),similar_mass(NHMAX,NHMAX)
        double precision Mh(NHMAX),Mh_avg(NHMAX),delta_Mh
	double precision Mhmin_SMcalc,Mhmax_SMcalc

#include "tables-tev_cb.h"
	Mhmin_SMcalc=65d0
	Mhmax_SMcalc=335d0

	do i=1,NHIGGS
	R_sigma_HW(i)=0d0
	R_sigma_HZ(i)=0d0
	R_sigma_H(i)=0d0
	R_sigma_Hb(i)=0d0
	enddo

	  do n=1,NHIGGS

	   if( (Mh(n) .ge. Mhmin_SMcalc) .and. (Mh(n) .le. Mhmax_SMcalc) ) then

	    R_sigma_HW(n)=R_hjWp_udb(n)*tev_rHWpm_udb(Mh(n)) 
     &		+R_hjWp_csb(n)*tev_rHWpm_csb(Mh(n)) 
     &		+R_hjWm_dub(n)*tev_rHWpm_dub(Mh(n)) 
     &		+R_hjWm_scb(n)*tev_rHWpm_scb(Mh(n))
	    R_sigma_HZ(n)=R_hjZ_ddb(n)*tev_rHZ_ddb(Mh(n)) 
     &		+R_hjZ_uub(n)*tev_rHZ_uub(Mh(n)) 
     &		+R_hjZ_ssb(n)*tev_rHZ_ssb(Mh(n)) 
     &		+R_hjZ_ccb(n)*tev_rHZ_ccb(Mh(n)) 
     &		+R_hjZ_bbb(n)*tev_rHZ_bbb(Mh(n))
	    R_sigma_H(n)=R_hj_gg(n)*tev_rH_gg(Mh(n)) 
     &		+R_hj_bb(n)*tev_rH_bb(Mh(n))
	    R_sigma_Hb(n)=R_hjb_bg(n)

	    endif

	   enddo

	end


!******************************************************************
	subroutine check_SM_likeness(table_applicable, tablenum, Mh,
     &		R_sigma_HW,R_sigma_HZ,R_sigma_H,R_sigma_H_VBF,
     &		R_BRhjbbbar,R_BRhjtautau,R_BRhjWW,R_BRhjgamgam)
!******************************************************************
	implicit none
	logical table_applicable
	integer max_n_prod_channels,max_n_decays
	integer n_prod_channels,n_decays,n,m,tablenum,i

	parameter(max_n_prod_channels=4)
	parameter(max_n_decays=4)

	double precision my_maxval

	double precision R_sigma(max_n_prod_channels),R_BR(max_n_decays)
        double precision rel_dR_sigma(max_n_prod_channels),rel_dR_BR(max_n_decays)
	double precision rel_delta_R_sig_BR(max_n_decays*max_n_prod_channels)
	double precision R_sigma_mean,R_BR_mean,eps
	
	double precision Mh,R_sigma_HW,R_sigma_HZ,R_sigma_H,R_sigma_H_VBF
        double precision R_BRhjbbbar,R_BRhjWW,R_BRhjtautau,R_BRhjgamgam

! The models is exactly SM-like, if for all x and y the ratios 
!	sigma(x)_model * BR_model(y)
!      ------------------------------ = R_sigma(x) * R_BR(y) 
!         sigma(x)_SM * BR_SM(y)
! are equal to the same constant.

	eps=EPS_SM_LIKENESS
! eps: tolerable maximum of individual deviations of R_sigma(x)*R_BR(y) 
!      from the mean. Eventually, this should be done Mh-dependent.
!      Suggested default value: 0.02
!      Value set in tables-nhmax.F .

! All tables, except tables 
! CDF9891, D05586, D00808_1266, CDF0802_0432, CDF8961_D05536, CDF9248, CDF9290_D05645, CDF9465_D05754, 
! CDF9674, D05858, CDF9887, D05871, D00901_1887, D00712_0598, CDF9713_D05889, CDF9897
! describe single search channel limits.
	if((tablenum .ne. CDF9891) .and. (tablenum .ne. CDF8961_D05536) 
     &		.and. (tablenum .ne. D05586) .and. (tablenum .ne. D00808_1266) 
     &		.and. (tablenum .ne. CDF0802_0432) .and. (tablenum .ne. CDF9248) 
     &		.and. (tablenum .ne. CDF9290_D05645) .and. (tablenum .ne. CDF9465_D05754) 
     &		.and. (tablenum .ne. CDF9674) .and. (tablenum .ne. D05858)
     & 		.and. (tablenum .ne. CDF9887) .and. (tablenum .ne. D05871)
     &		.and. (tablenum .ne. D00901_1887) .and. (tablenum .ne. D00712_0598)
     &          .and. (tablenum .ne. CDF9713_D05889) .and. (tablenum .ne. CDF9897) ) then
		table_applicable = .TRUE.
	else

#ifdef DEBUG
	write(*,*)'SM likeness check, Mh=',Mh
#endif

c>> 17.4.2009 probably not needed
c	 R_sigma(1)=R_sigma_HW
c         R_sigma(2)=R_sigma_HZ
c         R_sigma(3)=R_sigma_H
c	 R_sigma(4)=R_sigma_H_VBF
c 
c 	 R_BR(1)=R_BRhjbbbar
c         R_BR(2)=R_BRhjWW
c	 R_BR(3)=R_BRhjtautau
c         R_BR(4)=R_BRhjgamgam
c<<

	 do n=1,max_n_prod_channels,1
	 do m=1,max_n_decays,1
	 rel_delta_R_sig_BR(m+(n-1)*max_n_decays)=0d0
	 enddo
	 enddo

!**********
         if((tablenum .eq. CDF8961_D05536) .or. (tablenum .eq. D00712_0598)) then
	 
	  n_prod_channels=3
	  n_decays=2

	 R_sigma(1)=R_sigma_HW
         R_sigma(2)=R_sigma_HZ
         R_sigma(3)=R_sigma_H
 
 	 R_BR(1)=R_BRhjbbbar
         R_BR(2)=R_BRhjWW


	  R_sigma_mean=0d0
	  do n=1,n_prod_channels
	     R_sigma_mean=R_sigma_mean+R_sigma(n)/dble(n_prod_channels)
	  enddo
	  R_BR_mean=0d0
	  do n=1,n_decays
             R_BR_mean=R_BR_mean+R_BR(n)/dble(n_decays)
          enddo
	 
	  do n=1,n_prod_channels
	    rel_dR_sigma(n)=(R_sigma(n)-R_sigma_mean)/R_sigma_mean
	  enddo
          do n=1,n_decays
	    rel_dR_BR(n)=(R_BR(n)-R_BR_mean)/R_BR_mean
          enddo
 
          do n=1,n_prod_channels
	  do m=1,n_decays
	    rel_delta_R_sig_BR(m+(n-1)*n_decays)= 
     &	  	 rel_dR_sigma(n) 
     &	 	+rel_dR_BR(m) 
     &	 	+rel_dR_sigma(n)*rel_dR_BR(m)
	    rel_delta_R_sig_BR(m+(n-1)*n_decays)=dabs(rel_delta_R_sig_BR(m+(n-1)*n_decays))
          enddo
          enddo
	 
	  
          if(my_maxval(rel_delta_R_sig_BR,max_n_decays*max_n_prod_channels) 
     &	  .gt. eps) then
	 	table_applicable=.FALSE.
#ifdef DEBUG
 	 	write(*,*)'table ',tablenum,' does not apply to model point'
#endif
 	  else
 	 	table_applicable=.TRUE.
#ifdef DEBUG
 	 	write(*,*)'table ',tablenum,' applies to model point'
#endif
 	  endif


!**********
	 elseif((tablenum .eq. CDF9891) .or. (tablenum .eq. D05586) 
     & 		.or. (tablenum .eq. D00808_1266) .or. (tablenum .eq. CDF0802_0432)) then

! Combination of only 2 production channels: HW and HZ production 
! and one decay channel: H -> bb	

! Thus:
	  R_sigma_mean=(R_sigma_HW+R_sigma_HZ)/2d0
! Distance from mean is the same for both c.s. ratios.
	  if(R_sigma_mean .eq. 0d0) then
		table_applicable=.FALSE.
!		write(*,*)'table ',tablenum,' does not apply to model point'
	  else
           rel_dR_sigma(1)=(R_sigma_HW-R_sigma_mean)/R_sigma_mean	
! Single channel decay, so rel_dR_BR = 0.
	   rel_delta_R_sig_BR(1)=dabs(rel_dR_sigma(1))

           if(rel_delta_R_sig_BR(1) .gt. eps) then
	 	table_applicable=.FALSE.
#ifdef DEBUG
 	 	write(*,*)'table ',tablenum,' does not apply to model point'
#endif
 	   else
 	 	table_applicable=.TRUE.
#ifdef DEBUG
 	 	write(*,*)'table ',tablenum,' applies to model point'
#endif
 	   endif
	  endif


!**********
         elseif((tablenum .eq. CDF9248) .or. (tablenum .eq. D05858)
     &       	.or. (tablenum .eq. CDF9887) .or. (tablenum .eq. D05871)
     & 		.or. (tablenum .eq. D00901_1887)) then
	 
	  n_prod_channels=4
	  n_decays=1

	 R_sigma(1)=R_sigma_HW
         R_sigma(2)=R_sigma_HZ
         R_sigma(3)=R_sigma_H
	 R_sigma(4)=R_sigma_H_VBF
 
	  R_sigma_mean=0d0
	  do n=1,n_prod_channels
	     R_sigma_mean=R_sigma_mean+R_sigma(n)/dble(n_prod_channels)
	  enddo

! Definition of R_BR_mean not needed. Check applies to all single decay mode
! analyses with thes above combination of production channels, 
! e.g. H -> tau tau and H -> gamma gamma.
	 
	  do n=1,n_prod_channels
	    rel_dR_sigma(n)=(R_sigma(n)-R_sigma_mean)/R_sigma_mean
	  enddo
! Single channel decay, so rel_dR_BR = 0.
	    rel_dR_BR(1)=0d0

          do n=1,n_prod_channels
	  do m=1,n_decays
	    rel_delta_R_sig_BR(m+(n-1)*n_decays)= 
     &	  	 rel_dR_sigma(n) 
     &	 	+rel_dR_BR(m) 
     &	 	+rel_dR_sigma(n)*rel_dR_BR(m)
          enddo
          enddo
 
          if(my_maxval(rel_delta_R_sig_BR,max_n_decays*max_n_prod_channels) 
     &	  .gt. eps) then
	 	table_applicable=.FALSE.
#ifdef DEBUG
 	 	write(*,*)'table ',tablenum,' does not apply to model point'
#endif
 	  else
 	 	table_applicable=.TRUE.
#ifdef DEBUG
 	 	write(*,*)'table ',tablenum,' applies to model point'
#endif
 	  endif


!**********
         elseif(tablenum .eq. CDF9290_D05645) then
	 
	  n_prod_channels=4
	  n_decays=4

	 R_sigma(1)=R_sigma_HW
         R_sigma(2)=R_sigma_HZ
         R_sigma(3)=R_sigma_H
	 R_sigma(4)=R_sigma_H_VBF
 
 	 R_BR(1)=R_BRhjbbbar
         R_BR(2)=R_BRhjWW
	 R_BR(3)=R_BRhjtautau
         R_BR(4)=R_BRhjgamgam

	  R_sigma_mean=0d0
	  do n=1,n_prod_channels
	     R_sigma_mean=R_sigma_mean+R_sigma(n)/dble(n_prod_channels)
	  enddo
	  R_BR_mean=0d0
	  do n=1,n_decays
             R_BR_mean=R_BR_mean+R_BR(n)/dble(n_decays)
          enddo
	 
	  do n=1,n_prod_channels
	    rel_dR_sigma(n)=(R_sigma(n)-R_sigma_mean)/R_sigma_mean
	  enddo
          do n=1,n_decays
	    rel_dR_BR(n)=(R_BR(n)-R_BR_mean)/R_BR_mean
          enddo
 
          do n=1,n_prod_channels
	  do m=1,n_decays
	    rel_delta_R_sig_BR(m+(n-1)*n_decays)= 
     &	  	 rel_dR_sigma(n) 
     &	 	+rel_dR_BR(m) 
     &	 	+rel_dR_sigma(n)*rel_dR_BR(m)
          enddo
          enddo
 
          if(my_maxval(rel_delta_R_sig_BR,max_n_decays*max_n_prod_channels) 
     &	  .gt. eps) then
	 	table_applicable=.FALSE.
#ifdef DEBUG
 	 	write(*,*)'table CDF9290_D05645 does not apply to model point'
#endif
 	  else
 	 	table_applicable=.TRUE.
#ifdef DEBUG
 	 	write(*,*)'table CDF9290_D05645 applies to model point'
#endif
 	  endif


!**********
         elseif((tablenum .eq. CDF9674) .or. (tablenum .eq. CDF9897)) then
	 
	  n_prod_channels=4
	  n_decays=3

	 R_sigma(1)=R_sigma_HW
         R_sigma(2)=R_sigma_HZ
         R_sigma(3)=R_sigma_H
	 R_sigma(4)=R_sigma_H_VBF
 
 	 R_BR(1)=R_BRhjbbbar
         R_BR(2)=R_BRhjWW
	 R_BR(3)=R_BRhjtautau

	  R_sigma_mean=0d0
	  do n=1,n_prod_channels
	     R_sigma_mean=R_sigma_mean+R_sigma(n)/dble(n_prod_channels)
	  enddo
	  R_BR_mean=0d0
	  do n=1,n_decays
             R_BR_mean=R_BR_mean+R_BR(n)/dble(n_decays)
          enddo
	 
	  do n=1,n_prod_channels
	    rel_dR_sigma(n)=(R_sigma(n)-R_sigma_mean)/R_sigma_mean
	  enddo
          do n=1,n_decays
	    rel_dR_BR(n)=(R_BR(n)-R_BR_mean)/R_BR_mean
          enddo
 
          do n=1,n_prod_channels
	  do m=1,n_decays
	    rel_delta_R_sig_BR(m+(n-1)*n_decays)= 
     &	  	 rel_dR_sigma(n) 
     &	 	+rel_dR_BR(m) 
     &	 	+rel_dR_sigma(n)*rel_dR_BR(m)
          enddo
          enddo
 
          if(my_maxval(rel_delta_R_sig_BR,max_n_decays*max_n_prod_channels) 
     &	  .gt. eps) then
	 	table_applicable=.FALSE.
#ifdef DEBUG
 	 	write(*,*)'table ',tablenum,' does not apply to model point'
#endif
 	  else
 	 	table_applicable=.TRUE.
#ifdef DEBUG
 	 	write(*,*)'table ',tablenum,' applies to model point'
#endif
 	  endif


!**********
         elseif(tablenum .eq. CDF9465_D05754) then
	 
	  n_prod_channels=4
	  n_decays=1

	  R_sigma(1)=R_sigma_H
          R_sigma(2)=R_sigma_HW
	  R_sigma(3)=R_sigma_HZ
          R_sigma(4)=R_sigma_H_VBF
 
          R_BR(1)=R_BRhjWW

	  R_sigma_mean=0d0
	  do n=1,n_prod_channels
	     R_sigma_mean=R_sigma_mean+R_sigma(n)/dble(n_prod_channels)
	  enddo
	  R_BR_mean=R_BR(1)
	  rel_dR_BR(1)=0d0
	 
	  do n=1,n_prod_channels
	    rel_dR_sigma(n)=(R_sigma(n)-R_sigma_mean)/R_sigma_mean
	  enddo

 
          do n=1,n_prod_channels
	  do m=1,n_decays
	    rel_delta_R_sig_BR(m+(n-1)*n_decays)= 
     &	  	 rel_dR_sigma(n) 
     &	 	+rel_dR_BR(m) 
     &	 	+rel_dR_sigma(n)*rel_dR_BR(m)
	    rel_delta_R_sig_BR(m+(n-1)*n_decays)=dabs(rel_delta_R_sig_BR(m+(n-1)*n_decays))
          enddo
          enddo
	 
	  
          if(my_maxval(rel_delta_R_sig_BR,max_n_decays*max_n_prod_channels) 
     &	  .gt. eps) then
	 	table_applicable=.FALSE.
#ifdef DEBUG
 	 	write(*,*)'table  CDF9465_D05754 does not apply to model point'
#endif
 	  else
 	 	table_applicable=.TRUE.
#ifdef DEBUG
 	 	write(*,*)'table  CDF9465_D05754 applies to model point'
#endif
 	  endif


!**********
         elseif(tablenum .eq. CDF9713_D05889) then
! CDF9713_D05889 : 
! p p-bar 
! -> HW -> 	b b l nu, 					not yet needed
!	   	b b-bar tau nu_tau, 				not yet needed
!		W+ W- W -> l+- nu l+- nu			x
! -> HV -> 	j j tau tau/tau tau j j (for m_H <= 145 GeV), 	not yet needed
!		b b + miss. E_T,				not yet needed
! 		b b j j 	(V=W,Z)				not yet needed
! -> HZ -> b b l l, 						not yet needed
! -> H via single H/HW/HZ/VBF + 0/1/2 jets, H -> W+ W- -> l+ nu l- nu,		x
! -> Htt -> b b t t (for m_H <= 155 GeV)			not yet needed
! -> H via single H/WH/ZH/VBF, 
!	H -> tau tau, 						not yet needed
!	H -> gamma gamma					not yet needed
	 
	  n_prod_channels=4
	  n_decays=1

	 R_sigma(1)=R_sigma_HW
         R_sigma(2)=R_sigma_HZ
         R_sigma(3)=R_sigma_H
	 R_sigma(4)=R_sigma_H_VBF
 
 	 R_BR(1)=R_BRhjWW

	  R_sigma_mean=0d0
	  do n=1,n_prod_channels
	     R_sigma_mean=R_sigma_mean+R_sigma(n)/dble(n_prod_channels)
	  enddo

	  R_BR_mean=R_BR(1)

	 
	  do n=1,n_prod_channels
	    rel_dR_sigma(n)=(R_sigma(n)-R_sigma_mean)/R_sigma_mean
	  enddo
	  rel_dR_BR(1)=0d0

 
          do n=1,n_prod_channels
	  do m=1,n_decays
	    rel_delta_R_sig_BR(m+(n-1)*n_decays)= 
     &	  	 rel_dR_sigma(n) 
     &	 	+rel_dR_BR(m) 
     &	 	+rel_dR_sigma(n)*rel_dR_BR(m)
          enddo
          enddo
 
          if(my_maxval(rel_delta_R_sig_BR,max_n_decays*max_n_prod_channels) 
     &	  .gt. eps) then
	 	table_applicable=.FALSE.
#ifdef DEBUG
 	 	write(*,*)'table CDF9713_D05889 does not apply to model point'
#endif
 	  else
 	 	table_applicable=.TRUE.
#ifdef DEBUG
 	 	write(*,*)'table CDF9713_D05889 applies to model point'
#endif
 	  endif


	 endif

	endif
	

	end



!******************************************************************
	integer function TEV_reaction_number(NHIGGS,channel)
!******************************************************************
	implicit none
	integer NHIGGS,channel

	if(channel .le. TEV_NSCHN*NHIGGS) then
	  TEV_reaction_number=channel-int((channel-1)/NHIGGS)*NHIGGS
	else
	  TEV_reaction_number=channel-TEV_NSCHN*NHIGGS-int((channel-TEV_NSCHN*NHIGGS-1)/NHIGGS**2)*NHIGGS**2
	endif

	end



!******************************************************************
	subroutine get_1d_TEVtable_data(selec,tabsize,Mhmin,Mhmax,sep)
!******************************************************************
	implicit none
	integer tabsize,selec
        double precision sep
	double precision Mhmin,Mhmax

#include "tables-tev_cb.h"

	if(selec .eq. CDF0908_3534) then
	 Mhmin=Mhmin_CDF_ZH_llbb_pub_2
	 Mhmax=Mhmax_CDF_ZH_llbb_pub_2
         sep=sep_CDF_ZH_llbb_pub_2
	elseif(selec .eq. CDF9889) then
         Mhmin=Mhmin_CDF_ZH_llbb_4
         Mhmax=Mhmax_CDF_ZH_llbb_4
         sep=sep_CDF_ZH_llbb_4
        elseif(selec .eq. CDF9891) then
         Mhmin=Mhmin_CDF_VH_bbmET_4
         Mhmax=Mhmax_CDF_VH_bbmET_4
         sep=sep_CDF_VH_bbmET_4
        elseif(selec .eq. CDF0802_0432) then
         Mhmin=Mhmin_CDF_VH_bbmET_pub
         Mhmax=Mhmax_CDF_VH_bbmET_pub
         sep=sep_CDF_VH_bbmET_pub
        elseif(selec .eq. D05876) then
         Mhmin=Mhmin_D0_ZH_llbb_3
         Mhmax=Mhmax_D0_ZH_llbb_3
         sep=sep_D0_ZH_llbb_3
        elseif(selec .eq. D05586) then
         Mhmin=Mhmin_D0_VH_bbmET
         Mhmax=Mhmax_D0_VH_bbmET
         sep=sep_D0_VH_bbmET
        elseif(selec .eq. D00808_1266) then
         Mhmin=Mhmin_D0_VH_bbmET_pub
         Mhmax=Mhmax_D0_VH_bbmET_pub
         sep=sep_D0_VH_bbmET_pub
        elseif(selec .eq. CDF8961_D05536) then
         Mhmin=Mhmin_CDF_D0_SM
         Mhmax=Mhmax_CDF_D0_SM
         sep=sep_CDF_D0_SM
        elseif(selec .eq. CDF9284) then
         Mhmin=Mhmin_CDF_bH_bb
         Mhmax=Mhmax_CDF_bH_bb
         sep=sep_CDF_bH_bb
        elseif(selec .eq. D05972) then
         Mhmin=Mhmin_D0_WH_lnubb_3
         Mhmax=Mhmax_D0_WH_lnubb_3
         sep=sep_D0_WH_lnubb_3
	elseif(selec .eq. CDFwhlnubb_090814) then
	 Mhmin=Mhmin_CDF_WH_lnubb_5
	 Mhmax=Mhmax_CDF_WH_lnubb_5
         sep=sep_CDF_WH_lnubb_5
        elseif(selec .eq. D05726) then
         Mhmin=Mhmin_D0_bH_bb_tb_limit_2
         Mhmax=Mhmax_D0_bH_bb_tb_limit_2
         sep=sep_D0_bH_bb_tb_limit_2
        elseif(selec .eq. D00805_3556) then
         Mhmin=Mhmin_D0_bH_bb_2
         Mhmax=Mhmax_D0_bH_bb_2
         sep=sep_D0_bH_bb_2
        elseif(selec .eq. D00805_2491) then
         Mhmin=Mhmin_D0_ppH_tautau_pub
         Mhmax=Mhmax_D0_ppH_tautau_pub
         sep=sep_D0_ppH_tautau_pub
        elseif(selec .eq. CDF0906_1014) then
         Mhmin=Mhmin_CDF_ppH_tautau_pub
         Mhmax=Mhmax_CDF_ppH_tautau_pub
         sep=sep_CDF_ppH_tautau_pub
        elseif(selec .eq. D05873) then
         Mhmin=Mhmin_D0_WH_WWW_2
         Mhmax=Mhmax_D0_WH_WWW_2
         sep=sep_D0_WH_WWW_2
        elseif(selec .eq. CDF7307v3) then
         Mhmin=Mhmin_CDF_WH_WWW_2
         Mhmax=Mhmax_CDF_WH_WWW_2
         sep=sep_CDF_WH_WWW_2
        elseif(selec .eq. D05858) then
         Mhmin=Mhmin_D0_ppH_gamgam_4
         Mhmax=Mhmax_D0_ppH_gamgam_4
         sep=sep_D0_ppH_gamgam_4
!temporary removed         elseif(selec .eq. D00803_1514) then
!temporary removed          Mhmin=Mhmin_D0_ppH_gamgam_pub
!temporary removed          Mhmax=Mhmax_D0_ppH_gamgam_pub
!temporary removed          sep=sep_D0_ppH_gamgam_pub
        elseif(selec .eq. D00901_1887) then
         Mhmin=Mhmin_D0_ppH_gamgam_pub_2
         Mhmax=Mhmax_D0_ppH_gamgam_pub_2
         sep=sep_D0_ppH_gamgam_pub_2
        elseif(selec .eq. CDF9248) then
         Mhmin=Mhmin_CDF_SM_tautau
         Mhmax=Mhmax_CDF_SM_tautau
         sep=sep_CDF_SM_tautau
        elseif(selec .eq. CDF9290_D05645) then
         Mhmin=Mhmin_CDF_D0_SM_2
         Mhmax=Mhmax_CDF_D0_SM_2
         sep=sep_CDF_D0_SM_2
        elseif(selec .eq. D05757) then
         Mhmin=Mhmin_D0_ppH_WW_ll_2
         Mhmax=Mhmax_D0_ppH_WW_ll_2
         sep=sep_D0_ppH_WW_ll_2
        elseif(selec .eq. CDF0809_3930) then
         Mhmin=Mhmin_CDF_ppH_WW_4
         Mhmax=Mhmax_CDF_ppH_WW_4
         sep=sep_CDF_ppH_WW_4
        elseif(selec .eq. CDF9465_D05754) then
         Mhmin=Mhmin_CDF_D0_SM_3
         Mhmax=Mhmax_CDF_D0_SM_3
         sep=sep_CDF_D0_SM_3
        elseif(selec .eq. CDF9674) then
         Mhmin=Mhmin_CDF_SM
         Mhmax=Mhmax_CDF_SM
         sep=sep_CDF_SM
        elseif(selec .eq. CDF9897) then
         Mhmin=Mhmin_CDF_SM_2
         Mhmax=Mhmax_CDF_SM_2
         sep=sep_CDF_SM_2
        elseif(selec .eq. D00712_0598) then
         Mhmin=Mhmin_D0_SM_pub
         Mhmax=Mhmax_D0_SM_pub
         sep=sep_D0_SM_pub
        elseif(selec .eq. CDF9887) then
         Mhmin=Mhmin_CDF_ppH_WW_6
         Mhmax=Mhmax_CDF_ppH_WW_6
         sep=sep_CDF_ppH_WW_6
        elseif(selec .eq. D05871) then
         Mhmin=Mhmin_D0_ppH_WW_ll_3
         Mhmax=Mhmax_D0_ppH_WW_ll_3
         sep=sep_D0_ppH_WW_ll_3
        elseif(selec .eq. D00808_1970) then
         Mhmin=Mhmin_D0_WH_lnubb_pub
         Mhmax=Mhmax_D0_WH_lnubb_pub
         sep=sep_D0_WH_lnubb_pub
        elseif(selec .eq. CDF0906_5613) then
         Mhmin=Mhmin_CDF_WH_lnubb_pub_2
         Mhmax=Mhmax_CDF_WH_lnubb_pub_2
         sep=sep_CDF_WH_lnubb_pub_2
        elseif(selec .eq. D00811_0024) then
         Mhmin=Mhmin_D0_bH_tautau_pub
         Mhmax=Mhmax_D0_bH_tautau_pub
         sep=sep_D0_bH_tautau_pub
        elseif(selec .eq. D05985) then
         Mhmin=Mhmin_D0_bH_tautau
         Mhmax=Mhmax_D0_bH_tautau
         sep=sep_D0_bH_tautau
        elseif(selec .eq. D05740) then
         Mhmin=Mhmin_D0_ppH_tautau
         Mhmax=Mhmax_D0_ppH_tautau
         sep=sep_D0_ppH_tautau
        elseif(selec .eq. CDF9888_D05980) then
         Mhmin=Mhmin_CDF_D0_ppH_tautau
         Mhmax=Mhmax_CDF_D0_ppH_tautau
         sep=sep_CDF_D0_ppH_tautau
        elseif(selec .eq. CDF9713_D05889) then
         Mhmin=Mhmin_CDF_D0_SM_4
         Mhmax=Mhmax_CDF_D0_SM_4
         sep=sep_CDF_D0_SM_4
	endif

	tabsize=nint((Mhmax-Mhmin)/sep)+1

	end


!******************************************************************
	subroutine check_TEV_bound1(Mh1, 
     &	  R_sig_BR,ratio,selec,tabsize,Mhmin,Mhmax,sep,pred)
!******************************************************************
	implicit none
	logical pred
	integer i,ilow,ihigh,selec,tabsize
	double precision Mh1,Mh1_test,Mh1_below,R_sig_BR
        double precision sep
	double precision Mhmin,Mhmax,ratio
	double precision Mh1dummy,Mh2dummy,S95obs,S95pred	
        double precision S95_above,S95_below,S95_interpol
        double precision S95(tabsize)

#include "tables-tev_cb.h"


	if(selec .eq. CDF0908_3534) then
	  if(pred) then
		call eval_1d_arr(1,tabsize,S95,CDF_ZH_llbb_pub_2_pred)
	  else	
		call eval_1d_arr(1,tabsize,S95,CDF_ZH_llbb_pub_2_obs)
	  endif
	elseif(selec .eq. CDF9889) then
          if(pred) then
		call eval_1d_arr(1,tabsize,S95,CDF_ZH_llbb_4_pred)
          else
		call eval_1d_arr(1,tabsize,S95,CDF_ZH_llbb_4_obs)
          endif
        elseif(selec .eq. CDF9891) then
          if(pred) then
		call eval_1d_arr(1,tabsize,S95,CDF_VH_bbmET_4_pred)
          else
		call eval_1d_arr(1,tabsize,S95,CDF_VH_bbmET_4_obs)
          endif
        elseif(selec .eq. CDF0802_0432) then
          if(pred) then
		call eval_1d_arr(1,tabsize,S95,CDF_VH_bbmET_pub_pred)
          else
		call eval_1d_arr(1,tabsize,S95,CDF_VH_bbmET_pub_obs)
          endif
        elseif(selec .eq. D05876) then
          if(pred) then
		call eval_1d_arr(1,tabsize,S95,D0_ZH_llbb_3_pred)
          else
		call eval_1d_arr(1,tabsize,S95,D0_ZH_llbb_3_obs)
          endif
        elseif(selec .eq. D05586) then
          if(pred) then
		call eval_1d_arr(1,tabsize,S95,D0_VH_bbmET_pred)
          else
		call eval_1d_arr(1,tabsize,S95,D0_VH_bbmET_obs)
          endif
        elseif(selec .eq. D00808_1266) then
          if(pred) then
		call eval_1d_arr(1,tabsize,S95,D0_VH_bbmET_pub_pred)
          else
		call eval_1d_arr(1,tabsize,S95,D0_VH_bbmET_pub_obs)
          endif
        elseif(selec .eq. CDF8961_D05536) then
          if(pred) then
		call eval_1d_arr(1,tabsize,S95,CDF_D0_SM_pred)
          else
		call eval_1d_arr(1,tabsize,S95,CDF_D0_SM_obs)
          endif
        elseif(selec .eq. CDF9284) then
          if(pred) then
		call eval_1d_arr(1,tabsize,S95,CDF_bH_bb_pred)
          else
		call eval_1d_arr(1,tabsize,S95,CDF_bH_bb_obs)
          endif
        elseif(selec .eq. D05972) then
          if(pred) then
		call eval_1d_arr(1,tabsize,S95,D0_WH_lnubb_3_pred)
          else
		call eval_1d_arr(1,tabsize,S95,D0_WH_lnubb_3_obs)
          endif
	elseif(selec .eq. CDFwhlnubb_090814) then
	  if(pred) then
		call eval_1d_arr(1,tabsize,S95,CDF_WH_lnubb_5_pred)
	  else	
		call eval_1d_arr(1,tabsize,S95,CDF_WH_lnubb_5_obs)
	  endif
        elseif(selec .eq. D05726) then
          if(pred) then
		call eval_1d_arr(1,tabsize,S95,D0_bH_bb_tb_limit_2_pred)
          else
		call eval_1d_arr(1,tabsize,S95,D0_bH_bb_tb_limit_2_obs)
          endif
        elseif(selec .eq. D00805_3556) then
          if(pred) then
		call eval_1d_arr(1,tabsize,S95,D0_bH_bb_2_pred)
          else
		call eval_1d_arr(1,tabsize,S95,D0_bH_bb_2_obs)
          endif
        elseif(selec .eq. D00805_2491) then
          if(pred) then
		call eval_1d_arr(1,tabsize,S95,D0_ppH_tautau_pub_pred)
          else
		call eval_1d_arr(1,tabsize,S95,D0_ppH_tautau_pub_obs)
          endif
        elseif(selec .eq. CDF0906_1014) then
          if(pred) then
		call eval_1d_arr(1,tabsize,S95,CDF_ppH_tautau_pub_pred)
          else
		call eval_1d_arr(1,tabsize,S95,CDF_ppH_tautau_pub_obs)
          endif
        elseif(selec .eq. D05873) then
          if(pred) then
		call eval_1d_arr(1,tabsize,S95,D0_WH_WWW_2_pred)
          else
		call eval_1d_arr(1,tabsize,S95,D0_WH_WWW_2_obs)
          endif
        elseif(selec .eq. CDF7307v3) then
          if(pred) then
		call eval_1d_arr(1,tabsize,S95,CDF_WH_WWW_2_pred)
          else
		call eval_1d_arr(1,tabsize,S95,CDF_WH_WWW_2_obs)
          endif
        elseif(selec .eq. D05858) then
          if(pred) then
		call eval_1d_arr(1,tabsize,S95,D0_ppH_gamgam_4_pred)
          else
		call eval_1d_arr(1,tabsize,S95,D0_ppH_gamgam_4_obs)
          endif
!temporary removed         elseif(selec .eq. D00803_1514) then
!temporary removed           if(pred) then
!temporary removed 		call eval_1d_arr(1,tabsize,S95,D0_ppH_gamgam_pub_pred)
!temporary removed           else
!temporary removed 		call eval_1d_arr(1,tabsize,S95,D0_ppH_gamgam_pub_obs)
!temporary removed           endif
        elseif(selec .eq. D00901_1887) then
          if(pred) then
		call eval_1d_arr(1,tabsize,S95,D0_ppH_gamgam_pub_2_pred)
          else
		call eval_1d_arr(1,tabsize,S95,D0_ppH_gamgam_pub_2_obs)
          endif
        elseif(selec .eq. CDF9248) then
          if(pred) then
		call eval_1d_arr(1,tabsize,S95,CDF_SM_tautau_pred)
          else
		call eval_1d_arr(1,tabsize,S95,CDF_SM_tautau_obs)
          endif
        elseif(selec .eq. CDF9290_D05645) then
          if(pred) then
		call eval_1d_arr(1,tabsize,S95,CDF_D0_SM_2_pred)
          else
		call eval_1d_arr(1,tabsize,S95,CDF_D0_SM_2_obs)
          endif
        elseif(selec .eq. D05757) then
          if(pred) then
		call eval_1d_arr(1,tabsize,S95,D0_ppH_WW_ll_2_pred)
          else
		call eval_1d_arr(1,tabsize,S95,D0_ppH_WW_ll_2_obs)
          endif
        elseif(selec .eq. CDF0809_3930) then
          if(pred) then
		call eval_1d_arr(1,tabsize,S95,CDF_ppH_WW_4_pred)
          else
		call eval_1d_arr(1,tabsize,S95,CDF_ppH_WW_4_obs)
          endif
        elseif(selec .eq. CDF9465_D05754) then
          if(pred) then
		call eval_1d_arr(1,tabsize,S95,CDF_D0_SM_3_pred)
          else
		call eval_1d_arr(1,tabsize,S95,CDF_D0_SM_3_obs)
          endif
        elseif(selec .eq. CDF9674) then
          if(pred) then
		call eval_1d_arr(1,tabsize,S95,CDF_SM_pred)
          else
		call eval_1d_arr(1,tabsize,S95,CDF_SM_obs)
          endif
        elseif(selec .eq. CDF9897) then
          if(pred) then
		call eval_1d_arr(1,tabsize,S95,CDF_SM_2_pred)
          else
		call eval_1d_arr(1,tabsize,S95,CDF_SM_2_obs)
          endif
        elseif(selec .eq. D00712_0598) then
          if(pred) then
		call eval_1d_arr(1,tabsize,S95,D0_SM_pub_pred)
          else
		call eval_1d_arr(1,tabsize,S95,D0_SM_pub_obs)
          endif
        elseif(selec .eq. CDF9887) then
          if(pred) then
		call eval_1d_arr(1,tabsize,S95,CDF_ppH_WW_6_pred)
          else
		call eval_1d_arr(1,tabsize,S95,CDF_ppH_WW_6_obs)
          endif
        elseif(selec .eq. D05871) then
          if(pred) then
		call eval_1d_arr(1,tabsize,S95,D0_ppH_WW_ll_3_pred)
          else
		call eval_1d_arr(1,tabsize,S95,D0_ppH_WW_ll_3_obs)
          endif
        elseif(selec .eq. D00808_1970) then
          if(pred) then
		call eval_1d_arr(1,tabsize,S95,D0_WH_lnubb_pub_pred)
          else
		call eval_1d_arr(1,tabsize,S95,D0_WH_lnubb_pub_obs)
          endif
        elseif(selec .eq. CDF0906_5613) then
          if(pred) then
		call eval_1d_arr(1,tabsize,S95,CDF_WH_lnubb_pub_2_pred)
          else
		call eval_1d_arr(1,tabsize,S95,CDF_WH_lnubb_pub_2_obs)
          endif
        elseif(selec .eq. D00811_0024) then
          if(pred) then
		call eval_1d_arr(1,tabsize,S95,D0_bH_tautau_pub_pred)
          else
		call eval_1d_arr(1,tabsize,S95,D0_bH_tautau_pub_obs)
          endif
        elseif(selec .eq. D05985) then
          if(pred) then
		call eval_1d_arr(1,tabsize,S95,D0_bH_tautau_pred)
          else
		call eval_1d_arr(1,tabsize,S95,D0_bH_tautau_obs)
          endif
        elseif(selec .eq. D05740) then
          if(pred) then
                call eval_1d_arr(1,tabsize,S95,D0_ppH_tautau_pred)
          else
                call eval_1d_arr(1,tabsize,S95,D0_ppH_tautau_obs)
          endif
        elseif(selec .eq. CDF9888_D05980) then
          if(pred) then
                call eval_1d_arr(1,tabsize,S95,CDF_D0_ppH_tautau_pred)
          else
                call eval_1d_arr(1,tabsize,S95,CDF_D0_ppH_tautau_obs)
          endif
        elseif(selec .eq. CDF9713_D05889) then
          if(pred) then
		call eval_1d_arr(1,tabsize,S95,CDF_D0_SM_4_pred)
          else
		call eval_1d_arr(1,tabsize,S95,CDF_D0_SM_4_obs)
          endif

	endif

        ratio=0.0D0
	
! check if mass is within mass range of table:
	if((Mh1 .gt. Mhmax) .or. (Mh1 .lt. Mhmin)) then
         ratio=0D0
	else	
	
! linear interpolate between points in table
	
	ilow=int((Mh1-Mhmin)/sep)+1
! check if interpolation is needed
!	The expression in the followin if-statement should be mathematically zero.
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
        	 	 
! R_sig_BR is the ratio (sigma*BRs)_model/(sigma*BRs)_ref.
! R_sig_BR for the channel with the highest R_sig_BR/S95_interpol_predicted
! has to be smaller then S95_interpol_observed in order to yield an
! allowed scenario. 
! The latter is equivalent to ratio=R_sig_BR/S95_interpol_observed < 1.

	 ratio=R_sig_BR/S95_interpol		
	  
	endif 

	end



