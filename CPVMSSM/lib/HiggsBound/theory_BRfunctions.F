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

!******************************************************
	subroutine initialize_SM_decay_data()
!******************************************************
	implicit none

!        integer lnblnk
!        external lnblnk
	integer i

	character*20 dummy
	character*200 dir
	double precision num
	
	include 'SM_decay_data_cb.h'
	
	open(44,file=INSTDIR//"/SMBRs/br.sm1_HDecay")
	open(45,file=INSTDIR//"/SMBRs/br.sm2_HDecay")
	
! read 1st three lines and ignore them
	do i=1,3,1
	  read(44,*) dummy
	  read(45,*) dummy
	enddo
	
	do i=1,4993
	  read(44,*) hdecay_tab(i,1),hdecay_tab(i,2),hdecay_tab(i,3),
     & 	    hdecay_tab(i,4),hdecay_tab(i,5),hdecay_tab(i,6),hdecay_tab(i,7)
	  read(45,*) num,hdecay_tab(i,8),hdecay_tab(i,9),hdecay_tab(i,10),
     &	    hdecay_tab(i,11),hdecay_tab(i,12),hdecay_tab(i,13)
	enddo
	
!	write(*,*) (hdecay_tab(1,i), i=1,13)
!	write(*,*) (hdecay_tab(4993,i), i=1,13)
	
	close(44)
	close(45)
	
	end
	
!****************************************************** 
       double precision function SMGamma_h(Mh)
!******************************************************
       implicit none
       integer ilow
       double precision Mh,Mhmin,Mhmax,sep,Gamma_interpol
       double precision  mh_below,Gamma_below,Gamma_above
       
       include 'SM_decay_data_cb.h'
       
       Mhmin=hdecay_tab(1,1)
       Mhmax=hdecay_tab(4993,1)
       if((Mh .ge. Mhmin) .and. (Mh .le. Mhmax)) then
         sep=hdecay_tab(2,1)-hdecay_tab(1,1)
       
         ilow=int((Mh-Mhmin)/sep)+1
       
         if( ((Mh-Mhmin)/sep+1d0-dble(ilow)) .lt. 1d-7) then
       	   Gamma_interpol=hdecay_tab(ilow,13)
         else
           Mh_below=dble(ilow-1)*sep+Mhmin
           Gamma_below=hdecay_tab(ilow,13)
           Gamma_above=hdecay_tab(ilow+1,13)
           Gamma_interpol=Gamma_below +(Gamma_above-Gamma_below)/sep*(MH-Mh_below)         
         endif
       else
         write(*,*)'function SMGamma_h(Mh) out of range'
         stop
       endif
       
       SMGamma_h=Gamma_interpol
       end 


!****************************************************** 
       double precision function SMBR_HWW(Mh)
!******************************************************
       implicit none
       integer ilow
       double precision Mh,Mhmin,Mhmax,sep,BRinterpol
       double precision  mh_below,BRbelow,BRabove
       
       include 'SM_decay_data_cb.h'
       
       Mhmin=hdecay_tab(1,1)
       Mhmax=hdecay_tab(4993,1)
       if((Mh .ge. Mhmin) .and. (Mh .le. Mhmax)) then
         sep=hdecay_tab(2,1)-hdecay_tab(1,1)
         
         ilow=int((Mh-Mhmin)/sep)+1
         
         if( ((Mh-Mhmin)/sep+1d0-dble(ilow)) .lt. 1d-7) then
         	 BRinterpol=hdecay_tab(ilow,11)
         else
           Mh_below=dble(ilow-1)*sep+Mhmin
           BRbelow=hdecay_tab(ilow,11)
           BRabove=hdecay_tab(ilow+1,11)
           BRinterpol=BRbelow +(BRabove-BRbelow)/sep*(MH-Mh_below)
         endif
       else
         write(*,*)'function SMBR_HWW(Mh) out of range'
         stop
       endif
                                
       
       SMBR_HWW=BRinterpol
       end 


!******************************************************
       double precision function SMBR_Hbb(Mh)
!******************************************************
       implicit none
       integer ilow
       double precision Mh,Mhmin,Mhmax,sep,BRinterpol
       double precision mh_below,BRbelow,BRabove
       
       include 'SM_decay_data_cb.h'
       
       Mhmin=hdecay_tab(1,1)
       Mhmax=hdecay_tab(4993,1)
       if((Mh .ge. Mhmin) .and. (Mh .le. Mhmax)) then
         sep=hdecay_tab(2,1)-hdecay_tab(1,1)
         
         ilow=int((Mh-Mhmin)/sep)+1
         
         if( ((Mh-Mhmin)/sep+1d0-dble(ilow)) .lt. 1d-7) then
         	 BRinterpol=hdecay_tab(ilow,2)
         else
           Mh_below=dble(ilow-1)*sep+Mhmin
           BRbelow=hdecay_tab(ilow,2)
           BRabove=hdecay_tab(ilow+1,2)
           BRinterpol=BRbelow +(BRabove-BRbelow)/sep*(MH-Mh_below)
         endif
       else
         write(*,*)'function SMBR_Hbb(Mh) out of range'
         stop
       endif
                                 
       
       SMBR_Hbb=BRinterpol
       end 

 
!******************************************************
       double precision function SMBR_Htautau(Mh)
!******************************************************
       implicit none
       integer ilow
       double precision Mh,Mhmin,Mhmax,sep,BRinterpol
       double precision  mh_below,BRbelow,BRabove
       
       include 'SM_decay_data_cb.h'
       
       Mhmin=hdecay_tab(1,1)
       Mhmax=hdecay_tab(4993,1)
       if((Mh .ge. Mhmin) .and. (Mh .le. Mhmax)) then
         sep=hdecay_tab(2,1)-hdecay_tab(1,1)
         
         ilow=int((Mh-Mhmin)/sep)+1
         
         if( ((Mh-Mhmin)/sep+1d0-dble(ilow)) .lt. 1d-7) then
  	   BRinterpol=hdecay_tab(ilow,3)
         else
           Mh_below=dble(ilow-1)*sep+Mhmin
           BRbelow=hdecay_tab(ilow,3)
           BRabove=hdecay_tab(ilow+1,3)
           BRinterpol=BRbelow +(BRabove-BRbelow)/sep*(MH-Mh_below)
         endif
       else
         write(*,*)'function SMBR_Htautau(Mh) out of range'
         stop
       endif
       
       SMBR_Htautau=BRinterpol       
       end 

       
!******************************************************
       double precision function SMBR_Hgamgam(Mh)
!******************************************************
       implicit none
       integer ilow
       double precision Mh,Mhmin,Mhmax,sep,BRinterpol
       double precision  mh_below,BRbelow,BRabove
       
       include 'SM_decay_data_cb.h'
       
       Mhmin=hdecay_tab(1,1)
       Mhmax=hdecay_tab(4993,1)
       if((Mh .ge. Mhmin) .and. (Mh .le. Mhmax)) then
         sep=hdecay_tab(2,1)-hdecay_tab(1,1)
         
         ilow=int((Mh-Mhmin)/sep)+1
         
	 if( ((Mh-Mhmin)/sep+1d0-dble(ilow)) .lt. 1d-7) then
           BRinterpol=hdecay_tab(ilow,9)
         else
           Mh_below=dble(ilow-1)*sep+Mhmin
           BRbelow=hdecay_tab(ilow,9)
           BRabove=hdecay_tab(ilow+1,9)
           BRinterpol=BRbelow +(BRabove-BRbelow)/sep*(MH-Mh_below)
         endif
       else
         write(*,*)'function SMBR_Hgamgam(Mh) out of range'
         stop
       endif
       
       SMBR_Hgamgam=BRinterpol       
       end 

! not looked for in Higgs search analyses

!******************************************************
       double precision function SMBR_Hgg(Mh)
!******************************************************
       implicit none
       integer ilow
       double precision Mh,Mhmin,Mhmax,sep,BRinterpol
       double precision  mh_below,BRbelow,BRabove
       
       include 'SM_decay_data_cb.h'
       
       Mhmin=hdecay_tab(1,1)
       Mhmax=hdecay_tab(4993,1)
       if((Mh .ge. Mhmin) .and. (Mh .le. Mhmax)) then
         sep=hdecay_tab(2,1)-hdecay_tab(1,1)
         
         ilow=int((Mh-Mhmin)/sep)+1
         
         if( ((Mh-Mhmin)/sep+1d0-dble(ilow)) .lt. 1d-7) then
	   BRinterpol=hdecay_tab(ilow,8)
         else
           Mh_below=dble(ilow-1)*sep+Mhmin
           BRbelow=hdecay_tab(ilow,8)
           BRabove=hdecay_tab(ilow+1,8)
           BRinterpol=BRbelow +(BRabove-BRbelow)/sep*(MH-Mh_below)
         endif
       else
         write(*,*)'function SMBR_Hgg(Mh) out of range'
         stop
       endif
       
       SMBR_Hgg=BRinterpol       
       end 


!******************************************************
       double precision function SMBR_HZgam(Mh)
!******************************************************
       implicit none
       integer ilow
       double precision Mh,Mhmin,Mhmax,sep,BRinterpol
       double precision  mh_below,BRbelow,BRabove
       
       include 'SM_decay_data_cb.h'
       
       Mhmin=hdecay_tab(1,1)
       Mhmax=hdecay_tab(4993,1)
       if((Mh .ge. Mhmin) .and. (Mh .le. Mhmax)) then
         sep=hdecay_tab(2,1)-hdecay_tab(1,1)
         
         ilow=int((Mh-Mhmin)/sep)+1
         
         if( ((Mh-Mhmin)/sep+1d0-dble(ilow)) .lt. 1d-7) then
	   BRinterpol=hdecay_tab(ilow,10)
         else
           Mh_below=dble(ilow-1)*sep+Mhmin
           BRbelow=hdecay_tab(ilow,10)
           BRabove=hdecay_tab(ilow+1,10)
           BRinterpol=BRbelow +(BRabove-BRbelow)/sep*(MH-Mh_below)
         endif
       else
         write(*,*)'function SMBR_HZgam(Mh) out of range'
         stop
       endif
       
       SMBR_HZgam=BRinterpol       
       end 


!******************************************************
       double precision function SMBR_HZZ(Mh)
!******************************************************
       implicit none
       integer ilow
       double precision Mh,Mhmin,Mhmax,sep,BRinterpol
       double precision  mh_below,BRbelow,BRabove
       
       include 'SM_decay_data_cb.h'
       
       Mhmin=hdecay_tab(1,1)
       Mhmax=hdecay_tab(4993,1)
       if((Mh .ge. Mhmin) .and. (Mh .le. Mhmax)) then
         sep=hdecay_tab(2,1)-hdecay_tab(1,1)
         
         ilow=int((Mh-Mhmin)/sep)+1
         
         if( ((Mh-Mhmin)/sep+1d0-dble(ilow)) .lt. 1d-7) then
	   BRinterpol=hdecay_tab(ilow,12)
         else
           Mh_below=dble(ilow-1)*sep+Mhmin
           BRbelow=hdecay_tab(ilow,12)
           BRabove=hdecay_tab(ilow+1,12)
           BRinterpol=BRbelow +(BRabove-BRbelow)/sep*(MH-Mh_below)
         endif
       else
         write(*,*)'function SMBR_HZZ(Mh) out of range'
         stop
       endif
       
       SMBR_HZZ=BRinterpol       
       end 


!******************************************************
       double precision function SMBR_Htoptop(Mh)
!******************************************************
       implicit none
       integer ilow
       double precision Mh,Mhmin,Mhmax,sep,BRinterpol
       double precision  mh_below,BRbelow,BRabove
       
       include 'SM_decay_data_cb.h'
       
       Mhmin=hdecay_tab(1,1)
       Mhmax=hdecay_tab(4993,1)
       if((Mh .ge. Mhmin) .and. (Mh .le. Mhmax)) then
         sep=hdecay_tab(2,1)-hdecay_tab(1,1)
         
         ilow=int((Mh-Mhmin)/sep)+1
         
         if( ((Mh-Mhmin)/sep+1d0-dble(ilow)) .lt. 1d-7) then
	   BRinterpol=hdecay_tab(ilow,7)
         else
           Mh_below=dble(ilow-1)*sep+Mhmin
           BRbelow=hdecay_tab(ilow,7)
           BRabove=hdecay_tab(ilow+1,7)
           BRinterpol=BRbelow +(BRabove-BRbelow)/sep*(MH-Mh_below)
         endif
       else
         write(*,*)'function SMBR_Htoptop(Mh) out of range'
         stop
       endif
       
       SMBR_Htoptop=BRinterpol       
       end 


!******************************************************
       double precision function SMBR_Hcc(Mh)
!******************************************************
       implicit none
       integer ilow
       double precision Mh,Mhmin,Mhmax,sep,BRinterpol
       double precision  mh_below,BRbelow,BRabove
       
       include 'SM_decay_data_cb.h'
       
       Mhmin=hdecay_tab(1,1)
       Mhmax=hdecay_tab(4993,1)
       if((Mh .ge. Mhmin) .and. (Mh .le. Mhmax)) then
         sep=hdecay_tab(2,1)-hdecay_tab(1,1)
         
         ilow=int((Mh-Mhmin)/sep)+1
         
         if( ((Mh-Mhmin)/sep+1d0-dble(ilow)) .lt. 1d-7) then
           BRinterpol=hdecay_tab(ilow,6)
         else
           Mh_below=dble(ilow-1)*sep+Mhmin
           BRbelow=hdecay_tab(ilow,6)
           BRabove=hdecay_tab(ilow+1,6)
           BRinterpol=BRbelow +(BRabove-BRbelow)/sep*(MH-Mh_below)
         endif
       else
         write(*,*)'function SMBR_Hcc(Mh) out of range'
         stop
       endif
       
       SMBR_Hcc=BRinterpol       
       end 




!******************************************************
       double precision function SMBR_Hss(Mh)
!******************************************************
       implicit none
       integer ilow
       double precision Mh,Mhmin,Mhmax,sep,BRinterpol
       double precision  mh_below,BRbelow,BRabove
       
       include 'SM_decay_data_cb.h'
       
       Mhmin=hdecay_tab(1,1)
       Mhmax=hdecay_tab(4993,1)
       if((Mh .ge. Mhmin) .and. (Mh .le. Mhmax)) then       
         sep=hdecay_tab(2,1)-hdecay_tab(1,1)
         
         ilow=int((Mh-Mhmin)/sep)+1
         
         if( ((Mh-Mhmin)/sep+1d0-dble(ilow)) .lt. 1d-7) then
	   BRinterpol=hdecay_tab(ilow,5)
         else
           Mh_below=dble(ilow-1)*sep+Mhmin
           BRbelow=hdecay_tab(ilow,5)
           BRabove=hdecay_tab(ilow+1,5)
           BRinterpol=BRbelow +(BRabove-BRbelow)/sep*(MH-Mh_below)
         endif
       else
         write(*,*)'function SMBR_Hss(Mh) out of range'
         stop
       endif
       
       SMBR_Hss=BRinterpol       
       end 



!******************************************************
       double precision function SMBR_Hmumu(Mh)
!******************************************************
       implicit none
       integer ilow
       double precision Mh,Mhmin,Mhmax,sep,BRinterpol
       double precision  mh_below,BRbelow,BRabove
       
       include 'SM_decay_data_cb.h'
       
       Mhmin=hdecay_tab(1,1)
       Mhmax=hdecay_tab(4993,1)
       if((Mh .ge. Mhmin) .and. (Mh .le. Mhmax)) then
         sep=hdecay_tab(2,1)-hdecay_tab(1,1)
         
         ilow=int((Mh-Mhmin)/sep)+1
         
         if( ((Mh-Mhmin)/sep+1d0-dble(ilow)) .lt. 1d-7) then
	   BRinterpol=hdecay_tab(ilow,4)
         else
           Mh_below=dble(ilow-1)*sep+Mhmin
           BRbelow=hdecay_tab(ilow,4)
           BRabove=hdecay_tab(ilow+1,4)
           BRinterpol=BRbelow +(BRabove-BRbelow)/sep*(MH-Mh_below)
         endif
       else
         write(*,*)'function SMBR_Hmumu(Mh) out of range'
         stop
       endif       
       
       SMBR_Hmumu=BRinterpol       
       end 


