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

!>> input.F now included in main *.F-file
!! Definitions of NHMAX, TEV_NSCHN, TEV_NDCHN, LEP_NSCHN, LEP_NDCHN :
!#include "tables-nhmax.F"
!<<


*************************************************************************
	subroutine initialize_filestatus()
*************************************************************************
	implicit none
	logical open1,open2,open3,open4,open5,open6,open7,open8,open9
	logical open10,open11,open12,open13,file_exists

        common /fstat_MH_GammaTot/ open1
	common /fstat_BR_1H/ open2
	common /fstat_BR_2H/ open3
	common /fstat_LEP_HZ_CS_ratios/ open4
	common /fstat_LEP_2H_CS_ratios/ open5
        common /fstat_TEV_H_0jet_partCS_ratios/ open6
        common /fstat_TEV_H_1jet_partCS_ratios/ open7
        common /fstat_TEV_H_vbf_hadCS_ratios/ open8
        common /fstat_TEV_HW_partCS_ratios/ open9
        common /fstat_TEV_HZ_partCS_ratios/ open10
	common /fstat_TEV_1H_hadCS_ratios/ open11
        common /fstat_effC/ open12
        common /fstat_additional/ open13,file_exists

	  open1 = .FALSE.
	  open2 = .FALSE.
	  open3 = .FALSE.
	  open4 = .FALSE.
          open5 = .FALSE.
          open6 = .FALSE.
          open7 = .FALSE.
          open8 = .FALSE.
          open9 = .FALSE.
          open10 = .FALSE.
          open11 = .FALSE.
          open12 = .FALSE.
          open13 = .FALSE.
	  file_exists=.TRUE.

	end

**** masses and widths

*************************************************************************
	subroutine read_MH_GammaTot(NHIGGS,fileprefix,prefix,eof,
     &		id,Mh,Gamma_h)
*************************************************************************
	implicit none
	logical use_default,eof,open,prefix
	common /fstat_MH_GammaTot/ open
	character*150 filename,filename2
        character*100 fileprefix 
	integer id,NHIGGS,i,j
        integer lnblnk
!        external lnblnk

	double precision Mh(NHMAX),Gamma_h(NHMAX)
	
	if(prefix) then
		filename=fileprefix(:lnblnk(fileprefix))//'MH_GammaTot.dat'
	else
		filename='MH_GammaTot.dat'
	endif
	if(.not. open) then
	  open(10,file=filename,status='old')
	  open = .TRUE.
	endif

	read(10,*,err=91,end=92) 
     & 		id,(Mh(i),i=1,NHIGGS),(Gamma_h(i), i=1,NHIGGS)
	return



 91    write(*,*)'Error reading file: ',filename
	stop
 92	eof=.TRUE.
	close(10)


	end


**** branching ratios 

*************************************************************************
	subroutine read_BR_1H(NHIGGS,fileprefix,prefix,eof,
     &		id,BRhjbbbar,BRhjtautau,BRhjWW,BRhjgamgam)
*************************************************************************
	implicit none
	logical use_default,eof,open,prefix
	common /fstat_BR_1H/ open
	character*150 filename
        character*100 fileprefix 
	integer id,NHIGGS,i
        integer lnblnk
!        external lnblnk


	double precision BRhjWW(NHMAX),BRhjgamgam(NHMAX)
        double precision BRhjbbbar(NHMAX),BRhjtautau(NHMAX)

	if(prefix) then
		filename=fileprefix(:lnblnk(fileprefix))//'BR_1H.dat'
	else
		filename='BR_1H.dat'
	endif
	if(.not. open) then 
	  open(11,file=filename,status='old')
          open = .TRUE.
        endif

	read(11,*,err=101,end=102) 
     & 		id,(BRhjbbbar(i),i=1,NHIGGS),
     &		(BRhjtautau(i),i=1,NHIGGS),
     &		(BRhjWW(i),i=1,NHIGGS),
     &		(BRhjgamgam(i),i=1,NHIGGS)

	return
	
 101    write(*,*)'Error reading file: ',filename
	stop
 102	eof=.TRUE.
	close(11)


	end


*************************************************************************
	subroutine read_BR_2H(NHIGGS,fileprefix,prefix,eof,id,BRhjhihi)
*************************************************************************
	implicit none
	logical use_default,eof,open,prefix
	common /fstat_BR_2H/ open
	character*150 filename
        character*100 fileprefix 
	integer id,NHIGGS
        integer lnblnk
!        external lnblnk

	integer i,j,k,l
        double precision BRhjhihi(NHMAX,NHMAX)
	double precision BRhjhihi_compact(NHMAX**2-NHMAX)



	if(prefix) then
		filename=fileprefix(:lnblnk(fileprefix))//'BR_2H.dat'
	else
		filename='BR_2H.dat'
	endif

	if(NHIGGS .ne. 1) then

	if(.not. open) then
	  open(12,file=filename,status='old')
          open = .TRUE.
        endif

	read(12,*,err=201,end=202) id,(BRhjhihi_compact(i),i=1,NHIGGS*(NHIGGS-1))

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
c	  write(*,*) ( BRhjhihi(j,l), l=1,NHIGGS )
	enddo
	
	else

	write(*,*)'NHIGGS=1 : file :',filename,' not needed in this case.'
	BRhjhihi(1,1)=0d0
	endif

	return
	
 201    write(*,*)'Error reading file: ',filename
	stop
 202	eof=.TRUE.
        close(12)


	end




*** LEP input

c LEP c.s. ratios
*************************************************************************
	subroutine read_LEP_HZ_CS_ratios(NHIGGS,fileprefix,prefix,eof,
     &		id,R_hjZ_ee)
*************************************************************************
	implicit none
	logical use_default,eof,open,prefix
	common /fstat_LEP_HZ_CS_ratios/ open
	character*150 filename
        character*100 fileprefix 
	integer id,NHIGGS,i
        integer lnblnk
!        external lnblnk

	double precision R_hjZ_ee(NHMAX)

	if(prefix) then
		filename=fileprefix(:lnblnk(fileprefix))//'LEP_HZ_CS_ratios.dat'
	else
		filename='LEP_HZ_CS_ratios.dat'
	endif
	if(.not. open) then
	  open(13,file=filename,status='old')
          open = .TRUE.
        endif

	read(13,*,err=301,end=302) id,(R_hjZ_ee(i),i=1,NHIGGS)
	return
	
 301    write(*,*)'Error reading file: ',filename
	stop
 302	eof=.TRUE.
	close(13)


	end


*************************************************************************
	subroutine read_LEP_2H_CS_ratios(NHIGGS,fileprefix,prefix,eof,
     &		id,R_hjhi_ee)
*************************************************************************
	implicit none
	logical use_default,eof,open,prefix
	common /fstat_LEP_2H_CS_ratios/ open
	character*150 filename
        character*100 fileprefix 
	integer id,NHIGGS
        integer lnblnk
!        external lnblnk

	integer i,j,k,l
        double precision R_hjhi_ee(NHMAX,NHMAX)
        double precision R_hjhi_ee_compact((NHMAX*(NHMAX+1))/2)

	if(prefix) then
		filename=fileprefix(:lnblnk(fileprefix))//'LEP_2H_CS_ratios.dat'
	else
		filename='LEP_2H_CS_ratios.dat'
	endif
	if(.not. open) then
	  open(14,file=filename,status='old')
          open = .TRUE.
        endif

	read(14,*,err=401,end=402) id,
     &		(R_hjhi_ee_compact(i),i=1,(NHIGGS*(NHIGGS+1))/2)


	k=0
	do j=1,NHIGGS
	  do i=1,NHIGGS
	    if(i .le. j) then 
	      k=k+1
	      R_hjhi_ee(j,i)=R_hjhi_ee_compact(k)
	      if(i .ne. j) R_hjhi_ee(i,j)=R_hjhi_ee(j,i)
	    endif
	  enddo
c	  write(*,*) ( R_hjhi_ee(j,l), l=1,NHIGGS )
	enddo
	return
	
 401    write(*,*)'Error reading file: ',filename
	stop
 402	eof=.TRUE.
        close(14)
 	

	end




*** TEV input

c TEV partonic c.s. ratios
*************************************************************************
	subroutine read_TEV_H_0jet_partCS_ratios(NHIGGS,fileprefix,prefix,eof,
     &		id,R_hj_gg,R_hj_bb)
*************************************************************************
	implicit none
	logical use_default,eof,open,prefix
	common /fstat_H_0jet_partCS_ratios/ open
	character*150 filename
        character*100 fileprefix 
	integer id,NHIGGS,i
        integer lnblnk
!        external lnblnk

        double precision R_hj_gg(NHMAX),R_hj_bb(NHMAX)

	if(prefix) then
		filename=fileprefix(:lnblnk(fileprefix))//'TEV_H_0jet_partCS_ratios.dat'
	else
		filename='TEV_H_0jet_partCS_ratios.dat'
	endif
	if(.not. open) then
	  open(15,file=filename,status='old')
          open = .TRUE.
        endif

	read(15,*,err=501,end=502) id,(R_hj_gg(i),i=1,NHIGGS),
     &		(R_hj_bb(i),i=1,NHIGGS)
	return
	
 501    write(*,*)'Error reading file: ',filename
	stop
 502	eof=.TRUE.
	close(15)

	end

*************************************************************************
	subroutine read_TEV_H_1jet_partCS_ratios(NHIGGS,fileprefix,prefix,eof,
     &		id,R_hjb_bg)
*************************************************************************
	implicit none
	logical use_default,eof,open,prefix
	common /fstat_TEV_H_1jet_partCS_ratios/ open
	character*150 filename
        character*100 fileprefix 
	integer id,NHIGGS,i
        integer lnblnk
!        external lnblnk

        double precision R_hjb_bg(NHMAX)

	if(prefix) then
		filename=fileprefix(:lnblnk(fileprefix))//'TEV_H_1jet_partCS_ratios.dat'
	else
		filename='TEV_H_1jet_partCS_ratios.dat'
	endif
	if(.not. open) then
	  open(16,file=filename,status='old')
          open = .TRUE.
        endif

	read(16,*,err=601,end=602) id,(R_hjb_bg(i),i=1,NHIGGS)
	return
	
 601    write(*,*)'Error reading file: ',filename
	stop
 602	eof=.TRUE.
	close(16)

	end

*************************************************************************
	subroutine read_TEV_H_vbf_hadCS_ratios(NHIGGS,fileprefix,prefix,eof,
     &		id,R_sigma_H_VBF)
*************************************************************************
	implicit none
	logical use_default,eof,open,prefix
	common /fstat_TEV_H_vbf_hadCS_ratios/ open
	character*150 filename
        character*100 fileprefix 
	integer id,NHIGGS,i
        integer lnblnk
!        external lnblnk

        double precision R_sigma_H_VBF(NHMAX)

	if(prefix) then
		filename=fileprefix(:lnblnk(fileprefix))//'TEV_H_vbf_hadCS_ratios.dat'
	else
		filename='TEV_H_vbf_hadCS_ratios.dat'
	endif
	if(.not. open) then
	  open(17,file=filename,status='old')
          open = .TRUE.
        endif

	read(17,*,err=701,end=702) id,(R_sigma_H_VBF(i),i=1,NHIGGS)
	return
	
 701    write(*,*)'Error reading file: ',filename
	stop
 702	eof=.TRUE.
	close(17)
 	
	end

*************************************************************************
	subroutine read_TEV_HW_partCS_ratios(NHIGGS,fileprefix,prefix,eof,
     &		id,R_hjWp_udb,R_hjWp_csb,R_hjWm_dub,R_hjWm_scb)
*************************************************************************
	implicit none
	logical use_default,eof,open,prefix
	common /fstat_TEV_HW_partCS_ratios/ open
	character*150 filename
        character*100 fileprefix 
	integer id,NHIGGS,i
        integer lnblnk
!        external lnblnk

        double precision R_hjWp_udb(NHMAX),R_hjWp_csb(NHMAX),
     &                   R_hjWm_dub(NHMAX),R_hjWm_scb(NHMAX)

	if(prefix) then
		filename=fileprefix(:lnblnk(fileprefix))//'TEV_HW_partCS_ratios.dat'
	else
		filename='TEV_HW_partCS_ratios.dat'
	endif
	if(.not. open) then
	  open(18,file=filename,status='old')
          open = .TRUE.
        endif

	read(18,*,err=801,end=802) id,
     &          (R_hjWp_udb(i),i=1,NHIGGS),
     &		(R_hjWp_csb(i),i=1,NHIGGS),
     &		(R_hjWm_dub(i),i=1,NHIGGS),
     &		(R_hjWm_scb(i),i=1,NHIGGS)
	return
	
 801    write(*,*)'Error reading file: ',filename
	stop
 802	eof=.TRUE.
	close(18)

	end

*************************************************************************
	subroutine read_TEV_HZ_partCS_ratios(NHIGGS,fileprefix,prefix,eof,
     &		id,R_hjZ_ddb,R_hjZ_uub,R_hjZ_ssb,R_hjZ_ccb,R_hjZ_bbb)
*************************************************************************
	implicit none
	logical use_default,eof,open,prefix
	common /fstat_TEV_HZ_partCS_ratios/ open
	character*150 filename
        character*100 fileprefix 
	integer id,NHIGGS,i
        integer lnblnk
!        external lnblnk

        double precision R_hjZ_ddb(NHMAX),R_hjZ_uub(NHMAX),
     &                   R_hjZ_ssb(NHMAX),R_hjZ_ccb(NHMAX),R_hjZ_bbb(NHMAX)

	if(prefix) then
		filename=fileprefix(:lnblnk(fileprefix))//'TEV_HZ_partCS_ratios.dat'
	else
		filename='TEV_HZ_partCS_ratios.dat'
	endif
	if(.not. open) then
	  open(19,file=filename,status='old')
          open = .TRUE.
        endif

	read(19,*,err=901,end=902) id,
     &          (R_hjZ_ddb(i),i=1,NHIGGS),
     &		(R_hjZ_uub(i),i=1,NHIGGS),
     &		(R_hjZ_ssb(i),i=1,NHIGGS),
     &		(R_hjZ_ccb(i),i=1,NHIGGS),
     &		(R_hjZ_bbb(i),i=1,NHIGGS)
	return
	
 901    write(*,*)'Error reading file: ',filename
	stop
 902	eof=.TRUE.
	close(19)

	end


c hadronic c.s. ratios
*************************************************************************
	subroutine read_TEV_1H_hadCS_ratios(NHIGGS,fileprefix,prefix,eof,
     &		id,R_sigma_H,R_sigma_Hb,R_sigma_HW,R_sigma_HZ,R_sigma_H_VBF)
*************************************************************************
	implicit none
	logical use_default,eof,open,prefix
	common /fstat_TEV_1H_hadCS_ratios/ open
	character*150 filename
        character*100 fileprefix 
	integer id,NHIGGS,i
        integer lnblnk
!        external lnblnk

	double precision R_sigma_H(NHMAX),R_sigma_Hb(NHMAX)
	double precision R_sigma_HW(NHMAX),R_sigma_HZ(NHMAX)
	double precision R_sigma_H_VBF(NHMAX)

	if(prefix) then
		filename=fileprefix(:lnblnk(fileprefix))//'TEV_1H_hadCS_ratios.dat'
	else
		filename='TEV_1H_hadCS_ratios.dat'
	endif
	if(.not. open) then
	  open(20,file=filename,status='old')
          open = .TRUE.
        endif

	read(20,*,err=2001,end=2002) id,(R_sigma_H(i),i=1,NHIGGS),
     &		(R_sigma_Hb(i),i=1,NHIGGS),
     &		(R_sigma_HW(i),i=1,NHIGGS),
     &		(R_sigma_HZ(i),i=1,NHIGGS),
     &		(R_sigma_H_VBF(i),i=1,NHIGGS)
	return
	
 2001    write(*,*)'Error reading file: ',filename
	stop
 2002	eof=.TRUE.
	close(20)

	end


c couplings squared ratios
*************************************************************************
	subroutine read_couplings_sq_ratios(NHIGGS,fileprefix,prefix,eof,
     &  	id, g_hjbb_sq, g_hjtautau_sq, g_hjWW_sq, g_hjZZ_sq,
     &		g_hjgamgam_sq, g_hjgg_sq, g_hjhiZ_sq)
*************************************************************************
	implicit none
	logical use_default,eof,open,prefix
	common /fstat_effC/ open
	character*150 filename
        character*100 fileprefix 
	integer id,i,j,k,NHIGGS
        integer lnblnk
!        external lnblnk

	double precision g_hjbb_sq(NHMAX), g_hjtautau_sq(NHMAX), 
     &   	g_hjWW_sq(NHMAX), g_hjZZ_sq(NHMAX),
     &          g_hjgamgam_sq(NHMAX), g_hjgg_sq(NHMAX), 
     &		g_hjhiZ_sq(NHMAX,NHMAX)
	double precision g_hjhiZ_sq_compact((NHMAX*(NHMAX+1))/2)



	if(prefix) then
		filename=fileprefix(:lnblnk(fileprefix))//'effC.dat'
	else
		filename='effC.dat'
	endif
	if(.not. open) then
	  open(21,file=filename,status='old')
          open = .TRUE.
        endif

	read(21,*,err=2101,end=2102) id,
     &		(g_hjbb_sq(i),i=1,NHIGGS), 
     &          (g_hjtautau_sq(i),i=1,NHIGGS),
     &          (g_hjWW_sq(i),i=1,NHIGGS), 
     &          (g_hjZZ_sq(i),i=1,NHIGGS), 
     &          (g_hjgamgam_sq(i),i=1,NHIGGS), 
     &          (g_hjgg_sq(i),i=1,NHIGGS),
     &          (g_hjhiZ_sq_compact(i),i=1,(NHIGGS*(NHIGGS+1))/2)



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
	return
	
 2101    write(*,*)'Error reading file: ',filename
	stop
 2102	eof=.TRUE.
	close(21)

	end

*** additional data (optional)

*************************************************************************
	subroutine read_additional_data(NHIGGS,fileprefix,prefix,eof,
     &		id,add_dat,columns)
*************************************************************************
	implicit none
	logical use_default,eof,open,file_exists,prefix
	common /fstat_additional/ open,file_exists
	character*150 filename
        character*100 fileprefix 
	integer id,NHIGGS,i,j
        integer lnblnk
!        external lnblnk

	integer length,columns,thecolumns
	common /column_info/ thecolumns
        parameter(length=MAX_ADDITIONAL_COLUMNS*10)
        character(LEN=length) line
	double precision add_dat(MAX_ADDITIONAL_COLUMNS)

	if(.not. file_exists) then
	 	thecolumns=0
		columns=thecolumns
	else

	if(prefix) then
		filename=fileprefix(:lnblnk(fileprefix))//'additional.dat'
	else
		filename='additional.dat'
	endif

	if(.not. open) then
		open(22,file=filename,status='old',err=93)
		read(22,'(A)',err=91) line
        	call determine_columns_from_line(line,columns)
        	close(22)
		columns=columns-1
		thecolumns=columns
        	write(*,*)'# of additional data columns:',columns

		open(22,file=filename,status='old')
		open = .TRUE.
	else
		columns=thecolumns
	endif
	
	read(22,*,err=91,end=92) 
     & 		id,(add_dat(i),i=1,columns)

c check	write(*,*)'read add:id=',id
c check	write(*,*)'read add:add_dat(i):',(add_dat(i),i=1,columns)
	endif

	return


 91    write(*,*)'Error reading file: ',filename
	stop
 92	eof=.TRUE.
	close(22)
        return
 93     file_exists = .FALSE.
	thecolumns=0
	columns=thecolumns
        return



	end


*************************************************************
        subroutine determine_columns_from_line(line,columns)
*************************************************************
        implicit none
        logical wasblank
        integer i,columns,length
        parameter(length=MAX_ADDITIONAL_COLUMNS*10)
        character(LEN=length) line

        wasblank=.FALSE.
        columns=0
        do i=1,length
c distinguish between blank or tab character and other characters
          if((line(i:i) .ne. ' ') .and. (ICHAR(line(i:i)) .ne. 9)) then
            if(.not. wasblank) columns=columns+1
            wasblank=.TRUE.
          else
            wasblank=.FALSE.
          endif
          enddo

        end




*************************************************
        subroutine output_vector(NHIGGS,m)
*************************************************
        implicit none
        integer i,NHIGGS
        double precision m(NHMAX)

#ifdef G77VERSION
	if(NHIGGS .eq. 1) write(*,901) ( m(i), i=1,NHIGGS,1 )
        if(NHIGGS .eq. 2) write(*,902) ( m(i), i=1,NHIGGS,1 )
        if(NHIGGS .eq. 3) write(*,903) ( m(i), i=1,NHIGGS,1 )
        if(NHIGGS .eq. 4) write(*,904) ( m(i), i=1,NHIGGS,1 )
        if(NHIGGS .eq. 5) write(*,905) ( m(i), i=1,NHIGGS,1 )
        if(NHIGGS .eq. 6) write(*,906) ( m(i), i=1,NHIGGS,1 )
        if(NHIGGS .eq. 7) write(*,907) ( m(i), i=1,NHIGGS,1 )
        if(NHIGGS .eq. 8) write(*,908) ( m(i), i=1,NHIGGS,1 )
        if(NHIGGS .eq. 9) write(*,909) ( m(i), i=1,NHIGGS,1 )
        if(NHIGGS .eq. 10) write(*,910) ( m(i), i=1,NHIGGS,1 )
 901	FORMAT(G12.5)
 902    FORMAT(2G12.5)
 903    FORMAT(3G12.5)
 904    FORMAT(4G12.5)
 905    FORMAT(5G12.5)
 906    FORMAT(6G12.5)
 907    FORMAT(7G12.5)
 908    FORMAT(8G12.5)
 909    FORMAT(9G12.5)
 910    FORMAT(10G12.5)
#else
	write(*,333) ( m(i), i=1,NHIGGS,1 )
 333	format (<NHIGGS>G12.5)
#endif

        end
	
*************************************************
        subroutine output_matrix(NHIGGS,m)
*************************************************
        implicit none
        integer i,j,NHIGGS
        double precision m(NHMAX,NHMAX)

        do i=1,NHIGGS,1
#ifdef G77VERSION
	if(NHIGGS .eq. 1) write(*,901) ( m(i,j), j=1,NHIGGS,1 )
        if(NHIGGS .eq. 2) write(*,902) ( m(i,j), j=1,NHIGGS,1 ) 
        if(NHIGGS .eq. 3) write(*,903) ( m(i,j), j=1,NHIGGS,1 )
        if(NHIGGS .eq. 4) write(*,904) ( m(i,j), j=1,NHIGGS,1 )
        if(NHIGGS .eq. 5) write(*,905) ( m(i,j), j=1,NHIGGS,1 )
        if(NHIGGS .eq. 6) write(*,906) ( m(i,j), j=1,NHIGGS,1 )
        if(NHIGGS .eq. 7) write(*,907) ( m(i,j), j=1,NHIGGS,1 )
        if(NHIGGS .eq. 8) write(*,908) ( m(i,j), j=1,NHIGGS,1 )
        if(NHIGGS .eq. 9) write(*,909) ( m(i,j), j=1,NHIGGS,1 )
        if(NHIGGS .eq. 10) write(*,910) ( m(i,j), j=1,NHIGGS,1 )
 901	FORMAT(G12.5)
 902    FORMAT(2G12.5)
 903    FORMAT(3G12.5)
 904    FORMAT(4G12.5)
 905    FORMAT(5G12.5)
 906    FORMAT(6G12.5)
 907    FORMAT(7G12.5)
 908    FORMAT(8G12.5)
 909    FORMAT(9G12.5)
 910    FORMAT(10G12.5)
#else
	write(*,334) ( m(i,j), j=1,NHIGGS,1 )
 334	format (<NHIGGS>G12.5)
#endif
	enddo	

        end

!del *************************************************
!del         subroutine put_zeros(dim,f)
!del *************************************************
!del         implicit none
!del         integer dim,i,j
!del         double precision f(dim,dim)
!del 
!del         do i=1,dim,1
!del           do j=1,dim,1
!del             f(i,j)=0d0
!del           enddo
!del         enddo
!del 
!del         end

