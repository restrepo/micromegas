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

!******************************************************************
	integer function log2int(logvar)
!******************************************************************
	implicit none
	logical logvar

	if(logvar) then
		log2int=1
	else
		log2int=0
	endif

	end
	
!******************************************************************
	double precision function log2dble(logvar)
!******************************************************************
	implicit none
	logical logvar

	if(logvar) then
		log2dble=1d0
	else
		log2dble=0d0
	endif

	end
	

*************************************************
        subroutine eval_1d_arr(minindex,maxindex,f,g)
*************************************************
        implicit none
        integer minindex,maxindex,i
        double precision f(minindex:maxindex),g(minindex:maxindex)

        do i=minindex,maxindex,1
	  f(i)=g(i)
        enddo

        end


*************************************************
        subroutine eval_2d_arr(minind1,maxind1,minind2,maxind2,f,g)
*************************************************
        implicit none
        integer minind1,maxind1,minind2,maxind2,i,j
        double precision f(minind1:maxind1,minind2:maxind2)
        double precision g(minind1:maxind1,minind2:maxind2)

        do i=minind1,maxind1,1
        do j=minind2,maxind2,1
        f(i,j)=g(i,j)
        enddo
        enddo

        end



************************************************************************
	integer function my_maxloc(array,length)
************************************************************************
c determine the lowest index-position of the maximum of an 1-dim array of
c double precision numbers with index-range 1:length.
	implicit none
	integer length,i,loc
	double precision array(length),max

c start-value for max is the value of the array at the 
c highest index position loc=length
	max=array(length)
	loc=length

c do (length-1) 2-by-2 number comparisons and each time change my_maxloc to the 
c maxium of the two numbers
	do i=1,length-1
	  if(array(length-i) .ge. max) then
	    max = array(length-i)
	    loc=length-i
	  endif
	enddo
	
	my_maxloc=loc

	end
	
	
************************************************************************
	double precision function my_maxval(array,length)
************************************************************************
c determine the maximum of an 1-dim array of double precision numbers
c with index-range 1:length.
	implicit none
	integer length,i
	double precision array(length),themax

	themax=array(1)

	do i=2,length
	    themax = dmax1(themax,array(i))
	enddo
	
	my_maxval=themax

	end
