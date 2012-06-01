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
	program SM_vs_fourth_generation_model
!******************************************************
	implicit none
	integer nH,HBresult,chan,ncombined,i
	real*8 obsratio,mass
	real*8 SMGammaTotal,SMGamma_hgg,Gamma_hgg,SMGamma_h,
     &	       SMBR_Hgg,FGMGammaTotal

	real*8 Mh,GammaTotal,g2hjbb,g2hjtautau,
     &         g2hjWW,g2hjZZ,g2hjgaga,g2hjgg,
     &         g2hjhiZ,BR_hjhihi


	nH=1
	call initialize_HiggsBounds(nH,'LandT')
c	call initialize_HiggsBounds(nH,'onlyT')
	open(10,file='example-SM-results.dat')
	open(20,file='example-4thGen-results.dat')

	do i=90,200,5
          Mh=dble(i)
	  SMGammaTotal=SMGamma_h(Mh)
	  SMGamma_hgg=SMBR_Hgg(Mh)*SMGammaTotal
	  Gamma_hgg=9d0*SMBR_Hgg(Mh)*SMGammaTotal
	  FGMGammaTotal=SMGammaTotal-SMGamma_hgg+Gamma_hgg

! SM 
	  GammaTotal=SMGammaTotal
	  g2hjbb=1d0
	  g2hjtautau=1d0
	  g2hjWW=1d0
	  g2hjZZ=1d0
	  g2hjgaga=1d0
	  g2hjgg=1d0
	  g2hjhiZ=0d0
     	  BR_hjhihi=0d0

	  call run_HiggsBounds_effC(nH,Mh,SMGammaTotal,
     &        g2hjbb,g2hjtautau,g2hjWW,g2hjZZ,g2hjgaga,g2hjgg,g2hjhiZ,
     &        BR_hjhihi,HBresult,chan,obsratio,ncombined)
	  write(10,*) Mh,HBresult,chan,obsratio

! 4th Generation Model
	  GammaTotal=FGMGammaTotal
	  g2hjgaga=1d0/9d0*FGMGammaTotal/SMGammaTotal
	  g2hjgg=9d0

	  call run_HiggsBounds_effC(nH,Mh,GammaTotal,
     &        g2hjbb,g2hjtautau,g2hjWW,g2hjZZ,g2hjgaga,g2hjgg,g2hjhiZ,
     &        BR_hjhihi,HBresult,chan,obsratio,ncombined)
	  write(20,*) Mh,HBresult,chan,obsratio
	enddo

	close(10)
        close(20)

	end
