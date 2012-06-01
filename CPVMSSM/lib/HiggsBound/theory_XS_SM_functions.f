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

! /!\ :  All cross sections have been in units of femtobarn and are now given
!        in units of pikobarn. Hence, there is always an explicit division by 1000d0.

! functions are fitted to data 
! downloaded from http://maltoni.home.cern.ch/maltoni/TeV4LHC/bbh-tev_ed.dat
! on Thursday 10th April 2008
! data has range 100 to 300 GeV
 
!****************************************************** 
       function SMCS_tev_pp_qq_HW(x)
!******************************************************
       implicit none
       double precision  SMCS_tev_pp_qq_HW,x 
       double precision  a0,a0p5,a1,a2
       
       a0 = 5.7514105496046D0
       a0p5 = -0.375021739545092D0
       a1 = 0.0049451487167627D0
       a2 = -3.77008582179264D-06
       
       if(x.lt.65.0D0)then
        write(*,*) 'function SMCS_tev_pp_qq_HW might not be a good fit (m_H < 65 GeV)' 
       elseif(x.le.335.0D0)then
        SMCS_tev_pp_qq_HW=10.0D0**(a0+a0p5*x**0.5D0+a1*x+a2*x**2.0D0)
       else
        write(*,*) 'function SMCS_tev_pp_qq_HW might not be a good fit (m_H > 335 GeV)' 
       endif 

       SMCS_tev_pp_qq_HW=SMCS_tev_pp_qq_HW/1000d0
       
       end


!****************************************************** 
       function SMCS_tev_pp_qq_HZ(x)
!****************************************************** 
       implicit none
       double precision  SMCS_tev_pp_qq_HZ,x 
       double precision  a0,a0p5,a1,a2
 
       a0 = 5.29935340004443D0
       a0p5 = -0.351677660532052D0
       a1 = 0.0047848452802514D0
       a2 = -3.82425969474559D-06 
       
       if(x.lt.65.0D0)then
        write(*,*) 'function SMCS_tev_pp_qq_HZ might not be a good fit (m_H < 65 GeV)' 
       elseif(x.le.335.0D0)then
        SMCS_tev_pp_qq_HZ=10.0D0**(a0+a0p5*x**0.5D0+a1*x+a2*x**2.0D0)
       else
        write(*,*) 'function SMCS_tev_pp_qq_HZ might not be a good fit (m_H > 335 GeV)' 
       endif 
 
       SMCS_tev_pp_qq_HZ=SMCS_tev_pp_qq_HZ/1000d0
        
       end 


!****************************************************** 
       function SMCS_tev_pp_gg_H(x)
!******************************************************
       implicit none
       double precision  SMCS_tev_pp_gg_H,x 
       double precision  a0,a0p5,a1,a2
 
       a0 = 5.59682783597183D0
       a0p5 = -0.244216706673437D0
       a1 = 0.000365613425058581D0
       a2 = 2.66122261164927D-06  
	
       if(x.lt.65.0D0)then
        write(*,*) 'function SMCS_tev_pp_gg_H might not be a good fit (m_H < 65 GeV)' 
       elseif(x.le.335.0D0)then
        SMCS_tev_pp_gg_H=10.0D0**(a0+a0p5*x**0.5D0+a1*x+a2*x**2.0D0)
       else
        write(*,*) 'function SMCS_tev_pp_gg_H might not be a good fit (m_H > 335 GeV)' 
       endif 
 
       SMCS_tev_pp_gg_H=SMCS_tev_pp_gg_H/1000d0 
       
       end 


!****************************************************** 
       function SMCS_tev_pp_bb_H(x)
!******************************************************
       implicit none
       double precision  SMCS_tev_pp_bb_H,x 
       double precision  a0,a0p5,a1,a2
 
       a0 = 5.41583328209568D0
       a0p5 = -0.453323023285831D0
       a1 = 0.00514220061294974D0
       a2 = -3.31488355831377D-06
	
       if(x.lt.65.0D0)then
        write(*,*) 'function SMCS_tev_pp_bb_H might not be a good fit (m_H < 65 GeV)' 
       elseif(x.le.335.0D0)then
        SMCS_tev_pp_bb_H=10.0D0**(a0+a0p5*x**0.5D0+a1*x+a2*x**2.0D0)
       else
        write(*,*) 'function SMCS_tev_pp_bb_H might not be a good fit (m_H > 335 GeV)' 
       endif  

       SMCS_tev_pp_bb_H=SMCS_tev_pp_bb_H/1000d0
 
       end 



!****************************************************** 
       function SMCS_tev_pp_vbf_H(x)
!******************************************************
       implicit none
       double precision  SMCS_tev_pp_vbf_H,x 
       double precision  a0,a0p5,a1,a2
 
       a0 = 3.03330688577339D0
       a0p5 = -0.061167148821396D0
       a1 = -0.00424395931917914D0
       a2 = 7.67964289500027D-08
	
       if(x.lt.65.0D0)then
        write(*,*)'function SMCS_tev_pp_vbf_H might not be a good fit (m_H < 65 GeV)' 
       elseif(x.le.335.0D0)then
        SMCS_tev_pp_vbf_H=10.0D0**(a0+a0p5*x**0.5D0+a1*x+a2*x**2.0D0)
       else
        write(*,*)'function SMCS_tev_pp_vbf_H might not be a good fit (m_H > 335 GeV)' 
       endif  

       SMCS_tev_pp_vbf_H=SMCS_tev_pp_vbf_H/1000d0
 
       end 

