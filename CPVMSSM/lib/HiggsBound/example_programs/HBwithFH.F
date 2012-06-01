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

* This is modified version of the demo.F file which is supplied with 
* FeynHiggs (vers 2.6.5, downloaded 05.01.09, www.feynhiggs.de)
* This file is part of the HiggsBounds distribution.

* The phrase 'currently' as used in this file means 'correct
* on the date 05.01.09'.
* See www.feynhiggs.de for the latest version of FeynHiggs
* and the corresponding documentation.

	program HBwithFHdemo
	implicit none

	integer error

c used by FHHiggsCorr:
	double precision MHiggs(4)
	double complex SAeff, UHiggs(3,3), ZHiggs(3,3)

c used by FHSelectUZ:
	integer uzint, uzext

c used by FHCouplings:
#include "FHCouplings.h"
	integer excl
	double complex couplings(ncouplings)
	double complex couplingsms(ncouplingsms)
	double precision gammas(ngammas)
	double precision gammasms(ngammasms)

c used by FHHiggsProd:
	double precision prodxs(nprodxs)       
       
c used by FHRetrieveSMPara:
        double precision invAlfa, AlfasMZ, GF
        double precision MS, MC, MB, MW, MZ
        double precision CKMlambda, CKMA, CKMrho, CKMeta

c used by initialize_HiggsBounds and run_HiggsBounds_part
c HB input:
        integer nHiggs        
        character*5 whichanalyses
        double precision Mh(3),                            
     &          CS_lep_hjZ_ratio(3),                              
     &          CS_lep_hjhi_ratio_nHbynH(3,3),                  
     &          CS_tev_gg_hj_ratio(3),CS_tev_bb_hj_ratio(3),   
     &          CS_tev_bg_hjb_ratio(3),                        
     &          CS_tev_ud_hjWp_ratio(3),CS_tev_cs_hjWp_ratio(3),
     &          CS_tev_ud_hjWm_ratio(3),CS_tev_cs_hjWm_ratio(3), 
     &          CS_tev_dd_hjZ_ratio(3),CS_tev_uu_hjZ_ratio(3),  
     &          CS_tev_ss_hjZ_ratio(3),CS_tev_cc_hjZ_ratio(3),  
     &          CS_tev_bb_hjZ_ratio(3),                       
     &          CS_tev_pp_vbf_ratio(3),                        
     &          BR_hjbb(3),BR_hjtautau(3),                      
     &          BR_hjWW(3),BR_hjgaga(3),                        
     &          BR_hjhihi_nHbynH(3,3)    

c HB output:
        integer HBresult,chan,ncombined  
        double precision obsratio
c misc:
        integer i,j
        double precision norm,CW2,Pi
        double precision
     &          g2hjbb(3),g2hjWW(3),g2hjZZ(3),                
     &          g2hjgg(3),g2hjhiZ_nHbynH(3,3)

        Pi = 3.1415926535897932384626433832795029D0

*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
c The number of neutral Higgs bosons in the MSSM is 3, therefore set
        nHiggs=3

c The string 'whichanalyses' determines which subset of experimental 
c results are used.
c In this example, we've used the option 'LandT',
c which instructs HiggsBounds to use tables of results
c from both LEP and the Tevatron (i.e. the full set of
c results supplied with HiggsBounds).
        whichanalyses='LandT'

c The subroutine initialize_HiggsBounds reads in all necessary
c tables etc.
c It must be called before run_HiggsBounds_part subroutine.

        print*,'calling initialize_HiggsBounds'

        call initialize_HiggsBounds(nHiggs,whichanalyses)

c If you would like to perform scans over variables, the subroutine
c initialize_HiggsBounds (and finish_HiggsBounds) should be called
c outside the do-loops in order to save time.
*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
c calls to FH subroutines:

c these subroutines are also contained in this file
c (they will call FHSetFlags, FHSetPara and FHGetPara)

        print*,'calling FHSetFlags'
	call setFlags
        print*,'calling FHSetPara'
	call setPara
        print*,'calling FHGetPara'
	call getPara

        print*,'calling FHHiggsCorr'

	call FHHiggsCorr(error, MHiggs, SAeff, UHiggs, ZHiggs)
	if( error .ne. 0 ) stop

c NOTE: we are setting uzint=uzext
        uzint=2
        uzext=2
        call FHSelectUZ(error, uzint, uzext)
	if( error .ne. 0 ) stop

	call FHCouplings(error,
     &    couplings, couplingsms, gammas, gammasms, excl)
	if( error .ne. 0 ) stop
	
        call FHHiggsProd(error, prodxs)
	if( error .ne. 0 ) stop 

        call FHRetrieveSMPara(error,
     &    invAlfa, AlfasMZ, GF,
     &    MS, MC, MB, MW, MZ,
     &    CKMlambda, CKMA, CKMrho, CKMeta)

*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *      
c Set variables needed by HiggsBounds (using results from FeynHiggs).
c See HiggsBounds documentation for definition of variables used
c as arguments to run_HiggsBounds_part and FeynHiggs 
c documentation for all other variables.

c Note: It is slightly more accurate to use the subroutine run_HiggsBounds_part 
c rather than the subroutine run_HiggsBounds_effC because the SM branching ratios
c used internally in HiggsBounds (from HDecay) are not identical to the SM branching
c ratios used in FeynHiggs

        do i=1,3

         Mh(i)=MHiggs(i) 

         BR_hjbb(i)     = BR(H0FF(i,4,3))
         BR_hjtautau(i) = BR(H0FF(i,2,3))
         BR_hjWW(i)     = BR(H0VV(i,4))
         BR_hjgaga(i)   = BR(H0VV(i,1))

         if(GammaSM(H0FF(i,4,3)).le.0.0D0)then
          g2hjbb(i)=0.0D0
         else
          g2hjbb(i)=Gamma(H0FF(i,4,3))/GammaSM(H0FF(i,4,3))
         endif

c Note that this is currently equivalent to
c         g2hjbb(i)= bbhTev(i)/bbhTevSM(i)
c         g2hjbb(i)= bbhLHC(i)/bbhLHCSM(i)
c         g2hjbb(i)= btagbhTev(i)/btagbhTevSM(i)
c         g2hjbb(i)= btagbhLHC(i)/btagbhLHCSM(i)
c as long as MH>80 GeV

         CS_tev_bg_hjb_ratio(i) = g2hjbb(i)
         CS_tev_bb_hj_ratio(i)  = g2hjbb(i)

         g2hjWW(i)= dble(  Coupling(H0VV(i,4)) 
     &                   / CouplingSM(H0VV(i,4)) )**2.0D0
     &            + dimag( Coupling(H0VV(i,4)) 
     &                   / CouplingSM(H0VV(i,4)) )**2.0D0
c Note that this is currently equivalent to
c         g2hjWW(i)= WhTev(i)/WhTevSM(i)
c         g2hjWW(i)= WhLHC(i)/WhLHCSM(i)
c	  g2hjWW(i)= qqhTev(i)/qqhTevSM(i)
c	  g2hjWW(i)= qqhLHC(i)/qqhLHCSM(i)
c as long as MH>80 GeV and uzint=uzext

         g2hjZZ(i)= dble(  Coupling(H0VV(i,3)) 
     &                   / CouplingSM(H0VV(i,3)) )**2.0D0
     &            + dimag( Coupling(H0VV(i,3)) 
     &                   / CouplingSM(H0VV(i,3)) )**2.0D0
c Note that this is currently equivalent to
c         g2hjZZ(i)= ZhTev(i)/ZhTevSM(i)
c         g2hjZZ(i)= ZhLHC(i)/ZhLHCSM(i)
c as long as MH>80 GeV and uzint=uzext
c It is also equivalent to g2hjWW(i)

         CS_lep_hjZ_ratio(i)        = g2hjZZ(i)
         CS_tev_dd_hjZ_ratio(i)     = g2hjZZ(i)
         CS_tev_uu_hjZ_ratio(i)     = g2hjZZ(i)
         CS_tev_ss_hjZ_ratio(i)     = g2hjZZ(i)
         CS_tev_cc_hjZ_ratio(i)     = g2hjZZ(i)
         CS_tev_bb_hjZ_ratio(i)     = g2hjZZ(i)
         CS_tev_ud_hjWp_ratio(i)    = g2hjZZ(i)
         CS_tev_cs_hjWp_ratio(i)    = g2hjZZ(i)
         CS_tev_ud_hjWm_ratio(i)    = g2hjZZ(i)
         CS_tev_cs_hjWm_ratio(i)    = g2hjZZ(i)
         CS_tev_pp_vbf_ratio(i)     = g2hjZZ(i)
         
         if(Mh(i).gt.80.0D0)then
          g2hjgg(i)= gghTev(i)/gghTevSM(i)
         else
c Currently, this result will only be used in conjunction 
c with data from hep-ex/0803.1514v1, which has a lower
c limit of 70 GeV.
          if(GammaSM(H0VV(i,5)).le.0.0D0)then
           g2hjgg(i)=0.0D0
          else
           g2hjgg(i)= Gamma(H0VV(i,5))/GammaSM(H0VV(i,5))
          endif
         endif

c Note that gghTev(i)/gghTevSM(i) is *not* currently equivalent to 
c         g2hjgg(i)= Gamma(H0VV(i,5))/GammaSM(H0VV(i,5))
c or
c         g2hjgg(i)= gghLHC(i)/gghLHCSM(i)
c since the k-factors are different

          CS_tev_gg_hj_ratio(i) = g2hjgg(i)

        enddo
  
	if( invAlfa .gt. 0 ) then
          CW2 = (MW/MZ)**2.0D0
          norm = 1/invAlfa*Pi/CW2/(1.0D0-CW2)
	else
          norm=GF*sqrt(2.0D0)*MZ**2.0D0
	endif

        do j=1,3
          do i=1,3   
            g2hjhiZ_nHbynH(j,i)= (
     &          dble(  Coupling(H0HV(j,i)) )**2.0D0
     &        + dimag( Coupling(H0HV(j,i)) )**2.0D0
     &          )
     &       /norm
            
            CS_lep_hjhi_ratio_nHbynH(j,i) = g2hjhiZ_nHbynH(j,i)

            BR_hjhihi_nHbynH(j,i)=BR(H0HH(j,i,i))
          enddo
        enddo
       
*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
c calls to run_HiggsBounds_part subroutine

        print*,'calling run_HiggsBounds_part'
	
        call run_HiggsBounds_part(nHiggs,Mh,                  
     &          CS_lep_hjZ_ratio,                             
     &          CS_lep_hjhi_ratio_nHbynH,                     
     &          CS_tev_gg_hj_ratio,CS_tev_bb_hj_ratio,        
     &          CS_tev_bg_hjb_ratio,                          
     &          CS_tev_ud_hjWp_ratio,CS_tev_cs_hjWp_ratio,     
     &          CS_tev_ud_hjWm_ratio,CS_tev_cs_hjWm_ratio,     
     &          CS_tev_dd_hjZ_ratio,CS_tev_uu_hjZ_ratio,      
     &          CS_tev_ss_hjZ_ratio,CS_tev_cc_hjZ_ratio,      
     &          CS_tev_bb_hjZ_ratio,                          
     &          CS_tev_pp_vbf_ratio,                          
     &          BR_hjbb,BR_hjtautau,                          
     &          BR_hjWW,BR_hjgaga,                             
     &          BR_hjhihi_nHbynH,                             
     &          HBresult,chan,                                
     &          obsratio, ncombined                          )



        write(*,*)' '      
        write(*,*)'*************    HiggsBounds Results  **************'
        write(*,*)' '    
        write(*,*)'Is this parameter point excluded by either LEP'
        write(*,*)'or Tevatron data?'         
        write(*,*) HBresult, ',  where'
        write(*,*)'               0 = yes, it is excluded'
        write(*,*)'               1 = no, it has not been excluded'
        write(*,*)'              -1 = invalid parameter set'    
        write(*,*)' '    
        write(*,*)'The process with the highest statistical sensitivity'
        write(*,*)'is'
        write(*,*) chan,'(see Key.dat)'
        write(*,*)'This process has a theoretical rate vs. limit of'
        write(*,*) obsratio
        write(*,*)' '    
        write(*,*)'The number of Higgs which have contributed to the'
        write(*,*)'theoretical rate of this process was'
        write(*,*) ncombined
        write(*,*)' '    
        write(*,*)'See HiggsBounds documentation for more information.'
        write(*,*)'****************************************************'
        write(*,*)' '    
       
*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
c deallocates arrays used by HiggsBounds:

        call finish_HiggsBounds
        
	end


************************************************************************

	subroutine setFlags
	implicit none

	integer mssmpart, fieldren, tanbren, higgsmix, p2approx
	integer looplevel, runningMT, tlBotResum, tlCplxApprox
	parameter (mssmpart = 4)
	parameter (fieldren = 0)
	parameter (tanbren = 0)
	parameter (higgsmix = 2)
	parameter (p2approx = 0)
	parameter (looplevel = 2)
	parameter (runningMT = 3)
	parameter (tlBotResum = 1)
	parameter (tlCplxApprox = 3)

	integer error

	call FHSetFlags(error, mssmpart, fieldren, tanbren,
     &    higgsmix, p2approx, looplevel,
     &    runningMT, tlBotResum, tlCplxApprox)
	if( error .ne. 0 ) stop
	end

************************************************************************

	subroutine setPara
	implicit none

	double precision invAlfa, AlfasMZ, GF
	double precision MS, MC, MB, MZ, MW
	double precision CKMlambda, CKMA, CKMrho, CKMeta
	parameter (invAlfa = -1)
	parameter (AlfasMZ = -1)
	parameter (GF = -1)
	parameter (MS = -1)
	parameter (MC = -1)
	parameter (MB = -1)
	parameter (MW = -1)
	parameter (MZ = -1)
	parameter (CKMlambda = -1)
	parameter (CKMA = -1)
	parameter (CKMrho = -1)
	parameter (CKMeta = -1)

	double precision MT, TB, MA0, MHp
	parameter (MT = 170.9D0)
	parameter (TB = 5)
	parameter (MA0 = 250)
	parameter (MHp = -1)

	double precision MSusy
	double precision M3SL, M2SL, M1SL
	double precision M3SE, M2SE, M1SE
	double precision M3SQ, M2SQ, M1SQ
	double precision M3SU, M2SU, M1SU
	double precision M3SD, M2SD, M1SD
	parameter (MSusy = 1000)
	parameter (M3SL = MSusy)
	parameter (M2SL = M3SL)
	parameter (M1SL = M2SL)
	parameter (M3SE = MSusy)
	parameter (M2SE = M3SE)
	parameter (M1SE = M2SE)
	parameter (M3SQ = MSusy)
	parameter (M2SQ = M3SQ)
	parameter (M1SQ = M2SQ)
	parameter (M3SU = MSusy)
	parameter (M2SU = M3SU)
	parameter (M1SU = M2SU)
	parameter (M3SD = MSusy)
	parameter (M2SD = M3SD)
	parameter (M1SD = M2SD)

	double complex Atau, At, Ab
	double complex Amu, Ac, As
	double complex Ae, Au, Ad
	parameter (At = 2000)
	parameter (Ab = At)
	parameter (Atau = At)
	parameter (Ac = At)
	parameter (As = Ab)
	parameter (Amu = Atau)
	parameter (Au = Ac)
	parameter (Ad = As)
	parameter (Ae = Amu)

	double complex MUE, M_1, M_2, M_3
	parameter (MUE = 200)
	parameter (M_1 = 0)
	parameter (M_2 = 500)
	parameter (M_3 = 800)

	double precision Qtau, Qt, Qb
	parameter (Qtau = 0)
	parameter (Qt = 0)
	parameter (Qb = 0)

	double precision scalefactor
	parameter (scalefactor = 1)

	integer error

	call FHSetSMPara(error,
     &    invAlfa, AlfasMZ, GF,
     &    MS, MC, MB, MW, MZ,
     &    CKMlambda, CKMA, CKMrho, CKMeta)
	if( error .ne. 0 ) stop

	call FHSetPara(error, scalefactor,
     &    MT, TB, MA0, MHp,
     &    M3SL, M3SE, M3SQ, M3SU, M3SD,
     &    M2SL, M2SE, M2SQ, M2SU, M2SD,
     &    M1SL, M1SE, M1SQ, M1SU, M1SD,
     &    MUE,
     &    Atau, At, Ab,
     &    Amu, Ac, As,
     &    Ae, Au, Ad,
     &    M_1, M_2, M_3,
     &    Qtau, Qt, Qb)
        if( error .ne. 0 ) stop
	end

************************************************************************

	subroutine setSLHA(filename)
	implicit none
	character*(*) filename

#include "SLHA.h"

	double precision scalefactor
	parameter (scalefactor = 1)

	integer error
	double complex slhadata(nslhadata)

	call SLHARead(error, slhadata, filename, 1)
	if( error .ne. 0 ) stop

	call FHSetSLHA(error, slhadata, scalefactor)
	if( error .ne. 0 ) stop
	end

************************************************************************

	subroutine getPara
	implicit none

	integer nmfv
	double precision MSf(2,4,3), MASf(6,2), MCha(2), MNeu(4)
	double complex USf(2,2,4,3), UASf(6,6,2)
	double complex UCha(2,2), VCha(2,2), ZNeu(4,4)
	double complex DeltaMB
	double precision MGl

	integer error

	call FHGetPara(error, nmfv, MSf, USf, MASf, UASf,
     &    MCha, UCha, VCha, MNeu, ZNeu, DeltaMB, MGl)
	if( error .ne. 0 ) stop

** print some sample output:
*	print *, "MCha =", MCha
*	print *, "MNeu =", MNeu
	end


