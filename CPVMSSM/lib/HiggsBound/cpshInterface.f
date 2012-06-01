      subroutine applyHiggBounds(HMASS_H, GAMBRN, MHC_H, NHC_H, MBMT_H,
     .  HBresult, obsratio,chan,ncombined)
************************************************************************
* CPsuperH2 with Higgsbounds
* adapted from HBwithCPsuperH by SK 
* 19 June 2009
************************************************************************
c!      IMPLICIT REAL*8(A-H,M,O-Z)
      IMPLICIT NONE
* HB input:
      REAL*8     HMASS_H(3)
      REAL*8 GAMBRN(101,3,3)   ! 101 = IFLAG_H(20)+IFLAG_H(21)+1 =NMNH
*                                      ISMN       =ISUSYN        = 50
      COMPLEX*16 NHC_H(100,3)  ! 100 = NCMAX

      real*8 MHC_H,MBMT_H
     
* HB output:
      integer chan,ncombined  
      real*8 obsratio
      integer HBresult
************************************
      
*-----------------------------------------------------------------------
*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *   *    
* used by initialize_HiggsBounds and run_HiggsBounds_part
* HB input:
        integer nHiggs        
        character*5 whichexperiment
        double precision Mh(3),                            
     &          CS_lep_hjZ_ratio(3),                              
     &          CS_lep_hjhi_ratio_nHbynH(3,3),                  
     &          CS_tev_gg_hj_ratio_dummy(3),
     &          CS_tev_bb_hj_ratio(3),   
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


* misc:
        integer i,j
        double precision betasq
        double precision
     &          g2hjVV(3),g2hjbb(3),               
     &          g2hjhiZ_nHbynH(3,3)

c The number of neutral Higgs bosons in the MSSM is 3, therefore set
        nHiggs=3
c The string 'whichexperiment' determines which subset of experimental 
c results are used.
c In this example, we've used the option 'onlyL',
c which instructs HiggsBounds to use tables of results
c from LEP only (i.e. no Tevatron results).
        whichexperiment='onlyL'

c The subroutine initialize_HiggsBounds reads in all necessary
c tables etc.
c It must be called before calling the run_HiggsBounds_part subroutine.

        call initialize_HiggsBounds(nHiggs,whichexperiment) 

c If you would like to perform scans over variables, the subroutine
c initialize_HiggsBounds (and finish_HiggsBounds) should be called
c outside the do-loops in order to save time.
*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
*-----------------------------------------------------------------------
* back to cpsuperh
*=======================================================================
*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *      
c Set variables needed by HiggsBounds (using results from CPsuperH).
c See HiggsBounds documentation for definition of variables used
c as arguments to run_HiggsBounds_part and CPsuperH 
c documentation for all other variables.

c Note: It is slightly more accurate to use the subroutine run_HiggsBounds_part 
c rather than the subroutine run_HiggsBounds_effC because the SM branching ratios
c used internally in HiggsBounds (from HDecay) are not identical to the SM branching
c ratios used in CPsuperH

        do i=1,3

         Mh(i)=HMASS_H(i)         
         BR_hjbb(i)     = GAMBRN(6,2,i)
         BR_hjtautau(i) = GAMBRN(3,2,i) 
         BR_hjWW(i)     = GAMBRN(10,2,i) 
         BR_hjgaga(i)   = GAMBRN(17,2,i) 
         
         betasq=1.0D0-4.0D0*(MBMT_H/HMASS_H(i))**2.0D0
         g2hjbb(i)=     
     &      abs(NHC_H(17,i))**2.0D0   
     &    + abs(NHC_H(18,i))**2.0D0/betasq

         CS_tev_bg_hjb_ratio(i) = g2hjbb(i)
         CS_tev_bb_hj_ratio(i)  = g2hjbb(i)

         g2hjVV(i)= abs(NHC_H(70,i))**2.0D0

         CS_lep_hjZ_ratio(i)        = g2hjVV(i)
         CS_tev_dd_hjZ_ratio(i)     = g2hjVV(i)
         CS_tev_uu_hjZ_ratio(i)     = g2hjVV(i)
         CS_tev_ss_hjZ_ratio(i)     = g2hjVV(i)
         CS_tev_cc_hjZ_ratio(i)     = g2hjVV(i)
         CS_tev_bb_hjZ_ratio(i)     = g2hjVV(i)
         CS_tev_ud_hjWp_ratio(i)    = g2hjVV(i)
         CS_tev_cs_hjWp_ratio(i)    = g2hjVV(i)
         CS_tev_ud_hjWm_ratio(i)    = g2hjVV(i)
         CS_tev_cs_hjWm_ratio(i)    = g2hjVV(i)
         CS_tev_pp_vbf_ratio(i)     = g2hjVV(i)
c ------------------------------------------------------------------
c We set CS_tev_gg_hj_ratio to zero because it will not be required 
c by the option 'onlyL'
         CS_tev_gg_hj_ratio_dummy(i) = 0.0D0
c ------------------------------------------------------------------

         BR_hjhihi_nHbynH(i,1)=GAMBRN(14,3,i)
         BR_hjhihi_nHbynH(i,2)=GAMBRN(16,3,i)         
         BR_hjhihi_nHbynH(i,3)=0.0D0

        enddo
        
        do j=1,3    
         do i=1,3 
          if(i.lt.j)then          
           g2hjhiZ_nHbynH(j,i)=g2hjVV(6-j-i)
           g2hjhiZ_nHbynH(i,j)=g2hjhiZ_nHbynH(j,i)
          else
           g2hjhiZ_nHbynH(j,i)=0.0D0           
          endif    
         enddo
        enddo 

        do j=1,3    
         do i=1,3        
          CS_lep_hjhi_ratio_nHbynH(j,i) = g2hjhiZ_nHbynH(j,i)
         enddo
        enddo 

*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
c call to run_HiggsBounds_part subroutine:

        call run_HiggsBounds_part(nHiggs,Mh,                  
     &          CS_lep_hjZ_ratio,                             
     &          CS_lep_hjhi_ratio_nHbynH,                     
     &          CS_tev_gg_hj_ratio_dummy,CS_tev_bb_hj_ratio,        
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
     &          obsratio, ncombined )
*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
*=======================================================================
*
      END

