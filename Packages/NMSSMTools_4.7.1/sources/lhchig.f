      SUBROUTINE LHCHIG(PAR,PROB)

*   Subroutine to check LHC constraints

      IMPLICIT NONE

      INTEGER I,J

      DOUBLE PRECISION PAR(*),PROB(*),SIG(3,10)
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2),CMASS
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5),BRSS(5)
      DOUBLE PRECISION BRCC(5),BRBB(5),BRTT(5),BRWW(3),BRZZ(3)
      DOUBLE PRECISION BRGG(5),BRZG(5),BRHHH(4),BRHAA(3,3)
      DOUBLE PRECISION BRHCHC(3),BRHAZ(3,2),BRAHA(3),BRAHZ(2,3)
      DOUBLE PRECISION BRHCW(5),BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5)
      DOUBLE PRECISION CU(5),CD(5),CV(3),CJ(5),CG(5)
      DOUBLE PRECISION BRJJSM,BREESM,BRMMSM,BRLLSM,BRSSSM,BRCCSM
      DOUBLE PRECISION BRBBSM,BRTTSM,BRWWSM,BRZZSM,BRGGSM,BRZGSM
      DOUBLE PRECISION brtopbw,brtopbh,brtopneutrstop(5,2),LHC_TBH
      DOUBLE PRECISION HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC
      DOUBLE PRECISION HCBRBT,HCBRWH(5),HCBRWHT,HCBRNC(5,2)
      DOUBLE PRECISION HCBRSQ(5),HCBRSL(3),HCBRSUSY,HCWIDTH

      COMMON/BRN/BRJJ,BREE,BRMM,BRLL,BRSS,BRCC,BRBB,BRTT,BRWW,
     .      BRZZ,BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     .      BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     .      BRSUSY,WIDTH
      COMMON/REDCOUP/CU,CD,CV,CJ,CG
      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/BR_top2body/brtopbw,brtopbh,brtopneutrstop
      COMMON/BRC/HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC,
     .       HCBRBT,HCBRWH,HCBRWHT,HCBRNC,HCBRSQ,HCBRSL,
     .       HCBRSUSY,HCWIDTH
      COMMON/LHCSIG/SIG

* Loop over H1, H2, H3

      DO I=1,3


       DO J=1,10
        SIG(I,J)=0d0
       ENDDO

       CALL HDECAY(SMASS(I),BRJJSM,BREESM,BRMMSM,BRLLSM,BRSSSM,
     .      BRCCSM,BRBBSM,BRTTSM,BRWWSM,BRZZSM,BRGGSM,BRZGSM)

*   H -> tautau
* VBF/VH
       IF(BRLLSM.NE.0d0)SIG(I,1)=CV(I)**2*BRLL(I)/BRLLSM
* ggF
       IF(BRLLSM.NE.0d0)SIG(I,2)=CJ(I)**2*BRLL(I)/BRLLSM
       
*   H -> bb
* VBF/VH
       IF(BRBBSM.NE.0d0)SIG(I,3)=CV(I)**2*BRBB(I)/BRBBSM
* ttH
       IF(BRBBSM.NE.0d0)SIG(I,4)=CU(I)**2*BRBB(I)/BRBBSM

*   H -> ZZ/WW
* VBF/VH
       IF(BRZZSM.NE.0d0)SIG(I,5)=CV(I)**2*BRZZ(I)/BRZZSM
* ggF
       IF(BRZZSM.NE.0d0)SIG(I,6)=CJ(I)**2*BRZZ(I)/BRZZSM
       
*   H -> gammagamma
* VBF/VH
       IF(BRGGSM.NE.0d0)SIG(I,7)=CV(I)**2*BRGG(I)/BRGGSM
* ggF
       IF(BRGGSM.NE.0d0)SIG(I,8)=CJ(I)**2*BRGG(I)/BRGGSM

*   H -> invisible
* VBF/VH
       SIG(I,9)=CV(I)**2*BRNEU(I,1,1)
* ggF
       SIG(I,10)=CJ(I)**2*BRNEU(I,1,1)

      ENDDO

* Bound on Br(t->bH+)*BR(H+->tau nu)

      PROB(45)=DDIM(brtopbh*HCBRL/LHC_TBH(CMASS),1d0)

      END


      DOUBLE PRECISION FUNCTION LHC_TBH(M)

* ATLAS constraints on BR(t->bH+)*BR(H+->taunu), ATLAS-CONF-2011-151 tab.5

      IMPLICIT NONE
      INTEGER I,N
      PARAMETER(N=8)
      DOUBLE PRECISION X(N),Y(N),M

      DATA X/90d0,100d0,110d0,120d0,130d0,140d0,150d0,160d0/ 
      DATA Y/.104d0,.098d0,.095d0,.077d0,.066d0,.071d0,.052d0,.141d0/ 

      LHC_TBH=1d9
      DO I=1,N-1
       IF((M.GE.X(I)).AND.(M.LE.X(I+1)))THEN
        LHC_TBH=(Y(I)+(Y(I+1)-Y(I))*(M-X(I))/(X(I+1)-X(I)))
        RETURN
       ENDIF
      ENDDO

      END


      SUBROUTINE Higgs_CHI2(PAR,PROB)

*      PROB(46) =/= 0  No Higgs in the MHmin-MHmax GeV range
*      PROB(47) =/= 0  chi2gam > chi2max
*      PROB(48) =/= 0  chi2bb > chi2max
*      PROB(49) =/= 0  chi2zz > chi2max

      IMPLICIT NONE
      INTEGER GMUFLAG,HFLAG
      DOUBLE PRECISION PAR(*),PROB(*),SIG(3,10),D1,D2
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2),CMASS
      DOUBLE PRECISION chi2gam,chi2bb,chi2zz,chi2max,MHmin,MHmax
      DOUBLE PRECISION agg,bgg,cgg,mugcengg,muvcengg
      DOUBLE PRECISION abb,bbb,cbb,mugcenbb,muvcenbb
      DOUBLE PRECISION azz,bzz,czz,mugcenzz,muvcenzz

      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/LHCSIG/SIG
      COMMON/HIGGSFIT/MHmin,MHmax,chi2max,chi2gam,chi2bb,chi2zz
      COMMON/GMUFLAG/GMUFLAG,HFLAG

* adding linearly 1 GeV exp. + 2 GeV theor. errors
      MHmin=125.1d0-3d0
      MHmax=125.1d0+3d0
      chi2max=6.18d0

* From 1409.1588:

c Chi^2 from gammagamma:
      agg=18.26d0
      bgg=2.84d0
      cgg=5.08d0
      mugcengg=1.25d0
      muvcengg=1.09d0

c Chi^2 from bb/tautau: (For ggF: use SIG(1,2) from H->tautau)
      abb=6.29d0
      bbb=2.62d0
      cbb=14.86d0
      mugcenbb=0.83d0
      muvcenbb=1.40d0

c Chi^2 from ZZ/WW:
      azz=39.07d0
      bzz=4.68d0
      czz=6.52d0
      mugcenzz=1.03d0
      muvcenzz=1.12d0

      IF(HFLAG.EQ.2)THEN
       D1=1d99
      ELSE
       D1=DDIM(SMASS(1)/MHMAX,1d0)-DDIM(1d0,SMASS(1)/MHMIN)
      ENDIF
      IF(HFLAG.EQ.1)THEN
       D2=1d99
      ELSE
       D2=DDIM(SMASS(2)/MHMAX,1d0)-DDIM(1d0,SMASS(2)/MHMIN)
      ENDIF

      IF(D1.EQ.0d0 .and. D2.EQ.0d0)THEN
       IF(DABS(SMASS(1)-SMASS(2)).LE.3.D0) THEN
        chi2gam=agg*(SIG(1,8)+SIG(2,8)-mugcengg)**2 
     .     +cgg*(SIG(1,7)+SIG(2,7)-muvcengg)**2
     .     +2d0*bgg*(SIG(1,8)+SIG(2,8)-mugcengg)
     .       *(SIG(1,7)+SIG(2,7)-muvcengg)
        chi2zz=azz*(SIG(1,6)+SIG(2,6)-mugcenzz)**2 
     .    +czz*(SIG(1,5)+SIG(2,5)-muvcenzz)**2
     .    +2d0*bzz*(SIG(1,6)+SIG(2,6)-mugcenzz)
     .       *(SIG(1,5)+SIG(2,5)-muvcenzz)
       ELSE
        chi2gam=MIN(
     .    agg*(SIG(1,8)-mugcengg)**2 +cgg*(SIG(1,7)-muvcengg)**2
     .    +2d0*bgg*(SIG(1,8)-mugcengg)*(SIG(1,7)-muvcengg),
     .    agg*(SIG(2,8)-mugcengg)**2+cgg*(SIG(2,7)-muvcengg)**2
     .    +2d0*bgg*(SIG(2,8)-mugcengg)*(SIG(2,7)-muvcengg))
        chi2zz=MIN(
     .    azz*(SIG(1,6)-mugcenzz)**2+czz*(SIG(1,5)-muvcenzz)**2
     .    +2d0*bzz*(SIG(1,6)-mugcenzz)*(SIG(1,5)-muvcenzz),
     .    azz*(SIG(2,6)-mugcenzz)**2+czz*(SIG(2,5)-muvcenzz)**2
     .    +2d0*bzz*(SIG(2,6)-mugcenzz)*(SIG(2,5)-muvcenzz))
       ENDIF
       chi2bb=abb*(SIG(1,2)+SIG(2,2)-mugcenbb)**2 
     .    +cbb*(SIG(1,3)+SIG(2,3)-muvcenbb)**2
     .    +2d0*bbb*(SIG(1,2)+SIG(2,2)-mugcenbb)
     .       *(SIG(1,3)+SIG(2,3)-muvcenbb)
      ELSEIF(D1.EQ.0d0)THEN
       chi2gam=agg*(SIG(1,8)-mugcengg)**2 +cgg*(SIG(1,7)-muvcengg)**2
     .    +2d0*bgg*(SIG(1,8)-mugcengg)*(SIG(1,7)-muvcengg)
       chi2bb=abb*(SIG(1,2)-mugcenbb)**2+cbb*(SIG(1,3)-muvcenbb)**2
     .    +2d0*bbb*(SIG(1,2)-mugcenbb)*(SIG(1,3)-muvcenbb)
       chi2zz=azz*(SIG(1,6)-mugcenzz)**2+czz*(SIG(1,5)-muvcenzz)**2
     .    +2d0*bzz*(SIG(1,6)-mugcenzz)*(SIG(1,5)-muvcenzz)
      ELSEIF(D2.EQ.0d0)THEN
       chi2gam=agg*(SIG(2,8)-mugcengg)**2+cgg*(SIG(2,7)-muvcengg)**2
     .    +2d0*bgg*(SIG(2,8)-mugcengg)*(SIG(2,7)-muvcengg)
       chi2bb=abb*(SIG(2,2)-mugcenbb)**2+cbb*(SIG(2,3)-muvcenbb)**2
     .    +2d0*bbb*(SIG(2,2)-mugcenbb)*(SIG(2,3)-muvcenbb)
       chi2zz=azz*(SIG(2,6)-mugcenzz)**2+czz*(SIG(2,5)-muvcenzz)**2
     .    +2d0*bzz*(SIG(2,6)-mugcenzz)*(SIG(2,5)-muvcenzz)
      ELSE
       chi2gam=0d0
       chi2bb=0d0
       chi2zz=0d0
       IF(DABS(D1).LT.DABS(D2))THEN
        PROB(46)=DABS(D1)
       chi2gam=agg*(SIG(1,8)-mugcengg)**2 +cgg*(SIG(1,7)-muvcengg)**2
     .    +2d0*bgg*(SIG(1,8)-mugcengg)*(SIG(1,7)-muvcengg)
       chi2bb=abb*(SIG(1,2)-mugcenbb)**2+cbb*(SIG(1,3)-muvcenbb)**2
     .    +2d0*bbb*(SIG(1,2)-mugcenbb)*(SIG(1,3)-muvcenbb)
       chi2zz=azz*(SIG(1,6)-mugcenzz)**2+czz*(SIG(1,5)-muvcenzz)**2
     .    +2d0*bzz*(SIG(1,6)-mugcenzz)*(SIG(1,5)-muvcenzz)
       ELSE
        PROB(46)=DABS(D2)
       chi2gam=agg*(SIG(2,8)-mugcengg)**2+cgg*(SIG(2,7)-muvcengg)**2
     .    +2d0*bgg*(SIG(2,8)-mugcengg)*(SIG(2,7)-muvcengg)
       chi2bb=abb*(SIG(2,2)-mugcenbb)**2+cbb*(SIG(2,3)-muvcenbb)**2
     .    +2d0*bbb*(SIG(2,2)-mugcenbb)*(SIG(2,3)-muvcenbb)
       chi2zz=azz*(SIG(2,6)-mugcenzz)**2+czz*(SIG(2,5)-muvcenzz)**2
     .    +2d0*bzz*(SIG(2,6)-mugcenzz)*(SIG(2,5)-muvcenzz)
       ENDIF
      ENDIF

      PROB(47)=DDIM(chi2gam/chi2MAX,1d0)
      PROB(48)=DDIM(chi2bb/chi2MAX,1d0)
      PROB(49)=DDIM(chi2zz/chi2MAX,1d0)
      
      END
