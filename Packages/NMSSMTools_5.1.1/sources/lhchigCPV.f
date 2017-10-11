      SUBROUTINE LHC_HIGGS_CPV(PAR,PROB)

***********************************************************************
*   Subroutine to check Higgs LHC constraints:
*      PROB(45) =/= 0: Bound on Br(t->bH+)*BR(H+->tau nu)
*      PROB(46) =/= 0:  No Higgs in the MHmin-MHmax GeV range
*      PROB(47) =/= 0:  chi2gam > chi2max
*      PROB(48) =/= 0:  chi2bb > chi2max
*      PROB(49) =/= 0:  chi2zz > chi2max
*      PROB(51) =/= 0: excluded by ggF/bb->H/A->tautau
*      PROB(52) =/= 0: Excluded H_125->AA->4mu
***********************************************************************

      IMPLICIT NONE

      CHARACTER*256 FILENAME,EXPCON_PATH,catpath

      INTEGER I,I1,J,J1,NX,NX13,NY,JBAR,JBARbb,K,K1,K2,K3
      PARAMETER(NX=18,NX13=10,NY=6)

      DOUBLE PRECISION PAR(*),PROB(*)
      DOUBLE PRECISION HMAS(NX),XSM(NX),XSMbb(NX),LCMS(NX),LCMSbb(NX)
      DOUBLE PRECISION LATLAS(NX),LATLASbb(NX)
      DOUBLE PRECISION HMAS13(NX13),XSM13(NX13),XSMbb13(NX13)
      DOUBLE PRECISION LATLAS13(NX13),LATLASbb13(NX13)
      DOUBLE PRECISION MH(5),XSMH(5),LCMSH(5),SIG1(5),LATLASH(5)
      DOUBLE PRECISION XSMHbb(5),LCMSHbb(5),SIGbb(5),LATLASHbb(5)
      DOUBLE PRECISION DEL,SIGTOT,MBAR,LCMSMB,LATLASMB
      DOUBLE PRECISION SIGTOTbb,MBARbb,LCMSMBbb,LATLASMBbb
      DOUBLE PRECISION X1(NY),LSIGBR(NY),LIMIT,LHC_TBH
      DOUBLE PRECISION D1,D2,CJ2,BRHTOAA,MHcen,masstest
      DOUBLE PRECISION agg,bgg,cgg,mugcengg,muvcengg,SSIG(6)
      DOUBLE PRECISION abb,bbb,cbb,mugcenbb,muvcenbb
      DOUBLE PRECISION azz,bzz,czz,mugcenzz,muvcenzz
      DOUBLE PRECISION BRJJSM,BREESM,BRMMSM,BRLLSM,BRCCSM
      DOUBLE PRECISION BRBBSM,BRTTSM,BRWWSM,BRZZSM,BRGGSM,BRZGSM
      DOUBLE PRECISION ggHgg(100,2),dummy(100,4),SMXS(100,2)

      DOUBLE PRECISION MHC,XC(2,2),MH0(5),XH(5,5),MA2
      DOUBLE PRECISION CU(5),CUP(5),CD(5),CDP(5),CB(5),CBP(5),CJ(5),
     . CJP(5),CI(5),CG(5),CGP(5),CV(5),CZG(5),CZGP(5)
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5),
     . BRCC(5),BRBB(5),BRTT(5),BRWW(5),BRZZ(5),BRGG(5),BRZG(5)
      DOUBLE PRECISION BRHHH(5,10),BRHCHC(5),BRHAZ(5,4),BRHCW(5),
     . BRHIGGS(5)
      DOUBLE PRECISION BRNEU(5,5,5),BRCHA(5,3),BRHSQ(5,10),BRHSL(5,7),
     . BRSUSY(5)
      DOUBLE PRECISION brtopbw,brtopbh,brtopneutrstop(5,2)
      DOUBLE PRECISION HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC,HCBRBT
      DOUBLE PRECISION MHmin,MHmax,chi2max,chi2gam,chi2bb,chi2zz
      DOUBLE PRECISION SIG(5,10),muH2

      COMMON/HISPEC/MHC,XC,MH0,XH,MA2
      COMMON/HNSMCOUP/CU,CUP,CD,CDP,CB,CBP,CJ,CJP,CI,CG,CGP,CV,CZG,CZGP
      COMMON/HNSMBR/BRJJ,BREE,BRMM,BRLL,BRCC,BRBB,BRTT,BRWW,BRZZ,
     . BRGG,BRZG
      COMMON/HNHIBR/BRHHH,BRHCHC,BRHAZ,BRHCW,BRHIGGS
      COMMON/HNSUSYBR/BRNEU,BRCHA,BRHSQ,BRHSL,BRSUSY
      COMMON/BR_top2body/brtopbw,brtopbh,brtopneutrstop
      COMMON/HCSMBR/HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC,HCBRBT
      COMMON/HIGGSFIT/MHmin,MHmax,chi2max,chi2gam,chi2bb,chi2zz
      COMMON/LHCSIGCPV/SIG
      COMMON/HIGGSMS/muH2

***********************************************************************
* 8 TeV limits
      DATA HMAS/90d0,100d0,120d0,140d0,160d0,180d0,200d0,250d0,300d0,
     . 350d0,400d0,450d0,500d0,600d0,700d0,800d0,900d0,1000d0/

* SM Higgs ggF prod. cross sect. at 8 TeV from 
* https://twiki.cern.ch/twiki/bin/view/LHCPhysics/
* CERNYellowReportPageAt8TeV#gluon_gluon_Fusion_Process
      DATA XSM/36.32d0,29.68d0,20.8d0,15.42d0,11.96d0,8.98d0,7.081d0,
     . 4.783d0,3.594d0,3.401d0,2.921d0,2.002d0,1.283d0,.523d0,.229d0,
     . .1097d0,.0571d0,.032d0/

* SM Higgs bbH prod. cross sect. at 8 TeV from 
* https://twiki.cern.ch/twiki/bin/view/LHCPhysics/
* /CrossSectionsFigures#MSSM_WG_plots (estimated)
      DATA XSMbb/0.56d0,0.42d0,0.25d0,0.15d0,8.7d-2,5.1d-2,3.2d-2,
     . 1.2d-2,4.5d-3,2.8d-3,1.4d-3,8.2d-4,4.9d-4,2.6d-4,1.2d-4,
     . 5.2d-5,2.4d-5,1.2d-5/

* Upper limit on ggF->H->tautau (8 TeV) from CMS-PAS-HIG-13-021, Tab. 7
      DATA LCMS/50.2d0,31.3d0,7.38d0,2.27d0,.845d0,.549d0,.517d0,.315d0,
     . .15d0,.112d0,.103d0,.607d-1,.385d-1,.193d-1,.143d-1,.115d-1,
     . .923d-2,.865d-2/

* Upper limit on bbH->tautau (8 TeV) from CMS-PAS-HIG-13-021, Tab. 8
      DATA LCMSbb/6.03d0,4.14d0,1.76d0,1.25d0,.814d0,.659d0,.553d0,
     . .217d0,
     . .975d-1,.638d-1,.613d-1,.431d-1,.320d-1,.203d-1,.173d-1,.166d-1,
     . .146d-1,.133d-1/

* Upper limit on ggF->H->tautau (8 TeV) from ATLAS-CONF-2014-049, Fig. 7
      DATA LATLAS/29.1d0,24.0d0,5.25d0,2.02d0,1.39d0,1.00d0,.794d0,
     .  .281d0,.127d0,.112d0,.773d-1,.400d-1,.240d-1,.177d-1,.127d-1,
     .  .993d-2,.840d-2,.735d-2/

* Upper limit on bbH->tautau (8 TeV) from ATLAS-CONF-2014-049, Fig. 7
      DATA LATLASbb/6.32d0,6.32d0,2.73d0,1.27d0,.966d0,.606d0,.393d0,
     .  .305d0,.116d0,.101d0,.656d-1,.363d-1,.238d-1,.159d-1,.117d-1,
     .  .943d-2,.785d-2,.716d-2/

***********************************************************************
* 13 TeV limits
      DATA HMAS13/200d0,300d0,350d0,400d0,500d0,600d0,700d0,800d0,
     . 1000d0,1200d0/

* SM Higgs ggF prod. cross sect. at 13 TeV from 
* https://twiki.cern.ch/twiki/bin/view/LHCPhysics/
      DATA XSM13/1.812d1,9.823d0,1.025d1,1.000d1,4.538d0,2.006d0,
     . 9.235d-1,4.491d-1,1.233d-1,3.958d-2/

* SM Higgs bbH prod. cross sect. at 13 TeV from 
* https://twiki.cern.ch/twiki/bin/view/LHCPhysics/
      DATA XSMbb13/1.000d-1,2.180d-2,1.174d-2,6.731d-3,2.547d-3,
     . 1.101d-3,5.251d-4,2.692d-4,8.258d-5,2.950d-5/

* Upper limit on ggF->H->tautau (13 TeV) from ATLAS-CONF-2016-085, Tab. 4
      DATA LATLAS13/2.1d0,5.9d-1,2.1d-1,1.5d-1,6.3d-2,3.0d-2,
     . 2.2d-2,1.8d-2,1.5d-2,1.4d-2/

* Upper limit on bbH->tautau (13 TeV) from ATLAS-CONF-2016-085, Tab. 5
      DATA LATLASbb13/2.0d0,5.7d-1,3.3d-1,2.7d-1,1.4d-1,5.4d-2,
     . 3.6d-2,2.5d-2,1.6d-2,1.3d-2/

***********************************************************************

      DATA X1/.25d0,.5d0,.75d0,1d0,2d0,3d0/
      DATA LSIGBR/3.8d-5,3.6d-5,4.1d-5,4.2d-5,4.4d-5,4.6d-5/

***********************************************************************

c       I- Constraints from ggF/bb->H/A->tautau

* Loop over 5 Higgses
      DO I=1,5
        MH(I)=dsqrt(MH0(I))

* 8 TeV limits
        J=1
        DOWHILE(HMAS(J).LE.MH(I) .AND. J.LT.NX)
          J=J+1
        ENDDO
        IF(J.GE.2 .AND. MH(I).LT.HMAS(NX)) THEN

* SM Higgs ggF prod. cross sect.:
          XSMH(I)=XSM(J-1)+(MH(I)-HMAS(J-1))/
     .      (HMAS(J)-HMAS(J-1))*(XSM(J)-XSM(J-1))
* ggF Signal cross section*BR:
          SIG1(I)=(CJ(I)**2+CJP(I)**2)*BRLL(I)*XSMH(I)

* CMS ggF limit:
          LCMSH(I)=LCMS(J-1)+(MH(I)-HMAS(J-1))/
     .      (HMAS(J)-HMAS(J-1))*(LCMS(J)-LCMS(J-1))
          PROB(51)=PROB(51)+DDIM(1d0,LCMSH(I)/SIG1(I))

* ATLAS ggF limit:
          LATLASH(I)=LATLAS(J-1)+(MH(I)-HMAS(J-1))/
     .      (HMAS(J)-HMAS(J-1))*(LATLAS(J)-LATLAS(J-1))
* Correct for jump in Fig.7 at MA=200 GeV: J=8, 
* modif. LATLAS(J-1)=LATLAS(7)=.96D0 and not .794d0:
          IF(J.EQ.8) THEN
            LATLASH(I)=.96D0+(MH(I)-HMAS(J-1))/
     .        (HMAS(J)-HMAS(J-1))*(LATLAS(J)-.96D0)
          ENDIF
          PROB(51)=PROB(51)+DDIM(1D0,LATLASH(I)/SIG1(I))

* SM Higgs bbH prod. cross sect.:
          XSMHbb(I)=XSMbb(J-1)+(MH(I)-HMAS(J-1))/
     .      (HMAS(J)-HMAS(J-1))*(XSMbb(J)-XSMbb(J-1))
* bbH Signal cross section*BR:
          SIGbb(I)=(CB(I)**2+CBP(I)**2)*XSMHbb(I)*BRLL(I)

* CMS Hbb limit:
          LCMSHbb(I)=LCMSbb(J-1)+(MH(I)-HMAS(J-1))/
     .      (HMAS(J)-HMAS(J-1))*(LCMSbb(J)-LCMSbb(J-1))
          PROB(51)=PROB(51)+DDIM(1d0,LCMSHbb(I)/SIGbb(I))

* ATLAS Hbb limit:
          LATLASHbb(I)=LATLASbb(J-1)+(MH(I)-HMAS(J-1))/
     .      (HMAS(J)-HMAS(J-1))*(LATLASbb(J)-LATLASbb(J-1))
* Correct for jump in Fig.7 at MA=200 GeV: J=8, 
* modif. LATLASbb(J-1)=LATLASbb(7)=.858D0 and not .393d0:
          IF(J.EQ.8) THEN
            LATLASHbb(I)=.858d0+(MH(I)-HMAS(J-1))/
     .        (HMAS(J)-HMAS(J-1))*(LATLASbb(J)-.858d0)
          ENDIF
          PROB(51)=PROB(51)+DDIM(1d0,LATLASHbb(I)/SIGbb(I))

* Combine signal rates of 2 Higgses
          DO I1=1,I-1
            J1=1
            DOWHILE(HMAS(J1).LE.MH(I1) .AND. J1.LT.NX)
              J1=J1+1
            ENDDO
            IF(MH(I1).GE.HMAS(1) .AND. MH(I1).LT.HMAS(NX)) THEN

* Average masses weighted by the signal rates (MBAR for ggF, MBARbb for bbH):
             MBAR=(SIG1(I)*MH(I)+SIG1(I1)*MH(I1))/(SIG1(I)+SIG1(I1))
             JBAR=1
             DOWHILE(HMAS(JBAR).LE.MBAR.AND.JBAR.LT.NX)
               JBAR=JBAR+1
             ENDDO
             MBARbb=(SIGbb(I)*MH(I)
     .         +SIGbb(I1)*MH(I1))/(SIGbb(I)+SIGbb(I1))
             JBARbb=1
             DOWHILE(HMAS(JBARbb).LE.MBARbb.AND.JBARbb.LT.NX)
               JBARbb=JBARbb+1
             ENDDO
* DEL=mass difference divided by a (small) resolution squared:
* [DEL < 1 only if |MH(I)-MH(I1)| < (MH(I)+MH(I1))/15d0;
*  otherwise the combined signal rate is small]
             DEL=((MH(I)-MH(I1))/(MH(I)+MH(I1))*15d0)**2

* Estimate of the combined ggF signal rates:
             SIGTOT=SIG1(I)+SIG1(I1)
     .         -SIG1(I)*SIG1(I1)*DEL/(SIG1(I)+SIG1(I1))
* Continue only if SIGTOT > 0 and 90<MBAR<1000
*      and |MH(I)-MH(I1)|/MBAR<0.20:
             IF(SIGTOT.GT.0D0.AND.MBAR.GE.HMAS(1).AND.
     .         MBAR.LT.HMAS(NX).AND.
     .         dabs(MH(I)-MH(I1))/MBAR.LE.0.20d0) THEN

* CMS ggF limit at MBAR:
               LCMSMB=LCMS(JBAR-1)+(MBAR-HMAS(JBAR-1))/
     .          (HMAS(JBAR)-HMAS(JBAR-1))*(LCMS(JBAR)-LCMS(JBAR-1))
               PROB(51)=PROB(51)+DDIM(1d0,LCMSMB/SIGTOT)

* ATLAS ggF limit at MBAR:
               LATLASMB=LATLAS(JBAR-1)+(MBAR-HMAS(JBAR-1))/
     .          (HMAS(JBAR)-HMAS(JBAR-1))*(LATLAS(JBAR)-LATLAS(JBAR-1))
* Correct for jump in Fig.7 at MA=200 GeV: JBAR=8, 
* modif. LATLAS(JBAR-1)=LATLAS(7)=.96D0 and not .794d0:
               IF(J.EQ.8) THEN
                 LATLASMB=.96D0+(MBAR-HMAS(JBAR-1))/
     .            (HMAS(JBAR)-HMAS(JBAR-1))*(LATLAS(JBAR)-.96D0)
               ENDIF
               PROB(51)=PROB(51)+DDIM(1D0,LATLASMB/SIGTOT)
             ENDIF

* Estimate of the combined bbH signal rates:
             SIGTOTbb=SIGbb(I)+SIGbb(I1)
     .         -SIGbb(I)*SIGbb(I1)*DEL/(SIGbb(I)+SIGbb(I1))
* Continue only if SIGTOTbb > 0 and 90<MBARbb<1000
*      and |MH(I)-MH(I1)|/MBARbb<0.20:
             IF(SIGTOTbb.GT.0D0.AND.MBARbb.GE.HMAS(1).AND.
     .         MBARbb.LT.HMAS(NX).AND.
     .         dabs(MH(I)-MH(I1))/MBARbb.LE.0.20d0) THEN

* CMS bbH limit at MBARbb:
               LCMSMBbb=LCMSbb(JBARbb-1)+(MBARbb-HMAS(JBARbb-1))/
     .           (HMAS(JBARbb)-HMAS(JBARbb-1))*
     .           (LCMSbb(JBARbb)-LCMSbb(JBARbb-1))
               PROB(51)=PROB(51)+DDIM(1d0,LCMSMBbb/SIGTOTbb)

* ATLAS bbH limit at MBARbb:
               LATLASMBbb=LATLASbb(JBARbb-1)+(MBARbb-HMAS(JBARbb-1))/
     .           (HMAS(JBARbb)-HMAS(JBARbb-1))*
     .           (LATLASbb(JBARbb)-LATLASbb(JBARbb-1))
* Correct for jump in Fig.7 at MA=200 GeV: JBARbb=8, 
* modif. LATLASbb(JBARbb-1)=LATLASbb(7)=.858D0 and not .393d0:
               IF(J.EQ.8) THEN
                 LATLASMBbb=.858D0+(MBARbb-HMAS(JBARbb-1))/
     .             (HMAS(JBARbb)-HMAS(JBARbb-1))*
     .             (LATLASbb(JBARbb)-.858D0)
               ENDIF
               PROB(51)=PROB(51)+DDIM(1D0,LATLASMBbb/SIGTOTbb)
             ENDIF
            ENDIF
          ENDDO
        ENDIF

* 13 TeV limits
        J=1
        DOWHILE(HMAS13(J).LE.MH(I) .AND. J.LT.NX13)
          J=J+1
        ENDDO
        IF(J.GE.2 .AND. MH(I).LT.HMAS13(NX13)) THEN

* SM Higgs ggF prod. cross sect.:
          XSMH(I)=XSM13(J-1)+(MH(I)-HMAS13(J-1))/
     .      (HMAS13(J)-HMAS13(J-1))*(XSM13(J)-XSM13(J-1))
* ggF Signal cross section*BR:
          SIG1(I)=(CJ(I)**2+CJP(I)**2)*BRLL(I)*XSMH(I)

* ATLAS ggF limit:
          LATLASH(I)=LATLAS13(J-1)+(MH(I)-HMAS13(J-1))/
     .      (HMAS13(J)-HMAS13(J-1))*(LATLAS13(J)-LATLAS13(J-1))
          PROB(51)=PROB(51)+DDIM(1D0,LATLASH(I)/SIG1(I))

* SM Higgs bbH prod. cross sect.:
          XSMHbb(I)=XSMbb13(J-1)+(MH(I)-HMAS13(J-1))/
     .      (HMAS13(J)-HMAS13(J-1))*(XSMbb13(J)-XSMbb13(J-1))
* bbH Signal cross section*BR:
          SIGbb(I)=(CB(I)**2+CBP(I)**2)**2*XSMHbb(I)*BRLL(I)

* ATLAS Hbb limit:
          LATLASHbb(I)=LATLASbb13(J-1)+(MH(I)-HMAS13(J-1))/
     .      (HMAS13(J)-HMAS13(J-1))*(LATLASbb13(J)-LATLASbb13(J-1))
          PROB(51)=PROB(51)+DDIM(1d0,LATLASHbb(I)/SIGbb(I))

* Combine signal rates of 2 Higgses
          DO I1=1,I-1
            J1=1
            DOWHILE(HMAS13(J1).LE.MH(I1) .AND. J1.LT.NX13)
              J1=J1+1
            ENDDO
            IF(MH(I1).GE.HMAS13(1) .AND. MH(I1).LT.HMAS13(NX13))THEN

* Average masses weighted by the signal rates (MBAR for ggF, MBARbb for bbH):
             MBAR=(SIG1(I)*MH(I)+SIG1(I1)*MH(I1))/(SIG1(I)+SIG1(I1))
             JBAR=1
             DOWHILE(HMAS13(JBAR).LE.MBAR .AND. JBAR.LT.NX13)
               JBAR=JBAR+1
             ENDDO
             MBARbb=(SIGbb(I)*MH(I)
     .         +SIGbb(I1)*MH(I1))/(SIGbb(I)+SIGbb(I1))
             JBARbb=1
             DOWHILE(HMAS13(JBARbb).LE.MBARbb .AND. JBARbb.LT.NX13)
               JBARbb=JBARbb+1
             ENDDO
* DEL=mass difference divided by a (small) resolution squared:
* [DEL < 1 only if |MH(I)-MH(I1)| < (MH(I)+MH(I1))/15d0;
*  otherwise the combined signal rate is small]
             DEL=((MH(I)-MH(I1))/(MH(I)+MH(I1))*15d0)**2

* Estimate of the combined ggF signal rates:
             SIGTOT=SIG1(I)+SIG1(I1)
     .         -SIG1(I)*SIG1(I1)*DEL/(SIG1(I)+SIG1(I1))
* Continue only if SIGTOT > 0 and 200<MBAR<1200
*      and |MH(I)-MH(I1)|/MBAR<0.20:
             IF(SIGTOT.GT.0D0.AND.MBAR.GE.HMAS13(1).AND.
     .         MBAR.LT.HMAS13(NX13).AND.
     .         dabs(MH(I)-MH(I1))/MBAR.LE.0.20d0) THEN

* ATLAS ggF limit at MBAR:
               LATLASMB=LATLAS13(JBAR-1)+(MBAR-HMAS13(JBAR-1))/
     .          (HMAS13(JBAR)-HMAS13(JBAR-1))*
     .          (LATLAS13(JBAR)-LATLAS13(JBAR-1))
               PROB(51)=PROB(51)+DDIM(1D0,LATLASMB/SIGTOT)
             ENDIF

* Estimate of the combined bbH signal rates:
             SIGTOTbb=SIGbb(I)+SIGbb(I1)
     .         -SIGbb(I)*SIGbb(I1)*DEL/(SIGbb(I)+SIGbb(I1))
* Continue only if SIGTOTbb > 0 and 200<MBARbb<1200
*      and |MH(I)-MH(I1)|/MBARbb<0.20:
             IF(SIGTOTbb.GT.0D0.AND.MBARbb.GE.HMAS13(1).AND.
     .         MBARbb.LT.HMAS13(NX13).AND.
     .         dabs(MH(I)-MH(I1))/MBARbb.LE.0.20d0) THEN

* ATLAS bbH limit at MBARbb:
               LATLASMBbb=LATLASbb13(JBARbb-1)+
     .           (MBARbb-HMAS13(JBARbb-1))/
     .           (HMAS13(JBARbb)-HMAS13(JBARbb-1))*
     .           (LATLASbb13(JBARbb)-LATLASbb13(JBARbb-1))
               PROB(51)=PROB(51)+DDIM(1D0,LATLASMBbb/SIGTOTbb)
             ENDIF
            ENDIF
          ENDDO
        ENDIF
      ENDDO


c        II- Bound on Br(t->bH+)*BR(H+->tau nu)

      PROB(45)=DDIM(brtopbh*HCBRL/LHC_TBH(dsqrt(MHC)),1d0)


c        III- Constraints from ggF->H/A->gamgam from ATLAS-CONF-2014-031, M_H/A < 122

*   The EXPCON_PATH variable is set:
      CALL getenv('EXPCON_PATH',EXPCON_PATH)
      if(EXPCON_PATH.eq.' ')  EXPCON_PATH='../EXPCON'

* Read ATLAS upper limit
* ggHgg(I,1): Higgs mass
* ggHgg(I,2): upper limit in fb
      FILENAME=catpath(EXPCON_PATH,'ggHgg.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1801 READ(11,*,END=1802,ERR=2)(dummy(I,J),J=1,4)
      ggHgg(I,1)=dummy(I,1)
      ggHgg(I,2)=dummy(I,4)
      I=I+1
      GOTO 1801
 1802 CLOSE(11)

* Read SM Higgs ggF production cross section (60-122 GeV)
* SMXS(I,1): Higgs mass
* SMXS.65-122.dat:  ggF production cross section in pb
* SMXS(I,2): ggF production cross section in fb
      FILENAME=catpath(EXPCON_PATH,'SMXS.65-122.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1811 READ(11,*,END=1812,ERR=2)(SMXS(I,J),J=1,2)
      SMXS(I,2)=1d3*SMXS(I,2)
      I=I+1
      GOTO 1811
 1812 CLOSE(11)

* Loop over 5 Higgses
      DO I=1,5
        J=1
        DOWHILE(SMXS(J,1).LE.dsqrt(MH0(I)) .AND. J.LT.70)
          J=J+1
        ENDDO

       IF(J.GE.2 .AND. dsqrt(MH0(I)).LT.122d0) THEN
        XSMH(I)=0d0
        LATLASH(I)=1d10
* SM Higgs ggF prod. cross sect.:
          XSMH(I)=SMXS(J-1,2)+(dsqrt(MH0(I))-SMXS(J-1,1))/
     .      (SMXS(J,1)-SMXS(J-1,1))*(SMXS(J,2)-SMXS(J-1,2))
* ggF Signal cross section*BR(H->gamgam):
          SIG(I,8)=(CJ(I)**2+CJP(I)**2)*BRGG(I)*XSMH(I)
* ATLAS limit:
          LATLASH(I)=ggHgg(J-1,2)+(dsqrt(MH0(I))-SMXS(J-1,1))/
     .      (SMXS(J,1)-SMXS(J-1,1))*(ggHgg(J,2)-ggHgg(J-1,2))
          PROB(53)=PROB(53)+DDIM(1d0,LATLASH(I)/SIG(I,8))

        ENDIF
        ENDDO


c       IV- Higgs test at ~125 GeV

c      1) Data

* adding linearly 1 GeV exp. + 2 GeV theor. errors
      MHmin=DSQRT(muH2)-3d0
      MHmax=DSQRT(muH2)+3d0
      chi2max=6.18d0

* From J. Bernon with run-II data (sep. 2016):

c Chi^2 from gammagamma:
      agg=17.47d0
      bgg=3.17d0
      cgg=6.26d0
      mugcengg=1.18d0
      muvcengg=1.07d0

c Chi^2 from bb/tautau:
      abb=4.81d0
      bbb=2.69d0
      cbb=21.57d0
      mugcenbb=1.27d0
      muvcenbb=0.86d0

c Chi^2 from ZZ/WW:
      azz=36.40d0
      bzz=4.56d0
      czz=8.16d0
      mugcenzz=1.11d0
      muvcenzz=1.37d0

c      2) Higgs signals / SM

      DO I=1,5

       DO J=1,10
        SIG(I,J)=0d0
       ENDDO

       CALL HDECAY(dsqrt(MH0(I)),BRJJSM,BREESM,BRMMSM,BRLLSM,
     .      BRCCSM,BRBBSM,BRTTSM,BRWWSM,BRZZSM,BRGGSM,BRZGSM)

*   H -> tautau
* VBF/VH
       IF(BRLLSM.NE.0d0)SIG(I,1)=CV(I)**2*BRLL(I)/BRLLSM
* ggF
       IF(BRLLSM.NE.0d0)SIG(I,2)=(CJ(I)**2+CJP(I)**2)*BRLL(I)/BRLLSM
       
*   H -> bb
* VBF/VH
       IF(BRBBSM.NE.0d0)SIG(I,3)=CV(I)**2*BRBB(I)/BRBBSM
* ttH
       IF(BRBBSM.NE.0d0)SIG(I,4)=(CU(I)**2+CUP(I)**2)*BRBB(I)/BRBBSM

*   H -> ZZ/WW
* VBF/VH
       IF(BRZZSM.NE.0d0)SIG(I,5)=CV(I)**2*BRZZ(I)/BRZZSM
* ggF
       IF(BRZZSM.NE.0d0)SIG(I,6)=(CJ(I)**2+CJP(I)**2)*BRZZ(I)/BRZZSM
       
*   H -> gammagamma
* VBF/VH
       IF(BRGGSM.NE.0d0)SIG(I,7)=CV(I)**2*BRGG(I)/BRGGSM
* ggF
       IF(BRGGSM.NE.0d0)SIG(I,8)=(CJ(I)**2+CJP(I)**2)*BRGG(I)/BRGGSM

*   H -> invisible
* VBF/VH
       SIG(I,9)=CV(I)**2*BRNEU(I,1,1)
* ggF
       SIG(I,10)=(CJ(I)**2+CJP(I)**2)*BRNEU(I,1,1)

      ENDDO

c      3) Chi^2 test of the couplings

      MHcen=(MHmin+MHmax)/2d0
      D1=(MHmax-MHmin)/2d0
      masstest=2d0
      K=0
      K1=0
      K2=0
      K3=0
      J=1

      DO I=1,6
       SSIG(I)=0d0
      ENDDO

      DO I=1,5
      masstest=min(masstest,dabs(dsqrt(MH0(I))-MHcen)/D1)
       IF(dabs(dsqrt(MH0(I))-MHcen).le.D1)THEN
        if(K.eq.2)then
         K3=I
         K=3
        endif
        if(K.eq.1)then
         K2=I
         K=2
        endif
        if(K.eq.0)then
         K1=I
         K=1
        endif
       ELSE
        if(K.eq.0.and.I.gt.1)then
         if(dabs(dsqrt(MH0(I))-MHcen)
     .  .lt.dabs(dsqrt(MH0(I-1))-MHcen))J=I
        endif
       ENDIF
      ENDDO

      IF(K.eq.0)THEN
       chi2gam=agg*(SIG(J,1)-mugcengg)**2+cgg*(SIG(J,2)-muvcengg)**2
     .     +2d0*bgg*(SIG(J,1)-mugcengg)*(SIG(J,2)-muvcengg)
       chi2bb=abb*(SIG(J,3)-mugcenbb)**2+cbb*(SIG(J,4)-muvcenbb)**2
     .    +2d0*bbb*(SIG(J,3)-mugcenbb)*(SIG(J,4)-muvcenbb)
       chi2zz=azz*(SIG(J,5)-mugcenzz)**2+czz*(SIG(J,6)-muvcenzz)**2
     .    +2d0*bzz*(SIG(J,5)-mugcenzz)*(SIG(J,6)-muvcenzz)
      ENDIF

      IF(K.eq.1)THEN
       chi2gam=agg*(SIG(K1,1)-mugcengg)**2+cgg*(SIG(K1,2)-muvcengg)**2
     .     +2d0*bgg*(SIG(K1,1)-mugcengg)*(SIG(K1,2)-muvcengg)
       chi2bb=abb*(SIG(K1,3)-mugcenbb)**2+cbb*(SIG(K1,4)-muvcenbb)**2
     .    +2d0*bbb*(SIG(K1,3)-mugcenbb)*(SIG(K1,4)-muvcenbb)
       chi2zz=azz*(SIG(K1,5)-mugcenzz)**2+czz*(SIG(K1,6)-muvcenzz)**2
     .    +2d0*bzz*(SIG(K1,5)-mugcenzz)*(SIG(K1,6)-muvcenzz)
      ENDIF

      IF(K.eq.2)THEN
       chi2gam=min(agg*(SIG(K1,1)-mugcengg)**2
     .                +cgg*(SIG(K1,2)-muvcengg)**2
     .     +2d0*bgg*(SIG(K1,1)-mugcengg)*(SIG(K1,2)-muvcengg),
     .     agg*(SIG(K2,1)-mugcengg)**2+cgg*(SIG(K2,2)-muvcengg)**2
     .     +2d0*bgg*(SIG(K2,1)-mugcengg)*(SIG(K2,2)-muvcengg))
       chi2zz=min(azz*(SIG(K1,5)-mugcenzz)**2
     .                +czz*(SIG(K1,6)-muvcenzz)**2
     .    +2d0*bzz*(SIG(K1,5)-mugcenzz)*(SIG(K1,6)-muvcenzz),
     .      azz*(SIG(K2,5)-mugcenzz)**2+czz*(SIG(K2,6)-muvcenzz)**2
     .    +2d0*bzz*(SIG(K2,5)-mugcenzz)*(SIG(K2,6)-muvcenzz))

        SSIG(1)=SIG(K1,8)+SIG(K2,8)
        SSIG(2)=SIG(K1,7)+SIG(K2,7)
        SSIG(3)=SIG(K1,2)+SIG(K2,2)
        SSIG(4)=SIG(K1,3)+SIG(K2,3)
        SSIG(5)=SIG(K1,6)+SIG(K2,6)
        SSIG(6)=SIG(K1,5)+SIG(K2,5)

        if(dabs(dsqrt(MH0(K1))-dsqrt(MH0(K2))).le.3d0)then
       chi2gam=agg*(SSIG(1)-mugcengg)**2+cgg*(SSIG(2)-muvcengg)**2
     .     +2d0*bgg*(SSIG(1)-mugcengg)*(SSIG(2)-muvcengg)
       chi2zz=azz*(SSIG(5)-mugcenzz)**2+czz*(SSIG(6)-muvcenzz)**2
     .    +2d0*bzz*(SSIG(5)-mugcenzz)*(SSIG(6)-muvcenzz)
        endif

       chi2bb=abb*(SSIG(3)-mugcenbb)**2+cbb*(SSIG(4)-muvcenbb)**2
     .    +2d0*bbb*(SSIG(3)-mugcenbb)*(SSIG(4)-muvcenbb)
      ENDIF

      IF(K.eq.3)THEN
       chi2gam=min(agg*(SIG(K1,1)-mugcengg)**2
     .                 +cgg*(SIG(K1,2)-muvcengg)**2
     .     +2d0*bgg*(SIG(K1,1)-mugcengg)*(SIG(K1,2)-muvcengg),
     .     agg*(SIG(K2,1)-mugcengg)**2+cgg*(SIG(K2,2)-muvcengg)**2
     .     +2d0*bgg*(SIG(K2,1)-mugcengg)*(SIG(K2,2)-muvcengg),
     .     agg*(SIG(K3,1)-mugcengg)**2+cgg*(SIG(K3,2)-muvcengg)**2
     .     +2d0*bgg*(SIG(K3,1)-mugcengg)*(SIG(K3,2)-muvcengg))
       chi2zz=min(azz*(SIG(K1,5)-mugcenzz)**2
     .                 +czz*(SIG(K1,6)-muvcenzz)**2
     .    +2d0*bzz*(SIG(K1,5)-mugcenzz)*(SIG(K1,6)-muvcenzz),
     .      azz*(SIG(K2,5)-mugcenzz)**2+czz*(SIG(K2,6)-muvcenzz)**2
     .    +2d0*bzz*(SIG(K2,5)-mugcenzz)*(SIG(K2,6)-muvcenzz),
     .      azz*(SIG(K3,5)-mugcenzz)**2+czz*(SIG(K3,6)-muvcenzz)**2
     .    +2d0*bzz*(SIG(K3,5)-mugcenzz)*(SIG(K3,6)-muvcenzz))

        SSIG(1)=SIG(K1,8)+SIG(K2,8)+SIG(K3,8)
        SSIG(2)=SIG(K1,7)+SIG(K2,7)+SIG(K3,7)
        SSIG(3)=SIG(K1,2)+SIG(K2,2)+SIG(K3,2)
        SSIG(4)=SIG(K1,3)+SIG(K2,3)+SIG(K3,3)
        SSIG(5)=SIG(K1,6)+SIG(K2,6)+SIG(K3,6)
        SSIG(6)=SIG(K1,5)+SIG(K2,5)+SIG(K3,5)

        if(dabs(dsqrt(MH0(K1))-dsqrt(MH0(K3))).le.3d0)then
       chi2gam=agg*(SSIG(1)-mugcengg)**2+cgg*(SSIG(2)-muvcengg)**2
     .     +2d0*bgg*(SSIG(1)-mugcengg)*(SSIG(2)-muvcengg)
       chi2zz=azz*(SSIG(5)-mugcenzz)**2+czz*(SSIG(6)-muvcenzz)**2
     .    +2d0*bzz*(SSIG(5)-mugcenzz)*(SSIG(6)-muvcenzz)
        endif

       chi2bb=abb*(SSIG(3)-mugcenbb)**2+cbb*(SSIG(4)-muvcenbb)**2
     .    +2d0*bbb*(SSIG(3)-mugcenbb)*(SSIG(4)-muvcenbb)

        if(dabs(dsqrt(MH0(K1))-dsqrt(MH0(K3))).gt.3d0)then
         if(dabs(dsqrt(MH0(K1))-dsqrt(MH0(K2))).lt.3d0.and.
     .       dabs(dsqrt(MH0(K2))-dsqrt(MH0(K3))).gt.3d0)then
         SSIG(1)=SIG(K1,8)+SIG(K2,8)
         SSIG(2)=SIG(K1,7)+SIG(K2,7)
         SSIG(3)=SIG(K1,2)+SIG(K2,2)
         SSIG(4)=SIG(K1,3)+SIG(K2,3)
         SSIG(5)=SIG(K1,6)+SIG(K2,6)
         SSIG(6)=SIG(K1,5)+SIG(K2,5)
         chi2gam=min(agg*(SSIG(1)-mugcengg)**2
     .                       +cgg*(SSIG(2)-muvcengg)**2
     .     +2d0*bgg*(SSIG(1)-mugcengg)*(SSIG(2)-muvcengg),
     .     agg*(SIG(K3,1)-mugcengg)**2+cgg*(SIG(K3,2)-muvcengg)**2
     .     +2d0*bgg*(SIG(K3,1)-mugcengg)*(SIG(K3,2)-muvcengg))
         chi2zz=min(azz*(SSIG(5)-mugcenzz)**2
     .                       +czz*(SSIG(6)-muvcenzz)**2
     .    +2d0*bzz*(SSIG(5)-mugcenzz)*(SSIG(6)-muvcenzz),
     .      azz*(SIG(K3,5)-mugcenzz)**2+czz*(SIG(K3,6)-muvcenzz)**2
     .    +2d0*bzz*(SIG(K3,5)-mugcenzz)*(SIG(K3,6)-muvcenzz))
         endif
         if(dabs(dsqrt(MH0(K1))-dsqrt(MH0(K2))).gt.3d0.and.
     .       dabs(dsqrt(MH0(K2))-dsqrt(MH0(K3))).lt.3d0)then
         SSIG(1)=SIG(K3,8)+SIG(K2,8)
         SSIG(2)=SIG(K3,7)+SIG(K2,7)
         SSIG(3)=SIG(K3,2)+SIG(K2,2)
         SSIG(4)=SIG(K3,3)+SIG(K2,3)
         SSIG(5)=SIG(K3,6)+SIG(K2,6)
         SSIG(6)=SIG(K3,5)+SIG(K2,5)
         chi2gam=min(agg*(SSIG(1)-mugcengg)**2
     .                              +cgg*(SSIG(2)-muvcengg)**2
     .     +2d0*bgg*(SSIG(1)-mugcengg)*(SSIG(2)-muvcengg),
     .     agg*(SIG(K1,1)-mugcengg)**2+cgg*(SIG(K1,2)-muvcengg)**2
     .     +2d0*bgg*(SIG(K1,1)-mugcengg)*(SIG(K1,2)-muvcengg))
         chi2zz=min(azz*(SSIG(5)-mugcenzz)**2
     .                              +czz*(SSIG(6)-muvcenzz)**2
     .    +2d0*bzz*(SSIG(5)-mugcenzz)*(SSIG(6)-muvcenzz),
     .      azz*(SIG(K1,5)-mugcenzz)**2+czz*(SIG(K1,6)-muvcenzz)**2
     .    +2d0*bzz*(SIG(K1,5)-mugcenzz)*(SIG(K1,6)-muvcenzz))
         endif
         if(dabs(dsqrt(MH0(K1))-dsqrt(MH0(K2))).lt.3d0.and.
     .       dabs(dsqrt(MH0(K2))-dsqrt(MH0(K3))).lt.3d0)then
         SSIG(1)=SIG(K1,8)+SIG(K2,8)
         SSIG(2)=SIG(K1,7)+SIG(K2,7)
         SSIG(3)=SIG(K1,2)+SIG(K2,2)
         SSIG(4)=SIG(K1,3)+SIG(K2,3)
         SSIG(5)=SIG(K1,6)+SIG(K2,6)
         SSIG(6)=SIG(K1,5)+SIG(K2,5)
         chi2gam=agg*(SSIG(1)-mugcengg)**2+cgg*(SSIG(2)-muvcengg)**2
     .     +2d0*bgg*(SSIG(1)-mugcengg)*(SSIG(2)-muvcengg)
         chi2zz=azz*(SSIG(5)-mugcenzz)**2+czz*(SSIG(6)-muvcenzz)**2
     .    +2d0*bzz*(SSIG(5)-mugcenzz)*(SSIG(6)-muvcenzz)
         SSIG(1)=SIG(K3,8)+SIG(K2,8)
         SSIG(2)=SIG(K3,7)+SIG(K2,7)
         SSIG(3)=SIG(K3,2)+SIG(K2,2)
         SSIG(4)=SIG(K3,3)+SIG(K2,3)
         SSIG(5)=SIG(K3,6)+SIG(K2,6)
         SSIG(6)=SIG(K3,5)+SIG(K2,5)
         chi2gam=min(agg*(SSIG(1)-mugcengg)**2
     .                  +cgg*(SSIG(2)-muvcengg)**2
     .     +2d0*bgg*(SSIG(1)-mugcengg)*(SSIG(2)-muvcengg),chi2gam)
         chi2zz=min(azz*(SSIG(5)-mugcenzz)**2
     .                  +czz*(SSIG(6)-muvcenzz)**2
     .    +2d0*bzz*(SSIG(5)-mugcenzz)*(SSIG(6)-muvcenzz),chi2zz)
         endif
        endif

      ENDIF

      IF(masstest.gt.1d0)PROB(46)=masstest

      PROB(47)=DDIM(chi2gam/chi2MAX,1d0)
      PROB(48)=DDIM(chi2bb/chi2MAX,1d0)
      PROB(49)=DDIM(chi2zz/chi2MAX,1d0)


c        V- CMS constraints on h[125]->2A->4mu

      PROB(52)=0d0

      DO I1=2,5

      CJ2=0d0
      BRHTOAA=0d0

      IF(dabs(dsqrt(MH0(I1))-MHcen).le.D1)THEN
       CJ2=CJ(I1)**2+CJP(I1)**2
       BRHTOAA=BRHHH(I1,1)
      ENDIF

      I=1
      DOWHILE(dsqrt(MH0(1)).GE.X1(1).AND.X1(I).LE.dsqrt(MH0(1))
     .                                 .AND.I.LT.NY)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND.dsqrt(MH0(1)).LT.X1(NY))THEN
       LIMIT=LSIGBR(I-1)
     . +(dsqrt(MH0(1))-X1(I-1))/(X1(I)-X1(I-1))*(LSIGBR(I)-LSIGBR(I-1))
       PROB(52)=PROB(52)+DDIM(CJ2*BRHTOAA*BRMM(1)**2/LIMIT,1d0)
      ENDIF

      ENDDO

      RETURN

*   Error catch

 1    WRITE(*,*)"Cannot find the file ",FILENAME
      STOP

 2    WRITE(*,*)"Read error in the file ",FILENAME
      STOP

      END


      SUBROUTINE LHC_HSMAA_LEPTONS_CPV(PROB)

*     PROB(52) =/= 0: excluded 

      IMPLICIT NONE

      CHARACTER*256 FILENAME,EXPCON_PATH,catpath

      INTEGER NHAATAUS1,NHAATAUS2,NHAABMU
      INTEGER NHAAMUS1,NHAAMUS2,NHAAMUS3
      INTEGER NHAAMUS4,I,J,K,NM
      PARAMETER(NM=14)

      DOUBLE PRECISION MHC,XC(2,2),MH0(5),XH(5,5),MA2
      DOUBLE PRECISION PROB(*)
      DOUBLE PRECISION D1,D2,CJ2BRHTOAA,MH,MHcen,LIMIT
      DOUBLE PRECISION MHSM(NM),SMXS_8TeV(NM),SMXS_125_8TeV,SMXS_J_8TeV
      DOUBLE PRECISION LOWBOUND,HIGHBOUND

      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5),
     . BRCC(5),BRBB(5),BRTT(5),BRWW(5),BRZZ(5),BRGG(5),BRZG(5)
      DOUBLE PRECISION BRHHH(5,10),BRHCHC(5),BRHAZ(5,4),BRHCW(5),
     . BRHIGGS(5)
      DOUBLE PRECISION CU(5),CUP(5),CD(5),CDP(5),CB(5),CBP(5),CJ(5),
     . CJP(5),CI(5),CG(5),CGP(5),CV(5),CZG(5),CZGP(5)
      DOUBLE PRECISION MHmin,MHmax,chi2max,chi2gam,chi2bb,chi2zz
      DOUBLE PRECISION HAATAUS1(92,2),HAATAUS2(92,2),HAABMU(76,2)
      DOUBLE PRECISION HAAMUS1(20,2),HAAMUS2(12,2),HAAMUS3(7,2)
      DOUBLE PRECISION HAAMUS4(67,4)

      COMMON/HISPEC/MHC,XC,MH0,XH,MA2
      COMMON/HNSMBR/BRJJ,BREE,BRMM,BRLL,BRCC,BRBB,BRTT,BRWW,BRZZ,
     . BRGG,BRZG
      COMMON/HNHIBR/BRHHH,BRHCHC,BRHAZ,BRHCW,BRHIGGS
      COMMON/HNSMCOUP/CU,CUP,CD,CDP,CB,CBP,CJ,CJP,CI,CG,CGP,CV,CZG,CZGP
      COMMON/HIGGSFIT/MHmin,MHmax,chi2max,chi2gam,chi2bb,chi2zz
      COMMON/LHCHAA/HAATAUS1,HAATAUS2,HAABMU,
     .      HAAMUS1,HAAMUS2,HAAMUS3,HAAMUS4,
     .      NHAATAUS1,NHAATAUS2,NHAABMU,
     .      NHAAMUS1,NHAAMUS2,NHAAMUS3,NHAAMUS4

      DATA MHSM/85d0,90d0,95d0,100d0,105d0,110d0,115d0,120d0,125d0,
     . 130d0,135d0,140d0,145d0,150d0/
      DATA SMXS_8TeV/3.940d1,3.526d1,3.175d1,2.873d1,2.611d1,2.383d1,
     . 2.184d1,2.008d1,1.851d1,1.712d1,1.587d1,1.475d1,1.375d1,1.284d1/           ! in pb

* Determining the SM-like Higgs and its couplings / BR / XS

      CJ2BRHTOAA=0d0
      MHcen=(MHmin+MHmax)/2d0
      MH=0d0
      D2=0d0

      DO I=2,5
       D1=DDIM(dsqrt(MH0(I))/MHMAX,1d0)-DDIM(1d0,dsqrt(MH0(I))/MHMIN)
       IF(D1.EQ.0d0)THEN
        CJ2BRHTOAA=CJ2BRHTOAA+(CJ(I)**2+CJP(I)**2)*BRHHH(I,1)
        MH=MH+(CJ(I)**2+CJP(I)**2)*dsqrt(MH0(I))
        D2=D2+CJ(I)**2+CJP(I)**2
       ENDIF
      ENDDO
      MH=MH/Max(D2,1d-10)
      D1=DDIM(MH/MHMAX,1d0)-DDIM(1d0,MH/MHMIN)
      IF(D1.ne.0d0)THEN
       MH=MHcen
      ENDIF

      SMXS_125_8TeV=0d0
      I=1
      DOWHILE(MH.GE.MHSM(1) .AND. MHSM(I).LE.MH .AND. I.LT.NM)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MH.LT.MHSM(NM))THEN
       SMXS_125_8TeV=SMXS_8TeV(I-1)
     .    +(MH-MHSM(I-1))/(MHSM(I)-MHSM(I-1))
     .         *(SMXS_8TeV(I)-SMXS_8TeV(I-1))
      ENDIF

* Constraints from ggF->HSM->AA->2mu2tau from 1505.01609 (ATLAS), M_A < 50GeV
* BR(A->2mu) converted in BR(A->2tau), normalized to SM XS

      I=1
      DOWHILE(dsqrt(MH0(1)).GE.HAATAUS1(1,1)
     .        .AND. HAATAUS1(I,1).LE.dsqrt(MH0(1))
     .        .AND. I.LT.NHAATAUS1)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. dsqrt(MH0(1)).LT.HAATAUS1(NHAATAUS1,1))THEN
       LIMIT=HAATAUS1(I-1,2)
     .    +(dsqrt(MH0(1))-HAATAUS1(I-1,1))
     .         /(HAATAUS1(I,1)-HAATAUS1(I-1,1))
     .         *(HAATAUS1(I,2)-HAATAUS1(I-1,2))
       PROB(52)=PROB(52)+DDIM(CJ2BRHTOAA*BRLL(1)**2/LIMIT,1d0)
      ENDIF


* Constraints from ggF->HSM->AA->4tau from 1510.06534 (CMS), 4GeV < M_A < 8GeV

      I=1
      DOWHILE(dsqrt(MH0(1)).GE.HAATAUS2(1,1)
     .        .AND. HAATAUS2(I,1).LE.dsqrt(MH0(1))
     .        .AND. I.LT.NHAATAUS2)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. dsqrt(MH0(1)).LT.HAATAUS2(NHAATAUS2,1))THEN
       LIMIT=HAATAUS2(I-1,2)
     .    +(dsqrt(MH0(1))-HAATAUS2(I-1,1))
     .         /(HAATAUS2(I,1)-HAATAUS2(I-1,1))
     .         *(HAATAUS2(I,2)-HAATAUS2(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(SMXS_125_8TeV*CJ2BRHTOAA*BRLL(1)**2/LIMIT,1d0)
      ENDIF


* Constraints from ggF->HSM->AA->2mu2b from CMS-PAS-HIG-14-041, 25GeV < M_A < 63GeV
* normalized to SM XS

      I=1
      DOWHILE(dsqrt(MH0(1)).GE.HAABMU(1,1)
     .        .AND. HAABMU(I,1).LE.dsqrt(MH0(1))
     .        .AND. I.LT.NHAABMU)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. dsqrt(MH0(1)).LT.HAABMU(NHAABMU,1))THEN
       LIMIT=HAABMU(I-1,2)
     .    +(dsqrt(MH0(1))-HAABMU(I-1,1))
     .         /(HAABMU(I,1)-HAABMU(I-1,1))
     .         *(HAABMU(I,2)-HAABMU(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(CJ2BRHTOAA*BRMM(1)*BRBB(1)/LIMIT,1d0)
      ENDIF


* Constraints from ggF->HSM->AA->2mu2tau from CMS-PAS-HIG-15-011, 19GeV < M_A < 57GeV
* BR(A->2tau) converted in BR(A->2mu), normalized to SM XS

      I=1
      DOWHILE(dsqrt(MH0(1)).GE.HAAMUS1(1,1)
     .        .AND. HAAMUS1(I,1).LE.dsqrt(MH0(1))
     .        .AND. I.LT.NHAAMUS1)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. dsqrt(MH0(1)).LT.HAAMUS1(NHAAMUS1,1))THEN
       LIMIT=HAAMUS1(I-1,2)
     .    +(dsqrt(MH0(1))-HAAMUS1(I-1,1))
     .         /(HAAMUS1(I,1)-HAAMUS1(I-1,1))
     .         *(HAAMUS1(I,2)-HAAMUS1(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(CJ2BRHTOAA*BRMM(1)**2/LIMIT,1d0)
      ENDIF


* Constraints from ggF->HSM->AA->4tau from CMS-PAS-HIG-14-022, 5GeV < M_A < 14GeV
* BR(A->2tau) converted in BR(A->2mu), as shown in Fig7 of CMS-PAS-HIG-15-011,
* normalized to SM XS

      I=1
      DOWHILE(dsqrt(MH0(1)).GE.HAAMUS2(1,1)
     .        .AND. HAAMUS2(I,1).LE.dsqrt(MH0(1))
     .        .AND. I.LT.NHAAMUS2)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. dsqrt(MH0(1)).LT.HAAMUS2(NHAAMUS2,1))THEN
       LIMIT=HAAMUS2(I-1,2)
     .    +(dsqrt(MH0(1))-HAAMUS2(I-1,1))
     .         /(HAAMUS2(I,1)-HAAMUS2(I-1,1))
     .         *(HAAMUS2(I,2)-HAAMUS2(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(CJ2BRHTOAA*BRMM(1)**2/LIMIT,1d0)
      ENDIF


* Constraints from ggF->HSM->AA->4tau from CMS-PAS-HIG-14-019, 3GeV < M_A < 8GeV
* BR(A->2tau) converted in BR(A->2mu), as shown in Fig7 of CMS-PAS-HIG-15-011,
* normalized to SM XS

      I=1
      DOWHILE(dsqrt(MH0(1)).GE.HAAMUS3(1,1)
     .        .AND. HAAMUS3(I,1).LE.dsqrt(MH0(1))
     .        .AND. I.LT.NHAAMUS3)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. dsqrt(MH0(1)).LT.HAAMUS3(NHAAMUS3,1))THEN
       LIMIT=HAAMUS3(I-1,2)
     .    +(dsqrt(MH0(1))-HAAMUS3(I-1,1))
     .         /(HAAMUS3(I,1)-HAAMUS3(I-1,1))
     .         *(HAAMUS3(I,2)-HAAMUS3(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(CJ2BRHTOAA*BRMM(1)**2/LIMIT,1d0)
      ENDIF


* Constraints from ggF->H->AA->4mu from 1506.00424, 0.25GeV < M_A < 3.55GeV, 85GeV < m_H < 150GeV

      CJ2BRHTOAA=0d0

        DO J=2,5

      IF(J.ge.3 .and. dabs(dsqrt(MH0(J))-dsqrt(MH0(J-1))).lt.3d0)THEN
       IF(J.ge.4 .and. dabs(dsqrt(MH0(J))-dsqrt(MH0(J-2))).lt.3d0)THEN
        MH=((CJ(J-2)**2+CJP(J-2)**2)*dsqrt(MH0(J-2))
     .    +(CJ(J-1)**2+CJP(J-1)**2)*dsqrt(MH0(J-1))
     .    +(CJ(J)**2+CJP(J)**2)*dsqrt(MH0(J)))
     .  /Max(CJ(J)**2+CJ(J-1)**2+CJ(J-2)**2+CJP(J)**2+CJP(J-1)**2
     .                                           +CJP(J-2)**2,1d-10)
        CJ2BRHTOAA=CJ2BRHTOAA+(CJ(J)**2+CJP(J)**2)*BRHHH(J,1)
       ELSE
        MH=((CJ(J-1)**2+CJP(J-1)**2)*dsqrt(MH0(J-1))
     .    +(CJ(J)**2+CJP(J)**2)*dsqrt(MH0(J)))
     .            /Max(CJ(J)**2+CJ(J-1)**2+CJP(J)**2+CJP(J-1)**2,1d-10)
        CJ2BRHTOAA=CJ2BRHTOAA+(CJ(J)**2+CJP(J)**2)*BRHHH(J,1)
       ENDIF
      ELSE
       MH=dsqrt(MH0(J))
       CJ2BRHTOAA=(CJ(J)**2+CJP(J)**2)*BRHHH(J,1)
      ENDIF

      I=1
      DOWHILE(MH.GE.MHSM(1)
     .        .AND.MHSM(I).LE.MH.AND.I.LT.NM)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MH.LT.MHSM(14))THEN
       SMXS_J_8TeV=SMXS_8TeV(I-1)
     .    +(MH-MHSM(I-1))/(MHSM(I)-MHSM(I-1))
     .         *(SMXS_8TeV(I)-SMXS_8TeV(I-1))
       SMXS_J_8TeV=SMXS_J_8TeV*1d3         ! converting in fb
      ENDIF


      I=1
      DOWHILE(MH.GE.HAAMUS4(1,1)
     .        .AND. HAAMUS4(I,1).LE.MH .AND. I.LT.NHAAMUS4)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MH.LT.HAAMUS4(NHAAMUS4,1))THEN
       LOWBOUND=100d0
       HIGHBOUND=100d0
       IF(dsqrt(MH0(1)).ge.0.25d0 .and. dsqrt(MH0(1)).lt.2d0)THEN
        LOWBOUND=HAAMUS4(I-1,2)+(dsqrt(MH0(1))-0.25d0)/(2d0-0.25d0)
     .         *(HAAMUS4(I-1,3)-HAAMUS4(I-1,2))
        HIGHBOUND=HAAMUS4(I,2)+(dsqrt(MH0(1))-0.25d0)/(2d0-0.25d0)
     .         *(HAAMUS4(I,3)-HAAMUS4(I,2))
       ELSEIF(dsqrt(MH0(1)).ge.2d0 .and. dsqrt(MH0(1)).lt.3.55d0)THEN
        LOWBOUND=HAAMUS4(I-1,3)+(dsqrt(MH0(1))-2d0)/(3.55d0-2d0)
     .         *(HAAMUS4(I-1,4)-HAAMUS4(I-1,3))
        HIGHBOUND=HAAMUS4(I,3)+(dsqrt(MH0(1))-2d0)/(3.55d0-2d0)
     .         *(HAAMUS4(I,4)-HAAMUS4(I,3))
       ENDIF

       LIMIT=LOWBOUND
     .    +(MH-HAAMUS4(I-1,1))/(HAAMUS4(I,1)-HAAMUS4(I-1,1))
     .         *(HIGHBOUND-LOWBOUND)
       PROB(52)=PROB(52)
     .          +DDIM(SMXS_J_8TeV*CJ2BRHTOAA*BRMM(1)**2/LIMIT,1d0)
      ENDIF

        ENDDO

      RETURN

      END
