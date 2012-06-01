      SUBROUTINE FILLGAMBR2(NSMIN,NSSIN,SMPARA_H,SSPARA_H,NFLAG,IFLAG_H
     . ,MCH,HMASS_H,NCMAX,NHC_H,SHC_H,CHC_H,STMASS_H
     . ,STMIX_H,SBMASS_H,SBMIX_H,STAUMASS_H,STAUMIX_H,SNU3MASS_H
     . ,MC_H,MN_H,NMNH,GAMBRN,NMCH,GAMBRC)
************************************************************************
*
*JSL 10/Jun/2009: Input arrays are extended to include SMPARA_H(NSMIN)
*                 and SSPARA_H(NSSIN)
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*-----------------------------------------------------------------------
*+CDE HC_ COMMON BLOCKS:
      COMMON /HC_SMPARA/ AEM_H,ASMZ_H,MZ_H,SW_H,ME_H,MMU_H,MTAU_H,MDMT_H
     .                  ,MSMT_H,MBMT_H,MUMT_H,MCMT_H,MTPOLE_H,GAMW_H
     .                  ,GAMZ_H,EEM_H,ASMT_H,CW_H,TW_H,MW_H,GW_H,GP_H
     .                  ,V_H,GF_H,MTMT_H
*
      COMMON /HC_RSUSYPARA/ TB_H,CB_H,SB_H,MQ3_H,MU3_H,MD3_H,ML3_H,ME3_H
*
      COMPLEX*16 MU_H,M1_H,M2_H,M3_H,AT_H,AB_H,ATAU_H
      COMMON /HC_CSUSYPARA/ MU_H,M1_H,M2_H,M3_H,AT_H,AB_H,ATAU_H
*
*NEW COMMON BLOCKS for V2
*
      REAL*8     RAUX_H(999)
      COMPLEX*16 CAUX_H(999)
      COMMON /HC_RAUX/ RAUX_H
      COMMON /HC_CAUX/ CAUX_H
      DATA NAUX/999/
*-----------------------------------------------------------------------
* Input Array:
      REAL*8     SMPARA_H(NSMIN),SSPARA_H(NSSIN)
      INTEGER*8  IFLAG_H(NFLAG)
      REAL*8     HMASS_H(3)
      COMPLEX*16 NHC_H(NCMAX,3)
      REAL*8     SHC_H(NCMAX)
      COMPLEX*16 CHC_H(NCMAX)
      REAL*8     STMASS_H(2),SBMASS_H(2),STAUMASS_H(2),SNU3MASS_H
      COMPLEX*16 STMIX_H(2,2),SBMIX_H(2,2),STAUMIX_H(2,2)
      REAL*8     MC_H(2),MN_H(4)
*-----------------------------------------------------------------------
* Output Array:
      REAL*8 GAMBRN(NMNH,3,3)   
      REAL*8 GAMBRC(NMCH,3)      
*-----------------------------------------------------------------------
* For integration
      COMMON /HC_BODE/ EPSV,OMEGAI,OMEGAJ,XUP,XDW
      EXTERNAL FVVS,FHVS
*-----------------------------------------------------------------------
* Local :
      INTEGER*8 ISMN,ISUSYN,ISMC,ISUSYC
      COMPLEX*16 GHV,GSS
      COMPLEX*16 XI
      COMPLEX*16 GF_UD,GS_UD,GP_UD
      COMPLEX*16 GF_CS,GS_CS,GP_CS
      COMPLEX*16 GF_TB,GS_TB,GP_TB
      COMPLEX*16 GF1,GS1,GP1
      COMPLEX*16 GF2,GS2,GP2
      COMPLEX*16 GF3,GS3,GP3
      COMPLEX*16 GF4,GS4,GP4
* Radiative corrections to Htt Hbb Yukawa couplings
      COMPLEX*16 HB_H,HT_H,CKB_H,CKT_H
      COMPLEX*16 RB_H,RT_H,CKBB_H,CKBT_H
*JSL 10/Jun/2009: Including threshold corrections
      COMPLEX*16 CKD_H,CKS_H
* Higgs-photon-photon
      COMPLEX*16 SPH,PPH
      COMPLEX*16 SPP(15),PPP(7)
*-----------------------------------------------------------------------
      ISMN   = IFLAG_H(20)
      ISUSYN = IFLAG_H(21)
      ISMC   = IFLAG_H(22)
      ISUSYC = IFLAG_H(23)
*      print*,ismn,isusyn,ismc,isusyc
      PI=2.D0*DASIN(1.D0)
*-----------------------------------------------------------------------
*
* << NEUTRAL HIGGS BOSON DECAYS INTO SM PARTICLES >>

      DO IH=1,3
      MH    = HMASS_H(IH)
      SQRTS = MH
*=======================================================================
*---> running alpha_s and b-quark mass at SQRTS=MH scale
*     : SQRTS > MS^pole assumed
*-----------------------------------------------------------------------
      PI      = 2.D0*DASIN(1.D0)
      B3      = (11.D0-2.D0/3.D0*3.D0)/4.D0/PI
      B4      = (11.D0-2.D0/3.D0*4.D0)/4.D0/PI
      B5      = (11.D0-2.D0/3.D0*5.D0)/4.D0/PI
      B6      = (11.D0-2.D0/3.D0*6.D0)/4.D0/PI
*
      AS_MT   = ASMT_H
      AS_MZ   = ASMZ_H
      AS_MB   = RAUX_H(3)
      AS_MC   = RAUX_H(6)
      MB_POLE = RAUX_H(1)
      MC_POLE = RAUX_H(4)
*Quark masses at mb^pole and mc^pole
      MT_MB = MTMT_H*(AS_MB/AS_MT)**(1.D0/B5/PI)
      MB_MB = MBMT_H*(AS_MB/AS_MT)**(1.D0/B5/PI)
      MC_MB = MCMT_H*(AS_MB/AS_MT)**(1.D0/B5/PI)
      MS_MB = MSMT_H*(AS_MB/AS_MT)**(1.D0/B5/PI)
      MU_MB = MUMT_H*(AS_MB/AS_MT)**(1.D0/B5/PI)
      MD_MB = MDMT_H*(AS_MB/AS_MT)**(1.D0/B5/PI)
*      print*,'at mb^pole:',mt_mb,mb_mb,mc_mb,ms_mb,mu_mb,md_mb
*      print*,'at mb^pole:',mb_mb,' ?= ',raux_h(2)
      MT_MC = MT_MB *(AS_MC/AS_MB)**(1.D0/B4/PI)
      MB_MC = MB_MB *(AS_MC/AS_MB)**(1.D0/B4/PI)
      MC_MC = MC_MB *(AS_MC/AS_MB)**(1.D0/B4/PI)
      MS_MC = MS_MB *(AS_MC/AS_MB)**(1.D0/B4/PI)
      MU_MC = MU_MB *(AS_MC/AS_MB)**(1.D0/B4/PI)
      MD_MC = MD_MB *(AS_MC/AS_MB)**(1.D0/B4/PI)
*      print*,'at mc^pole:',mt_mc,mb_mc,mc_mc,ms_mc,mu_mc,md_mc
*      print*,'at mc^pole:',mc_mc,' ?= ',raux_h(5)
*-----
*AS(SQRTS)
*  mt^pole < ss
      IF(SQRTS.GT.MTPOLE_H) THEN
       AS_S = AS_MT/(1.D0+B6*AS_MT*DLOG(SQRTS**2/MTPOLE_H**2))
*  mb^pole < ss <=mt^pole
      ELSEIF(SQRTS.LE.MTPOLE_H .AND. SQRTS.GT.MB_POLE ) THEN
       AS_S = AS_MZ/(1.D0+B5*AS_MZ*DLOG(SQRTS**2/MZ_H**2))
*  mc^pole < ss <=mb^pole
      ELSEIF(SQRTS.LE.MB_POLE .AND. SQRTS.GT.MC_POLE  ) THEN
       AS_S = AS_MB/(1.D0+B4*AS_MB*DLOG(SQRTS**2/MB_POLE**2))
*            ss <=mc^pole
      ELSEIF(SQRTS.LE.MC_POLE) THEN
       AS_S = AS_MC/(1.D0+B3*AS_MC*DLOG(SQRTS**2/MC_POLE**2))
      ELSE
       print*,'SQRTS = ',sqrts,' is out of range !!!'
       STOP
      ENDIF
*MQ(SQRTS)
*  mt^pole < ss
      IF(SQRTS.GT.MTPOLE_H) THEN
       MT_S = MTMT_H*(AS_S/AS_MT)**(1.D0/B6/PI)
       MB_S = MBMT_H*(AS_S/AS_MT)**(1.D0/B6/PI)
       MC_S = MCMT_H*(AS_S/AS_MT)**(1.D0/B6/PI)
       MS_S = MSMT_H*(AS_S/AS_MT)**(1.D0/B6/PI)
       MU_S = MUMT_H*(AS_S/AS_MT)**(1.D0/B6/PI)
       MD_S = MDMT_H*(AS_S/AS_MT)**(1.D0/B6/PI)
*  mb^pole < ss <=mt^pole
      ELSEIF(SQRTS.LE.MTPOLE_H .AND. SQRTS.GT.MB_POLE ) THEN
       MT_S = MTMT_H*(AS_S/AS_MT)**(1.D0/B5/PI)
       MB_S = MBMT_H*(AS_S/AS_MT)**(1.D0/B5/PI)
       MC_S = MCMT_H*(AS_S/AS_MT)**(1.D0/B5/PI)
       MS_S = MSMT_H*(AS_S/AS_MT)**(1.D0/B5/PI)
       MU_S = MUMT_H*(AS_S/AS_MT)**(1.D0/B5/PI)
       MD_S = MDMT_H*(AS_S/AS_MT)**(1.D0/B5/PI)
*  mc^pole < ss <=mb^pole
      ELSEIF(SQRTS.LE.MB_POLE .AND. SQRTS.GT.MC_POLE  ) THEN
       MT_S = MT_MB *(AS_S/AS_MB)**(1.D0/B4/PI)
       MB_S = MB_MB *(AS_S/AS_MB)**(1.D0/B4/PI)
       MC_S = MC_MB *(AS_S/AS_MB)**(1.D0/B4/PI)
       MS_S = MS_MB *(AS_S/AS_MB)**(1.D0/B4/PI)
       MU_S = MU_MB *(AS_S/AS_MB)**(1.D0/B4/PI)
       MD_S = MD_MB *(AS_S/AS_MB)**(1.D0/B4/PI)
*            ss <=mc^pole
      ELSEIF(SQRTS.LE.MC_POLE) THEN
       MT_S = MT_MC *(AS_S/AS_MC)**(1.D0/B3/PI)
       MB_S = MB_MC *(AS_S/AS_MC)**(1.D0/B3/PI)
       MC_S = MC_MC *(AS_S/AS_MC)**(1.D0/B3/PI)
       MS_S = MS_MC *(AS_S/AS_MC)**(1.D0/B3/PI)
       MU_S = MU_MC *(AS_S/AS_MC)**(1.D0/B3/PI)
       MD_S = MD_MC *(AS_S/AS_MC)**(1.D0/B3/PI)
      ELSE
       print*,'SQRTS = ',sqrts,' is out of range !!!'
       STOP
      ENDIF
*
*      print*,'SQRTS,AS(SQRTS) =',SQRTS,AS_S
*      print*,' > MQ(SQRTS)    :',MT_S,MB_S,MC_S,MS_S
*-----------------------------------------------------------------------
      ASMH=AS_S
      MTMH=MT_S
      MBMH=MB_S
      MCMH=MC_S
      MSMH=MS_S
      MUMH=MU_S
      MDMH=MD_S
*=======================================================================
*
*---> H_IH -> e+ e-     [IM= 1]
      IM   = 1
      IFF  = 1
      MJ   = ME_H
      MK   = ME_H
      CF   = 1.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*---> H_IH -> mu+ mu-   [IM= 2]
      IM   = 2
      IFF  = 4
      MJ   = MMU_H
      MK   = MMU_H
      CF   = 1.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*---> H_IH -> tau+ tau- [IM= 3]
      IM   = 3
      IFF  = 7
      MJ   = MTAU_H
      MK   = MTAU_H
      CF   = 1.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*---> H_IH -> d d       [IM= 4]
      IM   = 4
      IFF  = 10
      MJ   = MDMH
      MK   = MDMH
      CF   = 3.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,MDMH/MDMT_H*NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*K-factor[See, for example, hep-ph/0305101, Eq.(6) and (7)]
      GAMBRN(IM,1,IH)=(1.D0+5.67D0*ASMH/PI)*GAMBRN(IM,1,IH)
*---> H_IH -> s s       [IM= 5]
      IM   = 5
      IFF  = 13
      MJ   = MSMH
      MK   = MSMH
      CF   = 3.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,MSMH/MSMT_H*NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*K-factor[See, for example, hep-ph/0305101, Eq.(6) and (7)]
      GAMBRN(IM,1,IH)=(1.D0+5.67D0*ASMH/PI)*GAMBRN(IM,1,IH)
*---> H_IH -> b b       [IM= 6]
      IM   = 6
      IFF  = 16
      MJ   = MBMH
      MK   = MBMH
      CF   = 3.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,MBMH/MBMT_H*NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*      print*,'SAME?',GAMBRN(IM,1,IH)
*recalculate taking into account running b-quark mass
*      IF((MH-2.D0*MBMH).GT.0.D0) THEN
*       BETABB=DSQRT(1.D0-4.D0*MBMH**2/MH**2)
*      ELSE
*       BETABB=0.D0
*      ENDIF
*      GBB=3.D0*(GW_H*MBMH/2.D0/MW_H)**2*MH*BETABB/8.D0/PI
*     .   *(BETABB**2*CDABS(NHC_H(IFF+1,IH))**2
*     .              +CDABS(NHC_H(IFF+2,IH))**2)
*      print*,'SAME?',GBB
*K-factor[See, for example, hep-ph/0305101, Eq.(6) and (7)]
      GAMBRN(IM,1,IH)=(1.D0+5.67D0*ASMH/PI)*GAMBRN(IM,1,IH)
*---> H_IH -> u u       [IM= 7]
      IM   = 7
      IFF  = 19
      MJ   = MUMH
      MK   = MUMH
      CF   = 3.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,MUMH/MUMT_H*NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*K-factor[See, for example, hep-ph/0305101, Eq.(6) and (7)]
      GAMBRN(IM,1,IH)=(1.D0+5.67D0*ASMH/PI)*GAMBRN(IM,1,IH)
*---> H_IH -> c c       [IM= 8]
      IM   = 8
      IFF  = 22
      MJ   = MCMH
      MK   = MCMH
      CF   = 3.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,MCMH/MCMT_H*NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*K-factor[See, for example, hep-ph/0305101, Eq.(6) and (7)]
      GAMBRN(IM,1,IH)=(1.D0+5.67D0*ASMH/PI)*GAMBRN(IM,1,IH)
*---> H_IH -> t t       [IM= 9]
      IM   = 9
      IFF  = 25
      MJ   = MTMH
      MK   = MTMH
      CF   = 3.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,MTMH/MTMT_H*NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*K-factor[See, for example, hep-ph/0305101, Eq.(6) and (7)]
      GAMBRN(IM,1,IH)=(1.D0+5.67D0*ASMH/PI)*GAMBRN(IM,1,IH)
* 
*      print*,'H',ih,' -> e   e   : ',gambrn(1,1,ih)
*      print*,'H',ih,' -> mu  mu  : ',gambrn(2,1,ih)
*      print*,'H',ih,' -> tau tau : ',gambrn(3,1,ih)
*      print*,'H',ih,' -> d   d   : ',gambrn(4,1,ih)
*      print*,'H',ih,' -> s   s   : ',gambrn(5,1,ih)
*      print*,'H',ih,' -> b   b   : ',gambrn(6,1,ih)
*      print*,'H',ih,' -> u   u   : ',gambrn(7,1,ih)
*      print*,'H',ih,' -> c   c   : ',gambrn(8,1,ih)
*      print*,'H',ih,' -> t   t   : ',gambrn(9,1,ih)
*
*---> H_IH -> W W       [IM=10]
      IM  = 10
      IVV = 70
      DV  = 2.D0
      CALL HVV(GF_H,DV,NHC_H(IVV,IH),MW_H,HMASS_H(IH),GAMBRN(IM,1,IH))
* into W + W^* 
      DVVS=2.D0
      IF( (HMASS_H(IH).GT.MW_H) .AND. 
     .    (HMASS_H(IH).LT.(2.D0*MW_H+GAMW_H)) ) THEN
       IF(DABS(HMASS_H(IH)-2.D0*MW_H).LT.GAMW_H) THEN
        DVVS=2.D0-(HMASS_H(IH)-2.D0*MW_H+GAMW_H)/2.D0/GAMW_H
       ENDIF
      EPSV=GAMW_H/MW_H
      OMEGAI=HMASS_H(IH)**2/MW_H**2
      OMEGAJ=0.D0
      XUP=(DSQRT(OMEGAI)-1.D0)**2
      XDW=0.D0
      NSTEP=500
      CALL BODE(FVVS,0.D0,1.D0,NSTEP,RES)
      GAMBRN(IM,1,IH)=GF_H*HMASS_H(IH)**3*DV*DVVS
     . *CDABS(NHC_H(IVV,IH))**2
     . /16.D0/DSQRT(2.D0)/PI**2*EPSV/OMEGAI**3*RES
      ENDIF
*
*---> H_IH -> Z Z       [IM=11]
      IM  = 11
      IVV = 70
      DV  = 1.D0
      CALL HVV(GF_H,DV,NHC_H(IVV,IH),MZ_H,HMASS_H(IH),GAMBRN(IM,1,IH))
* into Z + Z^* 
      DVVS=2.D0
      IF( (HMASS_H(IH).GT.MZ_H) .AND. 
     .    (HMASS_H(IH).LT.(2.D0*MZ_H+GAMZ_H)) ) THEN
       IF(DABS(HMASS_H(IH)-2.D0*MZ_H).LT.GAMZ_H) THEN
        DVVS=2.D0-(HMASS_H(IH)-2.D0*MZ_H+GAMZ_H)/2.D0/GAMZ_H
       ENDIF
      EPSV=GAMZ_H/MZ_H
      OMEGAI=HMASS_H(IH)**2/MZ_H**2
      OMEGAJ=0.D0
      XUP=(DSQRT(OMEGAI)-1.D0)**2
      XDW=0.D0
      NSTEP=500
      CALL BODE(FVVS,0.D0,1.D0,NSTEP,RES)
      GAMBRN(IM,1,IH)=GF_H*HMASS_H(IH)**3*DV*DVVS
     . *CDABS(NHC_H(IVV,IH))**2
     . /16.D0/DSQRT(2.D0)/PI**2*EPSV/OMEGAI**3*RES
      ENDIF
*      print*,'H',ih,' -> W   W   : ',gambrn(10,1,ih)
*      print*,'H',ih,' -> Z   Z   : ',gambrn(11,1,ih)
*---> H_IH -> H1 Z      [IM=12]
      IM  = 12
      IF(IH.EQ.1) GHV=DCMPLX(0.D0,0.D0) ! Neglecting overall sign
      IF(IH.EQ.2) GHV=NHC_H(70,3)       
      IF(IH.EQ.3) GHV=NHC_H(70,2)
      CALL HHV(GF_H,GHV,HMASS_H(IH),HMASS_H(1),MZ_H,GAMBRN(IM,1,IH))
* H_IH -> H1 Z*
      IF( (HMASS_H(IH).GT.HMASS_H(1)) .AND. 
     .    (HMASS_H(IH).LT.(HMASS_H(1)+MZ_H+5.D0*GAMZ_H)) ) THEN
      EPSV=GAMZ_H/MZ_H
      OMEGAI=HMASS_H(IH)**2/MZ_H**2
      OMEGAJ=HMASS_H(1)**2/MZ_H**2
      XUP=(DSQRT(OMEGAI)-DSQRT(OMEGAJ))**2
      XDW=0.D0
      NSTEP=500
      CALL BODE(FHVS,0.D0,1.D0,NSTEP,RES)
      GAMBRN(IM,1,IH)=GF_H*HMASS_H(IH)**3*CDABS(GHV)**2
     . /8.D0/SQRT(2.D0)/PI**2*EPSV/OMEGAI**3*RES
      ENDIF
*---> H_IH -> H2 Z      [IM=13]
      IM  = 13
      IF(IH.EQ.1) GHV=NHC_H(70,3)       ! Neglecting overall sign
      IF(IH.EQ.2) GHV=DCMPLX(0.D0,0.D0) 
      IF(IH.EQ.3) GHV=NHC_H(70,1)       
      CALL HHV(GF_H,GHV,HMASS_H(IH),HMASS_H(2),MZ_H,GAMBRN(IM,1,IH))
* H_IH -> H2 Z*
      IF( (HMASS_H(IH).GT.HMASS_H(2)) .AND. 
     .    (HMASS_H(IH).LT.(HMASS_H(2)+MZ_H+5.D0*GAMZ_H)) ) THEN
      EPSV=GAMZ_H/MZ_H
      OMEGAI=HMASS_H(IH)**2/MZ_H**2
      OMEGAJ=HMASS_H(2)**2/MZ_H**2
      XUP=(DSQRT(OMEGAI)-DSQRT(OMEGAJ))**2
      XDW=0.D0
      NSTEP=500
      CALL BODE(FHVS,0.D0,1.D0,NSTEP,RES)
      GAMBRN(IM,1,IH)=GF_H*HMASS_H(IH)**3*CDABS(GHV)**2
     . /8.D0/SQRT(2.D0)/PI**2*EPSV/OMEGAI**3*RES
      ENDIF
*      print*,'H',ih,' -> H1  Z   : ',gambrn(12,1,ih)
*      print*,'H',ih,' -> H2  Z   : ',gambrn(13,1,ih)
*---> H_IH -> H1 H1     [IM=14]
      IM   = 14
      SYMF = 2.D0
      IF(IH.EQ.1) GSS=DCMPLX(SHC_H(10),0.D0)
      IF(IH.EQ.2) GSS=DCMPLX(SHC_H(9),0.D0)
      IF(IH.EQ.3) GSS=DCMPLX(SHC_H(6),0.D0)
      CALL HSS(SYMF,V_H,GSS,HMASS_H(IH),HMASS_H(1),HMASS_H(1)
     .,GAMBRN(IM,1,IH))
*---> H_IH -> H1 H2     [IM=15]
      IM   = 15
      SYMF = 1.D0
      IF(IH.EQ.1) GSS=DCMPLX(SHC_H(9),0.D0)
      IF(IH.EQ.2) GSS=DCMPLX(SHC_H(8),0.D0)
      IF(IH.EQ.3) GSS=DCMPLX(SHC_H(5),0.D0)
      CALL HSS(SYMF,V_H,GSS,HMASS_H(IH),HMASS_H(1),HMASS_H(2)
     .,GAMBRN(IM,1,IH))
*---> H_IH -> H2 H2     [IM=16]
      IM   = 16
      SYMF = 2.D0
      IF(IH.EQ.1) GSS=DCMPLX(SHC_H(8),0.D0)
      IF(IH.EQ.2) GSS=DCMPLX(SHC_H(7),0.D0)
      IF(IH.EQ.3) GSS=DCMPLX(SHC_H(4),0.D0)
      CALL HSS(SYMF,V_H,GSS,HMASS_H(IH),HMASS_H(2),HMASS_H(2)
     .,GAMBRN(IM,1,IH))
*      print*,'H',ih,' -> H1  H1  : ',gambrn(14,1,ih)
*      print*,'H',ih,' -> H1  H2  : ',gambrn(15,1,ih)
*      print*,'H',ih,' -> H2  H2  : ',gambrn(16,1,ih)
*---> H_IH -> P P       [IM=17]
      CALL HPP(IH,MC_H,MCH,HMASS_H,STMASS_H,SBMASS_H,STAUMASS_H
     .        ,NCMAX,NHC_H,SPH,PPH,SPP,PPP)
*      print*,'---> FILLGAMBR'
*      print*,'H-P-P with IH = ',ih
*      print*,'bottom,top     ',spp(1),spp(2)
*      print*,'charm,tau      ',spp(3),spp(4)
*      print*,'c.ino11,c.ino22',spp(5),spp(6)
*      print*,'stop11,stop22  ',spp(7),spp(8)
*      print*,'sbot11,sbot22  ',spp(9),spp(10)
*      print*,'ww,c.hc.h      ',spp(11),spp(12)
*      print*,'stau11,stau22  ',spp(13),spp(14)
*      print*,'S(H-p-p)       ',nhc_h(88,ih)
*      print*,'zero?          ',nhc_h(88,ih)-spp(15)
*      print*,'zero?          ',nhc_h(88,ih)-sph
*      print*,'bottom,top     ',ppp(1),ppp(2)
*      print*,'charm,tau      ',ppp(3),ppp(4)
*      print*,'c.ino11,c.ino22',ppp(5),ppp(6)
*      print*,'P(H-p-p)       ',nhc_h(89,ih)
*      print*,'zero?',nhc_h(89,ih)-ppp(7)
*      print*,'zero?',nhc_h(89,ih)-pph
*
      DJT = -ASMH/PI
      DJSQ= 8.D0*ASMH/3.D0/PI
*      print*,DJT,DJSQ
*      print*,SPH
      SPH=SPH+DJT*SPP(2)+DJSQ*(SPP(7)+SPP(8)+SPP(9)+SPP(10))
*      print*,SPH
*
      GAMBRN(17,1,IH)=HMASS_H(IH)**3*AEM_H**2/256.D0/PI**3/V_H**2
     . *(CDABS(SPH)**2+CDABS(PPH)**2)
* ---> H_IH -> G G       [IM=18]
* running alpha_s(m_h) effect and K factor included
      IF(                        HMASS_H(IH).GT.MTMH) XNF=6.D0
      IF(HMASS_H(IH).LE.MTMH.AND.HMASS_H(IH).GT.MBMH) XNF=5.D0
      IF(HMASS_H(IH).LE.MBMH.AND.HMASS_H(IH).GT.MCMH) XNF=4.D0
      IF(HMASS_H(IH).LE.MCMH.AND.HMASS_H(IH).GT.MSMH) XNF=3.D0
      IF(HMASS_H(IH).LE.MSMH.AND.HMASS_H(IH).GT.MDMH) XNF=2.D0
      IF(HMASS_H(IH).LE.MDMH                        ) XNF=1.D0
      CKHG=1.D0+ASMH/PI*(95.D0/4.D0-7.D0/6.D0*XNF)
      CKAG=1.D0+ASMH/PI*(97.D0/4.D0-7.D0/6.D0*XNF)
      GAMBRN(18,1,IH)=HMASS_H(IH)**3*ASMH**2/32.D0/PI**3/V_H**2
     . *(CKHG*CDABS(NHC_H(84,IH))**2+CKAG*CDABS(NHC_H(85,IH))**2)
*      print*,IH,MH,ASMH,MBMH
*      print*,'H',ih,' -> P  P  : ',gambrn(17,1,ih)
*      print*,'H',ih,' -> G  G  : ',gambrn(18,1,ih)
*-----------------------------------------------------------------------
*
* << NEUTRAL HIGGS BOSON DECAYS INTO SUSY PARTICLES >>
 
*---> H_IH -> N1 N1       [IM=ISMN+1]
      IM   = ISMN+1
      IFF  = 28
      MJ   = MN_H(1)
      MK   = MN_H(1)
      CF   = 1.D0  ! Color Factor
      SYMF = 2.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 1.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*---> H_IH -> N1 N2       [IM=ISMN+2]
      IM   = ISMN+2
      IFF  = 40
      MJ   = MN_H(1)
      MK   = MN_H(2)
      CF   = 1.D0  ! Color Factor
      SYMF = 2.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*---> H_IH -> N1 N3       [IM=ISMN+3]
      IM   = ISMN+3
      IFF  = 43
      MJ   = MN_H(1)
      MK   = MN_H(3)
      CF   = 1.D0  ! Color Factor
      SYMF = 2.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*---> H_IH -> N1 N4       [IM=ISMN+4]
      IM   = ISMN+4
      IFF  = 46
      MJ   = MN_H(1)
      MK   = MN_H(4)
      CF   = 1.D0  ! Color Factor
      SYMF = 2.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*---> H_IH -> N2 N2       [IM=ISMN+5]
      IM   = ISMN+5
      IFF  = 31
      MJ   = MN_H(2)
      MK   = MN_H(2)
      CF   = 1.D0  ! Color Factor
      SYMF = 2.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 1.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*---> H_IH -> N2 N3       [IM=ISMN+6]
      IM   = ISMN+6
      IFF  = 49
      MJ   = MN_H(2)
      MK   = MN_H(3)
      CF   = 1.D0  ! Color Factor
      SYMF = 2.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*---> H_IH -> N2 N4       [IM=ISMN+7]
      IM   = ISMN+7
      IFF  = 52
      MJ   = MN_H(2)
      MK   = MN_H(4)
      CF   = 1.D0  ! Color Factor
      SYMF = 2.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*---> H_IH -> N3 N3       [IM=ISMN+8]
      IM   = ISMN+8
      IFF  = 34
      MJ   = MN_H(3)
      MK   = MN_H(3)
      CF   = 1.D0  ! Color Factor
      SYMF = 2.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 1.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*---> H_IH -> N3 N4       [IM=ISMN+9]
      IM   = ISMN+9
      IFF  = 55
      MJ   = MN_H(3)
      MK   = MN_H(4)
      CF   = 1.D0  ! Color Factor
      SYMF = 2.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*---> H_IH -> N4 N4       [IM=ISMN+10]
      IM   = ISMN+10
      IFF  = 37
      MJ   = MN_H(4)
      MK   = MN_H(4)
      CF   = 1.D0  ! Color Factor
      SYMF = 2.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 1.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*
*      print*,'H',IH,' -> N1 N1 : ',gambrn(ismn+1,1,ih)
*      print*,'H',IH,' -> N1 N2 : ',gambrn(ismn+2,1,ih)
*      print*,'H',IH,' -> N1 N3 : ',gambrn(ismn+3,1,ih)
*      print*,'H',IH,' -> N1 N4 : ',gambrn(ismn+4,1,ih)
*      print*,'H',IH,' -> N2 N2 : ',gambrn(ismn+5,1,ih)
*      print*,'H',IH,' -> N2 N3 : ',gambrn(ismn+6,1,ih)
*      print*,'H',IH,' -> N2 N4 : ',gambrn(ismn+7,1,ih)
*      print*,'H',IH,' -> N3 N3 : ',gambrn(ismn+8,1,ih)
*      print*,'H',IH,' -> N3 N4 : ',gambrn(ismn+9,1,ih)
*      print*,'H',IH,' -> N4 N4 : ',gambrn(ismn+10,1,ih)
*
*---> H_IH -> C1+ C1-       [IM=ISMN+11]
      IM   = ISMN+11
      IFF  = 58
      MJ   = MC_H(1)
      MK   = MC_H(1)
      CF   = 1.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*---> H_IH -> C1+ C2-       [IM=ISMN+11]
      IM   = ISMN+12
      IFF  = 64
      MJ   = MC_H(1)
      MK   = MC_H(2)
      CF   = 1.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*---> H_IH -> C2+ C1-       [IM=ISMN+11]
      IM   = ISMN+13
      IFF  = 61
      MJ   = MC_H(2)
      MK   = MC_H(1)
      CF   = 1.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*---> H_IH -> C2+ C2-       [IM=ISMN+11]
      IM   = ISMN+14
      IFF  = 67
      MJ   = MC_H(2)
      MK   = MC_H(2)
      CF   = 1.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,NHC_H(IFF,IH),NHC_H(IFF+1,IH)
     .,NHC_H(IFF+2,IH),HMASS_H(IH),MJ,MK,GAMBRN(IM,1,IH))
*
*      print*,'H',IH,' -> C1+ C1- : ',gambrn(ismn+11,1,ih)
*      print*,'H',IH,' -> C1+ C2- : ',gambrn(ismn+12,1,ih)
*      print*,'H',IH,' -> C2+ C1- : ',gambrn(ismn+13,1,ih)
*      print*,'H',IH,' -> C2+ C2- : ',gambrn(ismn+14,1,ih)
*
*---> H_IH -> ST1* ST1     [IM=ISMN+15]
      IM   = ISMN+15
      ISS  = 71
      SYMF = 3.D0 ! Color factor
      CALL HSS(SYMF,V_H,NHC_H(ISS,IH),HMASS_H(IH)
     .,STMASS_H(1),STMASS_H(1),GAMBRN(IM,1,IH))
*---> H_IH -> ST1* ST2     [IM=ISMN+16]
      IM   = ISMN+16
      ISS  = 73
      SYMF = 3.D0 ! Color factor
      CALL HSS(SYMF,V_H,NHC_H(ISS,IH),HMASS_H(IH)
     .,STMASS_H(1),STMASS_H(2),GAMBRN(IM,1,IH))
*---> H_IH -> ST2* ST1     [IM=ISMN+17]
      IM   = ISMN+17
      ISS  = 72
      SYMF = 3.D0 ! Color factor
      CALL HSS(SYMF,V_H,NHC_H(ISS,IH),HMASS_H(IH)
     .,STMASS_H(2),STMASS_H(1),GAMBRN(IM,1,IH))
*---> H_IH -> ST2* ST2     [IM=ISMN+18]
      IM   = ISMN+18
      ISS  = 74
      SYMF = 3.D0 ! Color factor
      CALL HSS(SYMF,V_H,NHC_H(ISS,IH),HMASS_H(IH)
     .,STMASS_H(2),STMASS_H(2),GAMBRN(IM,1,IH))
*---> H_IH -> SB1* SB1     [IM=ISMN+19]
      IM   = ISMN+19
      ISS  = 75
      SYMF = 3.D0 ! Color factor
      CALL HSS(SYMF,V_H,NHC_H(ISS,IH),HMASS_H(IH)
     .,SBMASS_H(1),SBMASS_H(1),GAMBRN(IM,1,IH))
*---> H_IH -> SB1* SB2     [IM=ISMN+20]
      IM   = ISMN+20
      ISS  = 77
      SYMF = 3.D0 ! Color factor
      CALL HSS(SYMF,V_H,NHC_H(ISS,IH),HMASS_H(IH)
     .,SBMASS_H(1),SBMASS_H(2),GAMBRN(IM,1,IH))
*---> H_IH -> SB2* SB1     [IM=ISMN+21]
      IM   = ISMN+21
      ISS  = 76
      SYMF = 3.D0 ! Color factor
      CALL HSS(SYMF,V_H,NHC_H(ISS,IH),HMASS_H(IH)
     .,SBMASS_H(2),SBMASS_H(1),GAMBRN(IM,1,IH))
*---> H_IH -> ST2* ST2     [IM=ISMN+22]
      IM   = ISMN+22
      ISS  = 78
      SYMF = 3.D0 ! Color factor
      CALL HSS(SYMF,V_H,NHC_H(ISS,IH),HMASS_H(IH)
     .,SBMASS_H(2),SBMASS_H(2),GAMBRN(IM,1,IH))
*JSL[01/SEP/05] Bug related STAUMASS_H fixed : Thanks to G. Belanger and S. Pukhov
*---> H_IH -> STAU1* STAU1 [IM=ISMN+23]
      IM   = ISMN+23
      ISS  = 79
      SYMF = 1.D0 ! Color factor
      CALL HSS(SYMF,V_H,NHC_H(ISS,IH),HMASS_H(IH)
     .,STAUMASS_H(1),STAUMASS_H(1),GAMBRN(IM,1,IH))
*---> H_IH -> STAU1* STAU2 [IM=ISMN+24]
      IM   = ISMN+24
      ISS  = 80
      SYMF = 1.D0 ! Color factor
      CALL HSS(SYMF,V_H,NHC_H(ISS,IH),HMASS_H(IH)
     .,STAUMASS_H(1),STAUMASS_H(2),GAMBRN(IM,1,IH))
*---> H_IH -> STAU2* STAU1 [IM=ISMN+25]
      IM   = ISMN+25
      ISS  = 81
      SYMF = 1.D0 ! Color factor
      CALL HSS(SYMF,V_H,NHC_H(ISS,IH),HMASS_H(IH)
     .,STAUMASS_H(2),STAUMASS_H(1),GAMBRN(IM,1,IH))
*---> H_IH -> STAU2* STAU2 [IM=ISMN+26]
      IM   = ISMN+26
      ISS  = 82
      SYMF = 1.D0 ! Color factor
      CALL HSS(SYMF,V_H,NHC_H(ISS,IH),HMASS_H(IH)
     .,STAUMASS_H(2),STAUMASS_H(2),GAMBRN(IM,1,IH))
*---> H_IH -> SNU3* SNU3   [IM=ISMN+27]
      IM   = ISMN+27
      ISS  = 83
      SYMF = 1.D0 ! Color factor
      CALL HSS(SYMF,V_H,NHC_H(ISS,IH),HMASS_H(IH)
     .,SNU3MASS_H,SNU3MASS_H,GAMBRN(IM,1,IH))
*
*      print*,'H',IH,' -> ST1* ST1 : ',gambrn(ismn+15,1,ih)
*      print*,'H',IH,' -> ST1* ST2 : ',gambrn(ismn+16,1,ih)
*      print*,'H',IH,' -> ST2* ST1 : ',gambrn(ismn+17,1,ih)
*      print*,'H',IH,' -> ST2* ST2 : ',gambrn(ismn+18,1,ih)
*      print*,'H',IH,' -> SB1* SB1 : ',gambrn(ismn+19,1,ih)
*      print*,'H',IH,' -> SB1* SB2 : ',gambrn(ismn+20,1,ih)
*      print*,'H',IH,' -> SB2* SB1 : ',gambrn(ismn+21,1,ih)
*      print*,'H',IH,' -> SB2* SB2 : ',gambrn(ismn+22,1,ih)
*
      ENDDO ! IH
*-----------------------------------------------------------------------
* 
* << BRANCHING RATIOS OF NEUTRAL HIGGS BOSON >>
*
      DO IH=1,3
*
       GAMBRN(ISMN,1,IH)=0.D0
       DO IM=1,ISMN-1
       GAMBRN(ISMN,1,IH)=GAMBRN(ISMN,1,IH)+GAMBRN(IM,1,IH)
       ENDDO
*
       IXX=ISMN+ISUSYN
       GAMBRN(IXX,1,IH)=0.D0
       DO IM=ISMN+1,IXX-1
       GAMBRN(IXX,1,IH)=GAMBRN(IXX,1,IH)+GAMBRN(IM,1,IH)
       ENDDO
*
       GAMBRN(NMNH,1,IH)=GAMBRN(ISMN,1,IH)+GAMBRN(IXX,1,IH)
*
       DO IM=1,ISMN+ISUSYN+1
        GAMBRN(IM,2,IH)=GAMBRN(IM,1,IH)/GAMBRN(ISMN,1,IH)
        GAMBRN(IM,3,IH)=GAMBRN(IM,1,IH)/GAMBRN(NMNH,1,IH)
       ENDDO
*
      ENDDO ! IH
*-----------------------------------------------------------------------
      IF(IFLAG_H(6).EQ.1) CALL DUMP_NHDCY(ISMN,ISUSYN,NMNH,GAMBRN,1)
      IF(IFLAG_H(6).EQ.2) CALL DUMP_NHDCY(ISMN,ISUSYN,NMNH,GAMBRN,2)
      IF(IFLAG_H(6).EQ.3) CALL DUMP_NHDCY(ISMN,ISUSYN,NMNH,GAMBRN,3)
      IF(IFLAG_H(6).EQ.5) CALL DUMP_NHDCY(ISMN,ISUSYN,NMNH,GAMBRN,5)
*-----------------------------------------------------------------------
*
* << CHARGED HIGGS BOSON DECAYS INTO SM PARTICLES >>
*

*---> running alpha_s and b-quark mass at Charged Higgs Mass scale
*     : MH > MB^pole assumed
      B5      = (11.D0-2.D0/3.D0*5.D0)/4.D0/PI
      B6      = (11.D0-2.D0/3.D0*6.D0)/4.D0/PI
      ASMZ    = ASMZ_H
      ASMT    = ASMT_H
      MH      = MCH
      IF(MH.LE.MTPOLE_H) THEN                               ! AS(MH)
       ASMH   = ASMZ/(1.D0+B5*ASMZ*DLOG(MH**2/MZ_H**2))
      ELSE
       ASMH   = ASMT/(1.D0+B6*ASMT*DLOG(MH**2/MTPOLE_H**2))
      ENDIF
      IF(MH.LE.MTPOLE_H) THEN                               ! MQ(MH)
       MTMH   = MTMT_H*(ASMH/ASMT)**(1.D0/B5/PI)
       MBMH   = MBMT_H*(ASMH/ASMT)**(1.D0/B5/PI)
       MCMH   = MCMT_H*(ASMH/ASMT)**(1.D0/B5/PI)
       MSMH   = MSMT_H*(ASMH/ASMT)**(1.D0/B5/PI)
       MUMH   = MUMT_H*(ASMH/ASMT)**(1.D0/B5/PI)
       MDMH   = MDMT_H*(ASMH/ASMT)**(1.D0/B5/PI)
      ELSE
       MTMH   = MTMT_H*(ASMH/ASMT)**(1.D0/B6/PI)
       MBMH   = MBMT_H*(ASMH/ASMT)**(1.D0/B6/PI)
       MCMH   = MCMT_H*(ASMH/ASMT)**(1.D0/B6/PI)
       MSMH   = MSMT_H*(ASMH/ASMT)**(1.D0/B6/PI)
       MUMH   = MUMT_H*(ASMH/ASMT)**(1.D0/B6/PI)
       MDMH   = MDMT_H*(ASMH/ASMT)**(1.D0/B6/PI)
      ENDIF
*      print*,'MCH,AS(MCH),=',MCH,ASMH,MTMH,MBMH,MCMH,MSMH,MUMH,MDMH
*
*JSL 10/Jun/2009: Including threshold corrections
*       R_123Q=SSPARA_H(22)
*       R_123D=SSPARA_H(24)
*       MQ12  =R_123Q*MQ3_H
*       MD12  =R_123D*MD3_H
*       Q12SQ =DMAX1(MQ12**2,MD12**2)
*       AS_M12=ASMT/(1.D0+B6*ASMT*DLOG(Q12SQ/MTPOLE_H**2))
*       CKS_H=2.D0*AS_M12/3.D0/PI*DCONJG(MU_H*M3_H)
*     .           *F_I(MD12**2,MQ12**2,CDABS(M3_H)**2)
       CKS_H=CAUX_H(11)  ! from FILLCOUPL2
       CKD_H=CKS_H
*       print*,'>>> FILLGAMBR2: cks_h',cks_h
*---> CH+ -> e+ nu     [IM= 1]
      IM   = 1
      IFF  = 1
      MJ   = ME_H
      MK   = 0.D0
      CF   = 1.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,CHC_H(IFF),CHC_H(IFF+1)
     .,CHC_H(IFF+2),MCH,MJ,MK,GAMBRC(IM,1))
*---> CH+ -> mu+ nu     [IM= 2]
      IM   = 2
      IFF  = 4
      MJ   = MMU_H
      MK   = 0.D0
      CF   = 1.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,CHC_H(IFF),CHC_H(IFF+1)
     .,CHC_H(IFF+2),MCH,MJ,MK,GAMBRC(IM,1))
*---> CH+ -> tau+ nu     [IM= 3]
      IM   = 3
      IFF  = 7
      MJ   = MTAU_H
      MK   = 0.D0
      CF   = 1.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,CHC_H(IFF),CHC_H(IFF+1)
     .,CHC_H(IFF+2),MCH,MJ,MK,GAMBRC(IM,1))
*---> CH+ -> u d         [IM= 4] 
      IM   = 4
      IFF  = 10
      GF_UD=DCMPLX(-GW_H*MUMH/DSQRT(2.D0)/MW_H,0.D0)
*      GS_UD=DCMPLX((1.D0/TB_H+MDMH/MUMH*TB_H)/2.D0,0.D0)
*      GP_UD=DCMPLX(0.D0,(1.D0/TB_H-MDMH/MUMH*TB_H)/2.D0)
*JSL 10/Jun/2009: Including threshold corrections
      GS_UD=(1.D0/TB_H+TB_H/(1.D0+DCONJG(CKD_H)*TB_H)*MDMH/MUMH)/2.D0
      GP_UD=(1.D0/TB_H-TB_H/(1.D0+DCONJG(CKD_H)*TB_H)*MDMH/MUMH)*XI/2.D0
*       print*,'FILLGAMBR2',gs_ud,chc_h(11)
*       print*,'FILLGAMBR2',gp_ud,chc_h(12)
      MJ   = MUMH
      MK   = MDMH
      CF   = 3.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,GF_UD,GS_UD,GP_UD,MCH,MJ,MK,GAMBRC(IM,1))
*K-factor
      GAMBRC(IM,1)=(1.D0+5.67D0*ASMH/PI)*GAMBRC(IM,1)
*---> CH+ -> c s         [IM= 5]
      IM   = 5
      IFF  = 13
      GF_CS=DCMPLX(-GW_H*MCMH/DSQRT(2.D0)/MW_H,0.D0)
*      GS_CS=DCMPLX((1.D0/TB_H+MSMH/MCMH*TB_H)/2.D0,0.D0)
*      GP_CS=DCMPLX(0.D0,(1.D0/TB_H-MSMH/MCMH*TB_H)/2.D0)
*JSL 10/Jun/2009: Including threshold corrections
      GS_CS=(1.D0/TB_H+TB_H/(1.D0+DCONJG(CKS_H)*TB_H)*MSMH/MCMH)/2.D0
      GP_CS=(1.D0/TB_H-TB_H/(1.D0+DCONJG(CKS_H)*TB_H)*MSMH/MCMH)*XI/2.D0
*       print*,'FILLGAMBR2',gs_cs,chc_h(14)
*       print*,'FILLGAMBR2',gp_cs,chc_h(15)
*       print*,'FILLGAMBR2',msmh,mcmh,msmh/mcmh
      MJ   = MCMH
      MK   = MSMH
      CF   = 3.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
      CALL HFF(CF,SYMF,DJK,GF_CS,GS_CS,GP_CS,MCH,MJ,MK,GAMBRC(IM,1))
*K-factor
      GAMBRC(IM,1)=(1.D0+5.67D0*ASMH/PI)*GAMBRC(IM,1)
*---> CH+ -> t b         [IM= 6]
      CALL RADNHTB(NFLAG,IFLAG_H,SBMASS_H,STMASS_H
     .            ,HB_H,HT_H,CKB_H,CKT_H)
*       print*,'FILLGAMBR2:CKB_H',ckb_h
      CALL RADCHTB(NFLAG,IFLAG_H,SBMIX_H,SBMASS_H,STMIX_H,STMASS_H
     .            ,RB_H,RT_H,CKBB_H,CKBT_H)
*      print*,'>> FILLGAMBR'
*      print*,DSQRT(2.D0)*MBMT_H/V_H/CB_H,DSQRT(2.D0)*MTMT_H/V_H/SB_H
*      print*,abs(hb_h),abs(ht_h)
*      print*,hb_h,ht_h
*      print*,ckb_h,ckt_h
*      print*,rb_h,rt_h
*      print*,ckbb_h,ckbt_h
      IM   = 6
      IFF  = 16
      XI   =DCMPLX(0.D0,1.D0)
      GF_TB=DCMPLX(-GW_H*MTMH/DSQRT(2.D0)/MW_H,0.D0)
      GS_TB=((1.D0/TB_H*(1.D0+RT_H)-CKBT_H)/(1.D0+CKT_H/TB_H)
     .          +(TB_H*(1.D0+DCONJG(RB_H))-DCONJG(CKBB_H))
     .           /(1.D0+DCONJG(CKB_H)*TB_H)*MBMH/MTMH)/2.D0
      GP_TB=((1.D0/TB_H*(1.D0+RT_H)-CKBT_H)/(1.D0+CKT_H/TB_H)
     .          -(TB_H*(1.D0+DCONJG(RB_H))-DCONJG(CKBB_H))
     .           /(1.D0+DCONJG(CKB_H)*TB_H)*MBMH/MTMH)*XI/2.D0
*      print*,mtmh,mtmt_h,mtmh/mtmt_h
*      print*,mbmh,mbmt_h,mbmh/mbmt_h
*      print*,chc_h(16),gf_tb,chc_h(16)*mtmh/mtmt_h
*      print*,chc_h(17),gs_tb
*      print*,chc_h(18),gp_tb
      MJ   = MTMH
      MK   = MBMH
      CF   = 3.D0  ! Color Factor
      SYMF = 1.D0  ! Sym. Factor [1=Dirac, 2=Majorana]
      DJK  = 0.D0  ! 1 For identical Majorana particles
*JSLEE 31/MAR/2005 : Three body decay has been included CH -> t* b -> W b b
*JSLEE 10/JUN/2009 : For LO QCD corrections, see, for example, hep-ph/9910339
*      print*,'delta_0',
*     . (1.D0-MW_H**2/MTPOLE_H**2)**2*(1.D0+2.D0*MW_H**2/MTPOLE_H**2) 
      GTWB=GW_H**2*MTPOLE_H**3/64.D0/PI/MW_H**2
     .    *( (1.D0-MW_H**2/MTPOLE_H**2)**2
     .      *(1.D0+2.D0*MW_H**2/MTPOLE_H**2) 
     .      -2.20D0*ASMT_H/PI )
      IF(MTPOLE_H.GT.MCH+MBMH) THEN
       GTHB=CDABS(GF_TB)**2*MTPOLE_H/16.D0/PI
     .    *(CDABS(GS_TB)**2+CDABS(GP_TB)**2)
     .    *(1.D0-MCH**2/MTPOLE_H**2)**2
      ELSE
       GTHB=0.D0
      ENDIF
*JSLEE 10/JUN/2009 : Store top-quark decay widths
       RAUX_H(50)=GTWB
       RAUX_H(51)=GTHB
       RAUX_H(52)=GTWB/(GTWB+GTHB)
       RAUX_H(53)=GTHB/(GTWB+GTHB)
*Top quark decay width
*      GTOP=GTWB+GTHB 
*H^\pm-b loop does not contrubute to the absorptive part of top self-energy
      GTOP=GTWB
*      print*,gtwb,gthb,gtop
      NSTEP=500
*
      IF(MCH.LE.MW_H+2.D0*MK) THEN
       GAMBRC(IM,1)=0.D0
      ELSEIF(MCH.LE.MJ+MK-2.D0) THEN
       CALL CHWBB(NSTEP,GW_H,GF_TB,GS_TB,GP_TB
     .           ,MCH,TB_H,MJ,MW_H,MK,GTOP,GCHWBB)
       GAMBRC(IM,1)=GCHWBB
      ELSEIF(MCH.LE.MJ+MK+2.D0) THEN
*
        MCH1=MJ+MK-3.D0
        ASMH1=ASMZ/(1.D0+B5*ASMZ*DLOG(MCH1**2/MZ_H**2))
        MT1  = MTMT_H*(ASMH1/ASMT)**(1.D0/B5/PI)
        MB1  = MBMT_H*(ASMH1/ASMT)**(1.D0/B5/PI)
        GF1=DCMPLX(-GW_H*MT1/DSQRT(2.D0)/MW_H,0.D0)
        GS1=((1.D0/TB_H*(1.D0+RT_H)-CKBT_H)/(1.D0+CKT_H/TB_H)
     .          +(TB_H*(1.D0+DCONJG(RB_H))-DCONJG(CKBB_H))
     .           /(1.D0+DCONJG(CKB_H)*TB_H)*MB1/MT1)/2.D0
        GP1=((1.D0/TB_H*(1.D0+RT_H)-CKBT_H)/(1.D0+CKT_H/TB_H)
     .          -(TB_H*(1.D0+DCONJG(RB_H))-DCONJG(CKBB_H))
     .           /(1.D0+DCONJG(CKB_H)*TB_H)*MB1/MT1)*XI/2.D0
       CALL CHWBB(NSTEP,GW_H,GF1,GS1,GP1
     .           ,MCH1,TB_H,MT1,MW_H,MB1,GTOP,GCHWBB1)
*
        MCH2=MJ+MK-2.D0
        ASMH2=ASMZ/(1.D0+B5*ASMZ*DLOG(MCH2**2/MZ_H**2))
        MT2  = MTMT_H*(ASMH2/ASMT)**(1.D0/B5/PI)
        MB2  = MBMT_H*(ASMH2/ASMT)**(1.D0/B5/PI)
        GF2=DCMPLX(-GW_H*MT2/DSQRT(2.D0)/MW_H,0.D0)
        GS2=((1.D0/TB_H*(1.D0+RT_H)-CKBT_H)/(1.D0+CKT_H/TB_H)
     .          +(TB_H*(1.D0+DCONJG(RB_H))-DCONJG(CKBB_H))
     .           /(1.D0+DCONJG(CKB_H)*TB_H)*MB2/MT2)/2.D0
        GP2=((1.D0/TB_H*(1.D0+RT_H)-CKBT_H)/(1.D0+CKT_H/TB_H)
     .          -(TB_H*(1.D0+DCONJG(RB_H))-DCONJG(CKBB_H))
     .           /(1.D0+DCONJG(CKB_H)*TB_H)*MB2/MT2)*XI/2.D0
       CALL CHWBB(NSTEP,GW_H,GF2,GS2,GP2
     .           ,MCH2,TB_H,MT2,MW_H,MB2,GTOP,GCHWBB2)
*
        MCH3=MJ+MK+2.D0
        ASMH3=ASMZ/(1.D0+B5*ASMZ*DLOG(MCH3**2/MZ_H**2))
        MT3  = MTMT_H*(ASMH3/ASMT)**(1.D0/B5/PI)
        MB3  = MBMT_H*(ASMH3/ASMT)**(1.D0/B5/PI)
        GF3=DCMPLX(-GW_H*MT3/DSQRT(2.D0)/MW_H,0.D0)
        GS3=((1.D0/TB_H*(1.D0+RT_H)-CKBT_H)/(1.D0+CKT_H/TB_H)
     .          +(TB_H*(1.D0+DCONJG(RB_H))-DCONJG(CKBB_H))
     .           /(1.D0+DCONJG(CKB_H)*TB_H)*MB3/MT3)/2.D0
        GP3=((1.D0/TB_H*(1.D0+RT_H)-CKBT_H)/(1.D0+CKT_H/TB_H)
     .          -(TB_H*(1.D0+DCONJG(RB_H))-DCONJG(CKBB_H))
     .           /(1.D0+DCONJG(CKB_H)*TB_H)*MB3/MT3)*XI/2.D0
       CALL HFF(CF,SYMF,DJK,GF3,GS3,GP3,MCH3,MT3,MB3,GCHWBB3)
        GCHWBB3=(1.D0+5.67D0*ASMH3/PI)*GCHWBB3
*
        MCH4=MJ+MK+3.D0
        ASMH4=ASMZ/(1.D0+B5*ASMZ*DLOG(MCH4**2/MZ_H**2))
        MT4  = MTMT_H*(ASMH4/ASMT)**(1.D0/B5/PI)
        MB4  = MBMT_H*(ASMH4/ASMT)**(1.D0/B5/PI)
        GF4=DCMPLX(-GW_H*MT4/DSQRT(2.D0)/MW_H,0.D0)
        GS4=((1.D0/TB_H*(1.D0+RT_H)-CKBT_H)/(1.D0+CKT_H/TB_H)
     .          +(TB_H*(1.D0+DCONJG(RB_H))-DCONJG(CKBB_H))
     .           /(1.D0+DCONJG(CKB_H)*TB_H)*MB4/MT4)/2.D0
        GP4=((1.D0/TB_H*(1.D0+RT_H)-CKBT_H)/(1.D0+CKT_H/TB_H)
     .          -(TB_H*(1.D0+DCONJG(RB_H))-DCONJG(CKBB_H))
     .           /(1.D0+DCONJG(CKB_H)*TB_H)*MB4/MT4)*XI/2.D0
       CALL HFF(CF,SYMF,DJK,GF4,GS4,GP4,MCH4,MT4,MB4,GCHWBB4)
        GCHWBB4=(1.D0+5.67D0*ASMH4/PI)*GCHWBB4
*
*       print*,mch,mj,mk,gtop
*       print*,mj+mk-3.d0,gchwbb1
*       print*,mj+mk-2.d0,gchwbb2
*       print*,mj+mk+2.d0,gchwbb3
*       print*,mj+mk+3.d0,gchwbb4
*Linear extrapolation
       DMCH=MCH-(MJ+MK-2.D0)
       GAMBRC(IM,1)=GCHWBB2+(GCHWBB3-GCHWBB2)*DMCH/4.D0
      ELSE
       CALL HFF(CF,SYMF,DJK,GF_TB,GS_TB,GP_TB,MCH,MJ,MK,GAMBRC(IM,1))
       GAMBRC(IM,1)=(1.D0+5.67D0*ASMH/PI)*GAMBRC(IM,1)
      ENDIF
*      print*,mch,mj,mk,gtop
*      print*,mch,gambrc(im,1)
*      MJ   = MTMT_H
*      MK   = MBMT_H
*      CALL HFF(CF,SYMF,DJK,CHC_H(IFF),CHC_H(IFF+1)
*     .,CHC_H(IFF+2),MCH,MJ,MK,GAMBRC(IM,1))
*      print*,gambrc(im,1)
*
*      print*,'CH+ -> e   nu  : ',gambrc(1,1)
*      print*,'CH+ -> mu  nu  : ',gambrc(2,1)
*      print*,'CH+ -> tau nu  : ',gambrc(3,1)
*      print*,'CH+ -> u   d   : ',gambrc(4,1)
*      print*,'CH+ -> c   s   : ',gambrc(5,1)
*      print*,'CH+ -> t   b   : ',gambrc(6,1)
*---> CH+ -> H1 W        [IM=7]
      IM  = 7
      IHV = 87
      CALL HHV(GF_H,NHC_H(IHV,1),MCH,HMASS_H(1),MW_H,GAMBRC(IM,1))
* CH+ -> H1 W*
      IF( (MCH.GT.HMASS_H(1)) .AND. 
     .    (MCH.LT.(HMASS_H(1)+MW_H+5.D0*GAMW_H)) ) THEN
      EPSV=GAMW_H/MW_H
      OMEGAI=MCH**2/MW_H**2
      OMEGAJ=HMASS_H(1)**2/MW_H**2
      XUP=(DSQRT(OMEGAI)-DSQRT(OMEGAJ))**2
      XDW=0.D0
      NSTEP=500
      CALL BODE(FHVS,0.D0,1.D0,NSTEP,RES)
      GAMBRC(IM,1)=GW_H**2*MCH*CDABS(NHC_H(IHV,1))**2 
     . /64.D0/PI**2*EPSV*(MW_H/MCH)**4*RES
      ENDIF
*---> CH+ -> H2 W        [IM=8]
      IM  = 8
      IHV = 87
      CALL HHV(GF_H,NHC_H(IHV,2),MCH,HMASS_H(2),MW_H,GAMBRC(IM,1))
* CH+ -> H2 W*
      IF( (MCH.GT.HMASS_H(2)) .AND. 
     .    (MCH.LT.(HMASS_H(2)+MW_H+5.D0*GAMW_H)) ) THEN
      EPSV=GAMW_H/MW_H
      OMEGAI=MCH**2/MW_H**2
      OMEGAJ=HMASS_H(2)**2/MW_H**2
      XUP=(DSQRT(OMEGAI)-DSQRT(OMEGAJ))**2
      XDW=0.D0
      NSTEP=500
      CALL BODE(FHVS,0.D0,1.D0,NSTEP,RES)
      GAMBRC(IM,1)=GW_H**2*MCH*CDABS(NHC_H(IHV,2))**2 
     . /64.D0/PI**2*EPSV*(MW_H/MCH)**4*RES
      ENDIF
*JSL: Added on Jun.06.2008
*---> CH+ -> H3 W        [IM=9]
      IM  = 9
      IHV = 87
      CALL HHV(GF_H,NHC_H(IHV,3),MCH,HMASS_H(3),MW_H,GAMBRC(IM,1))
* CH+ -> H3 W*
      IF( (MCH.GT.HMASS_H(3)) .AND.
     .    (MCH.LT.(HMASS_H(3)+MW_H+5.D0*GAMW_H)) ) THEN
      EPSV=GAMW_H/MW_H
      OMEGAI=MCH**2/MW_H**2
      OMEGAJ=HMASS_H(3)**2/MW_H**2
      XUP=(DSQRT(OMEGAI)-DSQRT(OMEGAJ))**2
      XDW=0.D0
      NSTEP=500
      CALL BODE(FHVS,0.D0,1.D0,NSTEP,RES)
      GAMBRC(IM,1)=GW_H**2*MCH*CDABS(NHC_H(IHV,3))**2
     . /64.D0/PI**2*EPSV*(MW_H/MCH)**4*RES
      ENDIF
*      print*,'CH+ -> H1   W   : ',gambrc(7,1)
*      print*,'CH+ -> H2   W   : ',gambrc(8,1)
*      print*,'CH+ -> H3   W   : ',gambrc(9,1)
*-----------------------------------------------------------------------
*
* << CHARGED HIGGS BOSON DECAYS INTO SUSY PARTICLES >>
 
*---> CH+ -> N1 C1+      [IM=ISMC+1]
      IM   = ISMC+1
      IFF  = 19
      MJ   = MN_H(1)
      MK   = MC_H(1)
      CALL HFF(1.D0,1.D0,0.D0,CHC_H(IFF),CHC_H(IFF+1)
     .,CHC_H(IFF+2),MCH,MJ,MK,GAMBRC(IM,1))
*---> CH+ -> N2 C1+      [IM=ISMC+2]
      IM   = ISMC+2
      IFF  = 25
      MJ   = MN_H(2)
      MK   = MC_H(1)
      CALL HFF(1.D0,1.D0,0.D0,CHC_H(IFF),CHC_H(IFF+1)
     .,CHC_H(IFF+2),MCH,MJ,MK,GAMBRC(IM,1))
*---> CH+ -> N3 C1+      [IM=ISMC+3]
      IM   = ISMC+3
      IFF  = 31
      MJ   = MN_H(3)
      MK   = MC_H(1)
      CALL HFF(1.D0,1.D0,0.D0,CHC_H(IFF),CHC_H(IFF+1)
     .,CHC_H(IFF+2),MCH,MJ,MK,GAMBRC(IM,1))
*---> CH+ -> N4 C1+      [IM=ISMC+4]
      IM   = ISMC+4
      IFF  = 37
      MJ   = MN_H(4)
      MK   = MC_H(1)
      CALL HFF(1.D0,1.D0,0.D0,CHC_H(IFF),CHC_H(IFF+1)
     .,CHC_H(IFF+2),MCH,MJ,MK,GAMBRC(IM,1))
*---> CH+ -> N1 C2+      [IM=ISMC+5]
      IM   = ISMC+5
      IFF  = 22
      MJ   = MN_H(1)
      MK   = MC_H(2)
      CALL HFF(1.D0,1.D0,0.D0,CHC_H(IFF),CHC_H(IFF+1)
     .,CHC_H(IFF+2),MCH,MJ,MK,GAMBRC(IM,1))
*---> CH+ -> N2 C2+      [IM=ISMC+6]
      IM   = ISMC+6
      IFF  = 28
      MJ   = MN_H(2)
      MK   = MC_H(2)
      CALL HFF(1.D0,1.D0,0.D0,CHC_H(IFF),CHC_H(IFF+1)
     .,CHC_H(IFF+2),MCH,MJ,MK,GAMBRC(IM,1))
*---> CH+ -> N3 C2+      [IM=ISMC+7]
      IM   = ISMC+7
      IFF  = 34
      MJ   = MN_H(3)
      MK   = MC_H(2)
      CALL HFF(1.D0,1.D0,0.D0,CHC_H(IFF),CHC_H(IFF+1)
     .,CHC_H(IFF+2),MCH,MJ,MK,GAMBRC(IM,1))
*---> CH+ -> N4 C2+      [IM=ISMC+8]
      IM   = ISMC+8
      IFF  = 40
      MJ   = MN_H(4)
      MK   = MC_H(2)
      CALL HFF(1.D0,1.D0,0.D0,CHC_H(IFF),CHC_H(IFF+1)
     .,CHC_H(IFF+2),MCH,MJ,MK,GAMBRC(IM,1))
*---> CH+ -> stop1 sbottom1*   [IM=ISMC+9]
      IM   = ISMC+9
      ISS  = 43
      SYMF = 3.D0 ! Color factor
      CALL HSS(SYMF,V_H,CHC_H(ISS),MCH
     .,STMASS_H(1),SBMASS_H(1),GAMBRC(IM,1))
*---> CH+ -> stop1 sbottom2*   [IM=ISMC+10]
      IM   = ISMC+10
      ISS  = 44
      SYMF = 3.D0 ! Color factor
      CALL HSS(SYMF,V_H,CHC_H(ISS),MCH
     .,STMASS_H(1),SBMASS_H(2),GAMBRC(IM,1))
*---> CH+ -> stop2 sbottom1*   [IM=ISMC+11]
      IM   = ISMC+11
      ISS  = 45
      SYMF = 3.D0 ! Color factor
      CALL HSS(SYMF,V_H,CHC_H(ISS),MCH
     .,STMASS_H(2),SBMASS_H(1),GAMBRC(IM,1))
*---> CH+ -> stop2 sbottom2*   [IM=ISMC+12]
      IM   = ISMC+12
      ISS  = 46
      SYMF = 3.D0 ! Color factor
      CALL HSS(SYMF,V_H,CHC_H(ISS),MCH
     .,STMASS_H(2),SBMASS_H(2),GAMBRC(IM,1))
*---> CH+ -> snu3 stau1*       [IM=ISMC+13]
      IM   = ISMC+13
      ISS  = 47
      SYMF = 1.D0 ! Color factor
      CALL HSS(SYMF,V_H,CHC_H(ISS),MCH
     .,SNU3MASS_H,STAUMASS_H(1),GAMBRC(IM,1))
*---> CH+ -> snu3 stau2*       [IM=ISMC+14]
      IM   = ISMC+14
      ISS  = 48
      SYMF = 1.D0 ! Color factor
      CALL HSS(SYMF,V_H,CHC_H(ISS),MCH
     .,SNU3MASS_H,STAUMASS_H(2),GAMBRC(IM,1))

*      print*,'CH+ -> N1   C1+   : ',gambrc(ismc+1,1)
*      print*,'CH+ -> N2   C1+   : ',gambrc(ismc+2,1)
*      print*,'CH+ -> N3   C1+   : ',gambrc(ismc+3,1)
*      print*,'CH+ -> N4   C1+   : ',gambrc(ismc+4,1)
*      print*,'CH+ -> N1   C2+   : ',gambrc(ismc+5,1)
*      print*,'CH+ -> N2   C2+   : ',gambrc(ismc+6,1)
*      print*,'CH+ -> N3   C2+   : ',gambrc(ismc+7,1)
*      print*,'CH+ -> N4   C2+   : ',gambrc(ismc+8,1)
*
*-----------------------------------------------------------------------
* 
* << BRANCHING RATIOS OF CHARGED HIGGS BOSON >>
*
       GAMBRC(ISMC,1)=0.D0
       DO IM=1,ISMC-1
       GAMBRC(ISMC,1)=GAMBRC(ISMC,1)+GAMBRC(IM,1)
       ENDDO
*
       IXX=ISMC+ISUSYC
       GAMBRC(IXX,1)=0.D0
       DO IM=ISMC+1,IXX-1
       GAMBRC(IXX,1)=GAMBRC(IXX,1)+GAMBRC(IM,1)
       ENDDO
*
       GAMBRC(NMCH,1)=GAMBRC(ISMC,1)+GAMBRC(IXX,1)
*
       DO IM=1,ISMC+ISUSYC+1
        GAMBRC(IM,2)=GAMBRC(IM,1)/GAMBRC(ISMC,1)
        GAMBRC(IM,3)=GAMBRC(IM,1)/GAMBRC(NMCH,1)
       ENDDO
*
*-----------------------------------------------------------------------
      IF(IFLAG_H(6).EQ.4.OR.IFLAG_H(6).EQ.5) 
     . CALL DUMP_CHDCY(ISMC,ISUSYC,NMCH,GAMBRC)
*-----------------------------------------------------------------------
      RETURN
      END

      SUBROUTINE DUMP_CHDCY(ISMC,ISUSYC,NMCH,GAMBRC)
************************************************************************
*
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      INTEGER*8 ISMC,ISUSYC
      REAL*8 GAMBRC(NMCH,3)   
*
      print*,'---------------------------------------------------------'
      print*,'Charged Higgs Boson Decays with '
     .,'ISMC = ',ISMC,' : ISUSYC = ',ISUSYC
      print*,'---------------------------------------------------------'
      print*,'DECAY MODE      [IM]   WIDTH[GeV]  BR[SM]      BR[TOTAL]'
      print*,'---------------------------------------------------------'
      WRITE(*,1)GAMBRC(1,1),GAMBRC(1,2),GAMBRC(1,3)
      WRITE(*,2)GAMBRC(2,1),GAMBRC(2,2),GAMBRC(2,3)
      WRITE(*,3)GAMBRC(3,1),GAMBRC(3,2),GAMBRC(3,3)
      WRITE(*,4)GAMBRC(4,1),GAMBRC(4,2),GAMBRC(4,3)
      WRITE(*,5)GAMBRC(5,1),GAMBRC(5,2),GAMBRC(5,3)
      WRITE(*,6)GAMBRC(6,1),GAMBRC(6,2),GAMBRC(6,3)
      WRITE(*,7)GAMBRC(7,1),GAMBRC(7,2),GAMBRC(7,3)
      WRITE(*,8)GAMBRC(8,1),GAMBRC(8,2),GAMBRC(8,3)
      WRITE(*,9)GAMBRC(9,1),GAMBRC(9,2),GAMBRC(9,3)
      WRITE(*,25)GAMBRC(25,1),GAMBRC(25,2),GAMBRC(25,3)
      WRITE(*,26)GAMBRC(26,1),GAMBRC(26,2),GAMBRC(26,3)
      WRITE(*,27)GAMBRC(27,1),GAMBRC(27,2),GAMBRC(27,3)
      WRITE(*,28)GAMBRC(28,1),GAMBRC(28,2),GAMBRC(28,3)
      WRITE(*,29)GAMBRC(29,1),GAMBRC(29,2),GAMBRC(29,3)
      WRITE(*,30)GAMBRC(30,1),GAMBRC(30,2),GAMBRC(30,3)
      WRITE(*,31)GAMBRC(31,1),GAMBRC(31,2),GAMBRC(31,3)
      WRITE(*,32)GAMBRC(32,1),GAMBRC(32,2),GAMBRC(32,3)
      WRITE(*,33)GAMBRC(33,1),GAMBRC(33,2),GAMBRC(33,3)
      WRITE(*,34)GAMBRC(34,1),GAMBRC(34,2),GAMBRC(34,3)
      WRITE(*,35)GAMBRC(35,1),GAMBRC(35,2),GAMBRC(35,3)
      WRITE(*,36)GAMBRC(36,1),GAMBRC(36,2),GAMBRC(36,3)
      WRITE(*,37)GAMBRC(37,1),GAMBRC(37,2),GAMBRC(37,3)
      WRITE(*,38)GAMBRC(38,1),GAMBRC(38,2),GAMBRC(38,3)
      WRITE(*,39)GAMBRC(39,1),GAMBRC(39,2),GAMBRC(39,3)
      WRITE(*,50)GAMBRC(50,1),GAMBRC(50,2),GAMBRC(50,3)
      WRITE(*,51)GAMBRC(51,1),GAMBRC(51,2),GAMBRC(51,3)
      print*,' '
      print*,'* Note : WIDTH=GAMBRC(IM,1), BR[SM]   =GAMBRC(IM,2) '
      print*,'                      and    BR[TOTAL]=GAMBRC(IM,3) '
      print*,'---------------------------------------------------------'
*
  1   FORMAT(1X,'CH+ -> e+   nu  [ 1]:',3(2X,E10.4))
  2   FORMAT(1X,'CH+ -> mu+  nu  [ 2]:',3(2X,E10.4))
  3   FORMAT(1X,'CH+ -> tau+ nu  [ 3]:',3(2X,E10.4))
  4   FORMAT(1X,'CH+ -> u    d   [ 4]:',3(2X,E10.4))
  5   FORMAT(1X,'CH+ -> c    s   [ 5]:',3(2X,E10.4))
  6   FORMAT(1X,'CH+ -> t    b   [ 6]:',3(2X,E10.4))
  7   FORMAT(1X,'CH+ -> H1   W   [ 7]:',3(2X,E10.4))
  8   FORMAT(1X,'CH+ -> H2   W   [ 8]:',3(2X,E10.4))
  9   FORMAT(1X,'CH+ -> H3   W   [ 9]:',3(2X,E10.4))
 25   FORMAT(1X,'CH+ TOTAL(SM)   [25]:',3(2X,E10.4))
 26   FORMAT(1X,'CH+ -> N1   C1+ [26]:',3(2X,E10.4))
 27   FORMAT(1X,'CH+ -> N2   C1+ [27]:',3(2X,E10.4))
 28   FORMAT(1X,'CH+ -> N3   C1+ [28]:',3(2X,E10.4))
 29   FORMAT(1X,'CH+ -> N4   C1+ [29]:',3(2X,E10.4))
 30   FORMAT(1X,'CH+ -> N1   C2+ [30]:',3(2X,E10.4))
 31   FORMAT(1X,'CH+ -> N2   C2+ [31]:',3(2X,E10.4))
 32   FORMAT(1X,'CH+ -> N3   C2+ [32]:',3(2X,E10.4))
 33   FORMAT(1X,'CH+ -> N4   C2+ [33]:',3(2X,E10.4))
 34   FORMAT(1X,'CH+ -> ST1 SB1* [34]:',3(2X,E10.4))
 35   FORMAT(1X,'CH+ -> ST1 SB2* [35]:',3(2X,E10.4))
 36   FORMAT(1X,'CH+ -> ST2 SB1* [36]:',3(2X,E10.4))
 37   FORMAT(1X,'CH+ -> ST2 SB2* [37]:',3(2X,E10.4))
 38   FORMAT(1X,'CH+ ->SNU3 STA1*[38]:',3(2X,E10.4))
 39   FORMAT(1X,'CH+ ->SNU3 STA2*[39]:',3(2X,E10.4))
 50   FORMAT(1X,'CH+ TOTAL(SUSY) [50]:',3(2X,E10.4))
 51   FORMAT(1X,'CH+ TOTAL       [51]:',3(2X,E10.4))
*
      RETURN
      END

      SUBROUTINE DUMP_NHDCY(ISMN,ISUSYN,NMNH,GAMBRN,IPRI)
************************************************************************
*
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      INTEGER*8 ISMN,ISUSYN
      REAL*8 GAMBRN(NMNH,3,3)   
*
      IF(IPRI.LE.3) IH=IPRI
      IC = 0
 999  CONTINUE
      IC = IC+1
      IF(IPRI.EQ.5) IH=IC
*
      print*,'---------------------------------------------------------'
      print*,'Neutral Higgs Boson Decays with '
     .,'ISMN = ',ISMN,' : ISUSYN = ',ISUSYN
      print*,'---------------------------------------------------------'
      print*,'DECAY MODE    [ IM]   WIDTH[GeV]  BR[SM]      BR[TOTAL]'
      print*,'---------------------------------------------------------'
       WRITE(*,1)IH,GAMBRN(1,1,ih),GAMBRN(1,2,ih),GAMBRN(1,3,ih)
       WRITE(*,2)IH,GAMBRN(2,1,ih),GAMBRN(2,2,ih),GAMBRN(2,3,ih)
       WRITE(*,3)IH,GAMBRN(3,1,ih),GAMBRN(3,2,ih),GAMBRN(3,3,ih)
       WRITE(*,4)IH,GAMBRN(4,1,ih),GAMBRN(4,2,ih),GAMBRN(4,3,ih)
       WRITE(*,5)IH,GAMBRN(5,1,ih),GAMBRN(5,2,ih),GAMBRN(5,3,ih)
       WRITE(*,6)IH,GAMBRN(6,1,ih),GAMBRN(6,2,ih),GAMBRN(6,3,ih)
       WRITE(*,7)IH,GAMBRN(7,1,ih),GAMBRN(7,2,ih),GAMBRN(7,3,ih)
       WRITE(*,8)IH,GAMBRN(8,1,ih),GAMBRN(8,2,ih),GAMBRN(8,3,ih)
       WRITE(*,9)IH,GAMBRN(9,1,ih),GAMBRN(9,2,ih),GAMBRN(9,3,ih)
       WRITE(*,10)IH,GAMBRN(10,1,ih),GAMBRN(10,2,ih),GAMBRN(10,3,ih)
       WRITE(*,11)IH,GAMBRN(11,1,ih),GAMBRN(11,2,ih),GAMBRN(11,3,ih)
       WRITE(*,12)IH,GAMBRN(12,1,ih),GAMBRN(12,2,ih),GAMBRN(12,3,ih)
       WRITE(*,13)IH,GAMBRN(13,1,ih),GAMBRN(13,2,ih),GAMBRN(13,3,ih)
       WRITE(*,14)IH,GAMBRN(14,1,ih),GAMBRN(14,2,ih),GAMBRN(14,3,ih)
       WRITE(*,15)IH,GAMBRN(15,1,ih),GAMBRN(15,2,ih),GAMBRN(15,3,ih)
       WRITE(*,16)IH,GAMBRN(16,1,ih),GAMBRN(16,2,ih),GAMBRN(16,3,ih)
       WRITE(*,17)IH,GAMBRN(17,1,ih),GAMBRN(17,2,ih),GAMBRN(17,3,ih)
       WRITE(*,18)IH,GAMBRN(18,1,ih),GAMBRN(18,2,ih),GAMBRN(18,3,ih)
       WRITE(*,50)IH,GAMBRN(50,1,ih),GAMBRN(50,2,ih),GAMBRN(50,3,ih)
       WRITE(*,51)IH,GAMBRN(51,1,ih),GAMBRN(51,2,ih),GAMBRN(51,3,ih)
       WRITE(*,52)IH,GAMBRN(52,1,ih),GAMBRN(52,2,ih),GAMBRN(52,3,ih)
       WRITE(*,53)IH,GAMBRN(53,1,ih),GAMBRN(53,2,ih),GAMBRN(53,3,ih)
       WRITE(*,54)IH,GAMBRN(54,1,ih),GAMBRN(54,2,ih),GAMBRN(54,3,ih)
       WRITE(*,55)IH,GAMBRN(55,1,ih),GAMBRN(55,2,ih),GAMBRN(55,3,ih)
       WRITE(*,56)IH,GAMBRN(56,1,ih),GAMBRN(56,2,ih),GAMBRN(56,3,ih)
       WRITE(*,57)IH,GAMBRN(57,1,ih),GAMBRN(57,2,ih),GAMBRN(57,3,ih)
       WRITE(*,58)IH,GAMBRN(58,1,ih),GAMBRN(58,2,ih),GAMBRN(58,3,ih)
       WRITE(*,59)IH,GAMBRN(59,1,ih),GAMBRN(59,2,ih),GAMBRN(59,3,ih)
       WRITE(*,60)IH,GAMBRN(60,1,ih),GAMBRN(60,2,ih),GAMBRN(60,3,ih)
       WRITE(*,61)IH,GAMBRN(61,1,ih),GAMBRN(61,2,ih),GAMBRN(61,3,ih)
       WRITE(*,62)IH,GAMBRN(62,1,ih),GAMBRN(62,2,ih),GAMBRN(62,3,ih)
       WRITE(*,63)IH,GAMBRN(63,1,ih),GAMBRN(63,2,ih),GAMBRN(63,3,ih)
       WRITE(*,64)IH,GAMBRN(64,1,ih),GAMBRN(64,2,ih),GAMBRN(64,3,ih)
       WRITE(*,65)IH,GAMBRN(65,1,ih),GAMBRN(65,2,ih),GAMBRN(65,3,ih)
       WRITE(*,66)IH,GAMBRN(66,1,ih),GAMBRN(66,2,ih),GAMBRN(66,3,ih)
       WRITE(*,67)IH,GAMBRN(67,1,ih),GAMBRN(67,2,ih),GAMBRN(67,3,ih)
       WRITE(*,68)IH,GAMBRN(68,1,ih),GAMBRN(68,2,ih),GAMBRN(68,3,ih)
       WRITE(*,69)IH,GAMBRN(69,1,ih),GAMBRN(69,2,ih),GAMBRN(69,3,ih)
       WRITE(*,70)IH,GAMBRN(70,1,ih),GAMBRN(70,2,ih),GAMBRN(70,3,ih)
       WRITE(*,71)IH,GAMBRN(71,1,ih),GAMBRN(71,2,ih),GAMBRN(71,3,ih)
       WRITE(*,72)IH,GAMBRN(72,1,ih),GAMBRN(72,2,ih),GAMBRN(72,3,ih)
       WRITE(*,73)IH,GAMBRN(73,1,ih),GAMBRN(73,2,ih),GAMBRN(73,3,ih)
       WRITE(*,74)IH,GAMBRN(74,1,ih),GAMBRN(74,2,ih),GAMBRN(74,3,ih)
       WRITE(*,75)IH,GAMBRN(75,1,ih),GAMBRN(75,2,ih),GAMBRN(75,3,ih)
       WRITE(*,76)IH,GAMBRN(76,1,ih),GAMBRN(76,2,ih),GAMBRN(76,3,ih)
       WRITE(*,77)IH,GAMBRN(77,1,ih),GAMBRN(77,2,ih),GAMBRN(77,3,ih)
       WRITE(*,100)IH,GAMBRN(100,1,ih),GAMBRN(100,2,ih),GAMBRN(100,3,ih)
       WRITE(*,101)IH,GAMBRN(101,1,ih),GAMBRN(101,2,ih),GAMBRN(101,3,ih)
      print*,' '
      IF(IH.EQ.1) THEN
      print*,'* Note : WIDTH=GAMBRN(IM,1,1), BR[SM]   =GAMBRN(IM,2,1) '
      print*,'                      and      BR[TOTAL]=GAMBRN(IM,3,1) '
      ELSEIF(IH.EQ.2) THEN
      print*,'* Note : WIDTH=GAMBRN(IM,1,2), BR[SM]   =GAMBRN(IM,2,2) '
      print*,'                      and      BR[TOTAL]=GAMBRN(IM,3,2) '
      ELSEIF(IH.EQ.3) THEN
      print*,'* Note : WIDTH=GAMBRN(IM,1,3), BR[SM]   =GAMBRN(IM,2,3) '
      print*,'                      and      BR[TOTAL]=GAMBRN(IM,3,3) '
      ENDIF
      IF(IPRI.EQ.5) THEN
       IF(IC.EQ.1.OR.IC.EQ.2) GOTO 999
      ENDIF
      print*,'---------------------------------------------------------'
*
  1   FORMAT(1X,'H',I1,' -> e    e  [  1]:',3(2X,E10.4))
  2   FORMAT(1X,'H',I1,' -> mu   mu [  2]:',3(2X,E10.4))
  3   FORMAT(1X,'H',I1,' -> tau  tau[  3]:',3(2X,E10.4))
  4   FORMAT(1X,'H',I1,' -> d    d  [  4]:',3(2X,E10.4))
  5   FORMAT(1X,'H',I1,' -> s    s  [  5]:',3(2X,E10.4))
  6   FORMAT(1X,'H',I1,' -> b    b  [  6]:',3(2X,E10.4))
  7   FORMAT(1X,'H',I1,' -> u    u  [  7]:',3(2X,E10.4))
  8   FORMAT(1X,'H',I1,' -> c    c  [  8]:',3(2X,E10.4))
  9   FORMAT(1X,'H',I1,' -> t    t  [  9]:',3(2X,E10.4))
 10   FORMAT(1X,'H',I1,' -> W    W  [ 10]:',3(2X,E10.4))
 11   FORMAT(1X,'H',I1,' -> Z    Z  [ 11]:',3(2X,E10.4))
 12   FORMAT(1X,'H',I1,' -> H1   Z  [ 12]:',3(2X,E10.4))
 13   FORMAT(1X,'H',I1,' -> H2   Z  [ 13]:',3(2X,E10.4))
 14   FORMAT(1X,'H',I1,' -> H1   H1 [ 14]:',3(2X,E10.4))
 15   FORMAT(1X,'H',I1,' -> H1   H2 [ 15]:',3(2X,E10.4))
 16   FORMAT(1X,'H',I1,' -> H2   H2 [ 16]:',3(2X,E10.4))
 17   FORMAT(1X,'H',I1,' -> ph   ph [ 17]:',3(2X,E10.4))
 18   FORMAT(1X,'H',I1,' -> gl   gl [ 18]:',3(2X,E10.4))
 50   FORMAT(1X,'H',I1,' TOTAL(SM)  [ 50]:',3(2X,E10.4))
 51   FORMAT(1X,'H',I1,' -> N1   N1 [ 51]:',3(2X,E10.4))
 52   FORMAT(1X,'H',I1,' -> N1   N2 [ 52]:',3(2X,E10.4))
 53   FORMAT(1X,'H',I1,' -> N1   N3 [ 53]:',3(2X,E10.4))
 54   FORMAT(1X,'H',I1,' -> N1   N4 [ 54]:',3(2X,E10.4))
 55   FORMAT(1X,'H',I1,' -> N2   N2 [ 55]:',3(2X,E10.4))
 56   FORMAT(1X,'H',I1,' -> N2   N3 [ 56]:',3(2X,E10.4))
 57   FORMAT(1X,'H',I1,' -> N2   N4 [ 57]:',3(2X,E10.4))
 58   FORMAT(1X,'H',I1,' -> N3   N3 [ 58]:',3(2X,E10.4))
 59   FORMAT(1X,'H',I1,' -> N3   N4 [ 59]:',3(2X,E10.4))
 60   FORMAT(1X,'H',I1,' -> N4   N4 [ 60]:',3(2X,E10.4))
 61   FORMAT(1X,'H',I1,' -> C1+  C1-[ 61]:',3(2X,E10.4))
 62   FORMAT(1X,'H',I1,' -> C1+  C2-[ 62]:',3(2X,E10.4))
 63   FORMAT(1X,'H',I1,' -> C2+  C1-[ 63]:',3(2X,E10.4))
 64   FORMAT(1X,'H',I1,' -> C2+  C2-[ 64]:',3(2X,E10.4))
 65   FORMAT(1X,'H',I1,' -> ST1* ST1[ 65]:',3(2X,E10.4))
 66   FORMAT(1X,'H',I1,' -> ST1* ST2[ 66]:',3(2X,E10.4))
 67   FORMAT(1X,'H',I1,' -> ST2* ST1[ 67]:',3(2X,E10.4))
 68   FORMAT(1X,'H',I1,' -> ST2* ST2[ 68]:',3(2X,E10.4))
 69   FORMAT(1X,'H',I1,' -> SB1* SB1[ 69]:',3(2X,E10.4))
 70   FORMAT(1X,'H',I1,' -> SB1* SB2[ 70]:',3(2X,E10.4))
 71   FORMAT(1X,'H',I1,' -> SB2* SB1[ 71]:',3(2X,E10.4))
 72   FORMAT(1X,'H',I1,' -> SB2* SB2[ 72]:',3(2X,E10.4))
 73   FORMAT(1X,'H',I1,' ->STA1*STA1[ 73]:',3(2X,E10.4))
 74   FORMAT(1X,'H',I1,' ->STA1*STA2[ 74]:',3(2X,E10.4))
 75   FORMAT(1X,'H',I1,' ->STA2*STA1[ 75]:',3(2X,E10.4))
 76   FORMAT(1X,'H',I1,' ->STA2*STA2[ 76]:',3(2X,E10.4))
 77   FORMAT(1X,'H',I1,' ->SNU3*SNU3[ 77]:',3(2X,E10.4))
100   FORMAT(1X,'H',I1,' TOTAL(SUSY)[100]:',3(2X,E10.4))
101   FORMAT(1X,'H',I1,' TOTAL      [101]:',3(2X,E10.4))
      RETURN
      END

      SUBROUTINE HSS(SF,V,G,MH,MHJ,MHK,GAM)
************************************************************************
*
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      COMPLEX*16 G
*
      PI=2.D0*DASIN(1.D0)
*
      XKJ=MHJ**2/MH**2
      XKK=MHK**2/MH**2
      XLAM=(1.D0-XKJ-XKK)**2-4.D0*XKJ*XKK
      IF(MH.GT.(MHJ+MHK) .AND. XLAM.GT.0.D0) THEN
       GAM=SF*V**2*CDABS(G)**2/16.D0/PI/MH*DSQRT(XLAM)
      ELSE
       GAM=0.D0
      ENDIF
*
      RETURN
      END

      SUBROUTINE HHV(GF,GV,MH,MHJ,MV,GAM)
************************************************************************
*
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      COMPLEX*16 GV
*
      PI=2.D0*DASIN(1.D0)
*
      XKV=MV**2/MH**2
      XKJ=MHJ**2/MH**2
      XLAM=(1.D0-XKJ-XKV)**2-4.D0*XKJ*XKV
      IF(MH.GT.(MHJ+MV) .AND. XLAM.GT.0.D0) THEN
       GAM=GF*MH**3/8.D0/DSQRT(2.D0)/PI*CDABS(GV)**2*XLAM*DSQRT(XLAM)
      ELSE
       GAM=0.D0
      ENDIF
*
      RETURN
      END

      SUBROUTINE HVV(GF,DV,GV,MV,MH,GAM)
************************************************************************
*
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      COMPLEX*16 GV
*
      PI=2.D0*DASIN(1.D0)
*
      XKV=MV**2/MH**2
      BETASQ=1.D0-4.D0*XKV
      IF(BETASQ.GE.0.D0) THEN
       BETA=DSQRT(BETASQ)
       GAM=GF*CDABS(GV)**2*MH**3*DV/16.D0/DSQRT(2.D0)/PI*BETA
     .    *(1.D0-4.D0*XKV+12.D0*XKV**2)
      ELSE
       GAM=0.D0
      ENDIF
*
      RETURN
      END

      SUBROUTINE HFF(CF,SF,DJK,GF,GS,GP,MH,MJ,MK,GAM)
************************************************************************
*
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      COMPLEX*16 GF,GS,GP
*
      PI=2.D0*DASIN(1.D0)
*
      XKJ=MJ**2/MH**2
      XKK=MK**2/MH**2
      XLAM=(1.D0-XKJ-XKK)**2-4.D0*XKJ*XKK
      IF(XLAM.LE.0.D0 .OR. MH.LE.(MJ+MK)) THEN
       GAM=0.D0
       RETURN
      ELSE
       GAM=SF**2/(1.D0+DJK)*CF*CDABS(GF)**2*MH*DSQRT(XLAM)*
     .    ((1.D0-XKJ-XKK)*(CDABS(GS)**2+CDABS(GP)**2)-
     .     2.D0*DSQRT(XKJ*XKK)*(CDABS(GS)**2-CDABS(GP)**2))/8.D0/PI
*
*       IF(MJ.EQ.MK) THEN
*         BETA=DSQRT(1.D0-4.D0*XKJ)
*         GAMP=SF**2/(1.D0+DJK)*CF*CDABS(GF)**2*MH*
*     .        BETA*(BETA**2*CDABS(GS)**2+CDABS(GP)**2)/8.D0/PI
*         print*,'same ?',gam,gamp
*       ENDIF
*
      ENDIF
*
      RETURN
      END

      SUBROUTINE BODE(FBODE,XIN,XOUT,NSTEP,YINT)
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%                                                                      %
c%     BODE makes one-dimensional numerical integration of a            % 
c%     F(x), using Bode's rule:                                         %
c%                                                                      %
c%     Integral_x0^x2 F(x) dx  =   h/598752 * [ 16067 (f_0 + f_10)      %
c%     + 106300 (f_1 + f_9) - 48525 (f_2 + f_8) + 272400 (f_3 +f_7)     %
c%     - 260550 (f_4 + f_6) + 427368 f_5 ]                              %
c%                                                                      %
c%    FBODE defines the integrand  F(x)                                 %
c%    XIN  is the lower limit of the integral                           %
c%    XOUT is the upper limit of the integral                           %
c%    NSTEP determines the total number of steps                        %
c%    YINT contains the result of the integration                       %
c%                                                                      %
c%                                                                      %
c%    NOTES: 1. We always assume that XOUT > XIN                        %
c%           2. The integrand F(x) must be defined as an external       %
C%              function in the main program                            %
c%                                                                      %
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      IMPLICIT REAL*8(A,B,D-H,O-Z), COMPLEX*16(C)
C
CCC
C
      YINT = 0.D0
*      DELTA = (XOUT - XIN)/DREAL(NSTEP)
      DELTA = (XOUT - XIN)/DBLE(NSTEP)
      X10 = XIN
999   CONTINUE
      X0 = X10
       X1 = X10 + DELTA/10.D0
        X2 = X10 + DELTA/5.D0
         X3 = X10 + 3.D0*DELTA/10.D0
          X4 = X10 + 2.D0*DELTA/5.D0
           X5 = X10 + DELTA/2.D0
            X6 = X10 + 3.D0*DELTA/5.D0
             X7 = X10 + 7.D0*DELTA/10.D0
              X8 = X10 + 4.D0*DELTA/5.D0
               X9 = X10 + 9.D0*DELTA/10.D0
                X10 = X10 + DELTA
      IF(X10.GE.XOUT) GOTO 9999
      F0 = FBODE(X0)
       F1 = FBODE(X1)
        F2 = FBODE(X2)
         F3 = FBODE(X3)
          F4 = FBODE(X4)
           F5 = FBODE(X5)
            F6 = FBODE(X6)
             F7 = FBODE(X7)
              F8 = FBODE(X8)
               F9 = FBODE(X9)
                F10 = FBODE(X10)
      YINT = DELTA*( 16067.D0*(F0+F10) + 106300.D0*(F1+F9)
     #- 48525.D0*(F2+F8) + 272400.D0*(F3+F7) - 260550.D0*(F4+F6)
     #+ 427368.D0*F5  ) / 598752.D0  + YINT
      GOTO 999
9999  CONTINUE
      RETURN
      END
C
      REAL*8 FUNCTION FVVS(R)
************************************************************************
*
* Here, we used the integration method:
*
* I    =               \int_{x-}^{x+} dx f(x) 
*      = {G(x+)-G(x-)} \int_{0}^{1}   dR f(x)/g(x)
*
*  where      f(x) = F(x)/[(x-a)^2+b^2] 
*             g(x) = b^2/[(x-a)^2+b^2]
*             G(x) = pi/2*b+b*arctg[(x-a)/b]
*               x  = G^{-1}[G(x-)+R*{G(x+)-G(x-)}]
*        G^{-1}(y) = a+b*tg[(y-pi/2*b)/b]
*
*  Note : a=1 and b=EPSV 
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
      COMMON /HC_BODE/ EPSV,OMEGAI,OMEGAJ,XUP,XDW
*
      PI=2.D0*DASIN(1.D0)
      GXUP=PI/2.D0*EPSV+EPSV*DATAN((XUP-1.D0)/EPSV)
      GXDW=PI/2.D0*EPSV+EPSV*DATAN((XDW-1.D0)/EPSV)
      Y=GXDW+R*(GXUP-GXDW)
      XX=1.D0+EPSV*DTAN((Y-PI/2.D0*EPSV)/EPSV)
      XLAM=(1.D0-OMEGAI-XX)**2-4.D0*OMEGAI*XX
      IF(XLAM.LT.0.D0) XLAM=0.D0
      F_XX=DSQRT(XLAM)*(XLAM+12.D0*XX)
      FVVS=(GXUP-GXDW)*F_XX/EPSV**2
*
      RETURN
      END

      REAL*8 FUNCTION FHVS(R)
************************************************************************
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
      COMMON /HC_BODE/ EPSV,OMEGAI,OMEGAJ,XUP,XDW
*
      PI=2.D0*DASIN(1.D0)
      GXUP=PI/2.D0*EPSV+EPSV*DATAN((XUP-1.D0)/EPSV)
      GXDW=PI/2.D0*EPSV+EPSV*DATAN((XDW-1.D0)/EPSV)
      Y=GXDW+R*(GXUP-GXDW)
      XX=1.D0+EPSV*DTAN((Y-PI/2.D0*EPSV)/EPSV)
      XLAM=(XX-OMEGAI-OMEGAJ)**2-4.D0*OMEGAI*OMEGAJ
      IF(XLAM.LT.0.D0) XLAM=0.D0
      F_XX=DSQRT(XLAM)**3
      FHVS=(GXUP-GXDW)*F_XX/EPSV**2
*
      RETURN
      END

      SUBROUTINE CHWBB(NX1,GW,GTB,GS,GP,MCH,TANB,MT,MW,MB,GTOP,GCHWBB)
************************************************************************
*
* Calculate the three body decay rate H+ -> t* b^bar -> W^+ b b^bar
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*-----------------------------------------------------------------------
* For integration
      COMMON /CHWBB_RBODE/ MCH_CHWBB,TANB_CHWBB,MT_CHWBB,MW_CHWBB
     .                    ,MB_CHWBB,GTOP_CHWBB
      COMPLEX*16 GS_CHWBB,GP_CHWBB
      COMMON /CHWBB_CBODE/ GS_CHWBB,GP_CHWBB
      EXTERNAL FX1
*Local
      COMPLEX*16 GTB,GS,GP
      PI=2.D0*DASIN(1.D0)
*-----------------------------------------------------------------------
      MCH_CHWBB  = MCH
      TANB_CHWBB = TANB
      MT_CHWBB   = MT
      MW_CHWBB   = MW
      MB_CHWBB   = MB
      GTOP_CHWBB = GTOP
*
      GS_CHWBB   = GS
      GP_CHWBB   = GP
*
      XKW = MW**2/MCH**2
*
      X1D=0.D0
      X1U=1.D0-MW**2/MCH**2
*
      CALL BODE(FX1,X1D,X1U,NX1,RES)
*
      GCHWBB=3.D0*GW**2*CDABS(GTB)**2*MCH/512.D0/PI**3*RES
*
      RETURN
      END

      REAL*8 FUNCTION FX1(X1)
************************************************************************
*
* It returns a function FX1 as function of x1 after x2 integration.
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*-----------------------------------------------------------------------
      COMMON /CHWBB_RBODE/ MCH_CHWBB,TANB_CHWBB,MT_CHWBB,MW_CHWBB
     .                    ,MB_CHWBB,GTOP_CHWBB
      COMPLEX*16 GS_CHWBB,GP_CHWBB
      COMMON /CHWBB_CBODE/ GS_CHWBB,GP_CHWBB
*Local
      COMPLEX*16 GS,GP,GL,GR
      COMPLEX*16 XI
*-----------------------------------------------------------------------
      XI   =DCMPLX(0.D0,1.D0)
*
      MCH  = MCH_CHWBB
      TANB = TANB_CHWBB
      MT   = MT_CHWBB
      MW   = MW_CHWBB
      MB   = MB_CHWBB
      GTOP = GTOP_CHWBB
      GS   = GS_CHWBB
      GP   = GP_CHWBB
* 
      XKT = MT**2/MCH**2
      XKW = MW**2/MCH**2
      XKB = MB**2/MCH**2
      XGT = GTOP**2/MCH**2
      GR0 = MB/MT*TANB
      GL0 = 1.D0/TANB
      GL  = GS-XI*GP
      GR  = GS+XI*GP
*      print*,'GR',gr0,dreal(gr),dimag(gr)
*      print*,'GL',gl0,dreal(gl),dimag(gl)
*
*The x2 intgration has been done by REDUCE:
*
      CALL GET_FN(XKT,XKW,XKB,XGT,X1,F0,F1,F2,F3)
*
      FXL2=-2.D0*XKB*XKT*F0+XKT*((1.D0-X1)*(F0-F1)/XKW+2.D0*X1*F0
     .     +2.D0*F1-3.D0*F0+2.D0*XKW*F0)
      FXR2=-2.D0*XKB**2*F0
     .     +XKB*((2.D0*XKW-2.D0*X1+3.D0)*F0
     .     +(-2.D0*F2-X1*F1+X1*F0+5.D0*F1-3.D0*F0)/XKW)
     .     +(F2+2.D0*X1*F1-2.D0*X1*F0-4.D0*F1+3.D0*F0-2.D0*XKW*F0)
     .     +(F3+X1*F2-3.D0*F2-2.D0*X1*F1+X1*F0+3.D0*F1-F0)/XKW
      FXLR=2.D0*XKB*F0-2.D0*XKW*F0-F1+F0+(F2-2.D0*F1+F0)/XKW
  
      FX1=CDABS(GL)**2*FXL2+CDABS(GR)**2*FXR2
     .   +2.D0*DSQRT(XKB*XKT)*DREAL(GL*DCONJG(GR))*FXLR
*
      RETURN
      END


      SUBROUTINE GET_FN(XKT,XKW,XKB,XGT,X1,F0,F1,F2,F3)
************************************************************************
*
*      /               x2^n
* Fn = | dx2  --------------------------
*      /      [(1-x2-xkt+xkb)^2+xkt*xgt]
*
*     x2^up =1-kw/(1-x1)   ; x2_down = 1-kw-x1
*
************************************************************************
      IMPLICIT REAL*8(A-H,M,O-Z)
*
      X2U=1.-XKW/(1.-X1)
      X2D=1.-XKW-X1
*-----------------------------------------------------------------------
*
* This is an out put of a REDUCE program 
*
*-----------------------------------------------------------------------
* f0
      f0=(sqrt(xkt)*sqrt(xgt)*(-atan((x2d-xkb+xkt-1.0)/(sqrt(xkt)*
     . sqrt(xgt)))+atan((x2u-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))))/(
     . xgt*xkt)
* f1
      f1=(-2.0*sqrt(xkt)*sqrt(xgt)*atan((x2d-xkb+xkt-1.0)/(sqrt(xkt)*
     . sqrt(xgt)))*xkb+2.0*sqrt(xkt)*sqrt(xgt)*atan((x2d-xkb+xkt-1.0)
     . /(sqrt(xkt)*sqrt(xgt)))*xkt-2.0*sqrt(xkt)*sqrt(xgt)*atan((x2d-
     . xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))+2.0*sqrt(xkt)*sqrt(xgt)*
     . atan((x2u-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xkb-2.0*sqrt(xkt
     . )*sqrt(xgt)*atan((x2u-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xkt+
     . 2.0*sqrt(xkt)*sqrt(xgt)*atan((x2u-xkb+xkt-1.0)/(sqrt(xkt)*sqrt
     . (xgt)))+log((x2u**2-2.0*x2u*xkb+2.0*x2u*xkt-2.0*x2u+xgt*xkt+
     . xkb**2-2.0*xkb*xkt+2.0*xkb+xkt**2-2.0*xkt+1.0)/(x2d**2-2.0*x2d
     . *xkb+2.0*x2d*xkt-2.0*x2d+xgt*xkt+xkb**2-2.0*xkb*xkt+2.0*xkb+
     . xkt**2-2.0*xkt+1.0))*xgt*xkt)/(2.0*xgt*xkt)
* f2
      ans2=log((x2d**2-2.0*x2d*xkb+2.0*x2d*xkt-2.0*x2d+xgt*xkt+xkb**2
     . -2.0*xkb*xkt+2.0*xkb+xkt**2-2.0*xkt+1.0)/(x2u**2-2.0*x2u*xkb+
     . 2.0*x2u*xkt-2.0*x2u+xgt*xkt+xkb**2-2.0*xkb*xkt+2.0*xkb+xkt**2-
     . 2.0*xkt+1.0))*xgt*xkt**2+log((x2u**2-2.0*x2u*xkb+2.0*x2u*xkt-
     . 2.0*x2u+xgt*xkt+xkb**2-2.0*xkb*xkt+2.0*xkb+xkt**2-2.0*xkt+1.0)
     . /(x2d**2-2.0*x2d*xkb+2.0*x2d*xkt-2.0*x2d+xgt*xkt+xkb**2-2.0*
     . xkb*xkt+2.0*xkb+xkt**2-2.0*xkt+1.0))*xgt*xkb*xkt+log((x2u**2-
     . 2.0*x2u*xkb+2.0*x2u*xkt-2.0*x2u+xgt*xkt+xkb**2-2.0*xkb*xkt+2.0
     . *xkb+xkt**2-2.0*xkt+1.0)/(x2d**2-2.0*x2d*xkb+2.0*x2d*xkt-2.0*
     . x2d+xgt*xkt+xkb**2-2.0*xkb*xkt+2.0*xkb+xkt**2-2.0*xkt+1.0))*
     . xgt*xkt-x2d*xgt*xkt+x2u*xgt*xkt
      ans1=sqrt(xkt)*sqrt(xgt)*atan((x2d-xkb+xkt-1.0)/(sqrt(xkt)*sqrt
     . (xgt)))*xgt*xkt-sqrt(xkt)*sqrt(xgt)*atan((x2d-xkb+xkt-1.0)/(
     . sqrt(xkt)*sqrt(xgt)))*xkb**2+2.0*sqrt(xkt)*sqrt(xgt)*atan((x2d
     . -xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xkb*xkt-2.0*sqrt(xkt)*
     . sqrt(xgt)*atan((x2d-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xkb-
     . sqrt(xkt)*sqrt(xgt)*atan((x2d-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt
     . )))*xkt**2+2.0*sqrt(xkt)*sqrt(xgt)*atan((x2d-xkb+xkt-1.0)/(
     . sqrt(xkt)*sqrt(xgt)))*xkt-sqrt(xkt)*sqrt(xgt)*atan((x2d-xkb+
     . xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))-sqrt(xkt)*sqrt(xgt)*atan((x2u-
     . xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xgt*xkt+sqrt(xkt)*sqrt(xgt
     . )*atan((x2u-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xkb**2-2.0*
     . sqrt(xkt)*sqrt(xgt)*atan((x2u-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt
     . )))*xkb*xkt+2.0*sqrt(xkt)*sqrt(xgt)*atan((x2u-xkb+xkt-1.0)/(
     . sqrt(xkt)*sqrt(xgt)))*xkb+sqrt(xkt)*sqrt(xgt)*atan((x2u-xkb+
     . xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xkt**2-2.0*sqrt(xkt)*sqrt(xgt)
     . *atan((x2u-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xkt+sqrt(xkt)*
     . sqrt(xgt)*atan((x2u-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))+ans2
      f2=ans1/(xgt*xkt)
* f3
      ans4=6.0*log((x2u**2-2.0*x2u*xkb+2.0*x2u*xkt-2.0*x2u+xgt*xkt+
     . xkb**2-2.0*xkb*xkt+2.0*xkb+xkt**2-2.0*xkt+1.0)/(x2d**2-2.0*x2d
     . *xkb+2.0*x2d*xkt-2.0*x2d+xgt*xkt+xkb**2-2.0*xkb*xkt+2.0*xkb+
     . xkt**2-2.0*xkt+1.0))*xgt*xkb*xkt+3.0*log((x2u**2-2.0*x2u*xkb+
     . 2.0*x2u*xkt-2.0*x2u+xgt*xkt+xkb**2-2.0*xkb*xkt+2.0*xkb+xkt**2-
     . 2.0*xkt+1.0)/(x2d**2-2.0*x2d*xkb+2.0*x2d*xkt-2.0*x2d+xgt*xkt+
     . xkb**2-2.0*xkb*xkt+2.0*xkb+xkt**2-2.0*xkt+1.0))*xgt*xkt**3+3.0
     . *log((x2u**2-2.0*x2u*xkb+2.0*x2u*xkt-2.0*x2u+xgt*xkt+xkb**2-
     . 2.0*xkb*xkt+2.0*xkb+xkt**2-2.0*xkt+1.0)/(x2d**2-2.0*x2d*xkb+
     . 2.0*x2d*xkt-2.0*x2d+xgt*xkt+xkb**2-2.0*xkb*xkt+2.0*xkb+xkt**2-
     . 2.0*xkt+1.0))*xgt*xkt-x2d**2*xgt*xkt-4.0*x2d*xgt*xkb*xkt+4.0*
     . x2d*xgt*xkt**2-4.0*x2d*xgt*xkt+x2u**2*xgt*xkt+4.0*x2u*xgt*xkb*
     . xkt-4.0*x2u*xgt*xkt**2+4.0*x2u*xgt*xkt
      ans3=log((x2d**2-2.0*x2d*xkb+2.0*x2d*xkt-2.0*x2d+xgt*xkt+xkb**2
     . -2.0*xkb*xkt+2.0*xkb+xkt**2-2.0*xkt+1.0)/(x2u**2-2.0*x2u*xkb+
     . 2.0*x2u*xkt-2.0*x2u+xgt*xkt+xkb**2-2.0*xkb*xkt+2.0*xkb+xkt**2-
     . 2.0*xkt+1.0))*xgt**2*xkt**2+6.0*log((x2d**2-2.0*x2d*xkb+2.0*
     . x2d*xkt-2.0*x2d+xgt*xkt+xkb**2-2.0*xkb*xkt+2.0*xkb+xkt**2-2.0*
     . xkt+1.0)/(x2u**2-2.0*x2u*xkb+2.0*x2u*xkt-2.0*x2u+xgt*xkt+xkb**
     . 2-2.0*xkb*xkt+2.0*xkb+xkt**2-2.0*xkt+1.0))*xgt*xkb*xkt**2+6.0*
     . log((x2d**2-2.0*x2d*xkb+2.0*x2d*xkt-2.0*x2d+xgt*xkt+xkb**2-2.0
     . *xkb*xkt+2.0*xkb+xkt**2-2.0*xkt+1.0)/(x2u**2-2.0*x2u*xkb+2.0*
     . x2u*xkt-2.0*x2u+xgt*xkt+xkb**2-2.0*xkb*xkt+2.0*xkb+xkt**2-2.0*
     . xkt+1.0))*xgt*xkt**2+3.0*log((x2u**2-2.0*x2u*xkb+2.0*x2u*xkt-
     . 2.0*x2u+xgt*xkt+xkb**2-2.0*xkb*xkt+2.0*xkb+xkt**2-2.0*xkt+1.0)
     . /(x2d**2-2.0*x2d*xkb+2.0*x2d*xkt-2.0*x2d+xgt*xkt+xkb**2-2.0*
     . xkb*xkt+2.0*xkb+xkt**2-2.0*xkt+1.0))*xgt*xkb**2*xkt+ans4
      ans2=-6.0*sqrt(xkt)*sqrt(xgt)*atan((x2u-xkb+xkt-1.0)/(sqrt(xkt)
     . *sqrt(xgt)))*xgt*xkb*xkt+6.0*sqrt(xkt)*sqrt(xgt)*atan((x2u-xkb
     . +xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xgt*xkt**2-6.0*sqrt(xkt)*sqrt
     . (xgt)*atan((x2u-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xgt*xkt+
     . 2.0*sqrt(xkt)*sqrt(xgt)*atan((x2u-xkb+xkt-1.0)/(sqrt(xkt)*sqrt
     . (xgt)))*xkb**3-6.0*sqrt(xkt)*sqrt(xgt)*atan((x2u-xkb+xkt-1.0)/
     . (sqrt(xkt)*sqrt(xgt)))*xkb**2*xkt+6.0*sqrt(xkt)*sqrt(xgt)*atan
     . ((x2u-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xkb**2+6.0*sqrt(xkt)
     . *sqrt(xgt)*atan((x2u-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xkb*
     . xkt**2-12.0*sqrt(xkt)*sqrt(xgt)*atan((x2u-xkb+xkt-1.0)/(sqrt(
     . xkt)*sqrt(xgt)))*xkb*xkt+6.0*sqrt(xkt)*sqrt(xgt)*atan((x2u-xkb
     . +xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xkb-2.0*sqrt(xkt)*sqrt(xgt)*
     . atan((x2u-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xkt**3+6.0*sqrt(
     . xkt)*sqrt(xgt)*atan((x2u-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*
     . xkt**2-6.0*sqrt(xkt)*sqrt(xgt)*atan((x2u-xkb+xkt-1.0)/(sqrt(
     . xkt)*sqrt(xgt)))*xkt+2.0*sqrt(xkt)*sqrt(xgt)*atan((x2u-xkb+xkt
     . -1.0)/(sqrt(xkt)*sqrt(xgt)))+ans3
      ans1=6.0*sqrt(xkt)*sqrt(xgt)*atan((x2d-xkb+xkt-1.0)/(sqrt(xkt)*
     . sqrt(xgt)))*xgt*xkb*xkt-6.0*sqrt(xkt)*sqrt(xgt)*atan((x2d-xkb+
     . xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xgt*xkt**2+6.0*sqrt(xkt)*sqrt(
     . xgt)*atan((x2d-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xgt*xkt-2.0
     . *sqrt(xkt)*sqrt(xgt)*atan((x2d-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(
     . xgt)))*xkb**3+6.0*sqrt(xkt)*sqrt(xgt)*atan((x2d-xkb+xkt-1.0)/(
     . sqrt(xkt)*sqrt(xgt)))*xkb**2*xkt-6.0*sqrt(xkt)*sqrt(xgt)*atan(
     . (x2d-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xkb**2-6.0*sqrt(xkt)*
     . sqrt(xgt)*atan((x2d-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xkb*
     . xkt**2+12.0*sqrt(xkt)*sqrt(xgt)*atan((x2d-xkb+xkt-1.0)/(sqrt(
     . xkt)*sqrt(xgt)))*xkb*xkt-6.0*sqrt(xkt)*sqrt(xgt)*atan((x2d-xkb
     . +xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xkb+2.0*sqrt(xkt)*sqrt(xgt)*
     . atan((x2d-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*xkt**3-6.0*sqrt(
     . xkt)*sqrt(xgt)*atan((x2d-xkb+xkt-1.0)/(sqrt(xkt)*sqrt(xgt)))*
     . xkt**2+6.0*sqrt(xkt)*sqrt(xgt)*atan((x2d-xkb+xkt-1.0)/(sqrt(
     . xkt)*sqrt(xgt)))*xkt-2.0*sqrt(xkt)*sqrt(xgt)*atan((x2d-xkb+xkt
     . -1.0)/(sqrt(xkt)*sqrt(xgt)))+ans2
      f3=ans1/(2.0*xgt*xkt)
*-----------------------------------------------------------------------
*END 
*-----------------------------------------------------------------------
      RETURN
      END
