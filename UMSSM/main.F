
C====== Modules =================
C   Keys to switch on
C   various modules of micrOMEGAs
C================================
   
#define MASSES_INFO
C    Display information about SUSY and Higgs masses
#define CONSTRAINTS
C    Display  deltarho, B and K observables, gmuon, Z1->inv,
C    check LEP mass limits and Zprime limits
C#define SLHAINPUT
C    Switch to use slha files as input instead of umssm.par or similar files
C#define HIGGSBOUNDS 
#define LILITH
C#define SMODELS

#define OMEGA
C    Calculate relic density and display contribution of
C                            individual channels 
#define INDIRECT_DETECTION
C     Compute spectra of gamma/positron/neutrinos
C     for DM annihilation; calculate <sigma*v> and
C     integrate gamma signal over DM galactic squared
C     density for given line of sight.
C#define RESET_FORMFACTORS
C     Modify default nucleus form factors, 
C     DM velocity distribution,
C     A-dependence of Fermi-dencity.
#define CDM_NUCLEON 
C     Calculate amplitudes and cross-sections for
C     CDM-mucleon collisions 
C#define CDM_NUCLEUS
C     Calculate number of events for 1kg*day
C     and recoil energy distibution for various nuclei

C================ end of Modules  ======

#define NEUTRINO   neutrino telescope

C#define DECAYS 
C     Calculates particle widths and branchings

#define CROSS_SECTIONS
C     Calculate cross sections and widths for
C     reactions specified by the user

C===== Options ========

C#define SHOWPLOTS
C     Display  graphical plots on the screen

#define CLEAN to crean intermediate files

C===== End of DEFINE  settings ===== 




      IMPLICIT NONE
C! external functions 
#include "../include/micromegas.fh"
#include "lib/pmodel.fh"
*      external iargc
*      integer iargc
C! local  Parameters   

      INTEGER fast/1/
      REAL*8  Beps/1.E-4/,cut/0.01/
C Variables      
      INTEGER err
      REAL*8 Omega,OmegaAppr,Xf
      character *100 argv, mess, buff
      character*20 cdmName
      integer spin2,charge3,cdim,PDGcdm
      integer narg 
      real*8 Emin,sigmaV,fi,dfi,Etest,cs,Qren,Qfact,pTmin
      real*8 Pcm,SMmev 
      real*8 SpA(250),Spe(250),SpP(250) 
      real*8 FluxA(250),FluxE(250),FluxP(250)
C  for CONSTRAINTS    
      real*8 constr0,constrM, constrP,csLim
C  for monoJet
      real*8 monojetCL
C for HIGGSBOUNDS/HIGGSSIGNALS
      integer  NH0/4/, NHch/1/
      real*8 HB_result,HB_obsratio,HS_observ,HS_chi2, HS_pval
      character*70 HB_chan
      character*7 HB_version, HS_version
C for SMODELS 
      integer result 
      real*8 Rvalue
      character*30 analysis,topology,smodels_version
C for LILITH 
      real*8 m2logL, m2logL_reference,pvalue,vPython1, vPython2
      integer exp_ndf,n_par,ndf
      character*100 call_lilith
      character*10 lilith_version
C for CDM_NUCLEON
      real*8  pA0(2),pA5(2),nA0(2),nA5(2),Nmass,SCcoeff
C  for CDM_NUCLEUS
      real*8  dNdE(200),rho,nEvents
C for DECAYS
      integer txtL(2)
      real*8 width   
      character*20 pname

C  for NEUTRINO 
      integer forSun
      real*8 IntFlux, rhoDet
      real*8 FluxN(NZ),FluxNb(NZ),FluxMu(NZ) 

C   for CROSS_SECTION
       character*10 next
       character*10 next_
       real*8 nextM 


C   to Force calculation of matrix elements in Unitary Gauge call with (1)
      call ForceUG(0)

C   to save the SLHA  input/output file  uncomment  
C*      call delFiles(0)


#ifdef SLHAINPUT
      write(*,*) '\n========= SLHA file input ========='
      narg=iargc()
      if(narg.eq.0) then 
        write(*,*)'The program needs one argument: name of SLHA file'
        write(*,*) 'Example: ./main UMSSM_inp.dat'
        stop
      endif
      call getarg(1,argv)
      write(*,*) 'Initial file: ',argv
     
      err=lesHinput(argv)
      if(err.ne.0) stop
#else
 
      if(iargc().eq.0) then
         write(*,*)'Correct usage: ./main  <file with parameters>'
         write(*,*)'Example: ./main  umssm.par' 
         stop
      endif

      call getarg(1,argv)
       
      err=readVar(argv)

#endif

      if(err.eq.-1) then
          write(*,*) 'Can not open the file'
          stop
      endif

      if(err.gt.0) then
          write(*,*)'Wrong file contents at line ',err
          stop
      endif

C  To print list of variables uncomment   
C      CALL printVar(6)

      err=sortOddParticles(cdmName)

      if(err.ne.0) then
        write(*,*) 'Can t calculate ', cdmName
        STOP
      endif

      PDGcdm=qNumbers(cdmName,spin2, charge3,  cdim)
      write(*,fmt="('Dark matter candidate is ', A5,'with spin=',
     >I1,'/2 mass=',1PE10.2)") cdmName, spin2,Mcdm

      if(charge3.ne.0) then 
         write(*,fmt="('DM has electric charge',I1,'/3')") charge3
         STOP
      endif

      if(cdim.ne.1) then 
         write(*,*) 'DM has color charge'
         STOP
      endif

      if(cdmName.ne.'~o1') then
          write(*,*) '~o1 is not LSP'
      else
          call o1Contents(6)
      endif

      err=umssmtools(PDGcdm)

      if(err.ne.0) then
        write(*,*) 'An error occurred running umssmtools.'
        STOP
      endif
      
C Get the corrected Higgs branching ratios from UMSSMTools :
      err=slhaRead('UMSSM_decay.dat',1)
      if(err.ne.0) stop

#ifdef MASSES_INFO
      write(*,*) '=== MASSES OF HIGGS AND SUSY PARTICLES: ===' 
      CALL printHiggs(6)
      CALL printMasses(6,1)
#endif


#ifdef CONSTRAINTS
      err=slhaRead('UMSSM_spectr.dat',0)
      if(err.ne.0) stop
      write(*,*)
      write(*,*) '==== Physical Constraints: ====='
      write(*,fmt='(A,1x, 1PE8.2)') 'deltartho=', deltarho()
      constr0=bsg(constrM,constrP)
      write(*,fmt="('b  -> s gamma     = ',1PE9.2,'(',E9.2,',',E9.2,')')") 
     >constr0,constrM, constrP

      constr0=bsmumu(constrM,constrP)
      write(*,fmt="('Bs -> mu+ mu-     = ',1PE9.2,'(',E9.2,',',E9.2,')')") 
     >constr0,constrM, constrP 
  
      constr0=btaunu(constrM,constrP)
      write(*,fmt="('B+ -> tau+ nu_tau = ',1PE9.2,'(',E9.2,',',E9.2,')')") 
     >constr0,constrM, constrP
  
      constr0=deltamd(constrM,constrP)
      write(*,fmt="('delta M_d         = ',1PE9.2,'(',E9.2,',',E9.2,') ps^-1')") 
     >constr0,constrM, constrP

      constr0=deltams(constrM,constrP)
      write(*,fmt="('delta M_s         = ',1PE9.2,'(',E9.2,',',E9.2,') ps^-1')") 
     >constr0,constrM, constrP

      constr0=gmuon(constrM,constrP)
      write(*,fmt="('gmuon/BSM         = ',1PE9.2,'(',E9.2,',',E9.2,')')") 
     >constr0,constrM, constrP

      constr0=bxislllow(constrM,constrP)
      write(*,fmt="('B  -> X_s l+ l- for low M_{l+l-}^2  = ',1PE9.2,'(',E9.2,',',E9.2,')')") 
     >constr0,constrM, constrP

      constr0=bxisllhigh(constrM,constrP)
      write(*,fmt="('B  -> X_s l+ l- for high M_{l+l-}^2 = ',1PE9.2,'(',E9.2,',',E9.2,')')") 
     >constr0,constrM, constrP

      constr0=bdg(constrM,constrP)
      write(*,fmt="('b  -> d gamma     = ',1PE9.2,'(',E9.2,',',E9.2,')')") 
     >constr0,constrM, constrP

      constr0=bdmumu(constrM,constrP)
      write(*,fmt="('Bd -> mu+ mu-     = ',1PE9.2,'(',E9.2,',',E9.2,')')") 
     >constr0,constrM, constrP

      constr0=bxisnunu(constrM,constrP)
      write(*,fmt="('B  -> Xs nu_L nubar_L = ',1PE9.2,'(',E9.2,',',E9.2,')')") 
     >constr0,constrM, constrP

      constr0=bpkpnunu(constrM,constrP)
      write(*,fmt="('B+ -> K+ nu_L nubar_L = ',1PE9.2,'(',E9.2,',',E9.2,')')") 
     >constr0,constrM, constrP

      constr0=bksnunu(constrM,constrP)
      write(*,fmt="('B  -> Ks nu_L nubar_L = ',1PE9.2,'(',E9.2,',',E9.2,')')") 
     >constr0,constrM, constrP

      constr0=rdtaul(constrM,constrP)
      write(*,fmt="('RD  = BR[B+ -> D  tau+ nu_tau]/BR[B+ -> D  l+ nu_l] = ',1PE9.2,'(',E9.2,',',E9.2,')')") 
     >constr0,constrM, constrP

      constr0=rdstaul(constrM,constrP)
      write(*,fmt="('RD* = BR[B+ -> D* tau+ nu_tau]/BR[B+ -> D* l+ nu_l] = ',1PE9.2,'(',E9.2,',',E9.2,')')") 
     >constr0,constrM, constrP

      constr0=kppipnunu(constrM,constrP)
      write(*,fmt="('K+ -> Pi+ nu_L nubar_L = ',1PE9.2,'(',E9.2,',',E9.2,')')") 
     >constr0,constrM, constrP

      constr0=klpi0nunu(constrM,constrP)
      write(*,fmt="('KL -> Pi0 nu_L nubar_L = ',1PE9.2,'(',E9.2,',',E9.2,')')") 
     >constr0,constrM, constrP

      constr0=deltamk(constrM,constrP)
      write(*,fmt="('delta M_K              = ',1PE9.2,'(',E9.2,',',E9.2,') ps^-1')") 
     >constr0,constrM, constrP

      constr0=epsk(constrM,constrP)
      write(*,fmt="('eps_K                  = ',1PE9.2,'(',E9.2,',',E9.2,')')") 
     >constr0,constrM, constrP


      if(MassLimits().eq.0)  write(*,*) 'LEP limits OK'
      if(zprimelimits().eq.0) write(*,*) 'LHC limits on new Zprime OK'
      if(Zinvisible().ne.0) write(*,*) "Excluded by Z->invisible"
      if(LspNlsp_LEP(csLim).ne.0) write(*,fmt='("Excluded by LEP  by e+,e- -> DM q qbar, Cross section=",E10.2,"pb")') csLim 
#endif

 
#if defined(HIGGSBOUNDS) || defined(HIGGSSIGNALS)
      NH0=hbBblocksMDL("HB.in",NHch)
      call system('echo BLOCK DMASS  >> HB.in')
      call system('echo "  25  2"  >> HB.in')
#include "../include/hBandS.inc_f"
#ifdef HIGGSBOUNDS
      write(*,fmt='("HB(",A, ")result=",1PE10.1,3x,"obsratio=", E10.1," channel=", A )') HB_version, HB_result,HB_obsratio,HB_chan
#endif
#ifdef HIGGSSIGNALS
      write(*,fmt='("HS(",A, "): Nobservables=",F5.0," chi^2 =", 1PE10.2," pval= ", E10.2)') HS_version,HS_observ,HS_chi2, HS_pval
#endif
#endif


#ifdef LILITH
      if(LilithMDL("Lilith_in.xml").gt.0) then 
         m2logL_reference=0
         n_par=0
#include "../include/Lilith.inc_f"
         if(ndf.ne.0)write(*,fmt='("LILITH(DB",A,"):  -2*log(L): ",F5.2,"; -2*log(L_reference): ",F5.2,"; ndf: ",I3"; p-value: ",F5.2)') lilith_version, m2logL,m2logL_reference,ndf,pvalue
      else 
         write(*,*) "LILITH: there is no Higgs candidate"
      endif     
#endif


#ifdef SMODELS
      write(*,*) "SMODELS:"
#include "../include/SMODELS.inc_f"
      if(result.ge.0) write(*,fmt='("highest R=",1PE10.2," from ",A," topology ",A)') Rvalue,analysis,topology
      if(result.eq.-1) write(*,*) 'no channels for testing.'
      if(result.eq.0)  write(*,*) 'not excluded.' 
      if(result.eq.1)  write(*,*) 'excluded.'
#endif

C Get the corrected Higgs branching ratios from UMSSMTools :
      err=slhaRead('UMSSM_decay.dat',1)
      if(err.ne.0) stop

#ifdef OMEGA

C to exclude processes with virtual W/Z in DM   annihilation
      call  setVVdecay(0,0)                                 

C to include processes with virtual W/Z  also  in co-annihilation
C     call setVVdecay(2,2)     

      write(*,*)
      write(*,*) '==== Calculation of relic density ====='

      Omega=darkOmega(Xf,fast,Beps)
      write(*,'(a,3x,1PE10.2)') 'Omega=',Omega

      if(Omega.gt.0) OmegaAppr=printChannels(Xf,cut,Beps,1,6)

C!  See an  example of direct access to annihilation channels    
c      do 555 i=1,20
c      err=omegach(i,weight,prtcls10_4)
c      if( weight.gt.0) then   
c      write(*,*) err, weight
c      write(*,*)prtcls10_4(1),prtcls10_4(2),prtcls10_4(3),prtcls10_4(4)
c555   contine     
c      endif
C    restore default VZdecay,VWdecay
       call setVVdecay(1,1)
#endif


#ifdef  INDIRECT_DETECTION
      write(*,*) 
      write(*,*) '=============== Indirect detection ==========='

C** Emin - energy cut in GeV
      Emin=1     
      Etest=Mcdm/2
C** fi angle of sight [rad], dfi - half of cone angle    
      fi=0.1
      dfi=M_PI/180
      sigmaV=calcSpectrum(4+2+1,SpA,SpE,SpP,NULL,NULL,NULL,err)    
C  Returns  sigma*v in cm^3/sec.
C  SpX - calculated spectra of annihilation.
C  Use SpectdNdE(E, SpX) to calculate energy distribution in 1/GeV units.
      write(*,*)
      write(*,fmt="('sigma*v=',1PE8.2,'[cm^3/s]')") sigmaV

      write(*,*) 'Photons'
      call gammaFluxTab(fi,dfi, sigmaV, SpA,  FluxA)
      write(*,fmt="('Photon flux  for angle of sight ',
     < F4.2,'[rad]')") fi
      write(*,fmt="('and spherical region described by cone with angle '
     >,F7.4,'[rad]' )") 2*dfi
#ifdef SHOWPLOTS      
      write(mess,fmt="('gamma flux[cm^2 s GeV]^{-1} at fi=',F4.2,
     > ' from cone with angle=',F4.2,'[rad]')")fi,dfi*2  
      err= displaySpectrum(mess,Emin,Mcdm,FluxA)
#endif
      write(*,fmt="('Photon flux =',1PE8.2,
     >'[cm^2 s GeV]^{-1} for E=',0PF5.1,'[GeV]')")
     > SpectdNdE(Etest, FluxA), Etest
      write(*,*) 'Positrons'
      call  posiFluxTab(Emin, sigmaV, SpE, FluxE)
#ifdef SHOWPLOTS      
      err= 
     >displaySpectrum('positron flux [cm^2 s sr GeV]^{-1}',
     >Emin,Mcdm,FluxE)
#endif
      SMmev=320 
      if(SMmev.gt.0) call solarModulation(SMmev,0.0005d0,FluxE,FluxE)
      write(*,fmt="('Positron flux  = ',1PE8.2,
     >'[cm^2 sr s GeV]^{-1} for E=',0PF6.1,'[GeV]')")
     > SpectdNdE(Etest, FluxE),  Etest     

      write(*,*) 'antiProtons'
      call  pbarFluxTab(Emin, sigmaV, SpP,  FluxP)

      if(SMmev.gt.0) call solarModulation(SMmev,1.d0,FluxP,FluxP)
#ifdef SHOWPLOTS      
      err= 
     >displaySpectrum('a-proton flux [cm^2 s sr GeV]^{-1}',Emin,
     >Mcdm,FluxP)
#endif
      write(*,fmt="('Anti proton flux  = ',1PE8.2,
     >'[cm^2 sr s GeV]^{-1} for E=',0PF6.1,'[GeV]')")
     > SpectdNdE(Etest, FluxP),  Etest
#endif

#ifdef RESET_FORMFACTORS

C   The user has approach to form factors  which specifies quark contents 
C   of  proton and nucleon via global parametes like
C      <Type>FF<Nucleon><q>
C   where <Type> can be "Scalar", "pVector", and "Sigma"; 
C         <Nucleon>     "P" or "N" for proton and neutron
C         <q>            "d", "u","s"
C
C   calcScalarQuarkFF( Mu/Md, Ms/Md, sigmaPiN[MeV], sigmaS[MeV])  
C   calculates and rewrites Scalar form factors

      write(*,fmt='(" protonFF (default)  d, u, s ", 1PE10.2,E10.2,E10.2)') ScalarFFPd,ScalarFFPu,ScalarFFPs
      write(*,fmt='(" neutronFF(default)  d, u, s ", 1PE10.2,E10.2,E10.2)') ScalarFFNd,ScalarFFNu,ScalarFFNs                          
      write(*,*) '===========Redefinition of form factors========='  
      
      call calcScalarQuarkFF(0.46d0,27.5d0,34.d0,42.d0)


C    To restore default form factors of  version 2  call
C      call  calcScalarQuarkFF(0.553d0,18.9d0,55.d0,243.5d0)


      write(*,fmt='(" protonFF (new)      d, u, s ", 1PE10.2,E10.2,E10.2)') ScalarFFPd,ScalarFFPu,ScalarFFPs                           
      write(*,fmt='(" neutronFF(new)      d, u, s ", 1PE10.2,E10.2,E10.2)') ScalarFFNd,ScalarFFNu,ScalarFFNs                    

#endif


#ifdef  CDM_NUCLEON
      write(*,*)
      write(*,*) '==== Calculation of CDM-nucleons amplitudes  ====='


      err=nucleonAmplitudes(CDM1, pA0,pA5,nA0,nA5)
      write(*,*)'CDM-nucleon amplitudes:'
      write(*,fmt="(' proton:  SI=',1PE9.2,' SD=',E9.2)")pA0(1),pA5(1)
      write(*,fmt="(' neutron: SI=',1PE9.2,' SD=',E9.2)")nA0(1),nA5(1) 

C  nucleon mass
      Nmass=0.939d0

      SCcoeff=4/M_PI*3.8937966E8
     >*(Nmass*Mcdm/(Nmass+ Mcdm))**2
      write(*,*) 'CDM-nucleon cross sections[pb]:'

      write(*,fmt="(' proton  SI=',1PE9.2,' SD=',E9.2)")
     >                           SCcoeff*pA0(1)**2,3*SCcoeff*pA5(1)**2
      write(*,fmt="(' neutron SI=',1PE9.2,' SD=',E9.2)")
     >                           SCcoeff*nA0(1)**2,3*SCcoeff*nA5(1)**2
#endif


#ifdef CDM_NUCLEUS
      write(*,*)
      write(*,*)'=========  Direct Detection ==========='

      write(*,*) 'Ge:'

      nEvents=nucleusRecoil(Maxwell,73,Z_Ge,J_Ge73,
     >SxxGe73,dNdE)
C  See '../include/micromegas.h' for description of arguments
C  Instead of Maxwell  one can use 'fDvDelta' Delta-function
C  velocity distribution.
      write(*,fmt="('Total number of events=',1PE9.2,' /day/kg')") 
     >nEvents

      write(*,fmt=
     >"(' Number of events in 10 - 50 KeV region=',1PE9.2,' /day/kg')") 
     >cutRecoilResult(dNdE,10.d0,50.d0) 


#ifdef SHOWPLOTS
       err=displayRecoilPlot(dNdE,
     >           'Distribution of recoil energy of 73Ge',0.d0,199.d0)
#endif
       write(*,*)
       write(*,*) '131Xe:'  
       nEvents=nucleusRecoil(Maxwell,131,Z_Xe,J_Xe131,
     > SxxXe131,dNdE)

      write(*,fmt=
     >"('Total number of events=',1PE9.2,' /day/kg')")nEvents
      write(*,fmt="(' Number of events in 10 - 50 KeV region=',
     >   1PE9.2, ' /day/kg')")cutRecoilResult(dNdE,1.d1,5.d1)

#ifdef SHOWPLOTS
      err=displayRecoilPlot(dNdE,
     >    'Distribution of recoil energy of 131Xe',0.d0,199.d0)
#endif

#endif

#ifdef NEUTRINO 
      forSun=1
      Emin=1
C rhoDet - detector density
      rhoDet=1     
      write(*,*)  
      write(*,fmt='("===============Neutrino Telescope=======  for  ",$)') 

      if(forSun.ne.0) then 
         write(*,*) "Sun"
      else 
         write(*,*) "Earth"  
      endif 

      err=neutrinoFlux(Maxwell,forSun, FluxN,FluxNb)
#ifdef SHOWPLOTS
      err=displaySpectrum("nu flux from Sun [1/Year/km^2/GeV]",Emin,Mcdm,FluxN)
      err=displaySpectrum("nu-bar from Sun  [1/Year/km^2/GeV]",Emin,Mcdm,FluxNb)
#endif
      write(*,fmt='("E>", 1PE8.2," GeV neutrino flux =", E8.2," [1/Year/km^2]")') Emin,spectrInfo(Emin,FluxN,NULL)
      write(*,fmt='("E>", 1PE8.2," GeV anti-neutrino flux =", E8.2," [1/Year/km^2]")') Emin, spectrInfo(Emin,FluxNb,NULL) 
 
C!  ICE CUBE
      if(forSun.eq.1) then 
         write(*,fmt='("IceCube22 exclusion confidence level = ",1PE8.2,"%")')  100*exLevIC22(FluxN,FluxNb,NULL)
      endif
C ======  Upward events
  
      call muonUpward(FluxN,FluxNb,FluxMu)
#ifdef SHOWPLOTS  
      err= displaySpectrum("Upward muons[1/Year/km^2/GeV]",Emin,Mcdm/2,FluxMu)
#endif

      write(*,fmt='("E>", 1PE8.2," GeV Upward muon flux =", E8.2," [1/Year/km^2]")') Emin,spectrInfo(Emin,FluxMu,NULL)
  
C ======  Contained events 
      call muonContained(FluxN,FluxNb,rhoDet, FluxMu)
#ifdef SHOWPLOTS  
      err=displaySpectrum("Contained  muons[1/Year/km^3/GeV]",Emin,Mcdm,FluxMu)
#endif
      write(*,fmt='("E>", 1PE8.2," Contained muon flux =", E8.2," [1/Year/km^3]")') Emin,spectrInfo(Emin,FluxMu,NULL)
#endif




#ifdef DECAYS

      write(*,*)
      write(*,*)'======Calculations of widths======='

      write(*,*) 
      pname = 'h1'
      width=pWidth(pname,txtL)
      write(*,fmt="(A8, ':  total width ',1PE8.2,'[GeV]')") pname, width
      write(*,*)'and Branchings:'
      call printTxtList(txtL,6)

      write(*,*) 
      pname = '~o2'
      width=pWidth(pname,txtL)
      write(*,fmt="(A8,':  total width ',1PE8.2,'[GeV]')")pname, width
      write(*,*)'and Branchings:'
      call printTxtList(txtL,6)

#endif

#ifdef CROSS_SECTIONS
      write(*,*)
      call nextOdd(next,1, nextM)
      if((next(1:1).ne.' ') .and. (nextM.lt.1000)) then 
         call antiparticle(next_,next)
         write(*,*)'======Calculations of cross sections======='
         write(*,*) ' Process p,p -> ', next,next_
         Pcm=6500 
         Qren=nextM
         Qfact=Qren
         pTmin=0
         cs=hCollider(Pcm,1,3,Qren, Qfact, next,next_,pTmin,1)
         write(*,fmt="('total cross section  ', 1PE8.2,'[pb]')") cs 
      endif    

#endif

#ifdef CLEAN
      call system("rm -f SM_decay.dat UMSSM*")
      call system("rm -f  nngg.out output.flha ")
      call system("rm -f HB.* HS.* hb.* hs.*  debug*  Key.dat")
      call system("rm -f Lilith_*   particles.py*")
      call system("rm -f smodels.in  smodels.log  smodels.out summary.*")
#endif 

      call killPlots
      END
