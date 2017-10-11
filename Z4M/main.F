
C====== Modules =================
C   Keys to switch on
C   various modules of micrOMEGAs
C================================
   
#define MASSES_INFO
C    Display information about SUSY and Higgs masses
#define CONSTRAINTS
C    Strongly depends on the model. 
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

C#define NEUTRINO
#define DECAYS
C#define CROSS_SECTIONS
C===== Options ========

C#define SHOWPLOTS
C     Display  graphical plots on the screen

#define CLEAN to crean intermediate files

C===== End of DEFINE  settings ===== 




      IMPLICIT NONE
C! external functions 
#include "../sources/micromegas.fh"
#include "lib/pmodel.fh"
*      external iargc
*      integer iargc
C! local  Parameters   

      INTEGER fast/1/
      REAL*8  Beps/1.E-4/,cut/0.01/
C Variables      
      INTEGER err
      REAL*8 Omega,OmegaAppr,Xf
      character *100 argv, mess
      character*20 cdmName
      integer spin2,charge3,cdim,PDGcdm
      real*8  Emin,sigmaV,fi,dfi,Etest
      real*8 SpA(250),Spe(250),SpP(250) 
      real*8 FluxA(250),FluxE(250),FluxP(250) 
C for CDM_NUCLEON
      real*8  pA0(2),pA5(2),nA0(2),nA5(2),Nmass,SCcoeff
C  for CDM_NUCLEUS
      real*8  dNdE(200),rho,nEvents

C  for NEUTRINO 
      integer forSun
      real*8 rhoDet
      real*8 FluxN(NZ),FluxNb(NZ),FluxMu(NZ) 

C for DECAYS
      integer*8  txtL
      real*8 width   
      character*20 pname

C For CROSS_SECTION
      real*8 Pcm,cosmin, cosmax,cs
      real*8 mass(4)
      character*6 name(4)
      integer*8 cc
      integer l,ntot,nin,nout

      

C   to save the SLHA  input/output file  uncomment  
C*      call delFiles(0)
 
      if(iargc().eq.0) then
         write(*,*)'Correct usage: ./main  <file with parameters>'
         write(*,*)'Example: ./main  data1.par' 
         stop
      endif

      call getarg(1,argv)                               
       
      err=readVar(argv)

      if(err.eq.-1) then
          write(*,*) 'Can not open the file'
          stop
      endif

      if(err.gt.0) then
          write(*,*)'Wrong file contents at line ',err
          stop
      endif

C  To print list od variables uncomment   
C      CALL printVar(6)

      err=sortOddParticles(cdmName)
      if(err.ne.0) then
        write(*,*) 'Can t calculate ', cdmName
        STOP
      endif 
      if(CDM1.ne.' ') then  
         PDGcdm=qNumbers(CDM1,spin2, charge3,cdim)
         write(*,fmt="('Dark matter candidate is ', A5,'with spin=', I1,'/2 mass=',1PE10.2)") CDM1, spin2,Mcdm1
         if(charge3.ne.0) then 
            write(*,fmt="('DM has electric charge',I1,'/3')") charge3
            STOP
         endif
         if(cdim.ne.1) then 
            write(*,fmt="('DM has color charge cdim=',I1,'/3')") cdim
            STOP
         endif
      endif
      if(CDM2.ne.' ') then
         PDGcdm=qNumbers(CDM2,spin2, charge3,cdim)
         write(*,fmt="('Dark matter candidate is ', A5,'with spin=', I1,'/2 mass=',1PE10.2)") CDM2, spin2,Mcdm2
         if(charge3.ne.0) then 
            write(*,fmt="('DM has electric charge',I1,'/3')") charge3
            STOP
         endif
         if(cdim.ne.1) then 
            write(*,fmt="('DM has color charge cdim=',I1,'/3')") cdim
            STOP
         endif
      endif


#ifdef MASSES_INFO  
      CALL printHiggs(6)
      write(*,*) '=== MASSES OF HIGGS AND ODD SECTOR: ===' 
      CALL printMasses(6,1)
#endif

#ifdef CONSTRAINTS

#endif



#ifdef OMEGA
      write(*,*)
      write(*,*) '==== Calculation of relic density ====='
      write(*,*) 'Dark Matter candidate is ', cdmName                                                                                                             
      Omega=darkOmega2(fast,Beps)
      write(*,'(a,3x,1PE10.2)') 'omega1=',Omega*(1-fracCDM2)
      write(*,'(a,3x,1PE10.2)') 'omega2=',Omega*fracCDM2 
      call displayFunc(vs1100F, Tend, Tstart,'vs1100') 
#endif


#ifdef  INDIRECT_DETECTION
      write(*,*) 
      write(*,*) '=============== Indirect detection ==========='

C** Emin - energy cut in GeV
      Emin=1     
      Etest=Mcdm/2
C** fi angle of sight [rad], dfi - half of cone angle    
      fi=0.1
      dfi=0.05
      sigmaV=calcSpectrum(4,SpA,SpE,SpP,NULL,NULL,NULL,err)    
C  Returns  sigma*v in cm^3/sec.
C  SpX - calculated spectra of annihilation.
C       Use SpectdNdE(E, SpX) to calculate energy distribution in  1/GeV units.
C       
C       First parameter 1-includes W/Z polarization
C                       2-includes gammas for 2->2+gamma
C                       4-print cross sections             


      write(*,*)
      write(*,fmt="('sigma*v=',1PE8.2,'[cm^3/s]')") sigmaV

      write(*,*) 'Photons'
      call gammaFluxTab(fi,dfi, sigmaV, SpA, FluxA)
      write(*,fmt="('Phonton flux  for angle of sight ',F4.2,'[rad]')") fi
      write(*,fmt="('and spherical region described by cone with angle ',F4.2,'[rad]' )") 2*dfi
#ifdef SHOWPLOTS      
      write(mess,fmt="('gamma flux[cm^2 s GeV]^{-1} at fi=',F4.2,' from cone with angle=',F4.2,'[rad]')")fi,dfi*2  
      err= displaySpectrum(mess,Emin,Mcdm,FluxA)
#endif
      write(*,fmt="('Photon flux =',1PE8.2,'[cm^2 s GeV]^{-1} for E=',0PF5.1,'[GeV]')")
     > SpectdNdE(Etest, FluxA), Etest
 
      write(*,*) 'Positrons'
      call  posiFluxTab(Emin, sigmaV, SpE,  FluxE)

#ifdef SHOWPLOTS      
      err=displaySpectrum('positron flux [cm^2 s sr GeV]^{-1}',
     >Emin,Mcdm,FluxE)
#endif
      write(*,fmt="('Positron flux  = ',1PE8.2,
     >'[cm^2 sr s GeV]^{-1} for E=',0PF6.1,'[GeV]')")
     > SpectdNdE(Etest, FluxE),  Etest     

      write(*,*) 'antiProtons'
      call  pbarFluxTab(Emin, sigmaV, SpP,  FluxP)

#ifdef SHOWPLOTS      
      err=displaySpectrum('a-proton flux [cm^2 s sr GeV]^{-1}',
     >Emin,Mcdm,FluxP)
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
C   calcScalarFF( Mu/Md, Ms/Md, sigmaPiN[MeV], sigma0[MeV])  
C   calculates and rewrites Scalar form factors
      write(*,*) 
      write(*,*) '===== RESET_FORMFACTOR ============'

      write(*,fmt='(" protonFF (default)  d, u, s ", 1PE10.2,E10.2,E10.2)') ScalarFFPd,ScalarFFPu,ScalarFFPs
      write(*,fmt='(" neutronFF(default)  d, u, s ", 1PE10.2,E10.2,E10.2)') ScalarFFNd,ScalarFFNu,ScalarFFNs                          

      write(*,*) '===========Redefinition of form factors========='  
      write(*,fmt='(" protonFF (new)      d, u, s ", 1PE10.2,E10.2,E10.2)') ScalarFFPd,ScalarFFPu,ScalarFFPs                           
      write(*,fmt='(" neutronFF(new)      d, u, s ", 1PE10.2,E10.2,E10.2)') ScalarFFNd,ScalarFFNu,ScalarFFNs                    

#endif


#ifdef  CDM_NUCLEON
      write(*,*)
      write(*,*) '==== Calculation of CDM-nucleons amplitudes  ====='   

      if(CDM1.ne.' ') then 
         err=nucleonAmplitudes(CDM1, pA0,pA5,nA0,nA5)
         write(*,*) CDM1,  ' -nucleon amplitudes:'
         write(*,fmt="(' proton:  SI=',1PE9.2,' SD=',E9.2)")pA0(1),pA5(1)
         write(*,fmt="(' neutron: SI=',1PE9.2,' SD=',E9.2)")nA0(1),nA5(1) 
      endif 

      if(CDM2.ne.' ') then 
         err=nucleonAmplitudes(CDM2, pA0,pA5,nA0,nA5)
         write(*,*) CDM2,  ' -nucleon amplitudes:'
         write(*,fmt="(' proton:  SI=',1PE9.2,' SD=',E9.2)")pA0(1),pA5(1)
         write(*,fmt="(' neutron: SI=',1PE9.2,' SD=',E9.2)")nA0(1),nA5(1) 
      endif 

 
C  nucleon mass
      Nmass=0.939d0

      SCcoeff=4/M_PI*3.8937966E8
     >*(Nmass*Mcdm/(Nmass+ Mcdm))**2
      write(*,*) cdmName , '-nucleon cross sections[pb]:'
    
      write(*,fmt="(' proton  SI=',1PE9.2,' SD=',E9.2)")
     >                           SCcoeff*pA0(1)**2,3*SCcoeff*pA5(1)**2
      write(*,fmt="(' neutron SI=',1PE9.2,' SD=',E9.2)")
     >                           SCcoeff*nA0(1)**2,3*SCcoeff*nA5(1)**2

      write(*,*) 'anti-',cdmName,  ' -nucleon amplitudes:'
      write(*,fmt="(' proton:  SI=',1PE9.2,' SD=',E9.2)")pA0(2),pA5(2)
      write(*,fmt="(' neutron: SI=',1PE9.2,' SD=',E9.2)")nA0(2),nA5(2) 


      write(*,*) 'anti-',cdmName , '- nucleon cross sections:'
    
      write(*,fmt="(' proton  SI=',1PE9.2,' SD=',E9.2)")
     >                           SCcoeff*pA0(2)**2,3*SCcoeff*pA5(2)**2
      write(*,fmt="(' neutron SI=',1PE9.2,' SD=',E9.2)")
     >                           SCcoeff*nA0(2)**2,3*SCcoeff*nA5(2)**2

#endif


#ifdef CDM_NUCLEUS
      write(*,*)

      write(*,*) '========= Direct Detection ==============='

      write(*,*) 'CDM-73Ge:'
      nEvents=nucleusRecoil(Maxwell,73,Z_Ge,J_Ge73,SxxGe73,dNdE)
      write(*,fmt="('73Ge: Total number of events=',1PE10.2,
     >' /day/kg')") nEvents
      write(*,fmt="(' Number of events in 10 - 50 KeV region=',
     >1PE10.2,' /day/kg')") cutRecoilResult(dNdE,10.D0,50.D0)

#ifdef SHOWPLOTS
       err=displayRecoilPlot(dNdE,
     >               'Distribution of recoil energy of 73Ge',0,199)
#endif
       write(*,*) 
       write(*,*) 'CDM-131Xe:' 
       nEvents=nucleusRecoil(Maxwell,131,Z_Xe,J_Xe131,SxxXe131,dNdE)
      write(*,fmt="('131Xe: Total number of events=',1PE10.2,
     >' /day/kg')") nEvents
      write(*,fmt="(' Number of events in 10 - 50 KeV region=',
     >1PE10.2,' /day/kg')") cutRecoilResult(dNdE,10.D0,50.D0)

#ifdef SHOWPLOTS
      err=displayRecoilPlot(dNdE,
     >    'Distribution of recoil energy of 131Xe',0,199);
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
      write(*,fmt='("E>", 1PE8.2," GeV anti-neutrino flux =", E8.2," [1/Year/km^2]")') Emin,spectrInfo(Emin,FluxNb,NULL) 

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

      call pdg2name(23,pname)
      if(pname .ne. ' ') then
         width=pWidth(pname,txtL)
         write(*,fmt="(A5,':  total width ',1PE8.2,'[GeV]')") pname,width
         write(*,*)'and Branchings:'
         call printTxtList(txtL,6)
      endif 

      call pdg2name(25,pname)
      if(pname .ne. ' ') then
         width=pWidth(pname,txtL)
         write(*,fmt="(A5,' :  total width ',1PE8.2,'[GeV]')") pname,width
         write(*,*)'and Branchings:'
         call printTxtList(txtL,6)
      endif 

#endif

#ifdef CROSS_SECTIONS
      write(*,*)
      write(*,*)'======Calculations of cross sections======='

      Pcm=500
      cosmin=-0.99
      cosmax=0.99
      write(*,fmt="('Process e,E->2*x at Pcm= ',F6.1,'[GeV]')") Pcm
      call setHelicities(0.5D0,-0.5D0)
      cc= newProcess('e%,E%->2*x')
      if(cc.ne.0 ) then
        call procInfo1(cc,ntot,nin,nout)
        do l=1,ntot
          cs= cs22(cc,l,Pcm,cosmin,cosmax,err)
          if((err.eq.0) .and. (cs.ne.0)) then 
            call procInfo2(cc,l,name,mass)
            write(*,fmt='(A3,1x,A3,A4,A3,1x,A3,$)')
     >      name(1),name(2),' -> ', name(3),name(4)
            write(*,fmt="(1PE8.2,'[pb]')") cs
          endif 
        enddo
      endif     
#endif


#ifdef CLEAN

#endif

      call killPlots
      END
