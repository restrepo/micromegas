      SUBROUTINE EWDAT_CPV(PAR,PROB)

      IMPLICIT NONE

      DOUBLE PRECISION PAR(*),PROB(*)
      DOUBLE PRECISION Pi,aux,NMB0,NMB22,NMH,CEW0,BEW1,FSF1,DEW0,DEW27
      DOUBLE PRECISION PiZMZ,PIZO,PiWMW,PiWO,Drho,rho2SM,dVB

      DOUBLE PRECISION GF,MZ,MW,g1,g2,alSMZ,S2TW,ALEMMZ
      DOUBLE PRECISION mt,mb,mtau,mmu,mel,MS,MC,MBP,MPI,MSTRANGE
      DOUBLE PRECISION tanb,cosb,sinb,vu,vd
      DOUBLE PRECISION MSL3,MSE3,MSL1,MSE1,ATAU,AMU
      DOUBLE PRECISION MSQ1,MSU1,MSD1
      DOUBLE PRECISION MSQ3,MSU3,MSD3,AT,AB
      DOUBLE PRECISION mur,M1r,M2r,msi

      COMMON/EWPAR/GF,MZ,MW,g1,g2,alSMZ,S2TW,ALEMMZ
      COMMON/SMFERM/mt,mb,mtau,mmu,mel,MS,MC,MBP,MPI,MSTRANGE
      COMMON/TBPAR/tanb,cosb,sinb,vu,vd
      COMMON/SLEPPAR/MSL3,MSE3,MSL1,MSE1,ATAU,AMU
      COMMON/SQUPAR/MSQ1,MSU1,MSD1
      COMMON/RADCOR2/MSQ3,MSU3,MSD3,AT,AB
      COMMON/GAUGINOPAR/mur,M1r,M2r,msi


      Pi=4d0*datan(1d0)

c      1) DRbar fine-structure constant

      aux=0.06887d0-1d0/(137.036d0*2d0*Pi)
     . *(-7d0*dlog(MW/MZ)+16d0/9d0*dlog(mt/MZ)+dlog(PAR(23)/MZ)/3d0
     .  +4d0/9d0*(2d0*dlog(MSQ1/MZ**2)+2d0*dlog(MSU1/MZ**2)
     .                    +dlog(MSQ3/MZ**2)+dlog(MSU3/MZ**2))/2d0
     .  +1d0/9d0*(2d0*dlog(MSQ1/MZ**2)+2d0*dlog(MSD1/MZ**2)
     .                    +dlog(MSQ3/MZ**2)+dlog(MSD3/MZ**2))/2d0
     .  +1d0/3d0*(2d0*dlog(MSL1/MZ**2)+2d0*dlog(MSE1/MZ**2)
     .                    +dlog(MSL3/MZ**2)+dlog(MSE3**2/MZ))/2d0
     .  +4d0/3d0*(dlog(M2r/MZ)+dlog(mur/MZ)))

      ALEMMZ=1d0/137.036d0/(1d0-aux)


c      2) EW gauge boson self-energies

      PiZMZ=MZ**2*NMB0(MZ**2,125d0**2,MZ**2,MZ**2)
     . -2d0*NMB22(MZ**2,125d0**2,MZ**2,MZ**2)
     . -2d0*NMB22(MZ**2,PAR(23)**2,PAR(23)**2,MZ**2)
     . -(1d0-2d0*S2TW)**2*(NMB22(MZ**2,MW**2,MW**2,MZ**2)
     .                 +NMB22(MZ**2,PAR(23)**2,PAR(23)**2,MZ**2))
     . -2d0*(1d0-S2TW)**2*(2d0*MZ**2+MW**2-MZ**2*S2TW**2/(1-S2TW))
     .                              *NMB0(MZ**2,MW**2,MW**2,MZ**2)
     . -8d0*(1d0-S2TW)**2*NMB22(MZ**2,MW**2,MW**2,MZ**2)
     . +3d0*((1d0/2d0-2d0/3d0*S2TW)**2+(2d0/3d0*S2TW)**2)
     .  *(2d0*NMH(MZ**2,0d0,0d0,MZ**2)+NMH(MZ**2,mt**2,mt**2,MZ**2))
     . +12d0*(2d0/3d0*S2TW)*(1d0/2d0-2d0/3d0*S2TW)*mt**2
     .                      *NMB0(MZ**2,mt**2,mt**2,MZ**2)
     . +3d0*((-1d0/2d0+S2TW/3d0)**2+(S2TW/3d0)**2)
     .  *(2d0*NMH(MZ**2,0d0,0d0,MZ**2)+NMH(MZ**2,mb**2,mb**2,MZ**2))
     . +12d0*(-1d0/3d0*S2TW)*(-1d0/2d0+1d0/3d0*S2TW)*mb**2
     .                      *NMB0(MZ**2,mb**2,mb**2,MZ**2)
     . +3d0*(1d0/2d0)**2*NMH(MZ**2,0d0,0d0,MZ**2)
     . +((-1d0/2d0+S2TW)**2+(S2TW)**2+(1d0/2d0)**2)
     .  *(2d0*NMH(MZ**2,0d0,0d0,MZ**2)+NMH(MZ**2,mtau**2,mtau**2,MZ**2))
     . +4d0*(-S2TW)*(-1d0/2d0+S2TW)*mtau**2
     .                      *NMB0(MZ**2,mtau**2,mtau**2,MZ**2)
     . -12d0*((1d0/2d0-2d0/3d0*S2TW)**2+(-1d0/2d0+S2TW/3d0)**2)
     .  *(2d0*NMB22(MZ**2,MSQ1,MSQ1,MZ**2)
     .                            +NMB22(MZ**2,MSQ3,MSQ3,MZ**2))
     . -12d0*(2d0/3d0*S2TW)**2
     .  *(2d0*NMB22(MZ**2,MSU1,MSU1,MZ**2)
     .                            +NMB22(MZ**2,MSU3,MSU3,MZ**2))
     . -12d0*(-1d0/3d0*S2TW)**2
     .  *(2d0*NMB22(MZ**2,MSD1,MSD1,MZ**2)
     .                            +NMB22(MZ**2,MSD3,MSD3,MZ**2))
     . -4d0*((1d0/2d0)**2+(-1d0/2d0+S2TW)**2)
     .  *(2d0*NMB22(MZ**2,MSL1,MSL1,MZ**2)
     .                            +NMB22(MZ**2,MSL3,MSL3,MZ**2))
     . -4d0*(-S2TW)**2
     .  *(2d0*NMB22(MZ**2,MSE1,MSE1,MZ**2)
     .                            +NMB22(MZ**2,MSE3,MSE3,MZ**2))
     . +(NMH(MZ**2,mur**2,mur**2,MZ**2)
     .                +2d0*mur**2*NMB0(MZ**2,mur**2,mur**2,MZ**2))/2d0
     . +2d0*(1d0-S2TW)**2*(NMH(MZ**2,M2r**2,M2r**2,MZ**2)
     .                +2d0*M2r**2*NMB0(MZ**2,M2r**2,M2r**2,MZ**2))
     . +(1d0-2d0*S2TW)**2*(NMH(MZ**2,mur**2,mur**2,MZ**2)
     .                +2d0*mur**2*NMB0(MZ**2,mur**2,mur**2,MZ**2))/2d0

      PiZMZ=(g1+g2)*PiZMZ/16d0/Pi**2

      PiZO=-MZ**2*FSF1(125d0**2,MZ**2,MZ**2)
     . -2d0*NMB22(0d0,125d0**2,MZ**2,MZ**2)
     . -2d0*NMB22(0d0,PAR(23)**2,PAR(23)**2,MZ**2)
     . -(1d0-2d0*S2TW)**2*(NMB22(0d0,MW**2,MW**2,MZ**2)
     .                 +NMB22(0d0,PAR(23)**2,PAR(23)**2,MZ**2))
     . +2d0*(1d0-S2TW)**2*(MW**2-MZ**2*S2TW**2/(1-S2TW))
     .                              *FSF1(MW**2,MW**2,MZ**2)
     . -8d0*(1d0-S2TW)**2*NMB22(0d0,MW**2,MW**2,MZ**2)
     . +3d0*((1d0/2d0-2d0/3d0*S2TW)**2+(2d0/3d0*S2TW)**2)
     .  *(2d0*0d0+NMH(0d0,mt**2,mt**2,MZ**2))
     . -12d0*(2d0/3d0*S2TW)*(1d0/2d0-2d0/3d0*S2TW)*mt**2
     .                      *FSF1(mt**2,mt**2,MZ**2)
     . +3d0*((-1d0/2d0+S2TW/3d0)**2+(S2TW/3d0)**2)
     .  *(2d0*0d0+NMH(0d0,mb**2,mb**2,MZ**2))
     . -12d0*(-1d0/3d0*S2TW)*(-1d0/2d0+1d0/3d0*S2TW)*mb**2
     .                      *FSF1(mb**2,mb**2,MZ**2)
     . +3d0*(1d0/2d0)**2*0d0
     . +((-1d0/2d0+S2TW)**2+(S2TW)**2+(1d0/2d0)**2)
     .  *(2d0*0d0+NMH(0d0,mtau**2,mtau**2,MZ**2))
     . -4d0*(-S2TW)*(-1d0/2d0+S2TW)*mtau**2
     .                      *FSF1(mtau**2,mtau**2,MZ**2)
     . -12d0*((1d0/2d0-2d0/3d0*S2TW)**2+(-1d0/2d0+S2TW/3d0)**2)
     .  *(2d0*NMB22(0d0,MSQ1,MSQ1,MZ**2)
     .                            +NMB22(0d0,MSQ3,MSQ3,MZ**2))
     . -12d0*(2d0/3d0*S2TW)**2
     .  *(2d0*NMB22(0d0,MSU1,MSU1,MZ**2)
     .                            +NMB22(0d0,MSU3,MSU3,MZ**2))
     . -12d0*(-1d0/3d0*S2TW)**2
     .  *(2d0*NMB22(0d0,MSD1,MSD1,MZ**2)
     .                            +NMB22(0d0,MSD3,MSD3,MZ**2))
     . -4d0*((1d0/2d0)**2+(-1d0/2d0+S2TW)**2)
     .  *(2d0*NMB22(0d0,MSL1,MSL1,MZ**2)
     .                            +NMB22(0d0,MSL3,MSL3,MZ**2))
     . -4d0*(-S2TW)**2
     .  *(2d0*NMB22(0d0,MSE1,MSE1,MZ**2)
     .                            +NMB22(0d0,MSE3,MSE3,MZ**2))
     . +(NMH(0d0,mur**2,mur**2,MZ**2)
     .                -2d0*mur**2*FSF1(mur**2,mur**2,MZ**2))/2d0
     . +2d0*(1d0-S2TW)**2*(NMH(0d0,M2r**2,M2r**2,MZ**2)
     .                -2d0*M2r**2*FSF1(M2r**2,M2r**2,MZ**2))
     . +(1d0-2d0*S2TW)**2*(NMH(0d0,mur**2,mur**2,MZ**2)
     .                -2d0*mur**2*FSF1(mur**2,mur**2,MZ**2))/2d0

      PiZO=(g1+g2)*PiZO/16d0/Pi**2

      PiWMW=MW**2*NMB0(MW**2,125d0**2,MW**2,MZ**2)
     . -NMB22(MW**2,125d0**2,MW**2,MZ**2)
     . -NMB22(MW**2,MZ**2,MW**2,MZ**2)
     . -2d0*NMB22(MW**2,PAR(23)**2,PAR(23)**2,MZ**2)
     . -((4d0*MW**2+MW**2+MZ**2)*(1d0-S2TW)-MZ**2*S2TW**2)
     .                        *NMB0(MW**2,MZ**2,MW**2,MZ**2)
     . -8d0*(1d0-S2TW)*NMB22(MW**2,MZ**2,MW**2,MZ**2)
     . -4d0*(2d0*NMB22(MW**2,MW**2,0d0,MZ**2)
     .                    +MW**2*NMB0(MW**2,MW**2,0d0,MZ**2))
     . +4d0*NMH(MW**2,0d0,0d0,MZ**2)
     . +3d0/2d0*NMH(MW**2,mt**2,mb**2,MZ**2)
     . +NMH(MW**2,mtau**2,0d0,MZ**2)/2d0
     . -6d0*(2d0*NMB22(MW**2,MSQ1,MSQ1,MZ**2)
     .             +NMB22(MW**2,MSQ3,MSQ3,MZ**2))
     . -2d0*(2d0*NMB22(MW**2,MSL1,MSL1,MZ**2)
     .             +NMB22(MW**2,MSL3,MSL3,MZ**2))
     . +2d0*(NMH(MW**2,M2r**2,M2r**2,MZ**2)
     .       +2d0*M2r**2*NMB0(MW**2,M2r**2,M2r**2,MZ**2))
     . +2d0*(NMH(MW**2,mur**2,mur**2,MZ**2)
     .       +mur**2*NMB0(MW**2,mur**2,mur**2,MZ**2))

      PiWMW=g2*PiWMW/16d0/Pi**2

      PiWO=-MW**2*FSF1(125d0**2,MW**2,MZ**2)
     . -NMB22(0d0,125d0**2,MW**2,MZ**2)
     . -NMB22(0d0,MZ**2,MW**2,MZ**2)
     . -2d0*NMB22(0d0,PAR(23)**2,PAR(23)**2,MZ**2)
     . +((MW**2+MZ**2)*(1d0-S2TW)-MZ**2*S2TW**2)
     .                        *FSF1(MZ**2,MW**2,MZ**2)
     . -8d0*(1d0-S2TW)*NMB22(0d0,MZ**2,MW**2,MZ**2)
     . -4d0*(2d0*NMB22(0d0,MW**2,0d0,MZ**2)
     .                    -MW**2*FSF1(MW**2,0d0,MZ**2))
     . +3d0/2d0*NMH(0d0,mt**2,mb**2,MZ**2)
     . +NMH(0d0,mtau**2,0d0,MZ**2)/2d0
     . -6d0*(2d0*NMB22(0d0,MSQ1,MSQ1,MZ**2)
     .             +NMB22(0d0,MSQ3,MSQ3,MZ**2))
     . -2d0*(2d0*NMB22(0d0,MSL1,MSL1,MZ**2)
     .             +NMB22(0d0,MSL3,MSL3,MZ**2))
     . +2d0*(NMH(0d0,M2r**2,M2r**2,MZ**2)
     .       -2d0*M2r**2*FSF1(M2r**2,M2r**2,MZ**2))
     . +2d0*(NMH(0d0,mur**2,mur**2,MZ**2)
     .       -mur**2*FSF1(mur**2,mur**2,MZ**2))

      PiWO=g2*PiWO/16d0/Pi**2


c      3) Custodial symmetry parameter

      Drho=(PiZMZ/MZ**2-PiWMW/MW**2)/(1d0+PiZMZ/MZ**2)
     . +ALEMMZ*ALSMZ/4d0/Pi**2/S2TW*(-2.145d0*(mt/MZ)**2
     .   +1.262d0*dlog(mt/MZ)-2.24d0-0.55*(MZ/mt))
     . +3*GF**2*mt**4/128d0/Pi**4*rho2SM(125d0/mt)


c      4) Vertex / Box contributions to the mu-nu_mu-e-nu_e eff. coupling

      dVB=ALEMMZ/4.d0/Pi/S2TW/(1d0-Drho)*(6d0+dlog(MW**2/MZ**2)
     . /(1d0-MW**2/MZ**2)*(1d0+5d0/2d0*MW**2/MZ**2
     . -(1d0-MW**2/MZ**2)*(5d0-3d0/2d0*MW**2/MZ**2/(1-S2TW))))

      aux=g1/2d0*BEW1(M1r**2,MSL1,MZ**2)
     . +3d0*g2/2d0*BEW1(M2r**2,MSL1,MZ**2)
     . +2d0*g2*M2r**2*CEW0(MSL1,M2r**2,M2r**2)
     . +g2*(FSF1(M2r**2,M2r**2,MZ**2)
     .          -MSL1*CEW0(MSL1,M2r**2,M2r**2)+1d0/2d0)
     . -g1/2d0/dsqrt(2d0)*(FSF1(MSL1,MSL1,MZ**2)
     .          -M1r**2*CEW0(M1r**2,MSL1,MSL1)-1d0/2d0)
     . +g2/2d0/dsqrt(2d0)*(FSF1(MSL1,MSL1,MZ**2)
     .          -M2r**2*CEW0(M2r**2,MSL1,MSL1)-1d0/2d0)

      dVB=dVB-2d0*aux/16d0/Pi**2

      aux=g1*g2/2d0*M1r*M2r*DEW0(MSL1,M1r**2,M2r**2)
     . -(g2*M2r)**2/2d0*DEW0(MSL1,M2r**2,M2r**2)
     . +g1*g2*DEW27(MSL1,M1r**2,M2r**2)
     . +g2**2*DEW27(MSL1,M2r**2,M2r**2)

      dVB=dVB-S2TW*(1d0-S2TW)/2d0/Pi/ALEMMZ*MZ**2/16d0/Pi**2*aux


c      5) EW vev's and couplings
c      print*,vu,vd,g1,g2
c      vu=vu/dsqrt(1d0-PiWO/MW**2-dVB)
c      vd=vd/dsqrt(1d0-PiWO/MW**2-dVB)

c      g2=2d0*MW**2/(vu**2+vd**2)*(1d0+PiWMW/MW**2)
c      g1=2d0*MZ**2/(vu**2+vd**2)*(1d0+PiZMZ/MZ**2)-g2
c      print*,vu,vd,g1,g2

      RETURN
      END

**********************************************************************

      DOUBLE PRECISION function NMB22(x,y,z,t)

c            ->NMB22(p^2,m1^2,m2^2,Q^2)
      
      IMPLICIT NONE
      DOUBLE PRECISION x,y,z,t,aux,NMB0,NMA0,FSF1

      IF(x.ge.1d-4)THEN
      aux=y+z-x/3d0-(x-2d0*y-2d0*z)*NMB0(x,y,z,t)/2d0+NMA0(y,t)/2d0
     .      +NMA0(z,t)/2d0+(z-y)/2d0/x*(NMA0(z,t)-NMA0(y,t)
     .                                     -(z-y)*NMB0(x,y,z,t))
      ELSE
       IF(Max(y,z).ge.1d-4)THEN
      IF(min(y,z).ge.1d-4)THEN
      IF(dabs(z-y).ge.1d-4)THEN
      aux=y+z-(y+z)*FSF1(y,z,t)+NMA0(y,t)/2d0
     .      +NMA0(z,t)/2d0+y*z*dlog(z/y)/2d0/(z-y)-(y+z)/4d0
      ELSE
      aux=y+z-(y+z)*FSF1(y,z,t)+NMA0(y,t)/2d0+NMA0(z,t)/2d0
      ENDIF
      ELSE
       aux=-Max(y,z)/4d0
      ENDIF
       ELSE
      aux=0d0
       ENDIF
      ENDIF

      NMB22=aux/6d0-(NMA0(y,t)+NMA0(z,t))/4d0

      RETURN
      END

**********************************************************************

      DOUBLE PRECISION function NMH(x,y,z,t)

c            ->NMH(p^2,m1^2,m2^2,Q^2)
      
      IMPLICIT NONE
      DOUBLE PRECISION x,y,z,t,aux,NMB0,NMA0,NMB22,FSF1

      IF(x.ge.1d-4)THEN
      aux=(x-y-z)*NMB0(x,y,z,t)
      ELSE
       IF(Max(y,z).ge.1d-4)THEN
      aux=(y+z)*FSF1(y,z,t)
       ELSE
      aux=0d0
       ENDIF
      ENDIF

      NMH=4d0*NMB22(x,y,z,t)+aux

      RETURN
      END

**********************************************************************

      DOUBLE PRECISION function rho2SM(x)
      
c            ->rho2SM(x)

      IMPLICIT NONE
      DOUBLE PRECISION x,Pi

      Pi=4d0*datan(1d0)

      rho2SM=19d0-33d0*x/2d0+43d0*x**2/12d0+7d0*x**3/120d0
     . -Pi*dsqrt(x)*(4d0-3d0/2d0*x+3d0/32d0*x**2+x**3/256d0)
     . -Pi**2*(2d0-2d0*x+x**2/2d0)-dlog(x)*(3d0*x-x**2/2d0)

      RETURN
      END

**********************************************************************

      DOUBLE PRECISION function BEW1(x,y,z)
      
c            ->BEW1(m1^2,m2^2,Q^2)

      IMPLICIT NONE
      DOUBLE PRECISION x,y,z,aux

      IF(min(x,y).ge.1d-4)THEN
      IF(dabs(x-y).ge.1d-4)THEN
       aux=1d0+dlog(z/y)+(x/(x-y))**2*dlog(y/x)+(x+y)/(x-y)/2d0
      ELSE
       aux=dlog(z/x)
      ENDIF
      ELSE
       IF(x.ge.1d-4)aux=dlog(z/x)+3d0/2d0
       IF(y.ge.1d-4)aux=dlog(z/y)+1d0/2d0
      ENDIF

      BEW1=aux/2d0

      RETURN
      END

**********************************************************************

      DOUBLE PRECISION function CEW0(x,y,z)
      
c            ->CEW0(m1^2,m2^2,m3^2)

      IMPLICIT NONE
      DOUBLE PRECISION x,y,z,aux

      IF(min(x,y,z).ge.1d-4)THEN
      IF(dabs(y-z).ge.1d-4)THEN
       IF(dabs(x-y).ge.1d-4)THEN
      IF(dabs(x-z).ge.1d-4)THEN
      aux=(y*dlog(y/x)/(x-y)-z*dlog(z/x)/(x-z))/(y-z)
      ELSE
      aux=(y*dlog(x/y)-x+y)/(x-y)**2
      ENDIF
       ELSE
      IF(dabs(x-z).ge.1d-4)THEN
      aux=(z*dlog(x/z)-x+z)/(x-z)**2
      ELSE
      aux=-1/2d0/x
        ENDIF
       ENDIF
      ELSE
       IF(dabs(x-y).ge.1d-4)THEN
      aux=(x*dlog(y/x)+x-y)/(x-y)**2
       ELSE
      aux=-1/2d0/x
       ENDIF
      ENDIF
      ELSE
       IF(x.ge.1d-4.and.y.ge.1d-4)THEN
      IF(dabs(y-x).ge.1d-4)THEN
      aux=dlog(y/x)/(x-y)
      ELSE
      aux=-1/x
      ENDIF
       ENDIF
       IF(x.ge.1d-4.and.y.ge.1d-4)THEN
      IF(dabs(y-x).ge.1d-4)THEN
      aux=dlog(y/x)/(x-y)
      ELSE
      aux=-1/x
      ENDIF
       ENDIF
       IF(x.ge.1d-4.and.z.ge.1d-4)THEN
      IF(dabs(z-x).ge.1d-4)THEN
      aux=dlog(z/x)/(x-z)
      ELSE
      aux=-1/x
      ENDIF
       ENDIF
       IF(y.ge.1d-4.and.z.ge.1d-4)THEN
      IF(dabs(z-y).ge.1d-4)THEN
      aux=dlog(z/y)/(y-z)
      ELSE
      aux=-1/y
      ENDIF
       ENDIF
      ENDIF

      CEW0=aux

      RETURN
      END

**********************************************************************

      DOUBLE PRECISION function DEW0(x,y,z)
      
c            ->DEW0(m1^2,m2^2,m3^2)

      IMPLICIT NONE
      DOUBLE PRECISION x,y,z,aux

      IF(min(x,y,z).ge.1d-4)THEN
       IF(dabs(y-z).ge.1d-4)THEN
      IF(dabs(x-y).ge.1d-4)THEN
       IF(dabs(x-z).ge.1d-4)THEN
        aux=-(y*(x-y+x*dlog(y/x))*(x-z)**2
     . -z*(x-z+x*dlog(z/x))*(x-y)**2)/x/(y-z)/(x-y)**2/(x-z)**2
       ELSE
        aux=(x**2-y**2+2d0*x*y*dlog(y/x))/2d0/x/(x-y)**3
       ENDIF
      ELSE
       IF(dabs(x-z).ge.1d-4)THEN
        aux=(x**2-z**2+2d0*x*z*dlog(z/x))/2d0/x/(x-z)**3
       ELSE
        aux=1d0/6d0/x**2
       ENDIF
      ENDIF
       IF(dabs(x-y).ge.1d-4)THEN
        aux=-(2d0*(x-y)+(x+y)*dlog(x/y))/(x-y)**3
       ELSE
        aux=1d0/6d0/x**2
       ENDIF
       ENDIF
      ELSE
       IF(x.ge.1d-4)THEN
      IF(max(y,z).ge.1d-4)THEN
         aux=(-x+Max(y,z)-x*dlog(Max(y,z)/x))/x/(x-Max(y,z))**2
      ENDIF
       ENDIF
      ENDIF

      DEW0=aux

      RETURN
      END

**********************************************************************

      DOUBLE PRECISION function DEW27(x,y,z)
      
c            ->DEW27(m1^2,m2^2,m3^2)

      IMPLICIT NONE
      DOUBLE PRECISION x,y,z,aux

      IF(min(x,y,z).ge.1d-4)THEN
       IF(dabs(y-z).ge.1d-4)THEN
      IF(dabs(x-y).ge.1d-4)THEN
       IF(dabs(x-z).ge.1d-4)THEN
        aux=-(y**2*(x-z)**2*dlog(y/x)-z**2*(x-y)**2*dlog(z/x)
     .        +x*(x-y)*(x-z)*(y-z))/4d0/(y-z)/(x-y)**2/(x-z)**2
       ELSE
        aux=-(x**2-4d0*x*y+3d0*y**2-2d0*y**2*dlog(y/x))
     .          /8d0/(x-y)**3
       ENDIF
      ELSE
       IF(dabs(x-z).ge.1d-4)THEN
        aux=-(x**2-4d0*x*z+3d0*z**2-2d0*z**2*dlog(z/x))
     .          /8d0/(x-z)**3
       ELSE
        aux=-1d0/12d0/x
       ENDIF
      ENDIF
       IF(dabs(x-y).ge.1d-4)THEN
        aux=(-x**2+y**2-2d0*x*y*dlog(y/x))/4d0/(x-y)**3
       ELSE
        aux=-1d0/12d0/x
       ENDIF
       ENDIF
      ELSE
       IF(x.ge.1d-4)THEN
      IF(Max(y,z).ge.1d-4)THEN
         aux=(x-Max(y,z)+Max(y,z)*dlog(Max(y,z)/x))/4d0
     .          /(x-Max(y,z))**2
      ELSE
       aux=-1/4d0/x
      ENDIF
       ELSE
      IF(min(y,z).ge.1d-4)THEN
       IF(dabs(y-z).ge.1d-4)THEN
        aux=-dlog(y/z)/4d0/(y-z)
       ELSE
        aux=-1/4d0/y
       ENDIF
      ENDIF
       ENDIF
      ENDIF

      DEW27=aux

      RETURN
      END
