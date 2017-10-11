
      implicit none

      character*20  fInput 
      character*20  fOutput
      real *8 v,vcsnngg,vcsnngz,vcsgg,vcsgz,zero
      integer nArgs,err ,slharead, fail
      
      call getarg(1,fInput)
      call getarg(2,fOutput)      
      
      OPEN(UNIT=77,FILE=fInput,STATUS='OLD')     
      call ModelVarIni(77)
      close(77)
      call ModelConstIni(fail)
      v=0D0
      vcsgg=vcsnngg(v)
      vcsgz=vcsnngz(v)      

      OPEN(UNIT=77,FILE=fOutput,STATUS='UNKNOWN')   
      write(77,fmt='("BLOCK lGamma  # AZ and AA cross sections")')
      write(77,fmt='(A6, 1PE10.4,A20)') '  1   ', vcsgz,'# ~o1,~o1->A,Z [pb]'
      write(77,fmt='(A6, 1PE10.4,A20)') '  2   ', vcsgg,'# ~o1,~o1->A,A [pb]'
      close(77)
      end 
