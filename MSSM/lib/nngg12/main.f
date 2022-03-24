
      implicit none

      character*100  fInput 
      character*100  fOutput
      character*100  buff
      real *8 v,vcsnngg,vcsnngz,vcsgg,vcsgz,MW,SW,EE,Mb
      integer  ModelConstIni,Nch
      integer nArgs,err ,slharead
      external ModelConstIni
      
      MW=80.385
      SW=0.487
      EE=0.3123 
      Mb=3      
      nArgs=iargc()
C      write(*,*) nArgs
      if(nArgs.lt.2) stop 
     >'2 or 3 arguments  expected for lGamma.exe: in/out files, GI key'     
      call getarg(1,fInput) 
      call getarg(2,fOutput)
      if(nArgs.ge.3) then 
         call getarg(3,buff)
         read(buff,*) MW
      endif    

      if(nArgs.ge.4) then 
         call getarg(4,buff)
         read(buff,*) SW
      endif    

      if(nArgs.ge.5) then 
         call getarg(5,buff)
         read(buff,*) EE
      endif    

      if(nArgs.ge.6) then 
         call getarg(6,buff)
         read(buff,*) Mb
      endif    
      
      err=slharead(fInput,1)
      if(err.ne.0)  stop 'Problem in reading  input file'
      Nch=ModelConstIni(MW,SW,EE,Mb)
      v=0D0   
      vcsgg=vcsnngg(v)
      if(Nch.gt.1) then    
        vcsgz=vcsnngz(v)
      else 
        vcsgz=0
      endif          
 
      OPEN(UNIT=78,FILE=fOutput,STATUS='UNKNOWN')   
      write(78,fmt='("BLOCK lGamma  # AZ and AA cross sections")')
      write(78,fmt='(A6, 1PE10.4,A20)') '  1   ', vcsgz,'# ~o1,~o1->A,Z [pb]'
      write(78,fmt='(A6, 1PE10.4,A20)') '  2   ', vcsgg,'# ~o1,~o1->A,A [pb]'
      close(78)
      end 
