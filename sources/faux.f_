      subroutine fortreread(Nch,fname)
      integer Nch,Nrd,i
      character*(*) fname
      character*200 buff
      character*10 fff
      logical open

      
      Nrd=12
4     Nrd=Nrd+1 
      INQUIRE(UNIT=Nrd,OPENED=open)
      if(open) goto 4       

      open(Nrd,FILE=fname,STATUS='OLD')
9     continue
        read(Nrd,fmt='(a200)',end=11) buff
        do i=200,1,-1
          if(buff(i:i) .ne. ' ') goto 10
        enddo
10      continue
        if(i.eq.0) i=1
        if(i.lt.10) then
           write(fff,fmt='(A2,I1,A1)')  '(a',i,')'
        else 
          if(i.lt.100) then 
            write(fff,fmt='(A2,I2,A1)') '(a',i,')' 
          else 
            write(fff,fmt='(A2,I3,A1)') '(a',i,')'
          endif
        endif    
        write(Nch,fmt=fff) buff
      goto 9
11    close(Nrd)
      end  
