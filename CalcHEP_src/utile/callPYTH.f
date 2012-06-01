      IMPLICIT NONE 

      INTEGER IMSS
      REAL*8 RMSS
      COMMON/PYMSSM/IMSS(0:99),RMSS(0:99)
      SAVE /PYMSSM/
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,IDUP(MAXNUP),
     &ISTUP(MAXNUP),MOTHUP(2,MAXNUP),ICOLUP(2,MAXNUP),PUP(5,MAXNUP),
     &VTIMUP(MAXNUP),SPINUP(MAXNUP)
       
      integer SCANDIR 
      INTEGER N,NSUB,NEV,MAXEVENTS
      real*8 cs
      CHARACTER  slhaFile*200 

      NSUB=SCANDIR('.')
      write(*,*) 'directory ./ ', NSUB, ' subprocesses found' 
C!  You can include  events stored in different  direstories
C!  Just add  SCANDIR commands.   


      CALL EVENTSTAT(CS, MAXEVENTS)
      write(*,*) 'total cross section ', CS
      write(*,*) 'Max number of events ', MAXEVENTS
C! Pogram can generate as may events as you ask.   
C! But statistical indepemdent will be only  'MAXEVENTS' ones.
      IF(CS.EQ.0) then 
        write(*,*) 'It looks like no generated events found'
        stop
      ENDIF

C!   number of events requested
       NEV=100
  
*  SUSY input if it needs, or empty string
      slhaFile='slha_1.txt'
      slhaFile=' '
      if(slhaFile.ne.' ') then
        IMSS(1)=11
        IMSS(21)=21
        OPEN(21,FILE=slhaFile,STATUS='OLD',ERR=100)
      else 
        IMSS(1)=1
      endif

      IF(NEV.GT.MAXEVENTS) THEN
        write(*,*) 'Only ', MAXEVENTS, ' events available'
      ENDIF

      CALL PYINIT('USER',' ',' ',0d0)
      DO 200 N=1,NEV
         CALL PYEVNT
*...<<<<<<< USER's routines for analysis should be call here >>>>>> 
c      call pylist(5)
      call pylist(7)    ! This is an example of user's routine:
      call pylist(1)    ! This is an example of user's routine:
c        call PYHEPC(1)      ! Convert PYJETS event to HEPEVT format.
*...<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
200   CONTINUE 
201   CALL closeEvents
      CLOSE(21)
      CALL PYSTAT(1)
100   continue
      END
C ============= Auxilary routins ========== 
      REAL *8 FUNCTION ID2MASS(ID)
      IMPLICIT NONE
      INTEGER PYCOMP,ID,KCHG
      REAL*8 PMAS,PARF,VCKM
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      ID2MASS=PMAS(PYCOMP(ID),1)      
c      ID2MASS=0
      END

      SUBROUTINE UPINIT
      IMPLICIT NONE
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,IDUP(MAXNUP),
     &ISTUP(MAXNUP),MOTHUP(2,MAXNUP),ICOLUP(2,MAXNUP),PUP(5,MAXNUP),
     &VTIMUP(MAXNUP),SPINUP(MAXNUP)

      INTEGER MAXPUP
      PARAMETER (MAXPUP=100)
      INTEGER IDBMUP,PDFGUP,PDFSUP,IDWTUP,NPRUP,LPRUP
      DOUBLE PRECISION EBMUP,XSECUP,XERRUP,XMAXUP
      COMMON/HEPRUP/IDBMUP(2),EBMUP(2),PDFGUP(2),PDFSUP(2),
     &IDWTUP,NPRUP,XSECUP(MAXPUP),XERRUP(MAXPUP),XMAXUP(MAXPUP),
     &LPRUP(MAXPUP)
      INTEGER I
      SAVE
* Only one subprocess is supported
      NPRUP= 1
      LPRUP(1)=1

c SHG defaults are used 
      DO 100 I=1,2
      PDFSUP(I)=-1
100   PDFGUP(I)=-1

      AQEDUP=-1
      AQCDUP=-1

      IDWTUP = 3
      XERRUP(1) = 0.
      XMAXUP(1) = 1.

      RETURN
      END
