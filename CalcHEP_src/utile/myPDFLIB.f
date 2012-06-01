
      SUBROUTINE PDFSET(PARM,VAL)
      CHARACTER*20 PARM(3)
      REAL*8 VAL(3)
      PARAMETER (NPTYMX = 4, NGRMAX = 9, NSETMX = 100)
      COMMON/W505120/ NPGSMX(NPTYMX,NGRMAX)
      CHARACTER*8 SFNAME(NPTYMX,NGRMAX,NSETMX)
      COMMON/W505110/ SFNAME
      COMMON/W50511/  NPTYPE,NGROUP,NSET,MODE,NFL,LO,TMAS
      INTEGER I,J,K 
      
      if(PARM(1).EQ.'Init0') then
C  Initialization 
        DO I=1,NPTYMX
        DO J=1,NGRMAX
           NPGSMX(I,J)=0
        DO K=1,NSETMX
          SFNAME(I,J,K)='        '
        ENDDO
        ENDDO
        ENDDO
C  Now we add one distribution
C  CalcHEP works only with type=1 distributions, 
C thus the first parameter should be 1   
        NPGSMX(1,1)=1
        SFNAME(1,1,1)='MYpdf_1'
CC To add second one 
*       NPGSMX(1,1)=2
*       SFNAME(1,1,2)='MYpdf_2'
      ELSE
C  Save numeration of distribution                      
         NPTYPE=VAL(1)
         NGROUP=VAL(2)
         NSET=VAL(3)
      ENDIF  
      END

      SUBROUTINE STRUCTM(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GL)
C
C  *********************************************************************
C  *                                                                   *
C  *   Main steering routine for all sets of structure functions       *
C  *                                                                   *
C  *                                                                   *
C  *   Input:    X     = x value of parton                             *
C  *             SCALE = QCD scale in GeV                              *
C  *                                                                   *
C  *   Output:   UPV   = up valence quark                              *
C  *             DNV   = down valence quark                            *
C  *             USEA  = sea (up_bar)                                  *
C  *             DSEA  = sea (down_bar)                                *
C  *             STR   = strange quark                                 *
C  *             CHM   = charm quark                                   *
C  *             BOT   = bottom quark                                  *
C  *             TOP   = top quark                                     *
C  *             GL    = gluon                                         *
C  *********************************************************************
C   If you have several sets, the function is characterized by     
      COMMON/W50511/NPTYPE,NGROUP,NSET,MODE,NFL,LO,TMAS 
      REAL*8 X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GL
C   Return values are has to present x*dN/dx
C   a stupid example:  
      UPV   = x
      DNV   = 0.5*x
      USEA  = 0.02
      DSEA  = 0.02
      STR   = 0.02
      CHM   = 0.02
      BOT   = 0.02
      TOP   = 0.
      GL    = 0.02
      END 

      REAL*8 FUNCTION ALPHAS2(SCALE)
      REAL*8 SCALE  
C   If you have several sets, the function is characterized by     
      COMMON/W50511/NPTYPE,NGROUP,NSET,MODE,NFL,LO,TMAS 

      ALPHAS2=0.1 
      
      END 
