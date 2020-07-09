*DECK PSPFIL
      SUBROUTINE PSPFIL(ISPSP,JSPSP,NAMPSP,NPAGE)
C
C---------------------------  PSPFIL  ---------------------------------
C
C  1- PROGRAMME STATISTICS:
C      NAME     : PSPFIL
C      USE      : PSP FILE ANLYSIS
C      AUTHOR   : G.MARLEAU
C      CREATED  : 99-01-21
C
C      MODIFICATION LOG
C      --------------------------------------------------------------
C      | DATE AND INITIALS  | MOTIVATIONS
C      --------------------------------------------------------------
C      | 99-01-21 G.M.      | IF OLD PSP FILE : TEST  COMPATIBILITY
C      |                    | IF NEW PSP FILE : INITIALIZE
C      ______________________________________________________________
C
C  2- ROUTINE PARAMETERS:
C    INPUT/OUTPUT
C      ISPSP    : PSP FILE UNIT                          I
C      JSPSP    : PSP FILE MODE                          I
C                 = 0 NEW
C                 = 1 UPDATE
C      ISPSP    : PSP FILE UNIT                          I
C      NAMPSP   : PSP FILE NAME                          C*12
C      NPAGE    : PAGE NUMBER                            I
C  3- ADDITIONAL COMMENTS
C      THIS ROUTINE BYPASSES THE FOLLOWING PSPLOT ROUTINES
C        NEWDEV -> OPENS THE POSTSCRIPT OUTPUT FILE
C        PSINIT -> INITIALIZE THE POSTSCRIPT FILE
C        CHOPIT -> PREPARE A NEW POSTSCRIPT PAGE FOR AN
C                  OLD POSTSCRIPT FILE
C
C---------------------------   PSPFIL  --------------------------------
C
      IMPLICIT         NONE
      INTEGER          IOUT
      CHARACTER        NAMSBR*6,PROGNM*6
      PARAMETER       (IOUT=6,NAMSBR='PSPFIL',PROGNM='DRAGON')
C----
C  ROUTINE PARAMETERS
C----
      INTEGER          JSPSP,ISPSP,NPAGE
      CHARACTER        NAMPSP*12
C----
C  LOCAL VARIABLES
C----
      INTEGER          IRL,IDR,ILINE,IPF,IPN
      CHARACTER        CMDSTR*132,CFMT*16
      REAL             XYPOS(2)
      NPAGE=0
      IF(JSPSP .EQ. 1) THEN
C----
C  TEST IF ADEQUATE DRAGON PS FILE TYPE
C----
        DO 100 IRL=1,3
          READ(ISPSP,'(A132)') CMDSTR
 100    CONTINUE
        READ(ISPSP,'(A132)') CMDSTR
        IDR=INDEX(CMDSTR,PROGNM)
        IF(IDR .EQ. 0) CALL XABORT(NAMSBR//
     >    ': NOT A DRAGON GENERATED POSTSCRIPT FILE')
        ILINE=0
        IPF=1
C----
C  LOCATE LAST PAGE NUMBER
C----
 110    CONTINUE
          READ(ISPSP,'(A132)',END=115) CMDSTR
          IPN=INDEX(CMDSTR,'%%Page')
          IF(IPN .NE. 0) THEN
            IPN=INDEX(CMDSTR,' ')
            IPF=INDEX(CMDSTR(IPN+1:132),' ')-1
            CFMT=' '
            WRITE(CFMT,'(2H(I,I1,1H))') IPF
            READ(CMDSTR(IPN+1:IPN+IPF),CFMT) NPAGE
          ENDIF
          GO TO 110
 115    CONTINUE
        BACKSPACE ISPSP
C----
C  SET NEXT PAGE NUMBER AND PREPARE FOR OUTPUT
C----
        NPAGE=NPAGE+1
        XYPOS(1)=0.5
        XYPOS(2)=0.5
        CALL PSPAGE(ISPSP,NPAGE,XYPOS)
      ELSE
        CALL PSHEAD(ISPSP,NAMPSP,PROGNM)
        NPAGE=NPAGE+1
        XYPOS(1)=0.5
        XYPOS(2)=0.5
        CALL PSPAGE(ISPSP,NPAGE,XYPOS)
      ENDIF
      RETURN
      END
