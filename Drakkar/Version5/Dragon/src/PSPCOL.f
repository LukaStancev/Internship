*DECK PSPCOL
      SUBROUTINE PSPCOL(ITCOL,NCOL,ICOL,RGB)
C
C---------------------------  PSPCOL  ---------------------------------
C
C  1- PROGRAMME STATISTICS:
C      NAME     : PSPCOL
C      USE      : PICK A COLOR NUMBER FROM A N-COLOR SET
C      AUTHOR   : G.MARLEAU
C      CREATED  : 99-01-21
C
C      MODIFICATION LOG
C      --------------------------------------------------------------
C      | DATE AND INITIALS  | MOTIVATIONS
C      --------------------------------------------------------------
C      | 99-01-21 G.M.      | PICK RGB COLOR IN A DISTRIBUTED
C      |                    | FASHION
C      ______________________________________________________________
C
C  2- ROUTINE PARAMETERS:
C    INPUT
C      ITCOL    : TYPE OF COLOR SET                      I
C                 = 1 GRAY
C                 = 2 RGB
C                 = 3 CMYK
C                 = 4 HSB
C      NCOL     : MAXIMUM NUMBER OF COLOR IN SET         I
C      ICOL     : REQUESTED COLOR NUMBER                 I
C    OUTPUT
C      RGB      : COLOR INTENSITY                        R(4)
C                 FOR GRAY USE ONLY RGB(1)
C                 FOR RGB USE ONLY RGB(1),RGB(2),RGB(3)
C                 FOR CMYK USE ALL
C                 FOR HSB USE ONLY RGB(1),RGB(2),RGB(3)
C
C---------------------------   PSPCOL  --------------------------------
C
      IMPLICIT         NONE
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='PSPCOL')
C----
C  ROUTINE PARAMETERS
C----
      INTEGER          ITCOL,NCOL,ICOL
      REAL             RGB(4)
C----
C  LOCAL PARAMETERS
C----
      INTEGER          IDC,JCOL
      REAL             DELCOL,DELSAT,DELBLK
C----
C  LOCAL VARIABLES
C----
      IF(ITCOL .EQ. 4) THEN
        RGB(4)=0.0
        IF(ICOL .LE. 0 ) THEN
          RGB(1)=0.0
          RGB(2)=0.0
          RGB(3)=1.0
        ELSE
          DELCOL=0.6667/FLOAT(NCOL-1)
          DELSAT=0.5/FLOAT(NCOL-1)
          DELBLK=0.5/FLOAT(NCOL-1)
          JCOL=ICOL-1
          RGB(1)=0.6667-DELCOL*FLOAT(JCOL)
          RGB(2)=0.5+DELSAT*FLOAT(JCOL)
          RGB(3)=0.5+DELBLK*FLOAT(JCOL)
        ENDIF
      ELSE IF(ITCOL .EQ. 3) THEN
        RGB(4)=0.0
        IF(ICOL .LE. 0 ) THEN
          RGB(1)=0.0
          RGB(2)=0.0
          RGB(3)=0.0
        ELSE
          IF     (NCOL .LE.       8) THEN
            IDC=2
          ELSE IF(NCOL .LE.      64) THEN
            IDC=4
          ELSE IF(NCOL .LE.     512) THEN
            IDC=8
          ELSE IF(NCOL .LE.    4096) THEN
            IDC=16
          ELSE IF(NCOL .LE.   32768) THEN
            IDC=32
          ELSE IF(NCOL .LE.  262144) THEN
            IDC=64
          ELSE
            IDC=128
          ENDIF
          JCOL=ICOL-1
          DELCOL=1.0/FLOAT(IDC)
          RGB(1)=1.0-DELCOL*FLOAT(MOD(JCOL,IDC))
          JCOL=JCOL/IDC
          RGB(2)=1.0-DELCOL*FLOAT(MOD(JCOL,IDC))
          JCOL=JCOL/IDC
          RGB(3)=1.0-DELCOL*FLOAT(MOD(JCOL,IDC))
        ENDIF
      ELSE IF(ITCOL .EQ. 2) THEN
        RGB(4)=0.0
        IF(ICOL .LE. 0 ) THEN
          RGB(1)=1.0
          RGB(2)=1.0
          RGB(3)=1.0
        ELSE
          IF     (NCOL .LE.       8) THEN
            IDC=2
          ELSE IF(NCOL .LE.      64) THEN
            IDC=4
          ELSE IF(NCOL .LE.     512) THEN
            IDC=8
          ELSE IF(NCOL .LE.    4096) THEN
            IDC=16
          ELSE IF(NCOL .LE.   32768) THEN
            IDC=32
          ELSE IF(NCOL .LE.  262144) THEN
            IDC=64
          ELSE
            IDC=128
          ENDIF
          JCOL=ICOL-1
          DELCOL=1.0/FLOAT(IDC)
          RGB(1)=1.0-DELCOL*FLOAT(MOD(JCOL,IDC)+1)
          JCOL=JCOL/IDC
          RGB(2)=1.0-DELCOL*FLOAT(MOD(JCOL,IDC)+1)
          JCOL=JCOL/IDC
          RGB(3)=1.0-DELCOL*FLOAT(MOD(JCOL,IDC)+1)
        ENDIF
      ELSE
        IF(ICOL .LE. 0 ) THEN
          RGB(1)=0.0
          RGB(2)=0.0
          RGB(3)=0.0
        ELSE
          IF     (NCOL .LE.       8) THEN
            IDC=8
          ELSE IF(NCOL .LE.      64) THEN
            IDC=64
          ELSE IF(NCOL .LE.     512) THEN
            IDC=512
          ELSE IF(NCOL .LE.    4096) THEN
            IDC=4096
          ELSE IF(NCOL .LE.   32768) THEN
            IDC=32768
          ELSE
            IDC=262144
          ENDIF
          JCOL=ICOL-1
          DELCOL=1.0/FLOAT(IDC)
          RGB(1)=1.0-DELCOL*FLOAT(MOD(JCOL,IDC))
          RGB(2)=RGB(1)
          RGB(3)=RGB(1)
        ENDIF
      ENDIF
      RETURN
      END
