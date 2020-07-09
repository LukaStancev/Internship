*DECK T16LST
      SUBROUTINE T16LST(IFT16 )
C
C----
C  1- PROGRAMME STATISTICS:
C      NAME     : T16DIM
C      USE      : PRINT TAPE16 CONTENTS
C      AUTHOR   : G.MARLEAU
C      CREATED  : 1999/10/22
C      REF      : EPM  IGE-244 REV.1
C                 EACL RC-1176 (COG-94-52)
C
C      MODIFICATION LOG
C      --------------------------------------------------------------
C      | DATE AND INITIALS  | MOTIVATIONS
C      --------------------------------------------------------------
C      | 1999/10/22 G.M.    | LIST THE CONTENTS OF A TAPE16 FILE
C      --------------------------------------------------------------
C
C  2- ROUTINE PARAMETERS:
C    INPUT
C      IFT16  : TAPE16 FILE UNIT                         I
C  3- ROUTINES CALLED
C    SPECIFIC T16CPO ROUTINES
C      T16FND : FIND A TAPE16 RECORD
C               EQUIVALENT TO FIND FUNCTION
C               IN APPENDIX E OF EACL RC-1176
C
C----
C
      IMPLICIT         NONE
      INTEGER          IFT16
C----
C  T16 KEYS
C----
      CHARACTER        TKEY1*10,TKEY2*10
      INTEGER          NKEY,IOPT,NBE
C----
C  LOCAL VARIABLES
C----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='T16LST')
      INTEGER          IPRINT
C----
C  LIST TAPE16 RECORDS AFTER REWINDING
C----
      WRITE(IOUT,6000) NAMSBR
      IPRINT=10000
      IOPT=1
      NKEY=1
      TKEY1='          '
      TKEY2=TKEY2
      CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1  ,TKEY2  ,
     >            NBE)
      WRITE(IOUT,6001)
      RETURN
C----
C  PRINT FORMAT
C----
 6000 FORMAT( 1X, 'PRINTING CONTENTS OF TAPE16 FILE USING ',A6)
 6001 FORMAT( 1X, 'END OF TAPE16 FILE REACHED')
      END
