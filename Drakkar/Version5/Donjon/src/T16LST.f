*DECK T16LST
      SUBROUTINE T16LST(IFT16 )
*
*----
*  1- PROGRAMME STATISTICS:
*      NAME     : T16DIM
*
*Purpose:
*  PRINT TAPE16 CONTENTS
*
*Author(s): 
* G.MARLEAU
*
*      CREATED  : 1999/10/22
*      REF      : EPM  IGE-244 REV.1
*                 EACL RC-1176 (COG-94-52)
*
*      MODIFICATION LOG
*      --------------------------------------------------------------
*      | DATE AND INITIALS  | MOTIVATIONS
*      --------------------------------------------------------------
*      | 1999/10/22 G.M.    | LIST THE CONTENTS OF A TAPE16 FILE
*      --------------------------------------------------------------
*
*  2- ROUTINE PARAMETERS:
*Parameters: input
* IFT16   TAPE16 FILE UNIT                         I
*
*  3- ROUTINES CALLED
*    SPECIFIC T16CPO ROUTINES
*      T16FND : FIND A TAPE16 RECORD
*               EQUIVALENT TO FIND FUNCTION
*               IN APPENDIX E OF EACL RC-1176
*
*----
*
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
