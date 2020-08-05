*DECK T16REC
      SUBROUTINE T16REC(IFT16 ,IPRINT,INEXTR)
*
*----
*  1- PROGRAMME STATISTICS:
*      NAME     : T16REC
*
*Purpose:
*  LOCATE NEXT SET OF RECORDS
*
*Author(s): 
* G.MARLEAU
*
*      CREATED  : 1999/10/21
*      REF      : IGE-244 REV.1
*
*      MODIFICATION LOG
*      --------------------------------------------------------------
*      | DATE AND INITIALS  | MOTIVATIONS
*      --------------------------------------------------------------
*      | 1999/12/17 G.M.    | EXTRACTED FROM T16FLX
*      --------------------------------------------------------------
*
*  2- ROUTINE PARAMETERS:
*Parameters: input
* IFT16   TAPE16 FILE UNIT                         I
* IPRINT  PRINT LEVEL                              I
*         =   0 NO PRINT
*         >=  1 PRINT PROCESSING OPTIONS READ
* INEXTR  NEXT RECORD TO READ                      I
*
*  3- ROUTINES CALLED
*    SPECIFIC T16CPO ROUTINES
*      T16FND : FIND A TAPE16 RECORD
*               EQUIVALENT TO FIND FUNCTION
*               IN APPENDIX E OF EACL RC-1176
*    UTILITIES ROUTINES
*      XABORT : ABORT ROUTINE
*      XDRSET : VECTOR INITIALIZATION ROUTINE
*
*----
*
      IMPLICIT         NONE
      INTEGER          IFT16,IPRINT,INEXTR
C----
C  T16 PARAMETERS
C----
      INTEGER          MAXKEY
      PARAMETER       (MAXKEY=3)
      CHARACTER        TKEY1(MAXKEY)*10,TKEY2(MAXKEY)*10,
     >                 RKEY1*10,RKEY2*10
      INTEGER          NKEY,IOPT,NBE,NSKIPR,ISKIPR
C----
C  LOCAL VARIABLES
C----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='T16REC')
C----
C  IF NEXT SET OF RECORD INEXTR <= LAST SET OF RECORD
C  REWIND AND SKIP FIRST INEXTR-1 SETS OF RECORDS
C----
      REWIND(IFT16)
      NSKIPR=INEXTR
      TKEY1(1)='MTR       '
      TKEY2(1)='FEWGROUPS '
      NKEY=1
      IOPT=-1
      DO 100 ISKIPR=1,NSKIPR
        CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1 ,TKEY2 ,
     >              NBE   )
        IF(NBE .EQ. -1) THEN
          WRITE(IOUT,9000) NAMSBR,TKEY1(1),TKEY2(1),INEXTR
          CALL XABORT(NAMSBR//
     >    ': INVALID RECORD NUMBER ON TAPE16')
        ENDIF
        READ(IFT16) RKEY1,RKEY2,NBE
 100  CONTINUE
      RETURN
C----
C  ABORT FORMAT
C----
 9000 FORMAT(1X,A6,1X,7('*'),' ERROR ',7('*')/
     >       8X,I6,' TAPE16 RECORD WITH KEYS =',2(A10,2X),
     >       'NOT FOUND'/
     >       8X,21('*'))
      END
