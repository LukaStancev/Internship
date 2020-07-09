*DECK T16REC
      SUBROUTINE T16REC(IFT16 ,IPRINT,INEXTR)
C
C----
C  1- PROGRAMME STATISTICS:
C      NAME     : T16REC
C      USE      : LOCATE NEXT SET OF RECORDS
C      AUTHOR   : G.MARLEAU
C      CREATED  : 1999/10/21
C      REF      : IGE-244 REV.1
C
C      MODIFICATION LOG
C      --------------------------------------------------------------
C      | DATE AND INITIALS  | MOTIVATIONS
C      --------------------------------------------------------------
C      | 1999/12/17 G.M.    | EXTRACTED FROM T16FLX
C      --------------------------------------------------------------
C
C  2- ROUTINE PARAMETERS:
C    INPUT
C      IFT16  : TAPE16 FILE UNIT                         I
C      IPRINT : PRINT LEVEL                              I
C               =   0 NO PRINT
C               >=  1 PRINT PROCESSING OPTIONS READ
C      INEXTR : NEXT RECORD TO READ                      I
C
C  3- ROUTINES CALLED
C    SPECIFIC T16CPO ROUTINES
C      T16FND : FIND A TAPE16 RECORD
C               EQUIVALENT TO FIND FUNCTION
C               IN APPENDIX E OF EACL RC-1176
C    UTILITIES ROUTINES
C      XABORT : ABORT ROUTINE
C      XDRSET : VECTOR INITIALIZATION ROUTINE
C
C----
C
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
