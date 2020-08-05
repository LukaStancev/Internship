*DECK T16FND
      SUBROUTINE T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1 ,TKEY2 ,
     >                  NBELEM)
*
*----
*  1- PROGRAMME STATISTICS:
*      NAME     : T16DIM
*
*Purpose:
*  FIND NEXT RECORD IDENTIFIED BY
*                 TKEY1 AND TKEY2 ON TAPE16 FILE
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
*      | 1999/10/22 G.M.    | FIND NEXT RECORD IDENTIFIED BY
*      |                    | TKEY1 AND TKEY2 ON TAPE16 FILE
*      |                    | CAN ALSO BE USED TO PRINT THE LIST
*      |                    | OF RECORDS ON TAPE16 IF IPRINT=10000
*      |                    | AND TKEY1=TKEY2=' '
*      --------------------------------------------------------------
*      | 1999/11/08 GM      | MULTIPLE SET OF RECORDS KEYS PERMITTED
*      |                    | (TKEY1(1),TKEY2(1)) MUST APPEAR
*      |                    | BEFORE (TKEY1(IK),TKEY2(IK)) IK=2,NKEY
*      --------------------------------------------------------------
*
*  2- ROUTINE PARAMETERS:
*Parameters: input
* IFT16   TAPE16 FILE UNIT                         I
* IPRINT  PRINT LEVEL                              I
*         <  100 NO PRINT
*         >= 100 PRINT RECORD TO READ
*         >= 10000 PRINT ALL RECORD READ TO REACH
*               REQUESTED RECORD
* IOPT    PROCESSING OPTION                        I
*         IOPT =-1 START AT CURRENT POSITION
*                  AND READ UNTIL END OF FILE
*                  NO BACKSPACE BEFORE RETURN
*         IOPT = 0 START AT CURRENT POSITION
*                  AND READ UNTIL END OF FILE
*                  BACKSPACE BEFORE RETURN
*         IOPT = 1 REWIND FILE BEFORE READING
*                  AND READ UNTIL END OF FILE
*         IOPT = 2 REWIND FROM CURRENT POSITION
*                  TO PREVIOUS POSITION WITH
*                  A REWIND AT END OF FILE
* NKEY    NUMBER OF KEYS                           I
*         = 1 -> SEARCH FOR TKEY1(1),TKEY2(1)
*                UNTIL END OF FILE REACHED
*         > 1 -> SEARCH FOR TKEY1(1),TKEY2(1)
*                UNTIL (TKEY1(IK),TKEY2(IK),IK=2,NKEY)
*                OR END OF FILE REACHED
* TKEY1   MAIN RECORD KEY REQUIRED                 C(NKEY)*10
* TKEY2   SUB RECORD KEY REQUIRED                  C(NKEY)*10
*
*Parameters: output
* NBELEM  NUMBER OF ELEMENTS ON RECORD FOUND       I
*         NBELEM < -1 RECORD NOT FOUND BEFORE
*                     ALTERNATIVE KEY -NBELEM
*                     REACHED
*         NBELEM = -1 RECORD NOT FOUND BEFORE
*                     END OF FILE
*         NBELEM =  0 TOTAL NUMBER OF ELEMENTS ON
*                     RECORD
*
*----
*
      IMPLICIT         NONE
      INTEGER          IFT16,IPRINT,IOPT,NKEY,NBELEM
      CHARACTER        TKEY1(NKEY)*10,TKEY2(NKEY)*10
C----
C  LOCAL VARIABLES
C----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='T16FND')
      CHARACTER        RKEY1*10,RKEY2*10
      INTEGER          NBE,IEND,IKEY
C----
C  REWIND FILE FIRST IF IOPT=1
C----
      IF(IPRINT .GE. 100) THEN
        IF(IPRINT .LT. 10000) THEN
          WRITE(6,6000) TKEY1(1),TKEY2(1)
        ENDIF
      ENDIF
      IEND=1
      IF(IOPT .EQ. 1) THEN
        REWIND(IFT16)
      ELSE IF (IOPT .EQ. 2) THEN
        IEND=0
      ENDIF
C----
C  LOOP FOR READ
C----
 100  CONTINUE
      READ(IFT16,END=105) RKEY1,RKEY2,NBE
      IF(IPRINT .GE. 10000) THEN
        WRITE(6,6003) RKEY1,RKEY2,NBE
      ENDIF
      IF(RKEY1 .EQ. TKEY1(1) .AND.
     >   RKEY2 .EQ. TKEY2(1)       ) THEN
C----
C  KEYS FOUND BACKSPACE AND RETURN
C----
        NBELEM=NBE
        IF(IOPT .GE. 0) BACKSPACE(IFT16)
        IF(IPRINT .GE. 100) THEN
          WRITE(6,6001) RKEY1,RKEY2,NBELEM
        ENDIF
        RETURN
      ELSE IF(NKEY .GE. 2) THEN
        DO 110 IKEY=2,NKEY
          IF(RKEY1 .EQ. TKEY1(IKEY) .AND.
     >       RKEY2 .EQ. TKEY2(IKEY)       ) THEN
            NBELEM=-IKEY
            IF(IOPT .GE. 0) BACKSPACE(IFT16)
            IF(IPRINT .GE. 100) THEN
              WRITE(6,6004) RKEY1,RKEY2,NBE,
     >                      TKEY1(1),TKEY2(1)
            ENDIF
            RETURN
          ENDIF
 110    CONTINUE
      ENDIF
C----
C  KEYS NOT FOUND READ NEXT RECORD
C----
      GO TO 100
C----
C  END OF FILE REACHED
C----
 105  CONTINUE
      IF(IEND .EQ. 0) THEN
C----
C  REWIND FILE AND CONTINUE READ
C----
        IEND=1
        REWIND(IFT16)
        GO TO 100
      ENDIF
C----
C  RECORD ABSENT, RETURN
C----
      NBELEM=-1
      IF(IPRINT .GE. 100) THEN
        IF(IPRINT .LT. 10000) THEN
          WRITE(6,6002) TKEY1(1),TKEY2(1)
        ENDIF
      ENDIF
      RETURN
C----
C  PRINT FORMAT
C----
 6000 FORMAT( 1X, 'FIND T16 RECORD = ',2(A10,2X))
 6001 FORMAT( 1X, '     T16 RECORD = ',2(A10,2X),I10,
     >        1X,'FOUND')
 6002 FORMAT( 1X, '     T16 RECORD = ',2(A10,2X),10X,
     >        1X,'NOT FOUND')
 6003 FORMAT(11X,'T16 RECORD READ = ',2(A10,2X),I10)
 6004 FORMAT( 1X,'T16 STOP RECORD = ',2(A10,2X),I10,
     >        1X,'FOUND BEFORE RECORD = ',2(A10,2X))
      END
