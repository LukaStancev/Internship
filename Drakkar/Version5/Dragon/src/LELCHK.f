*DECK LELCHK
      LOGICAL FUNCTION LELCHK(  NSOLD,  NVOLD, VOLOLD, MATOLD,  ICOLD,
     >                          NSNEW,  NVNEW, VOLNEW, MATNEW,  ICNEW,
     >                           IPRT )
************************************************************************
*                                                                      *
*           NAME: LELCHK                                               *
*      COMPONENT: EXCELL                                               *
*          LEVEL: 2 (CALLED BY 'EXCELT')                               *
*        VERSION: 1.0                                                  *
*       CREATION: 91/08                                                *
*       MODIFIED: 00/03 (R.R.) DECLARE ALL VARIABLE TYPES              *
*         AUTHOR: R. ROY                                               *
*                                                                      *
*     SUBROUTINE: THIS ROUTINE WILL CHECK COMPATIBILITY BETWEEN DATA   *
*                 IN THE OLD TRACKING FILE AND IN THE NEW GEOMETRY.    *
*                 THIS ROUTINE DOES NOT STOP THE EXECUTION.            *
*                                                                      *
*--------+-------------- V A R I A B L E S -------------+--+-----------*
*  NAME  /                  DESCRIPTION                 /IO/MOD(DIMENS)*
*--------+----------------------------------------------+--+-----------*
* NSOLD  / # OF SURFACES             IN TRACKING FILE.  /I./INT        *
* NVOLD  / # OF ZONES                IN TRACKING FILE.  /I./INT        *
* VOLOLD / VOLUMES & SURFACES        IN TRACKING FILE.  /I./REL(-NS:NV)*
* MATOLD / #ING OF SURFACES & ZONES  IN TRACKING FILE.  /I./INT(-NS:NV)*
* ICOLD  / INDEX OF B.C.             IN TRACKING FILE.  /I./INT(6)     *
* NSNEW  / # OF SURFACES             IN NEW GEOMETRY.   /I./INT        *
* NVNEW  / # OF ZONES                IN NEW GEOMETRY.   /I./INT        *
* VOLNEW / VOLUMES & SURFACES        IN NEW GEOMETRY.   /I./REL(-NS:NV)*
* MATNEW / #ING OF SURFACES & ZONES  IN NEW GEOMETRY.   /I./INT(-NS:NV)*
* ICNEW  / INDEX OF B.C.             IN NEW GEOMETRY.   /I./INT(6)     *
* IPRT   / PRINTING LEVEL ( 0: NO PRINT)                /I./INT        *
************************************************************************
      IMPLICIT    NONE
C
      INTEGER     NSOLD,NVOLD,MATOLD(-NSOLD:NVOLD),ICOLD(6),IPRT,IOUT,
     >            NSNEW,NVNEW,MATNEW(-NSOLD:NVOLD),ICNEW(6),IR,NERROC
      REAL        VOLOLD(-NSOLD:NVOLD),VOLNEW(-NSNEW:NVNEW),
     >            ZERO,HUND,EMAX
      PARAMETER ( IOUT=6, ZERO=0.0, HUND=100.0, EMAX=1.E-5 )
      LELCHK= .TRUE.
C
C1.1) CHECK # OF ZONES ------------------------------------------------
      IF( NVOLD.NE.NVNEW )THEN
         IF( IPRT.GT.0 )THEN
            WRITE(IOUT,'(40H *** INCONSISTENT # OF ZONES            )')
         ENDIF
         LELCHK=.FALSE.
         GO TO 999
      ENDIF
C
C1.2) CHECK # OF FACES ------------------------------------------------
      IF( NSOLD.NE.NSNEW )THEN
         IF( IPRT.GT.0 )THEN
            WRITE(IOUT,'(40H *** INCONSISTENT # OF FACES            )')
         ENDIF
         LELCHK=.FALSE.
         GO TO 999
      ENDIF
C
C1.3) CHECK CONSISTENCY OF INDEX *ICODE* ------------------------------
      DO 10 IR= 1, 6
         IF( ICOLD(IR).NE.ICNEW(IR) )THEN
            IF( IPRT.GT.0 )THEN
               WRITE(IOUT,'(9H   ICODE(,I1,3H)= ,I6,5H(WAS ,I6,1H))')
     >                                  IR,      ICNEW(IR), ICOLD(IR)
            ENDIF
            IF( ICOLD(IR).LE.0.OR.ICNEW(IR).LE.0 )THEN
               LELCHK=.FALSE.
               GO TO 999
            ENDIF
         ENDIF
   10 CONTINUE
C
C1.4) CHECK IF SOME FACES HAVE ICODE=0 --------------------------------
      DO 20 IR= -NSOLD, -1
         IF( ICNEW(-MATNEW(IR)).EQ.0 )THEN
            IF( IPRT.GT.0 )THEN
               WRITE(IOUT,'(9H    FACE(,I1,3H)= ,I6,12H HAS ICODE=0 )')
     >                                 -IR,      MATNEW(IR)
            ENDIF
            LELCHK=.FALSE.
            GO TO 999
         ENDIF
   20 CONTINUE
C
C2)   CHECK CONSISTENCY OF VECTORS *VOLSUR* AND *MATALB* --------------
      NERROC= 0
      DO 30 IR= -NSOLD, NVOLD
         IF( VOLOLD(IR)-VOLNEW(IR).GT.ZERO )THEN
            NERROC= NERROC+1
            IF( IR.EQ.0 ) GO TO 30
            LELCHK= LELCHK.AND.
     >              ABS((VOLNEW(IR)-VOLOLD(IR))/VOLOLD(IR)).LE.EMAX
         ENDIF
         IF( MATOLD(IR).NE.MATNEW(IR) )THEN
            NERROC= NERROC+1
            IF( IR.LE.0 ) LELCHK= .FALSE.
         ENDIF
   30 CONTINUE
      IF( IPRT.GT.0 )THEN
         WRITE(IOUT,'(1H )')
         IF( NERROC.EQ.0 )THEN
            WRITE(IOUT,'(60H ECHO = >>> CONSISTENCY BETWEEN '//
     >                 'TRACKING FILE AND GEOMETRY                 /)')
         ELSE
            WRITE(IOUT,'(60H ECHO = >>> WARNING: INCONSISTENT '//
     >                 'TRACKING FILE                              /)')
            DO 40 IR= -NSOLD, NVOLD
               IF( IR.EQ.0 ) GO TO 40
               IF( VOLOLD(IR)-VOLNEW(IR).GT.ZERO )THEN
               IF( IR.LE.0 )THEN
                  WRITE(IOUT,'(15H ERROR ON FACE(,I4,3H)= ,F10.7,1H%)')
     >                      -IR,HUND*(VOLNEW(IR)-VOLOLD(IR))/VOLOLD(IR)
               ELSE
                  WRITE(IOUT,'(15H ERROR ON ZONE(,I4,3H)= ,F10.7,1H%)')
     >                       IR,HUND*(VOLNEW(IR)-VOLOLD(IR))/VOLOLD(IR)
               ENDIF
               ENDIF
               IF( MATOLD(IR).NE.MATNEW(IR) )THEN
               IF( IR.LE.0 )THEN
                  WRITE(IOUT,'(9H    FACE(,I1,3H)= ,I6,5H(WAS ,I6,1H))')
     >                                    -IR,     MATNEW(IR),MATOLD(IR)
               ELSE
                  WRITE(IOUT,'(9H MIXTURE(,I1,3H)= ,I6,5H(WAS ,I6,1H))')
     >                                     IR,     MATNEW(IR),MATOLD(IR)
               ENDIF
               ENDIF
   40       CONTINUE
         ENDIF
      ENDIF
C
  999 RETURN
      END
