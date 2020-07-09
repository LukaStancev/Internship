*DECK XELCMP
      SUBROUTINE XELCMP(     NS,     NV,  VOLIN,  MATIN,  MRGIN,
     >                    NSOUT,  NVOUT, VOLOUT, MATOUT,  ITGEO,  ICODE)
************************************************************************
*                                                                      *
*           NAME: XELCMP                                               *
*      COMPONENT: EXCELL                                               *
*          LEVEL: 2 (CALLED BY 'EXCELT')                               *
*        VERSION: 1.0                                                  *
*       CREATION: 91/07                                                *
*       MODIFIED: 97/11 (G.M.) ELIMINATE XABORT FOR MRGIN()=0          *
*                 00/03 (R.R.) DECLARE ALL VARIABLE TYPES              *
*         AUTHOR: ROBERT ROY                                           *
*                                                                      *
*     SUBROUTINE: THIS ROUTINE WILL MERGE VOLUMES AND SURFACES         *
*                 AND RECOMPUTE THE NUMBER OF SURFACES AND VOLUMES.    *
*                                                                      *
*--------+-------------- V A R I A B L E S -------------+--+-----------*
*  NAME  /                  DESCRIPTION                 /IO/MOD(DIMENS)*
*--------+----------------------------------------------+--+-----------*
* NS     / # OF SURFACES           BEFORE MERGING.      /I./INT        *
* NV     / # OF ZONES              BEFORE MERGING.      /I./INT        *
* VOLIN  / VOLUMES & SURFACES      BEFORE MERGING.      /I./INT(-NS:NV)*
* MATIN  / #ING OF SUFACES & ZONES BEFORE MERGING.      /I./INT(-NS:NV)*
* MRGIN  / MERGING INDEX.                               /I./INT(-NS:NV)*
* NSOUT  / # OF SURFACES            AFTER MERGING.      /.O/INT        *
* NVOUT  / # OF ZONES               AFTER MERGING.      /.O/INT        *
* VOLOUT / VOLUMES & SURFACES       AFTER MERGING.      /.O/INT(*)     *
* MATOUT / #ING OF SUFACES & ZONES  AFTER MERGING.      /.O/INT(*)     *
* ITGEO  / KIND OF GEOMETRY(0,1,2,3).                   /I./INT        *
* ICODE  / INDEX OF BOUNDARY CONDITIONS.                /I./INT(6)     *
************************************************************************
      IMPLICIT      NONE
C
      INTEGER       NS,NV,NSOUT,NVOUT,ITGEO,IVS,IMR,ICNT,I0,IOUT,
     >              IR,JR,MATMRG,LESOIR,
     >              MATIN(-NS:NV),MRGIN(-NS:NV),MATOUT(*),ICODE(6)
      REAL          VOLIN(-NS:NV),VOLOUT(*),ZERO
      CHARACTER*4   CORIEN(0:3,-6:0)
      PARAMETER    ( ZERO= 0.0, IOUT=6 )
      DATA         ((CORIEN(JR,IR),IR=-6,0),JR=0,3)
     >       / ' O6 ',' O5 ',' O4 ',' O3 ',' O2 ',' O1 ','    ',
     >         ' Z+ ',' Z- ','****','****',' R+ ','****','    ',
     >         ' Z+ ',' Z- ','****','****','****','HBC ','    ',
     >         ' Z+ ',' Z- ',' Y+ ',' Y- ',' X+ ',' X- ','    '/
C
C     FIND NSOUT AND NVOUT & INITIALIZE VOLOUT AND MATOUT
      NSOUT= 0
      NVOUT= 0
      DO 10 IVS= -NS, NV
         VOLOUT(IVS+NS+1)= ZERO
         MATOUT(IVS+NS+1)= 0
         IF( IVS.GT.0 )THEN
            IF( MRGIN(IVS).LT.0 )THEN
               CALL XABORT( 'XELCMP: 1.INCOMPATIBLE MERGE INDEX' )
            ENDIF
         ELSEIF( IVS.LT.0 )THEN
            IF( MRGIN(IVS).GT.0 )THEN
               CALL XABORT( 'XELCMP: 2.INCOMPATIBLE MERGE INDEX' )
            ENDIF
         ELSE
            IF( MRGIN(IVS).NE.0 )THEN
               WRITE(IOUT,*) 'XELCMP: *KEYMRG* VECTOR IS:', MRGIN
               CALL XABORT( 'XELCMP: 3.INCOMPATIBLE MERGE INDEX' )
            ENDIF
            IF( VOLIN(IVS).NE.0.0 )THEN
               WRITE(IOUT,*) 'XELCMP: *VOLSUR* VECTOR IS:', VOLIN
               CALL XABORT( 'XELCMP: 4. VOLSUR(0).NE.0 ON TRACK-FILE' )
            ENDIF
            IF( MATIN(IVS).NE.0 )THEN
               WRITE(IOUT,*) 'XELCMP: *MATALB* VECTOR IS:', MATIN
               CALL XABORT( 'XELCMP: 5. MATALB(0).NE.0 ON TRACK-FILE' )
            ENDIF
         ENDIF
         NSOUT= MIN(NSOUT,MRGIN(IVS))
         NVOUT= MAX(NVOUT,MRGIN(IVS))
   10 CONTINUE
      NSOUT= -NSOUT
C
C     ALL VALUES MUST BE PRESENT BETWEEN -NSOUT AND NVOUT IN MRGIN(*)
C     BUT WITH THE SAME MATIN(*) NUMBER FOR MERGED ZONES.
C     NEW(97/11): 0 MEANS REGION IS REMOVED
      DO 30 IMR= -NSOUT, NVOUT
         ICNT= 0
         DO 20 IVS= -NS, NV
            IF( ICNT.EQ.0 ) MATMRG= MATIN(IVS)
            IF( MRGIN(IVS).EQ.IMR )THEN
               ICNT= ICNT+1
               IF( MATMRG.NE.MATIN(IVS) )THEN
                  LESOIR= MATIN(IVS)
                  IF( IVS.GE.0 )THEN
C
C                    FOR MERGING ZONES, ABORT IF NOT SAME *MATALB*
                     WRITE(IOUT,*) '*** ABORT *** ATTEMPT TO MERGE '//
     >                          'MIX ',MATMRG,' WITH MIX ',
     >                                 LESOIR,' IN ZONE #',IVS
                     CALL XABORT( 'XELCMP: 6.INCOMPATIBLE MERGE INDEX' )
                  ELSE
C
C                    FOR MERGING FACES, ABORT IF NOT SAME *ICODE*
                     IF( ICODE(-MATMRG).NE.ICODE(-LESOIR) )THEN
                        WRITE(IOUT,*) '*** ABORT *** ATTEMPT TO MERGE ',
     >                                ' FACE ',-IVS,
     >                             '( ',CORIEN(ITGEO,MATMRG),',ICODE=',
     >                                  ICODE(-MATMRG),') WITH A FACE ',
     >                             '( ',CORIEN(ITGEO,LESOIR),',ICODE=',
     >                                  ICODE(-LESOIR),'). '
                     CALL XABORT( 'XELCMP: 7.INCOMPATIBLE MERGE INDEX' )
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
   20    CONTINUE
         IF( ICNT.EQ.0 )THEN
            CALL XABORT( 'XELCMP: 8.MISSING VALUES IN THE MERGE INDEX' )
         ENDIF
   30 CONTINUE
C
C     COMPUTE VOLOUT AND MATOUT VALUES
      I0= 1 + NSOUT
      DO 40 IVS= -NS, NV
         VOLOUT(I0+MRGIN(IVS))= VOLOUT(I0+MRGIN(IVS))+VOLIN(IVS)
         MATOUT(I0+MRGIN(IVS))= MATIN(IVS)
   40 CONTINUE
C
      RETURN
      END
