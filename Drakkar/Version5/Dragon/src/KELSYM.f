*DECK KELSYM
      FUNCTION KELSYM(   IPRT,   NDIM,   MAXDO,  NSURO,  NVOLO,
     >                 IDLGEO,  INDEXO, MATGEO, KEYSYM )
************************************************************************
*                                                                      *
*           NAME: KELSYM                                               *
*      COMPONENT: EXCELL                                               *
*        VERSION: 1.0                                                  *
*       CREATION: 90/01                                                *
*       MODIFIED: 00/03 (R.R.) DECLARE ALL VARIABLE TYPES              *
*         AUTHOR: ROBERT ROY                                           *
*                                                                      *
*       FUNCTION: THIS FUNCTION WILL GENERATE THE VECTOR 'KEYSYM'      *
*                 FOR A BLOCK                                          *
*                                                                      *
*          LEVEL: 4 (CALLED BY 'XELTRP')                               *
*                                                                      *
*--------+-------------- V A R I A B L E S -------------+--+-----------*
*  NAME  /                  DESCRIPTION                 /IO/MOD(DIMENS)*
*--------+----------------------------------------------+--+-----------*
* KELSYM / FUNCTION = # OF SURFACES WITH SYMMETRIC.     /.O/INT        *
* IPRT   / INTERMEDIATE PRINTING LEVEL                  /I./INT        *
* NDIM   / # OF DIMENSIONS (2 OR 3)                     /I./INT        *
* MAXDO  / MAX INDEX VALUES FOR ALL AXES (RECT/CYL).    /I./INT(NC3MAX)*
* NSURO  / # OF SURFACES FOR A SPECIFIC GEOMETRY.       /I./INT        *
* NVOLO  / # OF ZONES FOR A SPECIFIC GEOMETRY.          /I./INT        *
* IDLGEO / SPECIFIC POSITION FOR A GEOMETRY.            /I./INT        *
* INDEXO / COORDINATES FOR ZONES & SURFACES OF A CELL.  /I./INT(4,*)   *
* MATGEO / MATERIAL #S CORRESPONDING TO GEOMETRIES.     /I./INT(NGIDL) *
* KEYSYM / SYMMETRY #S CORRESPONDING TO GEOMETRIES.     /.O/INT(NGIDL) *
************************************************************************
C
      IMPLICIT        NONE
C
      INTEGER         KELSYM, IPRT, NDIM, NSURO, NVOLO, IDLGEO
      INTEGER         MAXDO(*),INDEXO(4,*),KEYSYM(*),MATGEO(*)
C
      INTEGER         ICUR(4), I, J, IVS, MAXPRC, MAXSUI, ISYM, IND
      LOGICAL         SWITCH
      INTEGER         IOUT
      PARAMETER     ( IOUT=6 )
C
      IND(I)= IDLGEO + I
C
      DO 5 IVS= 0, NVOLO
         KEYSYM(IND(IVS))= 0
    5 CONTINUE
      KELSYM= 0
C
C     LOCATES THE SYMMETRIC SURFACE TO EACH SURFACE
      DO 50 IVS = NSURO, -1
         IF( MATGEO(IND(IVS)).EQ.0 )GO TO 51
         MAXPRC= 0
         DO 10 J = 1, 4
            ICUR(J)= INDEXO(J,IND(IVS))
C
C           FIND THE SYMMETRIC SURFACE BY CHANGING END-FACE
            IF( J.LE.NDIM )THEN
               MAXSUI= MAXDO(J)
               IF( ICUR(J).EQ.MAXPRC)THEN
                  ICUR(J)= MAXSUI
               ELSEIF( ICUR(J).EQ.MAXSUI)THEN
                  ICUR(J)= MAXPRC
               ENDIF
               MAXPRC= MAXSUI
            ENDIF
C
C           THE SENTINEL VALUE IS IVS=0
            INDEXO(J,IND(0))= ICUR(J)
   10    CONTINUE
         ISYM= NSURO
   20       SWITCH= .TRUE.
            DO 30 J    = 1, 4
               SWITCH= SWITCH .AND. ICUR(J).EQ.INDEXO(J,IND(ISYM))
   30       CONTINUE
            IF( SWITCH )GO TO 40
            ISYM= ISYM + 1
         GO TO 20
   40    KEYSYM(IND(IVS))= ISYM
         IF( IPRT.GE.10 )THEN
            WRITE(IOUT,'(22H SURFACE SYMMETRIC TO ,I6,4H IS ,I6)')
     >                                         -IVS,     -ISYM
         ENDIF
         IF( ISYM.NE.0 ) KELSYM=KELSYM-1
   51    CONTINUE
   50 CONTINUE
C
C     RESET SENTINEL INDEXO(J,IND(0)) FOR SUBSEQUENT USES
      DO 60 J= 1, 4
         INDEXO(J,IND(0))= 0
   60 CONTINUE
C
      RETURN
      END
