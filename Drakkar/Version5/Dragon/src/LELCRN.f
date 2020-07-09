*DECK LELCRN
      FUNCTION LELCRN( CENTEC, RAYONC, X, Y)
************************************************************************
*                                                                      *
*           NAME: LELCRN                                               *
*      COMPONENT: EXCELL                                               *
*          LEVEL: 5 (CALLED BY 'KELRNG' & 'XELVOL')                    *
*        VERSION: 1.0                                                  *
*       CREATION: 90/02                                                *
*       MODIFIED: 00/03 (R.R.) DECLARE ALL VARIABLE TYPES              *
*         AUTHOR: ROBERT ROY                                           *
*                                                                      *
*     SUBROUTINE: THIS ROUTINE WILL DECIDE IF THE CROWN INTERSECT      *
*                 A RECTANGULAR MESH.                                  *
*                                                                      *
*--------+-------------- V A R I A B L E S -------------+--+-----------*
*  NAME  /                  DESCRIPTION                 /IO/MOD(DIMENS)*
*--------+----------------------------------------------+--+-----------*
* CENTEC / COORDIANTES OF CENTER.                       /I./DBL(2)     *
* RAYONC / INNER AND OUTER RADIUS**2 OF THE CROWN.      /I./DBL(2)     *
* X      / X OF THE SQUARE.                             /I./DBL(2)     *
* Y      / Y OF THE SQUARE.                             /I./DBL(2)     *
* LELCRN / .T. IF INTERSECTION EXISTS.                  /.O/LOGICAL    *
************************************************************************
C
      IMPLICIT NONE
      LOGICAL  LELCRN
C
      DOUBLE PRECISION CENTEC(2), RAYONC(2), X(2), Y(2), R
      INTEGER NBEXT, NBINT, IX, IY
C
      NBEXT=0
      NBINT=0
      DO 10 IX=1, 2
      DO 10 IY=1, 2
         R= (X(IX)-CENTEC(1))*(X(IX)-CENTEC(1))
     >    + (Y(IY)-CENTEC(2))*(Y(IY)-CENTEC(2))
         IF( R.LE.RAYONC(1) ) NBINT= NBINT+1
         IF( R.GE.RAYONC(2) ) NBEXT= NBEXT+1
   10 CONTINUE
      IF( NBINT.EQ.4 )THEN
C
C        RECTANGLE IS CONTAINED INSIDE THE INTERNAL RADIUS
         LELCRN=.FALSE.
      ELSEIF( NBEXT.EQ.4 )THEN
         IF( Y(1).LT.CENTEC(2).AND.CENTEC(2).LT.Y(2) )THEN
            IF( CENTEC(1).LT.X(1) )THEN
               LELCRN= (X(1)-CENTEC(1))*(X(1)-CENTEC(1)).LT.RAYONC(2)
            ELSEIF( X(2).LT.CENTEC(1) )THEN
               LELCRN= (X(2)-CENTEC(1))*(X(2)-CENTEC(1)).LT.RAYONC(2)
            ELSE
               LELCRN=.TRUE.
            ENDIF
         ELSEIF( X(1).LT.CENTEC(1).AND.CENTEC(1).LT.X(2) )THEN
            IF( CENTEC(2).LT.Y(1) )THEN
               LELCRN= (Y(1)-CENTEC(2))*(Y(1)-CENTEC(2)).LT.RAYONC(2)
            ELSEIF( Y(2).LT.CENTEC(2) )THEN
               LELCRN= (Y(2)-CENTEC(2))*(Y(2)-CENTEC(2)).LT.RAYONC(2)
            ELSE
               LELCRN=.TRUE.
            ENDIF
         ELSE
            LELCRN=.FALSE.
         ENDIF
      ELSE
         LELCRN=.TRUE.
      ENDIF
C
      RETURN
      END
