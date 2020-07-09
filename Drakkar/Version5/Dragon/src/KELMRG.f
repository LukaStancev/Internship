*DECK KELMRG
      FUNCTION    KELMRG(IPGEOM, NSURO, NVOLO, IDLGEO, MATGEO)
************************************************************************
*                                                                      *
*           NAME: KELMRG                                               *
*      COMPONENT: EXCELL                                               *
*        VERSION: 1.0                                                  *
*       CREATION: 90/04                                                *
*       MODIFIED: 00/03 (R.R.) DECLARE ALL VARIABLE TYPES              *
*         AUTHOR: ROBERT ROY                                           *
*                                                                      *
*       FUNCTION: THIS FUNCTION WILL MERGE ZONES                       *
*                 FOR A HETEROGENEOUS BLOCK                            *
*                                                                      *
*          LEVEL: 4 (CALLED BY 'XELTRP')                               *
*                                                                      *
*--------+-------------- V A R I A B L E S -------------+--+-----------*
*  NAME  /                  DESCRIPTION                 /IO/MOD(DIMENS)*
*--------+----------------------------------------------+--+-----------*
* KELMRG / FUNCTION = # OF SURFACES & ZONES RENUMBERED. /.O/INT        *
* IPGEOM / POINTER TO THE GEOMETRY (L_GEOM)             /I./INT        *
* NSURO  / # OF SURFACES FOR A SPECIFIC GEOMETRY.       /I./INT        *
* NVOLO  / # OF ZONES FOR A SPECIFIC GEOMETRY.          /I./INT        *
* IDLGEO / SPECIFIC POSITION FOR A GEOMETRY.            /I./INT        *
* MATGEO / #ING OF ZONES AND SURFACES FOR ALL GEOMETRIES/.O/INT(*)     *
************************************************************************
      USE                 GANLIB
      IMPLICIT            NONE
      LOGICAL             SWONCE
      TYPE(C_PTR)         IPGEOM
      INTEGER             KELMRG,NSURO,NVOLO,IDLGEO,MATGEO(*)
      INTEGER             IOUT, IND, MATMIN, MATMAX, ITYLCM, IMRG, JMRG,
     >                    ILEN, I
      PARAMETER         ( IOUT=6 )
C
      IND(I)= IDLGEO + I
      CALL LCMLEN(IPGEOM, 'MERGE', ILEN, ITYLCM)
      IF( ILEN.EQ.0 )THEN
         KELMRG= NVOLO - NSURO + 1
      ELSE
         IF( ILEN.GT.NVOLO )
     >      CALL XABORT('KELMRG: MERGING HAS TOO MANY ZONES' )
         CALL LCMGET(IPGEOM, 'MERGE', MATGEO(IND(1)) )
         MATMIN= 100000000
         MATMAX=-100000000
         DO 10 IMRG= 1, ILEN
            IF( MATGEO(IND(IMRG)).LT.MATMIN) MATMIN= MATGEO(IND(IMRG))
            IF( MATGEO(IND(IMRG)).GT.MATMAX) MATMAX= MATGEO(IND(IMRG))
   10    CONTINUE
         IF( MATMIN.NE.1 )
     >      CALL XABORT('KELMRG: NO FIRST MERGING ZONE' )
         DO 30 JMRG= MATMIN, MATMAX
            SWONCE= .FALSE.
            DO 20 IMRG= 1, ILEN
               SWONCE= SWONCE.OR.(MATGEO(IND(IMRG)).EQ.JMRG)
   20       CONTINUE
            IF( .NOT.SWONCE )THEN
               WRITE(IOUT,*) 'WHERE IS MERGE REGION NO.', JMRG
               CALL XABORT('KELMRG: ERROR IN MERGE NUMBERING' )
            ENDIF
   30    CONTINUE
         KELMRG= MATMAX - NSURO + 1
      ENDIF
C
      RETURN
      END
