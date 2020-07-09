*DECK XELPRC
      SUBROUTINE XELPRC (IPGEOM,GEONAM,NDIM,NNCYL,NNSUR,NNVOL,NAXREM)
************************************************************************
*                                                                      *
*           NAME: XELPRC                                               *
*      COMPONENT: EXCELL                                               *
*          LEVEL: 4 (CALLED BY 'XELDCL')                               *
*        VERSION: 1.0                                                  *
*       CREATION: 89/12                                                *
*       MODIFIED: 97/11 (G.M.) ELIMINATE SPLIT>2 ABORT                 *
*                 00/03 (R.R.) DECLARE ALL VARIABLE TYPES              *
*         AUTHOR: ROBERT ROY                                           *
*                                                                      *
*     SUBROUTINE: THIS ROUTINE READS A CELL GEOMETRY ON LCM AND        *
*                 CHECK IF THE GEOMETRY IS ACCEPTABLE FOR "EXCELL".    *
*                                                                      *
*--------+-------------- V A R I A B L E S -------------+--+-----------*
*  NAME  /                  DESCRIPTION                 /IO/MOD(DIMENS)*
*--------+----------------------------------------------+--+-----------*
* IPGEOM / POINTER TO THE GEOMETRY (L_GEOM)             /I./INT        *
* GEONAM / GEOMETRY NAME                                /I./CAR*12     *
* NDIM   / # OF DIMENSIONS ( 2 OR 3)                    /I./INT        *
* NNCYL  / # OF CYLINDERS IN THE GEOMETRY               /.O/INT        *
* NNSUR  / # OF SURFACES                                /.O/INT        *
* NNVOL  / # OF VOLUMES                                 /.O/INT        *
* NAXREM / MAX # OF COORDINATES TO SPECIFY THAT CELL    /.O/INT        *
************************************************************************
C
      USE          GANLIB
      IMPLICIT     NONE
C
C     DECLARE      DUMMY ARGUMENTS
      TYPE(C_PTR)  IPGEOM 
      INTEGER      NDIM, NNCYL, NNSUR, NNVOL, NAXREM
      CHARACTER*12 GEONAM
C
C     DECLARE      LOCAL VARIABLES
      INTEGER      NLCM, NIXS, NIST, NSTATE, MAXSPL
      PARAMETER  ( NLCM=26, NIXS=11, NIST=2, NSTATE=40, MAXSPL=100 )
      CHARACTER*12 LCMNM(NLCM)
      INTEGER      LNLCM(NLCM),INVLCM(NIXS),INVSTA(NIST),
     >             ISTATE(NSTATE),ISPLT(MAXSPL)
      INTEGER      ILCM, IIXS, IIST, ITYPE, LR, LX, LY, LZ, ISPLIT,
     >             JX, JY, JZ, JR, JL, ILEN, ITYLCM
C
      DATA INVLCM/  6, 11, 12, 14,        16, 17, 18, 19,
     >             20, 21, 22 /
      DATA INVSTA/ 8, 12 /
      DATA LCMNM /  'MIX',  'MESHX',  'MESHY',   'MESHZ',  'RADIUS',
     >             'SIDE', 'SPLITX', 'SPLITY',  'SPLITZ',  'SPLITR',
     >             'CELL',  'COORD',  'MERGE',    'TURN', 'CLUSTER',
     >             'NPIN',   'RPIN',   'APIN',   'BIHET',  'POURCE',
     >           'PROCEL',   'IHEX',  'NCODE',   'ZCODE',   'ICODE',
     >           'CENTER'/
C
      DO 10 ILCM= 1, NLCM
         CALL LCMLEN(IPGEOM,LCMNM(ILCM),LNLCM(ILCM),ITYLCM)
   10 CONTINUE
C
C     ELIMINATES THE INVALID OPTIONS
      DO 20 IIXS= 1, NIXS
        IF( LNLCM(INVLCM(IIXS)).NE.0 )
     >     CALL XABORT( 'XELPRC:*'//GEONAM//'* IS '//
     >                  'NOT A VALID CELL GEOMETRY FOR EXCELL'//
     >                  ' (LCM BLOCK *'//LCMNM(INVLCM(IIXS))//'*)')
   20 CONTINUE
      CALL LCMLEN(IPGEOM,'STATE-VECTOR',ILEN,ITYLCM)
      IF(ILEN .LT. 1 .OR. ILEN .GT. NSTATE )
     >   CALL XABORT( 'XELPRC: GEOMETRY HAS INVALID STATE VECTOR')
      CALL LCMGET(IPGEOM,'STATE-VECTOR',ISTATE)
      DO 30 IIST= 1, NIST
        IF( ISTATE(INVSTA(IIST)).NE.0 )
     >     CALL XABORT( 'XELPRC: INVALID GEOMETRY FOR EXCELL')
   30 CONTINUE
C
      ITYPE=  ISTATE(1)
      LR=     ISTATE(2)
      LX=     MAX(1,ISTATE(3))
      LY=     MAX(1,ISTATE(4))
      LZ=     MAX(1,ISTATE(5))
      NNVOL=  ISTATE(6)
      ISPLIT= ISTATE(11)
C
C     GET THE SPLITTING INFORMATION, AND COMPUTE JR, JX, JY, JZ VALUES
      IF( ISPLIT.GT.0 )THEN
         JR= 0
         JX= 0
         JY= 0
         JZ= 0
         CALL LCMLEN(IPGEOM,'SPLITR',ILEN,ITYLCM)
         IF( ILEN.GT.MAXSPL )THEN
            CALL XABORT('XELPRC: SPLITR OVERFLOW')
         ELSEIF( ILEN.EQ.0 )THEN
            JR= LR
         ELSEIF( ILEN.NE.LR )THEN
            CALL XABORT( 'XELPRC: R-SPLITTING NOT ACCEPTED' )
         ELSE
            CALL LCMGET(IPGEOM,'SPLITR',ISPLT)
            JR=   0
            DO 15 JL= 1, ILEN
               JR= JR + ABS(ISPLT(JL))
   15       CONTINUE
         ENDIF
         CALL LCMLEN(IPGEOM,'SPLITX',ILEN,ITYLCM)
         IF( ILEN.GT.MAXSPL )THEN
            CALL XABORT('XELPRC: SPLITX OVERFLOW')
         ELSEIF( ILEN.EQ.0 )THEN
            JX= LX
         ELSEIF( ILEN.NE.LX )THEN
            CALL XABORT( 'XELPRC: X-SPLITTING NOT ACCEPTED' )
         ELSE
            CALL LCMGET(IPGEOM,'SPLITX',ISPLT)
            JX=   0
            DO 25 JL= 1, ILEN
               JX= JX + ISPLT(JL)
   25       CONTINUE
         ENDIF
         CALL LCMLEN(IPGEOM,'SPLITY',ILEN,ITYLCM)
         IF( ILEN.GT.MAXSPL )THEN
            CALL XABORT('XELPRC: SPLITY OVERFLOW')
         ELSEIF( ILEN.EQ.0 )THEN
            JY= LY
         ELSEIF( ILEN.NE.LY )THEN
            CALL XABORT( 'XELPRC: Y-SPLITTING NOT ACCEPTED' )
         ELSE
            CALL LCMGET(IPGEOM,'SPLITY',ISPLT)
            JY=   0
            DO 35 JL= 1, ILEN
               JY= JY + ISPLT(JL)
   35       CONTINUE
         ENDIF
         CALL LCMLEN(IPGEOM,'SPLITZ',ILEN,ITYLCM)
         IF(ILEN.GT.MAXSPL) CALL XABORT('XELPRC: SPLITZ OVERFLOW')
         IF( ILEN.EQ.0 )THEN
            JZ= LZ
         ELSEIF( ILEN.NE.LZ )THEN
            CALL XABORT( 'XELPRC: Z-SPLITTING NOT ACCEPTED' )
         ELSE
            JZ= 0
            CALL LCMGET(IPGEOM,'SPLITZ',ISPLT)
            DO 45 JL= 1, ILEN
               JZ= JZ + ISPLT(JL)
   45       CONTINUE
         ENDIF
      ELSE
         JR= LR
         JX= LX
         JY= LY
         JZ= LZ
      ENDIF
C
      IF( ITYPE.EQ.0 )THEN
C
C        VIRTUAL ELEMENT
         NNVOL= 0
         NNCYL= 0
         NNSUR= 0
         NAXREM= 0
      ELSE
         IF( NDIM.EQ.2 )THEN
            NNSUR= 2 * (JX+JY)
            NNVOL=  JX*JY
            IF( ITYPE.EQ.5 )THEN
C              FOR *CAR2D* GEOMETRY
C
               NNCYL=  0
C
C              X-AXIS:JX+1, Y-AXIS:JY+1, Z-AXIS:2
               NAXREM= JX+JY+4
            ELSEIF( ITYPE.EQ.3 )THEN
C              FOR *TUBE* GEOMETRY
C
               NNCYL=  1
               IF( JX.NE.1 .OR. JY.NE.1 )THEN
                  CALL XABORT( 'XELPRC: FOR TUBE, PLEASE NO XY SPLIT')
               ENDIF
               NNVOL=  NNVOL+JX*JY*JR
C
C              X-AXIS:JX+1, Y-AXIS:JY+1, Z-AXIS:2, R-AXIS:JR+3
               NAXREM= JX+JY+JR+7
            ELSEIF( ITYPE.EQ.20 )THEN
C              FOR *CARCEL* GEOMETRY
C
               NNCYL=  1
               NNVOL=  NNVOL+JX*JY*JR
C
C              X-AXIS:JX+1, Y-AXIS:JY+1, Z-AXIS:2, R-AXIS:JR+3
               NAXREM= JX+JY+JR+7
            ELSE
               CALL XABORT('XELPRC: INVALID CELL GEOMETRY FOR EXCELL=>'
     >                     //GEONAM(1:12) )
            ENDIF
         ELSE
            NNSUR=  2 * (JX*JY+JX*JZ+JY*JZ )
            NNVOL=  JX*JY*JZ
            IF( ITYPE.EQ.7 )THEN
C              FOR *CAR3D* GEOMETRY
C
               NNCYL=  0
C
C              X-AXIS:JX+1, Y-AXIS:JY+1, Z-AXIS:JZ+1
               NAXREM= JX+JY+JZ+3
            ELSEIF( ITYPE.EQ. 6 .OR. ITYPE.EQ.21 .OR.
     >              ITYPE.EQ.22 .OR. ITYPE.EQ.23 )THEN
C              FOR *TUBEZ*, *CARCELX*, *CARCELY* OR *CARCELZ* GEOMETRY
C
              NNCYL= 1
              IF( ITYPE.EQ.6 )THEN
                 IF( JX.NE.1 .OR. JY.NE.1 ) THEN
                    CALL XABORT('XELPRC: FOR TUBEZ, PLEASE NO XY SPLIT')
                 ENDIF
              ELSEIF( ITYPE.EQ.23 )THEN
                 NNSUR= NNSUR+2*JR*JX*JY
              ELSEIF( ITYPE.EQ.22 )THEN
                 NNSUR= NNSUR+2*JR*JX*JZ
              ELSEIF( ITYPE.EQ.21 )THEN
                 NNSUR= NNSUR+2*JR*JY*JZ
              ENDIF
              NNVOL= NNVOL+JR*JX*JY*JZ
C
C             X-AXIS:JX+1, Y-AXIS:JY+1, Z-AXIS:JZ+1, R-AXIS:JR+3
              NAXREM= JX+JY+JZ+JR+6
            ELSE
            CALL XABORT( 'XELPRC: INVALID CELL GEOMETRY FOR EXCELL=>'//
     >                   GEONAM(1:12) )
            ENDIF
         ENDIF
      ENDIF
C
      RETURN
      END
