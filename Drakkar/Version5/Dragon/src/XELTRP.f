*DECK XELTRP
      SUBROUTINE XELTRP( IPGEOM,  NGIDL,   NDIM, NGEOME, L1CELL,
     >                   NTOTCO, NEXTGE,  MAXRO,   IPRT,  CELLG,
     >                    NSURO,  NVOLO, IDLDIM, IDLGEO, KEYTRN,
     >                    MAXDO,  MINDO, ICORDO, RMESHO, IDLREM,
     >                   INDEXO,  VOLSO, MATGEO)
************************************************************************
*                                                                      *
*           NAME: XELTRP                                               *
*      COMPONENT: EXCELL                                               *
*          LEVEL: 3 (CALLED BY 'XELTRK')                               *
*        VERSION: 1.0                                                  *
*       CREATION: 87/01                                                *
*       MODIFIED: 97/11 (G.M.) ELIMINATE CONHERENCE TESTS              *
*                 00/03 (R.R.) DECLARE ALL VARIABLE TYPES              *
*         AUTHOR: ROBERT ROY                                           *
*                                                                      *
*     SUBROUTINE: THIS ROUTINE WILL PREPARE TRACKING BY PRODUCING      *
*                 THE REQUIRED NUMBERING AND CALCULATE VOLUMES AND     *
*                 SURFACES.                                            *
*                                                                      *
*--------+-------------- V A R I A B L E S -------------+--+-----------*
*  NAME  /                  DESCRIPTION                 /IO/MOD(DIMENS)*
*--------+----------------------------------------------+--+-----------*
* IPGEOM / POINTER TO THE GEOMETRY (L_GEOM)             /I./INT        *
* NGIDL  / LENGHT OF GEOMETRIC NUMBERING.               /I./INT        *
* NDIM   / # OF DIMENSIONS (2 OR 3).                    /I./INT        *
* NGEOME / # OF GEOMETRIES.                             /I./INT        *
* L1CELL / TO INDICATE IF THERE IS JUST 1 CELL.         /I./LOG        *
* NEXTGE / RECTANGULAR(0)/CIRCULAR(1) BOUNDARY.         /I./INT        *
* NTOTCO / TOT NUMBER OF CYLINDERS IN ALL GEOMETRIES.   /I./INT        *
* MAXRO  / MAX NUMBER OF REAL MESH VALUES IN 'RMESHO'.  /I./INT        *
* IPRT   / INTERMEDIATE PRINTING LEVEL FOR OUTPUT.      /I./INT        *
* CELLG  / TO KEEP GEOMETY NAMES.                       /I./C*4(3NGEOM)*
* NSURO  / # OF SURFACES OF EACH GEOMETRY.              /I./INT(NGEOME)*
* NVOLO  / # OF ZONES OF EACH GEOMETRY.                 /I./INT(NGEOME)*
* IDLDIM / POSITION OF EACH GEOMETRY IN CYLINDER #ING.  /I./INT(NGEOME)*
* IDLGEO / POSITION OF EACH GEOMETRY IN THE             /I./INT(NGEOME)*
*        /            GEOMETRY NUMBERING SCHEME.        /  /           *
* KEYTRN / TURN # OF EACH GEOMETRY.                     /I./INT(NGEOME)*
* MAXDO  / MAX INDEX VALUES FOR ALL AXES (RECT/CYL).    /.O/INT(NTOTCO)*
* MINDO  / MIN INDEX VALUES FOR ALL AXES (RECT/CYL).    /.O/INT(NTOTCO)*
* ICORDO / PRINCIPAL AXES DIRECTION (X/Y/Z) FOR MESHES. /.O/INT(NTOTCO)*
* RMESHO / REAL MESH VALUES (RECT/CYL).                 /.O/REL(MAXRO) *
* IDLREM / POSITION OF MESH VALUES PER GEOMETRY.        /.O/INT(NGEOME)*
* INDEXO / INDEX FOR SEARCH IN 'RMESHO'.                /.O/I(4*NGIDL) *
* VOLSO  / VOLUMES & SURFACES FOR EACH GEOMETRY.        /.O/REL(NGIDL) *
* MATGEO / MATERIAL #S CORRESPONDING TO GEOMETRIES.     /.O/INT(NGIDL) *
*--------+---------------- R O U T I N E S -------------+--+-----------*
*  NAME  /                  DESCRIPTION                                *
*--------+-------------------------------------------------------------*
* KELRNG / TO READ THE USER NUMBERING OF ZONES.                        *
* KELMRG / TO READ THE USER MERGING OF ZONES.                          *
* KELSYM / TO ESTABLISH THE SURFACE SYMMETRIES.                        *
* XELGRD / TO READ GEOMETRIC INPUT.                                    *
* XELVOL / TO COMPUTE VOLUME IN 3D GEOMETRIES.                         *
************************************************************************
C
      USE               GANLIB
      IMPLICIT          NONE
C
      TYPE(C_PTR)       IPGEOM 
      INTEGER           NGIDL, NDIM, NGEOME, NTOTCO, NEXTGE, MAXRO, IPRT
      INTEGER           MAXDO(NTOTCO), MINDO(NTOTCO),   ICORDO(NTOTCO),
     >                  MATGEO(NGIDL), CELLG(3*NGEOME),
     >                  NSURO(NGEOME),  NVOLO(NGEOME), IDLDIM(NGEOME),
     >                  IDLGEO(NGEOME), IDLREM(NGEOME), KEYTRN(NGEOME),
     >                  INDEXO(4,NGIDL)
      REAL              RMESHO(MAXRO), VOLSO(NGIDL)
C
      INTEGER           NSTATE, IOUT, MAXTUR
      PARAMETER       ( NSTATE=40, IOUT=6, MAXTUR=12 )
      INTEGER           ISTATE(NSTATE)
      INTEGER           NTOTRM, NGEO, NTC, ITURN, NC, NCPC, NVSP1,
     >                  NO, NSYM, MAXC, KELRNG, KELMRG, KELSYM
      LOGICAL           L1CELL
      CHARACTER         CNAMEG*12, CTURN(2*MAXTUR)*2
C----
C  ALLOCATABLE ARRAYS
C----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: KEYSYM
C----
C  DATA STATEMENTS
C----
      DATA       CTURN / ' A',' B',' C',' D',' E',' F',' G',' H',
     >                   ' I',' J',' K',' L',
     >                   '-A','-B','-C','-D','-E','-F','-G','-H',
     >                   '-I','-J','-K','-L' /
C----
C  SCRATCH STORAGE ALLOCATION
C   KEYSYM: symmetry key giving the symmetric surface
C----
      ALLOCATE(KEYSYM(NGIDL))
C
C     LOOP OVER ALL GEOMETRIES
      NTOTRM= 0
      DO 90 NGEO= 1, NGEOME
         NTC= IDLDIM(NGEO)+1
         ITURN= KEYTRN(NGEO)
         WRITE( CNAMEG( 1: 4),'(A4)') CELLG(3*NGEO-2)
         WRITE( CNAMEG( 5: 8),'(A4)') CELLG(3*NGEO-1)
         WRITE( CNAMEG( 9:12),'(A4)') CELLG(3*NGEO  )
         IF( .NOT.L1CELL ) CALL LCMSIX(IPGEOM, CNAMEG, 1)
         CALL XDISET(ISTATE,NSTATE,0)
         CALL LCMGET(IPGEOM, 'STATE-VECTOR', ISTATE)
         IF( ISTATE(1).GE.20.OR.ISTATE(1).EQ.3.OR.ISTATE(1).EQ.6 )THEN
            NC= 1
         ELSE
            NC= 0
         ENDIF
         IF( IPRT.GT.1 )THEN
            WRITE(IOUT,'(1H )')
            IF    ( NC.EQ.0 )THEN
              WRITE(IOUT,'(/27H NUMBERING PHYSICAL CELL # ,I8/6H  >>> ,
     >                     A12,6H /ROT ,A2,13H GEOMETRY <<<,
     >                     13H    (WITH  NO,11H CYLINDER ) )')
     >                                 NGEO,        CNAMEG,CTURN(ITURN)
            ELSEIF( NC.EQ.1 )THEN
              WRITE(IOUT,'(/27H NUMBERING PHYSICAL CELL # ,I8/6H  >>> ,
     >                     A12,6H /ROT ,A2,13H GEOMETRY <<<,
     >                     13H    (WITH ONE,11H CYLINDER ) )')
     >                                 NGEO,        CNAMEG,CTURN(ITURN)
            ELSE
              WRITE(IOUT,'(/27H NUMBERING PHYSICAL CELL # ,I8/6H  >>> ,
     >                     A12,6H /ROT ,A2,13H GEOMETRY <<<,
     >                     10H    (WITH ,I3,11H CYLINDERS) )')
     >                               NGEO, CNAMEG, CTURN(ITURN), NC
            ENDIF
         ENDIF
         NCPC  = NC + 3
         NVSP1 = NVOLO(NGEO) - NSURO(NGEO) + 1
C
C        LOOKING TO THE GEOMETRY
         CALL XELGRD( IPGEOM, IPRT, NDIM, NEXTGE, ITURN,
     >                MAXC, RMESHO(NTOTRM+1),
     >                MINDO(NTC), MAXDO(NTC), ICORDO(NTC))
C
C        RENUMBER
         NO=   KELRNG(IPRT, NDIM, NEXTGE, NCPC,
     >                MINDO(NTC), MAXDO(NTC), ICORDO(NTC),
     >                NSURO(NGEO), NVOLO(NGEO), IDLGEO(NGEO),
     >                MAXC, RMESHO(NTOTRM+1), MATGEO, VOLSO, INDEXO)
C
C        MERGE
         NO= KELMRG(IPGEOM,NSURO(NGEO),NVOLO(NGEO),IDLGEO(NGEO),MATGEO)
         IF( NO.NE.NVSP1 )THEN
            IF( IPRT.GT.1 )THEN
               WRITE(IOUT,'(1H )')
               WRITE(IOUT,'(22H     MERGE INTO   >>> ,I8,
     >                  13H  ZONES   <<<)')
     >                       NO+NSURO(NGEO)-1
            ENDIF
         ENDIF
C
C        ESTABLISH NECESSARY SYMMETRIES
         NSYM= KELSYM( IPRT, NDIM, MAXDO(NTC), NSURO(NGEO), NVOLO(NGEO),
     >                 IDLGEO(NGEO), INDEXO, MATGEO,KEYSYM)
C
C        COMPUTE VOLUMES
         CALL XELVOL( IPRT, NDIM, NEXTGE, NCPC,
     >                MINDO(NTC), MAXDO(NTC), ICORDO(NTC),
     >                NSURO(NGEO), NVOLO(NGEO), IDLGEO(NGEO),INDEXO,
     >                MAXC, RMESHO(NTOTRM+1), MATGEO, VOLSO )
         IDLREM(NGEO)= NTOTRM
         NTOTRM= NTOTRM + MAXC
         IF( .NOT.L1CELL ) CALL LCMSIX(IPGEOM, ' ', 2 )
   90 CONTINUE
      IF( NTOTRM.GT.MAXRO )THEN
         CALL XABORT( 'XELTRP : INCREASE MAXREM => SEE DEVELOPPER')
      ENDIF
C----
C  SCRATCH STORAGE DEALLOCATION
C----
      DEALLOCATE(KEYSYM)
C
      RETURN
      END
