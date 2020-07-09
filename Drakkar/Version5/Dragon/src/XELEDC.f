*DECK XELEDC
      SUBROUTINE XELEDC(   NDIM, MAXGRI, NGEOME, NTOTCO, NTYPES,
     >                   NBLOCK,  NUNKO,
     >                    NSURO,  NVOLO,  MINDO,  MAXDO,
     >                   ICORDO, IDLDIM, KEYGEO,
     >                   KEYTYP, IDLBLK, KEYINT,
     >                   NTOTCL,   MAXR,   NSUR,   NVOL, KEYCYL )
************************************************************************
*                                                                      *
*           NAME: XELEDC                                               *
*      COMPONENT: EXCELL                                               *
*          LEVEL: 3 (CALLED BY 'XELTRK')                               *
*        VERSION: 1.0                                                  *
*       CREATION: 90/08                                                *
*       MODIFIED: 00/03 (R.R.) DECLARE ALL VARIABLE TYPES              *
*         AUTHOR: ROBERT ROY                                           *
*                                                                      *
*     SUBROUTINE: THIS ROUTINE WILL ASSOCIATE ALL BLOCKS OF A PROBLEM  *
*                 TO ONLY ONE GEOMETRY & WILL GENERATE THE 4           *
*                 USEFUL INTEGER VALUES THAT WILL DESCRIBE THE PROBLEM *
*                 IN ITS EXACT GEOMETRIC DESCRIPTION.                  *
*                                                                      *
*--------+-------------- V A R I A B L E S -------------+--+-----------*
*  NAME  /                  DESCRIPTION                 /IO/MOD(DIMENS)*
*--------+----------------------------------------------+--+-----------*
* NDIM   / # OF DIMENSIONS.                             /I./INT        *
* MAXGRI / # OF GRID CELL IN X/Y/Z DIRECTIONS           /I./INT(3)     *
* NGEOME / # OF GEOMETRIES.                             /I./INT        *
* NTOTCO / TOT # OF CYLINDERS IN ALL GEOMETRIES         /I./INT        *
* NTYPES / # OF TYPES.                                  /I./INT        *
* NBLOCK / # OF BLOCKS.                                 /I./INT        *
* NUNKO  / # OF UNKNOWNS.                               /I./INT        *
* NSURO  / # OF SURFACES OF EACH GEOMETRY.              /I./INT(NTYPES)*
* NVOLO  / # OF ZONES OF EACH GEOMETRY.                 /I./INT(NTYPES)*
* MINDO  / MIN INDEX IN THE REMESH ARRAY.               /I./INT(NTOTCO)*
* MAXDO  / MIN INDEX IN THE REMESH ARRAY.               /I./INT(NTOTCO)*
* ICORDO / COORDINATE   FOR REMESH ARRAY.               /I./INT(NTOTCO)*
* IDLDIM / POSITION OF EACH GEOEMTRY IN CYLINDERS #ING. /I./INT(NGEOME)*
* KEYGEO / GEOMETRIC KEY FOR EACH TYPE.                 /I./INT(NTYPES)*
* KEYTYP / TYPE KEY FOR EACH BLOCK.                     /I./INT(NBLOCK)*
* IDLBLK / POSITION OF EACH BLOCK IN NUMBERING SCHEME.  /I./INT(NBLOCK)*
* KEYINT / #ING OF CELL INTERFACES.                     /I./INT(NUNKO )*
* NTOTCL / TOT # OF CYLINDERS IN EXACT GEOMETRY.        /.O/INT        *
* MAXR   / LENGHT TO STOCK REAL ABSCISSAE.              /.O/INT        *
* NSUR   / # OF SURFACES OF EXACT GEOMETRY (NEGATIVE).  /.O/INT        *
* NVOL   / # OF ZONES OF EXACT GEOMETRY.                /.O/INT        *
* KEYCYL / INDEX OF CYLINDERS BY BLOCK.                 /.O/INT(NBLOCK)*
************************************************************************
C
      IMPLICIT     NONE
C
      INTEGER      NDIM, NGEOME, NTOTCO, NTYPES, NBLOCK, NUNKO,
     >             NTOTCL, MAXR, NSUR, NVOL
      INTEGER      MAXGRI(3),      NSURO(NTYPES),  NVOLO(NTYPES),
     >             MINDO(NTOTCO),  MAXDO(NTOTCO),  ICORDO(NTOTCO),
     >             IDLDIM(NTYPES), KEYGEO(NTYPES),
     >             KEYTYP(NBLOCK), IDLBLK(NBLOCK), KEYCYL(NBLOCK),
     >             KEYINT( NUNKO)
C
      INTEGER      ICUR(3), IBLK, N, ICX, ITYP, IGEO, IDLD, MDMIN,
     >             NP1, NP2, IP1, IP2, IP3, NC, NSUX, NVOX, IVX
      INTEGER      NUMBLK, I, K
C
      NUMBLK(I,K)= I + IDLBLK(K)
C
      DO 5 IBLK= 1, NBLOCK
         KEYCYL(IBLK)= 0
    5 CONTINUE
C
C     DETERMINE: NTOTCL & MAXR
C.1)  RECONSTRUCT CARTESIAN MESH
      MAXR= 0
      NTOTCL= 3
      ICUR(1)= 1
      ICUR(2)= 1
      ICUR(3)= 1
      DO 30 N= 1, 3
C
C        SCANNING CELLS ON THE AXIS #N
         DO 20 ICX= 1, MAXGRI(N)
            ICUR(N)= ICX
            IF( NDIM.EQ.2 )THEN
               IBLK= MAXGRI(1) * (ICUR(2) - 1) + ICUR(1)
            ELSE
               IBLK= MAXGRI(1)*(MAXGRI(2)*ICUR(3)+ICUR(2)-
     >                          MAXGRI(2))+ICUR(1)-MAXGRI(1)
            ENDIF
            ITYP= KEYTYP(IBLK)
            IF( ITYP.EQ.0 ) GO TO 20
            IGEO= KEYGEO(ITYP)
            IDLD= IDLDIM(IGEO)
            MAXR= MAXR + (MAXDO(IDLD+N)-MINDO(IDLD+N))
   20    CONTINUE
         ICUR(N)= 1
         MAXR= MAXR+1
   30 CONTINUE
C
C.2)  RECONSTRUCT INFORMATIONS FOR CYLINDRICAL MESH
      IF( NDIM.EQ.2 )THEN
         MDMIN= 3
      ELSE
         MDMIN= 1
      ENDIF
      DO 130 N= MDMIN, 3
         ICUR(N)= 1
         NP1= MOD(N  ,3) + 1
         NP2= MOD(N+1,3) + 1
         DO 120 IP2= 1, MAXGRI(NP2)
         DO 110 IP1= 1, MAXGRI(NP1)
            ICUR(NP1)= IP1
            ICUR(NP2)= IP2
            IF( NDIM.EQ.2 )THEN
               IBLK= MAXGRI(1) * (ICUR(2) - 1) + ICUR(1)
            ELSE
               IBLK= MAXGRI(1)*(MAXGRI(2)*ICUR(3)+ICUR(2)-
     >                          MAXGRI(2))+ICUR(1)-MAXGRI(1)
            ENDIF
            ITYP= KEYTYP(IBLK)
            IF( ITYP.EQ.0 ) GO TO 110
            IGEO= KEYGEO(ITYP)
            IDLD= IDLDIM(IGEO)
            IF( IGEO.NE.NGEOME )THEN
               NC= IDLDIM(IGEO+1)-IDLD-3
            ELSE
               NC= NTOTCO-IDLD-3
            ENDIF
            IF( NC.EQ.1 )THEN
               IF( ICORDO(IDLD+4).EQ.N )THEN
                  NTOTCL= NTOTCL+1
                  MAXR= MAXR + 3 + (MAXDO(IDLD+4)-MINDO(IDLD+4))
                  DO 105 IP3= 1, MAXGRI(N)
                     ICUR(N)= IP3
                     IF( NDIM.EQ.2 )THEN
                        IBLK= MAXGRI(1) * (ICUR(2) - 1) + ICUR(1)
                     ELSE
                        IBLK= MAXGRI(1)*(MAXGRI(2)*ICUR(3)+ICUR(2)-
     >                                   MAXGRI(2))+ICUR(1)-MAXGRI(1)
                     ENDIF
                     KEYCYL(IBLK)= NTOTCL
  105             CONTINUE
                  ICUR(N)= 1
               ENDIF
            ENDIF
  110    CONTINUE
  120    CONTINUE
  130 CONTINUE
C
C     DETERMINE: NSUR & NVOL
      NSUR= 0
      NVOL= 0
      DO 230 IBLK= 1,NBLOCK
         ITYP= KEYTYP(IBLK)
         IF( ITYP.EQ.0 ) THEN
            CALL XABORT( '*** XELEDC: EXACT VOID CELL NOT ALLOWED')
         ENDIF
         IGEO= KEYGEO(ITYP)
         NSUX= NSURO(IGEO)
         NVOX= NVOLO(IGEO)
         DO 220 IVX= NSUX, NVOX
            IF( IVX.LT.0 )THEN
               IF( KEYINT(NUMBLK(IVX,IBLK)).EQ.0 ) NSUR= NSUR-1
            ELSEIF( IVX.GT.0 )THEN
               NVOL= NVOL + 1
            ENDIF
  220    CONTINUE
  230 CONTINUE
C
      RETURN
      END
