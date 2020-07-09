*DECK XELETR
      SUBROUTINE XELETR(   IPRT,   NDIM, MAXGRI, NGEOME, NTOTCO, NTYPES,
     >                    NTIDL, NBLOCK,   NSUR,   NVOL, NTOTCL,  NUNKO,
     >                    NSURO,  NVOLO,  MINDO,  MAXDO, ICORDO, IDLDIM,
     >                   IDLGEO, KEYGEO, IDLTYP, KEYTYP, IDLBLK, KEYCYL,
     >                   RMESHO, IDLREM, INDEXO,  VOLSO, MATGEO, KEYINT,
     >                   MATTYP, REMESH, MINDIM, MAXDIM,  ICORD, VOLSUR,
     >                   KEYMRG,  INDEX, INCELL, MATALB,  NSURC,  NVOLC)
************************************************************************
*                                                                      *
*           NAME: XELETR                                               *
*      COMPONENT: EXCELL                                               *
*          LEVEL: 3 (CALLED BY 'XELTRK')                               *
*        VERSION: 1.0                                                  *
*       CREATION: 90/08                                                *
*       MODIFIED: 97/11 (G.M.) INTRODUCE MATGEO()=0 CASE               *
*                 00/03 (R.R.) DECLARE ALL VARIABLE TYPES              *
*         AUTHOR: ROBERT ROY                                           *
*                                                                      *
*     SUBROUTINE: THIS ROUTINE WILL PREPARE TRACKING BY PRODUCING      *
*                 THE REQUIRED NUMBERING AND RECALCULATE MESH FOR      *
*                 AN EXACT GEOMETRY TREATMENT.                         *
*                                                                      *
*--------+-------------- V A R I A B L E S -------------+--+-----------*
*  NAME  /                  DESCRIPTION                 /IO/MOD(DIMENS)*
*--------+----------------------------------------------+--+-----------*
* IPRT   / INTERMEDIATE PRINTING LEVEL FOR OUTPUT.      /I./INT        *
* NDIM   / # OF DIMENSIONS (2 OR 3).                    /I./INT        *
* MAXGRI / # OF BLOCKS IN X/Y/Z DIRECTIONS.             /I./INT(3)     *
* NGEOME / # OF GEOMETRIES.                             /I./INT        *
* NTOTCO / TOT NUMBER OF CYLINDERS IN ALL GEOMETRIES.   /I./INT        *
* NTYPES / # OF CELL TYPES.                          .  /I./INT        *
* NTIDL  / LENGHT OF TYPE NUMBERING.                    /I./INT        *
* NBLOCK / # OF BLOCKS.                                 /I./INT        *
* NSUR   / # OF SURFACES.                               /I./INT        *
* NVOL   / # OF ZONES.                                  /I./INT        *
* NTOTCL / TOT # OF CYLINDERS IN EXACT GEOMETRY.        /I./INT        *
* NUNKO  / OLD # OF UNKNOWNS.                           /I./INT        *
* NSURO  / # OF SURFACES OF EACH GEOMETRY.              /I./INT(NGEOME)*
* NVOLO  / # OF ZONES OF EACH GEOMETRY.                 /I./INT(NGEOME)*
* MINDO  / MIN INDEX VALUES FOR ALL AXES (RECT/CYL).    /I./INT(NTOTCO)*
* MAXDO  / MAX INDEX VALUES FOR ALL AXES (RECT/CYL).    /I./INT(NTOTCO)*
* ICORDO / PRINCIPAL AXES DIRECTION (X/Y/Z) FOR MESHES. /I./INT(NTOTCO)*
* IDLDIM / POSITION OF EACH GEOMETRY IN CYLINDER #ING.  /I./INT(NGEOME)*
* IDLGEO / POSITION OF EACH GEOMETRY IN THE             /I./INT(NGEOME)*
*        /            GEOMETRY NUMBERING SCHEME.        /  /           *
* KEYGEO / GEOMETRIC KEY FOR EACH TYPE.                 /I./INT(NTYPES)*
* IDLTYP / POSITION OF EACH TYPE IN NUMBERING SCHEME.   /I./INT(NTYPES)*
* KEYTYP / TYPE KEY FOR EACH BLOCK.                     /I./INT(NBLOCK)*
* IDLBLK / POSITION OF EACH BLOCK IN NUMBERING SCHEME.  /I./INT(NBLOCK)*
* KEYCYL / INDEX OF CYLINDERS BY BLOCK.                 /I./INT(NBLOCK)*
* RMESHO / REAL MESH VALUES (RECT/CYL).                 /I./REL(*     )*
* IDLREM / POSITION OF MESH VALUES PER GEOMETRY.        /I./INT(NGEOME)*
* INDEXO / INDEX FOR SEARCH IN 'RMESHO'.                /I./I(4*NGIDL) *
* VOLSO  / VOLUMES & SURFACES FOR EACH GEOMETRY.        /I./REL(NGIDL) *
* MATGEO / MATERIAL #S CORRESPONDING TO GEOMETRIES.     /I./INT(NGIDL) *
* KEYINT / INTERFACE KEY (GIVING THE CONNECTED SURFACE)./I./INT(NUNKO )*
* MATTYP / MATERIAL #S FOR ZONES OF EVERY TYPE.         /I./INT(NTIDL) *
* REMESH / REAL MESH VALUES (RECT/CYL).                 /.O/REL(*     )*
* MINDIM / MIN INDEX VALUES FOR ALL AXES (RECT/CYL).    /.O/INT(NTOTCL)*
* MAXDIM / MAX INDEX VALUES FOR ALL AXES (RECT/CYL).    /.O/INT(NTOTCL)*
* ICORD  / PRINCIPAL AXES DIRECTION (X/Y/Z) FOR MESHES. /.O/INT(NTOTCL)*
* VOLSUR / VOLUME-SURFACE VECTOR OF EXACT GEOMETRY.     /.O/REL(*    ) *
* KEYMRG / MERGING VECTOR        OF EXACT GEOMETRY.     /.O/INT(*    ) *
* INDEX  / #ING OF SURFACES & ZONES.                    /.O/INT(4,*  ) *
* INCELL / BLOCK   #ING.                                /.O/INT(*)     *
* MATALB / MATERIAL TYPES.                              /.O/INT(*)     *
* NSURC  / # OF COMPRESSED SURFACES.                    /.O/INT        *
* NVOLC  / # OF COMPRESSED ZONES.                       /.O/INT        *
************************************************************************
C
      IMPLICIT           NONE
C
      INTEGER              IPRT,   NDIM, NGEOME, NTOTCO, NTYPES,
     >                    NTIDL, NBLOCK,   NSUR,   NVOL, NTOTCL,  NUNKO,
     >                    NSURC,  NVOLC
      INTEGER            MAXGRI(3),
     >                    MAXDO(NTOTCO),  MINDO(NTOTCO), ICORDO(NTOTCO),
     >                    NSURO(NGEOME),  NVOLO(NGEOME), IDLDIM(NGEOME),
     >                   IDLGEO(NGEOME), IDLREM(NGEOME), KEYGEO(NTYPES),
     >                   IDLTYP(NTYPES),
     >                   KEYTYP(NBLOCK), IDLBLK(NBLOCK), KEYCYL(NBLOCK),
     >                   INDEXO(4,*), MATGEO(*),
     >                   KEYINT(NUNKO), MATTYP(NTIDL),
     >                   MINDIM(NTOTCL), MAXDIM(NTOTCL), ICORD(NTOTCL),
     >                   INDEX(4,*), KEYMRG(*), MATALB(*), INCELL(*)
      REAL               RMESHO(*), REMESH(*), VOLSO(*), VOLSUR(*)
C
      INTEGER            ICUR(4)
      INTEGER            NUNK, IDLGE2, IG2, I4, N, ICX, IREM, I, J, K,
     >                   IBLK, ITYP, IGEO, IDLD, IDLR, MINABS, MAXABS,
     >                   J1, NTOTCX, MDMIN, NP1, NP2, IP1, IP2, IP3,
     >                   IOLD, ISU2, IVO2, ICREM, MINP1, MINP2, MINC,
     >                   ICYL, IDLTYX, IDLGEX, NO, NC, IMYG, IVSN,
     >                   IVS, IKREM
      REAL               RSTART, RMINUS, XP1, XP2
      CHARACTER          TEMESH(4)*8
      INTEGER            IOUT
      PARAMETER        ( IOUT=6 )
      INTEGER            NUMBLK, KL
      DATA        TEMESH / 'X', 'Y', 'Z', 'C' /
C
      NUMBLK(I,K)= I + IDLBLK(K)
C
C     INITIALIZE: NO INTERFACE & PUT INDEXES TO 0.
      NUNK = NVOL + 1 - NSUR
      IDLGE2=       1 - NSUR
      DO 5 IG2= 1, NUNK
         VOLSUR(IG2)= 0.0
         KEYMRG(IG2)= 0
         MATALB(IG2)= 0
         INCELL(IG2)= 0
      DO 5 I4= 1,4
         INDEX(I4,IG2)= 0
    5 CONTINUE
C
      IF( IPRT.GE.1 )THEN
          WRITE(IOUT,'(1H )')
          WRITE(IOUT,'(/24H  ====> GLOBAL MESHING   )')
      ENDIF
C
C     RECONSTRUCT CARTESIAN MESH
      J= 0
      ICUR(1)= 1
      ICUR(2)= 1
      ICUR(3)= 1
      IDLD=0
      DO 30 N= 1, 3
         RSTART= 0.0
         ICORD(N)= N
         MINDIM(N)= J+1
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
            IGEO= KEYGEO(ITYP)
            IDLD= IDLDIM(IGEO)
            IDLR= IDLREM(IGEO)
            MINABS= MINDO(IDLD+N)
            MAXABS= MAXDO(IDLD+N)
            RMINUS= RSTART - RMESHO(IDLR+MINABS)
            DO 10 IREM= MINABS, MAXABS-1
               J= J+1
               REMESH(J)= RMESHO(IDLR+IREM)+RMINUS
   10       CONTINUE
            ICUR(N)= 1
            RSTART= RMESHO(IDLR+MAXABS)+RMINUS
   20    CONTINUE
         J= J+1
         REMESH(J)= RSTART
         MAXDIM(N)= J
         IF( IPRT.GE.1.AND.N.LE.NDIM )THEN
            WRITE(IOUT,'(8X,A1,14H-COORDINATES: /(9X,5(1X,F13.6)))')
     >               TEMESH(N), (REMESH(J1),J1=MINDIM(N),MAXDIM(N))
         ENDIF
   30 CONTINUE
      NTOTCX= 3
C
C     RECONSTRUCT CYLINDRICAL MESH
      IF( NDIM.EQ.2 )THEN
         MDMIN= 3
      ELSE
         MDMIN= 1
      ENDIF
      DO 130 N= MDMIN, 3
         ICUR(N)= 1
         NP1= MOD(N  ,3) + 1
         NP2= MOD(N+1,3) + 1
C
C        (XP1,XP2) ARE COORDINATES AT BEGINNING OF BLOCK (IP1,IP2)
         XP2= 0.0
         DO 120 IP2= 1, MAXGRI(NP2)
         XP1= 0.0
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
            IF( ITYP.EQ.0 ) GO TO 105
            IGEO= KEYGEO(ITYP)
            IDLD= IDLDIM(IGEO)
            IDLR= IDLREM(IGEO)
            IF( IGEO.NE.NGEOME )THEN
               NC= IDLDIM(IGEO+1)-IDLD-3
            ELSE
               NC= NTOTCO-IDLD-3
            ENDIF
            IF( NC.EQ.1 )THEN
              IF( ICORDO(IDLD+4).EQ.N )THEN
                NTOTCX= NTOTCX+1
                IF( NTOTCX.GT.NTOTCL )
     >             CALL XABORT( '** XELETR: TOO MANY CYLINDERS' )
                MINP1 = MINDO(IDLD+NP1)
                MINP2 = MINDO(IDLD+NP2)
                MINC  = MINDO(IDLD+4)
                ICORD(NTOTCX)= ICORDO(IDLD+4)
C
C               RECENTER CYLINDERS
                REMESH(J+1)= RMESHO(IDLR+MINC-2)-RMESHO(IDLR+MINP1)+XP1
                REMESH(J+2)= RMESHO(IDLR+MINC-1)-RMESHO(IDLR+MINP2)+XP2
                J= J+2
                MINDIM(NTOTCX)= J+1
                DO 95 IREM= MINC, MAXDO(IDLD+4)
                   J=J+1
                   REMESH(J)= RMESHO(IDLR+IREM)
   95           CONTINUE
                MAXDIM(NTOTCX)= J
                IF( IPRT.GE.1 )THEN
                   WRITE(IOUT,'(13H        CELL(,I8,1H,,I8,1H,,I8,1H),
     >                     3H  (,A1,1H,,A1,10H)- CENTRE: ,
     >                     2H (,2(1X,F13.6),1H) )')
     >                      ICUR(1), ICUR(2), ICUR(3),
     >                      TEMESH(MOD(ICORD(NTOTCX)   ,3)+1),
     >                      TEMESH(MOD(ICORD(NTOTCX)+1,3)+1),
     >                      REMESH(MINDIM(NTOTCX)-2),
     >                      REMESH(MINDIM(NTOTCX)-1)
                   IF( NDIM.EQ.3 )THEN
                     WRITE(IOUT,'(24X,A1,8H-RADII: /(25X,5(1X,F13.6)))')
     >               TEMESH(ICORD(NTOTCX)),
     >               (SQRT(REMESH(J1)),J1=MINDIM(NTOTCX),MAXDIM(NTOTCX))
                   ELSE
                     WRITE(IOUT,'(26X,7HRADII: /(26X,5(1X,F13.6)))')
     >               (SQRT(REMESH(J1)),J1=MINDIM(NTOTCX),MAXDIM(NTOTCX))
                   ENDIF
                ENDIF
              ENDIF
            ENDIF
  105       CONTINUE
            IOLD= ICUR(NP2)
            ICUR(NP2)=   1
            IF( NDIM.EQ.2 )THEN
               IBLK= MAXGRI(1) * (ICUR(2) - 1) + ICUR(1)
            ELSE
               IBLK= MAXGRI(1)*(MAXGRI(2)*ICUR(3)+ICUR(2)-
     >                          MAXGRI(2))+ICUR(1)-MAXGRI(1)
            ENDIF
            ITYP= KEYTYP(IBLK)
            IGEO= KEYGEO(ITYP)
            IDLD= IDLDIM(IGEO)
            IDLR= IDLREM(IGEO)
            MINABS= MINDO(IDLD+NP1)
            MAXABS= MAXDO(IDLD+NP1)
            XP1= XP1 + (RMESHO(IDLR+MAXABS)-RMESHO(IDLR+MINABS))
            ICUR(NP2)= IOLD
  110    CONTINUE
            ICUR(NP1)=   1
            IF( NDIM.EQ.2 )THEN
               IBLK= MAXGRI(1) * (ICUR(2) - 1) + ICUR(1)
            ELSE
               IBLK= MAXGRI(1)*(MAXGRI(2)*ICUR(3)+ICUR(2)-
     >                          MAXGRI(2))+ICUR(1)-MAXGRI(1)
            ENDIF
            ITYP= KEYTYP(IBLK)
            IGEO= KEYGEO(ITYP)
            IDLD= IDLDIM(IGEO)
            IDLR= IDLREM(IGEO)
            MINABS= MINDO(IDLD+NP2)
            MAXABS= MAXDO(IDLD+NP2)
            XP2= XP2 + (RMESHO(IDLR+MAXABS)-RMESHO(IDLR+MINABS))
  120    CONTINUE
  130 CONTINUE
C
C     REESTABLISH INDEXING OF ALL UNKNOWNS
C     NOW, *ICUR()* IS THE INCREMENT FOR CARTESIAN CELL MESHING
      ISU2= 0
      IVO2= 0
      ICREM = 0
      ICUR(3)= 0
      DO 230 IP3= 1,MAXGRI(3)
      ICUR(2)= 0
      DO 220 IP2= 1,MAXGRI(2)
      ICUR(1)= 0
      DO 210 IP1= 1,MAXGRI(1)
         IF( NDIM.EQ.2 )THEN
            IBLK= MAXGRI(1)*(IP2-1)+IP1
         ELSE
            IBLK= MAXGRI(1)*(MAXGRI(2)*IP3+IP2-MAXGRI(2))+IP1-MAXGRI(1)
         ENDIF
         ITYP=   KEYTYP(IBLK)
         ICYL=   KEYCYL(IBLK)
         IGEO=   KEYGEO(ITYP)
         IDLTYX= IDLTYP(ITYP)
         IDLD=   IDLDIM(IGEO)
         IDLGEX= IDLGEO(IGEO)
         IKREM = ICREM
         DO 200 IVS= 1, NVOLO(IGEO)
            NO= NUMBLK(IVS, IBLK)
            IMYG=0
            IF( KEYINT(NO).NE.0 ) GO TO 200
            IMYG=MATGEO(IDLGEX+IVS)
            IVO2= IVO2 + 1
            IVSN= IVO2
            IF( IMYG.GE.0 )THEN
               IF( IVO2.GT.NVOL )
     >           CALL XABORT( '** XELETR: TOO MANY ZONES' )
               KEYMRG( IDLGE2+IVSN)= IMYG+ICREM
               VOLSUR( IDLGE2+IVSN)= VOLSO( IDLGEX+IVS)
               MATALB( IDLGE2+IVSN)= MATTYP( IDLTYX+IVS)
               INDEX(1,IDLGE2+IVSN)= INDEXO(1,IDLGEX+IVS)+ICUR(1)
     >                               + (MINDIM(1)-MINDO(IDLD+1))
               INDEX(2,IDLGE2+IVSN)= INDEXO(2,IDLGEX+IVS)+ICUR(2)
     >                               + (MINDIM(2)-MINDO(IDLD+2))
               INDEX(3,IDLGE2+IVSN)= INDEXO(3,IDLGEX+IVS)+ICUR(3)
     >                             + (MINDIM(3)-MINDO(IDLD+3))
               INDEX(4,IDLGE2+IVSN)= 0
               INCELL( IDLGE2+IVSN)= IBLK
               IF( ICYL.NE.0 )THEN
                  IF( INDEXO(4,IDLGEX+IVS).NE.MAXDO(IDLD+4) )THEN
C                    IF WE ARE INSIDE THE CYLINDER:
                     INDEX(4,IDLGE2+IVSN)=  INDEXO(4,IDLGEX+IVS)
     >                                   + (MINDIM(ICYL)-MINDO(IDLD+4))
                  ENDIF
               ENDIF
              IKREM=IKREM+1
            ELSE
               KEYMRG( IDLGE2+IVSN)= 0
               INCELL( IDLGE2+IVSN)= IBLK
            ENDIF
  200    CONTINUE
         ICREM=IKREM
         DO 400 IVS= -1,NSURO(IGEO),-1
            NO= NUMBLK(IVS, IBLK)
            IMYG=0
            IF( KEYINT(NO).NE.0 ) GO TO 400
            IMYG=MATGEO(IDLGEX+IVS)
            IF( IMYG.LT.0 )THEN
               ISU2= ISU2 - 1
               IVSN= ISU2
               IF( ISU2.LT. NSUR )
     >            CALL XABORT( '** XELETR: TOO MANY SURFACES' )
               KEYMRG( IDLGE2+IVSN)= IVSN
               VOLSUR( IDLGE2+IVSN)= VOLSO( IDLGEX+IVS)
               MATALB( IDLGE2+IVSN)= MATTYP( IDLTYX+IVS)
               INDEX(1,IDLGE2+IVSN)= INDEXO(1,IDLGEX+IVS)+ICUR(1)
     >                               + (MINDIM(1)-MINDO(IDLD+1))
               INDEX(2,IDLGE2+IVSN)= INDEXO(2,IDLGEX+IVS)+ICUR(2)
     >                               + (MINDIM(2)-MINDO(IDLD+2))
               INDEX(3,IDLGE2+IVSN)= INDEXO(3,IDLGEX+IVS)+ICUR(3)
     >                             + (MINDIM(3)-MINDO(IDLD+3))
               INDEX(4,IDLGE2+IVSN)= 0
               INCELL( IDLGE2+IVSN)= IBLK
               IF( ICYL.NE.0 )THEN
                  IF( INDEXO(4,IDLGEX+IVS).NE.MAXDO(IDLD+4) )THEN
C                    IF WE ARE INSIDE THE CYLINDER:
                       INDEX(4,IDLGE2+IVSN)=  INDEXO(4,IDLGEX+IVS)
     >                                   + (MINDIM(ICYL)-MINDO(IDLD+4))
                  ENDIF
               ENDIF
            ENDIF
  400    CONTINUE
         ICUR(1)= ICUR(1) + (MAXDO(IDLD+1)-MINDO(IDLD+1))
  210 CONTINUE
         ICUR(2)= ICUR(2) + (MAXDO(IDLD+2)-MINDO(IDLD+2))
  220 CONTINUE
         ICUR(3)= ICUR(3) + (MAXDO(IDLD+3)-MINDO(IDLD+3))
  230 CONTINUE
C----
C  REMOVE ZONES AND SURFACES WITH VANISHING VOLSUR
C---- 
      IVS=0
      DO 410 ISU2=NSUR,-1
        IF(VOLSUR(ISU2-NSUR+1) .GT. 0.0) THEN
          IVS=IVS+1
          VOLSUR(IVS)=VOLSUR(ISU2-NSUR+1)
          MATALB(IVS)=MATALB(ISU2-NSUR+1)
          KEYMRG(IVS)=KEYMRG(ISU2-NSUR+1)
          INCELL(IVS)=INCELL(ISU2-NSUR+1)
          DO 411 J1=1,4
            INDEX(J1,IVS)=INDEX(J1,ISU2-NSUR+1)
 411      CONTINUE
        ENDIF
 410  CONTINUE
      NSURC=-IVS
      IVS=IVS+1                 
      VOLSUR(IVS)=0.0
      MATALB(IVS)=0
      KEYMRG(IVS)=0
      INCELL(IVS)=0
      DO 420 J1=1,4
        INDEX(J1,IVS)=0
 420  CONTINUE
      DO 430 IVO2=1,NVOL
        IF(VOLSUR(IVO2-NSUR+1) .GT. 0.0) THEN
          IVS=IVS+1               
          VOLSUR(IVS)=VOLSUR(IVO2-NSUR+1)
          MATALB(IVS)=MATALB(IVO2-NSUR+1)
          KEYMRG(IVS)=KEYMRG(IVO2-NSUR+1)
          INCELL(IVS)=INCELL(IVO2-NSUR+1)
          DO 431 J1=1,4
            INDEX(J1,IVS)=INDEX(J1,IVO2-NSUR+1)
 431      CONTINUE
        ENDIF
 430  CONTINUE
      NVOLC=IVS+NSURC-1
      KL=1-NSURC
      IF( IPRT.GE.5 )THEN
         WRITE(IOUT,'(1H )')
         WRITE(IOUT,'(/13H RENUMBERING ,I8,13H VOLUMES AND ,'//
     >           'I8,10H SURFACES.)') NVOL,-NSUR
         WRITE(IOUT,'(20H CARTESIAN MESH      ,7HMINDIM=,3I8)')
     >                  (MINDIM(J1),J1=1,3)
         WRITE(IOUT,'(20X                     ,7HMAXDIM=,3I8)')
     >                  (MAXDIM(J1),J1=1,3)
         IF( NTOTCL.GT.3 )THEN
            DO 540 J1= 4, NTOTCL
               WRITE(IOUT,'(10H CYLINDER ,I8,6X     ,7HMINDIM=,15X,I8)')
     >                                 J1-3,MINDIM(J1)
               WRITE(IOUT,'(20X                     ,7HMAXDIM=,15X,I8)')
     >                                      MAXDIM(J1)
  540       CONTINUE
         ENDIF
         DO 550 IVS= NSURC, NVOLC
             IF( IVS.LT.0 )THEN
                IF( KEYMRG(IVS+KL) .EQ. 0 )THEN
                   WRITE(IOUT,'(8H KEYMRG(,I8,2H)=,I8,
     >             7H INDEX=,4I8,7H BLOCK=,I8,9H SURFACE=,F20.7,
     >             17H ABSENT FROM CELL)')
     >             IVS,KEYMRG(IVS+KL),(INDEX(J1,IVS+KL),J1=1,4),
     >             INCELL(IVS+KL),4.*VOLSUR(IVS+KL)
                ELSE
                   WRITE(IOUT,'(8H KEYMRG(,I8,2H)=,I8,
     >             7H INDEX=,4I8,7H BLOCK=,I8,9H SURFACE=,F20.7)')
     >             IVS,KEYMRG(IVS+KL),(INDEX(J1,IVS+KL),J1=1,4),
     >             INCELL(IVS+KL),4.*VOLSUR(IVS+KL)
                ENDIF
             ELSE
                IF( KEYMRG(IVS+KL) .EQ. 0 )THEN
                   WRITE(IOUT,'(8H KEYMRG(,I8,2H)=,I8,
     >             7H INDEX=,4I8,7H BLOCK=,I8,9H VOLUME= ,F20.7,
     >             17H ABSENT FROM CELL)')
     >             IVS,KEYMRG(IVS+KL),(INDEX(J1,IVS+KL),J1=1,4),
     >             INCELL(IVS+KL),VOLSUR(IVS+KL)
                ELSE
                   WRITE(IOUT,'(8H KEYMRG(,I8,2H)=,I8,
     >             7H INDEX=,4I8,7H BLOCK=,I8,9H VOLUME= ,F20.7)')
     >             IVS,KEYMRG(IVS+KL),(INDEX(J1,IVS+KL),J1=1,4),
     >             INCELL(IVS+KL),VOLSUR(IVS+KL)
                ENDIF
             ENDIF
  550    CONTINUE                              
      ENDIF
      RETURN
      END
