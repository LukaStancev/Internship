*DECK XELBIN
      SUBROUTINE XELBIN( IPGEOM,   NDIM, NGEOME, L1CELL, NTYPES,  NGIDL,
     >                    NTIDL, NBLOCK, MAXGRI,  NUNKO,   IPRT,  CELLG,
     >                    NSURO,  NVOLO, IDLGEO, MATGEO, KEYGEO, IDLTYP,
     >                   IDLBLK, KEYTYP, MATTYP, KEYINT)
************************************************************************
*                                                                      *
*           NAME: XELBIN                                               *
*      COMPONENT: EXCELL                                               *
*          LEVEL: 3 (CALLED BY 'XELTRK')                               *
*        VERSION: 1.0                                                  *
*       CREATION: 87/01                                                *
*       MODIFIED: 00/03 (R.R.) DECLARE ALL VARIABLE TYPES              *
*         AUTHOR: ROBERT ROY                                           *
*                                                                      *
*     SUBROUTINE: THIS ROUTINE WILL IDENTIFY EVERY ZONE OF EVERY TYPE  *
*                 TO ITS MATERIAL; IT WILL ALSO INTERFACE ALL INTERNAL *
*                 SURFACES FOR CELLS PRESENT IN THE SUPERCELL.         *
*                                                                      *
*--------+-------------- V A R I A B L E S -------------+--+-----------*
*  NAME  /                  DESCRIPTION                 /IO/MOD(DIMENS)*
*--------+----------------------------------------------+--+-----------*
* IPGEOM / POINTER TO THE GEOMETRY (L_GEOM)             /I./INT        *
* NDIM   / # OF DIMENSIONS (2 OR 3).                    /I./INT        *
* NGEOME / # OF GEOMETRIES.                             /I./INT        *
* L1CELL / .TRUE. IF ONLY ONE CELL.                     /I./LOG        *
* NTYPES / # OF TYPES.                                  /I./INT        *
* NGIDL  / LENGHT OF GEOMETRIC NUMBERING.               /I./INT        *
* NTIDL  / LENGHT OF TYPE NUMBERING.                    /I./INT        *
* NBLOCK / # OF BLOCKS.                                 /I./INT        *
* MAXGRI / # OF CELLS ALONG EACH AXIS.                  /I./INT(3)     *
* NUNKO  / OLD # OF UNKNOWNS.                           /I./INT        *
* IPRT   / INTERMEDIATE PRINTING LEVEL FOR OUTPUT.      /I./INT        *
* CELLG  / TO KEEP GEOMETRY  NAMES.                     /I./C*4(3*NGEOM*
* NSURO  / # OF SURFACES OF EACH GEOMETRY.              /I./INT(NGEOME)*
* NVOLO  / # OF ZONES OF EACH GEOMETRY.                 /I./INT(NGEOME)*
* IDLGEO / POSITION OF EACH GEOMETRY IN THE             /I./INT(NGEOME)*
*        /            GEOMETRY NUMBERING SCHEME.        /  /           *
* MATGEO / MATERIAL #S CORRESPONDING TO GEOMETRIES.     /I./INT(NGIDL) *
* KEYGEO / GEOMETRIC KEY FOR EACH TYPE.                 /I./INT(NTYPES)*
* IDLTYP / POSITION OF EACH TYPE IN NUMBERING SCHEME.   /I./INT(NTYPES)*
* IDLBLK / POSITION OF EACH BLOCK IN NUMBERING SCHEME.  /I./INT(NBLOCK)*
* KEYTYP / TYPE KEY FOR EACH BLOCK.                     /I./INT(NBLOCK)*
* MATTYP / MATERIAL #S FOR ZONES OF EVERY TYPE.         /.O/INT(NTIDL) *
* KEYINT / INTERFACE KEY (GIVING THE CONNECTED SURFACE)./.O/INT(NUNKO )*
************************************************************************
C
      USE                GANLIB
      IMPLICIT           NONE
C
      TYPE(C_PTR)        IPGEOM 
      INTEGER            NDIM, NGEOME, NTYPES, NGIDL, NTIDL, NBLOCK,
     >                   NUNKO, IPRT
      INTEGER             NSURO(NGEOME),  NVOLO(NGEOME), IDLGEO(NGEOME),
     >                   MATGEO( NGIDL), KEYGEO(NTYPES), IDLTYP(NTYPES),
     >                   MATTYP( NTIDL), KEYTYP(NBLOCK), IDLBLK(NBLOCK),
     >                   KEYINT(NUNKO ), MAXGRI(NDIM)  , CELLG(3*NTYPES)
C
      INTEGER            ILO(3,2), NO(2), KTYP(2),
     >                     KMAT(2),    KSUR(2), KABSO(2), KSID(2),
     >                   ICOORD(3),   NCODE(6)
      CHARACTER          GEOCEL*12, TEDATA*12, TEMESH(4)*7
      LOGICAL            SWKILL, L1CELL, LL1, LL2
      INTEGER            NSTATE, IOUT, MAXSPL
      PARAMETER        ( NSTATE=40, IOUT=6, MAXSPL=100 )
      INTEGER            ISTATE(NSTATE),ISPLT(MAXSPL)
      INTEGER            NUMGEO, NUMTYP, NUMBLK, I, K
      INTEGER            NBMD
      INTEGER            IMYT, IUNK, ITYP, IMYG, IGEO, NSUX, NVOX, ICYL,
     >                   ICX, ICY, ICZ, LR, LX, LY, LZ, KOLD, ITYPG,
     >                   ISUR, IX, IY, IZ, IOFF, KNEW, ILEN, ITYLCM,
     >                   ISX, ISY, ISZ, ISR, KIOFX, KIOFY, KIOFZ,
     >                   J0, J1, J2, JC, JR, IP0, IP1, IP2, N, NP1, NP2,
     >                   K0, K1, K2, K3, KR, IBLK, ISUX
      EQUIVALENCE      ( ICOORD(1),LX ),(ICOORD(2),LY),(ICOORD(3),LZ )
      DATA        TEMESH / 'X', 'Y', 'Z', 'R'/
C
      NUMGEO(I,K)= I + IDLGEO(K)
      NUMTYP(I,K)= I + IDLTYP(K)
      NUMBLK(I,K)= I + IDLBLK(K)
C
      SWKILL= .FALSE.
      LL1= .FALSE.
      LL2= .FALSE.
      DO 10 IMYT= 1,  NTIDL
          MATTYP(IMYT)=  0
   10 CONTINUE
      DO 20 IUNK= 1, NUNKO
          KEYINT(IUNK)=  0
   20 CONTINUE
      DO 40 ITYP= 1, NTYPES
         IGEO  = KEYGEO( ITYP   )
         NVOX  = NVOLO( IGEO )
         NSUX  = NSURO( IGEO ) 
         IF( .NOT.L1CELL )THEN
            WRITE(GEOCEL( 1: 4), '(A4)') CELLG(3*ITYP-2)
            WRITE(GEOCEL( 5: 8), '(A4)') CELLG(3*ITYP-1)
            WRITE(GEOCEL( 9:12), '(A4)') CELLG(3*ITYP  )
            CALL LCMSIX(IPGEOM, GEOCEL, 1)
         ELSE
            CALL LCMGET(IPGEOM,'NCODE', NCODE)
            LL1=((NCODE(2).EQ.3).AND.(NCODE(3).EQ.3))
            LL2=((NCODE(1).EQ.3).AND.(NCODE(4).EQ.3))
         ENDIF
         CALL XDISET(ISTATE,NSTATE,0)
         CALL LCMGET(IPGEOM,'STATE-VECTOR',ISTATE)
         ITYPG= ISTATE(1)
         IF( ITYPG.EQ.20) THEN
C           FOR *CARCEL* GEOMETRIES
            ICYL= 1
            ICX=  1
            ICY=  2
            ICZ=  3
         ELSEIF(ITYPG.EQ.3.OR.ITYPG.EQ.6 )THEN
C           FOR *CARCEL*, *TUBE* OR *TUBEZ* GEOMETRIES
            ICYL= 1
            ICX=  1
            ICY=  2
            ICZ=  3
            IF( LL1.OR.LL2 )THEN
               CALL XABORT( 'XELBIN: DIAGONAL SYMETRIES NOT POSSIBLE')
            ENDIF
         ELSEIF( ITYPG.GT.20 )THEN
C           FOR *CARCELX*, *CARCELY* OR *CARCELZ*
            ICYL= 1
            ICZ= ITYPG-20
            ICX= MOD(ICZ  , 3) + 1
            ICY= MOD(ICZ+1, 3) + 1
         ELSE
C           FOR *CAR2D* OR *CAR3D*
            ICYL= 0
            ICX=  1
            ICY=  2
            ICZ=  3
         ENDIF
         LR=    ISTATE(2)
         LX=    MAX(1,ISTATE(3))
         LY=    MAX(1,ISTATE(4))
         LZ=    MAX(1,ISTATE(5))
         KOLD=  ISTATE(6)
         DO 30 ISUR= NSUX, -1
            MATTYP(NUMTYP(ISUR,ITYP))= MATGEO(NUMGEO(ISUR,IGEO))
   30    CONTINUE
C
C        GET MIXTURE NUMBERS
         CALL LCMLEN(IPGEOM, 'MIX', ILEN, ITYLCM)
         IF( ILEN.NE.KOLD )THEN
            WRITE(IOUT,*) 'LENGHT(MIX)=  ',ILEN
            WRITE(IOUT,*) '# OF VOLUMES= ',KOLD
            CALL LCMLIB(IPGEOM)
            CALL XABORT( 'XELBIN: INVALID NUMBER OF MIXTURES')
         ENDIF
         CALL LCMGET(IPGEOM,'MIX',MATTYP(NUMTYP(1,ITYP)))
C
C        IN THE CASE OF DIAGONAL SYMMETRY IN 'ONE-CELL'
C        CAR2D AND CAR3D  GEOMETRY UNFOLD MIXTURES 
C
      IF(ITYPG .LT. 20) THEN
         K3=ISTATE(6)
         NBMD=(LZ*LY*(LX+1))/2
         IF(K3 .EQ. NBMD) THEN
C----
C MIXTURE ENTERED IN DIAGONAL FORM
C----
         IF( LL1 )THEN
            DO 70 IZ=LZ,1,-1
            IOFF=(IZ-1)*LX*LY
            DO 70 IY=LY,1,-1
            DO 60 IX=LX,IY+1,-1
   60       MATTYP(NUMTYP(IOFF+(IY-1)*LX+IX,ITYP))=
     >                    MATTYP(NUMTYP(IOFF+(IX-1)*LY+IY,ITYP))
            DO 70 IX=IY,1,-1
            MATTYP(NUMTYP(IOFF+(IY-1)*LX+IX,ITYP))=
     >                    MATTYP(NUMTYP(K3,ITYP))
            K3=K3-1
   70       CONTINUE
            KOLD= LX*LY*LZ
         ELSEIF( LL2 )THEN
            DO 80 IZ=LZ,1,-1
            IOFF=(IZ-1)*LX*LY
            DO 80 IY=LY,1,-1
            DO 80 IX=LX,IY,-1
            MATTYP(NUMTYP(IOFF+(IY-1)*LX+IX,ITYP))=
     >                    MATTYP(NUMTYP(K3,ITYP))
            K3=K3-1
   80       CONTINUE
            DO 90 IZ=1,LZ
            IOFF=(IZ-1)*LX*LY
            DO 90 IY=1,LY
            DO 90 IX=1,IY-1
   90       MATTYP(NUMTYP(IOFF+(IY-1)*LX+IX,ITYP))=
     >               MATTYP(NUMTYP(IOFF+(IX-1)*LY+IY,ITYP))
            KOLD= LX*LY*LZ
         ENDIF
         ENDIF
      ENDIF
C
C        FOR THE PARTICULAR CASE OF *TUBE* OR *TUBEZ* GEOMETRIES
         IF( ITYPG.EQ.3.OR.ITYPG.EQ.6 )THEN
            DO 39 IZ= 1, LZ
               MATTYP(NUMTYP(KOLD+IZ,ITYP))= -2
   39       CONTINUE
            KOLD= KOLD+LZ
         ENDIF
C
C        FILL UP MATTYP ACCORDING TO SPLITTING VALUES.
         KNEW= NVOX
         ISR= 0
         ISX= 0
         ISY= 0
         ISZ= 0
         DO 303 K0= ICOORD(ICZ),1,-1
            KIOFZ= KOLD
            TEDATA= 'SPLIT'//TEMESH(ICZ)
            CALL LCMLEN(IPGEOM,TEDATA,ILEN,ITYLCM)
            IF( ILEN.GT.MAXSPL )THEN
               CALL XABORT('XELBIN: SPLIT OVERFLOW ('//TEDATA//')')
            ELSEIF( ILEN.EQ.0 )THEN
               ISZ= 1
            ELSE
               CALL LCMGET(IPGEOM,TEDATA,ISPLT)
               ISZ= ISPLT(K0)
            ENDIF
         DO 303 J0=ISZ,1,-1
            KOLD= KIOFZ
         DO 303 K1= ICOORD(ICY),1,-1
            KIOFY= KOLD
            TEDATA= 'SPLIT'//TEMESH(ICY)
            CALL LCMLEN(IPGEOM,TEDATA,ILEN,ITYLCM)
            IF( ILEN.GT.MAXSPL )THEN
               CALL XABORT('XELBIN: SPLIT OVERFLOW ('//TEDATA//')')
            ELSEIF( ILEN.EQ.0 )THEN
               ISY= 1
            ELSE
               CALL LCMGET(IPGEOM,TEDATA,ISPLT)
               ISY= ISPLT(K1)
            ENDIF
         DO 303 J1=ISY,1,-1
            KOLD= KIOFY
         DO 303 K2= ICOORD(ICX),1,-1
            KIOFX= KOLD
            TEDATA= 'SPLIT'//TEMESH(ICX)
            CALL LCMLEN(IPGEOM,TEDATA,ILEN,ITYLCM)
            IF( ILEN.GT.MAXSPL )THEN
               CALL XABORT('XELBIN: SPLIT OVERFLOW ('//TEDATA//')')
            ELSEIF( ILEN.EQ.0 )THEN
               ISX= 1
            ELSE
               CALL LCMGET(IPGEOM,TEDATA,ISPLT)
               ISX= ISPLT(K2)
            ENDIF
         DO 303 J2=ISX,1,-1
            KOLD= KIOFX
C           FOR RECTANGULAR OUTER REGIONS.
            IMYT= MATTYP(NUMTYP(KOLD,ITYP))
            MATTYP(NUMTYP(KNEW,ITYP))= IMYT
            KNEW= KNEW-1
            KOLD= KOLD-1
            IF( ICYL.EQ.1 )THEN
C              FOR CYLINDRICAL INNER REGIONS.
               DO 302 KR= LR,1,-1
                  TEDATA= 'SPLIT'//TEMESH(4)
                  CALL LCMLEN(IPGEOM,TEDATA,ILEN,ITYLCM)
                  IF( ILEN.GT.MAXSPL )THEN
                    CALL XABORT('XELBIN: SPLIT OVERFLOW ('//TEDATA//')')
                  ELSEIF( ILEN.EQ.0 )THEN
                     ISR= 1
                  ELSE
                     CALL LCMGET(IPGEOM,TEDATA,ISPLT)
                     ISR= ABS(ISPLT(KR))
                  ENDIF
                  IMYT= MATTYP(NUMTYP(KOLD,ITYP))
               DO 301 JR=ISR,1,-1
                  MATTYP(NUMTYP(KNEW,ITYP))= IMYT
                  KNEW= KNEW-1
  301          CONTINUE
                  KOLD= KOLD-1
  302          CONTINUE
            ENDIF
  303    CONTINUE
         IF( KNEW.NE.0 )THEN
            WRITE(IOUT,*) 'XELBIN: KNEW.NE.0 = PROBLEM WITH SPLITTING'
            SWKILL= .TRUE.
         ENDIF
         IF( KOLD.NE.0 )THEN
            WRITE(IOUT,*) 'XELBIN: KOLD.NE.0 = PROBLEM WITH SPLITTING'
            SWKILL= .TRUE.
         ENDIF
C
         IF( .NOT.L1CELL ) CALL LCMSIX(IPGEOM, ' ',    2)
   40 CONTINUE
C
C     RECOMPOSE INTERNAL SURFACES COUPLING (INTERFACES)
C        THIS ASSUMES THAT AN ORDERING OF SURFACES IS DONE
C        BECAUSE:  SIDE-BY-SIDE INTERFACES
C                  ARE SUPPOSED IN INCREASING POSITION.
      DO 220 N= 1, NDIM
C
C        DEFINITION OF THE SIDE NUMBER TO COUPLE.
         KSID(1)=  -2*N
         KSID(2)= (-2*N) + 1
         NP1   = MOD(N  ,NDIM) + 1
         IF( NDIM.EQ.3 )THEN
            NP2   = MOD(N+1,NDIM) + 1
            DO 110 IP1= 1, MAXGRI(NP1)
               ILO(NP1,1)= IP1
               ILO(NP1,2)= IP1
            DO 110 IP2= 1, MAXGRI(NP2)
               ILO(NP2,1)= IP2
               ILO(NP2,2)= IP2
            DO 110 IP0= 1, MAXGRI(N)-1
               ILO(N  ,1)= IP0
               ILO(N  ,2)= IP0 + 1
               DO 100  JC= 1, 2
                  NO(JC)= MAXGRI(1)*(MAXGRI(2)*ILO(3,JC)+ILO(2,JC)-
     >                               MAXGRI(2))+ILO(1,JC)-MAXGRI(1)
                  KTYP(JC)= KEYTYP( NO(JC) )
                  IF( KTYP(JC).EQ.0 ) GO TO 110
                  IGEO  = KEYGEO( KTYP(JC) )
C                 SEARCH FROM THE END
                  KSUR(JC)= NSURO(IGEO)
                  KMAT(JC)= MATTYP( NUMTYP(KSUR(JC),KTYP(JC)) )
  100          CONTINUE
C
C              ORDERING INTERFACING OF THE TWO BLOCKS.
  101          CONTINUE
                  IF( KMAT(1).EQ.KSID(1).AND.KMAT(2).EQ.KSID(2) )THEN
                     IF( KSUR(1).EQ.0 .OR. KSUR(2).EQ.0 ) GO TO 109
                     KABSO(1)= NUMBLK( KSUR(1),NO(1) )
                     KABSO(2)= NUMBLK( KSUR(2),NO(2) )
                     KEYINT( KABSO(1) )= KABSO(2)
                     KEYINT( KABSO(2) )= KABSO(1)
                     KSUR(1)= KSUR(1)+1
                     KSUR(2)= KSUR(2)+1
                  ELSE
                     IF( KMAT(1).NE.KSID(1) ) KSUR(1)= KSUR(1)+1
                     IF( KMAT(2).NE.KSID(2) ) KSUR(2)= KSUR(2)+1
                  ENDIF
                  IF( KSUR(1).NE.0 )THEN
                     KMAT(1)= MATTYP( NUMTYP(KSUR(1),KTYP(1)) )
                  ELSE
                     KMAT(1)= KSID(1)
                  ENDIF
                  IF( KSUR(2).NE.0 )THEN
                     KMAT(2)= MATTYP( NUMTYP(KSUR(2),KTYP(2)) )
                  ELSE
                     KMAT(2)= KSID(2)
                  ENDIF
               GO TO 101
  109          IF( KSUR(1).NE.0 .OR. KSUR(2).NE.0 )THEN
                  WRITE(IOUT,'(1H ,I8,4H OF ,I8,5H <=> ,I8,4H OF ,I8)')
     >                          KSUR(1),  NO(1),     KSUR(2),  NO(2)
                  SWKILL=.TRUE.
               ENDIF
  110       CONTINUE
         ELSEIF( NDIM.EQ.2 )THEN
            DO 210 IP1= 1, MAXGRI(NP1)
               ILO(NP1,1)= IP1
               ILO(NP1,2)= IP1
            DO 210 IP0= 1, MAXGRI(N)-1
               ILO(N  ,1)= IP0
               ILO(N  ,2)= IP0 + 1
               DO 200  JC= 1, 2
                  NO(JC)= MAXGRI(1) * (ILO(2,JC) - 1) + ILO(1,JC)
                  KTYP(JC)= KEYTYP( NO(JC) )
                  IF( KTYP(JC).EQ.0 ) GO TO 210
                  IGEO  = KEYGEO( KTYP(JC) )
C                 SEARCH FROM THE END
                  KSUR(JC)= NSURO(IGEO)
                  KMAT(JC)= MATTYP( NUMTYP(KSUR(JC),KTYP(JC)) ) 
  200          CONTINUE
C
C              ORDERING INTERFACING OF THE TWO BLOCKS.
  201          CONTINUE
                  IF( KMAT(1).EQ.KSID(1).AND.KMAT(2).EQ.KSID(2) )THEN
                     IF( KSUR(1).EQ.0 .OR. KSUR(2).EQ.0 ) GO TO 209
                     KABSO(1)= NUMBLK( KSUR(1),NO(1) )
                     KABSO(2)= NUMBLK( KSUR(2),NO(2) )
                     KEYINT( KABSO(1) )= KABSO(2)
                     KEYINT( KABSO(2) )= KABSO(1)
                     KSUR(1)= KSUR(1)+1
                     KSUR(2)= KSUR(2)+1
                  ELSE
                     IF( KMAT(1).NE.KSID(1) ) KSUR(1)= KSUR(1)+1
                     IF( KMAT(2).NE.KSID(2) ) KSUR(2)= KSUR(2)+1
                  ENDIF
                  IF( KSUR(1).NE.0 )THEN
                     KMAT(1)= MATTYP( NUMTYP(KSUR(1),KTYP(1)) )
                  ELSE
                     KMAT(1)= KSID(1)
                  ENDIF
                  IF( KSUR(2).NE.0 )THEN
                     KMAT(2)= MATTYP( NUMTYP(KSUR(2),KTYP(2)) )
                  ELSE
                     KMAT(2)= KSID(2)
                  ENDIF
               GO TO 201
  209          IF( KSUR(1).NE.0 .OR. KSUR(2).NE.0 )THEN
                  WRITE(IOUT,'(1H ,I8,4H OF ,I8,5H <=> ,I8,4H OF ,I8)')
     >                          KSUR(1),  NO(1),     KSUR(2),  NO(2)
                  SWKILL=.TRUE.
               ENDIF
  210       CONTINUE
         ELSE
            CALL LCMLIB(IPGEOM)
            CALL XABORT( 'XELBIN: *** FALSE NDIM VALUE')
         ENDIF
  220 CONTINUE
C
      IF( IPRT.GE.100 .OR. SWKILL )THEN
         IUNK= 0
         WRITE(IOUT,'(/40H       KEYINT       COUPLE    MATERIAL  )')
         DO 250 IBLK= 1, NBLOCK
            ITYP= KEYTYP(IBLK)
            IGEO= KEYGEO(ITYP)
            NVOX= NVOLO(IGEO)
            NSUX= NSURO(IGEO)
            DO 240 ISUX= NSUX, NVOX
               IUNK= IUNK+1
               IMYT= MATTYP( NUMTYP(ISUX,ITYP) )
               IF( ISUX.LT.0 )THEN
                  WRITE(IOUT,
     >            '(5H SUR(,I8,5H) => ,I8,4H OF ,    I8)')
     >                      IUNK,      KEYINT(IUNK), IMYT
               ELSEIF( ISUX.GT.0 )THEN
                  IMYG= MATGEO( NUMGEO(ISUX,IGEO) )
                  WRITE(IOUT,
     >            '(5H VOL(,I8,5H) => ,I8,4H OF ,    I8,1H(,I8,1H))')
     >                      IUNK,      KEYINT(IUNK), IMYT,  IMYG
               ENDIF
  240       CONTINUE
            WRITE(IOUT,'(/1X)')
  250    CONTINUE
      ENDIF
      IF( SWKILL ) CALL XABORT( 'XELBIN: IMPOSSIBLE TO INTERFACE')
C
      RETURN
      END
