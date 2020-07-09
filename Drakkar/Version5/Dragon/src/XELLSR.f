*DECK XELLSR
      SUBROUTINE XELLSR(  NDIM,    NCP,   NSUR, MAXREM, REMESH,
     >                   INDEL, MINDIM, MAXDIM, ICOORD,   ICUR,   INCR,
     >                  TRKORI, TRKDIR, TRKCUT,  NSCUT,  NCROS,
     >                  TOTLEN)
************************************************************************
*                                                                      *
*           NAME: XELLSR                                               *
*      COMPONENT: EXCELL                                               *
*          LEVEL: 4 (CALLED BY 'XELTI2' & 'XELTI3' & 'XELTS2')         *
*        VERSION: 1.0                                                  *
*       CREATION: 90/05                                                *
*       MODIFIED: 91/07 (R.R.)                                         *
*                 00/03 (R.R.) DECLARE ALL VARIABLE TYPES              *
*         AUTHOR: ROBERT ROY                                           *
*                                                                      *
*     SUBROUTINE: THIS ROUTINE WILL FIND THE BEGINNING AND ENDING      *
*                 SURFACES CROSSED BY A TRACK.                         *
*                                                                      *
*--------+-------------- V A R I A B L E S -------------+--+-----------*
*  NAME  /                  DESCRIPTION                 /IO/MOD(DIMENS)*
*--------+----------------------------------------------+--+-----------*
* NDIM   / # OF DIMENSION (2 OR 3).                     /I./INT        *
* NCP    / # OF CYLINDRES OF A TYPE + 3    (.LT.NC3MAX)./I./INT        *
* NSUR   / # OF SURFACES.                               /I./INT        *
* MAXREM / MAX NUMBER OF REAL MESH VALUES IN "REMESH".  /I./INT        *
* REMESH / REAL MESH VALUES (RECT/CYL).                 /I./REL(MAXREM)*
* INDEL  / #ING OF SURFACES & ZONES.                    /I./INT(4*NVS) *
* MINDIM / MIN INDEX VALUES FOR ALL AXES (RECT/CYL).    /I./INT(NC3MAX)*
* MAXDIM / MAX INDEX VALUES FOR ALL AXES (RECT/CYL).    /I./INT(NC3MAX)*
* ICOORD / PRINCIPAL AXES DIRECTION (X/Y/Z) FOR MESHES. /I./INT(NC3MAX)*
* ICUR   / CURRENT ZONAL LOCATION FOR A TRACK SEGMENT.  /../INT(NC3MAX)*
* INCR   / INCREMENT DIRECTION FOR NEXT TRACK SEGMENT.  /../INT(NC3MAX)*
* TRKORI / ORIGIN OF A TRACK.                           /I./REL(3)     *
* TRKDIR / DIRECTION OF A TRACK IN ALL AXES.            /I./REL(3)     *
* TRKCUT / POINTS WHERE TRACK CUT THE DOMAIN.           /.O/REL(3,2)   *
* NSCUT  / SURFACE   WHERE THE TRACK BEGINS/ENDS        /.O/INT(2)     *
* NCROS  / # OF SURFACE CROSSING                        /.O/INT        *
* TOTLEN / TOTAL LENGTH OF THE TRACK                    /.O/REL        *
************************************************************************
C
      IMPLICIT NONE
C
      INTEGER NDIM, NCP, NSUR, MAXREM, NCROS
      REAL    TRKCUT(3,2), REMESH(MAXREM), TRKDIR(3), TRKORI(3), TOTLEN,
     >        TKBEG1, TKBEG2, R2BEG
      INTEGER MINDIM(NCP), MAXDIM(NCP), ICUR(NCP), INCR(NCP),
     >        ICOORD(NCP), IFACUT(2), ISFCUT(2),
     >        IORD(4), NSCUT(2), INDEL(4,*)
      INTEGER     IOUT
      PARAMETER ( IOUT=6 )
      REAL        ANORM2, CENTRE, A, B, XYZP2, CONBEG, CONEND, CON,
     >            XYZP1
      INTEGER     I, J, NUM, NCRBEG, NCREND, NP1, NUMP1, NP2, NUMP2,
     >            N, NUMP0, K, IBEGIN, KELSUR, KWW, IDM
C
      ANORM2(A,B)= A*A + B*B
      CENTRE(I,J)= REMESH( MAXDIM(I-1) + J )
      NUM(J)= J + 1 - NSUR
      NUMP2= 0
      CALL XDISET(IFACUT,2,0)
      CALL XDISET(ISFCUT,2,0)
C
      IF( NDIM.EQ.2 )THEN
         NP2= 3
         NUMP2 = ICOORD(NP2)
         XYZP2= 0.0
      ENDIF
C
C     IF THERE ARE NO CYLINDER AT ALL
      NSCUT(1)= 0
      NSCUT(2)= 0
      NCRBEG= 0
      NCREND= 0
      CONBEG=+1.0E+36
      CONEND=-1.0E+36
C
C     FING BEGINNING AND ENDING POINTS OF THE TRACK
      DO 75 N   = 1, NDIM
         NUMP0 = ICOORD(N  )
         IF( INCR(NUMP0).EQ.0 ) GO TO 75
         NP1   = MOD(N   ,NDIM)+1
         NUMP1 = ICOORD(NP1)
         IF( NDIM.EQ.3 )THEN
            NP2   = MOD(N+1 ,NDIM)+1
            NUMP2 = ICOORD(NP2)
         ENDIF
         DO 70 IDM = MINDIM(N), MAXDIM(N), MAXDIM(N)-MINDIM(N)
            CON   = (REMESH(IDM) - TRKORI(NUMP0)) / TRKDIR(NUMP0)
            XYZP1 = TRKORI(NUMP1) + CON * TRKDIR(NUMP1)
            IF( XYZP1.LT.REMESH(MINDIM(NP1)).OR.
     >          XYZP1.GT.REMESH(MAXDIM(NP1))) GO TO 70
            IF( NDIM.EQ.3 )THEN
               XYZP2 = TRKORI(NUMP2) + CON * TRKDIR(NUMP2)
               IF( XYZP2.LT.REMESH(MINDIM(NP2)).OR.
     >             XYZP2.GT.REMESH(MAXDIM(NP2))) GO TO 70
            ENDIF
            IF( CON.LT.CONBEG )THEN
               NCRBEG=1
               NCREND=MAX(1,NCREND)
               IFACUT(1)= NUMP0
               ISFCUT(1)= IDM
               IF( IDM.EQ.MINDIM(N) ) ISFCUT(1)= ISFCUT(1)-1
               CONBEG=CON
               TRKCUT(NUMP0,1)= REMESH(IDM)
               TRKCUT(NUMP1,1)= XYZP1
               TRKCUT(NUMP2,1)= XYZP2
            ENDIF
            IF( CON.GT.CONEND )THEN
               NCREND=2
               NCRBEG=MIN(2,NCRBEG)
               IFACUT(2)= NUMP0
               ISFCUT(2)= IDM
               IF( IDM.EQ.MINDIM(N) ) ISFCUT(2)= ISFCUT(2)-1
               CONEND=CON
               TRKCUT(NUMP0,2)= REMESH(IDM)
               TRKCUT(NUMP1,2)= XYZP1
               TRKCUT(NUMP2,2)= XYZP2
            ENDIF
   70    CONTINUE
   75 CONTINUE
      NCROS = NCREND + NCRBEG
      TOTLEN= CONEND - CONBEG
      IF( NCROS.EQ.0 ) GO TO 1000
      NCROS = NCREND + 1 - NCRBEG
C
C     FIND BEGINNING AND ENDING SURFACE NUMBERS
      DO 900 K= NCRBEG, NCREND
         DO 90 I   = 1, NDIM
            N      = ICOORD(I)
            ICUR(I)= MINDIM(I)
            DO 80 J = MINDIM(I), MAXDIM(I)-1
               IF(TRKCUT(N,K).GE.REMESH(J)) ICUR(I)= J
   80       CONTINUE
   90    CONTINUE
         ICUR(IFACUT(K))= ISFCUT(K)
         IBEGIN= MAXDIM(3) + 3
         DO 110 I  = 4, NCP
            N     = ICOORD(I)
            NP1   = MOD(N  ,3) + 1
            NP2   = MOD(N+1,3) + 1
            TKBEG1= CENTRE(I,1) - TRKCUT(NP1,K)
            TKBEG2= CENTRE(I,2) - TRKCUT(NP2,K)
            R2BEG = ANORM2(TKBEG1,TKBEG2)
            ICUR(I)  = IBEGIN - 1
            DO 100 J  = IBEGIN, MAXDIM(I)
               IF( R2BEG    .GE. REMESH(J) )ICUR(I)= J
  100       CONTINUE
            IBEGIN= MAXDIM(I) + 3
  110    CONTINUE
C
C        FIND IORD(4) FOR LOCATION IN THE INDEX VECTOR
         DO 115 I= 1,NCP
            IORD(MIN(4,I))= ICUR(I)
            IF( I.GT.3.AND.ICUR(I).LT.MAXDIM(I)) GOTO 116
  115    CONTINUE
         IORD(4)= 0
  116    CONTINUE
C
C        FIND NSCUT=BEGINNING/ENDING SURFACE #S
         KELSUR= NSUR
         INDEL(1,NUM(0))= IORD(1)
         INDEL(2,NUM(0))= IORD(2)
         INDEL(3,NUM(0))= IORD(3)
         INDEL(4,NUM(0))= IORD(4)
  880    CONTINUE
            IF( IORD(1).EQ.INDEL(1,NUM(KELSUR)).AND.
     >          IORD(2).EQ.INDEL(2,NUM(KELSUR)).AND.
     >          IORD(3).EQ.INDEL(3,NUM(KELSUR)).AND.
     >          IORD(4).EQ.INDEL(4,NUM(KELSUR)) ) GO TO 890
            KELSUR= KELSUR + 1
         GO TO 880
  890    NSCUT(K)= KELSUR
         IF( KELSUR.EQ.0 )THEN
            WRITE(IOUT,*) '         BAD SURFACE IDENTIFICATION'
            WRITE(IOUT,*) ' NSCUT=', NSCUT(K)
            WRITE(IOUT,*) 'TRKCUT=', (TRKCUT(KWW,K),KWW=1,3)
            WRITE(IOUT,*) '  IORD=', IORD
            WRITE(IOUT,*) '  ICUR=', (ICUR(KWW),KWW=1,NCP)
            CALL XABORT('XELLSR: BAD SURFACE IDENTIFICATION')
         ENDIF
  900 CONTINUE
C
 1000 RETURN
      END
