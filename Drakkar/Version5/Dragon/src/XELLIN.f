*DECK XELLIN
      SUBROUTINE XELLIN(NDIM,NCP,MAXREM,REMESH,
     >                  NSUR,NVOL,INDEL,MINDIM,MAXDIM,
     >                  ICOORD,ICUR,INCR,TRKBEG,TRKEND,TRKDIR,
     >                  PROJC2,TOTLEN,
     >                  CONV,LINMAX,LENGHT,NUMERO,LINE)
************************************************************************
*                                                                      *
*           NAME: XELLIN                                               *
*      COMPONENT: EXCELL                                               *
*          LEVEL: 4 (CALLED BY 'XELTI2' & 'XELTI3' & 'XELTS2' )        *
*        VERSION: 1.0                                                  *
*       CREATION: 87/01                                                *
*       MODIFIED: 91/07 (R.R.)                                         *
*                 00/03 (R.R.) DECLARE ALL VARIABLE TYPES              *
*         AUTHOR: ROBERT ROY                                           *
*                                                                      *
*     SUBROUTINE: THIS ROUTINE WILL CONSTRUCT ONE TRACKING LINE        *
*                 WHICH CONSISTS OF TWO VECTORS:                       *
*                          (LENGHT(I),I=1,LINE) <= SEGMENT LENGTHS     *
*                          (NUMERO(I),I=1,LINE) <= REGION  NUMBERS     *
*                                                                      *
*--------+-------------- V A R I A B L E S -------------+--+-----------*
*  NAME  /                  DESCRIPTION                 /IO/MOD(DIMENS)*
*--------+----------------------------------------------+--+-----------*
* NDIM   / # OF DIMENSION (2 OR 3).                     /I./INT        *
* NCP    / # OF CYLINDRES OF A TYPE + 3    (.LT.NC3MAX)./I./INT        *
* MAXREM / MAX NUMBER OF REAL MESH VALUES IN "REMESH".  /I./INT        *
* REMESH / REAL MESH VALUES (RECT/CYL).                 /I./REL(MAXREM)*
* NSUR   / # OF SURFACES.                               /I./INT        *
* NVOL   / # OF ZONES.                                  /I./INT        *
* INDEL  / #ING OF SURFACES & ZONES.                    /I./INT(4*NVS) *
* MINDIM / MIN INDEX VALUES FOR ALL AXES (RECT/CYL).    /I./INT(NC3MAX)*
* MAXDIM / MAX INDEX VALUES FOR ALL AXES (RECT/CYL).    /I./INT(NC3MAX)*
* ICOORD / PRINCIPAL AXES DIRECTION (X/Y/Z) FOR MESHES. /I./INT(NC3MAX)*
* ICUR   / CURRENT ZONAL LOCATION FOR A TRACK SEGMENT.  /../INT(NC3MAX)*
* INCR   / INCREMENT DIRECTION FOR NEXT TRACK SEGMENT.  /../INT(NC3MAX)*
* TRKBEG / POSITION WHERE A TRACK BEGINS.               /../REL(NC3MAX)*
* TRKEND / POSITION WHERE A TRACK ENDS.                 /../REL(NC3MAX)*
* TRKDIR / DIRECTION OF A TRACK IN ALL AXES.            /I./REL(NC3MAX)*
* PROJC2 / PROJECTIONS OF "TRKDIR" ALONG TRACKED ANGLES./I./REL(3)     *
* TOTLEN / TOTAL LENGHT OF THE TRACK.                   /I./REL        *
* CONV   / SEGMENTS OF TRACKS.                          /../REL(NC3MAX)*
* LINMAX / MAX. # OF TRACK SEGMENTS IN A SINGLE TRACK.  /I./INT        *
* LENGHT / RELATIVE LENGHT OF EACH SEGMENT IN A TRACK.  /.O/REL(LINMAX)*
* NUMERO / REGION IDENTIFICATION OF EACH TRACK SEGMENT  /.O/INT(LINMAX)*
* LINE   / LENGHT OF THE TRACK                          /.O/INT        *
************************************************************************
C
      IMPLICIT          NONE
C
      INTEGER           NDIM, NCP, MAXREM, NSUR, NVOL, LINMAX, LINE
      REAL              TRKBEG(NCP), TRKDIR(NCP), CONV(NCP),
     >                  REMESH(MAXREM), TRKEND(*), PROJC2(*), TOTLEN
      DOUBLE PRECISION  LENGHT(LINMAX)
      INTEGER           MINDIM(NCP), MAXDIM(NCP), ICUR(NCP),
     >                  ICOORD(NCP), INCR(NCP), INDEL(4,*),
     >                  NUMERO(LINMAX)
C
      INTEGER           IORD(4), N, NP1, NP2, IBEGIN, IEND, KELVOL
      REAL              TKBEG1, TKBEG2, TKEND1, TKEND2, R2BEG, R2END
      DOUBLE PRECISION  CONVOK, PAT0, PAT1
      LOGICAL           BETWEN
      INTEGER           NEXT, NUM, I, J
      REAL              ANORM2, CENTRE, A, B
C
      ANORM2(A,B)= A*A + B*B
      CENTRE(I,J)= REMESH( MAXDIM(I-1) + J )
      NEXT(J)= ICUR(J) + MAX( 0, INCR(J) )
      NUM(J)= J + 1 - NSUR
C
C     IF THERE ARE NO CYLINDER AT ALL
      IEND=0
      DO 90 I   = 1, NDIM
         N     = ICOORD(I)
C
C        FIND BEGINNING VOLUME #
         ICUR(I)= MINDIM(I)
         DO 80 J = MINDIM(I), MAXDIM(I)-1
            IF(TRKBEG(N).GE.REMESH(J)) ICUR(I)= J
            IF(TRKEND(N).GE.REMESH(J)) IEND= J
   80    CONTINUE
         IF( ICUR(I).EQ.IEND )THEN
            CONV(I)= TOTLEN
         ELSE
            IF( INCR(N).EQ.0 )THEN
               CONV(I)= TOTLEN
            ELSE
               CONV(I)=(REMESH(NEXT(I))-TRKBEG(N))/TRKDIR(N)
            ENDIF
         ENDIF
   90 CONTINUE
C
      IBEGIN= MAXDIM(3) + 3
      DO 110 I  = 4, NCP
         N     = ICOORD(I)
         NP1   = MOD(N  ,3) + 1
         NP2   = MOD(N+1,3) + 1
         TKBEG1= CENTRE(I,1) - TRKBEG(NP1)
         TKEND1= CENTRE(I,1) - TRKEND(NP1)
         TKBEG2= CENTRE(I,2) - TRKBEG(NP2)
         TKEND2= CENTRE(I,2) - TRKEND(NP2)
         R2BEG = ANORM2(TKBEG1,TKBEG2)
         R2END = ANORM2(TKEND1,TKEND2)
         TRKBEG(I)= (TKBEG1*TRKDIR(NP1)+TKBEG2*TRKDIR(NP2))/PROJC2(N)
         BETWEN   =  0.0 .LT. TRKBEG(I) .AND. TRKBEG(I) .LT. TOTLEN
         TRKDIR(I)=  R2BEG - TRKBEG(I) * TRKBEG(I) * PROJC2(N)
         ICUR(I)  = IBEGIN - 1
         MINDIM(I)= IBEGIN - 1
         IEND     = IBEGIN - 1
         DO 100 J  = IBEGIN, MAXDIM(I)
            IF( R2BEG    .GE. REMESH(J) )ICUR(I)= J
            IF( TRKDIR(I).GE. REMESH(J) )MINDIM(I)= J
            IF( R2END    .GE. REMESH(J) )IEND   = J
  100    CONTINUE
         IBEGIN= MAXDIM(I) + 3
         IF( ICUR(I).EQ.MINDIM(I) .AND.
     >       IEND   .EQ.MINDIM(I)       )THEN
            CONV(I)=TOTLEN
         ELSE
            IF( (BETWEN .AND. ICUR(I).NE.MINDIM(I)) .OR.
     >                        ICUR(I).GT.IEND            )THEN
               INCR(I)=-1
            ELSE
               INCR(I)=+1
            ENDIF
            IF( NEXT(I).GT.MAXDIM(I) )THEN
               CONV(I)=TOTLEN
            ELSE
               CONV(I)= TRKBEG(I) + INCR(I) *
     >                              SQRT((REMESH(NEXT(I))-TRKDIR(I))
     >                            / PROJC2(ICOORD(I)))
            ENDIF
         ENDIF
  110 CONTINUE
C
C     VOLUME TRACKED
      LINE  = 0
  120 LINE  = LINE + 1
C
C        LOOKING FOR THE MINIMUM VALUE IN "CONVOK"
         CONVOK= TOTLEN
         DO 130 I= 1, NCP
            IF( CONV(I) .LT. CONVOK ) CONVOK= CONV(I)
  130    CONTINUE
         DO 135 I= 1, NCP
            IORD(MIN(4,I))= ICUR(I)
            IF(I.GT.3.AND.ICUR(I).LT.MAXDIM(I)) GOTO 136
  135    CONTINUE
         IORD(4)= 0
  136    CONTINUE
         KELVOL= NVOL
         INDEL(1,NUM(0))= IORD(1)
         INDEL(2,NUM(0))= IORD(2)
         INDEL(3,NUM(0))= IORD(3)
         INDEL(4,NUM(0))= IORD(4)
  885    CONTINUE
            IF( IORD(1).EQ.INDEL(1,NUM(KELVOL)).AND.
     >          IORD(2).EQ.INDEL(2,NUM(KELVOL)).AND.
     >          IORD(3).EQ.INDEL(3,NUM(KELVOL)).AND.
     >          IORD(4).EQ.INDEL(4,NUM(KELVOL)) ) GO TO 895
            KELVOL= KELVOL - 1
         GO TO 885
  895    IF(KELVOL.EQ.0) CALL XABORT('XELLIN: TRACKING FAILURE.')
         NUMERO(LINE)= KELVOL
         LENGHT(LINE)= CONVOK
C
C        IF "CONVOK" IS "TOTLEN" THE TRACKING IS FINISHED
         IF( CONVOK.EQ.TOTLEN ) GO TO 160
C
C        UPDATE WHERE THE MINIMUM VALUE "CONVOK" IS OBTAINED
         DO 140 I   = 1, NDIM
            IF( CONV(I) .NE. CONVOK ) GO TO 140
            ICUR(I)= ICUR(I) + INCR(I)
            IF( NEXT(I).GT.MAXDIM(I) .OR.
     >          NEXT(I).LT.MINDIM(I)      ) GO TO 160
            N      = ICOORD(I)
            CONV(I)= ( REMESH(NEXT(I)) - TRKBEG(N) ) / TRKDIR(N)
  140    CONTINUE
         DO 150 I   = 4, NCP
            IF( CONV(I) .NE. CONVOK ) GO TO 150
            ICUR(I)= ICUR(I) + INCR(I)
            IF( ICUR(I) .EQ. MINDIM(I) .AND.
     >          INCR(I) .EQ. -1              )THEN
               INCR(I)= +1
            ENDIF
            IF( NEXT(I) .GT. MAXDIM(I) .OR.
     >          NEXT(I) .LT. MINDIM(I)       )THEN
               CONV(I)= TOTLEN
            ELSE
               N = ICOORD(I)
               CONV(I)= TRKBEG(I) + INCR(I) *
     >                  SQRT( (REMESH(NEXT(I))-TRKDIR(I))/PROJC2(N) )
            ENDIF
  150    CONTINUE
C
C     GO TO NEXT COORDINATE
      IF( LINE .NE. LINMAX ) GO TO 120
      CALL XABORT('XELLIN: TOO MANY TRACKS')
  160 CONTINUE
C
C     TRANSFORM LOCAL COORDINATES TO PATH LENGTHS
      PAT0=   0.0D0
      DO 170 I= 1, LINE
         PAT1= LENGHT(I)
         LENGHT(I)= PAT1-PAT0
         PAT0= PAT1
  170 CONTINUE
C
      RETURN
      END
