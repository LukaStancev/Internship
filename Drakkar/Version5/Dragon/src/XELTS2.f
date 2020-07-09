*DECK XELTS2
      SUBROUTINE XELTS2( IPRT, IFTEMP,NANGLE,DENUSR,NCODE,
     >                   ANGLES,DENSTY,SWZERO,
     >                   NTOTCL,MAXREM,REMESH,SUBMAX,LINMAX,
     >                   NSUR,NVOL,MATALB,INDEX,MINDIM,
     >                   MAXDIM,ICOORD,INCR,ICUR,TRKBEG,CONV,TRKDIR,
     >                   LENGHT,NUMERO,DDENWT)
************************************************************************
*                                                                      *
*           NAME: XELTS2                                               *
*      COMPONENT: EXCELL                                               *
*          LEVEL: 3 (CALLED BY 'XELTRK')                               *
*        VERSION: 1.0                                                  *
*       CREATION: 90/12                                                *
*       MODIFIED: 91/07 (R.R.)                                         *
*                 98/02 (G.M.) INCLUDE PERIODIC BOUNDARY OPTION        *
*                 00/03 (R.R.) DECLARE ALL VARIABLE TYPES              *
*         AUTHOR: ROBERT ROY                                           *
*                                                                      *
*     SUBROUTINE: THIS ROUTINE WILL CONSTRUCT THE SEQUENTIAL TAPE      *
*                 THAT WILL CONTAIN TRACKS FOR SPECULAR B.C.           *
*                 2-D CALCULATION. TRACKS ARE STORED ON 'IFTEMP' FILE. *
*                 HERE, THE CYCLIC TRACKING IS USED.                   *
*                                                                      *
*     REFERENCE: 'A CYCLIC TRACKING PROCEDURE FOR CP CALCULATIONS      *
*                 IN 2-D LATTICES', R.ROY ET AL.,                      *
*                 CONF/ADVANCES IN MATH, COMP & REACTOR PHYSICS,       *
*                 PITTSBURGH/USA, V 1, P 2.2 4-1 (1991).               *
*                                                                      *
*--------+-------------- V A R I A B L E S -------------+--+-----------*
*  NAME  /                  DESCRIPTION                 /IO/MOD(DIMENS)*
*--------+----------------------------------------------+--+-----------*
* IPRT   / INTERMEDIATE PRINTING LEVEL FOR OUTPUT.      /I./INT        *
* IFTEMP / TRACKING FILE #.                             /I./INT        *
* NANGLE / # OF ANGLES USED IN THE TRACKING PROCESS.    /I./INT        *
* DENUSR / DENSITY OF TRACKS IN THE PLANE PERPENDICULAR /I./REL        *
*        /      TO THE TRACKING ANGLES.                 /  /           *
* NCODE  / TYPE OF BOUNDARY CONDITIONS.                 /I./INT(6)     *
* ANGLES / 3D ANGLE VALUES.                             /I./R(3*NANGLE)*
* DENSTY / DENSITY OF TRACKS ANGLE BY ANGLE.            /I./REL(NANGLE)*
* SWZERO / LOGICAL VALUE (.T.= USE 0 AND PI/2 ANGLES).  /I./LOG        *
* NTOTCL / # OF CYLINDRES OF A TYPE + 2.                /I./INT        *
* MAXREM / MAX NUMBER OF REAL MESH VALUES IN "REMESH".  /I./INT        *
* REMESH / REAL MESH VALUES (RECT/CYL).                 /I./REL(MAXREM)*
* SUBMAX / MAX. # OF SUBTRACKS IN A SINGLE TRACK.       /I./INT        *
* LINMAX / MAX. # OF TRACK SEGMENTS IN A SINGLE TRACK.  /I./INT        *
*   NSUR / # OF SURFACES.                               /I./INT        *
*   NVOL / # OF ZONES.                                  /I./INT        *
* MATALB / MATERIAL TYPES (FACES FOR SURFACES).         /I./INT(  NVS) *
*  INDEX / #ING OF SURFACES & ZONES.                    /I./INT(4*NVS) *
* MINDIM / MIN INDEX VALUES FOR ALL AXES (RECT/CYL).    /I./INT(NTOTCL)*
* MAXDIM / MAX INDEX VALUES FOR ALL AXES (RECT/CYL).    /I./INT(NTOTCL)*
*   ICUR / CURRENT ZONAL LOCATION FOR A TRACK SEGMENT.  /../INT(NTOTCL)*
*   INCR / INCREMENT DIRECTION FOR NEXT TRACK SEGMENT.  /../INT(NTOTCL)*
* TRKBEG / POSITION WHERE A TRACK BEGINS.               /../REL(NTOTCL)*
*   CONV / SEGMENTS OF TRACKS.                          /../REL(NTOTCL)*
* TRKDIR / DIRECTION OF A TRACK IN ALL AXES.            /../REL(NTOTCL)*
* LENGHT / RELATIVE LENGHT OF EACH SEGMENT IN A TRACK.  /../REL(LINMAX)*
* NUMERO / MATERIAL IDENTIFICATION OF EACH TRACK SEGMENT/../INT(LINMAX)*
* DDENWT / DENSITY OF TRACKS ANGLE BY ANGLE.            /I./D(NANGLE)  *
*--------+---------------- R O U T I N E S -------------+--+-----------*
*  NAME  /                  DESCRIPTION                                *
*--------+-------------------------------------------------------------*
* XELLSR / TO FIND BEGINNING/ENDING SURFACES OF A TRACK.               *
* XELLIN / TO COMPUTE TRACKS IN GEOMETRIES.                            *
* XELTSA / TO GET THE ANGLES IN PERIODIC TRACKING.                     *
* XELTSW / TO GET THE WEIGHTS IN PERIODIC TRACKING.                    *
************************************************************************
C
      IMPLICIT NONE
C
      INTEGER           IPRT, IFTEMP, NANGLE, NTOTCL, MAXREM, SUBMAX,
     >                  LINMAX, NSUR,   NVOL
      REAL              TRKBEG(NTOTCL), TRKDIR(NTOTCL), CONV(NTOTCL),
     >                  REMESH(MAXREM),DENUSR
      DOUBLE PRECISION  DENSTY(4*NANGLE),ANGLES(2,4*NANGLE),
     >                  LENGHT(LINMAX),DDENWT(NANGLE)
      INTEGER           MINDIM(NTOTCL),MAXDIM(NTOTCL),ICUR(NTOTCL),
     >                  ICOORD(NTOTCL),INCR(NTOTCL),NUMERO(LINMAX),
     >                  NCODE(6),MATALB(NSUR:NVOL),INDEX(4,*)
C
      REAL              TRKEND(2), PROJC2(3), TRKCUT(3,2),
     >                  TRKPTS(3), OLDBEG(2), OLDDIR(2), 
     >                  EPS, TOTLEN, ZERO, ONE
      DOUBLE PRECISION  WEIGHT, ANGTSA(2,2), ANGLE2(2), ABSC(2), DENS,
     >                  DP, RCIRC, PROJ, PMAX, PMIN, DEPART, DENLIN,
     >                  DRKORI(2), RONEPS
      INTEGER           IOUT, NSCUT(2),INDC(2),IPER(2),IREFL(2)
      INTEGER           NDIM, IPERG, IDEB, ISUM, ISTRID, IANG, ITX, ITY,
     >                  IDIM, NOTRAK, NANGLS, NSOLMX, IREF1, ITG,
     >                  NSCAN, ISCAN, NTRAC, NPOINT,
     >                  NDEBS, IX, IY, I2, NSGANG, NTTRK, LINACT,
     >                  LINNUS, LINUSD, NCROS, JINT, KINT, II, IZZ,
     >                  NUMANG, I, J, K, LINE, NSUB
      LOGICAL           SWZERO, SWZDIR(3), SWBNEW
      CHARACTER         TEDATA*13
      PARAMETER       ( EPS=1.E-5, ZERO= 0.0E0, ONE=1.0E0, IOUT=6 ) 
C----
C  ALLOCATABLE ARRAYS
C----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: KANGL
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: PTSANG,WGTANG,
     > DNSANG
C----
C  SCRATCH STORAGE ALLOCATION
C   PTSANG: cosines of angles
C   WGTANG: weights of angles
C   DNSANG: densities of each angle
C----
      ALLOCATE(PTSANG(NANGLE),WGTANG(NANGLE),DNSANG(NANGLE))
C
      NDIM= 2
C
C     SET FLAG FOR SURFACE CROSSING
C       IPER(1) = X-PERIOD
C       IPER(2) = Y-PERIOD
C     VALUES ARE
C       IPER(I) = 1 FOR PERIODIC BC
C       IPER(I) = 2 FOR OTHER BC
      IPER(1)=2
      IPER(2)=2
      IF( NCODE(1).EQ.4 .AND. NCODE(2).EQ.4 ) THEN
         IPER(1)=1
      ENDIF
      IF( NCODE(3).EQ.4 .AND. NCODE(4).EQ.4 ) THEN
         IPER(2)=1
      ENDIF
      IPERG=   MIN(IPER(1)*IPER(2),2)
      ABSC(1)= DBLE(REMESH(MAXDIM(1)))-DBLE(REMESH(MINDIM(1)))
      ABSC(2)= DBLE(REMESH(MAXDIM(2)))-DBLE(REMESH(MINDIM(2)))
      RCIRC=   SQRT(ABSC(1)**2 + ABSC(2)**2)
      ABSC(1)= ABSC(1)/RCIRC
      ABSC(2)= ABSC(2)/RCIRC
C
C     SET ITX THE NUMBER OF X CROSSING
C     SET ITY THE NUMBER OF Y CROSSING
C     FOR STANDARD TRACKING INCLUDING POINTS AT 0.0 AND PI/2
C       SWZERO =.TRUE.
C       SCAN ITX FROM 0 TO NANGLE-1
C       SCAN ITY FROM 0 TO NANGLE-1
C       ONLY VALID VALUE IS ITX+ITY = NANGLE
C     FOR MEDI TRACKING EXCLUDING POINT AT O.O AND PI/2
C       SWZERO = .FALSE.
C       SCAN ITX FROM 1 TO 2+NANGLE BY STEPS OF 2
C       SCAN ITY FROM 1 TO 2*NANGLE BY STEPS OF 2
C       ONLY VALID VALUES IS ITX+ITY=2*NANGLE
      IF( SWZERO )THEN
         ISUM= NANGLE-1
         IDEB= 0
         ISTRID=1
      ELSE
         ISUM= 2*NANGLE
         IDEB= 1
         ISTRID=2
      ENDIF
      ALLOCATE(KANGL(SUBMAX))
C
C     FIRST ANGLE INITIALIZATION FOR STORING ON TRACKING FILE
C        ANGTSA(1,1)=  COS(THETA) WRT X-DIRECTION
C        ANGTSA(1,2)=  SIN(THETA) WRT X-DIRECTION
C        ANGTSA(2,1)=  SIN(THETA) WRT X-DIRECTION
C                   =  COS(PI/2-THETA) WRT Y-DIRECTION
C        ANGTSA(2,2)= -COS(THETA) WRT Y-DIRECTION
C                   =  SIN(PI/2-THETA) WRT Y-DIRECTION
C     1) GET SUCCESSIVE ANGLES COSINE USING XELTSA
C        RANGE 0 <= THETA <= PI/2
C        CAN BE EXTENDED TO 0 <= THETA <= PI
C        USING CHANGE OF SIGN FOR ANGTSA(1,1)
C     2) COMPUTE ALL ANGULAR INTEGRATION WEIGHTS USING XELTSW
C     3) STORE ON TRACKING FILE
      IANG= 0
      DO 80 ITX= IDEB, ISUM, ISTRID
         INDC(1)= ITX
         ITY=ISUM-ITX
         INDC(2)= ITY
C
C        READ ANGLE BY ANGLE
         CALL XELTSA( NDIM, ABSC, INDC, DENS, ANGTSA)
         IANG= IANG+1
         DENLIN= DENS / RCIRC
C
C        FOR HORIZONTAL & VERTICAL ANGLES
C           TRAK DENSITY = ORIGINAL DENSITY
C        OTHERWISE
C           FIND RATIO BETWEEN ORIGINAL DENSITY AND MINIMUM DENSITY
C           TRACK DENSITY =  CLOSEST MULTIPLE OF MINIMUM TRACK DENSITY
         IF( ITX.EQ.0 .OR. ITY.EQ.0 )THEN
            DENLIN= DBLE(DENUSR)
            DNSANG(IANG)= DENLIN
         ELSE
            NTRAC= MAX(1,INT(DBLE(DENUSR)/DENLIN+0.5D0))
            DNSANG(IANG)= DBLE(NTRAC) * DENLIN
         ENDIF
         PTSANG(IANG)=   REAL(ANGTSA(1,1))
         ANGLES(1,IANG)= REAL(ANGTSA(1,1))
         ANGLES(2,IANG)= REAL(ANGTSA(2,1))
  80  CONTINUE
C
C     COMPUTE ALL ANGULAR INTEGRATION WEIGHTS
      CALL XELTSW( ABSC, NANGLE, PTSANG, WGTANG)
      DO 90 IANG= 1, NANGLE
         DENSTY(IANG)= 2.0/REAL(WGTANG(IANG))
         IF( IPRT.GT.2 )THEN
            WRITE(IOUT,1000) IANG, PTSANG(IANG), WGTANG(IANG),
     >                     DNSANG(IANG), WGTANG(IANG)/DNSANG(IANG)
         ENDIF
  90  CONTINUE
      DO 100 IANG=1,NANGLE
        ANGLES(1,2*NANGLE-IANG+1)=-ANGLES(1,IANG)
        ANGLES(2,2*NANGLE-IANG+1)=ANGLES(2,IANG)
        DENSTY(2*NANGLE-IANG+1)=DENSTY(IANG)
        DDENWT(2*NANGLE-IANG+1)=0.25D0*DBLE(WGTANG(IANG)/DNSANG(IANG))
 100  CONTINUE
      DO 110 IANG=1,2*NANGLE
        ANGLES(1,4*NANGLE-IANG+1)=ANGLES(1,IANG)
        ANGLES(2,4*NANGLE-IANG+1)=-ANGLES(2,IANG)
        DENSTY(4*NANGLE-IANG+1)=DENSTY(IANG)
        DDENWT(4*NANGLE-IANG+1)=0.25D0*DBLE(WGTANG(IANG)/DNSANG(IANG))
 110  CONTINUE
C
C     COPY ANGLES AND DENSITIES ON TEMPORARY TRACKING FILE
      WRITE(IFTEMP) ((ANGLES(IDIM,IANG),IDIM=1,NDIM),IANG=1,4*NANGLE)
      WRITE(IFTEMP)  (DENSTY(IANG)                  ,IANG=1,4*NANGLE)
C
C     PREPARE FOR TRACKING
      PROJC2(1)= ZERO
      PROJC2(2)= ZERO
      PROJC2(3)= ONE
      TRKBEG(3)= ZERO
      TRKDIR(3)= ZERO
      NOTRAK= 0
      NANGLS=NANGLE
      NSOLMX= 0
      IF( IPRT.GT.1 )THEN
         WRITE(IOUT,'(1H )')
         WRITE(IOUT,1001) NANGLE
         NSOLMX= MIN(9, NANGLE/10)
         IREF1=0
         WRITE(IOUT,1002) (IREF1, IZZ=0,NSOLMX)
         WRITE(IOUT,1002) (MOD(IZZ,10), IZZ=0,NSOLMX)
         TEDATA= '(1H+,TXXX,I1)'
      ENDIF
C
C     READ SUCCESSIVE ANGLES COSINE USING XELTSA FOR TRACKING
      IANG= 0
      DO 120 ITX= IDEB, ISUM, ISTRID
         INDC(1)= ITX
         ITY=ISUM-ITX
         INDC(2)= ITY
C
C        READ ANGLE BY ANGLE
         CALL XELTSA( NDIM, ABSC, INDC, DENS, ANGTSA)
         IANG= IANG+1
C
C        COMPUTE NUMBER OF SEGMENTS OF THE ANGLE
C          ITX = ISUM -> ANGLE = 0
C          ITY = ISUM -> ANGLE = PI/2
         NUMANG=1
         NSCAN = IPERG
         IF( ITX.EQ.ISUM )THEN
            NSCAN = IPER(1)
         ELSEIF( ITY.EQ.ISUM )THEN
            NSCAN = IPER(2)
         ELSE
            NUMANG= ISUM
            DO 130 ITG= MIN(ITX,ITY), 2, -1
               IF( (ITX .EQ. (ITX/ITG)*ITG) .AND.
     >             (ITY .EQ. (ITY/ITG)*ITG)       )THEN
                 NUMANG=NUMANG/ITG
                 GO TO 135
               ENDIF
 130        CONTINUE
         ENDIF
 135     CONTINUE
         NUMANG=NUMANG*NSCAN
C
C        IF NSCAN = 2
C           0 TO PI/2 AND PI/2 TO PI ARE SCANNED SIMULTANEOUSLY
C        IF NSCAN= 1
C           FIRST TREAT 0 TO PI/2
C           THEN TREAT PI/2 TO PI
         DO 140 ISCAN=2,NSCAN,-1
C
C         ZEROS ON THE COMPONENT OF THE DIRECTION
          SWZDIR(1)= ITX.EQ.0
          SWZDIR(2)= ITY.EQ.0
          DENLIN = DNSANG(IANG)
          DP   = 1.D0 / DENLIN
          NDEBS= 0
          IF( (IPRT .GT. 1) .AND. (ISCAN .EQ. 2) ) THEN
            IF( MOD(IANG,100) .EQ. 0 )THEN
              IREF1=IREF1+1
              NDEBS= NSOLMX+1
              NSOLMX=MIN(NDEBS+9, NANGLE/10)
              WRITE(IOUT,1002)(IREF1,IZZ=NDEBS,NSOLMX)
              WRITE(IOUT,1002)(MOD(IZZ,10),IZZ=NDEBS,NSOLMX)
            ELSE
              IF( (IPRT.GT.10000) .AND. (MOD(IANG,100).NE.0) )THEN
                WRITE(IOUT,1002) (IREF1,IZZ=NDEBS,NSOLMX)
                WRITE(IOUT,1002) (MOD(IZZ,10),IZZ=NDEBS,NSOLMX)
              ENDIF
              WRITE(TEDATA(7:9),'(I3.3)') MOD(IANG,100) + 2
              WRITE(IOUT,TEDATA) MOD(IANG,10)
            ENDIF
          ENDIF
          DO 150 I   = 1, 2
            TRKDIR(I)= REAL(ANGTSA(I,1))
            INCR(I)  =  1
            IF( SWZDIR(I) ) INCR(I)= 0
 150      CONTINUE
C
C         CUT LENGHT FOR UNIT PLANE VECTORS
C         PROJECT THE 4 CORNERS ON THE PLANE
C
C         DIRECTION OF TRACK IS IN O TO PI/2 REPRESENTS ROTATED +X-AXIS
C         DIRECTION TRACK NORMAL IS IN -PI/2 TO 0 REPRESENTED   +Y-AXIS
C           PMIN IS LOWEST POINT FROM PROJECTING CARTESIAN CELL
C           IN THIS FRAME OF REFERENCE OF TRACKING LINE
C           PMAX IS HIGHEST POINT FROM PROJECTING CARTESIAN CELL
C           IN THIS FRAME OF REFERENCE OF TRACKING LINE
C           STARTING POINT IN THIS FRAME OF REFERENCE OF TRACKING LINE I
C              PMIN-0.5*(MESH SPACING DP)
C           STARTING POINT IN ORIGINAL FRAME OF REFERENCE IS
C              X=PMIN*ANGTSA(1,2)
C              Y=PMIN*ANGTSA(2,2)
C           MESH SPACING IN THE ORIGINAL FRAME OF REFERENCE IS
C              DX=DP*ANGTSA(1,2)
C              DY=DP*ANGTSA(2,2)
          PMIN = +1.0D+50
          PMAX = -1.0D+50
          DO 160   IX  =MINDIM(1),MAXDIM(1),MAXDIM(1)-MINDIM(1)
            DO 161 IY  =MINDIM(2),MAXDIM(2),MAXDIM(2)-MINDIM(2)
              PROJ  = DBLE(REMESH(IX)) * ANGTSA(1,2)
     >              + DBLE(REMESH(IY)) * ANGTSA(2,2)
              PMIN=MIN(PMIN,PROJ)
              PMAX=MAX(PMAX,PROJ)
 161        CONTINUE
 160      CONTINUE
C
C         NEAREST INTEGER -1 OR +1 FOR SECURITY
          NPOINT     =NINT((PMAX-PMIN)*DENLIN)+1
          DEPART     =PMIN - 0.5D0 * DP
          DO 170 J   = 1, 2
            DRKORI(J)= DEPART    * ANGTSA(J,2)
            ANGLE2(J)= DP        * ANGTSA(J,2)
 170      CONTINUE
          IF(ISCAN .EQ. 1) THEN
            IF(ITX .EQ. 0) THEN
              TRKDIR(2)=-TRKDIR(2)
              INCR(2)  =-1
            ELSE IF(ITY.EQ.0) THEN
              TRKDIR(1)=-TRKDIR(1)
              INCR(1)  =-1
            ENDIF
          ELSE
          ENDIF
          SWBNEW= .TRUE.
          NSGANG= 0
          IREFL(1)=1
          IREFL(2)=1
          LINACT= 0
          DO 180 I2  = 1,NSCAN*NPOINT
            IF( SWBNEW )THEN
              NSUB=0
              NTTRK=NOTRAK+1
              LINACT= 1
              LINNUS= LINMAX
              IF(NSGANG .EQ. 0) THEN
                DO 190 J   = 1, 2
                  DRKORI(J)= DRKORI(J) + ANGLE2(J)
                  TRKPTS(J)= REAL(DRKORI(J))
 190            CONTINUE
              ENDIF
              IF(ISCAN .EQ. 1 .AND. ITX*ITY .NE. 0) THEN
C
C               LOCATE STARTUP POSITION ON SURFACES
C               IDENTICAL TO CASE WITH ISCAN=2
                TRKDIR(1)= REAL(ANGTSA(1,1))
                TRKDIR(2)= REAL(ANGTSA(2,1))
                CALL XELLSR( NDIM, NTOTCL,   NSUR, MAXREM, REMESH,
     >                      INDEX, MINDIM, MAXDIM, ICOORD,   ICUR,
     >                       INCR, TRKPTS, TRKDIR, TRKCUT,  NSCUT,
     >                      NCROS, TOTLEN)
                IF( NCROS .LT. 2 ) GO TO 185
                TRKPTS(1)=TRKCUT(1,1)
                TRKPTS(2)=TRKCUT(2,1)
                JINT = (1-MATALB(NSCUT(1)))/2
                TRKDIR(JINT)= -REAL(ANGTSA(JINT,1))
                INCR(JINT)  = -1
              ENDIF
            ENDIF
C
C           LOCATE EXTERNAL SURFACES CROSSED BY THIS TRACK
            CALL XELLSR( NDIM, NTOTCL, NSUR, MAXREM, REMESH,
     >                  INDEX, MINDIM, MAXDIM, ICOORD, ICUR, INCR,
     >                  TRKPTS, TRKDIR, TRKCUT, NSCUT, NCROS,
     >                  TOTLEN)
C
C           VALID TRACK ONLY IF 2 SURFACES ARE CROSSED
C           OTHERWISE DO NOT CONSIDER TRACK
            IF( NCROS .LT. 2 ) GO TO 185
            NSGANG= NSGANG+1
            DO 200 K= 1, NDIM
              TRKBEG(K)= TRKCUT(K,1)
              TRKEND(K)= TRKCUT(K,2)
 200        CONTINUE
            IF( SWBNEW )THEN
              DO 210 J= 1, 2
                OLDBEG(J)= TRKBEG(J)
                OLDDIR(J)= TRKDIR(J)
 210          CONTINUE
            ENDIF
C
C           SAVE INITIAL SURFACE CROSSED BY TRACK
C           SINCE THE INITIAL IS DOUBLED SET
C            LENGTH TO 0.5 TO TAKE THIS EFFECT INTO ACCOUNT
            LENGHT(LINACT)= 0.5D0
            NUMERO(LINACT)= NSCUT(1)
            LINACT= LINACT + 1
C
C           LOCATE ALL REGIONS CROSSED BY LINE
            NSUB=NSUB+1
            IF(NSUB.GT.SUBMAX) CALL XABORT('XELTS2: SUBMAX OVERFLOW.')
            KANGL(NSUB)=0
            DO II=1,4*NANGLE
              IF((DBLE(TRKDIR(1)).EQ.ANGLES(1,II)).AND.
     >           (DBLE(TRKDIR(2)).EQ.ANGLES(2,II))) THEN
                  KANGL(NSUB)=II
                  GO TO 215
              ENDIF
            ENDDO
            CALL XABORT('XELTS2: UNABLE TO FIND AN ANGULAR INDEX FOR A'
     >      //' SUBTRACK')
 215        CALL XELLIN(  NDIM, NTOTCL, MAXREM, REMESH,
     >                    NSUR, NVOL, INDEX, MINDIM, MAXDIM,
     >                  ICOORD, ICUR, INCR, TRKBEG, TRKEND, TRKDIR,
     >                  PROJC2, TOTLEN,
     >                    CONV, LINNUS, LENGHT(LINACT), NUMERO(LINACT),
     >                  LINUSD)
            LINACT= LINACT + LINUSD
            LINNUS= LINNUS - LINUSD
C
C           SAVE FINAL SURFACE CROSSED BY TRACK
C           SINCE THE FINAL SURFACES IS DOUBLED SET
C            LENGTH TO 0.5 TO TAKE THIS EFFECT INTO ACCOUNT
            LENGHT(LINACT)= 0.5D0
            NUMERO(LINACT)= NSCUT(2)
            LINACT= LINACT + 1
            LINNUS= LINNUS - 1
C
C           FIND INTERSECTION DIRECTION
            JINT = (1-MATALB(NSCUT(2)))/2
            KINT = MOD(JINT,2) +1
C
C           IF NCODE(J)=4
C             -> TRANSLATION FOR THE FACE
C                FOR LOWER FACE (TRKBEG(J)=REMESH(MINDIM(J)))
C                  TRKBEG -> REMESH(MAXDIM(J))
C                FOR UPPER FACE (TRKBEG(J)=REMESH(MAXDIM(J)))
C                  TRKBEG -> REMESH(MINDIM(J))
C           OTHERWISE
C             -> SPECULAR REFLECTION FOR THE FACE
C                  TRKDIR(J)= -TRKDIR(J)
C                  INCR(J)=-INCR(J)
            IF( IPER(JINT).EQ.1 )THEN
              IF( TRKEND(JINT).EQ.REMESH(MAXDIM(JINT)) )THEN
                TRKEND(JINT)= REMESH(MINDIM(JINT))
              ELSEIF( TRKEND(JINT).EQ.REMESH(MINDIM(JINT)) )THEN
                TRKEND(JINT)= REMESH(MAXDIM(JINT))
              ELSE
                CALL XABORT('XELTS2: TRANSLATION ERROR')
              ENDIF
            ELSE
              TRKDIR(JINT)= -TRKDIR(JINT)
              INCR(JINT)=   -INCR(JINT)
            ENDIF
            RONEPS= ZERO
            DO 220 J= 1, 2
               TRKBEG(J)= TRKEND(J)
               RONEPS= RONEPS + (TRKBEG(J)-OLDBEG(J))**2
     >                        + (TRKDIR(J)-OLDDIR(J))**2
               TRKPTS(J)= TRKBEG(J)
 220        CONTINUE
            RONEPS=RONEPS/(RCIRC*RCIRC)
            SWBNEW= NSGANG.EQ.NUMANG
            IF(SWBNEW) THEN
              NSGANG=0
              IF(RONEPS .GT. EPS) THEN
                WRITE(IOUT,9000) (OLDBEG(J),J=1,2),(TRKBEG(J),J=1,2),
     >                           (OLDDIR(J),J=1,2),(TRKDIR(J),J=1,2),
     >                            100.0*SQRT(RONEPS)
                CALL XABORT
     >            ('XELTS2: ERROR ON FINAL POSITION OR DIRECTION')
              ENDIF
C
C             RESET ORIGINAL ANGLES IF ROTATION CONSIDERED
              DO 230 II=1,2
                TRKDIR(II)= IREFL(II)*TRKDIR(II)
                INCR(II)  = IREFL(II)*INCR(II)
                IREFL(II) = 1
 230          CONTINUE
            ELSE IF(NSCAN .EQ. 2) THEN
              IF(RONEPS .LE. EPS) THEN
                IF(IPER(JINT) .EQ. 1) THEN
C
C                 ROTATE ANGLE IF STARTUP SURFACE PERIODIC
                  TRKDIR(JINT)= -TRKDIR(JINT)
                  INCR(JINT)=   -INCR(JINT)
                  IREFL(JINT)=-IREFL(JINT)
                  SWBNEW= .TRUE.
                ELSE IF(IPER(KINT) .EQ. 1) THEN
C
C                 IF NORMAL SURFACE IS PERIODIC
C                 LOCATE FIRST NORMAL SURFACE REACHED AND ROTATE ANGLE
                  DO 240 II=1,NSCAN*NPOINT
                    CALL XELLSR( NDIM, NTOTCL,   NSUR, MAXREM, REMESH,
     >                          INDEX, MINDIM, MAXDIM, ICOORD,   ICUR,
     >                           INCR, TRKPTS, TRKDIR, TRKCUT,  NSCUT,
     >                          NCROS, TOTLEN)
                    JINT = (1-MATALB(NSCUT(2)))/2
                    IF(JINT.EQ.KINT) GO TO 245
 240              CONTINUE
                  CALL XABORT
     >              ('XELTS2: CANNOT FIND AN INITIAL PERIODIC SURFACE')
 245              CONTINUE
                  TRKPTS(1)=TRKCUT(1,2)
                  TRKPTS(2)=TRKCUT(2,2)
                  TRKDIR(JINT)= -TRKDIR(JINT)
                  INCR(JINT)  = -INCR(JINT)
                  IREFL(JINT)=-IREFL(JINT)
                  SWBNEW= .TRUE.
                ENDIF
              ENDIF
            ENDIF
            IF( SWBNEW  )THEN
C
C             NOW, WRITE THE TRACK
              NOTRAK= NOTRAK + 1
              NTTRK = NOTRAK
              LINE= LINACT-1
              WEIGHT= 0.25D0*DBLE(WGTANG(IANG)/DNSANG(IANG))
              WRITE(IFTEMP) NSUB, LINE, WEIGHT,
     >                     (KANGL(II),II=1,NSUB),
     >                     (NUMERO(I),I=1,LINE),
     >                     (LENGHT(I),I=1,LINE)
              IF( IPRT.GT.    1000 ) THEN
                WRITE(IOUT,1010) NOTRAK,IANG,LINE
                WRITE(IOUT,1011) (LENGHT(I),NUMERO(I),I=1,LINE)
              ENDIF
            ENDIF
 185        CONTINUE
 180      CONTINUE
 140    CONTINUE
 120  CONTINUE
C
      IF( IPRT.GT.1 )THEN
         WRITE(IOUT,'(1H )')
         WRITE(IOUT,1020)  IFTEMP,NANGLE, DENUSR,NOTRAK
      ENDIF
C----
C  SCRATCH STORAGE DEALLOCATION
C----
      DEALLOCATE(KANGL,PTSANG,WGTANG,DNSANG)
      RETURN
C
 1000 FORMAT(1X,I4,': COS=',F15.10,' WGT=',F15.10,' DNS=',F15.10,
     >             ' WGT/DEN=',F15.10)
 1001 FORMAT(1X,'ECHO = ',I10,' SOLID ANGLES TO BE TRACKED ')
 1002 FORMAT(1X,10(I1,9X))
 1010 FORMAT(1X,'#',I10,' IANG=',I10,' LEN=',I10)
 1011 FORMAT(1P,(1X,E15.3,1X,I10))
 1020 FORMAT(1X,'ECHO OF TRACKING PROPERTIES FOR TAPE ',I2/
     >       10X,'NUMBER OF ANGLES       =',I10/
     >       10X,'   TRACK DENSITY       =',F10.3,' LINES/CM'/
     >       10X,'NUMBER OF TRACK STORED =',I10)
C
 9000 FORMAT(1X,'FINAL TRACK POSITION EXPECTED  = (',
     >       F15.10,',',F15.10,')'/
     >       1X,'FINAL TRACK POSITION FOUND     = (',
     >       F15.10,',',F15.10,')'/
     >       1X,'FINAL TRACK DIRECTION EXPECTED = (',
     >       F15.10,',',F15.10,')'/
     >       1X,'FINAL TRACK DIRECTION FOUND    = (',
     >       F15.10,',',F15.10,')'/
     >       1X,'RMS ERROR = ',F15.10,' %')
      END
