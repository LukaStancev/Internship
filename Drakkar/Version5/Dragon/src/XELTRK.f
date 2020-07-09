*DECK XELTRK
      SUBROUTINE XELTRK(IPTRK ,IPGEOM,GEONAM,IDISP ,IFTEMP,
     >                  IPRT  ,NDIM  ,ITOPT ,NV    ,NS    ,NANGL ,
     >                  ISYMM ,DENUSR,RCUTOF,MXSUB ,MXSEG ,ICODE ,
     >                  TITREC,INSB  ,IZ    ,LPRISM,NPRISM)
************************************************************************
*                                                                      *
*           NAME: XELTRK                                               *
*      COMPONENT: EXCELL                                               *
*          LEVEL: 2 (CALLED BY 'EXCELT')                               *
*        VERSION: 1.0                                                  *
*       CREATION: 93/04 (R.R.)                                         *
*       MODIFIED: 97/11 (G.M.) INTRODUCE OFFCENTER VARIATION           *
*                 00/03 (R.R.) DECLARE ALL VARIABLE TYPES              *
*         AUTHOR: ROBERT ROY                                           *
*                                                                      *
*     SUBROUTINE: THIS ROUTINE IS USED TO TREAT CARTESIAN ASSEMBLIES   *
*                 OF CELLS; THIS IS A TWO-STEP PROCESS:                *
*                 1) STUDY THE GEOMETRY TO GET VOLUMES AND MATERIALS   *
*                 2) PRODUCE TEMPORARY TRACKING FILE IF NECESSARY      *
*                                                                      *
*--------+-------------- V A R I A B L E S -------------+--+-----------*
*  NAME  /                  DESCRIPTION                 /IO/MOD(DIMENS)*
*--------+----------------------------------------------+--+-----------*
* IPTRK  / POINTER TO THE TRACKING (L_TRACK)            /I./TYPE(C_PTR)*
* IPGEOM / POINTER TO THE GEOMETRY (L_GEOM)             /I./TYPE(C_PTR)*
* GEONAM / GEOMETRY NAME                                /I./CAR*12     *
* IDISP  / STATUS OF TRACKING FILE (NEWT:IDISP.GT.0)    /IO/INT        *
*        /  CHANGED ONLY FOR *HALT* OPTION.             /  /           *
* IFTEMP / UNIT NUMBER ALLOCATED TO TEMPORARY FILE      /I./INT        *
* IPRT   / GEOMETRY PRINT LEVEL                         /I./INT        *
* NDIM   / # OF DIMENSIONS (2.OR.3)                     /I./INT        *
* ITOPT  / KIND OF TRACKING (0:ISOTROPIC;1:SPECULAR)    /IO/INT        *
* NV     / # OF ZONES IN THE ASSEMBLY.                  /.O/INT        *
* NS     / # OF SURFACES IN THE ASSEMBLY.               /.O/INT        *
* NANGL  / # OF ANGLES USED IN TEMPORARY TRACKING FILE. /.O/INT        *
* ISYMM  / SYMMETRY FACTOR.                             /IO/INT        *
* DENUSR / DENSITY OF TRACKS IN THE PLANE PERPENDICULAR /IO/REL        *
*        /      TO THE TRACKING ANGLES.                 /  /           *
* RCUTOF / CUTOF FOR CORNER TRACKING( 0.25 SUGGESTED )  /I./REL        *
* MXSUB  / MAXIMUM # OF SUBTRACKS IN A SINGLE TRACK.    /.O/INT        *
* MXSEG  / MAXIMUM # OF SEGMENTS IN A SINGLE TRACK.     /.O/INT        *
* ICODE  / INDEX FOR BOUNDARY CONDITIONS.               /.O/INT(6)     *
* TITREC / TITLE FOR THIS CASE                          /I./C*72       *
* INSB   / CONTROL ON VECTORIZATION                     /I./INT        *
* IZ     / PROJECTION AXIS FOR 3D PRISMATIC GEOMETRY    /I./INT        *
* LPRISM / FLAG FOR 3D PRISMATIC GEOMETRY               /I./LOG        *
* NPRISM / NUMER OF PLANS FOR A 3D PRISMATIC GEOMETRY   /.O/INT        *
*--------+---------------- R O U T I N E S -------------+--+-----------*
*  NAME  /                  DESCRIPTION                                *
*--------+-------------------------------------------------------------*
* AXGXEL / TO ANALYZE GEOMETRY                                         *
* XELGPR / TO PRINT A GEOMETRY REPRESENTATION                          *
* XELTI2 / TO TRACK IN 2D GEOMETRIES (REGULAR TRACKS)                  *
* XELTI3 / TO TRACK IN 3D GEOMETRIES (REGULAR TRACKS)                  *
* XELTS2 / TO TRACK IN 2D GEOMETRIES (CYCLIC TRACKS)                   *
************************************************************************
C
      USE             GANLIB
      IMPLICIT        NONE
      TYPE(C_PTR)     IPTRK ,IPGEOM
      INTEGER         IDISP ,IFTEMP,IPRT  ,NDIM  ,ITOPT ,NV    ,NS    ,
     >                NANGL ,ISYMM ,MXSUB ,MXSEG ,INSB  ,IZ    ,NPRISM
      LOGICAL         LPRISM
      REAL            DENUSR,RCUTOF
      CHARACTER       GEONAM*12, TITREC*72
      INTEGER         ICODE(6)
C
      INTEGER         NSTATE, IOUT
      PARAMETER     ( NSTATE=40, IOUT=6)
C----
C  ALLOCATABLE ARRAYS
C----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: KEYMRG,MATALB,MINDIM,MAXDIM,
     > ICORD,INDEX,ICUR,INCR,NUMERO
      REAL, ALLOCATABLE, DIMENSION(:) :: VOLSUR,REMESH,CONV,TRKBEG,
     > TRKDIR
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: ANGLES,DENSTY,
     > LENGHT,DDENWT
C----
C  LOCAL PARAMETERS
C---- 
      INTEGER         LTRK  ,NANGLE,SUBMAX,LINMAX
      INTEGER         NSUR  ,NVOL  ,NTOTCL,MAXR  ,NUNK  ,NEXTGE
      INTEGER         NTX   ,NTY   ,NTZ   ,NTR   ,ICL   ,NC    ,
     >                NALBG ,JJ    ,INDLEC,ITYLCM,ILENGT,INDATA,
     >                NCOR       
      INTEGER         NCODE(6), LCLSYM(3), ISTATE(NSTATE)
      LOGICAL         SWZERO
      REAL            ALBEDO(6), EXTKOP(NSTATE), CUTOFX, REDATA
      DOUBLE PRECISION DBLINP
      CHARACTER       CTISO*8, CTSPC*8,  CCORN*8, CSYMM*8, CMEDI*8,
     >                CHALT*8, CBLAN*8, TEDATA*8, CTRK*4, COMENT*80
      INTEGER         MXANGL
C
      SAVE     CTISO, CTSPC, CCORN, CHALT, CSYMM, CMEDI, CBLAN
      DATA     CTISO, CTSPC, CCORN, CHALT, CSYMM, CMEDI, CBLAN
     >       / 'TISO','TSPC','CORN','HALT','SYMM','MEDI','    ' /
C
      SWZERO=.TRUE.
      CALL XDISET(ISTATE,NSTATE,0)
      CALL XDRSET(EXTKOP,NSTATE,0.0)
      CALL LCMLEN(IPTRK,'STATE-VECTOR',ILENGT,ITYLCM)
      IF(ILENGT .LE. 0 .OR. ILENGT .GT. NSTATE) THEN
         LTRK  = 0
         NANGLE= 0
         ISYMM=1
         DENUSR= 0.0
         RCUTOF= 0.0
      ELSE
         CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
         CALL LCMGET(IPTRK,'EXCELTRACKOP',EXTKOP)
         LTRK=ISTATE(9)+1
         NANGLE=ISTATE(11)
         ISYMM=ISTATE(12)
         DENUSR=EXTKOP(2)
         RCUTOF=EXTKOP(3)
      ENDIF
      CUTOFX= 0.0
C
C. 1) READ ALL USER INPUT.
C
C     READ TRACKING PARAMETERS  LTRK= 1 : ISOTROPIC TRACKING   (TISO)
C                               LTRK= 2 : SPECULAR TRACKING    (TSPC)
      TEDATA= CBLAN
      IF( IDISP.GT.0 )THEN
         CALL REDGET( INDLEC, NANGLE, DENUSR, TEDATA,DBLINP)
         IF( ILENGT.NE.0 )THEN
            IF( INDLEC.EQ.3.AND.TEDATA(1:4).EQ.';' )THEN
               GO TO 10
            ENDIF
         ENDIF
         IF( TEDATA(1:4).EQ.CCORN(1:4) ) THEN
            CALL REDGET( INDLEC, INDATA, RCUTOF, TEDATA,DBLINP)
            CALL REDGET( INDLEC, NANGLE, DENUSR, TEDATA,DBLINP)
         ENDIF
         IF( INDLEC.NE.3 )THEN
            LTRK  = 1
            CALL REDGET( INDLEC, NANGLE, DENUSR, TEDATA,DBLINP)
            IF( INDLEC.EQ.3 )THEN
            CALL XABORT('XELTRK: *TISO* ASSUMED, PUT NANGLE & DENSTY' )
            ENDIF
            CALL REDGET( INDLEC, INDATA, REDATA, TEDATA,DBLINP)
         ELSE
            IF( TEDATA(1:4).EQ.CTISO(1:4) )THEN
               LTRK  = 1
            ELSEIF( TEDATA(1:4).EQ.CTSPC(1:4) )THEN
               LTRK  = 2
            ENDIF
            IF( LTRK.GT.0 )THEN
               CALL REDGET( INDLEC, NANGLE, DENUSR, TEDATA,DBLINP)
               IF( LTRK.EQ.2.AND.TEDATA(1:4).EQ.CMEDI(1:4) )THEN
                  SWZERO= .FALSE.
                  CALL REDGET( INDLEC, NANGLE, DENUSR, TEDATA,DBLINP)
               ENDIF
               CALL REDGET( INDLEC, NANGLE, DENUSR, TEDATA,DBLINP)
               CALL REDGET( INDLEC, INDATA, REDATA, TEDATA,DBLINP)
            ELSE
               IF( TEDATA(1:4).EQ.CHALT(1:4) )THEN
                  CALL REDGET( INDLEC, INDATA, REDATA, TEDATA,DBLINP)
                  IDISP = -2
               ELSE
                  CALL XABORT( 'XELTRK: *TISO*,*TSPC*,*HALT* HERE')
               ENDIF
            ENDIF
         ENDIF
         IF( INDLEC.NE.3 )THEN
            CALL XABORT( 'XELTRK: ; SYMM or NOSYMM PERMITTED' )
         ELSEIF( TEDATA(1:4) .EQ. 'NOSY'   )THEN
            ISYMM=0
            CALL REDGET( INDLEC, INDATA, REDATA, TEDATA,DBLINP)
         ELSEIF( TEDATA(1:4).EQ.CSYMM(1:4) )THEN
            CALL REDGET( INDLEC, INDATA, REDATA, TEDATA,DBLINP)
            IF(INDLEC.NE.1) CALL XABORT('XELTRK: INTEGER DATA EXPECTED')
            ISYMM=INDATA
            CALL REDGET( INDLEC, INDATA, REDATA, TEDATA,DBLINP)
         ENDIF
         IF( INDLEC.NE.3.OR.TEDATA(1:4).NE.';' )THEN
            CALL XABORT( 'XELTRK: ; IS SUPPOSED TO BE HERE' )
         ENDIF
         IF( LTRK.GT.0.AND.NANGLE.LE.1 )THEN
          CALL XABORT( 'XELTRK: INVALID NUMBER OF ANGLES (NANGLE < 2)' )
         ENDIF
         IF( LTRK.GT.0.AND.DENUSR.LE.0.0 )THEN
          CALL XABORT( 'XELTRK: INVALID DENSITY (DENSTY < 0.0) ' )
         ENDIF
      ENDIF
   10 CONTINUE
C----
C  PROCESS THE GEOMETRY
C----
      CALL AXGXEL(IPGEOM,IPTRK ,IPRT  ,GEONAM)
      IF(IPRT .GE. 1)THEN
         WRITE(IOUT,'(1H )')
         WRITE(IOUT,'(26H >>> GEOMETRY TREATED:    ,A12)') GEONAM
         WRITE(IOUT,'(26H >>> EXCELL TREATMENT <<< )')
         WRITE(IOUT,'(1H )')
      ENDIF
C----
C  SAVE EXCELL SPECIFIC TRACKING INFORMATION.
C----
      ISTATE(9)=LTRK-1
      ISTATE(11)=NANGLE
      ISTATE(12)=ISYMM
      CALL LCMPUT(IPTRK,'STATE-VECTOR',NSTATE,1,ISTATE)
      EXTKOP(1)=CUTOFX
      EXTKOP(2)=DENUSR
      EXTKOP(3)=RCUTOF
      CALL LCMPUT(IPTRK,'EXCELTRACKOP',NSTATE,2,EXTKOP)
C----
C  IF A PRISMATIC 3D TRACKING IS REQUESTED, 
C  CREATE 2D PROJECTED GEOMETRY ANALYSIS
C----
      IF (LPRISM) THEN
         CALL XELPR3(IPTRK,IZ,NPRISM)
         CALL LCMSIX(IPTRK,'PROJECTION  ',1)
         CALL LCMSIX(IPTRK,'EXCELL      ',1)
      ELSE
         CALL LCMSIX(IPTRK,'EXCELL      ',1)         
      ENDIF
C----
C  ALLOCATE GEOMETRIC STRUCTURES (SEE COMMON CEXGEO)
C     KEYMRG   : INTEGER MERGE VECTOR
C     VOLSUR   : REAL VOLUME-SURFACE VECTOR
C     MATALB   : INTEGER MATERIAL-FACE VECTOR
C----
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      NDIM     =ISTATE(1)
      NSUR    =-ISTATE(2)
      NVOL     =ISTATE(3)
      NTOTCL   =ISTATE(4)
      MAXR     =ISTATE(5)
      NUNK     =ISTATE(6)
      NEXTGE   =ISTATE(7)
      ALLOCATE(KEYMRG(NUNK),MATALB(NUNK),VOLSUR(NUNK))
      CALL LCMGET(IPTRK,'KEYMRG      ',KEYMRG)
      CALL LCMGET(IPTRK,'MATALB      ',MATALB)
      CALL LCMGET(IPTRK,'VOLSUR      ',VOLSUR)
      CALL LCMSIX(IPTRK,'EXCELL      ',2)
      CALL LCMGET(IPTRK,'ICODE       ',ICODE )
      CALL LCMGET(IPTRK,'NCODE       ',NCODE )
      CALL LCMGET(IPTRK,'ALBEDO      ',ALBEDO)
      NV=  NVOL
      NS= -NSUR
      ITOPT=LTRK-1
      NANGL=0
C----
C  EXCELL-TYPE VECTORIZATION - THE TRACKING FILE IS NOT COMPUTED
C----
      IF(INSB.EQ.2) RETURN
C----
C  Intrinsic symmetries used in geometry
C  Use these to simplify tracking unless 
C  NOSYMM tracking option activated
C----
      LCLSYM(1) =ISTATE(8)
      LCLSYM(2) =ISTATE(9)
      LCLSYM(3) =ISTATE(10)
      IF(ISYMM .NE. 0) THEN
        ISYMM=0
        IF(NDIM .EQ. 2) THEN
          IF(LCLSYM(1) .NE. 0) THEN
C----
C  X SYMMETRY
C----
            ISYMM=2
          ENDIF
          IF(LCLSYM(2) .NE. 0) THEN
            IF(ISYMM .EQ. 0) THEN
C----
C  Y SYMMETRY
C----
              ISYMM=4
            ELSE
C----
C  X AND Y SYMMETRY
C----
              ISYMM=8 
            ENDIF
          ENDIF
CC          IF(ISTATE(11) .NE. 0) THEN
CCC----
CCC  X-Y DIAGONAL SYMMETRY
CCC---- 
CC            IF(ISYMM .EQ. 0) THEN
CC              ISYMM=10
CC            ELSE
CC              ISYMM=12
CC            ENDIF
CC          ENDIF 
        ELSE
          IF(LCLSYM(1) .NE. 0) THEN
C----
C  X SYMMETRY
C----
            ISYMM=2
          ENDIF
          IF(LCLSYM(2) .NE. 0) THEN
            IF(ISYMM .EQ. 0) THEN
C----
C  Y SYMMETRY
C----
              ISYMM=4
            ELSE
C----
C  X AND Y SYMMETRY
C----
              ISYMM=8 
            ENDIF
          ENDIF 
          IF(LCLSYM(3) .NE. 0) THEN
C----
C  Z SYMMETRY
C----
            ISYMM=ISYMM+16
          ENDIF
        ENDIF
        IF(ISYMM .EQ. 0) ISYMM=1
      ENDIF
C----
C  READ THE GEOMETRY INFORMATION STORED ON IPTRK
C----      
      SUBMAX= 0
      LINMAX= 0
      CALL LCMSIX(IPTRK,'EXCELL      ',1)
      ALLOCATE(MINDIM(NTOTCL),MAXDIM(NTOTCL),ICORD(NTOTCL),
     > INDEX(4*NUNK))
      ALLOCATE(REMESH(MAXR))
      CALL LCMGET(IPTRK,'MINDIM      ',MINDIM)
      CALL LCMGET(IPTRK,'MAXDIM      ',MAXDIM)
      CALL LCMGET(IPTRK,'ICORD       ',ICORD )
      CALL LCMGET(IPTRK,'INDEX       ',INDEX )
      CALL LCMGET(IPTRK,'REMESH      ',REMESH)
      CALL LCMSIX(IPTRK,'EXCELL      ',2)
      IF (LPRISM) CALL LCMSIX(IPTRK,'PROJECTION  ',2)
C----
C  VERIFY SYMMETRY AND
C  STUDY TRACKING PARAMETERS. ARE THEY BASICALLY POSSIBLE ?
C----
      MXANGL=0
      IF(LTRK .EQ. 1)THEN
        NCOR= 1
        IF(NDIM .EQ. 2) THEN
          MXANGL=NANGLE
CC          IF(ISYMM .EQ. 12) THEN
CC            NANGL = NANGLE/4
          IF(ISYMM .GE. 2) THEN
            NANGL = (NANGLE+1)/2
          ELSE 
            NANGL = NANGLE
          ENDIF
          IF( RCUTOF.GT.0.0 ) NCOR= 2
        ELSE IF(NDIM .EQ. 3) THEN
          IF(MOD(NANGLE,2) .EQ. 1)THEN
            NANGLE=NANGLE+1
            WRITE(IOUT,'(/31H MESS = ONLY EVEN # EQN ANGLES )')
          ENDIF
          IF(NANGLE .GT. 16)THEN
            NANGLE=16
            WRITE(IOUT,'(/31H MESS = 16 IS MAX # EQN ANGLES )')
          ENDIF
          MXANGL=(NANGLE * (NANGLE+2)) / 2 
          IF(NEXTGE .EQ. 1) THEN
            NANGL = (NANGLE * (NANGLE+2)) / 8
          ELSE
            IF(ISYMM .EQ. 8 .OR. ISYMM .EQ. 24) THEN
              NANGL = (NANGLE * (NANGLE+2)) / 8 
            ELSE IF(ISYMM .EQ. 2  .OR. ISYMM .EQ. 4 .OR. 
     >              ISYMM .EQ. 18 .OR. ISYMM .EQ. 20 ) THEN
              NANGL = (NANGLE * (NANGLE+2)) / 4 
            ELSE
              NANGL = (NANGLE * (NANGLE+2)) / 2 
            ENDIF
          ENDIF
          IF(RCUTOF .GT. 0.0) NCOR= 4
        ENDIF
      ELSEIF( LTRK.EQ.2 )THEN
        NCOR  = 1
        IF(NANGLE .GT. 24) THEN
          NANGLE=30
        ELSE IF(NANGLE .GT. 20) THEN
          NANGLE=24
        ELSE IF(NANGLE .GT. 18) THEN
          NANGLE=20
        ELSE IF(NANGLE .GT. 14) THEN
          NANGLE=18
        ELSE IF(NANGLE .GT. 12) THEN
          NANGLE=14
        ELSE IF(NANGLE .GT. 8) THEN
          NANGLE=12
        ELSE 
          NANGLE=8
        ENDIF
        MXANGL=4*NANGLE
        ISTATE(11)=NANGLE
        CALL LCMPUT(IPTRK,'STATE-VECTOR',NSTATE,1,ISTATE)
        IF( NDIM.EQ.2 )THEN
           NANGL = 4*NANGLE
        ELSEIF( NDIM.EQ.3 )THEN
           CALL XABORT('XELTRK: *TSPC* NOT AVAILABLE FOR 3-D GEOMETRY')
        ENDIF
        CUTOFX= RCUTOF
      ENDIF
      IF(IPRT .GT. 1 .AND. NEXTGE .EQ. 0)THEN
C----
C  IF PRINT REQUIRED AND OVERALL CARTESIAN GEOMETRY
C  PRINT CARTESIAN REGION MAP 
C----
         NTX= MAXDIM(1)-MINDIM(1)
         NTY= MAXDIM(2)-MINDIM(2)
         NTZ= MAXDIM(3)-MINDIM(3)
         NTR=0
         DO 103 ICL=3,NTOTCL-1
            NTR= MAX(NTR,MAXDIM(ICL+1)-MINDIM(ICL+1)+1)
  103    CONTINUE
         CALL XELGPR(NDIM,NTX,NTY,NTZ,NTR,ISYMM,NSUR,NVOL,NTOTCL,
     >               MINDIM,MAXDIM,KEYMRG,INDEX,MATALB)
      ENDIF
C
C. 3) DO THE TRACKING OF THE EXACT GEOMETRY FOR *NEWT* OPTION.
      IF( IDISP.GT.0.AND.LTRK.NE.0 )THEN
         NC= NTOTCL - 3
         IF( IPRT.GE.1 )THEN
            WRITE(IOUT,'(1H )')
            IF    ( NC.EQ.0 )THEN
               WRITE(IOUT,'(/21H NOW, TRACKING   >>> ,A12,
     >                   13H GEOMETRY <<<,
     >                   13H    (WITH  NO,11H CYLINDER ) /)')
     >                                            GEONAM
            ELSEIF( NC.EQ.1 )THEN
               WRITE(IOUT,'(/21H NOW, TRACKING   >>> ,A12,
     >                   13H GEOMETRY <<<,
     >                   13H    (WITH ONE,11H CYLINDER ) /)')
     >                                            GEONAM
            ELSE
               WRITE(IOUT,'(/21H NOW, TRACKING   >>> ,A12,
     >                   13H GEOMETRY <<<,
     >                   10H    (WITH ,I3,11H CYLINDERS) /)')
     >                                            GEONAM, NC
            ENDIF
         ENDIF 
         ALLOCATE(ICUR(NTOTCL),INCR(NTOTCL))
         ALLOCATE(CONV(NTOTCL),TRKBEG(NTOTCL),TRKDIR(NTOTCL))
C
C  3.0)  WRITE FIRST RECORDS OF THE UNNORMALIZED TRACKING FILE
         IF( LTRK.EQ.1 )THEN
            SUBMAX= 1
            LINMAX= 2*NVOL + 10
         ELSE
C
C           REQUIRED CHANGE LINMAX FROM  2*NANGL*(2*NVOL + 8)
C                                    TO  2*NANGL*(2*NVOL + 16)
C           TO TAKE INTO ACCOUNT INITIAL AND FINAL SURFACE
C           FOR PERIODIC BC
            SUBMAX= NANGL
            LINMAX= 2*NANGL*(2*NVOL + 16)
         ENDIF
         CTRK  = '$TRK'
         NALBG = 6
         WRITE(IFTEMP) CTRK,5,0,0
         COMENT='CREATOR  : DRAGON'
         WRITE(IFTEMP) COMENT
         COMENT='MODULE   : XELTRK'
         WRITE(IFTEMP) COMENT
         COMENT='TYPE     : CARTESIAN'
         WRITE(IFTEMP) COMENT
         COMENT='GEOMETRY : '//GEONAM
         WRITE(IFTEMP) COMENT
         COMENT=TITREC
         WRITE(IFTEMP) COMENT
         WRITE(IFTEMP) NDIM,ITOPT,NV,NS,NALBG,NCOR,NANGL,SUBMAX,LINMAX
         WRITE(IFTEMP) (VOLSUR(JJ),JJ=1,NUNK)
         WRITE(IFTEMP) (MATALB(JJ),JJ=1,NUNK)
         WRITE(IFTEMP) ( ICODE(JJ),JJ=1,NALBG)
         WRITE(IFTEMP) (ALBEDO(JJ),JJ=1,NALBG)
         ALLOCATE(NUMERO(LINMAX))
         ALLOCATE(LENGHT(LINMAX),ANGLES(3*MXANGL),DENSTY(MXANGL),
     >            DDENWT(MXANGL))
         MXSUB= SUBMAX
         MXSEG= LINMAX
         IF( LTRK.EQ.1 )THEN
C
C  3.1)     THE REGULAR TRACKING
            IF( NDIM.EQ.3 )THEN
               CALL XELTI3(  IPRT,IFTEMP,NANGLE,DENUSR, ISYMM,ANGLES,
     >                     DENSTY,NTOTCL,NEXTGE,  MAXR,REMESH,LINMAX,
     >                     RCUTOF,  NSUR,  NVOL, INDEX,MINDIM,MAXDIM,
     >                      ICORD,  INCR,  ICUR,TRKBEG,  CONV,TRKDIR,
     >                   LENGHT, NUMERO ,DDENWT)
               CALL LCMPUT(IPTRK,'TrackingDirc',2*NDIM*NANGLE,4,ANGLES)
               CALL LCMPUT(IPTRK,'TrackingTrkW',2*NANGLE,4,DDENWT)
               CALL LCMPUT(IPTRK,'TrackingSpaD',2*NANGLE,4,DENSTY)
            ELSEIF( NDIM.EQ.2 )THEN
               ICUR(3)= MINDIM(3)
               CONV(3)= 1.0E+36
               CALL XELTI2(  IPRT,IFTEMP,NANGLE,DENUSR, ISYMM,ANGLES,
     >                     DENSTY,NTOTCL,  MAXR,REMESH,LINMAX,RCUTOF,
     >                       NSUR,  NVOL, INDEX,MINDIM,MAXDIM, ICORD,
     >                       INCR,ICUR,  TRKBEG,  CONV,TRKDIR,LENGHT,
     >                     NUMERO,DDENWT)
               CALL LCMPUT(IPTRK,'TrackingDirc',2*NDIM*NANGLE,4,ANGLES)
               CALL LCMPUT(IPTRK,'TrackingTrkW',2*NANGLE,4,DDENWT)
               CALL LCMPUT(IPTRK,'TrackingSpaD',2*NANGLE,4,DENSTY)
            ENDIF
         ELSEIF( LTRK.EQ.2 )THEN
C
C  3.2)     THE CYCLIC TRACKING.
            ICUR(3)= MINDIM(3)
            CONV(3)= 1.0E+36
            CALL XELTS2(  IPRT, IFTEMP, NANGLE, DENUSR, NCODE,ANGLES,
     >                  DENSTY, SWZERO, NTOTCL,   MAXR,REMESH,SUBMAX,
     >                  LINMAX,   NSUR,   NVOL, MATALB, INDEX,MINDIM,
     >                  MAXDIM,  ICORD,   INCR,   ICUR,TRKBEG,  CONV,
     >                  TRKDIR, LENGHT, NUMERO,DDENWT)
               CALL LCMPUT(IPTRK,'TrackingDirc',2*NDIM*NANGLE,4,ANGLES)
               CALL LCMPUT(IPTRK,'TrackingTrkW',2*NANGLE,4,DDENWT)
               CALL LCMPUT(IPTRK,'TrackingSpaD',2*NANGLE,4,DENSTY)
         ENDIF
         DEALLOCATE(DENSTY,ANGLES,LENGHT,NUMERO,TRKDIR,TRKBEG,CONV,INCR,
     >   ICUR,DDENWT)
      ENDIF
      DEALLOCATE(REMESH,INDEX,ICORD,MAXDIM,MINDIM)
C---
C  IF A PRISMATIC 3D TRACKING IS REQUESTED, 
C  ALLOCATE GEOMETRIC STRUCTURES (SEE COMMON CEXGEO) CORRESPONDING TO
C  THE 3D INITIAL GEOMETRY
C---
      IF (LPRISM) THEN
         CALL LCMSIX(IPTRK,'EXCELL      ',1)
         CALL XDISET(ISTATE,NSTATE,0)
         CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
         NDIM     =ISTATE(1)
         NSUR    =-ISTATE(2)
         NVOL     =ISTATE(3)
         CALL LCMSIX(IPTRK,'EXCELL      ',2)
         CALL LCMGET(IPTRK,'ICODE       ',ICODE )
         CALL LCMGET(IPTRK,'NCODE       ',NCODE )
         CALL LCMGET(IPTRK,'ALBEDO      ',ALBEDO)
         NV=  NVOL
         NS= -NSUR
      ENDIF
      DEALLOCATE(VOLSUR,MATALB,KEYMRG)
C
      RETURN
      END
