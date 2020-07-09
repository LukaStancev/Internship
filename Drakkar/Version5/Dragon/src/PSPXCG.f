*DECK PSPXCG
      SUBROUTINE PSPXCG(IPRINT,ISPSP ,ICOLR ,NBAN  ,NRT   ,MSROD ,
     >                  NSURX ,NSUR  ,NVOL  ,COTE  ,
     >                  RAN   ,NRODS ,RODS  ,RODR  ,NRINFO,NRODR ,
     >                  NXRI  ,KEYMRG,COLREG)
C
C-------------------------    PSPXCG    -------------------------------
C
C 1- SUBROUTINE STATISTICS:
C     NAME     : PSPXCG
C     USE      : GRAPHICS FOR 2-D CLUSTER GEOMETRY
C     MODIFIED : 99-03-15
C     AUTHOR   : G. MARLEAU
C
C 2- PARAMETERS:
C  INPUT
C     IPRINT : PRINT LEVEL                             I
C     ISPSP  : PSP FILE UNIT                           I
C     ICOLR  : COLOR SET USED                          I
C              = -4 FILL HSB WITH NO-CONTOUR
C              = -3 FILL CMYK WITH NO-CONTOUR
C              = -2 FILL RGB WITH NO-CONTOUR
C              = -1 FILL BW WITH NO-CONTOUR
C              =  0 NO FILL CONTOUR ONLY
C              =  1 FILL BW AND CONTOUR
C              =  2 FILL RGB AND CONTOUR
C              =  3 FILL CMYK AND CONTOUR
C              =  4 FILL HSB AND CONTOUR
C     NBAN    : NUMBER OF CONCENTRIC REGIONS           I
C     NRT     : NUMBER OF ROD TYPES                    I
C     MSROD   : MAXIMUM NUMBER OF SUBROD PER RODS      I
C     NSURX   : NUMBER OF SURFACES                     I
C     NSUR    : NUMBER OF SURFACES.                    I
C     NVOL    : NUMBER OF REGIONS                      I
C     COTE    : Y DIMENSION FOR RECTANGLE              R
C     RAN     : RADIUS/LATTICE SIDE OF REGION          R(NBAN)
C     NRODS   : INTEGER DESCRIPTION OF ROD  TYPE       I(3,NRT)
C               NRODS(1,IRT) = NUMBER OF ROD
C               NRODS(2,IRT) = NUMBER OF SUBRODS IN ROD
C               NRODS(3,IRT) = ASSOCIATED ANNULUS
C     RODS    : DESCRIPTION OF ROD OF A GIVEN TYPE     R(2,NRT)
C               RODS(1,IRT) = ROD CENTER RADIUS
C               RODS(2,IRT) = ANGLE POSITION OF ONE ROD
C     RODR    : SUBROD RADIUS                          R(MSROD,NRT)
C     NRINFO  : ANNULAR REGION CONTENT                 I(2,NBAN)
C               NRINFO(1,IAN) = NEW REGION NUMBER
C               NRINFO(2,IAN) = ASSOCIATED CLUSTER
C                             = 0 NO CLUSTER
C     NRODR    : SUBROD REGION                         I(NRT)
C     NXRI     : ANNULAR REGION CONTENT MULTI-ROD      I(NRT,NBAN)
C     KEYMRG   : MERGE INDEX                           I(NSUR:NVOL)
C     COLREG   : REGION COLOR                          I(4,NVOL)
C  3-INTERNAL PARAMETERS
C     IUNOUT   : OUT UNIT NUMBER = 6
C     PI       : NUMERICAL VALUE OF PI
C----------------------------------------------------------------------
C
      IMPLICIT         NONE
      INTEGER          IOUT,NPTS
      REAL             PI,DIMX,DIMY,WLINE
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NPTS=6,PI=3.1415926535897932,
     >                 DIMX=3.5,DIMY=3.5,WLINE=0.002,NAMSBR='PSPXCG')
C----
C  ROUTINE PARAMETERS
C----
      INTEGER          IPRINT,ISPSP,ICOLR,NBAN,NRT,MSROD,NSURX,
     >                 NSUR,NVOL
      INTEGER          NRODS(3,NRT),NRINFO(2,NBAN),NRODR(NRT),
     >                 NXRI(NRT,NBAN),KEYMRG(NSUR:NVOL)
      REAL             COTE,RAN(NBAN),RODS(2,NRT),
     >                 RODR(MSROD,NRT),COLREG(4,NVOL)
C----
C  LOCAL PARAMETERS
C----
      INTEGER          ICOL,ICONT,IVOL,IMRG,NTAN,IPT,IRT,
     >                 NPROD,NINRD,IROD,ISBR,IAN,NSEG,KRT,JRT
      REAL             XYPOS(2,NPTS),RADEQ,FACT,ANGD,ANGR(2),
     >                 DANGR,RPIN,RROD,XINT,ANGA,
     >                 WLFAC
      INTEGER          KFS,KFR,KSS,KSR
C----
C  INITIALIZE
C    ICOL  FOR COLOR   (NONE, BW, RGB)
C    ICONT FOR CONTOUR (WITH OR WITHOUT CONTOUR)
C----
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
      KFS=0
      KFR=0
      KSS=0
      KSR=0
      ICONT=1
      WLFAC=1.0
      ICOL=ABS(ICOLR)
      IF(ICOLR .EQ. 0) THEN
        WLFAC=2.5
      ELSE IF(ICOLR .LT. 0) THEN
        ICONT=0
      ELSE
        KFS=1
        KSR=1
      ENDIF
C----
C  LOCATE PEN AT CENTER OF CELL
C  DETERMINE DIMENSION OF GRAPH USING CELL LIMIT
C  FOR HEXAGONAL CELL PRINT HEXAGONAL REGION
C  FOR CARTESIAN CELL PRINT CARTESIAN REGION
C----
      XYPOS(1,1)=DIMX
      XYPOS(2,1)=DIMY
      CALL PSMOVE(ISPSP,XYPOS,-3)
      IF(NSURX.EQ.6) THEN
        RADEQ=RAN(NBAN)
        FACT=DIMX/RADEQ
        RADEQ=DIMX
        NTAN=NBAN-1
C----
C  POSITION OF POINTS DEFINING THE HEXAGONAL SHAPE TO FILL
C----
        ANGD=0.0
        DO 100 IPT=1,NSURX
          XYPOS(1,IPT)=COS(ANGD)*RADEQ
          XYPOS(2,IPT)=SIN(ANGD)*RADEQ
          ANGD=ANGD+PI/3.0
 100    CONTINUE
        IVOL=NRINFO(1,NBAN)
        IMRG=KEYMRG(IVOL)
C----
C  FILL IF REQUIRED
C----
        CALL PSDREG(ISPSP,NSURX,XYPOS)
        IF(ICOL. GT. 0) THEN
          CALL PSFILL(ISPSP,ICOL,COLREG(1,IVOL),KFS,KFR)
        ENDIF
C----
C  STROKE CONTOUR IF REQUIRED
C----
        IF(ICONT .EQ. 1) THEN
          CALL PSSTRK(ISPSP,WLFAC*WLINE,KSS,KSR)
        ENDIF
      ELSE IF(NSURX.EQ.4) THEN
        RADEQ=0.5*MAX(RAN(NBAN),COTE)
        FACT=DIMX/RADEQ
        NTAN=NBAN-1
        XYPOS(1,1)=FACT*RAN(NBAN)/2
        XYPOS(2,1)=FACT*COTE/2
        XYPOS(1,2)=-XYPOS(1,1)
        XYPOS(2,2)=XYPOS(2,1)
        XYPOS(1,3)=XYPOS(1,2)
        XYPOS(2,3)=-XYPOS(2,2)
        XYPOS(1,4)=XYPOS(1,1)
        XYPOS(2,4)=XYPOS(2,3)
        IVOL=NRINFO(1,NBAN)
        IMRG=KEYMRG(IVOL)
C----
C  FILL IF REQUIRED
C----
        CALL PSDREG(ISPSP,NSURX,XYPOS)
        IF(ICOL. GT. 0) THEN
          CALL PSFILL(ISPSP,ICOL,COLREG(1,IVOL),KFS,KFR)
        ENDIF
C----
C  STROKE CONTOUR IF REQUIRED
C----
        IF(ICONT .EQ. 1) THEN
          CALL PSSTRK(ISPSP,WLFAC*WLINE,KSS,KSR)
        ENDIF
      ELSE
        FACT=DIMX/RAN(NBAN)
        NTAN=NBAN
      ENDIF
C----
C  ANNULAR REGIONS
C----
      DO 110 IAN=NTAN,1,-1
        RADEQ=FACT*RAN(IAN)
        XYPOS(1,1)=0.0
        XYPOS(2,1)=0.0
        IVOL=NRINFO(1,IAN)
        IMRG=KEYMRG(IVOL)
C----
C  FILL IF REQUIRED
C----
        IF(ICOL. GT. 0) THEN
          CALL PSDCIR(ISPSP,XYPOS,RADEQ)
          CALL PSFILL(ISPSP,ICOL,COLREG(1,IVOL),0,0)
        ENDIF
C----
C  STROKE CONTOUR IF REQUIRED
C----
        IF(ICONT .EQ. 1) THEN
          IF(NRINFO(2,IAN) .NE. 0) THEN
            NSEG=0
            DO 111 KRT=NRINFO(2,IAN),1,-1
              JRT=NXRI(KRT,IAN)
              IF(JRT .GT. 1000000 .AND. JRT .LT. 3000000) THEN
                IRT=MOD(JRT,1000000)
                NSEG=NSEG+1
C----
C  IF ANNULAR REGION CUT BY PINS
C  DRAW ARC SEGMENT
C----
                NPROD=NRODS(1,IRT)
                NINRD=NRODS(2,IRT)
                DANGR=2.*PI/FLOAT(NPROD)
                ANGD=RODS(2,IRT)
                RROD=FACT*RODR(NINRD,IRT)
                RPIN=FACT*RODS(1,IRT)
C----
C  ANNULUS INTERSECT RODS
C  1) FIND X (XINT) AND Y (YINT) INTERSECTION
C     XINT=(RADEQ**2+RPIN**2-RROD**2)/(2*RPIN)
C     YINT=SQRT(RAN**2-XINT**2)
C  2) FIND OPENNING ANGLE FOR VOLUME LIMITED BY
C     ANNULUS (ANGA)
C     ANGA=ACOS(XINT/RADEQ)
C----
                XINT=(RADEQ**2+RPIN**2-RROD**2)
     >               /(2.0*RPIN)
                ANGA=ACOS(XINT/RADEQ)
                DO 112 IROD=1,NPROD
                  ANGR(1)=180.0*(ANGD+ANGA)/PI
                  ANGD=ANGD+DANGR
                  ANGR(2)=180.0*(ANGD-ANGA)/PI
                  CALL PSLINW(ISPSP,WLFAC*WLINE)
                  CALL PSSARC(ISPSP,XYPOS,RADEQ,ANGR)
 112            CONTINUE
              ENDIF
 111        CONTINUE
            IF(NSEG .EQ. 0) THEN
              CALL PSDCIR(ISPSP,XYPOS,RADEQ)
              CALL PSSTRK(ISPSP,WLFAC*WLINE,0,0)
            ENDIF
          ELSE
C----
C  IF ANNULAR REGION NOT CUT BY PINS
C  STROKE CIRCLES
C----
            CALL PSDCIR(ISPSP,XYPOS,RADEQ)
            CALL PSSTRK(ISPSP,WLFAC*WLINE,0,0)
          ENDIF
        ENDIF
 110  CONTINUE
C----
C  ROD CLUSTER
C----
      DO 120 IRT=NRT,1,-1
        NPROD=NRODS(1,IRT)
        NINRD=NRODS(2,IRT)
        DANGR=2.*PI/FLOAT(NPROD)
        ANGD=RODS(2,IRT)
        RPIN=FACT*RODS(1,IRT)
        DO 121 IROD=1,NPROD
          XYPOS(1,1)=RPIN*COS(ANGD)
          XYPOS(2,1)=RPIN*SIN(ANGD)
          ANGD=ANGD+DANGR
          DO 122 ISBR=NINRD,1,-1
            IVOL=NRODR(IRT)-NINRD+ISBR
            IMRG=KEYMRG(IVOL)
            RADEQ=FACT*RODR(ISBR,IRT)
C----
C  FILL IF REQUIRED
C----
            CALL PSDCIR(ISPSP,XYPOS,RADEQ)
            IF(ICOL. GT. 0) THEN
              CALL PSFILL(ISPSP,ICOL,COLREG(1,IVOL),KFS,KFR)
            ENDIF
C----
C STROKE IF REQUIRED
C----
            IF(ICONT .EQ. 1) THEN
              CALL PSSTRK(ISPSP,WLFAC*WLINE,KSS,KSR)
            ENDIF
 122      CONTINUE
 121    CONTINUE
 120  CONTINUE
      XYPOS(1,1)=-DIMX
      XYPOS(2,1)=-DIMY
      CALL PSMOVE(ISPSP,XYPOS,-3)
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
      END
