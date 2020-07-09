*DECK PSPXEL
      SUBROUTINE PSPXEL(IPRINT,ISPSP ,ICOLR ,NDIM  ,NSUR  ,NVOL  ,
     >                  NTOTCL,MAXR  ,
     >                  MINDIM,MAXDIM,KEYMRG,INDEX ,REMESH,COLREG)
C
C-------------------------    PSPXEL    -------------------------------
C
C 1- SUBROUTINE STATISTICS:
C     NAME     : PSPXEL
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
C     NDIM   : NUMBER OF DIMENSIONS                    I
C     NSUR   : NUMBER OF SURFACES.                     I
C     NVOL   : NUMBER OF ZONES.                        I
C     NTOTCL : NUMBER OF CYLINDERS                     I
C     MAXR   : DIMENSION OF REMESH VECTOR              I
C     MINDIM : MIN INDEX VALUES FOR AXES               I(NTOTCL)
C     MAXDIM : MAX INDEX VALUES FOR AXES               I(NTOTCL)
C     KEYMRG : MERGE INDEX                             I(NSUR:NVOL)
C     INDEX  : #ING OF SURFACES & ZONES.               I(4,NSUR:NVOL)
C     REMESH : MESHING                                 R(MAXR)
C     COLREG : REGION COLOR                            I(4,NVOL)
C----------------------------------------------------------------------
C
      IMPLICIT         NONE
      INTEGER          IOUT,NPTS,MXDIM,NXY,NINT
      CHARACTER        NAMSBR*6
      REAL             PI,DIMX,DIMY,WLINE
      PARAMETER       (IOUT=6,NPTS=4,MXDIM=3,NXY=2,NINT=16,
     >                 PI=3.1415926535897932,
     >                 DIMX=3.5,DIMY=3.5,WLINE=0.002,NAMSBR='PSPXEL')
C----
C  ROUTINE PARAMETERS
C----
      INTEGER          IPRINT,ISPSP,ICOLR,NDIM,
     >                 NSUR,NVOL,NTOTCL,MAXR
      INTEGER          MINDIM(NTOTCL),MAXDIM(NTOTCL),
     >                 KEYMRG(NSUR:NVOL),INDEX(4,NSUR:NVOL)
      REAL             REMESH(MAXR),COLREG(4,NVOL)
C----
C  LOCAL PARAMETERS
C----
      INTEGER          ICOL,ICONT,IDIR,IVOL,IMRG,
     >                 IX,IY,IR,ICL,NSEG,IORDER(NINT)
      REAL             WLFAC,RCIRC,OFFDIR(MXDIM),XYPOS(NXY,NPTS),
     >                 FACT,CENTER(NXY),RADANG(NXY,NINT)
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
C  COMPUTE THE CIRCUMSCRIBED RADIUS
C  THE COORDINATE FOR THE TRUE CENTER OF THE CELL
C----
      RCIRC= 0.0
      DO 100 IDIR=1,NDIM
        OFFDIR(IDIR)=0.5
     >              *(REMESH(MAXDIM(IDIR))+REMESH(MINDIM(IDIR)))
        RCIRC=MAX(RCIRC,
     >            0.5*(REMESH(MAXDIM(IDIR))-REMESH(MINDIM(IDIR))))
 100  CONTINUE
C----
C  LOCATE PEN AT CENTER OF CELL
C  DETERMINE DIMENSION OF GRAPH USING CELL LIMIT
C  FOR HEXAGONAL CELL PRINT HEXAGONAL REGION
C  FOR CARTESIAN CELL PRINT CARTESIAN REGION
C----
      XYPOS(1,1)=DIMX
      XYPOS(2,1)=DIMY
      CALL PSMOVE(ISPSP,XYPOS,-3)
      FACT=DIMX/RCIRC
C----
C  SCAN ALL REGIONS AND LOCATE POSITION
C  REGION NUMBER FROM INSIDE ANNULUS
C  TO EXTERIOR CARTESIAN
C----
      DO 110 IVOL=NVOL,1,-1
        IMRG=KEYMRG(IVOL)
        IF(IMRG .NE. 0) THEN
C----
C  CARTESIAN CELL POSITION IN X AND Y
C----
          IX=INDEX(1,IVOL)
          IY=INDEX(2,IVOL)
          XYPOS(1,1)=FACT*(REMESH(IX)-OFFDIR(1))
          XYPOS(2,1)=FACT*(REMESH(IY)-OFFDIR(2))
          XYPOS(1,2)=FACT*(REMESH(IX+1)-OFFDIR(1))
          XYPOS(2,2)=XYPOS(2,1)
          XYPOS(1,3)=XYPOS(1,2)
          XYPOS(2,3)=FACT*(REMESH(IY+1)-OFFDIR(2))
          XYPOS(1,4)=XYPOS(1,1)
          XYPOS(2,4)=XYPOS(2,3)
          IF(INDEX(4,IVOL) .EQ. 0) THEN
C----
C  CARTESIAN POSITION GEOMETRY LOCATED
C  COLOR AND TRACE IT
C----
            CALL PSDREG(ISPSP,NPTS,XYPOS)
            IF(ICOL. GT. 0) THEN
              CALL PSFILL(ISPSP,ICOL,COLREG(1,IVOL),KFS,KFR)
            ENDIF
            IF(ICONT.EQ.1) THEN
              CALL PSSTRK(ISPSP,WLFAC*WLINE,KSS,KSR)
            ENDIF
          ELSE
C----
C  CARTESIAN GEOMETRY CONTAINS ANNULAR SUBDIVISION
C  DETERMINE WHICH ANNULUS
C----
            DO 111 ICL=4,NTOTCL
              IR=INDEX(4,IVOL)
              IF( IR .GE. MINDIM(ICL)-1 .AND.
     >            IR .LT. MAXDIM(ICL)      ) THEN
C----
C  ANNULUS IS DETERMINED
C  LOCATE ANNULAR/CARTESIAN AND ORDER CARTESIAN POINTS
C  FOR GEOMETRY TRACING
C----
                CENTER(1)=FACT*(REMESH(MINDIM(ICL)-2)-OFFDIR(1))
                CENTER(2)=FACT*(REMESH(MINDIM(ICL)-1)-OFFDIR(2))
                RCIRC=FACT*SQRT(REMESH(IR+1))
                CALL PSPRAI(NINT,NPTS,XYPOS,CENTER,RCIRC,
     >                      NSEG,IORDER,RADANG)
C----
C  COLOR AND TRACE RESULT
C----
                CALL PSMOVE(ISPSP,CENTER,-3)
                CALL PSDRAI(ISPSP,NSEG,IORDER,CENTER,RADANG)
                IF(ICOL. GT. 0) THEN
                  CALL PSFILL(ISPSP,ICOL,COLREG(1,IVOL),KFS,KFR)
                ENDIF
                IF(ICONT.EQ.1) THEN
                  CALL PSSTRK(ISPSP,WLFAC*WLINE,KSS,KSR)
                ENDIF
                CENTER(1)=-CENTER(1)
                CENTER(2)=-CENTER(2)
                CALL PSMOVE(ISPSP,CENTER,-3)
                GO TO 115
              ENDIF
 111        CONTINUE
 115        CONTINUE
          ENDIF
        ENDIF
 110  CONTINUE
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
