*DECK XCWREC
      SUBROUTINE XCWREC(ANGD,SIDE,TRKPOS,LINTER,ROTPOS,INDS,IMS)
C
C-------------------------    XCWREC    -------------------------------
C
C 1- SUBROUTINE STATISTICS:
C     NAME     : XCWREC
C     USE      : TRACK OUTER RECTANGLE FOR 2-D CLUSTER
C     MODIFIED : 92-02-18
C     AUTHOR   : G. MARLEAU
C
C 2- PARAMETERS:
C  INPUT
C     ANGD     : TRACK DIRECTOR COSINES (COS(A),SIN(A)  D(2)
C     SIDE     : SIDE OF RECTANGLE                      D(2)
C     IMS      : SURFACE MERGE                          I(6)
C
C  INPUT/OUTPUT
C     TRKPOS   : ONE TRACK POINT AT INPUT  (*,1)        D(2,2)
C                TRACK ORIGIN AT OUTPUT    (*,1)
C                TRACK ORIGIN AT INPUT     (*,2)
C                TRACK END AT OUTPUT       (*,2)
C  OUTPUT
C     LINTER   : INTERSECTION LOGICAL                   L
C     ROTPOS   : POSITION WRT ROTATED AXIS              D(2,2)
C     INDS     : SURFACE OF INTERSECTION                I(2)
C
C----------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION ANGD(2),SIDE(2),TRKPOS(2,2),ROTPOS(2,2)
      INTEGER    IMS(6),INDS(2)
      LOGICAL    LINTER
C----
C  EQUATIONS FOR SIDES
C     SIDE 1: XR= SIDE(1)/2   (-SIDE(2)/2<=YR<=SIDE(2)/2)
C     SIDE 2: YR= SIDE(2)/2   (-SIDE(1)/2<=XR<=SIDE(1)/2)
C     SIDE 3: XR=-SIDE(1)/2   (-SIDE(2)/2<=YR<=SIDE(2)/2)
C     SIDE 4: YR=-SIDE(2)/2   (-SIDE(1)/2<=XR<=SIDE(1)/2)
C  TRACK EQUATION
C             YR=  TAN(ANGD)*(XR-TRKPOS(1,1))+TRKPOS(2,1)
C          OR XR=COTAN(ANGD)*(YR-TRKPOS(2,1))+TRKPOS(1,1)
C----
      YTOP=0.5D0*SIDE(2)
      XTOP=0.5D0*SIDE(1)
      LINTER=.FALSE.
      IF(ANGD(1).EQ.0.0D0) THEN
C----
C  TRACK PARALLEL TO Y
C  TRACK INTERCEPT SURFACE 4 AND 2
C----
        IF(ABS(TRKPOS(1,1)).LT.XTOP) THEN
          TRKPOS(1,2)=TRKPOS(1,1)
           IF(ANGD(2).LT.0.0) THEN
            INDS(2)=IMS(4)
            INDS(1)=IMS(2)
            TRKPOS(2,2)=-YTOP
            TRKPOS(2,1)=YTOP
          ELSE
            INDS(2)=IMS(2)
            INDS(1)=IMS(4)
            TRKPOS(2,2)=YTOP
            TRKPOS(2,1)=-YTOP
          ENDIF
          LINTER=.TRUE.
        ENDIF
      ELSE IF(ANGD(2).EQ.0.0D0) THEN
C----
C  TRACK PARALLEL TO X
C  TRACK INTERCEPT SURFACE 3 AND 1
C----
        IF(ABS(TRKPOS(2,1)).LT.YTOP) THEN
          TRKPOS(2,2)=TRKPOS(2,1)
          IF(ANGD(1).LT.0.0D0) THEN
            INDS(2)=IMS(3)
            INDS(1)=IMS(1)
            TRKPOS(1,2)=-XTOP
            TRKPOS(1,1)=XTOP
          ELSE
            INDS(2)=IMS(1)
            INDS(1)=IMS(3)
            TRKPOS(1,2)=XTOP
            TRKPOS(1,1)=-XTOP
          ENDIF
          LINTER=.TRUE.
        ENDIF
      ELSE
        NSEG=1
        COSAI=1.0/ANGD(1)
        SINAI=1.0/ANGD(2)
C----
C    SLOPEY=TAN(ANGD)
C    SLOPEX=COTAN(ANGD)
C    RINTY=TRKPOS(2,1)-SLOPEY*TRKPOS(1,1)
C    RINTX=TRKPOS(1,1)-SOLPEX*TRKPOS(2,1)
C----
        SLOPEY=ANGD(2)*COSAI
        SLOPEX=ANGD(1)*SINAI
        RINTY=TRKPOS(2,1)-SLOPEY*TRKPOS(1,1)
        RINTX=TRKPOS(1,1)-SLOPEX*TRKPOS(2,1)
C----
C  SURFACE 3: YR=RINTY-SLOPEY*XTOP
C            (-YTOP <=YR<= YTOP)
C----
        TRKPOS(2,NSEG)=RINTY-SLOPEY*XTOP
        IF( ABS(TRKPOS(2,NSEG)).LE.YTOP ) THEN
C----
C  TRACK INTERSEPT SURFACE 3
C----
          INDS(NSEG)=IMS(3)
          TRKPOS(1,NSEG)=-XTOP
          NSEG=NSEG+1
        ENDIF
C----
C  SURFACE 1: YR=RINTY+SLOPEY*XTOP
C            (-YTOP <=YR<= YTOP)
C----
        TRKPOS(2,NSEG)=RINTY+SLOPEY*XTOP
        IF( ABS(TRKPOS(2,NSEG)).LE.YTOP ) THEN
C----
C  TRACK INTERSEPT SURFACE 1
C----
          INDS(NSEG)=IMS(1)
          TRKPOS(1,NSEG)=XTOP
          IF(NSEG.EQ.2) GO TO 100
          NSEG=NSEG+1
        ENDIF
C----
C  SURFACE 4: XR=RINTX-SLOPEX*YTOP
C            (-XTOP <=XR<= XTOP)
C----
        TRKPOS(1,NSEG)=RINTX-SLOPEX*YTOP
        IF( ABS(TRKPOS(1,NSEG)).LE.XTOP ) THEN
C----
C  TRACK INTERSEPT SURFACE 4
C----
          INDS(NSEG)=IMS(4)
          TRKPOS(2,NSEG)=-YTOP
          IF(NSEG.EQ.2) GO TO 100
          NSEG=NSEG+1
        ENDIF
C----
C  SURFACE 2: XR=RINTX+SLOPEX*YTOP
C            (-XTOP <=XR<= XTOP)
C----
        TRKPOS(1,NSEG)=RINTX+SLOPEX*YTOP
        IF( ABS(TRKPOS(1,NSEG)).LE.XTOP ) THEN
C----
C  TRACK INTERSEPT SURFACE 2
C----
          INDS(NSEG)=IMS(2)
          TRKPOS(2,NSEG)=YTOP
          IF(NSEG.EQ.2) GO TO 100
          NSEG=NSEG+1
        ENDIF
 100    CONTINUE
        IF(NSEG.EQ.2) THEN
          LINTER=.TRUE.
C----
C  REORDER INTERSECTION POINTS FOR DIRECTION
C----
          IF(ANGD(1).LT.0.0D0) THEN
            IF(TRKPOS(1,1).GT.TRKPOS(1,2)) THEN
              TRKTMP=TRKPOS(1,2)
              TRKPOS(1,2)=TRKPOS(1,1)
              TRKPOS(1,1)=TRKTMP
              TRKTMP=TRKPOS(2,2)
              TRKPOS(2,2)=TRKPOS(2,1)
              TRKPOS(2,1)=TRKTMP
              INDT=INDS(2)
              INDS(2)=INDS(1)
              INDS(1)=INDT
            ENDIF
          ELSE
            IF(TRKPOS(1,1).GT.TRKPOS(1,2)) THEN
              TRKTMP=TRKPOS(1,2)
              TRKPOS(1,2)=TRKPOS(1,1)
              TRKPOS(1,1)=TRKTMP
              TRKTMP=TRKPOS(2,2)
              TRKPOS(2,2)=TRKPOS(2,1)
              TRKPOS(2,1)=TRKTMP
              INDT=INDS(2)
              INDS(2)=INDS(1)
              INDS(1)=INDT
            ENDIF
          ENDIF
          IF(ANGD(2).LT.0.0D0) THEN
            IF(TRKPOS(2,2).GT.TRKPOS(2,1)) THEN
              TRKTMP=TRKPOS(1,2)
              TRKPOS(1,2)=TRKPOS(1,1)
              TRKPOS(1,1)=TRKTMP
              TRKTMP=TRKPOS(2,2)
              TRKPOS(2,2)=TRKPOS(2,1)
              TRKPOS(2,1)=TRKTMP
              INDT=INDS(2)
              INDS(2)=INDS(1)
              INDS(1)=INDT
            ENDIF
          ELSE
            IF(TRKPOS(2,1).GT.TRKPOS(2,2)) THEN
              TRKTMP=TRKPOS(1,2)
              TRKPOS(1,2)=TRKPOS(1,1)
              TRKPOS(1,1)=TRKTMP
              TRKTMP=TRKPOS(2,2)
              TRKPOS(2,2)=TRKPOS(2,1)
              TRKPOS(2,1)=TRKTMP
              INDT=INDS(2)
              INDS(2)=INDS(1)
              INDS(1)=INDT
            ENDIF
          ENDIF
        ENDIF
      ENDIF
C----
C  ROTATE RECTANGLE BY ANGD
C----
      IF(LINTER) THEN
        DO 110 II=1,2
          ROTPOS(1,II)=ANGD(1)*TRKPOS(1,II)+ANGD(2)*TRKPOS(2,II)
          ROTPOS(2,II)=-ANGD(2)*TRKPOS(1,II)+ANGD(1)*TRKPOS(2,II)
 110    CONTINUE
      ENDIF
      RETURN
      END
