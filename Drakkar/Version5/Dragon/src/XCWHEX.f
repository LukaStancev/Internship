*DECK XCWHEX
      SUBROUTINE XCWHEX(ANGD,RADC,SIDE,LINTER,XPOS,INDS,IMS)
C
C-------------------------    XCWHEX    -------------------------------
C
C 1- SUBROUTINE STATISTICS:
C     NAME     : XCWHEX
C     USE      : TRACK OUTER HEXAGONE FOR 2-D CLUSTER GEOMETRY
C     MODIFIED : 91-07-26
C     AUTHOR   : G. MARLEAU
C
C 2- PARAMETERS:
C  INPUT
C     ANGD     : TRACK ANGLE                            R
C     RADC     : Y-POSITION OF TRACK > 0                R
C     SIDE     : SIDE OF HEXAGONE                       R
C     IMS      : SURFACE MERGE                          I(6)
C  OUTPUT
C     LINTER   : INTERSECTION LOGICAL                   L
C     XPOS     : POINTS OF INTERSECTION                 R(2)
C     INDS     : SURFACE OF INTERSECTION                I(2)
C
C  3-INTERNAL PARAMETERS
C     SQ3      : SQRT(3)
C     OSQ3     : 1.0/SQRT(3)
C----------------------------------------------------------------------
C
      PARAMETER (SQ3=1.73205080756887729,OSQ3=0.577350269189625795)
      INTEGER    IMS(6),INDS(2)
      LOGICAL    LINTER
      REAL       ANGD,RADC,SIDE,XPOS(2)
C----
C  EQUATIONS FOR SIDES OF HEXAGONE YR(XR)
C     SIDE 1: YR=-SQ3*XR+SQ3*SIDE   (          0 <=YR<= SQ3*SIDE/2)
C                                   (     SIDE/2 <=XR<= SIDE      )
C     SIDE 2: YR= SQ3*SIDE/2        (    -SIDE/2 <=XR<=  SIDE/2   )
C     SIDE 3: YR= SQ3*XR+SQ3*SIDE   (          0 <=YR<= SQ3*SIDE/2)
C                                   (      -SIDE <=XR<= -SIDE/2   )
C     SIDE 4: YR=-SQ3*XR-SQ3*SIDE   (-SQ3*SIDE/2 <=YR<= 0         )
C                                   (      -SIDE <=XR<= -SIDE/2   )
C     SIDE 5: YR=-SQ3*SIDE/2        (    -SIDE/2 <=XR<= SIDE/2    )
C     SIDE 6: YR= SQ3*XR-SQ3*SIDE   (-SQ3*SIDE/2 <=YR<= 0         )
C                                  (     SIDE/2 <=XR<= SIDE      )
C  EQUATIONS FOR SIDES OF HEXAGONE XR(YR)
C     SIDE 1: XR=-OSQ3*YR+SIDE      (          0 <=YR<= SQ3*SIDE/2)
C                                   (     SIDE/2 <=XR<= SIDE      )
C     SIDE 3: XR= OSQ3*YR-SIDE      (          0 <=YR<= SQ3*SIDE/2)
C                                   (      -SIDE <=XR<= -SIDE/2   )
C     SIDE 4: XR=-OSQ3*YR-SIDE      (-SQ3*SIDE/2 <=YR<= 0         )
C                                   (      -SIDE <=XR<= -SIDE/2   )
C     SIDE 6: XR= OSQ3*YR+SIDE      (-SQ3*SIDE/2 <=YR<= 0         )
C                                   (     SIDE/2 <=XR<= SIDE      )
C  TRACK EQUATION:
C             YR= SQ3*(SLOPEY*XR+RINTY)
C          OR XR= OSQ3*SLOPEX*YR-RINTX
C----
      YRINT=SQ3*SIDE
      YLIM=0.5*YRINT
      XLIM=0.5*SIDE
      SINA=SIN(ANGD)
      COSA=COS(ANGD)
      LINTER=.FALSE.
      IF(COSA.EQ.0.0) THEN
C----
C  TRACK PARALLEL TO Y
C----
        IF( RADC.LT. XLIM ) THEN
C----
C  TRACK INTERCEPT SURFACE 5 AND 2
C----
          IF(SINA.LT.0.0) THEN
            INDS(2)=IMS(5)
            INDS(1)=IMS(2)
          ELSE
            INDS(2)=IMS(2)
            INDS(1)=IMS(5)
          ENDIF
          XPOS(2)=YLIM
          XPOS(1)=-XPOS(2)
          LINTER=.TRUE.
        ELSE IF(RADC.LE.SIDE) THEN
C----
C  TRACK INTERCEPT SURFACE 3 AND 4 OR 6 AND 1
C----
          IF(SINA.LT.0.0) THEN
            INDS(2)=IMS(3)
            INDS(1)=IMS(4)
          ELSE
            INDS(2)=IMS(6)
            INDS(1)=IMS(1)
          ENDIF
          XPOS(2)=YRINT-SQ3*RADC
          XPOS(1)=-XPOS(2)
          LINTER=.TRUE.
        ENDIF
      ELSE IF(SINA.EQ.0.0) THEN
C----
C  TRACK PARALLEL TO X
C----
        IF(RADC.LE.YLIM ) THEN
C----
C  TRACK INTERCEPT SURFACE 6 AND 4
C----
          INDS(2)=IMS(4)
          INDS(1)=IMS(6)
          XPOS(2)= OSQ3*RADC+SIDE
          XPOS(1)=-XPOS(2)
          LINTER=.TRUE.
        ENDIF
      ELSE
        NSEG=0
        COSAI=1.0/COSA
        SINAI=1.0/SINA
        SLOPEY=OSQ3*SINA*COSAI
        SLOPEX=SQ3*COSA*SINAI
        RINTY=OSQ3*RADC*COSAI
        RINTX=RADC*SINAI
        XREF=RADC*COSAI*SINA
        OPSY=1.0/(1+SLOPEY)
        OMSY=1.0/(1-SLOPEY)
        XLSX=SLOPEX*XLIM
        SPRY=SIDE+RINTY
        SMRY=SIDE-RINTY
C----
C  SURFACE 1: XR=(SIDE-RINTY)/(1+SLOPEY)
C            (SIDE/2 <=XR<= SIDE)
C----
        XR=SMRY*OPSY
        IF( (XLIM.LE.XR) .AND. (XR.LE.SIDE) ) THEN
C----
C  TRACK INTERSEPT SURFACE 1
C----
          NSEG=NSEG+1
          INDS(NSEG)=IMS(1)
          XPOS(NSEG)=XR
        ENDIF
C----
C  SURFACE 2: XR= SLOPEX*SIDE/2-RINTX
C            (-SIDE/2 <=XR<= SIDE/2)
C----
        XR=XLSX-RINTX
        IF( ABS(XR).LE.XLIM ) THEN
C----
C  TRACK INTERSEPT SURFACE 2
C----
          NSEG=NSEG+1
          INDS(NSEG)=IMS(2)
          XPOS(NSEG)=XR
          IF(NSEG.EQ.2) GO TO 100
        ENDIF
C----
C  SURFACE 3: XR=-(SIDE-RINTY)/(1-SLOPEY)
C            (-SIDE <=XR<= -SIDE/2)
C----
        XR=-SMRY*OMSY
        IF( (-SIDE.LE.XR) .AND. (XR.LE.-XLIM) )THEN
C----
C  TRACK INTERSEPT SURFACE 3
C----
          NSEG=NSEG+1
          INDS(NSEG)=IMS(3)
          XPOS(NSEG)=XR
          IF(NSEG.EQ.2) GO TO 100
        ENDIF
C----
C  SURFACE 4: XR=-(SIDE+RINTY)/(1+SLOPEY)
C            (-SIDE <=XR<= -SIDE/2)
C----
        XR=-SPRY*OPSY
        IF( (-SIDE.LE.XR) .AND. (XR.LE.-XLIM) ) THEN
C----
C  TRACK INTERSEPT SURFACE 4
C----
          NSEG=NSEG+1
          INDS(NSEG)=IMS(4)
          XPOS(NSEG)=XR
          IF(NSEG.EQ.2) GO TO 100
        ENDIF
C----
C  SURFACE 5: XR=-SLOPEX*SIDE/2-RINTX
C            (-SIDE/2 <=XR<= SIDE/2)
C----
        XR=-XLSX-RINTX
        IF( ABS(XR).LE.XLIM ) THEN
C----
C  TRACK INTERSEPT SURFACE 5
C----
          NSEG=NSEG+1
          INDS(NSEG)=IMS(5)
          XPOS(NSEG)=XR
          IF(NSEG.EQ.2) GO TO 100
        ENDIF
C----
C  SURFACE 6: XR=(RINTY+SIDE)/(1-SLOPEY)
C            (SIDE/2 <=XR<= SIDE)
C----
        XR=SPRY*OMSY
        IF( (XLIM.LE.XR) .AND. (XR.LE.SIDE) ) THEN
C----
C  TRACK INTERSEPT SURFACE 6
C----
          NSEG=NSEG+1
          INDS(NSEG)=IMS(6)
          XPOS(NSEG)=XR
        ENDIF
 100    CONTINUE
        IF(NSEG.EQ.2) THEN
          LINTER=.TRUE.
C----
C  ROTATE HEXAGONE BY -ANGD
C----
          XPOS(1)=XREF+XPOS(1)*COSAI
          XPOS(2)=XREF+XPOS(2)*COSAI
          IF( XPOS(1).GT.XPOS(2) ) THEN
            INDT=INDS(2)
            INDS(2)=INDS(1)
            INDS(1)=INDT
            XPOST=XPOS(2)
            XPOS(2)=XPOS(1)
            XPOS(1)=XPOST
          ENDIF
        ENDIF
      ENDIF
      RETURN
      END
