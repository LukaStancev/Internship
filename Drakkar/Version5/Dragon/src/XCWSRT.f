*DECK XCWSRT
      SUBROUTINE XCWSRT(IPRT,MXSEG,SEGLEN,NRSEG,NNSEG,NTSEG)
C
C-------------------------    XCWSRT    -------------------------------
C
C 1- SUBROUTINE STATISTICS:
C     NAME     : XCWSRT
C     USE      : SORT REGION INTERSECTION BY POSITION
C     MODIFIED : 94-01-19
C     AUTHOR   : G. MARLEAU
C
C 2- PARAMETERS:
C  INPUT
C     IPRT     : PRINT LEVEL                           I
C     MXSEG    : CURRENT MAXIMUM TRACK LENGTH          I
C  INPUT/OUTPUT
C     SEGLEN   : LENGTH OF TRACK                       D(MXSEG)
C     NRSEG    : REGION CROSSED BY TRACK               I(MXSEG)
C     NNSEG    : REGION CROSSED BY TRACK  (LEFT)       I(MXSEG)
C  OUTPUT
C     NTSEG    : TOTAL NUMBER OF SEGMENTS              I
C
C  3-INTERNAL PARAMETERS
C     IUNOUT   : OUT UNIT NUMBER = 6
C
C----------------------------------------------------------------------
C
      PARAMETER (IUNOUT=6)
      INTEGER    IPRT,MXSEG,NRSEG(*),NNSEG(*),NTSEG
      DOUBLE PRECISION SEGLEN(*)
C----
C LOCAL VARIABLES
C----
      INTEGER REFNR,REFNN
      DOUBLE PRECISION REFSL
C----
C  REMOVE TERM WITH NRSEG<=0
C----
      NTSEG=0
      DO 100 IS=1,MXSEG-1
        IF(NRSEG(IS).GT.0) THEN
          NTSEG=NTSEG+1
          NRSEG(NTSEG)=NRSEG(IS)
          NNSEG(NTSEG)=NNSEG(IS)
          SEGLEN(NTSEG)=SEGLEN(IS)
        ENDIF
 100  CONTINUE
      NSEG=NTSEG+1
      NRSEG(NSEG)=NRSEG(MXSEG)
      NNSEG(NSEG)=NNSEG(MXSEG)
      SEGLEN(NSEG)=SEGLEN(MXSEG)
      IF(IPRT.GE.200) THEN
        WRITE(IUNOUT,6000)
        WRITE(IUNOUT,6010) (IIJJ,SEGLEN(IIJJ),NNSEG(IIJJ),
     >                      NRSEG(IIJJ),IIJJ=1,NSEG)
      ENDIF
C----
C  SORT FROM MINIMUM TO MAXIMUM
C----
      DO 110 IS=2,NSEG
        REFSL=SEGLEN(IS)
        REFNR=NRSEG(IS)
        REFNN=NNSEG(IS)
        DO 111 JS=IS-1,1,-1
          KS=JS
          IF(SEGLEN(JS).GT.REFSL) THEN
            SEGLEN(JS+1)=SEGLEN(JS)
            NRSEG(JS+1)=NRSEG(JS)
            NNSEG(JS+1)=NNSEG(JS)
          ELSE
            GO TO 112
          ENDIF
 111    CONTINUE
        KS=0
 112    CONTINUE
        SEGLEN(KS+1)=REFSL
        NRSEG(KS+1)=REFNR
        NNSEG(KS+1)=REFNN
 110  CONTINUE
      IF(IPRT.GE.200) THEN
        WRITE(IUNOUT,6001)
        WRITE(IUNOUT,6010) (IIJJ,SEGLEN(IIJJ),NNSEG(IIJJ),
     >                      NRSEG(IIJJ),IIJJ=1,NSEG)
      ENDIF
C----
C  CHECK FOR ROD INTERSECTION WITH ANNULUS OR
C  ANNULUS LOCATED BETWEEN ROD SETS
C----
      DO 120 IS=1,NSEG
        NTB=NRSEG(IS)
        NFB=NNSEG(IS)
        IF(NTB.GT.0) THEN
          IF(NTB.LT.NFB) THEN
            DO 121 JS=IS+1,NSEG
              NTE=NRSEG(JS)
              NFE=NNSEG(JS)
              IF((NTE.EQ.NFB).AND.(NFE.EQ.NTB)) GO TO 122
              IF(NTE.GT.NTB) THEN
                NRSEG(JS)=NTB
              ENDIF
              IF(ABS(NFE).GT.NTB) THEN
                IF(NFE.LT.0) THEN
                  NNSEG(JS)=-NTB
                ELSE
                  NNSEG(JS)=NTB
                ENDIF
              ENDIF
 121        CONTINUE
          ENDIF
        ENDIF
 122    CONTINUE
        IF(NFB.GT.0) THEN
          DO 123 JS=IS-1,1,-1
            NTE=NRSEG(JS)
            IF(NTE.GT.0) THEN
              IF(NFB.NE.NTE) THEN
                NRSEG(IS)=0
              ENDIF
              GO TO 124
            ENDIF
 123      CONTINUE
        ENDIF
 124    CONTINUE
 120  CONTINUE
C----
C  REMOVE NEW TERMS WITH NRSEG<=0
C----
      NTSEG=0
      DO 130 IS=1,NSEG-1
        IF(NRSEG(IS).GT.0) THEN
          NTSEG=NTSEG+1
          NRSEG(NTSEG)=NRSEG(IS)
          NNSEG(NTSEG)=NNSEG(IS)
          SEGLEN(NTSEG)=SEGLEN(IS)
        ENDIF
 130  CONTINUE
      NSEG=NTSEG+1
      NRSEG(NSEG)=NRSEG(MXSEG)
      NNSEG(NSEG)=NNSEG(MXSEG)
      SEGLEN(NSEG)=SEGLEN(MXSEG)
      RETURN
C----
C  FORMATS
C----
 6000 FORMAT(' COMPRESSED TRACKING FILE'/
     >5X,'NUMBER',7X,'POSITION',4X,'BEFORE',5X,'AFTER')
 6001 FORMAT(' SORTED TRACKING FILE'/
     >5X,'NUMBER',7X,'POSITION',4X,'BEFORE',5X,'AFTER')
 6010 FORMAT((1X,I10,F15.7,2I10))
      END
