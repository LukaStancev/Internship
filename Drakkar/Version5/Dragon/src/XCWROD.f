*DECK XCWROD
      SUBROUTINE XCWROD(NRIN,NRODS,NRODR,RODR,RODP,RADC,NFSEG,NLSEG,
     >                  SEGLEN,NRSEG,NNSEG)
C
C-------------------------    XCWROD    -------------------------------
C
C 1- SUBROUTINE STATISTICS:
C     NAME     : XCWROD
C     USE      : PERFORM ROD TRACKING FOR 2-D CLUSTER GEOMETRY
C     MODIFIED : 92-02-18
C     AUTHOR   : G. MARLEAU
C
C 2- PARAMETERS:
C  INPUT
C     NRIN     : CURRENT REGION NUMBER                 I
C     NRODS    : INTEGER DESCRIPTION OF ROD  TYPE      I(3)
C                NRODS(1) = NUMBER OF ROD
C                NRODS(2) = NUMBER OF SUBRODS IN ROD
C     NRODR    : SUBROD REGION                         I
C     RODR     : SUBROD RADIUS                         R(MSROD)
C     RODP     : ROD POSITION                          R(2,NRODS)
C                RODP(1,IRD) = X-POSITION
C                RODP(2,IRD) = Y-POSITION
C     RADC     : Y-POSITION OF TRACK                   D
C  OUTPUT
C     NFSEG    : INITIAL SEGMENT POSITION              I
C     NLSEG    : FINAL SEGMENT POSITION                I
C     SEGLEN   : LENGTH OF TRACK                       D(MRTRK)
C     NRSEG    : REGION CROSSED BY TRACK               I(MRTRK)
C     NNSEG    : REGION CROSSED BY TRACK  (LEFT)       I(MXSEG)
C
C----------------------------------------------------------------------
C
      INTEGER    NRIN,NRODS(2),NRODR,NFSEG,NLSEG,NRSEG(*),NNSEG(*)
      REAL       RODR(*),RODP(2,*)
      DOUBLE PRECISION SEGLEN(*),RADC,RADR,RADR2
C----
C  FILL IN SEGLEN FROM THE END STARTING WITH ROD FURTHER FROM
C  TRACK STARTING POINT UNTIL CENTER OF TRACK REACHED
C----
      NPROD=(NRODS(1)+3)/2
      NSBR=NRODS(2)
      IF(RADC.GE.0.0D0) THEN
        IPDEB=1
        IPFIN=NPROD
        IPSTP=1
        IMDEB=NPROD
        IMFIN=1
        IMSTP=-1
      ELSE
        RADR=RODP(2,1)-RADC
        IF(ABS(RADR).LT.RODR(NSBR)) THEN
          IPDEB=NRODS(1)+1
          IPFIN=MAX(2,NRODS(1)+1-NPROD)
        ELSE
          IPDEB=NRODS(1)
          IPFIN=MAX(1,NRODS(1)-NPROD)
        ENDIF
        IPSTP=-1
        IMDEB=IPFIN
        IMFIN=IPDEB
        IMSTP=1
      ENDIF
      NXSEG=NLSEG
      DO 100 IRZ=IPDEB,IPFIN,IPSTP
        IF(IRZ.EQ.NRODS(1)+1) THEN
          IRD=1
        ELSE
          IRD=IRZ
        ENDIF
        RADR=RODP(2,IRD)-RADC
        RADR2=RADR*RADR
        NREG=NRIN
        IF( ABS(RADR).LT.RODR(NSBR) ) THEN
C----
C  ROD INTERCEPS
C----
          XTRA=SQRT(RODR(NSBR)*RODR(NSBR)-REAL(RADR2))
          XLST=RODP(1,IRD)+XTRA
          XFST=RODP(1,IRD)-XTRA
          IF(XLST.LT.0.0) THEN
C----
C  CENTER OF TRACK REACHED/EXIT
C----
            GO TO 1000
          ELSE
C----
C  SET POINTERS TO SEGLEN VECTOR W.R.T. LAST POSITION FREE
C----
            NFLSEG=NXSEG-2*NSBR
            NLLSEG=NXSEG
            NXSEG=NFLSEG
          ENDIF
          SEGLEN(NLLSEG)=XLST
          NRSEG(NLLSEG)=NREG
          NNSEG(NFLSEG+1)=-NREG
          NLLSEG=NLLSEG-1
          NREG=NRODR
          NFLSEG=NFLSEG+1
          SEGLEN(NFLSEG)=XFST
          NRSEG(NFLSEG)=NREG
          NNSEG(NLLSEG+1)=-NREG
          DO 110 ISBR=NSBR-1,1,-1
            IF( ABS(RADR).LT.RODR(ISBR) ) THEN
C----
C  SUBROD INTERCEPS
C----
              XTRA=SQRT(RODR(ISBR)*RODR(ISBR)-REAL(RADR2))
              SEGLEN(NLLSEG)=RODP(1,IRD)+XTRA
              NRSEG(NLLSEG)=NREG
              NNSEG(NFLSEG+1)=-NREG
              NLLSEG=NLLSEG-1
              NREG=NREG-1
              NFLSEG=NFLSEG+1
              SEGLEN(NFLSEG)=RODP(1,IRD)-XTRA
              NRSEG(NFLSEG)=NREG
              NNSEG(NLLSEG+1)=-NREG
            ENDIF
 110      CONTINUE
        ENDIF
 100  CONTINUE
 1000 CONTINUE
      NLSEG=NXSEG
C----
C  FILL IN SEGLEN FROM THE BEGINNING STARTING WITH ROD CLOSEST FROM
C  TRACK STARTING POINT UNTIL CENTER OF TRACK REACHED
C----
      NXSEG=NFSEG
      DO 200 IRZ=IMDEB,IMFIN,IMSTP
        IF(IRZ.EQ.NRODS(1)+1) THEN
          IRD=1
        ELSE
          IRD=IRZ
        ENDIF
        RADR=RODP(2,IRD)-RADC
        RADR2=RADR*RADR
        NREG=NRIN
        IF( ABS(RADR).LT.RODR(NSBR) ) THEN
C----
C  ROD INTERCEPS
C----
          XTRA=SQRT(RODR(NSBR)*RODR(NSBR)-REAL(RADR2))
          XLST=RODP(1,IRD)+XTRA
          XFST=RODP(1,IRD)-XTRA
          IF(XLST.LT.0.0) THEN
C----
C  SET POINTERS TO SEGLEN VECTOR W.R.T. FIRST POSITION FREE
C----
            NLLSEG=NXSEG+2*NSBR
            NFLSEG=NXSEG
            NXSEG=NLLSEG
          ELSE
C----
C  CENTER OF TRACK REACHED/EXIT
C----
            GO TO 2000
          ENDIF
          SEGLEN(NLLSEG)=XLST
          NRSEG(NLLSEG)=NREG
          NNSEG(NFLSEG+1)=-NREG
          NLLSEG=NLLSEG-1
          NREG=NRODR
          NFLSEG=NFLSEG+1
          SEGLEN(NFLSEG)=XFST
          NRSEG(NFLSEG)=NREG
          NNSEG(NLLSEG+1)=-NREG
          DO 210 ISBR=NSBR-1,1,-1
            IF( ABS(RADR).LT.RODR(ISBR) ) THEN
C----
C  SUBROD INTERCEPS
C----
              XTRA=SQRT(RODR(ISBR)*RODR(ISBR)-REAL(RADR2))
              SEGLEN(NLLSEG)=RODP(1,IRD)+XTRA
              NRSEG(NLLSEG)=NREG
              NNSEG(NFLSEG+1)=-NREG
              NLLSEG=NLLSEG-1
              NREG=NREG-1
              NFLSEG=NFLSEG+1
              SEGLEN(NFLSEG)=RODP(1,IRD)-XTRA
              NRSEG(NFLSEG)=NREG
              NNSEG(NLLSEG+1)=-NREG
            ENDIF
 210      CONTINUE
        ENDIF
 200  CONTINUE
 2000 CONTINUE
      NFSEG=NXSEG
      RETURN
      END
