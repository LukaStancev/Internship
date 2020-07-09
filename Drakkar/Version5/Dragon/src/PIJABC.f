*DECK PIJABC
      SUBROUTINE PIJABC(NREG,NSOUT,NPRB,SIGTAL,MATRT,PROB,PSST,PSVT)
C
C-------------------------    PIJABC    -------------------------------
C
C 1- SUBROUTINE STATISTICS:
C     NAME     : PIJABC
C     LEVEL    : 2 (CALLED BY 'EXCELP')
C     USE      : THIS ROUTINE RECONSTRUCT COLLISION PROBABILITIES (CP)
C                FOR  ALL ZONES ELIMINATING SURFACES FROM THE SYSTEM.
C     MODIFIED : 91-05-14 (G.M.)
C                91-07-12 (R.R.)
C                94-05-12 (I.P.) ( PIJK )
C                98-02-11 (G.M.) ( PERIODIC B.C. )
C     AUTHOR   : R. ROY (87-05-01)
C
C 2- PARAMETERS:
C  INPUT
C     NREG    : # OF ZONES FOR GEOMETRY.             I
C     NSOUT   : # OF SURFACES FOR GEOMETRY.          I
C     NPRB    : NUMBER OF PROBABILITIES IN PROB I
C     SIGTAL  : ALBEDO-SIGT VECTOR                   R(-NSOUT:NREG)
C     MATRT   : REFLECTION/TRANSMISSION VECTOR       I(NSOUT)
C     PROB    : -CP- MATRIX FOR ALL TYPES.           D(NPRB)
C  OUTPUT
C     PROB    : -CP- MATRIX FOR ALL TYPES.           D(NPRB)
C     PSST    : WORKING AREA FOR PSST                D(NSOUT,NSOUT)
C               PSST=(A**(-1)-PSS)**(-1)
C     PSVT    : WORKING AREA FOR PSVT                D(NSOUT,NREG)
C               PSVT=PSST*PSV
C
C  3-EXTERNAL ROUTINE CALLED
C     ALINVD  : INVERSE A DOUBLE PRECISION MATRIX
C     XABORT  : DRAGON ABORT ROUTINE
C
C----------------------------------------------------------------------
C
      IMPLICIT         NONE
C----
C VARIABLES
C-----
      INTEGER          NREG,NSOUT,NPRB,NSP1,IPSS,ISUR,ISX,IS,JSX,JS,IER,
     1                 IV,JV,IVSI,IVVI,ISV,IVS,IVV
      INTEGER          MATRT(NSOUT)
C
      REAL             SIGTAL(-NSOUT:NREG)
      DOUBLE PRECISION PROB(NPRB),PSST(NSOUT,NSOUT),PSVT(NSOUT,NREG)
C----
C  EVALUATE MATRIX (A**(-1)-PSS)
C----
      NSP1=NSOUT+1
      IPSS=0
      ISUR=(NSOUT*NSP1)/2
      ISX=0
      DO 100 IS=-NSOUT,-1,1
        ISX=ISX+1
        JSX=0
        ISUR=ISUR+1
        DO 101 JS=-NSOUT,IS,1
          JSX=JSX+1
          IPSS=IPSS+1
          IF((SIGTAL(IS).EQ.0.0).OR.(SIGTAL(JS).EQ.0.0)) THEN
            PSST(ISX,JSX)= 0.0D0
          ELSE
            PSST(ISX,JSX)=-PROB(IPSS)
          ENDIF
          IF(JS.NE.IS) THEN
            PSST(JSX,ISX)=PSST(ISX,JSX)
          ENDIF
 101    CONTINUE
        IF(SIGTAL(IS) .EQ. 0.0)THEN
          PSST(ISX,ISX)=PROB(ISUR)
        ELSE
          JS=-MATRT(-IS)
          IF(JS .EQ. IS) THEN
            PSST(ISX,ISX)=PSST(ISX,ISX)+PROB(ISUR)/SIGTAL(IS)
          ELSE IF(JS .LT. IS) THEN
            JSX=NSOUT+JS+1
            PSST(ISX,JSX)=PSST(ISX,JSX)+PROB(ISUR)/SIGTAL(IS)
            PSST(JSX,ISX)=PSST(ISX,JSX)
          ENDIF
        ENDIF
 100  CONTINUE
C----
C  INVERSE MATRIX PSST=(A**(-1)-PSS)
C----
      CALL ALINVD(NSOUT,PSST,NSOUT,IER)
C----
C  CHECK IF INVERSE IS VALID
C----
      IF(IER .NE. 0 ) CALL XABORT
     >  ('PIJABC: IMPOSSIBLE TO INVERT PSS COUPLING MATRIX')
      IVSI=(NSP1*(NSP1+1))/2
      IVVI=IVSI+NSP1
      DO 110 IV=1,NREG
C----
C    PSVT(IS,IV)=SUM(JSS) PSST(ISS,JSS)*PSV(JSS,IV)
C----
        DO 111 IS=1,NSOUT
          PSVT(IS,IV)=0.0D0
 111    CONTINUE
        DO 120 IS=1,NSOUT
          DO 121 JS=1,NSOUT
            ISV=IVSI+JS
            PSVT(IS,IV)=PSVT(IS,IV)+PSST(IS,JS)*PROB(ISV)
 121      CONTINUE
 120    CONTINUE
        IVV=IVVI
        DO 130 JV=1,IV
          IVV=IVV+1
          ISV=0
          IVS=IVSI
          DO 131 IS=-NSOUT,-1,1
            ISV=ISV+1
            IVS=IVS+1
            IF(SIGTAL(IS).NE.0.0) THEN
              PROB(IVV)=PROB(IVV)+PROB(IVS)*PSVT(ISV,JV)
            ENDIF
 131      CONTINUE
 130    CONTINUE
        IVSI=IVSI+NSP1+IV
        IVVI=IVVI+NSP1+IV
 110  CONTINUE
      RETURN
      END
