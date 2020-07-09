*DECK PIJS3D
      SUBROUTINE PIJS3D(NREG,NSOUT,NSLINE,NSBG,WEIGHT,
     >                  RCUTOF,SIGTAL,NPSYS,
     >                  SEGLEN,NRSEG,
     >                  STAYIN,GOSOUT,DPR)
C
C-------------------------    PIJS3D    -------------------------------
C
C 1- SUBROUTINE STATISTICS:
C     NAME     : PIJS3D
C     LEVEL    : 2 (CALLED BY 'EXCELP')
C     USE      : INTEGRATION FOR GENERAL 3D SPECULAR  B.C. TRACK FILE
C     MODIFIED : 91/07/12 (R.R.)
C                98/02/10 (G.M.)
C                MODIFIED BECAUSE OF SURFACE DOUBLING IN XELS2D
C                INSERTED TO TAKE INTO ACCOUNT
C                PERIODIC BOUNDARY CONDITIONS
C     AUTHOR   : R. ROY
C     REFERENCE: 'A CYCLIC TRACKING PROCEDURE FOR CP CALCULATIONS
C                 IN 2-D LATTICES', R.ROY ET AL.,
C                 CONF/ADVANCES IN MATH, COMP & REACTOR PHYSICS,
C                 PITTSBURGH/USA, V 1, P 2.2 4-1 (1991).
C
C 2- PARAMETERS:
C  INPUT
C     NREG    : TOTAL NUMBER OF REGIONS                I
C     NSOUT   : NUMBER OF OUTER SURFACE                I
C     NSLINE  : number of segemnts on line             
C     NSBG    : NUMBER OF SUBGROUP                     I
C     WEIGHT  : line weight 
C     RCUTOF  : MFP CUT-OFF FACTOR (TRUNCATE LINES)    R
C     SIGTAL  : ALBEDO-CROSS SECTION VECTOR            R(-NSOUT:NREG,
C               IS=-NSOUT,-1; SIGTAL(IS)=ALBEDO(-IS)     NSBG)
C               IV=1,NREG ;   SIGTAL(IV)=SIGT(IV)
C     NPSYS   : NON-CONVERGED ENERGY GROUP INDICES.    I(NSBG)
C     SEGLEN  : LENGTH OF TRACK              D(NSLINE)
C     NRSEG   : REGION CROSSED BY TRACK      I(NSLINE)
C  OUTPUT
C     DPR     : COLLISION PROBABILITIES      D(-NSOUT:NREG,
C                                              -NSOUT:NREG,NSBG)
C               VALUES ARE COMPUTED ONLY FOR ONE DIRECTION P(I->J)
C  WORK
C     STAYIN  : STAY-IN ZONE PROBABILITY     D(NSLINE)
C     GOSOUT  : GOES-OUT ZONE PROBABILITY    D(NSLINE)
C
C-------------------------    PIJS3D    -------------------------------
C
      IMPLICIT         NONE
C----
C VARIABLES
C----
      INTEGER          NREG,NSOUT,NSLINE,NSBG
      INTEGER          NRSEG(NSLINE),NPSYS(NSBG)
      REAL             RCUTOF,SIGTAL(-NSOUT:NREG,NSBG)
      DOUBLE PRECISION WEIGHT,SEGLEN(NSLINE),STAYIN(NSLINE),
     >                 GOSOUT(NSLINE)
      DOUBLE PRECISION DPR(-NSOUT:NREG,-NSOUT:NREG,NSBG)
C----
C  Local variables
C----
      INTEGER          ISBG,IL,JL,NOIL,NOJL,ISODD,JSODD,IJDEL
      DOUBLE PRECISION TTOT,XSIL,OPATH,FINV,CUTOF
      REAL             ZERO,ONE,HALF
      PARAMETER       (ZERO=0.0E0, ONE=1.0E0, HALF=0.5E0 )
      REAL             SIXT,CUTEXP
      PARAMETER       (SIXT=HALF/3.0,CUTEXP=0.02)
      DOUBLE PRECISION EXSIL,XSIL2
      TTOT= ONE
      DO 2001 ISBG=1,NSBG
        IF(NPSYS(ISBG).EQ.0) GO TO 2001
C
C1.1)    CHANGE PATHS => GOSOUT AND STAYIN PATHS, INCLUDING ALBEDOS
C        ADD *PII* LOCAL NON-CYCLIC CONTRIBUTIONS
        ISODD=0
        DO 30 IL= 1, NSLINE
          NOIL  = NRSEG(IL)
          IF( NOIL.LT.0 )THEN
            IF(ISODD .EQ. 1) THEN
              ISODD=0
C----
C  FOR SURFACES:
C    OLD VERSION BEFORE SURFACE DOUBLING
C      GOSOUT= ALBEDO * SURFACE WEIGHT
C      WHERE ALL SURFACE WEIGHTS WERE 1.0
C    NEW VERSION WITH SURFACE DOUBLING
C      GOSOUT= ALBEDO
C    STAYIN = 1- ALBEDO * SURFACE WEIGHT
C    TTOT   = PRODUCT OF GOSOUT
C----
              GOSOUT(IL)= SIGTAL(NOIL,ISBG)
              STAYIN(IL)= ONE - GOSOUT(IL)
              TTOT= TTOT * GOSOUT(IL)
            ELSE
              ISODD=1
C----
C  FOR SURFACES:
C    OLD VERSION BEFORE SURFACE DOUBLING
C      GOSOUT= ALBEDO * SURFACE WEIGHT
C      WHERE ALL SURFACE WEIGHTS WERE 1.0
C    NEW VERSION WITH SURFACE DOUBLING
C      GOSOUT= ALBEDO
C    STAYIN = 1- ALBEDO * SURFACE WEIGHT
C    TTOT   = PRODUCT OF GOSOUT
C----
              GOSOUT(IL)= SIGTAL(NOIL,ISBG)
              STAYIN(IL)= ONE
            ENDIF
          ELSE
C----
C  FOR REGIONS
C  STAYIN = 1 -  EXP[ -CROSS SECTION * LENGTH OF NSLINE]
C  GOSOUT = 1 -  STAYIN
C  TTOT   = PRODUCT OF GOSOUT
C----
            XSIL  = SIGTAL(NOIL,ISBG)
            IF( XSIL .EQ. ZERO) THEN
              GOSOUT(IL)= ONE
              STAYIN(IL)= SEGLEN(IL)
              DPR(NOIL,NOIL,ISBG)= DPR(NOIL,NOIL,ISBG)
     >                           + HALF*WEIGHT*STAYIN(IL)*STAYIN(IL)
            ELSE IF( XSIL .LT. CUTEXP) THEN
              OPATH= SIGTAL(NOIL,ISBG)*SEGLEN(IL)
              XSIL2=OPATH*OPATH
              EXSIL=XSIL2*(HALF-SIXT*OPATH+XSIL2/24.0)
              STAYIN(IL)=OPATH-EXSIL
              GOSOUT(IL)= ONE - STAYIN(IL)
              TTOT= TTOT * GOSOUT(IL)
              DPR(NOIL,NOIL,ISBG)= DPR(NOIL,NOIL,ISBG) + WEIGHT*EXSIL
            ELSE
              OPATH= SIGTAL(NOIL,ISBG)*SEGLEN(IL)
              EXSIL= EXP(-OPATH)
              STAYIN(IL)= ONE - EXSIL
              GOSOUT(IL)= EXSIL
              TTOT= TTOT * GOSOUT(IL)
              DPR(NOIL,NOIL,ISBG)= DPR(NOIL,NOIL,ISBG)
     >                           + WEIGHT*(OPATH-STAYIN(IL))
            ENDIF
          ENDIF
 30     CONTINUE
C
C1.2)    COMPUTE CYCLIC FACTORS BY ANGLE
C        USING GLOBAL TRACK ATTENUATION: BETA(TOT)*EXP(-MFP(TOT))
         IF(TTOT .GE. ONE )THEN
           CALL XABORT( 'PIJS3D: ALBEDOS ARE NOT COMPATIBLE')
         ENDIF
         FINV= WEIGHT / (ONE-TTOT)
C
C1.3)    ADD *PIJ* CONTRIBUTIONS FOR FORWARD SOURCES
        ISODD=0
        DO 50 IL= 1, NSLINE
          NOIL  = NRSEG(IL)
          TTOT= FINV * STAYIN(IL)
          CUTOF= RCUTOF*TTOT
          IF( NOIL .LT. 0) THEN
            ISODD=MOD(ISODD+1,2)
            JSODD=ISODD
            DO 70 IJDEL= 1, NSLINE
              JL= MOD(IL+IJDEL-1,NSLINE) + 1
              NOJL=NRSEG(JL)
              IF( NOJL .LT. 0 ) THEN
                JSODD=MOD(JSODD+1,2)
                IF( ISODD.EQ.1 .AND. JSODD .EQ.0) THEN
                  DPR(NOJL,NOIL,ISBG)= DPR(NOJL,NOIL,ISBG)
     >                               + TTOT * STAYIN(JL)
                  TTOT= TTOT * GOSOUT(JL)
                  IF( TTOT.LE.CUTOF ) GO TO 55
                ENDIF
              ELSE IF(ISODD.EQ.1) THEN
                DPR(NOJL,NOIL,ISBG)= DPR(NOJL,NOIL,ISBG)
     >                             + TTOT * STAYIN(JL)
                TTOT= TTOT * GOSOUT(JL)
                IF( TTOT.LE.CUTOF ) GO TO 55
              ENDIF
 70         CONTINUE
          ELSE
            JSODD=ISODD
            DO 80 IJDEL= 1, NSLINE
              JL= MOD(IL+IJDEL-1,NSLINE) + 1
              NOJL=NRSEG(JL)
              IF( NOJL .LT. 0 ) THEN
                JSODD=MOD(JSODD+1,2)
                IF( JSODD .EQ.0) THEN
                  DPR(NOJL,NOIL,ISBG)= DPR(NOJL,NOIL,ISBG)
     >                               + TTOT * STAYIN(JL)
                  TTOT= TTOT * GOSOUT(JL)
                  IF( TTOT.LE.CUTOF ) GO TO 55
                ENDIF
              ELSE
                DPR(NOJL,NOIL,ISBG)= DPR(NOJL,NOIL,ISBG)
     >                             + TTOT* STAYIN(JL)
                TTOT= TTOT * GOSOUT(JL)
                IF( TTOT.LE.CUTOF ) GO TO 55
              ENDIF
 80         CONTINUE
          ENDIF
 55       CONTINUE
 50     CONTINUE
 2001 CONTINUE
      RETURN
      END
