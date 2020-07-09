*DECK PIJRHL
      SUBROUTINE PIJRHL(IPRT,NREG,NSOUT,SIGTAL,PROB)
C
C-------------------------    PIJRHL    -------------------------------
C
C 1- SUBROUTINE STATISTICS:
C     NAME     : PIJRHL
C     LEVEL    : 2 (CALLED BY 'EXCELP')
C     USE      : NON-LINEAR NORMALIZATION OF COLLISION PROBS
C     AUTHOR   : R. ROY (94-04-18)
C     MODIFIED : E. VARIN (97-01-16)
C     REFERENCE: 'NORMALIZATION TECHNIQUES FOR CP MATRICES',
C                 R.ROY AND G.MARLEAU,
C                 CONF/PHYSOR-90, MARSEILLE/FRANCE, V 2, P IX-40 (1990).
C                 'HELIOS: ANGULARLY DEPENDENT COLLISION PROBABILITIES'
C                 E.A. VILLARINO, R.J.J.STAMM'LER
C                 AND A.A.FERRI AND J.J.CASAL
C                 Nucl.Sci.Eng. 112,16-31, 1992.
C
C 2- PARAMETERS:
C  INPUT
C     NREG    : # OF ZONES FOR GEOMETRY.             I
C     NSNEG   : # OF SURFACES FOR GEOMETRY.          I
C     SIGTAL  : ALBEDO-SIGT VECTOR                   R(-NSOUT:NREG)
C     PROB    : -CP- MATRIX FOR ALL TYPES.           D(NPRB)
C               NPRB=(NSOUT+NREG+1)*(NSOUT+NREG+2)/2
C  OUTPUT
C     PROB    : -CP- MATRIX FOR ALL TYPES.           D(NPRB)
C
C  3-INTERNAL PARAMETERS
C     EPSCON  : CONVERGENCE CRITERIA = 1.0E-6
C     NITMAX  : MAXIMUM NUMBER OF ITERATIONS = 20
C
C-------------------------    PIJRHL    -------------------------------
C
      IMPLICIT   NONE
      INTEGER    IPRT,NREG,NSOUT,IUNOUT,NITMAX,NIT,IPRINT,
     >           IR,JR,IP,IPRB,IND,I,J,CPTLB,CPTAC,CTOT,
     >           NSURC,NSURM,NVOLC,NVOLM
      REAL       SIGTAL(-NSOUT:NREG)
      LOGICAL    NOTCON
      DOUBLE PRECISION PROB(*),NOM,DENOM,DMU,WFSPAD,WFSP,EPSCON,R1,R2,
     >           TOTCON,TMPCON
      PARAMETER (IUNOUT=6, IPRINT=10, EPSCON=1.0E-6, NITMAX=20)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: CHI
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: WEIG
C
C----- INTRINSIC FUNCTION FOR POSITION IN CONDENSE PIJ MATRIX
C
      IND(I,J)=(MAX(I+NSOUT+1,J+NSOUT+1)*
     >         (MAX(I+NSOUT+1,J+NSOUT+1)-1))/2
     >         +MIN(I+NSOUT+1,J+NSOUT+1)
C----
C  SCRATCH STORAGE ALLOCATION
C   WEIG : ADDITIVE WEIGHT
C----
      ALLOCATE(WEIG(-NSOUT:NREG,3),CHI(-NSOUT:NREG))
C
      NOTCON= .FALSE.
      CPTLB = 3
      CPTAC = 3
      CTOT = CPTAC+CPTLB
C
C     INITIALISATION OF WEIGHTS
      DO 60 IR=-NSOUT, NREG
         WEIG(IR,1)=0.0D0
         WEIG(IR,2)=0.5D0
         WEIG(IR,3)=0.5D0
   60 CONTINUE
      DO 50 IR=-NSOUT, NREG
         CHI(IR)= 1.0D0
         IF( IR.GE.0.AND.SIGTAL(IR).EQ.0.0D0 )THEN
            CHI(IR)= 0.0D0
         ENDIF
   50 CONTINUE
C
C
C     MAIN ITERATION LOOP
      IF(IPRT.GT.2) WRITE(IUNOUT,'(A24)')
     >       'ITER.     MU      ERROR '
      DO 110 NIT=1,NITMAX
C
         DO 220 IR= -NSOUT, NREG
            WFSPAD = PROB(IND(IR,0))
     >                   + CHI(IR)*PROB(IND(IR,IR))*WEIG(IR,3)
            WFSP = CHI(IR)*PROB(IND(IR,IR))
            DO 200 JR=-NSOUT, NREG
               WFSPAD = WFSPAD - CHI(JR)*WEIG(JR,3)*PROB(IND(IR,JR))
               WFSP = WFSP + CHI(JR)*PROB(IND(IR,JR))
  200       CONTINUE
            WEIG(IR,3) = WFSPAD / WFSP
  220    CONTINUE
C
C        ACCELERATION TECHNIQUE
         IF(  MOD(NIT-1,CTOT).GE.CPTAC )THEN
            NOM   = 0.0D0
            DENOM = 0.0D0
            DO 10 IR=-NSOUT, NREG
               R1= WEIG(IR,2) - WEIG(IR,1)
               R2= WEIG(IR,3) - WEIG(IR,2)
               NOM = NOM + R1*(R2-R1)
               DENOM = DENOM + (R2-R1)*(R2-R1)
   10       CONTINUE
            IF(DENOM.EQ.0.0D0) THEN
              DMU = 1.0D0
            ELSE
              DMU = - NOM / DENOM
            ENDIF
            IF( DMU.GT.10.0D0 .OR. DMU.LT.0.0D0 )CALL XABORT('PIJRHL: '
     >         //'PROBLEM OF ACCELERATION')
            DO 20 IR=-NSOUT, NREG
               WEIG(IR,3) = WEIG(IR,2) + DMU *
     >                           (WEIG(IR,3) - WEIG(IR,2))
               WEIG(IR,2) = WEIG(IR,1) + DMU *
     >                           (WEIG(IR,2) - WEIG(IR,1))
   20       CONTINUE
         ELSE
            DMU = 1.0D0
         ENDIF
C
C        CALCULATIONS OF SQUARE DISTANCE BETWEEN 2 ITERATIONS
C        AND UPDATING THE SOLUTION
         TOTCON = 0.0D0
         DO 100 IR=-NSOUT, NREG
            TMPCON=ABS(WEIG(IR,3)-WEIG(IR,2))/WEIG(IR,3)
            TOTCON=MAX(TMPCON,TOTCON)
            WEIG(IR,1)= WEIG(IR,2)
            WEIG(IR,2)= WEIG(IR,3)
  100    CONTINUE
         IF( IPRT.GT.2 ) WRITE(IUNOUT,'(I3,F9.5,E15.7)') NIT,DMU,TOTCON
C
C        CONVERGENCE TEST
         IF( TOTCON.LT.EPSCON )GO TO 120
C
  110 CONTINUE
      NOTCON=.TRUE.
      WRITE(IUNOUT,'(35H PIJRHL: WEIGHTS NOT CONVERGED          )')
  120 CONTINUE
C
C     RENORMALIZE "PIJ" SYMMETRIC MATRIX
      IPRB = 0
      DO 240 IR   = -NSOUT, NREG
         DO 230 JR= -NSOUT, IR
            IPRB= IPRB+1
            IF( IR.NE.0.AND.JR.NE.0 )THEN
                PROB(IPRB)=PROB(IPRB)*(WEIG(IR,1)+WEIG(JR,1))
            ENDIF
  230    CONTINUE
  240 CONTINUE
C
C     PRINT WEIGHT FACTORS IF THERE IS A PROBLEM...
      IF( NOTCON .OR. IPRT.GE.IPRINT )THEN
         WRITE(IUNOUT,'(30H0 SURFACE WEIGHTS FACTORS                /)')
         NSURC = -1
         DO 300 IP  = 1, (9 +NSOUT) / 10
            NSURM= MAX( -NSOUT, NSURC-9 )
            WRITE(IUNOUT,'(10X,10( A5,    I6)/)')
     >                     (' SUR ',-IR,IR= NSURC, NSURM, -1)
            WRITE(IUNOUT,'(10H WEIGHT   ,10F11.5)')
     >                   (WEIG(IR,1),IR=NSURC,NSURM,-1)
            NSURC = NSURC - 10
 300     CONTINUE
         WRITE(IUNOUT,'(30H0  VOLUME WEIGHTS FACTORS                /)')
         NVOLC =  1
         DO 310 IP  = 1, (9 + NREG) / 10
            NVOLM= MIN( NREG, NVOLC+9 )
            WRITE(IUNOUT,'(10X,10( A5 ,  I6)/)')
     >                  (' VOL ',IR,IR=NVOLC,NVOLM, 1)
            WRITE(IUNOUT,'(10H WEIGHT   ,10F11.5)')
     >                   (WEIG(IR,1),IR=NVOLC,NVOLM, 1)
            NVOLC = NVOLC + 10
 310     CONTINUE
      ENDIF
C----
C  SCRATCH STORAGE DEALLOCATION
C----
      DEALLOCATE(CHI,WEIG)
      RETURN
      END
