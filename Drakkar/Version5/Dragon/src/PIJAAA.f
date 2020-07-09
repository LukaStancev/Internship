*DECK PIJAAA
      SUBROUTINE PIJAAA(NREG,NSOUT,SIGTAL,PROB,PSVT,PROBS)
C
C-------------------------    PIJAAA    -------------------------------
C
C 1- SUBROUTINE STATISTICS:
C     NAME     : PIJAAA
C     LEVEL    : 2 (CALLED BY 'EXCELP')
C     USE      : THIS ROUTINE CALCULATES DIRECTIONAL" COLLISION
C                PROBABILITIES FOR ALL ZONES ELIMINATING
C                SURFACES FROM THE SYSTEM:
C                PIJK"=PIJK+PISK*((1.-PSS)**(-1.))*PSJ
C     IMPORTANT: THIS ROUTINE HAS TO BE AFTER PIJABC(PROB) AND
C                BEFORE PIJABC(PROBX)
C     AUTHOR   : I.PETROVIC (94-05-18)
C
C 2- PARAMETERS:
C  INPUT
C     NREG    : # OF ZONES FOR GEOMETRY.             I
C     NSOUT   : # OF SURFACES FOR GEOMETRY.          I
C     NPRB    : NUMBER OF PROBABILITIES IN PROB      I
C     SIGTAL  : ALBEDO-SIGT VECTOR                   R(-NSOUT:NREG)
C     PROB    : DIRECTIONAL CP MATRIX FOR ALL TYPES  D(NPRB)
C     PSVT    : PSST MATRIX                          D(NSOUT,NREG)
c               PSVT=(A**(-1)-PSS)**(-1)*PSV
C  OUTPUT
C     PROBS   : DIRECTIONAL" CP MATRIX               R(3*NNREG)
C
C
C----------------------------------------------------------------------
C
      IMPLICIT NONE
C----
C INTERFACE VARIABLES
C----
      INTEGER          NREG,NSOUT
      REAL             SIGTAL(-NSOUT:NREG),PROBS(*)
      DOUBLE PRECISION PROB(*),PSVT(NSOUT,NREG)
C----
C LOCAL VARIABLES
C----
      INTEGER NSP1,IVSI,IDPSV,IV,IPRL,IPRU,JV,ISV,IPSV,IVS,ISU
C
      NSP1=NSOUT+1
      IVSI=(NSP1*(NSP1+1))/2
      IDPSV=IVSI
      DO 100 IV=1,NREG
        IPRL=NREG*(IV-1)+1
        IPRU=IV
        IPSV=IDPSV
        DO 110 JV=1,IV
          ISV=0
          IVS=IVSI
          DO 120 ISU=-NSOUT,-1,1
            ISV=ISV+1
            IVS=IVS+1
            IPSV=IPSV+1
            IF(SIGTAL(ISU).NE.0.0) THEN
              PROBS(IPRL)=PROBS(IPRL)+REAL(PROB(IVS)*PSVT(ISV,JV))
              IF(IPRL.NE.IPRU) THEN
                PROBS(IPRU)=PROBS(IPRU)+REAL(PROB(IPSV)*PSVT(ISV,IV))
              ENDIF
            ENDIF
 120      CONTINUE
          IPSV=IPSV+JV+1
          IPRL=IPRL+1
          IPRU=NREG*JV+IV
 110    CONTINUE
        IVSI=IVSI+NSP1+IV
 100  CONTINUE
C
      RETURN
      END
