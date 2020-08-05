*DECK PIJD2R
      SUBROUTINE PIJD2R(NREG,NSOUT,PROB,FACTOR,LPIJK,NELPIJ,N2PROB,PIJ)
C
C-------------------------    PIJD2R    -------------------------------
C
C 1- SUBROUTINE STATISTICS:
C     NAME     : PIJD2R
C     LEVEL    : 2 (CALLED BY 'EXCELP')
C     USE      : CHARGE PIJ MATRICES IN THE DRAGON SYMMETRIZED FORMAT
C     AUTHOR   : A. HEBERT
C
C 2- PARAMETERS:
C  INPUT
C     NREG    : TOTAL NUMBER OF REGIONS                I
C     NSOUT   : NUMBER OF OUTER SURFACE                I
C     PROB    : COLLISION PROBABILITIES                D(NPLEN)
C               NPLEN=(NREG+NSOUT+2)*(NREG+NSOUT+1)/2
C     FACTOR  : ONE OVER TOTAL XS                      R(NREG)
C     LPIJK   : PIJK FLAG                              L
C     NELPIJ  : NUMBER OF TERMS IN PIJ                 I
C     N2PROB  : NUMBER OF TERMS IN PROB                I
C
C  OUTPUT
C     PIJ     : SYMMETRIC PROBABILITY MATRIX           R(NELPIJ)
C
C-------------------------    PIJD2R    -------------------------------
C
      IMPLICIT NONE
C----
C VARIABLES
C----
      INTEGER          NREG,NSOUT,NELPIJ,N2PROB,IUN,JUN,KPRB,IVV
      DOUBLE PRECISION PROB(N2PROB)
      REAL             FACTOR(NREG),PIJ(NELPIJ),COEF
      LOGICAL          LPIJK
C----
C  STORE IN SYMMETRIC FORMAT
C----
      IVV=0
      COEF=1.0
      IF(LPIJK) COEF=1.5
      KPRB=(NSOUT+1)*(NSOUT+2)/2+NSOUT+1
      DO 20 IUN=1,NREG
         DO 10 JUN=1,IUN
            KPRB=KPRB+1
            IVV=IVV+1
            PIJ(IVV)=COEF*REAL(PROB(KPRB))*FACTOR(IUN)*FACTOR(JUN)
   10    CONTINUE
         KPRB=KPRB+NSOUT+1
   20 CONTINUE
C
      RETURN
      END
