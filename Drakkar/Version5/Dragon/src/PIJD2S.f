*DECK PIJD2S
      SUBROUTINE PIJD2S(NREG,NSOUT,PROB,PROBKS)
C
C-------------------------    PIJD2S    -------------------------------
C
C 1- SUBROUTINE STATISTICS:
C     NAME     : PIJD2S
C     LEVEL    : 2 (CALLED BY 'EXCELP')
C     USE      : CHARGE PROBKS MATRICES IN THE DRAGON SQUARE FORMAT
C     AUTHOR   : A. HEBERT
C
C 2- PARAMETERS:
C  INPUT
C     NREG    : TOTAL NUMBER OF REGIONS                I
C     NSOUT   : NUMBER OF OUTER SURFACE                I
C     PROB    : COLLISION PROBABILITIES                D(NPLEN)
C               NPLEN=(NREG+NSOUT+2)*(NREG+NSOUT+1)/2
C
C  OUTPUT
C     PROBKS  : SQUARE PROBABILITY MATRIX              R(NREG*NREG)
C
C-------------------------    PIJD2S    -------------------------------
C
      IMPLICIT NONE
C----
C VARIABLES
C----
      INTEGER          NREG,NSOUT,KPRB,IIU,IIL,II,JJ
      DOUBLE PRECISION PROB(((NSOUT+NREG+2)*(NSOUT+NREG+1))/2)
      REAL             PROBKS(NREG*NREG)
C----
C  STORE IN SQUARE FORMAT
C----
      KPRB=(NSOUT+1)*(NSOUT+2)/2+NSOUT+1
      DO 20 JJ=1,NREG
         IIU=JJ
         IIL=(JJ-1)*NREG+1
         DO 10 II=1,JJ
            KPRB=KPRB+1
            PROBKS(IIL)=REAL(PROB(KPRB))
            PROBKS(IIU)=PROBKS(IIL)
            IIU=JJ+II*NREG
            IIL=IIL+1
   10    CONTINUE
         KPRB=KPRB+NSOUT+1
   20 CONTINUE
C
      RETURN
      END
