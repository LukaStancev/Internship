*DECK PIJCMP
      SUBROUTINE PIJCMP(NREG,NSOUT,NCOR,DPR,VOLSUR,LPIJK,PROB)
C
C-------------------------    PIJCMP    -------------------------------
C
C 1- SUBROUTINE STATISTICS:
C     NAME     : PIJCMP
C     LEVEL    : 2 (CALLED BY 'EXCELP')
C     USE      : COMPRESSION OF PIJ MATRICES IN SYMMETRIC FORMAT
C     AUTHOR   : R. ROY (94-05-03)
C
C 2- PARAMETERS:
C  INPUT
C     NREG    : TOTAL NUMBER OF REGIONS                I
C     NSOUT   : NUMBER OF OUTER SURFACE                I
C     NCOR    : MAXIMUM NUMBER OF CORNERS              I
C     DPR     : COLLISION PROBABILITIES      D(-NSOUT:NREG,-NSOUT:NREG)
C     VOLSUR  : VOLUMES                                R(-NSOUT:NREG)
C     LPIJK   : PIJK FLAG                              L
C
C  OUTPUT
C     PROB    : COMPRESS PROBABILITY MATRIX            D(NPLEN)
C               NPLEN=(NREG+NSOUT+2)*(NREG+NSOUT+1)/2
C               IND(I,J)=MAX(I+NSOUT+1,J+NSOUT+1)
C                       *(MAX(I+NSOUT+1,J+NSOUT+1)-1)/2
C                       +MIN(I+NSOUT+1,J+NSOUT+1)
C         IS=-NSOUT,-1; JS=-NSOUT,IS; I=IND(IS,JS)
C           PROB(I)=VOLSUR(IS)*PSS(IS,JS)
C         IV=1,NREG; JS=-NSOUT,-1;    I=IND(IV,JS)
C           SIGT(IV).GT.0.0
C             PROB(I)=SIGT(IV)*VOLSUR(IV)*PVS(IV,JS)
C           SIGT(IV).EQ.0.0
C             PROB(I)=VOLSUR(IV)*PVS(IV,JS)
C         IV=1,NREG; JV=1,IV;       I=IND(IV,JV)
C           SIGT(IV).GT.0.0 AND SIGT(JV).GT.0.0
C             PROB(I)=SIGT(IV)*SIGT(JV)*VOLSUR(IV)*PVV(IV,JV)
C           SIGT(IV).GT.0.0 AND SIGT(JV).EQ.0.0
C             PROB(I)=SIGT(IV)*VOLSUR(IV)*PVV(IV,JV)
C           SIGT(IV).EQ.0.0 AND SIGT(JV).GT.0.0
C             PROB(I)=SIGT(JV)*VOLSUR(IV)*PVV(IV,JV)
C           SIGT(IV).EQ.0.0 AND SIGT(JV).EQ.0.0
C             PROB(I)=VOLSUR(IV)*PVV(IV,JV)
C
C-------------------------    PIJCMP    -------------------------------
C
      IMPLICIT NONE
C----
C VARIABLES
C----
      INTEGER          NREG,NSOUT,NCOR,IPR,IL,JL,IVOL,IUN
      DOUBLE PRECISION DPR(-NSOUT:NREG,-NSOUT:NREG),PROB(*),ZCOR,ZCOR1
      REAL             VOLSUR(-NSOUT:NREG),COEF
      LOGICAL          LPIJK
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: WORK
C----
C  SCRATCH STORAGE ALLOCATION
C----
      ALLOCATE(WORK((NREG+NSOUT+2)*(NREG+NSOUT+1)/2))
C----
C  SYMMETRIZE AND STORE IN PROB
C----
      IPR= 0
      DO  150 IL    = -NSOUT, NREG
         DO  160 JL = -NSOUT, IL
            IPR= IPR+1
            WORK(IPR)= DPR(IL,JL) + DPR(JL,IL)
  160    CONTINUE
  150 CONTINUE
      IF( NCOR.EQ.1 )THEN
         IPR= 0
         DO  250 IL    = -NSOUT, NREG
            DO  260 JL = -NSOUT, IL
               IPR= IPR+1
               PROB(IPR)= WORK(IPR)
  260       CONTINUE
  250    CONTINUE
      ELSE
         IPR= 0
         ZCOR1= 1.0D0/DBLE(NCOR)
         ZCOR=  1.0D0/DBLE(NCOR*NCOR)
         DO  251 IL    = -NSOUT, NREG
            IF( IL.GT.0 ) ZCOR= ZCOR1
            DO  261 JL = -NSOUT, IL
               IPR= IPR+1
               IF( JL.GT.0 ) ZCOR= 1.0D0
               PROB(IPR)= ZCOR * WORK(IPR)
  261       CONTINUE
  251    CONTINUE
      ENDIF
C----
C  CHARGE VOLUMES IN THE PROB MATRIX
C----
      COEF=1.0
      IVOL= NSOUT*(NSOUT+1)/2
      DO 300 IUN= -NSOUT, NREG
         IF( IUN.LE.0 )THEN
            IVOL= IVOL+1
            IF(LPIJK) COEF= 3./4.
         ELSE
            IVOL= IVOL+NSOUT+IUN
            IF(LPIJK) COEF= 2./3.
         ENDIF
         IF( PROB(IVOL).NE.0.0 )THEN
            CALL XABORT( 'PIJCMP: UNEXPECTED VALUE IN PROB MATRIX' )
         ENDIF
         PROB(IVOL) = VOLSUR(IUN)*COEF
  300 CONTINUE
C----
C  SCRATCH STORAGE DEALLOCATION
C----
      DEALLOCATE(WORK)
      RETURN
      END
