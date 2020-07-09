*DECK QIJCMP
      SUBROUTINE QIJCMP(NREG,NSOUT,NPIJ,NGRP,NCOR,VOLSUR,SIGTAL,DPR,
     > NPSYS)
C
C-------------------------    QIJCMP    -------------------------------
C
C 1- SUBROUTINE STATISTICS:
C     NAME     : QIJCMP
C     LEVEL    : 2 (CALLED BY 'XL3')
C     USE      : COMPRESSION OF PIJ MATRICES IN SYMMETRIC FORMAT
C     AUTHOR   : R. ROY (96-04-15)
C
C 2- PARAMETERS:
C  INPUT
C     NREG    : TOTAL NUMBER OF REGIONS                I
C     NSOUT   : NUMBER OF OUTER SURFACE                I
C     NPIJ    : NUMBER OF PROBABILITIES IN ONE GROUP   I
C     NGRP    : NUMBER OF ENERGY GROUPS                I
C     NCOR    : MAXIMUM NUMBER OF CORNERS              I
C     DPR     : COLLISION PROBABILITIES                D(NGRP,*)
C     NPSYS   : NON-CONVERGED ENERGY GROUP INDICES.    I(NGRP)
C
C  OUTPUT
C     DPR     : COMPRESS PROBABILITY MATRIX            D(NUN**2,NGRP)
C               NPLEN=(NREG-NSOUT+2)*(NREG-NSOUT+1)/2
C               IND(I,J)=MAX(I-NSOUT+1,J-NSOUT+1)
C                       *(MAX(I-NSOUT+1,J-NSOUT+1)-1)/2
C                       +MIN(I-NSOUT+1,J-NSOUT+1)
C         IS=NSOUT,-1; JS=NSOUT,IS; I=IND(IS,JS)
C           PROB(I)=VOLSUR(IS)*PSS(IS,JS)
C         IV=1,NREG; JS=NSOUT,-1;    I=IND(IV,JS)
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
C-------------------------    QIJCMP    -------------------------------
C
      IMPLICIT         NONE
      INTEGER          NREG,NSOUT,NPIJ,NGRP,NCOR,NPSYS(NGRP)
      INTEGER          IPR,IL,JL,IG,INDIJ
      REAL             VOLSUR(NSOUT:NREG),SIGTAL(NSOUT:NREG,NGRP),ZERO
      DOUBLE PRECISION DPR(NPIJ,NGRP),ZCOR,ZCOR1,DZERO
      PARAMETER      ( ZERO=0.0, DZERO=0.0D0 )
C----
C  SYMMETRIZE AND STORE IN PROB
C----
      INDIJ= 0
      DO 5 IL = 1, NREG-NSOUT+1
         INDIJ= INDIJ + IL
         DO 1 IG= 1, NGRP
           IF(NPSYS(IG).NE.0)
     >     DPR(INDIJ,IG)= DPR(INDIJ,IG) + DPR(INDIJ,IG)
    1    CONTINUE
    5 CONTINUE
      IF( NCOR.NE.1 )THEN
         ZCOR1= 1.0D0/DBLE(NCOR)
         ZCOR=  1.0D0/DBLE(NCOR*NCOR)
         INDIJ= 0
         DO 35 IL    = NSOUT, NREG
            IF( IL.GT.0 ) ZCOR= ZCOR1
            DO 25 JL = NSOUT, IL
               INDIJ= INDIJ + 1
               IF( JL.GT.0 ) ZCOR= 1.0D0
               DO 15 IG= 1, NGRP
                  IF(NPSYS(IG).NE.0)
     >            DPR(INDIJ,IG)= ZCOR * DPR(INDIJ,IG)
   15          CONTINUE
   25       CONTINUE
   35    CONTINUE
      ENDIF
      IPR=-((1-NSOUT)*NSOUT)/2
      DO 80 IL= NSOUT,-1
         IPR= IPR+1
         DO 70 IG= 1, NGRP
            IF(NPSYS(IG).NE.0) DPR(IPR,IG)= DBLE(VOLSUR(IL))
   70    CONTINUE
   80 CONTINUE
      IPR= IPR+1
      DO 90 IG= 1, NGRP
         DPR(IPR,IG)= DZERO
   90 CONTINUE
      DO 110 IL= 1,NREG
         IPR= IPR-NSOUT+IL
         DO 100 IG= 1, NGRP
            IF(NPSYS(IG).EQ.0) GO TO 100
            IF( SIGTAL(IL,IG).EQ.ZERO )THEN
               DPR(IPR,IG)= DBLE(VOLSUR(IL))
            ELSE
               DPR(IPR,IG)= DBLE(VOLSUR(IL)*SIGTAL(IL,IG))
            ENDIF
  100    CONTINUE
  110 CONTINUE
C
      RETURN
      END
