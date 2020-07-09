*DECK SYBJJ0
      SUBROUTINE SYBJJ0 (IPAS,NSUPCE,NPIJ,NUNKNO,EPSJ,FUNKNO,SUNKNO,
     1 IMPX,ISTAT,NMC,PROCEL,PIJW,PISW,PSJW,PSSW)
*
*-----------------------------------------------------------------------
*
*Purpose:
* compute the neutron flux and interface currents in a do-it-yourself
* geometry using the current iteration method.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPAS    total number of regions.
* NSUPCE  number of cells.
* NPIJ    length of cellwise scattering-reduced collision probability
*         matrices.
* EPSJ    stopping criterion for flux-current iterations.
* NUNKNO  total number of unknowns in vectors SUNKNO and FUNKNO.
* SUNKNO  input source vector.
* IMPX    print flag (equal to 0 for no print).
* ISTAT   statistical approximation flag (set with ISTAT=1).
* NMC     offset of the first volume in each cell.
* PROCEL  user supplied geometrical matrix.
* PIJW    cellwise scattering-reduced collision probability matrices.
* PISW    cellwise scattering-reduced escape probability matrices.
* PSJW    cellwise scattering-reduced collision probability matrices
*         for incoming neutrons.
* PSSW    cellwise scattering-reduced transmission probability matrices.
*
*Parameters: input/output
* FUNKNO  unknown vector.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IPAS,NSUPCE,NPIJ,NUNKNO,IMPX,ISTAT,NMC(NSUPCE+1)
      REAL EPSJ,FUNKNO(NUNKNO),SUNKNO(NUNKNO),PROCEL(NSUPCE,NSUPCE),
     1 PIJW(NPIJ),PISW(IPAS),PSJW(IPAS),PSSW(NSUPCE)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION S1,S2,S3,S4,S5,DET
      LOGICAL LOGTES
      PARAMETER (MAXIT=400,LACCFC=2,ICL1=3,ICL2=3)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, DIMENSION(:), POINTER :: INDPIJ
      REAL, DIMENSION(:), POINTER :: CIT0
      REAL, DIMENSION(:,:), POINTER :: CITR,AITR
      DOUBLE PRECISION, DIMENSION(:), POINTER :: WCURR
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(INDPIJ(NSUPCE))
      ALLOCATE(CITR(3,NSUPCE),CIT0(NSUPCE),AITR(2,NSUPCE))
      ALLOCATE(WCURR(NSUPCE))
*
      IPIJ=0
      DO 10 JKG=1,NSUPCE
      J2=NMC(JKG+1)-NMC(JKG)
      INDPIJ(JKG)=IPIJ
   10 IPIJ=IPIJ+J2*J2
*----
*  PROCESS STATISTICAL APPROXIMATION
*----
      IF(ISTAT.NE.0) THEN
         X1=0.0
         DO 20 IKK=1,NSUPCE
   20    X1=X1+PSSW(IKK)*PROCEL(1,IKK)
         X1=1.0/(1.0-X1)
         SSS=0.0
         DO 30 IKK=1,NSUPCE
         I1P=NMC(IKK)
         I2=NMC(IKK+1)-I1P
         DO 30 I=1,I2
   30    SSS=SSS+PROCEL(1,IKK)*X1*PSJW(I1P+I)*SUNKNO(I1P+I)
         IT3=1
         DO 40 IKK=1,NSUPCE
   40    CITR(IT3,IKK)=SSS
         GO TO 190
      ENDIF
*----
*  COMPUTE PSJW * Q(*) CONTRIBUTION
*----
      DO 50 IKK=1,NSUPCE
      CIT0(IKK)=0.0
      CITR(1,IKK)=FUNKNO(IPAS+IKK)
      DO 50 JKK=1,NSUPCE
      I1P=NMC(JKK)
      I2=NMC(JKK+1)-I1P
      DO 50 I=1,I2
   50 CIT0(IKK)=CIT0(IKK)+PROCEL(IKK,JKK)*PSJW(I1P+I)*SUNKNO(I1P+I)
*----
*  COMPUTE NORMALIZATION VECTOR WCURR
*----
      DO 60 JKK=1,NSUPCE
      WCURR(JKK)=1.0D0
      DO 60 IKK=1,NSUPCE
   60 WCURR(JKK)=WCURR(JKK)-PROCEL(IKK,JKK)*PSSW(JKK)
*
      ISTART=1
      TEST=0.0
      ITER=0
   70 ITER=ITER+1
      IF(ITER.GT.MAXIT) THEN
         WRITE(6,'(/47H SYBJJ0: *** WARNING *** MAXIMUM NUMBER OF ITER,
     1   15HATIONS REACHED.)')
         GO TO 190
      ENDIF
      IT3=MOD(ITER,3)+1
      IT2=MOD(ITER-1,3)+1
      IT1=MOD(ITER-2,3)+1
      DO 80 I=1,NSUPCE
   80 CITR(IT3,I)=CIT0(I)
*----
*  COMPUTE PSSW * J(-) CONTRIBUTION
*----
      DO 90 IKK=1,NSUPCE
      DO 90 JKK=1,NSUPCE
      PSS=PROCEL(IKK,JKK)*PSSW(JKK)
   90 CITR(IT3,IKK)=CITR(IT3,IKK)+PSS*CITR(IT2,JKK)
*----
*  NORMALIZATION
*----
      S1=0.0
      S2=0.0
      DO 100 I=1,NSUPCE
      S1=S1+WCURR(I)*CITR(IT3,I)
  100 S2=S2+CIT0(I)
      ZNORM=REAL(S2/S1)
      IF(ZNORM.LT.0.0) ZNORM=1.0
      DO 110 I=1,NSUPCE
  110 CITR(IT3,I)=CITR(IT3,I)*ZNORM
*----
*  ONE/TWO PARAMETER ACCELERATION
*----
      ALP=1.0
      BET=0.0
      LOGTES=(1+MOD(ITER-ISTART,ICL1+ICL2).GT.ICL1)
      IF(LOGTES) THEN
         DO 130 IKK=1,NSUPCE
         AITR(1,IKK)=CITR(IT3,IKK)-CITR(IT2,IKK)
         AITR(2,IKK)=CITR(IT2,IKK)-CITR(IT1,IKK)
         DO 130 JKK=1,NSUPCE
         PSS=PROCEL(IKK,JKK)*PSSW(JKK)
         AITR(1,IKK)=AITR(1,IKK)-PSS*(CITR(IT3,JKK)-CITR(IT2,JKK))
  130    AITR(2,IKK)=AITR(2,IKK)-PSS*(CITR(IT2,JKK)-CITR(IT1,JKK))
         IF((LACCFC.EQ.1).OR.(MOD(ITER-ISTART,ICL1+ICL2).EQ.ICL1)) THEN
            S1=0.0
            S2=0.0
            DO 140 I=1,NSUPCE
            S1=S1+(CITR(IT3,I)-CITR(IT2,I))*AITR(1,I)
  140       S2=S2+AITR(1,I)*AITR(1,I)
            IF(S2.EQ.0.0) THEN
               ISTART=ITER+1
            ELSE
               ALP=REAL(S1/S2)
               IF(ALP.LE.0.0) THEN
                  ISTART=ITER+1
                  ALP=1.0
               ENDIF
            ENDIF
            DO 150 I=1,NSUPCE
  150       CITR(IT3,I)=CITR(IT2,I)+ALP*(CITR(IT3,I)-CITR(IT2,I))
         ELSE IF(LACCFC.EQ.2) THEN
            S1=0.0
            S2=0.0
            S3=0.0
            S4=0.0
            S5=0.0
            DO 160 I=1,NSUPCE
            S1=S1+(CITR(IT3,I)-CITR(IT2,I))*AITR(1,I)
            S2=S2+AITR(1,I)*AITR(1,I)
            S3=S3+(CITR(IT3,I)-CITR(IT2,I))*AITR(2,I)
            S4=S4+AITR(1,I)*AITR(2,I)
  160       S5=S5+AITR(2,I)*AITR(2,I)
            DET=S2*S5-S4*S4
            IF(DET.EQ.0.0) THEN
               ISTART=ITER+1
            ELSE
               ALP=REAL((S5*S1-S4*S3)/DET)
               BET=REAL((S2*S3-S4*S1)/DET)
               IF(ALP.LE.0.0) THEN
                  ISTART=ITER+1
                  ALP=1.0
                  BET=0.0
               ENDIF
            ENDIF
            DO 170 I=1,NSUPCE
  170       CITR(IT3,I)=CITR(IT2,I)+ALP*(CITR(IT3,I)-CITR(IT2,I))+
     1      BET*(CITR(IT2,I)-CITR(IT1,I))
         ENDIF
      ENDIF
*----
*  CHECK THE CONVERGENCE ERROR
*----
      ERR1=0.0
      ERR2=0.0
      DO 180 I=1,NSUPCE
      ERR1=MAX(ERR1,ABS(CITR(IT3,I)-CITR(IT2,I)))
  180 ERR2=MAX(ERR2,ABS(CITR(IT3,I)))
      IF(IMPX.GT.3) WRITE(6,'(30H SYBJJ0: CURRENT ITERATION NB.,I4,
     1 7H ERROR=,1P,E10.3,5H OVER,E10.3,15H NORMALIZATION=,E10.3,
     2 14H ACCELERATION=,2E11.3,1H.)') ITER,ERR1,ERR2,ZNORM,ALP,
     3 BET/ALP
      IF(ITER.EQ.1) TEST=ERR1/ERR2
      IF((ITER.GT.20).AND.(ERR1/ERR2.GT.TEST)) CALL XABORT('SYBJJ0: '
     1 //'CONVERGENCE FAILURE.')
      IF(LOGTES.OR.(ERR1.GT.EPSJ*ERR2)) GO TO 70
      IF(IMPX.GT.2) WRITE(6,'(37H SYBJJ0: CURRENT CONVERGENCE AT ITERA,
     1 8HTION NB.,I4,7H ERROR=,1P,E10.3,5H OVER,E10.3,1H.)') ITER,ERR1,
     2 ERR2
*
  190 DO 200 I=1,IPAS
  200 FUNKNO(I)=0.0
      DO 210 I=1,NSUPCE
  210 FUNKNO(IPAS+I)=CITR(IT3,I)
*----
*  COMPUTE ( PISW * J(-) ) + ( PIJW * Q(*) ) CONTRIBUTION
*----
      DO 220 IKK=1,NSUPCE
      I1P=NMC(IKK)
      I2=NMC(IKK+1)-I1P
      DO 220 J=1,I2
      FUNKNO(I1P+J)=FUNKNO(I1P+J)+PISW(I1P+J)*FUNKNO(IPAS+IKK)
      DO 220 I=1,I2
      PIJ=PIJW(INDPIJ(IKK)+(I-1)*I2+J)
  220 FUNKNO(I1P+J)=FUNKNO(I1P+J)+PIJ*SUNKNO(I1P+I)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(WCURR)
      DEALLOCATE(AITR,CIT0,CITR)
      DEALLOCATE(INDPIJ)
      RETURN
      END
