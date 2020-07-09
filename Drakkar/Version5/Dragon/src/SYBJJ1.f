*DECK SYBJJ1
      SUBROUTINE SYBJJ1 (IPAS,NMCEL,NMERGE,NGEN,NPIJ,NPIS,EPSJ,NUNKNO,
     1 FUNKNO,SUNKNO,IMPX,NCOUR,XX,YY,NMC,IFR,ALB,INUM,IGEN,PIJW,PISW,
     2 PSJW,PSSW)
*
*-----------------------------------------------------------------------
*
*Purpose:
* compute the neutron flux and interface currents in a 2-D Cartesian
* or hexagonal assembly using the current iteration method with Roth
* approximation.
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
* NMCEL   total number of cells in the domain.
* NMERGE  total number of merged cells for which specific values
*         of the neutron flux and reactions rates are required.
*         Many cells with different position in the domain can
*         be merged before the neutron flux calculation if they
*         own the same generating cell. Equal to the total number
*         of distinct out-currents (NMERGE.le.NMCEL).
* NGEN    total number of generating cells. A generating cell is
*         defined by its material and dimensions, irrespective of
*         its position in the domain (NGEN.le.NMERGE).
* NPIJ    size of cellwise scattering-reduced collision probability
*         matrices.
* NPIS    size of cellwise scattering-reduced escape probability
*         matrices.
* EPSJ    stopping criterion for flux-current iterations.
* NUNKNO  total number of unknowns in vectors SUNKNO and FUNKNO.
* SUNKNO  input source vector.
* IMPX    print flag (equal to 0 for no print).
* NCOUR   number of incoming currents (=4: Cartesian lattice;
*         =6: hexagonal lattice).
* XX      X-thickness of the generating cells.
* YY      Y-thickness of the generating cells.
* NMC     offset of the first volume in each generating cell.
* IFR     index-number of in-currents.
* ALB     transmission/albedo associated with each in-current.
*         Note: IFR and ALB contains information to rebuild the
*         geometrical 'A' matrix.
* INUM    index-number of the merged cell associated to each cell.
* IGEN    index-number of the generating cell associated with each
*         merged cell.
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
      INTEGER IPAS,NMCEL,NMERGE,NGEN,NPIJ,NPIS,NUNKNO,IMPX,NCOUR,
     1 NMC(NGEN+1),IFR(NCOUR*NMCEL),INUM(NMCEL),IGEN(NMERGE)
      REAL EPSJ,FUNKNO(NUNKNO),SUNKNO(NUNKNO),XX(NGEN),YY(NGEN),
     1 ALB(NCOUR*NMCEL),PIJW(NPIJ),PISW(NPIS),PSJW(NPIS),PSSW(NGEN)
*----
*  LOCAL VARIABLES
*----
      REAL PIBB(6)
      DOUBLE PRECISION S1,S2,S3,S4,S5,DET
      LOGICAL LOGTES
      PARAMETER (MAXIT=400,LACCFC=2,ICL1=3,ICL2=3)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, DIMENSION(:), POINTER :: INDPIJ,INDNMC
      REAL, DIMENSION(:), POINTER :: CIT0
      REAL, DIMENSION(:,:), POINTER :: CITR,AITR
      DOUBLE PRECISION, DIMENSION(:), POINTER :: WCURR
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(INDPIJ(NGEN),INDNMC(NMERGE))
      ALLOCATE(CITR(3,NMERGE),CIT0(NMERGE),AITR(2,NMERGE))
      ALLOCATE(WCURR(NMERGE))
*
      IPIJ=0
      DO 10 JKG=1,NGEN
      J2=NMC(JKG+1)-NMC(JKG)
      INDPIJ(JKG)=IPIJ
   10 IPIJ=IPIJ+J2*J2
      KNMC=0
      DO 20 JKK=1,NMERGE
      JKG=IGEN(JKK)
      J2=NMC(JKG+1)-NMC(JKG)
      INDNMC(JKK)=KNMC
   20 KNMC=KNMC+J2
*
      DO 30 I=1,NMERGE
      WCURR(I)=1.0D0
      CIT0(I)=0.0
   30 CITR(1,I)=FUNKNO(IPAS+I)
*----
*  COMPUTE PSJW * Q(*) CONTRIBUTION
*----
      DO 40 IKK=1,NMERGE
      IKG=IGEN(IKK)
      I1P=NMC(IKG)
      I2=NMC(IKG+1)-I1P
      KNMC=INDNMC(IKK)
      DO 40 I=1,I2
   40 CIT0(IKK)=CIT0(IKK)+PSJW(I1P+I)*SUNKNO(KNMC+I)
*----
*  COMPUTE NORMALIZATION VECTOR WCURR
*----
      DO 60 ICEL=1,NMCEL
      IKK=INUM(ICEL)
      IS=NCOUR*(ICEL-1)
      IKG=IGEN(IKK)
      IF(NCOUR.EQ.4) THEN
         A=XX(IKG)
         B=YY(IKG)
         DEN1=2.0*(A+B)
         PIBB(1)=B/DEN1
         PIBB(2)=B/DEN1
         PIBB(3)=A/DEN1
         PIBB(4)=A/DEN1
      ELSE
         DO 50 JC=1,NCOUR
   50    PIBB(JC)=1.0/6.0
      ENDIF
      DO 60 JC=1,NCOUR
      J1=IFR(IS+JC)
   60 WCURR(J1)=WCURR(J1)-PSSW(IKG)*PIBB(JC)*ALB(IS+JC)
*
      ISTART=1
      TEST=0.0
      ITER=0
   70 ITER=ITER+1
      IF(ITER.GT.MAXIT) THEN
         WRITE(6,'(/47H SYBJJ1: *** WARNING *** MAXIMUM NUMBER OF ITER,
     1   15HATIONS REACHED.)')
         GO TO 190
      ENDIF
      IT3=MOD(ITER,3)+1
      IT2=MOD(ITER-1,3)+1
      IT1=MOD(ITER-2,3)+1
      DO 80 I=1,NMERGE
   80 CITR(IT3,I)=CIT0(I)
*----
*  COMPUTE PSSW * J(-) CONTRIBUTION
*----
      DO 90 ICEL=1,NMCEL
      IKK=INUM(ICEL)
      IS=NCOUR*(ICEL-1)
      IKG=IGEN(IKK)
      IF(NCOUR.EQ.4) THEN
         A=XX(IKG)
         B=YY(IKG)
         DEN1=2.0*(A+B)
         PIBB(1)=B/DEN1
         PIBB(2)=B/DEN1
         PIBB(3)=A/DEN1
         PIBB(4)=A/DEN1
      ELSE
         DO 85 JC=1,NCOUR
   85    PIBB(JC)=1.0/6.0
      ENDIF
      DO 90 JC=1,NCOUR
      J1=IFR(IS+JC)
      PSS=PSSW(IKG)*PIBB(JC)
   90 CITR(IT3,IKK)=CITR(IT3,IKK)+PSS*ALB(IS+JC)*CITR(IT2,J1)
*----
*  NORMALIZATION
*----
      S1=0.0
      S2=0.0
      DO 100 I=1,NMERGE
      S1=S1+WCURR(I)*CITR(IT3,I)
  100 S2=S2+CIT0(I)
      ZNORM=REAL(S2/S1)
      IF(ZNORM.LT.0.0) ZNORM=1.0
      DO 110 I=1,NMERGE
  110 CITR(IT3,I)=CITR(IT3,I)*ZNORM
*----
*  ONE/TWO PARAMETER ACCELERATION
*----
      ALP=1.0
      BET=0.0
      LOGTES=(1+MOD(ITER-ISTART,ICL1+ICL2).GT.ICL1)
      IF(LOGTES) THEN
         DO 120 I=1,NMERGE
         AITR(1,I)=CITR(IT3,I)-CITR(IT2,I)
  120    AITR(2,I)=CITR(IT2,I)-CITR(IT1,I)
         DO 130 ICEL=1,NMCEL
         IKK=INUM(ICEL)
         IS=NCOUR*(ICEL-1)
         IKG=IGEN(IKK)
         IF(NCOUR.EQ.4) THEN
            A=XX(IKG)
            B=YY(IKG)
            DEN1=2.0*(A+B)
            PIBB(1)=B/DEN1
            PIBB(2)=B/DEN1
            PIBB(3)=A/DEN1
            PIBB(4)=A/DEN1
         ELSE
            DO 125 JC=1,NCOUR
  125       PIBB(JC)=1.0/6.0
         ENDIF
         DO 130 JC=1,NCOUR
         J1=IFR(IS+JC)
         PSS=PSSW(IKG)*PIBB(JC)*ALB(IS+JC)
         AITR(1,IKK)=AITR(1,IKK)-PSS*(CITR(IT3,J1)-CITR(IT2,J1))
  130    AITR(2,IKK)=AITR(2,IKK)-PSS*(CITR(IT2,J1)-CITR(IT1,J1))
         IF((LACCFC.EQ.1).OR.(MOD(ITER-ISTART,ICL1+ICL2).EQ.ICL1)) THEN
            S1=0.0
            S2=0.0
            DO 140 I=1,NMERGE
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
            DO 150 I=1,NMERGE
  150       CITR(IT3,I)=CITR(IT2,I)+ALP*(CITR(IT3,I)-CITR(IT2,I))
         ELSE IF(LACCFC.EQ.2) THEN
            S1=0.0
            S2=0.0
            S3=0.0
            S4=0.0
            S5=0.0
            DO 160 I=1,NMERGE
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
            DO 170 I=1,NMERGE
  170       CITR(IT3,I)=CITR(IT2,I)+ALP*(CITR(IT3,I)-CITR(IT2,I))+
     1      BET*(CITR(IT2,I)-CITR(IT1,I))
         ENDIF
      ENDIF
*----
*  CHECK THE CONVERGENCE ERROR
*----
      ERR1=0.0
      ERR2=0.0
      DO 180 I=1,NMERGE
      ERR1=MAX(ERR1,ABS(CITR(IT3,I)-CITR(IT2,I)))
  180 ERR2=MAX(ERR2,ABS(CITR(IT3,I)))
      IF(IMPX.GT.3) WRITE(6,'(30H SYBJJ1: CURRENT ITERATION NB.,I4,
     1 7H ERROR=,1P,E10.3,5H OVER,E10.3,15H NORMALIZATION=,E10.3,
     2 14H ACCELERATION=,2E11.3,1H.)') ITER,ERR1,ERR2,ZNORM,ALP,
     3 BET/ALP
      IF(ITER.EQ.1) TEST=ERR1/ERR2
      IF((ITER.GT.20).AND.(ERR1/ERR2.GT.TEST)) CALL XABORT('SYBJJ1: '
     1 //'CONVERGENCE FAILURE.')
      IF(LOGTES.OR.(ERR1.GT.EPSJ*ERR2)) GO TO 70
      IF(IMPX.GT.2) WRITE(6,'(37H SYBJJ1: CURRENT CONVERGENCE AT ITERA,
     1 8HTION NB.,I4,7H ERROR=,1P,E10.3,5H OVER,E10.3,1H.)') ITER,ERR1,
     2 ERR2
*
  190 DO 200 I=1,IPAS
  200 FUNKNO(I)=0.0
      DO 210 I=1,NMERGE
  210 FUNKNO(IPAS+I)=CITR(IT3,I)
*----
*  COMPUTE PISW * J(-) CONTRIBUTION
*----
      DO 220 ICEL=1,NMCEL
      IKK=INUM(ICEL)
      IS=NCOUR*(ICEL-1)
      IKG=IGEN(IKK)
      IF(NCOUR.EQ.4) THEN
         A=XX(IKG)
         B=YY(IKG)
         DEN1=2.0*(A+B)
         PIBB(1)=B/DEN1
         PIBB(2)=B/DEN1
         PIBB(3)=A/DEN1
         PIBB(4)=A/DEN1
      ELSE
         DO 215 JC=1,NCOUR
  215    PIBB(JC)=1.0/6.0
      ENDIF
      I1P=NMC(IKG)
      I2=NMC(IKG+1)-I1P
      KNMC=INDNMC(IKK)
      DO 220 J=1,I2
      DO 220 JC=1,NCOUR
      J1=IFR(IS+JC)
      PIS=PISW(I1P+J)*PIBB(JC)
  220 FUNKNO(KNMC+J)=FUNKNO(KNMC+J)+PIS*ALB(IS+JC)*FUNKNO(IPAS+J1)
*----
*  COMPUTE PIJW * Q(*) CONTRIBUTION
*----
      DO 230 IKK=1,NMERGE
      IKG=IGEN(IKK)
      I2=NMC(IKG+1)-NMC(IKG)
      KNMC=INDNMC(IKK)
      DO 230 I=1,I2
      DO 230 J=1,I2
      PIJ=PIJW(INDPIJ(IKG)+(I-1)*I2+J)
  230 FUNKNO(KNMC+J)=FUNKNO(KNMC+J)+PIJ*SUNKNO(KNMC+I)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(WCURR)
      DEALLOCATE(AITR,CIT0,CITR)
      DEALLOCATE(INDNMC,INDPIJ)
      RETURN
      END
