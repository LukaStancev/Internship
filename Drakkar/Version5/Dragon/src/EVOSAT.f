*DECK EVOSAT
      SUBROUTINE EVOSAT(IMPX,MAXA,MAXB,MAXY,LOGY,NSAT,NVAR,KSAT,YST1,
     1 YSAT,MU1,IMA,NSUPF,NFISS,IDIRAC,KFISS,YSF,ADPL,BDPL,NSUPFG)
*
*-----------------------------------------------------------------------
*
*Purpose:
* lumping of the depletion matrix, fission yields, sources and initial
* conditions to take into account the saturation of depleting nuclides.
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
* IMPX    print parameter.
* MAXA    first dimension of matrices ADPL and AGAR.
* MAXB    first dimension of matrices BDPL, IMA and MU1.
* MAXY    second dimension of matrix YSF.
* LOGY    number of passes through EVOSAT:
*         first pass: update YSAT and YST1;
*         second pass: do not update YSAT and YST1.
* NSAT    number of saturating nuclides.
* NVAR    number of nuclides in the complete depletion chain.
* KSAT    position in chain of the saturating nuclides.
* YST1    number densities for all isotopes.
* NFISS   number of fissile isotopes producing fission products.
* IDIRAC  saturation model flag (=1 to use Dirac function contributions
*         in the saturating nuclide number densities.
* MU1     position of each diagonal element in vector ADPL.
* IMA     position of the first non-zero column element in vector ADPL.
* NSUPF   number of depleting fission products.
* KFISS   position in chain of the fissile isotopes.
* YSF     product of the fission yields and fission rates.
* ADPL    depletion matrix.
* BDPL    depletion source.
*
*Parameters: output
* NSUPFG  number of lumped depleting fission products.
* YST1    number densities of the non-saturated isotopes.
* YSAT    number densities of the saturating isotopes.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IMPX,MAXA,MAXB,MAXY,LOGY,NSAT,NVAR,KSAT(NSAT),MU1(MAXB),
     1 IMA(MAXB),NSUPF,NFISS,IDIRAC,KFISS(NFISS),NSUPFG
      REAL YST1(NVAR),YSAT(NSAT),YSF(NFISS,MAXY,LOGY),ADPL(MAXA,LOGY),
     1 BDPL(MAXB,LOGY)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(EPS=1.0E-5)
      CHARACTER HSMG*131
      LOGICAL LTEST
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: KEV,MGAR,IGAR
      REAL, ALLOCATABLE, DIMENSION(:) :: YSTG,AGAR,BGAR,GAR
      REAL, ALLOCATABLE, DIMENSION(:,:) :: A22,YSFG
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: A21,A12
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(KEV(NVAR),MGAR(NVAR-NSAT),IGAR(NVAR-NSAT))
      ALLOCATE(YSTG(NVAR-NSAT),A22(NSAT,NSAT),A21(NSAT,NVAR-NSAT,LOGY),
     1 A12(NVAR-NSAT,NSAT,LOGY),AGAR(MAXA),BGAR(NVAR-NSAT),
     2 YSFG(NFISS,NSUPF),GAR(NSAT))
*
      NSUPL=NVAR-NSUPF
      I0=0
      DO 40 I=1,NVAR
      DO 10 II=1,NSAT
      IF(I.EQ.KSAT(II)) GO TO 20
   10 CONTINUE
      I0=I0+1
      KEV(I)=I0
      GO TO 40
   20 DO 25 L=1,LOGY
      IF(ADPL(MU1(I),L).EQ.0.0) CALL XABORT('EVOSAT: ZERO DIAGONAL COM'
     1 //'PONENT FOR A SATURATING ISOTOPE.')
   25 CONTINUE
      DO 30 II=1,NFISS
      IF(I.EQ.KFISS(II)) CALL XABORT('EVOSAT: A FISSILE ISOTOPE IS SAT'
     1 //'URATING.')
   30 CONTINUE
      KEV(I)=0
   40 CONTINUE
      DO 50 I=1,NFISS
   50 KFISS(I)=KEV(KFISS(I))
*----
*  FIRST LOOP OVER LOGY
*----
      DO 270 L=1,LOGY
*----
*  COMPUTE MATRICES A22**-1, A21, AND A12
*----
      DO 90 II=1,NSAT
      I=KSAT(II)
      IMAM1=0
      IF(I.GT.1) IMAM1=IMA(I-1)
      DO 60 JJ=1,NSAT
      J=KSAT(JJ)
      IF((J.LE.I).AND.(J.GT.I+IMAM1-MU1(I))) THEN
         A22(II,JJ)=ADPL(MU1(I)-I+J,L)
      ELSE IF((I.LE.J).AND.(I.GE.J-IMA(J)+MU1(J))) THEN
         A22(II,JJ)=ADPL(MU1(J)+J-I,L)
      ELSE
         A22(II,JJ)=0.0
      ENDIF
   60 CONTINUE
      JMAM1=0
      DO 70 J=1,NVAR
      J0=KEV(J)
      IF(J0.EQ.0) GO TO 70
      IF((J.LE.I).AND.(J.GT.I+IMAM1-MU1(I))) THEN
         A21(II,J0,L)=ADPL(MU1(I)-I+J,L)
      ELSE IF((I.LE.J).AND.(I.GE.J-IMA(J)+MU1(J))) THEN
         A21(II,J0,L)=ADPL(MU1(J)+J-I,L)
      ELSE
         A21(II,J0,L)=0.0
      ENDIF
      IF((I.LE.J).AND.(I.GT.J+JMAM1-MU1(J))) THEN
         A12(J0,II,L)=ADPL(MU1(J)-J+I,L)
      ELSE IF((J.LE.I).AND.(J.GE.I-IMA(I)+MU1(I))) THEN
         A12(J0,II,L)=ADPL(MU1(I)+I-J,L)
      ELSE
         A12(J0,II,L)=0.0
      ENDIF
   70 JMAM1=IMA(J)
      IF(I.GT.NSUPL) THEN
         DO 80 K=1,NFISS
   80    A21(II,KFISS(K),L)=A21(II,KFISS(K),L)+YSF(K,I-NSUPL,L)
      ENDIF
   90 CONTINUE
      CALL ALINV(NSAT,A22,NSAT,IER)
      IF(IER.NE.0) CALL XABORT('EVOSAT: SINGULAR MATRIX.')
*----
*  COMPUTE VECTOR YSTG ANT YSAT
*----
      IF(L.EQ.1) THEN
*        BEGINNING-OF-STAGE DIRAC DELTA CONTRIBUTIONS:
         DO 100 I=1,NSAT
  100    YSAT(I)=YST1(KSAT(I))
         DO 110 I=1,NVAR
         IF(KEV(I).GT.0) YSTG(KEV(I))=YST1(I)
  110    CONTINUE
         IF(IDIRAC.EQ.0) THEN
            DO 120 I=1,NSAT
            GAR(I)=BDPL(KSAT(I),L)
            DO 120 J=1,NVAR-NSAT
  120       GAR(I)=GAR(I)+A21(I,J,L)*YSTG(J)
            DO 130 I=1,NSAT
            YSAT(I)=0.0
            DO 130 J=1,NSAT
  130       YSAT(I)=YSAT(I)-A22(I,J)*GAR(J)
            GO TO 220
         ENDIF
         ITER=0
  140    ITER=ITER+1
         IF(ITER.GT.50) CALL XABORT('EVOSAT: CONVERGENCE FAILURE.')
         DO 150 I=1,NSAT
         GAR(I)=BDPL(KSAT(I),L)
         DO 150 J=1,NVAR-NSAT
  150    GAR(I)=GAR(I)+A21(I,J,L)*YSTG(J)
         ERR1=0.0
         ERR2=0.0
         DO 170 I=1,NSAT
         ZCOMP=YSAT(I)
         YSAT(I)=0.0
         DO 160 J=1,NSAT
  160    YSAT(I)=YSAT(I)-A22(I,J)*GAR(J)
         ERR1=MAX(ERR1,ABS(ZCOMP-YSAT(I)))
  170    ERR2=MAX(ERR2,ABS(YSAT(I)))
         DO 180 I=1,NSAT
         GAR(I)=0.0
         DO 180 J=1,NSAT
  180    GAR(I)=GAR(I)-A22(I,J)*(YST1(KSAT(J))-YSAT(J))
         DO 190 I=1,NVAR
         IF(KEV(I).GT.0) YSTG(KEV(I))=YST1(I)
  190    CONTINUE
         DO 210 I=1,NVAR-NSAT
         DO 200 J=1,NSAT
  200    YSTG(I)=YSTG(I)+A12(I,J,L)*GAR(J)
  210    ERR2=MAX(ERR2,ABS(YSTG(I)))
         IF(ERR1.LE.EPS*ERR2) GO TO 220
         GO TO 140
      ENDIF
*----
*  COMPUTE MATRICES A21 AND BGAR
*----
  220 DO 230 I=1,NSAT
      GAR(I)=0.0
      DO 230 J=1,NSAT
  230 GAR(I)=GAR(I)-A22(I,J)*BDPL(KSAT(J),L)
      CALL XDRSET(BGAR,NVAR-NSAT,0.0)
      DO 240 I=1,NVAR
      IF(KEV(I).GT.0) BGAR(KEV(I))=BDPL(I,L)
  240 CONTINUE
      DO 250 I=1,NVAR-NSAT
      DO 250 J=1,NSAT
  250 BGAR(I)=BGAR(I)+A12(I,J,L)*GAR(J)
      DO 270 J=1,NVAR-NSAT
      BDPL(J,L)=BGAR(J)
      IF(L.EQ.1) YST1(J)=YSTG(J)
      DO 260 K=1,NSAT
  260 GAR(K)=A21(K,J,L)
      DO 270 I=1,NSAT
      A21(I,J,L)=0.0
      DO 270 K=1,NSAT
  270 A21(I,J,L)=A21(I,J,L)+A22(I,K)*GAR(K)
*----
*  DETERMINE THE PROFILE PATTERN OF THE LUMPED DEPLETION MATRIX.
*----
      NSUPLG=NSUPL
      DO 280 I=1,NVAR
      IF((KEV(I).EQ.0).AND.(I.LE.NSUPL)) NSUPLG=NSUPLG-1
  280 CONTINUE
      NSUPFG=NVAR-NSAT-NSUPLG
      CALL XDISET(MGAR,NVAR-NSAT,1)
      CALL XDISET(IGAR,NVAR-NSAT,1)
      IMAM1=0
      DO 300 I=1,NVAR
      IKEV=KEV(I)
      IF(IKEV.EQ.0) GO TO 300
      DO 290 J=1,NVAR
      JKEV=KEV(J)
      IF(JKEV.EQ.0) GO TO 290
      IF((J.LE.I).AND.(J.GT.I+IMAM1-MU1(I))) THEN
         MGAR(IKEV)=MAX(MGAR(IKEV),IKEV-JKEV+1)
      ELSE IF((I.LE.J).AND.(I.GE.J-IMA(J)+MU1(J))) THEN
         IGAR(JKEV)=MAX(IGAR(JKEV),JKEV-IKEV+1)
      ENDIF
  290 CONTINUE
  300 IMAM1=IMA(I)
      DO 330 J=1,NVAR-NSAT
      JIFI=0
      DO 310 IFI=1,NFISS
      IF(J.EQ.KFISS(IFI)) JIFI=IFI
  310 CONTINUE
      DO 330 I=1,NVAR-NSAT
      IF((I.GT.NSUPLG).AND.(JIFI.GT.0)) GO TO 330
      LTEST=.FALSE.
      DO 320 L=1,LOGY
      DO 320 K=1,NSAT
  320 LTEST=LTEST.OR.(A12(I,K,L)*A21(K,J,L).NE.0.0)
      IF(LTEST.AND.(J.LE.I)) THEN
         MGAR(I)=MAX(MGAR(I),I-J+1)
      ELSE IF(LTEST) THEN
         IGAR(J)=MAX(IGAR(J),J-I+1)
      ENDIF
  330 CONTINUE
      II=0
      DO 340 I=1,NVAR-NSAT
      II=II+MGAR(I)
      MGAR(I)=II
      II=II+IGAR(I)-1
  340 IGAR(I)=II
      IF(IMPX.GT.8) WRITE(6,'(/27H EVOSAT: REAL SIZE OF ADPL=,I9,3H AL,
     1 13HLOCATED SIZE=,I9,1H.)') IGAR(NVAR-NSAT),MAXA
      IF(IGAR(NVAR-NSAT).GT.MAXA) THEN
         WRITE(HSMG,'(24HEVOSAT: IGAR(NVAR-NSAT)=,I6,6H MAXA=,I6)')
     1   IGAR(NVAR-NSAT),MAXA
         CALL XABORT(HSMG)
      ENDIF
*----
*  SECOND LOOP OVER LOGY
*----
      DO 540 L=1,LOGY
*----
*  COMPUTE MATRIX AGAR AND YIELDS YSFG.
*----
      CALL XDRSET(AGAR,IGAR(NVAR-NSAT),0.0)
      IMAM1=0
      DO 440 I=1,NVAR
      IKEV=KEV(I)
      IF(IKEV.EQ.0) GO TO 440
      DO 420 J=1,NVAR
      JKEV=KEV(J)
      IF(JKEV.EQ.0) GO TO 420
      IF((J.LE.I).AND.(J.GT.I+IMAM1-MU1(I))) THEN
         AGAR(MGAR(IKEV)-IKEV+JKEV)=ADPL(MU1(I)-I+J,L)
      ELSE IF((I.LE.J).AND.(I.GE.J-IMA(J)+MU1(J))) THEN
         AGAR(MGAR(JKEV)+JKEV-IKEV)=ADPL(MU1(J)+J-I,L)
      ENDIF
  420 CONTINUE
      IF(I.GT.NSUPL) THEN
         DO 430 K=1,NFISS
  430    YSFG(K,IKEV-NSUPLG)=YSF(K,I-NSUPL,L)
      ENDIF
  440 IMAM1=IMA(I)
      DO 480 J=1,NVAR-NSAT
      JIFI=0
      DO 450 IFI=1,NFISS
      IF(J.EQ.KFISS(IFI)) JIFI=IFI
  450 CONTINUE
      IMAM1=0
      DO 480 I=1,NVAR-NSAT
      IF((I.GT.NSUPLG).AND.(JIFI.GT.0)) GO TO 480
      IF((J.LE.I).AND.(J.GT.I+IMAM1-MGAR(I))) THEN
         DO 460 K=1,NSAT
  460    AGAR(MGAR(I)-I+J)=AGAR(MGAR(I)-I+J)-A12(I,K,L)*A21(K,J,L)
      ELSE IF((I.LE.J).AND.(I.GE.J-IGAR(J)+MGAR(J))) THEN
         DO 470 K=1,NSAT
  470    AGAR(MGAR(J)+J-I)=AGAR(MGAR(J)+J-I)-A12(I,K,L)*A21(K,J,L)
      ENDIF
  480 IMAM1=IGAR(I)
      DO 490 I=NSUPLG+1,NVAR-NSAT
      DO 490 IFI=1,NFISS
      J=KFISS(IFI)
      DO 490 K=1,NSAT
  490 YSFG(IFI,I-NSUPLG)=YSFG(IFI,I-NSUPLG)-A12(I,K,L)*A21(K,J,L)
*----
*  REPLACE THE ORIGINAL INFORMATION WITH THE LUMPED ONE
*----
      DO 520 I=1,IGAR(NVAR-NSAT)
  520 ADPL(I,L)=AGAR(I)
      DO 530 I=1,NFISS
      DO 530 J=1,NSUPFG
  530 YSF(I,J,L)=YSFG(I,J)
  540 CONTINUE
      DO 550 I=1,NVAR-NSAT
      IMA(I)=IGAR(I)
  550 MU1(I)=MGAR(I)
      RETURN
      END
