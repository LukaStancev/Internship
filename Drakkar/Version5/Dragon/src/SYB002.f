*DECK SYB002
      SUBROUTINE SYB002 (NGEN,NPIJ,NPIS,SIGT2,SIGW2,IMPX,NCOUR,IWIGN,
     1 IQUAD,XX,YY,NMC,RAYRE,MAIL,RZMAIL,PIJW,PISW,PSJW,PSSW)
*
*-----------------------------------------------------------------------
*
*Purpose:
* compute the cellwise scattering-reduced collision, escape and
* transmission probabilities in a 2-D Cartesian or hexagonal assembly
* with Roth approximation.
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
* NGEN    total number of generating cells.
* NPIJ    length of cellwise scattering-reduced collision probability
*         matrices.
* NPIS    length of cellwise scattering-reduced escape probability
*         matrices (NPIS=NMC(NGEN+1)).
* SIGT2   total macroscopic cross sections.
* SIGW2   P0 within-group scattering macroscopic cross sections.
* IMPX    print flag (equal to 0 for no print).
* NCOUR   number of currents surrounding the cells (=4: Cartesian
*         lattice; =6: hexagonal lattice).
* IWIGN   type of cylinderization.
* IQUAD   quadrature parameters.
* XX      X-thickness of the generating cells.
* YY      Y-thickness of the generating cells.
* NMC     offset of the first volume in each generating cell.
* RAYRE   radius of the tubes in each generating cell.
* MAIL    offset of the first tracking information in each generating
*         cell.
* RZMAIL  real tracking information.
*
*Parameters: output
* PIJW    cellwise scattering-reduced collision probability matrices.
* PISW    cellwise scattering-reduced escape probability matrices.
* PSJW    cellwise scattering-reduced collision probability matrices
*         for incoming neutrons.
* PSSW    cellwise scattering-reduced transmission probability
*         matrices.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NGEN,NPIJ,NPIS,IMPX,NCOUR,IWIGN,IQUAD(4),NMC(NGEN+1),
     1 MAIL(2,NGEN)
      REAL SIGT2(NPIS),SIGW2(NPIS),XX(NGEN),YY(NGEN),RAYRE(NPIS),
     1 RZMAIL(*),PIJW(NPIJ),PISW(NPIS),PSJW(NPIS),PSSW(NGEN)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (PI=3.141592654)
      LOGICAL LSKIP
      REAL, ALLOCATABLE, DIMENSION(:) :: PIS,PSJ,RAYR2,WORK
      REAL, ALLOCATABLE, DIMENSION(:,:) :: PP
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(PIS(NPIS),PSJ(NPIS))
*
      MR=IQUAD(4)
      IPIJ=0
      DO 220 JKG=1,NGEN
      J1=NMC(JKG)
      J2=NMC(JKG+1)-J1
*----
*  CYLINDERIZATION OPTIONS
*----
      A=XX(JKG)
      B=YY(JKG)
      IB=MAIL(2,JKG)
      RJ1=RAYRE(NMC(JKG+1))
      SCALE1=1.0
      SCALE2=1.0
      ROUT=0.0
      IF((NCOUR.EQ.4).AND.(IWIGN.EQ.1)) THEN
*        ASKEW CYLINDERIZATION CARTESIAN.
         RJ2=(A+B)/PI
         SCALE1=(A*B-PI*RJ1**2)/(PI*RJ2**2-PI*RJ1**2)
         ROUT=RJ2
      ELSE IF((NCOUR.EQ.4).AND.(IWIGN.EQ.2)) THEN
*        WIGNER CYLINDERIZATION CARTESIAN.
         ROUT=SQRT(A*B/PI)
      ELSE IF((NCOUR.EQ.4).AND.(IWIGN.EQ.3)) THEN
*        SANCHEZ CYLINDERIZATION CARTESIAN.
         SCALE2=SQRT(PI*A*B)/(A+B)
         ROUT=SQRT(A*B/PI)
      ELSE IF(IWIGN.EQ.1) THEN
*        ASKEW CYLINDERIZATION HEXAGONAL.
         RJ2=3.0*A/PI
         SCALE1=(1.5*SQRT(3.0)*A*A-PI*RJ1**2)/(PI*RJ2**2-PI*RJ1**2)
         ROUT=RJ2
      ELSE IF(IWIGN.EQ.2) THEN
*        WIGNER CYLINDERIZATION HEXAGONAL.
         ROUT=SQRT(1.5*SQRT(3.0)/PI)*A
      ELSE IF(IWIGN.EQ.3) THEN
*        SANCHEZ CYLINDERIZATION HEXAGONAL.
         SCALE2=SQRT(PI*SQRT(3.0)/6.0)
         ROUT=SQRT(1.5*SQRT(3.0)/PI)*A
      ENDIF
      IF(ROUT.LE.RJ1) CALL XABORT('SYB002: CYLINDERIZATION ERROR.')
*----
*  COMPUTE THE REDUCED COLLISION PROBABILITY MATRIX
*----
      SURFA=0.5*PI*ROUT
      ALLOCATE(PP(J2,J2),RAYR2(J2+1))
      DO 10 I=1,J2
   10 RAYR2(I)=RAYRE(J1+I)
      RAYR2(J2+1)=ROUT
      SIGT2(J1+J2)=SIGT2(J1+J2)*SCALE1
      CALL SYBALC(J2,J2,RAYR2,SIGT2(J1+1),MR,0.0,RZMAIL(IB),PP)
      PSS=0.0
      RJ1=0.0
      DO 30 I=1,J2
      PIS(I)=1.0
      RJ2=RAYR2(I+1)**2
      VV=PI*(RJ2-RJ1)
      DO 20 J=1,J2
   20 PIS(I)=PIS(I)-PP(I,J)*SIGT2(J1+J)
      PSS=PSS+PIS(I)*SIGT2(J1+I)*VV/SURFA
   30 RJ1=RJ2
      DEALLOCATE(RAYR2)
      PSS=1.0-SCALE2*PSS
      IF(IMPX.GE.8) CALL SYBPRX(1,J2,JKG,SIGT2(J1+1),PP,PIS,PSS)
      SIGT2(J1+J2)=SIGT2(J1+J2)/SCALE1
*----
*  COMPUTE THE REDUCED COLLISION PROBABILITY MATRIX FOR INCOMING
*  NEUTRONS
*----
      SURFA=(0.5*PI*ROUT)/SCALE2
      RJ1=0.0
      DO 40 I=1,J2-1
      RJ2=PI*RAYRE(J1+I+1)**2
      PSJ(I)=PIS(I)*(RJ2-RJ1)/SURFA
   40 RJ1=RJ2
      RJ2=PI*ROUT**2
      PSJ(J2)=PIS(J2)*(RJ2-RJ1)*SCALE1/SURFA
*----
*  CHECK IF SCATTERING REDUCTION IS REQUIRED
*----
      LSKIP=.TRUE.
      DO 70 I=1,J2
   70 LSKIP=LSKIP.AND.(SIGW2(J1+I).EQ.0.0)
*----
*  SCATTERING REDUCTION IF LSKIP=.FALSE.
*----
      IF(LSKIP) THEN
*        DO NOT PERFORM SCATTERING REDUCTION.
         DO 90 I=1,J2
         DO 80 J=1,J2-1
   80    PIJW(IPIJ+(J-1)*J2+I)=PP(I,J)
   90    PIJW(IPIJ+(J2-1)*J2+I)=PP(I,J2)*SCALE1
         DO 100 I=1,J2
         PISW(J1+I)=PIS(I)
  100    PSJW(J1+I)=PSJ(I)
         PSSW(JKG)=PSS
      ELSE
*        COMPUTE THE SCATTERING-REDUCED COLLISION AND ESCAPE MATRICES.
         DO 120 I=1,J2
         DO 110 J=1,J2-1
  110    PIJW(IPIJ+(J-1)*J2+I)=-PP(I,J)*SIGW2(J1+J)
         PIJW(IPIJ+(J2-1)*J2+I)=-PP(I,J2)*SIGW2(J1+J2)*SCALE1
  120    PIJW(IPIJ+(I-1)*J2+I)=1.0+PIJW(IPIJ+(I-1)*J2+I)
         CALL ALINV(J2,PIJW(IPIJ+1),J2,IER)
         IF(IER.NE.0) CALL XABORT('SYB002: SINGULAR MATRIX.')
         ALLOCATE(WORK(J2))
         DO 170 I=1,J2
         DO 130 K=1,J2
  130    WORK(K)=PIJW(IPIJ+(K-1)*J2+I)
         DO 150 J=1,J2-1
         WGAR=0.0
         DO 140 K=1,J2
  140    WGAR=WGAR+WORK(K)*PP(K,J)
  150    PIJW(IPIJ+(J-1)*J2+I)=WGAR
         WGAR=0.0
         DO 155 K=1,J2
  155    WGAR=WGAR+WORK(K)*PP(K,J2)
         PIJW(IPIJ+(J2-1)*J2+I)=WGAR*SCALE1
         WGAR=0.0
         DO 160 K=1,J2
  160    WGAR=WGAR+WORK(K)*PIS(K)
  170    PISW(J1+I)=WGAR
         DEALLOCATE(WORK)
*
*        COMPUTE THE SCATTERING-REDUCED COLLISION PROBABILITY MATRIX
*        FOR INCOMING NEUTRONS.
         DO 190 J=1,J2
         WGAR=PSJ(J)
         DO 180 K=1,J2
  180    WGAR=WGAR+PSJ(K)*SIGW2(J1+K)*PIJW(IPIJ+(J-1)*J2+K)
  190    PSJW(J1+J)=WGAR
*
*        COMPUTE THE SCATTERING-REDUCED TRANSMISSION PROBABILITY MATRIX.
         WGAR=PSS
         DO 200 K=1,J2
  200    WGAR=WGAR+PSJ(K)*SIGW2(J1+K)*PISW(J1+K)
         PSSW(JKG)=WGAR
      ENDIF
      DEALLOCATE(PP)
  220 IPIJ=IPIJ+J2*J2
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(PSJ,PIS)
      RETURN
      END
