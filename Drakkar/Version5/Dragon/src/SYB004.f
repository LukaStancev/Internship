*DECK SYB004
      SUBROUTINE SYB004 (NGEN,NPIJ,NPIS,NRAYRE,SIGT2,SIGW2,IMPX,NCOUR,
     1 IQUAD,XX,YY,LSECT,NMC,NMCR,RAYRE,MAIL,IZMAIL,RZMAIL,PIJW,PISW,
     2 PSJW,PSSW)
*
*-----------------------------------------------------------------------
*
*Purpose:
* compute the cellwise scattering-reduced collision, escape and
* transmission probabilities in a 2-D Cartesian or hexagonal assembly
* with DP-0 approximation.
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
* NRAYRE  size of array RAYRE (NRAYRE=NMCR(NGEN+1)).
* SIGT2   total macroscopic cross sections.
* SIGW2   P0 within-group scattering macroscopic cross sections.
* IMPX    print flag (equal to 0 for no print).
* NCOUR   number of currents surrounding the cells (=4: Cartesian
*         lattice; =6: hexagonal lattice).
* IQUAD   quadrature parameters.
* XX      X-thickness of the generating cells.
* YY      Y-thickness of the generating cells.
* LSECT   type of sectorization.
* NMC     offset of the first volume in each generating cell.
* NMCR    offset of the first radius in each generating cell.
* RAYRE   radius of the tubes in each generating cell.
* MAIL    offset of the first tracking information in each generating
*         cell.
* IZMAIL  integer tracking information.
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
      INTEGER NGEN,NPIJ,NPIS,NRAYRE,IMPX,NCOUR,IQUAD(4),LSECT(NGEN),
     1 NMC(NGEN+1),NMCR(NGEN+1),MAIL(2,NGEN),IZMAIL(*)
      REAL SIGT2(NPIS),SIGW2(NPIS),XX(NGEN),YY(NGEN),RAYRE(NRAYRE),
     1 RZMAIL(*),PIJW(NPIJ),PISW(NCOUR*NPIS),PSJW(NCOUR*NPIS),
     2 PSSW(NGEN*NCOUR*NCOUR)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (PI=3.141592654,SMALL=5.0E-3,SQRT3=1.732050807568877)
      LOGICAL LSKIP
      REAL PSS(36),SURFA(6),ALPA(64),PWA(64)
      REAL, ALLOCATABLE, DIMENSION(:) :: VOL,WORK
      REAL, ALLOCATABLE, DIMENSION(:,:) :: PIS,PSJ,PP
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(PIS(NCOUR,NPIS),PSJ(NCOUR,NPIS))
*
      IPIJ=0
      IPIS=0
      IPSS=0
      DO 220 JKG=1,NGEN
      J1=NMC(JKG)
      J2=NMC(JKG+1)-J1
      J1R=NMCR(JKG)
      J2R=NMCR(JKG+1)-J1R
      ALLOCATE(PP(J2,J2),VOL(J2))
*----
*  COMPUTE THE REDUCED COLLISION PROBABILITY MATRIX
*----
      A=XX(JKG)
      B=YY(JKG)
      IF((NCOUR.EQ.4).AND.(LSECT(JKG).NE.0)) THEN
*        SECTORIZED CARTESIAN CELL.
         IB1=MAIL(1,JKG)
         IB2=MAIL(2,JKG)
         IF(LSECT(JKG).EQ.-999) THEN
            NSECT=4
         ELSE IF((LSECT(JKG).EQ.-1).OR.(LSECT(JKG).EQ.-101)) THEN
            NSECT=8
         ELSE
            NSECT=4*MOD(ABS(LSECT(JKG)),100)
         ENDIF
         MNA4=4*IQUAD(1)
         CALL SYB4QG(IMPX,1,MNA4,J2R,NSECT,LSECT(JKG),J2,RZMAIL(IB2),
     1   IZMAIL(IB1),A,B,RAYRE(J1R+2),SIGT2(J1+1),SMALL,VOL,PP,PIS,PSS)
      ELSE IF(LSECT(JKG).NE.0) THEN
*        SECTORIZED HEXAGONAL CELL.
         IB1=MAIL(1,JKG)
         IB2=MAIL(2,JKG)
         NSECT=6
         MNA4=12*IQUAD(1)
         CALL SYB7QG(IMPX,1,MNA4,J2R,NSECT,LSECT(JKG),J2,RZMAIL(IB2),
     1   IZMAIL(IB1),A,RAYRE(J1R+2),SIGT2(J1+1),SMALL,VOL,PP,PIS,PSS)
      ELSE IF((NCOUR.EQ.4).AND.(J2.EQ.1)) THEN
         CALL ALGPT(IQUAD(3),-1.0,1.0,ALPA,PWA)
         CALL RECT1(IQUAD(3),A,B,SIGT2(J1+1),SMALL,PP,PIS,PSS,ALPA,PWA)
         VOL(1)=A*B
      ELSE IF(J2.EQ.1) THEN
         CALL ALGPT(IQUAD(3),-1.0,1.0,ALPA,PWA)
         CALL XHX2D0(IQUAD(3),ALPA,PWA,A,SIGT2(J1+1),SMALL,PP,PIS,PSS)
         VOL(1)=1.5*SQRT3*A*A
      ELSE
*        NON-SECTORIZED CARTESIAN OR HEXAGONAL CELL.
         IB1=MAIL(1,JKG)
         IB2=MAIL(2,JKG)
         CALL SYBUP0(RZMAIL(IB2),IZMAIL(IB1),NCOUR,J2,SIGT2(J1+1),SMALL,
     1   A,B,IMPX,VOL,PP,PIS,PSS)
      ENDIF
      IF(IMPX.GE.8) CALL SYBPRX(NCOUR,J2,JKG,SIGT2(J1+1),PP,PIS,PSS)
*----
*  COMPUTE THE REDUCED COLLISION PROBABILITY MATRIX FOR INCOMING
*  NEUTRONS
*----
      DO 60 I=1,J2
      IF(NCOUR.EQ.4) THEN
         SURFA(1)=0.25*B
         SURFA(2)=0.25*B
         SURFA(3)=0.25*A
         SURFA(4)=0.25*A
      ELSE
         DO 50 JC=1,6
   50    SURFA(JC)=0.25*A
      ENDIF
      DO 60 JC=1,NCOUR
   60 PSJ(JC,I)=PIS(JC,I)*VOL(I)/SURFA(JC)
      DEALLOCATE(VOL)
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
         DO 80 I=1,J2
         DO 80 J=1,J2
   80    PIJW(IPIJ+(J-1)*J2+I)=PP(I,J)
         DO 90 I=1,J2
         DO 90 JC=1,NCOUR
         PISW(IPIS+(JC-1)*J2+I)=PIS(JC,I)
   90    PSJW(IPIS+(I-1)*NCOUR+JC)=PSJ(JC,I)
         DO 100 IC=1,NCOUR
         DO 100 JC=1,NCOUR
  100    PSSW(IPSS+(JC-1)*NCOUR+IC)=PSS((JC-1)*NCOUR+IC)
      ELSE
*        COMPUTE THE SCATTERING-REDUCED COLLISION AND ESCAPE MATRICES.
         DO 120 I=1,J2
         DO 110 J=1,J2
  110    PIJW(IPIJ+(J-1)*J2+I)=-PP(I,J)*SIGW2(J1+J)
  120    PIJW(IPIJ+(I-1)*J2+I)=1.0+PIJW(IPIJ+(I-1)*J2+I)
         CALL ALINV(J2,PIJW(IPIJ+1),J2,IER)
         IF(IER.NE.0) CALL XABORT('SYB004: SINGULAR MATRIX.')
         ALLOCATE(WORK(J2))
         DO 170 I=1,J2
         DO 130 K=1,J2
  130    WORK(K)=PIJW(IPIJ+(K-1)*J2+I)
         DO 150 J=1,J2
         WGAR=0.0
         DO 140 K=1,J2
  140    WGAR=WGAR+WORK(K)*PP(K,J)
  150    PIJW(IPIJ+(J-1)*J2+I)=WGAR
         DO 170 JC=1,NCOUR
         WGAR=0.0
         DO 160 K=1,J2
  160    WGAR=WGAR+WORK(K)*PIS(JC,K)
  170    PISW(IPIS+(JC-1)*J2+I)=WGAR
         DEALLOCATE(WORK)
*
*        COMPUTE THE SCATTERING-REDUCED COLLISION PROBABILITY MATRIX
*        FOR INCOMING NEUTRONS.
         DO 190 IC=1,NCOUR
         DO 190 J=1,J2
         WGAR=PSJ(IC,J)
         DO 180 K=1,J2
  180    WGAR=WGAR+PSJ(IC,K)*SIGW2(J1+K)*PIJW(IPIJ+(J-1)*J2+K)
  190    PSJW(IPIS+(J-1)*NCOUR+IC)=WGAR
*
*        COMPUTE THE SCATTERING-REDUCED TRANSMISSION PROBABILITY MATRIX.
         DO 210 IC=1,NCOUR
         DO 210 JC=1,NCOUR
         WGAR=PSS((JC-1)*NCOUR+IC)
         DO 200 K=1,J2
  200    WGAR=WGAR+PSJ(IC,K)*SIGW2(J1+K)*PISW(IPIS+(JC-1)*J2+K)
  210    PSSW(IPSS+(JC-1)*NCOUR+IC)=WGAR
      ENDIF
      DEALLOCATE(PP)
      IPIJ=IPIJ+J2*J2
      IPIS=IPIS+J2*NCOUR
  220 IPSS=IPSS+NCOUR*NCOUR
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(PSJ,PIS)
      RETURN
      END
