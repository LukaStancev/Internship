*DECK SYBRX3
      SUBROUTINE SYBRX3 (MULTC,IPAS,NPIJ,NPIS,NRAYRE,SIGT,SIGW,P,IMPX,
     1 NCOUR,IWIGN,NMCEL,NMERGE,NGEN,IJAT,IQUAD,XX,YY,LSECT,NMC,NMCR,
     2 RAYRE,MAIL,IZMAIL,RZMAIL,IFR,ALB,INUM,MIX,DVX,IGEN)
*
*-----------------------------------------------------------------------
*
*Purpose:
* compute the global scattering-reduced collision probabilities in a
* 2-D Cartesian or hexagonal assembly using the interface current
* method with Roth x 4, Roth x 6, DP-0 or DP-1 approximation.
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
* MULTC   type of interface cuttent approximation:
*         =2: Roth x 4 or Roth x 6 approximation;
*         =3: DP-0 approximation; =4: DP-1 approximation.
* IPAS    total number of volumes.
* NPIJ    length of cellwise scattering-reduced collision probability
*         matrices.
* NPIS    length of cellwise scattering-reduced collision probability
*         matrices (NPIS=NMC(NGEN+1)).
* NRAYRE  size of array rayre (NRAYRE=NMCR(NGEN+1)).
* SIGT    total macroscopic cross sections.
* SIGW    P0 within-group scattering macroscopic cross sections.
* IMPX    print flag (equal to 0 for no print).
* NCOUR   number of currents surrounding the cells (=4 or 12: Cartesian
*         lattice; =6 or 18: hexagonal lattice).
* IWIGN   type of cylinderization if MULTC=2.
* IQUAD   quadrature parameters.
* NMCEL   total number of cells in the domain.
* IFR     index-number of in-currents.
* ALB     transmission/albedo associated with each in-current.
* NMERGE  total number of merged cells for which specific values
*         of the neutron flux and reactions rates are required.
*         Many cells with different position in the domain can
*         be merged before the neutron flux calculation if they
*         own the same generating cell. This allows some reduction
*         in cpu time and memory (NMERGE.le.NMCEL).
* IJAT    total number of distinct out-currents.
* INUM    index-number of the merged cell associated to each cell.
* MIX     index-number of out-currents.
* DVX     weight associated with each out-current.
*         Note: IFR, ALB, MIX and DVX contains information to rebuild
*         the geometrical 'A' matrix.
* NGEN    total number of generating cells. A generating cell is
*         defined by its material and its position in the domain
*         (NGEN.le.NMERGE).
* XX      X-thickness of the generating cells.
* YY      Y-thickness of the generating cells.
* LSECT   type of sectorization.
* NMC     offset of the first volume in each generating cell.
* NMCR    offset of the first radius in each generating cell.
* RAYRE   radius of the tubes in each generating cell.
* MAIL    offset of the first tracking information in each generatin
*         cell.
* IZMAIL  integer tracking information.
* RZMAIL  real tracking information.
* IGEN    index-number of the generating cell associated with each
*         merged cell.
*
*Parameters: output
* P       reduced collision probabilities.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER MULTC,IPAS,NPIJ,NPIS,NRAYRE,IMPX,NCOUR,IWIGN,NMCEL,
     1 NMERGE,NGEN,IJAT,IQUAD(4),LSECT(NGEN),NMC(NGEN+1),NMCR(NGEN+1),
     2 MAIL(2,NGEN),IZMAIL(*),IFR(NCOUR*NMCEL),INUM(NMCEL),
     3 MIX(NCOUR*NMERGE),IGEN(NMERGE)
      REAL SIGT(IPAS),SIGW(IPAS),P(IPAS,IPAS),XX(NGEN),YY(NGEN),
     1 RAYRE(NRAYRE),RZMAIL(*),ALB(NCOUR*NMCEL),DVX(NCOUR*NMERGE)
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: SIGT2,SIGW2,PIJW,PISW,PSJW,
     1 PSSW
      REAL, ALLOCATABLE, DIMENSION(:,:) :: PSSB
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(PSSB(IJAT,2*IJAT),SIGT2(IPAS),SIGW2(IPAS),PIJW(NPIJ),
     1 PISW(NCOUR*NPIS),PSJW(NCOUR*NPIS),PSSW(NGEN*NCOUR*NCOUR))
*
      DO 10 I=1,IPAS
      DO 10 J=1,IPAS
   10 P(I,J)=0.0
      I1=0
      DO 30 IKK=1,NMERGE
      IKG=IGEN(IKK)
      J1=NMC(IKG)
      I2=NMC(IKG+1)-J1
      DO 20 I=1,I2
      SIGT2(J1+I)=SIGT(I1+I)
   20 SIGW2(J1+I)=SIGW(I1+I)
   30 I1=I1+I2
*
      IF(MULTC.EQ.2) THEN
*        ROTH X 4 OR ROTH X 6 APPROXIMATION.
         CALL SYB003 (NGEN,NPIJ,NPIS,SIGT2,SIGW2,IMPX,NCOUR,IWIGN,IQUAD,
     1   XX,YY,NMC,RAYRE,MAIL,RZMAIL,PIJW,PISW,PSJW,PSSW)
      ELSE IF(MULTC.EQ.3) THEN
*        DP-0 APPROXIMATION.
         CALL SYB004 (NGEN,NPIJ,NPIS,NRAYRE,SIGT2,SIGW2,IMPX,NCOUR,
     1   IQUAD,XX,YY,LSECT,NMC,NMCR,RAYRE,MAIL,IZMAIL,RZMAIL,PIJW,PISW,
     2   PSJW,PSSW)
      ELSE IF(MULTC.EQ.4) THEN
*        DP-1 APPROXIMATION.
         CALL SYB005 (NGEN,NPIJ,NPIS,NRAYRE,SIGT2,SIGW2,IMPX,NCOUR,
     1   IQUAD,XX,YY,LSECT,NMC,NMCR,RAYRE,MAIL,IZMAIL,RZMAIL,PIJW,PISW,
     2   PSJW,PSSW)
      ELSE
         CALL XABORT('SYBRX3: UNKNOWN CP MODULE.')
      ENDIF
*
      IPIJ=0
      DO 60 JKG=1,NGEN
      J2=NMC(JKG+1)-NMC(JKG)
      I1=0
      DO 50 IKK=1,NMERGE
      IKG=IGEN(IKK)
      I2=NMC(IKG+1)-NMC(IKG)
      IF(IKG.EQ.JKG) THEN
         DO 40 J=1,J2
         DO 40 I=1,J2
   40    P(I1+I,I1+J)=PIJW(IPIJ+(J-1)*J2+I)
      ENDIF
   50 I1=I1+I2
   60 IPIJ=IPIJ+J2*J2
*----
*  COMPUTATION OF PSSB=A*(I-PSS*A)**-1
*----
      DO 80 I=1,IJAT
      DO 70 J=1,IJAT
      PSSB(I,J)=0.0
   70 PSSB(I,IJAT+J)=0.0
   80 PSSB(I,I)=1.0
      DO 90 ICEL=1,NMCEL
      IKK=INUM(ICEL)
      IT=NCOUR*(IKK-1)
      IS=NCOUR*(ICEL-1)
      IKG=IGEN(IKK)
      IPSS=(IKG-1)*NCOUR*NCOUR
      DO 90 JC=1,NCOUR
      J1=IFR(IS+JC)
      J2=MIX(IT+JC)
      ALBEDO=ALB(IS+JC)
      PSSB(J1,IJAT+J2)=PSSB(J1,IJAT+J2)+ALBEDO*DVX(IT+JC)
      DO 90 IC=1,NCOUR
      J2=MIX(IT+IC)
   90 PSSB(J1,J2)=PSSB(J1,J2)-PSSW(IPSS+(JC-1)*NCOUR+IC)*ALBEDO*
     1 DVX(IT+IC)
      CALL ALSB (IJAT,IJAT,PSSB,IER,IJAT)
      IF(IER.NE.0) CALL XABORT('SYBRX3: SINGULAR MATRIX.')
*----
*  COMPUTATION OF PISW*PSSB*PSJW
*----
      I1=0
      DO 120 IKK=1,NMERGE
      IKG=IGEN(IKK)
      I1P=NMC(IKG)
      I2=NMC(IKG+1)-I1P
      IT=NCOUR*(IKK-1)
      DO 110 I=1,I2
      DO 110 IC=1,NCOUR
      ICC=MIX(IT+IC)
      ZZZ=PISW(I1P*NCOUR+(IC-1)*I2+I)*SIGN(1.0,DVX(IT+IC))
      J1=0
      DO 110 JKK=1,NMERGE
      JKG=IGEN(JKK)
      J1P=NMC(JKG)
      J2=NMC(JKG+1)-J1P
      JT=NCOUR*(JKK-1)
      DO 100 J=1,J2
      DO 100 JC=1,NCOUR
      JCC=MIX(JT+JC)
      PBJ=PSJW(J1P*NCOUR+(J-1)*NCOUR+JC)
  100 P(I1+I,J1+J)=P(I1+I,J1+J)+ZZZ*DVX(JT+JC)*PSSB(JCC,IJAT+ICC)*PBJ
  110 J1=J1+J2
  120 I1=I1+I2
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(PSSW,PSJW,PISW,PIJW,SIGW2,SIGT2,PSSB)
      RETURN
      END
