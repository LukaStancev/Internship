*DECK SNFTH2
      SUBROUTINE SNFTH2(NUN,NGRP,NGEFF,IMPX,INCONV,INGIND,NHEX,IELEM,
     1 ISPLH,NMAT,NPQ,NSCT,MAT,VOL,TOTAL,QEXT,DU,DE,W,DB,DA,PL,CONNEC,
     2 IZGLOB,SIDE,CONFROM,FUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* perform one inner iteration for solving SN equations in 2D hexagonal
* geometry for the DISCONTINUOUS GALERKIN method. Energy-angle
* multithreading. Albedo boundary conditions.
*
*Copyright:
* Copyright (C) 2019 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. A. Calloo
*
*Parameters: input
* NUN     total number of unknowns in vector FUNKNO.
* NGRP    number of energy groups.
* NGEFF   dimension of arrays INCONV and INGIND.
* IMPX    print flag (equal to zero for no print).
* INCONV  energy group convergence flag (set to .false. if converged).
* INGIND  energy group index assign to 1:NGEFF arrays.
* NHEX    number of hexagons in X-Y plane.
* IELEM   measure of order of the spatial approximation polynomial:
*         =1: constant - only for HODD, classical diamond scheme
*             - default for HODD;
*         =2: linear - default for DG;
*         =3: parabolic;
*         =4: cubic - only for DG.
* NMAT    number of material mixtures.
* NPQ     number of SN directions in four octants (including zero-weight
*         directions).
* NSCT    maximum number of spherical harmonics moments of the flux.
* MAT     material mixture index in each region.
* VOL     volumes of each region.
* TOTAL   macroscopic total cross sections.
* NCODE   boundary condition indices.
* ZCODE   albedos.
* QEXT    Legendre components of the fixed source.
* LFIXUP  flag to enable negative flux fixup.
* DU      first direction cosines ($\mu$).
* DE      second direction cosines ($\eta$).
* W       weights.
* MRM     quadrature index.
* MRMY    quadrature index.
* DB      diamond-scheme parameter.
* DA      diamond-scheme parameter.
* PL      discrete values of the spherical harmonics corresponding
*         to the 2D SN quadrature.
* IZGLOB  hexagon sweep order depending on direction
* SIDE    side of an hexagon.
* CONNEC  connectivity matrix for flux swapping -- which lozenges is the
*         lozenge under consideration connected to; in order to pass the
*         flux along. This is dependent on direction
* CONFROM matrix for incoming flux -- which lozenges are feeding into
*         the lozenge under consideration. This is dependent on
*         direction
*
*Parameters: output
* FUNKNO  Legendre components of the flux and boundary fluxes.
*
*-----------------------------------------------------------------------
*
#if defined(_OPENMP)
      USE omp_lib
#endif
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NUN,NGRP,NGEFF,IMPX,INGIND(NGEFF),NHEX,IELEM,ISPLH,NMAT,
     1 NPQ,NSCT,IZGLOB(NHEX,6),MAT(ISPLH,ISPLH,3,NHEX),
     2 CONNEC(3,(NHEX*3)*2,6),CONFROM(2,3,6)
      LOGICAL INCONV(NGEFF)
      REAL SIDE,VOL(ISPLH,ISPLH,3,NHEX), TOTAL(0:NMAT,NGRP),
     1 QEXT(NUN,NGRP), DU(NPQ), DB(ISPLH,ISPLH,3,NHEX,NPQ),
     2 DA(ISPLH,ISPLH,3,NHEX,NPQ), PL(NSCT,NPQ), DE(NPQ), W(NPQ),
     3 FUNKNO(NUN,NGRP)
*----
*  LOCAL VARIABLES
*----
      INTEGER NPQD(6),IIND(6)
      DOUBLE PRECISION :: Q(IELEM**2), Q2(IELEM**2,(IELEM**2)+1),
     1 VT,CORNERQ(IELEM**2)
      PARAMETER(IUNOUT=6,RLOG=1.0E-8,PI=3.141592654)
      INTEGER :: ILOZSWP(3,6), IFROMI, IFROMJ
      REAL :: JAC(2,2,3), MUH, ETAH, AAA, BBB, CCC, DDD, MUHTEMP,
     1 ETAHTEMP
      DOUBLE PRECISION :: THETA, XNI(IELEM,ISPLH), XNJ(IELEM)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: INDANG
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: FLUX
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: FLUX_G
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: BFLUX
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(INDANG(NPQ,6))
      ALLOCATE(BFLUX(2,NHEX*3,ISPLH,IELEM),
     1 FLUX(IELEM**2,NSCT,3*ISPLH**2,NHEX))
      ALLOCATE(FLUX_G(IELEM**2,NSCT,3*ISPLH**2,NHEX,NGEFF))
*----
*  CONSTRUCT JACOBIAN MATRIX FOR EACH LOZENGE
*----
      ILOZSWP = RESHAPE((/ 3, 2, 1, 3, 1, 2, 1, 3, 2, 1, 2, 3, 2, 1,
     1   3, 2, 3, 1 /), SHAPE(ILOZSWP))
      JAC = RESHAPE((/ 1., -SQRT(3.), 1., SQRT(3.), 2., 0., 1.,
     1    SQRT(3.), 2., 0., -1., SQRT(3.) /), SHAPE(JAC))
      JAC = (SIDE/2.)*JAC
*----
*  PARAMETER VALIDATION.
*----
      L4=3*(IELEM**2)*(ISPLH**2)*NHEX*NSCT
      IF((IELEM.LT.1).OR.(IELEM.GT.3))
     1   CALL XABORT('SNFTH2: INVALID IELEM (DIAM) VALUE. '
     2   //'CHECK INPUT DATA FILE.')
      FLUX_G(:IELEM**2,:NSCT,:3*ISPLH**2,:NHEX,:NGEFF)=0.0D0
*----
*  SET OCTANT SWAPPING ORDER.
*----
      NPQD(:6)=0
      INDANG(:NPQ,:6)=0
      DO 10 M=1,NPQ
        VU=DU(M)
        VE=DE(M)
        THETA=0.0D0
        IF(VE.GT.0.0)THEN
          IF(VU.EQ.0.0)THEN
            THETA = PI/2
          ELSEIF(VU.GT.0.0)THEN
            THETA = ATAN(ABS(VE/VU))
          ELSEIF(VU.LT.0.0)THEN
            THETA = PI - ATAN(ABS(VE/VU))
          ENDIF
        ELSEIF(VE.LT.0.0)THEN
          IF(VU.EQ.0.0)THEN
            THETA = 3*PI/2
          ELSEIF(VU.LT.0.0)THEN
            THETA = PI + ATAN(ABS(VE/VU))
          ELSEIF(VU.GT.0.0)THEN
            THETA = 2.*PI - ATAN(ABS(VE/VU))
          ENDIF
        ENDIF
        IND=0
        IF((THETA.GT.0.0).AND.(THETA.LT.(PI/3.)))THEN
          IND=1
        ELSEIF((THETA.GT.(PI/3.)).AND.(THETA.LT.(2.*PI/3.)))THEN
          IND=2
        ELSEIF((THETA.GT.(2.*PI/3.)).AND.(THETA.LT.(PI)))THEN
          IND=3
        ELSEIF((THETA.GT.(PI)).AND.(THETA.LT.(4.*PI/3.)))THEN
          IND=4
        ELSEIF((THETA.GT.(4.*PI/3.)).AND.(THETA.LT.(5.*PI/3.)))THEN
          IND=5
        ELSEIF((THETA.GT.(5.*PI/3.)).AND.(THETA.LT.(2.*PI)))THEN
          IND=6
        ENDIF
        IIND(IND)=IND ! assume IIND(I)=I in hexagonal geometry
        NPQD(IND)=NPQD(IND)+1
        INDANG(NPQD(IND),IND)=M
   10 CONTINUE
*----
*  MAIN LOOP OVER OCTANTS.
*----
      DO 220 JND=1,6
      IND=IIND(JND)
*----
*  MAIN SWAPPING LOOPS FOR SN FLUX CALCULATION
*----
*$OMP  PARALLEL DO
*$OMP1 PRIVATE(ITID,FLUX,BFLUX,M,IG,XNI,XNJ,Q,Q2,IOF,IER,II,JJ,IEL,JEL)
*$OMP2 PRIVATE(I,J,K,IL,JL,VT,AAA,BBB,CCC,DDD,IHEXI,ILOZI,ISIDEI,IHEXJ)
*$OMP3 PRIVATE(ILOZJ,ISIDEJ,IFROMI,IFROMJ,INDEXI,INDEXJ,MUHTEMP)
*$OMP4 PRIVATE(ETAHTEMP,MUH,ETAH,I2,I3,J2,J3,C1,ILOZ,CORNERQ)
*$OMP5 SHARED(FUNKNO) REDUCTION(+:FLUX_G) COLLAPSE(2)
      DO 210 IGEFF=1,NGEFF
      DO 200 IPQD=1,NPQD(IND)
#if defined(_OPENMP)
        ITID=omp_get_thread_num()
#else
        ITID=0
#endif
      IF(IMPX.GT.5) WRITE(IUNOUT,400) ITID,IGEFF,IPQD
      IF(.NOT.INCONV(IGEFF)) GO TO 200
      M=INDANG(IPQD,IND)
      IG=INGIND(IGEFF)
      FLUX(:IELEM**2,:NSCT,:3*ISPLH**2,:NHEX)=0.0D0
      BFLUX(:2,:NHEX*3,:ISPLH,:IELEM)=0.0D0
      IF(W(M).EQ.0.0) GO TO 200
*----
*  LOOP OVER X- AND Y-DIRECTED AXES.
*----
* SECOND loop over hexagons
      DO 190 IZ=1,NHEX
      I=IZGLOB(IZ,IND)

* THIRD loop over lozenges
      DO 180 JJ=1,3
      J=ILOZSWP(JJ,IND)
*
      AAA = JAC(1,1,J)
      BBB = JAC(1,2,J)
      CCC = JAC(2,1,J)
      DDD = JAC(2,2,J)
*
      IHEXI  = CONNEC(1,((I-1)*3*2) + ((J-1)*2) +1,IND)
      ILOZI  = CONNEC(2,((I-1)*3*2) + ((J-1)*2) +1,IND)
      ISIDEI = CONNEC(3,((I-1)*3*2) + ((J-1)*2) +1,IND)
      IHEXJ  = CONNEC(1,((I-1)*3*2) + ((J-1)*2) +2,IND)
      ILOZJ  = CONNEC(2,((I-1)*3*2) + ((J-1)*2) +2,IND)
      ISIDEJ = CONNEC(3,((I-1)*3*2) + ((J-1)*2) +2,IND)
      IFROMI = CONFROM(1,J,IND)
      IFROMJ = CONFROM(2,J,IND)
      INDEXI = ((IHEXI-1)*3)+ILOZI
      INDEXJ = ((IHEXJ-1)*3)+ILOZJ
*
      DO 170 IL=1,ISPLH
      I2=IL
      IF(CONFROM(1,J,IND).EQ.3) I2=ISPLH+1-I2

      DO IEL=1,IELEM
         XNJ(IEL) = BFLUX(2,((I-1)*3)+J,I2,IEL)
      ENDDO

      DO 160 JL=1,ISPLH
      J2=JL
      IF(CONFROM(2,J,IND).EQ.4) J2=ISPLH+1-J2

      IF(IL.EQ.1)THEN
         DO IEL=1,IELEM
            XNI(IEL,J2) = BFLUX(1,((I-1)*3)+J,J2,IEL)
         ENDDO
      ENDIF
*
      MUHTEMP  =  DA(I2,J2,J,I,M)
      ETAHTEMP =  DB(I2,J2,J,I,M)
      MUH = (MUHTEMP*DDD) - (ETAHTEMP*BBB)
      ETAH = (-MUHTEMP*CCC) + (ETAHTEMP*AAA)
*
      IF(MAT(I2,J2,J,I).EQ.0) GO TO 180

*     -----------------------------------------------------
      DO 110 IEL=1,IELEM**2
      Q(IEL)=0.0D0
      DO 100 K=1,NSCT
      IOF=((((((I-1)*3+(J-1))*ISPLH+(J2-1))*ISPLH+(I2-1))*NSCT+
     1 (K-1))*IELEM**2)+IEL
      Q(IEL)=Q(IEL)+QEXT(IOF,IG)*PL(K,M)/(4.0D0*PI)
  100 CONTINUE
  110 CONTINUE
*     -----------------------------------------------------
*
      VT=VOL(I2,J2,J,I)*TOTAL(MAT(I2,J2,J,I),IG)
      CALL XDDSET(Q2,(IELEM**2)*((IELEM**2)+1),0.0D0)
*
*     -----------------------------------------------------
      IF(IELEM.EQ.2) THEN
*     alpha=1 , beta=1
      Q2(1,1) = VT/4 + ABS(ETAH)/4 + ABS(MUH)/4
      Q2(1,2) = (3**(1.0D0/2.0D0)*MUH)/4
      Q2(1,3) = (3**(1.0D0/2.0D0)*ETAH)/4
      Q2(1,4) = 0

*     alpha=2 , beta=1
      Q2(2,1) = -(3**(1.0D0/2.0D0)*MUH)/12
      Q2(2,2) = VT/12 + ABS(ETAH)/12 + ABS(MUH)/4
      Q2(2,3) = 0
      Q2(2,4) = (3**(1.0D0/2.0D0)*ETAH)/12

*     alpha=1 , beta=2
      Q2(3,1) = -(3**(1.0D0/2.0D0)*ETAH)/12
      Q2(3,2) = 0
      Q2(3,3) = VT/12 + ABS(ETAH)/4 + ABS(MUH)/12
      Q2(3,4) = (3**(1.0D0/2.0D0)*MUH)/12

*     alpha=2 , beta=2
      Q2(4,1) = 0
      Q2(4,2) = -(3**(1.0D0/2.0D0)*ETAH)/36
      Q2(4,3) = -(3**(1.0D0/2.0D0)*MUH)/36
      Q2(4,4) = VT/36 + ABS(ETAH)/12 + ABS(MUH)/12

*     source
      Q2(1,5) = Q(1)*VOL(I2,J2,J,I)/4.0D0
      Q2(2,5) = Q(2)*VOL(I2,J2,J,I)/12.0D0
      Q2(3,5) = Q(3)*VOL(I2,J2,J,I)/12.0D0
      Q2(4,5) = Q(4)*VOL(I2,J2,J,I)/36.0D0

*     jump terms (dependent on direction)
      Q2(1,5)= Q2(1,5)
     >   + (1.0D0/8.0D0)*(
     >   MUH*SIGN(1.0,MUH)*(XNI(1,J2)+XNI(2,J2)) +
     >   ETAH*SIGN(1.0,ETAH)*(XNJ(1)+XNJ(2))   )
      Q2(2,5)= Q2(2,5)
     >   -SQRT(3.0D0)*(1.0D0/72.0D0)*(
     >   3.0D0*MUH*(XNI(1,J2)+XNI(2,J2)) +
     >   ETAH*SIGN(1.0,ETAH)*(XNJ(1)-XNJ(2))   )
      Q2(3,5)= Q2(3,5)
     >   -SQRT(3.0D0)*(1.0D0/72.0D0)*(
     >   MUH*SIGN(1.0,MUH)*(XNI(1,J2)-XNI(2,J2)) +
     >   3.0D0*ETAH*(XNJ(1)+XNJ(2))   )
      Q2(4,5)= Q2(4,5)
     >   + (1.0D0/72.0D0)*(
     >   MUH*(XNI(1,J2)-XNI(2,J2)) +
     >   ETAH*(XNJ(1)-XNJ(2))   )
*** ---------------------------------------------------------------- ***
      ELSEIF(IELEM.EQ.3) THEN
      Q2(1,1) = VT/4 + ABS(ETAH)/12 + ABS(MUH)/12
      Q2(1,2) = (5*3**(0.5D0)*MUH)/12
      Q2(1,3) = -(5**(0.5D0)*(6*VT + 2*ABS(ETAH) -
     > 15*ABS(MUH)))/180
      Q2(1,4) = (5*3**(0.5D0)*ETAH)/12
      Q2(1,5) = 0
      Q2(1,6) = -(15**(0.5D0)*ETAH)/18
      Q2(1,7) = -(5**(0.5D0)*(6*VT - 15*ABS(ETAH) +
     > 2*ABS(MUH)))/180
      Q2(1,8) = -(15**(0.5D0)*MUH)/18
      Q2(1,9) = VT/45 - ABS(ETAH)/18 - ABS(MUH)/18

      Q2(2,1) = -(3**(0.5D0)*MUH)/12
      Q2(2,2) = VT/12 + ABS(ETAH)/36 + ABS(MUH)/4
      Q2(2,3) = (15**(0.5D0)*MUH)/12
      Q2(2,4) = 0
      Q2(2,5) = (5*3**(0.5D0)*ETAH)/36
      Q2(2,6) = 0
      Q2(2,7) = (15**(0.5D0)*MUH)/90
      Q2(2,8) = -(5**(0.5D0)*(2*VT - 5*ABS(ETAH) +
     > 6*ABS(MUH)))/180
      Q2(2,9) = -(3**(0.5D0)*MUH)/18

      Q2(3,1) = -(5**(0.5D0)*(3*VT + ABS(ETAH) -
     > 3*ABS(MUH)))/90
      Q2(3,2) = -(15**(0.5D0)*MUH)/10
      Q2(3,3) = VT/15 + ABS(ETAH)/45 + ABS(MUH)/6
      Q2(3,4) = -(15**(0.5D0)*ETAH)/18
      Q2(3,5) = 0
      Q2(3,6) = (3**(0.5D0)*ETAH)/9
      Q2(3,7) = VT/45 - ABS(ETAH)/18 - ABS(MUH)/45
      Q2(3,8) = (3**(0.5D0)*MUH)/15
      Q2(3,9) = -(5**(0.5D0)*(2*VT - 5*ABS(ETAH) +
     > 5*ABS(MUH)))/225

      Q2(4,1) = -(3**(0.5D0)*ETAH)/12
      Q2(4,2) = 0
      Q2(4,3) = (15**(0.5D0)*ETAH)/90
      Q2(4,4) = VT/12 + ABS(ETAH)/4 + ABS(MUH)/36
      Q2(4,5) = (5*3**(0.5D0)*MUH)/36
      Q2(4,6) = -(5**(0.5D0)*(2*VT + 6*ABS(ETAH) -
     > 5*ABS(MUH)))/180
      Q2(4,7) = (15**(0.5D0)*ETAH)/12
      Q2(4,8) = 0
      Q2(4,9) = -(3**(0.5D0)*ETAH)/18

      Q2(5,1) = 0
      Q2(5,2) = -(3**(0.5D0)*ETAH)/36
      Q2(5,3) = 0
      Q2(5,4) = -(3**(0.5D0)*MUH)/36
      Q2(5,5) = VT/36 + ABS(ETAH)/12 + ABS(MUH)/12
      Q2(5,6) = (15**(0.5D0)*MUH)/36
      Q2(5,7) = 0
      Q2(5,8) = (15**(0.5D0)*ETAH)/36
      Q2(5,9) = 0

      Q2(6,1) = (15**(0.5D0)*ETAH)/90
      Q2(6,2) = 0
      Q2(6,3) = -(3**(0.5D0)*ETAH)/45
      Q2(6,4) = -(5**(0.5D0)*(VT + 3*ABS(ETAH) -
     > ABS(MUH)))/90
      Q2(6,5) = -(15**(0.5D0)*MUH)/30
      Q2(6,6) = VT/45 + ABS(ETAH)/15 + ABS(MUH)/18
      Q2(6,7) = -(3**(0.5D0)*ETAH)/18
      Q2(6,8) = 0
      Q2(6,9) = (15**(0.5D0)*ETAH)/45

      Q2(7,1) = -(5**(0.5D0)*(3*VT - 3*ABS(ETAH) +
     > ABS(MUH)))/90
      Q2(7,2) = -(15**(0.5D0)*MUH)/18
      Q2(7,3) = VT/45 - ABS(ETAH)/45 - ABS(MUH)/18
      Q2(7,4) = -(15**(0.5D0)*ETAH)/10
      Q2(7,5) = 0
      Q2(7,6) = (3**(0.5D0)*ETAH)/15
      Q2(7,7) = VT/15 + ABS(ETAH)/6 + ABS(MUH)/45
      Q2(7,8) = (3**(0.5D0)*MUH)/9
      Q2(7,9) = -(5**(0.5D0)*(2*VT + 5*ABS(ETAH) -
     > 5*ABS(MUH)))/225

      Q2(8,1) = (15**(0.5D0)*MUH)/90
      Q2(8,2) = -(5**(0.5D0)*(VT - ABS(ETAH) +
     > 3*ABS(MUH)))/90
      Q2(8,3) = -(3**(0.5D0)*MUH)/18
      Q2(8,4) = 0
      Q2(8,5) = -(15**(0.5D0)*ETAH)/30
      Q2(8,6) = 0
      Q2(8,7) = -(3**(0.5D0)*MUH)/45
      Q2(8,8) = VT/45 + ABS(ETAH)/18 + ABS(MUH)/15
      Q2(8,9) = (15**(0.5D0)*MUH)/45

      Q2(9,1) = VT/45 - ABS(ETAH)/45 - ABS(MUH)/45
      Q2(9,2) = (3**(0.5D0)*MUH)/15
      Q2(9,3) = -(5**(0.5D0)*(2*VT - 2*ABS(ETAH) +
     > 5*ABS(MUH)))/225
      Q2(9,4) = (3**(0.5D0)*ETAH)/15
      Q2(9,5) = 0
      Q2(9,6) = -(2*15**(0.5D0)*ETAH)/75
      Q2(9,7) = -(5**(0.5D0)*(2*VT + 5*ABS(ETAH) -
     > 2*ABS(MUH)))/225
      Q2(9,8) = -(2*15**(0.5D0)*MUH)/75
      Q2(9,9) = (4*VT)/225 + (2*ABS(ETAH))/45 +
     > (2*ABS(MUH))/45

      Q2(1,10) = (((45*Q(01) + 4*Q(09) - 6*5**(0.5D0)*Q(03) -
     > 6*5**(0.5D0)*Q(07)))/180)*VOL(I2,J2,J,I)
      Q2(2,10) = (((15*Q(02) - 2*5**(0.5D0)*Q(08)))/180)*VOL(I2,J2,J,I)
      Q2(3,10) = (((30*Q(03) + 10*Q(07) - 15*5**(0.5D0)*Q(01) -
     > 4*5**(0.5D0)*Q(09)))/450)*VOL(I2,J2,J,I)
      Q2(4,10) = (((15*Q(04) - 2*5**(0.5D0)*Q(06)))/180)*VOL(I2,J2,J,I)
      Q2(5,10) = ((Q(05))/36)*VOL(I2,J2,J,I)
      Q2(6,10) = (((2*Q(06) - 5**(0.5D0)*Q(04)))/90)*VOL(I2,J2,J,I)
      Q2(7,10) = (((10*Q(03) + 30*Q(07) - 15*5**(0.5D0)*Q(01) -
     > 4*5**(0.5D0)*Q(09)))/450)*VOL(I2,J2,J,I)
      Q2(8,10) = (((2*Q(08) - 5**(0.5D0)*Q(02)))/90)*VOL(I2,J2,J,I)
      Q2(9,10) = (((900*Q(01) + 720*Q(09) - 360*5**(0.5D0)*Q(03) -
     > 360*5**(0.5D0)*Q(07)))/40500)*VOL(I2,J2,J,I)

      Q2(1,10) = Q2(1,10) +
     >   ((11*XNI(3,J2))/1080 + (17*XNI(2,J2))/270 + (11*XNI(1,J2))/
     >   1080)*MUH*SIGN(1.0,MUH) + ((11*XNJ(1))/1080 +
     >   (17*XNJ(2))/270 + (11*XNJ(3))/1080)*ETAH*SIGN(1.0,ETAH)
      Q2(2,10) = Q2(2,10) +
     >   (-(3**(0.5D0)*(11*XNI(3,J2) + 68*XNI(2,J2) + 11*XNI(1,J2)))/
     >   1080)*MUH + (-(3**(0.5D0)*(5*XNJ(1) - 5*XNJ(3)))/
     >   1080)*ETAH*SIGN(1.0,ETAH)
      Q2(3,10) = Q2(3,10) +
     >   ((5**(0.5D0)*(11*XNI(3,J2) + 68*XNI(2,J2) + 11*XNI(1,J2)))/
     >   2700)*MUH*SIGN(1.0,MUH) + (-(5**(0.5D0)*(XNJ(1) +
     >   28*XNJ(2) + XNJ(3)))/2700)*ETAH*SIGN(1.0,ETAH)
      Q2(4,10) = Q2(4,10) +
     >   ((3**(0.5D0)*(5*XNI(3,J2) - 5*XNI(1,J2)))/1080)*MUH*
     >   SIGN(1.0,MUH) + (-(3**(0.5D0)*(11*XNJ(1) + 68*XNJ(2) +
     >   11*XNJ(3)))/1080)*ETAH
      Q2(5,10) = Q2(5,10) +
     >   (XNI(1,J2)/72 - XNI(3,J2)/72)*MUH + (XNJ(1)/72 -
     >   XNJ(3)/72)*ETAH
      Q2(6,10) = Q2(6,10) +
     >   ((3**(0.5D0)*5**(0.5D0)*(5*XNI(3,J2) - 5*XNI(1,J2)))/2700)*
     >   MUH*SIGN(1.0,MUH) + ((3**(0.5D0)*5**(0.5D0)*
     >   (XNJ(1) + 28*XNJ(2) + XNJ(3)))/2700)*ETAH
      Q2(7,10) = Q2(7,10) +
     >   (-(5**(0.5D0)*(XNI(3,J2) + 28*XNI(2,J2) + XNI(1,J2)))/
     >   2700)*MUH*SIGN(1.0,MUH) + ((5**(0.5D0)*(11*XNJ(1) +
     >   68*XNJ(2) + 11*XNJ(3)))/2700)*ETAH*SIGN(1.0,ETAH)
      Q2(8,10) = Q2(8,10) +
     >   ((3**(0.5D0)*5**(0.5D0)*(XNI(3,J2) + 28*XNI(2,J2) +
     >   XNI(1,J2)))/2700)*MUH + (-(3**(0.5D0)*5**(0.5D0)*
     >   (5*XNJ(1) - 5*XNJ(3)))/2700)*ETAH*SIGN(1.0,ETAH)
      Q2(9,10) = Q2(9,10) +
     >   (- XNI(3,J2)/1350 - (14*XNI(2,J2))/675 - XNI(1,J2)/1350)*
     >   MUH*SIGN(1.0,MUH) + (- XNJ(1)/1350 - (14*XNJ(2))/
     >   675 - XNJ(3)/1350)*ETAH*SIGN(1.0,ETAH)
*** ---------------------------------------------------------------- ***
      ELSEIF(IELEM.EQ.4) THEN
      Q2(1,1) = (25*VT)/256 + (5*ABS(ETAH))/128 +
     >   (5*ABS(MUH))/128
      Q2(1,2) = (5*3**(0.5D0)*MUH)/32
      Q2(1,3) = -(5**(0.5D0)*(15*VT + 6*ABS(ETAH) -
     >   50*ABS(MUH)))/1280
      Q2(1,4) = (5*7**(0.5D0)*MUH)/128
      Q2(1,5) = (5*3**(0.5D0)*ETAH)/32
      Q2(1,6) = 0
      Q2(1,7) = -(3*15**(0.5D0)*ETAH)/160
      Q2(1,8) = 0
      Q2(1,9) = -(5**(0.5D0)*(15*VT - 50*ABS(ETAH) +
     >   6*ABS(MUH)))/1280
      Q2(1,10) = -(3*15**(0.5D0)*MUH)/160
      Q2(1,11) = (9*VT)/1280 - (3*ABS(ETAH))/128 -
     >   (3*ABS(MUH))/128
      Q2(1,12) = -(3*35**(0.5D0)*MUH)/640
      Q2(1,13) = (5*7**(0.5D0)*ETAH)/128
      Q2(1,14) = 0
      Q2(1,15) = -(3*35**(0.5D0)*ETAH)/640
      Q2(1,16) = 0

      Q2(2,1) = -(11*3**(0.5D0)*MUH)/384
      Q2(2,2) = (85*VT)/768 + (17*ABS(ETAH))/384 +
     >   (11*ABS(MUH))/128
      Q2(2,3) = (37*15**(0.5D0)*MUH)/192
      Q2(2,4) = -(3**(0.5D0)*7**(0.5D0)*(45*VT +
     >   18*ABS(ETAH) - 110*ABS(MUH)))/3840
      Q2(2,5) = 0
      Q2(2,6) = (17*3**(0.5D0)*ETAH)/96
      Q2(2,7) = 0
      Q2(2,8) = -(9*7**(0.5D0)*ETAH)/160
      Q2(2,9) = (11*15**(0.5D0)*MUH)/3200
      Q2(2,10) = -(5**(0.5D0)*(255*VT - 850*ABS(ETAH) +
     >   198*ABS(MUH)))/19200
      Q2(2,11) = -(37*3**(0.5D0)*MUH)/320
      Q2(2,12) = -(3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(30*ABS(ETAH) -
     >   9*VT + 22*ABS(MUH)))/6400
      Q2(2,13) = 0
      Q2(2,14) = (17*7**(0.5D0)*ETAH)/384
      Q2(2,15) = 0
      Q2(2,16) = -(21*3**(0.5D0)*ETAH)/640

      Q2(3,1) = -(3*5**(0.5D0)*(5*VT + 2*ABS(ETAH) -
     >   10*ABS(MUH)))/1280
      Q2(3,2) = -(3*15**(0.5D0)*MUH)/64
      Q2(3,3) = (9*VT)/256 + (9*ABS(ETAH))/640 +
     >   (15*ABS(MUH))/128
      Q2(3,4) = (3*35**(0.5D0)*MUH)/128
      Q2(3,5) = -(3*15**(0.5D0)*ETAH)/160
      Q2(3,6) = 0
      Q2(3,7) = (9*3**(0.5D0)*ETAH)/160
      Q2(3,8) = 0
      Q2(3,9) = (9*VT)/1280 - (3*ABS(ETAH))/128 -
     >   (9*ABS(MUH))/640
      Q2(3,10) = (9*3**(0.5D0)*MUH)/320
      Q2(3,11) = -(9*5**(0.5D0)*(3*VT - 10*ABS(ETAH) +
     >   10*ABS(MUH)))/6400
      Q2(3,12) = -(9*7**(0.5D0)*MUH)/640
      Q2(3,13) = -(3*35**(0.5D0)*ETAH)/640
      Q2(3,14) = 0
      Q2(3,15) = (9*7**(0.5D0)*ETAH)/640
      Q2(3,16) = 0

      Q2(4,1) = -(9*7**(0.5D0)*MUH)/896
      Q2(4,2) = -(3*3**(0.5D0)*7**(0.5D0)*(35*VT +
     >   14*ABS(ETAH) - 30*ABS(MUH)))/8960
      Q2(4,3) = -(9*35**(0.5D0)*MUH)/112
      Q2(4,4) = (81*VT)/1792 + (81*ABS(ETAH))/4480 +
     >   (9*ABS(MUH))/128
      Q2(4,5) = 0
      Q2(4,6) = -(9*7**(0.5D0)*ETAH)/160
      Q2(4,7) = 0
      Q2(4,8) = (81*3**(0.5D0)*ETAH)/1120
      Q2(4,9) = (27*35**(0.5D0)*MUH)/22400
      Q2(4,10) = -(3*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*
     >   (70*ABS(ETAH) - 21*VT + 18*ABS(MUH)))/44800
      Q2(4,11) = (27*7**(0.5D0)*MUH)/560
      Q2(4,12) = -(27*5**(0.5D0)*(9*VT - 30*ABS(ETAH) +
     >   14*ABS(MUH)))/44800
      Q2(4,13) = 0
      Q2(4,14) = -(21*3**(0.5D0)*ETAH)/640
      Q2(4,15) = 0
      Q2(4,16) = (81*7**(0.5D0)*ETAH)/4480

      Q2(5,1) = -(11*3**(0.5D0)*ETAH)/384
      Q2(5,2) = 0
      Q2(5,3) = (11*15**(0.5D0)*ETAH)/3200
      Q2(5,4) = 0
      Q2(5,5) = (85*VT)/768 + (11*ABS(ETAH))/128 +
     >   (17*ABS(MUH))/384
      Q2(5,6) = (17*3**(0.5D0)*MUH)/96
      Q2(5,7) = -(5**(0.5D0)*(255*VT + 198*ABS(ETAH) -
     >   850*ABS(MUH)))/19200
      Q2(5,8) = (17*7**(0.5D0)*MUH)/384
      Q2(5,9) = (37*15**(0.5D0)*ETAH)/192
      Q2(5,10) = 0
      Q2(5,11) = -(37*3**(0.5D0)*ETAH)/320
      Q2(5,12) = 0
      Q2(5,13) = -(3**(0.5D0)*7**(0.5D0)*(45*VT -
     >   110*ABS(ETAH) + 18*ABS(MUH)))/3840
      Q2(5,14) = -(9*7**(0.5D0)*MUH)/160
      Q2(5,15) = -(3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(22*ABS(ETAH) -
     >   9*VT + 30*ABS(MUH)))/6400
      Q2(5,16) = -(21*3**(0.5D0)*MUH)/640

      Q2(6,1) = 0
      Q2(6,2) = -(187*3**(0.5D0)*ETAH)/5760
      Q2(6,3) = 0
      Q2(6,4) = (33*7**(0.5D0)*ETAH)/3200
      Q2(6,5) = -(187*3**(0.5D0)*MUH)/5760
      Q2(6,6) = (289*VT)/2304 + (187*ABS(ETAH))/1920 +
     >   (187*ABS(MUH))/1920
      Q2(6,7) = (629*15**(0.5D0)*MUH)/2880
      Q2(6,8) = -(3**(0.5D0)*7**(0.5D0)*(765*VT +
     >   594*ABS(ETAH) - 1870*ABS(MUH)))/57600
      Q2(6,9) = 0
      Q2(6,10) = (629*15**(0.5D0)*ETAH)/2880
      Q2(6,11) = 0
      Q2(6,12) = -(111*35**(0.5D0)*ETAH)/1600
      Q2(6,13) = (33*7**(0.5D0)*MUH)/3200
      Q2(6,14) = -(3**(0.5D0)*7**(0.5D0)*(765*VT -
     >   1870*ABS(ETAH) + 594*ABS(MUH)))/57600
      Q2(6,15) = -(111*35**(0.5D0)*MUH)/1600
      Q2(6,16) = (189*VT)/6400 - (231*ABS(ETAH))/3200 -
     >   (231*ABS(MUH))/3200

      Q2(7,1) = (11*15**(0.5D0)*ETAH)/3200
      Q2(7,2) = 0
      Q2(7,3) = -(33*3**(0.5D0)*ETAH)/3200
      Q2(7,4) = 0
      Q2(7,5) = -(5**(0.5D0)*(85*VT + 66*ABS(ETAH) -
     >   170*ABS(MUH)))/6400
      Q2(7,6) = -(17*15**(0.5D0)*MUH)/320
      Q2(7,7) = (51*VT)/1280 + (99*ABS(ETAH))/3200 +
     >   (17*ABS(MUH))/128
      Q2(7,8) = (17*35**(0.5D0)*MUH)/640
      Q2(7,9) = -(37*3**(0.5D0)*ETAH)/320
      Q2(7,10) = 0
      Q2(7,11) = (111*15**(0.5D0)*ETAH)/1600
      Q2(7,12) = 0
      Q2(7,13) = -(3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*
     >   (22*ABS(ETAH) - 9*VT + 18*ABS(MUH)))/6400
      Q2(7,14) = (27*35**(0.5D0)*MUH)/1600
      Q2(7,15) = -(3*3**(0.5D0)*7**(0.5D0)*(9*VT -
     >   22*ABS(ETAH) + 30*ABS(MUH)))/6400
      Q2(7,16) = -(63*15**(0.5D0)*MUH)/3200

      Q2(8,1) = 0
      Q2(8,2) = (33*7**(0.5D0)*ETAH)/3200
      Q2(8,3) = 0
      Q2(8,4) = -(297*3**(0.5D0)*ETAH)/22400
      Q2(8,5) = -(51*7**(0.5D0)*MUH)/4480
      Q2(8,6) = -(3**(0.5D0)*7**(0.5D0)*(595*VT +
     >   462*ABS(ETAH) - 510*ABS(MUH)))/44800
      Q2(8,7) = -(51*35**(0.5D0)*MUH)/560
      Q2(8,8) = (459*VT)/8960 + (891*ABS(ETAH))/22400 +
     >   (51*ABS(MUH))/640
      Q2(8,9) = 0
      Q2(8,10) = -(111*35**(0.5D0)*ETAH)/1600
      Q2(8,11) = 0
      Q2(8,12) = (999*15**(0.5D0)*ETAH)/11200
      Q2(8,13) = (27*3**(0.5D0)*MUH)/3200
      Q2(8,14) = (189*VT)/6400 - (231*ABS(ETAH))/3200 -
     >   (81*ABS(MUH))/3200
      Q2(8,15) = (27*15**(0.5D0)*MUH)/400
      Q2(8,16) = -(27*3**(0.5D0)*7**(0.5D0)*(9*VT -
     >   22*ABS(ETAH) + 14*ABS(MUH)))/44800

      Q2(9,1) = -(3*5**(0.5D0)*(5*VT - 10*ABS(ETAH) +
     >   2*ABS(MUH)))/1280
      Q2(9,2) = -(3*15**(0.5D0)*MUH)/160
      Q2(9,3) = (9*VT)/1280 - (9*ABS(ETAH))/640 -
     >   (3*ABS(MUH))/128
      Q2(9,4) = -(3*35**(0.5D0)*MUH)/640
      Q2(9,5) = -(3*15**(0.5D0)*ETAH)/64
      Q2(9,6) = 0
      Q2(9,7) = (9*3**(0.5D0)*ETAH)/320
      Q2(9,8) = 0
      Q2(9,9) = (9*VT)/256 + (15*ABS(ETAH))/128 +
     >   (9*ABS(MUH))/640
      Q2(9,10) = (9*3**(0.5D0)*MUH)/160
      Q2(9,11) = -(9*5**(0.5D0)*(3*VT + 10*ABS(ETAH) -
     >   10*ABS(MUH)))/6400
      Q2(9,12) = (9*7**(0.5D0)*MUH)/640
      Q2(9,13) = (3*35**(0.5D0)*ETAH)/128
      Q2(9,14) = 0
      Q2(9,15) = -(9*7**(0.5D0)*ETAH)/640
      Q2(9,16) = 0

      Q2(10,1) = (11*15**(0.5D0)*MUH)/3200
      Q2(10,2) = -(5**(0.5D0)*(85*VT - 170*ABS(ETAH) +
     >   66*ABS(MUH)))/6400
      Q2(10,3) = -(37*3**(0.5D0)*MUH)/320
      Q2(10,4) = -(3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(18*ABS(ETAH) -
     >   9*VT + 22*ABS(MUH)))/6400
      Q2(10,5) = 0
      Q2(10,6) = -(17*15**(0.5D0)*ETAH)/320
      Q2(10,7) = 0
      Q2(10,8) = (27*35**(0.5D0)*ETAH)/1600
      Q2(10,9) = -(33*3**(0.5D0)*MUH)/3200
      Q2(10,10) = (51*VT)/1280 + (17*ABS(ETAH))/128 +
     >   (99*ABS(MUH))/3200
      Q2(10,11) = (111*15**(0.5D0)*MUH)/1600
      Q2(10,12) = -(3*3**(0.5D0)*7**(0.5D0)*(9*VT + 30*ABS(ETAH) -
     >   22*ABS(MUH)))/6400
      Q2(10,13) = 0
      Q2(10,14) = (17*35**(0.5D0)*ETAH)/640
      Q2(10,15) = 0
      Q2(10,16) = -(63*15**(0.5D0)*ETAH)/3200

      Q2(11,1) = (9*VT)/1280 - (9*ABS(ETAH))/640 -
     >   (9*ABS(MUH))/640
      Q2(11,2) = (9*3**(0.5D0)*MUH)/320
      Q2(11,3) = -(9*5**(0.5D0)*(3*VT - 6*ABS(ETAH) +
     >   10*ABS(MUH)))/6400
      Q2(11,4) = -(9*7**(0.5D0)*MUH)/640
      Q2(11,5) = (9*3**(0.5D0)*ETAH)/320
      Q2(11,6) = 0
      Q2(11,7) = -(27*15**(0.5D0)*ETAH)/1600
      Q2(11,8) = 0
      Q2(11,9) = -(9*5**(0.5D0)*(3*VT + 10*ABS(ETAH) -
     >   6*ABS(MUH)))/6400
      Q2(11,10) = -(27*15**(0.5D0)*MUH)/1600
      Q2(11,11) = (81*VT)/6400 + (27*ABS(ETAH))/640 +
     >   (27*ABS(MUH))/640
      Q2(11,12) = (27*35**(0.5D0)*MUH)/3200
      Q2(11,13) = -(9*7**(0.5D0)*ETAH)/640
      Q2(11,14) = 0
      Q2(11,15) = (27*35**(0.5D0)*ETAH)/3200
      Q2(11,16) = 0

      Q2(12,1) = (27*35**(0.5D0)*MUH)/22400
      Q2(12,2) = -(9*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*
     >   (14*ABS(ETAH) - 7*VT + 6*ABS(MUH)))/44800
      Q2(12,3) = (27*7**(0.5D0)*MUH)/560
      Q2(12,4) = -(27*5**(0.5D0)*(9*VT - 18*ABS(ETAH) +
     >   14*ABS(MUH)))/44800
      Q2(12,5) = 0
      Q2(12,6) = (27*35**(0.5D0)*ETAH)/1600
      Q2(12,7) = 0
      Q2(12,8) = -(243*15**(0.5D0)*ETAH)/11200
      Q2(12,9) = -(81*7**(0.5D0)*MUH)/22400
      Q2(12,10) = -(9*3**(0.5D0)*7**(0.5D0)*(21*VT +
     >   70*ABS(ETAH) - 18*ABS(MUH)))/44800
      Q2(12,11) = -(81*35**(0.5D0)*MUH)/2800
      Q2(12,12) = (729*VT)/44800 + (243*ABS(ETAH))/4480 +
     >   (81*ABS(MUH))/3200
      Q2(12,13) = 0
      Q2(12,14) = -(63*15**(0.5D0)*ETAH)/3200
      Q2(12,15) = 0
      Q2(12,16) = (243*35**(0.5D0)*ETAH)/22400

      Q2(13,1) = -(9*7**(0.5D0)*ETAH)/896
      Q2(13,2) = 0
      Q2(13,3) = (27*35**(0.5D0)*ETAH)/22400
      Q2(13,4) = 0
      Q2(13,5) = -(3*3**(0.5D0)*7**(0.5D0)*(35*VT -
     >   30*ABS(ETAH) + 14*ABS(MUH)))/8960
      Q2(13,6) = -(9*7**(0.5D0)*MUH)/160
      Q2(13,7) = -(3*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*
     >   (18*ABS(ETAH) - 21*VT + 70*ABS(MUH)))/44800
      Q2(13,8) = -(21*3**(0.5D0)*MUH)/640
      Q2(13,9) = -(9*35**(0.5D0)*ETAH)/112
      Q2(13,10) = 0
      Q2(13,11) = (27*7**(0.5D0)*ETAH)/560
      Q2(13,12) = 0
      Q2(13,13) = (81*VT)/1792 + (9*ABS(ETAH))/128 +
     >   (81*ABS(MUH))/4480
      Q2(13,14) = (81*3**(0.5D0)*MUH)/1120
      Q2(13,15) = -(27*5**(0.5D0)*(9*VT + 14*ABS(ETAH) -
     >   30*ABS(MUH)))/44800
      Q2(13,16) = (81*7**(0.5D0)*MUH)/4480

      Q2(14,1) = 0
      Q2(14,2) = -(51*7**(0.5D0)*ETAH)/4480
      Q2(14,3) = 0
      Q2(14,4) = (27*3**(0.5D0)*ETAH)/3200
      Q2(14,5) = (33*7**(0.5D0)*MUH)/3200
      Q2(14,6) = -(3**(0.5D0)*7**(0.5D0)*(595*VT -
     >   510*ABS(ETAH) + 462*ABS(MUH)))/44800
      Q2(14,7) = -(111*35**(0.5D0)*MUH)/1600
      Q2(14,8) = (189*VT)/6400 - (81*ABS(ETAH))/3200 -
     >   (231*ABS(MUH))/3200
      Q2(14,9) = 0
      Q2(14,10) = -(51*35**(0.5D0)*ETAH)/560
      Q2(14,11) = 0
      Q2(14,12) = (27*15**(0.5D0)*ETAH)/400
      Q2(14,13) = -(297*3**(0.5D0)*MUH)/22400
      Q2(14,14) = (459*VT)/8960 + (51*ABS(ETAH))/640 +
     >   (891*ABS(MUH))/22400
      Q2(14,15) = (999*15**(0.5D0)*MUH)/11200
      Q2(14,16) = -(27*3**(0.5D0)*7**(0.5D0)*(9*VT +
     >   14*ABS(ETAH) - 22*ABS(MUH)))/44800

      Q2(15,1) = (27*35**(0.5D0)*ETAH)/22400
      Q2(15,2) = 0
      Q2(15,3) = -(81*7**(0.5D0)*ETAH)/22400
      Q2(15,4) = 0
      Q2(15,5) = -(9*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*
     >   (6*ABS(ETAH) - 7*VT + 14*ABS(MUH)))/44800
      Q2(15,6) = (27*35**(0.5D0)*MUH)/1600
      Q2(15,7) = -(9*3**(0.5D0)*7**(0.5D0)*(21*VT -
     >   18*ABS(ETAH) + 70*ABS(MUH)))/44800
      Q2(15,8) = -(63*15**(0.5D0)*MUH)/3200
      Q2(15,9) = (27*7**(0.5D0)*ETAH)/560
      Q2(15,10) = 0
      Q2(15,11) = -(81*35**(0.5D0)*ETAH)/2800
      Q2(15,12) = 0
      Q2(15,13) = -(27*5**(0.5D0)*(9*VT + 14*ABS(ETAH) -
     >   18*ABS(MUH)))/44800
      Q2(15,14) = -(243*15**(0.5D0)*MUH)/11200
      Q2(15,15) = (729*VT)/44800 + (81*ABS(ETAH))/3200 +
     >   (243*ABS(MUH))/4480
      Q2(15,16) = (243*35**(0.5D0)*MUH)/22400

      Q2(16,1) = 0
      Q2(16,2) = (27*3**(0.5D0)*ETAH)/3200
      Q2(16,3) = 0
      Q2(16,4) = -(729*7**(0.5D0)*ETAH)/156800
      Q2(16,5) = (27*3**(0.5D0)*MUH)/3200
      Q2(16,6) = (189*VT)/6400 - (81*ABS(ETAH))/3200 -
     >   (81*ABS(MUH))/3200
      Q2(16,7) = (27*15**(0.5D0)*MUH)/400
      Q2(16,8) = -(27*3**(0.5D0)*7**(0.5D0)*(63*VT -
     >   54*ABS(ETAH) + 98*ABS(MUH)))/313600
      Q2(16,9) = 0
      Q2(16,10) = (27*15**(0.5D0)*ETAH)/400
      Q2(16,11) = 0
      Q2(16,12) = -(729*35**(0.5D0)*ETAH)/19600
      Q2(16,13) = -(729*7**(0.5D0)*MUH)/156800
      Q2(16,14) = -(27*3**(0.5D0)*7**(0.5D0)*(63*VT +
     >   98*ABS(ETAH) - 54*ABS(MUH)))/313600
      Q2(16,15) = -(729*35**(0.5D0)*MUH)/19600
      Q2(16,16) = (6561*VT)/313600 + (729*ABS(ETAH))/22400 +
     >   (729*ABS(MUH))/22400

      Q2(1,17) = (((2500*Q(01) + 180*Q(11) - 300*5**(0.5D0)*Q(03) -
     >   300*5**(0.5D0)*Q(09)))/25600)*VOL(I2,J2,J,I)
      Q2(2,17) = (((2125*Q(02) - 255*5**(0.5D0)*Q(10) -
     >   225*21**(0.5D0)*Q(04) + 27*105**(0.5D0)*Q(12)))/19200)*
     >   VOL(I2,J2,J,I)
      Q2(3,17) = (-(3*5**(0.5D0)*(500*Q(01) + 180*Q(11) -
     >   300*5**(0.5D0)*Q(03) - 60*5**(0.5D0)*Q(09)))/128000)*
     >   VOL(I2,J2,J,I)
      Q2(4,17) = ((3*(675*Q(04) - 81*5**(0.5D0)*Q(12) -
     >   175*21**(0.5D0)*Q(02) + 21*105**(0.5D0)*Q(10)))/44800)*
     >   VOL(I2,J2,J,I)
      Q2(5,17) = (((2125*Q(05) - 255*5**(0.5D0)*Q(07) -
     >   225*21**(0.5D0)*Q(13) + 27*105**(0.5D0)*Q(15)))/19200)*
     >   VOL(I2,J2,J,I)
      Q2(6,17) = (((7225*Q(06) + 1701*Q(16) - 765*21**(0.5D0)*Q(08) -
     >   765*21**(0.5D0)*Q(14)))/57600)*VOL(I2,J2,J,I)
      Q2(7,17) = (-(15**(0.5D0)*(8500*3**(0.5D0)*Q(05) -
     >   2700*7**(0.5D0)*Q(13) - 5100*15**(0.5D0)*Q(07) +
     >   1620*35**(0.5D0)*Q(15)))/1920000)*VOL(I2,J2,J,I)
      Q2(8,17) = (((2295*Q(08) + 1323*Q(14) - 595*21**(0.5D0)*Q(06) -
     >   243*21**(0.5D0)*Q(16)))/44800)*VOL(I2,J2,J,I)
      Q2(9,17) = (-(3*5**(0.5D0)*(500*Q(01) + 180*Q(11) -
     >   60*5**(0.5D0)*Q(03) - 300*5**(0.5D0)*Q(09)))/128000)*
     >   VOL(I2,J2,J,I)
      Q2(10,17) = (-(15**(0.5D0)*(8500*3**(0.5D0)*Q(02) -
     >   2700*7**(0.5D0)*Q(04) - 5100*15**(0.5D0)*Q(10) +
     >   1620*35**(0.5D0)*Q(12)))/1920000)*VOL(I2,J2,J,I)
      Q2(11,17) = ((9*(100*Q(01) + 180*Q(11) - 60*5**(0.5D0)*Q(03) -
     >   60*5**(0.5D0)*Q(09)))/128000)*VOL(I2,J2,J,I)
      Q2(12,17) = ((9*35**(0.5D0)*(4900*3**(0.5D0)*Q(02) -
     >   2700*7**(0.5D0)*Q(04) - 2940*15**(0.5D0)*Q(10) +
     >   1620*35**(0.5D0)*Q(12)))/31360000)*VOL(I2,J2,J,I)
      Q2(13,17) = ((3*(675*Q(13) - 81*5**(0.5D0)*Q(15) -
     >   175*21**(0.5D0)*Q(05) + 21*105**(0.5D0)*Q(07)))/44800)*
     >   VOL(I2,J2,J,I)
      Q2(14,17) = (((1323*Q(08) + 2295*Q(14) - 595*21**(0.5D0)*Q(06) -
     >   243*21**(0.5D0)*Q(16)))/44800)*VOL(I2,J2,J,I)
      Q2(15,17) = ((9*35**(0.5D0)*(4900*3**(0.5D0)*Q(05) -
     >   2700*7**(0.5D0)*Q(13) - 2940*15**(0.5D0)*Q(07) +
     >   1620*35**(0.5D0)*Q(15)))/31360000)*VOL(I2,J2,J,I)
      Q2(16,17) = ((9*(720300*Q(06) + 510300*Q(16) -
     >   132300*21**(0.5D0)*Q(08) - 132300*21**(0.5D0)*Q(14)))/
     >   219520000)*VOL(I2,J2,J,I)

      Q2(1,17) = Q2(1,17) +
     >   (XNI(4,J2)/320 + (21*XNI(3,J2))/1280 + (21*XNI(2,J2))/1280 +
     >   XNI(1,J2)/320)*MUH*SIGN(1.0,MUH) + (XNJ(1)/320 +
     >   (21*XNJ(2))/1280 + (21*XNJ(3))/1280 + XNJ(4)/320)*ETAH*
     >   SIGN(1.0,ETAH)
      Q2(2,17) = Q2(2,17)
     >   - (3**(0.5D0)*ETAH*SIGN(1.0,ETAH)*(173*XNJ(1) +
     >   756*XNJ(2) - 756*XNJ(3) - 173*XNJ(4)))/57600 - (3**(0.5D0)*
     >   MUH*(132*XNI(4,J2) + 693*XNI(3,J2) + 693*XNI(2,J2) +
     >   132*XNI(1,J2)))/57600
      Q2(3,17) = Q2(3,17) +
     >   (3*5**(0.5D0)*MUH*SIGN(1.0,MUH)*(4*XNI(4,J2) +
     >   21*XNI(3,J2) + 21*XNI(2,J2) + 4*XNI(1,J2)))/6400 +
     >   (3*5**(0.5D0)*ETAH*SIGN(1.0,ETAH)*(XNJ(1) - 6*XNJ(2) -
     >   6*XNJ(3) + XNJ(4)))/6400
      Q2(4,17) = Q2(4,17) +
     >   (3*7**(0.5D0)*ETAH*SIGN(1.0,ETAH)*(74*XNJ(1) +
     >   513*XNJ(2) - 513*XNJ(3) - 74*XNJ(4)))/313600 -
     >   (3*7**(0.5D0)*MUH*(84*XNI(4,J2) + 441*XNI(3,J2) +
     >   441*XNI(2,J2) + 84*XNI(1,J2)))/313600
      Q2(5,17) = Q2(5,17) +
     >   (3**(0.5D0)*MUH*SIGN(1.0,MUH)*(173*XNI(4,J2) +
     >   756*XNI(3,J2) - 756*XNI(2,J2) - 173*XNI(1,J2)))/57600 -
     >   (3**(0.5D0)*ETAH*(132*XNJ(1) + 693*XNJ(2) + 693*XNJ(3) +
     >   132*XNJ(4)))/57600
      Q2(6,17) = Q2(6,17) +
     >   ((231*XNI(2,J2))/8000 - (231*XNI(3,J2))/8000 -
     >   (1903*XNI(4,J2))/288000 + (1903*XNI(1,J2))/288000)*MUH +
     >   ((1903*XNJ(1))/288000 + (231*XNJ(2))/8000 -
     >   (231*XNJ(3))/8000 - (1903*XNJ(4))/288000)*ETAH
      Q2(7,17) = Q2(7,17) +
     >   (15**(0.5D0)*MUH*SIGN(1.0,MUH)*(173*XNI(4,J2) +
     >   756*XNI(3,J2) - 756*XNI(2,J2) - 173*XNI(1,J2)))/96000 -
     >   (15**(0.5D0)*ETAH*(33*XNJ(1) - 198*XNJ(2) - 198*XNJ(3) +
     >   33*XNJ(4)))/96000
      Q2(8,17) = Q2(8,17)
     >   - (21**(0.5D0)*ETAH*(814*XNJ(1) + 5643*XNJ(2) -
     >   5643*XNJ(3) - 814*XNJ(4)))/1568000 - (21**(0.5D0)*MUH
     >   *(1211*XNI(4,J2) + 5292*XNI(3,J2) - 5292*XNI(2,J2) -
     >   1211*XNI(1,J2)))/1568000
      Q2(9,17) = Q2(9,17) +
     >   (3*5**(0.5D0)*ETAH*SIGN(1.0,ETAH)*(4*XNJ(1) +
     >   21*XNJ(2) + 21*XNJ(3) + 4*XNJ(4)))/6400 + (3*5**(0.5D0)*
     >   MUH*SIGN(1.0,MUH)*(XNI(4,J2) - 6*XNI(3,J2) -
     >   6*XNI(2,J2) + XNI(1,J2)))/6400
      Q2(10,17) = Q2(10,17)
     >   - (15**(0.5D0)*ETAH*SIGN(1.0,ETAH)*(173*XNJ(1) +
     >   756*XNJ(2) - 756*XNJ(3) - 173*XNJ(4)))/96000 -
     >   (15**(0.5D0)*MUH*(33*XNI(4,J2) - 198*XNI(3,J2) -
     >   198*XNI(2,J2) + 33*XNI(1,J2)))/96000
      Q2(11,17) = Q2(11,17) +
     >   ((9*XNI(4,J2))/6400 - (27*XNI(3,J2))/3200 -
     >   (27*XNI(2,J2))/3200 + (9*XNI(1,J2))/6400)*MUH*
     >   SIGN(1.0,MUH) + ((9*XNJ(1))/6400 - (27*XNJ(2))/3200 -
     >   (27*XNJ(3))/3200 + (9*XNJ(4))/6400)*ETAH*SIGN(1.0,ETAH)
      Q2(12,17) = Q2(12,17) +
     >   (-(9*5**(0.5D0)*7**(0.5D0)*(21*XNI(4,J2) - 126*XNI(3,J2) -
     >   126*XNI(2,J2) + 21*XNI(1,J2)))/1568000)*MUH +
     >   ((9*5**(0.5D0)*7**(0.5D0)*(74*XNJ(1) + 513*XNJ(2) -
     >   513*XNJ(3) - 74*XNJ(4)))/1568000)*ETAH*SIGN(1.0,ETAH)
      Q2(13,17) = Q2(13,17)
     >   - (3*7**(0.5D0)*ETAH*(84*XNJ(1) + 441*XNJ(2) +
     >   441*XNJ(3) + 84*XNJ(4)))/313600 - (3*7**(0.5D0)*MUH*
     >   SIGN(1.0,MUH)*(74*XNI(4,J2) + 513*XNI(3,J2) -
     >   513*XNI(2,J2) - 74*XNI(1,J2)))/313600
      Q2(14,17) = Q2(14,17) +
     >   (21**(0.5D0)*ETAH*(1211*XNJ(1) + 5292*XNJ(2) -
     >   5292*XNJ(3) - 1211*XNJ(4)))/1568000 + (21**(0.5D0)*
     >   MUH*(814*XNI(4,J2) + 5643*XNI(3,J2) - 5643*XNI(2,J2) -
     >   814*XNI(1,J2)))/1568000
      Q2(15,17) = Q2(15,17) +
     >   (-(9*5**(0.5D0)*7**(0.5D0)*(74*XNI(4,J2) + 513*XNI(3,J2) -
     >   513*XNI(2,J2) - 74*XNI(1,J2)))/1568000)*MUH*
     >   SIGN(1.0,MUH) + (-(9*5**(0.5D0)*7**(0.5D0)*(21*XNJ(1) -
     >   126*XNJ(2) - 126*XNJ(3) + 21*XNJ(4)))/1568000)*ETAH
      Q2(16,17) = Q2(16,17) +
     >   ((999*XNI(4,J2))/784000 + (13851*XNI(3,J2))/1568000 -
     >   (13851*XNI(2,J2))/1568000 - (999*XNI(1,J2))/784000)*
     >   MUH + ((13851*XNJ(3))/1568000 - (13851*XNJ(2))/
     >   1568000 - (999*XNJ(1))/784000 + (999*XNJ(4))/784000)*ETAH
      ENDIF
*
      CALL ALSBD(IELEM**2,1,Q2,IER,IELEM**2)
      IF(IER.NE.0) CALL XABORT('SNFTH2: SINGULAR MATRIX.')
*
      IF(IELEM.EQ.2)THEN
      CORNERQ(1) = (Q2(01,IELEM**2+1) + 3*Q2(04,IELEM**2+1) -
     >   3**(0.5D0)*Q2(02,IELEM**2+1) - 3**(0.5D0)*Q2(03,IELEM**2+1))
      CORNERQ(2) = (Q2(01,IELEM**2+1) - 3*Q2(04,IELEM**2+1) +
     >   3**(0.5D0)*Q2(02,IELEM**2+1) - 3**(0.5D0)*Q2(03,IELEM**2+1))
      CORNERQ(3) = (Q2(01,IELEM**2+1) - 3*Q2(04,IELEM**2+1) -
     >   3**(0.5D0)*Q2(02,IELEM**2+1) + 3**(0.5D0)*Q2(03,IELEM**2+1))
      CORNERQ(4) = (Q2(01,IELEM**2+1) + 3*Q2(04,IELEM**2+1) +
     >   3**(0.5D0)*Q2(02,IELEM**2+1) + 3**(0.5D0)*Q2(03,IELEM**2+1))
      ELSEIF(IELEM.EQ.3)THEN
      CORNERQ(1) = (Q2(01,IELEM**2+1) + 3*Q2(05,IELEM**2+1) +
     >   5*Q2(09,IELEM**2+1) - 3**(0.5D0)*Q2(02,IELEM**2+1) -
     >   3**(0.5D0)*Q2(04,IELEM**2+1) + 5**(0.5D0)*Q2(03,IELEM**2+1) +
     >   5**(0.5D0)*Q2(07,IELEM**2+1) - 15**(0.5D0)*Q2(06,IELEM**2+1) -
     >   15**(0.5D0)*Q2(08,IELEM**2+1))
      CORNERQ(2) = (Q2(01,IELEM**2+1) - (5*Q2(09,IELEM**2+1))/2 -
     >   3**(0.5D0)*Q2(04,IELEM**2+1) -
     >   (5**(0.5D0)*Q2(03,IELEM**2+1))/2 +
     > 5**(0.5D0)*Q2(07,IELEM**2+1) +
     > (15**(0.5D0)*Q2(06,IELEM**2+1))/2)
      CORNERQ(3) = (Q2(01,IELEM**2+1) - 3*Q2(05,IELEM**2+1) +
     >   5*Q2(09,IELEM**2+1) + 3**(0.5D0)*Q2(02,IELEM**2+1) -
     >   3**(0.5D0)*Q2(04,IELEM**2+1) + 5**(0.5D0)*Q2(03,IELEM**2+1) +
     >   5**(0.5D0)*Q2(07,IELEM**2+1) - 15**(0.5D0)*Q2(06,IELEM**2+1) +
     >   15**(0.5D0)*Q2(08,IELEM**2+1))
      CORNERQ(4) = (Q2(01,IELEM**2+1) - (5*Q2(09,IELEM**2+1))/2 -
     >   3**(0.5D0)*Q2(02,IELEM**2+1) + 5**(0.5D0)*Q2(03,IELEM**2+1) -
     >   (5**(0.5D0)*Q2(07,IELEM**2+1))/2 +
     >   (15**(0.5D0)*Q2(08,IELEM**2+1))/2)
      CORNERQ(5) = (Q2(01,IELEM**2+1) + (5*Q2(09,IELEM**2+1))/4 -
     >   (5**(0.5D0)*Q2(03,IELEM**2+1))/2 -
     >   (5**(0.5D0)*Q2(07,IELEM**2+1))/2)
      CORNERQ(6) = (Q2(01,IELEM**2+1) - (5*Q2(09,IELEM**2+1))/2 +
     >   3**(0.5D0)*Q2(02,IELEM**2+1) + 5**(0.5D0)*Q2(03,IELEM**2+1) -
     >   (5**(0.5D0)*Q2(07,IELEM**2+1))/2 -
     >   (15**(0.5D0)*Q2(08,IELEM**2+1))/2)
      CORNERQ(7) = (Q2(01,IELEM**2+1) - 3*Q2(05,IELEM**2+1) +
     >   5*Q2(09,IELEM**2+1) - 3**(0.5D0)*Q2(02,IELEM**2+1) +
     >   3**(0.5D0)*Q2(04,IELEM**2+1) + 5**(0.5D0)*Q2(03,IELEM**2+1) +
     >   5**(0.5D0)*Q2(07,IELEM**2+1) + 15**(0.5D0)*Q2(06,IELEM**2+1) -
     >   15**(0.5D0)*Q2(08,IELEM**2+1))
      CORNERQ(8) = (Q2(01,IELEM**2+1) - (5*Q2(09,IELEM**2+1))/2 +
     >   3**(0.5D0)*Q2(04,IELEM**2+1) -
     >   (5**(0.5D0)*Q2(03,IELEM**2+1))/2 +
     >   5**(0.5D0)*Q2(07,IELEM**2+1) -
     >   (15**(0.5D0)*Q2(06,IELEM**2+1))/2)
      CORNERQ(9) = (Q2(01,IELEM**2+1) + 3*Q2(05,IELEM**2+1) +
     >   5*Q2(09,IELEM**2+1) + 3**(0.5D0)*Q2(02,IELEM**2+1) +
     >   3**(0.5D0)*Q2(04,IELEM**2+1) + 5**(0.5D0)*Q2(03,IELEM**2+1) +
     >   5**(0.5D0)*Q2(07,IELEM**2+1) + 15**(0.5D0)*Q2(06,IELEM**2+1) +
     >   15**(0.5D0)*Q2(08,IELEM**2+1))
      ELSEIF(IELEM.EQ.4)THEN
      CORNERQ(01) = (Q2(01,IELEM**2+1) + 3*Q2(06,IELEM**2+1) +
     >   5*Q2(11,IELEM**2+1) + 7*Q2(16,IELEM**2+1) -
     >   3**(0.5D0)*Q2(02,IELEM**2+1) - 3**(0.5D0)*Q2(05,IELEM**2+1) +
     >   5**(0.5D0)*Q2(03,IELEM**2+1) - 7**(0.5D0)*Q2(04,IELEM**2+1) +
     > 5**(0.5D0)*Q2(09,IELEM**2+1) - 7**(0.5D0)*Q2(13,IELEM**2+1) -
     > 15**(0.5D0)*Q2(07,IELEM**2+1) - 15**(0.5D0)*Q2(10,IELEM**2+1) +
     > 21**(0.5D0)*Q2(08,IELEM**2+1) + 21**(0.5D0)*Q2(14,IELEM**2+1) -
     > 35**(0.5D0)*Q2(12,IELEM**2+1) - 35**(0.5D0)*Q2(15,IELEM**2+1))
      CORNERQ(02) = (Q2(01,IELEM**2+1) + Q2(06,IELEM**2+1) -
     >   (5*Q2(11,IELEM**2+1))/3 - (77*Q2(16,IELEM**2+1))/27 -
     >   (3**(0.5D0)*Q2(02,IELEM**2+1))/3 - 3**(0.5D0)*
     >   Q2(05,IELEM**2+1) - (5**(0.5D0)*Q2(03,IELEM**2+1))/3 +
     >   (11*7**(0.5D0)*Q2(04,IELEM**2+1))/27 + 5**(0.5D0)*
     >   Q2(09,IELEM**2+1) - 7**(0.5D0)*Q2(13,IELEM**2+1) +
     >   (15**(0.5D0)*Q2(07,IELEM**2+1))/3 - (15**(0.5D0)*
     >   Q2(10,IELEM**2+1))/3 - (11*21**(0.5D0)*
     >   Q2(08,IELEM**2+1))/27 + (21**(0.5D0)*Q2(14,IELEM**2+1))/3 +
     >   (11*35**(0.5D0)*Q2(12,IELEM**2+1))/27 +
     >   (35**(0.5D0)*Q2(15,IELEM**2+1))/3)
      CORNERQ(03) = (Q2(01,IELEM**2+1) - Q2(06,IELEM**2+1) -
     >   (5*Q2(11,IELEM**2+1))/3 + (77*Q2(16,IELEM**2+1))/27 +
     >   (3**(0.5D0)*Q2(02,IELEM**2+1))/3 - 3**(0.5D0)*
     >   Q2(05,IELEM**2+1) - (5**(0.5D0)*Q2(03,IELEM**2+1))/3 -
     >   (11*7**(0.5D0)*Q2(04,IELEM**2+1))/27 + 5**(0.5D0)*
     >   Q2(09,IELEM**2+1) - 7**(0.5D0)*Q2(13,IELEM**2+1) +
     >   (15**(0.5D0)*Q2(07,IELEM**2+1))/3 + (15**(0.5D0)*
     >   Q2(10,IELEM**2+1))/3 + (11*21**(0.5D0)*
     >   Q2(08,IELEM**2+1))/27 - (21**(0.5D0)*Q2(14,IELEM**2+1))/3 -
     >   (11*35**(0.5D0)*Q2(12,IELEM**2+1))/27 + (35**(0.5D0)*
     >   Q2(15,IELEM**2+1))/3)
      CORNERQ(04) = (Q2(01,IELEM**2+1) - 3*Q2(06,IELEM**2+1) +
     >   5*Q2(11,IELEM**2+1) - 7*Q2(16,IELEM**2+1) +
     >   3**(0.5D0)*Q2(02,IELEM**2+1) - 3**(0.5D0)*Q2(05,IELEM**2+1) +
     >   5**(0.5D0)*Q2(03,IELEM**2+1) + 7**(0.5D0)*Q2(04,IELEM**2+1) +
     >   5**(0.5D0)*Q2(09,IELEM**2+1) - 7**(0.5D0)*Q2(13,IELEM**2+1) -
     >   15**(0.5D0)*Q2(07,IELEM**2+1) + 15**(0.5D0)*
     >   Q2(10,IELEM**2+1) - 21**(0.5D0)*Q2(08,IELEM**2+1) -
     >   21**(0.5D0)*Q2(14,IELEM**2+1) + 35**(0.5D0)*
     >   Q2(12,IELEM**2+1) - 35**(0.5D0)*Q2(15,IELEM**2+1))
      CORNERQ(05) = (Q2(01,IELEM**2+1) + Q2(06,IELEM**2+1) -
     >   (5*Q2(11,IELEM**2+1))/3 - (77*Q2(16,IELEM**2+1))/27 -
     >   3**(0.5D0)*Q2(02,IELEM**2+1) - (3**(0.5D0)*
     >   Q2(05,IELEM**2+1))/3 + 5**(0.5D0)*Q2(03,IELEM**2+1) -
     >   7**(0.5D0)*Q2(04,IELEM**2+1) - (5**(0.5D0)*
     >   Q2(09,IELEM**2+1))/3 + (11*7**(0.5D0)*
     >   Q2(13,IELEM**2+1))/27 - (15**(0.5D0)*Q2(07,IELEM**2+1))/3 +
     >   (15**(0.5D0)*Q2(10,IELEM**2+1))/3 + (21**(0.5D0)*
     >   Q2(08,IELEM**2+1))/3 - (11*21**(0.5D0)*
     >   Q2(14,IELEM**2+1))/27 + (35**(0.5D0)*Q2(12,IELEM**2+1))/3 +
     >   (11*35**(0.5D0)*Q2(15,IELEM**2+1))/27)
      CORNERQ(06) = (Q2(01,IELEM**2+1) + Q2(06,IELEM**2+1)/3 +
     >   (5*Q2(11,IELEM**2+1))/9 + (847*Q2(16,IELEM**2+1))/729 -
     >   (3**(0.5D0)*Q2(02,IELEM**2+1))/3 - (3**(0.5D0)*
     >   Q2(05,IELEM**2+1))/3 - (5**(0.5D0)*Q2(03,IELEM**2+1))/3 +
     >   (11*7**(0.5D0)*Q2(04,IELEM**2+1))/27 - (5**(0.5D0)*
     >   Q2(09,IELEM**2+1))/3 + (11*7**(0.5D0)*
     >   Q2(13,IELEM**2+1))/27 + (15**(0.5D0)*
     >   Q2(07,IELEM**2+1))/9 + (15**(0.5D0)*
     >   Q2(10,IELEM**2+1))/9 - (11*21**(0.5D0)*
     >   Q2(08,IELEM**2+1))/81 - (11*21**(0.5D0)*
     >   Q2(14,IELEM**2+1))/81 - (11*35**(0.5D0)*
     >   Q2(12,IELEM**2+1))/81 - (11*35**(0.5D0)*
     >   Q2(15,IELEM**2+1))/81)
      CORNERQ(07) = (Q2(01,IELEM**2+1) - Q2(06,IELEM**2+1)/3 +
     >   (5*Q2(11,IELEM**2+1))/9 - (847*Q2(16,IELEM**2+1))/729 +
     >   (3**(0.5D0)*Q2(02,IELEM**2+1))/3 -
     >   (3**(0.5D0)*Q2(05,IELEM**2+1))/3 -
     >   (5**(0.5D0)*Q2(03,IELEM**2+1))/3 -
     >   (11*7**(0.5D0)*Q2(04,IELEM**2+1))/27 -
     >   (5**(0.5D0)*Q2(09,IELEM**2+1))/3 +
     >   (11*7**(0.5D0)*Q2(13,IELEM**2+1))/27 +
     >   (15**(0.5D0)*Q2(07,IELEM**2+1))/9 -
     >   (15**(0.5D0)*Q2(10,IELEM**2+1))/9 +
     >   (11*21**(0.5D0)*Q2(08,IELEM**2+1))/81 +
     >   (11*21**(0.5D0)*Q2(14,IELEM**2+1))/81 +
     >   (11*35**(0.5D0)*Q2(12,IELEM**2+1))/81 -
     >   (11*35**(0.5D0)*Q2(15,IELEM**2+1))/81)
      CORNERQ(08) = (Q2(01,IELEM**2+1) - Q2(06,IELEM**2+1) -
     >   (5*Q2(11,IELEM**2+1))/3 + (77*Q2(16,IELEM**2+1))/27 +
     >   3**(0.5D0)*Q2(02,IELEM**2+1) - (3**(0.5D0)*
     >   Q2(05,IELEM**2+1))/3 + 5**(0.5D0)*Q2(03,IELEM**2+1) +
     >   7**(0.5D0)*Q2(04,IELEM**2+1) - (5**(0.5D0)*
     >   Q2(09,IELEM**2+1))/3 + (11*7**(0.5D0)*
     >   Q2(13,IELEM**2+1))/27 - (15**(0.5D0)*Q2(07,IELEM**2+1))/3 -
     >   (15**(0.5D0)*Q2(10,IELEM**2+1))/3 - (21**(0.5D0)*
     >   Q2(08,IELEM**2+1))/3 + (11*21**(0.5D0)*
     >   Q2(14,IELEM**2+1))/27 - (35**(0.5D0)*Q2(12,IELEM**2+1))/3 +
     >   (11*35**(0.5D0)*Q2(15,IELEM**2+1))/27)
      CORNERQ(09) = (Q2(01,IELEM**2+1) - Q2(06,IELEM**2+1) -
     >   (5*Q2(11,IELEM**2+1))/3 + (77*Q2(16,IELEM**2+1))/27 -
     >   3**(0.5D0)*Q2(02,IELEM**2+1) + (3**(0.5D0)*
     >   Q2(05,IELEM**2+1))/3 + 5**(0.5D0)*Q2(03,IELEM**2+1) -
     >   7**(0.5D0)*Q2(04,IELEM**2+1) - (5**(0.5D0)*
     >   Q2(09,IELEM**2+1))/3 - (11*7**(0.5D0)*
     >   Q2(13,IELEM**2+1))/27 + (15**(0.5D0)*Q2(07,IELEM**2+1))/3 +
     >   (15**(0.5D0)*Q2(10,IELEM**2+1))/3 - (21**(0.5D0)*
     > Q2(08,IELEM**2+1))/3 + (11*21**(0.5D0)*
     >   Q2(14,IELEM**2+1))/27 + (35**(0.5D0)*
     >   Q2(12,IELEM**2+1))/3 - (11*35**(0.5D0)*Q2(15,IELEM**2+1))/27)
      CORNERQ(10) = (Q2(01,IELEM**2+1) - Q2(06,IELEM**2+1)/3 +
     >   (5*Q2(11,IELEM**2+1))/9 - (847*Q2(16,IELEM**2+1))/729 -
     >   (3**(0.5D0)*Q2(02,IELEM**2+1))/3 + (3**(0.5D0)*
     >   Q2(05,IELEM**2+1))/3 - (5**(0.5D0)*Q2(03,IELEM**2+1))/3 +
     >   (11*7**(0.5D0)*Q2(04,IELEM**2+1))/27 - (5**(0.5D0)*
     >   Q2(09,IELEM**2+1))/3 - (11*7**(0.5D0)*
     >   Q2(13,IELEM**2+1))/27 - (15**(0.5D0)*
     >   Q2(07,IELEM**2+1))/9 + (15**(0.5D0)*Q2(10,IELEM**2+1))/9 +
     >   (11*21**(0.5D0)*Q2(08,IELEM**2+1))/81 + (11*21**(0.5D0)*
     >   Q2(14,IELEM**2+1))/81 - (11*35**(0.5D0)*
     >   Q2(12,IELEM**2+1))/81 + (11*35**(0.5D0)*
     >   Q2(15,IELEM**2+1))/81)
      CORNERQ(11) = (Q2(01,IELEM**2+1) + Q2(06,IELEM**2+1)/3 +
     >   (5*Q2(11,IELEM**2+1))/9 + (847*Q2(16,IELEM**2+1))/729 +
     >   (3**(0.5D0)*Q2(02,IELEM**2+1))/3 + (3**(0.5D0)*
     >   Q2(05,IELEM**2+1))/3 - (5**(0.5D0)*Q2(03,IELEM**2+1))/3 -
     >   (11*7**(0.5D0)*Q2(04,IELEM**2+1))/27 - (5**(0.5D0)*
     >   Q2(09,IELEM**2+1))/3 - (11*7**(0.5D0)*
     >   Q2(13,IELEM**2+1))/27 - (15**(0.5D0)*Q2(07,IELEM**2+1))/9 -
     >   (15**(0.5D0)*Q2(10,IELEM**2+1))/9 - (11*21**(0.5D0)*
     >   Q2(08,IELEM**2+1))/81 - (11*21**(0.5D0)*
     >   Q2(14,IELEM**2+1))/81 + (11*35**(0.5D0)*
     >   Q2(12,IELEM**2+1))/81 + (11*35**(0.5D0)*
     >   Q2(15,IELEM**2+1))/81)
      CORNERQ(12) = (Q2(01,IELEM**2+1) + Q2(06,IELEM**2+1) -
     >   (5*Q2(11,IELEM**2+1))/3 - (77*Q2(16,IELEM**2+1))/27 +
     >   3**(0.5D0)*Q2(02,IELEM**2+1) + (3**(0.5D0)*
     >   Q2(05,IELEM**2+1))/3 + 5**(0.5D0)*Q2(03,IELEM**2+1) +
     >   7**(0.5D0)*Q2(04,IELEM**2+1) - (5**(0.5D0)*
     >   Q2(09,IELEM**2+1))/3 - (11*7**(0.5D0)*
     >   Q2(13,IELEM**2+1))/27 + (15**(0.5D0)*
     >   Q2(07,IELEM**2+1))/3 - (15**(0.5D0)*Q2(10,IELEM**2+1))/3 +
     >   (21**(0.5D0)*Q2(08,IELEM**2+1))/3 -
     >   (11*21**(0.5D0)*Q2(14,IELEM**2+1))/27 -
     >   (35**(0.5D0)*Q2(12,IELEM**2+1))/3 -
     >   (11*35**(0.5D0)*Q2(15,IELEM**2+1))/27)
      CORNERQ(13) = (Q2(01,IELEM**2+1) - 3*Q2(06,IELEM**2+1) +
     >   5*Q2(11,IELEM**2+1) - 7*Q2(16,IELEM**2+1) -
     >   3**(0.5D0)*Q2(02,IELEM**2+1) + 3**(0.5D0)*Q2(05,IELEM**2+1) +
     >   5**(0.5D0)*Q2(03,IELEM**2+1) - 7**(0.5D0)*Q2(04,IELEM**2+1) +
     >   5**(0.5D0)*Q2(09,IELEM**2+1) + 7**(0.5D0)*Q2(13,IELEM**2+1) +
     >   15**(0.5D0)*Q2(07,IELEM**2+1) - 15**(0.5D0)*
     >   Q2(10,IELEM**2+1) - 21**(0.5D0)*Q2(08,IELEM**2+1) -
     >   21**(0.5D0)*Q2(14,IELEM**2+1) - 35**(0.5D0)*
     >   Q2(12,IELEM**2+1) + 35**(0.5D0)*Q2(15,IELEM**2+1))
      CORNERQ(14) = (Q2(01,IELEM**2+1) - Q2(06,IELEM**2+1) -
     >   (5*Q2(11,IELEM**2+1))/3 + (77*Q2(16,IELEM**2+1))/27 -
     >   (3**(0.5D0)*Q2(02,IELEM**2+1))/3 + 3**(0.5D0)*
     >   Q2(05,IELEM**2+1) - (5**(0.5D0)*Q2(03,IELEM**2+1))/3 +
     >   (11*7**(0.5D0)*Q2(04,IELEM**2+1))/27 + 5**(0.5D0)*
     >   Q2(09,IELEM**2+1) + 7**(0.5D0)*Q2(13,IELEM**2+1) -
     >   (15**(0.5D0)*Q2(07,IELEM**2+1))/3 - (15**(0.5D0)*
     >   Q2(10,IELEM**2+1))/3 + (11*21**(0.5D0)*
     >   Q2(08,IELEM**2+1))/27 - (21**(0.5D0)*
     >   Q2(14,IELEM**2+1))/3 + (11*35**(0.5D0)*
     >   Q2(12,IELEM**2+1))/27 - (35**(0.5D0)*Q2(15,IELEM**2+1))/3)
      CORNERQ(15) = (Q2(01,IELEM**2+1) + Q2(06,IELEM**2+1) -
     >   (5*Q2(11,IELEM**2+1))/3 - (77*Q2(16,IELEM**2+1))/27 +
     >   (3**(0.5D0)*Q2(02,IELEM**2+1))/3 + 3**(0.5D0)*
     >   Q2(05,IELEM**2+1) - (5**(0.5D0)*Q2(03,IELEM**2+1))/3 -
     >   (11*7**(0.5D0)*Q2(04,IELEM**2+1))/27 + 5**(0.5D0)*
     >   Q2(09,IELEM**2+1) + 7**(0.5D0)*Q2(13,IELEM**2+1) -
     >   (15**(0.5D0)*Q2(07,IELEM**2+1))/3 + (15**(0.5D0)*
     >   Q2(10,IELEM**2+1))/3 - (11*21**(0.5D0)*
     >   Q2(08,IELEM**2+1))/27 + (21**(0.5D0)*
     >   Q2(14,IELEM**2+1))/3 - (11*35**(0.5D0)*
     >   Q2(12,IELEM**2+1))/27 - (35**(0.5D0)*Q2(15,IELEM**2+1))/3)
      CORNERQ(16) = (Q2(01,IELEM**2+1) + 3*Q2(06,IELEM**2+1) +
     >   5*Q2(11,IELEM**2+1) + 7*Q2(16,IELEM**2+1) +
     >   3**(0.5D0)*Q2(02,IELEM**2+1) + 3**(0.5D0)*
     >   Q2(05,IELEM**2+1) + 5**(0.5D0)*Q2(03,IELEM**2+1) +
     >   7**(0.5D0)*Q2(04,IELEM**2+1) + 5**(0.5D0)*Q2(09,IELEM**2+1) +
     >   7**(0.5D0)*Q2(13,IELEM**2+1) + 15**(0.5D0)*Q2(07,IELEM**2+1) +
     >   15**(0.5D0)*Q2(10,IELEM**2+1) + 21**(0.5D0)*
     >   Q2(08,IELEM**2+1) + 21**(0.5D0)*Q2(14,IELEM**2+1) +
     >   35**(0.5D0)*Q2(12,IELEM**2+1) + 35**(0.5D0)*Q2(15,IELEM**2+1))
      ENDIF
*
      IF(JL.LT.ISPLH)THEN
      DO IEL=1,IELEM
         IF(IFROMJ.EQ.2)THEN
            XNJ(IEL) = CORNERQ(IELEM**2 - IELEM +IEL)
         ELSEIF(IFROMJ.EQ.4)THEN
            XNJ(IEL) = CORNERQ(IEL)
         ENDIF
      ENDDO
      ELSEIF((JL.EQ.ISPLH).AND.(IHEXJ.LE.NHEX))THEN
      DO INDEXE=1,IELEM
         IEL=INDEXE
         IELJ=INDEXE
         I3=I2
         IF((J.EQ.1).AND.(ILOZJ.EQ.3))THEN
            IELJ=IELEM - (INDEXE-1)
            I3=ISPLH+1 -I2
         ENDIF
*
         IF(IFROMJ.EQ.2)THEN
            BFLUX(ISIDEJ,INDEXJ,I3,IELJ) = CORNERQ(IELEM**2-IELEM+IEL)
         ELSEIF(IFROMJ.EQ.4)THEN
            BFLUX(ISIDEJ,INDEXJ,I3,IELJ) = CORNERQ(IEL)
         ENDIF
      ENDDO
      ENDIF
*
      IF(IL.LT.ISPLH)THEN
      DO IEL=1,IELEM
         IF(IFROMI.EQ.1)THEN
            XNI(IEL,J2) = CORNERQ(IEL*IELEM)
         ELSEIF(IFROMI.EQ.3)THEN
            XNI(IEL,J2) = CORNERQ((IELEM*(IEL-1))+1)
         ENDIF
      ENDDO
      ELSEIF((IL.EQ.ISPLH).AND.(IHEXI.LE.NHEX))THEN
      DO INDEXE=1,IELEM
         IEL=INDEXE
         IELI=INDEXE
         J3=J2
         IF((J.EQ.3).AND.(ILOZI.EQ.1))THEN
            IELI=IELEM - (INDEXE-1)
            J3=ISPLH+1-J2
         ENDIF
*
         IF(IFROMI.EQ.1)THEN
            BFLUX(ISIDEI,INDEXI,J3,IELI) = CORNERQ(IEL*IELEM)
         ELSEIF(IFROMI.EQ.3)THEN
            BFLUX(ISIDEI,INDEXI,J3,IELI) = CORNERQ((IELEM*(IEL-1))+1)
         ENDIF
      ENDDO
      ENDIF
*
      ILOZ=((J-1)*ISPLH+(J2-1))*ISPLH+I2
      DO 150 K=1,NSCT
      DO 140 IEL=1,IELEM**2
      FLUX(IEL,K,ILOZ,I) = FLUX(IEL,K,ILOZ,I) +
     1   2.0D0*W(M)*Q2(IEL,IELEM**2+1)*PL(K,M)
  140 CONTINUE
  150 CONTINUE
*--
  160 CONTINUE
  170 CONTINUE
  180 CONTINUE
  190 CONTINUE
      FLUX_G(:,:,:,:,IGEFF)=FLUX_G(:,:,:,:,IGEFF)+FLUX(:,:,:,:)
  200 CONTINUE
  210 CONTINUE
*$OMP END PARALLEL DO
  220 CONTINUE
      DO 230 IGEFF=1,NGEFF
        IF(.NOT.INCONV(IGEFF)) GO TO 230
        IG=INGIND(IGEFF)
        FUNKNO(:L4,IG)=
     1  RESHAPE(REAL(FLUX_G(:IELEM**2,:NSCT,:3*ISPLH**2,:NHEX,IGEFF)),
     2  (/ L4 /) )
  230 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(FLUX_G,FLUX,BFLUX,INDANG)
      RETURN
  400 FORMAT(16H SNFTH2: thread=,I8,12H --->(group=,I4,7H angle=,I4,1H))
      END
