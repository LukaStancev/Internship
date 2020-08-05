*DECK SNFT22
      SUBROUTINE SNFT22(NUN,NGRP,NGEFF,IMPX,INCONV,INGIND,LX,LY,IELEM,
     1 NMAT,NPQ,NSCT,MAT,VOL,TOTAL,NCODE,ZCODE,QEXT,LFIXUP,DU,DE,W,MRM,
     2 MRMY,DB,DA,PL,FUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* perform one inner iteration for solving SN equations in 2D Cartesian
* for the DISCONTINUOUS GALERKIN method. Energy-angle multithreading.
* Albedo boundary conditions.
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
* LX      number of meshes along X axis.
* LY      number of meshes along Y axis.
* IELEM   measure of order of the spatial approximation polynomial:
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
* FUNKNO  Legendre components of the flux and boundary fluxes.
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
      INTEGER NUN,NGRP,NGEFF,IMPX,INGIND(NGEFF),LX,LY,IELEM,NMAT,NPQ,
     1 NSCT,MAT(LX,LY),NCODE(4),MRM(NPQ),MRMY(NPQ)
      LOGICAL INCONV(NGEFF)
      REAL VOL(LX,LY),TOTAL(0:NMAT,NGRP),ZCODE(4),QEXT(NUN,NGRP),
     1 DU(NPQ),DE(NPQ),W(NPQ),DB(LX,NPQ),DA(LX,LY,NPQ),PL(NSCT,NPQ),
     2 FUNKNO(NUN,NGRP)
      LOGICAL LFIXUP
*----
*  LOCAL VARIABLES
*----
      INTEGER NPQD(4),IIND(4)
      DOUBLE PRECISION Q(IELEM**2),Q2(IELEM**2,(IELEM**2)+1),XNJ(IELEM),
     1 VT,CONST0,CONST1,CONST2,CONST3,CONST4,CONST5,CORNERQ(IELEM**2)
      PARAMETER(IUNOUT=6,RLOG=1.0E-8,PI=3.141592654)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: INDANG
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: FLUX
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: FLUX_G
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: XNI
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(INDANG(NPQ,4))
      ALLOCATE(XNI(IELEM,LY),FLUX(IELEM**2,NSCT,LX,LY))
      ALLOCATE(FLUX_G(IELEM**2,NSCT,LX,LY,NGEFF))
*----
*  DEFINITION OF CONSTANTS.
*----
      L4=IELEM*IELEM*LX*LY*NSCT
      CONST0=1.0D0/12.0D0
      CONST1=1.0D0/36.0D0
      CONST2=1.0D0/180.0D0
      CONST3=1.0D0/900.0D0
      CONST4=1.0D0/(44800.0D0*3.0D0)
      CONST5=1.0D0/2822400.0D0
*----
*  PARAMETER VALIDATION.
*----
      IF((IELEM.LT.1).OR.(IELEM.GT.4)) CALL XABORT('SNFT22: INVALID IE'
     1 //'LEM (DIAM) VALUE. CHECK INPUT DATA FILE.')
      FLUX_G(:IELEM**2,:NSCT,:LX,:LY,:NGEFF)=0.0D0
*----
*  SET OCTANT SWAPPING ORDER.
*----
      NPQD(:4)=0
      INDANG(:NPQ,:4)=0
      DO M=1,NPQ
        VU=DU(M)
        VE=DE(M)
        IF((VU.GE.0.0).AND.(VE.GE.0.0)) THEN
          IND=1
          JND=4
        ELSE IF((VU.LE.0.0).AND.(VE.GE.0.0)) THEN
          IND=2
          JND=3
        ELSE IF((VU.LE.0.0).AND.(VE.LE.0.0)) THEN
          IND=3
          JND=1
        ELSE
          IND=4
          JND=2
        ENDIF
        IIND(JND)=IND
        NPQD(IND)=NPQD(IND)+1
        INDANG(NPQD(IND),IND)=M
      ENDDO
*----
*  MAIN LOOP OVER OCTANTS.
*----
      DO 190 JND=1,4
      IND=IIND(JND)
*----
*  PRELIMINARY LOOPS FOR SETTING BOUNDARY CONDITIONS.
*----
*$OMP  PARALLEL DO
*$OMP1 PRIVATE(M,IG,WEIGHT,M1,E1,IOF,JOF,IEL,I,J)
*$OMP2 SHARED(FUNKNO) COLLAPSE(2)
      DO 70 IGEFF=1,NGEFF
      DO 60 IPQD=1,NPQD(IND)
      IF(.NOT.INCONV(IGEFF)) GO TO 60
      M=INDANG(IPQD,IND)
      IG=INGIND(IGEFF)
      WEIGHT=W(M)
      IF((NCODE(1).NE.4).AND.(WEIGHT.EQ.0.0)) THEN
        M1=MRM(M)
        DO 15 IEL=1,IELEM
        DO 10 J=1,LY
        IOF=((M-1)*LY+(J-1))*IELEM+IEL
        JOF=((M1-1)*LY+(J-1))*IELEM+IEL
        FUNKNO(L4+IOF,IG)=FUNKNO(L4+JOF,IG)
   10   CONTINUE
   15   CONTINUE
      ELSE IF((NCODE(1).NE.4).AND.(DU(M).GT.0.0)) THEN
        M1=MRM(M)
        DO 30 IEL=1,IELEM
        DO 20 J=1,LY
        IOF=((M-1)*LY+(J-1))*IELEM+IEL
        JOF=((M1-1)*LY+(J-1))*IELEM+IEL
        E1=FUNKNO(L4+IOF,IG)
        FUNKNO(L4+IOF,IG)=FUNKNO(L4+JOF,IG)
        FUNKNO(L4+JOF,IG)=E1
   20   CONTINUE
   30   CONTINUE
      ENDIF
      IF((NCODE(3).NE.4).AND.(DE(M).GT.0)) THEN
        M1=MRMY(M)
        DO 50 IEL=1,IELEM
        DO 40 I=1,LX
        IOF=((M-1)*LX+(I-1))*IELEM+IEL
        JOF=((M1-1)*LX+(I-1))*IELEM+IEL
        E1=FUNKNO(L4+IELEM*LY*NPQ+IOF,IG)
        FUNKNO(L4+IELEM*LY*NPQ+IOF,IG)=FUNKNO(L4+IELEM*LY*NPQ+JOF,IG)
        FUNKNO(L4+IELEM*LY*NPQ+JOF,IG)=E1
   40   CONTINUE
   50   CONTINUE
      ENDIF
   60 CONTINUE
   70 CONTINUE
*$OMP END PARALLEL DO
*----
*  MAIN SWAPPING LOOPS FOR SN FLUX CALCULATION
*----
*$OMP  PARALLEL DO
*$OMP1 PRIVATE(ITID,FLUX,M,IG,XNI,XNJ,Q,Q2,CORNERQ,IOF,IER,II,JJ,IEL)
*$OMP2 PRIVATE(JEL,I,J,K,VT) SHARED(FUNKNO) REDUCTION(+:FLUX_G)
*$OMP3 COLLAPSE(2)
      DO 180 IGEFF=1,NGEFF
      DO 170 IPQD=1,NPQD(IND)
#if defined(_OPENMP)
        ITID=omp_get_thread_num()
#else
        ITID=0
#endif
      IF(IMPX.GT.5) WRITE(IUNOUT,400) ITID,IGEFF,IPQD
      IF(.NOT.INCONV(IGEFF)) GO TO 170
      M=INDANG(IPQD,IND)
      IG=INGIND(IGEFF)
      FLUX(:IELEM**2,:NSCT,:LX,:LY)=0.0D0
      IF(W(M).EQ.0.0) GO TO 170
*----
*  LOOP OVER X- AND Y-DIRECTED AXES.
*----
      DO 155 II=1,LX
      I=II
      IF((IND.EQ.2).OR.(IND.EQ.3)) I=LX+1-I
      DO 100 IEL=1,IELEM
      IOF=(M-1)*IELEM*LX+(I-1)*IELEM+IEL
      IF((IND.EQ.1).OR.(IND.EQ.2)) THEN
         XNJ(IEL)=FUNKNO(L4+IELEM*LY*NPQ+IOF,IG)*ZCODE(3)
      ELSE
         XNJ(IEL)=FUNKNO(L4+IELEM*LY*NPQ+IOF,IG)*ZCODE(4)
      ENDIF
  100 CONTINUE
      DO 140 JJ=1,LY
      J=JJ
      IF((IND.EQ.3).OR.(IND.EQ.4)) J=LY+1-J
      DO 105 IEL=1,IELEM
      IF(II.EQ.1) THEN
         IOF=(M-1)*IELEM*LY+(J-1)*IELEM+IEL
         IF((IND.EQ.1).OR.(IND.EQ.4)) THEN
            XNI(IEL,J)=FUNKNO(L4+IOF,IG)*ZCODE(1)
         ELSE
            XNI(IEL,J)=FUNKNO(L4+IOF,IG)*ZCODE(2)
         ENDIF
      ENDIF
  105 CONTINUE
      IF(MAT(I,J).EQ.0) GO TO 140
      DO 115 IEL=1,IELEM**2
      Q(IEL)=0.0D0
      DO 110 K=1,NSCT
      IOF=((J-1)*LX*NSCT+(I-1)*NSCT+(K-1))*IELEM*IELEM+IEL
      Q(IEL)=Q(IEL)+QEXT(IOF,IG)*PL(K,M)/(4.0D0*PI)
  110 CONTINUE
  115 CONTINUE

      VT=VOL(I,J)*TOTAL(MAT(I,J),IG)
      CALL XDDSET(Q2,(IELEM**2)*((IELEM**2)+1),0.0D0)
      IF(IELEM.EQ.1) THEN
         Q2(1,1) = VT + ABS(DB(I,M)) + ABS(DA(I,J,M))
         Q2(1,2) = Q(1)*VOL(I,J)
         Q2(1,2) = Q2(1,2) + XNJ(1)*(SIGN(1.0,DE(M))*DB(I,M) +
     >      SIGN(1.0,DU(M))*DA(I,J,M))
      ELSEIF(IELEM.EQ.2) THEN
*     alpha=1 , beta=1
      Q2(1,1) = VT/4 + ABS(DB(I,M))/4 + ABS(DA(I,J,M))/4
      Q2(1,2) = (3**(1.0D0/2.0D0)*DA(I,J,M))/4
      Q2(1,3) = (3**(1.0D0/2.0D0)*DB(I,M))/4
      Q2(1,4) = 0.0D0

*     alpha=2 , beta=1
      Q2(2,1) = -(3**(1.0D0/2.0D0)*DA(I,J,M))/12
      Q2(2,2) = VT/12 + ABS(DB(I,M))/12 + ABS(DA(I,J,M))/4
      Q2(2,3) = 0.0D0
      Q2(2,4) = (3**(1.0D0/2.0D0)*DB(I,M))/12

*     alpha=1 , beta=2
      Q2(3,1) = -(3**(1.0D0/2.0D0)*DB(I,M))/12
      Q2(3,2) = 0.0D0
      Q2(3,3) = VT/12 + ABS(DB(I,M))/4 + ABS(DA(I,J,M))/12
      Q2(3,4) = (3**(1.0D0/2.0D0)*DA(I,J,M))/12

*     alpha=2 , beta=2
      Q2(4,1) = 0.0D0
      Q2(4,2) = -(3**(1.0D0/2.0D0)*DB(I,M))/36
      Q2(4,3) = -(3**(1.0D0/2.0D0)*DA(I,J,M))/36
      Q2(4,4) = VT/36 + ABS(DB(I,M))/12 + ABS(DA(I,J,M))/12

*     source
      Q2(1,5) = Q(1)*VOL(I,J)/4.0D0
      Q2(2,5) = Q(2)*VOL(I,J)/12.0D0
      Q2(3,5) = Q(3)*VOL(I,J)/12.0D0
      Q2(4,5) = Q(4)*VOL(I,J)/36.0D0

*     jump terms (dependent on direction)
      Q2(1,5)= Q2(1,5)
     >   + (1.0D0/8.0D0)*(
     >   DA(I,J,M)*SIGN(1.0,DU(M))*(XNI(1,J)+XNI(2,J)) +
     >   DB(I,M)*SIGN(1.0,DE(M))*(XNJ(1)+XNJ(2))   )
      Q2(2,5)= Q2(2,5)
     >   -SQRT(3.0D0)*(1.0D0/72.0D0)*(
     >   3.0D0*DA(I,J,M)*(XNI(1,J)+XNI(2,J)) +
     >   DB(I,M)*SIGN(1.0,DE(M))*(XNJ(1)-XNJ(2))   )
      Q2(3,5)= Q2(3,5)
     >   -SQRT(3.0D0)*(1.0D0/72.0D0)*(
     >   DA(I,J,M)*SIGN(1.0,DU(M))*(XNI(1,J)-XNI(2,J)) +
     >   3.0D0*DB(I,M)*(XNJ(1)+XNJ(2))   )
      Q2(4,5)= Q2(4,5)
     >   + (1.0D0/72.0D0)*(
     >   DA(I,J,M)*(XNI(1,J)-XNI(2,J)) +
     >   DB(I,M)*(XNJ(1)-XNJ(2))   )
*** ---------------------------------------------------------------- ***
      ELSEIF(IELEM.EQ.3) THEN
      Q2(1,1) = VT/4 + ABS(DB(I,M))/12 + ABS(DA(I,J,M))/12
      Q2(1,2) = (5*3**(0.5D0)*DA(I,J,M))/12
      Q2(1,3) = -(5**(0.5D0)*(6*VT + 2*ABS(DB(I,M)) -
     > 15*ABS(DA(I,J,M))))/180
      Q2(1,4) = (5*3**(0.5D0)*DB(I,M))/12
      Q2(1,5) = 0.0D0
      Q2(1,6) = -(15**(0.5D0)*DB(I,M))/18
      Q2(1,7) = -(5**(0.5D0)*(6*VT - 15*ABS(DB(I,M)) +
     > 2*ABS(DA(I,J,M))))/180
      Q2(1,8) = -(15**(0.5D0)*DA(I,J,M))/18
      Q2(1,9) = VT/45 - ABS(DB(I,M))/18 - ABS(DA(I,J,M))/18

      Q2(2,1) = -(3**(0.5D0)*DA(I,J,M))/12
      Q2(2,2) = VT/12 + ABS(DB(I,M))/36 + ABS(DA(I,J,M))/4
      Q2(2,3) = (15**(0.5D0)*DA(I,J,M))/12
      Q2(2,4) = 0.0D0
      Q2(2,5) = (5*3**(0.5D0)*DB(I,M))/36
      Q2(2,6) = 0.0D0
      Q2(2,7) = (15**(0.5D0)*DA(I,J,M))/90
      Q2(2,8) = -(5**(0.5D0)*(2*VT - 5*ABS(DB(I,M)) +
     > 6*ABS(DA(I,J,M))))/180
      Q2(2,9) = -(3**(0.5D0)*DA(I,J,M))/18

      Q2(3,1) = -(5**(0.5D0)*(3*VT + ABS(DB(I,M)) -
     > 3*ABS(DA(I,J,M))))/90
      Q2(3,2) = -(15**(0.5D0)*DA(I,J,M))/10
      Q2(3,3) = VT/15 + ABS(DB(I,M))/45 + ABS(DA(I,J,M))/6
      Q2(3,4) = -(15**(0.5D0)*DB(I,M))/18
      Q2(3,5) = 0.0D0
      Q2(3,6) = (3**(0.5D0)*DB(I,M))/9
      Q2(3,7) = VT/45 - ABS(DB(I,M))/18 - ABS(DA(I,J,M))/45
      Q2(3,8) = (3**(0.5D0)*DA(I,J,M))/15
      Q2(3,9) = -(5**(0.5D0)*(2*VT - 5*ABS(DB(I,M)) +
     > 5*ABS(DA(I,J,M))))/225

      Q2(4,1) = -(3**(0.5D0)*DB(I,M))/12
      Q2(4,2) = 0.0D0
      Q2(4,3) = (15**(0.5D0)*DB(I,M))/90
      Q2(4,4) = VT/12 + ABS(DB(I,M))/4 + ABS(DA(I,J,M))/36
      Q2(4,5) = (5*3**(0.5D0)*DA(I,J,M))/36
      Q2(4,6) = -(5**(0.5D0)*(2*VT + 6*ABS(DB(I,M)) -
     > 5*ABS(DA(I,J,M))))/180
      Q2(4,7) = (15**(0.5D0)*DB(I,M))/12
      Q2(4,8) = 0.0D0
      Q2(4,9) = -(3**(0.5D0)*DB(I,M))/18

      Q2(5,1) = 0.0D0
      Q2(5,2) = -(3**(0.5D0)*DB(I,M))/36
      Q2(5,3) = 0.0D0
      Q2(5,4) = -(3**(0.5D0)*DA(I,J,M))/36
      Q2(5,5) = VT/36 + ABS(DB(I,M))/12 + ABS(DA(I,J,M))/12
      Q2(5,6) = (15**(0.5D0)*DA(I,J,M))/36
      Q2(5,7) = 0.0D0
      Q2(5,8) = (15**(0.5D0)*DB(I,M))/36
      Q2(5,9) = 0.0D0

      Q2(6,1) = (15**(0.5D0)*DB(I,M))/90
      Q2(6,2) = 0.0D0
      Q2(6,3) = -(3**(0.5D0)*DB(I,M))/45
      Q2(6,4) = -(5**(0.5D0)*(VT + 3*ABS(DB(I,M)) -
     > ABS(DA(I,J,M))))/90
      Q2(6,5) = -(15**(0.5D0)*DA(I,J,M))/30
      Q2(6,6) = VT/45 + ABS(DB(I,M))/15 + ABS(DA(I,J,M))/18
      Q2(6,7) = -(3**(0.5D0)*DB(I,M))/18
      Q2(6,8) = 0.0D0
      Q2(6,9) = (15**(0.5D0)*DB(I,M))/45

      Q2(7,1) = -(5**(0.5D0)*(3*VT - 3*ABS(DB(I,M)) +
     > ABS(DA(I,J,M))))/90
      Q2(7,2) = -(15**(0.5D0)*DA(I,J,M))/18
      Q2(7,3) = VT/45 - ABS(DB(I,M))/45 - ABS(DA(I,J,M))/18
      Q2(7,4) = -(15**(0.5D0)*DB(I,M))/10
      Q2(7,5) = 0.0D0
      Q2(7,6) = (3**(0.5D0)*DB(I,M))/15
      Q2(7,7) = VT/15 + ABS(DB(I,M))/6 + ABS(DA(I,J,M))/45
      Q2(7,8) = (3**(0.5D0)*DA(I,J,M))/9
      Q2(7,9) = -(5**(0.5D0)*(2*VT + 5*ABS(DB(I,M)) -
     > 5*ABS(DA(I,J,M))))/225

      Q2(8,1) = (15**(0.5D0)*DA(I,J,M))/90
      Q2(8,2) = -(5**(0.5D0)*(VT - ABS(DB(I,M)) +
     > 3*ABS(DA(I,J,M))))/90
      Q2(8,3) = -(3**(0.5D0)*DA(I,J,M))/18
      Q2(8,4) = 0.0D0
      Q2(8,5) = -(15**(0.5D0)*DB(I,M))/30
      Q2(8,6) = 0.0D0
      Q2(8,7) = -(3**(0.5D0)*DA(I,J,M))/45
      Q2(8,8) = VT/45 + ABS(DB(I,M))/18 + ABS(DA(I,J,M))/15
      Q2(8,9) = (15**(0.5D0)*DA(I,J,M))/45

      Q2(9,1) = VT/45 - ABS(DB(I,M))/45 - ABS(DA(I,J,M))/45
      Q2(9,2) = (3**(0.5D0)*DA(I,J,M))/15
      Q2(9,3) = -(5**(0.5D0)*(2*VT - 2*ABS(DB(I,M)) +
     > 5*ABS(DA(I,J,M))))/225
      Q2(9,4) = (3**(0.5D0)*DB(I,M))/15
      Q2(9,5) = 0.0D0
      Q2(9,6) = -(2*15**(0.5D0)*DB(I,M))/75
      Q2(9,7) = -(5**(0.5D0)*(2*VT + 5*ABS(DB(I,M)) -
     > 2*ABS(DA(I,J,M))))/225
      Q2(9,8) = -(2*15**(0.5D0)*DA(I,J,M))/75
      Q2(9,9) = (4*VT)/225 + (2*ABS(DB(I,M)))/45 +
     > (2*ABS(DA(I,J,M)))/45

      Q2(1,10) = (((45*Q(01) + 4*Q(09) - 6*5**(0.5D0)*Q(03) -
     > 6*5**(0.5D0)*Q(07)))/180)*VOL(I,J)
      Q2(2,10) = (((15*Q(02) - 2*5**(0.5D0)*Q(08)))/180)*VOL(I,J)
      Q2(3,10) = (((30*Q(03) + 10*Q(07) - 15*5**(0.5D0)*Q(01) -
     > 4*5**(0.5D0)*Q(09)))/450)*VOL(I,J)
      Q2(4,10) = (((15*Q(04) - 2*5**(0.5D0)*Q(06)))/180)*VOL(I,J)
      Q2(5,10) = ((Q(05))/36)*VOL(I,J)
      Q2(6,10) = (((2*Q(06) - 5**(0.5D0)*Q(04)))/90)*VOL(I,J)
      Q2(7,10) = (((10*Q(03) + 30*Q(07) - 15*5**(0.5D0)*Q(01) -
     > 4*5**(0.5D0)*Q(09)))/450)*VOL(I,J)
      Q2(8,10) = (((2*Q(08) - 5**(0.5D0)*Q(02)))/90)*VOL(I,J)
      Q2(9,10) = (((900*Q(01) + 720*Q(09) - 360*5**(0.5D0)*Q(03) -
     > 360*5**(0.5D0)*Q(07)))/40500)*VOL(I,J)

      Q2(1,10) = Q2(1,10) +
     >   ((11*XNI(3,J))/1080 + (17*XNI(2,J))/270 + (11*XNI(1,J))/
     >   1080)*DA(I,J,M)*SIGN(1.0,DU(M)) + ((11*XNJ(1))/1080 +
     >   (17*XNJ(2))/270 + (11*XNJ(3))/1080)*DB(I,M)*SIGN(1.0,DE(M))
      Q2(2,10) = Q2(2,10) +
     >   (-(3**(0.5D0)*(11*XNI(3,J) + 68*XNI(2,J) + 11*XNI(1,J)))/
     >   1080)*DA(I,J,M) + (-(3**(0.5D0)*(5*XNJ(1) - 5*XNJ(3)))/
     >   1080)*DB(I,M)*SIGN(1.0,DE(M))
      Q2(3,10) = Q2(3,10) +
     >   ((5**(0.5D0)*(11*XNI(3,J) + 68*XNI(2,J) + 11*XNI(1,J)))/
     >   2700)*DA(I,J,M)*SIGN(1.0,DU(M)) + (-(5**(0.5D0)*(XNJ(1) +
     >   28*XNJ(2) + XNJ(3)))/2700)*DB(I,M)*SIGN(1.0,DE(M))
      Q2(4,10) = Q2(4,10) +
     >   ((3**(0.5D0)*(5*XNI(3,J) - 5*XNI(1,J)))/1080)*DA(I,J,M)*
     >   SIGN(1.0,DU(M)) + (-(3**(0.5D0)*(11*XNJ(1) + 68*XNJ(2) +
     >   11*XNJ(3)))/1080)*DB(I,M)
      Q2(5,10) = Q2(5,10) +
     >   (XNI(1,J)/72 - XNI(3,J)/72)*DA(I,J,M) + (XNJ(1)/72 -
     >   XNJ(3)/72)*DB(I,M)
      Q2(6,10) = Q2(6,10) +
     >   ((3**(0.5D0)*5**(0.5D0)*(5*XNI(3,J) - 5*XNI(1,J)))/2700)*
     >   DA(I,J,M)*SIGN(1.0,DU(M)) + ((3**(0.5D0)*5**(0.5D0)*
     >   (XNJ(1) + 28*XNJ(2) + XNJ(3)))/2700)*DB(I,M)
      Q2(7,10) = Q2(7,10) +
     >   (-(5**(0.5D0)*(XNI(3,J) + 28*XNI(2,J) + XNI(1,J)))/
     >   2700)*DA(I,J,M)*SIGN(1.0,DU(M)) + ((5**(0.5D0)*(11*XNJ(1) +
     >   68*XNJ(2) + 11*XNJ(3)))/2700)*DB(I,M)*SIGN(1.0,DE(M))
      Q2(8,10) = Q2(8,10) +
     >   ((3**(0.5D0)*5**(0.5D0)*(XNI(3,J) + 28*XNI(2,J) +
     >   XNI(1,J)))/2700)*DA(I,J,M) + (-(3**(0.5D0)*5**(0.5D0)*
     >   (5*XNJ(1) - 5*XNJ(3)))/2700)*DB(I,M)*SIGN(1.0,DE(M))
      Q2(9,10) = Q2(9,10) +
     >   (- XNI(3,J)/1350 - (14*XNI(2,J))/675 - XNI(1,J)/1350)*
     >   DA(I,J,M)*SIGN(1.0,DU(M)) + (- XNJ(1)/1350 - (14*XNJ(2))/
     >   675 - XNJ(3)/1350)*DB(I,M)*SIGN(1.0,DE(M))
*** ---------------------------------------------------------------- ***
      ELSEIF(IELEM.EQ.4) THEN
      Q2(1,1) = (25*VT)/256 + (5*ABS(DB(I,M)))/128 +
     >   (5*ABS(DA(I,J,M)))/128
      Q2(1,2) = (5*3**(0.5D0)*DA(I,J,M))/32
      Q2(1,3) = -(5**(0.5D0)*(15*VT + 6*ABS(DB(I,M)) -
     >   50*ABS(DA(I,J,M))))/1280
      Q2(1,4) = (5*7**(0.5D0)*DA(I,J,M))/128
      Q2(1,5) = (5*3**(0.5D0)*DB(I,M))/32
      Q2(1,6) = 0.0D0
      Q2(1,7) = -(3*15**(0.5D0)*DB(I,M))/160
      Q2(1,8) = 0.0D0
      Q2(1,9) = -(5**(0.5D0)*(15*VT - 50*ABS(DB(I,M)) +
     >   6*ABS(DA(I,J,M))))/1280
      Q2(1,10) = -(3*15**(0.5D0)*DA(I,J,M))/160
      Q2(1,11) = (9*VT)/1280 - (3*ABS(DB(I,M)))/128 -
     >   (3*ABS(DA(I,J,M)))/128
      Q2(1,12) = -(3*35**(0.5D0)*DA(I,J,M))/640
      Q2(1,13) = (5*7**(0.5D0)*DB(I,M))/128
      Q2(1,14) = 0.0D0
      Q2(1,15) = -(3*35**(0.5D0)*DB(I,M))/640
      Q2(1,16) = 0.0D0

      Q2(2,1) = -(11*3**(0.5D0)*DA(I,J,M))/384
      Q2(2,2) = (85*VT)/768 + (17*ABS(DB(I,M)))/384 +
     >   (11*ABS(DA(I,J,M)))/128
      Q2(2,3) = (37*15**(0.5D0)*DA(I,J,M))/192
      Q2(2,4) = -(3**(0.5D0)*7**(0.5D0)*(45*VT +
     >   18*ABS(DB(I,M)) - 110*ABS(DA(I,J,M))))/3840
      Q2(2,5) = 0.0D0
      Q2(2,6) = (17*3**(0.5D0)*DB(I,M))/96
      Q2(2,7) = 0.0D0
      Q2(2,8) = -(9*7**(0.5D0)*DB(I,M))/160
      Q2(2,9) = (11*15**(0.5D0)*DA(I,J,M))/3200
      Q2(2,10) = -(5**(0.5D0)*(255*VT - 850*ABS(DB(I,M)) +
     >   198*ABS(DA(I,J,M))))/19200
      Q2(2,11) = -(37*3**(0.5D0)*DA(I,J,M))/320
      Q2(2,12) = -(3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(30*ABS(DB(I,M)) -
     >   9*VT + 22*ABS(DA(I,J,M))))/6400
      Q2(2,13) = 0.0D0
      Q2(2,14) = (17*7**(0.5D0)*DB(I,M))/384
      Q2(2,15) = 0.0D0
      Q2(2,16) = -(21*3**(0.5D0)*DB(I,M))/640

      Q2(3,1) = -(3*5**(0.5D0)*(5*VT + 2*ABS(DB(I,M)) -
     >   10*ABS(DA(I,J,M))))/1280
      Q2(3,2) = -(3*15**(0.5D0)*DA(I,J,M))/64
      Q2(3,3) = (9*VT)/256 + (9*ABS(DB(I,M)))/640 +
     >   (15*ABS(DA(I,J,M)))/128
      Q2(3,4) = (3*35**(0.5D0)*DA(I,J,M))/128
      Q2(3,5) = -(3*15**(0.5D0)*DB(I,M))/160
      Q2(3,6) = 0.0D0
      Q2(3,7) = (9*3**(0.5D0)*DB(I,M))/160
      Q2(3,8) = 0.0D0
      Q2(3,9) = (9*VT)/1280 - (3*ABS(DB(I,M)))/128 -
     >   (9*ABS(DA(I,J,M)))/640
      Q2(3,10) = (9*3**(0.5D0)*DA(I,J,M))/320
      Q2(3,11) = -(9*5**(0.5D0)*(3*VT - 10*ABS(DB(I,M)) +
     >   10*ABS(DA(I,J,M))))/6400
      Q2(3,12) = -(9*7**(0.5D0)*DA(I,J,M))/640
      Q2(3,13) = -(3*35**(0.5D0)*DB(I,M))/640
      Q2(3,14) = 0.0D0
      Q2(3,15) = (9*7**(0.5D0)*DB(I,M))/640
      Q2(3,16) = 0.0D0

      Q2(4,1) = -(9*7**(0.5D0)*DA(I,J,M))/896
      Q2(4,2) = -(3*3**(0.5D0)*7**(0.5D0)*(35*VT +
     >   14*ABS(DB(I,M)) - 30*ABS(DA(I,J,M))))/8960
      Q2(4,3) = -(9*35**(0.5D0)*DA(I,J,M))/112
      Q2(4,4) = (81*VT)/1792 + (81*ABS(DB(I,M)))/4480 +
     >   (9*ABS(DA(I,J,M)))/128
      Q2(4,5) = 0.0D0
      Q2(4,6) = -(9*7**(0.5D0)*DB(I,M))/160
      Q2(4,7) = 0.0D0
      Q2(4,8) = (81*3**(0.5D0)*DB(I,M))/1120
      Q2(4,9) = (27*35**(0.5D0)*DA(I,J,M))/22400
      Q2(4,10) = -(3*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*
     >   (70*ABS(DB(I,M)) - 21*VT + 18*ABS(DA(I,J,M))))/44800
      Q2(4,11) = (27*7**(0.5D0)*DA(I,J,M))/560
      Q2(4,12) = -(27*5**(0.5D0)*(9*VT - 30*ABS(DB(I,M)) +
     >   14*ABS(DA(I,J,M))))/44800
      Q2(4,13) = 0.0D0
      Q2(4,14) = -(21*3**(0.5D0)*DB(I,M))/640
      Q2(4,15) = 0.0D0
      Q2(4,16) = (81*7**(0.5D0)*DB(I,M))/4480

      Q2(5,1) = -(11*3**(0.5D0)*DB(I,M))/384
      Q2(5,2) = 0.0D0
      Q2(5,3) = (11*15**(0.5D0)*DB(I,M))/3200
      Q2(5,4) = 0.0D0
      Q2(5,5) = (85*VT)/768 + (11*ABS(DB(I,M)))/128 +
     >   (17*ABS(DA(I,J,M)))/384
      Q2(5,6) = (17*3**(0.5D0)*DA(I,J,M))/96
      Q2(5,7) = -(5**(0.5D0)*(255*VT + 198*ABS(DB(I,M)) -
     >   850*ABS(DA(I,J,M))))/19200
      Q2(5,8) = (17*7**(0.5D0)*DA(I,J,M))/384
      Q2(5,9) = (37*15**(0.5D0)*DB(I,M))/192
      Q2(5,10) = 0.0D0
      Q2(5,11) = -(37*3**(0.5D0)*DB(I,M))/320
      Q2(5,12) = 0.0D0
      Q2(5,13) = -(3**(0.5D0)*7**(0.5D0)*(45*VT -
     >   110*ABS(DB(I,M)) + 18*ABS(DA(I,J,M))))/3840
      Q2(5,14) = -(9*7**(0.5D0)*DA(I,J,M))/160
      Q2(5,15) = -(3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(22*ABS(DB(I,M)) -
     >   9*VT + 30*ABS(DA(I,J,M))))/6400
      Q2(5,16) = -(21*3**(0.5D0)*DA(I,J,M))/640

      Q2(6,1) = 0.0D0
      Q2(6,2) = -(187*3**(0.5D0)*DB(I,M))/5760
      Q2(6,3) = 0.0D0
      Q2(6,4) = (33*7**(0.5D0)*DB(I,M))/3200
      Q2(6,5) = -(187*3**(0.5D0)*DA(I,J,M))/5760
      Q2(6,6) = (289*VT)/2304 + (187*ABS(DB(I,M)))/1920 +
     >   (187*ABS(DA(I,J,M)))/1920
      Q2(6,7) = (629*15**(0.5D0)*DA(I,J,M))/2880
      Q2(6,8) = -(3**(0.5D0)*7**(0.5D0)*(765*VT +
     >   594*ABS(DB(I,M)) - 1870*ABS(DA(I,J,M))))/57600
      Q2(6,9) = 0.0D0
      Q2(6,10) = (629*15**(0.5D0)*DB(I,M))/2880
      Q2(6,11) = 0.0D0
      Q2(6,12) = -(111*35**(0.5D0)*DB(I,M))/1600
      Q2(6,13) = (33*7**(0.5D0)*DA(I,J,M))/3200
      Q2(6,14) = -(3**(0.5D0)*7**(0.5D0)*(765*VT -
     >   1870*ABS(DB(I,M)) + 594*ABS(DA(I,J,M))))/57600
      Q2(6,15) = -(111*35**(0.5D0)*DA(I,J,M))/1600
      Q2(6,16) = (189*VT)/6400 - (231*ABS(DB(I,M)))/3200 -
     >   (231*ABS(DA(I,J,M)))/3200

      Q2(7,1) = (11*15**(0.5D0)*DB(I,M))/3200
      Q2(7,2) = 0.0D0
      Q2(7,3) = -(33*3**(0.5D0)*DB(I,M))/3200
      Q2(7,4) = 0.0D0
      Q2(7,5) = -(5**(0.5D0)*(85*VT + 66*ABS(DB(I,M)) -
     >   170*ABS(DA(I,J,M))))/6400
      Q2(7,6) = -(17*15**(0.5D0)*DA(I,J,M))/320
      Q2(7,7) = (51*VT)/1280 + (99*ABS(DB(I,M)))/3200 +
     >   (17*ABS(DA(I,J,M)))/128
      Q2(7,8) = (17*35**(0.5D0)*DA(I,J,M))/640
      Q2(7,9) = -(37*3**(0.5D0)*DB(I,M))/320
      Q2(7,10) = 0.0D0
      Q2(7,11) = (111*15**(0.5D0)*DB(I,M))/1600
      Q2(7,12) = 0.0D0
      Q2(7,13) = -(3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*
     >   (22*ABS(DB(I,M)) - 9*VT + 18*ABS(DA(I,J,M))))/6400
      Q2(7,14) = (27*35**(0.5D0)*DA(I,J,M))/1600
      Q2(7,15) = -(3*3**(0.5D0)*7**(0.5D0)*(9*VT -
     >   22*ABS(DB(I,M)) + 30*ABS(DA(I,J,M))))/6400
      Q2(7,16) = -(63*15**(0.5D0)*DA(I,J,M))/3200

      Q2(8,1) = 0.0D0
      Q2(8,2) = (33*7**(0.5D0)*DB(I,M))/3200
      Q2(8,3) = 0.0D0
      Q2(8,4) = -(297*3**(0.5D0)*DB(I,M))/22400
      Q2(8,5) = -(51*7**(0.5D0)*DA(I,J,M))/4480
      Q2(8,6) = -(3**(0.5D0)*7**(0.5D0)*(595*VT +
     >   462*ABS(DB(I,M)) - 510*ABS(DA(I,J,M))))/44800
      Q2(8,7) = -(51*35**(0.5D0)*DA(I,J,M))/560
      Q2(8,8) = (459*VT)/8960 + (891*ABS(DB(I,M)))/22400 +
     >   (51*ABS(DA(I,J,M)))/640
      Q2(8,9) = 0.0D0
      Q2(8,10) = -(111*35**(0.5D0)*DB(I,M))/1600
      Q2(8,11) = 0.0D0
      Q2(8,12) = (999*15**(0.5D0)*DB(I,M))/11200
      Q2(8,13) = (27*3**(0.5D0)*DA(I,J,M))/3200
      Q2(8,14) = (189*VT)/6400 - (231*ABS(DB(I,M)))/3200 -
     >   (81*ABS(DA(I,J,M)))/3200
      Q2(8,15) = (27*15**(0.5D0)*DA(I,J,M))/400
      Q2(8,16) = -(27*3**(0.5D0)*7**(0.5D0)*(9*VT -
     >   22*ABS(DB(I,M)) + 14*ABS(DA(I,J,M))))/44800

      Q2(9,1) = -(3*5**(0.5D0)*(5*VT - 10*ABS(DB(I,M)) +
     >   2*ABS(DA(I,J,M))))/1280
      Q2(9,2) = -(3*15**(0.5D0)*DA(I,J,M))/160
      Q2(9,3) = (9*VT)/1280 - (9*ABS(DB(I,M)))/640 -
     >   (3*ABS(DA(I,J,M)))/128
      Q2(9,4) = -(3*35**(0.5D0)*DA(I,J,M))/640
      Q2(9,5) = -(3*15**(0.5D0)*DB(I,M))/64
      Q2(9,6) = 0.0D0
      Q2(9,7) = (9*3**(0.5D0)*DB(I,M))/320
      Q2(9,8) = 0.0D0
      Q2(9,9) = (9*VT)/256 + (15*ABS(DB(I,M)))/128 +
     >   (9*ABS(DA(I,J,M)))/640
      Q2(9,10) = (9*3**(0.5D0)*DA(I,J,M))/160
      Q2(9,11) = -(9*5**(0.5D0)*(3*VT + 10*ABS(DB(I,M)) -
     >   10*ABS(DA(I,J,M))))/6400
      Q2(9,12) = (9*7**(0.5D0)*DA(I,J,M))/640
      Q2(9,13) = (3*35**(0.5D0)*DB(I,M))/128
      Q2(9,14) = 0.0D0
      Q2(9,15) = -(9*7**(0.5D0)*DB(I,M))/640
      Q2(9,16) = 0.0D0

      Q2(10,1) = (11*15**(0.5D0)*DA(I,J,M))/3200
      Q2(10,2) = -(5**(0.5D0)*(85*VT - 170*ABS(DB(I,M)) +
     >   66*ABS(DA(I,J,M))))/6400
      Q2(10,3) = -(37*3**(0.5D0)*DA(I,J,M))/320
      Q2(10,4) = -(3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(18*ABS(DB(I,M)) -
     >   9*VT + 22*ABS(DA(I,J,M))))/6400
      Q2(10,5) = 0.0D0
      Q2(10,6) = -(17*15**(0.5D0)*DB(I,M))/320
      Q2(10,7) = 0.0D0
      Q2(10,8) = (27*35**(0.5D0)*DB(I,M))/1600
      Q2(10,9) = -(33*3**(0.5D0)*DA(I,J,M))/3200
      Q2(10,10) = (51*VT)/1280 + (17*ABS(DB(I,M)))/128 +
     >   (99*ABS(DA(I,J,M)))/3200
      Q2(10,11) = (111*15**(0.5D0)*DA(I,J,M))/1600
      Q2(10,12) = -(3*3**(0.5D0)*7**(0.5D0)*(9*VT + 30*ABS(DB(I,M)) -
     >   22*ABS(DA(I,J,M))))/6400
      Q2(10,13) = 0.0D0
      Q2(10,14) = (17*35**(0.5D0)*DB(I,M))/640
      Q2(10,15) = 0.0D0
      Q2(10,16) = -(63*15**(0.5D0)*DB(I,M))/3200

      Q2(11,1) = (9*VT)/1280 - (9*ABS(DB(I,M)))/640 -
     >   (9*ABS(DA(I,J,M)))/640
      Q2(11,2) = (9*3**(0.5D0)*DA(I,J,M))/320
      Q2(11,3) = -(9*5**(0.5D0)*(3*VT - 6*ABS(DB(I,M)) +
     >   10*ABS(DA(I,J,M))))/6400
      Q2(11,4) = -(9*7**(0.5D0)*DA(I,J,M))/640
      Q2(11,5) = (9*3**(0.5D0)*DB(I,M))/320
      Q2(11,6) = 0.0D0
      Q2(11,7) = -(27*15**(0.5D0)*DB(I,M))/1600
      Q2(11,8) = 0.0D0
      Q2(11,9) = -(9*5**(0.5D0)*(3*VT + 10*ABS(DB(I,M)) -
     >   6*ABS(DA(I,J,M))))/6400
      Q2(11,10) = -(27*15**(0.5D0)*DA(I,J,M))/1600
      Q2(11,11) = (81*VT)/6400 + (27*ABS(DB(I,M)))/640 +
     >   (27*ABS(DA(I,J,M)))/640
      Q2(11,12) = (27*35**(0.5D0)*DA(I,J,M))/3200
      Q2(11,13) = -(9*7**(0.5D0)*DB(I,M))/640
      Q2(11,14) = 0.0D0
      Q2(11,15) = (27*35**(0.5D0)*DB(I,M))/3200
      Q2(11,16) = 0.0D0

      Q2(12,1) = (27*35**(0.5D0)*DA(I,J,M))/22400
      Q2(12,2) = -(9*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*
     >   (14*ABS(DB(I,M)) - 7*VT + 6*ABS(DA(I,J,M))))/44800
      Q2(12,3) = (27*7**(0.5D0)*DA(I,J,M))/560
      Q2(12,4) = -(27*5**(0.5D0)*(9*VT - 18*ABS(DB(I,M)) +
     >   14*ABS(DA(I,J,M))))/44800
      Q2(12,5) = 0.0D0
      Q2(12,6) = (27*35**(0.5D0)*DB(I,M))/1600
      Q2(12,7) = 0.0D0
      Q2(12,8) = -(243*15**(0.5D0)*DB(I,M))/11200
      Q2(12,9) = -(81*7**(0.5D0)*DA(I,J,M))/22400
      Q2(12,10) = -(9*3**(0.5D0)*7**(0.5D0)*(21*VT +
     >   70*ABS(DB(I,M)) - 18*ABS(DA(I,J,M))))/44800
      Q2(12,11) = -(81*35**(0.5D0)*DA(I,J,M))/2800
      Q2(12,12) = (729*VT)/44800 + (243*ABS(DB(I,M)))/4480 +
     >   (81*ABS(DA(I,J,M)))/3200
      Q2(12,13) = 0.0D0
      Q2(12,14) = -(63*15**(0.5D0)*DB(I,M))/3200
      Q2(12,15) = 0.0D0
      Q2(12,16) = (243*35**(0.5D0)*DB(I,M))/22400

      Q2(13,1) = -(9*7**(0.5D0)*DB(I,M))/896
      Q2(13,2) = 0.0D0
      Q2(13,3) = (27*35**(0.5D0)*DB(I,M))/22400
      Q2(13,4) = 0.0D0
      Q2(13,5) = -(3*3**(0.5D0)*7**(0.5D0)*(35*VT -
     >   30*ABS(DB(I,M)) + 14*ABS(DA(I,J,M))))/8960
      Q2(13,6) = -(9*7**(0.5D0)*DA(I,J,M))/160
      Q2(13,7) = -(3*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*
     >   (18*ABS(DB(I,M)) - 21*VT + 70*ABS(DA(I,J,M))))/44800
      Q2(13,8) = -(21*3**(0.5D0)*DA(I,J,M))/640
      Q2(13,9) = -(9*35**(0.5D0)*DB(I,M))/112
      Q2(13,10) = 0.0D0
      Q2(13,11) = (27*7**(0.5D0)*DB(I,M))/560
      Q2(13,12) = 0.0D0
      Q2(13,13) = (81*VT)/1792 + (9*ABS(DB(I,M)))/128 +
     >   (81*ABS(DA(I,J,M)))/4480
      Q2(13,14) = (81*3**(0.5D0)*DA(I,J,M))/1120
      Q2(13,15) = -(27*5**(0.5D0)*(9*VT + 14*ABS(DB(I,M)) -
     >   30*ABS(DA(I,J,M))))/44800
      Q2(13,16) = (81*7**(0.5D0)*DA(I,J,M))/4480

      Q2(14,1) = 0.0D0
      Q2(14,2) = -(51*7**(0.5D0)*DB(I,M))/4480
      Q2(14,3) = 0.0D0
      Q2(14,4) = (27*3**(0.5D0)*DB(I,M))/3200
      Q2(14,5) = (33*7**(0.5D0)*DA(I,J,M))/3200
      Q2(14,6) = -(3**(0.5D0)*7**(0.5D0)*(595*VT -
     >   510*ABS(DB(I,M)) + 462*ABS(DA(I,J,M))))/44800
      Q2(14,7) = -(111*35**(0.5D0)*DA(I,J,M))/1600
      Q2(14,8) = (189*VT)/6400 - (81*ABS(DB(I,M)))/3200 -
     >   (231*ABS(DA(I,J,M)))/3200
      Q2(14,9) = 0.0D0
      Q2(14,10) = -(51*35**(0.5D0)*DB(I,M))/560
      Q2(14,11) = 0.0D0
      Q2(14,12) = (27*15**(0.5D0)*DB(I,M))/400
      Q2(14,13) = -(297*3**(0.5D0)*DA(I,J,M))/22400
      Q2(14,14) = (459*VT)/8960 + (51*ABS(DB(I,M)))/640 +
     >   (891*ABS(DA(I,J,M)))/22400
      Q2(14,15) = (999*15**(0.5D0)*DA(I,J,M))/11200
      Q2(14,16) = -(27*3**(0.5D0)*7**(0.5D0)*(9*VT +
     >   14*ABS(DB(I,M)) - 22*ABS(DA(I,J,M))))/44800

      Q2(15,1) = (27*35**(0.5D0)*DB(I,M))/22400
      Q2(15,2) = 0.0D0
      Q2(15,3) = -(81*7**(0.5D0)*DB(I,M))/22400
      Q2(15,4) = 0.0D0
      Q2(15,5) = -(9*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*
     >   (6*ABS(DB(I,M)) - 7*VT + 14*ABS(DA(I,J,M))))/44800
      Q2(15,6) = (27*35**(0.5D0)*DA(I,J,M))/1600
      Q2(15,7) = -(9*3**(0.5D0)*7**(0.5D0)*(21*VT -
     >   18*ABS(DB(I,M)) + 70*ABS(DA(I,J,M))))/44800
      Q2(15,8) = -(63*15**(0.5D0)*DA(I,J,M))/3200
      Q2(15,9) = (27*7**(0.5D0)*DB(I,M))/560
      Q2(15,10) = 0.0D0
      Q2(15,11) = -(81*35**(0.5D0)*DB(I,M))/2800
      Q2(15,12) = 0.0D0
      Q2(15,13) = -(27*5**(0.5D0)*(9*VT + 14*ABS(DB(I,M)) -
     >   18*ABS(DA(I,J,M))))/44800
      Q2(15,14) = -(243*15**(0.5D0)*DA(I,J,M))/11200
      Q2(15,15) = (729*VT)/44800 + (81*ABS(DB(I,M)))/3200 +
     >   (243*ABS(DA(I,J,M)))/4480
      Q2(15,16) = (243*35**(0.5D0)*DA(I,J,M))/22400

      Q2(16,1) = 0.0D0
      Q2(16,2) = (27*3**(0.5D0)*DB(I,M))/3200
      Q2(16,3) = 0.0D0
      Q2(16,4) = -(729*7**(0.5D0)*DB(I,M))/156800
      Q2(16,5) = (27*3**(0.5D0)*DA(I,J,M))/3200
      Q2(16,6) = (189*VT)/6400 - (81*ABS(DB(I,M)))/3200 -
     >   (81*ABS(DA(I,J,M)))/3200
      Q2(16,7) = (27*15**(0.5D0)*DA(I,J,M))/400
      Q2(16,8) = -(27*3**(0.5D0)*7**(0.5D0)*(63*VT -
     >   54*ABS(DB(I,M)) + 98*ABS(DA(I,J,M))))/313600
      Q2(16,9) = 0.0D0
      Q2(16,10) = (27*15**(0.5D0)*DB(I,M))/400
      Q2(16,11) = 0.0D0
      Q2(16,12) = -(729*35**(0.5D0)*DB(I,M))/19600
      Q2(16,13) = -(729*7**(0.5D0)*DA(I,J,M))/156800
      Q2(16,14) = -(27*3**(0.5D0)*7**(0.5D0)*(63*VT +
     >   98*ABS(DB(I,M)) - 54*ABS(DA(I,J,M))))/313600
      Q2(16,15) = -(729*35**(0.5D0)*DA(I,J,M))/19600
      Q2(16,16) = (6561*VT)/313600 + (729*ABS(DB(I,M)))/22400 +
     >   (729*ABS(DA(I,J,M)))/22400

      Q2(1,17) = (((2500*Q(01) + 180*Q(11) - 300*5**(0.5D0)*Q(03) -
     >   300*5**(0.5D0)*Q(09)))/25600)*VOL(I,J)
      Q2(2,17) = (((2125*Q(02) - 255*5**(0.5D0)*Q(10) -
     >   225*21**(0.5D0)*Q(04) + 27*105**(0.5D0)*Q(12)))/19200)*
     >   VOL(I,J)
      Q2(3,17) = (-(3*5**(0.5D0)*(500*Q(01) + 180*Q(11) -
     >   300*5**(0.5D0)*Q(03) - 60*5**(0.5D0)*Q(09)))/128000)*
     >   VOL(I,J)
      Q2(4,17) = ((3*(675*Q(04) - 81*5**(0.5D0)*Q(12) -
     >   175*21**(0.5D0)*Q(02) + 21*105**(0.5D0)*Q(10)))/44800)*
     >   VOL(I,J)
      Q2(5,17) = (((2125*Q(05) - 255*5**(0.5D0)*Q(07) -
     >   225*21**(0.5D0)*Q(13) + 27*105**(0.5D0)*Q(15)))/19200)*
     >   VOL(I,J)
      Q2(6,17) = (((7225*Q(06) + 1701*Q(16) - 765*21**(0.5D0)*Q(08) -
     >   765*21**(0.5D0)*Q(14)))/57600)*VOL(I,J)
      Q2(7,17) = (-(15**(0.5D0)*(8500*3**(0.5D0)*Q(05) -
     >   2700*7**(0.5D0)*Q(13) - 5100*15**(0.5D0)*Q(07) +
     >   1620*35**(0.5D0)*Q(15)))/1920000)*VOL(I,J)
      Q2(8,17) = (((2295*Q(08) + 1323*Q(14) - 595*21**(0.5D0)*Q(06) -
     >   243*21**(0.5D0)*Q(16)))/44800)*VOL(I,J)
      Q2(9,17) = (-(3*5**(0.5D0)*(500*Q(01) + 180*Q(11) -
     >   60*5**(0.5D0)*Q(03) - 300*5**(0.5D0)*Q(09)))/128000)*
     >   VOL(I,J)
      Q2(10,17) = (-(15**(0.5D0)*(8500*3**(0.5D0)*Q(02) -
     >   2700*7**(0.5D0)*Q(04) - 5100*15**(0.5D0)*Q(10) +
     >   1620*35**(0.5D0)*Q(12)))/1920000)*VOL(I,J)
      Q2(11,17) = ((9*(100*Q(01) + 180*Q(11) - 60*5**(0.5D0)*Q(03) -
     >   60*5**(0.5D0)*Q(09)))/128000)*VOL(I,J)
      Q2(12,17) = ((9*35**(0.5D0)*(4900*3**(0.5D0)*Q(02) -
     >   2700*7**(0.5D0)*Q(04) - 2940*15**(0.5D0)*Q(10) +
     >   1620*35**(0.5D0)*Q(12)))/31360000)*VOL(I,J)
      Q2(13,17) = ((3*(675*Q(13) - 81*5**(0.5D0)*Q(15) -
     >   175*21**(0.5D0)*Q(05) + 21*105**(0.5D0)*Q(07)))/44800)*
     >   VOL(I,J)
      Q2(14,17) = (((1323*Q(08) + 2295*Q(14) - 595*21**(0.5D0)*Q(06) -
     >   243*21**(0.5D0)*Q(16)))/44800)*VOL(I,J)
      Q2(15,17) = ((9*35**(0.5D0)*(4900*3**(0.5D0)*Q(05) -
     >   2700*7**(0.5D0)*Q(13) - 2940*15**(0.5D0)*Q(07) +
     >   1620*35**(0.5D0)*Q(15)))/31360000)*VOL(I,J)
      Q2(16,17) = ((9*(720300*Q(06) + 510300*Q(16) -
     >   132300*21**(0.5D0)*Q(08) - 132300*21**(0.5D0)*Q(14)))/
     >   219520000)*VOL(I,J)

      Q2(1,17) = Q2(1,17) +
     >   (XNI(4,J)/320 + (21*XNI(3,J))/1280 + (21*XNI(2,J))/1280 +
     >   XNI(1,J)/320)*DA(I,J,M)*SIGN(1.0,DU(M)) + (XNJ(1)/320 +
     >   (21*XNJ(2))/1280 + (21*XNJ(3))/1280 + XNJ(4)/320)*DB(I,M)*
     >   SIGN(1.0,DE(M))
      Q2(2,17) = Q2(2,17)
     >   - (3**(0.5D0)*DB(I,M)*SIGN(1.0,DE(M))*(173*XNJ(1) +
     >   756*XNJ(2) - 756*XNJ(3) - 173*XNJ(4)))/57600 - (3**(0.5D0)*
     >   DA(I,J,M)*(132*XNI(4,J) + 693*XNI(3,J) + 693*XNI(2,J) +
     >   132*XNI(1,J)))/57600
      Q2(3,17) = Q2(3,17) +
     >   (3*5**(0.5D0)*DA(I,J,M)*SIGN(1.0,DU(M))*(4*XNI(4,J) +
     >   21*XNI(3,J) + 21*XNI(2,J) + 4*XNI(1,J)))/6400 +
     >   (3*5**(0.5D0)*DB(I,M)*SIGN(1.0,DE(M))*(XNJ(1) - 6*XNJ(2) -
     >   6*XNJ(3) + XNJ(4)))/6400
      Q2(4,17) = Q2(4,17) +
     >   (3*7**(0.5D0)*DB(I,M)*SIGN(1.0,DE(M))*(74*XNJ(1) +
     >   513*XNJ(2) - 513*XNJ(3) - 74*XNJ(4)))/313600 -
     >   (3*7**(0.5D0)*DA(I,J,M)*(84*XNI(4,J) + 441*XNI(3,J) +
     >   441*XNI(2,J) + 84*XNI(1,J)))/313600
      Q2(5,17) = Q2(5,17) +
     >   (3**(0.5D0)*DA(I,J,M)*SIGN(1.0,DU(M))*(173*XNI(4,J) +
     >   756*XNI(3,J) - 756*XNI(2,J) - 173*XNI(1,J)))/57600 -
     >   (3**(0.5D0)*DB(I,M)*(132*XNJ(1) + 693*XNJ(2) + 693*XNJ(3) +
     >   132*XNJ(4)))/57600
      Q2(6,17) = Q2(6,17) +
     >   ((231*XNI(2,J))/8000 - (231*XNI(3,J))/8000 -
     >   (1903*XNI(4,J))/288000 + (1903*XNI(1,J))/288000)*DA(I,J,M) +
     >   ((1903*XNJ(1))/288000 + (231*XNJ(2))/8000 -
     >   (231*XNJ(3))/8000 - (1903*XNJ(4))/288000)*DB(I,M)
      Q2(7,17) = Q2(7,17) +
     >   (15**(0.5D0)*DA(I,J,M)*SIGN(1.0,DU(M))*(173*XNI(4,J) +
     >   756*XNI(3,J) - 756*XNI(2,J) - 173*XNI(1,J)))/96000 -
     >   (15**(0.5D0)*DB(I,M)*(33*XNJ(1) - 198*XNJ(2) - 198*XNJ(3) +
     >   33*XNJ(4)))/96000
      Q2(8,17) = Q2(8,17)
     >   - (21**(0.5D0)*DB(I,M)*(814*XNJ(1) + 5643*XNJ(2) -
     >   5643*XNJ(3) - 814*XNJ(4)))/1568000 - (21**(0.5D0)*DA(I,J,M)
     >   *(1211*XNI(4,J) + 5292*XNI(3,J) - 5292*XNI(2,J) -
     >   1211*XNI(1,J)))/1568000
      Q2(9,17) = Q2(9,17) +
     >   (3*5**(0.5D0)*DB(I,M)*SIGN(1.0,DE(M))*(4*XNJ(1) +
     >   21*XNJ(2) + 21*XNJ(3) + 4*XNJ(4)))/6400 + (3*5**(0.5D0)*
     >   DA(I,J,M)*SIGN(1.0,DU(M))*(XNI(4,J) - 6*XNI(3,J) -
     >   6*XNI(2,J) + XNI(1,J)))/6400
      Q2(10,17) = Q2(10,17)
     >   - (15**(0.5D0)*DB(I,M)*SIGN(1.0,DE(M))*(173*XNJ(1) +
     >   756*XNJ(2) - 756*XNJ(3) - 173*XNJ(4)))/96000 -
     >   (15**(0.5D0)*DA(I,J,M)*(33*XNI(4,J) - 198*XNI(3,J) -
     >   198*XNI(2,J) + 33*XNI(1,J)))/96000
      Q2(11,17) = Q2(11,17) +
     >   ((9*XNI(4,J))/6400 - (27*XNI(3,J))/3200 -
     >   (27*XNI(2,J))/3200 + (9*XNI(1,J))/6400)*DA(I,J,M)*
     >   SIGN(1.0,DU(M)) + ((9*XNJ(1))/6400 - (27*XNJ(2))/3200 -
     >   (27*XNJ(3))/3200 + (9*XNJ(4))/6400)*DB(I,M)*SIGN(1.0,DE(M))
      Q2(12,17) = Q2(12,17) +
     >   (-(9*5**(0.5D0)*7**(0.5D0)*(21*XNI(4,J) - 126*XNI(3,J) -
     >   126*XNI(2,J) + 21*XNI(1,J)))/1568000)*DA(I,J,M) +
     >   ((9*5**(0.5D0)*7**(0.5D0)*(74*XNJ(1) + 513*XNJ(2) -
     >   513*XNJ(3) - 74*XNJ(4)))/1568000)*DB(I,M)*SIGN(1.0,DE(M))
      Q2(13,17) = Q2(13,17)
     >   - (3*7**(0.5D0)*DB(I,M)*(84*XNJ(1) + 441*XNJ(2) +
     >   441*XNJ(3) + 84*XNJ(4)))/313600 - (3*7**(0.5D0)*DA(I,J,M)*
     >   SIGN(1.0,DU(M))*(74*XNI(4,J) + 513*XNI(3,J) -
     >   513*XNI(2,J) - 74*XNI(1,J)))/313600
      Q2(14,17) = Q2(14,17) +
     >   (21**(0.5D0)*DB(I,M)*(1211*XNJ(1) + 5292*XNJ(2) -
     >   5292*XNJ(3) - 1211*XNJ(4)))/1568000 + (21**(0.5D0)*
     >   DA(I,J,M)*(814*XNI(4,J) + 5643*XNI(3,J) - 5643*XNI(2,J) -
     >   814*XNI(1,J)))/1568000
      Q2(15,17) = Q2(15,17) +
     >   (-(9*5**(0.5D0)*7**(0.5D0)*(74*XNI(4,J) + 513*XNI(3,J) -
     >   513*XNI(2,J) - 74*XNI(1,J)))/1568000)*DA(I,J,M)*
     >   SIGN(1.0,DU(M)) + (-(9*5**(0.5D0)*7**(0.5D0)*(21*XNJ(1) -
     >   126*XNJ(2) - 126*XNJ(3) + 21*XNJ(4)))/1568000)*DB(I,M)
      Q2(16,17) = Q2(16,17) +
     >   ((999*XNI(4,J))/784000 + (13851*XNI(3,J))/1568000 -
     >   (13851*XNI(2,J))/1568000 - (999*XNI(1,J))/784000)*
     >   DA(I,J,M) + ((13851*XNJ(3))/1568000 - (13851*XNJ(2))/
     >   1568000 - (999*XNJ(1))/784000 + (999*XNJ(4))/784000)*DB(I,M)
      ENDIF
*
      CALL ALSBD(IELEM**2,1,Q2,IER,IELEM**2)
      IF(IER.NE.0) CALL XABORT('SNFT22: SINGULAR MATRIX.')
*
      IF(IELEM.EQ.1)THEN
      CORNERQ(1) =  Q2(01,IELEM**2+1)
      ELSEIF(IELEM.EQ.2)THEN
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

      IF((IELEM.EQ.1).OR.(IELEM.EQ.2).OR.(IELEM.EQ.3)
     >   .OR.(IELEM.EQ.4)) THEN
         IF(LFIXUP.AND.(Q2(1,2).LE.RLOG)) Q2(1,2)=0.0
         DO IEL=1,IELEM
            IF((IND.EQ.1).OR.(IND.EQ.4)) THEN
               XNI(IEL,J)=CORNERQ(IEL*IELEM)
               IF(IND.EQ.1) THEN
                  XNJ(IEL)=CORNERQ(IELEM**2 - IELEM +IEL)
               ELSEIF(IND.EQ.4) THEN
                  XNJ(IEL)=CORNERQ(IEL)
               ENDIF
            ELSEIF((IND.EQ.2).OR.(IND.EQ.3)) THEN
               XNI(IEL,J)=CORNERQ((IELEM*(IEL-1))+1)
               IF(IND.EQ.2) THEN
                  XNJ(IEL)=CORNERQ(IELEM**2 - IELEM +IEL)
               ELSEIF(IND.EQ.3) THEN
                  XNJ(IEL)=CORNERQ(IEL)
               ENDIF
            ENDIF
         ENDDO
         IF(LFIXUP.AND.(XNI(1,J).LE.RLOG)) XNI(1,J)=0.0
         IF(LFIXUP.AND.(XNJ(1).LE.RLOG)) XNJ(1)=0.0
      ENDIF

      DO 135 K=1,NSCT
      DO 130 IEL=1,IELEM**2
      FLUX(IEL,K,I,J)=FLUX(IEL,K,I,J)+2.0D0*W(M)*Q2(IEL,IELEM**2+1)*
     1 PL(K,M)
  130 CONTINUE
  135 CONTINUE
*--
  140 CONTINUE
      DO 150 IEL=1,IELEM
      IOF=(M-1)*IELEM*LX+(I-1)*IELEM+IEL
      FUNKNO(L4+IELEM*LY*NPQ+IOF,IG)=REAL(XNJ(IEL))
  150 CONTINUE
*--
  155 CONTINUE
      DO 165 J=1,LY
      DO 160 IEL=1,IELEM
      IOF=(M-1)*IELEM*LY+(J-1)*IELEM+IEL
      FUNKNO(L4+IOF,IG)=REAL(XNI(IEL,J))
  160 CONTINUE
  165 CONTINUE
      FLUX_G(:,:,:,:,IGEFF)=FLUX_G(:,:,:,:,IGEFF)+FLUX(:,:,:,:)
  170 CONTINUE
  180 CONTINUE
*$OMP END PARALLEL DO
  190 CONTINUE
      DO 200 IGEFF=1,NGEFF
        IF(.NOT.INCONV(IGEFF)) GO TO 200
        IG=INGIND(IGEFF)
        FUNKNO(:L4,IG)=
     1  RESHAPE(REAL(FLUX_G(:IELEM**2,:NSCT,:LX,:LY,IGEFF)), (/ L4 /) )
  200 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(XNI,FLUX_G,FLUX,INDANG)
      RETURN
  400 FORMAT(16H SNFT22: thread=,I8,12H --->(group=,I4,7H angle=,I4,1H))
      END
