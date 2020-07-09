*DECK SNFTH4
      SUBROUTINE SNFTH4(NHEX,LZ,IELEM,ISPLH,SIDE,NMAT,NPQ,NSCT,MAT,VOL,
     1 NCODE,ZCODE,TOTAL,QEXT,DU,DE,DZ,W,MRMZ,DC,DB,DA,PL,FLUX,XNEK,
     2 CONNEC,IZGLOB,CONFROM)
*
*-----------------------------------------------------------------------
*
*Purpose:
* perform one inner iteration for solving SN equations in 3D hexagonal
* geometry for the DISCONTINUOUS GALERKIN method. VOID boundary
* conditions on sides, and, VOID or REFL boundary conditions on top and 
* bottom.
*
*Copyright:
* Copyright (C) 2020 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. A. Calloo
*
*Parameters: input
* NHEX    number of hexagons in X-Y plane.
* LZ      number of meshes along Z axis.
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
* XNEI    X-directed SN boundary fluxes.
* XNEJ    Y-directed SN boundary fluxes.
* IZGLOB  hexagon sweep order depending on direction
* CONNEC  connectivity matrix for flux swapping -- which lozenges is the
*         lozenge under consideration connected to; in order to pass the
*         flux along. This is dependent on direction
* CONFROM matrix for incoming flux -- which lozenges are feeding into
*         the lozenge under consideration. This is dependent on
*         direction
*
*Parameters: output
* FLUX    Legendre components of the flux.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NHEX, LZ, IELEM, ISPLH, NMAT, NPQ, NSCT,IZGLOB(NHEX,6),
     1 MAT(ISPLH,ISPLH,3,NHEX,LZ),NCODE(6),MRMZ(NPQ),
     2 CONNEC(3,(NHEX*3)*2,6),CONFROM(2,3,6)
      REAL SIDE,VOL(ISPLH,ISPLH,3,NHEX,LZ), ZCODE(6), TOTAL(0:NMAT),
     1 QEXT(IELEM**3,NSCT,ISPLH,ISPLH,3,NHEX,LZ), DU(NPQ), DE(NPQ),
     2 DZ(NPQ),W(NPQ),DC(ISPLH*ISPLH*3*NHEX,1,NPQ),
     3 DB(ISPLH*ISPLH*3*NHEX,LZ,NPQ), DA(1,LZ,NPQ),
     4 PL(NSCT,NPQ), FLUX(IELEM**3,NSCT,ISPLH,ISPLH,3,NHEX,LZ),
     5 XNEK(IELEM,IELEM,ISPLH,ISPLH,3,NHEX,NPQ)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(RLOG=1.0E-8,PI=3.141592654)
      DOUBLE PRECISION Q(IELEM**3), Q2(IELEM**3,(IELEM**3)+1),
     1   CORNERQ(IELEM**2,IELEM) 
      INTEGER :: ILOZSWP(3,6), IFROMI, IFROMJ
      REAL :: JAC(2,2,3), MUH, ETAH, XIH, AAA, BBB, CCC, DDD, MUHTEMP,
     1 ETAHTEMP
      DOUBLE PRECISION :: THETA, XNI(IELEM,IELEM,ISPLH,LZ), 
     1 XNJ(IELEM,IELEM,LZ),XNK(IELEM,IELEM)  
      DOUBLE PRECISION CONST0, CONST1, CONST2, CONST3, CONST4, CONST5,
     1 CONST6, CONST7, CONST8, CONST9, CONST10, CONST11
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: BFLUX
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(BFLUX(2,NHEX*3,ISPLH,IELEM,IELEM,LZ))
*----
*  CONSTRUCT JACOBIAN MATRIX FOR EACH LOZENGE
*----
      ILOZSWP = RESHAPE((/ 3, 2, 1, 3, 1, 2, 1, 3, 2, 1, 2, 3, 2, 1,
     1   3, 2, 3, 1 /), SHAPE(ILOZSWP))
      JAC = RESHAPE((/ 1., -SQRT(3.), 1., SQRT(3.), 2., 0., 1.,
     1    SQRT(3.), 2., 0., -1., SQRT(3.) /), SHAPE(JAC))
      JAC = (SIDE/2.)*JAC
*----
*  DEFINITION OF CONSTANTS.
*----
      CONST0 =  27648E05
      CONST1 = 150528E05
      CONST2 =  50176E05
      CONST3 = 351232E05
      CONST4 =   4608E06
      CONST5 = 225792E06
      CONST6 =  25088E06
      CONST7 = 175616E06
      CONST8 =8605184E06
      CONST9 =  37632E05
      CONST10=  87808E05
      CONST11= 614656E05
*----
*  PARAMETER VALIDATION.
*----
      IF((IELEM.LT.2).OR.(IELEM.GT.4))
     1   CALL XABORT('SNFTH4: INVALID IELEM (DIAM) VALUE. '
     2   //'CHECK INPUT DATA FILE.')
*----
*  MAIN LOOP OVER SN ANGLES.
*----
      CALL XDRSET(FLUX,IELEM*IELEM*IELEM*NHEX*(3*ISPLH**2)*NSCT*LZ,0.0)

      DO 210 M=1,NPQ
      WEIGHT=W(M)
      VU=DU(M)
      VE=DE(M)
      VZ=DZ(M)

      IF (NCODE(5).NE.4) THEN
         M1=MRMZ(M)
         DO 61 IEL=1,IELEM
         DO 60 JEL=1,IELEM
         IF (VZ.GT.0.0) THEN
           DO  73 I =1,NHEX
           DO  72 J =1,3
           DO  71 IL=1,ISPLH
           DO  70 JL=1,ISPLH
           E1=XNEK(IEL,JEL,IL,JL,J,I,M)
           XNEK(IEL,JEL,IL,JL,J,I,M)=XNEK(IEL,JEL,IL,JL,J,I,M1)
           XNEK(IEL,JEL,IL,JL,J,I,M1)=E1
  70       CONTINUE
  71       CONTINUE
  72       CONTINUE
  73       CONTINUE
         ENDIF
  60  CONTINUE
  61  CONTINUE
      ENDIF

      IF(WEIGHT.EQ.0.0) GO TO 210
*---- CALCULATE SWEEP ANGLE IN X-Y PLANE
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
*---- UNFOLD OCTANTS
      IND=0
      IND2=0
      IF(VZ.GE.0.0)THEN
         IF((THETA.GT.0.0).AND.(THETA.LT.(PI/3.)))THEN
            IND=1
            IND2=1
         ELSEIF((THETA.GT.(PI/3.)).AND.(THETA.LT.(2.*PI/3.)))THEN
            IND=2
            IND2=2
         ELSEIF((THETA.GT.(2.*PI/3.)).AND.(THETA.LT.(PI)))THEN
            IND=3
            IND2=3
         ELSEIF((THETA.GT.(PI)).AND.(THETA.LT.(4.*PI/3.)))THEN
            IND=4
            IND2=4
         ELSEIF((THETA.GT.(4.*PI/3.)).AND.(THETA.LT.(5.*PI/3.)))THEN
            IND=5
            IND2=5
         ELSEIF((THETA.GT.(5.*PI/3.)).AND.(THETA.LT.(2.*PI)))THEN
            IND=6
            IND2=6
         ENDIF
      ELSEIF(VZ.LT.0.0)THEN
         IF((THETA.GT.0.0).AND.(THETA.LT.(PI/3.)))THEN
            IND=7
            IND2=1
         ELSEIF((THETA.GT.(PI/3.)).AND.(THETA.LT.(2.*PI/3.)))THEN
            IND=8
            IND2=2
         ELSEIF((THETA.GT.(2.*PI/3.)).AND.(THETA.LT.(PI)))THEN
            IND=9
            IND2=3
         ELSEIF((THETA.GT.(PI)).AND.(THETA.LT.(4.*PI/3.)))THEN
            IND=10
            IND2=4
         ELSEIF((THETA.GT.(4.*PI/3.)).AND.(THETA.LT.(5.*PI/3.)))THEN
            IND=11
            IND2=5
         ELSEIF((THETA.GT.(5.*PI/3.)).AND.(THETA.LT.(2.*PI)))THEN
            IND=12
            IND2=6
         ENDIF
      ENDIF

      BFLUX(:,:,:,:,:,:) = 0.0D0
*
*----
*  LOOP OVER X-Y PLANE AND Z-DIRECTED AXIS.
*----
* SECOND loop over hexagons
      DO 190 II=1,NHEX
      I=IZGLOB(II,IND2)

* THIRD loop over lozenges
      DO 180 JJ=1,3
      J=ILOZSWP(JJ,IND2)
*
      AAA = JAC(1,1,J)
      BBB = JAC(1,2,J)
      CCC = JAC(2,1,J)
      DDD = JAC(2,2,J)
*
      IHEXI  = CONNEC(1,((I-1)*3*2) + ((J-1)*2) +1,IND2)
      ILOZI  = CONNEC(2,((I-1)*3*2) + ((J-1)*2) +1,IND2)
      ISIDEI = CONNEC(3,((I-1)*3*2) + ((J-1)*2) +1,IND2)
      IHEXJ  = CONNEC(1,((I-1)*3*2) + ((J-1)*2) +2,IND2)
      ILOZJ  = CONNEC(2,((I-1)*3*2) + ((J-1)*2) +2,IND2)
      ISIDEJ = CONNEC(3,((I-1)*3*2) + ((J-1)*2) +2,IND2)
      IFROMI = CONFROM(1,J,IND2)
      IFROMJ = CONFROM(2,J,IND2)
      INDEXI = ((IHEXI-1)*3)+ILOZI
      INDEXJ = ((IHEXJ-1)*3)+ILOZJ
*
      DO 170 IL=1,ISPLH
      I2=IL
      IF(CONFROM(1,J,IND2).EQ.3) I2=ISPLH+1-IL

      DO 160 JL=1,ISPLH
      J2=JL
      IF(CONFROM(2,J,IND2).EQ.4) J2=ISPLH+1-JL

      DO IEL=1,IELEM
      DO JEL=1,IELEM
         IF(IND.LT.7) THEN  
            XNK(IEL,JEL)=XNEK(IEL,JEL,I2,J2,J,I,M)*ZCODE(5)
         ELSE
            XNK(IEL,JEL)=XNEK(IEL,JEL,I2,J2,J,I,M)*ZCODE(6)
         ENDIF
      ENDDO 
      ENDDO

      DO 150 IZ=1,LZ
      K2= IZ
      IF(IND.GE.7) K2=LZ+1-IZ

      IF(JL.EQ.1)THEN
         DO IEL=1,IELEM
            DO JEL=1,IELEM
               XNJ(IEL,JEL,K2) = BFLUX(2,((I-1)*3)+J,I2,IEL,JEL,K2)
            ENDDO
         ENDDO
      ENDIF

      IF(IL.EQ.1)THEN
         DO IEL=1,IELEM
            DO JEL=1,IELEM
               XNI(IEL,JEL,J2,K2) = BFLUX(1,((I-1)*3)+J,J2,IEL,JEL,K2)
            ENDDO
         ENDDO
      ENDIF
*
      MUHTEMP  =  DA(1,K2,M)
      ETAHTEMP =  DB(1,K2,M)
      MUH  = (MUHTEMP*DDD) - (ETAHTEMP*BBB)
      ETAH = (-MUHTEMP*CCC) + (ETAHTEMP*AAA)
      XIH  = DC(1,1,M)
*
      IF(MAT(I2,J2,J,I,K2).EQ.0) GO TO 180
*
      DO 111 IEL=1,IELEM**3
      Q(IEL)=0.0
      DO 110 K=1,NSCT
      Q(IEL)=Q(IEL)+QEXT(IEL,K,I2,J2,J,I,K2)*PL(K,M)/(4.0*PI)
  110 CONTINUE
  111 CONTINUE
*
      VT=VOL(I2,J2,J,I,K2)*TOTAL(MAT(I2,J2,J,I,K2))
      CALL XDDSET(Q2,(IELEM**3)*((IELEM**3)+1),0.0D0)
*
      IF(IELEM.EQ.2) THEN
      Q2(1,1) = VT/8.0D0 + ABS(ETAH)/8.0D0 + 
     >   ABS(MUH)/8.0D0 + 
     >   ABS(XIH)/8.0D0
      Q2(1,2) = (3.0D0**(0.5D0)*MUH)/8.0D0
      Q2(1,3) = (3.0D0**(0.5D0)*ETAH)/8.0D0
      Q2(1,4) = 0.0D0
      Q2(1,5) = (3.0D0**(0.5D0)*XIH)/8.0D0
      Q2(1,6) = 0.0D0
      Q2(1,7) = 0.0D0
      Q2(1,8) = 0.0D0

      Q2(2,1) = -(3.0D0**(0.5D0)*MUH)/24.0D0
      Q2(2,2) = VT/24.0D0 + ABS(ETAH)/24.0D0 + 
     >   ABS(MUH)/8.0D0 + 
     >   ABS(XIH)/24.0D0
      Q2(2,3) = 0.0D0
      Q2(2,4) = (3.0D0**(0.5D0)*ETAH)/24.0D0
      Q2(2,5) = 0.0D0
      Q2(2,6) = (3.0D0**(0.5D0)*XIH)/24.0D0
      Q2(2,7) = 0.0D0
      Q2(2,8) = 0.0D0

      Q2(3,1) = -(3.0D0**(0.5D0)*ETAH)/24.0D0
      Q2(3,2) = 0.0D0
      Q2(3,3) = VT/24.0D0 + ABS(ETAH)/8.0D0 + 
     >   ABS(MUH)/24.0D0 + 
     >   ABS(XIH)/24.0D0
      Q2(3,4) = (3.0D0**(0.5D0)*MUH)/24.0D0
      Q2(3,5) = 0.0D0
      Q2(3,6) = 0.0D0
      Q2(3,7) = (3.0D0**(0.5D0)*XIH)/24.0D0
      Q2(3,8) = 0.0D0

      Q2(4,1) = 0.0D0
      Q2(4,2) = -(3.0D0**(0.5D0)*ETAH)/72.0D0
      Q2(4,3) = -(3.0D0**(0.5D0)*MUH)/72.0D0
      Q2(4,4) = VT/72.0D0 + ABS(ETAH)/24.0D0 + 
     >   ABS(MUH)/24.0D0 + 
     >   ABS(XIH)/72.0D0
      Q2(4,5) = 0.0D0
      Q2(4,6) = 0.0D0
      Q2(4,7) = 0.0D0
      Q2(4,8) = (3.0D0**(0.5D0)*XIH)/72.0D0

      Q2(5,1) = -(3.0D0**(0.5D0)*XIH)/24.0D0
      Q2(5,2) = 0.0D0
      Q2(5,3) = 0.0D0
      Q2(5,4) = 0.0D0
      Q2(5,5) = VT/24.0D0 + ABS(ETAH)/24.0D0 + 
     >   ABS(MUH)/24.0D0 + 
     >   ABS(XIH)/8.0D0
      Q2(5,6) = (3.0D0**(0.5D0)*MUH)/24.0D0
      Q2(5,7) = (3.0D0**(0.5D0)*ETAH)/24.0D0
      Q2(5,8) = 0.0D0

      Q2(6,1) = 0.0D0
      Q2(6,2) = -(3.0D0**(0.5D0)*XIH)/72.0D0
      Q2(6,3) = 0.0D0
      Q2(6,4) = 0.0D0
      Q2(6,5) = -(3.0D0**(0.5D0)*MUH)/72.0D0
      Q2(6,6) = VT/72.0D0 + ABS(ETAH)/72.0D0 + 
     >   ABS(MUH)/24.0D0 + 
     >   ABS(XIH)/24.0D0
      Q2(6,7) = 0.0D0
      Q2(6,8) = (3.0D0**(0.5D0)*ETAH)/72.0D0

      Q2(7,1) = 0.0D0
      Q2(7,2) = 0.0D0
      Q2(7,3) = -(3.0D0**(0.5D0)*XIH)/72.0D0
      Q2(7,4) = 0.0D0
      Q2(7,5) = -(3.0D0**(0.5D0)*ETAH)/72.0D0
      Q2(7,6) = 0.0D0
      Q2(7,7) = VT/72.0D0 + ABS(ETAH)/24.0D0 + 
     >   ABS(MUH)/72.0D0 + 
     >   ABS(XIH)/24.0D0
      Q2(7,8) = (3.0D0**(0.5D0)*MUH)/72.0D0

      Q2(8,1) = 0.0D0
      Q2(8,2) = 0.0D0
      Q2(8,3) = 0.0D0
      Q2(8,4) = -(3.0D0**(0.5D0)*XIH)/216.0D0
      Q2(8,5) = 0.0D0
      Q2(8,6) = -(3.0D0**(0.5D0)*ETAH)/216.0D0
      Q2(8,7) = -(3.0D0**(0.5D0)*MUH)/216.0D0
      Q2(8,8) = VT/216.0D0 + ABS(ETAH)/72.0D0 + 
     >   ABS(MUH)/72.0D0 + 
     >   ABS(XIH)/72.0D0

      Q2(1,9) = ((Q(01))/8.0D0)*VOL(I2,J2,J,I,K2) 
      Q2(2,9) = ((Q(02))/24.0D0)*VOL(I2,J2,J,I,K2) 
      Q2(3,9) = ((Q(03))/24.0D0)*VOL(I2,J2,J,I,K2) 
      Q2(4,9) = ((Q(04))/72.0D0)*VOL(I2,J2,J,I,K2) 
      Q2(5,9) = ((Q(05))/24.0D0)*VOL(I2,J2,J,I,K2) 
      Q2(6,9) = ((Q(06))/72.0D0)*VOL(I2,J2,J,I,K2) 
      Q2(7,9) = ((Q(07))/72.0D0)*VOL(I2,J2,J,I,K2) 
      Q2(8,9) = ((Q(08))/216.0D0)*VOL(I2,J2,J,I,K2) 

      Q2(1,9) = Q2(1,9) + 
     >   (XN
     >   I(2,1,J2,K2)/32.0D0 + XNI(1,1,J2,K2)/32.0D0 + 
     >   XNI(2,2,J2,K2)/32.0D0 + 
     >   XNI(1,2,J2,K2)/32.0D0)*MUH*SIGN(1.0,MUH) + 
     >   (XNJ(1,1,K2)/32.0D0 + 
     >   XNJ(2,1,K2)/32.0D0 + XNJ(1,2,K2)/32.0D0 + 
     >   XNJ(2,2,K2)/32.0D0)*ETAH*SIGN(1.0,ETAH) + 
     >   (XNK(1,2)/32.0D0 + XNK(1,1)/32.0D0 + XNK(2,1)/32.0D0 + 
     >   XNK(2,2)/32.0D0)*XIH*SIGN(1.0,XIH)
      Q2(2,9) = Q2(2,9) 
     >   - (3.0D0**(0.5D0)*MUH*(3.0D0*XNI(2,1,J2,K2) 
     >   + 3.0D0*XNI(1,1,J2,K2) + 
     >   3.0D0*XNI(2,2,J2,K2) + 3.0D0*XNI(1,2,J2,K2)))/288.0D0 - 
     >   (3.0D0**(0.5D0)*ETAH*SIGN(1.0,ETAH)*(XNJ(1,1,K2) - 
     >   XNJ(2,1,K2) + XNJ(1,2,K2) - XNJ(2,2,K2)))/288.0D0 - 
     >   (3.0D0**(0.5D0)*XIH*SIGN(1.0,XIH)*(XNK(1,2) + 
     >   XNK(1,1) - XNK(2,1) - XNK(2,2)))/288.0D0
      Q2(3,9) = Q2(3,9) + 
     >   (3.0D0**(0.5D0)*MUH*SIGN(1.0,MUH)*(XNI(2,1,J2,K2) - 
     >   XNI(1,1,J2,K2) + XNI(2,2,J2,K2) - XNI(1,2,J2,K2)))/288.0D0 - 
     >   (3.0D0**(0.5D0)*ETAH*(3.0D0*XNJ(1,1,K2) 
     >   + 3.0D0*XNJ(2,1,K2) + 
     >   3.0D0*XNJ(1,2,K2) + 3.0D0*XNJ(2,2,K2)))/288.0D0 + 
     >   (3.0D0**(0.5D0)*XIH*SIGN(1.0,XIH)*(XNK(1,2) - 
     >   XNK(1,1) - XNK(2,1) + XNK(2,2)))/288.0D0
      Q2(4,9) = Q2(4,9) + 
     >   (XNI(1,1,J2,K2)/288.0D0 - XNI(2,1,J2,K2)/288.0D0 - 
     >   XNI(2,2,J2,K2)/288.0D0 + 
     >   XNI(1,2,J2,K2)/288.0D0)*MUH + (XNJ(1,1,K2)/288.0D0 - 
     >   XNJ(2,1,K2)/288.0D0 + XNJ(1,2,K2)/288.0D0 - 
     >   XNJ(2,2,K2)/288.0D0)*ETAH + (XNK(1,1)/864.0D0 
     >   - XNK(1,2)/864.0D0 - 
     >   XNK(2,1)/864.0D0 + XNK(2,2)/864.0D0)*XIH*SIGN(1.0,XIH)
      Q2(5,9) = Q2(5,9) 
     >   - (3.0D0**(0.5D0)*XIH*(3.0D0*XNK(1,2) + 3.0D0*XNK(1,1) + 
     >   3.0D0*XNK(2,1) + 3.0D0*XNK(2,2)))/288.0D0 - 
     >   (3.0D0**(0.5D0)*ETAH*SIGN(1.0,ETAH)*(XNJ(1,1,K2) + 
     >   XNJ(2,1,K2) - XNJ(1,2,K2) - XNJ(2,2,K2)))/288.0D0 - 
     >   (3.0D0**(0.5D0)*MUH*SIGN(1.0,MUH)*(XNI(2,1,J2,K2) + 
     >   XNI(1,1,J2,K2) - XNI(2,2,J2,K2) - XNI(1,2,J2,K2)))/288.0D0
      Q2(6,9) = Q2(6,9) + 
     >   (XNI(2,1,J2,K2)/288.0D0 + XNI(1,1,J2,K2)/288.0D0 
     >   - XNI(2,2,J2,K2)/288.0D0 - 
     >   XNI(1,2,J2,K2)/288.0D0)*MUH + (XNJ(1,1,K2)/864.0D0 - 
     >   XNJ(2,1,K2)/864.0D0 - XNJ(1,2,K2)/864.0D0 + 
     >   XNJ(2,2,K2)/864.0D0)*ETAH*SIGN(1.0,ETAH) 
     >   + (XNK(1,2)/288.0D0 + 
     >   XNK(1,1)/288.0D0 - XNK(2,1)/288.0D0 
     >   - XNK(2,2)/288.0D0)*XIH 
      Q2(7,9) = Q2(7,9) + 
     >   (XNI(1,1,J2,K2)/864.0D0 - XNI(2,1,J2,K2)/864.0D0 
     >   + XNI(2,2,J2,K2)/864.0D0 - 
     >   XNI(1,2,J2,K2)/864.0D0)*MUH*SIGN(1.0,MUH) + 
     >   (XNJ(1,1,K2)/288.0D0 + XNJ(2,1,K2)/288.0D0 - 
     >   XNJ(1,2,K2)/288.0D0 - 
     >   XNJ(2,2,K2)/288.0D0)*ETAH + (XNK(1,1)/288.0D0 
     >   - XNK(1,2)/288.0D0 + 
     >   XNK(2,1)/288.0D0 - XNK(2,2)/288.0D0)*XIH 
      Q2(8,9) = Q2(8,9) + 
     >   (3.0D0**(0.5D0)*MUH*(XNI(2,1,J2,K2) - XNI(1,1,J2,K2) - 
     >   XNI(2,2,J2,K2) + XNI(1,2,J2,K2)))/2592.0D0 - 
     >   (3.0D0**(0.5D0)*ETAH*(XNJ(1,1,K2) - XNJ(2,1,K2) - 
     >   XNJ(1,2,K2) + XNJ(2,2,K2)))/2592.0D0 + 
     >   (3.0D0**(0.5D0)*XIH*(XNK(1,2) - XNK(1,1) + XNK(2,1) - 
     >   XNK(2,2)))/2592.0D0
*** ---------------------------------------------------------------- ***
*** ---------------------------------------------------------------- ***
*** ---------------------------------------------------------------- ***
      ELSE IF(IELEM.EQ.3) THEN
*** ---------------------------------------------------------------- ***
*** ---------------------------------------------------------------- ***
*** ---------------------------------------------------------------- ***
      Q2(1,1) = VT/8 + ABS(ETAH)/24 + ABS(MUH)/24 + 
     >   ABS(XIH)/24 
      Q2(1,2) = (5*3**(0.5D0)*MUH)/24 
      Q2(1,3) = -(5**(0.5D0)*(6*VT + 2*ABS(ETAH) - 
     >   15*ABS(MUH) + 2*ABS(XIH)))/360 
      Q2(1,4) = (5*3**(0.5D0)*ETAH)/24 
      Q2(1,5) = 0 
      Q2(1,6) = -(15**(0.5D0)*ETAH)/36 
      Q2(1,7) = -(5**(0.5D0)*(6*VT - 15*ABS(ETAH) + 
     >   2*ABS(MUH) + 2*ABS(XIH)))/360 
      Q2(1,8) = -(15**(0.5D0)*MUH)/36 
      Q2(1,9) = VT/90 - ABS(ETAH)/36 - ABS(MUH)/36 + 
     >   ABS(XIH)/270 
      Q2(1,10) = (5*3**(0.5D0)*XIH)/24 
      Q2(1,11) = 0 
      Q2(1,12) = -(15**(0.5D0)*XIH)/36 
      Q2(1,13) = 0 
      Q2(1,14) = 0 
      Q2(1,15) = 0 
      Q2(1,16) = -(15**(0.5D0)*XIH)/36 
      Q2(1,17) = 0 
      Q2(1,18) = (3**(0.5D0)*XIH)/54 
      Q2(1,19) = -(5**(0.5D0)*(6*VT + 2*ABS(ETAH) + 
     >   2*ABS(MUH) - 15*ABS(XIH)))/360 
      Q2(1,20) = -(15**(0.5D0)*MUH)/36 
      Q2(1,21) = VT/90 + ABS(ETAH)/270 - ABS(MUH)/36 - 
     >   ABS(XIH)/36 
      Q2(1,22) = -(15**(0.5D0)*ETAH)/36 
      Q2(1,23) = 0 
      Q2(1,24) = (3**(0.5D0)*ETAH)/54 
      Q2(1,25) = VT/90 - ABS(ETAH)/36 + ABS(MUH)/270 - 
     >   ABS(XIH)/36 
      Q2(1,26) = (3**(0.5D0)*MUH)/54 
      Q2(1,27) = (5**(0.5D0)*(5*ABS(ETAH) - 2*VT + 
     >   5*ABS(MUH) + 5*ABS(XIH)))/1350 
      
      Q2(2,1) = -(3**(0.5D0)*MUH)/24 
      Q2(2,2) = VT/24 + ABS(ETAH)/72 + ABS(MUH)/8 + 
     >   ABS(XIH)/72 
      Q2(2,3) = (15**(0.5D0)*MUH)/24 
      Q2(2,4) = 0 
      Q2(2,5) = (5*3**(0.5D0)*ETAH)/72 
      Q2(2,6) = 0 
      Q2(2,7) = (15**(0.5D0)*MUH)/180 
      Q2(2,8) = -(5**(0.5D0)*(6*VT - 15*ABS(ETAH) + 
     >   18*ABS(MUH) + 2*ABS(XIH)))/1080 
      Q2(2,9) = -(3**(0.5D0)*MUH)/36 
      Q2(2,10) = 0 
      Q2(2,11) = (5*3**(0.5D0)*XIH)/72 
      Q2(2,12) = 0 
      Q2(2,13) = 0 
      Q2(2,14) = 0 
      Q2(2,15) = 0 
      Q2(2,16) = 0 
      Q2(2,17) = -(15**(0.5D0)*XIH)/108 
      Q2(2,18) = 0 
      Q2(2,19) = (15**(0.5D0)*MUH)/180 
      Q2(2,20) = -(5**(0.5D0)*(6*VT + 2*ABS(ETAH) + 
     >   18*ABS(MUH) - 15*ABS(XIH)))/1080 
      Q2(2,21) = -(3**(0.5D0)*MUH)/36 
      Q2(2,22) = 0 
      Q2(2,23) = -(15**(0.5D0)*ETAH)/108 
      Q2(2,24) = 0 
      Q2(2,25) = -(3**(0.5D0)*MUH)/270 
      Q2(2,26) = VT/270 - ABS(ETAH)/108 + ABS(MUH)/90 - 
     >   ABS(XIH)/108 
      Q2(2,27) = (15**(0.5D0)*MUH)/270 
      
      Q2(3,1) = -(5**(0.5D0)*(3*VT + ABS(ETAH) - 
     >   3*ABS(MUH) + ABS(XIH)))/180 
      Q2(3,2) = -(15**(0.5D0)*MUH)/20 
      Q2(3,3) = VT/30 + ABS(ETAH)/90 + ABS(MUH)/12 + 
     >   ABS(XIH)/90 
      Q2(3,4) = -(15**(0.5D0)*ETAH)/36 
      Q2(3,5) = 0 
      Q2(3,6) = (3**(0.5D0)*ETAH)/18 
      Q2(3,7) = VT/90 - ABS(ETAH)/36 - ABS(MUH)/90 + 
     >   ABS(XIH)/270 
      Q2(3,8) = (3**(0.5D0)*MUH)/30 
      Q2(3,9) = -(5**(0.5D0)*(6*VT - 15*ABS(ETAH) + 
     >   15*ABS(MUH) + 2*ABS(XIH)))/1350 
      Q2(3,10) = -(15**(0.5D0)*XIH)/36 
      Q2(3,11) = 0 
      Q2(3,12) = (3**(0.5D0)*XIH)/18 
      Q2(3,13) = 0 
      Q2(3,14) = 0 
      Q2(3,15) = 0 
      Q2(3,16) = (3**(0.5D0)*XIH)/54 
      Q2(3,17) = 0 
      Q2(3,18) = -(15**(0.5D0)*XIH)/135 
      Q2(3,19) = VT/90 + ABS(ETAH)/270 - ABS(MUH)/90 - 
     >   ABS(XIH)/36 
      Q2(3,20) = (3**(0.5D0)*MUH)/30 
      Q2(3,21) = -(5**(0.5D0)*(6*VT + 2*ABS(ETAH) + 
     >   15*ABS(MUH) - 15*ABS(XIH)))/1350 
      Q2(3,22) = (3**(0.5D0)*ETAH)/54 
      Q2(3,23) = 0 
      Q2(3,24) = -(15**(0.5D0)*ETAH)/135 
      Q2(3,25) = (5**(0.5D0)*(5*ABS(ETAH) - 2*VT + 
     >   2*ABS(MUH) + 5*ABS(XIH)))/1350 
      Q2(3,26) = -(15**(0.5D0)*MUH)/225 
      Q2(3,27) = (2*VT)/675 - ABS(ETAH)/135 + 
     >   ABS(MUH)/135 - ABS(XIH)/135 
      
      Q2(4,1) = -(3**(0.5D0)*ETAH)/24 
      Q2(4,2) = 0 
      Q2(4,3) = (15**(0.5D0)*ETAH)/180 
      Q2(4,4) = VT/24 + ABS(ETAH)/8 + ABS(MUH)/72 + 
     >   ABS(XIH)/72 
      Q2(4,5) = (5*3**(0.5D0)*MUH)/72 
      Q2(4,6) = -(5**(0.5D0)*(6*VT + 18*ABS(ETAH) - 
     >   15*ABS(MUH) + 2*ABS(XIH)))/1080 
      Q2(4,7) = (15**(0.5D0)*ETAH)/24 
      Q2(4,8) = 0 
      Q2(4,9) = -(3**(0.5D0)*ETAH)/36 
      Q2(4,10) = 0 
      Q2(4,11) = 0 
      Q2(4,12) = 0 
      Q2(4,13) = (5*3**(0.5D0)*XIH)/72 
      Q2(4,14) = 0 
      Q2(4,15) = -(15**(0.5D0)*XIH)/108 
      Q2(4,16) = 0 
      Q2(4,17) = 0 
      Q2(4,18) = 0 
      Q2(4,19) = (15**(0.5D0)*ETAH)/180 
      Q2(4,20) = 0 
      Q2(4,21) = -(3**(0.5D0)*ETAH)/270 
      Q2(4,22) = -(5**(0.5D0)*(6*VT + 18*ABS(ETAH) + 
     >   2*ABS(MUH) - 15*ABS(XIH)))/1080 
      Q2(4,23) = -(15**(0.5D0)*MUH)/108 
      Q2(4,24) = VT/270 + ABS(ETAH)/90 - ABS(MUH)/108 - 
     >   ABS(XIH)/108 
      Q2(4,25) = -(3**(0.5D0)*ETAH)/36 
      Q2(4,26) = 0 
      Q2(4,27) = (15**(0.5D0)*ETAH)/270 
      
      Q2(5,1) = 0 
      Q2(5,2) = -(3**(0.5D0)*ETAH)/72 
      Q2(5,3) = 0 
      Q2(5,4) = -(3**(0.5D0)*MUH)/72 
      Q2(5,5) = VT/72 + ABS(ETAH)/24 + ABS(MUH)/24 + 
     >   ABS(XIH)/216 
      Q2(5,6) = (15**(0.5D0)*MUH)/72 
      Q2(5,7) = 0 
      Q2(5,8) = (15**(0.5D0)*ETAH)/72 
      Q2(5,9) = 0 
      Q2(5,10) = 0 
      Q2(5,11) = 0 
      Q2(5,12) = 0 
      Q2(5,13) = 0 
      Q2(5,14) = (5*3**(0.5D0)*XIH)/216 
      Q2(5,15) = 0 
      Q2(5,16) = 0 
      Q2(5,17) = 0 
      Q2(5,18) = 0 
      Q2(5,19) = 0 
      Q2(5,20) = (15**(0.5D0)*ETAH)/540 
      Q2(5,21) = 0 
      Q2(5,22) = (15**(0.5D0)*MUH)/540 
      Q2(5,23) = -(5**(0.5D0)*(2*VT + 6*ABS(ETAH) + 
     >   6*ABS(MUH) - 5*ABS(XIH)))/1080 
      Q2(5,24) = -(3**(0.5D0)*MUH)/108 
      Q2(5,25) = 0 
      Q2(5,26) = -(3**(0.5D0)*ETAH)/108 
      Q2(5,27) = 0 
      
      Q2(6,1) = (15**(0.5D0)*ETAH)/180 
      Q2(6,2) = 0 
      Q2(6,3) = -(3**(0.5D0)*ETAH)/90 
      Q2(6,4) = -(5**(0.5D0)*(3*VT + 9*ABS(ETAH) - 
     >   3*ABS(MUH) + ABS(XIH)))/540 
      Q2(6,5) = -(15**(0.5D0)*MUH)/60 
      Q2(6,6) = VT/90 + ABS(ETAH)/30 + ABS(MUH)/36 + 
     >   ABS(XIH)/270 
      Q2(6,7) = -(3**(0.5D0)*ETAH)/36 
      Q2(6,8) = 0 
      Q2(6,9) = (15**(0.5D0)*ETAH)/90 
      Q2(6,10) = 0 
      Q2(6,11) = 0 
      Q2(6,12) = 0 
      Q2(6,13) = -(15**(0.5D0)*XIH)/108 
      Q2(6,14) = 0 
      Q2(6,15) = (3**(0.5D0)*XIH)/54 
      Q2(6,16) = 0 
      Q2(6,17) = 0 
      Q2(6,18) = 0 
      Q2(6,19) = -(3**(0.5D0)*ETAH)/270 
      Q2(6,20) = 0 
      Q2(6,21) = (15**(0.5D0)*ETAH)/675 
      Q2(6,22) = VT/270 + ABS(ETAH)/90 - ABS(MUH)/270 - 
     >   ABS(XIH)/108 
      Q2(6,23) = (3**(0.5D0)*MUH)/90 
      Q2(6,24) = -(5**(0.5D0)*(2*VT + 6*ABS(ETAH) + 
     >   5*ABS(MUH) - 5*ABS(XIH)))/1350 
      Q2(6,25) = (15**(0.5D0)*ETAH)/270 
      Q2(6,26) = 0 
      Q2(6,27) = -(3**(0.5D0)*ETAH)/135 
      
      Q2(7,1) = -(5**(0.5D0)*(3*VT - 3*ABS(ETAH) + 
     >   ABS(MUH) + ABS(XIH)))/180 
      Q2(7,2) = -(15**(0.5D0)*MUH)/36 
      Q2(7,3) = VT/90 - ABS(ETAH)/90 - ABS(MUH)/36 + 
     >   ABS(XIH)/270 
      Q2(7,4) = -(15**(0.5D0)*ETAH)/20 
      Q2(7,5) = 0 
      Q2(7,6) = (3**(0.5D0)*ETAH)/30 
      Q2(7,7) = VT/30 + ABS(ETAH)/12 + ABS(MUH)/90 + 
     >   ABS(XIH)/90 
      Q2(7,8) = (3**(0.5D0)*MUH)/18 
      Q2(7,9) = -(5**(0.5D0)*(6*VT + 15*ABS(ETAH) - 
     >   15*ABS(MUH) + 2*ABS(XIH)))/1350 
      Q2(7,10) = -(15**(0.5D0)*XIH)/36 
      Q2(7,11) = 0 
      Q2(7,12) = (3**(0.5D0)*XIH)/54 
      Q2(7,13) = 0 
      Q2(7,14) = 0 
      Q2(7,15) = 0 
      Q2(7,16) = (3**(0.5D0)*XIH)/18 
      Q2(7,17) = 0 
      Q2(7,18) = -(15**(0.5D0)*XIH)/135 
      Q2(7,19) = VT/90 - ABS(ETAH)/90 + ABS(MUH)/270 - 
     >   ABS(XIH)/36 
      Q2(7,20) = (3**(0.5D0)*MUH)/54 
      Q2(7,21) = (5**(0.5D0)*(2*ABS(ETAH) - 2*VT + 
     >   5*ABS(MUH) + 5*ABS(XIH)))/1350 
      Q2(7,22) = (3**(0.5D0)*ETAH)/30 
      Q2(7,23) = 0 
      Q2(7,24) = -(15**(0.5D0)*ETAH)/225 
      Q2(7,25) = -(5**(0.5D0)*(6*VT + 15*ABS(ETAH) + 
     >   2*ABS(MUH) - 15*ABS(XIH)))/1350 
      Q2(7,26) = -(15**(0.5D0)*MUH)/135 
      Q2(7,27) = (2*VT)/675 + ABS(ETAH)/135 - 
     >   ABS(MUH)/135 - ABS(XIH)/135 
      
      Q2(8,1) = (15**(0.5D0)*MUH)/180 
      Q2(8,2) = -(5**(0.5D0)*(3*VT - 3*ABS(ETAH) + 
     >   9*ABS(MUH) + ABS(XIH)))/540 
      Q2(8,3) = -(3**(0.5D0)*MUH)/36 
      Q2(8,4) = 0 
      Q2(8,5) = -(15**(0.5D0)*ETAH)/60 
      Q2(8,6) = 0 
      Q2(8,7) = -(3**(0.5D0)*MUH)/90 
      Q2(8,8) = VT/90 + ABS(ETAH)/36 + ABS(MUH)/30 + 
     >   ABS(XIH)/270 
      Q2(8,9) = (15**(0.5D0)*MUH)/90 
      Q2(8,10) = 0 
      Q2(8,11) = -(15**(0.5D0)*XIH)/108 
      Q2(8,12) = 0 
      Q2(8,13) = 0 
      Q2(8,14) = 0 
      Q2(8,15) = 0 
      Q2(8,16) = 0 
      Q2(8,17) = (3**(0.5D0)*XIH)/54 
      Q2(8,18) = 0 
      Q2(8,19) = -(3**(0.5D0)*MUH)/270 
      Q2(8,20) = VT/270 - ABS(ETAH)/270 + ABS(MUH)/90 - 
     >   ABS(XIH)/108 
      Q2(8,21) = (15**(0.5D0)*MUH)/270 
      Q2(8,22) = 0 
      Q2(8,23) = (3**(0.5D0)*ETAH)/90 
      Q2(8,24) = 0 
      Q2(8,25) = (15**(0.5D0)*MUH)/675 
      Q2(8,26) = -(5**(0.5D0)*(2*VT + 5*ABS(ETAH) + 
     >   6*ABS(MUH) - 5*ABS(XIH)))/1350 
      Q2(8,27) = -(3**(0.5D0)*MUH)/135 
      
      Q2(9,1) = VT/90 - ABS(ETAH)/90 - ABS(MUH)/90 + 
     >   ABS(XIH)/270 
      Q2(9,2) = (3**(0.5D0)*MUH)/30 
      Q2(9,3) = -(5**(0.5D0)*(6*VT - 6*ABS(ETAH) + 
     >   15*ABS(MUH) + 2*ABS(XIH)))/1350 
      Q2(9,4) = (3**(0.5D0)*ETAH)/30 
      Q2(9,5) = 0 
      Q2(9,6) = -(15**(0.5D0)*ETAH)/75 
      Q2(9,7) = -(5**(0.5D0)*(6*VT + 15*ABS(ETAH) - 
     >   6*ABS(MUH) + 2*ABS(XIH)))/1350 
      Q2(9,8) = -(15**(0.5D0)*MUH)/75 
      Q2(9,9) = (2*VT)/225 + ABS(ETAH)/45 + ABS(MUH)/45 + 
     >   (2*ABS(XIH))/675 
      Q2(9,10) = (3**(0.5D0)*XIH)/54 
      Q2(9,11) = 0 
      Q2(9,12) = -(15**(0.5D0)*XIH)/135 
      Q2(9,13) = 0 
      Q2(9,14) = 0 
      Q2(9,15) = 0 
      Q2(9,16) = -(15**(0.5D0)*XIH)/135 
      Q2(9,17) = 0 
      Q2(9,18) = (2*3**(0.5D0)*XIH)/135 
      Q2(9,19) = (5**(0.5D0)*(2*ABS(ETAH) - 2*VT + 
     >   2*ABS(MUH) + 5*ABS(XIH)))/1350 
      Q2(9,20) = -(15**(0.5D0)*MUH)/225 
      Q2(9,21) = (2*VT)/675 - (2*ABS(ETAH))/675 + 
     >   ABS(MUH)/135 - ABS(XIH)/135 
      Q2(9,22) = -(15**(0.5D0)*ETAH)/225 
      Q2(9,23) = 0 
      Q2(9,24) = (2*3**(0.5D0)*ETAH)/225 
      Q2(9,25) = (2*VT)/675 + ABS(ETAH)/135 - 
     >   (2*ABS(MUH))/675 - ABS(XIH)/135 
      Q2(9,26) = (2*3**(0.5D0)*MUH)/225 
      Q2(9,27) = -(2*5**(0.5D0)*(2*VT + 5*ABS(ETAH) + 
     >   5*ABS(MUH) - 5*ABS(XIH)))/3375 
      
      Q2(10,1) = -(3**(0.5D0)*XIH)/24 
      Q2(10,2) = 0 
      Q2(10,3) = (15**(0.5D0)*XIH)/180 
      Q2(10,4) = 0 
      Q2(10,5) = 0 
      Q2(10,6) = 0 
      Q2(10,7) = (15**(0.5D0)*XIH)/180 
      Q2(10,8) = 0 
      Q2(10,9) = -(3**(0.5D0)*XIH)/270 
      Q2(10,10) = VT/24 + ABS(ETAH)/72 + ABS(MUH)/72 + 
     >   ABS(XIH)/8 
      Q2(10,11) = (5*3**(0.5D0)*MUH)/72 
      Q2(10,12) = -(5**(0.5D0)*(6*VT + 2*ABS(ETAH) - 
     >   15*ABS(MUH) + 18*ABS(XIH)))/1080 
      Q2(10,13) = (5*3**(0.5D0)*ETAH)/72 
      Q2(10,14) = 0 
      Q2(10,15) = -(15**(0.5D0)*ETAH)/108 
      Q2(10,16) = -(5**(0.5D0)*(6*VT - 15*ABS(ETAH) + 
     >   2*ABS(MUH) + 18*ABS(XIH)))/1080 
      Q2(10,17) = -(15**(0.5D0)*MUH)/108 
      Q2(10,18) = VT/270 - ABS(ETAH)/108 - ABS(MUH)/108 + 
     >   ABS(XIH)/90 
      Q2(10,19) = (15**(0.5D0)*XIH)/24 
      Q2(10,20) = 0 
      Q2(10,21) = -(3**(0.5D0)*XIH)/36 
      Q2(10,22) = 0 
      Q2(10,23) = 0 
      Q2(10,24) = 0 
      Q2(10,25) = -(3**(0.5D0)*XIH)/36 
      Q2(10,26) = 0 
      Q2(10,27) = (15**(0.5D0)*XIH)/270 
      
      Q2(11,1) = 0 
      Q2(11,2) = -(3**(0.5D0)*XIH)/72 
      Q2(11,3) = 0 
      Q2(11,4) = 0 
      Q2(11,5) = 0 
      Q2(11,6) = 0 
      Q2(11,7) = 0 
      Q2(11,8) = (15**(0.5D0)*XIH)/540 
      Q2(11,9) = 0 
      Q2(11,10) = -(3**(0.5D0)*MUH)/72 
      Q2(11,11) = VT/72 + ABS(ETAH)/216 + ABS(MUH)/24 + 
     >   ABS(XIH)/24 
      Q2(11,12) = (15**(0.5D0)*MUH)/72 
      Q2(11,13) = 0 
      Q2(11,14) = (5*3**(0.5D0)*ETAH)/216 
      Q2(11,15) = 0 
      Q2(11,16) = (15**(0.5D0)*MUH)/540 
      Q2(11,17) = -(5**(0.5D0)*(2*VT - 5*ABS(ETAH) + 
     >   6*ABS(MUH) + 6*ABS(XIH)))/1080 
      Q2(11,18) = -(3**(0.5D0)*MUH)/108 
      Q2(11,19) = 0 
      Q2(11,20) = (15**(0.5D0)*XIH)/72 
      Q2(11,21) = 0 
      Q2(11,22) = 0 
      Q2(11,23) = 0 
      Q2(11,24) = 0 
      Q2(11,25) = 0 
      Q2(11,26) = -(3**(0.5D0)*XIH)/108 
      Q2(11,27) = 0 
      
      Q2(12,1) = (15**(0.5D0)*XIH)/180 
      Q2(12,2) = 0 
      Q2(12,3) = -(3**(0.5D0)*XIH)/90 
      Q2(12,4) = 0 
      Q2(12,5) = 0 
      Q2(12,6) = 0 
      Q2(12,7) = -(3**(0.5D0)*XIH)/270 
      Q2(12,8) = 0 
      Q2(12,9) = (15**(0.5D0)*XIH)/675 
      Q2(12,10) = -(5**(0.5D0)*(3*VT + ABS(ETAH) - 
     >   3*ABS(MUH) + 9*ABS(XIH)))/540 
      Q2(12,11) = -(15**(0.5D0)*MUH)/60 
      Q2(12,12) = VT/90 + ABS(ETAH)/270 + ABS(MUH)/36 + 
     >   ABS(XIH)/30 
      Q2(12,13) = -(15**(0.5D0)*ETAH)/108 
      Q2(12,14) = 0 
      Q2(12,15) = (3**(0.5D0)*ETAH)/54 
      Q2(12,16) = VT/270 - ABS(ETAH)/108 - ABS(MUH)/270 + 
     >   ABS(XIH)/90 
      Q2(12,17) = (3**(0.5D0)*MUH)/90 
      Q2(12,18) = -(5**(0.5D0)*(2*VT - 5*ABS(ETAH) + 
     >   5*ABS(MUH) + 6*ABS(XIH)))/1350 
      Q2(12,19) = -(3**(0.5D0)*XIH)/36 
      Q2(12,20) = 0 
      Q2(12,21) = (15**(0.5D0)*XIH)/90 
      Q2(12,22) = 0 
      Q2(12,23) = 0 
      Q2(12,24) = 0 
      Q2(12,25) = (15**(0.5D0)*XIH)/270 
      Q2(12,26) = 0 
      Q2(12,27) = -(3**(0.5D0)*XIH)/135 
      
      Q2(13,1) = 0 
      Q2(13,2) = 0 
      Q2(13,3) = 0 
      Q2(13,4) = -(3**(0.5D0)*XIH)/72 
      Q2(13,5) = 0 
      Q2(13,6) = (15**(0.5D0)*XIH)/540 
      Q2(13,7) = 0 
      Q2(13,8) = 0 
      Q2(13,9) = 0 
      Q2(13,10) = -(3**(0.5D0)*ETAH)/72 
      Q2(13,11) = 0 
      Q2(13,12) = (15**(0.5D0)*ETAH)/540 
      Q2(13,13) = VT/72 + ABS(ETAH)/24 + ABS(MUH)/216 + 
     >   ABS(XIH)/24 
      Q2(13,14) = (5*3**(0.5D0)*MUH)/216 
      Q2(13,15) = -(5**(0.5D0)*(2*VT + 6*ABS(ETAH) - 
     >   5*ABS(MUH) + 6*ABS(XIH)))/1080 
      Q2(13,16) = (15**(0.5D0)*ETAH)/72 
      Q2(13,17) = 0 
      Q2(13,18) = -(3**(0.5D0)*ETAH)/108 
      Q2(13,19) = 0 
      Q2(13,20) = 0 
      Q2(13,21) = 0 
      Q2(13,22) = (15**(0.5D0)*XIH)/72 
      Q2(13,23) = 0 
      Q2(13,24) = -(3**(0.5D0)*XIH)/108 
      Q2(13,25) = 0 
      Q2(13,26) = 0 
      Q2(13,27) = 0 
      
      Q2(14,1) = 0 
      Q2(14,2) = 0 
      Q2(14,3) = 0 
      Q2(14,4) = 0 
      Q2(14,5) = -(3**(0.5D0)*XIH)/216 
      Q2(14,6) = 0 
      Q2(14,7) = 0 
      Q2(14,8) = 0 
      Q2(14,9) = 0 
      Q2(14,10) = 0 
      Q2(14,11) = -(3**(0.5D0)*ETAH)/216 
      Q2(14,12) = 0 
      Q2(14,13) = -(3**(0.5D0)*MUH)/216 
      Q2(14,14) = VT/216 + ABS(ETAH)/72 + ABS(MUH)/72 + 
     >   ABS(XIH)/72 
      Q2(14,15) = (15**(0.5D0)*MUH)/216 
      Q2(14,16) = 0 
      Q2(14,17) = (15**(0.5D0)*ETAH)/216 
      Q2(14,18) = 0 
      Q2(14,19) = 0 
      Q2(14,20) = 0 
      Q2(14,21) = 0 
      Q2(14,22) = 0 
      Q2(14,23) = (15**(0.5D0)*XIH)/216 
      Q2(14,24) = 0 
      Q2(14,25) = 0 
      Q2(14,26) = 0 
      Q2(14,27) = 0 
      
      Q2(15,1) = 0 
      Q2(15,2) = 0 
      Q2(15,3) = 0 
      Q2(15,4) = (15**(0.5D0)*XIH)/540 
      Q2(15,5) = 0 
      Q2(15,6) = -(3**(0.5D0)*XIH)/270 
      Q2(15,7) = 0 
      Q2(15,8) = 0 
      Q2(15,9) = 0 
      Q2(15,10) = (15**(0.5D0)*ETAH)/540 
      Q2(15,11) = 0 
      Q2(15,12) = -(3**(0.5D0)*ETAH)/270 
      Q2(15,13) = -(5**(0.5D0)*(VT + 3*ABS(ETAH) - 
     >   ABS(MUH) + 3*ABS(XIH)))/540 
      Q2(15,14) = -(15**(0.5D0)*MUH)/180 
      Q2(15,15) = VT/270 + ABS(ETAH)/90 + ABS(MUH)/108 + 
     >   ABS(XIH)/90 
      Q2(15,16) = -(3**(0.5D0)*ETAH)/108 
      Q2(15,17) = 0 
      Q2(15,18) = (15**(0.5D0)*ETAH)/270 
      Q2(15,19) = 0 
      Q2(15,20) = 0 
      Q2(15,21) = 0 
      Q2(15,22) = -(3**(0.5D0)*XIH)/108 
      Q2(15,23) = 0 
      Q2(15,24) = (15**(0.5D0)*XIH)/270 
      Q2(15,25) = 0 
      Q2(15,26) = 0 
      Q2(15,27) = 0 
      
      Q2(16,1) = (15**(0.5D0)*XIH)/180 
      Q2(16,2) = 0 
      Q2(16,3) = -(3**(0.5D0)*XIH)/270 
      Q2(16,4) = 0 
      Q2(16,5) = 0 
      Q2(16,6) = 0 
      Q2(16,7) = -(3**(0.5D0)*XIH)/90 
      Q2(16,8) = 0 
      Q2(16,9) = (15**(0.5D0)*XIH)/675 
      Q2(16,10) = -(5**(0.5D0)*(3*VT - 3*ABS(ETAH) + 
     >   ABS(MUH) + 9*ABS(XIH)))/540 
      Q2(16,11) = -(15**(0.5D0)*MUH)/108 
      Q2(16,12) = VT/270 - ABS(ETAH)/270 - ABS(MUH)/108 + 
     >   ABS(XIH)/90 
      Q2(16,13) = -(15**(0.5D0)*ETAH)/60 
      Q2(16,14) = 0 
      Q2(16,15) = (3**(0.5D0)*ETAH)/90 
      Q2(16,16) = VT/90 + ABS(ETAH)/36 + ABS(MUH)/270 + 
     >   ABS(XIH)/30 
      Q2(16,17) = (3**(0.5D0)*MUH)/54 
      Q2(16,18) = -(5**(0.5D0)*(2*VT + 5*ABS(ETAH) - 
     >   5*ABS(MUH) + 6*ABS(XIH)))/1350 
      Q2(16,19) = -(3**(0.5D0)*XIH)/36 
      Q2(16,20) = 0 
      Q2(16,21) = (15**(0.5D0)*XIH)/270 
      Q2(16,22) = 0 
      Q2(16,23) = 0 
      Q2(16,24) = 0 
      Q2(16,25) = (15**(0.5D0)*XIH)/90 
      Q2(16,26) = 0 
      Q2(16,27) = -(3**(0.5D0)*XIH)/135 
      
      Q2(17,1) = 0 
      Q2(17,2) = (15**(0.5D0)*XIH)/540 
      Q2(17,3) = 0 
      Q2(17,4) = 0 
      Q2(17,5) = 0 
      Q2(17,6) = 0 
      Q2(17,7) = 0 
      Q2(17,8) = -(3**(0.5D0)*XIH)/270 
      Q2(17,9) = 0 
      Q2(17,10) = (15**(0.5D0)*MUH)/540 
      Q2(17,11) = -(5**(0.5D0)*(VT - ABS(ETAH) + 
     >   3*ABS(MUH) + 3*ABS(XIH)))/540 
      Q2(17,12) = -(3**(0.5D0)*MUH)/108 
      Q2(17,13) = 0 
      Q2(17,14) = -(15**(0.5D0)*ETAH)/180 
      Q2(17,15) = 0 
      Q2(17,16) = -(3**(0.5D0)*MUH)/270 
      Q2(17,17) = VT/270 + ABS(ETAH)/108 + ABS(MUH)/90 + 
     >   ABS(XIH)/90 
      Q2(17,18) = (15**(0.5D0)*MUH)/270 
      Q2(17,19) = 0 
      Q2(17,20) = -(3**(0.5D0)*XIH)/108 
      Q2(17,21) = 0 
      Q2(17,22) = 0 
      Q2(17,23) = 0 
      Q2(17,24) = 0 
      Q2(17,25) = 0 
      Q2(17,26) = (15**(0.5D0)*XIH)/270 
      Q2(17,27) = 0 
      
      Q2(18,1) = -(3**(0.5D0)*XIH)/270 
      Q2(18,2) = 0 
      Q2(18,3) = (15**(0.5D0)*XIH)/675 
      Q2(18,4) = 0 
      Q2(18,5) = 0 
      Q2(18,6) = 0 
      Q2(18,7) = (15**(0.5D0)*XIH)/675 
      Q2(18,8) = 0 
      Q2(18,9) = -(2*3**(0.5D0)*XIH)/675 
      Q2(18,10) = VT/270 - ABS(ETAH)/270 - ABS(MUH)/270 + 
     >   ABS(XIH)/90 
      Q2(18,11) = (3**(0.5D0)*MUH)/90 
      Q2(18,12) = -(5**(0.5D0)*(2*VT - 2*ABS(ETAH) + 
     >   5*ABS(MUH) + 6*ABS(XIH)))/1350 
      Q2(18,13) = (3**(0.5D0)*ETAH)/90 
      Q2(18,14) = 0 
      Q2(18,15) = -(15**(0.5D0)*ETAH)/225 
      Q2(18,16) = -(5**(0.5D0)*(2*VT + 5*ABS(ETAH) - 
     >   2*ABS(MUH) + 6*ABS(XIH)))/1350 
      Q2(18,17) = -(15**(0.5D0)*MUH)/225 
      Q2(18,18) = (2*VT)/675 + ABS(ETAH)/135 + 
     >   ABS(MUH)/135 + (2*ABS(XIH))/225 
      Q2(18,19) = (15**(0.5D0)*XIH)/270 
      Q2(18,20) = 0 
      Q2(18,21) = -(3**(0.5D0)*XIH)/135 
      Q2(18,22) = 0 
      Q2(18,23) = 0 
      Q2(18,24) = 0 
      Q2(18,25) = -(3**(0.5D0)*XIH)/135 
      Q2(18,26) = 0 
      Q2(18,27) = (2*15**(0.5D0)*XIH)/675 
      
      Q2(19,1) = -(5**(0.5D0)*(3*VT + ABS(ETAH) + 
     >   ABS(MUH) - 3*ABS(XIH)))/180 
      Q2(19,2) = -(15**(0.5D0)*MUH)/36 
      Q2(19,3) = VT/90 + ABS(ETAH)/270 - ABS(MUH)/36 - 
     >   ABS(XIH)/90 
      Q2(19,4) = -(15**(0.5D0)*ETAH)/36 
      Q2(19,5) = 0 
      Q2(19,6) = (3**(0.5D0)*ETAH)/54 
      Q2(19,7) = VT/90 - ABS(ETAH)/36 + ABS(MUH)/270 - 
     >   ABS(XIH)/90 
      Q2(19,8) = (3**(0.5D0)*MUH)/54 
      Q2(19,9) = (5**(0.5D0)*(5*ABS(ETAH) - 2*VT + 
     >   5*ABS(MUH) + 2*ABS(XIH)))/1350 
      Q2(19,10) = -(15**(0.5D0)*XIH)/20 
      Q2(19,11) = 0 
      Q2(19,12) = (3**(0.5D0)*XIH)/30 
      Q2(19,13) = 0 
      Q2(19,14) = 0 
      Q2(19,15) = 0 
      Q2(19,16) = (3**(0.5D0)*XIH)/30 
      Q2(19,17) = 0 
      Q2(19,18) = -(15**(0.5D0)*XIH)/225 
      Q2(19,19) = VT/30 + ABS(ETAH)/90 + ABS(MUH)/90 + 
     >   ABS(XIH)/12 
      Q2(19,20) = (3**(0.5D0)*MUH)/18 
      Q2(19,21) = -(5**(0.5D0)*(6*VT + 2*ABS(ETAH) - 
     >   15*ABS(MUH) + 15*ABS(XIH)))/1350 
      Q2(19,22) = (3**(0.5D0)*ETAH)/18 
      Q2(19,23) = 0 
      Q2(19,24) = -(15**(0.5D0)*ETAH)/135 
      Q2(19,25) = -(5**(0.5D0)*(6*VT - 15*ABS(ETAH) + 
     >   2*ABS(MUH) + 15*ABS(XIH)))/1350 
      Q2(19,26) = -(15**(0.5D0)*MUH)/135 
      Q2(19,27) = (2*VT)/675 - ABS(ETAH)/135 - 
     >   ABS(MUH)/135 + ABS(XIH)/135 
      
      Q2(20,1) = (15**(0.5D0)*MUH)/180 
      Q2(20,2) = -(5**(0.5D0)*(3*VT + ABS(ETAH) + 
     >   9*ABS(MUH) - 3*ABS(XIH)))/540 
      Q2(20,3) = -(3**(0.5D0)*MUH)/36 
      Q2(20,4) = 0 
      Q2(20,5) = -(15**(0.5D0)*ETAH)/108 
      Q2(20,6) = 0 
      Q2(20,7) = -(3**(0.5D0)*MUH)/270 
      Q2(20,8) = VT/270 - ABS(ETAH)/108 + ABS(MUH)/90 - 
     >   ABS(XIH)/270 
      Q2(20,9) = (15**(0.5D0)*MUH)/270 
      Q2(20,10) = 0 
      Q2(20,11) = -(15**(0.5D0)*XIH)/60 
      Q2(20,12) = 0 
      Q2(20,13) = 0 
      Q2(20,14) = 0 
      Q2(20,15) = 0 
      Q2(20,16) = 0 
      Q2(20,17) = (3**(0.5D0)*XIH)/90 
      Q2(20,18) = 0 
      Q2(20,19) = -(3**(0.5D0)*MUH)/90 
      Q2(20,20) = VT/90 + ABS(ETAH)/270 + ABS(MUH)/30 + 
     >   ABS(XIH)/36 
      Q2(20,21) = (15**(0.5D0)*MUH)/90 
      Q2(20,22) = 0 
      Q2(20,23) = (3**(0.5D0)*ETAH)/54 
      Q2(20,24) = 0 
      Q2(20,25) = (15**(0.5D0)*MUH)/675 
      Q2(20,26) = -(5**(0.5D0)*(2*VT - 5*ABS(ETAH) + 
     >   6*ABS(MUH) + 5*ABS(XIH)))/1350 
      Q2(20,27) = -(3**(0.5D0)*MUH)/135 
      
      Q2(21,1) = VT/90 + ABS(ETAH)/270 - ABS(MUH)/90 - 
     >   ABS(XIH)/90 
      Q2(21,2) = (3**(0.5D0)*MUH)/30 
      Q2(21,3) = -(5**(0.5D0)*(6*VT + 2*ABS(ETAH) + 
     >   15*ABS(MUH) - 6*ABS(XIH)))/1350 
      Q2(21,4) = (3**(0.5D0)*ETAH)/54 
      Q2(21,5) = 0 
      Q2(21,6) = -(15**(0.5D0)*ETAH)/135 
      Q2(21,7) = (5**(0.5D0)*(5*ABS(ETAH) - 2*VT + 
     >   2*ABS(MUH) + 2*ABS(XIH)))/1350 
      Q2(21,8) = -(15**(0.5D0)*MUH)/225 
      Q2(21,9) = (2*VT)/675 - ABS(ETAH)/135 + 
     >   ABS(MUH)/135 - (2*ABS(XIH))/675 
      Q2(21,10) = (3**(0.5D0)*XIH)/30 
      Q2(21,11) = 0 
      Q2(21,12) = -(15**(0.5D0)*XIH)/75 
      Q2(21,13) = 0 
      Q2(21,14) = 0 
      Q2(21,15) = 0 
      Q2(21,16) = -(15**(0.5D0)*XIH)/225 
      Q2(21,17) = 0 
      Q2(21,18) = (2*3**(0.5D0)*XIH)/225 
      Q2(21,19) = -(5**(0.5D0)*(6*VT + 2*ABS(ETAH) - 
     >   6*ABS(MUH) + 15*ABS(XIH)))/1350 
      Q2(21,20) = -(15**(0.5D0)*MUH)/75 
      Q2(21,21) = (2*VT)/225 + (2*ABS(ETAH))/675 + 
     >   ABS(MUH)/45 + ABS(XIH)/45 
      Q2(21,22) = -(15**(0.5D0)*ETAH)/135 
      Q2(21,23) = 0 
      Q2(21,24) = (2*3**(0.5D0)*ETAH)/135 
      Q2(21,25) = (2*VT)/675 - ABS(ETAH)/135 - 
     >   (2*ABS(MUH))/675 + ABS(XIH)/135 
      Q2(21,26) = (2*3**(0.5D0)*MUH)/225 
      Q2(21,27) = -(2*5**(0.5D0)*(2*VT - 5*ABS(ETAH) + 
     >   5*ABS(MUH) + 5*ABS(XIH)))/3375 
      
      Q2(22,1) = (15**(0.5D0)*ETAH)/180 
      Q2(22,2) = 0 
      Q2(22,3) = -(3**(0.5D0)*ETAH)/270 
      Q2(22,4) = -(5**(0.5D0)*(3*VT + 9*ABS(ETAH) + 
     >   ABS(MUH) - 3*ABS(XIH)))/540 
      Q2(22,5) = -(15**(0.5D0)*MUH)/108 
      Q2(22,6) = VT/270 + ABS(ETAH)/90 - ABS(MUH)/108 - 
     >   ABS(XIH)/270 
      Q2(22,7) = -(3**(0.5D0)*ETAH)/36 
      Q2(22,8) = 0 
      Q2(22,9) = (15**(0.5D0)*ETAH)/270 
      Q2(22,10) = 0 
      Q2(22,11) = 0 
      Q2(22,12) = 0 
      Q2(22,13) = -(15**(0.5D0)*XIH)/60 
      Q2(22,14) = 0 
      Q2(22,15) = (3**(0.5D0)*XIH)/90 
      Q2(22,16) = 0 
      Q2(22,17) = 0 
      Q2(22,18) = 0 
      Q2(22,19) = -(3**(0.5D0)*ETAH)/90 
      Q2(22,20) = 0 
      Q2(22,21) = (15**(0.5D0)*ETAH)/675 
      Q2(22,22) = VT/90 + ABS(ETAH)/30 + ABS(MUH)/270 + 
     >   ABS(XIH)/36 
      Q2(22,23) = (3**(0.5D0)*MUH)/54 
      Q2(22,24) = -(5**(0.5D0)*(2*VT + 6*ABS(ETAH) - 
     >   5*ABS(MUH) + 5*ABS(XIH)))/1350 
      Q2(22,25) = (15**(0.5D0)*ETAH)/90 
      Q2(22,26) = 0 
      Q2(22,27) = -(3**(0.5D0)*ETAH)/135 
      
      Q2(23,1) = 0 
      Q2(23,2) = (15**(0.5D0)*ETAH)/540 
      Q2(23,3) = 0 
      Q2(23,4) = (15**(0.5D0)*MUH)/540 
      Q2(23,5) = -(5**(0.5D0)*(VT + 3*ABS(ETAH) + 
     >   3*ABS(MUH) - ABS(XIH)))/540 
      Q2(23,6) = -(3**(0.5D0)*MUH)/108 
      Q2(23,7) = 0 
      Q2(23,8) = -(3**(0.5D0)*ETAH)/108 
      Q2(23,9) = 0 
      Q2(23,10) = 0 
      Q2(23,11) = 0 
      Q2(23,12) = 0 
      Q2(23,13) = 0 
      Q2(23,14) = -(15**(0.5D0)*XIH)/180 
      Q2(23,15) = 0 
      Q2(23,16) = 0 
      Q2(23,17) = 0 
      Q2(23,18) = 0 
      Q2(23,19) = 0 
      Q2(23,20) = -(3**(0.5D0)*ETAH)/270 
      Q2(23,21) = 0 
      Q2(23,22) = -(3**(0.5D0)*MUH)/270 
      Q2(23,23) = VT/270 + ABS(ETAH)/90 + ABS(MUH)/90 + 
     >   ABS(XIH)/108 
      Q2(23,24) = (15**(0.5D0)*MUH)/270 
      Q2(23,25) = 0 
      Q2(23,26) = (15**(0.5D0)*ETAH)/270 
      Q2(23,27) = 0 
      
      Q2(24,1) = -(3**(0.5D0)*ETAH)/270 
      Q2(24,2) = 0 
      Q2(24,3) = (15**(0.5D0)*ETAH)/675 
      Q2(24,4) = VT/270 + ABS(ETAH)/90 - ABS(MUH)/270 - 
     >   ABS(XIH)/270 
      Q2(24,5) = (3**(0.5D0)*MUH)/90 
      Q2(24,6) = -(5**(0.5D0)*(2*VT + 6*ABS(ETAH) + 
     >   5*ABS(MUH) - 2*ABS(XIH)))/1350 
      Q2(24,7) = (15**(0.5D0)*ETAH)/270 
      Q2(24,8) = 0 
      Q2(24,9) = -(3**(0.5D0)*ETAH)/135 
      Q2(24,10) = 0 
      Q2(24,11) = 0 
      Q2(24,12) = 0 
      Q2(24,13) = (3**(0.5D0)*XIH)/90 
      Q2(24,14) = 0 
      Q2(24,15) = -(15**(0.5D0)*XIH)/225 
      Q2(24,16) = 0 
      Q2(24,17) = 0 
      Q2(24,18) = 0 
      Q2(24,19) = (15**(0.5D0)*ETAH)/675 
      Q2(24,20) = 0 
      Q2(24,21) = -(2*3**(0.5D0)*ETAH)/675 
      Q2(24,22) = -(5**(0.5D0)*(2*VT + 6*ABS(ETAH) - 
     >   2*ABS(MUH) + 5*ABS(XIH)))/1350 
      Q2(24,23) = -(15**(0.5D0)*MUH)/225 
      Q2(24,24) = (2*VT)/675 + (2*ABS(ETAH))/225 + 
     >   ABS(MUH)/135 + ABS(XIH)/135 
      Q2(24,25) = -(3**(0.5D0)*ETAH)/135 
      Q2(24,26) = 0 
      Q2(24,27) = (2*15**(0.5D0)*ETAH)/675 
      
      Q2(25,1) = VT/90 - ABS(ETAH)/90 + ABS(MUH)/270 - 
     >   ABS(XIH)/90 
      Q2(25,2) = (3**(0.5D0)*MUH)/54 
      Q2(25,3) = (5**(0.5D0)*(2*ABS(ETAH) - 2*VT + 
     >   5*ABS(MUH) + 2*ABS(XIH)))/1350 
      Q2(25,4) = (3**(0.5D0)*ETAH)/30 
      Q2(25,5) = 0 
      Q2(25,6) = -(15**(0.5D0)*ETAH)/225 
      Q2(25,7) = -(5**(0.5D0)*(6*VT + 15*ABS(ETAH) + 
     >   2*ABS(MUH) - 6*ABS(XIH)))/1350 
      Q2(25,8) = -(15**(0.5D0)*MUH)/135 
      Q2(25,9) = (2*VT)/675 + ABS(ETAH)/135 - 
     >   ABS(MUH)/135 - (2*ABS(XIH))/675 
      Q2(25,10) = (3**(0.5D0)*XIH)/30 
      Q2(25,11) = 0 
      Q2(25,12) = -(15**(0.5D0)*XIH)/225 
      Q2(25,13) = 0 
      Q2(25,14) = 0 
      Q2(25,15) = 0 
      Q2(25,16) = -(15**(0.5D0)*XIH)/75 
      Q2(25,17) = 0 
      Q2(25,18) = (2*3**(0.5D0)*XIH)/225 
      Q2(25,19) = -(5**(0.5D0)*(6*VT - 6*ABS(ETAH) + 
     >   2*ABS(MUH) + 15*ABS(XIH)))/1350 
      Q2(25,20) = -(15**(0.5D0)*MUH)/135 
      Q2(25,21) = (2*VT)/675 - (2*ABS(ETAH))/675 - 
     >   ABS(MUH)/135 + ABS(XIH)/135 
      Q2(25,22) = -(15**(0.5D0)*ETAH)/75 
      Q2(25,23) = 0 
      Q2(25,24) = (2*3**(0.5D0)*ETAH)/225 
      Q2(25,25) = (2*VT)/225 + ABS(ETAH)/45 + 
     >   (2*ABS(MUH))/675 + ABS(XIH)/45 
      Q2(25,26) = (2*3**(0.5D0)*MUH)/135 
      Q2(25,27) = -(2*5**(0.5D0)*(2*VT + 5*ABS(ETAH) - 
     >   5*ABS(MUH) + 5*ABS(XIH)))/3375 
      
      Q2(26,1) = -(3**(0.5D0)*MUH)/270 
      Q2(26,2) = VT/270 - ABS(ETAH)/270 + ABS(MUH)/90 - 
     >   ABS(XIH)/270 
      Q2(26,3) = (15**(0.5D0)*MUH)/270 
      Q2(26,4) = 0 
      Q2(26,5) = (3**(0.5D0)*ETAH)/90 
      Q2(26,6) = 0 
      Q2(26,7) = (15**(0.5D0)*MUH)/675 
      Q2(26,8) = -(5**(0.5D0)*(2*VT + 5*ABS(ETAH) + 
     >   6*ABS(MUH) - 2*ABS(XIH)))/1350 
      Q2(26,9) = -(3**(0.5D0)*MUH)/135 
      Q2(26,10) = 0 
      Q2(26,11) = (3**(0.5D0)*XIH)/90 
      Q2(26,12) = 0 
      Q2(26,13) = 0 
      Q2(26,14) = 0 
      Q2(26,15) = 0 
      Q2(26,16) = 0 
      Q2(26,17) = -(15**(0.5D0)*XIH)/225 
      Q2(26,18) = 0 
      Q2(26,19) = (15**(0.5D0)*MUH)/675 
      Q2(26,20) = -(5**(0.5D0)*(2*VT - 2*ABS(ETAH) + 
     >   6*ABS(MUH) + 5*ABS(XIH)))/1350 
      Q2(26,21) = -(3**(0.5D0)*MUH)/135 
      Q2(26,22) = 0 
      Q2(26,23) = -(15**(0.5D0)*ETAH)/225 
      Q2(26,24) = 0 
      Q2(26,25) = -(2*3**(0.5D0)*MUH)/675 
      Q2(26,26) = (2*VT)/675 + ABS(ETAH)/135 + 
     >   (2*ABS(MUH))/225 + ABS(XIH)/135 
      Q2(26,27) = (2*15**(0.5D0)*MUH)/675 
      
      Q2(27,1) = (5**(0.5D0)*(ABS(ETAH) - VT + 
     >   ABS(MUH) + ABS(XIH)))/675 
      Q2(27,2) = -(15**(0.5D0)*MUH)/225 
      Q2(27,3) = (2*VT)/675 - (2*ABS(ETAH))/675 + 
     >   ABS(MUH)/135 - (2*ABS(XIH))/675 
      Q2(27,4) = -(15**(0.5D0)*ETAH)/225 
      Q2(27,5) = 0 
      Q2(27,6) = (2*3**(0.5D0)*ETAH)/225 
      Q2(27,7) = (2*VT)/675 + ABS(ETAH)/135 - 
     >   (2*ABS(MUH))/675 - (2*ABS(XIH))/675 
      Q2(27,8) = (2*3**(0.5D0)*MUH)/225 
      Q2(27,9) = -(2*5**(0.5D0)*(2*VT + 5*ABS(ETAH) + 
     >   5*ABS(MUH) - 2*ABS(XIH)))/3375 
      Q2(27,10) = -(15**(0.5D0)*XIH)/225 
      Q2(27,11) = 0 
      Q2(27,12) = (2*3**(0.5D0)*XIH)/225 
      Q2(27,13) = 0 
      Q2(27,14) = 0 
      Q2(27,15) = 0 
      Q2(27,16) = (2*3**(0.5D0)*XIH)/225 
      Q2(27,17) = 0 
      Q2(27,18) = -(4*15**(0.5D0)*XIH)/1125 
      Q2(27,19) = (2*VT)/675 - (2*ABS(ETAH))/675 - 
     >   (2*ABS(MUH))/675 + ABS(XIH)/135 
      Q2(27,20) = (2*3**(0.5D0)*MUH)/225 
      Q2(27,21) = -(2*5**(0.5D0)*(2*VT - 2*ABS(ETAH) + 
     >   5*ABS(MUH) + 5*ABS(XIH)))/3375 
      Q2(27,22) = (2*3**(0.5D0)*ETAH)/225 
      Q2(27,23) = 0 
      Q2(27,24) = -(4*15**(0.5D0)*ETAH)/1125 
      Q2(27,25) = -(2*5**(0.5D0)*(2*VT + 5*ABS(ETAH) - 
     >   2*ABS(MUH) + 5*ABS(XIH)))/3375 
      Q2(27,26) = -(4*15**(0.5D0)*MUH)/1125 
      Q2(27,27) = (8*VT)/3375 + (4*ABS(ETAH))/675 + 
     >   (4*ABS(MUH))/675 + (4*ABS(XIH))/675 
      
      Q2(1,28) = (((675*Q(01) + 60*Q(09) + 60*Q(21) + 60*Q(25) - 
     >   90*5**(0.5D0)*Q(03) - 90*5**(0.5D0)*Q(07) - 90*5**(0.5D0)*
     >   Q(19) - 8*5**(0.5D0)*Q(27)))/5400)*VOL(I2,J2,J,I,K2)  
      Q2(2,28) = (((45*Q(02) + 4*Q(26) - 6*5**(0.5D0)*Q(08) - 
     >   6*5**(0.5D0)*Q(20)))/1080)*VOL(I2,J2,J,I,K2)  
      Q2(3,28) = (((90*Q(03) + 30*Q(07) + 30*Q(19) + 8*Q(27) - 
     >   45*5**(0.5D0)*Q(01) - 12*5**(0.5D0)*Q(09) - 12*5**(0.5D0)*
     >   Q(21) - 4*5**(0.5D0)*Q(25)))/2700)*VOL(I2,J2,J,I,K2)  
      Q2(4,28) = (((45*Q(04) + 4*Q(24) - 6*5**(0.5D0)*Q(06) - 
     >   6*5**(0.5D0)*Q(22)))/1080)*VOL(I2,J2,J,I,K2)  
      Q2(5,28) = (((15*Q(05) - 2*5**(0.5D0)*Q(23)))/1080)*
     >   VOL(I2,J2,J,I,K2)  
      Q2(6,28) = (((30*Q(06) + 10*Q(22) - 15*5**(0.5D0)*Q(04) - 
     >   4*5**(0.5D0)*Q(24)))/2700)*VOL(I2,J2,J,I,K2)  
      Q2(7,28) = (((30*Q(03) + 90*Q(07) + 30*Q(19) + 8*Q(27) - 
     >   45*5**(0.5D0)*Q(01) - 12*5**(0.5D0)*Q(09) - 4*5**(0.5D0)*
     >   Q(21) - 12*5**(0.5D0)*Q(25)))/2700)*VOL(I2,J2,J,I,K2)  
      Q2(8,28) = (((30*Q(08) + 10*Q(20) - 15*5**(0.5D0)*Q(02) - 
     >   4*5**(0.5D0)*Q(26)))/2700)*VOL(I2,J2,J,I,K2)  
      Q2(9,28) = (((81000*Q(01) + 64800*Q(09) + 21600*Q(21) + 
     >   21600*Q(25) - 32400*5**(0.5D0)*Q(03) - 32400*5**(0.5D0)*
     >   Q(07) - 10800*5**(0.5D0)*Q(19) - 8640*5**(0.5D0)*
     >   Q(27)))/7290000)*VOL(I2,J2,J,I,K2)  
      Q2(10,28) = (((45*Q(10) + 4*Q(18) - 6*5**(0.5D0)*Q(12) - 
     >   6*5**(0.5D0)*Q(16)))/1080)*VOL(I2,J2,J,I,K2)  
      Q2(11,28) = (((15*Q(11) - 2*5**(0.5D0)*Q(17)))/1080)*
     >   VOL(I2,J2,J,I,K2)  
      Q2(12,28) = (((30*Q(12) + 10*Q(16) - 15*5**(0.5D0)*Q(10) - 
     >   4*5**(0.5D0)*Q(18)))/2700)*VOL(I2,J2,J,I,K2)  
      Q2(13,28) = (((15*Q(13) - 2*5**(0.5D0)*Q(15)))/1080)*
     >   VOL(I2,J2,J,I,K2)  
      Q2(14,28) = ((Q(14))/216)*VOL(I2,J2,J,I,K2)  
      Q2(15,28) = (((2*Q(15) - 5**(0.5D0)*Q(13)))/540)*
     >   VOL(I2,J2,J,I,K2)  
      Q2(16,28) = (((10*Q(12) + 30*Q(16) - 15*5**(0.5D0)*Q(10) - 
     >   4*5**(0.5D0)*Q(18)))/2700)*VOL(I2,J2,J,I,K2)  
      Q2(17,28) = (((2*Q(17) - 5**(0.5D0)*Q(11)))/540)*
     >   VOL(I2,J2,J,I,K2)  
      Q2(18,28) = (((5*Q(10) + 4*Q(18) - 2*5**(0.5D0)*Q(12) - 
     >   2*5**(0.5D0)*Q(16)))/1350)*VOL(I2,J2,J,I,K2)  
      Q2(19,28) = (((30*Q(03) + 30*Q(07) + 90*Q(19) + 8*Q(27) - 
     >   45*5**(0.5D0)*Q(01) - 4*5**(0.5D0)*Q(09) - 12*5**(0.5D0)*
     >   Q(21) - 12*5**(0.5D0)*Q(25)))/2700)*VOL(I2,J2,J,I,K2)  
      Q2(20,28) = (((10*Q(08) + 30*Q(20) - 15*5**(0.5D0)*Q(02) - 
     >   4*5**(0.5D0)*Q(26)))/2700)*VOL(I2,J2,J,I,K2)  
      Q2(21,28) = (((81000*Q(01) + 21600*Q(09) + 64800*Q(21) + 
     >   21600*Q(25) - 32400*5**(0.5D0)*Q(03) - 10800*5**(0.5D0)*
     >   Q(07) - 32400*5**(0.5D0)*Q(19) - 8640*5**(0.5D0)*
     >   Q(27)))/7290000)*VOL(I2,J2,J,I,K2)  
      Q2(22,28) = (((10*Q(06) + 30*Q(22) - 15*5**(0.5D0)*Q(04) - 
     >   4*5**(0.5D0)*Q(24)))/2700)*VOL(I2,J2,J,I,K2)  
      Q2(23,28) = (((2*Q(23) - 5**(0.5D0)*Q(05)))/540)*
     >   VOL(I2,J2,J,I,K2)  
      Q2(24,28) = (((5*Q(04) + 4*Q(24) - 2*5**(0.5D0)*Q(06) - 
     >   2*5**(0.5D0)*Q(22)))/1350)*VOL(I2,J2,J,I,K2)  
      Q2(25,28) = (((81000*Q(01) + 21600*Q(09) + 21600*Q(21) + 
     >   64800*Q(25) - 10800*5**(0.5D0)*Q(03) - 32400*5**(0.5D0)*
     >   Q(07) - 32400*5**(0.5D0)*Q(19) - 8640*5**(0.5D0)*
     >   Q(27)))/7290000)*VOL(I2,J2,J,I,K2)  
      Q2(26,28) = (((5*Q(02) + 4*Q(26) - 2*5**(0.5D0)*Q(08) - 
     >   2*5**(0.5D0)*Q(20)))/1350)*VOL(I2,J2,J,I,K2)  
      Q2(27,28) = (((10*Q(03) + 10*Q(07) + 10*Q(19) + 8*Q(27) - 
     >   5*5**(0.5D0)*Q(01) - 4*5**(0.5D0)*Q(09) - 4*5**(0.5D0)*
     >   Q(21) - 4*5**(0.5D0)*Q(25)))/3375)*VOL(I2,J2,J,I,K2)  

      Q2(1,28) = Q2(1,28) + 
     >   ((121*XNI(1,3,J2,K2))/194400 + (121*XNI(3,1,J2,K2))/194400 + 
     >   (187*XNI(2,1,J2,K2))/48600 + (121*XNI(1,1,J2,K2))/194400 + 
     >   (187*XNI(3,2,J2,K2))/48600 + (289*XNI(2,2,J2,K2))/12150 + 
     >   (187*XNI(1,2,J2,K2))/48600 + (121*XNI(3,3,J2,K2))/194400 + 
     >   (187*XNI(2,3,J2,K2))/48600)*MUH*SIGN(1.0,MUH) + 
     >   ((121*XNJ(1,3,K2))/194400 + (187*XNJ(2,3,K2))/48600 + 
     >   (121*XNJ(3,3,K2))/194400 + (121*XNJ(1,1,K2))/194400 + 
     >   (187*XNJ(2,1,K2))/48600 + (121*XNJ(3,1,K2))/194400 + 
     >   (187*XNJ(1,2,K2))/48600 + (289*XNJ(2,2,K2))/12150 + 
     >   (187*XNJ(3,2,K2))/48600)*ETAH*SIGN(1.0,ETAH) + 
     >   ((121*XNK(1,3))/194400 + (187*XNK(1,2))/48600 + 
     >   (121*XNK(1,1))/194400 + (187*XNK(2,1))/48600 + 
     >   (121*XNK(3,1))/194400 + (187*XNK(3,2))/48600 + 
     >   (121*XNK(3,3))/194400 + (187*XNK(2,3))/48600 + 
     >   (289*XNK(2,2))/12150)*XIH*SIGN(1.0,XIH) 
      Q2(2,28) = Q2(2,28)
     >   - (3**(0.5D0)*ETAH*SIGN(1.0,ETAH)*(55*XNJ(1,3,K2) - 
     >   55*XNJ(3,3,K2) + 55*XNJ(1,1,K2) - 55*XNJ(3,1,K2) + 
     >   340*XNJ(1,2,K2) - 340*XNJ(3,2,K2)))/194400 - 
     >   (3**(0.5D0)*XIH*SIGN(1.0,XIH)*(55*XNK(1,3) + 
     >   340*XNK(1,2) + 55*XNK(1,1) - 55*XNK(3,1) - 340*XNK(3,2) - 
     >   55*XNK(3,3)))/194400 - (3**(0.5D0)*MUH*
     >   (121*XNI(1,3,J2,K2) + 121*XNI(3,1,J2,K2) + 748*XNI(2,1,J2,K2)+
     >   121*XNI(1,1,J2,K2) + 748*XNI(3,2,J2,K2) + 4624*XNI(2,2,J2,K2)+
     >   748*XNI(1,2,J2,K2) + 121*XNI(3,3,J2,K2) + 
     >   748*XNI(2,3,J2,K2)))/
     >   194400 
      Q2(3,28) = Q2(3,28) + 
     >   (5**(0.5D0)*MUH*SIGN(1.0,MUH)*(121*XNI(1,3,J2,K2) + 
     >   121*XNI(3,1,J2,K2) + 748*XNI(2,1,J2,K2) + 121*XNI(1,1,J2,K2) +
     >   748*XNI(3,2,J2,K2) + 4624*XNI(2,2,J2,K2) + 748*XNI(1,2,J2,K2) +
     >   121*XNI(3,3,J2,K2) + 748*XNI(2,3,J2,K2)))/486000 - 
     >   (5**(0.5D0)*ETAH*SIGN(1.0,ETAH)*(11*XNJ(1,3,K2) + 
     >   308*XNJ(2,3,K2) + 11*XNJ(3,3,K2) + 11*XNJ(1,1,K2) + 
     >   308*XNJ(2,1,K2) + 11*XNJ(3,1,K2) + 68*XNJ(1,2,K2) + 
     >   1904*XNJ(2,2,K2) + 68*XNJ(3,2,K2)))/486000 - (5**(0.5D0)*
     >   XIH*SIGN(1.0,XIH)*(11*XNK(1,3) + 68*XNK(1,2) + 
     >   11*XNK(1,1) + 308*XNK(2,1) + 11*XNK(3,1) + 68*XNK(3,2) + 
     >   11*XNK(3,3) + 308*XNK(2,3) + 1904*XNK(2,2)))/486000 
      Q2(4,28) = Q2(4,28) + 
     >   (3**(0.5D0)*XIH*SIGN(1.0,XIH)*(55*XNK(1,3) - 
     >   55*XNK(1,1) - 340*XNK(2,1) - 55*XNK(3,1) + 55*XNK(3,3) + 
     >   340*XNK(2,3)))/194400 - (3**(0.5D0)*MUH*
     >   SIGN(1.0,MUH)*(55*XNI(1,3,J2,K2) - 55*XNI(3,1,J2,K2) + 
     >   55*XNI(1,1,J2,K2) - 340*XNI(3,2,J2,K2) + 340*XNI(1,2,J2,K2) - 
     >   55*XNI(3,3,J2,K2)))/194400 - (3**(0.5D0)*ETAH*
     >   (121*XNJ(1,3,K2) + 748*XNJ(2,3,K2) + 121*XNJ(3,3,K2) + 
     >   121*XNJ(1,1,K2) + 748*XNJ(2,1,K2) + 121*XNJ(3,1,K2) + 
     >   748*XNJ(1,2,K2) + 4624*XNJ(2,2,K2) + 748*XNJ(3,2,K2)))/194400
      Q2(5,28) = Q2(5,28) + 
     >   ((11*XNI(1,3,J2,K2))/12960 - (11*XNI(3,1,J2,K2))/12960 + 
     >   (11*XNI(1,1,J2,K2))/12960 - (17*XNI(3,2,J2,K2))/3240 + 
     >   (17*XNI(1,2,J2,K2))/3240 - (11*XNI(3,3,J2,K2))/12960)*
     >   MUH + ((11*XNJ(1,3,K2))/12960 - (11*XNJ(3,3,K2))/
     >   12960 + (11*XNJ(1,1,K2))/12960 - (11*XNJ(3,1,K2))/12960 + 
     >   (17*XNJ(1,2,K2))/3240 - (17*XNJ(3,2,K2))/3240)*ETAH + 
     >   (XNK(1,1)/2592 - XNK(1,3)/2592 - XNK(3,1)/2592 + 
     >   XNK(3,3)/2592)*XIH*SIGN(1.0,XIH) 
      Q2(6,28) = Q2(6,28) + 
     >   (15**(0.5D0)*ETAH*(11*XNJ(1,3,K2) + 308*XNJ(2,3,K2) + 
     >   11*XNJ(3,3,K2) + 11*XNJ(1,1,K2) + 308*XNJ(2,1,K2) + 
     >   11*XNJ(3,1,K2) + 68*XNJ(1,2,K2) + 1904*XNJ(2,2,K2) + 
     >   68*XNJ(3,2,K2)))/486000 - (15**(0.5D0)*XIH*
     >   SIGN(1.0,XIH)*(5*XNK(1,3) - 5*XNK(1,1) - 140*XNK(2,1) - 
     >   5*XNK(3,1) + 5*XNK(3,3) + 140*XNK(2,3)))/486000 - 
     >   (15**(0.5D0)*MUH*SIGN(1.0,MUH)*(55*XNI(1,3,J2,K2) - 
     >   55*XNI(3,1,J2,K2) + 55*XNI(1,1,J2,K2) - 340*XNI(3,2,J2,K2) + 
     >   340*XNI(1,2,J2,K2) - 55*XNI(3,3,J2,K2)))/486000 
      Q2(7,28) = Q2(7,28) + 
     >   (5**(0.5D0)*ETAH*SIGN(1.0,ETAH)*(121*XNJ(1,3,K2) + 
     >   748*XNJ(2,3,K2) + 121*XNJ(3,3,K2) + 121*XNJ(1,1,K2) + 
     >   748*XNJ(2,1,K2) + 121*XNJ(3,1,K2) + 748*XNJ(1,2,K2) + 
     >   4624*XNJ(2,2,K2) + 748*XNJ(3,2,K2)))/486000 - (5**(0.5D0)*
     >   MUH*SIGN(1.0,MUH)*(11*XNI(1,3,J2,K2) + 11*
     >   XNI(3,1,J2,K2) + 308*XNI(2,1,J2,K2) + 11*XNI(1,1,J2,K2) + 
     >   68*XNI(3,2,J2,K2) + 1904*XNI(2,2,J2,K2) + 68*XNI(1,2,J2,K2) + 
     >   11*XNI(3,3,J2,K2) + 308*XNI(2,3,J2,K2)))/486000 - (5**(0.5D0)*
     >   XIH*SIGN(1.0,XIH)*(11*XNK(1,3) + 308*XNK(1,2) + 
     >   11*XNK(1,1) + 68*XNK(2,1) + 11*XNK(3,1) + 308*XNK(3,2) + 
     >   11*XNK(3,3) + 68*XNK(2,3) + 1904*XNK(2,2)))/486000 
      Q2(8,28) = Q2(8,28) + 
     >   (15**(0.5D0)*XIH*SIGN(1.0,XIH)*(5*XNK(1,3) + 
     >   140*XNK(1,2) + 5*XNK(1,1) - 5*XNK(3,1) - 140*XNK(3,2) - 
     >   5*XNK(3,3)))/486000 - (15**(0.5D0)*ETAH*
     >   SIGN(1.0,ETAH)*(55*XNJ(1,3,K2) - 55*XNJ(3,3,K2) + 
     >   55*XNJ(1,1,K2) - 55*XNJ(3,1,K2) + 340*XNJ(1,2,K2) - 
     >   340*XNJ(3,2,K2)))/486000 + (15**(0.5D0)*MUH*
     >   (11*XNI(1,3,J2,K2) + 11*XNI(3,1,J2,K2) + 308*XNI(2,1,J2,K2) + 
     >   11*XNI(1,1,J2,K2) + 68*XNI(3,2,J2,K2) + 1904*XNI(2,2,J2,K2) + 
     >   68*XNI(1,2,J2,K2) + 11*XNI(3,3,J2,K2) + 308*XNI(2,3,J2,K2)))/
     >   486000 
      Q2(9,28) = Q2(9,28) + 
     >   XIH*SIGN(1.0,XIH)*(XNK(1,3)/243000 + 
     >   (7*XNK(1,2))/60750 + XNK(1,1)/243000 + (7*XNK(2,1))/60750 + 
     >   XNK(3,1)/243000 + (7*XNK(3,2))/60750 + XNK(3,3)/243000 + 
     >   (7*XNK(2,3))/60750 + (98*XNK(2,2))/30375) - MUH*
     >   SIGN(1.0,MUH)*((11*XNI(1,3,J2,K2))/243000 + (11*
     >   XNI(3,1,J2,K2))/243000 + (77*XNI(2,1,J2,K2))/60750 + (11*
     >   XNI(1,1,J2,K2))/243000 + (17*XNI(3,2,J2,K2))/60750 + (238*
     >   XNI(2,2,J2,K2))/30375 + (17*XNI(1,2,J2,K2))/60750 + (11*
     >   XNI(3,3,J2,K2))/243000 + (77*XNI(2,3,J2,K2))/60750) - 
     >   ETAH*SIGN(1.0,ETAH)*((11*XNJ(1,3,K2))/243000 + 
     >   (77*XNJ(2,3,K2))/60750 + (11*XNJ(3,3,K2))/243000 + 
     >   (11*XNJ(1,1,K2))/243000 + (77*XNJ(2,1,K2))/60750 + 
     >   (11*XNJ(3,1,K2))/243000 + (17*XNJ(1,2,K2))/60750 + 
     >   (238*XNJ(2,2,K2))/30375 + (17*XNJ(3,2,K2))/60750) 
      Q2(10,28) = Q2(10,28) + 
     >   (3**(0.5D0)*ETAH*SIGN(1.0,ETAH)*(55*XNJ(1,3,K2) + 
     >   340*XNJ(2,3,K2) + 55*XNJ(3,3,K2) - 55*XNJ(1,1,K2) - 
     >   340*XNJ(2,1,K2) - 55*XNJ(3,1,K2)))/194400 + (3**(0.5D0)*
     >   MUH*SIGN(1.0,MUH)*(55*XNI(1,3,J2,K2) - 
     >   55*XNI(3,1,J2,K2) - 340*XNI(2,1,J2,K2) - 55*XNI(1,1,J2,K2) + 
     >   55*XNI(3,3,J2,K2) + 340*XNI(2,3,J2,K2)))/194400 - 
     >   (3**(0.5D0)*XIH*(121*XNK(1,3) + 748*XNK(1,2) + 
     >   121*XNK(1,1) + 748*XNK(2,1) + 121*XNK(3,1) + 748*XNK(3,2) + 
     >   121*XNK(3,3) + 748*XNK(2,3) + 4624*XNK(2,2)))/194400 
      Q2(11,28) = Q2(11,28) + 
     >   ((11*XNI(3,1,J2,K2))/12960 - (11*XNI(1,3,J2,K2))/12960 + 
     >   (17*XNI(2,1,J2,K2))/3240 + (11*XNI(1,1,J2,K2))/12960 - 
     >   (11*XNI(3,3,J2,K2))/12960 - (17*XNI(2,3,J2,K2))/3240)*MUH + 
     >   (XNJ(3,3,K2)/2592 - XNJ(1,3,K2)/2592 + XNJ(1,1,K2)/2592 - 
     >   XNJ(3,1,K2)/2592)*ETAH*SIGN(1.0,ETAH) + ((11*XNK(1,3))/
     >   12960 + (17*XNK(1,2))/3240 + (11*XNK(1,1))/12960 - 
     >   (11*XNK(3,1))/12960 - (17*XNK(3,2))/3240 - (11*XNK(3,3))/
     >   12960)*XIH 
      Q2(12,28) = Q2(12,28) + 
     >   (15**(0.5D0)*MUH*SIGN(1.0,MUH)*(55*XNI(1,3,J2,K2) - 
     >   55*XNI(3,1,J2,K2) - 340*XNI(2,1,J2,K2) - 55*XNI(1,1,J2,K2) + 
     >   55*XNI(3,3,J2,K2) + 340*XNI(2,3,J2,K2)))/486000 - 
     >   (15**(0.5D0)*
     >   ETAH*SIGN(1.0,ETAH)*(5*XNJ(1,3,K2) + 140*XNJ(2,3,K2) + 
     >   5*XNJ(3,3,K2) - 5*XNJ(1,1,K2) - 140*XNJ(2,1,K2) - 
     >   5*XNJ(3,1,K2)))/486000 + (15**(0.5D0)*XIH*(11*XNK(1,3) + 
     >   68*XNK(1,2) + 11*XNK(1,1) + 308*XNK(2,1) + 11*XNK(3,1) + 
     >   68*XNK(3,2) + 11*XNK(3,3) + 308*XNK(2,3) + 
     >   1904*XNK(2,2)))/486000 
      Q2(13,28) = Q2(13,28) + 
     >   (XNI(1,1,J2,K2)/2592 - XNI(3,1,J2,K2)/2592 - 
     >   XNI(1,3,J2,K2)/2592 + 
     >   XNI(3,3,J2,K2)/2592)*MUH*SIGN(1.0,MUH) + 
     >   ((11*XNJ(1,1,K2))/12960 - (17*XNJ(2,3,K2))/3240 - 
     >   (11*XNJ(3,3,K2))/12960 - (11*XNJ(1,3,K2))/12960 + 
     >   (17*XNJ(2,1,K2))/3240 + (11*XNJ(3,1,K2))/12960)*ETAH + 
     >   ((11*XNK(1,1))/12960 - (11*XNK(1,3))/12960 + 
     >   (17*XNK(2,1))/3240 + (11*XNK(3,1))/12960 - (11*XNK(3,3))/
     >   12960 - (17*XNK(2,3))/3240)*XIH 
      Q2(14,28) = Q2(14,28) + 
     >   (3**(0.5D0)*ETAH*(XNJ(1,3,K2) - XNJ(3,3,K2) - 
     >   XNJ(1,1,K2) + XNJ(3,1,K2)))/2592 + (3**(0.5D0)*MUH*
     >   (XNI(1,3,J2,K2) + XNI(3,1,J2,K2) - XNI(1,1,J2,K2) - 
     >   XNI(3,3,J2,K2)))/2592 + (3**(0.5D0)*XIH*(XNK(1,3) - 
     >   XNK(1,1) + XNK(3,1) - XNK(3,3)))/2592 
      Q2(15,28) = Q2(15,28) + 
     >   (5**(0.5D0)*ETAH*(XNJ(1,3,K2) + 28*XNJ(2,3,K2) + 
     >   XNJ(3,3,K2) - XNJ(1,1,K2) - 28*XNJ(2,1,K2) - XNJ(3,1,K2)))/
     >   32400 + (5**(0.5D0)*XIH*(XNK(1,3) - XNK(1,1) - 
     >   28*XNK(2,1) - XNK(3,1) + XNK(3,3) + 28*XNK(2,3)))/32400 - 
     >   (5**(0.5D0)*MUH*SIGN(1.0,MUH)*(5*XNI(1,3,J2,K2) + 
     >   5*XNI(3,1,J2,K2) - 5*XNI(1,1,J2,K2) - 5*XNI(3,3,J2,K2)))/
     >   32400 
      Q2(16,28) = Q2(16,28) + 
     >   (15**(0.5D0)*ETAH*SIGN(1.0,ETAH)*(55*XNJ(1,3,K2) + 
     >   340*XNJ(2,3,K2) + 55*XNJ(3,3,K2) - 55*XNJ(1,1,K2) - 
     >   340*XNJ(2,1,K2) - 55*XNJ(3,1,K2)))/486000 - (15**(0.5D0)*
     >   MUH*SIGN(1.0,MUH)*(5*XNI(1,3,J2,K2) - 5*XNI(3,1,J2,K2) - 
     >   140*XNI(2,1,J2,K2) - 5*XNI(1,1,J2,K2) + 5*XNI(3,3,J2,K2) + 
     >   140*XNI(2,3,J2,K2)))/486000 + (15**(0.5D0)*XIH*
     >   (11*XNK(1,3) + 308*XNK(1,2) + 11*XNK(1,1) + 68*XNK(2,1) + 
     >   11*XNK(3,1) + 308*XNK(3,2) + 11*XNK(3,3) + 68*XNK(2,3) + 
     >   1904*XNK(2,2)))/486000 
      Q2(17,28) = Q2(17,28) + 
     >   (5**(0.5D0)*MUH*(XNI(1,3,J2,K2) - XNI(3,1,J2,K2) - 
     >   28*XNI(2,1,J2,K2) - XNI(1,1,J2,K2) + XNI(3,3,J2,K2) + 
     >   28*XNI(2,3,J2,K2)))/32400 - (5**(0.5D0)*XIH*(XNK(1,3) + 
     >   28*XNK(1,2) + XNK(1,1) - XNK(3,1) - 28*XNK(3,2) - 
     >   XNK(3,3)))/32400 - (5**(0.5D0)*ETAH*SIGN(1.0,ETAH)*
     >   (5*XNJ(1,3,K2) - 5*XNJ(3,3,K2) - 5*XNJ(1,1,K2) + 
     >   5*XNJ(3,1,K2)))/32400 
      Q2(18,28) = Q2(18,28) 
     >   - (3**(0.5D0)*XIH*(XNK(1,3) + 28*XNK(1,2) + XNK(1,1) + 
     >   28*XNK(2,1) + XNK(3,1) + 28*XNK(3,2) + XNK(3,3) + 
     >   28*XNK(2,3) + 784*XNK(2,2)))/243000 - (3**(0.5D0)*
     >   ETAH*SIGN(1.0,ETAH)*(5*XNJ(1,3,K2) + 140*XNJ(2,3,K2) + 
     >   5*XNJ(3,3,K2) - 5*XNJ(1,1,K2) - 140*XNJ(2,1,K2) - 
     >   5*XNJ(3,1,K2)))/243000 - (3**(0.5D0)*MUH*
     >   SIGN(1.0,MUH)*(5*XNI(1,3,J2,K2) - 5*XNI(3,1,J2,K2) - 
     >   140*XNI(2,1,J2,K2) - 5*XNI(1,1,J2,K2) + 5*XNI(3,3,J2,K2) + 
     >   140*XNI(2,3,J2,K2)))/243000 
      Q2(19,28) = Q2(19,28) + 
     >   (5**(0.5D0)*XIH*SIGN(1.0,XIH)*(121*XNK(1,3) + 
     >   748*XNK(1,2) + 121*XNK(1,1) + 748*XNK(2,1) + 121*XNK(3,1) + 
     >   748*XNK(3,2) + 121*XNK(3,3) + 748*XNK(2,3) + 4624*
     >   XNK(2,2)))/486000 - (5**(0.5D0)*MUH*SIGN(1.0,MUH)*
     >   (11*XNI(1,3,J2,K2) + 11*XNI(3,1,J2,K2) + 68*XNI(2,1,J2,K2) + 
     >   11*XNI(1,1,J2,K2) + 308*XNI(3,2,J2,K2) + 1904*XNI(2,2,J2,K2) + 
     >   308*XNI(1,2,J2,K2) + 11*XNI(3,3,J2,K2) + 68*XNI(2,3,J2,K2)))/
     >   486000 - (5**(0.5D0)*ETAH*SIGN(1.0,ETAH)*(11*
     >   XNJ(1,3,K2) + 68*XNJ(2,3,K2) + 11*XNJ(3,3,K2) + 
     >   11*XNJ(1,1,K2) + 
     >   68*XNJ(2,1,K2) + 11*XNJ(3,1,K2) + 308*XNJ(1,2,K2) + 1904*
     >   XNJ(2,2,K2) + 308*XNJ(3,2,K2)))/486000 
      Q2(20,28) = Q2(20,28) + 
     >   (15**(0.5D0)*ETAH*SIGN(1.0,ETAH)*(5*XNJ(1,3,K2) - 
     >   5*XNJ(3,3,K2) + 5*XNJ(1,1,K2) - 5*XNJ(3,1,K2) + 140*
     >   XNJ(1,2,K2) - 140*XNJ(3,2,K2)))/486000 - (15**(0.5D0)*
     >   XIH*SIGN(1.0,XIH)*(55*XNK(1,3) + 340*XNK(1,2) + 
     >   55*XNK(1,1) - 55*XNK(3,1) - 340*XNK(3,2) - 55*XNK(3,3)))/
     >   486000 + (15**(0.5D0)*MUH*(11*XNI(1,3,J2,K2) + 
     >   11*XNI(3,1,J2,K2) + 68*XNI(2,1,J2,K2) + 11*XNI(1,1,J2,K2) + 
     >   308*XNI(3,2,J2,K2) + 1904*XNI(2,2,J2,K2) + 308*XNI(1,2,J2,K2) +
     >   11*XNI(3,3,J2,K2) + 68*XNI(2,3,J2,K2)))/486000 
      Q2(21,28) = Q2(21,28) + 
     >   ETAH*SIGN(1.0,ETAH)*(XNJ(1,3,K2)/243000 + 
     >   (7*XNJ(2,3,K2))/60750 + XNJ(3,3,K2)/243000 + XNJ(1,1,K2)/
     >   243000 + (7*XNJ(2,1,K2))/60750 + XNJ(3,1,K2)/243000 + 
     >   (7*XNJ(1,2,K2))/60750 + (98*XNJ(2,2,K2))/30375 + 
     >   (7*XNJ(3,2,K2))/60750) - MUH*SIGN(1.0,MUH)*
     >   ((11*XNI(1,3,J2,K2))/243000 + (11*XNI(3,1,J2,K2))/243000 + 
     >   (17*XNI(2,1,J2,K2))/60750 + (11*XNI(1,1,J2,K2))/243000 + 
     >   (77*XNI(3,2,J2,K2))/60750 + (238*XNI(2,2,J2,K2))/30375 + 
     >   (77*XNI(1,2,J2,K2))/60750 + (11*XNI(3,3,J2,K2))/243000 + 
     >   (17*XNI(2,3,J2,K2))/60750) - XIH*SIGN(1.0,XIH)*
     >   ((11*XNK(1,3))/243000 + (17*XNK(1,2))/60750 + 
     >   (11*XNK(1,1))/243000 + (77*XNK(2,1))/60750 + 
     >   (11*XNK(3,1))/243000 + (17*XNK(3,2))/60750 + 
     >   (11*XNK(3,3))/243000 + (77*XNK(2,3))/60750 + 
     >   (238*XNK(2,2))/30375) 
      Q2(22,28) = Q2(22,28) + 
     >   (15**(0.5D0)*MUH*SIGN(1.0,MUH)*(5*XNI(1,3,J2,K2) - 
     >   5*XNI(3,1,J2,K2) + 5*XNI(1,1,J2,K2) - 140*XNI(3,2,J2,K2) + 
     >   140*XNI(1,2,J2,K2) - 5*XNI(3,3,J2,K2)))/486000 + (15**(0.5D0)*
     >   XIH*SIGN(1.0,XIH)*(55*XNK(1,3) - 55*XNK(1,1) - 
     >   340*XNK(2,1) - 55*XNK(3,1) + 55*XNK(3,3) + 340*XNK(2,3)))/
     >   486000 + (15**(0.5D0)*ETAH*(11*XNJ(1,3,K2) + 
     >   68*XNJ(2,3,K2) + 11*XNJ(3,3,K2) + 11*XNJ(1,1,K2) + 
     >   68*XNJ(2,1,K2) + 11*XNJ(3,1,K2) + 308*XNJ(1,2,K2) + 
     >   1904*XNJ(2,2,K2) + 308*XNJ(3,2,K2)))/486000 
      Q2(23,28) = Q2(23,28) 
     >   - (5**(0.5D0)*ETAH*(XNJ(1,3,K2) - XNJ(3,3,K2) + 
     >   XNJ(1,1,K2) - XNJ(3,1,K2) + 28*XNJ(1,2,K2) - 28*XNJ(3,2,K2)))/
     >   32400 - (5**(0.5D0)*MUH*(XNI(1,3,J2,K2) - XNI(3,1,J2,K2) + 
     >   XNI(1,1,J2,K2) - 28*XNI(3,2,J2,K2) + 28*XNI(1,2,J2,K2) - 
     >   XNI(3,3,J2,K2)))/32400 - (5**(0.5D0)*XIH*
     >   SIGN(1.0,XIH)*(5*XNK(1,3) - 5*XNK(1,1) + 5*XNK(3,1) - 
     >   5*XNK(3,3)))/32400 
      Q2(24,28) = Q2(24,28) + 
     >   (3**(0.5D0)*MUH*SIGN(1.0,MUH)*(5*XNI(1,3,J2,K2) - 
     >   5*XNI(3,1,J2,K2) + 5*XNI(1,1,J2,K2) - 140*XNI(3,2,J2,K2) + 
     >   140*XNI(1,2,J2,K2) - 5*XNI(3,3,J2,K2)))/243000 - (3**(0.5D0)*
     >   ETAH*(XNJ(1,3,K2) + 28*XNJ(2,3,K2) + XNJ(3,3,K2) + 
     >   XNJ(1,1,K2) + 28*XNJ(2,1,K2) + XNJ(3,1,K2) + 28*XNJ(1,2,K2) + 
     >   784*XNJ(2,2,K2) + 28*XNJ(3,2,K2)))/243000 - (3**(0.5D0)*
     >   XIH*SIGN(1.0,XIH)*(5*XNK(1,3) - 5*XNK(1,1) - 
     >   140*XNK(2,1) - 5*XNK(3,1) + 5*XNK(3,3) + 140*XNK(2,3)))/243000 
      Q2(25,28) = Q2(25,28) + 
     >   MUH*SIGN(1.0,MUH)*(XNI(1,3,J2,K2)/243000 + 
     >   XNI(3,1,J2,K2)/243000 + (7*XNI(2,1,J2,K2))/60750 + 
     >   XNI(1,1,J2,K2)/243000 + (7*XNI(3,2,J2,K2))/60750 + 
     >   (98*XNI(2,2,J2,K2))/30375 + (7*XNI(1,2,J2,K2))/60750 + 
     >   XNI(3,3,J2,K2)/243000 + (7*XNI(2,3,J2,K2))/60750) - ETAH*
     >   SIGN(1.0,ETAH)*((11*XNJ(1,3,K2))/243000 + (17*XNJ(2,3,K2))/
     >   60750 + (11*XNJ(3,3,K2))/243000 + (11*XNJ(1,1,K2))/243000 + 
     >   (17*XNJ(2,1,K2))/60750 + (11*XNJ(3,1,K2))/243000 + 
     >   (77*XNJ(1,2,K2))/60750 + (238*XNJ(2,2,K2))/30375 + 
     >   (77*XNJ(3,2,K2))/60750) - XIH*SIGN(1.0,XIH)*
     >   ((11*XNK(1,3))/243000 + (77*XNK(1,2))/60750 + 
     >   (11*XNK(1,1))/243000 + (17*XNK(2,1))/60750 + (11*XNK(3,1))/
     >   243000 + (77*XNK(3,2))/60750 + (11*XNK(3,3))/243000 + 
     >   (17*XNK(2,3))/60750 + (238*XNK(2,2))/30375) 
      Q2(26,28) = Q2(26,28) + 
     >   (3**(0.5D0)*ETAH*SIGN(1.0,ETAH)*(5*XNJ(1,3,K2) - 
     >   5*XNJ(3,3,K2) + 5*XNJ(1,1,K2) - 5*XNJ(3,1,K2) + 
     >   140*XNJ(1,2,K2) -
     >   140*XNJ(3,2,K2)))/243000 - (3**(0.5D0)*MUH*
     >   (XNI(1,3,J2,K2) + XNI(3,1,J2,K2) + 28*XNI(2,1,J2,K2) + 
     >   XNI(1,1,J2,K2) + 28*XNI(3,2,J2,K2) + 784*XNI(2,2,J2,K2) + 
     >   28*XNI(1,2,J2,K2) + XNI(3,3,J2,K2) + 28*XNI(2,3,J2,K2)))/
     >   243000 + 
     >   (3**(0.5D0)*XIH*SIGN(1.0,XIH)*(5*XNK(1,3) + 
     >   140*XNK(1,2) + 5*XNK(1,1) - 5*XNK(3,1) - 140*XNK(3,2) - 
     >   5*XNK(3,3)))/243000 
      Q2(27,28) = Q2(27,28) + 
     >   (5**(0.5D0)*ETAH*SIGN(1.0,ETAH)*(XNJ(1,3,K2) + 
     >   28*XNJ(2,3,K2) + XNJ(3,3,K2) + XNJ(1,1,K2) + 28*XNJ(2,1,K2) + 
     >   XNJ(3,1,K2) + 28*XNJ(1,2,K2) + 784*XNJ(2,2,K2) + 
     >   28*XNJ(3,2,K2)))/607500 + (5**(0.5D0)*MUH*
     >   SIGN(1.0,MUH)*(XNI(1,3,J2,K2) + XNI(3,1,J2,K2) + 
     >   28*XNI(2,1,J2,K2) + XNI(1,1,J2,K2) + 28*XNI(3,2,J2,K2) + 
     >   784*XNI(2,2,J2,K2) + 28*XNI(1,2,J2,K2) + XNI(3,3,J2,K2) + 
     >   28*XNI(2,3,J2,K2)))/607500 + (5**(0.5D0)*XIH*
     >   SIGN(1.0,XIH)*(XNK(1,3) + 28*XNK(1,2) + XNK(1,1) + 
     >   28*XNK(2,1) + XNK(3,1) + 28*XNK(3,2) + XNK(3,3) + 
     >   28*XNK(2,3) + 784*XNK(2,2)))/607500 
*** ---------------------------------------------------------------- ***
*** ---------------------------------------------------------------- ***
*** ---------------------------------------------------------------- ***
      ELSE IF(IELEM.EQ.4) THEN
*** ---------------------------------------------------------------- ***
*** ---------------------------------------------------------------- ***
*** ---------------------------------------------------------------- ***
      Q2(1,1) = (125*VT)/4096 + (25*ABS(ETAH))/2048 + 
     >   (25*ABS(MUH))/2048 + (25*ABS(XIH))/2048 
      Q2(1,2) = (25*3**(0.5D0)*MUH)/512 
      Q2(1,3) = -(5**(0.5D0)*(15*VT + 6*ABS(ETAH) - 
     >   50*ABS(MUH) + 6*ABS(XIH)))/4096 
      Q2(1,4) = (25*7**(0.5D0)*MUH)/2048 
      Q2(1,5) = (25*3**(0.5D0)*ETAH)/512 
      Q2(1,6) = 0 
      Q2(1,7) = -(3*15**(0.5D0)*ETAH)/512 
      Q2(1,8) = 0 
      Q2(1,9) = -(5**(0.5D0)*(15*VT - 50*ABS(ETAH) + 
     >   6*ABS(MUH) + 6*ABS(XIH)))/4096 
      Q2(1,10) = -(3*15**(0.5D0)*MUH)/512 
      Q2(1,11) = (9*VT)/4096 - (15*ABS(ETAH))/2048 - 
     >   (15*ABS(MUH))/2048 + (9*ABS(XIH))/10240 
      Q2(1,12) = -(3*35**(0.5D0)*MUH)/2048 
      Q2(1,13) = (25*7**(0.5D0)*ETAH)/2048 
      Q2(1,14) = 0 
      Q2(1,15) = -(3*35**(0.5D0)*ETAH)/2048 
      Q2(1,16) = 0 
      Q2(1,17) = (25*3**(0.5D0)*XIH)/512 
      Q2(1,18) = 0 
      Q2(1,19) = -(3*15**(0.5D0)*XIH)/512 
      Q2(1,20) = 0 
      Q2(1,21) = 0 
      Q2(1,22) = 0 
      Q2(1,23) = 0 
      Q2(1,24) = 0 
      Q2(1,25) = -(3*15**(0.5D0)*XIH)/512 
      Q2(1,26) = 0 
      Q2(1,27) = (9*3**(0.5D0)*XIH)/2560 
      Q2(1,28) = 0 
      Q2(1,29) = 0 
      Q2(1,30) = 0 
      Q2(1,31) = 0 
      Q2(1,32) = 0 
      Q2(1,33) = -(5**(0.5D0)*(15*VT + 6*ABS(ETAH) + 
     >   6*ABS(MUH) - 50*ABS(XIH)))/4096 
      Q2(1,34) = -(3*15**(0.5D0)*MUH)/512 
      Q2(1,35) = (9*VT)/4096 + (9*ABS(ETAH))/10240 - 
     >   (15*ABS(MUH))/2048 - (15*ABS(XIH))/2048 
      Q2(1,36) = -(3*35**(0.5D0)*MUH)/2048 
      Q2(1,37) = -(3*15**(0.5D0)*ETAH)/512 
      Q2(1,38) = 0 
      Q2(1,39) = (9*3**(0.5D0)*ETAH)/2560 
      Q2(1,40) = 0 
      Q2(1,41) = (9*VT)/4096 - (15*ABS(ETAH))/2048 + 
     >   (9*ABS(MUH))/10240 - (15*ABS(XIH))/2048 
      Q2(1,42) = (9*3**(0.5D0)*MUH)/2560 
      Q2(1,43) = (9*5**(0.5D0)*(10*ABS(ETAH) - 3*VT + 
     >   10*ABS(MUH) + 10*ABS(XIH)))/102400 
      Q2(1,44) = (9*7**(0.5D0)*MUH)/10240 
      Q2(1,45) = -(3*35**(0.5D0)*ETAH)/2048 
      Q2(1,46) = 0 
      Q2(1,47) = (9*7**(0.5D0)*ETAH)/10240 
      Q2(1,48) = 0 
      Q2(1,49) = (25*7**(0.5D0)*XIH)/2048 
      Q2(1,50) = 0 
      Q2(1,51) = -(3*35**(0.5D0)*XIH)/2048 
      Q2(1,52) = 0 
      Q2(1,53) = 0 
      Q2(1,54) = 0 
      Q2(1,55) = 0 
      Q2(1,56) = 0 
      Q2(1,57) = -(3*35**(0.5D0)*XIH)/2048 
      Q2(1,58) = 0 
      Q2(1,59) = (9*7**(0.5D0)*XIH)/10240 
      Q2(1,60) = 0 
      Q2(1,61) = 0 
      Q2(1,62) = 0 
      Q2(1,63) = 0 
      Q2(1,64) = 0 

      Q2(2,1) = -(55*3**(0.5D0)*MUH)/6144 
      Q2(2,2) = (425*VT)/12288 + (85*ABS(ETAH))/6144 + 
     >   (55*ABS(MUH))/2048 + (85*ABS(XIH))/6144 
      Q2(2,3) = (185*15**(0.5D0)*MUH)/3072 
      Q2(2,4) = -(3**(0.5D0)*7**(0.5D0)*(45*VT + 18*ABS(ETAH) - 
     >   110*ABS(MUH) + 18*ABS(XIH)))/12288 
      Q2(2,5) = 0 
      Q2(2,6) = (85*3**(0.5D0)*ETAH)/1536 
      Q2(2,7) = 0 
      Q2(2,8) = -(9*7**(0.5D0)*ETAH)/512 
      Q2(2,9) = (11*15**(0.5D0)*MUH)/10240 
      Q2(2,10) = -(5**(0.5D0)*(255*VT - 850*ABS(ETAH) + 
     >   198*ABS(MUH) + 102*ABS(XIH)))/61440 
      Q2(2,11) = -(37*3**(0.5D0)*MUH)/1024 
      Q2(2,12) = (3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(45*VT - 
     >   150*ABS(ETAH) - 110*ABS(MUH) + 
     >   18*ABS(XIH)))/102400 
      Q2(2,13) = 0 
      Q2(2,14) = (85*7**(0.5D0)*ETAH)/6144 
      Q2(2,15) = 0 
      Q2(2,16) = -(21*3**(0.5D0)*ETAH)/2048 
      Q2(2,17) = 0 
      Q2(2,18) = (85*3**(0.5D0)*XIH)/1536 
      Q2(2,19) = 0 
      Q2(2,20) = -(9*7**(0.5D0)*XIH)/512 
      Q2(2,21) = 0 
      Q2(2,22) = 0 
      Q2(2,23) = 0 
      Q2(2,24) = 0 
      Q2(2,25) = 0 
      Q2(2,26) = -(17*15**(0.5D0)*XIH)/2560 
      Q2(2,27) = 0 
      Q2(2,28) = (27*35**(0.5D0)*XIH)/12800 
      Q2(2,29) = 0 
      Q2(2,30) = 0 
      Q2(2,31) = 0 
      Q2(2,32) = 0 
      Q2(2,33) = (11*15**(0.5D0)*MUH)/10240 
      Q2(2,34) = -(5**(0.5D0)*(255*VT + 102*ABS(ETAH) + 
     >   198*ABS(MUH) - 850*ABS(XIH)))/61440 
      Q2(2,35) = -(37*3**(0.5D0)*MUH)/1024 
      Q2(2,36) = (3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(45*VT + 
     >   18*ABS(ETAH) - 110*ABS(MUH) - 
     >   150*ABS(XIH)))/102400 
      Q2(2,37) = 0 
      Q2(2,38) = -(17*15**(0.5D0)*ETAH)/2560 
      Q2(2,39) = 0 
      Q2(2,40) = (27*35**(0.5D0)*ETAH)/12800 
      Q2(2,41) = -(33*3**(0.5D0)*MUH)/51200 
      Q2(2,42) = (51*VT)/20480 - (17*ABS(ETAH))/2048 + 
     >   (99*ABS(MUH))/51200 - (17*ABS(XIH))/2048 
      Q2(2,43) = (111*15**(0.5D0)*MUH)/25600 
      Q2(2,44) = (3*3**(0.5D0)*7**(0.5D0)*(30*ABS(ETAH) - 
     >   9*VT + 22*ABS(MUH) + 30*ABS(XIH)))/102400 
      Q2(2,45) = 0 
      Q2(2,46) = -(17*35**(0.5D0)*ETAH)/10240 
      Q2(2,47) = 0 
      Q2(2,48) = (63*15**(0.5D0)*ETAH)/51200 
      Q2(2,49) = 0 
      Q2(2,50) = (85*7**(0.5D0)*XIH)/6144 
      Q2(2,51) = 0 
      Q2(2,52) = -(21*3**(0.5D0)*XIH)/2048 
      Q2(2,53) = 0 
      Q2(2,54) = 0 
      Q2(2,55) = 0 
      Q2(2,56) = 0 
      Q2(2,57) = 0 
      Q2(2,58) = -(17*35**(0.5D0)*XIH)/10240 
      Q2(2,59) = 0 
      Q2(2,60) = (63*15**(0.5D0)*XIH)/51200 
      Q2(2,61) = 0 
      Q2(2,62) = 0 
      Q2(2,63) = 0 
      Q2(2,64) = 0 

      Q2(3,1) = -(3*5**(0.5D0)*(5*VT + 2*ABS(ETAH) - 
     >   10*ABS(MUH) + 2*ABS(XIH)))/4096 
      Q2(3,2) = -(15*15**(0.5D0)*MUH)/1024 
      Q2(3,3) = (45*VT)/4096 + (9*ABS(ETAH))/2048 + 
     >   (75*ABS(MUH))/2048 + (9*ABS(XIH))/2048 
      Q2(3,4) = (15*35**(0.5D0)*MUH)/2048 
      Q2(3,5) = -(3*15**(0.5D0)*ETAH)/512 
      Q2(3,6) = 0 
      Q2(3,7) = (9*3**(0.5D0)*ETAH)/512 
      Q2(3,8) = 0 
      Q2(3,9) = (9*VT)/4096 - (15*ABS(ETAH))/2048 - 
     >   (9*ABS(MUH))/2048 + (9*ABS(XIH))/10240 
      Q2(3,10) = (9*3**(0.5D0)*MUH)/1024 
      Q2(3,11) = -(9*5**(0.5D0)*(15*VT - 50*ABS(ETAH) + 
     >   50*ABS(MUH) + 6*ABS(XIH)))/102400 
      Q2(3,12) = -(9*7**(0.5D0)*MUH)/2048 
      Q2(3,13) = -(3*35**(0.5D0)*ETAH)/2048 
      Q2(3,14) = 0 
      Q2(3,15) = (9*7**(0.5D0)*ETAH)/2048 
      Q2(3,16) = 0 
      Q2(3,17) = -(3*15**(0.5D0)*XIH)/512 
      Q2(3,18) = 0 
      Q2(3,19) = (9*3**(0.5D0)*XIH)/512 
      Q2(3,20) = 0 
      Q2(3,21) = 0 
      Q2(3,22) = 0 
      Q2(3,23) = 0 
      Q2(3,24) = 0 
      Q2(3,25) = (9*3**(0.5D0)*XIH)/2560 
      Q2(3,26) = 0 
      Q2(3,27) = -(27*15**(0.5D0)*XIH)/12800 
      Q2(3,28) = 0 
      Q2(3,29) = 0 
      Q2(3,30) = 0 
      Q2(3,31) = 0 
      Q2(3,32) = 0 
      Q2(3,33) = (9*VT)/4096 + (9*ABS(ETAH))/10240 - 
     >   (9*ABS(MUH))/2048 - (15*ABS(XIH))/2048 
      Q2(3,34) = (9*3**(0.5D0)*MUH)/1024 
      Q2(3,35) = -(9*5**(0.5D0)*(15*VT + 6*ABS(ETAH) + 
     >   50*ABS(MUH) - 50*ABS(XIH)))/102400 
      Q2(3,36) = -(9*7**(0.5D0)*MUH)/2048 
      Q2(3,37) = (9*3**(0.5D0)*ETAH)/2560 
      Q2(3,38) = 0 
      Q2(3,39) = -(27*15**(0.5D0)*ETAH)/12800 
      Q2(3,40) = 0 
      Q2(3,41) = (9*5**(0.5D0)*(10*ABS(ETAH) - 3*VT + 
     >   6*ABS(MUH) + 10*ABS(XIH)))/102400 
      Q2(3,42) = -(27*15**(0.5D0)*MUH)/25600 
      Q2(3,43) = (81*VT)/102400 - (27*ABS(ETAH))/10240 + 
     >   (27*ABS(MUH))/10240 - (27*ABS(XIH))/10240 
      Q2(3,44) = (27*35**(0.5D0)*MUH)/51200 
      Q2(3,45) = (9*7**(0.5D0)*ETAH)/10240 
      Q2(3,46) = 0 
      Q2(3,47) = -(27*35**(0.5D0)*ETAH)/51200 
      Q2(3,48) = 0 
      Q2(3,49) = -(3*35**(0.5D0)*XIH)/2048 
      Q2(3,50) = 0 
      Q2(3,51) = (9*7**(0.5D0)*XIH)/2048 
      Q2(3,52) = 0 
      Q2(3,53) = 0 
      Q2(3,54) = 0 
      Q2(3,55) = 0 
      Q2(3,56) = 0 
      Q2(3,57) = (9*7**(0.5D0)*XIH)/10240 
      Q2(3,58) = 0 
      Q2(3,59) = -(27*35**(0.5D0)*XIH)/51200 
      Q2(3,60) = 0 
      Q2(3,61) = 0 
      Q2(3,62) = 0 
      Q2(3,63) = 0 
      Q2(3,64) = 0 

      Q2(4,1) = -(45*7**(0.5D0)*MUH)/14336 
      Q2(4,2) = -(3*3**(0.5D0)*7**(0.5D0)*(35*VT + 
     >   14*ABS(ETAH) - 30*ABS(MUH) + 
     >   14*ABS(XIH)))/28672 
      Q2(4,3) = -(45*35**(0.5D0)*MUH)/1792 
      Q2(4,4) = (405*VT)/28672 + (81*ABS(ETAH))/14336 + 
     >   (45*ABS(MUH))/2048 + (81*ABS(XIH))/14336 
      Q2(4,5) = 0 
      Q2(4,6) = -(9*7**(0.5D0)*ETAH)/512 
      Q2(4,7) = 0 
      Q2(4,8) = (81*3**(0.5D0)*ETAH)/3584 
      Q2(4,9) = (27*35**(0.5D0)*MUH)/71680 
      Q2(4,10) = (3*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(105*VT - 
     >   350*ABS(ETAH) - 90*ABS(MUH) + 
     >   42*ABS(XIH)))/716800 
      Q2(4,11) = (27*7**(0.5D0)*MUH)/1792 
      Q2(4,12) = -(27*5**(0.5D0)*(45*VT - 150*ABS(ETAH) + 
     >   70*ABS(MUH) + 18*ABS(XIH)))/716800 
      Q2(4,13) = 0 
      Q2(4,14) = -(21*3**(0.5D0)*ETAH)/2048 
      Q2(4,15) = 0 
      Q2(4,16) = (81*7**(0.5D0)*ETAH)/14336 
      Q2(4,17) = 0 
      Q2(4,18) = -(9*7**(0.5D0)*XIH)/512 
      Q2(4,19) = 0 
      Q2(4,20) = (81*3**(0.5D0)*XIH)/3584 
      Q2(4,21) = 0 
      Q2(4,22) = 0 
      Q2(4,23) = 0 
      Q2(4,24) = 0 
      Q2(4,25) = 0 
      Q2(4,26) = (27*35**(0.5D0)*XIH)/12800 
      Q2(4,27) = 0 
      Q2(4,28) = -(243*15**(0.5D0)*XIH)/89600 
      Q2(4,29) = 0 
      Q2(4,30) = 0 
      Q2(4,31) = 0 
      Q2(4,32) = 0 
      Q2(4,33) = (27*35**(0.5D0)*MUH)/71680 
      Q2(4,34) = (3*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(105*VT + 
     >   42*ABS(ETAH) - 90*ABS(MUH) - 
     >   350*ABS(XIH)))/716800 
      Q2(4,35) = (27*7**(0.5D0)*MUH)/1792 
      Q2(4,36) = -(27*5**(0.5D0)*(45*VT + 18*ABS(ETAH) + 
     >   70*ABS(MUH) - 150*ABS(XIH)))/716800 
      Q2(4,37) = 0 
      Q2(4,38) = (27*35**(0.5D0)*ETAH)/12800 
      Q2(4,39) = 0 
      Q2(4,40) = -(243*15**(0.5D0)*ETAH)/89600 
      Q2(4,41) = -(81*7**(0.5D0)*MUH)/358400 
      Q2(4,42) = (9*3**(0.5D0)*7**(0.5D0)*(70*ABS(ETAH) - 
     >   21*VT + 18*ABS(MUH) + 70*ABS(XIH)))/716800 
      Q2(4,43) = -(81*35**(0.5D0)*MUH)/44800 
      Q2(4,44) = (729*VT)/716800 - (243*ABS(ETAH))/71680 + 
     >   (81*ABS(MUH))/51200 - (243*ABS(XIH))/71680 
      Q2(4,45) = 0 
      Q2(4,46) = (63*15**(0.5D0)*ETAH)/51200 
      Q2(4,47) = 0 
      Q2(4,48) = -(243*35**(0.5D0)*ETAH)/358400 
      Q2(4,49) = 0 
      Q2(4,50) = -(21*3**(0.5D0)*XIH)/2048 
      Q2(4,51) = 0 
      Q2(4,52) = (81*7**(0.5D0)*XIH)/14336 
      Q2(4,53) = 0 
      Q2(4,54) = 0 
      Q2(4,55) = 0 
      Q2(4,56) = 0 
      Q2(4,57) = 0 
      Q2(4,58) = (63*15**(0.5D0)*XIH)/51200 
      Q2(4,59) = 0 
      Q2(4,60) = -(243*35**(0.5D0)*XIH)/358400 
      Q2(4,61) = 0 
      Q2(4,62) = 0 
      Q2(4,63) = 0 
      Q2(4,64) = 0 

      Q2(5,1) = -(55*3**(0.5D0)*ETAH)/6144 
      Q2(5,2) = 0 
      Q2(5,3) = (11*15**(0.5D0)*ETAH)/10240 
      Q2(5,4) = 0 
      Q2(5,5) = (425*VT)/12288 + (55*ABS(ETAH))/2048 + 
     >   (85*ABS(MUH))/6144 + (85*ABS(XIH))/6144 
      Q2(5,6) = (85*3**(0.5D0)*MUH)/1536 
      Q2(5,7) = -(5**(0.5D0)*(255*VT + 198*ABS(ETAH) - 
     >   850*ABS(MUH) + 102*ABS(XIH)))/61440 
      Q2(5,8) = (85*7**(0.5D0)*MUH)/6144 
      Q2(5,9) = (185*15**(0.5D0)*ETAH)/3072 
      Q2(5,10) = 0 
      Q2(5,11) = -(37*3**(0.5D0)*ETAH)/1024 
      Q2(5,12) = 0 
      Q2(5,13) = -(3**(0.5D0)*7**(0.5D0)*(45*VT - 
     >   110*ABS(ETAH) + 18*ABS(MUH) + 
     >   18*ABS(XIH)))/12288 
      Q2(5,14) = -(9*7**(0.5D0)*MUH)/512 
      Q2(5,15) = (3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(45*VT - 
     >   110*ABS(ETAH) - 150*ABS(MUH) + 
     >   18*ABS(XIH)))/102400 
      Q2(5,16) = -(21*3**(0.5D0)*MUH)/2048 
      Q2(5,17) = 0 
      Q2(5,18) = 0 
      Q2(5,19) = 0 
      Q2(5,20) = 0 
      Q2(5,21) = (85*3**(0.5D0)*XIH)/1536 
      Q2(5,22) = 0 
      Q2(5,23) = -(17*15**(0.5D0)*XIH)/2560 
      Q2(5,24) = 0 
      Q2(5,25) = 0 
      Q2(5,26) = 0 
      Q2(5,27) = 0 
      Q2(5,28) = 0 
      Q2(5,29) = -(9*7**(0.5D0)*XIH)/512 
      Q2(5,30) = 0 
      Q2(5,31) = (27*35**(0.5D0)*XIH)/12800 
      Q2(5,32) = 0 
      Q2(5,33) = (11*15**(0.5D0)*ETAH)/10240 
      Q2(5,34) = 0 
      Q2(5,35) = -(33*3**(0.5D0)*ETAH)/51200 
      Q2(5,36) = 0 
      Q2(5,37) = -(5**(0.5D0)*(255*VT + 198*ABS(ETAH) + 
     >   102*ABS(MUH) - 850*ABS(XIH)))/61440 
      Q2(5,38) = -(17*15**(0.5D0)*MUH)/2560 
      Q2(5,39) = (51*VT)/20480 + (99*ABS(ETAH))/51200 - 
     >   (17*ABS(MUH))/2048 - (17*ABS(XIH))/2048 
      Q2(5,40) = -(17*35**(0.5D0)*MUH)/10240 
      Q2(5,41) = -(37*3**(0.5D0)*ETAH)/1024 
      Q2(5,42) = 0 
      Q2(5,43) = (111*15**(0.5D0)*ETAH)/25600 
      Q2(5,44) = 0 
      Q2(5,45) = (3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(45*VT - 
     >   110*ABS(ETAH) + 18*ABS(MUH) - 
     >   150*ABS(XIH)))/102400 
      Q2(5,46) = (27*35**(0.5D0)*MUH)/12800 
      Q2(5,47) = (3*3**(0.5D0)*7**(0.5D0)*(22*ABS(ETAH) - 
     >   9*VT + 30*ABS(MUH) + 30*ABS(XIH)))/102400 
      Q2(5,48) = (63*15**(0.5D0)*MUH)/51200 
      Q2(5,49) = 0 
      Q2(5,50) = 0 
      Q2(5,51) = 0 
      Q2(5,52) = 0 
      Q2(5,53) = (85*7**(0.5D0)*XIH)/6144 
      Q2(5,54) = 0 
      Q2(5,55) = -(17*35**(0.5D0)*XIH)/10240 
      Q2(5,56) = 0 
      Q2(5,57) = 0 
      Q2(5,58) = 0 
      Q2(5,59) = 0 
      Q2(5,60) = 0 
      Q2(5,61) = -(21*3**(0.5D0)*XIH)/2048 
      Q2(5,62) = 0 
      Q2(5,63) = (63*15**(0.5D0)*XIH)/51200 
      Q2(5,64) = 0 

      Q2(6,1) = 0 
      Q2(6,2) = -(187*3**(0.5D0)*ETAH)/18432 
      Q2(6,3) = 0 
      Q2(6,4) = (33*7**(0.5D0)*ETAH)/10240 
      Q2(6,5) = -(187*3**(0.5D0)*MUH)/18432 
      Q2(6,6) = (1445*VT)/36864 + (187*ABS(ETAH))/6144 + 
     >   (187*ABS(MUH))/6144 + (289*ABS(XIH))/18432 
      Q2(6,7) = (629*15**(0.5D0)*MUH)/9216 
      Q2(6,8) = -(3**(0.5D0)*7**(0.5D0)*(765*VT + 
     >   594*ABS(ETAH) - 1870*ABS(MUH) + 
     >   306*ABS(XIH)))/184320 
      Q2(6,9) = 0 
      Q2(6,10) = (629*15**(0.5D0)*ETAH)/9216 
      Q2(6,11) = 0 
      Q2(6,12) = -(111*35**(0.5D0)*ETAH)/5120 
      Q2(6,13) = (33*7**(0.5D0)*MUH)/10240 
      Q2(6,14) = -(3**(0.5D0)*7**(0.5D0)*(765*VT - 
     >   1870*ABS(ETAH) + 594*ABS(MUH) + 
     >   306*ABS(XIH)))/184320 
      Q2(6,15) = -(111*35**(0.5D0)*MUH)/5120 
      Q2(6,16) = (189*VT)/20480 - (231*ABS(ETAH))/10240 - 
     >   (231*ABS(MUH))/10240 + (189*ABS(XIH))/51200 
      Q2(6,17) = 0 
      Q2(6,18) = 0 
      Q2(6,19) = 0 
      Q2(6,20) = 0 
      Q2(6,21) = 0 
      Q2(6,22) = (289*3**(0.5D0)*XIH)/4608 
      Q2(6,23) = 0 
      Q2(6,24) = -(51*7**(0.5D0)*XIH)/2560 
      Q2(6,25) = 0 
      Q2(6,26) = 0 
      Q2(6,27) = 0 
      Q2(6,28) = 0 
      Q2(6,29) = 0 
      Q2(6,30) = -(51*7**(0.5D0)*XIH)/2560 
      Q2(6,31) = 0 
      Q2(6,32) = (189*3**(0.5D0)*XIH)/12800 
      Q2(6,33) = 0 
      Q2(6,34) = (187*15**(0.5D0)*ETAH)/153600 
      Q2(6,35) = 0 
      Q2(6,36) = -(99*35**(0.5D0)*ETAH)/256000 
      Q2(6,37) = (187*15**(0.5D0)*MUH)/153600 
      Q2(6,38) = -(17*5**(0.5D0)*(255*VT + 198*ABS(ETAH) + 
     >   198*ABS(MUH) - 850*ABS(XIH)))/921600 
      Q2(6,39) = -(629*3**(0.5D0)*MUH)/15360 
      Q2(6,40) = (3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(765*VT + 
     >   594*ABS(ETAH) - 1870*ABS(MUH) - 
     >   2550*ABS(XIH)))/1536000 
      Q2(6,41) = 0 
      Q2(6,42) = -(629*3**(0.5D0)*ETAH)/15360 
      Q2(6,43) = 0 
      Q2(6,44) = (333*7**(0.5D0)*ETAH)/25600 
      Q2(6,45) = -(99*35**(0.5D0)*MUH)/256000 
      Q2(6,46) = (3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(765*VT - 
     >   1870*ABS(ETAH) + 594*ABS(MUH) - 
     >   2550*ABS(XIH)))/1536000 
      Q2(6,47) = (333*7**(0.5D0)*MUH)/25600 
      Q2(6,48) = (63*5**(0.5D0)*(22*ABS(ETAH) - 9*VT + 
     >   22*ABS(MUH) + 30*ABS(XIH)))/512000 
      Q2(6,49) = 0 
      Q2(6,50) = 0 
      Q2(6,51) = 0 
      Q2(6,52) = 0 
      Q2(6,53) = 0 
      Q2(6,54) = (289*7**(0.5D0)*XIH)/18432 
      Q2(6,55) = 0 
      Q2(6,56) = -(119*3**(0.5D0)*XIH)/10240 
      Q2(6,57) = 0 
      Q2(6,58) = 0 
      Q2(6,59) = 0 
      Q2(6,60) = 0 
      Q2(6,61) = 0 
      Q2(6,62) = -(119*3**(0.5D0)*XIH)/10240 
      Q2(6,63) = 0 
      Q2(6,64) = (189*7**(0.5D0)*XIH)/51200 

      Q2(7,1) = (11*15**(0.5D0)*ETAH)/10240 
      Q2(7,2) = 0 
      Q2(7,3) = -(33*3**(0.5D0)*ETAH)/10240 
      Q2(7,4) = 0 
      Q2(7,5) = -(5**(0.5D0)*(85*VT + 66*ABS(ETAH) - 
     >   170*ABS(MUH) + 34*ABS(XIH)))/20480 
      Q2(7,6) = -(17*15**(0.5D0)*MUH)/1024 
      Q2(7,7) = (51*VT)/4096 + (99*ABS(ETAH))/10240 + 
     >   (85*ABS(MUH))/2048 + (51*ABS(XIH))/10240 
      Q2(7,8) = (17*35**(0.5D0)*MUH)/2048 
      Q2(7,9) = -(37*3**(0.5D0)*ETAH)/1024 
      Q2(7,10) = 0 
      Q2(7,11) = (111*15**(0.5D0)*ETAH)/5120 
      Q2(7,12) = 0 
      Q2(7,13) = (3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(45*VT - 
     >   110*ABS(ETAH) - 90*ABS(MUH) + 
     >   18*ABS(XIH)))/102400 
      Q2(7,14) = (27*35**(0.5D0)*MUH)/5120 
      Q2(7,15) = -(3*3**(0.5D0)*7**(0.5D0)*(45*VT - 
     >   110*ABS(ETAH) + 150*ABS(MUH) + 
     >   18*ABS(XIH)))/102400 
      Q2(7,16) = -(63*15**(0.5D0)*MUH)/10240 
      Q2(7,17) = 0 
      Q2(7,18) = 0 
      Q2(7,19) = 0 
      Q2(7,20) = 0 
      Q2(7,21) = -(17*15**(0.5D0)*XIH)/2560 
      Q2(7,22) = 0 
      Q2(7,23) = (51*3**(0.5D0)*XIH)/2560 
      Q2(7,24) = 0 
      Q2(7,25) = 0 
      Q2(7,26) = 0 
      Q2(7,27) = 0 
      Q2(7,28) = 0 
      Q2(7,29) = (27*35**(0.5D0)*XIH)/12800 
      Q2(7,30) = 0 
      Q2(7,31) = -(81*7**(0.5D0)*XIH)/12800 
      Q2(7,32) = 0 
      Q2(7,33) = -(33*3**(0.5D0)*ETAH)/51200 
      Q2(7,34) = 0 
      Q2(7,35) = (99*15**(0.5D0)*ETAH)/256000 
      Q2(7,36) = 0 
      Q2(7,37) = (51*VT)/20480 + (99*ABS(ETAH))/51200 - 
     >   (51*ABS(MUH))/10240 - (17*ABS(XIH))/2048 
      Q2(7,38) = (51*3**(0.5D0)*MUH)/5120 
      Q2(7,39) = -(3*5**(0.5D0)*(255*VT + 198*ABS(ETAH) + 
     >   850*ABS(MUH) - 850*ABS(XIH)))/512000 
      Q2(7,40) = -(51*7**(0.5D0)*MUH)/10240 
      Q2(7,41) = (111*15**(0.5D0)*ETAH)/25600 
      Q2(7,42) = 0 
      Q2(7,43) = -(333*3**(0.5D0)*ETAH)/25600 
      Q2(7,44) = 0 
      Q2(7,45) = (3*3**(0.5D0)*7**(0.5D0)*(22*ABS(ETAH) - 
     >   9*VT + 18*ABS(MUH) + 30*ABS(XIH)))/102400 
      Q2(7,46) = -(81*7**(0.5D0)*MUH)/25600 
      Q2(7,47) = (9*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(9*VT - 
     >   22*ABS(ETAH) + 30*ABS(MUH) - 
     >   30*ABS(XIH)))/512000 
      Q2(7,48) = (189*3**(0.5D0)*MUH)/51200 
      Q2(7,49) = 0 
      Q2(7,50) = 0 
      Q2(7,51) = 0 
      Q2(7,52) = 0 
      Q2(7,53) = -(17*35**(0.5D0)*XIH)/10240 
      Q2(7,54) = 0 
      Q2(7,55) = (51*7**(0.5D0)*XIH)/10240 
      Q2(7,56) = 0 
      Q2(7,57) = 0 
      Q2(7,58) = 0 
      Q2(7,59) = 0 
      Q2(7,60) = 0 
      Q2(7,61) = (63*15**(0.5D0)*XIH)/51200 
      Q2(7,62) = 0 
      Q2(7,63) = -(189*3**(0.5D0)*XIH)/51200 
      Q2(7,64) = 0 

      Q2(8,1) = 0 
      Q2(8,2) = (33*7**(0.5D0)*ETAH)/10240 
      Q2(8,3) = 0 
      Q2(8,4) = -(297*3**(0.5D0)*ETAH)/71680 
      Q2(8,5) = -(51*7**(0.5D0)*MUH)/14336 
      Q2(8,6) = -(3**(0.5D0)*7**(0.5D0)*(595*VT + 
     >   462*ABS(ETAH) - 510*ABS(MUH) + 
     >   238*ABS(XIH)))/143360 
      Q2(8,7) = -(51*35**(0.5D0)*MUH)/1792 
      Q2(8,8) = (459*VT)/28672 + (891*ABS(ETAH))/71680 + 
     >   (51*ABS(MUH))/2048 + (459*ABS(XIH))/71680 
      Q2(8,9) = 0 
      Q2(8,10) = -(111*35**(0.5D0)*ETAH)/5120 
      Q2(8,11) = 0 
      Q2(8,12) = (999*15**(0.5D0)*ETAH)/35840 
      Q2(8,13) = (27*3**(0.5D0)*MUH)/10240 
      Q2(8,14) = (189*VT)/20480 - (231*ABS(ETAH))/10240 - 
     >   (81*ABS(MUH))/10240 + (189*ABS(XIH))/51200 
      Q2(8,15) = (27*15**(0.5D0)*MUH)/1280 
      Q2(8,16) = -(27*3**(0.5D0)*7**(0.5D0)*(45*VT - 
     >   110*ABS(ETAH) + 70*ABS(MUH) + 
     >   18*ABS(XIH)))/716800 
      Q2(8,17) = 0 
      Q2(8,18) = 0 
      Q2(8,19) = 0 
      Q2(8,20) = 0 
      Q2(8,21) = 0 
      Q2(8,22) = -(51*7**(0.5D0)*XIH)/2560 
      Q2(8,23) = 0 
      Q2(8,24) = (459*3**(0.5D0)*XIH)/17920 
      Q2(8,25) = 0 
      Q2(8,26) = 0 
      Q2(8,27) = 0 
      Q2(8,28) = 0 
      Q2(8,29) = 0 
      Q2(8,30) = (189*3**(0.5D0)*XIH)/12800 
      Q2(8,31) = 0 
      Q2(8,32) = -(729*7**(0.5D0)*XIH)/89600 
      Q2(8,33) = 0 
      Q2(8,34) = -(99*35**(0.5D0)*ETAH)/256000 
      Q2(8,35) = 0 
      Q2(8,36) = (891*15**(0.5D0)*ETAH)/1792000 
      Q2(8,37) = (153*35**(0.5D0)*MUH)/358400 
      Q2(8,38) = (3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(1785*VT + 
     >   1386*ABS(ETAH) - 1530*ABS(MUH) - 
     >   5950*ABS(XIH)))/3584000 
      Q2(8,39) = (153*7**(0.5D0)*MUH)/8960 
      Q2(8,40) = -(9*5**(0.5D0)*(765*VT + 594*ABS(ETAH) + 
     >   1190*ABS(MUH) - 2550*ABS(XIH)))/3584000 
      Q2(8,41) = 0 
      Q2(8,42) = (333*7**(0.5D0)*ETAH)/25600 
      Q2(8,43) = 0 
      Q2(8,44) = -(2997*3**(0.5D0)*ETAH)/179200 
      Q2(8,45) = -(81*15**(0.5D0)*MUH)/256000 
      Q2(8,46) = (9*5**(0.5D0)*(154*ABS(ETAH) - 63*VT + 
     >   54*ABS(MUH) + 210*ABS(XIH)))/512000 
      Q2(8,47) = -(81*3**(0.5D0)*MUH)/6400 
      Q2(8,48) = (81*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(9*VT - 
     >   22*ABS(ETAH) + 14*ABS(MUH) - 
     >   30*ABS(XIH)))/3584000 
      Q2(8,49) = 0 
      Q2(8,50) = 0 
      Q2(8,51) = 0 
      Q2(8,52) = 0 
      Q2(8,53) = 0 
      Q2(8,54) = -(119*3**(0.5D0)*XIH)/10240 
      Q2(8,55) = 0 
      Q2(8,56) = (459*7**(0.5D0)*XIH)/71680 
      Q2(8,57) = 0 
      Q2(8,58) = 0 
      Q2(8,59) = 0 
      Q2(8,60) = 0 
      Q2(8,61) = 0 
      Q2(8,62) = (189*7**(0.5D0)*XIH)/51200 
      Q2(8,63) = 0 
      Q2(8,64) = -(243*3**(0.5D0)*XIH)/51200 

      Q2(9,1) = -(3*5**(0.5D0)*(5*VT - 10*ABS(ETAH) + 
     >   2*ABS(MUH) + 2*ABS(XIH)))/4096 
      Q2(9,2) = -(3*15**(0.5D0)*MUH)/512 
      Q2(9,3) = (9*VT)/4096 - (9*ABS(ETAH))/2048 - 
     >   (15*ABS(MUH))/2048 + (9*ABS(XIH))/10240 
      Q2(9,4) = -(3*35**(0.5D0)*MUH)/2048 
      Q2(9,5) = -(15*15**(0.5D0)*ETAH)/1024 
      Q2(9,6) = 0 
      Q2(9,7) = (9*3**(0.5D0)*ETAH)/1024 
      Q2(9,8) = 0 
      Q2(9,9) = (45*VT)/4096 + (75*ABS(ETAH))/2048 + 
     >   (9*ABS(MUH))/2048 + (9*ABS(XIH))/2048 
      Q2(9,10) = (9*3**(0.5D0)*MUH)/512 
      Q2(9,11) = -(9*5**(0.5D0)*(15*VT + 50*ABS(ETAH) - 
     >   50*ABS(MUH) + 6*ABS(XIH)))/102400 
      Q2(9,12) = (9*7**(0.5D0)*MUH)/2048 
      Q2(9,13) = (15*35**(0.5D0)*ETAH)/2048 
      Q2(9,14) = 0 
      Q2(9,15) = -(9*7**(0.5D0)*ETAH)/2048 
      Q2(9,16) = 0 
      Q2(9,17) = -(3*15**(0.5D0)*XIH)/512 
      Q2(9,18) = 0 
      Q2(9,19) = (9*3**(0.5D0)*XIH)/2560 
      Q2(9,20) = 0 
      Q2(9,21) = 0 
      Q2(9,22) = 0 
      Q2(9,23) = 0 
      Q2(9,24) = 0 
      Q2(9,25) = (9*3**(0.5D0)*XIH)/512 
      Q2(9,26) = 0 
      Q2(9,27) = -(27*15**(0.5D0)*XIH)/12800 
      Q2(9,28) = 0 
      Q2(9,29) = 0 
      Q2(9,30) = 0 
      Q2(9,31) = 0 
      Q2(9,32) = 0 
      Q2(9,33) = (9*VT)/4096 - (9*ABS(ETAH))/2048 + 
     >   (9*ABS(MUH))/10240 - (15*ABS(XIH))/2048 
      Q2(9,34) = (9*3**(0.5D0)*MUH)/2560 
      Q2(9,35) = (9*5**(0.5D0)*(6*ABS(ETAH) - 3*VT + 
     >   10*ABS(MUH) + 10*ABS(XIH)))/102400 
      Q2(9,36) = (9*7**(0.5D0)*MUH)/10240 
      Q2(9,37) = (9*3**(0.5D0)*ETAH)/1024 
      Q2(9,38) = 0 
      Q2(9,39) = -(27*15**(0.5D0)*ETAH)/25600 
      Q2(9,40) = 0 
      Q2(9,41) = -(9*5**(0.5D0)*(15*VT + 50*ABS(ETAH) + 
     >   6*ABS(MUH) - 50*ABS(XIH)))/102400 
      Q2(9,42) = -(27*15**(0.5D0)*MUH)/12800 
      Q2(9,43) = (81*VT)/102400 + (27*ABS(ETAH))/10240 - 
     >   (27*ABS(MUH))/10240 - (27*ABS(XIH))/10240 
      Q2(9,44) = -(27*35**(0.5D0)*MUH)/51200 
      Q2(9,45) = -(9*7**(0.5D0)*ETAH)/2048 
      Q2(9,46) = 0 
      Q2(9,47) = (27*35**(0.5D0)*ETAH)/51200 
      Q2(9,48) = 0 
      Q2(9,49) = -(3*35**(0.5D0)*XIH)/2048 
      Q2(9,50) = 0 
      Q2(9,51) = (9*7**(0.5D0)*XIH)/10240 
      Q2(9,52) = 0 
      Q2(9,53) = 0 
      Q2(9,54) = 0 
      Q2(9,55) = 0 
      Q2(9,56) = 0 
      Q2(9,57) = (9*7**(0.5D0)*XIH)/2048 
      Q2(9,58) = 0 
      Q2(9,59) = -(27*35**(0.5D0)*XIH)/51200 
      Q2(9,60) = 0 
      Q2(9,61) = 0 
      Q2(9,62) = 0 
      Q2(9,63) = 0 
      Q2(9,64) = 0 

      Q2(10,1) = (11*15**(0.5D0)*MUH)/10240 
      Q2(10,2) = -(5**(0.5D0)*(85*VT - 170*ABS(ETAH) + 
     >   66*ABS(MUH) + 34*ABS(XIH)))/20480 
      Q2(10,3) = -(37*3**(0.5D0)*MUH)/1024 
      Q2(10,4) = (3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(45*VT - 
     >   90*ABS(ETAH) - 110*ABS(MUH) + 
     >   18*ABS(XIH)))/102400 
      Q2(10,5) = 0 
      Q2(10,6) = -(17*15**(0.5D0)*ETAH)/1024 
      Q2(10,7) = 0 
      Q2(10,8) = (27*35**(0.5D0)*ETAH)/5120 
      Q2(10,9) = -(33*3**(0.5D0)*MUH)/10240 
      Q2(10,10) = (51*VT)/4096 + (85*ABS(ETAH))/2048 + 
     >   (99*ABS(MUH))/10240 + (51*ABS(XIH))/10240 
      Q2(10,11) = (111*15**(0.5D0)*MUH)/5120 
      Q2(10,12) = -(3*3**(0.5D0)*7**(0.5D0)*(45*VT + 
     >   150*ABS(ETAH) - 110*ABS(MUH) + 
     >   18*ABS(XIH)))/102400 
      Q2(10,13) = 0 
      Q2(10,14) = (17*35**(0.5D0)*ETAH)/2048 
      Q2(10,15) = 0 
      Q2(10,16) = -(63*15**(0.5D0)*ETAH)/10240 
      Q2(10,17) = 0 
      Q2(10,18) = -(17*15**(0.5D0)*XIH)/2560 
      Q2(10,19) = 0 
      Q2(10,20) = (27*35**(0.5D0)*XIH)/12800 
      Q2(10,21) = 0 
      Q2(10,22) = 0 
      Q2(10,23) = 0 
      Q2(10,24) = 0 
      Q2(10,25) = 0 
      Q2(10,26) = (51*3**(0.5D0)*XIH)/2560 
      Q2(10,27) = 0 
      Q2(10,28) = -(81*7**(0.5D0)*XIH)/12800 
      Q2(10,29) = 0 
      Q2(10,30) = 0 
      Q2(10,31) = 0 
      Q2(10,32) = 0 
      Q2(10,33) = -(33*3**(0.5D0)*MUH)/51200 
      Q2(10,34) = (51*VT)/20480 - (51*ABS(ETAH))/10240 + 
     >   (99*ABS(MUH))/51200 - (17*ABS(XIH))/2048 
      Q2(10,35) = (111*15**(0.5D0)*MUH)/25600 
      Q2(10,36) = (3*3**(0.5D0)*7**(0.5D0)*(18*ABS(ETAH) - 
     >   9*VT + 22*ABS(MUH) + 30*ABS(XIH)))/102400 
      Q2(10,37) = 0 
      Q2(10,38) = (51*3**(0.5D0)*ETAH)/5120 
      Q2(10,39) = 0 
      Q2(10,40) = -(81*7**(0.5D0)*ETAH)/25600 
      Q2(10,41) = (99*15**(0.5D0)*MUH)/256000 
      Q2(10,42) = -(3*5**(0.5D0)*(255*VT + 850*ABS(ETAH) + 
     >   198*ABS(MUH) - 850*ABS(XIH)))/512000 
      Q2(10,43) = -(333*3**(0.5D0)*MUH)/25600 
      Q2(10,44) = (9*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(9*VT + 
     >   30*ABS(ETAH) - 22*ABS(MUH) - 
     >   30*ABS(XIH)))/512000 
      Q2(10,45) = 0 
      Q2(10,46) = -(51*7**(0.5D0)*ETAH)/10240 
      Q2(10,47) = 0 
      Q2(10,48) = (189*3**(0.5D0)*ETAH)/51200 
      Q2(10,49) = 0 
      Q2(10,50) = -(17*35**(0.5D0)*XIH)/10240 
      Q2(10,51) = 0 
      Q2(10,52) = (63*15**(0.5D0)*XIH)/51200 
      Q2(10,53) = 0 
      Q2(10,54) = 0 
      Q2(10,55) = 0 
      Q2(10,56) = 0 
      Q2(10,57) = 0 
      Q2(10,58) = (51*7**(0.5D0)*XIH)/10240 
      Q2(10,59) = 0 
      Q2(10,60) = -(189*3**(0.5D0)*XIH)/51200 
      Q2(10,61) = 0 
      Q2(10,62) = 0 
      Q2(10,63) = 0 
      Q2(10,64) = 0 

      Q2(11,1) = (9*VT)/4096 - (9*ABS(ETAH))/2048 - 
     >   (9*ABS(MUH))/2048 + (9*ABS(XIH))/10240 
      Q2(11,2) = (9*3**(0.5D0)*MUH)/1024 
      Q2(11,3) = -(9*5**(0.5D0)*(15*VT - 30*ABS(ETAH) + 
     >   50*ABS(MUH) + 6*ABS(XIH)))/102400 
      Q2(11,4) = -(9*7**(0.5D0)*MUH)/2048 
      Q2(11,5) = (9*3**(0.5D0)*ETAH)/1024 
      Q2(11,6) = 0 
      Q2(11,7) = -(27*15**(0.5D0)*ETAH)/5120 
      Q2(11,8) = 0 
      Q2(11,9) = -(9*5**(0.5D0)*(15*VT + 50*ABS(ETAH) - 
     >   30*ABS(MUH) + 6*ABS(XIH)))/102400 
      Q2(11,10) = -(27*15**(0.5D0)*MUH)/5120 
      Q2(11,11) = (81*VT)/20480 + (27*ABS(ETAH))/2048 + 
     >   (27*ABS(MUH))/2048 + (81*ABS(XIH))/51200 
      Q2(11,12) = (27*35**(0.5D0)*MUH)/10240 
      Q2(11,13) = -(9*7**(0.5D0)*ETAH)/2048 
      Q2(11,14) = 0 
      Q2(11,15) = (27*35**(0.5D0)*ETAH)/10240 
      Q2(11,16) = 0 
      Q2(11,17) = (9*3**(0.5D0)*XIH)/2560 
      Q2(11,18) = 0 
      Q2(11,19) = -(27*15**(0.5D0)*XIH)/12800 
      Q2(11,20) = 0 
      Q2(11,21) = 0 
      Q2(11,22) = 0 
      Q2(11,23) = 0 
      Q2(11,24) = 0 
      Q2(11,25) = -(27*15**(0.5D0)*XIH)/12800 
      Q2(11,26) = 0 
      Q2(11,27) = (81*3**(0.5D0)*XIH)/12800 
      Q2(11,28) = 0 
      Q2(11,29) = 0 
      Q2(11,30) = 0 
      Q2(11,31) = 0 
      Q2(11,32) = 0 
      Q2(11,33) = (9*5**(0.5D0)*(6*ABS(ETAH) - 3*VT + 
     >   6*ABS(MUH) + 10*ABS(XIH)))/102400 
      Q2(11,34) = -(27*15**(0.5D0)*MUH)/25600 
      Q2(11,35) = (81*VT)/102400 - (81*ABS(ETAH))/51200 + 
     >   (27*ABS(MUH))/10240 - (27*ABS(XIH))/10240 
      Q2(11,36) = (27*35**(0.5D0)*MUH)/51200 
      Q2(11,37) = -(27*15**(0.5D0)*ETAH)/25600 
      Q2(11,38) = 0 
      Q2(11,39) = (81*3**(0.5D0)*ETAH)/25600 
      Q2(11,40) = 0 
      Q2(11,41) = (81*VT)/102400 + (27*ABS(ETAH))/10240 - 
     >   (81*ABS(MUH))/51200 - (27*ABS(XIH))/10240 
      Q2(11,42) = (81*3**(0.5D0)*MUH)/25600 
      Q2(11,43) = -(81*5**(0.5D0)*(3*VT + 10*ABS(ETAH) + 
     >   10*ABS(MUH) - 10*ABS(XIH)))/512000 
      Q2(11,44) = -(81*7**(0.5D0)*MUH)/51200 
      Q2(11,45) = (27*35**(0.5D0)*ETAH)/51200 
      Q2(11,46) = 0 
      Q2(11,47) = -(81*7**(0.5D0)*ETAH)/51200 
      Q2(11,48) = 0 
      Q2(11,49) = (9*7**(0.5D0)*XIH)/10240 
      Q2(11,50) = 0 
      Q2(11,51) = -(27*35**(0.5D0)*XIH)/51200 
      Q2(11,52) = 0 
      Q2(11,53) = 0 
      Q2(11,54) = 0 
      Q2(11,55) = 0 
      Q2(11,56) = 0 
      Q2(11,57) = -(27*35**(0.5D0)*XIH)/51200 
      Q2(11,58) = 0 
      Q2(11,59) = (81*7**(0.5D0)*XIH)/51200 
      Q2(11,60) = 0 
      Q2(11,61) = 0 
      Q2(11,62) = 0 
      Q2(11,63) = 0 
      Q2(11,64) = 0 

      Q2(12,1) = (27*35**(0.5D0)*MUH)/71680 
      Q2(12,2) = (9*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(35*VT - 
     >   70*ABS(ETAH) - 30*ABS(MUH) + 
     >   14*ABS(XIH)))/716800 
      Q2(12,3) = (27*7**(0.5D0)*MUH)/1792 
      Q2(12,4) = -(27*5**(0.5D0)*(45*VT - 90*ABS(ETAH) + 
     >   70*ABS(MUH) + 18*ABS(XIH)))/716800 
      Q2(12,5) = 0 
      Q2(12,6) = (27*35**(0.5D0)*ETAH)/5120 
      Q2(12,7) = 0 
      Q2(12,8) = -(243*15**(0.5D0)*ETAH)/35840 
      Q2(12,9) = -(81*7**(0.5D0)*MUH)/71680 
      Q2(12,10) = -(9*3**(0.5D0)*7**(0.5D0)*(105*VT + 
     >   350*ABS(ETAH) - 90*ABS(MUH) + 
     >   42*ABS(XIH)))/716800 
      Q2(12,11) = -(81*35**(0.5D0)*MUH)/8960 
      Q2(12,12) = (729*VT)/143360 + (243*ABS(ETAH))/14336 + 
     >   (81*ABS(MUH))/10240 + (729*ABS(XIH))/358400 
      Q2(12,13) = 0 
      Q2(12,14) = -(63*15**(0.5D0)*ETAH)/10240 
      Q2(12,15) = 0 
      Q2(12,16) = (243*35**(0.5D0)*ETAH)/71680 
      Q2(12,17) = 0 
      Q2(12,18) = (27*35**(0.5D0)*XIH)/12800 
      Q2(12,19) = 0 
      Q2(12,20) = -(243*15**(0.5D0)*XIH)/89600 
      Q2(12,21) = 0 
      Q2(12,22) = 0 
      Q2(12,23) = 0 
      Q2(12,24) = 0 
      Q2(12,25) = 0 
      Q2(12,26) = -(81*7**(0.5D0)*XIH)/12800 
      Q2(12,27) = 0 
      Q2(12,28) = (729*3**(0.5D0)*XIH)/89600 
      Q2(12,29) = 0 
      Q2(12,30) = 0 
      Q2(12,31) = 0 
      Q2(12,32) = 0 
      Q2(12,33) = -(81*7**(0.5D0)*MUH)/358400 
      Q2(12,34) = (9*3**(0.5D0)*7**(0.5D0)*(42*ABS(ETAH) - 
     >   21*VT + 18*ABS(MUH) + 70*ABS(XIH)))/716800 
      Q2(12,35) = -(81*35**(0.5D0)*MUH)/44800 
      Q2(12,36) = (729*VT)/716800 - (729*ABS(ETAH))/358400 + 
     >   (81*ABS(MUH))/51200 - (243*ABS(XIH))/71680 
      Q2(12,37) = 0 
      Q2(12,38) = -(81*7**(0.5D0)*ETAH)/25600 
      Q2(12,39) = 0 
      Q2(12,40) = (729*3**(0.5D0)*ETAH)/179200 
      Q2(12,41) = (243*35**(0.5D0)*MUH)/1792000 
      Q2(12,42) = (27*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(21*VT + 
     >   70*ABS(ETAH) - 18*ABS(MUH) - 
     >   70*ABS(XIH)))/3584000 
      Q2(12,43) = (243*7**(0.5D0)*MUH)/44800 
      Q2(12,44) = -(243*5**(0.5D0)*(9*VT + 30*ABS(ETAH) + 
     >   14*ABS(MUH) - 30*ABS(XIH)))/3584000 
      Q2(12,45) = 0 
      Q2(12,46) = (189*3**(0.5D0)*ETAH)/51200 
      Q2(12,47) = 0 
      Q2(12,48) = -(729*7**(0.5D0)*ETAH)/358400 
      Q2(12,49) = 0 
      Q2(12,50) = (63*15**(0.5D0)*XIH)/51200 
      Q2(12,51) = 0 
      Q2(12,52) = -(243*35**(0.5D0)*XIH)/358400 
      Q2(12,53) = 0 
      Q2(12,54) = 0 
      Q2(12,55) = 0 
      Q2(12,56) = 0 
      Q2(12,57) = 0 
      Q2(12,58) = -(189*3**(0.5D0)*XIH)/51200 
      Q2(12,59) = 0 
      Q2(12,60) = (729*7**(0.5D0)*XIH)/358400 
      Q2(12,61) = 0 
      Q2(12,62) = 0 
      Q2(12,63) = 0 
      Q2(12,64) = 0 

      Q2(13,1) = -(45*7**(0.5D0)*ETAH)/14336 
      Q2(13,2) = 0 
      Q2(13,3) = (27*35**(0.5D0)*ETAH)/71680 
      Q2(13,4) = 0 
      Q2(13,5) = -(3*3**(0.5D0)*7**(0.5D0)*(35*VT - 
     >   30*ABS(ETAH) + 14*ABS(MUH) + 
     >   14*ABS(XIH)))/28672 
      Q2(13,6) = -(9*7**(0.5D0)*MUH)/512 
      Q2(13,7) = (3*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(105*VT - 
     >   90*ABS(ETAH) - 350*ABS(MUH) + 
     >   42*ABS(XIH)))/716800 
      Q2(13,8) = -(21*3**(0.5D0)*MUH)/2048 
      Q2(13,9) = -(45*35**(0.5D0)*ETAH)/1792 
      Q2(13,10) = 0 
      Q2(13,11) = (27*7**(0.5D0)*ETAH)/1792 
      Q2(13,12) = 0 
      Q2(13,13) = (405*VT)/28672 + (45*ABS(ETAH))/2048 + 
     >   (81*ABS(MUH))/14336 + (81*ABS(XIH))/14336 
      Q2(13,14) = (81*3**(0.5D0)*MUH)/3584 
      Q2(13,15) = -(27*5**(0.5D0)*(45*VT + 70*ABS(ETAH) - 
     >   150*ABS(MUH) + 18*ABS(XIH)))/716800 
      Q2(13,16) = (81*7**(0.5D0)*MUH)/14336 
      Q2(13,17) = 0 
      Q2(13,18) = 0 
      Q2(13,19) = 0 
      Q2(13,20) = 0 
      Q2(13,21) = -(9*7**(0.5D0)*XIH)/512 
      Q2(13,22) = 0 
      Q2(13,23) = (27*35**(0.5D0)*XIH)/12800 
      Q2(13,24) = 0 
      Q2(13,25) = 0 
      Q2(13,26) = 0 
      Q2(13,27) = 0 
      Q2(13,28) = 0 
      Q2(13,29) = (81*3**(0.5D0)*XIH)/3584 
      Q2(13,30) = 0 
      Q2(13,31) = -(243*15**(0.5D0)*XIH)/89600 
      Q2(13,32) = 0 
      Q2(13,33) = (27*35**(0.5D0)*ETAH)/71680 
      Q2(13,34) = 0 
      Q2(13,35) = -(81*7**(0.5D0)*ETAH)/358400 
      Q2(13,36) = 0 
      Q2(13,37) = (3*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(105*VT - 
     >   90*ABS(ETAH) + 42*ABS(MUH) - 
     >   350*ABS(XIH)))/716800 
      Q2(13,38) = (27*35**(0.5D0)*MUH)/12800 
      Q2(13,39) = (9*3**(0.5D0)*7**(0.5D0)*(18*ABS(ETAH) - 
     >   21*VT + 70*ABS(MUH) + 70*ABS(XIH)))/716800 
      Q2(13,40) = (63*15**(0.5D0)*MUH)/51200 
      Q2(13,41) = (27*7**(0.5D0)*ETAH)/1792 
      Q2(13,42) = 0 
      Q2(13,43) = -(81*35**(0.5D0)*ETAH)/44800 
      Q2(13,44) = 0 
      Q2(13,45) = -(27*5**(0.5D0)*(45*VT + 70*ABS(ETAH) + 
     >   18*ABS(MUH) - 150*ABS(XIH)))/716800 
      Q2(13,46) = -(243*15**(0.5D0)*MUH)/89600 
      Q2(13,47) = (729*VT)/716800 + (81*ABS(ETAH))/51200 - 
     >   (243*ABS(MUH))/71680 - (243*ABS(XIH))/71680 
      Q2(13,48) = -(243*35**(0.5D0)*MUH)/358400 
      Q2(13,49) = 0 
      Q2(13,50) = 0 
      Q2(13,51) = 0 
      Q2(13,52) = 0 
      Q2(13,53) = -(21*3**(0.5D0)*XIH)/2048 
      Q2(13,54) = 0 
      Q2(13,55) = (63*15**(0.5D0)*XIH)/51200 
      Q2(13,56) = 0 
      Q2(13,57) = 0 
      Q2(13,58) = 0 
      Q2(13,59) = 0 
      Q2(13,60) = 0 
      Q2(13,61) = (81*7**(0.5D0)*XIH)/14336 
      Q2(13,62) = 0 
      Q2(13,63) = -(243*35**(0.5D0)*XIH)/358400 
      Q2(13,64) = 0 

      Q2(14,1) = 0 
      Q2(14,2) = -(51*7**(0.5D0)*ETAH)/14336 
      Q2(14,3) = 0 
      Q2(14,4) = (27*3**(0.5D0)*ETAH)/10240 
      Q2(14,5) = (33*7**(0.5D0)*MUH)/10240 
      Q2(14,6) = -(3**(0.5D0)*7**(0.5D0)*(595*VT - 
     >   510*ABS(ETAH) + 462*ABS(MUH) + 
     >   238*ABS(XIH)))/143360 
      Q2(14,7) = -(111*35**(0.5D0)*MUH)/5120 
      Q2(14,8) = (189*VT)/20480 - (81*ABS(ETAH))/10240 - 
     >   (231*ABS(MUH))/10240 + (189*ABS(XIH))/51200 
      Q2(14,9) = 0 
      Q2(14,10) = -(51*35**(0.5D0)*ETAH)/1792 
      Q2(14,11) = 0 
      Q2(14,12) = (27*15**(0.5D0)*ETAH)/1280 
      Q2(14,13) = -(297*3**(0.5D0)*MUH)/71680 
      Q2(14,14) = (459*VT)/28672 + (51*ABS(ETAH))/2048 + 
     >   (891*ABS(MUH))/71680 + (459*ABS(XIH))/71680 
      Q2(14,15) = (999*15**(0.5D0)*MUH)/35840 
      Q2(14,16) = -(27*3**(0.5D0)*7**(0.5D0)*(45*VT + 
     >   70*ABS(ETAH) - 110*ABS(MUH) + 
     >   18*ABS(XIH)))/716800 
      Q2(14,17) = 0 
      Q2(14,18) = 0 
      Q2(14,19) = 0 
      Q2(14,20) = 0 
      Q2(14,21) = 0 
      Q2(14,22) = -(51*7**(0.5D0)*XIH)/2560 
      Q2(14,23) = 0 
      Q2(14,24) = (189*3**(0.5D0)*XIH)/12800 
      Q2(14,25) = 0 
      Q2(14,26) = 0 
      Q2(14,27) = 0 
      Q2(14,28) = 0 
      Q2(14,29) = 0 
      Q2(14,30) = (459*3**(0.5D0)*XIH)/17920 
      Q2(14,31) = 0 
      Q2(14,32) = -(729*7**(0.5D0)*XIH)/89600 
      Q2(14,33) = 0 
      Q2(14,34) = (153*35**(0.5D0)*ETAH)/358400 
      Q2(14,35) = 0 
      Q2(14,36) = -(81*15**(0.5D0)*ETAH)/256000 
      Q2(14,37) = -(99*35**(0.5D0)*MUH)/256000 
      Q2(14,38) = (3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(1785*VT - 
     >   1530*ABS(ETAH) + 1386*ABS(MUH) - 
     >   5950*ABS(XIH)))/3584000 
      Q2(14,39) = (333*7**(0.5D0)*MUH)/25600 
      Q2(14,40) = (9*5**(0.5D0)*(54*ABS(ETAH) - 63*VT + 
     >   154*ABS(MUH) + 210*ABS(XIH)))/512000 
      Q2(14,41) = 0 
      Q2(14,42) = (153*7**(0.5D0)*ETAH)/8960 
      Q2(14,43) = 0 
      Q2(14,44) = -(81*3**(0.5D0)*ETAH)/6400 
      Q2(14,45) = (891*15**(0.5D0)*MUH)/1792000 
      Q2(14,46) = -(9*5**(0.5D0)*(765*VT + 1190*ABS(ETAH) + 
     >   594*ABS(MUH) - 2550*ABS(XIH)))/3584000 
      Q2(14,47) = -(2997*3**(0.5D0)*MUH)/179200 
      Q2(14,48) = (81*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(9*VT + 
     >   14*ABS(ETAH) - 22*ABS(MUH) - 
     >   30*ABS(XIH)))/3584000 
      Q2(14,49) = 0 
      Q2(14,50) = 0 
      Q2(14,51) = 0 
      Q2(14,52) = 0 
      Q2(14,53) = 0 
      Q2(14,54) = -(119*3**(0.5D0)*XIH)/10240 
      Q2(14,55) = 0 
      Q2(14,56) = (189*7**(0.5D0)*XIH)/51200 
      Q2(14,57) = 0 
      Q2(14,58) = 0 
      Q2(14,59) = 0 
      Q2(14,60) = 0 
      Q2(14,61) = 0 
      Q2(14,62) = (459*7**(0.5D0)*XIH)/71680 
      Q2(14,63) = 0 
      Q2(14,64) = -(243*3**(0.5D0)*XIH)/51200 

      Q2(15,1) = (27*35**(0.5D0)*ETAH)/71680 
      Q2(15,2) = 0 
      Q2(15,3) = -(81*7**(0.5D0)*ETAH)/71680 
      Q2(15,4) = 0 
      Q2(15,5) = (9*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(35*VT - 
     >   30*ABS(ETAH) - 70*ABS(MUH) + 
     >   14*ABS(XIH)))/716800 
      Q2(15,6) = (27*35**(0.5D0)*MUH)/5120 
      Q2(15,7) = -(9*3**(0.5D0)*7**(0.5D0)*(105*VT - 
     >   90*ABS(ETAH) + 350*ABS(MUH) + 
     >   42*ABS(XIH)))/716800 
      Q2(15,8) = -(63*15**(0.5D0)*MUH)/10240 
      Q2(15,9) = (27*7**(0.5D0)*ETAH)/1792 
      Q2(15,10) = 0 
      Q2(15,11) = -(81*35**(0.5D0)*ETAH)/8960 
      Q2(15,12) = 0 
      Q2(15,13) = -(27*5**(0.5D0)*(45*VT + 70*ABS(ETAH) - 
     >   90*ABS(MUH) + 18*ABS(XIH)))/716800 
      Q2(15,14) = -(243*15**(0.5D0)*MUH)/35840 
      Q2(15,15) = (729*VT)/143360 + (81*ABS(ETAH))/10240 + 
     >   (243*ABS(MUH))/14336 + (729*ABS(XIH))/358400 
      Q2(15,16) = (243*35**(0.5D0)*MUH)/71680 
      Q2(15,17) = 0 
      Q2(15,18) = 0 
      Q2(15,19) = 0 
      Q2(15,20) = 0 
      Q2(15,21) = (27*35**(0.5D0)*XIH)/12800 
      Q2(15,22) = 0 
      Q2(15,23) = -(81*7**(0.5D0)*XIH)/12800 
      Q2(15,24) = 0 
      Q2(15,25) = 0 
      Q2(15,26) = 0 
      Q2(15,27) = 0 
      Q2(15,28) = 0 
      Q2(15,29) = -(243*15**(0.5D0)*XIH)/89600 
      Q2(15,30) = 0 
      Q2(15,31) = (729*3**(0.5D0)*XIH)/89600 
      Q2(15,32) = 0 
      Q2(15,33) = -(81*7**(0.5D0)*ETAH)/358400 
      Q2(15,34) = 0 
      Q2(15,35) = (243*35**(0.5D0)*ETAH)/1792000 
      Q2(15,36) = 0 
      Q2(15,37) = (9*3**(0.5D0)*7**(0.5D0)*(18*ABS(ETAH) - 
     >   21*VT + 42*ABS(MUH) + 70*ABS(XIH)))/716800 
      Q2(15,38) = -(81*7**(0.5D0)*MUH)/25600 
      Q2(15,39) = (27*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(21*VT - 
     >   18*ABS(ETAH) + 70*ABS(MUH) - 
     >   70*ABS(XIH)))/3584000 
      Q2(15,40) = (189*3**(0.5D0)*MUH)/51200 
      Q2(15,41) = -(81*35**(0.5D0)*ETAH)/44800 
      Q2(15,42) = 0 
      Q2(15,43) = (243*7**(0.5D0)*ETAH)/44800 
      Q2(15,44) = 0 
      Q2(15,45) = (729*VT)/716800 + (81*ABS(ETAH))/51200 - 
     >   (729*ABS(MUH))/358400 - (243*ABS(XIH))/71680 
      Q2(15,46) = (729*3**(0.5D0)*MUH)/179200 
      Q2(15,47) = -(243*5**(0.5D0)*(9*VT + 14*ABS(ETAH) + 
     >   30*ABS(MUH) - 30*ABS(XIH)))/3584000 
      Q2(15,48) = -(729*7**(0.5D0)*MUH)/358400 
      Q2(15,49) = 0 
      Q2(15,50) = 0 
      Q2(15,51) = 0 
      Q2(15,52) = 0 
      Q2(15,53) = (63*15**(0.5D0)*XIH)/51200 
      Q2(15,54) = 0 
      Q2(15,55) = -(189*3**(0.5D0)*XIH)/51200 
      Q2(15,56) = 0 
      Q2(15,57) = 0 
      Q2(15,58) = 0 
      Q2(15,59) = 0 
      Q2(15,60) = 0 
      Q2(15,61) = -(243*35**(0.5D0)*XIH)/358400 
      Q2(15,62) = 0 
      Q2(15,63) = (729*7**(0.5D0)*XIH)/358400 
      Q2(15,64) = 0 

      Q2(16,1) = 0 
      Q2(16,2) = (27*3**(0.5D0)*ETAH)/10240 
      Q2(16,3) = 0 
      Q2(16,4) = -(729*7**(0.5D0)*ETAH)/501760 
      Q2(16,5) = (27*3**(0.5D0)*MUH)/10240 
      Q2(16,6) = (189*VT)/20480 - (81*ABS(ETAH))/10240 - 
     >   (81*ABS(MUH))/10240 + (189*ABS(XIH))/51200 
      Q2(16,7) = (27*15**(0.5D0)*MUH)/1280 
      Q2(16,8) = -(27*3**(0.5D0)*7**(0.5D0)*(315*VT - 
     >   270*ABS(ETAH) + 490*ABS(MUH) + 
     >   126*ABS(XIH)))/5017600 
      Q2(16,9) = 0 
      Q2(16,10) = (27*15**(0.5D0)*ETAH)/1280 
      Q2(16,11) = 0 
      Q2(16,12) = -(729*35**(0.5D0)*ETAH)/62720 
      Q2(16,13) = -(729*7**(0.5D0)*MUH)/501760 
      Q2(16,14) = -(27*3**(0.5D0)*7**(0.5D0)*(315*VT + 
     >   490*ABS(ETAH) - 270*ABS(MUH) + 
     >   126*ABS(XIH)))/5017600 
      Q2(16,15) = -(729*35**(0.5D0)*MUH)/62720 
      Q2(16,16) = (6561*VT)/1003520 + (729*ABS(ETAH))/71680 + 
     >   (729*ABS(MUH))/71680 + (6561*ABS(XIH))/2508800 
      Q2(16,17) = 0 
      Q2(16,18) = 0 
      Q2(16,19) = 0 
      Q2(16,20) = 0 
      Q2(16,21) = 0 
      Q2(16,22) = (189*3**(0.5D0)*XIH)/12800 
      Q2(16,23) = 0 
      Q2(16,24) = -(729*7**(0.5D0)*XIH)/89600 
      Q2(16,25) = 0 
      Q2(16,26) = 0 
      Q2(16,27) = 0 
      Q2(16,28) = 0 
      Q2(16,29) = 0 
      Q2(16,30) = -(729*7**(0.5D0)*XIH)/89600 
      Q2(16,31) = 0 
      Q2(16,32) = (6561*3**(0.5D0)*XIH)/627200 
      Q2(16,33) = 0 
      Q2(16,34) = -(81*15**(0.5D0)*ETAH)/256000 
      Q2(16,35) = 0 
      Q2(16,36) = (2187*35**(0.5D0)*ETAH)/12544000 
      Q2(16,37) = -(81*15**(0.5D0)*MUH)/256000 
      Q2(16,38) = (27*5**(0.5D0)*(18*ABS(ETAH) - 21*VT + 
     >   18*ABS(MUH) + 70*ABS(XIH)))/512000 
      Q2(16,39) = -(81*3**(0.5D0)*MUH)/6400 
      Q2(16,40) = (81*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(63*VT - 
     >   54*ABS(ETAH) + 98*ABS(MUH) - 
     >   210*ABS(XIH)))/25088000 
      Q2(16,41) = 0 
      Q2(16,42) = -(81*3**(0.5D0)*ETAH)/6400 
      Q2(16,43) = 0 
      Q2(16,44) = (2187*7**(0.5D0)*ETAH)/313600 
      Q2(16,45) = (2187*35**(0.5D0)*MUH)/12544000 
      Q2(16,46) = (81*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(63*VT + 
     >   98*ABS(ETAH) - 54*ABS(MUH) - 
     >   210*ABS(XIH)))/25088000 
      Q2(16,47) = (2187*7**(0.5D0)*MUH)/313600 
      Q2(16,48) = -(2187*5**(0.5D0)*(9*VT + 14*ABS(ETAH) + 
     >   14*ABS(MUH) - 30*ABS(XIH)))/25088000 
      Q2(16,49) = 0 
      Q2(16,50) = 0 
      Q2(16,51) = 0 
      Q2(16,52) = 0 
      Q2(16,53) = 0 
      Q2(16,54) = (189*7**(0.5D0)*XIH)/51200 
      Q2(16,55) = 0 
      Q2(16,56) = -(243*3**(0.5D0)*XIH)/51200 
      Q2(16,57) = 0 
      Q2(16,58) = 0 
      Q2(16,59) = 0 
      Q2(16,60) = 0 
      Q2(16,61) = 0 
      Q2(16,62) = -(243*3**(0.5D0)*XIH)/51200 
      Q2(16,63) = 0 
      Q2(16,64) = (6561*7**(0.5D0)*XIH)/2508800 

      Q2(17,1) = -(55*3**(0.5D0)*XIH)/6144 
      Q2(17,2) = 0 
      Q2(17,3) = (11*15**(0.5D0)*XIH)/10240 
      Q2(17,4) = 0 
      Q2(17,5) = 0 
      Q2(17,6) = 0 
      Q2(17,7) = 0 
      Q2(17,8) = 0 
      Q2(17,9) = (11*15**(0.5D0)*XIH)/10240 
      Q2(17,10) = 0 
      Q2(17,11) = -(33*3**(0.5D0)*XIH)/51200 
      Q2(17,12) = 0 
      Q2(17,13) = 0 
      Q2(17,14) = 0 
      Q2(17,15) = 0 
      Q2(17,16) = 0 
      Q2(17,17) = (425*VT)/12288 + (85*ABS(ETAH))/6144 + 
     >   (85*ABS(MUH))/6144 + (55*ABS(XIH))/2048 
      Q2(17,18) = (85*3**(0.5D0)*MUH)/1536 
      Q2(17,19) = -(5**(0.5D0)*(255*VT + 102*ABS(ETAH) - 
     >   850*ABS(MUH) + 198*ABS(XIH)))/61440 
      Q2(17,20) = (85*7**(0.5D0)*MUH)/6144 
      Q2(17,21) = (85*3**(0.5D0)*ETAH)/1536 
      Q2(17,22) = 0 
      Q2(17,23) = -(17*15**(0.5D0)*ETAH)/2560 
      Q2(17,24) = 0 
      Q2(17,25) = -(5**(0.5D0)*(255*VT - 850*ABS(ETAH) + 
     >   102*ABS(MUH) + 198*ABS(XIH)))/61440 
      Q2(17,26) = -(17*15**(0.5D0)*MUH)/2560 
      Q2(17,27) = (51*VT)/20480 - (17*ABS(ETAH))/2048 - 
     >   (17*ABS(MUH))/2048 + (99*ABS(XIH))/51200 
      Q2(17,28) = -(17*35**(0.5D0)*MUH)/10240 
      Q2(17,29) = (85*7**(0.5D0)*ETAH)/6144 
      Q2(17,30) = 0 
      Q2(17,31) = -(17*35**(0.5D0)*ETAH)/10240 
      Q2(17,32) = 0 
      Q2(17,33) = (185*15**(0.5D0)*XIH)/3072 
      Q2(17,34) = 0 
      Q2(17,35) = -(37*3**(0.5D0)*XIH)/1024 
      Q2(17,36) = 0 
      Q2(17,37) = 0 
      Q2(17,38) = 0 
      Q2(17,39) = 0 
      Q2(17,40) = 0 
      Q2(17,41) = -(37*3**(0.5D0)*XIH)/1024 
      Q2(17,42) = 0 
      Q2(17,43) = (111*15**(0.5D0)*XIH)/25600 
      Q2(17,44) = 0 
      Q2(17,45) = 0 
      Q2(17,46) = 0 
      Q2(17,47) = 0 
      Q2(17,48) = 0 
      Q2(17,49) = -(3**(0.5D0)*7**(0.5D0)*(45*VT + 
     >   18*ABS(ETAH) + 18*ABS(MUH) - 
     >   110*ABS(XIH)))/12288 
      Q2(17,50) = -(9*7**(0.5D0)*MUH)/512 
      Q2(17,51) = (3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(45*VT + 
     >   18*ABS(ETAH) - 150*ABS(MUH) - 
     >   110*ABS(XIH)))/102400 
      Q2(17,52) = -(21*3**(0.5D0)*MUH)/2048 
      Q2(17,53) = -(9*7**(0.5D0)*ETAH)/512 
      Q2(17,54) = 0 
      Q2(17,55) = (27*35**(0.5D0)*ETAH)/12800 
      Q2(17,56) = 0 
      Q2(17,57) = (3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(45*VT - 
     >   150*ABS(ETAH) + 18*ABS(MUH) - 
     >   110*ABS(XIH)))/102400 
      Q2(17,58) = (27*35**(0.5D0)*MUH)/12800 
      Q2(17,59) = (3*3**(0.5D0)*7**(0.5D0)*(30*ABS(ETAH) - 
     >   9*VT + 30*ABS(MUH) + 22*ABS(XIH)))/102400 
      Q2(17,60) = (63*15**(0.5D0)*MUH)/51200 
      Q2(17,61) = -(21*3**(0.5D0)*ETAH)/2048 
      Q2(17,62) = 0 
      Q2(17,63) = (63*15**(0.5D0)*ETAH)/51200 
      Q2(17,64) = 0 

      Q2(18,1) = 0 
      Q2(18,2) = -(187*3**(0.5D0)*XIH)/18432 
      Q2(18,3) = 0 
      Q2(18,4) = (33*7**(0.5D0)*XIH)/10240 
      Q2(18,5) = 0 
      Q2(18,6) = 0 
      Q2(18,7) = 0 
      Q2(18,8) = 0 
      Q2(18,9) = 0 
      Q2(18,10) = (187*15**(0.5D0)*XIH)/153600 
      Q2(18,11) = 0 
      Q2(18,12) = -(99*35**(0.5D0)*XIH)/256000 
      Q2(18,13) = 0 
      Q2(18,14) = 0 
      Q2(18,15) = 0 
      Q2(18,16) = 0 
      Q2(18,17) = -(187*3**(0.5D0)*MUH)/18432 
      Q2(18,18) = (1445*VT)/36864 + (289*ABS(ETAH))/18432 + 
     >   (187*ABS(MUH))/6144 + (187*ABS(XIH))/6144 
      Q2(18,19) = (629*15**(0.5D0)*MUH)/9216 
      Q2(18,20) = -(3**(0.5D0)*7**(0.5D0)*(765*VT + 
     >   306*ABS(ETAH) - 1870*ABS(MUH) + 
     >   594*ABS(XIH)))/184320 
      Q2(18,21) = 0 
      Q2(18,22) = (289*3**(0.5D0)*ETAH)/4608 
      Q2(18,23) = 0 
      Q2(18,24) = -(51*7**(0.5D0)*ETAH)/2560 
      Q2(18,25) = (187*15**(0.5D0)*MUH)/153600 
      Q2(18,26) = -(17*5**(0.5D0)*(255*VT - 850*ABS(ETAH) + 
     >   198*ABS(MUH) + 198*ABS(XIH)))/921600 
      Q2(18,27) = -(629*3**(0.5D0)*MUH)/15360 
      Q2(18,28) = (3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(765*VT - 
     >   2550*ABS(ETAH) - 1870*ABS(MUH) + 
     >   594*ABS(XIH)))/1536000 
      Q2(18,29) = 0 
      Q2(18,30) = (289*7**(0.5D0)*ETAH)/18432 
      Q2(18,31) = 0 
      Q2(18,32) = -(119*3**(0.5D0)*ETAH)/10240 
      Q2(18,33) = 0 
      Q2(18,34) = (629*15**(0.5D0)*XIH)/9216 
      Q2(18,35) = 0 
      Q2(18,36) = -(111*35**(0.5D0)*XIH)/5120 
      Q2(18,37) = 0 
      Q2(18,38) = 0 
      Q2(18,39) = 0 
      Q2(18,40) = 0 
      Q2(18,41) = 0 
      Q2(18,42) = -(629*3**(0.5D0)*XIH)/15360 
      Q2(18,43) = 0 
      Q2(18,44) = (333*7**(0.5D0)*XIH)/25600 
      Q2(18,45) = 0 
      Q2(18,46) = 0 
      Q2(18,47) = 0 
      Q2(18,48) = 0 
      Q2(18,49) = (33*7**(0.5D0)*MUH)/10240 
      Q2(18,50) = -(3**(0.5D0)*7**(0.5D0)*(765*VT + 
     >   306*ABS(ETAH) + 594*ABS(MUH) - 
     >   1870*ABS(XIH)))/184320 
      Q2(18,51) = -(111*35**(0.5D0)*MUH)/5120 
      Q2(18,52) = (189*VT)/20480 + (189*ABS(ETAH))/51200 - 
     >   (231*ABS(MUH))/10240 - (231*ABS(XIH))/10240 
      Q2(18,53) = 0 
      Q2(18,54) = -(51*7**(0.5D0)*ETAH)/2560 
      Q2(18,55) = 0 
      Q2(18,56) = (189*3**(0.5D0)*ETAH)/12800 
      Q2(18,57) = -(99*35**(0.5D0)*MUH)/256000 
      Q2(18,58) = (3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(765*VT - 
     >   2550*ABS(ETAH) + 594*ABS(MUH) - 
     >   1870*ABS(XIH)))/1536000 
      Q2(18,59) = (333*7**(0.5D0)*MUH)/25600 
      Q2(18,60) = (63*5**(0.5D0)*(30*ABS(ETAH) - 9*VT + 
     >   22*ABS(MUH) + 22*ABS(XIH)))/512000 
      Q2(18,61) = 0 
      Q2(18,62) = -(119*3**(0.5D0)*ETAH)/10240 
      Q2(18,63) = 0 
      Q2(18,64) = (189*7**(0.5D0)*ETAH)/51200 

      Q2(19,1) = (11*15**(0.5D0)*XIH)/10240 
      Q2(19,2) = 0 
      Q2(19,3) = -(33*3**(0.5D0)*XIH)/10240 
      Q2(19,4) = 0 
      Q2(19,5) = 0 
      Q2(19,6) = 0 
      Q2(19,7) = 0 
      Q2(19,8) = 0 
      Q2(19,9) = -(33*3**(0.5D0)*XIH)/51200 
      Q2(19,10) = 0 
      Q2(19,11) = (99*15**(0.5D0)*XIH)/256000 
      Q2(19,12) = 0 
      Q2(19,13) = 0 
      Q2(19,14) = 0 
      Q2(19,15) = 0 
      Q2(19,16) = 0 
      Q2(19,17) = -(5**(0.5D0)*(85*VT + 34*ABS(ETAH) - 
     >   170*ABS(MUH) + 66*ABS(XIH)))/20480 
      Q2(19,18) = -(17*15**(0.5D0)*MUH)/1024 
      Q2(19,19) = (51*VT)/4096 + (51*ABS(ETAH))/10240 + 
     >   (85*ABS(MUH))/2048 + (99*ABS(XIH))/10240 
      Q2(19,20) = (17*35**(0.5D0)*MUH)/2048 
      Q2(19,21) = -(17*15**(0.5D0)*ETAH)/2560 
      Q2(19,22) = 0 
      Q2(19,23) = (51*3**(0.5D0)*ETAH)/2560 
      Q2(19,24) = 0 
      Q2(19,25) = (51*VT)/20480 - (17*ABS(ETAH))/2048 - 
     >   (51*ABS(MUH))/10240 + (99*ABS(XIH))/51200 
      Q2(19,26) = (51*3**(0.5D0)*MUH)/5120 
      Q2(19,27) = -(3*5**(0.5D0)*(255*VT - 850*ABS(ETAH) + 
     >   850*ABS(MUH) + 198*ABS(XIH)))/512000 
      Q2(19,28) = -(51*7**(0.5D0)*MUH)/10240 
      Q2(19,29) = -(17*35**(0.5D0)*ETAH)/10240 
      Q2(19,30) = 0 
      Q2(19,31) = (51*7**(0.5D0)*ETAH)/10240 
      Q2(19,32) = 0 
      Q2(19,33) = -(37*3**(0.5D0)*XIH)/1024 
      Q2(19,34) = 0 
      Q2(19,35) = (111*15**(0.5D0)*XIH)/5120 
      Q2(19,36) = 0 
      Q2(19,37) = 0 
      Q2(19,38) = 0 
      Q2(19,39) = 0 
      Q2(19,40) = 0 
      Q2(19,41) = (111*15**(0.5D0)*XIH)/25600 
      Q2(19,42) = 0 
      Q2(19,43) = -(333*3**(0.5D0)*XIH)/25600 
      Q2(19,44) = 0 
      Q2(19,45) = 0 
      Q2(19,46) = 0 
      Q2(19,47) = 0 
      Q2(19,48) = 0 
      Q2(19,49) = (3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(45*VT + 
     >   18*ABS(ETAH) - 90*ABS(MUH) - 
     >   110*ABS(XIH)))/102400 
      Q2(19,50) = (27*35**(0.5D0)*MUH)/5120 
      Q2(19,51) = -(3*3**(0.5D0)*7**(0.5D0)*(45*VT + 
     >   18*ABS(ETAH) + 150*ABS(MUH) - 
     >   110*ABS(XIH)))/102400 
      Q2(19,52) = -(63*15**(0.5D0)*MUH)/10240 
      Q2(19,53) = (27*35**(0.5D0)*ETAH)/12800 
      Q2(19,54) = 0 
      Q2(19,55) = -(81*7**(0.5D0)*ETAH)/12800 
      Q2(19,56) = 0 
      Q2(19,57) = (3*3**(0.5D0)*7**(0.5D0)*(30*ABS(ETAH) - 
     >   9*VT + 18*ABS(MUH) + 22*ABS(XIH)))/102400 
      Q2(19,58) = -(81*7**(0.5D0)*MUH)/25600 
      Q2(19,59) = (9*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(9*VT - 
     >   30*ABS(ETAH) + 30*ABS(MUH) - 
     >   22*ABS(XIH)))/512000 
      Q2(19,60) = (189*3**(0.5D0)*MUH)/51200 
      Q2(19,61) = (63*15**(0.5D0)*ETAH)/51200 
      Q2(19,62) = 0 
      Q2(19,63) = -(189*3**(0.5D0)*ETAH)/51200 
      Q2(19,64) = 0 

      Q2(20,1) = 0 
      Q2(20,2) = (33*7**(0.5D0)*XIH)/10240 
      Q2(20,3) = 0 
      Q2(20,4) = -(297*3**(0.5D0)*XIH)/71680 
      Q2(20,5) = 0 
      Q2(20,6) = 0 
      Q2(20,7) = 0 
      Q2(20,8) = 0 
      Q2(20,9) = 0 
      Q2(20,10) = -(99*35**(0.5D0)*XIH)/256000 
      Q2(20,11) = 0 
      Q2(20,12) = (891*15**(0.5D0)*XIH)/1792000 
      Q2(20,13) = 0 
      Q2(20,14) = 0 
      Q2(20,15) = 0 
      Q2(20,16) = 0 
      Q2(20,17) = -(51*7**(0.5D0)*MUH)/14336 
      Q2(20,18) = -(3**(0.5D0)*7**(0.5D0)*(595*VT + 
     >   238*ABS(ETAH) - 510*ABS(MUH) + 
     >   462*ABS(XIH)))/143360 
      Q2(20,19) = -(51*35**(0.5D0)*MUH)/1792 
      Q2(20,20) = (459*VT)/28672 + (459*ABS(ETAH))/71680 + 
     >   (51*ABS(MUH))/2048 + (891*ABS(XIH))/71680 
      Q2(20,21) = 0 
      Q2(20,22) = -(51*7**(0.5D0)*ETAH)/2560 
      Q2(20,23) = 0 
      Q2(20,24) = (459*3**(0.5D0)*ETAH)/17920 
      Q2(20,25) = (153*35**(0.5D0)*MUH)/358400 
      Q2(20,26) = (3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(1785*VT - 
     >   5950*ABS(ETAH) - 1530*ABS(MUH) + 
     >   1386*ABS(XIH)))/3584000 
      Q2(20,27) = (153*7**(0.5D0)*MUH)/8960 
      Q2(20,28) = -(9*5**(0.5D0)*(765*VT - 2550*ABS(ETAH) + 
     >   1190*ABS(MUH) + 594*ABS(XIH)))/3584000 
      Q2(20,29) = 0 
      Q2(20,30) = -(119*3**(0.5D0)*ETAH)/10240 
      Q2(20,31) = 0 
      Q2(20,32) = (459*7**(0.5D0)*ETAH)/71680 
      Q2(20,33) = 0 
      Q2(20,34) = -(111*35**(0.5D0)*XIH)/5120 
      Q2(20,35) = 0 
      Q2(20,36) = (999*15**(0.5D0)*XIH)/35840 
      Q2(20,37) = 0 
      Q2(20,38) = 0 
      Q2(20,39) = 0 
      Q2(20,40) = 0 
      Q2(20,41) = 0 
      Q2(20,42) = (333*7**(0.5D0)*XIH)/25600 
      Q2(20,43) = 0 
      Q2(20,44) = -(2997*3**(0.5D0)*XIH)/179200 
      Q2(20,45) = 0 
      Q2(20,46) = 0 
      Q2(20,47) = 0 
      Q2(20,48) = 0 
      Q2(20,49) = (27*3**(0.5D0)*MUH)/10240 
      Q2(20,50) = (189*VT)/20480 + (189*ABS(ETAH))/51200 - 
     >   (81*ABS(MUH))/10240 - (231*ABS(XIH))/10240 
      Q2(20,51) = (27*15**(0.5D0)*MUH)/1280 
      Q2(20,52) = -(27*3**(0.5D0)*7**(0.5D0)*(45*VT + 
     >   18*ABS(ETAH) + 70*ABS(MUH) - 
     >   110*ABS(XIH)))/716800 
      Q2(20,53) = 0 
      Q2(20,54) = (189*3**(0.5D0)*ETAH)/12800 
      Q2(20,55) = 0 
      Q2(20,56) = -(729*7**(0.5D0)*ETAH)/89600 
      Q2(20,57) = -(81*15**(0.5D0)*MUH)/256000 
      Q2(20,58) = (9*5**(0.5D0)*(210*ABS(ETAH) - 63*VT + 
     >   54*ABS(MUH) + 154*ABS(XIH)))/512000 
      Q2(20,59) = -(81*3**(0.5D0)*MUH)/6400 
      Q2(20,60) = (81*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(9*VT - 
     >   30*ABS(ETAH) + 14*ABS(MUH) - 
     >   22*ABS(XIH)))/3584000 
      Q2(20,61) = 0 
      Q2(20,62) = (189*7**(0.5D0)*ETAH)/51200 
      Q2(20,63) = 0 
      Q2(20,64) = -(243*3**(0.5D0)*ETAH)/51200 

      Q2(21,1) = 0 
      Q2(21,2) = 0 
      Q2(21,3) = 0 
      Q2(21,4) = 0 
      Q2(21,5) = -(187*3**(0.5D0)*XIH)/18432 
      Q2(21,6) = 0 
      Q2(21,7) = (187*15**(0.5D0)*XIH)/153600 
      Q2(21,8) = 0 
      Q2(21,9) = 0 
      Q2(21,10) = 0 
      Q2(21,11) = 0 
      Q2(21,12) = 0 
      Q2(21,13) = (33*7**(0.5D0)*XIH)/10240 
      Q2(21,14) = 0 
      Q2(21,15) = -(99*35**(0.5D0)*XIH)/256000 
      Q2(21,16) = 0 
      Q2(21,17) = -(187*3**(0.5D0)*ETAH)/18432 
      Q2(21,18) = 0 
      Q2(21,19) = (187*15**(0.5D0)*ETAH)/153600 
      Q2(21,20) = 0 
      Q2(21,21) = (1445*VT)/36864 + (187*ABS(ETAH))/6144 + 
     >   (289*ABS(MUH))/18432 + (187*ABS(XIH))/6144 
      Q2(21,22) = (289*3**(0.5D0)*MUH)/4608 
      Q2(21,23) = -(17*5**(0.5D0)*(255*VT + 198*ABS(ETAH) - 
     >   850*ABS(MUH) + 198*ABS(XIH)))/921600 
      Q2(21,24) = (289*7**(0.5D0)*MUH)/18432 
      Q2(21,25) = (629*15**(0.5D0)*ETAH)/9216 
      Q2(21,26) = 0 
      Q2(21,27) = -(629*3**(0.5D0)*ETAH)/15360 
      Q2(21,28) = 0 
      Q2(21,29) = -(3**(0.5D0)*7**(0.5D0)*(765*VT - 
     >   1870*ABS(ETAH) + 306*ABS(MUH) + 
     >   594*ABS(XIH)))/184320 
      Q2(21,30) = -(51*7**(0.5D0)*MUH)/2560 
      Q2(21,31) = (3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(765*VT - 
     >   1870*ABS(ETAH) - 2550*ABS(MUH) + 
     >   594*ABS(XIH)))/1536000 
      Q2(21,32) = -(119*3**(0.5D0)*MUH)/10240 
      Q2(21,33) = 0 
      Q2(21,34) = 0 
      Q2(21,35) = 0 
      Q2(21,36) = 0 
      Q2(21,37) = (629*15**(0.5D0)*XIH)/9216 
      Q2(21,38) = 0 
      Q2(21,39) = -(629*3**(0.5D0)*XIH)/15360 
      Q2(21,40) = 0 
      Q2(21,41) = 0 
      Q2(21,42) = 0 
      Q2(21,43) = 0 
      Q2(21,44) = 0 
      Q2(21,45) = -(111*35**(0.5D0)*XIH)/5120 
      Q2(21,46) = 0 
      Q2(21,47) = (333*7**(0.5D0)*XIH)/25600 
      Q2(21,48) = 0 
      Q2(21,49) = (33*7**(0.5D0)*ETAH)/10240 
      Q2(21,50) = 0 
      Q2(21,51) = -(99*35**(0.5D0)*ETAH)/256000 
      Q2(21,52) = 0 
      Q2(21,53) = -(3**(0.5D0)*7**(0.5D0)*(765*VT + 
     >   594*ABS(ETAH) + 306*ABS(MUH) - 
     >   1870*ABS(XIH)))/184320 
      Q2(21,54) = -(51*7**(0.5D0)*MUH)/2560 
      Q2(21,55) = (3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(765*VT + 
     >   594*ABS(ETAH) - 2550*ABS(MUH) - 
     >   1870*ABS(XIH)))/1536000 
      Q2(21,56) = -(119*3**(0.5D0)*MUH)/10240 
      Q2(21,57) = -(111*35**(0.5D0)*ETAH)/5120 
      Q2(21,58) = 0 
      Q2(21,59) = (333*7**(0.5D0)*ETAH)/25600 
      Q2(21,60) = 0 
      Q2(21,61) = (189*VT)/20480 - (231*ABS(ETAH))/10240 + 
     >   (189*ABS(MUH))/51200 - (231*ABS(XIH))/10240 
      Q2(21,62) = (189*3**(0.5D0)*MUH)/12800 
      Q2(21,63) = (63*5**(0.5D0)*(22*ABS(ETAH) - 9*VT + 
     >   30*ABS(MUH) + 22*ABS(XIH)))/512000 
      Q2(21,64) = (189*7**(0.5D0)*MUH)/51200 

      Q2(22,1) = 0 
      Q2(22,2) = 0 
      Q2(22,3) = 0 
      Q2(22,4) = 0 
      Q2(22,5) = 0 
      Q2(22,6) = -(3179*3**(0.5D0)*XIH)/276480 
      Q2(22,7) = 0 
      Q2(22,8) = (187*7**(0.5D0)*XIH)/51200 
      Q2(22,9) = 0 
      Q2(22,10) = 0 
      Q2(22,11) = 0 
      Q2(22,12) = 0 
      Q2(22,13) = 0 
      Q2(22,14) = (187*7**(0.5D0)*XIH)/51200 
      Q2(22,15) = 0 
      Q2(22,16) = -(693*3**(0.5D0)*XIH)/256000 
      Q2(22,17) = 0 
      Q2(22,18) = -(3179*3**(0.5D0)*ETAH)/276480 
      Q2(22,19) = 0 
      Q2(22,20) = (187*7**(0.5D0)*ETAH)/51200 
      Q2(22,21) = -(3179*3**(0.5D0)*MUH)/276480 
      Q2(22,22) = (4913*VT)/110592 + (3179*ABS(ETAH))/92160 + 
     >   (3179*ABS(MUH))/92160 + (3179*ABS(XIH))/92160 
      Q2(22,23) = (10693*15**(0.5D0)*MUH)/138240 
      Q2(22,24) = -(17*3**(0.5D0)*7**(0.5D0)*(765*VT + 
     >   594*ABS(ETAH) - 1870*ABS(MUH) + 
     >   594*ABS(XIH)))/2764800 
      Q2(22,25) = 0 
      Q2(22,26) = (10693*15**(0.5D0)*ETAH)/138240 
      Q2(22,27) = 0 
      Q2(22,28) = -(629*35**(0.5D0)*ETAH)/25600 
      Q2(22,29) = (187*7**(0.5D0)*MUH)/51200 
      Q2(22,30) = -(17*3**(0.5D0)*7**(0.5D0)*(765*VT - 
     >   1870*ABS(ETAH) + 594*ABS(MUH) + 
     >   594*ABS(XIH)))/2764800 
      Q2(22,31) = -(629*35**(0.5D0)*MUH)/25600 
      Q2(22,32) = (1071*VT)/102400 - (1309*ABS(ETAH))/51200 - 
     >   (1309*ABS(MUH))/51200 + (2079*ABS(XIH))/256000 
      Q2(22,33) = 0 
      Q2(22,34) = 0 
      Q2(22,35) = 0 
      Q2(22,36) = 0 
      Q2(22,37) = 0 
      Q2(22,38) = (10693*15**(0.5D0)*XIH)/138240 
      Q2(22,39) = 0 
      Q2(22,40) = -(629*35**(0.5D0)*XIH)/25600 
      Q2(22,41) = 0 
      Q2(22,42) = 0 
      Q2(22,43) = 0 
      Q2(22,44) = 0 
      Q2(22,45) = 0 
      Q2(22,46) = -(629*35**(0.5D0)*XIH)/25600 
      Q2(22,47) = 0 
      Q2(22,48) = (2331*15**(0.5D0)*XIH)/128000 
      Q2(22,49) = 0 
      Q2(22,50) = (187*7**(0.5D0)*ETAH)/51200 
      Q2(22,51) = 0 
      Q2(22,52) = -(693*3**(0.5D0)*ETAH)/256000 
      Q2(22,53) = (187*7**(0.5D0)*MUH)/51200 
      Q2(22,54) = -(17*3**(0.5D0)*7**(0.5D0)*(765*VT + 
     >   594*ABS(ETAH) + 594*ABS(MUH) - 
     >   1870*ABS(XIH)))/2764800 
      Q2(22,55) = -(629*35**(0.5D0)*MUH)/25600 
      Q2(22,56) = (1071*VT)/102400 + (2079*ABS(ETAH))/256000 - 
     >   (1309*ABS(MUH))/51200 - (1309*ABS(XIH))/51200 
      Q2(22,57) = 0 
      Q2(22,58) = -(629*35**(0.5D0)*ETAH)/25600 
      Q2(22,59) = 0 
      Q2(22,60) = (2331*15**(0.5D0)*ETAH)/128000 
      Q2(22,61) = -(693*3**(0.5D0)*MUH)/256000 
      Q2(22,62) = (1071*VT)/102400 - (1309*ABS(ETAH))/51200 + 
     >   (2079*ABS(MUH))/256000 - (1309*ABS(XIH))/51200 
      Q2(22,63) = (2331*15**(0.5D0)*MUH)/128000 
      Q2(22,64) = (63*3**(0.5D0)*7**(0.5D0)*(22*ABS(ETAH) - 
     >   9*VT + 22*ABS(MUH) + 22*ABS(XIH)))/512000 

      Q2(23,1) = 0 
      Q2(23,2) = 0 
      Q2(23,3) = 0 
      Q2(23,4) = 0 
      Q2(23,5) = (187*15**(0.5D0)*XIH)/153600 
      Q2(23,6) = 0 
      Q2(23,7) = -(187*3**(0.5D0)*XIH)/51200 
      Q2(23,8) = 0 
      Q2(23,9) = 0 
      Q2(23,10) = 0 
      Q2(23,11) = 0 
      Q2(23,12) = 0 
      Q2(23,13) = -(99*35**(0.5D0)*XIH)/256000 
      Q2(23,14) = 0 
      Q2(23,15) = (297*7**(0.5D0)*XIH)/256000 
      Q2(23,16) = 0 
      Q2(23,17) = (187*15**(0.5D0)*ETAH)/153600 
      Q2(23,18) = 0 
      Q2(23,19) = -(187*3**(0.5D0)*ETAH)/51200 
      Q2(23,20) = 0 
      Q2(23,21) = -(17*5**(0.5D0)*(85*VT + 66*ABS(ETAH) - 
     >   170*ABS(MUH) + 66*ABS(XIH)))/307200 
      Q2(23,22) = -(289*15**(0.5D0)*MUH)/15360 
      Q2(23,23) = (289*VT)/20480 + (561*ABS(ETAH))/51200 + 
     >   (289*ABS(MUH))/6144 + (561*ABS(XIH))/51200 
      Q2(23,24) = (289*35**(0.5D0)*MUH)/30720 
      Q2(23,25) = -(629*3**(0.5D0)*ETAH)/15360 
      Q2(23,26) = 0 
      Q2(23,27) = (629*15**(0.5D0)*ETAH)/25600 
      Q2(23,28) = 0 
      Q2(23,29) = (3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(765*VT - 
     >   1870*ABS(ETAH) - 1530*ABS(MUH) + 
     >   594*ABS(XIH)))/1536000 
      Q2(23,30) = (153*35**(0.5D0)*MUH)/25600 
      Q2(23,31) = -(3**(0.5D0)*7**(0.5D0)*(765*VT - 
     >   1870*ABS(ETAH) + 2550*ABS(MUH) + 
     >   594*ABS(XIH)))/512000 
      Q2(23,32) = -(357*15**(0.5D0)*MUH)/51200 
      Q2(23,33) = 0 
      Q2(23,34) = 0 
      Q2(23,35) = 0 
      Q2(23,36) = 0 
      Q2(23,37) = -(629*3**(0.5D0)*XIH)/15360 
      Q2(23,38) = 0 
      Q2(23,39) = (629*15**(0.5D0)*XIH)/25600 
      Q2(23,40) = 0 
      Q2(23,41) = 0 
      Q2(23,42) = 0 
      Q2(23,43) = 0 
      Q2(23,44) = 0 
      Q2(23,45) = (333*7**(0.5D0)*XIH)/25600 
      Q2(23,46) = 0 
      Q2(23,47) = -(999*35**(0.5D0)*XIH)/128000 
      Q2(23,48) = 0 
      Q2(23,49) = -(99*35**(0.5D0)*ETAH)/256000 
      Q2(23,50) = 0 
      Q2(23,51) = (297*7**(0.5D0)*ETAH)/256000 
      Q2(23,52) = 0 
      Q2(23,53) = (3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(765*VT + 
     >   594*ABS(ETAH) - 1530*ABS(MUH) - 
     >   1870*ABS(XIH)))/1536000 
      Q2(23,54) = (153*35**(0.5D0)*MUH)/25600 
      Q2(23,55) = -(3**(0.5D0)*7**(0.5D0)*(765*VT + 
     >   594*ABS(ETAH) + 2550*ABS(MUH) - 
     >   1870*ABS(XIH)))/512000 
      Q2(23,56) = -(357*15**(0.5D0)*MUH)/51200 
      Q2(23,57) = (333*7**(0.5D0)*ETAH)/25600 
      Q2(23,58) = 0 
      Q2(23,59) = -(999*35**(0.5D0)*ETAH)/128000 
      Q2(23,60) = 0 
      Q2(23,61) = (63*5**(0.5D0)*(22*ABS(ETAH) - 9*VT + 
     >   18*ABS(MUH) + 22*ABS(XIH)))/512000 
      Q2(23,62) = -(567*15**(0.5D0)*MUH)/128000 
      Q2(23,63) = (1701*VT)/512000 - (2079*ABS(ETAH))/256000 + 
     >   (567*ABS(MUH))/51200 - (2079*ABS(XIH))/256000 
      Q2(23,64) = (567*35**(0.5D0)*MUH)/256000 

      Q2(24,1) = 0 
      Q2(24,2) = 0 
      Q2(24,3) = 0 
      Q2(24,4) = 0 
      Q2(24,5) = 0 
      Q2(24,6) = (187*7**(0.5D0)*XIH)/51200 
      Q2(24,7) = 0 
      Q2(24,8) = -(1683*3**(0.5D0)*XIH)/358400 
      Q2(24,9) = 0 
      Q2(24,10) = 0 
      Q2(24,11) = 0 
      Q2(24,12) = 0 
      Q2(24,13) = 0 
      Q2(24,14) = -(693*3**(0.5D0)*XIH)/256000 
      Q2(24,15) = 0 
      Q2(24,16) = (2673*7**(0.5D0)*XIH)/1792000 
      Q2(24,17) = 0 
      Q2(24,18) = (187*7**(0.5D0)*ETAH)/51200 
      Q2(24,19) = 0 
      Q2(24,20) = -(1683*3**(0.5D0)*ETAH)/358400 
      Q2(24,21) = -(289*7**(0.5D0)*MUH)/71680 
      Q2(24,22) = -(17*3**(0.5D0)*7**(0.5D0)*(595*VT + 
     >   462*ABS(ETAH) - 510*ABS(MUH) + 
     >   462*ABS(XIH)))/2150400 
      Q2(24,23) = -(289*35**(0.5D0)*MUH)/8960 
      Q2(24,24) = (2601*VT)/143360 + (5049*ABS(ETAH))/358400 + 
     >   (289*ABS(MUH))/10240 + (5049*ABS(XIH))/358400 
      Q2(24,25) = 0 
      Q2(24,26) = -(629*35**(0.5D0)*ETAH)/25600 
      Q2(24,27) = 0 
      Q2(24,28) = (5661*15**(0.5D0)*ETAH)/179200 
      Q2(24,29) = (153*3**(0.5D0)*MUH)/51200 
      Q2(24,30) = (1071*VT)/102400 - (1309*ABS(ETAH))/51200 - 
     >   (459*ABS(MUH))/51200 + (2079*ABS(XIH))/256000 
      Q2(24,31) = (153*15**(0.5D0)*MUH)/6400 
      Q2(24,32) = -(9*3**(0.5D0)*7**(0.5D0)*(765*VT - 
     >   1870*ABS(ETAH) + 1190*ABS(MUH) + 
     >   594*ABS(XIH)))/3584000 
      Q2(24,33) = 0 
      Q2(24,34) = 0 
      Q2(24,35) = 0 
      Q2(24,36) = 0 
      Q2(24,37) = 0 
      Q2(24,38) = -(629*35**(0.5D0)*XIH)/25600 
      Q2(24,39) = 0 
      Q2(24,40) = (5661*15**(0.5D0)*XIH)/179200 
      Q2(24,41) = 0 
      Q2(24,42) = 0 
      Q2(24,43) = 0 
      Q2(24,44) = 0 
      Q2(24,45) = 0 
      Q2(24,46) = (2331*15**(0.5D0)*XIH)/128000 
      Q2(24,47) = 0 
      Q2(24,48) = -(8991*35**(0.5D0)*XIH)/896000 
      Q2(24,49) = 0 
      Q2(24,50) = -(693*3**(0.5D0)*ETAH)/256000 
      Q2(24,51) = 0 
      Q2(24,52) = (2673*7**(0.5D0)*ETAH)/1792000 
      Q2(24,53) = (153*3**(0.5D0)*MUH)/51200 
      Q2(24,54) = (1071*VT)/102400 + (2079*ABS(ETAH))/256000 - 
     >   (459*ABS(MUH))/51200 - (1309*ABS(XIH))/51200 
      Q2(24,55) = (153*15**(0.5D0)*MUH)/6400 
      Q2(24,56) = -(9*3**(0.5D0)*7**(0.5D0)*(765*VT + 
     >   594*ABS(ETAH) + 1190*ABS(MUH) - 
     >   1870*ABS(XIH)))/3584000 
      Q2(24,57) = 0 
      Q2(24,58) = (2331*15**(0.5D0)*ETAH)/128000 
      Q2(24,59) = 0 
      Q2(24,60) = -(8991*35**(0.5D0)*ETAH)/896000 
      Q2(24,61) = -(243*7**(0.5D0)*MUH)/256000 
      Q2(24,62) = (9*3**(0.5D0)*7**(0.5D0)*(154*ABS(ETAH) - 
     >   63*VT + 54*ABS(MUH) + 154*ABS(XIH)))/512000 
      Q2(24,63) = -(243*35**(0.5D0)*MUH)/32000 
      Q2(24,64) = (2187*VT)/512000 - (2673*ABS(ETAH))/256000 + 
     >   (1701*ABS(MUH))/256000 - (2673*ABS(XIH))/256000 

      Q2(25,1) = (11*15**(0.5D0)*XIH)/10240 
      Q2(25,2) = 0 
      Q2(25,3) = -(33*3**(0.5D0)*XIH)/51200 
      Q2(25,4) = 0 
      Q2(25,5) = 0 
      Q2(25,6) = 0 
      Q2(25,7) = 0 
      Q2(25,8) = 0 
      Q2(25,9) = -(33*3**(0.5D0)*XIH)/10240 
      Q2(25,10) = 0 
      Q2(25,11) = (99*15**(0.5D0)*XIH)/256000 
      Q2(25,12) = 0 
      Q2(25,13) = 0 
      Q2(25,14) = 0 
      Q2(25,15) = 0 
      Q2(25,16) = 0 
      Q2(25,17) = -(5**(0.5D0)*(85*VT - 170*ABS(ETAH) + 
     >   34*ABS(MUH) + 66*ABS(XIH)))/20480 
      Q2(25,18) = -(17*15**(0.5D0)*MUH)/2560 
      Q2(25,19) = (51*VT)/20480 - (51*ABS(ETAH))/10240 - 
     >   (17*ABS(MUH))/2048 + (99*ABS(XIH))/51200 
      Q2(25,20) = -(17*35**(0.5D0)*MUH)/10240 
      Q2(25,21) = -(17*15**(0.5D0)*ETAH)/1024 
      Q2(25,22) = 0 
      Q2(25,23) = (51*3**(0.5D0)*ETAH)/5120 
      Q2(25,24) = 0 
      Q2(25,25) = (51*VT)/4096 + (85*ABS(ETAH))/2048 + 
     >   (51*ABS(MUH))/10240 + (99*ABS(XIH))/10240 
      Q2(25,26) = (51*3**(0.5D0)*MUH)/2560 
      Q2(25,27) = -(3*5**(0.5D0)*(255*VT + 850*ABS(ETAH) - 
     >   850*ABS(MUH) + 198*ABS(XIH)))/512000 
      Q2(25,28) = (51*7**(0.5D0)*MUH)/10240 
      Q2(25,29) = (17*35**(0.5D0)*ETAH)/2048 
      Q2(25,30) = 0 
      Q2(25,31) = -(51*7**(0.5D0)*ETAH)/10240 
      Q2(25,32) = 0 
      Q2(25,33) = -(37*3**(0.5D0)*XIH)/1024 
      Q2(25,34) = 0 
      Q2(25,35) = (111*15**(0.5D0)*XIH)/25600 
      Q2(25,36) = 0 
      Q2(25,37) = 0 
      Q2(25,38) = 0 
      Q2(25,39) = 0 
      Q2(25,40) = 0 
      Q2(25,41) = (111*15**(0.5D0)*XIH)/5120 
      Q2(25,42) = 0 
      Q2(25,43) = -(333*3**(0.5D0)*XIH)/25600 
      Q2(25,44) = 0 
      Q2(25,45) = 0 
      Q2(25,46) = 0 
      Q2(25,47) = 0 
      Q2(25,48) = 0 
      Q2(25,49) = (3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(45*VT - 
     >   90*ABS(ETAH) + 18*ABS(MUH) - 
     >   110*ABS(XIH)))/102400 
      Q2(25,50) = (27*35**(0.5D0)*MUH)/12800 
      Q2(25,51) = (3*3**(0.5D0)*7**(0.5D0)*(18*ABS(ETAH) - 
     >   9*VT + 30*ABS(MUH) + 22*ABS(XIH)))/102400 
      Q2(25,52) = (63*15**(0.5D0)*MUH)/51200 
      Q2(25,53) = (27*35**(0.5D0)*ETAH)/5120 
      Q2(25,54) = 0 
      Q2(25,55) = -(81*7**(0.5D0)*ETAH)/25600 
      Q2(25,56) = 0 
      Q2(25,57) = -(3*3**(0.5D0)*7**(0.5D0)*(45*VT + 
     >   150*ABS(ETAH) + 18*ABS(MUH) - 
     >   110*ABS(XIH)))/102400 
      Q2(25,58) = -(81*7**(0.5D0)*MUH)/12800 
      Q2(25,59) = (9*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(9*VT + 
     >   30*ABS(ETAH) - 30*ABS(MUH) - 
     >   22*ABS(XIH)))/512000 
      Q2(25,60) = -(189*3**(0.5D0)*MUH)/51200 
      Q2(25,61) = -(63*15**(0.5D0)*ETAH)/10240 
      Q2(25,62) = 0 
      Q2(25,63) = (189*3**(0.5D0)*ETAH)/51200 
      Q2(25,64) = 0 

      Q2(26,1) = 0 
      Q2(26,2) = (187*15**(0.5D0)*XIH)/153600 
      Q2(26,3) = 0 
      Q2(26,4) = -(99*35**(0.5D0)*XIH)/256000 
      Q2(26,5) = 0 
      Q2(26,6) = 0 
      Q2(26,7) = 0 
      Q2(26,8) = 0 
      Q2(26,9) = 0 
      Q2(26,10) = -(187*3**(0.5D0)*XIH)/51200 
      Q2(26,11) = 0 
      Q2(26,12) = (297*7**(0.5D0)*XIH)/256000 
      Q2(26,13) = 0 
      Q2(26,14) = 0 
      Q2(26,15) = 0 
      Q2(26,16) = 0 
      Q2(26,17) = (187*15**(0.5D0)*MUH)/153600 
      Q2(26,18) = -(17*5**(0.5D0)*(85*VT - 170*ABS(ETAH) + 
     >   66*ABS(MUH) + 66*ABS(XIH)))/307200 
      Q2(26,19) = -(629*3**(0.5D0)*MUH)/15360 
      Q2(26,20) = (3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(765*VT - 
     >   1530*ABS(ETAH) - 1870*ABS(MUH) + 
     >   594*ABS(XIH)))/1536000 
      Q2(26,21) = 0 
      Q2(26,22) = -(289*15**(0.5D0)*ETAH)/15360 
      Q2(26,23) = 0 
      Q2(26,24) = (153*35**(0.5D0)*ETAH)/25600 
      Q2(26,25) = -(187*3**(0.5D0)*MUH)/51200 
      Q2(26,26) = (289*VT)/20480 + (289*ABS(ETAH))/6144 + 
     >   (561*ABS(MUH))/51200 + (561*ABS(XIH))/51200 
      Q2(26,27) = (629*15**(0.5D0)*MUH)/25600 
      Q2(26,28) = -(3**(0.5D0)*7**(0.5D0)*(765*VT + 
     >   2550*ABS(ETAH) - 1870*ABS(MUH) + 
     >   594*ABS(XIH)))/512000 
      Q2(26,29) = 0 
      Q2(26,30) = (289*35**(0.5D0)*ETAH)/30720 
      Q2(26,31) = 0 
      Q2(26,32) = -(357*15**(0.5D0)*ETAH)/51200 
      Q2(26,33) = 0 
      Q2(26,34) = -(629*3**(0.5D0)*XIH)/15360 
      Q2(26,35) = 0 
      Q2(26,36) = (333*7**(0.5D0)*XIH)/25600 
      Q2(26,37) = 0 
      Q2(26,38) = 0 
      Q2(26,39) = 0 
      Q2(26,40) = 0 
      Q2(26,41) = 0 
      Q2(26,42) = (629*15**(0.5D0)*XIH)/25600 
      Q2(26,43) = 0 
      Q2(26,44) = -(999*35**(0.5D0)*XIH)/128000 
      Q2(26,45) = 0 
      Q2(26,46) = 0 
      Q2(26,47) = 0 
      Q2(26,48) = 0 
      Q2(26,49) = -(99*35**(0.5D0)*MUH)/256000 
      Q2(26,50) = (3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(765*VT - 
     >   1530*ABS(ETAH) + 594*ABS(MUH) - 
     >   1870*ABS(XIH)))/1536000 
      Q2(26,51) = (333*7**(0.5D0)*MUH)/25600 
      Q2(26,52) = (63*5**(0.5D0)*(18*ABS(ETAH) - 9*VT + 
     >   22*ABS(MUH) + 22*ABS(XIH)))/512000 
      Q2(26,53) = 0 
      Q2(26,54) = (153*35**(0.5D0)*ETAH)/25600 
      Q2(26,55) = 0 
      Q2(26,56) = -(567*15**(0.5D0)*ETAH)/128000 
      Q2(26,57) = (297*7**(0.5D0)*MUH)/256000 
      Q2(26,58) = -(3**(0.5D0)*7**(0.5D0)*(765*VT + 
     >   2550*ABS(ETAH) + 594*ABS(MUH) - 
     >   1870*ABS(XIH)))/512000 
      Q2(26,59) = -(999*35**(0.5D0)*MUH)/128000 
      Q2(26,60) = (1701*VT)/512000 + (567*ABS(ETAH))/51200 - 
     >   (2079*ABS(MUH))/256000 - (2079*ABS(XIH))/256000 
      Q2(26,61) = 0 
      Q2(26,62) = -(357*15**(0.5D0)*ETAH)/51200 
      Q2(26,63) = 0 
      Q2(26,64) = (567*35**(0.5D0)*ETAH)/256000 

      Q2(27,1) = -(33*3**(0.5D0)*XIH)/51200 
      Q2(27,2) = 0 
      Q2(27,3) = (99*15**(0.5D0)*XIH)/256000 
      Q2(27,4) = 0 
      Q2(27,5) = 0 
      Q2(27,6) = 0 
      Q2(27,7) = 0 
      Q2(27,8) = 0 
      Q2(27,9) = (99*15**(0.5D0)*XIH)/256000 
      Q2(27,10) = 0 
      Q2(27,11) = -(297*3**(0.5D0)*XIH)/256000 
      Q2(27,12) = 0 
      Q2(27,13) = 0 
      Q2(27,14) = 0 
      Q2(27,15) = 0 
      Q2(27,16) = 0 
      Q2(27,17) = (51*VT)/20480 - (51*ABS(ETAH))/10240 - 
     >   (51*ABS(MUH))/10240 + (99*ABS(XIH))/51200 
      Q2(27,18) = (51*3**(0.5D0)*MUH)/5120 
      Q2(27,19) = -(3*5**(0.5D0)*(255*VT - 510*ABS(ETAH) + 
     >   850*ABS(MUH) + 198*ABS(XIH)))/512000 
      Q2(27,20) = -(51*7**(0.5D0)*MUH)/10240 
      Q2(27,21) = (51*3**(0.5D0)*ETAH)/5120 
      Q2(27,22) = 0 
      Q2(27,23) = -(153*15**(0.5D0)*ETAH)/25600 
      Q2(27,24) = 0 
      Q2(27,25) = -(3*5**(0.5D0)*(255*VT + 850*ABS(ETAH) - 
     >   510*ABS(MUH) + 198*ABS(XIH)))/512000 
      Q2(27,26) = -(153*15**(0.5D0)*MUH)/25600 
      Q2(27,27) = (459*VT)/102400 + (153*ABS(ETAH))/10240 + 
     >   (153*ABS(MUH))/10240 + (891*ABS(XIH))/256000 
      Q2(27,28) = (153*35**(0.5D0)*MUH)/51200 
      Q2(27,29) = -(51*7**(0.5D0)*ETAH)/10240 
      Q2(27,30) = 0 
      Q2(27,31) = (153*35**(0.5D0)*ETAH)/51200 
      Q2(27,32) = 0 
      Q2(27,33) = (111*15**(0.5D0)*XIH)/25600 
      Q2(27,34) = 0 
      Q2(27,35) = -(333*3**(0.5D0)*XIH)/25600 
      Q2(27,36) = 0 
      Q2(27,37) = 0 
      Q2(27,38) = 0 
      Q2(27,39) = 0 
      Q2(27,40) = 0 
      Q2(27,41) = -(333*3**(0.5D0)*XIH)/25600 
      Q2(27,42) = 0 
      Q2(27,43) = (999*15**(0.5D0)*XIH)/128000 
      Q2(27,44) = 0 
      Q2(27,45) = 0 
      Q2(27,46) = 0 
      Q2(27,47) = 0 
      Q2(27,48) = 0 
      Q2(27,49) = (3*3**(0.5D0)*7**(0.5D0)*(18*ABS(ETAH) - 
     >   9*VT + 18*ABS(MUH) + 22*ABS(XIH)))/102400 
      Q2(27,50) = -(81*7**(0.5D0)*MUH)/25600 
      Q2(27,51) = (9*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(9*VT - 
     >   18*ABS(ETAH) + 30*ABS(MUH) - 
     >   22*ABS(XIH)))/512000 
      Q2(27,52) = (189*3**(0.5D0)*MUH)/51200 
      Q2(27,53) = -(81*7**(0.5D0)*ETAH)/25600 
      Q2(27,54) = 0 
      Q2(27,55) = (243*35**(0.5D0)*ETAH)/128000 
      Q2(27,56) = 0 
      Q2(27,57) = (9*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(9*VT + 
     >   30*ABS(ETAH) - 18*ABS(MUH) - 
     >   22*ABS(XIH)))/512000 
      Q2(27,58) = (243*35**(0.5D0)*MUH)/128000 
      Q2(27,59) = -(27*3**(0.5D0)*7**(0.5D0)*(9*VT + 
     >   30*ABS(ETAH) + 30*ABS(MUH) - 
     >   22*ABS(XIH)))/512000 
      Q2(27,60) = -(567*15**(0.5D0)*MUH)/256000 
      Q2(27,61) = (189*3**(0.5D0)*ETAH)/51200 
      Q2(27,62) = 0 
      Q2(27,63) = -(567*15**(0.5D0)*ETAH)/256000 
      Q2(27,64) = 0 

      Q2(28,1) = 0 
      Q2(28,2) = -(99*35**(0.5D0)*XIH)/256000 
      Q2(28,3) = 0 
      Q2(28,4) = (891*15**(0.5D0)*XIH)/1792000 
      Q2(28,5) = 0 
      Q2(28,6) = 0 
      Q2(28,7) = 0 
      Q2(28,8) = 0 
      Q2(28,9) = 0 
      Q2(28,10) = (297*7**(0.5D0)*XIH)/256000 
      Q2(28,11) = 0 
      Q2(28,12) = -(2673*3**(0.5D0)*XIH)/1792000 
      Q2(28,13) = 0 
      Q2(28,14) = 0 
      Q2(28,15) = 0 
      Q2(28,16) = 0 
      Q2(28,17) = (153*35**(0.5D0)*MUH)/358400 
      Q2(28,18) = (3*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(595*VT - 
     >   1190*ABS(ETAH) - 510*ABS(MUH) + 
     >   462*ABS(XIH)))/3584000 
      Q2(28,19) = (153*7**(0.5D0)*MUH)/8960 
      Q2(28,20) = -(9*5**(0.5D0)*(765*VT - 1530*ABS(ETAH) + 
     >   1190*ABS(MUH) + 594*ABS(XIH)))/3584000 
      Q2(28,21) = 0 
      Q2(28,22) = (153*35**(0.5D0)*ETAH)/25600 
      Q2(28,23) = 0 
      Q2(28,24) = -(1377*15**(0.5D0)*ETAH)/179200 
      Q2(28,25) = -(459*7**(0.5D0)*MUH)/358400 
      Q2(28,26) = -(3*3**(0.5D0)*7**(0.5D0)*(1785*VT + 
     >   5950*ABS(ETAH) - 1530*ABS(MUH) + 
     >   1386*ABS(XIH)))/3584000 
      Q2(28,27) = -(459*35**(0.5D0)*MUH)/44800 
      Q2(28,28) = (4131*VT)/716800 + (1377*ABS(ETAH))/71680 + 
     >   (459*ABS(MUH))/51200 + (8019*ABS(XIH))/1792000 
      Q2(28,29) = 0 
      Q2(28,30) = -(357*15**(0.5D0)*ETAH)/51200 
      Q2(28,31) = 0 
      Q2(28,32) = (1377*35**(0.5D0)*ETAH)/358400 
      Q2(28,33) = 0 
      Q2(28,34) = (333*7**(0.5D0)*XIH)/25600 
      Q2(28,35) = 0 
      Q2(28,36) = -(2997*3**(0.5D0)*XIH)/179200 
      Q2(28,37) = 0 
      Q2(28,38) = 0 
      Q2(28,39) = 0 
      Q2(28,40) = 0 
      Q2(28,41) = 0 
      Q2(28,42) = -(999*35**(0.5D0)*XIH)/128000 
      Q2(28,43) = 0 
      Q2(28,44) = (8991*15**(0.5D0)*XIH)/896000 
      Q2(28,45) = 0 
      Q2(28,46) = 0 
      Q2(28,47) = 0 
      Q2(28,48) = 0 
      Q2(28,49) = -(81*15**(0.5D0)*MUH)/256000 
      Q2(28,50) = (9*5**(0.5D0)*(126*ABS(ETAH) - 63*VT + 
     >   54*ABS(MUH) + 154*ABS(XIH)))/512000 
      Q2(28,51) = -(81*3**(0.5D0)*MUH)/6400 
      Q2(28,52) = (81*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(9*VT - 
     >   18*ABS(ETAH) + 14*ABS(MUH) - 
     >   22*ABS(XIH)))/3584000 
      Q2(28,53) = 0 
      Q2(28,54) = -(567*15**(0.5D0)*ETAH)/128000 
      Q2(28,55) = 0 
      Q2(28,56) = (2187*35**(0.5D0)*ETAH)/896000 
      Q2(28,57) = (243*3**(0.5D0)*MUH)/256000 
      Q2(28,58) = (1701*VT)/512000 + (567*ABS(ETAH))/51200 - 
     >   (729*ABS(MUH))/256000 - (2079*ABS(XIH))/256000 
      Q2(28,59) = (243*15**(0.5D0)*MUH)/32000 
      Q2(28,60) = -(243*3**(0.5D0)*7**(0.5D0)*(9*VT + 
     >   30*ABS(ETAH) + 14*ABS(MUH) - 
     >   22*ABS(XIH)))/3584000 
      Q2(28,61) = 0 
      Q2(28,62) = (567*35**(0.5D0)*ETAH)/256000 
      Q2(28,63) = 0 
      Q2(28,64) = -(729*15**(0.5D0)*ETAH)/256000 

      Q2(29,1) = 0 
      Q2(29,2) = 0 
      Q2(29,3) = 0 
      Q2(29,4) = 0 
      Q2(29,5) = (33*7**(0.5D0)*XIH)/10240 
      Q2(29,6) = 0 
      Q2(29,7) = -(99*35**(0.5D0)*XIH)/256000 
      Q2(29,8) = 0 
      Q2(29,9) = 0 
      Q2(29,10) = 0 
      Q2(29,11) = 0 
      Q2(29,12) = 0 
      Q2(29,13) = -(297*3**(0.5D0)*XIH)/71680 
      Q2(29,14) = 0 
      Q2(29,15) = (891*15**(0.5D0)*XIH)/1792000 
      Q2(29,16) = 0 
      Q2(29,17) = -(51*7**(0.5D0)*ETAH)/14336 
      Q2(29,18) = 0 
      Q2(29,19) = (153*35**(0.5D0)*ETAH)/358400 
      Q2(29,20) = 0 
      Q2(29,21) = -(3**(0.5D0)*7**(0.5D0)*(595*VT - 
     >   510*ABS(ETAH) + 238*ABS(MUH) + 
     >   462*ABS(XIH)))/143360 
      Q2(29,22) = -(51*7**(0.5D0)*MUH)/2560 
      Q2(29,23) = (3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(1785*VT - 
     >   1530*ABS(ETAH) - 5950*ABS(MUH) + 
     >   1386*ABS(XIH)))/3584000 
      Q2(29,24) = -(119*3**(0.5D0)*MUH)/10240 
      Q2(29,25) = -(51*35**(0.5D0)*ETAH)/1792 
      Q2(29,26) = 0 
      Q2(29,27) = (153*7**(0.5D0)*ETAH)/8960 
      Q2(29,28) = 0 
      Q2(29,29) = (459*VT)/28672 + (51*ABS(ETAH))/2048 + 
     >   (459*ABS(MUH))/71680 + (891*ABS(XIH))/71680 
      Q2(29,30) = (459*3**(0.5D0)*MUH)/17920 
      Q2(29,31) = -(9*5**(0.5D0)*(765*VT + 1190*ABS(ETAH) - 
     >   2550*ABS(MUH) + 594*ABS(XIH)))/3584000 
      Q2(29,32) = (459*7**(0.5D0)*MUH)/71680 
      Q2(29,33) = 0 
      Q2(29,34) = 0 
      Q2(29,35) = 0 
      Q2(29,36) = 0 
      Q2(29,37) = -(111*35**(0.5D0)*XIH)/5120 
      Q2(29,38) = 0 
      Q2(29,39) = (333*7**(0.5D0)*XIH)/25600 
      Q2(29,40) = 0 
      Q2(29,41) = 0 
      Q2(29,42) = 0 
      Q2(29,43) = 0 
      Q2(29,44) = 0 
      Q2(29,45) = (999*15**(0.5D0)*XIH)/35840 
      Q2(29,46) = 0 
      Q2(29,47) = -(2997*3**(0.5D0)*XIH)/179200 
      Q2(29,48) = 0 
      Q2(29,49) = (27*3**(0.5D0)*ETAH)/10240 
      Q2(29,50) = 0 
      Q2(29,51) = -(81*15**(0.5D0)*ETAH)/256000 
      Q2(29,52) = 0 
      Q2(29,53) = (189*VT)/20480 - (81*ABS(ETAH))/10240 + 
     >   (189*ABS(MUH))/51200 - (231*ABS(XIH))/10240 
      Q2(29,54) = (189*3**(0.5D0)*MUH)/12800 
      Q2(29,55) = (9*5**(0.5D0)*(54*ABS(ETAH) - 63*VT + 
     >   210*ABS(MUH) + 154*ABS(XIH)))/512000 
      Q2(29,56) = (189*7**(0.5D0)*MUH)/51200 
      Q2(29,57) = (27*15**(0.5D0)*ETAH)/1280 
      Q2(29,58) = 0 
      Q2(29,59) = -(81*3**(0.5D0)*ETAH)/6400 
      Q2(29,60) = 0 
      Q2(29,61) = -(27*3**(0.5D0)*7**(0.5D0)*(45*VT + 
     >   70*ABS(ETAH) + 18*ABS(MUH) - 
     >   110*ABS(XIH)))/716800 
      Q2(29,62) = -(729*7**(0.5D0)*MUH)/89600 
      Q2(29,63) = (81*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(9*VT + 
     >   14*ABS(ETAH) - 30*ABS(MUH) - 
     >   22*ABS(XIH)))/3584000 
      Q2(29,64) = -(243*3**(0.5D0)*MUH)/51200 

      Q2(30,1) = 0 
      Q2(30,2) = 0 
      Q2(30,3) = 0 
      Q2(30,4) = 0 
      Q2(30,5) = 0 
      Q2(30,6) = (187*7**(0.5D0)*XIH)/51200 
      Q2(30,7) = 0 
      Q2(30,8) = -(693*3**(0.5D0)*XIH)/256000 
      Q2(30,9) = 0 
      Q2(30,10) = 0 
      Q2(30,11) = 0 
      Q2(30,12) = 0 
      Q2(30,13) = 0 
      Q2(30,14) = -(1683*3**(0.5D0)*XIH)/358400 
      Q2(30,15) = 0 
      Q2(30,16) = (2673*7**(0.5D0)*XIH)/1792000 
      Q2(30,17) = 0 
      Q2(30,18) = -(289*7**(0.5D0)*ETAH)/71680 
      Q2(30,19) = 0 
      Q2(30,20) = (153*3**(0.5D0)*ETAH)/51200 
      Q2(30,21) = (187*7**(0.5D0)*MUH)/51200 
      Q2(30,22) = -(17*3**(0.5D0)*7**(0.5D0)*(595*VT - 
     >   510*ABS(ETAH) + 462*ABS(MUH) + 
     >   462*ABS(XIH)))/2150400 
      Q2(30,23) = -(629*35**(0.5D0)*MUH)/25600 
      Q2(30,24) = (1071*VT)/102400 - (459*ABS(ETAH))/51200 - 
     >   (1309*ABS(MUH))/51200 + (2079*ABS(XIH))/256000 
      Q2(30,25) = 0 
      Q2(30,26) = -(289*35**(0.5D0)*ETAH)/8960 
      Q2(30,27) = 0 
      Q2(30,28) = (153*15**(0.5D0)*ETAH)/6400 
      Q2(30,29) = -(1683*3**(0.5D0)*MUH)/358400 
      Q2(30,30) = (2601*VT)/143360 + (289*ABS(ETAH))/10240 + 
     >   (5049*ABS(MUH))/358400 + (5049*ABS(XIH))/358400 
      Q2(30,31) = (5661*15**(0.5D0)*MUH)/179200 
      Q2(30,32) = -(9*3**(0.5D0)*7**(0.5D0)*(765*VT + 
     >   1190*ABS(ETAH) - 1870*ABS(MUH) + 
     >   594*ABS(XIH)))/3584000 
      Q2(30,33) = 0 
      Q2(30,34) = 0 
      Q2(30,35) = 0 
      Q2(30,36) = 0 
      Q2(30,37) = 0 
      Q2(30,38) = -(629*35**(0.5D0)*XIH)/25600 
      Q2(30,39) = 0 
      Q2(30,40) = (2331*15**(0.5D0)*XIH)/128000 
      Q2(30,41) = 0 
      Q2(30,42) = 0 
      Q2(30,43) = 0 
      Q2(30,44) = 0 
      Q2(30,45) = 0 
      Q2(30,46) = (5661*15**(0.5D0)*XIH)/179200 
      Q2(30,47) = 0 
      Q2(30,48) = -(8991*35**(0.5D0)*XIH)/896000 
      Q2(30,49) = 0 
      Q2(30,50) = (153*3**(0.5D0)*ETAH)/51200 
      Q2(30,51) = 0 
      Q2(30,52) = -(243*7**(0.5D0)*ETAH)/256000 
      Q2(30,53) = -(693*3**(0.5D0)*MUH)/256000 
      Q2(30,54) = (1071*VT)/102400 - (459*ABS(ETAH))/51200 + 
     >   (2079*ABS(MUH))/256000 - (1309*ABS(XIH))/51200 
      Q2(30,55) = (2331*15**(0.5D0)*MUH)/128000 
      Q2(30,56) = (9*3**(0.5D0)*7**(0.5D0)*(54*ABS(ETAH) - 
     >   63*VT + 154*ABS(MUH) + 154*ABS(XIH)))/512000 
      Q2(30,57) = 0 
      Q2(30,58) = (153*15**(0.5D0)*ETAH)/6400 
      Q2(30,59) = 0 
      Q2(30,60) = -(243*35**(0.5D0)*ETAH)/32000 
      Q2(30,61) = (2673*7**(0.5D0)*MUH)/1792000 
      Q2(30,62) = -(9*3**(0.5D0)*7**(0.5D0)*(765*VT + 
     >   1190*ABS(ETAH) + 594*ABS(MUH) - 
     >   1870*ABS(XIH)))/3584000 
      Q2(30,63) = -(8991*35**(0.5D0)*MUH)/896000 
      Q2(30,64) = (2187*VT)/512000 + (1701*ABS(ETAH))/256000 - 
     >   (2673*ABS(MUH))/256000 - (2673*ABS(XIH))/256000 

      Q2(31,1) = 0 
      Q2(31,2) = 0 
      Q2(31,3) = 0 
      Q2(31,4) = 0 
      Q2(31,5) = -(99*35**(0.5D0)*XIH)/256000 
      Q2(31,6) = 0 
      Q2(31,7) = (297*7**(0.5D0)*XIH)/256000 
      Q2(31,8) = 0 
      Q2(31,9) = 0 
      Q2(31,10) = 0 
      Q2(31,11) = 0 
      Q2(31,12) = 0 
      Q2(31,13) = (891*15**(0.5D0)*XIH)/1792000 
      Q2(31,14) = 0 
      Q2(31,15) = -(2673*3**(0.5D0)*XIH)/1792000 
      Q2(31,16) = 0 
      Q2(31,17) = (153*35**(0.5D0)*ETAH)/358400 
      Q2(31,18) = 0 
      Q2(31,19) = -(459*7**(0.5D0)*ETAH)/358400 
      Q2(31,20) = 0 
      Q2(31,21) = (3*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(595*VT - 
     >   510*ABS(ETAH) - 1190*ABS(MUH) + 
     >   462*ABS(XIH)))/3584000 
      Q2(31,22) = (153*35**(0.5D0)*MUH)/25600 
      Q2(31,23) = -(3*3**(0.5D0)*7**(0.5D0)*(1785*VT - 
     >   1530*ABS(ETAH) + 5950*ABS(MUH) + 
     >   1386*ABS(XIH)))/3584000 
      Q2(31,24) = -(357*15**(0.5D0)*MUH)/51200 
      Q2(31,25) = (153*7**(0.5D0)*ETAH)/8960 
      Q2(31,26) = 0 
      Q2(31,27) = -(459*35**(0.5D0)*ETAH)/44800 
      Q2(31,28) = 0 
      Q2(31,29) = -(9*5**(0.5D0)*(765*VT + 1190*ABS(ETAH) - 
     >   1530*ABS(MUH) + 594*ABS(XIH)))/3584000 
      Q2(31,30) = -(1377*15**(0.5D0)*MUH)/179200 
      Q2(31,31) = (4131*VT)/716800 + (459*ABS(ETAH))/51200 + 
     >   (1377*ABS(MUH))/71680 + (8019*ABS(XIH))/1792000 
      Q2(31,32) = (1377*35**(0.5D0)*MUH)/358400 
      Q2(31,33) = 0 
      Q2(31,34) = 0 
      Q2(31,35) = 0 
      Q2(31,36) = 0 
      Q2(31,37) = (333*7**(0.5D0)*XIH)/25600 
      Q2(31,38) = 0 
      Q2(31,39) = -(999*35**(0.5D0)*XIH)/128000 
      Q2(31,40) = 0 
      Q2(31,41) = 0 
      Q2(31,42) = 0 
      Q2(31,43) = 0 
      Q2(31,44) = 0 
      Q2(31,45) = -(2997*3**(0.5D0)*XIH)/179200 
      Q2(31,46) = 0 
      Q2(31,47) = (8991*15**(0.5D0)*XIH)/896000 
      Q2(31,48) = 0 
      Q2(31,49) = -(81*15**(0.5D0)*ETAH)/256000 
      Q2(31,50) = 0 
      Q2(31,51) = (243*3**(0.5D0)*ETAH)/256000 
      Q2(31,52) = 0 
      Q2(31,53) = (9*5**(0.5D0)*(54*ABS(ETAH) - 63*VT + 
     >   126*ABS(MUH) + 154*ABS(XIH)))/512000 
      Q2(31,54) = -(567*15**(0.5D0)*MUH)/128000 
      Q2(31,55) = (1701*VT)/512000 - (729*ABS(ETAH))/256000 + 
     >   (567*ABS(MUH))/51200 - (2079*ABS(XIH))/256000 
      Q2(31,56) = (567*35**(0.5D0)*MUH)/256000 
      Q2(31,57) = -(81*3**(0.5D0)*ETAH)/6400 
      Q2(31,58) = 0 
      Q2(31,59) = (243*15**(0.5D0)*ETAH)/32000 
      Q2(31,60) = 0 
      Q2(31,61) = (81*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(9*VT + 
     >   14*ABS(ETAH) - 18*ABS(MUH) - 
     >   22*ABS(XIH)))/3584000 
      Q2(31,62) = (2187*35**(0.5D0)*MUH)/896000 
      Q2(31,63) = -(243*3**(0.5D0)*7**(0.5D0)*(9*VT + 
     >   14*ABS(ETAH) + 30*ABS(MUH) - 
     >   22*ABS(XIH)))/3584000 
      Q2(31,64) = -(729*15**(0.5D0)*MUH)/256000 

      Q2(32,1) = 0 
      Q2(32,2) = 0 
      Q2(32,3) = 0 
      Q2(32,4) = 0 
      Q2(32,5) = 0 
      Q2(32,6) = -(693*3**(0.5D0)*XIH)/256000 
      Q2(32,7) = 0 
      Q2(32,8) = (2673*7**(0.5D0)*XIH)/1792000 
      Q2(32,9) = 0 
      Q2(32,10) = 0 
      Q2(32,11) = 0 
      Q2(32,12) = 0 
      Q2(32,13) = 0 
      Q2(32,14) = (2673*7**(0.5D0)*XIH)/1792000 
      Q2(32,15) = 0 
      Q2(32,16) = -(24057*3**(0.5D0)*XIH)/12544000 
      Q2(32,17) = 0 
      Q2(32,18) = (153*3**(0.5D0)*ETAH)/51200 
      Q2(32,19) = 0 
      Q2(32,20) = -(4131*7**(0.5D0)*ETAH)/2508800 
      Q2(32,21) = (153*3**(0.5D0)*MUH)/51200 
      Q2(32,22) = (1071*VT)/102400 - (459*ABS(ETAH))/51200 - 
     >   (459*ABS(MUH))/51200 + (2079*ABS(XIH))/256000 
      Q2(32,23) = (153*15**(0.5D0)*MUH)/6400 
      Q2(32,24) = -(9*3**(0.5D0)*7**(0.5D0)*(5355*VT - 
     >   4590*ABS(ETAH) + 8330*ABS(MUH) + 
     >   4158*ABS(XIH)))/25088000 
      Q2(32,25) = 0 
      Q2(32,26) = (153*15**(0.5D0)*ETAH)/6400 
      Q2(32,27) = 0 
      Q2(32,28) = -(4131*35**(0.5D0)*ETAH)/313600 
      Q2(32,29) = -(4131*7**(0.5D0)*MUH)/2508800 
      Q2(32,30) = -(9*3**(0.5D0)*7**(0.5D0)*(5355*VT + 
     >   8330*ABS(ETAH) - 4590*ABS(MUH) + 
     >   4158*ABS(XIH)))/25088000 
      Q2(32,31) = -(4131*35**(0.5D0)*MUH)/313600 
      Q2(32,32) = (37179*VT)/5017600 + 
     >   (4131*ABS(ETAH))/358400 + 
     >   (4131*ABS(MUH))/358400 + 
     >   (72171*ABS(XIH))/12544000 
      Q2(32,33) = 0 
      Q2(32,34) = 0 
      Q2(32,35) = 0 
      Q2(32,36) = 0 
      Q2(32,37) = 0 
      Q2(32,38) = (2331*15**(0.5D0)*XIH)/128000 
      Q2(32,39) = 0 
      Q2(32,40) = -(8991*35**(0.5D0)*XIH)/896000 
      Q2(32,41) = 0 
      Q2(32,42) = 0 
      Q2(32,43) = 0 
      Q2(32,44) = 0 
      Q2(32,45) = 0 
      Q2(32,46) = -(8991*35**(0.5D0)*XIH)/896000 
      Q2(32,47) = 0 
      Q2(32,48) = (80919*15**(0.5D0)*XIH)/6272000 
      Q2(32,49) = 0 
      Q2(32,50) = -(243*7**(0.5D0)*ETAH)/256000 
      Q2(32,51) = 0 
      Q2(32,52) = (2187*3**(0.5D0)*ETAH)/1792000 
      Q2(32,53) = -(243*7**(0.5D0)*MUH)/256000 
      Q2(32,54) = (9*3**(0.5D0)*7**(0.5D0)*(54*ABS(ETAH) - 
     >   63*VT + 54*ABS(MUH) + 154*ABS(XIH)))/512000 
      Q2(32,55) = -(243*35**(0.5D0)*MUH)/32000 
      Q2(32,56) = (2187*VT)/512000 - (6561*ABS(ETAH))/1792000 + 
     >   (1701*ABS(MUH))/256000 - (2673*ABS(XIH))/256000 
      Q2(32,57) = 0 
      Q2(32,58) = -(243*35**(0.5D0)*ETAH)/32000 
      Q2(32,59) = 0 
      Q2(32,60) = (2187*15**(0.5D0)*ETAH)/224000 
      Q2(32,61) = (2187*3**(0.5D0)*MUH)/1792000 
      Q2(32,62) = (2187*VT)/512000 + (1701*ABS(ETAH))/256000 - 
     >   (6561*ABS(MUH))/1792000 - (2673*ABS(XIH))/256000 
      Q2(32,63) = (2187*15**(0.5D0)*MUH)/224000 
      Q2(32,64) = -(2187*3**(0.5D0)*7**(0.5D0)*(9*VT + 
     >   14*ABS(ETAH) + 14*ABS(MUH) - 
     >   22*ABS(XIH)))/25088000 

      Q2(33,1) = -(3*5**(0.5D0)*(5*VT + 2*ABS(ETAH) + 
     >   2*ABS(MUH) - 10*ABS(XIH)))/4096 
      Q2(33,2) = -(3*15**(0.5D0)*MUH)/512 
      Q2(33,3) = (9*VT)/4096 + (9*ABS(ETAH))/10240 - 
     >   (15*ABS(MUH))/2048 - (9*ABS(XIH))/2048 
      Q2(33,4) = -(3*35**(0.5D0)*MUH)/2048 
      Q2(33,5) = -(3*15**(0.5D0)*ETAH)/512 
      Q2(33,6) = 0 
      Q2(33,7) = (9*3**(0.5D0)*ETAH)/2560 
      Q2(33,8) = 0 
      Q2(33,9) = (9*VT)/4096 - (15*ABS(ETAH))/2048 + 
     >   (9*ABS(MUH))/10240 - (9*ABS(XIH))/2048 
      Q2(33,10) = (9*3**(0.5D0)*MUH)/2560 
      Q2(33,11) = (9*5**(0.5D0)*(10*ABS(ETAH) - 3*VT + 
     >   10*ABS(MUH) + 6*ABS(XIH)))/102400 
      Q2(33,12) = (9*7**(0.5D0)*MUH)/10240 
      Q2(33,13) = -(3*35**(0.5D0)*ETAH)/2048 
      Q2(33,14) = 0 
      Q2(33,15) = (9*7**(0.5D0)*ETAH)/10240 
      Q2(33,16) = 0 
      Q2(33,17) = -(15*15**(0.5D0)*XIH)/1024 
      Q2(33,18) = 0 
      Q2(33,19) = (9*3**(0.5D0)*XIH)/1024 
      Q2(33,20) = 0 
      Q2(33,21) = 0 
      Q2(33,22) = 0 
      Q2(33,23) = 0 
      Q2(33,24) = 0 
      Q2(33,25) = (9*3**(0.5D0)*XIH)/1024 
      Q2(33,26) = 0 
      Q2(33,27) = -(27*15**(0.5D0)*XIH)/25600 
      Q2(33,28) = 0 
      Q2(33,29) = 0 
      Q2(33,30) = 0 
      Q2(33,31) = 0 
      Q2(33,32) = 0 
      Q2(33,33) = (45*VT)/4096 + (9*ABS(ETAH))/2048 + 
     >   (9*ABS(MUH))/2048 + (75*ABS(XIH))/2048 
      Q2(33,34) = (9*3**(0.5D0)*MUH)/512 
      Q2(33,35) = -(9*5**(0.5D0)*(15*VT + 6*ABS(ETAH) - 
     >   50*ABS(MUH) + 50*ABS(XIH)))/102400 
      Q2(33,36) = (9*7**(0.5D0)*MUH)/2048 
      Q2(33,37) = (9*3**(0.5D0)*ETAH)/512 
      Q2(33,38) = 0 
      Q2(33,39) = -(27*15**(0.5D0)*ETAH)/12800 
      Q2(33,40) = 0 
      Q2(33,41) = -(9*5**(0.5D0)*(15*VT - 50*ABS(ETAH) + 
     >   6*ABS(MUH) + 50*ABS(XIH)))/102400 
      Q2(33,42) = -(27*15**(0.5D0)*MUH)/12800 
      Q2(33,43) = (81*VT)/102400 - (27*ABS(ETAH))/10240 - 
     >   (27*ABS(MUH))/10240 + (27*ABS(XIH))/10240 
      Q2(33,44) = -(27*35**(0.5D0)*MUH)/51200 
      Q2(33,45) = (9*7**(0.5D0)*ETAH)/2048 
      Q2(33,46) = 0 
      Q2(33,47) = -(27*35**(0.5D0)*ETAH)/51200 
      Q2(33,48) = 0 
      Q2(33,49) = (15*35**(0.5D0)*XIH)/2048 
      Q2(33,50) = 0 
      Q2(33,51) = -(9*7**(0.5D0)*XIH)/2048 
      Q2(33,52) = 0 
      Q2(33,53) = 0 
      Q2(33,54) = 0 
      Q2(33,55) = 0 
      Q2(33,56) = 0 
      Q2(33,57) = -(9*7**(0.5D0)*XIH)/2048 
      Q2(33,58) = 0 
      Q2(33,59) = (27*35**(0.5D0)*XIH)/51200 
      Q2(33,60) = 0 
      Q2(33,61) = 0 
      Q2(33,62) = 0 
      Q2(33,63) = 0 
      Q2(33,64) = 0 

      Q2(34,1) = (11*15**(0.5D0)*MUH)/10240 
      Q2(34,2) = -(5**(0.5D0)*(85*VT + 34*ABS(ETAH) + 
     >   66*ABS(MUH) - 170*ABS(XIH)))/20480 
      Q2(34,3) = -(37*3**(0.5D0)*MUH)/1024 
      Q2(34,4) = (3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(45*VT + 
     >   18*ABS(ETAH) - 110*ABS(MUH) - 
     >   90*ABS(XIH)))/102400 
      Q2(34,5) = 0 
      Q2(34,6) = -(17*15**(0.5D0)*ETAH)/2560 
      Q2(34,7) = 0 
      Q2(34,8) = (27*35**(0.5D0)*ETAH)/12800 
      Q2(34,9) = -(33*3**(0.5D0)*MUH)/51200 
      Q2(34,10) = (51*VT)/20480 - (17*ABS(ETAH))/2048 + 
     >   (99*ABS(MUH))/51200 - (51*ABS(XIH))/10240 
      Q2(34,11) = (111*15**(0.5D0)*MUH)/25600 
      Q2(34,12) = (3*3**(0.5D0)*7**(0.5D0)*(30*ABS(ETAH) - 
     >   9*VT + 22*ABS(MUH) + 18*ABS(XIH)))/102400 
      Q2(34,13) = 0 
      Q2(34,14) = -(17*35**(0.5D0)*ETAH)/10240 
      Q2(34,15) = 0 
      Q2(34,16) = (63*15**(0.5D0)*ETAH)/51200 
      Q2(34,17) = 0 
      Q2(34,18) = -(17*15**(0.5D0)*XIH)/1024 
      Q2(34,19) = 0 
      Q2(34,20) = (27*35**(0.5D0)*XIH)/5120 
      Q2(34,21) = 0 
      Q2(34,22) = 0 
      Q2(34,23) = 0 
      Q2(34,24) = 0 
      Q2(34,25) = 0 
      Q2(34,26) = (51*3**(0.5D0)*XIH)/5120 
      Q2(34,27) = 0 
      Q2(34,28) = -(81*7**(0.5D0)*XIH)/25600 
      Q2(34,29) = 0 
      Q2(34,30) = 0 
      Q2(34,31) = 0 
      Q2(34,32) = 0 
      Q2(34,33) = -(33*3**(0.5D0)*MUH)/10240 
      Q2(34,34) = (51*VT)/4096 + (51*ABS(ETAH))/10240 + 
     >   (99*ABS(MUH))/10240 + (85*ABS(XIH))/2048 
      Q2(34,35) = (111*15**(0.5D0)*MUH)/5120 
      Q2(34,36) = -(3*3**(0.5D0)*7**(0.5D0)*(45*VT + 
     >   18*ABS(ETAH) - 110*ABS(MUH) + 
     >   150*ABS(XIH)))/102400 
      Q2(34,37) = 0 
      Q2(34,38) = (51*3**(0.5D0)*ETAH)/2560 
      Q2(34,39) = 0 
      Q2(34,40) = -(81*7**(0.5D0)*ETAH)/12800 
      Q2(34,41) = (99*15**(0.5D0)*MUH)/256000 
      Q2(34,42) = -(3*5**(0.5D0)*(255*VT - 850*ABS(ETAH) + 
     >   198*ABS(MUH) + 850*ABS(XIH)))/512000 
      Q2(34,43) = -(333*3**(0.5D0)*MUH)/25600 
      Q2(34,44) = (9*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(9*VT - 
     >   30*ABS(ETAH) - 22*ABS(MUH) + 
     >   30*ABS(XIH)))/512000 
      Q2(34,45) = 0 
      Q2(34,46) = (51*7**(0.5D0)*ETAH)/10240 
      Q2(34,47) = 0 
      Q2(34,48) = -(189*3**(0.5D0)*ETAH)/51200 
      Q2(34,49) = 0 
      Q2(34,50) = (17*35**(0.5D0)*XIH)/2048 
      Q2(34,51) = 0 
      Q2(34,52) = -(63*15**(0.5D0)*XIH)/10240 
      Q2(34,53) = 0 
      Q2(34,54) = 0 
      Q2(34,55) = 0 
      Q2(34,56) = 0 
      Q2(34,57) = 0 
      Q2(34,58) = -(51*7**(0.5D0)*XIH)/10240 
      Q2(34,59) = 0 
      Q2(34,60) = (189*3**(0.5D0)*XIH)/51200 
      Q2(34,61) = 0 
      Q2(34,62) = 0 
      Q2(34,63) = 0 
      Q2(34,64) = 0 

      Q2(35,1) = (9*VT)/4096 + (9*ABS(ETAH))/10240 - 
     >   (9*ABS(MUH))/2048 - (9*ABS(XIH))/2048 
      Q2(35,2) = (9*3**(0.5D0)*MUH)/1024 
      Q2(35,3) = -(9*5**(0.5D0)*(15*VT + 6*ABS(ETAH) + 
     >   50*ABS(MUH) - 30*ABS(XIH)))/102400 
      Q2(35,4) = -(9*7**(0.5D0)*MUH)/2048 
      Q2(35,5) = (9*3**(0.5D0)*ETAH)/2560 
      Q2(35,6) = 0 
      Q2(35,7) = -(27*15**(0.5D0)*ETAH)/12800 
      Q2(35,8) = 0 
      Q2(35,9) = (9*5**(0.5D0)*(10*ABS(ETAH) - 3*VT + 
     >   6*ABS(MUH) + 6*ABS(XIH)))/102400 
      Q2(35,10) = -(27*15**(0.5D0)*MUH)/25600 
      Q2(35,11) = (81*VT)/102400 - (27*ABS(ETAH))/10240 + 
     >   (27*ABS(MUH))/10240 - (81*ABS(XIH))/51200 
      Q2(35,12) = (27*35**(0.5D0)*MUH)/51200 
      Q2(35,13) = (9*7**(0.5D0)*ETAH)/10240 
      Q2(35,14) = 0 
      Q2(35,15) = -(27*35**(0.5D0)*ETAH)/51200 
      Q2(35,16) = 0 
      Q2(35,17) = (9*3**(0.5D0)*XIH)/1024 
      Q2(35,18) = 0 
      Q2(35,19) = -(27*15**(0.5D0)*XIH)/5120 
      Q2(35,20) = 0 
      Q2(35,21) = 0 
      Q2(35,22) = 0 
      Q2(35,23) = 0 
      Q2(35,24) = 0 
      Q2(35,25) = -(27*15**(0.5D0)*XIH)/25600 
      Q2(35,26) = 0 
      Q2(35,27) = (81*3**(0.5D0)*XIH)/25600 
      Q2(35,28) = 0 
      Q2(35,29) = 0 
      Q2(35,30) = 0 
      Q2(35,31) = 0 
      Q2(35,32) = 0 
      Q2(35,33) = -(9*5**(0.5D0)*(15*VT + 6*ABS(ETAH) - 
     >   30*ABS(MUH) + 50*ABS(XIH)))/102400 
      Q2(35,34) = -(27*15**(0.5D0)*MUH)/5120 
      Q2(35,35) = (81*VT)/20480 + (81*ABS(ETAH))/51200 + 
     >   (27*ABS(MUH))/2048 + (27*ABS(XIH))/2048 
      Q2(35,36) = (27*35**(0.5D0)*MUH)/10240 
      Q2(35,37) = -(27*15**(0.5D0)*ETAH)/12800 
      Q2(35,38) = 0 
      Q2(35,39) = (81*3**(0.5D0)*ETAH)/12800 
      Q2(35,40) = 0 
      Q2(35,41) = (81*VT)/102400 - (27*ABS(ETAH))/10240 - 
     >   (81*ABS(MUH))/51200 + (27*ABS(XIH))/10240 
      Q2(35,42) = (81*3**(0.5D0)*MUH)/25600 
      Q2(35,43) = -(81*5**(0.5D0)*(3*VT - 10*ABS(ETAH) + 
     >   10*ABS(MUH) + 10*ABS(XIH)))/512000 
      Q2(35,44) = -(81*7**(0.5D0)*MUH)/51200 
      Q2(35,45) = -(27*35**(0.5D0)*ETAH)/51200 
      Q2(35,46) = 0 
      Q2(35,47) = (81*7**(0.5D0)*ETAH)/51200 
      Q2(35,48) = 0 
      Q2(35,49) = -(9*7**(0.5D0)*XIH)/2048 
      Q2(35,50) = 0 
      Q2(35,51) = (27*35**(0.5D0)*XIH)/10240 
      Q2(35,52) = 0 
      Q2(35,53) = 0 
      Q2(35,54) = 0 
      Q2(35,55) = 0 
      Q2(35,56) = 0 
      Q2(35,57) = (27*35**(0.5D0)*XIH)/51200 
      Q2(35,58) = 0 
      Q2(35,59) = -(81*7**(0.5D0)*XIH)/51200 
      Q2(35,60) = 0 
      Q2(35,61) = 0 
      Q2(35,62) = 0 
      Q2(35,63) = 0 
      Q2(35,64) = 0 

      Q2(36,1) = (27*35**(0.5D0)*MUH)/71680 
      Q2(36,2) = (9*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(35*VT + 
     >   14*ABS(ETAH) - 30*ABS(MUH) - 
     >   70*ABS(XIH)))/716800 
      Q2(36,3) = (27*7**(0.5D0)*MUH)/1792 
      Q2(36,4) = -(27*5**(0.5D0)*(45*VT + 18*ABS(ETAH) + 
     >   70*ABS(MUH) - 90*ABS(XIH)))/716800 
      Q2(36,5) = 0 
      Q2(36,6) = (27*35**(0.5D0)*ETAH)/12800 
      Q2(36,7) = 0 
      Q2(36,8) = -(243*15**(0.5D0)*ETAH)/89600 
      Q2(36,9) = -(81*7**(0.5D0)*MUH)/358400 
      Q2(36,10) = (9*3**(0.5D0)*7**(0.5D0)*(70*ABS(ETAH) - 
     >   21*VT + 18*ABS(MUH) + 42*ABS(XIH)))/716800 
      Q2(36,11) = -(81*35**(0.5D0)*MUH)/44800 
      Q2(36,12) = (729*VT)/716800 - (243*ABS(ETAH))/71680 + 
     >   (81*ABS(MUH))/51200 - (729*ABS(XIH))/358400 
      Q2(36,13) = 0 
      Q2(36,14) = (63*15**(0.5D0)*ETAH)/51200 
      Q2(36,15) = 0 
      Q2(36,16) = -(243*35**(0.5D0)*ETAH)/358400 
      Q2(36,17) = 0 
      Q2(36,18) = (27*35**(0.5D0)*XIH)/5120 
      Q2(36,19) = 0 
      Q2(36,20) = -(243*15**(0.5D0)*XIH)/35840 
      Q2(36,21) = 0 
      Q2(36,22) = 0 
      Q2(36,23) = 0 
      Q2(36,24) = 0 
      Q2(36,25) = 0 
      Q2(36,26) = -(81*7**(0.5D0)*XIH)/25600 
      Q2(36,27) = 0 
      Q2(36,28) = (729*3**(0.5D0)*XIH)/179200 
      Q2(36,29) = 0 
      Q2(36,30) = 0 
      Q2(36,31) = 0 
      Q2(36,32) = 0 
      Q2(36,33) = -(81*7**(0.5D0)*MUH)/71680 
      Q2(36,34) = -(9*3**(0.5D0)*7**(0.5D0)*(105*VT + 
     >   42*ABS(ETAH) - 90*ABS(MUH) + 
     >   350*ABS(XIH)))/716800 
      Q2(36,35) = -(81*35**(0.5D0)*MUH)/8960 
      Q2(36,36) = (729*VT)/143360 + (729*ABS(ETAH))/358400 + 
     >   (81*ABS(MUH))/10240 + (243*ABS(XIH))/14336 
      Q2(36,37) = 0 
      Q2(36,38) = -(81*7**(0.5D0)*ETAH)/12800 
      Q2(36,39) = 0 
      Q2(36,40) = (729*3**(0.5D0)*ETAH)/89600 
      Q2(36,41) = (243*35**(0.5D0)*MUH)/1792000 
      Q2(36,42) = (27*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(21*VT - 
     >   70*ABS(ETAH) - 18*ABS(MUH) + 
     >   70*ABS(XIH)))/3584000 
      Q2(36,43) = (243*7**(0.5D0)*MUH)/44800 
      Q2(36,44) = -(243*5**(0.5D0)*(9*VT - 30*ABS(ETAH) + 
     >   14*ABS(MUH) + 30*ABS(XIH)))/3584000 
      Q2(36,45) = 0 
      Q2(36,46) = -(189*3**(0.5D0)*ETAH)/51200 
      Q2(36,47) = 0 
      Q2(36,48) = (729*7**(0.5D0)*ETAH)/358400 
      Q2(36,49) = 0 
      Q2(36,50) = -(63*15**(0.5D0)*XIH)/10240 
      Q2(36,51) = 0 
      Q2(36,52) = (243*35**(0.5D0)*XIH)/71680 
      Q2(36,53) = 0 
      Q2(36,54) = 0 
      Q2(36,55) = 0 
      Q2(36,56) = 0 
      Q2(36,57) = 0 
      Q2(36,58) = (189*3**(0.5D0)*XIH)/51200 
      Q2(36,59) = 0 
      Q2(36,60) = -(729*7**(0.5D0)*XIH)/358400 
      Q2(36,61) = 0 
      Q2(36,62) = 0 
      Q2(36,63) = 0 
      Q2(36,64) = 0 

      Q2(37,1) = (11*15**(0.5D0)*ETAH)/10240 
      Q2(37,2) = 0 
      Q2(37,3) = -(33*3**(0.5D0)*ETAH)/51200 
      Q2(37,4) = 0 
      Q2(37,5) = -(5**(0.5D0)*(85*VT + 66*ABS(ETAH) + 
     >   34*ABS(MUH) - 170*ABS(XIH)))/20480 
      Q2(37,6) = -(17*15**(0.5D0)*MUH)/2560 
      Q2(37,7) = (51*VT)/20480 + (99*ABS(ETAH))/51200 - 
     >   (17*ABS(MUH))/2048 - (51*ABS(XIH))/10240 
      Q2(37,8) = -(17*35**(0.5D0)*MUH)/10240 
      Q2(37,9) = -(37*3**(0.5D0)*ETAH)/1024 
      Q2(37,10) = 0 
      Q2(37,11) = (111*15**(0.5D0)*ETAH)/25600 
      Q2(37,12) = 0 
      Q2(37,13) = (3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(45*VT - 
     >   110*ABS(ETAH) + 18*ABS(MUH) - 
     >   90*ABS(XIH)))/102400 
      Q2(37,14) = (27*35**(0.5D0)*MUH)/12800 
      Q2(37,15) = (3*3**(0.5D0)*7**(0.5D0)*(22*ABS(ETAH) - 
     >   9*VT + 30*ABS(MUH) + 18*ABS(XIH)))/102400 
      Q2(37,16) = (63*15**(0.5D0)*MUH)/51200 
      Q2(37,17) = 0 
      Q2(37,18) = 0 
      Q2(37,19) = 0 
      Q2(37,20) = 0 
      Q2(37,21) = -(17*15**(0.5D0)*XIH)/1024 
      Q2(37,22) = 0 
      Q2(37,23) = (51*3**(0.5D0)*XIH)/5120 
      Q2(37,24) = 0 
      Q2(37,25) = 0 
      Q2(37,26) = 0 
      Q2(37,27) = 0 
      Q2(37,28) = 0 
      Q2(37,29) = (27*35**(0.5D0)*XIH)/5120 
      Q2(37,30) = 0 
      Q2(37,31) = -(81*7**(0.5D0)*XIH)/25600 
      Q2(37,32) = 0 
      Q2(37,33) = -(33*3**(0.5D0)*ETAH)/10240 
      Q2(37,34) = 0 
      Q2(37,35) = (99*15**(0.5D0)*ETAH)/256000 
      Q2(37,36) = 0 
      Q2(37,37) = (51*VT)/4096 + (99*ABS(ETAH))/10240 + 
     >   (51*ABS(MUH))/10240 + (85*ABS(XIH))/2048 
      Q2(37,38) = (51*3**(0.5D0)*MUH)/2560 
      Q2(37,39) = -(3*5**(0.5D0)*(255*VT + 198*ABS(ETAH) - 
     >   850*ABS(MUH) + 850*ABS(XIH)))/512000 
      Q2(37,40) = (51*7**(0.5D0)*MUH)/10240 
      Q2(37,41) = (111*15**(0.5D0)*ETAH)/5120 
      Q2(37,42) = 0 
      Q2(37,43) = -(333*3**(0.5D0)*ETAH)/25600 
      Q2(37,44) = 0 
      Q2(37,45) = -(3*3**(0.5D0)*7**(0.5D0)*(45*VT - 
     >   110*ABS(ETAH) + 18*ABS(MUH) + 
     >   150*ABS(XIH)))/102400 
      Q2(37,46) = -(81*7**(0.5D0)*MUH)/12800 
      Q2(37,47) = (9*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(9*VT - 
     >   22*ABS(ETAH) - 30*ABS(MUH) + 
     >   30*ABS(XIH)))/512000 
      Q2(37,48) = -(189*3**(0.5D0)*MUH)/51200 
      Q2(37,49) = 0 
      Q2(37,50) = 0 
      Q2(37,51) = 0 
      Q2(37,52) = 0 
      Q2(37,53) = (17*35**(0.5D0)*XIH)/2048 
      Q2(37,54) = 0 
      Q2(37,55) = -(51*7**(0.5D0)*XIH)/10240 
      Q2(37,56) = 0 
      Q2(37,57) = 0 
      Q2(37,58) = 0 
      Q2(37,59) = 0 
      Q2(37,60) = 0 
      Q2(37,61) = -(63*15**(0.5D0)*XIH)/10240 
      Q2(37,62) = 0 
      Q2(37,63) = (189*3**(0.5D0)*XIH)/51200 
      Q2(37,64) = 0 

      Q2(38,1) = 0 
      Q2(38,2) = (187*15**(0.5D0)*ETAH)/153600 
      Q2(38,3) = 0 
      Q2(38,4) = -(99*35**(0.5D0)*ETAH)/256000 
      Q2(38,5) = (187*15**(0.5D0)*MUH)/153600 
      Q2(38,6) = -(17*5**(0.5D0)*(85*VT + 66*ABS(ETAH) + 
     >   66*ABS(MUH) - 170*ABS(XIH)))/307200 
      Q2(38,7) = -(629*3**(0.5D0)*MUH)/15360 
      Q2(38,8) = (3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(765*VT + 
     >   594*ABS(ETAH) - 1870*ABS(MUH) - 
     >   1530*ABS(XIH)))/1536000 
      Q2(38,9) = 0 
      Q2(38,10) = -(629*3**(0.5D0)*ETAH)/15360 
      Q2(38,11) = 0 
      Q2(38,12) = (333*7**(0.5D0)*ETAH)/25600 
      Q2(38,13) = -(99*35**(0.5D0)*MUH)/256000 
      Q2(38,14) = (3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(765*VT - 
     >   1870*ABS(ETAH) + 594*ABS(MUH) - 
     >   1530*ABS(XIH)))/1536000 
      Q2(38,15) = (333*7**(0.5D0)*MUH)/25600 
      Q2(38,16) = (63*5**(0.5D0)*(22*ABS(ETAH) - 9*VT + 
     >   22*ABS(MUH) + 18*ABS(XIH)))/512000 
      Q2(38,17) = 0 
      Q2(38,18) = 0 
      Q2(38,19) = 0 
      Q2(38,20) = 0 
      Q2(38,21) = 0 
      Q2(38,22) = -(289*15**(0.5D0)*XIH)/15360 
      Q2(38,23) = 0 
      Q2(38,24) = (153*35**(0.5D0)*XIH)/25600 
      Q2(38,25) = 0 
      Q2(38,26) = 0 
      Q2(38,27) = 0 
      Q2(38,28) = 0 
      Q2(38,29) = 0 
      Q2(38,30) = (153*35**(0.5D0)*XIH)/25600 
      Q2(38,31) = 0 
      Q2(38,32) = -(567*15**(0.5D0)*XIH)/128000 
      Q2(38,33) = 0 
      Q2(38,34) = -(187*3**(0.5D0)*ETAH)/51200 
      Q2(38,35) = 0 
      Q2(38,36) = (297*7**(0.5D0)*ETAH)/256000 
      Q2(38,37) = -(187*3**(0.5D0)*MUH)/51200 
      Q2(38,38) = (289*VT)/20480 + (561*ABS(ETAH))/51200 + 
     >   (561*ABS(MUH))/51200 + (289*ABS(XIH))/6144 
      Q2(38,39) = (629*15**(0.5D0)*MUH)/25600 
      Q2(38,40) = -(3**(0.5D0)*7**(0.5D0)*(765*VT + 
     >   594*ABS(ETAH) - 1870*ABS(MUH) + 
     >   2550*ABS(XIH)))/512000 
      Q2(38,41) = 0 
      Q2(38,42) = (629*15**(0.5D0)*ETAH)/25600 
      Q2(38,43) = 0 
      Q2(38,44) = -(999*35**(0.5D0)*ETAH)/128000 
      Q2(38,45) = (297*7**(0.5D0)*MUH)/256000 
      Q2(38,46) = -(3**(0.5D0)*7**(0.5D0)*(765*VT - 
     >   1870*ABS(ETAH) + 594*ABS(MUH) + 
     >   2550*ABS(XIH)))/512000 
      Q2(38,47) = -(999*35**(0.5D0)*MUH)/128000 
      Q2(38,48) = (1701*VT)/512000 - (2079*ABS(ETAH))/256000 - 
     >   (2079*ABS(MUH))/256000 + (567*ABS(XIH))/51200 
      Q2(38,49) = 0 
      Q2(38,50) = 0 
      Q2(38,51) = 0 
      Q2(38,52) = 0 
      Q2(38,53) = 0 
      Q2(38,54) = (289*35**(0.5D0)*XIH)/30720 
      Q2(38,55) = 0 
      Q2(38,56) = -(357*15**(0.5D0)*XIH)/51200 
      Q2(38,57) = 0 
      Q2(38,58) = 0 
      Q2(38,59) = 0 
      Q2(38,60) = 0 
      Q2(38,61) = 0 
      Q2(38,62) = -(357*15**(0.5D0)*XIH)/51200 
      Q2(38,63) = 0 
      Q2(38,64) = (567*35**(0.5D0)*XIH)/256000 

      Q2(39,1) = -(33*3**(0.5D0)*ETAH)/51200 
      Q2(39,2) = 0 
      Q2(39,3) = (99*15**(0.5D0)*ETAH)/256000 
      Q2(39,4) = 0 
      Q2(39,5) = (51*VT)/20480 + (99*ABS(ETAH))/51200 - 
     >   (51*ABS(MUH))/10240 - (51*ABS(XIH))/10240 
      Q2(39,6) = (51*3**(0.5D0)*MUH)/5120 
      Q2(39,7) = -(3*5**(0.5D0)*(255*VT + 198*ABS(ETAH) + 
     >   850*ABS(MUH) - 510*ABS(XIH)))/512000 
      Q2(39,8) = -(51*7**(0.5D0)*MUH)/10240 
      Q2(39,9) = (111*15**(0.5D0)*ETAH)/25600 
      Q2(39,10) = 0 
      Q2(39,11) = -(333*3**(0.5D0)*ETAH)/25600 
      Q2(39,12) = 0 
      Q2(39,13) = (3*3**(0.5D0)*7**(0.5D0)*(22*ABS(ETAH) - 
     >   9*VT + 18*ABS(MUH) + 18*ABS(XIH)))/102400 
      Q2(39,14) = -(81*7**(0.5D0)*MUH)/25600 
      Q2(39,15) = (9*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(9*VT - 
     >   22*ABS(ETAH) + 30*ABS(MUH) - 
     >   18*ABS(XIH)))/512000 
      Q2(39,16) = (189*3**(0.5D0)*MUH)/51200 
      Q2(39,17) = 0 
      Q2(39,18) = 0 
      Q2(39,19) = 0 
      Q2(39,20) = 0 
      Q2(39,21) = (51*3**(0.5D0)*XIH)/5120 
      Q2(39,22) = 0 
      Q2(39,23) = -(153*15**(0.5D0)*XIH)/25600 
      Q2(39,24) = 0 
      Q2(39,25) = 0 
      Q2(39,26) = 0 
      Q2(39,27) = 0 
      Q2(39,28) = 0 
      Q2(39,29) = -(81*7**(0.5D0)*XIH)/25600 
      Q2(39,30) = 0 
      Q2(39,31) = (243*35**(0.5D0)*XIH)/128000 
      Q2(39,32) = 0 
      Q2(39,33) = (99*15**(0.5D0)*ETAH)/256000 
      Q2(39,34) = 0 
      Q2(39,35) = -(297*3**(0.5D0)*ETAH)/256000 
      Q2(39,36) = 0 
      Q2(39,37) = -(3*5**(0.5D0)*(255*VT + 198*ABS(ETAH) - 
     >   510*ABS(MUH) + 850*ABS(XIH)))/512000 
      Q2(39,38) = -(153*15**(0.5D0)*MUH)/25600 
      Q2(39,39) = (459*VT)/102400 + (891*ABS(ETAH))/256000 + 
     >   (153*ABS(MUH))/10240 + (153*ABS(XIH))/10240 
      Q2(39,40) = (153*35**(0.5D0)*MUH)/51200 
      Q2(39,41) = -(333*3**(0.5D0)*ETAH)/25600 
      Q2(39,42) = 0 
      Q2(39,43) = (999*15**(0.5D0)*ETAH)/128000 
      Q2(39,44) = 0 
      Q2(39,45) = (9*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(9*VT - 
     >   22*ABS(ETAH) - 18*ABS(MUH) + 
     >   30*ABS(XIH)))/512000 
      Q2(39,46) = (243*35**(0.5D0)*MUH)/128000 
      Q2(39,47) = -(27*3**(0.5D0)*7**(0.5D0)*(9*VT - 
     >   22*ABS(ETAH) + 30*ABS(MUH) + 
     >   30*ABS(XIH)))/512000 
      Q2(39,48) = -(567*15**(0.5D0)*MUH)/256000 
      Q2(39,49) = 0 
      Q2(39,50) = 0 
      Q2(39,51) = 0 
      Q2(39,52) = 0 
      Q2(39,53) = -(51*7**(0.5D0)*XIH)/10240 
      Q2(39,54) = 0 
      Q2(39,55) = (153*35**(0.5D0)*XIH)/51200 
      Q2(39,56) = 0 
      Q2(39,57) = 0 
      Q2(39,58) = 0 
      Q2(39,59) = 0 
      Q2(39,60) = 0 
      Q2(39,61) = (189*3**(0.5D0)*XIH)/51200 
      Q2(39,62) = 0 
      Q2(39,63) = -(567*15**(0.5D0)*XIH)/256000 
      Q2(39,64) = 0 

      Q2(40,1) = 0 
      Q2(40,2) = -(99*35**(0.5D0)*ETAH)/256000 
      Q2(40,3) = 0 
      Q2(40,4) = (891*15**(0.5D0)*ETAH)/1792000 
      Q2(40,5) = (153*35**(0.5D0)*MUH)/358400 
      Q2(40,6) = (3*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(595*VT + 
     >   462*ABS(ETAH) - 510*ABS(MUH) - 
     >   1190*ABS(XIH)))/3584000 
      Q2(40,7) = (153*7**(0.5D0)*MUH)/8960 
      Q2(40,8) = -(9*5**(0.5D0)*(765*VT + 594*ABS(ETAH) + 
     >   1190*ABS(MUH) - 1530*ABS(XIH)))/3584000 
      Q2(40,9) = 0 
      Q2(40,10) = (333*7**(0.5D0)*ETAH)/25600 
      Q2(40,11) = 0 
      Q2(40,12) = -(2997*3**(0.5D0)*ETAH)/179200 
      Q2(40,13) = -(81*15**(0.5D0)*MUH)/256000 
      Q2(40,14) = (9*5**(0.5D0)*(154*ABS(ETAH) - 63*VT + 
     >   54*ABS(MUH) + 126*ABS(XIH)))/512000 
      Q2(40,15) = -(81*3**(0.5D0)*MUH)/6400 
      Q2(40,16) = (81*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(9*VT - 
     >   22*ABS(ETAH) + 14*ABS(MUH) - 
     >   18*ABS(XIH)))/3584000 
      Q2(40,17) = 0 
      Q2(40,18) = 0 
      Q2(40,19) = 0 
      Q2(40,20) = 0 
      Q2(40,21) = 0 
      Q2(40,22) = (153*35**(0.5D0)*XIH)/25600 
      Q2(40,23) = 0 
      Q2(40,24) = -(1377*15**(0.5D0)*XIH)/179200 
      Q2(40,25) = 0 
      Q2(40,26) = 0 
      Q2(40,27) = 0 
      Q2(40,28) = 0 
      Q2(40,29) = 0 
      Q2(40,30) = -(567*15**(0.5D0)*XIH)/128000 
      Q2(40,31) = 0 
      Q2(40,32) = (2187*35**(0.5D0)*XIH)/896000 
      Q2(40,33) = 0 
      Q2(40,34) = (297*7**(0.5D0)*ETAH)/256000 
      Q2(40,35) = 0 
      Q2(40,36) = -(2673*3**(0.5D0)*ETAH)/1792000 
      Q2(40,37) = -(459*7**(0.5D0)*MUH)/358400 
      Q2(40,38) = -(3*3**(0.5D0)*7**(0.5D0)*(1785*VT + 
     >   1386*ABS(ETAH) - 1530*ABS(MUH) + 
     >   5950*ABS(XIH)))/3584000 
      Q2(40,39) = -(459*35**(0.5D0)*MUH)/44800 
      Q2(40,40) = (4131*VT)/716800 + (8019*ABS(ETAH))/1792000 + 
     >   (459*ABS(MUH))/51200 + (1377*ABS(XIH))/71680 
      Q2(40,41) = 0 
      Q2(40,42) = -(999*35**(0.5D0)*ETAH)/128000 
      Q2(40,43) = 0 
      Q2(40,44) = (8991*15**(0.5D0)*ETAH)/896000 
      Q2(40,45) = (243*3**(0.5D0)*MUH)/256000 
      Q2(40,46) = (1701*VT)/512000 - (2079*ABS(ETAH))/256000 - 
     >   (729*ABS(MUH))/256000 + (567*ABS(XIH))/51200 
      Q2(40,47) = (243*15**(0.5D0)*MUH)/32000 
      Q2(40,48) = -(243*3**(0.5D0)*7**(0.5D0)*(9*VT - 
     >   22*ABS(ETAH) + 14*ABS(MUH) + 
     >   30*ABS(XIH)))/3584000 
      Q2(40,49) = 0 
      Q2(40,50) = 0 
      Q2(40,51) = 0 
      Q2(40,52) = 0 
      Q2(40,53) = 0 
      Q2(40,54) = -(357*15**(0.5D0)*XIH)/51200 
      Q2(40,55) = 0 
      Q2(40,56) = (1377*35**(0.5D0)*XIH)/358400 
      Q2(40,57) = 0 
      Q2(40,58) = 0 
      Q2(40,59) = 0 
      Q2(40,60) = 0 
      Q2(40,61) = 0 
      Q2(40,62) = (567*35**(0.5D0)*XIH)/256000 
      Q2(40,63) = 0 
      Q2(40,64) = -(729*15**(0.5D0)*XIH)/256000 

      Q2(41,1) = (9*VT)/4096 - (9*ABS(ETAH))/2048 + 
     >   (9*ABS(MUH))/10240 - (9*ABS(XIH))/2048 
      Q2(41,2) = (9*3**(0.5D0)*MUH)/2560 
      Q2(41,3) = (9*5**(0.5D0)*(6*ABS(ETAH) - 3*VT + 
     >   10*ABS(MUH) + 6*ABS(XIH)))/102400 
      Q2(41,4) = (9*7**(0.5D0)*MUH)/10240 
      Q2(41,5) = (9*3**(0.5D0)*ETAH)/1024 
      Q2(41,6) = 0 
      Q2(41,7) = -(27*15**(0.5D0)*ETAH)/25600 
      Q2(41,8) = 0 
      Q2(41,9) = -(9*5**(0.5D0)*(15*VT + 50*ABS(ETAH) + 
     >   6*ABS(MUH) - 30*ABS(XIH)))/102400 
      Q2(41,10) = -(27*15**(0.5D0)*MUH)/12800 
      Q2(41,11) = (81*VT)/102400 + (27*ABS(ETAH))/10240 - 
     >   (27*ABS(MUH))/10240 - (81*ABS(XIH))/51200 
      Q2(41,12) = -(27*35**(0.5D0)*MUH)/51200 
      Q2(41,13) = -(9*7**(0.5D0)*ETAH)/2048 
      Q2(41,14) = 0 
      Q2(41,15) = (27*35**(0.5D0)*ETAH)/51200 
      Q2(41,16) = 0 
      Q2(41,17) = (9*3**(0.5D0)*XIH)/1024 
      Q2(41,18) = 0 
      Q2(41,19) = -(27*15**(0.5D0)*XIH)/25600 
      Q2(41,20) = 0 
      Q2(41,21) = 0 
      Q2(41,22) = 0 
      Q2(41,23) = 0 
      Q2(41,24) = 0 
      Q2(41,25) = -(27*15**(0.5D0)*XIH)/5120 
      Q2(41,26) = 0 
      Q2(41,27) = (81*3**(0.5D0)*XIH)/25600 
      Q2(41,28) = 0 
      Q2(41,29) = 0 
      Q2(41,30) = 0 
      Q2(41,31) = 0 
      Q2(41,32) = 0 
      Q2(41,33) = -(9*5**(0.5D0)*(15*VT - 30*ABS(ETAH) + 
     >   6*ABS(MUH) + 50*ABS(XIH)))/102400 
      Q2(41,34) = -(27*15**(0.5D0)*MUH)/12800 
      Q2(41,35) = (81*VT)/102400 - (81*ABS(ETAH))/51200 - 
     >   (27*ABS(MUH))/10240 + (27*ABS(XIH))/10240 
      Q2(41,36) = -(27*35**(0.5D0)*MUH)/51200 
      Q2(41,37) = -(27*15**(0.5D0)*ETAH)/5120 
      Q2(41,38) = 0 
      Q2(41,39) = (81*3**(0.5D0)*ETAH)/25600 
      Q2(41,40) = 0 
      Q2(41,41) = (81*VT)/20480 + (27*ABS(ETAH))/2048 + 
     >   (81*ABS(MUH))/51200 + (27*ABS(XIH))/2048 
      Q2(41,42) = (81*3**(0.5D0)*MUH)/12800 
      Q2(41,43) = -(81*5**(0.5D0)*(3*VT + 10*ABS(ETAH) - 
     >   10*ABS(MUH) + 10*ABS(XIH)))/512000 
      Q2(41,44) = (81*7**(0.5D0)*MUH)/51200 
      Q2(41,45) = (27*35**(0.5D0)*ETAH)/10240 
      Q2(41,46) = 0 
      Q2(41,47) = -(81*7**(0.5D0)*ETAH)/51200 
      Q2(41,48) = 0 
      Q2(41,49) = -(9*7**(0.5D0)*XIH)/2048 
      Q2(41,50) = 0 
      Q2(41,51) = (27*35**(0.5D0)*XIH)/51200 
      Q2(41,52) = 0 
      Q2(41,53) = 0 
      Q2(41,54) = 0 
      Q2(41,55) = 0 
      Q2(41,56) = 0 
      Q2(41,57) = (27*35**(0.5D0)*XIH)/10240 
      Q2(41,58) = 0 
      Q2(41,59) = -(81*7**(0.5D0)*XIH)/51200 
      Q2(41,60) = 0 
      Q2(41,61) = 0 
      Q2(41,62) = 0 
      Q2(41,63) = 0 
      Q2(41,64) = 0 

      Q2(42,1) = -(33*3**(0.5D0)*MUH)/51200 
      Q2(42,2) = (51*VT)/20480 - (51*ABS(ETAH))/10240 + 
     >   (99*ABS(MUH))/51200 - (51*ABS(XIH))/10240 
      Q2(42,3) = (111*15**(0.5D0)*MUH)/25600 
      Q2(42,4) = (3*3**(0.5D0)*7**(0.5D0)*(18*ABS(ETAH) - 
     >   9*VT + 22*ABS(MUH) + 18*ABS(XIH)))/102400 
      Q2(42,5) = 0 
      Q2(42,6) = (51*3**(0.5D0)*ETAH)/5120 
      Q2(42,7) = 0 
      Q2(42,8) = -(81*7**(0.5D0)*ETAH)/25600 
      Q2(42,9) = (99*15**(0.5D0)*MUH)/256000 
      Q2(42,10) = -(3*5**(0.5D0)*(255*VT + 850*ABS(ETAH) + 
     >   198*ABS(MUH) - 510*ABS(XIH)))/512000 
      Q2(42,11) = -(333*3**(0.5D0)*MUH)/25600 
      Q2(42,12) = (9*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(9*VT + 
     >   30*ABS(ETAH) - 22*ABS(MUH) - 
     >   18*ABS(XIH)))/512000 
      Q2(42,13) = 0 
      Q2(42,14) = -(51*7**(0.5D0)*ETAH)/10240 
      Q2(42,15) = 0 
      Q2(42,16) = (189*3**(0.5D0)*ETAH)/51200 
      Q2(42,17) = 0 
      Q2(42,18) = (51*3**(0.5D0)*XIH)/5120 
      Q2(42,19) = 0 
      Q2(42,20) = -(81*7**(0.5D0)*XIH)/25600 
      Q2(42,21) = 0 
      Q2(42,22) = 0 
      Q2(42,23) = 0 
      Q2(42,24) = 0 
      Q2(42,25) = 0 
      Q2(42,26) = -(153*15**(0.5D0)*XIH)/25600 
      Q2(42,27) = 0 
      Q2(42,28) = (243*35**(0.5D0)*XIH)/128000 
      Q2(42,29) = 0 
      Q2(42,30) = 0 
      Q2(42,31) = 0 
      Q2(42,32) = 0 
      Q2(42,33) = (99*15**(0.5D0)*MUH)/256000 
      Q2(42,34) = -(3*5**(0.5D0)*(255*VT - 510*ABS(ETAH) + 
     >   198*ABS(MUH) + 850*ABS(XIH)))/512000 
      Q2(42,35) = -(333*3**(0.5D0)*MUH)/25600 
      Q2(42,36) = (9*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(9*VT - 
     >   18*ABS(ETAH) - 22*ABS(MUH) + 
     >   30*ABS(XIH)))/512000 
      Q2(42,37) = 0 
      Q2(42,38) = -(153*15**(0.5D0)*ETAH)/25600 
      Q2(42,39) = 0 
      Q2(42,40) = (243*35**(0.5D0)*ETAH)/128000 
      Q2(42,41) = -(297*3**(0.5D0)*MUH)/256000 
      Q2(42,42) = (459*VT)/102400 + (153*ABS(ETAH))/10240 + 
     >   (891*ABS(MUH))/256000 + (153*ABS(XIH))/10240 
      Q2(42,43) = (999*15**(0.5D0)*MUH)/128000 
      Q2(42,44) = -(27*3**(0.5D0)*7**(0.5D0)*(9*VT + 
     >   30*ABS(ETAH) - 22*ABS(MUH) + 
     >   30*ABS(XIH)))/512000 
      Q2(42,45) = 0 
      Q2(42,46) = (153*35**(0.5D0)*ETAH)/51200 
      Q2(42,47) = 0 
      Q2(42,48) = -(567*15**(0.5D0)*ETAH)/256000 
      Q2(42,49) = 0 
      Q2(42,50) = -(51*7**(0.5D0)*XIH)/10240 
      Q2(42,51) = 0 
      Q2(42,52) = (189*3**(0.5D0)*XIH)/51200 
      Q2(42,53) = 0 
      Q2(42,54) = 0 
      Q2(42,55) = 0 
      Q2(42,56) = 0 
      Q2(42,57) = 0 
      Q2(42,58) = (153*35**(0.5D0)*XIH)/51200 
      Q2(42,59) = 0 
      Q2(42,60) = -(567*15**(0.5D0)*XIH)/256000 
      Q2(42,61) = 0 
      Q2(42,62) = 0 
      Q2(42,63) = 0 
      Q2(42,64) = 0 

      Q2(43,1) = (27*5**(0.5D0)*(2*ABS(ETAH) - VT + 
     >   2*ABS(MUH) + 2*ABS(XIH)))/102400 
      Q2(43,2) = -(27*15**(0.5D0)*MUH)/25600 
      Q2(43,3) = (81*VT)/102400 - (81*ABS(ETAH))/51200 + 
     >   (27*ABS(MUH))/10240 - (81*ABS(XIH))/51200 
      Q2(43,4) = (27*35**(0.5D0)*MUH)/51200 
      Q2(43,5) = -(27*15**(0.5D0)*ETAH)/25600 
      Q2(43,6) = 0 
      Q2(43,7) = (81*3**(0.5D0)*ETAH)/25600 
      Q2(43,8) = 0 
      Q2(43,9) = (81*VT)/102400 + (27*ABS(ETAH))/10240 - 
     >   (81*ABS(MUH))/51200 - (81*ABS(XIH))/51200 
      Q2(43,10) = (81*3**(0.5D0)*MUH)/25600 
      Q2(43,11) = -(81*5**(0.5D0)*(3*VT + 10*ABS(ETAH) + 
     >   10*ABS(MUH) - 6*ABS(XIH)))/512000 
      Q2(43,12) = -(81*7**(0.5D0)*MUH)/51200 
      Q2(43,13) = (27*35**(0.5D0)*ETAH)/51200 
      Q2(43,14) = 0 
      Q2(43,15) = -(81*7**(0.5D0)*ETAH)/51200 
      Q2(43,16) = 0 
      Q2(43,17) = -(27*15**(0.5D0)*XIH)/25600 
      Q2(43,18) = 0 
      Q2(43,19) = (81*3**(0.5D0)*XIH)/25600 
      Q2(43,20) = 0 
      Q2(43,21) = 0 
      Q2(43,22) = 0 
      Q2(43,23) = 0 
      Q2(43,24) = 0 
      Q2(43,25) = (81*3**(0.5D0)*XIH)/25600 
      Q2(43,26) = 0 
      Q2(43,27) = -(243*15**(0.5D0)*XIH)/128000 
      Q2(43,28) = 0 
      Q2(43,29) = 0 
      Q2(43,30) = 0 
      Q2(43,31) = 0 
      Q2(43,32) = 0 
      Q2(43,33) = (81*VT)/102400 - (81*ABS(ETAH))/51200 - 
     >   (81*ABS(MUH))/51200 + (27*ABS(XIH))/10240 
      Q2(43,34) = (81*3**(0.5D0)*MUH)/25600 
      Q2(43,35) = -(81*5**(0.5D0)*(3*VT - 6*ABS(ETAH) + 
     >   10*ABS(MUH) + 10*ABS(XIH)))/512000 
      Q2(43,36) = -(81*7**(0.5D0)*MUH)/51200 
      Q2(43,37) = (81*3**(0.5D0)*ETAH)/25600 
      Q2(43,38) = 0 
      Q2(43,39) = -(243*15**(0.5D0)*ETAH)/128000 
      Q2(43,40) = 0 
      Q2(43,41) = -(81*5**(0.5D0)*(3*VT + 10*ABS(ETAH) - 
     >   6*ABS(MUH) + 10*ABS(XIH)))/512000 
      Q2(43,42) = -(243*15**(0.5D0)*MUH)/128000 
      Q2(43,43) = (729*VT)/512000 + (243*ABS(ETAH))/51200 + 
     >   (243*ABS(MUH))/51200 + (243*ABS(XIH))/51200 
      Q2(43,44) = (243*35**(0.5D0)*MUH)/256000 
      Q2(43,45) = -(81*7**(0.5D0)*ETAH)/51200 
      Q2(43,46) = 0 
      Q2(43,47) = (243*35**(0.5D0)*ETAH)/256000 
      Q2(43,48) = 0 
      Q2(43,49) = (27*35**(0.5D0)*XIH)/51200 
      Q2(43,50) = 0 
      Q2(43,51) = -(81*7**(0.5D0)*XIH)/51200 
      Q2(43,52) = 0 
      Q2(43,53) = 0 
      Q2(43,54) = 0 
      Q2(43,55) = 0 
      Q2(43,56) = 0 
      Q2(43,57) = -(81*7**(0.5D0)*XIH)/51200 
      Q2(43,58) = 0 
      Q2(43,59) = (243*35**(0.5D0)*XIH)/256000 
      Q2(43,60) = 0 
      Q2(43,61) = 0 
      Q2(43,62) = 0 
      Q2(43,63) = 0 
      Q2(43,64) = 0 

      Q2(44,1) = -(81*7**(0.5D0)*MUH)/358400 
      Q2(44,2) = (27*3**(0.5D0)*7**(0.5D0)*(14*ABS(ETAH) - 
     >   7*VT + 6*ABS(MUH) + 14*ABS(XIH)))/716800 
      Q2(44,3) = -(81*35**(0.5D0)*MUH)/44800 
      Q2(44,4) = (729*VT)/716800 - (729*ABS(ETAH))/358400 + 
     >   (81*ABS(MUH))/51200 - (729*ABS(XIH))/358400 
      Q2(44,5) = 0 
      Q2(44,6) = -(81*7**(0.5D0)*ETAH)/25600 
      Q2(44,7) = 0 
      Q2(44,8) = (729*3**(0.5D0)*ETAH)/179200 
      Q2(44,9) = (243*35**(0.5D0)*MUH)/1792000 
      Q2(44,10) = (27*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(21*VT + 
     >   70*ABS(ETAH) - 18*ABS(MUH) - 
     >   42*ABS(XIH)))/3584000 
      Q2(44,11) = (243*7**(0.5D0)*MUH)/44800 
      Q2(44,12) = -(243*5**(0.5D0)*(9*VT + 30*ABS(ETAH) + 
     >   14*ABS(MUH) - 18*ABS(XIH)))/3584000 
      Q2(44,13) = 0 
      Q2(44,14) = (189*3**(0.5D0)*ETAH)/51200 
      Q2(44,15) = 0 
      Q2(44,16) = -(729*7**(0.5D0)*ETAH)/358400 
      Q2(44,17) = 0 
      Q2(44,18) = -(81*7**(0.5D0)*XIH)/25600 
      Q2(44,19) = 0 
      Q2(44,20) = (729*3**(0.5D0)*XIH)/179200 
      Q2(44,21) = 0 
      Q2(44,22) = 0 
      Q2(44,23) = 0 
      Q2(44,24) = 0 
      Q2(44,25) = 0 
      Q2(44,26) = (243*35**(0.5D0)*XIH)/128000 
      Q2(44,27) = 0 
      Q2(44,28) = -(2187*15**(0.5D0)*XIH)/896000 
      Q2(44,29) = 0 
      Q2(44,30) = 0 
      Q2(44,31) = 0 
      Q2(44,32) = 0 
      Q2(44,33) = (243*35**(0.5D0)*MUH)/1792000 
      Q2(44,34) = (27*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(21*VT - 
     >   42*ABS(ETAH) - 18*ABS(MUH) + 
     >   70*ABS(XIH)))/3584000 
      Q2(44,35) = (243*7**(0.5D0)*MUH)/44800 
      Q2(44,36) = -(243*5**(0.5D0)*(9*VT - 18*ABS(ETAH) + 
     >   14*ABS(MUH) + 30*ABS(XIH)))/3584000 
      Q2(44,37) = 0 
      Q2(44,38) = (243*35**(0.5D0)*ETAH)/128000 
      Q2(44,39) = 0 
      Q2(44,40) = -(2187*15**(0.5D0)*ETAH)/896000 
      Q2(44,41) = -(729*7**(0.5D0)*MUH)/1792000 
      Q2(44,42) = -(81*3**(0.5D0)*7**(0.5D0)*(21*VT + 
     >   70*ABS(ETAH) - 18*ABS(MUH) + 
     >   70*ABS(XIH)))/3584000 
      Q2(44,43) = -(729*35**(0.5D0)*MUH)/224000 
      Q2(44,44) = (6561*VT)/3584000 + (2187*ABS(ETAH))/358400 + 
     >   (729*ABS(MUH))/256000 + (2187*ABS(XIH))/358400 
      Q2(44,45) = 0 
      Q2(44,46) = -(567*15**(0.5D0)*ETAH)/256000 
      Q2(44,47) = 0 
      Q2(44,48) = (2187*35**(0.5D0)*ETAH)/1792000 
      Q2(44,49) = 0 
      Q2(44,50) = (189*3**(0.5D0)*XIH)/51200 
      Q2(44,51) = 0 
      Q2(44,52) = -(729*7**(0.5D0)*XIH)/358400 
      Q2(44,53) = 0 
      Q2(44,54) = 0 
      Q2(44,55) = 0 
      Q2(44,56) = 0 
      Q2(44,57) = 0 
      Q2(44,58) = -(567*15**(0.5D0)*XIH)/256000 
      Q2(44,59) = 0 
      Q2(44,60) = (2187*35**(0.5D0)*XIH)/1792000 
      Q2(44,61) = 0 
      Q2(44,62) = 0 
      Q2(44,63) = 0 
      Q2(44,64) = 0 

      Q2(45,1) = (27*35**(0.5D0)*ETAH)/71680 
      Q2(45,2) = 0 
      Q2(45,3) = -(81*7**(0.5D0)*ETAH)/358400 
      Q2(45,4) = 0 
      Q2(45,5) = (9*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(35*VT - 
     >   30*ABS(ETAH) + 14*ABS(MUH) - 
     >   70*ABS(XIH)))/716800 
      Q2(45,6) = (27*35**(0.5D0)*MUH)/12800 
      Q2(45,7) = (9*3**(0.5D0)*7**(0.5D0)*(18*ABS(ETAH) - 
     >   21*VT + 70*ABS(MUH) + 42*ABS(XIH)))/716800 
      Q2(45,8) = (63*15**(0.5D0)*MUH)/51200 
      Q2(45,9) = (27*7**(0.5D0)*ETAH)/1792 
      Q2(45,10) = 0 
      Q2(45,11) = -(81*35**(0.5D0)*ETAH)/44800 
      Q2(45,12) = 0 
      Q2(45,13) = -(27*5**(0.5D0)*(45*VT + 70*ABS(ETAH) + 
     >   18*ABS(MUH) - 90*ABS(XIH)))/716800 
      Q2(45,14) = -(243*15**(0.5D0)*MUH)/89600 
      Q2(45,15) = (729*VT)/716800 + (81*ABS(ETAH))/51200 - 
     >   (243*ABS(MUH))/71680 - (729*ABS(XIH))/358400 
      Q2(45,16) = -(243*35**(0.5D0)*MUH)/358400 
      Q2(45,17) = 0 
      Q2(45,18) = 0 
      Q2(45,19) = 0 
      Q2(45,20) = 0 
      Q2(45,21) = (27*35**(0.5D0)*XIH)/5120 
      Q2(45,22) = 0 
      Q2(45,23) = -(81*7**(0.5D0)*XIH)/25600 
      Q2(45,24) = 0 
      Q2(45,25) = 0 
      Q2(45,26) = 0 
      Q2(45,27) = 0 
      Q2(45,28) = 0 
      Q2(45,29) = -(243*15**(0.5D0)*XIH)/35840 
      Q2(45,30) = 0 
      Q2(45,31) = (729*3**(0.5D0)*XIH)/179200 
      Q2(45,32) = 0 
      Q2(45,33) = -(81*7**(0.5D0)*ETAH)/71680 
      Q2(45,34) = 0 
      Q2(45,35) = (243*35**(0.5D0)*ETAH)/1792000 
      Q2(45,36) = 0 
      Q2(45,37) = -(9*3**(0.5D0)*7**(0.5D0)*(105*VT - 
     >   90*ABS(ETAH) + 42*ABS(MUH) + 
     >   350*ABS(XIH)))/716800 
      Q2(45,38) = -(81*7**(0.5D0)*MUH)/12800 
      Q2(45,39) = (27*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(21*VT - 
     >   18*ABS(ETAH) - 70*ABS(MUH) + 
     >   70*ABS(XIH)))/3584000 
      Q2(45,40) = -(189*3**(0.5D0)*MUH)/51200 
      Q2(45,41) = -(81*35**(0.5D0)*ETAH)/8960 
      Q2(45,42) = 0 
      Q2(45,43) = (243*7**(0.5D0)*ETAH)/44800 
      Q2(45,44) = 0 
      Q2(45,45) = (729*VT)/143360 + (81*ABS(ETAH))/10240 + 
     >   (729*ABS(MUH))/358400 + (243*ABS(XIH))/14336 
      Q2(45,46) = (729*3**(0.5D0)*MUH)/89600 
      Q2(45,47) = -(243*5**(0.5D0)*(9*VT + 14*ABS(ETAH) - 
     >   30*ABS(MUH) + 30*ABS(XIH)))/3584000 
      Q2(45,48) = (729*7**(0.5D0)*MUH)/358400 
      Q2(45,49) = 0 
      Q2(45,50) = 0 
      Q2(45,51) = 0 
      Q2(45,52) = 0 
      Q2(45,53) = -(63*15**(0.5D0)*XIH)/10240 
      Q2(45,54) = 0 
      Q2(45,55) = (189*3**(0.5D0)*XIH)/51200 
      Q2(45,56) = 0 
      Q2(45,57) = 0 
      Q2(45,58) = 0 
      Q2(45,59) = 0 
      Q2(45,60) = 0 
      Q2(45,61) = (243*35**(0.5D0)*XIH)/71680 
      Q2(45,62) = 0 
      Q2(45,63) = -(729*7**(0.5D0)*XIH)/358400 
      Q2(45,64) = 0 

      Q2(46,1) = 0 
      Q2(46,2) = (153*35**(0.5D0)*ETAH)/358400 
      Q2(46,3) = 0 
      Q2(46,4) = -(81*15**(0.5D0)*ETAH)/256000 
      Q2(46,5) = -(99*35**(0.5D0)*MUH)/256000 
      Q2(46,6) = (3*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(595*VT - 
     >   510*ABS(ETAH) + 462*ABS(MUH) - 
     >   1190*ABS(XIH)))/3584000 
      Q2(46,7) = (333*7**(0.5D0)*MUH)/25600 
      Q2(46,8) = (9*5**(0.5D0)*(54*ABS(ETAH) - 63*VT + 
     >   154*ABS(MUH) + 126*ABS(XIH)))/512000 
      Q2(46,9) = 0 
      Q2(46,10) = (153*7**(0.5D0)*ETAH)/8960 
      Q2(46,11) = 0 
      Q2(46,12) = -(81*3**(0.5D0)*ETAH)/6400 
      Q2(46,13) = (891*15**(0.5D0)*MUH)/1792000 
      Q2(46,14) = -(9*5**(0.5D0)*(765*VT + 1190*ABS(ETAH) + 
     >   594*ABS(MUH) - 1530*ABS(XIH)))/3584000 
      Q2(46,15) = -(2997*3**(0.5D0)*MUH)/179200 
      Q2(46,16) = (81*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(9*VT + 
     >   14*ABS(ETAH) - 22*ABS(MUH) - 
     >   18*ABS(XIH)))/3584000 
      Q2(46,17) = 0 
      Q2(46,18) = 0 
      Q2(46,19) = 0 
      Q2(46,20) = 0 
      Q2(46,21) = 0 
      Q2(46,22) = (153*35**(0.5D0)*XIH)/25600 
      Q2(46,23) = 0 
      Q2(46,24) = -(567*15**(0.5D0)*XIH)/128000 
      Q2(46,25) = 0 
      Q2(46,26) = 0 
      Q2(46,27) = 0 
      Q2(46,28) = 0 
      Q2(46,29) = 0 
      Q2(46,30) = -(1377*15**(0.5D0)*XIH)/179200 
      Q2(46,31) = 0 
      Q2(46,32) = (2187*35**(0.5D0)*XIH)/896000 
      Q2(46,33) = 0 
      Q2(46,34) = -(459*7**(0.5D0)*ETAH)/358400 
      Q2(46,35) = 0 
      Q2(46,36) = (243*3**(0.5D0)*ETAH)/256000 
      Q2(46,37) = (297*7**(0.5D0)*MUH)/256000 
      Q2(46,38) = -(3*3**(0.5D0)*7**(0.5D0)*(1785*VT - 
     >   1530*ABS(ETAH) + 1386*ABS(MUH) + 
     >   5950*ABS(XIH)))/3584000 
      Q2(46,39) = -(999*35**(0.5D0)*MUH)/128000 
      Q2(46,40) = (1701*VT)/512000 - (729*ABS(ETAH))/256000 - 
     >   (2079*ABS(MUH))/256000 + (567*ABS(XIH))/51200 
      Q2(46,41) = 0 
      Q2(46,42) = -(459*35**(0.5D0)*ETAH)/44800 
      Q2(46,43) = 0 
      Q2(46,44) = (243*15**(0.5D0)*ETAH)/32000 
      Q2(46,45) = -(2673*3**(0.5D0)*MUH)/1792000 
      Q2(46,46) = (4131*VT)/716800 + (459*ABS(ETAH))/51200 + 
     >   (8019*ABS(MUH))/1792000 + (1377*ABS(XIH))/71680 
      Q2(46,47) = (8991*15**(0.5D0)*MUH)/896000 
      Q2(46,48) = -(243*3**(0.5D0)*7**(0.5D0)*(9*VT + 
     >   14*ABS(ETAH) - 22*ABS(MUH) + 
     >   30*ABS(XIH)))/3584000 
      Q2(46,49) = 0 
      Q2(46,50) = 0 
      Q2(46,51) = 0 
      Q2(46,52) = 0 
      Q2(46,53) = 0 
      Q2(46,54) = -(357*15**(0.5D0)*XIH)/51200 
      Q2(46,55) = 0 
      Q2(46,56) = (567*35**(0.5D0)*XIH)/256000 
      Q2(46,57) = 0 
      Q2(46,58) = 0 
      Q2(46,59) = 0 
      Q2(46,60) = 0 
      Q2(46,61) = 0 
      Q2(46,62) = (1377*35**(0.5D0)*XIH)/358400 
      Q2(46,63) = 0 
      Q2(46,64) = -(729*15**(0.5D0)*XIH)/256000 

      Q2(47,1) = -(81*7**(0.5D0)*ETAH)/358400 
      Q2(47,2) = 0 
      Q2(47,3) = (243*35**(0.5D0)*ETAH)/1792000 
      Q2(47,4) = 0 
      Q2(47,5) = (27*3**(0.5D0)*7**(0.5D0)*(6*ABS(ETAH) - 
     >   7*VT + 14*ABS(MUH) + 14*ABS(XIH)))/716800 
      Q2(47,6) = -(81*7**(0.5D0)*MUH)/25600 
      Q2(47,7) = (27*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(21*VT - 
     >   18*ABS(ETAH) + 70*ABS(MUH) - 
     >   42*ABS(XIH)))/3584000 
      Q2(47,8) = (189*3**(0.5D0)*MUH)/51200 
      Q2(47,9) = -(81*35**(0.5D0)*ETAH)/44800 
      Q2(47,10) = 0 
      Q2(47,11) = (243*7**(0.5D0)*ETAH)/44800 
      Q2(47,12) = 0 
      Q2(47,13) = (729*VT)/716800 + (81*ABS(ETAH))/51200 - 
     >   (729*ABS(MUH))/358400 - (729*ABS(XIH))/358400 
      Q2(47,14) = (729*3**(0.5D0)*MUH)/179200 
      Q2(47,15) = -(243*5**(0.5D0)*(9*VT + 14*ABS(ETAH) + 
     >   30*ABS(MUH) - 18*ABS(XIH)))/3584000 
      Q2(47,16) = -(729*7**(0.5D0)*MUH)/358400 
      Q2(47,17) = 0 
      Q2(47,18) = 0 
      Q2(47,19) = 0 
      Q2(47,20) = 0 
      Q2(47,21) = -(81*7**(0.5D0)*XIH)/25600 
      Q2(47,22) = 0 
      Q2(47,23) = (243*35**(0.5D0)*XIH)/128000 
      Q2(47,24) = 0 
      Q2(47,25) = 0 
      Q2(47,26) = 0 
      Q2(47,27) = 0 
      Q2(47,28) = 0 
      Q2(47,29) = (729*3**(0.5D0)*XIH)/179200 
      Q2(47,30) = 0 
      Q2(47,31) = -(2187*15**(0.5D0)*XIH)/896000 
      Q2(47,32) = 0 
      Q2(47,33) = (243*35**(0.5D0)*ETAH)/1792000 
      Q2(47,34) = 0 
      Q2(47,35) = -(729*7**(0.5D0)*ETAH)/1792000 
      Q2(47,36) = 0 
      Q2(47,37) = (27*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(21*VT - 
     >   18*ABS(ETAH) - 42*ABS(MUH) + 
     >   70*ABS(XIH)))/3584000 
      Q2(47,38) = (243*35**(0.5D0)*MUH)/128000 
      Q2(47,39) = -(81*3**(0.5D0)*7**(0.5D0)*(21*VT - 
     >   18*ABS(ETAH) + 70*ABS(MUH) + 
     >   70*ABS(XIH)))/3584000 
      Q2(47,40) = -(567*15**(0.5D0)*MUH)/256000 
      Q2(47,41) = (243*7**(0.5D0)*ETAH)/44800 
      Q2(47,42) = 0 
      Q2(47,43) = -(729*35**(0.5D0)*ETAH)/224000 
      Q2(47,44) = 0 
      Q2(47,45) = -(243*5**(0.5D0)*(9*VT + 14*ABS(ETAH) - 
     >   18*ABS(MUH) + 30*ABS(XIH)))/3584000 
      Q2(47,46) = -(2187*15**(0.5D0)*MUH)/896000 
      Q2(47,47) = (6561*VT)/3584000 + (729*ABS(ETAH))/256000 + 
     >   (2187*ABS(MUH))/358400 + (2187*ABS(XIH))/358400 
      Q2(47,48) = (2187*35**(0.5D0)*MUH)/1792000 
      Q2(47,49) = 0 
      Q2(47,50) = 0 
      Q2(47,51) = 0 
      Q2(47,52) = 0 
      Q2(47,53) = (189*3**(0.5D0)*XIH)/51200 
      Q2(47,54) = 0 
      Q2(47,55) = -(567*15**(0.5D0)*XIH)/256000 
      Q2(47,56) = 0 
      Q2(47,57) = 0 
      Q2(47,58) = 0 
      Q2(47,59) = 0 
      Q2(47,60) = 0 
      Q2(47,61) = -(729*7**(0.5D0)*XIH)/358400 
      Q2(47,62) = 0 
      Q2(47,63) = (2187*35**(0.5D0)*XIH)/1792000 
      Q2(47,64) = 0 

      Q2(48,1) = 0 
      Q2(48,2) = -(81*15**(0.5D0)*ETAH)/256000 
      Q2(48,3) = 0 
      Q2(48,4) = (2187*35**(0.5D0)*ETAH)/12544000 
      Q2(48,5) = -(81*15**(0.5D0)*MUH)/256000 
      Q2(48,6) = (81*5**(0.5D0)*(6*ABS(ETAH) - 7*VT + 
     >   6*ABS(MUH) + 14*ABS(XIH)))/512000 
      Q2(48,7) = -(81*3**(0.5D0)*MUH)/6400 
      Q2(48,8) = (81*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(63*VT - 
     >   54*ABS(ETAH) + 98*ABS(MUH) - 
     >   126*ABS(XIH)))/25088000 
      Q2(48,9) = 0 
      Q2(48,10) = -(81*3**(0.5D0)*ETAH)/6400 
      Q2(48,11) = 0 
      Q2(48,12) = (2187*7**(0.5D0)*ETAH)/313600 
      Q2(48,13) = (2187*35**(0.5D0)*MUH)/12544000 
      Q2(48,14) = (81*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(63*VT + 
     >   98*ABS(ETAH) - 54*ABS(MUH) - 
     >   126*ABS(XIH)))/25088000 
      Q2(48,15) = (2187*7**(0.5D0)*MUH)/313600 
      Q2(48,16) = -(2187*5**(0.5D0)*(9*VT + 14*ABS(ETAH) + 
     >   14*ABS(MUH) - 18*ABS(XIH)))/25088000 
      Q2(48,17) = 0 
      Q2(48,18) = 0 
      Q2(48,19) = 0 
      Q2(48,20) = 0 
      Q2(48,21) = 0 
      Q2(48,22) = -(567*15**(0.5D0)*XIH)/128000 
      Q2(48,23) = 0 
      Q2(48,24) = (2187*35**(0.5D0)*XIH)/896000 
      Q2(48,25) = 0 
      Q2(48,26) = 0 
      Q2(48,27) = 0 
      Q2(48,28) = 0 
      Q2(48,29) = 0 
      Q2(48,30) = (2187*35**(0.5D0)*XIH)/896000 
      Q2(48,31) = 0 
      Q2(48,32) = -(19683*15**(0.5D0)*XIH)/6272000 
      Q2(48,33) = 0 
      Q2(48,34) = (243*3**(0.5D0)*ETAH)/256000 
      Q2(48,35) = 0 
      Q2(48,36) = -(6561*7**(0.5D0)*ETAH)/12544000 
      Q2(48,37) = (243*3**(0.5D0)*MUH)/256000 
      Q2(48,38) = (1701*VT)/512000 - (729*ABS(ETAH))/256000 - 
     >   (729*ABS(MUH))/256000 + (567*ABS(XIH))/51200 
      Q2(48,39) = (243*15**(0.5D0)*MUH)/32000 
      Q2(48,40) = -(243*3**(0.5D0)*7**(0.5D0)*(63*VT - 
     >   54*ABS(ETAH) + 98*ABS(MUH) + 
     >   210*ABS(XIH)))/25088000 
      Q2(48,41) = 0 
      Q2(48,42) = (243*15**(0.5D0)*ETAH)/32000 
      Q2(48,43) = 0 
      Q2(48,44) = -(6561*35**(0.5D0)*ETAH)/1568000 
      Q2(48,45) = -(6561*7**(0.5D0)*MUH)/12544000 
      Q2(48,46) = -(243*3**(0.5D0)*7**(0.5D0)*(63*VT + 
     >   98*ABS(ETAH) - 54*ABS(MUH) + 
     >   210*ABS(XIH)))/25088000 
      Q2(48,47) = -(6561*35**(0.5D0)*MUH)/1568000 
      Q2(48,48) = (59049*VT)/25088000 + 
     >   (6561*ABS(ETAH))/1792000 + 
     >   (6561*ABS(MUH))/1792000 + 
     >   (19683*ABS(XIH))/2508800 
      Q2(48,49) = 0 
      Q2(48,50) = 0 
      Q2(48,51) = 0 
      Q2(48,52) = 0 
      Q2(48,53) = 0 
      Q2(48,54) = (567*35**(0.5D0)*XIH)/256000 
      Q2(48,55) = 0 
      Q2(48,56) = -(729*15**(0.5D0)*XIH)/256000 
      Q2(48,57) = 0 
      Q2(48,58) = 0 
      Q2(48,59) = 0 
      Q2(48,60) = 0 
      Q2(48,61) = 0 
      Q2(48,62) = -(729*15**(0.5D0)*XIH)/256000 
      Q2(48,63) = 0 
      Q2(48,64) = (19683*35**(0.5D0)*XIH)/12544000 

      Q2(49,1) = -(45*7**(0.5D0)*XIH)/14336 
      Q2(49,2) = 0 
      Q2(49,3) = (27*35**(0.5D0)*XIH)/71680 
      Q2(49,4) = 0 
      Q2(49,5) = 0 
      Q2(49,6) = 0 
      Q2(49,7) = 0 
      Q2(49,8) = 0 
      Q2(49,9) = (27*35**(0.5D0)*XIH)/71680 
      Q2(49,10) = 0 
      Q2(49,11) = -(81*7**(0.5D0)*XIH)/358400 
      Q2(49,12) = 0 
      Q2(49,13) = 0 
      Q2(49,14) = 0 
      Q2(49,15) = 0 
      Q2(49,16) = 0 
      Q2(49,17) = -(3*3**(0.5D0)*7**(0.5D0)*(35*VT + 
     >   14*ABS(ETAH) + 14*ABS(MUH) - 
     >   30*ABS(XIH)))/28672 
      Q2(49,18) = -(9*7**(0.5D0)*MUH)/512 
      Q2(49,19) = (3*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(105*VT + 
     >   42*ABS(ETAH) - 350*ABS(MUH) - 
     >   90*ABS(XIH)))/716800 
      Q2(49,20) = -(21*3**(0.5D0)*MUH)/2048 
      Q2(49,21) = -(9*7**(0.5D0)*ETAH)/512 
      Q2(49,22) = 0 
      Q2(49,23) = (27*35**(0.5D0)*ETAH)/12800 
      Q2(49,24) = 0 
      Q2(49,25) = (3*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(105*VT - 
     >   350*ABS(ETAH) + 42*ABS(MUH) - 
     >   90*ABS(XIH)))/716800 
      Q2(49,26) = (27*35**(0.5D0)*MUH)/12800 
      Q2(49,27) = (9*3**(0.5D0)*7**(0.5D0)*(70*ABS(ETAH) - 
     >   21*VT + 70*ABS(MUH) + 18*ABS(XIH)))/716800 
      Q2(49,28) = (63*15**(0.5D0)*MUH)/51200 
      Q2(49,29) = -(21*3**(0.5D0)*ETAH)/2048 
      Q2(49,30) = 0 
      Q2(49,31) = (63*15**(0.5D0)*ETAH)/51200 
      Q2(49,32) = 0 
      Q2(49,33) = -(45*35**(0.5D0)*XIH)/1792 
      Q2(49,34) = 0 
      Q2(49,35) = (27*7**(0.5D0)*XIH)/1792 
      Q2(49,36) = 0 
      Q2(49,37) = 0 
      Q2(49,38) = 0 
      Q2(49,39) = 0 
      Q2(49,40) = 0 
      Q2(49,41) = (27*7**(0.5D0)*XIH)/1792 
      Q2(49,42) = 0 
      Q2(49,43) = -(81*35**(0.5D0)*XIH)/44800 
      Q2(49,44) = 0 
      Q2(49,45) = 0 
      Q2(49,46) = 0 
      Q2(49,47) = 0 
      Q2(49,48) = 0 
      Q2(49,49) = (405*VT)/28672 + (81*ABS(ETAH))/14336 + 
     >   (81*ABS(MUH))/14336 + (45*ABS(XIH))/2048 
      Q2(49,50) = (81*3**(0.5D0)*MUH)/3584 
      Q2(49,51) = -(27*5**(0.5D0)*(45*VT + 18*ABS(ETAH) - 
     >   150*ABS(MUH) + 70*ABS(XIH)))/716800 
      Q2(49,52) = (81*7**(0.5D0)*MUH)/14336 
      Q2(49,53) = (81*3**(0.5D0)*ETAH)/3584 
      Q2(49,54) = 0 
      Q2(49,55) = -(243*15**(0.5D0)*ETAH)/89600 
      Q2(49,56) = 0 
      Q2(49,57) = -(27*5**(0.5D0)*(45*VT - 150*ABS(ETAH) + 
     >   18*ABS(MUH) + 70*ABS(XIH)))/716800 
      Q2(49,58) = -(243*15**(0.5D0)*MUH)/89600 
      Q2(49,59) = (729*VT)/716800 - (243*ABS(ETAH))/71680 - 
     >   (243*ABS(MUH))/71680 + (81*ABS(XIH))/51200 
      Q2(49,60) = -(243*35**(0.5D0)*MUH)/358400 
      Q2(49,61) = (81*7**(0.5D0)*ETAH)/14336 
      Q2(49,62) = 0 
      Q2(49,63) = -(243*35**(0.5D0)*ETAH)/358400 
      Q2(49,64) = 0 

      Q2(50,1) = 0 
      Q2(50,2) = -(51*7**(0.5D0)*XIH)/14336 
      Q2(50,3) = 0 
      Q2(50,4) = (27*3**(0.5D0)*XIH)/10240 
      Q2(50,5) = 0 
      Q2(50,6) = 0 
      Q2(50,7) = 0 
      Q2(50,8) = 0 
      Q2(50,9) = 0 
      Q2(50,10) = (153*35**(0.5D0)*XIH)/358400 
      Q2(50,11) = 0 
      Q2(50,12) = -(81*15**(0.5D0)*XIH)/256000 
      Q2(50,13) = 0 
      Q2(50,14) = 0 
      Q2(50,15) = 0 
      Q2(50,16) = 0 
      Q2(50,17) = (33*7**(0.5D0)*MUH)/10240 
      Q2(50,18) = -(3**(0.5D0)*7**(0.5D0)*(595*VT + 
     >   238*ABS(ETAH) + 462*ABS(MUH) - 
     >   510*ABS(XIH)))/143360 
      Q2(50,19) = -(111*35**(0.5D0)*MUH)/5120 
      Q2(50,20) = (189*VT)/20480 + (189*ABS(ETAH))/51200 - 
     >   (231*ABS(MUH))/10240 - (81*ABS(XIH))/10240 
      Q2(50,21) = 0 
      Q2(50,22) = -(51*7**(0.5D0)*ETAH)/2560 
      Q2(50,23) = 0 
      Q2(50,24) = (189*3**(0.5D0)*ETAH)/12800 
      Q2(50,25) = -(99*35**(0.5D0)*MUH)/256000 
      Q2(50,26) = (3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(1785*VT - 
     >   5950*ABS(ETAH) + 1386*ABS(MUH) - 
     >   1530*ABS(XIH)))/3584000 
      Q2(50,27) = (333*7**(0.5D0)*MUH)/25600 
      Q2(50,28) = (9*5**(0.5D0)*(210*ABS(ETAH) - 63*VT + 
     >   154*ABS(MUH) + 54*ABS(XIH)))/512000 
      Q2(50,29) = 0 
      Q2(50,30) = -(119*3**(0.5D0)*ETAH)/10240 
      Q2(50,31) = 0 
      Q2(50,32) = (189*7**(0.5D0)*ETAH)/51200 
      Q2(50,33) = 0 
      Q2(50,34) = -(51*35**(0.5D0)*XIH)/1792 
      Q2(50,35) = 0 
      Q2(50,36) = (27*15**(0.5D0)*XIH)/1280 
      Q2(50,37) = 0 
      Q2(50,38) = 0 
      Q2(50,39) = 0 
      Q2(50,40) = 0 
      Q2(50,41) = 0 
      Q2(50,42) = (153*7**(0.5D0)*XIH)/8960 
      Q2(50,43) = 0 
      Q2(50,44) = -(81*3**(0.5D0)*XIH)/6400 
      Q2(50,45) = 0 
      Q2(50,46) = 0 
      Q2(50,47) = 0 
      Q2(50,48) = 0 
      Q2(50,49) = -(297*3**(0.5D0)*MUH)/71680 
      Q2(50,50) = (459*VT)/28672 + (459*ABS(ETAH))/71680 + 
     >   (891*ABS(MUH))/71680 + (51*ABS(XIH))/2048 
      Q2(50,51) = (999*15**(0.5D0)*MUH)/35840 
      Q2(50,52) = -(27*3**(0.5D0)*7**(0.5D0)*(45*VT + 
     >   18*ABS(ETAH) - 110*ABS(MUH) + 
     >   70*ABS(XIH)))/716800 
      Q2(50,53) = 0 
      Q2(50,54) = (459*3**(0.5D0)*ETAH)/17920 
      Q2(50,55) = 0 
      Q2(50,56) = -(729*7**(0.5D0)*ETAH)/89600 
      Q2(50,57) = (891*15**(0.5D0)*MUH)/1792000 
      Q2(50,58) = -(9*5**(0.5D0)*(765*VT - 2550*ABS(ETAH) + 
     >   594*ABS(MUH) + 1190*ABS(XIH)))/3584000 
      Q2(50,59) = -(2997*3**(0.5D0)*MUH)/179200 
      Q2(50,60) = (81*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(9*VT - 
     >   30*ABS(ETAH) - 22*ABS(MUH) + 
     >   14*ABS(XIH)))/3584000 
      Q2(50,61) = 0 
      Q2(50,62) = (459*7**(0.5D0)*ETAH)/71680 
      Q2(50,63) = 0 
      Q2(50,64) = -(243*3**(0.5D0)*ETAH)/51200 

      Q2(51,1) = (27*35**(0.5D0)*XIH)/71680 
      Q2(51,2) = 0 
      Q2(51,3) = -(81*7**(0.5D0)*XIH)/71680 
      Q2(51,4) = 0 
      Q2(51,5) = 0 
      Q2(51,6) = 0 
      Q2(51,7) = 0 
      Q2(51,8) = 0 
      Q2(51,9) = -(81*7**(0.5D0)*XIH)/358400 
      Q2(51,10) = 0 
      Q2(51,11) = (243*35**(0.5D0)*XIH)/1792000 
      Q2(51,12) = 0 
      Q2(51,13) = 0 
      Q2(51,14) = 0 
      Q2(51,15) = 0 
      Q2(51,16) = 0 
      Q2(51,17) = (9*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(35*VT + 
     >   14*ABS(ETAH) - 70*ABS(MUH) - 
     >   30*ABS(XIH)))/716800 
      Q2(51,18) = (27*35**(0.5D0)*MUH)/5120 
      Q2(51,19) = -(9*3**(0.5D0)*7**(0.5D0)*(105*VT + 
     >   42*ABS(ETAH) + 350*ABS(MUH) - 
     >   90*ABS(XIH)))/716800 
      Q2(51,20) = -(63*15**(0.5D0)*MUH)/10240 
      Q2(51,21) = (27*35**(0.5D0)*ETAH)/12800 
      Q2(51,22) = 0 
      Q2(51,23) = -(81*7**(0.5D0)*ETAH)/12800 
      Q2(51,24) = 0 
      Q2(51,25) = (9*3**(0.5D0)*7**(0.5D0)*(70*ABS(ETAH) - 
     >   21*VT + 42*ABS(MUH) + 18*ABS(XIH)))/716800 
      Q2(51,26) = -(81*7**(0.5D0)*MUH)/25600 
      Q2(51,27) = (27*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(21*VT - 
     >   70*ABS(ETAH) + 70*ABS(MUH) - 
     >   18*ABS(XIH)))/3584000 
      Q2(51,28) = (189*3**(0.5D0)*MUH)/51200 
      Q2(51,29) = (63*15**(0.5D0)*ETAH)/51200 
      Q2(51,30) = 0 
      Q2(51,31) = -(189*3**(0.5D0)*ETAH)/51200 
      Q2(51,32) = 0 
      Q2(51,33) = (27*7**(0.5D0)*XIH)/1792 
      Q2(51,34) = 0 
      Q2(51,35) = -(81*35**(0.5D0)*XIH)/8960 
      Q2(51,36) = 0 
      Q2(51,37) = 0 
      Q2(51,38) = 0 
      Q2(51,39) = 0 
      Q2(51,40) = 0 
      Q2(51,41) = -(81*35**(0.5D0)*XIH)/44800 
      Q2(51,42) = 0 
      Q2(51,43) = (243*7**(0.5D0)*XIH)/44800 
      Q2(51,44) = 0 
      Q2(51,45) = 0 
      Q2(51,46) = 0 
      Q2(51,47) = 0 
      Q2(51,48) = 0 
      Q2(51,49) = -(27*5**(0.5D0)*(45*VT + 18*ABS(ETAH) - 
     >   90*ABS(MUH) + 70*ABS(XIH)))/716800 
      Q2(51,50) = -(243*15**(0.5D0)*MUH)/35840 
      Q2(51,51) = (729*VT)/143360 + (729*ABS(ETAH))/358400 + 
     >   (243*ABS(MUH))/14336 + (81*ABS(XIH))/10240 
      Q2(51,52) = (243*35**(0.5D0)*MUH)/71680 
      Q2(51,53) = -(243*15**(0.5D0)*ETAH)/89600 
      Q2(51,54) = 0 
      Q2(51,55) = (729*3**(0.5D0)*ETAH)/89600 
      Q2(51,56) = 0 
      Q2(51,57) = (729*VT)/716800 - (243*ABS(ETAH))/71680 - 
     >   (729*ABS(MUH))/358400 + (81*ABS(XIH))/51200 
      Q2(51,58) = (729*3**(0.5D0)*MUH)/179200 
      Q2(51,59) = -(243*5**(0.5D0)*(9*VT - 30*ABS(ETAH) + 
     >   30*ABS(MUH) + 14*ABS(XIH)))/3584000 
      Q2(51,60) = -(729*7**(0.5D0)*MUH)/358400 
      Q2(51,61) = -(243*35**(0.5D0)*ETAH)/358400 
      Q2(51,62) = 0 
      Q2(51,63) = (729*7**(0.5D0)*ETAH)/358400 
      Q2(51,64) = 0 

      Q2(52,1) = 0 
      Q2(52,2) = (27*3**(0.5D0)*XIH)/10240 
      Q2(52,3) = 0 
      Q2(52,4) = -(729*7**(0.5D0)*XIH)/501760 
      Q2(52,5) = 0 
      Q2(52,6) = 0 
      Q2(52,7) = 0 
      Q2(52,8) = 0 
      Q2(52,9) = 0 
      Q2(52,10) = -(81*15**(0.5D0)*XIH)/256000 
      Q2(52,11) = 0 
      Q2(52,12) = (2187*35**(0.5D0)*XIH)/12544000 
      Q2(52,13) = 0 
      Q2(52,14) = 0 
      Q2(52,15) = 0 
      Q2(52,16) = 0 
      Q2(52,17) = (27*3**(0.5D0)*MUH)/10240 
      Q2(52,18) = (189*VT)/20480 + (189*ABS(ETAH))/51200 - 
     >   (81*ABS(MUH))/10240 - (81*ABS(XIH))/10240 
      Q2(52,19) = (27*15**(0.5D0)*MUH)/1280 
      Q2(52,20) = -(27*3**(0.5D0)*7**(0.5D0)*(315*VT + 
     >   126*ABS(ETAH) + 490*ABS(MUH) - 
     >   270*ABS(XIH)))/5017600 
      Q2(52,21) = 0 
      Q2(52,22) = (189*3**(0.5D0)*ETAH)/12800 
      Q2(52,23) = 0 
      Q2(52,24) = -(729*7**(0.5D0)*ETAH)/89600 
      Q2(52,25) = -(81*15**(0.5D0)*MUH)/256000 
      Q2(52,26) = (27*5**(0.5D0)*(70*ABS(ETAH) - 21*VT + 
     >   18*ABS(MUH) + 18*ABS(XIH)))/512000 
      Q2(52,27) = -(81*3**(0.5D0)*MUH)/6400 
      Q2(52,28) = (81*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(63*VT - 
     >   210*ABS(ETAH) + 98*ABS(MUH) - 
     >   54*ABS(XIH)))/25088000 
      Q2(52,29) = 0 
      Q2(52,30) = (189*7**(0.5D0)*ETAH)/51200 
      Q2(52,31) = 0 
      Q2(52,32) = -(243*3**(0.5D0)*ETAH)/51200 
      Q2(52,33) = 0 
      Q2(52,34) = (27*15**(0.5D0)*XIH)/1280 
      Q2(52,35) = 0 
      Q2(52,36) = -(729*35**(0.5D0)*XIH)/62720 
      Q2(52,37) = 0 
      Q2(52,38) = 0 
      Q2(52,39) = 0 
      Q2(52,40) = 0 
      Q2(52,41) = 0 
      Q2(52,42) = -(81*3**(0.5D0)*XIH)/6400 
      Q2(52,43) = 0 
      Q2(52,44) = (2187*7**(0.5D0)*XIH)/313600 
      Q2(52,45) = 0 
      Q2(52,46) = 0 
      Q2(52,47) = 0 
      Q2(52,48) = 0 
      Q2(52,49) = -(729*7**(0.5D0)*MUH)/501760 
      Q2(52,50) = -(27*3**(0.5D0)*7**(0.5D0)*(315*VT + 
     >   126*ABS(ETAH) - 270*ABS(MUH) + 
     >   490*ABS(XIH)))/5017600 
      Q2(52,51) = -(729*35**(0.5D0)*MUH)/62720 
      Q2(52,52) = (6561*VT)/1003520 + (6561*ABS(ETAH))/2508800 + 
     >   (729*ABS(MUH))/71680 + (729*ABS(XIH))/71680 
      Q2(52,53) = 0 
      Q2(52,54) = -(729*7**(0.5D0)*ETAH)/89600 
      Q2(52,55) = 0 
      Q2(52,56) = (6561*3**(0.5D0)*ETAH)/627200 
      Q2(52,57) = (2187*35**(0.5D0)*MUH)/12544000 
      Q2(52,58) = (81*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(63*VT - 
     >   210*ABS(ETAH) - 54*ABS(MUH) + 
     >   98*ABS(XIH)))/25088000 
      Q2(52,59) = (2187*7**(0.5D0)*MUH)/313600 
      Q2(52,60) = -(2187*5**(0.5D0)*(9*VT - 30*ABS(ETAH) + 
     >   14*ABS(MUH) + 14*ABS(XIH)))/25088000 
      Q2(52,61) = 0 
      Q2(52,62) = -(243*3**(0.5D0)*ETAH)/51200 
      Q2(52,63) = 0 
      Q2(52,64) = (6561*7**(0.5D0)*ETAH)/2508800 

      Q2(53,1) = 0 
      Q2(53,2) = 0 
      Q2(53,3) = 0 
      Q2(53,4) = 0 
      Q2(53,5) = -(51*7**(0.5D0)*XIH)/14336 
      Q2(53,6) = 0 
      Q2(53,7) = (153*35**(0.5D0)*XIH)/358400 
      Q2(53,8) = 0 
      Q2(53,9) = 0 
      Q2(53,10) = 0 
      Q2(53,11) = 0 
      Q2(53,12) = 0 
      Q2(53,13) = (27*3**(0.5D0)*XIH)/10240 
      Q2(53,14) = 0 
      Q2(53,15) = -(81*15**(0.5D0)*XIH)/256000 
      Q2(53,16) = 0 
      Q2(53,17) = (33*7**(0.5D0)*ETAH)/10240 
      Q2(53,18) = 0 
      Q2(53,19) = -(99*35**(0.5D0)*ETAH)/256000 
      Q2(53,20) = 0 
      Q2(53,21) = -(3**(0.5D0)*7**(0.5D0)*(595*VT + 
     >   462*ABS(ETAH) + 238*ABS(MUH) - 
     >   510*ABS(XIH)))/143360 
      Q2(53,22) = -(51*7**(0.5D0)*MUH)/2560 
      Q2(53,23) = (3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(1785*VT + 
     >   1386*ABS(ETAH) - 5950*ABS(MUH) - 
     >   1530*ABS(XIH)))/3584000 
      Q2(53,24) = -(119*3**(0.5D0)*MUH)/10240 
      Q2(53,25) = -(111*35**(0.5D0)*ETAH)/5120 
      Q2(53,26) = 0 
      Q2(53,27) = (333*7**(0.5D0)*ETAH)/25600 
      Q2(53,28) = 0 
      Q2(53,29) = (189*VT)/20480 - (231*ABS(ETAH))/10240 + 
     >   (189*ABS(MUH))/51200 - (81*ABS(XIH))/10240 
      Q2(53,30) = (189*3**(0.5D0)*MUH)/12800 
      Q2(53,31) = (9*5**(0.5D0)*(154*ABS(ETAH) - 63*VT + 
     >   210*ABS(MUH) + 54*ABS(XIH)))/512000 
      Q2(53,32) = (189*7**(0.5D0)*MUH)/51200 
      Q2(53,33) = 0 
      Q2(53,34) = 0 
      Q2(53,35) = 0 
      Q2(53,36) = 0 
      Q2(53,37) = -(51*35**(0.5D0)*XIH)/1792 
      Q2(53,38) = 0 
      Q2(53,39) = (153*7**(0.5D0)*XIH)/8960 
      Q2(53,40) = 0 
      Q2(53,41) = 0 
      Q2(53,42) = 0 
      Q2(53,43) = 0 
      Q2(53,44) = 0 
      Q2(53,45) = (27*15**(0.5D0)*XIH)/1280 
      Q2(53,46) = 0 
      Q2(53,47) = -(81*3**(0.5D0)*XIH)/6400 
      Q2(53,48) = 0 
      Q2(53,49) = -(297*3**(0.5D0)*ETAH)/71680 
      Q2(53,50) = 0 
      Q2(53,51) = (891*15**(0.5D0)*ETAH)/1792000 
      Q2(53,52) = 0 
      Q2(53,53) = (459*VT)/28672 + (891*ABS(ETAH))/71680 + 
     >   (459*ABS(MUH))/71680 + (51*ABS(XIH))/2048 
      Q2(53,54) = (459*3**(0.5D0)*MUH)/17920 
      Q2(53,55) = -(9*5**(0.5D0)*(765*VT + 594*ABS(ETAH) - 
     >   2550*ABS(MUH) + 1190*ABS(XIH)))/3584000 
      Q2(53,56) = (459*7**(0.5D0)*MUH)/71680 
      Q2(53,57) = (999*15**(0.5D0)*ETAH)/35840 
      Q2(53,58) = 0 
      Q2(53,59) = -(2997*3**(0.5D0)*ETAH)/179200 
      Q2(53,60) = 0 
      Q2(53,61) = -(27*3**(0.5D0)*7**(0.5D0)*(45*VT - 
     >   110*ABS(ETAH) + 18*ABS(MUH) + 
     >   70*ABS(XIH)))/716800 
      Q2(53,62) = -(729*7**(0.5D0)*MUH)/89600 
      Q2(53,63) = (81*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(9*VT - 
     >   22*ABS(ETAH) - 30*ABS(MUH) + 
     >   14*ABS(XIH)))/3584000 
      Q2(53,64) = -(243*3**(0.5D0)*MUH)/51200 

      Q2(54,1) = 0 
      Q2(54,2) = 0 
      Q2(54,3) = 0 
      Q2(54,4) = 0 
      Q2(54,5) = 0 
      Q2(54,6) = -(289*7**(0.5D0)*XIH)/71680 
      Q2(54,7) = 0 
      Q2(54,8) = (153*3**(0.5D0)*XIH)/51200 
      Q2(54,9) = 0 
      Q2(54,10) = 0 
      Q2(54,11) = 0 
      Q2(54,12) = 0 
      Q2(54,13) = 0 
      Q2(54,14) = (153*3**(0.5D0)*XIH)/51200 
      Q2(54,15) = 0 
      Q2(54,16) = -(243*7**(0.5D0)*XIH)/256000 
      Q2(54,17) = 0 
      Q2(54,18) = (187*7**(0.5D0)*ETAH)/51200 
      Q2(54,19) = 0 
      Q2(54,20) = -(693*3**(0.5D0)*ETAH)/256000 
      Q2(54,21) = (187*7**(0.5D0)*MUH)/51200 
      Q2(54,22) = -(17*3**(0.5D0)*7**(0.5D0)*(595*VT + 
     >   462*ABS(ETAH) + 462*ABS(MUH) - 
     >   510*ABS(XIH)))/2150400 
      Q2(54,23) = -(629*35**(0.5D0)*MUH)/25600 
      Q2(54,24) = (1071*VT)/102400 + (2079*ABS(ETAH))/256000 - 
     >   (1309*ABS(MUH))/51200 - (459*ABS(XIH))/51200 
      Q2(54,25) = 0 
      Q2(54,26) = -(629*35**(0.5D0)*ETAH)/25600 
      Q2(54,27) = 0 
      Q2(54,28) = (2331*15**(0.5D0)*ETAH)/128000 
      Q2(54,29) = -(693*3**(0.5D0)*MUH)/256000 
      Q2(54,30) = (1071*VT)/102400 - (1309*ABS(ETAH))/51200 + 
     >   (2079*ABS(MUH))/256000 - (459*ABS(XIH))/51200 
      Q2(54,31) = (2331*15**(0.5D0)*MUH)/128000 
      Q2(54,32) = (9*3**(0.5D0)*7**(0.5D0)*(154*ABS(ETAH) - 
     >   63*VT + 154*ABS(MUH) + 54*ABS(XIH)))/512000 
      Q2(54,33) = 0 
      Q2(54,34) = 0 
      Q2(54,35) = 0 
      Q2(54,36) = 0 
      Q2(54,37) = 0 
      Q2(54,38) = -(289*35**(0.5D0)*XIH)/8960 
      Q2(54,39) = 0 
      Q2(54,40) = (153*15**(0.5D0)*XIH)/6400 
      Q2(54,41) = 0 
      Q2(54,42) = 0 
      Q2(54,43) = 0 
      Q2(54,44) = 0 
      Q2(54,45) = 0 
      Q2(54,46) = (153*15**(0.5D0)*XIH)/6400 
      Q2(54,47) = 0 
      Q2(54,48) = -(243*35**(0.5D0)*XIH)/32000 
      Q2(54,49) = 0 
      Q2(54,50) = -(1683*3**(0.5D0)*ETAH)/358400 
      Q2(54,51) = 0 
      Q2(54,52) = (2673*7**(0.5D0)*ETAH)/1792000 
      Q2(54,53) = -(1683*3**(0.5D0)*MUH)/358400 
      Q2(54,54) = (2601*VT)/143360 + (5049*ABS(ETAH))/358400 + 
     >   (5049*ABS(MUH))/358400 + (289*ABS(XIH))/10240 
      Q2(54,55) = (5661*15**(0.5D0)*MUH)/179200 
      Q2(54,56) = -(9*3**(0.5D0)*7**(0.5D0)*(765*VT + 
     >   594*ABS(ETAH) - 1870*ABS(MUH) + 
     >   1190*ABS(XIH)))/3584000 
      Q2(54,57) = 0 
      Q2(54,58) = (5661*15**(0.5D0)*ETAH)/179200 
      Q2(54,59) = 0 
      Q2(54,60) = -(8991*35**(0.5D0)*ETAH)/896000 
      Q2(54,61) = (2673*7**(0.5D0)*MUH)/1792000 
      Q2(54,62) = -(9*3**(0.5D0)*7**(0.5D0)*(765*VT - 
     >   1870*ABS(ETAH) + 594*ABS(MUH) + 
     >   1190*ABS(XIH)))/3584000 
      Q2(54,63) = -(8991*35**(0.5D0)*MUH)/896000 
      Q2(54,64) = (2187*VT)/512000 - (2673*ABS(ETAH))/256000 - 
     >   (2673*ABS(MUH))/256000 + (1701*ABS(XIH))/256000 

      Q2(55,1) = 0 
      Q2(55,2) = 0 
      Q2(55,3) = 0 
      Q2(55,4) = 0 
      Q2(55,5) = (153*35**(0.5D0)*XIH)/358400 
      Q2(55,6) = 0 
      Q2(55,7) = -(459*7**(0.5D0)*XIH)/358400 
      Q2(55,8) = 0 
      Q2(55,9) = 0 
      Q2(55,10) = 0 
      Q2(55,11) = 0 
      Q2(55,12) = 0 
      Q2(55,13) = -(81*15**(0.5D0)*XIH)/256000 
      Q2(55,14) = 0 
      Q2(55,15) = (243*3**(0.5D0)*XIH)/256000 
      Q2(55,16) = 0 
      Q2(55,17) = -(99*35**(0.5D0)*ETAH)/256000 
      Q2(55,18) = 0 
      Q2(55,19) = (297*7**(0.5D0)*ETAH)/256000 
      Q2(55,20) = 0 
      Q2(55,21) = (3*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(595*VT + 
     >   462*ABS(ETAH) - 1190*ABS(MUH) - 
     >   510*ABS(XIH)))/3584000 
      Q2(55,22) = (153*35**(0.5D0)*MUH)/25600 
      Q2(55,23) = -(3*3**(0.5D0)*7**(0.5D0)*(1785*VT + 
     >   1386*ABS(ETAH) + 5950*ABS(MUH) - 
     >   1530*ABS(XIH)))/3584000 
      Q2(55,24) = -(357*15**(0.5D0)*MUH)/51200 
      Q2(55,25) = (333*7**(0.5D0)*ETAH)/25600 
      Q2(55,26) = 0 
      Q2(55,27) = -(999*35**(0.5D0)*ETAH)/128000 
      Q2(55,28) = 0 
      Q2(55,29) = (9*5**(0.5D0)*(154*ABS(ETAH) - 63*VT + 
     >   126*ABS(MUH) + 54*ABS(XIH)))/512000 
      Q2(55,30) = -(567*15**(0.5D0)*MUH)/128000 
      Q2(55,31) = (1701*VT)/512000 - (2079*ABS(ETAH))/256000 + 
     >   (567*ABS(MUH))/51200 - (729*ABS(XIH))/256000 
      Q2(55,32) = (567*35**(0.5D0)*MUH)/256000 
      Q2(55,33) = 0 
      Q2(55,34) = 0 
      Q2(55,35) = 0 
      Q2(55,36) = 0 
      Q2(55,37) = (153*7**(0.5D0)*XIH)/8960 
      Q2(55,38) = 0 
      Q2(55,39) = -(459*35**(0.5D0)*XIH)/44800 
      Q2(55,40) = 0 
      Q2(55,41) = 0 
      Q2(55,42) = 0 
      Q2(55,43) = 0 
      Q2(55,44) = 0 
      Q2(55,45) = -(81*3**(0.5D0)*XIH)/6400 
      Q2(55,46) = 0 
      Q2(55,47) = (243*15**(0.5D0)*XIH)/32000 
      Q2(55,48) = 0 
      Q2(55,49) = (891*15**(0.5D0)*ETAH)/1792000 
      Q2(55,50) = 0 
      Q2(55,51) = -(2673*3**(0.5D0)*ETAH)/1792000 
      Q2(55,52) = 0 
      Q2(55,53) = -(9*5**(0.5D0)*(765*VT + 594*ABS(ETAH) - 
     >   1530*ABS(MUH) + 1190*ABS(XIH)))/3584000 
      Q2(55,54) = -(1377*15**(0.5D0)*MUH)/179200 
      Q2(55,55) = (4131*VT)/716800 + (8019*ABS(ETAH))/1792000 + 
     >   (1377*ABS(MUH))/71680 + (459*ABS(XIH))/51200 
      Q2(55,56) = (1377*35**(0.5D0)*MUH)/358400 
      Q2(55,57) = -(2997*3**(0.5D0)*ETAH)/179200 
      Q2(55,58) = 0 
      Q2(55,59) = (8991*15**(0.5D0)*ETAH)/896000 
      Q2(55,60) = 0 
      Q2(55,61) = (81*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(9*VT - 
     >   22*ABS(ETAH) - 18*ABS(MUH) + 
     >   14*ABS(XIH)))/3584000 
      Q2(55,62) = (2187*35**(0.5D0)*MUH)/896000 
      Q2(55,63) = -(243*3**(0.5D0)*7**(0.5D0)*(9*VT - 
     >   22*ABS(ETAH) + 30*ABS(MUH) + 
     >   14*ABS(XIH)))/3584000 
      Q2(55,64) = -(729*15**(0.5D0)*MUH)/256000 

      Q2(56,1) = 0 
      Q2(56,2) = 0 
      Q2(56,3) = 0 
      Q2(56,4) = 0 
      Q2(56,5) = 0 
      Q2(56,6) = (153*3**(0.5D0)*XIH)/51200 
      Q2(56,7) = 0 
      Q2(56,8) = -(4131*7**(0.5D0)*XIH)/2508800 
      Q2(56,9) = 0 
      Q2(56,10) = 0 
      Q2(56,11) = 0 
      Q2(56,12) = 0 
      Q2(56,13) = 0 
      Q2(56,14) = -(243*7**(0.5D0)*XIH)/256000 
      Q2(56,15) = 0 
      Q2(56,16) = (2187*3**(0.5D0)*XIH)/1792000 
      Q2(56,17) = 0 
      Q2(56,18) = -(693*3**(0.5D0)*ETAH)/256000 
      Q2(56,19) = 0 
      Q2(56,20) = (2673*7**(0.5D0)*ETAH)/1792000 
      Q2(56,21) = (153*3**(0.5D0)*MUH)/51200 
      Q2(56,22) = (1071*VT)/102400 + (2079*ABS(ETAH))/256000 - 
     >   (459*ABS(MUH))/51200 - (459*ABS(XIH))/51200 
      Q2(56,23) = (153*15**(0.5D0)*MUH)/6400 
      Q2(56,24) = -(9*3**(0.5D0)*7**(0.5D0)*(5355*VT + 
     >   4158*ABS(ETAH) + 8330*ABS(MUH) - 
     >   4590*ABS(XIH)))/25088000 
      Q2(56,25) = 0 
      Q2(56,26) = (2331*15**(0.5D0)*ETAH)/128000 
      Q2(56,27) = 0 
      Q2(56,28) = -(8991*35**(0.5D0)*ETAH)/896000 
      Q2(56,29) = -(243*7**(0.5D0)*MUH)/256000 
      Q2(56,30) = (9*3**(0.5D0)*7**(0.5D0)*(154*ABS(ETAH) - 
     >   63*VT + 54*ABS(MUH) + 54*ABS(XIH)))/512000 
      Q2(56,31) = -(243*35**(0.5D0)*MUH)/32000 
      Q2(56,32) = (2187*VT)/512000 - (2673*ABS(ETAH))/256000 + 
     >   (1701*ABS(MUH))/256000 - (6561*ABS(XIH))/1792000 
      Q2(56,33) = 0 
      Q2(56,34) = 0 
      Q2(56,35) = 0 
      Q2(56,36) = 0 
      Q2(56,37) = 0 
      Q2(56,38) = (153*15**(0.5D0)*XIH)/6400 
      Q2(56,39) = 0 
      Q2(56,40) = -(4131*35**(0.5D0)*XIH)/313600 
      Q2(56,41) = 0 
      Q2(56,42) = 0 
      Q2(56,43) = 0 
      Q2(56,44) = 0 
      Q2(56,45) = 0 
      Q2(56,46) = -(243*35**(0.5D0)*XIH)/32000 
      Q2(56,47) = 0 
      Q2(56,48) = (2187*15**(0.5D0)*XIH)/224000 
      Q2(56,49) = 0 
      Q2(56,50) = (2673*7**(0.5D0)*ETAH)/1792000 
      Q2(56,51) = 0 
      Q2(56,52) = -(24057*3**(0.5D0)*ETAH)/12544000 
      Q2(56,53) = -(4131*7**(0.5D0)*MUH)/2508800 
      Q2(56,54) = -(9*3**(0.5D0)*7**(0.5D0)*(5355*VT + 
     >   4158*ABS(ETAH) - 4590*ABS(MUH) + 
     >   8330*ABS(XIH)))/25088000 
      Q2(56,55) = -(4131*35**(0.5D0)*MUH)/313600 
      Q2(56,56) = (37179*VT)/5017600 + 
     >   (72171*ABS(ETAH))/12544000 + 
     >   (4131*ABS(MUH))/358400 + (4131*ABS(XIH))/358400 
      Q2(56,57) = 0 
      Q2(56,58) = -(8991*35**(0.5D0)*ETAH)/896000 
      Q2(56,59) = 0 
      Q2(56,60) = (80919*15**(0.5D0)*ETAH)/6272000 
      Q2(56,61) = (2187*3**(0.5D0)*MUH)/1792000 
      Q2(56,62) = (2187*VT)/512000 - (2673*ABS(ETAH))/256000 - 
     >   (6561*ABS(MUH))/1792000 + (1701*ABS(XIH))/256000 
      Q2(56,63) = (2187*15**(0.5D0)*MUH)/224000 
      Q2(56,64) = -(2187*3**(0.5D0)*7**(0.5D0)*(9*VT - 
     >   22*ABS(ETAH) + 14*ABS(MUH) + 
     >   14*ABS(XIH)))/25088000 

      Q2(57,1) = (27*35**(0.5D0)*XIH)/71680 
      Q2(57,2) = 0 
      Q2(57,3) = -(81*7**(0.5D0)*XIH)/358400 
      Q2(57,4) = 0 
      Q2(57,5) = 0 
      Q2(57,6) = 0 
      Q2(57,7) = 0 
      Q2(57,8) = 0 
      Q2(57,9) = -(81*7**(0.5D0)*XIH)/71680 
      Q2(57,10) = 0 
      Q2(57,11) = (243*35**(0.5D0)*XIH)/1792000 
      Q2(57,12) = 0 
      Q2(57,13) = 0 
      Q2(57,14) = 0 
      Q2(57,15) = 0 
      Q2(57,16) = 0 
      Q2(57,17) = (9*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(35*VT - 
     >   70*ABS(ETAH) + 14*ABS(MUH) - 
     >   30*ABS(XIH)))/716800 
      Q2(57,18) = (27*35**(0.5D0)*MUH)/12800 
      Q2(57,19) = (9*3**(0.5D0)*7**(0.5D0)*(42*ABS(ETAH) - 
     >   21*VT + 70*ABS(MUH) + 18*ABS(XIH)))/716800 
      Q2(57,20) = (63*15**(0.5D0)*MUH)/51200 
      Q2(57,21) = (27*35**(0.5D0)*ETAH)/5120 
      Q2(57,22) = 0 
      Q2(57,23) = -(81*7**(0.5D0)*ETAH)/25600 
      Q2(57,24) = 0 
      Q2(57,25) = -(9*3**(0.5D0)*7**(0.5D0)*(105*VT + 
     >   350*ABS(ETAH) + 42*ABS(MUH) - 
     >   90*ABS(XIH)))/716800 
      Q2(57,26) = -(81*7**(0.5D0)*MUH)/12800 
      Q2(57,27) = (27*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(21*VT + 
     >   70*ABS(ETAH) - 70*ABS(MUH) - 
     >   18*ABS(XIH)))/3584000 
      Q2(57,28) = -(189*3**(0.5D0)*MUH)/51200 
      Q2(57,29) = -(63*15**(0.5D0)*ETAH)/10240 
      Q2(57,30) = 0 
      Q2(57,31) = (189*3**(0.5D0)*ETAH)/51200 
      Q2(57,32) = 0 
      Q2(57,33) = (27*7**(0.5D0)*XIH)/1792 
      Q2(57,34) = 0 
      Q2(57,35) = -(81*35**(0.5D0)*XIH)/44800 
      Q2(57,36) = 0 
      Q2(57,37) = 0 
      Q2(57,38) = 0 
      Q2(57,39) = 0 
      Q2(57,40) = 0 
      Q2(57,41) = -(81*35**(0.5D0)*XIH)/8960 
      Q2(57,42) = 0 
      Q2(57,43) = (243*7**(0.5D0)*XIH)/44800 
      Q2(57,44) = 0 
      Q2(57,45) = 0 
      Q2(57,46) = 0 
      Q2(57,47) = 0 
      Q2(57,48) = 0 
      Q2(57,49) = -(27*5**(0.5D0)*(45*VT - 90*ABS(ETAH) + 
     >   18*ABS(MUH) + 70*ABS(XIH)))/716800 
      Q2(57,50) = -(243*15**(0.5D0)*MUH)/89600 
      Q2(57,51) = (729*VT)/716800 - (729*ABS(ETAH))/358400 - 
     >   (243*ABS(MUH))/71680 + (81*ABS(XIH))/51200 
      Q2(57,52) = -(243*35**(0.5D0)*MUH)/358400 
      Q2(57,53) = -(243*15**(0.5D0)*ETAH)/35840 
      Q2(57,54) = 0 
      Q2(57,55) = (729*3**(0.5D0)*ETAH)/179200 
      Q2(57,56) = 0 
      Q2(57,57) = (729*VT)/143360 + (243*ABS(ETAH))/14336 + 
     >   (729*ABS(MUH))/358400 + (81*ABS(XIH))/10240 
      Q2(57,58) = (729*3**(0.5D0)*MUH)/89600 
      Q2(57,59) = -(243*5**(0.5D0)*(9*VT + 30*ABS(ETAH) - 
     >   30*ABS(MUH) + 14*ABS(XIH)))/3584000 
      Q2(57,60) = (729*7**(0.5D0)*MUH)/358400 
      Q2(57,61) = (243*35**(0.5D0)*ETAH)/71680 
      Q2(57,62) = 0 
      Q2(57,63) = -(729*7**(0.5D0)*ETAH)/358400 
      Q2(57,64) = 0 

      Q2(58,1) = 0 
      Q2(58,2) = (153*35**(0.5D0)*XIH)/358400 
      Q2(58,3) = 0 
      Q2(58,4) = -(81*15**(0.5D0)*XIH)/256000 
      Q2(58,5) = 0 
      Q2(58,6) = 0 
      Q2(58,7) = 0 
      Q2(58,8) = 0 
      Q2(58,9) = 0 
      Q2(58,10) = -(459*7**(0.5D0)*XIH)/358400 
      Q2(58,11) = 0 
      Q2(58,12) = (243*3**(0.5D0)*XIH)/256000 
      Q2(58,13) = 0 
      Q2(58,14) = 0 
      Q2(58,15) = 0 
      Q2(58,16) = 0 
      Q2(58,17) = -(99*35**(0.5D0)*MUH)/256000 
      Q2(58,18) = (3*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(595*VT - 
     >   1190*ABS(ETAH) + 462*ABS(MUH) - 
     >   510*ABS(XIH)))/3584000 
      Q2(58,19) = (333*7**(0.5D0)*MUH)/25600 
      Q2(58,20) = (9*5**(0.5D0)*(126*ABS(ETAH) - 63*VT + 
     >   154*ABS(MUH) + 54*ABS(XIH)))/512000 
      Q2(58,21) = 0 
      Q2(58,22) = (153*35**(0.5D0)*ETAH)/25600 
      Q2(58,23) = 0 
      Q2(58,24) = -(567*15**(0.5D0)*ETAH)/128000 
      Q2(58,25) = (297*7**(0.5D0)*MUH)/256000 
      Q2(58,26) = -(3*3**(0.5D0)*7**(0.5D0)*(1785*VT + 
     >   5950*ABS(ETAH) + 1386*ABS(MUH) - 
     >   1530*ABS(XIH)))/3584000 
      Q2(58,27) = -(999*35**(0.5D0)*MUH)/128000 
      Q2(58,28) = (1701*VT)/512000 + (567*ABS(ETAH))/51200 - 
     >   (2079*ABS(MUH))/256000 - (729*ABS(XIH))/256000 
      Q2(58,29) = 0 
      Q2(58,30) = -(357*15**(0.5D0)*ETAH)/51200 
      Q2(58,31) = 0 
      Q2(58,32) = (567*35**(0.5D0)*ETAH)/256000 
      Q2(58,33) = 0 
      Q2(58,34) = (153*7**(0.5D0)*XIH)/8960 
      Q2(58,35) = 0 
      Q2(58,36) = -(81*3**(0.5D0)*XIH)/6400 
      Q2(58,37) = 0 
      Q2(58,38) = 0 
      Q2(58,39) = 0 
      Q2(58,40) = 0 
      Q2(58,41) = 0 
      Q2(58,42) = -(459*35**(0.5D0)*XIH)/44800 
      Q2(58,43) = 0 
      Q2(58,44) = (243*15**(0.5D0)*XIH)/32000 
      Q2(58,45) = 0 
      Q2(58,46) = 0 
      Q2(58,47) = 0 
      Q2(58,48) = 0 
      Q2(58,49) = (891*15**(0.5D0)*MUH)/1792000 
      Q2(58,50) = -(9*5**(0.5D0)*(765*VT - 1530*ABS(ETAH) + 
     >   594*ABS(MUH) + 1190*ABS(XIH)))/3584000 
      Q2(58,51) = -(2997*3**(0.5D0)*MUH)/179200 
      Q2(58,52) = (81*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(9*VT - 
     >   18*ABS(ETAH) - 22*ABS(MUH) + 
     >   14*ABS(XIH)))/3584000 
      Q2(58,53) = 0 
      Q2(58,54) = -(1377*15**(0.5D0)*ETAH)/179200 
      Q2(58,55) = 0 
      Q2(58,56) = (2187*35**(0.5D0)*ETAH)/896000 
      Q2(58,57) = -(2673*3**(0.5D0)*MUH)/1792000 
      Q2(58,58) = (4131*VT)/716800 + (1377*ABS(ETAH))/71680 + 
     >   (8019*ABS(MUH))/1792000 + (459*ABS(XIH))/51200 
      Q2(58,59) = (8991*15**(0.5D0)*MUH)/896000 
      Q2(58,60) = -(243*3**(0.5D0)*7**(0.5D0)*(9*VT + 
     >   30*ABS(ETAH) - 22*ABS(MUH) + 
     >   14*ABS(XIH)))/3584000 
      Q2(58,61) = 0 
      Q2(58,62) = (1377*35**(0.5D0)*ETAH)/358400 
      Q2(58,63) = 0 
      Q2(58,64) = -(729*15**(0.5D0)*ETAH)/256000 

      Q2(59,1) = -(81*7**(0.5D0)*XIH)/358400 
      Q2(59,2) = 0 
      Q2(59,3) = (243*35**(0.5D0)*XIH)/1792000 
      Q2(59,4) = 0 
      Q2(59,5) = 0 
      Q2(59,6) = 0 
      Q2(59,7) = 0 
      Q2(59,8) = 0 
      Q2(59,9) = (243*35**(0.5D0)*XIH)/1792000 
      Q2(59,10) = 0 
      Q2(59,11) = -(729*7**(0.5D0)*XIH)/1792000 
      Q2(59,12) = 0 
      Q2(59,13) = 0 
      Q2(59,14) = 0 
      Q2(59,15) = 0 
      Q2(59,16) = 0 
      Q2(59,17) = (27*3**(0.5D0)*7**(0.5D0)*(14*ABS(ETAH) - 
     >   7*VT + 14*ABS(MUH) + 6*ABS(XIH)))/716800 
      Q2(59,18) = -(81*7**(0.5D0)*MUH)/25600 
      Q2(59,19) = (27*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(21*VT - 
     >   42*ABS(ETAH) + 70*ABS(MUH) - 
     >   18*ABS(XIH)))/3584000 
      Q2(59,20) = (189*3**(0.5D0)*MUH)/51200 
      Q2(59,21) = -(81*7**(0.5D0)*ETAH)/25600 
      Q2(59,22) = 0 
      Q2(59,23) = (243*35**(0.5D0)*ETAH)/128000 
      Q2(59,24) = 0 
      Q2(59,25) = (27*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(21*VT + 
     >   70*ABS(ETAH) - 42*ABS(MUH) - 
     >   18*ABS(XIH)))/3584000 
      Q2(59,26) = (243*35**(0.5D0)*MUH)/128000 
      Q2(59,27) = -(81*3**(0.5D0)*7**(0.5D0)*(21*VT + 
     >   70*ABS(ETAH) + 70*ABS(MUH) - 
     >   18*ABS(XIH)))/3584000 
      Q2(59,28) = -(567*15**(0.5D0)*MUH)/256000 
      Q2(59,29) = (189*3**(0.5D0)*ETAH)/51200 
      Q2(59,30) = 0 
      Q2(59,31) = -(567*15**(0.5D0)*ETAH)/256000 
      Q2(59,32) = 0 
      Q2(59,33) = -(81*35**(0.5D0)*XIH)/44800 
      Q2(59,34) = 0 
      Q2(59,35) = (243*7**(0.5D0)*XIH)/44800 
      Q2(59,36) = 0 
      Q2(59,37) = 0 
      Q2(59,38) = 0 
      Q2(59,39) = 0 
      Q2(59,40) = 0 
      Q2(59,41) = (243*7**(0.5D0)*XIH)/44800 
      Q2(59,42) = 0 
      Q2(59,43) = -(729*35**(0.5D0)*XIH)/224000 
      Q2(59,44) = 0 
      Q2(59,45) = 0 
      Q2(59,46) = 0 
      Q2(59,47) = 0 
      Q2(59,48) = 0 
      Q2(59,49) = (729*VT)/716800 - (729*ABS(ETAH))/358400 - 
     >   (729*ABS(MUH))/358400 + (81*ABS(XIH))/51200 
      Q2(59,50) = (729*3**(0.5D0)*MUH)/179200 
      Q2(59,51) = -(243*5**(0.5D0)*(9*VT - 18*ABS(ETAH) + 
     >   30*ABS(MUH) + 14*ABS(XIH)))/3584000 
      Q2(59,52) = -(729*7**(0.5D0)*MUH)/358400 
      Q2(59,53) = (729*3**(0.5D0)*ETAH)/179200 
      Q2(59,54) = 0 
      Q2(59,55) = -(2187*15**(0.5D0)*ETAH)/896000 
      Q2(59,56) = 0 
      Q2(59,57) = -(243*5**(0.5D0)*(9*VT + 30*ABS(ETAH) - 
     >   18*ABS(MUH) + 14*ABS(XIH)))/3584000 
      Q2(59,58) = -(2187*15**(0.5D0)*MUH)/896000 
      Q2(59,59) = (6561*VT)/3584000 + (2187*ABS(ETAH))/358400 + 
     >   (2187*ABS(MUH))/358400 + (729*ABS(XIH))/256000 
      Q2(59,60) = (2187*35**(0.5D0)*MUH)/1792000 
      Q2(59,61) = -(729*7**(0.5D0)*ETAH)/358400 
      Q2(59,62) = 0 
      Q2(59,63) = (2187*35**(0.5D0)*ETAH)/1792000 
      Q2(59,64) = 0 

      Q2(60,1) = 0 
      Q2(60,2) = -(81*15**(0.5D0)*XIH)/256000 
      Q2(60,3) = 0 
      Q2(60,4) = (2187*35**(0.5D0)*XIH)/12544000 
      Q2(60,5) = 0 
      Q2(60,6) = 0 
      Q2(60,7) = 0 
      Q2(60,8) = 0 
      Q2(60,9) = 0 
      Q2(60,10) = (243*3**(0.5D0)*XIH)/256000 
      Q2(60,11) = 0 
      Q2(60,12) = -(6561*7**(0.5D0)*XIH)/12544000 
      Q2(60,13) = 0 
      Q2(60,14) = 0 
      Q2(60,15) = 0 
      Q2(60,16) = 0 
      Q2(60,17) = -(81*15**(0.5D0)*MUH)/256000 
      Q2(60,18) = (81*5**(0.5D0)*(14*ABS(ETAH) - 7*VT + 
     >   6*ABS(MUH) + 6*ABS(XIH)))/512000 
      Q2(60,19) = -(81*3**(0.5D0)*MUH)/6400 
      Q2(60,20) = (81*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(63*VT - 
     >   126*ABS(ETAH) + 98*ABS(MUH) - 
     >   54*ABS(XIH)))/25088000 
      Q2(60,21) = 0 
      Q2(60,22) = -(567*15**(0.5D0)*ETAH)/128000 
      Q2(60,23) = 0 
      Q2(60,24) = (2187*35**(0.5D0)*ETAH)/896000 
      Q2(60,25) = (243*3**(0.5D0)*MUH)/256000 
      Q2(60,26) = (1701*VT)/512000 + (567*ABS(ETAH))/51200 - 
     >   (729*ABS(MUH))/256000 - (729*ABS(XIH))/256000 
      Q2(60,27) = (243*15**(0.5D0)*MUH)/32000 
      Q2(60,28) = -(243*3**(0.5D0)*7**(0.5D0)*(63*VT + 
     >   210*ABS(ETAH) + 98*ABS(MUH) - 
     >   54*ABS(XIH)))/25088000 
      Q2(60,29) = 0 
      Q2(60,30) = (567*35**(0.5D0)*ETAH)/256000 
      Q2(60,31) = 0 
      Q2(60,32) = -(729*15**(0.5D0)*ETAH)/256000 
      Q2(60,33) = 0 
      Q2(60,34) = -(81*3**(0.5D0)*XIH)/6400 
      Q2(60,35) = 0 
      Q2(60,36) = (2187*7**(0.5D0)*XIH)/313600 
      Q2(60,37) = 0 
      Q2(60,38) = 0 
      Q2(60,39) = 0 
      Q2(60,40) = 0 
      Q2(60,41) = 0 
      Q2(60,42) = (243*15**(0.5D0)*XIH)/32000 
      Q2(60,43) = 0 
      Q2(60,44) = -(6561*35**(0.5D0)*XIH)/1568000 
      Q2(60,45) = 0 
      Q2(60,46) = 0 
      Q2(60,47) = 0 
      Q2(60,48) = 0 
      Q2(60,49) = (2187*35**(0.5D0)*MUH)/12544000 
      Q2(60,50) = (81*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(63*VT - 
     >   126*ABS(ETAH) - 54*ABS(MUH) + 
     >   98*ABS(XIH)))/25088000 
      Q2(60,51) = (2187*7**(0.5D0)*MUH)/313600 
      Q2(60,52) = -(2187*5**(0.5D0)*(9*VT - 18*ABS(ETAH) + 
     >   14*ABS(MUH) + 14*ABS(XIH)))/25088000 
      Q2(60,53) = 0 
      Q2(60,54) = (2187*35**(0.5D0)*ETAH)/896000 
      Q2(60,55) = 0 
      Q2(60,56) = -(19683*15**(0.5D0)*ETAH)/6272000 
      Q2(60,57) = -(6561*7**(0.5D0)*MUH)/12544000 
      Q2(60,58) = -(243*3**(0.5D0)*7**(0.5D0)*(63*VT + 
     >   210*ABS(ETAH) - 54*ABS(MUH) + 
     >   98*ABS(XIH)))/25088000 
      Q2(60,59) = -(6561*35**(0.5D0)*MUH)/1568000 
      Q2(60,60) = (59049*VT)/25088000 + 
     >   (19683*ABS(ETAH))/2508800 + 
     >   (6561*ABS(MUH))/1792000 + (6561*ABS(XIH))/1792000 
      Q2(60,61) = 0 
      Q2(60,62) = -(729*15**(0.5D0)*ETAH)/256000 
      Q2(60,63) = 0 
      Q2(60,64) = (19683*35**(0.5D0)*ETAH)/12544000 

      Q2(61,1) = 0 
      Q2(61,2) = 0 
      Q2(61,3) = 0 
      Q2(61,4) = 0 
      Q2(61,5) = (27*3**(0.5D0)*XIH)/10240 
      Q2(61,6) = 0 
      Q2(61,7) = -(81*15**(0.5D0)*XIH)/256000 
      Q2(61,8) = 0 
      Q2(61,9) = 0 
      Q2(61,10) = 0 
      Q2(61,11) = 0 
      Q2(61,12) = 0 
      Q2(61,13) = -(729*7**(0.5D0)*XIH)/501760 
      Q2(61,14) = 0 
      Q2(61,15) = (2187*35**(0.5D0)*XIH)/12544000 
      Q2(61,16) = 0 
      Q2(61,17) = (27*3**(0.5D0)*ETAH)/10240 
      Q2(61,18) = 0 
      Q2(61,19) = -(81*15**(0.5D0)*ETAH)/256000 
      Q2(61,20) = 0 
      Q2(61,21) = (189*VT)/20480 - (81*ABS(ETAH))/10240 + 
     >   (189*ABS(MUH))/51200 - (81*ABS(XIH))/10240 
      Q2(61,22) = (189*3**(0.5D0)*MUH)/12800 
      Q2(61,23) = (27*5**(0.5D0)*(18*ABS(ETAH) - 21*VT + 
     >   70*ABS(MUH) + 18*ABS(XIH)))/512000 
      Q2(61,24) = (189*7**(0.5D0)*MUH)/51200 
      Q2(61,25) = (27*15**(0.5D0)*ETAH)/1280 
      Q2(61,26) = 0 
      Q2(61,27) = -(81*3**(0.5D0)*ETAH)/6400 
      Q2(61,28) = 0 
      Q2(61,29) = -(27*3**(0.5D0)*7**(0.5D0)*(315*VT + 
     >   490*ABS(ETAH) + 126*ABS(MUH) - 
     >   270*ABS(XIH)))/5017600 
      Q2(61,30) = -(729*7**(0.5D0)*MUH)/89600 
      Q2(61,31) = (81*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(63*VT + 
     >   98*ABS(ETAH) - 210*ABS(MUH) - 
     >   54*ABS(XIH)))/25088000 
      Q2(61,32) = -(243*3**(0.5D0)*MUH)/51200 
      Q2(61,33) = 0 
      Q2(61,34) = 0 
      Q2(61,35) = 0 
      Q2(61,36) = 0 
      Q2(61,37) = (27*15**(0.5D0)*XIH)/1280 
      Q2(61,38) = 0 
      Q2(61,39) = -(81*3**(0.5D0)*XIH)/6400 
      Q2(61,40) = 0 
      Q2(61,41) = 0 
      Q2(61,42) = 0 
      Q2(61,43) = 0 
      Q2(61,44) = 0 
      Q2(61,45) = -(729*35**(0.5D0)*XIH)/62720 
      Q2(61,46) = 0 
      Q2(61,47) = (2187*7**(0.5D0)*XIH)/313600 
      Q2(61,48) = 0 
      Q2(61,49) = -(729*7**(0.5D0)*ETAH)/501760 
      Q2(61,50) = 0 
      Q2(61,51) = (2187*35**(0.5D0)*ETAH)/12544000 
      Q2(61,52) = 0 
      Q2(61,53) = -(27*3**(0.5D0)*7**(0.5D0)*(315*VT - 
     >   270*ABS(ETAH) + 126*ABS(MUH) + 
     >   490*ABS(XIH)))/5017600 
      Q2(61,54) = -(729*7**(0.5D0)*MUH)/89600 
      Q2(61,55) = (81*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(63*VT - 
     >   54*ABS(ETAH) - 210*ABS(MUH) + 
     >   98*ABS(XIH)))/25088000 
      Q2(61,56) = -(243*3**(0.5D0)*MUH)/51200 
      Q2(61,57) = -(729*35**(0.5D0)*ETAH)/62720 
      Q2(61,58) = 0 
      Q2(61,59) = (2187*7**(0.5D0)*ETAH)/313600 
      Q2(61,60) = 0 
      Q2(61,61) = (6561*VT)/1003520 + (729*ABS(ETAH))/71680 + 
     >   (6561*ABS(MUH))/2508800 + (729*ABS(XIH))/71680 
      Q2(61,62) = (6561*3**(0.5D0)*MUH)/627200 
      Q2(61,63) = -(2187*5**(0.5D0)*(9*VT + 14*ABS(ETAH) - 
     >   30*ABS(MUH) + 14*ABS(XIH)))/25088000 
      Q2(61,64) = (6561*7**(0.5D0)*MUH)/2508800 

      Q2(62,1) = 0 
      Q2(62,2) = 0 
      Q2(62,3) = 0 
      Q2(62,4) = 0 
      Q2(62,5) = 0 
      Q2(62,6) = (153*3**(0.5D0)*XIH)/51200 
      Q2(62,7) = 0 
      Q2(62,8) = -(243*7**(0.5D0)*XIH)/256000 
      Q2(62,9) = 0 
      Q2(62,10) = 0 
      Q2(62,11) = 0 
      Q2(62,12) = 0 
      Q2(62,13) = 0 
      Q2(62,14) = -(4131*7**(0.5D0)*XIH)/2508800 
      Q2(62,15) = 0 
      Q2(62,16) = (2187*3**(0.5D0)*XIH)/1792000 
      Q2(62,17) = 0 
      Q2(62,18) = (153*3**(0.5D0)*ETAH)/51200 
      Q2(62,19) = 0 
      Q2(62,20) = -(243*7**(0.5D0)*ETAH)/256000 
      Q2(62,21) = -(693*3**(0.5D0)*MUH)/256000 
      Q2(62,22) = (1071*VT)/102400 - (459*ABS(ETAH))/51200 + 
     >   (2079*ABS(MUH))/256000 - (459*ABS(XIH))/51200 
      Q2(62,23) = (2331*15**(0.5D0)*MUH)/128000 
      Q2(62,24) = (9*3**(0.5D0)*7**(0.5D0)*(54*ABS(ETAH) - 
     >   63*VT + 154*ABS(MUH) + 54*ABS(XIH)))/512000 
      Q2(62,25) = 0 
      Q2(62,26) = (153*15**(0.5D0)*ETAH)/6400 
      Q2(62,27) = 0 
      Q2(62,28) = -(243*35**(0.5D0)*ETAH)/32000 
      Q2(62,29) = (2673*7**(0.5D0)*MUH)/1792000 
      Q2(62,30) = -(9*3**(0.5D0)*7**(0.5D0)*(5355*VT + 
     >   8330*ABS(ETAH) + 4158*ABS(MUH) - 
     >   4590*ABS(XIH)))/25088000 
      Q2(62,31) = -(8991*35**(0.5D0)*MUH)/896000 
      Q2(62,32) = (2187*VT)/512000 + (1701*ABS(ETAH))/256000 - 
     >   (2673*ABS(MUH))/256000 - (6561*ABS(XIH))/1792000 
      Q2(62,33) = 0 
      Q2(62,34) = 0 
      Q2(62,35) = 0 
      Q2(62,36) = 0 
      Q2(62,37) = 0 
      Q2(62,38) = (153*15**(0.5D0)*XIH)/6400 
      Q2(62,39) = 0 
      Q2(62,40) = -(243*35**(0.5D0)*XIH)/32000 
      Q2(62,41) = 0 
      Q2(62,42) = 0 
      Q2(62,43) = 0 
      Q2(62,44) = 0 
      Q2(62,45) = 0 
      Q2(62,46) = -(4131*35**(0.5D0)*XIH)/313600 
      Q2(62,47) = 0 
      Q2(62,48) = (2187*15**(0.5D0)*XIH)/224000 
      Q2(62,49) = 0 
      Q2(62,50) = -(4131*7**(0.5D0)*ETAH)/2508800 
      Q2(62,51) = 0 
      Q2(62,52) = (2187*3**(0.5D0)*ETAH)/1792000 
      Q2(62,53) = (2673*7**(0.5D0)*MUH)/1792000 
      Q2(62,54) = -(9*3**(0.5D0)*7**(0.5D0)*(5355*VT - 
     >   4590*ABS(ETAH) + 4158*ABS(MUH) + 
     >   8330*ABS(XIH)))/25088000 
      Q2(62,55) = -(8991*35**(0.5D0)*MUH)/896000 
      Q2(62,56) = (2187*VT)/512000 - (6561*ABS(ETAH))/1792000 - 
     >   (2673*ABS(MUH))/256000 + (1701*ABS(XIH))/256000 
      Q2(62,57) = 0 
      Q2(62,58) = -(4131*35**(0.5D0)*ETAH)/313600 
      Q2(62,59) = 0 
      Q2(62,60) = (2187*15**(0.5D0)*ETAH)/224000 
      Q2(62,61) = -(24057*3**(0.5D0)*MUH)/12544000 
      Q2(62,62) = (37179*VT)/5017600 + (4131*ABS(ETAH))/358400 + 
     >   (72171*ABS(MUH))/12544000 + (4131*ABS(XIH))/358400 
      Q2(62,63) = (80919*15**(0.5D0)*MUH)/6272000 
      Q2(62,64) = -(2187*3**(0.5D0)*7**(0.5D0)*(9*VT + 
     >   14*ABS(ETAH) - 22*ABS(MUH) + 
     >   14*ABS(XIH)))/25088000 

      Q2(63,1) = 0 
      Q2(63,2) = 0 
      Q2(63,3) = 0 
      Q2(63,4) = 0 
      Q2(63,5) = -(81*15**(0.5D0)*XIH)/256000 
      Q2(63,6) = 0 
      Q2(63,7) = (243*3**(0.5D0)*XIH)/256000 
      Q2(63,8) = 0 
      Q2(63,9) = 0 
      Q2(63,10) = 0 
      Q2(63,11) = 0 
      Q2(63,12) = 0 
      Q2(63,13) = (2187*35**(0.5D0)*XIH)/12544000 
      Q2(63,14) = 0 
      Q2(63,15) = -(6561*7**(0.5D0)*XIH)/12544000 
      Q2(63,16) = 0 
      Q2(63,17) = -(81*15**(0.5D0)*ETAH)/256000 
      Q2(63,18) = 0 
      Q2(63,19) = (243*3**(0.5D0)*ETAH)/256000 
      Q2(63,20) = 0 
      Q2(63,21) = (81*5**(0.5D0)*(6*ABS(ETAH) - 7*VT + 
     >   14*ABS(MUH) + 6*ABS(XIH)))/512000 
      Q2(63,22) = -(567*15**(0.5D0)*MUH)/128000 
      Q2(63,23) = (1701*VT)/512000 - (729*ABS(ETAH))/256000 + 
     >   (567*ABS(MUH))/51200 - (729*ABS(XIH))/256000 
      Q2(63,24) = (567*35**(0.5D0)*MUH)/256000 
      Q2(63,25) = -(81*3**(0.5D0)*ETAH)/6400 
      Q2(63,26) = 0 
      Q2(63,27) = (243*15**(0.5D0)*ETAH)/32000 
      Q2(63,28) = 0 
      Q2(63,29) = (81*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(63*VT + 
     >   98*ABS(ETAH) - 126*ABS(MUH) - 
     >   54*ABS(XIH)))/25088000 
      Q2(63,30) = (2187*35**(0.5D0)*MUH)/896000 
      Q2(63,31) = -(243*3**(0.5D0)*7**(0.5D0)*(63*VT + 
     >   98*ABS(ETAH) + 210*ABS(MUH) - 
     >   54*ABS(XIH)))/25088000 
      Q2(63,32) = -(729*15**(0.5D0)*MUH)/256000 
      Q2(63,33) = 0 
      Q2(63,34) = 0 
      Q2(63,35) = 0 
      Q2(63,36) = 0 
      Q2(63,37) = -(81*3**(0.5D0)*XIH)/6400 
      Q2(63,38) = 0 
      Q2(63,39) = (243*15**(0.5D0)*XIH)/32000 
      Q2(63,40) = 0 
      Q2(63,41) = 0 
      Q2(63,42) = 0 
      Q2(63,43) = 0 
      Q2(63,44) = 0 
      Q2(63,45) = (2187*7**(0.5D0)*XIH)/313600 
      Q2(63,46) = 0 
      Q2(63,47) = -(6561*35**(0.5D0)*XIH)/1568000 
      Q2(63,48) = 0 
      Q2(63,49) = (2187*35**(0.5D0)*ETAH)/12544000 
      Q2(63,50) = 0 
      Q2(63,51) = -(6561*7**(0.5D0)*ETAH)/12544000 
      Q2(63,52) = 0 
      Q2(63,53) = (81*3**(0.5D0)*5**(0.5D0)*7**(0.5D0)*(63*VT - 
     >   54*ABS(ETAH) - 126*ABS(MUH) + 
     >   98*ABS(XIH)))/25088000 
      Q2(63,54) = (2187*35**(0.5D0)*MUH)/896000 
      Q2(63,55) = -(243*3**(0.5D0)*7**(0.5D0)*(63*VT - 
     >   54*ABS(ETAH) + 210*ABS(MUH) + 
     >   98*ABS(XIH)))/25088000 
      Q2(63,56) = -(729*15**(0.5D0)*MUH)/256000 
      Q2(63,57) = (2187*7**(0.5D0)*ETAH)/313600 
      Q2(63,58) = 0 
      Q2(63,59) = -(6561*35**(0.5D0)*ETAH)/1568000 
      Q2(63,60) = 0 
      Q2(63,61) = -(2187*5**(0.5D0)*(9*VT + 14*ABS(ETAH) - 
     >   18*ABS(MUH) + 14*ABS(XIH)))/25088000 
      Q2(63,62) = -(19683*15**(0.5D0)*MUH)/6272000 
      Q2(63,63) = (59049*VT)/25088000 + 
     >   (6561*ABS(ETAH))/1792000 + 
     >   (19683*ABS(MUH))/2508800 + (6561*ABS(XIH))/1792000 
      Q2(63,64) = (19683*35**(0.5D0)*MUH)/12544000 

      Q2(64,1) = 0 
      Q2(64,2) = 0 
      Q2(64,3) = 0 
      Q2(64,4) = 0 
      Q2(64,5) = 0 
      Q2(64,6) = -(243*7**(0.5D0)*XIH)/256000 
      Q2(64,7) = 0 
      Q2(64,8) = (2187*3**(0.5D0)*XIH)/1792000 
      Q2(64,9) = 0 
      Q2(64,10) = 0 
      Q2(64,11) = 0 
      Q2(64,12) = 0 
      Q2(64,13) = 0 
      Q2(64,14) = (2187*3**(0.5D0)*XIH)/1792000 
      Q2(64,15) = 0 
      Q2(64,16) = -(59049*7**(0.5D0)*XIH)/87808000 
      Q2(64,17) = 0 
      Q2(64,18) = -(243*7**(0.5D0)*ETAH)/256000 
      Q2(64,19) = 0 
      Q2(64,20) = (2187*3**(0.5D0)*ETAH)/1792000 
      Q2(64,21) = -(243*7**(0.5D0)*MUH)/256000 
      Q2(64,22) = (81*3**(0.5D0)*7**(0.5D0)*(6*ABS(ETAH) - 
     >   7*VT + 6*ABS(MUH) + 6*ABS(XIH)))/512000 
      Q2(64,23) = -(243*35**(0.5D0)*MUH)/32000 
      Q2(64,24) = (2187*VT)/512000 - (6561*ABS(ETAH))/1792000 + 
     >   (1701*ABS(MUH))/256000 - (6561*ABS(XIH))/1792000 
      Q2(64,25) = 0 
      Q2(64,26) = -(243*35**(0.5D0)*ETAH)/32000 
      Q2(64,27) = 0 
      Q2(64,28) = (2187*15**(0.5D0)*ETAH)/224000 
      Q2(64,29) = (2187*3**(0.5D0)*MUH)/1792000 
      Q2(64,30) = (2187*VT)/512000 + (1701*ABS(ETAH))/256000 - 
     >   (6561*ABS(MUH))/1792000 - (6561*ABS(XIH))/1792000 
      Q2(64,31) = (2187*15**(0.5D0)*MUH)/224000 
      Q2(64,32) = -(2187*3**(0.5D0)*7**(0.5D0)*(63*VT + 
     >   98*ABS(ETAH) + 98*ABS(MUH) - 
     >   54*ABS(XIH)))/175616000 
      Q2(64,33) = 0 
      Q2(64,34) = 0 
      Q2(64,35) = 0 
      Q2(64,36) = 0 
      Q2(64,37) = 0 
      Q2(64,38) = -(243*35**(0.5D0)*XIH)/32000 
      Q2(64,39) = 0 
      Q2(64,40) = (2187*15**(0.5D0)*XIH)/224000 
      Q2(64,41) = 0 
      Q2(64,42) = 0 
      Q2(64,43) = 0 
      Q2(64,44) = 0 
      Q2(64,45) = 0 
      Q2(64,46) = (2187*15**(0.5D0)*XIH)/224000 
      Q2(64,47) = 0 
      Q2(64,48) = -(59049*35**(0.5D0)*XIH)/10976000 
      Q2(64,49) = 0 
      Q2(64,50) = (2187*3**(0.5D0)*ETAH)/1792000 
      Q2(64,51) = 0 
      Q2(64,52) = -(59049*7**(0.5D0)*ETAH)/87808000 
      Q2(64,53) = (2187*3**(0.5D0)*MUH)/1792000 
      Q2(64,54) = (2187*VT)/512000 - (6561*ABS(ETAH))/1792000 - 
     >   (6561*ABS(MUH))/1792000 + (1701*ABS(XIH))/256000 
      Q2(64,55) = (2187*15**(0.5D0)*MUH)/224000 
      Q2(64,56) = -(2187*3**(0.5D0)*7**(0.5D0)*(63*VT - 
     >   54*ABS(ETAH) + 98*ABS(MUH) + 
     >   98*ABS(XIH)))/175616000 
      Q2(64,57) = 0 
      Q2(64,58) = (2187*15**(0.5D0)*ETAH)/224000 
      Q2(64,59) = 0 
      Q2(64,60) = -(59049*35**(0.5D0)*ETAH)/10976000 
      Q2(64,61) = -(59049*7**(0.5D0)*MUH)/87808000 
      Q2(64,62) = -(2187*3**(0.5D0)*7**(0.5D0)*(63*VT + 
     >   98*ABS(ETAH) - 54*ABS(MUH) + 
     >   98*ABS(XIH)))/175616000 
      Q2(64,63) = -(59049*35**(0.5D0)*MUH)/10976000 
      Q2(64,64) = (531441*VT)/175616000 + 
     >   (59049*ABS(ETAH))/12544000 + 
     >   (59049*ABS(MUH))/12544000 + 
     >   (59049*ABS(XIH))/12544000 

      Q2(1,65) = (((125000*Q(01) + 9000*Q(11) + 9000*Q(35) + 
     >   9000*Q(41) - 15000*5**(0.5D0)*Q(03) - 15000*5**(0.5D0)*
     >   Q(09) - 15000*5**(0.5D0)*Q(33) - 1080*5**(0.5D0)*Q(43)))/
     >   4096000)*VOL(I2,J2,J,I,K2) 
      Q2(2,65) = ((3**(0.5D0)*(2125000*3**(0.5D0)*Q(02) - 
     >   675000*7**(0.5D0)*Q(04) - 255000*15**(0.5D0)*Q(10) + 
     >   153000*3**(0.5D0)*Q(42) + 81000*35**(0.5D0)*Q(12) - 
     >   255000*15**(0.5D0)*Q(34) - 48600*7**(0.5D0)*Q(44) + 
     >   81000*35**(0.5D0)*Q(36)))/184320000)*VOL(I2,J2,J,I,K2) 
      Q2(3,65) = (-(3*5**(0.5D0)*(25000*Q(01) + 9000*Q(11) + 
     >   9000*Q(35) + 1800*Q(41) - 15000*5**(0.5D0)*Q(03) - 
     >   3000*5**(0.5D0)*Q(09) - 3000*5**(0.5D0)*Q(33) - 
     >   1080*5**(0.5D0)*Q(43)))/20480000)*VOL(I2,J2,J,I,K2) 
      Q2(4,65) = (-(3*7**(0.5D0)*(1225000*3**(0.5D0)*Q(02) - 
     >   675000*7**(0.5D0)*Q(04) - 147000*15**(0.5D0)*Q(10) + 
     >   88200*3**(0.5D0)*Q(42) + 81000*35**(0.5D0)*Q(12) - 
     >   147000*15**(0.5D0)*Q(34) - 48600*7**(0.5D0)*Q(44) + 
     >   81000*35**(0.5D0)*Q(36)))/1003520000)*VOL(I2,J2,J,I,K2) 
      Q2(5,65) = ((3**(0.5D0)*(2125000*3**(0.5D0)*Q(05) - 
     >   675000*7**(0.5D0)*Q(13) - 255000*15**(0.5D0)*Q(07) + 
     >   153000*3**(0.5D0)*Q(39) + 81000*35**(0.5D0)*Q(15) - 
     >   255000*15**(0.5D0)*Q(37) - 48600*7**(0.5D0)*Q(47) + 
     >   81000*35**(0.5D0)*Q(45)))/184320000)*VOL(I2,J2,J,I,K2) 
      Q2(6,65) = (((108375000*Q(06) + 25515000*Q(16) - 
     >   11475000*21**(0.5D0)*Q(08) - 11475000*21**(0.5D0)*Q(14) - 
     >   13005000*5**(0.5D0)*Q(38) - 3061800*5**(0.5D0)*Q(48) + 
     >   1377000*105**(0.5D0)*Q(40) + 1377000*105**(0.5D0)*Q(46)))/
     >   CONST0)*VOL(I2,J2,J,I,K2) 
      Q2(7,65) = (-(15**(0.5D0)*(425000*3**(0.5D0)*Q(05) - 
     >   135000*7**(0.5D0)*Q(13) - 255000*15**(0.5D0)*Q(07) + 
     >   153000*3**(0.5D0)*Q(39) + 81000*35**(0.5D0)*Q(15) - 
     >   51000*15**(0.5D0)*Q(37) - 48600*7**(0.5D0)*Q(47) + 
     >   16200*35**(0.5D0)*Q(45)))/307200000)*VOL(I2,J2,J,I,K2) 
      Q2(8,65) = (-(21**(0.5D0)*(62475000*Q(06) + 25515000*Q(16) - 
     >   11475000*21**(0.5D0)*Q(08) - 6615000*21**(0.5D0)*Q(14) - 
     >   7497000*5**(0.5D0)*Q(38) - 3061800*5**(0.5D0)*Q(48) + 
     >   1377000*105**(0.5D0)*Q(40) + 793800*105**(0.5D0)*Q(46)))/
     >   CONST1)*VOL(I2,J2,J,I,K2) 
      Q2(9,65) = (-(3*5**(0.5D0)*(25000*Q(01) + 9000*Q(11) + 
     >   1800*Q(35) + 9000*Q(41) - 3000*5**(0.5D0)*Q(03) - 
     >   15000*5**(0.5D0)*Q(09) - 3000*5**(0.5D0)*Q(33) - 
     >   1080*5**(0.5D0)*Q(43)))/20480000)*VOL(I2,J2,J,I,K2) 
      Q2(10,65) = (-(15**(0.5D0)*(425000*3**(0.5D0)*Q(02) - 
     >   135000*7**(0.5D0)*Q(04) - 255000*15**(0.5D0)*Q(10) + 
     >   153000*3**(0.5D0)*Q(42) + 81000*35**(0.5D0)*Q(12) - 
     >   51000*15**(0.5D0)*Q(34) - 48600*7**(0.5D0)*Q(44) + 
     >   16200*35**(0.5D0)*Q(36)))/307200000)*VOL(I2,J2,J,I,K2) 
      Q2(11,65) = ((9*(5000*Q(01) + 9000*Q(11) + 1800*Q(35) + 
     >   1800*Q(41) - 3000*5**(0.5D0)*Q(03) - 3000*5**(0.5D0)*
     >   Q(09) - 600*5**(0.5D0)*Q(33) - 1080*5**(0.5D0)*Q(43)))/
     >   20480000)*VOL(I2,J2,J,I,K2) 
      Q2(12,65) = ((9*35**(0.5D0)*(245000*3**(0.5D0)*Q(02) - 
     >   135000*7**(0.5D0)*Q(04) - 147000*15**(0.5D0)*Q(10) + 
     >   88200*3**(0.5D0)*Q(42) + 81000*35**(0.5D0)*Q(12) - 
     >   29400*15**(0.5D0)*Q(34) - 48600*7**(0.5D0)*Q(44) + 
     >   16200*35**(0.5D0)*Q(36)))/CONST2)*VOL(I2,J2,J,I,K2) 
      Q2(13,65) = (-(3*7**(0.5D0)*(1225000*3**(0.5D0)*Q(05) - 
     >   675000*7**(0.5D0)*Q(13) - 147000*15**(0.5D0)*Q(07) + 
     >   88200*3**(0.5D0)*Q(39) + 81000*35**(0.5D0)*Q(15) - 
     >   147000*15**(0.5D0)*Q(37) - 48600*7**(0.5D0)*Q(47) + 
     >   81000*35**(0.5D0)*Q(45)))/1003520000)*VOL(I2,J2,J,I,K2) 
      Q2(14,65) = (-(21**(0.5D0)*(62475000*Q(06) + 25515000*Q(16) - 
     >   6615000*21**(0.5D0)*Q(08) - 11475000*21**(0.5D0)*Q(14) - 
     >   7497000*5**(0.5D0)*Q(38) - 3061800*5**(0.5D0)*Q(48) + 
     >   793800*105**(0.5D0)*Q(40) + 1377000*105**(0.5D0)*Q(46)))/
     >   CONST1)*VOL(I2,J2,J,I,K2) 
      Q2(15,65) = ((9*35**(0.5D0)*(245000*3**(0.5D0)*Q(05) - 
     >   135000*7**(0.5D0)*Q(13) - 147000*15**(0.5D0)*Q(07) + 
     >   88200*3**(0.5D0)*Q(39) + 81000*35**(0.5D0)*Q(15) - 
     >   29400*15**(0.5D0)*Q(37) - 48600*7**(0.5D0)*Q(47) + 
     >   16200*35**(0.5D0)*Q(45)))/CONST2)*VOL(I2,J2,J,I,K2) 
      Q2(16,65) = ((9*(36015000*Q(06) + 25515000*Q(16) - 
     >   6615000*21**(0.5D0)*Q(08) - 6615000*21**(0.5D0)*Q(14) - 
     >   4321800*5**(0.5D0)*Q(38) - 3061800*5**(0.5D0)*Q(48) + 
     >   793800*105**(0.5D0)*Q(40) + 793800*105**(0.5D0)*Q(46)))/
     >   CONST3)*VOL(I2,J2,J,I,K2) 
      Q2(17,65) = ((3**(0.5D0)*(2125000*3**(0.5D0)*Q(17) + 
     >   153000*3**(0.5D0)*Q(27) - 255000*15**(0.5D0)*Q(19) - 
     >   255000*15**(0.5D0)*Q(25) - 675000*7**(0.5D0)*Q(49) - 
     >   48600*7**(0.5D0)*Q(59) + 81000*35**(0.5D0)*Q(51) + 
     >   81000*35**(0.5D0)*Q(57)))/184320000)*VOL(I2,J2,J,I,K2) 
      Q2(18,65) = (((108375000*Q(18) + 25515000*Q(52) - 
     >   13005000*5**(0.5D0)*Q(26) - 11475000*21**(0.5D0)*Q(20) - 
     >   3061800*5**(0.5D0)*Q(60) - 11475000*21**(0.5D0)*Q(50) + 
     >   1377000*105**(0.5D0)*Q(28) + 1377000*105**(0.5D0)*Q(58)))/
     >   CONST0)*VOL(I2,J2,J,I,K2) 
      Q2(19,65) = (-(15**(0.5D0)*(425000*3**(0.5D0)*Q(17) + 
     >   153000*3**(0.5D0)*Q(27) - 255000*15**(0.5D0)*Q(19) - 
     >   51000*15**(0.5D0)*Q(25) - 135000*7**(0.5D0)*Q(49) - 
     >   48600*7**(0.5D0)*Q(59) + 81000*35**(0.5D0)*Q(51) + 
     >   16200*35**(0.5D0)*Q(57)))/307200000)*VOL(I2,J2,J,I,K2) 
      Q2(20,65) = (-(21**(0.5D0)*(62475000*Q(18) + 25515000*Q(52) - 
     >   7497000*5**(0.5D0)*Q(26) - 11475000*21**(0.5D0)*Q(20) - 
     >   3061800*5**(0.5D0)*Q(60) - 6615000*21**(0.5D0)*Q(50) + 
     >   1377000*105**(0.5D0)*Q(28) + 793800*105**(0.5D0)*Q(58)))/
     >   CONST1)*VOL(I2,J2,J,I,K2) 
      Q2(21,65) = (((108375000*Q(21) + 25515000*Q(61) - 
     >   13005000*5**(0.5D0)*Q(23) - 11475000*21**(0.5D0)*Q(29) - 
     >   3061800*5**(0.5D0)*Q(63) - 11475000*21**(0.5D0)*Q(53) + 
     >   1377000*105**(0.5D0)*Q(31) + 1377000*105**(0.5D0)*Q(55)))/
     >   CONST0)*VOL(I2,J2,J,I,K2) 
      Q2(22,65) = (((614125*Q(22) + 144585*Q(32) + 144585*Q(56) + 
     >   144585*Q(62) - 65025*21**(0.5D0)*Q(24) - 65025*21**(0.5D0)*
     >   Q(30) - 65025*21**(0.5D0)*Q(54) - 15309*21**(0.5D0)*
     >   Q(64)))/13824000)*VOL(I2,J2,J,I,K2) 
      Q2(23,65) = (-(5**(0.5D0)*(21675000*Q(21) + 5103000*Q(61) - 
     >   13005000*5**(0.5D0)*Q(23) - 2295000*21**(0.5D0)*Q(29) - 
     >   3061800*5**(0.5D0)*Q(63) - 2295000*21**(0.5D0)*Q(53) + 
     >   1377000*105**(0.5D0)*Q(31) + 1377000*105**(0.5D0)*Q(55)))/
     >   CONST4)*VOL(I2,J2,J,I,K2) 
      Q2(24,65) = (-(7**(0.5D0)*(1062075000*3**(0.5D0)*Q(22) - 
     >   585225000*7**(0.5D0)*Q(24) + 433755000*3**(0.5D0)*Q(32) - 
     >   337365000*7**(0.5D0)*Q(30) + 433755000*3**(0.5D0)*Q(56) - 
     >   337365000*7**(0.5D0)*Q(54) + 250047000*3**(0.5D0)*Q(62) - 
     >   137781000*7**(0.5D0)*Q(64)))/CONST5)*VOL(I2,J2,J,I,K2) 
      Q2(25,65) = (-(15**(0.5D0)*(425000*3**(0.5D0)*Q(17) + 
     >   153000*3**(0.5D0)*Q(27) - 51000*15**(0.5D0)*Q(19) - 
     >   255000*15**(0.5D0)*Q(25) - 135000*7**(0.5D0)*Q(49) - 
     >   48600*7**(0.5D0)*Q(59) + 16200*35**(0.5D0)*Q(51) + 
     >   81000*35**(0.5D0)*Q(57)))/307200000)*VOL(I2,J2,J,I,K2) 
      Q2(26,65) = (-(5**(0.5D0)*(21675000*Q(18) + 5103000*Q(52) - 
     >   13005000*5**(0.5D0)*Q(26) - 2295000*21**(0.5D0)*Q(20) - 
     >   3061800*5**(0.5D0)*Q(60) - 2295000*21**(0.5D0)*Q(50) + 
     >   1377000*105**(0.5D0)*Q(28) + 1377000*105**(0.5D0)*Q(58)))/
     >   CONST4)*VOL(I2,J2,J,I,K2) 
      Q2(27,65) = ((3**(0.5D0)*(85000*3**(0.5D0)*Q(17) + 
     >   153000*3**(0.5D0)*Q(27) - 51000*15**(0.5D0)*Q(19) - 
     >   51000*15**(0.5D0)*Q(25) - 27000*7**(0.5D0)*Q(49) - 
     >   48600*7**(0.5D0)*Q(59) + 16200*35**(0.5D0)*Q(51) + 
     >   16200*35**(0.5D0)*Q(57)))/102400000)*VOL(I2,J2,J,I,K2) 
      Q2(28,65) = ((105**(0.5D0)*(12495000*Q(18) + 5103000*Q(52) - 
     >   7497000*5**(0.5D0)*Q(26) - 2295000*21**(0.5D0)*Q(20) - 
     >   3061800*5**(0.5D0)*Q(60) - 1323000*21**(0.5D0)*Q(50) + 
     >   1377000*105**(0.5D0)*Q(28) + 793800*105**(0.5D0)*Q(58)))/
     >   CONST6)*VOL(I2,J2,J,I,K2) 
      Q2(29,65) = (-(21**(0.5D0)*(62475000*Q(21) + 25515000*Q(61) - 
     >   7497000*5**(0.5D0)*Q(23) - 11475000*21**(0.5D0)*Q(29) - 
     >   3061800*5**(0.5D0)*Q(63) - 6615000*21**(0.5D0)*Q(53) + 
     >   1377000*105**(0.5D0)*Q(31) + 793800*105**(0.5D0)*Q(55)))/
     >   CONST1)*VOL(I2,J2,J,I,K2) 
      Q2(30,65) = (-(7**(0.5D0)*(1062075000*3**(0.5D0)*Q(22) - 
     >   337365000*7**(0.5D0)*Q(24) + 433755000*3**(0.5D0)*Q(32) - 
     >   585225000*7**(0.5D0)*Q(30) + 250047000*3**(0.5D0)*Q(56) - 
     >   337365000*7**(0.5D0)*Q(54) + 433755000*3**(0.5D0)*Q(62) - 
     >   137781000*7**(0.5D0)*Q(64)))/CONST5)*VOL(I2,J2,J,I,K2) 
      Q2(31,65) = ((105**(0.5D0)*(12495000*Q(21) + 5103000*Q(61) - 
     >   7497000*5**(0.5D0)*Q(23) - 2295000*21**(0.5D0)*Q(29) - 
     >   3061800*5**(0.5D0)*Q(63) - 1323000*21**(0.5D0)*Q(53) + 
     >   1377000*105**(0.5D0)*Q(31) + 793800*105**(0.5D0)*Q(55)))/
     >   CONST6)*VOL(I2,J2,J,I,K2) 
      Q2(32,65) = ((3**(0.5D0)*(612255000*3**(0.5D0)*Q(22) - 
     >   337365000*7**(0.5D0)*Q(24) + 433755000*3**(0.5D0)*Q(32) - 
     >   337365000*7**(0.5D0)*Q(30) + 250047000*3**(0.5D0)*Q(56) - 
     >   194481000*7**(0.5D0)*Q(54) + 250047000*3**(0.5D0)*Q(62) - 
     >   137781000*7**(0.5D0)*Q(64)))/CONST7)*VOL(I2,J2,J,I,K2) 
      Q2(33,65) = (-(3*5**(0.5D0)*(25000*Q(01) + 1800*Q(11) + 
     >   9000*Q(35) + 9000*Q(41) - 3000*5**(0.5D0)*Q(03) - 
     >   3000*5**(0.5D0)*Q(09) - 15000*5**(0.5D0)*Q(33) - 
     >   1080*5**(0.5D0)*Q(43)))/20480000)*VOL(I2,J2,J,I,K2) 
      Q2(34,65) = (-(15**(0.5D0)*(425000*3**(0.5D0)*Q(02) - 
     >   135000*7**(0.5D0)*Q(04) - 51000*15**(0.5D0)*Q(10) + 
     >   153000*3**(0.5D0)*Q(42) + 16200*35**(0.5D0)*Q(12) - 
     >   255000*15**(0.5D0)*Q(34) - 48600*7**(0.5D0)*Q(44) + 
     >   81000*35**(0.5D0)*Q(36)))/307200000)*VOL(I2,J2,J,I,K2) 
      Q2(35,65) = ((9*(5000*Q(01) + 1800*Q(11) + 9000*Q(35) + 
     >   1800*Q(41) - 3000*5**(0.5D0)*Q(03) - 600*5**(0.5D0)*
     >   Q(09) - 3000*5**(0.5D0)*Q(33) - 1080*5**(0.5D0)*
     >   Q(43)))/20480000)*VOL(I2,J2,J,I,K2) 
      Q2(36,65) = ((9*35**(0.5D0)*(245000*3**(0.5D0)*Q(02) - 
     >   135000*7**(0.5D0)*Q(04) - 29400*15**(0.5D0)*Q(10) + 
     >   88200*3**(0.5D0)*Q(42) + 16200*35**(0.5D0)*Q(12) - 
     >   147000*15**(0.5D0)*Q(34) - 48600*7**(0.5D0)*Q(44) + 
     >   81000*35**(0.5D0)*Q(36)))/CONST2)*VOL(I2,J2,J,I,K2) 
      Q2(37,65) = (-(15**(0.5D0)*(425000*3**(0.5D0)*Q(05) - 
     >   135000*7**(0.5D0)*Q(13) - 51000*15**(0.5D0)*Q(07) + 
     >   153000*3**(0.5D0)*Q(39) + 16200*35**(0.5D0)*Q(15) - 
     >   255000*15**(0.5D0)*Q(37) - 48600*7**(0.5D0)*Q(47) + 
     >   81000*35**(0.5D0)*Q(45)))/307200000)*VOL(I2,J2,J,I,K2) 
      Q2(38,65) = (-(5**(0.5D0)*(21675000*Q(06) + 5103000*Q(16) - 
     >   2295000*21**(0.5D0)*Q(08) - 2295000*21**(0.5D0)*Q(14) - 
     >   13005000*5**(0.5D0)*Q(38) - 3061800*5**(0.5D0)*Q(48) + 
     >   1377000*105**(0.5D0)*Q(40) + 1377000*105**(0.5D0)*
     >   Q(46)))/CONST4)*VOL(I2,J2,J,I,K2) 
      Q2(39,65) = ((3**(0.5D0)*(85000*3**(0.5D0)*Q(05) - 
     >   27000*7**(0.5D0)*Q(13) - 51000*15**(0.5D0)*Q(07) + 
     >   153000*3**(0.5D0)*Q(39) + 16200*35**(0.5D0)*Q(15) - 
     >   51000*15**(0.5D0)*Q(37) - 48600*7**(0.5D0)*Q(47) + 
     >   16200*35**(0.5D0)*Q(45)))/102400000)*VOL(I2,J2,J,I,K2) 
      Q2(40,65) = ((105**(0.5D0)*(12495000*Q(06) + 5103000*Q(16) - 
     >   2295000*21**(0.5D0)*Q(08) - 1323000*21**(0.5D0)*Q(14) - 
     >   7497000*5**(0.5D0)*Q(38) - 3061800*5**(0.5D0)*Q(48) + 
     >   1377000*105**(0.5D0)*Q(40) + 793800*105**(0.5D0)*
     >   Q(46)))/CONST6)*VOL(I2,J2,J,I,K2) 
      Q2(41,65) = ((9*(5000*Q(01) + 1800*Q(11) + 1800*Q(35) + 
     >   9000*Q(41) - 600*5**(0.5D0)*Q(03) - 3000*5**(0.5D0)*Q(09) - 
     >   3000*5**(0.5D0)*Q(33) - 1080*5**(0.5D0)*Q(43)))/20480000)*
     >   VOL(I2,J2,J,I,K2) 
      Q2(42,65) = ((3**(0.5D0)*(85000*3**(0.5D0)*Q(02) - 
     >   27000*7**(0.5D0)*Q(04) - 51000*15**(0.5D0)*Q(10) + 
     >   153000*3**(0.5D0)*Q(42) + 16200*35**(0.5D0)*Q(12) - 
     >   51000*15**(0.5D0)*Q(34) - 48600*7**(0.5D0)*Q(44) + 
     >   16200*35**(0.5D0)*Q(36)))/102400000)*VOL(I2,J2,J,I,K2) 
      Q2(43,65) = (-(27*5**(0.5D0)*(1000*Q(01) + 1800*Q(11) + 
     >   1800*Q(35) + 1800*Q(41) - 600*5**(0.5D0)*Q(03) - 
     >   600*5**(0.5D0)*Q(09) - 600*5**(0.5D0)*Q(33) - 
     >   1080*5**(0.5D0)*Q(43)))/102400000)*VOL(I2,J2,J,I,K2) 
      Q2(44,65) = (-(27*7**(0.5D0)*(49000*3**(0.5D0)*Q(02) - 
     >   27000*7**(0.5D0)*Q(04) - 29400*15**(0.5D0)*Q(10) + 
     >   88200*3**(0.5D0)*Q(42) + 16200*35**(0.5D0)*Q(12) - 
     >   29400*15**(0.5D0)*Q(34) - 48600*7**(0.5D0)*Q(44) + 
     >   16200*35**(0.5D0)*Q(36)))/CONST2)*VOL(I2,J2,J,I,K2) 
      Q2(45,65) = ((9*35**(0.5D0)*(245000*3**(0.5D0)*Q(05) - 
     >   135000*7**(0.5D0)*Q(13) - 29400*15**(0.5D0)*Q(07) + 
     >   88200*3**(0.5D0)*Q(39) + 16200*35**(0.5D0)*Q(15) - 
     >   147000*15**(0.5D0)*Q(37) - 48600*7**(0.5D0)*Q(47) + 
     >   81000*35**(0.5D0)*Q(45)))/CONST2)*VOL(I2,J2,J,I,K2) 
      Q2(46,65) = ((105**(0.5D0)*(12495000*Q(06) + 5103000*Q(16) - 
     >   1323000*21**(0.5D0)*Q(08) - 2295000*21**(0.5D0)*Q(14) - 
     >   7497000*5**(0.5D0)*Q(38) - 3061800*5**(0.5D0)*Q(48) + 
     >   793800*105**(0.5D0)*Q(40) + 1377000*105**(0.5D0)*
     >   Q(46)))/CONST6)*VOL(I2,J2,J,I,K2) 
      Q2(47,65) = (-(27*7**(0.5D0)*(49000*3**(0.5D0)*Q(05) - 
     >   27000*7**(0.5D0)*Q(13) - 29400*15**(0.5D0)*Q(07) + 
     >   88200*3**(0.5D0)*Q(39) + 16200*35**(0.5D0)*Q(15) - 
     >   29400*15**(0.5D0)*Q(37) - 48600*7**(0.5D0)*Q(47) + 
     >   16200*35**(0.5D0)*Q(45)))/CONST2)*VOL(I2,J2,J,I,K2) 
      Q2(48,65) = (-(27*5**(0.5D0)*(7203000*Q(06) + 
     >   5103000*Q(16) - 1323000*21**(0.5D0)*Q(08) - 
     >   1323000*21**(0.5D0)*Q(14) - 4321800*5**(0.5D0)*Q(38) - 
     >   3061800*5**(0.5D0)*Q(48) + 793800*105**(0.5D0)*Q(40) + 
     >   793800*105**(0.5D0)*Q(46)))/CONST7)*VOL(I2,J2,J,I,K2) 
      Q2(49,65) = (-(3*7**(0.5D0)*(1225000*3**(0.5D0)*Q(17) + 
     >   88200*3**(0.5D0)*Q(27) - 147000*15**(0.5D0)*Q(19) - 
     >   147000*15**(0.5D0)*Q(25) - 675000*7**(0.5D0)*Q(49) - 
     >   48600*7**(0.5D0)*Q(59) + 81000*35**(0.5D0)*Q(51) + 
     >   81000*35**(0.5D0)*Q(57)))/1003520000)*VOL(I2,J2,J,I,K2) 
      Q2(50,65) = (-(21**(0.5D0)*(62475000*Q(18) + 25515000*Q(52) - 
     >   7497000*5**(0.5D0)*Q(26) - 6615000*21**(0.5D0)*Q(20) - 
     >   3061800*5**(0.5D0)*Q(60) - 11475000*21**(0.5D0)*Q(50) + 
     >   793800*105**(0.5D0)*Q(28) + 1377000*105**(0.5D0)*
     >   Q(58)))/CONST1)*VOL(I2,J2,J,I,K2) 
      Q2(51,65) = ((9*35**(0.5D0)*(245000*3**(0.5D0)*Q(17) + 
     >   88200*3**(0.5D0)*Q(27) - 147000*15**(0.5D0)*Q(19) - 
     >   29400*15**(0.5D0)*Q(25) - 135000*7**(0.5D0)*Q(49) - 
     >   48600*7**(0.5D0)*Q(59) + 81000*35**(0.5D0)*Q(51) + 
     >   16200*35**(0.5D0)*Q(57)))/CONST2)*VOL(I2,J2,J,I,K2) 
      Q2(52,65) = ((9*(36015000*Q(18) + 25515000*Q(52) - 
     >   4321800*5**(0.5D0)*Q(26) - 6615000*21**(0.5D0)*Q(20) - 
     >   3061800*5**(0.5D0)*Q(60) - 6615000*21**(0.5D0)*Q(50) + 
     >   793800*105**(0.5D0)*Q(28) + 793800*105**(0.5D0)*
     >   Q(58)))/CONST3)*VOL(I2,J2,J,I,K2) 
      Q2(53,65) = (-(21**(0.5D0)*(62475000*Q(21) + 25515000*Q(61) - 
     >   7497000*5**(0.5D0)*Q(23) - 6615000*21**(0.5D0)*Q(29) - 
     >   3061800*5**(0.5D0)*Q(63) - 11475000*21**(0.5D0)*Q(53) + 
     >   793800*105**(0.5D0)*Q(31) + 1377000*105**(0.5D0)*
     >   Q(55)))/CONST1)*VOL(I2,J2,J,I,K2) 
      Q2(54,65) = (-(7**(0.5D0)*(1062075000*3**(0.5D0)*Q(22) - 
     >   337365000*7**(0.5D0)*Q(24) + 250047000*3**(0.5D0)*Q(32) - 
     >   337365000*7**(0.5D0)*Q(30) + 433755000*3**(0.5D0)*Q(56) - 
     >   585225000*7**(0.5D0)*Q(54) + 433755000*3**(0.5D0)*Q(62) - 
     >   137781000*7**(0.5D0)*Q(64)))/CONST5)*VOL(I2,J2,J,I,K2) 
      Q2(55,65) = ((105**(0.5D0)*(12495000*Q(21) + 5103000*Q(61) - 
     >   7497000*5**(0.5D0)*Q(23) - 1323000*21**(0.5D0)*Q(29) - 
     >   3061800*5**(0.5D0)*Q(63) - 2295000*21**(0.5D0)*Q(53) + 
     >   793800*105**(0.5D0)*Q(31) + 1377000*105**(0.5D0)*
     >   Q(55)))/CONST6)*VOL(I2,J2,J,I,K2) 
      Q2(56,65) = ((3**(0.5D0)*(612255000*3**(0.5D0)*Q(22) - 
     >   337365000*7**(0.5D0)*Q(24) + 250047000*3**(0.5D0)*Q(32) - 
     >   194481000*7**(0.5D0)*Q(30) + 433755000*3**(0.5D0)*Q(56) - 
     >   337365000*7**(0.5D0)*Q(54) + 250047000*3**(0.5D0)*Q(62) - 
     >   137781000*7**(0.5D0)*Q(64)))/CONST7)*VOL(I2,J2,J,I,K2) 
      Q2(57,65) = ((9*35**(0.5D0)*(245000*3**(0.5D0)*Q(17) + 
     >   88200*3**(0.5D0)*Q(27) - 29400*15**(0.5D0)*Q(19) - 
     >   147000*15**(0.5D0)*Q(25) - 135000*7**(0.5D0)*Q(49) - 
     >   48600*7**(0.5D0)*Q(59) + 16200*35**(0.5D0)*Q(51) + 
     >   81000*35**(0.5D0)*Q(57)))/CONST2)*VOL(I2,J2,J,I,K2) 
      Q2(58,65) = ((105**(0.5D0)*(12495000*Q(18) + 5103000*Q(52) - 
     >   7497000*5**(0.5D0)*Q(26) - 1323000*21**(0.5D0)*Q(20) - 
     >   3061800*5**(0.5D0)*Q(60) - 2295000*21**(0.5D0)*Q(50) + 
     >   793800*105**(0.5D0)*Q(28) + 1377000*105**(0.5D0)*
     >   Q(58)))/CONST6)*VOL(I2,J2,J,I,K2) 
      Q2(59,65) = (-(27*7**(0.5D0)*(49000*3**(0.5D0)*Q(17) + 
     >   88200*3**(0.5D0)*Q(27) - 29400*15**(0.5D0)*Q(19) - 
     >   29400*15**(0.5D0)*Q(25) - 27000*7**(0.5D0)*Q(49) - 
     >   48600*7**(0.5D0)*Q(59) + 16200*35**(0.5D0)*Q(51) + 
     >   16200*35**(0.5D0)*Q(57)))/CONST2)*VOL(I2,J2,J,I,K2) 
      Q2(60,65) = (-(27*5**(0.5D0)*(7203000*Q(18) + 5103000*Q(52) - 
     >   4321800*5**(0.5D0)*Q(26) - 1323000*21**(0.5D0)*Q(20) - 
     >   3061800*5**(0.5D0)*Q(60) - 1323000*21**(0.5D0)*Q(50) + 
     >   793800*105**(0.5D0)*Q(28) + 793800*105**(0.5D0)*
     >   Q(58)))/CONST7)*VOL(I2,J2,J,I,K2) 
      Q2(61,65) = ((9*(36015000*Q(21) + 25515000*Q(61) - 
     >   4321800*5**(0.5D0)*Q(23) - 6615000*21**(0.5D0)*Q(29) - 
     >   3061800*5**(0.5D0)*Q(63) - 6615000*21**(0.5D0)*Q(53) + 
     >   793800*105**(0.5D0)*Q(31) + 793800*105**(0.5D0)*
     >   Q(55)))/CONST3)*VOL(I2,J2,J,I,K2) 
      Q2(62,65) = ((3**(0.5D0)*(612255000*3**(0.5D0)*Q(22) - 
     >   194481000*7**(0.5D0)*Q(24) + 250047000*3**(0.5D0)*Q(32) - 
     >   337365000*7**(0.5D0)*Q(30) + 250047000*3**(0.5D0)*Q(56) - 
     >   337365000*7**(0.5D0)*Q(54) + 433755000*3**(0.5D0)*Q(62) - 
     >   137781000*7**(0.5D0)*Q(64)))/CONST7)*VOL(I2,J2,J,I,K2) 
      Q2(63,65) = (-(27*5**(0.5D0)*(7203000*Q(21) + 5103000*Q(61) - 
     >   4321800*5**(0.5D0)*Q(23) - 1323000*21**(0.5D0)*Q(29) - 
     >   3061800*5**(0.5D0)*Q(63) - 1323000*21**(0.5D0)*Q(53) + 
     >   793800*105**(0.5D0)*Q(31) + 793800*105**(0.5D0)*
     >   Q(55)))/CONST7)*VOL(I2,J2,J,I,K2) 
      Q2(64,65) = (-(27*7**(0.5D0)*(352947000*3**(0.5D0)*Q(22) - 
     >   194481000*7**(0.5D0)*Q(24) + 250047000*3**(0.5D0)*Q(32) - 
     >   194481000*7**(0.5D0)*Q(30) + 250047000*3**(0.5D0)*Q(56) - 
     >   194481000*7**(0.5D0)*Q(54) + 250047000*3**(0.5D0)*Q(62) - 
     >   137781000*7**(0.5D0)*Q(64)))/CONST8)*VOL(I2,J2,J,I,K2) 

      Q2(1,65) = Q2(1,65) + 
     >   ((21*XNI(4,3,J2,K2))/51200 + (441*XNI(3,3,J2,K2))/204800 + 
     >   (441*XNI(2,3,J2,K2))/204800 + (21*XNI(1,3,J2,K2))/51200 + 
     >   XNI(4,4,J2,K2)/12800 + (21*XNI(3,4,J2,K2))/51200 + 
     >   (21*XNI(2,4,J2,K2))/51200 + XNI(1,4,J2,K2)/12800 + 
     >   XNI(4,1,J2,K2)/12800 + (21*XNI(3,1,J2,K2))/51200 + 
     >   (21*XNI(2,1,J2,K2))/51200 + XNI(1,1,J2,K2)/12800 + 
     >   (21*XNI(4,2,J2,K2))/51200 + (441*XNI(3,2,J2,K2))/204800 + 
     >   (441*XNI(2,2,J2,K2))/204800 + (21*XNI(1,2,J2,K2))/51200)*
     >   MUH*SIGN(1.0,MUH) + (XNJ(1,4,K2)/12800 + 
     >   (21*XNJ(2,4,K2))/51200 + (21*XNJ(3,4,K2))/51200 + 
     >   XNJ(4,4,K2)/12800 + (21*XNJ(1,3,K2))/51200 + 
     >   (441*XNJ(2,3,K2))/204800 + (441*XNJ(3,3,K2))/204800 + 
     >   (21*XNJ(4,3,K2))/51200 + XNJ(1,1,K2)/12800 + 
     >   (21*XNJ(2,1,K2))/51200 + (21*XNJ(3,1,K2))/51200 + 
     >   XNJ(4,1,K2)/12800 + (21*XNJ(1,2,K2))/51200 + 
     >   (441*XNJ(2,2,K2))/204800 + (441*XNJ(3,2,K2))/204800 + 
     >   (21*XNJ(4,2,K2))/51200)*ETAH*SIGN(1.0,ETAH) + 
     >   ((441*XNK(3,2))/204800 + (441*XNK(3,3))/204800 + 
     >   XNK(1,4)/12800 + (21*XNK(1,3))/51200 + (21*XNK(1,2))/51200 + 
     >   XNK(1,1)/12800 + (21*XNK(2,1))/51200 + (21*XNK(3,1))/51200 + 
     >   XNK(4,1)/12800 + (21*XNK(4,2))/51200 + (21*XNK(4,3))/51200 + 
     >   XNK(4,4)/12800 + (21*XNK(3,4))/51200 + (21*XNK(2,4))/51200 + 
     >   (441*XNK(2,3))/204800 + (441*XNK(2,2))/204800)*XIH*
     >   SIGN(1.0,XIH) 
      Q2(2,65) = Q2(2,65) + 
     >   (3**(0.5D0)*XIH*SIGN(1.0,XIH)*(15876*XNK(3,2) + 
     >   15876*XNK(3,3) - 692*XNK(1,4) - 3633*XNK(1,3) - 
     >   3633*XNK(1,2) - 692*XNK(1,1) - 3024*XNK(2,1) + 
     >   3024*XNK(3,1) + 692*XNK(4,1) + 3633*XNK(4,2) + 
     >   3633*XNK(4,3) + 692*XNK(4,4) + 3024*XNK(3,4) - 
     >   3024*XNK(2,4) - 15876*XNK(2,3) - 15876*XNK(2,2)))/9216000 - 
     >   (3**(0.5D0)*MUH*(2772*XNI(4,3,J2,K2) + 
     >   14553*XNI(3,3,J2,K2) + 14553*XNI(2,3,J2,K2) + 2772*
     >   XNI(1,3,J2,K2) + 
     >   528*XNI(4,4,J2,K2) + 2772*XNI(3,4,J2,K2) + 2772*XNI(2,4,J2,K2)+
     >   528*XNI(1,4,J2,K2) + 528*XNI(4,1,J2,K2) + 2772*XNI(3,1,J2,K2) +
     >   2772*XNI(2,1,J2,K2) + 528*XNI(1,1,J2,K2) + 2772*XNI(4,2,J2,K2)+
     >   14553*XNI(3,2,J2,K2) + 14553*XNI(2,2,J2,K2) + 2772*
     >   XNI(1,2,J2,K2)))/9216000 - (3**(0.5D0)*ETAH*
     >   SIGN(1.0,ETAH)*(692*XNJ(1,4,K2) + 3024*XNJ(2,4,K2) - 
     >   3024*XNJ(3,4,K2) - 692*XNJ(4,4,K2) + 3633*XNJ(1,3,K2) + 
     >   15876*XNJ(2,3,K2) - 15876*XNJ(3,3,K2) - 3633*XNJ(4,3,K2) + 
     >   692*XNJ(1,1,K2) + 3024*XNJ(2,1,K2) - 3024*XNJ(3,1,K2) - 
     >   692*XNJ(4,1,K2) + 3633*XNJ(1,2,K2) + 15876*XNJ(2,2,K2) - 
     >   15876*XNJ(3,2,K2) - 3633*XNJ(4,2,K2)))/9216000 
      Q2(3,65) = Q2(3,65) + 
     >   (3*5**(0.5D0)*ETAH*SIGN(1.0,ETAH)*(4*XNJ(1,4,K2) - 
     >   24*XNJ(2,4,K2) - 24*XNJ(3,4,K2) + 4*XNJ(4,4,K2) + 
     >   21*XNJ(1,3,K2) - 126*XNJ(2,3,K2) - 126*XNJ(3,3,K2) + 
     >   21*XNJ(4,3,K2) + 4*XNJ(1,1,K2) - 24*XNJ(2,1,K2) - 
     >   24*XNJ(3,1,K2) + 4*XNJ(4,1,K2) + 21*XNJ(1,2,K2) - 
     >   126*XNJ(2,2,K2) - 126*XNJ(3,2,K2) + 21*XNJ(4,2,K2)))/1024000 +
     >   (3*5**(0.5D0)*MUH*SIGN(1.0,MUH)*(84*XNI(4,3,J2,K2) + 
     >   441*XNI(3,3,J2,K2) + 441*XNI(2,3,J2,K2) + 84*XNI(1,3,J2,K2) + 
     >   16*XNI(4,4,J2,K2) + 84*XNI(3,4,J2,K2) + 84*XNI(2,4,J2,K2) + 
     >   16*XNI(1,4,J2,K2) + 16*XNI(4,1,J2,K2) + 84*XNI(3,1,J2,K2) + 
     >   84*XNI(2,1,J2,K2) + 16*XNI(1,1,J2,K2) + 84*XNI(4,2,J2,K2) + 
     >   441*XNI(3,2,J2,K2) + 441*XNI(2,2,J2,K2) + 84*XNI(1,2,J2,K2)))/
     >   1024000 - (3*5**(0.5D0)*XIH*SIGN(1.0,XIH)*(126*
     >   XNK(3,2) + 126*XNK(3,3) - 4*XNK(1,4) - 21*XNK(1,3) - 
     >   21*XNK(1,2) - 4*XNK(1,1) + 24*XNK(2,1) + 24*XNK(3,1) - 
     >   4*XNK(4,1) - 21*XNK(4,2) - 21*XNK(4,3) - 4*XNK(4,4) + 
     >   24*XNK(3,4) + 24*XNK(2,4) + 126*XNK(2,3) + 126*XNK(2,2)))/
     >   1024000 
      Q2(4,65) = Q2(4,65) + 
     >   (3*7**(0.5D0)*ETAH*SIGN(1.0,ETAH)*(296*XNJ(1,4,K2) + 
     >   2052*XNJ(2,4,K2) - 2052*XNJ(3,4,K2) - 296*XNJ(4,4,K2) + 
     >   1554*XNJ(1,3,K2) + 10773*XNJ(2,3,K2) - 10773*XNJ(3,3,K2) - 
     >   1554*XNJ(4,3,K2) + 296*XNJ(1,1,K2) + 2052*XNJ(2,1,K2) - 
     >   2052*XNJ(3,1,K2) - 296*XNJ(4,1,K2) + 1554*XNJ(1,2,K2) + 
     >   10773*XNJ(2,2,K2) - 10773*XNJ(3,2,K2) - 1554*XNJ(4,2,K2)))/
     >   50176000 - (3*7**(0.5D0)*MUH*(1764*XNI(4,3,J2,K2) + 
     >   9261*XNI(3,3,J2,K2) + 9261*XNI(2,3,J2,K2) + 1764*
     >   XNI(1,3,J2,K2) + 
     >   336*XNI(4,4,J2,K2) + 1764*XNI(3,4,J2,K2) + 1764*XNI(2,4,J2,K2)+
     >   336*XNI(1,4,J2,K2) + 336*XNI(4,1,J2,K2) + 1764*XNI(3,1,J2,K2) +
     >   1764*XNI(2,1,J2,K2) + 336*XNI(1,1,J2,K2) + 1764*XNI(4,2,J2,K2)+
     >   9261*XNI(3,2,J2,K2) + 9261*XNI(2,2,J2,K2) + 1764*
     >   XNI(1,2,J2,K2)))/
     >   50176000 - (3*7**(0.5D0)*XIH*SIGN(1.0,XIH)*(10773*
     >   XNK(3,2) + 10773*XNK(3,3) - 296*XNK(1,4) - 1554*XNK(1,3) - 
     >   1554*XNK(1,2) - 296*XNK(1,1) - 2052*XNK(2,1) + 2052*
     >   XNK(3,1) + 296*XNK(4,1) + 1554*XNK(4,2) + 1554*XNK(4,3) + 
     >   296*XNK(4,4) + 2052*XNK(3,4) - 2052*XNK(2,4) - 10773*
     >   XNK(2,3) - 10773*XNK(2,2)))/50176000 
      Q2(5,65) = Q2(5,65) + 
     >   (3**(0.5D0)*MUH*SIGN(1.0,MUH)*(3633*XNI(4,3,J2,K2) + 
     >   15876*XNI(3,3,J2,K2) - 15876*XNI(2,3,J2,K2) - 3633*
     >   XNI(1,3,J2,K2) + 
     >   692*XNI(4,4,J2,K2) + 3024*XNI(3,4,J2,K2) - 3024*XNI(2,4,J2,K2)-
     >   692*XNI(1,4,J2,K2) + 692*XNI(4,1,J2,K2) + 3024*XNI(3,1,J2,K2)-
     >   3024*XNI(2,1,J2,K2) - 692*XNI(1,1,J2,K2) + 3633*XNI(4,2,J2,K2)+
     >   15876*XNI(3,2,J2,K2) - 15876*XNI(2,2,J2,K2) - 3633*
     >   XNI(1,2,J2,K2)))/9216000 - (3**(0.5D0)*ETAH*(528*
     >   XNJ(1,4,K2) + 2772*XNJ(2,4,K2) + 2772*XNJ(3,4,K2) + 528*
     >   XNJ(4,4,K2) + 2772*XNJ(1,3,K2) + 14553*XNJ(2,3,K2) + 14553*
     >   XNJ(3,3,K2) + 2772*XNJ(4,3,K2) + 528*XNJ(1,1,K2) + 2772*
     >   XNJ(2,1,K2) + 2772*XNJ(3,1,K2) + 528*XNJ(4,1,K2) + 2772*
     >   XNJ(1,2,K2) + 14553*XNJ(2,2,K2) + 14553*XNJ(3,2,K2) + 2772*
     >   XNJ(4,2,K2)))/9216000 - (3**(0.5D0)*XIH*SIGN(1.0,XIH)*
     >   (15876*XNK(3,2) - 15876*XNK(3,3) - 692*XNK(1,4) - 3024*
     >   XNK(1,3) + 3024*XNK(1,2) + 692*XNK(1,1) + 3633*XNK(2,1) + 
     >   3633*XNK(3,1) + 692*XNK(4,1) + 3024*XNK(4,2) - 3024*
     >   XNK(4,3) - 692*XNK(4,4) - 3633*XNK(3,4) - 3633*XNK(2,4) - 
     >   15876*XNK(2,3) + 15876*XNK(2,2)))/9216000 
      Q2(6,65) = Q2(6,65) + 
     >   ((4851*XNI(2,3,J2,K2))/1280000 - (4851*XNI(3,3,J2,K2))/1280000-
     >   (13321*XNI(4,3,J2,K2))/15360000 + (13321*XNI(1,3,J2,K2))/
     >   15360000 - (1903*XNI(4,4,J2,K2))/11520000 - (231*
     >   XNI(3,4,J2,K2))/320000 + (231*XNI(2,4,J2,K2))/320000 + (1903*
     >   XNI(1,4,J2,K2))/11520000 - (1903*XNI(4,1,J2,K2))/11520000 - 
     >   (231*XNI(3,1,J2,K2))/320000 + (231*XNI(2,1,J2,K2))/320000 + 
     >   (1903*XNI(1,1,J2,K2))/11520000 - (13321*XNI(4,2,J2,K2))/
     >   15360000 - (4851*XNI(3,2,J2,K2))/1280000 + (4851*
     >   XNI(2,2,J2,K2))/1280000 + (13321*XNI(1,2,J2,K2))/15360000)*
     >   MUH + ((1903*XNJ(1,4,K2))/11520000 + (231*XNJ(2,4,K2))/
     >   320000 - (231*XNJ(3,4,K2))/320000 - (1903*XNJ(4,4,K2))/
     >   11520000 + (13321*XNJ(1,3,K2))/15360000 + (4851*XNJ(2,3,K2))/
     >   1280000 - (4851*XNJ(3,3,K2))/1280000 - (13321*XNJ(4,3,K2))/
     >   15360000 + (1903*XNJ(1,1,K2))/11520000 + (231*XNJ(2,1,K2))/
     >   320000 - (231*XNJ(3,1,K2))/320000 - (1903*XNJ(4,1,K2))/
     >   11520000 + (13321*XNJ(1,2,K2))/15360000 + (4851*XNJ(2,2,K2))/
     >   1280000 - (4851*XNJ(3,2,K2))/1280000 - (13321*XNJ(4,2,K2))/
     >   15360000)*ETAH + ((1323*XNK(3,3))/320000 - (1323*
     >   XNK(3,2))/320000 - (29929*XNK(1,4))/138240000 - (1211*
     >   XNK(1,3))/1280000 + (1211*XNK(1,2))/1280000 + (29929*
     >   XNK(1,1))/138240000 + (1211*XNK(2,1))/1280000 - (1211*
     >   XNK(3,1))/1280000 - (29929*XNK(4,1))/138240000 - (1211*
     >   XNK(4,2))/1280000 + (1211*XNK(4,3))/1280000 + (29929*
     >   XNK(4,4))/138240000 + (1211*XNK(3,4))/1280000 - (1211*
     >   XNK(2,4))/1280000 - (1323*XNK(2,3))/320000 + (1323*
     >   XNK(2,2))/320000)*XIH*SIGN(1.0,XIH) 
      Q2(7,65) = Q2(7,65) + 
     >   (15**(0.5D0)*MUH*SIGN(1.0,MUH)*(3633*XNI(4,3,J2,K2) + 
     >   15876*XNI(3,3,J2,K2) - 15876*XNI(2,3,J2,K2) - 3633*
     >   XNI(1,3,J2,K2) + 
     >   692*XNI(4,4,J2,K2) + 3024*XNI(3,4,J2,K2) - 3024*XNI(2,4,J2,K2)-
     >   692*XNI(1,4,J2,K2) + 692*XNI(4,1,J2,K2) + 3024*XNI(3,1,J2,K2) -
     >   3024*XNI(2,1,J2,K2) - 692*XNI(1,1,J2,K2) + 3633*XNI(4,2,J2,K2)+
     >   15876*XNI(3,2,J2,K2) - 15876*XNI(2,2,J2,K2) - 3633*
     >   XNI(1,2,J2,K2)))/15360000 - (15**(0.5D0)*ETAH*(132*
     >   XNJ(1,4,K2) - 792*XNJ(2,4,K2) - 792*XNJ(3,4,K2) + 132*
     >   XNJ(4,4,K2) + 693*XNJ(1,3,K2) - 4158*XNJ(2,3,K2) - 4158*
     >   XNJ(3,3,K2) + 693*XNJ(4,3,K2) + 132*XNJ(1,1,K2) - 792*
     >   XNJ(2,1,K2) - 792*XNJ(3,1,K2) + 132*XNJ(4,1,K2) + 693*
     >   XNJ(1,2,K2) - 4158*XNJ(2,2,K2) - 4158*XNJ(3,2,K2) + 693*
     >   XNJ(4,2,K2)))/15360000 + (15**(0.5D0)*XIH*
     >   SIGN(1.0,XIH)*(4536*XNK(3,2) - 4536*XNK(3,3) + 
     >  173*XNK(1,4) + 756*XNK(1,3) - 756*XNK(1,2) - 173*XNK(1,1) + 
     >  1038*XNK(2,1) + 1038*XNK(3,1) - 173*XNK(4,1) - 756*XNK(4,2) + 
     >   756*XNK(4,3) + 173*XNK(4,4) - 1038*XNK(3,4) - 1038*XNK(2,4) - 
     >   4536*XNK(2,3) + 4536*XNK(2,2)))/15360000 
      Q2(8,65) = Q2(8,65) + 
     >   (21**(0.5D0)*XIH*SIGN(1.0,XIH)*(387828*XNK(3,2) - 
     >   387828*XNK(3,3) + 12802*XNK(1,4) + 55944*XNK(1,3) - 
     >   55944*XNK(1,2) - 12802*XNK(1,1) - 88749*XNK(2,1) + 
     >   88749*XNK(3,1) + 12802*XNK(4,1) + 55944*XNK(4,2) - 
     >   55944*XNK(4,3) - 12802*XNK(4,4) - 88749*XNK(3,4) + 
     >   88749*XNK(2,4) + 387828*XNK(2,3) - 387828*XNK(2,2)))/
     >   752640000 - (21**(0.5D0)*MUH*(76293*XNI(4,3,J2,K2) + 
     >   333396*XNI(3,3,J2,K2) - 333396*XNI(2,3,J2,K2) - 76293*
     >   XNI(1,3,J2,K2) + 14532*XNI(4,4,J2,K2) + 63504*XNI(3,4,J2,K2) -
     >   63504*XNI(2,4,J2,K2) - 14532*XNI(1,4,J2,K2) + 14532*
     >   XNI(4,1,J2,K2) + 63504*XNI(3,1,J2,K2) - 63504*XNI(2,1,J2,K2) -
     >   14532*XNI(1,1,J2,K2) + 76293*XNI(4,2,J2,K2) + 333396*
     >   XNI(3,2,J2,K2) - 333396*XNI(2,2,J2,K2) - 76293*
     >   XNI(1,2,J2,K2)))/
     >   752640000 - (21**(0.5D0)*ETAH*(9768*XNJ(1,4,K2) + 
     >   67716*XNJ(2,4,K2) - 67716*XNJ(3,4,K2) - 9768*XNJ(4,4,K2) + 
     >   51282*XNJ(1,3,K2) + 355509*XNJ(2,3,K2) - 355509*XNJ(3,3,K2) - 
     >   51282*XNJ(4,3,K2) + 9768*XNJ(1,1,K2) + 67716*XNJ(2,1,K2) - 
     >   67716*XNJ(3,1,K2) - 9768*XNJ(4,1,K2) + 51282*XNJ(1,2,K2) + 
     >   355509*XNJ(2,2,K2) - 355509*XNJ(3,2,K2) - 51282*
     >   XNJ(4,2,K2)))/
     >   752640000 
      Q2(9,65) = Q2(9,65) + 
     >   (3*5**(0.5D0)*ETAH*SIGN(1.0,ETAH)*(16*XNJ(1,4,K2) + 
     >   84*XNJ(2,4,K2) + 84*XNJ(3,4,K2) + 16*XNJ(4,4,K2) + 
     >   84*XNJ(1,3,K2) + 441*XNJ(2,3,K2) + 441*XNJ(3,3,K2) + 
     >   84*XNJ(4,3,K2) + 16*XNJ(1,1,K2) + 84*XNJ(2,1,K2) + 
     >   84*XNJ(3,1,K2) + 16*XNJ(4,1,K2) + 84*XNJ(1,2,K2) + 
     >   441*XNJ(2,2,K2) + 441*XNJ(3,2,K2) + 84*XNJ(4,2,K2)))/1024000 +
     >   (3*5**(0.5D0)*MUH*SIGN(1.0,MUH)*(21*XNI(4,3,J2,K2) - 
     >   126*XNI(3,3,J2,K2) - 126*XNI(2,3,J2,K2) + 21*XNI(1,3,J2,K2) + 
     >   4*XNI(4,4,J2,K2) - 24*XNI(3,4,J2,K2) - 24*XNI(2,4,J2,K2) + 
     >   4*XNI(1,4,J2,K2) + 4*XNI(4,1,J2,K2) - 24*XNI(3,1,J2,K2) - 
     >   24*XNI(2,1,J2,K2) + 4*XNI(1,1,J2,K2) + 21*XNI(4,2,J2,K2) - 
     >   126*XNI(3,2,J2,K2) - 126*XNI(2,2,J2,K2) + 21*XNI(1,2,J2,K2)))/
     >   1024000 - (3*5**(0.5D0)*XIH*SIGN(1.0,XIH)*(126*
     >   XNK(3,2) + 126*XNK(3,3) - 4*XNK(1,4) + 24*XNK(1,3) + 
     >   24*XNK(1,2) - 4*XNK(1,1) - 21*XNK(2,1) - 21*XNK(3,1) - 
     >   4*XNK(4,1) + 24*XNK(4,2) + 24*XNK(4,3) - 4*XNK(4,4) - 
     >   21*XNK(3,4) - 21*XNK(2,4) + 126*XNK(2,3) + 126*XNK(2,2)))/
     >   1024000 
      Q2(10,65) = Q2(10,65) 
     >   - (15**(0.5D0)*ETAH*SIGN(1.0,ETAH)*(692*XNJ(1,4,K2) + 
     >   3024*XNJ(2,4,K2) - 3024*XNJ(3,4,K2) - 692*XNJ(4,4,K2) + 
     >   3633*XNJ(1,3,K2) + 15876*XNJ(2,3,K2) - 15876*XNJ(3,3,K2) - 
     >   3633*XNJ(4,3,K2) + 692*XNJ(1,1,K2) + 3024*XNJ(2,1,K2) - 
     >   3024*XNJ(3,1,K2) - 692*XNJ(4,1,K2) + 3633*XNJ(1,2,K2) + 
     >   15876*XNJ(2,2,K2) - 15876*XNJ(3,2,K2) - 3633*XNJ(4,2,K2)))/
     >   15360000 - (15**(0.5D0)*MUH*(693*XNI(4,3,J2,K2) - 
     >   4158*XNI(3,3,J2,K2) - 4158*XNI(2,3,J2,K2) + 693*
     >   XNI(1,3,J2,K2) + 
     >   132*XNI(4,4,J2,K2) - 792*XNI(3,4,J2,K2) - 792*XNI(2,4,J2,K2) +
     >   132*XNI(1,4,J2,K2) + 132*XNI(4,1,J2,K2) - 792*XNI(3,1,J2,K2) -
     >   792*XNI(2,1,J2,K2) + 132*XNI(1,1,J2,K2) + 693*XNI(4,2,J2,K2) -
     >   4158*XNI(3,2,J2,K2) - 4158*XNI(2,2,J2,K2) + 693*
     >   XNI(1,2,J2,K2)))/
     >   15360000 - (15**(0.5D0)*XIH*SIGN(1.0,XIH)*(4536*
     >   XNK(3,2) + 4536*XNK(3,3) + 173*XNK(1,4) - 1038*XNK(1,3) - 
     >   1038*XNK(1,2) + 173*XNK(1,1) + 756*XNK(2,1) - 756*XNK(3,1) - 
     >   173*XNK(4,1) + 1038*XNK(4,2) + 1038*XNK(4,3) - 173*XNK(4,4) - 
     >   756*XNK(3,4) + 756*XNK(2,4) - 4536*XNK(2,3) - 4536*
     >   XNK(2,2)))/15360000 
      Q2(11,65) = Q2(11,65) + 
     >   ((189*XNI(4,3,J2,K2))/1024000 - (567*XNI(3,3,J2,K2))/512000 -
     >   (567*XNI(2,3,J2,K2))/512000 + (189*XNI(1,3,J2,K2))/1024000 + 
     >   (9*XNI(4,4,J2,K2))/256000 - (27*XNI(3,4,J2,K2))/128000 - 
     >   (27*XNI(2,4,J2,K2))/128000 + (9*XNI(1,4,J2,K2))/256000 + 
     >   (9*XNI(4,1,J2,K2))/256000 - (27*XNI(3,1,J2,K2))/128000 - 
     >   (27*XNI(2,1,J2,K2))/128000 + (9*XNI(1,1,J2,K2))/256000 + 
     >   (189*XNI(4,2,J2,K2))/1024000 - (567*XNI(3,2,J2,K2))/512000 - 
     >   (567*XNI(2,2,J2,K2))/512000 + (189*XNI(1,2,J2,K2))/1024000)*
     >   MUH*SIGN(1.0,MUH) + ((9*XNJ(1,4,K2))/256000 - (27*
     >   XNJ(2,4,K2))/128000 - (27*XNJ(3,4,K2))/128000 + (9*
     >   XNJ(4,4,K2))/256000 + (189*XNJ(1,3,K2))/1024000 - (567*
     >   XNJ(2,3,K2))/512000 - (567*XNJ(3,3,K2))/512000 + (189*
     >   XNJ(4,3,K2))/1024000 + (9*XNJ(1,1,K2))/256000 - (27*
     >   XNJ(2,1,K2))/128000 - (27*XNJ(3,1,K2))/128000 + (9*
     >   XNJ(4,1,K2))/256000 + (189*XNJ(1,2,K2))/1024000 - (567*
     >   XNJ(2,2,K2))/512000 - (567*XNJ(3,2,K2))/512000 + (189*
     >   XNJ(4,2,K2))/1024000)*ETAH*SIGN(1.0,ETAH) + ((81*
     >   XNK(3,2))/256000 + (81*XNK(3,3))/256000 + (9*XNK(1,4))/
     >   1024000 - (27*XNK(1,3))/512000 - (27*XNK(1,2))/512000 + 
     >   (9*XNK(1,1))/1024000 - (27*XNK(2,1))/512000 - (27*XNK(3,1))/
     >   512000 + (9*XNK(4,1))/1024000 - (27*XNK(4,2))/512000 - (27*
     >   XNK(4,3))/512000 + (9*XNK(4,4))/1024000 - (27*XNK(3,4))/
     >   512000 - (27*XNK(2,4))/512000 + (81*XNK(2,3))/256000 + 
     >   (81*XNK(2,2))/256000)*XIH*SIGN(1.0,XIH) 
      Q2(12,65) = Q2(12,65) + 
     >   (9*35**(0.5D0)*ETAH*SIGN(1.0,ETAH)*(296*XNJ(1,4,K2) + 
     >   2052*XNJ(2,4,K2) - 2052*XNJ(3,4,K2) - 296*XNJ(4,4,K2) + 
     >   1554*XNJ(1,3,K2) + 10773*XNJ(2,3,K2) - 10773*XNJ(3,3,K2) - 
     >   1554*XNJ(4,3,K2) + 296*XNJ(1,1,K2) + 2052*XNJ(2,1,K2) - 
     >   2052*XNJ(3,1,K2) - 296*XNJ(4,1,K2) + 1554*XNJ(1,2,K2) + 
     >   10773*XNJ(2,2,K2) - 10773*XNJ(3,2,K2) - 1554*XNJ(4,2,K2)))/
     >   250880000 - (9*35**(0.5D0)*MUH*(441*XNI(4,3,J2,K2) - 
     >   2646*XNI(3,3,J2,K2) - 2646*XNI(2,3,J2,K2) + 441*
     >   XNI(1,3,J2,K2) + 
     >   84*XNI(4,4,J2,K2) - 504*XNI(3,4,J2,K2) - 504*XNI(2,4,J2,K2) + 
     >   84*XNI(1,4,J2,K2) + 84*XNI(4,1,J2,K2) - 504*XNI(3,1,J2,K2) - 
     >   504*XNI(2,1,J2,K2) + 84*XNI(1,1,J2,K2) + 441*XNI(4,2,J2,K2) -
     >   2646*XNI(3,2,J2,K2) - 2646*XNI(2,2,J2,K2) + 441*
     >   XNI(1,2,J2,K2)))/
     >   250880000 + (9*35**(0.5D0)*XIH*SIGN(1.0,XIH)*(3078*
     >   XNK(3,2) + 3078*XNK(3,3) + 74*XNK(1,4) - 444*XNK(1,3) - 
     >   444*XNK(1,2) + 74*XNK(1,1) + 513*XNK(2,1) - 513*XNK(3,1) - 
     >   74*XNK(4,1) + 444*XNK(4,2) + 444*XNK(4,3) - 74*XNK(4,4) - 
     >   513*XNK(3,4) + 513*XNK(2,4) - 3078*XNK(2,3) - 3078*
     >   XNK(2,2)))/250880000 
      Q2(13,65) = Q2(13,65) + 
     >   (3*7**(0.5D0)*XIH*SIGN(1.0,XIH)*(10773*XNK(3,2) - 
     >   10773*XNK(3,3) - 296*XNK(1,4) - 2052*XNK(1,3) + 2052*
     >   XNK(1,2) + 296*XNK(1,1) + 1554*XNK(2,1) + 1554*XNK(3,1) + 
     >   296*XNK(4,1) + 2052*XNK(4,2) - 2052*XNK(4,3) - 296*
     >   XNK(4,4) - 1554*XNK(3,4) - 1554*XNK(2,4) - 10773*XNK(2,3) + 
     >   10773*XNK(2,2)))/50176000 - (3*7**(0.5D0)*MUH*
     >   SIGN(1.0,MUH)*(1554*XNI(4,3,J2,K2) + 10773*XNI(3,3,J2,K2) - 
     >   10773*XNI(2,3,J2,K2) - 1554*XNI(1,3,J2,K2) + 296*
     >   XNI(4,4,J2,K2)+
     >   2052*XNI(3,4,J2,K2) - 2052*XNI(2,4,J2,K2) - 296*XNI(1,4,J2,K2)+
     >   296*XNI(4,1,J2,K2) + 2052*XNI(3,1,J2,K2) - 2052*XNI(2,1,J2,K2)-
     >   296*XNI(1,1,J2,K2) + 1554*XNI(4,2,J2,K2) + 10773*
     >   XNI(3,2,J2,K2) - 
     >   10773*XNI(2,2,J2,K2) - 1554*XNI(1,2,J2,K2)))/50176000 - 
     >   (3*7**(0.5D0)*ETAH*(336*XNJ(1,4,K2) + 1764*XNJ(2,4,K2) + 
     >   1764*XNJ(3,4,K2) + 336*XNJ(4,4,K2) + 1764*XNJ(1,3,K2) + 
     >   9261*XNJ(2,3,K2) + 9261*XNJ(3,3,K2) + 1764*XNJ(4,3,K2) + 
     >   336*XNJ(1,1,K2) + 1764*XNJ(2,1,K2) + 1764*XNJ(3,1,K2) + 
     >   336*XNJ(4,1,K2) + 1764*XNJ(1,2,K2) + 9261*XNJ(2,2,K2) + 
     >   9261*XNJ(3,2,K2) + 1764*XNJ(4,2,K2)))/50176000 
      Q2(14,65) = Q2(14,65) + 
     >   (21**(0.5D0)*ETAH*(14532*XNJ(1,4,K2) + 63504*XNJ(2,4,K2) - 
     >   63504*XNJ(3,4,K2) - 14532*XNJ(4,4,K2) + 76293*XNJ(1,3,K2) + 
     >   333396*XNJ(2,3,K2) - 333396*XNJ(3,3,K2) - 76293*XNJ(4,3,K2)+
     >   14532*XNJ(1,1,K2) + 63504*XNJ(2,1,K2) - 63504*XNJ(3,1,K2) - 
     >   14532*XNJ(4,1,K2) + 76293*XNJ(1,2,K2) + 333396*XNJ(2,2,K2) - 
     >   333396*XNJ(3,2,K2) - 76293*XNJ(4,2,K2)))/752640000 + 
     >   (21**(0.5D0)*MUH*(51282*XNI(4,3,J2,K2) + 
     >   355509*XNI(3,3,J2,K2) - 355509*XNI(2,3,J2,K2) - 51282*
     >   XNI(1,3,J2,K2) + 9768*XNI(4,4,J2,K2) + 67716*XNI(3,4,J2,K2) - 
     >   67716*XNI(2,4,J2,K2) - 9768*XNI(1,4,J2,K2) + 9768*
     >   XNI(4,1,J2,K2) + 
     >   67716*XNI(3,1,J2,K2) - 67716*XNI(2,1,J2,K2) - 9768*
     >   XNI(1,1,J2,K2) + 51282*XNI(4,2,J2,K2) + 355509*XNI(3,2,J2,K2)-
     >   355509*XNI(2,2,J2,K2) - 51282*XNI(1,2,J2,K2)))/752640000 + 
     >   (21**(0.5D0)*XIH*SIGN(1.0,XIH)*(387828*XNK(3,2) - 
     >   387828*XNK(3,3) + 12802*XNK(1,4) + 88749*XNK(1,3) - 
     >   88749*XNK(1,2) - 12802*XNK(1,1) - 55944*XNK(2,1) + 
     >   55944*XNK(3,1) + 12802*XNK(4,1) + 88749*XNK(4,2) - 
     >   88749*XNK(4,3) - 12802*XNK(4,4) - 55944*XNK(3,4) + 
     >   55944*XNK(2,4) + 387828*XNK(2,3) - 387828*XNK(2,2)))/
     >   752640000 
      Q2(15,65) = Q2(15,65) 
     >   - (9*35**(0.5D0)*ETAH*(84*XNJ(1,4,K2) - 504*
     >   XNJ(2,4,K2) - 504*XNJ(3,4,K2) + 84*XNJ(4,4,K2) + 441*
     >   XNJ(1,3,K2) - 2646*XNJ(2,3,K2) - 2646*XNJ(3,3,K2) + 441*
     >   XNJ(4,3,K2) + 84*XNJ(1,1,K2) - 504*XNJ(2,1,K2) - 504*
     >   XNJ(3,1,K2) + 84*XNJ(4,1,K2) + 441*XNJ(1,2,K2) - 2646*
     >   XNJ(2,2,K2) - 2646*XNJ(3,2,K2) + 441*XNJ(4,2,K2)))/250880000 - 
     >   (9*35**(0.5D0)*MUH*SIGN(1.0,MUH)*(1554*
     >   XNI(4,3,J2,K2) + 10773*XNI(3,3,J2,K2) - 10773*XNI(2,3,J2,K2) - 
     >   1554*XNI(1,3,J2,K2) + 296*XNI(4,4,J2,K2) + 2052*XNI(3,4,J2,K2)-
     >   2052*XNI(2,4,J2,K2) - 296*XNI(1,4,J2,K2) + 296*XNI(4,1,J2,K2) + 
     >   2052*XNI(3,1,J2,K2) - 2052*XNI(2,1,J2,K2) - 296*XNI(1,1,J2,K2)+
     >   1554*XNI(4,2,J2,K2) + 10773*XNI(3,2,J2,K2) - 10773*
     >   XNI(2,2,J2,K2) - 1554*XNI(1,2,J2,K2)))/250880000 - 
     >   (9*35**(0.5D0)*XIH*SIGN(1.0,XIH)*(3078*XNK(3,2) - 
     >   3078*XNK(3,3) + 74*XNK(1,4) + 513*XNK(1,3) - 513*XNK(1,2) - 
     >   74*XNK(1,1) + 444*XNK(2,1) + 444*XNK(3,1) - 74*XNK(4,1) - 
     >   513*XNK(4,2) + 513*XNK(4,3) + 74*XNK(4,4) - 444*XNK(3,4) - 
     >   444*XNK(2,4) - 3078*XNK(2,3) + 3078*XNK(2,2)))/250880000 
      Q2(16,65) = Q2(16,65) + 
     >   ((2997*XNI(4,3,J2,K2))/17920000 + (41553*XNI(3,3,J2,K2))/
     >   35840000 - (41553*XNI(2,3,J2,K2))/35840000 - (2997*
     >   XNI(1,3,J2,K2))/17920000 + (999*XNI(4,4,J2,K2))/31360000 + 
     >   (13851*XNI(3,4,J2,K2))/62720000 - (13851*XNI(2,4,J2,K2))/
     >   62720000 - (999*XNI(1,4,J2,K2))/31360000 + (999*
     >   XNI(4,1,J2,K2))/31360000 + (13851*XNI(3,1,J2,K2))/62720000 - 
     >   (13851*XNI(2,1,J2,K2))/62720000 - (999*XNI(1,1,J2,K2))/
     >   31360000 + (2997*XNI(4,2,J2,K2))/17920000 + (41553*
     >   XNI(3,2,J2,K2))/35840000 - (41553*XNI(2,2,J2,K2))/35840000 - 
     >   (2997*XNI(1,2,J2,K2))/17920000)*MUH + ((13851*
     >   XNJ(3,4,K2))/62720000 - (13851*XNJ(2,4,K2))/62720000 - 
     >   (999*XNJ(1,4,K2))/31360000 + (999*XNJ(4,4,K2))/31360000 - 
     >   (2997*XNJ(1,3,K2))/17920000 - (41553*XNJ(2,3,K2))/35840000 + 
     >   (41553*XNJ(3,3,K2))/35840000 + (2997*XNJ(4,3,K2))/17920000 - 
     >   (999*XNJ(1,1,K2))/31360000 - (13851*XNJ(2,1,K2))/62720000 + 
     >   (13851*XNJ(3,1,K2))/62720000 + (999*XNJ(4,1,K2))/31360000 - 
     >   (2997*XNJ(1,2,K2))/17920000 - (41553*XNJ(2,2,K2))/35840000 + 
     >   (41553*XNJ(3,2,K2))/35840000 + (2997*XNJ(4,2,K2))/17920000)*
     >   ETAH + ((2368521*XNK(3,3))/1756160000 - (2368521*
     >   XNK(3,2))/1756160000 - (12321*XNK(1,4))/439040000 - 
     >   (170829*XNK(1,3))/878080000 + (170829*XNK(1,2))/878080000 + 
     >   (12321*XNK(1,1))/439040000 + (170829*XNK(2,1))/878080000 - 
     >   (170829*XNK(3,1))/878080000 - (12321*XNK(4,1))/439040000 - 
     >   (170829*XNK(4,2))/878080000 + (170829*XNK(4,3))/878080000 + 
     >   (12321*XNK(4,4))/439040000 + (170829*XNK(3,4))/878080000 - 
     >   (170829*XNK(2,4))/878080000 - (2368521*XNK(2,3))/1756160000 + 
     >   (2368521*XNK(2,2))/1756160000)*XIH*SIGN(1.0,XIH) 
      Q2(17,65) = Q2(17,65) + 
     >   (3**(0.5D0)*ETAH*SIGN(1.0,ETAH)*(692*XNJ(1,4,K2) + 
     >   3633*XNJ(2,4,K2) + 3633*XNJ(3,4,K2) + 692*XNJ(4,4,K2) + 
     >   3024*XNJ(1,3,K2) + 15876*XNJ(2,3,K2) + 15876*XNJ(3,3,K2) + 
     >   3024*XNJ(4,3,K2) - 692*XNJ(1,1,K2) - 3633*XNJ(2,1,K2) - 
     >   3633*XNJ(3,1,K2) - 692*XNJ(4,1,K2) - 3024*XNJ(1,2,K2) - 
     >   15876*XNJ(2,2,K2) - 15876*XNJ(3,2,K2) - 3024*XNJ(4,2,K2)))/
     >   9216000 + (3**(0.5D0)*MUH*SIGN(1.0,MUH)*(3024*
     >   XNI(4,3,J2,K2) + 15876*XNI(3,3,J2,K2) + 15876*XNI(2,3,J2,K2) +
     >   3024*XNI(1,3,J2,K2) + 692*XNI(4,4,J2,K2) + 3633*XNI(3,4,J2,K2)+ 
     >   3633*XNI(2,4,J2,K2) + 692*XNI(1,4,J2,K2) - 692*XNI(4,1,J2,K2) -
     >   3633*XNI(3,1,J2,K2) - 3633*XNI(2,1,J2,K2) - 692*XNI(1,1,J2,K2)-
     >   3024*XNI(4,2,J2,K2) - 15876*XNI(3,2,J2,K2) - 15876*
     >   XNI(2,2,J2,K2) - 
     >   3024*XNI(1,2,J2,K2)))/9216000 - (3**(0.5D0)*XIH*(14553*
     >   XNK(3,2) + 14553*XNK(3,3) + 528*XNK(1,4) + 2772*XNK(1,3) + 
     >   2772*XNK(1,2) + 528*XNK(1,1) + 2772*XNK(2,1) + 2772*
     >   XNK(3,1) + 528*XNK(4,1) + 2772*XNK(4,2) + 2772*XNK(4,3) + 
     >   528*XNK(4,4) + 2772*XNK(3,4) + 2772*XNK(2,4) + 14553*
     >   XNK(2,3) + 14553*XNK(2,2)))/9216000 
      Q2(18,65) = Q2(18,65) + 
     >   ((1903*XNI(4,1,J2,K2))/11520000 - (4851*XNI(3,3,J2,K2))/
     >   1280000 -
     >   (4851*XNI(2,3,J2,K2))/1280000 - (231*XNI(1,3,J2,K2))/320000 - 
     >   (1903*XNI(4,4,J2,K2))/11520000 - (13321*XNI(3,4,J2,K2))/
     >   15360000 - (13321*XNI(2,4,J2,K2))/15360000 - (1903*
     >   XNI(1,4,J2,K2))/11520000 - (231*XNI(4,3,J2,K2))/320000 + 
     >   (13321*XNI(3,1,J2,K2))/15360000 + (13321*XNI(2,1,J2,K2))/
     >   15360000 + (1903*XNI(1,1,J2,K2))/11520000 + (231*
     >   XNI(4,2,J2,K2))/320000 + (4851*XNI(3,2,J2,K2))/1280000 + 
     >   (4851*XNI(2,2,J2,K2))/1280000 + (231*XNI(1,2,J2,K2))/320000)*
     >   MUH + ((1211*XNJ(3,4,K2))/1280000 - (1211*XNJ(2,4,K2))/
     >   1280000 - (29929*XNJ(1,4,K2))/138240000 + (29929*XNJ(4,4,K2))/
     >   138240000 - (1211*XNJ(1,3,K2))/1280000 - (1323*XNJ(2,3,K2))/
     >   320000 + (1323*XNJ(3,3,K2))/320000 + (1211*XNJ(4,3,K2))/
     >   1280000 + (29929*XNJ(1,1,K2))/138240000 + (1211*XNJ(2,1,K2))/
     >   1280000 - (1211*XNJ(3,1,K2))/1280000 - (29929*XNJ(4,1,K2))/
     >   138240000 + (1211*XNJ(1,2,K2))/1280000 + (1323*XNJ(2,2,K2))/
     >   320000 - (1323*XNJ(3,2,K2))/320000 - (1211*XNJ(4,2,K2))/
     >   1280000)*ETAH*SIGN(1.0,ETAH) + ((1903*XNK(1,4))/
     >   11520000 - (4851*XNK(3,3))/1280000 - (4851*XNK(3,2))/
     >   1280000 + (13321*XNK(1,3))/15360000 + (13321*XNK(1,2))/
     >   15360000 + (1903*XNK(1,1))/11520000 + (231*XNK(2,1))/
     >   320000 - (231*XNK(3,1))/320000 - (1903*XNK(4,1))/11520000 - 
     >   (13321*XNK(4,2))/15360000 - (13321*XNK(4,3))/15360000 - 
     >   (1903*XNK(4,4))/11520000 - (231*XNK(3,4))/320000 + 
     >   (231*XNK(2,4))/320000 + (4851*XNK(2,3))/1280000 + 
     >   (4851*XNK(2,2))/1280000)*XIH 
      Q2(19,65) = Q2(19,65) + 
     >   (15**(0.5D0)*ETAH*SIGN(1.0,ETAH)*(173*XNJ(1,4,K2) - 
     >   1038*XNJ(2,4,K2) - 1038*XNJ(3,4,K2) + 173*XNJ(4,4,K2) + 
     >   756*XNJ(1,3,K2) - 4536*XNJ(2,3,K2) - 4536*XNJ(3,3,K2) + 
     >   756*XNJ(4,3,K2) - 173*XNJ(1,1,K2) + 1038*XNJ(2,1,K2) + 
     >   1038*XNJ(3,1,K2) - 173*XNJ(4,1,K2) - 756*XNJ(1,2,K2) + 
     >   4536*XNJ(2,2,K2) + 4536*XNJ(3,2,K2) - 756*XNJ(4,2,K2)))/
     >   15360000 + (15**(0.5D0)*MUH*SIGN(1.0,MUH)*(3024*
     >   XNI(4,3,J2,K2) + 15876*XNI(3,3,J2,K2) + 15876*XNI(2,3,J2,K2) + 
     >   3024*XNI(1,3,J2,K2) + 692*XNI(4,4,J2,K2) + 3633*XNI(3,4,J2,K2)+
     >   3633*XNI(2,4,J2,K2) + 692*XNI(1,4,J2,K2) - 692*XNI(4,1,J2,K2)- 
     >   3633*XNI(3,1,J2,K2) - 3633*XNI(2,1,J2,K2) - 692*XNI(1,1,J2,K2)-
     >   3024*XNI(4,2,J2,K2) - 15876*XNI(3,2,J2,K2) - 15876*
     >   XNI(2,2,J2,K2) - 
     >   3024*XNI(1,2,J2,K2)))/15360000 + (15**(0.5D0)*XIH*(4158*
     >   XNK(3,2) + 4158*XNK(3,3) - 132*XNK(1,4) - 693*XNK(1,3) - 
     >   693*XNK(1,2) - 132*XNK(1,1) + 792*XNK(2,1) + 792*XNK(3,1) - 
     >   132*XNK(4,1) - 693*XNK(4,2) - 693*XNK(4,3) - 132*XNK(4,4) + 
     >   792*XNK(3,4) + 792*XNK(2,4) + 4158*XNK(2,3) + 4158*
     >   XNK(2,2)))/15360000 
      Q2(20,65) = Q2(20,65) + 
     >   (21**(0.5D0)*ETAH*SIGN(1.0,ETAH)*(12802*XNJ(1,4,K2) + 
     >   88749*XNJ(2,4,K2) - 88749*XNJ(3,4,K2) - 12802*XNJ(4,4,K2) + 
     >   55944*XNJ(1,3,K2) + 387828*XNJ(2,3,K2) - 387828*XNJ(3,3,K2) -
     >   55944*XNJ(4,3,K2) - 12802*XNJ(1,1,K2) - 88749*XNJ(2,1,K2) + 
     >   88749*XNJ(3,1,K2) + 12802*XNJ(4,1,K2) - 55944*XNJ(1,2,K2) - 
     >   387828*XNJ(2,2,K2) + 387828*XNJ(3,2,K2) + 55944*XNJ(4,2,K2)))/
     >   752640000 - (21**(0.5D0)*MUH*(63504*XNI(4,3,J2,K2) + 
     >   333396*XNI(3,3,J2,K2) + 333396*XNI(2,3,J2,K2) + 63504*
     >   XNI(1,3,J2,K2) + 14532*XNI(4,4,J2,K2) + 76293*XNI(3,4,J2,K2) +
     >   76293*XNI(2,4,J2,K2) + 14532*XNI(1,4,J2,K2) - 14532*
     >   XNI(4,1,J2,K2) - 76293*XNI(3,1,J2,K2) - 76293*XNI(2,1,J2,K2) -
     >   14532*XNI(1,1,J2,K2) - 63504*XNI(4,2,J2,K2) - 333396*
     >   XNI(3,2,J2,K2) - 333396*XNI(2,2,J2,K2) - 63504*
     >   XNI(1,2,J2,K2)))/
     >   752640000 + (21**(0.5D0)*XIH*(355509*XNK(3,2) + 
     >   355509*XNK(3,3) - 9768*XNK(1,4) - 51282*XNK(1,3) - 
     >   51282*XNK(1,2) - 9768*XNK(1,1) - 67716*XNK(2,1) + 
     >   67716*XNK(3,1) + 9768*XNK(4,1) + 51282*XNK(4,2) + 
     >   51282*XNK(4,3) + 9768*XNK(4,4) + 67716*XNK(3,4) - 
     >   67716*XNK(2,4) - 355509*XNK(2,3) - 355509*XNK(2,2)))/
     >   752640000 
      Q2(21,65) = Q2(21,65) + 
     >   ((1211*XNI(4,3,J2,K2))/1280000 + (1323*XNI(3,3,J2,K2))/320000-
     >   (1323*XNI(2,3,J2,K2))/320000 - (1211*XNI(1,3,J2,K2))/1280000 + 
     >   (29929*XNI(4,4,J2,K2))/138240000 + (1211*XNI(3,4,J2,K2))/
     >   1280000 - (1211*XNI(2,4,J2,K2))/1280000 - (29929*
     >   XNI(1,4,J2,K2))/138240000 - (29929*XNI(4,1,J2,K2))/138240000 -
     >   (1211*XNI(3,1,J2,K2))/1280000 + (1211*XNI(2,1,J2,K2))/1280000+
     >   (29929*XNI(1,1,J2,K2))/138240000 - (1211*XNI(4,2,J2,K2))/
     >   1280000 - (1323*XNI(3,2,J2,K2))/320000 + (1323*XNI(2,2,J2,K2))/
     >   320000 + (1211*XNI(1,2,J2,K2))/1280000)*MUH*
     >   SIGN(1.0,MUH) + ((1903*XNJ(1,1,K2))/11520000 - (13321*
     >   XNJ(2,4,K2))/15360000 - (13321*XNJ(3,4,K2))/15360000 - 
     >   (1903*XNJ(4,4,K2))/11520000 - (231*XNJ(1,3,K2))/320000 - 
     >   (4851*XNJ(2,3,K2))/1280000 - (4851*XNJ(3,3,K2))/1280000 - 
     >   (231*XNJ(4,3,K2))/320000 - (1903*XNJ(1,4,K2))/11520000 + 
     >   (13321*XNJ(2,1,K2))/15360000 + (13321*XNJ(3,1,K2))/15360000 + 
     >   (1903*XNJ(4,1,K2))/11520000 + (231*XNJ(1,2,K2))/320000 + 
     >   (4851*XNJ(2,2,K2))/1280000 + (4851*XNJ(3,2,K2))/1280000 + 
     >   (231*XNJ(4,2,K2))/320000)*ETAH + ((4851*XNK(3,2))/
     >   1280000 - (4851*XNK(3,3))/1280000 - (1903*XNK(1,4))/
     >   11520000 - (231*XNK(1,3))/320000 + (231*XNK(1,2))/320000 + 
     >   (1903*XNK(1,1))/11520000 + (13321*XNK(2,1))/15360000 + 
     >   (13321*XNK(3,1))/15360000 + (1903*XNK(4,1))/11520000 + 
     >   (231*XNK(4,2))/320000 - (231*XNK(4,3))/320000 - 
     >   (1903*XNK(4,4))/11520000 - (13321*XNK(3,4))/15360000 - 
     >   (13321*XNK(2,4))/15360000 - (4851*XNK(2,3))/1280000 + 
     >   (4851*XNK(2,2))/1280000)*XIH 
      Q2(22,65) = Q2(22,65) + 
     >   (11*3**(0.5D0)*ETAH*(29929*XNJ(1,4,K2) + 130788*
     >   XNJ(2,4,K2) - 130788*XNJ(3,4,K2) - 29929*XNJ(4,4,K2) + 
     >   130788*XNJ(1,3,K2) + 571536*XNJ(2,3,K2) - 571536*XNJ(3,3,K2) -
     >   130788*XNJ(4,3,K2) - 29929*XNJ(1,1,K2) - 130788*XNJ(2,1,K2) + 
     >   130788*XNJ(3,1,K2) + 29929*XNJ(4,1,K2) - 130788*XNJ(1,2,K2) - 
     >   571536*XNJ(2,2,K2) + 571536*XNJ(3,2,K2) + 130788*XNJ(4,2,K2)))/
     >   2073600000 - (11*3**(0.5D0)*MUH*(130788*XNI(4,3,J2,K2) + 
     >   571536*XNI(3,3,J2,K2) - 571536*XNI(2,3,J2,K2) - 130788*
     >   XNI(1,3,J2,K2) + 29929*XNI(4,4,J2,K2) + 130788*XNI(3,4,J2,K2)-
     >   130788*XNI(2,4,J2,K2) - 29929*XNI(1,4,J2,K2) - 29929*
     >   XNI(4,1,J2,K2) - 130788*XNI(3,1,J2,K2) + 130788*XNI(2,1,J2,K2)+
     >   29929*XNI(1,1,J2,K2) - 130788*XNI(4,2,J2,K2) - 571536*
     >   XNI(3,2,J2,K2) + 571536*XNI(2,2,J2,K2) + 130788*
     >   XNI(1,2,J2,K2)))/
     >   2073600000 + (11*3**(0.5D0)*XIH*(571536*XNK(3,2) - 
     >   571536*XNK(3,3) + 29929*XNK(1,4) + 130788*XNK(1,3) - 
     >   130788*XNK(1,2) - 29929*XNK(1,1) - 130788*XNK(2,1) + 
     >   130788*XNK(3,1) + 29929*XNK(4,1) + 130788*XNK(4,2) - 
     >   130788*XNK(4,3) - 29929*XNK(4,4) - 130788*XNK(3,4) + 
     >   130788*XNK(2,4) + 571536*XNK(2,3) - 571536*XNK(2,2)))/
     >   2073600000 
      Q2(23,65) = Q2(23,65) + 
     >   (5**(0.5D0)*MUH*SIGN(1.0,MUH)*(130788*XNI(4,3,J2,K2) + 
     >   571536*XNI(3,3,J2,K2) - 571536*XNI(2,3,J2,K2) - 130788*
     >   XNI(1,3,J2,K2) + 29929*XNI(4,4,J2,K2) + 130788*XNI(3,4,J2,K2) - 
     >   130788*XNI(2,4,J2,K2) - 29929*XNI(1,4,J2,K2) - 29929*
     >   XNI(4,1,J2,K2) - 130788*XNI(3,1,J2,K2) + 130788*XNI(2,1,J2,K2)+
     >   29929*XNI(1,1,J2,K2) - 130788*XNI(4,2,J2,K2) - 571536*
     >   XNI(3,2,J2,K2) + 571536*XNI(2,2,J2,K2) + 130788*
     >   XNI(1,2,J2,K2)))/230400000 - (5**(0.5D0)*ETAH*(5709*
     >   XNJ(1,4,K2) - 34254*XNJ(2,4,K2) - 34254*XNJ(3,4,K2) + 5709*
     >   XNJ(4,4,K2) + 24948*XNJ(1,3,K2) - 149688*XNJ(2,3,K2) - 149688*
     >   XNJ(3,3,K2) + 24948*XNJ(4,3,K2) - 5709*XNJ(1,1,K2) + 34254*
     >   XNJ(2,1,K2) + 34254*XNJ(3,1,K2) - 5709*XNJ(4,1,K2) - 24948*
     >   XNJ(1,2,K2) + 149688*XNJ(2,2,K2) + 149688*XNJ(3,2,K2) - 24948*
     >   XNJ(4,2,K2)))/230400000 - (5**(0.5D0)*XIH*(149688*
     >   XNK(3,2) - 149688*XNK(3,3) + 5709*XNK(1,4) + 24948*XNK(1,3) - 
     >   24948*XNK(1,2) - 5709*XNK(1,1) + 34254*XNK(2,1) + 34254*
     >   XNK(3,1) - 5709*XNK(4,1) - 24948*XNK(4,2) + 24948*
     >   XNK(4,3) + 5709*XNK(4,4) - 34254*XNK(3,4) - 34254*XNK(2,4) - 
     >   149688*XNK(2,3) + 149688*XNK(2,2)))/230400000 
      Q2(24,65) = Q2(24,65) 
     >   - (7**(0.5D0)*ETAH*(140822*XNJ(1,4,K2) + 976239*
     >   XNJ(2,4,K2) - 976239*XNJ(3,4,K2) - 140822*XNJ(4,4,K2) + 
     >   615384*XNJ(1,3,K2) + 4266108*XNJ(2,3,K2) - 4266108*
     >   XNJ(3,3,K2) - 615384*XNJ(4,3,K2) - 140822*XNJ(1,1,K2) - 
     >   976239*XNJ(2,1,K2) + 976239*XNJ(3,1,K2) + 140822*XNJ(4,1,K2) - 
     >   615384*XNJ(1,2,K2) - 4266108*XNJ(2,2,K2) + 4266108*XNJ(3,2,K2)+
     >   615384*XNJ(4,2,K2)))/CONST9 - (7**(0.5D0)*MUH*
     >   (915516*XNI(4,3,J2,K2) + 4000752*XNI(3,3,J2,K2) - 4000752*
     >   XNI(2,3,J2,K2) - 915516*XNI(1,3,J2,K2) + 209503*XNI(4,4,J2,K2)+
     >   915516*XNI(3,4,J2,K2) - 915516*XNI(2,4,J2,K2) - 209503*
     >   XNI(1,4,J2,K2) - 209503*XNI(4,1,J2,K2) - 915516*XNI(3,1,J2,K2)+
     >   915516*XNI(2,1,J2,K2) + 209503*XNI(1,1,J2,K2) - 915516*
     >   XNI(4,2,J2,K2) - 4000752*XNI(3,2,J2,K2) + 4000752*
     >   XNI(2,2,J2,K2) + 
     >   915516*XNI(1,2,J2,K2)))/CONST9 - (7**(0.5D0)*XIH*
     >   (4266108*XNK(3,2) - 4266108*XNK(3,3) + 140822*XNK(1,4) + 
     >   615384*XNK(1,3) - 615384*XNK(1,2) - 140822*XNK(1,1) - 
     >   976239*XNK(2,1) + 976239*XNK(3,1) + 140822*XNK(4,1) + 
     >   615384*XNK(4,2) - 615384*XNK(4,3) - 140822*XNK(4,4) - 
     >   976239*XNK(3,4) + 976239*XNK(2,4) + 4266108*XNK(2,3) - 
     >   4266108*XNK(2,2)))/CONST9 
      Q2(25,65) = Q2(25,65) + 
     >   (15**(0.5D0)*ETAH*SIGN(1.0,ETAH)*(692*XNJ(1,4,K2) + 
     >   3633*XNJ(2,4,K2) + 3633*XNJ(3,4,K2) + 692*XNJ(4,4,K2) + 
     >   3024*XNJ(1,3,K2) + 15876*XNJ(2,3,K2) + 15876*XNJ(3,3,K2) + 
     >   3024*XNJ(4,3,K2) - 692*XNJ(1,1,K2) - 3633*XNJ(2,1,K2) - 
     >   3633*XNJ(3,1,K2) - 692*XNJ(4,1,K2) - 3024*XNJ(1,2,K2) - 
     >   15876*XNJ(2,2,K2) - 15876*XNJ(3,2,K2) - 3024*XNJ(4,2,K2)))/
     >   15360000 + (15**(0.5D0)*MUH*SIGN(1.0,MUH)*(756*
     >   XNI(4,3,J2,K2) - 4536*XNI(3,3,J2,K2) - 4536*XNI(2,3,J2,K2) + 
     >   756*XNI(1,3,J2,K2) + 173*XNI(4,4,J2,K2) - 1038*XNI(3,4,J2,K2) -
     >   1038*XNI(2,4,J2,K2) + 173*XNI(1,4,J2,K2) - 173*XNI(4,1,J2,K2) +
     >   1038*XNI(3,1,J2,K2) + 1038*XNI(2,1,J2,K2) - 173*XNI(1,1,J2,K2)-
     >   756*XNI(4,2,J2,K2) + 4536*XNI(3,2,J2,K2) + 4536*XNI(2,2,J2,K2)-
     >   756*XNI(1,2,J2,K2)))/15360000 + (15**(0.5D0)*XIH*
     >   (4158*XNK(3,2) + 4158*XNK(3,3) - 132*XNK(1,4) + 
     >   792*XNK(1,3) + 792*XNK(1,2) - 132*XNK(1,1) - 693*XNK(2,1) - 
     >   693*XNK(3,1) - 132*XNK(4,1) + 792*XNK(4,2) + 792*XNK(4,3) - 
     >   132*XNK(4,4) - 693*XNK(3,4) - 693*XNK(2,4) + 4158*XNK(2,3) + 
     >   4158*XNK(2,2)))/15360000 
      Q2(26,65) = Q2(26,65) + 
     >   (5**(0.5D0)*XIH*(149688*XNK(3,2) + 149688*XNK(3,3) + 
     >   5709*XNK(1,4) - 34254*XNK(1,3) - 34254*XNK(1,2) + 5709*
     >   XNK(1,1) + 24948*XNK(2,1) - 24948*XNK(3,1) - 5709*XNK(4,1) + 
     >   34254*XNK(4,2) + 34254*XNK(4,3) - 5709*XNK(4,4) - 24948*
     >   XNK(3,4) + 24948*XNK(2,4) - 149688*XNK(2,3) - 149688*
     >   XNK(2,2)))/230400000 - (5**(0.5D0)*MUH*(24948*
     >   XNI(4,3,J2,K2) - 149688*XNI(3,3,J2,K2) - 149688*
     >   XNI(2,3,J2,K2) + 
     >   24948*XNI(1,3,J2,K2) + 5709*XNI(4,4,J2,K2) - 34254*
     >   XNI(3,4,J2,K2) - 
     >   34254*XNI(2,4,J2,K2) + 5709*XNI(1,4,J2,K2) - 5709*
     >   XNI(4,1,J2,K2) + 
     >   34254*XNI(3,1,J2,K2) + 34254*XNI(2,1,J2,K2) - 5709*
     >   XNI(1,1,J2,K2) - 
     >   24948*XNI(4,2,J2,K2) + 149688*XNI(3,2,J2,K2) + 149688*
     >   XNI(2,2,J2,K2) - 24948*XNI(1,2,J2,K2)))/230400000 - 
     >   (5**(0.5D0)*ETAH*SIGN(1.0,ETAH)*(29929*XNJ(1,4,K2) + 
     >   130788*XNJ(2,4,K2) - 130788*XNJ(3,4,K2) - 29929*XNJ(4,4,K2) + 
     >   130788*XNJ(1,3,K2) + 571536*XNJ(2,3,K2) - 571536*XNJ(3,3,K2) - 
     >   130788*XNJ(4,3,K2) - 29929*XNJ(1,1,K2) - 130788*XNJ(2,1,K2) + 
     >   130788*XNJ(3,1,K2) + 29929*XNJ(4,1,K2) - 130788*XNJ(1,2,K2) - 
     >   571536*XNJ(2,2,K2) + 571536*XNJ(3,2,K2) + 130788*XNJ(4,2,K2)))/
     >   230400000 
      Q2(27,65) = Q2(27,65) + 
     >   (3**(0.5D0)*ETAH*SIGN(1.0,ETAH)*(173*XNJ(1,4,K2) - 
     >   1038*XNJ(2,4,K2) - 1038*XNJ(3,4,K2) + 173*XNJ(4,4,K2) + 
     >   756*XNJ(1,3,K2) - 4536*XNJ(2,3,K2) - 4536*XNJ(3,3,K2) + 
     >   756*XNJ(4,3,K2) - 173*XNJ(1,1,K2) + 1038*XNJ(2,1,K2) + 
     >   1038*XNJ(3,1,K2) - 173*XNJ(4,1,K2) - 756*XNJ(1,2,K2) + 
     >   4536*XNJ(2,2,K2) + 4536*XNJ(3,2,K2) - 756*XNJ(4,2,K2)))/
     >   5120000 + (3**(0.5D0)*MUH*SIGN(1.0,MUH)*(756*
     >   XNI(4,3,J2,K2) - 4536*XNI(3,3,J2,K2) - 4536*XNI(2,3,J2,K2) + 
     >   756*XNI(1,3,J2,K2) + 173*XNI(4,4,J2,K2) - 1038*
     >   XNI(3,4,J2,K2) - 
     >   1038*XNI(2,4,J2,K2) + 173*XNI(1,4,J2,K2) - 173*
     >   XNI(4,1,J2,K2) + 
     >   1038*XNI(3,1,J2,K2) + 1038*XNI(2,1,J2,K2) - 173*
     >   XNI(1,1,J2,K2) - 
     >   756*XNI(4,2,J2,K2) + 4536*XNI(3,2,J2,K2) + 4536*
     >   XNI(2,2,J2,K2) - 
     >   756*XNI(1,2,J2,K2)))/5120000 - (3**(0.5D0)*XIH*(1188*
     >   XNK(3,2) + 1188*XNK(3,3) + 33*XNK(1,4) - 198*XNK(1,3) - 
     >   198*XNK(1,2) + 33*XNK(1,1) - 198*XNK(2,1) - 198*XNK(3,1) + 
     >   33*XNK(4,1) - 198*XNK(4,2) - 198*XNK(4,3) + 33*XNK(4,4) - 
     >   198*XNK(3,4) - 198*XNK(2,4) + 1188*XNK(2,3) + 1188*
     >   XNK(2,2)))/5120000 
      Q2(28,65) = Q2(28,65) + 
     >   (105**(0.5D0)*ETAH*SIGN(1.0,ETAH)*(12802*XNJ(1,4,K2) + 
     >   88749*XNJ(2,4,K2) - 88749*XNJ(3,4,K2) - 12802*XNJ(4,4,K2) + 
     >   55944*XNJ(1,3,K2) + 387828*XNJ(2,3,K2) - 387828*XNJ(3,3,K2) - 
     >   55944*XNJ(4,3,K2) - 12802*XNJ(1,1,K2) - 88749*XNJ(2,1,K2) + 
     >   88749*XNJ(3,1,K2) + 12802*XNJ(4,1,K2) - 55944*XNJ(1,2,K2) - 
     >   387828*XNJ(2,2,K2) + 387828*XNJ(3,2,K2) + 55944*XNJ(4,2,K2)))/
     >   1254400000 - (105**(0.5D0)*MUH*(15876*XNI(4,3,J2,K2) - 
     >   95256*XNI(3,3,J2,K2) - 95256*XNI(2,3,J2,K2) + 15876*
     >   XNI(1,3,J2,K2) + 3633*XNI(4,4,J2,K2) - 21798*XNI(3,4,J2,K2) - 
     >   21798*XNI(2,4,J2,K2) + 3633*XNI(1,4,J2,K2) - 3633*
     >   XNI(4,1,J2,K2) + 
     >   21798*XNI(3,1,J2,K2) + 21798*XNI(2,1,J2,K2) - 3633*
     >   XNI(1,1,J2,K2) - 
     >   15876*XNI(4,2,J2,K2) + 95256*XNI(3,2,J2,K2) + 95256*
     >   XNI(2,2,J2,K2) - 15876*XNI(1,2,J2,K2)))/1254400000 - 
     >   (105**(0.5D0)*XIH*(101574*XNK(3,2) + 101574*XNK(3,3) + 
     >   2442*XNK(1,4) - 14652*XNK(1,3) - 14652*XNK(1,2) + 
     >   2442*XNK(1,1) + 16929*XNK(2,1) - 16929*XNK(3,1) - 
     >   2442*XNK(4,1) + 14652*XNK(4,2) + 14652*XNK(4,3) - 
     >   2442*XNK(4,4) - 16929*XNK(3,4) + 16929*XNK(2,4) - 
     >   101574*XNK(2,3) - 101574*XNK(2,2)))/1254400000 
      Q2(29,65) = Q2(29,65) 
     >   - (21**(0.5D0)*ETAH*(14532*XNJ(1,4,K2) + 76293*
     >   XNJ(2,4,K2) + 76293*XNJ(3,4,K2) + 14532*XNJ(4,4,K2) + 
     >   63504*XNJ(1,3,K2) + 333396*XNJ(2,3,K2) + 333396*XNJ(3,3,K2) + 
     >   63504*XNJ(4,3,K2) - 14532*XNJ(1,1,K2) - 76293*XNJ(2,1,K2) - 
     >   76293*XNJ(3,1,K2) - 14532*XNJ(4,1,K2) - 63504*XNJ(1,2,K2) - 
     >   333396*XNJ(2,2,K2) - 333396*XNJ(3,2,K2) - 63504*XNJ(4,2,K2)))/
     >   752640000 - (21**(0.5D0)*MUH*SIGN(1.0,MUH)*(55944*
     >   XNI(4,3,J2,K2) + 387828*XNI(3,3,J2,K2) - 387828*
     >   XNI(2,3,J2,K2) - 
     >   55944*XNI(1,3,J2,K2) + 12802*XNI(4,4,J2,K2) + 88749*
     >   XNI(3,4,J2,K2) - 88749*XNI(2,4,J2,K2) - 12802*XNI(1,4,J2,K2) - 
     >   12802*XNI(4,1,J2,K2) - 88749*XNI(3,1,J2,K2) + 88749*
     >   XNI(2,1,J2,K2) + 12802*XNI(1,1,J2,K2) - 55944*XNI(4,2,J2,K2) - 
     >   387828*XNI(3,2,J2,K2) + 387828*XNI(2,2,J2,K2) + 55944*
     >   XNI(1,2,J2,K2)))/752640000 - (21**(0.5D0)*XIH*
     >   (355509*XNK(3,2) - 355509*XNK(3,3) - 9768*XNK(1,4) - 
     >   67716*XNK(1,3) + 67716*XNK(1,2) + 9768*XNK(1,1) + 
     >   51282*XNK(2,1) + 51282*XNK(3,1) + 9768*XNK(4,1) + 
     >   67716*XNK(4,2) - 67716*XNK(4,3) - 9768*XNK(4,4) - 
     >   51282*XNK(3,4) - 51282*XNK(2,4) - 355509*XNK(2,3) + 
     >   355509*XNK(2,2)))/752640000 
      Q2(30,65) = Q2(30,65) + 
     >   (7**(0.5D0)*ETAH*(209503*XNJ(1,4,K2) + 915516*
     >   XNJ(2,4,K2) - 915516*XNJ(3,4,K2) - 209503*XNJ(4,4,K2) + 
     >   915516*XNJ(1,3,K2) + 4000752*XNJ(2,3,K2) - 4000752*
     >   XNJ(3,3,K2) - 915516*XNJ(4,3,K2) - 209503*XNJ(1,1,K2) - 
     >   915516*XNJ(2,1,K2) + 915516*XNJ(3,1,K2) + 209503*XNJ(4,1,K2) -
     >   915516*XNJ(1,2,K2) - 4000752*XNJ(2,2,K2) + 4000752*
     >   XNJ(3,2,K2) + 
     >   915516*XNJ(4,2,K2)))/CONST9 + (7**(0.5D0)*MUH*
     >   (615384*XNI(4,3,J2,K2) + 4266108*XNI(3,3,J2,K2) - 4266108*
     >   XNI(2,3,J2,K2) - 615384*XNI(1,3,J2,K2) + 140822*
     >   XNI(4,4,J2,K2) + 
     >   976239*XNI(3,4,J2,K2) - 976239*XNI(2,4,J2,K2) - 140822*
     >   XNI(1,4,J2,K2) - 140822*XNI(4,1,J2,K2) - 976239*
     >   XNI(3,1,J2,K2) + 
     >   976239*XNI(2,1,J2,K2) + 140822*XNI(1,1,J2,K2) - 615384*
     >   XNI(4,2,J2,K2) - 4266108*XNI(3,2,J2,K2) + 4266108*
     >   XNI(2,2,J2,K2) + 
     >   615384*XNI(1,2,J2,K2)))/CONST9 - (7**(0.5D0)*XIH*
     >   (4266108*XNK(3,2) - 4266108*XNK(3,3) + 140822*XNK(1,4) + 
     >   976239*XNK(1,3) - 976239*XNK(1,2) - 140822*XNK(1,1) - 
     >   615384*XNK(2,1) + 615384*XNK(3,1) + 140822*XNK(4,1) + 
     >   976239*XNK(4,2) - 976239*XNK(4,3) - 140822*XNK(4,4) - 
     >   615384*XNK(3,4) + 615384*XNK(2,4) + 4266108*XNK(2,3) - 
     >   4266108*XNK(2,2)))/CONST9 
      Q2(31,65) = Q2(31,65) + 
     >   (105**(0.5D0)*XIH*(101574*XNK(3,2) - 101574*XNK(3,3) + 
     >   2442*XNK(1,4) + 16929*XNK(1,3) - 16929*XNK(1,2) - 2442*
     >   XNK(1,1) + 14652*XNK(2,1) + 14652*XNK(3,1) - 2442*XNK(4,1) - 
     >   16929*XNK(4,2) + 16929*XNK(4,3) + 2442*XNK(4,4) - 14652*
     >   XNK(3,4) - 14652*XNK(2,4) - 101574*XNK(2,3) + 101574*
     >   XNK(2,2)))/1254400000 - (105**(0.5D0)*MUH*
     >   SIGN(1.0,MUH)*(55944*XNI(4,3,J2,K2) + 387828*XNI(3,3,J2,K2) - 
     >   387828*XNI(2,3,J2,K2) - 55944*XNI(1,3,J2,K2) + 12802*
     >   XNI(4,4,J2,K2) + 88749*XNI(3,4,J2,K2) - 88749*XNI(2,4,J2,K2) - 
     >   12802*XNI(1,4,J2,K2) - 12802*XNI(4,1,J2,K2) - 88749*
     >   XNI(3,1,J2,K2) + 88749*XNI(2,1,J2,K2) + 12802*XNI(1,1,J2,K2) - 
     >   55944*XNI(4,2,J2,K2) - 387828*XNI(3,2,J2,K2) + 387828*
     >   XNI(2,2,J2,K2) + 55944*XNI(1,2,J2,K2)))/1254400000 - 
     >   (105**(0.5D0)*ETAH*(3633*XNJ(1,4,K2) - 21798*
     >   XNJ(2,4,K2) - 21798*XNJ(3,4,K2) + 3633*XNJ(4,4,K2) + 15876*
     >   XNJ(1,3,K2) - 95256*XNJ(2,3,K2) - 95256*XNJ(3,3,K2) + 15876*
     >   XNJ(4,3,K2) - 3633*XNJ(1,1,K2) + 21798*XNJ(2,1,K2) + 21798*
     >   XNJ(3,1,K2) - 3633*XNJ(4,1,K2) - 15876*XNJ(1,2,K2) + 95256*
     >   XNJ(2,2,K2) + 95256*XNJ(3,2,K2) - 15876*XNJ(4,2,K2)))/
     >   1254400000 
      Q2(32,65) = Q2(32,65) + 
     >   (3*3**(0.5D0)*MUH*(391608*XNI(4,3,J2,K2) + 2714796*
     >   XNI(3,3,J2,K2) - 2714796*XNI(2,3,J2,K2) - 391608*
     >   XNI(1,3,J2,K2) + 
     >   89614*XNI(4,4,J2,K2) + 621243*XNI(3,4,J2,K2) - 621243*
     >   XNI(2,4,J2,K2) - 89614*XNI(1,4,J2,K2) - 89614*XNI(4,1,J2,K2) - 
     >   621243*XNI(3,1,J2,K2) + 621243*XNI(2,1,J2,K2) + 89614*
     >   XNI(1,1,J2,K2) - 391608*XNI(4,2,J2,K2) - 2714796*
     >   XNI(3,2,J2,K2) + 
     >   2714796*XNI(2,2,J2,K2) + 391608*XNI(1,2,J2,K2)))/CONST10 - 
     >   (3*3**(0.5D0)*ETAH*(89614*XNJ(1,4,K2) + 621243*
     >   XNJ(2,4,K2) - 621243*XNJ(3,4,K2) - 89614*XNJ(4,4,K2) + 
     >   391608*XNJ(1,3,K2) + 2714796*XNJ(2,3,K2) - 2714796*
     >   XNJ(3,3,K2) - 391608*XNJ(4,3,K2) - 89614*XNJ(1,1,K2) - 
     >   621243*XNJ(2,1,K2) + 621243*XNJ(3,1,K2) + 89614*XNJ(4,1,K2) - 
     >   391608*XNJ(1,2,K2) - 2714796*XNJ(2,2,K2) + 2714796*
     >   XNJ(3,2,K2) + 391608*XNJ(4,2,K2)))/CONST10 + 
     >   (3*3**(0.5D0)*XIH*(2894859*XNK(3,2) - 2894859*
     >   XNK(3,3) + 60236*XNK(1,4) + 417582*XNK(1,3) - 417582*
     >   XNK(1,2) - 60236*XNK(1,1) - 417582*XNK(2,1) + 417582*
     >   XNK(3,1) + 60236*XNK(4,1) + 417582*XNK(4,2) - 417582*
     >   XNK(4,3) - 60236*XNK(4,4) - 417582*XNK(3,4) + 417582*
     >   XNK(2,4) + 2894859*XNK(2,3) - 2894859*XNK(2,2)))/CONST10 
      Q2(33,65) = Q2(33,65) + 
     >   (3*5**(0.5D0)*ETAH*SIGN(1.0,ETAH)*(4*XNJ(1,4,K2) + 
     >   21*XNJ(2,4,K2) + 21*XNJ(3,4,K2) + 4*XNJ(4,4,K2) - 
     >   24*XNJ(1,3,K2) - 126*XNJ(2,3,K2) - 126*XNJ(3,3,K2) - 
     >   24*XNJ(4,3,K2) + 4*XNJ(1,1,K2) + 21*XNJ(2,1,K2) + 
     >   21*XNJ(3,1,K2) + 4*XNJ(4,1,K2) - 24*XNJ(1,2,K2) - 
     >   126*XNJ(2,2,K2) - 126*XNJ(3,2,K2) - 24*XNJ(4,2,K2)))/1024000 - 
     >   (3*5**(0.5D0)*MUH*SIGN(1.0,MUH)*(24*XNI(4,3,J2,K2) + 
     >   126*XNI(3,3,J2,K2) + 126*XNI(2,3,J2,K2) + 24*XNI(1,3,J2,K2) - 
     >   4*XNI(4,4,J2,K2) - 21*XNI(3,4,J2,K2) - 21*XNI(2,4,J2,K2) - 
     >   4*XNI(1,4,J2,K2) - 4*XNI(4,1,J2,K2) - 21*XNI(3,1,J2,K2) - 
     >   21*XNI(2,1,J2,K2) - 4*XNI(1,1,J2,K2) + 24*XNI(4,2,J2,K2) + 
     >   126*XNI(3,2,J2,K2) + 126*XNI(2,2,J2,K2) + 24*XNI(1,2,J2,K2)))/
     >   1024000 + (3*5**(0.5D0)*XIH*SIGN(1.0,XIH)*(441*
     >   XNK(3,2) + 441*XNK(3,3) + 16*XNK(1,4) + 84*XNK(1,3) + 
     >   84*XNK(1,2) + 16*XNK(1,1) + 84*XNK(2,1) + 84*XNK(3,1) + 
     >   16*XNK(4,1) + 84*XNK(4,2) + 84*XNK(4,3) + 16*XNK(4,4) + 
     >   84*XNK(3,4) + 84*XNK(2,4) + 441*XNK(2,3) + 441*XNK(2,2)))/
     >   1024000 
      Q2(34,65) = Q2(34,65) + 
     >   (15**(0.5D0)*MUH*(792*XNI(4,3,J2,K2) + 4158*
     >   XNI(3,3,J2,K2) + 4158*XNI(2,3,J2,K2) + 792*XNI(1,3,J2,K2) - 
     >   132*XNI(4,4,J2,K2) - 693*XNI(3,4,J2,K2) - 693*XNI(2,4,J2,K2) - 
     >   132*XNI(1,4,J2,K2) - 132*XNI(4,1,J2,K2) - 693*XNI(3,1,J2,K2) - 
     >   693*XNI(2,1,J2,K2) - 132*XNI(1,1,J2,K2) + 792*XNI(4,2,J2,K2) + 
     >   4158*XNI(3,2,J2,K2) + 4158*XNI(2,2,J2,K2) + 792*
     >   XNI(1,2,J2,K2)))/
     >   15360000 - (15**(0.5D0)*ETAH*SIGN(1.0,ETAH)*(173*
     >   XNJ(1,4,K2) + 756*XNJ(2,4,K2) - 756*XNJ(3,4,K2) - 173*
     >   XNJ(4,4,K2) - 1038*XNJ(1,3,K2) - 4536*XNJ(2,3,K2) + 4536*
     >   XNJ(3,3,K2) + 1038*XNJ(4,3,K2) + 173*XNJ(1,1,K2) + 756*
     >   XNJ(2,1,K2) - 756*XNJ(3,1,K2) - 173*XNJ(4,1,K2) - 1038*
     >   XNJ(1,2,K2) - 4536*XNJ(2,2,K2) + 4536*XNJ(3,2,K2) + 1038*
     >   XNJ(4,2,K2)))/15360000 + (15**(0.5D0)*XIH*
     >   SIGN(1.0,XIH)*(15876*XNK(3,2) + 15876*XNK(3,3) - 
     >   692*XNK(1,4) - 3633*XNK(1,3) - 3633*XNK(1,2) - 
     >   692*XNK(1,1) - 3024*XNK(2,1) + 3024*XNK(3,1) + 
     >   692*XNK(4,1) + 3633*XNK(4,2) + 3633*XNK(4,3) + 
     >   692*XNK(4,4) + 3024*XNK(3,4) - 3024*XNK(2,4) - 
     >   15876*XNK(2,3) - 15876*XNK(2,2)))/15360000 
      Q2(35,65) = Q2(35,65) + 
     >   ((9*XNI(4,4,J2,K2))/256000 - (567*XNI(3,3,J2,K2))/512000 - 
     >   (567*XNI(2,3,J2,K2))/512000 - (27*XNI(1,3,J2,K2))/128000 - 
     >   (27*XNI(4,3,J2,K2))/128000 + (189*XNI(3,4,J2,K2))/1024000 + 
     >   (189*XNI(2,4,J2,K2))/1024000 + (9*XNI(1,4,J2,K2))/256000 + 
     >   (9*XNI(4,1,J2,K2))/256000 + (189*XNI(3,1,J2,K2))/1024000 + 
     >   (189*XNI(2,1,J2,K2))/1024000 + (9*XNI(1,1,J2,K2))/256000 - 
     >   (27*XNI(4,2,J2,K2))/128000 - (567*XNI(3,2,J2,K2))/512000 - 
     >   (567*XNI(2,2,J2,K2))/512000 - (27*XNI(1,2,J2,K2))/128000)*
     >   MUH*SIGN(1.0,MUH) + ((9*XNJ(1,4,K2))/1024000 - 
     >   (27*XNJ(2,4,K2))/512000 - (27*XNJ(3,4,K2))/512000 + (9*
     >   XNJ(4,4,K2))/1024000 - (27*XNJ(1,3,K2))/512000 + 
     >   (81*XNJ(2,3,K2))/256000 + (81*XNJ(3,3,K2))/256000 - 
     >   (27*XNJ(4,3,K2))/512000 + (9*XNJ(1,1,K2))/1024000 - 
     >   (27*XNJ(2,1,K2))/512000 - (27*XNJ(3,1,K2))/512000 + 
     >   (9*XNJ(4,1,K2))/1024000 - (27*XNJ(1,2,K2))/512000 + 
     >   (81*XNJ(2,2,K2))/256000 + (81*XNJ(3,2,K2))/256000 - 
     >   (27*XNJ(4,2,K2))/512000)*ETAH*SIGN(1.0,ETAH) + 
     >   ((9*XNK(1,4))/256000 - (567*XNK(3,3))/512000 - (567*
     >   XNK(3,2))/512000 + (189*XNK(1,3))/1024000 + (189*XNK(1,2))/
     >   1024000 + (9*XNK(1,1))/256000 - (27*XNK(2,1))/128000 - 
     >   (27*XNK(3,1))/128000 + (9*XNK(4,1))/256000 + (189*XNK(4,2))/
     >   1024000 + (189*XNK(4,3))/1024000 + (9*XNK(4,4))/256000 - 
     >   (27*XNK(3,4))/128000 - (27*XNK(2,4))/128000 - (567*
     >   XNK(2,3))/512000 - (567*XNK(2,2))/512000)*XIH*
     >   SIGN(1.0,XIH) 
      Q2(36,65) = Q2(36,65) + 
     >   (9*35**(0.5D0)*ETAH*SIGN(1.0,ETAH)*(74*XNJ(1,4,K2) + 
     >   513*XNJ(2,4,K2) - 513*XNJ(3,4,K2) - 74*XNJ(4,4,K2) - 444*
     >   XNJ(1,3,K2) - 3078*XNJ(2,3,K2) + 3078*XNJ(3,3,K2) + 444*
     >   XNJ(4,3,K2) + 74*XNJ(1,1,K2) + 513*XNJ(2,1,K2) - 513*
     >   XNJ(3,1,K2) - 74*XNJ(4,1,K2) - 444*XNJ(1,2,K2) - 3078*
     >   XNJ(2,2,K2) + 3078*XNJ(3,2,K2) + 444*XNJ(4,2,K2)))/250880000 +
     >   (9*35**(0.5D0)*MUH*(504*XNI(4,3,J2,K2) + 2646*
     >   XNI(3,3,J2,K2) + 2646*XNI(2,3,J2,K2) + 504*XNI(1,3,J2,K2) - 
     >   84*XNI(4,4,J2,K2) - 441*XNI(3,4,J2,K2) - 441*XNI(2,4,J2,K2) - 
     >   84*XNI(1,4,J2,K2) - 84*XNI(4,1,J2,K2) - 441*XNI(3,1,J2,K2) - 
     >   441*XNI(2,1,J2,K2) - 84*XNI(1,1,J2,K2) + 504*XNI(4,2,J2,K2) + 
     >   2646*XNI(3,2,J2,K2) + 2646*XNI(2,2,J2,K2) + 504*
     >   XNI(1,2,J2,K2)))/
     >   250880000 - (9*35**(0.5D0)*XIH*SIGN(1.0,XIH)*
     >   (10773*XNK(3,2) + 10773*XNK(3,3) - 296*XNK(1,4) - 1554*
     >   XNK(1,3) - 1554*XNK(1,2) - 296*XNK(1,1) - 2052*XNK(2,1) + 
     >   2052*XNK(3,1) + 296*XNK(4,1) + 1554*XNK(4,2) + 1554*
     >   XNK(4,3) + 296*XNK(4,4) + 2052*XNK(3,4) - 2052*XNK(2,4) - 
     >   10773*XNK(2,3) - 10773*XNK(2,2)))/250880000 
      Q2(37,65) = Q2(37,65) 
     >   - (15**(0.5D0)*ETAH*(132*XNJ(1,4,K2) + 693*XNJ(2,4,K2) + 
     >   693*XNJ(3,4,K2) + 132*XNJ(4,4,K2) - 792*XNJ(1,3,K2) - 4158*
     >   XNJ(2,3,K2) - 4158*XNJ(3,3,K2) - 792*XNJ(4,3,K2) + 132*
     >   XNJ(1,1,K2) + 693*XNJ(2,1,K2) + 693*XNJ(3,1,K2) + 132*
     >   XNJ(4,1,K2) - 792*XNJ(1,2,K2) - 4158*XNJ(2,2,K2) - 4158*
     >   XNJ(3,2,K2) - 792*XNJ(4,2,K2)))/15360000 - (15**(0.5D0)*
     >   MUH*SIGN(1.0,MUH)*(1038*XNI(4,3,J2,K2) + 4536*
     >   XNI(3,3,J2,K2) - 4536*XNI(2,3,J2,K2) - 1038*XNI(1,3,J2,K2) - 
     >   173*XNI(4,4,J2,K2) - 756*XNI(3,4,J2,K2) + 756*XNI(2,4,J2,K2)+
     >   173*XNI(1,4,J2,K2) - 173*XNI(4,1,J2,K2) - 756*XNI(3,1,J2,K2)+
     >   756*XNI(2,1,J2,K2) + 173*XNI(1,1,J2,K2) + 1038*XNI(4,2,J2,K2)+
     >   4536*XNI(3,2,J2,K2) - 4536*XNI(2,2,J2,K2) - 1038*
     >   XNI(1,2,J2,K2)))/
     >   15360000 - (15**(0.5D0)*XIH*SIGN(1.0,XIH)*(15876*
     >   XNK(3,2) - 15876*XNK(3,3) - 692*XNK(1,4) - 3024*XNK(1,3) + 
     >   3024*XNK(1,2) + 692*XNK(1,1) + 3633*XNK(2,1) + 3633*
     >   XNK(3,1) + 692*XNK(4,1) + 3024*XNK(4,2) - 3024*XNK(4,3) - 
     >   692*XNK(4,4) - 3633*XNK(3,4) - 3633*XNK(2,4) - 15876*
     >   XNK(2,3) + 15876*XNK(2,2)))/15360000 
      Q2(38,65) = Q2(38,65) + 
     >   (5**(0.5D0)*ETAH*(5709*XNJ(1,4,K2) + 24948*XNJ(2,4,K2) - 
     >   24948*XNJ(3,4,K2) - 5709*XNJ(4,4,K2) - 34254*XNJ(1,3,K2) - 
     >   149688*XNJ(2,3,K2) + 149688*XNJ(3,3,K2) + 34254*XNJ(4,3,K2) + 
     >   5709*XNJ(1,1,K2) + 24948*XNJ(2,1,K2) - 24948*XNJ(3,1,K2) - 
     >   5709*XNJ(4,1,K2) - 34254*XNJ(1,2,K2) - 149688*XNJ(2,2,K2) + 
     >   149688*XNJ(3,2,K2) + 34254*XNJ(4,2,K2)))/230400000 + 
     >   (5**(0.5D0)*MUH*(34254*XNI(4,3,J2,K2) + 149688*
     >   XNI(3,3,J2,K2) - 149688*XNI(2,3,J2,K2) - 34254*XNI(1,3,J2,K2)-
     >   5709*XNI(4,4,J2,K2) - 24948*XNI(3,4,J2,K2) + 24948*
     >   XNI(2,4,J2,K2) + 5709*XNI(1,4,J2,K2) - 5709*XNI(4,1,J2,K2) - 
     >   24948*XNI(3,1,J2,K2) + 24948*XNI(2,1,J2,K2) + 5709*
     >   XNI(1,1,J2,K2) + 34254*XNI(4,2,J2,K2) + 149688*XNI(3,2,J2,K2)-
     >   149688*XNI(2,2,J2,K2) - 34254*XNI(1,2,J2,K2)))/230400000 - 
     >   (5**(0.5D0)*XIH*SIGN(1.0,XIH)*(571536*XNK(3,2) - 
     >   571536*XNK(3,3) + 29929*XNK(1,4) + 130788*XNK(1,3) - 
     >   130788*XNK(1,2) - 29929*XNK(1,1) - 130788*XNK(2,1) + 
     >   130788*XNK(3,1) + 29929*XNK(4,1) + 130788*XNK(4,2) - 
     >   130788*XNK(4,3) - 29929*XNK(4,4) - 130788*XNK(3,4) + 
     >   130788*XNK(2,4) + 571536*XNK(2,3) - 571536*XNK(2,2)))/
     >   230400000 
      Q2(39,65) = Q2(39,65) + 
     >   (3**(0.5D0)*XIH*SIGN(1.0,XIH)*(4536*XNK(3,2) - 
     >   4536*XNK(3,3) + 173*XNK(1,4) + 756*XNK(1,3) - 756*XNK(1,2) - 
     >   173*XNK(1,1) + 1038*XNK(2,1) + 1038*XNK(3,1) - 173*
     >   XNK(4,1) - 756*XNK(4,2) + 756*XNK(4,3) + 173*XNK(4,4) - 
     >   1038*XNK(3,4) - 1038*XNK(2,4) - 4536*XNK(2,3) + 4536*
     >   XNK(2,2)))/5120000 - (3**(0.5D0)*MUH*SIGN(1.0,MUH)*
     >   (1038*XNI(4,3,J2,K2) + 4536*XNI(3,3,J2,K2) - 4536*
     >   XNI(2,3,J2,K2) - 
     >   1038*XNI(1,3,J2,K2) - 173*XNI(4,4,J2,K2) - 756*XNI(3,4,J2,K2)+
     >   756*XNI(2,4,J2,K2) + 173*XNI(1,4,J2,K2) - 173*XNI(4,1,J2,K2) - 
     >   756*XNI(3,1,J2,K2) + 756*XNI(2,1,J2,K2) + 173*XNI(1,1,J2,K2) + 
     >   1038*XNI(4,2,J2,K2) + 4536*XNI(3,2,J2,K2) - 4536*
     >   XNI(2,2,J2,K2) - 1038*XNI(1,2,J2,K2)))/5120000 - (3**(0.5D0)*
     >   ETAH*(33*XNJ(1,4,K2) - 198*XNJ(2,4,K2) - 198*XNJ(3,4,K2) + 
     >   33*XNJ(4,4,K2) - 198*XNJ(1,3,K2) + 1188*XNJ(2,3,K2) + 1188*
     >   XNJ(3,3,K2) - 198*XNJ(4,3,K2) + 33*XNJ(1,1,K2) - 198*
     >   XNJ(2,1,K2) - 198*XNJ(3,1,K2) + 33*XNJ(4,1,K2) - 198*
     >   XNJ(1,2,K2) + 1188*XNJ(2,2,K2) + 1188*XNJ(3,2,K2) - 198*
     >   XNJ(4,2,K2)))/5120000 
      Q2(40,65) = Q2(40,65) + 
     >   (105**(0.5D0)*MUH*(21798*XNI(4,3,J2,K2) + 95256*
     >   XNI(3,3,J2,K2) - 95256*XNI(2,3,J2,K2) - 21798*XNI(1,3,J2,K2) - 
     >   3633*XNI(4,4,J2,K2) - 15876*XNI(3,4,J2,K2) + 15876*
     >   XNI(2,4,J2,K2) + 3633*XNI(1,4,J2,K2) - 3633*XNI(4,1,J2,K2) - 
     >   15876*XNI(3,1,J2,K2) + 15876*XNI(2,1,J2,K2) + 3633*
     >   XNI(1,1,J2,K2) + 21798*XNI(4,2,J2,K2) + 95256*XNI(3,2,J2,K2) - 
     >   95256*XNI(2,2,J2,K2) - 21798*XNI(1,2,J2,K2)))/1254400000 - 
     >   (105**(0.5D0)*ETAH*(2442*XNJ(1,4,K2) + 16929*
     >   XNJ(2,4,K2) - 16929*XNJ(3,4,K2) - 2442*XNJ(4,4,K2) - 14652*
     >   XNJ(1,3,K2) - 101574*XNJ(2,3,K2) + 101574*XNJ(3,3,K2) + 
     >   14652*XNJ(4,3,K2) + 2442*XNJ(1,1,K2) + 16929*XNJ(2,1,K2) - 
     >   16929*XNJ(3,1,K2) - 2442*XNJ(4,1,K2) - 14652*XNJ(1,2,K2) - 
     >   101574*XNJ(2,2,K2) + 101574*XNJ(3,2,K2) + 14652*XNJ(4,2,K2)))/
     >   1254400000 + (105**(0.5D0)*XIH*SIGN(1.0,XIH)*
     >   (387828*XNK(3,2) - 387828*XNK(3,3) + 12802*XNK(1,4) + 
     >   55944*XNK(1,3) - 55944*XNK(1,2) - 12802*XNK(1,1) - 
     >   88749*XNK(2,1) + 88749*XNK(3,1) + 12802*XNK(4,1) + 
     >   55944*XNK(4,2) - 55944*XNK(4,3) - 12802*XNK(4,4) - 
     >   88749*XNK(3,4) + 88749*XNK(2,4) + 387828*XNK(2,3) - 
     >   387828*XNK(2,2)))/1254400000 
      Q2(41,65) = Q2(41,65) + 
     >   ((81*XNI(3,3,J2,K2))/256000 - (27*XNI(4,3,J2,K2))/512000 + 
     >   (81*XNI(2,3,J2,K2))/256000 - (27*XNI(1,3,J2,K2))/512000 + 
     >   (9*XNI(4,4,J2,K2))/1024000 - (27*XNI(3,4,J2,K2))/512000 - 
     >   (27*XNI(2,4,J2,K2))/512000 + (9*XNI(1,4,J2,K2))/1024000 + 
     >   (9*XNI(4,1,J2,K2))/1024000 - (27*XNI(3,1,J2,K2))/512000 - 
     >   (27*XNI(2,1,J2,K2))/512000 + (9*XNI(1,1,J2,K2))/1024000 - 
     >   (27*XNI(4,2,J2,K2))/512000 + (81*XNI(3,2,J2,K2))/256000 + 
     >   (81*XNI(2,2,J2,K2))/256000 - (27*XNI(1,2,J2,K2))/512000)*
     >   MUH*SIGN(1.0,MUH) + ((9*XNJ(1,4,K2))/256000 + 
     >  (189*XNJ(2,4,K2))/1024000 + (189*XNJ(3,4,K2))/1024000 + 
     >   (9*XNJ(4,4,K2))/256000 - (27*XNJ(1,3,K2))/128000 - (567*
     >   XNJ(2,3,K2))/512000 - (567*XNJ(3,3,K2))/512000 - (27*
     >   XNJ(4,3,K2))/128000 + (9*XNJ(1,1,K2))/256000 + (189*
     >   XNJ(2,1,K2))/1024000 + (189*XNJ(3,1,K2))/1024000 + (9*
     >   XNJ(4,1,K2))/256000 - (27*XNJ(1,2,K2))/128000 - (567*
     >   XNJ(2,2,K2))/512000 - (567*XNJ(3,2,K2))/512000 - (27*
     >   XNJ(4,2,K2))/128000)*ETAH*SIGN(1.0,ETAH) + ((9*
     >   XNK(1,4))/256000 - (567*XNK(3,3))/512000 - (567*
     >   XNK(3,2))/512000 - (27*XNK(1,3))/128000 - (27*
     >   XNK(1,2))/128000 + (9*XNK(1,1))/256000 + (189*
     >   XNK(2,1))/1024000 + (189*XNK(3,1))/1024000 + (9*
     >   XNK(4,1))/256000 - (27*XNK(4,2))/128000 - (27*
     >   XNK(4,3))/128000 + (9*XNK(4,4))/256000 + (189*
     >   XNK(3,4))/1024000 + (189*XNK(2,4))/1024000 - (567*
     >   XNK(2,3))/512000 - (567*XNK(2,2))/512000)*XIH*
     >   SIGN(1.0,XIH) 
      Q2(42,65) = Q2(42,65) + 
     >   (3**(0.5D0)*MUH*(198*XNI(4,3,J2,K2) - 1188*XNI(3,3,J2,K2) - 
     >   1188*XNI(2,3,J2,K2) + 198*XNI(1,3,J2,K2) - 33*XNI(4,4,J2,K2) + 
     >   198*XNI(3,4,J2,K2) + 198*XNI(2,4,J2,K2) - 33*XNI(1,4,J2,K2) - 
     >   33*XNI(4,1,J2,K2) + 198*XNI(3,1,J2,K2) + 198*XNI(2,1,J2,K2) - 
     >   33*XNI(1,1,J2,K2) + 198*XNI(4,2,J2,K2) - 1188*XNI(3,2,J2,K2) - 
     >   1188*XNI(2,2,J2,K2) + 198*XNI(1,2,J2,K2)))/5120000 - 
     >   (3**(0.5D0)*ETAH*SIGN(1.0,ETAH)*(173*XNJ(1,4,K2) + 
     >   756*XNJ(2,4,K2) - 756*XNJ(3,4,K2) - 173*XNJ(4,4,K2) - 
     >   1038*XNJ(1,3,K2) - 4536*XNJ(2,3,K2) + 4536*XNJ(3,3,K2) + 
     >   1038*XNJ(4,3,K2) + 173*XNJ(1,1,K2) + 756*XNJ(2,1,K2) - 
     >   756*XNJ(3,1,K2) - 173*XNJ(4,1,K2) - 1038*XNJ(1,2,K2) - 
     >   4536*XNJ(2,2,K2) + 4536*XNJ(3,2,K2) + 1038*XNJ(4,2,K2)))/
     >   5120000 - (3**(0.5D0)*XIH*SIGN(1.0,XIH)*(4536*
     >   XNK(3,2) + 4536*XNK(3,3) + 173*XNK(1,4) - 1038*XNK(1,3) - 
     >   1038*XNK(1,2) + 173*XNK(1,1) + 756*XNK(2,1) - 756*XNK(3,1) - 
     >   173*XNK(4,1) + 1038*XNK(4,2) + 1038*XNK(4,3) - 173*
     >   XNK(4,4) - 756*XNK(3,4) + 756*XNK(2,4) - 4536*XNK(2,3) - 
     >   4536*XNK(2,2)))/5120000 
      Q2(43,65) = Q2(43,65) + 
     >   (27*5**(0.5D0)*ETAH*SIGN(1.0,ETAH)*(XNJ(1,4,K2) - 
     >   6*XNJ(2,4,K2) - 6*XNJ(3,4,K2) + XNJ(4,4,K2) - 6*XNJ(1,3,K2) + 
     >   36*XNJ(2,3,K2) + 36*XNJ(3,3,K2) - 6*XNJ(4,3,K2) + XNJ(1,1,K2) - 
     >   6*XNJ(2,1,K2) - 6*XNJ(3,1,K2) + XNJ(4,1,K2) - 6*XNJ(1,2,K2) + 
     >   36*XNJ(2,2,K2) + 36*XNJ(3,2,K2) - 6*XNJ(4,2,K2)))/5120000 - 
     >   (27*5**(0.5D0)*MUH*SIGN(1.0,MUH)*(6*XNI(4,3,J2,K2) - 
     >   36*XNI(3,3,J2,K2) - 36*XNI(2,3,J2,K2) + 6*XNI(1,3,J2,K2) - 
     >   XNI(4,4,J2,K2) + 6*XNI(3,4,J2,K2) + 6*XNI(2,4,J2,K2) - 
     >   XNI(1,4,J2,K2) - XNI(4,1,J2,K2) + 6*XNI(3,1,J2,K2) + 
     >   6*XNI(2,1,J2,K2) - XNI(1,1,J2,K2) + 6*XNI(4,2,J2,K2) - 
     >   36*XNI(3,2,J2,K2) - 36*XNI(2,2,J2,K2) + 6*XNI(1,2,J2,K2)))/
     >   5120000 + (27*5**(0.5D0)*XIH*SIGN(1.0,XIH)*
     >   (36*XNK(3,2) + 36*XNK(3,3) + XNK(1,4) - 6*XNK(1,3) - 
     >   6*XNK(1,2) + XNK(1,1) - 6*XNK(2,1) - 6*XNK(3,1) + 
     >   XNK(4,1) - 6*XNK(4,2) - 6*XNK(4,3) + XNK(4,4) - 
     >   6*XNK(3,4) - 6*XNK(2,4) + 36*XNK(2,3) + 36*XNK(2,2)))/
     >   5120000 
      Q2(44,65) = Q2(44,65) + 
     >   (27*7**(0.5D0)*ETAH*SIGN(1.0,ETAH)*(74*XNJ(1,4,K2) + 
     >   513*XNJ(2,4,K2) - 513*XNJ(3,4,K2) - 74*XNJ(4,4,K2) - 
     >   444*XNJ(1,3,K2) - 3078*XNJ(2,3,K2) + 3078*XNJ(3,3,K2) + 
     >   444*XNJ(4,3,K2) + 74*XNJ(1,1,K2) + 513*XNJ(2,1,K2) - 513*
     >   XNJ(3,1,K2) - 74*XNJ(4,1,K2) - 444*XNJ(1,2,K2) - 3078*
     >   XNJ(2,2,K2) + 3078*XNJ(3,2,K2) + 444*XNJ(4,2,K2)))/250880000 +
     >   (27*7**(0.5D0)*MUH*(126*XNI(4,3,J2,K2) - 756*
     >   XNI(3,3,J2,K2) - 756*XNI(2,3,J2,K2) + 126*XNI(1,3,J2,K2) - 
     >   21*XNI(4,4,J2,K2) + 126*XNI(3,4,J2,K2) + 126*XNI(2,4,J2,K2) - 
     >   21*XNI(1,4,J2,K2) - 21*XNI(4,1,J2,K2) + 126*XNI(3,1,J2,K2) + 
     >   126*XNI(2,1,J2,K2) - 21*XNI(1,1,J2,K2) + 126*XNI(4,2,J2,K2) - 
     >   756*XNI(3,2,J2,K2) - 756*XNI(2,2,J2,K2) + 126*XNI(1,2,J2,K2)))/
     >   250880000 + (27*7**(0.5D0)*XIH*SIGN(1.0,XIH)*(3078*
     >   XNK(3,2) + 3078*XNK(3,3) + 74*XNK(1,4) - 444*XNK(1,3) - 
     >   444*XNK(1,2) + 74*XNK(1,1) + 513*XNK(2,1) - 513*XNK(3,1) - 
     >   74*XNK(4,1) + 444*XNK(4,2) + 444*XNK(4,3) - 74*XNK(4,4) - 
     >   513*XNK(3,4) + 513*XNK(2,4) - 3078*XNK(2,3) - 3078*
     >   XNK(2,2)))/250880000 
      Q2(45,65) = Q2(45,65) + 
     >   (9*35**(0.5D0)*MUH*SIGN(1.0,MUH)*(444*XNI(4,3,J2,K2) + 
     >   3078*XNI(3,3,J2,K2) - 3078*XNI(2,3,J2,K2) - 444*XNI(1,3,J2,K2)-
     >   74*XNI(4,4,J2,K2) - 513*XNI(3,4,J2,K2) + 513*XNI(2,4,J2,K2) + 
     >   74*XNI(1,4,J2,K2) - 74*XNI(4,1,J2,K2) - 513*XNI(3,1,J2,K2) + 
     >   513*XNI(2,1,J2,K2) + 74*XNI(1,1,J2,K2) + 444*XNI(4,2,J2,K2) + 
     >   3078*XNI(3,2,J2,K2) - 3078*XNI(2,2,J2,K2) - 444*
     >   XNI(1,2,J2,K2)))/
     >   250880000 - (9*35**(0.5D0)*ETAH*(84*XNJ(1,4,K2) + 441*
     >   XNJ(2,4,K2) + 441*XNJ(3,4,K2) + 84*XNJ(4,4,K2) - 504*
     >   XNJ(1,3,K2) - 2646*XNJ(2,3,K2) - 2646*XNJ(3,3,K2) - 504*
     >   XNJ(4,3,K2) + 84*XNJ(1,1,K2) + 441*XNJ(2,1,K2) + 441*
     >   XNJ(3,1,K2) + 84*XNJ(4,1,K2) - 504*XNJ(1,2,K2) - 2646*
     >   XNJ(2,2,K2) - 2646*XNJ(3,2,K2) - 504*XNJ(4,2,K2)))/
     >   250880000 + (9*35**(0.5D0)*XIH*SIGN(1.0,XIH)*
     >   (10773*XNK(3,2) - 10773*XNK(3,3) - 296*XNK(1,4) - 2052*
     >   XNK(1,3) + 2052*XNK(1,2) + 296*XNK(1,1) + 1554*XNK(2,1) + 
     >   1554*XNK(3,1) + 296*XNK(4,1) + 2052*XNK(4,2) - 2052*
     >   XNK(4,3) - 296*XNK(4,4) - 1554*XNK(3,4) - 1554*XNK(2,4) - 
     >   10773*XNK(2,3) + 10773*XNK(2,2)))/250880000 
      Q2(46,65) = Q2(46,65) + 
     >   (105**(0.5D0)*ETAH*(3633*XNJ(1,4,K2) + 15876*XNJ(2,4,K2) - 
     >   15876*XNJ(3,4,K2) - 3633*XNJ(4,4,K2) - 21798*XNJ(1,3,K2) - 
     >   95256*XNJ(2,3,K2) + 95256*XNJ(3,3,K2) + 21798*XNJ(4,3,K2) + 
     >   3633*XNJ(1,1,K2) + 15876*XNJ(2,1,K2) - 15876*XNJ(3,1,K2) - 
     >   3633*XNJ(4,1,K2) - 21798*XNJ(1,2,K2) - 95256*XNJ(2,2,K2) + 
     >   95256*XNJ(3,2,K2) + 21798*XNJ(4,2,K2)))/1254400000 - 
     >   (105**(0.5D0)*MUH*(14652*XNI(4,3,J2,K2) + 101574*
     >   XNI(3,3,J2,K2) - 101574*XNI(2,3,J2,K2) - 14652*XNI(1,3,J2,K2)-
     >   2442*XNI(4,4,J2,K2) - 16929*XNI(3,4,J2,K2) + 16929*
     >   XNI(2,4,J2,K2) + 
     >   2442*XNI(1,4,J2,K2) - 2442*XNI(4,1,J2,K2) - 16929*
     >   XNI(3,1,J2,K2) + 
     >   16929*XNI(2,1,J2,K2) + 2442*XNI(1,1,J2,K2) + 14652*
     >   XNI(4,2,J2,K2) + 
     >   101574*XNI(3,2,J2,K2) - 101574*XNI(2,2,J2,K2) - 14652*
     >   XNI(1,2,J2,K2)))/1254400000 + (105**(0.5D0)*XIH*
     >   SIGN(1.0,XIH)*(387828*XNK(3,2) - 387828*XNK(3,3) + 
     >   12802*XNK(1,4) + 88749*XNK(1,3) - 88749*XNK(1,2) - 
     >   12802*XNK(1,1) - 55944*XNK(2,1) + 55944*XNK(3,1) + 
     >   12802*XNK(4,1) + 88749*XNK(4,2) - 88749*XNK(4,3) - 
     >   12802*XNK(4,4) - 55944*XNK(3,4) + 55944*XNK(2,4) + 
     >   387828*XNK(2,3) - 387828*XNK(2,2)))/1254400000 
      Q2(47,65) = Q2(47,65) + 
     >   (27*7**(0.5D0)*MUH*SIGN(1.0,MUH)*(444*XNI(4,3,J2,K2) + 
     >   3078*XNI(3,3,J2,K2) - 3078*XNI(2,3,J2,K2) - 444*
     >   XNI(1,3,J2,K2) - 
     >   74*XNI(4,4,J2,K2) - 513*XNI(3,4,J2,K2) + 513*XNI(2,4,J2,K2) + 
     >   74*XNI(1,4,J2,K2) - 74*XNI(4,1,J2,K2) - 513*XNI(3,1,J2,K2) + 
     >   513*XNI(2,1,J2,K2) + 74*XNI(1,1,J2,K2) + 444*XNI(4,2,J2,K2) + 
     >   3078*XNI(3,2,J2,K2) - 3078*XNI(2,2,J2,K2) - 444*
     >   XNI(1,2,J2,K2)))/
     >   250880000 - (27*7**(0.5D0)*ETAH*(21*XNJ(1,4,K2) - 
     >   126*XNJ(2,4,K2) - 126*XNJ(3,4,K2) + 21*XNJ(4,4,K2) - 
     >   126*XNJ(1,3,K2) + 756*XNJ(2,3,K2) + 756*XNJ(3,3,K2) - 
     >   126*XNJ(4,3,K2) + 21*XNJ(1,1,K2) - 126*XNJ(2,1,K2) - 
     >   126*XNJ(3,1,K2) + 21*XNJ(4,1,K2) - 126*XNJ(1,2,K2) + 
     >   756*XNJ(2,2,K2) + 756*XNJ(3,2,K2) - 126*XNJ(4,2,K2)))/
     >   250880000 - (27*7**(0.5D0)*XIH*SIGN(1.0,XIH)*
     >   (3078*XNK(3,2) - 3078*XNK(3,3) + 74*XNK(1,4) + 513*
     >   XNK(1,3) - 513*XNK(1,2) - 74*XNK(1,1) + 444*XNK(2,1) + 
     >   444*XNK(3,1) - 74*XNK(4,1) - 513*XNK(4,2) + 513*XNK(4,3) + 
     >   74*XNK(4,4) - 444*XNK(3,4) - 444*XNK(2,4) - 3078*XNK(2,3) + 
     >   3078*XNK(2,2)))/250880000 
      Q2(48,65) = Q2(48,65) 
     >   - (27*5**(0.5D0)*ETAH*(1554*XNJ(1,4,K2) + 10773*
     >   XNJ(2,4,K2) - 10773*XNJ(3,4,K2) - 1554*XNJ(4,4,K2) - 9324*
     >   XNJ(1,3,K2) - 64638*XNJ(2,3,K2) + 64638*XNJ(3,3,K2) + 9324*
     >   XNJ(4,3,K2) + 1554*XNJ(1,1,K2) + 10773*XNJ(2,1,K2) - 10773*
     >   XNJ(3,1,K2) - 1554*XNJ(4,1,K2) - 9324*XNJ(1,2,K2) - 64638*
     >   XNJ(2,2,K2) + 64638*XNJ(3,2,K2) + 9324*XNJ(4,2,K2)))/
     >   CONST10 - (27*5**(0.5D0)*MUH*(9324*XNI(4,3,J2,K2) + 
     >   64638*XNI(3,3,J2,K2) - 64638*XNI(2,3,J2,K2) - 9324*
     >   XNI(1,3,J2,K2) - 1554*XNI(4,4,J2,K2) - 10773*XNI(3,4,J2,K2) + 
     >   10773*XNI(2,4,J2,K2) + 1554*XNI(1,4,J2,K2) - 1554*
     >   XNI(4,1,J2,K2) - 10773*XNI(3,1,J2,K2) + 10773*XNI(2,1,J2,K2) + 
     >   1554*XNI(1,1,J2,K2) + 9324*XNI(4,2,J2,K2) + 64638*
     >   XNI(3,2,J2,K2) - 64638*XNI(2,2,J2,K2) - 9324*XNI(1,2,J2,K2)))/
     >   CONST10 - (27*5**(0.5D0)*XIH*SIGN(1.0,XIH)*
     >   (263169*XNK(3,2) - 263169*XNK(3,3) + 5476*XNK(1,4) + 
     >   37962*XNK(1,3) - 37962*XNK(1,2) - 5476*XNK(1,1) - 
     >   37962*XNK(2,1) + 37962*XNK(3,1) + 5476*XNK(4,1) + 
     >   37962*XNK(4,2) - 37962*XNK(4,3) - 5476*XNK(4,4) - 
     >   37962*XNK(3,4) + 37962*XNK(2,4) + 263169*XNK(2,3) - 
     >   263169*XNK(2,2)))/CONST10 
      Q2(49,65) = Q2(49,65) 
     >   - (3*7**(0.5D0)*ETAH*SIGN(1.0,ETAH)*(296*XNJ(1,4,K2) + 
     >   1554*XNJ(2,4,K2) + 1554*XNJ(3,4,K2) + 296*XNJ(4,4,K2) + 
     >   2052*XNJ(1,3,K2) + 10773*XNJ(2,3,K2) + 10773*XNJ(3,3,K2) + 
     >   2052*XNJ(4,3,K2) - 296*XNJ(1,1,K2) - 1554*XNJ(2,1,K2) - 
     >   1554*XNJ(3,1,K2) - 296*XNJ(4,1,K2) - 2052*XNJ(1,2,K2) - 
     >   10773*XNJ(2,2,K2) - 10773*XNJ(3,2,K2) - 2052*XNJ(4,2,K2)))/
     >   50176000 - (3*7**(0.5D0)*MUH*SIGN(1.0,MUH)*
     >   (2052*XNI(4,3,J2,K2) + 10773*XNI(3,3,J2,K2) + 10773*
     >   XNI(2,3,J2,K2) + 2052*XNI(1,3,J2,K2) + 296*XNI(4,4,J2,K2) + 
     >   1554*XNI(3,4,J2,K2) + 1554*XNI(2,4,J2,K2) + 296*XNI(1,4,J2,K2)-
     >   296*XNI(4,1,J2,K2) - 1554*XNI(3,1,J2,K2) - 1554*XNI(2,1,J2,K2)-
     >   296*XNI(1,1,J2,K2) - 2052*XNI(4,2,J2,K2) - 10773*
     >   XNI(3,2,J2,K2) -
     >   10773*XNI(2,2,J2,K2) - 2052*XNI(1,2,J2,K2)))/50176000 - 
     >   (3*7**(0.5D0)*XIH*(9261*XNK(3,2) + 9261*XNK(3,3) + 
     >   336*XNK(1,4) + 1764*XNK(1,3) + 1764*XNK(1,2) + 336*
     >   XNK(1,1) + 1764*XNK(2,1) + 1764*XNK(3,1) + 336*XNK(4,1) + 
     >   1764*XNK(4,2) + 1764*XNK(4,3) + 336*XNK(4,4) + 1764*
     >   XNK(3,4) + 1764*XNK(2,4) + 9261*XNK(2,3) + 9261*
     >   XNK(2,2)))/50176000 
      Q2(50,65) = Q2(50,65) + 
     >   (21**(0.5D0)*ETAH*SIGN(1.0,ETAH)*(12802*XNJ(1,4,K2) + 
     >   55944*XNJ(2,4,K2) - 55944*XNJ(3,4,K2) - 12802*XNJ(4,4,K2) + 
     >   88749*XNJ(1,3,K2) + 387828*XNJ(2,3,K2) - 387828*XNJ(3,3,K2) - 
     >   88749*XNJ(4,3,K2) - 12802*XNJ(1,1,K2) - 55944*XNJ(2,1,K2) + 
     >   55944*XNJ(3,1,K2) + 12802*XNJ(4,1,K2) - 88749*XNJ(1,2,K2) - 
     >   387828*XNJ(2,2,K2) + 387828*XNJ(3,2,K2) + 88749*XNJ(4,2,K2)))/
     >   752640000 + (21**(0.5D0)*MUH*(67716*XNI(4,3,J2,K2) + 
     >   355509*XNI(3,3,J2,K2) + 355509*XNI(2,3,J2,K2) + 67716*
     >   XNI(1,3,J2,K2) + 9768*XNI(4,4,J2,K2) + 51282*XNI(3,4,J2,K2) + 
     >   51282*XNI(2,4,J2,K2) + 9768*XNI(1,4,J2,K2) - 9768*
     >   XNI(4,1,J2,K2) - 
     >   51282*XNI(3,1,J2,K2) - 51282*XNI(2,1,J2,K2) - 9768*
     >   XNI(1,1,J2,K2) - 67716*XNI(4,2,J2,K2) - 355509*XNI(3,2,J2,K2) -
     >   355509*XNI(2,2,J2,K2) - 67716*XNI(1,2,J2,K2)))/752640000 - 
     >   (21**(0.5D0)*XIH*(333396*XNK(3,2) + 333396*
     >   XNK(3,3) - 14532*XNK(1,4) - 76293*XNK(1,3) - 76293*
     >   XNK(1,2) - 14532*XNK(1,1) - 63504*XNK(2,1) + 63504*
     >   XNK(3,1) + 14532*XNK(4,1) + 76293*XNK(4,2) + 76293*
     >   XNK(4,3) + 14532*XNK(4,4) + 63504*XNK(3,4) - 63504*
     >   XNK(2,4) - 333396*XNK(2,3) - 333396*XNK(2,2)))/752640000 
      Q2(51,65) = Q2(51,65) + 
     >   (9*35**(0.5D0)*XIH*(2646*XNK(3,2) + 2646*XNK(3,3) - 
     >   84*XNK(1,4) - 441*XNK(1,3) - 441*XNK(1,2) - 84*XNK(1,1) + 
     >   504*XNK(2,1) + 504*XNK(3,1) - 84*XNK(4,1) - 441*XNK(4,2) - 
     >   441*XNK(4,3) - 84*XNK(4,4) + 504*XNK(3,4) + 504*XNK(2,4) + 
     >   2646*XNK(2,3) + 2646*XNK(2,2)))/250880000 - (9*35**(0.5D0)*
     >   MUH*SIGN(1.0,MUH)*(2052*XNI(4,3,J2,K2) + 10773*
     >   XNI(3,3,J2,K2) + 10773*XNI(2,3,J2,K2) + 2052*XNI(1,3,J2,K2) + 
     >   296*XNI(4,4,J2,K2) + 1554*XNI(3,4,J2,K2) + 1554*XNI(2,4,J2,K2)+
     >   296*XNI(1,4,J2,K2) - 296*XNI(4,1,J2,K2) - 1554*XNI(3,1,J2,K2)-
     >   1554*XNI(2,1,J2,K2) - 296*XNI(1,1,J2,K2) - 2052*XNI(4,2,J2,K2)-
     >   10773*XNI(3,2,J2,K2) - 10773*XNI(2,2,J2,K2) - 2052*
     >   XNI(1,2,J2,K2)))/250880000 - (9*35**(0.5D0)*ETAH*
     >   SIGN(1.0,ETAH)*(74*XNJ(1,4,K2) - 444*XNJ(2,4,K2) - 444*
     >   XNJ(3,4,K2) + 74*XNJ(4,4,K2) + 513*XNJ(1,3,K2) - 3078*
     >   XNJ(2,3,K2) - 3078*XNJ(3,3,K2) + 513*XNJ(4,3,K2) - 74*
     >   XNJ(1,1,K2) + 444*XNJ(2,1,K2) + 444*XNJ(3,1,K2) - 74*
     >   XNJ(4,1,K2) - 513*XNJ(1,2,K2) + 3078*XNJ(2,2,K2) + 3078*
     >   XNJ(3,2,K2) - 513*XNJ(4,2,K2)))/250880000 
      Q2(52,65) = Q2(52,65) + 
     >   ((13851*XNI(4,3,J2,K2))/62720000 + (41553*XNI(3,3,J2,K2))/
     >   35840000 + (41553*XNI(2,3,J2,K2))/35840000 + (13851*
     >   XNI(1,3,J2,K2))/62720000 + (999*XNI(4,4,J2,K2))/31360000 + 
     >   (2997*XNI(3,4,J2,K2))/17920000 + (2997*XNI(2,4,J2,K2))/
     >   17920000 + (999*XNI(1,4,J2,K2))/31360000 - (999*
     >   XNI(4,1,J2,K2))/31360000 - (2997*XNI(3,1,J2,K2))/17920000 - 
     >   (2997*XNI(2,1,J2,K2))/17920000 - (999*XNI(1,1,J2,K2))/31360000-
     >   (13851*XNI(4,2,J2,K2))/62720000 - (41553*XNI(3,2,J2,K2))/
     >   35840000 - (41553*XNI(2,2,J2,K2))/35840000 - (13851*
     >   XNI(1,2,J2,K2))/62720000)*MUH + ((170829*XNJ(3,4,K2))/
     >   878080000 - (170829*XNJ(2,4,K2))/878080000 - (12321*
     >   XNJ(1,4,K2))/439040000 + (12321*XNJ(4,4,K2))/439040000 - 
     >   (170829*XNJ(1,3,K2))/878080000 - (2368521*XNJ(2,3,K2))/
     >   1756160000 + (2368521*XNJ(3,3,K2))/1756160000 + (170829*
     >   XNJ(4,3,K2))/878080000 + (12321*XNJ(1,1,K2))/439040000 + 
     >   (170829*XNJ(2,1,K2))/878080000 - (170829*XNJ(3,1,K2))/
     >   878080000 - (12321*XNJ(4,1,K2))/439040000 + (170829*
     >   XNJ(1,2,K2))/878080000 + (2368521*XNJ(2,2,K2))/1756160000 - 
     >   (2368521*XNJ(3,2,K2))/1756160000 - (170829*XNJ(4,2,K2))/
     >   878080000)*ETAH*SIGN(1.0,ETAH) + ((41553*XNK(3,2))/
     >   35840000 + (41553*XNK(3,3))/35840000 - (999*XNK(1,4))/
     >   31360000 - (2997*XNK(1,3))/17920000 - (2997*XNK(1,2))/
     >   17920000 - (999*XNK(1,1))/31360000 - (13851*XNK(2,1))/
     >   62720000 + (13851*XNK(3,1))/62720000 + (999*XNK(4,1))/
     >   31360000 + (2997*XNK(4,2))/17920000 + (2997*XNK(4,3))/
     >   17920000 + (999*XNK(4,4))/31360000 + (13851*XNK(3,4))/
     >   62720000 - (13851*XNK(2,4))/62720000 - (41553*XNK(2,3))/
     >   35840000 - (41553*XNK(2,2))/35840000)*XIH 
      Q2(53,65) = Q2(53,65) + 
     >   (21**(0.5D0)*ETAH*(9768*XNJ(1,4,K2) + 51282*XNJ(2,4,K2) + 
     >   51282*XNJ(3,4,K2) + 9768*XNJ(4,4,K2) + 67716*XNJ(1,3,K2) + 
     >   355509*XNJ(2,3,K2) + 355509*XNJ(3,3,K2) + 67716*XNJ(4,3,K2)-
     >   9768*XNJ(1,1,K2) - 51282*XNJ(2,1,K2) - 51282*XNJ(3,1,K2) - 
     >   9768*XNJ(4,1,K2) - 67716*XNJ(1,2,K2) - 355509*XNJ(2,2,K2) - 
     >   355509*XNJ(3,2,K2) - 67716*XNJ(4,2,K2)))/752640000 - 
     >   (21**(0.5D0)*MUH*SIGN(1.0,MUH)*(88749*XNI(4,3,J2,K2) + 
     >   387828*XNI(3,3,J2,K2) - 387828*XNI(2,3,J2,K2) - 88749*
     >   XNI(1,3,J2,K2) + 12802*XNI(4,4,J2,K2) + 55944*XNI(3,4,J2,K2)-
     >   55944*XNI(2,4,J2,K2) - 12802*XNI(1,4,J2,K2) - 12802*
     >   XNI(4,1,J2,K2) - 55944*XNI(3,1,J2,K2) + 55944*XNI(2,1,J2,K2)+
     >   12802*XNI(1,1,J2,K2) - 88749*XNI(4,2,J2,K2) - 387828*
     >   XNI(3,2,J2,K2) + 387828*XNI(2,2,J2,K2) + 88749*
     >   XNI(1,2,J2,K2)))/
     >   752640000 + (21**(0.5D0)*XIH*(333396*XNK(3,2) - 
     >   333396*XNK(3,3) - 14532*XNK(1,4) - 63504*XNK(1,3) + 
     >   63504*XNK(1,2) + 14532*XNK(1,1) + 76293*XNK(2,1) + 
     >   76293*XNK(3,1) + 14532*XNK(4,1) + 63504*XNK(4,2) - 
     >   63504*XNK(4,3) - 14532*XNK(4,4) - 76293*XNK(3,4) - 
     >   76293*XNK(2,4) - 333396*XNK(2,3) + 333396*XNK(2,2)))/
     >   752640000 
      Q2(54,65) = Q2(54,65) + 
     >   (7**(0.5D0)*MUH*(976239*XNI(4,3,J2,K2) + 4266108*
     >   XNI(3,3,J2,K2) - 4266108*XNI(2,3,J2,K2) - 976239*
     >   XNI(1,3,J2,K2) + 
     >   140822*XNI(4,4,J2,K2) + 615384*XNI(3,4,J2,K2) - 615384*
     >   XNI(2,4,J2,K2) - 140822*XNI(1,4,J2,K2) - 140822*
     >   XNI(4,1,J2,K2) - 
     >   615384*XNI(3,1,J2,K2) + 615384*XNI(2,1,J2,K2) + 140822*
     >   XNI(1,1,J2,K2) - 976239*XNI(4,2,J2,K2) - 4266108*
     >   XNI(3,2,J2,K2) + 4266108*XNI(2,2,J2,K2) + 976239*
     >   XNI(1,2,J2,K2)))/CONST9 - (7**(0.5D0)*ETAH*
     >   (140822*XNJ(1,4,K2) + 615384*XNJ(2,4,K2) - 615384*
     >   XNJ(3,4,K2) - 140822*XNJ(4,4,K2) + 976239*XNJ(1,3,K2) + 
     >   4266108*XNJ(2,3,K2) - 4266108*XNJ(3,3,K2) - 976239*
     >   XNJ(4,3,K2) - 140822*XNJ(1,1,K2) - 615384*XNJ(2,1,K2) + 
     >   615384*XNJ(3,1,K2) + 140822*XNJ(4,1,K2) - 976239*
     >   XNJ(1,2,K2) - 
     >   4266108*XNJ(2,2,K2) + 4266108*XNJ(3,2,K2) + 976239*
     >   XNJ(4,2,K2)))/CONST9 + (7**(0.5D0)*XIH*
     >   (4000752*XNK(3,2) - 4000752*XNK(3,3) + 209503*XNK(1,4) + 
     >   915516*XNK(1,3) - 915516*XNK(1,2) - 209503*XNK(1,1) - 
     >   915516*XNK(2,1) + 915516*XNK(3,1) + 209503*XNK(4,1) + 
     >   915516*XNK(4,2) - 915516*XNK(4,3) - 209503*XNK(4,4) - 
     >   915516*XNK(3,4) + 915516*XNK(2,4) + 4000752*XNK(2,3) - 
     >   4000752*XNK(2,2)))/CONST9 
      Q2(55,65) = Q2(55,65) + 
     >   (105**(0.5D0)*ETAH*(2442*XNJ(1,4,K2) - 14652*
     >   XNJ(2,4,K2) - 14652*XNJ(3,4,K2) + 2442*XNJ(4,4,K2) + 
     >   16929*XNJ(1,3,K2) - 101574*XNJ(2,3,K2) - 101574*XNJ(3,3,K2) + 
     >   16929*XNJ(4,3,K2) - 2442*XNJ(1,1,K2) + 14652*XNJ(2,1,K2) + 
     >   14652*XNJ(3,1,K2) - 2442*XNJ(4,1,K2) - 16929*XNJ(1,2,K2) + 
     >   101574*XNJ(2,2,K2) + 101574*XNJ(3,2,K2) - 16929*XNJ(4,2,K2)))/
     >   1254400000 - (105**(0.5D0)*MUH*SIGN(1.0,MUH)*
     >   (88749*XNI(4,3,J2,K2) + 387828*XNI(3,3,J2,K2) - 387828*
     >   XNI(2,3,J2,K2) - 88749*XNI(1,3,J2,K2) + 12802*XNI(4,4,J2,K2) + 
     >   55944*XNI(3,4,J2,K2) - 55944*XNI(2,4,J2,K2) - 12802*
     >   XNI(1,4,J2,K2) - 12802*XNI(4,1,J2,K2) - 55944*XNI(3,1,J2,K2) + 
     >   55944*XNI(2,1,J2,K2) + 12802*XNI(1,1,J2,K2) - 88749*
     >   XNI(4,2,J2,K2) - 387828*XNI(3,2,J2,K2) + 387828*
     >   XNI(2,2,J2,K2) + 
     >   88749*XNI(1,2,J2,K2)))/1254400000 - (105**(0.5D0)*XIH*
     >   (95256*XNK(3,2) - 95256*XNK(3,3) + 3633*XNK(1,4) + 15876*
     >   XNK(1,3) - 15876*XNK(1,2) - 3633*XNK(1,1) + 21798*XNK(2,1) + 
     >   21798*XNK(3,1) - 3633*XNK(4,1) - 15876*XNK(4,2) + 15876*
     >   XNK(4,3) + 3633*XNK(4,4) - 21798*XNK(3,4) - 21798*
     >   XNK(2,4) - 95256*XNK(2,3) + 95256*XNK(2,2)))/1254400000 
      Q2(56,65) = Q2(56,65) + 
     >   (3*3**(0.5D0)*ETAH*(60236*XNJ(1,4,K2) + 417582*
     >   XNJ(2,4,K2) - 417582*XNJ(3,4,K2) - 60236*XNJ(4,4,K2) + 
     >   417582*XNJ(1,3,K2) + 2894859*XNJ(2,3,K2) - 2894859*
     >   XNJ(3,3,K2) - 417582*XNJ(4,3,K2) - 60236*XNJ(1,1,K2) - 
     >   417582*XNJ(2,1,K2) + 417582*XNJ(3,1,K2) + 60236*
     >   XNJ(4,1,K2) - 417582*XNJ(1,2,K2) - 2894859*XNJ(2,2,K2) + 
     >   2894859*XNJ(3,2,K2) + 417582*XNJ(4,2,K2)))/CONST10 + 
     >   (3*3**(0.5D0)*MUH*(621243*XNI(4,3,J2,K2) + 2714796*
     >   XNI(3,3,J2,K2) - 2714796*XNI(2,3,J2,K2) - 621243*
     >   XNI(1,3,J2,K2) + 
     >   89614*XNI(4,4,J2,K2) + 391608*XNI(3,4,J2,K2) - 391608*
     >   XNI(2,4,J2,K2) - 89614*XNI(1,4,J2,K2) - 89614*XNI(4,1,J2,K2) -
     >   391608*XNI(3,1,J2,K2) + 391608*XNI(2,1,J2,K2) + 89614*
     >   XNI(1,1,J2,K2) - 621243*XNI(4,2,J2,K2) - 2714796*
     >   XNI(3,2,J2,K2) + 
     >   2714796*XNI(2,2,J2,K2) + 621243*XNI(1,2,J2,K2)))/CONST10 - 
     >   (3*3**(0.5D0)*XIH*(2714796*XNK(3,2) - 2714796*
     >   XNK(3,3) + 89614*XNK(1,4) + 391608*XNK(1,3) - 391608*
     >   XNK(1,2) - 89614*XNK(1,1) - 621243*XNK(2,1) + 621243*
     >   XNK(3,1) + 89614*XNK(4,1) + 391608*XNK(4,2) - 391608*
     >   XNK(4,3) - 89614*XNK(4,4) - 621243*XNK(3,4) + 621243*
     >   XNK(2,4) + 2714796*XNK(2,3) - 2714796*XNK(2,2)))/
     >   CONST10 
      Q2(57,65) = Q2(57,65) + 
     >   (9*35**(0.5D0)*XIH*(2646*XNK(3,2) + 2646*XNK(3,3) - 
     >   84*XNK(1,4) + 504*XNK(1,3) + 504*XNK(1,2) - 84*XNK(1,1) - 
     >   441*XNK(2,1) - 441*XNK(3,1) - 84*XNK(4,1) + 504*XNK(4,2) + 
     >   504*XNK(4,3) - 84*XNK(4,4) - 441*XNK(3,4) - 441*XNK(2,4) + 
     >   2646*XNK(2,3) + 2646*XNK(2,2)))/250880000 - (9*35**(0.5D0)*
     >   MUH*SIGN(1.0,MUH)*(513*XNI(4,3,J2,K2) - 3078*
     >   XNI(3,3,J2,K2) - 3078*XNI(2,3,J2,K2) + 513*XNI(1,3,J2,K2) + 
     >   74*XNI(4,4,J2,K2) - 444*XNI(3,4,J2,K2) - 444*XNI(2,4,J2,K2) + 
     >   74*XNI(1,4,J2,K2) - 74*XNI(4,1,J2,K2) + 444*XNI(3,1,J2,K2) + 
     >   444*XNI(2,1,J2,K2) - 74*XNI(1,1,J2,K2) - 513*XNI(4,2,J2,K2) + 
     >   3078*XNI(3,2,J2,K2) + 3078*XNI(2,2,J2,K2) - 513*
     >   XNI(1,2,J2,K2)))/
     >   250880000 - (9*35**(0.5D0)*ETAH*SIGN(1.0,ETAH)*(296*
     >   XNJ(1,4,K2) + 1554*XNJ(2,4,K2) + 1554*XNJ(3,4,K2) + 296*
     >   XNJ(4,4,K2) + 2052*XNJ(1,3,K2) + 10773*XNJ(2,3,K2) + 10773*
     >   XNJ(3,3,K2) + 2052*XNJ(4,3,K2) - 296*XNJ(1,1,K2) - 1554*
     >   XNJ(2,1,K2) - 1554*XNJ(3,1,K2) - 296*XNJ(4,1,K2) - 2052*
     >   XNJ(1,2,K2) - 10773*XNJ(2,2,K2) - 10773*XNJ(3,2,K2) - 2052*
     >   XNJ(4,2,K2)))/250880000 
      Q2(58,65) = Q2(58,65) + 
     >   (105**(0.5D0)*ETAH*SIGN(1.0,ETAH)*(12802*XNJ(1,4,K2) + 
     >   55944*XNJ(2,4,K2) - 55944*XNJ(3,4,K2) - 12802*XNJ(4,4,K2) + 
     >   88749*XNJ(1,3,K2) + 387828*XNJ(2,3,K2) - 387828*XNJ(3,3,K2) - 
     >   88749*XNJ(4,3,K2) - 12802*XNJ(1,1,K2) - 55944*XNJ(2,1,K2) + 
     >   55944*XNJ(3,1,K2) + 12802*XNJ(4,1,K2) - 88749*XNJ(1,2,K2) - 
     >   387828*XNJ(2,2,K2) + 387828*XNJ(3,2,K2) + 88749*XNJ(4,2,K2)))/
     >   1254400000 + (105**(0.5D0)*MUH*(16929*XNI(4,3,J2,K2) - 
     >   101574*XNI(3,3,J2,K2) - 101574*XNI(2,3,J2,K2) + 16929*
     >   XNI(1,3,J2,K2) + 2442*XNI(4,4,J2,K2) - 14652*XNI(3,4,J2,K2) - 
     >   14652*XNI(2,4,J2,K2) + 2442*XNI(1,4,J2,K2) - 2442*
     >   XNI(4,1,J2,K2) + 
     >   14652*XNI(3,1,J2,K2) + 14652*XNI(2,1,J2,K2) - 2442*
     >   XNI(1,1,J2,K2) - 16929*XNI(4,2,J2,K2) + 101574*XNI(3,2,J2,K2)+
     >   101574*XNI(2,2,J2,K2) - 16929*XNI(1,2,J2,K2)))/1254400000 + 
     >   (105**(0.5D0)*XIH*(95256*XNK(3,2) + 95256*XNK(3,3) + 
     >   3633*XNK(1,4) - 21798*XNK(1,3) - 21798*XNK(1,2) + 
     >   3633*XNK(1,1) + 15876*XNK(2,1) - 15876*XNK(3,1) - 
     >   3633*XNK(4,1) + 21798*XNK(4,2) + 21798*XNK(4,3) - 
     >   3633*XNK(4,4) - 15876*XNK(3,4) + 15876*XNK(2,4) - 
     >   95256*XNK(2,3) - 95256*XNK(2,2)))/1254400000 
      Q2(59,65) = Q2(59,65) 
     >   - (27*7**(0.5D0)*ETAH*SIGN(1.0,ETAH)*(74*
     >   XNJ(1,4,K2) - 444*XNJ(2,4,K2) - 444*XNJ(3,4,K2) + 74*
     >   XNJ(4,4,K2) + 513*XNJ(1,3,K2) - 3078*XNJ(2,3,K2) - 3078*
     >   XNJ(3,3,K2) + 513*XNJ(4,3,K2) - 74*XNJ(1,1,K2) + 444*
     >   XNJ(2,1,K2) + 444*XNJ(3,1,K2) - 74*XNJ(4,1,K2) - 513*
     >   XNJ(1,2,K2) + 3078*XNJ(2,2,K2) + 3078*XNJ(3,2,K2) - 513*
     >   XNJ(4,2,K2)))/250880000 - (27*7**(0.5D0)*MUH*
     >   SIGN(1.0,MUH)*(513*XNI(4,3,J2,K2) - 3078*XNI(3,3,J2,K2) - 
     >   3078*XNI(2,3,J2,K2) + 513*XNI(1,3,J2,K2) + 74*XNI(4,4,J2,K2)-
     >   444*XNI(3,4,J2,K2) - 444*XNI(2,4,J2,K2) + 74*XNI(1,4,J2,K2) - 
     >   74*XNI(4,1,J2,K2) + 444*XNI(3,1,J2,K2) + 444*XNI(2,1,J2,K2) - 
     >   74*XNI(1,1,J2,K2) - 513*XNI(4,2,J2,K2) + 3078*XNI(3,2,J2,K2)+
     >   3078*XNI(2,2,J2,K2) - 513*XNI(1,2,J2,K2)))/250880000 - 
     >   (27*7**(0.5D0)*XIH*(756*XNK(3,2) + 756*XNK(3,3) + 
     >   21*XNK(1,4) - 126*XNK(1,3) - 126*XNK(1,2) + 21*XNK(1,1) - 
     >   126*XNK(2,1) - 126*XNK(3,1) + 21*XNK(4,1) - 126*XNK(4,2) - 
     >   126*XNK(4,3) + 21*XNK(4,4) - 126*XNK(3,4) - 126*XNK(2,4) + 
     >   756*XNK(2,3) + 756*XNK(2,2)))/250880000 
      Q2(60,65) = Q2(60,65) + 
     >   (27*5**(0.5D0)*MUH*(10773*XNI(4,3,J2,K2) - 64638*
     >   XNI(3,3,J2,K2) - 64638*XNI(2,3,J2,K2) + 10773*XNI(1,3,J2,K2) + 
     >   1554*XNI(4,4,J2,K2) - 9324*XNI(3,4,J2,K2) - 
     >   9324*XNI(2,4,J2,K2) + 
     >   1554*XNI(1,4,J2,K2) - 1554*XNI(4,1,J2,K2) + 
     >   9324*XNI(3,1,J2,K2) + 
     >   9324*XNI(2,1,J2,K2) - 1554*XNI(1,1,J2,K2) - 
     >   10773*XNI(4,2,J2,K2) + 
     >   64638*XNI(3,2,J2,K2) + 64638*XNI(2,2,J2,K2) - 10773*
     >   XNI(1,2,J2,K2)))/CONST10 - (27*5**(0.5D0)*ETAH*
     >   SIGN(1.0,ETAH)*(5476*XNJ(1,4,K2) + 37962*XNJ(2,4,K2) - 
     >   37962*XNJ(3,4,K2) - 5476*XNJ(4,4,K2) + 37962*XNJ(1,3,K2) + 
     >   263169*XNJ(2,3,K2) - 263169*XNJ(3,3,K2) - 37962*XNJ(4,3,K2) -
     >   5476*XNJ(1,1,K2) - 37962*XNJ(2,1,K2) + 37962*XNJ(3,1,K2) + 
     >   5476*XNJ(4,1,K2) - 37962*XNJ(1,2,K2) - 263169*XNJ(2,2,K2) + 
     >   263169*XNJ(3,2,K2) + 37962*XNJ(4,2,K2)))/CONST10 - 
     >   (27*5**(0.5D0)*XIH*(64638*XNK(3,2) + 64638*XNK(3,3) + 
     >   1554*XNK(1,4) - 9324*XNK(1,3) - 9324*XNK(1,2) + 1554*
     >   XNK(1,1) + 10773*XNK(2,1) - 10773*XNK(3,1) - 1554*
     >   XNK(4,1) + 9324*XNK(4,2) + 9324*XNK(4,3) - 1554*
     >   XNK(4,4) - 10773*XNK(3,4) + 10773*XNK(2,4) - 64638*
     >   XNK(2,3) - 64638*XNK(2,2)))/CONST10 
      Q2(61,65) = Q2(61,65) + 
     >   ((170829*XNI(4,3,J2,K2))/878080000 + (2368521*
     >   XNI(3,3,J2,K2))/1756160000 - (2368521*XNI(2,3,J2,K2))/
     >   1756160000 - (170829*XNI(1,3,J2,K2))/878080000 + (12321*
     >   XNI(4,4,J2,K2))/439040000 + (170829*XNI(3,4,J2,K2))/
     >   878080000 - (170829*XNI(2,4,J2,K2))/878080000 - 
     >   (12321*XNI(1,4,J2,K2))/439040000 - (12321*XNI(4,1,J2,K2))/
     >   439040000 - (170829*XNI(3,1,J2,K2))/878080000 + (170829*
     >   XNI(2,1,J2,K2))/878080000 + (12321*XNI(1,1,J2,K2))/439040000 -
     >   (170829*XNI(4,2,J2,K2))/878080000 - (2368521*XNI(3,2,J2,K2))/
     >   1756160000 + (2368521*XNI(2,2,J2,K2))/1756160000 + (170829*
     >   XNI(1,2,J2,K2))/878080000)*MUH*SIGN(1.0,MUH) + 
     >   ((999*XNJ(1,4,K2))/31360000 + (2997*XNJ(2,4,K2))/17920000 + 
     >   (2997*XNJ(3,4,K2))/17920000 + (999*XNJ(4,4,K2))/31360000 + 
     >   (13851*XNJ(1,3,K2))/62720000 + (41553*XNJ(2,3,K2))/35840000 + 
     >   (41553*XNJ(3,3,K2))/35840000 + (13851*XNJ(4,3,K2))/62720000 - 
     >   (999*XNJ(1,1,K2))/31360000 - (2997*XNJ(2,1,K2))/17920000 - 
     >   (2997*XNJ(3,1,K2))/17920000 - (999*XNJ(4,1,K2))/31360000 - 
     >   (13851*XNJ(1,2,K2))/62720000 - (41553*XNJ(2,2,K2))/35840000 -
     >   (41553*XNJ(3,2,K2))/35840000 - (13851*XNJ(4,2,K2))/62720000)*
     >   ETAH + ((41553*XNK(3,3))/35840000 - (41553*XNK(3,2))/
     >   35840000 + (999*XNK(1,4))/31360000 + (13851*XNK(1,3))/
     >   62720000 - (13851*XNK(1,2))/62720000 - (999*XNK(1,1))/
     >   31360000 - (2997*XNK(2,1))/17920000 - (2997*XNK(3,1))/
     >   17920000 - (999*XNK(4,1))/31360000 - (13851*XNK(4,2))/
     >   62720000 + (13851*XNK(4,3))/62720000 + (999*XNK(4,4))/
     >   31360000 + (2997*XNK(3,4))/17920000 + (2997*XNK(2,4))/
     >   17920000 + (41553*XNK(2,3))/35840000 - (41553*XNK(2,2))/
     >   35840000)*XIH 
      Q2(62,65) = Q2(62,65) 
     >   - (3*3**(0.5D0)*ETAH*(89614*XNJ(1,4,K2) + 391608*
     >   XNJ(2,4,K2) - 391608*XNJ(3,4,K2) - 89614*XNJ(4,4,K2) + 
     >   621243*XNJ(1,3,K2) + 2714796*XNJ(2,3,K2) - 2714796*
     >   XNJ(3,3,K2) - 621243*XNJ(4,3,K2) - 89614*XNJ(1,1,K2) - 
     >   391608*XNJ(2,1,K2) + 391608*XNJ(3,1,K2) + 89614*XNJ(4,1,K2) - 
     >   621243*XNJ(1,2,K2) - 2714796*XNJ(2,2,K2) + 2714796*
     >   XNJ(3,2,K2) + 621243*XNJ(4,2,K2)))/CONST10 - 
     >   (3*3**(0.5D0)*MUH*(417582*XNI(4,3,J2,K2) + 2894859*
     >   XNI(3,3,J2,K2) - 2894859*XNI(2,3,J2,K2) - 417582*
     >   XNI(1,3,J2,K2) + 
     >   60236*XNI(4,4,J2,K2) + 417582*XNI(3,4,J2,K2) - 417582*
     >   XNI(2,4,J2,K2) - 60236*XNI(1,4,J2,K2) - 60236*XNI(4,1,J2,K2) -
     >   417582*XNI(3,1,J2,K2) + 417582*XNI(2,1,J2,K2) + 60236*
     >   XNI(1,1,J2,K2) - 417582*XNI(4,2,J2,K2) - 2894859*
     >   XNI(3,2,J2,K2) + 2894859*XNI(2,2,J2,K2) + 417582*
     >   XNI(1,2,J2,K2)))/CONST10 - (3*3**(0.5D0)*XIH*
     >   (2714796*XNK(3,2) - 2714796*XNK(3,3) + 89614*XNK(1,4) + 
     >   621243*XNK(1,3) - 621243*XNK(1,2) - 89614*XNK(1,1) - 
     >   391608*XNK(2,1) + 391608*XNK(3,1) + 89614*XNK(4,1) + 
     >   621243*XNK(4,2) - 621243*XNK(4,3) - 89614*XNK(4,4) - 
     >   391608*XNK(3,4) + 391608*XNK(2,4) + 2714796*XNK(2,3) - 
     >   2714796*XNK(2,2)))/CONST10 
      Q2(63,65) = Q2(63,65) + 
     >   (27*5**(0.5D0)*ETAH*(1554*XNJ(1,4,K2) - 9324*
     >   XNJ(2,4,K2) - 9324*XNJ(3,4,K2) + 1554*XNJ(4,4,K2) + 10773*
     >   XNJ(1,3,K2) - 64638*XNJ(2,3,K2) - 64638*XNJ(3,3,K2) + 10773*
     >   XNJ(4,3,K2) - 1554*XNJ(1,1,K2) + 9324*XNJ(2,1,K2) + 9324*
     >   XNJ(3,1,K2) - 1554*XNJ(4,1,K2) - 10773*XNJ(1,2,K2) + 64638*
     >   XNJ(2,2,K2) + 64638*XNJ(3,2,K2) - 10773*XNJ(4,2,K2)))/
     >   CONST10 + (27*5**(0.5D0)*MUH*SIGN(1.0,MUH)*
     >   (37962*XNI(4,3,J2,K2) + 263169*XNI(3,3,J2,K2) - 263169*
     >   XNI(2,3,J2,K2) - 37962*XNI(1,3,J2,K2) + 5476*XNI(4,4,J2,K2) + 
     >   37962*XNI(3,4,J2,K2) - 37962*XNI(2,4,J2,K2) - 5476*
     >   XNI(1,4,J2,K2) - 5476*XNI(4,1,J2,K2) - 37962*XNI(3,1,J2,K2) + 
     >   37962*XNI(2,1,J2,K2) + 5476*XNI(1,1,J2,K2) - 37962*
     >   XNI(4,2,J2,K2) - 263169*XNI(3,2,J2,K2) + 263169*
     >   XNI(2,2,J2,K2) + 
     >   37962*XNI(1,2,J2,K2)))/CONST10 + (27*5**(0.5D0)*
     >   XIH*(64638*XNK(3,2) - 64638*XNK(3,3) + 1554*
     >   XNK(1,4) + 10773*XNK(1,3) - 10773*XNK(1,2) - 1554*
     >   XNK(1,1) + 9324*XNK(2,1) + 9324*XNK(3,1) - 1554*
     >   XNK(4,1) - 10773*XNK(4,2) + 10773*XNK(4,3) + 1554*
     >   XNK(4,4) - 9324*XNK(3,4) - 9324*XNK(2,4) - 64638*
     >   XNK(2,3) + 64638*XNK(2,2)))/CONST10 
      Q2(64,65) = Q2(64,65) + 
     >   (81*7**(0.5D0)*ETAH*(5476*XNJ(1,4,K2) + 37962*
     >   XNJ(2,4,K2) - 37962*XNJ(3,4,K2) - 5476*XNJ(4,4,K2) + 37962*
     >   XNJ(1,3,K2) + 263169*XNJ(2,3,K2) - 263169*XNJ(3,3,K2) - 37962*
     >   XNJ(4,3,K2) - 5476*XNJ(1,1,K2) - 37962*XNJ(2,1,K2) + 37962*
     >   XNJ(3,1,K2) + 5476*XNJ(4,1,K2) - 37962*XNJ(1,2,K2) - 263169*
     >   XNJ(2,2,K2) + 263169*XNJ(3,2,K2) + 37962*XNJ(4,2,K2)))/
     >   CONST11 - (81*7**(0.5D0)*MUH*(37962*
     >   XNI(4,3,J2,K2) + 263169*XNI(3,3,J2,K2) - 263169*XNI(2,3,J2,K2)-
     >   37962*XNI(1,3,J2,K2) + 5476*XNI(4,4,J2,K2) + 37962*
     >   XNI(3,4,J2,K2) - 37962*XNI(2,4,J2,K2) - 5476*XNI(1,4,J2,K2) - 
     >   5476*XNI(4,1,J2,K2) - 37962*XNI(3,1,J2,K2) + 37962*
     >   XNI(2,1,J2,K2) + 5476*XNI(1,1,J2,K2) - 37962*XNI(4,2,J2,K2) - 
     >   263169*XNI(3,2,J2,K2) + 263169*XNI(2,2,J2,K2) + 37962*
     >   XNI(1,2,J2,K2)))/CONST11 + (81*7**(0.5D0)*XIH*
     >   (263169*XNK(3,2) - 263169*XNK(3,3) + 5476*XNK(1,4) + 
     >   37962*XNK(1,3) - 37962*XNK(1,2) - 5476*XNK(1,1) - 37962*
     >   XNK(2,1) + 37962*XNK(3,1) + 5476*XNK(4,1) + 37962*
     >   XNK(4,2) - 37962*XNK(4,3) - 5476*XNK(4,4) - 37962*
     >   XNK(3,4) + 37962*XNK(2,4) + 263169*XNK(2,3) - 263169*
     >   XNK(2,2)))/CONST11 
*** ---------------------------------------------------------------- ***
*** ---------------------------------------------------------------- ***
*** ---------------------------------------------------------------- ***
      ENDIF
*
      CALL ALSBD(IELEM**3,1,Q2,IER,IELEM**3)
      IF(IER.NE.0) CALL XABORT('SNFTH4: SINGULAR MATRIX.')
*
      IF(IELEM.EQ.2)THEN
*** ---------------------------------------------------------------- ***
*** ---------------------------------------------------------------- ***
*** ---------------------------------------------------------------- ***
      CORNERQ(01,1) = (Q2(01,IELEM**3+1) + 3*Q2(04,IELEM**3+1) + 
     >   3*Q2(06,IELEM**3+1) + 3*Q2(07,IELEM**3+1) - 
     >   3**(0.5D0)*Q2(02,IELEM**3+1) - 3**(0.5D0)*Q2(03,IELEM**3+1) - 
     >   3**(0.5D0)*Q2(05,IELEM**3+1) - 3*3**(0.5D0)*Q2(08,IELEM**3+1)) 
      CORNERQ(02,1) = (Q2(01,IELEM**3+1) - 3*Q2(04,IELEM**3+1) - 
     >   3*Q2(06,IELEM**3+1) + 3*Q2(07,IELEM**3+1) + 
     >   3**(0.5D0)*Q2(02,IELEM**3+1) - 3**(0.5D0)*Q2(03,IELEM**3+1) - 
     >   3**(0.5D0)*Q2(05,IELEM**3+1) + 3*3**(0.5D0)*Q2(08,IELEM**3+1)) 
      CORNERQ(03,1) = (Q2(01,IELEM**3+1) - 3*Q2(04,IELEM**3+1) + 
     >   3*Q2(06,IELEM**3+1) - 3*Q2(07,IELEM**3+1) - 
     >   3**(0.5D0)*Q2(02,IELEM**3+1) + 3**(0.5D0)*Q2(03,IELEM**3+1) - 
     >   3**(0.5D0)*Q2(05,IELEM**3+1) + 3*3**(0.5D0)*Q2(08,IELEM**3+1)) 
      CORNERQ(04,1) = (Q2(01,IELEM**3+1) + 3*Q2(04,IELEM**3+1) - 
     >   3*Q2(06,IELEM**3+1) - 3*Q2(07,IELEM**3+1) + 
     >   3**(0.5D0)*Q2(02,IELEM**3+1) + 3**(0.5D0)*Q2(03,IELEM**3+1) - 
     >   3**(0.5D0)*Q2(05,IELEM**3+1) - 3*3**(0.5D0)*Q2(08,IELEM**3+1)) 
      CORNERQ(01,2) = (Q2(01,IELEM**3+1) + 3*Q2(04,IELEM**3+1) - 
     >   3*Q2(06,IELEM**3+1) - 3*Q2(07,IELEM**3+1) - 
     >   3**(0.5D0)*Q2(02,IELEM**3+1) - 3**(0.5D0)*Q2(03,IELEM**3+1) + 
     >   3**(0.5D0)*Q2(05,IELEM**3+1) + 3*3**(0.5D0)*Q2(08,IELEM**3+1)) 
      CORNERQ(02,2) = (Q2(01,IELEM**3+1) - 3*Q2(04,IELEM**3+1) + 
     >   3*Q2(06,IELEM**3+1) - 3*Q2(07,IELEM**3+1) + 
     >   3**(0.5D0)*Q2(02,IELEM**3+1) - 3**(0.5D0)*Q2(03,IELEM**3+1) + 
     >   3**(0.5D0)*Q2(05,IELEM**3+1) - 3*3**(0.5D0)*Q2(08,IELEM**3+1)) 
      CORNERQ(03,2) = (Q2(01,IELEM**3+1) - 3*Q2(04,IELEM**3+1) - 
     >   3*Q2(06,IELEM**3+1) + 3*Q2(07,IELEM**3+1) - 
     >   3**(0.5D0)*Q2(02,IELEM**3+1) + 3**(0.5D0)*Q2(03,IELEM**3+1) + 
     >   3**(0.5D0)*Q2(05,IELEM**3+1) - 3*3**(0.5D0)*Q2(08,IELEM**3+1)) 
      CORNERQ(04,2) = (Q2(01,IELEM**3+1) + 3*Q2(04,IELEM**3+1) + 
     >   3*Q2(06,IELEM**3+1) + 3*Q2(07,IELEM**3+1) + 
     >   3**(0.5D0)*Q2(02,IELEM**3+1) + 3**(0.5D0)*Q2(03,IELEM**3+1) + 
     >   3**(0.5D0)*Q2(05,IELEM**3+1) + 3*3**(0.5D0)*Q2(08,IELEM**3+1)) 
*** ---------------------------------------------------------------- ***
*** ---------------------------------------------------------------- ***
*** ---------------------------------------------------------------- ***
      ELSEIF(IELEM.EQ.3)THEN
*** ---------------------------------------------------------------- ***
*** ---------------------------------------------------------------- ***
*** ---------------------------------------------------------------- ***
      CORNERQ(1,1) = (Q2(01,IELEM**3+1) + 3*Q2(05,IELEM**3+1) + 
     >   5*Q2(09,IELEM**3+1) + 3*Q2(11,IELEM**3+1) + 
     >   3*Q2(13,IELEM**3+1) + 5*Q2(21,IELEM**3+1) + 
     >   5*Q2(25,IELEM**3+1) - 3**(0.5D0)*Q2(02,IELEM**3+1) - 
     >   3**(0.5D0)*Q2(04,IELEM**3+1) + 
     >   5**(0.5D0)*Q2(03,IELEM**3+1) + 
     >   5**(0.5D0)*Q2(07,IELEM**3+1) - 
     >   3**(0.5D0)*Q2(10,IELEM**3+1) - 
     >   3*3**(0.5D0)*Q2(14,IELEM**3+1) + 
     >   3*5**(0.5D0)*Q2(15,IELEM**3+1) - 
     >   5*3**(0.5D0)*Q2(18,IELEM**3+1) - 
     >   15**(0.5D0)*Q2(06,IELEM**3+1) + 
     >   3*5**(0.5D0)*Q2(17,IELEM**3+1) - 
     >   15**(0.5D0)*Q2(08,IELEM**3+1) + 
     >   5**(0.5D0)*Q2(19,IELEM**3+1) - 
     >   5*3**(0.5D0)*Q2(24,IELEM**3+1) - 
     >   15**(0.5D0)*Q2(12,IELEM**3+1) + 
     >   3*5**(0.5D0)*Q2(23,IELEM**3+1) - 
     >   5*3**(0.5D0)*Q2(26,IELEM**3+1) - 
     >   15**(0.5D0)*Q2(16,IELEM**3+1) + 
     >   5*5**(0.5D0)*Q2(27,IELEM**3+1) - 
     >   15**(0.5D0)*Q2(20,IELEM**3+1) - 
     >   15**(0.5D0)*Q2(22,IELEM**3+1)) 
      CORNERQ(2,1) = (Q2(01,IELEM**3+1) - 
     >   (5*Q2(09,IELEM**3+1))/2 + 3*Q2(13,IELEM**3+1) - 
     >   (5*Q2(21,IELEM**3+1))/2 + 5*Q2(25,IELEM**3+1) - 
     >   3**(0.5D0)*Q2(04,IELEM**3+1) - 
     >   (5**(0.5D0)*Q2(03,IELEM**3+1))/2 + 
     >   5**(0.5D0)*Q2(07,IELEM**3+1) - 
     >   3**(0.5D0)*Q2(10,IELEM**3+1) - 
     >   (3*5**(0.5D0)*Q2(15,IELEM**3+1))/2 + 
     >   (5*3**(0.5D0)*Q2(18,IELEM**3+1))/2 + 
     >   (15**(0.5D0)*Q2(06,IELEM**3+1))/2 + 
     >   5**(0.5D0)*Q2(19,IELEM**3+1) + 
     >   (5*3**(0.5D0)*Q2(24,IELEM**3+1))/2 + 
     >   (15**(0.5D0)*Q2(12,IELEM**3+1))/2 - 
     >   15**(0.5D0)*Q2(16,IELEM**3+1) - 
     >   (5*5**(0.5D0)*Q2(27,IELEM**3+1))/2 - 
     >   15**(0.5D0)*Q2(22,IELEM**3+1)) 
      CORNERQ(3,1) = (Q2(01,IELEM**3+1) - 3*Q2(05,IELEM**3+1) + 
     >   5*Q2(09,IELEM**3+1) - 3*Q2(11,IELEM**3+1) + 
     >   3*Q2(13,IELEM**3+1) + 5*Q2(21,IELEM**3+1) + 
     >   5*Q2(25,IELEM**3+1) + 3**(0.5D0)*Q2(02,IELEM**3+1) - 
     >   3**(0.5D0)*Q2(04,IELEM**3+1) + 5**(0.5D0)*Q2(03,IELEM**3+1) + 
     >   5**(0.5D0)*Q2(07,IELEM**3+1) - 3**(0.5D0)*Q2(10,IELEM**3+1) + 
     >   3*3**(0.5D0)*Q2(14,IELEM**3+1) + 
     >   3*5**(0.5D0)*Q2(15,IELEM**3+1) - 
     >   5*3**(0.5D0)*Q2(18,IELEM**3+1) - 
     >   15**(0.5D0)*Q2(06,IELEM**3+1) - 
     >   3*5**(0.5D0)*Q2(17,IELEM**3+1) + 
     >   15**(0.5D0)*Q2(08,IELEM**3+1) + 
     >   5**(0.5D0)*Q2(19,IELEM**3+1) - 
     >   5*3**(0.5D0)*Q2(24,IELEM**3+1) - 
     >   15**(0.5D0)*Q2(12,IELEM**3+1) - 
     >   3*5**(0.5D0)*Q2(23,IELEM**3+1) + 
     >   5*3**(0.5D0)*Q2(26,IELEM**3+1) - 
     >   15**(0.5D0)*Q2(16,IELEM**3+1) + 
     >   5*5**(0.5D0)*Q2(27,IELEM**3+1) + 
     >   15**(0.5D0)*Q2(20,IELEM**3+1) - 
     >   15**(0.5D0)*Q2(22,IELEM**3+1)) 
      CORNERQ(4,1) = (Q2(01,IELEM**3+1) - (5*Q2(09,IELEM**3+1))/2 + 
     >   3*Q2(11,IELEM**3+1) + 5*Q2(21,IELEM**3+1) - 
     >   (5*Q2(25,IELEM**3+1))/2 - 3**(0.5D0)*Q2(02,IELEM**3+1) + 
     >   5**(0.5D0)*Q2(03,IELEM**3+1) - 
     >   (5**(0.5D0)*Q2(07,IELEM**3+1))/2 - 
     >   3**(0.5D0)*Q2(10,IELEM**3+1) + 
     >   (5*3**(0.5D0)*Q2(18,IELEM**3+1))/2 - 
     >   (3*5**(0.5D0)*Q2(17,IELEM**3+1))/2 + 
     >   (15**(0.5D0)*Q2(08,IELEM**3+1))/2 + 
     >   5**(0.5D0)*Q2(19,IELEM**3+1) - 
     >   15**(0.5D0)*Q2(12,IELEM**3+1) + 
     >   (5*3**(0.5D0)*Q2(26,IELEM**3+1))/2 + 
     >   (15**(0.5D0)*Q2(16,IELEM**3+1))/2 - 
     >   (5*5**(0.5D0)*Q2(27,IELEM**3+1))/2 - 
     >   15**(0.5D0)*Q2(20,IELEM**3+1)) 
      CORNERQ(5,1) = (Q2(01,IELEM**3+1) + 
     >   (5*Q2(09,IELEM**3+1))/4 - (5*Q2(21,IELEM**3+1))/2 - 
     >   (5*Q2(25,IELEM**3+1))/2 - (5**(0.5D0)*Q2(03,IELEM**3+1))/2 - 
     >   (5**(0.5D0)*Q2(07,IELEM**3+1))/2 - 
     >   3**(0.5D0)*Q2(10,IELEM**3+1) - 
     >   (5*3**(0.5D0)*Q2(18,IELEM**3+1))/4 + 
     >   5**(0.5D0)*Q2(19,IELEM**3+1) + 
     >   (15**(0.5D0)*Q2(12,IELEM**3+1))/2 + 
     >   (15**(0.5D0)*Q2(16,IELEM**3+1))/2 + 
     >   (5*5**(0.5D0)*Q2(27,IELEM**3+1))/4) 
      CORNERQ(6,1) = (Q2(01,IELEM**3+1) - (5*Q2(09,IELEM**3+1))/2 - 
     >   3*Q2(11,IELEM**3+1) + 5*Q2(21,IELEM**3+1) - 
     >   (5*Q2(25,IELEM**3+1))/2 + 3**(0.5D0)*Q2(02,IELEM**3+1) + 
     >   5**(0.5D0)*Q2(03,IELEM**3+1) - 
     >   (5**(0.5D0)*Q2(07,IELEM**3+1))/2 - 
     >   3**(0.5D0)*Q2(10,IELEM**3+1) + 
     >   (5*3**(0.5D0)*Q2(18,IELEM**3+1))/2 + 
     >   (3*5**(0.5D0)*Q2(17,IELEM**3+1))/2 - 
     >   (15**(0.5D0)*Q2(08,IELEM**3+1))/2 + 
     >   5**(0.5D0)*Q2(19,IELEM**3+1) - 
     >   15**(0.5D0)*Q2(12,IELEM**3+1) - 
     >   (5*3**(0.5D0)*Q2(26,IELEM**3+1))/2 + 
     >   (15**(0.5D0)*Q2(16,IELEM**3+1))/2 - 
     >   (5*5**(0.5D0)*Q2(27,IELEM**3+1))/2 + 
     >   15**(0.5D0)*Q2(20,IELEM**3+1)) 
      CORNERQ(7,1) = (Q2(01,IELEM**3+1) - 3*Q2(05,IELEM**3+1) + 
     >   5*Q2(09,IELEM**3+1) + 3*Q2(11,IELEM**3+1) - 
     >   3*Q2(13,IELEM**3+1) + 5*Q2(21,IELEM**3+1) + 
     >   5*Q2(25,IELEM**3+1) - 3**(0.5D0)*Q2(02,IELEM**3+1) + 
     >   3**(0.5D0)*Q2(04,IELEM**3+1) + 
     >   5**(0.5D0)*Q2(03,IELEM**3+1) + 
     >   5**(0.5D0)*Q2(07,IELEM**3+1) - 
     >   3**(0.5D0)*Q2(10,IELEM**3+1) + 
     >   3*3**(0.5D0)*Q2(14,IELEM**3+1) - 
     >   3*5**(0.5D0)*Q2(15,IELEM**3+1) - 
     >   5*3**(0.5D0)*Q2(18,IELEM**3+1) + 
     >   15**(0.5D0)*Q2(06,IELEM**3+1) + 
     >   3*5**(0.5D0)*Q2(17,IELEM**3+1) - 
     >   15**(0.5D0)*Q2(08,IELEM**3+1) + 
     >   5**(0.5D0)*Q2(19,IELEM**3+1) + 
     >   5*3**(0.5D0)*Q2(24,IELEM**3+1) - 
     >   15**(0.5D0)*Q2(12,IELEM**3+1) - 
     >   3*5**(0.5D0)*Q2(23,IELEM**3+1) - 
     >   5*3**(0.5D0)*Q2(26,IELEM**3+1) - 
     >   15**(0.5D0)*Q2(16,IELEM**3+1) + 
     >   5*5**(0.5D0)*Q2(27,IELEM**3+1) - 
     >   15**(0.5D0)*Q2(20,IELEM**3+1) + 
     >   15**(0.5D0)*Q2(22,IELEM**3+1)) 
      CORNERQ(8,1) = (Q2(01,IELEM**3+1) - 
     >   (5*Q2(09,IELEM**3+1))/2 - 3*Q2(13,IELEM**3+1) - 
     >   (5*Q2(21,IELEM**3+1))/2 + 5*Q2(25,IELEM**3+1) + 
     >   3**(0.5D0)*Q2(04,IELEM**3+1) - 
     >   (5**(0.5D0)*Q2(03,IELEM**3+1))/2 + 
     >   5**(0.5D0)*Q2(07,IELEM**3+1) - 
     >   3**(0.5D0)*Q2(10,IELEM**3+1) + 
     >   (3*5**(0.5D0)*Q2(15,IELEM**3+1))/2 + 
     >   (5*3**(0.5D0)*Q2(18,IELEM**3+1))/2 - 
     >   (15**(0.5D0)*Q2(06,IELEM**3+1))/2 + 
     >   5**(0.5D0)*Q2(19,IELEM**3+1) - 
     >   (5*3**(0.5D0)*Q2(24,IELEM**3+1))/2 + 
     >   (15**(0.5D0)*Q2(12,IELEM**3+1))/2 - 
     >   15**(0.5D0)*Q2(16,IELEM**3+1) - 
     >   (5*5**(0.5D0)*Q2(27,IELEM**3+1))/2 + 
     >   15**(0.5D0)*Q2(22,IELEM**3+1)) 
      CORNERQ(9,1) = (Q2(01,IELEM**3+1) + 
     >   3*Q2(05,IELEM**3+1) + 5*Q2(09,IELEM**3+1) - 
     >   3*Q2(11,IELEM**3+1) - 3*Q2(13,IELEM**3+1) + 
     >   5*Q2(21,IELEM**3+1) + 5*Q2(25,IELEM**3+1) + 
     >   3**(0.5D0)*Q2(02,IELEM**3+1) + 
     >   3**(0.5D0)*Q2(04,IELEM**3+1) + 
     >   5**(0.5D0)*Q2(03,IELEM**3+1) + 
     >   5**(0.5D0)*Q2(07,IELEM**3+1) - 
     >   3**(0.5D0)*Q2(10,IELEM**3+1) - 
     >   3*3**(0.5D0)*Q2(14,IELEM**3+1) - 
     >   3*5**(0.5D0)*Q2(15,IELEM**3+1) - 
     >   5*3**(0.5D0)*Q2(18,IELEM**3+1) + 
     >   15**(0.5D0)*Q2(06,IELEM**3+1) - 
     >   3*5**(0.5D0)*Q2(17,IELEM**3+1) + 
     >   15**(0.5D0)*Q2(08,IELEM**3+1) + 
     >   5**(0.5D0)*Q2(19,IELEM**3+1) + 
     >   5*3**(0.5D0)*Q2(24,IELEM**3+1) - 
     >   15**(0.5D0)*Q2(12,IELEM**3+1) + 
     >   3*5**(0.5D0)*Q2(23,IELEM**3+1) + 
     >   5*3**(0.5D0)*Q2(26,IELEM**3+1) - 
     >   15**(0.5D0)*Q2(16,IELEM**3+1) + 
     >   5*5**(0.5D0)*Q2(27,IELEM**3+1) + 
     >   15**(0.5D0)*Q2(20,IELEM**3+1) + 
     >   15**(0.5D0)*Q2(22,IELEM**3+1)) 
      CORNERQ(1,2) = (Q2(01,IELEM**3+1) + 3*Q2(05,IELEM**3+1) + 
     >   5*Q2(09,IELEM**3+1) - (5*Q2(21,IELEM**3+1))/2 - 
     >   (5*Q2(25,IELEM**3+1))/2 - 3**(0.5D0)*Q2(02,IELEM**3+1) - 
     >   3**(0.5D0)*Q2(04,IELEM**3+1) + 
     >   5**(0.5D0)*Q2(03,IELEM**3+1) + 
     >   5**(0.5D0)*Q2(07,IELEM**3+1) - 15**(0.5D0)*Q2(06,IELEM**3+1) - 
     >   15**(0.5D0)*Q2(08,IELEM**3+1) - 
     >   (5**(0.5D0)*Q2(19,IELEM**3+1))/2 + 
     >   (5*3**(0.5D0)*Q2(24,IELEM**3+1))/2 - 
     >   (3*5**(0.5D0)*Q2(23,IELEM**3+1))/2 + (5*3**(0.5D0)*
     >   Q2(26,IELEM**3+1))/2 - (5*5**(0.5D0)*Q2(27,IELEM**3+1))/2 + 
     >   (15**(0.5D0)*Q2(20,IELEM**3+1))/2 + (15**(0.5D0)*
     >   Q2(22,IELEM**3+1))/2) 
      CORNERQ(2,2) = (Q2(01,IELEM**3+1) - (5*Q2(09,IELEM**3+1))/2 + 
     >   (5*Q2(21,IELEM**3+1))/4 - (5*Q2(25,IELEM**3+1))/2 - 
     >   3**(0.5D0)*Q2(04,IELEM**3+1) - (5**(0.5D0)*
     >   Q2(03,IELEM**3+1))/2 + 5**(0.5D0)*Q2(07,IELEM**3+1) + 
     >   (15**(0.5D0)*Q2(06,IELEM**3+1))/2 - (5**(0.5D0)*
     >   Q2(19,IELEM**3+1))/2 - (5*3**(0.5D0)*Q2(24,IELEM**3+1))/4 + 
     >   (5*5**(0.5D0)*Q2(27,IELEM**3+1))/4 + (15**(0.5D0)*
     >   Q2(22,IELEM**3+1))/2) 
      CORNERQ(3,2) = (Q2(01,IELEM**3+1) - 3*Q2(05,IELEM**3+1) + 
     >   5*Q2(09,IELEM**3+1) - (5*Q2(21,IELEM**3+1))/2 - 
     >   (5*Q2(25,IELEM**3+1))/2 + 3**(0.5D0)*Q2(02,IELEM**3+1) - 
     >   3**(0.5D0)*Q2(04,IELEM**3+1) + 5**(0.5D0)*Q2(03,IELEM**3+1) + 
     >   5**(0.5D0)*Q2(07,IELEM**3+1) - 15**(0.5D0)*
     >   Q2(06,IELEM**3+1) + 15**(0.5D0)*Q2(08,IELEM**3+1) - 
     >   (5**(0.5D0)*Q2(19,IELEM**3+1))/2 + (5*3**(0.5D0)*
     >   Q2(24,IELEM**3+1))/2 + (3*5**(0.5D0)*
     >   Q2(23,IELEM**3+1))/2 - (5*3**(0.5D0)*Q2(26,IELEM**3+1))/2 - 
     >   (5*5**(0.5D0)*Q2(27,IELEM**3+1))/2 - (15**(0.5D0)*
     >   Q2(20,IELEM**3+1))/2 + (15**(0.5D0)*Q2(22,IELEM**3+1))/2) 
      CORNERQ(4,2) = (Q2(01,IELEM**3+1) - (5*Q2(09,IELEM**3+1))/2 - 
     >   (5*Q2(21,IELEM**3+1))/2 + (5*Q2(25,IELEM**3+1))/4 - 
     >   3**(0.5D0)*Q2(02,IELEM**3+1) + 5**(0.5D0)*Q2(03,IELEM**3+1) - 
     >   (5**(0.5D0)*Q2(07,IELEM**3+1))/2 + (15**(0.5D0)*
     >   Q2(08,IELEM**3+1))/2 - (5**(0.5D0)*Q2(19,IELEM**3+1))/2 - 
     >   (5*3**(0.5D0)*Q2(26,IELEM**3+1))/4 + (5*5**(0.5D0)*
     >   Q2(27,IELEM**3+1))/4 + (15**(0.5D0)*Q2(20,IELEM**3+1))/2) 
      CORNERQ(5,2) = (Q2(01,IELEM**3+1) + (5*Q2(09,IELEM**3+1))/4 + 
     >   (5*Q2(21,IELEM**3+1))/4 + (5*Q2(25,IELEM**3+1))/4 - 
     >   (5**(0.5D0)*Q2(03,IELEM**3+1))/2 - (5**(0.5D0)*
     >   Q2(07,IELEM**3+1))/2 - (5**(0.5D0)*
     >   Q2(19,IELEM**3+1))/2 - (5*5**(0.5D0)*Q2(27,IELEM**3+1))/8) 
      CORNERQ(6,2) = (Q2(01,IELEM**3+1) - (5*Q2(09,IELEM**3+1))/2 - 
     >   (5*Q2(21,IELEM**3+1))/2 + (5*Q2(25,IELEM**3+1))/4 + 
     >   3**(0.5D0)*Q2(02,IELEM**3+1) + 5**(0.5D0)*Q2(03,IELEM**3+1) - 
     >   (5**(0.5D0)*Q2(07,IELEM**3+1))/2 - (15**(0.5D0)*
     >   Q2(08,IELEM**3+1))/2 - (5**(0.5D0)*Q2(19,IELEM**3+1))/2 + 
     >   (5*3**(0.5D0)*Q2(26,IELEM**3+1))/4 + (5*5**(0.5D0)*
     >   Q2(27,IELEM**3+1))/4 - (15**(0.5D0)*Q2(20,IELEM**3+1))/2) 
      CORNERQ(7,2) = (Q2(01,IELEM**3+1) - 3*Q2(05,IELEM**3+1) + 
     >   5*Q2(09,IELEM**3+1) - (5*Q2(21,IELEM**3+1))/2 - 
     >   (5*Q2(25,IELEM**3+1))/2 - 3**(0.5D0)*Q2(02,IELEM**3+1) + 
     >   3**(0.5D0)*Q2(04,IELEM**3+1) + 5**(0.5D0)*Q2(03,IELEM**3+1) + 
     >   5**(0.5D0)*Q2(07,IELEM**3+1) + 15**(0.5D0)*
     >   Q2(06,IELEM**3+1) - 15**(0.5D0)*Q2(08,IELEM**3+1) - 
     >   (5**(0.5D0)*Q2(19,IELEM**3+1))/2 - (5*3**(0.5D0)*
     >   Q2(24,IELEM**3+1))/2 + (3*5**(0.5D0)*Q2(23,IELEM**3+1))/2 + 
     >   (5*3**(0.5D0)*Q2(26,IELEM**3+1))/2 - (5*5**(0.5D0)*
     >   Q2(27,IELEM**3+1))/2 + (15**(0.5D0)*Q2(20,IELEM**3+1))/2 - 
     >   (15**(0.5D0)*Q2(22,IELEM**3+1))/2) 
      CORNERQ(8,2) = (Q2(01,IELEM**3+1) - (5*Q2(09,IELEM**3+1))/2 + 
     >   (5*Q2(21,IELEM**3+1))/4 - (5*Q2(25,IELEM**3+1))/2 + 
     >   3**(0.5D0)*Q2(04,IELEM**3+1) - (5**(0.5D0)*
     >   Q2(03,IELEM**3+1))/2 + 5**(0.5D0)*Q2(07,IELEM**3+1) - 
     >   (15**(0.5D0)*Q2(06,IELEM**3+1))/2 - (5**(0.5D0)*
     >   Q2(19,IELEM**3+1))/2 + (5*3**(0.5D0)*
     >   Q2(24,IELEM**3+1))/4 + (5*5**(0.5D0)*
     >   Q2(27,IELEM**3+1))/4 - (15**(0.5D0)*Q2(22,IELEM**3+1))/2) 
      CORNERQ(9,2) = (Q2(01,IELEM**3+1) + 3*Q2(05,IELEM**3+1) + 
     >   5*Q2(09,IELEM**3+1) - (5*Q2(21,IELEM**3+1))/2 - 
     >   (5*Q2(25,IELEM**3+1))/2 + 3**(0.5D0)*Q2(02,IELEM**3+1) + 
     >   3**(0.5D0)*Q2(04,IELEM**3+1) + 5**(0.5D0)*Q2(03,IELEM**3+1) + 
     >   5**(0.5D0)*Q2(07,IELEM**3+1) + 15**(0.5D0)*Q2(06,IELEM**3+1) + 
     >   15**(0.5D0)*Q2(08,IELEM**3+1) - (5**(0.5D0)*
     >   Q2(19,IELEM**3+1))/2 - (5*3**(0.5D0)*Q2(24,IELEM**3+1))/2 - 
     >   (3*5**(0.5D0)*Q2(23,IELEM**3+1))/2 - (5*3**(0.5D0)*
     >   Q2(26,IELEM**3+1))/2 - (5*5**(0.5D0)*Q2(27,IELEM**3+1))/2 - 
     >   (15**(0.5D0)*Q2(20,IELEM**3+1))/2 - (15**(0.5D0)*
     >   Q2(22,IELEM**3+1))/2) 
      CORNERQ(1,3) = (Q2(01,IELEM**3+1) + 3*Q2(05,IELEM**3+1) + 
     >   5*Q2(09,IELEM**3+1) - 3*Q2(11,IELEM**3+1) - 3*
     >   Q2(13,IELEM**3+1) + 5*Q2(21,IELEM**3+1) + 5*
     >   Q2(25,IELEM**3+1) - 3**(0.5D0)*Q2(02,IELEM**3+1) - 
     >   3**(0.5D0)*Q2(04,IELEM**3+1) + 5**(0.5D0)*Q2(03,IELEM**3+1) + 
     >   5**(0.5D0)*Q2(07,IELEM**3+1) + 3**(0.5D0)*Q2(10,IELEM**3+1) + 
     >   3*3**(0.5D0)*Q2(14,IELEM**3+1) - 3*5**(0.5D0)*
     >   Q2(15,IELEM**3+1) + 5*3**(0.5D0)*Q2(18,IELEM**3+1) - 
     >   15**(0.5D0)*Q2(06,IELEM**3+1) - 3*5**(0.5D0)*
     >   Q2(17,IELEM**3+1) - 15**(0.5D0)*Q2(08,IELEM**3+1) + 
     >   5**(0.5D0)*Q2(19,IELEM**3+1) - 5*3**(0.5D0)*
     >   Q2(24,IELEM**3+1) + 15**(0.5D0)*Q2(12,IELEM**3+1) + 
     >   3*5**(0.5D0)*Q2(23,IELEM**3+1) - 5*3**(0.5D0)*
     >   Q2(26,IELEM**3+1) + 15**(0.5D0)*Q2(16,IELEM**3+1) + 
     >   5*5**(0.5D0)*Q2(27,IELEM**3+1) - 15**(0.5D0)*
     >   Q2(20,IELEM**3+1) - 15**(0.5D0)*Q2(22,IELEM**3+1)) 
      CORNERQ(2,3) = (Q2(01,IELEM**3+1) - (5*Q2(09,IELEM**3+1))/2 - 
     >   3*Q2(13,IELEM**3+1) - (5*Q2(21,IELEM**3+1))/2 + 
     >   5*Q2(25,IELEM**3+1) - 3**(0.5D0)*Q2(04,IELEM**3+1) - 
     >   (5**(0.5D0)*Q2(03,IELEM**3+1))/2 + 5**(0.5D0)*
     >   Q2(07,IELEM**3+1) + 3**(0.5D0)*Q2(10,IELEM**3+1) + 
     >   (3*5**(0.5D0)*Q2(15,IELEM**3+1))/2 - (5*3**(0.5D0)*
     >   Q2(18,IELEM**3+1))/2 + (15**(0.5D0)*Q2(06,IELEM**3+1))/2 + 
     >   5**(0.5D0)*Q2(19,IELEM**3+1) + (5*3**(0.5D0)*
     >   Q2(24,IELEM**3+1))/2 - (15**(0.5D0)*Q2(12,IELEM**3+1))/2 + 
     >   15**(0.5D0)*Q2(16,IELEM**3+1) - (5*5**(0.5D0)*
     >   Q2(27,IELEM**3+1))/2 - 15**(0.5D0)*Q2(22,IELEM**3+1)) 
      CORNERQ(3,3) = (Q2(01,IELEM**3+1) - 3*Q2(05,IELEM**3+1) + 
     >   5*Q2(09,IELEM**3+1) + 3*Q2(11,IELEM**3+1) - 3*
     >   Q2(13,IELEM**3+1) + 5*Q2(21,IELEM**3+1) + 5*
     >   Q2(25,IELEM**3+1) + 3**(0.5D0)*Q2(02,IELEM**3+1) - 
     >   3**(0.5D0)*Q2(04,IELEM**3+1) + 5**(0.5D0)*
     >   Q2(03,IELEM**3+1) + 5**(0.5D0)*Q2(07,IELEM**3+1) + 
     >   3**(0.5D0)*Q2(10,IELEM**3+1) - 3*3**(0.5D0)*
     >   Q2(14,IELEM**3+1) - 3*5**(0.5D0)*Q2(15,IELEM**3+1) + 
     >   5*3**(0.5D0)*Q2(18,IELEM**3+1) - 15**(0.5D0)*
     >   Q2(06,IELEM**3+1) + 3*5**(0.5D0)*Q2(17,IELEM**3+1) + 
     >   15**(0.5D0)*Q2(08,IELEM**3+1) + 5**(0.5D0)*
     >   Q2(19,IELEM**3+1) - 5*3**(0.5D0)*Q2(24,IELEM**3+1) + 
     >   15**(0.5D0)*Q2(12,IELEM**3+1) - 3*5**(0.5D0)*
     >   Q2(23,IELEM**3+1) + 5*3**(0.5D0)*Q2(26,IELEM**3+1) + 
     >   15**(0.5D0)*Q2(16,IELEM**3+1) + 5*5**(0.5D0)*
     >   Q2(27,IELEM**3+1) + 15**(0.5D0)*Q2(20,IELEM**3+1) - 
     >   15**(0.5D0)*Q2(22,IELEM**3+1)) 
      CORNERQ(4,3) = (Q2(01,IELEM**3+1) - (5*Q2(09,IELEM**3+1))/2 - 
     >   3*Q2(11,IELEM**3+1) + 5*Q2(21,IELEM**3+1) - (5*
     >   Q2(25,IELEM**3+1))/2 - 3**(0.5D0)*Q2(02,IELEM**3+1) + 
     >   5**(0.5D0)*Q2(03,IELEM**3+1) - (5**(0.5D0)*
     >   Q2(07,IELEM**3+1))/2 + 3**(0.5D0)*Q2(10,IELEM**3+1) - 
     >   (5*3**(0.5D0)*Q2(18,IELEM**3+1))/2 + (3*5**(0.5D0)*
     >   Q2(17,IELEM**3+1))/2 + (15**(0.5D0)*
     >   Q2(08,IELEM**3+1))/2 + 5**(0.5D0)*Q2(19,IELEM**3+1) + 
     >   15**(0.5D0)*Q2(12,IELEM**3+1) + (5*3**(0.5D0)*
     >   Q2(26,IELEM**3+1))/2 - (15**(0.5D0)*
     >   Q2(16,IELEM**3+1))/2 - (5*5**(0.5D0)*
     >   Q2(27,IELEM**3+1))/2 - 15**(0.5D0)*Q2(20,IELEM**3+1)) 
      CORNERQ(5,3) = (Q2(01,IELEM**3+1) + (5*Q2(09,IELEM**3+1))/4 - 
     >   (5*Q2(21,IELEM**3+1))/2 - (5*Q2(25,IELEM**3+1))/2 - 
     >   (5**(0.5D0)*Q2(03,IELEM**3+1))/2 - (5**(0.5D0)*
     >   Q2(07,IELEM**3+1))/2 + 3**(0.5D0)*Q2(10,IELEM**3+1) + 
     >   (5*3**(0.5D0)*Q2(18,IELEM**3+1))/4 + 5**(0.5D0)*
     >   Q2(19,IELEM**3+1) - (15**(0.5D0)*Q2(12,IELEM**3+1))/2 - 
     >   (15**(0.5D0)*Q2(16,IELEM**3+1))/2 + (5*5**(0.5D0)*
     >   Q2(27,IELEM**3+1))/4) 
      CORNERQ(6,3) = (Q2(01,IELEM**3+1) - (5*Q2(09,IELEM**3+1))/2 + 
     >   3*Q2(11,IELEM**3+1) + 5*Q2(21,IELEM**3+1) - 
     >   (5*Q2(25,IELEM**3+1))/2 + 3**(0.5D0)*Q2(02,IELEM**3+1) + 
     >   5**(0.5D0)*Q2(03,IELEM**3+1) - (5**(0.5D0)*
     >   Q2(07,IELEM**3+1))/2 + 3**(0.5D0)*Q2(10,IELEM**3+1) - 
     >   (5*3**(0.5D0)*Q2(18,IELEM**3+1))/2 - (3*5**(0.5D0)*
     >   Q2(17,IELEM**3+1))/2 - (15**(0.5D0)*Q2(08,IELEM**3+1))/2 + 
     >   5**(0.5D0)*Q2(19,IELEM**3+1) + 15**(0.5D0)*
     >   Q2(12,IELEM**3+1) - (5*3**(0.5D0)*Q2(26,IELEM**3+1))/2 - 
     >   (15**(0.5D0)*Q2(16,IELEM**3+1))/2 - (5*5**(0.5D0)*
     >   Q2(27,IELEM**3+1))/2 + 15**(0.5D0)*Q2(20,IELEM**3+1)) 
      CORNERQ(7,3) = (Q2(01,IELEM**3+1) - 3*Q2(05,IELEM**3+1) + 
     >   5*Q2(09,IELEM**3+1) - 3*Q2(11,IELEM**3+1) + 
     >   3*Q2(13,IELEM**3+1) + 5*Q2(21,IELEM**3+1) + 
     >   5*Q2(25,IELEM**3+1) - 3**(0.5D0)*Q2(02,IELEM**3+1) + 
     >   3**(0.5D0)*Q2(04,IELEM**3+1) + 5**(0.5D0)*Q2(03,IELEM**3+1) + 
     >   5**(0.5D0)*Q2(07,IELEM**3+1) + 3**(0.5D0)*Q2(10,IELEM**3+1) - 
     >   3*3**(0.5D0)*Q2(14,IELEM**3+1) + 3*5**(0.5D0)*
     >   Q2(15,IELEM**3+1) + 5*3**(0.5D0)*Q2(18,IELEM**3+1) + 
     >   15**(0.5D0)*Q2(06,IELEM**3+1) - 3*5**(0.5D0)*
     >   Q2(17,IELEM**3+1) - 15**(0.5D0)*Q2(08,IELEM**3+1) + 
     >   5**(0.5D0)*Q2(19,IELEM**3+1) + 5*3**(0.5D0)*
     >   Q2(24,IELEM**3+1) + 15**(0.5D0)*Q2(12,IELEM**3+1) - 
     >   3*5**(0.5D0)*Q2(23,IELEM**3+1) - 5*3**(0.5D0)*
     >   Q2(26,IELEM**3+1) + 15**(0.5D0)*Q2(16,IELEM**3+1) + 
     >   5*5**(0.5D0)*Q2(27,IELEM**3+1) - 15**(0.5D0)*
     >   Q2(20,IELEM**3+1) + 15**(0.5D0)*Q2(22,IELEM**3+1)) 
      CORNERQ(8,3) = (Q2(01,IELEM**3+1) - (5*Q2(09,IELEM**3+1))/2 + 
     >   3*Q2(13,IELEM**3+1) - (5*Q2(21,IELEM**3+1))/2 + 
     >   5*Q2(25,IELEM**3+1) + 3**(0.5D0)*Q2(04,IELEM**3+1) - 
     >   (5**(0.5D0)*Q2(03,IELEM**3+1))/2 + 5**(0.5D0)*
     >   Q2(07,IELEM**3+1) + 3**(0.5D0)*Q2(10,IELEM**3+1) - 
     >   (3*5**(0.5D0)*Q2(15,IELEM**3+1))/2 - (5*3**(0.5D0)*
     >   Q2(18,IELEM**3+1))/2 - (15**(0.5D0)*Q2(06,IELEM**3+1))/2 + 
     >   5**(0.5D0)*Q2(19,IELEM**3+1) - (5*3**(0.5D0)*
     >   Q2(24,IELEM**3+1))/2 - (15**(0.5D0)*Q2(12,IELEM**3+1))/2 + 
     >   15**(0.5D0)*Q2(16,IELEM**3+1) - (5*5**(0.5D0)*
     >   Q2(27,IELEM**3+1))/2 + 15**(0.5D0)*Q2(22,IELEM**3+1)) 
      CORNERQ(9,3) = (Q2(01,IELEM**3+1) + 3*Q2(05,IELEM**3+1) + 
     >   5*Q2(09,IELEM**3+1) + 3*Q2(11,IELEM**3+1) + 
     >   3*Q2(13,IELEM**3+1) + 5*Q2(21,IELEM**3+1) + 
     >   5*Q2(25,IELEM**3+1) + 3**(0.5D0)*Q2(02,IELEM**3+1) + 
     >   3**(0.5D0)*Q2(04,IELEM**3+1) + 5**(0.5D0)*Q2(03,IELEM**3+1) + 
     >   5**(0.5D0)*Q2(07,IELEM**3+1) + 3**(0.5D0)*Q2(10,IELEM**3+1) + 
     >   3*3**(0.5D0)*Q2(14,IELEM**3+1) + 3*5**(0.5D0)*
     >   Q2(15,IELEM**3+1) + 5*3**(0.5D0)*Q2(18,IELEM**3+1) + 
     >   15**(0.5D0)*Q2(06,IELEM**3+1) + 3*5**(0.5D0)*
     >   Q2(17,IELEM**3+1) + 15**(0.5D0)*Q2(08,IELEM**3+1) + 
     >   5**(0.5D0)*Q2(19,IELEM**3+1) + 5*3**(0.5D0)*
     >   Q2(24,IELEM**3+1) + 15**(0.5D0)*Q2(12,IELEM**3+1) + 
     >   3*5**(0.5D0)*Q2(23,IELEM**3+1) + 5*3**(0.5D0)*
     >   Q2(26,IELEM**3+1) + 15**(0.5D0)*Q2(16,IELEM**3+1) + 
     >   5*5**(0.5D0)*Q2(27,IELEM**3+1) + 15**(0.5D0)*
     >   Q2(20,IELEM**3+1) + 15**(0.5D0)*Q2(22,IELEM**3+1)) 
*** ---------------------------------------------------------------- ***
*** ---------------------------------------------------------------- ***
*** ---------------------------------------------------------------- ***
      ELSEIF(IELEM.EQ.4)THEN
*** ---------------------------------------------------------------- ***
*** ---------------------------------------------------------------- ***
*** ---------------------------------------------------------------- ***
      CORNERQ(01,1) = (Q2(01,IELEM**3+1) + 3*Q2(06,IELEM**3+1) + 5*
     >   Q2(11,IELEM**3+1) + 7*Q2(16,IELEM**3+1) + 3*Q2(18,IELEM**3+1)
     >    + 3*Q2(21,IELEM**3+1) + 5*Q2(35,IELEM**3+1) + 5*Q2(41,IELEM*
     >   *3+1) + 7*Q2(52,IELEM**3+1) + 7*Q2(61,IELEM**3+1) - 3**(0.5D0
     >   )*Q2(02,IELEM**3+1) - 3**(0.5D0)*Q2(05,IELEM**3+1) + 5**(0.5D
     >   0)*Q2(03,IELEM**3+1) - 7**(0.5D0)*Q2(04,IELEM**3+1) + 5**(0.5
     >   D0)*Q2(09,IELEM**3+1) - 3**(0.5D0)*Q2(17,IELEM**3+1) - 7**(0.
     >   5D0)*Q2(13,IELEM**3+1) - 15**(0.5D0)*Q2(07,IELEM**3+1) - 3*3*
     >   *(0.5D0)*Q2(22,IELEM**3+1) - 15**(0.5D0)*Q2(10,IELEM**3+1) + 
     >   3*5**(0.5D0)*Q2(23,IELEM**3+1) + 21**(0.5D0)*Q2(08,IELEM**3+1
     >   ) - 5*3**(0.5D0)*Q2(27,IELEM**3+1) + 3*5**(0.5D0)*Q2(26,IELEM
     >   **3+1) - 3*7**(0.5D0)*Q2(24,IELEM**3+1) - 15**(0.5D0)*Q2(19,I
     >   ELEM**3+1) - 7*3**(0.5D0)*Q2(32,IELEM**3+1) + 21**(0.5D0)*Q2(
     >   14,IELEM**3+1) - 3*7**(0.5D0)*Q2(30,IELEM**3+1) + 5**(0.5D0)*
     >   Q2(33,IELEM**3+1) - 15**(0.5D0)*Q2(25,IELEM**3+1) + 21**(0.5D
     >   0)*Q2(20,IELEM**3+1) - 5*3**(0.5D0)*Q2(39,IELEM**3+1) + 3*5**
     >   (0.5D0)*Q2(38,IELEM**3+1) - 5*3**(0.5D0)*Q2(42,IELEM**3+1) - 
     >   35**(0.5D0)*Q2(12,IELEM**3+1) + 5*5**(0.5D0)*Q2(43,IELEM**3+1
     >   ) - 15**(0.5D0)*Q2(34,IELEM**3+1) + 21**(0.5D0)*Q2(29,IELEM**
     >   3+1) - 35**(0.5D0)*Q2(15,IELEM**3+1) - 5*7**(0.5D0)*Q2(44,IEL
     >   EM**3+1) - 15**(0.5D0)*Q2(37,IELEM**3+1) + 7*5**(0.5D0)*Q2(48
     >   ,IELEM**3+1) - 5*7**(0.5D0)*Q2(47,IELEM**3+1) - 7**(0.5D0)*Q2
     >   (49,IELEM**3+1) - 7*3**(0.5D0)*Q2(56,IELEM**3+1) - 3*7**(0.5D
     >   0)*Q2(54,IELEM**3+1) - 7*3**(0.5D0)*Q2(62,IELEM**3+1) + 7*5**
     >   (0.5D0)*Q2(60,IELEM**3+1) - 5*7**(0.5D0)*Q2(59,IELEM**3+1) + 
     >   7*5**(0.5D0)*Q2(63,IELEM**3+1) - 7*7**(0.5D0)*Q2(64,IELEM**3+
     >   1) + 21**(0.5D0)*Q2(50,IELEM**3+1) - 35**(0.5D0)*Q2(36,IELEM*
     >   *3+1) + 21**(0.5D0)*Q2(53,IELEM**3+1) - 35**(0.5D0)*Q2(45,IEL
     >   EM**3+1) - 35**(0.5D0)*Q2(51,IELEM**3+1) - 35**(0.5D0)*Q2(57,
     >   IELEM**3+1) + 105**(0.5D0)*Q2(28,IELEM**3+1) + 105**(0.5D0)*Q
     >   2(31,IELEM**3+1) + 105**(0.5D0)*Q2(40,IELEM**3+1) + 105**(0.5
     >   D0)*Q2(46,IELEM**3+1) + 105**(0.5D0)*Q2(55,IELEM**3+1) + 105*
     >   *(0.5D0)*Q2(58,IELEM**3+1)) 
      CORNERQ(02,1) = (Q2(01,IELEM**3+1) + Q2(06,IELEM**3+1) - (5*Q
     >   2(11,IELEM**3+1))/3 - (77*Q2(16,IELEM**3+1))/27 + Q2(18,IELEM
     >   **3+1) + 3*Q2(21,IELEM**3+1) - (5*Q2(35,IELEM**3+1))/3 + 5*Q2
     >   (41,IELEM**3+1) - (77*Q2(52,IELEM**3+1))/27 + 7*Q2(61,IELEM**
     >   3+1) - (3**(0.5D0)*Q2(02,IELEM**3+1))/3 - 3**(0.5D0)*Q2(05,IE
     >   LEM**3+1) - (5**(0.5D0)*Q2(03,IELEM**3+1))/3 + (11*7**(0.5D0)
     >   *Q2(04,IELEM**3+1))/27 + 5**(0.5D0)*Q2(09,IELEM**3+1) - 3**(0
     >   .5D0)*Q2(17,IELEM**3+1) - 7**(0.5D0)*Q2(13,IELEM**3+1) + (15*
     >   *(0.5D0)*Q2(07,IELEM**3+1))/3 - 3**(0.5D0)*Q2(22,IELEM**3+1) 
     >   - (15**(0.5D0)*Q2(10,IELEM**3+1))/3 - 5**(0.5D0)*Q2(23,IELEM*
     >   *3+1) - (11*21**(0.5D0)*Q2(08,IELEM**3+1))/27 + (5*3**(0.5D0)
     >   *Q2(27,IELEM**3+1))/3 + 5**(0.5D0)*Q2(26,IELEM**3+1) + (11*7*
     >   *(0.5D0)*Q2(24,IELEM**3+1))/9 + (15**(0.5D0)*Q2(19,IELEM**3+1
     >   ))/3 + (77*3**(0.5D0)*Q2(32,IELEM**3+1))/27 + (21**(0.5D0)*Q2
     >   (14,IELEM**3+1))/3 - 7**(0.5D0)*Q2(30,IELEM**3+1) + 5**(0.5D0
     >   )*Q2(33,IELEM**3+1) - 15**(0.5D0)*Q2(25,IELEM**3+1) - (11*21*
     >   *(0.5D0)*Q2(20,IELEM**3+1))/27 + (5*3**(0.5D0)*Q2(39,IELEM**3
     >   +1))/3 + 5**(0.5D0)*Q2(38,IELEM**3+1) - (5*3**(0.5D0)*Q2(42,I
     >   ELEM**3+1))/3 + (11*35**(0.5D0)*Q2(12,IELEM**3+1))/27 - (5*5*
     >   *(0.5D0)*Q2(43,IELEM**3+1))/3 - (15**(0.5D0)*Q2(34,IELEM**3+1
     >   ))/3 + 21**(0.5D0)*Q2(29,IELEM**3+1) + (35**(0.5D0)*Q2(15,IEL
     >   EM**3+1))/3 + (55*7**(0.5D0)*Q2(44,IELEM**3+1))/27 - 15**(0.5
     >   D0)*Q2(37,IELEM**3+1) - (77*5**(0.5D0)*Q2(48,IELEM**3+1))/27 
     >   + (5*7**(0.5D0)*Q2(47,IELEM**3+1))/3 - 7**(0.5D0)*Q2(49,IELEM
     >   **3+1) + (77*3**(0.5D0)*Q2(56,IELEM**3+1))/27 - 7**(0.5D0)*Q2
     >   (54,IELEM**3+1) - (7*3**(0.5D0)*Q2(62,IELEM**3+1))/3 - (77*5*
     >   *(0.5D0)*Q2(60,IELEM**3+1))/27 + (5*7**(0.5D0)*Q2(59,IELEM**3
     >   +1))/3 - (7*5**(0.5D0)*Q2(63,IELEM**3+1))/3 + (77*7**(0.5D0)*
     >   Q2(64,IELEM**3+1))/27 + (21**(0.5D0)*Q2(50,IELEM**3+1))/3 + (
     >   11*35**(0.5D0)*Q2(36,IELEM**3+1))/27 + 21**(0.5D0)*Q2(53,IELE
     >   M**3+1) - 35**(0.5D0)*Q2(45,IELEM**3+1) + (35**(0.5D0)*Q2(51,
     >   IELEM**3+1))/3 - 35**(0.5D0)*Q2(57,IELEM**3+1) - (11*105**(0.
     >   5D0)*Q2(28,IELEM**3+1))/27 - (105**(0.5D0)*Q2(31,IELEM**3+1))
     >   /3 - (11*105**(0.5D0)*Q2(40,IELEM**3+1))/27 + (105**(0.5D0)*Q
     >   2(46,IELEM**3+1))/3 - (105**(0.5D0)*Q2(55,IELEM**3+1))/3 + (1
     >   05**(0.5D0)*Q2(58,IELEM**3+1))/3) 
      CORNERQ(03,1) = (Q2(01,IELEM**3+1) - Q2(06,IELEM**3+1) - (5*Q
     >   2(11,IELEM**3+1))/3 + (77*Q2(16,IELEM**3+1))/27 - Q2(18,IELEM
     >   **3+1) + 3*Q2(21,IELEM**3+1) - (5*Q2(35,IELEM**3+1))/3 + 5*Q2
     >   (41,IELEM**3+1) + (77*Q2(52,IELEM**3+1))/27 + 7*Q2(61,IELEM**
     >   3+1) + (3**(0.5D0)*Q2(02,IELEM**3+1))/3 - 3**(0.5D0)*Q2(05,IE
     >   LEM**3+1) - (5**(0.5D0)*Q2(03,IELEM**3+1))/3 - (11*7**(0.5D0)
     >   *Q2(04,IELEM**3+1))/27 + 5**(0.5D0)*Q2(09,IELEM**3+1) - 3**(0
     >   .5D0)*Q2(17,IELEM**3+1) - 7**(0.5D0)*Q2(13,IELEM**3+1) + (15*
     >   *(0.5D0)*Q2(07,IELEM**3+1))/3 + 3**(0.5D0)*Q2(22,IELEM**3+1) 
     >   + (15**(0.5D0)*Q2(10,IELEM**3+1))/3 - 5**(0.5D0)*Q2(23,IELEM*
     >   *3+1) + (11*21**(0.5D0)*Q2(08,IELEM**3+1))/27 + (5*3**(0.5D0)
     >   *Q2(27,IELEM**3+1))/3 - 5**(0.5D0)*Q2(26,IELEM**3+1) - (11*7*
     >   *(0.5D0)*Q2(24,IELEM**3+1))/9 + (15**(0.5D0)*Q2(19,IELEM**3+1
     >   ))/3 - (77*3**(0.5D0)*Q2(32,IELEM**3+1))/27 - (21**(0.5D0)*Q2
     >   (14,IELEM**3+1))/3 + 7**(0.5D0)*Q2(30,IELEM**3+1) + 5**(0.5D0
     >   )*Q2(33,IELEM**3+1) - 15**(0.5D0)*Q2(25,IELEM**3+1) + (11*21*
     >   *(0.5D0)*Q2(20,IELEM**3+1))/27 + (5*3**(0.5D0)*Q2(39,IELEM**3
     >   +1))/3 - 5**(0.5D0)*Q2(38,IELEM**3+1) + (5*3**(0.5D0)*Q2(42,I
     >   ELEM**3+1))/3 - (11*35**(0.5D0)*Q2(12,IELEM**3+1))/27 - (5*5*
     >   *(0.5D0)*Q2(43,IELEM**3+1))/3 + (15**(0.5D0)*Q2(34,IELEM**3+1
     >   ))/3 + 21**(0.5D0)*Q2(29,IELEM**3+1) + (35**(0.5D0)*Q2(15,IEL
     >   EM**3+1))/3 - (55*7**(0.5D0)*Q2(44,IELEM**3+1))/27 - 15**(0.5
     >   D0)*Q2(37,IELEM**3+1) + (77*5**(0.5D0)*Q2(48,IELEM**3+1))/27 
     >   + (5*7**(0.5D0)*Q2(47,IELEM**3+1))/3 - 7**(0.5D0)*Q2(49,IELEM
     >   **3+1) - (77*3**(0.5D0)*Q2(56,IELEM**3+1))/27 + 7**(0.5D0)*Q2
     >   (54,IELEM**3+1) + (7*3**(0.5D0)*Q2(62,IELEM**3+1))/3 + (77*5*
     >   *(0.5D0)*Q2(60,IELEM**3+1))/27 + (5*7**(0.5D0)*Q2(59,IELEM**3
     >   +1))/3 - (7*5**(0.5D0)*Q2(63,IELEM**3+1))/3 - (77*7**(0.5D0)*
     >   Q2(64,IELEM**3+1))/27 - (21**(0.5D0)*Q2(50,IELEM**3+1))/3 - (
     >   11*35**(0.5D0)*Q2(36,IELEM**3+1))/27 + 21**(0.5D0)*Q2(53,IELE
     >   M**3+1) - 35**(0.5D0)*Q2(45,IELEM**3+1) + (35**(0.5D0)*Q2(51,
     >   IELEM**3+1))/3 - 35**(0.5D0)*Q2(57,IELEM**3+1) + (11*105**(0.
     >   5D0)*Q2(28,IELEM**3+1))/27 - (105**(0.5D0)*Q2(31,IELEM**3+1))
     >   /3 + (11*105**(0.5D0)*Q2(40,IELEM**3+1))/27 - (105**(0.5D0)*Q
     >   2(46,IELEM**3+1))/3 - (105**(0.5D0)*Q2(55,IELEM**3+1))/3 - (1
     >   05**(0.5D0)*Q2(58,IELEM**3+1))/3) 
      CORNERQ(04,1) = (Q2(01,IELEM**3+1) - 3*Q2(06,IELEM**3+1) + 5*
     >   Q2(11,IELEM**3+1) - 7*Q2(16,IELEM**3+1) - 3*Q2(18,IELEM**3+1)
     >    + 3*Q2(21,IELEM**3+1) + 5*Q2(35,IELEM**3+1) + 5*Q2(41,IELEM*
     >   *3+1) - 7*Q2(52,IELEM**3+1) + 7*Q2(61,IELEM**3+1) + 3**(0.5D0
     >   )*Q2(02,IELEM**3+1) - 3**(0.5D0)*Q2(05,IELEM**3+1) + 5**(0.5D
     >   0)*Q2(03,IELEM**3+1) + 7**(0.5D0)*Q2(04,IELEM**3+1) + 5**(0.5
     >   D0)*Q2(09,IELEM**3+1) - 3**(0.5D0)*Q2(17,IELEM**3+1) - 7**(0.
     >   5D0)*Q2(13,IELEM**3+1) - 15**(0.5D0)*Q2(07,IELEM**3+1) + 3*3*
     >   *(0.5D0)*Q2(22,IELEM**3+1) + 15**(0.5D0)*Q2(10,IELEM**3+1) + 
     >   3*5**(0.5D0)*Q2(23,IELEM**3+1) - 21**(0.5D0)*Q2(08,IELEM**3+1
     >   ) - 5*3**(0.5D0)*Q2(27,IELEM**3+1) - 3*5**(0.5D0)*Q2(26,IELEM
     >   **3+1) + 3*7**(0.5D0)*Q2(24,IELEM**3+1) - 15**(0.5D0)*Q2(19,I
     >   ELEM**3+1) + 7*3**(0.5D0)*Q2(32,IELEM**3+1) - 21**(0.5D0)*Q2(
     >   14,IELEM**3+1) + 3*7**(0.5D0)*Q2(30,IELEM**3+1) + 5**(0.5D0)*
     >   Q2(33,IELEM**3+1) - 15**(0.5D0)*Q2(25,IELEM**3+1) - 21**(0.5D
     >   0)*Q2(20,IELEM**3+1) - 5*3**(0.5D0)*Q2(39,IELEM**3+1) - 3*5**
     >   (0.5D0)*Q2(38,IELEM**3+1) + 5*3**(0.5D0)*Q2(42,IELEM**3+1) + 
     >   35**(0.5D0)*Q2(12,IELEM**3+1) + 5*5**(0.5D0)*Q2(43,IELEM**3+1
     >   ) + 15**(0.5D0)*Q2(34,IELEM**3+1) + 21**(0.5D0)*Q2(29,IELEM**
     >   3+1) - 35**(0.5D0)*Q2(15,IELEM**3+1) + 5*7**(0.5D0)*Q2(44,IEL
     >   EM**3+1) - 15**(0.5D0)*Q2(37,IELEM**3+1) - 7*5**(0.5D0)*Q2(48
     >   ,IELEM**3+1) - 5*7**(0.5D0)*Q2(47,IELEM**3+1) - 7**(0.5D0)*Q2
     >   (49,IELEM**3+1) + 7*3**(0.5D0)*Q2(56,IELEM**3+1) + 3*7**(0.5D
     >   0)*Q2(54,IELEM**3+1) + 7*3**(0.5D0)*Q2(62,IELEM**3+1) - 7*5**
     >   (0.5D0)*Q2(60,IELEM**3+1) - 5*7**(0.5D0)*Q2(59,IELEM**3+1) + 
     >   7*5**(0.5D0)*Q2(63,IELEM**3+1) + 7*7**(0.5D0)*Q2(64,IELEM**3+
     >   1) - 21**(0.5D0)*Q2(50,IELEM**3+1) + 35**(0.5D0)*Q2(36,IELEM*
     >   *3+1) + 21**(0.5D0)*Q2(53,IELEM**3+1) - 35**(0.5D0)*Q2(45,IEL
     >   EM**3+1) - 35**(0.5D0)*Q2(51,IELEM**3+1) - 35**(0.5D0)*Q2(57,
     >   IELEM**3+1) - 105**(0.5D0)*Q2(28,IELEM**3+1) + 105**(0.5D0)*Q
     >   2(31,IELEM**3+1) - 105**(0.5D0)*Q2(40,IELEM**3+1) - 105**(0.5
     >   D0)*Q2(46,IELEM**3+1) + 105**(0.5D0)*Q2(55,IELEM**3+1) - 105*
     >   *(0.5D0)*Q2(58,IELEM**3+1)) 
      CORNERQ(05,1) = (Q2(01,IELEM**3+1) + Q2(06,IELEM**3+1) - (5*Q
     >   2(11,IELEM**3+1))/3 - (77*Q2(16,IELEM**3+1))/27 + 3*Q2(18,IEL
     >   EM**3+1) + Q2(21,IELEM**3+1) + 5*Q2(35,IELEM**3+1) - (5*Q2(41
     >   ,IELEM**3+1))/3 + 7*Q2(52,IELEM**3+1) - (77*Q2(61,IELEM**3+1)
     >   )/27 - 3**(0.5D0)*Q2(02,IELEM**3+1) - (3**(0.5D0)*Q2(05,IELEM
     >   **3+1))/3 + 5**(0.5D0)*Q2(03,IELEM**3+1) - 7**(0.5D0)*Q2(04,I
     >   ELEM**3+1) - (5**(0.5D0)*Q2(09,IELEM**3+1))/3 - 3**(0.5D0)*Q2
     >   (17,IELEM**3+1) + (11*7**(0.5D0)*Q2(13,IELEM**3+1))/27 - (15*
     >   *(0.5D0)*Q2(07,IELEM**3+1))/3 - 3**(0.5D0)*Q2(22,IELEM**3+1) 
     >   + (15**(0.5D0)*Q2(10,IELEM**3+1))/3 + 5**(0.5D0)*Q2(23,IELEM*
     >   *3+1) + (21**(0.5D0)*Q2(08,IELEM**3+1))/3 + (5*3**(0.5D0)*Q2(
     >   27,IELEM**3+1))/3 - 5**(0.5D0)*Q2(26,IELEM**3+1) - 7**(0.5D0)
     >   *Q2(24,IELEM**3+1) - 15**(0.5D0)*Q2(19,IELEM**3+1) + (77*3**(
     >   0.5D0)*Q2(32,IELEM**3+1))/27 - (11*21**(0.5D0)*Q2(14,IELEM**3
     >   +1))/27 + (11*7**(0.5D0)*Q2(30,IELEM**3+1))/9 + 5**(0.5D0)*Q2
     >   (33,IELEM**3+1) + (15**(0.5D0)*Q2(25,IELEM**3+1))/3 + 21**(0.
     >   5D0)*Q2(20,IELEM**3+1) - (5*3**(0.5D0)*Q2(39,IELEM**3+1))/3 +
     >    5**(0.5D0)*Q2(38,IELEM**3+1) + (5*3**(0.5D0)*Q2(42,IELEM**3+
     >   1))/3 + (35**(0.5D0)*Q2(12,IELEM**3+1))/3 - (5*5**(0.5D0)*Q2(
     >   43,IELEM**3+1))/3 - 15**(0.5D0)*Q2(34,IELEM**3+1) - (11*21**(
     >   0.5D0)*Q2(29,IELEM**3+1))/27 + (11*35**(0.5D0)*Q2(15,IELEM**3
     >   +1))/27 + (5*7**(0.5D0)*Q2(44,IELEM**3+1))/3 - (15**(0.5D0)*Q
     >   2(37,IELEM**3+1))/3 - (77*5**(0.5D0)*Q2(48,IELEM**3+1))/27 + 
     >   (55*7**(0.5D0)*Q2(47,IELEM**3+1))/27 - 7**(0.5D0)*Q2(49,IELEM
     >   **3+1) - (7*3**(0.5D0)*Q2(56,IELEM**3+1))/3 - 7**(0.5D0)*Q2(5
     >   4,IELEM**3+1) + (77*3**(0.5D0)*Q2(62,IELEM**3+1))/27 - (7*5**
     >   (0.5D0)*Q2(60,IELEM**3+1))/3 + (5*7**(0.5D0)*Q2(59,IELEM**3+1
     >   ))/3 - (77*5**(0.5D0)*Q2(63,IELEM**3+1))/27 + (77*7**(0.5D0)*
     >   Q2(64,IELEM**3+1))/27 + 21**(0.5D0)*Q2(50,IELEM**3+1) - 35**(
     >   0.5D0)*Q2(36,IELEM**3+1) + (21**(0.5D0)*Q2(53,IELEM**3+1))/3 
     >   + (11*35**(0.5D0)*Q2(45,IELEM**3+1))/27 - 35**(0.5D0)*Q2(51,I
     >   ELEM**3+1) + (35**(0.5D0)*Q2(57,IELEM**3+1))/3 - (105**(0.5D0
     >   )*Q2(28,IELEM**3+1))/3 - (11*105**(0.5D0)*Q2(31,IELEM**3+1))/
     >   27 + (105**(0.5D0)*Q2(40,IELEM**3+1))/3 - (11*105**(0.5D0)*Q2
     >   (46,IELEM**3+1))/27 + (105**(0.5D0)*Q2(55,IELEM**3+1))/3 - (1
     >   05**(0.5D0)*Q2(58,IELEM**3+1))/3) 
      CORNERQ(06,1) = (Q2(01,IELEM**3+1) + Q2(06,IELEM**3+1)/3 + (5
     >   *Q2(11,IELEM**3+1))/9 + (847*Q2(16,IELEM**3+1))/729 + Q2(18,I
     >   ELEM**3+1) + Q2(21,IELEM**3+1) - (5*Q2(35,IELEM**3+1))/3 - (5
     >   *Q2(41,IELEM**3+1))/3 - (77*Q2(52,IELEM**3+1))/27 - (77*Q2(61
     >   ,IELEM**3+1))/27 - (3**(0.5D0)*Q2(02,IELEM**3+1))/3 - (3**(0.
     >   5D0)*Q2(05,IELEM**3+1))/3 - (5**(0.5D0)*Q2(03,IELEM**3+1))/3 
     >   + (11*7**(0.5D0)*Q2(04,IELEM**3+1))/27 - (5**(0.5D0)*Q2(09,IE
     >   LEM**3+1))/3 - 3**(0.5D0)*Q2(17,IELEM**3+1) + (11*7**(0.5D0)*
     >   Q2(13,IELEM**3+1))/27 + (15**(0.5D0)*Q2(07,IELEM**3+1))/9 - (
     >   3**(0.5D0)*Q2(22,IELEM**3+1))/3 + (15**(0.5D0)*Q2(10,IELEM**3
     >   +1))/9 - (5**(0.5D0)*Q2(23,IELEM**3+1))/3 - (11*21**(0.5D0)*Q
     >   2(08,IELEM**3+1))/81 - (5*3**(0.5D0)*Q2(27,IELEM**3+1))/9 - (
     >   5**(0.5D0)*Q2(26,IELEM**3+1))/3 + (11*7**(0.5D0)*Q2(24,IELEM*
     >   *3+1))/27 + (15**(0.5D0)*Q2(19,IELEM**3+1))/3 - (847*3**(0.5D
     >   0)*Q2(32,IELEM**3+1))/729 - (11*21**(0.5D0)*Q2(14,IELEM**3+1)
     >   )/81 + (11*7**(0.5D0)*Q2(30,IELEM**3+1))/27 + 5**(0.5D0)*Q2(3
     >   3,IELEM**3+1) + (15**(0.5D0)*Q2(25,IELEM**3+1))/3 - (11*21**(
     >   0.5D0)*Q2(20,IELEM**3+1))/27 + (5*3**(0.5D0)*Q2(39,IELEM**3+1
     >   ))/9 + (5**(0.5D0)*Q2(38,IELEM**3+1))/3 + (5*3**(0.5D0)*Q2(42
     >   ,IELEM**3+1))/9 - (11*35**(0.5D0)*Q2(12,IELEM**3+1))/81 + (5*
     >   5**(0.5D0)*Q2(43,IELEM**3+1))/9 - (15**(0.5D0)*Q2(34,IELEM**3
     >   +1))/3 - (11*21**(0.5D0)*Q2(29,IELEM**3+1))/27 - (11*35**(0.5
     >   D0)*Q2(15,IELEM**3+1))/81 - (55*7**(0.5D0)*Q2(44,IELEM**3+1))
     >   /81 - (15**(0.5D0)*Q2(37,IELEM**3+1))/3 + (847*5**(0.5D0)*Q2(
     >   48,IELEM**3+1))/729 - (55*7**(0.5D0)*Q2(47,IELEM**3+1))/81 - 
     >   7**(0.5D0)*Q2(49,IELEM**3+1) + (77*3**(0.5D0)*Q2(56,IELEM**3+
     >   1))/81 - (7**(0.5D0)*Q2(54,IELEM**3+1))/3 + (77*3**(0.5D0)*Q2
     >   (62,IELEM**3+1))/81 + (77*5**(0.5D0)*Q2(60,IELEM**3+1))/81 - 
     >   (5*7**(0.5D0)*Q2(59,IELEM**3+1))/9 + (77*5**(0.5D0)*Q2(63,IEL
     >   EM**3+1))/81 - (847*7**(0.5D0)*Q2(64,IELEM**3+1))/729 + (21**
     >   (0.5D0)*Q2(50,IELEM**3+1))/3 + (11*35**(0.5D0)*Q2(36,IELEM**3
     >   +1))/27 + (21**(0.5D0)*Q2(53,IELEM**3+1))/3 + (11*35**(0.5D0)
     >   *Q2(45,IELEM**3+1))/27 + (35**(0.5D0)*Q2(51,IELEM**3+1))/3 + 
     >   (35**(0.5D0)*Q2(57,IELEM**3+1))/3 + (11*105**(0.5D0)*Q2(28,IE
     >   LEM**3+1))/81 + (11*105**(0.5D0)*Q2(31,IELEM**3+1))/81 - (11*
     >   105**(0.5D0)*Q2(40,IELEM**3+1))/81 - (11*105**(0.5D0)*Q2(46,I
     >   ELEM**3+1))/81 - (105**(0.5D0)*Q2(55,IELEM**3+1))/9 - (105**(
     >   0.5D0)*Q2(58,IELEM**3+1))/9) 
      CORNERQ(07,1) = (Q2(01,IELEM**3+1) - Q2(06,IELEM**3+1)/3 + (5
     >   *Q2(11,IELEM**3+1))/9 - (847*Q2(16,IELEM**3+1))/729 - Q2(18,I
     >   ELEM**3+1) + Q2(21,IELEM**3+1) - (5*Q2(35,IELEM**3+1))/3 - (5
     >   *Q2(41,IELEM**3+1))/3 + (77*Q2(52,IELEM**3+1))/27 - (77*Q2(61
     >   ,IELEM**3+1))/27 + (3**(0.5D0)*Q2(02,IELEM**3+1))/3 - (3**(0.
     >   5D0)*Q2(05,IELEM**3+1))/3 - (5**(0.5D0)*Q2(03,IELEM**3+1))/3 
     >   - (11*7**(0.5D0)*Q2(04,IELEM**3+1))/27 - (5**(0.5D0)*Q2(09,IE
     >   LEM**3+1))/3 - 3**(0.5D0)*Q2(17,IELEM**3+1) + (11*7**(0.5D0)*
     >   Q2(13,IELEM**3+1))/27 + (15**(0.5D0)*Q2(07,IELEM**3+1))/9 + (
     >   3**(0.5D0)*Q2(22,IELEM**3+1))/3 - (15**(0.5D0)*Q2(10,IELEM**3
     >   +1))/9 - (5**(0.5D0)*Q2(23,IELEM**3+1))/3 + (11*21**(0.5D0)*Q
     >   2(08,IELEM**3+1))/81 - (5*3**(0.5D0)*Q2(27,IELEM**3+1))/9 + (
     >   5**(0.5D0)*Q2(26,IELEM**3+1))/3 - (11*7**(0.5D0)*Q2(24,IELEM*
     >   *3+1))/27 + (15**(0.5D0)*Q2(19,IELEM**3+1))/3 + (847*3**(0.5D
     >   0)*Q2(32,IELEM**3+1))/729 + (11*21**(0.5D0)*Q2(14,IELEM**3+1)
     >   )/81 - (11*7**(0.5D0)*Q2(30,IELEM**3+1))/27 + 5**(0.5D0)*Q2(3
     >   3,IELEM**3+1) + (15**(0.5D0)*Q2(25,IELEM**3+1))/3 + (11*21**(
     >   0.5D0)*Q2(20,IELEM**3+1))/27 + (5*3**(0.5D0)*Q2(39,IELEM**3+1
     >   ))/9 - (5**(0.5D0)*Q2(38,IELEM**3+1))/3 - (5*3**(0.5D0)*Q2(42
     >   ,IELEM**3+1))/9 + (11*35**(0.5D0)*Q2(12,IELEM**3+1))/81 + (5*
     >   5**(0.5D0)*Q2(43,IELEM**3+1))/9 + (15**(0.5D0)*Q2(34,IELEM**3
     >   +1))/3 - (11*21**(0.5D0)*Q2(29,IELEM**3+1))/27 - (11*35**(0.5
     >   D0)*Q2(15,IELEM**3+1))/81 + (55*7**(0.5D0)*Q2(44,IELEM**3+1))
     >   /81 - (15**(0.5D0)*Q2(37,IELEM**3+1))/3 - (847*5**(0.5D0)*Q2(
     >   48,IELEM**3+1))/729 - (55*7**(0.5D0)*Q2(47,IELEM**3+1))/81 - 
     >   7**(0.5D0)*Q2(49,IELEM**3+1) - (77*3**(0.5D0)*Q2(56,IELEM**3+
     >   1))/81 + (7**(0.5D0)*Q2(54,IELEM**3+1))/3 - (77*3**(0.5D0)*Q2
     >   (62,IELEM**3+1))/81 - (77*5**(0.5D0)*Q2(60,IELEM**3+1))/81 - 
     >   (5*7**(0.5D0)*Q2(59,IELEM**3+1))/9 + (77*5**(0.5D0)*Q2(63,IEL
     >   EM**3+1))/81 + (847*7**(0.5D0)*Q2(64,IELEM**3+1))/729 - (21**
     >   (0.5D0)*Q2(50,IELEM**3+1))/3 - (11*35**(0.5D0)*Q2(36,IELEM**3
     >   +1))/27 + (21**(0.5D0)*Q2(53,IELEM**3+1))/3 + (11*35**(0.5D0)
     >   *Q2(45,IELEM**3+1))/27 + (35**(0.5D0)*Q2(51,IELEM**3+1))/3 + 
     >   (35**(0.5D0)*Q2(57,IELEM**3+1))/3 - (11*105**(0.5D0)*Q2(28,IE
     >   LEM**3+1))/81 + (11*105**(0.5D0)*Q2(31,IELEM**3+1))/81 + (11*
     >   105**(0.5D0)*Q2(40,IELEM**3+1))/81 + (11*105**(0.5D0)*Q2(46,I
     >   ELEM**3+1))/81 - (105**(0.5D0)*Q2(55,IELEM**3+1))/9 + (105**(
     >   0.5D0)*Q2(58,IELEM**3+1))/9) 
      CORNERQ(08,1) = (Q2(01,IELEM**3+1) - Q2(06,IELEM**3+1) - (5*Q
     >   2(11,IELEM**3+1))/3 + (77*Q2(16,IELEM**3+1))/27 - 3*Q2(18,IEL
     >   EM**3+1) + Q2(21,IELEM**3+1) + 5*Q2(35,IELEM**3+1) - (5*Q2(41
     >   ,IELEM**3+1))/3 - 7*Q2(52,IELEM**3+1) - (77*Q2(61,IELEM**3+1)
     >   )/27 + 3**(0.5D0)*Q2(02,IELEM**3+1) - (3**(0.5D0)*Q2(05,IELEM
     >   **3+1))/3 + 5**(0.5D0)*Q2(03,IELEM**3+1) + 7**(0.5D0)*Q2(04,I
     >   ELEM**3+1) - (5**(0.5D0)*Q2(09,IELEM**3+1))/3 - 3**(0.5D0)*Q2
     >   (17,IELEM**3+1) + (11*7**(0.5D0)*Q2(13,IELEM**3+1))/27 - (15*
     >   *(0.5D0)*Q2(07,IELEM**3+1))/3 + 3**(0.5D0)*Q2(22,IELEM**3+1) 
     >   - (15**(0.5D0)*Q2(10,IELEM**3+1))/3 + 5**(0.5D0)*Q2(23,IELEM*
     >   *3+1) - (21**(0.5D0)*Q2(08,IELEM**3+1))/3 + (5*3**(0.5D0)*Q2(
     >   27,IELEM**3+1))/3 + 5**(0.5D0)*Q2(26,IELEM**3+1) + 7**(0.5D0)
     >   *Q2(24,IELEM**3+1) - 15**(0.5D0)*Q2(19,IELEM**3+1) - (77*3**(
     >   0.5D0)*Q2(32,IELEM**3+1))/27 + (11*21**(0.5D0)*Q2(14,IELEM**3
     >   +1))/27 - (11*7**(0.5D0)*Q2(30,IELEM**3+1))/9 + 5**(0.5D0)*Q2
     >   (33,IELEM**3+1) + (15**(0.5D0)*Q2(25,IELEM**3+1))/3 - 21**(0.
     >   5D0)*Q2(20,IELEM**3+1) - (5*3**(0.5D0)*Q2(39,IELEM**3+1))/3 -
     >    5**(0.5D0)*Q2(38,IELEM**3+1) - (5*3**(0.5D0)*Q2(42,IELEM**3+
     >   1))/3 - (35**(0.5D0)*Q2(12,IELEM**3+1))/3 - (5*5**(0.5D0)*Q2(
     >   43,IELEM**3+1))/3 + 15**(0.5D0)*Q2(34,IELEM**3+1) - (11*21**(
     >   0.5D0)*Q2(29,IELEM**3+1))/27 + (11*35**(0.5D0)*Q2(15,IELEM**3
     >   +1))/27 - (5*7**(0.5D0)*Q2(44,IELEM**3+1))/3 - (15**(0.5D0)*Q
     >   2(37,IELEM**3+1))/3 + (77*5**(0.5D0)*Q2(48,IELEM**3+1))/27 + 
     >   (55*7**(0.5D0)*Q2(47,IELEM**3+1))/27 - 7**(0.5D0)*Q2(49,IELEM
     >   **3+1) + (7*3**(0.5D0)*Q2(56,IELEM**3+1))/3 + 7**(0.5D0)*Q2(5
     >   4,IELEM**3+1) - (77*3**(0.5D0)*Q2(62,IELEM**3+1))/27 + (7*5**
     >   (0.5D0)*Q2(60,IELEM**3+1))/3 + (5*7**(0.5D0)*Q2(59,IELEM**3+1
     >   ))/3 - (77*5**(0.5D0)*Q2(63,IELEM**3+1))/27 - (77*7**(0.5D0)*
     >   Q2(64,IELEM**3+1))/27 - 21**(0.5D0)*Q2(50,IELEM**3+1) + 35**(
     >   0.5D0)*Q2(36,IELEM**3+1) + (21**(0.5D0)*Q2(53,IELEM**3+1))/3 
     >   + (11*35**(0.5D0)*Q2(45,IELEM**3+1))/27 - 35**(0.5D0)*Q2(51,I
     >   ELEM**3+1) + (35**(0.5D0)*Q2(57,IELEM**3+1))/3 + (105**(0.5D0
     >   )*Q2(28,IELEM**3+1))/3 - (11*105**(0.5D0)*Q2(31,IELEM**3+1))/
     >   27 - (105**(0.5D0)*Q2(40,IELEM**3+1))/3 + (11*105**(0.5D0)*Q2
     >   (46,IELEM**3+1))/27 + (105**(0.5D0)*Q2(55,IELEM**3+1))/3 + (1
     >   05**(0.5D0)*Q2(58,IELEM**3+1))/3) 
      CORNERQ(09,1) = (Q2(01,IELEM**3+1) - Q2(06,IELEM**3+1) - (5*Q
     >   2(11,IELEM**3+1))/3 + (77*Q2(16,IELEM**3+1))/27 + 3*Q2(18,IEL
     >   EM**3+1) - Q2(21,IELEM**3+1) + 5*Q2(35,IELEM**3+1) - (5*Q2(41
     >   ,IELEM**3+1))/3 + 7*Q2(52,IELEM**3+1) + (77*Q2(61,IELEM**3+1)
     >   )/27 - 3**(0.5D0)*Q2(02,IELEM**3+1) + (3**(0.5D0)*Q2(05,IELEM
     >   **3+1))/3 + 5**(0.5D0)*Q2(03,IELEM**3+1) - 7**(0.5D0)*Q2(04,I
     >   ELEM**3+1) - (5**(0.5D0)*Q2(09,IELEM**3+1))/3 - 3**(0.5D0)*Q2
     >   (17,IELEM**3+1) - (11*7**(0.5D0)*Q2(13,IELEM**3+1))/27 + (15*
     >   *(0.5D0)*Q2(07,IELEM**3+1))/3 + 3**(0.5D0)*Q2(22,IELEM**3+1) 
     >   + (15**(0.5D0)*Q2(10,IELEM**3+1))/3 - 5**(0.5D0)*Q2(23,IELEM*
     >   *3+1) - (21**(0.5D0)*Q2(08,IELEM**3+1))/3 + (5*3**(0.5D0)*Q2(
     >   27,IELEM**3+1))/3 - 5**(0.5D0)*Q2(26,IELEM**3+1) + 7**(0.5D0)
     >   *Q2(24,IELEM**3+1) - 15**(0.5D0)*Q2(19,IELEM**3+1) - (77*3**(
     >   0.5D0)*Q2(32,IELEM**3+1))/27 + (11*21**(0.5D0)*Q2(14,IELEM**3
     >   +1))/27 - (11*7**(0.5D0)*Q2(30,IELEM**3+1))/9 + 5**(0.5D0)*Q2
     >   (33,IELEM**3+1) + (15**(0.5D0)*Q2(25,IELEM**3+1))/3 + 21**(0.
     >   5D0)*Q2(20,IELEM**3+1) + (5*3**(0.5D0)*Q2(39,IELEM**3+1))/3 -
     >    5**(0.5D0)*Q2(38,IELEM**3+1) + (5*3**(0.5D0)*Q2(42,IELEM**3+
     >   1))/3 + (35**(0.5D0)*Q2(12,IELEM**3+1))/3 - (5*5**(0.5D0)*Q2(
     >   43,IELEM**3+1))/3 - 15**(0.5D0)*Q2(34,IELEM**3+1) + (11*21**(
     >   0.5D0)*Q2(29,IELEM**3+1))/27 - (11*35**(0.5D0)*Q2(15,IELEM**3
     >   +1))/27 + (5*7**(0.5D0)*Q2(44,IELEM**3+1))/3 + (15**(0.5D0)*Q
     >   2(37,IELEM**3+1))/3 + (77*5**(0.5D0)*Q2(48,IELEM**3+1))/27 - 
     >   (55*7**(0.5D0)*Q2(47,IELEM**3+1))/27 - 7**(0.5D0)*Q2(49,IELEM
     >   **3+1) + (7*3**(0.5D0)*Q2(56,IELEM**3+1))/3 + 7**(0.5D0)*Q2(5
     >   4,IELEM**3+1) - (77*3**(0.5D0)*Q2(62,IELEM**3+1))/27 - (7*5**
     >   (0.5D0)*Q2(60,IELEM**3+1))/3 + (5*7**(0.5D0)*Q2(59,IELEM**3+1
     >   ))/3 + (77*5**(0.5D0)*Q2(63,IELEM**3+1))/27 - (77*7**(0.5D0)*
     >   Q2(64,IELEM**3+1))/27 + 21**(0.5D0)*Q2(50,IELEM**3+1) - 35**(
     >   0.5D0)*Q2(36,IELEM**3+1) - (21**(0.5D0)*Q2(53,IELEM**3+1))/3 
     >   - (11*35**(0.5D0)*Q2(45,IELEM**3+1))/27 - 35**(0.5D0)*Q2(51,I
     >   ELEM**3+1) + (35**(0.5D0)*Q2(57,IELEM**3+1))/3 - (105**(0.5D0
     >   )*Q2(28,IELEM**3+1))/3 + (11*105**(0.5D0)*Q2(31,IELEM**3+1))/
     >   27 - (105**(0.5D0)*Q2(40,IELEM**3+1))/3 + (11*105**(0.5D0)*Q2
     >   (46,IELEM**3+1))/27 - (105**(0.5D0)*Q2(55,IELEM**3+1))/3 - (1
     >   05**(0.5D0)*Q2(58,IELEM**3+1))/3) 
      CORNERQ(10,1) = (Q2(01,IELEM**3+1) - Q2(06,IELEM**3+1)/3 + (5
     >   *Q2(11,IELEM**3+1))/9 - (847*Q2(16,IELEM**3+1))/729 + Q2(18,I
     >   ELEM**3+1) - Q2(21,IELEM**3+1) - (5*Q2(35,IELEM**3+1))/3 - (5
     >   *Q2(41,IELEM**3+1))/3 - (77*Q2(52,IELEM**3+1))/27 + (77*Q2(61
     >   ,IELEM**3+1))/27 - (3**(0.5D0)*Q2(02,IELEM**3+1))/3 + (3**(0.
     >   5D0)*Q2(05,IELEM**3+1))/3 - (5**(0.5D0)*Q2(03,IELEM**3+1))/3 
     >   + (11*7**(0.5D0)*Q2(04,IELEM**3+1))/27 - (5**(0.5D0)*Q2(09,IE
     >   LEM**3+1))/3 - 3**(0.5D0)*Q2(17,IELEM**3+1) - (11*7**(0.5D0)*
     >   Q2(13,IELEM**3+1))/27 - (15**(0.5D0)*Q2(07,IELEM**3+1))/9 + (
     >   3**(0.5D0)*Q2(22,IELEM**3+1))/3 + (15**(0.5D0)*Q2(10,IELEM**3
     >   +1))/9 + (5**(0.5D0)*Q2(23,IELEM**3+1))/3 + (11*21**(0.5D0)*Q
     >   2(08,IELEM**3+1))/81 - (5*3**(0.5D0)*Q2(27,IELEM**3+1))/9 - (
     >   5**(0.5D0)*Q2(26,IELEM**3+1))/3 - (11*7**(0.5D0)*Q2(24,IELEM*
     >   *3+1))/27 + (15**(0.5D0)*Q2(19,IELEM**3+1))/3 + (847*3**(0.5D
     >   0)*Q2(32,IELEM**3+1))/729 + (11*21**(0.5D0)*Q2(14,IELEM**3+1)
     >   )/81 - (11*7**(0.5D0)*Q2(30,IELEM**3+1))/27 + 5**(0.5D0)*Q2(3
     >   3,IELEM**3+1) + (15**(0.5D0)*Q2(25,IELEM**3+1))/3 - (11*21**(
     >   0.5D0)*Q2(20,IELEM**3+1))/27 - (5*3**(0.5D0)*Q2(39,IELEM**3+1
     >   ))/9 - (5**(0.5D0)*Q2(38,IELEM**3+1))/3 + (5*3**(0.5D0)*Q2(42
     >   ,IELEM**3+1))/9 - (11*35**(0.5D0)*Q2(12,IELEM**3+1))/81 + (5*
     >   5**(0.5D0)*Q2(43,IELEM**3+1))/9 - (15**(0.5D0)*Q2(34,IELEM**3
     >   +1))/3 + (11*21**(0.5D0)*Q2(29,IELEM**3+1))/27 + (11*35**(0.5
     >   D0)*Q2(15,IELEM**3+1))/81 - (55*7**(0.5D0)*Q2(44,IELEM**3+1))
     >   /81 + (15**(0.5D0)*Q2(37,IELEM**3+1))/3 - (847*5**(0.5D0)*Q2(
     >   48,IELEM**3+1))/729 + (55*7**(0.5D0)*Q2(47,IELEM**3+1))/81 - 
     >   7**(0.5D0)*Q2(49,IELEM**3+1) - (77*3**(0.5D0)*Q2(56,IELEM**3+
     >   1))/81 + (7**(0.5D0)*Q2(54,IELEM**3+1))/3 - (77*3**(0.5D0)*Q2
     >   (62,IELEM**3+1))/81 + (77*5**(0.5D0)*Q2(60,IELEM**3+1))/81 - 
     >   (5*7**(0.5D0)*Q2(59,IELEM**3+1))/9 - (77*5**(0.5D0)*Q2(63,IEL
     >   EM**3+1))/81 + (847*7**(0.5D0)*Q2(64,IELEM**3+1))/729 + (21**
     >   (0.5D0)*Q2(50,IELEM**3+1))/3 + (11*35**(0.5D0)*Q2(36,IELEM**3
     >   +1))/27 - (21**(0.5D0)*Q2(53,IELEM**3+1))/3 - (11*35**(0.5D0)
     >   *Q2(45,IELEM**3+1))/27 + (35**(0.5D0)*Q2(51,IELEM**3+1))/3 + 
     >   (35**(0.5D0)*Q2(57,IELEM**3+1))/3 + (11*105**(0.5D0)*Q2(28,IE
     >   LEM**3+1))/81 - (11*105**(0.5D0)*Q2(31,IELEM**3+1))/81 + (11*
     >   105**(0.5D0)*Q2(40,IELEM**3+1))/81 + (11*105**(0.5D0)*Q2(46,I
     >   ELEM**3+1))/81 + (105**(0.5D0)*Q2(55,IELEM**3+1))/9 - (105**(
     >   0.5D0)*Q2(58,IELEM**3+1))/9) 
      CORNERQ(11,1) = (Q2(01,IELEM**3+1) + Q2(06,IELEM**3+1)/3 + (5
     >   *Q2(11,IELEM**3+1))/9 + (847*Q2(16,IELEM**3+1))/729 - Q2(18,I
     >   ELEM**3+1) - Q2(21,IELEM**3+1) - (5*Q2(35,IELEM**3+1))/3 - (5
     >   *Q2(41,IELEM**3+1))/3 + (77*Q2(52,IELEM**3+1))/27 + (77*Q2(61
     >   ,IELEM**3+1))/27 + (3**(0.5D0)*Q2(02,IELEM**3+1))/3 + (3**(0.
     >   5D0)*Q2(05,IELEM**3+1))/3 - (5**(0.5D0)*Q2(03,IELEM**3+1))/3 
     >   - (11*7**(0.5D0)*Q2(04,IELEM**3+1))/27 - (5**(0.5D0)*Q2(09,IE
     >   LEM**3+1))/3 - 3**(0.5D0)*Q2(17,IELEM**3+1) - (11*7**(0.5D0)*
     >   Q2(13,IELEM**3+1))/27 - (15**(0.5D0)*Q2(07,IELEM**3+1))/9 - (
     >   3**(0.5D0)*Q2(22,IELEM**3+1))/3 - (15**(0.5D0)*Q2(10,IELEM**3
     >   +1))/9 + (5**(0.5D0)*Q2(23,IELEM**3+1))/3 - (11*21**(0.5D0)*Q
     >   2(08,IELEM**3+1))/81 - (5*3**(0.5D0)*Q2(27,IELEM**3+1))/9 + (
     >   5**(0.5D0)*Q2(26,IELEM**3+1))/3 + (11*7**(0.5D0)*Q2(24,IELEM*
     >   *3+1))/27 + (15**(0.5D0)*Q2(19,IELEM**3+1))/3 - (847*3**(0.5D
     >   0)*Q2(32,IELEM**3+1))/729 - (11*21**(0.5D0)*Q2(14,IELEM**3+1)
     >   )/81 + (11*7**(0.5D0)*Q2(30,IELEM**3+1))/27 + 5**(0.5D0)*Q2(3
     >   3,IELEM**3+1) + (15**(0.5D0)*Q2(25,IELEM**3+1))/3 + (11*21**(
     >   0.5D0)*Q2(20,IELEM**3+1))/27 - (5*3**(0.5D0)*Q2(39,IELEM**3+1
     >   ))/9 + (5**(0.5D0)*Q2(38,IELEM**3+1))/3 - (5*3**(0.5D0)*Q2(42
     >   ,IELEM**3+1))/9 + (11*35**(0.5D0)*Q2(12,IELEM**3+1))/81 + (5*
     >   5**(0.5D0)*Q2(43,IELEM**3+1))/9 + (15**(0.5D0)*Q2(34,IELEM**3
     >   +1))/3 + (11*21**(0.5D0)*Q2(29,IELEM**3+1))/27 + (11*35**(0.5
     >   D0)*Q2(15,IELEM**3+1))/81 + (55*7**(0.5D0)*Q2(44,IELEM**3+1))
     >   /81 + (15**(0.5D0)*Q2(37,IELEM**3+1))/3 + (847*5**(0.5D0)*Q2(
     >   48,IELEM**3+1))/729 + (55*7**(0.5D0)*Q2(47,IELEM**3+1))/81 - 
     >   7**(0.5D0)*Q2(49,IELEM**3+1) + (77*3**(0.5D0)*Q2(56,IELEM**3+
     >   1))/81 - (7**(0.5D0)*Q2(54,IELEM**3+1))/3 + (77*3**(0.5D0)*Q2
     >   (62,IELEM**3+1))/81 - (77*5**(0.5D0)*Q2(60,IELEM**3+1))/81 - 
     >   (5*7**(0.5D0)*Q2(59,IELEM**3+1))/9 - (77*5**(0.5D0)*Q2(63,IEL
     >   EM**3+1))/81 - (847*7**(0.5D0)*Q2(64,IELEM**3+1))/729 - (21**
     >   (0.5D0)*Q2(50,IELEM**3+1))/3 - (11*35**(0.5D0)*Q2(36,IELEM**3
     >   +1))/27 - (21**(0.5D0)*Q2(53,IELEM**3+1))/3 - (11*35**(0.5D0)
     >   *Q2(45,IELEM**3+1))/27 + (35**(0.5D0)*Q2(51,IELEM**3+1))/3 + 
     >   (35**(0.5D0)*Q2(57,IELEM**3+1))/3 - (11*105**(0.5D0)*Q2(28,IE
     >   LEM**3+1))/81 - (11*105**(0.5D0)*Q2(31,IELEM**3+1))/81 - (11*
     >   105**(0.5D0)*Q2(40,IELEM**3+1))/81 - (11*105**(0.5D0)*Q2(46,I
     >   ELEM**3+1))/81 + (105**(0.5D0)*Q2(55,IELEM**3+1))/9 + (105**(
     >   0.5D0)*Q2(58,IELEM**3+1))/9) 
      CORNERQ(12,1) = (Q2(01,IELEM**3+1) + Q2(06,IELEM**3+1) - (5*Q
     >   2(11,IELEM**3+1))/3 - (77*Q2(16,IELEM**3+1))/27 - 3*Q2(18,IEL
     >   EM**3+1) - Q2(21,IELEM**3+1) + 5*Q2(35,IELEM**3+1) - (5*Q2(41
     >   ,IELEM**3+1))/3 - 7*Q2(52,IELEM**3+1) + (77*Q2(61,IELEM**3+1)
     >   )/27 + 3**(0.5D0)*Q2(02,IELEM**3+1) + (3**(0.5D0)*Q2(05,IELEM
     >   **3+1))/3 + 5**(0.5D0)*Q2(03,IELEM**3+1) + 7**(0.5D0)*Q2(04,I
     >   ELEM**3+1) - (5**(0.5D0)*Q2(09,IELEM**3+1))/3 - 3**(0.5D0)*Q2
     >   (17,IELEM**3+1) - (11*7**(0.5D0)*Q2(13,IELEM**3+1))/27 + (15*
     >   *(0.5D0)*Q2(07,IELEM**3+1))/3 - 3**(0.5D0)*Q2(22,IELEM**3+1) 
     >   - (15**(0.5D0)*Q2(10,IELEM**3+1))/3 - 5**(0.5D0)*Q2(23,IELEM*
     >   *3+1) + (21**(0.5D0)*Q2(08,IELEM**3+1))/3 + (5*3**(0.5D0)*Q2(
     >   27,IELEM**3+1))/3 + 5**(0.5D0)*Q2(26,IELEM**3+1) - 7**(0.5D0)
     >   *Q2(24,IELEM**3+1) - 15**(0.5D0)*Q2(19,IELEM**3+1) + (77*3**(
     >   0.5D0)*Q2(32,IELEM**3+1))/27 - (11*21**(0.5D0)*Q2(14,IELEM**3
     >   +1))/27 + (11*7**(0.5D0)*Q2(30,IELEM**3+1))/9 + 5**(0.5D0)*Q2
     >   (33,IELEM**3+1) + (15**(0.5D0)*Q2(25,IELEM**3+1))/3 - 21**(0.
     >   5D0)*Q2(20,IELEM**3+1) + (5*3**(0.5D0)*Q2(39,IELEM**3+1))/3 +
     >    5**(0.5D0)*Q2(38,IELEM**3+1) - (5*3**(0.5D0)*Q2(42,IELEM**3+
     >   1))/3 - (35**(0.5D0)*Q2(12,IELEM**3+1))/3 - (5*5**(0.5D0)*Q2(
     >   43,IELEM**3+1))/3 + 15**(0.5D0)*Q2(34,IELEM**3+1) + (11*21**(
     >   0.5D0)*Q2(29,IELEM**3+1))/27 - (11*35**(0.5D0)*Q2(15,IELEM**3
     >   +1))/27 - (5*7**(0.5D0)*Q2(44,IELEM**3+1))/3 + (15**(0.5D0)*Q
     >   2(37,IELEM**3+1))/3 - (77*5**(0.5D0)*Q2(48,IELEM**3+1))/27 - 
     >   (55*7**(0.5D0)*Q2(47,IELEM**3+1))/27 - 7**(0.5D0)*Q2(49,IELEM
     >   **3+1) - (7*3**(0.5D0)*Q2(56,IELEM**3+1))/3 - 7**(0.5D0)*Q2(5
     >   4,IELEM**3+1) + (77*3**(0.5D0)*Q2(62,IELEM**3+1))/27 + (7*5**
     >   (0.5D0)*Q2(60,IELEM**3+1))/3 + (5*7**(0.5D0)*Q2(59,IELEM**3+1
     >   ))/3 + (77*5**(0.5D0)*Q2(63,IELEM**3+1))/27 + (77*7**(0.5D0)*
     >   Q2(64,IELEM**3+1))/27 - 21**(0.5D0)*Q2(50,IELEM**3+1) + 35**(
     >   0.5D0)*Q2(36,IELEM**3+1) - (21**(0.5D0)*Q2(53,IELEM**3+1))/3 
     >   - (11*35**(0.5D0)*Q2(45,IELEM**3+1))/27 - 35**(0.5D0)*Q2(51,I
     >   ELEM**3+1) + (35**(0.5D0)*Q2(57,IELEM**3+1))/3 + (105**(0.5D0
     >   )*Q2(28,IELEM**3+1))/3 + (11*105**(0.5D0)*Q2(31,IELEM**3+1))/
     >   27 + (105**(0.5D0)*Q2(40,IELEM**3+1))/3 - (11*105**(0.5D0)*Q2
     >   (46,IELEM**3+1))/27 - (105**(0.5D0)*Q2(55,IELEM**3+1))/3 + (1
     >   05**(0.5D0)*Q2(58,IELEM**3+1))/3) 
      CORNERQ(13,1) = (Q2(01,IELEM**3+1) - 3*Q2(06,IELEM**3+1) + 5*
     >   Q2(11,IELEM**3+1) - 7*Q2(16,IELEM**3+1) + 3*Q2(18,IELEM**3+1)
     >    - 3*Q2(21,IELEM**3+1) + 5*Q2(35,IELEM**3+1) + 5*Q2(41,IELEM*
     >   *3+1) + 7*Q2(52,IELEM**3+1) - 7*Q2(61,IELEM**3+1) - 3**(0.5D0
     >   )*Q2(02,IELEM**3+1) + 3**(0.5D0)*Q2(05,IELEM**3+1) + 5**(0.5D
     >   0)*Q2(03,IELEM**3+1) - 7**(0.5D0)*Q2(04,IELEM**3+1) + 5**(0.5
     >   D0)*Q2(09,IELEM**3+1) - 3**(0.5D0)*Q2(17,IELEM**3+1) + 7**(0.
     >   5D0)*Q2(13,IELEM**3+1) + 15**(0.5D0)*Q2(07,IELEM**3+1) + 3*3*
     >   *(0.5D0)*Q2(22,IELEM**3+1) - 15**(0.5D0)*Q2(10,IELEM**3+1) - 
     >   3*5**(0.5D0)*Q2(23,IELEM**3+1) - 21**(0.5D0)*Q2(08,IELEM**3+1
     >   ) - 5*3**(0.5D0)*Q2(27,IELEM**3+1) + 3*5**(0.5D0)*Q2(26,IELEM
     >   **3+1) + 3*7**(0.5D0)*Q2(24,IELEM**3+1) - 15**(0.5D0)*Q2(19,I
     >   ELEM**3+1) + 7*3**(0.5D0)*Q2(32,IELEM**3+1) - 21**(0.5D0)*Q2(
     >   14,IELEM**3+1) + 3*7**(0.5D0)*Q2(30,IELEM**3+1) + 5**(0.5D0)*
     >   Q2(33,IELEM**3+1) - 15**(0.5D0)*Q2(25,IELEM**3+1) + 21**(0.5D
     >   0)*Q2(20,IELEM**3+1) + 5*3**(0.5D0)*Q2(39,IELEM**3+1) - 3*5**
     >   (0.5D0)*Q2(38,IELEM**3+1) - 5*3**(0.5D0)*Q2(42,IELEM**3+1) - 
     >   35**(0.5D0)*Q2(12,IELEM**3+1) + 5*5**(0.5D0)*Q2(43,IELEM**3+1
     >   ) - 15**(0.5D0)*Q2(34,IELEM**3+1) - 21**(0.5D0)*Q2(29,IELEM**
     >   3+1) + 35**(0.5D0)*Q2(15,IELEM**3+1) - 5*7**(0.5D0)*Q2(44,IEL
     >   EM**3+1) + 15**(0.5D0)*Q2(37,IELEM**3+1) - 7*5**(0.5D0)*Q2(48
     >   ,IELEM**3+1) + 5*7**(0.5D0)*Q2(47,IELEM**3+1) - 7**(0.5D0)*Q2
     >   (49,IELEM**3+1) + 7*3**(0.5D0)*Q2(56,IELEM**3+1) + 3*7**(0.5D
     >   0)*Q2(54,IELEM**3+1) + 7*3**(0.5D0)*Q2(62,IELEM**3+1) + 7*5**
     >   (0.5D0)*Q2(60,IELEM**3+1) - 5*7**(0.5D0)*Q2(59,IELEM**3+1) - 
     >   7*5**(0.5D0)*Q2(63,IELEM**3+1) + 7*7**(0.5D0)*Q2(64,IELEM**3+
     >   1) + 21**(0.5D0)*Q2(50,IELEM**3+1) - 35**(0.5D0)*Q2(36,IELEM*
     >   *3+1) - 21**(0.5D0)*Q2(53,IELEM**3+1) + 35**(0.5D0)*Q2(45,IEL
     >   EM**3+1) - 35**(0.5D0)*Q2(51,IELEM**3+1) - 35**(0.5D0)*Q2(57,
     >   IELEM**3+1) + 105**(0.5D0)*Q2(28,IELEM**3+1) - 105**(0.5D0)*Q
     >   2(31,IELEM**3+1) - 105**(0.5D0)*Q2(40,IELEM**3+1) - 105**(0.5
     >   D0)*Q2(46,IELEM**3+1) - 105**(0.5D0)*Q2(55,IELEM**3+1) + 105*
     >   *(0.5D0)*Q2(58,IELEM**3+1)) 
      CORNERQ(14,1) = (Q2(01,IELEM**3+1) - Q2(06,IELEM**3+1) - (5*Q
     >   2(11,IELEM**3+1))/3 + (77*Q2(16,IELEM**3+1))/27 + Q2(18,IELEM
     >   **3+1) - 3*Q2(21,IELEM**3+1) - (5*Q2(35,IELEM**3+1))/3 + 5*Q2
     >   (41,IELEM**3+1) - (77*Q2(52,IELEM**3+1))/27 - 7*Q2(61,IELEM**
     >   3+1) - (3**(0.5D0)*Q2(02,IELEM**3+1))/3 + 3**(0.5D0)*Q2(05,IE
     >   LEM**3+1) - (5**(0.5D0)*Q2(03,IELEM**3+1))/3 + (11*7**(0.5D0)
     >   *Q2(04,IELEM**3+1))/27 + 5**(0.5D0)*Q2(09,IELEM**3+1) - 3**(0
     >   .5D0)*Q2(17,IELEM**3+1) + 7**(0.5D0)*Q2(13,IELEM**3+1) - (15*
     >   *(0.5D0)*Q2(07,IELEM**3+1))/3 + 3**(0.5D0)*Q2(22,IELEM**3+1) 
     >   - (15**(0.5D0)*Q2(10,IELEM**3+1))/3 + 5**(0.5D0)*Q2(23,IELEM*
     >   *3+1) + (11*21**(0.5D0)*Q2(08,IELEM**3+1))/27 + (5*3**(0.5D0)
     >   *Q2(27,IELEM**3+1))/3 + 5**(0.5D0)*Q2(26,IELEM**3+1) - (11*7*
     >   *(0.5D0)*Q2(24,IELEM**3+1))/9 + (15**(0.5D0)*Q2(19,IELEM**3+1
     >   ))/3 - (77*3**(0.5D0)*Q2(32,IELEM**3+1))/27 - (21**(0.5D0)*Q2
     >   (14,IELEM**3+1))/3 + 7**(0.5D0)*Q2(30,IELEM**3+1) + 5**(0.5D0
     >   )*Q2(33,IELEM**3+1) - 15**(0.5D0)*Q2(25,IELEM**3+1) - (11*21*
     >   *(0.5D0)*Q2(20,IELEM**3+1))/27 - (5*3**(0.5D0)*Q2(39,IELEM**3
     >   +1))/3 - 5**(0.5D0)*Q2(38,IELEM**3+1) - (5*3**(0.5D0)*Q2(42,I
     >   ELEM**3+1))/3 + (11*35**(0.5D0)*Q2(12,IELEM**3+1))/27 - (5*5*
     >   *(0.5D0)*Q2(43,IELEM**3+1))/3 - (15**(0.5D0)*Q2(34,IELEM**3+1
     >   ))/3 - 21**(0.5D0)*Q2(29,IELEM**3+1) - (35**(0.5D0)*Q2(15,IEL
     >   EM**3+1))/3 + (55*7**(0.5D0)*Q2(44,IELEM**3+1))/27 + 15**(0.5
     >   D0)*Q2(37,IELEM**3+1) + (77*5**(0.5D0)*Q2(48,IELEM**3+1))/27 
     >   - (5*7**(0.5D0)*Q2(47,IELEM**3+1))/3 - 7**(0.5D0)*Q2(49,IELEM
     >   **3+1) - (77*3**(0.5D0)*Q2(56,IELEM**3+1))/27 + 7**(0.5D0)*Q2
     >   (54,IELEM**3+1) + (7*3**(0.5D0)*Q2(62,IELEM**3+1))/3 - (77*5*
     >   *(0.5D0)*Q2(60,IELEM**3+1))/27 + (5*7**(0.5D0)*Q2(59,IELEM**3
     >   +1))/3 + (7*5**(0.5D0)*Q2(63,IELEM**3+1))/3 - (77*7**(0.5D0)*
     >   Q2(64,IELEM**3+1))/27 + (21**(0.5D0)*Q2(50,IELEM**3+1))/3 + (
     >   11*35**(0.5D0)*Q2(36,IELEM**3+1))/27 - 21**(0.5D0)*Q2(53,IELE
     >   M**3+1) + 35**(0.5D0)*Q2(45,IELEM**3+1) + (35**(0.5D0)*Q2(51,
     >   IELEM**3+1))/3 - 35**(0.5D0)*Q2(57,IELEM**3+1) - (11*105**(0.
     >   5D0)*Q2(28,IELEM**3+1))/27 + (105**(0.5D0)*Q2(31,IELEM**3+1))
     >   /3 + (11*105**(0.5D0)*Q2(40,IELEM**3+1))/27 - (105**(0.5D0)*Q
     >   2(46,IELEM**3+1))/3 + (105**(0.5D0)*Q2(55,IELEM**3+1))/3 + (1
     >   05**(0.5D0)*Q2(58,IELEM**3+1))/3) 
      CORNERQ(15,1) = (Q2(01,IELEM**3+1) + Q2(06,IELEM**3+1) - (5*Q
     >   2(11,IELEM**3+1))/3 - (77*Q2(16,IELEM**3+1))/27 - Q2(18,IELEM
     >   **3+1) - 3*Q2(21,IELEM**3+1) - (5*Q2(35,IELEM**3+1))/3 + 5*Q2
     >   (41,IELEM**3+1) + (77*Q2(52,IELEM**3+1))/27 - 7*Q2(61,IELEM**
     >   3+1) + (3**(0.5D0)*Q2(02,IELEM**3+1))/3 + 3**(0.5D0)*Q2(05,IE
     >   LEM**3+1) - (5**(0.5D0)*Q2(03,IELEM**3+1))/3 - (11*7**(0.5D0)
     >   *Q2(04,IELEM**3+1))/27 + 5**(0.5D0)*Q2(09,IELEM**3+1) - 3**(0
     >   .5D0)*Q2(17,IELEM**3+1) + 7**(0.5D0)*Q2(13,IELEM**3+1) - (15*
     >   *(0.5D0)*Q2(07,IELEM**3+1))/3 - 3**(0.5D0)*Q2(22,IELEM**3+1) 
     >   + (15**(0.5D0)*Q2(10,IELEM**3+1))/3 + 5**(0.5D0)*Q2(23,IELEM*
     >   *3+1) - (11*21**(0.5D0)*Q2(08,IELEM**3+1))/27 + (5*3**(0.5D0)
     >   *Q2(27,IELEM**3+1))/3 - 5**(0.5D0)*Q2(26,IELEM**3+1) + (11*7*
     >   *(0.5D0)*Q2(24,IELEM**3+1))/9 + (15**(0.5D0)*Q2(19,IELEM**3+1
     >   ))/3 + (77*3**(0.5D0)*Q2(32,IELEM**3+1))/27 + (21**(0.5D0)*Q2
     >   (14,IELEM**3+1))/3 - 7**(0.5D0)*Q2(30,IELEM**3+1) + 5**(0.5D0
     >   )*Q2(33,IELEM**3+1) - 15**(0.5D0)*Q2(25,IELEM**3+1) + (11*21*
     >   *(0.5D0)*Q2(20,IELEM**3+1))/27 - (5*3**(0.5D0)*Q2(39,IELEM**3
     >   +1))/3 + 5**(0.5D0)*Q2(38,IELEM**3+1) + (5*3**(0.5D0)*Q2(42,I
     >   ELEM**3+1))/3 - (11*35**(0.5D0)*Q2(12,IELEM**3+1))/27 - (5*5*
     >   *(0.5D0)*Q2(43,IELEM**3+1))/3 + (15**(0.5D0)*Q2(34,IELEM**3+1
     >   ))/3 - 21**(0.5D0)*Q2(29,IELEM**3+1) - (35**(0.5D0)*Q2(15,IEL
     >   EM**3+1))/3 - (55*7**(0.5D0)*Q2(44,IELEM**3+1))/27 + 15**(0.5
     >   D0)*Q2(37,IELEM**3+1) - (77*5**(0.5D0)*Q2(48,IELEM**3+1))/27 
     >   - (5*7**(0.5D0)*Q2(47,IELEM**3+1))/3 - 7**(0.5D0)*Q2(49,IELEM
     >   **3+1) + (77*3**(0.5D0)*Q2(56,IELEM**3+1))/27 - 7**(0.5D0)*Q2
     >   (54,IELEM**3+1) - (7*3**(0.5D0)*Q2(62,IELEM**3+1))/3 + (77*5*
     >   *(0.5D0)*Q2(60,IELEM**3+1))/27 + (5*7**(0.5D0)*Q2(59,IELEM**3
     >   +1))/3 + (7*5**(0.5D0)*Q2(63,IELEM**3+1))/3 + (77*7**(0.5D0)*
     >   Q2(64,IELEM**3+1))/27 - (21**(0.5D0)*Q2(50,IELEM**3+1))/3 - (
     >   11*35**(0.5D0)*Q2(36,IELEM**3+1))/27 - 21**(0.5D0)*Q2(53,IELE
     >   M**3+1) + 35**(0.5D0)*Q2(45,IELEM**3+1) + (35**(0.5D0)*Q2(51,
     >   IELEM**3+1))/3 - 35**(0.5D0)*Q2(57,IELEM**3+1) + (11*105**(0.
     >   5D0)*Q2(28,IELEM**3+1))/27 + (105**(0.5D0)*Q2(31,IELEM**3+1))
     >   /3 - (11*105**(0.5D0)*Q2(40,IELEM**3+1))/27 + (105**(0.5D0)*Q
     >   2(46,IELEM**3+1))/3 + (105**(0.5D0)*Q2(55,IELEM**3+1))/3 - (1
     >   05**(0.5D0)*Q2(58,IELEM**3+1))/3) 
      CORNERQ(16,1) = (Q2(01,IELEM**3+1) + 3*Q2(06,IELEM**3+1) + 5*
     >   Q2(11,IELEM**3+1) + 7*Q2(16,IELEM**3+1) - 3*Q2(18,IELEM**3+1)
     >    - 3*Q2(21,IELEM**3+1) + 5*Q2(35,IELEM**3+1) + 5*Q2(41,IELEM*
     >   *3+1) - 7*Q2(52,IELEM**3+1) - 7*Q2(61,IELEM**3+1) + 3**(0.5D0
     >   )*Q2(02,IELEM**3+1) + 3**(0.5D0)*Q2(05,IELEM**3+1) + 5**(0.5D
     >   0)*Q2(03,IELEM**3+1) + 7**(0.5D0)*Q2(04,IELEM**3+1) + 5**(0.5
     >   D0)*Q2(09,IELEM**3+1) - 3**(0.5D0)*Q2(17,IELEM**3+1) + 7**(0.
     >   5D0)*Q2(13,IELEM**3+1) + 15**(0.5D0)*Q2(07,IELEM**3+1) - 3*3*
     >   *(0.5D0)*Q2(22,IELEM**3+1) + 15**(0.5D0)*Q2(10,IELEM**3+1) - 
     >   3*5**(0.5D0)*Q2(23,IELEM**3+1) + 21**(0.5D0)*Q2(08,IELEM**3+1
     >   ) - 5*3**(0.5D0)*Q2(27,IELEM**3+1) - 3*5**(0.5D0)*Q2(26,IELEM
     >   **3+1) - 3*7**(0.5D0)*Q2(24,IELEM**3+1) - 15**(0.5D0)*Q2(19,I
     >   ELEM**3+1) - 7*3**(0.5D0)*Q2(32,IELEM**3+1) + 21**(0.5D0)*Q2(
     >   14,IELEM**3+1) - 3*7**(0.5D0)*Q2(30,IELEM**3+1) + 5**(0.5D0)*
     >   Q2(33,IELEM**3+1) - 15**(0.5D0)*Q2(25,IELEM**3+1) - 21**(0.5D
     >   0)*Q2(20,IELEM**3+1) + 5*3**(0.5D0)*Q2(39,IELEM**3+1) + 3*5**
     >   (0.5D0)*Q2(38,IELEM**3+1) + 5*3**(0.5D0)*Q2(42,IELEM**3+1) + 
     >   35**(0.5D0)*Q2(12,IELEM**3+1) + 5*5**(0.5D0)*Q2(43,IELEM**3+1
     >   ) + 15**(0.5D0)*Q2(34,IELEM**3+1) - 21**(0.5D0)*Q2(29,IELEM**
     >   3+1) + 35**(0.5D0)*Q2(15,IELEM**3+1) + 5*7**(0.5D0)*Q2(44,IEL
     >   EM**3+1) + 15**(0.5D0)*Q2(37,IELEM**3+1) + 7*5**(0.5D0)*Q2(48
     >   ,IELEM**3+1) + 5*7**(0.5D0)*Q2(47,IELEM**3+1) - 7**(0.5D0)*Q2
     >   (49,IELEM**3+1) - 7*3**(0.5D0)*Q2(56,IELEM**3+1) - 3*7**(0.5D
     >   0)*Q2(54,IELEM**3+1) - 7*3**(0.5D0)*Q2(62,IELEM**3+1) - 7*5**
     >   (0.5D0)*Q2(60,IELEM**3+1) - 5*7**(0.5D0)*Q2(59,IELEM**3+1) - 
     >   7*5**(0.5D0)*Q2(63,IELEM**3+1) - 7*7**(0.5D0)*Q2(64,IELEM**3+
     >   1) - 21**(0.5D0)*Q2(50,IELEM**3+1) + 35**(0.5D0)*Q2(36,IELEM*
     >   *3+1) - 21**(0.5D0)*Q2(53,IELEM**3+1) + 35**(0.5D0)*Q2(45,IEL
     >   EM**3+1) - 35**(0.5D0)*Q2(51,IELEM**3+1) - 35**(0.5D0)*Q2(57,
     >   IELEM**3+1) - 105**(0.5D0)*Q2(28,IELEM**3+1) - 105**(0.5D0)*Q
     >   2(31,IELEM**3+1) + 105**(0.5D0)*Q2(40,IELEM**3+1) + 105**(0.5
     >   D0)*Q2(46,IELEM**3+1) - 105**(0.5D0)*Q2(55,IELEM**3+1) - 105*
     >   *(0.5D0)*Q2(58,IELEM**3+1)) 
      CORNERQ(01,2) = (Q2(01,IELEM**3+1) + 3*Q2(06,IELEM**3+1) + 5*
     >   Q2(11,IELEM**3+1) + 7*Q2(16,IELEM**3+1) + Q2(18,IELEM**3+1) +
     >    Q2(21,IELEM**3+1) - (5*Q2(35,IELEM**3+1))/3 - (5*Q2(41,IELEM
     >   **3+1))/3 - (77*Q2(52,IELEM**3+1))/27 - (77*Q2(61,IELEM**3+1)
     >   )/27 - 3**(0.5D0)*Q2(02,IELEM**3+1) - 3**(0.5D0)*Q2(05,IELEM*
     >   *3+1) + 5**(0.5D0)*Q2(03,IELEM**3+1) - 7**(0.5D0)*Q2(04,IELEM
     >   **3+1) + 5**(0.5D0)*Q2(09,IELEM**3+1) - (3**(0.5D0)*Q2(17,IEL
     >   EM**3+1))/3 - 7**(0.5D0)*Q2(13,IELEM**3+1) - 15**(0.5D0)*Q2(0
     >   7,IELEM**3+1) - 3**(0.5D0)*Q2(22,IELEM**3+1) - 15**(0.5D0)*Q2
     >   (10,IELEM**3+1) + 5**(0.5D0)*Q2(23,IELEM**3+1) + 21**(0.5D0)*
     >   Q2(08,IELEM**3+1) - (5*3**(0.5D0)*Q2(27,IELEM**3+1))/3 + 5**(
     >   0.5D0)*Q2(26,IELEM**3+1) - 7**(0.5D0)*Q2(24,IELEM**3+1) - (15
     >   **(0.5D0)*Q2(19,IELEM**3+1))/3 - (7*3**(0.5D0)*Q2(32,IELEM**3
     >   +1))/3 + 21**(0.5D0)*Q2(14,IELEM**3+1) - 7**(0.5D0)*Q2(30,IEL
     >   EM**3+1) - (5**(0.5D0)*Q2(33,IELEM**3+1))/3 - (15**(0.5D0)*Q2
     >   (25,IELEM**3+1))/3 + (21**(0.5D0)*Q2(20,IELEM**3+1))/3 + (5*3
     >   **(0.5D0)*Q2(39,IELEM**3+1))/3 - 5**(0.5D0)*Q2(38,IELEM**3+1)
     >    + (5*3**(0.5D0)*Q2(42,IELEM**3+1))/3 - 35**(0.5D0)*Q2(12,IEL
     >   EM**3+1) - (5*5**(0.5D0)*Q2(43,IELEM**3+1))/3 + (15**(0.5D0)*
     >   Q2(34,IELEM**3+1))/3 + (21**(0.5D0)*Q2(29,IELEM**3+1))/3 - 35
     >   **(0.5D0)*Q2(15,IELEM**3+1) + (5*7**(0.5D0)*Q2(44,IELEM**3+1)
     >   )/3 + (15**(0.5D0)*Q2(37,IELEM**3+1))/3 - (7*5**(0.5D0)*Q2(48
     >   ,IELEM**3+1))/3 + (5*7**(0.5D0)*Q2(47,IELEM**3+1))/3 + (11*7*
     >   *(0.5D0)*Q2(49,IELEM**3+1))/27 + (77*3**(0.5D0)*Q2(56,IELEM**
     >   3+1))/27 + (11*7**(0.5D0)*Q2(54,IELEM**3+1))/9 + (77*3**(0.5D
     >   0)*Q2(62,IELEM**3+1))/27 - (77*5**(0.5D0)*Q2(60,IELEM**3+1))/
     >   27 + (55*7**(0.5D0)*Q2(59,IELEM**3+1))/27 - (77*5**(0.5D0)*Q2
     >   (63,IELEM**3+1))/27 + (77*7**(0.5D0)*Q2(64,IELEM**3+1))/27 - 
     >   (11*21**(0.5D0)*Q2(50,IELEM**3+1))/27 + (35**(0.5D0)*Q2(36,IE
     >   LEM**3+1))/3 - (11*21**(0.5D0)*Q2(53,IELEM**3+1))/27 + (35**(
     >   0.5D0)*Q2(45,IELEM**3+1))/3 + (11*35**(0.5D0)*Q2(51,IELEM**3+
     >   1))/27 + (11*35**(0.5D0)*Q2(57,IELEM**3+1))/27 + (105**(0.5D0
     >   )*Q2(28,IELEM**3+1))/3 + (105**(0.5D0)*Q2(31,IELEM**3+1))/3 -
     >    (105**(0.5D0)*Q2(40,IELEM**3+1))/3 - (105**(0.5D0)*Q2(46,IEL
     >   EM**3+1))/3 - (11*105**(0.5D0)*Q2(55,IELEM**3+1))/27 - (11*10
     >   5**(0.5D0)*Q2(58,IELEM**3+1))/27) 
      CORNERQ(02,2) = (Q2(01,IELEM**3+1) + Q2(06,IELEM**3+1) - (5*Q
     >   2(11,IELEM**3+1))/3 - (77*Q2(16,IELEM**3+1))/27 + Q2(18,IELEM
     >   **3+1)/3 + Q2(21,IELEM**3+1) + (5*Q2(35,IELEM**3+1))/9 - (5*Q
     >   2(41,IELEM**3+1))/3 + (847*Q2(52,IELEM**3+1))/729 - (77*Q2(61
     >   ,IELEM**3+1))/27 - (3**(0.5D0)*Q2(02,IELEM**3+1))/3 - 3**(0.5
     >   D0)*Q2(05,IELEM**3+1) - (5**(0.5D0)*Q2(03,IELEM**3+1))/3 + (1
     >   1*7**(0.5D0)*Q2(04,IELEM**3+1))/27 + 5**(0.5D0)*Q2(09,IELEM**
     >   3+1) - (3**(0.5D0)*Q2(17,IELEM**3+1))/3 - 7**(0.5D0)*Q2(13,IE
     >   LEM**3+1) + (15**(0.5D0)*Q2(07,IELEM**3+1))/3 - (3**(0.5D0)*Q
     >   2(22,IELEM**3+1))/3 - (15**(0.5D0)*Q2(10,IELEM**3+1))/3 - (5*
     >   *(0.5D0)*Q2(23,IELEM**3+1))/3 - (11*21**(0.5D0)*Q2(08,IELEM**
     >   3+1))/27 + (5*3**(0.5D0)*Q2(27,IELEM**3+1))/9 + (5**(0.5D0)*Q
     >   2(26,IELEM**3+1))/3 + (11*7**(0.5D0)*Q2(24,IELEM**3+1))/27 + 
     >   (15**(0.5D0)*Q2(19,IELEM**3+1))/9 + (77*3**(0.5D0)*Q2(32,IELE
     >   M**3+1))/81 + (21**(0.5D0)*Q2(14,IELEM**3+1))/3 - (7**(0.5D0)
     >   *Q2(30,IELEM**3+1))/3 - (5**(0.5D0)*Q2(33,IELEM**3+1))/3 - (1
     >   5**(0.5D0)*Q2(25,IELEM**3+1))/3 - (11*21**(0.5D0)*Q2(20,IELEM
     >   **3+1))/81 - (5*3**(0.5D0)*Q2(39,IELEM**3+1))/9 - (5**(0.5D0)
     >   *Q2(38,IELEM**3+1))/3 + (5*3**(0.5D0)*Q2(42,IELEM**3+1))/9 + 
     >   (11*35**(0.5D0)*Q2(12,IELEM**3+1))/27 + (5*5**(0.5D0)*Q2(43,I
     >   ELEM**3+1))/9 + (15**(0.5D0)*Q2(34,IELEM**3+1))/9 + (21**(0.5
     >   D0)*Q2(29,IELEM**3+1))/3 + (35**(0.5D0)*Q2(15,IELEM**3+1))/3 
     >   - (55*7**(0.5D0)*Q2(44,IELEM**3+1))/81 + (15**(0.5D0)*Q2(37,I
     >   ELEM**3+1))/3 + (77*5**(0.5D0)*Q2(48,IELEM**3+1))/81 - (5*7**
     >   (0.5D0)*Q2(47,IELEM**3+1))/9 + (11*7**(0.5D0)*Q2(49,IELEM**3+
     >   1))/27 - (847*3**(0.5D0)*Q2(56,IELEM**3+1))/729 + (11*7**(0.5
     >   D0)*Q2(54,IELEM**3+1))/27 + (77*3**(0.5D0)*Q2(62,IELEM**3+1))
     >   /81 + (847*5**(0.5D0)*Q2(60,IELEM**3+1))/729 - (55*7**(0.5D0)
     >   *Q2(59,IELEM**3+1))/81 + (77*5**(0.5D0)*Q2(63,IELEM**3+1))/81
     >    - (847*7**(0.5D0)*Q2(64,IELEM**3+1))/729 - (11*21**(0.5D0)*Q
     >   2(50,IELEM**3+1))/81 - (11*35**(0.5D0)*Q2(36,IELEM**3+1))/81 
     >   - (11*21**(0.5D0)*Q2(53,IELEM**3+1))/27 + (35**(0.5D0)*Q2(45,
     >   IELEM**3+1))/3 - (11*35**(0.5D0)*Q2(51,IELEM**3+1))/81 + (11*
     >   35**(0.5D0)*Q2(57,IELEM**3+1))/27 - (11*105**(0.5D0)*Q2(28,IE
     >   LEM**3+1))/81 - (105**(0.5D0)*Q2(31,IELEM**3+1))/9 + (11*105*
     >   *(0.5D0)*Q2(40,IELEM**3+1))/81 - (105**(0.5D0)*Q2(46,IELEM**3
     >   +1))/9 + (11*105**(0.5D0)*Q2(55,IELEM**3+1))/81 - (11*105**(0
     >   .5D0)*Q2(58,IELEM**3+1))/81) 
      CORNERQ(03,2) = (Q2(01,IELEM**3+1) - Q2(06,IELEM**3+1) - (5*Q
     >   2(11,IELEM**3+1))/3 + (77*Q2(16,IELEM**3+1))/27 - Q2(18,IELEM
     >   **3+1)/3 + Q2(21,IELEM**3+1) + (5*Q2(35,IELEM**3+1))/9 - (5*Q
     >   2(41,IELEM**3+1))/3 - (847*Q2(52,IELEM**3+1))/729 - (77*Q2(61
     >   ,IELEM**3+1))/27 + (3**(0.5D0)*Q2(02,IELEM**3+1))/3 - 3**(0.5
     >   D0)*Q2(05,IELEM**3+1) - (5**(0.5D0)*Q2(03,IELEM**3+1))/3 - (1
     >   1*7**(0.5D0)*Q2(04,IELEM**3+1))/27 + 5**(0.5D0)*Q2(09,IELEM**
     >   3+1) - (3**(0.5D0)*Q2(17,IELEM**3+1))/3 - 7**(0.5D0)*Q2(13,IE
     >   LEM**3+1) + (15**(0.5D0)*Q2(07,IELEM**3+1))/3 + (3**(0.5D0)*Q
     >   2(22,IELEM**3+1))/3 + (15**(0.5D0)*Q2(10,IELEM**3+1))/3 - (5*
     >   *(0.5D0)*Q2(23,IELEM**3+1))/3 + (11*21**(0.5D0)*Q2(08,IELEM**
     >   3+1))/27 + (5*3**(0.5D0)*Q2(27,IELEM**3+1))/9 - (5**(0.5D0)*Q
     >   2(26,IELEM**3+1))/3 - (11*7**(0.5D0)*Q2(24,IELEM**3+1))/27 + 
     >   (15**(0.5D0)*Q2(19,IELEM**3+1))/9 - (77*3**(0.5D0)*Q2(32,IELE
     >   M**3+1))/81 - (21**(0.5D0)*Q2(14,IELEM**3+1))/3 + (7**(0.5D0)
     >   *Q2(30,IELEM**3+1))/3 - (5**(0.5D0)*Q2(33,IELEM**3+1))/3 - (1
     >   5**(0.5D0)*Q2(25,IELEM**3+1))/3 + (11*21**(0.5D0)*Q2(20,IELEM
     >   **3+1))/81 - (5*3**(0.5D0)*Q2(39,IELEM**3+1))/9 + (5**(0.5D0)
     >   *Q2(38,IELEM**3+1))/3 - (5*3**(0.5D0)*Q2(42,IELEM**3+1))/9 - 
     >   (11*35**(0.5D0)*Q2(12,IELEM**3+1))/27 + (5*5**(0.5D0)*Q2(43,I
     >   ELEM**3+1))/9 - (15**(0.5D0)*Q2(34,IELEM**3+1))/9 + (21**(0.5
     >   D0)*Q2(29,IELEM**3+1))/3 + (35**(0.5D0)*Q2(15,IELEM**3+1))/3 
     >   + (55*7**(0.5D0)*Q2(44,IELEM**3+1))/81 + (15**(0.5D0)*Q2(37,I
     >   ELEM**3+1))/3 - (77*5**(0.5D0)*Q2(48,IELEM**3+1))/81 - (5*7**
     >   (0.5D0)*Q2(47,IELEM**3+1))/9 + (11*7**(0.5D0)*Q2(49,IELEM**3+
     >   1))/27 + (847*3**(0.5D0)*Q2(56,IELEM**3+1))/729 - (11*7**(0.5
     >   D0)*Q2(54,IELEM**3+1))/27 - (77*3**(0.5D0)*Q2(62,IELEM**3+1))
     >   /81 - (847*5**(0.5D0)*Q2(60,IELEM**3+1))/729 - (55*7**(0.5D0)
     >   *Q2(59,IELEM**3+1))/81 + (77*5**(0.5D0)*Q2(63,IELEM**3+1))/81
     >    + (847*7**(0.5D0)*Q2(64,IELEM**3+1))/729 + (11*21**(0.5D0)*Q
     >   2(50,IELEM**3+1))/81 + (11*35**(0.5D0)*Q2(36,IELEM**3+1))/81 
     >   - (11*21**(0.5D0)*Q2(53,IELEM**3+1))/27 + (35**(0.5D0)*Q2(45,
     >   IELEM**3+1))/3 - (11*35**(0.5D0)*Q2(51,IELEM**3+1))/81 + (11*
     >   35**(0.5D0)*Q2(57,IELEM**3+1))/27 + (11*105**(0.5D0)*Q2(28,IE
     >   LEM**3+1))/81 - (105**(0.5D0)*Q2(31,IELEM**3+1))/9 - (11*105*
     >   *(0.5D0)*Q2(40,IELEM**3+1))/81 + (105**(0.5D0)*Q2(46,IELEM**3
     >   +1))/9 + (11*105**(0.5D0)*Q2(55,IELEM**3+1))/81 + (11*105**(0
     >   .5D0)*Q2(58,IELEM**3+1))/81) 
      CORNERQ(04,2) = (Q2(01,IELEM**3+1) - 3*Q2(06,IELEM**3+1) + 5*
     >   Q2(11,IELEM**3+1) - 7*Q2(16,IELEM**3+1) - Q2(18,IELEM**3+1) +
     >    Q2(21,IELEM**3+1) - (5*Q2(35,IELEM**3+1))/3 - (5*Q2(41,IELEM
     >   **3+1))/3 + (77*Q2(52,IELEM**3+1))/27 - (77*Q2(61,IELEM**3+1)
     >   )/27 + 3**(0.5D0)*Q2(02,IELEM**3+1) - 3**(0.5D0)*Q2(05,IELEM*
     >   *3+1) + 5**(0.5D0)*Q2(03,IELEM**3+1) + 7**(0.5D0)*Q2(04,IELEM
     >   **3+1) + 5**(0.5D0)*Q2(09,IELEM**3+1) - (3**(0.5D0)*Q2(17,IEL
     >   EM**3+1))/3 - 7**(0.5D0)*Q2(13,IELEM**3+1) - 15**(0.5D0)*Q2(0
     >   7,IELEM**3+1) + 3**(0.5D0)*Q2(22,IELEM**3+1) + 15**(0.5D0)*Q2
     >   (10,IELEM**3+1) + 5**(0.5D0)*Q2(23,IELEM**3+1) - 21**(0.5D0)*
     >   Q2(08,IELEM**3+1) - (5*3**(0.5D0)*Q2(27,IELEM**3+1))/3 - 5**(
     >   0.5D0)*Q2(26,IELEM**3+1) + 7**(0.5D0)*Q2(24,IELEM**3+1) - (15
     >   **(0.5D0)*Q2(19,IELEM**3+1))/3 + (7*3**(0.5D0)*Q2(32,IELEM**3
     >   +1))/3 - 21**(0.5D0)*Q2(14,IELEM**3+1) + 7**(0.5D0)*Q2(30,IEL
     >   EM**3+1) - (5**(0.5D0)*Q2(33,IELEM**3+1))/3 - (15**(0.5D0)*Q2
     >   (25,IELEM**3+1))/3 - (21**(0.5D0)*Q2(20,IELEM**3+1))/3 + (5*3
     >   **(0.5D0)*Q2(39,IELEM**3+1))/3 + 5**(0.5D0)*Q2(38,IELEM**3+1)
     >    - (5*3**(0.5D0)*Q2(42,IELEM**3+1))/3 + 35**(0.5D0)*Q2(12,IEL
     >   EM**3+1) - (5*5**(0.5D0)*Q2(43,IELEM**3+1))/3 - (15**(0.5D0)*
     >   Q2(34,IELEM**3+1))/3 + (21**(0.5D0)*Q2(29,IELEM**3+1))/3 - 35
     >   **(0.5D0)*Q2(15,IELEM**3+1) - (5*7**(0.5D0)*Q2(44,IELEM**3+1)
     >   )/3 + (15**(0.5D0)*Q2(37,IELEM**3+1))/3 + (7*5**(0.5D0)*Q2(48
     >   ,IELEM**3+1))/3 + (5*7**(0.5D0)*Q2(47,IELEM**3+1))/3 + (11*7*
     >   *(0.5D0)*Q2(49,IELEM**3+1))/27 - (77*3**(0.5D0)*Q2(56,IELEM**
     >   3+1))/27 - (11*7**(0.5D0)*Q2(54,IELEM**3+1))/9 - (77*3**(0.5D
     >   0)*Q2(62,IELEM**3+1))/27 + (77*5**(0.5D0)*Q2(60,IELEM**3+1))/
     >   27 + (55*7**(0.5D0)*Q2(59,IELEM**3+1))/27 - (77*5**(0.5D0)*Q2
     >   (63,IELEM**3+1))/27 - (77*7**(0.5D0)*Q2(64,IELEM**3+1))/27 + 
     >   (11*21**(0.5D0)*Q2(50,IELEM**3+1))/27 - (35**(0.5D0)*Q2(36,IE
     >   LEM**3+1))/3 - (11*21**(0.5D0)*Q2(53,IELEM**3+1))/27 + (35**(
     >   0.5D0)*Q2(45,IELEM**3+1))/3 + (11*35**(0.5D0)*Q2(51,IELEM**3+
     >   1))/27 + (11*35**(0.5D0)*Q2(57,IELEM**3+1))/27 - (105**(0.5D0
     >   )*Q2(28,IELEM**3+1))/3 + (105**(0.5D0)*Q2(31,IELEM**3+1))/3 +
     >    (105**(0.5D0)*Q2(40,IELEM**3+1))/3 + (105**(0.5D0)*Q2(46,IEL
     >   EM**3+1))/3 - (11*105**(0.5D0)*Q2(55,IELEM**3+1))/27 + (11*10
     >   5**(0.5D0)*Q2(58,IELEM**3+1))/27) 
      CORNERQ(05,2) = (Q2(01,IELEM**3+1) + Q2(06,IELEM**3+1) - (5*Q
     >   2(11,IELEM**3+1))/3 - (77*Q2(16,IELEM**3+1))/27 + Q2(18,IELEM
     >   **3+1) + Q2(21,IELEM**3+1)/3 - (5*Q2(35,IELEM**3+1))/3 + (5*Q
     >   2(41,IELEM**3+1))/9 - (77*Q2(52,IELEM**3+1))/27 + (847*Q2(61,
     >   IELEM**3+1))/729 - 3**(0.5D0)*Q2(02,IELEM**3+1) - (3**(0.5D0)
     >   *Q2(05,IELEM**3+1))/3 + 5**(0.5D0)*Q2(03,IELEM**3+1) - 7**(0.
     >   5D0)*Q2(04,IELEM**3+1) - (5**(0.5D0)*Q2(09,IELEM**3+1))/3 - (
     >   3**(0.5D0)*Q2(17,IELEM**3+1))/3 + (11*7**(0.5D0)*Q2(13,IELEM*
     >   *3+1))/27 - (15**(0.5D0)*Q2(07,IELEM**3+1))/3 - (3**(0.5D0)*Q
     >   2(22,IELEM**3+1))/3 + (15**(0.5D0)*Q2(10,IELEM**3+1))/3 + (5*
     >   *(0.5D0)*Q2(23,IELEM**3+1))/3 + (21**(0.5D0)*Q2(08,IELEM**3+1
     >   ))/3 + (5*3**(0.5D0)*Q2(27,IELEM**3+1))/9 - (5**(0.5D0)*Q2(26
     >   ,IELEM**3+1))/3 - (7**(0.5D0)*Q2(24,IELEM**3+1))/3 - (15**(0.
     >   5D0)*Q2(19,IELEM**3+1))/3 + (77*3**(0.5D0)*Q2(32,IELEM**3+1))
     >   /81 - (11*21**(0.5D0)*Q2(14,IELEM**3+1))/27 + (11*7**(0.5D0)*
     >   Q2(30,IELEM**3+1))/27 - (5**(0.5D0)*Q2(33,IELEM**3+1))/3 + (1
     >   5**(0.5D0)*Q2(25,IELEM**3+1))/9 + (21**(0.5D0)*Q2(20,IELEM**3
     >   +1))/3 + (5*3**(0.5D0)*Q2(39,IELEM**3+1))/9 - (5**(0.5D0)*Q2(
     >   38,IELEM**3+1))/3 - (5*3**(0.5D0)*Q2(42,IELEM**3+1))/9 + (35*
     >   *(0.5D0)*Q2(12,IELEM**3+1))/3 + (5*5**(0.5D0)*Q2(43,IELEM**3+
     >   1))/9 + (15**(0.5D0)*Q2(34,IELEM**3+1))/3 - (11*21**(0.5D0)*Q
     >   2(29,IELEM**3+1))/81 + (11*35**(0.5D0)*Q2(15,IELEM**3+1))/27 
     >   - (5*7**(0.5D0)*Q2(44,IELEM**3+1))/9 + (15**(0.5D0)*Q2(37,IEL
     >   EM**3+1))/9 + (77*5**(0.5D0)*Q2(48,IELEM**3+1))/81 - (55*7**(
     >   0.5D0)*Q2(47,IELEM**3+1))/81 + (11*7**(0.5D0)*Q2(49,IELEM**3+
     >   1))/27 + (77*3**(0.5D0)*Q2(56,IELEM**3+1))/81 + (11*7**(0.5D0
     >   )*Q2(54,IELEM**3+1))/27 - (847*3**(0.5D0)*Q2(62,IELEM**3+1))/
     >   729 + (77*5**(0.5D0)*Q2(60,IELEM**3+1))/81 - (55*7**(0.5D0)*Q
     >   2(59,IELEM**3+1))/81 + (847*5**(0.5D0)*Q2(63,IELEM**3+1))/729
     >    - (847*7**(0.5D0)*Q2(64,IELEM**3+1))/729 - (11*21**(0.5D0)*Q
     >   2(50,IELEM**3+1))/27 + (35**(0.5D0)*Q2(36,IELEM**3+1))/3 - (1
     >   1*21**(0.5D0)*Q2(53,IELEM**3+1))/81 - (11*35**(0.5D0)*Q2(45,I
     >   ELEM**3+1))/81 + (11*35**(0.5D0)*Q2(51,IELEM**3+1))/27 - (11*
     >   35**(0.5D0)*Q2(57,IELEM**3+1))/81 - (105**(0.5D0)*Q2(28,IELEM
     >   **3+1))/9 - (11*105**(0.5D0)*Q2(31,IELEM**3+1))/81 - (105**(0
     >   .5D0)*Q2(40,IELEM**3+1))/9 + (11*105**(0.5D0)*Q2(46,IELEM**3+
     >   1))/81 - (11*105**(0.5D0)*Q2(55,IELEM**3+1))/81 + (11*105**(0
     >   .5D0)*Q2(58,IELEM**3+1))/81) 
      CORNERQ(06,2) = (Q2(01,IELEM**3+1) + Q2(06,IELEM**3+1)/3 + (5
     >   *Q2(11,IELEM**3+1))/9 + (847*Q2(16,IELEM**3+1))/729 + Q2(18,I
     >   ELEM**3+1)/3 + Q2(21,IELEM**3+1)/3 + (5*Q2(35,IELEM**3+1))/9 
     >   + (5*Q2(41,IELEM**3+1))/9 + (847*Q2(52,IELEM**3+1))/729 + (84
     >   7*Q2(61,IELEM**3+1))/729 - (3**(0.5D0)*Q2(02,IELEM**3+1))/3 -
     >    (3**(0.5D0)*Q2(05,IELEM**3+1))/3 - (5**(0.5D0)*Q2(03,IELEM**
     >   3+1))/3 + (11*7**(0.5D0)*Q2(04,IELEM**3+1))/27 - (5**(0.5D0)*
     >   Q2(09,IELEM**3+1))/3 - (3**(0.5D0)*Q2(17,IELEM**3+1))/3 + (11
     >   *7**(0.5D0)*Q2(13,IELEM**3+1))/27 + (15**(0.5D0)*Q2(07,IELEM*
     >   *3+1))/9 - (3**(0.5D0)*Q2(22,IELEM**3+1))/9 + (15**(0.5D0)*Q2
     >   (10,IELEM**3+1))/9 - (5**(0.5D0)*Q2(23,IELEM**3+1))/9 - (11*2
     >   1**(0.5D0)*Q2(08,IELEM**3+1))/81 - (5*3**(0.5D0)*Q2(27,IELEM*
     >   *3+1))/27 - (5**(0.5D0)*Q2(26,IELEM**3+1))/9 + (11*7**(0.5D0)
     >   *Q2(24,IELEM**3+1))/81 + (15**(0.5D0)*Q2(19,IELEM**3+1))/9 - 
     >   (847*3**(0.5D0)*Q2(32,IELEM**3+1))/2187 - (11*21**(0.5D0)*Q2(
     >   14,IELEM**3+1))/81 + (11*7**(0.5D0)*Q2(30,IELEM**3+1))/81 - (
     >   5**(0.5D0)*Q2(33,IELEM**3+1))/3 + (15**(0.5D0)*Q2(25,IELEM**3
     >   +1))/9 - (11*21**(0.5D0)*Q2(20,IELEM**3+1))/81 - (5*3**(0.5D0
     >   )*Q2(39,IELEM**3+1))/27 - (5**(0.5D0)*Q2(38,IELEM**3+1))/9 - 
     >   (5*3**(0.5D0)*Q2(42,IELEM**3+1))/27 - (11*35**(0.5D0)*Q2(12,I
     >   ELEM**3+1))/81 - (5*5**(0.5D0)*Q2(43,IELEM**3+1))/27 + (15**(
     >   0.5D0)*Q2(34,IELEM**3+1))/9 - (11*21**(0.5D0)*Q2(29,IELEM**3+
     >   1))/81 - (11*35**(0.5D0)*Q2(15,IELEM**3+1))/81 + (55*7**(0.5D
     >   0)*Q2(44,IELEM**3+1))/243 + (15**(0.5D0)*Q2(37,IELEM**3+1))/9
     >    - (847*5**(0.5D0)*Q2(48,IELEM**3+1))/2187 + (55*7**(0.5D0)*Q
     >   2(47,IELEM**3+1))/243 + (11*7**(0.5D0)*Q2(49,IELEM**3+1))/27 
     >   - (847*3**(0.5D0)*Q2(56,IELEM**3+1))/2187 + (11*7**(0.5D0)*Q2
     >   (54,IELEM**3+1))/81 - (847*3**(0.5D0)*Q2(62,IELEM**3+1))/2187
     >    - (847*5**(0.5D0)*Q2(60,IELEM**3+1))/2187 + (55*7**(0.5D0)*Q
     >   2(59,IELEM**3+1))/243 - (847*5**(0.5D0)*Q2(63,IELEM**3+1))/21
     >   87 + (9317*7**(0.5D0)*Q2(64,IELEM**3+1))/19683 - (11*21**(0.5
     >   D0)*Q2(50,IELEM**3+1))/81 - (11*35**(0.5D0)*Q2(36,IELEM**3+1)
     >   )/81 - (11*21**(0.5D0)*Q2(53,IELEM**3+1))/81 - (11*35**(0.5D0
     >   )*Q2(45,IELEM**3+1))/81 - (11*35**(0.5D0)*Q2(51,IELEM**3+1))/
     >   81 - (11*35**(0.5D0)*Q2(57,IELEM**3+1))/81 + (11*105**(0.5D0)
     >   *Q2(28,IELEM**3+1))/243 + (11*105**(0.5D0)*Q2(31,IELEM**3+1))
     >   /243 + (11*105**(0.5D0)*Q2(40,IELEM**3+1))/243 + (11*105**(0.
     >   5D0)*Q2(46,IELEM**3+1))/243 + (11*105**(0.5D0)*Q2(55,IELEM**3
     >   +1))/243 + (11*105**(0.5D0)*Q2(58,IELEM**3+1))/243) 
      CORNERQ(07,2) = (Q2(01,IELEM**3+1) - Q2(06,IELEM**3+1)/3 + (5
     >   *Q2(11,IELEM**3+1))/9 - (847*Q2(16,IELEM**3+1))/729 - Q2(18,I
     >   ELEM**3+1)/3 + Q2(21,IELEM**3+1)/3 + (5*Q2(35,IELEM**3+1))/9 
     >   + (5*Q2(41,IELEM**3+1))/9 - (847*Q2(52,IELEM**3+1))/729 + (84
     >   7*Q2(61,IELEM**3+1))/729 + (3**(0.5D0)*Q2(02,IELEM**3+1))/3 -
     >    (3**(0.5D0)*Q2(05,IELEM**3+1))/3 - (5**(0.5D0)*Q2(03,IELEM**
     >   3+1))/3 - (11*7**(0.5D0)*Q2(04,IELEM**3+1))/27 - (5**(0.5D0)*
     >   Q2(09,IELEM**3+1))/3 - (3**(0.5D0)*Q2(17,IELEM**3+1))/3 + (11
     >   *7**(0.5D0)*Q2(13,IELEM**3+1))/27 + (15**(0.5D0)*Q2(07,IELEM*
     >   *3+1))/9 + (3**(0.5D0)*Q2(22,IELEM**3+1))/9 - (15**(0.5D0)*Q2
     >   (10,IELEM**3+1))/9 - (5**(0.5D0)*Q2(23,IELEM**3+1))/9 + (11*2
     >   1**(0.5D0)*Q2(08,IELEM**3+1))/81 - (5*3**(0.5D0)*Q2(27,IELEM*
     >   *3+1))/27 + (5**(0.5D0)*Q2(26,IELEM**3+1))/9 - (11*7**(0.5D0)
     >   *Q2(24,IELEM**3+1))/81 + (15**(0.5D0)*Q2(19,IELEM**3+1))/9 + 
     >   (847*3**(0.5D0)*Q2(32,IELEM**3+1))/2187 + (11*21**(0.5D0)*Q2(
     >   14,IELEM**3+1))/81 - (11*7**(0.5D0)*Q2(30,IELEM**3+1))/81 - (
     >   5**(0.5D0)*Q2(33,IELEM**3+1))/3 + (15**(0.5D0)*Q2(25,IELEM**3
     >   +1))/9 + (11*21**(0.5D0)*Q2(20,IELEM**3+1))/81 - (5*3**(0.5D0
     >   )*Q2(39,IELEM**3+1))/27 + (5**(0.5D0)*Q2(38,IELEM**3+1))/9 + 
     >   (5*3**(0.5D0)*Q2(42,IELEM**3+1))/27 + (11*35**(0.5D0)*Q2(12,I
     >   ELEM**3+1))/81 - (5*5**(0.5D0)*Q2(43,IELEM**3+1))/27 - (15**(
     >   0.5D0)*Q2(34,IELEM**3+1))/9 - (11*21**(0.5D0)*Q2(29,IELEM**3+
     >   1))/81 - (11*35**(0.5D0)*Q2(15,IELEM**3+1))/81 - (55*7**(0.5D
     >   0)*Q2(44,IELEM**3+1))/243 + (15**(0.5D0)*Q2(37,IELEM**3+1))/9
     >    + (847*5**(0.5D0)*Q2(48,IELEM**3+1))/2187 + (55*7**(0.5D0)*Q
     >   2(47,IELEM**3+1))/243 + (11*7**(0.5D0)*Q2(49,IELEM**3+1))/27 
     >   + (847*3**(0.5D0)*Q2(56,IELEM**3+1))/2187 - (11*7**(0.5D0)*Q2
     >   (54,IELEM**3+1))/81 + (847*3**(0.5D0)*Q2(62,IELEM**3+1))/2187
     >    + (847*5**(0.5D0)*Q2(60,IELEM**3+1))/2187 + (55*7**(0.5D0)*Q
     >   2(59,IELEM**3+1))/243 - (847*5**(0.5D0)*Q2(63,IELEM**3+1))/21
     >   87 - (9317*7**(0.5D0)*Q2(64,IELEM**3+1))/19683 + (11*21**(0.5
     >   D0)*Q2(50,IELEM**3+1))/81 + (11*35**(0.5D0)*Q2(36,IELEM**3+1)
     >   )/81 - (11*21**(0.5D0)*Q2(53,IELEM**3+1))/81 - (11*35**(0.5D0
     >   )*Q2(45,IELEM**3+1))/81 - (11*35**(0.5D0)*Q2(51,IELEM**3+1))/
     >   81 - (11*35**(0.5D0)*Q2(57,IELEM**3+1))/81 - (11*105**(0.5D0)
     >   *Q2(28,IELEM**3+1))/243 + (11*105**(0.5D0)*Q2(31,IELEM**3+1))
     >   /243 - (11*105**(0.5D0)*Q2(40,IELEM**3+1))/243 - (11*105**(0.
     >   5D0)*Q2(46,IELEM**3+1))/243 + (11*105**(0.5D0)*Q2(55,IELEM**3
     >   +1))/243 - (11*105**(0.5D0)*Q2(58,IELEM**3+1))/243) 
      CORNERQ(08,2) = (Q2(01,IELEM**3+1) - Q2(06,IELEM**3+1) - (5*Q
     >   2(11,IELEM**3+1))/3 + (77*Q2(16,IELEM**3+1))/27 - Q2(18,IELEM
     >   **3+1) + Q2(21,IELEM**3+1)/3 - (5*Q2(35,IELEM**3+1))/3 + (5*Q
     >   2(41,IELEM**3+1))/9 + (77*Q2(52,IELEM**3+1))/27 + (847*Q2(61,
     >   IELEM**3+1))/729 + 3**(0.5D0)*Q2(02,IELEM**3+1) - (3**(0.5D0)
     >   *Q2(05,IELEM**3+1))/3 + 5**(0.5D0)*Q2(03,IELEM**3+1) + 7**(0.
     >   5D0)*Q2(04,IELEM**3+1) - (5**(0.5D0)*Q2(09,IELEM**3+1))/3 - (
     >   3**(0.5D0)*Q2(17,IELEM**3+1))/3 + (11*7**(0.5D0)*Q2(13,IELEM*
     >   *3+1))/27 - (15**(0.5D0)*Q2(07,IELEM**3+1))/3 + (3**(0.5D0)*Q
     >   2(22,IELEM**3+1))/3 - (15**(0.5D0)*Q2(10,IELEM**3+1))/3 + (5*
     >   *(0.5D0)*Q2(23,IELEM**3+1))/3 - (21**(0.5D0)*Q2(08,IELEM**3+1
     >   ))/3 + (5*3**(0.5D0)*Q2(27,IELEM**3+1))/9 + (5**(0.5D0)*Q2(26
     >   ,IELEM**3+1))/3 + (7**(0.5D0)*Q2(24,IELEM**3+1))/3 - (15**(0.
     >   5D0)*Q2(19,IELEM**3+1))/3 - (77*3**(0.5D0)*Q2(32,IELEM**3+1))
     >   /81 + (11*21**(0.5D0)*Q2(14,IELEM**3+1))/27 - (11*7**(0.5D0)*
     >   Q2(30,IELEM**3+1))/27 - (5**(0.5D0)*Q2(33,IELEM**3+1))/3 + (1
     >   5**(0.5D0)*Q2(25,IELEM**3+1))/9 - (21**(0.5D0)*Q2(20,IELEM**3
     >   +1))/3 + (5*3**(0.5D0)*Q2(39,IELEM**3+1))/9 + (5**(0.5D0)*Q2(
     >   38,IELEM**3+1))/3 + (5*3**(0.5D0)*Q2(42,IELEM**3+1))/9 - (35*
     >   *(0.5D0)*Q2(12,IELEM**3+1))/3 + (5*5**(0.5D0)*Q2(43,IELEM**3+
     >   1))/9 - (15**(0.5D0)*Q2(34,IELEM**3+1))/3 - (11*21**(0.5D0)*Q
     >   2(29,IELEM**3+1))/81 + (11*35**(0.5D0)*Q2(15,IELEM**3+1))/27 
     >   + (5*7**(0.5D0)*Q2(44,IELEM**3+1))/9 + (15**(0.5D0)*Q2(37,IEL
     >   EM**3+1))/9 - (77*5**(0.5D0)*Q2(48,IELEM**3+1))/81 - (55*7**(
     >   0.5D0)*Q2(47,IELEM**3+1))/81 + (11*7**(0.5D0)*Q2(49,IELEM**3+
     >   1))/27 - (77*3**(0.5D0)*Q2(56,IELEM**3+1))/81 - (11*7**(0.5D0
     >   )*Q2(54,IELEM**3+1))/27 + (847*3**(0.5D0)*Q2(62,IELEM**3+1))/
     >   729 - (77*5**(0.5D0)*Q2(60,IELEM**3+1))/81 - (55*7**(0.5D0)*Q
     >   2(59,IELEM**3+1))/81 + (847*5**(0.5D0)*Q2(63,IELEM**3+1))/729
     >    + (847*7**(0.5D0)*Q2(64,IELEM**3+1))/729 + (11*21**(0.5D0)*Q
     >   2(50,IELEM**3+1))/27 - (35**(0.5D0)*Q2(36,IELEM**3+1))/3 - (1
     >   1*21**(0.5D0)*Q2(53,IELEM**3+1))/81 - (11*35**(0.5D0)*Q2(45,I
     >   ELEM**3+1))/81 + (11*35**(0.5D0)*Q2(51,IELEM**3+1))/27 - (11*
     >   35**(0.5D0)*Q2(57,IELEM**3+1))/81 + (105**(0.5D0)*Q2(28,IELEM
     >   **3+1))/9 - (11*105**(0.5D0)*Q2(31,IELEM**3+1))/81 + (105**(0
     >   .5D0)*Q2(40,IELEM**3+1))/9 - (11*105**(0.5D0)*Q2(46,IELEM**3+
     >   1))/81 - (11*105**(0.5D0)*Q2(55,IELEM**3+1))/81 - (11*105**(0
     >   .5D0)*Q2(58,IELEM**3+1))/81) 
      CORNERQ(09,2) = (Q2(01,IELEM**3+1) - Q2(06,IELEM**3+1) - (5*Q
     >   2(11,IELEM**3+1))/3 + (77*Q2(16,IELEM**3+1))/27 + Q2(18,IELEM
     >   **3+1) - Q2(21,IELEM**3+1)/3 - (5*Q2(35,IELEM**3+1))/3 + (5*Q
     >   2(41,IELEM**3+1))/9 - (77*Q2(52,IELEM**3+1))/27 - (847*Q2(61,
     >   IELEM**3+1))/729 - 3**(0.5D0)*Q2(02,IELEM**3+1) + (3**(0.5D0)
     >   *Q2(05,IELEM**3+1))/3 + 5**(0.5D0)*Q2(03,IELEM**3+1) - 7**(0.
     >   5D0)*Q2(04,IELEM**3+1) - (5**(0.5D0)*Q2(09,IELEM**3+1))/3 - (
     >   3**(0.5D0)*Q2(17,IELEM**3+1))/3 - (11*7**(0.5D0)*Q2(13,IELEM*
     >   *3+1))/27 + (15**(0.5D0)*Q2(07,IELEM**3+1))/3 + (3**(0.5D0)*Q
     >   2(22,IELEM**3+1))/3 + (15**(0.5D0)*Q2(10,IELEM**3+1))/3 - (5*
     >   *(0.5D0)*Q2(23,IELEM**3+1))/3 - (21**(0.5D0)*Q2(08,IELEM**3+1
     >   ))/3 + (5*3**(0.5D0)*Q2(27,IELEM**3+1))/9 - (5**(0.5D0)*Q2(26
     >   ,IELEM**3+1))/3 + (7**(0.5D0)*Q2(24,IELEM**3+1))/3 - (15**(0.
     >   5D0)*Q2(19,IELEM**3+1))/3 - (77*3**(0.5D0)*Q2(32,IELEM**3+1))
     >   /81 + (11*21**(0.5D0)*Q2(14,IELEM**3+1))/27 - (11*7**(0.5D0)*
     >   Q2(30,IELEM**3+1))/27 - (5**(0.5D0)*Q2(33,IELEM**3+1))/3 + (1
     >   5**(0.5D0)*Q2(25,IELEM**3+1))/9 + (21**(0.5D0)*Q2(20,IELEM**3
     >   +1))/3 - (5*3**(0.5D0)*Q2(39,IELEM**3+1))/9 + (5**(0.5D0)*Q2(
     >   38,IELEM**3+1))/3 - (5*3**(0.5D0)*Q2(42,IELEM**3+1))/9 + (35*
     >   *(0.5D0)*Q2(12,IELEM**3+1))/3 + (5*5**(0.5D0)*Q2(43,IELEM**3+
     >   1))/9 + (15**(0.5D0)*Q2(34,IELEM**3+1))/3 + (11*21**(0.5D0)*Q
     >   2(29,IELEM**3+1))/81 - (11*35**(0.5D0)*Q2(15,IELEM**3+1))/27 
     >   - (5*7**(0.5D0)*Q2(44,IELEM**3+1))/9 - (15**(0.5D0)*Q2(37,IEL
     >   EM**3+1))/9 - (77*5**(0.5D0)*Q2(48,IELEM**3+1))/81 + (55*7**(
     >   0.5D0)*Q2(47,IELEM**3+1))/81 + (11*7**(0.5D0)*Q2(49,IELEM**3+
     >   1))/27 - (77*3**(0.5D0)*Q2(56,IELEM**3+1))/81 - (11*7**(0.5D0
     >   )*Q2(54,IELEM**3+1))/27 + (847*3**(0.5D0)*Q2(62,IELEM**3+1))/
     >   729 + (77*5**(0.5D0)*Q2(60,IELEM**3+1))/81 - (55*7**(0.5D0)*Q
     >   2(59,IELEM**3+1))/81 - (847*5**(0.5D0)*Q2(63,IELEM**3+1))/729
     >    + (847*7**(0.5D0)*Q2(64,IELEM**3+1))/729 - (11*21**(0.5D0)*Q
     >   2(50,IELEM**3+1))/27 + (35**(0.5D0)*Q2(36,IELEM**3+1))/3 + (1
     >   1*21**(0.5D0)*Q2(53,IELEM**3+1))/81 + (11*35**(0.5D0)*Q2(45,I
     >   ELEM**3+1))/81 + (11*35**(0.5D0)*Q2(51,IELEM**3+1))/27 - (11*
     >   35**(0.5D0)*Q2(57,IELEM**3+1))/81 - (105**(0.5D0)*Q2(28,IELEM
     >   **3+1))/9 + (11*105**(0.5D0)*Q2(31,IELEM**3+1))/81 + (105**(0
     >   .5D0)*Q2(40,IELEM**3+1))/9 - (11*105**(0.5D0)*Q2(46,IELEM**3+
     >   1))/81 + (11*105**(0.5D0)*Q2(55,IELEM**3+1))/81 + (11*105**(0
     >   .5D0)*Q2(58,IELEM**3+1))/81) 
      CORNERQ(10,2) = (Q2(01,IELEM**3+1) - Q2(06,IELEM**3+1)/3 + (5
     >   *Q2(11,IELEM**3+1))/9 - (847*Q2(16,IELEM**3+1))/729 + Q2(18,I
     >   ELEM**3+1)/3 - Q2(21,IELEM**3+1)/3 + (5*Q2(35,IELEM**3+1))/9 
     >   + (5*Q2(41,IELEM**3+1))/9 + (847*Q2(52,IELEM**3+1))/729 - (84
     >   7*Q2(61,IELEM**3+1))/729 - (3**(0.5D0)*Q2(02,IELEM**3+1))/3 +
     >    (3**(0.5D0)*Q2(05,IELEM**3+1))/3 - (5**(0.5D0)*Q2(03,IELEM**
     >   3+1))/3 + (11*7**(0.5D0)*Q2(04,IELEM**3+1))/27 - (5**(0.5D0)*
     >   Q2(09,IELEM**3+1))/3 - (3**(0.5D0)*Q2(17,IELEM**3+1))/3 - (11
     >   *7**(0.5D0)*Q2(13,IELEM**3+1))/27 - (15**(0.5D0)*Q2(07,IELEM*
     >   *3+1))/9 + (3**(0.5D0)*Q2(22,IELEM**3+1))/9 + (15**(0.5D0)*Q2
     >   (10,IELEM**3+1))/9 + (5**(0.5D0)*Q2(23,IELEM**3+1))/9 + (11*2
     >   1**(0.5D0)*Q2(08,IELEM**3+1))/81 - (5*3**(0.5D0)*Q2(27,IELEM*
     >   *3+1))/27 - (5**(0.5D0)*Q2(26,IELEM**3+1))/9 - (11*7**(0.5D0)
     >   *Q2(24,IELEM**3+1))/81 + (15**(0.5D0)*Q2(19,IELEM**3+1))/9 + 
     >   (847*3**(0.5D0)*Q2(32,IELEM**3+1))/2187 + (11*21**(0.5D0)*Q2(
     >   14,IELEM**3+1))/81 - (11*7**(0.5D0)*Q2(30,IELEM**3+1))/81 - (
     >   5**(0.5D0)*Q2(33,IELEM**3+1))/3 + (15**(0.5D0)*Q2(25,IELEM**3
     >   +1))/9 - (11*21**(0.5D0)*Q2(20,IELEM**3+1))/81 + (5*3**(0.5D0
     >   )*Q2(39,IELEM**3+1))/27 + (5**(0.5D0)*Q2(38,IELEM**3+1))/9 - 
     >   (5*3**(0.5D0)*Q2(42,IELEM**3+1))/27 - (11*35**(0.5D0)*Q2(12,I
     >   ELEM**3+1))/81 - (5*5**(0.5D0)*Q2(43,IELEM**3+1))/27 + (15**(
     >   0.5D0)*Q2(34,IELEM**3+1))/9 + (11*21**(0.5D0)*Q2(29,IELEM**3+
     >   1))/81 + (11*35**(0.5D0)*Q2(15,IELEM**3+1))/81 + (55*7**(0.5D
     >   0)*Q2(44,IELEM**3+1))/243 - (15**(0.5D0)*Q2(37,IELEM**3+1))/9
     >    + (847*5**(0.5D0)*Q2(48,IELEM**3+1))/2187 - (55*7**(0.5D0)*Q
     >   2(47,IELEM**3+1))/243 + (11*7**(0.5D0)*Q2(49,IELEM**3+1))/27 
     >   + (847*3**(0.5D0)*Q2(56,IELEM**3+1))/2187 - (11*7**(0.5D0)*Q2
     >   (54,IELEM**3+1))/81 + (847*3**(0.5D0)*Q2(62,IELEM**3+1))/2187
     >    - (847*5**(0.5D0)*Q2(60,IELEM**3+1))/2187 + (55*7**(0.5D0)*Q
     >   2(59,IELEM**3+1))/243 + (847*5**(0.5D0)*Q2(63,IELEM**3+1))/21
     >   87 - (9317*7**(0.5D0)*Q2(64,IELEM**3+1))/19683 - (11*21**(0.5
     >   D0)*Q2(50,IELEM**3+1))/81 - (11*35**(0.5D0)*Q2(36,IELEM**3+1)
     >   )/81 + (11*21**(0.5D0)*Q2(53,IELEM**3+1))/81 + (11*35**(0.5D0
     >   )*Q2(45,IELEM**3+1))/81 - (11*35**(0.5D0)*Q2(51,IELEM**3+1))/
     >   81 - (11*35**(0.5D0)*Q2(57,IELEM**3+1))/81 + (11*105**(0.5D0)
     >   *Q2(28,IELEM**3+1))/243 - (11*105**(0.5D0)*Q2(31,IELEM**3+1))
     >   /243 - (11*105**(0.5D0)*Q2(40,IELEM**3+1))/243 - (11*105**(0.
     >   5D0)*Q2(46,IELEM**3+1))/243 - (11*105**(0.5D0)*Q2(55,IELEM**3
     >   +1))/243 + (11*105**(0.5D0)*Q2(58,IELEM**3+1))/243) 
      CORNERQ(11,2) = (Q2(01,IELEM**3+1) + Q2(06,IELEM**3+1)/3 + (5
     >   *Q2(11,IELEM**3+1))/9 + (847*Q2(16,IELEM**3+1))/729 - Q2(18,I
     >   ELEM**3+1)/3 - Q2(21,IELEM**3+1)/3 + (5*Q2(35,IELEM**3+1))/9 
     >   + (5*Q2(41,IELEM**3+1))/9 - (847*Q2(52,IELEM**3+1))/729 - (84
     >   7*Q2(61,IELEM**3+1))/729 + (3**(0.5D0)*Q2(02,IELEM**3+1))/3 +
     >    (3**(0.5D0)*Q2(05,IELEM**3+1))/3 - (5**(0.5D0)*Q2(03,IELEM**
     >   3+1))/3 - (11*7**(0.5D0)*Q2(04,IELEM**3+1))/27 - (5**(0.5D0)*
     >   Q2(09,IELEM**3+1))/3 - (3**(0.5D0)*Q2(17,IELEM**3+1))/3 - (11
     >   *7**(0.5D0)*Q2(13,IELEM**3+1))/27 - (15**(0.5D0)*Q2(07,IELEM*
     >   *3+1))/9 - (3**(0.5D0)*Q2(22,IELEM**3+1))/9 - (15**(0.5D0)*Q2
     >   (10,IELEM**3+1))/9 + (5**(0.5D0)*Q2(23,IELEM**3+1))/9 - (11*2
     >   1**(0.5D0)*Q2(08,IELEM**3+1))/81 - (5*3**(0.5D0)*Q2(27,IELEM*
     >   *3+1))/27 + (5**(0.5D0)*Q2(26,IELEM**3+1))/9 + (11*7**(0.5D0)
     >   *Q2(24,IELEM**3+1))/81 + (15**(0.5D0)*Q2(19,IELEM**3+1))/9 - 
     >   (847*3**(0.5D0)*Q2(32,IELEM**3+1))/2187 - (11*21**(0.5D0)*Q2(
     >   14,IELEM**3+1))/81 + (11*7**(0.5D0)*Q2(30,IELEM**3+1))/81 - (
     >   5**(0.5D0)*Q2(33,IELEM**3+1))/3 + (15**(0.5D0)*Q2(25,IELEM**3
     >   +1))/9 + (11*21**(0.5D0)*Q2(20,IELEM**3+1))/81 + (5*3**(0.5D0
     >   )*Q2(39,IELEM**3+1))/27 - (5**(0.5D0)*Q2(38,IELEM**3+1))/9 + 
     >   (5*3**(0.5D0)*Q2(42,IELEM**3+1))/27 + (11*35**(0.5D0)*Q2(12,I
     >   ELEM**3+1))/81 - (5*5**(0.5D0)*Q2(43,IELEM**3+1))/27 - (15**(
     >   0.5D0)*Q2(34,IELEM**3+1))/9 + (11*21**(0.5D0)*Q2(29,IELEM**3+
     >   1))/81 + (11*35**(0.5D0)*Q2(15,IELEM**3+1))/81 - (55*7**(0.5D
     >   0)*Q2(44,IELEM**3+1))/243 - (15**(0.5D0)*Q2(37,IELEM**3+1))/9
     >    - (847*5**(0.5D0)*Q2(48,IELEM**3+1))/2187 - (55*7**(0.5D0)*Q
     >   2(47,IELEM**3+1))/243 + (11*7**(0.5D0)*Q2(49,IELEM**3+1))/27 
     >   - (847*3**(0.5D0)*Q2(56,IELEM**3+1))/2187 + (11*7**(0.5D0)*Q2
     >   (54,IELEM**3+1))/81 - (847*3**(0.5D0)*Q2(62,IELEM**3+1))/2187
     >    + (847*5**(0.5D0)*Q2(60,IELEM**3+1))/2187 + (55*7**(0.5D0)*Q
     >   2(59,IELEM**3+1))/243 + (847*5**(0.5D0)*Q2(63,IELEM**3+1))/21
     >   87 + (9317*7**(0.5D0)*Q2(64,IELEM**3+1))/19683 + (11*21**(0.5
     >   D0)*Q2(50,IELEM**3+1))/81 + (11*35**(0.5D0)*Q2(36,IELEM**3+1)
     >   )/81 + (11*21**(0.5D0)*Q2(53,IELEM**3+1))/81 + (11*35**(0.5D0
     >   )*Q2(45,IELEM**3+1))/81 - (11*35**(0.5D0)*Q2(51,IELEM**3+1))/
     >   81 - (11*35**(0.5D0)*Q2(57,IELEM**3+1))/81 - (11*105**(0.5D0)
     >   *Q2(28,IELEM**3+1))/243 - (11*105**(0.5D0)*Q2(31,IELEM**3+1))
     >   /243 + (11*105**(0.5D0)*Q2(40,IELEM**3+1))/243 + (11*105**(0.
     >   5D0)*Q2(46,IELEM**3+1))/243 - (11*105**(0.5D0)*Q2(55,IELEM**3
     >   +1))/243 - (11*105**(0.5D0)*Q2(58,IELEM**3+1))/243) 
      CORNERQ(12,2) = (Q2(01,IELEM**3+1) + Q2(06,IELEM**3+1) - (5*Q
     >   2(11,IELEM**3+1))/3 - (77*Q2(16,IELEM**3+1))/27 - Q2(18,IELEM
     >   **3+1) - Q2(21,IELEM**3+1)/3 - (5*Q2(35,IELEM**3+1))/3 + (5*Q
     >   2(41,IELEM**3+1))/9 + (77*Q2(52,IELEM**3+1))/27 - (847*Q2(61,
     >   IELEM**3+1))/729 + 3**(0.5D0)*Q2(02,IELEM**3+1) + (3**(0.5D0)
     >   *Q2(05,IELEM**3+1))/3 + 5**(0.5D0)*Q2(03,IELEM**3+1) + 7**(0.
     >   5D0)*Q2(04,IELEM**3+1) - (5**(0.5D0)*Q2(09,IELEM**3+1))/3 - (
     >   3**(0.5D0)*Q2(17,IELEM**3+1))/3 - (11*7**(0.5D0)*Q2(13,IELEM*
     >   *3+1))/27 + (15**(0.5D0)*Q2(07,IELEM**3+1))/3 - (3**(0.5D0)*Q
     >   2(22,IELEM**3+1))/3 - (15**(0.5D0)*Q2(10,IELEM**3+1))/3 - (5*
     >   *(0.5D0)*Q2(23,IELEM**3+1))/3 + (21**(0.5D0)*Q2(08,IELEM**3+1
     >   ))/3 + (5*3**(0.5D0)*Q2(27,IELEM**3+1))/9 + (5**(0.5D0)*Q2(26
     >   ,IELEM**3+1))/3 - (7**(0.5D0)*Q2(24,IELEM**3+1))/3 - (15**(0.
     >   5D0)*Q2(19,IELEM**3+1))/3 + (77*3**(0.5D0)*Q2(32,IELEM**3+1))
     >   /81 - (11*21**(0.5D0)*Q2(14,IELEM**3+1))/27 + (11*7**(0.5D0)*
     >   Q2(30,IELEM**3+1))/27 - (5**(0.5D0)*Q2(33,IELEM**3+1))/3 + (1
     >   5**(0.5D0)*Q2(25,IELEM**3+1))/9 - (21**(0.5D0)*Q2(20,IELEM**3
     >   +1))/3 - (5*3**(0.5D0)*Q2(39,IELEM**3+1))/9 - (5**(0.5D0)*Q2(
     >   38,IELEM**3+1))/3 + (5*3**(0.5D0)*Q2(42,IELEM**3+1))/9 - (35*
     >   *(0.5D0)*Q2(12,IELEM**3+1))/3 + (5*5**(0.5D0)*Q2(43,IELEM**3+
     >   1))/9 - (15**(0.5D0)*Q2(34,IELEM**3+1))/3 + (11*21**(0.5D0)*Q
     >   2(29,IELEM**3+1))/81 - (11*35**(0.5D0)*Q2(15,IELEM**3+1))/27 
     >   + (5*7**(0.5D0)*Q2(44,IELEM**3+1))/9 - (15**(0.5D0)*Q2(37,IEL
     >   EM**3+1))/9 + (77*5**(0.5D0)*Q2(48,IELEM**3+1))/81 + (55*7**(
     >   0.5D0)*Q2(47,IELEM**3+1))/81 + (11*7**(0.5D0)*Q2(49,IELEM**3+
     >   1))/27 + (77*3**(0.5D0)*Q2(56,IELEM**3+1))/81 + (11*7**(0.5D0
     >   )*Q2(54,IELEM**3+1))/27 - (847*3**(0.5D0)*Q2(62,IELEM**3+1))/
     >   729 - (77*5**(0.5D0)*Q2(60,IELEM**3+1))/81 - (55*7**(0.5D0)*Q
     >   2(59,IELEM**3+1))/81 - (847*5**(0.5D0)*Q2(63,IELEM**3+1))/729
     >    - (847*7**(0.5D0)*Q2(64,IELEM**3+1))/729 + (11*21**(0.5D0)*Q
     >   2(50,IELEM**3+1))/27 - (35**(0.5D0)*Q2(36,IELEM**3+1))/3 + (1
     >   1*21**(0.5D0)*Q2(53,IELEM**3+1))/81 + (11*35**(0.5D0)*Q2(45,I
     >   ELEM**3+1))/81 + (11*35**(0.5D0)*Q2(51,IELEM**3+1))/27 - (11*
     >   35**(0.5D0)*Q2(57,IELEM**3+1))/81 + (105**(0.5D0)*Q2(28,IELEM
     >   **3+1))/9 + (11*105**(0.5D0)*Q2(31,IELEM**3+1))/81 - (105**(0
     >   .5D0)*Q2(40,IELEM**3+1))/9 + (11*105**(0.5D0)*Q2(46,IELEM**3+
     >   1))/81 + (11*105**(0.5D0)*Q2(55,IELEM**3+1))/81 - (11*105**(0
     >   .5D0)*Q2(58,IELEM**3+1))/81) 
      CORNERQ(13,2) = (Q2(01,IELEM**3+1) - 3*Q2(06,IELEM**3+1) + 5*
     >   Q2(11,IELEM**3+1) - 7*Q2(16,IELEM**3+1) + Q2(18,IELEM**3+1) -
     >    Q2(21,IELEM**3+1) - (5*Q2(35,IELEM**3+1))/3 - (5*Q2(41,IELEM
     >   **3+1))/3 - (77*Q2(52,IELEM**3+1))/27 + (77*Q2(61,IELEM**3+1)
     >   )/27 - 3**(0.5D0)*Q2(02,IELEM**3+1) + 3**(0.5D0)*Q2(05,IELEM*
     >   *3+1) + 5**(0.5D0)*Q2(03,IELEM**3+1) - 7**(0.5D0)*Q2(04,IELEM
     >   **3+1) + 5**(0.5D0)*Q2(09,IELEM**3+1) - (3**(0.5D0)*Q2(17,IEL
     >   EM**3+1))/3 + 7**(0.5D0)*Q2(13,IELEM**3+1) + 15**(0.5D0)*Q2(0
     >   7,IELEM**3+1) + 3**(0.5D0)*Q2(22,IELEM**3+1) - 15**(0.5D0)*Q2
     >   (10,IELEM**3+1) - 5**(0.5D0)*Q2(23,IELEM**3+1) - 21**(0.5D0)*
     >   Q2(08,IELEM**3+1) - (5*3**(0.5D0)*Q2(27,IELEM**3+1))/3 + 5**(
     >   0.5D0)*Q2(26,IELEM**3+1) + 7**(0.5D0)*Q2(24,IELEM**3+1) - (15
     >   **(0.5D0)*Q2(19,IELEM**3+1))/3 + (7*3**(0.5D0)*Q2(32,IELEM**3
     >   +1))/3 - 21**(0.5D0)*Q2(14,IELEM**3+1) + 7**(0.5D0)*Q2(30,IEL
     >   EM**3+1) - (5**(0.5D0)*Q2(33,IELEM**3+1))/3 - (15**(0.5D0)*Q2
     >   (25,IELEM**3+1))/3 + (21**(0.5D0)*Q2(20,IELEM**3+1))/3 - (5*3
     >   **(0.5D0)*Q2(39,IELEM**3+1))/3 + 5**(0.5D0)*Q2(38,IELEM**3+1)
     >    + (5*3**(0.5D0)*Q2(42,IELEM**3+1))/3 - 35**(0.5D0)*Q2(12,IEL
     >   EM**3+1) - (5*5**(0.5D0)*Q2(43,IELEM**3+1))/3 + (15**(0.5D0)*
     >   Q2(34,IELEM**3+1))/3 - (21**(0.5D0)*Q2(29,IELEM**3+1))/3 + 35
     >   **(0.5D0)*Q2(15,IELEM**3+1) + (5*7**(0.5D0)*Q2(44,IELEM**3+1)
     >   )/3 - (15**(0.5D0)*Q2(37,IELEM**3+1))/3 + (7*5**(0.5D0)*Q2(48
     >   ,IELEM**3+1))/3 - (5*7**(0.5D0)*Q2(47,IELEM**3+1))/3 + (11*7*
     >   *(0.5D0)*Q2(49,IELEM**3+1))/27 - (77*3**(0.5D0)*Q2(56,IELEM**
     >   3+1))/27 - (11*7**(0.5D0)*Q2(54,IELEM**3+1))/9 - (77*3**(0.5D
     >   0)*Q2(62,IELEM**3+1))/27 - (77*5**(0.5D0)*Q2(60,IELEM**3+1))/
     >   27 + (55*7**(0.5D0)*Q2(59,IELEM**3+1))/27 + (77*5**(0.5D0)*Q2
     >   (63,IELEM**3+1))/27 - (77*7**(0.5D0)*Q2(64,IELEM**3+1))/27 - 
     >   (11*21**(0.5D0)*Q2(50,IELEM**3+1))/27 + (35**(0.5D0)*Q2(36,IE
     >   LEM**3+1))/3 + (11*21**(0.5D0)*Q2(53,IELEM**3+1))/27 - (35**(
     >   0.5D0)*Q2(45,IELEM**3+1))/3 + (11*35**(0.5D0)*Q2(51,IELEM**3+
     >   1))/27 + (11*35**(0.5D0)*Q2(57,IELEM**3+1))/27 + (105**(0.5D0
     >   )*Q2(28,IELEM**3+1))/3 - (105**(0.5D0)*Q2(31,IELEM**3+1))/3 +
     >    (105**(0.5D0)*Q2(40,IELEM**3+1))/3 + (105**(0.5D0)*Q2(46,IEL
     >   EM**3+1))/3 + (11*105**(0.5D0)*Q2(55,IELEM**3+1))/27 - (11*10
     >   5**(0.5D0)*Q2(58,IELEM**3+1))/27) 
      CORNERQ(14,2) = (Q2(01,IELEM**3+1) - Q2(06,IELEM**3+1) - (5*Q
     >   2(11,IELEM**3+1))/3 + (77*Q2(16,IELEM**3+1))/27 + Q2(18,IELEM
     >   **3+1)/3 - Q2(21,IELEM**3+1) + (5*Q2(35,IELEM**3+1))/9 - (5*Q
     >   2(41,IELEM**3+1))/3 + (847*Q2(52,IELEM**3+1))/729 + (77*Q2(61
     >   ,IELEM**3+1))/27 - (3**(0.5D0)*Q2(02,IELEM**3+1))/3 + 3**(0.5
     >   D0)*Q2(05,IELEM**3+1) - (5**(0.5D0)*Q2(03,IELEM**3+1))/3 + (1
     >   1*7**(0.5D0)*Q2(04,IELEM**3+1))/27 + 5**(0.5D0)*Q2(09,IELEM**
     >   3+1) - (3**(0.5D0)*Q2(17,IELEM**3+1))/3 + 7**(0.5D0)*Q2(13,IE
     >   LEM**3+1) - (15**(0.5D0)*Q2(07,IELEM**3+1))/3 + (3**(0.5D0)*Q
     >   2(22,IELEM**3+1))/3 - (15**(0.5D0)*Q2(10,IELEM**3+1))/3 + (5*
     >   *(0.5D0)*Q2(23,IELEM**3+1))/3 + (11*21**(0.5D0)*Q2(08,IELEM**
     >   3+1))/27 + (5*3**(0.5D0)*Q2(27,IELEM**3+1))/9 + (5**(0.5D0)*Q
     >   2(26,IELEM**3+1))/3 - (11*7**(0.5D0)*Q2(24,IELEM**3+1))/27 + 
     >   (15**(0.5D0)*Q2(19,IELEM**3+1))/9 - (77*3**(0.5D0)*Q2(32,IELE
     >   M**3+1))/81 - (21**(0.5D0)*Q2(14,IELEM**3+1))/3 + (7**(0.5D0)
     >   *Q2(30,IELEM**3+1))/3 - (5**(0.5D0)*Q2(33,IELEM**3+1))/3 - (1
     >   5**(0.5D0)*Q2(25,IELEM**3+1))/3 - (11*21**(0.5D0)*Q2(20,IELEM
     >   **3+1))/81 + (5*3**(0.5D0)*Q2(39,IELEM**3+1))/9 + (5**(0.5D0)
     >   *Q2(38,IELEM**3+1))/3 + (5*3**(0.5D0)*Q2(42,IELEM**3+1))/9 + 
     >   (11*35**(0.5D0)*Q2(12,IELEM**3+1))/27 + (5*5**(0.5D0)*Q2(43,I
     >   ELEM**3+1))/9 + (15**(0.5D0)*Q2(34,IELEM**3+1))/9 - (21**(0.5
     >   D0)*Q2(29,IELEM**3+1))/3 - (35**(0.5D0)*Q2(15,IELEM**3+1))/3 
     >   - (55*7**(0.5D0)*Q2(44,IELEM**3+1))/81 - (15**(0.5D0)*Q2(37,I
     >   ELEM**3+1))/3 - (77*5**(0.5D0)*Q2(48,IELEM**3+1))/81 + (5*7**
     >   (0.5D0)*Q2(47,IELEM**3+1))/9 + (11*7**(0.5D0)*Q2(49,IELEM**3+
     >   1))/27 + (847*3**(0.5D0)*Q2(56,IELEM**3+1))/729 - (11*7**(0.5
     >   D0)*Q2(54,IELEM**3+1))/27 - (77*3**(0.5D0)*Q2(62,IELEM**3+1))
     >   /81 + (847*5**(0.5D0)*Q2(60,IELEM**3+1))/729 - (55*7**(0.5D0)
     >   *Q2(59,IELEM**3+1))/81 - (77*5**(0.5D0)*Q2(63,IELEM**3+1))/81
     >    + (847*7**(0.5D0)*Q2(64,IELEM**3+1))/729 - (11*21**(0.5D0)*Q
     >   2(50,IELEM**3+1))/81 - (11*35**(0.5D0)*Q2(36,IELEM**3+1))/81 
     >   + (11*21**(0.5D0)*Q2(53,IELEM**3+1))/27 - (35**(0.5D0)*Q2(45,
     >   IELEM**3+1))/3 - (11*35**(0.5D0)*Q2(51,IELEM**3+1))/81 + (11*
     >   35**(0.5D0)*Q2(57,IELEM**3+1))/27 - (11*105**(0.5D0)*Q2(28,IE
     >   LEM**3+1))/81 + (105**(0.5D0)*Q2(31,IELEM**3+1))/9 - (11*105*
     >   *(0.5D0)*Q2(40,IELEM**3+1))/81 + (105**(0.5D0)*Q2(46,IELEM**3
     >   +1))/9 - (11*105**(0.5D0)*Q2(55,IELEM**3+1))/81 - (11*105**(0
     >   .5D0)*Q2(58,IELEM**3+1))/81) 
      CORNERQ(15,2) = (Q2(01,IELEM**3+1) + Q2(06,IELEM**3+1) - (5*Q
     >   2(11,IELEM**3+1))/3 - (77*Q2(16,IELEM**3+1))/27 - Q2(18,IELEM
     >   **3+1)/3 - Q2(21,IELEM**3+1) + (5*Q2(35,IELEM**3+1))/9 - (5*Q
     >   2(41,IELEM**3+1))/3 - (847*Q2(52,IELEM**3+1))/729 + (77*Q2(61
     >   ,IELEM**3+1))/27 + (3**(0.5D0)*Q2(02,IELEM**3+1))/3 + 3**(0.5
     >   D0)*Q2(05,IELEM**3+1) - (5**(0.5D0)*Q2(03,IELEM**3+1))/3 - (1
     >   1*7**(0.5D0)*Q2(04,IELEM**3+1))/27 + 5**(0.5D0)*Q2(09,IELEM**
     >   3+1) - (3**(0.5D0)*Q2(17,IELEM**3+1))/3 + 7**(0.5D0)*Q2(13,IE
     >   LEM**3+1) - (15**(0.5D0)*Q2(07,IELEM**3+1))/3 - (3**(0.5D0)*Q
     >   2(22,IELEM**3+1))/3 + (15**(0.5D0)*Q2(10,IELEM**3+1))/3 + (5*
     >   *(0.5D0)*Q2(23,IELEM**3+1))/3 - (11*21**(0.5D0)*Q2(08,IELEM**
     >   3+1))/27 + (5*3**(0.5D0)*Q2(27,IELEM**3+1))/9 - (5**(0.5D0)*Q
     >   2(26,IELEM**3+1))/3 + (11*7**(0.5D0)*Q2(24,IELEM**3+1))/27 + 
     >   (15**(0.5D0)*Q2(19,IELEM**3+1))/9 + (77*3**(0.5D0)*Q2(32,IELE
     >   M**3+1))/81 + (21**(0.5D0)*Q2(14,IELEM**3+1))/3 - (7**(0.5D0)
     >   *Q2(30,IELEM**3+1))/3 - (5**(0.5D0)*Q2(33,IELEM**3+1))/3 - (1
     >   5**(0.5D0)*Q2(25,IELEM**3+1))/3 + (11*21**(0.5D0)*Q2(20,IELEM
     >   **3+1))/81 + (5*3**(0.5D0)*Q2(39,IELEM**3+1))/9 - (5**(0.5D0)
     >   *Q2(38,IELEM**3+1))/3 - (5*3**(0.5D0)*Q2(42,IELEM**3+1))/9 - 
     >   (11*35**(0.5D0)*Q2(12,IELEM**3+1))/27 + (5*5**(0.5D0)*Q2(43,I
     >   ELEM**3+1))/9 - (15**(0.5D0)*Q2(34,IELEM**3+1))/9 - (21**(0.5
     >   D0)*Q2(29,IELEM**3+1))/3 - (35**(0.5D0)*Q2(15,IELEM**3+1))/3 
     >   + (55*7**(0.5D0)*Q2(44,IELEM**3+1))/81 - (15**(0.5D0)*Q2(37,I
     >   ELEM**3+1))/3 + (77*5**(0.5D0)*Q2(48,IELEM**3+1))/81 + (5*7**
     >   (0.5D0)*Q2(47,IELEM**3+1))/9 + (11*7**(0.5D0)*Q2(49,IELEM**3+
     >   1))/27 - (847*3**(0.5D0)*Q2(56,IELEM**3+1))/729 + (11*7**(0.5
     >   D0)*Q2(54,IELEM**3+1))/27 + (77*3**(0.5D0)*Q2(62,IELEM**3+1))
     >   /81 - (847*5**(0.5D0)*Q2(60,IELEM**3+1))/729 - (55*7**(0.5D0)
     >   *Q2(59,IELEM**3+1))/81 - (77*5**(0.5D0)*Q2(63,IELEM**3+1))/81
     >    - (847*7**(0.5D0)*Q2(64,IELEM**3+1))/729 + (11*21**(0.5D0)*Q
     >   2(50,IELEM**3+1))/81 + (11*35**(0.5D0)*Q2(36,IELEM**3+1))/81 
     >   + (11*21**(0.5D0)*Q2(53,IELEM**3+1))/27 - (35**(0.5D0)*Q2(45,
     >   IELEM**3+1))/3 - (11*35**(0.5D0)*Q2(51,IELEM**3+1))/81 + (11*
     >   35**(0.5D0)*Q2(57,IELEM**3+1))/27 + (11*105**(0.5D0)*Q2(28,IE
     >   LEM**3+1))/81 + (105**(0.5D0)*Q2(31,IELEM**3+1))/9 + (11*105*
     >   *(0.5D0)*Q2(40,IELEM**3+1))/81 - (105**(0.5D0)*Q2(46,IELEM**3
     >   +1))/9 - (11*105**(0.5D0)*Q2(55,IELEM**3+1))/81 + (11*105**(0
     >   .5D0)*Q2(58,IELEM**3+1))/81) 
      CORNERQ(16,2) = (Q2(01,IELEM**3+1) + 3*Q2(06,IELEM**3+1) + 5*
     >   Q2(11,IELEM**3+1) + 7*Q2(16,IELEM**3+1) - Q2(18,IELEM**3+1) -
     >    Q2(21,IELEM**3+1) - (5*Q2(35,IELEM**3+1))/3 - (5*Q2(41,IELEM
     >   **3+1))/3 + (77*Q2(52,IELEM**3+1))/27 + (77*Q2(61,IELEM**3+1)
     >   )/27 + 3**(0.5D0)*Q2(02,IELEM**3+1) + 3**(0.5D0)*Q2(05,IELEM*
     >   *3+1) + 5**(0.5D0)*Q2(03,IELEM**3+1) + 7**(0.5D0)*Q2(04,IELEM
     >   **3+1) + 5**(0.5D0)*Q2(09,IELEM**3+1) - (3**(0.5D0)*Q2(17,IEL
     >   EM**3+1))/3 + 7**(0.5D0)*Q2(13,IELEM**3+1) + 15**(0.5D0)*Q2(0
     >   7,IELEM**3+1) - 3**(0.5D0)*Q2(22,IELEM**3+1) + 15**(0.5D0)*Q2
     >   (10,IELEM**3+1) - 5**(0.5D0)*Q2(23,IELEM**3+1) + 21**(0.5D0)*
     >   Q2(08,IELEM**3+1) - (5*3**(0.5D0)*Q2(27,IELEM**3+1))/3 - 5**(
     >   0.5D0)*Q2(26,IELEM**3+1) - 7**(0.5D0)*Q2(24,IELEM**3+1) - (15
     >   **(0.5D0)*Q2(19,IELEM**3+1))/3 - (7*3**(0.5D0)*Q2(32,IELEM**3
     >   +1))/3 + 21**(0.5D0)*Q2(14,IELEM**3+1) - 7**(0.5D0)*Q2(30,IEL
     >   EM**3+1) - (5**(0.5D0)*Q2(33,IELEM**3+1))/3 - (15**(0.5D0)*Q2
     >   (25,IELEM**3+1))/3 - (21**(0.5D0)*Q2(20,IELEM**3+1))/3 - (5*3
     >   **(0.5D0)*Q2(39,IELEM**3+1))/3 - 5**(0.5D0)*Q2(38,IELEM**3+1)
     >    - (5*3**(0.5D0)*Q2(42,IELEM**3+1))/3 + 35**(0.5D0)*Q2(12,IEL
     >   EM**3+1) - (5*5**(0.5D0)*Q2(43,IELEM**3+1))/3 - (15**(0.5D0)*
     >   Q2(34,IELEM**3+1))/3 - (21**(0.5D0)*Q2(29,IELEM**3+1))/3 + 35
     >   **(0.5D0)*Q2(15,IELEM**3+1) - (5*7**(0.5D0)*Q2(44,IELEM**3+1)
     >   )/3 - (15**(0.5D0)*Q2(37,IELEM**3+1))/3 - (7*5**(0.5D0)*Q2(48
     >   ,IELEM**3+1))/3 - (5*7**(0.5D0)*Q2(47,IELEM**3+1))/3 + (11*7*
     >   *(0.5D0)*Q2(49,IELEM**3+1))/27 + (77*3**(0.5D0)*Q2(56,IELEM**
     >   3+1))/27 + (11*7**(0.5D0)*Q2(54,IELEM**3+1))/9 + (77*3**(0.5D
     >   0)*Q2(62,IELEM**3+1))/27 + (77*5**(0.5D0)*Q2(60,IELEM**3+1))/
     >   27 + (55*7**(0.5D0)*Q2(59,IELEM**3+1))/27 + (77*5**(0.5D0)*Q2
     >   (63,IELEM**3+1))/27 + (77*7**(0.5D0)*Q2(64,IELEM**3+1))/27 + 
     >   (11*21**(0.5D0)*Q2(50,IELEM**3+1))/27 - (35**(0.5D0)*Q2(36,IE
     >   LEM**3+1))/3 + (11*21**(0.5D0)*Q2(53,IELEM**3+1))/27 - (35**(
     >   0.5D0)*Q2(45,IELEM**3+1))/3 + (11*35**(0.5D0)*Q2(51,IELEM**3+
     >   1))/27 + (11*35**(0.5D0)*Q2(57,IELEM**3+1))/27 - (105**(0.5D0
     >   )*Q2(28,IELEM**3+1))/3 - (105**(0.5D0)*Q2(31,IELEM**3+1))/3 -
     >    (105**(0.5D0)*Q2(40,IELEM**3+1))/3 - (105**(0.5D0)*Q2(46,IEL
     >   EM**3+1))/3 + (11*105**(0.5D0)*Q2(55,IELEM**3+1))/27 + (11*10
     >   5**(0.5D0)*Q2(58,IELEM**3+1))/27) 
      CORNERQ(01,3) = (Q2(01,IELEM**3+1) + 3*Q2(06,IELEM**3+1) + 5*
     >   Q2(11,IELEM**3+1) + 7*Q2(16,IELEM**3+1) - Q2(18,IELEM**3+1) -
     >    Q2(21,IELEM**3+1) - (5*Q2(35,IELEM**3+1))/3 - (5*Q2(41,IELEM
     >   **3+1))/3 + (77*Q2(52,IELEM**3+1))/27 + (77*Q2(61,IELEM**3+1)
     >   )/27 - 3**(0.5D0)*Q2(02,IELEM**3+1) - 3**(0.5D0)*Q2(05,IELEM*
     >   *3+1) + 5**(0.5D0)*Q2(03,IELEM**3+1) - 7**(0.5D0)*Q2(04,IELEM
     >   **3+1) + 5**(0.5D0)*Q2(09,IELEM**3+1) + (3**(0.5D0)*Q2(17,IEL
     >   EM**3+1))/3 - 7**(0.5D0)*Q2(13,IELEM**3+1) - 15**(0.5D0)*Q2(0
     >   7,IELEM**3+1) + 3**(0.5D0)*Q2(22,IELEM**3+1) - 15**(0.5D0)*Q2
     >   (10,IELEM**3+1) - 5**(0.5D0)*Q2(23,IELEM**3+1) + 21**(0.5D0)*
     >   Q2(08,IELEM**3+1) + (5*3**(0.5D0)*Q2(27,IELEM**3+1))/3 - 5**(
     >   0.5D0)*Q2(26,IELEM**3+1) + 7**(0.5D0)*Q2(24,IELEM**3+1) + (15
     >   **(0.5D0)*Q2(19,IELEM**3+1))/3 + (7*3**(0.5D0)*Q2(32,IELEM**3
     >   +1))/3 + 21**(0.5D0)*Q2(14,IELEM**3+1) + 7**(0.5D0)*Q2(30,IEL
     >   EM**3+1) - (5**(0.5D0)*Q2(33,IELEM**3+1))/3 + (15**(0.5D0)*Q2
     >   (25,IELEM**3+1))/3 - (21**(0.5D0)*Q2(20,IELEM**3+1))/3 + (5*3
     >   **(0.5D0)*Q2(39,IELEM**3+1))/3 - 5**(0.5D0)*Q2(38,IELEM**3+1)
     >    + (5*3**(0.5D0)*Q2(42,IELEM**3+1))/3 - 35**(0.5D0)*Q2(12,IEL
     >   EM**3+1) - (5*5**(0.5D0)*Q2(43,IELEM**3+1))/3 + (15**(0.5D0)*
     >   Q2(34,IELEM**3+1))/3 - (21**(0.5D0)*Q2(29,IELEM**3+1))/3 - 35
     >   **(0.5D0)*Q2(15,IELEM**3+1) + (5*7**(0.5D0)*Q2(44,IELEM**3+1)
     >   )/3 + (15**(0.5D0)*Q2(37,IELEM**3+1))/3 - (7*5**(0.5D0)*Q2(48
     >   ,IELEM**3+1))/3 + (5*7**(0.5D0)*Q2(47,IELEM**3+1))/3 - (11*7*
     >   *(0.5D0)*Q2(49,IELEM**3+1))/27 - (77*3**(0.5D0)*Q2(56,IELEM**
     >   3+1))/27 - (11*7**(0.5D0)*Q2(54,IELEM**3+1))/9 - (77*3**(0.5D
     >   0)*Q2(62,IELEM**3+1))/27 + (77*5**(0.5D0)*Q2(60,IELEM**3+1))/
     >   27 - (55*7**(0.5D0)*Q2(59,IELEM**3+1))/27 + (77*5**(0.5D0)*Q2
     >   (63,IELEM**3+1))/27 - (77*7**(0.5D0)*Q2(64,IELEM**3+1))/27 + 
     >   (11*21**(0.5D0)*Q2(50,IELEM**3+1))/27 + (35**(0.5D0)*Q2(36,IE
     >   LEM**3+1))/3 + (11*21**(0.5D0)*Q2(53,IELEM**3+1))/27 + (35**(
     >   0.5D0)*Q2(45,IELEM**3+1))/3 - (11*35**(0.5D0)*Q2(51,IELEM**3+
     >   1))/27 - (11*35**(0.5D0)*Q2(57,IELEM**3+1))/27 - (105**(0.5D0
     >   )*Q2(28,IELEM**3+1))/3 - (105**(0.5D0)*Q2(31,IELEM**3+1))/3 -
     >    (105**(0.5D0)*Q2(40,IELEM**3+1))/3 - (105**(0.5D0)*Q2(46,IEL
     >   EM**3+1))/3 + (11*105**(0.5D0)*Q2(55,IELEM**3+1))/27 + (11*10
     >   5**(0.5D0)*Q2(58,IELEM**3+1))/27) 
      CORNERQ(02,3) = (Q2(01,IELEM**3+1) + Q2(06,IELEM**3+1) - (5*Q
     >   2(11,IELEM**3+1))/3 - (77*Q2(16,IELEM**3+1))/27 - Q2(18,IELEM
     >   **3+1)/3 - Q2(21,IELEM**3+1) + (5*Q2(35,IELEM**3+1))/9 - (5*Q
     >   2(41,IELEM**3+1))/3 - (847*Q2(52,IELEM**3+1))/729 + (77*Q2(61
     >   ,IELEM**3+1))/27 - (3**(0.5D0)*Q2(02,IELEM**3+1))/3 - 3**(0.5
     >   D0)*Q2(05,IELEM**3+1) - (5**(0.5D0)*Q2(03,IELEM**3+1))/3 + (1
     >   1*7**(0.5D0)*Q2(04,IELEM**3+1))/27 + 5**(0.5D0)*Q2(09,IELEM**
     >   3+1) + (3**(0.5D0)*Q2(17,IELEM**3+1))/3 - 7**(0.5D0)*Q2(13,IE
     >   LEM**3+1) + (15**(0.5D0)*Q2(07,IELEM**3+1))/3 + (3**(0.5D0)*Q
     >   2(22,IELEM**3+1))/3 - (15**(0.5D0)*Q2(10,IELEM**3+1))/3 + (5*
     >   *(0.5D0)*Q2(23,IELEM**3+1))/3 - (11*21**(0.5D0)*Q2(08,IELEM**
     >   3+1))/27 - (5*3**(0.5D0)*Q2(27,IELEM**3+1))/9 - (5**(0.5D0)*Q
     >   2(26,IELEM**3+1))/3 - (11*7**(0.5D0)*Q2(24,IELEM**3+1))/27 - 
     >   (15**(0.5D0)*Q2(19,IELEM**3+1))/9 - (77*3**(0.5D0)*Q2(32,IELE
     >   M**3+1))/81 + (21**(0.5D0)*Q2(14,IELEM**3+1))/3 + (7**(0.5D0)
     >   *Q2(30,IELEM**3+1))/3 - (5**(0.5D0)*Q2(33,IELEM**3+1))/3 + (1
     >   5**(0.5D0)*Q2(25,IELEM**3+1))/3 + (11*21**(0.5D0)*Q2(20,IELEM
     >   **3+1))/81 - (5*3**(0.5D0)*Q2(39,IELEM**3+1))/9 - (5**(0.5D0)
     >   *Q2(38,IELEM**3+1))/3 + (5*3**(0.5D0)*Q2(42,IELEM**3+1))/9 + 
     >   (11*35**(0.5D0)*Q2(12,IELEM**3+1))/27 + (5*5**(0.5D0)*Q2(43,I
     >   ELEM**3+1))/9 + (15**(0.5D0)*Q2(34,IELEM**3+1))/9 - (21**(0.5
     >   D0)*Q2(29,IELEM**3+1))/3 + (35**(0.5D0)*Q2(15,IELEM**3+1))/3 
     >   - (55*7**(0.5D0)*Q2(44,IELEM**3+1))/81 + (15**(0.5D0)*Q2(37,I
     >   ELEM**3+1))/3 + (77*5**(0.5D0)*Q2(48,IELEM**3+1))/81 - (5*7**
     >   (0.5D0)*Q2(47,IELEM**3+1))/9 - (11*7**(0.5D0)*Q2(49,IELEM**3+
     >   1))/27 + (847*3**(0.5D0)*Q2(56,IELEM**3+1))/729 - (11*7**(0.5
     >   D0)*Q2(54,IELEM**3+1))/27 - (77*3**(0.5D0)*Q2(62,IELEM**3+1))
     >   /81 - (847*5**(0.5D0)*Q2(60,IELEM**3+1))/729 + (55*7**(0.5D0)
     >   *Q2(59,IELEM**3+1))/81 - (77*5**(0.5D0)*Q2(63,IELEM**3+1))/81
     >    + (847*7**(0.5D0)*Q2(64,IELEM**3+1))/729 + (11*21**(0.5D0)*Q
     >   2(50,IELEM**3+1))/81 - (11*35**(0.5D0)*Q2(36,IELEM**3+1))/81 
     >   + (11*21**(0.5D0)*Q2(53,IELEM**3+1))/27 + (35**(0.5D0)*Q2(45,
     >   IELEM**3+1))/3 + (11*35**(0.5D0)*Q2(51,IELEM**3+1))/81 - (11*
     >   35**(0.5D0)*Q2(57,IELEM**3+1))/27 + (11*105**(0.5D0)*Q2(28,IE
     >   LEM**3+1))/81 + (105**(0.5D0)*Q2(31,IELEM**3+1))/9 + (11*105*
     >   *(0.5D0)*Q2(40,IELEM**3+1))/81 - (105**(0.5D0)*Q2(46,IELEM**3
     >   +1))/9 - (11*105**(0.5D0)*Q2(55,IELEM**3+1))/81 + (11*105**(0
     >   .5D0)*Q2(58,IELEM**3+1))/81) 
      CORNERQ(03,3) = (Q2(01,IELEM**3+1) - Q2(06,IELEM**3+1) - (5*Q
     >   2(11,IELEM**3+1))/3 + (77*Q2(16,IELEM**3+1))/27 + Q2(18,IELEM
     >   **3+1)/3 - Q2(21,IELEM**3+1) + (5*Q2(35,IELEM**3+1))/9 - (5*Q
     >   2(41,IELEM**3+1))/3 + (847*Q2(52,IELEM**3+1))/729 + (77*Q2(61
     >   ,IELEM**3+1))/27 + (3**(0.5D0)*Q2(02,IELEM**3+1))/3 - 3**(0.5
     >   D0)*Q2(05,IELEM**3+1) - (5**(0.5D0)*Q2(03,IELEM**3+1))/3 - (1
     >   1*7**(0.5D0)*Q2(04,IELEM**3+1))/27 + 5**(0.5D0)*Q2(09,IELEM**
     >   3+1) + (3**(0.5D0)*Q2(17,IELEM**3+1))/3 - 7**(0.5D0)*Q2(13,IE
     >   LEM**3+1) + (15**(0.5D0)*Q2(07,IELEM**3+1))/3 - (3**(0.5D0)*Q
     >   2(22,IELEM**3+1))/3 + (15**(0.5D0)*Q2(10,IELEM**3+1))/3 + (5*
     >   *(0.5D0)*Q2(23,IELEM**3+1))/3 + (11*21**(0.5D0)*Q2(08,IELEM**
     >   3+1))/27 - (5*3**(0.5D0)*Q2(27,IELEM**3+1))/9 + (5**(0.5D0)*Q
     >   2(26,IELEM**3+1))/3 + (11*7**(0.5D0)*Q2(24,IELEM**3+1))/27 - 
     >   (15**(0.5D0)*Q2(19,IELEM**3+1))/9 + (77*3**(0.5D0)*Q2(32,IELE
     >   M**3+1))/81 - (21**(0.5D0)*Q2(14,IELEM**3+1))/3 - (7**(0.5D0)
     >   *Q2(30,IELEM**3+1))/3 - (5**(0.5D0)*Q2(33,IELEM**3+1))/3 + (1
     >   5**(0.5D0)*Q2(25,IELEM**3+1))/3 - (11*21**(0.5D0)*Q2(20,IELEM
     >   **3+1))/81 - (5*3**(0.5D0)*Q2(39,IELEM**3+1))/9 + (5**(0.5D0)
     >   *Q2(38,IELEM**3+1))/3 - (5*3**(0.5D0)*Q2(42,IELEM**3+1))/9 - 
     >   (11*35**(0.5D0)*Q2(12,IELEM**3+1))/27 + (5*5**(0.5D0)*Q2(43,I
     >   ELEM**3+1))/9 - (15**(0.5D0)*Q2(34,IELEM**3+1))/9 - (21**(0.5
     >   D0)*Q2(29,IELEM**3+1))/3 + (35**(0.5D0)*Q2(15,IELEM**3+1))/3 
     >   + (55*7**(0.5D0)*Q2(44,IELEM**3+1))/81 + (15**(0.5D0)*Q2(37,I
     >   ELEM**3+1))/3 - (77*5**(0.5D0)*Q2(48,IELEM**3+1))/81 - (5*7**
     >   (0.5D0)*Q2(47,IELEM**3+1))/9 - (11*7**(0.5D0)*Q2(49,IELEM**3+
     >   1))/27 - (847*3**(0.5D0)*Q2(56,IELEM**3+1))/729 + (11*7**(0.5
     >   D0)*Q2(54,IELEM**3+1))/27 + (77*3**(0.5D0)*Q2(62,IELEM**3+1))
     >   /81 + (847*5**(0.5D0)*Q2(60,IELEM**3+1))/729 + (55*7**(0.5D0)
     >   *Q2(59,IELEM**3+1))/81 - (77*5**(0.5D0)*Q2(63,IELEM**3+1))/81
     >    - (847*7**(0.5D0)*Q2(64,IELEM**3+1))/729 - (11*21**(0.5D0)*Q
     >   2(50,IELEM**3+1))/81 + (11*35**(0.5D0)*Q2(36,IELEM**3+1))/81 
     >   + (11*21**(0.5D0)*Q2(53,IELEM**3+1))/27 + (35**(0.5D0)*Q2(45,
     >   IELEM**3+1))/3 + (11*35**(0.5D0)*Q2(51,IELEM**3+1))/81 - (11*
     >   35**(0.5D0)*Q2(57,IELEM**3+1))/27 - (11*105**(0.5D0)*Q2(28,IE
     >   LEM**3+1))/81 + (105**(0.5D0)*Q2(31,IELEM**3+1))/9 - (11*105*
     >   *(0.5D0)*Q2(40,IELEM**3+1))/81 + (105**(0.5D0)*Q2(46,IELEM**3
     >   +1))/9 - (11*105**(0.5D0)*Q2(55,IELEM**3+1))/81 - (11*105**(0
     >   .5D0)*Q2(58,IELEM**3+1))/81) 
      CORNERQ(04,3) = (Q2(01,IELEM**3+1) - 3*Q2(06,IELEM**3+1) + 5*
     >   Q2(11,IELEM**3+1) - 7*Q2(16,IELEM**3+1) + Q2(18,IELEM**3+1) -
     >    Q2(21,IELEM**3+1) - (5*Q2(35,IELEM**3+1))/3 - (5*Q2(41,IELEM
     >   **3+1))/3 - (77*Q2(52,IELEM**3+1))/27 + (77*Q2(61,IELEM**3+1)
     >   )/27 + 3**(0.5D0)*Q2(02,IELEM**3+1) - 3**(0.5D0)*Q2(05,IELEM*
     >   *3+1) + 5**(0.5D0)*Q2(03,IELEM**3+1) + 7**(0.5D0)*Q2(04,IELEM
     >   **3+1) + 5**(0.5D0)*Q2(09,IELEM**3+1) + (3**(0.5D0)*Q2(17,IEL
     >   EM**3+1))/3 - 7**(0.5D0)*Q2(13,IELEM**3+1) - 15**(0.5D0)*Q2(0
     >   7,IELEM**3+1) - 3**(0.5D0)*Q2(22,IELEM**3+1) + 15**(0.5D0)*Q2
     >   (10,IELEM**3+1) - 5**(0.5D0)*Q2(23,IELEM**3+1) - 21**(0.5D0)*
     >   Q2(08,IELEM**3+1) + (5*3**(0.5D0)*Q2(27,IELEM**3+1))/3 + 5**(
     >   0.5D0)*Q2(26,IELEM**3+1) - 7**(0.5D0)*Q2(24,IELEM**3+1) + (15
     >   **(0.5D0)*Q2(19,IELEM**3+1))/3 - (7*3**(0.5D0)*Q2(32,IELEM**3
     >   +1))/3 - 21**(0.5D0)*Q2(14,IELEM**3+1) - 7**(0.5D0)*Q2(30,IEL
     >   EM**3+1) - (5**(0.5D0)*Q2(33,IELEM**3+1))/3 + (15**(0.5D0)*Q2
     >   (25,IELEM**3+1))/3 + (21**(0.5D0)*Q2(20,IELEM**3+1))/3 + (5*3
     >   **(0.5D0)*Q2(39,IELEM**3+1))/3 + 5**(0.5D0)*Q2(38,IELEM**3+1)
     >    - (5*3**(0.5D0)*Q2(42,IELEM**3+1))/3 + 35**(0.5D0)*Q2(12,IEL
     >   EM**3+1) - (5*5**(0.5D0)*Q2(43,IELEM**3+1))/3 - (15**(0.5D0)*
     >   Q2(34,IELEM**3+1))/3 - (21**(0.5D0)*Q2(29,IELEM**3+1))/3 - 35
     >   **(0.5D0)*Q2(15,IELEM**3+1) - (5*7**(0.5D0)*Q2(44,IELEM**3+1)
     >   )/3 + (15**(0.5D0)*Q2(37,IELEM**3+1))/3 + (7*5**(0.5D0)*Q2(48
     >   ,IELEM**3+1))/3 + (5*7**(0.5D0)*Q2(47,IELEM**3+1))/3 - (11*7*
     >   *(0.5D0)*Q2(49,IELEM**3+1))/27 + (77*3**(0.5D0)*Q2(56,IELEM**
     >   3+1))/27 + (11*7**(0.5D0)*Q2(54,IELEM**3+1))/9 + (77*3**(0.5D
     >   0)*Q2(62,IELEM**3+1))/27 - (77*5**(0.5D0)*Q2(60,IELEM**3+1))/
     >   27 - (55*7**(0.5D0)*Q2(59,IELEM**3+1))/27 + (77*5**(0.5D0)*Q2
     >   (63,IELEM**3+1))/27 + (77*7**(0.5D0)*Q2(64,IELEM**3+1))/27 - 
     >   (11*21**(0.5D0)*Q2(50,IELEM**3+1))/27 - (35**(0.5D0)*Q2(36,IE
     >   LEM**3+1))/3 + (11*21**(0.5D0)*Q2(53,IELEM**3+1))/27 + (35**(
     >   0.5D0)*Q2(45,IELEM**3+1))/3 - (11*35**(0.5D0)*Q2(51,IELEM**3+
     >   1))/27 - (11*35**(0.5D0)*Q2(57,IELEM**3+1))/27 + (105**(0.5D0
     >   )*Q2(28,IELEM**3+1))/3 - (105**(0.5D0)*Q2(31,IELEM**3+1))/3 +
     >    (105**(0.5D0)*Q2(40,IELEM**3+1))/3 + (105**(0.5D0)*Q2(46,IEL
     >   EM**3+1))/3 + (11*105**(0.5D0)*Q2(55,IELEM**3+1))/27 - (11*10
     >   5**(0.5D0)*Q2(58,IELEM**3+1))/27) 
      CORNERQ(05,3) = (Q2(01,IELEM**3+1) + Q2(06,IELEM**3+1) - (5*Q
     >   2(11,IELEM**3+1))/3 - (77*Q2(16,IELEM**3+1))/27 - Q2(18,IELEM
     >   **3+1) - Q2(21,IELEM**3+1)/3 - (5*Q2(35,IELEM**3+1))/3 + (5*Q
     >   2(41,IELEM**3+1))/9 + (77*Q2(52,IELEM**3+1))/27 - (847*Q2(61,
     >   IELEM**3+1))/729 - 3**(0.5D0)*Q2(02,IELEM**3+1) - (3**(0.5D0)
     >   *Q2(05,IELEM**3+1))/3 + 5**(0.5D0)*Q2(03,IELEM**3+1) - 7**(0.
     >   5D0)*Q2(04,IELEM**3+1) - (5**(0.5D0)*Q2(09,IELEM**3+1))/3 + (
     >   3**(0.5D0)*Q2(17,IELEM**3+1))/3 + (11*7**(0.5D0)*Q2(13,IELEM*
     >   *3+1))/27 - (15**(0.5D0)*Q2(07,IELEM**3+1))/3 + (3**(0.5D0)*Q
     >   2(22,IELEM**3+1))/3 + (15**(0.5D0)*Q2(10,IELEM**3+1))/3 - (5*
     >   *(0.5D0)*Q2(23,IELEM**3+1))/3 + (21**(0.5D0)*Q2(08,IELEM**3+1
     >   ))/3 - (5*3**(0.5D0)*Q2(27,IELEM**3+1))/9 + (5**(0.5D0)*Q2(26
     >   ,IELEM**3+1))/3 + (7**(0.5D0)*Q2(24,IELEM**3+1))/3 + (15**(0.
     >   5D0)*Q2(19,IELEM**3+1))/3 - (77*3**(0.5D0)*Q2(32,IELEM**3+1))
     >   /81 - (11*21**(0.5D0)*Q2(14,IELEM**3+1))/27 - (11*7**(0.5D0)*
     >   Q2(30,IELEM**3+1))/27 - (5**(0.5D0)*Q2(33,IELEM**3+1))/3 - (1
     >   5**(0.5D0)*Q2(25,IELEM**3+1))/9 - (21**(0.5D0)*Q2(20,IELEM**3
     >   +1))/3 + (5*3**(0.5D0)*Q2(39,IELEM**3+1))/9 - (5**(0.5D0)*Q2(
     >   38,IELEM**3+1))/3 - (5*3**(0.5D0)*Q2(42,IELEM**3+1))/9 + (35*
     >   *(0.5D0)*Q2(12,IELEM**3+1))/3 + (5*5**(0.5D0)*Q2(43,IELEM**3+
     >   1))/9 + (15**(0.5D0)*Q2(34,IELEM**3+1))/3 + (11*21**(0.5D0)*Q
     >   2(29,IELEM**3+1))/81 + (11*35**(0.5D0)*Q2(15,IELEM**3+1))/27 
     >   - (5*7**(0.5D0)*Q2(44,IELEM**3+1))/9 + (15**(0.5D0)*Q2(37,IEL
     >   EM**3+1))/9 + (77*5**(0.5D0)*Q2(48,IELEM**3+1))/81 - (55*7**(
     >   0.5D0)*Q2(47,IELEM**3+1))/81 - (11*7**(0.5D0)*Q2(49,IELEM**3+
     >   1))/27 - (77*3**(0.5D0)*Q2(56,IELEM**3+1))/81 - (11*7**(0.5D0
     >   )*Q2(54,IELEM**3+1))/27 + (847*3**(0.5D0)*Q2(62,IELEM**3+1))/
     >   729 - (77*5**(0.5D0)*Q2(60,IELEM**3+1))/81 + (55*7**(0.5D0)*Q
     >   2(59,IELEM**3+1))/81 - (847*5**(0.5D0)*Q2(63,IELEM**3+1))/729
     >    + (847*7**(0.5D0)*Q2(64,IELEM**3+1))/729 + (11*21**(0.5D0)*Q
     >   2(50,IELEM**3+1))/27 + (35**(0.5D0)*Q2(36,IELEM**3+1))/3 + (1
     >   1*21**(0.5D0)*Q2(53,IELEM**3+1))/81 - (11*35**(0.5D0)*Q2(45,I
     >   ELEM**3+1))/81 - (11*35**(0.5D0)*Q2(51,IELEM**3+1))/27 + (11*
     >   35**(0.5D0)*Q2(57,IELEM**3+1))/81 + (105**(0.5D0)*Q2(28,IELEM
     >   **3+1))/9 + (11*105**(0.5D0)*Q2(31,IELEM**3+1))/81 - (105**(0
     >   .5D0)*Q2(40,IELEM**3+1))/9 + (11*105**(0.5D0)*Q2(46,IELEM**3+
     >   1))/81 + (11*105**(0.5D0)*Q2(55,IELEM**3+1))/81 - (11*105**(0
     >   .5D0)*Q2(58,IELEM**3+1))/81) 
      CORNERQ(06,3) = (Q2(01,IELEM**3+1) + Q2(06,IELEM**3+1)/3 + (5
     >   *Q2(11,IELEM**3+1))/9 + (847*Q2(16,IELEM**3+1))/729 - Q2(18,I
     >   ELEM**3+1)/3 - Q2(21,IELEM**3+1)/3 + (5*Q2(35,IELEM**3+1))/9 
     >   + (5*Q2(41,IELEM**3+1))/9 - (847*Q2(52,IELEM**3+1))/729 - (84
     >   7*Q2(61,IELEM**3+1))/729 - (3**(0.5D0)*Q2(02,IELEM**3+1))/3 -
     >    (3**(0.5D0)*Q2(05,IELEM**3+1))/3 - (5**(0.5D0)*Q2(03,IELEM**
     >   3+1))/3 + (11*7**(0.5D0)*Q2(04,IELEM**3+1))/27 - (5**(0.5D0)*
     >   Q2(09,IELEM**3+1))/3 + (3**(0.5D0)*Q2(17,IELEM**3+1))/3 + (11
     >   *7**(0.5D0)*Q2(13,IELEM**3+1))/27 + (15**(0.5D0)*Q2(07,IELEM*
     >   *3+1))/9 + (3**(0.5D0)*Q2(22,IELEM**3+1))/9 + (15**(0.5D0)*Q2
     >   (10,IELEM**3+1))/9 + (5**(0.5D0)*Q2(23,IELEM**3+1))/9 - (11*2
     >   1**(0.5D0)*Q2(08,IELEM**3+1))/81 + (5*3**(0.5D0)*Q2(27,IELEM*
     >   *3+1))/27 + (5**(0.5D0)*Q2(26,IELEM**3+1))/9 - (11*7**(0.5D0)
     >   *Q2(24,IELEM**3+1))/81 - (15**(0.5D0)*Q2(19,IELEM**3+1))/9 + 
     >   (847*3**(0.5D0)*Q2(32,IELEM**3+1))/2187 - (11*21**(0.5D0)*Q2(
     >   14,IELEM**3+1))/81 - (11*7**(0.5D0)*Q2(30,IELEM**3+1))/81 - (
     >   5**(0.5D0)*Q2(33,IELEM**3+1))/3 - (15**(0.5D0)*Q2(25,IELEM**3
     >   +1))/9 + (11*21**(0.5D0)*Q2(20,IELEM**3+1))/81 - (5*3**(0.5D0
     >   )*Q2(39,IELEM**3+1))/27 - (5**(0.5D0)*Q2(38,IELEM**3+1))/9 - 
     >   (5*3**(0.5D0)*Q2(42,IELEM**3+1))/27 - (11*35**(0.5D0)*Q2(12,I
     >   ELEM**3+1))/81 - (5*5**(0.5D0)*Q2(43,IELEM**3+1))/27 + (15**(
     >   0.5D0)*Q2(34,IELEM**3+1))/9 + (11*21**(0.5D0)*Q2(29,IELEM**3+
     >   1))/81 - (11*35**(0.5D0)*Q2(15,IELEM**3+1))/81 + (55*7**(0.5D
     >   0)*Q2(44,IELEM**3+1))/243 + (15**(0.5D0)*Q2(37,IELEM**3+1))/9
     >    - (847*5**(0.5D0)*Q2(48,IELEM**3+1))/2187 + (55*7**(0.5D0)*Q
     >   2(47,IELEM**3+1))/243 - (11*7**(0.5D0)*Q2(49,IELEM**3+1))/27 
     >   + (847*3**(0.5D0)*Q2(56,IELEM**3+1))/2187 - (11*7**(0.5D0)*Q2
     >   (54,IELEM**3+1))/81 + (847*3**(0.5D0)*Q2(62,IELEM**3+1))/2187
     >    + (847*5**(0.5D0)*Q2(60,IELEM**3+1))/2187 - (55*7**(0.5D0)*Q
     >   2(59,IELEM**3+1))/243 + (847*5**(0.5D0)*Q2(63,IELEM**3+1))/21
     >   87 - (9317*7**(0.5D0)*Q2(64,IELEM**3+1))/19683 + (11*21**(0.5
     >   D0)*Q2(50,IELEM**3+1))/81 - (11*35**(0.5D0)*Q2(36,IELEM**3+1)
     >   )/81 + (11*21**(0.5D0)*Q2(53,IELEM**3+1))/81 - (11*35**(0.5D0
     >   )*Q2(45,IELEM**3+1))/81 + (11*35**(0.5D0)*Q2(51,IELEM**3+1))/
     >   81 + (11*35**(0.5D0)*Q2(57,IELEM**3+1))/81 - (11*105**(0.5D0)
     >   *Q2(28,IELEM**3+1))/243 - (11*105**(0.5D0)*Q2(31,IELEM**3+1))
     >   /243 + (11*105**(0.5D0)*Q2(40,IELEM**3+1))/243 + (11*105**(0.
     >   5D0)*Q2(46,IELEM**3+1))/243 - (11*105**(0.5D0)*Q2(55,IELEM**3
     >   +1))/243 - (11*105**(0.5D0)*Q2(58,IELEM**3+1))/243) 
      CORNERQ(07,3) = (Q2(01,IELEM**3+1) - Q2(06,IELEM**3+1)/3 + (5
     >   *Q2(11,IELEM**3+1))/9 - (847*Q2(16,IELEM**3+1))/729 + Q2(18,I
     >   ELEM**3+1)/3 - Q2(21,IELEM**3+1)/3 + (5*Q2(35,IELEM**3+1))/9 
     >   + (5*Q2(41,IELEM**3+1))/9 + (847*Q2(52,IELEM**3+1))/729 - (84
     >   7*Q2(61,IELEM**3+1))/729 + (3**(0.5D0)*Q2(02,IELEM**3+1))/3 -
     >    (3**(0.5D0)*Q2(05,IELEM**3+1))/3 - (5**(0.5D0)*Q2(03,IELEM**
     >   3+1))/3 - (11*7**(0.5D0)*Q2(04,IELEM**3+1))/27 - (5**(0.5D0)*
     >   Q2(09,IELEM**3+1))/3 + (3**(0.5D0)*Q2(17,IELEM**3+1))/3 + (11
     >   *7**(0.5D0)*Q2(13,IELEM**3+1))/27 + (15**(0.5D0)*Q2(07,IELEM*
     >   *3+1))/9 - (3**(0.5D0)*Q2(22,IELEM**3+1))/9 - (15**(0.5D0)*Q2
     >   (10,IELEM**3+1))/9 + (5**(0.5D0)*Q2(23,IELEM**3+1))/9 + (11*2
     >   1**(0.5D0)*Q2(08,IELEM**3+1))/81 + (5*3**(0.5D0)*Q2(27,IELEM*
     >   *3+1))/27 - (5**(0.5D0)*Q2(26,IELEM**3+1))/9 + (11*7**(0.5D0)
     >   *Q2(24,IELEM**3+1))/81 - (15**(0.5D0)*Q2(19,IELEM**3+1))/9 - 
     >   (847*3**(0.5D0)*Q2(32,IELEM**3+1))/2187 + (11*21**(0.5D0)*Q2(
     >   14,IELEM**3+1))/81 + (11*7**(0.5D0)*Q2(30,IELEM**3+1))/81 - (
     >   5**(0.5D0)*Q2(33,IELEM**3+1))/3 - (15**(0.5D0)*Q2(25,IELEM**3
     >   +1))/9 - (11*21**(0.5D0)*Q2(20,IELEM**3+1))/81 - (5*3**(0.5D0
     >   )*Q2(39,IELEM**3+1))/27 + (5**(0.5D0)*Q2(38,IELEM**3+1))/9 + 
     >   (5*3**(0.5D0)*Q2(42,IELEM**3+1))/27 + (11*35**(0.5D0)*Q2(12,I
     >   ELEM**3+1))/81 - (5*5**(0.5D0)*Q2(43,IELEM**3+1))/27 - (15**(
     >   0.5D0)*Q2(34,IELEM**3+1))/9 + (11*21**(0.5D0)*Q2(29,IELEM**3+
     >   1))/81 - (11*35**(0.5D0)*Q2(15,IELEM**3+1))/81 - (55*7**(0.5D
     >   0)*Q2(44,IELEM**3+1))/243 + (15**(0.5D0)*Q2(37,IELEM**3+1))/9
     >    + (847*5**(0.5D0)*Q2(48,IELEM**3+1))/2187 + (55*7**(0.5D0)*Q
     >   2(47,IELEM**3+1))/243 - (11*7**(0.5D0)*Q2(49,IELEM**3+1))/27 
     >   - (847*3**(0.5D0)*Q2(56,IELEM**3+1))/2187 + (11*7**(0.5D0)*Q2
     >   (54,IELEM**3+1))/81 - (847*3**(0.5D0)*Q2(62,IELEM**3+1))/2187
     >    - (847*5**(0.5D0)*Q2(60,IELEM**3+1))/2187 - (55*7**(0.5D0)*Q
     >   2(59,IELEM**3+1))/243 + (847*5**(0.5D0)*Q2(63,IELEM**3+1))/21
     >   87 + (9317*7**(0.5D0)*Q2(64,IELEM**3+1))/19683 - (11*21**(0.5
     >   D0)*Q2(50,IELEM**3+1))/81 + (11*35**(0.5D0)*Q2(36,IELEM**3+1)
     >   )/81 + (11*21**(0.5D0)*Q2(53,IELEM**3+1))/81 - (11*35**(0.5D0
     >   )*Q2(45,IELEM**3+1))/81 + (11*35**(0.5D0)*Q2(51,IELEM**3+1))/
     >   81 + (11*35**(0.5D0)*Q2(57,IELEM**3+1))/81 + (11*105**(0.5D0)
     >   *Q2(28,IELEM**3+1))/243 - (11*105**(0.5D0)*Q2(31,IELEM**3+1))
     >   /243 - (11*105**(0.5D0)*Q2(40,IELEM**3+1))/243 - (11*105**(0.
     >   5D0)*Q2(46,IELEM**3+1))/243 - (11*105**(0.5D0)*Q2(55,IELEM**3
     >   +1))/243 + (11*105**(0.5D0)*Q2(58,IELEM**3+1))/243) 
      CORNERQ(08,3) = (Q2(01,IELEM**3+1) - Q2(06,IELEM**3+1) - (5*Q
     >   2(11,IELEM**3+1))/3 + (77*Q2(16,IELEM**3+1))/27 + Q2(18,IELEM
     >   **3+1) - Q2(21,IELEM**3+1)/3 - (5*Q2(35,IELEM**3+1))/3 + (5*Q
     >   2(41,IELEM**3+1))/9 - (77*Q2(52,IELEM**3+1))/27 - (847*Q2(61,
     >   IELEM**3+1))/729 + 3**(0.5D0)*Q2(02,IELEM**3+1) - (3**(0.5D0)
     >   *Q2(05,IELEM**3+1))/3 + 5**(0.5D0)*Q2(03,IELEM**3+1) + 7**(0.
     >   5D0)*Q2(04,IELEM**3+1) - (5**(0.5D0)*Q2(09,IELEM**3+1))/3 + (
     >   3**(0.5D0)*Q2(17,IELEM**3+1))/3 + (11*7**(0.5D0)*Q2(13,IELEM*
     >   *3+1))/27 - (15**(0.5D0)*Q2(07,IELEM**3+1))/3 - (3**(0.5D0)*Q
     >   2(22,IELEM**3+1))/3 - (15**(0.5D0)*Q2(10,IELEM**3+1))/3 - (5*
     >   *(0.5D0)*Q2(23,IELEM**3+1))/3 - (21**(0.5D0)*Q2(08,IELEM**3+1
     >   ))/3 - (5*3**(0.5D0)*Q2(27,IELEM**3+1))/9 - (5**(0.5D0)*Q2(26
     >   ,IELEM**3+1))/3 - (7**(0.5D0)*Q2(24,IELEM**3+1))/3 + (15**(0.
     >   5D0)*Q2(19,IELEM**3+1))/3 + (77*3**(0.5D0)*Q2(32,IELEM**3+1))
     >   /81 + (11*21**(0.5D0)*Q2(14,IELEM**3+1))/27 + (11*7**(0.5D0)*
     >   Q2(30,IELEM**3+1))/27 - (5**(0.5D0)*Q2(33,IELEM**3+1))/3 - (1
     >   5**(0.5D0)*Q2(25,IELEM**3+1))/9 + (21**(0.5D0)*Q2(20,IELEM**3
     >   +1))/3 + (5*3**(0.5D0)*Q2(39,IELEM**3+1))/9 + (5**(0.5D0)*Q2(
     >   38,IELEM**3+1))/3 + (5*3**(0.5D0)*Q2(42,IELEM**3+1))/9 - (35*
     >   *(0.5D0)*Q2(12,IELEM**3+1))/3 + (5*5**(0.5D0)*Q2(43,IELEM**3+
     >   1))/9 - (15**(0.5D0)*Q2(34,IELEM**3+1))/3 + (11*21**(0.5D0)*Q
     >   2(29,IELEM**3+1))/81 + (11*35**(0.5D0)*Q2(15,IELEM**3+1))/27 
     >   + (5*7**(0.5D0)*Q2(44,IELEM**3+1))/9 + (15**(0.5D0)*Q2(37,IEL
     >   EM**3+1))/9 - (77*5**(0.5D0)*Q2(48,IELEM**3+1))/81 - (55*7**(
     >   0.5D0)*Q2(47,IELEM**3+1))/81 - (11*7**(0.5D0)*Q2(49,IELEM**3+
     >   1))/27 + (77*3**(0.5D0)*Q2(56,IELEM**3+1))/81 + (11*7**(0.5D0
     >   )*Q2(54,IELEM**3+1))/27 - (847*3**(0.5D0)*Q2(62,IELEM**3+1))/
     >   729 + (77*5**(0.5D0)*Q2(60,IELEM**3+1))/81 + (55*7**(0.5D0)*Q
     >   2(59,IELEM**3+1))/81 - (847*5**(0.5D0)*Q2(63,IELEM**3+1))/729
     >    - (847*7**(0.5D0)*Q2(64,IELEM**3+1))/729 - (11*21**(0.5D0)*Q
     >   2(50,IELEM**3+1))/27 - (35**(0.5D0)*Q2(36,IELEM**3+1))/3 + (1
     >   1*21**(0.5D0)*Q2(53,IELEM**3+1))/81 - (11*35**(0.5D0)*Q2(45,I
     >   ELEM**3+1))/81 - (11*35**(0.5D0)*Q2(51,IELEM**3+1))/27 + (11*
     >   35**(0.5D0)*Q2(57,IELEM**3+1))/81 - (105**(0.5D0)*Q2(28,IELEM
     >   **3+1))/9 + (11*105**(0.5D0)*Q2(31,IELEM**3+1))/81 + (105**(0
     >   .5D0)*Q2(40,IELEM**3+1))/9 - (11*105**(0.5D0)*Q2(46,IELEM**3+
     >   1))/81 + (11*105**(0.5D0)*Q2(55,IELEM**3+1))/81 + (11*105**(0
     >   .5D0)*Q2(58,IELEM**3+1))/81) 
      CORNERQ(09,3) = (Q2(01,IELEM**3+1) - Q2(06,IELEM**3+1) - (5*Q
     >   2(11,IELEM**3+1))/3 + (77*Q2(16,IELEM**3+1))/27 - Q2(18,IELEM
     >   **3+1) + Q2(21,IELEM**3+1)/3 - (5*Q2(35,IELEM**3+1))/3 + (5*Q
     >   2(41,IELEM**3+1))/9 + (77*Q2(52,IELEM**3+1))/27 + (847*Q2(61,
     >   IELEM**3+1))/729 - 3**(0.5D0)*Q2(02,IELEM**3+1) + (3**(0.5D0)
     >   *Q2(05,IELEM**3+1))/3 + 5**(0.5D0)*Q2(03,IELEM**3+1) - 7**(0.
     >   5D0)*Q2(04,IELEM**3+1) - (5**(0.5D0)*Q2(09,IELEM**3+1))/3 + (
     >   3**(0.5D0)*Q2(17,IELEM**3+1))/3 - (11*7**(0.5D0)*Q2(13,IELEM*
     >   *3+1))/27 + (15**(0.5D0)*Q2(07,IELEM**3+1))/3 - (3**(0.5D0)*Q
     >   2(22,IELEM**3+1))/3 + (15**(0.5D0)*Q2(10,IELEM**3+1))/3 + (5*
     >   *(0.5D0)*Q2(23,IELEM**3+1))/3 - (21**(0.5D0)*Q2(08,IELEM**3+1
     >   ))/3 - (5*3**(0.5D0)*Q2(27,IELEM**3+1))/9 + (5**(0.5D0)*Q2(26
     >   ,IELEM**3+1))/3 - (7**(0.5D0)*Q2(24,IELEM**3+1))/3 + (15**(0.
     >   5D0)*Q2(19,IELEM**3+1))/3 + (77*3**(0.5D0)*Q2(32,IELEM**3+1))
     >   /81 + (11*21**(0.5D0)*Q2(14,IELEM**3+1))/27 + (11*7**(0.5D0)*
     >   Q2(30,IELEM**3+1))/27 - (5**(0.5D0)*Q2(33,IELEM**3+1))/3 - (1
     >   5**(0.5D0)*Q2(25,IELEM**3+1))/9 - (21**(0.5D0)*Q2(20,IELEM**3
     >   +1))/3 - (5*3**(0.5D0)*Q2(39,IELEM**3+1))/9 + (5**(0.5D0)*Q2(
     >   38,IELEM**3+1))/3 - (5*3**(0.5D0)*Q2(42,IELEM**3+1))/9 + (35*
     >   *(0.5D0)*Q2(12,IELEM**3+1))/3 + (5*5**(0.5D0)*Q2(43,IELEM**3+
     >   1))/9 + (15**(0.5D0)*Q2(34,IELEM**3+1))/3 - (11*21**(0.5D0)*Q
     >   2(29,IELEM**3+1))/81 - (11*35**(0.5D0)*Q2(15,IELEM**3+1))/27 
     >   - (5*7**(0.5D0)*Q2(44,IELEM**3+1))/9 - (15**(0.5D0)*Q2(37,IEL
     >   EM**3+1))/9 - (77*5**(0.5D0)*Q2(48,IELEM**3+1))/81 + (55*7**(
     >   0.5D0)*Q2(47,IELEM**3+1))/81 - (11*7**(0.5D0)*Q2(49,IELEM**3+
     >   1))/27 + (77*3**(0.5D0)*Q2(56,IELEM**3+1))/81 + (11*7**(0.5D0
     >   )*Q2(54,IELEM**3+1))/27 - (847*3**(0.5D0)*Q2(62,IELEM**3+1))/
     >   729 - (77*5**(0.5D0)*Q2(60,IELEM**3+1))/81 + (55*7**(0.5D0)*Q
     >   2(59,IELEM**3+1))/81 + (847*5**(0.5D0)*Q2(63,IELEM**3+1))/729
     >    - (847*7**(0.5D0)*Q2(64,IELEM**3+1))/729 + (11*21**(0.5D0)*Q
     >   2(50,IELEM**3+1))/27 + (35**(0.5D0)*Q2(36,IELEM**3+1))/3 - (1
     >   1*21**(0.5D0)*Q2(53,IELEM**3+1))/81 + (11*35**(0.5D0)*Q2(45,I
     >   ELEM**3+1))/81 - (11*35**(0.5D0)*Q2(51,IELEM**3+1))/27 + (11*
     >   35**(0.5D0)*Q2(57,IELEM**3+1))/81 + (105**(0.5D0)*Q2(28,IELEM
     >   **3+1))/9 - (11*105**(0.5D0)*Q2(31,IELEM**3+1))/81 + (105**(0
     >   .5D0)*Q2(40,IELEM**3+1))/9 - (11*105**(0.5D0)*Q2(46,IELEM**3+
     >   1))/81 - (11*105**(0.5D0)*Q2(55,IELEM**3+1))/81 - (11*105**(0
     >   .5D0)*Q2(58,IELEM**3+1))/81) 
      CORNERQ(10,3) = (Q2(01,IELEM**3+1) - Q2(06,IELEM**3+1)/3 + (5
     >   *Q2(11,IELEM**3+1))/9 - (847*Q2(16,IELEM**3+1))/729 - Q2(18,I
     >   ELEM**3+1)/3 + Q2(21,IELEM**3+1)/3 + (5*Q2(35,IELEM**3+1))/9 
     >   + (5*Q2(41,IELEM**3+1))/9 - (847*Q2(52,IELEM**3+1))/729 + (84
     >   7*Q2(61,IELEM**3+1))/729 - (3**(0.5D0)*Q2(02,IELEM**3+1))/3 +
     >    (3**(0.5D0)*Q2(05,IELEM**3+1))/3 - (5**(0.5D0)*Q2(03,IELEM**
     >   3+1))/3 + (11*7**(0.5D0)*Q2(04,IELEM**3+1))/27 - (5**(0.5D0)*
     >   Q2(09,IELEM**3+1))/3 + (3**(0.5D0)*Q2(17,IELEM**3+1))/3 - (11
     >   *7**(0.5D0)*Q2(13,IELEM**3+1))/27 - (15**(0.5D0)*Q2(07,IELEM*
     >   *3+1))/9 - (3**(0.5D0)*Q2(22,IELEM**3+1))/9 + (15**(0.5D0)*Q2
     >   (10,IELEM**3+1))/9 - (5**(0.5D0)*Q2(23,IELEM**3+1))/9 + (11*2
     >   1**(0.5D0)*Q2(08,IELEM**3+1))/81 + (5*3**(0.5D0)*Q2(27,IELEM*
     >   *3+1))/27 + (5**(0.5D0)*Q2(26,IELEM**3+1))/9 + (11*7**(0.5D0)
     >   *Q2(24,IELEM**3+1))/81 - (15**(0.5D0)*Q2(19,IELEM**3+1))/9 - 
     >   (847*3**(0.5D0)*Q2(32,IELEM**3+1))/2187 + (11*21**(0.5D0)*Q2(
     >   14,IELEM**3+1))/81 + (11*7**(0.5D0)*Q2(30,IELEM**3+1))/81 - (
     >   5**(0.5D0)*Q2(33,IELEM**3+1))/3 - (15**(0.5D0)*Q2(25,IELEM**3
     >   +1))/9 + (11*21**(0.5D0)*Q2(20,IELEM**3+1))/81 + (5*3**(0.5D0
     >   )*Q2(39,IELEM**3+1))/27 + (5**(0.5D0)*Q2(38,IELEM**3+1))/9 - 
     >   (5*3**(0.5D0)*Q2(42,IELEM**3+1))/27 - (11*35**(0.5D0)*Q2(12,I
     >   ELEM**3+1))/81 - (5*5**(0.5D0)*Q2(43,IELEM**3+1))/27 + (15**(
     >   0.5D0)*Q2(34,IELEM**3+1))/9 - (11*21**(0.5D0)*Q2(29,IELEM**3+
     >   1))/81 + (11*35**(0.5D0)*Q2(15,IELEM**3+1))/81 + (55*7**(0.5D
     >   0)*Q2(44,IELEM**3+1))/243 - (15**(0.5D0)*Q2(37,IELEM**3+1))/9
     >    + (847*5**(0.5D0)*Q2(48,IELEM**3+1))/2187 - (55*7**(0.5D0)*Q
     >   2(47,IELEM**3+1))/243 - (11*7**(0.5D0)*Q2(49,IELEM**3+1))/27 
     >   - (847*3**(0.5D0)*Q2(56,IELEM**3+1))/2187 + (11*7**(0.5D0)*Q2
     >   (54,IELEM**3+1))/81 - (847*3**(0.5D0)*Q2(62,IELEM**3+1))/2187
     >    + (847*5**(0.5D0)*Q2(60,IELEM**3+1))/2187 - (55*7**(0.5D0)*Q
     >   2(59,IELEM**3+1))/243 - (847*5**(0.5D0)*Q2(63,IELEM**3+1))/21
     >   87 + (9317*7**(0.5D0)*Q2(64,IELEM**3+1))/19683 + (11*21**(0.5
     >   D0)*Q2(50,IELEM**3+1))/81 - (11*35**(0.5D0)*Q2(36,IELEM**3+1)
     >   )/81 - (11*21**(0.5D0)*Q2(53,IELEM**3+1))/81 + (11*35**(0.5D0
     >   )*Q2(45,IELEM**3+1))/81 + (11*35**(0.5D0)*Q2(51,IELEM**3+1))/
     >   81 + (11*35**(0.5D0)*Q2(57,IELEM**3+1))/81 - (11*105**(0.5D0)
     >   *Q2(28,IELEM**3+1))/243 + (11*105**(0.5D0)*Q2(31,IELEM**3+1))
     >   /243 - (11*105**(0.5D0)*Q2(40,IELEM**3+1))/243 - (11*105**(0.
     >   5D0)*Q2(46,IELEM**3+1))/243 + (11*105**(0.5D0)*Q2(55,IELEM**3
     >   +1))/243 - (11*105**(0.5D0)*Q2(58,IELEM**3+1))/243) 
      CORNERQ(11,3) = (Q2(01,IELEM**3+1) + Q2(06,IELEM**3+1)/3 + (5
     >   *Q2(11,IELEM**3+1))/9 + (847*Q2(16,IELEM**3+1))/729 + Q2(18,I
     >   ELEM**3+1)/3 + Q2(21,IELEM**3+1)/3 + (5*Q2(35,IELEM**3+1))/9 
     >   + (5*Q2(41,IELEM**3+1))/9 + (847*Q2(52,IELEM**3+1))/729 + (84
     >   7*Q2(61,IELEM**3+1))/729 + (3**(0.5D0)*Q2(02,IELEM**3+1))/3 +
     >    (3**(0.5D0)*Q2(05,IELEM**3+1))/3 - (5**(0.5D0)*Q2(03,IELEM**
     >   3+1))/3 - (11*7**(0.5D0)*Q2(04,IELEM**3+1))/27 - (5**(0.5D0)*
     >   Q2(09,IELEM**3+1))/3 + (3**(0.5D0)*Q2(17,IELEM**3+1))/3 - (11
     >   *7**(0.5D0)*Q2(13,IELEM**3+1))/27 - (15**(0.5D0)*Q2(07,IELEM*
     >   *3+1))/9 + (3**(0.5D0)*Q2(22,IELEM**3+1))/9 - (15**(0.5D0)*Q2
     >   (10,IELEM**3+1))/9 - (5**(0.5D0)*Q2(23,IELEM**3+1))/9 - (11*2
     >   1**(0.5D0)*Q2(08,IELEM**3+1))/81 + (5*3**(0.5D0)*Q2(27,IELEM*
     >   *3+1))/27 - (5**(0.5D0)*Q2(26,IELEM**3+1))/9 - (11*7**(0.5D0)
     >   *Q2(24,IELEM**3+1))/81 - (15**(0.5D0)*Q2(19,IELEM**3+1))/9 + 
     >   (847*3**(0.5D0)*Q2(32,IELEM**3+1))/2187 - (11*21**(0.5D0)*Q2(
     >   14,IELEM**3+1))/81 - (11*7**(0.5D0)*Q2(30,IELEM**3+1))/81 - (
     >   5**(0.5D0)*Q2(33,IELEM**3+1))/3 - (15**(0.5D0)*Q2(25,IELEM**3
     >   +1))/9 - (11*21**(0.5D0)*Q2(20,IELEM**3+1))/81 + (5*3**(0.5D0
     >   )*Q2(39,IELEM**3+1))/27 - (5**(0.5D0)*Q2(38,IELEM**3+1))/9 + 
     >   (5*3**(0.5D0)*Q2(42,IELEM**3+1))/27 + (11*35**(0.5D0)*Q2(12,I
     >   ELEM**3+1))/81 - (5*5**(0.5D0)*Q2(43,IELEM**3+1))/27 - (15**(
     >   0.5D0)*Q2(34,IELEM**3+1))/9 - (11*21**(0.5D0)*Q2(29,IELEM**3+
     >   1))/81 + (11*35**(0.5D0)*Q2(15,IELEM**3+1))/81 - (55*7**(0.5D
     >   0)*Q2(44,IELEM**3+1))/243 - (15**(0.5D0)*Q2(37,IELEM**3+1))/9
     >    - (847*5**(0.5D0)*Q2(48,IELEM**3+1))/2187 - (55*7**(0.5D0)*Q
     >   2(47,IELEM**3+1))/243 - (11*7**(0.5D0)*Q2(49,IELEM**3+1))/27 
     >   + (847*3**(0.5D0)*Q2(56,IELEM**3+1))/2187 - (11*7**(0.5D0)*Q2
     >   (54,IELEM**3+1))/81 + (847*3**(0.5D0)*Q2(62,IELEM**3+1))/2187
     >    - (847*5**(0.5D0)*Q2(60,IELEM**3+1))/2187 - (55*7**(0.5D0)*Q
     >   2(59,IELEM**3+1))/243 - (847*5**(0.5D0)*Q2(63,IELEM**3+1))/21
     >   87 - (9317*7**(0.5D0)*Q2(64,IELEM**3+1))/19683 - (11*21**(0.5
     >   D0)*Q2(50,IELEM**3+1))/81 + (11*35**(0.5D0)*Q2(36,IELEM**3+1)
     >   )/81 - (11*21**(0.5D0)*Q2(53,IELEM**3+1))/81 + (11*35**(0.5D0
     >   )*Q2(45,IELEM**3+1))/81 + (11*35**(0.5D0)*Q2(51,IELEM**3+1))/
     >   81 + (11*35**(0.5D0)*Q2(57,IELEM**3+1))/81 + (11*105**(0.5D0)
     >   *Q2(28,IELEM**3+1))/243 + (11*105**(0.5D0)*Q2(31,IELEM**3+1))
     >   /243 + (11*105**(0.5D0)*Q2(40,IELEM**3+1))/243 + (11*105**(0.
     >   5D0)*Q2(46,IELEM**3+1))/243 + (11*105**(0.5D0)*Q2(55,IELEM**3
     >   +1))/243 + (11*105**(0.5D0)*Q2(58,IELEM**3+1))/243) 
      CORNERQ(12,3) = (Q2(01,IELEM**3+1) + Q2(06,IELEM**3+1) - (5*Q
     >   2(11,IELEM**3+1))/3 - (77*Q2(16,IELEM**3+1))/27 + Q2(18,IELEM
     >   **3+1) + Q2(21,IELEM**3+1)/3 - (5*Q2(35,IELEM**3+1))/3 + (5*Q
     >   2(41,IELEM**3+1))/9 - (77*Q2(52,IELEM**3+1))/27 + (847*Q2(61,
     >   IELEM**3+1))/729 + 3**(0.5D0)*Q2(02,IELEM**3+1) + (3**(0.5D0)
     >   *Q2(05,IELEM**3+1))/3 + 5**(0.5D0)*Q2(03,IELEM**3+1) + 7**(0.
     >   5D0)*Q2(04,IELEM**3+1) - (5**(0.5D0)*Q2(09,IELEM**3+1))/3 + (
     >   3**(0.5D0)*Q2(17,IELEM**3+1))/3 - (11*7**(0.5D0)*Q2(13,IELEM*
     >   *3+1))/27 + (15**(0.5D0)*Q2(07,IELEM**3+1))/3 + (3**(0.5D0)*Q
     >   2(22,IELEM**3+1))/3 - (15**(0.5D0)*Q2(10,IELEM**3+1))/3 + (5*
     >   *(0.5D0)*Q2(23,IELEM**3+1))/3 + (21**(0.5D0)*Q2(08,IELEM**3+1
     >   ))/3 - (5*3**(0.5D0)*Q2(27,IELEM**3+1))/9 - (5**(0.5D0)*Q2(26
     >   ,IELEM**3+1))/3 + (7**(0.5D0)*Q2(24,IELEM**3+1))/3 + (15**(0.
     >   5D0)*Q2(19,IELEM**3+1))/3 - (77*3**(0.5D0)*Q2(32,IELEM**3+1))
     >   /81 - (11*21**(0.5D0)*Q2(14,IELEM**3+1))/27 - (11*7**(0.5D0)*
     >   Q2(30,IELEM**3+1))/27 - (5**(0.5D0)*Q2(33,IELEM**3+1))/3 - (1
     >   5**(0.5D0)*Q2(25,IELEM**3+1))/9 + (21**(0.5D0)*Q2(20,IELEM**3
     >   +1))/3 - (5*3**(0.5D0)*Q2(39,IELEM**3+1))/9 - (5**(0.5D0)*Q2(
     >   38,IELEM**3+1))/3 + (5*3**(0.5D0)*Q2(42,IELEM**3+1))/9 - (35*
     >   *(0.5D0)*Q2(12,IELEM**3+1))/3 + (5*5**(0.5D0)*Q2(43,IELEM**3+
     >   1))/9 - (15**(0.5D0)*Q2(34,IELEM**3+1))/3 - (11*21**(0.5D0)*Q
     >   2(29,IELEM**3+1))/81 - (11*35**(0.5D0)*Q2(15,IELEM**3+1))/27 
     >   + (5*7**(0.5D0)*Q2(44,IELEM**3+1))/9 - (15**(0.5D0)*Q2(37,IEL
     >   EM**3+1))/9 + (77*5**(0.5D0)*Q2(48,IELEM**3+1))/81 + (55*7**(
     >   0.5D0)*Q2(47,IELEM**3+1))/81 - (11*7**(0.5D0)*Q2(49,IELEM**3+
     >   1))/27 - (77*3**(0.5D0)*Q2(56,IELEM**3+1))/81 - (11*7**(0.5D0
     >   )*Q2(54,IELEM**3+1))/27 + (847*3**(0.5D0)*Q2(62,IELEM**3+1))/
     >   729 + (77*5**(0.5D0)*Q2(60,IELEM**3+1))/81 + (55*7**(0.5D0)*Q
     >   2(59,IELEM**3+1))/81 + (847*5**(0.5D0)*Q2(63,IELEM**3+1))/729
     >    + (847*7**(0.5D0)*Q2(64,IELEM**3+1))/729 - (11*21**(0.5D0)*Q
     >   2(50,IELEM**3+1))/27 - (35**(0.5D0)*Q2(36,IELEM**3+1))/3 - (1
     >   1*21**(0.5D0)*Q2(53,IELEM**3+1))/81 + (11*35**(0.5D0)*Q2(45,I
     >   ELEM**3+1))/81 - (11*35**(0.5D0)*Q2(51,IELEM**3+1))/27 + (11*
     >   35**(0.5D0)*Q2(57,IELEM**3+1))/81 - (105**(0.5D0)*Q2(28,IELEM
     >   **3+1))/9 - (11*105**(0.5D0)*Q2(31,IELEM**3+1))/81 - (105**(0
     >   .5D0)*Q2(40,IELEM**3+1))/9 + (11*105**(0.5D0)*Q2(46,IELEM**3+
     >   1))/81 - (11*105**(0.5D0)*Q2(55,IELEM**3+1))/81 + (11*105**(0
     >   .5D0)*Q2(58,IELEM**3+1))/81) 
      CORNERQ(13,3) = (Q2(01,IELEM**3+1) - 3*Q2(06,IELEM**3+1) + 5*
     >   Q2(11,IELEM**3+1) - 7*Q2(16,IELEM**3+1) - Q2(18,IELEM**3+1) +
     >    Q2(21,IELEM**3+1) - (5*Q2(35,IELEM**3+1))/3 - (5*Q2(41,IELEM
     >   **3+1))/3 + (77*Q2(52,IELEM**3+1))/27 - (77*Q2(61,IELEM**3+1)
     >   )/27 - 3**(0.5D0)*Q2(02,IELEM**3+1) + 3**(0.5D0)*Q2(05,IELEM*
     >   *3+1) + 5**(0.5D0)*Q2(03,IELEM**3+1) - 7**(0.5D0)*Q2(04,IELEM
     >   **3+1) + 5**(0.5D0)*Q2(09,IELEM**3+1) + (3**(0.5D0)*Q2(17,IEL
     >   EM**3+1))/3 + 7**(0.5D0)*Q2(13,IELEM**3+1) + 15**(0.5D0)*Q2(0
     >   7,IELEM**3+1) - 3**(0.5D0)*Q2(22,IELEM**3+1) - 15**(0.5D0)*Q2
     >   (10,IELEM**3+1) + 5**(0.5D0)*Q2(23,IELEM**3+1) - 21**(0.5D0)*
     >   Q2(08,IELEM**3+1) + (5*3**(0.5D0)*Q2(27,IELEM**3+1))/3 - 5**(
     >   0.5D0)*Q2(26,IELEM**3+1) - 7**(0.5D0)*Q2(24,IELEM**3+1) + (15
     >   **(0.5D0)*Q2(19,IELEM**3+1))/3 - (7*3**(0.5D0)*Q2(32,IELEM**3
     >   +1))/3 - 21**(0.5D0)*Q2(14,IELEM**3+1) - 7**(0.5D0)*Q2(30,IEL
     >   EM**3+1) - (5**(0.5D0)*Q2(33,IELEM**3+1))/3 + (15**(0.5D0)*Q2
     >   (25,IELEM**3+1))/3 - (21**(0.5D0)*Q2(20,IELEM**3+1))/3 - (5*3
     >   **(0.5D0)*Q2(39,IELEM**3+1))/3 + 5**(0.5D0)*Q2(38,IELEM**3+1)
     >    + (5*3**(0.5D0)*Q2(42,IELEM**3+1))/3 - 35**(0.5D0)*Q2(12,IEL
     >   EM**3+1) - (5*5**(0.5D0)*Q2(43,IELEM**3+1))/3 + (15**(0.5D0)*
     >   Q2(34,IELEM**3+1))/3 + (21**(0.5D0)*Q2(29,IELEM**3+1))/3 + 35
     >   **(0.5D0)*Q2(15,IELEM**3+1) + (5*7**(0.5D0)*Q2(44,IELEM**3+1)
     >   )/3 - (15**(0.5D0)*Q2(37,IELEM**3+1))/3 + (7*5**(0.5D0)*Q2(48
     >   ,IELEM**3+1))/3 - (5*7**(0.5D0)*Q2(47,IELEM**3+1))/3 - (11*7*
     >   *(0.5D0)*Q2(49,IELEM**3+1))/27 + (77*3**(0.5D0)*Q2(56,IELEM**
     >   3+1))/27 + (11*7**(0.5D0)*Q2(54,IELEM**3+1))/9 + (77*3**(0.5D
     >   0)*Q2(62,IELEM**3+1))/27 + (77*5**(0.5D0)*Q2(60,IELEM**3+1))/
     >   27 - (55*7**(0.5D0)*Q2(59,IELEM**3+1))/27 - (77*5**(0.5D0)*Q2
     >   (63,IELEM**3+1))/27 + (77*7**(0.5D0)*Q2(64,IELEM**3+1))/27 + 
     >   (11*21**(0.5D0)*Q2(50,IELEM**3+1))/27 + (35**(0.5D0)*Q2(36,IE
     >   LEM**3+1))/3 - (11*21**(0.5D0)*Q2(53,IELEM**3+1))/27 - (35**(
     >   0.5D0)*Q2(45,IELEM**3+1))/3 - (11*35**(0.5D0)*Q2(51,IELEM**3+
     >   1))/27 - (11*35**(0.5D0)*Q2(57,IELEM**3+1))/27 - (105**(0.5D0
     >   )*Q2(28,IELEM**3+1))/3 + (105**(0.5D0)*Q2(31,IELEM**3+1))/3 +
     >    (105**(0.5D0)*Q2(40,IELEM**3+1))/3 + (105**(0.5D0)*Q2(46,IEL
     >   EM**3+1))/3 - (11*105**(0.5D0)*Q2(55,IELEM**3+1))/27 + (11*10
     >   5**(0.5D0)*Q2(58,IELEM**3+1))/27) 
      CORNERQ(14,3) = (Q2(01,IELEM**3+1) - Q2(06,IELEM**3+1) - (5*Q
     >   2(11,IELEM**3+1))/3 + (77*Q2(16,IELEM**3+1))/27 - Q2(18,IELEM
     >   **3+1)/3 + Q2(21,IELEM**3+1) + (5*Q2(35,IELEM**3+1))/9 - (5*Q
     >   2(41,IELEM**3+1))/3 - (847*Q2(52,IELEM**3+1))/729 - (77*Q2(61
     >   ,IELEM**3+1))/27 - (3**(0.5D0)*Q2(02,IELEM**3+1))/3 + 3**(0.5
     >   D0)*Q2(05,IELEM**3+1) - (5**(0.5D0)*Q2(03,IELEM**3+1))/3 + (1
     >   1*7**(0.5D0)*Q2(04,IELEM**3+1))/27 + 5**(0.5D0)*Q2(09,IELEM**
     >   3+1) + (3**(0.5D0)*Q2(17,IELEM**3+1))/3 + 7**(0.5D0)*Q2(13,IE
     >   LEM**3+1) - (15**(0.5D0)*Q2(07,IELEM**3+1))/3 - (3**(0.5D0)*Q
     >   2(22,IELEM**3+1))/3 - (15**(0.5D0)*Q2(10,IELEM**3+1))/3 - (5*
     >   *(0.5D0)*Q2(23,IELEM**3+1))/3 + (11*21**(0.5D0)*Q2(08,IELEM**
     >   3+1))/27 - (5*3**(0.5D0)*Q2(27,IELEM**3+1))/9 - (5**(0.5D0)*Q
     >   2(26,IELEM**3+1))/3 + (11*7**(0.5D0)*Q2(24,IELEM**3+1))/27 - 
     >   (15**(0.5D0)*Q2(19,IELEM**3+1))/9 + (77*3**(0.5D0)*Q2(32,IELE
     >   M**3+1))/81 - (21**(0.5D0)*Q2(14,IELEM**3+1))/3 - (7**(0.5D0)
     >   *Q2(30,IELEM**3+1))/3 - (5**(0.5D0)*Q2(33,IELEM**3+1))/3 + (1
     >   5**(0.5D0)*Q2(25,IELEM**3+1))/3 + (11*21**(0.5D0)*Q2(20,IELEM
     >   **3+1))/81 + (5*3**(0.5D0)*Q2(39,IELEM**3+1))/9 + (5**(0.5D0)
     >   *Q2(38,IELEM**3+1))/3 + (5*3**(0.5D0)*Q2(42,IELEM**3+1))/9 + 
     >   (11*35**(0.5D0)*Q2(12,IELEM**3+1))/27 + (5*5**(0.5D0)*Q2(43,I
     >   ELEM**3+1))/9 + (15**(0.5D0)*Q2(34,IELEM**3+1))/9 + (21**(0.5
     >   D0)*Q2(29,IELEM**3+1))/3 - (35**(0.5D0)*Q2(15,IELEM**3+1))/3 
     >   - (55*7**(0.5D0)*Q2(44,IELEM**3+1))/81 - (15**(0.5D0)*Q2(37,I
     >   ELEM**3+1))/3 - (77*5**(0.5D0)*Q2(48,IELEM**3+1))/81 + (5*7**
     >   (0.5D0)*Q2(47,IELEM**3+1))/9 - (11*7**(0.5D0)*Q2(49,IELEM**3+
     >   1))/27 - (847*3**(0.5D0)*Q2(56,IELEM**3+1))/729 + (11*7**(0.5
     >   D0)*Q2(54,IELEM**3+1))/27 + (77*3**(0.5D0)*Q2(62,IELEM**3+1))
     >   /81 - (847*5**(0.5D0)*Q2(60,IELEM**3+1))/729 + (55*7**(0.5D0)
     >   *Q2(59,IELEM**3+1))/81 + (77*5**(0.5D0)*Q2(63,IELEM**3+1))/81
     >    - (847*7**(0.5D0)*Q2(64,IELEM**3+1))/729 + (11*21**(0.5D0)*Q
     >   2(50,IELEM**3+1))/81 - (11*35**(0.5D0)*Q2(36,IELEM**3+1))/81 
     >   - (11*21**(0.5D0)*Q2(53,IELEM**3+1))/27 - (35**(0.5D0)*Q2(45,
     >   IELEM**3+1))/3 + (11*35**(0.5D0)*Q2(51,IELEM**3+1))/81 - (11*
     >   35**(0.5D0)*Q2(57,IELEM**3+1))/27 + (11*105**(0.5D0)*Q2(28,IE
     >   LEM**3+1))/81 - (105**(0.5D0)*Q2(31,IELEM**3+1))/9 - (11*105*
     >   *(0.5D0)*Q2(40,IELEM**3+1))/81 + (105**(0.5D0)*Q2(46,IELEM**3
     >   +1))/9 + (11*105**(0.5D0)*Q2(55,IELEM**3+1))/81 + (11*105**(0
     >   .5D0)*Q2(58,IELEM**3+1))/81) 
      CORNERQ(15,3) = (Q2(01,IELEM**3+1) + Q2(06,IELEM**3+1) - (5*Q
     >   2(11,IELEM**3+1))/3 - (77*Q2(16,IELEM**3+1))/27 + Q2(18,IELEM
     >   **3+1)/3 + Q2(21,IELEM**3+1) + (5*Q2(35,IELEM**3+1))/9 - (5*Q
     >   2(41,IELEM**3+1))/3 + (847*Q2(52,IELEM**3+1))/729 - (77*Q2(61
     >   ,IELEM**3+1))/27 + (3**(0.5D0)*Q2(02,IELEM**3+1))/3 + 3**(0.5
     >   D0)*Q2(05,IELEM**3+1) - (5**(0.5D0)*Q2(03,IELEM**3+1))/3 - (1
     >   1*7**(0.5D0)*Q2(04,IELEM**3+1))/27 + 5**(0.5D0)*Q2(09,IELEM**
     >   3+1) + (3**(0.5D0)*Q2(17,IELEM**3+1))/3 + 7**(0.5D0)*Q2(13,IE
     >   LEM**3+1) - (15**(0.5D0)*Q2(07,IELEM**3+1))/3 + (3**(0.5D0)*Q
     >   2(22,IELEM**3+1))/3 + (15**(0.5D0)*Q2(10,IELEM**3+1))/3 - (5*
     >   *(0.5D0)*Q2(23,IELEM**3+1))/3 - (11*21**(0.5D0)*Q2(08,IELEM**
     >   3+1))/27 - (5*3**(0.5D0)*Q2(27,IELEM**3+1))/9 + (5**(0.5D0)*Q
     >   2(26,IELEM**3+1))/3 - (11*7**(0.5D0)*Q2(24,IELEM**3+1))/27 - 
     >   (15**(0.5D0)*Q2(19,IELEM**3+1))/9 - (77*3**(0.5D0)*Q2(32,IELE
     >   M**3+1))/81 + (21**(0.5D0)*Q2(14,IELEM**3+1))/3 + (7**(0.5D0)
     >   *Q2(30,IELEM**3+1))/3 - (5**(0.5D0)*Q2(33,IELEM**3+1))/3 + (1
     >   5**(0.5D0)*Q2(25,IELEM**3+1))/3 - (11*21**(0.5D0)*Q2(20,IELEM
     >   **3+1))/81 + (5*3**(0.5D0)*Q2(39,IELEM**3+1))/9 - (5**(0.5D0)
     >   *Q2(38,IELEM**3+1))/3 - (5*3**(0.5D0)*Q2(42,IELEM**3+1))/9 - 
     >   (11*35**(0.5D0)*Q2(12,IELEM**3+1))/27 + (5*5**(0.5D0)*Q2(43,I
     >   ELEM**3+1))/9 - (15**(0.5D0)*Q2(34,IELEM**3+1))/9 + (21**(0.5
     >   D0)*Q2(29,IELEM**3+1))/3 - (35**(0.5D0)*Q2(15,IELEM**3+1))/3 
     >   + (55*7**(0.5D0)*Q2(44,IELEM**3+1))/81 - (15**(0.5D0)*Q2(37,I
     >   ELEM**3+1))/3 + (77*5**(0.5D0)*Q2(48,IELEM**3+1))/81 + (5*7**
     >   (0.5D0)*Q2(47,IELEM**3+1))/9 - (11*7**(0.5D0)*Q2(49,IELEM**3+
     >   1))/27 + (847*3**(0.5D0)*Q2(56,IELEM**3+1))/729 - (11*7**(0.5
     >   D0)*Q2(54,IELEM**3+1))/27 - (77*3**(0.5D0)*Q2(62,IELEM**3+1))
     >   /81 + (847*5**(0.5D0)*Q2(60,IELEM**3+1))/729 + (55*7**(0.5D0)
     >   *Q2(59,IELEM**3+1))/81 + (77*5**(0.5D0)*Q2(63,IELEM**3+1))/81
     >    + (847*7**(0.5D0)*Q2(64,IELEM**3+1))/729 - (11*21**(0.5D0)*Q
     >   2(50,IELEM**3+1))/81 + (11*35**(0.5D0)*Q2(36,IELEM**3+1))/81 
     >   - (11*21**(0.5D0)*Q2(53,IELEM**3+1))/27 - (35**(0.5D0)*Q2(45,
     >   IELEM**3+1))/3 + (11*35**(0.5D0)*Q2(51,IELEM**3+1))/81 - (11*
     >   35**(0.5D0)*Q2(57,IELEM**3+1))/27 - (11*105**(0.5D0)*Q2(28,IE
     >   LEM**3+1))/81 - (105**(0.5D0)*Q2(31,IELEM**3+1))/9 + (11*105*
     >   *(0.5D0)*Q2(40,IELEM**3+1))/81 - (105**(0.5D0)*Q2(46,IELEM**3
     >   +1))/9 + (11*105**(0.5D0)*Q2(55,IELEM**3+1))/81 - (11*105**(0
     >   .5D0)*Q2(58,IELEM**3+1))/81) 
      CORNERQ(16,3) = (Q2(01,IELEM**3+1) + 3*Q2(06,IELEM**3+1) + 5*
     >   Q2(11,IELEM**3+1) + 7*Q2(16,IELEM**3+1) + Q2(18,IELEM**3+1) +
     >    Q2(21,IELEM**3+1) - (5*Q2(35,IELEM**3+1))/3 - (5*Q2(41,IELEM
     >   **3+1))/3 - (77*Q2(52,IELEM**3+1))/27 - (77*Q2(61,IELEM**3+1)
     >   )/27 + 3**(0.5D0)*Q2(02,IELEM**3+1) + 3**(0.5D0)*Q2(05,IELEM*
     >   *3+1) + 5**(0.5D0)*Q2(03,IELEM**3+1) + 7**(0.5D0)*Q2(04,IELEM
     >   **3+1) + 5**(0.5D0)*Q2(09,IELEM**3+1) + (3**(0.5D0)*Q2(17,IEL
     >   EM**3+1))/3 + 7**(0.5D0)*Q2(13,IELEM**3+1) + 15**(0.5D0)*Q2(0
     >   7,IELEM**3+1) + 3**(0.5D0)*Q2(22,IELEM**3+1) + 15**(0.5D0)*Q2
     >   (10,IELEM**3+1) + 5**(0.5D0)*Q2(23,IELEM**3+1) + 21**(0.5D0)*
     >   Q2(08,IELEM**3+1) + (5*3**(0.5D0)*Q2(27,IELEM**3+1))/3 + 5**(
     >   0.5D0)*Q2(26,IELEM**3+1) + 7**(0.5D0)*Q2(24,IELEM**3+1) + (15
     >   **(0.5D0)*Q2(19,IELEM**3+1))/3 + (7*3**(0.5D0)*Q2(32,IELEM**3
     >   +1))/3 + 21**(0.5D0)*Q2(14,IELEM**3+1) + 7**(0.5D0)*Q2(30,IEL
     >   EM**3+1) - (5**(0.5D0)*Q2(33,IELEM**3+1))/3 + (15**(0.5D0)*Q2
     >   (25,IELEM**3+1))/3 + (21**(0.5D0)*Q2(20,IELEM**3+1))/3 - (5*3
     >   **(0.5D0)*Q2(39,IELEM**3+1))/3 - 5**(0.5D0)*Q2(38,IELEM**3+1)
     >    - (5*3**(0.5D0)*Q2(42,IELEM**3+1))/3 + 35**(0.5D0)*Q2(12,IEL
     >   EM**3+1) - (5*5**(0.5D0)*Q2(43,IELEM**3+1))/3 - (15**(0.5D0)*
     >   Q2(34,IELEM**3+1))/3 + (21**(0.5D0)*Q2(29,IELEM**3+1))/3 + 35
     >   **(0.5D0)*Q2(15,IELEM**3+1) - (5*7**(0.5D0)*Q2(44,IELEM**3+1)
     >   )/3 - (15**(0.5D0)*Q2(37,IELEM**3+1))/3 - (7*5**(0.5D0)*Q2(48
     >   ,IELEM**3+1))/3 - (5*7**(0.5D0)*Q2(47,IELEM**3+1))/3 - (11*7*
     >   *(0.5D0)*Q2(49,IELEM**3+1))/27 - (77*3**(0.5D0)*Q2(56,IELEM**
     >   3+1))/27 - (11*7**(0.5D0)*Q2(54,IELEM**3+1))/9 - (77*3**(0.5D
     >   0)*Q2(62,IELEM**3+1))/27 - (77*5**(0.5D0)*Q2(60,IELEM**3+1))/
     >   27 - (55*7**(0.5D0)*Q2(59,IELEM**3+1))/27 - (77*5**(0.5D0)*Q2
     >   (63,IELEM**3+1))/27 - (77*7**(0.5D0)*Q2(64,IELEM**3+1))/27 - 
     >   (11*21**(0.5D0)*Q2(50,IELEM**3+1))/27 - (35**(0.5D0)*Q2(36,IE
     >   LEM**3+1))/3 - (11*21**(0.5D0)*Q2(53,IELEM**3+1))/27 - (35**(
     >   0.5D0)*Q2(45,IELEM**3+1))/3 - (11*35**(0.5D0)*Q2(51,IELEM**3+
     >   1))/27 - (11*35**(0.5D0)*Q2(57,IELEM**3+1))/27 + (105**(0.5D0
     >   )*Q2(28,IELEM**3+1))/3 + (105**(0.5D0)*Q2(31,IELEM**3+1))/3 -
     >    (105**(0.5D0)*Q2(40,IELEM**3+1))/3 - (105**(0.5D0)*Q2(46,IEL
     >   EM**3+1))/3 - (11*105**(0.5D0)*Q2(55,IELEM**3+1))/27 - (11*10
     >   5**(0.5D0)*Q2(58,IELEM**3+1))/27) 
      CORNERQ(01,4) = (Q2(01,IELEM**3+1) + 3*Q2(06,IELEM**3+1) + 5*
     >   Q2(11,IELEM**3+1) + 7*Q2(16,IELEM**3+1) - 3*Q2(18,IELEM**3+1)
     >    - 3*Q2(21,IELEM**3+1) + 5*Q2(35,IELEM**3+1) + 5*Q2(41,IELEM*
     >   *3+1) - 7*Q2(52,IELEM**3+1) - 7*Q2(61,IELEM**3+1) - 3**(0.5D0
     >   )*Q2(02,IELEM**3+1) - 3**(0.5D0)*Q2(05,IELEM**3+1) + 5**(0.5D
     >   0)*Q2(03,IELEM**3+1) - 7**(0.5D0)*Q2(04,IELEM**3+1) + 5**(0.5
     >   D0)*Q2(09,IELEM**3+1) + 3**(0.5D0)*Q2(17,IELEM**3+1) - 7**(0.
     >   5D0)*Q2(13,IELEM**3+1) - 15**(0.5D0)*Q2(07,IELEM**3+1) + 3*3*
     >   *(0.5D0)*Q2(22,IELEM**3+1) - 15**(0.5D0)*Q2(10,IELEM**3+1) - 
     >   3*5**(0.5D0)*Q2(23,IELEM**3+1) + 21**(0.5D0)*Q2(08,IELEM**3+1
     >   ) + 5*3**(0.5D0)*Q2(27,IELEM**3+1) - 3*5**(0.5D0)*Q2(26,IELEM
     >   **3+1) + 3*7**(0.5D0)*Q2(24,IELEM**3+1) + 15**(0.5D0)*Q2(19,I
     >   ELEM**3+1) + 7*3**(0.5D0)*Q2(32,IELEM**3+1) + 21**(0.5D0)*Q2(
     >   14,IELEM**3+1) + 3*7**(0.5D0)*Q2(30,IELEM**3+1) + 5**(0.5D0)*
     >   Q2(33,IELEM**3+1) + 15**(0.5D0)*Q2(25,IELEM**3+1) - 21**(0.5D
     >   0)*Q2(20,IELEM**3+1) - 5*3**(0.5D0)*Q2(39,IELEM**3+1) + 3*5**
     >   (0.5D0)*Q2(38,IELEM**3+1) - 5*3**(0.5D0)*Q2(42,IELEM**3+1) - 
     >   35**(0.5D0)*Q2(12,IELEM**3+1) + 5*5**(0.5D0)*Q2(43,IELEM**3+1
     >   ) - 15**(0.5D0)*Q2(34,IELEM**3+1) - 21**(0.5D0)*Q2(29,IELEM**
     >   3+1) - 35**(0.5D0)*Q2(15,IELEM**3+1) - 5*7**(0.5D0)*Q2(44,IEL
     >   EM**3+1) - 15**(0.5D0)*Q2(37,IELEM**3+1) + 7*5**(0.5D0)*Q2(48
     >   ,IELEM**3+1) - 5*7**(0.5D0)*Q2(47,IELEM**3+1) + 7**(0.5D0)*Q2
     >   (49,IELEM**3+1) + 7*3**(0.5D0)*Q2(56,IELEM**3+1) + 3*7**(0.5D
     >   0)*Q2(54,IELEM**3+1) + 7*3**(0.5D0)*Q2(62,IELEM**3+1) - 7*5**
     >   (0.5D0)*Q2(60,IELEM**3+1) + 5*7**(0.5D0)*Q2(59,IELEM**3+1) - 
     >   7*5**(0.5D0)*Q2(63,IELEM**3+1) + 7*7**(0.5D0)*Q2(64,IELEM**3+
     >   1) - 21**(0.5D0)*Q2(50,IELEM**3+1) - 35**(0.5D0)*Q2(36,IELEM*
     >   *3+1) - 21**(0.5D0)*Q2(53,IELEM**3+1) - 35**(0.5D0)*Q2(45,IEL
     >   EM**3+1) + 35**(0.5D0)*Q2(51,IELEM**3+1) + 35**(0.5D0)*Q2(57,
     >   IELEM**3+1) - 105**(0.5D0)*Q2(28,IELEM**3+1) - 105**(0.5D0)*Q
     >   2(31,IELEM**3+1) + 105**(0.5D0)*Q2(40,IELEM**3+1) + 105**(0.5
     >   D0)*Q2(46,IELEM**3+1) - 105**(0.5D0)*Q2(55,IELEM**3+1) - 105*
     >   *(0.5D0)*Q2(58,IELEM**3+1)) 
      CORNERQ(02,4) = (Q2(01,IELEM**3+1) + Q2(06,IELEM**3+1) - (5*Q
     >   2(11,IELEM**3+1))/3 - (77*Q2(16,IELEM**3+1))/27 - Q2(18,IELEM
     >   **3+1) - 3*Q2(21,IELEM**3+1) - (5*Q2(35,IELEM**3+1))/3 + 5*Q2
     >   (41,IELEM**3+1) + (77*Q2(52,IELEM**3+1))/27 - 7*Q2(61,IELEM**
     >   3+1) - (3**(0.5D0)*Q2(02,IELEM**3+1))/3 - 3**(0.5D0)*Q2(05,IE
     >   LEM**3+1) - (5**(0.5D0)*Q2(03,IELEM**3+1))/3 + (11*7**(0.5D0)
     >   *Q2(04,IELEM**3+1))/27 + 5**(0.5D0)*Q2(09,IELEM**3+1) + 3**(0
     >   .5D0)*Q2(17,IELEM**3+1) - 7**(0.5D0)*Q2(13,IELEM**3+1) + (15*
     >   *(0.5D0)*Q2(07,IELEM**3+1))/3 + 3**(0.5D0)*Q2(22,IELEM**3+1) 
     >   - (15**(0.5D0)*Q2(10,IELEM**3+1))/3 + 5**(0.5D0)*Q2(23,IELEM*
     >   *3+1) - (11*21**(0.5D0)*Q2(08,IELEM**3+1))/27 - (5*3**(0.5D0)
     >   *Q2(27,IELEM**3+1))/3 - 5**(0.5D0)*Q2(26,IELEM**3+1) - (11*7*
     >   *(0.5D0)*Q2(24,IELEM**3+1))/9 - (15**(0.5D0)*Q2(19,IELEM**3+1
     >   ))/3 - (77*3**(0.5D0)*Q2(32,IELEM**3+1))/27 + (21**(0.5D0)*Q2
     >   (14,IELEM**3+1))/3 + 7**(0.5D0)*Q2(30,IELEM**3+1) + 5**(0.5D0
     >   )*Q2(33,IELEM**3+1) + 15**(0.5D0)*Q2(25,IELEM**3+1) + (11*21*
     >   *(0.5D0)*Q2(20,IELEM**3+1))/27 + (5*3**(0.5D0)*Q2(39,IELEM**3
     >   +1))/3 + 5**(0.5D0)*Q2(38,IELEM**3+1) - (5*3**(0.5D0)*Q2(42,I
     >   ELEM**3+1))/3 + (11*35**(0.5D0)*Q2(12,IELEM**3+1))/27 - (5*5*
     >   *(0.5D0)*Q2(43,IELEM**3+1))/3 - (15**(0.5D0)*Q2(34,IELEM**3+1
     >   ))/3 - 21**(0.5D0)*Q2(29,IELEM**3+1) + (35**(0.5D0)*Q2(15,IEL
     >   EM**3+1))/3 + (55*7**(0.5D0)*Q2(44,IELEM**3+1))/27 - 15**(0.5
     >   D0)*Q2(37,IELEM**3+1) - (77*5**(0.5D0)*Q2(48,IELEM**3+1))/27 
     >   + (5*7**(0.5D0)*Q2(47,IELEM**3+1))/3 + 7**(0.5D0)*Q2(49,IELEM
     >   **3+1) - (77*3**(0.5D0)*Q2(56,IELEM**3+1))/27 + 7**(0.5D0)*Q2
     >   (54,IELEM**3+1) + (7*3**(0.5D0)*Q2(62,IELEM**3+1))/3 + (77*5*
     >   *(0.5D0)*Q2(60,IELEM**3+1))/27 - (5*7**(0.5D0)*Q2(59,IELEM**3
     >   +1))/3 + (7*5**(0.5D0)*Q2(63,IELEM**3+1))/3 - (77*7**(0.5D0)*
     >   Q2(64,IELEM**3+1))/27 - (21**(0.5D0)*Q2(50,IELEM**3+1))/3 + (
     >   11*35**(0.5D0)*Q2(36,IELEM**3+1))/27 - 21**(0.5D0)*Q2(53,IELE
     >   M**3+1) - 35**(0.5D0)*Q2(45,IELEM**3+1) - (35**(0.5D0)*Q2(51,
     >   IELEM**3+1))/3 + 35**(0.5D0)*Q2(57,IELEM**3+1) + (11*105**(0.
     >   5D0)*Q2(28,IELEM**3+1))/27 + (105**(0.5D0)*Q2(31,IELEM**3+1))
     >   /3 - (11*105**(0.5D0)*Q2(40,IELEM**3+1))/27 + (105**(0.5D0)*Q
     >   2(46,IELEM**3+1))/3 + (105**(0.5D0)*Q2(55,IELEM**3+1))/3 - (1
     >   05**(0.5D0)*Q2(58,IELEM**3+1))/3) 
      CORNERQ(03,4) = (Q2(01,IELEM**3+1) - Q2(06,IELEM**3+1) - (5*Q
     >   2(11,IELEM**3+1))/3 + (77*Q2(16,IELEM**3+1))/27 + Q2(18,IELEM
     >   **3+1) - 3*Q2(21,IELEM**3+1) - (5*Q2(35,IELEM**3+1))/3 + 5*Q2
     >   (41,IELEM**3+1) - (77*Q2(52,IELEM**3+1))/27 - 7*Q2(61,IELEM**
     >   3+1) + (3**(0.5D0)*Q2(02,IELEM**3+1))/3 - 3**(0.5D0)*Q2(05,IE
     >   LEM**3+1) - (5**(0.5D0)*Q2(03,IELEM**3+1))/3 - (11*7**(0.5D0)
     >   *Q2(04,IELEM**3+1))/27 + 5**(0.5D0)*Q2(09,IELEM**3+1) + 3**(0
     >   .5D0)*Q2(17,IELEM**3+1) - 7**(0.5D0)*Q2(13,IELEM**3+1) + (15*
     >   *(0.5D0)*Q2(07,IELEM**3+1))/3 - 3**(0.5D0)*Q2(22,IELEM**3+1) 
     >   + (15**(0.5D0)*Q2(10,IELEM**3+1))/3 + 5**(0.5D0)*Q2(23,IELEM*
     >   *3+1) + (11*21**(0.5D0)*Q2(08,IELEM**3+1))/27 - (5*3**(0.5D0)
     >   *Q2(27,IELEM**3+1))/3 + 5**(0.5D0)*Q2(26,IELEM**3+1) + (11*7*
     >   *(0.5D0)*Q2(24,IELEM**3+1))/9 - (15**(0.5D0)*Q2(19,IELEM**3+1
     >   ))/3 + (77*3**(0.5D0)*Q2(32,IELEM**3+1))/27 - (21**(0.5D0)*Q2
     >   (14,IELEM**3+1))/3 - 7**(0.5D0)*Q2(30,IELEM**3+1) + 5**(0.5D0
     >   )*Q2(33,IELEM**3+1) + 15**(0.5D0)*Q2(25,IELEM**3+1) - (11*21*
     >   *(0.5D0)*Q2(20,IELEM**3+1))/27 + (5*3**(0.5D0)*Q2(39,IELEM**3
     >   +1))/3 - 5**(0.5D0)*Q2(38,IELEM**3+1) + (5*3**(0.5D0)*Q2(42,I
     >   ELEM**3+1))/3 - (11*35**(0.5D0)*Q2(12,IELEM**3+1))/27 - (5*5*
     >   *(0.5D0)*Q2(43,IELEM**3+1))/3 + (15**(0.5D0)*Q2(34,IELEM**3+1
     >   ))/3 - 21**(0.5D0)*Q2(29,IELEM**3+1) + (35**(0.5D0)*Q2(15,IEL
     >   EM**3+1))/3 - (55*7**(0.5D0)*Q2(44,IELEM**3+1))/27 - 15**(0.5
     >   D0)*Q2(37,IELEM**3+1) + (77*5**(0.5D0)*Q2(48,IELEM**3+1))/27 
     >   + (5*7**(0.5D0)*Q2(47,IELEM**3+1))/3 + 7**(0.5D0)*Q2(49,IELEM
     >   **3+1) + (77*3**(0.5D0)*Q2(56,IELEM**3+1))/27 - 7**(0.5D0)*Q2
     >   (54,IELEM**3+1) - (7*3**(0.5D0)*Q2(62,IELEM**3+1))/3 - (77*5*
     >   *(0.5D0)*Q2(60,IELEM**3+1))/27 - (5*7**(0.5D0)*Q2(59,IELEM**3
     >   +1))/3 + (7*5**(0.5D0)*Q2(63,IELEM**3+1))/3 + (77*7**(0.5D0)*
     >   Q2(64,IELEM**3+1))/27 + (21**(0.5D0)*Q2(50,IELEM**3+1))/3 - (
     >   11*35**(0.5D0)*Q2(36,IELEM**3+1))/27 - 21**(0.5D0)*Q2(53,IELE
     >   M**3+1) - 35**(0.5D0)*Q2(45,IELEM**3+1) - (35**(0.5D0)*Q2(51,
     >   IELEM**3+1))/3 + 35**(0.5D0)*Q2(57,IELEM**3+1) - (11*105**(0.
     >   5D0)*Q2(28,IELEM**3+1))/27 + (105**(0.5D0)*Q2(31,IELEM**3+1))
     >   /3 + (11*105**(0.5D0)*Q2(40,IELEM**3+1))/27 - (105**(0.5D0)*Q
     >   2(46,IELEM**3+1))/3 + (105**(0.5D0)*Q2(55,IELEM**3+1))/3 + (1
     >   05**(0.5D0)*Q2(58,IELEM**3+1))/3) 
      CORNERQ(04,4) = (Q2(01,IELEM**3+1) - 3*Q2(06,IELEM**3+1) + 5*
     >   Q2(11,IELEM**3+1) - 7*Q2(16,IELEM**3+1) + 3*Q2(18,IELEM**3+1)
     >    - 3*Q2(21,IELEM**3+1) + 5*Q2(35,IELEM**3+1) + 5*Q2(41,IELEM*
     >   *3+1) + 7*Q2(52,IELEM**3+1) - 7*Q2(61,IELEM**3+1) + 3**(0.5D0
     >   )*Q2(02,IELEM**3+1) - 3**(0.5D0)*Q2(05,IELEM**3+1) + 5**(0.5D
     >   0)*Q2(03,IELEM**3+1) + 7**(0.5D0)*Q2(04,IELEM**3+1) + 5**(0.5
     >   D0)*Q2(09,IELEM**3+1) + 3**(0.5D0)*Q2(17,IELEM**3+1) - 7**(0.
     >   5D0)*Q2(13,IELEM**3+1) - 15**(0.5D0)*Q2(07,IELEM**3+1) - 3*3*
     >   *(0.5D0)*Q2(22,IELEM**3+1) + 15**(0.5D0)*Q2(10,IELEM**3+1) - 
     >   3*5**(0.5D0)*Q2(23,IELEM**3+1) - 21**(0.5D0)*Q2(08,IELEM**3+1
     >   ) + 5*3**(0.5D0)*Q2(27,IELEM**3+1) + 3*5**(0.5D0)*Q2(26,IELEM
     >   **3+1) - 3*7**(0.5D0)*Q2(24,IELEM**3+1) + 15**(0.5D0)*Q2(19,I
     >   ELEM**3+1) - 7*3**(0.5D0)*Q2(32,IELEM**3+1) - 21**(0.5D0)*Q2(
     >   14,IELEM**3+1) - 3*7**(0.5D0)*Q2(30,IELEM**3+1) + 5**(0.5D0)*
     >   Q2(33,IELEM**3+1) + 15**(0.5D0)*Q2(25,IELEM**3+1) + 21**(0.5D
     >   0)*Q2(20,IELEM**3+1) - 5*3**(0.5D0)*Q2(39,IELEM**3+1) - 3*5**
     >   (0.5D0)*Q2(38,IELEM**3+1) + 5*3**(0.5D0)*Q2(42,IELEM**3+1) + 
     >   35**(0.5D0)*Q2(12,IELEM**3+1) + 5*5**(0.5D0)*Q2(43,IELEM**3+1
     >   ) + 15**(0.5D0)*Q2(34,IELEM**3+1) - 21**(0.5D0)*Q2(29,IELEM**
     >   3+1) - 35**(0.5D0)*Q2(15,IELEM**3+1) + 5*7**(0.5D0)*Q2(44,IEL
     >   EM**3+1) - 15**(0.5D0)*Q2(37,IELEM**3+1) - 7*5**(0.5D0)*Q2(48
     >   ,IELEM**3+1) - 5*7**(0.5D0)*Q2(47,IELEM**3+1) + 7**(0.5D0)*Q2
     >   (49,IELEM**3+1) - 7*3**(0.5D0)*Q2(56,IELEM**3+1) - 3*7**(0.5D
     >   0)*Q2(54,IELEM**3+1) - 7*3**(0.5D0)*Q2(62,IELEM**3+1) + 7*5**
     >   (0.5D0)*Q2(60,IELEM**3+1) + 5*7**(0.5D0)*Q2(59,IELEM**3+1) - 
     >   7*5**(0.5D0)*Q2(63,IELEM**3+1) - 7*7**(0.5D0)*Q2(64,IELEM**3+
     >   1) + 21**(0.5D0)*Q2(50,IELEM**3+1) + 35**(0.5D0)*Q2(36,IELEM*
     >   *3+1) - 21**(0.5D0)*Q2(53,IELEM**3+1) - 35**(0.5D0)*Q2(45,IEL
     >   EM**3+1) + 35**(0.5D0)*Q2(51,IELEM**3+1) + 35**(0.5D0)*Q2(57,
     >   IELEM**3+1) + 105**(0.5D0)*Q2(28,IELEM**3+1) - 105**(0.5D0)*Q
     >   2(31,IELEM**3+1) - 105**(0.5D0)*Q2(40,IELEM**3+1) - 105**(0.5
     >   D0)*Q2(46,IELEM**3+1) - 105**(0.5D0)*Q2(55,IELEM**3+1) + 105*
     >   *(0.5D0)*Q2(58,IELEM**3+1)) 
      CORNERQ(05,4) = (Q2(01,IELEM**3+1) + Q2(06,IELEM**3+1) - (5*Q
     >   2(11,IELEM**3+1))/3 - (77*Q2(16,IELEM**3+1))/27 - 3*Q2(18,IEL
     >   EM**3+1) - Q2(21,IELEM**3+1) + 5*Q2(35,IELEM**3+1) - (5*Q2(41
     >   ,IELEM**3+1))/3 - 7*Q2(52,IELEM**3+1) + (77*Q2(61,IELEM**3+1)
     >   )/27 - 3**(0.5D0)*Q2(02,IELEM**3+1) - (3**(0.5D0)*Q2(05,IELEM
     >   **3+1))/3 + 5**(0.5D0)*Q2(03,IELEM**3+1) - 7**(0.5D0)*Q2(04,I
     >   ELEM**3+1) - (5**(0.5D0)*Q2(09,IELEM**3+1))/3 + 3**(0.5D0)*Q2
     >   (17,IELEM**3+1) + (11*7**(0.5D0)*Q2(13,IELEM**3+1))/27 - (15*
     >   *(0.5D0)*Q2(07,IELEM**3+1))/3 + 3**(0.5D0)*Q2(22,IELEM**3+1) 
     >   + (15**(0.5D0)*Q2(10,IELEM**3+1))/3 - 5**(0.5D0)*Q2(23,IELEM*
     >   *3+1) + (21**(0.5D0)*Q2(08,IELEM**3+1))/3 - (5*3**(0.5D0)*Q2(
     >   27,IELEM**3+1))/3 + 5**(0.5D0)*Q2(26,IELEM**3+1) + 7**(0.5D0)
     >   *Q2(24,IELEM**3+1) + 15**(0.5D0)*Q2(19,IELEM**3+1) - (77*3**(
     >   0.5D0)*Q2(32,IELEM**3+1))/27 - (11*21**(0.5D0)*Q2(14,IELEM**3
     >   +1))/27 - (11*7**(0.5D0)*Q2(30,IELEM**3+1))/9 + 5**(0.5D0)*Q2
     >   (33,IELEM**3+1) - (15**(0.5D0)*Q2(25,IELEM**3+1))/3 - 21**(0.
     >   5D0)*Q2(20,IELEM**3+1) - (5*3**(0.5D0)*Q2(39,IELEM**3+1))/3 +
     >    5**(0.5D0)*Q2(38,IELEM**3+1) + (5*3**(0.5D0)*Q2(42,IELEM**3+
     >   1))/3 + (35**(0.5D0)*Q2(12,IELEM**3+1))/3 - (5*5**(0.5D0)*Q2(
     >   43,IELEM**3+1))/3 - 15**(0.5D0)*Q2(34,IELEM**3+1) + (11*21**(
     >   0.5D0)*Q2(29,IELEM**3+1))/27 + (11*35**(0.5D0)*Q2(15,IELEM**3
     >   +1))/27 + (5*7**(0.5D0)*Q2(44,IELEM**3+1))/3 - (15**(0.5D0)*Q
     >   2(37,IELEM**3+1))/3 - (77*5**(0.5D0)*Q2(48,IELEM**3+1))/27 + 
     >   (55*7**(0.5D0)*Q2(47,IELEM**3+1))/27 + 7**(0.5D0)*Q2(49,IELEM
     >   **3+1) + (7*3**(0.5D0)*Q2(56,IELEM**3+1))/3 + 7**(0.5D0)*Q2(5
     >   4,IELEM**3+1) - (77*3**(0.5D0)*Q2(62,IELEM**3+1))/27 + (7*5**
     >   (0.5D0)*Q2(60,IELEM**3+1))/3 - (5*7**(0.5D0)*Q2(59,IELEM**3+1
     >   ))/3 + (77*5**(0.5D0)*Q2(63,IELEM**3+1))/27 - (77*7**(0.5D0)*
     >   Q2(64,IELEM**3+1))/27 - 21**(0.5D0)*Q2(50,IELEM**3+1) - 35**(
     >   0.5D0)*Q2(36,IELEM**3+1) - (21**(0.5D0)*Q2(53,IELEM**3+1))/3 
     >   + (11*35**(0.5D0)*Q2(45,IELEM**3+1))/27 + 35**(0.5D0)*Q2(51,I
     >   ELEM**3+1) - (35**(0.5D0)*Q2(57,IELEM**3+1))/3 + (105**(0.5D0
     >   )*Q2(28,IELEM**3+1))/3 + (11*105**(0.5D0)*Q2(31,IELEM**3+1))/
     >   27 + (105**(0.5D0)*Q2(40,IELEM**3+1))/3 - (11*105**(0.5D0)*Q2
     >   (46,IELEM**3+1))/27 - (105**(0.5D0)*Q2(55,IELEM**3+1))/3 + (1
     >   05**(0.5D0)*Q2(58,IELEM**3+1))/3) 
      CORNERQ(06,4) = (Q2(01,IELEM**3+1) + Q2(06,IELEM**3+1)/3 + (5
     >   *Q2(11,IELEM**3+1))/9 + (847*Q2(16,IELEM**3+1))/729 - Q2(18,I
     >   ELEM**3+1) - Q2(21,IELEM**3+1) - (5*Q2(35,IELEM**3+1))/3 - (5
     >   *Q2(41,IELEM**3+1))/3 + (77*Q2(52,IELEM**3+1))/27 + (77*Q2(61
     >   ,IELEM**3+1))/27 - (3**(0.5D0)*Q2(02,IELEM**3+1))/3 - (3**(0.
     >   5D0)*Q2(05,IELEM**3+1))/3 - (5**(0.5D0)*Q2(03,IELEM**3+1))/3 
     >   + (11*7**(0.5D0)*Q2(04,IELEM**3+1))/27 - (5**(0.5D0)*Q2(09,IE
     >   LEM**3+1))/3 + 3**(0.5D0)*Q2(17,IELEM**3+1) + (11*7**(0.5D0)*
     >   Q2(13,IELEM**3+1))/27 + (15**(0.5D0)*Q2(07,IELEM**3+1))/9 + (
     >   3**(0.5D0)*Q2(22,IELEM**3+1))/3 + (15**(0.5D0)*Q2(10,IELEM**3
     >   +1))/9 + (5**(0.5D0)*Q2(23,IELEM**3+1))/3 - (11*21**(0.5D0)*Q
     >   2(08,IELEM**3+1))/81 + (5*3**(0.5D0)*Q2(27,IELEM**3+1))/9 + (
     >   5**(0.5D0)*Q2(26,IELEM**3+1))/3 - (11*7**(0.5D0)*Q2(24,IELEM*
     >   *3+1))/27 - (15**(0.5D0)*Q2(19,IELEM**3+1))/3 + (847*3**(0.5D
     >   0)*Q2(32,IELEM**3+1))/729 - (11*21**(0.5D0)*Q2(14,IELEM**3+1)
     >   )/81 - (11*7**(0.5D0)*Q2(30,IELEM**3+1))/27 + 5**(0.5D0)*Q2(3
     >   3,IELEM**3+1) - (15**(0.5D0)*Q2(25,IELEM**3+1))/3 + (11*21**(
     >   0.5D0)*Q2(20,IELEM**3+1))/27 + (5*3**(0.5D0)*Q2(39,IELEM**3+1
     >   ))/9 + (5**(0.5D0)*Q2(38,IELEM**3+1))/3 + (5*3**(0.5D0)*Q2(42
     >   ,IELEM**3+1))/9 - (11*35**(0.5D0)*Q2(12,IELEM**3+1))/81 + (5*
     >   5**(0.5D0)*Q2(43,IELEM**3+1))/9 - (15**(0.5D0)*Q2(34,IELEM**3
     >   +1))/3 + (11*21**(0.5D0)*Q2(29,IELEM**3+1))/27 - (11*35**(0.5
     >   D0)*Q2(15,IELEM**3+1))/81 - (55*7**(0.5D0)*Q2(44,IELEM**3+1))
     >   /81 - (15**(0.5D0)*Q2(37,IELEM**3+1))/3 + (847*5**(0.5D0)*Q2(
     >   48,IELEM**3+1))/729 - (55*7**(0.5D0)*Q2(47,IELEM**3+1))/81 + 
     >   7**(0.5D0)*Q2(49,IELEM**3+1) - (77*3**(0.5D0)*Q2(56,IELEM**3+
     >   1))/81 + (7**(0.5D0)*Q2(54,IELEM**3+1))/3 - (77*3**(0.5D0)*Q2
     >   (62,IELEM**3+1))/81 - (77*5**(0.5D0)*Q2(60,IELEM**3+1))/81 + 
     >   (5*7**(0.5D0)*Q2(59,IELEM**3+1))/9 - (77*5**(0.5D0)*Q2(63,IEL
     >   EM**3+1))/81 + (847*7**(0.5D0)*Q2(64,IELEM**3+1))/729 - (21**
     >   (0.5D0)*Q2(50,IELEM**3+1))/3 + (11*35**(0.5D0)*Q2(36,IELEM**3
     >   +1))/27 - (21**(0.5D0)*Q2(53,IELEM**3+1))/3 + (11*35**(0.5D0)
     >   *Q2(45,IELEM**3+1))/27 - (35**(0.5D0)*Q2(51,IELEM**3+1))/3 - 
     >   (35**(0.5D0)*Q2(57,IELEM**3+1))/3 - (11*105**(0.5D0)*Q2(28,IE
     >   LEM**3+1))/81 - (11*105**(0.5D0)*Q2(31,IELEM**3+1))/81 - (11*
     >   105**(0.5D0)*Q2(40,IELEM**3+1))/81 - (11*105**(0.5D0)*Q2(46,I
     >   ELEM**3+1))/81 + (105**(0.5D0)*Q2(55,IELEM**3+1))/9 + (105**(
     >   0.5D0)*Q2(58,IELEM**3+1))/9) 
      CORNERQ(07,4) = (Q2(01,IELEM**3+1) - Q2(06,IELEM**3+1)/3 + (5
     >   *Q2(11,IELEM**3+1))/9 - (847*Q2(16,IELEM**3+1))/729 + Q2(18,I
     >   ELEM**3+1) - Q2(21,IELEM**3+1) - (5*Q2(35,IELEM**3+1))/3 - (5
     >   *Q2(41,IELEM**3+1))/3 - (77*Q2(52,IELEM**3+1))/27 + (77*Q2(61
     >   ,IELEM**3+1))/27 + (3**(0.5D0)*Q2(02,IELEM**3+1))/3 - (3**(0.
     >   5D0)*Q2(05,IELEM**3+1))/3 - (5**(0.5D0)*Q2(03,IELEM**3+1))/3 
     >   - (11*7**(0.5D0)*Q2(04,IELEM**3+1))/27 - (5**(0.5D0)*Q2(09,IE
     >   LEM**3+1))/3 + 3**(0.5D0)*Q2(17,IELEM**3+1) + (11*7**(0.5D0)*
     >   Q2(13,IELEM**3+1))/27 + (15**(0.5D0)*Q2(07,IELEM**3+1))/9 - (
     >   3**(0.5D0)*Q2(22,IELEM**3+1))/3 - (15**(0.5D0)*Q2(10,IELEM**3
     >   +1))/9 + (5**(0.5D0)*Q2(23,IELEM**3+1))/3 + (11*21**(0.5D0)*Q
     >   2(08,IELEM**3+1))/81 + (5*3**(0.5D0)*Q2(27,IELEM**3+1))/9 - (
     >   5**(0.5D0)*Q2(26,IELEM**3+1))/3 + (11*7**(0.5D0)*Q2(24,IELEM*
     >   *3+1))/27 - (15**(0.5D0)*Q2(19,IELEM**3+1))/3 - (847*3**(0.5D
     >   0)*Q2(32,IELEM**3+1))/729 + (11*21**(0.5D0)*Q2(14,IELEM**3+1)
     >   )/81 + (11*7**(0.5D0)*Q2(30,IELEM**3+1))/27 + 5**(0.5D0)*Q2(3
     >   3,IELEM**3+1) - (15**(0.5D0)*Q2(25,IELEM**3+1))/3 - (11*21**(
     >   0.5D0)*Q2(20,IELEM**3+1))/27 + (5*3**(0.5D0)*Q2(39,IELEM**3+1
     >   ))/9 - (5**(0.5D0)*Q2(38,IELEM**3+1))/3 - (5*3**(0.5D0)*Q2(42
     >   ,IELEM**3+1))/9 + (11*35**(0.5D0)*Q2(12,IELEM**3+1))/81 + (5*
     >   5**(0.5D0)*Q2(43,IELEM**3+1))/9 + (15**(0.5D0)*Q2(34,IELEM**3
     >   +1))/3 + (11*21**(0.5D0)*Q2(29,IELEM**3+1))/27 - (11*35**(0.5
     >   D0)*Q2(15,IELEM**3+1))/81 + (55*7**(0.5D0)*Q2(44,IELEM**3+1))
     >   /81 - (15**(0.5D0)*Q2(37,IELEM**3+1))/3 - (847*5**(0.5D0)*Q2(
     >   48,IELEM**3+1))/729 - (55*7**(0.5D0)*Q2(47,IELEM**3+1))/81 + 
     >   7**(0.5D0)*Q2(49,IELEM**3+1) + (77*3**(0.5D0)*Q2(56,IELEM**3+
     >   1))/81 - (7**(0.5D0)*Q2(54,IELEM**3+1))/3 + (77*3**(0.5D0)*Q2
     >   (62,IELEM**3+1))/81 + (77*5**(0.5D0)*Q2(60,IELEM**3+1))/81 + 
     >   (5*7**(0.5D0)*Q2(59,IELEM**3+1))/9 - (77*5**(0.5D0)*Q2(63,IEL
     >   EM**3+1))/81 - (847*7**(0.5D0)*Q2(64,IELEM**3+1))/729 + (21**
     >   (0.5D0)*Q2(50,IELEM**3+1))/3 - (11*35**(0.5D0)*Q2(36,IELEM**3
     >   +1))/27 - (21**(0.5D0)*Q2(53,IELEM**3+1))/3 + (11*35**(0.5D0)
     >   *Q2(45,IELEM**3+1))/27 - (35**(0.5D0)*Q2(51,IELEM**3+1))/3 - 
     >   (35**(0.5D0)*Q2(57,IELEM**3+1))/3 + (11*105**(0.5D0)*Q2(28,IE
     >   LEM**3+1))/81 - (11*105**(0.5D0)*Q2(31,IELEM**3+1))/81 + (11*
     >   105**(0.5D0)*Q2(40,IELEM**3+1))/81 + (11*105**(0.5D0)*Q2(46,I
     >   ELEM**3+1))/81 + (105**(0.5D0)*Q2(55,IELEM**3+1))/9 - (105**(
     >   0.5D0)*Q2(58,IELEM**3+1))/9) 
      CORNERQ(08,4) = (Q2(01,IELEM**3+1) - Q2(06,IELEM**3+1) - (5*Q
     >   2(11,IELEM**3+1))/3 + (77*Q2(16,IELEM**3+1))/27 + 3*Q2(18,IEL
     >   EM**3+1) - Q2(21,IELEM**3+1) + 5*Q2(35,IELEM**3+1) - (5*Q2(41
     >   ,IELEM**3+1))/3 + 7*Q2(52,IELEM**3+1) + (77*Q2(61,IELEM**3+1)
     >   )/27 + 3**(0.5D0)*Q2(02,IELEM**3+1) - (3**(0.5D0)*Q2(05,IELEM
     >   **3+1))/3 + 5**(0.5D0)*Q2(03,IELEM**3+1) + 7**(0.5D0)*Q2(04,I
     >   ELEM**3+1) - (5**(0.5D0)*Q2(09,IELEM**3+1))/3 + 3**(0.5D0)*Q2
     >   (17,IELEM**3+1) + (11*7**(0.5D0)*Q2(13,IELEM**3+1))/27 - (15*
     >   *(0.5D0)*Q2(07,IELEM**3+1))/3 - 3**(0.5D0)*Q2(22,IELEM**3+1) 
     >   - (15**(0.5D0)*Q2(10,IELEM**3+1))/3 - 5**(0.5D0)*Q2(23,IELEM*
     >   *3+1) - (21**(0.5D0)*Q2(08,IELEM**3+1))/3 - (5*3**(0.5D0)*Q2(
     >   27,IELEM**3+1))/3 - 5**(0.5D0)*Q2(26,IELEM**3+1) - 7**(0.5D0)
     >   *Q2(24,IELEM**3+1) + 15**(0.5D0)*Q2(19,IELEM**3+1) + (77*3**(
     >   0.5D0)*Q2(32,IELEM**3+1))/27 + (11*21**(0.5D0)*Q2(14,IELEM**3
     >   +1))/27 + (11*7**(0.5D0)*Q2(30,IELEM**3+1))/9 + 5**(0.5D0)*Q2
     >   (33,IELEM**3+1) - (15**(0.5D0)*Q2(25,IELEM**3+1))/3 + 21**(0.
     >   5D0)*Q2(20,IELEM**3+1) - (5*3**(0.5D0)*Q2(39,IELEM**3+1))/3 -
     >    5**(0.5D0)*Q2(38,IELEM**3+1) - (5*3**(0.5D0)*Q2(42,IELEM**3+
     >   1))/3 - (35**(0.5D0)*Q2(12,IELEM**3+1))/3 - (5*5**(0.5D0)*Q2(
     >   43,IELEM**3+1))/3 + 15**(0.5D0)*Q2(34,IELEM**3+1) + (11*21**(
     >   0.5D0)*Q2(29,IELEM**3+1))/27 + (11*35**(0.5D0)*Q2(15,IELEM**3
     >   +1))/27 - (5*7**(0.5D0)*Q2(44,IELEM**3+1))/3 - (15**(0.5D0)*Q
     >   2(37,IELEM**3+1))/3 + (77*5**(0.5D0)*Q2(48,IELEM**3+1))/27 + 
     >   (55*7**(0.5D0)*Q2(47,IELEM**3+1))/27 + 7**(0.5D0)*Q2(49,IELEM
     >   **3+1) - (7*3**(0.5D0)*Q2(56,IELEM**3+1))/3 - 7**(0.5D0)*Q2(5
     >   4,IELEM**3+1) + (77*3**(0.5D0)*Q2(62,IELEM**3+1))/27 - (7*5**
     >   (0.5D0)*Q2(60,IELEM**3+1))/3 - (5*7**(0.5D0)*Q2(59,IELEM**3+1
     >   ))/3 + (77*5**(0.5D0)*Q2(63,IELEM**3+1))/27 + (77*7**(0.5D0)*
     >   Q2(64,IELEM**3+1))/27 + 21**(0.5D0)*Q2(50,IELEM**3+1) + 35**(
     >   0.5D0)*Q2(36,IELEM**3+1) - (21**(0.5D0)*Q2(53,IELEM**3+1))/3 
     >   + (11*35**(0.5D0)*Q2(45,IELEM**3+1))/27 + 35**(0.5D0)*Q2(51,I
     >   ELEM**3+1) - (35**(0.5D0)*Q2(57,IELEM**3+1))/3 - (105**(0.5D0
     >   )*Q2(28,IELEM**3+1))/3 + (11*105**(0.5D0)*Q2(31,IELEM**3+1))/
     >   27 - (105**(0.5D0)*Q2(40,IELEM**3+1))/3 + (11*105**(0.5D0)*Q2
     >   (46,IELEM**3+1))/27 - (105**(0.5D0)*Q2(55,IELEM**3+1))/3 - (1
     >   05**(0.5D0)*Q2(58,IELEM**3+1))/3) 
      CORNERQ(09,4) = (Q2(01,IELEM**3+1) - Q2(06,IELEM**3+1) - (5*Q
     >   2(11,IELEM**3+1))/3 + (77*Q2(16,IELEM**3+1))/27 - 3*Q2(18,IEL
     >   EM**3+1) + Q2(21,IELEM**3+1) + 5*Q2(35,IELEM**3+1) - (5*Q2(41
     >   ,IELEM**3+1))/3 - 7*Q2(52,IELEM**3+1) - (77*Q2(61,IELEM**3+1)
     >   )/27 - 3**(0.5D0)*Q2(02,IELEM**3+1) + (3**(0.5D0)*Q2(05,IELEM
     >   **3+1))/3 + 5**(0.5D0)*Q2(03,IELEM**3+1) - 7**(0.5D0)*Q2(04,I
     >   ELEM**3+1) - (5**(0.5D0)*Q2(09,IELEM**3+1))/3 + 3**(0.5D0)*Q2
     >   (17,IELEM**3+1) - (11*7**(0.5D0)*Q2(13,IELEM**3+1))/27 + (15*
     >   *(0.5D0)*Q2(07,IELEM**3+1))/3 - 3**(0.5D0)*Q2(22,IELEM**3+1) 
     >   + (15**(0.5D0)*Q2(10,IELEM**3+1))/3 + 5**(0.5D0)*Q2(23,IELEM*
     >   *3+1) - (21**(0.5D0)*Q2(08,IELEM**3+1))/3 - (5*3**(0.5D0)*Q2(
     >   27,IELEM**3+1))/3 + 5**(0.5D0)*Q2(26,IELEM**3+1) - 7**(0.5D0)
     >   *Q2(24,IELEM**3+1) + 15**(0.5D0)*Q2(19,IELEM**3+1) + (77*3**(
     >   0.5D0)*Q2(32,IELEM**3+1))/27 + (11*21**(0.5D0)*Q2(14,IELEM**3
     >   +1))/27 + (11*7**(0.5D0)*Q2(30,IELEM**3+1))/9 + 5**(0.5D0)*Q2
     >   (33,IELEM**3+1) - (15**(0.5D0)*Q2(25,IELEM**3+1))/3 - 21**(0.
     >   5D0)*Q2(20,IELEM**3+1) + (5*3**(0.5D0)*Q2(39,IELEM**3+1))/3 -
     >    5**(0.5D0)*Q2(38,IELEM**3+1) + (5*3**(0.5D0)*Q2(42,IELEM**3+
     >   1))/3 + (35**(0.5D0)*Q2(12,IELEM**3+1))/3 - (5*5**(0.5D0)*Q2(
     >   43,IELEM**3+1))/3 - 15**(0.5D0)*Q2(34,IELEM**3+1) - (11*21**(
     >   0.5D0)*Q2(29,IELEM**3+1))/27 - (11*35**(0.5D0)*Q2(15,IELEM**3
     >   +1))/27 + (5*7**(0.5D0)*Q2(44,IELEM**3+1))/3 + (15**(0.5D0)*Q
     >   2(37,IELEM**3+1))/3 + (77*5**(0.5D0)*Q2(48,IELEM**3+1))/27 - 
     >   (55*7**(0.5D0)*Q2(47,IELEM**3+1))/27 + 7**(0.5D0)*Q2(49,IELEM
     >   **3+1) - (7*3**(0.5D0)*Q2(56,IELEM**3+1))/3 - 7**(0.5D0)*Q2(5
     >   4,IELEM**3+1) + (77*3**(0.5D0)*Q2(62,IELEM**3+1))/27 + (7*5**
     >   (0.5D0)*Q2(60,IELEM**3+1))/3 - (5*7**(0.5D0)*Q2(59,IELEM**3+1
     >   ))/3 - (77*5**(0.5D0)*Q2(63,IELEM**3+1))/27 + (77*7**(0.5D0)*
     >   Q2(64,IELEM**3+1))/27 - 21**(0.5D0)*Q2(50,IELEM**3+1) - 35**(
     >   0.5D0)*Q2(36,IELEM**3+1) + (21**(0.5D0)*Q2(53,IELEM**3+1))/3 
     >   - (11*35**(0.5D0)*Q2(45,IELEM**3+1))/27 + 35**(0.5D0)*Q2(51,I
     >   ELEM**3+1) - (35**(0.5D0)*Q2(57,IELEM**3+1))/3 + (105**(0.5D0
     >   )*Q2(28,IELEM**3+1))/3 - (11*105**(0.5D0)*Q2(31,IELEM**3+1))/
     >   27 - (105**(0.5D0)*Q2(40,IELEM**3+1))/3 + (11*105**(0.5D0)*Q2
     >   (46,IELEM**3+1))/27 + (105**(0.5D0)*Q2(55,IELEM**3+1))/3 + (1
     >   05**(0.5D0)*Q2(58,IELEM**3+1))/3) 
      CORNERQ(10,4) = (Q2(01,IELEM**3+1) - Q2(06,IELEM**3+1)/3 + (5
     >   *Q2(11,IELEM**3+1))/9 - (847*Q2(16,IELEM**3+1))/729 - Q2(18,I
     >   ELEM**3+1) + Q2(21,IELEM**3+1) - (5*Q2(35,IELEM**3+1))/3 - (5
     >   *Q2(41,IELEM**3+1))/3 + (77*Q2(52,IELEM**3+1))/27 - (77*Q2(61
     >   ,IELEM**3+1))/27 - (3**(0.5D0)*Q2(02,IELEM**3+1))/3 + (3**(0.
     >   5D0)*Q2(05,IELEM**3+1))/3 - (5**(0.5D0)*Q2(03,IELEM**3+1))/3 
     >   + (11*7**(0.5D0)*Q2(04,IELEM**3+1))/27 - (5**(0.5D0)*Q2(09,IE
     >   LEM**3+1))/3 + 3**(0.5D0)*Q2(17,IELEM**3+1) - (11*7**(0.5D0)*
     >   Q2(13,IELEM**3+1))/27 - (15**(0.5D0)*Q2(07,IELEM**3+1))/9 - (
     >   3**(0.5D0)*Q2(22,IELEM**3+1))/3 + (15**(0.5D0)*Q2(10,IELEM**3
     >   +1))/9 - (5**(0.5D0)*Q2(23,IELEM**3+1))/3 + (11*21**(0.5D0)*Q
     >   2(08,IELEM**3+1))/81 + (5*3**(0.5D0)*Q2(27,IELEM**3+1))/9 + (
     >   5**(0.5D0)*Q2(26,IELEM**3+1))/3 + (11*7**(0.5D0)*Q2(24,IELEM*
     >   *3+1))/27 - (15**(0.5D0)*Q2(19,IELEM**3+1))/3 - (847*3**(0.5D
     >   0)*Q2(32,IELEM**3+1))/729 + (11*21**(0.5D0)*Q2(14,IELEM**3+1)
     >   )/81 + (11*7**(0.5D0)*Q2(30,IELEM**3+1))/27 + 5**(0.5D0)*Q2(3
     >   3,IELEM**3+1) - (15**(0.5D0)*Q2(25,IELEM**3+1))/3 + (11*21**(
     >   0.5D0)*Q2(20,IELEM**3+1))/27 - (5*3**(0.5D0)*Q2(39,IELEM**3+1
     >   ))/9 - (5**(0.5D0)*Q2(38,IELEM**3+1))/3 + (5*3**(0.5D0)*Q2(42
     >   ,IELEM**3+1))/9 - (11*35**(0.5D0)*Q2(12,IELEM**3+1))/81 + (5*
     >   5**(0.5D0)*Q2(43,IELEM**3+1))/9 - (15**(0.5D0)*Q2(34,IELEM**3
     >   +1))/3 - (11*21**(0.5D0)*Q2(29,IELEM**3+1))/27 + (11*35**(0.5
     >   D0)*Q2(15,IELEM**3+1))/81 - (55*7**(0.5D0)*Q2(44,IELEM**3+1))
     >   /81 + (15**(0.5D0)*Q2(37,IELEM**3+1))/3 - (847*5**(0.5D0)*Q2(
     >   48,IELEM**3+1))/729 + (55*7**(0.5D0)*Q2(47,IELEM**3+1))/81 + 
     >   7**(0.5D0)*Q2(49,IELEM**3+1) + (77*3**(0.5D0)*Q2(56,IELEM**3+
     >   1))/81 - (7**(0.5D0)*Q2(54,IELEM**3+1))/3 + (77*3**(0.5D0)*Q2
     >   (62,IELEM**3+1))/81 - (77*5**(0.5D0)*Q2(60,IELEM**3+1))/81 + 
     >   (5*7**(0.5D0)*Q2(59,IELEM**3+1))/9 + (77*5**(0.5D0)*Q2(63,IEL
     >   EM**3+1))/81 - (847*7**(0.5D0)*Q2(64,IELEM**3+1))/729 - (21**
     >   (0.5D0)*Q2(50,IELEM**3+1))/3 + (11*35**(0.5D0)*Q2(36,IELEM**3
     >   +1))/27 + (21**(0.5D0)*Q2(53,IELEM**3+1))/3 - (11*35**(0.5D0)
     >   *Q2(45,IELEM**3+1))/27 - (35**(0.5D0)*Q2(51,IELEM**3+1))/3 - 
     >   (35**(0.5D0)*Q2(57,IELEM**3+1))/3 - (11*105**(0.5D0)*Q2(28,IE
     >   LEM**3+1))/81 + (11*105**(0.5D0)*Q2(31,IELEM**3+1))/81 + (11*
     >   105**(0.5D0)*Q2(40,IELEM**3+1))/81 + (11*105**(0.5D0)*Q2(46,I
     >   ELEM**3+1))/81 - (105**(0.5D0)*Q2(55,IELEM**3+1))/9 + (105**(
     >   0.5D0)*Q2(58,IELEM**3+1))/9) 
      CORNERQ(11,4) = (Q2(01,IELEM**3+1) + Q2(06,IELEM**3+1)/3 + (5
     >   *Q2(11,IELEM**3+1))/9 + (847*Q2(16,IELEM**3+1))/729 + Q2(18,I
     >   ELEM**3+1) + Q2(21,IELEM**3+1) - (5*Q2(35,IELEM**3+1))/3 - (5
     >   *Q2(41,IELEM**3+1))/3 - (77*Q2(52,IELEM**3+1))/27 - (77*Q2(61
     >   ,IELEM**3+1))/27 + (3**(0.5D0)*Q2(02,IELEM**3+1))/3 + (3**(0.
     >   5D0)*Q2(05,IELEM**3+1))/3 - (5**(0.5D0)*Q2(03,IELEM**3+1))/3 
     >   - (11*7**(0.5D0)*Q2(04,IELEM**3+1))/27 - (5**(0.5D0)*Q2(09,IE
     >   LEM**3+1))/3 + 3**(0.5D0)*Q2(17,IELEM**3+1) - (11*7**(0.5D0)*
     >   Q2(13,IELEM**3+1))/27 - (15**(0.5D0)*Q2(07,IELEM**3+1))/9 + (
     >   3**(0.5D0)*Q2(22,IELEM**3+1))/3 - (15**(0.5D0)*Q2(10,IELEM**3
     >   +1))/9 - (5**(0.5D0)*Q2(23,IELEM**3+1))/3 - (11*21**(0.5D0)*Q
     >   2(08,IELEM**3+1))/81 + (5*3**(0.5D0)*Q2(27,IELEM**3+1))/9 - (
     >   5**(0.5D0)*Q2(26,IELEM**3+1))/3 - (11*7**(0.5D0)*Q2(24,IELEM*
     >   *3+1))/27 - (15**(0.5D0)*Q2(19,IELEM**3+1))/3 + (847*3**(0.5D
     >   0)*Q2(32,IELEM**3+1))/729 - (11*21**(0.5D0)*Q2(14,IELEM**3+1)
     >   )/81 - (11*7**(0.5D0)*Q2(30,IELEM**3+1))/27 + 5**(0.5D0)*Q2(3
     >   3,IELEM**3+1) - (15**(0.5D0)*Q2(25,IELEM**3+1))/3 - (11*21**(
     >   0.5D0)*Q2(20,IELEM**3+1))/27 - (5*3**(0.5D0)*Q2(39,IELEM**3+1
     >   ))/9 + (5**(0.5D0)*Q2(38,IELEM**3+1))/3 - (5*3**(0.5D0)*Q2(42
     >   ,IELEM**3+1))/9 + (11*35**(0.5D0)*Q2(12,IELEM**3+1))/81 + (5*
     >   5**(0.5D0)*Q2(43,IELEM**3+1))/9 + (15**(0.5D0)*Q2(34,IELEM**3
     >   +1))/3 - (11*21**(0.5D0)*Q2(29,IELEM**3+1))/27 + (11*35**(0.5
     >   D0)*Q2(15,IELEM**3+1))/81 + (55*7**(0.5D0)*Q2(44,IELEM**3+1))
     >   /81 + (15**(0.5D0)*Q2(37,IELEM**3+1))/3 + (847*5**(0.5D0)*Q2(
     >   48,IELEM**3+1))/729 + (55*7**(0.5D0)*Q2(47,IELEM**3+1))/81 + 
     >   7**(0.5D0)*Q2(49,IELEM**3+1) - (77*3**(0.5D0)*Q2(56,IELEM**3+
     >   1))/81 + (7**(0.5D0)*Q2(54,IELEM**3+1))/3 - (77*3**(0.5D0)*Q2
     >   (62,IELEM**3+1))/81 + (77*5**(0.5D0)*Q2(60,IELEM**3+1))/81 + 
     >   (5*7**(0.5D0)*Q2(59,IELEM**3+1))/9 + (77*5**(0.5D0)*Q2(63,IEL
     >   EM**3+1))/81 + (847*7**(0.5D0)*Q2(64,IELEM**3+1))/729 + (21**
     >   (0.5D0)*Q2(50,IELEM**3+1))/3 - (11*35**(0.5D0)*Q2(36,IELEM**3
     >   +1))/27 + (21**(0.5D0)*Q2(53,IELEM**3+1))/3 - (11*35**(0.5D0)
     >   *Q2(45,IELEM**3+1))/27 - (35**(0.5D0)*Q2(51,IELEM**3+1))/3 - 
     >   (35**(0.5D0)*Q2(57,IELEM**3+1))/3 + (11*105**(0.5D0)*Q2(28,IE
     >   LEM**3+1))/81 + (11*105**(0.5D0)*Q2(31,IELEM**3+1))/81 - (11*
     >   105**(0.5D0)*Q2(40,IELEM**3+1))/81 - (11*105**(0.5D0)*Q2(46,I
     >   ELEM**3+1))/81 - (105**(0.5D0)*Q2(55,IELEM**3+1))/9 - (105**(
     >   0.5D0)*Q2(58,IELEM**3+1))/9) 
      CORNERQ(12,4) = (Q2(01,IELEM**3+1) + Q2(06,IELEM**3+1) - (5*Q
     >   2(11,IELEM**3+1))/3 - (77*Q2(16,IELEM**3+1))/27 + 3*Q2(18,IEL
     >   EM**3+1) + Q2(21,IELEM**3+1) + 5*Q2(35,IELEM**3+1) - (5*Q2(41
     >   ,IELEM**3+1))/3 + 7*Q2(52,IELEM**3+1) - (77*Q2(61,IELEM**3+1)
     >   )/27 + 3**(0.5D0)*Q2(02,IELEM**3+1) + (3**(0.5D0)*Q2(05,IELEM
     >   **3+1))/3 + 5**(0.5D0)*Q2(03,IELEM**3+1) + 7**(0.5D0)*Q2(04,I
     >   ELEM**3+1) - (5**(0.5D0)*Q2(09,IELEM**3+1))/3 + 3**(0.5D0)*Q2
     >   (17,IELEM**3+1) - (11*7**(0.5D0)*Q2(13,IELEM**3+1))/27 + (15*
     >   *(0.5D0)*Q2(07,IELEM**3+1))/3 + 3**(0.5D0)*Q2(22,IELEM**3+1) 
     >   - (15**(0.5D0)*Q2(10,IELEM**3+1))/3 + 5**(0.5D0)*Q2(23,IELEM*
     >   *3+1) + (21**(0.5D0)*Q2(08,IELEM**3+1))/3 - (5*3**(0.5D0)*Q2(
     >   27,IELEM**3+1))/3 - 5**(0.5D0)*Q2(26,IELEM**3+1) + 7**(0.5D0)
     >   *Q2(24,IELEM**3+1) + 15**(0.5D0)*Q2(19,IELEM**3+1) - (77*3**(
     >   0.5D0)*Q2(32,IELEM**3+1))/27 - (11*21**(0.5D0)*Q2(14,IELEM**3
     >   +1))/27 - (11*7**(0.5D0)*Q2(30,IELEM**3+1))/9 + 5**(0.5D0)*Q2
     >   (33,IELEM**3+1) - (15**(0.5D0)*Q2(25,IELEM**3+1))/3 + 21**(0.
     >   5D0)*Q2(20,IELEM**3+1) + (5*3**(0.5D0)*Q2(39,IELEM**3+1))/3 +
     >    5**(0.5D0)*Q2(38,IELEM**3+1) - (5*3**(0.5D0)*Q2(42,IELEM**3+
     >   1))/3 - (35**(0.5D0)*Q2(12,IELEM**3+1))/3 - (5*5**(0.5D0)*Q2(
     >   43,IELEM**3+1))/3 + 15**(0.5D0)*Q2(34,IELEM**3+1) - (11*21**(
     >   0.5D0)*Q2(29,IELEM**3+1))/27 - (11*35**(0.5D0)*Q2(15,IELEM**3
     >   +1))/27 - (5*7**(0.5D0)*Q2(44,IELEM**3+1))/3 + (15**(0.5D0)*Q
     >   2(37,IELEM**3+1))/3 - (77*5**(0.5D0)*Q2(48,IELEM**3+1))/27 - 
     >   (55*7**(0.5D0)*Q2(47,IELEM**3+1))/27 + 7**(0.5D0)*Q2(49,IELEM
     >   **3+1) + (7*3**(0.5D0)*Q2(56,IELEM**3+1))/3 + 7**(0.5D0)*Q2(5
     >   4,IELEM**3+1) - (77*3**(0.5D0)*Q2(62,IELEM**3+1))/27 - (7*5**
     >   (0.5D0)*Q2(60,IELEM**3+1))/3 - (5*7**(0.5D0)*Q2(59,IELEM**3+1
     >   ))/3 - (77*5**(0.5D0)*Q2(63,IELEM**3+1))/27 - (77*7**(0.5D0)*
     >   Q2(64,IELEM**3+1))/27 + 21**(0.5D0)*Q2(50,IELEM**3+1) + 35**(
     >   0.5D0)*Q2(36,IELEM**3+1) + (21**(0.5D0)*Q2(53,IELEM**3+1))/3 
     >   - (11*35**(0.5D0)*Q2(45,IELEM**3+1))/27 + 35**(0.5D0)*Q2(51,I
     >   ELEM**3+1) - (35**(0.5D0)*Q2(57,IELEM**3+1))/3 - (105**(0.5D0
     >   )*Q2(28,IELEM**3+1))/3 - (11*105**(0.5D0)*Q2(31,IELEM**3+1))/
     >   27 + (105**(0.5D0)*Q2(40,IELEM**3+1))/3 - (11*105**(0.5D0)*Q2
     >   (46,IELEM**3+1))/27 + (105**(0.5D0)*Q2(55,IELEM**3+1))/3 - (1
     >   05**(0.5D0)*Q2(58,IELEM**3+1))/3) 
      CORNERQ(13,4) = (Q2(01,IELEM**3+1) - 3*Q2(06,IELEM**3+1) + 5*
     >   Q2(11,IELEM**3+1) - 7*Q2(16,IELEM**3+1) - 3*Q2(18,IELEM**3+1)
     >    + 3*Q2(21,IELEM**3+1) + 5*Q2(35,IELEM**3+1) + 5*Q2(41,IELEM*
     >   *3+1) - 7*Q2(52,IELEM**3+1) + 7*Q2(61,IELEM**3+1) - 3**(0.5D0
     >   )*Q2(02,IELEM**3+1) + 3**(0.5D0)*Q2(05,IELEM**3+1) + 5**(0.5D
     >   0)*Q2(03,IELEM**3+1) - 7**(0.5D0)*Q2(04,IELEM**3+1) + 5**(0.5
     >   D0)*Q2(09,IELEM**3+1) + 3**(0.5D0)*Q2(17,IELEM**3+1) + 7**(0.
     >   5D0)*Q2(13,IELEM**3+1) + 15**(0.5D0)*Q2(07,IELEM**3+1) - 3*3*
     >   *(0.5D0)*Q2(22,IELEM**3+1) - 15**(0.5D0)*Q2(10,IELEM**3+1) + 
     >   3*5**(0.5D0)*Q2(23,IELEM**3+1) - 21**(0.5D0)*Q2(08,IELEM**3+1
     >   ) + 5*3**(0.5D0)*Q2(27,IELEM**3+1) - 3*5**(0.5D0)*Q2(26,IELEM
     >   **3+1) - 3*7**(0.5D0)*Q2(24,IELEM**3+1) + 15**(0.5D0)*Q2(19,I
     >   ELEM**3+1) - 7*3**(0.5D0)*Q2(32,IELEM**3+1) - 21**(0.5D0)*Q2(
     >   14,IELEM**3+1) - 3*7**(0.5D0)*Q2(30,IELEM**3+1) + 5**(0.5D0)*
     >   Q2(33,IELEM**3+1) + 15**(0.5D0)*Q2(25,IELEM**3+1) - 21**(0.5D
     >   0)*Q2(20,IELEM**3+1) + 5*3**(0.5D0)*Q2(39,IELEM**3+1) - 3*5**
     >   (0.5D0)*Q2(38,IELEM**3+1) - 5*3**(0.5D0)*Q2(42,IELEM**3+1) - 
     >   35**(0.5D0)*Q2(12,IELEM**3+1) + 5*5**(0.5D0)*Q2(43,IELEM**3+1
     >   ) - 15**(0.5D0)*Q2(34,IELEM**3+1) + 21**(0.5D0)*Q2(29,IELEM**
     >   3+1) + 35**(0.5D0)*Q2(15,IELEM**3+1) - 5*7**(0.5D0)*Q2(44,IEL
     >   EM**3+1) + 15**(0.5D0)*Q2(37,IELEM**3+1) - 7*5**(0.5D0)*Q2(48
     >   ,IELEM**3+1) + 5*7**(0.5D0)*Q2(47,IELEM**3+1) + 7**(0.5D0)*Q2
     >   (49,IELEM**3+1) - 7*3**(0.5D0)*Q2(56,IELEM**3+1) - 3*7**(0.5D
     >   0)*Q2(54,IELEM**3+1) - 7*3**(0.5D0)*Q2(62,IELEM**3+1) - 7*5**
     >   (0.5D0)*Q2(60,IELEM**3+1) + 5*7**(0.5D0)*Q2(59,IELEM**3+1) + 
     >   7*5**(0.5D0)*Q2(63,IELEM**3+1) - 7*7**(0.5D0)*Q2(64,IELEM**3+
     >   1) - 21**(0.5D0)*Q2(50,IELEM**3+1) - 35**(0.5D0)*Q2(36,IELEM*
     >   *3+1) + 21**(0.5D0)*Q2(53,IELEM**3+1) + 35**(0.5D0)*Q2(45,IEL
     >   EM**3+1) + 35**(0.5D0)*Q2(51,IELEM**3+1) + 35**(0.5D0)*Q2(57,
     >   IELEM**3+1) - 105**(0.5D0)*Q2(28,IELEM**3+1) + 105**(0.5D0)*Q
     >   2(31,IELEM**3+1) - 105**(0.5D0)*Q2(40,IELEM**3+1) - 105**(0.5
     >   D0)*Q2(46,IELEM**3+1) + 105**(0.5D0)*Q2(55,IELEM**3+1) - 105*
     >   *(0.5D0)*Q2(58,IELEM**3+1)) 
      CORNERQ(14,4) = (Q2(01,IELEM**3+1) - Q2(06,IELEM**3+1) - (5*Q
     >   2(11,IELEM**3+1))/3 + (77*Q2(16,IELEM**3+1))/27 - Q2(18,IELEM
     >   **3+1) + 3*Q2(21,IELEM**3+1) - (5*Q2(35,IELEM**3+1))/3 + 5*Q2
     >   (41,IELEM**3+1) + (77*Q2(52,IELEM**3+1))/27 + 7*Q2(61,IELEM**
     >   3+1) - (3**(0.5D0)*Q2(02,IELEM**3+1))/3 + 3**(0.5D0)*Q2(05,IE
     >   LEM**3+1) - (5**(0.5D0)*Q2(03,IELEM**3+1))/3 + (11*7**(0.5D0)
     >   *Q2(04,IELEM**3+1))/27 + 5**(0.5D0)*Q2(09,IELEM**3+1) + 3**(0
     >   .5D0)*Q2(17,IELEM**3+1) + 7**(0.5D0)*Q2(13,IELEM**3+1) - (15*
     >   *(0.5D0)*Q2(07,IELEM**3+1))/3 - 3**(0.5D0)*Q2(22,IELEM**3+1) 
     >   - (15**(0.5D0)*Q2(10,IELEM**3+1))/3 - 5**(0.5D0)*Q2(23,IELEM*
     >   *3+1) + (11*21**(0.5D0)*Q2(08,IELEM**3+1))/27 - (5*3**(0.5D0)
     >   *Q2(27,IELEM**3+1))/3 - 5**(0.5D0)*Q2(26,IELEM**3+1) + (11*7*
     >   *(0.5D0)*Q2(24,IELEM**3+1))/9 - (15**(0.5D0)*Q2(19,IELEM**3+1
     >   ))/3 + (77*3**(0.5D0)*Q2(32,IELEM**3+1))/27 - (21**(0.5D0)*Q2
     >   (14,IELEM**3+1))/3 - 7**(0.5D0)*Q2(30,IELEM**3+1) + 5**(0.5D0
     >   )*Q2(33,IELEM**3+1) + 15**(0.5D0)*Q2(25,IELEM**3+1) + (11*21*
     >   *(0.5D0)*Q2(20,IELEM**3+1))/27 - (5*3**(0.5D0)*Q2(39,IELEM**3
     >   +1))/3 - 5**(0.5D0)*Q2(38,IELEM**3+1) - (5*3**(0.5D0)*Q2(42,I
     >   ELEM**3+1))/3 + (11*35**(0.5D0)*Q2(12,IELEM**3+1))/27 - (5*5*
     >   *(0.5D0)*Q2(43,IELEM**3+1))/3 - (15**(0.5D0)*Q2(34,IELEM**3+1
     >   ))/3 + 21**(0.5D0)*Q2(29,IELEM**3+1) - (35**(0.5D0)*Q2(15,IEL
     >   EM**3+1))/3 + (55*7**(0.5D0)*Q2(44,IELEM**3+1))/27 + 15**(0.5
     >   D0)*Q2(37,IELEM**3+1) + (77*5**(0.5D0)*Q2(48,IELEM**3+1))/27 
     >   - (5*7**(0.5D0)*Q2(47,IELEM**3+1))/3 + 7**(0.5D0)*Q2(49,IELEM
     >   **3+1) + (77*3**(0.5D0)*Q2(56,IELEM**3+1))/27 - 7**(0.5D0)*Q2
     >   (54,IELEM**3+1) - (7*3**(0.5D0)*Q2(62,IELEM**3+1))/3 + (77*5*
     >   *(0.5D0)*Q2(60,IELEM**3+1))/27 - (5*7**(0.5D0)*Q2(59,IELEM**3
     >   +1))/3 - (7*5**(0.5D0)*Q2(63,IELEM**3+1))/3 + (77*7**(0.5D0)*
     >   Q2(64,IELEM**3+1))/27 - (21**(0.5D0)*Q2(50,IELEM**3+1))/3 + (
     >   11*35**(0.5D0)*Q2(36,IELEM**3+1))/27 + 21**(0.5D0)*Q2(53,IELE
     >   M**3+1) + 35**(0.5D0)*Q2(45,IELEM**3+1) - (35**(0.5D0)*Q2(51,
     >   IELEM**3+1))/3 + 35**(0.5D0)*Q2(57,IELEM**3+1) + (11*105**(0.
     >   5D0)*Q2(28,IELEM**3+1))/27 - (105**(0.5D0)*Q2(31,IELEM**3+1))
     >   /3 + (11*105**(0.5D0)*Q2(40,IELEM**3+1))/27 - (105**(0.5D0)*Q
     >   2(46,IELEM**3+1))/3 - (105**(0.5D0)*Q2(55,IELEM**3+1))/3 - (1
     >   05**(0.5D0)*Q2(58,IELEM**3+1))/3) 
      CORNERQ(15,4) = (Q2(01,IELEM**3+1) + Q2(06,IELEM**3+1) - (5*Q
     >   2(11,IELEM**3+1))/3 - (77*Q2(16,IELEM**3+1))/27 + Q2(18,IELEM
     >   **3+1) + 3*Q2(21,IELEM**3+1) - (5*Q2(35,IELEM**3+1))/3 + 5*Q2
     >   (41,IELEM**3+1) - (77*Q2(52,IELEM**3+1))/27 + 7*Q2(61,IELEM**
     >   3+1) + (3**(0.5D0)*Q2(02,IELEM**3+1))/3 + 3**(0.5D0)*Q2(05,IE
     >   LEM**3+1) - (5**(0.5D0)*Q2(03,IELEM**3+1))/3 - (11*7**(0.5D0)
     >   *Q2(04,IELEM**3+1))/27 + 5**(0.5D0)*Q2(09,IELEM**3+1) + 3**(0
     >   .5D0)*Q2(17,IELEM**3+1) + 7**(0.5D0)*Q2(13,IELEM**3+1) - (15*
     >   *(0.5D0)*Q2(07,IELEM**3+1))/3 + 3**(0.5D0)*Q2(22,IELEM**3+1) 
     >   + (15**(0.5D0)*Q2(10,IELEM**3+1))/3 - 5**(0.5D0)*Q2(23,IELEM*
     >   *3+1) - (11*21**(0.5D0)*Q2(08,IELEM**3+1))/27 - (5*3**(0.5D0)
     >   *Q2(27,IELEM**3+1))/3 + 5**(0.5D0)*Q2(26,IELEM**3+1) - (11*7*
     >   *(0.5D0)*Q2(24,IELEM**3+1))/9 - (15**(0.5D0)*Q2(19,IELEM**3+1
     >   ))/3 - (77*3**(0.5D0)*Q2(32,IELEM**3+1))/27 + (21**(0.5D0)*Q2
     >   (14,IELEM**3+1))/3 + 7**(0.5D0)*Q2(30,IELEM**3+1) + 5**(0.5D0
     >   )*Q2(33,IELEM**3+1) + 15**(0.5D0)*Q2(25,IELEM**3+1) - (11*21*
     >   *(0.5D0)*Q2(20,IELEM**3+1))/27 - (5*3**(0.5D0)*Q2(39,IELEM**3
     >   +1))/3 + 5**(0.5D0)*Q2(38,IELEM**3+1) + (5*3**(0.5D0)*Q2(42,I
     >   ELEM**3+1))/3 - (11*35**(0.5D0)*Q2(12,IELEM**3+1))/27 - (5*5*
     >   *(0.5D0)*Q2(43,IELEM**3+1))/3 + (15**(0.5D0)*Q2(34,IELEM**3+1
     >   ))/3 + 21**(0.5D0)*Q2(29,IELEM**3+1) - (35**(0.5D0)*Q2(15,IEL
     >   EM**3+1))/3 - (55*7**(0.5D0)*Q2(44,IELEM**3+1))/27 + 15**(0.5
     >   D0)*Q2(37,IELEM**3+1) - (77*5**(0.5D0)*Q2(48,IELEM**3+1))/27 
     >   - (5*7**(0.5D0)*Q2(47,IELEM**3+1))/3 + 7**(0.5D0)*Q2(49,IELEM
     >   **3+1) - (77*3**(0.5D0)*Q2(56,IELEM**3+1))/27 + 7**(0.5D0)*Q2
     >   (54,IELEM**3+1) + (7*3**(0.5D0)*Q2(62,IELEM**3+1))/3 - (77*5*
     >   *(0.5D0)*Q2(60,IELEM**3+1))/27 - (5*7**(0.5D0)*Q2(59,IELEM**3
     >   +1))/3 - (7*5**(0.5D0)*Q2(63,IELEM**3+1))/3 - (77*7**(0.5D0)*
     >   Q2(64,IELEM**3+1))/27 + (21**(0.5D0)*Q2(50,IELEM**3+1))/3 - (
     >   11*35**(0.5D0)*Q2(36,IELEM**3+1))/27 + 21**(0.5D0)*Q2(53,IELE
     >   M**3+1) + 35**(0.5D0)*Q2(45,IELEM**3+1) - (35**(0.5D0)*Q2(51,
     >   IELEM**3+1))/3 + 35**(0.5D0)*Q2(57,IELEM**3+1) - (11*105**(0.
     >   5D0)*Q2(28,IELEM**3+1))/27 - (105**(0.5D0)*Q2(31,IELEM**3+1))
     >   /3 - (11*105**(0.5D0)*Q2(40,IELEM**3+1))/27 + (105**(0.5D0)*Q
     >   2(46,IELEM**3+1))/3 - (105**(0.5D0)*Q2(55,IELEM**3+1))/3 + (1
     >   05**(0.5D0)*Q2(58,IELEM**3+1))/3) 
      CORNERQ(16,4) = (Q2(01,IELEM**3+1) + 3*Q2(06,IELEM**3+1) + 5*
     >   Q2(11,IELEM**3+1) + 7*Q2(16,IELEM**3+1) + 3*Q2(18,IELEM**3+1)
     >    + 3*Q2(21,IELEM**3+1) + 5*Q2(35,IELEM**3+1) + 5*Q2(41,IELEM*
     >   *3+1) + 7*Q2(52,IELEM**3+1) + 7*Q2(61,IELEM**3+1) + 3**(0.5D0
     >   )*Q2(02,IELEM**3+1) + 3**(0.5D0)*Q2(05,IELEM**3+1) + 5**(0.5D
     >   0)*Q2(03,IELEM**3+1) + 7**(0.5D0)*Q2(04,IELEM**3+1) + 5**(0.5
     >   D0)*Q2(09,IELEM**3+1) + 3**(0.5D0)*Q2(17,IELEM**3+1) + 7**(0.
     >   5D0)*Q2(13,IELEM**3+1) + 15**(0.5D0)*Q2(07,IELEM**3+1) + 3*3*
     >   *(0.5D0)*Q2(22,IELEM**3+1) + 15**(0.5D0)*Q2(10,IELEM**3+1) + 
     >   3*5**(0.5D0)*Q2(23,IELEM**3+1) + 21**(0.5D0)*Q2(08,IELEM**3+1
     >   ) + 5*3**(0.5D0)*Q2(27,IELEM**3+1) + 3*5**(0.5D0)*Q2(26,IELEM
     >   **3+1) + 3*7**(0.5D0)*Q2(24,IELEM**3+1) + 15**(0.5D0)*Q2(19,I
     >   ELEM**3+1) + 7*3**(0.5D0)*Q2(32,IELEM**3+1) + 21**(0.5D0)*Q2(
     >   14,IELEM**3+1) + 3*7**(0.5D0)*Q2(30,IELEM**3+1) + 5**(0.5D0)*
     >   Q2(33,IELEM**3+1) + 15**(0.5D0)*Q2(25,IELEM**3+1) + 21**(0.5D
     >   0)*Q2(20,IELEM**3+1) + 5*3**(0.5D0)*Q2(39,IELEM**3+1) + 3*5**
     >   (0.5D0)*Q2(38,IELEM**3+1) + 5*3**(0.5D0)*Q2(42,IELEM**3+1) + 
     >   35**(0.5D0)*Q2(12,IELEM**3+1) + 5*5**(0.5D0)*Q2(43,IELEM**3+1
     >   ) + 15**(0.5D0)*Q2(34,IELEM**3+1) + 21**(0.5D0)*Q2(29,IELEM**
     >   3+1) + 35**(0.5D0)*Q2(15,IELEM**3+1) + 5*7**(0.5D0)*Q2(44,IEL
     >   EM**3+1) + 15**(0.5D0)*Q2(37,IELEM**3+1) + 7*5**(0.5D0)*Q2(48
     >   ,IELEM**3+1) + 5*7**(0.5D0)*Q2(47,IELEM**3+1) + 7**(0.5D0)*Q2
     >   (49,IELEM**3+1) + 7*3**(0.5D0)*Q2(56,IELEM**3+1) + 3*7**(0.5D
     >   0)*Q2(54,IELEM**3+1) + 7*3**(0.5D0)*Q2(62,IELEM**3+1) + 7*5**
     >   (0.5D0)*Q2(60,IELEM**3+1) + 5*7**(0.5D0)*Q2(59,IELEM**3+1) + 
     >   7*5**(0.5D0)*Q2(63,IELEM**3+1) + 7*7**(0.5D0)*Q2(64,IELEM**3+
     >   1) + 21**(0.5D0)*Q2(50,IELEM**3+1) + 35**(0.5D0)*Q2(36,IELEM*
     >   *3+1) + 21**(0.5D0)*Q2(53,IELEM**3+1) + 35**(0.5D0)*Q2(45,IEL
     >   EM**3+1) + 35**(0.5D0)*Q2(51,IELEM**3+1) + 35**(0.5D0)*Q2(57,
     >   IELEM**3+1) + 105**(0.5D0)*Q2(28,IELEM**3+1) + 105**(0.5D0)*Q
     >   2(31,IELEM**3+1) + 105**(0.5D0)*Q2(40,IELEM**3+1) + 105**(0.5
     >   D0)*Q2(46,IELEM**3+1) + 105**(0.5D0)*Q2(55,IELEM**3+1) + 105*
     >   *(0.5D0)*Q2(58,IELEM**3+1)) 
*** ---------------------------------------------------------------- ***
*** ---------------------------------------------------------------- ***
*** ---------------------------------------------------------------- ***
      ENDIF
*
      DO JEL=1,IELEM
      DO IEL=1,IELEM
         IF(IND.LE.6)THEN
            XNK(IEL,JEL)=CORNERQ(IEL + (JEL-1)*IELEM,IELEM)
         ELSEIF(IND.GE.7)THEN
            XNK(IEL,JEL)=CORNERQ(IEL + (JEL-1)*IELEM,1)
         ENDIF
      ENDDO
      ENDDO
*
      IF(JL.LT.ISPLH)THEN
         DO JEL=1,IELEM
         DO IEL=1,IELEM
            IF(IFROMJ.EQ.2)THEN
               XNJ(IEL,JEL,K2)=CORNERQ(IELEM**2 - IELEM +IEL,JEL)
            ELSEIF(IFROMJ.EQ.4)THEN
               XNJ(IEL,JEL,K2)=CORNERQ(IEL,JEL)
            ENDIF
         ENDDO
         ENDDO
      ELSEIF((JL.EQ.ISPLH).AND.(IHEXJ.LE.NHEX))THEN
         DO JEL=1,IELEM
         DO INDEXE=1,IELEM
            IEL=INDEXE
            IELJ=INDEXE
            I3=I2
            IF((J.EQ.1).AND.(ILOZJ.EQ.3))THEN
               IELJ=IELEM - (INDEXE-1)
               I3=ISPLH+1 -I2
            ENDIF
            IF(IFROMJ.EQ.2)THEN
               BFLUX(ISIDEJ,INDEXJ,I3,IELJ,JEL,K2) = CORNERQ
     1            (IELEM**2-IELEM+IEL,JEL)
            ELSEIF(IFROMJ.EQ.4)THEN
               BFLUX(ISIDEJ,INDEXJ,I3,IELJ,JEL,K2) = CORNERQ(IEL,JEL)
            ENDIF
         ENDDO
         ENDDO
      ENDIF
*
      IF(IL.LT.ISPLH)THEN
         DO JEL=1,IELEM
         DO IEL=1,IELEM
            IF(IFROMI.EQ.1)THEN
               XNI(IEL,JEL,J2,K2)=CORNERQ(IEL*IELEM,JEL)
            ELSEIF(IFROMI.EQ.3)THEN
               XNI(IEL,JEL,J2,K2)=CORNERQ((IELEM*(IEL-1))+1,JEL)
            ENDIF
         ENDDO
         ENDDO
      ELSEIF((IL.EQ.ISPLH).AND.(IHEXI.LE.NHEX))THEN
         DO JEL=1,IELEM
         DO INDEXE=1,IELEM
            IEL=INDEXE
            IELI=INDEXE
            J3=J2
            IF((J.EQ.3).AND.(ILOZI.EQ.1))THEN
               IELI=IELEM - (INDEXE-1)
               J3=ISPLH+1-J2
            ENDIF
            IF(IFROMI.EQ.1)THEN
               BFLUX(ISIDEI,INDEXI,J3,IELI,JEL,K2)=CORNERQ
     1            (IEL*IELEM,JEL)
            ELSEIF(IFROMI.EQ.3)THEN
               BFLUX(ISIDEI,INDEXI,J3,IELI,JEL,K2)=CORNERQ
     1            ((IELEM*(IEL-1))+1,JEL)
            ENDIF
         ENDDO
         ENDDO
      ENDIF
*
      DO 131 K=1,NSCT
      DO 130 IEL=1,IELEM**3
      FLUX(IEL,K,I2,J2,J,I,K2) = FLUX(IEL,K,I2,J2,J,I,K2) +
     1   W(M)*REAL(Q2(IEL,IELEM**3+1))*PL(K,M)
  130 CONTINUE
  131 CONTINUE
*
  150 CONTINUE

      DO 156 IEL=1,IELEM
      DO 155 JEL=1,IELEM
      XNEK(IEL,JEL,I2,J2,J,I,M)=REAL(XNK(IEL,JEL))
  155 CONTINUE
  156 CONTINUE

  160 CONTINUE
  170 CONTINUE
  180 CONTINUE
  190 CONTINUE
  210 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(BFLUX)
      RETURN
      END
