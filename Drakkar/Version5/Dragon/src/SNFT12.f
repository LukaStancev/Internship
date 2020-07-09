*DECK SNFT12
      SUBROUTINE SNFT12(LX,LY,IELEM,NMAT,NPQ,NSCT,MAT,VOL,TOTAL,NCODE,
     1 ZCODE,QEXT,LFIXUP,DU,DE,W,MRM,MRMY,DB,DA,PL,FLUX,XNEI,XNEJ)
*
*-----------------------------------------------------------------------
*
*Purpose:
* perform one inner iteration for solving SN equations in 2D Cartesian
* geometry for the HODD method. Albedo boundary conditions.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* LX      number of meshes along X axis.
* LY      number of meshes along Y axis.
* IELEM   measure of order of the spatial approximation polynomial:
*         =1: constant - classical diamond scheme - default for HODD;
*         =2: linear;
*         =3: parabolic.
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
*
*Parameters: output
* FLUX    Legendre components of the flux.
* XNEI    X-directed SN boundary fluxes.
* XNEJ    Y-directed SN boundary fluxes.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER LX,LY,IELEM,NMAT,NPQ,NSCT,MAT(LX,LY),NCODE(4),MRM(NPQ),
     1 MRMY(NPQ)
      REAL VOL(LX,LY),TOTAL(0:NMAT),ZCODE(4),QEXT(IELEM**2,NSCT,LX,LY),
     1 DU(NPQ),DE(NPQ),W(NPQ),DB(LX,NPQ),DA(LX,LY,NPQ),PL(NSCT,NPQ),
     2 FLUX(IELEM**2,NSCT,LX,LY),XNEI(IELEM,LY,NPQ),XNEJ(IELEM,LX,NPQ)
      LOGICAL LFIXUP
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION Q(9),XNJ(3),Q2(9,10),CONST0,CONST1,CONST2
      PARAMETER(RLOG=1.0E-8,PI=3.141592654)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: XNI
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(XNI(IELEM,LY))
*----
*  DEFINITION OF CONSTANTS.
*----
      CONST0=2.0D0*DSQRT(3.0D0)
      CONST1=2.0D0*DSQRT(5.0D0)
      CONST2=2.0D0*DSQRT(15.0D0)
*----
*  PARAMETER VALIDATION.
*----
      IF(IELEM.GT.3) CALL XABORT('SNFT12: INVALID IELEM (DIAM) VALUE. '
     1 //'CHECK INPUT DATA FILE.')
*----
*  MAIN LOOP OVER SN ANGLES.
*----
      CALL XDRSET(FLUX,IELEM*IELEM*LX*LY*NSCT,0.0)
      DO 170 M=1,NPQ
      WEIGHT=W(M)
      VU=DU(M)
      VE=DE(M)
      IF(NCODE(1).NE.4) THEN
         M1=MRM(M)
         DO 30 IEL=1,IELEM
         IF(WEIGHT.EQ.0.0) THEN
            DO 10 J=1,LY
            XNEI(IEL,J,M)=XNEI(IEL,J,M1)
   10       CONTINUE
         ELSE IF(VU.GT.0.0) THEN
            DO 20 J=1,LY
            E1=XNEI(IEL,J,M)
            XNEI(IEL,J,M)=XNEI(IEL,J,M1)
            XNEI(IEL,J,M1)=E1
   20       CONTINUE
         ENDIF
   30    CONTINUE
      ENDIF
      IF(NCODE(3).NE.4) THEN
         M1=MRMY(M)
         DO 50 IEL=1,IELEM
         IF(VE.GT.0) THEN
            DO 40 I=1,LX
            E1=XNEJ(IEL,I,M)
            XNEJ(IEL,I,M)=XNEJ(IEL,I,M1)
            XNEJ(IEL,I,M1)=E1
   40       CONTINUE
         ENDIF
   50    CONTINUE
      ENDIF
      IF(WEIGHT.EQ.0.0) GO TO 170
      IF(VE.GT.0.0) GOTO 70
      IF(VU.GT.0.0) GOTO 60
      IND=3
      GOTO 90
   60 IND=4
      GOTO 90
   70 IF(VU.GT.0.0) GOTO 80
      IND=2
      GOTO 90
   80 IND=1
*----
*  LOOP OVER X- AND Y-DIRECTED AXES.
*----
   90 DO 155 IZ=1,LX
      I=IZ
      IF((IND.EQ.2).OR.(IND.EQ.3)) I=LX+1-IZ
      DO 100 IEL=1,IELEM
      IF((IND.EQ.1).OR.(IND.EQ.2)) THEN
         XNJ(IEL)=XNEJ(IEL,I,M)*ZCODE(3)
      ELSE
         XNJ(IEL)=XNEJ(IEL,I,M)*ZCODE(4)
      ENDIF
  100 CONTINUE
      DO 140 JJ=1,LY
      J=JJ
      IF((IND.EQ.3).OR.(IND.EQ.4)) J=LY+1-JJ
      DO 105 IEL=1,IELEM
      IF(IZ.EQ.1) THEN
         IF((IND.EQ.1).OR.(IND.EQ.4)) THEN
            XNI(IEL,J)=XNEI(IEL,J,M)*ZCODE(1)
         ELSE
            XNI(IEL,J)=XNEI(IEL,J,M)*ZCODE(2)
         ENDIF
      ENDIF
  105 CONTINUE
      IF(MAT(I,J).EQ.0) GO TO 140
      DO 115 IEL=1,IELEM**2
      Q(IEL)=0.0
      DO 110 K=1,NSCT
      Q(IEL)=Q(IEL)+QEXT(IEL,K,I,J)*PL(K,M)/(4.0*PI)
  110 CONTINUE
  115 CONTINUE
      VT=VOL(I,J)*TOTAL(MAT(I,J))
      CALL XDDSET(Q2,90,0.0D0)
      IF(IELEM.EQ.1) THEN
         Q2(1,1)=2.0D0*ABS(DA(I,J,M))+2.0D0*ABS(DB(I,M))+VT
         Q2(1,2)=2.0D0*ABS(DA(I,J,M))*XNI(1,J)+2.0D0*ABS(DB(I,M))
     1           *XNJ(1)+VOL(I,J)*Q(1)
      ELSE IF(IELEM.EQ.2) THEN
         Q2(1,1)=VT
         Q2(2,1)=CONST0*DA(I,J,M)
         Q2(2,2)=-VT-6.0D0*ABS(DA(I,J,M))
         Q2(3,1)=CONST0*DB(I,M)
         Q2(3,3)=-VT-6.0D0*ABS(DB(I,M))
         Q2(4,2)=-CONST0*DB(I,M)
         Q2(4,3)=-CONST0*DA(I,J,M)
         Q2(4,4)=VT+6.0D0*ABS(DA(I,J,M))+6.0D0*ABS(DB(I,M))
*        ------
         Q2(1,5)=VOL(I,J)*Q(1)
         Q2(2,5)=-VOL(I,J)*Q(2)+CONST0*DA(I,J,M)*XNI(1,J)
         Q2(3,5)=-VOL(I,J)*Q(3)+CONST0*DB(I,M)*XNJ(1)
         Q2(4,5)=VOL(I,J)*Q(4)-CONST0*DA(I,J,M)*XNI(2,J)-CONST0*
     1           DB(I,M)*XNJ(2)
      ELSE IF(IELEM.EQ.3) THEN
         Q2(1,1)=VT+2.0D0*ABS(DA(I,J,M))+2.0D0*ABS(DB(I,M))
         Q2(2,2)=-VT-2.0D0*ABS(DB(I,M))
         Q2(3,1)=CONST1*ABS(DA(I,J,M))
         Q2(3,2)=-CONST2*DA(I,J,M)
         Q2(3,3)=VT+1.0D1*ABS(DA(I,J,M))+2.0D0*ABS(DB(I,M))
         Q2(4,4)=-VT-2.0D0*ABS(DA(I,J,M))
         Q2(5,5)=VT
         Q2(6,4)=-CONST1*ABS(DA(I,J,M))
         Q2(6,5)=CONST2*DA(I,J,M)
         Q2(6,6)=-VT-1.0D1*ABS(DA(I,J,M))
         Q2(7,1)=CONST1*ABS(DB(I,M))
         Q2(7,4)=-CONST2*DB(I,M)
         Q2(7,7)=VT+2.0D0*ABS(DA(I,J,M))+1.0D1*ABS(DB(I,M))
         Q2(8,2)=-CONST1*ABS(DB(I,M))
         Q2(8,5)=CONST2*DB(I,M)
         Q2(8,8)=-VT-1.0D1*ABS(DB(I,M))
         Q2(9,3)=CONST1*ABS(DB(I,M))
         Q2(9,6)=-CONST2*DB(I,M)
         Q2(9,7)=CONST1*ABS(DA(I,J,M))
         Q2(9,8)=-CONST2*DA(I,J,M)
         Q2(9,9)=VT+1.0D1*ABS(DA(I,J,M))+1.0D1*ABS(DB(I,M))
*        ------
         Q2(1,10)=VOL(I,J)*Q(1)+2.0D0*ABS(DA(I,J,M))*XNI(1,J)+2.0D0*
     1            ABS(DB(I,M))*XNJ(1)
         Q2(2,10)=-VOL(I,J)*Q(2)-2.0D0*ABS(DB(I,M))*XNJ(2)
         Q2(3,10)=VOL(I,J)*Q(3)+CONST1*ABS(DA(I,J,M))*XNI(1,J)+2.0D0*
     1            ABS(DB(I,M))*XNJ(3)
         Q2(4,10)=-VOL(I,J)*Q(4)-2.0D0*ABS(DA(I,J,M))*XNI(2,J)
         Q2(5,10)=VOL(I,J)*Q(5)
         Q2(6,10)=-VOL(I,J)*Q(6)-CONST1*ABS(DA(I,J,M))*XNI(2,J)
         Q2(7,10)=VOL(I,J)*Q(7)+2.0D0*ABS(DA(I,J,M))*XNI(3,J)+CONST1*
     1            ABS(DB(I,M))*XNJ(1)
         Q2(8,10)=-VOL(I,J)*Q(8)-CONST1*ABS(DB(I,M))*XNJ(2)
         Q2(9,10)=VOL(I,J)*Q(9)+CONST1*ABS(DA(I,J,M))*XNI(3,J)+CONST1*
     1            ABS(DB(I,M))*XNJ(3)
      ENDIF
      DO 125 IEL=1,IELEM**2
      DO 120 JEL=IEL+1,IELEM**2
      Q2(IEL,JEL)=Q2(JEL,IEL)
  120 CONTINUE
  125 CONTINUE
      CALL ALSBD(IELEM**2,1,Q2,IER,9)
      IF(IER.NE.0) CALL XABORT('SNFT12: SINGULAR MATRIX.')
      IF(IELEM.EQ.1) THEN
         IF(LFIXUP.AND.(Q2(1,2).LE.RLOG)) Q2(1,2)=0.0
         XNI(1,J)=2.0D0*Q2(1,2)-XNI(1,J)
         XNJ(1)=2.0D0*Q2(1,2)-XNJ(1)
         IF(LFIXUP.AND.(XNI(1,J).LE.RLOG)) XNI(1,J)=0.0
         IF(LFIXUP.AND.(XNJ(1).LE.RLOG)) XNJ(1)=0.0
      ELSE IF(IELEM.EQ.2) THEN
         XNI(1,J)=XNI(1,J)+SIGN(1.0,DU(M))*CONST0*Q2(2,5)
         XNI(2,J)=XNI(2,J)+SIGN(1.0,DU(M))*CONST0*Q2(4,5)
         XNJ(1)=XNJ(1)+SIGN(1.0,DE(M))*CONST0*Q2(3,5)
         XNJ(2)=XNJ(2)+SIGN(1.0,DE(M))*CONST0*Q2(4,5)
      ELSE IF(IELEM.EQ.3) THEN
         XNI(1,J)=2.0D0*Q2(1,10)+CONST1*Q2(3,10)-XNI(1,J)
         XNI(2,J)=2.0D0*Q2(4,10)+CONST1*Q2(6,10)-XNI(2,J)
         XNI(3,J)=2.0D0*Q2(7,10)+CONST1*Q2(9,10)-XNI(3,J)
         XNJ(1)=2.0D0*Q2(1,10)+CONST1*Q2(7,10)-XNJ(1)
         XNJ(2)=2.0D0*Q2(2,10)+CONST1*Q2(8,10)-XNJ(2)
         XNJ(3)=2.0D0*Q2(3,10)+CONST1*Q2(9,10)-XNJ(3)
      ENDIF
      DO 135 K=1,NSCT
      DO 130 IEL=1,IELEM**2
      FLUX(IEL,K,I,J)=FLUX(IEL,K,I,J)+2.0*W(M)*REAL(Q2(IEL,IELEM**2+1))*
     1 PL(K,M)
  130 CONTINUE
  135 CONTINUE
  140 CONTINUE
      DO 150 IEL=1,IELEM
      XNEJ(IEL,I,M)=REAL(XNJ(IEL))
  150 CONTINUE
  155 CONTINUE
      DO 165 J=1,LY
      DO 160 IEL=1,IELEM
      XNEI(IEL,J,M)=REAL(XNI(IEL,J))
  160 ENDDO
  165 ENDDO
  170 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(XNI)
      RETURN
      END
