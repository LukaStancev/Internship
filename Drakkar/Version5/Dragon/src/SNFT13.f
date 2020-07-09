*DECK SNFT13
      SUBROUTINE SNFT13(LX,LY,LZ,IELEM,NMAT,NPQ,NSCT,MAT,VOL,TOTAL,
     1 NCODE,ZCODE,QEXT,LFIXUP,DU,DE,DZ,W,MRMX,MRMY,MRMZ,DC,DB,DA,PL,
     2 FLUX,XNEI,XNEJ,XNEK)
*
*-----------------------------------------------------------------------
*
*Purpose:
* perform one inner iteration for solving SN equations in 3D Cartesian 
* geometry for the HODD method. Albedo boundary conditions.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): N. Martin
*
*Parameters: input
* LX      number of meshes along X axis.
* LY      number of meshes along Y axis.
* LZ      number of meshes along Z axis.
* IELEM   degree of the spatial approximation (=1: diamond scheme).
* NMAT    number of material mixtures.
* NPQ     number of SN directions in height octants.
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
* DZ      third direction cosines ($\xi$).
* W       weights.
* MRMX    quadrature index.
* MRMY    quadrature index.
* MRMZ    quadrature index.
* DC      diamond-scheme parameter.
* DB      diamond-scheme parameter.
* DA      diamond-scheme parameter.
* PL      discrete values of the spherical harmonics corresponding
*         to the 3D SN quadrature.
* XNEI    X-directed SN boundary fluxes.
* XNEJ    Y-directed SN boundary fluxes.
* XNEK    Z-directed SN boundary fluxes.
*
*Parameters: output
* FLUX    Legendre components of the flux.
* XNEI    X-directed SN boundary fluxes.
* XNEJ    Y-directed SN boundary fluxes.
* XNEK    Z-directed SN boundary fluxes.
*
*-----------------------------------------------------------------------
*
      INTEGER LX,LY,LZ,IELEM,NMAT,NPQ,NSCT,MAT(LX,LY,LZ),
     1 NCODE(6),MRMX(NPQ),MRMY(NPQ),MRMZ(NPQ)
      REAL VOL(LX,LY,LZ),TOTAL(0:NMAT),ZCODE(6),
     1 QEXT(IELEM**3,NSCT,LX,LY,LZ),
     2 DU(NPQ),DE(NPQ),DZ(NPQ),W(NPQ),DC(LX,LY,NPQ),DB(LX,LZ,NPQ),
     3 DA(LY,LZ,NPQ), PL(NSCT,NPQ),FLUX(IELEM**3,NSCT,LX,LY,LZ),
     4 XNEI(IELEM,IELEM,LY,LZ,NPQ), XNEJ(IELEM,IELEM,LX,LZ,NPQ),
     5 XNEK(IELEM,IELEM,LX,LY,NPQ) 
      LOGICAL LFIXUP
*----
*  LOCAL VARIABLES
*----
      PARAMETER(RLOG=1.0E-8,PI=3.141592654)
      DOUBLE PRECISION CONST0,CONST1,CONST2,Q(IELEM**3),
     >   Q2(IELEM**3,IELEM**3 +1),XNK(IELEM,IELEM)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: XNI
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: XNJ
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(XNI(IELEM,IELEM,LY,LZ),XNJ(IELEM,IELEM,LZ))
*
      CONST0=2.0D0*DSQRT(3.0D0)
      CONST1=2.0D0*DSQRT(5.0D0)
      CONST2=2.0D0*DSQRT(15.0D0) 
*----
*  PARAMETER VALIDATION.
*----
      IF((IELEM.LT.1).OR.(IELEM.GT.3))
     1   CALL XABORT('SNFTH1: INVALID IELEM (DIAM) VALUE. '
     2   //'CHECK INPUT DATA FILE.')   
*----
*  MAIN LOOP OVER SN ANGLES.
*----
      CALL XDRSET(FLUX,IELEM*IELEM*IELEM*LX*LY*LZ*NSCT,0.0)
*
      DO 170 M=1,NPQ
      WEIGHT=W(M)
      VU=DU(M)
      VE=DE(M)
      VZ=DZ(M)
*
      IF(NCODE(1).NE.4) THEN
         M1=MRMX(M)
         DO 31 IEL=1,IELEM
         DO 30 JEL=1,IELEM           
         IF(VU.GT.0.0) THEN
             DO 26 J=1,LY
             DO 25 K=1,LZ
             E1=XNEI(IEL,JEL,J,K,M)
             XNEI(IEL,JEL,J,K,M)=XNEI(IEL,JEL,J,K,M1)
             XNEI(IEL,JEL,J,K,M1)=E1
   25        CONTINUE
   26        CONTINUE
         ENDIF
   30    CONTINUE
   31    CONTINUE
      ENDIF
      IF(NCODE(3).NE.4) THEN
         M1=MRMY(M)
         DO 51 IEL=1,IELEM
         DO 50 JEL=1,IELEM
         IF(VE.GT.0.0) THEN
            DO 41 I=1,LX
            DO 40 K=1,LZ
            E1=XNEJ(IEL,JEL,I,K,M)
            XNEJ(IEL,JEL,I,K,M)=XNEJ(IEL,JEL,I,K,M1)
            XNEJ(IEL,JEL,I,K,M1)=E1
   40       CONTINUE
   41       CONTINUE
         ENDIF
   50    CONTINUE
   51    CONTINUE
      ENDIF
      IF (NCODE(5).NE.4) THEN
         M1=MRMZ(M)
         DO 61 IEL=1,IELEM
         DO 60 JEL=1,IELEM
         IF (VZ.GT.0.0) THEN
           DO  71 I=1,LX
           DO  70 J=1,LY
           E1=XNEK(IEL,JEL,I,J,M)
           XNEK(IEL,JEL,I,J,M)=XNEK(IEL,JEL,I,J,M1)
           XNEK(IEL,JEL,I,J,M1)=E1
  70       CONTINUE
  71       CONTINUE
         ENDIF
  60  CONTINUE
  61  CONTINUE
      ENDIF
*---- UNFOLD OCTANTS
      IND=0
      IF((VU.GE.0.0).AND.(VE.GE.0.0).AND.(VZ.GE.0.0)) THEN
      IND=1
      ELSE IF((VU.LE.0.0).AND.(VE.GE.0.0).AND.(VZ.GE.0.0)) THEN
      IND=2
      ELSE IF((VU.LE.0.0).AND.(VE.LE.0.0).AND.(VZ.GE.0.0)) THEN
      IND=3
      ELSE IF((VU.GE.0.0).AND.(VE.LE.0.0).AND.(VZ.GE.0.0)) THEN
      IND=4
      ELSE IF((VU.GE.0.0).AND.(VE.GE.0.0).AND.(VZ.LE.0.0)) THEN
      IND=5
      ELSE IF((VU.LE.0.0).AND.(VE.GE.0.0).AND.(VZ.LE.0.0)) THEN   
      IND=6
      ELSE IF((VU.LE.0.0).AND.(VE.LE.0.0).AND.(VZ.LE.0.0)) THEN
      IND=7
      ELSE IF((VU.GE.0.0).AND.(VE.LE.0.0).AND.(VZ.LE.0.0)) THEN
      IND=8
      ENDIF   
*----
*  LOOP OVER X-,Y- AND Z-DIRECTED AXES.
*----     
      DO 163 IX=1,LX
      I=IX
      IF((IND.EQ.2).OR.(IND.EQ.3).OR.(IND.EQ.6).OR.(IND.EQ.7)) I=LX+1-IX
      DO 152 IY=1,LY
      J=IY
      IF((IND.EQ.3).OR.(IND.EQ.4).OR.(IND.EQ.7).OR.(IND.EQ.8)) J=LY+1-IY
      DO 101 IEL=1,IELEM
      DO 100 JEL=1,IELEM
         IF((IND.EQ.1).OR.(IND.EQ.2).OR.(IND.EQ.3).OR.(IND.EQ.4)) THEN  
            XNK(IEL,JEL)=XNEK(IEL,JEL,I,J,M)*ZCODE(5)
         ELSE
            XNK(IEL,JEL)=XNEK(IEL,JEL,I,J,M)*ZCODE(6)
         ENDIF
  100 CONTINUE    
  101 CONTINUE    
      DO 140 IZ=1,LZ     
      K=IZ
      IF((IND.EQ.5).OR.(IND.EQ.6).OR.(IND.EQ.7).OR.(IND.EQ.8)) K=LZ+1-IZ
      DO 111 IEL=1,IELEM
      DO 110 JEL=1,IELEM
         IF(IY.EQ.1) THEN
           IF((IND.EQ.1).OR.(IND.EQ.2).OR.(IND.EQ.5).OR.(IND.EQ.6)) THEN
           XNJ(IEL,JEL,K)=XNEJ(IEL,JEL,I,K,M)*ZCODE(3)
           ELSE
           XNJ(IEL,JEL,K)=XNEJ(IEL,JEL,I,K,M)*ZCODE(4)
           ENDIF
         ENDIF
  110 CONTINUE
  111 CONTINUE
      DO 121 IEL=1,IELEM
      DO 120 JEL=1,IELEM
         IF(IX.EQ.1) THEN
           IF((IND.EQ.1).OR.(IND.EQ.4).OR.(IND.EQ.5).OR.(IND.EQ.8)) THEN
           XNI(IEL,JEL,J,K)=XNEI(IEL,JEL,J,K,M)*ZCODE(1)
           ELSE
           XNI(IEL,JEL,J,K)=XNEI(IEL,JEL,J,K,M)*ZCODE(2)
           ENDIF
         ENDIF    
  120 CONTINUE 
  121 CONTINUE 
      IF(MAT(I,J,K).EQ.0) GO TO 140
      DO 131 IEL=1,IELEM**3
      Q(IEL)=0.0
      DO 130 L=1,NSCT    
      Q(IEL)=Q(IEL)+QEXT(IEL,L,I,J,K)*PL(L,M)/(4.0*PI)
  130 CONTINUE
  131 CONTINUE
      VT=VOL(I,J,K)*TOTAL(MAT(I,J,K))
      CALL XDDSET(Q2,(IELEM**3)*((IELEM**3)+1),0.0D0)
*---------------------------------------
      IF(IELEM.EQ.1) THEN  
      Q2(1,1)=2.0D0*ABS(DA(J,K,M))+2.0D0*ABS(DB(I,K,M))+
     1        2.0D0*ABS(DC(I,J,M))+VT 
      Q2(1,2)=2.0D0*ABS(DA(J,K,M))*XNI(1,1,J,K)+
     1        2.0D0*ABS(DB(I,K,M))*XNJ(1,1,K)+
     2        2.0D0*ABS(DC(I,J,M))*XNK(1,1)+VOL(I,J,K)*Q(1)    
      ELSE IF(IELEM.EQ.2) THEN
      Q2(1,1)=VT
      Q2(1,2)=CONST0*DA(J,K,M) 
      Q2(1,3)=CONST0*DB(I,K,M)
      Q2(1,4)=CONST0*DC(I,J,M)
      Q2(2,2)=-VT-6.0D0*ABS(DA(J,K,M))    
      Q2(2,5)=-CONST0*DB(I,K,M)
      Q2(2,6)=-CONST0*DC(I,J,M)
      Q2(3,3)=-VT-6.0D0*ABS(DB(I,K,M))   
      Q2(3,5)=-CONST0*DA(J,K,M)    
      Q2(3,7)=-CONST0*DC(I,J,M)
      Q2(4,4)=-VT-6.0D0*ABS(DC(I,J,M))    
      Q2(4,6)=-CONST0*DA(J,K,M)
      Q2(4,7)=-CONST0*DB(I,K,M)
      Q2(5,5)=VT+6.0D0*ABS(DA(J,K,M))+6.0D0*ABS(DB(I,K,M))    
      Q2(5,8)=CONST0*DC(I,J,M) 
      Q2(6,6)=VT+6.0D0*ABS(DA(J,K,M))+6.0D0*ABS(DC(I,J,M))    
      Q2(6,8)=CONST0*DB(I,K,M) 
      Q2(7,7)=VT+6.0D0*ABS(DB(I,K,M))+6.0D0*ABS(DC(I,J,M))
      Q2(7,8)=CONST0*DA(J,K,M) 
      Q2(8,8)=-VT-6.0D0*ABS(DA(J,K,M))-6.0D0*ABS(DB(I,K,M))
     1        -6.0D0*ABS(DC(I,J,M))
*---------------
      Q2(1,9)=VOL(I,J,K)*Q(1)
      Q2(2,9)=-VOL(I,J,K)*Q(2)+CONST0*DA(J,K,M)*XNI(1,1,J,K)
      Q2(3,9)=-VOL(I,J,K)*Q(3)+CONST0*DB(I,K,M)*XNJ(1,1,K)
      Q2(4,9)=-VOL(I,J,K)*Q(4)+CONST0*DC(I,J,M)*XNK(1,1)
      Q2(5,9)=VOL(I,J,K)*Q(5)-CONST0*DA(J,K,M)*XNI(2,1,J,K)
     1       -CONST0*DB(I,K,M)*XNJ(2,1,K)
      Q2(6,9)=VOL(I,J,K)*Q(6)-CONST0*DA(J,K,M)*XNI(1,2,J,K)
     1       -CONST0*DC(I,J,M)*XNK(2,1)
      Q2(7,9)=VOL(I,J,K)*Q(7)-CONST0*DB(I,K,M)*XNJ(1,2,K)
     1       -CONST0*DC(I,J,M)*XNK(1,2)
      Q2(8,9)=-VOL(I,J,K)*Q(8)+CONST0*DA(J,K,M)*XNI(2,2,J,K)
     1       +CONST0*DB(I,K,M)*XNJ(2,2,K)
     2       +CONST0*DC(I,J,M)*XNK(2,2)  

      ELSE IF (IELEM.EQ.3) THEN
* UPPER DIAGONAL TERMS
      Q2(1,1)=-VT-2.0D0*ABS(DA(J,K,M))-2.0D0*ABS(DB(I,K,M))-
     1        2.0D0*ABS(DC(I,J,M))  
      Q2(1,9)=-CONST1*ABS(DA(J,K,M))
      Q2(1,10)=-CONST1*ABS(DB(I,K,M))
      Q2(1,11)=-CONST1*ABS(DC(I,J,M))
      Q2(2,2)=VT+2.0D0*ABS(DB(I,K,M))+2.0D0*ABS(DC(I,J,M))
      Q2(2,9)=CONST2*DA(J,K,M)
      Q2(2,16)=CONST1*ABS(DC(I,J,M))
      Q2(2,17)=CONST1*ABS(DB(I,K,M))
      Q2(3,3)=VT+2.0D0*ABS(DA(J,K,M))+2.0D0*ABS(DC(I,J,M))
      Q2(3,10)=CONST2*DB(I,K,M)
      Q2(3,12)=CONST1*ABS(DA(J,K,M))
      Q2(3,15)=CONST1*ABS(DC(I,J,M))
      Q2(4,4)=VT+2.0D0*ABS(DA(J,K,M))+2.0D0*ABS(DB(I,K,M))
      Q2(4,11)=CONST2*DC(I,J,M)
      Q2(4,13)=CONST1*ABS(DA(J,K,M))
      Q2(4,19)=CONST1*ABS(DB(I,K,M))
      Q2(5,5)=-VT-2.0D0*ABS(DC(I,J,M))
      Q2(5,12)=-CONST2*DA(J,K,M)
      Q2(5,17)=-CONST2*DB(I,K,M)
      Q2(5,20)=-CONST1*ABS(DC(I,J,M))
      Q2(6,6)=-VT-2.0D0*ABS(DA(J,K,M))
      Q2(6,14)=-CONST1*ABS(DA(J,K,M))
      Q2(6,15)=-CONST2*DC(I,J,M)
      Q2(6,19)=-CONST2*DB(I,K,M)
      Q2(7,7)=-VT-2.0D0*ABS(DB(I,K,M))
      Q2(7,13)=-CONST2*DA(J,K,M)
      Q2(7,16)=-CONST2*DC(I,J,M)
      Q2(7,18)=-CONST1*ABS(DB(I,K,M))
      Q2(8,8)=VT
      Q2(8,14)=CONST2*DA(J,K,M)
      Q2(8,18)=CONST2*DB(I,K,M)
      Q2(8,20)=CONST2*DC(I,J,M)
      Q2(9,9)=-VT-1.0D1*ABS(DA(J,K,M))-2.0D0*ABS(DB(I,K,M))-
     1         2.0D0*ABS(DC(I,J,M))
      Q2(9,21)=-CONST1*ABS(DB(I,K,M))
      Q2(9,23)=-CONST1*ABS(DC(I,J,M))
      Q2(10,10)=-VT-2.0D0*ABS(DA(J,K,M))-1.0D1*ABS(DB(I,K,M))-
     1        2.0D0*ABS(DC(I,J,M))
      Q2(10,21)=-CONST1*ABS(DA(J,K,M))
      Q2(10,22)=-CONST1*ABS(DC(I,J,M))
      Q2(11,11)=-VT-2.0D0*ABS(DA(J,K,M))-2.0D0*ABS(DB(I,K,M))-
     1        1.0D1*ABS(DC(I,J,M))
      Q2(11,22)=-CONST1*ABS(DB(I,K,M))
      Q2(11,23)=-CONST1*ABS(DA(J,K,M))
      Q2(12,12)=VT+1.0D1*ABS(DA(J,K,M))+2.0D0*ABS(DC(I,J,M))
      Q2(12,21)=CONST2*DB(I,K,M)
      Q2(12,25)=CONST1*ABS(DC(I,J,M))
      Q2(13,13)=VT+1.0D1*ABS(DA(J,K,M))+2.0D0*ABS(DB(I,K,M))
      Q2(13,23)=CONST2*DC(I,J,M)
      Q2(13,26)=CONST1*ABS(DB(I,K,M))
      Q2(14,14)=-VT-1.0D1*ABS(DA(J,K,M))
      Q2(14,25)=-CONST2*DC(I,J,M)
      Q2(14,26)=-CONST2*DB(I,K,M)
      Q2(15,15)=VT+2.0D0*ABS(DA(J,K,M))+1.0D1*ABS(DC(I,J,M))
      Q2(15,22)=CONST2*DB(I,K,M)
      Q2(15,25)=CONST1*ABS(DA(J,K,M))
      Q2(16,16)=VT+2.0D0*ABS(DB(I,K,M))+1.0D1*ABS(DC(I,J,M))
      Q2(16,23)=CONST2*DA(J,K,M)
      Q2(16,24)=CONST1*ABS(DB(I,K,M))
      Q2(17,17)=VT+1.0D1*ABS(DB(I,K,M))+2.0D0*ABS(DC(I,J,M))
      Q2(17,21)=CONST2*DA(J,K,M)
      Q2(17,24)=CONST1*ABS(DC(I,J,M))
      Q2(18,18)=-VT-1.0D1*ABS(DB(I,K,M))
      Q2(18,24)=-CONST2*DC(I,J,M)
      Q2(18,26)=-CONST2*DA(J,K,M)
      Q2(19,19)=VT+2.0D0*ABS(DA(J,K,M))+1.0D1*ABS(DB(I,K,M))
      Q2(19,22)=CONST2*DC(I,J,M)
      Q2(19,26)=CONST1*ABS(DA(J,K,M))
      Q2(20,20)=-VT-1.0D1*ABS(DC(I,J,M))
      Q2(20,24)=-CONST2*DB(I,K,M)
      Q2(20,25)=-CONST2*DA(J,K,M)
      Q2(21,21)=-VT-1.0D1*ABS(DA(J,K,M))-1.0D1*ABS(DB(I,K,M))-
     1        2.0D0*ABS(DC(I,J,M))
      Q2(21,27)=-CONST1*ABS(DC(I,J,M))
      Q2(22,22)=-VT-1.0D1*ABS(DB(I,K,M))-1.0D1*ABS(DC(I,J,M))-
     1        2.0D0*ABS(DA(J,K,M))
      Q2(22,27)=-CONST1*ABS(DA(J,K,M))
      Q2(23,23)=-VT-1.0D1*ABS(DA(J,K,M))-2.0D0*ABS(DB(I,K,M))-
     1        1.0D1*ABS(DC(I,J,M))
      Q2(23,27)=-CONST1*ABS(DB(I,K,M))
      Q2(24,24)=VT+1.0D1*ABS(DB(I,K,M))+1.0D1*ABS(DC(I,J,M))
      Q2(24,27)=CONST2*DA(J,K,M)
      Q2(25,25)=VT+1.0D1*ABS(DA(J,K,M))+1.0D1*ABS(DC(I,J,M))
      Q2(25,27)=CONST2*DB(I,K,M)
      Q2(26,26)=VT+1.0D1*ABS(DA(J,K,M))+1.0D1*ABS(DB(I,K,M))
      Q2(26,27)=CONST2*DC(I,J,M)
      Q2(27,27)=-VT-1.0D1*ABS(DA(J,K,M))-1.0D1*ABS(DB(I,K,M))-
     1        1.0D1*ABS(DC(I,J,M))    
*---------------------
      Q2(1,28)=-VOL(I,J,K)*Q(1)-2.0D0*ABS(DA(J,K,M))*XNI(1,1,J,K)-
     1        2.0D0*ABS(DB(I,K,M))*XNJ(1,1,K)-
     2        2.0D0*ABS(DC(I,J,M))*XNK(1,1)
      Q2(2,28)=VOL(I,J,K)*Q(2)+2.0D0*ABS(DB(I,K,M))*XNJ(2,1,K)+
     1        2.0D0*ABS(DC(I,J,M))*XNK(2,1)
      Q2(3,28)=VOL(I,J,K)*Q(3)+2.0D0*ABS(DA(J,K,M))*XNI(2,1,J,K)+
     1        2.0D0*ABS(DC(I,J,M))*XNK(1,2)
      Q2(4,28)=VOL(I,J,K)*Q(4)+2.0D0*ABS(DA(J,K,M))*XNI(1,2,J,K)+
     1        2.0D0*ABS(DB(I,K,M))*XNJ(1,2,K)
      Q2(5,28)=-VOL(I,J,K)*Q(5)-2.0D0*ABS(DC(I,J,M))*XNK(2,2)
      Q2(6,28)=-VOL(I,J,K)*Q(6)-2.0D0*ABS(DA(J,K,M))*XNI(2,2,J,K)
      Q2(7,28)=-VOL(I,J,K)*Q(7)-2.0D0*ABS(DB(I,K,M))*XNJ(2,2,K)
      Q2(8,28)=VOL(I,J,K)*Q(8)
      Q2(9,28)=-VOL(I,J,K)*Q(9)-CONST1*ABS(DA(J,K,M))*XNI(1,1,J,K)-
     1        2.0D0*ABS(DB(I,K,M))*XNJ(3,1,K)-
     2        2.0D0*ABS(DC(I,J,M))*XNK(3,1)
      Q2(10,28)=-VOL(I,J,K)*Q(10)-CONST1*ABS(DB(I,K,M))*XNJ(1,1,K)-
     1        2.0D0*ABS(DA(J,K,M))*XNI(3,1,J,K)-
     2        2.0D0*ABS(DC(I,J,M))*XNK(1,3)
      Q2(11,28)=-VOL(I,J,K)*Q(11)-2.0D0*ABS(DA(J,K,M))*XNI(1,3,J,K)-
     1        2.0D0*ABS(DB(I,K,M))*XNJ(1,3,K)-
     2        CONST1*ABS(DC(I,J,M))*XNK(1,1)
      Q2(12,28)=VOL(I,J,K)*Q(12)+CONST1*ABS(DA(J,K,M))*XNI(2,1,J,K)+
     1        2.0D0*ABS(DC(I,J,M))*XNK(3,2)
      Q2(13,28)=VOL(I,J,K)*Q(13)+2.0D0*ABS(DB(I,K,M))*XNJ(3,2,K)+
     1        CONST1*ABS(DA(J,K,M))*XNI(1,2,J,K)
      Q2(14,28)=-VOL(I,J,K)*Q(14)-CONST1*ABS(DA(J,K,M))*XNI(2,2,J,K)
      Q2(15,28)=VOL(I,J,K)*Q(15)+2.0D0*ABS(DA(J,K,M))*XNI(2,3,J,K)+
     1        CONST1*ABS(DC(I,J,M))*XNK(1,2)
      Q2(16,28)=VOL(I,J,K)*Q(16)+2.0D0*ABS(DB(I,K,M))*XNJ(2,3,K)+
     1        CONST1*ABS(DC(I,J,M))*XNK(2,1)
      Q2(17,28)=VOL(I,J,K)*Q(17)+CONST1*ABS(DB(I,K,M))*XNJ(2,1,K)+
     1        CONST1*ABS(DC(I,J,M))*XNK(2,3)
      Q2(18,28)=-VOL(I,J,K)*Q(18)-CONST1*ABS(DB(I,K,M))*XNJ(2,2,K)
      Q2(19,28)=VOL(I,J,K)*Q(19)+2.0D0*ABS(DA(J,K,M))*XNI(3,2,J,K)+
     1        CONST1*ABS(DB(I,K,M))*XNJ(1,2,K)
      Q2(20,28)=-VOL(I,J,K)*Q(20)-CONST1*ABS(DC(I,J,M))*XNK(2,2)
      Q2(21,28)=-VOL(I,J,K)*Q(21)-
     1        CONST1*ABS(DA(J,K,M))*XNI(3,1,J,K)-
     2        CONST1*ABS(DB(I,K,M))*XNJ(3,1,K)-
     3        2.0D0*CONST1*ABS(DC(I,J,M))*XNK(3,3)
      Q2(22,28)=-VOL(I,J,K)*Q(22)-2.0D0*ABS(DA(J,K,M))*XNI(3,3,J,K)
     1        -CONST1*ABS(DB(I,K,M))*XNJ(1,3,K)
     2        -CONST1*ABS(DC(I,J,M))*XNK(1,3)
      Q2(23,28)=-VOL(I,J,K)*Q(23)-
     1        CONST1*ABS(DA(J,K,M))*XNI(1,3,J,K)-
     2        2.0D0*ABS(DB(I,K,M))*XNJ(3,3,K)-
     3        CONST1*ABS(DC(I,J,M))*XNK(3,1)
      Q2(24,28)=VOL(I,J,K)*Q(24)+CONST1*ABS(DB(I,K,M))*XNJ(2,3,K)+
     1        CONST1*ABS(DC(I,J,M))*XNK(2,3)
      Q2(25,28)=VOL(I,J,K)*Q(25)+CONST1*ABS(DA(J,K,M))*XNI(2,3,J,K)+
     1        CONST1*ABS(DC(I,J,M))*XNK(3,2)
      Q2(26,28)=VOL(I,J,K)*Q(26)+CONST1*ABS(DA(J,K,M))*XNI(3,2,J,K)+
     1        CONST1*ABS(DB(I,K,M))*XNJ(3,2,K)
      Q2(27,28)=-VOL(I,J,K)*Q(27)-
     1        CONST1*ABS(DA(J,K,M))*XNI(3,3,J,K)-
     2        CONST1*ABS(DB(I,K,M))*XNJ(3,3,K)-
     3        CONST1*ABS(DC(I,J,M))*XNK(3,3)
      ENDIF
      DO 191 IEL=1,IELEM**3
      DO 190 JEL=IEL+1,IELEM**3     
      Q2(JEL,IEL)=Q2(IEL,JEL)         
 190  CONTINUE  
 191  CONTINUE  
      CALL ALSBD(IELEM**3,1,Q2,IER,IELEM**3)
      IF(IER.NE.0) CALL XABORT('SNFT13: SINGULAR MATRIX.')    
      IF(IELEM.EQ.1) THEN
      IF(LFIXUP.AND.(Q2(1,2).LE.RLOG)) Q2(1,2)=0.0
      XNI(1,1,J,K)=2.0D0*Q2(1,2)-XNI(1,1,J,K)
      XNJ(1,1,K)=2.0D0*Q2(1,2)-XNJ(1,1,K)
      XNK(1,1)=2.0D0*Q2(1,2)-XNK(1,1)      
      IF(LFIXUP.AND.(XNI(1,1,J,K).LE.RLOG)) XNI(1,1,J,K)=0.0
      IF(LFIXUP.AND.(XNJ(1,1,K).LE.RLOG)) XNJ(1,1,K)=0.0
      IF(LFIXUP.AND.(XNK(1,1).LE.RLOG)) XNK(1,1)=0.0   
      ELSE IF(IELEM.EQ.2) THEN  
      XNI(1,1,J,K)=XNI(1,1,J,K)+SIGN(1.0,DU(M))*CONST0*Q2(2,9)
      XNI(1,2,J,K)=XNI(1,2,J,K)+SIGN(1.0,DU(M))*CONST0*Q2(6,9)
      XNI(2,1,J,K)=XNI(2,1,J,K)+SIGN(1.0,DU(M))*CONST0*Q2(5,9)
      XNI(2,2,J,K)=XNI(2,2,J,K)+SIGN(1.0,DU(M))*CONST0*Q2(8,9)
      XNJ(1,1,K)=XNJ(1,1,K)+SIGN(1.0,DE(M))*CONST0*Q2(3,9)
      XNJ(1,2,K)=XNJ(1,2,K)+SIGN(1.0,DE(M))*CONST0*Q2(7,9)
      XNJ(2,1,K)=XNJ(2,1,K)+SIGN(1.0,DE(M))*CONST0*Q2(5,9)
      XNJ(2,2,K)=XNJ(2,2,K)+SIGN(1.0,DE(M))*CONST0*Q2(8,9)  
      XNK(1,1)=XNK(1,1)+SIGN(1.0,DZ(M))*CONST0*Q2(4,9)
      XNK(1,2)=XNK(1,2)+SIGN(1.0,DZ(M))*CONST0*Q2(7,9)
      XNK(2,1)=XNK(2,1)+SIGN(1.0,DZ(M))*CONST0*Q2(6,9)
      XNK(2,2)=XNK(2,2)+SIGN(1.0,DZ(M))*CONST0*Q2(8,9)
      ELSE IF(IELEM.EQ.3) THEN
      XNI(1,1,J,K)=2.0D0*Q2(1,28)+CONST1*Q2(9,28)-XNI(1,1,J,K)
      XNI(2,1,J,K)=2.0D0*Q2(3,28)+CONST1*Q2(12,28)-XNI(2,1,J,K)
      XNI(1,2,J,K)=2.0D0*Q2(4,28)+CONST1*Q2(13,28)-XNI(1,2,J,K)
      XNI(2,2,J,K)=2.0D0*Q2(6,28)+CONST1*Q2(14,28)-XNI(2,2,J,K)
      XNI(1,3,J,K)=2.0D0*Q2(11,28)+CONST1*Q2(23,28)-XNI(1,3,J,K)
      XNI(3,1,J,K)=2.0D0*Q2(10,28)+CONST1*Q2(21,28)-XNI(3,1,J,K)
      XNI(3,2,J,K)=2.0D0*Q2(19,28)+CONST1*Q2(26,28)-XNI(3,2,J,K)
      XNI(2,3,J,K)=2.0D0*Q2(15,28)+CONST1*Q2(25,28)-XNI(2,3,J,K)
      XNI(3,3,J,K)=2.0D0*Q2(22,28)+CONST1*Q2(27,28)-XNI(3,3,J,K)     
      XNJ(1,1,K)=2.0D0*Q2(1,28)+CONST1*Q2(10,28)-XNJ(1,1,K)
      XNJ(2,1,K)=2.0D0*Q2(2,28)+CONST1*Q2(17,28)-XNJ(2,1,K)
      XNJ(1,2,K)=2.0D0*Q2(4,28)+CONST1*Q2(19,28)-XNJ(1,2,K)
      XNJ(2,2,K)=2.0D0*Q2(7,28)+CONST1*Q2(18,28)-XNJ(2,2,K)
      XNJ(3,1,K)=2.0D0*Q2(9,28)+CONST1*Q2(21,28)-XNJ(3,1,K)
      XNJ(1,3,K)=2.0D0*Q2(11,28)+CONST1*Q2(22,28)-XNJ(1,3,K)
      XNJ(3,2,K)=2.0D0*Q2(13,28)+CONST1*Q2(26,28)-XNJ(3,2,K)
      XNJ(2,3,K)=2.0D0*Q2(16,28)+CONST1*Q2(24,28)-XNJ(2,3,K)
      XNJ(3,3,K)=2.0D0*Q2(23,28)+CONST1*Q2(27,28)-XNJ(3,3,K)        
      XNK(1,1)=2.0D0*Q2(1,28)+CONST1*Q2(11,28)-XNK(1,1)
      XNK(2,1)=2.0D0*Q2(2,28)+CONST1*Q2(16,28)-XNK(2,1)
      XNK(1,2)=2.0D0*Q2(3,28)+CONST1*Q2(15,28)-XNK(1,2)
      XNK(2,2)=2.0D0*Q2(5,28)+CONST1*Q2(20,28)-XNK(2,2)
      XNK(3,1)=2.0D0*Q2(9,28)+CONST1*Q2(23,28)-XNK(3,1)
      XNK(1,3)=2.0D0*Q2(10,28)+CONST1*Q2(22,28)-XNK(1,3)
      XNK(3,2)=2.0D0*Q2(12,28)+CONST1*Q2(25,28)-XNK(3,2)
      XNK(2,3)=2.0D0*Q2(17,28)+CONST1*Q2(24,28)-XNK(2,3)
      XNK(3,3)=2.0D0*Q2(21,28)+CONST1*Q2(27,28)-XNK(3,3)
      ENDIF
      DO 201 L=1,NSCT
      DO 200 IEL=1,IELEM**3
      FLUX(IEL,L,I,J,K)=FLUX(IEL,L,I,J,K)+
     1  W(M)*REAL(Q2(IEL,IELEM**3+1))*PL(L,M)
  200 CONTINUE
  201 CONTINUE
*
  140 CONTINUE
      DO 151 IEL=1,IELEM
      DO 150 JEL=1,IELEM
      XNEK(IEL,JEL,I,J,M)=REAL(XNK(IEL,JEL))
  150 CONTINUE
  151 CONTINUE
*
  152 CONTINUE
      DO 162 IEL=1,IELEM
      DO 161 JEL=1,IELEM
      DO 160 K=1,LZ
      XNEJ(IEL,JEL,I,K,M)=REAL(XNJ(IEL,JEL,K))
  160 CONTINUE
  161 CONTINUE
  162 CONTINUE
*
  163 CONTINUE
      DO 183 J=1,LY
      DO 182 K=1,LZ
      DO 181 IEL=1,IELEM
      DO 180 JEL=1,IELEM
      XNEI(IEL,JEL,J,K,M)=REAL(XNI(IEL,JEL,J,K))
  180 CONTINUE
  181 CONTINUE
  182 CONTINUE
  183 CONTINUE
*
  170 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(XNJ,XNI)
      RETURN
      END
