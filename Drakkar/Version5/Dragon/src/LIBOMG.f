*DECK LIBOMG
      SUBROUTINE LIBOMG(MX,IX,X,MY,IY,Y,DCM,OMEG)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the correlated weight matrix preserving a matrix of moments.
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
* MX      number of base points in the first variable.
* IX      order of the first moment of the first variable. We must
*         have 1-MX <= IX <= 0 (order 0 moment is always preserved).
* X       base points in the first variable.
* MY      number of base points in the second variable.
* IY      order of the first moment of the second variable. We must
*         have 1-MY <= IY <= 0 (order 0 moment is always preserved).
* Y       base points in the second variable.
* DCM     comoments.
*
*Parameters: output
* OMEG    correlated weight matrix.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      REAL X(MX),Y(MY),OMEG(MX,MY)
      DOUBLE PRECISION DCM(MX,MY)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (MAXNOR=20)
      DOUBLE PRECISION DD,DAUX,WORK,DDA(0:MAXNOR),PROD1(MAXNOR,MAXNOR),
     1 PROD2(MAXNOR,MAXNOR)
*
      IF(MX.GT.MAXNOR) CALL XABORT('LIBOMG: STORAGE OVERFLOW(1).')
      IF(MY.GT.MAXNOR) CALL XABORT('LIBOMG: STORAGE OVERFLOW(2).')
      IF((1-MX.GT.IX).OR.(IX.GT.0)) CALL XABORT('LIBOMG: INCONSISTENT '
     1 //'VALUE OF IX.')
      IF((1-MY.GT.IY).OR.(IY.GT.0)) CALL XABORT('LIBOMG: INCONSISTENT '
     1 //'VALUE OF IY.')
*
      DO 10 I=1,MX
      DO 10 J=1,MY
   10 PROD1(I,J)=0.0D0
      DO 50 I=1,MX
      DAUX=DBLE(X(I))
      DDA(0)=1.0D0
      J0=0
      DO 30 J=1,MX
      IF(J.EQ.I) GO TO 30
      J0=J0+1
      DDA(J0)=DDA(J0-1)
      DO 20 K=1,J0-1
   20 DDA(J0-K)=DDA(J0-K-1)-DDA(J0-K)*DBLE(X(J))
      DDA(0)=-DDA(0)*DBLE(X(J))
   30 CONTINUE
      DD=1.0D0
      DO 40 J=1,MX
      IF(J.NE.I) DD=DD*(DBLE(X(J))-DAUX)
   40 CONTINUE
      WORK=((-1.0D0)**(MX-1))*DAUX**(-IX)/DD
      DO 50 J=1,MY
      DO 50 K=1,MX
   50 PROD1(I,J)=PROD1(I,J)+WORK*DDA(K-1)*DCM(K,J)
*
      DO 60 I=1,MX
      DO 60 J=1,MY
   60 PROD2(I,J)=0.0D0
      DO 100 I=1,MY
      DAUX=DBLE(Y(I))
      DDA(0)=1.0D0
      J0=0
      DO 80 J=1,MY
      IF(J.EQ.I) GO TO 80
      J0=J0+1
      DDA(J0)=DDA(J0-1)
      DO 70 K=1,J0-1
   70 DDA(J0-K)=DDA(J0-K-1)-DDA(J0-K)*DBLE(Y(J))
      DDA(0)=-DDA(0)*DBLE(Y(J))
   80 CONTINUE
      DD=1.0D0
      DO 90 J=1,MY
      IF(J.NE.I) DD=DD*(DBLE(Y(J))-DAUX)
   90 CONTINUE
      WORK=((-1.0D0)**(MY-1))*DAUX**(-IY)/DD
      DO 100 J=1,MX
      DO 100 K=1,MY
  100 PROD2(J,I)=PROD2(J,I)+WORK*DDA(K-1)*PROD1(J,K)
*
      DO 120 I=1,MX
      DO 120 J=1,MY
  120 OMEG(I,J)=REAL(PROD2(I,J))
      RETURN
      END
