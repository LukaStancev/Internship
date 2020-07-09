*DECK LIBPRI
      SUBROUTINE LIBPRI(MAXTRA,DELI,AWR,IALTER,IL,N,PRI)
*
*-----------------------------------------------------------------------
*
*Purpose:
* compute the PRI array for various Legendre orders using Gaussian
* integration.
*
*Copyright:
* Copyright (C) 2003 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* MAXTRA  allocated dimension of array PRI.
* DELI    elementary lethargy width of the equi-width lethargy mesh.
* AWR     mass ratio for current isotope.
* IALTER  type of approximation (=0: use exponentials; =1: use Taylor
*         expansions).
* IL      Legendre order (=0: isotropic kernel).
*
*Parameters: output
* N       related to the exact dimension of array PRI.
* PRI     array containing the slowing-down probabilities defined on
*         an equi-width lethargy mesh.
*
*-----------------------------------------------------------------------
*
* Reference: M. Grandotto-Bietolli, "AUTOSECOL, un calcul automatique
*   de l'autoprotection des resonances des isotopes lourds",
*   Commissariat a l'Energie Atomique, Note CEA-N-1961, Mars 1977.
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER MAXTRA,IALTER,IL,N
      REAL DELI,AWR,PRI(MAXTRA)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NGPT=6,MAXNL=50)
      DOUBLE PRECISION AWRB,ALP,T0,FACT
      REAL UI(NGPT),WI(NGPT),UJ(NGPT),WJ(NGPT)
      REAL POLY(0:MAXNL),CALC(0:MAXNL,0:2)
      ZMU(AWR,U)=0.5*(AWR+1.0)*EXP(-0.5*U)-0.5*(AWR-1.0)*EXP(0.5*U)
*----
*  COMPUTE THE LEGENDRE POLYNOMIAL OF ORDER IL.
*----
      IF(IL.GT.MAXNL) CALL XABORT('LIBPRI: IL OVERFLOW.')
      IF(IL.EQ.0) THEN
         POLY(0)=1.0
      ELSE IF(IL.EQ.1) THEN
         POLY(0)=0.0
         POLY(1)=1.0
      ELSE
         CALL XDRSET(CALC(0,0),IL+1,0.0)
         CALL XDRSET(CALC(0,1),IL+1,0.0)
         CALC(0,0)=1.0
         CALC(1,1)=1.0
         DO 10 J=2,IL
         DO 10 I=0,IL
         T0=-REAL(J-1)*CALC(I,MOD(J-2,3))
         IF(I.GT.0) T0=T0+(2.0*REAL(J-1)+1.0)*CALC(I-1,MOD(J-1,3))
   10    CALC(I,MOD(J,3))=REAL(T0)/REAL(J)
         DO 20 I=0,IL
   20    POLY(I)=CALC(I,MOD(IL,3))
      ENDIF
*
      AWRB=AWR
      IF(AWR.LT.1.0001) AWRB=1.0001
      ALP=((AWRB-1.D0)/(AWRB+1.D0))**2
      REPS=REAL(-DLOG(ALP))
      N=INT(REPS/DELI)
      IF(N+2.GT.MAXTRA) CALL XABORT('LIBPRI: MAXTRA IS TOO SMALL.')
*
      IF(N.EQ.0) THEN
*        COMPUTE PRI(1).
         PRI(1)=0.0
         CALL ALGPT(NGPT,0.0,DELI-REPS,UI,WI)
         DO 40 I=1,NGPT
         CALL ALGPT(NGPT,UI(I),UI(I)+REPS,UJ,WJ)
         DO 40 J=1,NGPT
         FACT=POLY(0)
         T0=1.0D0
         DO 30 K=1,IL
         T0=T0*ZMU(AWR,UJ(J)-UI(I))
   30    FACT=FACT+POLY(K)*T0
         IF(IALTER.EQ.0) THEN
            PRI(1)=PRI(1)+WI(I)*WJ(J)*EXP(UI(I)-UJ(J))*REAL(FACT)
         ELSE
            PRI(1)=PRI(1)+WI(I)*WJ(J)*REAL(FACT)
         ENDIF
   40    CONTINUE
         CALL ALGPT(NGPT,DELI-REPS,DELI,UI,WI)
         DO 50 I=1,NGPT
         CALL ALGPT(NGPT,UI(I),DELI,UJ,WJ)
         DO 50 J=1,NGPT
         FACT=POLY(0)
         T0=1.0D0
         DO 45 K=1,IL
         T0=T0*ZMU(AWR,UJ(J)-UI(I))
   45    FACT=FACT+POLY(K)*T0
         IF(IALTER.EQ.0) THEN
            PRI(1)=PRI(1)+WI(I)*WJ(J)*EXP(UI(I)-UJ(J))*REAL(FACT)
         ELSE
            PRI(1)=PRI(1)+WI(I)*WJ(J)*REAL(FACT)
         ENDIF
   50    CONTINUE
*
*        COMPUTE PRI(2).
         PRI(2)=0.0
         CALL ALGPT(NGPT,DELI-REPS,DELI,UI,WI)
         DO 60 I=1,NGPT
         CALL ALGPT(NGPT,DELI,UI(I)+REPS,UJ,WJ)
         DO 60 J=1,NGPT
         FACT=POLY(0)
         T0=1.0D0
         DO 55 K=1,IL
         T0=T0*ZMU(AWR,UJ(J)-UI(I))
   55    FACT=FACT+POLY(K)*T0
         IF(IALTER.EQ.0) THEN
            PRI(2)=PRI(2)+WI(I)*WJ(J)*EXP(UI(I)-UJ(J))*REAL(FACT)
         ELSE
            PRI(2)=PRI(2)+WI(I)*WJ(J)*REAL(FACT)
         ENDIF
   60    CONTINUE
      ELSE
*        COMPUTE PRI(1).
         PRI(1)=0.0
         CALL ALGPT(NGPT,0.0,DELI,UI,WI)
         DO 70 I=1,NGPT
         CALL ALGPT(NGPT,REAL(UI(I)),DELI,UJ,WJ)
         DO 70 J=1,NGPT
         FACT=POLY(0)
         T0=1.0D0
         DO 65 K=1,IL
         T0=T0*ZMU(AWR,UJ(J)-UI(I))
   65    FACT=FACT+POLY(K)*T0
         IF(IALTER.EQ.0) THEN
            PRI(1)=PRI(1)+WI(I)*WJ(J)*EXP(UI(I)-UJ(J))*REAL(FACT)
         ELSE
            PRI(1)=PRI(1)+WI(I)*WJ(J)*REAL(FACT)
         ENDIF
   70    CONTINUE
*
*        COMPUTE PRI(L) FOR L=2,N.
         CALL ALGPT(NGPT,0.0,DELI,UI,WI)
         DO 80 L=2,N
         PRI(L)=0.0
         DO 80 I=1,NGPT
         CALL ALGPT(NGPT,REAL(L-1)*DELI,REAL(L)*DELI,UJ,WJ)
         DO 80 J=1,NGPT
         FACT=POLY(0)
         T0=1.0D0
         DO 75 K=1,IL
         T0=T0*ZMU(AWR,UJ(J)-UI(I))
   75    FACT=FACT+POLY(K)*T0
         IF(IALTER.EQ.0) THEN
            PRI(L)=PRI(L)+WI(I)*WJ(J)*EXP(UI(I)-UJ(J))*REAL(FACT)
         ELSE
            PRI(L)=PRI(L)+WI(I)*WJ(J)*REAL(FACT)
         ENDIF
   80    CONTINUE
*
*        COMPUTE PRI(N+1).
         PRI(N+1)=0.0
         CALL ALGPT(NGPT,0.0,REAL(N+1)*DELI-REPS,UI,WI)
         DO 90 I=1,NGPT
         CALL ALGPT(NGPT,REAL(N)*DELI,UI(I)+REPS,UJ,WJ)
         DO 90 J=1,NGPT
         FACT=POLY(0)
         T0=1.0D0
         DO 85 K=1,IL
         T0=T0*ZMU(AWR,UJ(J)-UI(I))
   85    FACT=FACT+POLY(K)*T0
         IF(IALTER.EQ.0) THEN
            PRI(N+1)=PRI(N+1)+WI(I)*WJ(J)*EXP(UI(I)-UJ(J))*REAL(FACT)
         ELSE
            PRI(N+1)=PRI(N+1)+WI(I)*WJ(J)*REAL(FACT)
         ENDIF
   90    CONTINUE
         CALL ALGPT(NGPT,REAL(N+1)*DELI-REPS,DELI,UI,WI)
         DO 100 I=1,NGPT
         CALL ALGPT(NGPT,REAL(N)*DELI,REAL(N+1)*DELI,UJ,WJ)
         DO 100 J=1,NGPT
         FACT=POLY(0)
         T0=1.0D0
         DO 95 K=1,IL
         T0=T0*ZMU(AWR,UJ(J)-UI(I))
   95    FACT=FACT+POLY(K)*T0
         IF(IALTER.EQ.0) THEN
            PRI(N+1)=PRI(N+1)+WI(I)*WJ(J)*EXP(UI(I)-UJ(J))*REAL(FACT)
         ELSE
            PRI(N+1)=PRI(N+1)+WI(I)*WJ(J)*REAL(FACT)
         ENDIF
  100    CONTINUE
*
*        COMPUTE PRI(N+2).
         PRI(N+2)=0.0
         CALL ALGPT(NGPT,REAL(N+1)*DELI-REPS,DELI,UI,WI)
         DO 110 I=1,NGPT
         CALL ALGPT(NGPT,REAL(N+1)*DELI,UI(I)+REPS,UJ,WJ)
         DO 110 J=1,NGPT
         FACT=POLY(0)
         T0=1.0D0
         DO 105 K=1,IL
         T0=T0*ZMU(AWR,UJ(J)-UI(I))
  105    FACT=FACT+POLY(K)*T0
         IF(IALTER.EQ.0) THEN
            PRI(N+2)=PRI(N+2)+WI(I)*WJ(J)*EXP(UI(I)-UJ(J))*REAL(FACT)
         ELSE
            PRI(N+2)=PRI(N+2)+WI(I)*WJ(J)*REAL(FACT)
         ENDIF
  110    CONTINUE
      ENDIF
      IF(IALTER.EQ.0) THEN
         DO 120 I=1,N+2
  120    PRI(I)=PRI(I)/DELI/REAL(1.0D0-ALP)
      ELSE
         DO 130 I=1,N+2
  130    PRI(I)=PRI(I)/DELI/REPS
      ENDIF
      RETURN
      END
