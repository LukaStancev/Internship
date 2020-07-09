*DECK SYBRHL
      SUBROUTINE SYBRHL(IPRT,NSOUT,NREG,G,PROB)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Stamm'ler normalisation of collision, escape and transmission
* probabilities.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Parameters: input
* IPRT    print parameter (equal to zero for no print).
* NSOUT   number of surfaces.
* NREG    number of regions.
* G       surface and volume array (first NSOUT values: surface/4;
*         last NREG values: volume * total macroscopic cross section).
* PROB    collision, escape and transmission probabilities.
*
*Parameters: output
* G       renormalization factors.
* PROB    normalized collision, escape and transmission probabilities.
*
*Reference:
* 'Normalization techniques for CP matrices', conf/Physor-90,
* Marseille/France, v 2, p ix-40 (1990).
* 'Helios: Angularly dependent collision probabilities' E.A. Villarino,
* R.J.J.Stamm'ler and A.A.Ferri and J.J.Casal, Nucl.Sci.Eng. 112,16-31,
* 1992.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IPRT,NSOUT,NREG
      REAL G(NSOUT+NREG),PROB((NSOUT+NREG)*(NSOUT*NREG+1)/2)
*----
*  LOCAL VARIABLES
*----
      INTEGER CPTLB,CPTAC,CTOT
      PARAMETER (CPTLB=3, CPTAC=3, CTOT=CPTAC+CPTLB, IUNOUT=6,
     1 EPSCON=1.0E-6, NITMAX=20)
      DOUBLE PRECISION WFSPAD,WFSP,NOM,DENOM,DMU
      LOGICAL NOTCON
      REAL, ALLOCATABLE, DIMENSION(:,:) :: WEIG
*
      INDPOS(I,J)=MAX(I,J)*(MAX(I,J)-1)/2+MIN(I,J)
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(WEIG(NSOUT+NREG,3))
*----
*  INITIALISATION OF WEIGHTS
*----
      NOTCON=.FALSE.
      DO 10 IR=1,NSOUT+NREG
      WEIG(IR,1)=0.0
      WEIG(IR,2)=0.5
   10 WEIG(IR,3)=0.5
*----
*  MAIN ITERATION LOOP
*----
      IF(IPRT.GT.8) THEN
         WRITE(IUNOUT,'(/30H SYBRHL: NORMALIZATION FACTORS)')
         WRITE(IUNOUT,'(1X,A24)')  'ITER.     MU      ERROR '
      ENDIF
      NIT=0
   20 NIT=NIT+1
      IF(NIT.GT.NITMAX) THEN
         NOTCON=.TRUE.
         WRITE(IUNOUT,'(31H SYBRHL: WEIGHTS NOT CONVERGED.)')
         GO TO 80
      ENDIF
      DO 40 IR=1,NSOUT+NREG
      WFSPAD=G(IR)+PROB(INDPOS(IR,IR))*WEIG(IR,3)
      WFSP=PROB(INDPOS(IR,IR))
      DO 30 JR=1,NSOUT+NREG
      WFSPAD=WFSPAD-WEIG(JR,3)*PROB(INDPOS(IR,JR))
   30 WFSP=WFSP+PROB(INDPOS(IR,JR))
   40 WEIG(IR,3)=REAL(WFSPAD/WFSP)
*----
*  ACCELERATION BY RESIDUAL MINIMIZATION
*----
      IF(MOD(NIT-1,CTOT).GE.CPTAC) THEN
         NOM   = 0.0D0
         DENOM = 0.0D0
         DO 50 IR=1,NSOUT+NREG
         R1= WEIG(IR,2) - WEIG(IR,1)
         R2= WEIG(IR,3) - WEIG(IR,2)
         NOM = NOM + R1*(R2-R1)
   50    DENOM = DENOM + (R2-R1)*(R2-R1)
         IF(DENOM.EQ.0.0D0) THEN
           DMU=1.0D0
         ELSE
           DMU=-NOM/DENOM
         ENDIF
         ZMU=REAL(DMU)
         IF(ZMU.GT.10.0 .OR. ZMU.LT.0.0) THEN
            IF( IPRT.GT.2 ) WRITE(IUNOUT,'(I3,1P,G12.4,A)') NIT,ZMU,
     >      ' =MU / SYBRHL: NON ACCELERATION'
            ZMU=1.0
         ENDIF
         DO 60 IR=1,NSOUT+NREG
         WEIG(IR,3)=WEIG(IR,2)+ZMU*(WEIG(IR,3)-WEIG(IR,2))
   60    WEIG(IR,2)=WEIG(IR,1)+ZMU*(WEIG(IR,2)-WEIG(IR,1))
      ELSE
         ZMU = 1.0
      ENDIF
*----
*  CALCULATIONS OF SQUARE DISTANCE BETWEEN 2 ITERATIONS AND UPDATING
*  OF THE SOLUTION
*----
      TOTCON = 0.0
      DO 70 IR=1,NSOUT+NREG
      TMPCON=ABS(WEIG(IR,3)-WEIG(IR,2))/WEIG(IR,3)
      TOTCON=MAX(TMPCON,TOTCON)
      WEIG(IR,1)=WEIG(IR,2)
   70 WEIG(IR,2)=WEIG(IR,3)
      IF(IPRT.GT.8) WRITE(IUNOUT,'(I3,F9.5,E15.7)') NIT,ZMU,TOTCON
*----
*  CONVERGENCE TEST
*----
      IF(TOTCON.LT.EPSCON) GO TO 80
      GO TO 20
*----
*  RENORMALIZE "PIJ" SYMMETRIC MATRIX
*----
   80 IPRB=0
      DO 90 IR=1,NSOUT+NREG
      G(IR)=WEIG(IR,1)
      DO 90 JR=1,IR
      IPRB=IPRB+1
   90 PROB(IPRB)=PROB(IPRB)*(WEIG(IR,1)+WEIG(JR,1))
*----
*  PRINT WEIGHT FACTORS IF THERE IS A PROBLEM
*----
      IF(NOTCON .OR. (IPRT.GE.15)) THEN
         WRITE(IUNOUT,'(24H SURFACE WEIGHTS FACTORS/)')
         NSURC=1
         DO 100 IP =1,(9+NSOUT)/10
         NSURM=MIN(NSOUT,NSURC+9)
         WRITE(IUNOUT,'(10X,10(A5,I6)/)')
     >                 (' SUR ',IR,IR= NSURC, NSURM)
         WRITE(IUNOUT,'(10H WEIGHT   ,10F11.5)')
     >                 (WEIG(IR,1),IR=NSURC,NSURM)
 100     NSURC=NSURC+10
         WRITE(IUNOUT,'(24H  VOLUME WEIGHTS FACTORS/)')
         NVOLC=NSOUT+1
         DO 110 IP=1,(9+NREG)/10
         NVOLM=MIN(NSOUT+NREG,NVOLC+9)
         WRITE(IUNOUT,'(10X,10(A5,I6)/)')
     >                 (' VOL ',IR,IR=NVOLC-NSOUT,NVOLM-NSOUT)
         WRITE(IUNOUT,'(10H WEIGHT   ,10F11.5)')
     >                 (WEIG(IR,1),IR=NVOLC,NVOLM)
 110     NVOLC=NVOLC+10
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(WEIG)
      RETURN
      END
