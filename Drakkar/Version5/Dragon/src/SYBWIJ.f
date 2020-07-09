*DECK SYBWIJ
      SUBROUTINE SYBWIJ (NREG,MAXPTS,SIGW,PIJ)
*
*-----------------------------------------------------------------------
*
*Purpose:
* scattering reduction for collision probabilities.
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
* NREG    total number of regions.
* MAXPTS  first dimension of matrix PIJ.
* SIGW    P0 within-group scattering macroscopic cross sections
*         ordered by volume.
* PIJ     reduced collision probability matrix.
*
*Parameters: output
*  PIJ    scattering-reduced collision probability matrix.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NREG,MAXPTS
      REAL SIGW(NREG),PIJ(MAXPTS,NREG)
*----
*  LOCAL VARIABLES
*----
      REAL, ALLOCATABLE, DIMENSION(:,:) :: WIJ
*
      ALLOCATE(WIJ(NREG,2*NREG))
      DO 20 I=1,NREG
      DO 10 J=1,NREG
      WIJ(I,NREG+J)=PIJ(I,J)
   10 WIJ(I,J)=-PIJ(I,J)*SIGW(J)
   20 WIJ(I,I)=1.0+WIJ(I,I)
      CALL ALSB(NREG,NREG,WIJ,IER,NREG)
      IF(IER.NE.0) CALL XABORT('SYBWIJ: SINGULAR MATRIX.')
      DO 30 J=1,NREG
      DO 30 I=1,NREG
   30 PIJ(I,J)=WIJ(I,NREG+J)
      DEALLOCATE(WIJ)
      RETURN
      END
