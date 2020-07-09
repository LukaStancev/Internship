*DECK CHAB04
      SUBROUTINE CHAB04(IPLIB,IMPX,IRHS,NGRP,NLEG,IMOD,IL,IGM,IGP,
     1 VALUE,DELTA,FMULT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* modify scattering information in a microlib or in a draglib.
*
*Copyright:
* Copyright (C) 2011 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPLIB   LCM pointer to the Microlib or Draglib.
* IMPX    print index.
* IRHS    type of IPLIB: =1: Microlib; =2: Draglib.
* NGRP    number of energy groups.
* NLEG    max Legendre order of scattering anisotropy (1=isotropic,
*         etc.).
* IMOD    type of modification: =1,2: replace the value; =3: increase by
*         VALUE; =4: multiply by VALUE.
* IL      Legendre order under consideration.
* IGM     first energy group to modify.
* IGP     last energy group to modify.
* VALUE   value array used in scattering modification operation.
*
*Parameters: output
* DELTA   difference in scattering cross section value.
* FMULT   multiplicative modification factors.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB
      INTEGER IMPX,IRHS,NGRP,NLEG,IMOD,IL,IGM,IGP
      REAL VALUE(NGRP),DELTA(NGRP),FMULT(NGRP)
*----
*  LOCAL VARIABLES
*----
      CHARACTER CM*2
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NJJ,IJJ,ITYPRO
      REAL, ALLOCATABLE, DIMENSION(:) :: GAR1,XSSCM
      REAL, ALLOCATABLE, DIMENSION(:,:) :: GAR2
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(NJJ(NGRP),IJJ(NGRP))
      ALLOCATE(GAR1(NGRP),GAR2(NGRP,NGRP))
*      
      CALL XDRSET(DELTA,NGRP,0.0)
      CALL XDRSET(FMULT,NGRP,1.0)
*----
*  RECOVER SCATTERING INFORMATION
*----
      IF(IL.GE.NLEG) CALL XABORT('CHAB04: LEGENDRE INDEX OVERFLOW.')
      IF(IRHS.EQ.1) THEN
         ALLOCATE(ITYPRO(NLEG))
         CALL XDRLGS(IPLIB,-1,IMPX,IL,IL,1,NGRP,GAR1,GAR2,ITYPRO)
         DEALLOCATE(ITYPRO)
      ELSE IF(IRHS.EQ.2) THEN
         ALLOCATE(XSSCM(NGRP*NGRP))
         WRITE (CM,'(I2.2)') IL
         CALL XDRSET(GAR1,NGRP,0.0)
         CALL XDRSET(GAR2,NGRP*NGRP,0.0)
         CALL LCMGET(IPLIB,'NJJS'//CM,NJJ)
         CALL LCMGET(IPLIB,'IJJS'//CM,IJJ)
         LENGT=0
         DO 10 I=1,NGRP
   10    LENGT=LENGT+NJJ(I)
         CALL XDRSET(XSSCM,LENGT,0.0)
         CALL LCMGET(IPLIB,'SCAT'//CM,XSSCM)
         IGAR=0
         DO 20 IG2=1,NGRP
         DO 20 IG1=IJJ(IG2),IJJ(IG2)-NJJ(IG2)+1,-1
         IGAR=IGAR+1
         GAR2(IG2,IG1)=XSSCM(IGAR)
   20    GAR1(IG1)=GAR1(IG1)+GAR2(IG2,IG1)
         DEALLOCATE(XSSCM)
      ENDIF
*----
*  MODIFY SCATTERING INFORMATION
*----
      DO 40 IG2=IGM,IGP
      IF(GAR1(IG2).NE.0.0) THEN
         DO 30 IG1=1,NGRP  
   30    GAR2(IG1,IG2)=GAR2(IG1,IG2)/GAR1(IG2)
      ENDIF
      IF((IMOD.EQ.1).OR.(IMOD.EQ.2)) THEN
         IF(GAR1(IG2).EQ.0.0) THEN
            FMULT(IG2)=1.0
         ELSE
            FMULT(IG2)=VALUE(IG2)/GAR1(IG2)
         ENDIF
         DELTA(IG2)=VALUE(IG2)-GAR1(IG2)
         GAR1(IG2)=VALUE(IG2)
      ELSE IF(IMOD.EQ.3) THEN
         IF(GAR1(IG2).EQ.0.0) THEN
            FMULT(IG2)=1.0
         ELSE
            FMULT(IG2)=1.0+VALUE(IG2)/GAR1(IG2)
         ENDIF
         DELTA(IG2)=VALUE(IG2)
         GAR1(IG2)=GAR1(IG2)+VALUE(IG2)
      ELSE IF(IMOD.EQ.4) THEN
         FMULT(IG2)=VALUE(IG2)
         DELTA(IG2)=GAR1(IG2)*(VALUE(IG2)-1.0)
         GAR1(IG2)=GAR1(IG2)*VALUE(IG2)
      ENDIF
      DO 40 IG1=1,NGRP
   40 GAR2(IG1,IG2)=GAR2(IG1,IG2)*GAR1(IG2)
*----
*  SAVE SCATTERING INFORMATION
*----
      IF(IRHS.EQ.1) THEN
         ALLOCATE(ITYPRO(NLEG))
         CALL XDRLGS(IPLIB,1,IMPX,IL,IL,1,NGRP,GAR1,GAR2,ITYPRO)
         DEALLOCATE(ITYPRO)
      ELSE IF(IRHS.EQ.2) THEN
         ALLOCATE(XSSCM(NGRP*NGRP))
         IGAR=0
         DO 50 IG2=1,NGRP
         DO 50 IG1=IJJ(IG2),IJJ(IG2)-NJJ(IG2)+1,-1
         IGAR=IGAR+1
   50    XSSCM(IGAR)=GAR2(IG2,IG1)
         CALL LCMPUT(IPLIB,'SCAT'//CM,IGAR,2,XSSCM)
         DEALLOCATE(XSSCM)
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(GAR2,GAR1)
      DEALLOCATE(IJJ,NJJ)
      RETURN
      END
