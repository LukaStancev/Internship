*DECK PNSH
      FUNCTION PNSH(L,M,ZMU,ETA,XI)
*
*-----------------------------------------------------------------------
*
*Purpose:
* return the real spherical harmonics corresponding to a set of
* direction cosines.
*
*Copyright:
* Copyright (C) 2004 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* L       Legendre order.
* M       azimuthal order.
* ZMU     X-directed direction cosine.
* ETA     Y-directed direction cosine.
* XI      Z-directed direction cosine.
*
*Parameters: output
* PNSH    value of the spherical harmonics.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER L,M
      REAL ZMU,ETA,XI
*
      TEST=ZMU*ZMU+ETA*ETA+XI*XI
      IF(ABS(TEST-1.0).GT.1.0E-5) THEN
         CALL XABORT('PNSH: INVALID DIRECTION COSINES.')
      ENDIF
      PNSH=0.0
      IF((L.EQ.0).AND.(M.EQ.0)) THEN
         PNSH=1.0
      ELSE IF((L.EQ.1).AND.(M.EQ.-1)) THEN
         PNSH=XI
      ELSE IF((L.EQ.1).AND.(M.EQ.0)) THEN
         PNSH=ZMU
      ELSE IF((L.EQ.1).AND.(M.EQ.1)) THEN
         PNSH=ETA
      ELSE IF((L.EQ.2).AND.(M.EQ.-2)) THEN
         PNSH=SQRT(3.0)*ETA*XI
      ELSE IF((L.EQ.2).AND.(M.EQ.-1)) THEN
         PNSH=SQRT(3.0)*ZMU*XI
      ELSE IF((L.EQ.2).AND.(M.EQ.0)) THEN
         PNSH=0.5*(3.0*ZMU*ZMU-1.0)
      ELSE IF((L.EQ.2).AND.(M.EQ.1)) THEN
         PNSH=SQRT(3.0)*ZMU*ETA
      ELSE IF((L.EQ.2).AND.(M.EQ.2)) THEN
         PNSH=0.5*SQRT(3.0)*(ETA*ETA-XI*XI)
      ELSE IF((L.EQ.3).AND.(M.EQ.-3)) THEN
         PNSH=SQRT(5./8.)*XI*(3.0*ETA*ETA-XI*XI)
      ELSE IF((L.EQ.3).AND.(M.EQ.-2)) THEN
         PNSH=SQRT(15.0)*ETA*XI*ZMU
      ELSE IF((L.EQ.3).AND.(M.EQ.-1)) THEN
         PNSH=SQRT(3./8.)*XI*(5.0*ZMU*ZMU-1.0)
      ELSE IF((L.EQ.3).AND.(M.EQ.0)) THEN
         PNSH=0.5*ZMU*(5.0*ZMU*ZMU-3.0)
      ELSE IF((L.EQ.3).AND.(M.EQ.1)) THEN
         PNSH=SQRT(3./8.)*ETA*(5.0*ZMU*ZMU-1.0)
      ELSE IF((L.EQ.3).AND.(M.EQ.2)) THEN
         PNSH=SQRT(15.0/4.0)*ZMU*(ETA*ETA-XI*XI)
      ELSE IF((L.EQ.3).AND.(M.EQ.3)) THEN
         PNSH=SQRT(5./8.)*ETA*(ETA*ETA-3.0*XI*XI)
      ELSE IF((L.EQ.4).AND.(M.EQ.-4)) THEN
         PNSH=0.5*SQRT(35.)*ETA*XI*(ETA*ETA-XI*XI)
      ELSE IF((L.EQ.4).AND.(M.EQ.-3)) THEN
         PNSH=0.5*SQRT(0.5*35.)*ZMU*XI*(3.*ETA*ETA-XI*XI)
      ELSE IF((L.EQ.4).AND.(M.EQ.-2)) THEN
         PNSH=SQRT(5.)*(21.*ZMU*ZMU-3.)*ETA*XI/6.
      ELSE IF((L.EQ.4).AND.(M.EQ.-1)) THEN
         PNSH=0.5*SQRT(2.5)*ZMU*XI*(7.*ZMU*ZMU-3.)
      ELSE IF((L.EQ.4).AND.(M.EQ.0)) THEN
         PNSH=(35.*ZMU**4-30.*ZMU*ZMU+3.)/8.
      ELSE IF((L.EQ.4).AND.(M.EQ.1)) THEN
         PNSH=0.5*SQRT(2.5)*ZMU*ETA*(7.*ZMU*ZMU-3.)
      ELSE IF((L.EQ.4).AND.(M.EQ.2)) THEN
         PNSH=SQRT(5.)*(21.*ZMU*ZMU-3.)*(ETA*ETA-XI*XI)/12.
      ELSE IF((L.EQ.4).AND.(M.EQ.3)) THEN
         PNSH=0.5*SQRT(0.5*35.)*ZMU*ETA*(ETA*ETA-3.*XI*XI)
      ELSE IF((L.EQ.4).AND.(M.EQ.4)) THEN
         PNSH=SQRT(35.)*(ETA**4-6.*(ETA*XI)**2+XI**4)/8.
      ELSE IF((L.EQ.5).AND.(M.EQ.-5)) THEN
         PNSH=21.*XI*(5.*ETA**4-10.*(ETA*XI)**2+XI**4)/(8.*SQRT(14.))
      ELSE IF((L.EQ.5).AND.(M.EQ.-4)) THEN
         PNSH=0.5*105.*ZMU*ETA*XI*(ETA*ETA-XI*XI)/SQRT(35.)
      ELSE IF((L.EQ.5).AND.(M.EQ.-3)) THEN
         PNSH=35.*(9*ZMU*ZMU-1.)*XI*(3.*ETA*ETA-XI*XI)/(8.*SQRT(70.))
      ELSE IF((L.EQ.5).AND.(M.EQ.-2)) THEN
         PNSH=0.5*SQRT(105.)*ZMU*(3.*ZMU*ZMU-1.)*ETA*XI
      ELSE IF((L.EQ.5).AND.(M.EQ.-1)) THEN
         PNSH=SQRT(15.)*XI*(21.*ZMU**4-14.*ZMU*ZMU+1.)/8.
      ELSE IF((L.EQ.5).AND.(M.EQ.0)) THEN
         PNSH=ZMU*(63.*ZMU**4-70.*ZMU*ZMU+15.)/8.
      ELSE IF((L.EQ.5).AND.(M.EQ.1)) THEN
         PNSH=SQRT(15.)*ETA*(21.*ZMU**4-14.*ZMU*ZMU+1.)/8.
      ELSE IF((L.EQ.5).AND.(M.EQ.2)) THEN
         PNSH=0.25*SQRT(105.)*ZMU*(3.*ZMU*ZMU-1.)*(ETA*ETA-XI*XI)
      ELSE IF((L.EQ.5).AND.(M.EQ.3)) THEN
         PNSH=35.*(9*ZMU*ZMU-1.)*ETA*(ETA*ETA-3.*XI*XI)/(8.*SQRT(70.))
      ELSE IF((L.EQ.5).AND.(M.EQ.4)) THEN
         PNSH=105.*ZMU*(ETA**4-6.*(ETA*XI)**2+XI**4)/(8.*SQRT(35.))
      ELSE IF((L.EQ.5).AND.(M.EQ.5)) THEN
         PNSH=21.*ETA*(ETA**4-10.*(ETA*XI)**2+5.*XI**4)/(8.*SQRT(14.))
      ELSE
         CALL XABORT('PNSH: LEGENDRE ORDER NOT AVAILABLE.')
      ENDIF
      RETURN
      END
