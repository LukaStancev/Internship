*DECK LIBXS3
      SUBROUTINE LIBXS3 (NAMFIL,NGRO,ENERG_PTR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* recover energy group information from an APOLIB-XSM library
*
*Copyright:
* Copyright (C) 2014 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* NAMFIL  name of the APOLIB-XSM file
*
*Parameters: output
* NGRO       number of energy groups
* ENERG_PTR  C_PTR pointer of the energy mesh limit array
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  Subroutine arguments
*----
      INTEGER NGRO
      CHARACTER NAMFIL*(*)
      TYPE(C_PTR) ENERG_PTR
*----
*  Local variables
*----
      TYPE(C_PTR) IPAP
      REAL, POINTER, DIMENSION(:) :: ENERG
*
      CALL LCMOP(IPAP,NAMFIL,2,2,0)
      CALL LCMSIX(IPAP,'PMAIL',1)
      CALL LCMLEN(IPAP,'E',NV,ITYLCM)
      NGRO=NV-1
      ENERG_PTR=LCMARA(NGRO+1)
      CALL C_F_POINTER(ENERG_PTR,ENERG,(/ NGRO+1 /))
      CALL LCMGET(IPAP,'E',ENERG)
      CALL LCMSIX(IPAP,' ',2)
      CALL LCMCL(IPAP,1)
      RETURN
      END
