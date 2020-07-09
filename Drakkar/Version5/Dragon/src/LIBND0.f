*DECK LIBND0
      SUBROUTINE LIBND0 (NAMFIL,NGRO,ENERG_PTR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* recover energy group information from a NDAS library
*
*Copyright:
* Copyright (C) 2006 Atomic Energy of Canada Limited and Ecole
* Polytechnique de Montreal
*
*Author(s): A. Hebert
*
*Parameters: input
* NAMFIL  name of the NDAS file
*
*Parameters: output
* NGRO       number of energy groups
* ENERG_PTR  C_PTR pointer of the energy mesh limit array
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  Subroutine arguments
*----
      INTEGER NGRO
      CHARACTER NAMFIL*(*)
      TYPE(C_PTR) ENERG_PTR
*----
*  Local variables
*----
      INTEGER IERR,HEADER(16)
      REAL, POINTER, DIMENSION(:) :: ENERG
*
      CALL XSDOPN(NAMFIL,IERR)
      IF(IERR.NE.0) CALL XABORT('LIBND0: XSDOPN could not open Library'
     >  //' files')
      CALL XSDBLD(6001,HEADER,IERR)
      IF(IERR.NE.0) CALL XABORT('LIBND0: XSDBLD could not read library'
     > //' parameters')
      NGRO=HEADER(2)
      ENERG_PTR=LCMARA(NGRO+1)
      CALL C_F_POINTER(ENERG_PTR,ENERG,(/ NGRO+1 /))
      CALL XSDBLD(5019,ENERG,IERR)
      IF(IERR.NE.0) CALL XABORT('LIBND0: XSDBLD could not read energy '
     >  //'group limits')
      CALL XSDCL()
      RETURN
      END
