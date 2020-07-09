*DECK SNEST
      SUBROUTINE SNEST (IPTRK,NREG,NUN,MAT,KEYFLX,FUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* rearrange SPn flux in the Sn order so that SPn can be used to 
* initialise Sn calculation.
*
*Copyright:
* Copyright (C) 2020 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. A. Calloo
*
*Parameters: input
* IPTRK   pointer to the tracking (L_TRACK signature).
* NREG    total number of regions for which specific values of the
*         neutron flux and reactions rates are required.
* NUN     total number of unknowns in vectors SUNKNO and FUNKNO.
* MAT     index-number of the mixture type assigned to each volume.
* KEYFLX  position of averaged flux elements in FUNKNO vector.
* FUNKNO  SN unknown vector.
*
*Parameters: output
* FUNKNO  SN unknown vector.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK
      INTEGER     NREG,NUN,MAT(NREG),KEYFLX(NREG) 
      REAL        FUNKNO(NUN)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (IUNOUT=6,NSTATE=40)
      INTEGER IPAR(NSTATE),NLOZH,SPLTL,SBMSH,REM,NUNSA,ISPLH
*
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IDL,TMPKEY
      REAL, ALLOCATABLE, DIMENSION(:) :: FUNSA
*
*----
*  RECOVER TRACKING INFORMATION
*----
      CALL LCMSIX(IPTRK,'DSA',1)
      CALL LCMGET(IPTRK,'STATE-VECTOR',IPAR)
      IF(NREG.NE.IPAR(1)) CALL XABORT('SNEST: INVALID NREG ON LCM.')
      NUNSA=IPAR(2)
      ITYPE=IPAR(6)
      ISPLH=1
      IF(ITYPE.EQ.7) THEN
        IELEM=IPAR(9)
        LX=IPAR(14)
        LY=IPAR(15)
        LZ=IPAR(16)
      ELSE
        IELEM=IPAR(8)
        LX=IPAR(12)
        LY=IPAR(13)
        LZ=0
      ENDIF
      IF(ITYPE.EQ.8) ISPLH=IPAR(10)
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IDL(NREG),FUNSA(NUNSA))
*----
*  RECOVER OR BUILD TRACKING INFORMATION
*----
      NLEG=0
      IF(ITYPE.EQ.2) THEN
         NLEG=IELEM
      ELSE IF((ITYPE.EQ.5).OR.(ITYPE.EQ.8)) THEN
         NLEG=IELEM*IELEM
      ELSE IF(ITYPE.EQ.7) THEN
         NLEG=IELEM*IELEM*IELEM
      ELSE
         CALL XABORT('SNEST: TYPE OF DISCRETIZATION NOT IMPLEMENTED.')
      ENDIF
      CALL LCMGET(IPTRK,'KEYFLX',IDL)
*----
*  REBUILD KEYFLX FOR 2D HEXAGONAL CASE
*----
      ! NLOZH - num. of loz. per hexagon
      ! SBMSH - num. of submeshes per lozenge (integer)
      ! SPLTL - split of the lozenge (ISPLH)
      IF(ITYPE.EQ.8)THEN
         ALLOCATE(TMPKEY(NREG))
         IND = 0
         JND = 0
         NLOZH  = 3*ISPLH**2
         NHEX   = NREG/NLOZH 
         SBMSH  = NLOZH/3
         SPLTL  = ISPLH
         DO IH=1,NHEX
            DO IM=1,SBMSH
               REM=MOD(IM-1,SPLTL)
               IF((REM.EQ.0).AND.(SBMSH.NE.1))THEN
                  JND = (IH-1)*NLOZH + SBMSH - (IM/SPLTL)
               ELSEIF((REM.NE.0).AND.(SBMSH.NE.1))THEN
                  JND = JND - (SBMSH*3) - SPLTL
               ENDIF
               DO ILZ=1,3
                  IND = (IH-1)*NLOZH + (IM-1)*3 + (ILZ-1) + 1
                  IF(SBMSH.EQ.1) JND = IND
                  TMPKEY(IND) = KEYFLX(JND)
                  JND = JND + SBMSH
               ENDDO
            ENDDO
         ENDDO
         KEYFLX(:) = TMPKEY(:)
         DEALLOCATE(TMPKEY)
      ELSEIF(ITYPE.EQ.9)THEN
         CALL XABORT('SNEST: 3D HEXAGONAL NOT IMPLEMENTED YET.')
      ENDIF
*----
*  REARRANGE SPN FLUX IN SN ORDER
*----
      CALL XDRSET(FUNSA,NUNSA,0.0)
      FUNSA(1:NUNSA) = FUNKNO(1:NUNSA)
      CALL XDRSET(FUNKNO,NUN,0.0)
      
      CALL LCMSIX(IPTRK,' ',2)

      DO 100 IR=1,NREG
      IF(MAT(IR).LE.0) GO TO 100
      DO  90 IEL=1,NLEG
      IND=KEYFLX(IR)+IEL-1
      JND=IDL(IR)+IEL-1
      FUNKNO(IND)=FUNSA(JND)
   90 CONTINUE
  100 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(FUNSA,IDL)
      RETURN
      END
