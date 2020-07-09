*DECK DEVINI
      SUBROUTINE DEVINI(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* read specification for the rod-devices, create a device object.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): D. Sekki
*
*Parameters: input/output
* NENTRY  number of LCM objects or files used by the operator.
* HENTRY  name of each LCM object or file:
*         HENTRY(1): create type(L_DEVICE);
*         HENTRY(2): modification type(L_MATEX).
* IENTRY  type of each LCM object or file:
*         =1 LCM memory object; =2 XSM file; =3 sequential binary file;
*         =4 sequential ascii file.
* JENTRY  access of each LCM object or file:
*         =0 the LCM object or file is created;
*         =1 the LCM object or file is open for modifications;
*         =2 the LCM object or file is open in read-only mode.
* KENTRY  LCM object address or file unit number.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER      NENTRY,IENTRY(NENTRY),JENTRY(NENTRY)
      TYPE(C_PTR)  KENTRY(NENTRY)
      CHARACTER    HENTRY(NENTRY)*12
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40)
      CHARACTER HSIGN*12,TEXT12*12
      INTEGER ISTATE(NSTATE)
      REAL LIMIT(6)
      TYPE(C_PTR) IPDEV,IPMTX
      REAL, ALLOCATABLE, DIMENSION(:) :: XXX,YYY,ZZZ
*----
*  PARAMETER VALIDATION
*----
      IF(NENTRY.NE.2)CALL XABORT('@DEVINI: TWO PARAMETERS EXPECTED')
      TEXT12=HENTRY(1)
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2))CALL XABORT('@DEVI'
     1 //'NI: LCM OBJECT FOR L_DEVICE EXPECTED ('//TEXT12//').')
      IF(JENTRY(1).NE.0)CALL XABORT('@DEVINI: CREATE MODE EXPECTE'
     1 //'D FOR L_DEVICE.')
      HSIGN='L_DEVICE'
      CALL LCMPTC(KENTRY(1),'SIGNATURE',12,1,HSIGN)
      IPDEV=KENTRY(1)
      TEXT12=HENTRY(2)
      IF((IENTRY(2).NE.1).AND.(IENTRY(2).NE.2))CALL XABORT('@DEVI'
     1 //'NI: LCM OBJECT FOR L_MATEX EXPECTED ('//TEXT12//').')
      CALL LCMGTC(KENTRY(2),'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_MATEX')CALL XABORT('@DEVINI: MISSING L_MATEX.')
      IF(JENTRY(2).NE.1)CALL XABORT('@DEVINI: MODIFICATION MODE E'
     1 //'XPECTED FOR L_MATEX.')
      IPMTX=KENTRY(2)
*----
*  RECOVER INFORMATION
*----
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(IPMTX,'STATE-VECTOR',ISTATE)
      IGEO=ISTATE(6)
      IF(IGEO.NE.7)CALL XABORT('@DEVINI: ON'
     1 //'LY 3D-CARTESIAN GEOMETRY ALLOWED.')
      NMIX=ISTATE(2)
      NTOT=ISTATE(5)
      LX=ISTATE(8)
      LY=ISTATE(9)
      LZ=ISTATE(10)
*     CORE LIMITS ALONG X-AXIS
      ALLOCATE(XXX(LX+1))
      CALL XDRSET(XXX,LX+1,0.)
      CALL LCMGET(IPMTX,'MESHX',XXX)
      LIMIT(1)=XXX(1)
      LIMIT(2)=XXX(LX+1)
      DEALLOCATE(XXX)
*     CORE LIMITS ALONG Y-AXIS
      ALLOCATE(YYY(LY+1))
      CALL XDRSET(YYY,LY+1,0.)
      CALL LCMGET(IPMTX,'MESHY',YYY)
      LIMIT(3)=YYY(1)
      LIMIT(4)=YYY(LY+1)
      DEALLOCATE(YYY)
*     CORE LIMITS ALONG Z-AXIS
      ALLOCATE(ZZZ(LZ+1))
      CALL XDRSET(ZZZ,LZ+1,0.)
      CALL LCMGET(IPMTX,'MESHZ',ZZZ)
      LIMIT(5)=ZZZ(1)
      LIMIT(6)=ZZZ(LZ+1)
      DEALLOCATE(ZZZ)
*     READ ROD-DEVICES DATA
      CALL DEVDRV(IPDEV,IPMTX,IGEO,NMIX,NTOT,LIMIT)
      RETURN
      END
