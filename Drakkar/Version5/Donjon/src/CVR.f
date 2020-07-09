*DECK CVR
      SUBROUTINE CVR(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* perform reordering of fuel-regions properties in the reactor core,
* according to the specified voiding pattern.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): D. Sekki
*
*Parameters: input/output
* NENTRY  number of LCM objects or files used by the operator.
* HENTRY  name of each LCM object or file:
*         HENTRY(1): create type(L_MAP);
*         HENTRY(2): read-only type(L_MAP).
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
      CHARACTER TEXT*12,HSIGN*12
      INTEGER ISTATE(NSTATE),IGST(NSTATE)
      DOUBLE PRECISION DFLOT
      TYPE(C_PTR) IPMAP,JPMAP
*----
*  PARAMETER VALIDATION
*----
      IF(NENTRY.NE.2)CALL XABORT('@CVR: TWO PARAMETERS EXPECTED.')
      TEXT=HENTRY(1)
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2))CALL XABORT('@CVR:'
     1 //' LCM OBJECT EXPECTED AT LHS ('//TEXT//').')
      IF(JENTRY(1).NE.0)CALL XABORT('@CVR: FUEL MAP OBJECT IN CRE'
     1 //'ATE MODE EXPECTED AT LHS ('//TEXT//').')
      TEXT=HENTRY(2)
      IF((IENTRY(2).NE.1).AND.(IENTRY(2).NE.2))CALL XABORT('@CVR:'
     1 //' LCM OBJECT EXPECTED AT RHS ('//TEXT//').')
      IF(JENTRY(2).NE.2)CALL XABORT('@CVR: FUEL MAP OBJECT IN REA'
     1 //'D-ONLY MODE EXPECTED AT RHS ('//TEXT//').')
      CALL LCMGTC(KENTRY(2),'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_MAP')THEN
        TEXT=HENTRY(2)
        CALL XABORT('@CVR: SIGNATURE OF '//TEXT//' IS '//HSIGN//
     1  '. L_MAP EXPECTED.')
      ENDIF
      IPMAP=KENTRY(1)
      CALL LCMEQU(KENTRY(2),IPMAP)
*----
*  RECOVER INFORMATION
*----
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(IPMAP,'STATE-VECTOR',ISTATE)
      NB=ISTATE(1)
      NCH=ISTATE(2)
      NFUEL=ISTATE(7)
      NPARM=ISTATE(8)
*     FUEL-MAP GEOMETRY
      JPMAP=LCMGID(IPMAP,'GEOMAP')
      CALL XDISET(IGST,NSTATE,0)
      CALL LCMGET(JPMAP,'STATE-VECTOR',IGST)
      IF(IGST(1).NE.7)CALL XABORT('@CVR: ONLY 3-D CART'
     1 //'ESIAN GEOMETRY ALLOWED.')
      NX=IGST(3)
      NY=IGST(4)
      NZ=IGST(5)
*     PRINTING INDEX
      CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(ITYP.NE.3)CALL XABORT('@CVR: CHARACTER DATA EXPECTED.')
      IF(TEXT.NE.'EDIT')CALL XABORT('@CVR: KEYWORD EDIT EXPECTED.')
      CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(ITYP.NE.1)CALL XABORT('@CVR: INTEGER FOR EDIT EXPECTED.')
      IMPX=MAX(0,NITMA)
*     READ INPUT DATA
      CALL CVRDRV(IPMAP,NCH,NB,NFUEL,NPARM,NX,NY,NZ,NVD,IVD,IMPX)
*     UPDATE STATE-VECTOR
      ISTATE(10)=NVD
      ISTATE(11)=IVD
      CALL LCMPUT(IPMAP,'STATE-VECTOR',NSTATE,1,ISTATE)
      IF(IMPX.GT.1) CALL LCMLIB(IPMAP)
      RETURN
      END
