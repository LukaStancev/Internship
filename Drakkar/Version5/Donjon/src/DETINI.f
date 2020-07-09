*DECK DETINI
      SUBROUTINE DETINI(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
*  reads detector information and stores them 
*
* Copyright (C) 2010 Ecole Polytechnique de Montreal.
*
*Author(s): J. Koclas, E. Varin, M. Guyot
*
*Parameters: input/output
* NENTRY  number of LCM objects or files used by the operator.
* HENTRY  name of each LCM object or file:
*         HENTRY(1): modification type(L_MAP);
*         HENTRY(2): read-only type(L_POWER).
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
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER      NENTRY,IENTRY(NENTRY),JENTRY(NENTRY)
      TYPE(C_PTR)  KENTRY(NENTRY)
      CHARACTER    HENTRY(NENTRY)*12
*----
*  LOCAL VARIABLES
*----
      INTEGER  NSTATE
      PARAMETER (NSTATE=40)
      CHARACTER TEXT*12,HSIGN*12
      INTEGER ISTATE(NSTATE),NGRP,NDETOT,IPRT,IHEX,ITYP,NITMA
      REAL FLOT
      DOUBLE PRECISION DFLOT      
      LOGICAL LHEX,LDET,LENTRY
      TYPE(C_PTR) IPDET
*----
*  PARAMETER VALIDATION
*----
      NDETOT = 0
      NGRP = 0
      LENTRY=.FALSE.
      CALL XDISET(ISTATE,NSTATE,0)
*
      IF(NENTRY.NE.1) CALL XABORT('@DETINI: PARAMETER EXPECTED.')
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2)) CALL XABORT('@D'
     + //'ETINI: LINKED LIST OR XSM FILE EXPECTED AT LHS.')
      IF((JENTRY(1).NE.0).AND.(JENTRY(1).NE.1)) CALL XABORT('@D'
     + //'ETINI: CREATE OR MODIFICATION MODE EXPECTED.')
*
      IPDET=KENTRY(1)
      IF(JENTRY(1).EQ.1) THEN
        TEXT=HENTRY(1)
        LENTRY=.TRUE.
        CALL LCMGTC(KENTRY(1),'SIGNATURE',12,1,HSIGN)
        IF(HSIGN.NE.'L_DETECT')CALL XABORT('@DETINI: L_DETECT'
     +   //' OBJECT IS EXPECTED (OBJECT='//TEXT//')')
        CALL LCMGET(IPDET,'STATE-VECTOR',ISTATE)
        NGRP = ISTATE(1)
        NDETOT = ISTATE(2)
      ENDIF
*----
*  READ INPUT DATA
*----
      IPRT = 0
      LHEX = .FALSE.
      LDET= .FALSE.
      IHEX = 0
   10 CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(ITYP.NE.3) CALL XABORT('@DETINI: CHARACTER DATA'
     +   //' EXPECTED(1).')
      IF(TEXT.EQ.'EDIT') THEN
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@DETINI: INTEGER DATA EXPECTED(1).')
        IPRT=MAX(0,NITMA)
      ELSEIF(TEXT.EQ.'HEXZ')THEN
        LHEX=.TRUE.
      ELSEIF(TEXT.EQ.'NGRP')THEN
        CALL REDGET(ITYP,NGRP,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@DETINI: INTEGER DATA EXPECTED(2).')
        IF(JENTRY(1).EQ.1) THEN
          CALL XABORT('@DETINI: ENERGY GROUP NUMBER REQUIRED ONLY AT'
     +      //' CREATION OF L_DETECT OBJECT')
        ENDIF
      ELSEIF(TEXT.EQ.'TYPE')THEN
        CALL DETDRV(IPDET,NGRP,IPRT,LHEX,NDETOT,LENTRY)
      ELSEIF(TEXT.EQ.';')THEN
        LDET=.TRUE.
      ELSE
        CALL XABORT('@DETINI: INVALID KEYWORD '//TEXT)
      ENDIF
      IF(.NOT.LDET) GOTO 10
*----
*  STATE-VECTOR STORAGE
*----
      IF(JENTRY(1).EQ.0) THEN
        HSIGN='L_DETECT'
        CALL LCMSIX(IPDET,' ',0)
        CALL LCMPTC(IPDET,'SIGNATURE',12,1,HSIGN)
      ENDIF
      CALL XDISET(ISTATE,NSTATE,0)
      ISTATE(1)=NGRP
      ISTATE(2)=NDETOT
      IF(LHEX) ISTATE(3)=1
      CALL LCMPUT(IPDET,'STATE-VECTOR',NSTATE,1,ISTATE)
      IF(IPRT.GT.2) CALL LCMLIB(IPDET)
      RETURN
      END
