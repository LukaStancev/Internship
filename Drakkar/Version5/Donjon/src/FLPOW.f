*DECK FLPOW
      SUBROUTINE FLPOW(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* compute and print power and flux distributions over the reactor core.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): D. Sekki
*
*  Modified 15/07/10 : Creation of L_FLUX object to be used by
*  module DETECT:, M. Guyot
*  
*
*Parameters: input/output
* NENTRY  number of LCM objects or files used by the operator.
* HENTRY  name of each LCM object or file:
*         HENTRY(1): create type(L_POWER);
*         HENTRY(2): optional read-only type(L_POWER);
*         HENTRY(3): read-only type(L_FLUX) or type(L_KINET) ;
*         HENTRY(4): read-only type(L_TRACK);
*         HENTRY(5): optional read-only type(L_MAP);
*         HENTRY(6): optional read-only type(L_MATEX);
*         HENTRY(7): optional read-only type(L_MACROLIB).
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
      CHARACTER HSIGN*12,TEXT*12
      LOGICAL LNEW,LMAP,LFLX,LRAT,LPOW,LFSTH,LFLU,LNRM,LBUN
      DOUBLE PRECISION DFLOT
      TYPE(C_PTR) IPPOW,IPFLX,IPKIN,IPTRK,IPMTX,IPMAP,IPMAC,IPNFX
*----
*  PARAMETER VALIDATION
*----
      LFLU=.FALSE.
      IF(NENTRY.LT.4)CALL XABORT('@FLPOW: PARAMETER EXPECTED.')
      TEXT=HENTRY(1)
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2))CALL XABORT('@FLPOW'
     1 //': LCM OBJECT FOR L_POWER EXPECTED AT LHS ('//TEXT//').')
      IF(JENTRY(1).NE.0)CALL XABORT('@FLPOW: CREATE MODE FOR L_POW'
     1 //'ER EXPECTED AT LHS ('//TEXT//').')
      IPPOW=KENTRY(1)
      IF(JENTRY(2).EQ.0)THEN
        LFLU=.TRUE.
        IPNFX=KENTRY(2)
      ENDIF
      IPFLX=C_NULL_PTR
      IPKIN=C_NULL_PTR
      IPTRK=C_NULL_PTR
      IPMTX=C_NULL_PTR
      IPMAP=C_NULL_PTR
      IPMAC=C_NULL_PTR
      LNEW=.FALSE.
      JMOD=0
      IF(LFLU)THEN
        NRHS=3
      ELSE
        NRHS=2
        IPNFX=C_NULL_PTR
      ENDIF
      DO 10 IEN=NRHS,NENTRY
      IF((IENTRY(IEN).NE.1).AND.(IENTRY(IEN).NE.2))CALL XABORT('@F'
     1 //'LPOW: LCM OBJECT EXPECTED AT THE RHS.')
      CALL LCMGTC(KENTRY(IEN),'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.EQ.'L_POWER')THEN
        IF(JENTRY(IEN).NE.2)CALL XABORT('@FLPOW: READ-ONLY MODE EXPE'
     1  //'CTED FOR THE L_POWER OBJECT AT RHS.')
        IF(LNEW)CALL XABORT('@FLPOW: L_POWER ALREADY DEFINED AT RHS.')
        CALL LCMEQU(KENTRY(IEN),IPPOW)
        LNEW=.TRUE.
      ELSEIF(HSIGN.EQ.'L_MATEX')THEN
        IF(JENTRY(IEN).NE.2)CALL XABORT('@FLPOW: READ-ONLY MODE EXPE'
     1  //'CTED FOR THE L_MATEX OBJECT AT RHS.')
        IF(.NOT.C_ASSOCIATED(IPMTX))THEN
          IPMTX=KENTRY(IEN)
        ELSE
          CALL XABORT('@FLPOW: L_MATEX ALREADY DEFINED.')
        ENDIF
      ELSEIF(HSIGN.EQ.'L_FLUX')THEN
        IF(JENTRY(IEN).NE.2)CALL XABORT('@FLPOW: READ-ONLY MODE EXPE'
     1  //'CTED FOR THE L_FLUX OBJECT AT RHS.')
        IF(.NOT.C_ASSOCIATED(IPFLX))THEN
          IPFLX=KENTRY(IEN)
        ELSE
          CALL XABORT('@FLPOW: L_FLUX ALREADY DEFINED.')
        ENDIF
      ELSEIF(HSIGN.EQ.'L_KINET')THEN
        IF(JENTRY(IEN).NE.2)CALL XABORT('@FLPOW: READ-ONLY MODE EXPE'
     1  //'CTED FOR THE L_KINET OBJECT AT RHS.')
        IF(.NOT.C_ASSOCIATED(IPKIN))THEN
          IPKIN=KENTRY(IEN)
        ELSE
          CALL XABORT('@FLPOW: L_KINET ALREADY DEFINED.')
        ENDIF
      ELSEIF(HSIGN.EQ.'L_TRACK')THEN
        IF(JENTRY(IEN).NE.2)CALL XABORT('@FLPOW: READ-ONLY MODE EXPE'
     1  //'CTED FOR THE L_TRACK OBJECT AT RHS.')
        IF(.NOT.C_ASSOCIATED(IPTRK))THEN
          IPTRK=KENTRY(IEN)
        ELSE
          CALL XABORT('@FLPOW: L_TRACK ALREADY DEFINED.')
        ENDIF
      ELSEIF(HSIGN.EQ.'L_MACROLIB')THEN
        IF(JENTRY(IEN).NE.2)CALL XABORT('@FLPOW: READ-ONLY MODE EXPE'
     1  //'CTED FOR THE L_MACROLIB OBJECT AT RHS.')
        IF(.NOT.C_ASSOCIATED(IPMAC))THEN
          IPMAC=KENTRY(IEN)
        ELSE
          CALL XABORT('@FLPOW: L_MACROLIB ALREADY DEFINED.')
        ENDIF
      ELSEIF(HSIGN.EQ.'L_MAP')THEN
        IF(JENTRY(IEN).EQ.1) JMOD=1
        IF(.NOT.C_ASSOCIATED(IPMAP))THEN
          IPMAP=KENTRY(IEN)
        ELSE
          CALL XABORT('@FLPOW: L_MAP ALREADY DEFINED.')
        ENDIF
      ENDIF
   10 CONTINUE
      IF((.NOT.C_ASSOCIATED(IPFLX)).AND.(.NOT.C_ASSOCIATED(IPKIN))) THEN
         CALL XABORT('@FLPOW: MISSING L_FLUX OR L_KINET OBJECT.')
      ELSE IF((C_ASSOCIATED(IPFLX)).AND.(C_ASSOCIATED(IPKIN))) THEN
         CALL XABORT('@FLPOW: L_FLUX AND L_KINET OBJECTS BOTH DEFINED.')
      ELSE IF(.NOT.C_ASSOCIATED(IPTRK)) THEN
         CALL XABORT('@FLPOW: MISSING L_TRACK OBJECT.')
      ELSE IF((C_ASSOCIATED(IPMAP)).AND.(.NOT.C_ASSOCIATED(IPMTX))) THEN
         CALL XABORT('@FLPOW: MISSING L_MATEX OBJECT.')
      ELSE IF((.NOT.C_ASSOCIATED(IPMAP)).AND.(C_ASSOCIATED(IPMTX))) THEN
         CALL XABORT('@FLPOW: MISSING L_MAP OBJECT.')
      ELSE IF((.NOT.C_ASSOCIATED(IPMTX)).AND.
     1        (.NOT.C_ASSOCIATED(IPMAC))) THEN
         CALL XABORT('@FLPOW: MISSING L_MATEX OR L_MACROLIB OBJECT.')
      ELSE IF((.NOT.C_ASSOCIATED(IPMAP)).AND.
     1        (.NOT.C_ASSOCIATED(IPMAC))) THEN
         CALL XABORT('@FLPOW: MISSING L_MAP OR L_MACROLIB OBJECT.')
      ENDIF
*----
*  READ KEYWORD
*----
      IMPX=1
      PTOT=0.0
      LFSTH=.FALSE.
      LNRM=.FALSE.
      LBUN=.FALSE.
      FSTH=0.0
      LFLX=.FALSE.
      LPOW=.FALSE.
      LMAP=.FALSE.
      LRAT=.FALSE.
   20 CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(ITYP.EQ.10) GO TO 40
   30 IF(ITYP.NE.3)CALL XABORT('@FLPOW: CHARACTER DATA EXPECTED.')
      IF(TEXT.EQ.'EDIT') THEN
*       PRINTING INDEX
        CALL REDGET(ITYP,IMPX,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@FLPOW: INTEGER DATA EXPECTED.')
      ELSE IF(TEXT.EQ.'P-NEW') THEN
        IF(.NOT.LNEW)CALL XABORT('@FLPOW: MISSING READ-ONLY L_POWER'
     1  //' OBJECT AT RHS.')
      ELSE IF(TEXT.EQ.'PTOT') THEN
        IF(LNEW)CALL XABORT('@FLPOW: ONLY ONE L_POWER OBJECT IN CRE'
     1  //'ATE MODE EXPECTED WITH PTOT OPTION.')
        CALL REDGET(ITYP,NITMA,PTOT,TEXT,DFLOT)
        IF(ITYP.NE.2)CALL XABORT('@FLPOW: REAL FOR PTOT EXPECTED.')
        IF(PTOT.LE.0.)CALL XABORT('@FLPOW: INVALID VALUE PTOT < 0.')
      ELSE IF(TEXT.EQ.'FSTH') THEN
        CALL REDGET(ITYP,NITMA,FSTH,TEXT,DFLOT)
        IF(ITYP.NE.2)CALL XABORT('@FLPOW: REAL DATA EXPECTED FOR FSTH.')
        IF((FSTH.GT.1.0).OR.(FSTH.LE.0.0)) CALL XABORT('@FLPOW: FSTH '
     1  //'SHOULD BE BETWEEN 0.0 AND 1.0.')
        LFSTH=.TRUE.
      ELSE IF(TEXT.EQ.'NORM') THEN
        LNRM=.TRUE.
      ELSE IF(TEXT.EQ.'BUND') THEN
        IF(.NOT.C_ASSOCIATED(IPMAP)) CALL XABORT('@FLPOW: NO RHS FUELM'
     1  //'AP DEFINED.')
        LBUN=.TRUE.
      ELSE IF(TEXT.EQ.'PRINT') THEN
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(TEXT.EQ.'MAP')THEN
          IF(.NOT.C_ASSOCIATED(IPMAP))CALL XABORT('@FLPOW: INVALID KEY'
     1    //'WORD MAP. MISSING L_MAP OBJECT FOR PRINT.')
          LMAP=.TRUE.
        ELSEIF(TEXT.EQ.'ALL')THEN
          LFLX=.TRUE.
          LPOW=.TRUE.
          IF(C_ASSOCIATED(IPMAP))LMAP=.TRUE.
          LRAT=.TRUE.
        ELSEIF(TEXT.EQ.'DISTR')THEN
          CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
          IF(ITYP.NE.3)CALL XABORT('@FLPOW: CHARACTER DATA EXPECTED AF'
     1    //'TER DISTR.')
          IF(TEXT.EQ.'FLUX')THEN
            IF(LFLX)CALL XABORT('@FLPOW: KEYWORD FLUX ALREADY READ.')
            LFLX=.TRUE.
          ELSEIF(TEXT.EQ.'POWER')THEN
            IF(LPOW)CALL XABORT('@FLPOW: KEYWORD POWER ALREADY READ.')
            LPOW=.TRUE.
          ELSEIF(TEXT.EQ.'RATIO')THEN
            IF(LRAT)CALL XABORT('@FLPOW: KEYWORD RATIO ALREADY READ.')
            LRAT=.TRUE.
          ELSE
            GO TO 30
          ENDIF
        ELSE
          CALL XABORT('@FLPOW: KEYWORD MAP/DISTR/ALL EXPECTED.')
        ENDIF
      ELSE IF(TEXT.EQ.'INIT') THEN
        IF(JENTRY(IEN).EQ.1) JMOD=2
      ELSE IF(TEXT.EQ.';') THEN
        GO TO 40
      ELSE
        CALL XABORT('@FLPOW: INVALID KEYWORD '//TEXT//'.')
      ENDIF
      GO TO 20
*----
*  CHECK CONSISTENCY
*----
   40 IF(LMAP) THEN
        IF(.NOT.C_ASSOCIATED(IPMAP)) THEN
           CALL XABORT('@FLPOW: MISSING L_MAP OBJECT.')
        ELSE IF(.NOT.C_ASSOCIATED(IPMTX)) THEN
           CALL XABORT('@FLPOW: MISSING L_MATEX OBJECT.')
        ELSE IF(.NOT.C_ASSOCIATED(IPMAC)) THEN
           CALL XABORT('@FLPOW: MISSING L_MACROLIB OBJECT.')
        ENDIF
      ENDIF
*----
*  PERFORM CALCULATION
*----
      CALL FLPDRV(IPPOW,IPNFX,IPFLX,IPKIN,IPTRK,IPMTX,IPMAP,IPMAC,PTOT,
     1 LNEW,LMAP,JMOD,LFLX,LPOW,LRAT,IMPX,FSTH,LFSTH,LFLU,LBUN,LNRM)
      RETURN
      END