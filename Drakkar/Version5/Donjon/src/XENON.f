*DECK XENON
      SUBROUTINE XENON(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* computing the Xenon distribution
*
*Copyright:
* Copyright (C) 2010 Ecole Polytechnique de Montreal.
*
*Author(s): M. Guyot
*
*Parameters: input/output
* NENTRY  number of LCM objects or files used by the operator.
* HENTRY  name of each LCM object or file:
*         HENTRY(1): create type(L_MACROLIB);
*         HENTRY(2): modification type(L_MATEX);
*         HENTRY(3): read-only type(L_MACROLIB);
*         HENTRY(4): read-only type(L_MACROLIB) (optional).
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
      PARAMETER(NSTATE=40,IOUT=6)
      CHARACTER HSIGN*12,TEXT*12
      INTEGER ISTATE(NSTATE),ITYP,NITMA
      REAL FLOT
      DOUBLE PRECISION DFLOT
      LOGICAL LINI
      TYPE(C_PTR) IPLIB,IPPOW
      REAL, ALLOCATABLE, DIMENSION(:) :: XEN
*----
*  PARAMETER VALIDATION
*----
      IPLIB=C_NULL_PTR
      IPPOW=C_NULL_PTR
      IF((NENTRY.NE.1).AND.(NENTRY.NE.2)) 
     1  CALL XABORT('@XENON: 1 OR 2 PARAMETERS EXPECTED.')
      DO I=1,NENTRY
         IF((IENTRY(I).NE.1).AND.(IENTRY(I).NE.2))
     1   CALL XABORT('@XENON: LCM OBJECT EXPECTED AT LHS')
      ENDDO
      IF(JENTRY(1).NE.1)CALL XABORT('@XENON: MODIFICATION MODE EXPECTED'
     1 //' FOR L_LIBRARY.')
      IF(NENTRY.EQ.2) THEN
        IF(JENTRY(2).NE.2)CALL XABORT('@XENON: READ-ONLY MODE EXPECTED'
     1   //' FOR L_POWER AT LHS.')
      ENDIF
      DO IEN=1,NENTRY
        CALL LCMGTC(KENTRY(IEN),'SIGNATURE',12,1,HSIGN)
*       L_LIBRARY
        IF(HSIGN.EQ.'L_LIBRARY')THEN
          IPLIB=KENTRY(IEN)
*       L_POWER
        ELSEIF(HSIGN.EQ.'L_POWER')THEN
          IPPOW=KENTRY(IEN)
        ELSE
          TEXT=HENTRY(IEN)
          CALL XABORT('@XENON: SIGNATURE OF '//TEXT//' IS '//HSIGN//
     1    '. L_LIBRARY OR L_POWER EXPECTED.')
        ENDIF
      ENDDO
*----
*  RECOVER INFORMATION
*----
*     L_LIBRARY
      CALL LCMSIX(IPLIB,' ',0)
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(IPLIB,'STATE-VECTOR',ISTATE)
      MAXMIX=ISTATE(1)
      NBISO=ISTATE(2)
      NGRP=ISTATE(3)
      NMIX=ISTATE(14)
*     L_POWER
      IF(C_ASSOCIATED(IPPOW)) THEN
        CALL LCMSIX(IPPOW,' ',0)
        CALL XDISET(ISTATE,NSTATE,0)
        CALL LCMGET(IPPOW,'STATE-VECTOR',ISTATE)
        IF(ISTATE(1).NE.NGRP)CALL XABORT('@XENON: DIFFERENT NGR'
     1   //'P NUMBER IN L_LIBRARY AND L_POWER OBJECT.')
        NCH=ISTATE(6)
        NB=ISTATE(7)
        IF(NCH*NB.NE.NMIX)CALL XABORT('@XENON: DIFFERENT '
     1   //'MIXTURE NUMBER IN L_LIBRARY AND L_POWER OBJECT.')
      ENDIF
*----
*  READ INPUT DATA
*----
      IPRT=0
      LINI=.FALSE.
*     READ KEYWORD
   10 CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(ITYP.NE.3)CALL XABORT('@XENON: CHARACTER DATA EXPECTED(1).')
      IF(TEXT.EQ.'EDIT')THEN
*       PRINTING INDEX
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@XENON: INTEGER DATA EXPECTED.')
        IPRT=MAX(0,NITMA)
        GOTO 10
      ELSEIF(TEXT.EQ.'INIT')THEN
        LINI=.TRUE.
        GOTO 10
      ELSEIF(TEXT.EQ.';')THEN
        GOTO 20
      ELSE
*       KEYWORD DOES NOT MATCH
        CALL XABORT('@XENON: WRONG KEYWORD:'//TEXT//'.')
      ENDIF

   20 IF((.NOT.C_ASSOCIATED(IPPOW)).AND.(.NOT.LINI)) THEN
         CALL XABORT('@XENON: L_POWER OBJECT REQUIRED .')
      ENDIF
      ALLOCATE(XEN(NMIX))
*----
*  COMPUTE THE VALUE OF THE XENON CONCENTRATIONS
*----
      IF(.NOT.LINI) THEN
        CALL XENCAL(IPLIB,IPPOW,NB,NCH,NGRP,NMIX,NBISO,XEN)
      ELSE
        CALL XDRSET(XEN,NMIX,0.0)
      ENDIF
*----
*  PUT THE CONCENTRATIONS IN THE LIBRARY AND COMPUTE NEW XS
*----
      CALL XENLIB(IPLIB,MAXMIX,NMIX,NBISO,NGRP,XEN)
      DEALLOCATE(XEN)
      RETURN
      END
