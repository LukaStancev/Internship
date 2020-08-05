*DECK AFMDRV
      SUBROUTINE AFMDRV (KENTRY,NENTRY,NPARM,ITYPE,NBURN,NGRP,NISO,ISC,
     1 MNPS,NL,ILEAK,NTYP,NBCH,NCCO,NCZO,NUT,CTITRE,LMCR,IXYZ,MMIX,MSFT,
     2 NISM)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Driver to generate a macrolib using fbm
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): 
* M.T. Sissaoui
*Update(s):
*  E. Varin 28/03/00, B. Dionne 26/02/01, 
*  A. Lagarrigue 30/07/05
*  A. Hebert 11/11/11 (remove table support)
*
*Parameters: input
* KENTRY  address of the LCM objects
* NENTRY  number of LCM objects
* NPARM   number of parameters in L_MAP object
* ITYPE   creation/modification flag for output macrolib
* NBURN   number of burnup steps
* NGRP    1+number of energy groups
* NISO    number of extracted isotopes
* ISC     type of cross-section calculation (=1: time average;
*         =2: instantaneous; =3: homogeneous)
* MNPS    number of shifts + 2
* NL      number of legendre orders (=1 for isotropic scattering)
* ILEAK   type of leakage
* NTYP
* NBCH    number of bundles per channel
* NCCO    number of channels in the core
* NCZO    number of combustion zones
* NUT     number of fuel types
* CTITRE  character*72 title
* LMCR    if true, create a macrolib containing only one non-zero
*         mixture
* IXYZ    type of diffusion coefficient (=0: isotropic; =1: directional)
* MMIX    number of mixtures in the output macrolib
* MSFT    second dimension of BSFT and PSFT
* NISM
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) KENTRY(NENTRY)
      INTEGER NPARM,ITYPE,NBURN,NGRP,NISO,ISC,MNPS,NL,ILEAK,NTYP,NBCH,
     1 NCCO,NCZO,NUT,IXYZ,MMIX,MSFT,NISM
      CHARACTER*72 CTITRE
      LOGICAL LMCR
*----
*  LOCAL VARIABLES
*----
      CHARACTER TEXTR*12,CM*2,TEXT4*5,HMICRO*12,TEXTB*12,TEXTD*12
      TYPE(C_PTR) IPMACX,JPMAC,KPMAC,IPFBM,IPMAP,JPMAP,KPMAP
      DOUBLE PRECISION DFLOTT,XCOF(3)
      REAL  STORE,RLOC(7)
      LOGICAL LNOMP,LTAV,LXENON,LSAM,LNEP,LXEREF,LNEREF,LTFUEL,LDRAH,
     1 LTCOOL,LDCOOL,LPWF,LINI
      CHARACTER PNAME*12,PARKEY*12
      INTEGER, DIMENSION(:), ALLOCATABLE :: IPOS,IJ,IZONE,IWORK,NJ,
     1 HISO,JTAB,INDEX,KTYP,ISFT,ITEXTR
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: IJJ,NJJ
      REAL, DIMENSION(:), ALLOCATABLE :: VOL,ENER,WORK,BURBG,BURED,
     1 POWER,PW,BRH,XSIGF,XSIGX,XFLUN,PDCOOL,PTCOOL,PTFUEL,SSCAT
      REAL, DIMENSION(:,:), ALLOCATABLE :: XBURN,OVERV,SIGS,FLUX,CHI,
     1 DIFFX,DIFFY,DIFFZ,FLUAV,BFLUX,BSFT,PSFT
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: SIGMA,SIGAV,DENSITB,HXEN1,
     1 HXEN2,HSAM1,HSAM2,HNEP1,HNEP2,CPW1B,CPW2B,FLUXB,CHIB,OVERVB
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: SCAT,SCATAV
      REAL, DIMENSION(:,:,:,:,:), ALLOCATABLE :: SMACB,XBORB,XXENB,
     1 XT1FB,XT2FB,XT1CB,XT2CB,XT1MB,XT2MB,XD1CB,XD2CB,XD1MB,XD2MB,
     2 XSMB,XNP9B,XMFDB,XMMDB,XPF1B,XPF2B,XPF1LB,XPF2LB,XPURB
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(SIGMA(MMIX,NGRP,NTYP),IJJ(MMIX,NL,NGRP),VOL(MMIX),
     1 NJJ(MMIX,NL,NGRP),XBURN(NBURN,NUT),OVERV(MMIX,NGRP),
     2 SIGS(MMIX,NGRP),FLUX(MMIX,NGRP),CHI(MMIX,NGRP),ENER(NGRP+1),
     3 IPOS(MMIX),SCAT(MMIX,NL,NGRP,NGRP),DIFFX(MMIX,NGRP),
     4 DIFFY(MMIX,NGRP),DIFFZ(MMIX,NGRP),IJ(NGRP),WORK(MMIX*NGRP*NBURN),
     5 IZONE(NCCO),BURBG(MMIX),BURED(MMIX),POWER(MMIX),
     6 FLUAV(NBURN,NGRP),SIGAV(NBURN,NGRP,NTYP),IWORK(MMIX*NGRP),
     7 SCATAV(NBURN,NL,NGRP,NGRP),PW(MNPS),BRH(MNPS),NJ(NGRP),
     8 BFLUX(NGRP,MMIX),DENSITB(NISO,NBURN,NUT),HISO(3*NISM),
     9 HXEN1(2,NBURN,NUT),HXEN2(2,NBURN,NUT),HSAM1(2,NBURN,NUT),
     1 HSAM2(2,NBURN,NUT),HNEP1(2,NBURN,NUT),HNEP2(2,NBURN,NUT),
     2 CPW1B(2,NBURN,NUT),CPW2B(2,NBURN,NUT),FLUXB(NGRP,NBURN,NUT),
     3 JTAB(NISO),CHIB(NGRP,NBURN,NUT),OVERVB(NGRP,NBURN,NUT),
     4 INDEX(MMIX),KTYP(NUT),XSIGF(NGRP),XSIGX(NGRP),XFLUN(NGRP),
     5 BSFT(MMIX,MSFT),PSFT(MMIX,MSFT),ISFT(MMIX),PDCOOL(MMIX),
     6 PTCOOL(MMIX),PTFUEL(MMIX),ITEXTR(3*NUT))
       ALLOCATE(SMACB(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     1          XBORB(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     2          XXENB(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     3          XT1FB(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     4          XT2FB(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     5          XT1CB(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     6          XT2CB(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     7          XT1MB(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     8          XT2MB(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     9          XD1CB(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     1          XD2CB(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     2          XD1MB(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     3          XD2MB(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     4          XSMB(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     5          XNP9B(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     6          XMFDB(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     7          XMMDB(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     8          XPF1B(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     9          XPF2B(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     1          XPF1LB(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     2          XPF2LB(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     3          XPURB(NGRP*NGRP,NTYP,NISO,NBURN,NUT))
*
      IPMACX=KENTRY(1)
      IPFBM=KENTRY(2)
      IF( .NOT.LMCR )IPMAP=KENTRY(3)
      CALL LCMLIB(IPMAP)
*---------------------------------------------------------------*
*      SET THE DEFAULT OPTIONS
      LNOMP=.FALSE.
      LTAV=.FALSE.
      LXENON=.FALSE.
      LSAM=.FALSE.
      LNEP=.FALSE.
      LXEREF=.FALSE.
      LNEREF=.FALSE.
      LTFUEL=.FALSE.
      LDRAH =.FALSE.
      LTCOOL=.FALSE.
      LDCOOL=.FALSE.
      LPWF=.TRUE.
      ILBFLU=0
      IMPX=0
      IXENO=0
      ISAMA=0
      INEPT=0
      IPROF2=0
      LINI=.FALSE.
      ILEAK=0
      PWREF=0.0
      DMR=0.0
      DCR=0.0
      NTM=0
* Set burnup interpolation method 
*   (default 0 for lagrangian interpolation)
*   (1 for linear)
      ILIN=0 
*     SET  HERMITE INTERPOLATION FOR TIME-AVERAGE CALCULATION
      ITM=3
*---------------------------------------------------------------*
*      MX IS THE MAXIMUN MIXTURE NUMBER
      MX=NBCH*NCCO
*---------------------------------------------------------------*
*     CHECK THE PARAMETERS
      IF(MX.EQ.0) CALL XABORT('AFMDRV: ZERO NUMBER OF MIXTURES.')
      IF(NGRP.EQ.0) CALL XABORT('AFMDRV: ZERO NUMBER OF GROUPS.')
      IF(NBURN.EQ.0) CALL XABORT('AFMDRV: ZERO NUMBER OF BURNUPS.')
*---------------------------------------------------------------*
*     INITIALISATION OF THE MATRICES
      NG2=NGRP*NGRP
      DO 81 IGR=1,NG2
       DO 81 IN=1,NUT
        DO 81 I=1,NBURN
         DO 81 ITY=1,NTYP
          DO 81 ISO=1,NISO
            XBORB(IGR,ITY,ISO,I,IN)=0.0
            XPURB(IGR,ITY,ISO,I,IN)=0.0
            XXENB(IGR,ITY,ISO,I,IN)=0.0
            XT1FB(IGR,ITY,ISO,I,IN)=0.0
            XT2FB(IGR,ITY,ISO,I,IN)=0.0
            XT1CB(IGR,ITY,ISO,I,IN)=0.0
            XT2CB(IGR,ITY,ISO,I,IN)=0.0
            XT1MB(IGR,ITY,ISO,I,IN)=0.0
            XT2MB(IGR,ITY,ISO,I,IN)=0.0
            XD1CB(IGR,ITY,ISO,I,IN)=0.0
            XD2CB(IGR,ITY,ISO,I,IN)=0.0
            XD1MB(IGR,ITY,ISO,I,IN)=0.0
            XD2MB(IGR,ITY,ISO,I,IN)=0.0
            XSMB(IGR,ITY,ISO,I,IN)=0.0
            XNP9B(IGR,ITY,ISO,I,IN)=0.0
            XMFDB(IGR,ITY,ISO,I,IN)=0.0
            XMMDB(IGR,ITY,ISO,I,IN)=0.0
            XPF1B(IGR,ITY,ISO,I,IN)=0.0
            XPF2B(IGR,ITY,ISO,I,IN)=0.0
            XPF1LB(IGR,ITY,ISO,I,IN)=0.0
            XPF2LB(IGR,ITY,ISO,I,IN)=0.0
            SMACB(IGR,ITY,ISO,I,IN)=0.0
   81 CONTINUE
*
      DO 10 IGR=1,NGRP
        DO 10 IMX=1,MX
           DIFFX(IMX,IGR)=0.0
           DIFFY(IMX,IGR)=0.0
           DIFFZ(IMX,IGR)=0.0
           FLUX(IMX,IGR)=0.0
           OVERV(IMX,IGR)=0.0
           CHI(IMX,IGR)=0.0
           DO 16 IL=1,NL
             DO 17 JGR=1,NGRP
               SCAT(IMX,IL,IGR,JGR)=0.0
   17        CONTINUE
             IJJ(IMX,IL,IGR)=IGR
             NJJ(IMX,IL,IGR)=1
   16      CONTINUE
           DO 11 ITYP=1,NTYP
             SIGMA(IMX,IGR,ITYP)=0.0
   11      CONTINUE
   10 CONTINUE
C
      DO 34 IBR=1,NBURN
        DO 34 IGR=1,NGRP
          FLUAV(IBR,IGR)=0.0
          DO 31 ITYP=1,NTYP
           SIGAV(IBR,IGR,ITYP)=0.0
   31     CONTINUE
          DO 32 JGR=1,NGRP
            DO 32 IL=1,NL
              SCATAV(IBR,IL,IGR,JGR)=0.0
   32     CONTINUE
   34 CONTINUE
* INITIALISATION OF THE HISTORY COEFFICIENT
      DO 82 IBR=1,NBURN
       DO 82 IN=1,NUT
         DO 82 I=1,2
           CPW1B(I,IBR,IN)=0.0
           CPW2B(I,IBR,IN)=0.0
           HXEN1(I,IBR,IN)=0.0
           HXEN2(I,IBR,IN)=0.0
           HSAM1(I,IBR,IN)=0.0
           HSAM2(I,IBR,IN)=0.0
           HNEP1(I,IBR,IN)=0.0
           HNEP2(I,IBR,IN)=0.0
   82 CONTINUE
*---------------------------------------------------------------*
* READ AN OPTION KEY WORD
  800 CALL REDGET (INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
      IF(INDIC.NE.3) CALL XABORT('AFMDRV: CHARACTER DATA EXPECTED.')
      IF(TEXT4.EQ.'EDIT') THEN
* READ THE PRINT INDEX.
         CALL REDGET(INDIC,IMPX,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('AFMDRV: INTEGER DATA EXPECTED.')
      ELSE IF(TEXT4.EQ.'REFT') THEN
         DO 111 IN=1,NUT
           CALL REDGET(INDIC,KTYP(IN),FLOTT,TEXT4,DFLOTT)
           IF(INDIC.NE.1) CALL XABORT('AFMDRV: INTEGER DATA EXPECTED.')
           CALL REDGET (INDIC,NITMA,FLOTT,TEXTR,DFLOTT)
           IF(INDIC.NE.3)
     1         CALL XABORT('AFMDRV: CHARACTER DATA EXPECTED.')
           READ(TEXTR,'(3A4)') (ITEXTR((IN-1)*3+I),I=1,3)
 111     CONTINUE
         IF(LMCR .AND. KTYP(1).GT.MX)
     +            CALL XABORT('AFMDRV: INVALID INDEX NUMBER.')
C
* CHECK THE NAME OF THE DIRECTORY
         WRITE(TEXTR,'(3A4)') (ITEXTR(I1),I1=1,3)
         CALL LCMLEN(IPFBM,TEXTR,ILENGT,ITYLCM)
         IF(ILENGT.EQ.0) THEN
           CALL XABORT('AFMDRV: UNABLE TO FIND '//TEXTR//' .')
         ENDIF
* RECOVER THE REFERENCE LOCAL PARAMETERS VALUES
         CALL LCMSIX(IPFBM,TEXTR,1)
         CALL LCMSIX(IPFBM,'INFO-NOMINA',1)
         CALL LCMLEN(IPFBM,'NOMINALP',ILP,ITYLCM)
         IF(ILP.GT.0) THEN
           CALL LCMGET(IPFBM,'NOMINALP',RLOC)
           CALL LCMGET(IPFBM,'NOMINALN',HISO)
           DO 888 I=1,ILP
             WRITE(HMICRO,'(3A4)') (HISO((I-1)*3+IH),IH=1,3)
             IF(HMICRO.EQ.'PW') PWREF=RLOC(I)
             IF(HMICRO.EQ.'TCOOL') TCR=RLOC(I)
             IF(HMICRO.EQ.'TMOD') TMR=RLOC(I)
             IF(HMICRO.EQ.'TFUEL') TFR=RLOC(I)
             IF(HMICRO.EQ.'RHOC') DCR=RLOC(I)
             IF(HMICRO.EQ.'RHOM') DMR=RLOC(I)
             IF(HMICRO.EQ.'PUR') XIR=RLOC(I)
 888       CONTINUE
         ENDIF
         CALL LCMSIX(IPFBM,' ',2)
         CALL LCMSIX(IPFBM,' ',2)
* REFERENCE PARAMETER VALUES
         PFIX=PWREF
         AW=15.9994 +2*(1-XIR)*1.0079 +2*XIR*2.014101
         PH=2*1.0079/AW
         PD=2*2.014101/AW
* INITIALISATION OF PERTURBED PARAMETER
         TF=TFR
         TC=TCR
         TM=TMR
         DC=1.0
         DM=1.0
         XI=XIR
         BOR=0.0
         SM=0.0
         RNP9=0.0
         XEN=0.0
*
         DO 15 IMX=1,MX
           POWER(IMX)=PWREF
           ISFT(IMX)=0
           BURBG(IMX)=0.0
           BURED(IMX)=0.0
           VOL(IMX)=0.0
           PDCOOL(IMX)=DCR
           PTCOOL(IMX)=TCR
           PTFUEL(IMX)=TFR
   15    CONTINUE
*        RECOVER THE TEMERATURE AND DENSITY PROFILES
         IF( (.NOT.LMCR).AND.(NPARM.GT.0) ) THEN
            JPMAP=LCMGID(IPMAP,'PARAM')
            DO 18 IPARM=1,NPARM
            KPMAP=LCMGIL(JPMAP,IPARM)
            CALL LCMGTC(KPMAP,'P-NAME',12,1,PNAME)
            CALL LCMGTC(KPMAP,'PARKEY',12,1,PARKEY)
            CALL LCMGET(KPMAP,'P-TYPE',IPTYPE)
            IF(IPTYPE.EQ.1) THEN
               CALL LCMGET(KPMAP,'P-VALUE',FLOTT)
            ELSE IF(IPTYPE.EQ.2) THEN
               CALL LCMLEN(KPMAP,'P-VALUE',NITMA,ITYLCM)
               IF(NITMA.NE.MX) CALL XABORT('@AFMDRV: INVALID LENGTH FO'
     1         //'R P-VALUE.')
            ENDIF
            IF(PNAME.EQ.'T-COOL') THEN
               WRITE(6,716) PNAME,PARKEY
               IF(IPTYPE.EQ.1) THEN
                  CALL XDRSET(PTCOOL,MX,FLOTT)
               ELSE IF(IPTYPE.EQ.2) THEN
                  CALL LCMGET(KPMAP,'P-VALUE',PTCOOL)
               ENDIF
            ELSE IF(PNAME.EQ.'D-COOL') THEN
               WRITE(6,716) PNAME,PARKEY
               IF(IPTYPE.EQ.1) THEN
                  CALL XDRSET(PDCOOL,MX,FLOTT)
               ELSE IF(IPTYPE.EQ.2) THEN
                  CALL LCMGET(KPMAP,'P-VALUE',PDCOOL)
               ENDIF
            ELSE IF(PNAME.EQ.'T-FUEL') THEN
               WRITE(6,716) PNAME,PARKEY
               IF(IPTYPE.EQ.1) THEN
                  CALL XDRSET(PTFUEL,MX,FLOTT)
               ELSE IF(IPTYPE.EQ.2) THEN
                  CALL LCMGET(KPMAP,'P-VALUE',PTFUEL)
               ENDIF
            ENDIF
   18       CONTINUE
         ENDIF
*
         CALL XDRSET(PW,MNPS,PWREF)
         CALL XDRSET(BRH,MNPS,0.0)
         CALL XDRSET(POWER,MX,PWREF)
*
      ELSE IF(TEXT4.EQ.'TFUEL') THEN
        CALL REDGET (INDIC,NITMA,TFU,TEXT4,DFLOTT)
        LTFUEL = .TRUE.
        IF(INDIC.NE.2) CALL XABORT('AFMDRV: REAL DATA EXPECTED.')
*
      ELSE IF(TEXT4.EQ.'TCOOL') THEN
         CALL REDGET (INDIC,NITMA,TCU,TEXT4,DFLOTT)
         IF(INDIC.NE.2) CALL XABORT('AFMDRV: REAL DATA EXPECTED.')
          LTCOOL = .TRUE.
          CALL XDRSET(PTCOOL,MX,TCU)
*
      ELSE IF(TEXT4.EQ.'TMOD') THEN
        CALL REDGET (INDIC,NITMA,TM,TEXT4,DFLOTT)
        IF(INDIC.NE.2) CALL XABORT('AFMDRV: REAL DATA EXPECTED.')
*
      ELSE IF(TEXT4.EQ.'RDCL') THEN
        CALL REDGET (INDIC,NITMA,DCU,TEXT4,DFLOTT)
        IF(INDIC.NE.2) CALL XABORT('AFMDRV: REAL DATA EXPECTED.')
          LDCOOL = .TRUE.
          CALL XDRSET(PDCOOL,MX,DCU)
*
      ELSE IF(TEXT4.EQ.'RDMD') THEN
        CALL REDGET (INDIC,NITMA,DM,TEXT4,DFLOTT)
        IF(INDIC.NE.2) CALL XABORT('AFMDRV: REAL DATA EXPECTED.')
        DM=DM/DMR
*
      ELSE IF(TEXT4.EQ.'BORON') THEN
        CALL REDGET (INDIC,NITMA,BOR,TEXT4,DFLOTT)
        IF(INDIC.NE.2) CALL XABORT('AFMDRV: REAL DATA EXPECTED.')
*
*  ppm eq 10**-6, NO CONSISTENCY WITH CFC CONCENTRATIONS
*  NEED TO ADD A COEFFICIENT TO FIT THE DATA (BREF should be 0.0ppm)
*
        BOR=BOR*1.E-6
*
      ELSE IF(TEXT4.EQ.'PUR') THEN
        CALL REDGET (INDIC,NITMA,XI,TEXT4,DFLOTT)
        IF(INDIC.NE.2) CALL XABORT('AFMDRV: REAL DATA EXPECTED.')
        XI=XI*1.0E-02
*
      ELSE IF(TEXT4.EQ.'FIXP') THEN
        CALL REDGET (INDIC,NITMA,PFIX,TEXT4,DFLOTT)
        IF(INDIC.EQ.2) THEN
          LNOMP=.TRUE.
        ELSE IF(TEXT4.EQ.'INIT') THEN
          LINI=.TRUE.
        ELSE
          CALL XABORT('AFMDRV: "INIT" or REAL DATA EXPECTED.')
        ENDIF
*
      ELSE IF(TEXT4.EQ.'IMET') THEN
        CALL REDGET(INDIC,ITM,FLOTT,TEXT4,DFLOTT)
        IF(INDIC.NE.1) CALL XABORT('AFMDRV: INTEGER DATA EXPECTED.')
*
      ELSE IF(TEXT4.EQ.'XENON') THEN
        LXENON=.TRUE.
        CALL REDGET (INDIC,NITMA,FXEN,TEXT4,DFLOTT)
        IF(INDIC.NE.2) CALL XABORT('AFMDRV: REAL DATA EXPECTED.')
*
      ELSE IF(TEXT4.EQ.'XEREF') THEN
        LXEREF=.TRUE.
*
      ELSE IF(TEXT4.EQ.'DRAH') THEN
         LDRAH=.TRUE.
*
      ELSE IF(TEXT4.EQ.'SAM') THEN
        LSAM=.TRUE.
        CALL REDGET (INDIC,NITMA,FSAM,TEXT4,DFLOTT)
        IF(INDIC.NE.2) CALL XABORT('AFMDRV: REAL DATA EXPECTED.')
*
      ELSE IF(TEXT4.EQ.'NEP') THEN
        LNEP=.TRUE.
        CALL REDGET (INDIC,NITMA,FNEP,TEXT4,DFLOTT)
        IF(INDIC.NE.2) CALL XABORT('AFMDRV: REAL DATA EXPECTED.')
*
      ELSE IF(TEXT4.EQ.'NREF') THEN
        LNEREF=.TRUE.
*
      ELSE IF(TEXT4.EQ.'BURN') THEN
        IF(LMCR) THEN
          CALL REDGET (INDIC,NITMA,FBUR,TEXT4,DFLOTT)
          IF(INDIC.NE.2) CALL XABORT('AFMDRV: REAL DATA EXPECTED.')
        ELSE
          CALL XABORT('AFMDRV: INVALID KEYWORD BURN.')
        ENDIF
*
      ELSE IF(TEXT4.EQ.'NPWF') THEN
         LPWF=.FALSE.
      ELSE IF(TEXT4.EQ.'PWF') THEN
         LPWF=.TRUE.
      ELSE IF(TEXT4.EQ.'BLIN') THEN
         ILIN=1
      ELSE IF(TEXT4.EQ.';') THEN
         GO TO 160
      ELSE
        CALL XABORT('AFMDRV: '//TEXT4//' IS AN INVALID KEY-WORD.')
      ENDIF
      GO TO 800
*     EQUIVALENT MODERATOR DENSITY FOR THE REFERENCE PURITY
  160 DXI = XI - XIR
* pas de modification de densite selon la purete D2O
*      DM=DM/(1.0+DXI*(PD-PH))
*---------------------------------------------------------------*
* RECOVER NEUTRONICS PARAMETRES
      WRITE(TEXTR,'(3A4)') (ITEXTR(I1),I1=1,3)
      CALL LCMSIX(IPFBM,TEXTR,1)
      CALL LCMGET(IPFBM,'VOLUME',VOL(1))
      CALL LCMGET(IPFBM,'ENERGY',ENER)
      CALL LCMGET(IPFBM,'HITAB',HISO)
      CALL LCMGET(IPFBM,'JTAB',JTAB)
      CALL LCMSIX(IPFBM,' ',2)
      DO 112 IN=1,NUT
        WRITE(TEXTR,'(3A4)') (ITEXTR((IN-1)*3+I1),I1=1,3)
        CALL LCMSIX(IPFBM,TEXTR,1)
        CALL LCMGET(IPFBM,'BURNUP',XBURN(1,IN))
*     RECOVER THE EXISTING DATABASE.
* RECOVER THE HISTORY COEFFICIENTS
        DO 910 I = 1,NBURN
          WRITE(TEXTB,'(4HBURN,4X,I4)') I
          CALL LCMSIX(IPFBM,TEXTB,1)
*
          IF(JTAB(1).EQ.1) THEN
            CALL LCMSIX(IPFBM,'HISTORY',1)
            CALL LCMGET(IPFBM,'PHIL1',CPW1B(1,I,IN))
            CALL LCMGET(IPFBM,'PHIS1',CPW1B(2,I,IN))
            CALL LCMGET(IPFBM,'PHIL2',CPW2B(1,I,IN))
            CALL LCMGET(IPFBM,'PHIS2',CPW2B(2,I,IN))
            CALL LCMLEN(IPFBM,'PHISX1',IHISTO,ITYLCM)
            IF(IHISTO.GT.0) THEN
              CALL LCMGET(IPFBM,'PHILX1',HXEN1(1,I,IN))
              CALL LCMGET(IPFBM,'PHISX1',HXEN1(2,I,IN))
              CALL LCMGET(IPFBM,'PHILX2',HXEN2(1,I,IN))
              CALL LCMGET(IPFBM,'PHISX2',HXEN2(2,I,IN))
C
              CALL LCMGET(IPFBM,'PHILS1',HSAM1(1,I,IN))
              CALL LCMGET(IPFBM,'PHISS1',HSAM1(2,I,IN))
              CALL LCMGET(IPFBM,'PHILS2',HSAM2(1,I,IN))
              CALL LCMGET(IPFBM,'PHISS2',HSAM2(2,I,IN))
C
              CALL LCMGET(IPFBM,'PHILN1',HNEP1(1,I,IN))
              CALL LCMGET(IPFBM,'PHISN1',HNEP1(2,I,IN))
              CALL LCMGET(IPFBM,'PHILN2',HNEP2(1,I,IN))
              CALL LCMGET(IPFBM,'PHISN2',HNEP2(2,I,IN))
            ENDIF
            CALL LCMSIX(IPFBM,' ',2)
          ENDIF
*
          CALL LCMGET(IPFBM,'FLUX-INTG',FLUXB(1,I,IN))
          CALL LCMGET(IPFBM,'OVERV',OVERVB(1,I,IN))
          CALL LCMGET(IPFBM,'ISOTOPESDENS',DENSITB(1,I,IN))
* COMPUTE DELTA-CONCENTRATION
          DO 49 ISO=1,NISO
            WRITE(HMICRO,'(3A4)') (HISO((ISO-1)*3+IH),IH=1,3)
            CALL LCMSIX(IPFBM,HMICRO,1)
            IF(JTAB(1).EQ.1) THEN
              IF((HMICRO.EQ.'XE135').OR.(HMICRO.EQ.'Xe135')) IXENO=ISO
              IF((HMICRO.EQ.'SM149').OR.(HMICRO.EQ.'Sm149')) ISAMA=ISO
              IF((HMICRO.EQ.'NP239').OR.(HMICRO.EQ.'Np239')) INEPT=ISO
              IF(HMICRO.EQ.'MACR ')
     1             CALL LCMGET(IPFBM,'CHI',CHIB(1,I,IN))
            ENDIF
* RECOVER MACROSCOPIC X-SECTIONS
            NTM=4+2*IXYZ
            DO 100 ITY=1,NTM
              IF(ITY.EQ.1) THEN
                IF(IXYZ.EQ.0)  THEN
                  TEXTD = 'STRD'
                ELSE IF(IXYZ.EQ.1) THEN
                  TEXTD = 'STRD X'
                ENDIF
              ENDIF
              IF(ITY.EQ.2) TEXTD = 'ABS'
              IF(ITY.EQ.3) TEXTD = 'NUSIGF'
              IF(ITY.EQ.4) TEXTD = 'H-FACTORS'
              IF(ITY.EQ.5) TEXTD = 'STRD Y'
              IF(ITY.EQ.6) TEXTD = 'STRD Z'
              CALL LCMLEN(IPFBM,TEXTD,ILENG,ITYXSM)
*
              IF(ILENG.NE.0) THEN
                CALL LCMSIX(IPFBM,TEXTD,1)
                CALL LCMGET(IPFBM,'REF',SMACB(1,ITY,ISO,I,IN))
                CALL LCMGET(IPFBM,'BOR',XBORB(1,ITY,ISO,I,IN))
                CALL LCMGET(IPFBM,'PUR',XPURB(1,ITY,ISO,I,IN))
                CALL LCMGET(IPFBM,'T1M',XT1MB(1,ITY,ISO,I,IN))
                CALL LCMGET(IPFBM,'T2M',XT2MB(1,ITY,ISO,I,IN))
                CALL LCMGET(IPFBM,'D1M',XD1MB(1,ITY,ISO,I,IN))
                CALL LCMGET(IPFBM,'D2M',XD2MB(1,ITY,ISO,I,IN))
                IF(JTAB(1).EQ.1) THEN
                  CALL LCMLEN(IPFBM,'XEN',ILENGX,ITYXSM)
                  IF(ILENGX.GT.0)
     +               CALL LCMGET(IPFBM,'XEN',XXENB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'T1F',XT1FB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'T2F',XT2FB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'T1C',XT1CB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'T2C',XT2CB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'D1C',XD1CB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'D2C',XD2CB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'SM149',XSMB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'NP239',XNP9B(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'MIXFD',XMFDB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'MIXMD',XMMDB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'FPCH1',XPF1B(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'FPCL1',XPF1LB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'FPCH2',XPF2B(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'FPCL2',XPF2LB(1,ITY,ISO,I,IN))
                ENDIF
*
                CALL LCMSIX(IPFBM,' ',2)
              ENDIF
 100        CONTINUE
*
            CALL LCMLEN(IPFBM,'NFTOT',ILNF,ITYXSM)
            IF(ILNF.NE.0) THEN
              CALL LCMGET(IPFBM,'NFTOT',SMACB(1,NTM+1,ISO,I,IN))
            ENDIF
            CALL LCMSIX(IPFBM,' ',2)
  49      CONTINUE
*     SCATTERING CROSS-SECTIONS
          DO 161 ISO=1,NISO
            WRITE(HMICRO,'(3A4)') (HISO((ISO-1)*3+IH),IH=1,3)
            CALL LCMLEN(IPFBM,HMICRO,ILENG,ITYLCM)
            IF(ILENG.EQ.0) GO TO 160
            CALL LCMSIX(IPFBM,HMICRO,1)
C            DO 150 IL=1,NL
              IL=1
              ITY=NTM+1+IL
              LTST=0
              WRITE (CM,'(I2.2)') IL-1
              CALL LCMLEN(IPFBM,'SCAT'//CM,ILENG,ITYXSM)
              IF(ILENG.NE.0) THEN
                LTST=1
              ELSE
                WRITE (CM,'(I2)') IL-1
                CALL LCMLEN(IPFBM,'SCAT'//CM,ILENG,ITYXSM)
                IF(ILENG.NE.0) THEN
                  LTST=2
                ENDIF
              ENDIF
              IF (LTST.GE.1) THEN
                CALL LCMSIX(IPFBM,'SCAT'//CM,1)
                IF(HMICRO.EQ.'MACR') THEN
                  IF (LTST.EQ.1) THEN
                    CALL LCMGET(IPFBM,'NJJS',NJ)
                    CALL LCMGET(IPFBM,'IJJS',IJ)
                  ELSE
                    CALL LCMGET(IPFBM,'NJJ',NJ)
                    CALL LCMGET(IPFBM,'IJJ',IJ)
                  ENDIF
                ENDIF
                CALL LCMGET(IPFBM,'REF',SMACB(1,ITY,ISO,I,IN))
                CALL LCMGET(IPFBM,'BOR',XBORB(1,ITY,ISO,I,IN))
                CALL LCMGET(IPFBM,'PUR',XPURB(1,ITY,ISO,I,IN))
                CALL LCMGET(IPFBM,'T1M',XT1MB(1,ITY,ISO,I,IN))
                CALL LCMGET(IPFBM,'T2M',XT2MB(1,ITY,ISO,I,IN))
                CALL LCMGET(IPFBM,'D1M',XD1MB(1,ITY,ISO,I,IN))
                CALL LCMGET(IPFBM,'D2M',XD2MB(1,ITY,ISO,I,IN))
                IF(JTAB(1).EQ.1) THEN
                  CALL LCMLEN(IPFBM,'XEN',ILENG,ITYXSM)
                  IF(ILENG.GT.0) THEN
                    CALL LCMGET(IPFBM,'XEN',XXENB(1,ITY,ISO,I,IN))
                  ENDIF
                  CALL LCMGET(IPFBM,'T1F',XT1FB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'T2F',XT2FB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'T1C',XT1CB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'T2C',XT2CB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'D1C',XD1CB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'D2C',XD2CB(1,ITY,ISO,I,IN))
*
                  CALL LCMGET(IPFBM,'SM149',XSMB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'NP239',XNP9B(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'MIXFD',XMFDB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'MIXMD',XMMDB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'FPCH1',XPF1B(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'FPCL1',XPF1LB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'FPCH2',XPF2B(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'FPCL2',XPF2LB(1,ITY,ISO,I,IN))
                ENDIF
*
                CALL LCMSIX(IPFBM,' ',2)
              ENDIF
C  150       CONTINUE
            CALL LCMSIX(IPFBM,' ',2)
  161     CONTINUE
C
          CALL LCMSIX(IPFBM,' ',2)
  910   CONTINUE
        CALL LCMSIX(IPFBM,' ',2)
 112  CONTINUE
      IF(JTAB(1).EQ.1) THEN
        IF(IXENO.EQ.0) CALL XABORT('NO XE135 FOUND ')
        IF(ISAMA.EQ.0) CALL XABORT('NO SM149 FOUND ')
        IF(INEPT.EQ.0) CALL XABORT('NO NP239 FOUND ')
      ENDIF
*     END OF THE RECOVERING PROCESS
*---------------------------------------------------------------*
*        ISC INDICATE THE TYPE OF CROSS-SECTION CALCULATION
*        ISC=1 ; TIME AVERAGE CALCULATION
*        ISC=2 ; INSTANTANEOUS CALCULATION
*        ISC=3 ; HOMOGENEOUS CALCULATION
*---------------------------------------------------------------*
*
      IF(ISC.EQ.0) THEN
        CALL XABORT('AFMDRV: TIMAV/INSTANT BURNUP TREATMENT NOT SET')
      ELSE IF(ISC.EQ.1) THEN
*       Time-averaged calculation
        WRITE(6,700)
        MMIX=NBCH*NCCO
        LTAV=.TRUE.
        CALL LCMGET(IPMAP,'FLMIX',INDEX)
        CALL LCMGET(IPMAP,'BURN-BEG',BURBG)
        CALL LCMGET(IPMAP,'BURN-END',BURED)
        CALL LCMLEN(IPMAP,'BUND-PW',ILPW,ITYLCM)
          IF((ILPW.NE.0).AND.LPWF) THEN
            IF(IMPX.GE.1) WRITE(6,702)
            IF(.NOT.LINI) THEN
              CALL LCMGET(IPMAP,'BUND-PW',POWER)
            ELSE 
              CALL LCMLEN(IPMAP,'BUND-PW-INI',ILPW,ITYLCM)
              IF(ILPW.NE.0) THEN
                CALL LCMGET(IPMAP,'BUND-PW-INI',POWER)
              ELSE
                CALL XABORT('AFMDRV: NO INITIAL POWER IN L_MAP')
              ENDIF
            ENDIF
          ELSE
            CALL XDRSET(POWER,MMIX,PWREF)
          ENDIF
        CALL LCMLEN(IPMAP,'FLUX-AV',ILBFLU,ITYLCM)
        IF(ILBFLU.NE.0) THEN
          IF(IMPX.GE.1) WRITE(6,703)
          CALL LCMGET(IPMAP,'FLUX-AV',WORK)
          DO 13 IGR=1,NGRP
           DO 13 IBF=1,MMIX
             IIBF=MMIX*(IGR-1)+IBF
             BFLUX(IGR,IBF)=WORK(IIBF)
 13       CONTINUE
        ENDIF
      ELSE IF(ISC.EQ.2) THEN
*       Instantaneous calculation
        IF(LMCR) THEN
          MMIX=NBCH*NCCO
          CALL XDRSET(POWER,MMIX,PWREF)
        ELSE
          WRITE(6,701)
          MMIX=NBCH*NCCO
          CALL LCMGET(IPMAP,'FLMIX',INDEX)
          CALL LCMGET(IPMAP,'BURN-INST',BURBG)
          CALL LCMLEN(IPMAP,'BUND-PW',ILPW,ITYLCM)
          IF((ILPW.NE.0).AND.LPWF) THEN
            IF(IMPX.GE.1) WRITE(6,702)
            IF(.NOT.LINI) THEN
              CALL LCMGET(IPMAP,'BUND-PW',POWER)
            ELSE 
              CALL LCMLEN(IPMAP,'BUND-PW-INI',ILPW,ITYLCM)
              IF(ILPW.NE.0) THEN
                CALL LCMGET(IPMAP,'BUND-PW-INI',POWER)
              ELSE
                CALL XABORT('AFMDRV: NO INITIAL POWER IN L_MAP')
              ENDIF
            ENDIF
          ELSE
            CALL XDRSET(POWER,MMIX,PWREF)
          ENDIF
          CALL LCMLEN(IPMAP,'FLUX-AV',ILBFLU,ITYLCM)
          IF(ILBFLU.NE.0) THEN
            IF(IMPX.GE.1) WRITE(6,703)
            CALL LCMGET(IPMAP,'FLUX-AV',WORK)
            DO 12 IGR=1,NGRP
             DO 12 IBF=1,MMIX
               IIBF=MMIX*(IGR-1)+IBF
               BFLUX(IGR,IBF)=WORK(IIBF)
 12         CONTINUE
          ENDIF
* RECOVER THE SHIFT INFORMATION
          IF(MNPS.GT.2) THEN
            IF(IMPX.GE.1) WRITE(6,704)
            CALL LCMGET(IPMAP,'ISHIFT',ISFT)
            DO 96 IS=1,MNPS-2
              WRITE (CM,'(I2)') IS
              CALL LCMGET(IPMAP,'BSHIFT'//CM,BSFT(1,IS))
              CALL LCMGET(IPMAP,'PSHIFT'//CM,PSFT(1,IS))
 96         CONTINUE
          ENDIF
        ENDIF
      ELSE IF(ISC.EQ.3) THEN
*       Homogeneous calculation
        MMIX=NCZO
        LTAV=.TRUE.
        CALL LCMGET(IPMAP,'B-ZONE',IZONE)
        CALL LCMGET(IPMAP,'FLMIX',INDEX)
        CALL LCMGET(IPMAP,'BURN-AVG',BURED)
      ENDIF
*---------------------------------------------------------------*
      IF(IMPX.GE.1) THEN
         IF(LNOMP) WRITE(6,705) PFIX
         IF(LXENON) WRITE(6,706) FXEN
         IF(LSAM) WRITE(6,710) FSAM
         IF(LNEP) WRITE(6,711) FNEP
         IF(LXEREF) WRITE(6,712)
         IF(LNEREF) WRITE(6,713)
         IF(LTFUEL) WRITE(6,714) TFU
         IF(IHISTO.GT.0.AND.LDRAH) WRITE(6,715)
         IF(LTCOOL) WRITE(6,717) TCU
         IF(LDCOOL) WRITE(6,718) DCU
      ENDIF
*---------------------------------------------------------------*
* MIXTURE SHIFT
      IF(LMCR) THEN
        MXSH=MMIX
        CALL XDRSET(VOL,MMIX,VOL(1))
      ELSE
        MXSH=1
      ENDIF
*---------------------------------------------------------------*
*        LOOP OVER THE MIXTURES
      DO 302 NMIX=MXSH,MMIX
        TC=PTCOOL(NMIX)
        DC=PDCOOL(NMIX)/DCR
        IF(LMCR) THEN
          NPS=2
          IDF=1
        ELSE
          VOL(NMIX)=VOL(1)
          NPS=ISFT(NMIX)+2
          KDF=0
          DO 113 IN=1,NUT
            IF(INDEX(NMIX).EQ.KTYP(IN)) THEN
              IDF=IN
              KDF=1
            ENDIF
  113     CONTINUE
          IF(KDF.EQ.0) CALL XABORT('AFMDRV: WRONG NUMBER OF INDEX')
        ENDIF
* IF TIME AVERAGE CALCULATION:
* EVALUATION OF THE BURNUPS STEPS EMBEDED IN THE INTEGRATION
        IF(LTAV) THEN
          XBMIN=BURBG(NMIX)
          XBMAX=BURED(NMIX)
* TIME AVERAGE BURNUP LOCALISATION
          CALL AFMLOC(NBURN,NTP,XBMAX,XBMIN,XBURN(1,IDF),
     1                IMAX,IMIN,XCOF,ILIN)
*       LAGRANGE METHOD (TIME-AVERAGE)
          IMINR=IMIN
          IMAXR=ABS(IMAX)
*       SPLINE OR HERMITE METHOD (TIME-AVERAGE)
          IF(ITM.EQ.2.OR.ITM.EQ.3) THEN
            IMINR=1
            IMAXR=NBURN
          ENDIF
*
        ELSE
          IMINR=1
          IMAXR=1
        ENDIF
C
        DO 175 JR=IMINR,IMAXR
         IF(LTAV) THEN
           IRAV=JR
           NPS=2
         ELSE
           IF(NPS.GT.2) THEN
             DO 75 K=2,NPS-1
               IS=K-1
               BRH(K)=BSFT(NMIX,IS)
 75          CONTINUE
           ENDIF
           IF(LMCR) THEN
             BRH(NPS)=FBUR
             IF(JTAB(1).EQ.0) BRH(NPS)=0.0
           ELSE
             BRH(NPS)=BURBG(NMIX)
           ENDIF
         ENDIF
*
         IF(LNOMP) THEN
           DO 76 K=2,NPS
             PW(K)=PFIX
 76        CONTINUE
         ELSE
           IF(NPS.GT.2) THEN
             DO 95 K=2,NPS-1
               IS=K-1
               PW(K)=PSFT(NMIX,IS)
 95          CONTINUE
           ENDIF
           PW(NPS)=POWER(NMIX)
         ENDIF
*        D. Rozon 'Introduction a la Cinetique des Reacteur Nucleaires'
*        Edition E.P., 1992. (p.217) or 1998 (p.185)
*        PW is assumed to be in kW.
         IF(IPROF2.GT.0) THEN
           TF = PTFUEL(NMIX)
         ELSE
           TF= TC + 0.476*PW(NPS) + 2.267*PW(NPS)*PW(NPS)*1.0E-04
         ENDIF
C INITIAL CONCENTRATIONS
         ZXREF=0.0
         SM=0.0
         ZRNP9=0.0
*      IF FUEL
         IF(JTAB(1).EQ.1) THEN
* BURNUP LOCALISATION FOR XENON AND FISSION X-SECTION INTERPOLATION
           IF(LTAV) THEN
             XIFL=XBURN(IRAV,IDF)
             IMAXX=IRAV
             IMINX=IRAV
             XCOF(1)=1.0D0
             XCOF(2)=0.0D0
             XCOF(3)=0.0D0
           ELSE
             XIFL=BRH(NPS)
             CALL AFMLOC(NBURN,NTP,BRH(NPS),BRH(NPS),XBURN(1,IDF),
     1                 IMAXX,IMINX,XCOF,ILIN)
           ENDIF
*
           DO 92 IGR = 1,NGRP
             XSIGX(IGR)=0.0
             XFLUN(IGR)=0.0
             XSIGF(IGR)=0.0
  92       CONTINUE
*       INTERPOLATION OF THE CONCENTRATION
*
           IIX=0
           DO 90 I = IMINX,IMAXX
             IIX=IIX+1
             RXCOF=REAL(XCOF(IIX))
             ZXREF=DENSITB(IXENO,I,IDF)*RXCOF  +ZXREF
             XEN=ZXREF
             SM=DENSITB(ISAMA,I,IDF)*RXCOF     +SM
             ZRNP9=DENSITB(INEPT,I,IDF)*RXCOF  +ZRNP9
             RNP9=ZRNP9
*
             DO 90 IGR=1,NGRP
               XSIGX(IGR)=SMACB(IGR,2,IXENO,I,IDF)*RXCOF
     1                    + XSIGX(IGR)
               XFLUN(IGR)=FLUXB(IGR,I,IDF)*RXCOF     + XFLUN(IGR)
               XSIGF(IGR)=SMACB(IGR,5,1,I,IDF)*RXCOF + XSIGF(IGR)
  90       CONTINUE
           IF(LDRAH.AND.IHISTO.GT.0) THEN
             IF(PW(NPS).GT.PWREF) THEN
               XPW=ALOG(PW(NPS)/PW(1))
               XPWM=1.0/PW(NPS)-1.0/PW(1)
               IFH=1
             ELSE
               XPW=PW(NPS)-PW(1)
               XPWM=(PW(NPS)-PW(1))**2
               IFH=2
             ENDIF
C
             XEN  =ZXREF
             RNP9 =ZRNP9
             IIX=0
             DO 290 I = IMINX,IMAXX
               IIX=IIX+1
               RXCOF=REAL(XCOF(IIX))
*            COMPUTE XENON-SAMRIUM-NEPTUNIUM CONCENTRATION USING DRAGON
               XEN  =XEN   +HXEN1(IFH,I,IDF)*XPW*RXCOF+
     1                      HXEN2(IFH,I,IDF)*XPWM*RXCOF
               SM  =SM     +HSAM1(IFH,I,IDF)*XPW*RXCOF+
     1                      HSAM2(IFH,I,IDF)*XPWM*RXCOF
               RNP9 =RNP9  +HNEP1(IFH,I,IDF)*XPW*RXCOF+
     1                      HNEP2(IFH,I,IDF)*XPWM*RXCOF
 290         CONTINUE
           ELSE IF(ILBFLU.NE.0.AND.XIFL.NE.0.0) THEN
*          COMPUTE THE XENON AND NEPTUNIUM CONCENTRATIONS
              CALL AFMXNC(NGRP,XSIGX,XSIGF,BFLUX(1,NMIX),
     1                   XEN,RNP9,XFLUN)
           ENDIF
* COMPUTE THE XENON AND NEPTUNIUM CONCENTRATIONS
           IF(LXENON) XEN=FXEN
           IF(LSAM) SM=FSAM
           IF(LNEP) RNP9=FNEP
           IF(LXEREF) XEN=ZXREF
           IF(LNEREF) RNP9=ZRNP9
           IF(LTFUEL) THEN
!       fuel temperature as input
             TF=TFU
!       reference fuel temperature 
           ELSEIF(LMCR) THEN
             TF=TFR
           ENDIF
         ENDIF
*---------------------------------------------------------------*
* XSECTION CALCULATION
*---------------------------------------------------------------*
         CALL AFMCPT(KENTRY,NBURN,NGRP,NISO,
     1   NL,IMPX,SMACB,XBORB,XPURB,XXENB,XT1FB,XT2FB,XT1CB,
     1   XT2CB,XT1MB,XT2MB,XD1CB,XD2CB,XD1MB,XD2MB,
     1   XSMB,XNP9B,XMFDB,XMMDB,XPF1B,XPF2B,XPF1LB,XPF2LB,
     1   DENSITB,CPW1B,CPW2B,FLUXB,OVERVB,CHIB,
     1   IJ,NJ,HISO,CTITRE,
     1   NMIX,SIGMA,NTYP,TF,TC,TM,DC,DM,BOR,XEN,SM,RNP9,XI,
     1   TFR,TCR,TMR,XIR,OVERV,FLUX,CHI,SCAT,MX,NPS,PW,BRH,
     1   XBURN,LTAV,IRAV,IDF,JTAB,IXYZ,ILIN)
*---------------------------------------------------------------*
*
         DO 102 IGR=1,NGRP
          FLUAV(JR,IGR)=FLUX(NMIX,IGR)
          DO 102 ITY=1,NTM+1
            SIGAV(JR,IGR,ITY)=SIGMA(NMIX,IGR,ITY)
 102     CONTINUE
         IL =1
         DO 103 IGR=1,NGRP
            DO 103 JGR=1,NGRP
              SCATAV(JR,IL,JGR,IGR)=SCAT(NMIX,IL,JGR,IGR)
 103     CONTINUE
 175    CONTINUE
        IF(LTAV) THEN
* COMPUTE  TIME AVERAGED X-SECTIONS
          DO 101 IGR=1,NGRP
            CALL AFMTAV(NBURN,ITM,XBMAX,XBMIN,FLUAV(1,IGR),IMIN,IMAX,
     1      XBURN,FLUX(NMIX,IGR))
            DO 101 ITY=1,NTM+1
              CALL AFMTAV(NBURN,ITM,XBMAX,XBMIN,SIGAV(1,IGR,ITY),IMIN,
     1        IMAX,XBURN,SIGMA(NMIX,IGR,ITY))
 101      CONTINUE
*
          DO 104 IGR=1,NGRP
           DO 104 JGR=1,NGRP
               IL=1
               CALL AFMTAV(NBURN,ITM,XBMAX,XBMIN,SCATAV(1,IL,IGR,JGR),
     1         IMIN,IMAX,XBURN,SCAT(NMIX,IL,IGR,JGR))
 104      CONTINUE
*
        ENDIF
* COMPUTE DIRECTIONAL DIFFUSION COEFFICIENTS FROM STRD
*  X-SECTIONS.
        IF(IXYZ.EQ.0) THEN
          DO 170 IGR=1,NGRP
            DIFFX(NMIX,IGR)=1.0/(3.0*SIGMA(NMIX,IGR,1))
 170      CONTINUE
          ILEAK=1
        ELSE IF(IXYZ.EQ.1) THEN
          DO 171 IGR=1,NGRP
            DIFFX(NMIX,IGR)=1.0/(3.0*SIGMA(NMIX,IGR,1))
            DIFFY(NMIX,IGR)=1.0/(3.0*SIGMA(NMIX,IGR,5))
            DIFFZ(NMIX,IGR)=1.0/(3.0*SIGMA(NMIX,IGR,6))
            ILEAK=2
 171      CONTINUE
        ENDIF
*
        IL=1
        DO 980 IGR=1,NGRP
            NJJ(NMIX,IL,IGR)=NJ(IGR)
            IJJ(NMIX,IL,IGR)=IJ(IGR)
            IF(LMCR) THEN
              DO 98 NI=1,MMIX
                NJJ(NI,IL,IGR)=NJ(IGR)
                IJJ(NI,IL,IGR)=IJ(IGR)
  98          CONTINUE
            ENDIF
 980    CONTINUE
* MIX LOOP
  302 CONTINUE
*
      IF(LTAV) THEN
        IF(IMPX.GE.1.AND.ITM.EQ.1) WRITE(6,707)
        IF(IMPX.GE.1.AND.ITM.EQ.2) WRITE(6,708)
        IF(IMPX.GE.1.AND.ITM.EQ.3) WRITE(6,709)
      ENDIF
*---------------------------------------------------------------*
*        DECOMPRESS BURN ZONE FOR ALL THE BUNDLES
      IF(ISC.EQ.3) THEN
       MMIX=NBCH*NCCO
       DO 303 IGR=1,NGRP
        DO 304 IZ=1,NCZO
            WORK(IZ)=DIFFX(IZ,IGR)
 304    CONTINUE
        DO 305 IC=1,NCCO
          DO 305 IB=1,NBCH
            ICB=NBCH*(IC-1)+IB
            DIFFX(ICB,IGR)=WORK(IZONE(IC))
 305    CONTINUE
*
        IF(ILEAK.EQ.2) THEN
          DO 314 IZ=1,NCZO
            WORK(IZ)=DIFFY(IZ,IGR)
 314      CONTINUE
          DO 315 IC=1,NCCO
           DO 315 IB=1,NBCH
             ICB=NBCH*(IC-1)+IB
             DIFFY(ICB,IGR)=WORK(IZONE(IC))
 315      CONTINUE
*
          DO 324 IZ=1,NCZO
            WORK(IZ)=DIFFZ(IZ,IGR)
 324      CONTINUE
          DO 325 IC=1,NCCO
           DO 325 IB=1,NBCH
             ICB=NBCH*(IC-1)+IB
             DIFFZ(ICB,IGR)=WORK(IZONE(IC))
 325      CONTINUE
        ENDIF
*
        DO 336 ITY=2,NTM+1
         DO 334 IZ=1,NCZO
           WORK(IZ)=SIGMA(IZ,IGR,ITY)
 334     CONTINUE
         DO 335 IC=1,NCCO
          DO 335 IB=1,NBCH
            ICB=NBCH*(IC-1)+IB
            SIGMA(ICB,IGR,ITY)=WORK(IZONE(IC))
 335     CONTINUE
 336    CONTINUE
*
        DO 354 IZ=1,NCZO
          WORK(IZ)=FLUX(IZ,IGR)
 354    CONTINUE
        DO 355 IC=1,NCCO
         DO 355 IB=1,NBCH
           ICB=NBCH*(IC-1)+IB
           FLUX(ICB,IGR)=WORK(IZONE(IC))
 355    CONTINUE
*
        DO 364 IZ=1,NCZO
          WORK(IZ)=OVERV(IZ,IGR)
 364    CONTINUE
        DO 365 IC=1,NCCO
         DO 365 IB=1,NBCH
           ICB=NBCH*(IC-1)+IB
           OVERV(ICB,IGR)=WORK(IZONE(IC))
 365    CONTINUE
*
        DO 374 IZ=1,NCZO
          WORK(IZ)=CHI(IZ,IGR)
 374    CONTINUE
        DO 375 IC=1,NCCO
         DO 375 IB=1,NBCH
           ICB=NBCH*(IC-1)+IB
           CHI(ICB,IGR)=WORK(IZONE(IC))
 375    CONTINUE
*
        IL=1
        DO 377 JGR=1,NGRP
           DO 378 IZ=1,NCZO
             WORK(IZ)=SCAT(IZ,IL,IGR,JGR)
 378       CONTINUE
           DO 379 IC=1,NCCO
            DO 379 IB=1,NBCH
              ICB=NBCH*(IC-1)+IB
              SCAT(ICB,IL,IGR,JGR)=WORK(IZONE(IC))
 379       CONTINUE
 377     CONTINUE
*
         DO 384 IZ=1,NCZO
           IWORK(IZ)=NJJ(IZ,IL,IGR)
 384     CONTINUE
         DO 385 IC=1,NCCO
          DO 385 IB=1,NBCH
            ICB=NBCH*(IC-1)+IB
            NJJ(ICB,IL,IGR)=IWORK(IZONE(IC))
 385     CONTINUE
*
         DO 394 IZ=1,NCZO
           IWORK(IZ)=IJJ(IZ,IL,IGR)
 394     CONTINUE
         DO 395 IC=1,NCCO
          DO 395 IB=1,NBCH
            ICB=NBCH*(IC-1)+IB
            IJJ(ICB,IL,IGR)=IWORK(IZONE(IC))
 395     CONTINUE
*
 303   CONTINUE
*
       DO 344 IZ=1,NCZO
         WORK(IZ)=VOL(IZ)
 344   CONTINUE
       DO 345 IC=1,NCCO
        DO 345 IB=1,NBCH
          ICB=NBCH*(IC-1)+IB
          VOL(ICB)=WORK(IZONE(IC))
 345   CONTINUE
*
      ENDIF
*---
* STORE MACROLIB INFORMATIONS
*---
      IF(ITYPE.EQ.0)THEN
        CALL LCMPUT(IPMACX,'VOLUME',MMIX,2,VOL)
        CALL LCMPUT(IPMACX,'ENERGY',NGRP+1,2,ENER)
      ENDIF
*
      IF(LMCR) THEN
        STORE=VOL(MMIX)
        VOL(MMIX)= 0.0
*  MACROLIB EN MODIFICATION
        IF(ITYPE.NE.0) THEN
           CALL LCMGET(IPMACX,'VOLUME',VOL)
        ENDIF
        VOL(KTYP(1)) = STORE
        CALL LCMPUT(IPMACX,'VOLUME',MMIX,2,VOL)
        JPMAC=LCMLID(IPMACX,'GROUP',NGRP)
        DO 211 JGR=1,NGRP
          KPMAC=LCMDIL(JPMAC,JGR)
          STORE=SIGMA(MMIX,JGR,2)
          SIGMA(MMIX,JGR,2) = 0.0
*  MACROLIB EN MODIFICATION
          IF(ITYPE.NE.0) THEN
            CALL LCMGET(KPMAC,'NTOT0',SIGMA(1,JGR,2))
          ENDIF
          SIGMA(KTYP(1),JGR,2) = STORE
*
          STORE=OVERV(MMIX,JGR)
          OVERV(MMIX,JGR) = 0.0
*  MACROLIB EN MODIFICATION
          IF(ITYPE.NE.0) THEN
            CALL LCMGET(KPMAC,'OVERV',OVERV(1,JGR))
          ENDIF
          OVERV(KTYP(1),JGR) = STORE
*
          STORE=DIFFX(MMIX,JGR)
          DIFFX(MMIX,JGR) = 0.0
*  MACROLIB EN MODIFICATION
          IF(ITYPE.NE.0) THEN
             CALL LCMGET(KPMAC,'DIFFX',DIFFX(1,JGR))
          ENDIF
          DIFFX(KTYP(1),JGR) = STORE
*
          IF(ILEAK.EQ.2) THEN
            STORE=DIFFY(MMIX,JGR)
            DIFFY(MMIX,JGR) = 0.0
            IF(ITYPE.NE.0) THEN
               CALL LCMGET(KPMAC,'DIFFY',DIFFY(1,JGR))
            ENDIF
            DIFFY(KTYP(1),JGR) = STORE
*
            STORE=DIFFZ(MMIX,JGR)
            DIFFZ(MMIX,JGR) = 0.0
            IF(ITYPE.NE.0) THEN
               CALL LCMGET(KPMAC,'DIFFZ',DIFFZ(1,JGR))
            ENDIF
            DIFFZ(KTYP(1),JGR) = STORE
          ENDIF
*
          STORE = FLUX(MMIX,JGR)
          FLUX(MMIX,JGR) = 0.0
          IF(ITYPE.NE.0) THEN
             CALL LCMGET(KPMAC,'FLUX-INTG',FLUX(1,JGR))
          ENDIF
          FLUX(KTYP(1),JGR) = STORE
*
          IF(JTAB(1).EQ.1 .OR. ITYPE.NE.0) THEN
            STORE = CHI(MMIX,JGR)
            CHI(MMIX,JGR) = 0.0
            IF(ITYPE.NE.0) THEN
              CALL LCMGET(KPMAC,'CHI',CHI(1,JGR))
            ENDIF
            CHI(KTYP(1),JGR) = STORE
*
            STORE=SIGMA(MMIX,JGR,3)
            SIGMA(MMIX,JGR,3) = 0.0
            IF(ITYPE.NE.0) THEN
              CALL LCMGET(KPMAC,'NUSIGF',SIGMA(1,JGR,3))
            ENDIF
            SIGMA(KTYP(1),JGR,3) = STORE
*
            STORE=SIGMA(MMIX,JGR,5)
            SIGMA(MMIX,JGR,5) = 0.0
            IF(ITYPE.NE.0) THEN
              CALL LCMGET(KPMAC,'NFTOT',SIGMA(1,JGR,5))
            ENDIF
            SIGMA(KTYP(1),JGR,5) = STORE
*
            STORE=SIGMA(MMIX,JGR,4)
            SIGMA(MMIX,JGR,4) = 0.0
            IF(ITYPE.NE.0) THEN
              CALL LCMGET(KPMAC,'H-FACTOR',SIGMA(1,JGR,4))
            ENDIF
            SIGMA(KTYP(1),JGR,4) = STORE
*
          ENDIF
*
          IL=1
          ALLOCATE(SSCAT(NGRP))
          DO 212 IGR=1,NGRP
            SSCAT(IGR)= SCAT(MMIX,IL,IGR,JGR)
            SCAT(MMIX,IL,IGR,JGR) = 0.0
 212      CONTINUE
          IF(ITYPE.NE.0) THEN
!!  ATTENTION isotropy is supposed
!!
            IL=1
            WRITE (CM,'(I2.2)') IL-1
            CALL LCMGET(KPMAC,'SCAT'//CM,WORK)
            CALL LCMGET(KPMAC,'NJJS'//CM,NJJ(1,IL,JGR))
            CALL LCMGET(KPMAC,'IJJS'//CM,IJJ(1,IL,JGR))
            CALL LCMGET(KPMAC,'IPOS'//CM,IPOS)
            DO 213 IBM=1,MMIX
               IJJ0=IJJ(IBM,IL,JGR)
               IPOSDE = IPOS(IBM)
               DO 213 IGR=IJJ0,IJJ0-NJJ(IBM,IL,JGR)+1,-1
                SCAT(IBM,IL,IGR,JGR)=WORK(IPOSDE)
                IPOSDE=IPOSDE+1
 213        CONTINUE
          ENDIF
*
          DO 214 IGR=1,NGRP
            SCAT(KTYP(1),IL,IGR,JGR) = SSCAT(IGR)
 214      CONTINUE
          DEALLOCATE(SSCAT)
  211   CONTINUE
      ENDIF
*
      DO 131 IX=1,MMIX
       DO 130 JGR=1,NGRP
        DO 130 IL=1,NL
          IGMIN=JGR
          IGMAX=JGR
          DO 120 IGR=NGRP,1,-1
           IF (SCAT(IX,IL,IGR,JGR).NE.0.0) THEN
             IGMIN=MIN(IGMIN,IGR)
             IGMAX=MAX(IGMAX,IGR)
           ENDIF
  120     CONTINUE
          IJJ(IX,IL,JGR)=IGMAX
          NJJ(IX,IL,JGR)=IGMAX-IGMIN+1
  130  CONTINUE
  131 CONTINUE
*
      CALL XDRSET(SIGS,MMIX*NGRP,0.0)
      JPMAC=LCMLID(IPMACX,'GROUP',NGRP)
      DO 210 JGR=1,NGRP
        KPMAC=LCMDIL(JPMAC,JGR)
        CALL LCMPUT(KPMAC,'NTOT0',MMIX,2,SIGMA(1,JGR,2))
        CALL LCMPUT(KPMAC,'OVERV',MMIX,2,OVERV(1,JGR))
        IF(ILEAK.EQ.1) THEN
          CALL LCMPUT(KPMAC,'DIFF',MMIX,2,DIFFX(1,JGR))
        ELSE IF(ILEAK.EQ.2) THEN
          CALL LCMPUT(KPMAC,'DIFFX',MMIX,2,DIFFX(1,JGR))
          CALL LCMPUT(KPMAC,'DIFFY',MMIX,2,DIFFY(1,JGR))
          CALL LCMPUT(KPMAC,'DIFFZ',MMIX,2,DIFFZ(1,JGR))
        ENDIF
        CALL LCMPUT(KPMAC,'FLUX-INTG',MMIX,2,FLUX(1,JGR))
        IF(JTAB(1).EQ.1 .OR. ITYPE.NE.0) THEN
          CALL LCMPUT(KPMAC,'CHI   ',MMIX,2,CHI(1,JGR))
          CALL LCMPUT(KPMAC,'NUSIGF   ',MMIX,2,SIGMA(1,JGR,3))
          CALL LCMPUT(KPMAC,'H-FACTOR',MMIX,2,SIGMA(1,JGR,4))
          CALL LCMPUT(KPMAC,'NFTOT',MMIX,2,SIGMA(1,JGR,5))
        ENDIF
*
        IL=1
        WRITE (CM,'(I2.2)') IL-1
        IPOSDE=0
        DO 190 IX=1,MMIX
          IPOS(IX)=IPOSDE+1
          DO 190 IGR=IJJ(IX,IL,JGR),IJJ(IX,IL,JGR)-NJJ(IX,IL,JGR)+1,-1
             IPOSDE=IPOSDE+1
             WORK(IPOSDE)=SCAT(IX,IL,IGR,JGR)
             SIGS(IX,IGR)=SIGS(IX,IGR)+ SCAT(IX,IL,IGR,JGR)
  190   CONTINUE
*
        CALL LCMPUT(KPMAC,'SCAT'//CM,IPOSDE,2,WORK)
        CALL LCMPUT(KPMAC,'IPOS'//CM,MMIX,1,IPOS)
        CALL LCMPUT(KPMAC,'NJJS'//CM,MMIX,1,NJJ(1,IL,JGR))
        CALL LCMPUT(KPMAC,'IJJS'//CM,MMIX,1,IJJ(1,IL,JGR))
        CALL LCMPUT(KPMAC,'SIGW'//CM,MMIX,2,SCAT(1,IL,JGR,JGR))
  210 CONTINUE
      DO 220 JGR=1,NGRP
        KPMAC=LCMDIL(JPMAC,JGR)
        IL=1
        WRITE (CM,'(I2.2)') IL-1
        CALL LCMPUT(KPMAC,'SIGS'//CM,MMIX,2,SIGS(1,JGR))
  220 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(XPURB,XPF2LB,XPF1LB,XPF2B,XPF1B,XMMDB,XMFDB,XNP9B,
     1 XSMB,XD2MB,XD1MB,XD2CB,XD1CB,XT2MB,XT1MB,XT2CB,XT1CB,XT2FB,XT1FB,
     2 XXENB,XBORB,SMACB)
      DEALLOCATE(ITEXTR,PTFUEL,PTCOOL,PDCOOL,ISFT,PSFT,BSFT,XFLUN,XSIGX,
     1 XSIGF,KTYP,INDEX,OVERVB,CHIB,JTAB,FLUXB,CPW2B,CPW1B,HNEP2,HNEP1,
     2 HSAM2,HSAM1,HXEN2,HXEN1,HISO,DENSITB,BFLUX,NJ,BRH,PW,SCATAV,
     3 IWORK,SIGAV,FLUAV,POWER,BURED,BURBG,IZONE,WORK,IJ,DIFFZ,DIFFY,
     4 DIFFX,SCAT,IPOS,ENER,CHI,FLUX,SIGS,OVERV,XBURN,NJJ,VOL,IJJ,SIGMA)
      RETURN
*
  700 FORMAT(/' AFMDRV: THE CROSS SECTIONS ARE GENERATED FOR A',
     1 ' TIME AVERAGE CALCULATION.')
  701 FORMAT(/' AFMDRV: THE CROSS SECTIONS ARE GENERATED FOR A',
     1 ' SNAPSHOT CALCULATION.')
  702 FORMAT(/' AFMDRV: POWER ARE RECOVERED FROM L_MAP.')
  703 FORMAT(/' AFMDRV: FLUX  ARE RECOVERED FROM L_MAP.')
  704 FORMAT(/' AFMDRV: BUNDLES POWER SHIFT ARE CORRECTED.')
  705 FORMAT(/' AFMDRV: BUNDLES POWER = ',F12.2,1X,'KW IS FIXED',
     1 ' BY THE USER.')
  706 FORMAT(/' AFMDRV: BUNDLES XENON = ',E14.8,1X,'IS FIXED',
     1 ' BY THE USER.')
  707 FORMAT(/' AFMDRV: LAGRANGE INTERPOLATION IS USED TO COMPUTE',
     1 ' TIME AVERAGED CROSS SECTIONS.')
  708 FORMAT(/' AFMDRV: SPLINE 3 INTERPOLATION IS USED TO COMPUTE',
     1 ' TIME AVERAGED CROSS SECTIONS.')
  709 FORMAT(/' AFMDRV: HERMITE 3 INTERPOLATION IS USED TO COMPUT',
     1 'E TIME AVERAGED CROSS SECTIONS.')
  710 FORMAT(/' AFMDRV: BUNDLES SAMARIUM = ',E14.8,1X,'IS FIXED',
     1 ' BY THE USER.')
  711 FORMAT(/' AFMDRV: BUNDLES NEPTUNIUM = ',E14.8,1X,'IS FIXED',
     1 ' BY THE USER.')
  712 FORMAT(/' AFMDRV: NOMINAL XENON IS USED.')
  713 FORMAT(/' AFMDRV: NOMINAL NEPTUNIUM IS USED.')
  714 FORMAT(/' AFMDRV: BUNDLES TFUEL = ',F12.2,1X,'K IS FIXED',
     1 ' BY THE USER.')
  715 FORMAT(/' AFMDRV: DRAGON CONCENTRATIONS ARE USED (XE135'
     1 //' NP239, SM149).')
  716 FORMAT(/' AFMDRV: ',A12,' PROFILES ARE RECOVERED FROM L_MAP.',
     1 ' PARKEY=',A12)
  717 FORMAT(/' AFMDRV: BUNDLES COOL. TEMP. TCOOL = ',F12.2,1X,
     1 'K IS FIXED BY THE USER.')
  718 FORMAT(/' AFMDRV: BUNDLES COOL. DENSITY RDCL = ',F12.9,1X,
     1 'K IS FIXED BY THE USER.')
      END
