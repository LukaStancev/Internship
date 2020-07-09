*DECK CFCFBM
      SUBROUTINE CFCFBM (TEXT1,TEXT2,IPLISU,IPLISD,IPFBM,NGRP,NBUM,NISM,
     1 NBURN,NISO,HISO,NL,IPRINT,TOTAL,ZNUG,DIFFX,DIFFY,DIFFZ,H,SCAT,
     1 MIJ,MNJ,TMREF,SMREF,DMREFX,DMREFY,DMREFZ,TOTAF,ZNUF,DXF,DYF,DZF,
     1 HF,SCATF,WORK3,REFC,TMICR,SMICR,DMICRX,DMICRY,DMICRZ,DELTA,
     1 DENSIT,TFR,TCR,TMR,XIR,TEXT,TEXTR,NB,FMICR,HMICR,FMREF,HMREF,
     1 JTAB,MIXP,V,EFJ,NXS,IXYZ,NBPARA,DBPARA)
*
*-----------------------------------------------------------------------
*
*Purpose:
* compute and store FBM coefficients.
*
*Copyright:
* Copyright (C) 1996 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): M. T. Sissaoui
*
*Parameters: input
*  IPLISU   address of the Compo object.
*  IPLISD   address of the Compo object.
*  IPFBM    address of the feedback data dase.
*  NISO     1+number of extracted isotopes.
*  TEXT1    name of the first feedback coefficient.
*  TEXT2    name of the second feedback coefficient.
*  TEXTR    name of the record.
*  TFR      reference fuel temperature.
*  TCR      reference coolant temperature.
*  TMR      reference moderator temperature.
*  NGRP     number of energy groups.
*  NB       number of feedback coefficient per parameter.
*  NL       number of Legendre orders (=1 for isotropic scattering).
*  IPRINT   print parameter. Equal to zero for no print.
*  HISO     Hollerith name information for extracted isotopes.
*  DENSIT   number densities.
*  REFC     reference number densities of the parameter
*  TOTAL    reference total macroscopic x-sections.
*  ZNUG     reference nu * fission macroscopic x-sections.
*  DIFF     reference x-directed diffusion coefficients.
*  H        reference H-FACTORS (kappa * fission mac. x-sect.).
*  SCAT     reference scattering macroscopic x-sections.
*  TMREF    reference total microscopic x-sections.
*  DMREF    reference mic. x-directed diffusion coefficients.
*  SMREF    reference scattering microscopic x-sections.
*  FMREF    reference nu * fission microscopic x-sections.
*  HMREF    reference microscopic H-FACTORS.
*  TOTAF    feedback total macroscopic x-sections.
*  ZNUF     feedback nu * fission macroscopic x-sections.
*  DIF      feedback x-directed diffusion coefficients.
*  HF       feedback H-FACTORS (kappa * fission mac. x-sect.).
*  SCATF    feedback scattering macroscopic x-sections.
*  TMICR    feedback total microscopic x-sections.
*  DMICR    feedback microscipic x-directed diffusion coefficients.
*  SMICR    feedback scattering microscopic x-sections.
*  NBPARA   Number of parameters for FBM
*  DBPARA   Values of parameters for FBM
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER TEXT1*8,TEXT2*8,TEXT(2)*12,TEXTR*12
      TYPE(C_PTR) IPLISU,IPLISD,IPFBM
      INTEGER NGRP,NBUM,NISM,NBURN,NISO,HISO(3*NISO),NL,IPRINT,
     1 MIJ(NGRP),MNJ(NGRP),NB,JTAB(NISO),MIXP,NXS,IXYZ,NBPARA
      REAL TOTAL(NGRP,NBURN),ZNUG(NGRP,NBURN),DIFFX(NGRP,NBURN),
     1 DIFFY(NGRP,NBURN),DIFFZ(NGRP,NBURN),H(NGRP,NBURN),
     2 SCAT(NBURN,NL,NGRP,NGRP),TMREF(NGRP,NBUM,NISO),
     3 SMREF(NISM,NBUM,NL,NGRP,NGRP),DMREFX(NGRP,NBUM,NISO),
     4 DMREFY(NGRP,NBUM,NISO),DMREFZ(NGRP,NBUM,NISO),
     5 TOTAF(NGRP,NBUM,NB),ZNUF(NGRP,NBUM,NB),DXF(NGRP,NBUM,NB),
     6 DYF(NGRP,NBUM,NB),DZF(NGRP,NBUM,NB),HF(NGRP,NBUM,NB),
     7 SCATF(NB,NBUM,NL,NGRP,NGRP),WORK3(NGRP*NGRP),REFC(NBUM,NISO),
     8 TMICR(NGRP,NISM,NBUM,NB),SMICR(NB,NISM,NBUM,NL,NGRP,NGRP),
     9 DMICRX(NGRP,NISM,NBUM,NB),DMICRY(NGRP,NISM,NBUM,NB),
     1 DMICRZ(NGRP,NISM,NBUM,NB),DELTA(NBUM,2),DENSIT(NISO),
     2 TFR,TCR,TMR,XIR,FMICR(NGRP,NISM,NBUM,NB),
     3 HMICR(NGRP,NISM,NBUM,NB),FMREF(NGRP,NBUM,NISO),
     4 HMREF(NGRP,NBUM,NISO),V(NBUM,8,NB),EFJ(NISO),DBPARA(NBPARA)
*----
*  LOCAL PARAMETERS
*----
      TYPE(C_PTR)   IPLIST
      INTEGER       IOUT
      PARAMETER    (IOUT=6)
      CHARACTER HMICRO*12,CM*2,TEXTB*12,TMIX(8)*8,HSMG*131
      LOGICAL LOGI,LOHIS
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IJJ,NJJ,IXS
      SAVE TMIX
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IJJ(NGRP),NJJ(NGRP),IXS(NXS))
C-----
C     PARAMETER VALUES ( TEMPERETURES, POWER AND PURITY )
C-----
      TFU=DBPARA(8)
      TCU=DBPARA(9)
      TMU=DBPARA(15)
      TFD=DBPARA(16)
      TCD=DBPARA(17)
      TMD=DBPARA(18)
      PWU=DBPARA(10)
      PWD=DBPARA(13)
      XI=DBPARA(14)
C-----
C     SET ALL THE VARIABLES TO ZERO
C-----
C
C     REAL VARIABLE
C
       DEL=0.0
       PV1U=0.0
       PV2U=0.0
       PV2UB=0.0
       PV1D=0.0
       PV2D=0.0
       PV2DB=0.0
C
      DO 600 III=1,8
        TMIX(III)=' '
 600  CONTINUE
      CALL XDRSET(V,8*NBUM*2,0.0)
      CALL XDRSET(DELTA,2*NBUM,0.0)
      CALL XDRSET(TOTAF,2*NBUM*NGRP,0.0)
      CALL XDRSET(ZNUF,2*NBUM*NGRP,0.0)
      CALL XDRSET(HF,2*NBUM*NGRP,0.0)
      CALL XDRSET(DXF,2*NBUM*NGRP,0.0)
      CALL XDRSET(DYF,2*NBUM*NGRP,0.0)
      CALL XDRSET(DZF,2*NBUM*NGRP,0.0)
      CALL XDRSET(TMICR,2*NBUM*NGRP*NISO,0.0)
      CALL XDRSET(FMICR,2*NBUM*NGRP*NISO,0.0)
      CALL XDRSET(HMICR,2*NBUM*NGRP*NISO,0.0)
      CALL XDRSET(DMICRX,2*NBUM*NGRP*NISO,0.0)
      CALL XDRSET(DMICRY,2*NBUM*NGRP*NISO,0.0)
      CALL XDRSET(DMICRZ,2*NBUM*NGRP*NISO,0.0)
      CALL XDRSET(SCATF,2*NBUM*NGRP*NGRP*NL,0.0)
      CALL XDRSET(SMICR,2*NBUM*NGRP*NGRP*NL*NISO,0.0)
C
      DO 10 IGR=1,NGRP
      IJJ(IGR)=IGR
      NJJ(IGR)=1
   10 CONTINUE
C
C     LOGICAL VARIABLE
C
      LOHIS=.FALSE.
C
C         INITIAL UNIT NUMBER
C
          IPLIST=IPLISU
C----------------------------------------------------------------------C
C-----
C        RECOVER NEUTRONICS PARAMETRES
C-----
         DO 900 J=1,NB
         IF(J.EQ.2) IPLIST=IPLISD
         CALL LCMSIX(IPLIST,TEXT(J),1)
         CALL LCMGET(IPLIST,'ISOTOPESNAME',HISO)
C
         I=1
         WRITE(TEXTB,'(4HBURN,4X,I4)') I
         CALL LCMSIX(IPLIST,TEXTB,1)
         IXYZF=0
         DO 90 ISO=1,NISO
         WRITE(HMICRO,'(3A4)') (HISO((ISO-1)*3+IH),IH=1,3)
         CALL LCMSIX(IPLIST,HMICRO,1)
         CALL LCMGET(IPLIST,'XS-SAVED',IXS)
         CALL LCMGET(IPLIST,'SCAT-SAVED',IXS(21))
         IF(IXS(18).EQ.1) IXYZF=1
         CALL LCMSIX(IPLIST,' ',2)
 90      CONTINUE
         IF(IXYZF.NE.IXYZ) THEN
         WRITE(HSMG,
     >    '(15HXS_SAVED(18) = ,I5,17HREF XS_SAVED(18)=,I5)')
     >      IXS(18),IXYZ
         CALL XABORT('CFCFBM:  INCONSISTENT NB OF FLAGS '
     1   //TEXT(J)//' IS '//HSMG//' ')
         ENDIF
         CALL LCMSIX(IPLIST,' ',2)
C
         DO 21 I=1,NBURN
         WRITE(TEXTB,'(4HBURN,4X,I4)') I
         CALL LCMSIX(IPLIST,TEXTB,1)
         CALL LCMGET(IPLIST,'ISOTOPESDENS',DENSIT)
         CALL LCMGET(IPLIST,'ISOTOPES-EFJ',EFJ)
         IF(DENSIT(1).NE.1.0) CALL XABORT('CFCFBM: DENSIT(1).NE.1.')
C
C        RECOVER FEEDBACK MACROSCOPIC  X-SECTIONS.
C
      CALL LCMSIX(IPLIST,'MACR',1)
      CALL LCMGET(IPLIST,'XS-SAVED',IXS)
      CALL LCMGET(IPLIST,'SCAT-SAVED',IXS(21))
      IF(IXS(1).EQ.1) CALL LCMGET(IPLIST,'TOTAL',TOTAF(1,I,J))
      IF(IXS(3).EQ.1) CALL LCMGET(IPLIST,'NUSIGF',ZNUF(1,I,J))
      IF(IXS(4).EQ.1) THEN
          CALL LCMGET(IPLIST,'NFTOT',HF(1,I,J))
          DO 11 IGR=1,NGRP
          HF(IGR,I,J)=HF(IGR,I,J)*EFJ(1)
  11      CONTINUE
         ENDIF
         IL=1
         WRITE (CM,'(I2.2)') IL-1
         IF(IXS(20+IL).EQ.1) THEN
         CALL LCMGET(IPLIST,'SIGS'//CM,WORK3)
         DO 231 IGR=1,NGRP
         TOTAF(IGR,I,J)= TOTAF(IGR,I,J)-WORK3(IGR)
 231     CONTINUE
         ENDIF
C
         IF(IXS(17).EQ.1) CALL LCMGET(IPLIST,'STRD  ',DXF(1,I,J))
         IF(IXS(18).EQ.1) CALL LCMGET(IPLIST,'STRD X',DXF(1,I,J))
         IF(IXS(19).EQ.1) CALL LCMGET(IPLIST,'STRD Y',DYF(1,I,J))
         IF(IXS(20).EQ.1) CALL LCMGET(IPLIST,'STRD Z',DZF(1,I,J))
C
      CALL LCMSIX(IPLIST,' ',2)
C
C     RECOVER FEEDBACK  DENSITIES.
C
      DO 40 ISO=2,NISO
      WRITE(HMICRO,'(3A4)') (HISO((ISO-1)*3+IH),IH=1,3)
      IF(HMICRO.EQ.'BMOD') THEN
        IF(TEXT1.EQ.'BOR') THEN
          DELTA(I,1)=DENSIT(ISO)- REFC(I,ISO)
          DELTA(I,2)=0.0
        ENDIF
      ELSE IF(HMICRO.EQ.'XE135') THEN
       IF(TEXT1.EQ.'XEN') THEN
         DELTA(I,1)=DENSIT(ISO)- REFC(I,ISO)
         DELTA(I,2)=0.0
         IF(I.EQ.1) DELTA(I,1)=0.0
       ELSE IF(TEXT1.EQ.'FPCH1'.OR.TEXT1.EQ.'FPCL1') THEN
         V(I,1,J)=DENSIT(ISO)- REFC(I,ISO)
         TMIX(1)='XEN'
       ENDIF
      ELSE IF(HMICRO.EQ.'SM149') THEN
       IF(TEXT1.EQ.'SM149') THEN
         DELTA(I,1)=DENSIT(ISO)- REFC(I,ISO)
         DELTA(I,2)=0.0
         IF(I.EQ.1) DELTA(I,1)=0.0
       ELSE IF(TEXT1.EQ.'FPCH1'.OR.TEXT1.EQ.'FPCL1') THEN
         V(I,2,J)=DENSIT(ISO)- REFC(I,ISO)
         TMIX(2)='SM149'
       ENDIF
      ELSE IF(HMICRO.EQ.'NP239') THEN
       IF(TEXT1.EQ.'NP239') THEN
         DELTA(I,1)=DENSIT(ISO)- REFC(I,ISO)
         DELTA(I,2)=0.0
         IF(I.EQ.1) DELTA(I,1)=0.0
       ELSE IF(TEXT1.EQ.'FPCH1'.OR.TEXT1.EQ.'FPCL1') THEN
         V(I,3,J)=DENSIT(ISO)- REFC(I,ISO)
         TMIX(3)='NP239'
       ENDIF
      ELSE IF(HMICRO.EQ.'FPC') THEN
         IF(TEXT1.EQ.'FPCH1'.OR.TEXT1.EQ.'FPCL1') THEN
         DELTA(I,J)=DENSIT(ISO)- REFC(I,ISO)
         ENDIF
      ELSE IF(HMICRO.EQ.'CWAT') THEN
        IF(TEXT1.EQ.'D1C') THEN
          IF(J.EQ.1) THEN
          PV1U=DENSIT(ISO)- REFC(I,ISO)
          PV2U=PV1U*PV1U
          PV2UB=PV1U*PV1U
          DELTA(I,J)=PV1U
          ELSE
          PV1D=DENSIT(ISO)- REFC(I,ISO)
          PV2D=PV1D*PV1D
          PV2DB=PV1D*PV1D
          DELTA(I,J)=PV1D
          ENDIF
        ENDIF
        IF(TEXT1.EQ.'MIXFD'.OR.TEXT1.EQ.'MIXMD') THEN
        DELTA(I,1)=DENSIT(ISO)- REFC(I,ISO)
        DELTA(I,2)=0.0
        V(I,1,J)=DENSIT(ISO)- REFC(I,ISO)
        V(I,2,J)=V(I,1,J)*V(I,1,J)
        TMIX(1)='D1C'
        TMIX(2)='D2C'
        ENDIF
      ELSE IF(HMICRO.EQ.'MWAT') THEN
        IF(TEXT1.EQ.'D1M') THEN
          IF(J.EQ.1) THEN
          PV1U=ALOG(DENSIT(ISO)/REFC(I,ISO))
          PV2U=1.0/DENSIT(ISO) - 1.0/REFC(I,ISO)
          PV2UB=DENSIT(ISO)- REFC(I,ISO)
          DELTA(I,J)=PV2UB
          ELSE
          PV1D=ALOG(DENSIT(ISO)/REFC(I,ISO))
          PV2D=1.0/DENSIT(ISO) - 1.0/REFC(I,ISO)
          PV2DB=DENSIT(ISO)- REFC(I,ISO)
          DELTA(I,J)=PV2DB
          ENDIF
        ELSE IF(TEXT1.EQ.'PUR') THEN
          DELTA(I,1)=(XI-XIR)*REFC(I,ISO)
          DELTA(I,2)=0.0
        ENDIF
      ENDIF
C
C     RECOVER FEEDBACK MICROSCOPIC  X-SECTIONS.
C
      CALL LCMSIX(IPLIST,HMICRO,1)
      CALL LCMGET(IPLIST,'XS-SAVED',IXS)
      CALL LCMGET(IPLIST,'SCAT-SAVED',IXS(21))
      IF(IXS(1).EQ.1) CALL LCMGET(IPLIST,'TOTAL',TMICR(1,ISO,I,J))
      IF(IXS(3).EQ.1) CALL LCMGET(IPLIST,'NUSIGF',FMICR(1,ISO,I,J))
      IF(IXS(4).EQ.1) THEN
          CALL LCMGET(IPLIST,'NFTOT',HMICR(1,ISO,I,J))
          DO 25 IGR=1,NGRP
          HMICR(IGR,ISO,I,J)=HMICR(IGR,ISO,I,J)*EFJ(ISO)
  25      CONTINUE
         ENDIF
         IL=1
         WRITE (CM,'(I2.2)') IL-1
         IF(IXS(20+IL).EQ.1) THEN
         CALL LCMGET(IPLIST,'SIGS'//CM,WORK3)
         DO 233 IGR=1,NGRP
         TMICR(IGR,ISO,I,J)= TMICR(IGR,ISO,I,J)-WORK3(IGR)
 233     CONTINUE
         ENDIF
C
C
      IF(IXS(17).EQ.1) CALL LCMGET(IPLIST,'STRD  ',DMICRX(1,ISO,I,J))
      IF(IXS(18).EQ.1) CALL LCMGET(IPLIST,'STRD X',DMICRX(1,ISO,I,J))
      IF(IXS(19).EQ.1) CALL LCMGET(IPLIST,'STRD Y',DMICRY(1,ISO,I,J))
      IF(IXS(20).EQ.1) CALL LCMGET(IPLIST,'STRD Z',DMICRZ(1,ISO,I,J))
C
C     ADD THE CONTRIBUTION OF MIC. X-SECT. IN MAC. X-S
C
      DO 20 IGR=1,NGRP
       TOTAF(IGR,I,J)=TOTAF(IGR,I,J)+DENSIT(ISO)*TMICR(IGR,ISO,I,J)
       DXF(IGR,I,J) =DXF(IGR,I,J)  +DENSIT(ISO)*DMICRX(IGR,ISO,I,J)
       DYF(IGR,I,J) =DYF(IGR,I,J)  +DENSIT(ISO)*DMICRY(IGR,ISO,I,J)
       DZF(IGR,I,J) =DZF(IGR,I,J)  +DENSIT(ISO)*DMICRZ(IGR,ISO,I,J)
 20   CONTINUE
      IF(JTAB(ISO).EQ.1) THEN
       DO 30 IGR=1,NGRP
       ZNUF(IGR,I,J)=ZNUF(IGR,I,J)+DENSIT(ISO)*FMICR(IGR,ISO,I,J)
       HF(IGR,I,J)  =HF(IGR,I,J)  +DENSIT(ISO)*HMICR(IGR,ISO,I,J)
 30   CONTINUE
      ENDIF
      CALL LCMSIX(IPLIST,' ',2)
   40 CONTINUE
C
C     RECOVER MACROSCOPIC  SCATTERING X-SECTIONS.
C
      CALL LCMSIX(IPLIST,'MACR',1)
      CALL LCMGET(IPLIST,'XS-SAVED',IXS)
      CALL LCMGET(IPLIST,'SCAT-SAVED',IXS(21))
      IL=1
      WRITE (CM,'(I2.2)') IL-1
      IF(IXS(20+IL).EQ.1) THEN
         CALL LCMGET(IPLIST,'SCAT'//CM,WORK3)
         CALL LCMGET(IPLIST,'NJJS'//CM,NJJ)
         CALL LCMGET(IPLIST,'IJJS'//CM,IJJ)
         IGAR=0
         DO 120 JGR=1,NGRP
         DO 120 IGR=IJJ(JGR),IJJ(JGR)-NJJ(JGR)+1,-1
         IGAR=IGAR+1
         SCATF(J,I,IL,IGR,JGR)=WORK3(IGAR)
  120    CONTINUE
      ENDIF
C
      CALL LCMSIX(IPLIST,' ',2)
C
C     RECOVER MICROSCOPIC CONTRIBUTIONS OF SCATTERING X-SECTIONS.
C
      DO 160 ISO=2,NISO
      WRITE(HMICRO,'(3A4)') (HISO((ISO-1)*3+IH),IH=1,3)
      CALL LCMSIX(IPLIST,HMICRO,1)
      CALL LCMGET(IPLIST,'XS-SAVED',IXS)
      CALL LCMGET(IPLIST,'SCAT-SAVED',IXS(21))
      IL=1
      WRITE (CM,'(I2.2)') IL-1
      IF(IXS(20+IL).EQ.1) THEN
         CALL LCMGET(IPLIST,'SCAT'//CM,WORK3)
         CALL LCMGET(IPLIST,'NJJS'//CM,NJJ)
         CALL LCMGET(IPLIST,'IJJS'//CM,IJJ)
         IGAR=0
         DO 140 JGR=1,NGRP
         DO 140 IGR=IJJ(JGR),IJJ(JGR)-NJJ(JGR)+1,-1
         IGAR=IGAR+1
         SMICR(J,ISO,I,IL,IGR,JGR)=WORK3(IGAR)
         SCATF(J,I,IL,IGR,JGR)=SCATF(J,I,IL,IGR,JGR)+
     1   DENSIT(ISO)*WORK3(IGAR)
  140 CONTINUE
      ENDIF
C
      CALL LCMSIX(IPLIST,' ',2)
  160 CONTINUE
C      GOING UP FOR BURN
       CALL LCMSIX(IPLIST,' ',2)
 21       CONTINUE
       CALL LCMSIX(IPLIST,' ',2)
 900  CONTINUE
C----------------------------------------------------------------------C
C                                                                      C
C      BEGIN THE COEFFICIENTS CALCULATION                              C
C                                                                      C
C----------------------------------------------------------------------C
       DT=0.0
       IF(TEXT1.EQ.'T1F') THEN
       PV1U=SQRT(TFU)-SQRT(TFR)
       PV2U=TFU-TFR
       PV2UB=PV2U
       PV1D=SQRT(TFD)-SQRT(TFR)
       PV2D=TFD-TFR
       PV2DB=PV2D
       ELSE IF(TEXT1.EQ.'T1C') THEN
       PV1U=ALOG(TCU/TCR)
       PV2U=1.0/TCU - 1.0/TCR
       PV2UB=PV2U
       PV1D=ALOG(TCD/TCR)
       PV2D=1.0/TCD - 1.0/TCR
       PV2DB=PV2D
       ELSE IF(TEXT1.EQ.'T1M') THEN
       PV1U=ALOG(TMU/TMR)
       PV2U=1.0/TMU - 1.0/TMR
       PV2UB=PV2U
       PV1D=ALOG(TMD/TMR)
       PV2D=1.0/TMD - 1.0/TMR
       PV2DB=PV2D
       ELSE IF(TEXT1.EQ.'MIXMD') THEN
       DT=ALOG(TCU/TCR)
       DO 152 I=1,NBURN
       V(I,3,1)=ALOG(TCU/TCR)
       V(I,4,1)=1.0/TCU - 1.0/TCR
 152   CONTINUE
       TMIX(1)='D1C'
       TMIX(2)='D2C'
       TMIX(3)='T1C'
       TMIX(4)='T2C'
       ELSE IF(TEXT1.EQ.'MIXFD') THEN
       DT=SQRT(TFU)-SQRT(TFR)
       DO 151 I=1,NBURN
       V(I,3,1)=SQRT(TFU)-SQRT(TFR)
       V(I,4,1)=TFU-TFR
 151   CONTINUE
       TMIX(1)='D1C'
       TMIX(2)='D2C'
       TMIX(3)='T1F'
       TMIX(4)='T2F'
       ENDIF
C
C----------------------------------------------------------------------C
C
C      COMPUTE DELTA SIGMA
C
       DO 801 I=1,NBURN
       DO 801 J=1,NB
       DO 807 IGR=1,NGRP
       TOTAF(IGR,I,J)=TOTAF(IGR,I,J)-TOTAL(IGR,I)
        DXF(IGR,I,J)=DXF(IGR,I,J)-DIFFX(IGR,I)
        DYF(IGR,I,J)=DYF(IGR,I,J)-DIFFY(IGR,I)
        DZF(IGR,I,J)=DZF(IGR,I,J)-DIFFZ(IGR,I)
       ZNUF(IGR,I,J)=ZNUF(IGR,I,J)-ZNUG(IGR,I)
       HF(IGR,I,J)=HF(IGR,I,J)-H(IGR,I)
       DO 812 ISO=2,NISO
        TMICR(IGR,ISO,I,J) =TMICR(IGR,ISO,I,J) - TMREF(IGR,I,ISO)
        DMICRX(IGR,ISO,I,J)=DMICRX(IGR,ISO,I,J)- DMREFX(IGR,I,ISO)
        DMICRY(IGR,ISO,I,J)=DMICRY(IGR,ISO,I,J)- DMREFY(IGR,I,ISO)
        DMICRZ(IGR,ISO,I,J)=DMICRZ(IGR,ISO,I,J)- DMREFZ(IGR,I,ISO)
       IF(JTAB(ISO).EQ.1) THEN
       FMICR(IGR,ISO,I,J)=FMICR(IGR,ISO,I,J)- FMREF(IGR,I,ISO)
       HMICR(IGR,ISO,I,J)=HMICR(IGR,ISO,I,J)- HMREF(IGR,I,ISO)
       ENDIF
 812   CONTINUE
       IL=1
       DO 813 JGR=1,NGRP
       SCATF(J,I,IL,IGR,JGR)=SCATF(J,I,IL,IGR,JGR)-
     1 SCAT(I,IL,IGR,JGR)
       DO 813 ISO=2,NISO
       SMICR(J,ISO,I,IL,IGR,JGR)=SMICR(J,ISO,I,IL,IGR,JGR)-
     1  SMREF(ISO,I,IL,IGR,JGR)
 813   CONTINUE
 807   CONTINUE
C
C       CORRECTION OF MACRO. X-SECTIONS
C
        LOGI=.FALSE.
        DO 872 ISO=2,NISO
        WRITE(HMICRO,'(3A4)') (HISO((ISO-1)*3+IH),IH=1,3)
        IF(HMICRO.EQ.'BMOD'.AND.TEXT1.EQ.'BOR') THEN
        LOGI =.TRUE.
        DELC=DELTA(I,J)
        ELSE IF(HMICRO.EQ.'XE135'.AND.TEXT1.EQ.'XEN') THEN
        LOGI =.TRUE.
        DELC=DELTA(I,J)
        ELSE IF(HMICRO.EQ.'SM149'.AND.TEXT1.EQ.'SM149') THEN
        LOGI =.TRUE.
        DELC=DELTA(I,J)
        ELSE IF(HMICRO.EQ.'NP239'.AND.TEXT1.EQ.'NP239') THEN
        LOGI =.TRUE.
        DELC=DELTA(I,J)
        ELSE IF(HMICRO.EQ.'CWAT'.AND.TEXT1.EQ.'D1C') THEN
        LOGI =.TRUE.
        DELC=DELTA(I,J)
        ELSE IF(HMICRO.EQ.'MWAT'.AND.TEXT1.EQ.'D1M') THEN
        LOGI =.TRUE.
        DELC=DELTA(I,J)
        ELSE IF(HMICRO.EQ.'CWAT'.AND.TEXT1.EQ.'MIXFD') THEN
        LOGI =.TRUE.
        DELC=DELTA(I,J)
        ELSE IF(HMICRO.EQ.'CWAT'.AND.TEXT1.EQ.'MIXMD') THEN
        LOGI =.TRUE.
        DELC=DELTA(I,J)
        ELSE IF(HMICRO.EQ.'XE135'.AND.TEXT1.EQ.'FPCH1') THEN
        LOGI =.TRUE.
        DELC=V(I,1,J)
        ELSE IF(HMICRO.EQ.'XE135'.AND.TEXT1.EQ.'FPCL1') THEN
        LOGI =.TRUE.
        DELC=V(I,1,J)
        ELSE IF(HMICRO.EQ.'SM149'.AND.TEXT1.EQ.'FPCH1') THEN
        LOGI =.TRUE.
        DELC=V(I,2,J)
        ELSE IF(HMICRO.EQ.'SM149'.AND.TEXT1.EQ.'FPCL1') THEN
        LOGI =.TRUE.
        DELC=V(I,2,J)
        ELSE IF(HMICRO.EQ.'NP239'.AND.TEXT1.EQ.'FPCH1') THEN
        LOGI =.TRUE.
        DELC=V(I,3,J)
        ELSE IF(HMICRO.EQ.'NP239'.AND.TEXT1.EQ.'FPCL1') THEN
        LOGI =.TRUE.
        DELC=V(I,3,J)
        ELSE IF(HMICRO.EQ.'FPC'.AND.TEXT1.EQ.'FPCH1') THEN
        LOGI =.TRUE.
        DELC=DELTA(I,J)
        ELSE IF(HMICRO.EQ.'FPC'.AND.TEXT1.EQ.'FPCL1') THEN
        LOGI =.TRUE.
        DELC=DELTA(I,J)
        ENDIF
C
        IF(LOGI) THEN
        DO 808 IGR=1,NGRP
        TOTAF(IGR,I,J)=TOTAF(IGR,I,J)-TMICR(IGR,ISO,I,J)*DELC
         DXF(IGR,I,J)=DXF(IGR,I,J)-DMICRX(IGR,ISO,I,J)*DELC
         DYF(IGR,I,J)=DYF(IGR,I,J)-DMICRY(IGR,ISO,I,J)*DELC
         DZF(IGR,I,J)=DZF(IGR,I,J)-DMICRZ(IGR,ISO,I,J)*DELC
        IF(JTAB(ISO).EQ.1) THEN
        ZNUF(IGR,I,J)= ZNUF(IGR,I,J)-FMICR(IGR,ISO,I,J)*DELC
        HF(IGR,I,J)=HF(IGR,I,J)-HMICR(IGR,ISO,I,J)*DELC
        ENDIF
        IL=1
        DO 873 JGR=1,NGRP
        SCATF(J,I,IL,IGR,JGR)=SCATF(J,I,IL,IGR,JGR)-
     1  SMICR(J,ISO,I,IL,IGR,JGR)*DELC
 873   CONTINUE
 808   CONTINUE
        ENDIF
       LOGI=.FALSE.
 872   CONTINUE
 801   CONTINUE
C----------------------------------------------------------------------C
C     'MIXMD' AND 'MIXFD'
C     TAKE OFF THE INDIVIDUAL VARIATION CONTRIBUTION OF:
C     FUEL     TEMPERATURE
C     COOLANT  TEMPERATURE
C     COOLANT  DENSITY
C----------------------------
C     'FPCH1' AND 'FPCL1'
C     XENON     CONCENTRATION
C     SAMARIUM  CONCENTRATION
C     NEPTUNIUM CONCENTRATION
C
       IF(MIXP.EQ.1) THEN
       IF(TEXT1.EQ.'MIXMD'.OR.TEXT1.EQ.'MIXFD') NCOR=4
       IF(TEXT1.EQ.'FPCH1'.OR.TEXT1.EQ.'FPCL1') NCOR=3
       IF(ABS(IPRINT) .GT. 5) THEN
         WRITE(IOUT,6000) TEXT1,NCOR,(TMIX(II),II=1,NCOR)
       ENDIF 
       CALL LCMSIX(IPFBM,TEXTR,1)
       DO 711 I=1,NBURN
        WRITE(TEXTB,'(4HBURN,4X,I4)') I
         CALL LCMSIX(IPFBM,TEXTB,1)
         CALL LCMSIX(IPFBM,'MACR',1)
         CALL LCMSIX(IPFBM,'ABS',1)
         DO 332 II=1,NCOR
         DO 332 J=1,NB
         CALL LCMGET(IPFBM,TMIX(II),WORK3)
         DO 331 IGR=1,NGRP
 331     TOTAF(IGR,I,J)= TOTAF(IGR,I,J)-WORK3(IGR)*V(I,II,J)
 332   CONTINUE
       CALL LCMSIX(IPFBM,' ',2)
C
       IF(IXYZ.EQ.0) THEN
       CALL LCMSIX(IPFBM,'STRD  ',1)
       DO 333 II=1,NCOR
       DO 333 J=1,NB
       CALL LCMGET(IPFBM,TMIX(II),WORK3)
       DO 334 IGR=1,NGRP
 334    DXF(IGR,I,J)=  DXF(IGR,I,J)-WORK3(IGR)*V(I,II,J)
 333   CONTINUE
      CALL LCMSIX(IPFBM,' ',2)
C
       ELSE IF(IXYZ.EQ.1) THEN
       CALL LCMSIX(IPFBM,'STRD X',1)
       DO 50 II=1,NCOR
       DO 50 J=1,NB
       CALL LCMGET(IPFBM,TMIX(II),WORK3)
       DO 51 IGR=1,NGRP
 51    DXF(IGR,I,J)=  DXF(IGR,I,J)-WORK3(IGR)*V(I,II,J)
 50   CONTINUE
      CALL LCMSIX(IPFBM,' ',2)
C
       CALL LCMSIX(IPFBM,'STRD Y',1)
       DO 52 II=1,NCOR
       DO 52 J=1,NB
       CALL LCMGET(IPFBM,TMIX(II),WORK3)
       DO 53 IGR=1,NGRP
 53    DYF(IGR,I,J)=  DYF(IGR,I,J)-WORK3(IGR)*V(I,II,J)
 52   CONTINUE
      CALL LCMSIX(IPFBM,' ',2)
C
       CALL LCMSIX(IPFBM,'STRD Z',1)
       DO 54 II=1,NCOR
       DO 54 J=1,NB
       CALL LCMGET(IPFBM,TMIX(II),WORK3)
       DO 55 IGR=1,NGRP
 55    DZF(IGR,I,J)=  DZF(IGR,I,J)-WORK3(IGR)*V(I,II,J)
 54   CONTINUE
      CALL LCMSIX(IPFBM,' ',2)
      ENDIF
C
      IF(JTAB(1).EQ.1) THEN
      CALL LCMSIX(IPFBM,'NUSIGF',1)
      DO 335 II=1,NCOR
      DO 335 J=1,NB
       CALL LCMGET(IPFBM,TMIX(II),WORK3)
       DO 336 IGR=1,NGRP
 336   ZNUF(IGR,I,J)= ZNUF(IGR,I,J)-WORK3(IGR)*V(I,II,J)
 335   CONTINUE
      CALL LCMSIX(IPFBM,' ',2)
C
      CALL LCMSIX(IPFBM,'H-FACTORS',1)
       DO 337 II=1,NCOR
       DO 337 J=1,NB
       CALL LCMGET(IPFBM,TMIX(II),WORK3)
       DO 338 IGR=1,NGRP
 338   HF(IGR,I,J)= HF(IGR,I,J)-WORK3(IGR)*V(I,II,J)
 337   CONTINUE
      CALL LCMSIX(IPFBM,' ',2)
      ENDIF
C
      IL=1
      WRITE (CM,'(I2)') IL-1
      CALL LCMSIX(IPFBM,'SCAT'//CM,1)
      CALL LCMGET(IPFBM,'NJJ',NJJ)
      CALL LCMGET(IPFBM,'IJJ',IJJ)
      DO 340 II=1,NCOR
      DO 340 J=1,NB
        CALL LCMGET(IPFBM,TMIX(II),WORK3)
        IGAR=0
        DO 341 JGR=1,NGRP
        DO 341 IGR=IJJ(JGR),IJJ(JGR)-NJJ(JGR)+1,-1
          IGAR=IGAR+1
          SCATF(J,I,IL,IGR,JGR)=SCATF(J,I,IL,IGR,JGR)-
     1    WORK3(IGAR)*V(I,II,J)
 341    CONTINUE
 340  CONTINUE
      CALL LCMSIX(IPFBM,' ',2)
C
C      GO UP FOR MACR
       CALL LCMSIX(IPFBM,' ',2)
C
C     MICROSCOPIC X-SECTION CORRECTION
C
      DO 360 ISO=2,NISO
      WRITE(HMICRO,'(3A4)') (HISO((ISO-1)*3+IH),IH=1,3)
      CALL LCMSIX(IPFBM,HMICRO,1)
C
      CALL LCMSIX(IPFBM,'ABS',1)
       DO 342 II=1,NCOR
       DO 342 J=1,NB
       CALL LCMGET(IPFBM,TMIX(II),WORK3)
       DO 343 IGR=1,NGRP
 343   TMICR(IGR,ISO,I,J)=TMICR(IGR,ISO,I,J)-
     1 WORK3(IGR)*V(I,II,J)
 342   CONTINUE
      CALL LCMSIX(IPFBM,' ',2)
C
      IF(IXYZ.EQ.0) THEN
      CALL LCMSIX(IPFBM,'STRD  ',1)
      DO 344 II=1,NCOR
      DO 344 J=1,NB
       CALL LCMGET(IPFBM,TMIX(II),WORK3)
       DO 345 IGR=1,NGRP
 345   DMICRX(IGR,ISO,I,J)=DMICRX(IGR,ISO,I,J)-
     1 WORK3(IGR)*V(I,II,J)
 344   CONTINUE
       CALL LCMSIX(IPFBM,' ',2)
C
       ELSE IF(IXYZ.EQ.1) THEN
       CALL LCMSIX(IPFBM,'STRD X',1)
       DO 60 II=1,NCOR
       DO 60 J=1,NB
       CALL LCMGET(IPFBM,TMIX(II),WORK3)
       DO 61 IGR=1,NGRP
 61    DMICRX(IGR,ISO,I,J)=DMICRX(IGR,ISO,I,J)-
     1 WORK3(IGR)*V(I,II,J)
 60   CONTINUE
       CALL LCMSIX(IPFBM,' ',2)
C
       CALL LCMSIX(IPFBM,'STRD Y',1)
       DO 62 II=1,NCOR
       DO 62 J=1,NB
       CALL LCMGET(IPFBM,TMIX(II),WORK3)
       DO 63 IGR=1,NGRP
 63    DMICRY(IGR,ISO,I,J)=DMICRY(IGR,ISO,I,J)-
     1 WORK3(IGR)*V(I,II,J)
 62   CONTINUE
       CALL LCMSIX(IPFBM,' ',2)
C
       CALL LCMSIX(IPFBM,'STRD Z',1)
       DO 64 II=1,NCOR
       DO 64 J=1,NB
       CALL LCMGET(IPFBM,TMIX(II),WORK3)
       DO 65 IGR=1,NGRP
 65    DMICRZ(IGR,ISO,I,J)=DMICRZ(IGR,ISO,I,J)-
     1 WORK3(IGR)*V(I,II,J)
 64   CONTINUE
       CALL LCMSIX(IPFBM,' ',2)
       ENDIF
C
      IF(JTAB(ISO).EQ.1) THEN
      CALL LCMSIX(IPFBM,'NUSIGF',1)
      DO 346 II=1,NCOR
      DO 346 J=1,NB
       CALL LCMGET(IPFBM,TMIX(II),WORK3)
       DO 347 IGR=1,NGRP
 347   FMICR(IGR,ISO,I,J)=FMICR(IGR,ISO,I,J)-
     1 WORK3(IGR)*V(I,II,J)
 346   CONTINUE
      CALL LCMSIX(IPFBM,' ',2)
C
      CALL LCMSIX(IPFBM,'H-FACTORS',1)
      DO 348 II=1,NCOR
      DO 348 J=1,NB
      CALL LCMGET(IPFBM,TMIX(II),WORK3)
      DO 349 IGR=1,NGRP
 349  HMICR(IGR,ISO,I,J)= HMICR(IGR,ISO,I,J)-
     1 WORK3(IGR)*V(I,II,J)
 348  CONTINUE
      CALL LCMSIX(IPFBM,' ',2)
      ENDIF
C
      IL=1
      WRITE (CM,'(I2)') IL-1
      CALL LCMSIX(IPFBM,'SCAT'//CM,1)
       DO 740 II=1,NCOR
       DO 740 J=1,NB
      CALL LCMGET(IPFBM,TMIX(II),WORK3)
      IGAR=0
      DO 741 JGR=1,NGRP
      DO 741 IGR=MIJ(JGR),MIJ(JGR)-MNJ(JGR)+1,-1
      IGAR=IGAR+1
       SMICR(J,ISO,I,IL,IGR,JGR)=SMICR(J,ISO,I,IL,IGR,JGR)-
     1 WORK3(IGAR)*V(I,II,J)
 741   CONTINUE
 740   CONTINUE
       CALL LCMSIX(IPFBM,' ',2)
C
C
      CALL LCMSIX(IPFBM,' ',2)
 360  CONTINUE
C     GO UP FOR BURNUP
        CALL LCMSIX(IPFBM,' ',2)
 711   CONTINUE
       CALL LCMSIX(IPFBM,' ',2)
       ENDIF
C
C      END OF INDIVIDUAL CORRECTION
C----------------------------------------------------------------------C
C-----
C     INVERT THE FEEDBACK FORMULAS
C-----
       DO 811 I=1,NBURN
       DO 811 IGR=1,NGRP
       IF(NB.EQ.1) THEN
C
C      ONLY ONE COEFFICIENT IS REQUIRED (NB=1) FOR:
C      BORON CONCENTRATION
C      XENON CONCENTRATION
C      SAMARIUM CONCENTRATION
C      NEPTUNIUM CONCENTRATION
C      MODERATOR PURITY
C
       IF(TEXT1.EQ.'MIXMD'.AND.IGR.EQ.1) DELTA(I,1)=DELTA(I,1)*DT
       IF(TEXT1.EQ.'MIXFD'.AND.IGR.EQ.1) DELTA(I,1)=DELTA(I,1)*DT
        IF(DELTA(I,1).NE.0.0) THEN
        TOTAF(IGR,I,1)=TOTAF(IGR,I,1)/DELTA(I,1)
         DXF(IGR,I,1)=DXF(IGR,I,1)/DELTA(I,1)
         DYF(IGR,I,1)=DYF(IGR,I,1)/DELTA(I,1)
         DZF(IGR,I,1)=DZF(IGR,I,1)/DELTA(I,1)
        ZNUF(IGR,I,1)=ZNUF(IGR,I,1)/DELTA(I,1)
        HF(IGR,I,1)=HF(IGR,I,1)/DELTA(I,1)
        DO 802 ISO=2,NISO
        TMICR(IGR,ISO,I,1)=TMICR(IGR,ISO,I,1)/DELTA(I,1)
         DMICRX(IGR,ISO,I,1)=DMICRX(IGR,ISO,I,1)/DELTA(I,1)
         DMICRY(IGR,ISO,I,1)=DMICRY(IGR,ISO,I,1)/DELTA(I,1)
         DMICRZ(IGR,ISO,I,1)=DMICRZ(IGR,ISO,I,1)/DELTA(I,1)
         IF(JTAB(ISO).EQ.1) THEN
         FMICR(IGR,ISO,I,1)=FMICR(IGR,ISO,I,1)/DELTA(I,1)
         HMICR(IGR,ISO,I,1)=HMICR(IGR,ISO,I,1)/DELTA(I,1)
         ENDIF
 802    CONTINUE
        IL=1
        DO 803 JGR=1,NGRP
        SCATF(1,I,IL,IGR,JGR)=SCATF(1,I,IL,IGR,JGR)/DELTA(I,1)
        DO 803 ISO=2,NISO
       SMICR(1,ISO,I,IL,IGR,JGR)=SMICR(1,ISO,I,IL,IGR,JGR)/DELTA(I,1)
 803    CONTINUE
        ELSE
         TOTAF(IGR,I,1)=0.0
         DXF(IGR,I,1) =0.0
         DYF(IGR,I,1) =0.0
         DZF(IGR,I,1) =0.0
C
        ZNUF(IGR,I,1) =0.0
        HF(IGR,I,1)  =0.0
        DO 552 ISO=2,NISO
         TMICR(IGR,ISO,I,1) =0.0
         DMICRX(IGR,ISO,I,1)=0.0
         DMICRY(IGR,ISO,I,1)=0.0
         DMICRZ(IGR,ISO,I,1)=0.0
         IF(JTAB(ISO).EQ.1) THEN
         FMICR(IGR,ISO,I,1)=0.0
         HMICR(IGR,ISO,I,1)=0.0
         ENDIF
 552    CONTINUE
        DO 553 IL=1,NL
        DO 553 JGR=1,NGRP
        SCATF(1,I,IL,IGR,JGR)=0.0
        DO 553 ISO=2,NISO
        SMICR(1,ISO,I,IL,IGR,JGR)=0.0
 553    CONTINUE
        ENDIF
C
       ELSE IF(NB.EQ.2) THEN
C
C      INVERT THE FEEDBACK FORMULAS
C      TWO FBM COEFFICIENTS ARE COMPUTED
C      TEMPERATURES
C      DENSITIES
C      POWER LEVEL
C
          IF(TEXT1.EQ.'FPCH1'.OR.TEXT1.EQ.'FPCL1') THEN
          PV1U=DELTA(I,1)
          PV2U=PV1U*PV1U
          PV2UB=PV2U
          PV1D=DELTA(I,2)
          PV2D=PV1D*PV1D
          PV2DB=PV2D
          ENDIF
C
       TX=PV2U*PV1D - PV2D*PV1U
       TXB=PV2UB*PV1D - PV2DB*PV1U
C
       IF(TX.NE.0.0.AND.TXB.NE.0.0) THEN
       TOTAF(IGR,I,2)=(TOTAF(IGR,I,1)*PV1D-TOTAF(IGR,I,2)*PV1U)/TX
       TOTAF(IGR,I,1)=(TOTAF(IGR,I,1) - TOTAF(IGR,I,2)*PV2U)/PV1U
C
       DXF(IGR,I,2)=(DXF(IGR,I,1)*PV1D -DXF(IGR,I,2)*PV1U)/TXB
       DXF(IGR,I,1)=(DXF(IGR,I,1) - DXF(IGR,I,2)*PV2UB)/PV1U
       DYF(IGR,I,2)=(DYF(IGR,I,1)*PV1D -DYF(IGR,I,2)*PV1U)/TXB
       DYF(IGR,I,1)=(DYF(IGR,I,1) - DYF(IGR,I,2)*PV2UB)/PV1U
       DZF(IGR,I,2)=(DZF(IGR,I,1)*PV1D -DZF(IGR,I,2)*PV1U)/TXB
       DZF(IGR,I,1)=(DZF(IGR,I,1) - DZF(IGR,I,2)*PV2UB)/PV1U
C
       ZNUF(IGR,I,2)=(ZNUF(IGR,I,1)*PV1D - ZNUF(IGR,I,2)*PV1U)/TX
       ZNUF(IGR,I,1)=(ZNUF(IGR,I,1) - ZNUF(IGR,I,2)*PV2U)/PV1U
C
       HF(IGR,I,2)=(HF(IGR,I,1)*PV1D - HF(IGR,I,2)*PV1U)/TX
       HF(IGR,I,1)=(HF(IGR,I,1) - HF(IGR,I,2)*PV2U)/PV1U
C
       DO 814 ISO=2,NISO
       TMICR(IGR,ISO,I,2)=(TMICR(IGR,ISO,I,1)*PV1D -
     1 TMICR(IGR,ISO,I,2)*PV1U)/TX
       TMICR(IGR,ISO,I,1)=(TMICR(IGR,ISO,I,1) -
     1 TMICR(IGR,ISO,I,2)*PV2U)/PV1U
C
       DMICRX(IGR,ISO,I,2)=(DMICRX(IGR,ISO,I,1)*PV1D -
     1 DMICRX(IGR,ISO,I,2)*PV1U)/TX
       DMICRX(IGR,ISO,I,1)=(DMICRX(IGR,ISO,I,1) -
     1 DMICRX(IGR,ISO,I,2)*PV2U)/PV1U
       DMICRY(IGR,ISO,I,2)=(DMICRY(IGR,ISO,I,1)*PV1D -
     1 DMICRY(IGR,ISO,I,2)*PV1U)/TX
       DMICRY(IGR,ISO,I,1)=(DMICRY(IGR,ISO,I,1) -
     1 DMICRY(IGR,ISO,I,2)*PV2U)/PV1U
       DMICRZ(IGR,ISO,I,2)=(DMICRZ(IGR,ISO,I,1)*PV1D -
     1 DMICRZ(IGR,ISO,I,2)*PV1U)/TX
       DMICRZ(IGR,ISO,I,1)=(DMICRZ(IGR,ISO,I,1) -
     1 DMICRZ(IGR,ISO,I,2)*PV2U)/PV1U
C
        IF(JTAB(ISO).EQ.1) THEN
        FMICR(IGR,ISO,I,2)=(FMICR(IGR,ISO,I,1)*PV1D -
     1  FMICR(IGR,ISO,I,2)*PV1U)/TX
        FMICR(IGR,ISO,I,1)=(FMICR(IGR,ISO,I,1) -
     1  FMICR(IGR,ISO,I,2)*PV2U)/PV1U
C
        HMICR(IGR,ISO,I,2)=(HMICR(IGR,ISO,I,1)*PV1D -
     1  HMICR(IGR,ISO,I,2)*PV1U)/TX
        HMICR(IGR,ISO,I,1)=(HMICR(IGR,ISO,I,1) -
     1  HMICR(IGR,ISO,I,2)*PV2U)/PV1U
        ENDIF
 814   CONTINUE
C
       IL=1
       DO 815 JGR=1,NGRP
       SCATF(2,I,IL,IGR,JGR)=(SCATF(1,I,IL,IGR,JGR)*PV1D-
     1 SCATF(2,I,IL,IGR,JGR)*PV1U)/TXB
       SCATF(1,I,IL,IGR,JGR)=(SCATF(1,I,IL,IGR,JGR) -
     1 SCATF(2,I,IL,IGR,JGR)*PV2UB)/PV1U
       DO 815 ISO=2,NISO
       SMICR(2,ISO,I,IL,IGR,JGR)=(SMICR(1,ISO,I,IL,IGR,JGR)*PV1D-
     1 SMICR(2,ISO,I,IL,IGR,JGR)*PV1U)/TX
       SMICR(1,ISO,I,IL,IGR,JGR)=(SMICR(1,ISO,I,IL,IGR,JGR) -
     1 SMICR(2,ISO,I,IL,IGR,JGR)*PV2U)/PV1U
 815   CONTINUE
         ELSE
         DO 876 J=1,NB
         TOTAF(IGR,I,J)=0.0
         DXF(IGR,I,J)  =0.0
         DYF(IGR,I,J)  =0.0
         DZF(IGR,I,J)  =0.0
         ZNUF(IGR,I,J)=0.0
         HF(IGR,I,J)=0.0
         DO 874 ISO=2,NISO
          TMICR(IGR,ISO,I,J) =0.0
          DMICRX(IGR,ISO,I,J)=0.0
          DMICRY(IGR,ISO,I,J)=0.0
          DMICRZ(IGR,ISO,I,J)=0.0
         IF(JTAB(ISO).EQ.1) THEN
         FMICR(IGR,ISO,I,J)=0.0
         HMICR(IGR,ISO,I,J)=0.0
         ENDIF
 874     CONTINUE
         DO 875 IL=1,NL
         DO 875 JGR=1,NGRP
         SCATF(J,I,IL,IGR,JGR)=0.0
         DO 875 ISO=2,NISO
         SMICR(J,ISO,I,IL,IGR,JGR)=0.0
 875     CONTINUE
 876     CONTINUE
         ENDIF
       ENDIF
 811   CONTINUE
C
C      ALL NOMINAL NEUTRONICS CONSTANTS ARE ALREDY STORED
C-----
C      STORING PROGRAM FOR THE FEEDBACK COEFFICIENTS.
C-----
         CALL LCMSIX(IPFBM,TEXTR,1)
         DO 391 I=1,NBURN
         WRITE(TEXTB,'(4HBURN,4X,I4)') I
         CALL LCMSIX(IPFBM,TEXTB,1)
         CALL LCMSIX(IPFBM,'MACR',1)
         CALL LCMSIX(IPFBM,'ABS',1)
         CALL LCMPUT(IPFBM,TEXT1,NGRP,2,TOTAF(1,I,1))
         CALL LCMSIX(IPFBM,' ',2)
C
      IF(IXYZ.EQ.0) THEN
       CALL LCMSIX(IPFBM,'STRD  ',1)
       CALL LCMPUT(IPFBM,TEXT1,NGRP,2,DXF(1,I,1))
       CALL LCMSIX(IPFBM,' ',2)
C
      ELSE IF(IXYZ.EQ.1) THEN
       CALL LCMSIX(IPFBM,'STRD X',1)
       CALL LCMPUT(IPFBM,TEXT1,NGRP,2,DXF(1,I,1))
       CALL LCMSIX(IPFBM,' ',2)
       CALL LCMSIX(IPFBM,'STRD Y',1)
       CALL LCMPUT(IPFBM,TEXT1,NGRP,2,DYF(1,I,1))
       CALL LCMSIX(IPFBM,' ',2)
       CALL LCMSIX(IPFBM,'STRD Z',1)
       CALL LCMPUT(IPFBM,TEXT1,NGRP,2,DZF(1,I,1))
       CALL LCMSIX(IPFBM,' ',2)
      ENDIF
C
      IF(JTAB(1).EQ.1) THEN
         CALL LCMSIX(IPFBM,'NUSIGF',1)
         CALL LCMPUT(IPFBM,TEXT1,NGRP,2,ZNUF(1,I,1))
         CALL LCMSIX(IPFBM,' ',2)
C
         CALL LCMSIX(IPFBM,'H-FACTORS',1)
         CALL LCMPUT(IPFBM,TEXT1,NGRP,2,HF(1,I,1))
         CALL LCMSIX(IPFBM,' ',2)
       ENDIF
C
         IF(NB.EQ.2) THEN
         CALL LCMSIX(IPFBM,'ABS',1)
         CALL LCMPUT(IPFBM,TEXT2,NGRP,2,TOTAF(1,I,2))
         CALL LCMSIX(IPFBM,' ',2)
C
         IF(IXYZ.EQ.0) THEN
          CALL LCMSIX(IPFBM,'STRD  ',1)
          CALL LCMPUT(IPFBM,TEXT2,NGRP,2,DXF(1,I,2))
          CALL LCMSIX(IPFBM,' ',2)
C
         ELSE IF(IXYZ.EQ.1) THEN
          CALL LCMSIX(IPFBM,'STRD X',1)
          CALL LCMPUT(IPFBM,TEXT2,NGRP,2,DXF(1,I,2))
          CALL LCMSIX(IPFBM,' ',2)
          CALL LCMSIX(IPFBM,'STRD Y',1)
          CALL LCMPUT(IPFBM,TEXT2,NGRP,2,DYF(1,I,2))
          CALL LCMSIX(IPFBM,' ',2)
          CALL LCMSIX(IPFBM,'STRD Z',1)
          CALL LCMPUT(IPFBM,TEXT2,NGRP,2,DZF(1,I,2))
          CALL LCMSIX(IPFBM,' ',2)
         ENDIF
C
          IF(JTAB(1).EQ.1) THEN
          CALL LCMSIX(IPFBM,'NUSIGF',1)
          CALL LCMPUT(IPFBM,TEXT2,NGRP,2,ZNUF(1,I,2))
          CALL LCMSIX(IPFBM,' ',2)
C
          CALL LCMSIX(IPFBM,'H-FACTORS',1)
          CALL LCMPUT(IPFBM,TEXT2,NGRP,2,HF(1,I,2))
          CALL LCMSIX(IPFBM,' ',2)
          ENDIF
         ENDIF
C
      IL=1
      WRITE (CM,'(I2)') IL-1
      CALL LCMSIX(IPFBM,'SCAT'//CM,1)
      CALL LCMLEN(IPFBM,'REF',ILENG,ITYXSM)
      IF(ILENG.GT.0) THEN
      IGAR=0
      DO 190 JGR=1,NGRP
      DO 190 IGR=MIJ(JGR),MIJ(JGR)-MNJ(JGR)+1,-1
      IGAR=IGAR+1
  190 WORK3(IGAR)=SCATF(1,I,IL,IGR,JGR)
      CALL LCMPUT(IPFBM,TEXT1,IGAR,2,WORK3)
       IF(NB.EQ.2) THEN
       IGAR=0
       DO 290 JGR=1,NGRP
       DO 290 IGR=MIJ(JGR),MIJ(JGR)-MNJ(JGR)+1,-1
       IGAR=IGAR+1
  290  WORK3(IGAR)=SCATF(2,I,IL,IGR,JGR)
       CALL LCMPUT(IPFBM,TEXT2,IGAR,2,WORK3)
       ENDIF
      ENDIF
      CALL LCMSIX(IPFBM,' ',2)
C
C
C    GO UP FOR MACR
        CALL LCMSIX(IPFBM,' ',2)
C-----
C      STORE MICROSCOPIC INFONFORMATION
C-----
      DO 401 ISO=2,NISO
      WRITE(HMICRO,'(3A4)') (HISO((ISO-1)*3+IH),IH=1,3)
      CALL LCMLEN(IPFBM,HMICRO,ILENG,ITYLCM)
      CALL LCMSIX(IPFBM,HMICRO,1)
C
      CALL LCMSIX(IPFBM,'ABS',1)
      CALL LCMPUT(IPFBM,TEXT1,NGRP,2,TMICR(1,ISO,I,1))
      CALL LCMSIX(IPFBM,' ',2)
C
      IF(IXYZ.EQ.0) THEN
       CALL LCMSIX(IPFBM,'STRD  ',1)
       CALL LCMPUT(IPFBM,TEXT1,NGRP,2,DMICRX(1,ISO,I,1))
       CALL LCMSIX(IPFBM,' ',2)
C
      ELSE IF(IXYZ.EQ.1) THEN
       CALL LCMSIX(IPFBM,'STRD X',1)
       CALL LCMPUT(IPFBM,TEXT1,NGRP,2,DMICRX(1,ISO,I,1))
       CALL LCMSIX(IPFBM,' ',2)
       CALL LCMSIX(IPFBM,'STRD Y',1)
       CALL LCMPUT(IPFBM,TEXT1,NGRP,2,DMICRY(1,ISO,I,1))
       CALL LCMSIX(IPFBM,' ',2)
       CALL LCMSIX(IPFBM,'STRD Z',1)
       CALL LCMPUT(IPFBM,TEXT1,NGRP,2,DMICRZ(1,ISO,I,1))
       CALL LCMSIX(IPFBM,' ',2)
      ENDIF
C
      IF(JTAB(ISO).EQ.1) THEN
      CALL LCMSIX(IPFBM,'NUSIGF',1)
      CALL LCMPUT(IPFBM,TEXT1,NGRP,2,FMICR(1,ISO,I,1))
      CALL LCMSIX(IPFBM,' ',2)
C
      CALL LCMSIX(IPFBM,'H-FACTORS',1)
      CALL LCMPUT(IPFBM,TEXT1,NGRP,2,HMICR(1,ISO,I,1))
      CALL LCMSIX(IPFBM,' ',2)
      ENDIF
C
       IF(NB.EQ.2) THEN
       CALL LCMSIX(IPFBM,'ABS',1)
       CALL LCMPUT(IPFBM,TEXT2,NGRP,2,TMICR(1,ISO,I,2))
       CALL LCMSIX(IPFBM,' ',2)
C
      IF(IXYZ.EQ.0) THEN
       CALL LCMSIX(IPFBM,'STRD  ',1)
       CALL LCMPUT(IPFBM,TEXT2,NGRP,2,DMICRX(1,ISO,I,2))
       CALL LCMSIX(IPFBM,' ',2)
C
      ELSE IF(IXYZ.EQ.1) THEN
       CALL LCMSIX(IPFBM,'STRD X',1)
       CALL LCMPUT(IPFBM,TEXT2,NGRP,2,DMICRX(1,ISO,I,2))
       CALL LCMSIX(IPFBM,' ',2)
       CALL LCMSIX(IPFBM,'STRD Y',1)
       CALL LCMPUT(IPFBM,TEXT2,NGRP,2,DMICRY(1,ISO,I,2))
       CALL LCMSIX(IPFBM,' ',2)
       CALL LCMSIX(IPFBM,'STRD Z',1)
       CALL LCMPUT(IPFBM,TEXT2,NGRP,2,DMICRZ(1,ISO,I,2))
       CALL LCMSIX(IPFBM,' ',2)
      ENDIF
C
       IF(JTAB(ISO).EQ.1) THEN
       CALL LCMSIX(IPFBM,'NUSIGF',1)
       CALL LCMPUT(IPFBM,TEXT2,NGRP,2,FMICR(1,ISO,I,2))
       CALL LCMSIX(IPFBM,' ',2)
C
       CALL LCMSIX(IPFBM,'H-FACTORS',1)
       CALL LCMPUT(IPFBM,TEXT2,NGRP,2,HMICR(1,ISO,I,2))
       CALL LCMSIX(IPFBM,' ',2)
       ENDIF
       ENDIF
C
      IL=1
      WRITE (CM,'(I2)') IL-1
      CALL LCMSIX(IPFBM,'SCAT'//CM,1)
      CALL LCMLEN(IPFBM,'REF',ILENG,ITYXSM)
      IF(ILENG.GT.0) THEN
      IGAR=0
      DO 191 JGR=1,NGRP
      DO 191 IGR=MIJ(JGR),MIJ(JGR)-MNJ(JGR)+1,-1
      IGAR=IGAR+1
  191 WORK3(IGAR)=SMICR(1,ISO,I,IL,IGR,JGR)
      CALL LCMPUT(IPFBM,TEXT1,IGAR,2,WORK3)
       IF(NB.EQ.2) THEN
       IGAR=0
       DO 291 JGR=1,NGRP
       DO 291 IGR=MIJ(JGR),MIJ(JGR)-MNJ(JGR)+1,-1
       IGAR=IGAR+1
  291  WORK3(IGAR)=SMICR(2,ISO,I,IL,IGR,JGR)
       CALL LCMPUT(IPFBM,TEXT2,IGAR,2,WORK3)
       ENDIF
      ENDIF
      CALL LCMSIX(IPFBM,' ',2)
C
C     GO UP FOR MICR
        CALL LCMSIX(IPFBM,' ',2)
  401  CONTINUE
C     GO UP FOR BURN
        CALL LCMSIX(IPFBM,' ',2)
 391    CONTINUE
C
        CALL LCMSIX(IPFBM,' ',2)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(IXS,NJJ,IJJ)
C
      RETURN
6000  FORMAT(' Keyword =',A8,' ncor =',i4,' Param =',8(2X,A8))
      END