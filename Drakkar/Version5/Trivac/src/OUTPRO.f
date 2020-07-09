*DECK OUTPRO
      SUBROUTINE OUTPRO (IPMAC1,IPMAC2,NBMIX,NL,NBFIS,NGRP,NEL,NUN,
     1 NALBP,NZS,MAT,VOL,IDL,EVECT,ADECT,IHOM,IMPX,MAXALC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* perform direct-adjoint homogenization into NZS regions based on
* averaged fluxes contained in EVECT and adjoint fluxes contained in
* ADECT. Create an output extended macrolib containing homogenized
* volumes, integrated fluxes and cross sections.
*
*Copyright:
* Copyright (C) 2018 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPMAC1  L_MACROLIB pointer to the input macrolib.
* IPMAC2  L_MACROLIB pointer to the output extended macrolib.
* NBMIX   number of material mixtures.
* NL      scattering anisotropy.
* NBFIS   number of fissionable isotopes.
* NGRP    total number of energy groups.
* NEL     number of finite elements.
* NUN     number of unknowns per energy group.
* NALBP   number of physical albedos.
* NZS     number of homogenized regions so that NZS=max(IHOM(i)).
* MAT     index-number of the mixture type assigned to each volume.
* VOL     volumes.
* IDL     position of the average flux component associated with
*         each volume.
* EVECT   flux unknowns.
* ADECT   adjoint flux unknowns.
* IHOM    homogenized index assigned to each element.
* IMPX    print parameter (equal to zero for no print).
* MAXALC  maximum of NBMIX and NZS.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAC1,IPMAC2
      PARAMETER(NREAC=10)
      INTEGER NBMIX,NL,NBFIS,NGRP,NEL,NUN,NALBP,NZS,MAT(NEL),IDL(NEL),
     1 IHOM(NEL),IMPX,MAXALC
      REAL VOL(NEL),EVECT(NUN,NGRP),ADECT(NUN,NGRP)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPMAC1,KPMAC1,JPMAC2,KPMAC2
      PARAMETER(NSTATE=40)
      CHARACTER HREAC(NREAC)*12,TEXT12*12,SUFF*2
      INTEGER IDATA(NSTATE)
      LOGICAL LNUSIG
      INTEGER, DIMENSION(:), ALLOCATABLE :: IJJ,NJJ,IPOS
      REAL, DIMENSION(:), ALLOCATABLE :: WORK,RATE,GAR
      REAL, DIMENSION(:,:), ALLOCATABLE :: FLINT,AFLINT,SCAT,CHI,ZUFIS,
     1 OUTR,ALBP
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ACCUM
      DATA HREAC/'NTOT0','SIGW00','NUSIGF','NFTOT','H-FACTOR',
     1 'OVERV','DIFF','DIFFX','DIFFY','DIFFZ'/
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IJJ(MAXALC),NJJ(MAXALC),IPOS(MAXALC))
      ALLOCATE(WORK(MAXALC*NGRP),RATE(NZS),FLINT(NZS,NGRP),
     1 AFLINT(NZS,NGRP),GAR(NGRP),SCAT(NZS,NGRP),CHI(NBMIX,NBFIS),
     2 ZUFIS(NBMIX,NBFIS),OUTR(NZS+1,NREAC+3))
      ALLOCATE(ACCUM(NZS,NBFIS))
*
      CALL XDISET(IDATA,NSTATE,0)
      PRINT *,'process with OUTPRO'
*----
*  DIRECT FLUX CALCULATION.
*----
      CALL XDRSET(WORK,NZS,0.0)
      CALL XDRSET(FLINT,NZS*NGRP,0.0)
      DO 20 K=1,NEL
      IBM=IHOM(K)
      IPFL=IDL(K)
      IF((IBM.NE.0).AND.(MAT(K).NE.0).AND.(IPFL.NE.0)) THEN
         WORK(IBM)=WORK(IBM)+VOL(K)
         DO 10 IGR=1,NGRP
         FLINT(IBM,IGR)=FLINT(IBM,IGR)+EVECT(IPFL,IGR)*VOL(K)
   10    CONTINUE
      ENDIF
   20 CONTINUE
      CALL LCMPUT(IPMAC2,'VOLUME',NZS,2,WORK)
      DO 25 IBM=1,NZS
      OUTR(IBM,NREAC+1)=WORK(IBM)
   25 CONTINUE
*----
*  ADJOINT FLUX CALCULATION.
*----
      CALL XDRSET(AFLINT,NZS*NGRP,0.0)
      DO 35 K=1,NEL
      IBM=IHOM(K)
      IPFL=IDL(K)
      IF((IBM.NE.0).AND.(MAT(K).NE.0).AND.(IPFL.NE.0)) THEN
        DO 30 IGR=1,NGRP
        AFLINT(IBM,IGR)=AFLINT(IBM,IGR)+ADECT(IPFL,IGR)*VOL(K)
   30   CONTINUE
      ENDIF
   35 CONTINUE
      DO 45 IGR=1,NGRP
      DO 40 IBM=1,NZS
      AFLINT(IBM,IGR)=AFLINT(IBM,IGR)/WORK(IBM)
   40 CONTINUE
   45 CONTINUE
*----
*  FISSION RATE CALCULATION.
*----
      IF(IMPX.GT.0) WRITE(6,'(/35H OUTPRO: REACTION RATE CALCULATION.)')
      JPMAC1=LCMGID(IPMAC1,'GROUP')
      JPMAC2=LCMLID(IPMAC2,'GROUP',NGRP)
      IF(NBFIS.GT.0) THEN
         DO 55 IFISS=1,NBFIS
         DO 50 IBM=1,NZS
         ACCUM(IBM,IFISS)=0.0D0
   50    CONTINUE
   55    CONTINUE
         DO 62 IGR=1,NGRP
         KPMAC1=LCMGIL(JPMAC1,IGR)
         CALL LCMGET(KPMAC1,'NUSIGF',ZUFIS)
         DO 61 IFISS=1,NBFIS
         DO 60 K=1,NEL
         IBM=IHOM(K)
         L=MAT(K)
         IPFL=IDL(K)
         IF((IBM.NE.0).AND.(L.NE.0).AND.(IPFL.NE.0)) THEN
            ACCUM(IBM,IFISS)=ACCUM(IBM,IFISS)+EVECT(IPFL,IGR)*VOL(K)*
     1      ZUFIS(L,IFISS)
         ENDIF
   60    CONTINUE
   61    CONTINUE
   62    CONTINUE
      ENDIF
*----
*  LOOP OVER ENERGY GROUP LIST.
*----
      DO 310 IGR=1,NGRP
      KPMAC1=LCMGIL(JPMAC1,IGR)
      KPMAC2=LCMDIL(JPMAC2,IGR)
      CALL LCMPUT(KPMAC2,'FLUX-INTG',NZS,2,FLINT(1,IGR))
      CALL LCMPUT(KPMAC2,'NWAT0',NZS,2,AFLINT(1,IGR))
      DO 70 IBM=1,NZS
      OUTR(IBM,NREAC+2)=FLINT(IBM,IGR)
   70 CONTINUE
*----
*  CROSS SECTIONS HOMOGENIZATION.
*----
      LNUSIG=.FALSE.
      DO 130 IREAC=1,NREAC
      DO 80 IBM=1,NZS
      OUTR(IBM,IREAC)=0.0
   80 CONTINUE
      CALL LCMLEN(KPMAC1,HREAC(IREAC),LENGT,ITYLCM)
      IF((HREAC(IREAC).EQ.'H-FACTOR').AND.(LENGT.EQ.0)) THEN
         WRITE(6,'(/46H OUTPRO: *** WARNING *** NO H-FACTOR FOUND ON ,
     1   25HLCM. USE NU*SIGF INSTEAD.)')
         LNUSIG=.TRUE.
         GO TO 130
      ELSE IF(HREAC(IREAC).EQ.'NUSIGF') THEN
         GO TO 130
      ELSE IF(HREAC(IREAC).EQ.'SIGW00') THEN
         GO TO 130
      ELSE
         TEXT12=HREAC(IREAC)
      ENDIF
      IF(LENGT.GT.0) THEN
         IF(LENGT.GT.NBMIX) CALL XABORT('OUTPRO: INVALID LENGTH FOR '//
     1   HREAC(IREAC)//' CROSS SECTIONS.')
         CALL LCMGET(KPMAC1,TEXT12,WORK)
         DO 90 IBM=1,NZS
         RATE(IBM)=0.0
   90    CONTINUE
         DO 100 K=1,NEL
         IBM=IHOM(K)
         L=MAT(K)
         IPFL=IDL(K)
         IF((IBM.NE.0).AND.(L.NE.0).AND.(IPFL.NE.0)) THEN
            RATE(IBM)=RATE(IBM)+ADECT(IPFL,IGR)*EVECT(IPFL,IGR)*VOL(K)*
     1      WORK(L)
         ENDIF
  100    CONTINUE
         DO 120 IBM=1,NZS
         OUTR(IBM,IREAC)=RATE(IBM)
         IF(RATE(IBM).NE.0.0) RATE(IBM)=RATE(IBM)/(AFLINT(IBM,IGR)*
     1   FLINT(IBM,IGR))
  120    CONTINUE
         CALL LCMPUT(KPMAC2,HREAC(IREAC),NZS,2,RATE)
         IF(HREAC(IREAC).EQ.'DIFF') THEN
            IDATA(9)=1
         ELSE IF(HREAC(IREAC).EQ.'DIFFX') THEN
            IDATA(9)=2
         ENDIF
      ENDIF
  130 CONTINUE
*----
*  SCATTERING MATRIX HOMOGENIZATION.
*----
      DO 215 IL=1,NL
      WRITE(SUFF,'(I2.2)') IL-1
      CALL LCMLEN(KPMAC1,'NJJS'//SUFF,LENGT,ITYLCM)
      IF(LENGT.GT.0) THEN
         IF(LENGT.GT.NBMIX) CALL XABORT('OUTPRO: INVALID LENGTH FOR '//
     1   'SCATTERING CROSS SECTIONS.')
         DO 145 JGR=1,NGRP
         DO 140 IBM=1,NZS
        SCAT(IBM,JGR)=0.0
  140    CONTINUE
  145    CONTINUE
         CALL LCMGET(KPMAC1,'NJJS'//SUFF,NJJ)
         CALL LCMGET(KPMAC1,'IJJS'//SUFF,IJJ)
         CALL LCMGET(KPMAC1,'IPOS'//SUFF,IPOS)
         CALL LCMGET(KPMAC1,'SCAT'//SUFF,WORK)
         IPOSDE=0
         DO 180 K=1,NEL
         IBM=IHOM(K)
         L=MAT(K)
         IPFL=IDL(K)
         IF((IBM.NE.0).AND.(L.NE.0).AND.(IPFL.NE.0)) THEN
            DO 150 JGR=1,NGRP
            GAR(JGR)=0.0
  150       CONTINUE
            IPOSDE=IPOS(L)-1
            DO 160 JGR=IJJ(L),IJJ(L)-NJJ(L)+1,-1
            IPOSDE=IPOSDE+1
            GAR(JGR)=WORK(IPOSDE)
  160       CONTINUE
            DO 170 JGR=1,NGRP
            SCAT(IBM,JGR)=SCAT(IBM,JGR)+ADECT(IPFL,IGR)*EVECT(IPFL,JGR)*
     1      VOL(K)*GAR(JGR)
  170       CONTINUE
         ENDIF
  180    CONTINUE
         DO 200 IBM=1,NZS
         IGMIN=IGR
         IGMAX=IGR
         DO 190 JGR=NGRP,1,-1
         IF(SCAT(IBM,JGR).NE.0.0) THEN
            IGMIN=MIN(IGMIN,JGR)
            IGMAX=MAX(IGMAX,JGR)
            SCAT(IBM,JGR)=SCAT(IBM,JGR)/(AFLINT(IBM,IGR)*FLINT(IBM,JGR))
         ENDIF
  190    CONTINUE
         IJJ(IBM)=IGMAX
         NJJ(IBM)=IGMAX-IGMIN+1
  200    CONTINUE
         IPOSDE=0
         DO 211 IBM=1,NZS
         IPOS(IBM)=IPOSDE+1
         DO 210 JGR=IJJ(IBM),IJJ(IBM)-NJJ(IBM)+1,-1
         IPOSDE=IPOSDE+1
         WORK(IPOSDE)=SCAT(IBM,JGR)
  210    CONTINUE
  211    CONTINUE
         CALL LCMPUT(KPMAC2,'SCAT'//SUFF,IPOSDE,2,WORK)
         CALL LCMPUT(KPMAC2,'IPOS'//SUFF,NZS,1,IPOS)
         CALL LCMPUT(KPMAC2,'NJJS'//SUFF,NZS,1,NJJ)
         CALL LCMPUT(KPMAC2,'IJJS'//SUFF,NZS,1,IJJ)
         CALL LCMPUT(KPMAC2,'SIGW'//SUFF,NZS,2,SCAT(1,IGR))
         IF(IL.EQ.1) OUTR(:NZS,2)=SCAT(:NZS,IGR)*OUTR(:NZS,NREAC+2)
      ENDIF
  215 CONTINUE
*----
*  FISSION SPECTRUM AND NUSIGF HOMOGENIZATION.
*----
      CALL XDRSET(OUTR(1,NREAC+3),NZS,0.0)
      IF(NBFIS.GT.0) THEN
         CALL LCMLEN(KPMAC1,'NUSIGF',LENGT,ITYLCM)
         IF(LENGT.NE.NBMIX*NBFIS) CALL XABORT('OUTPRO: INVALID LENGTH '
     1   //'FOR FISSION SPECTRUM.')
         CALL LCMGET(KPMAC1,'NUSIGF',ZUFIS)
         CALL XDRSET(RATE,NZS,0.0)
         CALL LCMLEN(KPMAC1,'CHI',LENGT,ITYLCM)
         IF(LENGT.EQ.0) THEN
            IF(IGR.EQ.1) THEN
               CALL XDRSET(RATE,NZS,1.0)
               CALL XDRSET(OUTR(1,NREAC+3),NZS,1.0)
            ENDIF
         ELSE
            CALL LCMGET(KPMAC1,'CHI',CHI)
            CALL XDRSET(SCAT(1,1),NZS,0.0)
            DO 230 K=1,NEL
            IBM=IHOM(K)
            L=MAT(K)
            IPFL=IDL(K)
            IF((IBM.NE.0).AND.(L.NE.0).AND.(IPFL.NE.0)) THEN
               DO 220 IFISS=1,NBFIS
               RATE(IBM)=RATE(IBM)+CHI(L,IFISS)*ADECT(IPFL,IGR)*
     1         REAL(ACCUM(IBM,IFISS))
               SCAT(IBM,1)=SCAT(IBM,1)+REAL(ACCUM(IBM,IFISS))
  220          CONTINUE
            ENDIF
  230       CONTINUE
            DO 240 IBM=1,NZS
            IF(SCAT(IBM,1).NE.0.0) THEN
               RATE(IBM)=RATE(IBM)/(SCAT(IBM,1)*AFLINT(IBM,IGR))
               OUTR(IBM,NREAC+3)=RATE(IBM)
            ENDIF
  240       CONTINUE
         ENDIF
         CALL LCMPUT(KPMAC2,'CHI',NZS,2,RATE)
         DO 250 IBM=1,NZS
         RATE(IBM)=0.0
  250    CONTINUE
         DO 265 IFISS=1,NBFIS
         DO 260 K=1,NEL
         IBM=IHOM(K)
         L=MAT(K)
         IPFL=IDL(K)
         IF((IBM.NE.0).AND.(L.NE.0).AND.(IPFL.NE.0)) THEN
            RATE(IBM)=RATE(IBM)+EVECT(IPFL,IGR)*VOL(K)*ZUFIS(L,IFISS)
         ENDIF
  260    CONTINUE
  265    CONTINUE
         DO 270 IBM=1,NZS
         OUTR(IBM,3)=RATE(IBM)
         IF(RATE(IBM).NE.0.0) RATE(IBM)=RATE(IBM)/FLINT(IBM,IGR)
  270    CONTINUE
         CALL LCMPUT(KPMAC2,'NUSIGF',NZS,2,RATE)
         IF(LNUSIG) CALL LCMPUT(KPMAC2,'H-FACTOR',NZS,2,RATE)
      ENDIF
*----
*  PHYSICAL ALBEDOS.
*----
      IF(NALBP.GT.0) THEN
         ALLOCATE(ALBP(NALBP,NGRP))
         CALL LCMGET(IPMAC1,'ALBEDO',ALBP)
         CALL LCMPUT(IPMAC2,'ALBEDO',NALBP*NGRP,2,ALBP)
         DEALLOCATE(ALBP)
      ENDIF
*----
*  PRINT THE REACTION RATES:
*----
      IF(IMPX.GT.0) THEN
         DO 280 I=1,NREAC+2
         OUTR(NZS+1,I)=0.0
  280    CONTINUE
         WRITE(6,500) IGR,'VOLUME     ','FLUX-INTG  ',(HREAC(I),I=1,6),
     1   'CHI        '
         DO 300 IBM=1,NZS
         DO 290 I=1,NREAC+2
         OUTR(NZS+1,I)=OUTR(NZS+1,I)+OUTR(IBM,I)
  290    CONTINUE
         WRITE(6,510) IBM,OUTR(IBM,NREAC+1),OUTR(IBM,NREAC+2),
     1   (OUTR(IBM,I),I=1,6),OUTR(IBM,NREAC+3)
  300    CONTINUE
         WRITE(6,520) OUTR(NZS+1,NREAC+1),OUTR(NZS+1,NREAC+2),
     1   (OUTR(NZS+1,I),I=1,6)
      ENDIF
  310 CONTINUE
      IDATA(1)=NGRP
      IDATA(2)=NZS
      IDATA(3)=NL
      IDATA(4)=1
      IDATA(8)=NALBP
      IDATA(15)=1
      CALL LCMPUT(IPMAC2,'STATE-VECTOR',NSTATE,1,IDATA)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(ACCUM)
      DEALLOCATE(WORK,RATE,AFLINT,FLINT,GAR,SCAT,CHI,ZUFIS,OUTR)
      DEALLOCATE(IJJ,NJJ,IPOS)
      RETURN
*
  500 FORMAT(/' G R O U P   : ',I3/1X,'IHOM',9A14)
  510 FORMAT(1X,I4,1P,9E14.5)
  520 FORMAT(/5H  SUM,1P,8E14.5)
      END
