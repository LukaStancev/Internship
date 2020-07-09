*DECK OUTDRV
      SUBROUTINE OUTDRV (IPGEOM,IPMAC1,IPFLUX,IPMAC2,MAXNEL,NBMIX,NL,
     1 NBFIS,NGRP,NEL,NUN,NALBP,HTRACK,IELEM,ICOL,MAT,VOL,IDL,TITR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* driver for the port-treatment of reactor calculation results.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A Hebert
*
*Parameters: input
* IPGEOM  L_GEOM pointer to the geometry.
* IPMAC1  L_MACROLIB pointer to the nuclear properties.
* IPFLUX  L_FLUX pointer to the solution.
* IPMAC2  L_MACROLIB pointer to the edition information.
* MAXNEL  maximum number of finite elements.
* NBMIX   number of material mixtures.
* NL      scattering anisotropy.
* NBFIS   number of fissionable isotopes.
* NGRP    total number of energy groups.
* NEL     total number of finite elements.
* NUN     total number of unknowns per group.
* NALBP   number of physical albedos.
* HTRACK  type of tracking (equal to 'BIVAC' or 'TRIVAC').
* IELEM   degree of the Lagrangian finite elements:
* ICOL    type of quadrature used to integrate the mass matrix
* MAT     index-number of the mixture type assigned to each volume.
* VOL     volumes.
* IDL     position of the average flux component associated with
*         each volume.
* TITR    title.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPGEOM,IPMAC1,IPMAC2,IPFLUX
      CHARACTER TITR*72,HTRACK*12
      INTEGER MAXNEL,NBMIX,NL,NBFIS,NGRP,NEL,NUN,NALBP,IELEM,ICOL,
     1 MAT(NEL),IDL(NEL)
      REAL VOL(NEL)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPMAC1,KPMAC1
      CHARACTER TEXT4*4
      DOUBLE PRECISION DFLOTT,ZNORM
      INTEGER, DIMENSION(:), ALLOCATABLE :: IHOM
      REAL, DIMENSION(:), ALLOCATABLE :: SGD
      REAL, DIMENSION(:,:), ALLOCATABLE :: EVECT,ADECT,ZUFIS
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IHOM(NEL),EVECT(NUN,NGRP),SGD(NBMIX))
*
      TKR=0.0
      IMPX=1
      IADJ=0
      LMOD=0
      CALL KDRCPU(TK1)
*----
*  RECOVER THE K-EFFECTIVE.
*----
      CALL LCMLEN(IPFLUX,'K-EFFECTIVE',ILEN,ITYLCM)
      IF(ILEN.GT.0) THEN
         CALL LCMGET(IPFLUX,'K-EFFECTIVE',FKEFF)
         CALL LCMPUT(IPMAC2,'K-EFFECTIVE',1,2,FKEFF)
      ENDIF
*
   30 CALL REDGET (INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
      IF(INDIC.NE.3) CALL XABORT('OUTDRV: CHARACTER DATA EXPECTED.')
*
      IF(TEXT4.EQ.'POWR') THEN
*        NORMALIZATION TO A GIVEN FISSION POWER.
         CALL REDGET (INDIC,NITMA,POWER,TEXT4,DFLOTT)
         IF(INDIC.NE.2) CALL XABORT('OUTDRV: REAL DATA EXPECTED.')
*        NORMALIZATION FACTOR FOR THE DIRECT FLUX.
         CALL OUTFLX(IPFLUX,0,NGRP,NUN,LMOD,IMPX,EVECT)
         ZNORM=0.0D0
         JPMAC1=LCMGID(IPMAC1,'GROUP')
         DO 45 IGR=1,NGRP
         KPMAC1=LCMGIL(JPMAC1,IGR)
         CALL LCMLEN(KPMAC1,'H-FACTOR',LENGT,ITYLCM)
         IF(LENGT.GT.0) THEN
            CALL LCMGET(KPMAC1,'H-FACTOR',SGD)
         ELSE
            WRITE(6,'(/43H OUTDRV: *** WARNING *** NO H-FACTOR FOUND ,
     1      28HON LCM. USE NU*SIGF INSTEAD.)')
            ALLOCATE(ZUFIS(NBMIX,NBFIS))
            CALL XDRSET(SGD,NBMIX,0.0)
            CALL LCMGET(KPMAC1,'NUSIGF',ZUFIS)
            DO IBM=1,NBMIX
              DO IFISS=1,NBFIS
                SGD(IBM)=SGD(IBM)+ZUFIS(IBM,IFISS)
              ENDDO
            ENDDO
            DEALLOCATE(ZUFIS)
         ENDIF
         DO 40 K=1,NEL
         L=MAT(K)
         IF((L.EQ.0).OR.(IDL(K).EQ.0)) GO TO 40
         ZNORM=ZNORM+EVECT(IDL(K),IGR)*VOL(K)*SGD(L)
   40    CONTINUE
   45    CONTINUE
         ZNORM=POWER/ZNORM
         WRITE(6,300) ' DIRECT',ZNORM
         DO 55 IGR=1,NGRP
         DO 50 I=1,NUN
         EVECT(I,IGR)=EVECT(I,IGR)*REAL(ZNORM)
   50    CONTINUE
   55    CONTINUE
      ELSE IF(TEXT4.EQ.'EDIT') THEN
         CALL REDGET(INDIC,IMPX,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('OUT: INTEGER DATA EXPECTED.')
      ELSE IF(TEXT4.EQ.'MODE') THEN
         CALL REDGET(INDIC,LMOD,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('OUT: INTEGER DATA EXPECTED.')
      ELSE IF(TEXT4.EQ.'DIRE') THEN
         IADJ=0
      ELSE IF(TEXT4.EQ.'PROD') THEN
         IADJ=1
      ELSE IF(TEXT4.EQ.'INTG') THEN
*        COMPUTE AND DISPLAY THE MACRO-ZONE REACTION RATES.
*        READ THE MACRO-ZONES DEFINITION.
         CALL OUTHOM (MAXNEL,IPGEOM,IMPX,NEL,IELEM,ICOL,HTRACK,MAT,NZS,
     1   IHOM)
         IF(NZS.GT.NEL) CALL XABORT('OUTDRV: INVALID VALUE OF NZS.')
         MAXALC=MAX(NBMIX,NZS)
         CALL OUTFLX(IPFLUX,0,NGRP,NUN,LMOD,IMPX,EVECT)
         IF(IMPX.GT.0) WRITE(6,320) TITR
         IF(IADJ.EQ.0) THEN
           CALL OUTAUX(IPMAC1,IPMAC2,NBMIX,NL,NBFIS,NGRP,NEL,NUN,NALBP,
     1     NZS,MAT,VOL,IDL,EVECT,IHOM,IMPX,MAXALC)
         ELSE IF(IADJ.EQ.1) THEN
           ALLOCATE(ADECT(NUN,NGRP))
           CALL OUTFLX(IPFLUX,1,NGRP,NUN,LMOD,IMPX,ADECT)
           CALL OUTPRO(IPMAC1,IPMAC2,NBMIX,NL,NBFIS,NGRP,NEL,NUN,NALBP,
     1     NZS,MAT,VOL,IDL,EVECT,ADECT,IHOM,IMPX,MAXALC)
           DEALLOCATE(ADECT)
         ENDIF
      ELSE IF(TEXT4.EQ.';') THEN
         CALL KDRCPU(TK2)
         TKR=TK2-TK1
         WRITE(6,310) TKR
         GO TO 60
      ELSE
         CALL XABORT('OUTDRV: '//TEXT4//' IS AN INVALID KEY WORD.')
      ENDIF
      GO TO 30
*----
*  SCRATCH STORAGE DEALLOCATION
*----
   60 DEALLOCATE(IHOM,EVECT,SGD)
      RETURN
*
  300 FORMAT(/9H OUTDRV: ,A7,28H FLUX NORMALIZATION FACTOR =,1P,E13.5)
  310 FORMAT(/49H OUTDRV: CPU TIME FOR REACTION RATE CALCULATION =,F7.3)
  320 FORMAT(/12H OUTDRV: ***,A72,3H***)
      END
