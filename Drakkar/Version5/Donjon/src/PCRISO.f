*DECK PCRISO
      SUBROUTINE PCRISO(IPLIB,KPTMP,HNAME,JSO,NCAL,NGRP,NL,NED,HVECT,
     1 NDEL,IMPX,TERP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* recover nuclear data from a single isotopic directory.
*
*Copyright:
* Copyright (C) 2019 Ecole Polytechnique de Montreal
*
*Author(s): A. Hebert
*
*Parameters: input
* IPLIB   address of the microlib LCM object.
* KPTMP   address of the 'CALCULATIONS' list.
* HNAME   character*12 name of the PMAXS isotope been processed.
* JSO     index of the PMAXS isotope been processed.
* NCAL    number of elementary calculations in the PMAXS file.
* NGRP    number of energy groups.
* NL      number of Legendre orders.
* NED     number of extra vector edits.
* HVECT   character names of the extra vector edits.
* NDEL    number of delayed precursor groups.
* IMPX    print parameter (equal to zero for no print).
* TERP    interpolation weights.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB,KPTMP
      INTEGER JSO,NCAL,NGRP,NL,NED,NDEL,IMPX
      REAL TERP(NCAL)
      CHARACTER HNAME*12,HVECT(NED)*(*)
*----
*  LOCAL VARIABLES
*----
      INTEGER, PARAMETER::IOUT=6
      REAL TAUXFI, TAUXF, WEIGHT
      INTEGER ICAL, IDEL, IED, IG1, IG2, IG, ILENG, IL, ITYLCM, J,
     & LENGTH, MAXH
      LOGICAL LWD
      CHARACTER CM*2,HMAKE(100)*12,TEXT12*12
      TYPE(C_PTR) LPTMP,MPTMP,NPTMP
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ITYPR
      REAL, ALLOCATABLE, DIMENSION(:) :: WDLA
      REAL, ALLOCATABLE, DIMENSION(:,:) :: GAR1,GAR2
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: WSCA1,WSCA2
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(ITYPR(NL))
      ALLOCATE(GAR1(NGRP,10+NL+NED+2*NDEL),WSCA1(NGRP,NGRP,NL),
     1 GAR2(NGRP,10+NL+NED+2*NDEL),WSCA2(NGRP,NGRP,NL),WDLA(NDEL))
*----
*  RECOVER GENERIC ISOTOPIC DATA FROM THE PMAXS FILE
*----
      LWD=.FALSE.
      DO 10 ICAL=1,NCAL
      WEIGHT=TERP(ICAL)
      IF(WEIGHT.EQ.0.0) GO TO 10
      LPTMP=LCMGIL(KPTMP,ICAL)
      CALL LCMLEN(LPTMP,'ISOTOPESLIST',LENGTH,ITYLCM)
      IF(LENGTH.EQ.0) GO TO 10
      MPTMP=LCMGID(LPTMP,'ISOTOPESLIST')
      CALL LCMLEL(MPTMP,JSO,ILENG,ITYLCM)
      IF(ILENG.EQ.0) GO TO 10
      NPTMP=LCMGIL(MPTMP,JSO)
      CALL LCMGTC(NPTMP,'ALIAS',12,1,TEXT12)
      IF(TEXT12(:8).NE.HNAME(:8)) GO TO 10
      CALL LCMLEN(NPTMP,'LAMBDA-D',LENGTH,ITYLCM)
      LWD=(LENGTH.EQ.NDEL).AND.(NDEL.GT.0)
      IF(LWD) CALL LCMGET(NPTMP,'LAMBDA-D',WDLA)
      GO TO 15
   10 CONTINUE
      CALL XABORT('PCRISO: UNABLE TO FIND A DIRECTORY FOR ISOTOPE '//
     1 HNAME//'.')
*----
*  LOOP OVER ELEMENTARY CALCULATIONS
*----
   15 MAXH=10+NL+NED+2*NDEL
      IF(MAXH+NL.GT.100) CALL XABORT('PCRISO: STATIC STORAGE EXCEEDED')
      DO J=1,MAXH+NL
         HMAKE(J)=' '
      ENDDO
      GAR2(:NGRP,:MAXH)=0.0
      WSCA2(:NGRP,:NGRP,:NL)=0.0
      TAUXFI=0.0
      DO 120 ICAL=1,NCAL
      WEIGHT=TERP(ICAL)
      IF(WEIGHT.EQ.0.0) GO TO 120
      LPTMP=LCMGIL(KPTMP,ICAL)
      IF(IMPX.GT.4) THEN
         WRITE(IOUT,'(34H PCRISO: PMAXS ACCESS FOR ISOTOPE ,A,6H AND C,
     1   10HALCULATION,I5,1H.)') HNAME,ICAL
         IF(IMPX.GT.50) CALL LCMLIB(LPTMP)
      ENDIF
      MPTMP=LCMGID(LPTMP,'ISOTOPESLIST')
      CALL LCMLEL(MPTMP,JSO,ILENG,ITYLCM)
      IF(ILENG.EQ.0) GO TO 120
      NPTMP=LCMGIL(MPTMP,JSO)
*----
*  RECOVER CALCULATION-SPECIFIC ISOTOPIC DATA FROM THE PMAXS FILE
*----
      CALL LCMLEN(NPTMP,'NWT0',LENGTH,ITYLCM)
      IF(LENGTH.EQ.NGRP) THEN
         CALL LCMGET(NPTMP,'NWT0',GAR1(1,1))
         HMAKE(1)='NWT0'
      ENDIF
      CALL LCMLEN(NPTMP,'NWT1',LENGTH,ITYLCM)
      IF(LENGTH.EQ.NGRP) THEN
         CALL LCMGET(NPTMP,'NWT1',GAR1(1,2))
         HMAKE(2)='NWT1'
      ENDIF
      CALL XDRLGS(NPTMP,-1,IMPX,0,NL-1,1,NGRP,GAR1(1,3),WSCA1,ITYPR)
      DO IL=0,NL-1
         IF(ITYPR(IL+1).NE.0) THEN
            WRITE (CM,'(I2.2)') IL
            HMAKE(3+IL)='SIGS'//CM
         ENDIF
      ENDDO
      CALL LCMGET(NPTMP,'NTOT0',GAR1(1,3+NL))
      HMAKE(3+NL)='NTOT0'
      CALL LCMLEN(NPTMP,'NTOT1',LENGTH,ITYLCM)
      IF(LENGTH.EQ.NGRP) THEN
         CALL LCMGET(NPTMP,'NTOT1',GAR1(1,4+NL))
         HMAKE(4+NL)='NTOT1'
      ENDIF
      CALL LCMLEN(NPTMP,'NUSIGF',LENGTH,ITYLCM)
      IF(LENGTH.EQ.NGRP) THEN
         CALL LCMGET(NPTMP,'NUSIGF',GAR1(1,5+NL))
         HMAKE(5+NL)='NUSIGF'
         CALL LCMGET(NPTMP,'CHI',GAR1(1,MAXH-NDEL-1))
         HMAKE(MAXH-NDEL-1)='CHI'
      ENDIF
      IF(NDEL.GT.0) THEN
         WRITE(TEXT12,'(6HNUSIGF,I2.2)') NDEL
         CALL LCMLEN(NPTMP,TEXT12,LENGTH,ITYLCM)
         IF(LENGTH.EQ.NGRP) THEN
            DO IDEL=1,NDEL
               WRITE(TEXT12,'(6HNUSIGF,I2.2)') IDEL
               CALL LCMGET(NPTMP,TEXT12,GAR1(1,MAXH-2*NDEL-2+IDEL))
               HMAKE(MAXH-2*NDEL-2+IDEL)=TEXT12
            ENDDO
         ENDIF
         WRITE(TEXT12,'(3HCHI,I2.2)') NDEL
         CALL LCMLEN(NPTMP,TEXT12,LENGTH,ITYLCM)
         IF(LENGTH.EQ.NGRP) THEN
            DO IDEL=1,NDEL
               WRITE(TEXT12,'(3HCHI,I2.2)') IDEL
               CALL LCMGET(NPTMP,TEXT12,GAR1(1,MAXH-NDEL-1+IDEL))
               HMAKE(MAXH-NDEL-1+IDEL)=TEXT12
            ENDDO
         ENDIF
      ENDIF
      CALL LCMLEN(NPTMP,'H-FACTOR',LENGTH,ITYLCM)
      IF(LENGTH.EQ.NGRP) THEN
         CALL LCMGET(NPTMP,'H-FACTOR',GAR1(1,MAXH-2*NDEL-4))
         HMAKE(MAXH-2*NDEL-4)='H-FACTOR'
      ENDIF
      CALL LCMLEN(NPTMP,'OVERV',LENGTH,ITYLCM)
      IF(LENGTH.EQ.NGRP) THEN
         CALL LCMGET(NPTMP,'OVERV',GAR1(1,MAXH-2*NDEL-3))
         HMAKE(MAXH-2*NDEL-3)='OVERV'
      ENDIF
      CALL LCMLEN(NPTMP,'TRANC',LENGTH,ITYLCM)
      IF(LENGTH.EQ.NGRP) THEN
         CALL LCMGET(NPTMP,'TRANC',GAR1(1,MAXH-2*NDEL-2))
         HMAKE(MAXH-2*NDEL-2)='TRANC'
      ENDIF
      DO IED=1,NED
         CALL LCMLEN(NPTMP,HVECT(IED),LENGTH,ITYLCM)
         IF((LENGTH.GT.0).AND.(HVECT(IED).NE.'TRANC')) THEN
            CALL LCMGET(NPTMP,HVECT(IED),GAR1(1,5+NL+IED))
            HMAKE(5+NL+IED)=HVECT(IED)
         ENDIF
      ENDDO
      CALL LCMLEN(NPTMP,'STRD',LENGTH,ITYLCM)
      IF(LENGTH.EQ.NGRP) THEN
         CALL LCMGET(NPTMP,'STRD',GAR1(1,MAXH))
         HMAKE(MAXH)='STRD'
      ENDIF
*----
*  COMPUTE FISSION RATE FOR A SINGLE ELEMENTARY CALCULATION
*----
      TAUXF=0.0
      IF(HMAKE(5+NL).EQ.'NUSIGF') THEN
        DO IG=1,NGRP
           TAUXF=TAUXF+GAR1(IG,5+NL)*GAR1(IG,1)
        ENDDO
        TAUXFI=TAUXFI+WEIGHT*TAUXF
      ENDIF
*----
*  ADD CONTRIBUTIONS FROM A SINGLE ELEMENTARY CALCULATION
*----
      DO J=1,MAXH
         IF((HMAKE(J).NE.' ').AND.(HMAKE(J)(:4).NE.'SIGS')) THEN
            DO IG=1,NGRP
               GAR2(IG,J)=GAR2(IG,J)+WEIGHT*GAR1(IG,J)
            ENDDO
         ENDIF
      ENDDO
      DO IL=1,NL
         ITYPR(IL)=0
         IF(HMAKE(MAXH+IL).NE.' ') ITYPR(IL)=1
         DO IG2=1,NGRP
            GAR2(IG2,2+IL)=GAR2(IG2,2+IL)+WEIGHT*GAR1(IG2,2+IL)
            DO IG1=1,NGRP
               WSCA2(IG1,IG2,IL)=WSCA2(IG1,IG2,IL)+WEIGHT*
     1         WSCA1(IG1,IG2,IL)
            ENDDO
         ENDDO
      ENDDO
  120 CONTINUE
*----
*  SAVE ISOTOPIC DATA IN THE MICROLIB
*----
      CALL LCMPTC(IPLIB,'ALIAS',12,1,HNAME)
      IF(LWD) CALL LCMPUT(IPLIB,'LAMBDA-D',NDEL,2,WDLA)
      DO J=1,MAXH
         IF((HMAKE(J).NE.' ').AND.(HMAKE(J)(:4).NE.'SIGS')) THEN
            CALL LCMPUT(IPLIB,HMAKE(J),NGRP,2,GAR2(1,J))
         ENDIF
      ENDDO
      CALL XDRLGS(IPLIB,1,IMPX,0,NL-1,1,NGRP,GAR2(1,3),WSCA2,ITYPR)
      IF(IMPX.GT.50) CALL LCMLIB(IPLIB)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(WDLA,WSCA2,GAR2,WSCA1,GAR1)
      DEALLOCATE(ITYPR)
      RETURN
      END