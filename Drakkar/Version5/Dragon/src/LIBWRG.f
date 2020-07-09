*DECK LIBWRG
      SUBROUTINE LIBWRG(IUNIT,NTYP,NGR,NRTOT,MAXTEM,MAXDIL,NSRES,RID,
     >                  NTM,NDI,RTMP,RDIL,RESI)
C
C------------------------------  LIBWRG  ------------------------------
C
C  PROGRAMME STATISTICS:
C     NAME     : LIBWRG
C     ENTRY    : LIBWRG
C     USE      : READ RESONANCE INFORMATION FROM WIMS-D4 LIBRARY
C     MODIFIED : 97-01-30
C     AUTHOR   : G. MARLEAU
C
C  ROUTINE PARAMETERS:
C   INPUT
C     IUNIT  : WIMS-D4 READ UNIT                    I
C     NTYP   : NUMBER OF RESONANCE TABLES PER ISOTOPES (2 or 3) I
C     NGR    : NUMBER OF RESONANCE GROUPS           I
C     NRTOT  : NUMBER OF RESONANCE SETS             I
C     MAXTEM : MAX NB TEMPERATURE                   I
C     MAXDIL : MAX NB DILUTIONS                     I
C     NSRES  : NB OF RESONANCE SET                  I
C     RID    : RESONANCE ID                         I(NRTOT)
C     NTM    : NUMBER OF TEMPERATURES               I(NTYP,NRTOT,
C                                                     NGR)
C     NDI    : NUMBER OF DILUTIONS                  I(NTYP,NRTOT,
C                                                     NGR)
C     RTMP   : RESONANCE TEMPERATURE                R(MAXTEM,
C                                                NTYP,NRTOT,NGR)
C     RDIL   : RESONANCE DILUTION                   R(MAXDIL,
C                                             NTYP,NRTOT,NGR)
C     RESI   : RESONANCE INTEGRALS                  R(MAXDIL,
C                                      MAXTEM,NTYP,NRTOT,NGR)
C
C------------------------------  LIBWRG  ------------------------------
C
      IMPLICIT NONE
C----
C PARAMETERS
C----
      INTEGER    IOUT
      PARAMETER (IOUT=6)
C----
C INTERFACE PARAMETERS
C----
      INTEGER    IUNIT,NTYP,NGR,NRTOT,MAXTEM,MAXDIL
      INTEGER    NTM(NTYP,NRTOT,NGR),NDI(NTYP,NRTOT,NGR)
C
      REAL       RID(NRTOT),RTMP(MAXTEM,NTYP,NRTOT,NGR),
     1           RDIL(MAXDIL,NTYP,NRTOT,NGR),
     2           RESI(MAXDIL,MAXTEM,NTYP,NRTOT,NGR)
C----
C LOCAL VARIABLES
C----
      INTEGER    IGR,NSRES,ISRES,IPREV,IRS,M1,M2,IT,ID,ISR,ITYP,
     1           NTIS
      REAL       XIDR,ENDR
C----
C ALLOCATABLE ARRAYS
C----
      REAL, ALLOCATABLE, DIMENSION(:) :: TMPT,DILT
      REAL, ALLOCATABLE, DIMENSION(:,:) :: REST
C----
C  SCRATCH STORAGE ALLOCATION
C     TMPT   : TEMPERATURE
C     DILT   : DILUTION
C     REST   : RESONANCE INTEGRALS
C----
      ALLOCATE(TMPT(MAXTEM),DILT(MAXDIL),REST(MAXDIL,MAXTEM))
C----
C  SCAN OVER RESONANCE GROUPS
C----
      NSRES=0
      ISRES=0
      DO 100 IGR=1,NGR
        IPREV=0
C----
C  SCAN OVER RESONANCE SETS + 1
C  AND READ RESONANCE INFO
C----
        DO 110 IRS=1,NTYP*NRTOT+1
          READ(IUNIT) XIDR,M1,M2,
     >     (TMPT(IT),IT=1,M1),(DILT(ID),ID=1,M2),
     >    ((REST(ID,IT),ID=1,M2),IT=1,M1)
          IF(XIDR.EQ.0.0) GO TO 115
          IF((M1.EQ.0).AND.(M2.EQ.0)) GO TO 110
          DO 120 ISR=1,NSRES
            IF(XIDR.EQ.RID(ISR)) THEN
              ISRES=ISR
              GO TO 125
            ENDIF
 120      CONTINUE
          NSRES=NSRES+1
          IF(NSRES.GT.NRTOT) THEN
            CALL XABORT('LIBWRG: TO MANY RESONANCE SET')
          ENDIF
          ISRES=NSRES
          IPREV=0
          RID(ISRES)=XIDR
 125      CONTINUE
          IF(ISRES.NE.IPREV) THEN
            ITYP=1
            IPREV=ISRES
          ELSE IF((ISRES.EQ.IPREV).AND.(ITYP.EQ.1)) THEN
            ITYP=2
          ELSE IF((ISRES.EQ.IPREV).AND.(ITYP.EQ.2)) THEN
            ITYP=3
            IPREV=0
          ENDIF
          NTIS=NTM(ITYP,ISRES,IGR)
          IF(NTIS.GT.0) THEN
            WRITE(IOUT,9000) IGR,ISRES,ITYP,XIDR
            CALL XABORT('LIBWRG: DUPLICATE RESONANCE SET')
          ENDIF
C----
C  SAVE RESONANCE INFORMATION FOR THIS SET
C----
          NTM(ITYP,ISRES,IGR)=M1
          NDI(ITYP,ISRES,IGR)=M2
          DO 130 IT=1,M1
            RTMP(IT,ITYP,ISRES,IGR)=TMPT(IT)
 130      CONTINUE
          DO 131 ID=1,M2
            RDIL(ID,ITYP,ISRES,IGR)=DILT(ID)
 131      CONTINUE
          DO 140 IT=1,M1
            DO 141 ID=1,M2
              RESI(ID,IT,ITYP,ISRES,IGR)=REST(ID,IT)
 141        CONTINUE
 140      CONTINUE
 110    CONTINUE
 115    CONTINUE
        IF(NTYP.EQ.2) READ(IUNIT) ENDR
 100  CONTINUE
C----
C  SCRATCH STORAGE DEALLOCATION
C----
      DEALLOCATE(REST,DILT,TMPT)
      RETURN
C----
C  FORMAT
C----
 9000 FORMAT(' LIBWRG ERROR - WIMS-D4 DUPLICATE RESONANCE SET'/
     >       ' RESONANCE GROUP = ',I10/
     >       '   RESONANCE SET = ',I10/
     >       '   INTEGRAL TYPE = ',I10/
     >       '    RESONANCE ID = ',F20.5)
      END
