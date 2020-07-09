*DECK LIBWRP
      SUBROUTINE LIBWRP(IPRINT,NTYP,NGR,NRTOT,MAXTEM,MAXDIL,IGR,IRES,
     >                  ITYP,DSIGPL,NTM,NDI,RTMP,RDIL,RESI,NTMPR,NDILR,
     >                  TMPT,DILT,REST)
C
C------------------------------  LIBWRP  ------------------------------
C
C  PROGRAMME STATISTICS:
C     NAME     : LIBWRP
C     ENTRY    : LIBWRP
C     USE      : PREPARE WIMS-D4 RESONANCE DATA
C     MODIFIED : 97-08-17
C     AUTHOR   : G. MARLEAU
C
C  ROUTINE PARAMETERS:
C   INPUT
C     IPRINT : PRINT FLAG                           I
C     NTYP   : NUMBER OF RESONANCE TABLES PER ISOTOPES (2 or 3) I
C     NGR    : NUMBER OF RESONANCE GROUPS           I
C     NRTOT  : MAMINUM NUMBER OF RESONANT ISOTOPES  I
C     MAXTEM : MAMINUM NUMBER OF TEMPERATURE        I
C     MAXDIL : MAMINUM NUMBER OF DILUTIONS          I
C     IGR    : RESONANCE GROUP NUMBER               I
C     IRES   : RESONANCE ISOTOPE SET                I
C     ITYP   : XS TYPE                              I
C     DSIGPL : BACKGROUND XS                        R
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
C   OUTPUT
C     NTMPR  : NUMBER OF LOCAL TEMPERATURES         I
C     NDILR  : NUMBER OF LOCAL DILUTIONS            I
C     TMPT   : WORK TEMPERATURE                     R(MAXTEM)
C     DILT   : WORK DILUTION                        R(MAXDIL)
C     REST   : WORK RESONANCE INTEGRALS             R(MAXDIL*
C                                                     MAXTEM)
C
C------------------------------  LIBWRP  ------------------------------
C
      IMPLICIT NONE
C----
C PARAMETERS
C----
      INTEGER   IOUT
      CHARACTER NAMSBR*6
      PARAMETER (IOUT=6,NAMSBR='LIBWRP')
C----
C INTERFACE VARIABLES
C----
      INTEGER   IPRINT,NTYP,NGR,NRTOT,MAXTEM,MAXDIL,IGR,IRES,ITYP,
     1          NTMPR,NDILR
      INTEGER   NTM(NTYP,NRTOT,NGR),NDI(NTYP,NRTOT,NGR)
      REAL      DSIGPL
      REAL      RTMP(MAXTEM,NTYP,NRTOT,NGR),RDIL(MAXDIL,NTYP,NRTOT,NGR),
     1          RESI(MAXDIL,MAXTEM,NTYP,NRTOT,NGR),TMPT(MAXTEM),
     2          DILT(MAXDIL),REST(MAXDIL*MAXTEM)
C----
C LOCAL VARIABLES
C----
      INTEGER   ITT,IT,IPOS,ID
      REAL      XDIL
C
C----
C
      NTMPR=NTM(ITYP,IRES,IGR)
      NDILR=NDI(ITYP,IRES,IGR)
      IF(ABS(IPRINT) .GE. 100) THEN
        WRITE(IOUT,6010) NAMSBR
        WRITE(IOUT,6000)
        WRITE(IOUT,6002) (RTMP(ITT,ITYP,IRES,IGR),ITT=1,NTMPR)
        WRITE(IOUT,6001)
        WRITE(IOUT,6002) (RDIL(ITT,ITYP,IRES,IGR),ITT=1,NDILR)
      ENDIF
      DO 100 IT=1,NTMPR
        TMPT(IT)=SQRT(RTMP(IT,ITYP,IRES,IGR))
 100  CONTINUE
      DO 110 ID=1,NDILR
        XDIL=RDIL(ID,ITYP,IRES,IGR)-DSIGPL
        IF(XDIL.GT.0.0) THEN
          DILT(ID)=SQRT(XDIL)
        ELSE
          DILT(ID)=0.0
        ENDIF
 110  CONTINUE
      IPOS=0
      DO 120 IT=1,NTMPR
        DO 121 ID=1,NDILR
          IPOS=IPOS+1
          REST(IPOS)=RESI(ID,IT,ITYP,IRES,IGR)
 121    CONTINUE
 120  CONTINUE
      IF(ABS(IPRINT) .GE. 100) THEN
        WRITE(IOUT,6011) NAMSBR
      ENDIF
      RETURN
C----
C  FORMAT
C----
 6000 FORMAT('   RESONANCE TEMPERATURE TABULATION = ')
 6001 FORMAT('   RESONANCE DILUTIONS TABULATION   = ')
 6002 FORMAT(1P,5E15.7)
 6010 FORMAT('(* Output from --',A6,'-- follows ')
 6011 FORMAT('   Output from --',A6,'-- completed *)')
      END
