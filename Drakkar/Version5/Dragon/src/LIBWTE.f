*DECK LIBWTE
      SUBROUTINE LIBWTE(IACT,ITXS,NGROUP,NGTHER,NTMP,NF,TERP,SCAT,
     >                  SIGS,XSNG,SIGF,XSFI,TRAN,TMPXS,TMPSC)
C
C------------------------------  LIBWTE  ------------------------------
C
C  PROGRAMME STATISTICS:
C     NAME     : LIBWTE
C     ENTRY    : LIBWTE
C     USE      : PERFORM TEMPERATURE INTERPOLATION FOR
C                WIMS-AECL OR WIMS-D4 XS
C     MODIFIED : 97-08-17
C     AUTHOR   : G. MARLEAU
C
C  ROUTINE PARAMETERS:
C   INPUT
C     IACT   : ACTION                               I
C              = 1 INITIALIZE BEFORE ADDING
C              = 2 ONLY ADD
C     ITXS   : TYPE                                 I
C              = 1 ALL CROSS SECTIONS
C              = 2 ONLY SCATTERING
C     NGROUP : NUMBER OF GROUPS                     I
C     NGTHER : NUMBER OF THERMAL GROUPS             I
C     NTMP   : NUMBER OF TEMPERATURES               I
C     NF     : FISSILE FLAG                         I
C     TERP   : TEMPERATURE COEFFICIENTS             D(NTMP)
C   INPUT/OUTPUT
C     SCAT   : COMPLETE SCATTERING MATRIX           R(NGROUP,
C              SCAT(JG,IG) IS FROM IG TO JG           NGROUP)
C     SIGS   : TOTAL SCATTERING OUT OF GROUP        R(NGROUP)
C     XSNG   : NG XS                                R(NGROUP)
C     SIGF   : NU*FISSION XS                        R(NGROUP)
C     XSFI   : FISSION XS                           R(NGROUP)
C     TRAN   : TRANSPORT XS                         R(NGROUP)
C   WORK
C     TMPXS  : TEMPERATURE DEPENDENT VECT XS        R(NGTHER,5,NTMP)
C     TMPSC  : TEMPERATURE DEPENDENT SCAT XS        R(NGROUP,
C                                               NGROUP,NTMP)
C
C------------------------------  LIBWTE  ------------------------------
C
      IMPLICIT NONE
C----
C INTERFACE VARIABLES
C----
      INTEGER          IACT,ITXS,NGROUP,NGTHER,NTMP,NF
      DOUBLE PRECISION TERP(NTMP)
      REAL             SCAT(NGROUP,NGROUP),SIGS(NGROUP),
     1                 XSNG(NGROUP),SIGF(NGROUP),XSFI(NGROUP),
     2                 TRAN(NGROUP),TMPXS(NGROUP,5,NTMP),
     3                 TMPSC(NGROUP,NGROUP,NTMP)
C----
C LOCAL VARIABLES
C----
      INTEGER          IGF,ITM,IGD,NGD
      REAL             RTERP
C----
C  INITIALIZED IF REQUIRED
C----
      NGD=NGROUP-NGTHER+1
      IF(IACT.EQ.1) THEN
        IF(ITXS.EQ.1) THEN
          CALL XDRSET(XSNG(NGD),NGTHER,0.0)
          CALL XDRSET(TRAN(NGD),NGTHER,0.0)
          IF(NF.GT.1) THEN
            CALL XDRSET(SIGF(NGD),NGTHER,0.0)
            CALL XDRSET(XSFI(NGD),NGTHER,0.0)
          ENDIF
        ENDIF
        IF(ITXS.GE.1) THEN
          CALL XDRSET(SIGS(NGD),NGTHER,0.0)
          DO 110 IGD=NGD,NGROUP
            CALL XDRSET(SCAT(1,IGD),NGROUP,0.0)
 110      CONTINUE
        ENDIF
      ENDIF
C----
C  INTERPOLATE STANDARD CROSS SECTIONS IN TEMPERATURE
C----
      IF(ITXS.EQ.1) THEN
        DO 120 ITM=1,NTMP
          RTERP=REAL(TERP(ITM))
          IF(RTERP.NE.0.0) THEN
            DO 121 IGD=NGD,NGROUP
              TRAN(IGD)=TRAN(IGD)+RTERP*TMPXS(IGD,1,ITM)
              XSNG(IGD)=XSNG(IGD)+RTERP*TMPXS(IGD,2,ITM)
              IF(NF.GT.1) THEN
                SIGF(IGD)=SIGF(IGD)+RTERP*TMPXS(IGD,3,ITM)
                XSFI(IGD)=XSFI(IGD)+RTERP*TMPXS(IGD,4,ITM)
              ENDIF
 121        CONTINUE
          ENDIF
 120    CONTINUE
      ENDIF
C----
C  INTERPOLATE SCATTERING CROSS SECTIONS IN TEMPERATURE
C----
      IF(ITXS.GE.1) THEN
        DO 130 ITM=1,NTMP
          RTERP=REAL(TERP(ITM))
          IF(RTERP.NE.0.0D0) THEN
            DO 131 IGD=NGD,NGROUP
              SIGS(IGD)=SIGS(IGD)+RTERP*TMPXS(IGD,5,ITM)
              DO 132 IGF=1,NGROUP
                SCAT(IGF,IGD)=SCAT(IGF,IGD)+RTERP*TMPSC(IGF,IGD,ITM)
 132          CONTINUE
 131        CONTINUE
          ENDIF
 130    CONTINUE
      ENDIF
      RETURN
      END
