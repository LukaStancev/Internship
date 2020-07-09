*DECK LIBWRE
      SUBROUTINE LIBWRE(NTYP,IPRINT,ITLIB ,NGROUP,NL,IGRF,IGRL,NGR,
     >                  SCAT,SIGS,TOTAL,XSNG,SIGF,XSFI,XNU,DELTA,
     >                  DIL,DLJ,XSOUT,XSCOR,DSIGPL)
C
C-------------------------- LIBWRE   ----------------------------------
C
C  PROGRAMME STATISTICS:
C     NAME     : LIBWRE
C     ENTRY    : LIBWRE
C     USE      : RESONANCE INTEGRAL TEMP AND DIL INTERPOLATION
C     MODIFIED : 93-08-17
C     AUTHOR   : G. MARLEAU
C
C  ROUTINE PARAMETERS:
C   INPUT
C     NTYP   : TYPES OF SELF SHIELDING RATES        I
C              = 1 ONLY ABSORPTION
C              = 2 ABSORPTION+FISSION
C              = 3 ABSORPTION+FISSION+SCATT
C     IPRINT : PRINT FLAG                           I
C     ITLIB  : TYPES OF LIBRARY                     I
C              = 1 WIMS-AECL
C              = 2 WIMS-D4
C     NGROUP : NUMBER OF GROUPS                     I
C     NL     : NUMBER OF LEGENDRE SCATTERING ORDER  I
C     IGRF   : FIRST RESONANCE GROUP TO TREAT       I
C     IGRL   : LAST RESONANCE GROUP TO TREAT        I
C     NGR    : NUMBER OF RESONANCE GROUPS           I
C     SCAT   : COMPLETE SCATTERING MATRIX           R(NGROUP,
C              SCAT(JG,IG) IS FROM IG TO JG          NGROUP,NL)
C     SIGS   : TOTAL SCATTERING OUT OF GROUP        R(NGROUP)
C     TOTAL  : TOTAL XS                             R(NGROUP)
C     XSNG   : NG XS                                R(NGROUP)
C     SIGF   : NU*FISSION XS                        R(NGROUP)
C     XSFI   : FISSION XS                           R(NGROUP)
C     XNU    : 1/NU                                 R(NGROUP)
C     DELTA  : LETHARGY                             R(NGROUP)
C     DIL    : STANDARD DILUTION                    R(NGROUP)
C     DLJ    : LIVOLANT-JEANPIERRE DILUTION         R(NGROUP)
C     XSOUT  : RESONNANCES INTEGRALS                R(NGROUP,7)
C     XSCOR  : TOTAL CORRECTION                     R(4)
C     DSIGPL : POTENTIAL XS TIMES G-C PARAMETERS    R(NGR)
C
C-------------------------- LIBWRE   ----------------------------------
C
      IMPLICIT NONE
      INTEGER     IOUT
      CHARACTER   NAMSBR*6
      PARAMETER  (IOUT=6,NAMSBR='LIBWRE')
C----
C INTEFACE VRAIABLES
C----
      INTEGER    NTYP,IPRINT,ITLIB,NGROUP,NL,IGRF,IGRL,NGR
      REAL       SCAT(NGROUP,NGROUP,NL),SIGS(NGROUP,NL),TOTAL(NGROUP),
     1           XSNG(NGROUP),SIGF(NGROUP),XSFI(NGROUP),
     2           XNU(NGROUP),DELTA(NGROUP),DIL(NGROUP),
     3           DLJ(NGROUP),XSOUT(NGROUP,7),XSCOR(4),DSIGPL(NGR)
C----
C LOCAL VARIABLES
C----
      INTEGER          IGRR,JG,IG1,IG2,IL
      REAL             XSF,DDIL,DDLJ
      DOUBLE PRECISION XNUMER,XDENOM
C
C-----
C
      IF(ABS(IPRINT) .GE. 100) THEN
        WRITE(IOUT,6000) NAMSBR 
      ENDIF
      DO 100 IGRR=IGRF,IGRL
        DDLJ=DLJ(IGRR)
        DDIL=DIL(IGRR)
        IF(ABS(IPRINT) .GE. 1000) THEN
          WRITE(IOUT,*) 'Potential XS*GC parameter ',DSIGPL(IGRR-IGRF+1)
        ENDIF
        IF(NTYP.EQ.3.AND.XSCOR(3).GT.0.0) THEN
C----
C  COMPUTE FLUX
C  SCATTERING IS SELF-SHIELDED
C----
          IF(IGRR .EQ. IGRF) THEN
            IF(ABS(IPRINT) .GE. 100) THEN
              WRITE(IOUT,6020)  
            ENDIF
          ENDIF
          IF(ITLIB.EQ.1) THEN
            XNUMER=DBLE(DDLJ-XSOUT(IGRR,1))
            DO 110 JG=1,IGRR-1
              XNUMER=XNUMER+DBLE(
     >        XSOUT(JG,4)*DELTA(JG)*SCAT(IGRR,JG,1)/DELTA(IGRR))
 110        CONTINUE
            XDENOM=DBLE(DDLJ+SIGS(IGRR,1)-SCAT(IGRR,IGRR,1))
            XSOUT(IGRR,4)=REAL(XNUMER/XDENOM)
            XNUMER=DBLE(DDIL-XSOUT(IGRR,1))
            DO 115 JG=1,IGRR-1
             XNUMER=XNUMER+DBLE(
     >        XSOUT(JG,5)*DELTA(JG)*SCAT(IGRR,JG,1)/DELTA(IGRR))
 115        CONTINUE
            XDENOM=DBLE(DDIL+SIGS(IGRR,1)-SCAT(IGRR,IGRR,1))
            XSOUT(IGRR,5)=REAL(XNUMER/XDENOM)
          ELSE
            XSOUT(IGRR,4)=(DDLJ-XSOUT(IGRR,1))/DDLJ
            XSOUT(IGRR,5)=(DDIL-XSOUT(IGRR,1))/DDIL
          ENDIF
          IF(ABS(IPRINT) .GE. 1000) THEN
            WRITE(IOUT,*) 'Flux 2 ',XSOUT(IGRR,4)
          ENDIF
          IF(XNU(IGRR).NE.0.0) THEN
             XSOUT(IGRR,6)=XSOUT(IGRR,2)/XNU(IGRR)
          ELSE
             XSOUT(IGRR,6)=0.0
          ENDIF
          XSOUT(IGRR,7)=XSOUT(IGRR,1)-XSOUT(IGRR,6)
          IF(ABS(IPRINT) .GE. 1000) THEN
            WRITE(IOUT,*) 'Total ',XSOUT(IGRR,1),XSOUT(IGRR,3),
     >                 XSOUT(IGRR,1)+XSOUT(IGRR,3)
          ENDIF
          XSOUT(IGRR,1)=XSOUT(IGRR,1)+XSOUT(IGRR,3)
        ELSE IF(XSCOR(1).GT.0.0) THEN
C----
C  COMPUTE FLUX AND DRAGLIB FLUX
C  SCATTERING IS NOT SELF-SHIELDED
C----
          IF(IGRR .EQ. IGRF) THEN
            IF(ABS(IPRINT) .GE. 100) THEN
              WRITE(IOUT,6021)  
            ENDIF
          ENDIF
          IF(ITLIB.EQ.1) THEN
            XNUMER=DBLE(DDLJ-XSOUT(IGRR,1))
            DO 120 JG=1,IGRR-1
              XNUMER=XNUMER+DBLE(
     >          XSOUT(JG,4)*DELTA(JG)*SCAT(IGRR,JG,1)/DELTA(IGRR))
 120        CONTINUE
            XDENOM=DBLE(DDLJ+SIGS(IGRR,1)-SCAT(IGRR,IGRR,1))
            XSOUT(IGRR,4)=REAL(XNUMER/XDENOM)
            IF(ABS(IPRINT) .GE. 1000) THEN
              WRITE(IOUT,*) 'Flux 1 ',XSOUT(IGRR,4),XNUMER,XDENOM
            ENDIF
            XNUMER=DBLE(DDIL-XSOUT(IGRR,1))
            DO 130 JG=1,IGRR-1
              XNUMER=XNUMER+DBLE(
     >          XSOUT(JG,5)*DELTA(JG)*SCAT(IGRR,JG,1)/DELTA(IGRR))
 130        CONTINUE
            XDENOM=DBLE(DDIL+SIGS(IGRR,1)-SCAT(IGRR,IGRR,1))
            XSOUT(IGRR,5)=REAL(XNUMER/XDENOM)
          ELSE
            XSOUT(IGRR,4)=(DDLJ-XSOUT(IGRR,1))/DDLJ
            XSOUT(IGRR,5)=(DDIL-XSOUT(IGRR,1))/DDIL
            IF(ABS(IPRINT) .GE. 1000) THEN
              WRITE(IOUT,*) 'Flux 2 ',XSOUT(IGRR,4)
            ENDIF
          ENDIF
          XSOUT(IGRR,3)=SIGS(IGRR,1)*XSOUT(IGRR,4)
          IF(NTYP.LT.2) THEN
            XSOUT(IGRR,2)=SIGF(IGRR)*XSOUT(IGRR,4)
          ENDIF
          IF(XNU(IGRR).NE.0.0) THEN
            XSOUT(IGRR,6)=XSOUT(IGRR,2)/XNU(IGRR)
          ELSE
            XSOUT(IGRR,6)=0.0
          ENDIF
          IF(ABS(IPRINT) .GE. 1000) THEN
            WRITE(IOUT,*) 'Total ',XSOUT(IGRR,1),
     >      SIGS(IGRR,1)*XSOUT(IGRR,5),
     >      XSOUT(IGRR,1)+SIGS(IGRR,1)*XSOUT(IGRR,5)
          ENDIF
          XSOUT(IGRR,1)=SIGS(IGRR,1)*XSOUT(IGRR,5)+XSOUT(IGRR,1)
          XSOUT(IGRR,7)=XSOUT(IGRR,1)-XSOUT(IGRR,3)-XSOUT(IGRR,6)
        ELSE IF(XSOUT(IGRR,4).NE.0.0) THEN
          IF(IGRR .EQ. IGRF) THEN
            IF(ABS(IPRINT) .GE. 100) THEN
              WRITE(IOUT,6022)  
            ENDIF
          ENDIF
          IF(XNU(IGRR).NE.0.0) THEN
            XSOUT(IGRR,6)=SIGF(IGRR)*XSOUT(IGRR,4)/XNU(IGRR)
          ELSE
            XSOUT(IGRR,6)=0.0
          ENDIF
          XSOUT(IGRR,1)=TOTAL(IGRR)*XSOUT(IGRR,4)
          XSOUT(IGRR,2)=SIGF(IGRR)*XSOUT(IGRR,4)
          XSOUT(IGRR,3)=SIGS(IGRR,1)*XSOUT(IGRR,4)
          XSOUT(IGRR,7)=XSOUT(IGRR,1)-XSOUT(IGRR,3)-XSOUT(IGRR,6)
          IF(ABS(IPRINT) .GE. 1000) THEN
            WRITE(IOUT,*) 'Total ',TOTAL(IGRR)*XSOUT(IGRR,4)
          ENDIF
        ELSE
          IF(IGRR .EQ. IGRF) THEN
            IF(ABS(IPRINT) .GE. 100) THEN
              WRITE(IOUT,6023)  
            ENDIF
          ENDIF
          XSOUT(IGRR,4)=1.0
          IF(XNU(IGRR).NE.0.0) THEN
            XSOUT(IGRR,6)=SIGF(IGRR)/XNU(IGRR)
          ELSE
            XSOUT(IGRR,6)=0.0
          ENDIF
          XSOUT(IGRR,1)=TOTAL(IGRR)
          XSOUT(IGRR,2)=SIGF(IGRR)
          XSOUT(IGRR,3)=SIGS(IGRR,1)
          XSOUT(IGRR,7)=XSOUT(IGRR,1)-XSOUT(IGRR,3)-XSOUT(IGRR,6)
          IF(ABS(IPRINT) .GE. 1000) THEN
            WRITE(IOUT,*) 'Total ',TOTAL(IGRR)
          ENDIF
        ENDIF
 100  CONTINUE
      DO 141 IG2=IGRF,IGRL
        XSF=XSOUT(IG2,3)/(XSOUT(IG2,4)*SIGS(IG2,1))
        DO 140 IL=1,NL
        DO 140 IG1=1,NGROUP
          SCAT(IG1,IG2,IL)=XSF*SCAT(IG1,IG2,IL)
 140    CONTINUE
 141  CONTINUE
      IF(ABS(IPRINT) .GE. 100) THEN
        WRITE(IOUT,6010) 'FLUX            '
        WRITE(IOUT,6011) (XSOUT(IG1,4),IG1=IGRF,IGRL)
        WRITE(IOUT,6010) 'TOTAL RATE      '
        WRITE(IOUT,6011) (XSOUT(IG1,1),IG1=IGRF,IGRL)
        WRITE(IOUT,6010) 'TOTAL XS        '
        WRITE(IOUT,6011) (XSOUT(IG1,1)/XSOUT(IG1,4),IG1=IGRF,IGRL)
        WRITE(IOUT,6010) 'FISSION RATE    '
        WRITE(IOUT,6011) (XSOUT(IG1,2),IG1=IGRF,IGRL)
        WRITE(IOUT,6010) 'FISSION XS      '
        WRITE(IOUT,6011) (XSOUT(IG1,2)/XSOUT(IG1,4),IG1=IGRF,IGRL)
        WRITE(IOUT,6010) 'SCATTERING RATE '
        WRITE(IOUT,6011) (XSOUT(IG1,3),IG1=IGRF,IGRL)
        WRITE(IOUT,6010) 'SCATTERING XS   '
        WRITE(IOUT,6011) (XSOUT(IG1,3)/XSOUT(IG1,4),IG1=IGRF,IGRL)
        WRITE(IOUT,6010) 'NG RATE         '
        WRITE(IOUT,6011) (XSOUT(IG1,7),IG1=IGRF,IGRL)
        WRITE(IOUT,6010) 'NG XS           '
        WRITE(IOUT,6011) (XSOUT(IG1,7)/XSOUT(IG1,4),IG1=IGRF,IGRL)
        WRITE(IOUT,6001) NAMSBR 
      ENDIF
      DO 150 IG1=IGRF,IGRL
        TOTAL(IG1)=XSOUT(IG1,1)/XSOUT(IG1,4)
        SIGF(IG1)=XSOUT(IG1,2)/XSOUT(IG1,4)
        XSF=XSOUT(IG1,3)/(XSOUT(IG1,4)*SIGS(IG1,1))
        SIGS(IG1,1)=XSOUT(IG1,3)/XSOUT(IG1,4)
        DO IL=2,NL
          SIGS(IG1,IL)=XSF*SIGS(IG1,IL)
        ENDDO
        XSFI(IG1)=XSOUT(IG1,6)/XSOUT(IG1,4)
        XSNG(IG1)=XSOUT(IG1,7)/XSOUT(IG1,4)
 150  CONTINUE
C----
C  RETURN LIBWRE
C----
      RETURN
C----
C  FORMAT
C----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6010 FORMAT('   Record  = ',A16)
 6011 FORMAT(1P,5E15.7)
 6020 FORMAT(' Flux computed with self shielded scattering XS.')
 6021 FORMAT(' Flux computed without self shielded scattering XS.')
 6022 FORMAT(' Flux tabulated.')
 6023 FORMAT(' Flux initialized to unity.')
      END
