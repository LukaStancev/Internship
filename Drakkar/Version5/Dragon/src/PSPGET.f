*DECK PSPGET
      SUBROUTINE PSPGET(IPRINT,ITYPE,ICOLR,NGROUP,NGT,ICOND)
C
C---------------------------  PSPGET  ---------------------------------
C
C  1- PROGRAMME STATISTICS:
C      NAME     : PSPGET
C      USE      : READ PSP OPTIONS PARAMETERS
C      AUTHOR   : G.MARLEAU
C      CREATED  : 99-01-21
C
C      MODIFICATION LOG
C      --------------------------------------------------------------
C      | DATE AND INITIALS  | MOTIVATIONS
C      --------------------------------------------------------------
C      | 99-01-21 G.M.      | READ EDIT LEVEL AND PROCESSING OPTIONS
C      ______________________________________________________________
C
C  2- ROUTINE PARAMETERS:
C    OUTPUT
C     IPRINT : PRINT LEVEL                            I
C     ITYPE  : TYPE OF GRAPHIC                        I
C              =  0 COLOR PER REGION NUMBER
C                 1 COLOR PER MATERIAL
C                 2 COLOR FOR FLUX (ONE GROUP)
C                 3 COLOR FOR FLUX (MULTIGROUP)
C                 4 COLOR PER MATERIAL FOR HOMOGENIZATION (HMIX)
C                 5 COLOR FOR MODE (ONE GROUP)
C                 6 COLOR FOR MODE (MULTIGROUP)
C     ICOLR  : COLOR SET USED                         I
C              = -4 FILL HSB WITH NO-CONTOUR
C              = -3 FILL CMYK WITH NO-CONTOUR
C              = -2 FILL RGB WITH NO-CONTOUR
C              = -1 FILL BW WITH NO-CONTOUR
C              =  0 NO FILL CONTOUR ONLY
C              =  1 FILL BW AND CONTOUR
C              =  2 FILL RGB AND CONTOUR
C              =  3 FILL CMYK AND CONTOUR
C              =  4 FILL HSB AND CONTOUR
C     NGROUP : number of groups for flux
C     NGT    : number of condensed groups for flux
C     ICOND  : upper group condensation limit
C  3- READ FORMAT
C     [ EDIT iprint ]
C     [ FILL  { NONE | GRAY | RGB | CMYK | HSB }  [ NOCONTOUR ]  ]
C     [ TYPE  { REGION | MIXTURE | FLUX | HMIX | 
C               MGFLUX (icond(i),i=1,ngt) }  ]
C       ;
C     DEFAULT:
C       IPRINT = 1 -> EDIT  1
C       ITYPE  = 0 -> PER REGION NUMBER
C       ICOLR  = 4 -> FILL HSB WITH CONTOUR
C
C---------------------------   PSPGET  --------------------------------
C
      IMPLICIT         NONE
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='PSPGET')
C----
C  ROUTINE PARAMETERS
C----
      INTEGER          IPRINT,ITYPE,ICOLR,NGROUP,NGT
      INTEGER          ICOND(NGROUP)
C----
C  REDGET INPUT VARIABLES
C----
      INTEGER          ITYPLU,INTLIR
      CHARACTER        CARLIR*12
      REAL             REALIR
      DOUBLE PRECISION DBLLIR
C----
C  LOCAL PARAMETERS
C----
      INTEGER          ICOL,ITY,ICONT,IGT
C----
C  READ OPTIONS
C----
      IPRINT=1
      ICOL=4
      ITY=0
      ICONT=1
 100  CONTINUE
      CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
 101  CONTINUE
      IF(ITYPLU .EQ. 10) THEN
        GO TO 105
      ELSE IF(ITYPLU .NE. 3) THEN
        CALL XABORT(NAMSBR//': ERROR -> CHARACTER VARIABLE EXPECTED')
      ENDIF
      IF(CARLIR(1:1) .EQ. ';' ) THEN
        GO TO 105
      ELSE IF(CARLIR .EQ. 'EDIT' ) THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 1 ) GO TO 101
        IPRINT=INTLIR
      ELSE IF(CARLIR(1:4) .EQ. 'FILL' ) THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 3 ) GO TO 101
        IF(CARLIR .EQ. 'NONE') THEN
          ICOL=0
        ELSE IF(CARLIR .EQ. 'GRAY') THEN
          ICOL=1
        ELSE IF(CARLIR .EQ. 'RGB') THEN
          ICOL=2
        ELSE IF(CARLIR .EQ. 'CMYK') THEN
          ICOL=3
        ELSE IF(CARLIR .EQ. 'HSB') THEN
          ICOL=4
        ELSE
          CALL XABORT(NAMSBR//': ILEGAL FILL KEYWORD '//CARLIR//
     >    'KEYWORD EXPECTED: NONE, GRAY, RGB, CMYK, HSB')
        ENDIF
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 3 ) GO TO 101
        IF(CARLIR(1:4) .EQ. 'NOCO') THEN
          ICONT=0
        ELSE
          GO TO 101
        ENDIF
      ELSE IF(CARLIR(1:4) .EQ. 'TYPE' ) THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(CARLIR(1:4) .EQ. 'REGI') THEN
          ITY=0
        ELSE IF(CARLIR(1:4) .EQ. 'MIXT') THEN
          ITY=1
        ELSE IF(CARLIR(1:4) .EQ. 'FLUX') THEN
          ITY=2
          NGT=1
          ICOND(NGT)=NGROUP
        ELSE IF(CARLIR(1:4) .EQ. 'MODE') THEN
          ITY=5
          NGT=1
          ICOND(NGT)=NGROUP
        ELSE IF(CARLIR(1:4) .EQ. 'MGFL') THEN
          ITY=3
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU .NE. 1) THEN
            NGT=NGROUP
            DO IGT=1,NGT
              ICOND(IGT)=IGT
            ENDDO
            GO TO 101
          ENDIF
          NGT=0
          DO IGT=1,NGROUP
            NGT=NGT+1
            IF(INTLIR .LT. 1 .OR. INTLIR .GT. NGROUP)
     >CALL XABORT(NAMSBR//': illegal group condensation number')
            IF(IGT .GT. 1) THEN
              IF(INTLIR .LE. ICOND(IGT-1))
     >CALL XABORT(NAMSBR//': group numbers must be increasing')
            ENDIF
            ICOND(IGT)=INTLIR
            CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
            IF(ITYPLU .NE. 1 ) THEN
              IF(ICOND(IGT) .NE. NGROUP) THEN
                NGT=NGT+1
                ICOND(NGT)=NGROUP
               ENDIF
              GO TO 101
            ENDIF
          ENDDO
        ELSE IF(CARLIR(1:4) .EQ. 'MGMD') THEN
          ITY=6
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU .NE. 1) THEN
            NGT=NGROUP
            DO IGT=1,NGT
              ICOND(IGT)=IGT
            ENDDO
            GO TO 101
          ENDIF
          NGT=0
          DO IGT=1,NGROUP
            NGT=NGT+1
            IF(INTLIR .LT. 1 .OR. INTLIR .GT. NGROUP)
     >CALL XABORT(NAMSBR//': illegal group condensation number')
            IF(IGT .GT. 1) THEN
              IF(INTLIR .LE. ICOND(IGT-1))
     >CALL XABORT(NAMSBR//': group numbers must be increasing')
            ENDIF
            ICOND(IGT)=INTLIR
            CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
            IF(ITYPLU .NE. 1 ) THEN
              IF(ICOND(IGT) .NE. NGROUP) THEN
                NGT=NGT+1
                ICOND(NGT)=NGROUP
               ENDIF
              GO TO 101
            ENDIF
          ENDDO
        ELSE IF(CARLIR(1:4) .EQ. 'HMIX') THEN
          ITY=4
        ELSE
          CALL XABORT(NAMSBR//': ILEGAL TYPE KEYWORD '//CARLIR//
     >    'KEYWORD EXPECTED: REGION, MIXTURE, FLUX, MGFLUX')
        ENDIF
      ELSE
C----
C  INVALID OPTION
C----
        CALL XABORT(NAMSBR//': ILEGAL MAIN KEYWORD '//CARLIR//
     >  'KEYWORD EXPECTED: FILL, TYPE, EDIT OR ; ')
      ENDIF
      GO TO 100
 105  CONTINUE
C----
C  TEST READ OPTIONS
C  IF FILL = NONE (ICOL = 0) IMPOSE CONTOUR
C----
      IF(ICONT .EQ. 0) THEN
        ICOLR=-ICOL
      ELSE
        ICOLR=ICOL
      ENDIF
      ITYPE=ITY
C----
C  PRINT ECHO OF PSP OPTIONS THAT WILL BE USED
C----
      IF(IPRINT .GE. 1 ) THEN
        WRITE(IOUT,6000) IPRINT,ICOL,ITY,ICONT
      ENDIF
C----
C  RETURN
C----
      RETURN
C----
C  FORMAT
C----
 6000 FORMAT(' ------  PSP EXECUTION OPTIONS --------'/
     >       ' PRINT LEVEL = ',I8                    /
     >       ' COLOR       = ',I8                    /
     >       ' TYPE        = ',I8                    /
     >       ' CONTOUR     = ',I8                    /
     >       ' --------------------------------------')
      END
