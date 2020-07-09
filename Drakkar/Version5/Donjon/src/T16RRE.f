*DECK T16RRE
      SUBROUTINE T16RRE(IFT16 ,IPRINT,NGCCPO,NGMTR ,IFGMTR,NVXSR ,
     >                  NMXSR ,IMIREG,VELMTR,B2CRI ,BRNIRR,
     >                  FLXINT,OVERV,RECXSV,RECXSM,RECTMP,RECSCA)
C
C----
C  1- PROGRAMME STATISTICS:
C      NAME     : T16RRE
C      USE      : READ TAPE16 REGION CROSS SECTIONS
C                 AT A SPECIFIC BURNUP
C      AUTHOR   : G.MARLEAU
C      CREATED  : 1999/10/21
C      REF      : IGE-244 REV.1
C
C      MODIFICATION LOG
C      --------------------------------------------------------------
C      | DATE AND INITIALS  | MOTIVATIONS
C      --------------------------------------------------------------
C      | 1999/10/21 G.M.    | READ TAPE16 CROSS SECTION
C      |                    | ASSOCIATED WITH A REGION
C      |                    | AT A SPECIFIC BURNUP
C      --------------------------------------------------------------
C
C  2- ROUTINE PARAMETERS:
C    INPUT
C      IFT16  : TAPE16 FILE UNIT                         I
C      IPRINT : PRINT LEVEL                              I
C               =   0 NO PRINT
C               >=  1 PRINT PROCESSING OPTIONS READ
C      NGCCPO : NUMBER OF FINAL CONDENSED GROUPS         I
C      NGMTR  : NUMBER OF MAIN TRANSPORT GROUP           I
C      IFGMTR : CPO FEW GROUP IDENTIFIER                 I(NGCCPO)
C               WITH RESPECT TO MTR GROUPS
C      NVXSR  : NUMBER OF VECTORIAL XS                   I
C      NMXSR  : NUMBER OF MATRIX XS                      I
C      IMIREG : REGION TO CONSIDER FOR MIXTURE           I
C      VELMTR : AVERAGE VELOCITY IN MAIN GROUPS          R(NGMTR)
C      B2CRI  : CRITICAL BUCKLINGS                       R(3)
C    OUTPUT
C      BRNIRR : BURNUP IRRADIATION ENERGY                R(3)
C      FLXINT : VOLUME INTEGRATED FLUX                   R(NGCCPO)
C      OVERV  : 1/V CROSS SECTION                        R(NGCCPO)
C      RECXSV : VECTOR CROSS SECTIONS RECORDS            R(NGCCPO,
C                                                     NVXSR+NMXSR)
C      RECXSM : MATRIX CROSS SECTIONS RECORDS            R(NGCCPO,
C               FORMAT OF RECXSM IS                  NGCCPO,NMXSR)
C               RECXSM(IGTO,IGFROM,IL) REPRESENT
C               SCATTERING CROSS SECTION
C               FROM GROUP "IGFROM" TO GROUP "IGTO"
C               FOR ANISOTROPY LEVEL IL
C    WORK
C      RECTMP : VECTOR CROSS SECTIONS RECORDS            R(NGMTR,4)
C      RECSCA : SCATT CROSS SECTIONS RECORDS             R(NGMTR,
C                                                          NGMTR)
C
C  3- ROUTINES CALLED
C    SPECIFIC T16CPO ROUTINES
C      T16FND : FIND A TAPE16 RECORD
C               EQUIVALENT TO FIND FUNCTION
C               IN APPENDIX E OF EACL RC-1176
C    UTILITIES ROUTINES
C      XABORT : ABORT ROUTINE
C      XDRSET : VECTOR INITIALIZATION ROUTINE
C
C----
C
      IMPLICIT         NONE
      INTEGER          IFT16,IPRINT,NGCCPO,NGMTR,NVXSR,NMXSR,IMIREG
      INTEGER          IFGMTR(NGCCPO)
      REAL             VELMTR(NGMTR),B2CRI(3),BRNIRR(3),
     >                 FLXINT(NGCCPO),OVERV(NGCCPO),
     >                 RECXSV(NGCCPO,NVXSR+NMXSR),
     >                 RECXSM(NGCCPO,NGCCPO,NMXSR),
     >                 RECTMP(NGMTR,4),RECSCA(NGMTR,NGMTR)
C----
C  T16 PARAMETERS
C----
      INTEGER          MAXKEY
      PARAMETER       (MAXKEY=2)
      CHARACTER        TKEY1(MAXKEY)*10,TKEY2(MAXKEY)*10,
     >                 RKEY1*10,RKEY2*10
      INTEGER          NKEY,IOPT,NBE,NID,NJD
C----
C  LOCAL VARIABLES
C  WSMEV FACTOR TO TRANSFORM MEV IN JOULES (WS)
C----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      REAL             WSMEV
      PARAMETER       (IOUT=6,NAMSBR='T16RRE',WSMEV=1.602189E-13)
      INTEGER          IREG,IGR,IGC,IGD,IGF,JGR,JGC,JGD,JGF,
     >                 NREGON
      REAL             VOLUME,BRNTMP(3),RTIME
C----
C  INITIALIZE CROSS SECTION VECTORS
C----
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
      CALL XDRSET(RECXSV,NGCCPO*(NVXSR+NMXSR),0.0)
      CALL XDRSET(RECXSM,NGCCPO*NGCCPO*NMXSR,0.0)
C----
C  LOCATE NEXT REGION DIMENSIONS RECORD
C  AND READ NREGON
C----
      IOPT=0
      TKEY1(1)='REGION    '
      TKEY2(1)='DIMENSIONS'
      NKEY=1
      CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1 ,TKEY2 ,
     >            NBE   )
      IF( NBE .NE. 2 ) CALL XABORT(NAMSBR//
     >': CANNOT FIND '//TKEY1(1)//' '//TKEY2(1))
      READ(IFT16) RKEY1,RKEY2,NBE,NREGON
      TKEY1(2)='CELLAV    '
      TKEY2(2)='NGROUPS   '
      NKEY=2
      DO 100 IREG=1,NREGON
C----
C  REGIONAL FLUX
C----
        TKEY1(1)='REGION    '
        TKEY2(1)='FLUX      '
        CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1 ,TKEY2 ,
     >              NBE   )
        IF( NBE .NE. 3+NGMTR ) CALL XABORT(NAMSBR//
     >  ': CANNOT FIND '//TKEY1(1)//' '//TKEY2(1))
        IF(IMIREG .EQ. IREG) THEN
          READ(IFT16) RKEY1,RKEY2,NBE,NID,NJD,VOLUME,
     >               (RECTMP(IGR,1),IGR=1,NGMTR)
          IF(IPRINT .GE. 100) THEN
            WRITE(IOUT,6100) TKEY2(1)
            WRITE(IOUT,6110) (RECTMP(IGR,1),IGR=1,NGMTR)
          ENDIF
C----
C  TREAT ALL CONDENSED GROUPS
C----
          TKEY1(1)='REGION    '
          TKEY2(1)='SIGMAS    '
          IGF=0
          DO 110 IGC=1,NGCCPO
            IGD=IGF+1
            IGF=IFGMTR(IGC)
C----
C  FLUX AND 1/V CROSS SECTION CONDENSATION
C----
            DO 111 IGR=IGD,IGF
              FLXINT(IGC)=FLXINT(IGC)+RECTMP(IGR,1)
              OVERV(IGC)=OVERV(IGC)+RECTMP(IGR,1)/VELMTR(IGR)
 111        CONTINUE
            IF(FLXINT(IGC) .NE. 0.0) THEN
              OVERV(IGC)=OVERV(IGC)/FLXINT(IGC)
              DO 112 IGR=IGD,IGF
                RECTMP(IGR,1)=RECTMP(IGR,1)/FLXINT(IGC)
 112          CONTINUE
              FLXINT(IGC)=FLXINT(IGC)*VOLUME
            ENDIF
C----
C  LOOP OBER MTR GROUP ASSOCIATED WITH CPO GROUPS
C----
            DO 120 IGR=IGD,IGF
C----
C  READ CROSS SECTIONS
C----
              CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1 ,TKEY2 ,
     >                    NBE   )
              IF( NBE .NE. 3+NGMTR ) CALL XABORT(NAMSBR//
     >        ': CANNOT FIND '//TKEY1(1)//' '//TKEY2(1))
              READ(IFT16) RKEY1,RKEY2,NBE,
     >                    RECTMP(IGR,4),RECTMP(IGR,3),RECTMP(IGR,2),
     >                    (RECSCA(IGR,JGR),JGR=1,NGMTR)
              IF(IPRINT .GE. 100) THEN
                WRITE(IOUT,6101) TKEY2(1),IGR
                WRITE(IOUT,6110)
     >            RECTMP(IGR,4),RECTMP(IGR,3),RECTMP(IGR,2),
     >           (RECSCA(IGR,JGR),JGR=1,NGMTR)
              ENDIF
C----
C  ABSORPTION, NU-FISSION AND TRANSPORT SECTION CONDENSATION
C----
              RECXSV(IGC, 2)=RECXSV(IGC, 2)
     >                      +RECTMP(IGR,2)*RECTMP(IGR,1)
              RECXSV(IGC, 3)=RECXSV(IGC, 3)
     >                      +RECTMP(IGR,3)*RECTMP(IGR,1)
              RECXSV(IGC,15)=RECXSV(IGC,15)
     >                      +RECTMP(IGR,4)*RECTMP(IGR,1)
C----
C  SCATTERING SECTION CONDENSATION
C----
              JGF=0
              DO 121 JGC=1,NGCCPO
                JGD=JGF+1
                JGF=IFGMTR(JGC)
                DO 122 JGR=JGD,JGF
                  RECXSM(JGC,IGC,1)=RECXSM(JGC,IGC,1)
     >                             +RECSCA(IGR,JGR)*RECTMP(IGR,1)
                  RECXSV(IGC,21)=RECXSV(IGC,21)
     >                          +RECSCA(IGR,JGR)*RECTMP(IGR,1)
 122            CONTINUE
 121          CONTINUE
 120        CONTINUE
C----
C  TOTAL AND TRANSPORT CORRECTION
C----
            RECXSV(IGC,1)=RECXSV(IGC,15)+RECXSV(IGC,21)
            RECXSV(IGC,2)=RECXSV(IGC,1)-RECXSV(IGC,2)
 110      CONTINUE
          TKEY1(1)='REGION    '
          TKEY2(1)='DIFFUSION '
          CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1 ,TKEY2 ,
     >              NBE   )
          IF( NBE .EQ. 2*NGMTR ) THEN
            READ(IFT16) RKEY1,RKEY2,NBE,
     >                 (RECTMP(IGR,2),IGR=1,NGMTR),
     >                 (RECTMP(IGR,3),IGR=1,NGMTR)
            IF(IPRINT .GE. 100) THEN
              WRITE(IOUT,6100) TKEY2(1)
              WRITE(IOUT,6110) (RECTMP(IGR,2),IGR=1,NGMTR)
              WRITE(IOUT,6110) (RECTMP(IGR,3),IGR=1,NGMTR)
            ENDIF
C----
C  CONDENSE DIFFUSION COEFFICIENTS
C  COMPUTE STRD=1/3*DIFF
C----
            IGF=0
            DO 140 IGC=1,NGCCPO
              IGD=IGF+1
              IGF=IFGMTR(IGC)
              DO 141 IGR=IGD,IGF
                RECXSV(IGC,17)=RECXSV(IGC,17)+RECTMP(IGR,1)
     >             *(B2CRI(1)*RECTMP(IGR,2)+B2CRI(2)*RECTMP(IGR,3))
                RECXSV(IGC,18)=RECXSV(IGC,18)
     >                        +RECTMP(IGR,1)*RECTMP(IGR,2)
                RECXSV(IGC,19)=RECXSV(IGC,19)
     >                        +RECTMP(IGR,1)*RECTMP(IGR,2)
                RECXSV(IGC,20)=RECXSV(IGC,20)
     >                        +RECTMP(IGR,1)*RECTMP(IGR,3)
 141          CONTINUE
              IF(RECXSV(IGC,17) .EQ. 0.0 .OR.
     >           RECXSV(IGC,18) .EQ. 0.0 .OR.
     >           RECXSV(IGC,19) .EQ. 0.0 .OR.
     >           RECXSV(IGC,19) .EQ. 0.0 ) THEN
                RECXSV(IGC,17)=RECXSV(IGC,1)-RECXSV(IGC,2)
                RECXSV(IGC,18)=0.0
                RECXSV(IGC,19)=0.0
                RECXSV(IGC,20)=0.0
              ELSE
                RECXSV(IGC,17)=1.0/(3.0*RECXSV(IGC,17))
                RECXSV(IGC,18)=1.0/(3.0*RECXSV(IGC,18))
                RECXSV(IGC,19)=1.0/(3.0*RECXSV(IGC,19))
                RECXSV(IGC,20)=1.0/(3.0*RECXSV(IGC,20))
              ENDIF
 140        CONTINUE
          ELSE
            DO 142 IGC=1,NGCCPO
              RECXSV(IGC,17)=1.0/(3.0*(RECXSV(IGC,1)-RECXSV(IGC,2)))
              RECXSV(IGC,18)=RECXSV(IGC,17)
              RECXSV(IGC,19)=RECXSV(IGC,17)
              RECXSV(IGC,20)=RECXSV(IGC,17)
 142        CONTINUE
          ENDIF
          GO TO 105
        ELSE
          READ(IFT16) RKEY1,RKEY2,NBE
        ENDIF
 100  CONTINUE
 105  CONTINUE
C----
C  LOCATE NEXT CELLAV RECORD
C----
      TKEY1(2)='BEGIN     '
      TKEY2(2)='LEAKAGE   '
C----
C  READ FISSION SPECTRUM
C----
      TKEY1(1)='CELLAV    '
      TKEY2(1)='FISSPECT  '
      CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1 ,TKEY2 ,
     >            NBE   )
      IF( NBE .NE. NGMTR ) CALL XABORT(NAMSBR//
     >': CANNOT FIND '//TKEY1(1)//' '//TKEY2(1))
      READ(IFT16) RKEY1,RKEY2,NBE,(RECTMP(IGR,4),IGR=1,NGMTR)
C----
C  CONDENSE FISSION SPECTRUM OVER CPO GROUPS
C----
      IGF=0
      DO 150 IGC=1,NGCCPO
        IGD=IGF+1
        IGF=IFGMTR(IGC)
        DO 151 IGR=IGD,IGF
          RECXSV(IGC, 5)=RECXSV(IGC,5)+RECTMP(IGR,4)
 151    CONTINUE
 150  CONTINUE
C----
C  BURNUP INFORMATION
C----
      TKEY1(2)='MTR       '
      TKEY2(2)='FEWGROUPS '
      TKEY1(1)='CELLAV    '
      TKEY2(1)='AVG-ENERGY'
      CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1 ,TKEY2 ,
     >            NBE   )
      IF( NBE .EQ. 5 ) THEN
        READ(IFT16) RKEY1,RKEY2,NBE,RTIME,
     >              BRNTMP(3),BRNTMP(1),BRNTMP(2)
        IF(IPRINT .GE. 10) THEN
          WRITE(IOUT,6010) RTIME,BRNTMP(3),BRNTMP(1),BRNTMP(2)
        ENDIF
        BRNIRR(1)=BRNTMP(1)
        BRNIRR(2)=BRNTMP(2)
        BRNIRR(3)=WSMEV*BRNTMP(3)
      ENDIF
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6001)
      ENDIF
      RETURN
C----
C  PRINT FORMAT
C----
 6000 FORMAT(1X,5('*'),' OUTPUT FROM ',A6,1X,5('*'))
 6001 FORMAT(1X,30('*'))
 6010 FORMAT(6X,'BURNUP IRRADIATION '/1P,
     >       6X,'TIME    (DAYS)     = ',E10.3/
     >       6X,'ENERGY  (MEV)      = ',E10.3/
     >       6X,'BURNUP  (MWD/T)    = ',E10.3/
     >       6X,'IRRADIATION (N/KB) = ',E10.3)
 6100 FORMAT(6X,'CELLAV MAIN TRANSPORT GROUP ',A10)
 6101 FORMAT(6X,'CELLAV MAIN TRANSPORT GROUP ',A10,
     >       6X,'GROUP  =',I10)
 6110 FORMAT(1P,10(2X,E10.3))
      END
