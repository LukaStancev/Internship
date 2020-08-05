*DECK T16FLX
      SUBROUTINE T16FLX(IFT16 ,IPRINT,NGCCPO,NGMTR ,NMATZ ,MTRMSH,
     >                  IFGMTR,VELMTR,MXFGET,IMIREG,VOLUME,B2CRI ,
     >                  FLXINT,FLXDIS,OVERV ,KMSPEC,MATMSH,VQLE  ,
     >                  PHI   )
*
*----
*  1- PROGRAMME STATISTICS:
*      NAME     : T16FLX
*
*Purpose:
*  READ MAIN TRANSPORT GROUP FLUX AND COMPUTE
*                 INTEGRATED FLUX, FLUX DISADVANTAGE FACTOR
*                 AND I/V CROSS SECTION
*
*Author(s): 
* G.MARLEAU
*
*      CREATED  : 1999/10/21
*      REF      : IGE-244 REV.1
*
*      MODIFICATION LOG
*      --------------------------------------------------------------
*      | DATE AND INITIALS  | MOTIVATIONS
*      --------------------------------------------------------------
*      | 1999/10/21 G.M.    | READ MAIN TRANSPORT GROUP FLUX
*      |                    | COMPUTE INTEGRATED FLUX,
*      |                    | FLUX DISADVANTAGE FACTOR
*      |                    | AND I/V CROSS SECTION
*      --------------------------------------------------------------
*
*  2- ROUTINE PARAMETERS:
*Parameters: input
* IFT16   TAPE16 FILE UNIT                         I
* IPRINT  PRINT LEVEL                              I
*         =   0 NO PRINT
*         >=  1 PRINT PROCESSING OPTIONS READ
* NGCCPO  NUMBER OF FINAL CONDENSED GROUPS         I
* NGMTR   NUMBER OF MAIN TRANSPORT GROUP           I
* NMATZ   NUMBER OF MIXTURES                       I
* MTRMSH  NUMBER OF MAIN TRANSPORT MESH POINTS     I
* IFGMTR  CPO FEW GROUP IDENTIFIER                 I(NGCCPO)
*         WITH RESPECT TO MTR GROUPS
* VELMTR  AVERAGE VELOCITY IN MAIN GROUPS          R(NGMTR)
* MXFGET  MAXIMUM DIMENSION OF FLUX VECTOR         I
* IMIREG  MIXTURE UPDATE IDENTIFIER                I
*         =  0 DO NOT UPDATE
*         = -1 UPDATE USING CELLAV INFORMATION
*         >  0 UPDATE USING SPECIFIED REGION NUMBER
*
*Parameters: output
* VOLUME  TOTAL VOLUME                             R
* B2CRI   CRITICAL BUCKLINGS                       R(3)
* FLXINT  VOLUME INTEGRATED FLUX                   R(NGCCPO)
* FLXDIS  FLUX DISADVANTAGE FACTOR                 R(NGCCPO)
* OVERV   1/V CROSS SECTION                        R(NGCCPO)
*
*Parameters: work
* KMSPEC  MATERIAL TYPE                            I(NMATZ)
* MATMSH  MATERIAL AT EACH MESH POINT              I(MTRMSH)
* VQLE    VOLUME OF EACH MESH POINT                R(MTRMSH)
* PHI     MULTIGROUP FLUX AT EACH MESH POINT       R(MXFGET)
*
*  3- ROUTINES CALLED
*    SPECIFIC T16CPO ROUTINES
*      T16FND : FIND A TAPE16 RECORD
*               EQUIVALENT TO FIND FUNCTION
*               IN APPENDIX E OF EACL RC-1176
*    UTILITIES ROUTINES
*      XABORT : ABORT ROUTINE
*      XDRSET : VECTOR INITIALIZATION ROUTINE
*
*----
*
      IMPLICIT         NONE
      INTEGER          IFT16,IPRINT,NGCCPO,NGMTR,NMATZ,MTRMSH,
     >                 MXFGET,IMIREG
      INTEGER          IFGMTR(NGCCPO),KMSPEC(NMATZ),MATMSH(MTRMSH)
      REAL             VELMTR(NGMTR),VOLUME,B2CRI(3),FLXINT(NGCCPO),
     >                 FLXDIS(NGCCPO),OVERV(NGCCPO),VQLE(MTRMSH),
     >                 PHI(MXFGET)
C----
C  T16 PARAMETERS
C----
      INTEGER          MAXKEY
      PARAMETER       (MAXKEY=3)
      CHARACTER        TKEY1(MAXKEY)*10,TKEY2(MAXKEY)*10,
     >                 RKEY1*10,RKEY2*10
      INTEGER          NKEY,IOPT,NBE,NID,IR
      REAL             RID
C----
C  LOCAL VARIABLES
C----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='T16FLX')
      INTEGER          IGR,IGC,IGD,IGF,IMIX,ITRFL,IBUCK
      REAL             B2INI(3)
C----
C  SET END RECORDS FOR THIS SEARCH
C----
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
      IOPT=0
      CALL XDRSET(FLXINT,NGCCPO              ,0.0)
      CALL XDRSET(FLXDIS,NGCCPO              ,0.0)
      CALL XDRSET(OVERV ,NGCCPO              ,0.0)
C----
C  MTRFLX RECORDS
C----
      NKEY=2
      TKEY1(2)='REGION    '
      TKEY2(2)='DESCRIPTON'
      TKEY1(1)='MTRFLX    '
      TKEY2(1)='FLUX      '
      CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1 ,TKEY2 ,
     >            NBE   )
      ITRFL=0
      IF(NBE .GT. 0 ) THEN
        ITRFL=1
      ELSE IF( NBE .LT. -1 ) THEN
        READ(IFT16) RKEY1,RKEY2,NBE
        IF(IMIREG .GT. 0) THEN
          TKEY1(2)='CELLAV    '
          TKEY2(2)='NGROUPS   '
          TKEY1(1)='REGION    '
          TKEY2(1)='FLUX      '
          DO 100 IR=1,IMIREG-1
            CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1 ,TKEY2 ,
     >                  NBE   )
            IF( NBE .LE. 0 ) CALL XABORT(NAMSBR//
     >      ': REGION FLUX NOT AVAILABLE')
            READ(IFT16) RKEY1,RKEY2
 100      CONTINUE
          CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1 ,TKEY2 ,
     >                NBE   )
          IF( NBE .GT. 0 ) ITRFL=2
        ELSE IF(IMIREG .LT. 0) THEN
          TKEY1(2)='CELLAV    '
          TKEY2(2)='K         '
          TKEY1(1)='CELLAV    '
          TKEY2(1)='FLUX      '
          CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1 ,TKEY2 ,
     >                NBE   )
          IF( NBE .GT. 0 ) ITRFL=3
        ENDIF
      ENDIF
      IF( ITRFL .EQ. 0 ) THEN
         CALL XABORT(NAMSBR//
     >   ': CANNOT FIND '//TKEY1(1)//' '//TKEY2(1)//' OR '//
     >   TKEY1(2)//' '//TKEY2(2))
      ELSE IF(ITRFL .EQ. 1) THEN
C----
C  USE MTRFLX
C  1) CONDENSE AND HOMOGENIZE FLUX
C  2) COMPUTE FLUX DISADVANTAGE FACTOR
C  3) COMPUTE VOLUME
C  4) COMPUTE OVERV
C----
        READ(IFT16) RKEY1,RKEY2,NBE,NID,NID,
     >   (MATMSH(IR),VQLE(IR),
     >   (PHI(IGR+(IR-1)*NGMTR),IGR=1,NGMTR),IR=1,MTRMSH)
        VOLUME=0.0
        DO 110 IR=1,MTRMSH
          IGF=0
          VOLUME=VOLUME+VQLE(IR)
          IF(IPRINT .GE. 100) THEN
            WRITE(IOUT,6100) IR,VQLE(IR)
            WRITE(IOUT,6110)(PHI(IGR+(IR-1)*NGMTR),IGR=1,NGMTR)
          ENDIF
          DO 111 IGC=1,NGCCPO
            IGD=IGF+1
            IGF=IFGMTR(IGC)
            DO 112 IGR=IGD,IGF
              FLXINT(IGC)=FLXINT(IGC)+PHI(IGR+(IR-1)*NGMTR)*VQLE(IR)
              OVERV(IGC)=OVERV(IGC)
     >                  +PHI(IGR+(IR-1)*NGMTR)*VQLE(IR)/VELMTR(IGR)
 112        CONTINUE
 111      CONTINUE
          IMIX=MATMSH(IR)
          IF(KMSPEC(IMIX) .EQ. 1) THEN
            IGF=0
            DO 113 IGC=1,NGCCPO
              IGD=IGF+1
              IGF=IFGMTR(IGC)
              DO 114 IGR=IGD,IGF
                FLXDIS(IGC)=FLXDIS(IGC)+PHI(IGR+(IR-1)*NGMTR)*VQLE(IR)
 114          CONTINUE
 113        CONTINUE
          ENDIF
 110    CONTINUE
      ELSE IF(ITRFL .EQ. 2) THEN
C----
C  USE REGION FLUX
C  1) CONDENSE AND HOMOGENIZE FLUX
C  2) COMPUTE FLUX DISADVANTAGE FACTOR
C  4) COMPUTE OVERV
C----
        READ(IFT16) RKEY1,RKEY2,NBE,NID,NID,VOLUME,
     >    (PHI(IGR),IGR=1,NGMTR)
        IR=IMIREG
        IGF=0
        IF(IPRINT .GE. 100) THEN
          WRITE(IOUT,6100) IR,VOLUME
          WRITE(IOUT,6110)(PHI(IGR),IGR=1,NGMTR)
        ENDIF
        DO 120 IGC=1,NGCCPO
          IGD=IGF+1
          IGF=IFGMTR(IGC)
          DO 121 IGR=IGD,IGF
            FLXINT(IGC)=FLXINT(IGC)+PHI(IGR)*VOLUME
            OVERV(IGC)=OVERV(IGC)
     >              +(PHI(IGR)*VOLUME)/VELMTR(IGR)
 121      CONTINUE
          FLXDIS(IGC)=FLXINT(IGC)
 120    CONTINUE
      ELSE
C----
C  USE CELLAV FLUX
C  1) CONDENSE AND HOMOGENIZE FLUX
C  2) COMPUTE FLUX DISADVANTAGE FACTOR
C  3) COMPUTE VOLUME
C  4) COMPUTE OVERV
C----
        VOLUME=1.0
        READ(IFT16) RKEY1,RKEY2,NBE,
     >   (PHI(IGR),IGR=1,NGMTR)
        IF(IPRINT .GE. 100) THEN
          WRITE(IOUT,6101)
          WRITE(IOUT,6110)(PHI(IGR),IGR=1,NGMTR)
        ENDIF
        IGF=0
        DO 130 IGC=1,NGCCPO
          IGD=IGF+1
          IGF=IFGMTR(IGC)
          DO 131 IGR=IGD,IGF
            FLXINT(IGC)=FLXINT(IGC)+PHI(IGR)
            OVERV(IGC)=OVERV(IGC)+PHI(IGR)/VELMTR(IGR)
 131      CONTINUE
          FLXDIS(IGC)=FLXINT(IGC)
 130    CONTINUE
      ENDIF
      DO 140 IGC=1,NGCCPO
        FLXDIS(IGC)=FLXDIS(IGC)/FLXINT(IGC)
        OVERV(IGC)=OVERV(IGC)/FLXINT(IGC)
 140  CONTINUE
C----
C  RADIAL AND AXIAL DIFFUSION COEFFICIENTS
C  AND BUCKLING
C----
      TKEY1(2)='CELLAV    '
      TKEY2(2)='K         '
      TKEY1(1)='CELLAV    '
      TKEY2(1)='DIFFUSION '
      CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1 ,TKEY2 ,
     >            NBE   )
      IF( NBE .NE. 5*NGMTR+5 ) CALL XABORT(NAMSBR//
     >': CANNOT FIND '//TKEY1(1)//' '//TKEY2(1))
      READ(IFT16) RKEY1,RKEY2,NBE,(NID,IR=1,3),
     >           (RID,IGR=1,NGMTR),
     >           (RID,IGR=1,NGMTR),
     >           (RID,IGR=1,NGMTR),
     >           (B2INI(IR),IR=1,2)
      TKEY1(1)='CELLAV    '
      TKEY2(1)='CRITICALB '
      CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1 ,TKEY2 ,
     >            NBE   )
      IF( NBE .EQ. 2*NGMTR+4 ) THEN
        READ(IFT16) RKEY1,RKEY2,NBE,IBUCK,
     >             (B2CRI(IR),IR=1,3)
        IF(IBUCK .EQ. 2) THEN
          B2CRI(3)=B2INI(1)+B2CRI(2)
          B2CRI(1)=B2INI(1)/B2CRI(3)
          B2CRI(2)=B2CRI(2)/B2CRI(3)
        ELSE IF(IBUCK .EQ. 3) THEN
          B2CRI(3)=B2CRI(1)+B2INI(2)
          B2CRI(1)=B2CRI(1)/B2CRI(3)
          B2CRI(2)=B2INI(2)/B2CRI(3)
        ELSE
          B2CRI(1)=B2CRI(1)/B2CRI(3)
          B2CRI(2)=B2CRI(2)/B2CRI(3)
        ENDIF
      ELSE IF(NBE .EQ. -2) THEN
        B2CRI(3)=B2INI(1)+B2INI(2)
        B2CRI(1)=B2INI(1)/B2CRI(3)
        B2CRI(2)=B2INI(2)/B2CRI(3)
      ELSE
        CALL XABORT(NAMSBR//
     >': CANNOT FIND '//TKEY1(2)//' '//TKEY2(2))
      ENDIF
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6010) (FLXINT(IGC),IGC=1,NGCCPO)
        WRITE(IOUT,6011) (FLXDIS(IGC),IGC=1,NGCCPO)
        WRITE(IOUT,6012) (OVERV(IGC),IGC=1,NGCCPO)
        WRITE(IOUT,6013) (B2CRI(IR),IR=1,3)
        WRITE(IOUT,6001)
      ENDIF
      RETURN
C----
C  PRINT FORMAT
C----
 6000 FORMAT(1X,5('*'),' OUTPUT FROM ',A6,1X,5('*'))
 6001 FORMAT(1X,30('*'))
 6010 FORMAT(6X,'INTEGRATED FLUXES'/
     >1P,10(2X,E10.3))
 6011 FORMAT(6X,'FLUX DISADVANTAGE FACTORS'/
     >1P,10(2X,E10.3))
 6012 FORMAT(6X,'1/V '/
     >1P,10(2X,E10.3))
 6013 FORMAT(6X,'CRITICAL BUCKLINGS'/
     >1P,3(2X,E10.3))
 6100 FORMAT(6X,'MAIN TRANSPORT GROUP FLUX IN REGION  = ',I10,
     >       5X,'OF VOLUME = ',1P,E10.3)
 6101 FORMAT(6X,'CELLAV MAIN TRANSPORT GROUP FLUX  ')
 6110 FORMAT(1P,10(2X,E10.3))
      END
