*DECK PSPLEG
      SUBROUTINE PSPLEG(IPRINT,ISPSP ,ITYPE ,ICOLR ,NSUR  ,NVOL  ,
     >                  NAMLEG,NUNKNO,FLUX  ,NREGT ,
     >                  MATALB,KEYMRG,KEYFLX,COLREG)
C
C-------------------------    PSPLEG    -------------------------------
C
C 1- SUBROUTINE STATISTICS:
C     NAME     : PSPLEG
C     USE      : ASSOCIATE A COLOR TO A REGION AND PRINT LEGEND
C     MODIFIED : 99-03-15
C     AUTHOR   : G. MARLEAU
C
C 2- PARAMETERS:
C  INPUT
C     IPRINT : PRINT LEVEL                            I
C     ISPSP  : PSP FILE UNIT                          I
C     ITYPE  : TYPE OF GRAPHIC                        I
C              =  0 COLOR PER REGION NUMBER
C                 1 COLOR PER MATERIAL
C                 2 COLOR FOR FLUX (ONE GROUP)
C                 3 COLOR FOR FLUX (MULTIGROUP)
C                 4 COLOR PER MATERIAL FOR HOMOGENIZATION (HMIX)
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
C     NSUR    : -NUMBER OF OUTER SURFACE              I
C     NVOL    : MAXIMUM NUMBER OF REGIONS             I
C     NAMLEG  : LEGEND NAME                           C*24
C     NUNKNO  : NUMBER OF UNKNOW                      I
C     FLUX    : UNKNOWN VECTOR                        R(NUNKNO)
C     NREGT   : DIMENSION OF KEYFLX VECTOR            I
C     MATALB  : ALBEDO-MATERIAL OF REGIONS            I(NSUR:NVOL)
C     KEYMRG  : MERGE INDEX                           I(NSUR:NVOL)
C     KEYFLX  : FLUX LOCATION                         I(NREGT)
C     COLREG  : REGION COLOR                          I(4,NVOL)
C----------------------------------------------------------------------
C
      IMPLICIT         NONE
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      REAL             WLINE
      PARAMETER       (IOUT=6,WLINE=0.002,NAMSBR='PSPLEG')
C----
C  ROUTINE PARAMETERS
C----
      INTEGER          IPRINT,ISPSP,ITYPE,ICOLR,NSUR,NVOL,
     >                 NUNKNO,NREGT
      INTEGER          MATALB(NSUR:NVOL),KEYMRG(NSUR:NVOL),
     >                 KEYFLX(NVOL)
      REAL             FLUX(NUNKNO),COLREG(4,NVOL)
      REAL             COLTMP(4)
      CHARACTER        NAMLEG*24
C----
C  LOCAL PARAMETERS
C----
      CHARACTER        COLNAM*4,LEGTXT*48,FLXTXT*80
      INTEGER          MXMIX,MREG,IVOL,IMX,IRG,ICOLA,
     >                 ILEG,IFRM,MXCOL,ICOLF,IKEY
      INTEGER          KMX,ICT
      REAL             XYPOS(2),POSL,POSB,DELX,DELY,DELXC,DELYC,
     >                 XYPTS(2,4),FLXMIN,FLXMAX,DELFLX,COLFLX(4)
      INTEGER          KFS,KFR,KSS,KSR
C----
C  INITIALIZE LEGEND
C----
      KFS=0
      KFR=0
      KSS=0
      KSR=0
      ICOLA=ABS(ICOLR)
      IF(ICOLA .GT. 0) THEN
        KFS=1
        KSR=1
      ENDIF
      IF(ICOLA .GE. 2) THEN
        LEGTXT='Color by '//NAMLEG
      ELSE
        LEGTXT='Graylevel by '//NAMLEG
      ENDIF
      ILEG=1
      IF(IPRINT .LE. 0) THEN
        ILEG=0
      ENDIF
C----
C  GENERATE RANDOM COLOR
C  FOR RGB USE ALL THREE COLORS
C  FOR BW USE ONLY FIRST COLOR
C  SKIP FOR NONE
C----
      IF(ICOLA .GT. 0) THEN
        POSL=0.0
        POSB=10.0
        XYPOS(1)=POSL
        XYPOS(2)=POSB
        IF(ILEG .EQ. 1) THEN
          CALL PSTEXT(ISPSP,6,'Legend',
     >      XYPOS,0.1,0,0.0)
        ENDIF
        IF(ITYPE .EQ. 0) THEN
C----
C  COMPUTE NUMBER OF REGIONS AFTER MERGE
C----
          MREG=0
          DO 100 IVOL=1,NVOL
            MREG=MAX(MREG,KEYMRG(IVOL))
 100      CONTINUE
C----
C  GENERATE ONE COLOR PER REGION
C----
          POSB=POSB-0.2
          XYPOS(2)=POSB
          IF(ILEG .EQ. 1) THEN
            CALL PSTEXT(ISPSP,48,LEGTXT,XYPOS,0.1,0,0.0)
          ENDIF
          POSB=POSB-0.2
          IF(MREG .GT. 10000) THEN
            ILEG=0
          ENDIF
          DELX=0.2
          DELY=DELX/2.0
          DELXC=DELY
          DELYC=DELXC/4.0
          DO 110 IRG=1,MREG
            IFRM=0
            IF(MOD(IRG-1,30) .EQ. 0 .AND. ILEG .EQ. 1) THEN
              POSB=POSB-DELY
            ENDIF
            DO 111 IVOL=1,NVOL
              IF(KEYMRG(IVOL) .EQ. IRG) THEN
                CALL PSPCOL(ICOLA,MREG,IRG,COLREG(1,IVOL))
                IF(IFRM .EQ. 0 .AND. ILEG .EQ.1) THEN
                  IFRM=IFRM+1
                  POSL=MOD(IRG-1,30)*DELX
                  XYPTS(1,1)=POSL
                  XYPTS(2,1)=POSB
                  XYPTS(1,2)=POSL+DELX
                  XYPTS(2,2)=POSB
                  XYPTS(1,3)=POSL+DELX
                  XYPTS(2,3)=POSB+DELY
                  XYPTS(1,4)=POSL
                  XYPTS(2,4)=POSB+DELY
                  CALL PSDREG(ISPSP,4,XYPTS)
                  IF(ICOLA .GT. 0) THEN
                    CALL PSFILL(ISPSP,ICOLA,COLREG(1,IVOL),KFS,KFR)
                  ENDIF
                  CALL PSSTRK(ISPSP,WLINE,KSS,KSR)
                  WRITE(COLNAM,'(I4)') IRG
                  XYPOS(1)=POSL+DELXC
                  XYPOS(2)=POSB+DELYC
                  CALL PSTEXT(ISPSP,4,COLNAM,XYPOS,0.05,1,0.0)
                ENDIF
              ENDIF
 111        CONTINUE
 110      CONTINUE
        ELSE IF(ITYPE .EQ. 1 .OR. ITYPE .EQ. 4) THEN
C----
C  COMPUTE NUMBER OF MIXTURES
C----
          MXMIX=0
          DO 120 IVOL=1,NVOL
            MXMIX=MAX(MXMIX,MATALB(IVOL))
 120      CONTINUE
          POSB=POSB-0.2
          XYPOS(2)=POSB
          IF(ILEG .EQ. 1) THEN
            CALL PSTEXT(ISPSP,32,LEGTXT,XYPOS,0.1,0,0.0)
          ENDIF
          POSB=POSB-0.2
          IF(MXMIX .GT. 10000) THEN
            ILEG=0
          ENDIF
          KMX=0
          DELX=0.2
          DELY=DELX/2.0
          DELXC=DELY
          DELYC=DELXC/4.0
C----
C  GENERATE ONE COLOR PER MIXTURE
C----
          DO 130 IMX=0,MXMIX
            KMX=KMX+1
            IFRM=0
            IF(MOD(KMX-1,30).EQ.0 .AND. ILEG .EQ. 1) THEN
              POSB=POSB-DELY
            ENDIF
            CALL PSPCOL(ICOLA,MXMIX,IMX,COLTMP)
            IF (ILEG.EQ.1) THEN
              POSL=MOD(KMX-1,30)*DELX
              XYPTS(1,1)=POSL
              XYPTS(2,1)=POSB
              XYPTS(1,2)=POSL+DELX
              XYPTS(2,2)=POSB
              XYPTS(1,3)=POSL+DELX
              XYPTS(2,3)=POSB+DELY
              XYPTS(1,4)=POSL
              XYPTS(2,4)=POSB+DELY
              CALL PSDREG(ISPSP,4,XYPTS)
              IF(ICOLA .GT. 0) THEN
                 CALL PSFILL(ISPSP,ICOLA,COLTMP,KFS,KFR)
              ENDIF
              CALL PSSTRK(ISPSP,WLINE,KSS,KSR)
              WRITE(COLNAM,'(I4)') IMX
              XYPOS(1)=POSL+DELXC
              XYPOS(2)=POSB+DELYC
              CALL PSTEXT(ISPSP,4,COLNAM,XYPOS,0.05,1,0.0)
            ENDIF
C----
C  ASSOCIATE MIXTURE COLOR WITH REGION
C----
            DO 131 IVOL=1,NVOL
              IF(MATALB(IVOL) .EQ. IMX) THEN
                DO 132 ICT=1,4
                  COLREG(ICT,IVOL)=COLTMP(ICT)
 132            CONTINUE
              ENDIF
 131        CONTINUE
 130      CONTINUE
        ELSE IF(ITYPE .EQ. 2 .OR. ITYPE .EQ. 3 .OR.
     >          ITYPE .EQ. 5 .OR. ITYPE .EQ. 6) THEN
C----
C  COMPUTE NUMBER OF REGIONS AFTER MERGE
C----
          POSB=POSB-0.2
          XYPOS(2)=POSB
          IF(ILEG .EQ. 1) THEN
            CALL PSTEXT(ISPSP,32,LEGTXT,XYPOS,0.1,0,0.0)
          ENDIF
          POSB=POSB-0.2
C----
C  FIND MAXIMUM AND MINIMUM FLUX
C----
          FLXMAX=FLUX(KEYFLX(1))
          FLXMIN=FLUX(KEYFLX(1))
          DO 150 IRG=2,NREGT
            IKEY=KEYFLX(IRG)
            FLXMAX=MAX(FLXMAX,FLUX(IKEY))
            FLXMIN=MIN(FLXMIN,FLUX(IKEY))
 150      CONTINUE
          MXCOL=20
          DELFLX=(FLXMAX-FLXMIN)/REAL(MXCOL)
          WRITE(FLXTXT,5000) FLXMIN,DELFLX,FLXMIN,DELFLX
          XYPOS(2)=POSB
          IF(ILEG .EQ. 1) THEN
            CALL PSTEXT(ISPSP,80,FLXTXT,XYPOS,0.1,0,0.0)
          ENDIF
          POSB=POSB-0.2
          DELX=0.2
          DELY=DELX/2.0
          DELXC=DELY
          DELYC=DELXC/4.0
C----
C  GENERATE ONE COLOR PER FLUX LEVEL
C  COLOR I IS GIVEN BY:
C  I=MIN(INT((FLUX-FLXMIN)/DELFLX)+1,MXCOL)
C----
          POSB=POSB-DELY
          DO 160 ICOLF=1,MXCOL
            CALL PSPCOL(ICOLA,MXCOL,ICOLF,COLFLX)
            POSL=MOD(ICOLF-1,30)*DELX
            XYPTS(1,1)=POSL
            XYPTS(2,1)=POSB
            XYPTS(1,2)=POSL+DELX
            XYPTS(2,2)=POSB
            XYPTS(1,3)=POSL+DELX
            XYPTS(2,3)=POSB+DELY
            XYPTS(1,4)=POSL
            XYPTS(2,4)=POSB+DELY
            CALL PSDREG(ISPSP,4,XYPTS)
            IF(ICOLA .GT. 0) THEN
              CALL PSFILL(ISPSP,ICOLA,COLFLX,KFS,KFR)
            ENDIF
            CALL PSSTRK(ISPSP,WLINE,KSS,KSR)
            WRITE(COLNAM,'(I4)') ICOLF
            XYPOS(1)=POSL+DELXC
            XYPOS(2)=POSB+DELYC
            CALL PSTEXT(ISPSP,4,COLNAM,XYPOS,0.05,1,0.0)
 160      CONTINUE
          DO 170 IRG=1,NREGT
            IKEY=KEYFLX(IRG)
            ICOLF=INT((FLUX(IKEY)-FLXMIN)/DELFLX)+1
            ICOLF=MIN(ICOLF,MXCOL)
            DO 171 IVOL=1,NVOL
              IF(KEYMRG(IVOL) .EQ. IRG) THEN
                CALL PSPCOL(ICOLA,MXCOL,ICOLF,COLREG(1,IVOL))
              ENDIF
 171        CONTINUE
 170      CONTINUE
        ENDIF
      ENDIF
      RETURN
C----
C  FORMAT
C----
 5000 FORMAT(1P,E9.2,'+(i-1)*',E9.2,
     >       ' < Flux(i) <= ',E9.2,'+i*',E9.2)
      END
