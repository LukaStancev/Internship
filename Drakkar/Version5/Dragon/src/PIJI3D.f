*DECK PIJI3D
      SUBROUTINE PIJI3D(NREG,NSOUT,NSLINE,NCOR,NSBG,
     >                  SWVOID,SIGTAL,NPSYS,WEIGHT,
     >                  SEGLEN,NRSEG,
     >                  STAYIN,GOSOUT,DPR)
C
C-------------------------    PIJI3D    -------------------------------
C
C 1- SUBROUTINE STATISTICS:
C     NAME     : PIJI3D
C     LEVEL    : 2 (CALLED BY 'EXCELP')
C     USE      : INTEGRATION FOR GENERAL 3D ISOTROPIC B.C. TRACKING
C     MODIFIED : 91/07/12 (R.R.)
C     AUTHOR   : R. ROY
C
C 2- PARAMETERS:
C  INPUT
C     NREG    : TOTAL NUMBER OF REGIONS                I
C     NSOUT   : NUMBER OF OUTER SURFACE                I
C     NSLINE  : number of segemnts on line             
C     NCOR    : MAXIMUM NUMBER OF CORNERS              I
C     NSBG    : NUMBER OF SUBGROUP                     I
C     SWVOID  : FLAG TO INDICATE IF THERE ARE VOIDS    L
C     SIGTAL  : ALBEDO-CROSS SECTION VECTOR            R(-NSOUT:NREG,
C               IS=-NSOUT,-1; SIGTAL(IS)=ALBEDO(-IS)     NSBG)
C               IV=1,NREG ;   SIGTAL(IV)=SIGT(IV)
C     NPSYS   : NON-CONVERGED ENERGY GROUP INDICES.    I(NSBG)
C     SEGLEN  : LENGTH OF TRACK              D(NSLINE)
C     NRSEG   : REGION CROSSED BY TRACK      I(NSLINE)
C  OUTPUT
C     DPR     : COLLISION PROBABILITIES      D(-NSOUT:NREG,
C                                              -NSOUT:NREG,NSBG)
C               VALUES ARE COMPUTED ONLY FOR ONE DIRECTION P(I->J)
C  WORK
C     STAYIN  : STAY-IN ZONE PROBABILITY     D(NSLINE)
C     GOSOUT  : GOES-OUT ZONE PROBABILITY    D(NSLINE)
C
C-------------------------    PIJI3D    -------------------------------
C
      IMPLICIT         NONE
C----
C VARIABLES
C----
      INTEGER          NREG,NSOUT,NSLINE,NCOR,NSBG
      INTEGER          NRSEG(NSLINE),NPSYS(NSBG)
      LOGICAL          SWVOID
      REAL             SIGTAL(-NSOUT:NREG,NSBG)
      DOUBLE PRECISION WEIGHT,SEGLEN(NSLINE),STAYIN(NSLINE),
     >                 GOSOUT(NSLINE)
      DOUBLE PRECISION DPR(-NSOUT:NREG,-NSOUT:NREG,NSBG)
C----
C  Local variables
C----
      INTEGER          IL,JL,NOIL,ISBG
      REAL             ZERO, ONE, HALF
      DOUBLE PRECISION XSIL, PRODUC, DSCBEG, DSCEND, ZCOR, ZCOR2
      INTEGER          ICSEG,JCSEG,ISD,ISF
      PARAMETER       (ZERO=0.0E0, ONE=1.0E0, HALF=0.5E0 )
      REAL             SIXT,CUTEXP
      PARAMETER       (SIXT=HALF/3.0,CUTEXP=0.02)
      DOUBLE PRECISION EXSIL,XSIL2
      DO 2001 ISBG=1,NSBG
      IF(NPSYS(ISBG).EQ.0) GO TO 2001
*----
*  Process track required
*----
      IF( NCOR.EQ.1 )THEN
C
C1)   ONLY ONE EXTERNAL SURFACE AT END --------------------------------
      ISD=NRSEG(1)
      ISF=NRSEG(NSLINE)
      IF( SWVOID )THEN
         PRODUC= WEIGHT
C        PII CALCULATION AND ESCAPE
         DO 40 IL = 1,NSLINE-2
            ICSEG=IL+1
            NOIL  = NRSEG(ICSEG)
            XSIL  = SIGTAL(NOIL,ISBG)*SEGLEN(ICSEG)
            IF( XSIL.EQ.ZERO )THEN
               GOSOUT(IL)= ONE
               STAYIN(IL)= SEGLEN(ICSEG)
               DPR(NOIL,NOIL,ISBG)= DPR(NOIL,NOIL,ISBG)
     >                       + HALF*WEIGHT*SEGLEN(ICSEG)*SEGLEN(ICSEG)
            ELSE IF(XSIL .LT. CUTEXP) THEN
               XSIL2=XSIL*XSIL
               EXSIL=XSIL2*(HALF-SIXT*XSIL)
               STAYIN(IL)=XSIL-EXSIL
               GOSOUT(IL)=ONE-STAYIN(IL)
               PRODUC= PRODUC * GOSOUT(IL)
               DPR(NOIL,NOIL,ISBG)= DPR(NOIL,NOIL,ISBG)
     >          + WEIGHT*EXSIL
            ELSE
               EXSIL=EXP( - XSIL )
               STAYIN(IL)= ONE - EXSIL
               GOSOUT(IL)= EXSIL
               PRODUC= PRODUC * GOSOUT(IL)
               DPR(NOIL,NOIL,ISBG)= DPR(NOIL,NOIL,ISBG)
     >                       + WEIGHT*(XSIL-STAYIN(IL))
            ENDIF
   40    CONTINUE
C        PIJ CALCULATION
         DSCBEG= WEIGHT
         DO 60 IL = 1, NSLINE-2
            ICSEG=IL+1
            NOIL  = NRSEG(ICSEG)
            DSCEND= WEIGHT*STAYIN(IL)
            DO 50 JL  = IL+1, NSLINE-2
               JCSEG=JL+1
               DPR(NRSEG(JCSEG),NOIL,ISBG)= 
     >         DPR(NRSEG(JCSEG),NOIL,ISBG)+ STAYIN(JL)*DSCEND
               DSCEND= DSCEND*GOSOUT(JL)
   50       CONTINUE
C           PIS CALCULATION
            DPR(ISD,NOIL,ISBG)= DPR(ISD,NOIL,ISBG)+DSCBEG*STAYIN(IL)
            DPR(ISF,NOIL,ISBG)= DPR(ISF,NOIL,ISBG)+DSCEND
            DSCBEG= DSCBEG * GOSOUT(IL)
   60    CONTINUE
C        PSS CALCULATION
         DPR(ISD,ISF,ISBG)= DPR(ISD,ISF,ISBG) + PRODUC
      ELSE
C
C1.2) NO VOID REGION
         PRODUC= WEIGHT
C        PII CALCULATION AND ESCAPE
         DO 140 IL = 1,NSLINE-2
            ICSEG=IL+1
            NOIL  = NRSEG(ICSEG)
            XSIL  = SIGTAL(NOIL,ISBG)*SEGLEN(ICSEG)
            IF(XSIL .LT. CUTEXP) THEN
               XSIL2=XSIL*XSIL
               EXSIL=XSIL2*(HALF-SIXT*XSIL)
               STAYIN(IL)=XSIL-EXSIL
               GOSOUT(IL)=ONE-STAYIN(IL)
               PRODUC= PRODUC * GOSOUT(IL)
               DPR(NOIL,NOIL,ISBG)= DPR(NOIL,NOIL,ISBG)
     >          + WEIGHT*EXSIL
            ELSE
               EXSIL=EXP( - XSIL )
               STAYIN(IL)= ONE - EXSIL
               GOSOUT(IL)= EXSIL
               PRODUC= PRODUC * GOSOUT(IL)
               DPR(NOIL,NOIL,ISBG)= DPR(NOIL,NOIL,ISBG)
     >                       + WEIGHT*(XSIL-STAYIN(IL))
            ENDIF
  140    CONTINUE
C        PIJ CALCULATION
         DSCBEG= WEIGHT
         DO 160 IL = 1, NSLINE-2
            ICSEG=IL+1
            NOIL  = NRSEG(ICSEG)
            DSCEND= WEIGHT*STAYIN(IL)
            DO 150 JL  = IL+1, NSLINE-2
               JCSEG=JL+1
               DPR(NRSEG(JCSEG),NOIL,ISBG)= 
     >         DPR(NRSEG(JCSEG),NOIL,ISBG)+ STAYIN(JL)*DSCEND
               DSCEND= DSCEND*GOSOUT(JL)
  150       CONTINUE
C           PIS CALCULATION
            DPR(ISD,NOIL,ISBG)= DPR(ISD,NOIL,ISBG)+DSCBEG*STAYIN(IL)
            DPR(ISF,NOIL,ISBG)= DPR(ISF,NOIL,ISBG)+DSCEND
            DSCBEG= DSCBEG * GOSOUT(IL)
  160    CONTINUE
C        PSS CALCULATION
         DPR(ISD,ISF,ISBG)= DPR(ISD,ISF,ISBG) + PRODUC
      ENDIF
      ELSE
C
C2)   MORE THAN ONE SURFACE PER LINE ----------------------------------
      ZCOR= 1./FLOAT(NCOR)
      ZCOR2= ZCOR*ZCOR
      IF( SWVOID )THEN
C
C2.1) VOIDS ARE POSSIBLE
         PRODUC= WEIGHT*ZCOR2
C        PII CALCULATION AND ESCAPE
         DO 240 IL = 1,NSLINE-2*NCOR
            ICSEG=IL+NCOR
            NOIL  = NRSEG(ICSEG)
            XSIL  = SIGTAL(NOIL,ISBG)*SEGLEN(ICSEG)
            IF( XSIL.EQ.ZERO )THEN
               GOSOUT(IL)= ONE
               STAYIN(IL)= SEGLEN(ICSEG)
               DPR(NOIL,NOIL,ISBG)= DPR(NOIL,NOIL,ISBG)
     >                       + HALF*WEIGHT*SEGLEN(ICSEG)*SEGLEN(ICSEG)
            ELSE IF(XSIL .LT. CUTEXP) THEN
               XSIL2=XSIL*XSIL
               STAYIN(IL)=XSIL-XSIL2*(HALF-SIXT*XSIL)
               GOSOUT(IL)=ONE-STAYIN(IL)
               PRODUC= PRODUC * GOSOUT(IL)
               DPR(NOIL,NOIL,ISBG)= DPR(NOIL,NOIL,ISBG)
     >          + WEIGHT*XSIL2*(HALF-SIXT*XSIL)
            ELSE
               GOSOUT(IL)= EXP( - XSIL )
               STAYIN(IL)= (ONE - GOSOUT(IL))
               PRODUC= PRODUC * GOSOUT(IL)
               DPR(NOIL,NOIL,ISBG)= DPR(NOIL,NOIL,ISBG)
     >                       + WEIGHT*(XSIL-STAYIN(IL))
            ENDIF
  240    CONTINUE
C        PIJ CALCULATION
         DSCBEG= WEIGHT*ZCOR
         DO 260 IL = 1, NSLINE-2*NCOR
            ICSEG=IL+NCOR
            NOIL  = NRSEG(ICSEG)
            DSCEND= WEIGHT*STAYIN(IL)
            DO 250 JL  = IL+1, NSLINE-2*NCOR
               JCSEG=JL+NCOR
               DPR(NRSEG(JCSEG),NOIL,ISBG)= 
     >         DPR(NRSEG(JCSEG),NOIL,ISBG)+ STAYIN(JL)*DSCEND
               DSCEND= DSCEND*GOSOUT(JL)
  250       CONTINUE
C           PIS CALCULATION
            DO 261 JL = 1, NCOR
               ISD=NRSEG(JL)
               ISF=NRSEG(NSLINE-NCOR+JL)
               DPR(ISD,NOIL,ISBG)= DPR(ISD,NOIL,ISBG)+DSCBEG*STAYIN(IL)
               DPR(ISF,NOIL,ISBG)= DPR(ISF,NOIL,ISBG)+DSCEND*ZCOR
  261       CONTINUE
            DSCBEG= DSCBEG*GOSOUT(IL)
  260    CONTINUE
C        PSS CALCULATION
         DO 265 IL = 1, NCOR
         ISD=NRSEG(IL)
         DO 265 JL = 1, NCOR
            ISF=NRSEG(NSLINE-NCOR+JL)
            DPR(ISD,ISF,ISBG)= DPR(ISD,ISF,ISBG) + PRODUC
  265    CONTINUE
      ELSE
C
C2.2) NO VOID REGION
         PRODUC= WEIGHT*ZCOR2
C        PII CALCULATION AND ESCAPE
         DO 340 IL = 1,NSLINE-2*NCOR
            ICSEG=IL+NCOR
            NOIL  = NRSEG(ICSEG)
            XSIL  = SIGTAL(NOIL,ISBG)*SEGLEN(ICSEG)
            IF(XSIL .LT. CUTEXP) THEN
               XSIL2=XSIL*XSIL
               STAYIN(IL)=XSIL-XSIL2*(HALF-SIXT*XSIL)
               GOSOUT(IL)=ONE-STAYIN(IL)
               PRODUC= PRODUC * GOSOUT(IL)
               DPR(NOIL,NOIL,ISBG)= DPR(NOIL,NOIL,ISBG)
     >          + WEIGHT*XSIL2*(HALF-SIXT*XSIL)
            ELSE
               GOSOUT(IL)= EXP( - XSIL )
               STAYIN(IL)= (ONE - GOSOUT(IL))
               PRODUC= PRODUC * GOSOUT(IL)
               DPR(NOIL,NOIL,ISBG)= DPR(NOIL,NOIL,ISBG)
     >                       + WEIGHT*(XSIL-STAYIN(IL))
            ENDIF
  340    CONTINUE
C        PIJ CALCULATION
         DSCBEG= WEIGHT*ZCOR
         DO 360 IL = 1, NSLINE-2*NCOR
            ICSEG=IL+NCOR
            NOIL  = NRSEG(ICSEG)
            DSCEND= WEIGHT*STAYIN(IL)
            DO 350 JL  = IL+1, NSLINE-2*NCOR
               JCSEG=JL+NCOR
               DPR(NRSEG(JCSEG),NOIL,ISBG)=
     >         DPR(NRSEG(JCSEG),NOIL,ISBG)+ STAYIN(JL)*DSCEND
               DSCEND= DSCEND*GOSOUT(JL)
  350       CONTINUE
C           PIS CALCULATION
            DO 361 JL = 1, NCOR
               ISD=NRSEG(JL)
               ISF=NRSEG(NSLINE-NCOR+JL)
               DPR(ISD,NOIL,ISBG)= DPR(ISD,NOIL,ISBG)+DSCBEG*STAYIN(IL)
               DPR(ISF,NOIL,ISBG)= DPR(ISF,NOIL,ISBG)+DSCEND*ZCOR
  361       CONTINUE
            DSCBEG= DSCBEG * GOSOUT(IL)
  360    CONTINUE
C        PSS CALCULATION
         DO 365 IL = 1, NCOR
         ISD=NRSEG(IL)
         DO 365 JL = 1, NCOR
            ISF=NRSEG(NSLINE-NCOR+JL)
            DPR(ISD,ISF,ISBG)= DPR(ISD,ISF,ISBG) + PRODUC
  365    CONTINUE
      ENDIF
      ENDIF
 2001 CONTINUE
      RETURN
      END
