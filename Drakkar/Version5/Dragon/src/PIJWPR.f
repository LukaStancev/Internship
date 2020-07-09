*DECK PIJWPR
      SUBROUTINE PIJWPR(LOPT,NREG,NSOUT,SIGTAL,PROB,SIGVOL,MSYM)
C
C-------------------------    PIJWPR    -------------------------------
C
C 1- SUBROUTINE STATISTICS:
C     NAME     : PIJWPR
C     LEVEL    : 2 (CALLED BY 'EXCELP')
C     USE      : PRINT-OUT FOR PROBABILITY MATRICES
C     MODIFIED : 94-05-25 (I.P.)
C     MODIFIED : 91-07-12 (R.R.)
C     AUTHOR   : R. ROY
C
C 2- PARAMETERS:
C  INPUT
C     LOPT    : PRINT-OUT FORM                         I
C               LOPT.LE.0
C                PRINT ALL  PSS/PVS/PVV
C               LOPT.GT.0
C                PRINT ONLY PVV
C     MSYM    : 1  SYMMETRIC TRINNG. MATRIX
C               0  NON-SYMMETRIC FULL MATRIX
C     NREG    : TOTAL NUMBER OF REGIONS                I
C     NSOUT   : NUMBER OF OUTER SURFACE                I
C     SIGTAL  : ALBEDO-CROSS SECTION VECTOR            R(-NSOUT:NREG)
C               IS=-NSOUT,-1; SIGTAL(IS)=ALBEDO(-IS)
C               IV=1,NREG ;   SIGTAL(IV)=SIGT(IV)
C     PROB    : COMPRESS PROBABILITY MATRIX            D(NPLEN)
C               NPLEN=(NREG+NSOUT+2)*(NREG+NSOUT+1)/2
C               IND(I,J)=MAX(I+NSOUT+1,J+NSOUT+1)
C                       *(MAX(I+NSOUT+1,J+NSOUT+1)-1)/2
C                       +MIN(I+NSOUT+1,J+NSOUT+1)
C
C         JS=-NSOUT,-1;               I=IND(0,JS)
C           PROB(I)=VOLSUR(JS)
C         IV=1,NREG;                  I=IND(IV,0)
C           SIGT(IV).GT.0.0
C             PROB(I)=SIGT(IV)*VOLSUR(IV)
C           SIGT(IV).EQ.0.0
C             PROB(I)=VOLSUR(IV)
C
C         IS=-NSOUT,-1; JS=-NSOUT,IS; I=IND(IS,JS)
C           PROB(I)=VOLSUR(IS)*PSS(IS,JS)
C         IV=1,NREG; JS=-NSOUT,-1;    I=IND(IV,JS)
C           SIGT(IV).GT.0.0
C             PROB(I)=SIGT(IV)*VOLSUR(IV)*PVS(IV,JS)
C           SIGT(IV).EQ.0.0
C             PROB(I)=VOLSUR(IV)*PVS(IV,JS)
C         IV=1,NREG; JV=1,IV;         I=IND(IV,JV)
C           SIGT(IV).GT.0.0 AND SIGT(JV).GT.0.0
C             PROB(I)=SIGT(IV)*SIGT(JV)*VOLSUR(IV)*PVV(IV,JV)
C           SIGT(IV).GT.0.0 AND SIGT(JV).EQ.0.0
C             PROB(I)=SIGT(IV)*VOLSUR(IV)*PVV(IV,JV)
C           SIGT(IV).EQ.0.0 AND SIGT(JV).GT.0.0
C             PROB(I)=SIGT(JV)*VOLSUR(IV)*PVV(IV,JV)
C           SIGT(IV).EQ.0.0 AND SIGT(JV).EQ.0.0
C             PROB(I)=VOLSUR(IV)*PVV(IV,JV)
C
C-------------------------    PIJWPR    -------------------------------
C
      IMPLICIT            NONE
      INTEGER           IUNOUT,   LOPT,   NREG,  NSOUT,   MSYM,
     >                     IND,      I,      J,   NSUR,   NVOL,
     >                   NSURC,  NSURM,  NVOLC,  NVOLM,     IP,     IR,
     >                      JR,    III
      PARAMETER        (IUNOUT=6)
      REAL              SIGTAL(-NSOUT:NREG), BILANP(10), XSJR, WPR,
     >                  VPR(10), SIGVOL(NREG), COF
      DOUBLE PRECISION  PROB(*)
C
      IND(I,J) = MAX(I+NSOUT+1,J+NSOUT+1)*(MAX(I+NSOUT+1,J+NSOUT+1)-1)/2
     >         + MIN(I+NSOUT+1,J+NSOUT+1)
C
      WPR(I,J)= REAL(PROB( IND(I,J) ) / PROB( IND(I,0) ))
C
CNOTE:
C     IF( SIGT(I).NE.0.0 )THEN
C        PROB(IND(I,0)= SIGT(I) * VOLSUR(I)
C     ELSE
C        PROB(IND(I,0)= VOLSUR(I)
C     ENDIF
C
      NSUR= -NSOUT
      NVOL=  NREG
C
      WRITE(IUNOUT,'(24H REGIONAL CROSS SECTIONS)')
      WRITE(IUNOUT,'(5(1X,6HREGION,5X,16HCROSS SECTION   ))')
      WRITE(IUNOUT,'(5(1X,I6,3X,E15.7))') (JR,SIGTAL(JR),JR=1,NVOL)
      IF(MSYM .EQ. 0) THEN
        WRITE(IUNOUT,'(5(1X,6HREGION,5X,16HVOLUMES         ))')
        WRITE(IUNOUT,'(5(1X,I6,3X,E15.7))') 
     >  (JR,SIGVOL(JR),JR=1,NVOL)
      ELSE
        WRITE(IUNOUT,'(5(1X,6HREGION,5X,16HSURFACE/VOLUMES ))')
        WRITE(IUNOUT,'(5(1X,I6,3X,E15.7))') 
     >  (JR,PROB(IND(JR,0)),JR=-NSOUT,NVOL)
      ENDIF
      IF( LOPT.LE.0 )THEN
         NSURC = -1
         DO 40 IP  = 1, (9 - NSUR) / 10
            NSURM= MAX( NSUR, NSURC-9 )
            WRITE(IUNOUT,'(30H0  SURFACE CONSERVATION LAWS: ,
     >                     31H( P.S<-S + P.V<-S = 1 + ERR.S ) ,
     >                     31H FOR XS.TOTAL=0, REDUCED P.V<-S ,
     >                     11H IS PRINTED )')
            WRITE(IUNOUT,'(1X,8H(P.S<-S),1X,10( A5,    I6,:)/)')
     >               (' SUR ',-IR,IR= NSURC, NSURM, -1)
            DO 10 IR  =NSURC, NSURM, -1
               BILANP(IR-NSURM+1)= 0.0
   10       CONTINUE
            DO 20 JR = -1,  NSUR, -1
               WRITE(IUNOUT,'(5H SUR ,I4,1X,10F11.8)')
     >                          -JR, (WPR(IR,JR),IR=NSURC,NSURM,-1)
               DO 20 IR  = NSURC, NSURM, -1
                  BILANP(IR-NSURM+1)=  BILANP(IR-NSURM+1)
     >                              +  WPR(IR,JR)
   20       CONTINUE
            WRITE(IUNOUT,'(1X,8H(P.V<-S) )')
            DO 30 JR  =  1,  NVOL,  1
               IF( SIGTAL(JR).EQ.0.0 ) THEN
                  WRITE(IUNOUT,'(5H VOL ,I4,1X,3H 0*,10(F8.5,:,3H 0*))')
     >                      JR,(WPR(IR,JR),IR=NSURC,NSURM,-1)
                  XSJR= 0.0
               ELSE
                  WRITE(IUNOUT,'(5H VOL ,I4,1X,10F11.8)')
     >                      JR,(WPR(IR,JR),IR=NSURC,NSURM,-1)
                  XSJR= 1.0
               ENDIF
               DO 30 IR  = NSURC, NSURM, -1
                  BILANP(IR-NSURM+1)= BILANP(IR-NSURM+1)
     >                              + XSJR * WPR(IR,JR)
   30       CONTINUE
            WRITE(IUNOUT,'(1H )')
            WRITE(IUNOUT,'(5H SUM ,5X,10F11.8)')
     >                          (BILANP(IR-NSURM+1),IR=NSURC,NSURM,-1)
            NSURC = NSURC - 10
   40    CONTINUE
      ENDIF
      NVOLC =  1
      DO 90 IP  = 1, (9 + NVOL) / 10
         NVOLM= MIN( NVOL, NVOLC+9 )
         IF( LOPT.LE.0 )THEN
            WRITE(IUNOUT,'(30H0  VOLUME  CONSERVATION LAWS: ,
     >                     31H( P.S<-V + P.V<-V = 1 + ERR.V ) ,
     >                     31H FOR XS.TOTAL=0, REDUCED P.V<-V ,
     >                     11H IS PRINTED )')
         ELSE
            WRITE(IUNOUT,'(30H0  VOLUME  CONSERVATION LAWS: ,
     >                     32H( SUM OF P.V<-V = 1 + ESCAPE.V ) ,
     >                     31H FOR XS.TOTAL=0, REDUCED P.V<-V ,
     >                     11H IS PRINTED )')
         ENDIF
         DO 50 IR  = NVOLC, NVOLM,   1
            BILANP(IR-NVOLC+1)= 0.0
   50    CONTINUE
         IF( LOPT.LE.0 )THEN
            WRITE(IUNOUT,'(1X,8H(P.S<-V),1X,10( A5 ,  I6,:)/)')
     >                    (' VOL ',IR,IR=NVOLC,NVOLM, 1)
            DO 60 JR = -1,  NSUR, -1
               WRITE(IUNOUT,'(5H SUR ,I4,1X,10F11.8)')
     >         -JR, (WPR(IR,JR),IR=NVOLC,NVOLM, 1)
            DO 60 IR  = NVOLC, NVOLM,   1
               BILANP(IR-NVOLC+1)= BILANP(IR-NVOLC+1)
     >                           + WPR(IR,JR)
   60       CONTINUE
            WRITE(IUNOUT,'(1X,8H(P.V<-V) )')
         ELSE
            WRITE(IUNOUT,'(1X,8H(P.V<-V),1X,10( A5 ,  I6,:)/)')
     >                    (' VOL ',IR,IR=NVOLC,NVOLM, 1)
         ENDIF
C
         IF(LOPT.GT.0.AND.MSYM.EQ.0)THEN
C
C        PRINTING OF PIJK" FULL MATRIX
C
         COF=1.5
         DO 70 JR =  1,  NVOL,  1
            IF( SIGTAL(JR).EQ.0.0 )THEN
              DO 75 IR=NVOLC, NVOLM, 1
               III=JR+NREG*(IR-1)
               VPR(IR-NVOLC+1)=COF*REAL(PROB(III))/SIGVOL(IR)
   75         CONTINUE
               WRITE(IUNOUT,'(5H VOL ,I4,1X,3H 0*,10(F8.5,:,3H 0*))')
     >         JR,(VPR(IR-NVOLC+1),IR=NVOLC,NVOLM,1)
            ELSE
              DO 76 IR=NVOLC, NVOLM, 1
               III=JR+NREG*(IR-1)
               VPR(IR-NVOLC+1)=COF*REAL(PROB(III))/SIGVOL(IR)
               BILANP(IR-NVOLC+1)=BILANP(IR-NVOLC+1)+VPR(IR-NVOLC+1)
   76         CONTINUE
               WRITE(IUNOUT,'(5H VOL ,I4,1X,10F11.8)')
     >         JR,(VPR(IR-NVOLC+1),IR=NVOLC,NVOLM,1)
            ENDIF
   70    CONTINUE
         WRITE(IUNOUT,'(1H )')
         WRITE(IUNOUT,'(5H SUM ,5X,10F11.8)')
     >                         (BILANP(IR-NVOLC+1),IR=NVOLC,NVOLM, 1)
C
         ELSE
C
         DO 80 JR =  1,  NVOL,  1
            IF( SIGTAL(JR).EQ.0.0 )THEN
               WRITE(IUNOUT,'(5H VOL ,I4,1X,3H 0*,10(F8.5,:,3H 0*))')
     >         JR,(WPR(IR,JR),IR=NVOLC,NVOLM,1)
               XSJR= 0.0
            ELSE
               WRITE(IUNOUT,'(5H VOL ,I4,1X,10F11.8)')
     >         JR,(WPR(IR,JR),IR=NVOLC,NVOLM,1)
               XSJR= 1.0
            ENDIF
         DO 80 IR  = NVOLC, NVOLM,   1
            BILANP(IR-NVOLC+1)= BILANP(IR-NVOLC+1)
     >                        + XSJR*WPR(IR,JR)
   80    CONTINUE
         WRITE(IUNOUT,'(1H )')
         WRITE(IUNOUT,'(5H SUM ,5X,10F11.8)')
     >                         (BILANP(IR-NVOLC+1),IR=NVOLC,NVOLM, 1)
         ENDIF
         NVOLC = NVOLC + 10
   90 CONTINUE
C
      RETURN
      END
