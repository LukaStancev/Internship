*DECK XL3NTR
      SUBROUTINE XL3NTR(  IPRT,   NDIM,   ISPEC, NS, NV, NORE,
     >                   VOLIN,  MRGIN,   MATIN,
     >                   NANGL, VOLTRK,  DENSTY )
************************************************************************
*                                                                      *
*           NAME: XL3NTR                                               *
*      COMPONENT: EXCELL                                               *
*          LEVEL: 2 (CALLED BY 'EXCELT')                               *
*        VERSION: 1.0                                                  *
*       CREATION: 91/07                                                *
*       MODIFIED: 00/03 (R.R.) DECLARE ALL VARIABLE TYPES              *
*         AUTHOR: ROBERT ROY                                           *
*                                                                      *
*     SUBROUTINE: THIS ROUTINE WILL COMPUTE RENORMALIZED TRACKS        *
*                 TO OBTAIN TRUE VOLUME VALUES. THE FILE "IFOLD"       *
*                 CONTAINS THE OLD TRACKS; THE FILE "IFTRAK" WILL      *
*                 CONTAIN THE NORMALIZED TRACKS.                       *
*                                                                      *
*           NOTE: THE FILES "IFOLD" AND "IFTRAK" ARE SUPPOSED TO BE:   *
*                 1) CONNECTED AND OPENED;                             *
*                 2) PLACED FOR ACCESSING THE FIRST RECORD (REWIND).   *
*                                                                      *
*--------+-------------- V A R I A B L E S -------------+--+-----------*
*  NAME  /                  DESCRIPTION                 /IO/MOD(DIMENS)*
*--------+----------------------------------------------+--+-----------*
* IPRT   / INTERMEDIATE PRINTING LEVEL FOR PRINOUT.     /I./INT        *
* NDIM   / # OF DIMENSIONS (2D OR 3D).                  /I./INT        *
* ISPEC  / KIND OF TRACKING (0:ISOTROPIC;1:SPECULAR)    /I./INT        *
* NS     / # OF SURFACES BEFORE MERGING.                /I./INT        *
* NV     / # OF ZONES BEFORE MERGING.                   /I./INT        *
* NORE   / TRACK NORMALIZATION (-1:YES 1:NO)            /I./I          *
* VOLIN  / VOLUMES & SURFACES BEFORE MERGING.           /I./REL(-NS:NV)*
* MRGIN  / MERGING INDEX.                               /I./INT(-NS:NV)*
* MATIN  / MATERIAL #S BEFORE MERGING.                  /I./INT(-NS:NV)*
* NANGL  / # OF ANGLES TO RENORMALIZE TRACKS BY ANGLE.  /I./INT        *
* VOLTRK / VOLUMES & SURFACES AS COMPUTED BY TRACKING.  /../R*8(-NS:NV,*
*        /                                              /../   0:NANGL)*
* DENSTY / WEIGHTS BY ANGLE.                            /../REL(NANGL) *
************************************************************************
      IMPLICIT           NONE
C
      INTEGER            NDIM,NS,NV,NANGL,IPRT,IANG,IP,IR,ISPEC,ITGEO,
     >                   IVS,IVSC,MNSUR,MXVOL,NANG2,IOUT,NORE,
     >                   NSURC,NSURM,NVOLC,NVOLM,MRGIN(-NS:NV),
     >                   MATIN(-NS:NV),NTMP,JR
      REAL               VOLIN(-NS:NV),
     >                   DENSTY(NANGL),
     >                   ERRSUR,ERRVOL,ERRVM,ERRSM,TMPERR(10)
      DOUBLE PRECISION   VOLTRK(-NS:NV,0:NANGL),APRSUR,APRVOL,
     >                   TOTVOL,TOTSUR,ZERO,ONE,TWO,FOUR,HALF,QUART,
     >                   HUND,PI,FACVOL,FACSUR
      CHARACTER          CORIEN(0:3,-6:-1)*4
      PARAMETER        ( PI=3.14159265358979323846D0, IOUT=6,
     >                   ZERO=0.D0, ONE=1.D0, TWO=2.D0, FOUR=4.D0,
     >                   HUND=1.D2, HALF=0.5D0, QUART=0.25D0, ITGEO=3 )
      DATA         ((CORIEN(JR,IR),IR=-6,-1),JR=0,3)
     >             / ' 6  ',' 5  ',' 4  ',' 3  ',' 2  ',' 1  ',
     >               ' Z+ ',' Z- ','****','****',' R+ ','****',
     >               ' Z+ ',' Z- ','****','****','****','HBC ',
     >               ' Z+ ',' Z- ',' Y+ ',' Y- ',' X+ ',' X- ' /
C
      FACVOL= TWO
      FACSUR= ONE
      IF( ISPEC.EQ.0 )THEN
         IF( NDIM.EQ.2 )THEN
            FACSUR= QUART*PI
         ELSEIF( NDIM.EQ.3 )THEN
            FACSUR= ONE
         ENDIF
      ELSEIF( ISPEC.EQ.1 )THEN
         IF( NDIM.EQ.2 )THEN
            FACSUR= HALF*PI
         ELSEIF( NDIM.EQ.3 )THEN
            FACSUR= ONE
         ENDIF
      ENDIF
      DO 47 IVS=  -NS, NV
        DO 46 IANG= 1, NANGL
           VOLTRK(IVS,0)= VOLTRK(IVS,0) + VOLTRK(IVS,IANG)
           VOLTRK(IVS,IANG)= VOLTRK(IVS,IANG)*DENSTY(IANG)
           IF( VOLTRK(IVS,IANG).NE.ZERO )THEN
C
C             CONVERT INTO NORMALIZATION FACTORS
              VOLTRK(IVS,IANG)= VOLIN(IVS)/VOLTRK(IVS,IANG)
           ELSE
              VOLTRK(IVS,IANG)= ONE
           ENDIF
   46   CONTINUE
   47 CONTINUE
C
C     COMPUTE ERRORS FOR CONSERVATION LAWS
      TOTSUR=ZERO
      APRSUR=ZERO
      TOTVOL=ZERO
      APRVOL=ZERO
      ERRSM=0.0
      ERRVM=0.0
      IVSC=0
      DO 50 IVS= -NS, NV
        IF( VOLTRK(IVS,0).EQ.ZERO.AND.VOLIN(IVS).GT.0.0)THEN
           IVSC= IVS
        ENDIF
        IF( IVS.LT.0 )THEN
           VOLTRK(IVS,0)= REAL(FACSUR)*VOLTRK(IVS,0)
           IF(VOLIN(IVS).NE.0.0) THEN
             ERRSM=MAX(ERRSM,
     >         REAL(100.0*ABS(1.0-VOLTRK(IVS,0)/VOLIN(IVS))))
           ENDIF
           TOTSUR=TOTSUR+VOLIN(IVS)
           APRSUR=APRSUR+VOLTRK(IVS,0)
        ELSEIF( IVS.GT.0 )THEN
           VOLTRK(IVS,0)= FACVOL*VOLTRK(IVS,0)
           TOTVOL=TOTVOL+VOLIN(IVS)
           APRVOL=APRVOL+VOLTRK(IVS,0)
           IF(VOLIN(IVS).NE.0.0) THEN
             ERRVM=MAX(ERRVM,
     >         REAL(100.0*ABS(1.0-VOLTRK(IVS,0)/VOLIN(IVS))))
           ENDIF
        ENDIF
   50 CONTINUE
      ERRSUR=100.*REAL(1.0-APRSUR/TOTSUR)
      ERRVOL=100.*REAL(1.0-APRVOL/TOTVOL)
      IF( IPRT.GT.1 )THEN
         MNSUR = -NS
         MXVOL =  NV
         NSURC = -1
         WRITE(IOUT,'(1H )')
         WRITE(IOUT,7000) ERRSUR,ERRSM
         DO 80 IP  = 1, (9 - MNSUR) / 10
            NSURM= MAX( MNSUR, NSURC-9 )
            WRITE(IOUT,'(10X,10(A5,I6))')(' FACE',-IR,IR=NSURC,NSURM,-1)
            WRITE(IOUT,'(8H SURFACE,2X,1P,10E11.4)')
     >                              (4.*VOLIN(IR),IR=NSURC,NSURM,-1)
            WRITE(IOUT,'(8H SIDE   ,2X,10(A4,7X))')
     >               (CORIEN(ITGEO,MATIN(IR)),IR=NSURC,NSURM,-1)
            WRITE(IOUT,'(8H APPROX ,2X,1P,10E11.4)')
     >                         (FOUR*VOLTRK(IR,0),IR=NSURC,NSURM,-1)
            NTMP=0
            DO 81  IR=NSURC,NSURM,-1
              NTMP=NTMP+1
              IF(VOLIN(IR).NE.0.0) THEN
                TMPERR(NTMP)=REAL(HUND-HUND*VOLTRK(IR,0)/VOLIN(IR))
              ELSE
                TMPERR(NTMP)=0.0
              ENDIF
  81        CONTINUE
            WRITE(IOUT,'(8H ERR(%) ,2X,10F11.5)')
     >              (TMPERR(IR),IR=1,NTMP)
            WRITE(IOUT,'(9H MERGE TO,1X,10(A5,I6))')
     >               (' FACE',-MRGIN(IR),IR=NSURC,NSURM,-1)
            WRITE(IOUT,'(1H )')
            NSURC = NSURC - 10
   80    CONTINUE
         NVOLC= 1
         WRITE(IOUT,'(1H )')
         WRITE(IOUT,7001) ERRVOL,ERRVM
         DO 90 IP  = 1, (9 + MXVOL) / 10
            NVOLM= MIN( MXVOL, NVOLC+9 )
            WRITE(IOUT,'(10X,10(A5,I6))') (' ZONE',IR,IR=NVOLC,NVOLM)
            WRITE(IOUT,'(8H VOLUME ,2X,1P,10E11.4)')
     >                                  (VOLIN(IR),IR=NVOLC,NVOLM)
            WRITE(IOUT,'(9H MIXTURE ,1X,10(A5,I6))')
     >                         (' MIX ', MATIN(IR),IR=NVOLC,NVOLM)
            WRITE(IOUT,'(8H APPROX ,2X,1P,10E11.4)')
     >                               (VOLTRK(IR,0),IR=NVOLC,NVOLM)
            NTMP=0
            DO 91  IR= NVOLC,NVOLM
              NTMP=NTMP+1
              IF(VOLIN(IR).NE.0.0) THEN
                TMPERR(NTMP)=REAL(HUND-HUND*VOLTRK(IR,0)/VOLIN(IR))
              ELSE
                TMPERR(NTMP)=0.0
              ENDIF
  91        CONTINUE
            WRITE(IOUT,'(8H ERR(%) ,2X,10F11.5)')
     >              (TMPERR(IR),IR=1,NTMP)
            WRITE(IOUT,'(9H MERGE TO,1X,10(A5,I6))')
     >                        (' ZONE',MRGIN(IR),IR=NVOLC,NVOLM)
            WRITE(IOUT,'(1H )')
            NVOLC = NVOLC + 10
   90    CONTINUE
         IF( IPRT.GT.5 )THEN
            NVOLC= 1
            NANG2= NANGL+2
            WRITE(IOUT,'(1H )')
            IF( NORE.EQ.-1 )THEN
               WRITE(IOUT,7002)
            ELSE IF( NORE.EQ.1 )THEN
               WRITE(IOUT,7003)
            ELSE
               CALL XABORT('XL3NTR: INVALID NORMALIZATION OPTION.')
            ENDIF
            DO 110 IP  = 1, (9 + MXVOL) / 10
               NVOLM= MIN( MXVOL, NVOLC+9 )
               WRITE(IOUT,'(10X,10(A5,I6))') (' VOL ',IR,IR=NVOLC,NVOLM)
               DO 100 IANG= 1, NANGL
                  WRITE(IOUT,'(4H ANG,I4 ,2X,1P,10E11.4)')
     >            IANG, (VOLTRK(IR,IANG),IR=NVOLC,NVOLM)
  100          CONTINUE
               WRITE(IOUT,'(1H )')
               NVOLC = NVOLC + 10
  110       CONTINUE
         ENDIF
      ENDIF
      IF( IVSC.NE.0 )THEN
         WRITE(IOUT,*) ' VOLUME # ',IVSC,' NOT TRACKED'
         WRITE(IOUT,*) ' USE FINER TRACKING'
         CALL XABORT( 'XL3NTR: CHECK NUMBERING OR USE FINER TRACKING')
      ENDIF
C
      RETURN
 7000 FORMAT(/' TRACKING ERRORS ON SURFACE   AVERAGE ERROR: ',F10.4,
     >        ' % ',5X,'MAXIMUM ERROR: ',F10.4,' % (BEFORE MERGE)')
 7001 FORMAT( ' TRACKING ERRORS ON VOLUME    AVERAGE ERROR: ',F10.4,
     >        ' % ',5X,'MAXIMUM ERROR: ',F10.4,' % (BEFORE MERGE)')
 7002 FORMAT(/' ANGLE-BY-ANGLE RENORMALIZATION FACTORS: '/)
 7003 FORMAT(/' ANGLE-BY-ANGLE RENORMALIZATION FACTORS(**NOT USED): '/)
      END
