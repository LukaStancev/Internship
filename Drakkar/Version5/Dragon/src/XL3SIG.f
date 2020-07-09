*DECK XL3SIG
      SUBROUTINE XL3SIG(  NGRT,  NBMIX,  XSSIGT,ALBEDO,  NPSYS,
     >                    NGRP,     NS,     NR, MATALB,    VOL,
     >                  SIGTAL, SIGVOL, SWVOID, SWNZBC)
************************************************************************
*                                                                      *
*           NAME: XL3SIG                                               *
*      COMPONENT: EXCELL                                               *
*          LEVEL: 4 (CALLED BY 'XL3TRK')                               *
*        VERSION: 1.0                                                  *
*       CREATION: 96/09 (R.R.)                                         *
*       MODIFIED: 00/03 (R.R.) DECLARE ALL VARIABLE TYPES              *
*         AUTHOR: ROBERT ROY                                           *
*                                                                      *
*     SUBROUTINE: THIS ROUTINE IS USED TO UNFOLD CROSS-SECTION DATA    *
*                 WHICH BECOMES AVAILABLE BY SUBSET OF GROUPS          *
*                                                                      *
*--------+-------------- V A R I A B L E S -------------+--+-----------*
*  NAME  /                  DESCRIPTION                 /IO/MOD(DIMENS)*
*--------+----------------------------------------------+--+-----------*
* NGRT   / TOTAL NUMBER OF GROUPS                       /I./INT        *
* NBMIX  / # OF MIXTURES IN THE MACROLIB                /I./INT        *
* XSSIGT / TOTAL XS FOR MIXTURES IN THE MACROLIB        /I./REL(NBMIX  *
*        /                                              /I./  +1,NGRT) *
* ALBEDO / GEOMETRIC ALBEDOS.                           /I./REL(6)     *
* NPSYS  / GROUP MASKS                                  /I./INT(NGRP)  *
* NGRP   / NUMBER OF GROUPS                             /I./INT        *
* NS     / # OF SURFACES IN THE ASSEMBLY.               /I./INT        *
* NR     / # OF ZONES IN THE ASSEMBLY.                  /I./INT        *
* MATALB / MATERIAL #S FOR ZONES IN THE SUPERCELL       /I./INT(NS:NR) *
* VOL    / VOLUMES                                      /I./INT        *
* SIGTAL / TOTAL XS & ALBEDOS BY REGION & SURFACE       /.O/REL(NS:NR, *
*        /                                              /  /     NGRP) *
* SIGVOL / VOLUME TIMES TOTAL XS BY REGION              /.O/REL(NR,    *
*        /                                              /  /     NGRP) *
* SWVOID / LOGICAL SWITCH (.TRUE. IF VOID REGIONS)      /.O/LOG        *
* SWNZBC / LOGICAL SWITCH (.TRUE. IF NON-ZERO B.C.)     /.O/LOG        *
************************************************************************
      IMPLICIT   NONE
C
      INTEGER    NGRT,NBMIX,NGRP,NS,NR,NPSYS(NGRP)
      REAL       XSSIGT(0:NBMIX,NGRT),VOL(NR),SIGTAL(NS:NR,NGRP),
     >           SIGVOL(NR,NGRP),ALBEDO(6)
      INTEGER    MATALB(NS:NR)
      INTEGER    IUN,JG,I
      LOGICAL    SWVOID,SWNZBC
      REAL       ZERO
      INTEGER    IOUT
      PARAMETER (ZERO=0.0,IOUT=6)
C
      SWVOID= .FALSE.
      SWNZBC= .FALSE.
C
      IF( NS.GT.0 ) CALL XABORT('XL3SIG: # OF SURFACES IS > 0')
      IF( NR.LT.0 ) CALL XABORT('XL3SIG: # OF REGIONS  IS < 0')
C
      DO 10 IUN= NS, NR
      DO 20 JG= 1, NGRP
         IF(NPSYS(JG).EQ.0) GO TO 20
         IF( IUN.LT.0 )THEN
            SIGTAL(IUN,JG)= ALBEDO(-MATALB(IUN))
            SWNZBC=SWNZBC.OR.(ALBEDO(-MATALB(IUN)).NE.ZERO)
         ELSEIF( IUN.EQ.0 )THEN
            SIGTAL(IUN,JG)= ZERO
         ELSE
            IF( MATALB(IUN).LT.0.OR.MATALB(IUN).GT.NBMIX)THEN
               WRITE(IOUT,*) 'NBMIX=',NBMIX
               WRITE(IOUT,*) 'IG/NGRT=',JG,NGRT
               WRITE(IOUT,*) 'XSSIGT=',(XSSIGT(I,JG),I=0,NBMIX)
               WRITE(IOUT,*) 'MATALB<=0 =',(MATALB(I),I=NS,0)
               WRITE(IOUT,*) 'MATALB >0 =',(MATALB(I),I=1,NR)
               CALL XABORT('XL3SIG: INVALID NUMBER OF MIXTURES')
            ENDIF
            SIGTAL(IUN,JG)= XSSIGT(MATALB(IUN),JG)
            SIGVOL(IUN,JG)= SIGTAL(IUN,JG)*VOL(IUN)
         ENDIF
   20 CONTINUE
   10 CONTINUE
C
      RETURN
      END
