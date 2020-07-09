*DECK SNDSA
      SUBROUTINE SNDSA (KPSYS,INCONV,INGIND,IPTRK,IMPX,NGRP,NGEFF,NREG,
     1 NBMIX,NUN,ISCAT,MAT,VOL,KEYFLX,NUNSA,ZCODE,FUNOLD,FUNKNO,NHEX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* perform a synthetic acceleration using BIVAC (2D) or TRIVAC (3D)    
* for the discrete ordinates (SN) method using P1 approximation.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert and N. Martin and A. A. Calloo
*
*Parameters: input
* KPSYS   pointer to the assembly matrices. KPSYS is an array of
*         directories.
* INCONV  energy group convergence flag (set to .false. if converged).
* INGIND  energy group index assign to 1:NGEFF arrays.
* IPTRK   pointer to the tracking (L_TRACK signature).
* IMPX    print flag (equal to zero for no print).
* NGRP    number of energy groups.
* NGEFF   dimension of arrays KPSYS, INCONV and INGIND.
* NREG    total number of regions for which specific values of the
*         neutron flux and reactions rates are required.
* NBMIX   number of mixtures.
* NUN     total number of unknowns in vectors SUNKNO and FUNKNO.
* ISCAT   anisotropy of one-speed sources in SN method.
* MAT     index-number of the mixture type assigned to each volume.
* VOL     volumes.
* KEYFLX  position of averaged flux elements in FUNKNO vector.
* SUNKNO  input source vector.
* NUNSA   number of unknowns in BIVAC/TRIVAC.
* ZCODE   albedos.
* FUNOLD  SN unknown vector at iteration kappa.
* FUNKNO  SN unknown vector at iteration kappa+1/2.
*
*Parameters: output
* FUNKNO  SN unknown vector at iteration kappa+1 (with DSA correction).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) KPSYS(NGEFF),IPTRK
      INTEGER     NGEFF,INGIND(NGEFF),IMPX,NGRP,NREG,NBMIX,NUN,ISCAT,
     >            MAT(NREG),KEYFLX(NREG),NUNSA,NHEX
      LOGICAL     INCONV(NGEFF)
      REAL        VOL(NREG),ZCODE(6),FUNOLD(NUN,NGRP),FUNKNO(NUN,NGRP)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (IUNOUT=6,NSTATE=40,PI=3.141592654)
      TYPE(C_PTR) JPSYS 
      INTEGER IPAR(NSTATE),NLOZH,SPLTL,REM,SBMSH
*
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IDL,TMPKEY
      REAL, ALLOCATABLE, DIMENSION(:) :: SGAS
      REAL, ALLOCATABLE, DIMENSION(:,:) :: FUNSA,SUNSA
*
      TYPE(C_PTR) DU_PTR,DE_PTR,W_PTR,DZ_PTR
      REAL, POINTER, DIMENSION(:) :: DU,DE,W,DZ
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IDL(NREG),FUNSA(NUNSA,NGRP),SUNSA(NUNSA,NGRP))
*----
*  RECOVER TRACKING INFORMATION
*----
      CALL LCMSIX(IPTRK,'DSA',1)
      CALL LCMGET(IPTRK,'STATE-VECTOR',IPAR)
      IF(NREG.NE.IPAR(1)) CALL XABORT('SNDSA: INVALID NREG ON LCM.')
      IF(NUNSA.NE.IPAR(2)) CALL XABORT('SNDSA: INVALID NUN ON LCM.')
      ITYPE=IPAR(6)
      ISPLH=1
      IF(ITYPE.EQ.7) THEN
        IELEM=IPAR(9)
        LX=IPAR(14)
        LY=IPAR(15)
        LZ=IPAR(16)
      ELSE
        IELEM=IPAR(8)
        LX=IPAR(12)
        LY=IPAR(13)
        LZ=0
      ENDIF
      NSCT=0
      NLEG=0
      IF(ITYPE.EQ.2) THEN
         NSCT=ISCAT
         NLEG=IELEM
      ELSE IF((ITYPE.EQ.5).OR.(ITYPE.EQ.8)) THEN
         NSCT=ISCAT*(ISCAT+1)/2
         NLEG=IELEM*IELEM
      ELSE IF((ITYPE.EQ.7).OR.(ITYPE.EQ.9)) THEN
         NSCT=(ISCAT)**2
         NLEG=IELEM*IELEM*IELEM
      ELSE
         CALL XABORT('SNDSA: TYPE OF DISCRETIZATION NOT IMPLEMENTED.')
      ENDIF
      IF(ITYPE.EQ.8) ISPLH=IPAR(10)
      CALL LCMGET(IPTRK,'KEYFLX',IDL)
*----
*  LOOP OVER ENERGY GROUPS.
*----
      CALL XDRSET(SUNSA,NUNSA*NGRP,0.0)
      DO 30 II=1,NGEFF
      IF(.NOT.INCONV(II)) GO TO 30
      JPSYS=KPSYS(II)
      IG=INGIND(II)
      IF(IMPX.GT.1) WRITE(IUNOUT,'(/24H SNDSA: PROCESSING GROUP,I5,
     1 6H WITH ,A,1H.)') IG,'SN/DSA'
*----
*  RECOVER WITHIN-GROUP SCATTERING CROSS SECTION.
*----
      CALL LCMLEN(JPSYS,'DRAGON-TXSC',ILONG,ITYLCM)
      IF(ILONG.NE.NBMIX+1) CALL XABORT('SNDSA: INVALID TXSC LENGTH.')
      CALL LCMLEN(JPSYS,'DRAGON-S0XSC',ILONG,ITYLCM)
      NANI=ILONG/(NBMIX+1)
      ALLOCATE(SGAS(ILONG))
      CALL LCMGET(JPSYS,'DRAGON-S0XSC',SGAS)
*----
*  REBUILD KEYFLX FOR HEXAGONAL CASE
*----
      ! NLOZH - num. of loz. per hexagon
      ! SBMSH - num. of submeshes per lozenge (integer)
      ! SPLTL - split of the lozenge (ISPLH)
      IF(ITYPE.EQ.8)THEN
         ALLOCATE(TMPKEY(NREG))
         IND = 0
         JND = 0
         NLOZH  = 3*ISPLH**2
         SBMSH  = NLOZH/3
         SPLTL  = ISPLH
         DO IH=1,NHEX
            DO IM=1,SBMSH
               REM=MOD(IM-1,SPLTL)
               IF((REM.EQ.0).AND.(SBMSH.NE.1))THEN
                  JND = (IH-1)*NLOZH + SBMSH - (IM/SPLTL)
               ELSEIF((REM.NE.0).AND.(SBMSH.NE.1))THEN
                  JND = JND - (SBMSH*3) - SPLTL
               ENDIF
               DO ILZ=1,3
                  IND = (IH-1)*NLOZH + (IM-1)*3 + (ILZ-1) + 1
                  IF(SBMSH.EQ.1) JND = IND
                  TMPKEY(IND) = KEYFLX(JND)
                  JND = JND + SBMSH
               ENDDO
            ENDDO
         ENDDO
         KEYFLX(:) = TMPKEY(:)
         DEALLOCATE(TMPKEY)
      ELSEIF(ITYPE.EQ.9)THEN
         CALL XABORT('SNDSA: 3D HEXAGONAL NOT IMPLEMENTED YET.')
      ENDIF
*----
*  COMPUTE THE SOURCE OF THE DSA EQUATION.
*----
      DO 20 IR=1,NREG
      IBM=MAT(IR)
      IF(IBM.LE.0) GO TO 20
      SIGS=SGAS(IBM+1)
      DO 10 IEL=1,NLEG
      IND=KEYFLX(IR)+IEL-1
      JND=IDL(IR)+IEL-1
      SUNSA(JND,IG)=SUNSA(JND,IG)+SIGS*(FUNKNO(IND,IG)-FUNOLD(IND,IG))
  10  CONTINUE       
  20  CONTINUE      
      DEALLOCATE(SGAS)
  30  CONTINUE
*----
*  SOLVE THE DSA EQUATION USING A P1 METHOD.
*----
      CALL XDRSET(FUNSA,NUNSA*NGRP,0.0)
      IF(ITYPE.EQ.7) THEN  
        CALL TRIFLV(KPSYS,INCONV,INGIND,IPTRK,IMPX,NGRP,NGEFF,NREG,
     >  NUNSA,MAT,VOL,IDL,FUNSA,SUNSA)
      ELSE
        CALL PNFLV(KPSYS,INCONV,INGIND,IPTRK,IMPX,NGRP,NGEFF,NREG,
     >  NBMIX,NUNSA,MAT,VOL,IDL,FUNSA,SUNSA)
      ENDIF
      CALL LCMSIX(IPTRK,' ',2)
*----
*  LOOP OVER ENERGY GROUPS.
*----
      DO 400 II=1,NGEFF
      IF(.NOT.INCONV(II)) GO TO 400
      IG=INGIND(II)
*----
*  UPGRADE THE SN SURFACE FLUX IN 2D CARTESIAN CASES
*----
      IF(ITYPE.EQ.5) THEN  
         CALL LCMLEN(IPTRK,'DU',NPQ,ITYLCM)
         CALL LCMGPD(IPTRK,'DU',DU_PTR)
         CALL LCMGPD(IPTRK,'DE',DE_PTR)
         CALL LCMGPD(IPTRK,'W',W_PTR)
         CALL C_F_POINTER(DU_PTR,DU,(/ NPQ /))
         CALL C_F_POINTER(DE_PTR,DE,(/ NPQ /))
         CALL C_F_POINTER(W_PTR,W,(/ NPQ /))
         IR=0
         DO 161 K1=1,LY
         DO 160 K2=1,LX
         IR=IR+1
         IF(MAT(IR).EQ.0) GO TO 160
         IF(VOL(IR).EQ.0.0) GO TO 150
*******XNEI-
         IF((K2.EQ.1).AND.(ZCODE(1).NE.0.0)) THEN
            DO 55 J0=1,IELEM
            FSN=0.0
            DFPN=0.0
            SG=1.0
            DO 40 I0=1,IELEM
            IND0=KEYFLX(IR)+(J0-1)*IELEM+I0-1
            IND1=IDL(IR)+(J0-1)*IELEM+I0-1
            FSN=FSN+SG*SQRT(REAL(2*I0-1))*FUNKNO(IND0,IG)
            DFPN=DFPN+SG*SQRT(REAL(2*I0-1))*FUNSA(IND1,IG)
            SG=-SG
   40       CONTINUE
            BSN=0.0
            DO 45 M=1,NPQ
            IF(DU(M).LT.0.0) THEN
               IND0=LY*LX*NSCT*IELEM**2+(M-1)*LY*IELEM+(K1-1)*IELEM+J0
               BSN=BSN+2.0*W(M)*FUNKNO(IND0,IG)*(1.0+ZCODE(1))
            ENDIF
   45       CONTINUE
            IF(FSN.NE.0.0) DFPN=DFPN*BSN/FSN
            DO 50 M=1,NPQ
            IF(DU(M).LT.0.0) THEN
               IND0=LY*LX*NSCT*IELEM**2+(M-1)*LY*IELEM+(K1-1)*IELEM+J0
               FUNKNO(IND0,IG)=FUNKNO(IND0,IG)+DFPN/(4.0*PI)
            ENDIF
   50       CONTINUE
   55       CONTINUE
         ENDIF
*******XNEI+
         IF((K2.EQ.LX).AND.(ZCODE(2).NE.0.0)) THEN
            DO 75 J0=1,IELEM
            FSN=0.0
            DFPN=0.0
            DO 60 I0=1,IELEM
            IND0=KEYFLX(IR)+(J0-1)*IELEM+I0-1
            IND1=IDL(IR)+(J0-1)*IELEM+I0-1
            FSN=FSN+SQRT(REAL(2*I0-1))*FUNKNO(IND0,IG)
            DFPN=DFPN+SQRT(REAL(2*I0-1))*FUNSA(IND1,IG)
   60       CONTINUE
            BSN=0.0
            DO 65 M=1,NPQ
            IF(DU(M).GT.0.0) THEN
               IND0=LY*LX*NSCT*IELEM**2+(M-1)*LY*IELEM+(K1-1)*IELEM+J0
               BSN=BSN+2.0*W(M)*FUNKNO(IND0,IG)*(1.0+ZCODE(2))
            ENDIF
   65       CONTINUE
            IF(FSN.NE.0.0) DFPN=DFPN*BSN/FSN
            DO 70 M=1,NPQ
            IF(DU(M).GT.0.0) THEN
               IND0=LY*LX*NSCT*IELEM**2+(M-1)*LY*IELEM+(K1-1)*IELEM+J0
               FUNKNO(IND0,IG)=FUNKNO(IND0,IG)+DFPN/(4.0*PI)
            ENDIF
   70       CONTINUE
   75       CONTINUE
         ENDIF
******XNEJ-
         IF((K1.EQ.1).AND.(ZCODE(3).NE.0.0)) THEN
            DO 95 I0=1,IELEM
            FSN=0.0
            DFPN=0.0
            SG=1.0
            DO 80 J0=1,IELEM
            IND0=KEYFLX(IR)+(I0-1)*IELEM+J0-1
            IND1=IDL(IR)+(I0-1)*IELEM+J0-1
            FSN=FSN+SG*SQRT(REAL(2*J0-1))*FUNKNO(IND0,IG)
            DFPN=DFPN+SG*SQRT(REAL(2*J0-1))*FUNSA(IND1,IG)
            SG=-SG
   80       CONTINUE
            BSN=0.0
            DO 85 M=1,NPQ
            IF(DE(M).LT.0.0) THEN
               IND0=LY*LX*NSCT*IELEM**2+IELEM*LY*NPQ+(M-1)*LX*IELEM+
     >         (K2-1)*IELEM+I0
               BSN=BSN+2.0*W(M)*FUNKNO(IND0,IG)*(1.0+ZCODE(3))
            ENDIF
   85       CONTINUE
            IF(FSN.NE.0.0) DFPN=DFPN*BSN/FSN
            DO 90 M=1,NPQ
            IF(DE(M).LT.0.0) THEN
               IND0=LY*LX*NSCT*IELEM**2+IELEM*LY*NPQ+(M-1)*LX*IELEM+
     >         (K2-1)*IELEM+I0
               FUNKNO(IND0,IG)=FUNKNO(IND0,IG)+DFPN/(4.0*PI)
            ENDIF
   90       CONTINUE
   95       CONTINUE
          ENDIF
*****XNEJ+
         IF((K1.EQ.LY).AND.(ZCODE(4).NE.0.0)) THEN
            DO 115 I0=1,IELEM
            FSN=0.0
            DFPN=0.0
            DO 100 J0=1,IELEM
            IND0=KEYFLX(IR)+(I0-1)*IELEM+J0-1
            IND1=IDL(IR)+(I0-1)*IELEM+J0-1
            FSN=FSN+SQRT(REAL(2*J0-1))*FUNKNO(IND0,IG)
            DFPN=DFPN+SQRT(REAL(2*J0-1))*FUNSA(IND1,IG)
  100       CONTINUE
            BSN=0.0
            DO 105 M=1,NPQ
            IF(DE(M).GT.0.0) THEN
               IND0=LY*LX*NSCT*IELEM**2+IELEM*LY*NPQ+(M-1)*LX*IELEM+
     >         (K2-1)*IELEM+I0
               BSN=BSN+2.0*W(M)*FUNKNO(IND0,IG)*(1.0+ZCODE(4))
            ENDIF
  105       CONTINUE
            IF(FSN.NE.0.0) DFPN=DFPN*BSN/FSN
            DO 110 M=1,NPQ
            IF(DE(M).GT.0.0) THEN
               IND0=LY*LX*NSCT*IELEM**2+IELEM*LY*NPQ+(M-1)*LX*IELEM+
     >         (K2-1)*IELEM+I0
               FUNKNO(IND0,IG)=FUNKNO(IND0,IG)+DFPN/(4.0*PI)
            ENDIF
  110       CONTINUE
  115       CONTINUE
         ENDIF
  150    CONTINUE
  160    CONTINUE   
  161    CONTINUE   
*----
*  UPGRADE THE SN SURFACE FLUX IN 3D CARTESIAN CASES
*----
      ELSE IF(ITYPE.EQ.7) THEN   
         CALL LCMLEN(IPTRK,'DU',NPQ,ITYLCM)
         CALL LCMGPD(IPTRK,'DU',DU_PTR)
         CALL LCMGPD(IPTRK,'DE',DE_PTR)
         CALL LCMGPD(IPTRK,'DZ',DZ_PTR)
         CALL LCMGPD(IPTRK,'W',W_PTR)         
         CALL C_F_POINTER(DU_PTR,DU,(/ NPQ /))
         CALL C_F_POINTER(DE_PTR,DE,(/ NPQ /))
         CALL C_F_POINTER(W_PTR,W,(/ NPQ /))
         CALL C_F_POINTER(DZ_PTR,DZ,(/ NPQ /))
         IR=0
         DO 182 K3=1,LZ
         DO 181 K2=1,LY
         DO 180 K1=1,LX
         IR=IR+1
         IF(MAT(IR).EQ.0) GO TO 180
         IF(VOL(IR).EQ.0.0) GO TO 180
******** XNEI-
         IF((K1.EQ.1).AND.(ZCODE(1).NE.0.0)) THEN 
            DO 202 J0=1,IELEM
            DO 201 K0=1,IELEM           
            FSN=0.0
            DFPN=0.0
            SG=1.0
            DO 210 I0=1,IELEM
            IND0=KEYFLX(IR)+(J0-1)*IELEM**2+(K0-1)*IELEM+I0-1
            IND1=IDL(IR)+(J0-1)*IELEM**2+(K0-1)*IELEM+I0-1
            FSN=FSN+SG*SQRT(REAL(2*I0-1))*FUNKNO(IND0,IG)
            DFPN=DFPN+SG*SQRT(REAL(2*I0-1))*FUNSA(IND1,IG)
            SG=-SG
  210       CONTINUE
            BSN=0.0           
            DO 220 M=1,NPQ
            IF(DU(M).LT.0.0) THEN
              IND0=LZ*LY*LX*NSCT*IELEM**3+(M-1)*LY*LZ*IELEM**2+
     >        LY*(K3-1)*IELEM**2+(K2-1)*IELEM**2+(J0-1)*IELEM+K0
              BSN=BSN+W(M)*FUNKNO(IND0,IG)*(1.0+ZCODE(1))
            ENDIF
  220       CONTINUE
            IF(FSN.NE.0.0) DFPN=DFPN*BSN/FSN
            DO 200 M=1,NPQ
            IF(DU(M).LT.0.0) THEN
              IND0=LZ*LY*LX*NSCT*IELEM**3+(M-1)*LY*LZ*IELEM**2+
     >        LY*(K3-1)*IELEM**2+(K2-1)*IELEM**2+(J0-1)*IELEM+K0
              FUNKNO(IND0,IG)=FUNKNO(IND0,IG)+DFPN/(4.0*PI)
            ENDIF
  200       CONTINUE
  201       CONTINUE
  202       CONTINUE
         ENDIF
******** XNEI+
         IF((K1.EQ.LX).AND.(ZCODE(2).NE.0.0)) THEN            
            DO 232 J0=1,IELEM
            DO 231 K0=1,IELEM
            FSN=0.0
            DFPN=0.0
            DO 240 I0=1,IELEM
            IND0=KEYFLX(IR)+(J0-1)*IELEM**2+(K0-1)*IELEM+I0-1
            IND1=IDL(IR)+(J0-1)*IELEM**2+(K0-1)*IELEM+I0-1   
            FSN=FSN+SQRT(REAL(2*I0-1))*FUNKNO(IND0,IG)
            DFPN=DFPN+SQRT(REAL(2*I0-1))*FUNSA(IND1,IG)
  240       CONTINUE
            BSN=0.0
            DO 250 M=1,NPQ
            IF(DU(M).GT.0.0) THEN
              IND0=LZ*LY*LX*NSCT*IELEM**3+(M-1)*LY*LZ*IELEM**2+
     >        LY*(K3-1)*IELEM**2+(K2-1)*IELEM**2+(J0-1)*IELEM+K0
              BSN=BSN+W(M)*FUNKNO(IND0,IG)*(1.0+ZCODE(2))
            ENDIF
  250       CONTINUE
            IF(FSN.NE.0.0) DFPN=DFPN*BSN/FSN
            DO 230 M=1,NPQ
            IF(DU(M).GT.0.0) THEN
              IND0=LZ*LY*LX*NSCT*IELEM**3+(M-1)*LY*LZ*IELEM**2+
     >        LY*(K3-1)*IELEM**2+(K2-1)*IELEM**2+(J0-1)*IELEM+K0
              FUNKNO(IND0,IG)=FUNKNO(IND0,IG)+DFPN/(4.0*PI)
            ENDIF
  230       CONTINUE
  231       CONTINUE
  232       CONTINUE
         ENDIF
***********XNEJ-
         IF((K2.EQ.1).AND.(ZCODE(3).NE.0.0)) THEN                    
            DO 262 I0=1,IELEM
            DO 261 K0=1,IELEM
            FSN=0.0
            DFPN=0.0
            SG=1.0
            DO 270 J0=1,IELEM
            IND0=KEYFLX(IR)+(I0-1)*IELEM*IELEM+(K0-1)*IELEM+J0-1
            IND1=IDL(IR)+(I0-1)*IELEM*IELEM+(K0-1)*IELEM+J0-1 
            FSN=FSN+SG*SQRT(REAL(2*J0-1))*FUNKNO(IND0,IG)
            DFPN=DFPN+SG*SQRT(REAL(2*J0-1))*FUNSA(IND1,IG)
            SG=-SG
  270       CONTINUE
            BSN=0.0
            DO 280 M=1,NPQ
            IF(DE(M).LT.0.0) THEN
              IND0=LZ*LY*LX*NSCT*IELEM**3+(M-1)*LX*LZ*IELEM**2+
     >        LX*(K3-1)*IELEM**2+(K1-1)*IELEM**2+(I0-1)*IELEM+K0
     >        +LY*LZ*NPQ*IELEM**2
              BSN=BSN+W(M)*FUNKNO(IND0,IG)*(1.0+ZCODE(3))
            ENDIF
  280       CONTINUE
            IF(FSN.NE.0.0) DFPN=DFPN*BSN/FSN
            DO 260 M=1,NPQ
            IF(DE(M).LT.0.0) THEN
              IND0=LZ*LY*LX*NSCT*IELEM**3+(M-1)*LX*LZ*IELEM**2+
     >        LX*(K3-1)*IELEM**2+(K1-1)*IELEM**2+(I0-1)*IELEM+K0
     >        +LY*LZ*NPQ*IELEM**2
              FUNKNO(IND0,IG)=FUNKNO(IND0,IG)+DFPN/(4.0*PI)
            ENDIF
  260       CONTINUE
  261       CONTINUE
  262       CONTINUE
         ENDIF
*******XNEJ +
         IF((K2.EQ.LY).AND.(ZCODE(4).NE.0.0)) THEN
            DO 292 I0=1,IELEM
            DO 291 K0=1,IELEM
            FSN=0.0
            DFPN=0.0
            DO 300 J0=1,IELEM
            IND0=KEYFLX(IR)+(I0-1)*IELEM*IELEM+(K0-1)*IELEM+J0-1
            IND1=IDL(IR)+(I0-1)*IELEM*IELEM+(K0-1)*IELEM+J0-1  
            FSN=FSN+SQRT(REAL(2*J0-1))*FUNKNO(IND0,IG)
            DFPN=DFPN+SQRT(REAL(2*J0-1))*FUNSA(IND1,IG)
  300       CONTINUE
            BSN=0.0
            DO 310 M=1,NPQ
            IF(DE(M).GT.0.0) THEN
              IND0=LZ*LY*LX*NSCT*IELEM**3+(M-1)*LX*LZ*IELEM**2+
     >        LX*(K3-1)*IELEM**2+(K1-1)*IELEM**2+(I0-1)*IELEM+K0
     >        +LY*LZ*NPQ*IELEM**2
              BSN=BSN+W(M)*FUNKNO(IND0,IG)*(1.0+ZCODE(4))
            ENDIF
  310       CONTINUE
            IF(FSN.NE.0.0) DFPN=DFPN*BSN/FSN
            DO 290 M=1,NPQ
            IF(DE(M).GT.0.0) THEN      
              IND0=LZ*LY*LX*NSCT*IELEM**3+(M-1)*LX*LZ*IELEM**2+
     >        LX*(K3-1)*IELEM**2+(K1-1)*IELEM**2+(I0-1)*IELEM+K0
     >        +LY*LZ*NPQ*IELEM**2
              FUNKNO(IND0,IG)=FUNKNO(IND0,IG)+DFPN/(4.0*PI)
            ENDIF
  290       CONTINUE
  291       CONTINUE
  292       CONTINUE
         ENDIF
********* XNEK -
         IF((K3.EQ.1).AND.(ZCODE(5).NE.0.0)) THEN
             DO 322 I0=1,IELEM
             DO 321 J0=1,IELEM
             FSN=0.0
             DFPN=0.0
             SG=1.0
             DO 330 K0=1,IELEM
             IND0=KEYFLX(IR)+(I0-1)*IELEM**2+(J0-1)*IELEM+K0-1
             IND1=IDL(IR)+(I0-1)*IELEM**2+(J0-1)*IELEM+K0-1
             FSN=FSN+SG*SQRT(REAL(2*K0-1))*FUNKNO(IND0,IG)
             DFPN=DFPN+SG*SQRT(REAL(2*K0-1))*FUNSA(IND1,IG)
             SG=-SG
  330        CONTINUE
             BSN=0.0
             DO 340 M=1,NPQ
             IF(DZ(M).LT.0.0) THEN
               IND0=LZ*LY*LX*NSCT*IELEM**3+LY*LZ*NPQ*IELEM**2+
     >         LX*LZ*NPQ*IELEM**2+(M-1)*LX*LY*IELEM**2+
     >         LX*(K2-1)*IELEM**2+(K1-1)*IELEM**2+(I0-1)*IELEM+J0
               BSN=BSN+W(M)*FUNKNO(IND0,IG)*(1.0+ZCODE(5))
             ENDIF
  340        CONTINUE
             IF(FSN.NE.0.0) DFPN=DFPN*BSN/FSN
             DO 320 M=1,NPQ
             IF(DZ(M).LT.0.0) THEN
               IND0=LZ*LY*LX*NSCT*IELEM**3+LY*LZ*NPQ*IELEM**2+
     >         LX*LZ*NPQ*IELEM**2+(M-1)*LX*LY*IELEM**2+
     >         LX*(K2-1)*IELEM**2+(K1-1)*IELEM**2+(I0-1)*IELEM+J0
               FUNKNO(IND0,IG)=FUNKNO(IND0,IG)+DFPN/(4.0*PI)
             ENDIF
  320        CONTINUE
  321        CONTINUE
  322        CONTINUE
         ENDIF
********** XNEK +
         IF((K3.EQ.LZ).AND.(ZCODE(6).NE.0.0)) THEN
             DO 352 I0=1,IELEM
             DO 351 J0=1,IELEM
             FSN=0.0
             DFPN=0.0
             DO 360 K0=1,IELEM                      
             IND0=KEYFLX(IR)+(I0-1)*IELEM**2+(J0-1)*IELEM+K0-1
             IND1=IDL(IR)+(I0-1)*IELEM**2+(J0-1)*IELEM+K0-1
             FSN=FSN+SQRT(REAL(2*K0-1))*FUNKNO(IND0,IG)
             DFPN=DFPN+SQRT(REAL(2*K0-1))*FUNSA(IND1,IG)
  360        CONTINUE
             BSN=0.0
             DO 370 M=1,NPQ
             IF(DZ(M).GT.0.0) THEN        
               IND0=LZ*LY*LX*NSCT*IELEM**3+LY*LZ*NPQ*IELEM**2+
     >         LX*LZ*NPQ*IELEM**2+(M-1)*LX*LY*IELEM**2+
     >         LX*(K2-1)*IELEM**2+(K1-1)*IELEM**2+(I0-1)*IELEM+J0
               BSN=BSN+W(M)*FUNKNO(IND0,IG)*(1.0+ZCODE(6))
             ENDIF
  370        CONTINUE
             IF(FSN.NE.0.0) DFPN=DFPN*BSN/FSN
             DO 350 M=1,NPQ
             IF(DZ(M).GT.0.0) THEN                      
               IND0=LZ*LY*LX*NSCT*IELEM**3+LY*LZ*NPQ*IELEM**2+
     >         LX*LZ*NPQ*IELEM**2+(M-1)*LX*LY*IELEM**2+
     >         LX*(K2-1)*IELEM**2+(K1-1)*IELEM**2+(I0-1)*IELEM+J0
               FUNKNO(IND0,IG)=FUNKNO(IND0,IG)+DFPN/(4.0*PI)
             ENDIF
  350        CONTINUE
  351        CONTINUE
  352        CONTINUE
         ENDIF
  180    CONTINUE
  181    CONTINUE
  182    CONTINUE
      ENDIF
*----
*  UPGRADE THE SOLUTION OF THE SN EQUATION.
*----
      DO 170 IR=1,NREG
      IF(MAT(IR).LE.0) GO TO 170
      DO 165 IEL=1,NLEG
      IND=KEYFLX(IR)+IEL-1
      JND=IDL(IR)+IEL-1
      OLD=FUNKNO(IND,IG)
      FUNKNO(IND,IG)=FUNKNO(IND,IG)+FUNSA(JND,IG)
 165  CONTINUE
 170  CONTINUE
 400  CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(SUNSA,FUNSA,IDL)
      RETURN
      END
