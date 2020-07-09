*DECK LIBTR2
      SUBROUTINE LIBTR2 (IPLIB,NAMFIL,NGRO,NBISO,NL,ISONAM,ISONRF,
     1 IPISO,ICOHNA,IINCNA,NTFG,TN,SN,SB,MASKI,NED,HVECT,ITIME,IMPX,
     2 NGF,NGFR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* transcription of the usefull interpolated microscopic cross section
* data from matxs to lcm data structures. Use matxs format from NJOY-91
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPLIB   pointer to the lattice microscopic cross section library
*         (L_LIBRARY signature)
* NAMFIL  name of the MATXS library file
* NGRO    number of energy groups
* NBISO   number of isotopes present in the calculation domain
* NL      number of Legendre orders required in the calculation
*         NL=1 or higher
* ISONAM  alias name of isotopes
* ISONRF  library name of isotopes
* IPISO   pointer array towards microlib isotopes
* ICOHNA  hcoh name
* IINCNA  hinc name
* NTFG    number of thermal groups where the thermal inelastic
*         correction is applied
* TN      temperature of each isotope
* SN      dilution cross section in each energy group of each
*         isotope. A value of 1.0E10 is used for infinite dilution
* SB      dilution cross section as used by livolant and jeanpierre
*         normalization
* MASKI   isotopic mask. Isotope with index I is processed if
*         MASKI(I)=.true.
* NED     number of extra vector edits from matxs
* HVECT   matxs names of the extra vector edits
*          MATXS reserved names:
*          NWT0/NWT1:    p0/p1 library weight function
*          NTOT0/NTOT1:  p0/p1 neutron total cross sections
*          NELAS:  neutron elastic scattering cross section
*          NINEL:  neutron inelastic scattering cross section
*          NG:     radiative capture cross section
*          NFTOT:  total fission cross section
*          NUDEL:  number of delayed secondary neutrons (nu-d)
*          CHID:  delayed fission spectrum
*          NF/NNF/N2NF/N3NF:  nu * partial fission cross sections
*          N2N/N3N/N4N:  (n,2n),(n,3n),(n,4n) cross sections
* ITIME   MATXS type of fission spectrum:
*         =1 steady-state; =2 prompt
* IMPX    print flag
*
*Parameters: output
* NGF     number of fast groups without self-shielding
* NGFR    number of fast and resonance groups
*
*Reference:
*  R. E. Macfarlane, TRANSX 2: A code for interfacing matxs cross-
*  section libraries to nuclear transport codes, Los Alamos National
*  Laboratory, Report LA-12312-MS, New Mexico, July 1992.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT CHARACTER*6 (H)
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER*(*) HVECT(NED),NAMFIL
      TYPE(C_PTR) IPLIB,IPISO(NBISO)
      INTEGER NGRO,NBISO,NL,ISONAM(3,NBISO),ISONRF(3,NBISO),
     1 ICOHNA(2,NBISO),IINCNA(2,NBISO),NTFG(NBISO),NED,ITIME,IMPX,
     2 NGF,NGFR
      LOGICAL MASKI(NBISO)
      REAL TN(NBISO),SN(NGRO,NBISO),SB(NGRO,NBISO)
*----
*  LOCAL VARIABLES
*----
      CHARACTER FORM*4,HSMG*131,HNISOR*12,HINC*6,HCOH*6,README*88,
     1 HNAMIS*12,TEXT12*12
      CHARACTER HN*8,HU*8,HS*8
      PARAMETER (MULT=2,IOUT=6,FORM='(A6)',MAXA=10000)
      TYPE(C_PTR) KPLIB
      LOGICAL LSUBM1,LTHERM,LTIME,LTERP
      DOUBLE PRECISION XHA(MAXA/2)
      REAL A(MAXA)
      INTEGER IA(MAXA),IHGAR(22)
      CHARACTER*6 HGAR(18)
      EQUIVALENCE (A(1),IA(1),XHA(1))
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IPR,ITYPRO
      REAL, ALLOCATABLE, DIMENSION(:) :: SFIS,SAVE,CHI,SIGF,TOTAL,FLUX,
     1 VECT,GAR,XS,TERP,TEMP,SIGZ
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SIGS
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: SCAT,XSMAT
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: LOGIED
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IPR(NBISO),ITYPRO(NL))
      ALLOCATE(SFIS(NGRO),SAVE(NGRO),CHI(NGRO),SIGS(NGRO,NL),SIGF(NGRO),
     1 TOTAL(NGRO),SCAT(NGRO,NGRO,NL),FLUX(NGRO),VECT(NGRO),GAR(NGRO),
     2 XSMAT(NGRO,NGRO,NL))
      ALLOCATE(LOGIED(NED))
*
      NGF=NGRO+1
      NGFR=0
      DO 100 I=1,NBISO
  100 IPR(I)=0
      IF(IMPX.GT.0) WRITE (IOUT,890) NAMFIL
      NIN=KDROPN(NAMFIL,2,2,0)
      IF(NIN.LE.0) THEN
         WRITE (HSMG,'(36HLIBTR2: UNABLE TO OPEN LIBRARY FILE ,A8,1H.)')
     1   NAMFIL
         CALL XABORT(HSMG)
      ENDIF
*----
*  INITIALIZE MATXS LIBRARY
*----
      NWDS=1+3*MULT
      IREC=1
*     --FILE IDENTIFICATION-----------
      CALL XDREED (NIN,IREC,A(1),NWDS)
*     --------------------------------
      WRITE(HN,'(A8)') XHA(1)
      WRITE(HU,'(A8)') XHA(2)
      WRITE(HS,'(A8)') XHA(3)
      IVER=IA(1+3*MULT)
      IF(IMPX.GT.0) WRITE (IOUT,935) HN,HU,HS,IVER
*
      NWDS=6
      IREC=2
*     --FILE CONTROL------------------
      CALL XDREED (NIN,IREC,A(1),NWDS)
*     --------------------------------
      NPART=IA(1)
      NTYPE=IA(2)
      NHOLL=IA(3)
      NMAT=IA(4)
      MAXW=IA(5)
*
      NWDS=NHOLL*MULT
      IF(NWDS.GT.MAXA)
     1 CALL XABORT('LIBTR2: INSUFFICIENT VALUE OF MAXA(1).')
      IREC=3
*     --HOLLERITH IDENTIFICATION------
      CALL XDREED (NIN,IREC,A(1),NWDS)
*     --------------------------------
      WRITE(README(9:),'(6H FROM ,9A8)') (XHA(I),I=1,MIN(NHOLL,9))
      IF(IMPX.GT.0) WRITE (IOUT,'(1X,9A8)') (XHA(I),I=1,MIN(NHOLL,9))
*
      NWDS=(NPART+NTYPE+NMAT)*MULT+2*NTYPE+NPART+2*NMAT
      IF(NWDS.GT.MAXA)
     1 CALL XABORT('LIBTR2: INSUFFICIENT VALUE OF MAXA(2).')
      IREC=4
*     --FILE DATA---------------------
      CALL XDREED (NIN,IREC,A(1),NWDS)
*     --------------------------------
      IF((NWDS/2)*2.NE.NWDS) NWDS=NWDS+1
      L2=1+NWDS
      L2H=1+NWDS/MULT
*----
*  CHECK GROUP STRUCTURES
*----
      NEX1=(NPART+NTYPE+NMAT)*MULT
      DO 120 I=1,NPART
      WRITE(HPRT,FORM) XHA(I)
      CALL LIBCOV(HPRT)
      NG=IA(NEX1+I)
      IF((HPRT.EQ.'N').AND.(NG.NE.NGRO))
     1 CALL XABORT('LIBTR2: INCONSISTENT GROUP STRUCTURES.')
      NWDS=NG+1
      ALLOCATE(XS(NWDS))
      IREC=IREC+1
*     --GROUP STRUCTURE-------------
      CALL XDREED (NIN,IREC,XS,NWDS)
*     ------------------------------
      IF(HPRT.EQ.'N') THEN
*        ENERGY BOUND IN EACH GROUP (IN EV):
         CALL LCMPUT(IPLIB,'ENERGY',NGRO+1,2,XS)
         DO 110 J=1,NGRO
  110    VECT(J)=LOG(XS(J)/XS(J+1))
         CALL LCMPUT(IPLIB,'DELTAU',NGRO,2,VECT)
      ENDIF
      DEALLOCATE(XS)
  120 CONTINUE
*----
*  READ THROUGH MATXS FILE AND ACCUMULATE CROSS SECTIONS FOR THIS RANGE
*  MATS, LEGENDRE ORDERS, AND GROUPS.
*
*  ***MATERIAL/ISOTOPE LOOP***
*----
      IRZM=IREC+1
      DO 840 IM=1,NMAT
  130 CNORM=0.0
      DO 150 IG1=1,NGRO
      CHI(IG1)=0.0
      SIGF(IG1)=0.0
      TOTAL(IG1)=0.0
      DO 150 IL=1,NL
      SIGS(IG1,IL)=0.0
      DO 150 IG2=1,NGRO
  150 SCAT(IG1,IG2,IL)=0.0
      WRITE (HMAT,FORM) XHA(NPART+NTYPE+IM)
      DO 160 IMX=1,NBISO
      IF(MASKI(IMX)) THEN
         IMT=IMX
         WRITE(HNAMIS,'(3A4)') (ISONAM(ITC,IMX),ITC=1,3)
         WRITE(HNISOR,'(3A4)') (ISONRF(ITC,IMX),ITC=1,3)
         WRITE(HCOH,'(A4,A2)') (ICOHNA(ITC,IMX),ITC=1,2)
         WRITE(HINC,'(A4,A2)') (IINCNA(ITC,IMX),ITC=1,2)
         CALL LIBCOV(HCOH)
         CALL LIBCOV(HINC)
         IF((HMAT.EQ.HNISOR(:6)).AND.(IPR(IMX).EQ.0)) GO TO 170
      ENDIF
  160 CONTINUE
      GO TO 840
*----
*  RECOVER THE MATERIAL CONTROL
*----
  170 IPR(IMT)=1
      DO 174 IED=1,NED
  174 LOGIED(IED)=.FALSE.
      KPLIB=IPISO(IMT) ! set IMT-th isotope
      LOC=(NPART+NTYPE+NMAT)*MULT+NPART+2*NTYPE+IM
      NSUB=IA(LOC)
      LOCM=IA(LOC+NMAT)
      IREC=LOCM+IRZM
      NWDS=MULT+1+6*NSUB
      IF(L2+NWDS-1.GT.MAXA)
     1 CALL XABORT('LIBTR2: INSUFFICIENT VALUE OF MAXA(3).')
*     --MATERIAL CONTROL---------------
      CALL XDREED (NIN,IREC,A(L2),NWDS)
*     ---------------------------------
*     MASS RATIO OF EACH MATERIAL/ISOTOPE IN THE CALCULATION DOMAIN:
      AWR=A(L2+MULT)
      IF((NWDS/2)*2.NE.NWDS) NWDS=NWDS+1
      L3=L2+NWDS
      L3H=L2H+NWDS/MULT
      ALLOCATE(TERP(NSUB*NGRO),TEMP(NSUB),SIGZ(NSUB))
      DO 175 ISUBM=1,NSUB
      TEMP(ISUBM)=A(L2+MULT+6*(ISUBM-1)+1)
  175 SIGZ(ISUBM)=A(L2+MULT+6*(ISUBM-1)+2)
      NSUB0=0
      DO 185 ITYPE=1,NTYPE
      N0=0
      DO 180 ISUBM=1,NSUB
      IF(IA(L2+MULT+6*(ISUBM-1)+3).EQ.ITYPE) N0=N0+1
  180 CONTINUE
      CALL LIBTE2(NGRO,N0,TEMP(NSUB0+1),SIGZ(NSUB0+1),TN(IMT),
     1 SN(1,IMT),TERP(NSUB0*NGRO+1))
  185 NSUB0=NSUB0+N0
      IF(NSUB0.NE.NSUB) CALL XABORT('LIBTR2: DATA TYPE FAILURE.')
      DEALLOCATE(SIGZ,TEMP)
*----
*  TEMPERATURE AND BACKGROUND LOOP
*----
      IOLDTY=0
      L5=0
      DNORM=0.0
      DO 720 ISUBM=1,NSUB
      LOC=L2+MULT+6*(ISUBM-1)
      TMAT=A(LOC+1)
      SMAT=A(LOC+2)
      ITYPE=IA(LOC+3)
      LSUBM1=(ITYPE.NE.IOLDTY)
      IOLDTY=ITYPE
      N1D=IA(LOC+4)
      N2D=IA(LOC+5)
      LOCS=IA(LOC+6)
      LOCG=(NPART+NTYPE+NMAT)*MULT
      JINP=IA(LOCG+NPART+ITYPE)
      NING=IA(LOCG+JINP)
      JOUTP=IA(LOCG+NPART+NTYPE+ITYPE)
      NOUTG=IA(LOCG+JOUTP)
      WRITE(HPRT,FORM) XHA(JINP)
      CALL LIBCOV(HPRT)
      WRITE(HTYPE,FORM) XHA(NPART+ITYPE)
      CALL LIBCOV(HTYPE)
      IF(IMPX.GT.6) WRITE(IOUT,860) ISUBM,HPRT,HTYPE,HMAT
      IF((HTYPE.NE.'NSCAT').AND.(HTYPE.NE.'NTHERM')) GO TO 720
      IF(.NOT.LSUBM1) THEN
         LTERP=.TRUE.
         DO 190 IK=1,NGRO
  190    LTERP=LTERP.AND.(TERP(NGRO*(ISUBM-1)+IK).EQ.0.0)
         IF(LTERP) GO TO 720
      ENDIF
*----
*  PROCESS THIS SUBMATERIAL
*----
      JREC=IREC+LOCS
      IF(N1D.EQ.0) GO TO 460
      NWDS=(2+MULT)*N1D
      IF(L3+NWDS-1.GT.MAXA)
     1 CALL XABORT('LIBTR2: INSUFFICIENT VALUE OF MAXA(4).')
      JREC=JREC+1
*     --VECTOR CONTROL-----------------
      CALL XDREED (NIN,JREC,A(L3),NWDS)
*     ---------------------------------
      NEX1=L3-1+MULT*N1D
      NEX2=NEX1+N1D
      IF(LSUBM1.AND.(IMPX.GT.4)) THEN
         WRITE (IOUT,870) HTYPE,HMAT,(XHA(L3H+IR-1),IR=1,N1D)
      ENDIF
*----
*  VECTOR PARTIALS
*----
      IF(LSUBM1) THEN
         DO 210 KG=1,NGRO
         SFIS(KG)=0.0
  210    SAVE(KG)=0.0
      ENDIF
*----
*  LOOP OVER REACTIONS
*----
      IRMAX=0
      DO 455 IR=1,N1D
      IF(IR.GT.IRMAX) THEN
*        MANY VECTORS (REACTIONS) ARE STORED IN VECTOR BLOCK.
         NWDS=0
         IJ0=IRMAX+1
         DO 220 IJ=IJ0,N1D
         NW=IA(NEX2+IJ)-IA(NEX1+IJ)+1
         IF(NWDS+NW.GE.MAXW) GO TO 230
         IRMAX=IRMAX+1
  220    NWDS=NWDS+NW
  230    IF(NWDS.EQ.0) CALL XABORT('LIBTR2: MAXW IS TOO SMALL.')
         ALLOCATE(XS(NWDS))
         JREC=JREC+1
*        --VECTOR BLOCK----------------
         CALL XDREED (NIN,JREC,XS,NWDS)
*        ------------------------------
         L5=0
      ENDIF
      WRITE(HVPS,FORM) XHA(L3H-1+IR)
      CALL LIBCOV(HVPS)
      NK=IA(NEX2+IR)-IA(NEX1+IR)+1
*----
*  SAVE REQUIRED EXTRA EDIT.
*----
      DO 260 I=1,NED
      TEXT12=HVECT(I)
      CALL LIBCOV(TEXT12)
      IF(HVPS.EQ.TEXT12) THEN
         IF(LSUBM1) THEN
            DO 240 IK=1,NGRO
  240       VECT(IK)=0.0
         ELSE
            CALL LCMGET(KPLIB,HVECT(I),VECT)
         ENDIF
         DO 250 IK=1,NK
         IF(XS(L5+IK).EQ.0.0) GO TO 250
         JJ=IA(NEX1+IR)+IK-1
         TERPZ=1.0
         IF(.NOT.LSUBM1) TERPZ=TERP(NGRO*(ISUBM-1)+JJ)
         VECT(JJ)=VECT(JJ)+TERPZ*XS(L5+IK)
  250    CONTINUE
         LOGIED(I)=.TRUE.
         CALL LCMPUT(KPLIB,HVECT(I),NGRO,2,VECT)
         GO TO 270
      ENDIF
  260 CONTINUE
*----
*  SAVE MODEL WEIGHT FUNCTIONS
*----
  270 IF((HTYPE.EQ.'NSCAT').AND.(HVPS.EQ.'NWT0').AND.LSUBM1) THEN
         DO 280 IK=1,NK
         JJ=IA(NEX1+IR)+IK-1
  280    FLUX(JJ)=XS(L5+IK)
         GO TO 450
      ENDIF
      IF((HTYPE.EQ.'NTHERM').AND.(HVPS.NE.HINC).AND.
     1   (HVPS.NE.HCOH)) GO TO 450
*----
*  LOOP OVER GROUPS
*----
      DO 440 IK=1,NK
      IF(XS(L5+IK).EQ.0.0) GO TO 440
      JJ=IA(NEX1+IR)+IK-1
      LTHERM=(JJ.GE.NGRO-NTFG(IMT)+1)
      LTIME=(ITIME.EQ.1)
*----
*  INTERPOLATION FACTOR
*----
      TERPZ=1.0
      IF(.NOT.LSUBM1) TERPZ=TERP(NGRO*(ISUBM-1)+JJ)
      IF((SMAT.LT.0.9E10).AND.(ABS(XS(L5+IK)).GT.1.0E-6).AND.
     1 (.NOT.LSUBM1).AND.(HVPS.EQ.'NTOT0')) THEN
         NGF=MIN(NGF,JJ-1)
         NGFR=MAX(NGFR,JJ)
      ENDIF
      IF(ABS(TERPZ).LT.1.0E-3) GO TO 440
      ADD=TERPZ*XS(L5+IK)
*
      IF(HVPS.EQ.'NTOT0') THEN
*        TOTAL XSEC
         TOTAL(JJ)=TOTAL(JJ)+ADD
      ELSE IF((.NOT.LSUBM1).AND.(HVPS.EQ.'NFTOT')) THEN
*        FISSION CROSS SECTION
         SIGF(JJ)=SIGF(JJ)+ADD*SAVE(JJ)
      ELSE IF(LSUBM1.AND.(HVPS.EQ.'NFTOT')) THEN
         SFIS(JJ)=SFIS(JJ)+ADD
      ELSE IF(LSUBM1.AND.LTIME.AND.(HVPS.EQ.'NUDEL')) THEN
*        DELAYED FISSION
         SIGF(JJ)=SIGF(JJ)+ADD*SFIS(JJ)
         SAVE(JJ)=SAVE(JJ)+SFIS(JJ)*ADD
         IF(IK.EQ.1) DNORM=0.0
         DNORM=DNORM+ADD*SFIS(JJ)*FLUX(JJ)
      ELSE IF(LSUBM1.AND.LTIME.AND.(HVPS.EQ.'CHID')) THEN
*        DELAYED FISSION
         IF(DNORM.EQ.0.0) THEN
            WRITE (HSMG,1060) HMAT
            CALL XABORT(HSMG)
         ENDIF
         ADDD=DNORM*XS(L5+IK)
         CNORM=CNORM+ADDD
         CHI(JJ)=CHI(JJ)+ADDD
      ENDIF
  440 CONTINUE
*
* END OF REACTION LOOP
  450 L5=L5+NK
      IF(L5.EQ.NWDS) DEALLOCATE(XS)
  455 CONTINUE
*----
*  RECOVER SCATTERING MATRIX CONTROL INFORMATION.
*----
  460 IF(N2D.EQ.0) GO TO 720
      DO 700 K=1,N2D
      NWDS=MULT+2+2*NOUTG
      IF(L3+NWDS-1.GT.MAXA)
     1 CALL XABORT('LIBTR2: INSUFFICIENT VALUE OF MAXA(5).')
      JREC=JREC+1
*     --MATRIX CONTROL-----------------
      CALL XDREED (NIN,JREC,A(L3),NWDS)
*     ---------------------------------
      LORD=IA(L3+MULT)
      IF(LORD.EQ.0) GO TO 700
      WRITE(HMTX,FORM) XHA(L3H)
      HMTX2=HMTX
      CALL LIBCOV(HMTX)
      LN=L3+MULT+1
      LG=LN+NOUTG
      IFISN=0
      IF(HTYPE.EQ.'NSCAT'.AND.(HMTX.EQ.'NF'.OR.HMTX.EQ.'NNF'
     1   .OR.HMTX.EQ.'N2NF'.OR.HMTX.EQ.'N3NF')) IFISN=1
      IF(HTYPE.EQ.'NSCAT'.AND.HMTX.EQ.'NFTOT') IFISN=2
      JCONST=IA(L3+MULT+1)
*----
*  RECOVER A NEW SCATTERING MATRIX SUB-BLOCK.
*----
      IF(NING.NE.NOUTG) CALL XABORT('LIBTR2: ONLY (N,N) ALLOWED.')
      DO 465 IL=1,NL
      DO 465 JJ=1,NOUTG
      DO 465 JJP=1,NING
  465 XSMAT(JJP,JJ,IL)=0.0
      NOUMAX=0
  470 NWDS=0
      NOUMIN=NOUMAX+1
      DO 475 JJ=NOUMIN,NOUTG
      NW=IA(LN+JJ)*LORD
      IF(NWDS+NW.GE.MAXW) GO TO 480
      NOUMAX=NOUMAX+1
  475 NWDS=NWDS+NW
  480 IF(NWDS.GT.0) THEN
         ALLOCATE(XS(NWDS))
         JREC=JREC+1
*        --MATRIX SUB-BLOCK------------
         CALL XDREED (NIN,JREC,XS,NWDS)
*        ------------------------------
         L5=0
      ELSE
         GO TO 520
      ENDIF
      DO 510 JJ=NOUMIN,NOUMAX
      NP=IA(LN+JJ)
      IF(NP.EQ.0) GO TO 510
      DO 500 IL=1,LORD
      IF(IL.GT.NL) GO TO 500
      DO 490 IP=1,NP
      JJP=IA(LG+JJ)-IP+1
  490 XSMAT(JJP,JJ,IL)=XS(L5+IP+NP*(IL-1))
  500 CONTINUE
  510 L5=L5+NP*LORD
      DEALLOCATE(XS)
  520 IF(NOUMAX.LT.NOUTG) GO TO 470
      IF(JCONST.NE.0) THEN
         IF(LORD.GT.1) CALL XABORT('LIBTR2: INVALID DATA ON MATXS2.')
         NWDS=NOUTG+JCONST
         ALLOCATE(XS(NWDS))
         JREC=JREC+1
*        --CONSTANT SUB-BLOCK----------
         CALL XDREED (NIN,JREC,XS,NWDS)
*        ------------------------------
         L5=0
         DO 530 JJ=1,NOUTG
         SPEC=XS(L5+JJ)
         JJP0=NING-JCONST+1
         DO 530 JJP=JJP0,NING
  530    XSMAT(JJP,JJ,1)=XSMAT(JJP,JJ,1)+SPEC*XS(L5+NOUTG+JJP-JJP0+1)
         DEALLOCATE(XS)
      ENDIF
*----
*  STORE DESIRED CROSS SECTIONS
*----
      IF(HTYPE.EQ.'NTHERM'.AND.HMTX.NE.HINC.AND.
     1   HMTX.NE.HCOH) GO TO 670
*----
*  LOOP OVER SINK, ORDER, SOURCE
*----
      DO 650 JJ=1,NOUTG
      DO 650 IL=1,LORD
      IF(IL.GT.NL) GO TO 650
      DO 640 JJP=1,NING
      XSNOW=XSMAT(JJP,JJ,IL)
      IF(XSNOW.EQ.0.) GO TO 640
*----
*  INTERPOLATION FACTOR
*----
      TERPZ=1.0
      IF(.NOT.LSUBM1) TERPZ=TERP(NGRO*(ISUBM-1)+JJP)
      IF(ABS(TERPZ).LT.1.0E-3) GO TO 640
      XSEC=TERPZ*XSNOW
*----
*  CHECK FOR SCATTERING AND FISSION MATRICES
*----
      IF(IFISN.EQ.0) THEN
*        THERMAL CORRECTION TO SCATTERING MATRIX
         IF((HMTX.EQ.'NELAS').AND.(JJP.GE.NGRO-NTFG(IMT)+1)) THEN
            IF(IL.EQ.1) TOTAL(JJP)=TOTAL(JJP)-XSEC
            GO TO 640
         ENDIF
         IF(((HMTX.EQ.HINC).OR.(HMTX.EQ.HCOH)).AND.(JJP.LT.
     1   NGRO-NTFG(IMT)+1)) GO TO 640
*        TOTAL SCATTERING MATRIX
*        SCAT(SECONDARY,PRIMARY,ORDER+1)
         SCAT(JJ,JJP,IL)=SCAT(JJ,JJP,IL)+XSEC
*        TOTAL XS AND TOTAL SCATTERING VECTOR
         SIGS(JJP,IL)=SIGS(JJP,IL)+XSEC
         IF((IL.EQ.1).AND.(JJP.GE.NGRO-NTFG(IMT)+1)) THEN
            TOTAL(JJP)=TOTAL(JJP)+XSEC
         ENDIF
      ELSE IF((IL.EQ.1).AND.(IFISN.NE.0).AND.(HTYPE.EQ.'NSCAT')) THEN
*        FISSION VECTORS
         SIGF(JJP)=SIGF(JJP)+XSEC
         CNORM=CNORM+XSEC*FLUX(JJP)
         CHI(JJ)=CHI(JJ)+XSEC*FLUX(JJP)
      ENDIF
  640 CONTINUE
  650 CONTINUE
*----
*  ACCUMULATE FISSION NUBAR
*----
  670 IF(LSUBM1.AND.(IFISN.NE.0).AND.(HTYPE.EQ.'NSCAT')) THEN
         DO 680 JJ=1,NOUTG
         DO 680 JJP=1,NING
  680    SAVE(JJP)=SAVE(JJP)+XSMAT(JJP,JJ,1)
      ENDIF
*
      HGAR(MOD(K-1,18)+1)=HMTX2
      IF((K.EQ.1).AND.LSUBM1.AND.(IMPX.GT.4)) THEN
         WRITE (IOUT,880) HTYPE,HMAT
      ENDIF
      IF((MOD(K-1,18).EQ.17).AND.LSUBM1.AND.(IMPX.GT.4)) THEN
         WRITE (IOUT,885) (HGAR(I)//' ',I=1,18)
      ELSE IF((K.EQ.N2D).AND.LSUBM1.AND.(IMPX.GT.4)) THEN
         WRITE (IOUT,885) (HGAR(I)//' ',I=1,MOD(N2D-1,18)+1)
      ENDIF
  700 CONTINUE
*----
*  SAVE FISSION NU FOR SHIELDING TERMS
*----
      IF(LSUBM1.AND.(HTYPE.EQ.'NSCAT')) THEN
         DO 710 JJ=1,NGRO
         IF(SFIS(JJ).EQ.0) GO TO 710
         SAVE(JJ)=SAVE(JJ)/SFIS(JJ)
  710    CONTINUE
      ENDIF
*
* END OF SUBMATERIAL LOOP.
  720 CONTINUE
      DEALLOCATE(TERP)
*----
*  PRINT FINAL FLUX COMPONENTS
*----
      IF((IMPX.GT.6).AND.MASKI(IMT)) THEN
         SUM=0.0
         DO 730 JJ=1,NGRO
  730    SUM=SUM+FLUX(JJ)
         WRITE(IOUT,927) HNAMIS,SUM
         WRITE(IOUT,928) (FLUX(I),I=1,NGRO)
      ENDIF
*----
*  PERFORM LIVOLANT-JEANPIERRE NORMALIZATION AND SAVE CROSS SECTION
*  INFORMATION ON LCM.
*----
      DO 740 I=1,NGRO
      IF((SN(I,IMT).NE.SB(I,IMT)).AND.(SN(I,IMT).LT.1.0E10)) THEN
         VECT(I)=1.0/(1.0+(TOTAL(I)-SIGS(I,1))*(1.0/SN(I,IMT)
     1   -1.0/SB(I,IMT)))
      ELSE
         VECT(I)=1.0
      ENDIF
      IF(SN(I,IMT).LT.1.0E10) THEN
         FLUX(I)=SN(I,IMT)/(SN(I,IMT)+TOTAL(I)-SIGS(I,1))/VECT(I)
      ELSE
         FLUX(I)=1.0
      ENDIF
  740 TOTAL(I)=TOTAL(I)*VECT(I)
      IF(IMPX.GT.5) THEN
         WRITE(IOUT,920) HNAMIS
         WRITE(IOUT,928) (VECT(I),I=1,NGRO)
      ENDIF
      CALL LCMPUT(KPLIB,'NTOT0',NGRO,2,TOTAL)
      CALL LCMPUT(KPLIB,'NWT0',NGRO,2,FLUX)
      DO 750 IL=0,NL-1
      DO 750 IG2=1,NGRO
      FACTOR=VECT(IG2)
      SIGS(IG2,IL+1)=SIGS(IG2,IL+1)*FACTOR
      DO 750 IG1=1,NGRO
  750 SCAT(IG1,IG2,IL+1)=SCAT(IG1,IG2,IL+1)*FACTOR
*
      CALL XDRLGS(KPLIB,1,0,0,NL-1,1,NGRO,SIGS,SCAT,ITYPRO)
*
      DO 810 IED=1,NED
      TEXT12=HVECT(IED)
      CALL LIBCOV(TEXT12)
      IF(LOGIED(IED).AND.(TEXT12(:3).NE.'CHI')
     1              .AND.(TEXT12(:2).NE.'NU')
     2              .AND.(TEXT12.NE.'NTOT0')
     3              .AND.(TEXT12(:3).NE.'NWT')) THEN
         CALL LCMGET(KPLIB,HVECT(IED),GAR)
         DO 800 I=1,NGRO
  800    GAR(I)=GAR(I)*VECT(I)
         CALL LCMPUT(KPLIB,HVECT(IED),NGRO,2,GAR)
      ENDIF
  810 CONTINUE
*
      IF(CNORM.NE.0.0) THEN
*        FISSION SOURCE NORMALIZATION
         DO 820 JJ=1,NGRO
         CHI(JJ)=CHI(JJ)/CNORM
  820    SIGF(JJ)=SIGF(JJ)*VECT(JJ)
         CALL LCMPUT(KPLIB,'NUSIGF',NGRO,2,SIGF)
         CALL LCMPUT(KPLIB,'CHI',NGRO,2,CHI)
      ENDIF
      CALL LCMPTC(KPLIB,'ALIAS',12,1,HNAMIS)
      CALL LCMPUT(KPLIB,'AWR',1,2,AWR)
      WRITE(README(:8),'(A8)') HNAMIS(1:8)
      READ(README,'(22A4)') (IHGAR(I),I=1,22)
      CALL LCMPUT(KPLIB,'README',22,3,IHGAR)
      GO TO 130
*----
*  END OF MATERIAL/ISOTOPE LOOP.
*----
  840 CONTINUE
*----
*  CLOSE MATXS FILE.
*----
      IER=KDRCLS(NIN,1)
      IF(IER.LT.0) THEN
         WRITE (HSMG,'(37HLIBTR2: UNABLE TO CLOSE LIBRARY FILE ,A8,1H.
     1   )') NAMFIL
         CALL XABORT(HSMG)
      ENDIF
*----
*  CHECK IF ALL NBISO ISOTOPES HAVE BEEN PROCESSED.
*----
      NISOT=0
      DO 850 IMT=1,NBISO
      IF(MASKI(IMT)) THEN
         IF(IPR(IMT).EQ.0) THEN
            WRITE (IOUT,910) (ISONAM(ITC,IMT),ITC=1,3),NAMFIL
            NISOT=NISOT+1
         ENDIF
      ENDIF
  850 CONTINUE
      IF(NISOT.GT.0) CALL XABORT('LIBTR2: MISSING ISOTOPES')
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(LOGIED)
      DEALLOCATE(XSMAT,GAR,VECT,FLUX,SCAT,TOTAL,SIGF,SIGS,CHI,SAVE,SFIS)
      DEALLOCATE(ITYPRO,IPR)
      RETURN
*
  860 FORMAT(/31H LIBTR2: PROCESSING SUBMATERIAL,I5,5X,12HINCIDENT PAR,
     1 6HTICLE=,A6,5X,10HDATA TYPE=,A6,5X,9HMATERIAL=,A6)
  870 FORMAT(/52H AVAILABLE IDENTIFIERS OF REACTION VECTORS FOR TYPE ,
     1 A6,14H AND MATERIAL ,A6,1H:/(1X,18A7))
  880 FORMAT(/53H AVAILABLE IDENTIFIERS OF REACTION MATRICES FOR TYPE ,
     1 A6,14H AND MATERIAL ,A6,1H:)
  885 FORMAT(1X,18A7)
  890 FORMAT(/33H PROCESSING MATXS2 LIBRARY NAMED ,A8,1H.)
  910 FORMAT(/27H LIBTR2: MATERIAL/ISOTOPE ',3A4,16H' IS MISSING ON ,
     1 16HMATXS FILE NAME ,A8,1H.)
  920 FORMAT(/40H L-J NORMALIZATION FACTORS FOR MATERIAL ,A12)
  927 FORMAT(/19H FLUX FOR MATERIAL ,A12,7H   SUM=,1P,E12.5)
  928 FORMAT(1X,1P,10E12.4)
  935 FORMAT(/17H MATXS2 FILE ID: ,3A8,6H VERS ,I2)
 1060 FORMAT(35HLIBTR2: DNORM MISSING FOR MATERIAL ,A6,1H.)
      END
