*DECK KINSRC
      SUBROUTINE KINSRC(IPTRK,IPSYS,CMOD,IMPX,IFL,IPR,IEXP,NGR,NBM,
     1 NBFIS,NDG,ITY,LL4,NUN,NUP,PDC,TTF,TTP,DT,OVR,CHI,CHD,SGF,SGD,
     2 OMEGA,EVECT,PC,SRC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* compute the space-time kinetics source for a known neutron flux.
*
*Copyright:
* Copyright (C) 2010 Ecole Polytechnique de Montreal.
*
*Author(s): A. Hebert
*
*Parameters: input
* IPTRK  pointer to L_TRACK object.
* IPSYS  pointer to L_SYSTEM object.
* CMOD   name of the assembly door (BIVAC or TRIVAC).
* IMPX   print parameter (equal to zero for no print).
* IFL    integration scheme for fluxes: =1 implicit;
*        =2 Crank-Nicholson; =3 theta.
* IPR    integration scheme for precursors: =1 implicit;
*        =2 Crank-Nicholson; =3 theta; =4 exponential.
* IEXP   exponential transformation flag (=1 to activate).
* NGR    number of energy groups.
* NBM    number of material mixtures.
* NBFIS  number of fissile isotopes.
* NDG    number of delayed-neutron groups.
* ITY    type of solution: =1: classical Bivac/diffusion;
*        =2: classical Trivac/diffusion; =3 Raviart-Thomas in
*        Trivac/diffusion; =11: Bivac/SPN; =13 Trivac/SPN.
* LL4    order of the system matrices.
* NUN    total number of unknowns per energy group.
* NUP    total number of precursor unknowns per precursor group.
* PDC    precursor decay constants.
* TTF    value of theta-parameter for fluxes.
* TTP    value of theta-parameter for precursors.
* DT     current time increment.
* OVR    reciprocal neutron velocities/DT.
* CHI    steady-state fission spectrum.
* CHD    delayed fission spectrum
* SGF    nu*fission macroscopic x-sections/keff.
* SGD    delayed nu*fission macroscopic x-sections/keff.
* OMEGA  exponential transformation parameter.
* EVECT  neutron flux.
* PC     precursor concentrations.
*
*Parameters: output
* SRC    space-time kinetics source.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,IPSYS
      INTEGER IMPX,IFL,IPR,IEXP,NGR,NBM,NBFIS,NDG,ITY,LL4,NUN,NUP
      REAL PDC(NDG),TTF,TTP,DT,OVR(NBM,NGR),CHI(NBM,NBFIS,NGR),
     1 CHD(NBM,NBFIS,NGR,NDG),SGF(NBM,NBFIS,NGR),SGD(NBM,NBFIS,NGR,NDG),
     2 OMEGA(NBM,NGR),EVECT(NUN,NGR),PC(NUP,NDG,NBFIS)
      DOUBLE PRECISION SRC(NUN,NGR)
      CHARACTER CMOD*12
*----
*  LOCAL VARIABLES
*----
      PARAMETER(IOS=6)
      DOUBLE PRECISION DTF,DTP,DARG,DK,DSM
      LOGICAL LFIS
      CHARACTER TEXT12*12
      REAL, DIMENSION(:), ALLOCATABLE :: WORK1,WORK2,CHEXP
      REAL, DIMENSION(:), POINTER :: AGAR
      TYPE(C_PTR) AGAR_PTR
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(WORK1(LL4),WORK2(NBM),CHEXP(NBM))
*
      DTF=9999.0D0
      DTP=9999.0D0
      IF(IFL.EQ.1)THEN
        DTF=1.0D0
      ELSEIF(IFL.EQ.2)THEN
        DTF=0.5D0
      ELSEIF(IFL.EQ.3)THEN
        DTF=DBLE(TTF)
      ENDIF
      IF(IPR.EQ.2)THEN
        DTP=0.5D0
      ELSEIF(IPR.EQ.3)THEN
        DTP=DBLE(TTP)
      ENDIF
*
      IF(IMPX.GT.0)WRITE(IOS,1001)CMOD
      CALL XDDSET(SRC,NUN*NGR,0.0D0)
      DO 110 IGR=1,NGR
      WRITE(TEXT12,'(1HA,2I3.3)') IGR,IGR
      CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,EVECT(1,IGR),WORK1)
      DO 10 IND=1,LL4
      SRC(IND,IGR)=-(1.0D0-DTF)*WORK1(IND)
   10 CONTINUE
      DO 20 JGR=1,NGR
      IF(JGR.EQ.IGR) GO TO 20
      WRITE(TEXT12,'(1HA,2I3.3)') IGR,JGR
      CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
      IF(ILONG.EQ.0) GO TO 20
      IF((CMOD.EQ.'BIVAC').OR.(ITY.EQ.13))THEN
        CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,EVECT(1,JGR),WORK1)
        DO 14 IND=1,LL4
        SRC(IND,IGR)=SRC(IND,IGR)+(1.0D0-DTF)*WORK1(IND)
   14   CONTINUE
      ELSEIF(CMOD.EQ.'TRIVAC')THEN
        CALL LCMGPD(IPSYS,TEXT12,AGAR_PTR)
        CALL C_F_POINTER(AGAR_PTR,AGAR,(/ ILONG /))
        DO 15 IND=1,ILONG
        SRC(IND,IGR)=SRC(IND,IGR)+(1.0D0-DTF)*AGAR(IND)*EVECT(IND,JGR)
   15   CONTINUE
      ENDIF
   20 CONTINUE
*----
*  PRECURSOR CONTRIBUTION
*----
      DO 93 IFIS=1,NBFIS
      DO 50 IDG=1,NDG
      DARG=PDC(IDG)*DT
      IF(IPR.EQ.1)THEN
        DK=1.0D0/(1.0D0+DARG)
      ELSEIF(IPR.EQ.4)THEN
        DK=DEXP(-DARG)
      ELSE
        DK=(1.0D0-(1.0D0-DTP)*DARG)/(1.0D0+DTP*DARG)
      ENDIF
      DSM=1.0D0-DTF+DTF*DK
      LFIS=.FALSE.
      DO 30 IBM=1,NBM
      LFIS=LFIS.OR.(CHD(IBM,IFIS,IGR,IDG).NE.0.0)
   30 CONTINUE
      IF(LFIS) THEN
        DO 35 IBM=1,NBM
        IF(IEXP.EQ.0) THEN
          CHEXP(IBM)=CHD(IBM,IFIS,IGR,IDG)
        ELSE
*         exponential transformation
          CHEXP(IBM)=CHD(IBM,IFIS,IGR,IDG)*EXP(-OMEGA(IBM,IGR)*DT)
        ENDIF
   35   CONTINUE
        IF(CMOD.EQ.'BIVAC')THEN
          CALL KINBLM(IPTRK,NBM,LL4,CHEXP,PC(1,IDG,IFIS),WORK1)
        ELSEIF(CMOD.EQ.'TRIVAC')THEN
          CALL KINTLM(IPTRK,NBM,LL4,CHEXP,PC(1,IDG,IFIS),WORK1)
        ENDIF
        DO 40 IND=1,LL4
        SRC(IND,IGR)=SRC(IND,IGR)+PDC(IDG)*DSM*WORK1(IND)
   40   CONTINUE
      ENDIF
   50 CONTINUE
*----
*  FISSION CONTRIBUTION
*----
      DO 92 JGR=1,NGR
      DO 60 IBM=1,NBM
      WORK2(IBM)=CHI(IBM,IFIS,IGR)*SGF(IBM,IFIS,JGR)
   60 CONTINUE
      IF(CMOD.EQ.'BIVAC')THEN
        CALL KINBLM(IPTRK,NBM,LL4,WORK2,EVECT(1,JGR),WORK1)
      ELSEIF(CMOD.EQ.'TRIVAC')THEN
        CALL KINTLM(IPTRK,NBM,LL4,WORK2,EVECT(1,JGR),WORK1)
      ENDIF
      DO 70 IND=1,LL4
      SRC(IND,IGR)=SRC(IND,IGR)+(1.0D0-DTF)*WORK1(IND)
   70 CONTINUE
      DO 91 IDG=1,NDG
      DARG=PDC(IDG)*DT
      IF(IPR.EQ.1)THEN
        DK=0.0D0
      ELSEIF(IPR.EQ.4)THEN
        DK=(1.0D0-DEXP(-DARG))/DARG-DEXP(-DARG)
      ELSE
        DK=(1.0D0-DTP)*DARG/(1.0D0+DTP*DARG)
      ENDIF
      DSM=1.0D0-DTF-DTF*DK
      DO 80 IBM=1,NBM
      WORK2(IBM)=CHD(IBM,IFIS,IGR,IDG)*SGD(IBM,IFIS,JGR,IDG)
   80 CONTINUE
      IF(CMOD.EQ.'BIVAC')THEN
        CALL KINBLM(IPTRK,NBM,LL4,WORK2,EVECT(1,JGR),WORK1)
      ELSEIF(CMOD.EQ.'TRIVAC')THEN
        CALL KINTLM(IPTRK,NBM,LL4,WORK2,EVECT(1,JGR),WORK1)
      ENDIF
      DO 90 IND=1,LL4
      SRC(IND,IGR)=SRC(IND,IGR)-DSM*WORK1(IND)
   90 CONTINUE
   91 CONTINUE
   92 CONTINUE
   93 CONTINUE
*----
*  1/V CONTRIBUTION
*----
      IF(CMOD.EQ.'BIVAC')THEN
        CALL KINBLM(IPTRK,NBM,LL4,OVR(1,IGR),EVECT(1,IGR),WORK1)
      ELSEIF(CMOD.EQ.'TRIVAC')THEN
        CALL KINTLM(IPTRK,NBM,LL4,OVR(1,IGR),EVECT(1,IGR),WORK1)
      ENDIF
      DO 100 IND=1,LL4
      SRC(IND,IGR)=SRC(IND,IGR)+WORK1(IND)
  100 CONTINUE
  110 CONTINUE
*----
*  EDITION
*----
      IF(IMPX.GT.5) THEN
        WRITE(IOS,1002)
        DO 120 IGR=1,NGR
        WRITE(IOS,1003) IGR,(SRC(IND,IGR),IND=1,LL4)
  120   CONTINUE
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(CHEXP,WORK2,WORK1)
      RETURN
*
 1001 FORMAT(/1X,'COMPUTING THE SPACE-TIME KINETICS SOURCE VECTOR',
     1 1X,'ACCORDING TO THE TRACKING TYPE: ',A6/)
 1002 FORMAT(/1X,'=> COMPUTED SPACE-TIME KINETICS SOURCE VECTOR')
 1003 FORMAT(/15H NEUTRON GROUP=,I5/(1P,8D14.5))
      END
