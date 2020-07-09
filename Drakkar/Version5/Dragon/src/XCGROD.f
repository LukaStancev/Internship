*DECK XCGROD
      SUBROUTINE XCGROD(NRT,MSROD,NRODS,RODS,MATROD,RODR)
C
C----------------------------------------------------------------------
C
C 1-  SUBROUTINE STATISTICS:
C
C          NAME      -> XCGROD
C          USE       -> CHECK GEOMETRY AND REORDER ROD CLUSTERS IF
C                       NECESSARY
C          DATE      -> 13-01-1994
C          AUTHOR    -> G. MARLEAU
C
C 2-  PARAMETERS:
C
C INPUT
C  NRT     : NUMBER OF ROD TYPES                         I
C  MSROD   : MAXIMUM NUMBER OF SUBRODS PER RODS          I
C INPUT/OUTPUT
C  NRODS   : INTEGER DESCRIPTION OF ROD OF A GIVEN TYPE  I(3,NRT)
C            NRODS(1,IRT) = NUMBER OF ROD
C            NRODS(2,IRT) = NUMBER OF SUBRODS IN ROD
C            NRODS(3,IRT) = FIRST CONCENTRIC REGION
C  RODS    : REAL DESCRIPTION OF ROD OF A GIVEN TYPE     R(2,NRT)
C            RODS(1,IRT) = ROD CENTER RADIUS
C            RODS(2,IRT) = ANGULAR POSITION OF FIRST ROD
C  MATROD  : TYPE OF MATERIAL FOR EACH SUBROD            I(MSROD,NRT)
C  RODR    : SUBROD RADIUS                               R(MSROD,NRT)
C
C----------------------------------------------------------------------
C
      INTEGER    IOUT
      REAL       PI
      PARAMETER (IOUT=6,PI=3.1415926535898)
      INTEGER    NRT,NRODS(3,NRT),MATROD(MSROD,NRT)
      REAL       RODS(2,NRT),RODR(MSROD,NRT)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IORD
C----
C  SCRATCH STORAGE ALLOCATION
C   IORD    : NEW ROD CLUSTER ORDER                       I(NRT)
C----
      ALLOCATE(IORD(NRT))
C----
C  CLASSIFY ROD CLUSTER BY INCREASING DISTANCE OF CENTER AND ANGLE
C----
      DO 100 IRT=1,NRT
        IORD(IRT)=IRT
 100  CONTINUE
      DO 110 IRT=2,NRT
        REFR=RODS(1,IRT)
        REFA=RODS(2,IRT)
        IPOS=IORD(IRT)
        DO 111 JRT=IRT-1,1,-1
          KRT=JRT
          IF(RODS(1,JRT).GT.REFR) THEN
            RODS(1,JRT+1)=RODS(1,JRT)
            RODS(2,JRT+1)=RODS(2,JRT)
            IORD(JRT+1)=IORD(JRT)
          ELSE IF(RODS(1,JRT).EQ.REFR) THEN
            IPOS=-IPOS
            GO TO 112
          ELSE
            GO TO 112
          ENDIF
 111    CONTINUE
        KRT=0
 112    CONTINUE
        RODS(1,KRT+1)=REFR
        RODS(2,KRT+1)=REFA
        IORD(KRT+1)=IPOS
        IF(IPOS.LT.0) THEN
          DO 113 JRT=KRT,1,-1
            LRT=JRT
            IF((RODS(2,JRT).GT.REFA).AND.
     >         (RODS(1,JRT).EQ.REFR)) THEN
              RODS(1,JRT+1)=RODS(1,JRT)
              RODS(2,JRT+1)=RODS(2,JRT)
              IORD(JRT+1)=IORD(JRT)
            ELSE
              GO TO 114
            ENDIF
 113      CONTINUE
          LRT=0
 114      CONTINUE
          RODS(1,LRT+1)=REFR
          RODS(2,LRT+1)=REFA
          IORD(LRT+1)=-IPOS
        ENDIF
 110  CONTINUE
C----
C  REORDER REMAINING VECTORS NRODS,MATROD,RODR
C----
      DO 140 IRT=1,NRT
        JRT=IORD(IRT)
        IF(JRT.NE.IRT) THEN
          DO 141 IX=1,3
            NNR=NRODS(IX,IRT)
            NRODS(IX,IRT)=NRODS(IX,JRT)
            NRODS(IX,JRT)=NNR
 141      CONTINUE
          DO 142 IS=1,MSROD
            MATT=MATROD(IS,IRT)
            MATROD(IS,IRT)=MATROD(IS,JRT)
            MATROD(IS,JRT)=MATT
            RROD=RODR(IS,IRT)
            RODR(IS,IRT)=RODR(IS,JRT)
            RODR(IS,JRT)=RROD
 142      CONTINUE
          DO 143 KRT=IRT+1,NRT
            IF(IORD(KRT).EQ.IRT) THEN
              IORD(KRT)=JRT
              IORD(IRT)=IRT
              GO TO 144
            ENDIF
 143      CONTINUE
 144      CONTINUE
        ENDIF
 140  CONTINUE
C----
C  FIND IF ROD OVERLAPP
C----
      DO 150 IRT=1,NRT
        NRDB=NRODS(1,IRT)
        NSBRB=NRODS(2,IRT)
        RODRB=RODR(NSBRB,IRT)
        RODRB2=RODRB*RODRB
        RDPB=RODS(1,IRT)
        XBOT=RDPB-RODRB
        DANGB=2.*PI/FLOAT(NRDB)
        ANGB=RODS(2,IRT)
C----
C  CHECK FOR ROD OVERLAPP INSIDE EACH CLUSTER
C----
        IF(NRDB.GT.1) THEN
          IF(RODRB.GT.RDPB) THEN
            WRITE(IOUT,'(1X,24HROD OVERLAP IN CLUSTER =,I10)') IRT
            CALL XABORT('XCGROD: ROD OVERLAP IN A CLUSTER')
          ELSE
            ANGMIN=2.*ASIN(RODRB/RDPB)
            IF(DANGB.LE.ANGMIN) THEN
              WRITE(IOUT,'(1X,24HROD OVERLAP IN CLUSTER =,I10)') IRT
              CALL XABORT('XCGROD: ROD OVERLAP IN A CLUSTER')
            ENDIF
          ENDIF
        ENDIF
C----
C  CHECK FOR ROD OVERLAPP BETWEEN DIFFERENT CLUSTERS
C----
        DO 151 JRT=IRT-1,1,-1
          NRDT=NRODS(1,JRT)
          NSBRT=NRODS(2,JRT)
          RODRT=RODR(NSBRT,JRT)
          RODRT2=RODRT*RODRT
          RDPT=RODS(1,JRT)
          XTOP=RDPT+RODRT
          DANGT=2.*PI/FLOAT(NRDT)
          ANGT=RODS(2,JRT)
C----
C  NO OVERLAPP
C----
          IF(XTOP.LT.XBOT) GO TO 152
C----
C  SOME OVERLAPP POSSIBLE TEST FOR INTERSECTION
C----
          ANG1=ANGB
          DO 160 IA1=1,NRDB
C----
C  FIND POSITION OF ROD (X0,Y0)
C----
            X01=RDPB*COS(ANG1)
            Y01=RDPB*SIN(ANG1)
            RRX=RODRB2-X01*X01
            RRY=RODRB2-Y01*Y01
            XY=X01*Y01
            RR1=(RRX-Y01*Y01)
            ANG2=ANGT
            DO 161 IA2=1,NRDT
              X02=RDPT*COS(ANG2)
              Y02=RDPT*SIN(ANG2)
              RR2=(RODRT2-X02*X02-Y02*Y02)
C----
C  CHECK FOR ROD INSIDE ROD
C----
              DELX=X02-X01
              DELY=Y02-Y01
              DIST=SQRT(DELX**2+DELY**2)
              IF(DIST.LT.RODRT+RODRB) THEN
                WRITE(IOUT,'(1X,25HROD OVERLAP IN CLUSTERS =,2I10)')
     >                 IRT,JRT
                CALL XABORT('XCGROD: ROD OVERLAP IN 2 CLUSTERS')
              ENDIF
C----
C  FIND IF CIRCLES
C  (X-X01)**2+(Y-Y01)**2=RODRB*2
C  (X-X02)**2+(Y-Y02)**2=RODRT*2
C  INTERSECT
C----
              IF(X02.NE.X01) THEN
                CCR=1./DELX
                BBR=-DELY*CCR
                AAR=0.5*CCR*(RR1-RR2)
                ARGSQ=AAR*(2.*X01-2.*BBR*Y01-AAR)
     >               +BBR*(BBR*RRY+2.*XY)+RRX
              ELSE
                CCR=1./DELY
                BBR=-DELX*CCR
                AAR=0.5*CCR*(RR1-RR2)
                ARGSQ=AAR*(2.*Y01-2.*BBR*X01-AAR)
     >               +BBR*(BBR*RRX+2.*XY)+RRY
              ENDIF
              IF(ARGSQ.GE.0.0) THEN
                WRITE(IOUT,'(1X,25HROD OVERLAP IN CLUSTERS =,2I10)')
     >                   IRT,JRT
                CALL XABORT('XCGROD: ROD OVERLAP IN 2 CLUSTERS')
              ENDIF
              ANG2=ANG2+DANGT
 161        CONTINUE
            ANG1=ANG1+DANGB
 160      CONTINUE
 151    CONTINUE
 152    CONTINUE
 150  CONTINUE
C----
C  SCRATCH STORAGE DEALLOCATION
C----
      DEALLOCATE(IORD)
C----
C  RETURN
C----
      RETURN
      END
