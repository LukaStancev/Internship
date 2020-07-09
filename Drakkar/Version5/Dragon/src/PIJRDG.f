*DECK PIJRDG
      SUBROUTINE PIJRDG(IPRT,NREG,NSOUT,SIGTAL,PROB)
C
C-------------------------    PIJRDG    -------------------------------
C
C 1- SUBROUTINE STATISTICS:
C     NAME     : PIJRDG
C     LEVEL    : 2 (CALLED BY 'EXCELP')
C     USE      : DIAGONAL NORMALIZATION OF COLLISION PROBS (CP)
C     MODIFIED : 91-07-12 (R.R.)
C     MODIFIED : 91-05-17 (G.M.)
C     AUTHOR   : R. ROY (89-06-01)
C     REFERENCE: 'NORMALIZATION TECHNIQUES FOR CP MATRICES',
C                 R.ROY AND G.MARLEAU,
C                 CONF/PHYSOR-90, MARSEILLE/FRANCE, V 2, P IX-40 (1990).
C
C 2- PARAMETERS:
C  INPUT
C     IPRT    : PRINT LEVEL                          I
C     NREG    : # OF ZONES FOR GEOMETRY.             I
C     NSOUT   : # OF SURFACES FOR GEOMETRY.          I
C     SIGTAL  : ALBEDO-SIGT VECTOR                   R(-NSOUT:NREG)
C     PROB    : -CP- MATRIX FOR ALL TYPES.           D(NPRB)
C               NPRB=(NSOUT+NREG+1)*(NSOUT+NREG+2)/2
C  OUTPUT
C     PROB    : -CP- MATRIX FOR ALL TYPES.           D(NPRB)
C
C-------------------------    PIJRDG    -------------------------------
C
      IMPLICIT   NONE
      INTEGER    IPRT,NREG,NSOUT,IR,JR,IPRB,IPRF,IUNK,JUNK,IVOL
      REAL       SIGTAL(-NSOUT:NREG)
      INTEGER    IUNOUT,IPRINT
      PARAMETER (IUNOUT=6, IPRINT=4)
      DOUBLE PRECISION PROB(*),BILAN
      IPRB= 0
      IUNK= 0
      IVOL= NSOUT*(NSOUT+1)/2
C
C     RENORMALIZE ALL DIAGONAL ELEMENTS OF MATRIX *PROB*
      IF( IPRT.GE.IPRINT )THEN
         WRITE(IUNOUT,9000) 'PIJRDG'
      ENDIF
      DO 100 IR = -NSOUT, NREG
         IUNK= IUNK+1
         BILAN=0.0
         DO 10 JR= -NSOUT, IR-1
            IPRB= IPRB+1
            IF( JR.LT.0.OR.SIGTAL(JR).GT.0.0 )THEN
               BILAN=BILAN + PROB(IPRB)
            ENDIF
   10    CONTINUE
         IPRB= IPRB+1
         IPRF= IPRB
         JUNK= IUNK
         DO 20 JR=  IR+1 , NREG
            IPRF= IPRF+JUNK
            JUNK= JUNK+1
            IF( JR.LT.0.OR.SIGTAL(JR).GT.0.0 )THEN
               BILAN=BILAN + PROB(IPRF)
            ENDIF
   20    CONTINUE
         IF( IR.LT.0 )THEN
            IVOL= IVOL+1
            PROB(IPRB)= PROB(IVOL)-BILAN
            IF( IPRT.GE.IPRINT )THEN
              WRITE(IUNOUT,9001) -IR,BILAN
            ENDIF
         ELSEIF( IR.GT.0 )THEN
            IVOL= IVOL+IUNK-1
            IF( SIGTAL(IR).GT.0.0 )THEN
C
C              VOIDS ARE NOT BE RENORMALIZED
               PROB(IPRB)= PROB(IVOL)-BILAN
              IF( IPRT.GE.IPRINT )THEN
                WRITE(IUNOUT,9002) IR,BILAN
              ENDIF
            ENDIF
         ELSE
            IVOL= IVOL+1
         ENDIF
 100  CONTINUE
C
9000  FORMAT('Diagonal correction factors for CP in ',A6)
9001  FORMAT('Surface ',I10,5X,E15.6)
9002  FORMAT('Region  ',I10,5X,E15.6)
      RETURN
      END
