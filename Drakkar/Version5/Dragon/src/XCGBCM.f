*DECK XCGBCM
      SUBROUTINE XCGBCM(IPTRK,NSOUT,NCODE,MATRT)
C
C----------------------------------------------------------------------
C
C 1-  SUBROUTINE STATISTICS:
C
C          NAME      -> XCGBCM
C          USE       -> BUILT BOUNDARY CONDITION MATRIX FOR
C                       REFLECTION AND TRANSMISSION
C          DATE      -> 16-02-1998
C          AUTHOR    -> G. MARLEAU
C
C 2-  PARAMETERS:
C
C INPUT
C  IPTRK   : POINTER TO THE TRACKING FILE                I
C  NSOUT   : NUMBER OF OUTER SURFACE                     I
C  NCODE   : ALBEDO TYPE                                 I(6)
C OUTPUT
C  MATRT   : BC MATRIX FOR REFLECTION/TRANSMISSION       I(NSOUT)
C
C----------------------------------------------------------------------
C
      USE GANLIB
      PARAMETER (NMCOD=6)
      TYPE(C_PTR) IPTRK
      INTEGER NSOUT,NCODE(NMCOD),MATRT(NSOUT),ISOUT
C----
C  INITIALIZE MATRT TO REFLECTION
C----
      DO 100 ISOUT=1,NSOUT
        MATRT(ISOUT)=ISOUT
 100  CONTINUE
C----
C  FOR CARTESIAN CELL LOOK AT PERIODIC BOUNDARY CONDITIONS
C  AND SET TRANSMISSION MATRIX
C----
      IF(NSOUT.EQ.4) THEN
        IF((NCODE(1) .EQ. 4) .AND. (NCODE(2) .EQ.4)) THEN
          MATRT(1)=3
          MATRT(3)=1
        ENDIF
        IF((NCODE(3) .EQ. 4) .AND. (NCODE(4) .EQ.4)) THEN
          MATRT(2)=4
          MATRT(4)=2
        ENDIF
      ENDIF
      CALL LCMPUT(IPTRK,'BC-REFL+TRAN',NSOUT,1,MATRT)
      RETURN
      END
