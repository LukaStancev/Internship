*DECK EDIHMX
      SUBROUTINE EDIHMX(IPTRK,NREGIO,NMERGE,IMERGE)
C
C---------------------------  EDIHMX  ---------------------------------
C
C  1- PROGRAMME STATISTICS:
C      NAME     : EDIHMX
C      USE      : FIND MERGE VECTOR IN IPTRK FOR MERGE BY HMIX
C      MODIFIED : 2001/10/30 (G.M)
C      AUTHOR   : G.MARLEAU
C
C  2- ROUTINE PARAMETERS:
C      IPTRK    : CALCULATION TRACKING DATA STRUCTURE
C                 ***> INTEGER IPTRKI
C      NREGIO   : NUMBER OF REGIONS
C                 ***> INTEGER NREGIO
C      NMERGE   : FINAL NUMBER OF MERGED REGIONS
C                 ***> INTEGER NMERGE
C      IMERGE   : MERGED REGIONS POSITION
C                 ***> INTEGER IMERGE(NREGIO)
C
C---------------------------   EDIHMX  --------------------------------
C
      USE         GANLIB
      IMPLICIT    NONE
      INTEGER     IOUT,NSTATE
      CHARACTER   NAMSBR*6
      PARAMETER  (IOUT=6,NSTATE=40,NAMSBR='EDIHMX')
C----
C  ROUTINE PARAMETERS
C----
      TYPE(C_PTR) IPTRK
      INTEGER     NREGIO
      INTEGER     NMERGE
      INTEGER     IMERGE(NREGIO)
C----
C  LOCAL PARAMETERS
C----
      INTEGER     IMRGLN,IMRGTY,IREG
      INTEGER     ISTATE(NSTATE)
C----
C  IMERGE is HOMMATCOD or MATCOD
C----
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      CALL LCMLEN(IPTRK,'HOMMATCOD   ',IMRGLN,IMRGTY)
      IF(IMRGLN .EQ. 0) THEN
        WRITE(IOUT,8000) NAMSBR
        CALL LCMGET(IPTRK,'MATCOD      ',IMERGE)
      ELSE
        CALL LCMGET(IPTRK,'HOMMATCOD   ',IMERGE)
      ENDIF
C----
C  Check for double heterogeneity
C----
      IF(ISTATE(40).EQ.1) THEN
        CALL EDIBHX (NREGIO,IPTRK,IMRGLN,IMERGE)
      ENDIF
      IF(IMRGLN.NE.NREGIO) CALL XABORT('EDIHMX: bad nb of regions')
C----
C  Compute number of merged regions
C----
      NMERGE=0
      DO IREG=1,NREGIO
        NMERGE=MAX(NMERGE,IMERGE(IREG))
      ENDDO
      RETURN
C----
C  WARNING FORMAT
C----
 8000 FORMAT('***** Warning in routine - ',A6,' - *****'/
     >'No HMIX data in GEO: for NXT tracking file'/
     >'Homogenize using MIX instead of HMIX')
      END
