*DECK FMT
      SUBROUTINE FMT(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To store and retreive information from binary and ascii files.
*
*Copyright:
* Copyright (C) 2009 Ecole Polytechnique de Montreal.
*
*Author(s):
* G. Marleau
*
*Update(s):
* None.
*
*Reference:
*
*Parameters: input
* NENTRY  number of data structures transfered to this module.
* HENTRY  name of the data structures.
* IENTRY  data structure type where:
*         \begin{itemize}
*         \item IENTRY=1 for LCM memory object;
*         \item IENTRY=2 for XSM file;
*         \item IENTRY=3 for sequential binary file;
*         \item IENTRY=4 for sequential ASCII file.
*         \end{itemize}
*JENTRY   access permission for the data structure where:
*         \begin{itemize}
*         \item JENTRY=0 for a data structure in creation mode;
*         \item JENTRY=1 for a data structure in modifications mode;
*         \item JENTRY=2 for a data structure in read-only mode.
*         \end{itemize}
* KENTRY  data structure pointer.
*
*Comments:
* Instructions for the use of the FMT: module:
*   [[ OutFiles ]] := FMT: [[ InDds ]] :: (FMTget) ;
* or
*   [[ UpdDds ]] := FMT: [[ UpdDds ]] [[ Infiles ]] :: (FMTget) ;
* where
*     OutFiles : sequential binary/ASCII output files.
*     UpdDds   : Data structures to update.
*     InFiles  : sequential binary/ASCII input files.
*     InDds    : Input data structure.
*     (FMTget) : Processing options
*                (read from input using the FMTGET routine).
*
*-----------------------------------------------------------------------
*
      USE              GANLIB
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          NENTRY,IENTRY(NENTRY),JENTRY(NENTRY)
      TYPE(C_PTR)      KENTRY(NENTRY)
      CHARACTER        HENTRY(NENTRY)*12
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='FMT   ')
      INTEGER          ILCMUP,ILCMDN,MXFIL,MXOPT
      PARAMETER       (ILCMUP=1,ILCMDN=2,MXFIL=20,MXOPT=20)
      INTEGER          NSTATE
      PARAMETER       (NSTATE=40)
*----
*  Local variables
*----
      CHARACTER*12     SENTRY(MXFIL)
      INTEGER          IEN
      CHARACTER        HSIGN*12
      INTEGER          IPRINT,NOPT,IOPT(MXOPT)
*----
*  Validate entry parameters
*----
      IF(NENTRY .GT. MXFIL) CALL XABORT(NAMSBR//
     >  ': Too many files or data structures for this module.')
*----
*  Scan data structure to determine signature (input or update)
*----
      DO IEN=1,NENTRY
        SENTRY(IEN)='            '
        IF(IENTRY(IEN) .EQ. 1 .OR. IENTRY(IEN) .EQ. 2) THEN
          IF(JENTRY(IEN) .NE. 0) THEN
            CALL LCMGTC(KENTRY(IEN),'SIGNATURE',12,1,HSIGN)
            SENTRY(IEN)=HSIGN
          ENDIF
        ENDIF
      ENDDO
*----
*  Recover processing option
*----
      NOPT=MXOPT
      CALL XDISET(IOPT,NOPT,0)
      CALL FMTGET(IPRINT,NOPT,IOPT)
*----
*  Process files
*----
      IF(IOPT(1) .EQ. 1) THEN
*----
*  SUS3D format
*----
        CALL FMTSUS(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY,SENTRY,
     >              IPRINT,NOPT,IOPT)
      ELSE IF(IOPT(1) .EQ. 2) THEN
*----
*  DIRFLX format
*----
        CALL FMTDFL(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY,SENTRY,
     >              IPRINT)
      ELSE IF(IOPT(1) .EQ. 3) THEN
*----
*  BURNUP format
*----
        CALL FMTBRN(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY,SENTRY,
     >              IPRINT,NOPT,IOPT)
      ENDIF
*----
*  Processing finished, return
*----
      RETURN
*----
*  Warning formats
*----
      END
