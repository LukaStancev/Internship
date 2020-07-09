*DECK NXTHRS
      FUNCTION NXTHRS(ITRCUR,ISYM)
*
*----------
*
*Purpose:
* Find new DRAGON TURN factor after an hexagonal symmetry
* is applied on an old DRAGON TURN factor.
*
*Copyright:
* Copyright (C) 2004 Ecole Polytechnique de Montreal.
*
*Author(s): G. Marleau.
*
*Update(s):
* 2005/10/290: Validated with TSTTRS.
*
*Parameters: input
* ITRCUR  initial turn factor.
* ISYMM   symmetry to consider where
*         \begin{itemize}
*         \item \verb|ISYMM|=-1 indicates $Z$ reflection symmetry;
*         \item \verb|ISYMM|=i indicates $H$ reflection symmetry;
*         \end{itemize}
*
*Parameters: output
* NXTHRS  turn factor after symmetry is applied.
*
*----
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          ITRCUR,ISYM
      INTEGER          NXTHRS
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTHRS')
      INTEGER          MAXTUR
      PARAMETER       (MAXTUR=12)
*----
*  Local variables
*----
      INTEGER          ISZ,ICTP,IRTR
*----
*  Test input data
*  ITRCUR can be 1-8 or 13-20
*  ISYM can be 1-4
*----
      IF((ITRCUR .LE.  0) .OR.
     >   (ITRCUR .GE.  9 .AND. ITRCUR .LE. 12) .OR.
     >   (ITRCUR .GE. 21)) CALL XABORT(NAMSBR//
     >  ': Invalid TURN')
      IF(ISYM .LT. -1 .OR.
     >   ISYM .EQ.  0 .OR.
     >   ISYM .GT. 4 ) CALL XABORT(NAMSBR//
     >  ': Invalid symmetry')
*----
*  Find current symmetry factor in Z (ISZ) and current turn
*  number in plane X-Y (ICTP)
*----
      ISZ=((ITRCUR-1)/MAXTUR)
      ICTP=MOD(ITRCUR-1,MAXTUR)+1
      IF(ISYM .EQ. -1) THEN
*----
*  Z symmetry
*----
        ISZ=(1-ISZ)
      ELSE
*----
*  X, X-Y AND Y symmetry
*----
        IRTR=((ICTP-1)/4)*4
        ICTP=MOD(4-ICTP+IRTR+ISYM,4)+5-IRTR
      ENDIF
      NXTHRS=ICTP+MAXTUR*ISZ
      RETURN
      END
