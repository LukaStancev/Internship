*DECK MSTCDI
      SUBROUTINE MSTCDI(IPSTR,ACSTR,IPRINT,NBLOCK,MYDIR,BLNAM,BLTYP,
     1                  BLLEN)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Create a new directory and move in a structure according to a defined
* directory name.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Le Tellier
*
*Parameters: input
* ACSTR   structure access.
* IPRINT  level of print index.
* NBLOCK  number of existing block in the directory.
* MYDIR   name of the directory to be created/moved in.
* BLNAM   names of these blocks.
* BLTYP   types of these blocks.
* BLLEN   lengths of these blocks.
*
*Parameters: input/output
* IPSTR   entering/leaving directory address.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPSTR
      INTEGER ACSTR,IPRINT,NBLOCK,BLNAM(3,NBLOCK+1),BLTYP(NBLOCK+1),
     1 BLLEN(NBLOCK+1)
      CHARACTER*12 MYDIR
*----
*  LOCAL VARIABLES
*----
      INTEGER IOUT,NSTATE
      PARAMETER (IOUT=6,NSTATE=40)
      INTEGER ISTATE(NSTATE),KCHAR(3),JJ,ILEN,ITYP
      LOGICAL EXIST
*----      
* PERFORM CD RELATED ACTIONS
*----
      IF (MYDIR.EQ.'..') THEN
*     GOING TO FATHER DIR
         IF (IPRINT.GT.2)
     1        WRITE(IOUT,*) 'MSTCDI: GOING TO FATHER DIR'
         CALL LCMSIX(IPSTR,' ',2)
      ELSE
*     GOING TO SON DIR
         EXIST=.FALSE.
         IF (NBLOCK.NE.0) THEN
*        IS THIS SON DIR ALREADY PART OF THE STRUCTURE ?
            READ(MYDIR,'(3A4)') (KCHAR(JJ),JJ=1,3)
            DO JJ=1,NBLOCK
               IF((BLNAM(1,JJ).EQ.KCHAR(1)).AND.
     1            (BLNAM(2,JJ).EQ.KCHAR(2)).AND.
     2            (BLNAM(3,JJ).EQ.KCHAR(3))) THEN
                  IF (BLLEN(JJ).NE.-1)
     1          CALL XABORT('MSTCDI: '//MYDIR//' IS AN EXISTING BLOCK.')
                  EXIST=.TRUE.
                  GOTO 10
               ENDIF
            ENDDO
 10         CONTINUE
         ENDIF
         IF (EXIST) THEN
*        YES: 
            IF (IPRINT.GT.2)
     1         WRITE(IOUT,*) 'MSTCDI: ENTERING EXISTING DIR '//MYDIR
         ELSE
*        NO:
            CALL LCMLEN(IPSTR,MYDIR,ILEN,ITYP)
            IF (ILEN.NE.0) THEN
*              IT IS ASSUMED THAT THIS IS AN EXTERNAL STRUCTURE FROM WHICH INFORMATION CAN BE RETRIEVED
               EXIST=.TRUE.
               GOTO 20
            ENDIF
            IF (ACSTR.EQ.2)
     1        CALL XABORT('MSTCDI: CANNOT CREATE DIR IN READ-ONLY MODE')
            IF (IPRINT.GT.2)
     1         WRITE(IOUT,*) 'MSTCDI: CREATING DIR '//MYDIR
            CALL LCMGET(IPSTR,'STATE-VECTOR',ISTATE)
            ISTATE(40)=ISTATE(40)+1
            READ(MYDIR,'(3A4)') (BLNAM(JJ,NBLOCK+1),JJ=1,3)
            BLTYP(NBLOCK+1)=0
            BLLEN(NBLOCK+1)=-1
            CALL LCMPUT(IPSTR,'STATE-VECTOR',NSTATE,1,ISTATE)
            CALL LCMPUT(IPSTR,'REC-NAMES',3*(NBLOCK+1),3,BLNAM)
            CALL LCMPUT(IPSTR,'REC-TYPES',(NBLOCK+1),1,BLTYP)
            CALL LCMPUT(IPSTR,'REC-LENGTHS',(NBLOCK+1),1,BLLEN)
         ENDIF
 20      CALL LCMSIX(IPSTR,MYDIR,1)
         IF (.NOT.EXIST) THEN
            ISTATE(:NSTATE)=0
            CALL LCMPUT(IPSTR,'STATE-VECTOR',NSTATE,1,ISTATE)
         ENDIF
      ENDIF
      
      RETURN
      END
