*DECK KDRVER
      SUBROUTINE KDRVER(REV,DATE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To extract CVS or SVN version and production date.
*
*Copyright:
* Copyright (C) 2006 Ecole Polytechnique de Montreal
*
*Author(s): A. Hebert
*
*Parameters: output
* REV     revision character identification
* DATE    revision character date 
*
*-----------------------------------------------------------------------
*
      CHARACTER REV*48,DATE*64
*
      REV='Version 5.0.6 ($Revision: 1800 $)'
      DATE='$Date: 2020-02-01 14:49:32 -0500 (Sat, 01 Feb 2020) $'
      IF(REV(22:).EQ.'ion$)') THEN
*        CVS or SVN keyword expansion not performed
         REV='Version 5.0.6'
         DATE='February 1, 2020'
      ENDIF
      RETURN
      END
