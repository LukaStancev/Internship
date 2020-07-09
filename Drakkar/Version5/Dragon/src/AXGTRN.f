      FUNCTION           AXGTRN(ITRCUR)
*----
*  Associate to TURN number a DRAGON name
*----
      IMPLICIT           NONE
*----
*  Local parameters
*----
      INTEGER            MAXTUR
      CHARACTER          NAMSBR*6
      PARAMETER         (MAXTUR=12,NAMSBR='AXGTRN')
*----
*  Routine input and output variables 
*----     
      INTEGER            ITRCUR
      CHARACTER          AXGTRN*(*)
*----
*  local variables
*----
      CHARACTER*2        CTURN(2*MAXTUR)
      SAVE               CTURN
*----
*  DEFINITION OF TURNS
*----
      DATA CTURN        /' A',' B',' C',' D',' E',' F',' G',' H',
     >                   ' I',' J',' K',' L',
     >                   '-A','-B','-C','-D','-E','-F','-G','-H',
     >                   '-I','-J','-K','-L'/
      IF(ITRCUR .LE. 0 .OR. ITRCUR .GT. 2*MAXTUR) CALL XABORT(NAMSBR//
     >  ': INVALID TURN NUMBER')
      AXGTRN=CTURN(ITRCUR)
      RETURN
      END
