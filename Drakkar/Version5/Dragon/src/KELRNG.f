*DECK KELRNG
      FUNCTION KELRNG(  IPRT,   NDIM, NEXTGE,   NCPC,  MINDO,  MAXDO,
     >                ICORDO,  NSURO,  NVOLO, IDLGEO,
     >                  MAXC, RMESHO, MATGEO,  VOLSO, INDEXO )
************************************************************************
*                                                                      *
*           NAME: KELRNG                                               *
*      COMPONENT: EXCELL                                               *
*        VERSION: 1.0                                                  *
*       CREATION: 90/01                                                *
*       MODIFIED: 00/03 (R.R.) DECLARE ALL VARIABLE TYPES              *
*         AUTHOR: ROBERT ROY                                           *
*                                                                      *
*       FUNCTION: THIS FUNCTION WILL RENUMBER ALL ZONES AND SURFACES   *
*                 FOR A BLOCK BY THE COORDINATE (RECT/CYL) VALUES.     *
*                                                                      *
*          LEVEL: 4 (CALLED BY 'XELTRP')                               *
*                                                                      *
*--------+-------------- V A R I A B L E S -------------+--+-----------*
*  NAME  /                  DESCRIPTION                 /IO/MOD(DIMENS)*
*--------+----------------------------------------------+--+-----------*
* KELRNG / FUNCTION = # OF SURFACES & ZONES RENUMBERED. /.O/INT        *
* IPRT   / INTERMEDIATE PRINTING LEVEL FOR OUTPUT.      /I./INT        *
* NDIM   / # OF DIMENSIONS                              |I./INT        *
* NEXTGE / RECTANGULAR(0)/CIRCULAR(1) BOUNDARY.         |I./INT        *
* NCPC   / # OF CYLINDRES IN A TYPE + 3    (.LT.NC3MAX).|I./INT        *
* MINDO  / MIN INDEX VALUES FOR ALL AXES (RECT/CYL).    |I./INT(NC3MAX)*
* MAXDO  / MAX INDEX VALUES FOR ALL AXES (RECT/CYL).    |I./INT(NC3MAX)*
* ICORDO / PRINCIPAL AXIS DIRECTIONS (X/Y/Z) MESHES.    /../INT(NC3MAX)*
* NSURO  / # OF SURFACES FOR A SPECIFIC GEOMETRY.       /I./INT        *
* NVOLO  / # OF ZONES FOR A SPECIFIC GEOMETRY.          /I./INT        *
* IDLGEO / SPECIFIC POSITION FOR A GEOMETRY.            /I./INT        *
* MAXC   / DIMENSION OF RMESHO().                       /I./INT        *
* RMESHO / REAL MESH VALUES (RECT/CYL).                 /I./REL(MAXC  )*
* MATGEO / MATERIAL #S CORRESPONDING TO GEOMETRIES.     /I./INT(NGIDL) *
* VOLSO  / VOLUMES & SURFACES FOR EACH GEOMETRY.        /I./REL(NGIDL) *
* INDEXO / COORDINATES FOR ZONES & SURFACES OF A CELL.  /.O/INT(4,*)   *
*--------+--------------- R O U T I N E S --------------+--+-----------*
*  NAME  /                  DESCRIPTION                                *
*--------+-------------------------------------------------------------*
* LELCRN / TO DECIDE CROWN/RECTANGLE INTERSECTIONS.                    *
************************************************************************
C
      IMPLICIT        NONE
C
      INTEGER         KELRNG, IPRT, NDIM, NEXTGE, NCPC, MAXC,
     >                NSURO, NVOLO, IDLGEO
      INTEGER         ICUR(4), MINDO(NCPC), MAXDO(NCPC),
     >                ICORDO(NCPC), INDEXO(4,*), MATGEO(*), IXYZ(3),
     >                MINT(3), MAXT(3), INCT(3), JINT(3), JAXT(3)
      REAL            VOLSO(*), RMESHO(MAXC)
      DOUBLE PRECISION  RECT(2,3), RAC(2), CEC(2), RAYMAX, RAYXY
      LOGICAL         LELCRN
C
      INTEGER         NSU, NVO, IDLGE, NCP, I, J, KSUR, ISUX, IVOX,
     >                INOX, ICX, ICY, ICZ, IDX, IDY, IDZ, JX, JY, JZ,
     >                IMAT, IMATX, IMATY, IMATZ, IMATYZ, IMATR, NBEXT,
     >                JCP, IX, IY, JRAY, NO, NESUR, NEVOL
C
      INTEGER         IOUT, IND
      PARAMETER     ( IOUT=6 )
C
      IND(I)= IDLGE + I
C
      NSU   =  NSURO
      NVO   =  NVOLO
      IDLGE = IDLGEO
      NCP   =  NCPC 
C
C     INITIALISATION OF INDEX AND VARIOUS THINGS
      DO 20 I= NSU, NVO
         MATGEO(IND(I))=0
         VOLSO(IND(I))=0.0
         DO 10 J= 1, 4
            INDEXO(J,IND(I))= 0
   10    CONTINUE
   20 CONTINUE
      KSUR= MOD(NDIM+1,3)
      DO 25 I= 1, 3
         IXYZ(I)= ABS(ICORDO(I))
         JINT(I)= MINDO(IXYZ(I))
         JAXT(I)= MAXDO(IXYZ(I))+1
         IF( ICORDO(I).GT.0 )THEN
            IF( I.EQ.3 )THEN
               MINT(I)= MINDO(IXYZ(I))+1-KSUR
               MAXT(I)= MAXDO(IXYZ(I))+KSUR
            ELSE
               MINT(I)= MINDO(IXYZ(I))
               MAXT(I)= MAXDO(IXYZ(I))+1
            ENDIF
            INCT(I)= +1
         ELSE
            IF( I.EQ.3 )THEN
               MINT(I)= MAXDO(IXYZ(I))+KSUR
               MAXT(I)= MINDO(IXYZ(I))+1-KSUR
            ELSE
               MINT(I)= MAXDO(IXYZ(I))+1
               MAXT(I)= MINDO(IXYZ(I))
            ENDIF
            INCT(I)= -1
         ENDIF
   25 CONTINUE
C
      KELRNG= 0
      ISUX= 0
      IVOX= 0
      INOX= 0
C
C     NUMBER ZONES & SURFACES
      IF( NCP.LT.4 )THEN
C        THERE ARE NO CYLINDER AT ALL
         J= 3
         ICUR(4)= 0
         ICZ= 3
      ELSE
         J= 4
         CEC(1)= DBLE(RMESHO(MINDO(J)-2))
         CEC(2)= DBLE(RMESHO(MINDO(J)-1))
         ICZ= ICORDO(J)
      ENDIF
C
C     AXIS ORDER IN TRUE GEOMETRY
      ICX= MOD(ICZ  , 3) + 1
      ICY= MOD(ICZ+1, 3) + 1
C
C     AXIS ORDER FOR NUMBERING PROCESS
      IDX= IXYZ(ICX)
      IDY= IXYZ(ICY)
      IDZ= IXYZ(ICZ)
C
C     LOOP OVER ALL "ICZ,ICY,ICX" ZONES, THEN RADIUS
      DO 260 JZ= MINT(ICZ), MAXT(ICZ), INCT(ICZ)
         ICUR(IDZ)= JZ-1
         IF( JZ.NE.JINT(ICZ).AND.JZ.NE.JAXT(ICZ) )THEN
            IMATZ= 0
         ELSE
            IMATZ=   - 2*ICZ
            IF(    (INCT(ICZ).EQ.+1.AND.JZ.EQ.MINT(ICZ))
     >         .OR.(INCT(ICZ).EQ.-1.AND.JZ.EQ.MAXT(ICZ)) )
     >              IMATZ= IMATZ+1
         ENDIF
         DO 250 JY= MINT(ICY), MAXT(ICY), INCT(ICY)
            RECT(1,IDY)= DBLE(RMESHO(MAX(JINT(ICY)  ,JY-1)))
            RECT(2,IDY)= DBLE(RMESHO(MIN(JAXT(ICY)-1,JY  )))
            ICUR(IDY)= JY-1
            IF( JY.NE.JINT(ICY).AND.JY.NE.JAXT(ICY) )THEN
               IMATY= 0
            ELSE
               IMATY= -2*IDY
               IF(    (INCT(ICY).EQ.+1.AND.JY.EQ.MINT(ICY))
     >            .OR.(INCT(ICY).EQ.-1.AND.JY.EQ.MAXT(ICY)) )
     >                 IMATY= IMATY+1
            ENDIF
C
C           TO EXCLUDE LINES
            IF( IMATY*IMATZ .NE. 0 ) GO TO 250
            IMATYZ= IMATY + IMATZ
            DO 240 JX= MINT(ICX), MAXT(ICX), INCT(ICX)
               RECT(1,IDX)= DBLE(RMESHO(MAX(JINT(ICX)  ,JX-1)))
               RECT(2,IDX)= DBLE(RMESHO(MIN(JAXT(ICX)-1,JX  )))
               ICUR(IDX)= JX-1
               IF( JX.NE.JINT(ICX).AND.JX.NE.JAXT(ICX) )THEN
                  IMATX= 0
               ELSE
                  IMATX= -2*IDX
                  IF(    (INCT(ICX).EQ.+1.AND.JX.EQ.MINT(ICX))
     >               .OR.(INCT(ICX).EQ.-1.AND.JX.EQ.MAXT(ICX)) )
     >                    IMATX= IMATX+1
               ENDIF
C
C              TO EXCLUDE SINGLE POINTS
               IF( IMATYZ*IMATX .NE. 0 ) GO TO 240
               IMAT= IMATYZ + IMATX
               NBEXT=1
               IF( NCP.GT.3 )THEN
                  IMATR= IMAT
                  RAC(1)= 0.0D0
                  DO 230 JRAY= MINDO(J), MAXDO(J)
                     RAC(2)= DBLE(RMESHO(JRAY))
                     ICUR(4)= JRAY-1
                     IF(LELCRN(CEC,RAC,RECT(1,ICX),RECT(1,ICY)))THEN
                        IF( IMAT.EQ.0 )THEN
C                          ZONE NUMBERING
                           IVOX= IVOX + 1
                           INOX= INOX + 1
                           NO=INOX
                           IMATR= IVOX
                        ELSE
C                          SURFACE NUMBERING
                           ISUX= ISUX - 1
                           NO= ISUX
                        ENDIF
C                       IDENTIFY FACE AND CHARGE THE ZONE OR SURFACE NO
                        MATGEO(IND(NO))= IMATR
                        DO 220 JCP= 1, 4
                           INDEXO(JCP,IND(NO))= ICUR(JCP)
  220                   CONTINUE
                     ELSE
                        IF( IMAT.EQ.0 )THEN
C                          ZONE NUMBERING
                           INOX= INOX + 1
C                       IDENTIFY FACE AND CHARGE THE ZONE OR SURFACE NO
                           MATGEO(IND(INOX))= -1
C                        ELSE
C                          ISUX=ISUX-1
                        ENDIF
                     ENDIF
                     RAC(1)= RAC(2)
  230             CONTINUE
                  ICUR(4)= MAXDO(J)
                  RAYMAX= DBLE(RMESHO(MAXDO(J)))
                  NBEXT=0
                  DO 232 IX= 1, 2
                  DO 232 IY= 1, 2
                  RAYXY= (RECT(IX,ICX)-CEC(1))*(RECT(IX,ICX)-CEC(1))
     >                 + (RECT(IY,ICY)-CEC(2))*(RECT(IY,ICY)-CEC(2))
                    IF( RAYXY.GE.RAYMAX ) NBEXT= NBEXT + 1
  232             CONTINUE
               ENDIF
               IF( NBEXT.EQ.0 )THEN
C
C                 NUMBER 'INSIDE' OF CYLINDER
                  IF( NEXTGE.EQ.0 )THEN
C
C                    CONSIDER ONLY FOR OVERALL CARTESIAN GEOMETRY
C                    SET IMAT TO -1 TO IDENTIFY REGION EXTRACTED
                     IF( IMAT.EQ.0 )THEN
C
C                       ZONE NUMBERING
                        INOX= INOX + 1
                        MATGEO(IND(INOX))= -1
                     ENDIF
                  ENDIF
               ELSE
C
C                 NUMBER 'OUTSIDE' OF CYLINDER
                  IF( IMAT.EQ.0 )THEN
C
C                    ZONE NUMBERING
                     IVOX= IVOX + 1
                     INOX= INOX + 1
                     IMAT= IVOX
                     NO  = INOX
                  ELSE
C
C                    SURFACE NUMBERING
                     ISUX= ISUX - 1
                     NO= ISUX
                  ENDIF
C
C                 IDENTIFY FACE AND CHARGE THE ZONE OR SURFACE NO
                  MATGEO(IND(NO))= IMAT
                  DO 235 JCP= 1, 4
                     INDEXO(JCP,IND(NO))= ICUR(JCP)
  235             CONTINUE
               ENDIF
  240       CONTINUE
  250    CONTINUE
  260 CONTINUE
C
      KELRNG= IVOX - ISUX + 1
C
      IF( IPRT.GT.5 )THEN
         NESUR= 0
         NEVOL= 0
         DO 549 I= NSU, NVO
            IF( I.GT.0.AND.MATGEO(IND(I)).LT.0 ) NEVOL= NEVOL+1
            IF( I.LT.0.AND.MATGEO(IND(I)).EQ.0 ) NESUR= NESUR-1
  549    CONTINUE
         WRITE(IOUT,'(/13H   NUMBERING ,I8,13H VOLUMES AND ,'//
     >           'I8,10H SURFACES.)') NVO-NEVOL,-NSU+NESUR
         WRITE(IOUT,'(17X,7HMINDIM=,10I8)') (MINDO(J),J=1,NCP)
         WRITE(IOUT,'(17X,7HMAXDIM=,10I8)') (MAXDO(J),J=1,NCP)
C
         DO 550 I= NSU-NESUR, NVO
             WRITE(IOUT,'(8H MATGEO(,I8,2H)=,I6,7H INDEX=,4I8)')
     >                    I, MATGEO(IND(I)), (INDEXO(J,IND(I)),J=1,4)
  550    CONTINUE
      ENDIF
C
      RETURN
      END
