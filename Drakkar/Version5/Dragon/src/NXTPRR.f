*DECK NXTPRR
      FUNCTION NXTPRR(XYREC1,XYREC2,XYRECI)
*
*----------
*
*Purpose:
* Find the rectangle representing the intersection of two rectangles.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal.
*
*Author(s): G. Marleau.
*
*Reference:
*  G. Marleau,
*  \textsl{New Geometries Processing in DRAGON: The NXT: Module},
*  Report IGE-260, \'{E}cole Polytechnique de Montr\'{e}al,
*  Montr\'{e}al, 2005.
*
*Parameters: input
* XYREC1  spatial description of the first rectangle with
*         XYREC1(1) for left face, XYREC1(2) for right face,
*         XYREC1(3) for bottom face and XYREC1(4) for top face
*         positions.
* XYREC2  spatial description of the second rectangle with
*         XYREC2(1) for left face, XYREC2(2) for right face,
*         XYREC2(3) for bottom face and XYREC2(4) for top face
*         positions.
*Parameters: output
* NXTPRR  type of intersection between rectangles where
*         \begin{itemize}
*         \item \verb|NXTPRR|=0 means that there is no intersection
*         between the two regions;
*         \item \verb|NXTPRR|=1 means that there is an intersection between
*         between the two regions.
*         \end{itemize}
* XYRECI  spatial description of the intersection rectangle with
*         XYRECI(1) for left face, XYRECI(2) for right face,
*         XYRECI(3) for bottom face and XYRECI(4) for top face
*         positions.
*----------
*
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          NXTPRR
      DOUBLE PRECISION XYREC1(4),XYREC2(4)
      DOUBLE PRECISION XYRECI(4)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTPRR')
*----
*  Local variables
*----
      INTEGER          IFACE
*----
*  Find position of surface of intersection
*  -> for left and bottom faces, maximum $X$ and $Y$ location.
*  -> for right and top faces, minimum $X$ and $Y$ location.
*----
      NXTPRR=1
      DO 100 IFACE=1,3,2
        XYRECI(IFACE)=MAX(XYREC2(IFACE),XYREC1(IFACE))
 100  CONTINUE
      DO 101 IFACE=2,4,2
        XYRECI(IFACE)=MIN(XYREC2(IFACE),XYREC1(IFACE))
 101  CONTINUE
*----
*  Test if intersection is valid
*----
      IF(XYRECI(1) .GE. XYRECI(2) .OR.
     >   XYRECI(3) .GE. XYRECI(4)) NXTPRR=0
      RETURN
      END
