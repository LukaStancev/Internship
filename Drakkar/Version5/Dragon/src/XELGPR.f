*DECK XELGPR
      SUBROUTINE XELGPR(  NDIM,   NTX,   NTY,   NTZ,  NTR,ISYMM,
     >                    NSUR,  NVOL,NTOTCL,MINDIM,MAXDIM,
     >                  KEYMRG, INDEX,MATALB)
C
C--------------------------    XELGPR    -------------------------------
C
C 1- PROGRAMME STATISTICS:
C           NAME: XELGPR
C      COMPONENT: EXCELL GEOMETRY ANALYSIS
C      CALLED BY: XELTRK
C        CALLING: ------
C        VERSION: 1.0
C       CREATION: 97/10/31
C         AUTHOR: G. MARLEAU
C            USE: PRINTS A SEMI-GRAPHICAL REPRESENTATION
C                 OF THE GEOMETRY COMPUTE ANNULAR SURFACE
C 2- INPUT
C   RANN   : ANNULAR RADIUS                                 R
C   PLANE  : CARTESIAN PLANE LOCATION                       R
C   NTX    : NUMBER OF X-MESH                               I
C   NTY    : NUMBER OF Y-MESH                               I
C   NTZ    : NUMBER OF Z-MESH                               I
C   NTR    : NUMBER OF R-MESH                               I
C   ISYMM  : FLAG FOR INTRINSIC SYMMETRY                    I        
C            ISYMM = 2 -> reflection plane normal to X axis
C            ISYMM = 4 -> reflection plane normal to Y axis
C            ISYMM = 8  -> reflection plane normal to X and Y axis  
C            ISYMM =16 -> reflection plane normal to Z axis
C            ISYMM =18 -> reflection plane normal to X and Z axis
C            ISYMM =20 -> reflection plane normal to Y and Z axis
C            ISYMM =24 -> reflection plane normal to X, Y and Z axis  
C   NSUR   : # OF SURFACES.                                 I
C   NVOL   : # OF ZONES.                                    I
C   NTOTCL : TOT # OF CYLINDERS IN EXACT GEOMETRY.          I
C   MINDIM : MIN INDEX VALUES FOR ALL AXES (RECT/CYL).      I(NTOTCL)
C   MAXDIM : MAX INDEX VALUES FOR ALL AXES (RECT/CYL)       I(NTOTCL)
C   KEYMRG : MERGING VECTOR        OF EXACT GEOMETRY.      I(NSUR:NVOL)
C   INDEX  : #ING OF SURFACES & ZONES.                   I(4,NSUR:NVOL)
C   MATALB : MATERIAL/ALBEDO                               I(NSUR:NVOL)
C
C--------------------------    XELGPR    -------------------------------
C
      IMPLICIT             NONE
C
      INTEGER              NDIM,   NTX,   NTY,   NTZ,   NTR,ISYMM,
     >                     NSUR,  NVOL,NTOTCL,
     >                     MINDIM(NTOTCL),
     >                     MAXDIM(NTOTCL),
     >                     KEYMRG(NSUR:NVOL),
     >                     INDEX(4,NSUR:NVOL),
     >                     MATALB(NSUR:NVOL),
     >                     NTC,IOUT
      PARAMETER          ( NTC=4,IOUT=6 )
      CHARACTER            CABS*16,CNON*16,CNAM*16
      CHARACTER            FMTB*24,FMTVS*24,FMTE*24
C
      INTEGER              MAXZ, MINZ, MAXY, MINY, MAXX, MINX, LNFMT,
     >                     ITRZ, NTRZ, IZ, IY, IX, ISURZ, ISURY, ISURX,
     >                     ISURG, ITC, IVS, ICL, IR
      INTEGER              IPX,IPY,IPZ,IPPZ
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: NAMNUM
C----
C  SCRATCH STORAGE ALLOCATION
C----
      ALLOCATE(NAMNUM(4,NTR+1,0:NTX+1,0:NTY+1,0:NTZ+1))
C
      CABS='          ABSENT'
      CNON='                '
C----
C  COMPUTE MEMORY SIZE REQUIRED
C----
      MAXZ=MAXDIM(3)
      MINZ=MINDIM(3)-1
      MAXY=MAXDIM(2)
      MINY=MINDIM(2)-1
      MAXX=MAXDIM(1)
      MINX=MINDIM(1)-1
      LNFMT=MAXX-MINX+1
      WRITE(FMTB ,5000) LNFMT*18-2
      WRITE(FMTVS,5001) LNFMT
      WRITE(FMTE ,5002) LNFMT*18-2
      ITRZ=0
      NTRZ=NTZ+1
      IF(NDIM .EQ. 2) THEN
        ITRZ=1
        NTRZ=1
      ENDIF
C----
C  INITIALIZE NAMNUM
C----
      DO 100 IZ=ITRZ,NTRZ
        IF(IZ .EQ. 0 .OR. IZ .EQ. NTZ+1) THEN
          ISURZ=1
        ELSE
          ISURZ=0
        ENDIF
        DO 101 IY=0,NTY+1
          IF(IY .EQ. 0 .OR. IY .EQ. NTY+1) THEN
            ISURY=1
          ELSE
            ISURY=0
          ENDIF
          DO 102 IX=0,NTX+1
            IF(IX .EQ. 0 .OR. IX .EQ. NTX+1) THEN
              ISURX=1
            ELSE
              ISURX=0
            ENDIF
C----
C  DETERMINE IF SURFACE REPRESENTS A LINE OR CORNER
C  FOR SURFACE
C----
            ISURG=ISURX*ISURY+ISURX*ISURZ+ISURY*ISURZ
            DO 103 IR=1,NTR+1
              IF(ISURG.EQ.0) THEN
C----
C REGION REQUIRED
C INITIALIZED TO ABSENT
C----
                READ(CABS,5010) (NAMNUM(ITC,IR,IX,IY,IZ),ITC=1,NTC)
              ELSE
C----
C  REGION NOT REQUIRED
C  INITIALIZE TO BLANK
C----
                READ(CNON,5010) (NAMNUM(ITC,IR,IX,IY,IZ),ITC=1,NTC)
              ENDIF
 103        CONTINUE
 102      CONTINUE
 101    CONTINUE
 100  CONTINUE
C----
C  SCAN ALL SURFACE AND REGIONS AND LOCATE POSITION
C  STORE ADEQUATE REGION NUMVER IN NAMNUM
C----
      DO 110 IVS=NSUR,NVOL
        IF(KEYMRG(IVS) .NE. 0) THEN
C----
C  POSITION IN X, Y AND Z LOCATED
C----
          IX=INDEX(1,IVS)-MINX
          IY=INDEX(2,IVS)-MINY
          IZ=INDEX(3,IVS)-MINZ
          IF(INDEX(4,IVS) .EQ. 0) THEN
C----
C  CARTESIAN POSITION
C  STORE AT LOCATION NTR+1
C----
            IR=NTR+1
            WRITE(CNAM,5011) MATALB(IVS),KEYMRG(IVS)
            READ(CNAM,5010)
     >        (NAMNUM(ITC,IR,IX,IY,IZ),ITC=1,NTC)
          ELSE
C----
C  ANNULAR POSITION
C  DETERMINE WHICH ANNULUS
C----
            DO 111 ICL=4,NTOTCL
              IF( INDEX(4,IVS) .GE. MINDIM(ICL)-1 .AND.
     >            INDEX(4,IVS) .LT. MAXDIM(ICL)      ) THEN
                IR=INDEX(4,IVS)-MINDIM(ICL)+2
C----
C  ANNULAR POSITION
C  STORE AT LOCATION IR
C----
                WRITE(CNAM,5011) MATALB(IVS),KEYMRG(IVS)
                READ(CNAM,5010)
     >            (NAMNUM(ITC,IR,IX,IY,IZ),ITC=1,NTC)
                GO TO 115
              ENDIF
 111        CONTINUE
 115        CONTINUE
          ENDIF
        ENDIF
 110  CONTINUE
C----
C PRINT HEADER
C----
      WRITE(IOUT,6000)
C----
C  PRINT NAMNUM MATRIX
C----
      IPZ=NTRZ
      IPY=0
      IPX=0
      IF(ISYMM .GE. 16) THEN
C----
C  Z SYMMETRY
C----
        IPZ=(NTZ+1)/2
        WRITE(IOUT,6010)
      ENDIF
      IF(ISYMM .EQ. 8 .OR. ISYMM .EQ. 24) THEN
C----
C  X AND Y SYMMETRY
C----
        IPX=NTX/2+1
        IPY=NTY/2+1
        WRITE(IOUT,6011)
      ELSE IF(ISYMM .EQ. 4 .OR. ISYMM .EQ. 20) THEN
C----
C  Y SYMMETRY
C----
        IPY=NTY/2+1
        WRITE(IOUT,6012)
      ELSE IF(ISYMM .EQ. 2 .OR. ISYMM .EQ. 18) THEN
C----
C  X SYMMETRY
C----
        IPX=NTX/2+1
        WRITE(IOUT,6013)
      ENDIF
C----
C Start test print
C      write(IOUT,7000) isymm,ntx,ipx,nty,ipy,ntz,ipz
C 7000 format(1x,'Test print:'/
C     >       1x,'Symmetry factor = ',i10/
C     >       1x,'ntx,ipx =',2i10/
C     >       1x,'nty,ipy =',2i10/
C     >       1x,'ntz,ipz =',2i10/
C     >       1x,'keymrg follows')
C      write(IOUT,7001) (ir,keymrg(ir),ir=-1,nsur,-1)
C      write(IOUT,7001) (ir,keymrg(ir),ir=1,nvol)
C 7001 format(10i10)
C Finish test print
C----
      DO 140 IZ=NTRZ,ITRZ,-1
        IPPZ=1
        IF(NDIM .EQ. 3) THEN
          IF(IZ .LE. IPZ) THEN
            IF(IZ .EQ. 0) THEN
              WRITE(IOUT,6001)
            ELSE IF(IZ .EQ. NTZ+1) THEN
              WRITE(IOUT,6002)
            ELSE
              WRITE(IOUT,6003) IZ
            ENDIF
          ELSE
            IPPZ=0
          ENDIF
        ELSE
          WRITE(IOUT,6004)
        ENDIF
        IF(IPPZ .EQ. 1) THEN 
          DO 141 IY=NTY+1,0,-1
            IF(IY .GE. IPY) THEN
              WRITE(IOUT,FMTB)
              DO 142 IR=NTR+1,1,-1
                WRITE(IOUT,FMTVS)
     <           ((NAMNUM(ITC,IR,IX,IY,IZ),ITC=1,NTC),IX=IPX,NTX+1)
 142          CONTINUE
            ENDIF
 141      CONTINUE
        ENDIF
        WRITE(IOUT,FMTE)
 140  CONTINUE
C----
C  SCRATCH STORAGE DEALLOCATION
C----
      DEALLOCATE(NAMNUM)
      RETURN
C----
C  FORMATS TO CREATE FORMATS
C----
 5000 FORMAT('(2X,',I10,'(1H-) )   ')
 5001 FORMAT('(   ',I10,'(2X,4A4)) ')
 5002 FORMAT('(2X,',I10,'(1H-)/)   ')
 5010 FORMAT(4A4)
 5011 FORMAT('(',I6,') ',I7)
C----
C  OTHER PRINT FORMATS
C----
 6000 FORMAT(//' PRINTING GEOMETRY DESCRIPTION BY PLANES '/
     >         ' ---- NOTATION USED:'/
     >10X,'NEGATIVE INTEGERS REPRESENT SURFACES'/
     >10X,'POSITIVE INTEGERS REPRESENT REGIONS'/
     >10X,'ABSENT MEANS THAT THE REGION OR SURFACE DOES NOT EXIST'/
     >10X,'FIRST LINE REPRESENTS REGION OR VOLUME IN CARTESIAN MESH'/
     >10X,'ADDITIONAL LINES REPRESENT REGION OR SURFACE IN ',
     >    'RADIAL MESH (OUTER TO INNER)'/
     >10X,'FOR 3-D MODEL, START WITH TOP Z-SURFACE ',
     >    'THEN GO DOWN ALONG Z-AXIS AND FINISH BY BOTTOM Z-SURFACE'/
     >10X,'FOR 2-D X-Y PLANE FIRST LINE IS FOR TOP Y-SURFACE ',
     >    'THEN GO DOWN ALONG Y-AXIS AND FINISH BY BOTTOM Y-SURFACE'/
     >10X,'FOR A LINE FIRST POINT IS FOR LEFT X-SURFACE ',
     >    'THEN INCREASE ALONG X-AXIS AND FINISH BY RIGHT X-SURFACE'/
     >10X,'MATERIAL AND ABLEDO NUMBERS ARE IN PARENTHESIS'/)
 6001 FORMAT(/' X-Y MESH ON BOTTOM Z-SURFACE')
 6002 FORMAT(/' X-Y MESH ON TOP Z-SURFACE')
 6003 FORMAT(/' X-Y MESH IN Z-PLANE = ',I10)
 6004 FORMAT(/' X-Y MESH')
 6010 FORMAT(/' GEOMETRY HAS CENTRAL Z SYMMETRY '/
     >        ' ONLY BOTTOM-Z PLANES  PRINTED')
 6011 FORMAT(/' GEOMETRY HAS CENTRAL X AND Y SYMMETRY '/
     >        ' ONLY TOP-Y RIGHT-X REGIONS PRINTED')
 6012 FORMAT(/' GEOMETRY HAS CENTRAL Y SYMMETRY '/
     >        ' ONLY TOP-Y REGIONS PRINTED')
 6013 FORMAT(/' GEOMETRY HAS CENTRAL X SYMMETRY '/
     >        ' ONLY RIGHT-X REGIONS PRINTED')
      END
