#if defined (VAR_STAR2)
      INTEGER*2 MJWOP
#endif
#if defined (VAR_STAR4)
      INTEGER*4 MJWOP
#endif
      LOGICAL QRREST
C     MAXWOP = maximum number of rotations (dimension of JWOP)
C     MAXWOP IS DEFINED IN INFVAR WHICH HAS TO PRECEDE
C     THE CALL OF QRINF
      COMMON /QRINF/  MZVAR(8), MZYVAR(8), MZCONF(8),
     *                MZYCON(8), MZWOPT(8), MZYWOP(8), MSYMA ,MSYMB ,
     *                MSYMC, MZYVMX, MZWOPH(8),MZVARH(8),
     *                QRREST
