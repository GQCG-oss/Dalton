! FILE: blocks.h
! Used in hermit routines and routines referencing to hermit routines
! MBSISH has been added for multiple basis sets (WK/UniKA/31-10-2002).

      REAL*8   CENTSH

      INTEGER  MAXSHL, NLRGBL, NSMLBL, NHKTSH, KHKTSH, KCKTSH,
     &         ISTBSH, NUCOSH, NORBSH, NSTRSH, NCNTSH, NSETSH,
     &         IORBSB, NRCSH,  LCLASH, NO2INT, NLRBL,  ISYMBL,
     &         NSYMBL, MBSISH

      LOGICAL BIGVEC, SEGMEN, SEGMSH, SPHRSH

      INTEGER BLOCKSlast
      COMMON /BLOCKS/ CENTSH(MXSHEL,3),
     &                MAXSHL, NLRGBL, NSMLBL,
     &                NHKTSH(MXSHEL), KHKTSH(MXSHEL), KCKTSH(MXSHEL),
     &                ISTBSH(MXSHEL), NUCOSH(MXSHEL), NORBSH(MXSHEL),
     &                NSTRSH(MXSHEL), NCNTSH(MXSHEL), NSETSH(MXSHEL,2),
     &                IORBSB(0:MXCORB-1), NRCSH(MXSHEL), LCLASH(MXSHEL),
     &                NLRBL,ISYMBL(MXSHEL,8), NSYMBL, MBSISH(MXSHEL),
     &                BIGVEC, SEGMEN, SEGMSH(MXSHEL), SPHRSH(MXCORB),
     &   BLOCKSlast !  Very important:
      !  Always keep BLOCKSlast as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.
! end of blocks.h
