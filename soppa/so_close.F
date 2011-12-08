C
C  /* Deck so_close */
      SUBROUTINE SO_CLOSE(LUFILE,FILENM,STATUS)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, May 1996
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Close the file LUFILE.
C
#include <implicit.h>
#include <soppinf.h>
C
      CHARACTER*7 FILENM
      CHARACTER*4 STATUS
C
      DTIME      = SECOND()
C
      CALL WCLOSE2(LUFILE,FILENM,STATUS)
C
      DTIME     = SECOND()  - DTIME
      SOORWC(4) = SOORWC(4) + DTIME
C
      RETURN
C
      END
