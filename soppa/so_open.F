#include <single.h>
C
C  /* Deck so_open */
      SUBROUTINE SO_OPEN(LUFILE,FILENM,LARRAY)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, May 1996
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Open and rewind the file LUFILE.
C             
#include <implicit.h>
#include <soppinf.h>
C
      CHARACTER*7 FILENM
C
      DTIME      = SECOND()
C
      IF (LARRAY .EQ. 0) RETURN
C
      LUFILE = 0
      CALL WOPEN2(LUFILE,FILENM,64,0)
C
      DTIME     = SECOND()  - DTIME
      SOORWC(3) = SOORWC(3) + DTIME
C
      RETURN
C
      END
