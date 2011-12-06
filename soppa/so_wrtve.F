C
C-----------------------------------------------------------------------
C  /* Deck so_wrtve */
      SUBROUTINE SO_WRTVE(ARRAY,LARRAY,ISYMTR,LABEL,FREQ,LU)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Andrea Ligabue, december 2003
C
C     Write solution or residual vector on file 
C     and all the needed information related to that vector
C
C     ARRAY      solution/residual vecotr
C     LARRAY     ARRAY length
C     ISYMTR     symmetry of the property
C     LABEL      property label
C     FREQ       frequency of the solution/residual
C     LU         file where we have to write the property
C
#include<implicit.h>
C
      DIMENSION ARRAY(LARRAY)
C
      CHARACTER*8 LABEL,LSTARS
      DATA LSTARS /'********'/
C
      WRITE(LU) LSTARS,LABEL,ISYMTR,LARRAY,FREQ
      WRITE(LU) ARRAY
C
      RETURN
C
      END
