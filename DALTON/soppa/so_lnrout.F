C
C  /* Deck so_lnrout */
      SUBROUTINE SO_LNROUT(POLDD,POLDL,POLDA,POLVL,POLVV,FOVIBG,
     &                     WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Stephan P. A. Sauer, February 1999
C     Rasmus Faber, October 2015, rewrite for merge with modern Dalton
C
C     PURPOSE: Driver for printing final output for the RPA, 
C              SOPPA, and SOPPA(CCSD) polarizabilities.
C
#include "implicit.h"
#include "priunit.h"
#include "soppinf.h"
#include "cbilnr.h"
C
      INTEGER, PARAMETER :: NMODEL = 4
C         
      DIMENSION POLDD(2,3,3,NFRVAL,NMODEL)!,POLDQ(2,3,3,3,NFRVAL,4)
      DIMENSION POLDL(2,3,3,NFRVAL,NMODEL),POLDA(2,3,3,NFRVAL,NMODEL)
C These are new, check dimension
      DIMENSION POLVL(2,3,3,NFRVAL,NMODEL),POLVV(2,3,3,NFRVAL,NMODEL)
      DIMENSION WORK(LWORK)

      LOGICAL   JOBS(4)
      CHARACTER(LEN=11), PARAMETER, DIMENSION(NMODEL) :: METHOD_NAMES = 
     &     (/ 'RPA        ', 'SOPPA      ', 'SOPPA(CCSD)', 
     &        'DCRPA      '  /)

C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_LNROUT')
C      
      JOBS = (/ AORPA, AOSOP, AOSOC, DCRPA /)
C
C===============================================
C     LOOP OVER METHODS
C===============================================
C
      DO I = 1, NMODEL         
         IF(.NOT.JOBS(I) ) CYCLE 
C
C-----------------------------------------------
C     Write calculated properties to output.
C-----------------------------------------------
C
         WRITE(LUPRI,9000)
         WRITE(LUPRI,'(30X,A,A)') TRIM(METHOD_NAMES(I)),' results:'
         WRITE(LUPRI,9001)
C
         CALL LNROUT (POLDD(1,1,1,1,I),POLDL(1,1,1,1,I),
     &                POLDA(1,1,1,1,I),POLVL(1,1,1,1,I),
     &                POLVV(1,1,1,1,I),FOVIBG,
     &                IPRSOP,WORK,LWORK)
C
      ENDDO
C
C===============================================
C     Write calculated RPA properties to output.
C===============================================
C
C      IF (AORPA) THEN
C
C         WRITE(LUPRI,9000)
C         WRITE(LUPRI,'(30X,A)') ' RPA results:'
C         WRITE(LUPRI,9001)
C
C         CALL LNROUT (POLDD(1,1,1,1,1),POLDL(1,1,1,1,1),
C    &                POLDA(1,1,1,1,1),POLVL(1,1,1,1,1),
C    &                POLVV(1,1,1,1,1),FOVIBG,
C    &                IPRSOP,WORK,LWORK)
C
C     END IF
C
C=================================================
C     Write calculated SOPPA properties to output.
C=================================================
C
C     IF (AOSOP) THEN
C
C        WRITE(LUPRI,9000)
C        WRITE(LUPRI,'(30X,A)') ' SOPPA results:'
C        WRITE(LUPRI,9001)
C
C        CALL LNROUT (POLDD(1,1,1,1,2),POLDL(1,1,1,1,2),
C    &                POLDA(1,1,1,1,2),IPRSOP)
C
C     END IF
C
C=======================================================
C     Write calculated SOPPA(CCSD) properties to output.
C=======================================================
C
C     IF (AOSOC) THEN
C
C        WRITE(LUPRI,9000)
C        WRITE(LUPRI,'(30X,A)') ' SOPPA(CCSD) results:'
C        WRITE(LUPRI,9001)
C
C        CALL LNROUT (POLDD(1,1,1,1,3),POLDL(1,1,1,1,3),
C    &                POLDA(1,1,1,1,3),IPRSOP)
C
C     END IF
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL FLSHFO(LUPRI)
C
      CALL QEXIT('SO_LNROUT')
C
      RETURN
C
 9000 FORMAT(//' =========================================',
     &       '=============================')
 9001 FORMAT(' -----------------------------------------',
     &       '-----------------------------')
      END
