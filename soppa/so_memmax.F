#include <single.h>
C
C  /* Deck so_memmax */
      SUBROUTINE SO_MEMMAX(SUBNM,MEMOFREE)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Stephan P. A. Sauer, 25 November 2003
C
C     PURPOSE: Collect and print information about memory usage in SOPPA
C              calculation.
C
#include <implicit.h>
#include <priunit.h>
C
#include <soppinf.h>
C
C
      CHARACTER*(*) SUBNM
      INTEGER       MEMOFREE
      PARAMETER (WORTOMB = 131072.D0)
C
      CHARACTER*56 LINE1,LINE2
      REAL*8       SOMEMMB
c      write (lupri,*) 'in memmax :',subnm,' : ',memofree
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_MEMMAX')
C
      LINE1 = '|=====================================================|'
      LINE2 = '|-----------------------------------------------------|'
C
C-------------------------
C     Initializes routine.
C-------------------------
C
      IF (SUBNM(1:5) .EQ. 'START') THEN
C
         DO I = 1, LSOSUB
            SOSUBNM(I) = '                '
         END DO
C
         ISOSUB = 1
C
         LWTOTAL = MEMOFREE
C
         LNM = LEN(SUBNM)
         SOSUBNM(ISOSUB) = SUBNM(1:LNM)
         SOMEMO(ISOSUB)  = LWTOTAL - MEMOFREE
C
C--------------------------------
C     Print memory statistics.
C--------------------------------
C
      ELSE IF (SUBNM(1:10) .EQ. 'STATISTICS') THEN
C
         WRITE(LUPRI,'(1X)')
         WRITE(LUPRI,'(1X,A)') LINE1
         WRITE(LUPRI,'(1X,A)')
     &    '|     Memory statistics for AO-SOPPA subroutines      |'
         WRITE(LUPRI,'(1X,A)') LINE1
         WRITE(LUPRI,'(1X,A)') 
     &    '|   Routine                  in Words         in MB   |'
         WRITE(LUPRI,'(1X,A)') LINE2
C
         DO I = 1, ISOSUB
C     
            SOMEMMB = SOMEMO(I)/WORTOMB
            WRITE(LUPRI,'(1X,A3,A16,6X,I12,6X,F8.1,A4)')
     &          '|  ',SOSUBNM(I),SOMEMO(I),SOMEMMB,'  |'
C
         END DO
         WRITE(LUPRI,'(1X,A)') LINE1
C
C----------------------------------------------------------------------
C        Compares memory requirements with previous maximum.
C        If the new memory requirement is larger, adds it to the array.
C----------------------------------------------------------------------
C
      ELSE
C
         IF ( (LWTOTAL - MEMOFREE) .GT. SOMEMO(ISOSUB) ) THEN
C     
            ISOSUB = ISOSUB + 1
C
            IF (ISOSUB .GT. LSOSUB) THEN
C
               WRITE(LUPRI,'(A,/,A,A)') 
     &             'ERROR: ISOSUB exceeds LSOSUB in SO_MEMMAX.',
     &             ' Recompile with larger value of LSOSUB ',
     &             ' in soppinf.h!' 
               CALL QUIT('ERROR:  dimension to small in SO_MEMMAX')
C
            END IF
C 
            LNM = LEN(SUBNM)
            SOSUBNM(ISOSUB) = SUBNM(1:LNM)
            SOMEMO(ISOSUB)  = LWTOTAL - MEMOFREE
C     
         END IF
C
      END IF
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_MEMMAX')
C
C
      RETURN
      END
