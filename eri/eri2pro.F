!
!...   Copyright (c) 2011 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2011 (2011), see http://daltonprogram.org"
!...
!...   This source code is provided under a written licence and may be
!...   used, copied, transmitted, or stored only in accord with that
!...   written licence.
!...
!...   In particular, no part of the source code or compiled modules may
!...   be distributed outside the research group of the licence holder.
!...   This means also that persons (e.g. post-docs) leaving the research
!...   group of the licence holder may not take any part of Dalton,
!...   including modified files, with him/her, unless that person has
!...   obtained his/her own licence.
!...
!...   For further information, including how to get a licence, see:
!...      http://daltonprogram.org
!
!
#include <single.h>
C
C  /* Deck eripro */
      SUBROUTINE ERIPRO(FMAT,DMAT,NDMT,ISYMDM,IFCTYP,DOGDER,DOBDER,
     &                  IPRFCK,CCFBT,INDXBT,WORK,LWORK)
#include <implicit.h>
#include <iratdef.h>
#include <priunit.h>
#include <mxcent.h>
#include <dummy.h>
#include <aovec.h>
#include <maxaqn.h>
#include <maxorb.h>
C
      LOGICAL DOBDER, DOGDER
      DIMENSION FMAT(*), DMAT(*), WORK(LWORK), IFCTYP(NDMT),
     &          ISYMDM(NDMT), CCFBT(*), INDXBT(*)
#include <ccom.h>
#include <cbieri.h>
#include <ericom.h>
#include <erithr.h>
#include <erimem.h>
#include <aobtch.h>
#include <odbtch.h>
#include <symmet.h>
#include <infpar.h>
C
      IF (SLAVE) THEN
         IPRINT = 0
      ELSE
         CALL TIMER('START ',TIMSTR,TIMEND)
C
C        Initialization in ER2INI
C
         CALL ER2INI
C
         IPRINT = MAX(IPRERI,IPRFCK)
      END IF
C
      THRSH  = MAX(THRS,1.00D-15)
      NDMAT  = NDMT
C
      GDER = DOGDER
      BDER = DOBDER 
      CCRUN = .FALSE.
      EXPERI = .TRUE.
      UNDIFF = .FALSE.
      WRTINT = .FALSE.
C
      IF (.NOT. SLAVE) THEN
         IF (NDMT .GT. 1) THEN
            WRITE (LUPRI,'(/A,I5,/A)')
     &           ' NDMT .gt. 1 in ERIPRO.', NDMT,
     &           ' Process aborted'
            CALL QUIT('NDMT .gt. 1 in ERIPRO')
         END IF
      END IF
C
C     Memory
C
      MEMOK  = .TRUE.
      MEMADD = 0
      MODAB  = 0
      MODCD  = 0
C
C     AO batches
C     ==========
C
      CALL SETAOB(CCFBT,INDXBT,WORK,LWORK,IPRINT)
C
C     OD batches
C     ==========
C
C     This subroutine returns several arrays for each electron
C     starting at addresses K????1 and K????2. These are to be
C     transferred to ODCDRV.
C
      CALL ODBCHS(KODCL1,KODCL2,
     &            KODBC1,KODBC2,KRDBC1,KRDBC2,
     &            KODPP1,KODPP2,KRDPP1,KRDPP2,
     &            KFREE,LFREE,CCFBT,WORK,
     &            LWORK,IPRINT)
C
      IF (IPRINT .GT. 2) THEN
         WRITE (LUPRI,'(2(/,2X,A,I10))')
     &      ' Memory requirements for ODBCHS:',LWORK - LFREE,
     &      ' Memory left for ODCDRV:        ',LFREE
      END IF
C
      ICALL = 0
      CALL GETDST(ICALL,ICALL,IPRINT)
C
C     Select integrals to be calculated
C     =================================
C
      CALL PICKAO(IPRINT)
C
C     Information about distributions
C     ===============================
C
      CALL ERIDSI(INDXBT,IPRINT)
C
      KLAST = KFREE
      LWRK  = LFREE
C
C     Calculate integrals
C     ===================
C
      IF (SLAVE) THEN
         CALL ODCDRV(WORK(KODCL1),WORK(KODCL2),
     &               WORK(KODBC1),WORK(KODBC2),
     &               WORK(KRDBC1),WORK(KRDBC2),
     &               WORK(KODPP1),WORK(KODPP2),
     &               WORK(KRDPP1),WORK(KRDPP2),
     &               FMAT,DMAT,IDUMMY,IDUMMY,DUMMY,IDUMMY,CCFBT,INDXBT,
     &               WORK(KLAST),LWRK,IPRFCK)
      ELSE
         IF (.NOT.INTSKP) THEN
            CALL ODCDRV(WORK(KODCL1),WORK(KODCL2),
     &                  WORK(KODBC1),WORK(KODBC2),
     &                  WORK(KRDBC1),WORK(KRDBC2),
     &                  WORK(KODPP1),WORK(KODPP2),
     &                  WORK(KRDPP1),WORK(KRDPP2),
     &                  FMAT,DMAT,IDUMMY,IDUMMY,DUMMY,IDUMMY,CCFBT,
     &                  INDXBT,WORK(KLAST),LWRK,IPRINT)
C
C           Error message in case of insufficient memory
C
            IF (.NOT.MEMOK) THEN
               WRITE (LUPRI,'(//,1X,A,3(/,1X,A,I10))')
     &            ' Not enough memory for this run of ERIPRO.',
     &            ' Available memory in ERIPRO:',LWORK,
     &            ' Required memory for ERIPRO:',LWORK + MEMADD,
     &            ' Increase memory (LWORK) by:',MEMADD
               WRITE (LUPRI,'(/,1X,A,2I5)')
     &            ' Memory requirements largest for OD classes :',
     &              MODAB,MODCD
               CALL QUIT('Insufficient memory in ERIPRO.')
            END IF
         END IF
C
C        Print densities 
C        ===============
C
         IF (IPRINT.GT.4) THEN
            CALL HEADER('Density matrix in ERIPRO',-1)
            KSTR = 1
            DO I = 1, NDMT
               WRITE (LUPRI,'(//,1X,A,I3)') ' Density matrix No.',I
               CALL OUTPUT(DMAT(KSTR),1,NBASE,1,NBASE,NBASE,
     &                     NBASE,1,LUPRI)
               KSTR = NBASE*NBASE
            END DO
         END IF
C
         CALL TIMER('ERIPRO',TIMSTR,TIMEND)
         CALL FLSHFO(LUPRI)
      END IF
C
      RETURN
      END
