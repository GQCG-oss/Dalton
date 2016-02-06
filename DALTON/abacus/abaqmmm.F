!
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2016 (2015), see http://daltonprogram.org"
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
C FILE: abacus/abaqmmm.F
C*******************************************************************
C/* deck qmmmfirst */
      SUBROUTINE QMMMFIRST(VEC2,IPRINT,INTPRI,WORK,LWORK)
C
c This routine adds the MM environmental contribution 
c to the 1st order (NCOMP=3) magnetic field perturbation due to the use of
c London orbitals
c Lund, Nov. 2008, J. Kongsted
c Based on QM3FIRST by K. Ruud
c
c VEC2 (OUTPUT)  : contribution to be added to perturbed Fock matrix
c NCOMP (INPUT)  : number of independent perturbation (3)

#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "orgcom.h"
#include "nuclei.h"
#include "qm3.h"
#include "qmmm.h"
#include "mmtimes.h"
#include "inforb.h"
#include "infpar.h"
C
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL,LOCDEB
      LOGICAL LCOMB
      DIMENSION WORK(*),VEC2(*)
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)

      CALL QENTER('QMMMFIRST')
      IF (MMTIME) DTIME = SECOND()

      LOCDEB = .FALSE.

C     If we have BOTH permanent dipoles AND induced dipoles then these
C     can be added together and we can treat the effect of those
C     simultanously. Then no seperate loop for dipoles and
C     polarizabilities but only one combined loop.

      LCOMB = .FALSE.
      IF ( (IPOLTP .GT. 0) .AND. (NMULT .GE. 1) ) LCOMB = .TRUE.
C     SPLNMR overrules these settings   
      IF (SPLNMR) LCOMB = .FALSE.

C     We don't have these integrals yet !     
      IF (NMULT .GE. 2) THEN
        CALL QUIT('Quadrupole contr. not implemented in QMMMFIRST')     
      ENDIF

      OBKPX = DIPORG(1)
      OBKPY = DIPORG(2)
      OBKPZ = DIPORG(3)

      CALL DZERO(VEC2,3*N2BASX)

C     1. Construct the charge contribution for all MM centers

      IF (NMULT .GE. 0) THEN

#if defined(VAR_MPI)
        IF (NODTOT.GE.1)THEN
           CALL QMMMFIRST_M1(VEC2,IPRINT,INTPRI,WORK(1),LWORK)
        ELSE
#endif
        KMAT = 1
        KLAST = KMAT + 3*N2BASX
        LWRK = LWORK - KLAST + 1

        IF (LWRK .LT. 0) CALL ERRWRK('QMMMFIRST',-KLAST,LWORK)

        CALL DZERO(WORK(KMAT),3*N2BASX)

        DO 100 I = 1,MMCENT

           IF (ABS(MUL0MM(I)) .LE. THRMM) GOTO 100

           KPATOM = 0
           NOSIM = 3
           TOFILE = .FALSE.
           TRIMAT = .FALSE.
           EXP1VL = .FALSE.
           DIPORG(1) = MMCORD(1,I)
           DIPORG(2) = MMCORD(2,I)
           DIPORG(3) = MMCORD(3,I)

           CALL DZERO(WORK(KMAT),3*N2BASX)

           CALL GET1IN(WORK(KMAT),'PCMBSOL',NOSIM,WORK(KLAST),
     &                 LWRK,LABINT,INTREP,INTADR,I,TOFILE,
     &                 KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,INTPRI)

           IF ( (INTPRI .GE. 12) .OR. (LOCDEB) ) THEN
              WRITE (LUPRI,'(/A,I3,A)')
     *             ' N(',I,')_ao matrix (Bx): '
              CALL OUTPUT(WORK(KMAT),1,NBAST,1,NBAST,
     &                    NBAST,NBAST,1,LUPRI)
              WRITE (LUPRI,'(/A,I3,A)')
     *             ' N(',I,')_ao matrix (By): '
              CALL OUTPUT(WORK(KMAT+N2BASX),1,NBAST,1,NBAST,
     &                    NBAST,NBAST,1,LUPRI)
              WRITE (LUPRI,'(/A,I3,A)')
     *             ' N(',I,')_ao matrix (Bz): '
              CALL OUTPUT(WORK(KMAT+2*N2BASX),1,NBAST,1,NBAST,
     &                    NBAST,NBAST,1,LUPRI)
           END IF

           FAC1 = -1.0D0*MUL0MM(I)
           CALL DAXPY(3*N2BASX,FAC1,WORK(KMAT),1,VEC2,1)

  100   CONTINUE

#if defined(VAR_MPI)
        ENDIF ! (NODTOT.GE.1) ... ELSE 
#endif
      ENDIF !( NMULT .GE. 0)

C     2. Construct the permanent dipole contribution for all MM centers

      IF ( (NMULT .GE. 1) .AND. (.NOT. LCOMB) ) THEN

        IF (LOCDEB) THEN
          WRITE(LUPRI,*) 'London gradient from permanent dipoles'
        ENDIF

#if defined(VAR_MPI)
        IF (NODTOT .GE. 1) THEN
           CALL QMMMFIRST_M2(VEC2,IPRINT,INTPRI,WORK(1),LWORK)
        ELSE
#endif
        KMAT = 1
        KLAST = KMAT + 9*N2BASX
        LWRK = LWORK - KLAST + 1

        IF (LWRK .LT. 0) CALL ERRWRK('QMMMFIRST',-KLAST,LWORK)

        CALL DZERO(WORK(KMAT),9*NNBASX)

        DO 200 I = 1,MMCENT

          DNORM2 = MUL1MM(1,I)**2+MUL1MM(2,I)**2+MUL1MM(3,I)**2
          DNORM = SQRT(DNORM2)

          IF (ABS(DNORM) .LE. THRMM) GOTO 200

          CALL DZERO(WORK(KMAT),9*NNBASX)

          KPATOM = 0
          NOSIM = 9
          TOFILE = .FALSE.
          TRIMAT = .FALSE.
          EXP1VL = .FALSE.
          DIPORG(1) = MMCORD(1,I)
          DIPORG(2) = MMCORD(2,I)
          DIPORG(3) = MMCORD(3,I)

          CALL GET1IN(WORK(KMAT),'EFIELB1',NOSIM,WORK(KLAST),
     &                LWRK,LABINT,INTREP,INTADR,I,TOFILE,
     &                KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,INTPRI)
          IF ( (INTPRI .GE. 15) .OR. (LOCDEB) ) THEN
             WRITE (LUPRI,'(/A)') ' Efield_ao Ex Bx matrix:'
             CALL OUTPUT(WORK(KMAT),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ex By matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ex Bz matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*2),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ey Bx matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*3),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ey By matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*4),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ey Bz matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*5),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ez Bx matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*6),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ez By matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*7),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ez Bz matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*8),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
          END IF

          FACX = -1.0D0*MUL1MM(1,I)
          FACY = -1.0D0*MUL1MM(2,I)
          FACZ = -1.0D0*MUL1MM(3,I)

          CALL DAXPY(N2BASX,FACX,WORK(KMAT),1,VEC2,1)
          CALL DAXPY(N2BASX,FACX,WORK(KMAT+N2BASX),1,
     &               VEC2(1+N2BASX),1)
          CALL DAXPY(N2BASX,FACX,WORK(KMAT+2*N2BASX),1,
     &               VEC2(1+2*N2BASX),1)
          CALL DAXPY(N2BASX,FACY,WORK(KMAT+3*N2BASX),1,
     &               VEC2,1)
          CALL DAXPY(N2BASX,FACY,WORK(KMAT+4*N2BASX),1,
     &               VEC2(1+N2BASX),1)
          CALL DAXPY(N2BASX,FACY,WORK(KMAT+5*N2BASX),1,
     &               VEC2(1+2*N2BASX),1)
          CALL DAXPY(N2BASX,FACZ,WORK(KMAT+6*N2BASX),1,
     &               VEC2,1)
          CALL DAXPY(N2BASX,FACZ,WORK(KMAT+7*N2BASX),1,
     &               VEC2(1+N2BASX),1)
          CALL DAXPY(N2BASX,FACZ,WORK(KMAT+8*N2BASX),1,
     &               VEC2(1+2*N2BASX),1)

  200   CONTINUE

#if defined(VAR_MPI)
        ENDIF ! (NODTOT.GE.1) ... ELSE 
#endif
      ENDIF !( NMULT .GE. 1)

C     3. Construct the induced dipole contribution for all MM centers
C        NOTE: This could be done more effectively by adding the
C        permanent dipoles and induced dipoles and then performing the
C        London orbital contribution. See the LCOMB option below.

      IF ( (IPOLTP .GT. 0) .AND. (.NOT. LCOMB) ) THEN

        IF (LOCDEB) THEN
          WRITE(LUPRI,*) 'London gradient from induced dipoles'
        ENDIF

        KMAT    = 1
        KINDMOM = KMAT    + 9*N2BASX ! INDMOM contains the ind. dips.
        KLAST   = KINDMOM + 3*NNZAL
        LWRK    = LWORK   - KLAST + 1 

        IF (LWRK .LT. 0) CALL ERRWRK('QMMMFIRST',-KLAST,LWORK)

        CALL DZERO(WORK(KMAT),9*NNBASX)
        CALL DZERO(WORK(KINDMOM),3*NNZAL)

        IF (LOCDEB) THEN
           WRITE(LUPRI,*)
           WRITE(LUPRI,*) 'Ind. dips read from file in QMMMFIRST'
           WRITE(LUPRI,*)
        ENDIF 

        CALL GET_FROM_FILE_1('INDUCED_DIPOLES',NNZAL,WORK(KINDMOM))

#if defined(VAR_MPI)
        IF (NODTOT.GE.1) THEN
           CALL QMMMFIRST_M3(WORK(KINDMOM),VEC2,
     &          IPRINT,INTPRI,WORK(KLAST),LWRK)
        ELSE
#endif
        IL = 0 ! Counter for induced dipole moments

        DO 300 I = 1,MMCENT

          IF (ZEROAL(I) .EQ. -1) GOTO 300

          CALL DZERO(WORK(KMAT),9*NNBASX)

          KPATOM = 0
          NOSIM = 9
          TOFILE = .FALSE.
          TRIMAT = .FALSE.
          EXP1VL = .FALSE.
          DIPORG(1) = MMCORD(1,I)
          DIPORG(2) = MMCORD(2,I)
          DIPORG(3) = MMCORD(3,I)

          CALL GET1IN(WORK(KMAT),'EFIELB1',NOSIM,WORK(KLAST),
     &                LWRK,LABINT,INTREP,INTADR,I,TOFILE,
     &                KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,INTPRI)
          IF ( (INTPRI .GE. 15) .OR. (LOCDEB) ) THEN
             WRITE (LUPRI,'(/A)') ' Efield_ao Ex Bx matrix:'
             CALL OUTPUT(WORK(KMAT),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ex By matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ex Bz matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*2),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ey Bx matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*3),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ey By matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*4),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ey Bz matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*5),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ez Bx matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*6),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ez By matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*7),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ez Bz matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*8),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
          END IF

          FACX = -1.0D0*WORK(KINDMOM+IL+0)
          FACY = -1.0D0*WORK(KINDMOM+IL+1)
          FACZ = -1.0D0*WORK(KINDMOM+IL+2)

          CALL DAXPY(N2BASX,FACX,WORK(KMAT),1,VEC2,1)
          CALL DAXPY(N2BASX,FACX,WORK(KMAT+N2BASX),1,
     &               VEC2(1+N2BASX),1)
          CALL DAXPY(N2BASX,FACX,WORK(KMAT+2*N2BASX),1,
     &               VEC2(1+2*N2BASX),1)
          CALL DAXPY(N2BASX,FACY,WORK(KMAT+3*N2BASX),1,
     &               VEC2,1)
          CALL DAXPY(N2BASX,FACY,WORK(KMAT+4*N2BASX),1,
     &               VEC2(1+N2BASX),1)
          CALL DAXPY(N2BASX,FACY,WORK(KMAT+5*N2BASX),1,
     &               VEC2(1+2*N2BASX),1)
          CALL DAXPY(N2BASX,FACZ,WORK(KMAT+6*N2BASX),1,
     &               VEC2,1)
          CALL DAXPY(N2BASX,FACZ,WORK(KMAT+7*N2BASX),1,
     &               VEC2(1+N2BASX),1)
          CALL DAXPY(N2BASX,FACZ,WORK(KMAT+8*N2BASX),1,
     &               VEC2(1+2*N2BASX),1)

          IL = IL+3

  300   CONTINUE

#if defined(VAR_MPI)
        ENDIF ! (NODTOT.GE.1) ... ELSE 
#endif
      ENDIF !(IPOLTP .GT. 0)

C     3a. If both permanent dipoles AND polarizabilities construct TOTAL 
C         dipole contribution for all MM centers

      IF (LCOMB)  THEN

        IF (LOCDEB) THEN
          WRITE(LUPRI,*) 'London gradient from total dipoles'
        ENDIF

        KMAT    = 1
        KINDMOM = KMAT    + 9*N2BASX ! INDMOM contains the ind. dips.
        KTOTMOM = KINDMOM + 3*NNZAL  ! TOTMOM contains the TOTAL dips.
        KLAST   = KTOTMOM + 3*MMCENT
        LWRK    = LWORK   - KLAST + 1

        IF (LWRK .LT. 0) CALL ERRWRK('QMMMFIRST',-KLAST,LWORK)

        CALL DZERO(WORK(KMAT),9*NNBASX)
        CALL DZERO(WORK(KINDMOM),3*NNZAL)
        CALL DZERO(WORK(KTOTMOM),3*MMCENT)

        IF (LOCDEB) THEN
           WRITE(LUPRI,*)
           WRITE(LUPRI,*) 'Ind. dips read from file in QMMMFIRST'
           WRITE(LUPRI,*)
        ENDIF 

        CALL GET_FROM_FILE_1('INDUCED_DIPOLES',NNZAL,WORK(KINDMOM))

        IL = 0
        IJ = 0
        DO 400 I = 1,MMCENT
          IF (ZEROAL(I) .EQ. -1) THEN
            WORK(KTOTMOM+IL+0) = MUL1MM(1,I)
            WORK(KTOTMOM+IL+1) = MUL1MM(2,I)
            WORK(KTOTMOM+IL+2) = MUL1MM(3,I)
          ELSE IF (ZEROAL(I) .EQ. 1) THEN 
            WORK(KTOTMOM+IL+0) = MUL1MM(1,I) + WORK(KINDMOM+IJ+0)
            WORK(KTOTMOM+IL+1) = MUL1MM(2,I) + WORK(KINDMOM+IJ+1)
            WORK(KTOTMOM+IL+2) = MUL1MM(3,I) + WORK(KINDMOM+IJ+2)
            IJ=IJ+3
          ENDIF
          IL=IL+3
  400   CONTINUE

#if defined(VAR_MPI)
        IF (NODTOT.GE.1)THEN
C     If we are doing a parallel calculation
           CALL QMMMFIRST_M4(WORK(KTOTMOM),VEC2,
     &          IPRINT,INTPRI,WORK(KLAST),LWRK)
        ELSE
#endif
Cahs        IL = 0 ! Counter for dipole moments

        DO 500 I = 1,MMCENT
          IL = 3*(I - 1)
C         See if the total dipole moment at this site is zero
          DNORM2 = WORK(KTOTMOM+IL+0)*WORK(KTOTMOM+IL+0)
     &           + WORK(KTOTMOM+IL+1)*WORK(KTOTMOM+IL+1)
     &           + WORK(KTOTMOM+IL+2)*WORK(KTOTMOM+IL+2)
          DNORM = SQRT(DNORM2)
          IF (ABS(DNORM) .LE. THRMM) GOTO 500

          CALL DZERO(WORK(KMAT),9*NNBASX)

          KPATOM = 0
          NOSIM = 9
          TOFILE = .FALSE.
          TRIMAT = .FALSE.
          EXP1VL = .FALSE.
          DIPORG(1) = MMCORD(1,I)
          DIPORG(2) = MMCORD(2,I)
          DIPORG(3) = MMCORD(3,I)

          CALL GET1IN(WORK(KMAT),'EFIELB1',NOSIM,WORK(KLAST),
     &                LWRK,LABINT,INTREP,INTADR,I,TOFILE,
     &                KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,INTPRI)

          IF ( (INTPRI .GE. 15) .OR. (LOCDEB) ) THEN
             WRITE (LUPRI,'(/A)') ' Efield_ao Ex Bx matrix:'
             CALL OUTPUT(WORK(KMAT),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ex By matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ex Bz matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*2),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ey Bx matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*3),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ey By matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*4),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ey Bz matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*5),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ez Bx matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*6),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ez By matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*7),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ez Bz matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*8),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
          END IF

          FACX = -1.0D0*WORK(KTOTMOM+IL+0)
          FACY = -1.0D0*WORK(KTOTMOM+IL+1)
          FACZ = -1.0D0*WORK(KTOTMOM+IL+2)

          CALL DAXPY(N2BASX,FACX,WORK(KMAT),1,VEC2,1)
          CALL DAXPY(N2BASX,FACX,WORK(KMAT+N2BASX),1,
     &               VEC2(1+N2BASX),1)
          CALL DAXPY(N2BASX,FACX,WORK(KMAT+2*N2BASX),1,
     &               VEC2(1+2*N2BASX),1)
          CALL DAXPY(N2BASX,FACY,WORK(KMAT+3*N2BASX),1,
     &               VEC2,1)
          CALL DAXPY(N2BASX,FACY,WORK(KMAT+4*N2BASX),1,
     &               VEC2(1+N2BASX),1)
          CALL DAXPY(N2BASX,FACY,WORK(KMAT+5*N2BASX),1,
     &               VEC2(1+2*N2BASX),1)
          CALL DAXPY(N2BASX,FACZ,WORK(KMAT+6*N2BASX),1,
     &               VEC2,1)
          CALL DAXPY(N2BASX,FACZ,WORK(KMAT+7*N2BASX),1,
     &               VEC2(1+N2BASX),1)
          CALL DAXPY(N2BASX,FACZ,WORK(KMAT+8*N2BASX),1,
     &               VEC2(1+2*N2BASX),1)

Cahs          IL = IL+3

  500   CONTINUE

#if defined(VAR_MPI)
        ENDIF ! (NODTOT.GE.1) ... ELSE 
#endif
      ENDIF !(LCOMB)

C     **print out section**
C
      IF (IPRINT .GT. 5) THEN
         CALL AROUND(
     &        'First order qmmm contr. to gradient in QMMMFIRST')
         WRITE (LUPRI,'(2X,A)') 'X coordinate'
         CALL OUTPUT(VEC2(1),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         WRITE (LUPRI,'(2X,A)') 'Y coordinate'
         CALL OUTPUT(VEC2(1+N2BASX),1,NBAST,1,NBAST,NBAST,NBAST,1,
     &        LUPRI)
         WRITE (LUPRI,'(2X,A)') 'Z coordinate'
         CALL OUTPUT(VEC2(1+2*N2BASX),1,NBAST,1,NBAST,
     &        NBAST,NBAST,1,LUPRI)
      END IF
C
      DIPORG(1) = OBKPX
      DIPORG(2) = OBKPY
      DIPORG(3) = OBKPZ
C
      IF (MMTIME) THEN
        DTIME = SECOND() - DTIME
        TMMFIRST = TMMFIRST + DTIME
      ENDIF
C
      CALL QEXIT('QMMMFIRST')
      RETURN
      END
C*****************************************************************************
C/* deck qmmmb2 */
      SUBROUTINE QMMMB2(EXPVAL,DENMAT,WORK,LWORK,IPRINT)
C
c This routine adds the /MM contribution 
c to the 2nd-order magnetic field perturbation due to the use of
c London orbitals. Works in combination with new QMMM code.
c Lund, Nov. 2008, J. Kongsted, based on QM3B3 by K. Ruud.
c
c EXPVAL (OUTPUT): contribution to be added to diamagnetic magnetizability
c
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "orgcom.h"
#include "nuclei.h"
#include "qm3.h"
#include "qmmm.h"
#include "mmtimes.h"
#include "inforb.h"
#include "infpar.h"
C
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL
      LOGICAL LOCDEB, LCOMB
      DIMENSION WORK(*), EXPVAL(6), DENMAT(*)
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
C

      CALL QENTER('QMMMB2')
      IF (MMTIME) DTIME = SECOND()

      LOCDEB = .FALSE.

C     If we have BOTH permanent dipoles AND induced dipoles then these
C     can be added together and we can treat the effect of those
C     simultanously. Then no seperate loop for dipoles and
C     polarizabilities but only one combined loop.

      LCOMB = .FALSE.
      IF ( (IPOLTP .GT. 0) .AND. (NMULT .GE. 1) ) LCOMB = .TRUE.

C     SPLNMR overrules these settings   
      IF (SPLNMR) LCOMB = .FALSE.

C     We don't have these integrals yet !     
      IF (NMULT .GE. 2) THEN
        CALL QUIT('Quadrupole contr. not implemented in QMMMB2')
      ENDIF

      OBKPX = DIPORG(1)
      OBKPY = DIPORG(2)
      OBKPZ = DIPORG(3)

      CALL DZERO(EXPVAL,6)

C     1. Construct the charge contribution for all MM centers

      IF (NMULT .GE. 0) THEN

#if defined(VAR_MPI)
        IF (NODTOT .GE. 1) THEN
           CALL QMMMB2_M1(EXPVAL,IPRINT,DENMAT,WORK(1),LWORK)
        ELSE
#endif
        KTMP  = 1
        KLAST = KTMP  + 6
        LWRK  = LWORK - KLAST + 1

        IF (LWRK .LT. 0) CALL ERRWRK('QMMMB2',-KLAST,LWORK)

        CALL DZERO(WORK(KTMP),6)

        DO 100 I = 1,MMCENT

          IF (ABS(MUL0MM(I)) .LE. THRMM) GOTO 100

          CALL DZERO(WORK(KTMP),6)

          KPATOM = 0
          NOSIM = 0
          TOFILE = .FALSE.
          TRIMAT = .FALSE.
          EXP1VL = .TRUE.
          DIPORG(1) = MMCORD(1,I)
          DIPORG(2) = MMCORD(2,I)
          DIPORG(3) = MMCORD(3,I)

          CALL GET1IN(DUMMY,'PCMB2SL',NOSIM,WORK(KLAST),
     &                LWRK,LABINT,INTREP,INTADR,I,TOFILE,
     &                KPATOM,TRIMAT,WORK(KTMP),EXP1VL,DENMAT,
     &                IPRINT)

          IF (IPRINT .GE. 4) THEN
             WRITE (LUPRI,'(/A,I3,A)')
     *       ' N(',I,') QMMM correction matrix '
             CALL OUTPAK(WORK(KTMP),3,1,LUPRI)
          END IF

          FAC1 = -1.0D0*MUL0MM(I)
          CALL DAXPY(6,FAC1,WORK(KTMP),1,EXPVAL,1)

  100   CONTINUE

#if defined(VAR_MPI)
        ENDIF ! (NODTOT.GE.1) ... ELSE 
#endif
      ENDIF !( NMULT .GE. 0)

C     2. Construct the permanent dipole contribution for all MM centers

      IF ( (NMULT .GE. 1) .AND. (.NOT. LCOMB) ) THEN

        IF (LOCDEB) THEN
          WRITE(LUPRI,*) 'London gradient from permanent dipoles'
        ENDIF

#if defined(VAR_MPI)
        IF (NODTOT .GE. 1) THEN
           CALL QMMMB2_M2(EXPVAL,IPRINT,DENMAT,WORK(1),LWORK)
        ELSE
#endif
        KTMP  = 1
        KLAST = KTMP  + 18
        LWRK  = LWORK - KLAST + 1

        IF (LWRK .LT. 0) CALL ERRWRK('QMMMB2',-KLAST,LWORK)

        CALL DZERO(WORK(KTMP),18)

        DO 200 I = 1,MMCENT

          DNORM2 = MUL1MM(1,I)**2+MUL1MM(2,I)**2+MUL1MM(3,I)**2
          DNORM = SQRT(DNORM2)

          IF (ABS(DNORM) .LE. THRMM) GOTO 200

          CALL DZERO(WORK(KTMP),18)

          KPATOM = 0
          NOSIM = 0
          TOFILE = .FALSE.
          TRIMAT = .FALSE.
          EXP1VL = .TRUE.
          DIPORG(1) = MMCORD(1,I)
          DIPORG(2) = MMCORD(2,I)
          DIPORG(3) = MMCORD(3,I)
          CALL GET1IN(DUMMY,'EFIELB2',NOSIM,WORK(KLAST),
     &                LWRK,LABINT,INTREP,INTADR,I,TOFILE,
     &                KPATOM,TRIMAT,WORK(KTMP),EXP1VL,DENMAT,
     &                IPRINT)

          IF (IPRINT .GE. 4) THEN
             WRITE (LUPRI,'(/A)') ' Efield Ex 2.order matrix:'
             CALL OUTPAK(WORK(KTMP),3,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield Ey 2.order matrix:'
             CALL OUTPAK(WORK(KTMP+6),3,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield Ez 2.order matrix:'
             CALL OUTPAK(WORK(KTMP+12),3,1,LUPRI)
          END IF

          FACX = -1.0D0*MUL1MM(1,I)
          FACY = -1.0D0*MUL1MM(2,I)
          FACZ = -1.0D0*MUL1MM(3,I)

          DO I2 = 1,6
             EXPVAL(I2) = EXPVAL(I2) + WORK(KTMP+I2-1)*FACX +
     &                    WORK(KTMP+I2+5)*FACY + WORK(KTMP+I2+11)*FACZ
          END DO

  200   CONTINUE

#if defined(VAR_MPI)
        ENDIF ! (NODTOT.GE.1) ... ELSE 
#endif
      ENDIF !( NMULT .GE. 1)

C     3. Construct the induced dipole contribution for all MM centers
C        NOTE: This could be done more effectively by adding the
C        permanent dipoles and induced dipoles and then performing the
C        London orbital contribution. See the LCOMB option below.

      IF ( (IPOLTP .GT. 0) .AND. (.NOT. LCOMB) ) THEN

        IF (LOCDEB) THEN
          WRITE(LUPRI,*) 'London gradient from induced dipoles'
        ENDIF

        KTMP    = 1
        KINDMOM = KTMP    + 18        !INDMOM contains the ind. dips.
        KLAST   = KINDMOM + 3*NNZAL
        LWRK    = LWORK   - KLAST + 1

        IF (LWRK .LT. 0) CALL ERRWRK('QMMMB2',-KLAST,LWORK)

        CALL DZERO(WORK(KTMP),18)
        CALL DZERO(WORK(KINDMOM),3*NNZAL)

        IF (LOCDEB) THEN
           WRITE(LUPRI,*)
           WRITE(LUPRI,*) 'Ind. dips read from file in QMMMFB2'
           WRITE(LUPRI,*)
        ENDIF

        CALL GET_FROM_FILE_1('INDUCED_DIPOLES',NNZAL,WORK(KINDMOM))

#if defined(VAR_MPI)
        IF (NODTOT .GE. 1) THEN
           CALL QMMMB2_M3(WORK(KINDMOM),EXPVAL,IPRINT,DENMAT,
     &                    WORK(KLAST),LWRK)
        ELSE
#endif
        IL = 0 ! Counter for induced dipole moments

        DO 300 I = 1,MMCENT

          IF (ZEROAL(I) .EQ. -1) GOTO 300

          CALL DZERO(WORK(KTMP),18)

          KPATOM = 0
          NOSIM = 0
          TOFILE = .FALSE.
          TRIMAT = .FALSE.
          EXP1VL = .TRUE.
          DIPORG(1) = MMCORD(1,I)
          DIPORG(2) = MMCORD(2,I)
          DIPORG(3) = MMCORD(3,I)
          CALL GET1IN(DUMMY,'EFIELB2',NOSIM,WORK(KLAST),
     &                LWRK,LABINT,INTREP,INTADR,I,TOFILE,
     &                KPATOM,TRIMAT,WORK(KTMP),EXP1VL,DENMAT,
     &                IPRINT)

          IF (IPRINT .GE. 4) THEN
             WRITE (LUPRI,'(/A)') ' Efield Ex 2.order matrix:'
             CALL OUTPAK(WORK(KTMP),3,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield Ey 2.order matrix:'
             CALL OUTPAK(WORK(KTMP+6),3,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield Ez 2.order matrix:'
             CALL OUTPAK(WORK(KTMP+12),3,1,LUPRI)
          END IF

          FACX = -1.0D0*WORK(KINDMOM+IL+0)
          FACY = -1.0D0*WORK(KINDMOM+IL+1)
          FACZ = -1.0D0*WORK(KINDMOM+IL+2)

          DO I2 = 1,6
             EXPVAL(I2) = EXPVAL(I2) + WORK(KTMP+I2-1)*FACX +
     &                    WORK(KTMP+I2+5)*FACY + WORK(KTMP+I2+11)*FACZ
          END DO

          IL = IL+3

  300   CONTINUE

#if defined(VAR_MPI)
        ENDIF ! (NODTOT.GE.1) ... ELSE 
#endif
      ENDIF !(IPOLTP .GT. 0)

C     3a. If both permanent dipoles AND polarizabilities construct TOTAL 
C         dipole contribution for all MM centers

      IF (LCOMB)  THEN

        IF (LOCDEB) THEN
          WRITE(LUPRI,*) 'London gradient from total dipoles'
        ENDIF

        KTMP    = 1
        KINDMOM = KTMP    + 18       ! INDMOM contains the ind. dips.
        KTOTMOM = KINDMOM + 3*NNZAL  ! TOTMOM contains the TOTAL dips.
        KLAST   = KTOTMOM + 3*MMCENT
        LWRK    = LWORK   - KLAST +1

        IF (LWRK .LT. 0) CALL ERRWRK('QMMMB2',-KLAST,LWORK)

        CALL DZERO(WORK(KTMP),18)
        CALL DZERO(WORK(KINDMOM),3*NNZAL)
        CALL DZERO(WORK(KTOTMOM),3*MMCENT)

        IF (LOCDEB) THEN
           WRITE(LUPRI,*)
           WRITE(LUPRI,*) 'Ind. dips read from file in QMMMB2'
           WRITE(LUPRI,*)
        ENDIF

        CALL GET_FROM_FILE_1('INDUCED_DIPOLES',NNZAL,WORK(KINDMOM))

        IL = 0
        IJ = 0
        DO 400 I = 1,MMCENT
          IF (ZEROAL(I) .EQ. -1) THEN
            WORK(KTOTMOM+IL+0) = MUL1MM(1,I)
            WORK(KTOTMOM+IL+1) = MUL1MM(2,I)
            WORK(KTOTMOM+IL+2) = MUL1MM(3,I)
          ELSE IF (ZEROAL(I) .EQ. 1) THEN
            WORK(KTOTMOM+IL+0) = MUL1MM(1,I) + WORK(KINDMOM+IJ+0)
            WORK(KTOTMOM+IL+1) = MUL1MM(2,I) + WORK(KINDMOM+IJ+1)
            WORK(KTOTMOM+IL+2) = MUL1MM(3,I) + WORK(KINDMOM+IJ+2)
            IJ=IJ+3
          ENDIF
          IL=IL+3
  400   CONTINUE

#if defined(VAR_MPI)
        IF (NODTOT .GE. 1) THEN
           CALL QMMMB2_M4(WORK(KTOTMOM),EXPVAL,IPRINT,DENMAT,
     &                    WORK(KLAST),LWRK)
        ELSE
#endif
Cahs        IL = 0 ! Counter for dipole moments

        DO 500 I = 1,MMCENT
          IL = 3*(I - 1)
C         See if the total dipole moment at this site is zero
          DNORM2 = WORK(KTOTMOM+IL+0)*WORK(KTOTMOM+IL+0)
     &           + WORK(KTOTMOM+IL+1)*WORK(KTOTMOM+IL+1)
     &           + WORK(KTOTMOM+IL+2)*WORK(KTOTMOM+IL+2)
          DNORM = SQRT(DNORM2)
          IF (ABS(DNORM) .LE. THRMM) GOTO 500

          CALL DZERO(WORK(KTMP),18)

          KPATOM = 0
          NOSIM = 0
          TOFILE = .FALSE.
          TRIMAT = .FALSE.
          EXP1VL = .TRUE.
          DIPORG(1) = MMCORD(1,I)
          DIPORG(2) = MMCORD(2,I)
          DIPORG(3) = MMCORD(3,I)
          CALL GET1IN(DUMMY,'EFIELB2',NOSIM,WORK(KLAST),
     &                LWRK,LABINT,INTREP,INTADR,I,TOFILE,
     &                KPATOM,TRIMAT,WORK(KTMP),EXP1VL,DENMAT,
     &                IPRINT)

          IF (IPRINT .GE. 4) THEN
             WRITE (LUPRI,'(/A)') ' Efield Ex 2.order matrix:'
             CALL OUTPAK(WORK(KTMP),3,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield Ey 2.order matrix:'
             CALL OUTPAK(WORK(KTMP+6),3,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield Ez 2.order matrix:'
             CALL OUTPAK(WORK(KTMP+12),3,1,LUPRI)
          END IF

          FACX = -1.0D0*WORK(KTOTMOM+IL+0)
          FACY = -1.0D0*WORK(KTOTMOM+IL+1)
          FACZ = -1.0D0*WORK(KTOTMOM+IL+2)

          DO I2 = 1,6
             EXPVAL(I2) = EXPVAL(I2) + WORK(KTMP+I2-1)*FACX +
     &                    WORK(KTMP+I2+5)*FACY + WORK(KTMP+I2+11)*FACZ
          END DO

Cahs          IL = IL+3

  500   CONTINUE

#if defined(VAR_MPI)
        ENDIF ! (NODTOT.GE.1) ... ELSE 
#endif
      ENDIF !(LCOMB)

C     Put back the dipole origin

      DIPORG(1) = OBKPX
      DIPORG(2) = OBKPY
      DIPORG(3) = OBKPZ
C
      IF (MMTIME) THEN
        DTIME = SECOND() - DTIME
        TMMB2 = TMMB2 + DTIME
      ENDIF
C
      CALL QEXIT('QMMMB2')
      RETURN
      END
#if defined(VAR_MPI)
C****************************************************************************
C     Parallel section, AHS
C****************************************************************************
C/* deck qmmmfirst_m1 */
      SUBROUTINE QMMMFIRST_M1(VEC2,IPRINT,INTPRI,WORK,LWORK)
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "orgcom.h"
#include "nuclei.h"
#include "qm3.h"
#include "qmmm.h"
#include "inforb.h"
#include "inftap.h"
#include "magone.h"
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"
C defined parallel calculation types  
#include "iprtyp.h"

      DIMENSION WORK(*),VEC2(*)

      CALL QENTER('QMMMFIRST_M1')

      KVEC2 = 1
      JVEC2 = KVEC2 + 3*N2BASX
      KLAST = JVEC2 + 3*N2BASX
      LWRK  = LWORK - KLAST + 1

      IF (LWRK .LT. 0) CALL ERRWRK('QMMMFIRST_M1',-KLAST,LWORK)

      IPRTYP = QMMMFIRST_1_WORK

      CALL MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPRINT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(N2BASX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NNBASX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MMCENT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(IPOLTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MXMMIT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MXMMDI,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPQMMM,1,'INTEGER',MASTER)

      CALL MPIXBCAST(MMCORD,3*MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(MUL0MM,MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(QMCOM,3,'DOUBLE',MASTER)

      CALL MPIXBCAST(ICMO,8,'INTEGER',MASTER)
      CALL MPIXBCAST(FORQM3,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(MMMAT,1,'LOGICAL',MASTER)

      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(CORD,3*MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(CHARGE,MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)

      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NPATOM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPATOM,MXCENT,'INTEGER',MASTER)

      DO I = 1,MMCENT
         IWHO = -1
         CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
         CALL MPIXSEND(I,1,'INTEGER',NWHO,MPTAG2)
      ENDDO

C     Send end message to all slaves
      LEND = -1
      DO ISLAVE = 1, NODTOT
         IWHO = -1
         CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
         CALL MPIXSEND(LEND, 1, 'INTEGER', NWHO, MPTAG2)
      END DO

C     Collect data from all slaves
      CALL DZERO(WORK(KVEC2),3*N2BASX)
      CALL DZERO(WORK(JVEC2),3*N2BASX)

      CALL MPI_REDUCE(WORK(KVEC2),WORK(JVEC2),3*N2BASX,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)

      FAC = 1.0D0
      CALL DAXPY(3*N2BASX,FAC,WORK(JVEC2),1,VEC2(1),1)


      CALL QEXIT('QMMMFIRST_M1')
      RETURN
      END
C****************************************************************************
C/* deck qmmmfirst_s1 */
      SUBROUTINE QMMMFIRST_S1(WORK,LWORK,IPRINT)
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "orgcom.h"
#include "nuclei.h"
#include "qm3.h"
#include "qmmm.h"
#include "inforb.h"
#include "inftap.h"
#include "gnrinf.h"
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"

      DIMENSION WORK(LWORK)
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)

      CALL QENTER('QMMMFIRST_S1')

      QMMM = .TRUE.

      CALL MPIXBCAST(N2BASX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NNBASX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MMCENT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(IPOLTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MXMMIT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MXMMDI,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPQMMM,1,'INTEGER',MASTER)

      CALL MPIXBCAST(MMCORD,3*MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(MUL0MM,MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(QMCOM,3,'DOUBLE',MASTER)

      CALL MPIXBCAST(ICMO,8,'INTEGER',MASTER)
      CALL MPIXBCAST(FORQM3,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(MMMAT,1,'LOGICAL',MASTER)

      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(CORD,3*MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(CHARGE,MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)

      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NPATOM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPATOM,MXCENT,'INTEGER',MASTER)

      KVEC2 = 1
      KLAST = KVEC2 + 3*N2BASX
      LWRK =  LWORK - KLAST + 1

      IF (LWRK .LT. 0) CALL ERRWRK('QMMMFIRST_S1',-KLAST,LWORK)

      CALL DZERO(WORK(KVEC2),3*N2BASX)

 100  CONTINUE

      CALL MPIXSEND(MYNUM, 1, 'INTEGER', MASTER, MPTAG1)
      CALL MPIXRECV(I, 1, 'INTEGER', MASTER, MPTAG2)

      IF (I .GE. 1) THEN

         IF (ABS(MUL0MM(I)) .LE. THRMM) GOTO 100

         KPATOM = 0
         NOSIM = 3
         TOFILE = .FALSE.
         TRIMAT = .FALSE.
         EXP1VL = .FALSE.
         DIPORG(1) = MMCORD(1,I)
         DIPORG(2) = MMCORD(2,I)
         DIPORG(3) = MMCORD(3,I)

         KMAT   = KLAST
         KLAST2 = KMAT + 3*N2BASX
         LWRK2  = LWORK - KLAST2 + 1

         IF (LWRK2 .LT. 0) CALL ERRWRK('QMMMFIRST_S1',-KLAST2,LWORK)

         CALL DZERO(WORK(KMAT),3*NNBASX)

         CALL GET1IN(WORK(KMAT),'PCMBSOL',NOSIM,WORK(KLAST2),
     &                 LWRK2,LABINT,INTREP,INTADR,I,TOFILE,
     &                 KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRINT)

         IF (IPRINT .GE. 12) THEN
            WRITE (LUPRI,'(/A,I3,A)')
     *             ' N(',I,')_ao matrix (Bx): '
            CALL OUTPUT(WORK(KMAT),1,NBAST,1,NBAST,
     &                    NBAST,NBAST,1,LUPRI)
            WRITE (LUPRI,'(/A,I3,A)')
     *             ' N(',I,')_ao matrix (By): '
            CALL OUTPUT(WORK(KMAT+N2BASX),1,NBAST,1,NBAST,
     &                    NBAST,NBAST,1,LUPRI)
            WRITE (LUPRI,'(/A,I3,A)')
     *             ' N(',I,')_ao matrix (Bz): '
            CALL OUTPUT(WORK(KMAT+2*N2BASX),1,NBAST,1,NBAST,
     &                    NBAST,NBAST,1,LUPRI)
         END IF

         FAC1 = -1.0D0*MUL0MM(I)
         CALL DAXPY(3*N2BASX,FAC1,WORK(KMAT),1,WORK(KVEC2),1)
         GO TO 100
      ENDIF

      CALL MPI_REDUCE(WORK(KVEC2),MPI_IN_PLACE,3*N2BASX,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)

      CALL QEXIT('QMMMFIRST_S1')
      RETURN
      END

C/* deck qmmmfirst_m2 */
      SUBROUTINE QMMMFIRST_M2(VEC2,IPRINT,INTPRI,WORK,LWORK)
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "orgcom.h"
#include "nuclei.h"
#include "qm3.h"
#include "qmmm.h"
#include "inforb.h"
#include "inftap.h"
#include "magone.h"
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"
C defined parallel calculation types  
#include "iprtyp.h"

      DIMENSION WORK(*),VEC2(*)

      CALL QENTER('QMMMFIRST_M2')

      KVEC2 = 1
      JVEC2 = KVEC2 + 3*N2BASX
      KLAST = JVEC2 + 3*N2BASX
      LWRK  = LWORK - KLAST + 1

      IF (LWRK .LT. 0) CALL ERRWRK('QMMMFIRST_M2',-KLAST,LWORK)

      IPRTYP = QMMMFIRST_2_WORK

      CALL MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPRINT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(N2BASX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NNBASX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MMCENT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(IPOLTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MXMMIT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MXMMDI,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPQMMM,1,'INTEGER',MASTER)

      CALL MPIXBCAST(MMCORD,3*MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(MUL1MM,3*MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(QMCOM,3,'DOUBLE',MASTER)

      CALL MPIXBCAST(ICMO,8,'INTEGER',MASTER)
      CALL MPIXBCAST(FORQM3,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(MMMAT,1,'LOGICAL',MASTER)

      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(CORD,3*MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(CHARGE,MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)

      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NPATOM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPATOM,MXCENT,'INTEGER',MASTER)

      DO I = 1,MMCENT
         IWHO = -1
         CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
         CALL MPIXSEND(I,1,'INTEGER',NWHO,MPTAG2)
      ENDDO

C     Send end message to all slaves
      LEND = -1
      DO ISLAVE = 1, NODTOT
         IWHO = -1
         CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
         CALL MPIXSEND(LEND, 1, 'INTEGER', NWHO, MPTAG2)
      END DO

C     Collect data from all slaves
      CALL DZERO(WORK(KVEC2),3*N2BASX)
      CALL DZERO(WORK(JVEC2),3*N2BASX)

      CALL MPI_REDUCE(WORK(KVEC2),WORK(JVEC2),3*N2BASX,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)

      FAC = 1.0D0
      CALL DAXPY(3*N2BASX,FAC,WORK(JVEC2),1,VEC2,1)


      CALL QEXIT('QMMMFIRST_M2')
      RETURN
      END
C****************************************************************************
C/* deck qmmmfirst_S2 */
      SUBROUTINE QMMMFIRST_S2(WORK,LWORK,INTPRI)
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "orgcom.h"
#include "nuclei.h"
#include "qm3.h"
#include "qmmm.h"
#include "inforb.h"
#include "inftap.h"
#include "gnrinf.h"
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"

      DIMENSION WORK(LWORK)
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)

      CALL QENTER('QMMMFIRST_S2')

      QMMM = .TRUE.

      CALL MPIXBCAST(N2BASX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NNBASX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MMCENT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(IPOLTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MXMMIT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MXMMDI,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPQMMM,1,'INTEGER',MASTER)

      CALL MPIXBCAST(MMCORD,3*MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(MUL1MM,3*MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(QMCOM,3,'DOUBLE',MASTER)

      CALL MPIXBCAST(ICMO,8,'INTEGER',MASTER)
      CALL MPIXBCAST(FORQM3,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(MMMAT,1,'LOGICAL',MASTER)

      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(CORD,3*MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(CHARGE,MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)

      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NPATOM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPATOM,MXCENT,'INTEGER',MASTER)

      KVEC2 = 1
      KMAT  = KVEC2 + 3*N2BASX
      KLAST = KMAT + 9*N2BASX
      LWRK  =  LWORK - KLAST + 1

      IF (LWRK .LT. 0) CALL ERRWRK('QMMMFIRST_S2',-KLAST,LWORK)

      CALL DZERO(WORK(KVEC2),3*N2BASX)

 100  CONTINUE

      CALL MPIXSEND(MYNUM, 1, 'INTEGER', MASTER, MPTAG1)
      CALL MPIXRECV(I, 1, 'INTEGER', MASTER, MPTAG2)

      IF (I .GE. 1) THEN
          DNORM2 = MUL1MM(1,I)**2+MUL1MM(2,I)**2+MUL1MM(3,I)**2
          DNORM = SQRT(DNORM2)

          IF (ABS(DNORM) .LE. THRMM) GOTO 100

          CALL DZERO(WORK(KMAT),9*NNBASX)

          KPATOM = 0
          NOSIM = 9
          TOFILE = .FALSE.
          TRIMAT = .FALSE.
          EXP1VL = .FALSE.
          DIPORG(1) = MMCORD(1,I)
          DIPORG(2) = MMCORD(2,I)
          DIPORG(3) = MMCORD(3,I)

          CALL GET1IN(WORK(KMAT),'EFIELB1',NOSIM,WORK(KLAST),
     &                LWRK,LABINT,INTREP,INTADR,I,TOFILE,
     &                KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,INTPRI)
          IF (INTPRI .GE. 15) THEN
             WRITE (LUPRI,'(/A)') ' Efield_ao Ex Bx matrix:'
             CALL OUTPUT(WORK(KMAT),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ex By matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ex Bz matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*2),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ey Bx matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*3),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ey By matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*4),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ey Bz matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*5),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ez Bx matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*6),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ez By matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*7),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ez Bz matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*8),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
          END IF

          FACX = -1.0D0*MUL1MM(1,I)
          FACY = -1.0D0*MUL1MM(2,I)
          FACZ = -1.0D0*MUL1MM(3,I)

          CALL DAXPY(N2BASX,FACX,WORK(KMAT),1,WORK(KVEC2),1)
          CALL DAXPY(N2BASX,FACX,WORK(KMAT+N2BASX),1,
     &               WORK(KVEC2+N2BASX),1)
          CALL DAXPY(N2BASX,FACX,WORK(KMAT+2*N2BASX),1,
     &               WORK(KVEC2+2*N2BASX),1)
          CALL DAXPY(N2BASX,FACY,WORK(KMAT+3*N2BASX),1,
     &               WORK(KVEC2),1)
          CALL DAXPY(N2BASX,FACY,WORK(KMAT+4*N2BASX),1,
     &               WORK(KVEC2+N2BASX),1)
          CALL DAXPY(N2BASX,FACY,WORK(KMAT+5*N2BASX),1,
     &               WORK(KVEC2+2*N2BASX),1)
          CALL DAXPY(N2BASX,FACZ,WORK(KMAT+6*N2BASX),1,
     &               WORK(KVEC2),1)
          CALL DAXPY(N2BASX,FACZ,WORK(KMAT+7*N2BASX),1,
     &               WORK(KVEC2+N2BASX),1)
          CALL DAXPY(N2BASX,FACZ,WORK(KMAT+8*N2BASX),1,
     &               WORK(KVEC2+2*N2BASX),1)

         GO TO 100
      ENDIF

      CALL MPI_REDUCE(WORK(KVEC2),MPI_IN_PLACE,3*N2BASX,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)

      CALL QEXIT('QMMMFIRST_S2')
      RETURN
      END

C/* deck qmmmfirst_m3 */
      SUBROUTINE QMMMFIRST_M3(AINDMOM,VEC2,IPRINT,INTPRI,WORK,LWORK)
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "orgcom.h"
#include "nuclei.h"
#include "qm3.h"
#include "qmmm.h"
#include "inforb.h"
#include "inftap.h"
#include "magone.h"
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"
C defined parallel calculation types  
#include "iprtyp.h"

      DIMENSION WORK(*),VEC2(*),AINDMOM(3*NNZAL)

      CALL QENTER('QMMMFIRST_M3')

      KVEC2 = 1
      JVEC2 = KVEC2 + 3*N2BASX
      KLAST = JVEC2 + 3*N2BASX
      LWRK  = LWORK - KLAST + 1

      IF (LWRK .LT. 0) CALL ERRWRK('QMMMFIRST_M3',-KLAST,LWORK)

      IPRTYP = QMMMFIRST_3_WORK

      CALL MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(INTPRI,1,'INTEGER',MASTER)

      CALL MPIXBCAST(N2BASX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NNBASX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MMCENT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(IPOLTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MXMMIT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MXMMDI,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPQMMM,1,'INTEGER',MASTER)

      CALL MPIXBCAST(MMCORD,3*MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(MUL0MM,MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(QMCOM,3,'DOUBLE',MASTER)

      CALL MPIXBCAST(ICMO,8,'INTEGER',MASTER)
      CALL MPIXBCAST(FORQM3,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(MMMAT,1,'LOGICAL',MASTER)

      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NNZAL,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(CORD,3*MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(CHARGE,MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)

      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NPATOM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPATOM,MXCENT,'INTEGER',MASTER)

      CALL MPIXBCAST(AINDMOM,3*NNZAL,'DOUBLE',MASTER)

      IL = 0 ! Counter for induced dipole moments
      DO 300 I = 1,MMCENT
         IF (ZEROAL(I) .EQ. -1) GOTO 300
         IWHO = -1
         CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
         CALL MPIXSEND(I,   1,'INTEGER',NWHO,MPTAG2)
         CALL MPIXSEND(IL,  1,'INTEGER',NWHO,MPTAG2)
         IL = IL+3
 300  CONTINUE

C     Send end message to all slaves
      LEND = -1
      DO ISLAVE = 1, NODTOT
         IWHO = -1
         CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
         CALL MPIXSEND(LEND,1,'INTEGER',NWHO,MPTAG2)
         CALL MPIXSEND(LEND,1,'INTEGER',NWHO,MPTAG2)
      END DO

C     Collect data from all slaves
      CALL DZERO(WORK(KVEC2),3*N2BASX)
      CALL DZERO(WORK(JVEC2),3*N2BASX)

      CALL MPI_REDUCE(WORK(KVEC2),WORK(JVEC2),3*N2BASX,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)

      FAC = 1.0D0
      CALL DAXPY(3*N2BASX,FAC,WORK(JVEC2),1,VEC2,1)


      CALL QEXIT('QMMMFIRST_M3')
      RETURN
      END
C****************************************************************************
C/* deck qmmmfirst_s3 */
      SUBROUTINE QMMMFIRST_S3(WORK,LWORK,INTPRI)
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "orgcom.h"
#include "nuclei.h"
#include "qm3.h"
#include "qmmm.h"
#include "inforb.h"
#include "inftap.h"
#include "gnrinf.h"
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"

      DIMENSION WORK(LWORK)
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)

      CALL QENTER('QMMMFIRST_S3')

      QMMM = .TRUE.

      CALL MPIXBCAST(N2BASX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NNBASX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MMCENT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(IPOLTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MXMMIT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MXMMDI,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPQMMM,1,'INTEGER',MASTER)

      CALL MPIXBCAST(MMCORD,3*MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(MUL0MM,MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(QMCOM,3,'DOUBLE',MASTER)

      CALL MPIXBCAST(ICMO,8,'INTEGER',MASTER)
      CALL MPIXBCAST(FORQM3,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(MMMAT,1,'LOGICAL',MASTER)

      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NNZAL,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(CORD,3*MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(CHARGE,MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)

      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NPATOM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPATOM,MXCENT,'INTEGER',MASTER)

      KVEC2   = 1
      KMAT    = KVEC2 + 3*N2BASX
      KINDMOM = KMAT  + 9*N2BASX
      KLAST   = KINDMOM + 3*NNZAL
      LWRK    = LWORK - KLAST + 1

      IF (LWRK .LT. 0) CALL ERRWRK('QMMMFIRST_S3',-KLAST,LWORK)

      CALL MPIXBCAST(WORK(KINDMOM),3*NNZAL,'DOUBLE',MASTER)

      CALL DZERO(WORK(KVEC2),3*N2BASX)

 100  CONTINUE

      CALL MPIXSEND(MYNUM, 1, 'INTEGER', MASTER, MPTAG1)
      CALL MPIXRECV(I, 1, 'INTEGER', MASTER, MPTAG2)
      CALL MPIXRECV(IL, 1, 'INTEGER', MASTER, MPTAG2)

      IF (I .GE. 1) THEN

          CALL DZERO(WORK(KMAT),9*NNBASX)

          KPATOM = 0
          NOSIM = 9
          TOFILE = .FALSE.
          TRIMAT = .FALSE.
          EXP1VL = .FALSE.
          DIPORG(1) = MMCORD(1,I)
          DIPORG(2) = MMCORD(2,I)
          DIPORG(3) = MMCORD(3,I)

          CALL GET1IN(WORK(KMAT),'EFIELB1',NOSIM,WORK(KLAST),
     &                LWRK,LABINT,INTREP,INTADR,I,TOFILE,
     &                KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,INTPRI)
          IF (INTPRI .GE. 15) THEN
             WRITE (LUPRI,'(/A)') ' Efield_ao Ex Bx matrix:'
             CALL OUTPUT(WORK(KMAT),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ex By matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ex Bz matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*2),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ey Bx matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*3),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ey By matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*4),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ey Bz matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*5),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ez Bx matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*6),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ez By matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*7),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ez Bz matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*8),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
          END IF

          FACX = -1.0D0*WORK(KINDMOM+IL+0)
          FACY = -1.0D0*WORK(KINDMOM+IL+1)
          FACZ = -1.0D0*WORK(KINDMOM+IL+2)

          CALL DAXPY(N2BASX,FACX,WORK(KMAT),1,WORK(KVEC2),1)
          CALL DAXPY(N2BASX,FACX,WORK(KMAT+N2BASX),1,
     &               WORK(KVEC2+N2BASX),1)
          CALL DAXPY(N2BASX,FACX,WORK(KMAT+2*N2BASX),1,
     &               WORK(KVEC2+2*N2BASX),1)
          CALL DAXPY(N2BASX,FACY,WORK(KMAT+3*N2BASX),1,
     &               WORK(KVEC2),1)
          CALL DAXPY(N2BASX,FACY,WORK(KMAT+4*N2BASX),1,
     &               WORK(KVEC2+N2BASX),1)
          CALL DAXPY(N2BASX,FACY,WORK(KMAT+5*N2BASX),1,
     &               WORK(KVEC2+2*N2BASX),1)
          CALL DAXPY(N2BASX,FACZ,WORK(KMAT+6*N2BASX),1,
     &               WORK(KVEC2),1)
          CALL DAXPY(N2BASX,FACZ,WORK(KMAT+7*N2BASX),1,
     &               WORK(KVEC2+N2BASX),1)
          CALL DAXPY(N2BASX,FACZ,WORK(KMAT+8*N2BASX),1,
     &               WORK(KVEC2+2*N2BASX),1)

         GO TO 100
      ENDIF

      CALL MPI_REDUCE(WORK(KVEC2),MPI_IN_PLACE,3*N2BASX,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)

      CALL QEXIT('QMMMFIRST_S3')
      RETURN
      END

C/* deck qmmmfirst_m4 */
      SUBROUTINE QMMMFIRST_M4(TOTMOM,VEC2,IPRINT,INTPRI,WORK,LWORK)
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "orgcom.h"
#include "nuclei.h"
#include "qm3.h"
#include "qmmm.h"
#include "inforb.h"
#include "inftap.h"
#include "magone.h"
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"
C defined parallel calculation types  
#include "iprtyp.h"

      DIMENSION WORK(LWORK),VEC2(3*N2BASX),TOTMOM(3*MMCENT)

      CALL QENTER('QMMMFIRST_M4')

      KVEC2 = 1
      JVEC2 = KVEC2 + 3*N2BASX
      KLAST = JVEC2 + 3*N2BASX
      LWRK  = LWORK - KLAST + 1

      IF (LWRK .LT. 0) CALL ERRWRK('QMMMFIRST_M4',-KLAST,LWORK)

      IPRTYP = QMMMFIRST_4_WORK

      CALL MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPRINT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(N2BASX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NNBASX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MMCENT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(IPOLTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MXMMIT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MXMMDI,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPQMMM,1,'INTEGER',MASTER)

      CALL MPIXBCAST(MMCORD,3*MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(MUL0MM,MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(QMCOM,3,'DOUBLE',MASTER)

      CALL MPIXBCAST(ICMO,8,'INTEGER',MASTER)
      CALL MPIXBCAST(FORQM3,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(MMMAT,1,'LOGICAL',MASTER)

      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(CORD,3*MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(CHARGE,MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)

      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NPATOM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPATOM,MXCENT,'INTEGER',MASTER)

      CALL MPIXBCAST(TOTMOM,3*MMCENT,'DOUBLE',MASTER)

      DO I = 1,MMCENT
         IWHO = -1
         CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
         CALL MPIXSEND(I,1,'INTEGER',NWHO,MPTAG2)
      ENDDO

C     Send end message to all slaves
      LEND = -1
      DO ISLAVE = 1, NODTOT
         IWHO = -1
         CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
         CALL MPIXSEND(LEND, 1, 'INTEGER', NWHO, MPTAG2)
      END DO

C     Collect data from all slaves
      CALL DZERO(WORK(KVEC2),3*N2BASX)
      CALL DZERO(WORK(JVEC2),3*N2BASX)

      CALL MPI_REDUCE(WORK(KVEC2),WORK(JVEC2),3*N2BASX,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)

      FAC = 1.0D0
      CALL DAXPY(3*N2BASX,FAC,WORK(JVEC2),1,VEC2,1)


      CALL QEXIT('QMMMFIRST_M4')
      RETURN
      END
C****************************************************************************
C/* deck qmmmfirst_s4 */
      SUBROUTINE QMMMFIRST_S4(WORK,LWORK,INTPRI)
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "orgcom.h"
C nuclei.h: NCTOT, CORD, CHARGE, NUCIND, NUCDEP, 
#include "nuclei.h"
C qm3.h: QMCOM, FORQM3, ISYTP
#include "qm3.h"
C qmmm.h: MMCENT,IPOLTP,MXMMIT,MXMMDI,IPQMMM,MMCORD, MMMAT
#include "qmmm.h"
C inforb.h: N2BASX, NNBASX, ICMO, NBAST
#include "inforb.h"
C inftap.h: LUPROP
#include "inftap.h"
#include "gnrinf.h"
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
C cbiher.h: NPATOM, IPATOM
#include "cbiher.h"

      DIMENSION WORK(LWORK)
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)

      CALL QENTER('QMMMFIRST_S4')

      QMMM = .TRUE.

      CALL MPIXBCAST(N2BASX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NNBASX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MMCENT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(IPOLTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MXMMIT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MXMMDI,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPQMMM,1,'INTEGER',MASTER)

      CALL MPIXBCAST(MMCORD,3*MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(MUL0MM,MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(QMCOM,3,'DOUBLE',MASTER)

      CALL MPIXBCAST(ICMO,8,'INTEGER',MASTER)
      CALL MPIXBCAST(FORQM3,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(MMMAT,1,'LOGICAL',MASTER)

      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(CORD,3*MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(CHARGE,MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)

      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NPATOM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPATOM,MXCENT,'INTEGER',MASTER)

      KMAT    = 1
      KVEC2   = KMAT + 9*N2BASX
      KTOTMOM = KVEC2 + 3*N2BASX
      KLAST   = KTOTMOM + 3*MMCENT
      LWRK    = LWORK - KLAST + 1

      IF (LWRK .LT. 0) CALL ERRWRK('QMMMFIRST_S4',-KLAST,LWORK)

      CALL MPIXBCAST(WORK(KTOTMOM),3*MMCENT,'DOUBLE',MASTER)

      CALL DZERO(WORK(KVEC2),3*N2BASX)

 100  CONTINUE

      CALL MPIXSEND(MYNUM, 1, 'INTEGER', MASTER, MPTAG1)
      CALL MPIXRECV(I, 1, 'INTEGER', MASTER, MPTAG2)

      IF (I .GE. 1) THEN
          IL = 3*(I - 1)
C         See if the total dipole moment at this site is zero
          DNORM2 = WORK(KTOTMOM+IL+0)*WORK(KTOTMOM+IL+0)
     &           + WORK(KTOTMOM+IL+1)*WORK(KTOTMOM+IL+1)
     &           + WORK(KTOTMOM+IL+2)*WORK(KTOTMOM+IL+2)
          DNORM = SQRT(DNORM2)
          IF (ABS(DNORM) .LE. THRMM) GOTO 100

          CALL DZERO(WORK(KMAT),9*NNBASX)

          KPATOM = 0
          NOSIM = 9
          TOFILE = .FALSE.
          TRIMAT = .FALSE.
          EXP1VL = .FALSE.
          DIPORG(1) = MMCORD(1,I)
          DIPORG(2) = MMCORD(2,I)
          DIPORG(3) = MMCORD(3,I)

          CALL GET1IN(WORK(KMAT),'EFIELB1',NOSIM,WORK(KLAST),
     &                LWRK,LABINT,INTREP,INTADR,I,TOFILE,
     &                KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,INTPRI)

          IF (INTPRI .GE. 15) THEN
             WRITE (LUPRI,'(/A)') ' Efield_ao Ex Bx matrix:'
             CALL OUTPUT(WORK(KMAT),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ex By matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ex Bz matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*2),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ey Bx matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*3),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ey By matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*4),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ey Bz matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*5),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ez Bx matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*6),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ez By matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*7),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ez Bz matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*8),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
          END IF

          FACX = -1.0D0*WORK(KTOTMOM+IL+0)
          FACY = -1.0D0*WORK(KTOTMOM+IL+1)
          FACZ = -1.0D0*WORK(KTOTMOM+IL+2)

          CALL DAXPY(N2BASX,FACX,WORK(KMAT),1,WORK(KVEC2),1)
          CALL DAXPY(N2BASX,FACX,WORK(KMAT+N2BASX),1,
     &               WORK(KVEC2+N2BASX),1)
          CALL DAXPY(N2BASX,FACX,WORK(KMAT+2*N2BASX),1,
     &               WORK(KVEC2+2*N2BASX),1)
          CALL DAXPY(N2BASX,FACY,WORK(KMAT+3*N2BASX),1,
     &               WORK(KVEC2),1)
          CALL DAXPY(N2BASX,FACY,WORK(KMAT+4*N2BASX),1,
     &               WORK(KVEC2+N2BASX),1)
          CALL DAXPY(N2BASX,FACY,WORK(KMAT+5*N2BASX),1,
     &               WORK(KVEC2+2*N2BASX),1)
          CALL DAXPY(N2BASX,FACZ,WORK(KMAT+6*N2BASX),1,
     &               WORK(KVEC2),1)
          CALL DAXPY(N2BASX,FACZ,WORK(KMAT+7*N2BASX),1,
     &               WORK(KVEC2+N2BASX),1)
          CALL DAXPY(N2BASX,FACZ,WORK(KMAT+8*N2BASX),1,
     &               WORK(KVEC2+2*N2BASX),1)

         GO TO 100
      ENDIF

      CALL MPI_REDUCE(WORK(KVEC2),MPI_IN_PLACE,3*N2BASX,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)

      CALL QEXIT('QMMMFIRST_S4')
      RETURN
      END
C****************************************************************************
C/* deck qmmmb2_m1 */
      SUBROUTINE QMMMB2_M1(EXPVAL,IPRINT,DENMAT,WORK,LWORK)
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "orgcom.h"
#include "nuclei.h"
#include "qm3.h"
#include "qmmm.h"
#include "inforb.h"
#include "inftap.h"
#include "magone.h"
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"
C defined parallel calculation types  
#include "iprtyp.h"

      DIMENSION WORK(LWORK),EXPVAL(6), DENMAT(NNBASX)

      CALL QENTER('QMMMB2_M1')

      KEXPVAL = 1
      JEXPVAL = KEXPVAL + 6
      KLAST   = JEXPVAL + 6
      LWRK    = LWORK - KLAST + 1

      IF (LWRK .LT. 0) CALL ERRWRK('QMMMB2_M1',-KLAST,LWORK)

      IPRTYP = QMMMB2_1_WORK

      CALL MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPRINT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(N2BASX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NNBASX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MMCENT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(IPOLTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MXMMIT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MXMMDI,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPQMMM,1,'INTEGER',MASTER)

      CALL MPIXBCAST(MMCORD,3*MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(MUL0MM,MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(QMCOM,3,'DOUBLE',MASTER)

      CALL MPIXBCAST(ICMO,8,'INTEGER',MASTER)
      CALL MPIXBCAST(FORQM3,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(MMMAT,1,'LOGICAL',MASTER)

      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(CORD,3*MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(CHARGE,MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)

      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NPATOM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPATOM,MXCENT,'INTEGER',MASTER)

      CALL MPIXBCAST(DENMAT,NNBASX,'DOUBLE',MASTER)

      DO I = 1,MMCENT
         IWHO = -1
         CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
         CALL MPIXSEND(I,1,'INTEGER',NWHO,MPTAG2)
      ENDDO

C     Send end message to all slaves
      LEND = -1
      DO ISLAVE = 1, NODTOT
         IWHO = -1
         CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
         CALL MPIXSEND(LEND, 1, 'INTEGER', NWHO, MPTAG2)
      END DO

C     Collect data from all slaves
      CALL DZERO(WORK(KEXPVAL),6)
      CALL DZERO(WORK(JEXPVAL),6)

      CALL MPI_REDUCE(WORK(KEXPVAL),WORK(JEXPVAL),6,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)

      FAC = 1.0D0
      CALL DAXPY(6,FAC,WORK(JEXPVAL),1,EXPVAL,1)


      CALL QEXIT('QMMMB2_M1')
      RETURN
      END
C****************************************************************************
C/* deck qmmmb2_s1 */
      SUBROUTINE QMMMB2_S1(WORK,LWORK,IPRINT)
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "orgcom.h"
#include "nuclei.h"
#include "qm3.h"
#include "qmmm.h"
#include "inforb.h"
#include "inftap.h"
#include "gnrinf.h"
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"

      DIMENSION WORK(LWORK)
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)

      CALL QENTER('QMMMB2_S1')

      QMMM = .TRUE.

      CALL MPIXBCAST(N2BASX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NNBASX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MMCENT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(IPOLTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MXMMIT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MXMMDI,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPQMMM,1,'INTEGER',MASTER)

      CALL MPIXBCAST(MMCORD,3*MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(MUL0MM,MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(QMCOM,3,'DOUBLE',MASTER)

      CALL MPIXBCAST(ICMO,8,'INTEGER',MASTER)
      CALL MPIXBCAST(FORQM3,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(MMMAT,1,'LOGICAL',MASTER)

      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(CORD,3*MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(CHARGE,MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)

      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NPATOM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPATOM,MXCENT,'INTEGER',MASTER)

      KEXPVAL = 1
      KTMP    = KEXPVAL + 6
      KDENMAT = KTMP + 6
      KLAST   = KDENMAT + NNBASX
      LWRK    = LWORK - KLAST + 1

      CALL MPIXBCAST(WORK(KDENMAT),NNBASX,'DOUBLE',MASTER)

      IF (LWRK .LT. 0) CALL ERRWRK('QMMMB2_S1',-KLAST,LWORK)

      CALL DZERO(WORK(KEXPVAL),6)

 100  CONTINUE

      CALL MPIXSEND(MYNUM, 1, 'INTEGER', MASTER, MPTAG1)
      CALL MPIXRECV(I, 1, 'INTEGER', MASTER, MPTAG2)

      IF (I .GE. 1) THEN

         IF (ABS(MUL0MM(I)) .LE. THRMM) GOTO 100

         CALL DZERO(WORK(KTMP),6)

         KPATOM = 0
         NOSIM = 0
         TOFILE = .FALSE.
         TRIMAT = .FALSE.
         EXP1VL = .TRUE.
         DIPORG(1) = MMCORD(1,I)
         DIPORG(2) = MMCORD(2,I)
         DIPORG(3) = MMCORD(3,I)

         CALL GET1IN(DUMMY,'PCMB2SL',NOSIM,WORK(KLAST),
     &                LWRK,LABINT,INTREP,INTADR,I,TOFILE,
     &                KPATOM,TRIMAT,WORK(KTMP),EXP1VL,
     &                WORK(KDENMAT),IPRINT)

         IF (IPRINT .GE. 4) THEN
            WRITE (LUPRI,'(/A,I3,A)')
     *       ' N(',I,') QMMM correction matrix '
            CALL OUTPAK(WORK(KTMP),3,1,LUPRI)
         END IF

         FAC1 = -1.0D0*MUL0MM(I)
         CALL DAXPY(6,FAC1,WORK(KTMP),1,WORK(KEXPVAL),1)

         GO TO 100
      ENDIF

      CALL MPI_REDUCE(WORK(KEXPVAL),MPI_IN_PLACE,6,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)

      CALL QEXIT('QMMMB2_S1')
      RETURN
      END
C****************************************************************************
C/* deck qmmmb2_m2 */
      SUBROUTINE QMMMB2_M2(EXPVAL,IPRINT,DENMAT,WORK,LWORK)
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "orgcom.h"
#include "nuclei.h"
#include "qm3.h"
#include "qmmm.h"
#include "inforb.h"
#include "inftap.h"
#include "magone.h"
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"
C defined parallel calculation types  
#include "iprtyp.h"

      DIMENSION WORK(*),EXPVAL(6), DENMAT(NNBASX)

      CALL QENTER('QMMMB2_M2')

      KEXPVAL = 1
      JEXPVAL = KEXPVAL + 6
      KLAST   = JEXPVAL + 6
      LWRK    = LWORK - KLAST + 1

      IF (LWRK .LT. 0) CALL ERRWRK('QMMMB2_M2',-KLAST,LWORK)

      IPRTYP = QMMMB2_2_WORK

      CALL MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPRINT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(N2BASX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NNBASX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MMCENT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(IPOLTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MXMMIT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MXMMDI,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPQMMM,1,'INTEGER',MASTER)

      CALL MPIXBCAST(MMCORD,3*MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(MUL1MM,3*MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(QMCOM,3,'DOUBLE',MASTER)

      CALL MPIXBCAST(ICMO,8,'INTEGER',MASTER)
      CALL MPIXBCAST(FORQM3,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(MMMAT,1,'LOGICAL',MASTER)

      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(CORD,3*MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(CHARGE,MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)

      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NPATOM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPATOM,MXCENT,'INTEGER',MASTER)

      CALL MPIXBCAST(DENMAT,NNBASX,'DOUBLE',MASTER)

      DO I = 1,MMCENT
         IWHO = -1
         CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
         CALL MPIXSEND(I,1,'INTEGER',NWHO,MPTAG2)
      ENDDO

C     Send end message to all slaves
      LEND = -1
      DO ISLAVE = 1, NODTOT
         IWHO = -1
         CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
         CALL MPIXSEND(LEND, 1, 'INTEGER', NWHO, MPTAG2)
      END DO

C     Collect data from all slaves
      CALL DZERO(WORK(KEXPVAL),6)
      CALL DZERO(WORK(JEXPVAL),6)

      CALL MPI_REDUCE(WORK(KEXPVAL),WORK(JEXPVAL),6,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)

      FAC = 1.0D0
      CALL DAXPY(6,FAC,WORK(JEXPVAL),1,EXPVAL,1)

      CALL QEXIT('QMMMB2_M2')
      RETURN
      END
C****************************************************************************
C/* deck qmmmb2_S2 */
      SUBROUTINE QMMMB2_S2(WORK,LWORK,IPRINT)
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "orgcom.h"
#include "nuclei.h"
#include "qm3.h"
#include "qmmm.h"
#include "inforb.h"
#include "inftap.h"
#include "gnrinf.h"
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"

      DIMENSION WORK(LWORK)
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)

      CALL QENTER('QMMMB2_S2')

      QMMM = .TRUE.

      CALL MPIXBCAST(N2BASX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NNBASX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MMCENT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(IPOLTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MXMMIT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MXMMDI,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPQMMM,1,'INTEGER',MASTER)

      CALL MPIXBCAST(MMCORD,3*MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(MUL1MM,3*MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(QMCOM,3,'DOUBLE',MASTER)

      CALL MPIXBCAST(ICMO,8,'INTEGER',MASTER)
      CALL MPIXBCAST(FORQM3,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(MMMAT,1,'LOGICAL',MASTER)

      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(CORD,3*MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(CHARGE,MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)

      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NPATOM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPATOM,MXCENT,'INTEGER',MASTER)

      KEXPVAL = 1
      KTMP    = KEXPVAL + 6
      KDENMAT = KTMP + 18
      KLAST   = KDENMAT + NNBASX
      LWRK    = LWORK - KLAST + 1

      IF (LWRK .LT. 0) CALL ERRWRK('QMMMB2_S2',-KLAST,LWORK)

      CALL MPIXBCAST(WORK(KDENMAT),NNBASX,'DOUBLE',MASTER)

      CALL DZERO(WORK(KEXPVAL),6)

 100  CONTINUE

      CALL MPIXSEND(MYNUM, 1, 'INTEGER', MASTER, MPTAG1)
      CALL MPIXRECV(I, 1, 'INTEGER', MASTER, MPTAG2)

      IF (I .GE. 1) THEN
          DNORM2 = MUL1MM(1,I)**2+MUL1MM(2,I)**2+MUL1MM(3,I)**2
          DNORM = SQRT(DNORM2)

          IF (ABS(DNORM) .LE. THRMM) GOTO 100

          CALL DZERO(WORK(KTMP),18)

          KPATOM = 0
          NOSIM = 0
          TOFILE = .FALSE.
          TRIMAT = .FALSE.
          EXP1VL = .TRUE.
          DIPORG(1) = MMCORD(1,I)
          DIPORG(2) = MMCORD(2,I)
          DIPORG(3) = MMCORD(3,I)
          CALL GET1IN(DUMMY,'EFIELB2',NOSIM,WORK(KLAST),
     &                LWRK,LABINT,INTREP,INTADR,I,TOFILE,
     &                KPATOM,TRIMAT,WORK(KTMP),EXP1VL,
     &                WORK(KDENMAT),IPRINT)

          IF (IPRINT .GE. 4) THEN
             WRITE (LUPRI,'(/A)') ' Efield Ex 2.order matrix:'
             CALL OUTPAK(WORK(KTMP),3,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield Ey 2.order matrix:'
             CALL OUTPAK(WORK(KTMP+6),3,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield Ez 2.order matrix:'
             CALL OUTPAK(WORK(KTMP+12),3,1,LUPRI)
          END IF

          FACX = -1.0D0*MUL1MM(1,I)
          FACY = -1.0D0*MUL1MM(2,I)
          FACZ = -1.0D0*MUL1MM(3,I)

          DO I2 = 1,6
             WORK(KEXPVAL+I2-1) = WORK(KEXPVAL+I2-1) + 
     &                            WORK(KTMP+I2-1)*FACX +
     &                            WORK(KTMP+I2+5)*FACY + 
     &                            WORK(KTMP+I2+11)*FACZ
          END DO

         GO TO 100
      ENDIF

      CALL MPI_REDUCE(WORK(KEXPVAL),MPI_IN_PLACE,6,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)

      CALL QEXIT('QMMMB2_S2')
      RETURN
      END
C****************************************************************************
C/* deck qmmmb2_m3 */
      SUBROUTINE QMMMB2_M3(AINDMOM,EXPVAL,IPRINT,DENMAT,WORK,LWORK)
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "orgcom.h"
#include "nuclei.h"
#include "qm3.h"
#include "qmmm.h"
#include "inforb.h"
#include "inftap.h"
#include "magone.h"
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"
C defined parallel calculation types  
#include "iprtyp.h"

      DIMENSION WORK(*),EXPVAL(6),AINDMOM(3*NNZAL), DENMAT(NNBASX)

      CALL QENTER('QMMMB2_M3')

      KEXPVAL = 1
      JEXPVAL = KEXPVAL + 6
      KLAST   = JEXPVAL + 6
      LWRK    = LWORK - KLAST + 1

      IF (LWRK .LT. 0) CALL ERRWRK('QMMMB2_M3',-KLAST,LWORK)

      IPRTYP = QMMMB2_3_WORK

      CALL MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPRINT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(N2BASX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NNBASX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MMCENT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(IPOLTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MXMMIT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MXMMDI,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPQMMM,1,'INTEGER',MASTER)

      CALL MPIXBCAST(MMCORD,3*MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(MUL0MM,MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(QMCOM,3,'DOUBLE',MASTER)

      CALL MPIXBCAST(ICMO,8,'INTEGER',MASTER)
      CALL MPIXBCAST(FORQM3,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(MMMAT,1,'LOGICAL',MASTER)

      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NNZAL,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(CORD,3*MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(CHARGE,MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)

      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NPATOM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPATOM,MXCENT,'INTEGER',MASTER)

      CALL MPIXBCAST(AINDMOM,3*NNZAL,'DOUBLE',MASTER)
      CALL MPIXBCAST(DENMAT,NNBASX,'DOUBLE',MASTER)

      IL = 0 ! Counter for induced dipole moments
      DO 300 I = 1,MMCENT
         IF (ZEROAL(I) .EQ. -1) GOTO 300
         IWHO = -1
         CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
         CALL MPIXSEND(I,   1,'INTEGER',NWHO,MPTAG2)
         CALL MPIXSEND(IL,  1,'INTEGER',NWHO,MPTAG2)
         IL = IL+3
 300  CONTINUE

C     Send end message to all slaves
      LEND = -1
      DO ISLAVE = 1, NODTOT
         IWHO = -1
         CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
         CALL MPIXSEND(LEND,1,'INTEGER',NWHO,MPTAG2)
         CALL MPIXSEND(LEND,1,'INTEGER',NWHO,MPTAG2)
      END DO

C     Collect data from all slaves
      CALL DZERO(WORK(KEXPVAL),6)
      CALL DZERO(WORK(JEXPVAL),6)

      CALL MPI_REDUCE(WORK(KEXPVAL),WORK(JEXPVAL),6,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)

      FAC = 1.0D0
      CALL DAXPY(6,FAC,WORK(JEXPVAL),1,EXPVAL,1)

      CALL QEXIT('QMMMB2_M3')
      RETURN
      END
C****************************************************************************
C/* deck qmmmb2_s3 */
      SUBROUTINE QMMMB2_S3(WORK,LWORK,IPRINT)
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "orgcom.h"
#include "nuclei.h"
#include "qm3.h"
#include "qmmm.h"
#include "inforb.h"
#include "inftap.h"
#include "gnrinf.h"
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"

      DIMENSION WORK(LWORK)
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)

      CALL QENTER('QMMMB2_S3')

      QMMM = .TRUE.

      CALL MPIXBCAST(N2BASX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NNBASX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MMCENT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(IPOLTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MXMMIT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MXMMDI,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPQMMM,1,'INTEGER',MASTER)

      CALL MPIXBCAST(MMCORD,3*MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(MUL0MM,MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(QMCOM,3,'DOUBLE',MASTER)

      CALL MPIXBCAST(ICMO,8,'INTEGER',MASTER)
      CALL MPIXBCAST(FORQM3,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(MMMAT,1,'LOGICAL',MASTER)

      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NNZAL,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(CORD,3*MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(CHARGE,MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)

      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NPATOM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPATOM,MXCENT,'INTEGER',MASTER)

      KEXPVAL   = 1
      KTMP      = KEXPVAL + 6
      KINDMOM   = KTMP  + 18
      KDENMAT   = KINDMOM + 3*NNZAL
      KLAST     = KDENMAT + NNBASX
      LWRK      = LWORK - KLAST + 1

      IF (LWRK .LT. 0) CALL ERRWRK('QMMMB2_S3',-KLAST,LWORK)

      CALL MPIXBCAST(WORK(KINDMOM),3*NNZAL,'DOUBLE',MASTER)
      CALL MPIXBCAST(WORK(KDENMAT),NNBASX,'DOUBLE',MASTER)

      CALL DZERO(WORK(KEXPVAL),6)

 100  CONTINUE

      CALL MPIXSEND(MYNUM, 1, 'INTEGER', MASTER, MPTAG1)
      CALL MPIXRECV(I, 1, 'INTEGER', MASTER, MPTAG2)
      CALL MPIXRECV(IL, 1, 'INTEGER', MASTER, MPTAG2)

      IF (I .GE. 1) THEN
          CALL DZERO(WORK(KTMP),18)

          KPATOM = 0
          NOSIM = 0
          TOFILE = .FALSE.
          TRIMAT = .FALSE.
          EXP1VL = .TRUE.
          DIPORG(1) = MMCORD(1,I)
          DIPORG(2) = MMCORD(2,I)
          DIPORG(3) = MMCORD(3,I)
          CALL GET1IN(DUMMY,'EFIELB2',NOSIM,WORK(KLAST),
     &                LWRK,LABINT,INTREP,INTADR,I,TOFILE,
     &                KPATOM,TRIMAT,WORK(KTMP),EXP1VL,
     &                WORK(KDENMAT),IPRINT)

          IF (IPRINT .GE. 4) THEN
             WRITE (LUPRI,'(/A)') ' Efield Ex 2.order matrix:'
             CALL OUTPAK(WORK(KTMP),3,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield Ey 2.order matrix:'
             CALL OUTPAK(WORK(KTMP+6),3,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield Ez 2.order matrix:'
             CALL OUTPAK(WORK(KTMP+12),3,1,LUPRI)
          END IF

          FACX = -1.0D0*WORK(KINDMOM+IL+0)
          FACY = -1.0D0*WORK(KINDMOM+IL+1)
          FACZ = -1.0D0*WORK(KINDMOM+IL+2)

          DO I2 = 1,6
             WORK(KEXPVAL+I2-1) = WORK(KEXPVAL+I2-1) + 
     &                            WORK(KTMP+I2-1)*FACX +
     &                            WORK(KTMP+I2+5)*FACY + 
     &                            WORK(KTMP+I2+11)*FACZ
          END DO

         GO TO 100
      ENDIF

      CALL MPI_REDUCE(WORK(KEXPVAL),MPI_IN_PLACE,6,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)

      CALL QEXIT('QMMMB2_S3')
      RETURN
      END
C****************************************************************************
C/* deck qmmmb2_m4 */
      SUBROUTINE QMMMB2_M4(TOTMOM,EXPVAL,IPRINT,DENMAT,WORK,LWORK)
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "orgcom.h"
#include "nuclei.h"
#include "qm3.h"
#include "qmmm.h"
#include "inforb.h"
#include "inftap.h"
#include "magone.h"
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"
C defined parallel calculation types  
#include "iprtyp.h"

      DIMENSION WORK(LWORK),EXPVAL(6),TOTMOM(3*MMCENT), DENMAT(NNBASX)

      CALL QENTER('QMMMB2_M4')

      KEXPVAL = 1
      JEXPVAL = KEXPVAL + 6
      KLAST   = JEXPVAL + 6
      LWRK    = LWORK - KLAST + 1

      IF (LWRK .LT. 0) CALL ERRWRK('QMMMB2_M4',-KLAST,LWORK)

      IPRTYP = QMMMB2_4_WORK

      CALL MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPRINT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(N2BASX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NNBASX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MMCENT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(IPOLTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MXMMIT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MXMMDI,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPQMMM,1,'INTEGER',MASTER)

      CALL MPIXBCAST(MMCORD,3*MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(MUL0MM,MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(QMCOM,3,'DOUBLE',MASTER)

      CALL MPIXBCAST(ICMO,8,'INTEGER',MASTER)
      CALL MPIXBCAST(FORQM3,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(MMMAT,1,'LOGICAL',MASTER)

      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(CORD,3*MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(CHARGE,MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)

      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NPATOM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPATOM,MXCENT,'INTEGER',MASTER)

      CALL MPIXBCAST(TOTMOM,3*MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(DENMAT,NNBASX,'DOUBLE',MASTER)

      DO I = 1,MMCENT
         IWHO = -1
         CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
         CALL MPIXSEND(I,1,'INTEGER',NWHO,MPTAG2)
      ENDDO

C     Send end message to all slaves
      LEND = -1
      DO ISLAVE = 1, NODTOT
         IWHO = -1
         CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
         CALL MPIXSEND(LEND, 1, 'INTEGER', NWHO, MPTAG2)
      END DO

C     Collect data from all slaves
      CALL DZERO(WORK(KEXPVAL),6)
      CALL DZERO(WORK(JEXPVAL),6)

      CALL MPI_REDUCE(WORK(KEXPVAL),WORK(JEXPVAL),6,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)

      FAC = 1.0D0
      CALL DAXPY(6,FAC,WORK(JEXPVAL),1,EXPVAL,1)

      CALL QEXIT('QMMMB2_M4')
      RETURN
      END
C****************************************************************************
C/* deck qmmmb2_s4 */
      SUBROUTINE QMMMB2_S4(WORK,LWORK,IPRINT)
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "orgcom.h"
C nuclei.h: NCTOT, CORD, CHARGE, NUCIND, NUCDEP, 
#include "nuclei.h"
C qm3.h: QMCOM, FORQM3, ISYTP
#include "qm3.h"
C qmmm.h: MMCENT,IPOLTP,MXMMIT,MXMMDI,IPQMMM,MMCORD, MMMAT
#include "qmmm.h"
C inforb.h: N2BASX, NNBASX, ICMO, NBAST
#include "inforb.h"
C inftap.h: LUPROP
#include "inftap.h"
#include "gnrinf.h"
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
C cbiher.h: NPATOM, IPATOM
#include "cbiher.h"

      DIMENSION WORK(LWORK)
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)

      CALL QENTER('QMMMB2_S4')

      QMMM = .TRUE.

      CALL MPIXBCAST(N2BASX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NNBASX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MMCENT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(IPOLTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MXMMIT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MXMMDI,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPQMMM,1,'INTEGER',MASTER)

      CALL MPIXBCAST(MMCORD,3*MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(MUL0MM,MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(QMCOM,3,'DOUBLE',MASTER)

      CALL MPIXBCAST(ICMO,8,'INTEGER',MASTER)
      CALL MPIXBCAST(FORQM3,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(MMMAT,1,'LOGICAL',MASTER)

      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(CORD,3*MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(CHARGE,MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)

      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NPATOM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPATOM,MXCENT,'INTEGER',MASTER)

      KTMP    = 1
      KEXPVAL = KTMP + 18
      KTOTMOM = KEXPVAL + 6
      KDENMAT = KTOTMOM + 3*MMCENT
      KLAST   = KDENMAT + NNBASX
      LWRK    = LWORK - KLAST + 1

      IF (LWRK .LT. 0) CALL ERRWRK('QMMMB2_S4',-KLAST,LWORK)

      CALL MPIXBCAST(WORK(KTOTMOM),3*MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(WORK(KDENMAT),NNBASX,'DOUBLE',MASTER)

      CALL DZERO(WORK(KEXPVAL),6)

 100  CONTINUE

      CALL MPIXSEND(MYNUM, 1, 'INTEGER', MASTER, MPTAG1)
      CALL MPIXRECV(I, 1, 'INTEGER', MASTER, MPTAG2)

      IF (I .GE. 1) THEN
          IL = 3*(I - 1)
C         See if the total dipole moment at this site is zero
          DNORM2 = WORK(KTOTMOM+IL+0)*WORK(KTOTMOM+IL+0)
     &           + WORK(KTOTMOM+IL+1)*WORK(KTOTMOM+IL+1)
     &           + WORK(KTOTMOM+IL+2)*WORK(KTOTMOM+IL+2)
          DNORM = SQRT(DNORM2)
          IF (ABS(DNORM) .LE. THRMM) GOTO 100

          CALL DZERO(WORK(KTMP),18)

          KPATOM = 0
          NOSIM = 0
          TOFILE = .FALSE.
          TRIMAT = .FALSE.
          EXP1VL = .TRUE.
          DIPORG(1) = MMCORD(1,I)
          DIPORG(2) = MMCORD(2,I)
          DIPORG(3) = MMCORD(3,I)
          CALL GET1IN(DUMMY,'EFIELB2',NOSIM,WORK(KLAST),
     &                LWRK,LABINT,INTREP,INTADR,I,TOFILE,
     &                KPATOM,TRIMAT,WORK(KTMP),EXP1VL,
     &                WORK(KDENMAT),IPRINT)

          IF (IPRINT .GE. 4) THEN
             WRITE (LUPRI,'(/A)') ' Efield Ex 2.order matrix:'
             CALL OUTPAK(WORK(KTMP),3,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield Ey 2.order matrix:'
             CALL OUTPAK(WORK(KTMP+6),3,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield Ez 2.order matrix:'
             CALL OUTPAK(WORK(KTMP+12),3,1,LUPRI)
          END IF

          FACX = -1.0D0*WORK(KTOTMOM+IL+0)
          FACY = -1.0D0*WORK(KTOTMOM+IL+1)
          FACZ = -1.0D0*WORK(KTOTMOM+IL+2)

          DO I2 = 1,6
             WORK(KEXPVAL+I2-1) = WORK(KEXPVAL+I2-1) + 
     &                            WORK(KTMP+I2-1)*FACX +
     &                            WORK(KTMP+I2+5)*FACY + 
     &                            WORK(KTMP+I2+11)*FACZ
          END DO

         GO TO 100
      ENDIF

      CALL MPI_REDUCE(WORK(KEXPVAL),MPI_IN_PLACE,6,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)

      CALL QEXIT('QMMMB2_S4')
      RETURN
      END
C****************************************************************************
#endif
