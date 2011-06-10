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
#include "single.h"
C
C  /* Deck pr1drv */
      SUBROUTINE PR1DRV(SOINT,NELMNT,WORK,LWORK,NPQUAD,LABINT,
     &                  INTTYP,INTREP,NOPTYP,NBAST,ANTI,IORDER,DOATOM,
     &                  INTADR,TRIANG,NATOM,SQUARE,IPRINT,DOINT,
     &                  EXPVAL,EXP1EL,DENMAT)
C
C   (I)  DOINT(2,2)     : which integral classes to do
C   (O)  EXPVAL(NOPTYP) : expectation values for each operator component
C   (I)  logical EXP1EL : calculate EXPVAL, using DENMAT
C   (I)  DENMAT(:)      : one-electron density matrix
C
#include "implicit.h"
#include "priunit.h"
C
      LOGICAL ANTI, DOATOM(*), TRIANG, SQUARE, DOINT(2,2), EXP1EL
      CHARACTER LABINT(NOPTYP)*8
      DIMENSION SOINT(NELMNT,NOPTYP), INTREP(NOPTYP),
     &          WORK(LWORK), INTADR(*), EXPVAL(NOPTYP), DENMAT(*)
C
      IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/)')
     &     ' <<<<<<<<<< OUTPUT FROM PR1DRV >>>>>>>>>> '
      IF (INTTYP .EQ. 12 .OR. INTTYP .EQ. 96) THEN
         NPOINT = NPQUAD
         NSHINT = 9*NATOM*NATOM
      ELSE IF (INTTYP .EQ. 52) THEN
         NPOINT = 1
         NSHINT = 63
      ELSE IF (INTTYP .EQ. 53) THEN
         NPOINT = 1
         NSHINT = 36
      ELSE IF (INTTYP .EQ. 54) THEN
         NPOINT = 1
         NSHINT = 60
      ELSE IF (INTTYP .EQ. 58) THEN
         NPOINT = NPQUAD
         NSHINT = 9*9*NATOM
      ELSE IF (INTTYP .EQ. 70) THEN
         NPOINT = 1
         NSHINT = 3*NOPTYP
      ELSE IF (INTTYP .EQ. 102) THEN
         NPOINT = 1
         NSHINT = 18
      ELSE
         NPOINT = 1
         NSHINT = NOPTYP
      END IF
C
      KWEIGH = 1
      KABSCI = KWEIGH + NPOINT
      KLAST  = KABSCI + NPOINT
      IF (KLAST .GT. LWORK) CALL STOPIT('PR1DRV',' ',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
      CALL PR1DR1(SOINT,WORK(KLAST),LWRK,NPOINT,LABINT,
     &            INTTYP,INTREP,NOPTYP,NBAST,NELMNT,ANTI,IORDER,DOATOM,
     &            WORK(KWEIGH),WORK(KABSCI),TRIANG,NATOM,INTADR,NSHINT,
     &            SQUARE,DOINT,EXPVAL,EXP1EL,DENMAT,IPRINT)
      RETURN
      END
C  /* Deck pr1dr1 */
      SUBROUTINE PR1DR1(SOINT,WORK,LWORK,NPOINT,LABINT,INTTYP,
     &                  INTREP,NOPTYP,NBAST,NELMNT,ANTI,IORDER,DOATOM,
     &                  WEIGHT,ABSCIS,TRIANG,NATOM,INTADR,NSHINT,SQUARE,
     &                  DOINT,EXPVAL,EXP1EL,DENMAT,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "pi.h"
#include "maxaqn.h"
      PARAMETER (D1 = 1.00 D00, D4 = 4.D0,D8 = 8.D0,DMTRA = .75D0)
C
      LOGICAL ANTI, DOATOM(*), HER, TRIANG, SQUARE, DONUC1, DOMOM1,
     &        DOINT(2,2), EXP1EL
C
      CHARACTER LABINT(NOPTYP)*8
C
      DIMENSION SOINT(NELMNT,NOPTYP), WORK(LWORK), INTREP(NOPTYP),
     &          WEIGHT(NPOINT), ABSCIS(NPOINT), INTADR(*),
     &          EXPVAL(NOPTYP), DENMAT(NELMNT)
C
#include "onecom.h"
#include "orgcom.h"
#include "lmns.h"
#include "ccom.h"
#include "nuclei.h"
#include "shells.h"
#include "symmet.h"
#include "symind.h"
#include "huckel.h"
C
#include "ibtfun.h"
C
      IF(THRS.EQ.0.D0) THRS = 1.D-15
      TOLS   = THRS*THRS
      TOLOG  = - LOG(TOLS)
C
C   JMAXD   max j in (d/dx)**j
C   JMAXM : max j in ( x  )**j
C   DONUC1: true  - do < (d/dx)**j A | ... | B>
C           false - do <           A | ... | (d/dx)**j B>
C   DOMOM1: true  - do < ( x  )**j A | ... | B>
C           false - do <           A | ... | ( x  )**j B>
C   (comment added Nov. 2002 / HJAaJ)
C
      IF (ABS(INTTYP) .EQ. 1 .OR. INTTYP .EQ. 45 .OR.
     &    INTTYP .EQ. 1001) THEN
         JMAXD  = 0
         JMAXM  = 0
         DONUC1 = .TRUE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ. 2) THEN
         JMAXD  = 0
         JMAXM  = 0
         DONUC1 = .TRUE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ. 3) THEN
         JMAXD  = 1
         JMAXM  = 0
         DONUC1 = .TRUE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ. 4) THEN
         JMAXD  = 0
         JMAXM  = 2
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ. 5) THEN
         JMAXD  = 1
         JMAXM  = 0
         DONUC1 = .TRUE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ. 6) THEN
         JMAXD  = 0
         JMAXM  = 2
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ. 7) THEN
         JMAXD  = 0
         JMAXM  = 2
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ. 8) THEN
         JMAXD  = 0
         JMAXM  = 0
         DONUC1 = .TRUE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ. 9) THEN
         JMAXD  = 0
         JMAXM  = 0
         DONUC1 = .TRUE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.10) THEN
         JMAXD  = 1
         JMAXM  = 0
         DONUC1 = .TRUE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.11) THEN
         JMAXD  = 0
         JMAXM  = 0
         DONUC1 = .TRUE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.12) THEN
         JMAXD  = 1
         JMAXM  = 0
         DONUC1 = .TRUE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.13) THEN
         JMAXD  = 0
         JMAXM  = 0
         DONUC1 = .TRUE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.14) THEN
         JMAXD  = 1
         JMAXM  = 0
         DONUC1 = .TRUE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.15) THEN
         JMAXD  = 0
         JMAXM  = 1
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.16) THEN
         JMAXD  = 0
         JMAXM  = 2
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.17) THEN
         JMAXD  = 1
         JMAXM  = 1
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.18) THEN
         JMAXD  = 1
         JMAXM  = 1
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.19) THEN
         JMAXD  = 2
         JMAXM  = 1
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.20) THEN
         JMAXD  = 2
         JMAXM  = 1
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.21) THEN
         JMAXD  = 2
         JMAXM  = 0
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.22) THEN
         JMAXD = 0
         JMAXM = 2
         DONUC1= .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.23) THEN
         JMAXD = 1
         JMAXM = 1
         DONUC1= .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.24) THEN
         JMAXD = 2
         JMAXM = 2
         DONUC1= .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.25) THEN
         JMAXD = 2
         JMAXM = 2
         DONUC1= .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.26) THEN
         JMAXD = 0
         JMAXM = 1
         DONUC1= .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.27) THEN
         JMAXD = 1
         JMAXM = 1
         DONUC1= .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.28) THEN
         JMAXD = 1
         JMAXM = 1
         DONUC1= .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.29) THEN
         JMAXD = 0
         JMAXM = 0
         DONUC1= .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.30) THEN
         JMAXD = 0
         JMAXM = 0
         DONUC1= .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.31) THEN
         JMAXD = 0
         JMAXM = 0
         DONUC1= .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.32) THEN
         JMAXD = 0
         JMAXM = 1
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.33) THEN
         JMAXD = 0
         JMAXM = 1
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.34) THEN
         JMAXD = 1
         JMAXM = 1
         DONUC1 = .TRUE.
         DOMOM1 = .FALSE.
      ELSE IF (INTTYP .EQ.35) THEN
         JMAXD = 0
         JMAXM = 0
         DONUC1 = .FALSE.
         DOMOM1 = .FALSE.
      ELSE IF (INTTYP .EQ.36) THEN
         JMAXD = 0
         JMAXM = 1
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.37) THEN
         JMAXD = 0
         JMAXM = 2
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.38) THEN
         JMAXD = 0
         JMAXM = 1
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.39) THEN
         JMAXD  = 0
         JMAXM  = 0
         DONUC1 = .TRUE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.40) THEN
         JMAXD = 4
         JMAXM = 0
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.41) THEN
         JMAXD = 0
         JMAXM = 0
         DONUC1 = .TRUE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.42) THEN
         JMAXD = 0
         JMAXM = 2
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.43) THEN
         JMAXD = 0
         JMAXM = 3
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.44) THEN
         JMAXD  = 1
         JMAXM  = 0
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
C     type 45 already defined together with type 1
      ELSE IF (INTTYP .EQ.46) THEN
         JMAXD  = 0
         JMAXM  = 1 + IORDER
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.47) THEN
         JMAXD  = 0
         JMAXM  = 2 + IORDER
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ. 48) THEN
         JMAXD  = 0
         JMAXM  = 1
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.49) THEN
         JMAXD  = 0
         JMAXM  = 1
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.50) THEN
         JMAXD  = 1
         JMAXM  = 0
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.51) THEN
         JMAXD  = 2
         JMAXM  = 1
         DONUC1 = .TRUE.
         DOMOM1 = .FALSE.
      ELSE IF (INTTYP .EQ.52) THEN
         JMAXD  = 1
         JMAXM  = 1
         DONUC1 = .TRUE.
         DOMOM1 = .FALSE.
      ELSE IF (INTTYP .EQ.53) THEN
         JMAXD  = 1
         JMAXM  = 2
         DONUC1 = .TRUE.
         DOMOM1 = .FALSE.
      ELSE IF (INTTYP .EQ. 54) THEN
         JMAXD  = 1
         JMAXM  = 3
         DONUC1 = .TRUE.
         DOMOM1 = .FALSE.
      ELSE IF (INTTYP .EQ.55) THEN
         JMAXD  = 1
         JMAXM  = 1
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.56) THEN
         JMAXD = 0
         JMAXM = 3
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.57) THEN
         JMAXD = 0
         JMAXM = 1
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.58) THEN
         JMAXD = 1
         JMAXM = 0
         DONUC1 = .TRUE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.59) THEN
         JMAXD  = 1
         JMAXM  = 0
         DONUC1 = .TRUE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.61) THEN
          JMAXD = 0
          JMAXM = 1
          DONUC1= .FALSE.
          DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.62) THEN
          JMAXD = 1
          JMAXM = 1
          DONUC1= .FALSE.
          DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.63) THEN
         JMAXD  = 1
         JMAXM  = 1
         DONUC1 = .TRUE.
         DOMOM1 = .FALSE.
      ELSE IF (INTTYP .EQ.64) THEN
         JMAXD  = 0
         JMAXM  = 0
         DONUC1 = .FALSE.
         DOMOM1 = .FALSE.
      ELSE IF (INTTYP .EQ.65) THEN
         JMAXD = 0
         JMAXM = 3
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.66) THEN
         JMAXD  = 0
         JMAXM  = 1
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.67) THEN
         JMAXD  = 0
         JMAXM  = 2
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.70) THEN
         JMAXD  = 3
         JMAXM  = 0
         DONUC1 = .TRUE.
         DOMOM1 = .FALSE.
      ELSE IF (INTTYP .EQ.71) THEN
         JMAXD  = 2
         JMAXM  = 0
         DONUC1 = .TRUE.
         DOMOM1 = .FALSE.
      ELSE IF (INTTYP .EQ.72) THEN
         JMAXD  = 2
         JMAXM  = 0
         DONUC1 = .FALSE.
         DOMOM1 = .FALSE.
      ELSE IF (INTTYP .EQ.73) THEN
         JMAXD  = 2
         JMAXM  = 0
         DONUC1 = .FALSE.
         DOMOM1 = .FALSE.
      ELSE IF (INTTYP .EQ.74) THEN
         JMAXD  = 2
         JMAXM  = 1
         DONUC1 = .FALSE.
         DOMOM1 = .FALSE.
      ELSE IF (INTTYP .EQ.75) THEN
         JMAXD  = 0
         JMAXM  = 0
         DONUC1 = .FALSE.
         DOMOM1 = .FALSE.
      ELSE IF (INTTYP .EQ.80) THEN
         JMAXD = 1
         JMAXM = 0
         DONUC1 = .TRUE.
         DOMOM1 = .TRUE.
cLig >> setting for INTTYP 81 and 82
      ELSE IF (INTTYP .EQ.81) THEN
         JMAXD  = 1
         JMAXM  = 2
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.82) THEN
         JMAXD  = 1
         JMAXM  = 1
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
cLig <<
      ELSE IF (INTTYP .EQ.83) THEN
         JMAXD  = 1
         JMAXM  = 1
         DONUC1 = .TRUE.
         DOMOM1 = .FALSE.
      ELSE IF (INTTYP .EQ.84) THEN
         JMAXD  = 2
         JMAXM  = 1
         DONUC1 = .TRUE.
         DOMOM1 = .FALSE.
      ELSE IF (INTTYP .EQ.91) THEN
         JMAXD  = 3
         JMAXM  = 1
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.92) THEN
         JMAXD  = 3
         JMAXM  = 0
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.93) THEN
         JMAXD = 2
         JMAXM = 1
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.94) THEN
         JMAXD  = 2
         JMAXM  = 0
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.95) THEN
         JMAXD  = 0
         JMAXM  = 0
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.96) THEN
         JMAXD  = 2
         JMAXM  = 1
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.97) THEN
         JMAXD = 2
         JMAXM = 1
         DONUC1= .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.98) THEN
         JMAXD  = 2
         JMAXM  = 0
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.99) THEN
         JMAXD  = 0
         JMAXM  = 1
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.100) THEN
         JMAXD  = 0
         JMAXM  = 2
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.101) THEN
         JMAXD = 1
         JMAXM = 2
         DONUC1= .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.102) THEN
         JMAXD  = 2
         JMAXM  = 1
         DONUC1 = .FALSE.
         DOMOM1 = .TRUE.
      ELSE IF (INTTYP .EQ.103) THEN
         JMAXD  = 2
         JMAXM  = 1
         DONUC1 = .TRUE.
         DOMOM1 = .FALSE.
      ELSE
         WRITE(LUPRI,*) 'PR1DR1 - Unknown INTTYP : ',INTTYP
         CALL QUIT('PR1DR1: Unknown INTTYP')
      END IF
C
      IF (IPRINT .GT. 5) THEN
         WRITE (LUPRI,'(A,I5)') ' NBAST  ', NBAST
         WRITE (LUPRI,'(A,I5)') ' NELMNT ', NELMNT
         WRITE (LUPRI,'(A,I5)') ' INTTYP ', INTTYP
         WRITE (LUPRI,'(A,I5)') ' NOPTYP ', NOPTYP
         WRITE (LUPRI,'(A,L5)') ' ANTI   ', ANTI
         WRITE (LUPRI,'(A,L5)') ' SQUARE ', SQUARE
         WRITE (LUPRI,'(A,(10I5))') ' INTREP ', (INTREP(I),I=1,NOPTYP)
         WRITE (LUPRI,'(A,(10(1X,A)))')' LABINT ',(LABINT(I),I=1,NOPTYP)
         WRITE (LUPRI,'(A,3F20.10)') ' ORIGIN ', (ORIGIN(I),I=1,3)
         WRITE (LUPRI,'(A,3F20.10)') ' GAUGE  ', (GAGORG(I),I=1,3)
         WRITE (LUPRI,'(A,3F20.10)') ' DIPOLE ', (DIPORG(I),I=1,3)
         WRITE (LUPRI,'(A,3F20.10)') ' CAVITY ', (CAVORG(I),I=1,3)
         WRITE (LUPRI,'(A,I5)') ' JMAXD  ',JMAXD
         WRITE (LUPRI,'(A,I5)') ' JMAXM  ',JMAXM
      END IF
C
      IF (.NOT. EXP1EL) CALL DZERO(SOINT,NELMNT*NOPTYP)
C
C     ****************************************************************
C     *** Prepare for Gaussian quadrature (diamagnetic spin-orbit) ***
C     ****************************************************************
C
      IF (INTTYP .EQ. 12 .OR. INTTYP .EQ. 58 .OR. INTTYP .EQ. 96) THEN
C        HER = NPOINT .EQ. 5 .OR. NPOINT .EQ. 6 .OR.
C    &         NPOINT .EQ. 8 .OR. NPOINT .EQ. 10
         HER = .FALSE.
         IF (HER) THEN
            CALL GAUHER(ABSCIS,WEIGHT,NPOINT)
            DO 100 IPOINT = 1, NPOINT
               EXPD = DUMMY
               ABSCIS(IPOINT) = ABSCIS(IPOINT)**2
               WEIGHT(IPOINT) = D4*WEIGHT(IPOINT)*EXPD/SQRTPI
  100       CONTINUE
         ELSE
            CALL GAULEG(-D1,D1,ABSCIS,WEIGHT,NPOINT)
            DO 110 IPOINT = 1, NPOINT
               D1PA   = D1 + ABSCIS(IPOINT)
               D1MA   = D1 - ABSCIS(IPOINT)
               ABSCIS(IPOINT) = (DMTRA*D1PA/D1MA)**2
               WEIGHT(IPOINT) = D8*(DMTRA**3)*WEIGHT(IPOINT)*(D1PA**2)
     &                                           /(SQRTPI*(D1MA**4))
 110        CONTINUE
         END IF
      END IF
C
C     ***********************************************************
C     ***** If expectation values, construct density matrix *****
C     ***********************************************************
C
      IF (EXP1EL) CALL DZERO(EXPVAL,NOPTYP)
C
C     ************************************************************
C     ***** Triangular loop over symmetry independent shells *****
C     ************************************************************
C
      IBCNT = 1
      IF (INTTYP .EQ. 45) IBCNT = 2
      IORBA = 0
      CALL IZERO(ISOFRA, 8)
      KMAXT = KMAX
      IDENA = 0
      DO 200 ISHELA = 1,KMAXT
         NHKTA = NHKT(ISHELA)
         KHKTA = KHKT(ISHELA)
         KCKTA = KCKT(ISHELA)
         ICA   = LCLASS(ISHELA)
         SPHRA = SPHR(ISHELA)
         CALL LMNVAL(NHKTA,KCKTA,LVALUA,MVALUA,NVALUA)
         NCENTA = NCENT(ISHELA)
         ICENTA = NUCNUM(NCENTA,1)
         MULA   = ISTBAO(ISHELA)
         MULTA  = MULT(MULA)
         NUCA   = NUCO(ISHELA)
         NUMCFA = NUMCF(ISHELA)
         JSTA   = JSTRT(ISHELA)
         CORAX  = CENT(ISHELA,1,1)
         CORAY  = CENT(ISHELA,2,1)
         CORAZ  = CENT(ISHELA,3,1)
C
C        Compute symmetry integral pointers for contributions
C        from this block.
C
         IF (.NOT.SQUARE) THEN
            DO 210 I = 1, 8
               ISOFRB(I) = 0
               DO 220 J = 1, MXAQN
                 INDFA(I,J) = -10 000 000
 220           CONTINUE
 210        CONTINUE
            DO 230 NA = 1, KHKTA
               DO 240 IREP = 0, MAXREP
               IF (IBTAND(MULA,IBTXOR(IREP,ISYMAO(NHKTA,NA))).EQ.0) THEN
                  ISOFRA(IREP+1)    = ISOFRA(IREP+1) + 1
                  INDFA (IREP+1,NA) = ISOFRA(IREP+1)
               END IF
 240           CONTINUE
 230        CONTINUE
            IF (IPRINT .GT. 20) THEN
               WRITE(LUPRI,'(A,I4)')' IA offsets for shell ',ISHELA
               DO 250 NA = 1,KHKTA
                  WRITE(LUPRI,'(8(1X,I5))') (INDFA(I,NA),I=1,MAXREP+1)
 250           CONTINUE
            END IF
         END IF
         IORBB = 0
         ISHMXB = ISHELA
         IF (SQUARE) ISHMXB = KMAXT
         IDENB0 = 0
         DO 300 ISHELB = 1, ISHMXB
            NHKTB = NHKT(ISHELB)
            KHKTB = KHKT(ISHELB)
            KCKTB = KCKT(ISHELB)
            ICB   = LCLASS(ISHELB)
            SPHRB = SPHR(ISHELB)
            CALL LMNVAL(NHKTB,KCKTB,LVALUB,MVALUB,NVALUB)
            NCENTB = NCENT(ISHELB)
            MULB   = ISTBAO(ISHELB)
            MULTB  = MULT(MULB)
            NUCB   = NUCO(ISHELB)
            NUMCFB = NUMCF(ISHELB)
            JSTB   = JSTRT(ISHELB)
            CORBX0 = CENT(ISHELB,1,IBCNT)
            CORBY0 = CENT(ISHELB,2,IBCNT)
            CORBZ0 = CENT(ISHELB,3,IBCNT)
            KHKTAB = KHKTA*KHKTB
            KCKTAB = KCKTA*KCKTB
            MAB    = IBTOR(MULA,MULB)
            KAB    = IBTAND(MULA,MULB)
            HKAB   = FMULT(KAB)
C
            SPHRAB = SPHRA .OR. SPHRB
C
C           Compute symmetry integral pointers.
C
            IF (.NOT.SQUARE) THEN
               DO 310 I = 1, 8
                  DO 320 J = 1, MXAQN
                    INDFB(I,J) = -10 000 000
 320              CONTINUE
 310           CONTINUE
               DO 330 NB = 1, KHKTB
               DO 340 IREP = 0, MAXREP
               IF (IBTAND(MULB,IBTXOR(IREP,ISYMAO(NHKTB,NB))).EQ.0) THEN
                  ISOFRB(IREP+1)    = ISOFRB(IREP+1) + 1
                  INDFB (IREP+1,NB) = ISOFRB(IREP+1)
               END IF
 340           CONTINUE
 330           CONTINUE
               IF (IPRINT .GT. 20) THEN
                  WRITE(LUPRI,'(A,I4)')' IB offsets for shell ',ISHELB
                  DO 350 NB = 1, KHKTB
                     WRITE(LUPRI,'(8(1X,I5))')(INDFB(I,NB),I=1,MAXREP+1)
 350              CONTINUE
               END IF
            END IF
            IF(.NOT.DOINT(ICA,ICB)) GOTO 300
            IF (IPRINT .GT. 05) WRITE (LUPRI, 1000) ISHELA, ISHELB
            IF (IPRINT .GE. 10) THEN
                WRITE (LUPRI,'(A,2I10)') ' NHKT   ', NHKTA, NHKTB
                WRITE (LUPRI,'(A,2I10)') ' KHKT   ', KHKTA, KHKTB
                WRITE (LUPRI,'(A,2I10)') ' KCKT   ', KCKTA, KCKTB
                WRITE (LUPRI,'(A,2I10)') ' NCENT  ', NCENTA, NCENTB
                WRITE (LUPRI,'(A,2I10)') ' ISTBAO ', MULA, MULB
                WRITE (LUPRI,'(A,2I10)') ' MULT   ', MULTA, MULTB
                WRITE (LUPRI,'(A,2I10)') ' NUC    ', NUCA, NUCB
                WRITE (LUPRI,'(A,2I10)') ' NUMCF  ', NUMCFA, NUMCFB
                WRITE (LUPRI,'(A,2I10)') ' JST    ', JSTA, JSTB
                WRITE (LUPRI,'(A,2F12.6)') ' CORAX    ', CORAX, CORBX0
                WRITE (LUPRI,'(A,2F12.6)') ' CORAY    ', CORAY, CORBY0
                WRITE (LUPRI,'(A,2F12.6)') ' CORAZ    ', CORAZ, CORBZ0
            END IF
C
C           ******************************
C           ***** Symmetry integrals *****
C           ******************************
C
            CALL PR1SOP(SOINT,DENMAT,EXPVAL,WORK,LWORK,
     &                  NPOINT,LABINT,INTTYP,INTREP,NOPTYP,NBAST,NELMNT,
     &                  ANTI,IORDER,DOATOM,WEIGHT,ABSCIS,TRIANG,NATOM,
     &                  INTADR,NSHINT,SQUARE,IPRINT,JMAXD,JMAXM,ISHELA,
     &                  ISHELB,IORBA,IORBB,KAB,CORBX0,CORBY0,CORBZ0,
     &                  TOLOG,TOLS,DONUC1,DOMOM1,MULTA,MULTB,EXP1EL,
     &                  IDENB0)
C
         IDENB0= IDENB0 + KHKTB*MULTB
  300    IORBB = IORBB  + KHKTB
         IORBA = IORBA  + KHKTA
         IDENA = IDENA  + KHKTA*MULTA
  200 CONTINUE
C
C     ****************************
C     ***** Unpack integrals *****
C     ****************************
C
      IF (.NOT.SQUARE .AND. .NOT. EXP1EL) THEN
         IF (NELMNT.GT.LWORK)CALL STOPIT('PR1DR1','SYMUPK',NELMNT,LWORK)
         DO 800 I = 1, NOPTYP
            CALL DCOPY(NELMNT,SOINT(1,I),1,WORK,1)
            CALL SYMUPK(WORK,SOINT(1,I),INTREP(I) + 1,NELMNT)
  800    CONTINUE
      END IF
      RETURN
 1000 FORMAT (//,2X,'***************************************',
     *         /,2X,'********** ISHELA/B =',I3,',',I3,' **********',
     *         /,2X,'***************************************',/)
      END
C  /* Deck pr1sop */
      SUBROUTINE PR1SOP(SOINT,DENMAT,EXPVAL,WORK,LWORK,NPOINT,LABINT,
     &                  INTTYP,INTREP,NOPTYP,NBAST,NELMNT,ANTI,IORDER,
     &                  DOATOM,WEIGHT,ABSCIS,TRIANG,NATOM,INTADR,NSHINT,
     &                  SQUARE,IPRINT,JMAXD,JMAXM,ISHELA,ISHELB,IORBA,
     &                  IORBB,KAB,CORBX0,CORBY0,CORBZ0,TOLOG,TOLS,
     &                  DONUC1,DOMOM1,MULTA,MULTB,EXP1EL,IDENB0)
!
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "mxcent.h"
#include "qm3.h"
#include "nuclei.h"
      LOGICAL ANTI, DOATOM(*), TRIANG, SQUARE, DONUC1, DOMOM1,
     &        EXP1EL
      CHARACTER LABINT(NOPTYP)*8
      DIMENSION SOINT(NELMNT,NOPTYP), WORK(LWORK), INTREP(NOPTYP),
     &          WEIGHT(NPOINT), ABSCIS(NPOINT), INTADR(*),
     &          DENMAT(NELMNT), EXPVAL(NOPTYP)
#include "onecom.h"
C
      INATM = NUCDEP
      IF (FORQM3 .AND. (INTTYP .EQ. 35 .OR. INTTYP .EQ. 29
     &               .OR. INTTYP .EQ. 30)) INATM = NCTOT
      KSHINT = 1
      KSYMC  = KSHINT + KCKTAB*NSHINT
      KCENTC = KSYMC  + (INATM + 1)/IRAT
      KNCNTC = KCENTC + (INATM + 1)/IRAT
      KFACTR = KNCNTC + (INATM + 1)/IRAT
      KCORCX = KFACTR + INATM
      KCORCY = KCORCX + INATM
      KCORCZ = KCORCY + INATM
      KGEXP  = KCORCZ + INATM
      KLAST  = KGEXP  + INATM
      IF (KLAST .GT. LWORK) CALL STOPIT('PR1SOP',' ',KLAST,LWORK)
      LWRK = LWORK - KLAST + 1
      CALL PR1SO1(SOINT,WORK(KSHINT),DENMAT,EXPVAL,WORK(KLAST),LWRK,
     &            NPOINT,LABINT,INTTYP,INTREP,NOPTYP,NBAST,NELMNT,ANTI,
     &            IORDER,DOATOM,WEIGHT,ABSCIS,TRIANG,NATOM,INTADR,
     &            NSHINT,SQUARE,IPRINT,WORK(KFACTR),JMAXD,JMAXM,ISHELA,
     &            ISHELB,IORBA,IORBB,KAB,CORBX0,CORBY0,CORBZ0,
     &            WORK(KCORCX),WORK(KCORCY),WORK(KCORCZ),WORK(KCENTC),
     &            WORK(KSYMC),WORK(KGEXP),TOLOG,TOLS,DONUC1,DOMOM1,
     &            MULTA,MULTB,WORK(KNCNTC),EXP1EL,IDENB0)
      RETURN
      END
C  /* Deck pr1so1 */
      SUBROUTINE PR1SO1(SOINT,SHLINT,DENMAT,EXPVAL,WORK,LWORK,NPOINT,
     &                  LABINT,INTTYP,INTREP,NOPTYP,NBAST,NELMNT,ANTI,
     &                  IORDER,DOATOM,WEIGHT,ABSCIS,TRIANG,NATOM,INTADR,
     &                  NSHINT,SQUARE,IPRINT,FACINT,JMAXD,JMAXM,ISHELA,
     &                  ISHELB,IORBA,IORBB,KAB,CORBX0,CORBY0,CORBZ0,
     &                  CORCX,CORCY,CORCZ,JCENTC,JSYMC,GEXP,TOLOG,TOLS,
     &                  DONUC1,DOMOM1,MULTA,MULTB,NCENTC,EXP1EL,IDENB0)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
C
      LOGICAL   FULMAT, ANTI, DOATOM(*), TRIANG, SQUARE, DONUC1,
     &          MULCHA, DOMOM1, HUCORB, EXP1EL
      CHARACTER LABINT(NOPTYP)*8
      DIMENSION SOINT(NELMNT,NOPTYP), SHLINT(KCKTAB,NSHINT),
     &          WORK(LWORK), INTREP(NOPTYP),
     &          WEIGHT(NPOINT), ABSCIS(NPOINT), FACINT(*),
     &          CORCX(*), CORCY(*), CORCZ(*), NCENTC(*),
     &          INTADR(*), JSYMC(*), JCENTC(*), GEXP(*),
     &          DENMAT(NELMNT),EXPVAL(NOPTYP)
C
      PARAMETER (D0 = 0.0D0)
C
#include "nuclei.h"
#include "onecom.h"
#include "orgcom.h"
#include "shells.h"
#include "huckel.h"
#include "symmet.h"
#include "qm3.h"
#include "gnrinf.h"
#include "elweak.h"
C
#include "ibtfun.h"
      XYZ(I,J) = PT(IBTAND(ISYMAX(I,1),J))
      ITRI(I,J) = MAX(I,J)*(MAX(I,J) - 1)/2 + MIN(I,J)
C
C     Initialization for Coulomb integrals
C
      NATOMC = 0
cLig <> added JMAX adjustment  for INTTYP = 82
      IF ((INTTYP .EQ. 5) .OR. (INTTYP .EQ. 10)
     &                    .OR. (INTTYP .EQ. 11)
     &                    .OR. (INTTYP .EQ. 12)
     &                    .OR. (INTTYP .EQ. 13)
     &                    .OR. (INTTYP .EQ. 19)
     &                    .OR. (INTTYP .EQ. 20)
     &                    .OR. (INTTYP .EQ. 24)
     &                    .OR. (INTTYP .EQ. 25)
     &                    .OR. (INTTYP .EQ. 26)
     &                    .OR. (INTTYP .EQ. 27)
     &                    .OR. (INTTYP .EQ. 28)
     &                    .OR. (INTTYP .EQ. 29)
     &                    .OR. (INTTYP .EQ. 30)
     &                    .OR. (INTTYP .EQ. 31)
     &                    .OR. (INTTYP .EQ. 35)
     &                    .OR. (INTTYP .EQ. 38)
     &                    .OR. (INTTYP .EQ. 48)
     &                    .OR. (INTTYP .EQ. 51)
     &                    .OR. (INTTYP .EQ. 57)
     &                    .OR. (INTTYP .EQ. 58)
     &                    .OR. (INTTYP .EQ. 61)
     &                    .OR. (INTTYP .EQ. 62)
     &                    .OR. (INTTYP .EQ. 63)
     &                    .OR. (INTTYP .EQ. 64)
     &                    .OR. (INTTYP .EQ. 66)
     &                    .OR. (INTTYP .EQ. 67)
     &                    .OR. (INTTYP .EQ. 70)
     &                    .OR. (INTTYP .EQ. 72)
     &                    .OR. (INTTYP .EQ. 73)
     &                    .OR. (INTTYP .EQ. 74)
     &                    .OR. (INTTYP .EQ. 75)
     &                    .OR. (INTTYP .EQ. 80)
     &                    .OR. (INTTYP .EQ. 82)
     &                    .OR. (INTTYP .EQ. 92)
     &                    .OR. (INTTYP .EQ. 93)
     &                    .OR. (INTTYP .EQ. 94)
     &                    .OR. (INTTYP .EQ. 96)
     &                    .OR. (INTTYP .EQ. 97)
     &                    .OR. (INTTYP .EQ. 99)
     &                    .OR. (INTTYP .EQ. 100)
     &   ) THEN

         JMAX = NHKTA + NHKTB - 2
         IF (INTTYP .EQ.  5) JMAX = JMAX + 1
         IF (INTTYP .EQ. 10) JMAX = JMAX + 1
         IF (INTTYP .EQ. 11) JMAX = JMAX + 2
         IF (INTTYP .EQ. 12) JMAX = JMAX + 2
         IF (INTTYP .EQ. 13) JMAX = JMAX + 2
         IF (INTTYP .EQ. 19) JMAX = JMAX + 1
         IF (INTTYP .EQ. 20) JMAX = JMAX + 1
         IF (INTTYP .EQ. 24) JMAX = JMAX + 2
         IF (INTTYP .EQ. 25) JMAX = JMAX + 2
         IF (INTTYP .EQ. 26) JMAX = JMAX + 2
         IF (INTTYP .EQ. 27) JMAX = JMAX + 2
         IF (INTTYP .EQ. 28) JMAX = JMAX + 2
         IF (INTTYP .EQ. 29) JMAX = JMAX + 1
         IF (INTTYP .EQ. 30) JMAX = JMAX + 2
         IF (INTTYP .EQ. 31) JMAX = JMAX + 2
         IF (INTTYP .EQ. 38) JMAX = JMAX + 2
         IF (INTTYP .EQ. 48) JMAX = JMAX + 2
         IF (INTTYP .EQ. 51) JMAX = JMAX + 1
         IF (INTTYP .EQ. 57) JMAX = JMAX + 2
         IF (INTTYP .EQ. 58) JMAX = JMAX + 2
         IF (INTTYP .EQ. 58) JMAX = JMAX + 1
         IF (INTTYP .EQ. 61) JMAX = JMAX + 2
         IF (INTTYP .EQ. 62) JMAX = JMAX + 2
         IF (INTTYP .EQ. 63) JMAX = JMAX + 2
         IF (INTTYP .EQ. 66) JMAX = JMAX + 1
         IF (INTTYP .EQ. 67) JMAX = JMAX + 2
         IF (INTTYP .EQ. 70) JMAX = JMAX + 1
         IF (INTTYP .EQ. 72) JMAX = JMAX + 2
         IF (INTTYP .EQ. 73) JMAX = JMAX + 2
         IF (INTTYP .EQ. 74) JMAX = JMAX + 2
         IF (INTTYP .EQ. 75) JMAX = JMAX + 1
         IF (INTTYP .EQ. 80) JMAX = JMAX + 3
cLig <> added JMAX = JMAX+2 for INTTYP = 82
         IF (INTTYP .EQ. 82) JMAX = JMAX + 2
         IF (INTTYP .EQ. 92) JMAX = JMAX + 6
         IF (INTTYP .EQ. 93) JMAX = JMAX + 4
         IF (INTTYP .EQ. 94) JMAX = JMAX + 4
         IF (INTTYP .EQ. 96) JMAX = JMAX + 4
         IF (INTTYP .EQ. 97) JMAX = JMAX + 6
         IF (INTTYP .EQ. 99) JMAX = JMAX + 2
         IF (INTTYP .EQ. 100) JMAX = JMAX + 3
         MULCHA = (INTTYP .EQ. 5) .OR. (INTTYP .EQ. 19)
     &                            .OR. (INTTYP .EQ. 20)
     &                            .OR. (INTTYP .EQ. 24)
     &                            .OR. (INTTYP .EQ. 25)
     &                            .OR. (INTTYP .EQ. 57)
     &                            .OR. (INTTYP .EQ. 61)
     &                            .OR. (INTTYP .EQ. 62)
     &                            .OR. (INTTYP .EQ. 63)
     &                            .OR. (INTTYP .EQ. 64)
     &                            .OR. (INTTYP .EQ. 70)
     &                            .OR. (INTTYP .EQ. 72)
     &                            .OR. (INTTYP .EQ. 73)
         ISTEPU = JMAX + 1
         ISTEPV = ISTEPU*ISTEPU
         NAHGTF = ISTEPU*ISTEPV
C
         IF (      INTTYP .EQ. 75 .OR. INTTYP .EQ. 66
     &        .OR. INTTYP .EQ. 67 .OR. INTTYP .EQ. 99
     &        .OR. INTTYP .EQ. 100
     &        .OR. ( (INTTYP .EQ. 29) .AND. RUNQM3 )
     &        .OR. ( (INTTYP .EQ. 30) .AND. RUNQM3 ) ) THEN
            NATOMC = 0
C
            DO 301 IATOMC = 1, NPOINT
               MULC = MAXREP + 1
               MABC = IBTOR(MULC,KAB)
               CORCX0 = DIPORG(1)
               CORCY0 = DIPORG(2)
               CORCZ0 = DIPORG(3)
               FACTOR = - FMULT(IBTAND(MULC,KAB))
               DO 311 ISYMOP = 0, MAXOPR
                  IF (IBTAND(ISYMOP,MABC) .EQ. 0) THEN
                     NATOMC = NATOMC + 1
                     JSYMC(NATOMC)  = ISYMOP
                     JCENTC(NATOMC) = IATOMC
                     CORCX(NATOMC)  = XYZ(1,ISYMOP)*CORCX0
                     CORCY(NATOMC)  = XYZ(2,ISYMOP)*CORCY0
                     CORCZ(NATOMC)  = XYZ(3,ISYMOP)*CORCZ0
                     FACINT(NATOMC) = FACTOR
                  END IF
 311           CONTINUE
 301        CONTINUE
         ELSE IF (INTTYP.NE.12 .AND. INTTYP.NE.58
     &      .AND. INTTYP.NE.96) THEN
            NATOMC = 0
            IMXATM = NUCIND
            IF (FORQM3 .AND.
     &         (INTTYP .EQ. 35 .OR. INTTYP .EQ. 29 .OR. INTTYP .EQ. 30))
     &         IMXATM = NCTOT
            DO 300 IATOMC = 1, IMXATM
            IF (DOATOM(IATOMC)) THEN
               MULC   = ISTBNU(IATOMC)
               MABC   = IBTOR(MULC,KAB)
               CORCX0 =  CORD(1,IATOMC)
               CORCY0 =  CORD(2,IATOMC)
               CORCZ0 =  CORD(3,IATOMC)
               FACTOR = - FMULT(IBTAND(MULC,KAB))
               IF (MULCHA) FACTOR = FACTOR*CHARGE(IATOMC)/HKAB
               IF ( FORQM3 .AND. (INTTYP .EQ. 35) ) THEN
                 IF (ISUBSY(IATOMC) .NE. 0) THEN
                   IF ( OLDTG ) THEN
                     FACTOR = - FACTOR*CHAOLD(IATOMC)/HKAB
                   ELSE
                     FACTOR = - FACTOR*CHARGE(IATOMC)/HKAB
                   END IF
                 END IF
               END IF
               IF (INTTYP .EQ. 80) THEN
                  FACTOR = FACTOR/HKAB*(
     &            DISOTP(IZATOM(IATOMC),ISOTOP(IATOMC),'NEUTRONS')
     &            -BGWEIN*CHARGE(IATOMC))
               END IF
               DO 310 ISYMOP = 0, MAXOPR
                  IF (IBTAND(ISYMOP,MABC) .EQ. 0) THEN
                     NATOMC = NATOMC + 1
                     JSYMC(NATOMC)  = ISYMOP
                     JCENTC(NATOMC) = IATOMC
                     CORCX(NATOMC)  = XYZ(1,ISYMOP)*CORCX0
                     CORCY(NATOMC)  = XYZ(2,ISYMOP)*CORCY0
                     CORCZ(NATOMC)  = XYZ(3,ISYMOP)*CORCZ0
                     GEXP(NATOMC)   = GNUEXP(IATOMC)
                     FACINT(NATOMC) = FACTOR
                     NCENTC(NATOMC) = NUCNUM(IATOMC,ISYMOP+1)
                  END IF
  310          CONTINUE
            END IF
  300       CONTINUE
         ELSE
C        ... here if INTTYP.EQ.12 .OR. INTTYP.EQ.58 .OR. INTTYP.EQ.96
            NATOMC = NUCDEP
C           ... hjaaj Oct07: must be sure allocation for AHGTF
C               in PR1PRM is big enough for PR1PRD
         END IF
      END IF
C
C     *****************************************
C     ***** Loop over symmetry operations *****
C     *****************************************
C
      IDENB = IDENB0 - KHKTB
      DO 100 ISYMOP = 0, MAXOPR
      IF(IBTAND(ISYMOP,MAB) .EQ. 0) THEN
         LDIAG  = (ISHELA .EQ. ISHELB) .AND. (.NOT.SQUARE)
         ONECEN = .FALSE.
         IDENB  = IDENB + KHKTB
         CORBX  = XYZ(1,ISYMOP)*CORBX0
         CORBY  = XYZ(2,ISYMOP)*CORBY0
         CORBZ  = XYZ(3,ISYMOP)*CORBZ0
         ICENTB = NUCNUM(NCENTB,ISYMOP+1)
         DIFABX = CORAX - CORBX
         DIFABY = CORAY - CORBY
         DIFABZ = CORAZ - CORBZ
         DISTAB = DIFABX*DIFABX + DIFABY*DIFABY + DIFABZ*DIFABZ
C
         IF (IPRINT .GT. 05) THEN
            WRITE (LUPRI, 1000) ISYMOP
            IF (IPRINT .GE. 10) WRITE (LUPRI,'(A,3F12.6)')
     &         ' CORB[XYZ] ', CORBX,CORBY,CORBZ
         END IF
C
C        Huckel / EWMO information matrix
C        ...  INTTYP .eq. -1 signals HUCKEL integrals
C
         IF (INTTYP .EQ. -1) THEN
            SCAL   = 1.0D0
            HUCFAC = 1.0D0
            IF (.NOT. EWMO) THEN
C     Scale overlap integrals for the Huckel matrix part
C hj-aug99: basis functions are normalized to MULT[AB],
C           and HUCMAT was scaled with the 1/sqrt(SCAL) factor.
C           Code modified so projection vectors are also properly scaled.
C
               IF (LCLASS(ISHELA) .EQ. 2) SCAL = DFLOAT(MULTA)
               IF (LCLASS(ISHELB) .EQ. 2) SCAL = DFLOAT(MULTB)*SCAL
            END IF
            IF (LCLASS(ISHELA) .EQ. 2 .AND. LCLASS(ISHELB) .EQ. 2) THEN
               IHUCA = IHUCPT(ISHELA)
               IHUCB = IHUCPT(ISHELB)
               IF (IHUCA .EQ. 0 .OR. IHUCB .EQ. 0) THEN
                  CALL QUIT('Program error IHUCPT(ISHELx) .eq. 0')
               ELSE IF (ISHELA .EQ. ISHELB .AND. DISTAB .LT. 1.D-3) THEN
                  HUCFAC = HUCEXC(IHUCA)
               ELSE IF (DISTAB .LT. 1.0D-3 .OR. EWMO) THEN
C                 ... no intercenter terms in Hamiltonian matrix in EWMO
                  HUCFAC = D0
               ELSE
                  HUCFAC = 0.50D0*HUCCNT*(HUCEXC(IHUCA) + HUCEXC(IHUCB))
               END IF
            END IF
            HUCFAC = HUCFAC/SQRT(SCAL)
         END IF
C
C        *****************************************
C        ***** Calculate Cartesian integrals *****
C        *****************************************
C
         CALL DZERO(SHLINT,KCKTAB*NSHINT)
         CALL PR1PRM(SHLINT,WORK,LWORK,LABINT,INTTYP,
     &               NOPTYP,NBAST,IORDER,DOATOM,IPRINT,TOLS,TOLOG,
     &               NATOMC,DISTAB,FACINT,CORCX,CORCY,CORCZ,JMAXD,
     &               JMAXM,DIFABX,DIFABY,DIFABZ,NPOINT,WEIGHT,ABSCIS,
     &               GEXP,TRIANG,KAB,DONUC1,DOMOM1,HUCFAC,NCENTC,NSHINT)
C
C        ****************************************
C        ***** Transform to spherical basis *****
C        ****************************************
C
         IF (SPHRAB) THEN
            CALL SPHRM1(SHLINT,SHLINT,NSHINT,WORK,LWORK,.FALSE.,IPRINT)
         END IF
C
C        *******************************************
C        ***** Expectation values of integrals *****
C        *******************************************
C
         ONECEN = ICENTA .EQ. ICENTB
         IF (EXP1EL) THEN
            CALL DZERO(WORK,NOPTYP)
            IF ((INTTYP .EQ. 26) .OR. (INTTYP .EQ. 27) .OR.
     &          (INTTYP .EQ. 28) .OR. (INTTYP .EQ. 38) .OR.
     &          (INTTYP .EQ. 48) .OR. (INTTYP .EQ. 74) .OR.
     &          (INTTYP .EQ. 82) .OR. (INTTYP .EQ. 93) .OR.
     &          (INTTYP .EQ. 97)) THEN
               CALL NSTTRA(SHLINT,DUMMY,FACINT,JSYMC,JCENTC,NATOMC,
     &                  ISYMOP,NELMNT,NOPTYP,ANTI,IORBA,IORBB,MULTA,
     &                  MULTB,NBAST,INTTYP,IPRINT,INTADR,EXPVAL,DENMAT,
     &                  EXP1EL)
            ELSE IF (INTTYP .EQ. 29) THEN
               CALL EFTRA(SHLINT,SOINT,FACINT,JSYMC,JCENTC,NATOMC,
     &                    ISYMOP,NELMNT,NOPTYP,ANTI,INTADR,IPRINT,
     &                    MULTA,MULTB,EXPVAL,DENMAT,EXP1EL)
            ELSE IF (INTTYP .EQ. 75) THEN
               CALL TESTRA(SHLINT,EXPVAL,FACINT,JSYMC,JCENTC,NATOMC,
     &                     ISYMOP,1,NOPTYP,ANTI,EXP1EL,DENMAT,MULTA,
     &                     MULTB,IPRINT)
            ELSE IF (INTTYP.EQ.12 .OR. INTTYP.EQ.58 .OR.
     &               INTTYP.EQ.96) THEN
               CALL DSOTRA(SHLINT,EXPVAL,ISYMOP,1,NOPTYP,ANTI,DOATOM,
     &                     KAB,TRIANG,NATOM,INTADR,NSHINT,INTTYP,EXP1EL,
     &                     DENMAT,MULTA,MULTB,IPRINT)
            ELSE IF (INTTYP .EQ. 30 .OR. INTTYP .EQ. 31) THEN
               CALL EFGTRA(SHLINT,EXPVAL,FACINT,JSYMC,JCENTC,NATOMC,
     &                     ISYMOP,1,NOPTYP,ANTI,INTADR,EXP1EL,
     &                     DENMAT,MULTA,MULTB,IPRINT)
            ELSE
               DO ICOMP = 1, NOPTYP
                  ICOUNT = 0
                  DO IORBA = IDENA + 1, IDENA + KHKTA
                     DO IORBB = IDENB + 1, IDENB + KHKTB
                        FAC = 1.0D0
                        IF (ONECEN.AND.LDIAG.AND.
     &                      IORBB.NE.IORBA) FAC=0.5D0
                        IORBAB = ITRI(IORBA,IORBB)
                        ICOUNT = ICOUNT + 1
                        WORK(ICOMP)=WORK(ICOMP)
     &                       +FAC*DENMAT(IORBAB)*SHLINT(ICOUNT,ICOMP)
                     END DO
                  END DO
               END DO
               DO ICOMP = 1, NOPTYP
                  IF (INTREP(ICOMP) .EQ. 0) EXPVAL(ICOMP) =
     &                                      EXPVAL(ICOMP) + WORK(ICOMP)
               END DO
            END IF
C
C        *************************************************
C        ***** Transform integrals to symmetry basis *****
C        *************************************************
C
         ELSE IF (INTTYP .EQ. 10 .OR. INTTYP .EQ. 92) THEN
            CALL PSOTRA(SHLINT,SOINT,FACINT,JSYMC,JCENTC,NATOMC,ISYMOP,
     &                  NELMNT,NOPTYP,ANTI,IORBA,IORBB,NBAST,INTTYP,
     &                  INTADR,IPRINT)
         ELSE IF (INTTYP .EQ. 11 .OR. INTTYP .EQ. 13 .OR.
     &            INTTYP .EQ. 94) THEN
            CALL SDTRA(SHLINT,SOINT,FACINT,JSYMC,JCENTC,NATOMC,ISYMOP,
     &                 NELMNT,NOPTYP,ANTI,IORBA,IORBB,NBAST,INTTYP,
     &                 INTADR,IPRINT)
         ELSE IF (INTTYP.EQ.12 .OR. INTTYP.EQ.58 .OR. INTTYP.EQ.96) THEN
            CALL DSOTRA(SHLINT,SOINT,ISYMOP,NELMNT,NOPTYP,ANTI,DOATOM,
     &                  KAB,TRIANG,NATOM,INTADR,NSHINT,INTTYP,EXP1EL,
     &                  DUMMY,MULTA,MULTB,IPRINT)
         ELSE IF (INTTYP.EQ.14 .OR. INTTYP.EQ.44 .OR. INTTYP.EQ.98) THEN
            CALL HDOTRA(SHLINT,SOINT,ISYMOP,NELMNT,NOPTYP,IORBA,IORBB,
     &                  NBAST,INTTYP,IPRINT)
cLig <> added the call to NSSTRA for INTTYP = 82
         ELSE IF ((INTTYP .EQ. 26) .OR. (INTTYP .EQ. 27) .OR.
     &            (INTTYP .EQ. 28) .OR. (INTTYP .EQ. 38) .OR.
     &            (INTTYP .EQ. 48) .OR. (INTTYP .EQ. 74) .OR.
     &            (INTTYP .EQ. 82) .OR. (INTTYP .EQ. 93) .OR.
     &            (INTTYP .EQ. 97)) THEN
            CALL NSTTRA(SHLINT,SOINT,FACINT,JSYMC,JCENTC,NATOMC,ISYMOP,
     &                  NELMNT,NOPTYP,ANTI,IORBA,IORBB,MULTA,MULTB,
     &                  NBAST,INTTYP,IPRINT,INTADR,EXPVAL,DUMMY,EXP1EL)
         ELSE IF (INTTYP .EQ. 29) THEN
            CALL EFTRA(SHLINT,SOINT,FACINT,JSYMC,JCENTC,NATOMC,
     &                 ISYMOP,NELMNT,NOPTYP,ANTI,INTADR,IPRINT,
     &                 MULTA,MULTB,EXPVAL,DENMAT,EXP1EL)
         ELSE IF (INTTYP .EQ. 30 .OR. INTTYP .EQ. 31) THEN
            CALL EFGTRA(SHLINT,SOINT,FACINT,JSYMC,JCENTC,NATOMC,ISYMOP,
     &                  NELMNT,NOPTYP,ANTI,INTADR,EXP1EL,DUMMY,
     &                  MULTA,MULTB,IPRINT)
         ELSE IF (INTTYP .EQ. 34) THEN
            CALL HDBTRA(SHLINT,SOINT,ISYMOP,NELMNT,NOPTYP,IORBA,IORBB,
     &                  NBAST,INTADR,IPRINT)
         ELSE IF (INTTYP .EQ. 35) THEN
            CALL NPETRA(SHLINT,SOINT,FACINT,JSYMC,JCENTC,NATOMC,ISYMOP,
     &                  NELMNT,NOPTYP,ANTI,INTADR,IPRINT)
         ELSE IF (INTTYP .EQ. 52) THEN
            CALL DPGTRA(SHLINT,SOINT,ISYMOP,NELMNT,NOPTYP,ANTI,
     &                  INTADR,IPRINT)
         ELSE IF (INTTYP .EQ. 53) THEN
            CALL QUGTRA(SHLINT,SOINT,ISYMOP,NELMNT,NOPTYP,ANTI,
     &                  INTADR,IPRINT)
         ELSE IF (INTTYP .EQ. 54) THEN
            CALL OCGTRA(SHLINT,SOINT,ISYMOP,NELMNT,NOPTYP,ANTI,
     &                  INTADR,IPRINT)
         ELSE IF (INTTYP .EQ. 59) THEN
            CALL S1HTRA(SHLINT,SOINT,ISYMOP,NELMNT,NOPTYP,ANTI,IPRINT)
         ELSE IF (INTTYP .EQ. 102) THEN
            CALL AMDTRA(SHLINT,SOINT,ISYMOP,NELMNT,NOPTYP,ANTI,
     &                  INTADR,IPRINT)
         ELSE IF (INTTYP .EQ. 70) THEN
            CALL G1HTRA(SHLINT,SOINT,FACINT,JSYMC,JCENTC,NCENTC,NATOMC,
     &                  ISYMOP,NELMNT,NOPTYP,ANTI,INTADR,IPRINT)
         ELSE IF (INTTYP .EQ. 75) THEN
            CALL TESTRA(SHLINT,SOINT,FACINT,JSYMC,JCENTC,NATOMC,ISYMOP,
     &                  NELMNT,NOPTYP,ANTI,EXP1EL,DUMMY,MULTA,MULTB,
     &                  IPRINT)
         ELSE IF (INTTYP .EQ. 103) THEN
            CALL DPHTRA(SHLINT,SOINT,ISYMOP,NELMNT,NOPTYP,ANTI,
     &                  INTADR,IPRINT)
         ELSE
            FULMAT = .TRUE.
            DO 200 I = 1, NOPTYP
               IF (SQUARE) THEN
                  CALL SYMSQR(SHLINT(1,I),SOINT(1,I),INTREP(I),ISYMOP,
     &                        IORBA,IORBB,-HKAB,NBAST,IPRINT)
               ELSE IF (INTREP(I) .EQ. 0) THEN
                  CALL SYM1S(SHLINT(1,I),SOINT(1,I),ISYMOP,
     &                       MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,-HKAB,
     &                       LDIAG,FULMAT,DUM,IDUM,IPRINT)
               ELSE
                  CALL SYM1N(SHLINT(1,I),SOINT(1,I),INTREP(I),ISYMOP,
     &                       MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,-HKAB,
     &                       LDIAG,FULMAT,ANTI,DUM,IDUM,IPRINT)
               END IF
  200       CONTINUE
         END IF
      END IF
  100 CONTINUE
      RETURN
 1000 FORMAT (//,2X,'***************************************',
     &         /,2X,'******** Symmetry operation ',I2,' ********',
     &         /,2X,'***************************************',/)
      END
C  /* Deck pr1prm */
      SUBROUTINE PR1PRM(SHLINT,WORK,LWORK,LABINT,INTTYP,
     &                  NOPTYP,NBAST,IORDER,DOATOM,IPRINT,TOLS,TOLOG,
     &                  NATOMC,DISTAB,FACINT,CORCX,CORCY,CORCZ,JMAXD,
     &                  JMAXM,DIFABX,DIFABY,DIFABZ,NPOINT,WEIGHT,ABSCIS,
     &                  GEXP,TRIANG,KAB,DONUC1,DOMOM1,HUCFAC,
     &                  NCENTC,NSHINT)
#include "implicit.h"
#include "priunit.h"
C
      LOGICAL DOATOM(*), TRIANG, DONUC1, DOMOM1, HUCORB
      CHARACTER LABINT(NOPTYP)*8
      DIMENSION WORK(LWORK), FACINT(NATOMC), CORCX(NATOMC),
     &          CORCY(NATOMC), CORCZ(NATOMC), SHLINT(KCKTAB,*), GEXP(*),
     &          NCENTC(*)
#include "onecom.h"
C
      CALL QENTER('PR1PRM')
      JMAXA  = NHKTA - 1
      IF (INTTYP .EQ. 42) JMAXA = JMAXA + 1
      IF (INTTYP .EQ. 43) JMAXA = JMAXA + 2
      IF (INTTYP .EQ. 46) JMAXA = JMAXA + 1
      IF (INTTYP .EQ. 47) JMAXA = JMAXA + 2
      IF (INTTYP .EQ. 65) JMAXA = JMAXA + 2
      JMAXB  = NHKTB - 1
      JMAXT  = JMAXA + JMAXB + JMAXD + JMAXM
      IF (INTTYP .EQ. 92) JMAXT = JMAXT + 3
C
      KFRWRK = 1
      KFREE  = KFRWRK
      LFREE  = LWORK
      LODC   = 3*(JMAXA+1)*(JMAXB+1)*(JMAXT+1)*(JMAXD+1)*(JMAXM+1)
      CALL MEMGET('REAL',KODC,  LODC,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KAHGTF,NAHGTF*(NATOMC + 1),WORK,KFREE,LFREE)
C
      IF (INTTYP .NE. 12 .AND. INTTYP .NE. 58 .AND. INTTYP .NE. 96) THEN
         CALL PR1PRA(WORK(KODC),JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,
     &               WORK(KAHGTF),SHLINT,WORK(KFREE),LFREE,
     &               LABINT,INTTYP,NOPTYP,NBAST,IORDER,DOATOM,IPRINT,
     &               TOLS,TOLOG,NATOMC,DISTAB,FACINT,CORCX,CORCY,CORCZ,
     &               DIFABX,DIFABY,DIFABZ,GEXP,DONUC1,DOMOM1,HUCFAC,
     &               NCENTC,NSHINT)
      ELSE
         CALL PR1PRD(WORK(KODC),JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,
     &               WORK(KAHGTF),SHLINT,WORK(KFREE),LFREE,
     &               NPOINT,LABINT,INTTYP,NOPTYP,NBAST,
     &               DOATOM,WEIGHT,ABSCIS,TRIANG,IPRINT,
     &               TOLS,TOLOG,NATOMC,DISTAB,FACINT,CORCX,CORCY,CORCZ,
     &               DIFABX,DIFABY,DIFABZ,KAB,DONUC1,DOMOM1)
      END IF
      CALL MEMREL('PR1PRM',WORK,1,KFRWRK,KFREE,LFREE)
      CALL QEXIT('PR1PRM')
      RETURN
      END
C  /* Deck pr1pra */
      SUBROUTINE PR1PRA(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT,
     &                  WORK,LWORK,LABINT,INTTYP,NOPTYP,
     &                  NBAST,IORDER,DOATOM,IPRINT,TOLS,TOLOG,NATOMC,
     &                  DISTAB,FACINT,CORCX,CORCY,CORCZ,DIFABX,DIFABY,
     &                  DIFABZ,GEXP,DONUC1,DOMOM1,HUCFAC,NCENTC,NSHINT)
#include "implicit.h"
#include "priunit.h"
#include "aovec.h"
#include "maxorb.h"
#include "mxcent.h"
#include "pi.h"
      PARAMETER (D1 = 1.00 D00, D3 = 3.00 D00, D3INV = D1/D3,
     &           D2INV = 0.50D0)
C
      LOGICAL DOATOM(*), DIFODC, DONUC1, DOMOM1, ADDNAI, HUCORB, ONCN
      LOGICAL KINODC
      CHARACTER LABINT(NOPTYP)*8
      DIMENSION WORK(LWORK), TEMP(3),
     &          FACINT(NATOMC), CORCX(NATOMC), CORCY(NATOMC),
     &          CORCZ(NATOMC), SHLINT(KCKTAB,*), AHGTF(*), GEXP(*),
     &          ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          NCENTC(*)
#include "huckel.h"
#include "nuclei.h"
#include "primit.h"
#include "elweak.h"
C
#include "onecom.h"
#include "orgcom.h"
C
C     DIFODC true if we have JMAXD .gt. 0 (take derivatives)
cLig <> set DIFODC true also for INTTYP = 81 and 82
      DIFODC = (INTTYP.EQ. 3) .OR. (INTTYP.EQ. 5) .OR. (INTTYP.EQ.10)
     &    .OR. (INTTYP.EQ.14) .OR. (INTTYP.EQ.17) .OR. (INTTYP.EQ.18)
     &    .OR. (INTTYP.EQ.19) .OR. (INTTYP.EQ.20) .OR. (INTTYP.EQ.21)
     &    .OR. (INTTYP.EQ.23) .OR. (INTTYP.EQ.24) .OR. (INTTYP.EQ.25)
     &    .OR. (INTTYP.EQ.27) .OR. (INTTYP.EQ.28) .OR. (INTTYP.EQ.34)
     &    .OR. (INTTYP.EQ.40) .OR. (INTTYP.EQ.44) .OR. (INTTYP.EQ.50)
     &    .OR. (INTTYP.EQ.51) .OR. (INTTYP.EQ.52) .OR. (INTTYP.EQ.53)
     &    .OR. (INTTYP.EQ.54) .OR. (INTTYP.EQ.55) .OR. (INTTYP.EQ.59)
     &    .OR. (INTTYP.EQ.61) .OR. (INTTYP.EQ.63) .OR. (INTTYP.EQ.66)
     &    .OR. (INTTYP.EQ.67) .OR. (INTTYP.EQ.70) .OR. (INTTYP.EQ.71)
     &    .OR. (INTTYP.EQ.72) .OR. (INTTYP.EQ.73)
     &    .OR. (INTTYP.EQ.74) .OR. (INTTYP.EQ.80) .OR. (INTTYP.EQ.81)
     &    .OR. (INTTYP.EQ.82) .OR. (INTTYP.EQ.83) .OR. (INTTYP.EQ.84)
     &    .OR. (INTTYP.EQ.91) .OR. (INTTYP.EQ.92) .OR. (INTTYP.EQ.93)
     &    .OR. (INTTYP.EQ.94) .OR. (INTTYP.EQ.97) .OR. (INTTYP.EQ.98)
     &    .OR. (INTTYP.EQ.99) .OR. (INTTYP.EQ.100).OR. (INTTYP.EQ.101)
     &    .OR. (INTTYP.EQ.102).OR. (INTTYP.EQ.103)
      ADDNAI = (INTTYP.EQ. 5) .OR. (INTTYP.EQ.19) .OR. (INTTYP.EQ.20)
     &    .OR. (INTTYP.EQ.24) .OR. (INTTYP.EQ.25) .OR. (INTTYP.EQ.51)
     &    .OR. (INTTYP.EQ.57) .OR. (INTTYP.EQ.61) .OR. (INTTYP.EQ.62)
     &    .OR. (INTTYP.EQ.63) .OR. (INTTYP.EQ.64) .OR. (INTTYP.EQ.72)
     &    .OR. (INTTYP.EQ.73) .OR. (INTTYP.EQ.80)
Ckr
Ckr   We may contemplate setting ADDNAI equal to true for EFBINT and EFB2IN
Ckr
C
      DO 100 IPRIMA = 1,NUCA
         JPRIMA = JSTA + IPRIMA
         CONTA = PRICCF(JPRIMA,NUMCFA)
         EXPA = PRIEXP(JPRIMA)
      DO 100 IPRIMB = 1,NUCB
         JPRIMB = JSTB + IPRIMB
         CONTB = PRICCF(JPRIMB,NUMCFB)
         EXPB = PRIEXP(JPRIMB)
         EXPP = EXPA + EXPB
         EXPPI = D1/EXPP
         EXPABQ = EXPA*EXPB*DISTAB*EXPPI
         IF (EXPABQ.GT.TOLOG) GO TO 200
         SAAB = CONTA*CONTB*EXP(-EXPABQ)
         ASAAB = ABS(SAAB)
         IF (ASAAB.LT.TOLS) GO TO 200
         SAAB13 = SIGN(ASAAB**D3INV,SAAB)
C
C        Calculate coordinates of product Gaussian P
C
         EXPAPI = EXPA*EXPPI
         EXPBPI = EXPB*EXPPI
         CORPX  = EXPAPI*CORAX + EXPBPI*CORBX
         CORPY  = EXPAPI*CORAY + EXPBPI*CORBY
         CORPZ  = EXPAPI*CORAZ + EXPBPI*CORBZ
C
C        *********************************************
C        ***** Overlap Distribution Coefficients *****
C        *********************************************
C
C        Expansion coefficients
C
         DIFPAX = CORPX - CORAX
         DIFPAY = CORPY - CORAY
         DIFPAZ = CORPZ - CORAZ
         DIFPBX = CORPX - CORBX
         DIFPBY = CORPY - CORBY
         DIFPBZ = CORPZ - CORBZ
C
C     Choose origin for expansions coefficients:
C     1) A certain nuclei
C     2) Dipole origin
C     3) Gauge origin
C     4) Cavity center ("molecular" origin)
C     5) Center P of overlap distribution
C
         IF ((INTTYP .EQ. 17) .OR. (INTTYP .EQ. 20) .OR.
     &       (INTTYP .EQ. 22) .OR. (INTTYP .EQ. 26)) THEN
            TEMP(1) = CORBX
            TEMP(2) = CORBY
            TEMP(3) = CORBZ
         ELSE IF ((INTTYP .EQ. 2)  .OR. (INTTYP .EQ. 4)  .OR.
     &            (INTTYP .EQ. 6)  .OR. (INTTYP .EQ. 7)  .OR.
     &            (INTTYP .EQ. 8)  .OR. (INTTYP .EQ. 29) .OR.
     &            (INTTYP .EQ. 30) .OR. (INTTYP .EQ. 31) .OR.
     &            (INTTYP .EQ. 56) .OR. (INTTYP .EQ. 53) .OR.
     &            (INTTYP .EQ. 54) .OR. (INTTYP .EQ. 83) .OR.
     &            (INTTYP .EQ. 84) ) THEN
            TEMP(1) = DIPORG(1)
            TEMP(2) = DIPORG(2)
            TEMP(3) = DIPORG(3)
cLig <> put the origin in GAGORG also for INTTYP = 81 and 82
         ELSE IF ((INTTYP .EQ. 18) .OR. (INTTYP .EQ. 37) .OR.
     &            (INTTYP .EQ. 38) .OR. (INTTYP .EQ. 57) .OR.
     &            (INTTYP .EQ. 61) .OR. (INTTYP .EQ. 81) .OR.
     &            (INTTYP .EQ. 82) .OR. (INTTYP .EQ. 91) .OR.
     &            (INTTYP .EQ. 93) .OR. (INTTYP .EQ. 97) .OR.
     &            (INTTYP .EQ.101) .OR. (INTTYP .EQ. 102)) THEN
            TEMP(1) = GAGORG(1)
            TEMP(2) = GAGORG(2)
            TEMP(3) = GAGORG(3)
         ELSE IF ((INTTYP .EQ. 46) .OR. (INTTYP .EQ. 47)) THEN
            TEMP(1) = CAVORG(1)
            TEMP(2) = CAVORG(2)
            TEMP(3) = CAVORG(3)
         ELSE IF (INTTYP .EQ. 48) THEN
C
C     We want to place the origin at the nucleus investigated
C
            KR_SET = 0
            DO 37 KR = 1, NUCIND
               IF (DOATOM(KR)) THEN
                  IF (KR_SET .NE. 0) THEN
                    CALL QUIT('More than one DOATOM true for INTTYP 48')
                  END IF
                  KR_SET = KR
                  TEMP(1) = CORD(1,KR)
                  TEMP(2) = CORD(2,KR)
                  TEMP(3) = CORD(3,KR)
               END IF
 37         CONTINUE
         ELSE
            TEMP(1) = ORIGIN(1)
            TEMP(2) = ORIGIN(2)
            TEMP(3) = ORIGIN(3)
         END IF
C
C        * NEXT IF-STATEMENT CHANGED TO INCLUDE DARWIN TERM*
C
         IF (INTTYP .NE. 9 .AND. INTTYP .NE. 41 .AND.
     &       INTTYP .NE. 95) THEN
Clf
            KINODC = .FALSE.
Clf
            CALL GETODC(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,DIFODC,
     &                  KINODC,ONECEN,EXPA,EXPB,IPRINT,SAAB13,EXPPI,
     &                  WORK,LWORK,CORPX,CORPY,CORPZ,DONUC1,DOMOM1,
     &                  TEMP,INTTYP)
         END IF
C
C        **********************************************
C        ***** Calculation of Hermitian integrals *****
C        **********************************************
C
C        Overlap integral
C
         SHGTF = SQRT(PI*EXPPI)
C
C        Nuclear attraction integrals for spin-orbit and London contribution
C        to angular momentum.
C
cLig <> added the call to DZERO also for INTTYP = 82
         IF ((INTTYP .EQ. 5) .OR. (INTTYP .EQ. 10)
     &                       .OR. (INTTYP .EQ. 11)
     &                       .OR. (INTTYP .EQ. 13)
     &                       .OR. (INTTYP .EQ. 19)
     &                       .OR. (INTTYP .EQ. 20)
     &                       .OR. (INTTYP .EQ. 24)
     &                       .OR. (INTTYP .EQ. 25)
     &                       .OR. (INTTYP .EQ. 26)
     &                       .OR. (INTTYP .EQ. 27)
     &                       .OR. (INTTYP .EQ. 28)
     &                       .OR. (INTTYP .EQ. 29)
     &                       .OR. (INTTYP .EQ. 30)
     &                       .OR. (INTTYP .EQ. 31)
     &                       .OR. (INTTYP .EQ. 35)
     &                       .OR. (INTTYP .EQ. 38)
     &                       .OR. (INTTYP .EQ. 48)
     &                       .OR. (INTTYP .EQ. 51)
     &                       .OR. (INTTYP .EQ. 57)
     &                       .OR. (INTTYP .EQ. 61)
     &                       .OR. (INTTYP .EQ. 62)
     &                       .OR. (INTTYP .EQ. 63)
     &                       .OR. (INTTYP .EQ. 64)
     &                       .OR. (INTTYP .EQ. 66)
     &                       .OR. (INTTYP .EQ. 67)
     &                       .OR. (INTTYP .EQ. 70)
     &                       .OR. (INTTYP .EQ. 72)
     &                       .OR. (INTTYP .EQ. 73)
     &                       .OR. (INTTYP .EQ. 74)
     &                       .OR. (INTTYP .EQ. 75)
     &                       .OR. (INTTYP .EQ. 80)
     &                       .OR. (INTTYP .EQ. 82)
     &                       .OR. (INTTYP .EQ. 92)
     &                       .OR. (INTTYP .EQ. 93)
     &                       .OR. (INTTYP .EQ. 94)
     &                       .OR. (INTTYP .EQ. 97)
     &                       .OR. (INTTYP .EQ. 99)
     &                       .OR. (INTTYP .EQ. 100)
     &    ) THEN
               CALL DZERO(AHGTF,(NATOMC + 1)*NAHGTF)
            IF (ADDNAI) THEN
               CALL DZERO(AHGTF,NAHGTF)
               IADR = 1 + NAHGTF
            ELSE
               IADR = 1
            END IF
            FACTOR = D1
            DO 300 IATOMC = 1, NATOMC
               IF (ADDNAI .OR. INTTYP.EQ.70) FACTOR = FACINT(IATOMC)
               DIFCPX = CORCX(IATOMC) - CORPX
               DIFCPY = CORCY(IATOMC) - CORPY
               DIFCPZ = CORCZ(IATOMC) - CORPZ
               IF (GNUEXP(IATOMC) .GT. 0.0D0) THEN
                  PEXP=EXPP*GNUEXP(IATOMC)/(EXPP+GNUEXP(IATOMC))
                  GNUFAC=FACTOR * (GNUEXP(IATOMC)/
     &                 (EXPP + GNUEXP(IATOMC)))**1.5D0
                  CALL HERNAI(AHGTF,JMAX,PEXP,DIFCPX,DIFCPY,DIFCPZ,
     &                        GNUFAC,IADR,ISTEPU,ISTEPV,NAHGTF,IPRINT)
               ELSE
                  CALL HERNAI(AHGTF,JMAX,EXPP,DIFCPX,DIFCPY,DIFCPZ,
     &                        FACTOR,IADR,ISTEPU,ISTEPV,NAHGTF,IPRINT)
               END IF

               IF (ADDNAI) THEN
                  CALL DAXPY(NAHGTF,D1,AHGTF(IADR),1,AHGTF(1),1)
               ELSE
                  IADR = IADR + NAHGTF
               END IF
  300       CONTINUE
         END IF
C
C        **********************************************
C        ***** Calculation of Cartesian integrals *****
C        **********************************************
C
C        Overlap integrals
C        -----------------
C
         IF (ABS(INTTYP) .EQ. 1 .OR. INTTYP .EQ. 45) THEN
            CALL OVLINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,HUCFAC,
     &                  SHGTF,SHLINT,INTTYP)
C          
C        POPOVLP /hjaaj Mar 2006
C        Modified Overlap integrals for population analysis
C        --------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 1001) THEN
            CALL POPOVLINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,
     &                     EXPBPI,SHGTF,SHLINT,INTTYP)
C
C        Dipole length integrals
C        -----------------------
C
         ELSE IF (INTTYP .EQ. 2) THEN
            CALL DPLINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,DIPORG)
C
C        Dipole velocity and half-derivative overlap integrals
C        -----------------------------------------------------
C
         ELSE IF ((INTTYP .EQ. 3) .OR. (INTTYP .EQ. 14) .OR.
     &            (INTTYP .EQ. 44).OR. (INTTYP .EQ. 59)) THEN
            CALL DPVINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT)
C
C        <d2/dx2>, <d2/dy2>, <d2/dz2> integrals for SQH2DO /hjaaj Oct 2004
C        -------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 98) THEN
            CALL DV2INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT)
C
C        Quadrupole integrals
C        --------------------
C
         ELSE IF (INTTYP .EQ. 4) THEN
            CALL QDPINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,EXPPI,INTTYP)
C
C        Spin-orbit integrals
C        --------------------
C
         ELSE IF (INTTYP .EQ. 5) THEN
            CALL SSOINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT)
C
C        Second moments integrals
C        ------------------------
C
         ELSE IF (INTTYP .EQ. 6) THEN
            CALL QDPINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,EXPPI,INTTYP)
C
C        Traceless theta quadrupole moments integrals
C        --------------------------------------------
C
         ELSE IF (INTTYP .EQ. 7) THEN
            CALL QDPINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,EXPPI,INTTYP)
C
C        Multipole moment integrals
C        --------------------------
C
         ELSE IF (INTTYP .EQ. 8) THEN
            IF (FMMORI) THEN
               DIPORG(1) = CORPX
               DIPORG(2) = CORPY
               DIPORG(3) = CORPZ
            END IF
            CALL MOMINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,EXPPI,DIPORG,IORDER,WORK,
     &                  LWORK)
C
C        Fermi contact integrals
C        -----------------------
C
         ELSE IF (INTTYP .EQ. 9) THEN
            CALL FRMINT(SHLINT,NOPTYP,EXPP,CORPX,CORPY,CORPZ,SAAB,
     &                  DOATOM)
C
C        Paramagnetic spin-orbit integrals
C        ---------------------------------
C
         ELSE IF (INTTYP .EQ. 10) THEN
            CALL PSOINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT,
     &                  NOPTYP,NATOMC)
C
C        Spin-dipole integrals
C        ---------------------
C
         ELSE IF (INTTYP .EQ. 11 .OR. INTTYP .EQ. 13) THEN
            CALL SDINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT,
     &                 NOPTYP,NATOMC,INTTYP,SAAB,EXPP,CORCX,CORCY,CORCZ,
     &                 CORPX,CORPY,CORPZ)
C
C        Magnetic derivatives of overlap matrix
C        --------------------------------------
C
         ELSE IF (INTTYP .EQ. 15) THEN
            CALL SM1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  DIFABX,DIFABY,DIFABZ)
C
C        Second order magnetic derivatives of overlap matrix
C        ---------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 16) THEN
            CALL SM2INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  EXPPI,DIFABX,DIFABY,DIFABZ)
C
C        Angular momentum around the nuclei
C        ----------------------------------
C
         ELSE IF (INTTYP .EQ. 17) THEN
            CALL AM1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  D1)
C
C        Angular momentum around the origin
C        ----------------------------------
C
         ELSE IF (INTTYP .EQ. 18) THEN
            CALL AM1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  D1)
C
C        Differentiated angular momentum around the origin
C        -------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 102) THEN
            CALL AMDINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  EXPA)
            DONUC1 = .TRUE.
            DOMOM1 = .FALSE.
            CALL GETODC(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,DIFODC,
     &                  .FALSE.,ONECEN,EXPA,EXPB,IPRINT,SAAB13,EXPPI,
     &                  WORK,LWORK,CORPX,CORPY,CORPZ,DONUC1,DOMOM1,
     &                  TEMP,INTTYP)
            CALL AMDINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,
     &                  SHLINT(1,10),EXPA)
            DONUC1 = .FALSE.
            DOMOM1 = .TRUE.
C
C        London orbital contribution to magnetic moment
C        ----------------------------------------------
C
         ELSE IF (INTTYP .EQ. 19 .OR. INTTYP .EQ. 66) THEN
            CALL MM1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  AHGTF,DIFABX,DIFABY,DIFABZ,NATOMC,INTTYP)
C
C        One-electron contribution to magnetic moment
C        --------------------------------------------
C
         ELSE IF (INTTYP .EQ. 20) THEN
            CALL AM1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  D2INV)
            CALL GETODC(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,DIFODC,
     &                  .FALSE.,ONECEN,EXPA,EXPB,IPRINT,SAAB13,EXPPI,
     &                  WORK,LWORK,CORPX,CORPY,CORPZ,DONUC1,DOMOM1,
     &                  ORIGIN,INTTYP)
            CALL MM1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  AHGTF,DIFABX,DIFABY,DIFABZ,NATOMC,INTTYP)
C
C        Kinetic energy integrals
C        ------------------------
C
         ELSE IF (INTTYP .EQ. 21) THEN
            CALL KININT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT)
C
C        Diamagnetic susceptiblity
C        -------------------------
C
         ELSE IF (INTTYP .EQ. 22) THEN
            CALL QDPINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,EXPPI,INTTYP)
C
C        Angular London orbital contribution to magnetic susceptibility
C        --------------------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 23) THEN
            CALL DSUSAN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,
     &                  SHLINT,CORPX,CORPY,CORPZ,DIFABX,DIFABY,DIFABZ)
C
C        London orbital contribution to magnetic susceptibility
C        ------------------------------------------------------
C
         ELSE IF ((INTTYP .EQ. 24).OR.(INTTYP .EQ. 67)) THEN
            CALL DSUSLH(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  AHGTF,DIFABX,DIFABY,DIFABZ,NATOMC,INTTYP)
C
C        Diamagnetic susceptibility
C        --------------------------
C
         ELSE IF (INTTYP .EQ. 25) THEN
            CALL DSUSAN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,
     &                  SHLINT,CORPX,CORPY,CORPZ,DIFABX,DIFABY,DIFABZ)
            CALL DSUSLH(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  AHGTF,DIFABX,DIFABY,DIFABZ,NATOMC,INTTYP)
            TEMP(1) = CORBX
            TEMP(2) = CORBY
            TEMP(3) = CORBZ
            CALL GETODC(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,DIFODC,
     &                  .FALSE.,ONECEN,EXPA,EXPB,IPRINT,SAAB13,EXPPI,
     &                  WORK,LWORK,CORPX,CORPY,CORPZ,DONUC1,DOMOM1,
     &                  TEMP,INTTYP)
            CALL QDPINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,EXPPI,INTTYP)
C
C        Nuclear shielding integrals
C        ---------------------------
C
         ELSE IF (INTTYP .EQ. 26) THEN
            CALL NSNLIN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,
     &                  SHLINT,NOPTYP,NATOMC)
C
C        London orbital contribution to nuclear shielding tensor
C        -------------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 27) THEN
            CALL NSLOIN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,
     &                  SHLINT,NOPTYP,NATOMC,DIFABX,DIFABY,DIFABZ)
C
C        Nuclear shielding tensor integrals
C        ----------------------------------
C
         ELSE IF (INTTYP .EQ. 28) THEN
            CALL NSLOIN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,
     &                  SHLINT,NOPTYP,NATOMC,DIFABX,DIFABY,DIFABZ)
            TEMP(1) = CORBX
            TEMP(2) = CORBY
            TEMP(3) = CORBZ
            CALL GETODC(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,DIFODC,
     &                  .FALSE.,ONECEN,EXPA,EXPB,IPRINT,SAAB13,EXPPI,
     &                  WORK,LWORK,CORPX,CORPY,CORPZ,DONUC1,DOMOM1,
     &                  TEMP,INTTYP)
            CALL NSNLIN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,
     &                  SHLINT,NOPTYP,NATOMC)
C
C        Electric field from the individual nuclei
C        -----------------------------------------
C
         ELSE IF (INTTYP .EQ. 29) THEN
            CALL EF1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,
     &                  SHLINT,NOPTYP,NATOMC)
C
C        Electric field gradient from the individual nuclei
C        --------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 30 .OR. INTTYP .EQ. 31) THEN
            CALL EFGINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,
     &                  SHLINT,NOPTYP,NATOMC,DISTAB)
C
C        Bra-differentiated overlap matrix with respect to magnetic field
C        ----------------------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 32) THEN
            DIFAOX = CORAX - GAGORG(1)
            DIFAOY = CORAY - GAGORG(2)
            DIFAOZ = CORAZ - GAGORG(3)
            CALL SM1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  DIFAOX,DIFAOY,DIFAOZ)
C
C        Ket-differentiated overlap matrix with respect to magnetic field
C        ----------------------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 33) THEN
            DIFOBX = GAGORG(1) - CORBX
            DIFOBY = GAGORG(2) - CORBY
            DIFOBZ = GAGORG(3) - CORBZ
            CALL SM1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  DIFOBX,DIFOBY,DIFOBZ)
C
C        Ket-differentiated HDO-integrals with respect to magnetic field
C        ---------------------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 34) THEN
            DIFOBX = GAGORG(1) - CORBX
            DIFOBY = GAGORG(2) - CORBY
            DIFOBZ = GAGORG(3) - CORBZ
            CALL HDBINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  DIFOBX,DIFOBY,DIFOBZ)
C
C        Potential energy from each individual nucleus
C        ---------------------------------------------
C
         ELSE IF (INTTYP .EQ. 35 .OR. INTTYP .EQ. 75) THEN
            CALL NPEINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,
     &                  SHLINT,NOPTYP,NATOMC)
C
C        Half B-differentiated overlap matrix
C        ------------------------------------
C
         ELSE IF (INTTYP .EQ. 36) THEN
            DIFOBX = CORBX - GAGORG(1)
            DIFOBY = CORBY - GAGORG(2)
            DIFOBZ = CORBZ - GAGORG(3)
            CALL SM1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  DIFOBX,DIFOBY,DIFOBZ)
            DIFAOX = CORAX - GAGORG(1)
            DIFAOY = CORAY - GAGORG(2)
            DIFAOZ = CORAZ - GAGORG(3)
            CALL SM1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  DIFAOX,DIFAOY,DIFAOZ)
C
C        Diamagnetic susceptiblity with common gauge origin
C        --------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 37) THEN
            CALL QDPINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,EXPPI,INTTYP)
C
C        Nuclear shielding integrals with common gauge origin
C        ----------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 38 .OR. INTTYP .EQ. 48) THEN
            CALL NSNLIN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,
     &                  SHLINT,NOPTYP,NATOMC)
C
C        Magnetic-field correction to spin-orbit integrals
C        -------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 57) THEN
            CALL SOFINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,
     &                  SHLINT,NOPTYP,NATOMC)
C
C        Cosine and sine integrals
C        -------------------------
C
         ELSE IF (INTTYP .EQ. 39) THEN
            CALL GOSINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,ORIGIN,EXPPI)
C
C        Mass velocity integrals
C        -----------------------
C
         ELSE IF (INTTYP .EQ. 40) THEN
            CALL MVLINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT)
C
C        Darwin term integrals
C        ---------------------
C
         ELSE IF (INTTYP .EQ. 41) THEN
            CALL DWNINT(SHLINT,EXPP,CORPX,CORPY,CORPZ,SAAB)
C
C        Electric field contribution to magnetic moment
C        ----------------------------------------------
C
         ELSE IF (INTTYP .EQ. 42) THEN
            CALL CM1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  EXPPI,DIFABX,DIFABY,DIFABZ)
C
C        Electric field contribution to diamagnetic magnetizability
C        ----------------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 43) THEN
            CALL CM2INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  EXPPI,DIFABX,DIFABY,DIFABZ)
C
C        First magnetic derivative of solvent integrals
C        ----------------------------------------------
C
         ELSE IF (INTTYP .EQ. 46) THEN
            CALL SL1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,DIFABX,DIFABY,DIFABZ,EXPPI,
     &                  IORDER,WORK,LWORK)
C
C        Second magnetic derivative of solvent integrals
C        -----------------------------------------------
C
         ELSE IF (INTTYP .EQ. 47) THEN
            CALL SL2INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,DIFABX,DIFABY,DIFABZ,EXPPI,
     &                  IORDER,WORK,LWORK)
C
C        First electric derivative of overlap integrals. Type A
C        ------------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 49) THEN
            CALL SE1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT)
C
C        First electric derivative of overlap integrals. Type B
C        ------------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 50) THEN
            CALL SE1INB(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,
     $           SHLINT, EXPA, EXPB)
C
C        First electric deriv. of 1-electron Ham. integrals
C        --------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 51) THEN
            CALL H1EINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,
     $           SHLINT,AHGTF,NATOMC)
C
C        Dipole gradient integrals
C        -------------------------
C
         ELSE IF (INTTYP .EQ. 52) THEN
            CALL DPGINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,
     &                  SHGTF,SHLINT(1,1))
            DONUC1 = .FALSE.
            DOMOM1 = .TRUE.
            CALL GETODC(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,DIFODC,
     &                  .FALSE.,ONECEN,EXPA,EXPB,IPRINT,SAAB13,EXPPI,
     &                  WORK,LWORK,CORPX,CORPY,CORPZ,DONUC1,DOMOM1,
     &                  TEMP,INTTYP)
            DONUC1 = .TRUE.
            DOMOM1 = .FALSE.
            CALL DPGINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,
     &                  SHGTF,SHLINT(1,10))
C
C        Quadrupole gradient integrals
C        -----------------------------
C
         ELSE IF (INTTYP .EQ. 53) THEN
            CALL QUGINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,
     &                  SHGTF,SHLINT(1,1))
            DONUC1 = .FALSE.
            DOMOM1 = .TRUE.
            CALL GETODC(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,DIFODC,
     &                  .FALSE.,ONECEN,EXPA,EXPB,IPRINT,SAAB13,EXPPI,
     &                  WORK,LWORK,CORPX,CORPY,CORPZ,DONUC1,DOMOM1,
     &                  TEMP,INTTYP)
            DONUC1 = .TRUE.
            DOMOM1 = .FALSE.
            CALL QUGINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,
     &                  SHGTF,SHLINT(1,19))
C
C        Octupole gradient integrals
C        ---------------------------
C
         ELSE IF (INTTYP .EQ. 54) THEN
            DSIGN = 1.0D0
            CALL OCGINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,
     &                  SHGTF,SHLINT(1,1))
            DONUC1 = .FALSE.
            DOMOM1 = .TRUE.
            CALL GETODC(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,DIFODC,
     &                  .FALSE.,ONECEN,EXPA,EXPB,IPRINT,SAAB13,EXPPI,
     &                  WORK,LWORK,CORPX,CORPY,CORPZ,DONUC1,DOMOM1,
     &                  TEMP,INTTYP)
            DONUC1 = .TRUE.
            DOMOM1 = .FALSE.
            CALL OCGINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,
     &                  SHGTF,SHLINT(1,31))
C
C        Rotational strength integrals
C        -----------------------------
C
         ELSE IF (INTTYP .EQ. 55) THEN
            CALL RSTINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT)
C
C        Third moments integrals
C        ------------------------
C
         ELSE IF (INTTYP .EQ. 56) THEN
            CALL THMINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT)
C
C        Diamagnetic one-electron spin-orbit integral without London
C        -----------------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 61) THEN
               CALL DELGIN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,
     &                  SHLINT,NOPTYP)
C
C        Diamagnetic one-electron spin-orbit integral with London
C        --------------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 62) THEN
               CALL DELGIN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,
     &                  SHLINT,NOPTYP)
C
C        Douglas-Kroll pVp integrals
C        ------------------------
C
         ELSE IF (INTTYP .EQ. 63) THEN
            CALL PVPINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT,
     &                  NATOMC,EXPB)
C
C        Potential energy integrals
C        ------------------------
C
         ELSE IF (INTTYP .EQ. 64) THEN
            CALL POTINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT)
C
C        Electric field contribution to magnetic moment
C        ----------------------------------------------
C
         ELSE IF (INTTYP .EQ. 65) THEN
            CALL QDBINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  EXPPI,DIFABX,DIFABY,DIFABZ)
C
C        Geometrical derivative of one-electron Hamiltonian
C        --------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 70) THEN
            CALL G1HINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  AHGTF,NATOMC,NOPTYP,NCENTC)
C
C        Overlap integrals for the small-component in DPT
C        ------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 71) THEN
            CALL RELOVL(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,
     &                  SHLINT)
C
C        Potential energy int for the small-component in DPT
C        ---------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 72) THEN
            CALL RELPO1(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,
     &                  SHLINT,AHGTF,NATOMC,EXPA,EXPB)
         ELSE IF (INTTYP .EQ. 73) THEN
            CALL RELPO2(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,
     &                  SHLINT,AHGTF,NATOMC,EXPA,EXPB)
C
C        PSO DPT-lookalike integrals
C        ---------------------------
C
         ELSE IF (INTTYP .EQ. 74) THEN
            CALL NSYINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT,
     &                  NOPTYP,NATOMC)
C
C        Parity violating electroweak interaction
C        ----------------------------------------
         ELSE IF (INTTYP .EQ. 80) THEN
            CALL CNTFPV(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT)
C
cLig >> added the call to RAMINT and RMIINT
C
C        (r-r')l' integrals
C        ------------------
C
         ELSE IF (INTTYP .EQ. 81) THEN
            CALL RAMINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT)
C
C        (r-r')l'/|r-R_I|**3 integrals
C        -----------------------------
C
         ELSE IF (INTTYP .EQ. 82) THEN
            CALL RMIINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,
     &                  SHLINT,NOPTYP,NATOMC)
cLig <<
C
C        DPT correction to dipole
C        ------------------------
C
         ELSE IF (INTTYP .EQ. 83) THEN
            CALL PRPINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,DIPORG)
C
C        DPT correction to dipole
C        ------------------------
C
         ELSE IF (INTTYP .EQ. 84) THEN
            CALL PXPINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,DIPORG)
C
C        Kinetic energy correction to orbital Zeeman
C        -------------------------------------------
C
         ELSE IF (INTTYP .EQ. 91) THEN
            CALL OZKEIN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT)
C
C        Kinetic energy correction to paramagnetic spin-orbit integrals
C        --------------------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 92) THEN
            CALL PSOKIN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT,
     &                  NOPTYP,NATOMC)
C
C        Kinetic energy correction to diamagnetic nuclear shieldings
C        -----------------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 93) THEN
            CALL NSKEIN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,
     &                  SHLINT,NOPTYP,NATOMC)
C
C        Kinetic energy correction to spin-dipole integrals
C        --------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 94) THEN
            CALL SDKINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT,
     &                 NOPTYP,NATOMC,INTTYP,SAAB,EXPP,EXPB,CORCX,CORCY,
     &                 CORCZ,CORPX,CORPY,CORPZ)
C
C        Kinetic energy correction to Fermi contact integrals
C        ----------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 95) THEN
            CALL FRMKIN(SHLINT,NOPTYP,EXPP,EXPB,CORPX,CORPY,CORPZ,SAAB,
     &                  DOATOM)
C
C        Orbital Zeeman correction to PSO integrals
C        ------------------------------------------
C
         ELSE IF (INTTYP .EQ. 97) THEN
            CALL PSOZIN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,
     &                  SHLINT,NOPTYP,NATOMC)
C
C        London orbital contribution to electric field at point in space
C        ---------------------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 99) THEN
            CALL EFBINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT,
     &                  DIFABX,DIFABY,DIFABZ)
C
C        London orbital contribution to electric field at point in space
C        ---------------------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 100) THEN
            CALL EFB2IN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT,
     &                  DIFABX,DIFABY,DIFABZ)
C
C        Angular London orbital contribution to magnetic susceptibility
C        --------------------------------------------------------------
C
         ELSE IF (INTTYP .EQ. 101) THEN
            CALL MQDPIN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT)
C
C        Dipole Hessian integrals
C        -------------------------
C
         ELSE IF (INTTYP .EQ. 103) THEN
            IF (.NOT. ONECEN) THEN
               CALL DPHINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,
     &                     SHGTF,SHLINT(1,1))
               DONUC1 = .FALSE.
               DOMOM1 = .TRUE.
               CALL GETODC(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,DIFODC,
     &                     .FALSE.,ONECEN,EXPA,EXPB,IPRINT,SAAB13,EXPPI,
     &                     WORK,LWORK,CORPX,CORPY,CORPZ,DONUC1,DOMOM1,
     &                     TEMP,INTTYP)
               DONUC1 = .TRUE.
               DOMOM1 = .FALSE.
               CALL DPHINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,
     &                     SHGTF,SHLINT(1,19))
C
            END IF
         END IF
  200 CONTINUE
  100 CONTINUE
C
C     Print
C
      IF (IPRINT .GT. 30) THEN
         DO 400 I = 1, NSHINT
            CALL AROUND('SHLINT for '//LABINT(I)//' in PR1PRA')
            CALL OUTPUT(SHLINT(1,I),1,KCKTA,1,KCKTB,KCKTA,KCKTB,1,LUPRI)
  400    CONTINUE
      END IF
      RETURN
      END
C  /* Deck pr1prd */
      SUBROUTINE PR1PRD(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT,
     &                  WORK,LWORK,NPOINT,LABINT,INTTYP,NOPTYP,NBAST,
     &                  DOATOM,WEIGHT,ABSCIS,TRIANG,IPRINT,TOLS,TOLOG,
     &                  NATOMC,DISTAB,FACINT,CORCX,CORCY,CORCZ,DIFABX,
     &                  DIFABY,DIFABZ,KAB,DONUC1,DOMOM1)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "nuclei.h"
      PARAMETER (D1 = 1.00 D00, D3 = 3.00 D00, D3INV = D1/D3)
C
      LOGICAL DOATOM(*), TRIANG, SAMECD, DONUC1, DOMOM1
      CHARACTER LABINT(NOPTYP)*8
      DIMENSION WORK(LWORK), WEIGHT(NPOINT), ABSCIS(NPOINT),
     &          SHLINT(KCKTAB,*), AHGTF(*),
     &          ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,3)
      DIMENSION FACINT(NATOMC),
     &          CORCX(NATOMC), CORCY(NATOMC), CORCZ(NATOMC)
C
#include "onecom.h"
#include "orgcom.h"
#include "primit.h"
#include "symmet.h"
C
#include "ibtfun.h"
      XYZ(I,J) = PT(IBTAND(ISYMAX(I,1),J))
C
      DO 100 IPRIMA = 1,NUCA
         JPRIMA = JSTA + IPRIMA
         CONTA  = PRICCF(JPRIMA,NUMCFA)
         EXPA   = PRIEXP(JPRIMA)
      DO 100 IPRIMB = 1,NUCB
         JPRIMB = JSTB + IPRIMB
         CONTB  = PRICCF(JPRIMB,NUMCFB)
         EXPB   = PRIEXP(JPRIMB)
C
C        Loop over fourth center, atom D
C        -------------------------------
C
         ISTDSO = 1
         DO 200 IATOMD = 1, NUCIND
         IF (DOATOM(IATOMD)) THEN
            MULD = ISTBNU(IATOMD)
            MABD = IBTOR(MULD,KAB)
            KABD = IBTAND(MULD,KAB)
            DO 300 ISYMD = 0, MAXOPR
            IF (IBTAND(ISYMD,MABD).EQ.0) THEN
               CORDX = XYZ(1,ISYMD)*CORD(1,IATOMD)
               CORDY = XYZ(2,ISYMD)*CORD(2,IATOMD)
               CORDZ = XYZ(3,ISYMD)*CORD(3,IATOMD)
               DISTAD = (CORAX - CORDX)**2 + (CORAY - CORDY)**2
     &                                     + (CORAZ - CORDZ)**2
               DISTBD = (CORBX - CORDX)**2 + (CORBY - CORDY)**2
     &                                     + (CORBZ - CORDZ)**2
C
C              Prepare for loop over atom C
C              ----------------------------
C
               JATOMC = 0
               MXATMC = NUCIND
               IF (TRIANG .AND. (INTTYP .EQ. 12 .OR. INTTYP .EQ. 96))
     &              MXATMC = IATOMD
               DO 400 IATOMC = 1, MXATMC
               IF (DOATOM(IATOMC)) THEN
                  SAMECD = IATOMC .EQ. IATOMD
                  MULC   = ISTBNU(IATOMC)
                  MABCD  = IBTOR(MULC,KABD)
                  FACTOR = - FMULT(IBTAND(MULC,KABD))
                  DO 500 ISYMC = 0, MAXOPR
                  IF (SAMECD .AND. ISYMC.EQ.ISYMD) GO TO 500
                  IF (IBTAND(ISYMC,MABCD) .EQ. 0) THEN
                     JATOMC = JATOMC + 1
                     CORCX(JATOMC) = XYZ(1,ISYMC)*CORD(1,IATOMC)
                     CORCY(JATOMC) = XYZ(2,ISYMC)*CORD(2,IATOMC)
                     CORCZ(JATOMC) = XYZ(3,ISYMC)*CORD(3,IATOMC)
                     FACINT(JATOMC)= FACTOR
                  END IF
  500             CONTINUE
               END IF
  400          CONTINUE
               IF (JATOMC .EQ. 0) GO TO 300
               IF (JATOMC .GT. NATOMC) THEN
                 WRITE(LUPRI,*) 'DSO error, JATOMC,NATOMC',JATOMC,NATOMC
                 CALL QUIT(
     &            'JATOMC.gt.NATOMC, insufficient allocation for AHGTF')
               END IF
C
C              Loop over quadrature points
C              ---------------------------
C
               DO 600 IPOINT = 1, NPOINT
                  EXPD   = ABSCIS(IPOINT)
                  EXPP   = EXPA + EXPB + EXPD
                  EXPPI  = D1/EXPP
                  EXPABQ = EXPPI*(EXPA*EXPB*DISTAB
     &                   + EXPA*EXPD*DISTAD + EXPB*EXPD*DISTBD)
               IF (EXPABQ.GT.TOLOG) GO TO 600
                  SAAB   = WEIGHT(IPOINT)*CONTA*CONTB
     &                                   *EXP(-EXPABQ)
                  ASAAB  = ABS(SAAB)
               IF (ASAAB.LT.TOLS) GO TO 600
                  SAAB13 = SIGN(ASAAB**D3INV,SAAB)
                  EXPAPI = EXPA*EXPPI
                  EXPBPI = EXPB*EXPPI
                  EXPCPI = EXPD*EXPPI
                  CORPX =EXPAPI*CORAX+EXPBPI*CORBX+EXPCPI*CORDX
                  CORPY =EXPAPI*CORAY+EXPBPI*CORBY+EXPCPI*CORDY
                  CORPZ =EXPAPI*CORAZ+EXPBPI*CORBZ+EXPCPI*CORDZ
C
C                 *********************************************
C                 ***** Overlap Distribution Coefficients *****
C                 *********************************************
C
                  IF (INTTYP .NE. 96 ) THEN
                     CALL GETODC(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,
     &                          .FALSE.,.FALSE.,.FALSE.,EXPA,EXPB,
     &                          IPRINT,SAAB13,EXPPI,WORK,LWORK,CORPX,
     &                          CORPY,CORPZ,DONUC1,DOMOM1,ORIGIN,INTTYP)
                  ELSE
                     CALL GETODC(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,
     &                          .TRUE.,.FALSE.,.FALSE.,EXPA,EXPB,
     &                          IPRINT,SAAB13,EXPPI,WORK,LWORK,CORPX,
     &                          CORPY,CORPZ,DONUC1,DOMOM1,ORIGIN,INTTYP)
                  END IF
C
                  DIFPDX = CORPX - CORDX
                  DIFPDY = CORPY - CORDY
                  DIFPDZ = CORPZ - CORDZ
                  CALL DSOODC(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,EXPPI,
     &                        DIFPDX,DIFPDY,DIFPDZ,IPRINT,INTTYP)
C
C                 ********************************************
C                 **** Calculation of Hermitian integrals ****
C                 ********************************************
C
                  DO 700 I = 1, JATOMC
                     FACTOR = FACINT(I)
                     DIFCPX = CORCX(I) - CORPX
                     DIFCPY = CORCY(I) - CORPY
                     DIFCPZ = CORCZ(I) - CORPZ
                     IADR = 1 + (I - 1)*NAHGTF
                     CALL HERNAI(AHGTF,JMAX,EXPP,DIFCPX,DIFCPY,DIFCPZ,
     &                           FACTOR,IADR,ISTEPU,ISTEPV,NAHGTF,
     &                           IPRINT)
  700             CONTINUE
C
C                 ********************************************
C                 **** Calculation of Cartesian integrals ****
C                 ********************************************
C
                  IF (INTTYP .EQ. 12 .OR. INTTYP .EQ. 58) THEN
                     CALL DSOINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,
     &                           AHGTF,SHLINT(1,ISTDSO),NOPTYP,JATOMC)
                  ELSE IF (INTTYP .EQ. 96) THEN
                     CALL DSOKIN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,
     &                           AHGTF,SHLINT(1,ISTDSO),NOPTYP,JATOMC)
                  END IF
  600          CONTINUE
               ISTDSO = ISTDSO + 9*JATOMC
            END IF
  300       CONTINUE
         END IF
  200    CONTINUE
C
C     End loop over primitives
C
  100 CONTINUE
C
C     Print
C
      IF (IPRINT .GT. 30) THEN
         DO 900 I = 1, NOPTYP
            CALL AROUND('SHLINT for '//LABINT(I)//' in PR1PRD')
            CALL OUTPUT(SHLINT(1,I),1,KCKTA,1,KCKTB,KCKTA,KCKTB,1,LUPRI)
  900    CONTINUE
      END IF
      RETURN
      END
C  /* Deck ssoint */
      SUBROUTINE SSOINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT)
C
C     Spatial one-electron spin-orbit integrals
C
C     tuh 23 Nov 1989
C
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          AHGTF(*), SHLINT(KCKTAB,3)
#include "onecom.h"
#include "lmns.h"
      INT = 1
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         MAXT = LVALA + LVALB + 1
         MAXU = MVALA + MVALB + 1
         MAXV = NVALA + NVALB + 1
         IADRAV = 1
         DO 200 IV = 0,MAXV
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            FV = ODC(NVALA,NVALB,IV,1,0,3)
            IADRAU = IADRAV
            DO 300 IU = 0,MAXU
               EU = ODC(MVALA,MVALB,IU,0,0,2)
               FU = ODC(MVALA,MVALB,IU,1,0,2)
               EE = EU*EV
               FE = FU*EV
               EF = EU*FV
               DO 400 IT = 0,MAXT
                  ET = ODC(LVALA,LVALB,IT,0,0,1)
                  FT = ODC(LVALA,LVALB,IT,1,0,1)
                  FEE = FT*EE
                  EFE = ET*FE
                  EEF = ET*EF
                  AH0T = AHGTF(IADRAU + IT + 1)
                  AH0U = AHGTF(IADRAU + IT + ISTEPU)
                  AH0V = AHGTF(IADRAU + IT + ISTEPV)
                  SHLINT(INT,1) = SHLINT(INT,1) + EFE*AH0V - EEF*AH0U
                  SHLINT(INT,2) = SHLINT(INT,2) + EEF*AH0T - FEE*AH0V
                  SHLINT(INT,3) = SHLINT(INT,3) + FEE*AH0U - EFE*AH0T
  400          CONTINUE
               IADRAU = IADRAU + ISTEPU
  300       CONTINUE
            IADRAV = IADRAV + ISTEPV
  200    CONTINUE
         INT = INT + 1
  100 CONTINUE
      RETURN
      END
C  /* Deck ovlint */
      SUBROUTINE OVLINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,HUCFAC,
     &                  SHGTF,SHLINT,INTTYP)
C
C     tuh 1989
C     revised Jan 2005 hjaaj for NOPTYP=2 for HUCKEL
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,2)
#include "onecom.h"
#include "lmns.h"
      LOGICAL DOHUCMAT
C
      DOHUCMAT = INTTYP .LT. 0
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
C
         SHLINT(INT,1) = SHLINT(INT,1) - SX0*SY0*SZ0
         IF (DOHUCMAT)
     &   SHLINT(INT,2) = SHLINT(INT,2) - SX0*SY0*SZ0*HUCFAC
  100 CONTINUE
      RETURN
      END
C  /* Deck popovlint */
      SUBROUTINE POPOVLINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,
     &                     EXPBPI,SHGTF,SHLINT,INTTYP)
C
C     Mar 2006 hjaaj -- POPOVLP type for population analysis
C     Based on subr. OVLINT
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,2)
#include <onecom.h>
#include <lmns.h>
C
      FAC = 2*EXPBPI
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
C
         SHLINT(INT,1) = SHLINT(INT,1) - SX0*SY0*SZ0*FAC
  100 CONTINUE
      RETURN
      END
C  /* Deck dplint */
      SUBROUTINE DPLINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,ORIGIN)
C
C     tuh 1989
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,3), ORIGIN(3)
#include "onecom.h"
#include "lmns.h"
      POX = CORPX - ORIGIN(1)
      POY = CORPY - ORIGIN(2)
      POZ = CORPZ - ORIGIN(3)
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         DX0 = POX*SX0
         DY0 = POY*SY0
         DZ0 = POZ*SZ0
         IF (LVALA+LVALB.GT.0) DX0 = DX0+ SHGTF*ODC(LVALA,LVALB,1,0,0,1)
         IF (MVALA+MVALB.GT.0) DY0 = DY0+ SHGTF*ODC(MVALA,MVALB,1,0,0,2)
         IF (NVALA+NVALB.GT.0) DZ0 = DZ0+ SHGTF*ODC(NVALA,NVALB,1,0,0,3)
C
         INT = INT + 1
         SHLINT(INT,1) = SHLINT(INT,1) - DX0*SY0*SZ0
         SHLINT(INT,2) = SHLINT(INT,2) - SX0*DY0*SZ0
         SHLINT(INT,3) = SHLINT(INT,3) - SX0*SY0*DZ0
  100 CONTINUE
      RETURN
      END
C  /* Deck prpint */
      SUBROUTINE PRPINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,ORIGIN)
C
C     Computation of pRp integrals: 1/4 < grad A | {x,y,z} | grad B >
C     DPT correction to dipole integrals (same as PXPINT)
C
C     (WK/UniKA/10-03-2004)
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D25 = 0.25D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,3), ORIGIN(3)
#include "onecom.h"
#include "lmns.h"
      POX = CORPX - ORIGIN(1)
      POY = CORPY - ORIGIN(2)
      POZ = CORPZ - ORIGIN(3)
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         SX2 = SHGTF*ODC(LVALA,LVALB,0,1,1,1)
         SY2 = SHGTF*ODC(MVALA,MVALB,0,1,1,2)
         SZ2 = SHGTF*ODC(NVALA,NVALB,0,1,1,3)
         DX0 = POX*SX0 + SHGTF*ODC(LVALA,LVALB,1,0,0,1)
         DY0 = POY*SY0 + SHGTF*ODC(MVALA,MVALB,1,0,0,2)
         DZ0 = POZ*SZ0 + SHGTF*ODC(NVALA,NVALB,1,0,0,3)
         DX2 = POX*SX2 + SHGTF*ODC(LVALA,LVALB,1,1,1,1)
         DY2 = POY*SY2 + SHGTF*ODC(MVALA,MVALB,1,1,1,2)
         DZ2 = POZ*SZ2 + SHGTF*ODC(NVALA,NVALB,1,1,1,3)
         INT = INT + 1
         XXX = (DX2*SY0*SZ0 + DX0*SY2*SZ0 + DX0*SY0*SZ2) * D25
         YYY = (SX2*DY0*SZ0 + SX0*DY2*SZ0 + SX0*DY0*SZ2) * D25
         ZZZ = (SX2*SY0*DZ0 + SX0*SY2*DZ0 + SX0*SY0*DZ2) * D25
         SHLINT(INT,1) = SHLINT(INT,1) - XXX
         SHLINT(INT,2) = SHLINT(INT,2) - YYY
         SHLINT(INT,3) = SHLINT(INT,3) - ZZZ
  100 CONTINUE
      RETURN
      END
C  /* Deck pxpint */
      SUBROUTINE PXPINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,ORIGIN)
C
C     Computation of pRp integrals: < grad A | {x,y,z} | grad B >
C     DPT correction to dipole integrals (same as PRPINT)
C
C     (WK/UniKA/09-03-2004)
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D25 = 0.25D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,3), ORIGIN(3)
#include "onecom.h"
#include "lmns.h"
      POX = CORPX - ORIGIN(1)
      POY = CORPY - ORIGIN(2)
      POZ = CORPZ - ORIGIN(3)
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SX1 = SHGTF*ODC(LVALA,LVALB,0,1,0,1)
         SX2 = SHGTF*ODC(LVALA,LVALB,0,2,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SY1 = SHGTF*ODC(MVALA,MVALB,0,1,0,2)
         SY2 = SHGTF*ODC(MVALA,MVALB,0,2,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         SZ1 = SHGTF*ODC(NVALA,NVALB,0,1,0,3)
         SZ2 = SHGTF*ODC(NVALA,NVALB,0,2,0,3)
         DX0 = POX*SX0 + SHGTF*ODC(LVALA,LVALB,1,0,0,1)
         DY0 = POY*SY0 + SHGTF*ODC(MVALA,MVALB,1,0,0,2)
         DZ0 = POZ*SZ0 + SHGTF*ODC(NVALA,NVALB,1,0,0,3)
         DX2 = POX*SX2 + SHGTF*ODC(LVALA,LVALB,1,2,0,1) - SX1
         DY2 = POY*SY2 + SHGTF*ODC(MVALA,MVALB,1,2,0,2) - SY1
         DZ2 = POZ*SZ2 + SHGTF*ODC(NVALA,NVALB,1,2,0,3) - SZ1
         INT = INT + 1
         XXX = (DX2*SY0*SZ0 + DX0*SY2*SZ0 + DX0*SY0*SZ2) * D25
         YYY = (SX2*DY0*SZ0 + SX0*DY2*SZ0 + SX0*DY0*SZ2) * D25
         ZZZ = (SX2*SY0*DZ0 + SX0*SY2*DZ0 + SX0*SY0*DZ2) * D25
         SHLINT(INT,1) = SHLINT(INT,1) + XXX
         SHLINT(INT,2) = SHLINT(INT,2) + YYY
         SHLINT(INT,3) = SHLINT(INT,3) + ZZZ
  100 CONTINUE
      RETURN
      END
C  /* Deck kinint */
      SUBROUTINE KININT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT)
C
C     K.Ruud, Nov. 1991. Based on TUH's subroutine CINT0, but modified
C     for use in HERPRO.
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (DP5 = 0.5 D00)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB)
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
C    **********************************************
C    ***** CALCULATE KINETIC ENERGY INTEGRALS *****
C    **********************************************
C
         X0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         Y0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         Z0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         X2 = SHGTF*ODC(LVALA,LVALB,0,2,0,1)
         Y2 = SHGTF*ODC(MVALA,MVALB,0,2,0,2)
         Z2 = SHGTF*ODC(NVALA,NVALB,0,2,0,3)
         INT = INT + 1
         SHLINT(INT) = SHLINT(INT)
     &               + DP5*(X2*Y0*Z0 + X0*Y2*Z0 + X0*Y0*Z2)
 100  CONTINUE
      RETURN
      END
C  /* Deck mvlint */
      SUBROUTINE MVLINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT)
C
C     Sheela Kirpekar  jan. 93 (based on KININT)
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB)
#include "onecom.h"
#include "lmns.h"
#include "codata.h"
      PARAMETER (ALPMVL = ALPHA2*0.125D0, ALPMV2 = ALPHA2*0.25D0)
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
C    **********************************************
C    ***** CALCULATE MASS VELOCITY INTEGRALS ******
C    **********************************************
C
         X0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         Y0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         Z0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         X2 = SHGTF*ODC(LVALA,LVALB,0,2,0,1)
         Y2 = SHGTF*ODC(MVALA,MVALB,0,2,0,2)
         Z2 = SHGTF*ODC(NVALA,NVALB,0,2,0,3)
         X4 = SHGTF*ODC(LVALA,LVALB,0,4,0,1)
         Y4 = SHGTF*ODC(MVALA,MVALB,0,4,0,2)
         Z4 = SHGTF*ODC(NVALA,NVALB,0,4,0,3)
         INT = INT + 1
         SHLINT(INT) = SHLINT(INT)
     &               + ALPMVL*(X4*Y0*Z0 + X0*Y4*Z0 + X0*Y0*Z4)
     &               + ALPMV2*(X2*Y2*Z0 + X2*Y0*Z2 + X0*Y2*Z2)
 100  CONTINUE
      RETURN
      END
C  /* Deck hdbint */
      SUBROUTINE HDBINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  DIFOBX,DIFOBY,DIFOBZ)
C
C     K.Ruud, Aug. 1992
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (DP5 = 0.5D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,9)
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1, KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         SX1 = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         SY1 = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         SZ1 = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
         DX0 = SHGTF*ODC(LVALA,LVALB,0,1,0,1)
         DY0 = SHGTF*ODC(MVALA,MVALB,0,1,0,2)
         DZ0 = SHGTF*ODC(NVALA,NVALB,0,1,0,3)
         DX1 = SHGTF*ODC(LVALA,LVALB,0,1,1,1)
         DY1 = SHGTF*ODC(MVALA,MVALB,0,1,1,2)
         DZ1 = SHGTF*ODC(NVALA,NVALB,0,1,1,3)
         INT = INT + 1
C
         SHLINT(INT,1) = SHLINT(INT,1) - DP5*(DIFOBY*SY0*SZ1
     &                                 -      DIFOBZ*SY1*SZ0)*DX0
         SHLINT(INT,2) = SHLINT(INT,2) - DP5*(DIFOBZ*DX1*SZ0
     &                                 -      DIFOBX*DX0*SZ1)*SY0
         SHLINT(INT,3) = SHLINT(INT,3) - DP5*(DIFOBX*DX0*SY1
     &                                 -      DIFOBY*DX1*SY0)*SZ0
         SHLINT(INT,4) = SHLINT(INT,4) - DP5*(DIFOBY*DY0*SZ1
     &                                 -      DIFOBZ*DY1*SZ0)*SX0
         SHLINT(INT,5) = SHLINT(INT,5) - DP5*(DIFOBZ*SX1*SZ0
     &                                 -      DIFOBX*SX0*SZ1)*DY0
         SHLINT(INT,6) = SHLINT(INT,6) - DP5*(DIFOBX*SX0*DY1
     &                                 -      DIFOBY*SX1*DY0)*SZ0
         SHLINT(INT,7) = SHLINT(INT,7) - DP5*(DIFOBY*SY0*DZ1
     &                                 -      DIFOBZ*SY1*DZ0)*SX0
         SHLINT(INT,8) = SHLINT(INT,8) - DP5*(DIFOBZ*SX1*DZ0
     &                                 -      DIFOBX*SX0*DZ1)*SY0
         SHLINT(INT,9) = SHLINT(INT,9) - DP5*(DIFOBX*SX0*SY1
     &                                 -      DIFOBY*SX1*SY0)*DZ0
 100  CONTINUE
      RETURN
      END
C  /* Deck dpgint */
      SUBROUTINE DPGINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT)
C
C     K.Ruud, July-97
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,9)
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1, KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         INT = INT + 1
C
         SX1 = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         SY1 = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         SZ1 = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
         DX0 = SHGTF*ODC(LVALA,LVALB,0,1,0,1)
         DY0 = SHGTF*ODC(MVALA,MVALB,0,1,0,2)
         DZ0 = SHGTF*ODC(NVALA,NVALB,0,1,0,3)
         DX1 = SHGTF*ODC(LVALA,LVALB,0,1,1,1)
         DY1 = SHGTF*ODC(MVALA,MVALB,0,1,1,2)
         DZ1 = SHGTF*ODC(NVALA,NVALB,0,1,1,3)
C
         SHLINT(INT,1) = SHLINT(INT,1) - DX1*SY0*SZ0
         SHLINT(INT,2) = SHLINT(INT,2) - DX0*SY1*SZ0
         SHLINT(INT,3) = SHLINT(INT,3) - DX0*SY0*SZ1
         SHLINT(INT,4) = SHLINT(INT,4) - SX1*DY0*SZ0
         SHLINT(INT,5) = SHLINT(INT,5) - SX0*DY1*SZ0
         SHLINT(INT,6) = SHLINT(INT,6) - SX0*DY0*SZ1
         SHLINT(INT,7) = SHLINT(INT,7) - SX1*SY0*DZ0
         SHLINT(INT,8) = SHLINT(INT,8) - SX0*SY1*DZ0
         SHLINT(INT,9) = SHLINT(INT,9) - SX0*SY0*DZ1
 100  CONTINUE
      RETURN
      END
C  /* Deck dphint */
      SUBROUTINE DPHINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT)
C
C     K.Ruud, July-07
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,18)
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1, KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         SX1 = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         SY1 = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         SZ1 = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
         DX1 = SHGTF*ODC(LVALA,LVALB,0,1,0,1)
         DY1 = SHGTF*ODC(MVALA,MVALB,0,1,0,2)
         DZ1 = SHGTF*ODC(NVALA,NVALB,0,1,0,3)
         DX2 = SHGTF*ODC(LVALA,LVALB,0,2,0,1)
         DY2 = SHGTF*ODC(MVALA,MVALB,0,2,0,2)
         DZ2 = SHGTF*ODC(NVALA,NVALB,0,2,0,3)
         EX2 = SHGTF*ODC(LVALA,LVALB,0,1,1,1)
         EY2 = SHGTF*ODC(MVALA,MVALB,0,1,1,2)
         EZ2 = SHGTF*ODC(NVALA,NVALB,0,1,1,3)
         EX3 = SHGTF*ODC(LVALA,LVALB,0,2,1,1)
         EY3 = SHGTF*ODC(MVALA,MVALB,0,2,1,2)
         EZ3 = SHGTF*ODC(NVALA,NVALB,0,2,1,3)
C
         INT = INT + 1
C
         SHLINT(INT, 1) = SHLINT(INT, 1) - EX3*SY0*SZ0
         SHLINT(INT, 2) = SHLINT(INT, 2) - DX2*SY1*SZ0
         SHLINT(INT, 3) = SHLINT(INT, 3) - DX2*SY0*SZ1
         SHLINT(INT, 4) = SHLINT(INT, 4) - EX2*DY1*SZ0
         SHLINT(INT, 5) = SHLINT(INT, 5) - DX1*EY2*SZ0
         SHLINT(INT, 6) = SHLINT(INT, 6) - DX1*DY1*SZ1
         SHLINT(INT, 7) = SHLINT(INT, 7) - EX2*SY0*DZ1
         SHLINT(INT, 8) = SHLINT(INT, 8) - DX1*SY1*DZ1
         SHLINT(INT, 9) = SHLINT(INT, 9) - DX1*SY0*EZ2
         SHLINT(INT,10) = SHLINT(INT,10) - SX1*DY2*SZ0
         SHLINT(INT,11) = SHLINT(INT,11) - SX0*EY3*SZ0
         SHLINT(INT,12) = SHLINT(INT,12) - SX0*DY2*SZ1
         SHLINT(INT,13) = SHLINT(INT,13) - SX1*DY1*DZ1
         SHLINT(INT,14) = SHLINT(INT,14) - SX0*EY2*DZ1
         SHLINT(INT,15) = SHLINT(INT,15) - SX0*DY1*EZ2
         SHLINT(INT,16) = SHLINT(INT,16) - SX1*SY0*DZ2
         SHLINT(INT,17) = SHLINT(INT,17) - SX0*SY1*DZ2
         SHLINT(INT,18) = SHLINT(INT,18) - SX0*SY0*EZ3
 100  CONTINUE
      RETURN
      END
C  /* Deck qugint */
      SUBROUTINE QUGINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT)
C
C     K.Ruud, October/November-97
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,18)
#include <onecom.h>
#include <lmns.h>
      INT = 0
      DO 100 ICOMPA = 1, KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         SX1 = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         SY1 = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         SZ1 = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
         SX2 = SHGTF*ODC(LVALA,LVALB,0,0,2,1)
         SY2 = SHGTF*ODC(MVALA,MVALB,0,0,2,2)
         SZ2 = SHGTF*ODC(NVALA,NVALB,0,0,2,3)
         DX0 = SHGTF*ODC(LVALA,LVALB,0,1,0,1)
         DY0 = SHGTF*ODC(MVALA,MVALB,0,1,0,2)
         DZ0 = SHGTF*ODC(NVALA,NVALB,0,1,0,3)
         DX1 = SHGTF*ODC(LVALA,LVALB,0,1,1,1)
         DY1 = SHGTF*ODC(MVALA,MVALB,0,1,1,2)
         DZ1 = SHGTF*ODC(NVALA,NVALB,0,1,1,3)
         DX2 = SHGTF*ODC(LVALA,LVALB,0,1,2,1)
         DY2 = SHGTF*ODC(MVALA,MVALB,0,1,2,2)
         DZ2 = SHGTF*ODC(NVALA,NVALB,0,1,2,3)
         INT = INT + 1
C
         SHLINT(INT,1) = SHLINT(INT,1) - DX2*SY0*SZ0
         SHLINT(INT,2) = SHLINT(INT,2) - DX1*SY1*SZ0
         SHLINT(INT,3) = SHLINT(INT,3) - DX1*SY0*SZ1
         SHLINT(INT,4) = SHLINT(INT,4) - DX0*SY2*SZ0
         SHLINT(INT,5) = SHLINT(INT,5) - DX0*SY1*SZ1
         SHLINT(INT,6) = SHLINT(INT,6) - DX0*SY0*SZ2
         SHLINT(INT,7) = SHLINT(INT,7) - SX2*DY0*SZ0
         SHLINT(INT,8) = SHLINT(INT,8) - SX1*DY1*SZ0
         SHLINT(INT,9) = SHLINT(INT,9) - SX1*DY0*SZ1
         SHLINT(INT,10) = SHLINT(INT,10) - SX0*DY2*SZ0
         SHLINT(INT,11) = SHLINT(INT,11) - SX0*DY1*SZ1
         SHLINT(INT,12) = SHLINT(INT,12) - SX0*DY0*SZ2
         SHLINT(INT,13) = SHLINT(INT,13) - SX2*SY0*DZ0
         SHLINT(INT,14) = SHLINT(INT,14) - SX1*SY1*DZ0
         SHLINT(INT,15) = SHLINT(INT,15) - SX1*SY0*DZ1
         SHLINT(INT,16) = SHLINT(INT,16) - SX0*SY2*DZ0
         SHLINT(INT,17) = SHLINT(INT,17) - SX0*SY1*DZ1
         SHLINT(INT,18) = SHLINT(INT,18) - SX0*SY0*DZ2
 100  CONTINUE
      RETURN
      END
C  /* Deck ocgint */
      SUBROUTINE OCGINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT)
C
C     K.Ruud, November-97
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,30)
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1, KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         SX1 = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         SY1 = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         SZ1 = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
         SX2 = SHGTF*ODC(LVALA,LVALB,0,0,2,1)
         SY2 = SHGTF*ODC(MVALA,MVALB,0,0,2,2)
         SZ2 = SHGTF*ODC(NVALA,NVALB,0,0,2,3)
         SX3 = SHGTF*ODC(LVALA,LVALB,0,0,3,1)
         SY3 = SHGTF*ODC(MVALA,MVALB,0,0,3,2)
         SZ3 = SHGTF*ODC(NVALA,NVALB,0,0,3,3)
         DX0 = SHGTF*ODC(LVALA,LVALB,0,1,0,1)
         DY0 = SHGTF*ODC(MVALA,MVALB,0,1,0,2)
         DZ0 = SHGTF*ODC(NVALA,NVALB,0,1,0,3)
         DX1 = SHGTF*ODC(LVALA,LVALB,0,1,1,1)
         DY1 = SHGTF*ODC(MVALA,MVALB,0,1,1,2)
         DZ1 = SHGTF*ODC(NVALA,NVALB,0,1,1,3)
         DX2 = SHGTF*ODC(LVALA,LVALB,0,1,2,1)
         DY2 = SHGTF*ODC(MVALA,MVALB,0,1,2,2)
         DZ2 = SHGTF*ODC(NVALA,NVALB,0,1,2,3)
         DX3 = SHGTF*ODC(LVALA,LVALB,0,1,3,1)
         DY3 = SHGTF*ODC(MVALA,MVALB,0,1,3,2)
         DZ3 = SHGTF*ODC(NVALA,NVALB,0,1,3,3)
         INT = INT + 1
C
         SHLINT(INT,1)  = SHLINT(INT,1)  - DX3*SY0*SZ0
         SHLINT(INT,2)  = SHLINT(INT,2)  - DX2*SY1*SZ0
         SHLINT(INT,3)  = SHLINT(INT,3)  - DX2*SY0*SZ1
         SHLINT(INT,4)  = SHLINT(INT,4)  - DX1*SY2*SZ0
         SHLINT(INT,5)  = SHLINT(INT,5)  - DX1*SY1*SZ1
         SHLINT(INT,6)  = SHLINT(INT,6)  - DX1*SY0*SZ2
         SHLINT(INT,7)  = SHLINT(INT,7)  - DX0*SY3*SZ0
         SHLINT(INT,8)  = SHLINT(INT,8)  - DX0*SY2*SZ1
         SHLINT(INT,9)  = SHLINT(INT,9)  - DX0*SY1*SZ2
         SHLINT(INT,10) = SHLINT(INT,10) - DX0*SY0*SZ3
         SHLINT(INT,11) = SHLINT(INT,11) - SX3*DY0*SZ0
         SHLINT(INT,12) = SHLINT(INT,12) - SX2*DY1*SZ0
         SHLINT(INT,13) = SHLINT(INT,13) - SX2*DY0*SZ1
         SHLINT(INT,14) = SHLINT(INT,14) - SX1*DY2*SZ0
         SHLINT(INT,15) = SHLINT(INT,15) - SX1*DY1*SZ1
         SHLINT(INT,16) = SHLINT(INT,16) - SX1*DY0*SZ2
         SHLINT(INT,17) = SHLINT(INT,17) - SX0*DY3*SZ0
         SHLINT(INT,18) = SHLINT(INT,18) - SX0*DY2*SZ1
         SHLINT(INT,19) = SHLINT(INT,19) - SX0*DY1*SZ2
         SHLINT(INT,20) = SHLINT(INT,20) - SX0*DY0*SZ3
         SHLINT(INT,21) = SHLINT(INT,21) - SX3*SY0*DZ0
         SHLINT(INT,22) = SHLINT(INT,22) - SX2*SY1*DZ0
         SHLINT(INT,23) = SHLINT(INT,23) - SX2*SY0*DZ1
         SHLINT(INT,24) = SHLINT(INT,24) - SX1*SY2*DZ0
         SHLINT(INT,25) = SHLINT(INT,25) - SX1*SY1*DZ1
         SHLINT(INT,26) = SHLINT(INT,26) - SX1*SY0*DZ2
         SHLINT(INT,27) = SHLINT(INT,27) - SX0*SY3*DZ0
         SHLINT(INT,28) = SHLINT(INT,28) - SX0*SY2*DZ1
         SHLINT(INT,29) = SHLINT(INT,29) - SX0*SY1*DZ2
         SHLINT(INT,30) = SHLINT(INT,30) - SX0*SY0*DZ3
 100  CONTINUE
      RETURN
      END
C  /* Deck sm1int */
      SUBROUTINE SM1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  DIFABX,DIFABY,DIFABZ)
C
C     Kenneth Ruud's first subroutine. Modified Dec. 1991, KR
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D2I = 0.5D00)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,3)
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         DX0 = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         DY0 = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         DZ0 = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
C
         DPLX = DX0*SY0*SZ0
         DPLY = SX0*DY0*SZ0
         DPLZ = SX0*SY0*DZ0
C
         SHLINT(INT,1) = SHLINT(INT,1) - D2I*(DIFABY*DPLZ - DIFABZ*DPLY)
         SHLINT(INT,2) = SHLINT(INT,2) - D2I*(DIFABZ*DPLX - DIFABX*DPLZ)
         SHLINT(INT,3) = SHLINT(INT,3) - D2I*(DIFABX*DPLY - DIFABY*DPLX)
 100  CONTINUE
      RETURN
      END
C  /* Deck sm2int */
      SUBROUTINE SM2INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  EXPPI,DIFABX,DIFABY,DIFABZ)
C
C     K. Ruud, Oct. 1991
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0, D4INV = 0.25D0)
      DIMENSION SHLINT(KCKTAB,6)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1, KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1,KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         DX1 = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         DY1 = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         DZ1 = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
         DX2 = SHGTF*ODC(LVALA,LVALB,0,0,2,1)
         DY2 = SHGTF*ODC(MVALA,MVALB,0,0,2,2)
         DZ2 = SHGTF*ODC(NVALA,NVALB,0,0,2,3)
C
         SHLINT(INT,1) = SHLINT(INT,1) - (D2*DIFABZ*DIFABY*DY1*DZ1*SX0 -
     &                   DIFABZ*DIFABZ*DY2*SX0*SZ0 -
     &                   DIFABY*DIFABY*DZ2*SX0*SY0)*D4INV
         SHLINT(INT,2) = SHLINT(INT,2) - (DIFABZ*DIFABZ*DX1*DY1*SZ0 -
     &                   DIFABY*DIFABZ*DZ1*DX1*SY0 -
     &                   DIFABZ*DIFABX*DY1*DZ1*SX0 +
     &                   DIFABY*DIFABX*DZ2*SX0*SY0)*D4INV
         SHLINT(INT,3) = SHLINT(INT,3) - (DIFABY*DIFABY*DX1*DZ1*SY0 -
     &                   DIFABX*DIFABY*DY1*DZ1*SX0 +
     &                   DIFABZ*DIFABX*DY2*SX0*SZ0 -
     &                   DIFABY*DIFABZ*DY1*DX1*SZ0)*D4INV
         SHLINT(INT,4) = SHLINT(INT,4) - (D2*DIFABZ*DIFABX*DX1*DZ1*SY0 -
     &                   DIFABZ*DIFABZ*DX2*SY0*SZ0 -
     &                   DIFABX*DIFABX*DZ2*SX0*SY0)*D4INV
         SHLINT(INT,5) = SHLINT(INT,5) - (DIFABZ*DIFABY*DX2*SY0*SZ0 -
     &                   DIFABZ*DIFABX*DX1*DY1*SZ0 -
     &                   DIFABY*DIFABX*DX1*DZ1*SY0 +
     &                   DIFABX*DIFABX*DY1*DZ1*SX0)*D4INV
         SHLINT(INT,6) = SHLINT(INT,6) - (D2*DIFABX*DIFABY*DX1*DY1*SZ0 -
     &                   DIFABY*DIFABY*DX2*SY0*SZ0 -
     &                   DIFABX*DIFABX*DY2*SX0*SZ0)*D4INV
 100  CONTINUE
      RETURN
      END
C  /* Deck cm1int */
      SUBROUTINE CM1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  EXPPI,DIFABX,DIFABY,DIFABZ)
C
C     K. Ruud, Aug.-93.
C     Rewritten in order to separate dipole origin, aug.-94 KR
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "efield.h"
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0, D2I = 1.0D0/D2)
      DIMENSION SHLINT(KCKTAB,6)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
#include "onecom.h"
#include "orgcom.h"
#include "lmns.h"
      INT = 0
      ADX = CORAX - DIPORG(1)
      ADY = CORAY - DIPORG(2)
      ADZ = CORAZ - DIPORG(3)
      DO 100 ICOMPA = 1, KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
C
      DO 100 ICOMPB = 1,KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         SX1 = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         SY1 = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         SZ1 = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
         SX2 = SHGTF*ODC(LVALA,LVALB,0,0,2,1)
         SY2 = SHGTF*ODC(MVALA,MVALB,0,0,2,2)
         SZ2 = SHGTF*ODC(NVALA,NVALB,0,0,2,3)
         DX1 = SHGTF*(ODC(LVALA + 1,LVALB,0,0,0,1) +
     &           ADX* ODC(LVALA    ,LVALB,0,0,0,1))
         DY1 = SHGTF*(ODC(MVALA + 1,MVALB,0,0,0,2) +
     &           ADY* ODC(MVALA    ,MVALB,0,0,0,2))
         DZ1 = SHGTF*(ODC(NVALA + 1,NVALB,0,0,0,3) +
     &           ADZ* ODC(NVALA    ,NVALB,0,0,0,3))
         DX2 = SHGTF*(ODC(LVALA + 1,LVALB,0,0,1,1) +
     &           ADX* ODC(LVALA    ,LVALB,0,0,1,1))
         DY2 = SHGTF*(ODC(MVALA + 1,MVALB,0,0,1,2) +
     &           ADY* ODC(MVALA    ,MVALB,0,0,1,2))
         DZ2 = SHGTF*(ODC(NVALA + 1,NVALB,0,0,1,3) +
     &           ADZ* ODC(NVALA    ,NVALB,0,0,1,3))
C
         IF (FIELD1 .EQ. 'X-FIELD') THEN
            DPLX = DX2*SY0*SZ0
            DPLY = DX1*SY1*SZ0
            DPLZ = DX1*SY0*SZ1
         ELSE IF (FIELD1 .EQ. 'Y-FIELD') THEN
            DPLX = SX1*DY1*SZ0
            DPLY = SX0*DY2*SZ0
            DPLZ = SX0*DY1*SZ1
         ELSE
            DPLX = SX1*SY0*DZ1
            DPLY = SX0*SY1*DZ1
            DPLZ = SX0*SY0*DZ2
         END IF
C
         SHLINT(INT,1) = SHLINT(INT,1) - D2I*(DIFABY*DPLZ - DIFABZ*DPLY)
         SHLINT(INT,2) = SHLINT(INT,2) - D2I*(DIFABZ*DPLX - DIFABX*DPLZ)
         SHLINT(INT,3) = SHLINT(INT,3) - D2I*(DIFABX*DPLY - DIFABY*DPLX)
 100  CONTINUE
      RETURN
      END
C  /* Deck cm2int */
      SUBROUTINE CM2INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  EXPPI,DIFABX,DIFABY,DIFABZ)
C
C     K. Ruud, Aug.-93
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "efield.h"
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0, D4INV = 0.25D0)
      DIMENSION SHLINT(KCKTAB,6)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
#include "onecom.h"
#include "orgcom.h"
#include "lmns.h"
      INT = 0
      ADX = CORAX - DIPORG(1)
      ADY = CORAY - DIPORG(2)
      ADZ = CORAZ - DIPORG(3)
      DO 100 ICOMPA = 1, KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1,KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         SX1 = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         SY1 = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         SZ1 = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
         SX2 = SHGTF*ODC(LVALA,LVALB,0,0,2,1)
         SY2 = SHGTF*ODC(MVALA,MVALB,0,0,2,2)
         SZ2 = SHGTF*ODC(NVALA,NVALB,0,0,2,3)
C
         DX1 = SHGTF*(ODC(LVALA + 1,LVALB,0,0,0,1) +
     &           ADX* ODC(LVALA    ,LVALB,0,0,0,1))
         DY1 = SHGTF*(ODC(MVALA + 1,MVALB,0,0,0,2) +
     &           ADY* ODC(MVALA    ,MVALB,0,0,0,2))
         DZ1 = SHGTF*(ODC(NVALA + 1,NVALB,0,0,0,3) +
     &           ADZ* ODC(NVALA    ,NVALB,0,0,0,3))
         DX2 = SHGTF*(ODC(LVALA + 1,LVALB,0,0,1,1) +
     &           ADX* ODC(LVALA    ,LVALB,0,0,1,1))
         DY2 = SHGTF*(ODC(MVALA + 1,MVALB,0,0,1,2) +
     &           ADY* ODC(MVALA    ,MVALB,0,0,1,2))
         DZ2 = SHGTF*(ODC(NVALA + 1,NVALB,0,0,1,3) +
     &           ADZ* ODC(NVALA    ,NVALB,0,0,1,3))
         DX3 = SHGTF*(ODC(LVALA + 1,LVALB,0,0,2,1) +
     &           ADX* ODC(LVALA    ,LVALB,0,0,2,1))
         DY3 = SHGTF*(ODC(MVALA + 1,MVALB,0,0,2,2) +
     &           ADY* ODC(MVALA    ,MVALB,0,0,2,2))
         DZ3 = SHGTF*(ODC(NVALA + 1,NVALB,0,0,2,3) +
     &           ADZ* ODC(NVALA    ,NVALB,0,0,2,3))
C
         IF (FIELD2 .EQ. 'X-FIELD') THEN
            DPLXX = DX3*SY0*SZ0
            DPLXY = DX2*SY1*SZ0
            DPLXZ = DX2*SY0*SZ1
            DPLYY = DX1*SY2*SZ0
            DPLYZ = DX1*SY1*SZ1
            DPLZZ = DX1*SY0*SZ2
         ELSE IF (FIELD2 .EQ. 'Y-FIELD') THEN
            DPLXX = SX2*DY1*SZ0
            DPLXY = SX1*DY2*SZ0
            DPLXZ = SX1*DY1*SZ1
            DPLYY = SX0*DY3*SZ0
            DPLYZ = SX0*DY2*SZ1
            DPLZZ = SX0*DY1*SZ2
         ELSE
            DPLXX = SX2*SY0*DZ1
            DPLXY = SX1*SY1*DZ1
            DPLXZ = SX1*SY0*DZ2
            DPLYY = SX0*SY2*DZ1
            DPLYZ = SX0*SY1*DZ2
            DPLZZ = SX0*SY0*DZ3
         END IF
C
         SHLINT(INT,1) = SHLINT(INT,1) - (D2*DIFABZ*DIFABY*DPLYZ
     &                                 -  DIFABZ*DIFABZ*DPLYY
     &                                 - DIFABY*DIFABY*DPLZZ)*D4INV
         SHLINT(INT,2) = SHLINT(INT,2) - (DIFABZ*DIFABZ*DPLXY -
     &                   DIFABY*DIFABZ*DPLXZ -
     &                   DIFABZ*DIFABX*DPLYZ +
     &                   DIFABY*DIFABX*DPLZZ)*D4INV
         SHLINT(INT,3) = SHLINT(INT,3) - (DIFABY*DIFABY*DPLXZ -
     &                   DIFABX*DIFABY*DPLYZ +
     &                   DIFABZ*DIFABX*DPLYY -
     &                   DIFABY*DIFABZ*DPLXY)*D4INV
         SHLINT(INT,4) = SHLINT(INT,4) - (D2*DIFABZ*DIFABX*DPLXZ -
     &                   DIFABZ*DIFABZ*DPLXX -
     &                   DIFABX*DIFABX*DPLZZ)*D4INV
         SHLINT(INT,5) = SHLINT(INT,5) - (DIFABZ*DIFABY*DPLXX -
     &                   DIFABZ*DIFABX*DPLXY -
     &                   DIFABY*DIFABX*DPLXZ +
     &                   DIFABX*DIFABX*DPLYZ)*D4INV
         SHLINT(INT,6) = SHLINT(INT,6) - (D2*DIFABX*DIFABY*DPLXY -
     &                   DIFABY*DIFABY*DPLXX -
     &                   DIFABX*DIFABX*DPLYY)*D4INV
 100  CONTINUE
      RETURN
      END
C  /* Deck qdbint */
      SUBROUTINE QDBINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  EXPPI,DIFABX,DIFABY,DIFABZ)
C
C     K. Ruud, February.-02.
C     Based on CM1INT
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "efield.h"
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0, D2I = 1.0D0/D2, D1P5 = 1.5D0,
     &           D1 = 1.0D0)
      DIMENSION SHLINT(KCKTAB,6)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
#include "onecom.h"
#include "orgcom.h"
#include "lmns.h"
      INT = 0
      ADX = CORAX - DIPORG(1)
      ADY = CORAY - DIPORG(2)
      ADZ = CORAZ - DIPORG(3)
      DO 100 ICOMPA = 1, KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
C
      DO 100 ICOMPB = 1,KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         SX1 = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         SY1 = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         SZ1 = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
         DX1 = SHGTF*(ODC(LVALA + 1,LVALB,0,0,0,1) +
     &           ADX* ODC(LVALA    ,LVALB,0,0,0,1))
         DY1 = SHGTF*(ODC(MVALA + 1,MVALB,0,0,0,2) +
     &           ADY* ODC(MVALA    ,MVALB,0,0,0,2))
         DZ1 = SHGTF*(ODC(NVALA + 1,NVALB,0,0,0,3) +
     &           ADZ* ODC(NVALA    ,NVALB,0,0,0,3))
         DX2 = SHGTF*(ODC(LVALA + 2,LVALB,0,0,0,1) +
     &        D2*ADX* ODC(LVALA + 1,LVALB,0,0,0,1) +
     &        ADX*ADX*ODC(LVALA    ,LVALB,0,0,0,1))
         DY2 = SHGTF*(ODC(MVALA + 2,MVALB,0,0,0,2) +
     &        D2*ADY* ODC(MVALA + 1,MVALB,0,0,0,2) +
     &        ADY*ADY*ODC(MVALA    ,MVALB,0,0,0,2))
         DZ2 = SHGTF*(ODC(NVALA + 2,NVALB,0,0,0,3) +
     &        D2*ADZ* ODC(NVALA + 1,NVALB,0,0,0,3) +
     &        ADZ*ADZ*ODC(NVALA    ,NVALB,0,0,0,3))
         SDX2 = SHGTF*(ODC(LVALA + 1,LVALB,0,0,1,1) +
     &            ADX* ODC(LVALA    ,LVALB,0,0,1,1))
         SDY2 = SHGTF*(ODC(MVALA + 1,MVALB,0,0,1,2) +
     &            ADY* ODC(MVALA    ,MVALB,0,0,1,2))
         SDZ2 = SHGTF*(ODC(NVALA + 1,NVALB,0,0,1,3) +
     &            ADZ* ODC(NVALA    ,NVALB,0,0,1,3))
         DX3 = SHGTF*(ODC(LVALA + 2,LVALB,0,0,1,1) +
     &        D2*ADX* ODC(LVALA + 1,LVALB,0,0,1,1) +
     &        ADX*ADX*ODC(LVALA    ,LVALB,0,0,1,1))
         DY3 = SHGTF*(ODC(MVALA + 2,MVALB,0,0,1,2) +
     &        D2*ADY* ODC(MVALA + 1,MVALB,0,0,1,2) +
     &        ADY*ADY*ODC(MVALA    ,MVALB,0,0,1,2))
         DZ3 = SHGTF*(ODC(NVALA + 2,NVALB,0,0,1,3) +
     &        D2*ADZ* ODC(NVALA + 1,NVALB,0,0,1,3) +
     &        ADZ*ADZ*ODC(NVALA    ,NVALB,0,0,1,3))
C
         IF (FIELD3(4:7) .EQ. 'FGRD') THEN
            FAC = D1P5
         ELSE
            FAC = D1
         END IF
         IF (FIELD3(1:2) .EQ. 'XX' .OR. FIELD3(1:2) .EQ. 'YY' .OR.
     &       FIELD3(1:2) .EQ. 'ZZ') THEN
            IF (FIELD3(4:7) .EQ. 'FGRD') THEN
               XR2 = D2I*(DX3*SY0*SZ0 + SX1*DY2*SZ0 + SX1*SY0*DZ2)
               YR2 = D2I*(DX2*SY1*SZ0 + SX0*DY3*SZ0 + SX0*SY1*DZ2)
               ZR2 = D2I*(DX2*SY0*SZ1 + SX0*DY2*SZ1 + SX0*SY0*DZ3)
            ELSE
               XR2 = D0
               YR2 = D0
               ZR2 = D0
            END IF
            IF (FIELD3(1:2) .EQ. 'XX') THEN
               DPLX = - XR2 + FAC*DX3*SY0*SZ0
               DPLY = - YR2 + FAC*DX2*SY1*SZ0
               DPLZ = - ZR2 + FAC*DX2*SY0*SZ1
            ELSE IF (FIELD3(1:2) .EQ. 'YY') THEN
               DPLX = - XR2 + FAC*SX1*DY2*SZ0
               DPLY = - YR2 + FAC*SX0*DY3*SZ0
               DPLZ = - ZR2 + FAC*SX0*DY2*SZ1
            ELSE
               DPLX = - XR2 + FAC*SX1*SY0*DZ2
               DPLY = - YR2 + FAC*SX0*SY1*DZ2
               DPLZ = - ZR2 + FAC*SX0*SY0*DZ3
            END IF
         ELSE IF (FIELD3(1:2) .EQ. 'XY') THEN
            DPLX = FAC*SDX2*DY1*SZ0
            DPLY = FAC*DX1*SDY2*SZ0
            DPLZ = FAC*DX1*DY1*SZ1
         ELSE IF (FIELD3(1:2) .EQ. 'XZ') THEN
            DPLX = FAC*SDX2*SY0*SZ1
            DPLY = FAC*DX1*SY1*DZ1
            DPLZ = FAC*DX1*SY0*SDZ2
         ELSE IF (FIELD3(1:2) .EQ. 'YZ') THEN
            DPLX = FAC*SX1*DY1*DZ1
            DPLY = FAC*SX0*SDY2*DZ1
            DPLZ = FAC*SX0*SY1*SDZ2
         END IF
C
         SHLINT(INT,1) = SHLINT(INT,1) - D2I*(DIFABY*DPLZ - DIFABZ*DPLY)
         SHLINT(INT,2) = SHLINT(INT,2) - D2I*(DIFABZ*DPLX - DIFABX*DPLZ)
         SHLINT(INT,3) = SHLINT(INT,3) - D2I*(DIFABX*DPLY - DIFABY*DPLX)
 100  CONTINUE
      RETURN
      END
C  /* Deck am1int */
      SUBROUTINE AM1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  FAC)
C
C     K.Ruud, Oct 1991, Modified Dec. 1991, KR
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,3)
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1,KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         SX1 = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         SY1 = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         SZ1 = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
         DX0 = SHGTF*ODC(LVALA,LVALB,0,1,0,1)
         DY0 = SHGTF*ODC(MVALA,MVALB,0,1,0,2)
         DZ0 = SHGTF*ODC(NVALA,NVALB,0,1,0,3)
C
         SHLINT(INT,1) = SHLINT(INT,1) - (SY1*DZ0 - DY0*SZ1)*SX0*FAC
         SHLINT(INT,2) = SHLINT(INT,2) - (SZ1*DX0 - SX1*DZ0)*SY0*FAC
         SHLINT(INT,3) = SHLINT(INT,3) - (SX1*DY0 - DX0*SY1)*SZ0*FAC
 100  CONTINUE
      RETURN
      END
C  /* Deck am1int */
      SUBROUTINE AMDINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  EXPA)
C
C     K.Ruud, Pisa July 2007
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,18)
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1,KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         SX1 = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         SY1 = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         SZ1 = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
         DX1 = SHGTF*ODC(LVALA,LVALB,0,1,0,1)
         DY1 = SHGTF*ODC(MVALA,MVALB,0,1,0,2)
         DZ1 = SHGTF*ODC(NVALA,NVALB,0,1,0,3)
         DX2 = SHGTF*ODC(LVALA,LVALB,0,2,0,1)
         DY2 = SHGTF*ODC(MVALA,MVALB,0,2,0,2)
         DZ2 = SHGTF*ODC(NVALA,NVALB,0,2,0,3)
         EX2 = SHGTF*ODC(LVALA,LVALB,0,1,1,1)
         EY2 = SHGTF*ODC(MVALA,MVALB,0,1,1,2)
         EZ2 = SHGTF*ODC(NVALA,NVALB,0,1,1,3)
C Diff.x direction on center 2
         SHLINT(INT,1) = SHLINT(INT,1) - (SY1*DZ1 - DY1*SZ1)*DX1
         SHLINT(INT,2) = SHLINT(INT,2) - (SZ1*DX2 - EX2*DZ1)*SY0
         SHLINT(INT,3) = SHLINT(INT,3) - (EX2*DY1 - DX2*SY1)*SZ0
C Diff.y direction on center 2
         SHLINT(INT,4) = SHLINT(INT,4) - (EY2*DZ1 - DY2*SZ1)*SX0
         SHLINT(INT,5) = SHLINT(INT,5) - (SZ1*DX1 - SX1*DZ1)*DY1
         SHLINT(INT,6) = SHLINT(INT,6) - (SX1*DY2 - DX1*EY2)*SZ0
C Diff.z direction on center 2
         SHLINT(INT,7) = SHLINT(INT,7) - (SY1*DZ2 - DY1*EZ2)*SX0
         SHLINT(INT,8) = SHLINT(INT,8) - (EZ2*DX1 - SX1*DZ2)*SY0
         SHLINT(INT,9) = SHLINT(INT,9) - (SX1*DY1 - DX1*SY1)*DZ1
 100  CONTINUE
      RETURN
      END
      SUBROUTINE OZKEIN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT)
C
C     K.Ruud and J.Vaara, Sep.01, based on AM1INT
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (DP25 = 0.25D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,3)
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1,KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         SX1 = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         SY1 = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         SZ1 = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
         DX0 = SHGTF*ODC(LVALA,LVALB,0,1,0,1)
         DY0 = SHGTF*ODC(MVALA,MVALB,0,1,0,2)
         DZ0 = SHGTF*ODC(NVALA,NVALB,0,1,0,3)
         DX1 = SHGTF*ODC(LVALA,LVALB,0,2,0,1)
         DY1 = SHGTF*ODC(MVALA,MVALB,0,2,0,2)
         DZ1 = SHGTF*ODC(NVALA,NVALB,0,2,0,3)
         DX2 = SHGTF*ODC(LVALA,LVALB,0,3,0,1)
         DY2 = SHGTF*ODC(MVALA,MVALB,0,3,0,2)
         DZ2 = SHGTF*ODC(NVALA,NVALB,0,3,0,3)
         SDX = SHGTF*ODC(LVALA,LVALB,0,2,1,1)
         SDY = SHGTF*ODC(MVALA,MVALB,0,2,1,2)
         SDZ = SHGTF*ODC(NVALA,NVALB,0,2,1,3)
C
         SHLINT(INT,1) = SHLINT(INT,1) - (DX1*SY1*DZ0 + SX0*SDY*DZ0
     &                 + SX0*SY1*DZ2 - DX1*DY0*SZ1 - SX0*DY2*SZ1
     &                 - SX0*DY0*SDZ)*DP25
         SHLINT(INT,2) = SHLINT(INT,2) - (DX2*SY0*SZ1 + DX0*DY1*SZ1
     &                 + DX0*SY0*SDZ - SDX*SY0*DZ0 - SX1*DY1*DZ0
     &                 - SX1*SY0*DZ2)*DP25
         SHLINT(INT,3) = SHLINT(INT,3) - (SDX*DY0*SZ0 + SX1*DY2*SZ0
     &                 + SX1*DY0*DZ1 - DX2*SY1*SZ0 - DX0*SDY*SZ0
     &                 - DX0*SY1*DZ1)*DP25
 100  CONTINUE
      RETURN
      END
C  /* Deck mm1int */
      SUBROUTINE MM1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  AHGTF,DIFABX,DIFABY,DIFABZ,NATOMC,INTTYP)
C
C     K.Ruud, Feb. 1992
C
C     INTTYP decides whether to calculate only the potential contribution
C     or also the kinetic energy contribution
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D4INV = 0.250D0, D2INV = 0.50D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,3), AHGTF(*)
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1,KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
         SX0  = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0  = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0  = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         SX1  = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         SY1  = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         SZ1  = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
         SX21 = SHGTF*ODC(LVALA,LVALB,0,2,1,1)
         SY21 = SHGTF*ODC(MVALA,MVALB,0,2,1,2)
         SZ21 = SHGTF*ODC(NVALA,NVALB,0,2,1,3)
         DX2  = SHGTF*ODC(LVALA,LVALB,0,2,0,1)
         DY2  = SHGTF*ODC(MVALA,MVALB,0,2,0,2)
         DZ2  = SHGTF*ODC(NVALA,NVALB,0,2,0,3)
C
C       Kinetic energy contribution
C
         IF (INTTYP .EQ. 19 .OR. INTTYP .EQ. 20) THEN
            SHLINT(INT,1) = SHLINT(INT,1)+D4INV*( DIFABY*(DX2*SY0*SZ1 +
     &                                       SX0*DY2*SZ1 + SX0*SY0*SZ21)
     &                                   - DIFABZ*(DX2*SY1*SZ0 +
     &                                     SX0*SY21*SZ0 + SX0*SY1*DZ2))
            SHLINT(INT,2) = SHLINT(INT,2)+D4INV*( DIFABZ*(SX21*SY0*SZ0 +
     &                                     SX1*DY2*SZ0 + SX1*SY0*DZ2)
     &                                   - DIFABX*(DX2*SY0*SZ1 +
     &                                     SX0*DY2*SZ1 + SX0*SY0*SZ21))
            SHLINT(INT,3) = SHLINT(INT,3)+D4INV*( DIFABX*(DX2*SY1*SZ0 +
     &                                     SX0*SY21*SZ0 + SX0*SY1*DZ2)
     &                                   - DIFABY*(SX21*SY0*SZ0 +
     &                                     SX1*DY2*SZ0 + SX1*SY0*DZ2))
         END IF
C
C     Nuclear attraction contribution
C
         MAXT = LVALA + LVALB
         MAXU = MVALA + MVALB
         MAXV = NVALA + NVALB
         IADRAV = 1
         DO 200 IV = 0, MAXV + 1
            IADRAU = IADRAV
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            FV = ODC(NVALA,NVALB,IV,0,1,3)
            IUMAX = MAXU + 1
            IF (IV .GT. MAXV) IUMAX = MAXU
            DO 300 IU = 0, IUMAX
               EU = ODC(MVALA,MVALB,IU,0,0,2)
               FU = ODC(MVALA,MVALB,IU,0,1,2)
               ITMAX = MAXT + 1
               IF ((IU .GT. MAXU) .OR. (IV .GT. MAXV)) ITMAX = MAXT
               DO 400 IT = 0, ITMAX
                  ET = ODC(LVALA,LVALB,IT,0,0,1)
                  FT = ODC(LVALA,LVALB,IT,0,1,1)
                  FX = FT*EU*EV
                  FY = ET*FU*EV
                  FZ = ET*EU*FV
                  ATUV = -D2INV*AHGTF(IADRAU + IT)
                  SHLINT(INT,1)=SHLINT(INT,1)+ATUV*(DIFABY*FZ-DIFABZ*FY)
                  SHLINT(INT,2)=SHLINT(INT,2)+ATUV*(DIFABZ*FX-DIFABX*FZ)
                  SHLINT(INT,3)=SHLINT(INT,3)+ATUV*(DIFABX*FY-DIFABY*FX)
 400           CONTINUE
               IADRAU = IADRAU + ISTEPU
 300        CONTINUE
            IADRAV = IADRAV + ISTEPV
 200     CONTINUE
 100  CONTINUE
      RETURN
      END
C  /* Deck pvpint */
      SUBROUTINE PVPINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT,
     &                  NATOMC,EXPB)
C
C     K.Ruud, September 2000
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB), AHGTF(*)
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1,KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
         MAXT = LVALA + LVALB
         MAXU = MVALA + MVALB
         MAXV = NVALA + NVALB
         IADRAV = 1
         DO 200 IV = 0, MAXV + 2
            IADRAU = IADRAV
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            FV = ODC(NVALA,NVALB,IV,1,1,3)
            IUMAX = MAXU + 2
            IF (IV .GT. MAXV + 1) THEN
               IUMAX = MAXU
            ELSE IF (IV .GT. MAXV) THEN
               IUMAX = MAXU + 1
            END IF
            DO 300 IU = 0, IUMAX
               EU = ODC(MVALA,MVALB,IU,0,0,2)
               FU = ODC(MVALA,MVALB,IU,1,1,2)
               ITMAX = MAXT + 2
               IF (((IU .GT. MAXU + 1) .OR. (IV .GT. MAXV + 1))
     &               .OR. ((IU .EQ. MAXU + 1) .AND. (IV .EQ. MAXV + 1)))
     &               THEN
                  ITMAX = MAXT
               ELSE IF ((IU .GT. MAXU) .OR. (IV .GT. MAXV)) THEN
                  ITMAX = MAXT + 1
               END IF
               DO 400 IT = 0, ITMAX
                  ET = ODC(LVALA,LVALB,IT,0,0,1)
                  FT = ODC(LVALA,LVALB,IT,1,1,1)
                  GXX = FT*EU*EV
                  GYY = ET*FU*EV
                  GZZ = ET*EU*FV
C
                  SHLINT(INT) = SHLINT(INT)
     &                        + (GXX + GYY + GZZ)*AHGTF(IADRAU+IT)
 400           CONTINUE
               IADRAU = IADRAU + ISTEPU
 300        CONTINUE
            IADRAV = IADRAV + ISTEPV
 200     CONTINUE
 100  CONTINUE
      RETURN
      END
C  /* Deck dsusan */
      SUBROUTINE DSUSAN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,
     &                  SHLINT,CORPX,CORPY,CORPZ,DIFABX,DIFABY,DIFABZ)
C
C     K.Ruud, Feb. 1992
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D2INV = 0.50D0, D4INV = 0.25D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,6), ORIGIN(3)
#include "onecom.h"
#include "lmns.h"
      INT = 0
      PBX = CORBX - CORPX
      PBY = CORBY - CORPY
      PBZ = CORBZ - CORPZ
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1,KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
         SX0  = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0  = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0  = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         SX1  = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         SY1  = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         SZ1  = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
         SXT0 = SHGTF*ODC(LVALA,LVALB,1,0,0,1)
         SYT0 = SHGTF*ODC(MVALA,MVALB,1,0,0,2)
         SZT0 = SHGTF*ODC(NVALA,NVALB,1,0,0,3)
         SXT1 = SHGTF*ODC(LVALA,LVALB,1,0,1,1)
         SYT1 = SHGTF*ODC(MVALA,MVALB,1,0,1,2)
         SZT1 = SHGTF*ODC(NVALA,NVALB,1,0,1,3)
         SX11 = SHGTF*ODC(LVALA,LVALB,0,1,1,1)
         SY11 = SHGTF*ODC(MVALA,MVALB,0,1,1,2)
         SZ11 = SHGTF*ODC(NVALA,NVALB,0,1,1,3)
         DX0  = SHGTF*ODC(LVALA,LVALB,0,1,0,1)
         DY0  = SHGTF*ODC(MVALA,MVALB,0,1,0,2)
         DZ0  = SHGTF*ODC(NVALA,NVALB,0,1,0,3)
C
         SX1N = SXT0 - PBX*SX0
         SY1N = SYT0 - PBY*SY0
         SZ1N = SZT0 - PBZ*SZ0
         SX2N = SXT1 - PBX*SX1
         SY2N = SYT1 - PBY*SY1
         SZ2N = SZT1 - PBZ*SZ1
C
         SHLINT(INT,1) = SHLINT(INT,1) + (DIFABY*SX0*
     &                                   (SZ11*SY1N - DY0*SZ2N )
     &                                 -  DIFABZ*SX0*
     &                                   (SY2N*DZ0  - SY11*SZ1N))*D2INV
         SHLINT(INT,2) = SHLINT(INT,2) + (DIFABY*SY0*
     &                                   (SZ2N*DX0  - SZ11*SX1N)
     &                                 - (DX0*SY1*SZ1N + SX1*DY0*SZ1N
     &                                 - (SX1*SY1N + SX1N*SY1)*DZ0)*
     &                                    DIFABZ
     &                                 -  DIFABX*SX0*
     &                                   (SZ11*SY1N - SZ2N*DY0))*D4INV
         SHLINT(INT,3) = SHLINT(INT,3) + (DIFABX*SX0*
     &                                   (SY2N*DZ0  - SY11*SZ1N )
     &                                 + (DY0*(SX1N*SZ1 + SX1*SZ1N)
     &                                 -  DX0*SY1N*SZ1 - SX1*SY1N*DZ0)*
     &                                    DIFABY
     &                                 -  DIFABZ*SZ0*
     &                                   (SY11*SX1N - DX0*SY2N))*D4INV
         SHLINT(INT,4) = SHLINT(INT,4) + (DIFABZ*SY0*
     &                                   (SZ1N*SX11  - SX2N*DZ0 )
     &                                 -  DIFABX*SY0*
     &                                   (SZ2N*DX0  - SX1N*SZ11))*D2INV
         SHLINT(INT,5) = SHLINT(INT,5) + (DIFABZ*SZ0*
     &                                   (SX2N*DY0  - SX11*SY1N)
     &                                 + (DX0*(SY1*SZ1N + SY1N*SZ1)
     &                                 -  SX1N*SY1*DZ0 - SX1N*DY0*SZ1)*
     &                                    DIFABX
     &                                 -  DIFABY*SY0*
     &                                   (SX11*SZ1N - SX2N*DZ0))*D4INV
         SHLINT(INT,6) = SHLINT(INT,6) + (DIFABX*SZ0*
     &                                   (SY11*SX1N - SY2N*DX0)
     &                                 -  DIFABY*SZ0*
     &                                   (SX2N*DY0  - SX11*SY1N))*D2INV
 100  CONTINUE
      RETURN
      END
C  /* Deck mqdpin */
      SUBROUTINE MQDPIN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT)
C
C     K.Ruud, Feb. 1992
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D3INV = 1.0D0/3.0D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,9)
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1,KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
         SX0  = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0  = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0  = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         SX1  = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         SY1  = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         SZ1  = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
         SX2  = SHGTF*ODC(LVALA,LVALB,0,0,2,1)
         SY2  = SHGTF*ODC(MVALA,MVALB,0,0,2,2)
         SZ2  = SHGTF*ODC(NVALA,NVALB,0,0,2,3)
         DX0  = SHGTF*ODC(LVALA,LVALB,0,1,0,1)
         DY0  = SHGTF*ODC(MVALA,MVALB,0,1,0,2)
         DZ0  = SHGTF*ODC(NVALA,NVALB,0,1,0,3)
         DX1  = SHGTF*ODC(LVALA,LVALB,0,1,1,1)
         DY1  = SHGTF*ODC(MVALA,MVALB,0,1,1,2)
         DZ1  = SHGTF*ODC(NVALA,NVALB,0,1,1,3)
C
         SHLINT(INT,1) = SHLINT(INT,1) - D3INV*SX1*(SY1*DZ0 - DY0*SZ1)
         SHLINT(INT,2) = SHLINT(INT,2) - D3INV*SY0*(DX1*SZ1 - SX2*DZ0)
         SHLINT(INT,3) = SHLINT(INT,3) - D3INV*SZ0*(SX2*DY0 - DX1*SY1)
         SHLINT(INT,4) = SHLINT(INT,4) - D3INV*SX0*(SY2*DZ0 - DY1*SZ1)
         SHLINT(INT,5) = SHLINT(INT,5) - D3INV*SY1*(DX0*SZ1 - SX1*DZ0)
         SHLINT(INT,6) = SHLINT(INT,6) - D3INV*SZ0*(SX1*DY1 - DX0*SY2)
         SHLINT(INT,7) = SHLINT(INT,7) - D3INV*SX0*(SY1*DZ1 - DY0*SZ2)
         SHLINT(INT,8) = SHLINT(INT,8) - D3INV*SY0*(DX0*SZ2 - SX1*DZ1)
         SHLINT(INT,9) = SHLINT(INT,9) - D3INV*SZ1*(SX1*DY0 - DX0*SY1)
 100  CONTINUE
      RETURN
      END
C  /* Deck dsuslh */
      SUBROUTINE DSUSLH(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  AHGTF,DIFABX,DIFABY,DIFABZ,NATOMC,INTTYP)
C
C     K.Ruud, Feb. 1992
C

#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D4INV = 0.250D0, D8INV = 0.1250D0, D2 = 2.0D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,6), AHGTF(*)
      INTEGER   INTTYP
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1,KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
C  .... when wished, calculate the kinetic energy contribution
C       Skip the Kinetic contribution if PCM type integral 67
C
       IF (INTTYP .NE. 67) THEN

         SX0  = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0  = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0  = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         SX1  = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         SY1  = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         SZ1  = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
         SX2  = SHGTF*ODC(LVALA,LVALB,0,0,2,1)
         SY2  = SHGTF*ODC(MVALA,MVALB,0,0,2,2)
         SZ2  = SHGTF*ODC(NVALA,NVALB,0,0,2,3)
         SX21 = SHGTF*ODC(LVALA,LVALB,0,2,1,1)
         SY21 = SHGTF*ODC(MVALA,MVALB,0,2,1,2)
         SZ21 = SHGTF*ODC(NVALA,NVALB,0,2,1,3)
         SX22 = SHGTF*ODC(LVALA,LVALB,0,2,2,1)
         SY22 = SHGTF*ODC(MVALA,MVALB,0,2,2,2)
         SZ22 = SHGTF*ODC(NVALA,NVALB,0,2,2,3)
         DX2  = SHGTF*ODC(LVALA,LVALB,0,2,0,1)
         DY2  = SHGTF*ODC(MVALA,MVALB,0,2,0,2)
         DZ2  = SHGTF*ODC(NVALA,NVALB,0,2,0,3)
         XX   = SX22*SY0*SZ0 + SX2*DY2*SZ0 + SX2*SY0*DZ2
         XY   = SX21*SY1*SZ0 + SX1*SY21*SZ0 + SX1*SY1*DZ2
         XZ   = SX21*SY0*SZ1 + SX1*DY2*SZ1 + SX1*SY0*SZ21
         YY   = DX2*SY2*SZ0 + SX0*SY22*SZ0 + SX0*SY2*DZ2
         YZ   = DX2*SY1*SZ1 + SX0*SY21*SZ1 + SX0*SY1*SZ21
         ZZ   = DX2*SY0*SZ2 + SX0*DY2*SZ2 + SX0*SY0*SZ22
C
C       Kinetic energy contribution
C
         SHLINT(INT,1) = SHLINT(INT,1) + (D2*DIFABY*DIFABZ*YZ
     &                                 -  DIFABY*DIFABY*ZZ
     &                                 -  DIFABZ*DIFABZ*YY)*D8INV
         SHLINT(INT,2) = SHLINT(INT,2) + (DIFABZ*DIFABZ*XY
     &                                 +  DIFABY*DIFABX*ZZ
     &                                 -  DIFABZ*DIFABX*YZ
     &                                 -  DIFABY*DIFABZ*XZ)*D8INV
         SHLINT(INT,3) = SHLINT(INT,3) + (DIFABX*DIFABZ*YY
     &                                 +  DIFABY*DIFABY*XZ
     &                                 -  DIFABY*DIFABZ*XY
     &                                 -  DIFABX*DIFABY*YZ)*D8INV
         SHLINT(INT,4) = SHLINT(INT,4) + (D2*DIFABX*DIFABZ*XZ
     &                                 -  DIFABX*DIFABX*ZZ
     &                                 -  DIFABZ*DIFABZ*XX)*D8INV
         SHLINT(INT,5) = SHLINT(INT,5) + (DIFABX*DIFABX*YZ
     &                                 +  DIFABY*DIFABZ*XX
     &                                 -  DIFABX*DIFABY*XZ
     &                                 -  DIFABX*DIFABZ*XY)*D8INV
         SHLINT(INT,6) = SHLINT(INT,6) + (DIFABX*DIFABY*XY*D2
     &                                 -  DIFABY*DIFABY*XX
     &                                 -  DIFABX*DIFABX*YY)*D8INV
      END IF
C
C     Nuclear attraction contribution
C

         MAXT = LVALA + LVALB
         MAXU = MVALA + MVALB
         MAXV = NVALA + NVALB
         IADRAV = 1
         DO 200 IV = 0, MAXV + 2
            IADRAU = IADRAV
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            FV = ODC(NVALA,NVALB,IV,0,1,3)
            GV = ODC(NVALA,NVALB,IV,0,2,3)
            IUMAX = MAXU + 2
            IF (IV .GT. MAXV + 1) THEN
               IUMAX = MAXU
            ELSE IF (IV .GT. MAXV) THEN
               IUMAX = MAXU + 1
            END IF
            DO 300 IU = 0, IUMAX
               EU = ODC(MVALA,MVALB,IU,0,0,2)
               FU = ODC(MVALA,MVALB,IU,0,1,2)
               GU = ODC(MVALA,MVALB,IU,0,2,2)
               ITMAX = MAXT + 2
               IF (((IU .GT. MAXU + 1) .OR. (IV .GT. MAXV + 1))
     &               .OR. ((IU .EQ. MAXU + 1) .AND. (IV .EQ. MAXV + 1)))
     &               THEN
                  ITMAX = MAXT
               ELSE IF ((IU .GT. MAXU) .OR. (IV .GT. MAXV)) THEN
                  ITMAX = MAXT + 1
               END IF
               DO 400 IT = 0, ITMAX
                  ET = ODC(LVALA,LVALB,IT,0,0,1)
                  FT = ODC(LVALA,LVALB,IT,0,1,1)
                  GT = ODC(LVALA,LVALB,IT,0,2,1)
                  FXX = GT*EU*EV
                  FXY = FT*FU*EV
                  FXZ = FT*EU*FV
                  FYY = ET*GU*EV
                  FYZ = ET*FU*FV
                  FZZ = ET*EU*GV
                  ATUV = -D4INV*AHGTF(IADRAU + IT)
                  SHLINT(INT,1)=SHLINT(INT,1)+ATUV*(DIFABY*DIFABZ*FYZ*D2
     &                                            - DIFABY*DIFABY*FZZ
     &                                            - DIFABZ*DIFABZ*FYY)
                  SHLINT(INT,2)=SHLINT(INT,2)+ATUV*(DIFABZ*DIFABZ*FXY
     &                                            + DIFABY*DIFABX*FZZ
     &                                            - DIFABZ*DIFABX*FYZ
     &                                            - DIFABY*DIFABZ*FXZ)
                  SHLINT(INT,3)=SHLINT(INT,3)+ATUV*(DIFABX*DIFABZ*FYY
     &                                            + DIFABY*DIFABY*FXZ
     &                                            - DIFABY*DIFABZ*FXY
     &                                            - DIFABX*DIFABY*FYZ)
                  SHLINT(INT,4)=SHLINT(INT,4)+ATUV*(DIFABX*DIFABZ*FXZ*D2
     &                                            - DIFABX*DIFABX*FZZ
     &                                            - DIFABZ*DIFABZ*FXX)
                  SHLINT(INT,5)=SHLINT(INT,5)+ATUV*(DIFABX*DIFABX*FYZ
     &                                            + DIFABY*DIFABZ*FXX
     &                                            - DIFABX*DIFABY*FXZ
     &                                            - DIFABX*DIFABZ*FXY)
                  SHLINT(INT,6)=SHLINT(INT,6)+ATUV*(DIFABX*DIFABY*FXY*D2
     &                                            - DIFABY*DIFABY*FXX
     &                                            - DIFABX*DIFABX*FYY)
 400           CONTINUE
               IADRAU = IADRAU + ISTEPU
 300        CONTINUE
            IADRAV = IADRAV + ISTEPV
 200     CONTINUE
 100  CONTINUE
      RETURN
      END
C  /* Deck delgin */
      SUBROUTINE DELGIN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,
     &                  SHLINT,NOPTYP)
C
C     Diamagnetic one-electron spin-orbit integrals
C
C     K. Ruud, June 1997
C
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"
      PARAMETER (D2INV = 0.50D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,NOPTYP), AHGTF(*)
#include "onecom.h"
#include "lmns.h"
      INT = 1
      DO 100 ICOMPA = 1, KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         MAXT = LVALA + LVALB + 1
         MAXU = MVALA + MVALB + 1
         MAXV = NVALA + NVALB + 1
         IADRAV = 1
         DO 200 IV = 0, MAXV
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            FV = ODC(NVALA,NVALB,IV,0,1,3)
            IADRAU = IADRAV
            DO 300 IU = 0, MAXU
               EU = ODC(MVALA,MVALB,IU,0,0,2)
               FU = ODC(MVALA,MVALB,IU,0,1,2)
               DO 400 IT = 0, MAXT
                  ET = ODC(LVALA,LVALB,IT,0,0,1)
                  FT = ODC(LVALA,LVALB,IT,0,1,1)
                  FEE = FT*EU*EV
                  EFE = ET*FU*EV
                  EEF = ET*EU*FV
                  AH0T = AHGTF(IADRAU + IT + 1)
                  AH0U = AHGTF(IADRAU + IT + ISTEPU)
                  AH0V = AHGTF(IADRAU + IT + ISTEPV)
C
                  SHLINT(INT,1) = SHLINT(INT,1) - D2INV*
     &                              (EFE*AH0U + EEF*AH0V)
                  SHLINT(INT,2) = SHLINT(INT,2) + D2INV*FEE*AH0U
                  SHLINT(INT,3) = SHLINT(INT,3) + D2INV*FEE*AH0V
                  SHLINT(INT,4) = SHLINT(INT,4) + D2INV*EFE*AH0T
                  SHLINT(INT,5) = SHLINT(INT,5) - D2INV*
     *                              (FEE*AH0T + EEF*AH0V)
                  SHLINT(INT,6) = SHLINT(INT,6) + D2INV*EFE*AH0V
                  SHLINT(INT,7) = SHLINT(INT,7) + D2INV*EEF*AH0T
                  SHLINT(INT,8) = SHLINT(INT,8) + D2INV*EEF*AH0U
                  SHLINT(INT,9) = SHLINT(INT,9) - D2INV*
     &                              (FEE*AH0T + EFE*AH0U)
 400           CONTINUE
               IADRAU = IADRAU + ISTEPU
 300        CONTINUE
            IADRAV = IADRAV + ISTEPV
 200     CONTINUE
         INT = INT + 1
 100  CONTINUE
      RETURN
      END
C  /* Deck nsnlin */
      SUBROUTINE NSNLIN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,
     &                  SHLINT,NOPTYP,NATOMC)
C
C     Nuclear shielding tensor integrals without London orbital contribution
C
C     K. Ruud, March 1992
C
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"
      PARAMETER (D2INV = 0.50D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,NOPTYP), AHGTF(*)
#include "onecom.h"
#include "lmns.h"
      INT = 1
      DO 100 ICOMPA = 1, KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         MAXT = LVALA + LVALB + 1
         MAXU = MVALA + MVALB + 1
         MAXV = NVALA + NVALB + 1
         IADRAV = 1
         DO 200 IV = 0, MAXV
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            FV = ODC(NVALA,NVALB,IV,0,1,3)
            IADRAU = IADRAV
            DO 300 IU = 0, MAXU
               EU = ODC(MVALA,MVALB,IU,0,0,2)
               FU = ODC(MVALA,MVALB,IU,0,1,2)
               DO 400 IT = 0, MAXT
                  ET = ODC(LVALA,LVALB,IT,0,0,1)
                  FT = ODC(LVALA,LVALB,IT,0,1,1)
                  FEE = FT*EU*EV
                  EFE = ET*FU*EV
                  EEF = ET*EU*FV
                  IOFF = 0
                  DO 500 IATOM = 1, NATOMC
                     IXX = 9*(IATOM - 1) + 1
                     IYX = 9*(IATOM - 1) + 2
                     IZX = 9*(IATOM - 1) + 3
                     IXY = 9*(IATOM - 1) + 4
                     IYY = 9*(IATOM - 1) + 5
                     IZY = 9*(IATOM - 1) + 6
                     IXZ = 9*(IATOM - 1) + 7
                     IYZ = 9*(IATOM - 1) + 8
                     IZZ = 9*(IATOM - 1) + 9
                     AH0T = AHGTF(IOFF + IADRAU + IT + 1)
                     AH0U = AHGTF(IOFF + IADRAU + IT + ISTEPU)
                     AH0V = AHGTF(IOFF + IADRAU + IT + ISTEPV)
C
                     SHLINT(INT,IXX) = SHLINT(INT,IXX) - D2INV*
     &                                 (EFE*AH0U + EEF*AH0V)
                     SHLINT(INT,IYX) = SHLINT(INT,IYX) + D2INV*FEE*AH0U
                     SHLINT(INT,IZX) = SHLINT(INT,IZX) + D2INV*FEE*AH0V
                     SHLINT(INT,IXY) = SHLINT(INT,IXY) + D2INV*EFE*AH0T
                     SHLINT(INT,IYY) = SHLINT(INT,IYY) - D2INV*
     *                                 (FEE*AH0T + EEF*AH0V)
                     SHLINT(INT,IZY) = SHLINT(INT,IZY) + D2INV*EFE*AH0V
                     SHLINT(INT,IXZ) = SHLINT(INT,IXZ) + D2INV*EEF*AH0T
                     SHLINT(INT,IYZ) = SHLINT(INT,IYZ) + D2INV*EEF*AH0U
                     SHLINT(INT,IZZ) = SHLINT(INT,IZZ) - D2INV*
     &                                 (FEE*AH0T + EFE*AH0U)
                     IOFF = IOFF + NAHGTF
 500              CONTINUE
 400           CONTINUE
               IADRAU = IADRAU + ISTEPU
 300        CONTINUE
            IADRAV = IADRAV + ISTEPV
 200     CONTINUE
         INT = INT + 1
 100  CONTINUE
      RETURN
      END
      SUBROUTINE NSKEIN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,
     &                  SHLINT,NOPTYP,NATOMC)
C
C     Kinetic energy correction to diamagnetic nuclear shielding tensor
C     integrals with common gauge origin
C
C     K. Ruud, Sep.2001
C
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"
      PARAMETER (DP75 = 0.75D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,NOPTYP), AHGTF(*)
#include "onecom.h"
#include "lmns.h"
      INT = 1
      DO 100 ICOMPA = 1, KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         MAXT = LVALA + LVALB + 3
         MAXU = MVALA + MVALB + 3
         MAXV = NVALA + NVALB + 3
         IADRAV = 1
         DO 200 IV = 0, MAXV
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            FV1 = ODC(NVALA,NVALB,IV,0,1,3)
            FV2 = ODC(NVALA,NVALB,IV,2,0,3)
            FV3 = ODC(NVALA,NVALB,IV,2,1,3)
            IADRAU = IADRAV
            DO 300 IU = 0,MAXU
               EU  = ODC(MVALA,MVALB,IU,0,0,2)
               FU1 = ODC(MVALA,MVALB,IU,0,1,2)
               FU2 = ODC(MVALA,MVALB,IU,2,0,2)
               FU3 = ODC(MVALA,MVALB,IU,2,1,2)
C
               EE   = EU*EV
               F1E  = FU1*EV
               F2E  = FU2*EV
               F3E  = FU3*EV
               EF1  = EU*FV1
               F2F1 = FU2*FV1
               F1F2 = FU1*FV2
               EF2  = EU*FV2
               EF3  = EU*FV3
               DO 400 IT = 0,MAXT
                  ET  = ODC(LVALA,LVALB,IT,0,0,1)
                  FT1 = ODC(LVALA,LVALB,IT,0,1,1)
                  FT2 = ODC(LVALA,LVALB,IT,2,0,1)
                  FT3 = ODC(LVALA,LVALB,IT,2,1,1)
C
                  F1F2E = FT1*F2E
                  F3EE  = FT3*EE
                  EF3E  = ET*F3E
                  F2F1E = FT2*F1E
                  F2EF1 = FT2*EF1
                  EF2F1 = ET*F2F1
                  EF1F2 = ET*F1F2
                  F1EF2 = FT1*EF2
                  EEF3  = ET*EF3
                  IOFF = 0
                  DO 500 IATOM = 1, NATOMC
                     IXX = 9*(IATOM - 1) + 1
                     IYX = 9*(IATOM - 1) + 2
                     IZX = 9*(IATOM - 1) + 3
                     IXY = 9*(IATOM - 1) + 4
                     IYY = 9*(IATOM - 1) + 5
                     IZY = 9*(IATOM - 1) + 6
                     IXZ = 9*(IATOM - 1) + 7
                     IYZ = 9*(IATOM - 1) + 8
                     IZZ = 9*(IATOM - 1) + 9
                     AH0T = AHGTF(IOFF + IADRAU + IT + 1)
                     AH0U = AHGTF(IOFF + IADRAU + IT + ISTEPU)
                     AH0V = AHGTF(IOFF + IADRAU + IT + ISTEPV)
C
                     SHLINT(INT,IXX) = SHLINT(INT,IXX) - DP75*
     &                                 (AH0U*(F2F1E + EF1F2 + EF3E) +
     &                                  AH0V*(F2EF1 + EF2F1 + EEF3))
                     SHLINT(INT,IYX) = SHLINT(INT,IYX) + DP75*AH0U*
     &                                 (F3EE + F1F2E + F1EF2)
                     SHLINT(INT,IZX) = SHLINT(INT,IZX) + DP75*AH0V*
     &                                 (F3EE + F1F2E + F1EF2)
                     SHLINT(INT,IXY) = SHLINT(INT,IXY) + DP75*AH0T*
     &                                 (F2F1E + EF3E + EF1F2)
                     SHLINT(INT,IYY) = SHLINT(INT,IYY) - DP75*
     &                                 (AH0T*(F3EE + F1F2E + F1EF2) +
     &                                  AH0V*(F2EF1 + EF2F1 + EEF3))
                     SHLINT(INT,IZY) = SHLINT(INT,IZY) + DP75*AH0V*
     &                                 (F2F1E + EF3E + EF1F2)
                     SHLINT(INT,IXZ) = SHLINT(INT,IXZ) + DP75*AH0T*
     &                                 (F2EF1 + EF2F1 + EEF3)
                     SHLINT(INT,IYZ) = SHLINT(INT,IYZ) + DP75*AH0U*
     &                                 (F2EF1 + EF2F1 + EEF3)
                     SHLINT(INT,IZZ) = SHLINT(INT,IZZ) - DP75*
     &                                 (AH0T*(F3EE + F1F2E + F1EF2) +
     &                                  AH0U*(F2F1E + EF1F2 + EF3E))
                     IOFF = IOFF + NAHGTF
 500              CONTINUE
 400           CONTINUE
               IADRAU = IADRAU + ISTEPU
 300        CONTINUE
            IADRAV = IADRAV + ISTEPV
 200     CONTINUE
         INT = INT + 1
 100  CONTINUE
      RETURN
      END
C  /* Deck sofint */
      SUBROUTINE SOFINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,
     &                  SHLINT,NOPTYP,NATOMC)
C
C     Magnetic-field corrections to spin-orbit integrals
C
C     K. Ruud, January 1998
C
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"
      PARAMETER (D2INV = 0.50D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,NOPTYP), AHGTF(*)
#include "onecom.h"
#include "lmns.h"
      INT = 1
      DO 100 ICOMPA = 1, KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         MAXT = LVALA + LVALB + 1
         MAXU = MVALA + MVALB + 1
         MAXV = NVALA + NVALB + 1
         IADRAV = 1
         DO 200 IV = 0, MAXV
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            FV = ODC(NVALA,NVALB,IV,0,1,3)
            IADRAU = IADRAV
            DO 300 IU = 0, MAXU
               EU = ODC(MVALA,MVALB,IU,0,0,2)
               FU = ODC(MVALA,MVALB,IU,0,1,2)
               DO 400 IT = 0, MAXT
                  ET = ODC(LVALA,LVALB,IT,0,0,1)
                  FT = ODC(LVALA,LVALB,IT,0,1,1)
                  FEE = FT*EU*EV
                  EFE = ET*FU*EV
                  EEF = ET*EU*FV
                  AH0T = AHGTF(IADRAU + IT + 1)
                  AH0U = AHGTF(IADRAU + IT + ISTEPU)
                  AH0V = AHGTF(IADRAU + IT + ISTEPV)
C
                  SHLINT(INT,1) = SHLINT(INT,1) - D2INV*
     &                              (EFE*AH0U + EEF*AH0V)
                  SHLINT(INT,2) = SHLINT(INT,2) + D2INV*FEE*AH0U
                  SHLINT(INT,3) = SHLINT(INT,3) + D2INV*FEE*AH0V
                  SHLINT(INT,4) = SHLINT(INT,4) + D2INV*EFE*AH0T
                  SHLINT(INT,5) = SHLINT(INT,5) - D2INV*
     *                              (FEE*AH0T + EEF*AH0V)
                  SHLINT(INT,6) = SHLINT(INT,6) + D2INV*EFE*AH0V
                  SHLINT(INT,7) = SHLINT(INT,7) + D2INV*EEF*AH0T
                  SHLINT(INT,8) = SHLINT(INT,8) + D2INV*EEF*AH0U
                  SHLINT(INT,9) = SHLINT(INT,9) - D2INV*
     &                              (FEE*AH0T + EFE*AH0U)
 400           CONTINUE
               IADRAU = IADRAU + ISTEPU
 300        CONTINUE
            IADRAV = IADRAV + ISTEPV
 200     CONTINUE
         INT = INT + 1
 100  CONTINUE
      RETURN
      END
C  /* Deck nsloin */
      SUBROUTINE NSLOIN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,
     &                  SHLINT,NOPTYP,NATOMC,DIFABX,DIFABY,DIFABZ)
C
C     London orbital contribution to nuclear shielding integrals
C
C     K. Ruud, March 1992
C
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"
      PARAMETER (D2INV = 0.50D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,NOPTYP), AHGTF(*)
#include "onecom.h"
#include "lmns.h"
      INT = 1
      DO 100 ICOMPA = 1, KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         MAXT = LVALA + LVALB + 2
         MAXU = MVALA + MVALB + 2
         MAXV = NVALA + NVALB + 2
         IADRAV = 1
         DO 200 IV = 0, MAXV
            DV = ODC(NVALA,NVALB,IV,1,1,3)
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            FV = ODC(NVALA,NVALB,IV,1,0,3)
            GV = ODC(NVALA,NVALB,IV,0,1,3)
            IADRAU = IADRAV
            DO 300 IU = 0, MAXU
               DU = ODC(MVALA,MVALB,IU,1,1,2)
               EU = ODC(MVALA,MVALB,IU,0,0,2)
               FU = ODC(MVALA,MVALB,IU,1,0,2)
               GU = ODC(MVALA,MVALB,IU,0,1,2)
               DO 400 IT = 0, MAXT
                  DT = ODC(LVALA,LVALB,IT,1,1,1)
                  ET = ODC(LVALA,LVALB,IT,0,0,1)
                  FT = ODC(LVALA,LVALB,IT,1,0,1)
                  GT = ODC(LVALA,LVALB,IT,0,1,1)
                  IOFF = 0
                  DO 500 IATOM = 1, NATOMC
                     IXX = 9*(IATOM - 1) + 1
                     IYX = 9*(IATOM - 1) + 2
                     IZX = 9*(IATOM - 1) + 3
                     IXY = 9*(IATOM - 1) + 4
                     IYY = 9*(IATOM - 1) + 5
                     IZY = 9*(IATOM - 1) + 6
                     IXZ = 9*(IATOM - 1) + 7
                     IYZ = 9*(IATOM - 1) + 8
                     IZZ = 9*(IATOM - 1) + 9
                     AH0T = AHGTF(IOFF + IADRAU + IT + 1)
                     AH0U = AHGTF(IOFF + IADRAU + IT + ISTEPU)
                     AH0V = AHGTF(IOFF + IADRAU + IT + ISTEPV)
                     SHLINT(INT,IXX) = SHLINT(INT,IXX) - (DIFABY*ET*
     &                                 (FU*GV*AH0V - EU*DV*AH0U)
     &                               + (GU*FV*AH0U - DU*EV*AH0V)*
     &                                 DIFABZ*ET)*D2INV
                     SHLINT(INT,IYX) = SHLINT(INT,IYX) - (DIFABZ*GT*
     &                                 (FU*EV*AH0V - EU*FV*AH0U)
     &                               + (EU*DV*AH0U - FU*GV*AH0V)*
     &                                 DIFABX*ET)*D2INV
                     SHLINT(INT,IZX) = SHLINT(INT,IZX) - (DIFABX*ET*
     &                                 (DU*EV*AH0V - GU*FV*AH0U)
     &                               + (EU*FV*AH0U - FU*EV*AH0V)*
     &                                 DIFABY*GT)*D2INV
                     SHLINT(INT,IXY) = SHLINT(INT,IXY) + (DIFABY*EU*
     &                                 (FT*GV*AH0V - ET*DV*AH0T)
     &                               + (ET*FV*AH0T - FT*EV*AH0V)*
     &                                 DIFABZ*GU)*D2INV
                     SHLINT(INT,IYY) = SHLINT(INT,IYY) + (DIFABZ*EU*
     &                                 (DT*EV*AH0V - GT*FV*AH0T)
     &                               + (ET*DV*AH0T - FT*GV*AH0V)*
     &                                 DIFABX*EU)*D2INV
                     SHLINT(INT,IZY) = SHLINT(INT,IZY) + (DIFABX*GU*
     &                                 (FT*EV*AH0V - ET*FV*AH0T)
     &                               + (GT*FV*AH0T - DT*EV*AH0V)*
     &                                 DIFABY*EU)*D2INV
                     SHLINT(INT,IXZ) = SHLINT(INT,IXZ) - (DIFABY*GV*
     &                                 (FT*EU*AH0U - ET*FU*AH0T)
     &                               + (ET*DU*AH0T - FT*GU*AH0U)*
     &                                 DIFABZ*EV)*D2INV
                     SHLINT(INT,IYZ) = SHLINT(INT,IYZ) - (DIFABZ*EV*
     &                                 (DT*EU*AH0U - GT*FU*AH0T)
     &                               + (ET*FU*AH0T - FT*EU*AH0U)*
     &                                 DIFABX*GV)*D2INV
                     SHLINT(INT,IZZ) = SHLINT(INT,IZZ) - (DIFABX*EV*
     &                                 (FT*GU*AH0U - ET*DU*AH0T)
     &                               + (GT*FU*AH0T - DT*EU*AH0U)*
     &                                 DIFABY*EV)*D2INV
                     IOFF = IOFF + NAHGTF
 500              CONTINUE
 400           CONTINUE
               IADRAU = IADRAU + ISTEPU
 300        CONTINUE
            IADRAV = IADRAV + ISTEPV
 200     CONTINUE
         INT = INT + 1
 100  CONTINUE
      RETURN
      END
C  /* Deck npeint */
      SUBROUTINE NPEINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT,
     &                  NOPTYP,NATOMC)
C
C     Potential energy from each individual nucleus IATOM in SHLINT(*,IATOM).
C
C     The total potential energy calculated in POTINT is the sum :
C         POTINT = sum(iatom=1,natomc) NPEINT(iatom)
C
C     K.Ruud, July 1992
C
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &              SHLINT(KCKTAB,NOPTYP), AHGTF(*)
#include "onecom.h"
#include "lmns.h"
C
      INT = 1
      DO 100 ICOMPA = 1, KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         MAXT = LVALA + LVALB
         MAXU = MVALA + MVALB
         MAXV = NVALA + NVALB
         IADRAV = 1
         DO 200 IV = 0, MAXV
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            IADRAU = IADRAV
            DO 300 IU = 0, MAXU
               EU = ODC(MVALA,MVALB,IU,0,0,2)
               DO 400 IT = 0, MAXT
                  ET = ODC(LVALA,LVALB,IT,0,0,1)
                  EEE = ET*EU*EV
                  IOFF = 0
                  DO 500 IATOM = 1, NATOMC
                     SHLINT(INT,IATOM) = SHLINT(INT,IATOM)
     &                                 + EEE*AHGTF(IOFF + IADRAU + IT)
                     IOFF = IOFF + NAHGTF
 500              CONTINUE
 400           CONTINUE
               IADRAU = IADRAU + ISTEPU
 300        CONTINUE
            IADRAV = IADRAV + ISTEPV
 200     CONTINUE
         INT = INT + 1
 100  CONTINUE
      RETURN
      END
C  /* Deck potint */
      SUBROUTINE POTINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT)
C
C     One-electron potential energy from all the nuclei
C     (use NPEINT for potential energy from each individual energy)
C
C     K.Ruud, September 2000
C
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &              SHLINT(KCKTAB), AHGTF(*)
#include "onecom.h"
#include "lmns.h"
C
      INT = 1
      DO 100 ICOMPA = 1, KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         MAXT = LVALA + LVALB
         MAXU = MVALA + MVALB
         MAXV = NVALA + NVALB
         IADRAV = 1
         DO 200 IV = 0, MAXV
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            IADRAU = IADRAV
            DO 300 IU = 0, MAXU
               EU = ODC(MVALA,MVALB,IU,0,0,2)
               DO 400 IT = 0, MAXT
                  ET = ODC(LVALA,LVALB,IT,0,0,1)
                  EEE = ET*EU*EV
                  SHLINT(INT) = SHLINT(INT)
     &                        + EEE*AHGTF(IADRAU + IT)
 400           CONTINUE
               IADRAU = IADRAU + ISTEPU
 300        CONTINUE
            IADRAV = IADRAV + ISTEPV
 200     CONTINUE
         INT = INT + 1
 100  CONTINUE
      RETURN
      END
CC  /* Deck efgint */
      SUBROUTINE EFGINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT,
     &                  NOPTYP,NATOMC,DISTAB)
C
C     Electric field gradient integrals
C
C     K.Ruud, June 1992
C
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxmom.h"
#include "maxorb.h"
      PARAMETER (D1 = 1.0D0, D3 = 3.0D0, D3INV = D1/D3,
     &           THRESH = 1.0D-10)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &              SHLINT(KCKTAB,NOPTYP), AHGTF(*)
#include "nuclei.h"
#include "symmet.h"
#include "onecom.h"
#include "lmns.h"
#include "ibtfun.h"
C
      INT = 1
      DO 100 ICOMPA = 1, KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         MAXT = LVALA + LVALB
         MAXU = MVALA + MVALB
         MAXV = NVALA + NVALB
         IADRAV = 1
         DO 200 IV = 0, MAXV
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            IADRAU = IADRAV
            DO 300 IU = 0, MAXU
               EU = ODC(MVALA,MVALB,IU,0,0,2)
               DO 400 IT = 0, MAXT
                  ET = ODC(LVALA,LVALB,IT,0,0,1)
                  EEE = ET*EU*EV
                  IOFF = 0
                  INTTYP = 0
                  DO 500 IATOM = 1, NATOMC
                     IADR0  = IOFF  + IADRAU + IT
                     IADRTT = IADR0 + 2
                     IADRTU = IADR0 + 1 + ISTEPU
                     IADRTV = IADR0 + 1 + ISTEPV
                     IADRUU = IADR0 + 2*ISTEPU
                     IADRUV = IADR0 + ISTEPU + ISTEPV
                     IADRVV = IADR0 + 2*ISTEPV
C
C     Possible numerical instabilities, as pointed out by Luuk Visscher,
C     still remain. K.Ruud-Dec96
C
                        SHLINT(INT,INTTYP+1) = SHLINT(INT,INTTYP+1)
     &                                  + EEE*D3INV*(2*AHGTF(IADRTT)
     &                                  - AHGTF(IADRUU) - AHGTF(IADRVV))
                        SHLINT(INT,INTTYP+4) = SHLINT(INT,INTTYP+4)
     &                                  + EEE*D3INV*(2*AHGTF(IADRUU)
     &                                  - AHGTF(IADRTT) - AHGTF(IADRVV))
                        SHLINT(INT,INTTYP+6) = SHLINT(INT,INTTYP+6)
     &                                  + EEE*D3INV*(2*AHGTF(IADRVV)
     &                                  - AHGTF(IADRTT) - AHGTF(IADRUU))
                     SHLINT(INT,INTTYP+2) = SHLINT(INT,INTTYP+2)
     &                                  + EEE*AHGTF(IADRTU)
                     SHLINT(INT,INTTYP+3) = SHLINT(INT,INTTYP+3)
     &                                  + EEE*AHGTF(IADRTV)
                     SHLINT(INT,INTTYP+5) = SHLINT(INT,INTTYP+5)
     &                                  + EEE*AHGTF(IADRUV)
                     IOFF = IOFF + NAHGTF
                     INTTYP = INTTYP + 6
 500              CONTINUE
 400           CONTINUE
               IADRAU = IADRAU + ISTEPU
 300        CONTINUE
            IADRAV = IADRAV + ISTEPV
 200     CONTINUE
         INT = INT + 1
 100  CONTINUE
      RETURN
      END
C  /* Deck ef1int */
      SUBROUTINE EF1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT,
     &                  NOPTYP,NATOMC)
C
C     Nuclear electric field integrals
C
C     K.Ruud, March 1992
C
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,NOPTYP), AHGTF(*)
#include "onecom.h"
#include "lmns.h"
      INT = 1
      DO 100 ICOMPA = 1, KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         MAXT = LVALA + LVALB
         MAXU = MVALA + MVALB
         MAXV = NVALA + NVALB
         IADRAV = 1
         DO 200 IV = 0, MAXV
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            IADRAU = IADRAV
            DO 300 IU = 0, MAXU
               EU = ODC(MVALA,MVALB,IU,0,0,2)
               DO 400 IT = 0, MAXT
                  ET = ODC(LVALA,LVALB,IT,0,0,1)
                  EEE = ET*EU*EV
                  IOFF = 0
                  DO 500 IATOM = 1, NATOMC
                     IX = 3*(IATOM - 1) + 1
                     IY = 3*(IATOM - 1) + 2
                     IZ = 3*(IATOM - 1) + 3
                     AH0T = AHGTF(IOFF + IADRAU + IT + 1)
                     AH0U = AHGTF(IOFF + IADRAU + IT + ISTEPU)
                     AH0V = AHGTF(IOFF + IADRAU + IT + ISTEPV)
                     SHLINT(INT,IX) = SHLINT(INT,IX) - EEE*AH0T
                     SHLINT(INT,IY) = SHLINT(INT,IY) - EEE*AH0U
                     SHLINT(INT,IZ) = SHLINT(INT,IZ) - EEE*AH0V
                     IOFF = IOFF + NAHGTF
 500              CONTINUE
 400           CONTINUE
               IADRAU = IADRAU + ISTEPU
 300        CONTINUE
            IADRAV = IADRAV + ISTEPV
 200     CONTINUE
         INT = INT + 1
 100  CONTINUE
      RETURN
      END
C  /* Deck efbint */
      SUBROUTINE EFBINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT,
     &                  DIFABX,DIFABY,DIFABZ)
C
C     London orbital correction to the electric field calculated at a
C     point (e.g. in a MM system) in space.
C
C     K.Ruud, Copenhagen January 2006
C
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"
      PARAMETER (DP5 = 0.50D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,9), AHGTF(*)
#include "onecom.h"
#include "lmns.h"
      INT = 1
      DO 100 ICOMPA = 1, KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         MAXT = LVALA + LVALB
         MAXU = MVALA + MVALB
         MAXV = NVALA + NVALB
         IADRAV = 1
         DO 200 IV = 0, MAXV + 1
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            FV = ODC(NVALA,NVALB,IV,0,1,3)
            IADRAU = IADRAV
            IUMAX = MAXU + 1
            IF (IV .GT. MAXV) IUMAX = MAXU
            DO 300 IU = 0, IUMAX
               EU = ODC(MVALA,MVALB,IU,0,0,2)
               FU = ODC(MVALA,MVALB,IU,0,1,2)
               ITMAX = MAXT + 1
               IF ((IU .GT. MAXU) .OR. (IV .GT. MAXV)) ITMAX = MAXT
               DO 400 IT = 0, ITMAX
                  ET = ODC(LVALA,LVALB,IT,0,0,1)
                  FT = ODC(LVALA,LVALB,IT,0,1,1)
                  FX = FT*EU*EV
                  FY = ET*FU*EV
                  FZ = ET*EU*FV
                  AH0T = DP5*AHGTF(IADRAU + IT + 1)
                  AH0U = DP5*AHGTF(IADRAU + IT + ISTEPU)
                  AH0V = DP5*AHGTF(IADRAU + IT + ISTEPV)
                  SHLINT(INT,1) = SHLINT(INT,1)
     &                           + AH0T*(DIFABY*FZ-DIFABZ*FY)
                  SHLINT(INT,2) = SHLINT(INT,2)
     &                           + AH0T*(DIFABZ*FX-DIFABX*FZ)
                  SHLINT(INT,3) = SHLINT(INT,3)
     &                           + AH0T*(DIFABX*FY-DIFABY*FX)
                  SHLINT(INT,4) = SHLINT(INT,4)
     &                           + AH0U*(DIFABY*FZ-DIFABZ*FY)
                  SHLINT(INT,5) = SHLINT(INT,5)
     &                           + AH0U*(DIFABZ*FX-DIFABX*FZ)
                  SHLINT(INT,6) = SHLINT(INT,6)
     &                           + AH0U*(DIFABX*FY-DIFABY*FX)
                  SHLINT(INT,7) = SHLINT(INT,7)
     &                           + AH0V*(DIFABY*FZ-DIFABZ*FY)
                  SHLINT(INT,8) = SHLINT(INT,8)
     &                           + AH0V*(DIFABZ*FX-DIFABX*FZ)
                  SHLINT(INT,9) = SHLINT(INT,9)
     &                           + AH0V*(DIFABX*FY-DIFABY*FX)
 400           CONTINUE
               IADRAU = IADRAU + ISTEPU
 300        CONTINUE
            IADRAV = IADRAV + ISTEPV
 200     CONTINUE
         INT = INT + 1
 100  CONTINUE
      RETURN
      END
C  /* Deck efb2in */
      SUBROUTINE EFB2IN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT,
     &                  DIFABX,DIFABY,DIFABZ)
C
C     2.nd-order London orbital correction to the electric field
C     calculated at a point (e.g. in a MM system) in space.
C
C     K.Ruud, Copenhagen January 2006
C
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"
      PARAMETER (DP25 = 0.25D0, D2 = 2.0D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,18), AHGTF(*)
#include "onecom.h"
#include "lmns.h"
      INT = 1
      DO 100 ICOMPA = 1, KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         MAXT = LVALA + LVALB
         MAXU = MVALA + MVALB
         MAXV = NVALA + NVALB
         IADRAV = 1
         DO 200 IV = 0, MAXV + 2
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            FV = ODC(NVALA,NVALB,IV,0,1,3)
            GV = ODC(NVALA,NVALB,IV,0,2,3)
            IADRAU = IADRAV
            IUMAX = MAXU + 2
            IF (IV .GT. MAXV + 1) THEN
               IUMAX = MAXU
            ELSE IF (IV .GT. MAXV) THEN
               IUMAX = MAXU + 1
            END IF
            DO 300 IU = 0, IUMAX
               EU = ODC(MVALA,MVALB,IU,0,0,2)
               FU = ODC(MVALA,MVALB,IU,0,1,2)
               GU = ODC(MVALA,MVALB,IU,0,2,2)
               ITMAX = MAXT + 2
               IF (((IU .GT. MAXU + 1) .OR. (IV .GT. MAXV + 1))
     &               .OR. ((IU .EQ. MAXU + 1) .AND. (IV .EQ. MAXV + 1)))
     &               THEN
                  ITMAX = MAXT
               ELSE IF ((IU .GT. MAXU) .OR. (IV .GT. MAXV)) THEN
                  ITMAX = MAXT + 1
               END IF
               DO 400 IT = 0, ITMAX
                  ET = ODC(LVALA,LVALB,IT,0,0,1)
                  FT = ODC(LVALA,LVALB,IT,0,1,1)
                  GT = ODC(LVALA,LVALB,IT,0,2,1)
                  FXX = GT*EU*EV
                  FXY = FT*FU*EV
                  FXZ = FT*EU*FV
                  FYY = ET*GU*EV
                  FYZ = ET*FU*FV
                  FZZ = ET*EU*GV
                  AH0T = DP25*AHGTF(IADRAU + IT + 1)
                  AH0U = DP25*AHGTF(IADRAU + IT + ISTEPU)
                  AH0V = DP25*AHGTF(IADRAU + IT + ISTEPV)
                SHLINT(INT,1)=SHLINT(INT,1)+AH0T*(DIFABY*DIFABZ*FYZ*D2
     &                                            - DIFABY*DIFABY*FZZ
     &                                            - DIFABZ*DIFABZ*FYY)
                SHLINT(INT,2)=SHLINT(INT,2)+AH0T*(DIFABZ*DIFABZ*FXY
     &                                            + DIFABY*DIFABX*FZZ
     &                                            - DIFABZ*DIFABX*FYZ
     &                                            - DIFABY*DIFABZ*FXZ)
                SHLINT(INT,3)=SHLINT(INT,3)+AH0T*(DIFABX*DIFABZ*FYY
     &                                            + DIFABY*DIFABY*FXZ
     &                                            - DIFABY*DIFABZ*FXY
     &                                            - DIFABX*DIFABY*FYZ)
                SHLINT(INT,4)=SHLINT(INT,4)+AH0T*(DIFABX*DIFABZ*FXZ*D2
     &                                            - DIFABX*DIFABX*FZZ
     &                                            - DIFABZ*DIFABZ*FXX)
                SHLINT(INT,5)=SHLINT(INT,5)+AH0T*(DIFABX*DIFABX*FYZ
     &                                            + DIFABY*DIFABZ*FXX
     &                                            - DIFABX*DIFABY*FXZ
     &                                            - DIFABX*DIFABZ*FXY)
                SHLINT(INT,6)=SHLINT(INT,6)+AH0T*(DIFABX*DIFABY*FXY*D2
     &                                            - DIFABY*DIFABY*FXX
     &                                            - DIFABX*DIFABX*FYY)
                SHLINT(INT,7)=SHLINT(INT,7)+AH0U*(DIFABY*DIFABZ*FYZ*D2
     &                                            - DIFABY*DIFABY*FZZ
     &                                            - DIFABZ*DIFABZ*FYY)
                SHLINT(INT,8)=SHLINT(INT,8)+AH0U*(DIFABZ*DIFABZ*FXY
     &                                            + DIFABY*DIFABX*FZZ
     &                                            - DIFABZ*DIFABX*FYZ
     &                                            - DIFABY*DIFABZ*FXZ)
                SHLINT(INT,9)=SHLINT(INT,9)+AH0U*(DIFABX*DIFABZ*FYY
     &                                            + DIFABY*DIFABY*FXZ
     &                                            - DIFABY*DIFABZ*FXY
     &                                            - DIFABX*DIFABY*FYZ)
                SHLINT(INT,10)=SHLINT(INT,10)+AH0U*(DIFABX*DIFABZ*FXZ*D2
     &                                            - DIFABX*DIFABX*FZZ
     &                                            - DIFABZ*DIFABZ*FXX)
                SHLINT(INT,11)=SHLINT(INT,11)+AH0U*(DIFABX*DIFABX*FYZ
     &                                            + DIFABY*DIFABZ*FXX
     &                                            - DIFABX*DIFABY*FXZ
     &                                            - DIFABX*DIFABZ*FXY)
                SHLINT(INT,12)=SHLINT(INT,12)+AH0U*(DIFABX*DIFABY*FXY*D2
     &                                            - DIFABY*DIFABY*FXX
     &                                            - DIFABX*DIFABX*FYY)
                SHLINT(INT,13)=SHLINT(INT,13)+AH0V*(DIFABY*DIFABZ*FYZ*D2
     &                                            - DIFABY*DIFABY*FZZ
     &                                            - DIFABZ*DIFABZ*FYY)
                SHLINT(INT,14)=SHLINT(INT,14)+AH0V*(DIFABZ*DIFABZ*FXY
     &                                            + DIFABY*DIFABX*FZZ
     &                                            - DIFABZ*DIFABX*FYZ
     &                                            - DIFABY*DIFABZ*FXZ)
                SHLINT(INT,15)=SHLINT(INT,15)+AH0V*(DIFABX*DIFABZ*FYY
     &                                            + DIFABY*DIFABY*FXZ
     &                                            - DIFABY*DIFABZ*FXY
     &                                            - DIFABX*DIFABY*FYZ)
                SHLINT(INT,16)=SHLINT(INT,16)+AH0V*(DIFABX*DIFABZ*FXZ*D2
     &                                            - DIFABX*DIFABX*FZZ
     &                                            - DIFABZ*DIFABZ*FXX)
                SHLINT(INT,17)=SHLINT(INT,17)+AH0V*(DIFABX*DIFABX*FYZ
     &                                            + DIFABY*DIFABZ*FXX
     &                                            - DIFABX*DIFABY*FXZ
     &                                            - DIFABX*DIFABZ*FXY)
                SHLINT(INT,18)=SHLINT(INT,18)+AH0V*(DIFABX*DIFABY*FXY*D2
     &                                            - DIFABY*DIFABY*FXX
     &                                            - DIFABX*DIFABX*FYY)
 400           CONTINUE
               IADRAU = IADRAU + ISTEPU
 300        CONTINUE
            IADRAV = IADRAV + ISTEPV
 200     CONTINUE
         INT = INT + 1
 100  CONTINUE
      RETURN
      END
C  /* Deck qdpint */
      SUBROUTINE QDPINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,EXPPI,INTTYP)
C
C     tuh 1989
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0, D1P5 = 1.5D0, D4 = 4.0D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,6)
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         DX1 = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         DY1 = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         DZ1 = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
         DX2 = SHGTF*ODC(LVALA,LVALB,0,0,2,1)
         DY2 = SHGTF*ODC(MVALA,MVALB,0,0,2,2)
         DZ2 = SHGTF*ODC(NVALA,NVALB,0,0,2,3)
C
         INT = INT + 1
         IF (INTTYP .EQ. 6) THEN
C           'SECMOM '
            SHLINT(INT,1) = SHLINT(INT,1) - DX2*SY0*SZ0
            SHLINT(INT,2) = SHLINT(INT,2) - DX1*DY1*SZ0
            SHLINT(INT,3) = SHLINT(INT,3) - DX1*SY0*DZ1
            SHLINT(INT,4) = SHLINT(INT,4) - SX0*DY2*SZ0
            SHLINT(INT,5) = SHLINT(INT,5) - SX0*DY1*DZ1
            SHLINT(INT,6) = SHLINT(INT,6) - SX0*SY0*DZ2
         ELSE IF (INTTYP .EQ. 7) THEN
C           'THETA  '
            X2 = DX2*SY0*SZ0
            Y2 = SX0*DY2*SZ0
            Z2 = SX0*SY0*DZ2
            R2 = (X2 + Y2 + Z2)/D2
            SHLINT(INT,1) = SHLINT(INT,1) - D1P5*X2 + R2
            SHLINT(INT,2) = SHLINT(INT,2) - D1P5*DX1*DY1*SZ0
            SHLINT(INT,3) = SHLINT(INT,3) - D1P5*DX1*SY0*DZ1
            SHLINT(INT,4) = SHLINT(INT,4) - D1P5*Y2 + R2
            SHLINT(INT,5) = SHLINT(INT,5) - D1P5*SX0*DY1*DZ1
            SHLINT(INT,6) = SHLINT(INT,6) - D1P5*Z2 + R2
         ELSE
C           'QUADRUP' and various diamagn. susceptibility integrals
            X2 = DX2*SY0*SZ0
            Y2 = SX0*DY2*SZ0
            Z2 = SX0*SY0*DZ2
            R2 = (X2 + Y2 + Z2)/D4
            SHLINT(INT,1) = SHLINT(INT,1) - R2 + X2/D4
            SHLINT(INT,2) = SHLINT(INT,2) + DX1*DY1*SZ0/D4
            SHLINT(INT,3) = SHLINT(INT,3) + DX1*SY0*DZ1/D4
            SHLINT(INT,4) = SHLINT(INT,4) - R2 + Y2/D4
            SHLINT(INT,5) = SHLINT(INT,5) + SX0*DY1*DZ1/D4
            SHLINT(INT,6) = SHLINT(INT,6) - R2 + Z2/D4
         END IF
  100 CONTINUE
      RETURN
      END
C  /* Deck thmint */
      SUBROUTINE THMINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT)
C
C     K.Ruud, Nov-97
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,10)
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         DX1 = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         DY1 = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         DZ1 = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
         DX2 = SHGTF*ODC(LVALA,LVALB,0,0,2,1)
         DY2 = SHGTF*ODC(MVALA,MVALB,0,0,2,2)
         DZ2 = SHGTF*ODC(NVALA,NVALB,0,0,2,3)
         DX3 = SHGTF*ODC(LVALA,LVALB,0,0,3,1)
         DY3 = SHGTF*ODC(MVALA,MVALB,0,0,3,2)
         DZ3 = SHGTF*ODC(NVALA,NVALB,0,0,3,3)
C
         INT = INT + 1
         SHLINT(INT,1)  = SHLINT(INT,1)  - DX3*SY0*SZ0
         SHLINT(INT,2)  = SHLINT(INT,2)  - DX2*DY1*SZ0
         SHLINT(INT,3)  = SHLINT(INT,3)  - DX2*SY0*DZ1
         SHLINT(INT,4)  = SHLINT(INT,4)  - DX1*DY2*SZ0
         SHLINT(INT,5)  = SHLINT(INT,5)  - DX1*DY1*DZ1
         SHLINT(INT,6)  = SHLINT(INT,6)  - DX1*SY0*DZ2
         SHLINT(INT,7)  = SHLINT(INT,7)  - SX0*DY3*SZ0
         SHLINT(INT,8)  = SHLINT(INT,8)  - SX0*DY2*DZ1
         SHLINT(INT,9)  = SHLINT(INT,9)  - SX0*DY1*DZ2
         SHLINT(INT,10) = SHLINT(INT,10) - SX0*SY0*DZ3
  100 CONTINUE
      RETURN
      END
C  /* Deck momint */
      SUBROUTINE MOMINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,EXPPI,ORIGIN,IORDER,WORK,
     &                  LWORK)
C
C     tuh 1989
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
      DIMENSION WORK(LWORK)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          ORIGIN(3), SHLINT(KCKTAB,(IORDER+1)*(IORDER+2)/2)
#include "onecom.h"
C
      KHMULT = 1
      KDX    = KHMULT + 3*((IORDER + 1)**2)
      KDY    = KDX    + IORDER + 1
      KDZ    = KDY    + IORDER + 1
      KLO    = KDZ    + IORDER + 1
      KMO    = KLO    + (((IORDER + 1)*(IORDER + 2)/2) + 1)/IRAT
      KNO    = KMO    + (((IORDER + 1)*(IORDER + 2)/2) + 1)/IRAT
      KLAST  = KNO    + (((IORDER + 1)*(IORDER + 2)/2) + 1)/IRAT
      IF (KLAST .GT. LWORK) CALL STOPIT('MOMINT',' ',KLAST,LWORK)
      CALL MOMIN1(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,CORPX,
     &            CORPY,CORPZ,EXPPI,ORIGIN,IORDER,WORK(KHMULT),
     &            WORK(KDX),WORK(KDY),WORK(KDZ),WORK(KLO),WORK(KMO),
     &            WORK(KNO))
      RETURN
      END
C  /* Deck momin1 */
      SUBROUTINE MOMIN1(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,EXPPI,ORIGIN,IORDER,HMULT,DX,
     &                  DY,DZ,LO,MO,NO)
C
C     tuh 1989
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          ORIGIN(3), SHLINT(KCKTAB,(IORDER+1)*(IORDER+2)/2),
     &          HMULT(3,0:IORDER,0:IORDER),
     &          DX(0:IORDER), DY(0:IORDER), DZ(0:IORDER),
     &          LO((IORDER + 1)*(IORDER + 2)/2),
     &          MO((IORDER + 1)*(IORDER + 2)/2),
     &          NO((IORDER + 1)*(IORDER + 2)/2)
#include "onecom.h"
#include "lmns.h"
C
C     Hermitian integrals
C     -------------------
C
      HMULT(1,0,0) = SHGTF
      HMULT(2,0,0) = SHGTF
      HMULT(3,0,0) = SHGTF
      DO 100 IO = 1, IORDER
         DO 200 IT = 0, IO
            HX = D0
            HY = D0
            HZ = D0
            IF (IT .GT. 0) THEN
               X_IT = IT
               HX = HX + X_IT*HMULT(1,IO-1,IT-1)
               HY = HY + X_IT*HMULT(2,IO-1,IT-1)
               HZ = HZ + X_IT*HMULT(3,IO-1,IT-1)
            END IF
            IF (IT .LE. IO-1) THEN
               HX = HX + (CORPX - ORIGIN(1))*HMULT(1,IO-1,IT)
               HY = HY + (CORPY - ORIGIN(2))*HMULT(2,IO-1,IT)
               HZ = HZ + (CORPZ - ORIGIN(3))*HMULT(3,IO-1,IT)
            END IF
            IF (IT .LE. IO-2) THEN
               HX = HX + DP5*EXPPI*HMULT(1,IO-1,IT+1)
               HY = HY + DP5*EXPPI*HMULT(2,IO-1,IT+1)
               HZ = HZ + DP5*EXPPI*HMULT(3,IO-1,IT+1)
            END IF
            HMULT(1,IO,IT) = HX
            HMULT(2,IO,IT) = HY
            HMULT(3,IO,IT) = HZ
  200    CONTINUE
  100 CONTINUE
C
C     Cartesian integrals
C
      INT = 0
      DO 300 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 300 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
C        One-dimensional integrals
C
         DO 400 IO = 0, IORDER
            DX(IO) = D0
            DY(IO) = D0
            DZ(IO) = D0
            DO 500 IT = 0, MIN(LVALA+LVALB,IO)
               DX(IO) = DX(IO)+ ODC(LVALA,LVALB,IT,0,0,1)*HMULT(1,IO,IT)
  500       CONTINUE
            DO 510 IU = 0, MIN(MVALA+MVALB,IO)
               DY(IO) = DY(IO)+ ODC(MVALA,MVALB,IU,0,0,2)*HMULT(2,IO,IU)
  510       CONTINUE
            DO 520 IV = 0, MIN(NVALA+NVALB,IO)
               DZ(IO) = DZ(IO)+ ODC(NVALA,NVALB,IV,0,0,3)*HMULT(3,IO,IV)
  520       CONTINUE
  400    CONTINUE
C
C        Three-dimensional integrals
C
         INT = INT + 1
         CALL LMNVAL(IORDER+1,(IORDER + 1)*(IORDER + 2)/2,LO,MO,NO)
         DO 600 I = 1, (IORDER + 1)*(IORDER + 2)/2
            SHLINT(INT,I) = SHLINT(INT,I)-DX(LO(I))*DY(MO(I))*DZ(NO(I))
  600    CONTINUE
  300 CONTINUE
      RETURN
      END
C  /* Deck dpvint */
      SUBROUTINE DPVINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT)
C
C     tuh Jan. 17 1990
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,3)
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         DX0 = SHGTF*ODC(LVALA,LVALB,0,1,0,1)
         DY0 = SHGTF*ODC(MVALA,MVALB,0,1,0,2)
         DZ0 = SHGTF*ODC(NVALA,NVALB,0,1,0,3)
         INT = INT + 1
         SHLINT(INT,1) = SHLINT(INT,1) - DX0*SY0*SZ0
         SHLINT(INT,2) = SHLINT(INT,2) - SX0*DY0*SZ0
         SHLINT(INT,3) = SHLINT(INT,3) - SX0*SY0*DZ0
  100 CONTINUE
      RETURN
      END
C  /* Deck dv2int */
      SUBROUTINE DV2INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT)
C
C     <d2/dx2>, <d2/dy2>, <d2/dz2>.
C
C     H.J.Aa.Jensen, Oct. 2004. Based on KININT and DPVINT.
C     For use with SQH2DO.
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (DP5 = 0.5 D00)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,3)
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
C    **********************************************
C    ***** CALCULATE second derivative INTEGRALS **
C    **********************************************
C
         X0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         Y0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         Z0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         X2 = SHGTF*ODC(LVALA,LVALB,0,2,0,1)
         Y2 = SHGTF*ODC(MVALA,MVALB,0,2,0,2)
         Z2 = SHGTF*ODC(NVALA,NVALB,0,2,0,3)
         INT = INT + 1
         SHLINT(INT,1) = SHLINT(INT,1) + X2*Y0*Z0
         SHLINT(INT,2) = SHLINT(INT,2) + X0*Y2*Z0
         SHLINT(INT,3) = SHLINT(INT,3) + X0*Y0*Z2
 100  CONTINUE
      RETURN
      END
C  /* Deck rstint */
      SUBROUTINE RSTINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT)
C
C     K.Ruud, October-97
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D2 = 2.0D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,6)
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         SX1 = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         SY1 = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         SZ1 = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
         DX0 = SHGTF*ODC(LVALA,LVALB,0,1,0,1)
         DY0 = SHGTF*ODC(MVALA,MVALB,0,1,0,2)
         DZ0 = SHGTF*ODC(NVALA,NVALB,0,1,0,3)
         DX1 = SHGTF*ODC(LVALA,LVALB,0,1,1,1)
         DY1 = SHGTF*ODC(MVALA,MVALB,0,1,1,2)
         DZ1 = SHGTF*ODC(NVALA,NVALB,0,1,1,3)
         INT = INT + 1
         OVLP = SX0*SY0*SZ0
         SHLINT(INT,1) = SHLINT(INT,1) + D2*DX1*SY0*SZ0 - OVLP
         SHLINT(INT,2) = SHLINT(INT,2) + (SX1*DY0 + DX0*SY1)*SZ0
         SHLINT(INT,3) = SHLINT(INT,3) + (SX1*DZ0 + DX0*SZ1)*SY0
         SHLINT(INT,4) = SHLINT(INT,4) + D2*SX0*DY1*SZ0 - OVLP
         SHLINT(INT,5) = SHLINT(INT,5) + SX0*(DY0*SZ1 + SY1*DZ0)
         SHLINT(INT,6) = SHLINT(INT,6) + D2*SX0*SY0*DZ1 - OVLP
  100 CONTINUE
      RETURN
      END
C  /* Deck frmint */
      SUBROUTINE FRMINT(SHLINT,NOPTYP,EXPP,XP,YP,ZP,SAAB,DOATOM)
C
C     tuh Feb 91
C
#include "implicit.h"
#include "priunit.h"
#include "pi.h"
#include "gfac.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "nuclei.h"
C
      PARAMETER (D0 = 0.0D0, D4 = 4.0D0, D3 = 3.0D0)
      PARAMETER (DFAC = D4*GFAC*PI/D3)
      LOGICAL DOATOM(*)
      DIMENSION SHLINT(KCKTAB,NOPTYP)
#include "onecom.h"
#include "lmns.h"
#include "symmet.h"
C
#include "ibtfun.h"
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
         INTTYP = 0
         INT = INT + 1
         DO 200 IREP = 0, MAXREP
            DO 300 IATOM = 1, NUCIND
            IF (DOATOM(IATOM)) THEN
            IF (IBTAND(IREP,ISTBNU(IATOM)).EQ.0) THEN
               INTTYP = INTTYP + 1
               FRM = D0
               DO 400 ISYMOP = 0, MAXOPR
               IF (IBTAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
                  XC = PT(IBTAND(ISYMAX(1,1),ISYMOP))*CORD(1,IATOM)
                  YC = PT(IBTAND(ISYMAX(2,1),ISYMOP))*CORD(2,IATOM)
                  ZC = PT(IBTAND(ISYMAX(3,1),ISYMOP))*CORD(3,IATOM)
                  RCP2 = (XC-XP)**2 + (YC-YP)**2 + (ZC-ZP)**2
                  FINT = PT(IBTAND(IREP,ISYMOP))*EXP(-EXPP*RCP2)
                  IF (LVALA .GT. 0) FINT = FINT*((XC - CORAX)**LVALA)
                  IF (MVALA .GT. 0) FINT = FINT*((YC - CORAY)**MVALA)
                  IF (NVALA .GT. 0) FINT = FINT*((ZC - CORAZ)**NVALA)
                  IF (LVALB .GT. 0) FINT = FINT*((XC - CORBX)**LVALB)
                  IF (MVALB .GT. 0) FINT = FINT*((YC - CORBY)**MVALB)
                  IF (NVALB .GT. 0) FINT = FINT*((ZC - CORBZ)**NVALB)
                  FRM = FRM + FINT
               END IF
  400          CONTINUE
               SHLINT(INT,INTTYP) = SHLINT(INT,INTTYP) - SAAB*FRM*DFAC
            END IF
            END IF
  300       CONTINUE
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
      SUBROUTINE FRMKIN(SHLINT,NOPTYP,EXPP,EXPB,XP,YP,ZP,SAAB,DOATOM)
C
C     K.Ruud, sep 2001
C
#include "implicit.h"
#include "priunit.h"
#include "pi.h"
#include "gfac.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
C
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0, D3 = 3.0D0, D4 = 4.0D0)
      PARAMETER (D1 = 1.0D0, DFAC = D2*GFAC*PI/D3)
      LOGICAL DOATOM(*)
      DIMENSION SHLINT(KCKTAB,NOPTYP)
#include "onecom.h"
#include "lmns.h"
#include "symmet.h"
#include "nuclei.h"
C
#include "ibtfun.h"
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
         INTTYP = 0
         INT = INT + 1
         DO 200 IREP = 0, MAXREP
            DO 300 IATOM = 1, NUCIND
            IF (DOATOM(IATOM)) THEN
            IF (IBTAND(IREP,ISTBNU(IATOM)).EQ.0) THEN
               INTTYP = INTTYP + 1
               FRM = D0
               DO 400 ISYMOP = 0, MAXOPR
               IF (IBTAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
                  XC = PT(IBTAND(ISYMAX(1,1),ISYMOP))*CORD(1,IATOM)
                  YC = PT(IBTAND(ISYMAX(2,1),ISYMOP))*CORD(2,IATOM)
                  ZC = PT(IBTAND(ISYMAX(3,1),ISYMOP))*CORD(3,IATOM)
                  RCP2 = (XC-XP)**2 + (YC-YP)**2 + (ZC-ZP)**2
                  FINT = PT(IBTAND(IREP,ISYMOP))*EXP(-EXPP*RCP2)
C
C     Contributions from center A
C
                  IF (LVALA .GT. 0) FINT = FINT*((XC - CORAX)**LVALA)
                  IF (MVALA .GT. 0) FINT = FINT*((YC - CORAY)**MVALA)
                  IF (NVALA .GT. 0) FINT = FINT*((ZC - CORAZ)**NVALA)
C
C     Undifferentiated contribution from center B
C
                  FINTYZ = D1
                  FINTXZ = D1
                  FINTXY = D1
                  IF (LVALB .GT. 0) THEN
                     FINTXY = FINTXY*((XC - CORBX)**LVALB)
                     FINTXZ = FINTXZ*((XC - CORBX)**LVALB)
                  END IF
                  IF (MVALB .GT. 0) THEN
                     FINTXY = FINTXY*((YC - CORBY)**MVALB)
                     FINTYZ = FINTYZ*((YC - CORBY)**MVALB)
                  END IF
                  IF (NVALB .GT. 0) THEN
                     FINTYZ = FINTYZ*((ZC - CORBZ)**NVALB)
                     FINTXZ = FINTXZ*((ZC - CORBZ)**NVALB)
                  END IF
C
C     Differentiated contribution from center B
C
                  IF (LVALB .GT. 0) THEN
                     XFINT = - D2*EXPB*((XC-CORBX)**LVALB)
                     XFINT = (2.0D0*LVALB + 1.0D0)*XFINT
                  ELSE
                     XFINT = - D2*EXPB
                  END IF
                  XFINT = XFINT
     &                  + D4*EXPB*EXPB*((XC - CORBX)**(LVALB + 2))
                  IF (LVALB .GT. 2) THEN
                    XFINT = XFINT +
     &                (LVALB*(LVALB-1.0D0)) * ((XC-CORBX)**(LVALB-2))
                  ELSE IF (LVALB .EQ. 2) THEN
                     XFINT = XFINT + D2
                  END IF
C
                  IF (MVALB .GT. 0) THEN
                     YFINT = - D2*EXPB*((YC-CORBY)**MVALB)
                     YFINT = (2.0D0*MVALB + 1.0D0)*YFINT
                  ELSE
                     YFINT = - D2*EXPB
                  END IF
                  YFINT = YFINT
     &                  + D4*EXPB*EXPB*((YC - CORBY)**(MVALB + 2))
                  IF (MVALB .GT. 2) THEN
                    YFINT = YFINT +
     &                (MVALB*(MVALB-1.0D0)) * ((YC-CORBY)**(MVALB-2))
                  ELSE IF (MVALB .EQ. 2) THEN
                     YFINT = YFINT + D2
                  END IF
C
                  IF (NVALB .GT. 0) THEN
                     ZFINT = - D2*EXPB*((ZC-CORBZ)**NVALB)
                     ZFINT = (2*NVALB + 1) * ZFINT
                  ELSE
                     ZFINT = - D2*EXPB
                  END IF
                  ZFINT = ZFINT
     &                  + D4*EXPB*EXPB*((ZC - CORBZ)**(NVALB + 2))
                  IF (NVALB .GT. 2) THEN
                     ZFINT = ZFINT +
     &                 (NVALB*(NVALB-1.0D0)) * ((ZC-CORBZ)**(NVALB-2))
                  ELSE IF (NVALB .EQ. 2) THEN
                     ZFINT = ZFINT + D2
                  END IF
                  FRM = FRM + FINT*(XFINT*FINTYZ + YFINT*FINTXZ
     &                            + ZFINT*FINTXY)
               END IF
  400          CONTINUE
               SHLINT(INT,INTTYP) = SHLINT(INT,INTTYP) - SAAB*FRM*DFAC
            END IF
            END IF
  300       CONTINUE
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck dwnint */
      SUBROUTINE DWNINT(SHLINT,EXPP,XP,YP,ZP,SAAB)
C
C 930309-Sheela Kirpekar: added MVLINT and DWNINT
C        (DWNINT based on FRMINT)
C

#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "codata.h"
C
      PARAMETER (D0 = 0.0D0)
      PARAMETER (DWNFAC = 0.5D0*PI*ALPHA2)
      DIMENSION SHLINT(KCKTAB)
#include "onecom.h"
#include "lmns.h"
#include "symmet.h"
#include "nuclei.h"
C
#include "ibtfun.h"
      DFAC = SAAB*DWNFAC
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
         INT = INT + 1
         DO 300 IATOM = 1, NUCIND
            FRM = D0
            DO 400 ISYMOP = 0, MAXOPR
            IF (IBTAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
               XC = PT(IBTAND(ISYMAX(1,1),ISYMOP))*CORD(1,IATOM)
               YC = PT(IBTAND(ISYMAX(2,1),ISYMOP))*CORD(2,IATOM)
               ZC = PT(IBTAND(ISYMAX(3,1),ISYMOP))*CORD(3,IATOM)
               RCP2 = (XC-XP)**2 + (YC-YP)**2 + (ZC-ZP)**2
               FINT = EXP(-EXPP*RCP2)
               IF (LVALA .GT. 0) FINT = FINT*((XC - CORAX)**LVALA)
               IF (MVALA .GT. 0) FINT = FINT*((YC - CORAY)**MVALA)
               IF (NVALA .GT. 0) FINT = FINT*((ZC - CORAZ)**NVALA)
               IF (LVALB .GT. 0) FINT = FINT*((XC - CORBX)**LVALB)
               IF (MVALB .GT. 0) FINT = FINT*((YC - CORBY)**MVALB)
               IF (NVALB .GT. 0) FINT = FINT*((ZC - CORBZ)**NVALB)
               FRM = FRM + FINT
            END IF
  400       CONTINUE
            SHLINT(INT) = SHLINT(INT) - FRM*DFAC*CHARGE(IATOM)
  300    CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck psoint */
      SUBROUTINE PSOINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT,
     &                  NOPTYP,NATOMC)
C
C     Paramagnetic spin-orbit integrals
C
C     tuh 8 Feb 1991
C
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,NOPTYP), AHGTF(*)
#include "onecom.h"
#include "lmns.h"
      INT = 1
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         MAXT = LVALA + LVALB + 1
         MAXU = MVALA + MVALB + 1
         MAXV = NVALA + NVALB + 1
         IADRAV = 1
         DO 200 IV = 0,MAXV
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            FV = ODC(NVALA,NVALB,IV,1,0,3)
            IADRAU = IADRAV
            DO 300 IU = 0,MAXU
               EU = ODC(MVALA,MVALB,IU,0,0,2)
               FU = ODC(MVALA,MVALB,IU,1,0,2)
               EE = EU*EV
               FE = FU*EV
               EF = EU*FV
               DO 400 IT = 0,MAXT
                  ET = ODC(LVALA,LVALB,IT,0,0,1)
                  FT = ODC(LVALA,LVALB,IT,1,0,1)
                  FEE = FT*EE
                  EFE = ET*FE
                  EEF = ET*EF
                  IOFF = 0
                  DO 500 IATOM = 1, NATOMC
                     IX = 3*(IATOM - 1) + 1
                     IY = 3*(IATOM - 1) + 2
                     IZ = 3*(IATOM - 1) + 3
                     AH0T = AHGTF(IOFF + IADRAU + IT + 1)
                     AH0U = AHGTF(IOFF + IADRAU + IT + ISTEPU)
                     AH0V = AHGTF(IOFF + IADRAU + IT + ISTEPV)
                     SHLINT(INT,IX) = SHLINT(INT,IX)+EFE*AH0V-EEF*AH0U
                     SHLINT(INT,IY) = SHLINT(INT,IY)+EEF*AH0T-FEE*AH0V
                     SHLINT(INT,IZ) = SHLINT(INT,IZ)+FEE*AH0U-EFE*AH0T
                     IOFF = IOFF + NAHGTF
  500             CONTINUE
  400          CONTINUE
               IADRAU = IADRAU + ISTEPU
  300       CONTINUE
            IADRAV = IADRAV + ISTEPV
  200    CONTINUE
         INT = INT + 1
  100 CONTINUE
      RETURN
      END
C
      SUBROUTINE PSOKIN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT,
     &                  NOPTYP,NATOMC)
C
C     Kinetic energy correction to paramagnetic spin-orbit integrals
C
C     K.Ruud, Sep.2001, based on PSOINT
C
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"
      PARAMETER (DP5 = 0.50D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,NOPTYP), AHGTF(*)
#include "onecom.h"
#include "lmns.h"
      INT = 1
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         MAXT = LVALA + LVALB + 3
         MAXU = MVALA + MVALB + 3
         MAXV = NVALA + NVALB + 3
         IADRAV = 1
         DO 200 IV = 0, MAXV
            EV  = ODC(NVALA,NVALB,IV,0,0,3)
            FV1 = ODC(NVALA,NVALB,IV,1,0,3)
            FV2 = ODC(NVALA,NVALB,IV,2,0,3)
            FV3 = ODC(NVALA,NVALB,IV,3,0,3)
            IADRAU = IADRAV
            DO 300 IU = 0, MAXU
               EU  = ODC(MVALA,MVALB,IU,0,0,2)
               FU1 = ODC(MVALA,MVALB,IU,1,0,2)
               FU2 = ODC(MVALA,MVALB,IU,2,0,2)
               FU3 = ODC(MVALA,MVALB,IU,3,0,2)
C
               EE   = EU*EV
               F1E  = FU1*EV
               F2E  = FU2*EV
               F3E  = FU3*EV
               EF1  = EU*FV1
               F2F1 = FU2*FV1
               F1F2 = FU1*FV2
               EF2  = EU*FV2
               EF3  = EU*FV3
               DO 400 IT = 0, MAXT
                  ET  = ODC(LVALA,LVALB,IT,0,0,1)
                  FT1 = ODC(LVALA,LVALB,IT,1,0,1)
                  FT2 = ODC(LVALA,LVALB,IT,2,0,1)
                  FT3 = ODC(LVALA,LVALB,IT,3,0,1)
C
                  F1F2E = FT1*F2E
                  F3EE  = FT3*EE
                  EF3E  = ET*F3E
                  F2F1E = FT2*F1E
                  F2EF1 = FT2*EF1
                  EF2F1 = ET*F2F1
                  EF1F2 = ET*F1F2
                  F1EF2 = FT1*EF2
                  EEF3  = ET*EF3
                  IOFF = 0
                  DO 500 IATOM = 1, NATOMC
                     IX = 3*(IATOM - 1) + 1
                     IY = 3*(IATOM - 1) + 2
                     IZ = 3*(IATOM - 1) + 3
                     AH0T = AHGTF(IOFF + IADRAU + IT + 1)
                     AH0U = AHGTF(IOFF + IADRAU + IT + ISTEPU)
                     AH0V = AHGTF(IOFF + IADRAU + IT + ISTEPV)
ckr                     write (2,*) 'f2ef1,ef2f1,eef3,ah0u',
ckr     &                    f2ef1,ef2f1,eef3,ah0u
ckr                     write (2,*) 'f2f1e,ef3e,ef1f2,ah0v',
ckr     &                    f2f1e,ef3e,ef1f2,ah0v
ckr                     write (2,*) 'ituv,iatom,istepu,istepv,ioff',
ckr     &                    lvala,mvala,nvala,lvalb,mvalb,nvalb,it,iu,iv,
ckr     &                    iatom,istepu,istepv,ioff
                     SHLINT(INT,IX) = SHLINT(INT,IX) + DP5*
     &                    ((F2EF1 + EF2F1 + EEF3)*AH0U
     &                   - (F2F1E + EF3E + EF1F2)*AH0V)
ckr                     write(2,*) 'int,ix',int,ix,shlint(int,ix)
                     SHLINT(INT,IY) = SHLINT(INT,IY) + DP5*
     &                    ((F3EE + F1F2E + F1EF2)*AH0V
     &                   - (F2EF1 + EF2F1 + EEF3)*AH0T)
                     SHLINT(INT,IZ) = SHLINT(INT,IZ) + DP5*
     &                    ((F2F1E + EF3E + EF1F2)*AH0T
     &                   - (F3EE + F1F2E + F1EF2)*AH0U)
                     IOFF = IOFF + NAHGTF
  500             CONTINUE
  400          CONTINUE
               IADRAU = IADRAU + ISTEPU
  300       CONTINUE
            IADRAV = IADRAV + ISTEPV
  200    CONTINUE
         INT = INT + 1
  100 CONTINUE
      RETURN
      END
      SUBROUTINE PSOZIN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,
     &                  SHLINT,NOPTYP,NATOMC)
C
C     Orbital Zeeman correction to PSO integrals
C
C     K. Ruud and P. Manninen, February 2004
C
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"
      PARAMETER (D2INV = 0.50D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,NOPTYP), AHGTF(*)
#include "onecom.h"
#include "lmns.h"
      INT = 1
      DO 100 ICOMPA = 1, KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         MAXT = LVALA + LVALB + 3
         MAXU = MVALA + MVALB + 3
         MAXV = NVALA + NVALB + 3
         IADRAV = 1
         DO 200 IV = 0, MAXV
            DV = ODC(NVALA,NVALB,IV,0,1,3)
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            FV = ODC(NVALA,NVALB,IV,1,0,3)
            GV = ODC(NVALA,NVALB,IV,1,1,3)
            HV = ODC(NVALA,NVALB,IV,2,0,3)
            IADRAU = IADRAV
            DO 300 IU = 0, MAXU
               DU = ODC(MVALA,MVALB,IU,0,1,2)
               EU = ODC(MVALA,MVALB,IU,0,0,2)
               FU = ODC(MVALA,MVALB,IU,1,0,2)
               GU = ODC(MVALA,MVALB,IU,1,1,2)
               HU = ODC(MVALA,MVALB,IU,2,0,2)
               DO 400 IT = 0, MAXT
                  DT = ODC(LVALA,LVALB,IT,0,1,1)
                  ET = ODC(LVALA,LVALB,IT,0,0,1)
                  FT = ODC(LVALA,LVALB,IT,1,0,1)
                  GT = ODC(LVALA,LVALB,IT,1,1,1)
                  HT = ODC(LVALA,LVALB,IT,2,0,1)
                  DHE = DT*HU*EV
                  HDE = HT*DU*EV
                  DEH = DT*EU*HV
                  HED = HT*EU*DV
                  EDH = ET*DU*HV
                  EHD = ET*HU*DV
                  FEE = FT*EU*EV
                  EFE = ET*FU*EV
                  EEF = ET*EU*FV
                  FGE = FT*GU*EV
                  GFE = GT*FU*EV
                  FEG = FT*EU*GV
                  GEF = GT*EU*FV
                  EFG = ET*FU*GV
                  EGF = ET*GU*FV
                  FFD = FT*FU*DV
                  FDF = FT*DU*FV
                  DFF = DT*FU*FV
                  IOFF = 0
                  DO 500 IATOM = 1, NATOMC
                     IXX = 9*(IATOM - 1) + 1
                     IYX = 9*(IATOM - 1) + 2
                     IZX = 9*(IATOM - 1) + 3
                     IXY = 9*(IATOM - 1) + 4
                     IYY = 9*(IATOM - 1) + 5
                     IZY = 9*(IATOM - 1) + 6
                     IXZ = 9*(IATOM - 1) + 7
                     IYZ = 9*(IATOM - 1) + 8
                     IZZ = 9*(IATOM - 1) + 9
                     AH0T = AHGTF(IOFF + IADRAU + IT + 1)
                     AH0U = AHGTF(IOFF + IADRAU + IT + ISTEPU)
                     AH0V = AHGTF(IOFF + IADRAU + IT + ISTEPV)
C
                     SHLINT(INT,IXX) = SHLINT(INT,IXX)
     &                               + ((EFE - EFG + EDH)*AH0U
     &                               + (EEF - EGF + EHD)*AH0V)
                     SHLINT(INT,IXY) = SHLINT(INT,IXY)
     &                               + ((EFG - EFE - EDH)*AH0T
     &                               + (FDF - FFD)*AH0V)
                     SHLINT(INT,IXZ) = SHLINT(INT,IXZ)
     &                               + ((EGF - EEF - EHD)*AH0T
     &                               + (FFD - FDF)*AH0U)
                     SHLINT(INT,IYX) = SHLINT(INT,IYX)
     &                               + ((FEG - FEE - DEH)*AH0U
     &                               + (DFF - FFD)*AH0V)
                     SHLINT(INT,IYY) = SHLINT(INT,IYY)
     &                               + ((EEF - GEF + HED)*AH0V
     &                               + (FEE - FEG + DEH)*AH0T)
                     SHLINT(INT,IYZ) = SHLINT(INT,IYZ)
     &                               + ((GEF - EEF - HED)*AH0U
     &                               + (FFD - DFF)*AH0T)
                     SHLINT(INT,IZX) = SHLINT(INT,IZX)
     &                               + ((FGE - FEE - DHE)*AH0V
     &                               + (DFF - FDF)*AH0U)
                     SHLINT(INT,IZY) = SHLINT(INT,IZY)
     &                               + ((GFE - EFE - HDE)*AH0V
     &                               + (FDF - DFF)*AH0T)
                     SHLINT(INT,IZZ) = SHLINT(INT,IZZ)
     &                               + ((FEE - FGE + DHE)*AH0T
     &                               + (EFE - GFE + HDE)*AH0U)
                     IOFF = IOFF + NAHGTF
ckr                  write(2,*) 'shlint(int,ixy)',shlint(int,ixy),ixy,int
 500              CONTINUE
 400           CONTINUE
               IADRAU = IADRAU + ISTEPU
 300        CONTINUE
            IADRAV = IADRAV + ISTEPV
 200     CONTINUE
         INT = INT + 1
 100  CONTINUE
      RETURN
      END
C  /* Deck g1hint */
      SUBROUTINE G1HINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  AHGTF,NATOMC,NOPTYP,NCENTC)
C
C     K.Ruud, May 1998
C     Rewritten for symmetry tuh 00
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
      PARAMETER (D2INV = 0.50D0)
      INTEGER XA, YA, ZA, XB, YB, ZB, XC, YC, ZC
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,NOPTYP,3), AHGTF(*), NCENTC(*)
#include "onecom.h"
#include "nuclei.h"
#include "lmns.h"
C
      XA = 3*(ICENTA - 1) + 1
      YA = 3*(ICENTA - 1) + 2
      ZA = 3*(ICENTA - 1) + 3
      XB = 3*(ICENTB - 1) + 1
      YB = 3*(ICENTB - 1) + 2
      ZB = 3*(ICENTB - 1) + 3
C
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1,KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
         SX0  = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0  = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0  = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         DX1  = SHGTF*ODC(LVALA,LVALB,0,1,0,1)
         DY1  = SHGTF*ODC(MVALA,MVALB,0,1,0,2)
         DZ1  = SHGTF*ODC(NVALA,NVALB,0,1,0,3)
         DX2  = SHGTF*ODC(LVALA,LVALB,0,2,0,1)
         DY2  = SHGTF*ODC(MVALA,MVALB,0,2,0,2)
         DZ2  = SHGTF*ODC(NVALA,NVALB,0,2,0,3)
         DX3  = SHGTF*ODC(LVALA,LVALB,0,3,0,1)
         DY3  = SHGTF*ODC(MVALA,MVALB,0,3,0,2)
         DZ3  = SHGTF*ODC(NVALA,NVALB,0,3,0,3)
C
C        Kinetic energy contribution
C
         TX = -D2INV*(DX3*SY0*SZ0 + DX1*DY2*SZ0 + DX1*SY0*DZ2)
         TY = -D2INV*(DX2*DY1*SZ0 + SX0*DY3*SZ0 + SX0*DY1*DZ2)
         TZ = -D2INV*(DX2*SY0*DZ1 + SX0*DY2*DZ1 + SX0*SY0*DZ3)
C
         SHLINT(INT,XA,1) = SHLINT(INT,XA,1) + TX
         SHLINT(INT,YA,1) = SHLINT(INT,YA,1) + TY
         SHLINT(INT,ZA,1) = SHLINT(INT,ZA,1) + TZ
         SHLINT(INT,XB,2) = SHLINT(INT,XB,2) - TX
         SHLINT(INT,YB,2) = SHLINT(INT,YB,2) - TY
         SHLINT(INT,ZB,2) = SHLINT(INT,ZB,2) - TZ
C
C        Nuclear attraction contribution
C
         MAXT = LVALA + LVALB + 1
         MAXU = MVALA + MVALB + 1
         MAXV = NVALA + NVALB + 1
         IADRAV = 1
         DO 200 IV = 0, MAXV
            IADRAU = IADRAV
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            FV = ODC(NVALA,NVALB,IV,1,0,3)
            IUMAX = MAXU
            DO 300 IU = 0, IUMAX
               EU = ODC(MVALA,MVALB,IU,0,0,2)
               FU = ODC(MVALA,MVALB,IU,1,0,2)
               ITMAX = MAXT
               DO 400 IT = 0, ITMAX
                  ET = ODC(LVALA,LVALB,IT,0,0,1)
                  FT = ODC(LVALA,LVALB,IT,1,0,1)
                  EE = ET*EU*EV
                  FX = FT*EU*EV
                  FY = ET*FU*EV
                  FZ = ET*EU*FV
C
                  IOFF = 0
                  DO 500 IATOM = 1, NATOMC
C
                     XC = 3*(NCENTC(IATOM) - 1) + 1
                     YC = 3*(NCENTC(IATOM) - 1) + 2
                     ZC = 3*(NCENTC(IATOM) - 1) + 3
C
                     ADERX =  FX*AHGTF(IOFF + IADRAU + IT)
                     ADERY =  FY*AHGTF(IOFF + IADRAU + IT)
                     ADERZ =  FZ*AHGTF(IOFF + IADRAU + IT)
                     CDERX = -EE*AHGTF(IOFF + IADRAU + IT + 1)
                     CDERY = -EE*AHGTF(IOFF + IADRAU + IT + ISTEPU)
                     CDERZ = -EE*AHGTF(IOFF + IADRAU + IT + ISTEPV)
C
                     SHLINT(INT,XA,1) = SHLINT(INT,XA,1) + ADERX
                     SHLINT(INT,YA,1) = SHLINT(INT,YA,1) + ADERY
                     SHLINT(INT,ZA,1) = SHLINT(INT,ZA,1) + ADERZ
                     SHLINT(INT,XB,2) = SHLINT(INT,XB,2) - ADERX - CDERX
                     SHLINT(INT,YB,2) = SHLINT(INT,YB,2) - ADERY - CDERY
                     SHLINT(INT,ZB,2) = SHLINT(INT,ZB,2) - ADERZ - CDERZ
                     SHLINT(INT,XC,3) = SHLINT(INT,XC,3) + CDERX
                     SHLINT(INT,YC,3) = SHLINT(INT,YC,3) + CDERY
                     SHLINT(INT,ZC,3) = SHLINT(INT,ZC,3) + CDERZ
C
                     IOFF = IOFF + NAHGTF
 500              CONTINUE
C
 400           CONTINUE
               IADRAU = IADRAU + ISTEPU
 300        CONTINUE
            IADRAV = IADRAV + ISTEPV
 200     CONTINUE
 100  CONTINUE
      RETURN
      END
C  /* Deck nsttra */
      SUBROUTINE NSTTRA(SHLINT,SOINT,FACINT,JSYMC,JCENTC,NATOMC,ISYMOP,
     &                  NELMNT,NOPTYP,ANTI,IORBA,IORBB,MULTA,MULTB,
     &                  NBAST,INTTYP,IPRINT,INTADR,EXPVAL,DMAT,EXP1EL)
C
#include "implicit.h"
#include "dummy.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
      LOGICAL FULMAT, ANTI, EXP1EL
      DIMENSION SHLINT(KCKTAB,NOPTYP), SOINT(NELMNT,NOPTYP),
     &          FACINT(NATOMC), JSYMC(NATOMC), JCENTC(NATOMC),
     &          INTADR(*)
      DIMENSION EXPVAL(*), DMAT(*)
#include "onecom.h"
#include "symmet.h"
#include "ibtfun.h"
      FULMAT = .TRUE.
      DO 100 IREPC = 0, MAXREP
         DO 200 IATOMC = 1, NATOMC
            ICENTC = JCENTC(IATOMC)
            ISYMC  = JSYMC(IATOMC)
            FACTOR = FACINT(IATOMC)
            DO 300 LCOOR = 1, 3
               ISCOOR = IPTCNT(3*(ICENTC - 1) + LCOOR,IREPC,2)
               IF (ISCOOR .GT. 0) THEN
                  DO 400 ICOOR = 1, 3
                     IREPCD = IBTXOR(IREPC,ISYMAX(ICOOR,2))
                     ISYMCA = IBTXOR(ISYMAX(ICOOR,2),ISYMAX(LCOOR,2))
cLig >> added how to set the symmetry for RPSO integrals
                     IF(INTTYP.EQ.82) THEN
                       IREPCD = IBTXOR(IREPC,ISYMAX(ICOOR,1))
                       ISYMCA = IBTXOR(ISYMAX(ICOOR,1),ISYMAX(LCOOR,2))
                     ENDIF
cLig <<
                     FACSYM = - FACTOR*PT(IBTAND(ISYMCA,ISYMC))
     &                                *PT(IBTAND(IREPCD,ISYMC))
                     ITYP = 9*(IATOMC - 1) + 3*(LCOOR - 1) + ICOOR
                     IADR = INTADR(3*(ISCOOR - 1) + ICOOR)
                     IF (EXP1EL) THEN
                        IF (IREPCD .EQ. 0) CALL SYM1EV(SHLINT(1,ITYP),
     &                                    DMAT,EXPVAL(IADR),FACSYM,
     &                                    MULTA,MULTB,IPRINT)
                     ELSE IF (INTTYP .EQ. 26 .OR. INTTYP .EQ. 27
     &                   .OR. INTTYP .EQ. 74
     &                   .OR. INTTYP .EQ. 82 .OR. INTTYP .EQ. 93
     &                   .OR. INTTYP .EQ. 97) THEN
cLig <> added the call to SYMSQR for INTTYP = 82
                        CALL SYMSQR(SHLINT(1,ITYP),SOINT(1,IADR),IREPCD,
     &                              ISYMOP,IORBA,IORBB,FACSYM,NBAST,
     &                              IPRINT)
                     ELSE IF (IREPCD .EQ. 0) THEN
                        CALL SYM1S(SHLINT(1,ITYP),SOINT(1,IADR),ISYMOP,
     &                             MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                             FACSYM,LDIAG,FULMAT,DUM,IDUM,IPRINT)
                     ELSE
                        CALL SYM1N(SHLINT(1,ITYP),SOINT(1,IADR),IREPCD,
     &                             ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,
     &                             KHKTB,FACSYM,LDIAG,FULMAT,ANTI,DUM,
     &                             IDUM,IPRINT)
                     END IF
 400              CONTINUE
               END IF
 300        CONTINUE
 200     CONTINUE
 100   CONTINUE
       RETURN
       END
C  /* Deck psotra */
      SUBROUTINE PSOTRA(SHLINT,SOINT,FACINT,JSYMC,JCENTC,NATOMC,ISYMOP,
     &                  NELMNT,NOPTYP,ANTI,IORBA,IORBB,NBAST,INTTYP,
     &                  INTADR,IPRINT)
#include "implicit.h"
#include "dummy.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
      LOGICAL FULMAT, ANTI
      DIMENSION SHLINT(KCKTAB,NOPTYP), SOINT(NELMNT,NOPTYP),
     &          FACINT(NATOMC), JSYMC(NATOMC), JCENTC(NATOMC),
     &          INTADR(*)
#include "onecom.h"
#include "symmet.h"
#include "ibtfun.h"
      FULMAT = .TRUE.
      DO 100 IREPC = 0, MAXREP
         DO 200 IATOMC = 1, NATOMC
            ICENTC = JCENTC(IATOMC)
            ISYMC  = JSYMC (IATOMC)
            FACTOR = FACINT(IATOMC)
            DO 300 ICOOR = 1, 3
               ISCOOR = IPTCNT(3*(ICENTC - 1) + ICOOR,IREPC,2)
               IF (ISCOOR .GT. 0) THEN
                  ISYMCR = ISYMAX(ICOOR,2)
                  FACSYM = - FACTOR*PT(IBTAND(ISYMCR,ISYMC))
     &                             *PT(IBTAND(IREPC ,ISYMC))
                  ITYP = 3*(IATOMC - 1) + ICOOR
                  IADR = INTADR(ISCOOR)
                  IF (INTTYP .EQ. 92) THEN
                     CALL SYMSQR(SHLINT(1,ITYP),SOINT(1,IADR),IREPC,
     &                           ISYMOP,IORBA,IORBB,FACSYM,NBAST,IPRINT)
                  ELSE
                     IF (IREPC .EQ. 0) THEN
                        CALL SYM1S(SHLINT(1,ITYP),SOINT(1,IADR),ISYMOP,
     &                             MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                              FACSYM,LDIAG,FULMAT,DUM,IDUM,IPRINT)
                     ELSE
                        CALL SYM1N(SHLINT(1,ITYP),SOINT(1,IADR),IREPC,
     &                             ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,
     &                             KHKTB,FACSYM,LDIAG,FULMAT,ANTI,DUM,
     &                             IDUM,IPRINT)
                     END IF
                  END IF
               END IF
  300       CONTINUE
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck hdbtra */
      SUBROUTINE HDBTRA(SHLINT,SOINT,ISYMOP,NELMNT,NOPTYP,IORBA,IORBB,
     &                  NBAST,INTADR,IPRINT)
#include "implicit.h"
#include "maxaqn.h"
#include "maxmom.h"
#include "mxcent.h"
#include "maxorb.h"
      INTEGER B, X
      DIMENSION SHLINT(KCKTAB,NOPTYP), SOINT(NELMNT,NOPTYP),
     &          INTADR(*)
#include "onecom.h"
#include "symmet.h"
#include "ibtfun.h"
      ICOMP = 0
      DO 100 X = 1, 3
         DO 200 B = 1, 3
            ICOMP = ICOMP + 1
            DO 300 IREPO = 0, MAXREP
               ISCOOR = IPTCNT(3*(NCENTA - 1) + X,IREPO,1)
               IF (ISCOOR .GT. 0) THEN
                  ITYP   = INTADR(3*(ISCOOR - 1) + B)
                  IREPCL = IBTXOR(IREPO,ISYMAX(B,2))
                  FACSYM = HKAB*PT(IBTAND(IREPCL,ISYMOP))
                  CALL SYMSQR(SHLINT(1,ICOMP),SOINT(1,ITYP),IREPCL,
     &                        ISYMOP,IORBA,IORBB,FACSYM,NBAST,IPRINT)
               END IF
 300        CONTINUE
 200     CONTINUE
 100  CONTINUE
      RETURN
      END
C  /* Deck dpgtra */
      SUBROUTINE DPGTRA(SHLINT,SOINT,ISYMOP,NELMNT,NOPTYP,ANTI,
     &                  INTADR,IPRINT)
#include "implicit.h"
#include "dummy.h"
#include "maxaqn.h"
#include "maxmom.h"
#include "mxcent.h"
#include "maxorb.h"
      INTEGER B, X
      LOGICAL FULMAT, ANTI
      DIMENSION SHLINT(KCKTAB,18), SOINT(NELMNT,NOPTYP),
     &          INTADR(*)
#include "onecom.h"
#include "symmet.h"
#include "ibtfun.h"
      FULMAT = .TRUE.
      ICOMP = 0

      ! loop over cartesian displacements of the atoms
      DO 100 X = 1, 3
         JCOMP = 0

         ! loop over components of the dipole operator
         DO 200 B = 1, 3
            ICOMP = ICOMP + 1
            JCOMP = JCOMP + 1
            IREPC  = ISYMAX(B,1) ! irrep of dipole component

            ! loop over possible irreps of symmetry coordinates
            DO IREP = 0, MAXREP
               ISCOOR = IPTCNT(3*(NCENTA - 1) + X,IREP,1)
               IF (ISCOOR .GT. 0) THEN
                  ITYP   = INTADR(3*(ISCOOR - 1) + JCOMP)
                  IREPO = IBTXOR(IREP,IREPC) ! irrep of operator
ckr                  FACSYM = HKAB*PT(IBTAND(IREPO,ISYMOP))
C                  FACSYM = HKAB*PT(IBTAND(IREPC,ISYMOP))
C     &                         *PT(IBTAND(IREP,ISYMOP))
ckr                  FACSYM = HKAB
                  FACSYM = HKAB
                  IF (IREPO .EQ. 0) THEN
                     CALL SYM1S(SHLINT(1,ICOMP),SOINT(1,ITYP),ISYMOP,
     &                          MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                          FACSYM,LDIAG,FULMAT,DUM,IDUM,IPRINT)
                  ELSE
                     CALL SYM1N(SHLINT(1,ICOMP),SOINT(1,ITYP),IREPO,
     &                          ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,
     &                          KHKTB,FACSYM,LDIAG,FULMAT,ANTI,DUM,
     &                          IDUM,IPRINT)
                  END IF
               END IF
               ISCOOR = IPTCNT(3*(NCENTB - 1) + X,IREP,1)
               IF (ISCOOR .GT. 0) THEN
                  ITYP   = INTADR(3*(ISCOOR - 1) + JCOMP)
                  IREPO = IBTXOR(IREP,IREPC) ! irrep of operator
ckr                  FACSYM = HKAB*PT(IBTAND(IREPO,ISYMOP))
C                  FACSYM = HKAB*PT(IBTAND(IREPC,ISYMOP))
C    &                          *PT(IBTAND(IREP,ISYMOP))
                   FACSYM = HKAB*PT(IBTAND(ISYMAX(X,1),ISYMOP))
     &                          *PT(IBTAND(IREP,ISYMOP))
                  IF (IREPO .EQ. 0) THEN
                     CALL SYM1S(SHLINT(1,ICOMP + 9),SOINT(1,ITYP),
     &                          ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,
     &                          KHKTB,FACSYM,LDIAG,FULMAT,DUM,IDUM,
     &                          IPRINT)
                  ELSE
                     CALL SYM1N(SHLINT(1,ICOMP + 9),SOINT(1,ITYP),
     &                          IREPO,ISYMOP,MULA,MULB,NHKTA,NHKTB,
     &                          KHKTA,KHKTB,FACSYM,LDIAG,FULMAT,ANTI,
     &                          DUM,IDUM,IPRINT)
                  END IF
               END IF
            END DO
 200     CONTINUE
 100  CONTINUE
      RETURN
      END
C  /* Deck dphtra */
      SUBROUTINE DPHTRA(SHLINT,SOINT,ISYMOP,NELMNT,NOPTYP,ANTI,
     &                  INTADR,IPRINT)
#include "implicit.h"
#include "dummy.h"
#include "maxaqn.h"
#include "maxmom.h"
#include "mxcent.h"
#include "maxorb.h"
      INTEGER B, X, Y
      LOGICAL FULMAT, ANTI
      DIMENSION SHLINT(KCKTAB,18), SOINT(NELMNT,NOPTYP),
     &          INTADR(*)
#include "onecom.h"
#include "nuclei.h"
! nuclei.h needed for NUCDEP
#include "symmet.h"
#include "ibtfun.h"
      FULMAT = .TRUE.
      ICOMP = 0

      ! loop over cartesian displacements of the atoms
      DO X = 1, 3
         DO Y = 1, X
            JCOMP = 0

         ! loop over components of the dipole operator
            DO 200 B = 1, 3
               ICOMP = ICOMP + 1
               JCOMP = JCOMP + 1
               IREPC  = ISYMAX(B,1) ! irrep of dipole component

            ! loop over possible irreps of symmetry coordinates
               DO IREP = 0, MAXREP
                  ISCOOR = IPTCNT(3*(NCENTA - 1) + X,IREP,1)
C
C     We first distribute the diagonal blocks of the Hessian-type matrix
C     Centre 1
C
                  IF (ISCOOR .GT. 0) THEN
                     ITYP   = INTADR(9*NUCDEP*(ISCOOR - 1)
     &                                       + ISCOOR + JCOMP)
                     IREPO = IREPC ! irrep of operator=irrep of dipole op.
                     FACSYM = HKAB
                     IF (IREPO .EQ. 0) THEN
                        CALL SYM1S(SHLINT(1,ICOMP),SOINT(1,ITYP),ISYMOP,
     &                             MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                             FACSYM,LDIAG,FULMAT,DUM,IDUM,IPRINT)
                     ELSE
                        CALL SYM1N(SHLINT(1,ICOMP),SOINT(1,ITYP),IREPO,
     &                             ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,
     &                             KHKTB,FACSYM,LDIAG,FULMAT,ANTI,DUM,
     &                             IDUM,IPRINT)
                     END IF
                  END IF
                  JSCOOR = IPTCNT(3*(NCENTB - 1) + X,IREP,1)
C     Centre 2
                  IF (JSCOOR .GT. 0) THEN
                     ITYP   = INTADR(9*NUCDEP*(JSCOOR - 1)
     &                                       + JSCOOR + JCOMP)
                     IREPO = IREPC ! irrep of operator=irrep of dipole op.
                     FACSYM = HKAB

ckr                  FACSYM = HKAB*PT(IBTAND(IREPO,ISYMOP))
C                  FACSYM = HKAB*PT(IBTAND(IREPC,ISYMOP))
C    &                          *PT(IBTAND(IREP,ISYMOP))
                     FACSYM = HKAB*PT(IBTAND(IBTXOR(ISYMAX(X,1),
     &                                ISYMAX(Y,1)),ISYMOP))
     &                            *PT(IBTAND(IREP,ISYMOP))
                     IF (IREPO .EQ. 0) THEN
                        CALL SYM1S(SHLINT(1,ICOMP + 18),SOINT(1,ITYP),
     &                             ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,
     &                             KHKTB,FACSYM,LDIAG,FULMAT,DUM,IDUM,
     &                             IPRINT)
                     ELSE
                        CALL SYM1N(SHLINT(1,ICOMP + 9),SOINT(1,ITYP),
     &                             IREPO,ISYMOP,MULA,MULB,NHKTA,NHKTB,
     &                             KHKTA,KHKTB,FACSYM,LDIAG,FULMAT,ANTI,
     &                             DUM,IDUM,IPRINT)
                     END IF
                  END IF
C
C     Then we do the off-diagonal blocks
C
Ckr                  DO JREP = 0, MAXREP
Ckr                     JSCOOR = IPTCNT(3*(NCENTB - 1) + Y,JREP,1)
Ckr                     IF (JSCOOR .GT. 0) THEN 
Cr                        IREPN = IBTXOR(IREP,JREP)
Ckr                        IREPO = IBTXOR(IREPN,IREPC)
               END DO
 200        CONTINUE
         END DO
      END DO
      RETURN
      END
C  /* Deck qugtra */
      SUBROUTINE QUGTRA(SHLINT,SOINT,ISYMOP,NELMNT,NOPTYP,ANTI,
     &                  INTADR,IPRINT)
#include "implicit.h"
#include "dummy.h"
#include "maxaqn.h"
#include "maxmom.h"
#include "mxcent.h"
#include "maxorb.h"
      INTEGER B, X, C
      LOGICAL FULMAT, ANTI
      DIMENSION SHLINT(KCKTAB,36), SOINT(NELMNT,NOPTYP),
     &          INTADR(*)
#include "onecom.h"
#include "symmet.h"
#include "ibtfun.h"
      FULMAT = .TRUE.
      ICOMP = 0
      DO 100 X = 1, 3
         JCOMP = 0
         DO 200 B = 1, 3
         DO 200 C = B, 3
            ICOMP = ICOMP + 1
            JCOMP = JCOMP + 1
            IREPC  = IBTXOR(ISYMAX(B,1),ISYMAX(C,1))
            DO IREPO = 0, MAXREP
               ISCOOR = IPTCNT(3*(NCENTA - 1) + X,IREPO,1)
               IF (ISCOOR .GT. 0) THEN
                  ITYP   = INTADR(6*(ISCOOR - 1) + JCOMP)
                  IREPCL = IBTXOR(IREPO,IREPC)
ckr                  FACSYM = HKAB*PT(IBTAND(IREPCL,ISYMOP))
                  FACSYM = HKAB*PT(IBTAND(IREPC,ISYMOP))
     &                         *PT(IBTAND(IREPO,ISYMOP))
ckr                  FACSYM = HKAB
                  IF (IREPO .EQ. 0) THEN
                     CALL SYM1S(SHLINT(1,ICOMP),SOINT(1,ITYP),ISYMOP,
     &                          MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                          FACSYM,LDIAG,FULMAT,DUM,IDUM,IPRINT)
                  ELSE
                     CALL SYM1N(SHLINT(1,ICOMP),SOINT(1,ITYP),IREPO,
     &                          ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,
     &                          KHKTB,FACSYM,LDIAG,FULMAT,ANTI,DUM,
     &                          IDUM,IPRINT)
                  END IF
               END IF
               ISCOOR = IPTCNT(3*(NCENTB - 1) + X,IREPO,1)
               IF (ISCOOR .GT. 0) THEN
                  ITYP   = INTADR(6*(ISCOOR - 1) + JCOMP)
                  IREPCL = IBTXOR(IREPO,IREPC)
ckr                  FACSYM = HKAB*PT(IBTAND(IREPCL,ISYMOP))
                  FACSYM = HKAB*PT(IBTAND(IREPC,ISYMOP))
     &                         *PT(IBTAND(IREPO,ISYMOP))
                  IF (IREPO .EQ. 0) THEN
                     CALL SYM1S(SHLINT(1,ICOMP + 18),SOINT(1,ITYP),
     &                          ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,
     &                          KHKTB,FACSYM,LDIAG,FULMAT,DUM,IDUM,
     &                          IPRINT)
                  ELSE
                     CALL SYM1N(SHLINT(1,ICOMP + 18),SOINT(1,ITYP),
     &                          IREPO,ISYMOP,MULA,MULB,NHKTA,NHKTB,
     &                          KHKTA,KHKTB,FACSYM,LDIAG,FULMAT,ANTI,
     &                          DUM,IDUM,IPRINT)
                  END IF
               END IF
            END DO
 200     CONTINUE
 100  CONTINUE
      RETURN
      END
C  /* Deck ocgtra */
      SUBROUTINE OCGTRA(SHLINT,SOINT,ISYMOP,NELMNT,NOPTYP,ANTI,
     &                  INTADR,IPRINT)
#include "implicit.h"
#include "dummy.h"
#include "maxaqn.h"
#include "maxmom.h"
#include "mxcent.h"
#include "maxorb.h"
      INTEGER B, X, C, D
      LOGICAL FULMAT, ANTI
      DIMENSION SHLINT(KCKTAB,60), SOINT(NELMNT,NOPTYP),
     &          INTADR(*)
#include "onecom.h"
#include "symmet.h"
#include "ibtfun.h"
      FULMAT = .TRUE.
      ICOMP = 0
      DO 100 X = 1, 3
         JCOMP = 0
         DO 200 B = 1, 3
         DO 200 C = B, 3
         DO 200 D = C, 3
            ICOMP = ICOMP + 1
            JCOMP = JCOMP + 1
            IREPB = IBTXOR(ISYMAX(B,1),ISYMAX(C,1))
            IREPC = IBTXOR(IREPB,ISYMAX(D,1))
            DO IREPO = 0, MAXREP
               ISCOOR = IPTCNT(3*(NCENTA - 1) + X,IREPO,1)
               IF (ISCOOR .GT. 0) THEN
                  ITYP   = INTADR(10*(ISCOOR - 1) + JCOMP)
                  IREPCL = IBTXOR(IREPO,IREPC)
                  FACSYM = HKAB*PT(IBTAND(IREPC,ISYMOP))
     &                         *PT(IBTAND(IREPO,ISYMOP))
                  IF (IREPO .EQ. 0) THEN
                     CALL SYM1S(SHLINT(1,ICOMP),SOINT(1,ITYP),ISYMOP,
     &                          MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                          FACSYM,LDIAG,FULMAT,DUM,IDUM,IPRINT)
                  ELSE
                     CALL SYM1N(SHLINT(1,ICOMP),SOINT(1,ITYP),IREPO,
     &                          ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,
     &                          KHKTB,FACSYM,LDIAG,FULMAT,ANTI,DUM,
     &                          IDUM,IPRINT)
                  END IF
               END IF
               ISCOOR = IPTCNT(3*(NCENTB - 1) + X,IREPO,1)
               IF (ISCOOR .GT. 0) THEN
                  ITYP   = INTADR(10*(ISCOOR - 1) + JCOMP)
                  IREPCL = IBTXOR(IREPO,IREPC)
                  FACSYM = HKAB*PT(IBTAND(IREPC,ISYMOP))
     &                         *PT(IBTAND(IREPO,ISYMOP))
                  IF (IREPO .EQ. 0) THEN
                     CALL SYM1S(SHLINT(1,ICOMP + 30),SOINT(1,ITYP),
     &                          ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,
     &                          KHKTB,FACSYM,LDIAG,FULMAT,DUM,IDUM,
     &                          IPRINT)
                  ELSE
                     CALL SYM1N(SHLINT(1,ICOMP + 30),SOINT(1,ITYP),
     &                          IREPO,ISYMOP,MULA,MULB,NHKTA,NHKTB,
     &                          KHKTA,KHKTB,FACSYM,LDIAG,FULMAT,ANTI,
     &                          DUM,IDUM,IPRINT)
                  END IF
               END IF
            END DO
 200     CONTINUE
 100  CONTINUE
      RETURN
      END
C  /* Deck hdotra */
      SUBROUTINE HDOTRA(SHLINT,SOINT,ISYMOP,NELMNT,NOPTYP,IORBA,IORBB,
     &                  NBAST,INTTYP,IPRINT)
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
      LOGICAL FULMAT, ANTI
      DIMENSION SHLINT(KCKTAB,NOPTYP), SOINT(NELMNT,NOPTYP)
#include "onecom.h"
#include "symmet.h"
#include "ibtfun.h"
      FULMAT = .TRUE.
      ANTI = .FALSE.
      DO 100 IREPO = 0, MAXREP
         DO 200 ICOOR = 1, 3
            IF (INTTYP .EQ. 14) THEN
               ISCOOR = IPTCNT(3*(NCENTA - 1) + ICOOR,IREPO,1)
               FAC = HKAB
            ELSE
C              ... INTTYP 44 and 98
               ISCOOR = IPTCNT(3*(NCENTB - 1) + ICOOR,IREPO,1)
               FAC = HKAB*PT(IBTAND(ISYMAX(ICOOR,1),ISYMOP))
     *                   *PT(IBTAND(IREPO, ISYMOP))
            END IF
            IF (ISCOOR .GT. 0) THEN
               CALL SYMSQR(SHLINT(1,ICOOR),SOINT(1,ISCOOR),IREPO,ISYMOP,
     &                     IORBA,IORBB,FAC,NBAST,IPRINT)
            END IF
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck symsqr */
      SUBROUTINE SYMSQR(AO,SO,IREPO,ISYMOP,IORBA,IORBB,FACTOR,NBAST,
     &                  IPRINT)
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
      DIMENSION AO(KHKTAB), SO(NBAST,NBAST)
#include "onecom.h"
#include "symmet.h"
#include "ibtfun.h"
      DO 100 IREPA = 0, MAXREP
         IREPB = IBTXOR(IREPO,IREPA)
         DO 200 NA = 1, KHKTA
         IF (IBTAND(MULA,IBTXOR(IREPA,ISYMAO(NHKTA,NA))).EQ.0) THEN
            NAT = KHKTB*(NA - 1)
            IA  = IPTSYM(IORBA + NA,IREPA)
            DO 300 NB = 1,KHKTB
            IF (IBTAND(MULB,IBTXOR(IREPB,ISYMAO(NHKTB,NB))).EQ.0) THEN
               IB  = IPTSYM(IORBB + NB,IREPB)
               FAC = PT(IBTAND(ISYMOP,IBTXOR(IREPB,ISYMAO(NHKTB,NB))))
               SO(IA,IB) = SO(IA,IB) + FAC*FACTOR*AO(NAT+NB)
            END IF
300         CONTINUE
         END IF
200      CONTINUE
100   CONTINUE
      RETURN
      END
C  /* Deck s1htra */
      SUBROUTINE S1HTRA(SHLINT,SOINT,ISYMOP,NELMNT,NOPTYP,ANTI,IPRINT)
C
C     tuh
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
      LOGICAL FULMAT, ANTI
      DIMENSION SHLINT(KCKTAB,3), SOINT(NELMNT,NOPTYP)
#include "onecom.h"
#include "symmet.h"
#include "ibtfun.h"
C
      FULMAT = .TRUE.
C
C     Loops over ireps and Cartesian coordinats
C
      DO 100 IREP = 0, MAXREP
      DO 100 ICOOR = 1, 3
C
C        Center A
C        ========
C
         ISCORA = IPTCNT(3*(NCENTA - 1) + ICOOR,IREP,1)
         IF (ISCORA .GT. 0) THEN
            FC = HKAB
            IF (IREP .EQ. 0) THEN
               CALL SYM1S(SHLINT(1,ICOOR),SOINT(1,ISCORA),
     &                    ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                    FC,LDIAG,FULMAT,DUM,IDUM,IPRINT)
            ELSE
               CALL SYM1N(SHLINT(1,ICOOR),SOINT(1,ISCORA),
     &                    IREP,ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                    FC,LDIAG,FULMAT,ANTI,DUM,IDUM,IPRINT)
            END IF
         END IF
C
C        Center B
C        ========
C
         ISCORB = IPTCNT(3*(NCENTB - 1) + ICOOR,IREP,1)
         IF (ISCORB .GT. 0) THEN
            FC = -HKAB*PT(IBTAND(ISYMAX(ICOOR,1),ISYMOP))
     &               *PT(IBTAND(IREP,ISYMOP))
            IF (IREP .EQ. 0) THEN
               CALL SYM1S(SHLINT(1,ICOOR),SOINT(1,ISCORB),
     &                    ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                    FC,LDIAG,FULMAT,DUM,IDUM,IPRINT)
            ELSE
               CALL SYM1N(SHLINT(1,ICOOR),SOINT(1,ISCORB),
     &                    IREP,ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                    FC,LDIAG,FULMAT,ANTI,DUM,IDUM,IPRINT)
            END IF
         END IF
C
  100 CONTINUE
C
      RETURN
      END
C  /* Deck amdtra */
      SUBROUTINE AMDTRA(SHLINT,SOINT,ISYMOP,NELMNT,NOPTYP,ANTI,
     &                  INTADR,IPRINT)
#include "implicit.h"
#include "dummy.h"
#include "maxaqn.h"
#include "maxmom.h"
#include "mxcent.h"
#include "maxorb.h"
      INTEGER B, X
      LOGICAL FULMAT, ANTI
      DIMENSION SHLINT(KCKTAB,18), SOINT(NELMNT,NOPTYP),
     &          INTADR(*)
#include "onecom.h"
#include "symmet.h"
#include "ibtfun.h"
      FULMAT = .TRUE.
      ICOMP = 0

      ! loop over cartesian displacements of the atoms
      DO 100 X = 1, 3
         JCOMP = 0

         ! loop over components of the dipole operator
         DO 200 B = 1, 3
            ICOMP = ICOMP + 1
            JCOMP = JCOMP + 1
            IREPC  = ISYMAX(B,2) ! irrep of angular moment component

            ! loop over possible irreps of symmetry coordinates
            DO IREP = 0, MAXREP
               ISCOOR = IPTCNT(3*(NCENTA - 1) + X,IREP,1)
               IF (ISCOOR .GT. 0) THEN
                  ITYP   = INTADR(3*(ISCOOR - 1) + JCOMP)
ckr                  write (2,*) 'center 1',iscoor,jcomp,ityp
                  IREPO = IBTXOR(IREP,IREPC) ! irrep of operator
ckr                  FACSYM = HKAB*PT(IBTAND(IREPO,ISYMOP))
C                  FACSYM = HKAB*PT(IBTAND(IREPC,ISYMOP))
C     &                         *PT(IBTAND(IREP,ISYMOP))
ckr                  FACSYM = HKAB
                  FACSYM = HKAB
                  IF (IREPO .EQ. 0) THEN
                     CALL SYM1S(SHLINT(1,ICOMP+9),SOINT(1,ITYP),ISYMOP,
     &                          MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                          FACSYM,LDIAG,FULMAT,DUM,IDUM,IPRINT)
                  ELSE
                     CALL SYM1N(SHLINT(1,ICOMP+9),SOINT(1,ITYP),IREPO,
     &                          ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,
     &                          KHKTB,FACSYM,LDIAG,FULMAT,ANTI,DUM,
     &                          IDUM,IPRINT)
                  END IF
               END IF
               ISCOOR = IPTCNT(3*(NCENTB - 1) + X,IREP,1)
               IF (ISCOOR .GT. 0) THEN
                  ITYP   = INTADR(3*(ISCOOR - 1) + JCOMP)
ckr                  write (2,*) 'center 2',iscoor,jcomp,ityp
                  IREPO = IBTXOR(IREP,IREPC) ! irrep of operator
ckr                  FACSYM = HKAB*PT(IBTAND(IREPO,ISYMOP))
C                  FACSYM = HKAB*PT(IBTAND(IREPC,ISYMOP))
C    &                          *PT(IBTAND(IREP,ISYMOP))
                   FACSYM = -HKAB*PT(IBTAND(ISYMAX(X,1),ISYMOP))
     &                          *PT(IBTAND(IREP,ISYMOP))
                  IF (IREPO .EQ. 0) THEN
                     CALL SYM1S(SHLINT(1,ICOMP),SOINT(1,ITYP),
     &                          ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,
     &                          KHKTB,FACSYM,LDIAG,FULMAT,DUM,IDUM,
     &                          IPRINT)
                  ELSE
                     CALL SYM1N(SHLINT(1,ICOMP),SOINT(1,ITYP),
     &                          IREPO,ISYMOP,MULA,MULB,NHKTA,NHKTB,
     &                          KHKTA,KHKTB,FACSYM,LDIAG,FULMAT,ANTI,
     &                          DUM,IDUM,IPRINT)
                  END IF
               END IF
            END DO
 200     CONTINUE
 100  CONTINUE
      RETURN
      END
CC  /* Deck amdtra */
C      SUBROUTINE AMDTRA(SHLINT,SOINT,ISYMOP,NELMNT,NOPTYP,IORBA,IORBB,
C     &                  NBAST,IPRINT)
CC
CC     K.Ruud, spring 2007
CC
C#include "implicit.h"
C#include "priunit.h"
C#include "dummy.h"
C#include "maxaqn.h"
C#include "mxcent.h"
C#include "maxorb.h"
C      LOGICAL FULMAT
C      DIMENSION SHLINT(KCKTAB,9), SOINT(NELMNT,NOPTYP)
C#include "onecom.h"
C#include "symmet.h"
C#include "ibtfun.h"
CC
C      FULMAT = .TRUE.
CC
CC     Loops over ireps and Cartesian coordinats
CC
C      DO IREP = 0, MAXREP
C         DO JBCOR = 1, 3
C            DO ICOOR = 1, 3
CC
CC        Center A
CC        ========
CC
C               IREP2 = IBTXOR(ISYMAX(JBCOR,2),IREP)
C               ICOMP = (ICOOR - 1)*3 + JBCOR
C               ISCORA = IPTCNT(3*(NCENTA - 1) + ICOOR,IREP,1)
Cckr               write (2,*) 'ncenta,icoor,iscora,irep',
Cckr     &              ncenta,icoor,iscora,jbcor
C               IF (ISCORA .GT. 0) THEN
C                  FC = HKAB
C                  ISCORT = 3*(ISCORA - 1) + JBCOR
Cckr                  WRITE (2,*) 'jbcor,icop,iscort',
Cckr     &                 jbcor,icomp,iscort,irep2,isymop
C                  CALL SYMSQR(SHLINT(1,ICOMP),SOINT(1,ISCORT),IREP2,
C     &                        ISYMOP,IORBA,IORBB,FC,NBAST,IPRINT)
C               END IF
CC
CC        Center B
CC        ========
CC
C               ISCORB = IPTCNT(3*(NCENTB - 1) + ICOOR,IREP,1)
CC               write (2,*) 'ncentb,icoor,iscorb,irep',
CC     &              ncentb,icoor,iscorb,irep
C               IF (ISCORB .GT. 0) THEN
CCkr Things to do
C                  FC = -HKAB*PT(IBTAND(ISYMAX(ICOOR,1),ISYMOP))
C     &                 *PT(IBTAND(IREP,ISYMOP))
C                  ISCORT = 3*(ISCORB - 1) + JBCOR
CC                  write (2,*) 'iscorb,jbcor,ncenta,jscort,icoor',
CC     &                 iscorb,jbcor,ncenta,iscort,icoor
C                  CALL SYMSQR(SHLINT(1,ICOMP),SOINT(1,ISCORT),IREP2,
C     &                        ISYMOP,IORBA,IORBB,FC,NBAST,IPRINT)
C               END IF
C            END DO
C         END DO
C      END DO
CC
C      RETURN
C      END
C  /* Deck h1dtra */
      SUBROUTINE G1HTRA(SHLINT,SOINT,FACINT,JSYMC,JCENTC,NCENTC,NATOMC,
     &                  ISYMOP,NELMNT,NOPTYP,ANTI,INTADR,IPRINT)
C
C     tuh
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
      LOGICAL FULMAT, ANTI
      DIMENSION SHLINT(KCKTAB,NOPTYP,3), SOINT(NELMNT,NOPTYP),
     &          FACINT(NATOMC), JSYMC(NATOMC), JCENTC(NATOMC),
     &          NCENTC(NATOMC),INTADR(*)
#include "onecom.h"
#include "symmet.h"
#include "ibtfun.h"
C
      FULMAT = .TRUE.
C
C     Loops over ireps and Cartesian coordinats
C
      DO 100 IREP = 0, MAXREP
      DO 100 ICOOR = 1, 3
C
C        Center A
C        ========
C
         ICOORA = 3*(ICENTA - 1) + ICOOR
         ISCORA = IPTCNT(3*(NCENTA - 1) + ICOOR,IREP,1)
         IF (ISCORA .GT. 0) THEN
            FC = HKAB
            IF (IREP .EQ. 0) THEN
               CALL SYM1S(SHLINT(1,ICOORA,1),SOINT(1,ISCORA),
     &                    ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                    FC,LDIAG,FULMAT,DUM,IDUM,IPRINT)
            ELSE
               CALL SYM1N(SHLINT(1,ICOORA,1),SOINT(1,ISCORA),
     &                    IREP,ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                    FC,LDIAG,FULMAT,ANTI,DUM,IDUM,IPRINT)
            END IF
         END IF
C
C        Center B
C        ========
C
         ICOORB = 3*(ICENTB - 1) + ICOOR
         ISCORB = IPTCNT(3*(NCENTB - 1) + ICOOR,IREP,1)
         IF (ISCORB .GT. 0) THEN
            FC = HKAB*PT(IBTAND(ISYMAX(ICOOR,1),ISYMOP))
     &               *PT(IBTAND(IREP,ISYMOP))
            IF (IREP .EQ. 0) THEN
               CALL SYM1S(SHLINT(1,ICOORB,2),SOINT(1,ISCORB),
     &                    ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                    FC,LDIAG,FULMAT,DUM,IDUM,IPRINT)
            ELSE
               CALL SYM1N(SHLINT(1,ICOORB,2),SOINT(1,ISCORB),
     &                    IREP,ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                    FC,LDIAG,FULMAT,ANTI,DUM,IDUM,IPRINT)
            END IF
         END IF
C
C        Center C
C        ========
C
         DO 200 IATOM = 1, NATOMC
            ICOORC = 3*(NCENTC(IATOM) - 1) + ICOOR
            ISCORC = IPTCNT(3*(JCENTC(IATOM) - 1) + ICOOR,IREP,1)
            IF (ISCORC .GT. 0) THEN
               FC = HKAB*PT(IBTAND(ISYMAX(ICOOR,1),JSYMC(IATOM)))
     &                  *PT(IBTAND(IREP,JSYMC(IATOM)))
               IF (IREP .EQ. 0) THEN
                  CALL SYM1S(SHLINT(1,ICOORC,3),SOINT(1,ISCORC),
     &                       ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                       FC,LDIAG,FULMAT,DUM,IDUM,IPRINT)
               ELSE
                  CALL SYM1N(SHLINT(1,ICOORC,3),SOINT(1,ISCORC),
     &                       IREP,ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,
     &                       KHKTB,FC,LDIAG,FULMAT,ANTI,DUM,IDUM,IPRINT)
               END IF
            END IF
  200    CONTINUE
C
  100 CONTINUE
C
      RETURN
      END
C  /* Deck sdint */
      SUBROUTINE SDINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT,
     &                 NOPTYP,NATOMC,INTTYP,SAAB,EXPP,CORCX,CORCY,CORCZ,
     &                 CORPX,CORPY,CORPZ)
C
C     tuh Feb 09 91 (use g_e factor 930726-hjaaj)
C
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "pi.h"
#include "gfac.h"
      PARAMETER (D2 = 2.0D0, D3 = 3.0D0)
      PARAMETER (DFAC11 = D2*GFAC*PI/D3, DFAC13 = D2*GFAC*PI)
      PARAMETER (DFACSD = GFAC / D2)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),AHGTF(*),
     &          CORCX(*), CORCY(*), CORCZ(*), SHLINT(KCKTAB,NOPTYP)
#include "onecom.h"
#include "lmns.h"
C
C        Correct Fermi contact contribution:
C        SDINT uses grad**2 (1/r) which gives -0.5 times
C        the Fermi contact contribution. This is removed
C        with DFAC11. IF SD+FC we add the Fermi contact
C        contribution to this correction, which gives DFAC13.
C
      IF (INTTYP .EQ. 11) THEN
         FACTOR = DFAC11
      ELSE
         FACTOR = DFAC13
      END IF
C
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
         MAXT = LVALA + LVALB
         MAXU = MVALA + MVALB
         MAXV = NVALA + NVALB
         INT = INT + 1
         DO 150 IATOM = 1,NATOMC
            I0 = 6*(IATOM - 1)
            XC = CORCX(IATOM)
            YC = CORCY(IATOM)
            ZC = CORCZ(IATOM)
            DIS = (XC-CORPX)**2 + (YC-CORPY)**2 + (ZC-CORPZ)**2
            DIR = SAAB*FACTOR*EXP(-EXPP*DIS)
            IF (LVALA.GT.0) DIR = DIR*((XC - CORAX)**LVALA)
            IF (MVALA.GT.0) DIR = DIR*((YC - CORAY)**MVALA)
            IF (NVALA.GT.0) DIR = DIR*((ZC - CORAZ)**NVALA)
            IF (LVALB.GT.0) DIR = DIR*((XC - CORBX)**LVALB)
            IF (MVALB.GT.0) DIR = DIR*((YC - CORBY)**MVALB)
            IF (NVALB.GT.0) DIR = DIR*((ZC - CORBZ)**NVALB)
            SHLINT(INT,I0+1) = SHLINT(INT,I0+1) + DIR
            SHLINT(INT,I0+3) = SHLINT(INT,I0+3) + DIR
            SHLINT(INT,I0+6) = SHLINT(INT,I0+6) + DIR
 150     CONTINUE
         IADRAV = 1
         DO 200 IV = 0,MAXV
            EV = DFACSD*ODC(NVALA,NVALB,IV,0,0,3)
            IADRAU = IADRAV
            DO 300 IU = 0,MAXU
               EE = ODC(MVALA,MVALB,IU,0,0,2)*EV
               DO 400 IT = 0,MAXT
                  EEE = ODC(LVALA,LVALB,IT,0,0,1)*EE
                  IADR00 = IADRAU + IT
                  IADR0T = IADR00 + 1
                  IADR0U = IADR00 + ISTEPU
                  IADR0V = IADR00 + ISTEPV
                  IADRTT = IADR0T + 1
                  IADRTU = IADR0T + ISTEPU
                  IADRTV = IADR0T + ISTEPV
                  IADRUU = IADR0U + ISTEPU
                  IADRUV = IADR0U + ISTEPV
                  IADRVV = IADR0V + ISTEPV
                  IADD = - NAHGTF
C
C                 ***** LOOP OVER NUCLEI *****
C
                  DO 500 IATOM = 1,NATOMC
                     I0 = 6*(IATOM - 1)
                     IADD = IADD + NAHGTF
                     SHLINT(INT,I0+1) = SHLINT(INT,I0+1)
     &                                + EEE*AHGTF(IADRTT + IADD)
                     SHLINT(INT,I0+2) = SHLINT(INT,I0+2)
     &                                + EEE*AHGTF(IADRTU + IADD)
                     SHLINT(INT,I0+3) = SHLINT(INT,I0+3)
     &                                + EEE*AHGTF(IADRUU + IADD)
                     SHLINT(INT,I0+4) = SHLINT(INT,I0+4)
     &                                + EEE*AHGTF(IADRTV + IADD)
                     SHLINT(INT,I0+5) = SHLINT(INT,I0+5)
     &                                + EEE*AHGTF(IADRUV + IADD)
                     SHLINT(INT,I0+6) = SHLINT(INT,I0+6)
     &                                + EEE*AHGTF(IADRVV + IADD)
  500             CONTINUE
  400          CONTINUE
               IADRAU = IADRAU + ISTEPU
  300       CONTINUE
            IADRAV = IADRAV + ISTEPV
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
      SUBROUTINE SDKINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT,
     &                 NOPTYP,NATOMC,INTTYP,SAAB,EXPP,EXPB,CORCX,CORCY,
     &                 CORCZ,CORPX,CORPY,CORPZ)
C
C     K.Ruud, sep.2001
C
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "pi.h"
#include "gfac.h"
      PARAMETER (D1 = 1.0D0, D2 = 2.0D0, D4 = 4.0D0, D3 = 3.0D0)
      PARAMETER (DFAC11 = GFAC*PI/D3)
      PARAMETER (DFACSD = GFAC / D4)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),AHGTF(*),
     &          CORCX(*), CORCY(*), CORCZ(*), SHLINT(KCKTAB,NOPTYP)
#include "onecom.h"
#include "lmns.h"

      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
         MAXT = LVALA + LVALB
         MAXU = MVALA + MVALB
         MAXV = NVALA + NVALB
         INT = INT + 1
C
C        Correct Fermi contact contribution:
C        SDINT uses grad**2 (1/r) which gives -0.5 times
C        the Fermi contact contribution. This is removed
C        with DFAC11. IF SD+FC we add the Fermi contact
C        contribution to this correction, which gives DFAC13.
C
         DO 150 IATOM = 1,NATOMC
            I0 = 6*(IATOM - 1)
            XC = CORCX(IATOM)
            YC = CORCY(IATOM)
            ZC = CORCZ(IATOM)
            DIS = (XC-CORPX)**2 + (YC-CORPY)**2 + (ZC-CORPZ)**2
            FINT = SAAB*DFAC11*EXP(-EXPP*DIS)
C
C     Contributions from center A
C
            IF (LVALA .GT. 0) FINT = FINT*((XC - CORAX)**LVALA)
            IF (MVALA .GT. 0) FINT = FINT*((YC - CORAY)**MVALA)
            IF (NVALA .GT. 0) FINT = FINT*((ZC - CORAZ)**NVALA)
C
C     Undifferentiated contribution from center B
C
            FINTYZ = D1
            FINTXZ = D1
            FINTXY = D1
            IF (LVALB .GT. 0) THEN
               FINTXY = FINTXY*((XC - CORBX)**LVALB)
               FINTXZ = FINTXZ*((XC - CORBX)**LVALB)
            END IF
            IF (MVALB .GT. 0) THEN
               FINTXY = FINTXY*((YC - CORBY)**MVALB)
               FINTYZ = FINTYZ*((YC - CORBY)**MVALB)
            END IF
            IF (NVALB .GT. 0) THEN
               FINTYZ = FINTYZ*((ZC - CORBZ)**NVALB)
               FINTXZ = FINTXZ*((ZC - CORBZ)**NVALB)
            END IF
C
C     Differentiated contribution from center B
C
            IF (LVALB .GT. 0) THEN
               XFINT = - D2*EXPB*((XC-CORBX)**LVALB)
               XFINT = (2*LVALB + 1) * XFINT
            ELSE
               XFINT = - D2*EXPB
            END IF
            XFINT = XFINT + D4*EXPB*EXPB*((XC - CORBX)**(LVALB + 2))
            IF (LVALB .GT. 2) THEN
               XFINT = XFINT +
     &              (LVALB*(LVALB-1.0D0)) * ((XC-CORBX)**(LVALB-2))
            ELSE IF (LVALB .EQ. 2) THEN
               XFINT = XFINT + D2
            END IF
C
            IF (MVALB .GT. 0) THEN
               YFINT = - D2*EXPB*((YC-CORBY)**MVALB)
               YFINT = (2*MVALB + 1) * YFINT
            ELSE
               YFINT = - D2*EXPB
            END IF
            YFINT = YFINT + D4*EXPB*EXPB*((YC - CORBY)**(MVALB + 2))
            IF (MVALB .GT. 2) THEN
               YFINT = YFINT +
     &              (MVALB*(MVALB-1.0D0)) * ((YC-CORBY)**(MVALB-2))
            ELSE IF (MVALB .EQ. 2) THEN
               YFINT = YFINT + D2
            END IF
C
            IF (NVALB .GT. 0) THEN
               ZFINT = - D2*EXPB*((ZC-CORBZ)**NVALB)
               ZFINT = (2*NVALB + 1) * ZFINT
            ELSE
               ZFINT = - D2*EXPB
            END IF
            ZFINT = ZFINT + D4*EXPB*EXPB*((ZC - CORBZ)**(NVALB + 2))
            IF (NVALB .GT. 2) THEN
               ZFINT = ZFINT +
     &              (NVALB*(NVALB-1.0D0)) * ((ZC-CORBZ)**(NVALB-2))
            ELSE IF (NVALB .EQ. 2) THEN
               ZFINT = ZFINT + D2
            END IF
            FRM = FINT*(XFINT*FINTYZ + YFINT*FINTXZ + ZFINT*FINTXY)
            SHLINT(INT,I0+1) = SHLINT(INT,I0+1) + FRM
            SHLINT(INT,I0+3) = SHLINT(INT,I0+3) + FRM
            SHLINT(INT,I0+6) = SHLINT(INT,I0+6) + FRM
 150     CONTINUE
         IADRAV = 1
         DO 200 IV = 0,MAXV + 2
            EV = DFACSD*ODC(NVALA,NVALB,IV,0,0,3)
            FV = DFACSD*ODC(NVALA,NVALB,IV,2,0,3)
            IADRAU = IADRAV
            DO 300 IU = 0,MAXU + 2
               EE = ODC(MVALA,MVALB,IU,0,0,2)*EV
               FE = ODC(MVALA,MVALB,IU,2,0,2)*EV
               EF = ODC(MVALA,MVALB,IU,0,0,2)*FV
               DO 400 IT = 0,MAXT + 2
                  FEE = ODC(LVALA,LVALB,IT,2,0,1)*EE
                  EFE = ODC(LVALA,LVALB,IT,0,0,1)*FE
                  EEF = ODC(LVALA,LVALB,IT,0,0,1)*EF
                  TOTEXP = FEE + EFE + EEF
                  IADR00 = IADRAU + IT
                  IADR0T = IADR00 + 1
                  IADR0U = IADR00 + ISTEPU
                  IADR0V = IADR00 + ISTEPV
                  IADRTT = IADR0T + 1
                  IADRTU = IADR0T + ISTEPU
                  IADRTV = IADR0T + ISTEPV
                  IADRUU = IADR0U + ISTEPU
                  IADRUV = IADR0U + ISTEPV
                  IADRVV = IADR0V + ISTEPV
                  IADD = - NAHGTF
C
C                 ***** LOOP OVER NUCLEI *****
C
                  DO 500 IATOM = 1,NATOMC
                     I0 = 6*(IATOM - 1)
                     IADD = IADD + NAHGTF
                     SHLINT(INT,I0+1) = SHLINT(INT,I0+1)
     &                     + TOTEXP*AHGTF(IADRTT + IADD)
                     SHLINT(INT,I0+2) = SHLINT(INT,I0+2)
     &                     + TOTEXP*AHGTF(IADRTU + IADD)
                     SHLINT(INT,I0+3) = SHLINT(INT,I0+3)
     &                     + TOTEXP*AHGTF(IADRUU + IADD)
                     SHLINT(INT,I0+4) = SHLINT(INT,I0+4)
     &                     + TOTEXP*AHGTF(IADRTV + IADD)
                     SHLINT(INT,I0+5) = SHLINT(INT,I0+5)
     &                     + TOTEXP*AHGTF(IADRUV + IADD)
                     SHLINT(INT,I0+6) = SHLINT(INT,I0+6)
     &                     + TOTEXP*AHGTF(IADRVV + IADD)
  500             CONTINUE
  400          CONTINUE
               IADRAU = IADRAU + ISTEPU
  300       CONTINUE
            IADRAV = IADRAV + ISTEPV
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck sdtra */
      SUBROUTINE SDTRA(SHLINT,SOINT,FACINT,JSYMC,JCENTC,NATOMC,ISYMOP,
     &                 NELMNT,NOPTYP,ANTI,IORBA,IORBB,NBAST,INTTYP,
     &                 INTADR,IPRINT)
#include "implicit.h"
#include "dummy.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
      LOGICAL FULMAT, ANTI
      DIMENSION SHLINT(KCKTAB,NOPTYP), SOINT(NELMNT,NOPTYP),
     &          INTADR(*), FACINT(NATOMC), JSYMC(NATOMC),
     &          JCENTC(NATOMC)
#include "onecom.h"
#include "symmet.h"
#include "ibtfun.h"
      FULMAT = .TRUE.
      DO 100 IREPC = 0, MAXREP
         DO 200 IATOMC = 1, NATOMC
            ICENTC =   JCENTC(IATOMC)
            ISYMC  =   JSYMC (IATOMC)
            FACTOR = - FACINT(IATOMC)
            DO 300 ICOOR1 = 1, 3
               ISCOR1 = IPTCNT(3*(ICENTC - 1) + ICOOR1,IREPC,2)
               IF (ISCOR1 .GT. 0) THEN
                  DO 400 ICOOR2 = 1, 3
                     ISCOOR = 3*(ISCOR1 - 1) + ICOOR2
                     IREPCD = IBTXOR(IREPC,ISYMAX(ICOOR2,2))
                     ISYMCR = IBTXOR(ISYMAX(ICOOR1,1),ISYMAX(ICOOR2,1))
                     FACSYM = FACTOR*PT(IBTAND(ISYMCR,ISYMC))
     &                              *PT(IBTAND(IREPCD,ISYMC))
                     IADR = INTADR(ISCOOR)
                     MXCR = MAX(ICOOR1,ICOOR2)
                     MNCR = MIN(ICOOR1,ICOOR2)
                     ITYP = 6*(IATOMC - 1) + MXCR*(MXCR - 1)/2 + MNCR
                     IF (INTTYP .EQ. 94) THEN
                        CALL SYMSQR(SHLINT(1,ITYP),SOINT(1,IADR),IREPCD,
     &                           ISYMOP,IORBA,IORBB,FACSYM,NBAST,IPRINT)
                     ELSE
                        IF (IREPCD .EQ. 0) THEN
                         CALL SYM1S(SHLINT(1,ITYP),SOINT(1,IADR),
     &                              ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,
     &                              KHKTB,FACSYM,LDIAG,FULMAT,DUM,IDUM,
     &                              IPRINT)
                        ELSE
                         CALL SYM1N(SHLINT(1,ITYP),SOINT(1,IADR),IREPCD,
     &                              ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,
     &                              KHKTB,FACSYM,LDIAG,FULMAT,ANTI,DUM,
     &                              IDUM,IPRINT)
                        END IF
                     END IF
  400             CONTINUE
               END IF
  300       CONTINUE
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck npetra */
      SUBROUTINE NPETRA(SHLINT,SOINT,FACINT,JSYMC,JCENTC,NATOMC,ISYMOP,
     &                  NELMNT,NOPTYP,ANTI,INTADR,IPRINT)
C
C     K.Ruud, July 92
C     Symmetry finally introduced Aug-01, kr
C
#include "implicit.h"
#include "dummy.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "nuclei.h"
#include "symmet.h"
      LOGICAL FULMAT, ANTI
      DIMENSION SHLINT(KCKTAB,NOPTYP), SOINT(NELMNT,NOPTYP),
     &          FACINT(NATOMC), JSYMC(NATOMC), JCENTC(NATOMC)
#include "onecom.h"
#include "ibtfun.h"
      FULMAT = .TRUE.
      ITYP = 0
      DO IREPC = 0, MAXREP
         DO IATOMC = 1, NATOMC
            ICENTC = JCENTC(IATOMC)
            ISYMC  = JSYMC (IATOMC)
            FACTOR = FACINT(IATOMC)
            IF (IBTAND(ISTBNU(ICENTC),IREPC) .EQ. 0) THEN
               FACSYM = FACTOR*PT(IBTAND(IREPC,ISYMC))
               IADR = IPTNUC(ICENTC,IREPC)
               IF (IREPC .EQ. 0) THEN
                  CALL SYM1S(SHLINT(1,IATOMC),SOINT(1,IADR),ISYMOP,
     &                       MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                       FACSYM,LDIAG,FULMAT,DUM,IDUM,IPRINT)
               ELSE
                  CALL SYM1N(SHLINT(1,IATOMC),SOINT(1,IADR),IREPC,
     &                       ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                       FACSYM,LDIAG,FULMAT,ANTI,DUM,IDUM,IPRINT)
               END IF
            END IF
         END DO
      END DO
      RETURN
      END
C  /* Deck testra */
      SUBROUTINE TESTRA(SHLINT,SOINT,FACINT,JSYMC,JCENTC,NATOMC,ISYMOP,
     &                  NELMNT,NOPTYP,ANTI,EXP1EL,DMAT,MULTA,MULTB,
     &                  IPRINT)
C
C     K.Ruud, June 2001
C
#include "implicit.h"
#include "dummy.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "symmet.h"
      LOGICAL FULMAT, ANTI, EXP1EL
      DIMENSION SHLINT(KCKTAB,NOPTYP), SOINT(NELMNT,NOPTYP),
     &          FACINT(NATOMC), JSYMC(NATOMC), JCENTC(NATOMC)
      DIMENSION DMAT(*)
#include "krcom.h"
#include "onecom.h"
#include "ibtfun.h"
      FULMAT = .TRUE.
      ITYP = 0
      DO IREPC = 0, MAXREP
         DO IATOMC = 1, NATOMC
            ICENTC = JCENTC(IATOMC)
            ISYMC  = JSYMC (IATOMC)
            FACTOR = FACINT(IATOMC)
            IF (IBTAND(MAXREP+1,IREPC) .EQ. 0) THEN
               FACSYM = FACTOR*PT(IBTAND(IREPC,ISYMC))
               IF (EXP1EL) THEN
                  IF (ISYMD .EQ. 0) THEN
                     CALL SYM1EV(SHLINT(1,IATOMC),DMAT,SOINT(1,IREPC+1),
     &                           -FACSYM,MULTA,MULTB,IPRINT)
                  ELSE
                     IF (ISYMD .EQ. IREPC)
     &                    CALL SYM1NC(SHLINT(1,IATOMC),DMAT,
     &                                SOINT(1,IREPC + 1),IREPC,ISYMOP,
     &                                MULA,MULB,NHKTA,NHKTB,KHKTA,
     &                                KHKTB,FACSYM,LDIAG,ANTI,IPRINT)
                  END IF
               ELSE IF (IREPC .EQ. 0) THEN
                  CALL SYM1S(SHLINT(1,IATOMC),SOINT(1,IREPC + 1),ISYMOP,
     &                       MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                       FACSYM,LDIAG,FULMAT,DUM,IDUM,IPRINT)
               ELSE
                  CALL SYM1N(SHLINT(1,IATOMC),SOINT(1,IREPC + 1),IREPC,
     &                       ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                       FACSYM,LDIAG,FULMAT,ANTI,DUM,IDUM,IPRINT)
               END IF
            END IF
         END DO
      END DO
      RETURN
      END
C  /* Deck eftra */
      SUBROUTINE EFTRA(SHLINT,SOINT,FACINT,JSYMC,JCENTC,NATOMC,
     &                  ISYMOP,NELMNT,NOPTYP,ANTI,INTADR,IPRINT,
     &                  MULTA,MULTB,EXPVAL,DMAT,EXP1EL)
C
C  KR, March 92
C     Symmetry properly treated october-95 (sic!)
C
#include "implicit.h"
#include "dummy.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
      LOGICAL FULMAT, ANTI, EXP1EL
      DIMENSION SHLINT(KCKTAB,NOPTYP), SOINT(NELMNT,NOPTYP),
     &          FACINT(NATOMC), JSYMC(NATOMC), JCENTC(NATOMC),
     &          INTADR(*)
      DIMENSION EXPVAL(*), DMAT(*)
#include "onecom.h"
#include "nuclei.h"
#include "symmet.h"
#include "ibtfun.h"
      FULMAT = .TRUE.
      DO 100 IREPC = 0, MAXREP
         DO 200 IATOMC = 1, NATOMC
            ICENTC = JCENTC(IATOMC)
            ISYMC  = JSYMC (IATOMC)
            FACTOR = FACINT(IATOMC)
            DO 300 ICOOR = 1, 3
               ISCOOR = IPTCNT(3*(ICENTC - 1) + ICOOR,IREPC,1)
               IF (ISCOOR .GT. 0) THEN
                  ISYMCR = ISYMAX(ICOOR,1)
                  FACSYM = - FACTOR*PT(IBTAND(ISYMCR,ISYMC))
     &                             *PT(IBTAND(IREPC ,ISYMC))
                  ITYP = 3*(IATOMC - 1) +ICOOR
                  IF (EXP1EL .AND. (IREPC .EQ. 0)) THEN
                     CALL SYM1EV(SHLINT(1,ITYP),DMAT,EXPVAL(ISCOOR),
     &                           FACSYM,MULTA,MULTB,IPRINT)
                  ELSE IF (IREPC .EQ. 0) THEN
                     CALL SYM1S(SHLINT(1,ITYP),SOINT(1,ISCOOR),ISYMOP,
     &                          MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                          FACSYM,LDIAG,FULMAT,DUM,IDUM,IPRINT)
                  ELSE
                     CALL SYM1N(SHLINT(1,ITYP),SOINT(1,ISCOOR),IREPC,
     &                          ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,
     &                          KHKTB,FACSYM,LDIAG,FULMAT,ANTI,DUM,
     &                          IDUM,IPRINT)
                  END IF
               END IF
 300        CONTINUE
 200     CONTINUE
 100  CONTINUE
      RETURN
      END
C  /* Deck efgtra */
      SUBROUTINE EFGTRA(SHLINT,SOINT,FACINT,JSYMC,JCENTC,NATOMC,ISYMOP,
     &                  NELMNT,NOPTYP,ANTI,IFGTBL,EXP1EL,DMAT,
     &                  MULTA,MULTB,IPRINT)
#include "implicit.h"
#include "dummy.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "nuclei.h"
#include "symmet.h"
      LOGICAL FULMAT, ANTI, EXP1EL
      DIMENSION SHLINT(KCKTAB,NOPTYP), SOINT(NELMNT,NOPTYP),
     &          FACINT(NATOMC), JSYMC(NATOMC), JCENTC(NATOMC),
     &          IFGTBL(NUCIND,6,0:MAXREP)
      DIMENSION DMAT(*)
#include "onecom.h"
#include "ibtfun.h"
C
      FULMAT = .TRUE.
      ITYP = 0
      DO 100 IATOMC = 1, NATOMC
         ICENTC = JCENTC(IATOMC)
         ISYMC  = JSYMC(IATOMC)
         FACTOR = FACINT(IATOMC)
         IJ = 0
         DO 200 ICOOR1 = 1, 3
         DO 200 ICOOR2 = ICOOR1, 3
            IJ = IJ + 1
            ITYP = ITYP + 1
            ISYMXY = IBTXOR(ISYMAX(ICOOR1,1),ISYMAX(ICOOR2,1))
            DO 300 IREPC = 0, MAXREP
              IF (IBTAND(ISTBNU(ICENTC),IBTXOR(IREPC,ISYMXY)).EQ.0) THEN
                  FACSYM = FACTOR*PT(IBTAND(ISYMXY,ISYMC))
     &                           *PT(IBTAND(IREPC ,ISYMC))
                  IADR = IFGTBL(ICENTC,IJ,IREPC)
                  IF (EXP1EL) THEN
                     IF (IREPC .EQ. 0) CALL SYM1EV(SHLINT(1,ITYP),DMAT,
     &                                      SOINT(1,IADR),FACSYM,
     &                                      MULTA,MULTB,IPRINT)
                  ELSE IF (IREPC .EQ. 0) THEN
                     CALL SYM1S(SHLINT(1,ITYP),SOINT(1,IADR),ISYMOP,
     &                          MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                          FACSYM,LDIAG,FULMAT,DUM,IDUM,IPRINT)
                  ELSE
                     CALL SYM1N(SHLINT(1,ITYP),SOINT(1,IADR),IREPC,
     &                          ISYMOP,MULA,MULB,NHKTA,NHKTB,KHKTA,
     &                          KHKTB,FACSYM,LDIAG,FULMAT,ANTI,DUM,
     &                          IDUM,IPRINT)
                  END IF
               END IF
 300        CONTINUE
 200     CONTINUE
 100   CONTINUE
       RETURN
       END
C  /* Deck dsotra */
      SUBROUTINE DSOTRA(SHLINT,SOINT,ISYMOP,NELMNT,NOPTYP,ANTI,DOATOM,
     &                  KAB,TRIANG,NATOM,INTADR,NSHINT,INTTYP,EXP1EL,
     &                  DMAT,MULTA,MULTB,IPRINT)
#include "implicit.h"
#include "dummy.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "nuclei.h"
      LOGICAL FULMAT, ANTI, DOATOM(*), SAMECD, TRIANG, TRICR, EXP1EL
      DIMENSION SHLINT(KCKTAB,NSHINT), SOINT(NELMNT,NOPTYP),
     &          INTADR(*), DMAT(*)
#include "onecom.h"
#include "symmet.h"
#include "ibtfun.h"
      FULMAT = .TRUE.
C
      DO 100 IREPO = 0, MAXREP
        ISTAO = 0
C
C       Atom D
C
        DO 110 IATOMD = 1, NUCIND
        IF (DOATOM(IATOMD)) THEN
          ISTABD = ISTBNU(IATOMD)
          DO 120 ISYMD = 0, MAXOPR
          IF (IBTAND(ISYMD,IBTOR(ISTABD,KAB)) .EQ. 0) THEN
            KABD = IBTAND(ISTABD,KAB)
C
C           Atom C
C
C           (First determine JATOMC)
C
            JATOMC = 0
            MXATMC = NUCIND
            IF (TRIANG .AND. INTTYP .EQ. 12) MXATMC = IATOMD
            DO 130 IATOMC = 1, MXATMC
            IF (DOATOM(IATOMC)) THEN
               SAMECD = IATOMC .EQ. IATOMD
               MABCD  = IBTOR(ISTBNU(IATOMC),KABD)
               DO 135 ISYMC = 0, MAXOPR
                  IF (SAMECD .AND. ISYMC.EQ.ISYMD) GO TO 135
                  IF (IBTAND(ISYMC,MABCD) .EQ. 0) JATOMC = JATOMC + 1
  135          CONTINUE
            END IF
  130       CONTINUE
            IF (JATOMC .GT. 0) THEN
              ITYP0 = ISTAO
              MXATMC = NUCIND
              IF (TRIANG .AND. INTTYP .EQ. 12) MXATMC = IATOMD
              DO 200 IATOMC = 1, MXATMC
              IF (DOATOM(IATOMC)) THEN
                SAMECD = IATOMC .EQ. IATOMD
                TRICR = (TRIANG .AND. INTTYP .EQ. 12) .AND. SAMECD
                DO 210 ISYMC = 0, MAXOPR
                IF (SAMECD .AND. ISYMC.EQ.ISYMD) GO TO 210
                IF (IBTAND(ISYMC,IBTOR(ISTBNU(IATOMC),KABD)).EQ.0) THEN
C
C                 Cartesian directions
C
                  DO 300 ICOORD = 1, 3
                    ISCORD = IPTCNT(3*(IATOMD - 1) + ICOORD,IREPO,2)
                    IF (ISCORD .GT. 0) THEN
                      FACD = - PT(IBTAND(ISYMAX(ICOORD,2),ISYMD))
     &                        *PT(IBTAND(IREPO ,ISYMD))
                      MXCRC = 3
                      IF (TRICR) MXCRC = ICOORD
                      DO 310 ICOORC = 1, MXCRC
                        ISCORC = IPTCNT(3*(IATOMC - 1) + ICOORC,IREPO,2)
                        IF (ISCORC .GT. 0) THEN
                          IF (INTTYP .EQ. 12) THEN
                             IF (TRIANG) THEN
                                MXCOR = MAX(ISCORC,ISCORD)
                                MNCOR = MIN(ISCORC,ISCORD)
                                ISCOOR = MXCOR*(MXCOR - 1)/2 + MNCOR
                             ELSE
                                ISCOOR = 3*NUCDEP*(ISCORD - 1) + ISCORC
                             END IF
                          ELSE
                             ISCOOR = 3*(ISCORD - 1) + ICOORC
                          END IF
                          IADR = INTADR(ISCOOR)
                          FACSYM=FACD*PT(IBTAND(ISYMAX(ICOORC,2),ISYMC))
     &                               *PT(IBTAND(IREPO ,ISYMC))
                          ITYP = ITYP0 + 3*(ICOORC - 1) + ICOORD
                          IF (INTTYP .EQ. 58) CALL DSCAL(KCKTAB,
     &                         CHARGE(IATOMC),SHLINT(1,ITYP),1)
                          IF (EXP1EL) THEN
                             CALL SYM1EV(SHLINT(1,ITYP),DMAT,
     &                                   SOINT(1,IADR),FACSYM,MULTA,
     &                                   MULTB,IPRINT)
                          ELSE
                             CALL SYM1S(SHLINT(1,ITYP),SOINT(1,IADR),
     &                                  ISYMOP,MULA,MULB,NHKTA,NHKTB,
     &                                  KHKTA,KHKTB,FACSYM,LDIAG,FULMAT,
     &                                  DUM,IDUM,IPRINT)
                          END IF
                        END IF
  310                 CONTINUE
                    END IF
  300             CONTINUE
                  ITYP0 = ITYP0 + 9
                END IF
  210           CONTINUE
              END IF
  200         CONTINUE
              ISTAO = ISTAO + 9*JATOMC
            END IF
          END IF
  120     CONTINUE
        END IF
  110   CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck dsoodc */
      SUBROUTINE DSOODC(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,EXPPI,DIFPDX,
     &                  DIFPDY,DIFPDZ,IPRINT,INTTYP)
C
C     TUH Feb 19 1991
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D0 = 0.D0, D2 = 2.D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)

      EXPPIH = EXPPI/D2
      DO 100 IA = 0, JMAXA
         DO 200 IB = 0, JMAXB
            DO 300 I = 0, IA + IB + 1
               FX = D0
               FY = D0
               FZ = D0
               X_IP1 = I + 1
               IF (I .LE. IA + IB) THEN
                  FX = FX + DIFPDX*ODC(IA,IB,I,0,0,1)
                  FY = FY + DIFPDY*ODC(IA,IB,I,0,0,2)
                  FZ = FZ + DIFPDZ*ODC(IA,IB,I,0,0,3)
               END IF
               IF (I .GT. 0) THEN
                  FX = FX + EXPPIH*ODC(IA,IB,I - 1,0,0,1)
                  FY = FY + EXPPIH*ODC(IA,IB,I - 1,0,0,2)
                  FZ = FZ + EXPPIH*ODC(IA,IB,I - 1,0,0,3)
               END IF
               IF (I .LT. IA + IB) THEN
                  FX = FX + X_IP1*ODC(IA,IB,I + 1,0,0,1)
                  FY = FY + X_IP1*ODC(IA,IB,I + 1,0,0,2)
                  FZ = FZ + X_IP1*ODC(IA,IB,I + 1,0,0,3)
               END IF
               IF (INTTYP .EQ. 96) THEN
                  ODC(IA,IB,I,0,1,1) = FX
                  ODC(IA,IB,I,0,1,2) = FY
                  ODC(IA,IB,I,0,1,3) = FZ
                  FX = D0
                  FY = D0
                  FZ = D0
                  IF (I .LE. IA + IB) THEN
                     FX = FX + DIFPDX*ODC(IA,IB,I,2,0,1)
                     FY = FY + DIFPDY*ODC(IA,IB,I,2,0,2)
                     FZ = FZ + DIFPDZ*ODC(IA,IB,I,2,0,3)
                  END IF
                  IF (I .GT. 0) THEN
                     FX = FX + EXPPIH*ODC(IA,IB,I - 1,2,0,1)
                     FY = FY + EXPPIH*ODC(IA,IB,I - 1,2,0,2)
                     FZ = FZ + EXPPIH*ODC(IA,IB,I - 1,2,0,3)
                  END IF
                  IF (I .LT. IA + IB) THEN
                     FX = FX + X_IP1*ODC(IA,IB,I + 1,2,0,1)
                     FY = FY + X_IP1*ODC(IA,IB,I + 1,2,0,2)
                     FZ = FZ + X_IP1*ODC(IA,IB,I + 1,2,0,3)
                  END IF
                  ODC(IA,IB,I,2,1,1) = FX
                  ODC(IA,IB,I,2,1,2) = FY
                  ODC(IA,IB,I,2,1,3) = FZ
               ELSE
                  ODC(IA,IB,I,1,0,1) = FX
                  ODC(IA,IB,I,1,0,2) = FY
                  ODC(IA,IB,I,1,0,3) = FZ
               END IF
  300       CONTINUE
  200    CONTINUE
  100 CONTINUE
      IF (IPRINT .GE. 20) THEN
         IF (INTTYP .EQ. 12 .OR. INTTYP .EQ. 58) THEN
            JMAXAB = (JMAXA + 1)*(JMAXB + 1)
            ISTP   = JMAXA + JMAXB + 2
            CALL TITLER('Output from DSOODC','*',103)
            CALL AROUND('x component')
            CALL OUTPUT(ODC(0,0,0,1,0,1),1,JMAXAB,1,ISTP,JMAXAB,ISTP,1,
     &                  LUPRI)
            CALL AROUND('y component')
            CALL OUTPUT(ODC(0,0,0,1,0,2),1,JMAXAB,1,ISTP,JMAXAB,ISTP,1,
     &                  LUPRI)
            CALL AROUND('z component')
            CALL OUTPUT(ODC(0,0,0,1,0,3),1,JMAXAB,1,ISTP,JMAXAB,ISTP,1,
     &                  LUPRI)
         ELSE IF (INTTYP .EQ. 96) THEN
            JMAXAB = (JMAXA + 1)*(JMAXB + 1)
            ISTP   = JMAXA + JMAXB + 2
            CALL TITLER('Output from DSOODC','*',103)
            CALL AROUND('x component')
            CALL OUTPUT(ODC(0,0,0,0,1,1),1,JMAXAB,1,ISTP,JMAXAB,ISTP,1,
     &                  LUPRI)
            CALL AROUND('y component')
            CALL OUTPUT(ODC(0,0,0,0,1,2),1,JMAXAB,1,ISTP,JMAXAB,ISTP,1,
     &                  LUPRI)
            CALL AROUND('z component')
            CALL OUTPUT(ODC(0,0,0,0,1,3),1,JMAXAB,1,ISTP,JMAXAB,ISTP,1,
     &                  LUPRI)
         END IF
      END IF
      RETURN
      END
C  /* Deck dsoint */
      SUBROUTINE DSOINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT,
     &                  NOPTYP,NATOMC)
C
C     tuh & ov Feb 09 91
C
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),AHGTF(*),
     &          SHLINT(KCKTAB,*)
#include "onecom.h"
#include "nuclei.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
         MAXT = LVALA + LVALB
         MAXU = MVALA + MVALB
         MAXV = NVALA + NVALB
         INT = INT + 1
         IADRAV = 1
         DO 200 IV = 0, MAXV + 1
            IADRAU = IADRAV
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            FV = ODC(NVALA,NVALB,IV,1,0,3)
            IUMAX = MAXU + 1
            IF (IV .GT. MAXV) IUMAX = MAXU
            DO 300 IU = 0, IUMAX
               EU = ODC(MVALA,MVALB,IU,0,0,2)
               FU = ODC(MVALA,MVALB,IU,1,0,2)
               ITMAX = MAXT + 1
               IF ((IU .GT. MAXU) .OR. (IV .GT. MAXV)) ITMAX = MAXT
               DO 400 IT = 0, ITMAX
                  ET = ODC(LVALA,LVALB,IT,0,0,1)
                  FT = ODC(LVALA,LVALB,IT,1,0,1)
                  FX = FT*EU*EV
                  FY = ET*FU*EV
                  FZ = ET*EU*FV
                  IADR0 = IADRAU + IT
                  IADRT = IADR0  + 1
                  IADRU = IADR0  + ISTEPU
                  IADRV = IADR0  + ISTEPV
C
C                 ***** LOOP OVER NUCLEI *****
C
                  DO 500 IATOM = 1, NATOMC
                     I0 = 9*(IATOM - 1)
                     SHLINT(INT,I0+1) = SHLINT(INT,I0+1)+FY*AHGTF(IADRU)
     &                                                  +FZ*AHGTF(IADRV)
                     SHLINT(INT,I0+2) = SHLINT(INT,I0+2)-FX*AHGTF(IADRU)
                     SHLINT(INT,I0+3) = SHLINT(INT,I0+3)-FX*AHGTF(IADRV)
                     SHLINT(INT,I0+4) = SHLINT(INT,I0+4)-FY*AHGTF(IADRT)
                     SHLINT(INT,I0+5) = SHLINT(INT,I0+5)+FX*AHGTF(IADRT)
     &                                                  +FZ*AHGTF(IADRV)
                     SHLINT(INT,I0+6) = SHLINT(INT,I0+6)-FY*AHGTF(IADRV)
                     SHLINT(INT,I0+7) = SHLINT(INT,I0+7)-FZ*AHGTF(IADRT)
                     SHLINT(INT,I0+8) = SHLINT(INT,I0+8)-FZ*AHGTF(IADRU)
                     SHLINT(INT,I0+9) = SHLINT(INT,I0+9)+FX*AHGTF(IADRT)
     &                                                  +FY*AHGTF(IADRU)
C                    SHLINT(INT,I0+1) = SHLINT(INT,I0+1)-FX*AHGTF(IADRT)
C                    SHLINT(INT,I0+5) = SHLINT(INT,I0+5)-FY*AHGTF(IADRU)
C                    SHLINT(INT,I0+9) = SHLINT(INT,I0+9)-FZ*AHGTF(IADRV)
                     IADRT = IADRT + NAHGTF
                     IADRU = IADRU + NAHGTF
                     IADRV = IADRV + NAHGTF
  500             CONTINUE
  400          CONTINUE
               IADRAU = IADRAU + ISTEPU
  300       CONTINUE
            IADRAV = IADRAV + ISTEPV
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
      SUBROUTINE DSOKIN(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT,
     &                  NOPTYP,NATOMC)
C
C     K.Ruud, Sep. 2001. Based on DSOINT
C
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),AHGTF(*),
     &          SHLINT(KCKTAB,*)
#include "onecom.h"
#include "nuclei.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
         MAXT = LVALA + LVALB
         MAXU = MVALA + MVALB
         MAXV = NVALA + NVALB
         INT = INT + 1
         IADRAV = 1
         DO 200 IV = 0, MAXV + 1
            IADRAU = IADRAV
            E1V = ODC(NVALA,NVALB,IV,0,0,3)
            F1V = ODC(NVALA,NVALB,IV,0,1,3)
            E2V = ODC(NVALA,NVALB,IV,2,0,3)
            F2V = ODC(NVALA,NVALB,IV,2,1,3)
            IUMAX = MAXU + 1
            IF (IV .GT. MAXV) IUMAX = MAXU
            DO 300 IU = 0, IUMAX
               E1U = ODC(MVALA,MVALB,IU,0,0,2)
               F1U = ODC(MVALA,MVALB,IU,0,1,2)
               E2U = ODC(MVALA,MVALB,IU,2,0,2)
               F2U = ODC(MVALA,MVALB,IU,2,1,2)
               ITMAX = MAXT + 1
               IF ((IU .GT. MAXU) .OR. (IV .GT. MAXV)) ITMAX = MAXT
               DO 400 IT = 0, ITMAX
                  E1T = ODC(LVALA,LVALB,IT,0,0,1)
                  F1T = ODC(LVALA,LVALB,IT,0,1,1)
                  E2T = ODC(LVALA,LVALB,IT,2,0,1)
                  F2T = ODC(LVALA,LVALB,IT,2,1,1)
                  FX = F2T*E1U*E1V + F1T*E2U*E1V + F1T*E1U*E2V
                  FY = E2T*F1U*E1V + E1T*F2U*E1V + E1T*F1U*E2V
                  FZ = E2T*E1U*F1V + E1T*E2U*F1V + E1T*E1U*F2V
                  IADR0 = IADRAU + IT
                  IADRT = IADR0  + 1
                  IADRU = IADR0  + ISTEPU
                  IADRV = IADR0  + ISTEPV
C
C                 ***** LOOP OVER NUCLEI *****
C
                  DO 500 IATOM = 1, NATOMC
                     I0 = 9*(IATOM - 1)
                     SHLINT(INT,I0+1) = SHLINT(INT,I0+1)+FY*AHGTF(IADRU)
     &                                                  +FZ*AHGTF(IADRV)
                     SHLINT(INT,I0+2) = SHLINT(INT,I0+2)-FX*AHGTF(IADRU)
                     SHLINT(INT,I0+3) = SHLINT(INT,I0+3)-FX*AHGTF(IADRV)
                     SHLINT(INT,I0+4) = SHLINT(INT,I0+4)-FY*AHGTF(IADRT)
                     SHLINT(INT,I0+5) = SHLINT(INT,I0+5)+FX*AHGTF(IADRT)
     &                                                  +FZ*AHGTF(IADRV)
                     SHLINT(INT,I0+6) = SHLINT(INT,I0+6)-FY*AHGTF(IADRV)
                     SHLINT(INT,I0+7) = SHLINT(INT,I0+7)-FZ*AHGTF(IADRT)
                     SHLINT(INT,I0+8) = SHLINT(INT,I0+8)-FZ*AHGTF(IADRU)
                     SHLINT(INT,I0+9) = SHLINT(INT,I0+9)+FX*AHGTF(IADRT)
     &                                                  +FY*AHGTF(IADRU)
                     IADRT = IADRT + NAHGTF
                     IADRU = IADRU + NAHGTF
                     IADRV = IADRV + NAHGTF
  500             CONTINUE
  400          CONTINUE
               IADRAU = IADRAU + ISTEPU
  300       CONTINUE
            IADRAV = IADRAV + ISTEPV
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck symupk */
      SUBROUTINE SYMUPK(AINT,WORK,ISYMO,NELMNT)
C
C  290689 Henrik Koch
C
C  Purpose:
C          Symmetry unpack the matrix AINT into WORK
C
C
#include "implicit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
C
      DIMENSION AINT(NELMNT),WORK(NELMNT),IOFFBL(8)
C
#include "symmet.h"
C
#include "ibtfun.h"
C
C     Initialize work array.
C
      CALL DZERO(WORK,NELMNT)
C
C     Setup offset array to the different subblocks.
C
      IOFF = 0
      DO 50 I = 1,MAXREP + 1
         IOFFBL(I) = IOFF
         IOFF = IOFF + NAOS(I)
   50 CONTINUE
C
      IF (ISYMO .EQ. 1) THEN
         IOFF2 = 0
         DO 100 ISYMA = 1,MAXREP + 1
            NDIMA = NAOS(ISYMA)
            IOFF1 = IOFFBL(ISYMA)
            DO 110 J = 1,NDIMA
               DO 120 K = J,NDIMA
                  KT  = IOFF1 + K
                  JT  = IOFF1 + J
                  KJT = KT*(KT-1)/2 + JT
                  WORK(KJT) = AINT(IOFF2 + K*(K-1)/2 + J)
  120          CONTINUE
  110       CONTINUE
            IOFF2 = IOFF2 + NDIMA*(NDIMA + 1)/2
  100    CONTINUE
      ELSE
#if defined (SYS_AIX)
CAIX 4.x compiler fix provided by T.Saue, Nov.12-97
         IREPO = ISYMO - 1
         DO 200 IREPB = 0, MAXREP
            ISYMB = IREPB + 1
            ISYMA = IBTXOR(IREPO,IREPB) + 1
#else
         DO 200 ISYMB = 1, MAXREP + 1
            ISYMA = IBTXOR(ISYMO - 1,ISYMB - 1) + 1
#endif
            IF (ISYMA .GT. ISYMB) THEN
               NDIMA = NAOS(ISYMA)
               NDIMB = NAOS(ISYMB)
               IF (NDIMA.GT.0 .AND. NDIMB.GT.0) THEN
                  IOFF = NPARNU(ISYMO,ISYMA)
                  DO 210 I = 1,NDIMB
                     DO 220 J = 1,NDIMA
                        JT  = IOFFBL(ISYMA) + J
                        JIT = JT*(JT - 1)/2 + IOFFBL(ISYMB) + I
                        WORK(JIT) = AINT(IOFF + I + (J - 1)*NDIMB)
  220                CONTINUE
  210             CONTINUE
               ENDIF
            ENDIF
  200    CONTINUE
      ENDIF
      RETURN
      END
C  /* Deck gosint */
      SUBROUTINE GOSINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,ORIGIN,EXPPI)
C
C     tuh 1993 -- calculate integrals over cos(kr) and sin(kr)
C                 for GOS : generalized oscillator strengths
C     hjaaj 1995 -- modified to cos(kr)/k and sin(kr)/k
C                   (if k .ne. 0)
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D0 = 0.0D0, DP25 = 0.25D0, D1 = 1.0D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,6), ORIGIN(3)
      INTEGER T, U, V
#include "onecom.h"
#include "lmns.h"
      DKX = ORIGIN(1)
      DKY = ORIGIN(2)
      DKZ = ORIGIN(3)
C begin hjaaj mod. 950705
      IF (DKX .NE. D0) THEN
         FX  = D1 / DKX
      ELSE
         FX  = D1
      END IF
      IF (DKY .NE. D0) THEN
         FY  = D1 / DKY
      ELSE
         FY  = D1
      END IF
      IF (DKZ .NE. D0) THEN
         FZ  = D1 / DKZ
      ELSE
         FZ  = D1
      END IF
C end hjaaj mod. 950705
      FACX = SHGTF/EXP(DP25*EXPPI*DKX**2)
      FACY = SHGTF/EXP(DP25*EXPPI*DKY**2)
      FACZ = SHGTF/EXP(DP25*EXPPI*DKZ**2)
      COSX = COS(CORPX*DKX)*FACX
      COSY = COS(CORPY*DKY)*FACY
      COSZ = COS(CORPZ*DKZ)*FACZ
      SINX = SIN(CORPX*DKX)*FACX
      SINY = SIN(CORPY*DKY)*FACY
      SINZ = SIN(CORPZ*DKZ)*FACZ
C
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         CX    = D0
         SX    = D0
         SGNC  = D1
         SGNS  = D1
         FACKX = D1
         DO 200 T = 0, LVALA + LVALB
            IF (MOD(T,2) .EQ. 0) THEN
               CX = CX + SGNC*FACKX*COSX*ODC(LVALA,LVALB,T,0,0,1)
               SX = SX + SGNS*FACKX*SINX*ODC(LVALA,LVALB,T,0,0,1)
               SGNC = - SGNC
            ELSE
               CX = CX + SGNC*FACKX*SINX*ODC(LVALA,LVALB,T,0,0,1)
               SX = SX + SGNS*FACKX*COSX*ODC(LVALA,LVALB,T,0,0,1)
               SGNS = - SGNS
            END IF
            FACKX = FACKX*DKX
  200    CONTINUE
C
         CY    = D0
         SY    = D0
         SGNC  = D1
         SGNS  = D1
         FACKY = D1
         DO 300 U = 0, MVALA + MVALB
            IF (MOD(U,2) .EQ. 0) THEN
               CY = CY + SGNC*FACKY*COSY*ODC(MVALA,MVALB,U,0,0,2)
               SY = SY + SGNS*FACKY*SINY*ODC(MVALA,MVALB,U,0,0,2)
               SGNC = - SGNC
            ELSE
               CY = CY + SGNC*FACKY*SINY*ODC(MVALA,MVALB,U,0,0,2)
               SY = SY + SGNS*FACKY*COSY*ODC(MVALA,MVALB,U,0,0,2)
               SGNS = - SGNS
            END IF
            FACKY = FACKY*DKY
  300    CONTINUE
C
         CZ    = D0
         SZ    = D0
         SGNC  = D1
         SGNS  = D1
         FACKZ = D1
         DO 400 V = 0, NVALA + NVALB
            IF (MOD(V,2) .EQ. 0) THEN
               CZ = CZ + SGNC*FACKZ*COSZ*ODC(NVALA,NVALB,V,0,0,3)
               SZ = SZ + SGNS*FACKZ*SINZ*ODC(NVALA,NVALB,V,0,0,3)
               SGNC = - SGNC
            ELSE
               CZ = CZ + SGNC*FACKZ*SINZ*ODC(NVALA,NVALB,V,0,0,3)
               SZ = SZ + SGNS*FACKZ*COSZ*ODC(NVALA,NVALB,V,0,0,3)
               SGNS = - SGNS
            END IF
            FACKZ = FACKZ*DKZ
  400    CONTINUE
C
         OX = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         OY = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         OZ = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
C
         INT = INT + 1
C 950705-hjaaj: next 6 lines multiplied by FX, FY, or FZ
         SHLINT(INT,1) = SHLINT(INT,1) - FX*CX*OY*OZ
         SHLINT(INT,2) = SHLINT(INT,2) - FY*OX*CY*OZ
         SHLINT(INT,3) = SHLINT(INT,3) - FZ*OX*OY*CZ
         SHLINT(INT,4) = SHLINT(INT,4) - FX*SX*OY*OZ
         SHLINT(INT,5) = SHLINT(INT,5) - FY*OX*SY*OZ
         SHLINT(INT,6) = SHLINT(INT,6) - FZ*OX*OY*SZ
  100 CONTINUE
      RETURN
      END
C  /* Deck sl1int */
      SUBROUTINE SL1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,DIFABX,DIFABY,DIFABZ,EXPPI,
     &                  IORDER,WORK,LWORK)
C
C     K.Ruud, 1994
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
      DIMENSION WORK(LWORK)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,3*(IORDER+1)*(IORDER+2)/2)
#include "onecom.h"
C
      KDX    = 1
      KDY    = KDX    + IORDER + 1
      KDZ    = KDY    + IORDER + 1
      KDX1   = KDZ    + IORDER + 1
      KDY1   = KDX1   + IORDER + 1
      KDZ1   = KDY1   + IORDER + 1
      KLO    = KDZ1   + IORDER + 1
      KMO    = KLO    + (((IORDER + 1)*(IORDER + 2)/2) + 1)/IRAT
      KNO    = KMO    + (((IORDER + 1)*(IORDER + 2)/2) + 1)/IRAT
      KLAST  = KNO    + (((IORDER + 1)*(IORDER + 2)/2) + 1)/IRAT
      IF (KLAST .GT. LWORK) CALL STOPIT('SL1INT',' ',KLAST,LWORK)
      CALL SL1IN1(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,CORPX,
     &            CORPY,CORPZ,DIFABX,DIFABY,DIFABZ,EXPPI,IORDER,
     &            WORK(KDX),WORK(KDY),WORK(KDZ),WORK(KDX1),
     &            WORK(KDY1),WORK(KDZ1),WORK(KLO),WORK(KMO),WORK(KNO))
      RETURN
      END
C  /* Deck sl1in1 */
      SUBROUTINE SL1IN1(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,DIFABX,DIFABY,DIFABZ,EXPPI,
     &                  IORDER,DX,DY,DZ,DX1,DY1,DZ1,LO,
     &                  MO,NO)
C
C     K.Ruud-94
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,3*(IORDER+1)*(IORDER+2)/2),
     &          DX(0:IORDER), DY(0:IORDER), DZ(0:IORDER),
     &          DX1(0:IORDER), DY1(0:IORDER), DZ1(0:IORDER),
     &          LO((IORDER + 1)*(IORDER + 2)/2),
     &          MO((IORDER + 1)*(IORDER + 2)/2),
     &          NO((IORDER + 1)*(IORDER + 2)/2)
#include "onecom.h"
#include "lmns.h"
#include "orgcom.h"
#include "mgsolt.h"
C
      INT = 0
      ADX = CORAX - ORIGIN(1)
      ADY = CORAY - ORIGIN(2)
      ADZ = CORAZ - ORIGIN(3)
C
C     Cartesian integrals
C
      DO 300 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 300 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
C        One-dimensional integrals
C
         DO 400 IO = 0, IORDER
            DX(IO)  = SHGTF*ODC(LVALA,LVALB,0,0,IO,1)
            DY(IO)  = SHGTF*ODC(MVALA,MVALB,0,0,IO,2)
            DZ(IO)  = SHGTF*ODC(NVALA,NVALB,0,0,IO,3)
C
            DX1(IO) = SHGTF*(ODC(LVALA + 1,LVALB,0,0,IO,1) +
     &                       ODC(LVALA    ,LVALB,0,0,IO,1)*ADX)
            DY1(IO) = SHGTF*(ODC(MVALA + 1,MVALB,0,0,IO,2) +
     &                       ODC(MVALA    ,MVALB,0,0,IO,2)*ADY)
            DZ1(IO) = SHGTF*(ODC(NVALA + 1,NVALB,0,0,IO,3) +
     &                       ODC(NVALA    ,NVALB,0,0,IO,3)*ADZ)
 400     CONTINUE
C
C        Three-dimensional integrals
C
         INT = INT + 1
         CALL LMNVAL(IORDER+1,(IORDER + 1)*(IORDER + 2)/2,LO,MO,NO)
         DO 600 I = 1, (IORDER + 1)*(IORDER + 2)/2
            IF (DOALL_MGSOLT) THEN
               IX = 3*(I - 1) + 1
               IY = 3*(I - 1) + 2
               IZ = 3*(I - 1) + 3
            ELSE
               IX = I
               IY = I
               IZ = I
            END IF
            IF (DOALL_MGSOLT .OR. ICOMP_MGSOLT .EQ. 1)
     &           SHLINT(INT,IX) = SHLINT(INT,IX) - DP5*DX(LO(I))
     &                           *(DIFABY*DY(MO(I))*DZ1(NO(I)) -
     &                            DIFABZ*DZ(NO(I))*DY1(MO(I)))
            IF (DOALL_MGSOLT .OR. ICOMP_MGSOLT .EQ. 2)
     &           SHLINT(INT,IY) = SHLINT(INT,IY) - DP5*DY(MO(I))
     &                           *(DIFABZ*DZ(NO(I))*DX1(LO(I)) -
     &                             DIFABX*DX(LO(I))*DZ1(NO(I)))
            IF (DOALL_MGSOLT .OR. ICOMP_MGSOLT .EQ. 3)
     &           SHLINT(INT,IZ) = SHLINT(INT,IZ) - DP5*DZ(NO(I))
     &                           *(DIFABX*DX(LO(I))*DY1(MO(I)) -
     &                            DIFABY*DY(MO(I))*DX1(LO(I)))
 600     CONTINUE
 300  CONTINUE
      RETURN
      END
C  /* Deck sl2int */
      SUBROUTINE SL2INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,DIFABX,DIFABY,DIFABZ,EXPPI,
     &                  IORDER,WORK,LWORK)
C
C     K.Ruud, 1994
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
      DIMENSION WORK(LWORK)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,6*(IORDER+1)*(IORDER+2)/2)
#include "onecom.h"
C
      KDX    = 1
      KDY    = KDX    + IORDER + 1
      KDZ    = KDY    + IORDER + 1
      KDX1   = KDZ    + IORDER + 1
      KDY1   = KDX1   + IORDER + 1
      KDZ1   = KDY1   + IORDER + 1
      KDX2   = KDZ1   + IORDER + 1
      KDY2   = KDX2   + IORDER + 1
      KDZ2   = KDY2   + IORDER + 1
      KLO    = KDZ2   + IORDER + 1
      KMO    = KLO    + (((IORDER + 1)*(IORDER + 2)/2) + 1)/IRAT
      KNO    = KMO    + (((IORDER + 1)*(IORDER + 2)/2) + 1)/IRAT
      KLAST  = KNO    + (((IORDER + 1)*(IORDER + 2)/2) + 1)/IRAT
      IF (KLAST .GT. LWORK) CALL STOPIT('SL2INT',' ',KLAST,LWORK)
      CALL SL2IN1(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,CORPX,
     &            CORPY,CORPZ,DIFABX,DIFABY,DIFABZ,EXPPI,IORDER,
     &            WORK(KDX),WORK(KDY),WORK(KDZ),WORK(KDX1),
     &            WORK(KDY1),WORK(KDZ1),WORK(KDX2),WORK(KDY2),
     &            WORK(KDZ2),WORK(KLO),WORK(KMO),WORK(KNO))
      RETURN
      END
C  /* Deck sl2in1 */
      SUBROUTINE SL2IN1(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT,
     &                  CORPX,CORPY,CORPZ,DIFABX,DIFABY,DIFABZ,EXPPI,
     &                  IORDER,DX,DY,DZ,DX1,DY1,DZ1,DX2,
     &                  DY2,DZ2,LO,MO,NO)
C
C     K.Ruud-94
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D0 = 0.0D0, DP25 = 0.25D0, D2=2.0D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,6*(IORDER+1)*(IORDER+2)/2),
     &          DX(0:IORDER), DY(0:IORDER), DZ(0:IORDER),
     &          DX1(0:IORDER), DY1(0:IORDER), DZ1(0:IORDER),
     &          DX2(0:IORDER), DY2(0:IORDER), DZ2(0:IORDER),
     &          LO((IORDER + 1)*(IORDER + 2)/2),
     &          MO((IORDER + 1)*(IORDER + 2)/2),
     &          NO((IORDER + 1)*(IORDER + 2)/2)
#include "onecom.h"
#include "lmns.h"
#include "orgcom.h"
#include "mgsolt.h"
C
      INT = 0
      ADX = CORAX - ORIGIN(1)
      ADY = CORAY - ORIGIN(2)
      ADZ = CORAZ - ORIGIN(3)
C
C     Cartesian integrals
C
      DO 300 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 300 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
C        One-dimensional integrals
C
         DO 400 IO = 0, IORDER
            DX(IO)  = SHGTF*ODC(LVALA,LVALB,0,0,IO,1)
            DY(IO)  = SHGTF*ODC(MVALA,MVALB,0,0,IO,2)
            DZ(IO)  = SHGTF*ODC(NVALA,NVALB,0,0,IO,3)
            DX1(IO) = SHGTF*(ODC(LVALA + 1,LVALB,0,0,IO,1) +
     &                   ADX*ODC(LVALA    ,LVALB,0,0,IO,1))
            DX2(IO) = SHGTF*(ODC(LVALA + 2,LVALB,0,0,IO,1) +
     &                D2*ADX*ODC(LVALA + 1,LVALB,0,0,IO,1) +
     &               ADX*ADX*ODC(LVALA    ,LVALB,0,0,IO,1))
            DY1(IO) = SHGTF*(ODC(MVALA + 1,MVALB,0,0,IO,2) +
     &                   ADY*ODC(MVALA    ,MVALB,0,0,IO,2))
            DY2(IO) = SHGTF*(ODC(MVALA + 2,MVALB,0,0,IO,2) +
     &                D2*ADY*ODC(MVALA + 1,MVALB,0,0,IO,2) +
     &               ADY*ADY*ODC(MVALA    ,MVALB,0,0,IO,2))
            DZ1(IO) = SHGTF*(ODC(NVALA + 1,NVALB,0,0,IO,3) +
     &                   ADZ*ODC(NVALA    ,NVALB,0,0,IO,3))
            DZ2(IO) = SHGTF*(ODC(NVALA + 2,NVALB,0,0,IO,3) +
     &                D2*ADZ*ODC(NVALA + 1,NVALB,0,0,IO,3) +
     &               ADZ*ADZ*ODC(NVALA    ,NVALB,0,0,IO,3))
  400    CONTINUE
C
C        Three-dimensional integrals
C
         INT = INT + 1
         CALL LMNVAL(IORDER+1,(IORDER + 1)*(IORDER + 2)/2,LO,MO,NO)
         DO 600 I = 1, (IORDER + 1)*(IORDER + 2)/2
            IF (DOALL_MGSOLT) THEN
               IXX = 6*(I - 1) + 1
               IXY = 6*(I - 1) + 2
               IXZ = 6*(I - 1) + 3
               IYY = 6*(I - 1) + 4
               IYZ = 6*(I - 1) + 5
               IZZ = 6*(I - 1) + 6
            ELSE
               IXX = I
               IXY = I
               IXZ = I
               IYY = I
               IYZ = I
               IZZ = I
            END IF
            SX0 = DX(LO(I))
            SY0 = DY(MO(I))
            SZ0 = DZ(NO(I))
            SX1 = DX1(LO(I))
            SY1 = DY1(MO(I))
            SZ1 = DZ1(NO(I))
            SX2 = DX2(LO(I))
            SY2 = DY2(MO(I))
            SZ2 = DZ2(NO(I))
            IF (DOALL_MGSOLT .OR. ICOMP_MGSOLT .EQ. 1)
     &         SHLINT(INT,IXX) = SHLINT(INT,IXX)
     &                         - (D2*DIFABZ*DIFABY*SY1*SZ1*SX0
     &                         - DIFABZ*DIFABZ*SY2*SX0*SZ0
     &                         - DIFABY*DIFABY*SZ2*SX0*SY0)*DP25
            IF (DOALL_MGSOLT .OR. ICOMP_MGSOLT .EQ. 2)
     &         SHLINT(INT,IXY) = SHLINT(INT,IXY)
     &                         - (DIFABZ*DIFABZ*SX1*SY1*SZ0
     &                         - DIFABY*DIFABZ*SZ1*SX1*SY0
     &                         - DIFABZ*DIFABX*SY1*SZ1*SX0
     &                         + DIFABY*DIFABX*SZ2*SX0*SY0)*DP25
            IF (DOALL_MGSOLT .OR. ICOMP_MGSOLT .EQ. 3)
     &         SHLINT(INT,IXZ) = SHLINT(INT,IXZ)
     &                         - (DIFABY*DIFABY*SX1*SZ1*SY0
     &                         - DIFABX*DIFABY*SY1*SZ1*SX0
     &                         + DIFABZ*DIFABX*SY2*SX0*SZ0
     &                         - DIFABY*DIFABZ*SY1*SX1*SZ0)*DP25
            IF (DOALL_MGSOLT .OR. ICOMP_MGSOLT .EQ. 4)
     &         SHLINT(INT,IYY) = SHLINT(INT,IYY)
     &                         - (D2*DIFABZ*DIFABX*SX1*SZ1*SY0
     &                         - DIFABZ*DIFABZ*SX2*SY0*SZ0
     &                         - DIFABX*DIFABX*SZ2*SX0*SY0)*DP25
            IF (DOALL_MGSOLT .OR. ICOMP_MGSOLT .EQ. 5)
     &         SHLINT(INT,IYZ) = SHLINT(INT,IYZ)
     &                         - (DIFABZ*DIFABY*SX2*SY0*SZ0
     &                         - DIFABZ*DIFABX*SX1*SY1*SZ0
     &                         - DIFABY*DIFABX*SX1*SZ1*SY0
     &                         + DIFABX*DIFABX*SY1*SZ1*SX0)*DP25
            IF (DOALL_MGSOLT .OR. ICOMP_MGSOLT .EQ. 6)
     &         SHLINT(INT,IZZ) = SHLINT(INT,IZZ)
     &                         - (D2*DIFABX*DIFABY*SX1*SY1*SZ0
     &                         - DIFABY*DIFABY*SX2*SY0*SZ0
     &                         - DIFABX*DIFABX*SY2*SX0*SZ0)*DP25
  600    CONTINUE
  300 CONTINUE
      RETURN
      END
C
C  /* Deck se1int */
      SUBROUTINE SE1INT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT)
C
C     H. Heiberg, -95
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,3)
#include "onecom.h"
#include "lmns.h"
      INT = 0
C
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         DX0 = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         DY0 = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         DZ0 = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
C
         DPLX = DX0*SY0*SZ0
         DPLY = SX0*DY0*SZ0
         DPLZ = SX0*SY0*DZ0
C
         SHLINT(INT,1) = SHLINT(INT,1) - DPLX
         SHLINT(INT,2) = SHLINT(INT,2) - DPLY
         SHLINT(INT,3) = SHLINT(INT,3) - DPLZ
 100  CONTINUE
      RETURN
      END
C
C  /* Deck se1inb */
      SUBROUTINE SE1INB(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,
     $           SHLINT, EXPA, EXPB)
C
C        H. Heiberg, -95
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,3)
#include "onecom.h"
#include "lmns.h"
      PARAMETER (D4I = 0.25D0)
C
      EXPAIQ = D4I / (EXPA * EXPA)
      EXPBIQ = D4I / (EXPB * EXPB)
      DIFEXP = EXPAIQ - EXPBIQ
C
      INT = 0
C
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
         SX0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         DX0 = DIFEXP*SHGTF*ODC(LVALA,LVALB,0,1,0,1)
         DY0 = DIFEXP*SHGTF*ODC(MVALA,MVALB,0,1,0,2)
         DZ0 = DIFEXP*SHGTF*ODC(NVALA,NVALB,0,1,0,3)
C
         DPLX = DX0*SY0*SZ0
         DPLY = SX0*DY0*SZ0
         DPLZ = SX0*SY0*DZ0
C
         SHLINT(INT,1) = SHLINT(INT,1) - DPLX
         SHLINT(INT,2) = SHLINT(INT,2) - DPLY
         SHLINT(INT,3) = SHLINT(INT,3) - DPLZ
 100  CONTINUE
      RETURN
      END
C
C  /* Deck h1eint */
      SUBROUTINE H1EINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,
     $           SHLINT,AHGTF,NATOMC)
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D4INV = 0.250D0, D2INV = 0.50D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,3), AHGTF(*)
#include "onecom.h"
#include "lmns.h"
      INT = 0
C
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1,KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
         SX0  = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0  = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0  = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         SX1  = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         SY1  = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         SZ1  = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
         DX2  = SHGTF*ODC(LVALA,LVALB,0,2,0,1)
         DY2  = SHGTF*ODC(MVALA,MVALB,0,2,0,2)
         DZ2  = SHGTF*ODC(NVALA,NVALB,0,2,0,3)
         DX21 = SHGTF*ODC(LVALA,LVALB,0,2,1,1)
         DY21 = SHGTF*ODC(MVALA,MVALB,0,2,1,2)
         DZ21 = SHGTF*ODC(NVALA,NVALB,0,2,1,3)
C
C       Kinetic energy contribution
C
         SHLINT(INT,1) = SHLINT(INT,1) - D2INV*( DX21*SY0*SZ0
     &                                         + SX1 *DY2*SZ0
     &                                         + SX1 *SY0*DZ2)
         SHLINT(INT,2) = SHLINT(INT,2) - D2INV*( DX2*SY1 *SZ0
     $                                         + SX0*DY21*SZ0
     $                                         + SX0*SY1 *DZ2)
         SHLINT(INT,3) = SHLINT(INT,3) - D2INV*( DX2*SY0*SZ1
     $                                         + SX0*DY2*SZ1
     $                                         + SX0*SY0*DZ21)
 100  CONTINUE
      RETURN
      END
C  /* Deck relvol */
      SUBROUTINE RELOVL(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT)
C ach
C     a. chr. hennum's  1. subroutine.
C
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D25=0.25D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,6)
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
         DO 100 ICOMPB  = 1, KCKTB
            LVALB = LVALUB(ICOMPB)
            MVALB = MVALUB(ICOMPB)
            NVALB = NVALUB(ICOMPB)
C
C     **********************************************
C     ***** DPT OVERLAP INTEGRALS *****
C     **********************************************
C
C            WRITE(LUPRI,*) '   ACHODCER : '
C            WRITE(LUPRI,*) 'KCKTA B' , KCKTA,KCKTB
C            WRITE(LUPRI,*) 'LVAL A B',LVALA,LVALB
C            WRITE(LUPRI,*) 'MVAL A B',MVALA,MVALB
C            WRITE(LUPRI,*) 'NVAL A B',NVALA,NVALB
C
            X0 = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
            Y0 = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
            Z0 = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
            X1 = SHGTF*ODC(LVALA,LVALB,0,1,0,1)
            Y1 = SHGTF*ODC(MVALA,MVALB,0,1,0,2)
            Z1 = SHGTF*ODC(NVALA,NVALB,0,1,0,3)
            X2 = SHGTF*ODC(LVALA,LVALB,0,2,0,1)
            Y2 = SHGTF*ODC(MVALA,MVALB,0,2,0,2)
            Z2 = SHGTF*ODC(NVALA,NVALB,0,2,0,3)
            INT = INT + 1
            SHLINT(INT,1) = SHLINT(INT,1) + D25*X2*Y0*Z0
            SHLINT(INT,2) = SHLINT(INT,2) + D25*X1*Y1*Z0
            SHLINT(INT,3) = SHLINT(INT,3) + D25*X1*Y0*Z1
            SHLINT(INT,4) = SHLINT(INT,4) + D25*X0*Y2*Z0
            SHLINT(INT,5) = SHLINT(INT,5) + D25*X0*Y1*Z1
            SHLINT(INT,6) = SHLINT(INT,6) + D25*X0*Y0*Z2
C ACH CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C            TOT=SHLINT(INT,1)+SHLINT(INT,4)+SHLINT(INT,6)+
C     %           +2.0D0*(SHLINT(INT,2)+SHLINT(INT,3)+SHLINT(INT,5))
C            WRITE(LUPRI,*) 'ACH TOT ',2.0D0*TOT
C            WRITE(LUPRI,*) '  ACHSHLINTer : '
C            WRITE(LUPRI,*) 'XX',X2*Y0*Z0
C            WRITE(LUPRI,*) 'YY',X0*Y2*Z0
C            WRITE(LUPRI,*) 'ZZ',X0*Y0*Z2
C            WRITE(LUPRI,*) 'XY',X1*Y1*Z0
C            WRITE(LUPRI,*) 'XZ',X1*Y0*Z1
C            WRITE(LUPRI,*) 'YZ',X0*Y1*Z1
C            WRITE(LUPRI,*)
 100  CONTINUE
      RETURN
      END
C
C
C
C
C  /* Deck relpo1 */
      SUBROUTINE RELPO1(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT
     &                  ,AHGTF,NATOMC,EXPA,EXPB)
C
C
C     a. chr. hennum's  2. subroutine.
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D25=0.25D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,6),AHGTF(*)
#include "onecom.h"
#include "lmns.h"
C
      MAXADD = 2
      INT = 0
C      TOT=0.0D0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1,KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         MAXT = LVALA + LVALB + MAXADD
         MAXU = MVALA + MVALB + MAXADD
         MAXV = NVALA + NVALB + MAXADD
         INT = INT + 1
         IADRAV = 1
         DO 200 IV = 0,MAXV
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            FV = ODC(NVALA,NVALB,IV,1,0,3)
            GV = ODC(NVALA,NVALB,IV,2,0,3)
            IADRAU = IADRAV
            DO 300 IU = 0,MAXU
               EU = ODC(MVALA,MVALB,IU,0,0,2)
               FU = ODC(MVALA,MVALB,IU,1,0,2)
               GU = ODC(MVALA,MVALB,IU,2,0,2)
               EE = EU*EV
               FE = FU*EV
               GE = GU*EV
               EF = EU*FV
               FF = FU*FV
               EG = EU*GV
               DO 400 IT = 0,MAXT
                  ET = ODC(LVALA,LVALB,IT,0,0,1)
                  FT = ODC(LVALA,LVALB,IT,1,0,1)
                  GT = ODC(LVALA,LVALB,IT,2,0,1)
                  EEE = ET*EE
                  FEE = FT*EE
                  EFE = ET*FE
                  EEF = ET*EF
                  FFE = FT*FE
                  FEF = FT*EF
                  EFF = ET*FF
                  GEE = GT*EE
                  EGE = ET*GE
                  EEG = ET*EG
                  IADR00 = IADRAU + IT
                  IADR0T = IADR00 + 1
                  IADR0U = IADR00 + ISTEPU
                  IADR0V = IADR00 + ISTEPV
C
                  IADD = - NAHGTF
C
C                 ***** Loop over nuclei *****
C
C                  DO 500 I = 1, NATOMC
C
C                    Pick up HGTF integrals
C
                     IADD = IADD + NAHGTF
                     AH00 = AHGTF(IADR00 + IADD)
                     AH0T = AHGTF(IADR0T + IADD)
                     AH0U = AHGTF(IADR0U + IADD)
                     AH0V = AHGTF(IADR0V + IADD)
C                     write(*,*) '1ACH',FEE,AH0U,EFE,AH0T
C
C     A-C differentiated integrals:
C     ach
C     Adjusted to closed shell CC-module. Won't work as
C     planned. Added Y and Z terms to the X-term to get
C     only one part.
C
                     SHLINT(INT,1)=SHLINT(INT,1)-D25*((FEE*AH0T
     &                    -GEE*AH00) +(EFE*AH0U
     &                    -EGE*AH00) +(EEF*AH0V
     &                    -EEG*AH00))
                     SHLINT(INT,2)=SHLINT(INT,2)-D25*(FEE*AH0U
     &                    +EFE*AH0T-2.0D0*FFE*AH00)
                     SHLINT(INT,3)=SHLINT(INT,3)-D25*(FEE*AH0V
     &                    +EEF*AH0T-2.0D0*FEF*AH00)
                     SHLINT(INT,4)=SHLINT(INT,4)-D25*(EFE*AH0U
     &                    -EGE*AH00)
                     SHLINT(INT,5)=SHLINT(INT,5)-D25*(EFE*AH0V
     &                    +EEF*AH0U-2.0D0*EFF*AH00)
                     SHLINT(INT,6)=SHLINT(INT,6)-D25*(EEF*AH0V
     &                    -EEG*AH00)
C  500             CONTINUE
  400          CONTINUE
               IADRAU = IADRAU + ISTEPU
  300       CONTINUE
            IADRAV = IADRAV + ISTEPV
  200    CONTINUE
  100 CONTINUE
C
C     ACH
C      WRITE(*,*) 'ACH : ', SHLINT(1,2)
      RETURN
      END
C
C  /* Deck relpo2 */
      SUBROUTINE RELPO2(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,SHLINT
     &                  ,AHGTF,NATOMC,EXPA,EXPB)
C     Based on CINT2
C
C     a. chr. hennum's  3. subroutine.
C     Anti-symmetric
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D25=0.25D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,6),AHGTF(*)
#include "onecom.h"
#include "lmns.h"
C
      MAXADD = 2
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1,KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         MAXT = LVALA + LVALB + MAXADD
         MAXU = MVALA + MVALB + MAXADD
         MAXV = NVALA + NVALB + MAXADD
         INT = INT + 1
         IADRAV = 1
         DO 200 IV = 0,MAXV
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            FV = ODC(NVALA,NVALB,IV,1,0,3)
            GV = ODC(NVALA,NVALB,IV,2,0,3)
            IADRAU = IADRAV
            DO 300 IU = 0,MAXU
               EU = ODC(MVALA,MVALB,IU,0,0,2)
               FU = ODC(MVALA,MVALB,IU,1,0,2)
               GU = ODC(MVALA,MVALB,IU,2,0,2)
               EE = EU*EV
               FE = FU*EV
               GE = GU*EV
               EF = EU*FV
               FF = FU*FV
               EG = EU*GV
               DO 400 IT = 0,MAXT
                  ET = ODC(LVALA,LVALB,IT,0,0,1)
                  FT = ODC(LVALA,LVALB,IT,1,0,1)
                  GT = ODC(LVALA,LVALB,IT,2,0,1)
                  EEE = ET*EE
                  FEE = FT*EE
                  EFE = ET*FE
                  EEF = ET*EF
                  FFE = FT*FE
                  FEF = FT*EF
                  EFF = ET*FF
                  GEE = GT*EE
                  EGE = ET*GE
                  EEG = ET*EG
                  IADR00 = IADRAU + IT
                  IADR0T = IADR00 + 1
                  IADR0U = IADR00 + ISTEPU
                  IADR0V = IADR00 + ISTEPV
C                  IADRTT = IADR0T + 1
C                  IADRTU = IADR0T + ISTEPU
C                  IADRTV = IADR0T + ISTEPV
C                  IADRUU = IADR0U + ISTEPU
C                  IADRUV = IADR0U + ISTEPV
C                  IADRVV = IADR0V + ISTEPV
C
                  IADD = - NAHGTF
C
C                 ***** Loop over nuclei *****
C
C                  DO 500 I = 1, NATOMC
C
C                    Pick up HGTF integrals
C
                     IADD = IADD + NAHGTF
                     AH00 = AHGTF(IADR00 + IADD)
                     AH0T = AHGTF(IADR0T + IADD)
                     AH0U = AHGTF(IADR0U + IADD)
                     AH0V = AHGTF(IADR0V + IADD)
C
C                     write(*,*) 'ACH',FEE,AH0U,EFE,AH0T
                     SHLINT(INT,1)=SHLINT(INT,1)-D25*(FEE*AH0U
     &                    -EFE*AH0T )
                     SHLINT(INT,2)=SHLINT(INT,2)-D25*(FEE*AH0V
     &                    -EEF*AH0T )
                     SHLINT(INT,3)=SHLINT(INT,3)-D25*(EFE*AH0V
     &                    -EEF*AH0U )
C  500             CONTINUE
  400          CONTINUE
               IADRAU = IADRAU + ISTEPU
  300       CONTINUE
            IADRAV = IADRAV + ISTEPV
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck nsyint */
      SUBROUTINE NSYINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT,
     &                  NOPTYP,NATOMC)
C
C     DPT (x/r^3 d/dx) (square)
C
C     ACh 2.3.98
C
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "priunit.h"
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,NOPTYP), AHGTF(*)
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         MAXT = LVALA + LVALB + 1
         MAXU = MVALA + MVALB + 1
         MAXV = NVALA + NVALB + 1
         INT=INT+1
         IADRAV = 1
         DO 200 IV = 0,MAXV
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            FV = ODC(NVALA,NVALB,IV,1,0,3)
            IADRAU = IADRAV
            DO 300 IU = 0,MAXU
               EU = ODC(MVALA,MVALB,IU,0,0,2)
               FU = ODC(MVALA,MVALB,IU,1,0,2)
               EE = EU*EV
               FE = FU*EV
               EF = EU*FV
               DO 400 IT = 0,MAXT
                  ET = ODC(LVALA,LVALB,IT,0,0,1)
                  FT = ODC(LVALA,LVALB,IT,1,0,1)
                  FEE = FT*EE
                  EFE = ET*FE
                  EEF = ET*EF
                  IOFF = 0
                  DO 500 IATOM = 1, NATOMC
                     IXX = 9*(IATOM - 1) + 1
                     IXY = 9*(IATOM - 1) + 2
                     IXZ = 9*(IATOM - 1) + 3
                     IYX = 9*(IATOM - 1) + 4
                     IYY = 9*(IATOM - 1) + 5
                     IYZ = 9*(IATOM - 1) + 6
                     IZX = 9*(IATOM - 1) + 7
                     IZY = 9*(IATOM - 1) + 8
                     IZZ = 9*(IATOM - 1) + 9
                     AH0T = AHGTF(IOFF + IADRAU + IT + 1)
                     AH0U = AHGTF(IOFF + IADRAU + IT + ISTEPU)
                     AH0V = AHGTF(IOFF + IADRAU + IT + ISTEPV)
                     SHLINT(INT,IXX) = SHLINT(INT,IXX)+ FEE * AH0T
                     SHLINT(INT,IXY) = SHLINT(INT,IXY)+ EFE * AH0T
                     SHLINT(INT,IXZ) = SHLINT(INT,IXZ)+ EEF * AH0T
                     SHLINT(INT,IYX) = SHLINT(INT,IYX)+ FEE * AH0U
                     SHLINT(INT,IYY) = SHLINT(INT,IYY)+ EFE * AH0U
                     SHLINT(INT,IYZ) = SHLINT(INT,IYZ)+ EEF * AH0U
                     SHLINT(INT,IZX) = SHLINT(INT,IZX)+ FEE * AH0V
                     SHLINT(INT,IZY) = SHLINT(INT,IZY)+ EFE * AH0V
                     SHLINT(INT,IZZ) = SHLINT(INT,IZZ)+ EEF * AH0V
                     IOFF = IOFF + NAHGTF
  500             CONTINUE
  400          CONTINUE
               IADRAU = IADRAU + ISTEPU
  300       CONTINUE
            IADRAV = IADRAV + ISTEPV
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck cntfpv */
      SUBROUTINE CNTFPV(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,SHLINT)
C
C     Parity violating electroweak interaction.
C
C     The array AHGTF contains Hermite integrals over a weighted sum
C     of nuclear attraction integrals . The weights are the nucleus-
C     dependent effective electroweak charges (enhancement factors) Q(Z):
C
C           Q(Z) = [N(neutrons) - (1 - 4*sin(theta_W)**2) * Z] ,
C
C     where theta_W is the Weinberg angle, Z the nuclear charge,
C     and N(neutrons) the number of neutrons in the nucleus.
C
C     See:
C     P. Lazzeretti and R. Zanasi, Chem. Phys. Lett. 279, 349 (1997).
C     A. Bakasov, T.-K. Ha, and M. Quack, J. Chem. Phys. 109, 7263 (1998).
C
C     Wim Klopper, Utrecht University, October 18, 1999
C
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"
      PARAMETER (D1 = 1D0, D2 = 2.D0, D16 = 16D0)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          AHGTF(*), SHLINT(KCKTAB,3)
#include "onecom.h"
#include "lmns.h"
      FFA = - D1 / (D16 * DATAN(D1))
      INT = 1
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         MAXT = LVALA + LVALB + 1
         MAXU = MVALA + MVALB + 1
         MAXV = NVALA + NVALB + 1
         IADRAV = 1
         DO 200 IV = 0,MAXV
            EV = FFA * ODC(NVALA,NVALB,IV,0,0,3)
            FV = FFA * ODC(NVALA,NVALB,IV,1,0,3)
            IADRAU = IADRAV
            DO 300 IU = 0,MAXU
               EU = ODC(MVALA,MVALB,IU,0,0,2)
               FU = ODC(MVALA,MVALB,IU,1,0,2)
               EE = EU*EV
               FE = FU*EV
               EF = EU*FV
               DO 400 IT = 0,MAXT
                  ET = ODC(LVALA,LVALB,IT,0,0,1)
                  FT = ODC(LVALA,LVALB,IT,1,0,1)
                  EEE = ET*EE
                  FEE = FT*EE
                  EFE = ET*FE
                  EEF = ET*EF
                  AH00 = (AHGTF(IADRAU + IT + 2 )
     &                 +  AHGTF(IADRAU + IT + 2 * ISTEPU)
     &                 +  AHGTF(IADRAU + IT + 2 * ISTEPV)) * D2
                  AH0T =  AHGTF(IADRAU + IT + 3)
     &                 +  AHGTF(IADRAU + IT + 2 * ISTEPU + 1)
     &                 +  AHGTF(IADRAU + IT + 2 * ISTEPV + 1)
                  AH0U =  AHGTF(IADRAU + IT + 2 + ISTEPU)
     &                 +  AHGTF(IADRAU + IT + 3 * ISTEPU)
     &                 +  AHGTF(IADRAU + IT + 2 * ISTEPV + ISTEPU)
                  AH0V =  AHGTF(IADRAU + IT + 2 + ISTEPV)
     &                 +  AHGTF(IADRAU + IT + 2 * ISTEPU + ISTEPV)
     &                 +  AHGTF(IADRAU + IT + 3 * ISTEPV)
                  SHLINT(INT,1) = SHLINT(INT,1) + FEE*AH00 - EEE*AH0T
                  SHLINT(INT,2) = SHLINT(INT,2) + EFE*AH00 - EEE*AH0U
                  SHLINT(INT,3) = SHLINT(INT,3) + EEF*AH00 - EEE*AH0V
  400          CONTINUE
               IADRAU = IADRAU + ISTEPU
  300       CONTINUE
            IADRAV = IADRAV + ISTEPV
  200    CONTINUE
         INT = INT + 1
  100 CONTINUE
      RETURN
      END
C
cLig >> added the RAMINT and RMIINT subroutines
C  /* Deck ramint */
      SUBROUTINE RAMINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,
     &                  SHLINT)
C
C     Calculation of the (r-r')l' integrals for determination of
C     magnetizability in an analytical way
C
C     A. Ligabue Nov. 1999
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,9)
#include "onecom.h"
#include "lmns.h"
      INT = 0
      DO 100 ICOMPA = 1, KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         INT = INT + 1
C
         SX0  = SHGTF*ODC(LVALA,LVALB,0,0,0,1)
         SY0  = SHGTF*ODC(MVALA,MVALB,0,0,0,2)
         SZ0  = SHGTF*ODC(NVALA,NVALB,0,0,0,3)
         SX1  = SHGTF*ODC(LVALA,LVALB,0,0,1,1)
         SY1  = SHGTF*ODC(MVALA,MVALB,0,0,1,2)
         SZ1  = SHGTF*ODC(NVALA,NVALB,0,0,1,3)
         SX2  = SHGTF*ODC(LVALA,LVALB,0,0,2,1)
         SY2  = SHGTF*ODC(MVALA,MVALB,0,0,2,2)
         SZ2  = SHGTF*ODC(NVALA,NVALB,0,0,2,3)
         SX11 = SHGTF*ODC(LVALA,LVALB,0,1,1,1)
         SY11 = SHGTF*ODC(MVALA,MVALB,0,1,1,2)
         SZ11 = SHGTF*ODC(NVALA,NVALB,0,1,1,3)
         DX0  = SHGTF*ODC(LVALA,LVALB,0,1,0,1)
         DY0  = SHGTF*ODC(MVALA,MVALB,0,1,0,2)
         DZ0  = SHGTF*ODC(NVALA,NVALB,0,1,0,3)
C
         SHLINT(INT,1) = SHLINT(INT,1) + SX1*(SY1*DZ0 - DY0*SZ1)
         SHLINT(INT,2) = SHLINT(INT,2) + SY0*(SZ1*SX11 - SX2*DZ0)
         SHLINT(INT,3) = SHLINT(INT,3) + SZ0*(SX2*DY0 - SY1*SX11)
         SHLINT(INT,4) = SHLINT(INT,4) + SX0*(SY2*DZ0 - SZ1*SY11)
         SHLINT(INT,5) = SHLINT(INT,5) + SY1*(SZ1*DX0 - DZ0*SX1)
         SHLINT(INT,6) = SHLINT(INT,6) + SZ0*(SX1*SY11 - SY2*DX0)
         SHLINT(INT,7) = SHLINT(INT,7) + SX0*(SY1*SZ11 - SZ2*DY0)
         SHLINT(INT,8) = SHLINT(INT,8) + SY0*(SZ2*DX0 - SX1*SZ11)
         SHLINT(INT,9) = SHLINT(INT,9) + SZ1*(SX1*DY0 - SY1*DX0)
C
 100  CONTINUE
      RETURN
      END
C
C  /* Deck rmiint */
C
      SUBROUTINE RMIINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AHGTF,
     &                  SHLINT,NOPTYP,NATOMC)
C
C     Calculation of the integrals (r-r')l'/|r-R_I|**3 for the
C     determination of the nuclear shieldings in an analytical way
C
C     A. Ligabue, Nov. 1999
C
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          SHLINT(KCKTAB,NOPTYP), AHGTF(*)
#include "onecom.h"
#include "lmns.h"
      INT = 1
      DO 100 ICOMPA = 1, KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
         MAXT = LVALA + LVALB + 2
         MAXU = MVALA + MVALB + 2
         MAXV = NVALA + NVALB + 2
         IADRAV = 1
         DO 200 IV = 0, MAXV
            S0Z = ODC(NVALA,NVALB,IV,0,0,3)
            D1Z = ODC(NVALA,NVALB,IV,1,1,3)
            D0Z = ODC(NVALA,NVALB,IV,1,0,3)
            S1Z = ODC(NVALA,NVALB,IV,0,1,3)
            IADRAU = IADRAV
            DO 300 IU = 0, MAXU
               S0Y = ODC(MVALA,MVALB,IU,0,0,2)
               D1Y = ODC(MVALA,MVALB,IU,1,1,2)
               D0Y = ODC(MVALA,MVALB,IU,1,0,2)
               S1Y = ODC(MVALA,MVALB,IU,0,1,2)
               DO 400 IT = 0, MAXT
                  S0X = ODC(LVALA,LVALB,IT,0,0,1)
                  D1X = ODC(LVALA,LVALB,IT,1,1,1)
                  D0X = ODC(LVALA,LVALB,IT,1,0,1)
                  S1X = ODC(LVALA,LVALB,IT,0,1,1)
                  IOFF = 0
                  DO 500 IATOM = 1, NATOMC
c                     IXX = 9*(IATOM - 1) + 1
c                     IXY = 9*(IATOM - 1) + 2
c                     IXZ = 9*(IATOM - 1) + 3
c                     IYX = 9*(IATOM - 1) + 4
c                     IYY = 9*(IATOM - 1) + 5
c                     IYZ = 9*(IATOM - 1) + 6
c                     IZX = 9*(IATOM - 1) + 7
c                     IZY = 9*(IATOM - 1) + 8
c                     IZZ = 9*(IATOM - 1) + 9
                     IXX = 9*(IATOM - 1) + 1
                     IYX = 9*(IATOM - 1) + 2
                     IZX = 9*(IATOM - 1) + 3
                     IXY = 9*(IATOM - 1) + 4
                     IYY = 9*(IATOM - 1) + 5
                     IZY = 9*(IATOM - 1) + 6
                     IXZ = 9*(IATOM - 1) + 7
                     IYZ = 9*(IATOM - 1) + 8
                     IZZ = 9*(IATOM - 1) + 9
                     AH0X = AHGTF(IOFF + IADRAU + IT + 1)
                     AH0Y = AHGTF(IOFF + IADRAU + IT + ISTEPU)
                     AH0Z = AHGTF(IOFF + IADRAU + IT + ISTEPV)
                     SHLINT(INT,IXX) = SHLINT(INT,IXX) + S1X*
     &                                 (S0Y*AH0Y*D0Z - D0Y*S0Z*AH0Z)
                     SHLINT(INT,IXY) = SHLINT(INT,IXY) + S0Y*
     &                                 (D1X*S0Z*AH0Z - S1X*AH0X*D0Z)
                     SHLINT(INT,IXZ) = SHLINT(INT,IXZ) + S0Z*
     &                                 (S1X*AH0X*D0Y - D1X*S0Y*AH0Y)
                     SHLINT(INT,IYX) = SHLINT(INT,IYX) + S0X*
     &                                 (S1Y*AH0Y*D0Z - D1Y*S0Z*AH0Z)
                     SHLINT(INT,IYY) = SHLINT(INT,IYY) + S1Y*
     &                                 (D0X*S0Z*AH0Z - S0X*AH0X*D0Z)
                     SHLINT(INT,IYZ) = SHLINT(INT,IYZ) + S0Z*
     &                                 (S0X*AH0X*D1Y - S1Y*AH0Y*D0X)
                     SHLINT(INT,IZX) = SHLINT(INT,IZX) + S0X*
     &                                 (S0Y*AH0Y*D1Z - S1Z*AH0Z*D0Y)
                     SHLINT(INT,IZY) = SHLINT(INT,IZY) + S0Y*
     &                                 (D0X*S1Z*AH0Z - S0X*AH0X*D1Z)
                     SHLINT(INT,IZZ) = SHLINT(INT,IZZ) + S1Z*
     &                                 (S0X*AH0X*D0Y - D0X*S0Y*AH0Y)

                     IOFF = IOFF + NAHGTF
 500              CONTINUE
 400           CONTINUE
               IADRAU = IADRAU + ISTEPU
 300        CONTINUE
            IADRAV = IADRAV + ISTEPV
 200     CONTINUE
         INT = INT + 1
 100  CONTINUE
      RETURN
      END
cLig <<
