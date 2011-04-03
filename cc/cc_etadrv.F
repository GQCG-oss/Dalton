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
c /* deck CC_ETADRV */
*=====================================================================*
      SUBROUTINE CC_ETADRV(TYPE,LABELV,ISYMS,ISTAT,EIGV,
     &                     ISYMO,FREQS,ICAU,NVEC,MAXVEC,
     &                     WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: calculate response eta vectors, used to build the
*             right-hand-side vectors for the lagrangian multipliers 
*             and as intermediates in the hyperpolarizability
*             and n-photon-transition matrix calculations
*
*             for excited states the X vectors are identical to the
*             rhs vectors for the left eigenvector response equations
*
*     implemented:  L:   ORDER = 2, 3
*                   LE:  ORDER = 1, 2
*                   CL:  ORDER = 2
*
*     Written by Christof Haettig april/june/july 1997.
*     extensions for Cauchy eta vectors in March 1998.
*     adapted for CC-R12 by Christian Neiss, june 2005
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <cclists.h>
#include <dummy.h>
#include <r12int.h>

* local parameters:
      CHARACTER*(19) MSGDBG
      PARAMETER (MSGDBG = '[debug] CC_ETADRV> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE. )

      CHARACTER TYPE*(*)

      INTEGER NVEC, MAXVEC, LWORK
      INTEGER ISYMO(MAXVEC,*), ICAU(MAXVEC,*)
      INTEGER ISYMS(MAXVEC,*), ISTAT(MAXVEC,*)

      CHARACTER*8 LABELV(MAXVEC,*)
      CHARACTER*3 APROXR12

#if defined (SYS_CRAY)
      REAL FREQS(MAXVEC,*), EIGV(MAXVEC,*)
      REAL WORK(LWORK)
      REAL ZERO
      REAL DDOT, XNORM, RDUM
#else
      DOUBLE PRECISION FREQS(MAXVEC,*), EIGV(MAXVEC,*)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION ZERO
      DOUBLE PRECISION DDOT, XNORM, RDUM
#endif
      PARAMETER (ZERO = 0.0d0)

      CHARACTER MODEL*(10)
      INTEGER MX0KTRAN, MX1GTRAN, MX2FTRAN, MX1FATRAN
      INTEGER MXTRAN, MX0GTRAN, MX1FTRAN, MX0FATRAN, MXXETRAN
      INTEGER K0KTRAN, K1GTRAN, K2FTRAN, K1FATRAN
      INTEGER K0GTRAN, K1FTRAN, K0FATRAN, KXETRAN
      INTEGER N0KTRAN, N1GTRAN, N2FTRAN, N1FATRAN
      INTEGER N0GTRAN, N1FTRAN, N0FATRAN, NEATRAN
      INTEGER IOPT, ISYM, IVEC, ORDER, MPERM, NSTAT, IOPTH, IOPTRW
      INTEGER KEND0, LEND0, KEND1, LEND1, LMAX1, LMAX2, KCHI1, KCHI2
      INTEGER KEND2, LEND2, KRHS1, KRHS2, IDUM
      INTEGER KCHIR12, LMAXR12, IOPTRWR12, MODLEN, KRHSR12

* external functions
      INTEGER ILSTSYM

*---------------------------------------------------------------------*
* check number of required eta/rhs vectors, if zero return immediatly:
*---------------------------------------------------------------------*
      IF (NVEC.EQ.0) RETURN

*---------------------------------------------------------------------*
* print header for eta/rhs vector section
*---------------------------------------------------------------------*
      WRITE (LUPRI,'(7(/1X,2A),/)')
     & '------------------------------------',
     &                               '-------------------------------',
     & '|                 OUTPUT FROM ETA/RH',    
     &                               'S VECTOR SECTION              |',
     & '------------------------------------',
     &                               '-------------------------------' 
      CALL FLSHFO(LUPRI)

*---------------------------------------------------------------------*
      IF (.NOT. (CCS .OR. CC2 .OR. CCSD) ) THEN
         CALL QUIT('CC_ETADRV called for unknown Coupled Cluster.')
      END IF

      IF (TYPE(1:3).EQ.'X1 ') THEN
        WRITE (LUPRI,*) 'X1 vectors not implemented in CC_ETADRV,'
        WRITE (LUPRI,*) 'routine CCRHSVEC should be used instead.'
        CALL QUIT('X1 vectors not implemented in CC_ETADRV.')
      ELSE IF (TYPE(1:2).EQ.'X2') THEN
        ORDER = 2
        NSTAT = 0
        MPERM = 2
      ELSE IF (TYPE(1:2).EQ.'X3') THEN
        ORDER = 3
        NSTAT = 0
        MPERM = 6
C     ELSE IF (TYPE(1:2).EQ.'X4') THEN
C       ORDER = 4
C       NSTAT = 0
C       MPERM = ??
      ELSE IF (TYPE(1:3).EQ.'EX1') THEN
        ORDER = 1
        NSTAT = 1
        MPERM = 1
      ELSE IF (TYPE(1:3).EQ.'EX2') THEN
        ORDER = 2
        NSTAT = 1
        MPERM = 2
        WRITE (LUPRI,*) 'warning: X vectors ',TYPE(1:3),
     &       ' not tested!!!.'
      ELSE IF (TYPE(1:3).EQ.'CX2') THEN
        ORDER = 2
        NSTAT = 0
        MPERM = 2
      ELSE
        WRITE (LUPRI,*) 'rhs vectors ',TYPE(1:2),' not implemented.'
        CALL QUIT('required rhs vectors not implemented.')
      END IF


* print some debug/info output
      IF (IPRINT .GT. 10) WRITE(LUPRI,*) 'CC_ETADRV Workspace:',LWORK
  
*---------------------------------------------------------------------*
* allocate & initialize work space for lists
*---------------------------------------------------------------------*

      MXTRAN  = MPERM * NVEC

      MX0KTRAN  = 5 * MXTRAN
      MX0GTRAN  = 4 * MXTRAN
      MX1GTRAN  = 4 * MXTRAN
      MX1FTRAN  = 3 * MXTRAN
      MX2FTRAN  = 3 * MXTRAN
      MX0FATRAN = 5 * MXTRAN
      MX1FATRAN = 5 * MXTRAN
      MXXETRAN  = MXDIM_XEVEC * MXTRAN

      K0KTRAN  = 1
      K0GTRAN  = K0KTRAN  + MX0KTRAN
      K1GTRAN  = K0GTRAN  + MX0GTRAN
      K1FTRAN  = K1GTRAN  + MX1GTRAN
      K2FTRAN  = K1FTRAN  + MX1FTRAN
      K0FATRAN = K2FTRAN  + MX2FTRAN
      K1FATRAN = K0FATRAN + MX0FATRAN
      KXETRAN  = K1FATRAN + MX1FATRAN
      KEND0    = KXETRAN  + MXXETRAN
      LEND0    = LWORK - KEND0

      IF (LEND0 .LT. 0 ) THEN
        CALL QUIT('Insufficient work space in CC_ETADRV.')
      END IF

*---------------------------------------------------------------------*
* set up lists for G, F and F{A} transformations and ETA{O} vectors:
*---------------------------------------------------------------------*
      CALL CC_ETA_SETUP(TYPE,NSTAT,ORDER,LABELV,ISTAT,EIGV,FREQS,ICAU,
     &                  NVEC, MAXVEC,  MXTRAN, 
     &                  WORK(K0KTRAN), N0KTRAN,
     &                  WORK(K0GTRAN), N0GTRAN,
     &                  WORK(K1GTRAN), N1GTRAN,
     &                  WORK(K1FTRAN), N1FTRAN,
     &                  WORK(K2FTRAN), N2FTRAN,
     &                  WORK(K0FATRAN),N0FATRAN,
     &                  WORK(K1FATRAN),N1FATRAN,
     &                  WORK(KXETRAN), NEATRAN  )

*---------------------------------------------------------------------*
* initialize ETA vector files:
*---------------------------------------------------------------------*
      LMAX1 = 0
      LMAX2 = 0
      LMAXR12 = 0
      DO ISYM = 1, NSYM
        LMAX1 = MAX(LMAX1,NT1AM(ISYM))
        LMAX2 = MAX(LMAX2,NT2AM(ISYM))
        IF (CCR12) LMAXR12 = MAX(LMAXR12,NTR12AM(ISYM))
      END DO

      KCHI1 = KEND0
      KCHI2 = KCHI1 + LMAX1
      KCHIR12 = KCHI2 + LMAXR12
      KEND1 = KCHIR12 + LMAX2
      LEND1 = LWORK - KEND1

      IF (LEND1 .LT. 0 ) THEN
        CALL QUIT('Insufficient work space in CC_ETADRV.')
      END IF

      CALL DZERO(WORK(KCHI1),LMAX1)
      IF (.NOT.CCS) CALL DZERO(WORK(KCHI2),LMAX2)
      IF (CCR12) CALL DZERO(WORK(KCHIR12),LMAXR12)

      IF (CCS) THEN
         MODEL  = 'CCS       '
         IOPTRW = 1
      ELSE IF (CC2) THEN
         MODEL  = 'CC2       '
         IOPTRW = 3
      ELSE IF (CCSD) THEN
         MODEL  = 'CCSD      '
         IOPTRW = 3
      ELSE
         CALL QUIT('Unknown coupled cluster model in CC_ETADRV.')
      END IF
      IF (CCR12) THEN
        APROXR12 = '   '
        IOPTRWR12 = 32
      END IF
      CALL CCSD_MODEL(MODEL,MODLEN,10,MODEL,10,APROXR12)

      DO IVEC = 1, NVEC
        ISYM = ILSTSYM(TYPE,IVEC)
        CALL CC_WRRSP(TYPE,IVEC,ISYM,IOPTRW,MODEL,IDUMMY,
     &                WORK(KCHI1),WORK(KCHI2),WORK(KEND1),LEND1)
        IF (CCR12) THEN
          CALL CC_WRRSP(TYPE,IVEC,ISYM,IOPTRWR12,MODEL,IDUMMY,
     &                  IDUMMY,WORK(KCHIR12),WORK(KEND1),LEND1)
        END IF
      END DO

*---------------------------------------------------------------------*
* calculate H matrix contributions:
*---------------------------------------------------------------------*
      IF (TYPE(1:2).EQ.'X3') THEN
        IOPTH = 4
        CALL CC_HMAT('L0','R1','R1','R1',TYPE,N0KTRAN, 0,
     &               WORK(K0KTRAN),IDUMMY,IDUMMY,
     &               WORK(KEND0), LEND0, IOPTH )
      END IF

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 
     &        MSGDBG, 'NORM^2 of ETA vectors after H matrix terms:'
        DO IVEC = 1, NVEC
          ISYM = ILSTSYM(TYPE,IVEC)
          CALL CC_RDRSP(TYPE,IVEC,ISYM,IOPTRW,MODEL,
     &                  WORK(KCHI1),WORK(KCHI2))
          IF (CCR12) CALL CC_RDRSP(TYPE,IVEC,ISYM,IOPTRWR12,MODEL,
     &                             DUMMY,WORK(KCHIR12))
          XNORM = DDOT(NT1AM(ISYM),WORK(KCHI1),1,WORK(KCHI1),1)
          IF (.NOT.CCS) 
     &     XNORM = XNORM+DDOT(NT2AM(ISYM),WORK(KCHI2),1,WORK(KCHI2),1)
          IF (CCR12) 
     &     XNORM = XNORM+DDOT(NTR12AM(ISYM),WORK(KCHIR12),1,
     &                        WORK(KCHIR12),1)
          WRITE (LUPRI,*) MSGDBG, IVEC,XNORM
        END DO
      END IF

*---------------------------------------------------------------------*
* calculate G matrix contributions:
*---------------------------------------------------------------------*
      IF (TYPE(1:2).EQ.'X2') THEN
        IOPT = 4
        CALL CC_GMATRIX('L0 ','R1 ','R1 ',TYPE,N0GTRAN, 0, 
     &                  WORK(K0GTRAN),IDUM,RDUM,WORK(KEND0),LEND0,IOPT)
      ELSE IF (TYPE(1:2).EQ.'X3') THEN
        IOPT = 4
        CALL CC_GMATRIX('L0 ','R2 ','R1 ',TYPE,N0GTRAN, 0, 
     &                  WORK(K0GTRAN),IDUM,RDUM,WORK(KEND0),LEND0,IOPT)
        IOPT = 4
        CALL CC_GMATRIX('L1 ','R1 ','R1 ',TYPE,N1GTRAN, 0, 
     &                  WORK(K1GTRAN),IDUM,RDUM,WORK(KEND0),LEND0,IOPT)
      ELSE IF (TYPE(1:3).EQ.'EX2') THEN
        IOPT = 4
        CALL CC_GMATRIX('LE ','R1 ','R1 ',TYPE,N0GTRAN, 0, 
     &                  WORK(K0GTRAN),IDUM,RDUM,WORK(KEND0),LEND0,IOPT)
      ELSE IF (TYPE(1:3).EQ.'CX2') THEN
        IOPT = 4
        CALL CC_GMATRIX('L0 ','RC ','RC ',TYPE,N0GTRAN, 0, 
     &                  WORK(K0GTRAN),IDUM,RDUM,WORK(KEND0),LEND0,IOPT)
      END IF

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) MSGDBG, 
     &        'NORM^2 of ETA vectors after G matrix terms:'
        DO IVEC = 1, NVEC
          ISYM = ILSTSYM(TYPE,IVEC)
          CALL CC_RDRSP(TYPE,IVEC,ISYM,IOPTRW,MODEL,
     &                  WORK(KCHI1),WORK(KCHI2))
          IF (CCR12.AND..NOT.(CCS.OR.CC2)) 
     &      CALL CC_RDRSP(TYPE,IVEC,ISYM,IOPTRWR12,MODEL,
     &                    DUMMY,WORK(KCHIR12))
          XNORM = DDOT(NT1AM(ISYM),WORK(KCHI1),1,WORK(KCHI1),1)
          IF (.NOT.CCS) 
     &     XNORM = XNORM+DDOT(NT2AM(ISYM),WORK(KCHI2),1,WORK(KCHI2),1)
          IF (CCR12.AND..NOT.(CCS.OR.CC2))
     &     XNORM = XNORM+DDOT(NTR12AM(ISYM),WORK(KCHIR12),1,
     &                        WORK(KCHIR12),1)
          WRITE (LUPRI,*) MSGDBG, IVEC,XNORM
        END DO
      END IF

*---------------------------------------------------------------------*
* calculate F matrix contributions:
*---------------------------------------------------------------------*
      IF (TYPE(1:2).EQ.'X2') THEN
        IOPT = 4
        CALL CC_FMATRIX(WORK(K1FTRAN),N1FTRAN,'L1 ','R1 ',IOPT,TYPE,
     &                  IDUM, RDUM, 0, WORK(KEND0), LEND0)
      ELSE IF (TYPE(1:3).EQ.'CX2') THEN
        IOPT = 4
        CALL CC_FMATRIX(WORK(K1FTRAN),N1FTRAN,'LC ','RC ',IOPT,TYPE,
     &                  IDUM, RDUM, 0, WORK(KEND0), LEND0)
      ELSE IF (TYPE(1:2).EQ.'X3') THEN
        IOPT = 4
        CALL CC_FMATRIX(WORK(K1FTRAN),N1FTRAN,'L1 ','R2 ',IOPT,TYPE,
     &                  IDUM, RDUM, 0, WORK(KEND0), LEND0)
        IOPT = 4
        CALL CC_FMATRIX(WORK(K2FTRAN),N2FTRAN,'L2 ','R1 ',IOPT,TYPE,
     &                  IDUM, RDUM, 0, WORK(KEND0), LEND0)
      ELSE IF (TYPE(1:3).EQ.'EX2') THEN
        IOPT = 4
        CALL CC_FMATRIX(WORK(K1FTRAN),N1FTRAN,'LE ','R2 ',IOPT,TYPE,
     &                  IDUM, RDUM, 0, WORK(KEND0), LEND0)
        IOPT = 4
        CALL CC_FMATRIX(WORK(K2FTRAN),N2FTRAN,'EL1','R1 ',IOPT,TYPE,
     &                  IDUM, RDUM, 0, WORK(KEND0), LEND0)
      ELSE IF (TYPE(1:3).EQ.'EX1') THEN
        IOPT = 4
        CALL CC_FMATRIX(WORK(K1FTRAN),N1FTRAN,'LE ','R1 ',IOPT,TYPE,
     &                  IDUM, RDUM, 0, WORK(KEND0), LEND0)
      END IF

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) MSGDBG, 
     &        'NORM^2 of ETA vectors after F matrix terms:'
        DO IVEC = 1, NVEC
          ISYM = ILSTSYM(TYPE,IVEC)
          CALL CC_RDRSP(TYPE,IVEC,ISYM,IOPTRW,MODEL,
     &                  WORK(KCHI1),WORK(KCHI2))
          IF (CCR12) CALL CC_RDRSP(TYPE,IVEC,ISYM,IOPTRWR12,MODEL,
     &                             DUMMY,WORK(KCHIR12))
          XNORM = DDOT(NT1AM(ISYM),WORK(KCHI1),1,WORK(KCHI1),1)
          IF (.NOT.CCS) 
     &     XNORM = XNORM+DDOT(NT2AM(ISYM),WORK(KCHI2),1,WORK(KCHI2),1)
          IF (CCR12)
     &     XNORM = XNORM+DDOT(NTR12AM(ISYM),WORK(KCHIR12),1,
     &                        WORK(KCHIR12),1)
          WRITE (LUPRI,*) MSGDBG, IVEC,XNORM
        END DO
      END IF

*---------------------------------------------------------------------*
* calculate F{O} matrix contributions:
*---------------------------------------------------------------------*
      IF (TYPE(1:2).EQ.'X2') THEN
        CALL CCQR_FADRV('L0 ','o1 ','R1 ',TYPE,N0FATRAN, 0,
     &                   WORK(K0FATRAN),IDUMMY,IDUMMY,
     &                   WORK(KEND0), LEND0, 'FILE' )
      ELSE IF (TYPE(1:3).EQ.'CX2') THEN
        CALL CCQR_FADRV('L0 ','o1 ','RC ',TYPE,N0FATRAN, 0,
     &                   WORK(K0FATRAN),IDUMMY,IDUMMY,
     &                   WORK(KEND0), LEND0, 'FILE' )
      ELSE IF (TYPE(1:2).EQ.'X3') THEN
        CALL CCQR_FADRV('L0 ','o1 ','R2 ',TYPE,N0FATRAN, 0,
     &                   WORK(K0FATRAN),IDUMMY,IDUMMY,
     &                   WORK(KEND0), LEND0, 'FILE' )
        CALL CCQR_FADRV('L1 ','o1 ','R1 ',TYPE,N1FATRAN, 0,
     &                   WORK(K1FATRAN),IDUMMY,IDUMMY,
     &                   WORK(KEND0), LEND0, 'FILE' )
      ELSE IF (TYPE(1:3).EQ.'EX2') THEN
        CALL CCQR_FADRV('LE ','o1 ','R1 ',TYPE,N0FATRAN, 0,
     &                   WORK(K0FATRAN),IDUMMY,IDUMMY,
     &                   WORK(KEND0), LEND0, 'FILE' )
      END IF

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) MSGDBG,
     &        'NORM^2 of ETA vectors after F{O} matrix terms:'
        DO IVEC = 1, NVEC
          ISYM = ILSTSYM(TYPE,IVEC)
          CALL CC_RDRSP(TYPE,IVEC,ISYM,IOPTRW,MODEL,
     &                  WORK(KCHI1),WORK(KCHI2))
          IF (CCR12) CALL CC_RDRSP(TYPE,IVEC,ISYM,IOPTRWR12,MODEL,
     &                             DUMMY,WORK(KCHIR12))
          XNORM = DDOT(NT1AM(ISYM),WORK(KCHI1),1,WORK(KCHI1),1)
          IF (.NOT.CCS) 
     &     XNORM = XNORM+DDOT(NT2AM(ISYM),WORK(KCHI2),1,WORK(KCHI2),1)
          IF (CCR12)
     &     XNORM = XNORM+DDOT(NTR12AM(ISYM),WORK(KCHIR12),1,
     &                        WORK(KCHIR12),1)
          WRITE (LUPRI,*) MSGDBG, IVEC,XNORM
        END DO
      END IF

*---------------------------------------------------------------------*
* calculate ETA{O} vector contributions:
*---------------------------------------------------------------------*
      IF (TYPE(1:2).EQ.'X2') THEN
        IOPT  = 4
        ORDER = 1
        CALL CC_XIETA(WORK(KXETRAN),NEATRAN,IOPT, ORDER, 'L1 ',
     &                'O1 ', IDUM, RDUM, TYPE, IDUM, RDUM,
     &                .FALSE.,0, WORK(KEND0),LEND0)
      ELSE IF (TYPE(1:3).EQ.'CX2') THEN
        IOPT  = 4
        ORDER = 1
        CALL CC_XIETA(WORK(KXETRAN),NEATRAN,IOPT, ORDER, 'LC ',
     &                'O1 ', IDUM, RDUM, TYPE, IDUM, RDUM,
     &                .FALSE.,0, WORK(KEND0),LEND0)
      ELSE IF (TYPE(1:2).EQ.'X3') THEN
        IOPT  = 4
        ORDER = 1
        CALL CC_XIETA(WORK(KXETRAN),NEATRAN,IOPT, ORDER, 'L2 ',
     &                'O1 ', IDUM, RDUM, TYPE, IDUM, RDUM,
     &                .FALSE.,0, WORK(KEND0),LEND0)
      ELSE IF (TYPE(1:3).EQ.'EX2') THEN
        IOPT  = 4
        ORDER = 1
        CALL CC_XIETA(WORK(KXETRAN),NEATRAN,IOPT, ORDER, 'EL1',
     &                'O1 ', IDUM, RDUM, TYPE, IDUM, RDUM,
     &                .FALSE.,0, WORK(KEND0),LEND0)
      ELSE IF (TYPE(1:3).EQ.'EX1') THEN
        IOPT  = 4
        ORDER = 1
        CALL CC_XIETA(WORK(KXETRAN),NEATRAN,IOPT, ORDER, 'LE ',
     &                'O1 ', IDUM, RDUM, TYPE, IDUM, RDUM,
     &                .FALSE.,0, WORK(KEND0),LEND0)
      END IF

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) MSGDBG,
     &        'NORM^2 of ETA vectors after ETA{O} vec. terms:'
        DO IVEC = 1, NVEC
          ISYM = ILSTSYM(TYPE,IVEC)
          CALL CC_RDRSP(TYPE,IVEC,ISYM,IOPTRW,MODEL,
     &                  WORK(KCHI1),WORK(KCHI2))
          IF (CCR12) CALL CC_RDRSP(TYPE,IVEC,ISYM,IOPTRWR12,MODEL,
     &                             DUMMY,WORK(KCHIR12))
          XNORM = DDOT(NT1AM(ISYM),WORK(KCHI1),1,WORK(KCHI1),1)
          IF (.NOT.CCS) 
     &     XNORM = XNORM+DDOT(NT2AM(ISYM),WORK(KCHI2),1,WORK(KCHI2),1)
          IF (CCR12)
     &     XNORM = XNORM+DDOT(NTR12AM(ISYM),WORK(KCHIR12),1,
     &                        WORK(KCHIR12),1)
          WRITE (LUPRI,*) MSGDBG, IVEC,XNORM
        END DO
      END IF
*---------------------------------------------------------------------*
* test (static) EX1 vectors by calculating the excited state FOP's
*---------------------------------------------------------------------*
      IF (LOCDBG .AND. TYPE(1:3).EQ.'EX1') THEN
        KRHS1 = KEND1
        KRHS2 = KRHS1 + LMAX1
        KEND2 = KRHS2 + LMAX2
        IF (CCR12) THEN
          KRHSR12 = KEND2
          KEND2   = KRHSR12 + LMAXR12
        END IF
        LEND2 = LWORK - KEND2

        IF (LEND2 .LT. 0 ) THEN
          CALL QUIT('Insufficient work space in CC_ETADRV.')
        END IF

        WRITE (LUPRI,*) MSGDBG, 'excited state first-order properties:'
        DO IVEC = 1, NVEC
        IF (ISYMO(IVEC,1).EQ.1 .AND. FREQS(IVEC,1).EQ.ZERO) THEN
          ISYM = ILSTSYM(TYPE,IVEC)
          CALL CC_RDRSP(TYPE,IVEC,ISYM,IOPTRW,MODEL,
     &                  WORK(KCHI1),WORK(KCHI2))
          IF (CCR12) CALL CC_RDRSP(TYPE,IVEC,ISYM,IOPTRWR12,MODEL,
     &                             DUMMY,WORK(KCHIR12))
          CALL CC_RDRSP('RE',ISTAT(IVEC,1),ISYMS(IVEC,1),IOPTRW,MODEL,
     &                  WORK(KRHS1),WORK(KRHS2))
          IF (CCR12) CALL CC_RDRSP('RE',ISTAT(IVEC,1),ISYMS(IVEC,1),
     &                             IOPTRWR12,MODEL,DUMMY,WORK(KRHSR12))
          XNORM = DDOT(NT1AM(ISYM),WORK(KCHI1),1,WORK(KRHS1),1)
          IF (.NOT. CCS)
     &     XNORM = XNORM+DDOT(NT2AM(ISYM),WORK(KCHI2),1,WORK(KRHS2),1)
          IF (CCR12) 
     &     XNORM = XNORM+DDOT(NTR12AM(ISYM),WORK(KCHIR12),1,
     &                        WORK(KRHSR12),1)
          WRITE (LUPRI,'(A,I3,2X,F12.8,2X,A,2X,F12.8)') MSGDBG,
     &              ISTAT(IVEC,1),EIGV(IVEC,1),LABELV(IVEC,1),XNORM
        ELSE
          WRITE (LUPRI,'(A,I3,2X,F12.8,2X,A,2X,F12.8)') MSGDBG,
     &              ISTAT(IVEC,1),EIGV(IVEC,1),LABELV(IVEC,1),ZERO
        END IF
        END DO

      END IF
*---------------------------------------------------------------------*
* that's it:
*---------------------------------------------------------------------*

      RETURN
      END

*=====================================================================*
*              END OF SUBROUTINE CC_ETADRV                            *
*=====================================================================*

c /* deck CC_ETA_SETUP */
*=====================================================================*
      SUBROUTINE CC_ETA_SETUP(TYPE,NSTAT,ORDER,LAB,
     &                        ISTAT,EIGV,FREQ,ICAU,
     &                        NVEC, MAXVEC, MXTRAN, 
     &                        I0KTRAN, N0KTRAN,
     &                        I0GTRAN, N0GTRAN,
     &                        I1GTRAN, N1GTRAN,
     &                        I1FTRAN, N1FTRAN,
     &                        I2FTRAN, N2FTRAN,
     &                        I0FATRAN,N0FATRAN,
     &                        I1FATRAN,N1FATRAN,
     &                        IXETRAN, NEATRAN  )
*---------------------------------------------------------------------*
*
*    Purpose: set up for CC_ETA section
*                - list of G matrix transformations 
*                - list of F matrix transformations 
*                - list of F{O} matrix transformations 
*                - list of ETA{O} vector calculations 
*
*     Written by Christof Haettig, april/june/july 1997.
*     extensions for Cauchy eta vectors in march 1998.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <cclists.h>

* local parameters:
      CHARACTER*(22) MSGDBG
      PARAMETER (MSGDBG = '[debug] CC_ETA_SETUP> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)


      INTEGER MXORD, MXORD2, MXORD3, MXSTAT
      PARAMETER (MXORD  = 4, MXSTAT = 2)
      PARAMETER (MXORD2 = MXORD *(MXORD-1)/2 )
      PARAMETER (MXORD3 = MXORD2*(MXORD-2)/3 )


      INTEGER MXTRAN, NSTAT, ORDER, MAXVEC, NVEC

      CHARACTER*(*) TYPE

      CHARACTER*(8) LAB(MAXVEC,*)
      INTEGER ISTAT(MAXVEC,*), ICAU(MAXVEC,*)

#if defined (SYS_CRAY)
      REAL FREQ(MAXVEC,*), EIGV(MAXVEC,*)
#else
      DOUBLE PRECISION FREQ(MAXVEC,*), EIGV(MAXVEC,*)
#endif

      INTEGER I0KTRAN(5,MXTRAN)
      INTEGER I0GTRAN(4,MXTRAN)
      INTEGER I1GTRAN(4,MXTRAN)
      INTEGER I1FTRAN(3,MXTRAN)
      INTEGER I2FTRAN(3,MXTRAN)
      INTEGER I0FATRAN(5,MXTRAN)
      INTEGER I1FATRAN(5,MXTRAN)
      INTEGER IXETRAN(MXDIM_XEVEC,MXTRAN)

      INTEGER N0KTRAN, N1GTRAN, N2FTRAN, N1FATRAN
      INTEGER          N0GTRAN, N1FTRAN, N0FATRAN, NEATRAN

      INTEGER IVEC, ISYML, ITRAN, I, IDX, IDXA, IDXB, IDXAB, IDXS

      INTEGER A, B, C, D
      PARAMETER (A = 1, B = 2, C = 3, D = 4)
      INTEGER AB, AC, AD, BC, BD, CD
      PARAMETER (AB = 1, AC = 2, BC = 3, AD = 4, BD = 5, CD = 6)
      INTEGER ABC, ABD, ACD, BCD 
      PARAMETER (ABC = 1, ABD = 2, ACD = 3, BCD = 4) 

      INTEGER NS2A, NS3A, NP3AB, NP4AB, NT4ABC
      PARAMETER (NS2A = 2, NS3A = 3, NP3AB = 3, NP4AB = 6, NT4ABC = 4)

      INTEGER ISA(NS3A), ISB(NS3A), ISC(NS3A)
      INTEGER IPAB(NP4AB), IPC(NP4AB), IPD(NP4AB), IPCD(NP4AB)
      INTEGER ITABC(NT4ABC), ITD(NT4ABC)

      DATA ISA  / A, B, C/
      DATA ISB  / B, A, A/
      DATA ISC  / C, C, B/

      DATA IPAB / AB, AC, BC, AD, BD, CD /
      DATA IPC  / C,  B,  A,  B,  A,  A  /
      DATA IPD  / D,  D,  D,  C,  C,  B  /
      DATA IPCD / CD, BD, AD, BC, AC, AB /

      DATA ITABC / ABC, ABD, ACD, BCD /
      DATA ITD   / D,   C,   B,   A   /


      INTEGER IL0
      PARAMETER (IL0 = 0)  ! index for zeroth-order zeta vector
      INTEGER IL1(MXORD),  IR1(MXORD), IOP(MXORD), ISYM(MXORD)
      INTEGER IL2(MXORD2), IR2(MXORD2)
      INTEGER IE0(MXSTAT), IE1(MXORD,MXSTAT), ISYMS(MXSTAT)
      INTEGER LEN

      CHARACTER CLASS*(5)


* external functions:
      INTEGER IROPER
      INTEGER IR1TAMP
      INTEGER IR2TAMP
      INTEGER IL1ZETA
      INTEGER IL2ZETA
      INTEGER IEL1AMP
      INTEGER IEL2AMP
      INTEGER ILRCAMP
      INTEGER ILC1AMP

*---------------------------------------------------------------------*
* initializations:
*---------------------------------------------------------------------*
      N0KTRAN  = 0
      N0GTRAN  = 0
      N1GTRAN  = 0
      N1FTRAN  = 0
      N2FTRAN  = 0
      N0FATRAN = 0
      N1FATRAN = 0
      NEATRAN  = 0

*---------------------------------------------------------------------*
* start loop over all requested ETA-vectors:
*---------------------------------------------------------------------*
 
      DO IVEC = 1, NVEC

* eigenvectors that contribute:
        IF (NSTAT.EQ.1) THEN
          DO IDXS = 1, NSTAT
            IE0(IDXS) = ISTAT(IVEC,IDXS)
          END DO
        END IF

* operators:
        IF (ORDER.GE.1) THEN
          DO IDXA = 1, ORDER
            IOP(IDXA) = IROPER(LAB(IVEC,IDXA),ISYML)
          END DO
        END IF

* operators and first-order vectors that contribute:
        IF (TYPE(1:1).EQ.'X' .AND. ORDER.GT.1) THEN
          DO IDXA = 1, ORDER
            IL1(IDXA) = IL1ZETA(LAB(IVEC,IDXA),.FALSE.,
     &                          FREQ(IVEC,IDXA),ISYML)
            IR1(IDXA) = IR1TAMP(LAB(IVEC,IDXA),.FALSE.,
     &                          FREQ(IVEC,IDXA),ISYML)
          END DO
        END IF
        IF (TYPE(1:2).EQ.'CX' .AND. ORDER.GT.1) THEN
          DO IDXA = 1, ORDER
            IR1(IDXA) = ILRCAMP(LAB(IVEC,IDXA),ICAU(IVEC,IDXA),ISYML)
            IL1(IDXA) = ILC1AMP(LAB(IVEC,IDXA),ICAU(IVEC,IDXA),ISYML)
          END DO
        END IF
        IF (TYPE(1:2).EQ.'EX' .AND. ORDER.GE.1) THEN
          DO IDXA = 1, ORDER
            IR1(IDXA) = IR1TAMP(LAB(IVEC,IDXA),.FALSE.,
     &                          FREQ(IVEC,IDXA),ISYML)
          END DO
          IF (ORDER.GT.1) THEN
            IE1(IDXA,1) = 
     &           IEL1AMP(ISTAT(IVEC,1),EIGV(IVEC,1),ISYMS(1),
     &                   LAB(IVEC,IDXA),FREQ(IVEC,IDXA),ISYM(IDXA),
     &                   .FALSE.,.FALSE.)
          END IF
        END IF
* second-order vectors that contribute:
      IF (ORDER.GT.2 .OR. (ORDER.GE.2 .AND. NSTAT.GE.1)) THEN

        IDXAB  = 0
        DO IDXB = 2, ORDER
        DO IDXA = 1, IDXB-1
         IDXAB = IDXAB + 1
         IR2(IDXAB) =
     &       IR2TAMP(LAB(IVEC,IDXA),.FALSE.,FREQ(IVEC,IDXA),ISYM(IDXA),
     &               LAB(IVEC,IDXB),.FALSE.,FREQ(IVEC,IDXB),ISYM(IDXB))
        END DO
        END DO

       IF (TYPE(1:2).NE.'EX') THEN
        IDXAB  = 0
        DO IDXB = 2, ORDER
        DO IDXA = 1, IDXB-1
         IDXAB = IDXAB + 1
         IL2(IDXAB) =
     &           IL2ZETA(LAB(IVEC,IDXA),FREQ(IVEC,IDXA),ISYM(IDXA),
     &                   LAB(IVEC,IDXB),FREQ(IVEC,IDXB),ISYM(IDXB))
        END DO
        END DO
       END IF

      END IF


*---------------------------------------------------------------------*
* set up list of H matrix transformations
*---------------------------------------------------------------------*
        IF (TYPE(1:2).EQ.'X3') THEN
          N0KTRAN = N0KTRAN + 1
          I0KTRAN(1,N0KTRAN) = IL0
          I0KTRAN(2,N0KTRAN) = IR1(A)
          I0KTRAN(3,N0KTRAN) = IR1(B)
          I0KTRAN(4,N0KTRAN) = IR1(C)
          I0KTRAN(5,N0KTRAN) = IVEC
        END IF
*---------------------------------------------------------------------*
* set up list of G matrix transformations
*---------------------------------------------------------------------*
        IF (TYPE(1:2).EQ.'X2') THEN
          N0GTRAN = N0GTRAN + 1
          I0GTRAN(1,N0GTRAN) = IL0
          I0GTRAN(2,N0GTRAN) = IR1(A)
          I0GTRAN(3,N0GTRAN) = IR1(B)
          I0GTRAN(4,N0GTRAN) = IVEC
        ELSE IF (TYPE(1:3).EQ.'CX2') THEN
          N0GTRAN = N0GTRAN + 1
          I0GTRAN(1,N0GTRAN) = IL0
          I0GTRAN(2,N0GTRAN) = IR1(A)
          I0GTRAN(3,N0GTRAN) = IR1(B)
          I0GTRAN(4,N0GTRAN) = IVEC
        ELSE IF (TYPE(1:2).EQ.'X3') THEN
          DO IDX = 1, NP3AB
            N0GTRAN = N0GTRAN + 1
            I0GTRAN(1,N0GTRAN) = IL0
            I0GTRAN(2,N0GTRAN) = IR2(IPAB(IDX))
            I0GTRAN(3,N0GTRAN) = IR1(IPC(IDX))
            I0GTRAN(4,N0GTRAN) = IVEC
          END DO

          DO IDX = 1, NS3A
            N1GTRAN = N1GTRAN + 1
            I1GTRAN(1,N1GTRAN) = IL1(ISA(IDX))
            I1GTRAN(2,N1GTRAN) = IR1(ISB(IDX))
            I1GTRAN(3,N1GTRAN) = IR1(ISC(IDX))
            I1GTRAN(4,N1GTRAN) = IVEC
          END DO
        ELSE IF (TYPE(1:3).EQ.'EX2') THEN
          N0GTRAN = N0GTRAN + 1
          I0GTRAN(1,N0GTRAN) = IE0(1)
          I0GTRAN(2,N0GTRAN) = IR1(A)
          I0GTRAN(3,N0GTRAN) = IR1(B)
          I0GTRAN(4,N0GTRAN) = IVEC
        END IF

*---------------------------------------------------------------------*
* set up list of F matrix transformations
*---------------------------------------------------------------------*
        IF (TYPE(1:2).EQ.'X2') THEN
          N1FTRAN = N1FTRAN + 1
          I1FTRAN(1,N1FTRAN) = IL1(A)
          I1FTRAN(2,N1FTRAN) = IR1(B)
          I1FTRAN(3,N1FTRAN) = IVEC

          N1FTRAN = N1FTRAN + 1
          I1FTRAN(1,N1FTRAN) = IL1(B)
          I1FTRAN(2,N1FTRAN) = IR1(A)
          I1FTRAN(3,N1FTRAN) = IVEC
        ELSE IF (TYPE(1:3).EQ.'CX2') THEN
          N1FTRAN = N1FTRAN + 1
          I1FTRAN(1,N1FTRAN) = IL1(A)
          I1FTRAN(2,N1FTRAN) = IR1(B)
          I1FTRAN(3,N1FTRAN) = IVEC

          N1FTRAN = N1FTRAN + 1
          I1FTRAN(1,N1FTRAN) = IL1(B)
          I1FTRAN(2,N1FTRAN) = IR1(A)
          I1FTRAN(3,N1FTRAN) = IVEC
        ELSE IF (TYPE(1:2).EQ.'X3') THEN
          DO IDX = 1, NP3AB
            N1FTRAN = N1FTRAN + 1
            I1FTRAN(1,N1FTRAN) = IL1(IPC(IDX))
            I1FTRAN(2,N1FTRAN) = IR2(IPAB(IDX))
            I1FTRAN(3,N1FTRAN) = IVEC
          END DO

          DO IDX = 1, NP3AB
            N2FTRAN = N2FTRAN + 1
            I2FTRAN(1,N2FTRAN) = IL2(IPAB(IDX))
            I2FTRAN(2,N2FTRAN) = IR1(IPC(IDX))
            I2FTRAN(3,N2FTRAN) = IVEC
          END DO
        ELSE IF (TYPE(1:3).EQ.'EX2') THEN
          N1FTRAN = N1FTRAN + 1
          I1FTRAN(1,N1FTRAN) = IE0(1)
          I1FTRAN(2,N1FTRAN) = IR2(AB)
          I1FTRAN(3,N1FTRAN) = IVEC

          N2FTRAN = N2FTRAN + 1
          I2FTRAN(1,N2FTRAN) = IE1(A,1)
          I2FTRAN(2,N2FTRAN) = IR1(B)
          I2FTRAN(3,N2FTRAN) = IVEC

          N2FTRAN = N2FTRAN + 1
          I2FTRAN(1,N2FTRAN) = IE1(B,1)
          I2FTRAN(2,N2FTRAN) = IR1(A)
          I2FTRAN(3,N2FTRAN) = IVEC
        ELSE IF (TYPE(1:3).EQ.'EX1') THEN
          N1FTRAN = N1FTRAN + 1
          I1FTRAN(1,N1FTRAN) = IE0(1)
          I1FTRAN(2,N1FTRAN) = IR1(A)
          I1FTRAN(3,N1FTRAN) = IVEC
        END IF

*---------------------------------------------------------------------*
* set up list of F{O} matrix transformations
*---------------------------------------------------------------------*
        IF (TYPE(1:2).EQ.'X2') THEN
          N0FATRAN = N0FATRAN + 1
          I0FATRAN(1,N0FATRAN) = IL0
          I0FATRAN(2,N0FATRAN) = IOP(A)
          I0FATRAN(3,N0FATRAN) = IR1(B)
          I0FATRAN(4,N0FATRAN) = IVEC
          I0FATRAN(5,N0FATRAN) = 0

          N0FATRAN = N0FATRAN + 1
          I0FATRAN(1,N0FATRAN) = IL0
          I0FATRAN(2,N0FATRAN) = IOP(B)
          I0FATRAN(3,N0FATRAN) = IR1(A)
          I0FATRAN(4,N0FATRAN) = IVEC
          I0FATRAN(5,N0FATRAN) = 0
        ELSE IF (TYPE(1:3).EQ.'CX2') THEN
          IF (ICAU(IVEC,A).EQ.0) THEN
            N0FATRAN = N0FATRAN + 1
            I0FATRAN(1,N0FATRAN) = IL0
            I0FATRAN(2,N0FATRAN) = IOP(A)
            I0FATRAN(3,N0FATRAN) = IR1(B)
            I0FATRAN(4,N0FATRAN) = IVEC
            I0FATRAN(5,N0FATRAN) = 0
          END IF

          IF (ICAU(IVEC,B).EQ.0) THEN
            N0FATRAN = N0FATRAN + 1
            I0FATRAN(1,N0FATRAN) = IL0
            I0FATRAN(2,N0FATRAN) = IOP(B)
            I0FATRAN(3,N0FATRAN) = IR1(A)
            I0FATRAN(4,N0FATRAN) = IVEC
            I0FATRAN(5,N0FATRAN) = 0
          END IF
        ELSE IF (TYPE(1:2).EQ.'X3') THEN
          DO IDX = 1, NP3AB
            N0FATRAN = N0FATRAN + 1
            I0FATRAN(1,N0FATRAN) = IL0
            I0FATRAN(2,N0FATRAN) = IOP(IPC(IDX))
            I0FATRAN(3,N0FATRAN) = IR2(IPAB(IDX))
            I0FATRAN(4,N0FATRAN) = IVEC
            I0FATRAN(5,N0FATRAN) = 0
          END DO

          DO IDX = 1, NP3AB
            N1FATRAN = N1FATRAN + 1
            I1FATRAN(1,N1FATRAN) = IL1(ISA(IDX))
            I1FATRAN(2,N1FATRAN) = IOP(ISB(IDX))
            I1FATRAN(3,N1FATRAN) = IR1(ISC(IDX))
            I1FATRAN(4,N1FATRAN) = IVEC
            I1FATRAN(5,N1FATRAN) = 0
            N1FATRAN = N1FATRAN + 1
            I1FATRAN(1,N1FATRAN) = IL1(ISA(IDX))
            I1FATRAN(2,N1FATRAN) = IOP(ISC(IDX))
            I1FATRAN(3,N1FATRAN) = IR1(ISB(IDX))
            I1FATRAN(4,N1FATRAN) = IVEC
            I1FATRAN(5,N1FATRAN) = 0
          END DO
        ELSE IF (TYPE(1:3).EQ.'EX2') THEN
          N0FATRAN = N0FATRAN + 1
          I0FATRAN(1,N0FATRAN) = IE0(1)
          I0FATRAN(2,N0FATRAN) = IOP(A)
          I0FATRAN(3,N0FATRAN) = IR1(B)
          I0FATRAN(4,N0FATRAN) = IVEC
          I0FATRAN(5,N0FATRAN) = 0

          N0FATRAN = N0FATRAN + 1
          I0FATRAN(1,N0FATRAN) = IE0(1)
          I0FATRAN(2,N0FATRAN) = IOP(B)
          I0FATRAN(3,N0FATRAN) = IR1(A)
          I0FATRAN(4,N0FATRAN) = IVEC
          I0FATRAN(5,N0FATRAN) = 0
        END IF

*---------------------------------------------------------------------*
* set up list of ETA{O} vector calculations:
*---------------------------------------------------------------------*
        IF (TYPE(1:2).EQ.'X2') THEN
          NEATRAN = NEATRAN + 1

          IXETRAN(1,NEATRAN) = IOP(B)
          IXETRAN(2,NEATRAN) = IL1(A)
          IXETRAN(3,NEATRAN) = -1
          IXETRAN(4,NEATRAN) = IVEC
          IXETRAN(5,NEATRAN) = 0
          IXETRAN(6,NEATRAN) = 0
          IXETRAN(7,NEATRAN) = 0
          IXETRAN(8,NEATRAN) = 0
  
          NEATRAN = NEATRAN + 1

          IXETRAN(1,NEATRAN) = IOP(A)
          IXETRAN(2,NEATRAN) = IL1(B)
          IXETRAN(3,NEATRAN) = -1
          IXETRAN(4,NEATRAN) = IVEC
          IXETRAN(5,NEATRAN) = 0
          IXETRAN(6,NEATRAN) = 0
          IXETRAN(7,NEATRAN) = 0
          IXETRAN(8,NEATRAN) = 0
        ELSE IF (TYPE(1:3).EQ.'CX2') THEN
          IF (ICAU(IVEC,B).EQ.0) THEN
            NEATRAN = NEATRAN + 1

            IXETRAN(1,NEATRAN) = IOP(B)
            IXETRAN(2,NEATRAN) = IL1(A)
            IXETRAN(3,NEATRAN) = -1
            IXETRAN(4,NEATRAN) = IVEC
            IXETRAN(5,NEATRAN) = 0
            IXETRAN(6,NEATRAN) = 0
            IXETRAN(7,NEATRAN) = 0
            IXETRAN(8,NEATRAN) = 0
          END IF
  
          IF (ICAU(IVEC,A).EQ.0) THEN
            NEATRAN = NEATRAN + 1

            IXETRAN(1,NEATRAN) = IOP(A)
            IXETRAN(2,NEATRAN) = IL1(B)
            IXETRAN(3,NEATRAN) = -1
            IXETRAN(4,NEATRAN) = IVEC
            IXETRAN(5,NEATRAN) = 0
            IXETRAN(6,NEATRAN) = 0
            IXETRAN(7,NEATRAN) = 0
            IXETRAN(8,NEATRAN) = 0
          END IF
        ELSE IF (TYPE(1:2).EQ.'X3') THEN
          DO IDX = 1, NP3AB
            NEATRAN = NEATRAN + 1

            IXETRAN(1,NEATRAN) = IOP(IPC(IDX))
            IXETRAN(2,NEATRAN) = IL2(IPAB(IDX))
            IXETRAN(3,NEATRAN) = -1
            IXETRAN(4,NEATRAN) = IVEC
            IXETRAN(5,NEATRAN) = 0
            IXETRAN(6,NEATRAN) = 0
            IXETRAN(7,NEATRAN) = 0
            IXETRAN(8,NEATRAN) = 0
          END DO
        ELSE IF (TYPE(1:3).EQ.'EX2') THEN
          NEATRAN = NEATRAN + 1

          IXETRAN(1,NEATRAN) = IOP(B)
          IXETRAN(2,NEATRAN) = IE1(A,1)
          IXETRAN(3,NEATRAN) = -1
          IXETRAN(4,NEATRAN) = IVEC
          IXETRAN(5,NEATRAN) = 0
          IXETRAN(6,NEATRAN) = 0
          IXETRAN(7,NEATRAN) = 0
          IXETRAN(8,NEATRAN) = 0
  
          NEATRAN = NEATRAN + 1

          IXETRAN(1,NEATRAN) = IOP(A)
          IXETRAN(2,NEATRAN) = IE1(B,1)
          IXETRAN(3,NEATRAN) = -1
          IXETRAN(4,NEATRAN) = IVEC
          IXETRAN(5,NEATRAN) = 0
          IXETRAN(6,NEATRAN) = 0
          IXETRAN(7,NEATRAN) = 0
          IXETRAN(8,NEATRAN) = 0
        ELSE IF (TYPE(1:3).EQ.'EX1') THEN
          NEATRAN = NEATRAN + 1

          IXETRAN(1,NEATRAN) = IOP(A)
          IXETRAN(2,NEATRAN) = IE0(1)
          IXETRAN(3,NEATRAN) = -1
          IXETRAN(4,NEATRAN) = IVEC
          IXETRAN(5,NEATRAN) = 0
          IXETRAN(6,NEATRAN) = 0
          IXETRAN(7,NEATRAN) = 0
          IXETRAN(8,NEATRAN) = 0
        END IF
  
*---------------------------------------------------------------------*
* end loop over all requested ETA vectors
*---------------------------------------------------------------------*
      END DO

*---------------------------------------------------------------------*
* print the lists: 
*---------------------------------------------------------------------*
* general statistics:
      IF (TYPE(1:1).EQ.'X') THEN 
          LEN = 2
          CLASS = ' eta '
      ELSE IF (TYPE(1:2).EQ.'CX') THEN
          LEN = 3
          CLASS = ' eta '
      ELSE IF (TYPE(1:2).EQ.'EX') THEN
          LEN = 3
          CLASS = ' rhs '
      ELSE
          LEN = 2
          CLASS = '     '
      END IF
      WRITE (LUPRI,'(/,/3X,A,I3,1X,3A)') 'For the requested',NVEC,
     &      TYPE(1:LEN),CLASS,' vectors'
      WRITE (LUPRI,'((8X,A,I3,A))') 
     &   ' - ',N0KTRAN,  ' H matrix transformations ',
     &   ' - ',N0GTRAN,  ' G matrix transformations ',
     &   ' - ',N1GTRAN,  ' generalized G matrix transformations ',
     &   ' - ',(N1FTRAN+N2FTRAN), 
     &                   ' generalized F matrix transformations ',
     &   ' - ',N0FATRAN, ' F{O} matrix transformations ',
     &   ' - ',N1FATRAN, ' generalized F{O} matrix transformations ',
     &   ' - ',NEATRAN,  ' generalized ETA{O} vector calculations ' 
      WRITE (LUPRI,'(3X,A,/,/)') 'will be performed.'


      IF (LOCDBG) THEN

* H matrix transformations:
        IF (N0KTRAN.GT.0) WRITE(LUPRI,*)
     &        'List of H matrix transformations:'
        DO ITRAN = 1, N0KTRAN
          WRITE(LUPRI,'(A,4I5,5X,(12I5,20X))') MSGDBG,
     &     (I0KTRAN(I,ITRAN),I=1,4)
        END DO
        WRITE (LUPRI,*)

* G matrix transformations:
        IF (N0GTRAN.GT.0) WRITE(LUPRI,*)
     &       'List of G matrix transformations:'
        DO ITRAN = 1, N0GTRAN
          WRITE(LUPRI,'(A,4I5,5X,(12I5,20X))') MSGDBG,
     &     (I0GTRAN(I,ITRAN),I=1,4)
        END DO
        WRITE (LUPRI,*)

        IF (N1GTRAN.GT.0) 
     &      WRITE (LUPRI,*) 'List of (T^1 C) matrix transformations:'
        DO ITRAN = 1, N1GTRAN
          WRITE(LUPRI,'(A,4I5,5X,(12I5,20X))') MSGDBG,
     &     (I1GTRAN(I,ITRAN),I=1,4)
        END DO
        WRITE (LUPRI,*)

* F matrix transformations:
        IF (N1FTRAN.GT.0)
     &      WRITE (LUPRI,*) 'List of (T^1 B) matrix transformations:'
        DO ITRAN = 1, N1FTRAN
          WRITE(LUPRI,'(A,3I5,5X,(12I5,20X))') MSGDBG,
     &     (I1FTRAN(I,ITRAN),I=1,3)
        END DO
        WRITE (LUPRI,*)

        IF (N2FTRAN.GT.0)
     &    WRITE (LUPRI,*) 'List of (T^2 B) matrix transformations:'
        DO ITRAN = 1, N2FTRAN
          WRITE(LUPRI,'(A,3I5,5X,(12I5,20X))') MSGDBG,
     &     (I2FTRAN(I,ITRAN),I=1,3)
        END DO
        WRITE (LUPRI,*)

* F{O} matrix transformations:
        IF (N0FATRAN.GT.0)
     &     WRITE (LUPRI,*) 'List of F{O} matrix transformations:'
        DO ITRAN = 1, N0FATRAN
          WRITE(LUPRI,'(A,4I5,5X,(12I5,20X))') MSGDBG,
     &     (I0FATRAN(I,ITRAN),I=1,4)
        END DO
        WRITE (LUPRI,*)

        IF (N1FATRAN.GT.0)
     &      WRITE (LUPRI,*) 'List of (T^1 B{O}) matrix transformations:'
        DO ITRAN = 1, N1FATRAN
          WRITE(LUPRI,'(A,4I5,5X,(12I5,20X))') MSGDBG,
     &     (I1FATRAN(I,ITRAN),I=1,4)
        END DO
        WRITE (LUPRI,*)

* ETA{O} vector calculations:
        IF (NEATRAN.GT.0)
     &    WRITE (LUPRI,*) 'List of (T^n A{O}) matrix transformations:'
        DO ITRAN = 1, NEATRAN
          WRITE(LUPRI,'(A,3I5,5X,(12I5,20X))') MSGDBG,
     &     (IXETRAN(I,ITRAN),I=1,4)
        END DO
        WRITE (LUPRI,*)

      END IF

      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CC_ETA_SETUP                         *
*---------------------------------------------------------------------*
