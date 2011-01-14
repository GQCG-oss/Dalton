C*******************************************************************************
C
C This file contains the core routines of the
C
C      Cholesky Decomposition Utility, Version 2.0
C
C written 2004/2005 by
C
C      Thomas Bondo Pedersen
C      Department of Theoretical Chemistry
C      Chemical Center, University of Lund
C      P.O. Box 124, S-221 00 Lund, Sweden.
C
C      Alfredo M. J. Sanchez de Meras
C      Institute of Molecular Science, Department of Physical Chemistry
C      University of Valencia
C      E-46100 Burjassot, Valencia, Spain.
C
C      Henrik Koch
C      Department of Physical Chemistry
C      Norwegian University of Science and Technology
C      N-7491 Trondheim, Norway.
C
C*******************************************************************************
C  /* Deck cho_drv */
      SUBROUTINE CHO_DRV(WORK,LWORK)
#include <implicit.h>
      DIMENSION WORK(LWORK)
#include <dummy.h>
C
C     Purpose: driver for the Cholesky decomposition of two-electron
C              integrals.
C
#include <cholesky.h>

      INTEGER ISYLST(8)

      LOGICAL LOCDBG, LCONV
      PARAMETER (LOCDBG = .FALSE.)

      CHARACTER*7 SECNAM
      PARAMETER (SECNAM = 'CHO_DRV')

C     Start overall timing.
C     ---------------------

      CALL CHO_TIMER(TCPU0,TWALL0)
      IF (LOCDBG) THEN
         WRITE(LUPRI,*) 'Entering ',SECNAM,': LWORK = ',LWORK
         CALL CHO_FLUSH(LUPRI)
      END IF

C     INITIALIZATION.
C     ===============

      ISEC = 1
      CALL CHO_TIMER(TIMSEC(1,ISEC),TIMSEC(3,ISEC))

C     Set some memory info (only strictly necessary for DALTON).
C     ----------------------------------------------------------

      KFREE = 1
      LFREE = LWORK

C     Initializations (incl. initialization of integral program).
C     Note that the initialization does some allocations, returning
C     pointers in KINFRD and KINFVC.
C     -------------------------------------------------------------

      KINFRD = -99999999
      KINFVC = KINFRD
      CALL CHO_INIT(WORK,LWORK,KINFRD,KINFVC,KFREE,LFREE)

C     Allocate memory for reduced set index arrays.
C     ---------------------------------------------

      NNSHL  = NSHELL*(NSHELL + 1)/2
      LRINDX = NSYM*NNSHL*3
      CALL CHO_MEM('iibstrsh','ALLO','INTE',KIIBST,LRINDX,
     &             KFREE,LFREE)
      CALL CHO_MEM('nnbstrsh','ALLO','INTE',KNNBST,LRINDX,
     &             KFREE,LFREE)
      CALL CHO_MEM('intmap','ALLO','INTE',KINMAP,NNSHL,
     &             KFREE,LFREE)

C     Initialize integral SP counter.
C     -------------------------------

      CALL CHO_INIMAP(WORK(KINMAP),NNSHL)

      CALL CHO_TIMER(TIMSEC(2,ISEC),TIMSEC(4,ISEC))
      CALL CHO_PRTTIM('Cholesky initialization',
     &                TIMSEC(2,ISEC),TIMSEC(1,ISEC),
     &                TIMSEC(4,ISEC),TIMSEC(3,ISEC),
     &                1)

C     GET DIAGONAL.
C     =============

      ISEC = 2
      CALL CHO_TIMER(TIMSEC(1,ISEC),TIMSEC(3,ISEC))
      IF (RSTDIA) THEN

C        Restart: first read index arrays.
C        ---------------------------------

         CALL CHO_RSTD_GETIND1(WORK(KIIBST),WORK(KNNBST),NSYM,NNSHL)

C        Allocate mapping arrays between reduced sets.
C        ---------------------------------------------

         LINDR = NNBSTRT(1)*3
         LINDS = NNBSTRT(1)
         CALL CHO_MEM('indred','ALLO','INTE',KINDR,LINDR,
     &                KFREE,LFREE)
         CALL CHO_MEM('indrsh','ALLO','INTE',KINDS,LINDS,
     &                KFREE,LFREE)

C        Read mapping arrays.
C        --------------------

         CALL CHO_RSTD_GETIND2(WORK(KINDR),WORK(KINDS),NNBSTRT(1))

C        Allocation: diagonal.
C        ---------------------

         NEEDR = 1
         NEEDI = 4*NEEDR

         CALL CHO_MEM('diarst','ALLO','REAL',KDIAG,NNBSTRT(1),
     &                KFREE,LFREE)
         CALL CHO_MEM('buf.2','ALLO','REAL',KBUF,NEEDR,
     &                KFREE,LFREE)
         CALL CHO_MEM('ibuf.2','ALLO','INTE',KIBUF,NEEDI,
     &                KFREE,LFREE)

         KREL = KBUF

C        Read diagonal.
C        --------------

         CALL CHO_GETDIAG1(WORK(KIIBST),WORK(KNNBST),WORK(KINDR),
     &                     WORK(KINDS),
     &                     WORK(KIBUF),WORK(KBUF),WORK(KDIAG),
     &                     NSYM,NNSHL,NNBSTRT(1),NEEDR,NDUMP)

         CALL CHO_MEM('buf.2','FLUSH','REAL',KREL,NEEDR,
     &                KFREE,LFREE)

      ELSE

C        Calculate diagonal and get 1st reduced set.
C        -------------------------------------------

         LSCR  = MX2SH
         NEEDR = LBUF + LSCR
         NEEDI = 4*LBUF

         CALL CHO_MEM('buf','ALLO','REAL',KREL,NEEDR,
     &                KFREE,LFREE)
         CALL CHO_MEM('ibuf','ALLO','INTE',KIBUF,NEEDI,
     &                KFREE,LFREE)

         KBUF  = KREL
         KSCR  = KBUF + LBUF

         NDUMP = 0

         CALL CHO_CALCDIAG(WORK(KBUF),WORK(KIBUF),LBUF,WORK(KSCR),LSCR,
     &                     WORK(KFREE),LFREE,WORK(KIIBST),WORK(KNNBST),
     &                     NSYM,NNSHL,NDUMP)

         CALL CHO_MEM('buf','FLUSH','REAL',KREL,NEEDR,
     &                KFREE,LFREE)

C        Allocate diagonal and mapping array between reduced sets.
C        Reallocate buffer.
C        ---------------------------------------------------------

         LINDR = NNBSTRT(1)*3
         LINDS = NNBSTRT(1)
         NEEDR = LBUF
         NEEDI = 4*LBUF

         CALL CHO_MEM('map','ALLO','INTE',KINDR,LINDR,
     &                KFREE,LFREE)
         CALL CHO_MEM('mapsh','ALLO','INTE',KINDS,LINDS,
     &                KFREE,LFREE)
         CALL CHO_MEM('dia','ALLO','REAL',KDIAG,NNBSTRT(1),
     &                KFREE,LFREE)
         CALL CHO_MEM('buf.2','ALLO','REAL',KBUF,NEEDR,
     &                KFREE,LFREE)
         CALL CHO_MEM('ibuf.2','ALLO','INTE',KIBUF,NEEDI,
     &                KFREE,LFREE)

         KREL = KBUF

C        Get diagonal in first reduced set.
C        ----------------------------------

         CALL CHO_GETDIAG1(WORK(KIIBST),WORK(KNNBST),WORK(KINDR),
     &                     WORK(KINDS),
     &                     WORK(KIBUF),WORK(KBUF),WORK(KDIAG),
     &                     NSYM,NNSHL,NNBSTRT(1),LBUF,NDUMP)

C        Write diagonal to disk.
C        -----------------------

         IOPT = 1
         CALL CHO_IODIAG(WORK(KDIAG),IOPT)

C        Deallocate back to and including buffer.
C        ----------------------------------------

         CALL CHO_MEM('buf','FLUSH','REAL',KREL,NEEDR,
     &                KFREE,LFREE)

      END IF

C     Allocate memory for iscratch array for reading vectors.
C     -------------------------------------------------------

      IF (CHO_IOVEC .EQ. 1) THEN
         LISCR = NNBSTR(1,1)
         DO ISYM = 2,NSYM
            LISCR = MAX(LISCR,NNBSTR(ISYM,1))
         END DO
      ELSE
         LISCR = 1
      END IF
      CALL CHO_MEM('ISCR','ALLO','INTE',KISCR,LISCR,KFREE,LFREE)

C     Update diagonal if restart, else just do analysis.
C     --------------------------------------------------

      LCONV = .FALSE.
      IF (RSTCHO) THEN
         CALL CHO_MEM('Cho.Rs1','MAX ','REAL',KWRK,LWRK,KFREE,LFREE)
         IF (LOCDBG) THEN
            WRITE(LUPRI,*) SECNAM,': restart diagonal:'
            DO ISYM = 1,NSYM
               ISYLST(ISYM) = ISYM
            END DO
            CALL CHO_PRTDIA(WORK(KDIAG),WORK(KIIBST),WORK(KNNBST),
     &                      WORK(KINDR),WORK(KINDS),ISYLST,
     &                      NSYM,NNSHL,NNBSTRT(1),NSYM,1)
         END IF
         CALL CHO_RESTART2(WORK(KIIBST),WORK(KNNBST),WORK(KINDR),
     &                     WORK(KINDS),WORK(KINFRD),WORK(KINFVC),
     &                     WORK(KISCR),WORK(KDIAG),WORK(KWRK),LWRK,
     &                     NSYM,NNSHL,NNBSTRT(1),MAXRED,MAXVEC,LISCR,
     &                     .FALSE.,LCONV)
         CALL CHO_MEM('Cho.Rs1','FREE','REAL',KWRK,LWRK,KFREE,LFREE)
         IPRTRED = 2  ! print flag for cho_prtred
      ELSE
         BIN1  = 1.0D2
         STEP  = 1.0D-1
         NBIN  = 18
         CALL CHO_ANADIA(WORK(KDIAG),BIN1,STEP,NBIN,.TRUE.)
         IPRTRED = 1  ! print flag for cho_prtred
      END IF

      CALL CHO_PRTRED(WORK(KNNBST),NSYM,NNSHL,IPRTRED)

      CALL CHO_TIMER(TIMSEC(2,ISEC),TIMSEC(4,ISEC))
      CALL CHO_PRTTIM('Cholesky diagonal setup',
     &                TIMSEC(2,ISEC),TIMSEC(1,ISEC),
     &                TIMSEC(4,ISEC),TIMSEC(3,ISEC),
     &                1)

C     DECOMPOSITION.
C     ==============

      ISEC = 3
      IF (LCONV) THEN

         IF (RSTCHO) THEN
            WRITE(LUPRI,'(///,10X,A,A,A,//)')
     &      '***** ',SECNAM,': restarted calculation converged. *****'
            CALL CHO_DZERO(TIMSEC(1,ISEC),4)
         ELSE
            WRITE(LUPRI,'(A,A)')
     &      SECNAM,': logical error: converged but not restart?!?!'
            CALL CHO_QUIT('Error in '//SECNAM,103)
         END IF

      ELSE

         CALL CHO_TIMER(TIMSEC(1,ISEC),TIMSEC(3,ISEC))

         KFR = 1
         LFR = LFREE
         CALL CHO_SETADDR(WORK(KINFRD),WORK(KINFVC),
     &                    WORK(KIIBST),WORK(KNNBST),
     &                    WORK(KINDR),WORK(KINDS),
     &                    MAXRED,MAXVEC,NSYM,NNSHL,NNBSTRT(1))
         CALL CHO_DECDRV(WORK(KINFRD),WORK(KINFVC),
     &                   WORK(KIIBST),WORK(KNNBST),WORK(KISCR),
     &                   WORK(KINDR),WORK(KINDS),WORK(KINMAP),
     &                   WORK(KDIAG),WORK(KFREE),LFREE,
     &                   MAXRED,MAXVEC,NSYM,NNSHL,NNBSTRT(1),
     &                   LISCR,KFR,LFR)

         CALL CHO_TIMER(TIMSEC(2,ISEC),TIMSEC(4,ISEC))
         CALL CHO_PRTTIM('Cholesky decomposition',
     &                   TIMSEC(2,ISEC),TIMSEC(1,ISEC),
     &                   TIMSEC(4,ISEC),TIMSEC(3,ISEC),
     &                   1)

      END IF

C     CHECK DIAGONAL.
C     ===============

      ISEC = 4
      IF (LCONV) THEN

         CALL CHO_DZERO(TIMSEC(1,ISEC),4)

      ELSE

         CALL CHO_TIMER(TIMSEC(1,ISEC),TIMSEC(3,ISEC))

         CALL CHO_MEM('Cho.Rst','MAX ','REAL',KWRK,LWRK,KFREE,LFREE)
         CALL CHO_RESTART2(WORK(KIIBST),WORK(KNNBST),WORK(KINDR),
     &                     WORK(KINDS),WORK(KINFRD),WORK(KINFVC),
     &                     WORK(KISCR),WORK(KDIAG),WORK(KWRK),LWRK,
     &                     NSYM,NNSHL,NNBSTRT(1),MAXRED,MAXVEC,LISCR,
     &                     .TRUE.,LCONV)
         CALL CHO_MEM('Cho.Rst','FREE','REAL',KWRK,LWRK,KFREE,LFREE)

         IF (.NOT. LCONV) THEN
            WRITE(LUPRI,'(A,A)')
     &      SECNAM,': decomposition failed!'
            CALL CHO_QUIT('Decomposition failed',104)
         END IF

         CALL CHO_TIMER(TIMSEC(2,ISEC),TIMSEC(4,ISEC))
         CALL CHO_PRTTIM('Cholesky diagonal check',
     &                   TIMSEC(2,ISEC),TIMSEC(1,ISEC),
     &                   TIMSEC(4,ISEC),TIMSEC(3,ISEC),
     &                   1)

      END IF

C     CHECK INTEGRALS.
C     ================

      ISEC = 5
      IF (CHO_INTCHK) THEN

         CALL CHO_TIMER(TIMSEC(1,ISEC),TIMSEC(3,ISEC))
         CALL CHO_DBGINT(WORK(KINFRD),WORK(KINFVC),
     &                   WORK(KIIBST),WORK(KNNBST),WORK(KISCR),
     &                   WORK(KINDR),WORK(KINDS),
     &                   MAXRED,MAXVEC,NSYM,NNSHL,NNBSTRT(1),LISCR)
         CALL CHO_TIMER(TIMSEC(2,ISEC),TIMSEC(4,ISEC))
         CALL CHO_PRTTIM('Cholesky integral check',
     &                   TIMSEC(2,ISEC),TIMSEC(1,ISEC),
     &                   TIMSEC(4,ISEC),TIMSEC(3,ISEC),
     &                   1)

      ELSE

         CALL CHO_DZERO(TIMSEC(1,ISEC),4)

      END IF

C     STATISTICS.
C     ===========

      ISEC = 6
      CALL CHO_TIMER(TIMSEC(1,ISEC),TIMSEC(3,ISEC))
      CALL CHO_STAT(WORK(KINFRD),WORK(KINFVC),
     &              WORK(KIIBST),WORK(KNNBST),
     &              WORK(KINDR),WORK(KINDS),WORK(KINMAP),
     &              MAXRED,MAXVEC,NSYM,NNSHL,NNBSTRT(1))
      CALL CHO_TIMER(TIMSEC(2,ISEC),TIMSEC(4,ISEC))
      CALL CHO_PRTTIM('Cholesky statistics',
     &                TIMSEC(2,ISEC),TIMSEC(1,ISEC),
     &                TIMSEC(4,ISEC),TIMSEC(3,ISEC),
     &                1)

C     Finalize decomposition.
C     -----------------------

      CALL CHO_FINAL()

C     Close restart file.
C     -------------------

      CALL CHO_OPENVR(2)

C     Print total timing.
C     -------------------

      CALL CHO_TIMER(TCPU1,TWALL1)
      CALL CHO_PRTTIM('Cholesky procedure',TCPU1,TCPU0,
     &                TWALL1,TWALL0,1)

      END
C  /* Deck cho_init */
      SUBROUTINE CHO_INIT(WORK,LWORK,KINFRD,KINFVC,KFREE,LFREE)
#include <implicit.h>
      DIMENSION WORK(LWORK)
C
C     Purpose: initializations.
C
#include <choorb.h>
#include <cholesky.h>
#include <chov2io.h>

      DIMENSION XXB(8)

      CHARACTER*1  LINE
      CHARACTER*8  SECNAM
      CHARACTER*17 STRING
      PARAMETER (LINE = '=', SECNAM = 'CHO_INIT',
     &           STRING = 'Information from ')

      PARAMETER (GBLIM = 2.147483648D9)

C     Integral program specific initializations
C     =========================================

      CALL CHO_DAL_INIT(WORK,LWORK)

C     Set screening damping.
C     ----------------------

      CALL CHO_SETDAMP()

C     Initialize timings etc.
C     -----------------------

      CALL CHO_DZERO(TDECDRV,2)
      CALL CHO_DZERO(TINTEG,2*NINTEG)
      CALL CHO_DZERO(TDECOM,2*NDECOM)
      CALL CHO_DZERO(VCSTOR,NSYM)
      CALL CHO_IZERO(ICHKQ,4*(NCHKQ+1))
      CALL CHO_IZERO(NVECRS1,NSYM)
      CALL CHO_IZERO(CHOV2I,3)

      DID_DECDRV = .FALSE.

      DIAMNZ = 0.0D0
      IABMNZ = 0
      NNZTOT = 0

C     Allocate memory for INFRED and INFVEC arrays.
C     In so doing, determine the max. #vectors and #reduced sets.
C     -----------------------------------------------------------

      IF (MAXRED.LT.1 .OR. MAXVEC.LT.1) THEN
         XXBMX = -1.0D8
         XXBT  = 0.0D0
         DO ISYM = 1,NSYM
            XXB(ISYM) = 0.0D0
            DO ISYMB = 1,NSYM
               ISYMA = MULD2H(ISYMB,ISYM)
               IF (ISYMA .EQ. ISYMB) THEN
                  XA = DBLE(NBAS(ISYMA))
                  XXB(ISYM) = XXB(ISYM) + XA*(XA + 1.0D0)/2.0D0 
               ELSE IF (ISYMA .GT. ISYMB) THEN
                  XA = DBLE(NBAS(ISYMA))
                  XB = DBLE(NBAS(ISYMB))
                  XXB(ISYM) = XXB(ISYM) + XA*XB
               END IF
            END DO
            XXBT  = XXBT + XXB(ISYM)     ! total diag. dim.
            XXBMX = MAX(XXBMX,XXB(ISYM)) ! max. diag. block
         END DO
         IF (MAXVEC .LT. 1) THEN
            NBSMAX = NBAS(1)
            DO ISYM = 2,NSYM
               NBSMAX = MAX(NBSMAX,NBAS(ISYM))
            END DO
            MAXVEC = 20*NBSMAX ! default max. #vectors
            IF (XXBMX .LT. GBLIM) THEN
               NNBMX  = INT(XXBMX)
               MAXVEC = MIN(MAXVEC,NNBMX) ! reset if less than default
            END IF
         END IF
         IF (MAXRED .LT. 1) THEN
            MAXRED = NSYM*MAXVEC ! default max. #red. sets
            IF (XXBT .LT. GBLIM) THEN
               NNBT   = INT(XXBT)
               MAXRED = MIN(MAXRED,NNBT) ! reset if less than default
            END IF
         END IF
      END IF

      IF (MAXRED.LT.1 .OR. MAXVEC.LT.1) THEN
         WRITE(LUPRI,*) SECNAM,': MAXRED = ',MAXRED
         WRITE(LUPRI,*) SECNAM,': MAXVEC = ',MAXVEC
         CALL CHO_QUIT('MAXRED/MAXVEC error in '//SECNAM,103)
      ELSE
         LINFRD = MAXRED
         LINFVC = MAXVEC*3*NSYM
         CALL CHO_MEM('INFRED','ALLO','INTE',KINFRD,LINFRD,KFREE,LFREE)
         CALL CHO_MEM('INFVEC','ALLO','INTE',KINFVC,LINFVC,KFREE,LFREE)
      END IF

C     Initialize INFRED, INFVEC, vector counter, etc.
C     Special handling depending on Cholesky restart.
C     -----------------------------------------------

      NR  = MAXRED
      NV1 = MAXVEC
      NV2 = 3
      NV3 = NSYM
      CALL CHO_INIT1(WORK(KINFRD),NR,WORK(KINFVC),NV1,NV2,NV3)

C     Print header and configuration.
C     -------------------------------

      CALL CHO_PRTHEAD(.FALSE.)

C     Check configuration.
C     --------------------

      NCONFL = 0
      CALL CHO_CHKCONF(NCONFL,.TRUE.)
      IF (CHKONLY) THEN
         WRITE(LUPRI,*) SECNAM,': ',NCONFL,' conflicts detected in ',
     &                  'Cholesky decomposition configuration.'
         CALL CHO_QUIT('End of configuration check in '//SECNAM,100)
      ELSE IF (NCONFL .NE. 0) THEN
         WRITE(LUPRI,*) SECNAM,': ',NCONFL,' conflicts detected in ',
     &                  'Cholesky decomposition configuration.'
         CALL CHO_QUIT('Configuration conflicts in '//SECNAM,105)
      END IF

C     Set screening mode.
C     -------------------

      IF (CHO_DECALG .EQ. 1) THEN
         MODE_SCREEN = 1  ! damped screening
      ELSE IF (CHO_DECALG.EQ.2) THEN
         MODE_SCREEN = 2  ! remove diagonals < THRCOM
      ELSE IF (CHO_DECALG.EQ.3) THEN
         MODE_SCREEN = 2  ! remove diagonals < THRCOM
         IF (THRDIAG .EQ. 0.0D0) THEN
            THRDIAG = THRCOM ! discard initial diagonals > THRCOM
         END IF
      ELSE
         MODE_SCREEN = 1
      END IF

C     Open files for vector and reduced set storage.
C     Open restart files.
C     ----------------------------------------------

      CALL CHO_OPENVR(1)

C     Print section.
C     --------------

      IF (IPRINT .GE. 0) THEN

         CALL CHO_HEAD(STRING//SECNAM,LINE,80,LUPRI)

         WRITE(LUPRI,'(/,2X,A,I10)')
     &   'Number of irreps     : ',NSYM
         WRITE(LUPRI,'(2X,A,I10)')
     &   'Number of SOs        : ',NBAST
         WRITE(LUPRI,'(2X,A,I10)')
     &   'Number of shells     : ',NSHELL
         WRITE(LUPRI,'(2X,A,I10)')
     &   'Number of shell pairs: ',NSHELL*(NSHELL+1)/2
         WRITE(LUPRI,'(2X,A,I10)')
     &   'Max. shell dimension : ',MXORSH
         WRITE(LUPRI,'(2X,A,I10)')
     &   'Max. shell pair dim. : ',MX2SH

         IF (IPRINT .GE. 2) THEN

C           Basis size info.
C           ----------------

            WRITE(LUPRI,'(/,2X,A,/,2X,A)')
     &      '  Symmetry        NBAS        IBAS',
     &      '----------------------------------'
            DO ISYM = 1,NSYM
               WRITE(LUPRI,'(2X,I10,2X,I10,2X,I10)')
     &         ISYM,NBAS(ISYM),IBAS(ISYM)
            END DO
            WRITE(LUPRI,'(2X,A)')
     &      '----------------------------------'

         END IF

         IF (IPRINT .GE. 25) THEN

C           Shell info.
C           -----------

            WRITE(LUPRI,'(/,2X,A,/,2X,A,/,2X,A)')
     &     '     Shell   Dimension    Symmetry   Dimension      Offset',
     &     '             (NBSTSH)                (NBASSH)     (IBASSH)',
     &     '----------------------------------------------------------'
            DO ISHL = 1,NSHELL
               DO ISYM = 1,NSYM
                  IF (ISYM .EQ. 1) THEN
                     WRITE(LUPRI,'(2X,I10,2X,I10,2X,I10,2X,I10,2X,I10)')
     &               ISHL,NBSTSH(ISHL),
     &               ISYM,NBASSH(ISYM,ISHL),IBASSH(ISYM,ISHL)
                  ELSE
                     WRITE(LUPRI,'(26X,I10,2X,I10,2X,I10)')
     &               ISYM,NBASSH(ISYM,ISHL),IBASSH(ISYM,ISHL)
                  END IF
               END DO
            END DO
            WRITE(LUPRI,'(2X,A)')
     &     '----------------------------------------------------------'

            WRITE(LUPRI,'(/,2X,A,/,2X,A,/,2X,A)')
     &      '    SO        SO    sym    Shell     Index ',
     &      ' (global) (reduced)      (ISOSHL)  (ISHLSO)',
     &      '-------------------------------------------'
            DO ISYM = 1,NSYM
               DO I = 1,NBAS(ISYM)
                  IA = IBAS(ISYM) + I
                  WRITE(LUPRI,'(2X,I9,1X,I9,1X,I3,1X,I9,1X,I9)')
     &                 IA,I,ISYM,ISOSHL(IA),ISHLSO(IA)
               END DO
            END DO
            WRITE(LUPRI,'(2X,A)')
     &      '-------------------------------------------'

         END IF

      END IF

      END
C  /* Deck cho_init1 */
      SUBROUTINE CHO_INIT1(INFRED,NR,INFVEC,NV1,NV2,NV3)
C
C     Purpose: initialize counter arrays.
C
#include <implicit.h>
      INTEGER INFRED(NR), INFVEC(NV1,NV2,NV3)
#include <cholesky.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CHO_INIT1')

      INTEGER CHO_ISUMELM

      IF (NR .NE. MAXRED) THEN
         CALL CHO_QUIT('Dimension problem 0 in '//SECNAM,101)
      END IF
      IF (NV1 .NE. MAXVEC) THEN
         CALL CHO_QUIT('Dimension problem 1 in '//SECNAM,101)
      END IF
      IF (NV2 .NE. 3) THEN
         CALL CHO_QUIT('Dimension problem 2 in '//SECNAM,101)
      END IF
      IF (NV3 .NE. NSYM) THEN
         CALL CHO_QUIT('Dimension problem 3 in '//SECNAM,101)
      END IF

      IF (RSTCHO) THEN

C        Read restart information.
C        -------------------------

         CALL CHO_GETRSTC(INFRED,NR,INFVEC,NV1,NV2,NV3)
         NUMCHT = CHO_ISUMELM(NUMCHO,NSYM)

      ELSE

C        Initialize vector info and counters.
C        ------------------------------------

         CALL CHO_IZERO(INFVEC,NV1*NV2*NV3)
         CALL CHO_IZERO(NUMCHO,NSYM)
         NUMCHT = 0

C        Initialize reduced set info.
C        ----------------------------

         CALL CHO_IZERO(INFRED,NR)

C        Initialize global integral pass counter.
C        ----------------------------------------

         XNPASS = 0

      END IF

      END
C  /* Deck cho_openvr */
      SUBROUTINE CHO_OPENVR(IOPT)
C
C     Purpose: open (IOPT=1) or close (IOPT=2) files for vector
C              and reduced set storage as well as restart files.
C
#include <implicit.h>
#include <cholesky.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_OPENVR')

      CALL CHO_DAL_OPENVR(IOPT)

      END
C  /* Deck cho_calcdiag */
      SUBROUTINE CHO_CALCDIAG(BUF,IBUF,LENBUF,SCR,LENSCR,WORK,LWORK,
     &                        IIBSTRSH,NNBSTRSH,MSYM,MMSHL,
     &                        NDUMP)
C
C     Purpose: shell-driven calculation of the integral diagonal and
C              setup of the first reduced set.
C
#include <implicit.h>
      DIMENSION BUF(LENBUF), SCR(LENSCR), WORK(LWORK)
      INTEGER   IBUF(4,LENBUF), IIBSTRSH(MSYM,MMSHL,3),
     &          NNBSTRSH(MSYM,MMSHL,3)
#include <cholesky.h>
#include <choorb.h>

      CHARACTER*12 SECNAM
      PARAMETER (SECNAM = 'CHO_CALCDIAG')

      PARAMETER (INFOD = 50)

      DIMENSION SCRMAX(8)

      INTEGER CHO_ISAOSH

      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J

C     Check dimensions.
C     -----------------

      IF (LENBUF .LT. LBUF) THEN
         WRITE(LUPRI,'(//,1X,A,A)') SECNAM,': LENBUF >= LBUF required!'
         WRITE(LUPRI,'(1X,A,I10)')    'LENBUF = ',LENBUF
         WRITE(LUPRI,'(1X,A,I10,/)')  'LBUF   = ',LBUF
         CALL CHO_QUIT('Buffer error in '//SECNAM,102)
      END IF
      LSCR = MX2SH
      IF (LENSCR .LT. LSCR) THEN
         WRITE(LUPRI,'(//,1X,A,A)') SECNAM,
     &   ': LENSCR >= MX2SH required!'
         WRITE(LUPRI,'(1X,A,I10)')    'LENSCR = ',LENSCR
         WRITE(LUPRI,'(1X,A,I10,/)')  'MX2SH  = ',LSCR
         CALL CHO_QUIT('Scratch space error in '//SECNAM,102)
      END IF
      IF (MSYM .NE. NSYM) THEN
         WRITE(LUPRI,'(//,1X,A,A)') SECNAM,': MSYM = NSYM required!'
         WRITE(LUPRI,'(1X,A,I10)')    'MSYM = ',MSYM
         WRITE(LUPRI,'(1X,A,I10,/)')  'NSYM = ',NSYM
         CALL CHO_QUIT('Symmetry error in '//SECNAM,102)
      END IF
      NNSHL = NSHELL*(NSHELL + 1)/2
      IF (MMSHL .NE. NNSHL) THEN
         WRITE(LUPRI,'(//,1X,A,A)') SECNAM,
     &   ': MMSHL = NSHELL*(NSHELL+1)/2 required!'
         WRITE(LUPRI,'(1X,A,I10)')    'MMSHL               = ',MMSHL
         WRITE(LUPRI,'(1X,A,I10,/)')  'NSHELL*(NSHELL+1)/2 = ',NNSHL
         CALL CHO_QUIT('Shell error in '//SECNAM,102)
      END IF

C     Open scratch files.
C     -------------------

      IUNIT = -1
      CALL CHO_OPEN(IUNIT,'CHO_DIASCR2')
      JUNIT = -1
      CALL CHO_OPEN(JUNIT,'CHO_DIASCR1')
      REWIND(IUNIT)
      REWIND(JUNIT)

C     Make JUNIT available outside this routine.
C     ------------------------------------------

      LUSCR = JUNIT

C     Initialize abs. max. diag. array.
C     ---------------------------------

      CALL CHO_DZERO(DIAMAX,NSYM)

C     Calculate diagonal in loop over shell-pairs.
C     --------------------------------------------

      XLDIAG = 0.0D0
      ICOUNT = 0
      NDUMP  = 0
      DO ISHLA = 1,NSHELL
         DO ISHLB = 1,ISHLA
 
            NUMA = NBSTSH(ISHLA)
            NUMB = NBSTSH(ISHLB)
            IF (ISHLA .EQ. ISHLB) THEN
               NUMAB = NUMA*(NUMA + 1)/2
            ELSE
               NUMAB = NUMA*NUMB
            END IF
 
            CALL CHO_DAL_GABN(ISHLA,ISHLB,SCR,WORK,LWORK)
 
            IF (IPRINT .GE. INFOD) THEN
               IF ((ISHLA.EQ.1) .AND. (ISHLB.EQ.1)) THEN
                  CALL CHO_HEAD(SECNAM//': Unscreened Diagonal',
     &                          '=',80,LUPRI)
               END IF
               WRITE(LUPRI,'(/,2X,A,I10,1X,I10,1X,I10)')
     &         'Diagonal shell block A,B,AB = ',ISHLA,ISHLB,
     &                                          INDEX(ISHLA,ISHLB)
               IF (ISHLA .EQ. ISHLB) THEN
                  CALL CHO_OUTPAK(SCR,NUMA,1,LUPRI)
               ELSE 
                  CALL CHO_OUTPUT(SCR,1,NUMA,1,NUMB,NUMA,NUMB,1,LUPRI)
               END IF                   
            END IF

            IF (ISHLA .EQ. ISHLB) THEN
               DO IA = 1,NUMA
                  ISYMA = CHO_ISAOSH(IA,ISHLA)
                  DO IB = 1,IA
                     ISYMB  = CHO_ISAOSH(IB,ISHLB)
                     ISYMAB = MULD2H(ISYMB,ISYMA) 
                     IAB    = INDEX(IA,IB)
                     DIAAB  = ABS(SCR(IAB))
                     IF (DIAAB .GT. THRDIAG) THEN
                        DIAMAX(ISYMAB) = MAX(DIAMAX(ISYMAB),DIAAB)
                        ICOUNT = ICOUNT + 1
                        BUF(ICOUNT)    = SCR(IAB)
                        IBUF(1,ICOUNT) = INDEX(ISHLA,ISHLB)
                        IBUF(2,ICOUNT) = IAB
                        IBUF(3,ICOUNT) = ISYMAB
                        IBUF(4,ICOUNT) = IAB
                        IF (ICOUNT .EQ. LBUF) THEN
                           CALL CHO_WRBUF(LBUF,BUF,IBUF,LBUF,IUNIT) 
                           XLDIAG = XLDIAG + DBLE(LBUF)
                           ICOUNT = 0
                           NDUMP  = NDUMP + 1
                        END IF
                     END IF
                  END DO
               END DO
            ELSE
               DO ISYMB = 1,NSYM
                  DO IBB = 1,NBASSH(ISYMB,ISHLB)
                     IB = IBASSH(ISYMB,ISHLB) + IBB
                     DO ISYMA = 1,NSYM
                        DO IAA = 1,NBASSH(ISYMA,ISHLA)
                           IA = IBASSH(ISYMA,ISHLA) + IAA
                           ISYMAB = MULD2H(ISYMA,ISYMB) 
                           IAB    = NUMA*(IB - 1) + IA
                           DIAAB  = ABS(SCR(IAB))
                           IF (DIAAB .GT. THRDIAG) THEN
                              DIAMAX(ISYMAB) = MAX(DIAMAX(ISYMAB),DIAAB)
                              ICOUNT = ICOUNT + 1
                              BUF(ICOUNT)    = SCR(IAB)
                              IBUF(1,ICOUNT) = INDEX(ISHLA,ISHLB)
                              IBUF(2,ICOUNT) = IAB
                              IBUF(3,ICOUNT) = ISYMAB
                              IBUF(4,ICOUNT) = IAB
                              IF (ICOUNT .EQ. LBUF) THEN
                                 CALL CHO_WRBUF(LBUF,BUF,IBUF,LBUF,
     &                                          IUNIT)
                                 XLDIAG = XLDIAG + DBLE(LBUF)
                                 ICOUNT = 0
                                 NDUMP  = NDUMP + 1
                              END IF
                           END IF
                        END DO
                     END DO
                  END DO
               END DO
            END IF
 
         END DO
      END DO
 
      IF (ICOUNT .GT. 0) THEN
         IF (ICOUNT .GT. LBUF) THEN
            CALL CHO_QUIT('Logical error in '//SECNAM,103)
         END IF
         CALL CHO_WRBUF(ICOUNT,BUF,IBUF,LBUF,IUNIT)
         XLDIAG = XLDIAG + DBLE(ICOUNT)
         ICOUNT = 0
         NDUMP  = NDUMP + 1
      END IF

      IF ((XLDIAG.LT.1.0D0) .OR. (NDUMP.LT.1)) THEN
         WRITE(LUPRI,*) SECNAM,' Error:'
         IF (XLDIAG .LT. 1.0D0) THEN
            WRITE(LUPRI,*) '   Diagonal dimension is ',XLDIAG
         END IF
         IF (NDUMP .LT. 1) THEN
            WRITE(LUPRI,*) '   Nothing written to disk; NDUMP = ',NDUMP
         END IF
         CALL CHO_QUIT('Error in '//SECNAM,104)
      END IF
      IF (IPRINT.GE.INFOD .OR. THRDIAG.GT.0.0D0) THEN
         CALL CHO_HEAD(SECNAM//': Diagonal Info','=',80,LUPRI)
         XXX   = DBLE(NBAST)
         XMDIA = XXX*(XXX + 1.0D0)/2.0D0
         XLDIA = XLDIAG
         SAVD  = 1.0D2*(XMDIA - XLDIA)/XMDIA
         WRITE(LUPRI,'(/,2X,A,1P,D15.6)')
     &   'Screening threshold for initial diagonal: ',THRDIAG
         WRITE(LUPRI,'(2X,A,F15.1,/,2X,A,F15.1)')
     &   'Dimension of unscreened initial diagonal: ',XMDIA,
     &   'Dimension of   screened initial diagonal: ',XLDIA
         WRITE(LUPRI,'(2X,A,7X,F8.3,A)')
     &   'Saving from screening                   : ',SAVD,'%'
         DO ISYM = 1,NSYM
            WRITE(LUPRI,'(2X,A,I2,12X,A,1P,D15.6)')
     &      'Maximum diagonal, symmetry',ISYM,': ',DIAMAX(ISYM)
         END DO
      END IF

C     Read through the file to get first reduced set.
C     -----------------------------------------------

      DEL1 = THRCOM*THRCOM/DAMP(1)
      DO ISYM = 1,NSYM
         IF (ABS(DIAMAX(ISYM)) .GT. 0.0D0) THEN
            SCRMAX(ISYM) = DEL1/DIAMAX(ISYM)
         ELSE
            SCRMAX(ISYM) = 1.0D40
         END IF
      END DO
      CALL CHO_IZERO(NNBSTRSH(1,1,1),MSYM*MMSHL)

      REWIND(IUNIT)
      REWIND(JUNIT)
      IF (CHO_USEABS) THEN
         DO IDUMP = 1,NDUMP
            CALL CHO_RDBUF(LENGTH,BUF,IBUF,LBUF,IUNIT)
            IF (IDUMP .EQ. NDUMP) THEN
               CALL CHO_CLOSE(IUNIT,'DELETE')
            END IF
            DO L = 1,LENGTH
               DIAGAB = BUF(L)
               ISYMAB = IBUF(3,L)
               IF (ABS(DIAGAB) .LT. SCRMAX(ISYMAB)) THEN
                  BUF(L)    = 0.0D0
                  IBUF(2,L) = -1
               ELSE
                  ISHLAB = IBUF(1,L)
                  NNBSTRSH(ISYMAB,ISHLAB,1) = NNBSTRSH(ISYMAB,ISHLAB,1)
     &                                      + 1
                  IBUF(2,L) = NNBSTRSH(ISYMAB,ISHLAB,1)
               END IF
            END DO
            CALL CHO_WRBUF(LENGTH,BUF,IBUF,LBUF,JUNIT)
         END DO
      ELSE
         DO IDUMP = 1,NDUMP
            CALL CHO_RDBUF(LENGTH,BUF,IBUF,LBUF,IUNIT)
            IF (IDUMP .EQ. NDUMP) THEN
               CALL CHO_CLOSE(IUNIT,'DELETE')
            END IF
            DO L = 1,LENGTH
               DIAGAB = BUF(L)
               ISYMAB = IBUF(3,L)
               IF (DIAGAB .LT. SCRMAX(ISYMAB)) THEN
                  BUF(L)    = 0.0D0
                  IBUF(2,L) = -1
               ELSE
                  ISHLAB = IBUF(1,L)
                  NNBSTRSH(ISYMAB,ISHLAB,1) = NNBSTRSH(ISYMAB,ISHLAB,1)
     &                                      + 1
                  IBUF(2,L) = NNBSTRSH(ISYMAB,ISHLAB,1)
               END IF
            END DO
            CALL CHO_WRBUF(LENGTH,BUF,IBUF,LBUF,JUNIT)
         END DO
      END IF

      CALL CHO_SETREDIND(IIBSTRSH,NNBSTRSH,MSYM,MMSHL,1)

      END
C  /* Deck cho_getdiag1 */
      SUBROUTINE CHO_GETDIAG1(IIBSTRSH,NNBSTRSH,INDRED,INDRSH,
     &                        IBUF,BUF,DIAG,
     &                        MSYM,MMSHL,MMBSTRT,LENBUF,NDUMP)
C
C     Purpose: read diagonal in first reduced set.
C
#include <implicit.h>
      INTEGER   IIBSTRSH(MSYM,MMSHL,3), NNBSTRSH(MSYM,MMSHL,3)
      INTEGER   INDRED(MMBSTRT,3), INDRSH(MMBSTRT),IBUF(4,LENBUF)
      DIMENSION DIAG(MMBSTRT), BUF(LENBUF)
#include <cholesky.h>

      CHARACTER*12 SECNAM
      PARAMETER (SECNAM = 'CHO_GETDIAG1')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.) 

      INTEGER ISYLST(8)

      PARAMETER (INFOD = 50)
      PARAMETER (TINY  = 1.0D-14)

C     Check dimensions.
C     -----------------

      IF (MSYM .NE. NSYM) THEN
         WRITE(LUPRI,'(//,1X,A,A)') SECNAM,': MSYM = NSYM required!'
         WRITE(LUPRI,'(1X,A,I10)')    'MSYM = ',MSYM
         WRITE(LUPRI,'(1X,A,I10,/)')  'NSYM = ',NSYM
         CALL CHO_QUIT('Symmetry error in '//SECNAM,102)
      END IF
      NNSHL = NSHELL*(NSHELL + 1)/2
      IF (MMSHL .NE. NNSHL) THEN
         WRITE(LUPRI,'(//,1X,A,A)') SECNAM,
     &   ': MMSHL = NSHELL*(NSHELL+1)/2 required!'
         WRITE(LUPRI,'(1X,A,I10)')    'MMSHL               = ',MMSHL
         WRITE(LUPRI,'(1X,A,I10,/)')  'NSHELL*(NSHELL+1)/2 = ',NNSHL
         CALL CHO_QUIT('Shell error in '//SECNAM,102)
      END IF
      IF (MMBSTRT .NE. NNBSTRT(1)) THEN
         WRITE(LUPRI,'(//,1X,A,A)') SECNAM,
     &   ': MMBSTRT = NNBSTRT(1) required!'
         WRITE(LUPRI,'(1X,A,I10)')    'MMBSTRT    = ',MMBSTRT
         WRITE(LUPRI,'(1X,A,I10,/)')  'NNBSTRT(1) = ',NNBSTRT(1)
         CALL CHO_QUIT('Dimension error in '//SECNAM,102)
      END IF
      IF (.NOT. RSTDIA) THEN
         IF (LENBUF .LT. LBUF) THEN
            WRITE(LUPRI,'(//,1X,A,A)') SECNAM,
     &                                 ': LENBUF >= LBUF required!'
            WRITE(LUPRI,'(1X,A,I10)')    'LENBUF = ',LENBUF
            WRITE(LUPRI,'(1X,A,I10,/)')  'LBUF   = ',LBUF
            CALL CHO_QUIT('Buffer error in '//SECNAM,102)
         END IF
      END IF

C     Read diagonal from file.
C     ------------------------

      IF (RSTDIA) THEN

         IOPT = 2
         CALL CHO_IODIAG(DIAG,IOPT)

      ELSE

         IUNIT = LUSCR
         LUSCR = -1
         REWIND(IUNIT)

         DO IDUMP = 1,NDUMP
            CALL CHO_RDBUF(LENGTH,BUF,IBUF,LBUF,IUNIT)
            IF (IDUMP .EQ. NDUMP) THEN
                CALL CHO_CLOSE(IUNIT,'DELETE')
            END IF
            DO L = 1,LENGTH
               IF (IBUF(2,L) .GT. 0) THEN
                  ISHLAB = IBUF(1,L)
                  ISYMAB = IBUF(3,L)
                  IAB    = IIBSTR(ISYMAB,1) + IIBSTRSH(ISYMAB,ISHLAB,1)
     &                   + IBUF(2,L)
                  DIAG(IAB) = BUF(L)
                  INDRSH(IAB) = ISHLAB
                  INDRED(IAB,1) = IBUF(4,L)
               END IF
            END DO
         END DO

      END IF

C     Copy info to current reduced set (IRED=2).
C     Also set up IRED=3 (although it should be redundant).
C     -----------------------------------------------------

      LEN2 = NSYM*NNSHL
      LEN3 = NSYM
      DO IRED = 2,3
         CALL ICOPY(LEN2,IIBSTRSH(1,1,1),1,IIBSTRSH(1,1,IRED),1)
         CALL ICOPY(LEN2,NNBSTRSH(1,1,1),1,NNBSTRSH(1,1,IRED),1)
         CALL ICOPY(LEN3,IIBSTR(1,1),1,IIBSTR(1,IRED),1)
         CALL ICOPY(LEN3,NNBSTR(1,1),1,NNBSTR(1,IRED),1)
         NNBSTRT(IRED) = NNBSTRT(1)
         DO IAB = 1,NNBSTRT(1)
            INDRED(IAB,IRED) = IAB
         END DO
      END DO

#if defined (CHODAL_DEBUG)
C     Set up array mapping from 1st reduced set to full storage.
C     ----------------------------------------------------------

      CALL CHO_DAL_SETRS1TOF(IIBSTRSH,NNBSTRSH,INDRED,
     &                       MSYM,MMSHL,NNBSTRT(1))
      WRITE(LUPRI,*) 'Initial diagonal: alpha,beta,diag'
      XSUM = 0.0D0
      XWSM = 0.0D0
      XXSM = 0.0D0
      DO I = 1,NNBSTRT(1)
         WRITE(LUPRI,*) IRS1TOF(1,I),IRS1TOF(2,I),DIAG(I)
         XSUM = XSUM + DIAG(I)
         XWSM = XWSM + DIAG(I)*DBLE(IRS1TOF(1,I))/DBLE(IRS1TOF(2,I))
         XXSM = XXSM + DIAG(I)*DBLE(IRS1TOF(2,I))/DBLE(IRS1TOF(1,I))
      END DO
      XNRM = SQRT(DDOT(NNBSTRT(1),DIAG,1,DIAG,1))
      WRITE(LUPRI,*) 'Sum  of diagonal elements: ',XSUM
      WRITE(LUPRI,*) 'Wsum of diagonal elements: ',XWSM
      WRITE(LUPRI,*) 'Xsum of diagonal elements: ',XXSM
      WRITE(LUPRI,*) 'Norm of diagonal elements: ',XNRM
#endif

C     Print.
C     ------

      IF (LOCDBG .OR. (IPRINT.GE.INFOD)) THEN
         DO ISYM = 1,NSYM
            ISYLST(ISYM) = ISYM
         END DO
         NSYLST = NSYM
         IRED = 1
         CALL CHO_PRTDIA(DIAG,IIBSTRSH,NNBSTRSH,INDRED,INDRSH,ISYLST,
     &                   MSYM,MMSHL,MMBSTRT,NSYLST,IRED)
         IF (LOCDBG) THEN
            IRED = 2
            CALL CHO_PRTDIA(DIAG,IIBSTRSH,NNBSTRSH,INDRED,INDRSH,ISYLST,
     &                      MSYM,MMSHL,MMBSTRT,NSYLST,IRED)
         END IF
      END IF

      END
C  /* Deck cho_anadia */
      SUBROUTINE CHO_ANADIA(DIAG,BIN1,STEP,NUMBIN,FULL)
C
C     Purpose: analyze diagonal (histogram).
C
#include <implicit.h>
      DIMENSION DIAG(*)
      LOGICAL   FULL
#include <cholesky.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_ANADIA')

      PARAMETER (MAXBIN = 50)
      DIMENSION BIN(MAXBIN)

      LOGICAL FOUND

C     Print header.
C     -------------

      CALL CHO_HEAD('Histogram of Diagonal Elements','=',80,LUPRI)

C     Set up size bins for analysis of diagonal.
C     ------------------------------------------

      IF (NUMBIN .LT. 1) THEN
         MBIN = MIN(10,MAXBIN)
         BINLOC = 1.0D2
         STPLOC = 1.0D-2
      ELSE
         MBIN = MIN(NUMBIN,MAXBIN)
         BINLOC = BIN1
         STPLOC = STEP
      END IF

      BIN(1) = BINLOC
      DO IBIN = 2,MBIN
         BIN(IBIN) = BIN(IBIN-1)*STPLOC
      END DO

C     Set smallest BIN according to full.
C     -----------------------------------

      IF (FULL) THEN
         NBIN = MBIN
      ELSE
         NBIN  = MBIN
         IBIN  = MBIN
         FOUND = .FALSE.
         DO WHILE ((IBIN.GT.1) .AND. (.NOT.FOUND))
            IBIN = IBIN - 1
            IF (THRCOM .GE. BIN(IBIN)) THEN
               NBIN = IBIN + 1
            ELSE
               FOUND = .TRUE.
            END IF
         END DO
      END IF

C     Histogram.
C     ----------

      CALL CHO_ANASIZE(DIAG,NNBSTRT(1),BIN,NBIN,LUPRI)

C     Count converged.
C     ----------------

      NCONV = 0
      DO IAB = 1,NNBSTRT(1)
         IF (DIAG(IAB) .LE. THRCOM) NCONV = NCONV + 1
      END DO
      WRITE(LUPRI,'(/,1X,A,I10,/,1X,A,I10)')
     & 'Converged  : ',NCONV,'Unconverged: ',NNBSTRT(1)-NCONV

C     Print total number of negative zeroed diagonal as well as the most
C     negative one.
C     ------------------------------------------------------------------

      WRITE(LUPRI,'(/,1X,A,5X,I10)')
     & 'Total number of zeroed negative diagonals: ',NNZTOT
      IF (NNZTOT .GT. 0) THEN
         IF (IABMNZ .LT. 1) THEN
            WRITE(LUPRI,'(1X,A)')
     &     'WARNING: most negative zeroed diagonal has not been stored!'
         ELSE
            WRITE(LUPRI,'(1X,A,1P,D15.6)')
     &      '- most negative zeroed diagonal          : ',DIAMNZ
         END IF
      END IF

      END
C  /* Deck cho_decdrv */
      SUBROUTINE CHO_DECDRV(INFRED,INFVEC,
     &                      IIBSTRSH,NNBSTRSH,ISCR,INDRED,INDRSH,INTMAP,
     &                      DIAG,WORK,LWORK,
     &                      MRED,MVEC,MSYM,MMSHL,MMBSTRT,LISCR,
     &                      KFRIN,LFRIN)
#include <implicit.h>
      INTEGER   IIBSTRSH(MSYM,MMSHL,3), NNBSTRSH(MSYM,MMSHL,3)
      INTEGER   ISCR(LISCR)
      INTEGER   INDRED(MMBSTRT,3), INDRSH(MMBSTRT), INTMAP(MMSHL)
      INTEGER   INFRED(MRED), INFVEC(MVEC,3,MSYM)
      DIMENSION DIAG(*), WORK(LWORK)
C
C     Purpose: driver for the decomposition of the two-electron integral
C              matrix based on the reduced diagonal.
C
#include <cholesky.h>

      INTEGER ISYLST(8)

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_DECDRV')

      CHARACTER*(8) FILSEL(8)
      DATA FILSEL /'CHOSEL_1','CHOSEL_2','CHOSEL_3','CHOSEL_4',
     &             'CHOSEL_5','CHOSEL_6','CHOSEL_7','CHOSEL_8'/

      LOGICAL CONV, LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

C     Start timing.
C     -------------

      CALL CHO_TIMER(TCPU1,TWALL1)

C     Set flag.
C     ---------

      DID_DECDRV = .TRUE.

C     Initializations and static settings.
C     IRED=2: points to current reduced set in index arrays.
C     ------------------------------------------------------

      IRED  = 2
      NNSHL = NSHELL*(NSHELL + 1)/2
      N2SHL = NSHELL*NSHELL
      CONV  = .FALSE.
      KFREE = KFRIN
      LFREE = LFRIN

C     Check dimensions.
C     -----------------

      IF (MSYM .NE. NSYM) THEN
         WRITE(LUPRI,'(//,1X,A,A)') SECNAM,': MSYM = NSYM required!'
         WRITE(LUPRI,'(1X,A,I10)')    'MSYM = ',MSYM
         WRITE(LUPRI,'(1X,A,I10,/)')  'NSYM = ',NSYM
         CALL CHO_QUIT('Symmetry error in '//SECNAM,102)
      END IF
      IF (MMSHL .NE. NNSHL) THEN
         WRITE(LUPRI,'(//,1X,A,A)') SECNAM,
     &   ': MMSHL = NSHELL*(NSHELL+1)/2 required!'
         WRITE(LUPRI,'(1X,A,I10)')    'MMSHL               = ',MMSHL
         WRITE(LUPRI,'(1X,A,I10,/)')  'NSHELL*(NSHELL+1)/2 = ',NNSHL
         CALL CHO_QUIT('Shell error in '//SECNAM,102)
      END IF
      IF (MMBSTRT .NE. NNBSTRT(1)) THEN
         WRITE(LUPRI,'(//,1X,A,A)') SECNAM,
     &   ': MMBSTRT = NNBSTRT(1) required!'
         WRITE(LUPRI,'(1X,A,I10)')    'MMBSTRT    = ',MMBSTRT
         WRITE(LUPRI,'(1X,A,I10,/)')  'NNBSTRT(1) = ',NNBSTRT(1)
         CALL CHO_QUIT('Dimension error in '//SECNAM,102)
      END IF

C     Allocate memory for shell pair based diagonal.
C     ----------------------------------------------

      CALL CHO_MEM('DIASH','ALLO','REAL',KDIASH,N2SHL,KFREE,LFREE)
      KFLUSH = KDIASH
      CALL CHO_MEM('ISYSH','ALLO','INTE',KISYSH,N2SHL,KFREE,LFREE)

C     Set first integral pass.
C     ------------------------

      NPOTSH = 0
      CALL CHO_SETPASS(DIAG,WORK(KDIASH),WORK(KISYSH),
     &                 INDRED,INDRSH,IRED,NSHELL,NNBSTRT(1),
     &                 CONV,NPOTSH)

C     Loop over integral passes. Continue until convergence or
C     until the max. number of integral passes has been reached.
C     To each integral pass there is associated a reduced set,
C     so the IPASS counter is also used as identifier of reduced
C     set during I/O.
C     ----------------------------------------------------------

      IPASS = XNPASS
      JPASS = 0
      MPASS = NNBSTRT(IRED)
      DO WHILE ((.NOT.CONV) .AND. (JPASS.LT.MPASS))

C        Update integral pass counter.
C        -----------------------------

         JPASS = JPASS + 1
         IPASS = IPASS + 1

C        Debug: print diagonal.
C        ----------------------

         IF (LOCDBG) THEN
            WRITE(LUPRI,*) SECNAM,': debug: diagonal before pass ',
     &                     IPASS
            DO ISYM = 1,NSYM
               ISYLST(ISYM) = ISYM
            END DO
            CALL CHO_PRTDIA(DIAG,IIBSTRSH,NNBSTRSH,INDRED,INDRSH,ISYLST,
     &                      NSYM,NNSHL,NNBSTRT(1),NSYM,IRED)
            WRITE(LUPRI,*)
            WRITE(LUPRI,*) SECNAM,': INFRED before pass ',IPASS
            WRITE(LUPRI,'(10I8)') (INFRED(I),I=1,MIN(IPASS,MAXRED))
            CALL CHO_FLUSH(LUPRI)
         END IF

C        Write index arrays for reduced set to disk
C        and update disk address.
C        ------------------------------------------

         IF (IPASS .GT. MAXRED) THEN
            IF (RSTCHO) THEN
               WRITE(LUPRI,*) SECNAM,': total integral pass ',IPASS
               WRITE(LUPRI,*) SECNAM,': local integral pass ',JPASS
            ELSE
               WRITE(LUPRI,*) SECNAM,': integral pass ',IPASS
            END IF
            WRITE(LUPRI,*) SECNAM,': max. allowed is ',MAXRED
            WRITE(LUPRI,*) SECNAM,': please increase max. allowed!'
            CALL CHO_QUIT('Too many integral passes in '//SECNAM,104)
         ELSE IF (IPASS .EQ. 1) THEN
            CALL CHO_PUTRED(INFRED,
     &                      NNBSTRSH(1,1,1),INDRED(1,1),INDRSH,
     &                      MRED,NSYM,NNSHL,NNBSTRT(1),IPASS)
            INFRED(IPASS+1) = INFRED(IPASS)
     &                      + NSYM*NNSHL + 2*NNBSTRT(1)
         ELSE IF (IPASS .EQ. MAXRED) THEN
            CALL CHO_PUTRED(INFRED,
     &                      NNBSTRSH(1,1,IRED),INDRED(1,IRED),INDRSH,
     &                      MRED,NSYM,NNSHL,NNBSTRT(1),IPASS)
         ELSE
            CALL CHO_PUTRED(INFRED,
     &                      NNBSTRSH(1,1,IRED),INDRED(1,IRED),INDRSH,
     &                      MRED,NSYM,NNSHL,NNBSTRT(1),IPASS)
            INFRED(IPASS+1) = INFRED(IPASS)
     &                      + NSYM*NNSHL + NNBSTRT(IRED)
         END IF

C        Open scratch files for qualified integral columns.
C        --------------------------------------------------

         DO ISYM = 1,NSYM
            LUSEL(ISYM) = -1
            CALL CHO_OPEN(LUSEL(ISYM),FILSEL(ISYM))
         END DO

C        Get integral columns on disk stored in current reduced set.
C        -----------------------------------------------------------

         NUM = 0

         KFR = 1
         LFR = LFREE
         CALL CHO_GETINT(DIAG,WORK(KDIASH),WORK(KISYSH),
     &                   INDRED,INDRSH,IIBSTRSH,NNBSTRSH,INTMAP,
     &                   NSYM,NSHELL,NNSHL,NNBSTRT(1),
     &                   WORK(KFREE),LFREE,NPOTSH,NUM,
     &                   KFR,LFR)
         CALL CHO_FLUSH(LUPRI)

C        Allocate max. work space for decomposition.
C        -------------------------------------------

         CALL CHO_MEM('MaxMem','MAX ','REAL',KWRK,LWRK,KFREE,LFREE)

C        Decompose the qualified integral columns.
C        -----------------------------------------

         CALL CHO_DECOM(DIAG,INFRED,INFVEC,INDRED,INDRSH,
     &                  IIBSTRSH,NNBSTRSH,ISCR,
     &                  MRED,MVEC,NSYM,NNSHL,NNBSTRT(1),LISCR,
     &                  WORK(KWRK),LWRK,IPASS,NUM)
         CALL CHO_FLUSH(LUPRI)

C        Write restart info to disk.
C        ---------------------------

         CALL CHO_WRRSTC(INFRED,MRED,INFVEC,MVEC,3,NSYM,
     &                   IPASS,INTMAP,NNSHL)

C        Free memory.
C        ------------

         CALL CHO_MEM('MaxMem','FREE','REAL',KWRK,LWRK,KFREE,LFREE)

C        Delete scratch files for qualified integral columns.
C        ----------------------------------------------------

         DO ISYM = 1,NSYM
            IF (LUSEL(ISYM) .GT. 0) THEN
               CALL CHO_CLOSE(LUSEL(ISYM),'DELETE')
            END IF
         END DO

C        Analyze diagonal.
C        -----------------

         BIN1 = 1.0D2
         STEP = 1.0D-1
         NBIN = 18
         CALL CHO_ANADIA(DIAG,BIN1,STEP,NBIN,.FALSE.)

C        Get next reduced set.
C        ---------------------

         CALL CHO_SETRED(DIAG,IIBSTRSH,NNBSTRSH,INDRED,INDRSH,
     &                   NSYM,NNBSTRT(1),NNSHL)
         CALL CHO_PRTRED(NNBSTRSH,NSYM,NNSHL,2)
         CALL CHO_FLUSH(LUPRI)

C        Check convergence and, if not converged, set next integral pass.
C        ----------------------------------------------------------------

         NPOTSH = 0
         CALL CHO_SETPASS(DIAG,WORK(KDIASH),WORK(KISYSH),
     &                    INDRED,INDRSH,IRED,NSHELL,
     &                    NNBSTRT(1),CONV,NPOTSH)

      END DO

C     Free memory for shell pair based diagonal.
C     ------------------------------------------

      CALL CHO_MEM('DIASH','FLUSH','REAL',KDIASH,N2SHL,KFREE,LFREE)

C     Set flag for statistics.
C     ------------------------

      DID_DECDRV = .TRUE.

C     Timing.
C     -------

      CALL CHO_TIMER(TCPU2,TWALL2)
      TDECDRV(1) = TCPU2  - TCPU1
      TDECDRV(2) = TWALL2 - TWALL1

      END
C  /* Deck cho_setpass */
      SUBROUTINE CHO_SETPASS(DIAG,DIASH,ISYSH,
     &                       INDRED,INDRSH,IRED,MSH,MMBSTRT,CONV,
     &                       NPOTSH)
C
C     Purpose: Check convergence and, if not converged, set up
C              integral pass.
C
#include <implicit.h>
      DIMENSION DIAG(*), DIASH(MSH,MSH)
      INTEGER   ISYSH(MSH,MSH), INDRED(MMBSTRT,3), INDRSH(MMBSTRT)
      LOGICAL   CONV
#include <cholesky.h>

C     Initialize the potential number of shell pairs that can contribute.
C     -------------------------------------------------------------------

      NPOTSH = 0

C     Find max. abs. diagonal in each symmetry and the global max.
C     ------------------------------------------------------------

      DGMAX = -1.0D15
      CALL CHO_MAXABSDIAG(DIAG,INDRED,INDRSH,IRED,MMBSTRT,DGMAX)

C     If not converged, set next integral pass.
C     -----------------------------------------

      CONV = DGMAX .LT. THRCOM
      IF (.NOT. CONV) THEN
         CALL CHO_SETMAXSHL(DIAG,DIASH,ISYSH,INDRED,INDRSH,
     &                      IRED,MSH,MMBSTRT)
         DO ISYM = 1,NSYM
            DIAMIN(ISYM) = MAX(DIAMAX(ISYM)*SPAN,THRCOM)
         END DO
         DO ISHLB = 1,NSHELL
            DO ISHLA = 1,ISHLB
               IF (DIASH(ISHLA,ISHLB) .GT. 0.0D0) NPOTSH = NPOTSH + 1
            END DO
         END DO
      END IF

      END
C  /* Deck cho_maxabsdiag */
      SUBROUTINE CHO_MAXABSDIAG(DIAG,INDRED,INDRSH,IRED,MMBSTRT,DGMAX)
C
C     Purpose: set max. abs. DIAG (reduced set IRED) in each symmetry, and
C              return the global max. abs. in DGMAX.
C
#include <implicit.h>
      DIMENSION DIAG(*)
      INTEGER   INDRED(MMBSTRT,3), INDRSH(MMBSTRT)
#include <cholesky.h>

      CHARACTER*14 SECNAM
      PARAMETER (SECNAM = 'CHO_MAXABSDIAG')

      INTEGER AB, AB1, AB2

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      IF (IRED .EQ. 1) THEN
         DO ISYM = 1,NSYM
            IF (NNBSTR(ISYM,IRED) .LT. 1) THEN
               DIAMAX(ISYM) = 0.0D0
            ELSE
               DIAMAX(ISYM) = ABS(DIAG(IIBSTR(ISYM,IRED)+1))
            END IF
            AB1 = IIBSTR(ISYM,IRED) + 2
            AB2 = IIBSTR(ISYM,IRED) + NNBSTR(ISYM,IRED)
            DO AB = AB1,AB2
               DIAMAX(ISYM) = MAX(DIAMAX(ISYM),ABS(DIAG(AB)))
            END DO
         END DO
      ELSE IF ((IRED.EQ.2) .OR. (IRED.EQ.3)) THEN
         DO ISYM = 1,NSYM
            IF (NNBSTR(ISYM,IRED) .LT. 1) THEN
               DIAMAX(ISYM) = 0.0D0
            ELSE
               AB = INDRED(IIBSTR(ISYM,IRED)+1,IRED)
               DIAMAX(ISYM) = ABS(DIAG(AB))
            END IF
            AB1 = IIBSTR(ISYM,IRED) + 2
            AB2 = IIBSTR(ISYM,IRED) + NNBSTR(ISYM,IRED)
            DO IAB = AB1,AB2
               AB = INDRED(IAB,IRED)
               DIAMAX(ISYM) = MAX(DIAMAX(ISYM),ABS(DIAG(AB)))
            END DO
         END DO
      ELSE
         WRITE(LUPRI,*) SECNAM,': unknown reduced set, IRED = ',IRED
         CALL CHO_QUIT('Unknown reduced set in '//SECNAM,104)
      END IF

      DGMAX = DIAMAX(1)
      DO ISYM = 2,NSYM
         DGMAX = MAX(DGMAX,DIAMAX(ISYM))
      END DO

      IF (LOCDBG) THEN
         WRITE(LUPRI,*) SECNAM,': in reduced set ',IRED,':'
         WRITE(LUPRI,*) 'DIAMAX = ',(DIAMAX(ISYM),ISYM=1,NSYM)
         WRITE(LUPRI,*) 'DGMAX  = ',DGMAX
      END IF

      END
C  /* Deck cho_setmaxshl */
      SUBROUTINE CHO_SETMAXSHL(DIAG,DIASH,ISYSH,
     &                         INDRED,INDRSH,IRED,MSH,MMBSTRT)
C
C     Purpose: set max. shell pair data for selection procedure.
C
#include <implicit.h>
      DIMENSION DIAG(*), DIASH(MSH,MSH)
      INTEGER   ISYSH(MSH,MSH), INDRED(MMBSTRT,3), INDRSH(MMBSTRT)
#include <cholesky.h>

      CHARACTER*13 SECNAM
      PARAMETER (SECNAM = 'CHO_SETMAXSHL')

      LOGICAL LOW, LOCDBG
      PARAMETER (LOW = .TRUE., LOCDBG = .FALSE.)

      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J

C     Initialize the largest diagonal in each shell pair.
C     ---------------------------------------------------

      CALL CHO_DZERO(DIASH,MSH*MSH)
      CALL CHO_IZERO(ISYSH,MSH*MSH)

C     Find largest diagonal in each shell pair. Loop only
C     over those that are included in the reduced set at hand.
C     --------------------------------------------------------
 
      IF (IRED .EQ. 1) THEN
         DO ISYMAB = 1,NSYM
            IAB1 = IIBSTR(ISYMAB,IRED) + 1
            IAB2 = IIBSTR(ISYMAB,IRED) + NNBSTR(ISYMAB,IRED)
            DO IAB = IAB1,IAB2
               ISHLAB = INDRSH(IAB)  ! get shell pair index
               CALL CHO_INVPCK(ISHLAB,ISHLA,ISHLB,LOW) ! invert index
               IF (DIASH(ISHLA,ISHLB) .GT. DIAG(IAB)) THEN
                  DIASH(ISHLA,ISHLB) = DIAG(IAB)
                  ISYSH(ISHLA,ISHLB) = ISYMAB
               END IF
            END DO
         END DO
      ELSE IF ((IRED.EQ.2) .OR. (IRED.EQ.3)) THEN
         DO ISYMAB = 1,NSYM
            JAB1 = IIBSTR(ISYMAB,IRED) + 1
            JAB2 = IIBSTR(ISYMAB,IRED) + NNBSTR(ISYMAB,IRED)
            DO JAB = JAB1,JAB2
               IAB    = INDRED(JAB,IRED)
               ISHLAB = INDRSH(IAB)
               CALL CHO_INVPCK(ISHLAB,ISHLA,ISHLB,LOW)
               IF (DIASH(ISHLA,ISHLB) .LT. DIAG(IAB)) THEN
                  DIASH(ISHLA,ISHLB) = DIAG(IAB)
                  ISYSH(ISHLA,ISHLB) = ISYMAB
               END IF
            END DO
         END DO
      ELSE
         WRITE(LUPRI,*) SECNAM,': unknown reduced set, IRED = ',IRED
         CALL CHO_QUIT('Unknown reduced set in '//SECNAM,104)
      END IF

C     Ensure symmetric DIASH and ISYSH.
C     ---------------------------------

      IF (LOCDBG) THEN
         WRITE(LUPRI,*) SECNAM,': array DIASHL before symmetrization:'
         CALL CHO_OUTPUT(DIASH,1,MSH,1,MSH,MSH,MSH,1,LUPRI)
      END IF

      IF (LOW) THEN
         DO ISHLA = 1,NSHELL
            DO ISHLB = 1,ISHLA
               DIASH(ISHLB,ISHLA) = DIASH(ISHLA,ISHLB)
               ISYSH(ISHLB,ISHLA) = ISYSH(ISHLA,ISHLB)
            END DO
         END DO
      ELSE
         DO ISHLB = 1,NSHELL
            DO ISHLA = 1,ISHLB
               DIASH(ISHLB,ISHLA) = DIASH(ISHLA,ISHLB)
               ISYSH(ISHLB,ISHLA) = ISYSH(ISHLA,ISHLB)
            END DO
         END DO
      END IF

      IF (LOCDBG) THEN
         WRITE(LUPRI,*) SECNAM,': array DIASHL after symmetrization:'
         CALL CHO_OUTPUT(DIASH,1,MSH,1,MSH,MSH,MSH,1,LUPRI)
      END IF

      END
C  /* Deck cho_getmaxshl */
      SUBROUTINE CHO_GETMAXSHL(DIASH,MSH,SMAX,ISHLA,ISHLB)
C
C     Purpose: Get max. shell pair and update DIASH.
C
#include <implicit.h>
      DIMENSION DIASH(MSH,MSH)
#include <cholesky.h>

      CHARACTER*13 SECNAM
      PARAMETER (SECNAM = 'CHO_GETMAXSHL')

      ISHLA = -1
      ISHLB = -1
      SMAX  = -1.0D8
      DO JSHLA = 1,NSHELL
         DO JSHLB = 1,JSHLA
            IF (DIASH(JSHLA,JSHLB) .GT. SMAX) THEN
               SMAX  = DIASH(JSHLA,JSHLB)
               ISHLA = JSHLA
               ISHLB = JSHLB
            END IF
         END DO
      END DO

      IF ((ISHLA.LT.1) .OR. (ISHLB.LT.1)) THEN
         CALL CHO_QUIT('Error in '//SECNAM,104)
      ELSE
         DIASH(ISHLA,ISHLB) = 0.0D0
      END IF

      END
C  /* Deck cho_getint */
      SUBROUTINE CHO_GETINT(DIAG,DIASH,ISYSH,INDRED,INDRSH,
     &                      IIBSTRSH,NNBSTRSH,
     &                      INTMAP,MSYM,MSH,MMSHL,MMBSTRT,
     &                      WORK,LWORK,NPOTSH,ICOUNT,
     &                      KFRIN,LFRIN)
#include <implicit.h>
      DIMENSION WORK(LWORK)
C
C     Purpose: get qualified integral columns for Cholesky decomposition.
C
C     DIASH(i,j): max. diagonal in shell pair i,j
C     NPOTSH    : the number of shell pairs that can be qualified.
C
#include <cholesky.h>

      DIMENSION DIAG(*), DIASH(MSH,MSH)
      INTEGER   ISYSH(MSH,MSH)
      INTEGER   INDRED(MMBSTRT,3), INDRSH(MMBSTRT)
      INTEGER   IIBSTRSH(MSYM,MMSHL,3), NNBSTRSH(MSYM,MMSHL,3)
      INTEGER   INTMAP(MMSHL)

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_GETINT')

      PARAMETER (ZERO = 0.0D0)

      LOGICAL DODECO, LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER CHO_ISUMELM

      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J

C     Initializations.
C     ----------------

      DO ISYM = 1,NSYM
         IF (LUSEL(ISYM) .GT. 0) REWIND(LUSEL(ISYM))
      END DO

      KFREE = KFRIN
      LFREE = LFRIN

      CALL CHO_IZERO(NQUAL,NSYM)
      ICOUNT = 0
      IF (MXSHPR .GT. 0) THEN
         MCOUNT = MIN(NPOTSH,MXSHPR)
      ELSE
         MCOUNT = NPOTSH
      END IF
      DODECO = .FALSE.

C     Shell pair qualification loop.
C     ------------------------------

      DO WHILE ((.NOT.DODECO) .AND. (ICOUNT.LT.MCOUNT))

C        Update shell pair counter.
C        --------------------------

         ICOUNT = ICOUNT + 1

C        Get shell pair corresponding to largest diagonal.
C        -------------------------------------------------

         CALL CHO_GETMAXSHL(DIASH,NSHELL,SMAX,ISHLA,ISHLB)
         ISYMAB = ISYSH(ISHLA,ISHLB)
         ISHLAB = INDEX(ISHLA,ISHLB)

         IF ((SMAX.EQ.ZERO) .OR. (ABS(SMAX).LT.DIAMIN(ISYMAB))) THEN

C           Diagonal too small to be qualified for decomposition.
C           -----------------------------------------------------

            IF (ICOUNT .EQ. 1) THEN
               WRITE(LUPRI,*) SECNAM,': no integrals calculated; ',
     &                        'unable to proceed to decomposition!'
               WRITE(LUPRI,*) 'Max. abs. diagonal for shell pair ',
     &                        ISHLA,', ',ISHLB,': ',ABS(SMAX)
               WRITE(LUPRI,*) 'Max. abs. diagonal allowed: ',
     &                        DIAMIN(ISYMAB),' (sym. ',ISYMAB,')'
               CALL CHO_QUIT('Fatal error in '//SECNAM,104)
            ELSE
               ICOUNT = ICOUNT - 1
               NSEL   = CHO_ISUMELM(NQUAL,NSYM)
               DODECO = NSEL .GT. 0
            END IF

         ELSE

C           Qualify diagonals within this shell pair.
C           -----------------------------------------

            CALL CHO_QUALIFY(DIAG,ISHLA,ISHLB,
     &                       IIBSTRSH,NNBSTRSH,INDRED,INDRSH,
     &                       MSYM,MMBSTRT,MMSHL)

C           Calculate integral columns; get qualified ones stored in
C           current reduced set; write these to disk on temporary
C           file(s).
C           --------------------------------------------------------

            NSEL   = CHO_ISUMELM(NQUAL,NSYM)
            NCOLAB = NSEL - CHO_ISUMELM(IOFFQ,NSYM)

            IF (NCOLAB .GT. 0) THEN

               INTMAP(ISHLAB) = INTMAP(ISHLAB) + 1
               WRITE(LUPRI,'(/,A,I5,1X,I5,A,I9,A)')
     &         'Calculating shell pair (**|',ISHLA,ISHLB,
     &         '):',NCOLAB,' columns have been qualified'
               WRITE(LUPRI,'(80A)') ('=',i=1,77)
               WRITE(LUPRI,'(A,I12)')
     &         'Number of calculations so far for this shell pair: ',
     &         INTMAP(ISHLAB)
               CALL CHO_FLUSH(LUPRI)

               CALL CHO_DAL_CALCINT(IIBSTRSH,NNBSTRSH,INDRED,INDRSH,
     &                              MSYM,MMSHL,MMBSTRT,ISHLA,ISHLB,
     &                              WORK,LWORK)

C              Enough integral columns for proceeding to decomposition?
C              --------------------------------------------------------

               DODECO = NSEL .GE. MINQUAL

            ELSE IF (NCOLAB .EQ. 0) THEN

               IF (NSEL .LT. 1) THEN
                  WRITE(LUPRI,*) SECNAM,': logical error: ',
     &                                  'unable to qualify diagonals'
                  WRITE(LUPRI,*) SECNAM,': NCOLAB = ',NCOLAB
                  WRITE(LUPRI,*) SECNAM,': NSEL   = ',NSEL
                  CALL CHO_QUIT('[0] Logical error in '//SECNAM,104)
               ELSE
                  ICOUNT = ICOUNT - 1
                  DODECO = .TRUE.
               END IF

            ELSE

               WRITE(LUPRI,*) SECNAM,': logical error: ',
     &                               'unable to qualify diagonals'
               WRITE(LUPRI,*) SECNAM,': NCOLAB = ',NCOLAB
               WRITE(LUPRI,*) SECNAM,': NSEL   = ',NSEL
               CALL CHO_QUIT('[1] Logical error in '//SECNAM,104)

            END IF

         END IF

      END DO

C     Test loop exit (we may have calculated all possible integrals, yet
C     NSEL < MINQUAL).
C     ------------------------------------------------------------------

      IF (.NOT. DODECO) THEN
         NSEL = CHO_ISUMELM(NQUAL,NSYM)
         IF (NSEL .LT. 1) THEN
            WRITE(LUPRI,*) SECNAM,': logical error: ',
     &                            'unable to qualify diagonals'
            WRITE(LUPRI,*) SECNAM,': Flag DODECO is ',DODECO
            WRITE(LUPRI,*) SECNAM,': NSEL    = ',NSEL
            WRITE(LUPRI,*) SECNAM,': ICOUNT  = ',ICOUNT
            WRITE(LUPRI,*) SECNAM,': MCOUNT  = ',MCOUNT
            WRITE(LUPRI,*) SECNAM,': NPOTSH  = ',NPOTSH
            WRITE(LUPRI,*) SECNAM,': MINQUAL = ',MINQUAL
            CALL CHO_QUIT('[2] Logical error in '//SECNAM,103)
         ELSE
            DODECO = .TRUE.
         END IF
      END IF

      END
C  /* Deck cho_setred */
      SUBROUTINE CHO_SETRED(DIAG,IIBSTRSH,NNBSTRSH,INDRED,INDRSH,
     &                      MSYM,MMBSTRT,MMSHL)
C
C     Purpose: set next reduced set. A copy of the previous set
C              is stored in IRED=3.
C
#include <implicit.h>
      DIMENSION DIAG(*)
      INTEGER   INDRED(MMBSTRT,3), INDRSH(MMBSTRT)
      INTEGER   IIBSTRSH(MSYM,MMSHL,3), NNBSTRSH(MSYM,MMSHL,3)
#include <cholesky.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_SETRED')

C     Debug print.
C     ------------

      IF (CHO_TRCNEG) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) SECNAM,
     &                  ': tracing of negative diagonals activated.'
         WRITE(LUPRI,*) SECNAM,': flag CHO_USEABS is ',CHO_USEABS
         WRITE(LUPRI,*) SECNAM,': checking for negative diagonals in ',
     &                  'first reduced set:'
         NNEG = 0
         DO ISYM = 1,NSYM
            JAB1 = IIBSTR(ISYM,1) + 1
            JAB2 = JAB1 + NNBSTR(ISYM,1) - 1
            INEG = 0
            DO JAB = JAB1,JAB2
               IF (DIAG(JAB) .LT. 0.0D0) THEN
                  INEG = INEG + 1
               END IF
            END DO
            NNEG = NNEG + INEG
            WRITE(LUPRI,*) SECNAM,': #negative in symmetry ',ISYM,
     &                     ': ',INEG
         END DO
         WRITE(LUPRI,*) SECNAM,': total #negative: ',NNEG
      END IF

C     Set number of shell pairs.
C     --------------------------

      NNSHL = NSHELL*(NSHELL + 1)/2

C     Copy index arrays.
C     ------------------

      CALL ICOPY(MMBSTRT,INDRED(1,2),1,INDRED(1,3),1)
      CALL ICOPY(MSYM*MMSHL,IIBSTRSH(1,1,2),1,IIBSTRSH(1,1,3),1)
      CALL ICOPY(MSYM*MMSHL,NNBSTRSH(1,1,2),1,NNBSTRSH(1,1,3),1)
      CALL ICOPY(MSYM,IIBSTR(1,2),1,IIBSTR(1,3),1)
      CALL ICOPY(MSYM,NNBSTR(1,2),1,NNBSTR(1,3),1)
      NNBSTRT(3) = NNBSTRT(2)

C     Re-initialize IRED=2.
C     ---------------------

      CALL CHO_IZERO(INDRED(1,2),MMBSTRT)
      CALL CHO_IZERO(IIBSTRSH(1,1,2),MSYM*MMSHL)
      CALL CHO_IZERO(NNBSTRSH(1,1,2),MSYM*MMSHL)
      CALL CHO_IZERO(IIBSTR(1,2),MSYM)
      CALL CHO_IZERO(NNBSTR(1,2),MSYM)
      NNBSTRT(2) = 0

C     Set new reduced set: mapping array and shell-pair count.
C     --------------------------------------------------------

      IF (SCDIAG) THEN  ! do screening

         IF (MODE_SCREEN .EQ. 1) THEN ! both conv. and unconv. included

            IF (CHO_USEABS) THEN ! neg. diag. might be included

               KAB = 0
               DO ISYM = 1,NSYM
                  IF (NNBSTR(ISYM,3) .GT. 0) THEN

                     JAB1 = IIBSTR(ISYM,3) + 1
                     JAB2 = JAB1 + NNBSTR(ISYM,3) - 1

                     IAB = INDRED(JAB1,3)
                     XM  = ABS(DIAG(IAB))
                     DO JAB = JAB1+1,JAB2
                        IAB = INDRED(JAB,3)
                        XM  = MAX(XM,ABS(DIAG(IAB)))
                     END DO

                     IF (XM .GT. THRCOM) THEN  ! only if not converged
                        DO JAB = JAB1,JAB2
                           IAB = INDRED(JAB,3)
                           TST = SQRT(ABS(DIAG(IAB))*XM)*DAMP(2)
                           IF (TST .GT. THRCOM) THEN
                              ISHLAB = INDRSH(IAB)
                              KAB    = KAB + 1
                              INDRED(KAB,2) = IAB
                              NNBSTRSH(ISYM,ISHLAB,2) =
     &                                       NNBSTRSH(ISYM,ISHLAB,2) + 1
                           END IF
                        END DO
                     END IF

                  END IF
               END DO

            ELSE ! neg. diag. excluded

               KAB = 0
               DO ISYM = 1,NSYM
                  IF (NNBSTR(ISYM,3) .GT. 0) THEN

                     JAB1 = IIBSTR(ISYM,3) + 1
                     JAB2 = JAB1 + NNBSTR(ISYM,3) - 1

                     IAB = INDRED(JAB1,3)
                     XM  = ABS(DIAG(IAB))
                     DO JAB = JAB1+1,JAB2
                        IAB = INDRED(JAB,3)
                        XM  = MAX(XM,ABS(DIAG(IAB)))
                     END DO

                     IF (XM .GT. THRCOM) THEN  ! only if not converged
                        DO JAB = JAB1,JAB2
                           IAB = INDRED(JAB,3)
                           IF (DIAG(IAB) .GT. 0.0D0) THEN ! neg => conv
                              TST = SQRT(DIAG(IAB)*XM)*DAMP(2)
                              IF (TST .GT. THRCOM) THEN
                                 ISHLAB = INDRSH(IAB)
                                 KAB    = KAB + 1
                                 INDRED(KAB,2) = IAB
                                 NNBSTRSH(ISYM,ISHLAB,2) =
     &                                       NNBSTRSH(ISYM,ISHLAB,2) + 1
                              END IF
                           END IF
                        END DO
                     END IF

                  END IF
               END DO

            END IF

         ELSE IF (MODE_SCREEN .EQ. 2) THEN ! only unconv. included

            IF (CHO_USEABS) THEN ! neg. diag. might be included

               KAB = 0
               DO ISYM = 1,NSYM
                  IF (NNBSTR(ISYM,3) .GT. 0) THEN

                     JAB1 = IIBSTR(ISYM,3) + 1
                     JAB2 = JAB1 + NNBSTR(ISYM,3) - 1

                     DO JAB = JAB1,JAB2
                        IAB = INDRED(JAB,3)
                        IF (ABS(DIAG(IAB)) .GT. THRCOM) THEN
                           ISHLAB = INDRSH(IAB)
                           KAB    = KAB + 1
                           INDRED(KAB,2) = IAB
                           NNBSTRSH(ISYM,ISHLAB,2) =
     &                                       NNBSTRSH(ISYM,ISHLAB,2) + 1
                        END IF
                     END DO

                  END IF

               END DO

            ELSE ! neg. diag. excluded

               KAB = 0
               DO ISYM = 1,NSYM
                  IF (NNBSTR(ISYM,3) .GT. 0) THEN

                     JAB1 = IIBSTR(ISYM,3) + 1
                     JAB2 = JAB1 + NNBSTR(ISYM,3) - 1

                     DO JAB = JAB1,JAB2
                        IAB = INDRED(JAB,3)
                        IF (DIAG(IAB) .GT. THRCOM) THEN
                           ISHLAB = INDRSH(IAB)
                           KAB    = KAB + 1
                           INDRED(KAB,2) = IAB
                           NNBSTRSH(ISYM,ISHLAB,2) =
     &                                       NNBSTRSH(ISYM,ISHLAB,2) + 1
                        END IF
                     END DO

                  END IF
               END DO

            END IF

         ELSE ! MODE_SCREEN out of bounds

            WRITE(LUPRI,*) SECNAM,': MODE_SCREEN = ',MODE_SCREEN
            CALL CHO_QUIT('MODE_SCREEN out of bounds in '//SECNAM,103)

         END IF

      ELSE ! no screening; remove zero diagonals and check convergence

         IF (CHO_USEABS) THEN

            KAB = 0
            DO ISYM = 1,NSYM
               IF (NNBSTR(ISYM,3) .GT. 0) THEN

                  JAB1 = IIBSTR(ISYM,3) + 1
                  JAB2 = JAB1 + NNBSTR(ISYM,3) - 1

                  IAB = INDRED(JAB1,3)
                  XM  = ABS(DIAG(IAB))
                  DO JAB = JAB1+1,JAB2
                     IAB = INDRED(JAB,3)
                     XM  = MAX(XM,ABS(DIAG(IAB)))
                  END DO

                  IF (XM .GT. THRCOM) THEN  ! only if not converged
                     DO JAB = JAB1,JAB2
                        IAB = INDRED(JAB,3)
                        IF (ABS(DIAG(IAB)) .NE. 0.0D0) THEN
                           ISHLAB = INDRSH(IAB)
                           KAB    = KAB + 1
                           INDRED(KAB,2) = IAB
                           NNBSTRSH(ISYM,ISHLAB,2) =
     &                                       NNBSTRSH(ISYM,ISHLAB,2) + 1
                        END IF
                     END DO
                  END IF

               END IF
            END DO

         ELSE

            KAB = 0
            DO ISYM = 1,NSYM
               IF (NNBSTR(ISYM,3) .GT. 0) THEN

                  JAB1 = IIBSTR(ISYM,3) + 1
                  JAB2 = JAB1 + NNBSTR(ISYM,3) - 1

                  IAB = INDRED(JAB1,3)
                  XM  = ABS(DIAG(IAB))
                  DO JAB = JAB1+1,JAB2
                     IAB = INDRED(JAB,3)
                     XM  = MAX(XM,ABS(DIAG(IAB)))
                  END DO

                  IF (XM .GT. THRCOM) THEN  ! only if not converged
                     DO JAB = JAB1,JAB2
                        IAB = INDRED(JAB,3)
                        IF (DIAG(IAB) .NE. 0.0D0) THEN
                           ISHLAB = INDRSH(IAB)
                           KAB    = KAB + 1
                           INDRED(KAB,2) = IAB
                           NNBSTRSH(ISYM,ISHLAB,2) =
     &                                       NNBSTRSH(ISYM,ISHLAB,2) + 1
                        END IF
                     END DO
                  END IF

               END IF
            END DO

         END IF

      END IF

C     Construct remaining reduced set index arrays.
C     ---------------------------------------------

      CALL CHO_SETREDIND(IIBSTRSH,NNBSTRSH,MSYM,MMSHL,2)

C     Debug print.
C     ------------

      IF (CHO_TRCNEG) THEN
         WRITE(LUPRI,*) SECNAM,': checking for negative diagonals in ',
     &                  'next reduced set:'
         NNEG = 0
         DO ISYM = 1,NSYM
            JAB1 = IIBSTR(ISYM,2) + 1
            JAB2 = JAB1 + NNBSTR(ISYM,2) - 1
            INEG = 0
            DO JAB = JAB1,JAB2
               IAB = INDRED(JAB,2)
               IF (DIAG(IAB) .LT. 0.0D0) THEN
                  INEG = INEG + 1
               END IF
            END DO
            NNEG = NNEG + INEG
            WRITE(LUPRI,*) SECNAM,': #negative in symmetry ',ISYM,
     &                     ': ',INEG
         END DO
         WRITE(LUPRI,*) SECNAM,': total #negative: ',NNEG
      END IF

      END
C  /* Deck cho_setredind */
      SUBROUTINE CHO_SETREDIND(IIBSTRSH,NNBSTRSH,MSYM,MMSHL,IRED)
C
C     Purpose: set index arrays for reduced set IRED. The counter
C              array NNBSTRSH must be set on entry.
C
#include <implicit.h>
      INTEGER IIBSTRSH(MSYM,MMSHL,3), NNBSTRSH(MSYM,MMSHL,3)
#include <cholesky.h>

      CHARACTER*13 SECNAM
      PARAMETER (SECNAM = 'CHO_SETREDIND')

      NNSHL = NSHELL*(NSHELL + 1)/2
      IF ((NNSHL.NE.MMSHL) .OR. (NSYM.NE.MSYM))
     & CALL CHO_QUIT('[1] Dimension error in '//SECNAM,104)

      J = IRED
      IF ((J.LT.1) .OR. (J.GT.3))
     & CALL CHO_QUIT('[2] Dimension error in '//SECNAM,104)

      NNBSTRT(J) = 0
      DO ISYM = 1,NSYM
         IIBSTRSH(ISYM,1,J) = 0
         NNBSTR(ISYM,J) = NNBSTRSH(ISYM,1,J)
         DO ISHLAB = 2,NNSHL
            IIBSTRSH(ISYM,ISHLAB,J) = NNBSTR(ISYM,J)
            NNBSTR(ISYM,J) = NNBSTR(ISYM,J) + NNBSTRSH(ISYM,ISHLAB,J)
         END DO
         IIBSTR(ISYM,J) = NNBSTRT(J)
         NNBSTRT(J) = NNBSTRT(J) + NNBSTR(ISYM,J)
      END DO

      END
C  /* Deck cho_qualify */
      SUBROUTINE CHO_QUALIFY(DIAG,ISHLA,ISHLB,
     &                       IIBSTRSH,NNBSTRSH,INDRED,INDRSH,
     &                       MSYM,MMBSTRT,MMSHL)
C
C     Purpose: qualify diagonal elements for decomposition in
C              current reduced set.
C
#include <implicit.h>
      DIMENSION DIAG(*)
      INTEGER   INDRED(MMBSTRT,3), INDRSH(MMBSTRT)
      INTEGER   IIBSTRSH(MSYM,MMSHL,3), NNBSTRSH(MSYM,MMSHL,3)
#include <cholesky.h>

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CHO_QUALIFY')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J

C     Copy counter to offset array.
C     -----------------------------

      CALL ICOPY(NSYM,NQUAL,1,IOFFQ,1)

C     Get shell pair index.
C     ---------------------

      ISHLAB = INDEX(ISHLA,ISHLB)

C     Qualify.
C     --------

      IF (IALQUA .EQ. 1) THEN  ! qualify until full
         DO ISYM = 1,NSYM
            I  = IIBSTR(ISYM,2) + IIBSTRSH(ISYM,ISHLAB,2)
            I2 = I + NNBSTRSH(ISYM,ISHLAB,2)
            DO WHILE ((I.LT.I2) .AND. (NQUAL(ISYM).LT.MAXQUAL))
               I = I + 1
               J = INDRED(I,2)
               IF (DIAG(J) .GE. DIAMIN(ISYM)) THEN
                  NQUAL(ISYM) = NQUAL(ISYM) + 1
                  IQUAB(NQUAL(ISYM),ISYM) = I
               END IF
            END DO
         END DO
      ELSE ! qualify until full, then largest
         DO ISYM = 1,NSYM
            NDIM = NNBSTRSH(ISYM,ISHLAB,2)
            NMAX = MAXQUAL - NQUAL(ISYM)
            IF ((NDIM.GT.0) .AND. (NMAX.GT.0)) THEN
               I1 = IIBSTR(ISYM,2) + IIBSTRSH(ISYM,ISHLAB,2) + 1
               I2 = I1 + NDIM - 1
               IF (NMAX .EQ. 1) THEN ! qualify the largest > DIAMIN
                  XMAX = DIAMIN(ISYM)
                  IMAX = -1
                  DO I = I1,I2
                     J = INDRED(I,2)
                     IF (DIAG(J) .GE. XMAX) THEN
                        XMAX = DIAG(J)
                        IMAX = I
                     END IF
                  END DO
                  IF (IMAX .GT. 0) THEN
                     NQUAL(ISYM) = NQUAL(ISYM) + 1
                     IQUAB(NQUAL(ISYM),ISYM) = IMAX
                  END IF
               ELSE ! full search
                  DO I = I1,I2
                     J = INDRED(I,2)
                     IF (DIAG(J) .GE. DIAMIN(ISYM)) THEN
                        IF (NQUAL(ISYM) .LT. MAXQUAL) THEN
                           NQUAL(ISYM) = NQUAL(ISYM) + 1
                           IQUAB(NQUAL(ISYM),ISYM) = I
                        ELSE IF (NQUAL(ISYM) .EQ. MAXQUAL) THEN
                           K1  = IOFFQ(ISYM) + 1
                           K2  = NQUAL(ISYM)
                           II1 = IQUAB(K1,ISYM)
                           JJ1 = INDRED(II1,2)
                           XMIN = DIAG(JJ1)
                           KKMN = K1
                           DO K = K1+1,K2 ! find min. among qualified
                              II = IQUAB(K,ISYM)
                              JJ = INDRED(II,2)
                              IF (DIAG(JJ) .LT. XMIN) THEN
                                 XMIN = DIAG(JJ)
                                 KKMN = K
                              END IF
                           END DO
                           IF (DIAG(J) .GT. XMIN) THEN ! replace
                              IQUAB(KKMN,ISYM) = I
                           END IF
                        ELSE
                           CALL CHO_QUIT('Logical error in '//SECNAM,
     &                                   104)
                        END IF
                     END IF
                  END DO
               END IF
            END IF
         END DO
      END IF

C     Debug: print.
C     -------------

      IF (LOCDBG) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) SECNAM,': qualified diagonals from shell-pair ',
     &                  ISHLA,ISHLB,':'
         WRITE(LUPRI,*) SECNAM,': algorithm: ',IALQUA
         DO ISYM = 1,NSYM
            NUM = NQUAL(ISYM) - IOFFQ(ISYM)
            WRITE(LUPRI,*)
            WRITE(LUPRI,*) 'Sym.,dimension,#qualified,threshold: ',
     &                     ISYM,NNBSTRSH(ISYM,ISHLAB,2),NUM,DIAMIN(ISYM)
            IF (NNBSTRSH(ISYM,ISHLAB,2) .GT. 0) THEN
               I1 = IIBSTR(ISYM,2) + IIBSTRSH(ISYM,ISHLAB,2) + 1
               I2 = I1 + NNBSTRSH(ISYM,ISHLAB,2) - 1
               WRITE(LUPRI,*) 'Diagonal (current reduced set):'
               WRITE(LUPRI,*) 'Index range, total: ',I1,I2
               WRITE(LUPRI,*) 'Index range, sym. : ',I1-IIBSTR(ISYM,2),
     &                                               I2-IIBSTR(ISYM,2)
               WRITE(LUPRI,'(5F15.8)') (DIAG(INDRED(I,2)), I=I1,I2)
               K1 = IOFFQ(ISYM) + 1
               K2 = NQUAL(ISYM)
               WRITE(LUPRI,*) 'Qualified diagonals:'
               WRITE(LUPRI,'(5F15.8)') (DIAG(INDRED(IQUAB(K,ISYM),2)),
     &                                  K = K1,K2)
            END IF
         END DO
         WRITE(LUPRI,*)
         WRITE(LUPRI,*)
      END IF

      END
C  /* Deck cho_qualify_old */
      SUBROUTINE CHO_QUALIFY_OLD(DIAG,ISHLA,ISHLB,ISYMAB,
     &                       IIBSTRSH,NNBSTRSH,INDRED,INDRSH,
     &                       MSYM,MMBSTRT,MMSHL)
C
C     Purpose: qualify diagonal elements for decomposition in
C              current reduced set (2). ISYMAB is the symmetry of the
C              largest diagonal in shell pair ISHLA,ISHLB.
C
#include <implicit.h>
      DIMENSION DIAG(*)
      INTEGER   INDRED(MMBSTRT,3), INDRSH(MMBSTRT)
      INTEGER   IIBSTRSH(MSYM,MMSHL,3), NNBSTRSH(MSYM,MMSHL,3)
#include <cholesky.h>

      CHARACTER*15 SECNAM
      PARAMETER (SECNAM = 'CHO_QUALIFY_OLD')

      INTEGER CHO_ISUMELM

      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J

C     Copy counter to offset array.
C     -----------------------------

      CALL ICOPY(NSYM,NQUAL,1,IOFFQ,1)

C     Get shell pair index.
C     ---------------------

      ISHLAB = INDEX(ISHLA,ISHLB)

C     To ensure that the largest diagonal in this shell pair
C     is indeed qualified, we start with symmetry ISYMAB.
C     ------------------------------------------------------

      ISYM = ISYMAB

C     If more diagonals are needed, qualify among the remaining
C     symmetries.
C     ---------------------------------------------------------

      NSEL = CHO_ISUMELM(NQUAL,NSYM)
      IF (NSEL .LT. MINQUAL) THEN
         DO ISYM = 1,NSYM
            IF (ISYM .NE. ISYMAB) THEN

               IAB  = IIBSTR(ISYM,2) + IIBSTRSH(ISYM,ISHLAB,2)
               IAB2 = IAB + NNBSTRSH(ISYM,ISHLAB,2)

               DO WHILE ((IAB.LT.IAB2).AND.(NQUAL(ISYM).LT.MAXQUAL))
                  IAB = IAB + 1
                  JAB = INDRED(IAB,2) ! get index in original red. set
                  IF (DIAG(JAB) .GE. DIAMIN(ISYM)) THEN
                     NQUAL(ISYM) = NQUAL(ISYM) + 1
                     IQUAB(NQUAL(ISYM),ISYM) = IAB
                  END IF
               END DO
 
            END IF
         END DO
      END IF

      END
C  /* Deck cho_decom */
      SUBROUTINE CHO_DECOM(DIAG,INFRED,INFVEC,INDRED,INDRSH,
     &                     IIBSTRSH,NNBSTRSH,ISCR,
     &                     MRED,MVEC,MSYM,MMSHL,MMBSTRT,LISCR,
     &                     WORK,LWORK,IPASS,NUM)
C
C     Purpose: calculate Cholesky vectors from qualified integral
C              columns (from disk).
C
#include <implicit.h>
      DIMENSION DIAG(*), WORK(LWORK)
      INTEGER   INFRED(MRED), INFVEC(MVEC,3,MSYM)
      INTEGER   INDRED(MMBSTRT,3), INDRSH(MMBSTRT)
      INTEGER   IIBSTRSH(MSYM,MMSHL,3), NNBSTRSH(MSYM,MMSHL,3)
      INTEGER   ISCR(LISCR)
#include <cholesky.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CHO_DECOM')

      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)

      LOGICAL LAST

      LENLIN = 0  ! to avoid compiler warnings...
      IF (IPRINT .GE. INF_PROGRESS) THEN
         CALL CHO_HEAD(SECNAM//
     &                 ': Decomposition of Qualified Diagonals','=',
     &                 80,LUPRI)
         WRITE(LUPRI,'(/,A,I5,A,I4,A)')
     &   'Integral pass number',IPASS,' (',NUM,
     &   ' shell pair distributions calculated)'
         WRITE(LUPRI,'(A,8I8)')
     &   '#Cholesky vec. : ',(NUMCHO(ISYM),ISYM=1,NSYM)
         WRITE(LUPRI,'(A,8I8)')
     &   '#qualifieds    : ',(NQUAL(ISYM),ISYM=1,NSYM)
         WRITE(LUPRI,'(A,8I8)')
     &   'Red.  dimension: ',(NNBSTR(ISYM,2),ISYM=1,NSYM)
         WRITE(LUPRI,'(A,8I8)')
     &   'Orig. dimension: ',(NNBSTR(ISYM,1),ISYM=1,NSYM)
         WRITE(LUPRI,'(/,A,/,A,A)')
     &   '           #Vectors             Treated Diagonal',
     &   'Sym.     Sym.     Total     Index     Before      After',
     &   '   Conv. Neg.   New Max'
         LENLIN = 79
         WRITE(LUPRI,'(80A)') ('-',I=1,LENLIN)
         CALL CHO_FLUSH(LUPRI)
      END IF

      DO ISYM = 1,NSYM

C        Cycle loop if nothing to do in this symmetry.
C        ---------------------------------------------

         IF ((NQUAL(ISYM).LT.1) .OR. (NNBSTR(ISYM,2).LT.1)) GO TO 100

C        Reserve space for qualified integral columns.
C        ---------------------------------------------

         LINT1 = NNBSTR(ISYM,2)*NQUAL(ISYM) ! integrals

         KINT1 = 1
         KEND0 = KINT1 + LINT1
         LWRK0 = LWORK - KEND0 + 1
         IF (LWRK0 .LE. 0) THEN
            CALL CHO_QUIT('[0] Insufficient memory in '//SECNAM,101)
         END IF

C        Determine size of Cholesky vector buffer.
C        -----------------------------------------

         NUMBUF = MIN(LWRK0/NNBSTR(ISYM,2),NQUAL(ISYM))
         IF (NUMBUF .LT. 1) THEN
            CALL CHO_QUIT('[1] Insufficient memory in '//SECNAM,101)
         END IF
         KCHO1 = KEND0
         KEND1 = KCHO1 + NNBSTR(ISYM,2)*NUMBUF
         LWRK1 = LWORK - KEND1 + 1
         IF (LWRK1 .LT. 0) THEN  ! should be redundant...
            CALL CHO_QUIT('Buffer allocation error in '//SECNAM,101)
         END IF

C        Read qualified integral columns.
C        --------------------------------

         CALL CHO_TIMER(C1,W1)
         REWIND(LUSEL(ISYM))
         DO I = 1,NQUAL(ISYM)
            KOFF = KINT1 + NNBSTR(ISYM,2)*(I - 1) - 1
            READ(LUSEL(ISYM)) (WORK(KOFF+J), J=1,NNBSTR(ISYM,2))
         END DO
         CALL CHO_TIMER(C2,W2)
         TDECOM(1,1) = TDECOM(1,1) + C2 - C1
         TDECOM(2,1) = TDECOM(2,1) + W2 - W1

C        Subtract contributions from previous vectors.
C        ---------------------------------------------

         CALL CHO_SUBTR(WORK(KINT1),INFRED,INFVEC,
     &                  INDRED,INDRSH,IIBSTRSH,NNBSTRSH,ISCR,
     &                  MRED,MVEC,MSYM,MMSHL,MMBSTRT,LISCR,
     &                  WORK(KEND0),LWRK0,ISYM)

C        Decompose in loop over qualified columns.
C        -----------------------------------------

         IVEC  = NUMCHO(ISYM)
         IVECT = NUMCHT
         IDUMP = 0
         DO ICHO = 1,NQUAL(ISYM)

C           Find max. diagonal among qualified.
C           -----------------------------------

            IAB  = 1
            IABG = INDRED(IQUAB(IAB,ISYM),2)
            XC   = DIAG(IABG)
            DO I = 2,NQUAL(ISYM)
               KAB = INDRED(IQUAB(I,ISYM),2)
               IF (DIAG(KAB) .GT. XC) THEN
                  IAB  = I
                  IABG = KAB
                  XC   = DIAG(KAB)
               END IF
            END DO

C           Decompose if max. diagonal is still qualified.
C           ----------------------------------------------

            LAST = (XC.LT.DIAMIN(ISYM)) .OR. (XC.LT.THRCOM)
            IF (.NOT. LAST) THEN

C              Offset to max. diagonal column.
C              -------------------------------

               KOFF0 = KINT1 + NNBSTR(ISYM,2)*(IAB - 1) - 1

C              Scale column corresponding to max. diagonal to obtain
C              the Cholesky vector.
C              -----------------------------------------------------

               FAC  = ONE/SQRT(XC)
               KOFF = KOFF0 + 1
               CALL DSCAL(NNBSTR(ISYM,2),FAC,WORK(KOFF),1)

C              Zero entries in Cholesky vector corresponding to zero
C              diagonals.
C              -----------------------------------------------------

               DO I = 1,NNBSTR(ISYM,2)
                  II = IIBSTR(ISYM,2) + I
                  JJ = INDRED(II,2)
                  IF (DIAG(JJ) .EQ. ZERO) THEN
                     KOFF = KOFF0 + I
                     WORK(KOFF) = ZERO
                  END IF
               END DO

C              Update diagonal.
C              ----------------

               DO I = 1,NNBSTR(ISYM,2)
                  II   = IIBSTR(ISYM,2) + I
                  JJ   = INDRED(II,2)
                  KOFF = KOFF0 + I
                  DIAG(JJ) = DIAG(JJ) - WORK(KOFF)*WORK(KOFF)
               END DO

C              Zero treated diagonal element and analyze updated diagonal.
C              -----------------------------------------------------------

               OLDIAG     = DIAG(IABG)
               DIAG(IABG) = ZERO
               CALL CHO_CHKDIA2(DIAG,INDRED,INDRSH,MMBSTRT,ISYM,
     &                          XMIN,XMAX,XM,NNEGT,NNEG,NCONV)

C              Update total number of zeroed negative diagonals.
C              -------------------------------------------------

               NNZTOT = NNZTOT + NNEG

C              Update DIAMIN from max. abs. diagonal element XM.
C              -------------------------------------------------

               DIAMIN(ISYM) = MAX(XM*SPAN,THRCOM)

C              Subtract this Cholesky vector from integrals. If
C              the corresponding diagonal element is zero, the
C              column will no longer be qualified and subtraction
C              can safely be skipped.
C              --------------------------------------------------

               KOFF1 = KOFF0 + 1
               DO I = 1,NQUAL(ISYM)
                  II = IQUAB(I,ISYM)
                  JJ = INDRED(II,2)
                  IF (DIAG(JJ) .NE. ZERO) THEN
                     KOFF2 = KINT1 + NNBSTR(ISYM,2)*(I - 1)
                     KOFF3 = KOFF0 + II - IIBSTR(ISYM,2)
                     FAC   = -WORK(KOFF3)
                     CALL DAXPY(NNBSTR(ISYM,2),FAC,WORK(KOFF1),1,
     &                                             WORK(KOFF2),1)
                  END IF
               END DO

C              Store Cholesky vector in buffer.
C              --------------------------------

               IDUMP = IDUMP + 1

               KOFF1 = KOFF0 + 1
               KOFF2 = KCHO1 + NNBSTR(ISYM,2)*(IDUMP - 1)
               CALL DCOPY(NNBSTR(ISYM,2),WORK(KOFF1),1,WORK(KOFF2),1)

C              Update Cholesky vector counters.
C              --------------------------------

               IVEC  = IVEC  + 1
               IVECT = IVECT + 1

C              Set info for this vector.
C              -------------------------

               IF (IVEC .GT. MAXVEC) THEN
                  WRITE(LUPRI,*) SECNAM,': too many Cholesky vectors!'
                  WRITE(LUPRI,*) SECNAM,': symmetry: ',ISYM
                  WRITE(LUPRI,*) SECNAM,': max. allowed is ',MAXVEC
                  WRITE(LUPRI,*) SECNAM,': please increase max. ',
     &                           'allowed'
                  CALL CHO_QUIT('Too many Cholesky vectors in '
     &                          //SECNAM,104)
               ELSE IF (IVEC .EQ. MAXVEC) THEN ! no set next disk addr.
                  INFVEC(IVEC,1,ISYM) = IABG   ! diag. index red. set 1
                  INFVEC(IVEC,2,ISYM) = IPASS  ! global red. set
               ELSE
                  INFVEC(IVEC,1,ISYM)   = IABG  ! diag. index red. set 1
                  INFVEC(IVEC,2,ISYM)   = IPASS ! global red. set
                  INFVEC(IVEC+1,3,ISYM) = INFVEC(IVEC,3,ISYM)
     &                                  + NNBSTR(ISYM,2) ! next addr.
               END IF

C              Print progress report.
C              ----------------------

               IF (IPRINT .GE. INF_PROGRESS) THEN
              WRITE(LUPRI,'(I3,3(1X,I9),2(1X,D11.3),2(1X,I4),1X,D11.3)')
     &            ISYM,IVEC,IVECT,IABG,XC,OLDIAG,NCONV,NNEG,XM
               END IF

            END IF

C           Dump vectors to disk when there is no more to be done, or
C           when the buffer is full.
C           ---------------------------------------------------------

            IF (LAST .OR. (IDUMP.EQ.NUMBUF)) THEN
               CALL CHO_TIMER(C1,W1)
               IVEC1 = NUMCHO(ISYM) + 1
               CALL CHO_PUTVEC(WORK(KCHO1),NNBSTR(ISYM,2),IDUMP,IVEC1,
     &                         ISYM,INFVEC,MVEC,MSYM)
               NUMCHO(ISYM) = NUMCHO(ISYM) + IDUMP
               NUMCHT       = NUMCHT       + IDUMP
               CALL CHO_TIMER(C2,W2)
               TDECOM(1,2) = TDECOM(1,2) + C2 - C1
               TDECOM(2,2) = TDECOM(2,2) + W2 - W1
               IF (LAST) THEN
                  GO TO 100   ! cycle symmetry loop
               ELSE
                  IVEC  = NUMCHO(ISYM)
                  IVECT = NUMCHT
                  IDUMP = 0
               END IF
            END IF

         END DO

C        Cycle point: go to next symmetry.
C        ---------------------------------

  100    CONTINUE
         IF (IPRINT .GE. INF_PROGRESS) CALL CHO_FLUSH(LUPRI)

      END DO

      IF (IPRINT .GE. INF_PROGRESS) THEN
         WRITE(LUPRI,'(80A)') ('-',I=1,LENLIN)
      END IF

      END
C  /* Deck cho_subtr */
      SUBROUTINE CHO_SUBTR(XINT,INFRED,INFVEC,
     &                     INDRED,INDRSH,IIBSTRSH,NNBSTRSH,ISCR,
     &                     MRED,MVEC,MSYM,MMSHL,MMBSTRT,LISCR,
     &                     WORK,LWORK,ISYM)
C
C     Purpose: subtract contributions from previous vectors
C              from the qualified integrals (in XINT).
C
#include <implicit.h>
      DIMENSION XINT(*), WORK(LWORK)
      INTEGER   INFRED(MRED), INFVEC(MVEC,3,MSYM)
      INTEGER   INDRED(MMBSTRT,3), INDRSH(MMBSTRT)
      INTEGER   IIBSTRSH(MSYM,MMSHL,3), NNBSTRSH(MSYM,MMSHL,3)
      INTEGER   ISCR(LISCR)
#include <cholesky.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CHO_SUBTR')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER CHO_LREAD

      PARAMETER (XMONE = -1.0D0, ONE = 1.0D0)

C     Return if nothing to do.
C     ------------------------

      IF (NUMCHO(ISYM) .LT. 1) THEN ! no prev. vectors.
         IF (LOCDBG) THEN
            WRITE(LUPRI,*) SECNAM,': nothing done because NUMCHO = ',
     &                     NUMCHO(ISYM),' (sym. ',ISYM,')'
         END IF
         RETURN
      ELSE IF (NNBSTR(ISYM,2) .LT. 1) THEN ! nothing to do (this sym.)
         IF (LOCDBG) THEN
            WRITE(LUPRI,*) SECNAM,': nothing done because NNBSTR = ',
     &                     NNBSTR(ISYM,2),' (sym. ',ISYM,')'
         END IF
         RETURN
      ELSE IF (NQUAL(ISYM) .LT. 1) THEN ! no qualifieds in this sym.
         IF (LOCDBG) THEN
            WRITE(LUPRI,*) SECNAM,': nothing done because NQUAL  = ',
     &                     NQUAL(ISYM),' (sym. ',ISYM,')'
         END IF
         RETURN
      END IF

C     Debug: read original diagonal and check that these elements are
C            included in the integrals.
C     ---------------------------------------------------------------

      IF (LOCDBG) THEN
         KDIAG = 1
         KEND  = KDIAG + NNBSTRT(1)
         LWRK  = LWORK - KEND + 1
         IF (LWRK .LT. 0) THEN
            WRITE(LUPRI,*) SECNAM,': diagonal/integral check skipped ',
     &                     'due to insufficient memory'
         ELSE
            TOL  = 1.0D-14
            NERR = 0
            CALL CHO_CHKINTO(XINT,WORK(KDIAG),INDRED,INDRSH,
     &                       IIBSTRSH,NNBSTRSH,
     &                       MSYM,MMSHL,MMBSTRT,ISYM,NERR,TOL,.TRUE.)
            IF (NERR .NE. 0) THEN
               WRITE(LUPRI,*) SECNAM,': ',NERR,' diagonal errors found!'
               WRITE(LUPRI,*) '          #tests: ',NQUAL(ISYM)
               WRITE(LUPRI,*) '          Printing integrals:'
               CALL CHO_OUTPUT(XINT,1,NNBSTR(ISYM,2),1,NQUAL(ISYM),
     &                         NNBSTR(ISYM,2),NQUAL(ISYM),1,LUPRI)
               CALL CHO_QUIT('Diagonal errors in '//SECNAM,104)
            ELSE
               WRITE(LUPRI,*) SECNAM,': comparison of qual. integrals ',
     &                     'and original diagonal: no errors !'
            END IF
         END IF
      END IF

C     Reserve space needed for reading previous vectors.
C     --------------------------------------------------

      LREAD = CHO_LREAD(ISYM,INFVEC,MVEC,MSYM,LWORK)
      IF (LREAD .LT. 1) THEN
         WRITE(LUPRI,*) SECNAM,': CHO_LREAD returned ',LREAD
         CALL CHO_QUIT('Memory error in '//SECNAM,101)
         LWRK1 = 0 ! to avoid compiler warnings
      ELSE
         LWRK1 = LWORK - LREAD
      END IF

C     Set up batch.
C     -------------

      MMEM = NNBSTR(ISYM,2) + NQUAL(ISYM)
      NVEC = MIN(LWRK1/MMEM,NUMCHO(ISYM))
      IF (NVEC .LT. 1) THEN
         CALL CHO_QUIT('Batch failure in '//SECNAM,101)
      END IF
      NBATCH = (NUMCHO(ISYM) - 1)/NVEC + 1

C     Start batch loop.
C     -----------------

      DO IBATCH = 1,NBATCH

         IF (IBATCH .EQ. NBATCH) THEN
            NUMV = NUMCHO(ISYM) - NVEC*(NBATCH - 1)
         ELSE
            NUMV = NVEC
         END IF
         IVEC1 = NVEC*(IBATCH - 1) + 1

C        Complete allocation.
C        --------------------

         KCHO1 = 1
         KCHO2 = KCHO1 + NNBSTR(ISYM,2)*NUMV
         KEND2 = KCHO2 + NQUAL(ISYM)*NUMV
         LWRK2 = LWORK - KEND2 + 1
         IF (LWRK2 .LT. LREAD) THEN
            CALL CHO_QUIT('Batch error in '//SECNAM,104)
         END IF

C        Read previous vectors.
C        ----------------------

         CALL CHO_TIMER(C1,W1)
         CALL CHO_GETVEC(WORK(KCHO1),NNBSTR(ISYM,2),NUMV,IVEC1,ISYM,
     &                   INFRED,INFVEC,INDRED,INDRSH,IIBSTRSH,NNBSTRSH,
     &                   ISCR,
     &                   MRED,MVEC,MSYM,MMBSTRT,MMSHL,MMBSTRT,
     &                   WORK(KEND2),LWRK2)
         CALL CHO_TIMER(C2,W2)
         TDECOM(1,2) = TDECOM(1,2) + C2 - C1
         TDECOM(2,2) = TDECOM(2,2) + W2 - W1

C        Copy out sub-blocks corresponding to qualified diagonals:
C        L({ab},#J)
C        ---------------------------------------------------------

         CALL CHO_TIMER(C1,W1)

         DO J = 1,NUMV
            DO IAB = 1,NQUAL(ISYM)
               KOFF1 = KCHO2 + NQUAL(ISYM)*(J - 1)
     &               + IAB - 1
               KOFF2 = KCHO1 + NNBSTR(ISYM,2)*(J - 1)
     &               + IQUAB(IAB,ISYM) - IIBSTR(ISYM,2) - 1
               WORK(KOFF1) = WORK(KOFF2)
            END DO
         END DO

C        Subtract:
C        (gd|{ab}) <- (gd|{ab}) - sum_J L(gd,#J) * L({ab},#J)
C        ----------------------------------------------------

         CALL DGEMM('N','T',NNBSTR(ISYM,2),NQUAL(ISYM),NUMV,
     &              XMONE,WORK(KCHO1),NNBSTR(ISYM,2),
     &              WORK(KCHO2),NQUAL(ISYM),ONE,XINT,NNBSTR(ISYM,2))

         CALL CHO_TIMER(C2,W2)
         TDECOM(1,3) = TDECOM(1,3) + C2 - C1
         TDECOM(2,3) = TDECOM(2,3) + W2 - W1

      END DO

      END
C  /* Deck cho_chkdia2 */
      SUBROUTINE CHO_CHKDIA2(DIAG,INDRED,INDRSH,MMBSTRT,ISYM,
     &                       XM,YM,ZM,NNEGT,NNEG,NCONV)
C
C     Purpose: 1) find min. in updated diagonal, XM (sym. ISYM only)
C              2) find max. in updated diagonal, YM (sym. ISYM only)
C              3) find abs. max. in updated diagonal, ZM (sym. ISYM only)
C              4) count #diagonals < 0.0D0, NNEGT
C              5) count #diagonals < THRNEG, NNEGT
C              6) count #screenable diagonals, NCONV
C
C     -- also: a) zero negative diagonals < THRNEG (from cholesky.h)
C              b) screen diagonal if requested (flag SCDIAG from cholesky.h)
C              c) Keep track of most negative zeroed diagonal.
C
#include <implicit.h>
      DIMENSION DIAG(*)
      INTEGER   INDRED(MMBSTRT,3), INDRSH(MMBSTRT)
#include <cholesky.h>

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CHO_CHKDIA2')

      PARAMETER (ZERO = 0.0D0)

C     Initialization.
C     ---------------

      NNEGT = 0
      NNEG  = 0
      NCONV = 0

      IF (NNBSTR(ISYM,2) .GT. 0) THEN
         JAB1 = IIBSTR(ISYM,2) + 1
         JAB2 = JAB1 + NNBSTR(ISYM,2) - 1
         XM   = DIAG(INDRED(JAB1,2))
         YM   = DIAG(INDRED(JAB1,2))
         ZM   = ABS(YM)
      ELSE
         XM = ZERO
         YM = ZERO
         ZM = ZERO
         RETURN
      END IF

C     Find min. and max. diagonal and zero too negative diagonals.
C     ------------------------------------------------------------

      DO JAB = JAB1,JAB2
         IAB = INDRED(JAB,2)  ! get address in first red. set
         XM  = MIN(XM,DIAG(IAB))
         YM  = MAX(YM,DIAG(IAB))
         IF (DIAG(IAB) .LT. ZERO) THEN
            NNEGT = NNEGT + 1
            IF (DIAG(IAB) .LT. THRNEG) THEN
               NNEG   = NNEG + 1
               IF (DIAG(IAB) .LT. TOONEG) THEN
                  WRITE(LUPRI,'(A,A,I12,1X,1P,D16.8)')
     &            SECNAM,': diagonal too negative: ',IAB,DIAG(IAB)
                  WRITE(LUPRI,'(A,A)')
     &            SECNAM,': shutting down Cholesky decomposition!'
                  CALL CHO_QUIT('Diagonal too negative in '//SECNAM,104)
               END IF
               IF (DIAG(IAB) .LT. WARNEG) THEN
                  WRITE(LUPRI,'(A,A,I12,1X,1P,D16.8,A)')
     &            SECNAM,': Negative diagonal: ',IAB,DIAG(IAB),
     &            ' (zeroed)'
               END IF
               IF (DIAG(IAB) .LT. DIAMNZ) THEN
                  DIAMNZ = DIAG(IAB)
                  IABMNZ = IAB
               END IF
               DIAG(IAB) = ZERO
            END IF
         END IF
      END DO
      ZM = MAX(ABS(XM),ABS(YM))

C     Screen diagonal (if requested) and count the screenables as
C     converged.
C     NOTE: the screening is actually identical to setting up
C           reduced sets. However, doing the screening here will
C           forces entries in later vectors of this integral pass
C           to have zero entries at screened diagonals.
C     -----------------------------------------------------------

      IF (SCDIAG) THEN
         DO JAB = JAB1,JAB2
            IAB = INDRED(JAB,2)  ! get address in first red. set
            TST = SQRT(ABS(DIAG(IAB))*ZM)*DAMP(2)
            IF (TST .LE. THRCOM) THEN
               NCONV     = NCONV + 1
               DIAG(IAB) = ZERO
            END IF
         END DO
      ELSE
         DO JAB = JAB1,JAB2
            IAB = INDRED(JAB,2)  ! get address in first red. set
            TST = SQRT(ABS(DIAG(IAB))*ZM)*DAMP(2)
            IF (TST .LE. THRCOM) THEN
               NCONV = NCONV + 1
            END IF
         END DO
      END IF

      END
C  /* Deck cho_restart2 */
      SUBROUTINE CHO_RESTART2(IIBSTRSH,NNBSTRSH,INDRED,INDRSH,
     &                        INFRED,INFVEC,ISCR,
     &                        DIAG,WORK,LWORK,
     &                        MSYM,MMSHL,MMBSTRT,MRED,MVEC,LISCR,
     &                        DSKDIA,LCONV)
C
C     Purpose: update and analyze diagonal for restart
C              (or check of decomposition). Index arrays
C              for first reduced set must be set up before
C              this routine is called. Reduced set 2, on the
C              other hand, is set up here.
C
#include <implicit.h>
      INTEGER   IIBSTRSH(MSYM,MMSHL,3), NNBSTRSH(MSYM,MMSHL,3)
      INTEGER   INDRED(MMBSTRT,3), INDRSH(MMBSTRT)
      INTEGER   INFRED(MRED), INFVEC(MVEC,3,MSYM), ISCR(LISCR)
      DIMENSION DIAG(*), WORK(LWORK)
      LOGICAL   DSKDIA, LCONV
#include <cholesky.h>

      CHARACTER*12 SECNAM
      PARAMETER (SECNAM = 'CHO_RESTART2')

      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0)

C     Read diagonal (in reduced set 1).
C     ---------------------------------

      IF (DSKDIA) THEN
         IOPT = 2
         CALL CHO_IODIAG(DIAG,IOPT)
      END IF

      WRITE(LUPRI,'(/,A,I10,/)')
     & 'Number of diagonal elements (1st reduced set): ',NNBSTRT(1)

C     Analyze diagonal before update.
C     -------------------------------

      BIN1 = 1.0D2
      STEP = 1.0D-1
      NBIN = 18
      CALL CHO_ANADIA(DIAG,BIN1,STEP,NBIN,.TRUE.)

C     Copy reduced set 1 to 2.
C     ------------------------

      CALL ICOPY(MSYM*MMSHL,IIBSTRSH(1,1,1),1,IIBSTRSH(1,1,2),1)
      CALL ICOPY(MSYM*MMSHL,NNBSTRSH(1,1,1),1,NNBSTRSH(1,1,2),1)
      CALL ICOPY(MSYM,IIBSTR(1,1),1,IIBSTR(1,2),1)
      CALL ICOPY(MSYM,NNBSTR(1,1),1,NNBSTR(1,2),1)
      NNBSTRT(2) = NNBSTRT(1)
      DO IAB = 1,NNBSTRT(1)
         INDRED(IAB,2) = IAB
      END DO

      IMXAB  = 0
      NCONVT = 0
      DO ISYM = 1,NSYM

         NDIM = NNBSTR(ISYM,2)
         NVEC = NUMCHO(ISYM)

         WRITE(LUPRI,'(//,A,I2)')
     &   'Check information, symmetry',ISYM
         WRITE(LUPRI,'(/,A,6X,I12)')
     &   'Dimension, 1st reduced set: ',NDIM
         WRITE(LUPRI,'(A,6X,I12)')
     &   'Number of Cholesky vectors: ',NVEC

         IF ((NVEC.GT.0) .AND. (NDIM.GT.0)) THEN

            KDIAG = 1
            KEND1 = KDIAG + NDIM
            LWRK1 = LWORK - KEND1 + 1

            IF (LWRK1 .LE. 0) THEN
               CALL CHO_QUIT('Insufficient memory in '//SECNAM,101)
            END IF

C           Calculate Cholesky diagonal.
C           ----------------------------

            CALL CHO_DIACHO(WORK(KDIAG),ISYM,INDRED,INDRSH,
     &                      IIBSTRSH,NNBSTRSH,INFRED,INFVEC,ISCR,
     &                      MSYM,MMSHL,MMBSTRT,MRED,MVEC,LISCR,
     &                      WORK(KEND1),LWRK1)

C           Update diagonal, find min. and max. error and save original
C           value.
C           -----------------------------------------------------------

            ERRMX = -1.0D10
            ERRMN =  1.0D10
            EXAMX = ZERO
            EXAMN = ZERO
            DO JAB = 1,NDIM
               IAB = INDRED(IIBSTR(ISYM,2)+JAB,2)
               SAV = DIAG(IAB)
               DIAG(IAB) = DIAG(IAB) - WORK(KDIAG+JAB-1)
               ERR = ABS(DIAG(IAB))
               IF (ERR .GT. ERRMX) THEN
                  ERRMX = ERR
                  EXAMX = SAV
                  IMXAB = IAB
               END IF
               IF (ERR .LT. ERRMN) THEN
                  ERRMN = ERR
                  EXAMN = SAV
               END IF
            END DO

C           Zero diagonals corresponding to Cholesky vectors,
C           i.e. those that have been explicitly treated.
C           -------------------------------------------------

            DO IVEC = 1,NVEC
               DIAG(INFVEC(IVEC,1,ISYM)) = ZERO
            END DO

C           Find min. and max. diagonals, zero too negative diagonals,
C           and screen (if requested).
C           ----------------------------------------------------------

            CALL CHO_CHKDIA2(DIAG,INDRED,INDRSH,MMBSTRT,ISYM,
     &                       XMIN,XMAX,XAMAX,NNEGT,NNEG,NSCR)

C           Count converged diagonals.
C           --------------------------

            NCONV = 0
            DO JAB = 1,NDIM
               IAB = INDRED(IIBSTR(ISYM,2)+JAB,2)
               IF (ABS(DIAG(IAB)) .LT. THRCOM) NCONV = NCONV + 1
            END DO

C           Calculate average and RMS error.
C           --------------------------------

            KOFF   = IIBSTR(ISYM,1) + 1
            XDIM   = DBLE(NDIM)
            RMSERR = SQRT(DDOT(NDIM,DIAG(KOFF),1,DIAG(KOFF),1)/XDIM)
            AVEERR = CHO_DSUMELM(DIAG(KOFF),NDIM)/XDIM

C           Print.
C           ------

            WRITE(LUPRI,'(A,1P,D18.8)')
     &      'Minimum diagonal          : ',XMIN
            WRITE(LUPRI,'(A,1P,D18.8)')
     &      'Maximum diagonal          : ',XMAX
            WRITE(LUPRI,'(A,1P,D18.8,1X,D18.8)')
     &      'Minimum absolute error    : ',ERRMN,EXAMN
            WRITE(LUPRI,'(A,1P,D18.8,1X,D18.8)')
     &      'Maximum absolute error    : ',ERRMX,EXAMX
            WRITE(LUPRI,'(A,1P,D18.8)')
     &      'Average error             : ',AVEERR
            WRITE(LUPRI,'(A,1P,D18.8)')
     &      'Root-mean-square error    : ',RMSERR
            WRITE(LUPRI,'(A,6X,I12)')
     &      'Converged diagonals       : ',NCONV
            WRITE(LUPRI,'(A,6X,I12)')
     &      'Unconverged diagonals     : ',NDIM-NCONV
            WRITE(LUPRI,'(A,6X,I12)')
     &      'Zeroed negative diagonals : ',NNEG
            IF (SCDIAG) THEN
               WRITE(LUPRI,'(A,6X,I12)')
     &         'Screened diagonals        : ',NSCR
            ELSE
               WRITE(LUPRI,'(A,6X,I12,A)')
     &         'Screenable diagonals      : ',NSCR,' (not screened)'
            END IF
            NCONVT = NCONVT + NCONV

         END IF

      END DO

C     Set reduced set.
C     ----------------

      CALL CHO_SETRED(DIAG,IIBSTRSH,NNBSTRSH,INDRED,INDRSH,
     &                MSYM,MMBSTRT,MMSHL)

C     Analyze updated diagonal.
C     -------------------------

      BIN1 = 1.0D2
      STEP = 1.0D-1
      NBIN = 18
      CALL CHO_ANADIA(DIAG,BIN1,STEP,NBIN,.FALSE.)

C     Set data for minimal integral checking.
C     ---------------------------------------

      IF (CHO_MINCHK) THEN
         IF (IMXAB .GT. 0) THEN
            ISHLAB = INDRSH(IMXAB)
            CALL CHO_INTCHK_REG('MAX ERR ',ISHLAB,ISHLAB)
         END IF
         NNSHL = NSHELL*(NSHELL + 1)/2
         JMXAB = 0
         I0AB  = 0
         XX    = 0.0D0
         DO ISHLAB = 1,NNSHL
            NTOT = 0
            DO ISYM = 1,NSYM
               IF (NNBSTRSH(ISYM,ISHLAB,1) .GT. 0) THEN
                  IAB1 = IIBSTR(ISYM,1) + IIBSTRSH(ISYM,ISHLAB,1) + 1
                  IAB2 = IAB1 + NNBSTRSH(ISYM,ISHLAB,1) - 1
                  DO IAB = IAB1,IAB2
                     IF (DIAG(IAB) .LT. XX) THEN
                        XX    = DIAG(IAB)
                        JMXAB = IAB
                     END IF
                  END DO
                  NTOT = NTOT + 1
               END IF
            END DO
            IF (NTOT .EQ. 0) THEN
               I0AB = ISHLAB
            END IF
         END DO
         IF (JMXAB .GT. 0) THEN
            JSHLAB = INDRSH(JMXAB)
            CALL CHO_INTCHK_REG('NEG DIAG',JSHLAB,JSHLAB)
            IF (IMXAB .GT. 0) THEN
               ISHLAB = INDRSH(IMXAB)
               CALL CHO_INTCHK_REG('MAX|NEG ',ISHLAB,JSHLAB)
            END IF
         END IF
         IF (I0AB .GT. 0) THEN
            JSHLAB = INDRSH(I0AB)
            CALL CHO_INTCHK_REG('EXCL RS1',JSHLAB,JSHLAB)
            IF (IMXAB .GT. 0) THEN
               ISHLAB = INDRSH(IMXAB)
               CALL CHO_INTCHK_REG('MAX|XRS1',ISHLAB,JSHLAB)
            END IF
         END IF
         IF (IABMNZ.GT.0 .AND. IABMNZ.LE.NNSHL) THEN
            JSHLAB = INDRSH(IABMNZ)
            CALL CHO_INTCHK_REG('NEG->ZER',JSHLAB,JSHLAB)
            IF (IMXAB .GT. 0) THEN
               ISHLAB = INDRSH(IMXAB)
               CALL CHO_INTCHK_REG('MAX|NEGZ',ISHLAB,JSHLAB)
            END IF
         END IF
      END IF

C     Set convergence flag.
C     ---------------------

      LCONV = NCONVT.EQ.NNBSTRT(1)

      END
C  /* Deck cho_diacho */
      SUBROUTINE CHO_DIACHO(DIAG,ISYM,INDRED,INDRSH,IIBSTRSH,NNBSTRSH,
     &                      INFRED,INFVEC,ISCR,
     &                      MSYM,MMSHL,MMBSTRT,MRED,MVEC,LISCR,
     &                      WORK,LWORK)
C
C     Purpose: calculate symmetry block ISYM of diagonal in
C              red. set 1. It is assumed that red. set 2 and
C              red. set 1 are identical!!
C
C     TODO/FIXME: operation count may be reduced by calc. directly in
C                 each reduced set. Requires batching within each
C                 reduced set individually....
C
#include <implicit.h>
      DIMENSION DIAG(*), WORK(LWORK)
      INTEGER   IIBSTRSH(MSYM,MMSHL,3), NNBSTRSH(MSYM,MMSHL,3)
      INTEGER   INDRED(MMBSTRT,3), INDRSH(MMBSTRT)
      INTEGER   INFRED(MRED), INFVEC(MVEC,3,MSYM)
      INTEGER   ISCR(LISCR)
#include <cholesky.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_DIACHO')

      INTEGER CHO_LREAD

C     Return if nothing to do.
C     ------------------------

      IF (NNBSTR(ISYM,1) .LT. 1) RETURN
      IF (NUMCHO(ISYM)   .LT. 1) RETURN

C     Check that red. set 1 and 2 are identical (superficial check).
C     --------------------------------------------------------------

      IF (NNBSTRT(2) .NE. NNBSTRT(1)) THEN
         CALL CHO_QUIT('Red. set assumption violated in '//SECNAM,104)
      END IF

C     Initialize result array.
C     ------------------------

      CALL CHO_DZERO(DIAG,NNBSTR(ISYM,1))

C     Reserve work space for read.
C     ----------------------------

      LREAD = CHO_LREAD(ISYM,INFVEC,MVEC,MSYM,LWORK)

      KREAD = 1
      KEND1 = KREAD + LREAD
      LWRK1 = LWORK - KEND1 + 1

      IF (LWRK1 .LE. 0) THEN
         CALL CHO_QUIT('Insufficient memory in '//SECNAM,101)
      END IF

C     Set up batch.
C     -------------

      MMEM = NNBSTR(ISYM,1)
      NVEC = MIN(LWRK1/MMEM,NUMCHO(ISYM))
      IF (NVEC .LT. 1) THEN
         CALL CHO_QUIT('Batch failure in '//SECNAM,101)
      END IF
      NBATCH = (NUMCHO(ISYM) - 1)/NVEC + 1

C     Start batch loop.
C     -----------------

      DO IBATCH = 1,NBATCH

         IF (IBATCH .EQ. NBATCH) THEN
            NUMV = NUMCHO(ISYM) - NVEC*(NBATCH - 1)
         ELSE
            NUMV = NVEC
         END IF
         IVEC1 = NVEC*(IBATCH - 1) + 1
         KCHO1 = KEND1

C        Read Cholesky vectors.
C        ----------------------

         CALL CHO_GETVEC(WORK(KCHO1),NNBSTR(ISYM,2),NUMV,IVEC1,ISYM,
     &                   INFRED,INFVEC,INDRED,INDRSH,
     &                   IIBSTRSH,NNBSTRSH,ISCR,
     &                   MRED,MVEC,MSYM,MMBSTRT,MMSHL,LISCR,
     &                   WORK(KREAD),LREAD)

C        Calculate contributions to diagonal.
C        ------------------------------------

         DO IVEC = 1,NUMV
            DO IAB = 1,NNBSTR(ISYM,1)
               KOFF = KCHO1 + NNBSTR(ISYM,1)*(IVEC - 1) + IAB - 1
               DIAG(IAB) = DIAG(IAB) + WORK(KOFF)*WORK(KOFF)
            END DO
         END DO

      END DO

      END
C  /* Deck cho_putred */
      SUBROUTINE CHO_PUTRED(INFRED,NNBSTRSH,INDRED,INDRSH,
     &                      MRED,MSYM,MMSHL,MMBSTRT,IPASS)
C
C     Purpose: write index arrays for current reduced set (reduced set
C              IPASS).
C
#include <implicit.h>
      INTEGER INFRED(MRED)
      INTEGER NNBSTRSH(MSYM*MMSHL), INDRED(MMBSTRT), INDRSH(MMBSTRT)
#include <cholesky.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_PUTRED')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

C     Test dimensions.
C     ----------------

      NNSHL = NSHELL*(NSHELL + 1)/2

      IF (MSYM .NE. NSYM) THEN
         CALL CHO_QUIT('NSYM error in '//SECNAM,104)
      END IF

      IF (MMSHL .NE. NNSHL) THEN
         CALL CHO_QUIT('NNSHL error in '//SECNAM,104)
      END IF

      IF (MMBSTRT .NE. NNBSTRT(1)) THEN
         CALL CHO_QUIT('NNBSTRT(1) error in '//SECNAM,104)
      END IF

      IF (MMBSTRT .LT. NNBSTRT(2)) THEN
         CALL CHO_QUIT('NNBSTRT(2) error in '//SECNAM,104)
      END IF

      IF ((IPASS.LT.1) .OR. (IPASS.GT.MAXRED)) THEN
         CALL CHO_QUIT('IPASS error in '//SECNAM,104)
      END IF

C     Get first address.
C     ------------------

      IADR1 = INFRED(IPASS)
      IF (IADR1 .LT. 1) THEN
         WRITE(LUPRI,*) SECNAM,': nonpositive address for reduced set ',
     &                  IPASS,': ',IADR1
         CALL CHO_QUIT('Error in '//SECNAM,104)
      END IF

      IF (LOCDBG) THEN
         WRITE(LUPRI,*) SECNAM,': putting reduced set ',IPASS,
     &                  ' at addr: ',IADR1,' on unit ',LURED
         IF (IPASS .EQ. 1) THEN
            WRITE(LUPRI,*) SECNAM,': putting INDRSH on unit ',LURSH
         END IF
         CALL CHO_FLUSH(LUPRI)
      END IF

C     Write index arrays.
C     -------------------

      IADR = IADR1
      WRITE(LURED,REC=IADR) (NNBSTRSH(I),I=1,NSYM*NNSHL),
     &                      (INDRED(J),J=1,NNBSTRT(1))
      IF (IPASS .EQ. 1) THEN
         REWIND(LURSH)
         WRITE(LURSH) (INDRSH(I),I=1,NNBSTRT(1))
      END IF

      END
C  /* Deck cho_getred */
      SUBROUTINE CHO_GETRED(INFRED,NNBSTRSH,INDRED,INDRSH,
     &                      MRED,MSYM,MMSHL,MMBSTRT,
     &                      IPASS,LRSH)
C
C     Purpose: read index arrays for current reduced set (reduced set
C              IPASS).
C
#include <implicit.h>
      INTEGER INFRED(MRED)
      INTEGER NNBSTRSH(MSYM*MMSHL), INDRED(MMBSTRT), INDRSH(MMBSTRT)
      LOGICAL LRSH
#include <cholesky.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_GETRED')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER CHO_ISUMELM

C     Test dimensions.
C     ----------------

      NNSHL = NSHELL*(NSHELL + 1)/2

      IF (MSYM .NE. NSYM) THEN
         CALL CHO_QUIT('NSYM error in '//SECNAM,104)
      END IF

      IF (MMSHL .NE. NNSHL) THEN
         CALL CHO_QUIT('NNSHL error in '//SECNAM,104)
      END IF

      IF (MMBSTRT .NE. NNBSTRT(1)) THEN
         CALL CHO_QUIT('NNBSTRT(1) error in '//SECNAM,104)
      END IF

      IF ((IPASS.LT.1) .OR. (IPASS.GT.MAXRED)) THEN
         CALL CHO_QUIT('IPASS error in '//SECNAM,104)
      END IF

C     Get first address.
C     ------------------

      IADR1 = INFRED(IPASS)
      IF (IADR1 .LT. 1) THEN
         WRITE(LUPRI,*) SECNAM,': nonpositive address for reduced set ',
     &                  IPASS,': ',IADR1
         CALL CHO_QUIT('Error in '//SECNAM,104)
      END IF

      IF (LOCDBG) THEN
         WRITE(LUPRI,*) SECNAM,': getting reduced set ',IPASS,
     &                  ' at addr: ',IADR1
      END IF

C     Read index arrays.
C     ------------------

      IADR = IADR1
      READ(LURED,REC=IADR) (NNBSTRSH(I),I=1,NSYM*NNSHL),
     &                     (INDRED(I),I=1,NNBSTRT(1))
      IF (LRSH) THEN
         REWIND(LURSH)
         READ(LURSH) (INDRSH(I),I=1,NNBSTRT(1))
      END IF

      END
C  /* Deck cho_putvec */
      SUBROUTINE CHO_PUTVEC(CHOVEC,LENVEC,NUMVEC,IVEC1,ISYM,
     &                      INFVEC,MVEC,MSYM)
C
C     Purpose: write Cholesky vectors IVEC=IVEC1,...,IVEC1+NUMVEC-1
C              of symmetry ISYM to file.
C
#include <implicit.h>
      DIMENSION CHOVEC(LENVEC,NUMVEC)
      INTEGER   INFVEC(MVEC,3,MSYM)
#include <cholesky.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_PUTVEC')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER*8 IADR

C     Return if no vectors.
C     ---------------------

      IF (NUMVEC .LT. 1) THEN
         IF (LOCDBG) THEN
            WRITE(LUPRI,*) SECNAM,': WARNING: no vectors in this call!'
            WRITE(LUPRI,*) SECNAM,': NUMVEC = ',NUMVEC
         END IF
         RETURN
      END IF

C     Check vector dimension: should be the same as current reduced
C     set.
C     -------------------------------------------------------------

      IF (LENVEC .NE. NNBSTR(ISYM,2)) THEN
         CALL CHO_QUIT('Illegal vector dimension in '//SECNAM,104)
      END IF
      IF (LENVEC .LT. 1) THEN
         IF (LOCDBG) THEN
            WRITE(LUPRI,*) SECNAM,': WARNING: negative vector dimension'
            WRITE(LUPRI,*) SECNAM,': LENVEC = ',LENVEC
         END IF
         RETURN
      END IF

C     Check symmetry.
C     ---------------

      IF ((ISYM.LT.1) .OR. (ISYM.GT.NSYM)) THEN
         WRITE(LUPRI,*) SECNAM,': symmetry out of bounds'
         WRITE(LUPRI,*) 'ISYM = ',ISYM
         CALL CHO_QUIT('Symmetry out of bounds in '//SECNAM,104)
      END IF

C     Check vector index.
C     -------------------

      IVEC2 = IVEC1 + NUMVEC - 1
      IADR2 = INFVEC(IVEC2,3,ISYM)
      IF ((IVEC1.LT.1) .OR. (IVEC1.GT.MAXVEC) .OR.
     &    (IVEC2.LT.1) .OR. (IVEC2.GT.MAXVEC)) THEN
         WRITE(LUPRI,*) SECNAM,': vector index out of bounds'
         WRITE(LUPRI,*) 'IVEC1 = ',IVEC1,' IVEC2 = ',IVEC2
         WRITE(LUPRI,*) '...must be between 1 and ',MAXVEC
         CALL CHO_QUIT('Vector index out of bounds in '
     &                 //SECNAM,104)
      ELSE IF (INFVEC(IVEC1,3,ISYM) .LT. 1) THEN
         WRITE(LUPRI,*) 'Error in ',SECNAM,':'
         WRITE(LUPRI,*) 'Illegal disk address for first vector: ',
     &                  INFVEC(IVEC1,3,ISYM)
         IF (INFVEC(IVEC1,3,ISYM) .LT. 0) THEN
            WRITE(LUPRI,*) '....is it an overflow?'
         END IF
         WRITE(LUPRI,*) 'IVEC1 = ',IVEC1,' ISYM = ',ISYM
         CALL CHO_QUIT('Illegal disk address in '//SECNAM,104)
      ELSE IF (IADR2 .LT. INFVEC(IVEC1,3,ISYM)) THEN
         WRITE(LUPRI,*) 'Error in ',SECNAM,':'
         WRITE(LUPRI,*) 'Illegal disk address for last vector: ',
     &                  IADR2
         IF (IADR2 .LT. 0) THEN
            WRITE(LUPRI,*) '....is it an overflow?'
         END IF
         WRITE(LUPRI,*) 'IVEC2 = ',IVEC2,' ISYM = ',ISYM
         CALL CHO_QUIT('Illegal disk address in '//SECNAM,104)
      END IF

C     Call the low-level I/O routines.
C     --------------------------------

      IADR = INFVEC(IVEC1,3,ISYM)
      LTOT = LENVEC*NUMVEC

      CALL PUTWA2U(LUNCHO(ISYM),FNVEC(ISYM),CHOVEC,IADR,LTOT)

C     Debug stuff.
C     ------------

      IF (LOCDBG) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) SECNAM,':'
         WRITE(LUPRI,*) 'Vectors ',IVEC1,' to ',IVEC1+NUMVEC-1,
     &                  ' of symmetry ',ISYM,' written to unit ',
     &                  LUNCHO(ISYM)
         WRITE(LUPRI,*) 'Vector dimension: ',LENVEC
         DO IVEC = 1,NUMVEC
            JVEC = IVEC1 + IVEC - 1
            JADR = INFVEC(JVEC,3,ISYM)
            XNRM = SQRT(DDOT(LENVEC,CHOVEC(1,IVEC),1,CHOVEC(1,IVEC),1))
            WRITE(LUPRI,*) 'Vector:',JVEC,' address: ',JADR,' norm: ',
     &                     XNRM
         END DO
      END IF

      END
C  /* Deck cho_chkinto */
      SUBROUTINE CHO_CHKINTO(XINT,DIAG,INDRED,INDRSH,IIBSTRSH,NNBSTRSH,
     &                       MSYM,MMSHL,MMBSTRT,ISYM,NERR,TOL,REPORT)
C
C     Purpose: check diagonals in qualified integral columns against
C              original diagonal (read in here).
C
#include <implicit.h>
      DIMENSION XINT(*), DIAG(*)
      INTEGER   INDRED(MMBSTRT,3), INDRSH(MMBSTRT)
      INTEGER   IIBSTRSH(MSYM,MMSHL,3), NNBSTRSH(MSYM,MMSHL,3)
      LOGICAL   REPORT
#include <cholesky.h>

c     REWIND(LUDIAG)
c     READ(LUDIAG) (DIAG(I), I=1,NNBSTRT(1))
      IOPT = 2
      CALL CHO_IODIAG(DIAG,IOPT)
      CALL CHO_CHKINT(XINT,DIAG,INDRED,INDRSH,IIBSTRSH,NNBSTRSH,
     &                MSYM,MMSHL,MMBSTRT,ISYM,NERR,TOL,REPORT)

      END
C  /* Deck cho_chkint */
      SUBROUTINE CHO_CHKINT(XINT,DIAG,INDRED,INDRSH,IIBSTRSH,NNBSTRSH,
     &                      MSYM,MMSHL,MMBSTRT,ISYM,NERR,TOL,REPORT)
C
C     Purpose: check diagonals in qualified integral columns.
C
#include <implicit.h>
      DIMENSION XINT(*), DIAG(*)
      INTEGER   INDRED(MMBSTRT,3), INDRSH(MMBSTRT)
      INTEGER   IIBSTRSH(MSYM,MMSHL,3), NNBSTRSH(MSYM,MMSHL,3)
      LOGICAL   REPORT
#include <cholesky.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_CHKINT')

      NERR = 0
      DO I = 1,NQUAL(ISYM)
         II = IQUAB(I,ISYM)
         JJ = INDRED(II,2)
         IK = II - IIBSTR(ISYM,2)
         KK = NNBSTR(ISYM,2)*(I - 1) + IK
         DF = DIAG(JJ) - XINT(KK)
         IF (ABS(DF) .GT. TOL) THEN
            NERR = NERR + 1
            IF (REPORT) THEN
               WRITE(LUPRI,*) SECNAM,': diag error: ',DIAG(JJ),XINT(KK)
               WRITE(LUPRI,*) '            diagonal elm    : ',JJ,
     &                        ' (rs1) ',II,' (rs2)'
               WRITE(LUPRI,*) '            integral row,col: ',IK,I
            END IF
         END IF
      END DO

      END
C  /* Deck cho_stat */
      SUBROUTINE CHO_STAT(INFRED,INFVEC,IIBSTRSH,NNBSTRSH,
     &                    INDRED,INDRSH,INTMAP,
     &                    MRED,MVEC,MSYM,MMSHL,MMBSTRT)
C
C     Purpose: print statistics from decomposition.
C
#include <implicit.h>
      INTEGER INFRED(MRED), INFVEC(MVEC,3,MSYM)
      INTEGER IIBSTRSH(MSYM,MMSHL,3), NNBSTRSH(MSYM,MMSHL,3)
      INTEGER INDRED(MMBSTRT,3), INDRSH(MMBSTRT)
      INTEGER INTMAP(MMSHL)
#include <cholesky.h>
#include <choorb.h>

      CHARACTER*8 SECNAM
      PARAMETER (SECNAM = 'CHO_STAT')

      CHARACTER*25 STRING
      CHARACTER*2  UNT

      REAL*8 XXBST(8), VCSTOR(8)

      LOGICAL DOCPCT, DOWPCT

      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J

C     Dimension test.
C     ---------------

      NNSHL = NSHELL*(NSHELL + 1)/2
      IF (MMSHL .NE. NNSHL) CALL CHO_QUIT('Dimension error in '//SECNAM,
     &                                    104)
      IF (NNSHL .LT. 1) RETURN

C     Overall header.
C     ---------------

      CALL CHO_HEAD('Cholesky Decomposition Statistics','=',80,LUPRI)
      IF (RSTDIA) THEN
         WRITE(LUPRI,'(/,A)')
     &   'Calculation restarted from diagonal on disk'
      END IF
      IF (RSTCHO) THEN
         IF (RSTDIA) THEN
            WRITE(LUPRI,'(A)')
     &      'Calculation restarted from Cholesky vectors on disk'
         ELSE
            WRITE(LUPRI,'(/,A)')
     &      'Calculation restarted from Cholesky vectors on disk'
         END IF
      END IF

C     Configuration.
C     --------------

      CALL CHO_HEAD('Configuration','-',80,LUPRI)
      CALL CHO_PRTHEAD(.TRUE.)

C     Vector statistics.
C     ------------------

      CALL CHO_HEAD('Vector statistics','-',80,LUPRI)
      WRITE(LUPRI,'(/,A,/,A,/,A)')
     & '              #Vectors                  N                ',
     & 'Sym.     Maximum    Actual %Saving (#Orbitals) #Vectors/N',
     & '---------------------------------------------------------'
      NTOT = 0
      DO ISYM = 1,NSYM
         NN = 0
         XXBST(ISYM) = 0.0D0
         DO JSYM = 1,NSYM
            KSYM = MULD2H(JSYM,ISYM)
            IF (JSYM .GT. KSYM) THEN
               NN = NN + NBAS(JSYM)*NBAS(KSYM)
               XXBST(ISYM) = XXBST(ISYM)
     &                     + DBLE(NBAS(JSYM))*DBLE(NBAS(KSYM))
            ELSE IF (JSYM .EQ. KSYM) THEN
               NN = NN + NBAS(JSYM)*(NBAS(JSYM) + 1)/2
               XXBST(ISYM) = XXBST(ISYM)
     &            + DBLE(NBAS(JSYM))*(DBLE(NBAS(JSYM)) + 1.0D0)/2.0D0
            END IF
         END DO
         IF (XXBST(ISYM) .GT. 0.0D0) THEN
            SAV = 1.0D2*(XXBST(ISYM)-DBLE(NUMCHO(ISYM)))/XXBST(ISYM)
         ELSE
            SAV = 9.0D9
         END IF
         IF (NBAS(ISYM) .NE. 0) THEN
            X = DBLE(NBAS(ISYM))
            WRITE(LUPRI,'(I3,4X,I9,1X,I9,1X,F7.2,1X,I7,5X,F10.3)')
     &      ISYM,NN,NUMCHO(ISYM),SAV,NBAS(ISYM),DBLE(NUMCHO(ISYM))/X
         ELSE
            WRITE(LUPRI,'(I3,4X,I9,1X,I9,1X,F7.2,1X,I7)')
     &      ISYM,NN,NUMCHO(ISYM),SAV,NBAS(ISYM)
         END IF
         NTOT = NTOT + NUMCHO(ISYM)
      END DO
      NN  = NBAST*(NBAST + 1)/2
      XX  = DBLE(NBAST)*(DBLE(NBAST) + 1.0D0)/2.0D0
      IF (XX .GT. 0.0D0) THEN
         SAV = 1.0D2*(XX-DBLE(NUMCHT))/XX
      ELSE
         SAV = 9.0D9
      END IF
      WRITE(LUPRI,'(A)')
     & '---------------------------------------------------------'
      IF (NBAST .NE. 0) THEN
         XT = DBLE(NBAST)
         WRITE(LUPRI,'(A,I9,1X,I9,1X,F7.2,1X,I7,5X,F10.3)')
     &   'Total: ',NN,NUMCHT,SAV,NBAST,DBLE(NUMCHT)/XT
      ELSE
         WRITE(LUPRI,'(A,I9,1X,I9,1X,F7.2,1X,I7)')
     &   'Total: ',NN,NUMCHT,SAV,NBAST
      END IF
      WRITE(LUPRI,'(A)')
     & '---------------------------------------------------------'
      IF (NTOT .NE. NUMCHT) THEN
         WRITE(LUPRI,'(A)')
     &   'WARNING: total number of vectors is wrong!!!'
      END IF

      CALL CHO_GETSTOR(INFRED,INFVEC,IIBSTRSH,NNBSTRSH,INDRED,
     &                 INDRSH,VCSTOR,MRED,MVEC,MSYM,MMSHL,MMBSTRT)

      WRITE(LUPRI,'(/,A,/,A,/,A)')
     & '                       %Saving relative to',
     & 'Sym.     Storage      1st Red. Set     Full',
     & '---------------------------------------------'
      VCTOT = 0.0D0
      X1TOT = 0.0D0
      XXTOT = 0.0D0
      DO ISYM = 1,NSYM
         X1RED = DBLE(NNBSTR(ISYM,1))*DBLE(NUMCHO(ISYM))
         IF (X1RED .GT. 0.0D0) THEN
            SAV1  = 1.0D2*(X1RED - VCSTOR(ISYM))/X1RED
         ELSE
            SAV1 = 9.0D9
         END IF
         XX = XXBST(ISYM)*DBLE(NUMCHO(ISYM))
         IF (XX .GT. 0.0D0) THEN
            SAV2  = 1.0D2*(XX - VCSTOR(ISYM))/XX
         ELSE
            SAV1 = 9.0D9
         END IF
         XGB = VCSTOR(ISYM)*8.0D0/1.024D3
         UNT = 'kb'
         IF (XGB .GT. 1.0D3) THEN
            XGB = XGB/1.024D3
            UNT = 'Mb'
            IF (XGB .GT. 1.0D3) THEN
               XGB = XGB/1.024D3
               UNT = 'Gb'
            END IF
         END IF
         WRITE(LUPRI,'(I2,4X,F10.3,1X,A,4X,F9.4,4X,F9.4)')
     &   ISYM,XGB,UNT,SAV1,SAV2
         VCTOT = VCTOT + VCSTOR(ISYM)
         X1TOT = X1TOT + X1RED
         XXTOT = XXTOT + XX
      END DO
      WRITE(LUPRI,'(A)')
     & '---------------------------------------------'
      IF (X1TOT .GT. 0.0D0) THEN
         SAV1  = 1.0D2*(X1TOT - VCTOT)/X1TOT
      ELSE
         SAV1 = 9.0D9
      END IF
      IF (XXTOT .GT. 0.0D0) THEN
         SAV2  = 1.0D2*(XXTOT - VCTOT)/XXTOT
      ELSE
         SAV2 = 9.0D9
      END IF
      XGB = VCTOT*8.0D0/1.024D3
      UNT = 'kb'
      IF (XGB .GT. 1.0D3) THEN
         XGB = XGB/1.024D3
         UNT = 'Mb'
         IF (XGB .GT. 1.0D3) THEN
            XGB = XGB/1.024D3
            UNT = 'Gb'
            IF (XGB .GT. 1.0D3) THEN
               XGB = XGB/1.024D3
               UNT = 'Tb'
            END IF
         END IF
      END IF
      WRITE(LUPRI,'(A6,F10.3,1X,A,4X,F9.4,4X,F9.4)')
     & 'Total:',XGB,UNT,SAV1,SAV2
      WRITE(LUPRI,'(A)')
     & '---------------------------------------------'

C     Integral statistics.
C     --------------------

      MAXCAL = 0
      NCAL   = 0
      NREP   = 0
      DO ISHLAB = 1,NNSHL
         MAXCAL = MAX(MAXCAL,INTMAP(ISHLAB))
         IF (INTMAP(ISHLAB) .GT. 0) THEN
            NCAL = NCAL + 1
            IF (INTMAP(ISHLAB) .GT. 1) NREP = NREP + 1
         END IF
      END DO

      XXSHL = DBLE(NNSHL)
      XCAL  = DBLE(NCAL)
      XREP  = DBLE(NREP)

      CALL CHO_HEAD('Integral statistics','-',80,LUPRI)
      WRITE(LUPRI,'(/,A,I10)')
     & '#Shells                  :',NSHELL
      WRITE(LUPRI,'(A,I10)')
     & '#Shell Pair Distributions:',NNSHL
      WRITE(LUPRI,'(A,I10,A,F8.3,A,A)')
     & '#Shell Pairs Calculated  :',NCAL,' (',XCAL*1.0D2/XXSHL,' %',
     & ' of total)'
      WRITE(LUPRI,'(A,I10,A,F8.3,A,A)')
     & '#Shell Pairs Repeated    :',NREP,' (',XREP*1.0D2/XCAL,' %',
     & ' of calculated)'

      WRITE(LUPRI,'(/,A,/,A)')
     & '#Calculations     #Shell Pairs   Percentage',
     & '-------------------------------------------'
      DO ICAL = 1,MAXCAL
         N = 0
         DO ISHLAB = 1,NNSHL
            IF (INTMAP(ISHLAB) .EQ. ICAL) N = N + 1
         END DO
         IF (N .GT. 0) THEN
            X = DBLE(N)
            WRITE(LUPRI,'(I12,6X,I12,5X,F8.3)') ICAL,N,X*1.0D2/XXSHL
         END IF
      END DO
      WRITE(LUPRI,'(A)')
     & '-------------------------------------------'

C     Section timings.
C     ----------------

      TCINI = TIMSEC(2,1) - TIMSEC(1,1)
      TWINI = TIMSEC(4,1) - TIMSEC(3,1)
      TCDIA = TIMSEC(2,2) - TIMSEC(1,2)
      TWDIA = TIMSEC(4,2) - TIMSEC(3,2)
      TCDEC = TIMSEC(2,3) - TIMSEC(1,3)
      TWDEC = TIMSEC(4,3) - TIMSEC(3,3)
      TCCHD = TIMSEC(2,4) - TIMSEC(1,4)
      TWCHD = TIMSEC(4,4) - TIMSEC(3,4)
      TCCHA = TIMSEC(2,5) - TIMSEC(1,5)
      TWCHA = TIMSEC(4,5) - TIMSEC(3,5)

      CALL CHO_HEAD('Section timings','-',80,LUPRI)
      WRITE(LUPRI,'(/,A,/,A,/,A)')
     &'                                 CPU time          Wall time',
     &'Section                      hours min. sec.    hours min. sec.',
     &'---------------------------------------------------------------'
      STRING = 'Initialization           '
      CALL CHO_CNVTIM(TCINI,IHC,IMC,SCC)
      CALL CHO_CNVTIM(TWINI,IHW,IMW,SCW)
      WRITE(LUPRI,'(A,1X,I8,2X,I2,1X,F5.1,1X,I8,2X,I2,1X,F5.1)')
     & STRING,IHC,IMC,SCC,IHW,IMW,SCW
      STRING = 'Diagonal setup           '
      CALL CHO_CNVTIM(TCDIA,IHC,IMC,SCC)
      CALL CHO_CNVTIM(TWDIA,IHW,IMW,SCW)
      WRITE(LUPRI,'(A,1X,I8,2X,I2,1X,F5.1,1X,I8,2X,I2,1X,F5.1)')
     & STRING,IHC,IMC,SCC,IHW,IMW,SCW
      STRING = 'Cholesky decomposition   '
      CALL CHO_CNVTIM(TCDEC,IHC,IMC,SCC)
      CALL CHO_CNVTIM(TWDEC,IHW,IMW,SCW)
      WRITE(LUPRI,'(A,1X,I8,2X,I2,1X,F5.1,1X,I8,2X,I2,1X,F5.1)')
     & STRING,IHC,IMC,SCC,IHW,IMW,SCW
      STRING = 'Diagonal check           '
      CALL CHO_CNVTIM(TCCHD,IHC,IMC,SCC)
      CALL CHO_CNVTIM(TWCHD,IHW,IMW,SCW)
      WRITE(LUPRI,'(A,1X,I8,2X,I2,1X,F5.1,1X,I8,2X,I2,1X,F5.1)')
     & STRING,IHC,IMC,SCC,IHW,IMW,SCW
      IF (CHO_INTCHK) THEN
         STRING = 'Integral check (debug)   '
         CALL CHO_CNVTIM(TCCHA,IHC,IMC,SCC)
         CALL CHO_CNVTIM(TWCHA,IHW,IMW,SCW)
         WRITE(LUPRI,'(A,1X,I8,2X,I2,1X,F5.1,1X,I8,2X,I2,1X,F5.1)')
     &   STRING,IHC,IMC,SCC,IHW,IMW,SCW
      END IF
      WRITE(LUPRI,'(A)')
     &'---------------------------------------------------------------'

C     Timing of decomposition driver.
C     -------------------------------

      IF (DID_DECDRV) THEN

         CMISC = TDECDRV(1)
         WMISC = TDECDRV(2)
         DO J = 1,NINTEG
            CMISC = CMISC - TINTEG(1,J)
            WMISC = WMISC - TINTEG(2,J)
         END DO
         DO J = 1,NDECOM
            CMISC = CMISC - TDECOM(1,J)
            WMISC = WMISC - TDECOM(2,J)
         END DO
         CPCT = -9.0D9
         CFAC = -9.0D9
         WPCT = -9.0D9
         WFAC = -9.0D9
         IF (TDECDRV(1) .GT. 0.0D0) THEN
            DOCPCT = .TRUE.
            CFAC   = 1.0D2/TDECDRV(1)
         ELSE
            DOCPCT = .FALSE.
         END IF
         IF (TDECDRV(2) .GT. 0.0D0) THEN
            DOWPCT = .TRUE.
            WFAC   = 1.0D2/TDECDRV(2)
         ELSE
            DOWPCT = .FALSE.
         END IF

         CALL CHO_HEAD('Timing of decomposition driver','-',80,LUPRI)
         WRITE(LUPRI,'(/,A,A,/,A,A)')
     &    'Task           Component           CPU (min.)     %',
     &    '   Wall (min.)     %',
     &    '---------------------------------------------------',
     &    '----------------------'
         CTIM = TINTEG(1,1)
         WTIM = TINTEG(2,1)
         IF (DOCPCT) CPCT = CTIM*CFAC
         IF (DOWPCT) WPCT = WTIM*WFAC
         WRITE(LUPRI,'(A,1X,F10.2,1X,F7.2,2X,F10.2,1X,F7.2)')
     &   'Integrals      calculation        ',
     &   CTIM/6.0D1,CPCT,WTIM/6.0D1,WPCT
         CTIM = TINTEG(1,2)
         WTIM = TINTEG(2,2)
         IF (DOCPCT) CPCT = CTIM*CFAC
         IF (DOWPCT) WPCT = WTIM*WFAC
         WRITE(LUPRI,'(A,1X,F10.2,1X,F7.2,2X,F10.2,1X,F7.2)')
     &   '               I/O, qualifieds    ',
     &   CTIM/6.0D1,CPCT,WTIM/6.0D1,WPCT
         CTIM = TDECOM(1,1)
         WTIM = TDECOM(2,1)
         IF (DOCPCT) CPCT = CTIM*CFAC
         IF (DOWPCT) WPCT = WTIM*WFAC
         WRITE(LUPRI,'(A,1X,F10.2,1X,F7.2,2X,F10.2,1X,F7.2)')
     &   'Decomposition  I/O, qualifieds    ',
     &   CTIM/6.0D1,CPCT,WTIM/6.0D1,WPCT
         CTIM = TDECOM(1,2)
         WTIM = TDECOM(2,2)
         IF (DOCPCT) CPCT = CTIM*CFAC
         IF (DOWPCT) WPCT = WTIM*WFAC
         WRITE(LUPRI,'(A,1X,F10.2,1X,F7.2,2X,F10.2,1X,F7.2)')
     &   '               I/O, vectors       ',
     &   CTIM/6.0D1,CPCT,WTIM/6.0D1,WPCT
         CTIM = TDECOM(1,3)
         WTIM = TDECOM(2,3)
         IF (DOCPCT) CPCT = CTIM*CFAC
         IF (DOWPCT) WPCT = WTIM*WFAC
         WRITE(LUPRI,'(A,1X,F10.2,1X,F7.2,2X,F10.2,1X,F7.2)')
     &   '               vector subtraction ',
     &   CTIM/6.0D1,CPCT,WTIM/6.0D1,WPCT
         CTIM = CMISC
         WTIM = WMISC
         IF (DOCPCT) CPCT = CTIM*CFAC
         IF (DOWPCT) WPCT = WTIM*WFAC
         WRITE(LUPRI,'(A,1X,F10.2,1X,F7.2,2X,F10.2,1X,F7.2)')
     &   'Misc.          qualification etc. ',
     &   CTIM/6.0D1,CPCT,WTIM/6.0D1,WPCT
         WRITE(LUPRI,'(A,A)')
     &   '---------------------------------------------------',
     &   '----------------------'
         WRITE(LUPRI,'(A,1X,F10.2,10X,F10.2)')
     &   'Total:                            ',
     &   TDECDRV(1)/6.0D1,TDECDRV(2)/6.0D1
         WRITE(LUPRI,'(A,A)')
     &   '---------------------------------------------------',
     &   '----------------------'

      END IF

      END
C  /* Deck cho_getstor */
      SUBROUTINE CHO_GETSTOR(INFRED,INFVEC,IIBSTRSH,NNBSTRSH,INDRED,
     &                       INDRSH,VCSTOR,MRED,MVEC,MSYM,MMSHL,MMBSTRT)
C
C     Purpose: get total vector storage (in words).
C
#include <implicit.h>
      INTEGER   INFRED(MRED), INFVEC(MVEC,3,MSYM)
      INTEGER   IIBSTRSH(MSYM,MMSHL,3), NNBSTRSH(MSYM,MMSHL,3)
      INTEGER   INDRED(MMBSTRT,3), INDRSH(MMBSTRT)
      DIMENSION VCSTOR(MSYM)
#include <cholesky.h>

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CHO_GETSTOR')

      DO ISYM = 1,NSYM
         IF (NUMCHO(ISYM) .GT. MAXVEC) THEN
            WRITE(LUPRI,*) SECNAM,': too many Cholesky vectors ',
     &                     'in symmetry ',ISYM,': ',NUMCHO(ISYM)
            CALL CHO_QUIT('Error in '//SECNAM,103)
            VCSTOR(ISYM) = 0.0D0   ! to avoid compiler warnings
         ELSE IF (NUMCHO(ISYM) .LT. 1) THEN
            VCSTOR(ISYM) = 0.0D0
         ELSE
            IRED = INFVEC(NUMCHO(ISYM),2,ISYM)
            JRED = 3
            CALL CHO_GETRED(INFRED,NNBSTRSH(1,1,JRED),INDRED(1,JRED),
     &                      INDRSH,MRED,MSYM,MMSHL,MMBSTRT,IRED,.FALSE.)
            CALL CHO_SETREDIND(IIBSTRSH,NNBSTRSH,MSYM,MMSHL,JRED)
            VCSTOR(ISYM) = INFVEC(NUMCHO(ISYM),3,ISYM)
     &                   + NNBSTR(ISYM,JRED) - 1
         END IF
      END DO

      END
C  /* Deck cho_dbgint */
      SUBROUTINE CHO_DBGINT(INFRED,INFVEC,
     &                      IIBSTRSH,NNBSTRSH,ISCR,INDRED,INDRSH,
     &                      MRED,MVEC,MSYM,MMSHL,MMBSTRT,LISCR)
C
C     Purpose: regenerate and check integrals as specified in input.
C
#include <implicit.h>
      INTEGER   INFRED(MRED), INFVEC(MVEC,3,MSYM), ISCR(LISCR)
      INTEGER   IIBSTRSH(MSYM,MMSHL,3), NNBSTRSH(MSYM,MMSHL,3)
      INTEGER   INDRED(MMBSTRT,3), INDRSH(MMBSTRT)
#include <cholesky.h>
#include <choorb.h>

      WRITE(LUPRI,*)
      WRITE(LUPRI,*)
      WRITE(LUPRI,*) '********************************************',
     &               '********************************'
      WRITE(LUPRI,*) '* NOTICE: integral check has not (yet) been ',
     &               'implemented. Check is skipped. *'
      WRITE(LUPRI,*) '********************************************',
     &               '********************************'
      WRITE(LUPRI,*)
      WRITE(LUPRI,*)

      END
C  /* Deck cho_dbgint_cho */
      SUBROUTINE CHO_DBGINT_CHO(XINT,NCD,NAB,WORK,LWORK,
     &                          INFRED,INFVEC,
     &                          IIBSTRSH,NNBSTRSH,ISCR,INDRED,INDRSH,
     &                          MRED,MVEC,MSYM,MMSHL,MMBSTRT,LISCR,
     &                          ERRMAX,ERRMIN,ERRRMS,NCMP,
     &                          ISHLC,ISHLD,ISHLA,ISHLB)
C
C     Purpose: calculate integrals in shell quadruple (CD|AB) from
C              Cholesky vectors on disk and compare to those in
C              XINT (for debugging).
C
C     NOTE: this is *only* for debugging.
C
#include <implicit.h>
      DIMENSION XINT(NCD,NAB), WORK(LWORK)
      INTEGER   INFRED(MRED), INFVEC(MVEC,3,MSYM)
      INTEGER   IIBSTRSH(MSYM,MMSHL,3), NNBSTRSH(MSYM,MMSHL,3)
      INTEGER   ISCR(LISCR)
      INTEGER   INDRED(MMBSTRT,3), INDRSH(MMBSTRT)
#include <cholesky.h>

      CHARACTER*14 SECNAM
      PARAMETER (SECNAM = 'CHO_DBGINT_CHO')

      INTEGER CHO_LREAD

      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J

C     Initializations.
C     ----------------

      ISHLCD = INDEX(ISHLC,ISHLD)
      ISHLAB = INDEX(ISHLA,ISHLB)
      ERRMAX = -1.0D12
      ERRMIN =  1.0D12
      ERRRMS = 0.0D0
      NCMP   = 0
      LCDABT = NCD*NAB

      IF (ISHLC .EQ. ISHLD) THEN
         NCDL = NBSTSH(ISHLC)*(NBSTSH(ISHLC) + 1)/2
      ELSE
         NCDL = NBSTSH(ISHLC)*NBSTSH(ISHLD)
      END IF
      IF (ISHLA .EQ. ISHLB) THEN
         NABL = NBSTSH(ISHLA)*(NBSTSH(ISHLA) + 1)/2
      ELSE
         NABL = NBSTSH(ISHLA)*NBSTSH(ISHLB)
      END IF
      IF (NCDL .GT. NCD) CALL CHO_QUIT('NCD error in '//SECNAM,104)
      IF (NABL .GT. NAB) CALL CHO_QUIT('NAB error in '//SECNAM,104)
      IF (NAB.LT.1 .OR. NCD.LT.1) RETURN

C     Get a copy of XINT.
C     -------------------

      KXINT = 1
      KEND0 = KXINT + LCDABT
      LWRK0 = LWORK - KEND0
      IF (LWRK0 .LE. 0) THEN
         CALL CHO_QUIT('Insufficient memory in '//SECNAM//' [0]',101)
      END IF

      CALL DCOPY(LCDABT,XINT,1,WORK(KXINT),1)

C     Start symmetry loop.
C     --------------------

      DO ISYM = 1,NSYM

         NUMCD = NNBSTRSH(ISYM,ISHLCD,2)
         NUMAB = NNBSTRSH(ISYM,ISHLAB,2)

         IF (NUMCD.GT.0 .AND. NUMAB.GT.0 .AND. NUMCHO(ISYM).GT.0) THEN

C           Allocate space for integrals and for Cholesky reading.
C           ------------------------------------------------------

            LENINT = NUMCD*NUMAB
            LREAD  = 2*NNBSTR(ISYM,1) + 1
            LVEC1  = NNBSTR(ISYM,2)

            KINT  = KEND0
            KREAD = KINT  + LENINT
            KVEC1 = KREAD + LREAD
            KEND1 = KVEC1 + LVEC1
            LWRK1 = LWORK - KEND1 + 1

            IF (LWRK1 .LE. 0) THEN
               CALL CHO_QUIT('Insufficient memory in '//SECNAM,104)
            END IF

C           Initialize integral array.
C           --------------------------

            CALL CHO_DZERO(WORK(KINT),LENINT)

C           Set up batch over Cholesky vectors.
C           -----------------------------------

            MINM = NUMCD + NUMAB
            NVEC = MIN(LWRK1/MINM,NUMCHO(ISYM))
            IF (NVEC .LT. 1) THEN
               CALL CHO_QUIT('Batch problem in '//SECNAM,104)
            END IF
            NBATCH = (NUMCHO(ISYM) - 1)/NVEC + 1

C           Start batch loop.
C           -----------------

            DO IBATCH = 1,NBATCH

               IF (IBATCH .EQ. NBATCH) THEN
                  NUMV = NUMCHO(ISYM) - NVEC*(NBATCH - 1)
               ELSE
                  NUMV = NVEC
               END IF
               JVEC1 = NVEC*(IBATCH - 1) + 1

               KCHOCD = KEND1
               KCHOAB = KCHOCD + NUMCD*NUMV
               KEND2  = KCHOAB + NUMAB*NUMV
               LWRK2  = LWORK  - KEND2 + 1

               IF (LWRK2 .LT. 0) THEN
                  CALL CHO_QUIT('Batch error in '//SECNAM,4)
               END IF

C              Read vectors.
C              -------------

               DO IVEC = 1,NUMV
                  JVEC = JVEC1 + IVEC - 1
                  CALL CHO_GETVEC(WORK(KVEC1),LVEC1,1,JVEC,ISYM,
     &                            INFRED,INFVEC,
     &                            INDRED,INDRSH,IIBSTRSH,NNBSTRSH,
     &                            ISCR,MRED,MVEC,
     &                            MSYM,MMBSTRT,MMSHL,LISCR,
     &                            WORK(KREAD),LREAD)
                  KOFF1 = KVEC1  + IIBSTRSH(ISYM,ISHLCD,2)
                  KOFF2 = KCHOCD + NUMCD*(IVEC - 1)
                  CALL DCOPY(NUMCD,WORK(KOFF1),1,WORK(KOFF2),1)
                  KOFF1 = KVEC1  + IIBSTRSH(ISYM,ISHLAB,2)
                  KOFF2 = KCHOAB + NUMAB*(IVEC - 1)
                  CALL DCOPY(NUMAB,WORK(KOFF1),1,WORK(KOFF2),1)
               END DO

C              Calculate contribution.
C              -----------------------

               CALL DGEMM('N','T',NUMCD,NUMAB,NUMV,
     &                    1.0D0,WORK(KCHOCD),NUMCD,WORK(KCHOAB),NUMAB,
     &                    1.0D0,WORK(KINT),NUMCD)

            END DO

C           Subtract contribution from full shell pair.
C           -------------------------------------------

            DO IAB = 1,NUMAB
               JAB = IIBSTR(ISYM,2) + IIBSTRSH(ISYM,ISHLAB,2) + IAB
               KAB = INDRED(INDRED(JAB,2),1)
               DO ICD = 1,NUMCD
                  JCD = IIBSTR(ISYM,2) + IIBSTRSH(ISYM,ISHLCD,2) + ICD
                  KCD = INDRED(INDRED(JCD,2),1)
                  ICDAB = KINT  + NUMCD*(IAB - 1) + ICD - 1
                  KCDAB = KXINT + NCD*(KAB - 1)   + KCD - 1
                  WORK(KCDAB) = WORK(KCDAB) - WORK(ICDAB)
               END DO
            END DO

         END IF

      END DO

C     Compare full shell pair.
C     ------------------------

      DO KAB = 1,NAB
         DO KCD = 1,NCD
            KCDAB = KXINT + NCD*(KAB - 1) + KCD - 1
            DIFF  = WORK(KCDAB)
            NCMP  = NCMP + 1
            IF (NCMP .EQ. 1) THEN
               ERRMAX = DIFF
               ERRMIN = DIFF
            ELSE
               IF (ABS(DIFF) .GT. ABS(ERRMAX)) THEN
                  ERRMAX = DIFF
               END IF
               IF (ABS(DIFF) .LT. ABS(ERRMIN)) THEN
                  ERRMIN = DIFF
               END IF
            END IF
            ERRRMS = ERRRMS + DIFF*DIFF
         END DO
      END DO

      END
C  /* Deck cho_inp */
      SUBROUTINE CHO_INP(DFONLY,LUNIT,LUOUT)
C
C     Purpose: If DFONLY, set defaults only.
C              Else, read and process input for Cholesky decomposition
C              from unit LUNIT. LUOUT is the unit of the output file
C              which is stored internally in the Cholesky program as
C              LUPRI (in cholesky.h).
C
#include <implicit.h>
      LOGICAL DFONLY
#include <cholesky.h>
#include <choscf.h>

      CHARACTER*7 SECNAM
      PARAMETER (SECNAM = 'CHO_INP')

      PARAMETER (NOPTION = 35, LOPTION = 50)  ! #options
      CHARACTER*(LOPTION) OPTION(NOPTION)

      DATA OPTION /'set decomposition threshold                       ',
     &             'set print level                                   ',
     &             'set buffer length for diagonal calculation        ',
     &             'set threshold for discarding initial diag. elem.  ',
     &             'set damping factor for first reduced set          ',
     &             'set damping factor for later reduced sets         ',
     &             'set span                                          ',
     &             'set minimum number of qualifieds for decomposition',
     &             'set maximum number of shell pair distributions    ',
     &             'turn on diag. screening (using damping factors)   ',
     &             'turn off diag. screening                          ',
     &             'set qualification algorithm                       ',
     &             'set threshold for zeroing negative diagonals      ',
     &             'set threshold for warning about negative diagonals',
     &             'set threshold for shutdown due to negative diag.  ',
     &             'check configuration only                          ',
     &             'check all integrals                               ',
     &             'diagonal restart                                  ',
     &             'decomposition restart with default restart model  ',
     &             'set restart model                                 ',
     &             'set maximum #qualifieds per symmetry              ',
     &             'set maximum #Cholesky vectors per symmetry        ',
     &             'set max. #reduced sets (i.e., integral passes)    ',
     &             'check specified #columns                          ',
     &             'minimal integral check                            ',
     &             'use abs. value to set up reduced sets (debug)     ',
     &             'do not use abs. value to set up red. sets (debug) ',
     &             'turn on tracing of negative diagonals (debug)     ',
     &             'set model for reading Cholesky vectors            ',
     &             'one-step decomposition                            ',
     &             'two-step decomposition                            ',
     &             'naive decomposition                               ',
     &             'set buffer length used by integral program        ',
     &             'decompose density in SCF: set threshold           ',
     &             'unused                                            '/


C     Set output unit.
C     ----------------

      LUPRI = LUOUT

C     Set defaults.
C     -------------

      CHO_DECALG = 1     ! default decomposition algorithm
      THRDEF  = 1.0D-8   ! default decomposition threshold
      THRCOM  = THRDEF   ! actual threshold
      IPRINT  = 0        ! print level
      LBUF    = 1000000  ! buffer length for initial diagonal
      THRDIAG = 0.0D0    ! threshold for initial diag. screening
      XLDIAG  = 0.0D0    ! just an initialization
      DAMP(1) = -1.0D9   ! damping for initial screening (neg.=>generic)
      DAMP(2) = -1.0D9   ! damping for later screenings (neg.=>generic)
      SPAN    = 1.0D-2   ! span factor
      MINQUAL = 50       ! min. #qualifieds needed for proceeding
      MAXQUAL = MIN(75,CHO_MAXQUAL) ! max. #qual. needed for proceeding
      MXSHPR  = 0        ! max. #sh. pairs for proceeding (0=>generic)
      SCDIAG  = .TRUE.   ! screen diagonal during decom.
      IALQUA  = 2        ! qualification algorithm
      THRNEG  = -1.0D-13 ! diag<THRNEG => diag=0
      WARNEG  = -1.0D-10 ! diag<WARNEG => diag=0, issue warning
      TOONEG  = -1.0D-8  ! diag<TOONEG => shutdown
      CHKONLY = .FALSE.  ! flag for "check input only"
      CHO_INTCHK = .FALSE. ! check integrals after decomposition
      CHO_MINCHK = .FALSE. ! minimal integral check
      RSTDIA  = .FALSE.  ! flag for diagonal restart
      RSTCHO  = .FALSE.  ! flag for decomposition restart
      MODRST  = -1       ! default restart model (if restart at all)
      MAXVEC  = 0        ! max. #Cholesky vectors (0=>generic)
      MAXRED  = 0        ! max. #reduced sets (0=>generic)
      NCOL_CHK= 0        ! #columns to check (0=>all if requested)
      CHO_USEABS = .TRUE.  ! use abs. value to set up reduced sets
      CHO_TRCNEG = .FALSE. ! tracing of negative diagonals
      CHO_IOVEC  = 1     ! i/o model used when reading Cholesky vectors
      LERIBUF = 250000   ! buffer length in integral program.

C     Set defaults for Cholesky SCF.
C     ------------------------------

      CALL SET_CHOSCFDEF()

C     Return if defaults only.
C     ------------------------

      IF (DFONLY) GO TO 2000

C     Loop through keyword input.
C     ---------------------------

 1000 CONTINUE

C        Read next keyword and get corresponding ID.
C        -------------------------------------------

         IDKEY = 0
         CALL CHO_DAL_GETKEY(LUNIT,OPTION,LOPTION,NOPTION,IDKEY,LUPRI)

         IF (IDKEY.GE.1 .AND. IDKEY.LE.NOPTION) THEN  ! key found

C           Branch for further processing.
C           ------------------------------

            GOTO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &            21,22,23,24,25,26,27,28,29,30,31,32,33,34,35),
     &           IDKEY

C           If this section is executed, it's most likely a bug.
C           ----------------------------------------------------

            WRITE(LUPRI,*) 'IDKEY = ',IDKEY,' is formally legal, but:'
            WRITE(LUPRI,*) 'Did you forget to change the ',
     &                     'computed goto in ',SECNAM,'?'
            CALL CHO_QUIT('Illegal address in computed GOTO in '
     &                    //SECNAM,105)

C           Read decomposition threshold.
C           -----------------------------

    1       CONTINUE
               READ(LUNIT,*) THRCOM
            GO TO 1000

C           Read print level.
C           -----------------

    2       CONTINUE
               READ(LUNIT,*) IPRINT
            GO TO 1000

C           Read buffer length for diagonal calculation.
C           --------------------------------------------

    3       CONTINUE
               READ(LUNIT,*) LBUF
            GO TO 1000

C           Read threshold for discarding initial diagonal elements.
C           --------------------------------------------------------

    4       CONTINUE
               READ(LUNIT,*) THRDIAG
            GO TO 1000

C           Read damping factor for 1st reduced set.
C           ----------------------------------------

    5       CONTINUE
               READ(LUNIT,*) DAMP(1)
            GO TO 1000

C           Read damping factor for later reduced set.
C           ------------------------------------------

    6       CONTINUE
               READ(LUNIT,*) DAMP(2)
            GO TO 1000

C           Read span factor.
C           -----------------

    7       CONTINUE
               READ(LUNIT,*) SPAN
            GO TO 1000

C           Read minimum number of qualifieds for decomposition.
C           ----------------------------------------------------

    8       CONTINUE
               READ(LUNIT,*) MINQUAL
            GO TO 1000

C           Read maximum number of shell pair distributions that should
C           be calculated before passing to decomposition.
C           -----------------------------------------------------------

    9       CONTINUE
               READ(LUNIT,*) MXSHPR
            GO TO 1000

C           Turn on screening.
C           ------------------

   10       CONTINUE
               SCDIAG = .TRUE.
            GO TO 1000

C           Turn off screening.
C           -------------------

   11       CONTINUE
               SCDIAG = .FALSE.
            GO TO 1000

C           Read algorithm ID for qualification procedure,
C           1: qualify from 1 to max.
C           2: qualify the largest qualifiables
C           ----------------------------------------------

   12       CONTINUE
               READ(LUNIT,*) IALQUA
            GO TO 1000

C           Read threshold for zeroing (small) negative diagonals.
C           Diag < THRNEG => Diag = 0.
C           ------------------------------------------------------

   13       CONTINUE
               READ(LUNIT,*) THRNEG
            GO TO 1000

C           Read threshold for warning about negative diagonals.
C           Diag < WARNEG => Diag = 0, but issue warning.
C           ----------------------------------------------------

   14       CONTINUE
               READ(LUNIT,*) WARNEG
            GO TO 1000

C           Read threshold for shutdown due to negative diagonals.
C           Diag < TOONEG => shutdown.
C           ------------------------------------------------------

   15       CONTINUE
               READ(LUNIT,*) TOONEG
            GO TO 1000

C           Check configuration only.
C           -------------------------

   16       CONTINUE
               CHKONLY = .TRUE.
            GO TO 1000

C           Check all integrals after decomposition.
C           ----------------------------------------

   17       CONTINUE
               CHO_INTCHK = .TRUE.
               CHO_MINCHK = .FALSE.
               NCOL_CHK   = 0  ! => all
            GO TO 1000

C           Restart diagonal.
C           -----------------

   18       CONTINUE
               RSTDIA = .TRUE.
            GO TO 1000

C           Restart decomposition using default model.
C           ------------------------------------------

   19       CONTINUE
               IF (.NOT. RSTCHO) THEN
                  RSTCHO = .TRUE.
                  MODRST = -1  ! use configuration from restart file
               END IF
            GO TO 1000

C           Read restart model and set restart.
C           -----------------------------------

   20       CONTINUE
               RSTCHO = .TRUE.
               READ(LUNIT,*) MODRST
               IF (MODRST .LT. 0) THEN
                  MODRST = -1  ! use configuration from restart file
               ELSE IF (MODRST .GT. 0) THEN
                  MODRST = 1   ! use configuration from input
               END IF
            GO TO 1000

C           Read max. #qualifieds per symmetry.
C           -----------------------------------

   21       CONTINUE
               READ(LUNIT,*) MAXQUAL
            GO TO 1000

C           Read max. #Cholesky vectors per symmetry.
C           -----------------------------------------

   22       CONTINUE
               READ(LUNIT,*) MAXVEC
            GO TO 1000

C           Read max. #reduced sets.
C           ------------------------

   23       CONTINUE
               READ(LUNIT,*) MAXRED
            GO TO 1000

C           Read number of integral columns to check.
C           -----------------------------------------

   24       CONTINUE
               CHO_INTCHK = .TRUE.
               CHO_MINCHK = .FALSE.
               READ(LUNIT,*) NCOL_CHK
            GO TO 1000

C           Minimal integral check.
C           -----------------------

   25       CONTINUE
               CHO_INTCHK = .TRUE.
               CHO_MINCHK = .TRUE.
            GO TO 1000

C           Use abs. value to set up reduced sets.
C           --------------------------------------

   26       CONTINUE
               CHO_USEABS = .TRUE.
            GO TO 1000

C           Do not use abs. value to set up reduced sets.
C           ---------------------------------------------

   27       CONTINUE
               CHO_USEABS = .FALSE.
            GO TO 1000

C           Turn on tracing of negative diagonals.
C           --------------------------------------

   28       CONTINUE
               CHO_TRCNEG = .TRUE.
            GO TO 1000

C           Read I/O model used for reading Cholesky vectors.
C           -------------------------------------------------

   29       CONTINUE
               READ(LUNIT,*) CHO_IOVEC
               IF (CHO_IOVEC .LT. 1) THEN
                  CHO_IOVEC = 1
               ELSE IF (CHO_IOVEC .GT. 1) THEN
                  CHO_IOVEC = 2
               END IF
            GO TO 1000

C           One-step decomposition.
C           -----------------------

   30       CONTINUE
               CHO_DECALG = 1
            GO TO 1000

C           Two-step decomposition.
C           -----------------------

   31       CONTINUE
               CHO_DECALG = 2
            GO TO 1000

C           Naive decomposition.
C           --------------------

   32       CONTINUE
               CHO_DECALG = 3
            GO TO 1000

C           Buffer length transferred to integral program.
C           ----------------------------------------------

   33       CONTINUE
               READ(LUNIT,*) LLL
               IF (LLL .GT. 0) THEN
                  LERIBUF = LLL
               END IF
            GO TO 1000

C           Decompose density in SCF; set threshold.
C           ----------------------------------------

   34       CONTINUE
               CCMODSK = .FALSE.
               READ(LUNIT,*) THRDC
            GO TO 1000

C           Unused.
C           -------

   35       CONTINUE
            GO TO 1000

C           Unused.
C           -------

         ELSE IF (IDKEY .EQ. NOPTION+1) THEN  ! end of input
            GO TO 2000
         ELSE IF (IDKEY .EQ. -5) THEN ! some internal error
            GO TO 1300
         ELSE IF (IDKEY .EQ. -1) THEN ! keyword error
            GO TO 1200
         ELSE  ! catch all out-of-bound errors
            GO TO 1100
         END IF

C     Normal exit point.
C     ------------------

 2000 RETURN

C     Error branches.
C     ---------------

 1100 CONTINUE ! IDKEY out of bounds
      WRITE(LUPRI,*) SECNAM,': IDKEY out of bounds: ',IDKEY
      GO TO 9999

 1200 CONTINUE ! keyword error
      WRITE(LUPRI,*) SECNAM,': keyword error detected, IDKEY = ',IDKEY
      GO TO 9999

 1300 CONTINUE ! some internal error
      WRITE(LUPRI,*) SECNAM,': internal error detected, IDKEY = ',IDKEY
      GO TO 9999

 9999 CALL CHO_QUIT('Error in '//SECNAM,105)
      END
C  /* Deck cho_chkconf */
      SUBROUTINE CHO_CHKCONF(NCONFL,VERBOSE)
C
C     Purpose: check configuration, return the number of errors NCONFL.
C
      IMPLICIT NONE
      INTEGER NCONFL
      LOGICAL VERBOSE
#include <cholesky.h>
#include <choorb.h>

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CHO_CHKCONF')

      LOGICAL REPORT
      INTEGER NNN, MMM, INEGRR
      REAL*8 XLBUF, XMBUF

C     Initialize.
C     -----------

      NCONFL = 0
      REPORT = VERBOSE

C     Check that output unit is appropriately set.
C     (Upper bound not checked, as it may change over time.)
C     ------------------------------------------------------

      IF (LUPRI .LT. 1) THEN
         NCONFL = NCONFL + 1
         REPORT = .FALSE.
      END IF

C     Check decomposition algorithm.
C     ------------------------------

      IF (CHO_DECALG.LT.1 .OR. CHO_DECALG.GT.3) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,*) 'Illegal decomposition algorithm, ',
     &                     'CHO_DECALG = ',CHO_DECALG
            WRITE(LUPRI,*) 'Assuming default algorithm, ',
     &                     'CHO_DECALG = 1'
         END IF
         CHO_DECALG = 1
      END IF
      IF (CHO_DECALG .EQ. 2) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,*) 'Two-step decomposition has not been ',
     &                     'implemented. Reverting to one-step.'
         END IF
         CHO_DECALG = 1
      END IF

C     Check max. number of Cholesky vectors and reduced sets.
C     -------------------------------------------------------

      IF (MAXVEC .LT. 1) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,*) 'Max. number of vectors < 1: ',MAXVEC
         END IF
         NCONFL = NCONFL + 1
      END IF
      IF (MAXRED .LT. 1) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,*) 'Max. number of reduced sets < 1: ',MAXRED
         END IF
         NCONFL = NCONFL + 1
      END IF

C     Check qualification algorithm.
C     ------------------------------

      IF (IALQUA .LT. 1) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,*) 'Qualification algorithm reset from ',
     &                     IALQUA,' to 1.'
         END IF
         IALQUA = 1
      ELSE IF (IALQUA .GT. NALQUA) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,*) 'Qualification algorithm reset from ',
     &                     IALQUA,' to ',NALQUA
         END IF
         IALQUA = NALQUA
      END IF

C     Decomposition threshold.
C     ------------------------

      IF (THRCOM .LT. 0.0D0) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,*) 'Decomposition threshold not positive: ',
     &                     THRCOM,' (apparent default value: ',THRDEF,
     &                     ')'
         END IF
         NCONFL = NCONFL + 1
      END IF

C     Buffer length.
C     --------------

      IF (LBUF .LT. 1) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,*) 'Buffer length < 1: ',LBUF
         END IF
         NCONFL = NCONFL + 1
      ELSE
         XLBUF = DBLE(LBUF)
         XMBUF = DBLE(NBAST)*(DBLE(NBAST)+1.0D0)/2.0D0
         IF (XLBUF .GT. XMBUF) THEN ! make sure LBUF is not too large
            LBUF = NINT(XMBUF)
         END IF
      END IF

C     Threshold for discarding elements of initial diagonal.
C     ------------------------------------------------------

      IF (THRDIAG .GT. THRCOM) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,*) 'Threshold for discarding initial diagonals',
     &                     ': ',THRDIAG
            WRITE(LUPRI,*) 'is larger than decomposition threshold: ',
     &                     THRCOM
         END IF
         NCONFL = NCONFL + 1
      END IF

C     Damping factors.
C     ----------------

      IF (DAMP(1) .LT. 1.0D0) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,*) 'First damping factor < 1: ',DAMP(1)
         END IF
         NCONFL = NCONFL + 1
      END IF
      IF (DAMP(2) .LT. 1.0D0) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,*) 'Second damping factor < 1: ',DAMP(2)
         END IF
         NCONFL = NCONFL + 1
      END IF

C     Span factor.
C     ------------

      IF (SPAN .GT. 1.0D0) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,*) 'Span factor > 1: ',SPAN
         END IF
         NCONFL = NCONFL + 1
      END IF

C     Max. #shell pairs allowed before proceeding to deco.
C     Max. #qualifieds allowed to proceed to decomposition.
C     Min. #qualifieds needed to proceed to decomposition.
C     -----------------------------------------------------

      NNN = NSHELL*(NSHELL + 1)/2
      IF (MXSHPR.LT.0 .OR. MXSHPR.GT.NNN) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,*) 'max. #shell pairs allowed ',
     &                     'before proceeding to deco. is: ',MXSHPR
            WRITE(LUPRI,*) 'resetting generic: 0'
         END IF
         MXSHPR = 0
      END IF
      IF (CHO_MAXQUAL .LT. 1) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,*) 'error: CHO_MAXQUAL = ',CHO_MAXQUAL
         END IF
         CALL CHO_QUIT('Install config error detected in '//SECNAM,103)
      END IF
      IF (MAXQUAL.LT.1 .OR. MAXQUAL.GT.CHO_MAXQUAL) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,*) 'max. #qualified allowed ',
     &                     'before proceeding to deco. is: ',MAXQUAL
            WRITE(LUPRI,*) 'resetting to max.: ',CHO_MAXQUAL
         END IF
         MAXQUAL = CHO_MAXQUAL
      END IF
      MMM = NSYM*MAXQUAL
      IF (MINQUAL.LT.1 .OR. MINQUAL.GT.MMM) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,*) 'max. #qualified needed ',
     &                     'before proceeding to deco. is: ',MINQUAL
            WRITE(LUPRI,*) 'resetting to max.: ',MMM
         END IF
         MINQUAL = MMM
      END IF

C     Handling of negative diagonals.
C     -------------------------------

      INEGRR = 0
      IF (THRNEG .GT. 0.0D0) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,*) 'Threshold for zeroing neg. diag. > 0: ',
     &                     THRNEG
         END IF
         INEGRR = INEGRR + 1
         NCONFL = NCONFL + 1
      END IF
      IF (WARNEG .GT. 0.0D0) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,*) 'Threshold for warning about neg. diag. ',
     &                     ' > 0: ',WARNEG
         END IF
         INEGRR = INEGRR + 1
         NCONFL = NCONFL + 1
      END IF
      IF (TOONEG .GT. 0.0D0) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,*) 'Threshold for shutdown due to neg. diag. ',
     &                     ' > 0: ',TOONEG
         END IF
         INEGRR = INEGRR + 1
         NCONFL = NCONFL + 1
      END IF
      IF (INEGRR .EQ. 0) THEN
         IF (THRNEG .LT. WARNEG) THEN
            IF (REPORT) THEN
               WRITE(LUPRI,*) 'Threshold for zeroing neg. diag. > ',
     &                       'threshold for warning about neg. diag.: ',
     &                        THRNEG,WARNEG
            END IF
            NCONFL = NCONFL + 1
         END IF
         IF (WARNEG .LT. TOONEG) THEN
            IF (REPORT) THEN
               WRITE(LUPRI,*)
     &                     'Threshold for warning about neg. diag. > ',
     &                     'threshold for shutdown due to neg. diag.: ',
     &                     WARNEG,TOONEG
            END IF
            NCONFL = NCONFL + 1
         END IF
      END IF

      END
C  /* Deck cho_iodiag */
      SUBROUTINE CHO_IODIAG(DIAG,IOPT)
C
C     Purpose: write/read a copy of diagonal to disk (1st reduced set).
C              The file is opened and closed here.
C              IOPT=1: write
C              IOPT=2: read
C
#include <implicit.h>
      DIMENSION DIAG(*)
#include <cholesky.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_IODIAG')

      CHARACTER*7 FNAME
      PARAMETER (FNAME = 'CHODIAG')

      LUNIT = -1
      IDUM  = -1
      IF (IOPT .EQ. 1) THEN
         CALL GPOPEN(LUNIT,FNAME,'UNKNOWN','SEQUENTIAL',
     &               'UNFORMATTED',IDUM,.FALSE.)
      ELSE IF (IOPT .EQ. 2) THEN
         CALL GPOPEN(LUNIT,FNAME,'OLD','SEQUENTIAL',
     &               'UNFORMATTED',IDUM,.FALSE.)
      ELSE
         CALL CHO_QUIT('IOPT error in '//SECNAM,104)
      END IF

C     I/O section.
C     ------------

      IF (IOPT .EQ. 1) THEN
         WRITE(LUNIT) (DIAG(I),I=1,NNBSTRT(1))
      ELSE IF (IOPT .EQ. 2) THEN
         READ(LUNIT) (DIAG(I),I=1,NNBSTRT(1))
      ELSE
         CALL CHO_QUIT('IOPT error in '//SECNAM,104)
      END IF

C     Close the file.
C     ---------------

      CALL GPCLOSE(LUNIT,'KEEP')

      END
C  /* Deck cho_rstd_getind1 */
      SUBROUTINE CHO_RSTD_GETIND1(IIBSTRSH,NNBSTRSH,MSYM,MMSHL)
C
C     Purpose: read some index arrays for diagonal restart.
C
#include <implicit.h>
      INTEGER IIBSTRSH(MSYM,MMSHL,3), NNBSTRSH(MSYM*MMSHL*3)
#include <cholesky.h>

      CHARACTER*16 SECNAM
      PARAMETER (SECNAM = 'CHO_RSTD_GETIND1')

C     Check dimensions.
C     -----------------

      NNSHL = NSHELL*(NSHELL + 1)/2
      IF (MSYM.NE.NSYM .OR. MMSHL.NE.NNSHL) THEN
         CALL CHO_QUIT('Dimension error in '//SECNAM,104)
      END IF

C     Read info.
C     ----------

      IADR = 1
      REWIND(LURED)
      READ(LURED,REC=IADR) (NNBSTRSH(I),I=1,NSYM*NNSHL)

C     Set up IIBSTRSH etc.
C     --------------------

      IRED = 1
      CALL CHO_SETREDIND(IIBSTRSH,NNBSTRSH,MSYM,MMSHL,IRED)

      END
C  /* Deck cho_rstd_getind2 */
      SUBROUTINE CHO_RSTD_GETIND2(INDRED,INDRSH,MMBSTRT)
C
C     Purpose: read mapping arrays for diagonal restart.
C
#include <implicit.h>
      INTEGER INDRED(MMBSTRT,3), INDRSH(MMBSTRT)
#include <cholesky.h>

      CHARACTER*16 SECNAM
      PARAMETER (SECNAM = 'CHO_RSTD_GETIND2')

C     Check dimension.
C     ----------------

      IF (MMBSTRT .NE. NNBSTRT(1)) THEN
         CALL CHO_QUIT('Dimension error in '//SECNAM,104)
      END IF

C     Read info.
C     ----------

      READ(LURED,REC=1) (IDUM,I=1,NSYM*NSHELL*(NSHELL+1)/2),
     &                  (INDRED(I,1),I=1,NNBSTRT(1))
      REWIND(LURSH)
      READ(LURSH) (INDRSH(I),I=1,NNBSTRT(1))

      END
C  /* Deck cho_getrstc */
      SUBROUTINE CHO_GETRSTC(INFRED,NR,INFVEC,NV1,NV2,NV3)
C
C     Purpose: read and check decomposition restart info.
C
#include <implicit.h>
      INTEGER INFRED(NR), INFVEC(NV1,NV2,NV3)
#include <cholesky.h>

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CHO_GETRSTC')

C     Read restart file, populating restart common block.
C     ---------------------------------------------------

      IFAIL = 0
      CALL CHO_RDRSTC(INFRED,NR,INFVEC,NV1,NV2,NV3,IFAIL)
      IF (IFAIL .NE. 0) THEN
         WRITE(LUPRI,'(A,A)')
     &   SECNAM,': error reading decomposition restart file.'
         WRITE(LUPRI,'(A,A,I10)')
     &   SECNAM,': return code from reading routine:',IFAIL
         CALL CHO_QUIT('Error reading decomposition restart file',104)
      END IF

C     Check system info.
C     ------------------

      IFAIL = 0
      CALL CHO_RSTMOL(IFAIL)
      IF (IFAIL .NE. 0) THEN
         WRITE(LUPRI,'(A,A)')
     &   SECNAM,': decomposition restart failure.'
         CALL CHO_QUIT('Decomposition restart failure in '//SECNAM,105)
      END IF

C     Check decomposition configuration.
C     ----------------------------------

      IFAIL = 0
      CALL CHO_RSTCNF(IFAIL)
      IF (IFAIL .NE. 0) THEN
         WRITE(LUPRI,'(A,A,I6,A)')
     &   SECNAM,':',IFAIL,' configuration discrepancies detected.'
         IF (MODRST .EQ. -1) THEN
            WRITE(LUPRI,'(A)')
     &      'Recovery: using configuration from restart file.'
            CALL CHO_RESETCNF
         ELSE IF (MODRST .EQ. 0) THEN
            WRITE(LUPRI,'(A)')
     &      'Recovery: none, program stops.'
            CALL CHO_QUIT('Restart configuration error in '//SECNAM,105)
         ELSE IF (MODRST .EQ. 1) THEN
            WRITE(LUPRI,'(A)')
     &      'Recovery: using input configuration.'
         ELSE
            WRITE(LUPRI,'(A,A,I6,A)')
     &      SECNAM,': restart model,',MODRST,', not recognized.'
            CALL CHO_QUIT('Error in '//SECNAM,103)
         END IF
      END IF

      END
C  /* Deck cho_wrrstc */
      SUBROUTINE CHO_WRRSTC(INFRED,NR,INFVEC,NV1,NV2,NV3,
     &                      IPASS,INTMAP,MMSHL)
C
C     Purpose: write decomposition restart info for integral pass IPASS.
C
C     NB!!!  The restart files are assumed open on entry.
C
#include <implicit.h>
      INTEGER INFRED(NR), INFVEC(NV1,NV2,NV3)
      INTEGER INTMAP(MMSHL)
#include <cholesky.h>
#include <choorb.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_WRRSTC')

      INTEGER*8 IADR

      NNSHL = NSHELL*(NSHELL + 1)/2
      IF (MMSHL .NE. NNSHL) THEN
         WRITE(LUPRI,'(A,A,I10,A,I10)')
     &   SECNAM,': MMSHL:',MMSHL,' expected:',NNSHL
         CALL CHO_QUIT('Dimension error -1 in '//SECNAM,104)
      END IF
      IF (NR .NE. MAXRED) THEN
         WRITE(LUPRI,'(A,A,I10,A,I10)')
     &   SECNAM,': NR:',NR,' expected:',MAXRED
         CALL CHO_QUIT('Dimension error 0 in '//SECNAM,104)
      END IF
      IF (NV1 .NE. MAXVEC) THEN
         WRITE(LUPRI,'(A,A,I10,A,I10)')
     &   SECNAM,': NV1:',NV1,' expected:',MAXVEC
         CALL CHO_QUIT('Dimension error 1 in '//SECNAM,104)
      END IF
      N2 = 3
      IF (NV2 .NE. N2) THEN
         WRITE(LUPRI,'(A,A,I10,A,I10)')
     &   SECNAM,': NV2:',NV2,' expected:',N2
         CALL CHO_QUIT('Dimension error 2 in '//SECNAM,104)
      END IF
      IF (NV3 .NE. NSYM) THEN
         WRITE(LUPRI,'(A,A,I10,A,I10)')
     &   SECNAM,': NV3:',NV3,' expected:',NSYM
         CALL CHO_QUIT('Dimension error 3 in '//SECNAM,104)
      END IF

      IF (SCDIAG) THEN
         ISCDIA = 1
      ELSE
         ISCDIA = 0
      END IF
      REWIND(LURST)
      WRITE(LURST) IPASS,NSYM,NSHELL,ISCDIA,
     &             THRCOM,THRDIAG,DAMP(1),DAMP(2),SPAN,
     &             THRNEG,WARNEG,TOONEG
      WRITE(LURST) (NBAS(ISYM),ISYM=1,NSYM)
      WRITE(LURST) (NUMCHO(ISYM),ISYM=1,NSYM)
      WRITE(LURST) (INFRED(I),I=1,IPASS)
      DO ISYM = 1,NSYM
         IF (NUMCHO(ISYM) .GT. 0) THEN
            DO J = 1,3
               WRITE(LURST) (INFVEC(I,J,ISYM),I=1,NUMCHO(ISYM))
            END DO
         END IF
      END DO
      REWIND(LUMAP)
      WRITE(LUMAP) (INTMAP(I),I=1,MMSHL)

      END
C  /* Deck cho_rdrstc */
      SUBROUTINE CHO_RDRSTC(INFRED,NR,INFVEC,NV1,NV2,NV3,IFAIL)
C
C     Purpose: read decomposition restart info and store in common
C              block. If IFAIL != 0 on exit, some error occurred and,
C              most likely, some of the restart info is not
C              defined/initialized.
C
C     NB!!!  The restart file is opened AND closed here...
C
#include <implicit.h>
      INTEGER INFRED(NR), INFVEC(NV1,NV2,NV3)
#include <cholesky.h>
#include <choorb.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_RDRSTC')

C     Set return code.
C     ----------------

      IFAIL = 0

C     Open decomposition restart file.
C     --------------------------------

      LURST = -1
      IDUM  = -1
      CALL GPOPEN(LURST,FRST,'OLD','SEQUENTIAL',
     &            'UNFORMATTED',IDUM,.FALSE.)

C     Read info.
C     ----------

      REWIND(LURST)
      READ(LURST) XNPASS,XNSYM,XNSHELL,ISCDIA,
     &            XTHRCOM,XTHRDIAG,XDAMP(1),XDAMP(2),XSPAN,
     &            XTHRNEG,XWARNEG,XTOONEG
      IF (XNSYM.LT.1 .OR. XNSYM.GT.8) THEN
         WRITE(LUPRI,'(A,A,I10)')
     &   SECNAM,': #irreps from restart file: ',XNSYM
         IFAIL = 1
         GO TO 100
      ELSE
         READ(LURST) (XNBAS(ISYM),ISYM=1,XNSYM)
         READ(LURST) (NUMCHO(ISYM),ISYM=1,XNSYM)
      END IF
      IF (ISCDIA .EQ. 0) THEN
         SCDIAG = .FALSE.
      ELSE IF (ISCDIA .EQ. 1) THEN
         SCDIAG = .TRUE.
      ELSE
         WRITE(LUPRI,'(A,A,I10)')
     &   SECNAM,': integer flag for screening not recognized:',ISCDIA
         IFAIL = 2
         GO TO 100
      END IF
      IF (XNPASS.LT.1 .OR. XNPASS.GT.MAXRED) THEN
         WRITE(LUPRI,'(A,A,I10)')
     &   SECNAM,': #reduced sets in restart:',XNPASS
         IFAIL = 3
         GO TO 100
      ELSE
         READ(LURST) (INFRED(I),I=1,XNPASS)
         IF (INFRED(1) .NE. 1) THEN
            WRITE(LUPRI,'(A,A,I10)')
     &      SECNAM,': disk address of 1st reduced set:',INFRED(1)
            IFAIL = 4
            GO TO 100
         END IF
         LREST = MAXRED - XNPASS
         IF (LREST .GT. 0) CALL CHO_IZERO(INFRED(XNPASS+1),LREST)
      END IF
      DO ISYM = 1,XNSYM
         IF (NUMCHO(ISYM).LT.0 .OR. NUMCHO(ISYM).GT.MAXVEC) THEN
            WRITE(LUPRI,'(A,A,I2,A,I10)')
     &      SECNAM,': #Cholesky vectors (sym.',ISYM,'): ',NUMCHO(ISYM)
            IFAIL = 5
            GO TO 100
         ELSE IF (NUMCHO(ISYM) .EQ. 0) THEN
            NDIM = MAXVEC*3
            CALL CHO_IZERO(INFVEC(1,1,ISYM),NDIM)
         ELSE
            DO J = 1,3
               READ(LURST) (INFVEC(I,J,ISYM),I=1,NUMCHO(ISYM))
               LREST = MAXVEC - NUMCHO(ISYM)
               IF (LREST .GT. 0) THEN
                  CALL CHO_IZERO(INFVEC(NUMCHO(ISYM)+1,J,ISYM),LREST)
               END IF
            END DO
         END IF
      END DO

C     Close restart file.
C     -------------------

  100 IF (IFAIL .NE. 0) THEN  ! failures jump to this point
         WRITE(LUPRI,'(A,A)')
     &   SECNAM,': refusing to read more restart info!'
      END IF

      CALL GPCLOSE(LURST,'KEEP')
      LURST = -1

      END
C  /* Deck cho_resetcnf */
      SUBROUTINE CHO_RESETCNF
C
C     Purpose: reset configuration of decomposition to that read from
C              restart file. Original configuration will be saved in
C              restart common block.
C
      IMPLICIT NONE
#include <cholesky.h>

      CALL CHO_DSWAP(THRCOM,XTHRCOM)
      CALL CHO_DSWAP(THRDIAG,XTHRDIAG)
      CALL CHO_DSWAP(DAMP(1),XDAMP(1))
      CALL CHO_DSWAP(DAMP(2),XDAMP(2))
      CALL CHO_DSWAP(SPAN,XSPAN)
      CALL CHO_DSWAP(THRNEG,XTHRNEG)
      CALL CHO_DSWAP(WARNEG,XWARNEG)
      CALL CHO_DSWAP(TOONEG,XTOONEG)
      CALL CHO_LSWAP(SCDIAG,XSCDIAG)

      END
C  /* Deck cho_lswap */
      SUBROUTINE CHO_LSWAP(X,Y)
C
C     Purpose: swap X and Y, logical version.
C
      IMPLICIT NONE
      LOGICAL  X, Y, TEMP

      TEMP = X
      X    = Y
      Y    = TEMP

      END
C  /* Deck cho_dswap */
      SUBROUTINE CHO_DSWAP(X,Y)
C
C     Purpose: swap X and Y, double precision version.
C
      IMPLICIT NONE
      REAL*8  X, Y, TEMP

      TEMP = X
      X    = Y
      Y    = TEMP

      END
C  /* Deck cho_rstmol */
      SUBROUTINE CHO_RSTMOL(NERR)
C
C     Purpose: check restart molecular info.
C
      IMPLICIT NONE
      INTEGER NERR
#include <cholesky.h>
#include <choorb.h>

      INTEGER ISYM

      NERR = 0

      IF (XNSYM .NE. NSYM) THEN
         WRITE(LUPRI,'(A,I3,A,I3)')
     &   'RESTART ERROR: #irreps from restart file:',XNSYM,
     &   ' Expected:',NSYM
         NERR = NERR + 1
      ELSE
         DO ISYM = 1,NSYM
            IF (XNBAS(ISYM) .NE. NBAS(ISYM)) THEN
               WRITE(LUPRI,'(A,I2,A,I6,A,I6)')
     &         'RESTART ERROR: #basis functions (sym.',ISYM,
     &         ') from restart file:',XNBAS(ISYM),
     &         ' Expected:',NBAS(ISYM)
               NERR = NERR + 1
            END IF
         END DO
      END IF

      IF (XNSHELL .NE. NSHELL) THEN
         WRITE(LUPRI,'(A,I3,A,I3)')
     &   'RESTART ERROR: #shells from restart file:',XNSHELL,
     &   ' Expected:',NSHELL
         NERR = NERR + 1
      END IF

      END
C  /* Deck cho_rstcnf */
      SUBROUTINE CHO_RSTCNF(NERR)
C
C     Purpose: check restart configuration info.
C
      IMPLICIT NONE
      INTEGER NERR
#include <cholesky.h>

      INTEGER I, J

      CHARACTER*3 SWITCH(2)
      DATA SWITCH /' ON','OFF'/

      REAL*8 ERRTOL, ERR
      PARAMETER (ERRTOL = 1.0D-14) 

      NERR = 0

      ERR = ABS(THRCOM - XTHRCOM)
      IF (ERR .GT. ERRTOL) THEN
         WRITE(LUPRI,'(A,D16.8,/,A,D16.8)')
     &   'RESTART: decomposition threshold from restart file: ',XTHRCOM,
     &   '         decomposition threshold from input       : ',THRCOM
         NERR = NERR + 1
      END IF

      ERR = ABS(THRDIAG - XTHRDIAG)
      IF (ERR .GT. ERRTOL) THEN
         WRITE(LUPRI,'(A,D16.8,/,A,D16.8)')
     &   'RESTART: init. diag. screening from restart file: ',XTHRDIAG,
     &   '         init. diag. screening from input       : ',THRDIAG
         NERR = NERR + 1
      END IF

      ERR = ABS(DAMP(1) - XDAMP(1))
      IF (ERR .GT. ERRTOL) THEN
         WRITE(LUPRI,'(A,D16.8,/,A,D16.8)')
     &   'RESTART: 1st screening damping from restart file: ',XDAMP(1),
     &   '         1st screening damping from input       : ',DAMP(1)
         NERR = NERR + 1
      END IF

      ERR = ABS(DAMP(2) - XDAMP(2))
      IF (ERR .GT. ERRTOL) THEN
         WRITE(LUPRI,'(A,D16.8,/,A,D16.8)')
     &   'RESTART: 2nd screening damping from restart file: ',XDAMP(2),
     &   '         2nd screening damping from input       : ',DAMP(2)
         NERR = NERR + 1
      END IF

      IF (SCDIAG .NEQV. XSCDIAG) THEN
         IF (XSCDIAG) THEN
            I = 1
         ELSE
            I = 2
         END IF
         IF (SCDIAG) THEN
            J = 1
         ELSE
            J = 2
         END IF
         WRITE(LUPRI,'(A,A,/,A,A)')
     &   'RESTART: diag. screening from restart file: ',SWITCH(I),
     &   '         diag. screening from input       : ',SWITCH(J)
         NERR = NERR + 1
      END IF

      ERR = ABS(THRNEG - XTHRNEG)
      IF (ERR .GT. ERRTOL) THEN
         WRITE(LUPRI,'(A,D16.8,/,A,D16.8)')
     &   'RESTART: neg. diag. threshold from restart file: ',XTHRNEG,
     &   '         neg. diag. threshold from input       : ',THRNEG
         NERR = NERR + 1
      END IF

      ERR = ABS(WARNEG - XWARNEG)
      IF (ERR .GT. ERRTOL) THEN
         WRITE(LUPRI,'(A,D16.8,/,A,D16.8)')
     &   'RESTART: neg. diag. warn thr. from restart file: ',XWARNEG,
     &   '         neg. diag. warn thr. from input       : ',WARNEG
         NERR = NERR + 1
      END IF

      ERR = ABS(TOONEG - XTOONEG)
      IF (ERR .GT. ERRTOL) THEN
         WRITE(LUPRI,'(A,D16.8,/,A,D16.8)')
     &   'RESTART: too neg. diag. thr. from restart file: ',XTOONEG,
     &   '         too neg. diag. thr. from input       : ',TOONEG
         NERR = NERR + 1
      END IF

      ERR = ABS(SPAN - XSPAN)
      IF (ERR .GT. ERRTOL) THEN
         WRITE(LUPRI,'(A,D16.8,/,A,D16.8)')
     &   'RESTART: span factor from restart file: ',XSPAN,
     &   '         span factor from input       : ',SPAN
         NERR = NERR + 1
      END IF

      END
C  /* Deck cho_inimap */
      SUBROUTINE CHO_INIMAP(INTMAP,MMSHL)
C
C     Purpose: initialize integral shell pair calculation map.
C
C     NB!!!!! file is assumed open (restart only)
C
#include <implicit.h>
      INTEGER INTMAP(MMSHL)
#include <cholesky.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_INIMAP')

      INTEGER*8 IADR

C     Check dimension.
C     ----------------

      NNSHL = NSHELL*(NSHELL + 1)/2
      IF (MMSHL .NE. NNSHL) THEN
         WRITE(LUPRI,'(A,A,I7,A,I7)')
     &   SECNAM,': dimension error: MMSHL = ',MMSHL,
     &   ' expected: ',NNSHL
         CALL CHO_QUIT('Error in '//SECNAM,104)
      END IF

C     If restart, read the map from disk.
C     Else, zero all elements.
C     -----------------------------------

      IF (RSTCHO) THEN
         REWIND(LUMAP)
         READ(LUMAP) (INTMAP(I),I=1,NNSHL)
      ELSE
         CALL CHO_IZERO(INTMAP,NNSHL)
      END IF

      END
C  /* Deck cho_setaddr */
      SUBROUTINE CHO_SETADDR(INFRED,INFVEC,IIBSTRSH,NNBSTRSH,INDRED,
     &                       INDRSH,MRED,MVEC,MSYM,MMSHL,MMBSTRT)
C
C     Purpose: set first disk addresses for reduced set info and
C              vectors.
C
#include <implicit.h>
      INTEGER INFRED(MRED), INFVEC(MVEC,3,MSYM)
      INTEGER IIBSTRSH(MSYM,MMSHL,3), NNBSTRSH(MSYM,MMSHL,3)
      INTEGER INDRED(MMBSTRT,3), INDRSH(MMBSTRT)
#include <cholesky.h>

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CHO_SETADDR')

C     Check dimensions.
C     -----------------

      NNSHL = NSHELL*(NSHELL + 1)/2
      IF (MRED .NE. MAXRED) THEN
         CALL CHO_QUIT('Error [1] in '//SECNAM,104)
      END IF
      IF (MVEC .NE. MAXVEC) THEN
         CALL CHO_QUIT('Error [2] in '//SECNAM,104)
      END IF
      IF (MSYM .NE. NSYM) THEN
         CALL CHO_QUIT('Error [3] in '//SECNAM,104)
      END IF
      IF (MMSHL .NE. NNSHL) THEN
         CALL CHO_QUIT('Error [4] in '//SECNAM,104)
      END IF
      IF (MMBSTRT .NE. NNBSTRT(1)) THEN
         CALL CHO_QUIT('Error [5] in '//SECNAM,104)
      END IF

C     Set addresses.
C     --------------

      IF (XNPASS .EQ. 0) THEN
         INFRED(1) = 1
         DO ISYM = 1,NSYM
            INFVEC(1,3,ISYM) = 1
         END DO
      ELSE IF (XNPASS .GT. 0) THEN
         IRED  = 3
         IPASS = XNPASS
         CALL CHO_GETRED(INFRED,NNBSTRSH(1,1,IRED),INDRED(1,IRED),
     &                   INDRSH,MRED,MSYM,MMSHL,MMBSTRT,IPASS,.FALSE.)
         CALL CHO_SETREDIND(IIBSTRSH,NNBSTRSH,MSYM,MMSHL,IRED)
         INFRED(IPASS+1) = INFRED(IPASS)
     &                   + NSYM*NNSHL + NNBSTRT(IRED)
         DO ISYM = 1,NSYM
            IF (NUMCHO(ISYM) .EQ. 0) THEN
               INFVEC(1,3,ISYM) = 1
            ELSE IF (NUMCHO(ISYM) .GT. 0) THEN
               INFVEC(NUMCHO(ISYM)+1,3,ISYM) =
     &         INFVEC(NUMCHO(ISYM),3,ISYM) + NNBSTR(ISYM,IRED)
            ELSE
               CALL CHO_QUIT('Dec. restart error in '//SECNAM,104)
            END IF
         END DO
      ELSE
         CALL CHO_QUIT('XNPASS is negative in '//SECNAM,104)
      END IF

      END
C  /* Deck cho_prttim */
      SUBROUTINE CHO_PRTTIM(SECTION,TCPU2,TCPU1,TWALL2,TWALL1,IOPT)
C
C     Purpose: print timing for a section.
C
      IMPLICIT NONE
      CHARACTER*(*) SECTION
      REAL*8        TCPU1, TCPU2, TWALL1, TWALL2
      INTEGER       IOPT
#include <cholesky.h>

      INTEGER IHRC, IMNC
      INTEGER IHRW, IMNW
      INTEGER LENSEC
      REAL*8  TCPUT, TWALLT, SECC, SECW
      CHARACTER*80 STRING

      TCPUT  = TCPU2  - TCPU1
      TWALLT = TWALL2 - TWALL1
      CALL CHO_CNVTIM(TCPUT,IHRC,IMNC,SECC)
      CALL CHO_CNVTIM(TWALLT,IHRW,IMNW,SECW)

      IF (IOPT .EQ. 1) THEN
         LENSEC = LEN(SECTION)
         WRITE(LUPRI,'(///,A,A,A)')
     &   '***** ',SECTION(1:LENSEC),' completed *****'
         WRITE(LUPRI,'(A,I8,A,I2,A,F6.2,A)')
     &   'Total CPU  time:',IHRC,' hours ',IMNC,' minutes ',SECC,
     &   ' seconds'
         WRITE(LUPRI,'(A,I8,A,I2,A,F6.2,A,//)')
     &   'Total wall time:',IHRW,' hours ',IMNW,' minutes ',SECW,
     &   ' seconds'
      ELSE IF (IOPT .EQ. 2) THEN
         LENSEC = MIN(LEN(SECTION),70)
         WRITE(STRING,'(A10,A)') 'Timing of ',SECTION(1:LENSEC)
         LENSEC = LENSEC + 10
         CALL CHO_HEAD(STRING(1:LENSEC),'=',80,LUPRI)
         WRITE(LUPRI,'(/,A,I8,A,I2,A,I2,A)')
     &   'Total CPU  time:',IHRC,' hours ',IMNC,' minutes ',SECC,
     &   ' seconds'
         WRITE(LUPRI,'(A,I8,A,I2,A,I2,A)')
     &   'Total wall time:',IHRW,' hours ',IMNW,' minutes ',SECW,
     &   ' seconds'
      ELSE
         WRITE(LUPRI,'(/,A,I8,A,I2,A,I2,A)')
     &   'Total CPU  time:',IHRC,' hours ',IMNC,' minutes ',SECC,
     &   ' seconds'
         WRITE(LUPRI,'(A,I8,A,I2,A,I2,A)')
     &   'Total wall time:',IHRW,' hours ',IMNW,' minutes ',SECW,
     &   ' seconds'
      END IF

      CALL CHO_FLUSH(LUPRI)

      END
C  /* Deck cho_cnvtim */
      SUBROUTINE CHO_CNVTIM(TIM,IHR,IMN,SEC)
C
C     Purpose: convert TIM to hours/minutes/seconds
C
      IMPLICIT NONE
      REAL*8  TIM
      INTEGER IHR, IMN
      REAL*8  XHR, XMN, SEC

      XHR = TIM/3.6D3
      IHR = INT(XHR)

      XMN = (TIM - 3.6D3*DBLE(IHR))/6.0D1
      IMN = INT(XMN)

      SEC = TIM - 3.6D3*DBLE(IHR) - 6.0D1*DBLE(IMN)

      END
C  /* Deck cho_prtred */
      SUBROUTINE CHO_PRTRED(NNBSTRSH,MSYM,MMSHL,IOPT)
C
C     Purpose: print information about reduced set.
C
#include <implicit.h>
      INTEGER NNBSTRSH(MSYM,MMSHL,3)
#include <choorb.h>
#include <cholesky.h>

      DIMENSION XBAS(8), XXBAS(8)

      INTEGER NSHP(2)
      LOGICAL CONTRIB(2)

      NNSHL = NSHELL*(NSHELL + 1)/2

      DO ISYM = 1,NSYM
         XBAS(ISYM) = DBLE(NBAS(ISYM))
      END DO

      XXBAST = 0.0D0
      DO ISYM = 1,NSYM
         XXBAS(ISYM) = 0.0D0
         DO ISYMB = 1,NSYM
            ISYMA = MULD2H(ISYMB,ISYM)
            IF (ISYMA .EQ. ISYMB) THEN
               XXBAS(ISYM) = XXBAS(ISYM)
     &                     + XBAS(ISYMA)*(XBAS(ISYMA) + 1.0D0)/2.0D0
            ELSE IF (ISYMA .GT. ISYMB) THEN
               XXBAS(ISYM) = XXBAS(ISYM) + XBAS(ISYMA)*XBAS(ISYMB)
            END IF
         END DO
         XXBAST = XXBAST + XXBAS(ISYM)
      END DO

      IF (IOPT .EQ. 1) THEN
         NRED = 1
      ELSE
         NRED = 2
      END IF

      DO IRED = 1,NRED
         NSHP(IRED) = 0
         DO ISHLAB = 1,NNSHL
            CONTRIB(IRED) = .FALSE.
            ISYM = 0
            DO WHILE ((.NOT.CONTRIB(IRED)) .AND. ISYM.LT.NSYM)
               ISYM = ISYM + 1
               IF (NNBSTRSH(ISYM,ISHLAB,IRED) .GT. 0) THEN
                  CONTRIB(IRED) = .TRUE.
               END IF
            END DO
            IF (CONTRIB(IRED)) THEN
               NSHP(IRED) = NSHP(IRED) + 1
            END IF
         END DO
      END DO

      CALL CHO_HEAD('Reduced Set Information','=',80,LUPRI)

      IF (NNSHL .EQ. 0) THEN
         PCT1 = 9.9D9
      ELSE
         PCT1 = 1.0D2*DBLE(NSHP(1))/DBLE(NNSHL)
      END IF
      IF (IOPT .EQ. 1) THEN
         WRITE(LUPRI,'(/,A,/,A)')
     &   'Sym.          Full   First Red. Set',
     &   '-----------------------------------'
         DO ISYM = 1,NSYM
            WRITE(LUPRI,'(I3,3X,F12.1,7X,I10)')
     &      ISYM,XXBAS(ISYM),NNBSTR(ISYM,1)
         END DO
         WRITE(LUPRI,'(A)')
     &   '-----------------------------------'
         WRITE(LUPRI,'(A,F12.1,7X,I10)')
     &   'Total:',XXBAST,NNBSTRT(1)
         WRITE(LUPRI,'(A)')
     &   '-----------------------------------'
         WRITE(LUPRI,'(/,A,I10,A,I10,A,F7.2,A)')
     &   'First Reduced Set:',NSHP(1),' of',NNSHL,
     &   ' shell pairs contribute (',PCT1,'%)'
      ELSE
         IF (NNSHL .EQ. 0) THEN
            PCT2 = 9.9D9
         ELSE
            PCT2 = 1.0D2*DBLE(NSHP(2))/DBLE(NNSHL)
         END IF
         WRITE(LUPRI,'(/,A,/,A,/,A)')
     &   '                          Reduced Set',
     &   'Sym.          Full      First    Current',
     &   '----------------------------------------'
         DO ISYM = 1,NSYM
            WRITE(LUPRI,'(I3,3X,F12.1,1X,I10,1X,I10)')
     &      ISYM,XXBAS(ISYM),NNBSTR(ISYM,1),NNBSTR(ISYM,2)
         END DO
         WRITE(LUPRI,'(A)')
     &   '----------------------------------------'
         WRITE(LUPRI,'(A,F12.1,1X,I10,1X,I10)')
     &   'Total:',XXBAST,NNBSTRT(1),NNBSTRT(2)
         WRITE(LUPRI,'(A)')
     &   '----------------------------------------'
         WRITE(LUPRI,'(/,A,I10,A,I10,A,F7.2,A)')
     &   'First Reduced Set:',NSHP(1),' of',NNSHL,
     &   ' shell pairs contribute (',PCT1,'%)'
         WRITE(LUPRI,'(A,I10,A,I10,A,F7.2,A)')
     &   'Curr. Reduced Set:',NSHP(2),' of',NNSHL,
     &   ' shell pairs contribute (',PCT2,'%)'
      END IF

      CALL CHO_FLUSH(LUPRI)

      END
C  /* Deck cho_intchk_reg */
      SUBROUTINE CHO_INTCHK_REG(LABEL,ISHLCD,ISHLAB)
C
C     Purpose: register a shell quadruple (CD|AB) for minimal integral
C              check using LABEL to keep track of its origin.
C
#include <implicit.h>
      CHARACTER*8 LABEL
#include <cholesky.h>

      CHARACTER*14 SECNAM
      PARAMETER (SECNAM = 'CHO_INTCHK_REG')

C     Check shell pair index.
C     -----------------------

      NNSHL = NSHELL*(NSHELL + 1)/2
      IF (ISHLCD.LT.1 .OR. ISHLCD.GT.NNSHL) THEN
         CALL CHO_QUIT('Shell index error 1 in '//SECNAM,103)
      END IF
      IF (ISHLAB.LT.1 .OR. ISHLAB.GT.NNSHL) THEN
         CALL CHO_QUIT('Shell index error 2 in '//SECNAM,103)
      END IF

C     Registration.
C     -------------

      CALL CHO_INVPCK(ISHLCD,ISHLC,ISHLD,.TRUE.)
      CALL CHO_INVPCK(ISHLAB,ISHLA,ISHLB,.TRUE.)
      CALL CHO_INTCHK_ID_OF(LABEL,ID,1)
      IF (ID.LT.1 .OR. ID.GT.NCHKQ) THEN
         ID = NCHKQ + 1 ! junk yard
      END IF
      ICHKQ(1,ID) = ISHLC
      ICHKQ(2,ID) = ISHLD
      ICHKQ(3,ID) = ISHLA
      ICHKQ(4,ID) = ISHLB

      END
C  /* Deck cho_intchk_id_of */
      SUBROUTINE CHO_INTCHK_ID_OF(LABEL,ID,IOPT)
C
C     Purpose: for minimal integral checking,
C              IOPT=-1 : return label corresponding to id ID.
C              else    : return index of shell quadruple corresponding
C                        to check label LABEL.
C
      IMPLICIT NONE
      CHARACTER*8 LABEL
      INTEGER ID, IOPT
#include <cholesky.h>

      INTEGER CHO_TABIND

      INTEGER     NTABLE
      PARAMETER   (NTABLE = 7)
      CHARACTER*8 TABLE(NTABLE)
      DATA TABLE  /'EXCL RS1','MAX|XRS1',
     &             'NEG DIAG','MAX|NEG ',
     &             'NEG->ZER','MAX|NEGZ',
     &             'MAX ERR '/

      IF (IOPT .EQ. -1) THEN
         IF (ID.LT.1 .OR. ID.GT.NTABLE) THEN
            LABEL = 'UNKNOWN '
         ELSE
            LABEL = TABLE(ID)
         END IF
      ELSE
         ID = CHO_TABIND(TABLE,8,NTABLE,' ',0,0,LABEL)
      END IF

      END
C  /* Deck cho_getvec */
      SUBROUTINE CHO_GETVEC(CHOVEC,LENVEC,NUMVEC,IVEC1,ISYM,
     &                      INFRED,INFVEC,
     &                      INDRED,INDRSH,IIBSTRSH,NNBSTRSH,
     &                      ISCR,
     &                      MRED,MVEC,MSYM,MMBSTRT,MMSHL,LISCR,
     &                      SCR,LSCR)
C
C     Purpose: read Cholesky vectors IVEC=IVEC1,....,IVEC1+NUMVEC-1
C              of symmetry ISYM from file. The vectors are returned
C              in the "current" reduced set. The algorithm used for
C              reading is taken from input (via cholesky.h header
C              file).
C
#include <implicit.h>
      DIMENSION CHOVEC(LENVEC,NUMVEC)
      INTEGER   INFRED(MRED), INFVEC(MVEC,3,MSYM)
      INTEGER   INDRED(MMBSTRT,3), INDRSH(MMBSTRT)
      INTEGER   IIBSTRSH(MSYM,MMSHL,3), NNBSTRSH(MSYM,MMSHL,3)
      INTEGER   ISCR(LISCR)
      DIMENSION SCR(LSCR)
#include <cholesky.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_GETVEC')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

C     Return if no vectors requested.
C     -------------------------------

      IF (NUMVEC .LT. 1) THEN
         IF (LOCDBG) THEN
            WRITE(LUPRI,*) SECNAM,': WARNING: no vectors in this call!'
            WRITE(LUPRI,*) SECNAM,': NUMVEC = ',NUMVEC
         END IF
         RETURN
      END IF

C     Check vector dimension: should be identical to current reduced
C     set. Check also symmetry and vector index.
C     --------------------------------------------------------------

      IF (LOCDBG) THEN
         IFAIL = 0
         IF (LENVEC.NE.NNBSTR(ISYM,2) .OR. LENVEC.LT.1) THEN
            WRITE(LUPRI,*) SECNAM,': illegal vector dimension:'
            WRITE(LUPRI,*) SECNAM,': LENVEC = ',LENVEC
            IFAIL = IFAIL + 1
         END IF
         IF (ISYM.LT.1 .OR. ISYM.GT.NSYM) THEN
            WRITE(LUPRI,*) SECNAM,': illegal symmetry input:'
            WRITE(LUPRI,*) SECNAM,': ISYM = ',ISYM
            IFAIL = IFAIL + 1
         END IF
         IVEC2 = IVEC1 + NUMVEC - 1
         IF (IVEC1.LT.1 .OR. IVEC1.GT.MAXVEC .OR.
     &       IVEC2.LT.1 .OR. IVEC2.GT.MAXVEC) THEN
            WRITE(LUPRI,*) SECNAM,': illegal vector index:'
            WRITE(LUPRI,*) SECNAM,': IVEC1,IVEC2 = ',IVEC1,IVEC2
            IFAIL = IFAIL + 1
         ELSE
            IF (INFVEC(IVEC1,3,ISYM) .LT. 1) THEN
               WRITE(LUPRI,*) SECNAM,': illegal first vector address:'
               WRITE(LUPRI,*) SECNAM,': address: ',INFVEC(IVEC1,3,ISYM)
               IFAIL = IFAIL + 1
            END IF
            IF (INFVEC(IVEC2,3,ISYM) .LT. 1) THEN
               WRITE(LUPRI,*) SECNAM,': illegal last vector address:'
               WRITE(LUPRI,*) SECNAM,': address: ',INFVEC(IVEC2,3,ISYM)
               IFAIL = IFAIL + 1
            END IF
         END IF
         IF (CHO_IOVEC .EQ. 1) THEN
            IF (LISCR .LT. NNBSTR(ISYM,2)) THEN
               WRITE(LUPRI,*) SECNAM,': insufficient iscratch:'
               WRITE(LUPRI,*) SECNAM,': LISCR = ',LISCR
               IFAIL = IFAIL + 1
            END IF
         END IF
         IF (IFAIL .NE. 0) THEN
            WRITE(LUPRI,*) SECNAM,': unable to continue!'
            CALL CHO_QUIT('Error in '//SECNAM,104)
         END IF
      END IF

C     Call reading routine.
C     ---------------------

      IF (CHO_IOVEC .EQ. 1) THEN
         IREDC = -1
         CALL CHO_GETVEC1(CHOVEC,LENVEC,NUMVEC,IVEC1,ISYM,
     &                    INFRED,INFVEC,
     &                    INDRED,INDRSH,IIBSTRSH,NNBSTRSH,
     &                    ISCR,
     &                    MRED,MVEC,MSYM,MMBSTRT,MMSHL,LISCR,
     &                    SCR,LSCR,IREDC)
      ELSE
         CALL CHO_GETVEC2(CHOVEC,LENVEC,NUMVEC,IVEC1,ISYM,
     &                    INFRED,INFVEC,
     &                    INDRED,INDRSH,IIBSTRSH,NNBSTRSH,
     &                    MRED,MVEC,MSYM,MMBSTRT,MMSHL,
     &                    SCR,LSCR)
      END IF

C     Debug print.
C     ------------

      IF (LOCDBG) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) SECNAM,':'
         WRITE(LUPRI,*) 'Vectors ',IVEC1,' to ',IVEC1+NUMVEC-1,
     &                  ' of symmetry ',ISYM,' read from unit ',
     &                  LUNCHO(ISYM)
         WRITE(LUPRI,*) 'Vector dimension: ',LENVEC,
     &                  ' (current reduced set)'
         WRITE(LUPRI,*) 'Algorithm: ',CHO_IOVEC
         DO IVEC = 1,NUMVEC
            JVEC = IVEC1 + IVEC - 1
            JADR = INFVEC(JVEC,3,ISYM)
            XNRM = SQRT(DDOT(LENVEC,CHOVEC(1,IVEC),1,CHOVEC(1,IVEC),1))
            WRITE(LUPRI,*) 'Vector:',JVEC,' address: ',JADR,' norm: ',
     &                     XNRM
         END DO
      END IF

      END
C  /* Deck cho_getvec1 */
      SUBROUTINE CHO_GETVEC1(CHOVEC,LENVEC,NUMVEC,IVEC1,ISYM,
     &                       INFRED,INFVEC,
     &                       INDRED,INDRSH,IIBSTRSH,NNBSTRSH,
     &                       ISCR,
     &                       MRED,MVEC,MSYM,MMBSTRT,MMSHL,LISCR,
     &                       SCR,LSCR,IREDC)
C
C     Purpose: read Cholesky vectors IVEC=IVEC1,....,IVEC1+NUMVEC-1
C              of symmetry ISYM from file. The vectors are returned
C              in the "current" reduced set. This routine attempts
C              to minimize gather/scatter operations and uses batched
C              reading to (hopefully) improve buffering.
C              IREDC identifies the reduced set stored at location 3 in
C              the index arrays (NNBSTR(*,3) etc.). If unknown or not
C              set, specify IREDC=-1. Note that IREDC is input and
C              output!
C
C     NOTE: the scratch array SCR(LSCR) is used to read vectors from
C           disk and should not be smaller than NNBSTR(ISYM,1)+1.
C
#include <implicit.h>
      DIMENSION CHOVEC(LENVEC,NUMVEC)
      INTEGER   INFRED(MRED), INFVEC(MVEC,3,MSYM)
      INTEGER   INDRED(MMBSTRT,3), INDRSH(MMBSTRT)
      INTEGER   IIBSTRSH(MSYM,MMSHL,3), NNBSTRSH(MSYM,MMSHL,3)
      INTEGER   ISCR(LISCR)
      DIMENSION SCR(LSCR)
#include <cholesky.h>

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CHO_GETVEC1')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER IOFF(0:1)

      INTEGER*8 IADR

C     Some initializations.
C     ---------------------

      NNSHL = NSHELL*(NSHELL + 1)/2
      IVEC2 = IVEC1 + NUMVEC - 1

      IF (LISCR .LT. NNBSTR(ISYM,2)) THEN
         CALL CHO_QUIT('LISCR too small in '//SECNAM,103)
      END IF

      KJUNK = 1
      KSCR  = KJUNK + 1
      LEFT  = LSCR  - KSCR + 1
      IF (LEFT .LE. 0) THEN
         CALL CHO_QUIT('Insufficient scratch space in '//SECNAM,101)
      END IF

      SCR(KJUNK) = 0.0D0
      IOFF(0)    = KJUNK

C     Get reduced sets of first and last vector.
C     ------------------------------------------

      IRED1 = INFVEC(IVEC1,2,ISYM)
      IRED2 = INFVEC(IVEC2,2,ISYM)

C     Loop through reduced sets to be read.
C     -------------------------------------

      KVEC1 = 1
      JVEC1 = IVEC1
      DO IRED = IRED1,IRED2

C        Count vectors in this reduced set.
C        ----------------------------------

         JNUM = 0
         JVEC = JVEC1 - 1
         DO WHILE (JVEC.LT.IVEC2)
            JVEC = JVEC + 1
            JRED = INFVEC(JVEC,2,ISYM)
            IF (JRED .EQ. IRED) THEN
               JNUM = JNUM + 1 ! increase counter
            ELSE
               JVEC = IVEC2 ! break while loop
            END IF
         END DO

C        Skip if this reduced set is empty.
C        ----------------------------------

         IF (JNUM .EQ. 0) GO TO 100

C        Check vector range.
C        -------------------

         IF (LOCDBG) THEN
            JVEC2 = JVEC1 + JNUM - 1
            IF (JVEC2 .GT. IVEC2) THEN
               WRITE(LUPRI,*) SECNAM,': IRED  = ',IRED
               WRITE(LUPRI,*) SECNAM,': JNUM  = ',JNUM
               WRITE(LUPRI,*) SECNAM,': JVEC1 = ',JVEC1
               WRITE(LUPRI,*) SECNAM,': JVEC2 = ',JVEC2
               CALL CHO_QUIT('Vector index error in '//SECNAM,103)
            END IF
         END IF

C        Read reduced set index arrays.
C        ------------------------------

         IF (IRED .NE. IREDC) THEN
            CALL CHO_GETRED(INFRED,NNBSTRSH(1,1,3),INDRED(1,3),INDRSH,
     &                      MRED,MSYM,MMSHL,MMBSTRT,IRED,.FALSE.)
            CALL CHO_SETREDIND(IIBSTRSH,NNBSTRSH,MSYM,MMSHL,3)
            IREDC = IRED
         END IF

C        If reduced sets are identical, simply read the vectors
C        directly into CHOVEC array and go to next reduced set.
C        ------------------------------------------------------

         IF (NNBSTR(ISYM,3) .EQ. NNBSTR(ISYM,2)) THEN
            IADR = INFVEC(JVEC1,3,ISYM)
            LTOT = NNBSTR(ISYM,2)*JNUM
            CALL GETWA2U(LUNCHO(ISYM),FNVEC(ISYM),CHOVEC(1,KVEC1),IADR,
     &                   LTOT)
            GO TO 100
         END IF

C        For the first reduced set, INDRED addresses into shell pair.
C        Therefore, we need to reset it for our present purposes.
C        (Only the relevant symmetry block is reset.)
C        ------------------------------------------------------------

         IF (IRED .EQ. 1) THEN
            IF (NNBSTRT(1) .NE. NNBSTRT(3)) THEN
               WRITE(LUPRI,*) 'Dimension problem: rs1 != rs3 for IRED=1'
               CALL CHO_QUIT('Dimension error in '//SECNAM,103)
            END IF
            I1 = IIBSTR(ISYM,3) + 1
            I2 = I1 + NNBSTR(ISYM,3) - 1
            DO I = I1,I2
               INDRED(I,3) = I
            END DO
         END IF

C        Set up batch over vectors to be read.
C        -------------------------------------

         MINL = NNBSTR(ISYM,3)
         IF (MINL .LT. 1) THEN
            NVEC = 0
         ELSE
            NVEC = MIN(LEFT/MINL,JNUM)
         END IF
         IF (NVEC .LT. 1) THEN
            WRITE(LUPRI,*) SECNAM,': insufficient scratch space:'
            WRITE(LUPRI,*) 'LEFT = ',LEFT
            WRITE(LUPRI,*) 'JNUM = ',JNUM
            WRITE(LUPRI,*) 'MINL = ',MINL
            WRITE(LUPRI,*) 'NVEC = ',NVEC
            WRITE(LUPRI,*) 'Input:'
            WRITE(LUPRI,*) 'IVEC1  = ',IVEC1
            WRITE(LUPRI,*) 'NUMVEC = ',NUMVEC
            WRITE(LUPRI,*) 'LENVEC = ',LENVEC
            WRITE(LUPRI,*) 'ISYM   = ',ISYM
            CALL CHO_QUIT('Insufficient scratch space in '//SECNAM,104)
            NBATCH = 0 ! to avoid compiler warnings
         ELSE
            NBATCH = (JNUM - 1)/NVEC + 1
         END IF

C        Set up mapping between reduced sets.
C        NOTE: this requires that the reduced sets are created
C        by simply popping elements of the reduced lists so that
C        the basic ordering is conserved between the sets.
C        -------------------------------------------------------

         CALL CHO_IZERO(ISCR,NNBSTR(ISYM,2))
         DO ISHLAB = 1,NNSHL
            N2 = NNBSTRSH(ISYM,ISHLAB,2)
            N3 = NNBSTRSH(ISYM,ISHLAB,3)
            IF (N2.GT.0 .AND. N3.GT.0) THEN
               IF (N2 .LT. N3) THEN
                  IAB1 = IIBSTRSH(ISYM,ISHLAB,2) + 1
                  IAB2 = IAB1 + N2 - 1
                  LAST = 0
                  DO IAB = IAB1,IAB2
                     JAB = INDRED(IIBSTR(ISYM,2)+IAB,2)
                     KAB = LAST
                     DO WHILE (KAB .LT. NNBSTRSH(ISYM,ISHLAB,3))
                        KAB = KAB + 1
                        LAB = IIBSTRSH(ISYM,ISHLAB,3) + KAB
                        MAB = INDRED(IIBSTR(ISYM,3)+LAB,3)
                        IF (MAB .EQ. JAB) THEN
                           ISCR(IAB) = LAB
                           LAST = KAB
                           KAB  = NNBSTRSH(ISYM,ISHLAB,3)
                        END IF
                     END DO
                  END DO
               ELSE
                  LAB1 = IIBSTRSH(ISYM,ISHLAB,3) + 1
                  LAB2 = LAB1 + N3 - 1
                  LAST = 0
                  DO LAB = LAB1,LAB2
                     MAB = INDRED(IIBSTR(ISYM,3)+LAB,3)
                     KAB = LAST
                     DO WHILE (KAB .LT. NNBSTRSH(ISYM,ISHLAB,2))
                        KAB = KAB + 1
                        IAB = IIBSTRSH(ISYM,ISHLAB,2) + KAB
                        JAB = INDRED(IIBSTR(ISYM,2)+IAB,2)
                        IF (JAB .EQ. MAB) THEN
                           ISCR(IAB) = LAB
                           LAST = KAB
                           KAB  = NNBSTRSH(ISYM,ISHLAB,2)
                        END IF
                     END DO
                  END DO
               END IF
            END IF
         END DO

C        Start batch loop.
C        -----------------

         DO IBATCH = 1,NBATCH

            IF (IBATCH .EQ. NBATCH) THEN
               NUMV = JNUM - NVEC*(NBATCH - 1)
            ELSE
               NUMV = NVEC
            END IF
            IBVEC1 = JVEC1 + NVEC*(IBATCH - 1)
            KBVEC1 = KVEC1 + NVEC*(IBATCH - 1)

C           Read vectors.
C           -------------

            IADR = INFVEC(IBVEC1,3,ISYM)
            LTOT = NNBSTR(ISYM,3)*NUMV
            CALL GETWA2U(LUNCHO(ISYM),FNVEC(ISYM),SCR(KSCR),IADR,LTOT)

C           Debug print.
C           ------------

            IF (LOCDBG) THEN
               WRITE(LUPRI,*)
               WRITE(LUPRI,*) SECNAM,':'
               WRITE(LUPRI,*) 'Vectors ',IBVEC1,' to ',IBVEC1+NUMV-1,
     &                        ' of symmetry ',ISYM,' read from unit ',
     &                        LUNCHO(ISYM)
               WRITE(LUPRI,*) 'Vector dimension: ',NNBSTR(ISYM,3),
     &                        ' (own reduced set)'
               DO IIII = 1,NUMV
                  JJJJ = IBVEC1 + IIII - 1
                  JADR = INFVEC(JJJJ,3,ISYM)
                  KKKK = KSCR + NNBSTR(ISYM,3)*(IIII-1)
                  XNRM = SQRT(DDOT(NNBSTR(ISYM,3),SCR(KKKK),1,
     &                                            SCR(KKKK),1))
                  WRITE(LUPRI,*)
     &            'Vector:',JJJJ,' address: ',JADR,' norm: ',
     &            XNRM
               END DO
               CALL CHO_FLUSH(LUPRI)
            END IF

C           Copy vectors into result array.
C           -------------------------------

            DO JVEC = 1,NUMV
               KVEC = KBVEC1 + JVEC - 1
               IOFF(1) = IOFF(0) + NNBSTR(ISYM,3)*(JVEC - 1)
               DO IAB = 1,NNBSTR(ISYM,2)
                  KOFF = IOFF(MIN(ISCR(IAB),1)) + ISCR(IAB)
                  CHOVEC(IAB,KVEC) = SCR(KOFF)
               END DO
            END DO

         END DO

C        Set next vector to be treated.
C        ------------------------------

  100    KVEC1 = KVEC1 + JNUM
         JVEC1 = JVEC1 + JNUM

      END DO

      END
C  /* Deck cho_getvec2 */
      SUBROUTINE CHO_GETVEC2(CHOVEC,LENVEC,NUMVEC,IVEC1,ISYM,
     &                       INFRED,INFVEC,
     &                       INDRED,INDRSH,IIBSTRSH,NNBSTRSH,
     &                       MRED,MVEC,MSYM,MMBSTRT,MMSHL,
     &                       SCR,LSCR)
C
C     Purpose: read Cholesky vectors IVEC=IVEC1,...,IVEC1+NUMVEC-1
C              of symmetry ISYM from file. The vectors are returned
C              in the "current" reduced set.
C
C     NOTE: array SCR(LSCR) is used for storing the vectors in the
C           red. set from disk and for a full first red. set vector.
C           Thus, to be certain that enough memory is available,
C           use LSCR = 2 x dimension of first reduced set.
C
C     TODO/FIXME: read more than 1 vector at a time...use batch based on
C                 LSCR similar to getvec1 routine.
C
#include <implicit.h>
      DIMENSION CHOVEC(LENVEC,NUMVEC)
      INTEGER   INFRED(MRED), INFVEC(MVEC,3,MSYM)
      INTEGER   INDRED(MMBSTRT,3), INDRSH(MMBSTRT)
      INTEGER   IIBSTRSH(MSYM,MMSHL,3), NNBSTRSH(MSYM,MMSHL,3)
      DIMENSION SCR(LSCR)
#include <cholesky.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_GETVEC')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER*8 IADR

C     Initialize output array.
C     ------------------------

      CALL CHO_DZERO(CHOVEC,LENVEC*NUMVEC)

C     Read reduced set index arrays for first vector.
C     -----------------------------------------------

      IRED = INFVEC(IVEC1,2,ISYM)
      CALL CHO_GETRED(INFRED,NNBSTRSH(1,1,3),INDRED(1,3),INDRSH,
     &                MRED,MSYM,MMSHL,MMBSTRT,IRED,.FALSE.)
      CALL CHO_SETREDIND(IIBSTRSH,NNBSTRSH,MSYM,MMSHL,3)
      KRED1 = 1
      KREAD = KRED1 + NNBSTR(ISYM,1)
      KEND1 = KREAD + NNBSTR(ISYM,3)
      LSCR1 = LSCR  - KEND1 + 1
      IF (LSCR1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Insufficient scratch space in ',SECNAM
         WRITE(LUPRI,*) 'Available: ',LSCR,'   Need: ',KEND1-1
         WRITE(LUPRI,*) '- needed for RED1: ',NNBSTR(ISYM,1)
         WRITE(LUPRI,*) '- needed for READ: ',NNBSTR(ISYM,3)
         CALL CHO_QUIT('[1] Insufficient scratch space in '//SECNAM,102)
      END IF

C     Read vectors and re-order into current reduced set via reduced
C     set 1.
C     NOTE: if the read vectors are already in red. set 1, don't resort.
C     ------------------------------------------------------------------

      DO JVEC = 1,NUMVEC

         IVEC = IVEC1 + JVEC - 1
         JRED = INFVEC(IVEC,2,ISYM)
         IF (JRED .NE. IRED) THEN   ! read new reduced set
            CALL CHO_GETRED(INFRED,NNBSTRSH(1,1,3),INDRED(1,3),INDRSH,
     &                      MRED,MSYM,MMSHL,MMBSTRT,JRED,.FALSE.)
            CALL CHO_SETREDIND(IIBSTRSH,NNBSTRSH,MSYM,MMSHL,3)
            KEND1 = KREAD + NNBSTR(ISYM,3)
            LSCR1 = LSCR  - KEND1 + 1
            IF (LSCR1 .LT. 0) THEN
               WRITE(LUPRI,*) 'Insufficient scratch space in ',SECNAM
               WRITE(LUPRI,*) 'Available: ',LSCR,'   Need: ',KEND1-1
               WRITE(LUPRI,*) '- needed for RED1: ',NNBSTR(ISYM,1)
               WRITE(LUPRI,*) '- needed for READ: ',NNBSTR(ISYM,3)
               CALL CHO_QUIT('[2] Insufficient scratch space in '
     &                       //SECNAM,102)
            END IF
            IRED = JRED
         END IF

         IADR = INFVEC(IVEC,3,ISYM)
         CALL GETWA2U(LUNCHO(ISYM),FNVEC(ISYM),SCR(KREAD),IADR,
     &                NNBSTR(ISYM,3))

         IF (LOCDBG) THEN
            IADR = INFVEC(IVEC,3,ISYM)
            XNRM = SQRT(DDOT(NNBSTR(ISYM,3),SCR(KREAD),1,SCR(KREAD),1))
            WRITE(LUPRI,*) SECNAM,': ',
     &                     'Vector:',IVEC,' address: ',
     &                     INFVEC(IVEC,3,ISYM),
     &                     ' norm: ',XNRM,' sym. ',ISYM,' red. set: ',
     &                     IRED,' dim.: ',NNBSTR(ISYM,3)
         END IF

         CALL CHO_DZERO(SCR(KRED1),NNBSTR(ISYM,1))
         IF (IRED .GT. 1) THEN
            DO JAB = 1,NNBSTR(ISYM,3)   ! sort into rs1 ordering
               KAB = IIBSTR(ISYM,3)  + JAB
               IAB = INDRED(KAB,3) - IIBSTR(ISYM,1)
               SCR(KRED1+IAB-1) = SCR(KREAD+JAB-1)
            END DO
            KREDU = KRED1   ! point rs2 sort to red1 resort
         ELSE IF (IRED .EQ. 1) THEN
            KREDU = KREAD   ! point rs2 sort to read (already red1)
         ELSE
            WRITE(LUPRI,*) SECNAM,': ERROR: IRED is negative: ',IRED
            CALL CHO_QUIT('Reduced set error in '//SECNAM,104)
            KREDU = -999999 ! just to avoid compiler warnings
         END IF

         DO JAB = 1,NNBSTR(ISYM,2)  ! sort into in rs2 ordering
            KAB = IIBSTR(ISYM,2)  + JAB
            IAB = INDRED(KAB,2) - IIBSTR(ISYM,1)
            CHOVEC(JAB,JVEC) = SCR(KREDU+IAB-1)
         END DO

      END DO

      END
C  /* Deck cho_lread */
      INTEGER FUNCTION CHO_LREAD(ISYM,INFVEC,MVEC,MSYM,LWORK)
C
C     Purpose: return a reasonable scratch space dimension for reading
C              previous vectors using cho_getvec.
C
      IMPLICIT NONE
      INTEGER ISYM, MVEC, MSYM, LWORK
      INTEGER INFVEC(MVEC,3,MSYM)
#include <cholesky.h>

      INTEGER JVEC, IRED, JRED, LEN1, LEN2, LEN3, LMIN

      INTEGER MNVECRS1
      PARAMETER (MNVECRS1 = 5)

      LMIN = 2*NNBSTR(ISYM,1)

      IF (CHO_IOVEC .EQ. 1) THEN
         IF (NVECRS1(ISYM).LT.1 .AND. NUMCHO(ISYM).GT.0) THEN
            NVECRS1(ISYM) = 1
            JVEC = 1
            IRED = INFVEC(JVEC,2,ISYM)
            DO WHILE (JVEC .LT. NUMCHO(ISYM))
               JVEC = JVEC + 1
               JRED = INFVEC(JVEC,2,ISYM)
               IF (JRED .EQ. IRED) THEN
                  NVECRS1(ISYM) = NVECRS1(ISYM) + 1
               ELSE
                  JVEC = NUMCHO(ISYM)
               END IF
            END DO
         END IF
         LEN1 = LWORK/3 - 1
         LEN2 = MAX(NVECRS1(ISYM),MNVECRS1)*NNBSTR(ISYM,1)
         LEN3 = MIN(LEN1,LEN2)
         CHO_LREAD = MAX(LEN3,LMIN) + 1
      ELSE
         CHO_LREAD = LMIN
      END IF

      END
C  /* Deck cho_flush */
      SUBROUTINE CHO_FLUSH(LUNIT)
C
C     Purpose: flush unit LUNIT.
C
      IMPLICIT NONE
      INTEGER LUNIT

      CALL FLSHFO(LUNIT)

      END
C  /* Deck cho_prthead */
      SUBROUTINE CHO_PRTHEAD(SKIPH)
C
C     Purpose: print Cholesky header.
C
#include <implicit.h>
      LOGICAL SKIPH
#include <cholesky.h>

      PARAMETER (NALG = 3)

      CHARACTER*7  INTPRO
      CHARACTER*12 ALGORITHM(0:NALG)
      CHARACTER*22 USED(2)

      DATA ALGORITHM /'     unknown','    one-step','    two-step',
     &                '       naive'/
      DATA USED /'(screening turned off)','(screening turned on) '/

      IF (LUPRI .LT. 1) THEN
         STOP ' *** LUPRI undefined in Cholesky decomposition *** '
      END IF

      IF (.NOT. SKIPH) THEN
         WRITE(LUPRI,'(//,80A)') ('*', I=1,80)
         WRITE(LUPRI,'(A,78X,A)') ('*',I=1,2)
         WRITE(LUPRI,'(A,10X,A,10X,A)')
     & '*','Cholesky Decomposition of Two-Electron Repulsion Integrals',
     &   '*'
         WRITE(LUPRI,'(A,78X,A)') ('*',I=1,2)
         WRITE(LUPRI,'(80A,//)') ('*', I=1,80)

         IF (RSTDIA) THEN
            WRITE(LUPRI,'(/,A)')
     &      '***** Using Restart Integral Diagonal *****'
         END IF
         IF (RSTCHO) THEN
            IF (RSTDIA) THEN
               WRITE(LUPRI,'(A)')
     &         '***** Using Restart Cholesky Vectors  *****'
            ELSE
               WRITE(LUPRI,'(/,A)')
     &         '***** Using Restart Cholesky Vectors  *****'
            END IF
         END IF
      END IF

      IF (SCDIAG) THEN
         IUSE = 2
      ELSE
         IUSE = 1
      END IF
      IF (CHO_DECALG.LT.1 .OR. CHO_DECALG.GT.NALG) THEN
         IALG = 0
      ELSE
         IALG = CHO_DECALG
      END IF
      IF (.NOT. SKIPH) CALL CHO_HEAD('Configuration','=',80,LUPRI)
      WRITE(LUPRI,'(/,A,A)')
     & 'Decomposition algorithm                   : ',ALGORITHM(IALG)
      WRITE(LUPRI,'(A,1P,D12.4)')
     & 'Decomposition threshold                   : ',THRCOM
      WRITE(LUPRI,'(A,1P,D12.4)')
     & 'Initial diagonal screening                : ',THRDIAG
      WRITE(LUPRI,'(A,1P,D12.4,1X,A)')
     & 'First  screening damping                  : ',DAMP(1),USED(IUSE)
      WRITE(LUPRI,'(A,1P,D12.4,1X,A)')
     & 'Second screening damping                  : ',DAMP(2),USED(IUSE)
      WRITE(LUPRI,'(A,1P,D12.4)')
     & 'Span factor                               : ',SPAN
      WRITE(LUPRI,'(A,I12)')
     & 'Max. #Cholesky vectors per symmetry       : ',MAXVEC
      WRITE(LUPRI,'(A,I12)')
     & 'Max. #reduced sets (i.e., integral passes): ',MAXRED
      WRITE(LUPRI,'(A,I12)')
     & 'Min. #qualified required for decomposition: ',MINQUAL
      WRITE(LUPRI,'(A,I12)')
     & 'Max. #qualified per symmetry              : ',MAXQUAL
      IF (MXSHPR .EQ. 0) THEN
         WRITE(LUPRI,'(A)')
     &   'Max. #shell pair allowed per integral pass:      generic'
      ELSE
         WRITE(LUPRI,'(A,I12)')
     &   'Max. #shell pair allowed per integral pass: ',MXSHPR
      END IF
      IF (IALQUA .EQ. 1) THEN
         WRITE(LUPRI,'(A)')
     &   'Qualification algorithm                   :   sequential'
      ELSE
         WRITE(LUPRI,'(A)')
     &   'Qualification algorithm                   :      sorting'
      END IF
      IF (CHO_IOVEC .EQ. 1) THEN
         WRITE(LUPRI,'(A)')
     &   'Algorithm for Cholesky vector I/O         : rs-to-rs map'
      ELSE
         WRITE(LUPRI,'(A)')
     &   'Algorithm for Cholesky vector I/O         : copy via rs1'
      END IF

      END
C  /* Deck cho_quit */
      SUBROUTINE CHO_QUIT(STRING,IERR)
C
C     Purpose: echo message STRING and abort execution.
C
#include <implicit.h>
      CHARACTER*(*) STRING
#include <cholesky.h>

      WRITE(LUPRI,'(//,1X,A)') '***'
      IF (IERR .EQ. 0) THEN
         WRITE(LUPRI,'(1X,A)')
     &   '*** Execution stopped by Cholesky Program'
      ELSE
         WRITE(LUPRI,'(1X,A)')    '*** Error in Cholesky Program'
         WRITE(LUPRI,'(1X,A,I5)') '*** Code:',IERR
      END IF
      WRITE(LUPRI,'(1X,A,//)')  '***'
      CALL CHO_FLUSH(LUPRI)
      CALL QUIT(STRING)

      END
C  /* Deck cho_mem */
      SUBROUTINE CHO_MEM(TEXT,JOB,TYPE,KPOINT,LENGTH,KFREE,LFREE)
C
C     Purpose: memory management.
C
      IMPLICIT NONE
      CHARACTER*(*) TEXT,JOB,TYPE
      INTEGER KPOINT,LENGTH,KFREE,LFREE

      CALL CHO_DAL_MEM(TEXT,JOB,TYPE,KPOINT,LENGTH,KFREE,LFREE)

      END
C  /* Deck cho_output */
      SUBROUTINE CHO_OUTPUT(AMATRX,ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM,
     &                      COLDIM,NCTL,LUPRI)
C
C     Purpose: print rectangular matrix.
C
      IMPLICIT NONE
      INTEGER ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM,COLDIM,NCTL,LUPRI
      REAL*8  AMATRX(ROWDIM,COLDIM)

      CALL OUTPUT(AMATRX,ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM,COLDIM,
     &            NCTL,LUPRI)

      END
C  /* Deck cho_outpak */
      SUBROUTINE CHO_OUTPAK(AMATRX,NROW,NCTL,LUPRI)
C
C     Purpose: print packed (lower triangle) matrix.
C
      IMPLICIT NONE
      REAL*8 AMATRX(*)
      INTEGER NROW, NCTL, LUPRI

      CALL OUTPAK(AMATRX,NROW,NCTL,LUPRI)

      END
C  /* Deck cho_dzero */
      SUBROUTINE CHO_DZERO(VEC,N)
C
C     Purpose: zero double precision vector.
C
      IMPLICIT NONE
      REAL*8  VEC(*)
      INTEGER I, N

      DO I = 1,N
         VEC(I) = 0.0D0
      END DO

      END
C  /* Deck cho_izero */
      SUBROUTINE CHO_IZERO(IVEC,N)
C
C     Purpose: zero integer vector.
C
      IMPLICIT NONE
      INTEGER IVEC(*)
      INTEGER I, N

      DO I = 1,N
         IVEC(I) = 0
      END DO

      END
C  /* Deck cho_isumelm */
      INTEGER FUNCTION CHO_ISUMELM(IVEC,N)
C
C     Purpose: sum elements of integer vector.
C
      IMPLICIT NONE
      INTEGER IVEC(*)
      INTEGER I, N

      CHO_ISUMELM = 0
      DO I = 1,N
         CHO_ISUMELM = CHO_ISUMELM + IVEC(I)
      END DO

      END
C  /* Deck cho_dsumelm */
      FUNCTION CHO_DSUMELM(VEC,N)
C
C     Purpose: sum elements of double precision vector.
C
      IMPLICIT NONE
      REAL*8  CHO_DSUMELM
      REAL*8  VEC(*)
      INTEGER I, N

      CHO_DSUMELM = 0.0D0
      DO I = 1,N
         CHO_DSUMELM = CHO_DSUMELM + VEC(I)
      END DO

      END
C  /* Deck cho_isao */
      INTEGER FUNCTION CHO_ISAO(IAO)
C
C     Purpose: return symmetry of AO number IAO (in global list).
C
C     TODO/FIXME: when we are comfortable with this routine
C                 the internal consistency checks can be omitted.
C
#include <implicit.h>
#include <cholesky.h>
#include <choorb.h>

      CHARACTER*8 SECNAM
      PARAMETER (SECNAM = 'CHO_ISAO')

      INTEGER CHO_IRANGE

      CHO_ISAO = -999999

      IF ((IAO.GT.NBAST) .OR. (IAO.LT.1)) THEN
         WRITE(LUPRI,'(//,1X,A,A,I10)')
     &   SECNAM,': AO index out of bounds: ',IAO
         WRITE(LUPRI,'(A,I10,A,/)')
     &   'Maximum possible: NBAST = ',NBAST,'(from common block)'
         IF ((NBAST.LT.1) .OR. (NBAST.GT.CHO_MAXBAS)) THEN
            CALL CHO_QUIT('Initialization error detected in '//SECNAM,
     &                    102)
         ELSE
            CALL CHO_QUIT('Internal error detected in '//SECNAM,
     &                    103)
         END IF
      ELSE
         CHO_ISAO = CHO_IRANGE(IAO,IBAS,NSYM,.FALSE.)
      END IF

      END
C  /* Deck cho_isaosh */
      INTEGER FUNCTION CHO_ISAOSH(IAO,ISHL)
C
C     Purpose: return symmetry of AO number IAO in shell ISHL.
C
C     TODO/FIXME: when we are comfortable with this routine
C                 the internal consistency checks can be omitted.
C
#include <implicit.h>
#include <cholesky.h>
#include <choorb.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_ISAOSH')

      INTEGER CHO_IRANGE

      CHO_ISAOSH = -999999

      IF ((ISHL.GT.NSHELL) .OR. (ISHL.LT.1)) THEN
         WRITE(LUPRI,'(//,1X,A,A,I10)')
     &   SECNAM,': shell index out of bounds: ',ISHL
         WRITE(LUPRI,'(A,I10,A,/)')
     &   'Maximum possible: NSHELL = ',NSHELL,'(from common block)'
         IF ((NSHELL.LT.1) .OR. (NSHELL.GT.CHO_MAXSHL)) THEN
            CALL CHO_QUIT('Initialization error detected in '//SECNAM,
     &                    102)
         ELSE
            CALL CHO_QUIT('Internal error detected in '//SECNAM,
     &                    103)
         END IF
      ELSE IF ((IAO.GT.NBSTSH(ISHL)) .OR. (IAO.LT.1)) THEN
         WRITE(LUPRI,'(//,1X,A,A,I10)')
     &   SECNAM,': AO index out of bounds: ',IAO,' shell: ',ISHL
         WRITE(LUPRI,'(A,I10,A,/)')
     &   'Maximum possible: NBSTSH(ISHL) = ',NBSTSH(ISHL),
     &   '(from common block)'
         IF ((NBSTSH(ISHL).LT.1) .OR. (NBSTSH(ISHL).GT.CHO_MAXBAS)) THEN
            CALL CHO_QUIT('Initialization error detected in '//SECNAM,
     &                    102)
         ELSE
            CALL CHO_QUIT('Internal error detected in '//SECNAM,
     &                    103)
         END IF
      ELSE
         CHO_ISAOSH = CHO_IRANGE(IAO,IBASSH(1,ISHL),NSYM,.FALSE.)
      END IF

      END
C  /* Deck cho_irange */
      INTEGER FUNCTION CHO_IRANGE(I,BIN,NBIN,LOWER)
C
C     Purpose: return range to which I belongs.
C              It is assumed that BIN is sorted
C              in ascending order so that:
C
C                I > BIN(NBIN)   ==> CHO_IRANGE = NBIN
C              Else
C                I > BIN(NBIN-1) ==> CHO_IRANGE = NBIN-1
C              Else
C                I > BIN(NBIN-2) ==> CHO_IRANGE = NBIN-2
C                ..
C                ..
C                ..
C              Else
C                I < BIN(1)      ==> CHO_IRANGE = 0
C
C              If NBIN <= 0      ==> CHO_IRANGE = -1
C
C              In case of degeneracies, the smallest bin is
C              returned if flag LOWER is set, i.e.,
C
C              I > BIN(IBIN) = BIN(IBIN-1)
C                            = ...
C                            = BIN(IBIN-N) ==> CHO_IRANGE = IBIN-N
C
C              else the largest bin (IBIN in the case above).
C
      IMPLICIT NONE
      INTEGER I, NBIN, BIN(*)
      LOGICAL LOWER
      INTEGER IBIN, JBIN

      CHO_IRANGE = 0

      IF (NBIN .LT. 1) THEN
         CHO_IRANGE = -1 
      ELSE
         IBIN = NBIN
         DO WHILE (IBIN .GT. 0)
            IF (I .GT. BIN(IBIN)) THEN
               CHO_IRANGE = IBIN
               IF (LOWER) THEN
                  JBIN = IBIN - 1
                  DO WHILE (JBIN.GT.0)
                     IF (BIN(JBIN) .EQ. BIN(IBIN)) THEN
                        CHO_IRANGE = JBIN
                        JBIN = JBIN - 1
                     ELSE
                        JBIN = 0
                     END IF
                  END DO
               END IF
               IBIN = 0
            ELSE
               IBIN = IBIN - 1
            END IF
         END DO
      END IF

      END
C  /* Deck cho_head */
      SUBROUTINE CHO_HEAD(STRING,LINE,LENMAX,LUNIT)
C
C     Purpose: print a header.
C
      IMPLICIT NONE
      CHARACTER*(*) STRING
      CHARACTER*1   LINE
      INTEGER       LENMAX, LUNIT

      INTEGER LENSTR, LENTOT, I

      LENSTR = LEN(STRING)
      LENTOT = MIN(LENSTR,LENMAX-2)
      IF (LENTOT .GT. 0) THEN
         WRITE(LUNIT,'(//,2X,A)') STRING(1:LENTOT)
         WRITE(LUNIT,'(2X,80A)') (LINE,I=1,LENTOT)
      ELSE
         WRITE(LUNIT,'(//,2X,A,/)') STRING(1:LENSTR)
      END IF

      END
C  /* Deck cho_anasize */
      SUBROUTINE CHO_ANASIZE(VEC,LVEC,BIN,LBIN,LUPRI)
C
C     Purpose: analyse vector (histogram).
C
#include <implicit.h>
      DIMENSION VEC(LVEC), BIN(LBIN)

      PARAMETER (ZERO = 0.0D0)

      PARAMETER (MBIN = 20)
      INTEGER   ICOUNT(MBIN)
      LOGICAL   FOUND

C     Return if nothing to do.
C     ------------------------

      IF ((LVEC.LT.1) .OR. (LBIN.LT.1)) RETURN

C     Ensure that BIN is in descending order.
C     ---------------------------------------

      IJOB = -1
      CALL CHO_ORDER(BIN,LBIN,IJOB)

C     Test that BIN is positive.
C     --------------------------

      IF (BIN(1) .LE. ZERO) RETURN

C     Analysis.
C     ---------

      NBIN = MIN(LBIN,MBIN)
      CALL CHO_IZERO(ICOUNT,NBIN)
      NLOW = 0
      NZER = 0
      NNEG = 0
      XNEG = ZERO

      DO I = 1,LVEC

         TEST = VEC(I)

         IF (TEST .LT. ZERO) THEN
            NNEG = NNEG + 1
            XNEG = MIN(XNEG,TEST)
         ELSE IF (TEST .EQ. ZERO) THEN
            NZER = NZER + 1
         END IF

         IBIN  = 0
         FOUND = .FALSE.
         DO WHILE ((.NOT.FOUND) .AND. (IBIN.LT.NBIN))
            IBIN = IBIN + 1
            IF (TEST .GE. BIN(IBIN)) THEN
               ICOUNT(IBIN) = ICOUNT(IBIN) + 1
               FOUND = .TRUE.
            END IF
         END DO
         IF (.NOT.FOUND) NLOW = NLOW + 1

      END DO

C     Print.
C     ------

      WRITE(LUPRI,'(/,1X,A,11X,1P,D10.4,A,I12)')
     & 'Larger than ',BIN(1),':',ICOUNT(1)
      DO IBIN = 2,NBIN
         WRITE(LUPRI,'(1X,A,1P,D10.4,A,D10.4,A,I12)')
     &   'Between ',BIN(IBIN-1),' and ',BIN(IBIN),':',ICOUNT(IBIN)
      END DO
      WRITE(LUPRI,'(1X,A,10X,1P,D10.4,A,I12)')
     & 'Smaller than ',BIN(NBIN),':',NLOW

      WRITE(LUPRI,'(/,1X,A,I12)')
     & 'Number of elements exactly 0.0D0 :',NZER
      WRITE(LUPRI,'(1X,A,I12)')
     & 'Number of negative elements      :',NNEG
      IF (NNEG .GT. 0) THEN
         WRITE(LUPRI,'(1X,A,1P,D12.4)')
     & ' - numerically largest           :',XNEG
      END IF

      END
C  /* Deck cho_order */
      SUBROUTINE CHO_ORDER(VEC,LVEC,IJOB)
C
C     Purpose: sort elements of VEC in
C              IJOB = -1: descending order
C              IJOB =  1: ascending  order
C              (all other IJOB values are ignored)
C
#include <implicit.h>
      DIMENSION VEC(LVEC)

      IF (IJOB .EQ. -1) THEN

         DO I = 1,LVEC-1
            VMAX = VEC(I)
            IMAX = I
            DO J = I+1,LVEC
               IF (VEC(J) .GT. VMAX) THEN
                  VMAX = VEC(J)
                  IMAX = J
               END IF
            END DO
            IF (IMAX .NE. I) THEN
               VEC(IMAX) = VEC(I)
               VEC(I)    = VMAX
            END IF
         END DO

      ELSE IF (IJOB .EQ. 1) THEN

         DO I = 1,LVEC-1
            VMIN = VEC(I)
            IMIN = I
            DO J = I+1,LVEC
               IF (VEC(J) .LT. VMIN) THEN
                  VMIN = VEC(J)
                  IMIN = J
               END IF
            END DO
            IF (IMIN .NE. I) THEN
               VEC(IMIN) = VEC(I)
               VEC(I)    = VMIN
            END IF
         END DO

      END IF

      END
C  /* Deck cho_invpck */
      SUBROUTINE CHO_INVPCK(IJ,I,J,LOW)
C
C     Purpose: invert triangular packing index ij to
C              rectangular indices i and j.
C              Flag LOW specifies packing convention:
C              LOW = T: i>=j
C              LOW = F: i<=j
C
C     TODO/FIXME: should be included as a statement function when
C                 needed!
C
#include <implicit.h>
      LOGICAL LOW

      PARAMETER (ONE = 1.0D0, TWO = 2.0D0, THREE = 3.0D0, EIGHT = 8.0D0)

      IF (IJ .GT. 0) THEN

         XX = EIGHT*IJ - THREE
         XI = (ONE + SQRT(XX))/TWO
         I  = INT(XI)
         J  = IJ - I*(I - 1)/2

         IF (.NOT. LOW) THEN
            ITMP = I
            I = J
            J = ITMP
         END IF

      ELSE

         I = -1
         J = -2

      END IF

      END
C  /* Deck cho_invrec */
      SUBROUTINE CHO_INVREC(AB,A,B,NA,NB)
C
C     Purpose: invert rectangular index AB.
C
      IMPLICIT NONE
      INTEGER AB,A,B,NA,NB

      INTEGER IA, IB, IAB, IAB1, IAB2

      A = -999999
      B = -999999

      DO IB = 1,NB
         IAB1 = NA*(IB - 1) + 1
         IAB2 = IAB1 + NA - 1
         IF (AB.GE.IAB1 .AND. AB.LE.IAB2) THEN
            IA   = 0
            IAB  = IAB1 - 1
            DO WHILE (IA.LT.NA .AND. AB.NE.IAB)
               IA  = IA  + 1
               IAB = IAB + 1
            END DO
            IF (AB .EQ. IAB) THEN
               B = IB
               A = IA
            END IF
            RETURN
         END IF
      END DO

      END
C  /* Deck cho_wrbuf */
      SUBROUTINE CHO_WRBUF(LENGTH,BUF,IBUF,LENBUF,IUNIT)
C
C     Purpose: write buffer to disk.
C
#include <implicit.h>
      DIMENSION BUF(LENBUF)
      INTEGER   IBUF(4,LENBUF)

      WRITE(IUNIT) LENGTH,BUF,IBUF

      END
C  /* Deck cho_rdbuf */
      SUBROUTINE CHO_RDBUF(LENGTH,BUF,IBUF,LENBUF,IUNIT)
C
C     Purpose: read buffer from disk.
C
#include <implicit.h>
      DIMENSION BUF(LENBUF)
      INTEGER   IBUF(4,LENBUF)

      READ(IUNIT) LENGTH,BUF,IBUF

      END
C  /* Deck cho_open */
      SUBROUTINE CHO_OPEN(LUNIT,FNAME)
C
C     Purpose: open sequential unformatted fortran file.
C
      IMPLICIT NONE
      INTEGER  LUNIT
      CHARACTER*(*) FNAME

      INTEGER  LOCUNT, ISEED

      INTEGER ISFREEUNIT

      IF (LUNIT .GT. 1) THEN
         LOCUNT = LUNIT
      ELSE
         LOCUNT = -1
      END IF
      CALL GPOPEN(LOCUNT,FNAME,'UNKNOWN','SEQUENTIAL',
     &            'UNFORMATTED',-1,.FALSE.)

      LUNIT = LOCUNT

      END
C  /* Deck cho_close */
      SUBROUTINE CHO_CLOSE(LUNIT,STAT)
C
C     Purpose: close sequential unformatted fortran file.
C
      IMPLICIT NONE
      INTEGER  LUNIT
      CHARACTER*(*) STAT

      CALL GPCLOSE(LUNIT,STAT)
      LUNIT = -1

      END
C  /* Deck cho_prtdia */
      SUBROUTINE CHO_PRTDIA(DIAG,IIBSTRSH,NNBSTRSH,INDRED,INDRSH,ISYLST,
     &                      MSYM,MMSHL,MMBSTRT,NSYLST,IRED)
C
C     Purpose: print requested symmetry block(s) of diagonal in
C              first (IRED=1) or current (IRED=2) reduced set.
C
#include <implicit.h>
      DIMENSION DIAG(*)
      INTEGER   IIBSTRSH(MSYM,MMSHL,3), NNBSTRSH(MSYM,MMSHL,3)
      INTEGER   INDRED(MMBSTRT,3), INDRSH(MMBSTRT)
      INTEGER   ISYLST(NSYLST)
#include <cholesky.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_PRTDIA')

C     Check dimensions.
C     -----------------

      IF (NSYLST .GT. NSYM) THEN
         WRITE(LUPRI,'(//,1X,A,A)') SECNAM,': NSYLST <= NSYM required!'
         WRITE(LUPRI,'(1X,A,I10)')    'NSYLST = ',NSYLST
         WRITE(LUPRI,'(1X,A,I10,/)')  'NSYM   = ',NSYM
         CALL CHO_QUIT('[0] Symmetry error in '//SECNAM,102)
      END IF
      IF (MSYM .NE. NSYM) THEN
         WRITE(LUPRI,'(//,1X,A,A)') SECNAM,': MSYM = NSYM required!'
         WRITE(LUPRI,'(1X,A,I10)')    'MSYM = ',MSYM
         WRITE(LUPRI,'(1X,A,I10,/)')  'NSYM = ',NSYM
         CALL CHO_QUIT('[1] Symmetry error in '//SECNAM,102)
      END IF
      NNSHL = NSHELL*(NSHELL + 1)/2
      IF (MMSHL .NE. NNSHL) THEN
         WRITE(LUPRI,'(//,1X,A,A)') SECNAM,
     &   ': MMSHL = NSHELL*(NSHELL+1)/2 required!'
         WRITE(LUPRI,'(1X,A,I10)')    'MMSHL               = ',MMSHL
         WRITE(LUPRI,'(1X,A,I10,/)')  'NSHELL*(NSHELL+1)/2 = ',NNSHL
         CALL CHO_QUIT('Shell error in '//SECNAM,102)
      END IF
      IF (MMBSTRT .NE. NNBSTRT(1)) THEN
         WRITE(LUPRI,'(//,1X,A,A)') SECNAM,
     &   ': MMBSTRT = NNBSTRT(1) required!'
         WRITE(LUPRI,'(1X,A,I10)')    'MMBSTRT    = ',MMBSTRT
         WRITE(LUPRI,'(1X,A,I10,/)')  'NNBSTRT(1) = ',NNBSTRT(1)
         CALL CHO_QUIT('Dimension error in '//SECNAM,102)
      END IF

C     Code for first or current reduced set.
C     --------------------------------------

      IF (IRED .EQ. 1) THEN
         CALL CHO_HEAD(SECNAM//': Diagonal in Original Reduced Set',
     &                 '=',80,LUPRI)
         DO ILST = 1,NSYLST
            ISYM = ISYLST(ILST)
            IF (ISYM.LT.1 .OR. ISYM.GT.NSYM) THEN
               WRITE(LUPRI,*) SECNAM,': element ',ILST,': ',ISYM,
     &                        ' of list ISYLST is out of bounds!'
               CALL CHO_QUIT('ISYLST input error in '//SECNAM,104)
            ELSE
               WRITE(LUPRI,'(/,A,I2)') 'Symmetry block:',ISYM
               WRITE(LUPRI,'(/,A,/,A)')
     &         '  Element Shell-Pair  SP Index         Diagonal',
     &         '-----------------------------------------------' 
               DO ISHLAB = 1,NNSHL
                  IAB1 = IIBSTR(ISYM,IRED) + IIBSTRSH(ISYM,ISHLAB,IRED)
     &                 + 1
                  IAB2 = IAB1 + NNBSTRSH(ISYM,ISHLAB,IRED) - 1
                  DO IAB = IAB1,IAB2
                     IF (INDRSH(IAB) .NE. ISHLAB) THEN
                        WRITE(LUPRI,*)
     &                  'Shell Pair error: INDRSH,ISHLAB',
     &                  INDRSH(IAB),ISHLAB
                        CALL CHO_QUIT('Shell-Pair error in '//SECNAM,
     &                                104)
                     ELSE
                        JAB = IAB
                        WRITE(LUPRI,'(I9,2X,I9,1X,I9,1X,1P,D16.8)')
     &                  JAB,ISHLAB,INDRED(IAB,IRED),DIAG(JAB)
                     END IF
                  END DO
               END DO
               WRITE(LUPRI,'(A)')
     &         '-----------------------------------------------' 
            END IF
         END DO
      ELSE IF (IRED .EQ. 2) THEN
         CALL CHO_HEAD(SECNAM//': Diagonal in Current Reduced Set',
     &                 '=',80,LUPRI)
         DO ILST = 1,NSYLST
            ISYM = ISYLST(ILST)
            IF (ISYM.LT.1 .OR. ISYM.GT.NSYM) THEN
               WRITE(LUPRI,*) SECNAM,': element ',ILST,': ',ISYM,
     &                        ' of list ISYLST is out of bounds!'
               CALL CHO_QUIT('ISYLST input error in '//SECNAM,104)
            ELSE
               WRITE(LUPRI,'(/,A,I2)') 'Symmetry block:',ISYM
               WRITE(LUPRI,'(/,A,A,/,A,A)')
     &         '  Element  RedSet 1 Shell-Pair  SP Index',
     &         '         Diagonal',
     &         '----------------------------------------',
     &         '-----------------' 
               DO ISHLAB = 1,NNSHL
                  IAB1 = IIBSTR(ISYM,IRED) + IIBSTRSH(ISYM,ISHLAB,IRED)
     &                 + 1
                  IAB2 = IAB1 + NNBSTRSH(ISYM,ISHLAB,IRED) - 1
                  DO IAB = IAB1,IAB2
                     JAB = INDRED(IAB,IRED)
                     IF (INDRSH(JAB) .NE. ISHLAB) THEN
                        WRITE(LUPRI,*)
     &                  'Shell Pair error: INDRSH,ISHLAB',
     &                  INDRSH(JAB),ISHLAB
                        CALL CHO_QUIT('Shell-Pair error in '//SECNAM,
     &                                104)
                     ELSE
                       WRITE(LUPRI,'(I9,1X,I9,2X,I9,1X,I9,1X,1P,D16.8)')
     &                  IAB,JAB,ISHLAB,INDRED(JAB,1),DIAG(JAB)
                     END IF
                  END DO
               END DO
               WRITE(LUPRI,'(A,A)')
     &         '----------------------------------------',
     &         '-----------------' 
            END IF
         END DO
      END IF

      CALL CHO_FLUSH(LUPRI)

      END
C  /* Deck cho_tabind */
      INTEGER FUNCTION CHO_TABIND(TABLE,LKEY,NTABLE,EOINP,LEOINP,NEOINP,
     &                            WORD)
C
C     Purpose: table lookup.
C              First, try to find WORD in TABLE. If success, return ID,
C              else, check if WORD is a special string
C              in EOINP (if any supplied). If success, return NTABLE+1,
C              else, return -1.
C
      IMPLICIT NONE
      INTEGER LKEY, NTABLE, LEOINP, NEOINP
      CHARACTER*(*) TABLE(NTABLE)  ! <-- character*(lkey)
      CHARACTER*(*) EOINP(NEOINP)  ! <-- character*(leoinp)
      CHARACTER*(*) WORD           ! <-- character*(lkey)

      INTEGER IJUMP, LCMP

C     Find entry.
C     -----------

      IF (LKEY.GT.0 .AND. NTABLE.GT.0) THEN
         IJUMP = 1
         DO WHILE (IJUMP.LE.NTABLE .AND. TABLE(IJUMP).NE.WORD)
            IJUMP = IJUMP + 1
         END DO
         IF (IJUMP .GT. NTABLE) THEN
            IF (LEOINP.GT.0 .AND. NEOINP.GT.0) THEN
               LCMP  = MIN(LEOINP,LKEY)
               IJUMP = 1
               DO WHILE (IJUMP.LE.NEOINP .AND.
     &                   EOINP(IJUMP)(1:LCMP).NE.WORD(1:LCMP))
                  IJUMP = IJUMP + 1
               END DO
               IF (IJUMP .GT. NEOINP) THEN
                  CHO_TABIND = -1
               ELSE
                  CHO_TABIND = NTABLE + 1
               END IF
            ELSE
               CHO_TABIND = -1
            END IF
         ELSE
            CHO_TABIND = IJUMP
         END IF
      ELSE
         CHO_TABIND = -1
      END IF

      END
C  /* Deck cho_timer */
      SUBROUTINE CHO_TIMER(TCPU,TWALL)
C
C     Purpose: Get CPU and Wall time clock.
C
      IMPLICIT NONE
      REAL*8 TCPU, TWALL

      CALL GETTIM(TCPU,TWALL)

      END
C  /* Deck cho_setdamp */
      SUBROUTINE CHO_SETDAMP()
C
C     Purpose: set screening damping, unless user-defined.
C
      IMPLICIT NONE
#include <cholesky.h>
      INTEGER I

      DO I = 1,2
         IF (DAMP(I) .LT. 0.0D0) THEN
            IF (THRCOM .GT. 9.99D-3) THEN ! >= 1.0D-2
               DAMP(I) = 1.0D7
            ELSE IF (THRCOM .GT. 9.99D-4) THEN ! >= 1.0D-3
               DAMP(I) = 1.0D6
            ELSE IF (THRCOM .GT. 9.99D-5) THEN ! >= 1.0D-4
               DAMP(I) = 1.0D5
            ELSE IF (THRCOM .GT. 9.99D-6) THEN ! >= 1.0D-5
               DAMP(I) = 1.0D4
            ELSE IF (THRCOM .GT. 9.99D-7) THEN ! >= 1.0D-6
               DAMP(I) = 1.0D3
            ELSE IF (THRCOM .GT. 9.99D-8) THEN ! >= 1.0D-7
               DAMP(I) = 1.0D2
            ELSE IF (THRCOM .GT. 9.99D-9) THEN ! >= 1.0D-8
               DAMP(I) = 1.0D1
            ELSE
               DAMP(I) = 1.0D0
            END IF
         END IF
      END DO

      END
C  /* Deck cho_final */
      SUBROUTINE CHO_FINAL()
C
C     Purpose: decomposition finalization.
C
      IMPLICIT NONE
#include <cholesky.h>

C     Transfer info from cholesky.h to ccdeco.h.
C     ------------------------------------------

      CALL CHO_DAL_FINAL(THRCOM,DAMP,SPAN,IIBSTR,NNBSTR,NUMCHO,
     &                   RSTDIA,RSTCHO,NSYM)

      END
C  /* Deck cho_rscopy */
      SUBROUTINE CHO_RSCOPY(IIBSTRSH,NNBSTRSH,INDRED,IRS1,IRS2,
     &                      MSYM,MMSHL,LMMBSTRT,NRS)
C
C     Purpose: copy red. set info from location IRS1 to IRS2.
C              Special action is taken with INDRED if IRS1=1 so that it
C              will point as expected for the "current" reduced set.
C
      IMPLICIT NONE
      INTEGER IRS1, IRS2, MSYM, MMSHL, LMMBSTRT, NRS
      INTEGER IIBSTRSH(MSYM,MMSHL,NRS), NNBSTRSH(MSYM,MMSHL,NRS)
      INTEGER INDRED(LMMBSTRT,NRS)
#include <cholesky.h>

      INTEGER IAB

      CALL ICOPY(MSYM*MMSHL,IIBSTRSH(1,1,IRS1),1,IIBSTRSH(1,1,IRS2),1)
      CALL ICOPY(MSYM*MMSHL,NNBSTRSH(1,1,IRS1),1,NNBSTRSH(1,1,IRS2),1)
      CALL ICOPY(MSYM,IIBSTR(1,IRS1),1,IIBSTR(1,IRS2),1)
      CALL ICOPY(MSYM,NNBSTR(1,IRS1),1,NNBSTR(1,IRS2),1)
      IF (IRS1 .EQ. 1) THEN
         DO IAB = 1,LMMBSTRT
            INDRED(IAB,IRS2) = IAB
         END DO
      ELSE
         CALL ICOPY(LMMBSTRT,INDRED(1,IRS1),1,INDRED(1,IRS2),1)
      END IF
      NNBSTRT(IRS2) = NNBSTRT(IRS1)

      END
