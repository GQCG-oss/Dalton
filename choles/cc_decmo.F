C  /* Deck cc_decmo */
      SUBROUTINE CC_DECMO(XDIANL,NTOVEC,NCOLUM,INDIA,THRDCM,SPADCM,
     &                    FILSAV,SKPDCM,MXQUAL,MXREAD,NDIM,IPRINT,IMAT,
     &                    NSYM,THRZER,FOCKD,WORK,LWORK)
C
C     Henrik Koch, Alfredo Sanchez, and Thomas Bondo Pedersen, July 2002.
C
C     tbp, august 2002: generalized version essentially identical to 
C                       the original CC_DECOMP2.
C                       To add a new IMAT:
C                       - plug in a new name in array OBJECT (this routine)
C                         and set up necesary type specifiers here and
C                         in CHO_MOP for the communications via direct-access
C                         files.
C                       - plug in the relevant routine calls in CC_DCMDIA
C                         and CC_DCMMAT.
C
C     tbp, september 2002: attempt to hold original set of MO Cholesky vectors
C                          in core, thus minimizing I/O.
C
C     Purpose: Cholesky decompose a positive definite matrix
C              specified through input variable IMAT:
C
C              IMAT = 1: (ia|jb) integrals constructed
C                        from L(ia,J) vectors.
C
C              IMAT = 2: (ai|bj) integrals constructed
C                        from L(ai,J) vectors.
C
C              IMAT = 3: CC2 doubles amplitudes constructed
C                        from L(ai,J) vectors.
C                        (scaled by -1 for positiveness)
C
C     Input:
C
C        INDIA  -- Scratch index array for largest diagonals.
C                  Dimension: MXQUAL.
C
C        THRDCM -- Threshold for decomposition.
C
C        SPADCM -- Span factor.
C
C        FILSAV -- File name for saving/reading information.
C
C        SKPDCM -- Logical flag for skipping decomposition.
C                  if .TRUE.: information is read from file FILSAV.
C
C        MXQUAL -- Maximum number of diagonals qualified.
C
C        MXREAD -- Maximum number of previous vectors read.
C
C        NDIM   -- Linear dimension in each symmetry (e.g. NT1AM()
C                  for IMAT = 1,2, or 3; included for facilitating
C                  future IMAT-implementations).
C                  Dimension: NSYM.
C
C        IPRINT -- Local print level.
C
C        NSYM   -- Number of symmetries. Should be as in ccorb.h
C                  (not tested).
C
C        THRZER -- Threshold for zeroing out diagonals.
C                  No screening for negative THRZER.
C                  Screening formula (NOTE the square root!!):
C                  SQRT(DIAG(i)*MAXDIAG) < THRZER => DIAG(i)=0
C
C        FOCKD  -- Canonical orbital energies (IMAT=3 only).
C
C        WORK   -- Work space.
C                  Dimension: LWORK.
C
C     Output:
C
C        XDIANL -- Array containing the min. error, max. error, and the
C                  rms error for the diagonal of each symmetry.
C                  Dimension: (3,NSYM).
C
C        NTOVEC -- Number of Cholesky vectors in each symmetry.
C                  Dimension: NSYM.
C
C        NCOLUM -- Number of columns actually calculated for decomposition.
C                  Dimension: NSYM.
C
C     Additional:
C
C        The diagonal and integral drivers, which are used as "plug-ins",
C        are assumed to recognize IMAT.
C
C        The Cholesky vectors are returned on disk on direct-access files
C        as set in routine CHO_MOP.
C
#include <implicit.h>
      DIMENSION XDIANL(3,NSYM), FOCKD(*), WORK(LWORK)
      INTEGER   NTOVEC(NSYM), NCOLUM(NSYM), INDIA(NSYM), NDIM(NSYM)
      LOGICAL   SKPDCM
      CHARACTER*(*) FILSAV
#include <priunit.h>
#include <iratdef.h>
#include <maxorb.h>
#include <ccdeco.h>

      LOGICAL LAST, ZEROUT

      INTEGER INCORE(8)

      CHARACTER*8 SECNAM
      PARAMETER (SECNAM = 'CC_DECMO')

      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
      PARAMETER (INFO = 0, INFC = 1, INFI = 2, INFD = 7)
      PARAMETER (IOPEN = -1, IKEEP = 1)
      PARAMETER (MXLIST = 16)

      PARAMETER (NMAT = 3)
      CHARACTER*7 OBJECT(NMAT)
      DATA OBJECT /'(ia|jb)','(ai|bj)','CC2 am.'/

C     Start timing.
C     -------------

      TIMTOT = SECOND()
      TIMDIA = ZERO
      TIMINT = ZERO

C     Test IMAT.
C     ----------

      IF ((IMAT.LE.0) .OR. (IMAT.GT.NMAT)) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Error in ',SECNAM,' - specifier IMAT not recognized:'
         WRITE(LUPRI,'(5X,A,I10)')
     &   'IMAT = ',IMAT
         WRITE(LUPRI,'(5X,A)')
     &   'Implemented specifiers:'
         DO JMAT = 1,NMAT
            WRITE(LUPRI,'(5X,I3,1X,A,A,A)') JMAT,'[',OBJECT(JMAT),']'
         ENDDO
         WRITE(LUPRI,*)
         CALL QUIT('Error in '//SECNAM)
      ENDIF

C     For skip option, read info from disk.
C     -------------------------------------

      IF (SKPDCM) THEN
         CALL RD_DECMOS(THRDCM,SPADCM,MXQUAL,MXREAD,NTOVEC,NCOLUM,
     &                  XDIANL,NSYM,THRZER,FILSAV,SECNAM)
      ENDIF

C     Set logical for zeroing out.
C     ----------------------------

      ZEROUT = THRZER .GE. ZERO

C     Print.
C     ------

      IF (IPRINT .GE. INFO) THEN
         CALL AROUND('Output from '//SECNAM)
         IF (SKPDCM) THEN
            WRITE(LUPRI,'(10X,A,A)')
     &      OBJECT(IMAT),' decomposition skipped.'
            WRITE(LUPRI,'(10X,A,A,A)')
     &      'Information read from file ',FILSAV,':'
         ELSE
            WRITE(LUPRI,'(10X,A,A,A)')
     &      'Cholesky decomposition of ',OBJECT(IMAT),':'
         ENDIF
         WRITE(LUPRI,'(10X,A,1P,D15.6)')
     &   'Decomposition threshold: ',THRDCM
         WRITE(LUPRI,'(10X,A,1P,D15.6)')
     &   'Span factor            : ',SPADCM
         IF (ZEROUT) THEN
            WRITE(LUPRI,'(10X,A,1P,D15.6)')
     &      'Threshold for zeroing  : ',THRZER
         ELSE
            WRITE(LUPRI,'(10X,A)')
     &      'No zeroing of diagonals.'
         ENDIF
         WRITE(LUPRI,'(10X,A,I10)')
     &   'Max. number of diagonals qualified: ',MXQUAL
         WRITE(LUPRI,'(10X,A,I10,/)')
     &   'Max. number of prev. vectors read : ',MXREAD
         CALL FLSHFO(LUPRI)
      ENDIF

C     Skip if requested.
C     ------------------

      IF (SKPDCM) GOTO 999

C     Set up stuff for analysis of diagonal.
C     --------------------------------------

      IF (IPRINT .GE. INFI) THEN

         KLIST = 1
         KEND0 = KLIST + MXLIST
         LWRK0 = LWORK - KEND0 - 1

         IF (LWRK0 .LT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Insufficient memory in ',SECNAM,' - allocation: list'
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need (more than): ',KEND0-1,
     &      'Available       : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

         WORK(KLIST) = 1.0D0
         DO ILIST = 2,MXLIST
            KOFF = KLIST + ILIST - 1
            WORK(KOFF) = WORK(KOFF-1)/10.0D0
         ENDDO

         NLIST = MXLIST
         DO ILIST = MXLIST-1,1,-1
            KOFF = KLIST + ILIST - 1
            IF (THRDCM .GE. WORK(KOFF)) THEN
               NLIST = ILIST + 1
            ELSE
               GO TO 10
            ENDIF
         ENDDO
   10    CONTINUE

      ELSE

         KEND0 = 1
         LWRK0 = LWORK
         KLIST = 1
         NLIST = 0

      ENDIF

C     Error test.
C     -----------

      IF ((MXQUAL.LE.0) .OR. (MXREAD.LE.0)) THEN
         WRITE(LUPRI,'(//,5X,A,A,A,A)')
     &   SECNAM,': unable to decompose ',OBJECT(IMAT),':'
         IF (MXQUAL .LE. 0) THEN
            WRITE(LUPRI,'(5X,A,I10)')
     &      'Number of qualified diagonals too small: MXQUAL = ',
     &      MXQUAL
         ENDIF
         IF (MXREAD .LE. 0) THEN
            WRITE(LUPRI,'(5X,A,I10)')
     &      'Number of prev. vectors read too small: MXREAD = ',
     &      MXREAD
         ENDIF
         CALL QUIT('Error in '//SECNAM)
      ENDIF

C     Set type identifiers for CHO_MOP.
C     ---------------------------------

      IF (IMAT .EQ. 1) THEN
         ITYMAT = 1
         ITYVEC = 2
      ELSE IF (IMAT .EQ. 2) THEN
         ITYMAT = 3
         ITYVEC = 5
      ELSE IF (IMAT .EQ. 3) THEN
         ITYMAT = 3
         ITYVEC = 6
      ENDIF

C     Loop over symmetries to decompose.
C     ----------------------------------

      DO ISYM = 1,NSYM

         TIMSYM = SECOND()

C        Print.
C        ------

         IF (IPRINT .GE. INFI) THEN
            WRITE(LUPRI,'(5X,A,A,A,I1,A)')
     &      'Cholesky decomposition of ',OBJECT(IMAT),
     &      ' symmetry block ',ISYM,':'
            IF (NDIM(ISYM) .LE. 0) THEN
               WRITE(LUPRI,'(5X,A,I10,/)')
     &         '- empty block: NDIM = ',NDIM(ISYM)
            ENDIF
            CALL FLSHFO(LUPRI)
         ENDIF

C        Initialization.
C        ---------------

         NTOVEC(ISYM) = 0
         NCOLUM(ISYM) = 0

C        Skip if no elements.
C        --------------------

         IF (NDIM(ISYM) .LE. 0) GOTO 100

C        Open files.
C        -----------

         CALL CHO_MOP(IOPEN,ITYMAT,ISYM,LUCHMO,1,1)
         CALL CHO_MOP(IOPEN,ITYVEC,ISYM,LUDECM,1,1)

C        Calculate 10% of the matrix to be decomposed,
C        although never more than MXQUAL.
C        ---------------------------------------------

         IF (NDIM(ISYM) .GT. 10) THEN
            NUMBER = NDIM(ISYM)/10
            MDECOM = MIN(MXQUAL,NUMBER)
         ELSE
            MDECOM = MIN(MXQUAL,NDIM(ISYM))
         ENDIF

C        Test if we can hold original set of MO vectors
C        in core for this symmetry.
C        N.B. CHANGE THIS IF YOU CHANGE ANY ALLOCATIONS
C             IN THIS ROUTINE!!!!
C        NEED = diagonal
C             + pointer array
C             + integral columns
C             + prev./new vectors
C             + scratch
C        ----------------------------------------------

         NEED = NDIM(ISYM)
     &        + (NDIM(ISYM) - 1)/IRAT + 1
     &        + NDIM(ISYM)*MDECOM
     &        + NDIM(ISYM)*MXREAD
     &        + MXREAD*MDECOM
         XNED = ONE*NEED

         XDIM  = ONE*NDIM(ISYM)
         XCHO  = ONE*NUMCHO(ISYM)
         XDEC  = ONE*MDECOM
         XWORK = ONE*LWORK

         KORIG = KEND0

         XORIG = ONE*KORIG
         XSORT = XORIG + XDIM*XCHO
         XENDT = XSORT + XCHO*XDEC
         XWRKT = XWORK - XENDT + ONE

         IF (XWRKT .GT. XNED) THEN
            INCORE(ISYM) = 1
            KSORT = KORIG + NDIM(ISYM)*NUMCHO(ISYM)
            KENDT = KSORT + NUMCHO(ISYM)*MDECOM
            TLINC = SECOND()
            CALL CHO_MOREAD(WORK(KORIG),NDIM(ISYM),NUMCHO(ISYM),1,
     &                      LUCHMO)
            TLINC = SECOND() - TLINC
            IF (IPRINT .GE. INFI) THEN
               WRITE(LUPRI,'(5X,A,/,5X,A,F10.2,A)')
     &        'Full set of original MO Cholesky vectors stored in core',
     &         'Time used for reading: ',TLINC,' seconds'
               CALL FLSHFO(LUPRI)
            ENDIF
         ELSE
            INCORE(ISYM) = 0
            KSORT = KORIG
            KENDT = KEND0
         ENDIF

C        Initialize record counter.
C        --------------------------

         IREC = 1

C        Get the diagonal.
C        -----------------

         KDIAG = KENDT
         KPOIN = KDIAG + NDIM(ISYM)
         KEND1 = KPOIN + (NDIM(ISYM) - 1)/IRAT + 1
         LWRK1 = LWORK - KEND1

         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Insufficient memory in ',SECNAM,' - allocation: diagonal'
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need (more than): ',KEND1,
     &      'Available       : ',LWORK
            IF (INCORE(ISYM) .EQ. 1) CALL CC_DECMO_E(LUPRI)
            CALL QUIT('Not enough space in '//SECNAM)
         ENDIF

         DTIME  = SECOND()
         CALL CC_DCMDIA(WORK(KDIAG),FOCKD,WORK(KORIG),WORK(KEND1),LWRK1,
     &                  ISYM,LUCHMO,IMAT,INCORE(ISYM))
         DTIME  = SECOND() - DTIME
         TIMDIA = TIMDIA   + DTIME

C        Print analysis of initial diagonal.
C        -----------------------------------

         IF (IPRINT .GE. INFI) THEN
            WRITE(LUPRI,'(/,5X,A,I1,A)')
     &      'Initial diagonal, symmetry ',ISYM,':'
            CALL CC_DCMANL(WORK(KDIAG),NDIM(ISYM),WORK(KLIST),NLIST)
         ENDIF

C        Start of iterative loop.
C        ------------------------

         ITER = 0
  200    CONTINUE

C        If converged, go to next symmetry.
C        (After closing files at line 150.)
C        ----------------------------------

         XM = ZERO
         DO I = 1,NDIM(ISYM)
            IF (WORK(KDIAG+I-1) .GT. XM) XM = WORK(KDIAG+I-1)
         ENDDO
C
         DIAMIN = XM * SPADCM
         IF (DIAMIN .LE. THRDCM) DIAMIN = THRDCM
C
         IF (ABS(XM) .LE. THRDCM) THEN
           CALL CC_DCMDIA(WORK(KDIAG),FOCKD,WORK(KORIG),
     &                    WORK(KEND1),LWRK1,
     &                    ISYM,LUCHMO,IMAT,INCORE(ISYM))
           IF ((NTOVEC(ISYM).LE.NUMCHO(ISYM)) .AND.
     &         (INCORE(ISYM).EQ.1)) THEN
              CALL CHO_MOREAD(WORK(KORIG),NDIM(ISYM),NTOVEC(ISYM),1,
     &                        LUDECM)
              JNCORE = 1
           ELSE
              JNCORE = 0
           ENDIF
           CALL CC_DCMXDI(XDIANL(1,ISYM),WORK(KDIAG),WORK(KORIG),
     &                    WORK(KEND1),LWRK1,NTOVEC(ISYM),NDIM(ISYM),
     &                    LUDECM,JNCORE)
           IF (IPRINT .GE. INFI) THEN
              WRITE(LUPRI,'(/,5X,A,I1,A)')
     &        'Final diagonal, symmetry ',ISYM,':'
              CALL CC_DCMANL(WORK(KDIAG),NDIM(ISYM),WORK(KLIST),NLIST)
           ENDIF
           IF (IPRINT .GE. INFC) THEN
              TIMSYM = SECOND() - TIMSYM
              WRITE(LUPRI,'(/,5X,A,A,A,I1,A)')
     &        'Cholesky decomposition of ',OBJECT(IMAT),
     &        ' completed for symmetry ',ISYM,':'
              IF (INCORE(ISYM) .EQ. 1) THEN
                 WRITE(LUPRI,'(5X,A)')
     &           'Original MO Cholesky vectors held in core.'
              ELSE
                 WRITE(LUPRI,'(5X,A)')
     &           'Original MO Cholesky vectors read in when needed.'
              ENDIF
              WRITE(LUPRI,'(5X,A,5X,I10)')
     &        'Total number of distributions     : ',NDIM(ISYM)
              WRITE(LUPRI,'(5X,A,5X,I10)')
     &        'Number of distributions calculated: ',NCOLUM(ISYM)
              WRITE(LUPRI,'(5X,A,5X,I10)')
     &        'Number of Cholesky vectors stored : ',NTOVEC(ISYM)
              IF (ABS(XDIANL(2,ISYM)) .GT. THRDCM) THEN
                 WRITE(LUPRI,'(5X,A,1P,D15.6,A)')
     &           'Largest diagonal                  : ',XDIANL(2,ISYM),
     &           '  *** WARNING ***'
              ELSE
                 WRITE(LUPRI,'(5X,A,1P,D15.6)')
     &           'Largest diagonal                  : ',XDIANL(2,ISYM)
              ENDIF
              WRITE(LUPRI,'(5X,A,5X,F10.2,A,/)')
     &        'Total time for this symmetry      : ',TIMSYM,' seconds'
              CALL FLSHFO(LUPRI)
           ENDIF
           GOTO 150
         ENDIF

C        Find MDECOM largest diagonal elements.
C        Leave pointers in INDIA.
C        --------------------------------------

         CALL CC_DCMLRG(WORK(KDIAG),MDECOM,INDIA,WORK(KPOIN),NDIM(ISYM))

C        Remove those that are already converged.
C        ----------------------------------------

         I = MDECOM + 1
   50    CONTINUE
         I = I - 1
         IF (MDECOM .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,I10,A,A)')
     &      'MDECOM = ',MDECOM,' in ',SECNAM
            WRITE(LUPRI,'(5X,A,A,A,/)')
     &      '- decomposition of ',OBJECT(IMAT),
     &      ' should be done ?!?!'
            CALL QUIT('Error in '//SECNAM)
         ENDIF
         IF (WORK(KDIAG+INDIA(I)-1) .LE. THRDCM) THEN
            MDECOM = MDECOM - 1
            GOTO 50
         ENDIF

         IF (IPRINT .GE. INFI) THEN
            KOFF1 = KDIAG + INDIA(MDECOM) - 1
            KOFF2 = KDIAG + INDIA(1)      - 1
            WRITE(LUPRI,'(5X,A,I10,A)')
     &      'Iteration ',ITER,':'
            WRITE(LUPRI,'(6X,A,I10,/,6X,A,1P,D15.6,1X,D15.6)')
     &      'Number of diagonals qualified  : ',MDECOM,
     &      'Smallest and largest among them: ',WORK(KOFF1),WORK(KOFF2)
            CALL FLSHFO(LUPRI)
         ENDIF

C        Update iteration counter.
C        -------------------------

         ITER = ITER + 1
         IF (ITER .GT. NDIM(ISYM)) THEN
            WRITE(LUPRI,*) SECNAM,': Iteration counter out of bounds.'
            WRITE(LUPRI,*) 'ITER = ',ITER,' > NDIM(ISYM) = ',NDIM(ISYM)
            CALL QUIT('Error in '//SECNAM)
         ENDIF

C        Allocation.
C        -----------

         KINT1 = KEND1
         KCHO1 = KINT1 + NDIM(ISYM)*MDECOM
         KCHO2 = KCHO1 + NDIM(ISYM)*MXREAD
         KEND2 = KCHO2 + MXREAD*MDECOM
         LWRK2 = LWORK - KEND2 + 1
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Insufficient memory in ',SECNAM,' - allocation: integrals'
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need     : ',KEND2-1,
     &      'Available: ',LWORK
            IF (INCORE(ISYM) .EQ. 1) CALL CC_DECMO_E(LUPRI)
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Calculate integral columns corresponding to indices in INDIA.
C        -------------------------------------------------------------

         DTIME  = SECOND()
         CALL CC_DCMMAT(WORK(KINT1),FOCKD,WORK(KORIG),WORK(KSORT),
     &                  WORK(KEND2),LWRK2,
     &                  INDIA,MDECOM,ISYM,LUCHMO,IMAT,INCORE(ISYM))
         DTIME  = SECOND() - DTIME
         TIMINT = TIMINT   + DTIME

         NCOLUM(ISYM) = NCOLUM(ISYM) + MDECOM

C        Batch over previous vectors.
C        ----------------------------

         NBATV = (NTOVEC(ISYM)-1)/MXREAD + 1

         IBATV2 = 0
         DO IBATV = 1,NBATV

            IBATV1 = IBATV2 + 1
            IBATV2 = IBATV2 + MXREAD
            IF (IBATV2 .GT. NTOVEC(ISYM)) IBATV2 = NTOVEC(ISYM)
            NUMVEC = IBATV2 - IBATV1 + 1
            IF (NUMVEC .EQ. 0) GOTO 300

            CALL CHO_MOREAD(WORK(KCHO1),NDIM(ISYM),NUMVEC,IBATV1,
     &                      LUDECM)

C           Prepare factor array.
C           ---------------------

            DO J = 1,MDECOM
               DO JJ = 1,NUMVEC
                  KOFF1 = KCHO1 + NDIM(ISYM)*(JJ-1) + INDIA(J) - 1
                  KOFF2 = KCHO2 + NUMVEC*(J-1) + JJ - 1
                  WORK(KOFF2) = WORK(KOFF1)
               ENDDO
            ENDDO

C           Subtraction.
C           ------------

            NTOT = MAX(NDIM(ISYM),1)

            CALL DGEMM('N','N',NDIM(ISYM),MDECOM,NUMVEC,
     &                 -ONE,WORK(KCHO1),NTOT,WORK(KCHO2),
     &                 NUMVEC,ONE,WORK(KINT1),NTOT)

         END DO
  300    CONTINUE

C        Decompose vectors in core.
C        --------------------------

         IDUMP = 0
         LAST  = .FALSE.
         DO 400 ICHO = 1,MDECOM

            XC = -ONE
            DO I = 1,MDECOM
               XXX = WORK(KDIAG + INDIA(I) - 1)
               IF (XXX .GT. XC) THEN
                  XC   = XXX
                  IC   = I
                  ICAI = INDIA(I)
               ENDIF
            ENDDO

            IF (IPRINT .GE. INFD) THEN
               WRITE(LUPRI,'(8X,A,I10,1X,1P,D15.6)')
     &         'Diagonal to be explicitly treated    : ',ICAI,XC
               CALL FLSHFO(LUPRI)
            ENDIF

            IF ((XC.LT.DIAMIN) .OR. (XC.LE.THRDCM)) THEN
               LAST = .TRUE.
               GOTO 500
            ENDIF

            XD   = ONE/SQRT(XC)
            KOFF = KINT1 + NDIM(ISYM)*(IC-1)
            CALL DSCAL(NDIM(ISYM),XD,WORK(KOFF),1)

            DO I = 1,NDIM(ISYM)
               KOFF = KINT1 + NDIM(ISYM)*(IC-1) + I - 1
               IF (WORK(KDIAG+I-1) .EQ. ZERO) WORK(KOFF) = ZERO
            ENDDO

C           Update diagonal.
C           ----------------

            XM   = ZERO
            NNEG = 0
            DO I = 1,NDIM(ISYM)
               KOFF = KINT1 + NDIM(ISYM)*(IC-1) + I - 1
               WORK(KDIAG+I-1) = WORK(KDIAG+I-1) - WORK(KOFF)*WORK(KOFF)
               IF (WORK(KDIAG+I-1) .GT. XM) XM = WORK(KDIAG+I-1)
               IF (WORK(KDIAG+I-1) .LT. ZERO) NNEG = NNEG + 1
            ENDDO

C           Set explicitly treated diagonal to zero.
C           ----------------------------------------

            DTST = DABS(WORK(KDIAG+ICAI-1))
            IF ((IPRINT.GE.INFD) .OR. (DTST.GT.THRDCM)) THEN
               WRITE(LUPRI,'(8X,A,I10,1X,1P,D15.6)')
     &         'Dia. explicitly treated (set to 0.D0): ',
     &         ICAI,WORK(KDIAG+ICAI-1)
               CALL FLSHFO(LUPRI)
            ENDIF

            WORK(KDIAG+ICAI-1) = ZERO

C           Zero out diagonals contributing less than THRZER
C           to integrals (if requested).
C           ------------------------------------------------

            IF (ZEROUT) THEN
               NZER = 0
               DO I = 1,NDIM(ISYM)
                  XCONTR = ABS(WORK(KDIAG+I-1)*XM)
                  XCONTR = DSQRT(XCONTR)
                  IF (XCONTR .LT. THRZER) THEN
                     NZER = NZER + 1
                     WORK(KDIAG+I-1) = ZERO
                  END IF
               ENDDO
            ENDIF

            IF (IPRINT .GE. INFD) THEN
               WRITE(LUPRI,'(8X,A,I10)')
     &         'Number of negative updated diagonals : ',NNEG
               IF (ZEROUT) THEN
                  WRITE(LUPRI,'(8X,A,I10)')
     &            'Number of  zeroed  updated diagonals : ',NZER
               ENDIF
               CALL FLSHFO(LUPRI)
            ENDIF

            DIAMIN = XM * SPADCM
            IF (DIAMIN .LE. THRDCM) DIAMIN = THRDCM

            KOFF1 = KINT1 + NDIM(ISYM)*(IC-1)
            DO I = 1,MDECOM
               II = INDIA(I)
               IF (WORK(KDIAG+II-1) .NE. ZERO) THEN

                  KOFF2 = KINT1 + NDIM(ISYM)*(I-1)

                  FACTOR = -WORK(KOFF1+II-1)

                  CALL DAXPY(NDIM(ISYM),FACTOR,
     *                       WORK(KOFF1),1,WORK(KOFF2),1)

               ENDIF
            ENDDO

C           Write cholesky vectors to disk.
C           -------------------------------

            IDUMP        = IDUMP + 1
            NTOVEC(ISYM) = NTOVEC(ISYM) + 1

            KOFF1 = KINT1 + NDIM(ISYM)*(IC-1)
            KOFF2 = KCHO1 + NDIM(ISYM)*(IDUMP-1)
            CALL DCOPY(NDIM(ISYM),WORK(KOFF1),1,WORK(KOFF2),1)

  500       CONTINUE

            IF ((IDUMP.EQ.MXREAD) .OR. (ICHO.EQ.MDECOM)
     *                            .OR. LAST) THEN

               JVEC1 = IREC
               CALL CHO_MOWRITE(WORK(KCHO1),NDIM(ISYM),IDUMP,JVEC1,
     &                          LUDECM)
               IREC  = IREC + IDUMP

               IF (IPRINT .GE. INFI) THEN
                  WRITE(LUPRI,'(6X,A,I10)')
     &            'Number of vectors written to disk    : ',IDUMP
                  WRITE(LUPRI,'(6X,A,I10)')
     &            'Accumulated number of vectors on disk: ',NTOVEC(ISYM)
                  CALL FLSHFO(LUPRI)
               ENDIF

               IF (LAST) GOTO 200  ! Go to next iteration
               IDUMP = 0
            ENDIF

  400    CONTINUE

C        Go to next iteration.
C        ---------------------

         GOTO 200

C        Close files and release units.
C        ------------------------------

  150    CALL CHO_MOP(IKEEP,ITYVEC,ISYM,LUDECM,1,1)
         CALL CHO_MOP(IKEEP,ITYMAT,ISYM,LUCHMO,1,1)

  100    CONTINUE

      ENDDO

C     Save information on disk.
C     -------------------------

      CALL WR_DECMOS(THRDCM,SPADCM,MXQUAL,MXREAD,NTOVEC,NCOLUM,
     &               XDIANL,NSYM,THRZER,FILSAV)

C     Print info.
C     -----------

  999 IF (IPRINT .GE. INFO) THEN

         TIMTOT = SECOND() - TIMTOT

         CALL HEADER('Summary of '//OBJECT(IMAT)//' decomposition',-1)

         NWRN  = 0
         IFAIL = 0
         INEG  = 0
         DO ISYM = 1,NSYM
            IF (DABS(XDIANL(1,ISYM)) .GT. THRDCM) INEG = INEG + 1
            IF (XDIANL(2,ISYM) .GT. THRDCM) IFAIL = IFAIL + 1
         ENDDO

         WRITE(LUPRI,'(5X,A,A)')
     &   'Symmetry        Min. Diag.           Max. Diag.',
     &   '            RMS Error   '
         WRITE(LUPRI,'(5X,A,A)')
     &   '---------------------------------------------',
     &   '-----------------------'
         DO ISYM = 1,NSYM
            WRITE(LUPRI,'(9X,I1,6X,1P,D15.6,6X,1P,D15.6,6X,1P,D15.6)')
     &      ISYM,XDIANL(1,ISYM),XDIANL(2,ISYM),XDIANL(3,ISYM)
         ENDDO
         IF ((IFAIL.EQ.0) .AND. (INEG.EQ.0)) THEN
            WRITE(LUPRI,'(5X,A,A,/)')
     &      '---------------------------------------------',
     &      '-----------------------'
         ELSE
            WRITE(LUPRI,'(5X,A,A)')
     &      '---------------------------------------------',
     &      '-----------------------'
            IF (IFAIL .NE. 0) THEN
               WRITE(LUPRI,'(5X,A,I2,A)')
     &         '**** WARNING **** Decomposition failed for',IFAIL,
     &         ' symmetries!'
               NWRN = NWRN + 1
            ENDIF
            IF (INEG .NE. 0) THEN
               WRITE(LUPRI,'(5X,A,I2,A)')
     &         '**** WARNING **** Large neg. diagonals for',INEG,
     &         ' symmetries!'
               NWRN = NWRN + 1
            ENDIF
            IF (ZEROUT) THEN
               WRITE(LUPRI,'(5X,A,1P,D15.6)')
     &         'The threshold for zeroing was: ',THRZER
               IF (THRZER .GT. THRDCM) THEN
                  WRITE(LUPRI,'(5X,A,A)')
     &            'Zeroing threshold is larger than decomposition',
     &            ' threshold!'
               ENDIF
               WRITE(LUPRI,'(5X,A,/)')
     &         'Program continues nevertheless...'
            ELSE
               WRITE(LUPRI,'(5X,A,/)')
     &         'No zeroing was performed - program stops!'
               CALL QUIT(SECNAM//': Decomposition failure.')
            ENDIF
         ENDIF

         WRITE(LUPRI,'(5X,A)')
     &   'Symmetry      NDIM        NTOVEC   %Saving(*)    NUMCHO'
         WRITE(LUPRI,'(5X,A)')
     &   '-------------------------------------------------------'
         XFC = 0.0D0
         XCT = 0.0D0
         XAT = 0.0D0
         DO ISYM = 1,NSYM
            IF (NDIM(ISYM) .NE. 0) THEN
               XAI = 1.0D0*NDIM(ISYM)
               XMO = 1.0D0*NTOVEC(ISYM)
               SAV = 100.0D0*(XAI-XMO)/XAI
               XAO = 1.0D0*NUMCHO(ISYM)
               XCL = 1.0D0*NCOLUM(ISYM)
               TMP = (XCL*XAO+XAI*XMO)/(XAI*XAO)
               XFC = XFC + TMP
               XCT = XCT + XAI*XCL
               XAT = XAT + XAI*XAI
            ELSE
               SAV = -3.0D10
            ENDIF
            WRITE(LUPRI,'(9X,I1,4X,I10,3X,I10,3X,F7.2,3X,I10)')
     &      ISYM,NDIM(ISYM),NTOVEC(ISYM),SAV,NUMCHO(ISYM)
         ENDDO
         WRITE(LUPRI,'(5X,A)')
     &   '-------------------------------------------------------'
         WRITE(LUPRI,'(5X,A,/)')
     &   '(*) (NDIM-NTOVEC)/NDIM'

         IF (IMAT .EQ. 1) THEN
            XSYM = 1.0D0*NSYM
            XFC  = XSYM/XFC
            IF (XFC .GE. 1.0D0) THEN
               WRITE(LUPRI,'(5X,A,F10.2,/)')
     &  'Estimated total MP2 operation count reduced by a factor of ',
     &         XFC
            ELSE
               XFC = 1.0D0/XFC
               WRITE(LUPRI,'(5X,A,F10.2,/)')
     &  'Estimated total MP2 operation count increased by a factor of ',
     &         XFC
            ENDIF
         ENDIF

         IF (.NOT. SKPDCM) THEN
            ICOUNT = 0
            DO ISYM = 1,NSYM
               IF (INCORE(ISYM) .EQ. 1) THEN
                  ICOUNT = ICOUNT + 1
                  INCORE(ICOUNT) = ISYM
               ENDIF
            ENDDO
            IF (ICOUNT .GT. 0) THEN
               WRITE(LUPRI,'(5X,A,A)')
     &         'Original MO Cholesky vectors held in core for',
     &         ' the following symmetries:'
               WRITE(LUPRI,'(5X,8I2)') (INCORE(I), I = 1,ICOUNT)
               WRITE(LUPRI,*)
            ENDIF
         ENDIF

         IF (XAT .NE. ZERO) THEN
            XCP = 100.0D0*XCT/XAT
         ELSE
            XCP = -3.0D10
         ENDIF
         IF (SKPDCM) THEN
            WRITE(LUPRI,'(5X,F7.2,A,/)')
     &      XCP,'% integrals calculated'
         ELSE
            CALL HEADER('Timing of '//OBJECT(IMAT)//' decomposition',-1)
            TIMDEC = TIMTOT - TIMDIA - TIMINT
            WRITE(LUPRI,'(5X,A,F10.2,A)')
     &      'Diagonal time     : ',TIMDIA,' seconds'
            WRITE(LUPRI,'(5X,A,F10.2,A,F7.2,A)')
     &      'Integral time     : ',TIMINT,' seconds (',XCP,
     &      '% integrals calculated)'
            WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &      'Decomposition time: ',TIMDEC,' seconds'
         ENDIF

         WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &   'Total time        : ',TIMTOT,' seconds'

         IF (NWRN .NE. 0) THEN
            WRITE(LUPRI,'(5X,A)')
     &      'Check output for warnings!'
         ENDIF
         WRITE(LUPRI,'(5X,A,A,/)')
     &   '- End of ',SECNAM

         CALL FLSHFO(LUPRI)

      ENDIF

      RETURN
      END
C  /* Deck cc_dcmlrg */
      SUBROUTINE CC_DCMLRG(DIAG,NUM,INDIA,ISCR,NDIM)
C
C     Adapted version of SO_SORT by K. L. Bak.
C
C     Purpose:
C        Find NUM largest elements in DIAG and leave pointers
C        to these elements in INDIA. NDIM is the linear dimension
C        of array DIAG.
C
#include <implicit.h>
      DIMENSION DIAG(NDIM)
      INTEGER INDIA(NUM),ISCR(NDIM)

C     Initialize scratch pointer.
C     ---------------------------

      DO I = 1,NDIM
         ISCR(I) = I
      ENDDO

C     Find the NUM largest elements of DIAG by
C     sorting indices in ISCR according to size of DIAG.
C     --------------------------------------------------

      DO J = 1,NUM
         DO I = NDIM,J+1,-1
            IF (DIAG(ISCR(I)) .GT. DIAG(ISCR(I-1))) THEN
               ITMP      = ISCR(I-1)
               ISCR(I-1) = ISCR(I)
               ISCR(I)   = ITMP
            ENDIF
         ENDDO
      ENDDO

C     Copy pointers to INDIA.
C     -----------------------

      DO I = 1,NUM
         INDIA(I) = ISCR(I)
      ENDDO

      RETURN
      END
C  /* Deck cc_dcmdia */
      SUBROUTINE CC_DCMDIA(DIAG,FOCKD,CHIA,WORK,LWORK,ISYM,LUCHMO,IMAT,
     &                     INCORE)
C
C     Thomas Bondo Pedersen, August 2002.
C
C     Purpose: Calculate diagonal for MO decomposition.
C
C     INCORE = 1: CHIA contains the Cholesky vectors for constructing
C                 the diagonal.
C
#include <implicit.h>
      DIMENSION DIAG(*), FOCKD(*), CHIA(*), WORK(LWORK)
#include <priunit.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_DCMDIA')

      IF ((IMAT.EQ.1) .OR. (IMAT.EQ.2)) THEN

C        (ia|jb) or (ai|bj) integral diagonals.
C        --------------------------------------

         IF (INCORE .EQ. 1) THEN
            CALL CHO_AIBJD_I(DIAG,CHIA,ISYM)
         ELSE
            CALL CHO_AIBJD(DIAG,WORK,LWORK,ISYM,LUCHMO)
         ENDIF

      ELSE IF (IMAT .EQ. 3) THEN

C        Minus the CC2 doubles amplitudes diagonal.
C        ------------------------------------------

         IF (INCORE .EQ. 1) THEN
            CALL CC2_CHOAMD_I(DIAG,FOCKD,CHIA,ISYM)
         ELSE
            CALL CC2_CHOAMD(DIAG,FOCKD,WORK,LWORK,ISYM,LUCHMO)
         ENDIF

      ELSE

C        Error.
C        ------

         WRITE(LUPRI,'(//,5X,A,I10,A,A)')
     &   'IMAT = ',IMAT,' option not implemented in ',SECNAM
         CALL QUIT('Unknown option in '//SECNAM)

      ENDIF

      RETURN
      END
C  /* Deck cc_dcmmat */
      SUBROUTINE CC_DCMMAT(XINT,FOCKD,CHIA,CHJB,WORK,LWORK,INDIA,MDECOM,
     &                     ISYM,LUCHMO,IMAT,INCORE)
C
C     Thomas Bondo Pedersen, August 2002.
C
C     Purpose: Calculate matrix (XINT) to be decomposed.
C
C     INCORE = 1: CHIA contains the full set of Cholesky vectors, and
C                 CHJB must be allocated to store the MDECOM rows of
C                 CHIA as specified in INDIA.
C
#include <implicit.h>
      DIMENSION XINT(*), FOCKD(*), CHIA(*), CHJB(*), WORK(LWORK)
      INTEGER INDIA(*)
#include <priunit.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_DCMMAT')

      IF ((IMAT.EQ.1) .OR. (IMAT.EQ.2)) THEN

C        (ia|jb) or (ai|bj) integrals (decided by which file
C        is connected to unit LUCHMO, or which is stored in
C        CHIA).
C        ---------------------------------------------------

         IF (INCORE .EQ. 1) THEN
            CALL CHO_AIBJ_I(XINT,CHIA,CHJB,INDIA,MDECOM,ISYM)
         ELSE
            CALL CHO_AIBJ(XINT,WORK,LWORK,INDIA,MDECOM,ISYM,LUCHMO)
         ENDIF

      ELSE IF (IMAT .EQ. 3) THEN

C        Minus the CC2 doubles amplitudes.
C        ---------------------------------

         IF (INCORE .EQ. 1) THEN
            CALL CC2_CHOAM_I(XINT,FOCKD,CHIA,CHJB,INDIA,MDECOM,ISYM)
         ELSE
            CALL CC2_CHOAM(XINT,FOCKD,WORK,LWORK,INDIA,MDECOM,ISYM,
     &                     LUCHMO)
         ENDIF

      ELSE

C        Error.
C        ------

         WRITE(LUPRI,'(//,5X,A,I10,A,A)')
     &   'IMAT = ',IMAT,' option not implemented in ',SECNAM
         CALL QUIT('Unknown option in '//SECNAM)

      ENDIF

      RETURN
      END
C  /* Deck cc_dcmanl */
      SUBROUTINE CC_DCMANL(DIAG,NDIAG,XLIST,NLIST)
C
C     Thomas Bondo Pedersen, September 2002.
C
C     Purpose: Print analysis of diagonal (histogram).
C
#include <implicit.h>
      DIMENSION DIAG(*), XLIST(*)
#include <priunit.h>

      PARAMETER (MXLIST = 20)
      INTEGER ICOUNT(MXLIST)

      PARAMETER (ZERO = 0.0D0)

C     Test input.
C     -----------

      IF (NDIAG .LE. 0) RETURN
      IF (NLIST .LE. 0) RETURN

C     Make sure that XLIST is in descending order.
C     --------------------------------------------

      IJOB = -1
      CALL ORDER3(XDUM,XLIST,IDUM,NLIST,0,IJOB)

C     Test that XLIST is positive.
C     ----------------------------

      IF (XLIST(1) .LE. ZERO) RETURN

C     Analyze.
C     --------

      NLTOT = MIN(NLIST,MXLIST)
      CALL IZERO(ICOUNT,NLTOT)
      NUMN  = 0
      NNEG  = 0
      XNEG  = ZERO

      DO IDIAG = 1,NDIAG

         TEST = DIAG(IDIAG)
         IF (TEST .LT. ZERO) THEN
            NNEG = NNEG + 1
            IF (TEST .LT. XNEG) XNEG = TEST
         ENDIF

         ILIST = 0
  100    CONTINUE

            ILIST = ILIST + 1
            IF (ILIST .GT. NLTOT) THEN
               NUMN = NUMN + 1
               GO TO 200
            ENDIF

            IF (TEST .GE. XLIST(ILIST)) THEN
               ICOUNT(ILIST) = ICOUNT(ILIST) + 1
               GO TO 200
            ENDIF

            GO TO 100

  200    CONTINUE

      ENDDO

C     Print.
C     ------

      WRITE(LUPRI,'(6X,A,11X,1P,D10.4,A,I10)')
     & 'Larger than ',XLIST(1),':',ICOUNT(1)
      DO ILIST = 2,NLTOT
         WRITE(LUPRI,'(6X,A,1P,D10.4,A,D10.4,A,I10)')
     &  'Between ',XLIST(ILIST-1),' and ',XLIST(ILIST),':',ICOUNT(ILIST)
      ENDDO
      WRITE(LUPRI,'(6X,A,10X,1P,D10.4,A,I10,/)')
     & 'Smaller than ',XLIST(NLTOT),':',NUMN

      IF (NNEG .GT. 0) THEN
         WRITE(LUPRI,'(6X,A,2X,I10,/,6X,A,1P,D12.4,/)')
     &   'Number of negative diagonals: ',NNEG,
     &   'Numerically largest         : ',XNEG
      ENDIF

      RETURN
      END
C  /* Deck cc_decmo_e */
      SUBROUTINE CC_DECMO_E(LUNIT)
#include <implicit.h>
      WRITE(LUNIT,'(5X,A,/)')
     & 'NOTICE: Original vectors held in core!'
      RETURN
      END
C  /* Deck cc_dcmcpai */
      SUBROUTINE CC_DCMCPAI(CHAI,CHBJ,LISTBJ,NUMBJ,ISYM)
C
C     Thomas Bondo Pedersen, September 2002.
C
C     Purpose: Copy sub-block CHBJ: L(J,bj) for bj in LISTBJ.
C
#include <implicit.h>
      DIMENSION CHAI(*), CHBJ(*)
      INTEGER LISTBJ(*)
#include <maxorb.h>
#include <ccdeco.h>
#include <ccsdsym.h>
#include <ccorb.h>

      INTEGER BJ

C      DO KBJ = 1,NUMBJ
C         BJ   = LISTBJ(KBJ)
C         KOFF = NUMCHO(ISYM)*(KBJ - 1) + 1
C         CALL DCOPY(NUMCHO(ISYM),CHAI(BJ),NT1AM(ISYM),
C     &                           CHBJ(KOFF),1)
C      ENDDO

      DO JVEC = 1,NUMCHO(ISYM)
         DO KBJ = 1,NUMBJ

            KOFF1 = NUMCHO(ISYM)*(KBJ - 1) + JVEC
            KOFF2 = NT1AM(ISYM)*(JVEC - 1) + LISTBJ(KBJ)

            CHBJ(KOFF1) = CHAI(KOFF2)

         ENDDO
      ENDDO

      RETURN
      END
C  /* Deck cc_dcmdxi */
      SUBROUTINE CC_DCMXDI(XDIANL,DIAG,CHOL,WORK,LWORK,NUMVEC,
     &                     NDIM,LUNIT,INCORE)
C
C     Thomas Bondo Pedersen, April 2003.
C
C     Purpose: "Poor man's check of MO decomposition"
C              On exit,
C              DIAG(i) <-- DIAG(i) - sum(J) L(i,J)**2, for i=1,NDIM and
C                                                          J=1,NUMVEC
C              XDIANL(1) = Min. diagonal in DIAG
C              XDIANL(2) = Max. diagonal in DIAG
C              XDIANL(3) = RMS error
C
C     If INCORE = 1: CHOL must contain the Cholesky vectors.
C
#include <implicit.h>
      DIMENSION XDIANL(3), DIAG(NDIM), CHOL(*), WORK(LWORK)
#include <priunit.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_DCMXDI')

C     Check that there's anything to do.
C     ==================================

      IF ((NDIM.LE.0) .OR. (NUMVEC.LE.0)) RETURN

C     Calculate DIAG.
C     ===============

      IF (INCORE .EQ. 1) THEN

C        Vectors available in CHOL array.
C        --------------------------------

         DO J = 1,NUMVEC
            KIJ = NDIM*(J - 1)
            DO I = 1,NDIM
               IJ = KIJ + I
               DIAG(I) = DIAG(I) - CHOL(IJ)*CHOL(IJ)
            ENDDO
         ENDDO

      ELSE

C        Vectors available on disk.
C        --------------------------

         NVEC = MIN(LWORK/NDIM,NUMVEC)
         IF (NVEC .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory for batch in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Minimum memory needed : ',NDIM,
     &      'Total memory available: ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF
         NBATCH = (NUMVEC - 1)/NVEC + 1

         DO IBATCH = 1,NBATCH

            NUMV = NVEC
            IF (IBATCH .EQ. NBATCH) THEN
               NUMV = NUMVEC - NVEC*(NBATCH - 1)
            ENDIF
            JVEC1 = NVEC*(IBATCH - 1) + 1

            CALL CHO_MOREAD(WORK,NDIM,NUMV,JVEC1,LUNIT)

            DO J = 1,NUMV
               KIJ = NDIM*(J - 1)
               DO I = 1,NDIM
                  IJ = KIJ + I
                  DIAG(I) = DIAG(I) - WORK(IJ)*WORK(IJ)
               ENDDO
            ENDDO

         ENDDO

      ENDIF

C     Calculate XDIANL: min, max, and RMS.
C     ------------------------------------

      XMIN = DIAG(1)
      XMAX = DIAG(1)
      DO I = 2,NDIM
         XMIN = MIN(XMIN,DIAG(I))
         XMAX = MAX(XMAX,DIAG(I))
      ENDDO

      XDIM = 1.0D0*NDIM
      RMS  = DDOT(NDIM,DIAG,1,DIAG,1)/XDIM

      XDIANL(1) = XMIN
      XDIANL(2) = XMAX
      XDIANL(3) = DSQRT(RMS)

      RETURN
      END
