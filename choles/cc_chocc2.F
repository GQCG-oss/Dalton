C  /* Deck cc_choecc2 */
      SUBROUTINE CC_CHOECC2(WORK,LWORK,ESCF,ECC2,RSPIM2)
C
C     Thomas Bondo Pedersen, August 2002.
C
C     Purpose:
C        Optimization of the CC2 wave function using Cholesky
C        decomposed two-electron integrals. The doubles amplitudes
C        are never stored in memory (unless memory actually allows
C        this).
C
C     Overall structure is based on CCSD_ENERGY by Henrik Koch.
C
C     IF (RSPIM2): Calculate the E intermediates for CCLR.
C
C     WARNING: NOCCIT has not been implemented: some essential intermediates
C              are assumed to be generated in the course of the calculation,
C              such as transformed Cholesky vectors and F(ia). These will be
C              needed if this is a response calculation. Thus, at least one
C              iteration must be carried out. (Unless, of course, you opt for
C              a rewriting of the code....)
C
#include <implicit.h>
      DIMENSION WORK(LWORK)
      LOGICAL RSPIM2
      COMMON /LUDIIS/ LUTDIS, LUSDIS
#include <maxorb.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <ccdeco.h>
#include <priunit.h>
#include <chocc2.h>
#include <cclr.h>
#include <dummy.h>
#include <ccfro.h>

      LOGICAL LCONVG,LXDIS

      CHARACTER*10 MODEL

      PARAMETER (LDIIST = 28, LDIISS = 29)

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CC_CHOECC2')

C     Start timing.
C     -------------

      CALL QENTER(SECNAM)
      CALL GETTIM(CSTR,WSTR)
      TIMTOT = SECOND()

C     Check that this is a Cholesky run.
C     ----------------------------------

      IF (.NOT. CHOINT) THEN
         WRITE(LUPRI,'(//,5X,A,A,A,/,5X,A,/)')
     &   'FATAL ERROR IN ',SECNAM,':',
     &   '- integrals must be Cholesky decomposed!'
         CALL QUIT('Fatal error in '//SECNAM)
      ENDIF

C     Set MODEL (for CC_WRRSP and CC_RDRSP).
C     --------------------------------------

      MODEL = 'CC2       '

C     Open DIIS file (CRAY only).
C     ---------------------------

#if defined (SYS_CRAY)
      CALL WOPEN('CC_DIIS',64,0,IERR)
      IF (IERR .NE. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A,/,5X,A,I10,/)')
     &   'Error opening file CC_DIIS in ',SECNAM,':',
     &   'Subroutine WOPEN returned IERR = ',IERR
         CALL QUIT('Error opening CC_DIIS in '//SECNAM)
      ENDIF
#endif

C     Print header.
C     -------------

      CALL HEADER('Output from '//SECNAM,-1)
      WRITE(LUPRI,'(5X,A)')
     & 'Cholesky based optimization of the CC2 wave function:'
      IF (CHOT2) THEN
         WRITE(LUPRI,'(5X,A)')
     &   'The doubles amplitudes will be separately decomposed.'
      ELSE IF (CHOMO2) THEN
         WRITE(LUPRI,'(5X,A)')
     &   'The (ai|bj) integral matrix will be separately decomposed.'
      ENDIF
      WRITE(LUPRI,'(5X,A,1P,D15.6,/,5X,A,1P,D15.6)')
     & 'Threshold for energy         : ',THRENR,
     & 'Threshold for vector function: ',THRVEC
      WRITE(LUPRI,'(5X,A,5X,I10,/)')
     & 'Maximum number of iterations : ',MAXITE

      CALL FLSHFO(LUPRI)

C     Calculate T1-independent MO transformations.
C     --------------------------------------------

      CALL CC2_TRF0(WORK,LWORK)

C     Allocation.
C     -----------

      KFOCKD = 1
      KT1AM  = KFOCKD + NORBTS
      KOMEG1 = KT1AM  + NT1AMX
      KEND1  = KOMEG1 + NT1AMX
      LWRK1  = LWORK  - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: T1'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need (more than): ',KEND1-1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Get SCF and orbital energies.
C     -----------------------------

      CALL CC2_GETSI(WORK(KFOCKD),WORK(KEND1),LWRK1,ESCF,POTNUC)

C     Set up initial amplitudes.
C     --------------------------

      CALL CC2_INITAM(WORK(KFOCKD),WORK(KT1AM),WORK(KOMEG1),
     &                WORK(KEND1),LWRK1,ECC2,IRST)


C     Start of iterative loop.
C     ------------------------

      ITER = 0
      EN1  = ECC2

      IF (IRST .EQ. 1) THEN
         WRITE(LUPRI,'(1X,A,I3,A,F23.16)')
     &   'Iter.',ITER,': Coupled cluster RSTAR energy : ',ECC2
         CALL FLSHFO(LUPRI)
      ENDIF

  100 CONTINUE

C        Update iteration counter.
C        -------------------------

         ITER = ITER + 1

         IF (ITER .GT. MAXITE) THEN
            WRITE(LUPRI,'(//,5X,A,A,I4,A,/)')
     &      SECNAM,': CC2 wave function not converged in ',MAXITE,
     &      ' iterations.'
            CALL QUIT(SECNAM//': CC2 wave function not converged')
         ENDIF

         IF (IPRINT .GT. 2) THEN
            WRITE(LUPRI,'(/,7X,A,I3)') 'Iteration no.:',ITER
            WRITE(LUPRI,'(7X,A,/)')    '-----------------'
         ENDIF

C        Calculate the CC2 vector function and energy.
C        ---------------------------------------------

c      lwsav = lwrk1
c      lwrk1 = 10000
c      write(lupri,*) 'calling rhs with lwork : ',LWRK1

         CALL CC2_CHORHS(WORK(KFOCKD),WORK(KT1AM),WORK(KOMEG1),
     &                   WORK(KEND1),LWRK1,T2NM,ECC2,ESCF)
         EN2 = ECC2

c      lwrk1 = lwsav

         IF (IPRINT .GE. 5) THEN
            T1NM = DDOT(NT1AMX,WORK(KT1AM),1,WORK(KT1AM),1)
            WRITE(LUPRI,'(7X,A,1P,D25.10)')
     &      'Norm of t1am   after ccvec: ',DSQRT(T1NM)
            WRITE(LUPRI,'(7X,A,1P,D25.10)')
     &      'Norm of t2am   after ccvec: ',DSQRT(T2NM)
            WRITE(LUPRI,'(7X,A,1P,D25.10)')
     &      'Total norm of amplitudes  : ',DSQRT(T1NM+T2NM)
         ENDIF

         O1NM = DSQRT(DDOT(NT1AMX,WORK(KOMEG1),1,WORK(KOMEG1),1))
         IF (IPRINT .GE. 3) THEN
            WRITE(LUPRI,'(7X,A,1P,D25.10)')
     &      'Norm of omega1 after ccvec: ',O1NM
            IF (IPRINT .GE. 25) THEN
               CALL HEADER('CC2 Vector Function',-1)
               DO ISYM = 1,NSYM
                  WRITE(LUPRI,'(10X,A,I1,A,/)')
     &            'Symmetry block ',ISYM,':'
                  NA = NVIR(ISYM)
                  NI = NRHF(ISYM)
                  IF ((NA.GT.0) .AND. (NI.GT.0)) THEN
                     KOFF = KOMEG1 + IT1AM(ISYM,ISYM)
                     CALL OUTPUT(WORK(KOFF),1,NA,1,NI,NA,NI,1,LUPRI)
                     WRITE(LUPRI,*)
                  ELSE
                     WRITE(LUPRI,'(10X,A,/)')
     &               ' - empty block'
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF ((ITER.EQ.1) .AND. (IRST.EQ.0)) THEN
            WRITE(LUPRI,'(1X,A,I3,A,F23.16)')
     &      'Iter.',ITER,': Coupled cluster MP2   energy : ',ECC2
         ELSE
            WRITE(LUPRI,'(1X,A,I3,A,F23.16)')
     &      'Iter.',ITER,': Coupled cluster CC2   energy : ',ECC2
         ENDIF
         IF (IPRINT .GE. 2) THEN
            WRITE(LUPRI,'(33X,A,9X,1P,D14.6)')
     &      'DeltaE : ',EN2-EN1
         END IF

         CALL FLSHFO(LUPRI)

C        Save information for this iteration.
C        ------------------------------------

         CALL CC2_SAVTAM(WORK(KT1AM),WORK(KOMEG1),ECC2)


C        Check convergence.
C        ------------------

         DELE   = DABS(EN2-EN1)
         LCONVG = (DELE.LE.THRENR) .AND. (O1NM.LE.THRVEC)
         IF (.NOT. LCONVG) THEN
            CALL CCSD_NXTAM(WORK(KT1AM),DUM1,DUM12,WORK(KOMEG1),
     &                      DUM2,DUM22,WORK(KFOCKD),.FALSE.,1,0.0D0)
            CALL CCSD_DIIS(WORK(KOMEG1),WORK(KT1AM),NT1AMX,ITER)
            EN1 = EN2
            GOTO 100
         ENDIF

C     CC2 wave function has converged.
C     --------------------------------

      WRITE(LUPRI,'(/,1X,A,D10.4,A,F16.9)')
     & 'CC2 energy converged to within ',THRENR,' is ',ECC2
      WRITE(LUPRI,'(1X,A,F16.9)')
     & '- the CC2 correlation energy contribution is ',ECC2-ESCF
      WRITE(LUPRI,'(1X,A,6X,1P,D14.4)')
     & 'Final norm of the CC2 vector function is ',O1NM 
      TIMTOT = SECOND() - TIMTOT
      WRITE(LUPRI,'(1X,A,7X,F14.2,A,/)')
     & 'Total time used for CC2 optimization is ',TIMTOT,' seconds'

      CALL FLSHFO(LUPRI)

C     Calculate and write Fock matrix
C     -------------------------------

      IF (RSPIM2) THEN
         CALL CHO_FCKH(WORK(KT1AM),WORK(KEND1),LWRK1)
      END IF

C     Close DIIS files.
C     -----------------

#if defined (SYS_CRAY)
      CALL WCLOSE('CC_DIIS',IERR)
      INFO = ISHELL('rm CC_DIIS')
#endif
#if !defined (SYS_CRAY)
      IF (LUTDIS .GT. 0) THEN
         INQUIRE(UNIT=LUTDIS,OPENED=LXDIS,ERR=121)
         IF (LXDIS) CALL GPCLOSE(LUTDIS,'DELETE')
      ENDIF
      IF (LUSDIS .GT. 0) THEN
         INQUIRE(UNIT=LUSDIS,OPENED=LXDIS,ERR=121)
         CALL GPCLOSE(LUSDIS,'DELETE')
      ENDIF
  121 CONTINUE
#endif

C     Print largest components of 0th order wave function.
C     ----------------------------------------------------

      IF (IPRINT .GT. 2) THEN
         CALL AROUND('Largest amplitudes in converged solution')
         CALL CC_PRAM(WORK(KT1AM),PT1,1)
      ENDIF

C     Save a copy on rsp-file system.
C     -------------------------------

      KT0AM = KT1AM + NT1AMX
      KEND1 = KT0AM + 2*NALLAI(1)
      LWRK1 = LWORK - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   'Insufficient memory for amplitude saving in ',SECNAM
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need     : ',KEND1-1,
     &   'Available: ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

      CALL DZERO(WORK(KT0AM),2*NALLAI(1))

      IOPT = 7
      CALL CC_WRRSP('R0',0,1,IOPT,MODEL,WORK(KT0AM),WORK(KT1AM),
     &              DUMMY,WORK(KEND1),LWRK1)

C     Calculate global E intermediates for response.
C     ----------------------------------------------

      IF (RSPIM2.AND.(.NOT.IMSKIP)) THEN

         RSPIM = RSPIM2

         WRITE(LUPRI,'(/)')
         CALL AROUND('Calculating Intermediates for CCLR')
         WRITE(LUPRI,'(/)')

         KEND1 = KT1AM + NT1AMX
         LWRK1 = LWORK - KEND1 + 1
         CALL CC_CHOEIM(WORK(KFOCKD),WORK(KT1AM),WORK(KEND1),LWRK1,
     &                  .TRUE.)

         IF (IPRINT .GT. 1) WRITE(LUPRI,'(/)')
         WRITE(LUPRI,'(12X,A)') 'E-intermediates calculated'

      ELSE IF (RSPIM2.AND.IMSKIP) THEN

         RSPIM = RSPIM2
         WRITE(LUPRI,'(12X,A)')
     &        'Intermediates assumed to be restart IM. '

      ENDIF

C     Delete files with Y intermediates.
C     (Might have been done already in cc_choeim).
C     --------------------------------------------

      DO ISYCHO = 1,NSYM
         CALL CC_CYIOP(-1,ISYCHO,0)
         CALL CC_CYIOP(0,ISYCHO,0)
      ENDDO

      CALL GETTIM(CEND,WEND)
      CTOT = CEND - CSTR
      WTOT = WEND - WSTR
      CALL TIMTXT('>>>> Total CPU  time used in '//SECNAM//':',CTOT,
     &            LUPRI)
      CALL TIMTXT('>>>> Total wall time used in '//SECNAM//':',WTOT,
     &            LUPRI)
      CALL QEXIT(SECNAM)
      RETURN
      END
C  /* Deck cc2_getsi */
      SUBROUTINE CC2_GETSI(FOCKD,WORK,LWORK,ESCF,POTNUC)
C
C     Thomas Bondo Pedersen, July 2002.
C
C     Purpose: Read Fock diagonal and SCF energy from SIRIUS
C              interface file.
C
#include <implicit.h>
      DIMENSION FOCKD(*), WORK(LWORK)
#include <dummy.h>

      CALL CHO_RDSIR(POTNUC,ESCF,FOCKD,DUMMY,WORK,LWORK,.TRUE.,.TRUE.,
     &               .FALSE.)

      RETURN
      END
C  /* Deck cc2_initam */
      SUBROUTINE CC2_INITAM(FOCKD,T1AM,OMEGA1,WORK,LWORK,ECC2,IRST)
C
C     Thomas Bondo Pedersen, August 2002.
C
C     Purpose: Set up initial Cholesky CC2 amplitudes.
C              Restart if requested and if possible.
C
C     On exit:
C
C        IRST = 0: no restart [if requested, it failed]
C        IRST = 1: "full" restart: T1AM, OMEGA1 and ECC2 read from disk.
C        IRST = 2: restart only with T1AM [as in conventional code].
C
C     The reason for trying to read not only T1AM but also OMEGA1 and ECC2
C     is that with conventional amplitude restart, the Cholesky CC2 energy
C     solver will spend 2 iterations before realizing that the amplitudes
C     are converged: the first iteration to get energy and vector function,
C     and convergence will be found after the second (if the initial ones
C     were in fact converged, obviously).
C
#include <implicit.h>
      DIMENSION FOCKD(*), T1AM(*), OMEGA1(*), WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>
#include <dummy.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CC2_INITAM')

      LOGICAL RSTART

      CHARACTER*10 FILRST, MODEL
      PARAMETER (FILRST = 'CHOCC2.RST')

      PARAMETER (ZERO = 0.0D0)

C     Set model for reading in CC_RDRSP.
C     ----------------------------------

      IF (CC2) THEN
         MODEL = 'CC2       '
      ELSE
         WRITE(LUPRI,'(//,5X,A,A)')
     &   SECNAM,': Error: only CC2 model recognized!'
         CALL QUIT('Error in '//SECNAM)
      ENDIF

C     Default is MP2 guess.
C     ---------------------

      IRST = 0

C     If requested, try to restart.
C     -----------------------------

      IF (CCRSTR) THEN

C        First try full restart.
C        -----------------------

         INQUIRE(FILE=FILRST,EXIST=RSTART,ERR=200)
         IF (RSTART) THEN
            IRST = 1
         ELSE
            IRST = 2
         ENDIF
         GO TO 100

  200    CONTINUE
         IF (IPRINT .GE. 1) THEN
            WRITE(LUPRI,'(5X,A,A)')
     &      SECNAM,': an error ocurred inquering file ',FILRST
            WRITE(LUPRI,'(5X,A,A)')
     &      SECNAM,': trying to restart from amplitudes only.'
         ENDIF
         IRST = 2

  100    CONTINUE

         IF (IRST .EQ. 1) THEN

C           FILRST exists, try to read data while checking that the
C           amount of data [NT1AMX] is consistent with current run.
C           If not, reset IRST to 2.
C           -------------------------------------------------------

            IFAIL = 0
            CALL CHO_RDRST(ECC2,T1AM,OMEGA1,.TRUE.,.TRUE.,.TRUE.,IFAIL)
            IF (IFAIL .EQ. 0) THEN
               IF (IPRINT .GE. 1) THEN
                  WRITE(LUPRI,'(5X,A,A,A,A)')
     &            SECNAM,': file ',FILRST,
     &            ' detected and accepted for full restart.'
               ENDIF
            ELSE
               IF (IPRINT .GE. 1) THEN
                  WRITE(LUPRI,'(5X,A,A,A,A)')
     &            SECNAM,': file ',FILRST,
     &            ' detected but rejected for restart.'
               ENDIF
               IRST = 2
            ENDIF

         ENDIF

      ENDIF

C     Action taken according to restart (and which) or not.
C     =====================================================

      IF (IRST .EQ. 1) THEN

C        T1AM and OMEGA1 (and ECC2) were successfully read from FILRST.
C        Generate new amplitudes by perturbation theory.
C        --------------------------------------------------------------

         ISYALF = 1
         CALL CCSD_NXTAM(T1AM,DUMMY,DUM12,OMEGA1,DUMMY,DUM22,FOCKD,
     &                   .FALSE.,ISYALF,0.0D0)
         CALL DCOPY(NT1AMX,OMEGA1,1,T1AM,1)

      ELSE IF (IRST .EQ. 2) THEN

C        For whatever reason, restart from FILRST failed. Now try
C        to get amplitudes from 'R0' file or, if that fails, the
C        original ("old restart") CCSD_TAM file.
C        --------------------------------------------------------

         ECC2 = ZERO

         ILST = 1
         ISYM = 1
         IOPT = 33
         CALL CC_RDRSP('R0 ',ILST,ISYM,IOPT,MODEL,T1AM,DUMMY)

         IF (IOPT .EQ. 33) THEN
            INQUIRE(FILE='CCSD_TAM',EXIST=RSTART,ERR=1000)
            GO TO 1100
 1000       RSTART = .FALSE.
 1100       CONTINUE
            IF (RSTART) THEN
              LUTAM = -1
              CALL GPOPEN(LUTAM,'CCSD_TAM','UNKNOWN',' ','UNFORMATTED',
     *                    IDUMMY,.FALSE.)
              REWIND(LUTAM)
              READ(LUTAM) (T1AM(I), I = 1,NT1AMX)
              CALL GPCLOSE(LUTAM,'KEEP')
              IF (IPRINT .GE. 1) THEN
                 WRITE(LUPRI,'(5X,A,A)')
     &           SECNAM,': "old" amplitude restart succeeded.'
              ENDIF
            ELSE
               IF (IPRINT .GE. 1) THEN
                  WRITE(LUPRI,'(5X,A,A)')
     &            SECNAM,': amplitude restart failed.'
                  WRITE(LUPRI,'(5X,A,A)')
     &            SECNAM,': setting up MP2 initial guess instead.'
               ENDIF
               CALL DZERO(T1AM,NT1AMX)
               IRST = 0
            ENDIF
         ELSE
            IF (IPRINT .GE. 1) THEN
               WRITE(LUPRI,'(5X,A,A)')
     &         SECNAM,': amplitude restart succeeded.'
            ENDIF
         ENDIF

      ELSE IF (IRST .EQ. 0) THEN

C        MP2 initial guess.
C        ------------------

         IF (IPRINT .GT. 10) THEN
            WRITE(LUPRI,'(5X,A,A)')
     &      SECNAM,': Setting up MP2 initial guess.'
         ENDIF

         ECC2 = ZERO
         CALL DZERO(T1AM,NT1AMX)

      ELSE

C        Logical error.
C        --------------

         WRITE(LUPRI,'(//,5X,A,A)')
     &   'Logical error ocurred in ',SECNAM
         WRITE(LUPRI,'(5X,A,I10,/)')
     &   'IRST out of bounds: ',IRST
         CALL QUIT('Logical error in '//SECNAM)

      ENDIF

C     Print.
C     ------

      IF (CCRSTR .AND. (IRST.EQ.0)) THEN
         WRITE(LUPRI,'(5X,A,A)')
     &   SECNAM,': Restart failed. Initial amplitudes are MP2 guess'
      ENDIF

      IF ((IPRINT.GT.100) .OR. DEBUG) THEN
         T1NRM = DSQRT(DDOT(NT1AMX,T1AM,1,T1AM,1))
         CALL AROUND('Initial CC2 amplitudes in '//SECNAM)
         WRITE(LUPRI,'(10X,A,I2)')
     &   'IRST =',IRST
         WRITE(LUPRI,'(10X,A,1P,D22.15,/)')
     &   'Norm of amplitudes: ',T1NRM
         CALL NOCC_PRT(T1AM,1,'AI  ')
      ENDIF

      RETURN
      END
C  /* Deck cc2_savtam */
      SUBROUTINE CC2_SAVTAM(T1AM,OMEGA1,ECC2)
C
C     Thomas Bondo Pedersen, August 2002.
C
C     Purpose: Save CC2 amplitudes, vector function, and energy
C              on disk.
C
#include <implicit.h>
      DIMENSION T1AM(*), OMEGA1(*)

      CALL CHO_WRRST(ECC2,T1AM,OMEGA1)

      RETURN
      END
C  /* Deck cc2_trf0 */
      SUBROUTINE CC2_TRF0(WORK,LWORK)
C
C     Thomas Bondo Pedersen, July 2002.
C
C     Purpose: T1-independent MO transformations of Cholesky
C              vectors; results are written to disk.
C
#include <implicit.h>
      DIMENSION WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      CHARACTER*8 SECNAM
      PARAMETER (SECNAM = 'CC2_TRF0')

C     Allocation.
C     -----------

      KCMO  = 1
      KEND1 = KCMO  + NLAMDS
      LWRK1 = LWORK - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: 1'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need (more than): ',KEND1-1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Read the MO coefficient matrix. Change from SIRIUS to CC ordering.
C     ------------------------------------------------------------------

      CALL CHO_RDSIR(DUM1,DUM2,DUM3,WORK(KCMO),WORK(KEND1),LWRK1,
     &               .FALSE.,.FALSE.,.TRUE.)

C     Transform Cholesky vectors: L(ia,J).
C     Result vectors are stored on disk.
C     ------------------------------------

      CALL CHO_TRFAI(WORK(KCMO),WORK(KEND1),LWRK1)

C     Transform one-electron integrals: h(ia).
C     Result is stored on disk.
C     ----------------------------------------
Ctbp: obsolete in new version
Ctbp  CALL CC2_ONETRF(WORK(KCMO),WORK(KEND1),LWRK1)
Ctbp
      RETURN
      END
C  /* Deck cc2_onetrf */
      SUBROUTINE CC2_ONETRF(CMO,WORK,LWORK)
C
C     Thomas Bondo Pedersen, July 2002.
C
C     Purpose: Partially transform AO one-electron integrals
C              to MO basis. Result array h(ia) stored on disk.
C
#include <implicit.h>
      DIMENSION CMO(*), WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CC2_ONETRF')

      PARAMETER (ITYP = 1, IOPEN = -1, IKEEP = 1)
      PARAMETER (ZERO = 0.00D0, ONE = 1.00D0)

C     Allocation.
C     -----------

      MAXALI = 0
      DO ISYM = 1,NSYM
         MAXALI = MAX(MAXALI,NBAS(ISYM)*NRHF(ISYM))
      ENDDO

      KONEL = 1
      KHIA  = KONEL + N2BST(1)
      KSCR1 = KHIA  + NT1AM(1)
      KEND1 = KSCR1 + MAXALI
      LWRK1 = LWORK - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: ONEL'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need     : ',KEND1-1,
     &   'Available: ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Read AO integrals.
C     ------------------

      CALL CCRHS_ONEAO(WORK(KONEL),WORK(KEND1),LWRK1)

C     Transform.
C     ----------

      DO ISYMI = 1,NSYM

         ISYMG = ISYMI
         ISYMD = ISYMG
         ISYMA = ISYMD

         NTOTG = MAX(NBAS(ISYMG),1)
         NTOTD = MAX(NBAS(ISYMD),1)
         NTOTA = MAX(NVIR(ISYMA),1)

         KOFF1 = KONEL + IAODIS(ISYMG,ISYMD)
         KOFF2 = IGLMRH(ISYMG,ISYMI) + 1
         KOFF3 = IGLMVI(ISYMD,ISYMA) + 1
         KOFF4 = KHIA + IT1AM(ISYMA,ISYMI)

         CALL DGEMM('T','N',NBAS(ISYMD),NRHF(ISYMI),NBAS(ISYMG),
     &              ONE,WORK(KOFF1),NTOTG,CMO(KOFF2),NTOTG,
     &              ZERO,WORK(KSCR1),NTOTD)

         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMD),
     &              ONE,CMO(KOFF3),NTOTD,WORK(KSCR1),NTOTD,
     &              ZERO,WORK(KOFF4),NTOTA)

      ENDDO

C     Write h(ia) to disk.
C     --------------------

      CALL ONEL_OP(IOPEN,ITYP,LUHIA)
      CALL CHO_MOWRITE(WORK(KHIA),NT1AM(1),1,1,LUHIA)
      CALL ONEL_OP(IKEEP,ITYP,LUHIA)

      RETURN
      END
C  /* Deck cc2_chorhs */
      SUBROUTINE CC2_CHORHS(FOCKD,T1AM,OMEGA1,WORK,LWORK,T2NRM,ECC2,
     &                      ESCF)
C
C     Thomas Bondo Pedersen, August 2002.
C
C     Purpose: Calculate the CC2 vector function from Cholesky
C              decomposed two-electron integrals. The CC2 energy
C              ECC2 is calculated as a byproduct; note that it
C              contains the SCF energy.
C
C     Notes:
C        (1) The T1-independent MO transformations are assumed
C            available on disk, i.e. L(ia,J).
C        (2) The T2 amplitudes may be treated in 3 different ways:
C            - by straightforward calculation from (ai|bj) integrals
C              represented by the transformed AO Cholesky vectors.
C            - by calculation from separately Cholesky decomposed
C              (ai|bj) integrals.
C            - by decomposition of the amplitudes themselves.
C
#include <implicit.h>
      DIMENSION FOCKD(*), T1AM(*), OMEGA1(*), WORK(LWORK)
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CC2_CHORHS')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      PARAMETER (INFTOT = 2, INFO = 3)
      PARAMETER (ZERO = 0.0D0)

C     Start timing this routine.
C     --------------------------

      TIMTOT = SECOND()

      IF (LOCDBG) THEN
         OESUM = ZERO
         DO ISYMI = 1,NSYM
            DO I = 1,NRHF(ISYMI)
               KOFFI = IRHF(ISYMI) + I
               OESUM = OESUM + FOCKD(KOFFI)
            ENDDO
         ENDDO
         WRITE(LUPRI,'(5X,A,A,F16.9)')
     &   SECNAM,': debug: Sum of occupied orbital energies: ',OESUM
      ENDIF

C     Initialization.
C     ---------------

      ECC2  = ESCF
      T2NRM = ZERO

C     Allocation.
C     -----------

      KTRACE = 1
      KLAMDP = KTRACE + NUMCHO(1)
      KLAMDH = KLAMDP + NGLMDT(1)
      KEND1  = KLAMDH + NGLMDT(1)
      LWRK1  = LWORK  - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: Lambda'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need (more than): ',KEND1-1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Calculate Lambda matrices.
C     --------------------------

      TIMLAM = SECOND()
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),T1AM,WORK(KEND1),LWRK1)
      TIMLAM = SECOND() - TIMLAM

C     Calculate T1-dependent MO transformations.
C     Calculate T1-dependent energy contribution.
C     -------------------------------------------

      TIMTRL = SECOND()
      ENERGY = ZERO
      CALL CC2_TRFL(T1AM,WORK(KLAMDP),WORK(KLAMDH),WORK(KTRACE),
     &              WORK(KEND1),LWRK1,ENERGY)
      ECC2   = ECC2 + ENERGY
      TIMTRL = SECOND() - TIMTRL

      IF (LOCDBG) THEN
         WRITE(LUPRI,'(5X,A,A,F16.9,/,5X,A,A,F16.9)')
     &   SECNAM,': debug: T1 contribution to energy: ',ENERGY,
     &   SECNAM,': debug: Accumulated CC2 energy   : ',ECC2
      ENDIF

C     Initialize result vector.
C     -------------------------

      CALL CC2_INIOME(FOCKD,T1AM,OMEGA1,1)

C     Calculate Y-intermediates and the I-term.
C     The Y-intermediates are stored on disk,
C     and the I-term is added into OMEGA1.
C     Calculate T2-dependent energy contribution.
C     -------------------------------------------

      TIMYI  = SECOND()
      ENERGY = ZERO
      CALL CHOCC2_YI(FOCKD,OMEGA1,WORK(KEND1),LWRK1,T2NRM,ENERGY)
      ECC2  = ECC2 + ENERGY
      TIMYI = SECOND() - TIMYI

      IF (LOCDBG) THEN
         OMNRM = DSQRT(DDOT(NT1AMX,OMEGA1,1,OMEGA1,1))
         WRITE(LUPRI,'(5X,A,A,F16.9,/,5X,A,A,F16.9)')
     &   SECNAM,': debug: T2 contribution to energy: ',ENERGY,
     &   SECNAM,': debug: Accumulated CC2 energy   : ',ECC2
         WRITE(LUPRI,'(5X,A,A,1P,D15.6)')
     &   SECNAM,': debug: Norm of Omega1 after CHOCC2_YI: ',OMNRM
      ENDIF

C     Calculate the JG-term.
C     N.B.: TRACE is scaled by 2 on exit from CHOCC2_JG.
C     --------------------------------------------------

      TIMJG = SECOND()
      CALL CHOCC2_JG(T1AM,OMEGA1,WORK(KLAMDP),WORK(KLAMDH),
     &               WORK(KTRACE),WORK(KEND1),LWRK1)
      TIMJG = SECOND() - TIMJG

      IF (LOCDBG) THEN
         OMNRM = DSQRT(DDOT(NT1AMX,OMEGA1,1,OMEGA1,1))
         WRITE(LUPRI,'(5X,A,A,1P,D15.6)')
     &   SECNAM,': debug: Norm of Omega1 after CHOCC2_JG: ',OMNRM
      ENDIF

C     Calculate the H term.
C     Pass full memory: trace and Lambda matrices no 
C     longer needed.
C     ----------------------------------------------

      TIMJH = SECOND()
      CALL CHOCC2_H(OMEGA1,WORK,LWORK)
      TIMJH = SECOND() - TIMJH

      IF (LOCDBG) THEN
         OMNRM = DSQRT(DDOT(NT1AMX,OMEGA1,1,OMEGA1,1))
         WRITE(LUPRI,'(5X,A,A,1P,D15.6)')
     &   SECNAM,': debug: Norm of Omega1 after CHOCC2_H: ',OMNRM
      ENDIF

C     Print section.
C     --------------

      IF (IPRINT .GT. INFTOT) THEN
         TIMTOT = SECOND() - TIMTOT
         WRITE(LUPRI,'(7X,A,F17.2,A)')
     &   'Time used in RHS - TOTAL  : ',TIMTOT,' seconds'
         IF (IPRINT .GE. INFO) THEN
            WRITE(LUPRI,'(7X,A,F17.2,A)')
     &      'Time used in LAMMAT       : ',TIMLAM,' seconds'
            WRITE(LUPRI,'(7X,A,F17.2,A)')
     &      'Time used in CC2_TRFL     : ',TIMTRL,' seconds'
            WRITE(LUPRI,'(7X,A,F17.2,A)')
     &      'Time used in CHOCC2_YI    : ',TIMYI,' seconds'
            WRITE(LUPRI,'(7X,A,F17.2,A)')
     &      'Time used in CHOCC2_JG    : ',TIMJG,' seconds'
            WRITE(LUPRI,'(7X,A,F17.2,A)')
     &      'Time used in CHOCC2_H     : ',TIMJH,' seconds'
         ENDIF
      ENDIF

      RETURN
      END
C  /* Deck cc2_trfl */
      SUBROUTINE CC2_TRFL(T1AM,XLAMDP,XLAMDH,TRACE,WORK,LWORK,ECC2)
C
C     Thomas Bondo Pedersen, August 2002.
C
C     Purpose: Calculate T1-dependent MO transformations for
C              Cholesky CC2.
C              The results are written to disk:
C              L(ij,J), L(ai,J), F(ia).
C              Calculate T1-dependent energy contribution.
C              On exit, the TRACE array contains the vector
C              TRACE(J) = sum(ai) L(ia,J) * T1AM(ai).
C
#include <implicit.h>
      DIMENSION T1AM(*), XLAMDP(*), XLAMDH(*), TRACE(*), WORK(LWORK)
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      CHARACTER*8 SECNAM
      PARAMETER (SECNAM = 'CC2_TRFL')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      PARAMETER (IOPEN = -1, IKEEP = 1)
      PARAMETER (ZERO = 0.00D0, ONE = 1.00D0, TWO = 2.00D0)

C     Calculate F(ia), L(ij,J), and L(ai,J).
C     Calculate T1-dependent energy contributions.
C     --------------------------------------------

      KFOCK = 1
      KEND1 = KFOCK + NT1AMX
      LWRK1 = LWORK - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: Fock'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need     : ',KEND1-1,
     &   'Available: ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

      CALL DZERO(WORK(KFOCK),NT1AMX)

      CALL CC2_TRF1(T1AM,XLAMDP,XLAMDH,TRACE,WORK(KFOCK),
     &              WORK(KEND1),LWRK1,ECC2)

      IF (LOCDBG) THEN
         FNORM = DSQRT(DDOT(NT1AMX,WORK(KFOCK),1,WORK(KFOCK),1))
         WRITE(LUPRI,'(5X,A,A,1P,D15.6)')
     &   SECNAM,': debug: Norm of F(ia) after CC2_TRF1: ',FNORM
      ENDIF

      CALL ONEL_OP(IOPEN,3,LUFOCK)
      CALL CHO_MOWRITE(WORK(KFOCK),NT1AMX,1,1,LUFOCK)
      CALL ONEL_OP(IKEEP,3,LUFOCK)

      RETURN
      END
C  /* Deck cc2_trf1 */
      SUBROUTINE CC2_TRF1(T1AM,XLAMDP,XLAMDH,TRACE,FOCK,WORK,LWORK,ECC2)
C
C     Thomas Bondo Pedersen, August 2002.
C
C     Purpose: Calculate T1-dependent MO transformations
C              F(ia), L(ji,J), and L(ai,J).
C              The Fock matrix must be initialized.
C              Calculate energy contribution from singles:
C
C              ECC2 = ECC2 + 2 * sum(J) [sum(ai) T1AM(ai) * L(ia,J)]**2
C
C                          -   sum(jiJ) [sum(a)  T1AM(aj) * L(ia,J)]
C
C                                     * [sum(b)  T1AM(bi) * L(jb,J)]
C
C              On exit, the TRACE array contains the vector
C
C              TRACE(J) = sum(ai) L(ia,J) * T1AM(ai).
C
#include <implicit.h>
      DIMENSION T1AM(*), XLAMDP(*), XLAMDH(*), TRACE(*), FOCK(*)
      DIMENSION WORK(LWORK)
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <chocc2.h>
#include <priunit.h>

      CHARACTER*8 SECNAM
      PARAMETER (SECNAM = 'CC2_TRF1')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      PARAMETER (IOPEN = -1, IKEEP = 1, IOPTR = 2)
      PARAMETER (ITYIA = 1, ITYAI = 3, ITYJI = 4)
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)

C     Initialize energy contributions.
C     --------------------------------

      ENRC = ZERO
      ENRX = ZERO

C     Read reduce index array.
C     ------------------------

      KIND1 = 1
      CALL CC_GETIND1(WORK(KIND1),LWORK,LIND1)
      KEND0 = KIND1 + LIND1
      LWRK0 = LWORK - KEND0 + 1

      IF (LWRK0 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: index'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need (more than): ',KEND0-1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Start Cholesky symmetry loop.
C     -----------------------------

      DO ISYCHO = 1,NSYM

         IF (NUMCHO(ISYCHO) .LE. 0) GO TO 999
         IF (N2BST(ISYCHO)  .LE. 0) GO TO 999

C        Open MO Cholesky files for this symmetry.
C        -----------------------------------------

         CALL CHO_MOP(IOPEN,ITYJI,ISYCHO,LUJI,1,1)
         CALL CHO_MOP(IOPEN,ITYAI,ISYCHO,LUAI,1,1)
         CALL CHO_MOP(IOPEN,ITYIA,ISYCHO,LUIA,1,1)

C        Calculate size of scratch space needed for read
C        and first half-transformation [L(gamma i)].
C        -----------------------------------------------

         MAXGI = 0
         DO ISYMI = 1,NSYM
            ISYMG = MULD2H(ISYMI,ISYCHO)
            MAXGI = MAX(MAXGI,NBAS(ISYMG)*NRHF(ISYMI))
         ENDDO
         LREAD = MEMRD(1,ISYCHO,IOPTR)
         LSCR1 = MAX(LREAD,MAXGI)

C        Allocation.
C        -----------

         KCHOL = KEND0
         KSCR1 = KCHOL + N2BST(ISYCHO)
         KEND1 = KSCR1 + LSCR1
         LWRK1 = LWORK - KEND1 + 1

         IF (LWRK1 .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Insufficient memory in ',SECNAM,' - allocation: Cholesky'
            MREQ = KEND1 + NT1AM(ISYCHO) + NMATIJ(ISYCHO) - 1
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need (minimum): ',MREQ,
     &      'Available     : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Set up batch.
C        -------------

         MINMEM = NT1AM(ISYCHO) + NMATIJ(ISYCHO)
         NVEC   = MIN(LWRK1/MINMEM,NUMCHO(ISYCHO))

         IF (NVEC .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory for batch in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,I10,/,5X,A,I10,/)')
     &      'Minimum memory required   : ',MINMEM,
     &      'Memory available for batch: ',LWRK1,
     &      'Memory available in total : ',LWORK,
     &      'Number of Cholesky vectors: ',NUMCHO(ISYCHO)
            CALL QUIT('Insufficient memory for batch in '//SECNAM)
         ENDIF

         NBATCH = (NUMCHO(ISYCHO) - 1)/NVEC + 1

C        Batch loop.
C        -----------

         DO IBATCH = 1,NBATCH

            NUMV = NVEC
            IF (IBATCH .EQ. NBATCH) THEN
               NUMV = NUMCHO(ISYCHO) - NVEC*(NBATCH - 1)
            ENDIF

            JVEC1 = NVEC*(IBATCH - 1) + 1

C           Complete allocation.
C           --------------------

            KCHJI = KEND1
            KCHAI = KCHJI + NMATIJ(ISYCHO)*NUMV

            IF (LOCDBG) THEN
               KEND2 = KCHAI + NT1AM(ISYCHO)*NUMV
               LWRK2 = LWORK - KEND2 + 1
               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,'(//,5X,A,A,A)')
     &            'Error in batching in ',SECNAM,':'
                  WRITE(LUPRI,*)
     &            'ISYCHO, NUMCHO: ',ISYCHO,NUMCHO(ISYCHO)
                  WRITE(LUPRI,*)
     &            'NVEC, NUMV    : ',NVEC,NUMV
                  WRITE(LUPRI,*)
     &            'IBATCH, NBATCH: ',IBATCH,NBATCH
                  WRITE(LUPRI,*)
     &            'LWRK1, MINMEM : ',LWRK1,MINMEM
                  WRITE(LUPRI,*)
     &            'LWORK, KEND1-1: ',LWORK,KEND1-1
                  CALL QUIT('Error in '//SECNAM)
               ENDIF
            ENDIF

C           Loop over vectors in this batch.        
C           --------------------------------

            DO LVEC = 1,NUMV

C              Read current vector.
C              --------------------

               JVEC = JVEC1 + LVEC - 1
               CALL CHO_READN(WORK(KCHOL),JVEC,1,WORK(KIND1),IDUM2,
     &                        ISYCHO,IOPTR,WORK(KSCR1),LSCR1)

C              Transform.
C              ----------

               DO ISYMI = 1,NSYM

                  NI = NRHF(ISYMI)
                  IF (NI .LE. 0) GOTO 998

                  ISYMA = MULD2H(ISYMI,ISYCHO)
                  ISYMJ = ISYMA
                  ISYMG = ISYMA
                  ISYMD = ISYMI

                  NG    = NBAS(ISYMG)
                  NTOTG = MAX(NG,1)
                  ND    = NBAS(ISYMD)
                  NTOTD = MAX(ND,1)
                  NA    = NVIR(ISYMA)
                  NTOTA = MAX(NA,1)
                  NJ    = NRHF(ISYMJ)
                  NTOTJ = MAX(NJ,1)

                  KOFF1 = KCHOL + IAODIS(ISYMG,ISYMD)
                  KOFF2 = IGLMRH(ISYMD,ISYMI) + 1

                  CALL DGEMM('N','N',NG,NI,ND,
     &                       ONE,WORK(KOFF1),NTOTG,XLAMDH(KOFF2),NTOTD,
     &                       ZERO,WORK(KSCR1),NTOTG)

                  KOFF1 = IGLMRH(ISYMG,ISYMJ) + 1
                  KOFF2 = KCHJI + NMATIJ(ISYCHO)*(LVEC - 1)
     &                  + IMATIJ(ISYMJ,ISYMI)

                  CALL DGEMM('T','N',NJ,NI,NG,
     &                       ONE,XLAMDP(KOFF1),NTOTG,WORK(KSCR1),NTOTG,
     &                       ZERO,WORK(KOFF2),NTOTJ)

                  KOFF1 = IGLMVI(ISYMG,ISYMA) + 1
                  KOFF2 = KCHAI + NT1AM(ISYCHO)*(LVEC - 1)
     &                  + IT1AM(ISYMA,ISYMI)

                  CALL DGEMM('T','N',NA,NI,NG,
     &                       ONE,XLAMDP(KOFF1),NTOTG,WORK(KSCR1),NTOTG,
     &                       ZERO,WORK(KOFF2),NTOTA)

  998             CONTINUE

               ENDDO

            ENDDO

C           Write this batch to disk.
C           -------------------------

            CALL CHO_MOWRITE(WORK(KCHJI),NMATIJ(ISYCHO),NUMV,JVEC1,LUJI)
            CALL CHO_MOWRITE(WORK(KCHAI),NT1AM(ISYCHO),NUMV,JVEC1,LUAI)

C           Read L(ia,J) into the space of L(ai,J).
C           ---------------------------------------

            CALL CHO_MOREAD(WORK(KCHAI),NT1AM(ISYCHO),NUMV,JVEC1,LUIA)

            IF (ISYCHO .EQ. 1) THEN

C              Calculate TRACE vector.
C              -----------------------

               NT1T = MAX(NT1AMX,1)

               CALL DGEMV('T',NT1AMX,NUMV,
     &                    ONE,WORK(KCHAI),NT1T,T1AM,1,
     &                    ZERO,TRACE(JVEC1),1)

C              Calculate Coulomb contribution to F(ia).
C              ----------------------------------------

               CALL DGEMV('N',NT1AMX,NUMV,
     &                    TWO,WORK(KCHAI),NT1T,TRACE(JVEC1),1,
     &                    ONE,FOCK,1)

            ENDIF

C           Calculate exchange contribution to F(ia) and energy.
C           First L(ji,J) vector is overwritten here!!!
C           ----------------------------------------------------

            KUMAT = KCHJI

            DO LVEC = 1,NUMV

C              Calculate U(ji) = sum(a) T1AM(aj) * L(ia,J).
C              --------------------------------------------

               DO ISYMI = 1,NSYM

                  ISYMJ = MULD2H(ISYMI,ISYCHO)
                  ISYMA = ISYMJ

                  KOFF1 = IT1AM(ISYMA,ISYMJ) + 1
                  KOFF2 = KCHAI + NT1AM(ISYCHO)*(LVEC - 1)
     &                  + IT1AM(ISYMA,ISYMI)
                  KOFF3 = KUMAT + IMATIJ(ISYMJ,ISYMI)

                  NA    = NVIR(ISYMA)
                  NTOTA = MAX(NA,1)
                  NJ    = NRHF(ISYMJ)
                  NTOTJ = MAX(NJ,1)
                  NI    = NRHF(ISYMI)

                  CALL DGEMM('T','N',NJ,NI,NA,
     &                       ONE,T1AM(KOFF1),NTOTA,WORK(KOFF2),NTOTA,
     &                       ZERO,WORK(KOFF3),NTOTJ)

               ENDDO

C              Calculate exchange contribution to F(ia).
C              -----------------------------------------

               DO ISYMK = 1,NSYM

                  ISYMI = MULD2H(ISYMK,ISYCHO)
                  ISYMA = ISYMI

                  KOFF1 = KCHAI + NT1AM(ISYCHO)*(LVEC - 1)
     &                  + IT1AM(ISYMA,ISYMK)
                  KOFF2 = KUMAT + IMATIJ(ISYMK,ISYMI)
                  KOFF3 = IT1AM(ISYMA,ISYMI) + 1

                  NA    = NVIR(ISYMA)
                  NTOTA = MAX(NA,1)
                  NI    = NRHF(ISYMI)
                  NK    = NRHF(ISYMK)
                  NTOTK = MAX(NK,1)

                  CALL DGEMM('N','N',NA,NI,NK,
     &                       XMONE,WORK(KOFF1),NTOTA,WORK(KOFF2),NTOTK,
     &                       ONE,FOCK(KOFF3),NTOTA)

               ENDDO

C              Calculate "exchange" energy contribution.
C              -----------------------------------------

               DO ISYMI = 1,NSYM

                  ISYMJ = MULD2H(ISYMI,ISYCHO)

                  IF (NRHF(ISYMJ) .GT. 0) THEN

                     DO I = 1,NRHF(ISYMI)

                        KOFF1 = KUMAT + IMATIJ(ISYMJ,ISYMI)
     &                        + NRHF(ISYMJ)*(I - 1)
                        KOFF2 = KUMAT + IMATIJ(ISYMI,ISYMJ)
     &                        + I - 1

                        ENRX = ENRX
     &                       + DDOT(NRHF(ISYMJ),WORK(KOFF1),1,
     &                                          WORK(KOFF2),NRHF(ISYMI))

                     ENDDO

                  ENDIF

               ENDDO

            ENDDO

         ENDDO

C        Close MO Cholesky files for this symmetry.
C        ------------------------------------------

         CALL CHO_MOP(IKEEP,ITYIA,ISYCHO,LUIA,1,1)
         CALL CHO_MOP(IKEEP,ITYAI,ISYCHO,LUAI,1,1)
         CALL CHO_MOP(IKEEP,ITYJI,ISYCHO,LUJI,1,1)

  999    CONTINUE

      ENDDO

C     Calculate Coulomb energy contribution and sum.
C     ----------------------------------------------

      ENRC = TWO*DDOT(NUMCHO(1),TRACE,1,TRACE,1)
      ECC2 = ECC2 + ENRC - ENRX

      IF (LOCDBG) THEN
         FNORM = DSQRT(DDOT(NT1AMX,FOCK,1,FOCK,1))
         TNORM = DSQRT(DDOT(NUMCHO(1),TRACE,1,TRACE,1))
         WRITE(LUPRI,'(5X,A,A,F16.9,/,5X,A,A,F16.9,/,5X,A,A,F16.9)')
     &   SECNAM,': debug: Coulomb contribution to energy: ',ENRC,
     &   SECNAM,': debug: Exch.   contribution to energy: ',-ENRX,
     &   SECNAM,': debug: Accumulated CC2 energy        : ',ECC2
         WRITE(LUPRI,'(5X,A,A,1P,D15.6)')
     &   SECNAM,': debug: Norm of F(ia): ',FNORM
         WRITE(LUPRI,'(5X,A,A,1P,D15.6)')
     &   SECNAM,': debug: Norm of TRACE: ',TNORM
      ENDIF

      RETURN
      END
C  /* Deck cc2_choamd */
      SUBROUTINE CC2_CHOAMD(DIAG,FOCKD,WORK,LWORK,ISYM,LUCHMO)
C
C     Thomas Bondo Pedersen, August 2002.
C
C     Purpose: Calculate (minus) the diagonal of
C              the CC2 amplitudes of symmetry ISYM. The calculation
C              uses the (ai|bj) integral representation
C              from file LUCHMO; this file is assumed
C              to be open and to contain the transformed
C              AO Cholesky vectors L(ai,J).
C
#include <implicit.h>
      DIMENSION DIAG(*), FOCKD(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>

      INTEGER AI

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CC2_CHOAMD')

      PARAMETER (INFO = 20)
      PARAMETER (TWO = 2.0D0)

C     Start timing.
C     -------------

      TIMTOT = SECOND()

C     Calculate the (ai|ai) integral diagonal.
C     ----------------------------------------

      TIMDIA = SECOND()
      CALL CHO_AIBJD(DIAG,WORK,LWORK,ISYM,LUCHMO)
      TIMDIA = SECOND() - TIMDIA

C     Divide by the orbital energy denominators.
C     t(ai,ai) = (ai|ai)/[2*e(a)-2*e(i)].
C     ------------------------------------------

      TIMDIV = SECOND()
      DO ISYMI = 1,NSYM

         ISYMA = MULD2H(ISYMI,ISYM)

         DO I = 1,NRHF(ISYMI)

            KOFFI = IRHF(ISYMI) + I

            DO A = 1,NVIR(ISYMA)

               KOFFA = IVIR(ISYMA) + A

               AI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A

               DENOM    = TWO*(FOCKD(KOFFA)-FOCKD(KOFFI))
               DIAG(AI) = DIAG(AI)/DENOM

            ENDDO

         ENDDO

      ENDDO
      TIMDIV = SECOND() - TIMDIV

C     Print.
C     ------

      IF (IPRINT .GE. INFO) THEN
         TIMTOT = SECOND() - TIMTOT
         XNRM = DSQRT(DDOT(NT1AM(ISYM),DIAG,1,DIAG,1))
         CALL HEADER('Information from '//SECNAM,-1)
         WRITE(LUPRI,'(5X,A,I1)')
     &   'CC2 amplitude diagonal calculated, symmetry block:',ISYM
         WRITE(LUPRI,'(5X,A,1P,D30.15,/)')
     &   'Norm of calculated diagonal: ',XNRM
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time for calculating (ai|ai) diagonal: ',TIMDIA,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time for orbital energy denominators : ',TIMDIV,' seconds'
         WRITE(LUPRI,'(5X,A)')
     &   '---------------------------------------------------------'
         WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &   'Total time                           : ',TIMTOT,' seconds'
      ENDIF

      RETURN
      END
C  /* Deck cc2_choamd_i */
      SUBROUTINE CC2_CHOAMD_I(DIAG,FOCKD,CHAI,ISYM)
C
C     Thomas Bondo Pedersen, August 2002.
C
C     Purpose: Calculate (minus) the diagonal of
C              the CC2 amplitudes of symmetry ISYM.
C              CHAI must contain the L(ai,J) vectors.
C
#include <implicit.h>
      DIMENSION DIAG(*), FOCKD(*), CHAI(*)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>

      INTEGER AI

      CHARACTER*12 SECNAM
      PARAMETER (SECNAM = 'CC2_CHOAMD_I')

      PARAMETER (INFO = 20)
      PARAMETER (TWO = 2.0D0)

C     Start timing.
C     -------------

      TIMTOT = SECOND()

C     Calculate the (ai|ai) integral diagonal.
C     ----------------------------------------

      TIMDIA = SECOND()
      CALL CHO_AIBJD_I(DIAG,CHAI,ISYM)
      TIMDIA = SECOND() - TIMDIA

C     Divide by the orbital energy denominators.
C     t(ai,ai) = (ai|ai)/[2*e(a)-2*e(i)].
C     ------------------------------------------

      TIMDIV = SECOND()
      DO ISYMI = 1,NSYM

         ISYMA = MULD2H(ISYMI,ISYM)

         DO I = 1,NRHF(ISYMI)

            KOFFI = IRHF(ISYMI) + I

            DO A = 1,NVIR(ISYMA)

               KOFFA = IVIR(ISYMA) + A

               AI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A

               DENOM    = TWO*(FOCKD(KOFFA)-FOCKD(KOFFI))
               DIAG(AI) = DIAG(AI)/DENOM

            ENDDO

         ENDDO

      ENDDO
      TIMDIV = SECOND() - TIMDIV

C     Print.
C     ------

      IF (IPRINT .GE. INFO) THEN
         TIMTOT = SECOND() - TIMTOT
         XNRM = DSQRT(DDOT(NT1AM(ISYM),DIAG,1,DIAG,1))
         CALL HEADER('Information from '//SECNAM,-1)
         WRITE(LUPRI,'(5X,A,I1)')
     &   'CC2 amplitude diagonal calculated, symmetry block:',ISYM
         WRITE(LUPRI,'(5X,A,1P,D30.15,/)')
     &   'Norm of calculated diagonal: ',XNRM
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time for calculating (ai|ai) diagonal: ',TIMDIA,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time for orbital energy denominators : ',TIMDIV,' seconds'
         WRITE(LUPRI,'(5X,A)')
     &   '---------------------------------------------------------'
         WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &   'Total time                           : ',TIMTOT,' seconds'
      ENDIF

      RETURN
      END
C  /* Deck cc2_choam */
      SUBROUTINE CC2_CHOAM(T2AM,FOCKD,WORK,LWORK,LISTBJ,NUMBJ,ISYMBJ,
     &                     LUCHMO)
C
C     Thomas Bondo Pedersen, August 2002.
C
C     Purpose:
C        Calculate minus CC2 amplitudes for a list, LISTBJ, of NUMBJ column
C        (bj) indices of symmetry ISYMBJ.
C
C     Notes:
C        (a) MO Cholesky vectors, L(ai,J), are assumed available on disk on
C            opened unit LUCHMO (in subroutine CHO_AIBJ).
C
#include <implicit.h>
      DIMENSION T2AM(*), FOCKD(*), WORK(LWORK)
      INTEGER LISTBJ(*)
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>

      INTEGER AI,BJ,AIBJ

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC2_CHOAM')

      PARAMETER (INFO = 20)
      PARAMETER (ZERO = 0.0D0)

C     Start timing.
C     -------------

      TIMTOT = SECOND()
      TIMINV = ZERO

C     Calculate (ai|#[bj]) integrals.
C     -------------------------------

      TIMINT = SECOND()
      CALL CHO_AIBJ(T2AM,WORK,LWORK,LISTBJ,NUMBJ,ISYMBJ,LUCHMO)
      TIMINT = SECOND() - TIMINT

C     Divide by the orbital energy denominators.
C     t(ai,bj) = (ai|bj)/[e(a)-e(i)+e(b)-e(j)].
C     ------------------------------------------

      TIMDIV = SECOND()
      ISYMAI = ISYMBJ
      DO IBJ = 1,NUMBJ

         DTIME  = SECOND()
         BJ = LISTBJ(IBJ)
         CALL CC2_INVBJ(BJ,ISYMBJ,B,ISYMB,J,ISYMJ)
         DTIME  = SECOND() - DTIME
         TIMINV = TIMINV   + DTIME

         KOFFJ = IRHF(ISYMJ) + J
         KOFFB = IVIR(ISYMB) + B
         DENBJ = FOCKD(KOFFB) - FOCKD(KOFFJ)

         DO ISYMI = 1,NSYM

            ISYMA = MULD2H(ISYMI,ISYMAI)

            DO I = 1,NRHF(ISYMI)

               KOFFI = IRHF(ISYMI) + I

               DO A = 1,NVIR(ISYMA)

                  KOFFA = IVIR(ISYMA) + A

                  AI   = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A
                  AIBJ = NT1AM(ISYMAI)*(IBJ - 1) + AI

                  DENOM      = DENBJ + FOCKD(KOFFA) - FOCKD(KOFFI)
                  T2AM(AIBJ) = T2AM(AIBJ)/DENOM

               ENDDO

            ENDDO

         ENDDO

      ENDDO
      TIMDIV = SECOND() - TIMDIV

C     Print.
C     ------

      IF (IPRINT .GE. INFO) THEN
         TIMTOT = SECOND() - TIMTOT
         LEN  = NT1AM(ISYMBJ)*NUMBJ
         XNRM = DSQRT(DDOT(LEN,T2AM,1,T2AM,1))
         CALL HEADER('Information from '//SECNAM,-1)
         WRITE(LUPRI,'(5X,A,I1)')
     &   'CC2 amplitudes calculated, symmetry block:',ISYMBJ
         WRITE(LUPRI,'(5X,A,I8,A)')
     &   'CC2 amplitudes calculated for the following ',NUMBJ,
     &   ' columns:'
         WRITE(LUPRI,'(10I8)') (LISTBJ(I),I=1,NUMBJ)
         WRITE(LUPRI,'(5X,A,1P,D30.15,/)')
     &   'Norm of calculated amplitudes: ',XNRM
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time for calculating (ai|bj) integrals: ',TIMINT,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time for orbital energy denominators  : ',TIMDIV,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   '  - of which CC2_INVBJ required       : ',TIMINV,' seconds'
         WRITE(LUPRI,'(5X,A)')
     &   '----------------------------------------------------------'
         WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &   'Total time                            : ',TIMTOT,' seconds'
      ENDIF

      RETURN
      END
C  /* Deck cc2_choam_i */
      SUBROUTINE CC2_CHOAM_I(T2AM,FOCKD,CHAI,CHBJ,LISTBJ,NUMBJ,ISYMBJ)
C
C     Thomas Bondo Pedersen, September 2002.
C
C     Purpose:
C        Calculate minus CC2 amplitudes for a list, LISTBJ, of NUMBJ column
C        (bj) indices of symmetry ISYMBJ.
C
C     Notes:
C        (a) MO Cholesky vectors, L(ai,J), are assumed stored in CHAI.
C
#include <implicit.h>
      DIMENSION T2AM(*), FOCKD(*), CHAI(*), CHBJ(*)
      INTEGER LISTBJ(*)
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>

      INTEGER AI,BJ,AIBJ

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CC2_CHOAM_I')

      PARAMETER (INFO = 20)
      PARAMETER (ZERO = 0.0D0)

C     Start timing.
C     -------------

      TIMTOT = SECOND()
      TIMINV = ZERO

C     Calculate (ai|#[bj]) integrals.
C     -------------------------------

      TIMINT = SECOND()
      CALL CHO_AIBJ_I(T2AM,CHAI,CHBJ,LISTBJ,NUMBJ,ISYMBJ)
      TIMINT = SECOND() - TIMINT

C     Divide by the orbital energy denominators.
C     t(ai,bj) = (ai|bj)/[e(a)-e(i)+e(b)-e(j)].
C     ------------------------------------------

      TIMDIV = SECOND()
      ISYMAI = ISYMBJ
      DO IBJ = 1,NUMBJ

         DTIME  = SECOND()
         BJ = LISTBJ(IBJ)
         CALL CC2_INVBJ(BJ,ISYMBJ,B,ISYMB,J,ISYMJ)
         DTIME  = SECOND() - DTIME
         TIMINV = TIMINV   + DTIME

         KOFFJ = IRHF(ISYMJ) + J
         KOFFB = IVIR(ISYMB) + B
         DENBJ = FOCKD(KOFFB) - FOCKD(KOFFJ)

         DO ISYMI = 1,NSYM

            ISYMA = MULD2H(ISYMI,ISYMAI)

            DO I = 1,NRHF(ISYMI)

               KOFFI = IRHF(ISYMI) + I

               DO A = 1,NVIR(ISYMA)

                  KOFFA = IVIR(ISYMA) + A

                  AI   = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A
                  AIBJ = NT1AM(ISYMAI)*(IBJ - 1) + AI

                  DENOM      = DENBJ + FOCKD(KOFFA) - FOCKD(KOFFI)
                  T2AM(AIBJ) = T2AM(AIBJ)/DENOM

               ENDDO

            ENDDO

         ENDDO

      ENDDO
      TIMDIV = SECOND() - TIMDIV

C     Print.
C     ------

      IF (IPRINT .GE. INFO) THEN
         TIMTOT = SECOND() - TIMTOT
         LEN  = NT1AM(ISYMBJ)*NUMBJ
         XNRM = DSQRT(DDOT(LEN,T2AM,1,T2AM,1))
         CALL HEADER('Information from '//SECNAM,-1)
         WRITE(LUPRI,'(5X,A,I1)')
     &   'CC2 amplitudes calculated, symmetry block:',ISYMBJ
         WRITE(LUPRI,'(5X,A,I8,A)')
     &   'CC2 amplitudes calculated for the following ',NUMBJ,
     &   ' columns:'
         WRITE(LUPRI,'(10I8)') (LISTBJ(I),I=1,NUMBJ)
         WRITE(LUPRI,'(5X,A,1P,D30.15,/)')
     &   'Norm of calculated amplitudes: ',XNRM
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time for calculating (ai|bj) integrals: ',TIMINT,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time for orbital energy denominators  : ',TIMDIV,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   '  - of which CC2_INVBJ required       : ',TIMINV,' seconds'
         WRITE(LUPRI,'(5X,A)')
     &   '----------------------------------------------------------'
         WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &   'Total time                            : ',TIMTOT,' seconds'
      ENDIF

      RETURN
      END
C  /* Deck cc2_invbj */
      SUBROUTINE CC2_INVBJ(BJ,ISYMBJ,B,ISYMB,J,ISYMJ)
C
C     Thomas Bondo Pedersen, August 2002.
C
C     Purpose: "Invert" the compund index BJ of symmetry ISYMBJ
C              to get B,ISYMB and J,ISYMJ.
C
#include <implicit.h>
      INTEGER BJ
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC2_INVBJ')

C     Check for errors.
C     -----------------

      IF ((BJ.LE.0) .OR. (BJ.GT.NT1AM(ISYMBJ))) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'BJ out of bounds in ',SECNAM,':'
         WRITE(LUPRI,'(5X,A,I10,A,I1)')
     &   'BJ = ',BJ,', ISYMBJ = ',ISYMBJ
         IF ((ISYMBJ.LE.0) .OR. (ISYMBJ.GT.NSYM)) THEN
            WRITE(LUPRI,'(5X,A)')
     &      '(The symmetry is out of bounds)'
         ENDIF
         WRITE(LUPRI,*)
         CALL QUIT('Error in '//SECNAM)
      ENDIF

C     Invert.
C     Detailed search only in relevant symmetry block.
C     ------------------------------------------------

      DO ISYMO = 1,NSYM

         ISYMV = MULD2H(ISYMO,ISYMBJ)

         IBJF = IT1AM(ISYMV,ISYMO) + 1
         IBJL = IBJF + NVIR(ISYMV)*NRHF(ISYMO) - 1

         IF ((BJ.GE.IBJF) .AND. (BJ.LE.IBJL)) THEN

            ICOUNT = IBJF - 1

            DO IJ = 1,NRHF(ISYMO)
               DO IB = 1,NVIR(ISYMV)

                  ICOUNT = ICOUNT + 1

                  IF (ICOUNT .EQ. BJ) THEN
                     J     = IJ
                     ISYMJ = ISYMO
                     B     = IB
                     ISYMB = ISYMV
                     GOTO 100
                  ENDIF

               ENDDO
            ENDDO

         ENDIF

      ENDDO

C     Debug: print.
C     -------------

  100 IF (LOCDBG) THEN
         WRITE(LUPRI,'(/,5X,A,A)')
     &   SECNAM,': Debug flag set. Testing bj -> b,j:'
         WRITE(LUPRI,'(5X,A,A,I10,1X,I1)')
     &   SECNAM,': Input : BJ,ISYMBJ: ',BJ,ISYMBJ
         WRITE(LUPRI,'(5X,A,A,I10,1X,I1)')
     &   SECNAM,': Output: B,ISYMB  : ',B,ISYMB
         WRITE(LUPRI,'(5X,A,A,I10,1X,I1)')
     &   SECNAM,': Output: J,ISYMJ  : ',J,ISYMJ
         IERRS = 0
         IERRC = 0
         ISYMVO = MULD2H(ISYMB,ISYMJ)
         IF (ISYMVO .NE. ISYMBJ) IERRS = 1
         IBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
         IF (IBJ .NE. BJ) IERRC = 1
         IERR = IERRS + IERRC
         IF (IERR .NE. 0) THEN
            WRITE(LUPRI,'(5X,A,A,I1,A)')
     &      SECNAM,': ',IERR,' error(s) detected.'
            IF (IERRS .NE. 0) THEN
               WRITE(LUPRI,'(5X,A,A)')
     &         SECNAM,': Symmetry error detected.'
            ENDIF
            IF (IERRC .NE. 0) THEN
               WRITE(LUPRI,'(5X,A,A)')
     &         SECNAM,': Index error detected.'
            ENDIF
            WRITE(LUPRI,*)
            CALL QUIT('Error detected in '//SECNAM)
         ELSE
            WRITE(LUPRI,'(5X,A,A,/)')
     &      SECNAM,': No errors detected.'
         ENDIF
      ENDIF

      RETURN
      END
C  /* Deck chocc2_yi */
      SUBROUTINE CHOCC2_YI(FOCKD,OMEGA1,WORK,LWORK,T2NRM,ECC2)
C
C     Thomas Bondo Pedersen, August 2002.
C
C     Purpose: Calculate Y-intermediates and the I-term:
C
C                 Y(ai,J) = sum(bj) s(ai,bj) * L(jb,J)
C
C                 OMEGA1(ai) = OMEGA1(ai)
C                            + sum(bj) s(ai,bj) * F(jb)
C
C              where s(ai,bj) = 2*t(ai,bj) - t(aj,bi), and
C              t(ai,bj) = -(ai|bj)/[e(a)-e(i)+e(b)-e(j)].
C
C              Calculate contribution to energy from doubles
C              amplitudes:
C
C                 ECC2 = ECC2 + sum(aiJ) Y(ai,J) * L(ia,J)
C
C     Notes:
C        (1) The doubles amplitudes are never stored in full, neither on
C            disk nor in core. There are two main routes, controlled by
C            flag CHOT2:
C            (a) calculation from (ai|bj) integrals (CHOT2=.FALSE.):
C                - represented by the transformed AO Cholesky vectors, OR
C                - represented by vectors from separately Cholesky
C                  decomposed (ai|bj) integrals.
C            (b) decomposition of the amplitudes themselves (CHOT2=.TRUE.).
C        (2) The norm of the packed amplitudes is returned in T2NRM,
C            partly for debugging and partly for mimicking the output
C            of the "standard" code.
C
#include <implicit.h>
      DIMENSION FOCKD(*), OMEGA1(*), WORK(LWORK)
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <chocc2.h>
#include <priunit.h>
#include <iratdef.h>
#include <dummy.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CHOCC2_YI')

      PARAMETER (ZERO = 0.0D0)

C     Decomposition as requested from input.
C     --------------------------------------

      IF (CHOT2 .OR. CHOMO2) THEN

C        Set threshold and span factor.
C        ------------------------------

         IF (THRCC2 .LT. ZERO) THRCC2 = THRCOM
         IF (SPACC2 .LT. ZERO) SPACC2 = SPAN

C        Set requested decomposition.
C        Amplitude deco. precedes integral deco.
C        ---------------------------------------

         IF (CHOT2) THEN
            IMAT = 3
         ELSE IF (CHOMO2) THEN
            IMAT = 2
         ENDIF

C        Scratch allocations.
C        --------------------

         KDIANL = 1
         KCOLUM = KDIANL + 3*NSYM
         KINDIA = KCOLUM + (NSYM - 1)/IRAT   + 1
         KEND1  = KINDIA + (MXDEC2 - 1)/IRAT + 1
         LWRK1  = LWORK  - KEND1

         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Insufficient memory in ',SECNAM,' - allocation: decomp.'
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need (more than): ',KEND1,
     &      'Available       : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Decompose.
C        ----------

         CALL CC_DECMO(WORK(KDIANL),NCHMO2,WORK(KCOLUM),WORK(KINDIA),
     &                 THRCC2,SPACC2,
     &                 FCC2S,.FALSE.,MXDEC2,NCHRD2,NT1AM,IPRINT,IMAT,
     &                 NSYM,THZCC2,FOCKD,WORK(KEND1),LWRK1)

      ENDIF

C     Calculate Y-intermediates, I-term and energy contribution.
C     ----------------------------------------------------------

      IF (IALCC2 .EQ. 1) THEN

         CALL CC2_CHOYI1(FOCKD,OMEGA1,WORK,LWORK,T2NRM,ECC2)

      ELSE IF (IALCC2 .EQ. 2) THEN

         KFCKIA = 1
         KDUM1  = KFCKIA + NT1AM(1)
         KDUM2  = KDUM1  + 1
         KEND1  = KDUM2  + 1
         LWRK1  = LWORK  - KEND1 + 1

         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Insufficient memory in ',SECNAM,
     &      ' - allocation: F(ia)'
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need (more than): ',KEND1-1,
     &      'Available       : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

         CALL ONEL_OP(-1,3,LUFIA)
         CALL CHO_MOREAD(WORK(KFCKIA),NT1AM(1),1,1,LUFIA)
         CALL ONEL_OP(1,3,LUFIA)

         CALL CC_CYIV0(FOCKD,OMEGA1,WORK(KFCKIA),WORK(KEND1),LWRK1,
     &                 T2NRM,T2CNM,.FALSE.,.FALSE.)
         CALL CC_CYIENR(WORK,LWORK,ECC2)

      ELSE

         WRITE(LUPRI,'(//,5X,A,A)')
     &   'Error in ',SECNAM
         WRITE(LUPRI,'(5X,A,I10,/)')
     &   'Algorithm option not recognized: IALCC2 = ',IALCC2
         CALL QUIT('Error in '//SECNAM)

      ENDIF

      RETURN
      END
C  /* Deck cc_cyiv0 */
      SUBROUTINE CC_CYIV0(FOCKD,OMEGA1,FIA,WORK,LWORK,T2NRM,T2CNM,
     &                    DELP1,DELP2)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Calculate Y intermediates and the I term for the CC2
C              vector function.
C
#include <implicit.h>
      DIMENSION FOCKD(*), OMEGA1(*), FIA(*), WORK(LWORK)
      LOGICAL DELP1, DELP2
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>
#include <chocc2.h>

      INTEGER NCHP12(8)

C     Set Cholesky vector type for amplitude assembly.
C     ------------------------------------------------

      FAC1   = 1.0D0
      FAC2   = 0.0D0
      NUMP12 = 1
      NUMUV  = 0
      ISYMP1 = 1
      ISYMP2 = ISYMP1
      SCD    = 1.0D0
      SCDG   = 1.0D0
      IOPTDN = 1
      FREQ   = 0.0D0
      IOPTCE = 1
      IF (CHOT2) THEN
         JTYP1  = 6
         IOPTDN = 0
         DO ISYCHO = 1,NSYM
            NCHP12(ISYCHO) = NCHMO2(ISYCHO)
         ENDDO
      ELSE IF (CHOMO2) THEN
         JTYP1 = 5
         DO ISYCHO = 1,NSYM
            NCHP12(ISYCHO) = NCHMO2(ISYCHO)
         ENDDO
      ELSE
         JTYP1 = 3
         DO ISYCHO = 1,NSYM
            NCHP12(ISYCHO) = NUMCHO(ISYCHO)
         ENDDO
      ENDIF
      JTYP2  = JTYP1

C     Set info for Y intermediates.
C     -----------------------------

      NUMZ  = 1
      JTYPZ = 1
      ISYMZ = 1
      JTYPY = 0

C     Calculate Y intermediates and I-term.
C     -------------------------------------

      CALL CC_CYI(JTYP1,ISYMP1,JTYP2,ISYMP2,NCHP12,FAC1,NUMP12,
     &            DUMMY,IDUMMY,DUMMY,IDUMMY,FAC2,NUMUV,
     &            SCD,SCDG,IOPTDN,FOCKD,FREQ,IOPTCE,
     &            JTYPZ,ISYMZ,NUMCHO,NUMZ,JTYPY,
     &            FIA,1,1,OMEGA1,
     &            WORK,LWORK,T2NRM,T2CNM,DELP1,DELP2)

      RETURN
      END
C  /* Deck cc2_choyi1 */
      SUBROUTINE CC2_CHOYI1(FOCKD,OMEGA1,WORK,LWORK,T2NRM,ECC2)
C
C     Thomas Bondo Pedersen, August 2002.
C
C     Purpose: Calculate Y-intermediates and the I-term.
C              Calculate energy contribution from doubles
C              amplitudes.
C
#include <implicit.h>
      DIMENSION FOCKD(*), OMEGA1(*), WORK(LWORK)
#include <maxorb.h>
#include <ccdeco.h>
#include <ccisvi.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <chocc2.h>

      INTEGER ABEG, AEND
      INTEGER LVIR(8), IOFA1(8), IX1AM(8,8), NX1AM(8), IX2SQ(8)

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CC2_CHOYI1')

      PARAMETER (INFO = 3)
      PARAMETER (IOPEN = -1, IKEEP = 1)
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0)
      PARAMETER (WTOMB = 7.62939453125D-6, D100 = 100.00D0)

C     Initilize timings.
C     ------------------

      TIMT = SECOND()
      TIMG = ZERO
      TIM2 = ZERO
      TIMI = ZERO
      TIMR = ZERO
      TIMY = ZERO

C     Initial check of memory.
C     We must be able to hold at least 1 virtual in each symmetry,
C     as well as at least 1 Cholesky vector in each symmetry.
C     Note: the minimal requirements for amplitude generation,
C     for the Y-intermediates, and for the I-term are identical.
C     ------------------------------------------------------------

      NEED = 0
      DO ISYMA = 1,NSYM
         NEEDG = 0
         DO ISYCHO = 1,NSYM
            ISYMI = MULD2H(ISYMA,ISYCHO)
            NTAMG = NT1AM(ISYCHO) + NRHF(ISYMI)
            NEEDG = MAX(NEEDG,NTAMG)
         ENDDO
         NEEDA = NCKI(ISYMA) + NEEDG
         NEED  = MAX(NEED,NEEDG)
      ENDDO
      IF (LWORK .LT. NEED) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   'Insufficient memory in ',SECNAM
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Minimum memory required: ',NEED,
     &   'Available memory       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Initializations.
C     ----------------

      MXMEMT = 0
      MXMEML = 0
      XLWORK = ONE*LWORK

C     Determine memory split.
C     -----------------------

      CALL MP2_MEMSPL2(LWORK,NUMCHO,LWRK)

C     Print.
C     ------

      IF (IPRINT .GE. INFO) THEN
         CALL AROUND('Output from '//SECNAM)
         XMB  = XLWORK*WTOMB
         LEFT = LWORK - LWRK
         WRITE(LUPRI,'(12X,A,/,12X,A,I10,A,F7.1,A)')
     &   'CC2 amplitudes will be calculated on the fly.',
     &   'Available memory: ',LWORK,' words (',XMB,' Mb).'
         WRITE(LUPRI,'(12X,A,1X,I10,A,/,12X,A,1X,I10,A)')
     &   'Maximum memory allowed for CC2 amplitudes: ',LWRK,' words.',
     &   'Maximum memory allowed for Cholesky vecs.: ',LEFT,' words.'
         WRITE(LUPRI,'(A)') ' '
         WRITE(LUPRI,'(5X,A,A,/,5X,A,A)')
     &   '   #a        First           Last         Mem. Usage',
     &   '        Time   ',
     &   '          (a, sym. a)     (a, sym. a)        (%)    ',
     &   '      (seconds)'
         WRITE(LUPRI,'(5X,A,A)')
     &   '----------------------------------------------------',
     &   '---------------'
      ENDIF

C     Start batched loop over a.
C     --------------------------

      ABEG  = 1
      NUMA  = 0
      NPASS = 0
  100 CONTINUE

         ABEG = ABEG + NUMA
         IF (ABEG .GT. NVIRT) GOTO 200

C        Time this batch.
C        ----------------

         TABAT = SECOND()

C        Find out how many a's can be treated.
C        -------------------------------------

         CALL GET_BTCH(LVIR,ABEG,NUMA,MEM,LWRK)

C        Check for errors.
C        -----------------

         AEND = ABEG + NUMA - 1
         IF (AEND .GT. NVIRT) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Batch error in ',SECNAM,' !'
            WRITE(LUPRI,'(5X,A)')
     &      '- dumping presumably most relevant info:'
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,I10)')
     &      'First A           : ',ABEG,
     &      'Last  A           : ',AEND,
     &      'Number of virtuals: ',NVIRT
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Available memory  : ',LWRK,
     &      'Needed    memory  : ',MEM
            CALL QUIT('Batch error in '//SECNAM)
         ENDIF

         IF (NUMA .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory for a-batch in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/)')
     &      'Available memory: ',LWRK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Get symmetries of first and last a in batch.
C        --------------------------------------------

         ISABEG = ISVI(ABEG)
         ISAEND = ISVI(AEND)

C        Complete allocation for amplitudes.
C        -----------------------------------

         KT2AM = 1
         KEND1 = KT2AM + MEM
         LWRK1 = LWORK - KEND1 + 1

         KLAST  = KEND1 - 1
         MXMEML = MAX(MXMEML,KLAST)
         IF (MXMEML .GT. LWORK) THEN
            ID = 1
            GOTO 2000
         ENDIF

C        Set up local index arrays.
C        --------------------------

         CALL IZERO(IOFA1,NSYM)

         IA1 = ABEG
         DO ISYMA = ISABEG,ISAEND
            IOFA1(ISYMA) = IA1 + NRHFT - IVIR(ISYMA)
            IA1 = IA1 + LVIR(ISYMA)
         ENDDO

         ICOUN2 =  0
         DO ISYMAI = 1,NSYM

            ICOUN1 = 0
            DO ISYMI = 1,NSYM
               ISYMA = MULD2H(ISYMI,ISYMAI)
               IX1AM(ISYMA,ISYMI) = ICOUN1
               ICOUN1 = ICOUN1 + LVIR(ISYMA)*NRHF(ISYMI)
            ENDDO
            NX1AM(ISYMAI) = ICOUN1

            IX2SQ(ISYMAI) = ICOUN2
            ICOUN2 = ICOUN2 + NT1AM(ISYMAI)*NX1AM(ISYMAI)

         ENDDO
         NX2SQ = ICOUN2

         IF (NX2SQ .NE. MEM) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Error detected in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,/)')
     &      'NX2SQ is calculated to be ',NX2SQ,
     &      'MEM determines alloc. as  ',MEM,
     &      'These numbers must be equal....'
            CALL QUIT('Error in '//SECNAM)
         ENDIF

C        Calculate the CC2 amplitudes.
C        Scale by -1 to get the right sign.
C        ----------------------------------

         DTIME = SECOND()
         CALL CC2_CHOAMG(WORK(KT2AM),FOCKD,WORK(KEND1),LWRK1,
     &                   LVIR,IOFA1,IX1AM,NX1AM,IX2SQ,NX2SQ,
     &                   T2NRM,.TRUE.,MEMLEN)
         CALL DSCAL(NX2SQ,XMONE,WORK(KT2AM),1)
         DTIME = SECOND() - DTIME
         TIMG  = TIMG     + DTIME

         KTEST  = KLAST + MEMLEN
         MXMEML = MAX(MXMEML,KTEST)
         IF (MXMEML .GT. LWORK) THEN
            ID = 2
            GOTO 2000
         ENDIF

C        Set up 2 Coulomb minus exchange (almost) in place.
C        --------------------------------------------------

         DTIME = SECOND()
         CALL CC2_TCME(WORK(KT2AM),WORK(KEND1),LWRK1,LVIR,IOFA1,IX1AM,
     &                 NX1AM,IX2SQ,NX2SQ,MEMLEN)
         DTIME = SECOND() - DTIME
         TIM2  = TIM2     + DTIME

         KTEST  = KLAST + MEMLEN
         MXMEML = MAX(MXMEML,KTEST)
         IF (MXMEML .GT. LWORK) THEN
            ID = 3
            GOTO 2000
         ENDIF

C        Calculate the I-term.
C        ---------------------

         DTIME = SECOND()
         CALL CC2_ITERM(OMEGA1,WORK(KT2AM),WORK(KEND1),LWRK1,LVIR,IOFA1,
     &                  IX1AM,NX1AM,IX2SQ,MEMLEN)
         DTIME = SECOND() - DTIME
         TIMI  = TIMI     + DTIME

         KTEST  = KLAST + MEMLEN
         MXMEML = MAX(MXMEML,KTEST)
         IF (MXMEML .GT. LWORK) THEN
            ID = 4
            GOTO 2000
         ENDIF

C        Calculate the Y-intermediates and energy contribution
C        from doubles amplitudes.
C        -----------------------------------------------------

         DTIME = SECOND()
         CALL CC2_Y(WORK(KT2AM),WORK(KEND1),LWRK1,LVIR,IOFA1,IX1AM,
     &              NX1AM,IX2SQ,MEMLEN,ECC2)
         DTIME = SECOND() - DTIME
         TIMY  = TIMY     + DTIME

         KTEST  = KLAST + MEMLEN
         MXMEML = MAX(MXMEML,KTEST)
         IF (MXMEML .GT. LWORK) THEN
            ID = 5
            GOTO 2000
         ENDIF

C        Print information from this a-batch.
C        ------------------------------------

         IF (IPRINT .GE. INFO) THEN
            TABAT = SECOND() - TABAT
            XMUSE = (D100*MXMEML)/XLWORK
            LAFST = IOFA1(ISABEG)
            LALST = IOFA1(ISAEND) + LVIR(ISAEND) - 1
         WRITE(LUPRI,'(5X,I6,4X,I6,1X,I1,8X,I6,1X,I1,9X,F6.2,8X,F10.2)')
     &      NUMA,LAFST,ISABEG,LALST,ISAEND,XMUSE,TABAT
         ENDIF

C        Update memory info.
C        -------------------

         MXMEMT = MAX(MXMEMT,MXMEML)
         MXMEML = 0

C        Update NPASS and go to next a-batch.
C        (Which will exit immediately if no more a's.)
C        ---------------------------------------------

         NPASS = NPASS + 1
         GO TO 100

C     Exit a-batch: calculation done.
C     -------------------------------

  200 CONTINUE

C     Print.
C     ------

      IF (IPRINT .GE. INFO) THEN
         WRITE(LUPRI,'(5X,A,A,/)')
     &   '----------------------------------------------------',
     &   '---------------'
         TIMT = SECOND() - TIMT
         CALL HEADER('Timing of '//SECNAM,-1)
         WRITE(LUPRI,'(5X,A,I10,/)')
     &   'Passes through Cholesky file(s)     : ',NPASS
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for generating amplitudes : ',TIMG,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for two Coulomb - exchange: ',TIM2,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for the I-term            : ',TIMI,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for the Y-intermediates   : ',TIMY,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   '--------------------------------------------------------'
         WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &   'Total time                          : ',TIMT,' seconds'
      ENDIF

C     Normal exit.
C     ------------

      RETURN

C     Memory error.
C     -------------

 2000 WRITE(LUPRI,'(//,5X,A,A,A)')
     & 'Allocation error detected in ',SECNAM,':'
      WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     & 'Last memory location used: ',MXMEML,
     & 'Memory available         : ',LWORK
      WRITE(LUPRI,'(5X,A,I3,/)')
     & 'Identifier: ',ID
      CALL QUIT('Allocation error in '//SECNAM)

      END
C  /* Deck cc2_y */
      SUBROUTINE CC2_Y(T2AM,WORK,LWORK,LVIR,IOFA1,IX1AM,NX1AM,IX2SQ,
     &                 MEMLEN,ECC2)
C
C     Thomas Bondo Pedersen, August 2002.
C
C     Purpose: Calculate the Y-intermediates and write them to disk.
C              Calculate contribution to energy from doubles.
C
C     Note: The Y-intermediate I/O section is probably crappy....
C
#include <implicit.h>
      DIMENSION T2AM(*), WORK(LWORK)
      INTEGER LVIR(8), IOFA1(8), IX1AM(8,8), NX1AM(8), IX2SQ(8)
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <chocc2.h>
#include <ccsdinp.h>
#include <priunit.h>

      INTEGER LROW, ICOL, ADDR

      CHARACTER*5 SECNAM
      PARAMETER (SECNAM = 'CC2_Y')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      PARAMETER (IOPEN = -1, IKEEP = 1, ITYP = 1, INFO = 5)
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)

C     Initialize timings.
C     -------------------

      TIMTOT = SECOND()
      TIMRMO = ZERO
      TIMYCL = ZERO
      TIMENR = ZERO
      TIMWRY = ZERO

C     Initialize MEMLEN.
C     ------------------

      MEMLEN = 0

C     Start Cholesky symmetry loop.
C     -----------------------------

      DO ISYCHO = 1,NSYM

         IF (NUMCHO(ISYCHO) .LE. 0) GOTO 1000
         IF (NT1AM(ISYCHO)  .LE. 0) GOTO 1000
         IF (NX1AM(ISYCHO)  .LE. 0) GOTO 1000

C        Open MO file for this symmetry [L(jb,J)].
C        -----------------------------------------

         DTIME  = SECOND()
         CALL CHO_MOP(IOPEN,ITYP,ISYCHO,LUCHMO,1,1)
         DTIME  = SECOND() - DTIME
         TIMRMO = TIMRMO   + DTIME

C        Open Y-intermediate file for this symmetry.
C        -------------------------------------------

         DTIME  = SECOND()
         CALL WOPEN2(LCC2YM,FCC2YM(ISYCHO),64,0)
         DTIME  = SECOND() - DTIME
         TIMWRY = TIMWRY   + DTIME

C        Set up batch.
C        -------------

         MINMEM = NT1AM(ISYCHO) + NX1AM(ISYCHO)
         NVEC   = MIN(LWORK/MINMEM,NUMCHO(ISYCHO))

         IF (NVEC .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory for batch in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need (pref. multiple) : ',MINMEM,
     &      'Total memory available: ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

         NBATCH = (NUMCHO(ISYCHO) - 1)/NVEC + 1

C        Batch loop.
C        -----------

         DO IBATCH = 1,NBATCH

            NUMV = NVEC
            IF (IBATCH .EQ. NBATCH) THEN
               NUMV = NUMCHO(ISYCHO) - NVEC*(NBATCH - 1)
            ENDIF

            JVEC1 = NVEC*(IBATCH - 1) + 1

C           Complete allocation.
C           --------------------

            KRESY = 1
            KCHOL = KRESY + NX1AM(ISYCHO)*NUMV
            KLAST = KCHOL + NT1AM(ISYCHO)*NUMV - 1

            MEMLEN = MAX(MEMLEN,KLAST)

C           Read vectors.
C           -------------

            DTIME  = SECOND()
            CALL CHO_MOREAD(WORK(KCHOL),NT1AM(ISYCHO),NUMV,JVEC1,LUCHMO)
            DTIME  = SECOND() - DTIME
            TIMRMO = TIMRMO   + DTIME

C           Calculate Y-intermediates.
C           --------------------------

            DTIME = SECOND()

            KOFFT = IX2SQ(ISYCHO) + 1

            NBJ = NT1AM(ISYCHO)
            NAI = NX1AM(ISYCHO)

            CALL DGEMM('T','N',NAI,NUMV,NBJ,
     &                 ONE,T2AM(KOFFT),NBJ,WORK(KCHOL),NBJ,
     &                 ZERO,WORK(KRESY),NAI)

            DTIME  = SECOND() - DTIME
            TIMYCL = TIMYCL   + DTIME

C           Calculate energy contribution from doubles.
C           -------------------------------------------

            DTIME  = SECOND()
            IF (NX1AM(ISYCHO) .EQ. NT1AM(ISYCHO)) THEN
               LTOT = NT1AM(ISYCHO)*NUMV
               ECC2 = ECC2 + DDOT(LTOT,WORK(KRESY),1,WORK(KCHOL),1)
            ELSE
               DO IVEC = 1,NUMV
                  DO ISYMI = 1,NSYM
                     ISYMA = MULD2H(ISYMI,ISYCHO)
                     IF (LVIR(ISYMA) .GT. 0) THEN
                        IF (LVIR(ISYMA) .EQ. NVIR(ISYMA)) THEN
                           KOFF1 = KRESY + NX1AM(ISYCHO)*(IVEC - 1)
     &                           + IX1AM(ISYMA,ISYMI)
                           KOFF2 = KCHOL + NT1AM(ISYCHO)*(IVEC - 1)
     &                           + IT1AM(ISYMA,ISYMI)
                           LENAI = NVIR(ISYMA)*NRHF(ISYMI)
                           ECC2  = ECC2
     &                           + DDOT(LENAI,WORK(KOFF1),1,
     &                                        WORK(KOFF2),1)
                        ELSE
                           DO I = 1,NRHF(ISYMI)
                              KOFF1 = KRESY + NX1AM(ISYCHO)*(IVEC - 1)
     &                              + IX1AM(ISYMA,ISYMI)
     &                              + LVIR(ISYMA)*(I - 1)
                              KOFF2 = KCHOL + NT1AM(ISYCHO)*(IVEC - 1)
     &                              + IT1AM(ISYMA,ISYMI)
     &                              + NVIR(ISYMA)*(I - 1)
     &                              + IOFA1(ISYMA) - 1
                              ECC2  = ECC2
     &                              + DDOT(LVIR(ISYMA),WORK(KOFF1),1,
     &                                                 WORK(KOFF2),1)
                           ENDDO
                        ENDIF
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF
            DTIME  = SECOND() - DTIME
            TIMENR = TIMENR   + DTIME

C           Save Y-intermediates on disk.
C           -----------------------------

            DTIME = SECOND()
            IF (NX1AM(ISYCHO) .EQ. NT1AM(ISYCHO)) THEN

C              Full vector; simply write.
C              --------------------------

               ICOL = JVEC1 - 1
               LROW = NT1AM(ISYCHO)
               ADDR = LROW*ICOL + 1
               LEN  = NT1AM(ISYCHO)*NUMV

               IF (LOCDBG) THEN
                  WRITE(LUPRI,'(5X,A,A)')
     &            SECNAM,': debug: calling PUTWA2:'
                  WRITE(LUPRI,'(7X,A,I1)')
     &            '- full amplitude array of sym. ',ISYCHO
               ENDIF

               CALL PUTWA2(LCC2YM,FCC2YM(ISYCHO),WORK(KRESY),ADDR,LEN)

            ELSE IF (NX1AM(ISYCHO) .LT. NT1AM(ISYCHO)) THEN

C              Part of vector; write to disk in small chunks.
C              ----------------------------------------------

               DO IVEC = 1,NUMV

                  JVEC = JVEC1 + IVEC - 1
                  ICOL = JVEC  - 1
                  LROW = NT1AM(ISYCHO)

                  DO ISYMI = 1,NSYM

                     ISYMA = MULD2H(ISYMI,ISYCHO)

                     IF (LVIR(ISYMA) .GT. 0) THEN

                        IF (LVIR(ISYMA) .EQ. NVIR(ISYMA)) THEN

                           KOFFY = KRESY
     &                           + NX1AM(ISYCHO)*(IVEC - 1)
     &                           + IX1AM(ISYMA,ISYMI)
                           ADDR  = LROW*ICOL
     &                           + IT1AM(ISYMA,ISYMI) + 1

                           IF (LOCDBG) THEN
                              WRITE(LUPRI,'(5X,A,A)')
     &                        SECNAM,': debug: calling PUTWA2:'
                              WRITE(LUPRI,'(7X,A,I1,1X,I1)')
     &                        '- full symmetry block a,i: ',ISYMA,ISYMI
                           ENDIF

c      write(lupri,*) 'LROW, ICOL: ',LROW,ICOL
c      write(lupri,*) 'LROW*ICOL : ',LROW*ICOL
c      write(lupri,*) 'IT1AM     : ',IT1AM(ISYMA,ISYMI)
c      write(lupri,*) 'ADDR expr.: ',LROW*ICOL+IT1AM(ISYMA,ISYMI)
c      write(lupri,*) 'ADDR pass.: ',ADDR

                           CALL PUTWA2(LCC2YM,FCC2YM(ISYCHO),
     &                                  WORK(KOFFY),ADDR,
     &                                  NVIR(ISYMA)*NRHF(ISYMI))

                        ELSE

                           DO I = 1,NRHF(ISYMI)

                              KOFFY = KRESY
     &                              + NX1AM(ISYCHO)*(IVEC - 1)
     &                              + IX1AM(ISYMA,ISYMI)
     &                              + LVIR(ISYMA)*(I - 1)
                              ADDR  = LROW*ICOL
     &                              + IT1AM(ISYMA,ISYMI)
     &                              + NVIR(ISYMA)*(I - 1)
     &                              + IOFA1(ISYMA)

                              IF (LOCDBG) THEN
                                 WRITE(LUPRI,'(5X,A,A)')
     &                           SECNAM,': debug: calling PUTWA2:'
                                 WRITE(LUPRI,'(7X,A,I1,1X,I1)')
     &                           '- symmetry block a,i: ',ISYMA,ISYMI
                                 WRITE(LUPRI,'(7X,A,I10,A,I10)')
     &                           '- writing ',LVIR(ISYMA),' for i = ',I
                              ENDIF

                              CALL PUTWA2(LCC2YM,FCC2YM(ISYCHO),
     &                                     WORK(KOFFY),ADDR,
     &                                     LVIR(ISYMA))

                           ENDDO

                        ENDIF

                     ENDIF

                  ENDDO

               ENDDO

            ELSE

C              Error.
C              ------

               WRITE(LUPRI,'(//,5X,A,A,A)')
     &         'Error detected in ',SECNAM,':'
               WRITE(LUPRI,'(5X,A,I10,A,I1,A)')
     &         'Row dimension of Y-intermediate is ',NX1AM(ISYCHO),
     &         ' (symmetry ',ISYCHO,')'
               WRITE(LUPRI,'(5X,A,I10)')
     &         'Maximum possible dimension is      ',NT1AM(ISYCHO)
               WRITE(LUPRI,'(5X,A,/)')
     &         'The error is probably caused by calling routine.'
               CALL QUIT('Error detected in '//SECNAM)

            ENDIF
            DTIME  = SECOND() - DTIME
            TIMWRY = TIMWRY   + DTIME

         ENDDO

C        Close Y-intermediate file for this symmetry.
C        --------------------------------------------

         DTIME  = SECOND()
         CALL WCLOSE2(LCC2YM,FCC2YM(ISYCHO),'KEEP')
         DTIME  = SECOND() - DTIME
         TIMWRY = TIMWRY   + DTIME

C        Close MO file for this symmetry [L(jb,J)].
C        ------------------------------------------

         DTIME  = SECOND()
         CALL CHO_MOP(IKEEP,ITYP,ISYCHO,LUCHMO,1,1)
         DTIME  = SECOND() - DTIME
         TIMRMO = TIMRMO   + DTIME

 1000    CONTINUE

      ENDDO

C     Print.
C     ------

      IF (IPRINT .GE. INFO) THEN
         TIMTOT = SECOND() - TIMTOT
         CALL HEADER('Timing of '//SECNAM,-1)
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for I/O of MO Cholesky vectors : ',TIMRMO,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for I/O of the  Y-intermediates: ',TIMWRY,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for calculating Y-intermediates: ',TIMYCL,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for calculating energy contr.  : ',TIMENR,' seconds'
         WRITE(LUPRI,'(5X,A)')
     &   '-------------------------------------------------------------'
         WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &   'Total time                               : ',TIMTOT,' seconds'
      ENDIF

      RETURN
      END
C  /* Deck cc2_iterm */
      SUBROUTINE CC2_ITERM(OMEGA1,T2AM,WORK,LWORK,LVIR,IOFA1,IX1AM,
     &                     NX1AM,IX2SQ,MEMLEN)
C
C     Thomas Bondo Pedersen, August 2002.
C
C     Purpose: Calculate the I-term from a batch of doubles amplitudes.
C
#include <implicit.h>
      DIMENSION OMEGA1(*), T2AM(*), WORK(LWORK)
      INTEGER LVIR(8), IOFA1(8), IX1AM(8,8), NX1AM(8), IX2SQ(8)
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC2_ITERM')

      PARAMETER (IOPEN = -1, IKEEP = 1, ITYP = 3)
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)

C     Initialize MEMLEN.
C     ------------------

      MEMLEN = 0

C     Test if anything to do.
C     -----------------------

      IF ((NX1AM(1).LE.0) .OR. (NT1AMX.LE.0)) RETURN

C     Allocation.
C     -----------

      KFOCK = 1
      KSCR1 = KFOCK + NT1AMX
      KEND1 = KSCR1 + NX1AM(1)
      LWRK1 = LWORK - KEND1 + 1

      KLAST  = KEND1 - 1
      MEMLEN = MAX(MEMLEN,KLAST)

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   'Insufficient memory in ',SECNAM
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need     : ',KEND1-1,
     &   'Available: ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Read Fock matrix block F(jb).
C     -----------------------------

      CALL ONEL_OP(IOPEN,ITYP,LUFOCK)
      CALL CHO_MOREAD(WORK(KFOCK),NT1AMX,1,1,LUFOCK)
      CALL ONEL_OP(IKEEP,ITYP,LUFOCK)

C     Calculate contribution to I-term in scratch space.
C     --------------------------------------------------

      KOFFT = IX2SQ(1) + 1
      CALL DGEMV('T',NT1AMX,NX1AM(1),ONE,T2AM(KOFFT),NT1AMX,
     &           WORK(KFOCK),1,ZERO,WORK(KSCR1),1)

C     Add into OMEGA1.
C     ----------------

      DO ISYMI = 1,NSYM

         ISYMA = ISYMI

         IF (LVIR(ISYMA) .LE. 0) GOTO 100

         DO I = 1,NRHF(ISYMI)

            KOFF1 = KSCR1 + IX1AM(ISYMA,ISYMI) + LVIR(ISYMA)*(I - 1)
            KOFF2 = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1)
     &            + IOFA1(ISYMA)

            CALL DAXPY(LVIR(ISYMA),ONE,WORK(KOFF1),1,OMEGA1(KOFF2),1)

         ENDDO

  100    CONTINUE

      ENDDO

      RETURN
      END
C  /* Deck cc2_tcme */
      SUBROUTINE CC2_TCME(T2AM,WORK,LWORK,LVIR,IOFB1,IX1AM,NX1AM,IX2SQ,
     &                    NX2SQ,MEMLEN)
C
C     Thomas Bondo Pedersen, August 2002.
C
C     Purpose: Set up two Coulomb minus exchange in place
C              for a batch of CC2 amplitudes, t(ai,#bj).
C              LWORK must be appr. 2*V.
C
C     The overall structure is taken from the standard "exchange"
C     CCSD_T2TP routine by A. Sanchez and H. Koch.
C
#include <implicit.h>
      DIMENSION T2AM(*), WORK(LWORK)
      INTEGER LVIR(8), IOFB1(8), IX1AM(8,8), NX1AM(8), IX2SQ(8)
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      CHARACTER*8 SECNAM
      PARAMETER (SECNAM = 'CC2_TCME')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      PARAMETER (XMONE = -1.0D0, TWO = 2.0D0)
      PARAMETER (TOL = 1.0D-15)

C     Initialize MEMLEN.
C     ------------------

      MEMLEN = 0

C     Debug: set up tcme in a straightforward fashion
C     using (potentially a lot of) work space.
C     -----------------------------------------------

      IF (LOCDBG) THEN
         MAXVIR = 0
         DO ISYM = 1,NSYM
            MAXVIR = MAX(MAXVIR,NVIR(ISYM))
         ENDDO
         KSCRT = 2*MAXVIR + 1
         KEND0 = KSCRT + NX2SQ
         LWRK0 = LWORK - KEND0 + 1
         KLAST = KEND0 - 1
         MEMLEN = MAX(MEMLEN,KLAST)
         IF (LWRK0 .LT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory for debug in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need     : ',KEND0-1,
     &      'Available: ',LWORK
            CALL QUIT('Insufficient memory for debug in '//SECNAM)
         ENDIF
         DO ISYMBJ = 1,NSYM
            ISYMAI = ISYMBJ
            DO ISYMJ = 1,NSYM
               ISYMB = MULD2H(ISYMJ,ISYMBJ)
               DO J = 1,NRHF(ISYMJ)
                  DO LB = 1,LVIR(ISYMB)
                     LBJ = IX1AM(ISYMB,ISYMJ) + LVIR(ISYMB)*(J - 1)
     &                   + LB
                     DO ISYMI = 1,NSYM
                        ISYMA  = MULD2H(ISYMI,ISYMAI)
                        ISYMAJ = MULD2H(ISYMA,ISYMJ)
                        ISYMBI = MULD2H(ISYMB,ISYMI)
                        DO I = 1,NRHF(ISYMI)
                           LBI = IX1AM(ISYMB,ISYMI)
     &                         + LVIR(ISYMB)*(I - 1) + LB
                           DO A = 1,NVIR(ISYMA)
                              NAI = IT1AM(ISYMA,ISYMI)
     &                            + NVIR(ISYMA)*(I - 1) + A
                              NAJ = IT1AM(ISYMA,ISYMJ)
     &                            + NVIR(ISYMA)*(J - 1) + A
                              NAIBJ = IX2SQ(ISYMBJ)
     &                              + NT1AM(ISYMAI)*(LBJ - 1) + NAI
                              NAJBI = IX2SQ(ISYMBI)
     &                              + NT1AM(ISYMAJ)*(LBI - 1) + NAJ
                              WORK(KSCRT+NAIBJ-1) = TWO*T2AM(NAIBJ)
     &                                            -     T2AM(NAJBI)
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDIF

C     Set up TCME.
C     ------------

      DO ISYMJ = 1,NSYM
         DO J = 1,NRHF(ISYMJ)
            DO ISYMB = 1,NSYM

               ISYMBJ = MULD2H(ISYMB,ISYMJ)
               ISYMAI = ISYMBJ

               DO LB = 1,LVIR(ISYMB)

                  LBJ = IX1AM(ISYMB,ISYMJ) + LVIR(ISYMB)*(J - 1) + LB

                  DO ISYMI = 1,ISYMJ

                     ISYMA  = MULD2H(ISYMI,ISYMAI)
                     ISYMAJ = MULD2H(ISYMA,ISYMJ)
                     ISYMBI = MULD2H(ISYMB,ISYMI)

                     KAIBJ = 1
                     KAJBI = KAIBJ + NVIR(ISYMA)
                     KEND1 = KAJBI + NVIR(ISYMA)
                     LWRK1 = LWORK - KEND1 + 1

                     KLAST  = KEND1 - 1
                     MEMLEN = MAX(MEMLEN,KLAST)

                     IF (LWRK1 .LT. 0) THEN
                      WRITE(LUPRI,'(//,5X,A,A,/,5X,A,I10,/,5X,A,I10,/)')
     &                  'Insufficient memory in ',SECNAM,
     &                  'Need     : ',KEND1-1,
     &                  'Available: ',LWORK
                        CALL QUIT('Insufficient memory in '//SECNAM)
                     ENDIF

                     IF (ISYMI .EQ. ISYMJ) THEN
                        NRHFI = J - 1
                     ELSE
                        NRHFI = NRHF(ISYMI)
                     ENDIF

                     DO I = 1,NRHFI

                        LBI = IX1AM(ISYMB,ISYMI) + LVIR(ISYMB)*(I - 1)
     &                      + LB

                        NAIBJ = IX2SQ(ISYMAI) + NT1AM(ISYMAI)*(LBJ - 1)
     &                        + IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1)
     &                        + 1
                        NAJBI = IX2SQ(ISYMAJ) + NT1AM(ISYMAJ)*(LBI - 1)
     &                        + IT1AM(ISYMA,ISYMJ) + NVIR(ISYMA)*(J - 1)
     &                        + 1

                        CALL DCOPY(NVIR(ISYMA),T2AM(NAIBJ),1,
     &                                         WORK(KAIBJ),1)
                        CALL DCOPY(NVIR(ISYMA),T2AM(NAJBI),1,
     &                                         WORK(KAJBI),1)
                        CALL DSCAL(NVIR(ISYMA),TWO,T2AM(NAIBJ),1)
                        CALL DSCAL(NVIR(ISYMA),TWO,T2AM(NAJBI),1)
                        CALL DAXPY(NVIR(ISYMA),XMONE,WORK(KAJBI),1,
     &                                               T2AM(NAIBJ),1)
                        CALL DAXPY(NVIR(ISYMA),XMONE,WORK(KAIBJ),1,
     &                                               T2AM(NAJBI),1)

                     ENDDO

                  ENDDO

               ENDDO

            ENDDO
         ENDDO
      ENDDO

C     Debug: compare results.
C     -----------------------

      IF (LOCDBG) THEN
         CALL AROUND('Diff. TCME Amplitudes in '//SECNAM)
         CALL DAXPY(NX2SQ,XMONE,T2AM,1,WORK(KSCRT),1)
         DNRM = DDOT(NX2SQ,WORK(KSCRT),1,WORK(KSCRT),1)
         XNX2 = 1.0D0*NX2SQ
         RMS  = DSQRT(DNRM/XNX2)
         WRITE(LUPRI,'(10X,A,I10,1X,1P,D14.6,/)')
     &   'Number of amplitudes and RMS error: ',NX2SQ,RMS
         IF (RMS .LE. TOL) GOTO 100
         WRITE(LUPRI,'(10X,A,1P,D14.6,A,/)')
     &   'Printing all differences larger than ',TOL,':'
         WRITE(LUPRI,'(1X,A)')
     &   '    A ISYMA     I ISYMI     B ISYMB     J ISYMJ     Amplitude'
         WRITE(LUPRI,'(1X,A)')
     &   '-------------------------------------------------------------'
         DO ISYMBJ = 1,NSYM
            ISYMAI = ISYMBJ
            DO ISYMJ = 1,NSYM
               ISYMB = MULD2H(ISYMJ,ISYMBJ)
               DO J = 1,NRHF(ISYMJ)
                  DO LB = 1,LVIR(ISYMB)
                     LBJ = IX1AM(ISYMB,ISYMJ) + LVIR(ISYMB)*(J - 1) + LB
                     B   = IOFB1(ISYMB) + LB - 1
                     DO ISYMI = 1,NSYM
                        ISYMA = MULD2H(ISYMI,ISYMAI)
                        DO I = 1,NRHF(ISYMI)
                           DO A = 1,NVIR(ISYMA)
                              NAI = IT1AM(ISYMA,ISYMI)
     &                            + NVIR(ISYMA)*(I - 1) + A
                              NAIBJ = KSCRT + IX2SQ(ISYMBJ)
     &                              + NT1AM(ISYMAI)*(LBJ - 1) + NAI - 1
                              IF (DABS(WORK(NAIBJ)) .GT. TOL) THEN
                                 WRITE(LUPRI,99)
     &                           A,ISYMA,I,ISYMI,B,ISYMB,J,ISYMJ,
     &                           WORK(NAIBJ)
                              ENDIF
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         WRITE(LUPRI,'(1X,A)')
     &   '-------------------------------------------------------------'
  100    CONTINUE
      ENDIF

      RETURN
   99 FORMAT(1X,I5,3X,I1,3X,I5,3X,I1,3X,I5,3X,I1,3X,I5,3X,I1,2X,1P,
     &       D14.6)
      END
C  /* Deck cc2_choamg */
      SUBROUTINE CC2_CHOAMG(T2AM,FOCKD,WORK,LWORK,LVIR,IOFB1,IX1AM,
     &                      NX1AM,IX2SQ,NX2SQ,T2NRM,CLNRM,MEMLN)
C
C     Thomas Bondo Pedersen, August 2002.
C
C     Purpose: Generate a batch of (minus) CC2 amplitudes t(ai,#bj).
C
C     Input:
C
C        FOCKD -- Orbital energies in CC ordering.
C                 Dummy if amplitudes are generated directly
C                 from T2 decomposition. (CHOT2 = .TRUE. in  chocc2.h).
C
C        WORK  -- Work space, dimension: LWORK.
C
C        LVIR  -- Number of b's in each symmetry.
C                 I.e., a "local" version of NVIR(*).
C
C        IOFB1 -- Offset to first b in each symmetry.
C                 I.e., IOFB1(ISYMB) returns the global index
C                 of the first b (in LVIR) of symmetry ISYMB.
C
C        IX1AM -- Offset to symmetry blocks of bj.
C                 I.e., a "local" version of IT1AM(*,*).
C
C        NX1AM -- Number of bj's in each symmetry.
C                 I.e., a "local" version of NT1AM(*).
C
C        IX2SQ -- Offset to symmetry blocks ai,bj in T2AM.
C                 I.e., a "local" version of IT2SQ(*,*)
C                 except that only 1 symmetry argument is
C                 used, as T2AM is tot. symmetric.
C
C        NX2SQ -- The total dimension of T2AM.
C                 I.e., a "local" version of NT2SQ(1).
C
C        CLNRM -- .TRUE.: calculate T2NRM.
C
C     Output:
C
C        T2AM  -- The amplitudes sorted according to IX2SQ.
C
C        T2NRM -- Norm of the PACKED T2AM array.
C                 (Must be initialized on input,
C                  only calculated if CLNRM = .TRUE.)
C
C        MEMLN -- Length of work space used (for statistics).
C
C     Notes:
C
C        There are currently 3 possible generations controlled
C        by input in chocc2.h:
C
C           (1) Generation from L(ai,J) vectors, i.e. transformed
C               AO Cholesky vectors (CHOMO2 = .FALSE. and CHOT2 = .FALSE.).
C
C           (2) Generation from separate decomposition of (ai|bj)
C               integrals (CHOMO2 = .TRUE. and CHOT2 = .FALSE.).
C
C           (3) Generation from separately decomposed CC2 amplitudes
C               (CHOT2 = .TRUE.)
C
#include <implicit.h>
      DIMENSION T2AM(*), FOCKD(*), WORK(LWORK)
      INTEGER LVIR(8), IOFB1(8), IX1AM(8,8), NX1AM(8), IX2SQ(8)
      LOGICAL CLNRM
#include <maxorb.h>
#include <ccdeco.h>
#include <chocc2.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>
#include <ccsdinp.h>

      INTEGER AI,BJ,AIBJ
      INTEGER NVECTOT(8)

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CC2_CHOAMG')

      PARAMETER (IOPEN = -1, IKEEP = 1, INFO = 3)
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)

C     Initialize timings.
C     -------------------

      TIMT = SECOND()
      TIMR = ZERO
      TIMS = ZERO
      TIMC = ZERO
      TIMD = ZERO
      TIMN = ZERO

C     Initialize MEMLN.
C     -----------------

      MEMLN = 0

C     Debug: print input.
C     -------------------

      IF (LOCDBG) THEN
         CALL HEADER('Input Variables Entering '//SECNAM,-1)
         WRITE(LUPRI,'(5X,A)')
     &   'LVIR - local NVIR array:'
         WRITE(LUPRI,'(8I10)') (LVIR(ISYM),ISYM=1,NSYM)
         WRITE(LUPRI,'(5X,A)')
     &   'IOFB1 - offsets to virtuals:'
         WRITE(LUPRI,'(8I10)') (IOFB1(ISYM),ISYM=1,NSYM)
         WRITE(LUPRI,'(5X,A)')
     &   'IX1AM - local IT1AM array:'
         DO JSYM = 1,NSYM
            WRITE(LUPRI,'(8I10)') (IX1AM(ISYM,JSYM),ISYM=1,NSYM)
         ENDDO
         WRITE(LUPRI,'(5X,A)')
     &   'NX1AM - local NT1AM array:'
         WRITE(LUPRI,'(8I10)') (NX1AM(ISYM),ISYM=1,NSYM)
         WRITE(LUPRI,'(5X,A)')
     &   'IX2SQ - local IT2SQ array:'
         WRITE(LUPRI,'(8I10)') (IX2SQ(ISYM),ISYM=1,NSYM)
         WRITE(LUPRI,'(5X,A,I10,/)')
     &   'NX2SQ - local dimension of T2AM: ',NX2SQ
      ENDIF

C     If nothing to do, return.
C     -------------------------

      IF (NX2SQ .LE. 0) THEN
         IF (LOCDBG) THEN
            WRITE(LUPRI,'(5X,A,A,I10,/)')
     &      SECNAM,' returns immediately: NX2SQ = ',NX2SQ
         ENDIF
         RETURN
      ENDIF

C     Initialize T2AM.
C     ----------------

      CALL DZERO(T2AM,NX2SQ)

C     Set number of vectors and the type of vectors in CHO_MOP.
C     ---------------------------------------------------------

      IF (CHOT2) THEN
         ITYP = 6
         DO ISYM = 1,NSYM
            NVECTOT(ISYM) = NCHMO2(ISYM)
         ENDDO
      ELSE IF (CHOMO2) THEN
         ITYP = 5
         DO ISYM = 1,NSYM
            NVECTOT(ISYM) = NCHMO2(ISYM)
         ENDDO
      ELSE
         ITYP = 3
         DO ISYM = 1,NSYM
            NVECTOT(ISYM) = NUMCHO(ISYM)
         ENDDO
      ENDIF

C     Start Cholesky symmetry loop.
C     -----------------------------

      DO ISYCHO = 1,NSYM

         ISYMAI = ISYCHO
         ISYMBJ = ISYMAI

C        If nothing to do, skip this symmetry.
C        -------------------------------------

         IF (NVECTOT(ISYCHO) .LE. 0) GOTO 1000
         IF (NT1AM(ISYMAI)   .LE. 0) GOTO 1000
         IF (NX1AM(ISYMBJ)   .LE. 0) GOTO 1000

C        Open file for this symmetry.
C        ----------------------------

         DTIME = SECOND()
         CALL CHO_MOP(IOPEN,ITYP,ISYCHO,LUCHMO,1,1)
         DTIME = SECOND() - DTIME
         TIMR  = TIMR     + DTIME

C        Set up batch.
C        -------------

         MINMEM = NT1AM(ISYMAI) + NX1AM(ISYMBJ)
         NVEC   = MIN(LWORK/MINMEM,NVECTOT(ISYCHO))

         IF (NVEC .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory for batch in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Minimum memory required: ',MINMEM,
     &      'Total memory available : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

         NBATCH = (NVECTOT(ISYCHO) - 1)/NVEC + 1

C        Batch loop.
C        -----------

         DO IBATCH = 1,NBATCH

            NUMV = NVEC
            IF (IBATCH .EQ. NBATCH) THEN
               NUMV = NVECTOT(ISYCHO) - NVEC*(NBATCH - 1)
            ENDIF

            JVEC1 = NVEC*(IBATCH - 1) + 1

C           Complete allocation.
C           --------------------

            KCHAI = 1
            KCHBJ = KCHAI + NT1AM(ISYMAI)*NUMV
            KLAST = KCHBJ + NX1AM(ISYMBJ)*NUMV - 1

            MEMLN = MAX(MEMLN,KLAST)

C           Read vectors.
C           -------------

            DTIME = SECOND()
            CALL CHO_MOREAD(WORK(KCHAI),NT1AM(ISYMAI),NUMV,JVEC1,LUCHMO)
            DTIME = SECOND() - DTIME
            TIMR  = TIMR     + DTIME

C           Copy out the #bj block.
C           -----------------------

            DTIME = SECOND()
            DO IVEC = 1,NUMV
               DO ISYMJ = 1,NSYM

                  ISYMB = MULD2H(ISYMJ,ISYMBJ)
                  IF (LVIR(ISYMB) .LE. 0) GOTO 999

                  DO J = 1,NRHF(ISYMJ)

                     KOFF1 = KCHAI + NT1AM(ISYCHO)*(IVEC - 1)
     &                     + IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1)
     &                     + IOFB1(ISYMB) - 1
                     KOFF2 = KCHBJ + NX1AM(ISYCHO)*(IVEC - 1)
     &                     + IX1AM(ISYMB,ISYMJ) + LVIR(ISYMB)*(J - 1)

                     CALL DCOPY(LVIR(ISYMB),WORK(KOFF1),1,
     &                                      WORK(KOFF2),1)

                  ENDDO

  999             CONTINUE

               ENDDO
            ENDDO
            DTIME = SECOND() - DTIME
            TIMS  = TIMS     + DTIME

C           Calculate the product.
C           ----------------------

            DTIME = SECOND()

            NTOAI = MAX(NT1AM(ISYMAI),1)
            NTOBJ = MAX(NX1AM(ISYMBJ),1)

            KOFFT = IX2SQ(ISYCHO) + 1

            CALL DGEMM('N','T',NT1AM(ISYMAI),NX1AM(ISYMBJ),NUMV,
     &                 ONE,WORK(KCHAI),NTOAI,WORK(KCHBJ),NTOBJ,
     &                 ONE,T2AM(KOFFT),NTOAI)

            DTIME = SECOND() - DTIME
            TIMC  = TIMC     + DTIME

         ENDDO

C        Close file for this symmetry.
C        -----------------------------

         DTIME = SECOND()
         CALL CHO_MOP(IKEEP,ITYP,ISYCHO,LUCHMO,1,1)
         DTIME = SECOND() - DTIME
         TIMR  = TIMR     + DTIME

 1000    CONTINUE

      ENDDO

C     If necessary, divide by orbital energy denominators.
C     ----------------------------------------------------

      IF (.NOT. CHOT2) THEN

         TIMD = SECOND()

         DO ISYMBJ = 1,NSYM

            IF (NX1AM(ISYMBJ) .LE. 0) GOTO 1002

            ISYMAI = ISYMBJ

            DO ISYMJ = 1,NSYM

               ISYMB = MULD2H(ISYMJ,ISYMBJ)

               IF (LVIR(ISYMB) .LE. 0) GOTO 1001

               DO J = 1,NRHF(ISYMJ)

                  KOFFJ = IRHF(ISYMJ) + J

                  DO LB = 1,LVIR(ISYMB)

                     B     = IOFB1(ISYMB) + LB - 1
                     KOFFB = IVIR(ISYMB)  + B

                     LBJ = IX1AM(ISYMB,ISYMJ) + LVIR(ISYMB)*(J - 1) + LB

                     DENBJ = FOCKD(KOFFB) - FOCKD(KOFFJ)

                     DO ISYMI = 1,NSYM

                        ISYMA = MULD2H(ISYMI,ISYMAI)

                        DO I = 1,NRHF(ISYMI)

                           KOFFI = IRHF(ISYMI) + I

                           DNBJI = DENBJ - FOCKD(KOFFI)

                           DO A = 1,NVIR(ISYMA)

                              KOFFA = IVIR(ISYMA) + A

                              AI = IT1AM(ISYMA,ISYMI)
     &                           + NVIR(ISYMA)*(I - 1) + A

                              AIBJ = IX2SQ(ISYMBJ)
     &                             + NT1AM(ISYMAI)*(LBJ - 1) + AI

                              DENOM = DNBJI + FOCKD(KOFFA)

                              T2AM(AIBJ) = T2AM(AIBJ)/DENOM

                           ENDDO

                        ENDDO

                     ENDDO

                  ENDDO

               ENDDO

 1001          CONTINUE

            ENDDO

 1002       CONTINUE

         ENDDO

         TIMD = SECOND() - TIMD

      ENDIF

C     If requested, calculate T2NRM.
C     ------------------------------

      IF (CLNRM) THEN

         TIMN = SECOND()

         DO ISYMBJ = 1,NSYM
            ISYMAI = ISYMBJ
            DO ISYMJ = 1,NSYM
               ISYMB = MULD2H(ISYMJ,ISYMBJ)
               IF (LVIR(ISYMB) .GT. 0) THEN
                  DO J = 1,NRHF(ISYMJ)
                     DO LB = 1,LVIR(ISYMB)
                        LBJ = IX1AM(ISYMB,ISYMJ) + LVIR(ISYMB)*(J - 1)
     &                      + LB
                        B   = IOFB1(ISYMB) + LB - 1
                        BJ  = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1)
     &                      + B
                        DO AI = 1,BJ
                           AIBJ  = IX2SQ(ISYMBJ)
     &                           + NT1AM(ISYMAI)*(LBJ - 1) + AI
                           T2NRM = T2NRM + T2AM(AIBJ)*T2AM(AIBJ)
                        ENDDO
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
         ENDDO

         TIMN = SECOND() - TIMN

      ENDIF

C     Print.
C     ------

      IF (IPRINT .GE. INFO) THEN
         TIMT = SECOND() - TIMT
         CALL HEADER('Information from '//SECNAM,-1)
         IF (CHOT2) THEN
            WRITE(LUPRI,'(5X,A,/)')
     &      'CC2 doubles amplitudes separately decomposed'
         ELSE IF (CHOMO2) THEN
            WRITE(LUPRI,'(5X,A,/)')
     &      'CC2 doubles amplitudes from separately decomposed (ai|bj)'
         ELSE
            WRITE(LUPRI,'(5X,A,/)')
     &      'CC2 doubles amplitudes calculated from (ai|bj)'
         ENDIF
         IF (CLNRM) THEN
            WRITE(LUPRI,'(5X,A,1P,D30.15)')
     &      'Accumulated amplitude norm: ',DSQRT(T2NRM)
         ENDIF
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for I/O of Cholesky vectors: ',TIMR,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for sort of MO vectors     : ',TIMS,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for multiplying vectors    : ',TIMC,' seconds'
         IF (.NOT. CHOT2) THEN
            WRITE(LUPRI,'(5X,A,F10.2,A)')
     &      'Time used for orbital denominators   : ',TIMD,' seconds'
         ENDIF
         IF (CLNRM) THEN
            WRITE(LUPRI,'(5X,A,F10.2,A)')
     &      'Time used for calculating T2 norm    : ',TIMN,' seconds'
         ENDIF
         WRITE(LUPRI,'(5X,A)')
     &   '---------------------------------------------------------'
         WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &   'Time used in total                   : ',TIMT,' seconds'
      ENDIF

      RETURN
      END
C  /* Deck chocc2_jg */
      SUBROUTINE CHOCC2_JG(T1AM,OMEGA1,XLAMDP,XLAMDH,TRACE,WORK,LWORK)
C
C     Thomas Bondo Pedersen, August 2002.
C
C     Purpose: Calculate the CC2 JG-term.
C
C     Notice: TRACE is scaled by 2 on exit !!!
C
#include <implicit.h>
      DIMENSION T1AM(*), OMEGA1(*), XLAMDP(*), XLAMDH(*), TRACE(*)
      DIMENSION WORK(LWORK)
#include <priunit.h>
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <chocc2.h>

      INTEGER IOFFC(8), IOFFZ(8), IOFFY(8)

      INTEGER LROW, ICOL, ADDR

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CHOCC2_JG')

      PARAMETER (IOPEN = -1, IKEEP = 1, ITYKI = 4, IOPTR = 2)
      PARAMETER (INFO = 5)
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)

C     Initialize timings.
C     -------------------

      TIMTOT = SECOND()
      TIMYIO = ZERO
      TIMCIO = ZERO
      TIMMIO = ZERO
      TIMYAD = ZERO
      TIMYBK = ZERO
      TIMYSO = ZERO
      TIMCSO = ZERO
      TIMOMA = ZERO
      TIMOMM = ZERO

C     Scale trace by 2.
C     -----------------

      CALL DSCAL(NUMCHO(1),TWO,TRACE,1)

C     Read reduce index array.
C     ------------------------

      DTIME = SECOND()

      KIND1 = 1
      CALL CC_GETIND1(WORK(KIND1),LWORK,LIND1)
      KEND0 = KIND1 + LIND1
      LWRK0 = LWORK - KEND0 + 1

      IF (LWRK0 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: index'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need (more than): ',KEND0-1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

      DTIME  = SECOND() - DTIME
      TIMCIO = TIMCIO   + DTIME

C     Allocation.
C     -----------

      KOMAO = KEND0
      KEND1 = KOMAO + NT1AO(1)
      LWRK1 = LWORK - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: Omega'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need (more than): ',KEND1-1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Initialize AO result vector.
C     ----------------------------

      CALL DZERO(WORK(KOMAO),NT1AO(1))

C     Start Cholesky symmetry loop.
C     -----------------------------

      DO ISYCHO = 1,NSYM

         IF (NUMCHO(ISYCHO) .LE. 0) GOTO 1000
         IF (N2BST(ISYCHO)  .LE. 0) GOTO 1000
         IF (NT1AM(ISYCHO)  .LE. 0) GOTO 1000
         IF (NT1AO(ISYCHO)  .LE. 0) GOTO 1000

C        Open file with Y-intermediates.
C        -------------------------------

         DTIME  = SECOND()
         CALL WOPEN2(LCC2YM,FCC2YM(ISYCHO),64,0)
         DTIME  = SECOND() - DTIME
         TIMYIO = TIMYIO   + DTIME

C        Open file with L(ki,J).
C        -----------------------

         DTIME  = SECOND()
         CALL CHO_MOP(IOPEN,ITYKI,ISYCHO,LUCHKI,1,1)
         DTIME  = SECOND() - DTIME
         TIMMIO = TIMMIO   + DTIME

C        Allocation.
C        -----------

         LREAD = MEMRD(1,ISYCHO,IOPTR)

         KSCRC = KEND1
         KREAD = KSCRC + N2BST(ISYCHO)
         KEND2 = KREAD + LREAD
         LWRK2 = LWORK - KEND2 + 1

         IF (LWRK2 .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Insufficient memory in ',SECNAM,' - allocation: Cholesky'
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need (more than): ',KEND2-1,
     &      'Available       : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Set up batch.
C        -------------

         LENC   = MAX(N2BST(ISYCHO),NT1AM(ISYCHO),NMATIJ(ISYCHO))
         LENZ   = MAX(NT1AO(ISYCHO),NT1AM(ISYCHO))
         MINMEM = LENC + LENZ
         NVEC   = MIN(LWRK2/MINMEM,NUMCHO(ISYCHO))

         IF (NVEC .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory for batch in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,I10,/)')
     &      'Minimum memory needed: ',MINMEM,
     &      'Available for batch  : ',LWRK2,
     &      'Available in total   : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

         NBATCH = (NUMCHO(ISYCHO) - 1)/NVEC + 1

C        Start batch loop.
C        -----------------

         DO IBATCH = 1,NBATCH

            NUMV = NVEC
            IF (IBATCH .EQ. NBATCH) THEN
               NUMV = NUMCHO(ISYCHO) - NVEC*(NBATCH - 1)
            ENDIF

C           Set first vector in this batch.
C           -------------------------------

            JVEC1 = NVEC*(IBATCH - 1) + 1

C           Set up index arrays IOFFC, IOFFY, and IOFFZ.
C           --------------------------------------------

            ICOUN1 = 0
            ICOUN2 = 0
            ICOUN3 = 0

            DO ISYMB = 1,NSYM

               ISYMA = MULD2H(ISYMB,ISYCHO)
               ISYMI = ISYMA
               ISYMC = ISYMB

               IOFFC(ISYMB) = ICOUN1
               IOFFZ(ISYMB) = ICOUN2
               IOFFY(ISYMC) = ICOUN3

               LENBJ  = NBAS(ISYMB)*NUMV
               LENCJ  = NVIR(ISYMC)*NUMV
               ICOUN1 = ICOUN1 + NBAS(ISYMA)*LENBJ
               ICOUN2 = ICOUN2 + LENBJ*NRHF(ISYMI)
               ICOUN3 = ICOUN3 + LENCJ*NRHF(ISYMI)

            ENDDO

C           Complete allocation.
C           --------------------

            KCHOL = KEND2
            KZMAT = KCHOL + LENC*NUMV

C           Read Y(ci,#J) in place of Z.
C           ----------------------------

            DTIME = SECOND()

            ICOL = JVEC1 - 1
            LROW = NT1AM(ISYCHO)
            ADDR = LROW*ICOL + 1
            LEN  = NT1AM(ISYCHO)*NUMV
            CALL GETWA2(LCC2YM,FCC2YM(ISYCHO),WORK(KZMAT),ADDR,LEN)

            DTIME  = SECOND() - DTIME
            TIMYIO = TIMYIO   + DTIME

C           Read L(ki,#J) in place of L(al,be,#J).
C           --------------------------------------

            DTIME  = SECOND()
            CALL CHO_MOREAD(WORK(KCHOL),NMATIJ(ISYCHO),NUMV,JVEC1,
     &                      LUCHKI)
            DTIME  = SECOND() - DTIME
            TIMMIO = TIMMIO   + DTIME

C           Subtract sum(k) T1AM(ck) * L(ki,#J) from Y(ci,#J).
C           --------------------------------------------------

            DTIME = SECOND()

            DO IVEC = 1,NUMV
               DO ISYMI = 1,NSYM

                  ISYMK = MULD2H(ISYMI,ISYCHO)
                  ISYMC = ISYMK

                  KOFF1 = IT1AM(ISYMC,ISYMK) + 1
                  KOFF2 = KCHOL + NMATIJ(ISYCHO)*(IVEC - 1)
     &                  + IMATIJ(ISYMK,ISYMI)
                  KOFF3 = KZMAT + NT1AM(ISYCHO)*(IVEC - 1)
     &                  + IT1AM(ISYMC,ISYMI)

                  NC    = NVIR(ISYMC) 
                  NTOTC = MAX(NC,1)
                  NI    = NRHF(ISYMI)
                  NK    = NRHF(ISYMK)
                  NTOTK = MAX(NK,1)

                  CALL DGEMM('N','N',NC,NI,NK,
     &                       XMONE,T1AM(KOFF1),NTOTC,WORK(KOFF2),NTOTK,
     &                       ONE,WORK(KOFF3),NTOTC)

               ENDDO
            ENDDO

            DTIME  = SECOND() - DTIME
            TIMYAD = TIMYAD   + DTIME

C           Reorder to obtain Y(c#J,i) in place of L(al,be,#J).
C           ---------------------------------------------------

            DTIME = SECOND()

            DO IVEC = 1,NUMV
               DO ISYMI = 1,NSYM

                  ISYMC = MULD2H(ISYMI,ISYCHO)

                  LENCJ = NVIR(ISYMC)*NUMV

                  DO I = 1,NRHF(ISYMI)

                     KOFF1 = KZMAT + NT1AM(ISYCHO)*(IVEC - 1)
     &                     + IT1AM(ISYMC,ISYMI) + NVIR(ISYMC)*(I - 1)
                     KOFF2 = KCHOL + IOFFY(ISYMC) + LENCJ*(I - 1)
     &                     + NVIR(ISYMC)*(IVEC - 1)

                     CALL DCOPY(NVIR(ISYMC),WORK(KOFF1),1,WORK(KOFF2),1)

                  ENDDO

               ENDDO
            ENDDO

            DTIME  = SECOND() - DTIME
            TIMYSO = TIMYSO   + DTIME

C           Calculate Z(be#J,i) by backtransforming virtual index of Y.
C           -----------------------------------------------------------

            DTIME = SECOND()
            DO ISYMI = 1,NSYM

               ISYMC = MULD2H(ISYMI,ISYCHO)
               ISYMB = ISYMC

               KOFF1 = IGLMVI(ISYMB,ISYMC) + 1
               KOFF2 = KCHOL + IOFFY(ISYMC)
               KOFF3 = KZMAT + IOFFZ(ISYMB)

               NTOTB = MAX(NBAS(ISYMB),1)
               NTOTC = MAX(NVIR(ISYMC),1)

               CALL DGEMM('N','N',
     &                    NBAS(ISYMB),NUMV*NRHF(ISYMI),NVIR(ISYMC),
     &                    ONE,XLAMDH(KOFF1),NTOTB,WORK(KOFF2),NTOTC,
     &                    ZERO,WORK(KOFF3),NTOTB)

            ENDDO
            DTIME  = SECOND() - DTIME
            TIMYBK = TIMYBK   + DTIME

C           Include Coulomb part from J-term.
C           N.B.: TRACE has already been scaled by 2 above.
C           -----------------------------------------------

            IF (ISYCHO .EQ. 1) THEN

               DTIME = SECOND()
               DO ISYMB = 1,NSYM

                  IF (NBAS(ISYMB) .GT. 0) THEN

                     ISYMI = ISYMB

                     LENBJ = NBAS(ISYMB)*NUMV

                     DO I = 1,NRHF(ISYMI)

                        KOFF2 = IGLMRH(ISYMB,ISYMI)
     &                        + NBAS(ISYMB)*(I - 1) + 1

                        DO IVEC = 1,NUMV

                           KOFF1 = JVEC1 + IVEC - 1
                           KOFF3 = KZMAT + IOFFZ(ISYMB)
     &                           + LENBJ*(I - 1)
     &                           + NBAS(ISYMB)*(IVEC - 1)

                           CALL DAXPY(NBAS(ISYMB),TRACE(KOFF1),
     &                                XLAMDH(KOFF2),1,WORK(KOFF3),1)

                        ENDDO

                     ENDDO

                  ENDIF

               ENDDO
               DTIME  = SECOND() - DTIME
               TIMYAD = TIMYAD   + DTIME

            ENDIF

C           Set up L(al,be#J).
C           ------------------

            DO IVEC = 1,NUMV

               DTIME = SECOND()
               JVEC  = JVEC1 + IVEC - 1
               CALL CHO_READN(WORK(KSCRC),JVEC,1,WORK(KIND1),IDUM2,
     &                        ISYCHO,IOPTR,WORK(KREAD),LREAD)
               DTIME  = SECOND() - DTIME
               TIMCIO = TIMCIO   + DTIME

               DTIME = SECOND()
               DO ISYMB = 1,NSYM

                  ISYMA = MULD2H(ISYMB,ISYCHO)

                  LENAB = NBAS(ISYMA)*NBAS(ISYMB)

                  KOFF1 = KSCRC + IAODIS(ISYMA,ISYMB)
                  KOFF2 = KCHOL + IOFFC(ISYMB) + LENAB*(IVEC - 1)

                  CALL DCOPY(LENAB,WORK(KOFF1),1,WORK(KOFF2),1)

               ENDDO
               DTIME  = SECOND() - DTIME
               TIMCSO = TIMCSO   + DTIME

            ENDDO

C           Calculate contribution in AO basis.
C           -----------------------------------

            DTIME = SECOND()
            DO ISYMB = 1,NSYM

               ISYMI = MULD2H(ISYMB,ISYCHO)
               ISYMA = ISYMI

               KOFF1 = KCHOL + IOFFC(ISYMB)
               KOFF2 = KZMAT + IOFFZ(ISYMB)
               KOFF3 = KOMAO + IT1AO(ISYMA,ISYMI)
 
               LENBJ = NBAS(ISYMB)*NUMV

               NTOTA = MAX(NBAS(ISYMA),1)
               NTOBJ = MAX(LENBJ,1)

               CALL DGEMM('N','N',NBAS(ISYMA),NRHF(ISYMI),LENBJ,
     &                    ONE,WORK(KOFF1),NTOTA,WORK(KOFF2),NTOBJ,
     &                    ONE,WORK(KOFF3),NTOTA)

            ENDDO
            DTIME  = SECOND() - DTIME
            TIMOMA = TIMOMA   + DTIME

         ENDDO

C        Close file with L(ki,J).
C        ------------------------

         DTIME  = SECOND()
         CALL CHO_MOP(IKEEP,ITYKI,ISYCHO,LUCHKI,1,1)
         DTIME  = SECOND() - DTIME
         TIMMIO = TIMMIO   + DTIME

C        Close file with Y-intermediates.
C        --------------------------------

         DTIME  = SECOND()
         CALL WCLOSE2(LCC2YM,FCC2YM(ISYCHO),'KEEP')
         DTIME  = SECOND() - DTIME
         TIMYIO = TIMYIO   + DTIME

 1000    CONTINUE

      ENDDO

C     Transform AO index to virtual.
C     ------------------------------

      DTIME = SECOND()
      DO ISYMI = 1,NSYM

         ISYMA = ISYMI
         ISYMG = ISYMA

         NTOTG = MAX(NBAS(ISYMG),1)
         NTOTA = MAX(NVIR(ISYMA),1)

         KOFF1 = IGLMVI(ISYMG,ISYMA) + 1
         KOFF2 = KOMAO + IT1AO(ISYMG,ISYMI)
         KOFF3 = IT1AM(ISYMA,ISYMI)  + 1

         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMG),
     &              ONE,XLAMDP(KOFF1),NTOTG,WORK(KOFF2),NTOTG,
     &              ONE,OMEGA1(KOFF3),NTOTA)

      ENDDO
      DTIME  = SECOND() - DTIME
      TIMOMM = TIMOMM   + DTIME

C     Print.
C     ------

      IF (IPRINT .GE. INFO) THEN
         TIMTOT = SECOND() - TIMTOT
         CALL HEADER('Timing of '//SECNAM,-1)
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for I/O of Y-intermediates  : ',TIMYIO,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for J-contributions to Z    : ',TIMYAD,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for reorderering Y          : ',TIMYSO,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for backtransformation of Y : ',TIMYBK,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for I/O of MO Cholesky vecs.: ',TIMMIO,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for I/O of AO Cholesky vecs.: ',TIMCIO,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for reorder. AO  Chol. vecs.: ',TIMCSO,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for calculating  AO Omega1  : ',TIMOMA,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for transforming AO Omega1  : ',TIMOMM,' seconds'
         WRITE(LUPRI,'(5X,A)')
     &   '----------------------------------------------------------'
         WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &   'Total time                            : ',TIMTOT,' seconds'
      ENDIF

      RETURN
      END
C  /* Deck chocc2_h */
      SUBROUTINE CHOCC2_H(OMEGA1,WORK,LWORK)
C
C     Thomas Bondo Pedersen, August 2002.
C
C     Purpose: Calculate the CC2 H-term.
C
#include <implicit.h>
      DIMENSION OMEGA1(*), WORK(LWORK)
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <chocc2.h>
#include <priunit.h>

      INTEGER IOFFX(8), IOFFC(8)

      INTEGER LROW, ICOL, ADDR

      CHARACTER*8 SECNAM
      PARAMETER (SECNAM = 'CHOCC2_H')

      PARAMETER (IOPEN = -1, IKEEP = 1, INFO = 5)
      PARAMETER (ITYKI = 4)
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0)

C     Initialize timings.
C     -------------------

      TIMTOT = SECOND()
      TIMCIO = ZERO
      TIMYIO = ZERO
      TIMXIM = ZERO
      TIMCAL = ZERO

C     Start Cholesky symmetry loop.
C     -----------------------------

      DO ISYCHO = 1,NSYM

         IF (NUMCHO(ISYCHO) .LE. 0) GOTO 1000
         IF (NT1AM(ISYCHO)  .LE. 0) GOTO 1000
         IF (NMATIJ(ISYCHO) .LE. 0) GOTO 1000

C        Open files for MO Cholesky vectors.
C        -----------------------------------

         DTIME  = SECOND()
         CALL CHO_MOP(IOPEN,ITYKI,ISYCHO,LUCHKI,1,1)
         DTIME  = SECOND() - DTIME
         TIMCIO = TIMCIO   + DTIME

C        Open file for Y-intermediates.
C        ------------------------------

         DTIME  = SECOND()
         CALL WOPEN2(LCC2YM,FCC2YM(ISYCHO),64,0)
         DTIME  = SECOND() - DTIME
         TIMYIO = TIMYIO   + DTIME

C        Allocation: scratch space for read.
C        -----------------------------------

         KSCR1 = 1
         KEND1 = KSCR1 + MAX(NT1AM(ISYCHO),NMATIJ(ISYCHO))
         LWRK1 = LWORK - KEND1 + 1

         IF (LWRK1 .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Insufficient memory in ',SECNAM,' - allocation: scratch'
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need (more than): ',KEND1-1,
     &      'Available       : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Set up batch.
C        -------------

         MINMEM = NT1AM(ISYCHO) + NMATIJ(ISYCHO)
         NVEC   = MIN(LWRK1/MINMEM,NUMCHO(ISYCHO))

         IF (NVEC .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory for batch in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,I10,/)')
     &      'Minimum memory required for batch: ',MINMEM,
     &      'Available for batch              : ',LWRK1,
     &      'Available in total               : ',LWORK
            CALL QUIT('Insufficient memory for batch in '//SECNAM)
         ENDIF

         NBATCH = (NUMCHO(ISYCHO) - 1)/NVEC + 1

C        Start batch loop.
C        -----------------

         DO IBATCH = 1,NBATCH

            NUMV = NVEC
            IF (IBATCH .EQ. NBATCH) THEN
               NUMV = NUMCHO(ISYCHO) - NVEC*(NBATCH - 1)
            ENDIF

            JVEC1 = NVEC*(IBATCH - 1) + 1

C           Complete allocation.
C           --------------------

            KYMAT = KEND1
            KCHKI = KYMAT + NT1AM(ISYCHO)*NUMV

C           Set up index arrays IOFFX and IOFFC.
C           ------------------------------------

            ICOUN1 = 0
            ICOUN2 = 0

            DO ISYMK = 1,NSYM

               ISYMA = MULD2H(ISYMK,ISYCHO)
               ISYMI = ISYMA

               IOFFX(ISYMK) = ICOUN1
               IOFFC(ISYMK) = ICOUN2

               LENKJ  = NRHF(ISYMK)*NUMV
               ICOUN1 = ICOUN1 + NVIR(ISYMA)*LENKJ
               ICOUN2 = ICOUN2 + LENKJ*NRHF(ISYMI)

            ENDDO

            DO IVEC = 1,NUMV

C              Set current vector.
C              -------------------

               JVEC = JVEC1 + IVEC - 1

C              Read L(ki,J).
C              -------------

               DTIME  = SECOND()
               CALL CHO_MOREAD(WORK(KSCR1),NMATIJ(ISYCHO),1,JVEC,LUCHKI)
               DTIME  = SECOND() - DTIME
               TIMCIO = TIMCIO   + DTIME

C              Reorder: L(k#J,i).
C              ------------------

               DTIME  = SECOND()
               DO ISYMI = 1,NSYM

                  ISYMK = MULD2H(ISYMI,ISYCHO)

                  IF (NRHF(ISYMK) .GT. 0) THEN
                     LENKJ = NRHF(ISYMK)*NUMV
                     DO I = 1,NRHF(ISYMI)

                        KOFF1 = KSCR1 + IMATIJ(ISYMK,ISYMI)
     &                        + NRHF(ISYMK)*(I - 1)
                        KOFF2 = KCHKI + IOFFC(ISYMK)
     &                        + LENKJ*(I - 1) + NRHF(ISYMK)*(IVEC - 1)

                        CALL DCOPY(NRHF(ISYMK),WORK(KOFF1),1,
     &                                         WORK(KOFF2),1)

                     ENDDO
                  ENDIF

               ENDDO
               DTIME  = SECOND() - DTIME
               TIMXIM = TIMXIM   + DTIME

C              Read Y(ak,J).
C              -------------

               DTIME  = SECOND()
               ICOL   = JVEC - 1
               LROW   = NT1AM(ISYCHO)
               ADDR   = LROW*ICOL + 1
               CALL GETWA2(LCC2YM,FCC2YM(ISYCHO),WORK(KSCR1),ADDR,
     &                      NT1AM(ISYCHO))
               DTIME  = SECOND() - DTIME
               TIMYIO = TIMYIO   + DTIME

C              Reorder: Y(a,k#J).
C              ------------------

               DTIME  = SECOND()
               DO ISYMK = 1,NSYM

                  ISYMA = MULD2H(ISYMK,ISYCHO)

                  LENAK = NVIR(ISYMA)*NRHF(ISYMK)

                  KOFF1 = KSCR1 + IT1AM(ISYMA,ISYMK)
                  KOFF2 = KYMAT + IOFFX(ISYMK)
     &                  + LENAK*(IVEC - 1)

                  CALL DCOPY(LENAK,WORK(KOFF1),1,WORK(KOFF2),1)

               ENDDO
               DTIME  = SECOND() - DTIME
               TIMXIM = TIMXIM   + DTIME

            ENDDO

C           Calculate contribution.
C           -----------------------

            DTIME  = SECOND()
            DO ISYMK = 1,NSYM

               ISYMA = MULD2H(ISYMK,ISYCHO)
               ISYMI = ISYMA

               KOFF1 = KYMAT + IOFFX(ISYMK)
               KOFF2 = KCHKI + IOFFC(ISYMK)
               KOFF3 = IT1AM(ISYMA,ISYMI) + 1

               LENKJ = NRHF(ISYMK)*NUMV
               NTOTA = MAX(NVIR(ISYMA),1)
               NTOKJ = MAX(LENKJ,1)

               CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),LENKJ,
     &                    XMONE,WORK(KOFF1),NTOTA,WORK(KOFF2),NTOKJ,
     &                    ONE,OMEGA1(KOFF3),NTOTA)

            ENDDO
            DTIME  = SECOND() - DTIME
            TIMCAL = TIMCAL   + DTIME

         ENDDO

C        Close file for Y-intermediates.
C        -------------------------------

         DTIME  = SECOND()
         CALL WCLOSE2(LCC2YM,FCC2YM(ISYCHO),'KEEP')
         DTIME  = SECOND() - DTIME
         TIMYIO = TIMYIO   + DTIME

C        Close file for MO Cholesky vectors.
C        -----------------------------------

         DTIME  = SECOND()
         CALL CHO_MOP(IKEEP,ITYKI,ISYCHO,LUCHKI,1,1)
         DTIME  = SECOND() - DTIME
         TIMCIO = TIMCIO   + DTIME

 1000    CONTINUE

      ENDDO

C     Print.
C     ------

      IF (IPRINT .GE. INFO) THEN
         TIMTOT = SECOND() - TIMTOT
         CALL HEADER('Timing of '//SECNAM,-1)
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for I/O of Y-intermediates   : ',TIMYIO,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for I/O of MO Chol. vectors  : ',TIMCIO,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for reordering factor arrays : ',TIMXIM,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for calculating contribution : ',TIMCAL,' seconds'
         WRITE(LUPRI,'(5X,A)')
     &   '-----------------------------------------------------------'
         WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &   'Total time                             : ',TIMTOT,' seconds'
      ENDIF

      RETURN
      END
C  /* Deck cc2_iniome */
      SUBROUTINE CC2_INIOME(FOCKD,T1AM,OMEGA1,ISYM)
C
C     Thomas Bondo Pedersen, August 2002.
C
C     Purpose: Initialize the CC2 vector function:
C
C              OMEGA1(ai) = [e(a) - e(i)] * T1AM(ai)
C
C
#include <implicit.h>
      DIMENSION FOCKD(*), T1AM(*), OMEGA1(*)
#include <ccorb.h>
#include <ccsdsym.h>

      INTEGER AI

      DO ISYMI = 1,NSYM

         ISYMA = MULD2H(ISYMI,ISYM)

         DO I = 1,NRHF(ISYMI)

            KOFFI = IRHF(ISYMI) + I

            DO A = 1,NVIR(ISYMA)

               KOFFA = IVIR(ISYMA) + A
               AI    = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A

               OMEGA1(AI) = (FOCKD(KOFFA) - FOCKD(KOFFI))*T1AM(AI)

            ENDDO

         ENDDO

      ENDDO

      RETURN
      END
C
C /* Deck cho_fckh */
      SUBROUTINE CHO_FCKH(T1AM,WORK,LWORK)
C
C     Write Fock matrix in CC_FCKH  needed in CC property calculation
C     Use Cholesky vectors
C
C     asm September 2005
C
#include <implicit.h>
#include <maxorb.h>
#include <priunit.h>
#include <dummy.h>
#include <inftap.h>
#include <ccorb.h>
#include <ccdeco.h>
#include <chocc2.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <ccfield.h>
C
      LOGICAL TINFO, IOPTPV
C
      DIMENSION T1AM(*)
      DIMENSION WORK(LWORK)
C
C
C---------------
C     Initialize
C---------------
C
      TINFO   = IPRINT .GE. 5
      ISYDEN = 1
C
C-----------------
C     Allocation.1
C-----------------
C
      KDENSI = 1
      KLAMDP = KDENSI + N2BST(ISYDEN)
      KLAMDH = KLAMDP + NLAMDT
      KEND0  = KLAMDH + NLAMDT
      LWRK0  = LWORK  - KEND0 + 1
C
      IF (LWRK0 .LT. 0) CALL QUIT('Not enough space in CHO_FCKH.0')
C
C------------------------------
C     Construct density matrix
C------------------------------
C
c     IOPT = 1
c     CALL CC_RDRSP('R0 ',1,1,IOPT,MODEL,WORK(KT1AM),DUMMY)
C
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),T1AM,WORK(KEND0),LWRK0)
C
      IC = 1
      CALL CC_AODENS(WORK(KLAMDP),WORK(KLAMDH),WORK(KDENSI),ISYDEN,
     &               IC,WORK(KEND0),LWRK0)
C
      IF (IPRINT .GT. 50) THEN
         CALL AROUND( 'Density matrix' )
         CALL CC_PRFCKAO(WORK(KDENSI),ISYDEN)
      END IF
C
C-----------------
C     Allocation.2
C-----------------
C
      KFOCK = KLAMDP
      KEND1 = KFOCK + N2BST(ISYDEN)
      LWRK1 = LWORK - KEND1 + 1
C
      IF (LWRK1 .LT. 0) CALL QUIT('Not enough space in CHO_FCKH.1')
C
C--------------------------
C     Construct Fock matrix
C--------------------------
C
      CALL DZERO(WORK(KFOCK),N2BST(ISYDEN))
C
      CALL CCRHS_ONEAO(WORK(KFOCK),WORK(KEND1),LWRK1)
C
      DO I = 1, NFIELD
         FIELD = EFIELD(I)
         CALL CC_ONEP(WORK(KFOCK),WORK(KEND1),LWRK1,FIELD,1,LFIELD(I))
      END DO
C
      IOPTPV = .TRUE.
      CALL SCF_CHOFCK3(WORK(KDENSI),WORK(KFOCK),WORK(KEND1),LWRK1,
     &                 ISYDEN,TINFO,IOPTPV)
C
C------------------
C     Write do disk
C------------------
C
      LUFCK = -1
      CALL GPOPEN(LUFCK,'CC_FCKH','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUFCK)
      WRITE(LUFCK)(WORK(KFOCK+I-1), I = 1,N2BST(ISYDEN))
      CALL GPCLOSE(LUFCK,'KEEP' )
C
      IF (IPRINT .GT. 50) THEN
         CALL AROUND( 'Fock AO matrix written to disk' )
         CALL CC_PRFCKAO(WORK(KFOCK),ISYDEN)
      END IF
C
      RETURN
      END

