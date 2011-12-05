C  /* Deck cc_chopden */
      SUBROUTINE CC_CHOPDEN(IPDEN,NPDEN,LISTL,LISTR,LISTD,
     &                      ORDER,ISIDE,FILDEN,DELDEN,WORK,LWORK)
C
C     Thomas Bondo Pedersen, April 2003.
C
C     Purpose: Driver routine for calculating perturbed one-electron
C              density matrices for Cholesky CC2.
C
C     Input:
C     ======
C
C     IPDEN(1,*) --- Indices of LISTL vectors.
C     IPDEN(2,*) --- Indices of LISTR vectors.
C     NPDEN      --- Number of perturbed density matrices requested.
C     LISTL      --- List specifier for "left" vectors.
C     LISTR      --- List specifier for "right" vectors.
C     LISTD      --- List specifier for output densities.
C     ORDER      --- Integer perturbation order.
C     ISIDE      --- Integer specifier of "side":
C                    -1: Left.
C                    +1: Right.
C     FILDEN     --- File for storage of MO density blocks.
C     DELDEN     --- Flag for deleting FILDEN before exit.
C
C     Output:
C     =======
C
C     IPDEN(3,*) --- Indices of output densities on disk, as specified
C                    through LISTD.
C
C     Note:
C     This routine is actually far less general than the argument list
C     suggests. Since the amplitudes are constructed on-the-fly, each
C     list (L or R) must be coded explicitly with the current structure
C     of the program.
C
#include <implicit.h>
      CHARACTER*(*) LISTL, LISTR, LISTD
      CHARACTER*10  FILDEN
      LOGICAL   DELDEN
      INTEGER   IPDEN(3,NPDEN)
      INTEGER   ORDER
      DIMENSION WORK(LWORK)
#include <maxorb.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccdeco.h>
#include <ccsdinp.h>
#include <priunit.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CC_CHOPDEN')

      LOGICAL LOCDBG, DENDEL
      PARAMETER (LOCDBG = .FALSE.)

      CHARACTER*10 MODEL

C     Return if nothing to do.
C     ------------------------

      IF (NPDEN .LE. 0) RETURN

C     Initialize error counter.
C     -------------------------

      JERR = 0

C     Check that this is Cholesky CC2.
C     --------------------------------

      IF (.NOT. CHOINT) THEN
         WRITE(LUPRI,'(5X,A,A)')
     &   SECNAM,': ERROR: integrals must be Cholesky decomposed!'
         JERR = JERR + 1
      ENDIF
      IF (CC2) THEN
         MODEL = 'CC2       '
      ELSE
         WRITE(LUPRI,'(5X,A,A)')
     &   SECNAM,': ERROR: only possible model is CC2!'
         JERR = JERR + 1
      ENDIF

C     Quit if any errors detected.
C     ----------------------------

      IF (JERR .NE. 0) THEN
         IF (JERR .EQ. 1) THEN
            WRITE(LUPRI,'(5X,A,A,I2,A,/)')
     &      SECNAM,':',JERR,' error detected!'
         ELSE
            WRITE(LUPRI,'(5X,A,A,I2,A,/)')
     &      SECNAM,':',JERR,' errors detected!'
         ENDIF
         CALL QUIT('Error in '//SECNAM)
      ENDIF

C     Call appropriate routine to do the job.
C     ---------------------------------------

      IF (LISTL(1:2) .EQ. 'L0') THEN

         IF ((ORDER.EQ.1) .AND. (ISIDE.EQ.1)) THEN

            IF (LOCDBG) THEN
               DENDEL = .FALSE.
            ELSE
               DENDEL = DELDEN
            ENDIF

            CALL CC_CHOR1D(IPDEN,NPDEN,LISTR,LISTD,MODEL,WORK,LWORK,
     &                     FILDEN,DENDEL,.FALSE.)

            IF (LOCDBG) THEN
                CALL CC_CHOR1DBGAO(IPDEN,NPDEN,LISTR,LISTD,MODEL,
     &                             WORK,LWORK,.FALSE.)
                CALL CC_CHOR1DBGMO(IPDEN,NPDEN,LISTR,LISTD,MODEL,
     &                             WORK,LWORK,FILDEN,DELDEN,.TRUE.)
            ENDIF

         ELSE

            WRITE(LUPRI,'(//,5X,A,A,I2,A,I2,A)')
     &      SECNAM,': ERROR: ORDER =',ORDER,' ISIDE =',ISIDE,
     &      ' not implemented!'
            CALL QUIT('Error in '//SECNAM)

         ENDIF

      ELSE

         WRITE(LUPRI,'(//,5X,A,A,A,A)')
     &   SECNAM,': ERROR: LISTL = ',LISTL(1:3),' not implemented!'
         CALL QUIT('Error in '//SECNAM)

      ENDIF

      RETURN
      END
C  /* Deck cc_chor1d */
      SUBROUTINE CC_CHOR1D(IPDEN,NPDEN,LISTR,LISTD,MODEL,WORK,LWORK,
     &                     FILDEN,DELDEN,SKPSNG)
C
C     Thomas Bondo Pedersen, April 2003.
C
C     Purpose: Calculate the right first-order CC2 density matrix using
C              Cholesky decomposed integrals. The algorithm is
C              quasi-T2-direct: tbar(#ai,dl), for all dl, is stored
C              temporarily on disk. The backtransformed AO density
C              is returned on disk as specified through LISTD; indices
C              on disk is returned in IPDEN(3,*).
C
C     Formulas:
C
C     D(jk) = - sum(b)   R(bj) * D0(bk)
C             - sum(aib) R(ai,bj) * tbar(ai,bk)
C     D(jb) =   sum(i)   D0(ji) * R(bi)
C             - sum(a)   R(aj) * D0(ab)
C             + sum(ai)  [2*R(ai,bj)-R(aj,bi)] * D0(ai)
C     D(cb) =   sum(i)   D0(ci) * R(bi)
C             + sum(aij) R(ai,bj) * tbar(ai,cj)
C
C     where tbar denotes 0th order multipliers and D0 the 0th order density
C     which must be available on disk (D0(ai) = tbar(ai)). The contributions
C     from singles R(ai) and 0th order densities may be skipped by setting
C     the SKPSNG flag.
C
C     The MO blocks are stored consecutively on file FILDEN (via long
C     integer crayio), which is deleted in the end if requested through
C     flag DELDEN. Thus, the ordering of this file is (for each amplitude
C     in LISTR as specified in IPDEN(2,*)):
C
C        D(jk),D(jb),D(cb).
C
#include <implicit.h>
      INTEGER IPDEN(3,NPDEN)
      CHARACTER*(*) LISTR, LISTD
      CHARACTER*10  MODEL, FILDEN
      DIMENSION WORK(LWORK)
      LOGICAL   DELDEN, SKPSNG
#include <maxorb.h>
#include <ccdeco.h>
#include <ccisvi.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <dccorb.h>
#include <dccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>
#include <dummy.h>
#include <ciarc.h>
#include <ccr1rsp.h>
#include <chocc2.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_CHOR1D')

      LOGICAL LOCDBG, DEL
      PARAMETER (LOCDBG = .FALSE.)

      LOGICAL LONLYAO
      DATA LONLYAO /.TRUE./

      INTEGER ADR

      INTEGER BBEG, BEND
      INTEGER ABEG, AEND
      INTEGER LAIJ(8)

      PARAMETER (INFO = 3)
      PARAMETER (KTYP1 = 1, KTYP2 = 2, KTYP3 = 3, KTYP4 = 4)
      PARAMETER (KTYP5 = 5, KTYP6 = 6, KTYP7 = 7, KTYP8 = 8)
      PARAMETER (KTYP9 = 9, KTYP10 = 10, KTYP11 = 11)
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (D100 = 1.0D2)

      PARAMETER (MXTYP = 2)
      DIMENSION SCD(MXTYP)
      INTEGER JTYP1(MXTYP), JTYP2(MXTYP), ISYMP1(MXTYP), ISYMP2(MXTYP)
      INTEGER NCHP12(8,MXTYP)
      LOGICAL DELP1(MXTYP), DELP2(MXTYP)

      CHARACTER*8 LABEL, FILMUL(8)

      DATA FILMUL /'SCRMUL_1','SCRMUL_2','SCRMUL_3','SCRMUL_4', 
     &             'SCRMUL_5','SCRMUL_6','SCRMUL_7','SCRMUL_8'/ 

      INTEGER ILSTSYM

C     Initialize timings.
C     -------------------

      TIMTOT = SECOND()
      TIMMUL = ZERO
      TIMRAM = ZERO
      TIMDEN = ZERO
      TIMMOL = ZERO
      TIMMOR = ZERO

C     Other initializations.
C     ----------------------

      MXMEMT = 0
      MXMEML = 0
      XLWORK = ONE*LWORK

C     Print.
C     ------

      IF (IPRINT .GE. INFO) THEN
         CALL AROUND('Output from '//SECNAM)
         WRITE(LUPRI,'(5X,A,A,1X,A,A,I10)')
     &   'Number of ',LISTD(1:3),MODEL,' densities requested:',NPDEN
         IF (SKPSNG) THEN
            WRITE(LUPRI,'(5X,A)')
     &      'Singles contributions NOT included!'
         ELSE
            WRITE(LUPRI,'(5X,A)')
     &      'Singles contributions included.'
         ENDIF
         WRITE(LUPRI,'(5X,A,I10,A)')
     &   'Available memory: ',LWORK,' words.'
         WRITE(LUPRI,'(5X,A)')
     &   'Left  amplitudes are 0th order multipliers.'
         WRITE(LUPRI,'(5X,A,A,A)')
     &   'Right amplitudes are of type ',LISTR(1:3),':'
         IF (LISTR(1:2) .EQ. 'R1') THEN
            WRITE(LUPRI,'(8X,A,/,8X,A)')
     &      'List #    Symmetry   Label       Frequency',
     &      '------------------------------------------'
            DO IVEC = 1,NPDEN
               ILSTNR = IPDEN(2,IVEC)
               ISYDEN = ILSTSYM(LISTR,ILSTNR)
               LABEL  = LRTLBL(ILSTNR)
               FREQ   = FRQLRT(ILSTNR)
               WRITE(LUPRI,'(8X,I6,4X,I4,6X,A8,2X,F12.5)')
     &         ILSTNR,ISYDEN,LABEL,FREQ
            ENDDO
            WRITE(LUPRI,'(8X,A)')
     &      '------------------------------------------'
            WRITE(LUPRI,'(5X,A,A)')
     &      'File used for storing MO density blocks: ',FILDEN
            IF (DELDEN) THEN
               WRITE(LUPRI,'(5X,A,A)')
     &         '- deleted on exit from ',SECNAM
            ENDIF
            WRITE(LUPRI,*)
         ELSE
            WRITE(LUPRI,'(5X,A)')
     &      'Sorry! This type has not yet been implemented!'
            CALL QUIT('Error in '//SECNAM)
         ENDIF
         CALL FLSHFO(LUPRI)
      ENDIF

C     Calculate singles/D0 contributions if requested, else
C     initialize densities on disk.
C     -----------------------------------------------------

      IF (SKPSNG) THEN
         CALL CC_CHODINI(IPDEN,NPDEN,LISTR,WORK,LWORK,FILDEN)
      ELSE
         TIMSNG = SECOND()
         CALL CC_CHOR1DS(IPDEN,NPDEN,LISTR,WORK,LWORK,MODEL,FILDEN)
         TIMSNG = SECOND() - TIMSNG
      ENDIF

C     Set up partial info for CC_CAA depending on LISTR.
C     --------------------------------------------------

      IF (LISTR(1:2) .EQ. 'R1') THEN
         FAC1      = 1.0D0
         FAC2      = 0.0D0
         SCDG      = 1.0D0
         IOPTDN    = 1
         IOPTCE    = 0
         NUMP12    = 2
         NUMUV     = 0
         JTYP1(1)  = KTYP6
         JTYP2(1)  = KTYP7
         ISYMP1(1) = 1
         DELP1(1)  = .FALSE.
         DELP2(1)  = .FALSE.
         SCD(1)    = 1.0D0
         JTYP1(2)  = KTYP3
         JTYP2(2)  = KTYP9
         ISYMP1(2) = 1
         DO ISYCHO = 1,NSYM
            NCHP12(ISYCHO,2) = NUMCHO(ISYCHO)
         ENDDO
         DELP1(2)  = .FALSE.
         DELP2(2)  = .FALSE.
         SCD(2)    = 1.0D0
      ELSE
         WRITE(LUPRI,'(//,5X,A,A,A,A)')
     &   SECNAM,': ERROR: LISTR = ',LISTR(1:3),' not implemented!'
         CALL QUIT('Error in '//SECNAM)
      ENDIF

C     Allocation.
C     -----------

      KNORM  = 1
      KFOCKD = KNORM  + NPDEN
      KFIA   = KFOCKD + NORBTS
      KLAMDP = KFIA   + NT1AM(1)
      KLAMDH = KLAMDP + NGLMDT(1)
      KTBAR  = KLAMDH + NGLMDT(1)
      KEND0  = KTBAR  + NT1AM(1)
      LWRK0  = LWORK  - KEND0 + 1

      KT1AM = KTBAR

      KSCR0 = KFIA
      LSCR0 = LWORK - KSCR0 + 1

      IF (LWRK0 .LT. 0) THEN
         CALL QUIT('Insufficient memory in '//SECNAM//' [0]')
      ENDIF

C     Initializations.
C     ----------------

      CALL DZERO(WORK(KNORM),NPDEN)
      TBNRM = ZERO

C     Open FILDEN.
C     ------------

      LUDEN = 97
      CALL WOPEN2(LUDEN,FILDEN,64,0)

C     Read canonical orbital energies.
C     --------------------------------

      CALL CHO_RDSIR(DUMMY,DUMMY,WORK(KFOCKD),DUMMY,WORK(KSCR0),LSCR0,
     &               .FALSE.,.TRUE.,.FALSE.)

C     Make sure that the 0th order doubles amplitudes have
C     been Cholesky decomposed.
C     ----------------------------------------------------

      IF (.NOT. CHOT2C) THEN
         CALL CC_CHOCIM1(WORK(KFOCKD),DUMMY,DUMMY,WORK(KSCR0),LSCR0,0,0)
      ENDIF

C     Read F(ia).
C     -----------

      CALL ONEL_OP(-1,3,LUFIA)
      CALL CHO_MOREAD(WORK(KFIA),NT1AM(1),1,1,LUFIA)
      CALL ONEL_OP(1,3,LUFIA)

C     Calculate Lambda matrices.
C     --------------------------

      IOPT = 1
      CALL CC_RDRSP('R0 ',1,1,IOPT,MODEL,WORK(KT1AM),DUMMY)
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),
     &            WORK(KEND0),LWRK0)

      IF (LOCDBG) THEN
         T1NM2 = DDOT(NT1AM(1),WORK(KT1AM),1,WORK(KT1AM),1)
         T1NRM = DSQRT(T1NM2)
         WRITE(LUPRI,*)
     &   SECNAM,' after reading T1AM vector'
         WRITE(LUPRI,*) '   Square norm T1AM: ',T1NM2
         WRITE(LUPRI,*) '   2-norm      T1AM: ',T1NRM
         CALL FLSHFO(LUPRI)
      ENDIF

C     Read 0th order multipliers.
C     ---------------------------

      IOPT = 1
      CALL CC_RDRSP('L0 ',1,1,IOPT,MODEL,WORK(KTBAR),DUMMY)

      IF (LOCDBG) THEN
         T1NM2 = DDOT(NT1AM(1),WORK(KTBAR),1,WORK(KTBAR),1)
         T1NRM = DSQRT(T1NM2)
         WRITE(LUPRI,*)
     &   SECNAM,' after reading TBAR1 vector'
         WRITE(LUPRI,*) '   Square norm TBAR1: ',T1NM2
         WRITE(LUPRI,*) '   2-norm      TBAR1: ',T1NRM
         CALL FLSHFO(LUPRI)
      ENDIF

C     Calculate Ltilde(ia,J) Cholesky vectors.
C     ----------------------------------------

      DTIME = SECOND()
      CALL CC_CHOTG(WORK(KLAMDP),WORK(KLAMDH),WORK(KTBAR),
     &              WORK(KEND0),LWRK0,1,1,-1)
      DTIME  = SECOND() - DTIME
      TIMMOL = TIMMOL   + DTIME

C     Add: L(ia,J) + Ltilde(ia,J), store in KTYP8 files.
C     --------------------------------------------------

      DTIME = SECOND()
      CALL CC_CHOADD(KTYP1,KTYP8,1,NUMCHO,WORK(KEND0),LWRK0,ONE,ONE)
      DTIME  = SECOND() - DTIME
      TIMMOL = TIMMOL   + DTIME

C     Allocate space for the largest possible density blocks.
C     Find also largest set of amplitudes.
C     -------------------------------------------------------

      ISYMAX = 0
      XMAXAM = -1.0D8
      MAXDEN = -999
      DO IVEC = 1,NPDEN
         ILSTNR = IPDEN(2,IVEC)
         ISYVEC = ILSTSYM(LISTR,ILSTNR)
         NDNTST = NT1AM(ISYVEC) + NMATIJ(ISYVEC) + NMATAB(ISYVEC)
         MAXDEN = MAX(MAXDEN,NDNTST)
         IF (XMAXAM .LT. XT2SQ(ISYVEC)) THEN
            ISYMAX = ISYVEC
            XMAXAM = XT2SQ(ISYVEC)
         ENDIF
      ENDDO
      IF (MAXDEN .LE. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,I10)')
     &   SECNAM,': ERROR: MAXDEN = ',MAXDEN
         CALL QUIT('Error in '//SECNAM)
      ENDIF
      IF ((ISYMAX.LE.0) .OR. (ISYMAX.GT.NSYM)) THEN
         WRITE(LUPRI,'(//,5X,A,A,I10)')
     &   SECNAM,': ERROR: ISYMAX = ',ISYMAX
         CALL QUIT('Error in '//SECNAM)
      ENDIF
      IF (XMAXAM .LE. ZERO) THEN
         WRITE(LUPRI,'(//,5X,A,A,I10)')
     &   SECNAM,': ERROR: XMAXAM = ',XMAXAM
         CALL QUIT('Error in '//SECNAM)
      ENDIF

      KDEN  = KEND0
      KEND1 = KDEN  + MAXDEN
      LWRK1 = LWORK - KEND1 + 1

      IF (LWRK1 .LE. 1) THEN
         CALL QUIT('Insufficient memory in '//SECNAM//' [1]')
      ENDIF

C     Split remaining work space in two equal parts:
C     one for 0th order multipliers, one for LISTR amplitudes.
C     --------------------------------------------------------

      LWRK = LWRK1/2

C     Start batch loop over a index.
C     ------------------------------

      IBATA = 0
      ABEG  = 1
      NUMA  = 0
  100 CONTINUE

         ABEG = ABEG + NUMA
         IF (ABEG .GT. NVIRT) GO TO 200

C        Time this batch.
C        ----------------

         TABAT = SECOND()

C        Update batch counter.
C        ---------------------

         IBATA = IBATA + 1

C        Find out how many a's can be treated, estimating
C        the number of b's along the way using ISYMAX as
C        LISTR amplitude symmetry.
C        ------------------------------------------------

         CALL CC_CYIBTB(LVIRA,ABEG,NUMA,LWRK,ISYMAX)

C        Check for errors.
C        -----------------

         AEND = ABEG + NUMA - 1
         IF (AEND .GT. NVIRT) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Batch error in ',SECNAM,' !'
            WRITE(LUPRI,'(5X,A)')
     &      '- dumping presumably most relevant info:'
            WRITE(LUPRI,'(5X,A,I10)')
     &      'a-batch number    : ',IBATA
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,I10)')
     &      'First A           : ',ABEG,
     &      'Last  A           : ',AEND,
     &      'Number of virtuals: ',NVIRT
            WRITE(LUPRI,'(5X,A,I10,/)')
     &      'Memory for batch  : ',LWRK
            CALL QUIT('Batch error in '//SECNAM)
         ENDIF

         IF (NUMA .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory for a-batch in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/)')
     &      'Memory for batch: ',LWRK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Get symmetries of first and last a in batch.
C        --------------------------------------------

         ISABEG = ISVI(ABEG)
         ISAEND = ISVI(AEND)

C        Set up index arrays.
C        --------------------

         CALL GET_INDVIR(ABEG,AEND,LVIRA,IOFA1,IX1AMA,NX1AMA)

         DO ISYAIJ = 1,NSYM
            LAIJ(ISYAIJ) = 0
            DO ISYMAI = 1,NSYM
               ISYMJ = MULD2H(ISYMAI,ISYAIJ)
               LAIJ(ISYAIJ) = LAIJ(ISYAIJ)
     &                      + NRHF(ISYMJ)*NX1AMA(ISYMAI)
            ENDDO
         ENDDO

C        Calculate 0th order multipliers tbar(#ai,dl) for all d
C        and store on disk (direct-access scratch files).
C        Note that the full work space (LWRK0) is passed here.
C        ------------------------------------------------------

         DTMM = SECOND()
         CALL CC_CHOMUL0(WORK(KFOCKD),WORK(KTBAR),WORK(KFIA),
     &                   WORK(KEND0),LWRK0,LAIJ,KTYP1,KTYP8,
     &                   NUMCHO,.FALSE.,.FALSE.,FILMUL,TBNRM,KLAST)
         DTMM = SECOND() - DTMM
         TIMMUL = TIMMUL + DTMM

         MUSED  = KEND0 + KLAST - 1
         MXMEMT = MAX(MXMEMT,MUSED)

C        Print.
C        ------

         IF (IPRINT .GE. INFO) THEN
            XUSED = (D100*MUSED)/XLWORK
            WRITE(LUPRI,'(5X,A,I6,A,/,5X,A)')
     &      'a-batch number ',IBATA,':',
     &      '======================'
            WRITE(LUPRI,'(5X,A,I6)')
     &      '#a   : ',NUMA
            WRITE(LUPRI,'(5X,A,I6,A,I1,/,5X,A,I6,A,I1)')
     &      'First a: ',IOFA1(ISABEG),'   symmetry: ',ISABEG,
     &      'Last  a: ',IOFA1(ISAEND)+LVIRA(ISAEND)-1,'   symmetry: ',
     &                                                  ISAEND
            WRITE(LUPRI,'(5X,A,F10.2,/,5X,A,1P,D22.15)')
     &      'Time for 0th order multipliers, tbar(#ai,dl): ',DTMM,
     &      'Accumulated norm of packed multipliers: ',DSQRT(TBNRM)
            WRITE(LUPRI,'(5X,A,F7.2,A,/)')
     &      'Memory used for assembling multipliers: ',XUSED,'%'
         ENDIF

C        Initialize address on FILDEN.
C        -----------------------------

         ADR = 1

C        Start loop over LISTR vectors.
C        ------------------------------

         DO IVEC = 1,NPDEN

C           Get vector/operator/frequency info.
C           -----------------------------------

            ILSTNR = IPDEN(2,IVEC)
            ISYDEN = ILSTSYM(LISTR,ILSTNR)
            LABEL  = LRTLBL(ILSTNR)
            FREQ   = FRQLRT(ILSTNR)

C           Read AO property integrals corresponding to operator LABEL.
C           Check info. Use work space from KEND0.
C           -----------------------------------------------------------

            KXIJ  = KEND0
            KXAB  = KXIJ  + NMATIJ(ISYDEN)
            KXAO  = KXAB  + NMATAB(ISYDEN)
            KEND2 = KXAO  + N2BST(ISYDEN)
            LWRK2 = LWORK - KEND2 + 1

            IF (LWRK2 .LT. 0) THEN
               CALL QUIT('Insufficient memory in '//SECNAM//' [2]')
            ENDIF 
            
            IERR = -1
            CALL CCPRPAO(LABEL,LONLYAO,WORK(KXAO),IRREP,IREAL,IERR,
     &                   WORK(KEND2),LWRK2)

            IF (IERR .GT. 0) THEN
               WRITE(LUPRI,'(//,5X,A,A)')
     &         SECNAM,': I/O error ocurred in CCPRPAO'
               WRITE(LUPRI,'(5X,A,I10,A,/,5X,A,2X,A8,/,5X,A,I10)')
     &         'Operator number  : ',ILSTNR,' (in LRTLBL list)',
     &         'Operator label   : ',LABEL,
     &         'Operator symmetry: ',ISYDEN
               WRITE(LUPRI,'(5X,A,I10,/)')
     &         'CCPRPAO returned IERR = ',IERR
               CALL QUIT('I/O error in CCPRPAO ('//SECNAM//')')
            ELSE IF (IERR .LT. 0) THEN
               WRITE(LUPRI,'(//,5X,A,A)')
     &         SECNAM,': CCPRPAO failed to determine operator symmetry'
               WRITE(LUPRI,'(5X,A,I10,A,/,5X,A,2X,A8,/,5X,A,I10)')
     &         'Operator number  : ',ILSTNR,' (in LRTLBL list)',
     &         'Operator label   : ',LABEL,
     &         'Operator symmetry: ',ISYDEN
               WRITE(LUPRI,'(5X,A,I10,/)')
     &         'CCPRPAO returned IERR = ',IERR
               CALL QUIT('Irrep. error in CCPRPAO ('//SECNAM//')')
            ENDIF

            IF (IRREP .NE. ISYDEN) THEN
               WRITE(LUPRI,'(//,5X,A,A)')
     &         'Symmetry mismatch in ',SECNAM
               WRITE(LUPRI,'(5X,A,I10,A,/,5X,A,2X,A8)')
     &         'Operator number  : ',ILSTNR,' (in LRTLBL list)',
     &         'Operator label   : ',LABEL
               WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &         'CCPRPAO returned: ',IRREP,
     &         '- expected      : ',ISYDEN
               CALL QUIT('Symmetry mismatch in '//SECNAM)
            ENDIF

C           Transform integrals: occupied and virtual blocks only.
C           ------------------------------------------------------

            CALL CC_CHOTRFOV(WORK(KLAMDP),WORK(KLAMDH),WORK(KXAO),
     &                       WORK(KXIJ),WORK(KXAB),WORK(KEND2),LWRK2,
     &                       1,1,ISYDEN)

C           Reset memory: AO integrals no longer needed.
C           --------------------------------------------

            KEND2 = KXAO
            LWRK2 = LWORK - KEND2 + 1

C           Transform Cholesky vectors corresponding to 0th order
C           amplitudes, incorporating the minus sign. Store result
C           in KTYP7 files.
C           ------------------------------------------------------

            DTIME = SECOND()
            CALL CC_CHOTRPRP(WORK(KXIJ),WORK(KXAB),WORK(KEND2),LWRK2,
     &                       ISYDEN,KTYP6,KTYP7,1,NCHMOC)
            DTIME = SECOND() - DTIME
            TIMMOR = TIMMOR  + DTIME

C           Re-allocate from KEND0.
C           -----------------------

            KLAM2 = KEND0
            KR1AM = KLAM2 + NGLMDT(ISYDEN)
            KEND3 = KR1AM + NT1AM(ISYDEN)
            LWRK3 = LWORK - KEND3 + 1

            IF (LWRK3 .LT. 0) THEN
               CALL QUIT('Insufficient memory in '//SECNAM//' [3]')
            ENDIF

C           Read R1 amplitudes.
C           -------------------

            IOPT = 1
            CALL CC_RDRSP(LISTR,ILSTNR,ISYDEN,IOPT,MODEL,WORK(KR1AM),
     &                    DUMMY)

            IF (LOCDBG) THEN
               R1NM2 = DDOT(NT1AM(ISYDEN),WORK(KR1AM),1,WORK(KR1AM),1)
               R1NRM = DSQRT(R1NM2)
               WRITE(LUPRI,*)
     &         SECNAM,' after reading R vector number: ',ILSTNR
               WRITE(LUPRI,*) '   symmetry        : ',ISYDEN
               WRITE(LUPRI,*) '   assoc. operator : ',LABEL
               WRITE(LUPRI,*) '   assoc. frequency: ',FREQ
               WRITE(LUPRI,*) '   Square norm R1AM: ',R1NM2
               WRITE(LUPRI,*) '   2-norm      R1AM: ',R1NRM
               CALL FLSHFO(LUPRI)
            ENDIF

C           Calculate perturbed lambda matrix.
C           ----------------------------------

            CALL CC_LAMPT(WORK(KLAMDP),WORK(KLAMDH),WORK(KR1AM),
     &                    WORK(KLAM2),1,ISYDEN,1)

C           Reset memory: R1 no longer needed.
C           ----------------------------------

            KEND3 = KR1AM
            LWRK3 = LWORK - KEND3 + 1

C           Transform Cholesky vectors: Lbar(ai,J) stored in KTYP9 files.
C           -------------------------------------------------------------

            DTIME = SECOND()
            CALL CC_CHOTG1(WORK(KLAMDH),WORK(KLAMDP),WORK(KLAM2),
     &                     WORK(KEND3),LWRK3,1,ISYDEN,KTYP9)
            DTIME = SECOND() - DTIME
            TIMMOR = TIMMOR  + DTIME

C           Make allocation for MO density matrix.
C           NB!! Ordering is important!
C           --------------------------------------

            LDEN = NMATIJ(ISYDEN) + NT1AM(ISYDEN) + NMATAB(ISYDEN)

            KDNJK = KEND0
            KDNJB = KDNJK + NMATIJ(ISYDEN)
            KDNCB = KDNJB + NT1AM(ISYDEN)
            KEND4 = KDNCB + NMATAB(ISYDEN)
            LWRK4 = LWORK - KEND4 + 1

            IF (LWRK4 .LE. 1) THEN  ! Should never happen....
               CALL QUIT('Insufficient memory in '//SECNAM//' [4]')
            ENDIF

C           Read appropriate density matrix.
C           --------------------------------

            CALL GETWA2(LUDEN,FILDEN,WORK(KDNJK),ADR,LDEN)

            IF (IPRINT .GE. INFO) THEN
               WRITE(LUPRI,'(7X,A,A,A,I10,A)')
     &         'Amplitudes: ',LISTR(1:3),' number',ILSTNR,':'
               WRITE(LUPRI,'(7X,A,A,/,7X,A,A)')
     &         '    #b   First     Last    Memory             CPU Time',
     &         ' (Seconds)',
     &         '       (b,sym.b) (b,sym.b) Usage (%)',
     &         '  Amplitudes  Density       Total'
               WRITE(LUPRI,'(7X,A,A)')
     &         '------------------------------------------------------',
     &         '---------------'
               CALL FLSHFO(LUPRI)
            ENDIF

C           Start batch loop over b.
C           ------------------------

            LSCR1 = LWRK4/2
            IBATB = 0
            BBEG  = 1
            NUMB  = 0
  300       CONTINUE

               BBEG = BBEG + NUMB
               IF (BBEG .GT. NVIRT) GO TO 400

C              Time this batch.
C              ----------------

               TBBAT = SECOND()

C              Update batch counter.
C              ---------------------

               IBATB = IBATB + 1

C              Find out how many b's can be treated.
C              -------------------------------------

               CALL CC_CYIBTA(LAIJ,NX1AMA,LVIRB,BBEG,NUMB,MEM,LSCR1,
     &                        ISYDEN)

C              Check for errors.
C              -----------------

               BEND = BBEG + NUMB - 1
               IF (BEND .GT. NVIRT) THEN
                  WRITE(LUPRI,'(//,5X,A,A,A)')
     &            'Batch error in ',SECNAM,' !'
                  WRITE(LUPRI,'(5X,A)')
     &            '- dumping presumably most relevant info:'
                  WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     &            'A-batch number    : ',IBATA,
     &            'B-batch number    : ',IBATB
                  WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     &            'First B           : ',BBEG,
     &            'Last  B           : ',BEND
                  WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,I10)')
     &            'First A           : ',ABEG,
     &            'Last  A           : ',AEND,
     &            'Number of virtuals: ',NVIRT
                  WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     &            'Available memory  : ',LSCR1,
     &            'Needed    memory  : ',MEM
                  WRITE(LUPRI,'(5X,A,I10,A,/,5X,A,2X,A8,/,5X,A,I10,/)')
     &            'Operator number  : ',ILSTNR,' (in LRTLBL list)',
     &            'Operator label   : ',LABEL,
     &            'Operator symmetry: ',ISYDEN
                  CALL QUIT('Batch error in '//SECNAM)
               ENDIF

               IF (NUMB .LE. 0) THEN
                  WRITE(LUPRI,'(//,5X,A,A)')
     &            'Insufficient memory for b-batch in ',SECNAM
                  WRITE(LUPRI,'(5X,A,I10,A,/,5X,A,2X,A8,/,5X,A,I10,/)')
     &            'Operator number  : ',ILSTNR,' (in LRTLBL list)',
     &            'Operator label   : ',LABEL,
     &            'Operator symmetry: ',ISYDEN
                  WRITE(LUPRI,'(5X,A,I10,/)')
     &            'Available memory: ',LSCR1
                  CALL QUIT('Insufficient memory in '//SECNAM)
               ENDIF

C              Get symmetries of first and last b in batch.
C              --------------------------------------------

               ISBBEG = ISVI(BBEG)
               ISBEND = ISVI(BEND)

C              Set up index arrays.
C              --------------------

               CALL GET_INDVIR(BBEG,BEND,LVIRB,IOFB1,IX1AMB,NX1AMB)

               CALL IZERO(IX2SQ,64)
               ICOUNT = 0
               DO ISYMBJ = 1,NSYM
                  ISYMAI = MULD2H(ISYMBJ,ISYDEN)
                  IX2SQ(ISYMAI,ISYMBJ) = ICOUNT
                  ICOUNT = ICOUNT + NX1AMA(ISYMAI)*NX1AMB(ISYMBJ)
               ENDDO
               NX2SQ = ICOUNT

               IF (NX2SQ .NE. MEM) THEN
                  WRITE(LUPRI,'(//,5X,A,A)')
     &            'Error detected in ',SECNAM
                  WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,/)')
     &            'NX2SQ is calculated to be ',NX2SQ,
     &            'MEM determines alloc. as  ',MEM,
     &            'These numbers must be equal....'
                  CALL QUIT('Error in '//SECNAM)
               ENDIF

C              Allocation for LISTR amplitude array.
C              -------------------------------------

               KR2AM = KEND4
               KEND5 = KR2AM + NX2SQ
               LWRK5 = LWORK - KEND5 + 1

C              Calculate amplitude batch.
C              --------------------------

               ISYMP2(1) = ISYDEN
               ISYMP2(2) = ISYDEN
               DO ISYCHO = 1,NSYM
                  NCHP12(ISYCHO,1) = NCHMOC(ISYCHO)
               ENDDO

               DTMR = SECOND()
               CALL CC_CAA(JTYP1,ISYMP1,JTYP2,ISYMP2,NCHP12,FAC1,NUMP12,
     &                     DUMMY,IDUMMY,DUMMY,IDUMMY,FAC2,NUMUV,
     &                     SCD,SCDG,IOPTDN,WORK(KFOCKD),FREQ,IOPTCE,
     &                     WORK(KR2AM),WORK(KEND5),LWRK5,X2NRM,X2CNM,
     &                     DELP1,DELP2,KLAST)
               DTMR = SECOND() - DTMR
               TIMRAM = TIMRAM + DTMR

               KOFF = KNORM + IVEC - 1
               WORK(KOFF) = WORK(KOFF) + X2NRM

               MUSED  = KEND5 + KLAST - 1
               MXMEML = MAX(MXMEML,MUSED)

C              Calculate contributions to density matrix blocks.
C              NB!!! The amplitudes are replaced here with 2CME!
C              -------------------------------------------------

               DTMD = SECOND()
               CALL CC_CHOR1D1(WORK(KDNJK),WORK(KDNJB),WORK(KDNCB),
     &                         WORK(KTBAR),WORK(KR2AM),
     &                         WORK(KEND5),LWRK5,1,ISYDEN,LAIJ,FILMUL,
     &                         KLAST)
               DTMD = SECOND() - DTMD
               TIMDEN = TIMDEN + DTMD

               MUSED  = KEND5 + KLAST - 1
               MXMEML = MAX(MXMEML,MUSED)

C              Print info from this b-batch.
C              -----------------------------

               IF (IPRINT .GE. INFO) THEN
                  TBBAT = SECOND() - TBBAT
                  XMUSE = (D100*MXMEML)/XLWORK
                  LBFST = IOFB1(ISBBEG)
                  LBLST = IOFB1(ISBEND) + LVIRB(ISBEND) - 1
                  WRITE(LUPRI,1) NUMB,LBFST,ISBBEG,LBLST,ISBEND,XMUSE,
     &                           DTMR,DTMD,TBBAT
    1             FORMAT(7X,I6,1X,I6,1X,I1,2X,I6,1X,I1,2X,F7.2,3X,F10.2,
     &                   1X,F10.2,1X,F10.2)
                  CALL FLSHFO(LUPRI)
               ENDIF

C              Update memory info.
C              -------------------

               MXMEMT = MAX(MXMEMT,MXMEML)
               MXMEML = 0

C              Go to next b-batch.
C              -------------------

               GO TO 300

C           Exit point for b-batch loop.
C           ----------------------------

  400       CONTINUE

C           Write density blocks to disk and update address.
C           ------------------------------------------------

            CALL PUTWA2(LUDEN,FILDEN,WORK(KDNJK),ADR,LDEN)
            ADR = ADR + LDEN

C           Print info for this amplitude.
C           ------------------------------

            IF (IPRINT .GE. INFO) THEN
               WRITE(LUPRI,'(7X,A,A)')
     &         '------------------------------------------------------',
     &         '---------------'
               X2NRM = DSQRT(WORK(KNORM+IVEC-1))
               WRITE(LUPRI,'(7X,A,1P,D22.15,/)')
     &         'Accumulated packed amplitude norm: ',X2NRM
               CALL FLSHFO(LUPRI)
            ENDIF

         ENDDO  ! end of amplitude/operator loop

C        Print info from this a-batch.
C        -----------------------------

         IF (IPRINT .GE. INFO) THEN
            TABAT  = SECOND() - TABAT
            DELTM  = SECOND() - TIMTOT
            XVTOT  = ONE*NVIRT
            PCTDON = (D100*AEND)/XVTOT
            WRITE(LUPRI,'(5X,A,F10.2,A,/,5X,A,F10.2,A,A,F6.2,A,/)')
     &      'Total time for this a-batch: ',TABAT,' seconds',
     &      'Total time used so far     : ',DELTM,' seconds',
     &      ' (',PCTDON,'% done)'
            CALL FLSHFO(LUPRI)
         ENDIF

C        Go to next a-batch.
C        -------------------

         GO TO 100

C     Exit point for a-batch loop.
C     ----------------------------

  200 CONTINUE

C     Clean up Cholesky vector files KTYP7, KTYP8, and KTYP9.
C     -------------------------------------------------------

      DO ISYCHO = 1,NSYM
         CALL CHO_MOP(-1,KTYP7,ISYCHO,LUCH,1,1)
         CALL CHO_MOP(0,KTYP7,ISYCHO,LUCH,1,1)
         CALL CHO_MOP(-1,KTYP8,ISYCHO,LUCH,1,1)
         CALL CHO_MOP(0,KTYP8,ISYCHO,LUCH,1,1)
         CALL CHO_MOP(-1,KTYP9,ISYCHO,LUCH,1,1)
         CALL CHO_MOP(0,KTYP9,ISYCHO,LUCH,1,1)
      ENDDO

C     Backtransform densities to AO basis and write to disk.
C     ------------------------------------------------------

      ADR = 1

      TIMBCK = SECOND()
      DO IDEN = 1,NPDEN

         ILSTNR = IPDEN(2,IDEN)
         ISYDEN = ILSTSYM(LISTR,ILSTNR)
c        LABEL  = LRTLBL(ILSTNR)
c        FREQ   = FRQLRT(ILSTNR)

         KDNJK = KEND0
         KDNJB = KDNJK + NMATIJ(ISYDEN)
         KDNCB = KDNJB + NT1AM(ISYDEN)
         KDNBJ = KDNCB + NMATAB(ISYDEN)
         KDNAO = KDNBJ + NT1AM(ISYDEN)
         KENDT = KDNAO + N2BST(ISYDEN)
         LWRKT = LWORK - KENDT + 1

         IF (LWRKT .LT. 0) THEN
            CALL QUIT('Insufficient memory in '//SECNAM//' [T]')
         ENDIF

         LDEN = NMATIJ(ISYDEN) + NT1AM(ISYDEN) + NMATAB(ISYDEN)
         CALL GETWA2(LUDEN,FILDEN,WORK(KDNJK),ADR,LDEN)
         CALL DZERO(WORK(KDNBJ),NT1AM(ISYDEN))
         ADR = ADR + LDEN

         CALL DZERO(WORK(KDNAO),N2BST(ISYDEN))
         CALL CC_DENAO(WORK(KDNAO),ISYDEN,WORK(KDNBJ),WORK(KDNCB),
     &                 WORK(KDNJK),WORK(KDNJB),ISYDEN,WORK(KLAMDP),
     &                 1,WORK(KLAMDH),1,WORK(KENDT),LWRKT)

         IPDEN(3,IDEN) = ILSTNR
         CALL CC_WRRSPD(LISTD,ILSTNR,ISYDEN,MODEL,.FALSE.,WORK(KDNAO),
     &                  WORK(KENDT),LWRKT)

      ENDDO
      TIMBCK = SECOND() - TIMBCK
      TIMDEN = TIMDEN   + TIMBCK

C     Close FILDEN, delete if requested.
C     ----------------------------------

      IF (DELDEN) THEN
         CALL WCLOSE2(LUDEN,FILDEN,'DELETE')
      ELSE
         CALL WCLOSE2(LUDEN,FILDEN,'KEEP')
      ENDIF

C     Print.
C     ------

      IF (IPRINT .GE. INFO) THEN
         TIMTOT = SECOND() - TIMTOT
         CALL HEADER('Timing of '//SECNAM,-1)
         IF (.NOT. SKPSNG) THEN
            WRITE(LUPRI,'(5X,A,F10.2,A)')
     &      'Time used for singles contributions   : ',TIMSNG,' seconds'
         ENDIF
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for 0th order multipliers   : ',TIMMUL,' seconds'
         WRITE(LUPRI,'(5X,A,A,A,F10.2,A)')
     &   'Time used for ',LISTR(1:3),' amplitudes          : ',TIMRAM,
     &   ' seconds'
         WRITE(LUPRI,'(5X,A,A,A,F10.2,A)')
     &   'Time used for ',LISTD(1:3),' AO densities        : ',TIMDEN,
     &   ' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for Chol. preparation, left : ',TIMMOL,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for Chol. preparation, right: ',TIMMOR,' seconds'
         WRITE(LUPRI,'(5X,A)')
     &   '----------------------------------------------------------'
         WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &   'Time used in total                    : ',TIMTOT,' seconds'
         CALL FLSHFO(LUPRI)
      ENDIF

      RETURN
      END
C  /* Deck cc_chomul0 */
      SUBROUTINE CC_CHOMUL0(FOCKD,TBAR1,FIA,WORK,LWORK,LAIJ,
     &                      ITY1,ITY2,NUMCHO,DEL1,DEL2,FILMUL,
     &                      TBNRM,MEMUSE)
C
C     Thomas Bondo Pedersen, April 2003.
C
C     Purpose: Calculate batch of 0th order multipliers tbar(#ai,bj) and
C              store on direct-access files FILMUL(*), each record
C              containing an #ai column of tbar.
C              Information about the a batch is taken from ciarc.h (except
C              for the LAIJ(*) array).
C              It is implicitly assumed that the necessary Cholesky
C              vectors are available on disk in files specified through
C              ITY1 and ITY2 (number of vectors in each symmetry in NUMCHO).
C              The Cholesky vectors are assumed total symmetric.
C
#include <implicit.h>
      DIMENSION FOCKD(*), TBAR1(*), FIA(*), WORK(LWORK)
      INTEGER   LAIJ(8), NUMCHO(8)
      LOGICAL   DEL1, DEL2
      CHARACTER*8 FILMUL(8)
#include <maxorb.h>
#include <ccisvi.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>
#include <ciarc.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CC_CHOMUL0')

      LOGICAL OLDDX

      INTEGER ABEG, AEND
      INTEGER BBEG, BEND

      PARAMETER (MXTYP = 1)
      DIMENSION SCD(MXTYP)
      INTEGER JTYP1(MXTYP), JTYP2(MXTYP), ISYMP1(MXTYP), ISYMP2(MXTYP)
      INTEGER NCHP12(8,MXTYP)
      LOGICAL DELP1(MXTYP), DELP2(MXTYP)

      PARAMETER (MXNMUV = 1)
      INTEGER ISYMU(MXNMUV), ISYMV(MXNMUV)

      PARAMETER (INFO = 3)
      PARAMETER (IOPEN = -1, IKEEP = 1)
      PARAMETER (ZERO = 0.00D0, ONE = 1.00D0, TWO = 2.00D0)
      PARAMETER (WTOMB = 7.62939453125D-6, D100 = 100.00D0)

C     Start timing.
C     -------------

      TIMT = SECOND()
      TIMA = ZERO
      TIMW = ZERO

C     Initializations.
C     ----------------

      MEMUSE = 0
      MXMEMT = -999
      XLWORK = ONE*LWORK

C     Set up info for amplitude assembler
C     -----------------------------------

      NUMP12 = 1
      NUMUV  = 1

      FAC1 = 1.0D0
      FAC2 = 1.0D0
      SCDG = 1.0D0

      IOPTDN = 1
      IOPTCE = 1

      NUMP12 = 1
      NUMUV  = 1

      JTYP1(1)  = ITY1
      JTYP2(1)  = ITY2
      ISYMP1(1) = 1
      ISYMP2(1) = 1
      DO ISYCHO = 1,NSYM
         NCHP12(ISYCHO,1) = NUMCHO(ISYCHO)
      ENDDO

      DELP1(1) = DEL1
      DELP2(1) = DEL2

      SCD(1) = 1.0D0

      ISYMU(1) = 1
      ISYMV(1) = 1

      FRQL0 = 0.0D0

C     For the purpose of the memory split, set control parameters
C     in ciarc.h.
C     -----------------------------------------------------------

      DO ISYCHO = 1,NSYM
         NTOVEC(ISYCHO) = NCHP12(ISYCHO,1)
      ENDDO
      ISYCH1 = ISYMP1(1)
      ISYCH2 = ISYMP2(1)
      ITYP1  = JTYP1(1)
      ITYP2  = JTYP2(1)

      LIDEN  = ITYP1 .EQ. ITYP2
      SYMTRZ = ITYP1 .NE. ITYP2
      CIAMIO = .FALSE.
      GETMNM = .FALSE.
      INXINT = .TRUE.

C     Split memory: one part for amplitudes, one for Cholesky vectors.
C     ----------------------------------------------------------------

      CALL CC_ML0MSP(LWORK,LWRK,LAIJ)
      LEFT = LWORK - LWRK

C     Check.
C     ------

      IF ((LWRK.LE.0) .OR. (LEFT.LE.0)) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,I10,/)')
     &   'Total memory available                  : ',LWORK,
     &   'Memory reserved for integrals/amplitudes: ',LWRK,
     &   'Memory reserved for Cholesky vectors    : ',LEFT
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Print.
C     ------

      IF (IPRINT .GE. INFO) THEN
         CALL AROUND('Output from '//SECNAM)
         XMB = XLWORK*WTOMB
         WRITE(LUPRI,'(5X,A,I10,A,F7.1,A)')
     &   'Available memory: ',LWORK,' words (',XMB,' Mb).'
         WRITE(LUPRI,'(5X,A,1X,I10,A)')
     &   'Maximum memory for integral intermediates: ',LWRK,' words.'
         CALL INV_BTCH(NUMA,ABEG,ISABEG,AEND,ISAEND,IOFA1,LVIRA)
         WRITE(LUPRI,'(5X,A,I10)')
     &   '#a: ',NUMA
         LAFST = IOFA1(ISABEG)
         LALST = IOFA1(ISAEND) + LVIRA(ISAEND) - 1
         WRITE(LUPRI,'(5X,A,I10,A,I2)')
     &   'First a: ',LAFST,' Symmetry:',ISABEG
         WRITE(LUPRI,'(5X,A,I10,A,I2)')
     &   'Last  a: ',LALST,' Symmetry:',ISAEND
         WRITE(LUPRI,'(A)') ' '
         WRITE(LUPRI,'(5X,A,A,/,5X,A,A)')
     &   '   #b        First           Last         Mem. Usage',
     &   '        Time   ',
     &   '          (b, sym. b)     (b, sym. b)        (%)    ',
     &   '      (seconds)'
         WRITE(LUPRI,'(5X,A,A)')
     &   '----------------------------------------------------',
     &   '---------------'
         CALL FLSHFO(LUPRI)
      ENDIF

C     Set integral symmetry.
C     ----------------------

      ISYINT = MULD2H(ISYCH1,ISYCH2)
      IF (ISYINT .NE. 1) CALL QUIT('Symmetry error in '//SECNAM)

C     Start batch loop over b.
C     ------------------------

      IBATB = 0
      BBEG  = 1
      NUMB  = 0
  100 CONTINUE

         BBEG = BBEG + NUMB
         IF (BBEG .GT. NVIRT) GO TO 200

C        Time this batch.
C        ----------------

         TBBAT = SECOND()

C        Update batch counter.
C        ---------------------

         IBATB = IBATB + 1

C        Find out how many b's can be treated.
C        -------------------------------------

         CALL CC_CYIBTA(LAIJ,NX1AMA,LVIRB,BBEG,NUMB,MEM,LWRK,ISYINT)

C        Check for errors.
C        -----------------

         BEND = BBEG + NUMB - 1
         IF (BEND .GT. NVIRT) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Batch error in ',SECNAM,' !'
            WRITE(LUPRI,'(5X,A)')
     &      '- dumping presumably most relevant info:'
            WRITE(LUPRI,'(5X,A,I10)')
     &      'B-batch number    : ',IBATB
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     &      'First B           : ',BBEG,
     &      'Last  B           : ',BEND
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Available memory  : ',LWRK,
     &      'Needed    memory  : ',MEM
            CALL QUIT('Batch error in '//SECNAM)
         ENDIF

         IF (NUMB .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory for B-batch in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10)')
     &      'Available memory: ',LWRK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Get symmetries of first and last b in batch.
C        --------------------------------------------

         ISBBEG = ISVI(BBEG)
         ISBEND = ISVI(BEND)

C        Set up index arrays.
C        --------------------

         CALL GET_INDVIR(BBEG,BEND,LVIRB,IOFB1,IX1AMB,NX1AMB)

         CALL IZERO(IX2SQ,64)
         ICOUNT = 0
         DO ISYMBJ = 1,NSYM
            ISYMAI = MULD2H(ISYMBJ,ISYINT)
            IX2SQ(ISYMAI,ISYMBJ) = ICOUNT
            ICOUNT = ICOUNT + NX1AMA(ISYMAI)*NX1AMB(ISYMBJ)
         ENDDO
         NX2SQ = ICOUNT

         IF (NX2SQ .NE. MEM) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Error detected in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,/)')
     &      'NX2SQ is calculated to be ',NX2SQ,
     &      'MEM determines alloc. as  ',MEM,
     &      'These numbers must be equal....'
            CALL QUIT('Error in '//SECNAM)
         ENDIF

C        Allocation for amplitudes.
C        --------------------------

         KTBAR = 1
         KEND1 = KTBAR + NX2SQ
         LWRK1 = LWORK - KEND1 + 1

C        Calculate amplitudes.
C        ---------------------

         DTIME = SECOND()
         CALL CC_CAA(JTYP1,ISYMP1,JTYP2,ISYMP2,NCHP12,FAC1,NUMP12,
     &               TBAR1,ISYMU,FIA,ISYMV,FAC2,NUMUV,
     &               SCD,SCDG,IOPTDN,FOCKD,FRQL0,IOPTCE,
     &               WORK(KTBAR),WORK(KEND1),LWRK1,X2NRM,X2CNM,
     &               DELP1,DELP2,KLAST)
         DTIME = SECOND() - DTIME
         TIMA  = TIMA     + DTIME

         TBNRM  = TBNRM + X2CNM
         MXMEML = KEND1 + KLAST - 1

C        Write amplitude batch to disk.
C        ------------------------------

         DTIME = SECOND()
         DO ISYMBJ = 1,NSYM

            ISYMAI = MULD2H(ISYMBJ,ISYINT)
            LENAI  = NX1AMA(ISYMAI)

            IF (LENAI .GT. 0) THEN

               LRLEN = 2*LENAI
               LUNIT = -1
               CALL GPOPEN(LUNIT,FILMUL(ISYMBJ),'UNKNOWN','DIRECT',
     &                     'UNFORMATTED',LRLEN,OLDDX)

               DO ISYMJ = 1,NSYM

                  ISYMB = MULD2H(ISYMJ,ISYMBJ)
                  LENB  = LVIRB(ISYMB)

                  DO J = 1,NRHF(ISYMJ)

                     NBJ1 = IT1AM(ISYMB,ISYMJ)
     &                    + NVIR(ISYMB)*(J - 1) + IOFB1(ISYMB)
                     LBJ1 = IX1AMB(ISYMB,ISYMJ)
     &                    + LVIRB(ISYMB)*(J - 1)
                     KOFF = KTBAR + IX2SQ(ISYMAI,ISYMBJ) + LENAI*LBJ1

                     CALL CHO_MOWRITE(WORK(KOFF),LENAI,LENB,NBJ1,LUNIT)

                  ENDDO

               ENDDO

               CALL GPCLOSE(LUNIT,'KEEP')

            ENDIF

         ENDDO
         DTIME = SECOND() - DTIME
         TIMW  = TIMW     + DTIME

C        Print information from this b-batch.
C        ------------------------------------

         IF (IPRINT .GE. INFO) THEN
            TBBAT = SECOND() - TBBAT
            XMUSE = (D100*MXMEML)/XLWORK
            LBFST = IOFB1(ISBBEG)
            LBLST = IOFB1(ISBEND) + LVIRB(ISBEND) - 1
         WRITE(LUPRI,'(5X,I6,4X,I6,1X,I1,8X,I6,1X,I1,9X,F6.2,8X,F10.2)')
     &      NUMB,LBFST,ISBBEG,LBLST,ISBEND,XMUSE,TBBAT
            CALL FLSHFO(LUPRI)
         ENDIF

C        Update memory info and go to next b-batch.
C        ------------------------------------------

         MXMEMT = MAX(MXMEMT,MXMEML)
         GO TO 100

C     Exit point for b-batch.
C     -----------------------

  200 CONTINUE

C     Set MEMUSE.
C     -----------

      MEMUSE = MXMEMT

C     Print.
C     ------

      IF (IPRINT .GE. INFO) THEN
         WRITE(LUPRI,'(5X,A,A)')
     &   '----------------------------------------------------',
     &   '---------------'
         WRITE(LUPRI,'(5X,A,1P,D22.15,/)')
     &   'Accumulated multiplier norm, packed: ',DSQRT(TBNRM)
         TIMT = SECOND() - TIMT
         CALL HEADER('Timing of '//SECNAM,-1)
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for multiplier assembly: ',TIMA,' seconds.'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for multiplier I/O     : ',TIMW,' seconds.'
         WRITE(LUPRI,'(5X,A)')
     &   '------------------------------------------------------'
         WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &   'Time used in total               : ',TIMT,' seconds.'
         CALL FLSHFO(LUPRI)
      ENDIF

      RETURN
      END
C  /* Deck cc_ml0msp */
      SUBROUTINE CC_ML0MSP(LWORK,LWRK,LAIJ)
C
C     Thomas Bondo Pedersen, April 2003.
C     - almost exact copy of CC_CYIMSP.
C
C     Purpose: Split memory for CC_CHOMUL0.
C              The routine should be flexible enough for use also with
C              other amplitude types than 0th order multipliers.
C              Note that this routine makes assumptions about the memory
C              allocations in CC_CIA and it is furthermore assumed that the
C              size of the a-batch is already known.
C
C     Input:
C
C        LWORK  : Work space available.
C
C        Logical control variables and # Cholesky vectors
C        in ciarc.h
C
C     Output:
C
C        LWRK   : Dimension of work space reserved for
C                 tbar(#ai,#bj) integral/amplitude batch.
C
#include <implicit.h>
      INTEGER LAIJ(8)
#include <ccorb.h>
#include <ccsdsym.h>
#include <dccorb.h>
#include <dccsdsym.h>
#include <priunit.h>
#include <ciarc.h>
#include <chocc2.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_ML0MSP')

      LOGICAL LALLA

      PARAMETER (MINDEF = 5, LFRAC = 5)
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (THREE = 3.0D0, FOUR = 4.0D0)
      PARAMETER (TINY = 1.0D-14)

C     Test LWORK.
C     ===========

      IF (LWORK .LE. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A,I10,/)')
     &   'Work space non-positive in ',SECNAM,': LWORK = ',LWORK
         CALL QUIT('Error in '//SECNAM)
      ENDIF
      XWORK  = ONE*LWORK
      ISYINT = MULD2H(ISYCH1,ISYCH2)

C     Set Cholesky weight.
C     ====================

      WCHOL = SPLITC

      IF (WCHOL .LE. ZERO) THEN
         WRITE(LUPRI,'(//,5X,A,A,A,/,5X,A,1P,D22.15,/)')
     &   'Error in ',SECNAM,':',
     &   'Cholesky weight in memory split, WCHOL = ',WCHOL
         CALL QUIT('Error in '//SECNAM)
      ENDIF

C     Calculate max. size of integral/amplitude matrix.
C     =================================================

      XDIM = ZERO
      DO ISYMBJ = 1,NSYM
         ISYMAI = MULD2H(ISYMBJ,ISYINT)
         XDIM   = XDIM + NX1AMA(ISYMAI)*XT1AM(ISYMBJ)
      ENDDO
      IF ((XDIM.LE.ZERO) .OR. (XDIM.GT.XT2SQ(ISYINT))) THEN
         WRITE(LUPRI,'(//,5X,A,A,1P,D22.15,/)')
     &   SECNAM,': ERROR: XDIM = ',XDIM
         CALL QUIT('Error in '//SECNAM)
      ENDIF

C     Check if we can store the entire amplitude array while still
C     holding 1/LFRAC of the Cholesky vectors in each symmetry.
C     ============================================================

      XLFT = XWORK - XDIM
      IF (XLFT .GT. ZERO) THEN

         DIFF  = DABS(XT2SQ(ISYINT) - XDIM)
         LALLA = DIFF .LT. TINY
         
         LWRK = INT(XDIM)
         LEFT = LWORK - LWRK

         IF (LIDEN) THEN
         
C           L1 = L2 (CC_CIA1 or CC_CIA2).
C           -----------------------------

            IF (LALLA) THEN
               DO ISYCHO = 1,NSYM
                  IF (NTOVEC(ISYCHO) .GT. LFRAC) THEN
                     MINVEC = NTOVEC(ISYCHO)/LFRAC
                  ELSE
                     MINVEC = MIN(NTOVEC(ISYCHO),MINDEF)
                  ENDIF
                  ISYMAI = MULD2H(ISYCHO,ISYCH1)
                  NEED   = NT1AM(ISYMAI)*MINVEC
                  IF (NEED .GT. LEFT) GO TO 100
               ENDDO
            ELSE
               DO ISYCHO = 1,NSYM
                  IF (NTOVEC(ISYCHO) .GT. LFRAC) THEN
                     MINVEC = NTOVEC(ISYCHO)/LFRAC
                  ELSE
                     MINVEC = MIN(NTOVEC(ISYCHO),MINDEF)
                  ENDIF
                  ISYMAI = MULD2H(ISYCHO,ISYCH1)
                  NEED = NT1AM(ISYMAI)*(2*MINVEC + 1)
                  IF (NEED .GT. LEFT) GO TO 100
               ENDDO
            ENDIF

         ELSE

C           L1 != L2 (CC_CIA3 or CC_CIA4_1 or CC_CIA4_3).
C           ---------------------------------------------

            IF (LALLA) THEN
               DO ISYCHO = 1,NSYM
                  IF (NTOVEC(ISYCHO) .GT. LFRAC) THEN
                     MINVEC = NTOVEC(ISYCHO)/LFRAC
                  ELSE
                     MINVEC = MIN(NTOVEC(ISYCHO),MINDEF)
                  ENDIF
                  ISYMAI = MULD2H(ISYCHO,ISYCH1)
                  ISYMBJ = MULD2H(ISYCHO,ISYCH2)
                  NEED   = (NT1AM(ISYMAI) + NT1AM(ISYMBJ))*MINVEC
                  IF (NEED .GT. LEFT) GO TO 100
               ENDDO
            ELSE
               IF (SYMTRZ .AND. CIAMIO) THEN
                  DO ISYCHO = 1,NSYM
                     IF (NTOVEC(ISYCHO) .GT. LFRAC) THEN
                        MINVEC = NTOVEC(ISYCHO)/LFRAC
                     ELSE
                        MINVEC = MIN(NTOVEC(ISYCHO),MINDEF)
                     ENDIF
                     ISYMAI = MULD2H(ISYCHO,ISYCH1)
                     ISYMBJ = MULD2H(ISYCHO,ISYCH2)
                     NEED   = (NT1AM(ISYMAI) + NT1AM(ISYMBJ))
     &                       *(2*MINVEC + 1)
                     IF (NEED .GT. LEFT) GO TO 100
                  ENDDO
               ELSE
                  DO ISYCHO = 1,NSYM
                     IF (NTOVEC(ISYCHO) .GT. LFRAC) THEN
                        MINVEC = NTOVEC(ISYCHO)/LFRAC
                     ELSE
                        MINVEC = MIN(NTOVEC(ISYCHO),MINDEF)
                     ENDIF
                     ISYMAI = MULD2H(ISYCHO,ISYCH1)
                     ISYMBJ = MULD2H(ISYCHO,ISYCH2)
                     NEED   = (NT1AM(ISYMAI) + NT1AM(ISYMBJ))
     &                       *(MINVEC + 1)
                     IF (NEED .GT. LEFT) GO TO 100
                  ENDDO
               ENDIF
            ENDIF

         ENDIF

C        If we reach this point, "full" calculation is possible.
C        -------------------------------------------------------

         RETURN

      ENDIF

C     Ressume processing here if batching is needed.
C     ----------------------------------------------

  100 CONTINUE

C     Calculate memory split IDIV according to relative
C     sizes of integrals and Cholesky vectors.
C     =================================================

      MAXAIJ = LAIJ(1)
      DO ISYM = 2,NSYM
         MAXAIJ = MAX(MAXAIJ,LAIJ(ISYM))
      ENDDO

      IF (LIDEN) THEN

C        L1 = L2 (CC_CIA2).
C        Integrals/amplitudes: XDIM
C        MO Cholesky vectors : 2*L1(ai,J) + one additional vector.
C        ---------------------------------------------------------

         XCHO = ZERO
         DO ISYCHO = 1,NSYM
            ISYMAI = MULD2H(ISYCHO,ISYCH1)
            XTST   = XT1AM(ISYMAI)*(2*NTOVEC(ISYCHO) + 1)
            IF (XCHO .LT. XTST) THEN
               XCHO = XTST
            ENDIF
         ENDDO

         IF (XCHO .LE. ZERO) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Error in ',SECNAM,':'
            WRITE(LUPRI,'(5X,A,1P,D30.15)')
     &      'Calculated max. dimension of Cholesky part: ',XCHO
            WRITE(LUPRI,'(5X,A)')
     &      'Number of Cholesky vectors supplied:'
            WRITE(LUPRI,'(8I10)')
     &      (NTOVEC(ISYCHO),ISYCHO=1,NSYM)
            CALL QUIT('Error in '//SECNAM)
         ENDIF

         XCHOF = WCHOL*XCHO
         IF (XDIM .GT. XCHOF) THEN
            ISPLIT = 1
            XFRAC  = XDIM/XCHOF
         ELSE
            ISPLIT = 2
            XFRAC  = XCHOF/XDIM
         ENDIF

         IFRAC = INT(XFRAC)
         IDIV  = MAX(IFRAC+1,2)

C        Split memory.
C        -------------

         LWRK = LWORK/IDIV
         IF (ISPLIT .EQ. 1) THEN
            LWRK = LWORK - LWRK
         ENDIF

C        Make sure that we do not reserve too much space
C        for Cholesky vectors.
C        -----------------------------------------------

         LEFT = LWORK - LWRK
         XLFT = ONE*LEFT
         IF (XLFT .GT. XCHO) THEN
            LEFT = INT(XCHO)
            LWRK = LWORK - LEFT
         ENDIF

C        Make sure that integral/amplitude batch can be held in
C        core for at least 1 virtual. This section should be
C        relevant only for work spaces in the neighborhood of
C        the absolute minimal required.
C        ------------------------------------------------------

         IF (LWRK .LT. MAXAIJ) THEN
            LWRK = MIN(LWORK,MAXAIJ)
         ENDIF

         LEFT = LWORK - LWRK
         DO ISYMB = 1,NSYM
            DO ISYCHO = 1,NSYM
               ISYMAI = MULD2H(ISYCHO,ISYCH1)
               ISYMBJ = ISYMAI
               ISYMJ  = MULD2H(ISYMBJ,ISYMB)
               NTST   = NT1AM(ISYMAI) + NX1AMA(ISYMAI) + NRHF(ISYMJ)
               IF (LEFT .LT. NTST) THEN
                  LEFT = MIN(LWORK,NTST)
                  LWRK = LWORK - LEFT
               ENDIF
            ENDDO
         ENDDO

      ELSE

         IF (SYMTRZ .AND. CIAMIO) THEN

C           L1 != L2, symmetrization with minmum I/O (CC_CIA4_1).
C           Integrals/amplitudes: XDIM.
C           MO Cholesky vectors : 2*(L1(ai,J) + L2(bj,J)) + 2 additional.
C           -------------------------------------------------------------

            XCHO = ZERO
            DO ISYCHO = 1,NSYM
               ISYMAI = MULD2H(ISYCHO,ISYCH1)
               ISYMBJ = MULD2H(ISYCHO,ISYCH2)
               XT1TOT = XT1AM(ISYMAI) + XT1AM(ISYMBJ)
               XTST   = XT1TOT*(2*NTOVEC(ISYCHO) + 1)
               IF (XCHO .LT. XTST) XCHO = XTST
            ENDDO

            IF (XCHO .LE. ZERO) THEN
               WRITE(LUPRI,'(//,5X,A,A,A)')
     &         'Error in ',SECNAM,':'
               WRITE(LUPRI,'(5X,A,1P,D30.15)')
     &         'Calculated max. dimension of Cholesky part: ',XCHO
               WRITE(LUPRI,'(5X,A)')
     &         'Number of Cholesky vectors supplied:'
               WRITE(LUPRI,'(8I10)')
     &         (NTOVEC(ISYCHO),ISYCHO=1,NSYM)
               CALL QUIT('Error in '//SECNAM)
            ENDIF

            XCHOF = WCHOL*XCHO
            IF (XDIM .GT. XCHOF) THEN
               ISPLIT = 1
               XFRAC  = XDIM/XCHOF
            ELSE
               ISPLIT = 2
               XFRAC  = XCHOF/XDIM
            ENDIF

            IFRAC = INT(XFRAC)
            IDIV  = MAX(IFRAC+1,2)

C           Split memory.
C           -------------

            LWRK = LWORK/IDIV
            IF (ISPLIT .EQ. 1) THEN
               LWRK = LWORK - LWRK
            ENDIF

C           Make sure that we do not reserve too much space
C           for Cholesky vectors.
C           -----------------------------------------------

            LEFT = LWORK - LWRK
            XLFT = ONE*LEFT
            IF (XLFT .GT. XCHO) THEN
               LEFT = INT(XCHO)
               LWRK = LWORK - LEFT
            ENDIF

C           Make sure that integral/amplitude batch can be held in
C           core for at least 1 virtual. This section should be
C           relevant only for work spaces in the neighborhood of
C           the absolute minimal required.
C           ------------------------------------------------------

            IF (LWRK .LT. MAXAIJ) THEN
               LWRK = MIN(LWORK,MAXAIJ)
            ENDIF

            LEFT = LWORK - LWRK
            DO ISYMB = 1,NSYM
               DO ISYCHO = 1,NSYM
                  ISYMAI = MULD2H(ISYCHO,ISYCH1)
                  ISYMBJ = MULD2H(ISYCHO,ISYCH2)
                  ISYMJ1 = MULD2H(ISYMBJ,ISYMB)
                  ISYMJ2 = MULD2H(ISYMAI,ISYMB)
                  NTST   = NT1AM(ISYMAI)  + NT1AM(ISYMBJ)
     &                   + NX1AMA(ISYMAI) + NRHF(ISYMJ1)
     &                   + NX1AMA(ISYMBJ) + NRHF(ISYMJ2)
                  IF (LEFT .LT. NTST) THEN
                     LEFT = MIN(LWORK,NTST)
                     LWRK = LWORK - LEFT
                  ENDIF
               ENDDO
            ENDDO

         ELSE

C           L1 != L2, no symmetrization, or symmetrization by 2 passes
C                     through Cholesky files (CC_CIA4_3).
C           Integrals/amplitudes: XDIM.
C           MO Cholesky vectors : L1(ai,J) + L2(bj,J) + 2 additional.
C           ----------------------------------------------------------

            XCHO = ZERO
            DO ISYCHO = 1,NSYM
               ISYMAI = MULD2H(ISYCHO,ISYCH1)
               ISYMBJ = MULD2H(ISYCHO,ISYCH2)
               XT1TOT = XT1AM(ISYMAI) + XT1AM(ISYMBJ)
               XTST   = XT1TOT*(NTOVEC(ISYCHO) + 1)
               IF (XCHO .LT. XTST) XCHO = XTST
            ENDDO

            IF (XCHO .LE. ZERO) THEN
               WRITE(LUPRI,'(//,5X,A,A,A)')
     &         'Error in ',SECNAM,':'
               WRITE(LUPRI,'(5X,A,1P,D30.15)')
     &         'Calculated max. dimension of Cholesky part: ',XCHO
               WRITE(LUPRI,'(5X,A)')
     &         'Number of Cholesky vectors supplied:'
               WRITE(LUPRI,'(8I10)')
     &         (NTOVEC(ISYCHO),ISYCHO=1,NSYM)
               WRITE(LUPRI,*)
               CALL QUIT('Error in '//SECNAM)
            ENDIF

            XCHOF = WCHOL*XCHO
            IF (XDIM .GT. XCHOF) THEN
               ISPLIT = 1
               XFRAC  = XDIM/XCHOF
            ELSE
               ISPLIT = 2
               XFRAC  = XCHOF/XDIM
            ENDIF

            IFRAC = INT(XFRAC)
            IDIV  = MAX(IFRAC+1,2)

C           Split memory.
C           -------------

            LWRK = LWORK/IDIV
            IF (ISPLIT .EQ. 1) THEN
               LWRK = LWORK - LWRK
            ENDIF

C           Make sure that we do not reserve too much space
C           for Cholesky vectors.
C           -----------------------------------------------

            LEFT = LWORK - LWRK
            XLFT = ONE*LEFT
            IF (XLFT .GT. XCHO) THEN
               LEFT = INT(XCHO)
               LWRK = LWORK - LEFT
            ENDIF

C           Make sure that integral/amplitude batch can be held in
C           core for at least 1 virtual. This section should be
C           relevant only for work spaces in the neighborhood of
C           the absolute minimal required.
C           ------------------------------------------------------

            IF (LWRK .LT. MAXAIJ) THEN
               LWRK = MIN(LWORK,MAXAIJ)
            ENDIF

            LEFT = LWORK - LWRK
            DO ISYMB = 1,NSYM
               DO ISYCHO = 1,NSYM
                  ISYMAI = MULD2H(ISYCHO,ISYCH1)
                  ISYMBJ = MULD2H(ISYCHO,ISYCH2)
                  ISYMJ  = MULD2H(ISYMBJ,ISYMB)
                  NTST   = NT1AM(ISYMAI)  + NT1AM(ISYMBJ)
     &                   + NX1AMA(ISYMAI) + NRHF(ISYMJ)
                  IF (LEFT .LT. NTST) THEN
                     LEFT = MIN(LWORK,NTST)
                     LWRK = LWORK - LEFT
                  ENDIF
               ENDDO
            ENDDO

         ENDIF

      ENDIF

      RETURN
      END
C  /* Deck inv_btch */
      SUBROUTINE INV_BTCH(NUMA,ABEG,ISABEG,AEND,ISAEND,IOFA1,LVIRA)
C
C     Thomas Bondo Pedersen, April 2003.
C
C     Purpose: Get start and end info of a virtual batch from IOFA1 and LVIRA.
C
#include <implicit.h>
      INTEGER ABEG, AEND
      INTEGER IOFA1(8), LVIRA(8)
#include <ccorb.h>
#include <ccsdsym.h>

C     Find first and last a.
C     ----------------------

      ISABEG = 0
      DO ISYMA = 1,NSYM
         IF (LVIRA(ISYMA) .GT. 0) THEN
            ABEG   = IOFA1(ISYMA) + IVIR(ISYMA) - NRHFT
            ISABEG = ISYMA
            GO TO 50
         ENDIF
      ENDDO
   50 IF ((ISABEG.LE.0) .OR. (ISABEG.GT.NSYM)) THEN
         CALL QUIT('Symmetry error [1] in INV_BTCH')
      ENDIF

      ISAEND = 0
      DO ISYMA = NSYM,1,-1
         IF (LVIRA(ISYMA) .GT. 0) THEN
            AEND   = IOFA1(ISYMA) + IVIR(ISYMA) - NRHFT + LVIRA(ISYMA)
     &             - 1
            ISAEND = ISYMA
            GO TO 60
         ENDIF
      ENDDO
   60 IF ((ISAEND.LE.0) .OR. (ISAEND.GT.NSYM)) THEN
         CALL QUIT('Symmetry error [2] in INV_BTCH')
      ENDIF

C     Calculate length of batch.
C     --------------------------

      NUMA = 0
      DO ISYMA = ISABEG,ISAEND
         NUMA = NUMA + LVIRA(ISYMA)
      ENDDO
      IF ((NUMA.LE.0) .OR. (NUMA.GT.NVIRT)) THEN
         CALL QUIT('Batch length error in INV_BTCH')
      ENDIF

      RETURN
      END
C  /* Deck cc_chor1d1 */
      SUBROUTINE CC_CHOR1D1(DENJK,DENJB,DENCB,XL1AM,XR2AM,WORK,LWORK,
     &                      ISYMXL,ISYMXR,LAIJ,FILMUL,MEM)
C
C     Thomas Bondo Pedersen, April 2003.
C
C     Purpose: Calculate right first-order CC2 density matrix blocks
C
C        D(jk)  = D(jk)  - sum(#ai#b) R(#ai,#bj) * L(#ai,#bk)
C
C        D(c#b) = D(c#b) + sum(#aij)  L(#ai,cj)  * R(#ai,#bj)
C
C        D(j#b) = D(j#b) + sum(#ai) [2*R(#ai,#bj)-R(#aj,#bi)] * L(#ai)
C
C     for batches of virtual orbitals. It is assumed that the density
C     matrices are passed in full storage (i.e. for all orbital indices).
C     Furthermore, the left doubles vector must be available on disk for
C     the current batch of a's [L(#ai,dl), stored on direct-access files
C     FILMUL(ISYMDL) with the dl-index labelling records].
C
C     Information about the a- and b-batches is taken from ciarc.h, except
C     array LAIJ(*) which should contain the length of the array #ai,j for
C     a given symmetry ISYAIJ, as needed for setting up a batched read of
C     the left doubles array from disk.
C
C     NB!!! On exit, R is replaced by 2CME amplitudes:
C           R(#ai,#bj) --> 2*R(#ai,#bj) - R(#aj,#bi)
C
C
C     TODO/FIXME: For small batch sizes in calculations with symmetry,
C                 too much of the left doubles array will be read. May in
C                 principle be fixed by reading the relevant L blocks when
C                 needed. This, however, leads to two reads for larger
C                 batch sizes (and, notably, in calculations without symmetry).
C
#include <implicit.h>
      DIMENSION DENJK(*), DENJB(*), DENCB(*), XL1AM(*), XR2AM(*)
      DIMENSION WORK(LWORK)
      INTEGER   LAIJ(8)
      CHARACTER*8 FILMUL(8)
#include <maxorb.h>
#include <ccisvi.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>
#include <ciarc.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CC_CHOR1D1')

      LOGICAL LOCDBG, OLDDX
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER DBEG, DEND
      INTEGER BBEG, BEND
      INTEGER BPBEG, BPEND

      INTEGER IOFD1(8), LVIRD(8), IX1AMD(8,8), NX1AMD(8)
      INTEGER IL2SQ(8,8), NL2SQ
      INTEGER LVIRBP(8), IOFBPB(8), IOFBPD(8)

      PARAMETER (INFO = 20)
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0)

C     Start timing.
C     -------------

      TIMT = SECOND()
      TIMR = ZERO
      TIMO = ZERO
      TIMV = ZERO

C     Set result symmetry.
C     --------------------

      ISYDEN = MULD2H(ISYMXL,ISYMXR)
      IF ((ISYDEN.LE.0) .OR. (ISYDEN.GT.NSYM)) THEN
         WRITE(LUPRI,'(//,5X,A,A,I10)')
     &   SECNAM,': ERROR: ISYDEN =',ISYDEN
         CALL QUIT('Error in '//SECNAM)
      ENDIF

C     Debug print.
C     ------------

      IF (LOCDBG) THEN
         XL1NRM = DDOT(NT1AM(ISYMXL),XL1AM,1,XL1AM,1)
         XR2NRM = ZERO
         CALL CC_CYINRM(XR2AM,ISYMXR,XR2NRM)
         XJKNRM = DDOT(NMATIJ(ISYDEN),DENJK,1,DENJK,1)
         XJBNRM = DDOT(NT1AM(ISYDEN),DENJB,1,DENJB,1)
         XCBNRM = DDOT(NMATAB(ISYDEN),DENCB,1,DENCB,1)
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) SECNAM,': Entering: Debug info:'
         WRITE(LUPRI,*) 'Symmetry of left  amplitudes: ',ISYMXL
         WRITE(LUPRI,*) 'Symmetry of right amplitudes: ',ISYMXR
         WRITE(LUPRI,*) 'Symmetry of density matrix  : ',ISYDEN
         WRITE(LUPRI,*) ' Sq-norm of left singles amplitudes: ',XL1NRM
         WRITE(LUPRI,*) ' 2--norm of left singles amplitudes: ',
     &                  DSQRT(XL1NRM)
         WRITE(LUPRI,*) ' Sq-norm of right doubl. amplitudes: ',XR2NRM,
     &                  ' (packed)'
         WRITE(LUPRI,*) ' 2--norm of right doubl. amplitudes: ',
     &                  DSQRT(XR2NRM),' (packed)'
         WRITE(LUPRI,*) ' Sq-norm of DENJK                  : ',XJKNRM
         WRITE(LUPRI,*) ' 2--norm of DENJK                  : ',
     &                  DSQRT(XJKNRM)
         WRITE(LUPRI,*) ' Sq-norm of DENJB                  : ',XJBNRM
         WRITE(LUPRI,*) ' 2--norm of DENJB                  : ',
     &                  DSQRT(XJBNRM)
         WRITE(LUPRI,*) ' Sq-norm of DENCB                  : ',XCBNRM
         WRITE(LUPRI,*) ' 2--norm of DENCB                  : ',
     &                  DSQRT(XCBNRM)
         CALL FLSHFO(LUPRI)
      ENDIF

C     Allocate scratch space for largest D(c,#b).
C     -------------------------------------------

      MAXCB = 0
      DO ISYMB = 1,NSYM
         ISYMC = MULD2H(ISYMB,ISYDEN)
         NEED  = NVIR(ISYMC)*LVIRB(ISYMB)
         MAXCB = MAX(MAXCB,NEED)
      ENDDO

      KDNCB = 1
      KEND1 = KDNCB + MAXCB
      LWRK1 = LWORK - KEND1 + 1

      IF (LWRK1 .LE. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   'Insufficient memory in ',SECNAM
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     &   'Need (more than): ',KEND1-1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Initialize memory info.
C     -----------------------

      MXMEMT = KEND1 - 1

C     Start batch loop over d in L(#ai,dl).
C     -------------------------------------

      IBATD = 0
      DBEG  = 1
      NUMD  = 0
  100 CONTINUE

         DBEG = DBEG + NUMD
         IF (DBEG .GT. NVIRT) GO TO 200

C        Time this batch.
C        ----------------

         TDBAT = SECOND()

C        Update batch counter.
C        ---------------------

         IBATD = IBATD + 1

C        Find out how many d's can be treated.
C        -------------------------------------

         CALL CC_CYIBTA(LAIJ,NX1AMA,LVIRD,DBEG,NUMD,MEM,LWRK1,ISYMXL)

C        Check for errors.
C        -----------------

         DEND = DBEG + NUMD - 1
         IF (DEND .GT. NVIRT) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Batch error in ',SECNAM,' !'
            WRITE(LUPRI,'(5X,A)')
     &      '- dumping presumably most relevant info:'
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,I10)')
     &      'First D           : ',DBEG,
     &      'Last  D           : ',DEND,
     &      'Number of virtuals: ',NVIRT
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     &      'Available memory  : ',LWRK1,
     &      'Needed    memory  : ',MEM
            CALL QUIT('Batch error in '//SECNAM)
         ENDIF

         IF (NUMD .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory for d-batch in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/)')
     &      'Available memory: ',LWRK1
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Initialize local timings.
C        -------------------------

         DTMR = ZERO
         DTMV = ZERO
         DTMO = ZERO

C        Get symmetries of first and last d in batch.
C        --------------------------------------------

         ISDBEG = ISVI(DBEG)
         ISDEND = ISVI(DEND)

C        Set up index arrays.
C        --------------------

         CALL GET_INDVIR(DBEG,DEND,LVIRD,IOFD1,IX1AMD,NX1AMD)

         CALL IZERO(IL2SQ,64)
         ICOUNT = 0
         DO ISYMDL = 1,NSYM
            ISYMAI = MULD2H(ISYMDL,ISYMXL)
            IL2SQ(ISYMAI,ISYMDL) = ICOUNT
            ICOUNT = ICOUNT + NX1AMA(ISYMAI)*NX1AMD(ISYMDL)
         ENDDO
         NL2SQ = ICOUNT

         IF (NL2SQ .NE. MEM) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Error detected in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,/)')
     &      'NL2SQ is calculated to be ',NL2SQ,
     &      'MEM determines alloc. as  ',MEM,
     &      'These numbers must be equal....'
            CALL QUIT('Error in '//SECNAM)
         ENDIF

C        Allocation.
C        -----------

         KL2AM = KEND1
         KEND2 = KL2AM + NL2SQ
         LWRK2 = LWORK - KEND2 + 1

         KLAST  = KEND2 - 1
         MXMEMT = MAX(MXMEMT,KLAST)

C        Find the union of d- and b-batches.
C        -----------------------------------

         CALL INV_BTCH(NUMB,BBEG,ISBBEG,BEND,ISBEND,IOFB1,LVIRB)

         BPBEG = MAX(BBEG,DBEG)
         BPEND = MIN(BEND,DEND)
         NUMBP = BPEND - BPBEG + 1
         NUMBP = MAX(NUMBP,0)

         IF (NUMBP .GT. 0) THEN
            CALL IZERO(LVIRBP,NSYM)
            CALL IZERO(IOFBPB,NSYM)
            CALL IZERO(IOFBPD,NSYM)
            DO IBP = BPBEG,BPEND
               ISYMBP = ISVI(IBP)
               LVIRBP(ISYMBP) = LVIRBP(ISYMBP) + 1
            ENDDO
            IBP1 = BPBEG
            DO ISYMBP = ISVI(BPBEG),ISVI(BPEND)
               IOFBP1 = IBP1 + NRHFT - IVIR(ISYMBP)
               IBPB   = ABS(IOFBP1-IOFB1(ISYMBP))
               IBPD   = ABS(IOFBP1-IOFD1(ISYMBP))
               IOFBPB(ISYMBP) = IBPB + 1
               IOFBPD(ISYMBP) = IBPD + 1
               IBP1   = IBP1 + LVIRBP(ISYMBP)
            ENDDO
         ENDIF

         IF (LOCDBG) THEN
            WRITE(LUPRI,*)
            ID1 = IOFD1(ISDBEG)
            ID2 = IOFD1(ISDEND) + LVIRD(ISDEND) - 1
            WRITE(LUPRI,*) SECNAM,': d-batch number ',IBATD
            WRITE(LUPRI,*) 'First d: ',DBEG,'(symmetry ',ISVI(DBEG),
     &                     ' number ',ID1,')'
            WRITE(LUPRI,*) 'Last  d: ',DEND,'(symmetry ',ISVI(DEND),
     &                     ' number ',ID2,')'
            IB1 = IOFB1(ISBBEG)
            IB2 = IOFB1(ISBEND) + LVIRB(ISBEND) - 1
            WRITE(LUPRI,*) 'First b: ',BBEG,'(symmetry ',ISBBEG,
     &                     ' number ',IB1,')'
            WRITE(LUPRI,*) 'Last  b: ',BEND,'(symmetry ',ISBEND,
     &                     ' number ',IB2,')'
            WRITE(LUPRI,*) 'Dimension of union: ',NUMBP
            IF (NUMBP .GT. 0) THEN
               WRITE(LUPRI,*) 'LVIRBP:'
               WRITE(LUPRI,'(8I10)') (LVIRBP(I),I=1,NSYM)
               WRITE(LUPRI,*) 'IOFBPB:'
               WRITE(LUPRI,'(8I10)') (IOFBPB(I),I=1,NSYM)
               WRITE(LUPRI,*) 'IOFBPD:'
               WRITE(LUPRI,'(8I10)') (IOFBPD(I),I=1,NSYM)
            ENDIF
            WRITE(LUPRI,*)
            CALL FLSHFO(LUPRI)
         ENDIF

C        Read left amplitudes: L(#ai,#dl).
C        ---------------------------------

         DTMR = SECOND()
         DO ISYMDL = 1,NSYM

            ISYMAI = MULD2H(ISYMDL,ISYMXL)
            LENAI  = NX1AMA(ISYMAI)

            IF (LENAI .GT. 0) THEN

               LRLEN = 2*NX1AMA(ISYMAI)
               LUNIT = -1
               CALL GPOPEN(LUNIT,FILMUL(ISYMDL),'OLD','DIRECT',
     &                     'UNFORMATTED',LRLEN,OLDDX)

               DO ISYML = 1,NSYM

                  ISYMD = MULD2H(ISYML,ISYMDL)
                  LEND  = LVIRD(ISYMD)

                  DO L = 1,NRHF(ISYML)

                     NDL1 = IT1AM(ISYMD,ISYML) + NVIR(ISYMD)*(L - 1)
     &                    + IOFD1(ISYMD)
                     LDL1 = IX1AMD(ISYMD,ISYML) + LEND*(L - 1)
                     KOFF = KL2AM + IL2SQ(ISYMAI,ISYMDL) + LENAI*LDL1

                     CALL CHO_MOREAD(WORK(KOFF),LENAI,LEND,NDL1,LUNIT)

                  ENDDO

               ENDDO

               CALL GPCLOSE(LUNIT,'KEEP')
 
            ENDIF

         ENDDO
         DTMR = SECOND() - DTMR
         TIMR = TIMR     + DTMR

C        Calculate occupied and virtual density blocks.
C        ==============================================

         DO ISYMDL = 1,NSYM

            ISYMAI = MULD2H(ISYMDL,ISYMXL)
            ISYMBJ = MULD2H(ISYMAI,ISYMXR)
            ISYMCJ = ISYMDL
            ISYMBK = ISYMDL

            DO ISYML = 1,NSYM

               ISYMD = MULD2H(ISYML,ISYMDL)

C              Virtual part (using scratch space):
C              D(#c,#b) = D(#c,#b)
C                       + sum(j) sum(#ai) L(#ai,#c;j) * R(#ai,#b;j).
C              -----------------------------------------------------

               DTIME = SECOND()

               ISYMC = ISYMD
               ISYMJ = ISYML
               ISYMB = MULD2H(ISYMJ,ISYMBJ)
               LENC  = LVIRD(ISYMC)
               LENB  = LVIRB(ISYMB)
               LENCB = LENC*LENB
               IF (LENCB .GT. 0) THEN

                  CALL DZERO(WORK(KDNCB),LENCB)

                  DO J = 1,NRHF(ISYMJ)

                     LCJ1 = IX1AMD(ISYMC,ISYMJ) + LENC*(J - 1)
                     LBJ1 = IX1AMB(ISYMB,ISYMJ) + LENB*(J - 1)

                     KOFFL = KL2AM + IL2SQ(ISYMAI,ISYMCJ)
     &                     + NX1AMA(ISYMAI)*LCJ1
                     KOFFR = IX2SQ(ISYMAI,ISYMBJ)
     &                     + NX1AMA(ISYMAI)*LBJ1 + 1

                     NTAI = MAX(NX1AMA(ISYMAI),1)

                     CALL DGEMM('T','N',LENC,LENB,NX1AMA(ISYMAI),
     &                          ONE,WORK(KOFFL),NTAI,XR2AM(KOFFR),NTAI,
     &                          ONE,WORK(KDNCB),LENC)

                  ENDDO

                  DO LB = 1,LENB
                     B     = IOFB1(ISYMB) + LB - 1
                     KOFF1 = KDNCB + LENC*(LB - 1)
                     KOFF2 = IMATAB(ISYMC,ISYMB)
     &                     + NVIR(ISYMC)*(B - 1) + IOFD1(ISYMC)
                     CALL DAXPY(LENC,ONE,WORK(KOFF1),1,DENCB(KOFF2),1)
                  ENDDO

               ENDIF

               DTIME = SECOND() - DTIME
               DTMV  = DTMV     + DTIME
               TIMV  = TIMV     + DTIME

C              Occupied part (only for #b' [union of #b and #d]):
C              D(jk) = D(jk) - sum(#ai,#b) R(#ai,#b';j) * L(#ai,#b';k).
C              --------------------------------------------------------

               IF (NUMBP .GT. 0) THEN

                  ISYMB = ISYMD
                  ISYMK = ISYML
                  ISYMJ = MULD2H(ISYMB,ISYMBJ)

                  IF (LVIRBP(ISYMB) .GT. 0) THEN

                     DTIME = SECOND()

                     LENAIB = NX1AMA(ISYMAI)*LVIRBP(ISYMB)

                     DO K = 1,NRHF(ISYMK)

                        LBK1  = IX1AMD(ISYMB,ISYMK)
     &                        + LVIRD(ISYMB)*(K - 1) + IOFBPD(ISYMB)
                        KOFFL = KL2AM + IL2SQ(ISYMAI,ISYMBK)
     &                        + NX1AMA(ISYMAI)*(LBK1 - 1)

                        DO J = 1,NRHF(ISYMJ)

                           LBJ1 = IX1AMB(ISYMB,ISYMJ)
     &                          + LVIRB(ISYMB)*(J - 1) + IOFBPB(ISYMB)
                           JK   = IMATIJ(ISYMJ,ISYMK)
     &                          + NRHF(ISYMJ)*(K - 1) + J

                           KOFFR = IX2SQ(ISYMAI,ISYMBJ)
     &                           + NX1AMA(ISYMAI)*(LBJ1 - 1) + 1

                           DENJK(JK) = DENJK(JK)
     &                               - DDOT(LENAIB,XR2AM(KOFFR),1,
     &                                             WORK(KOFFL),1)

                        ENDDO

                     ENDDO

                     DTIME = SECOND() - DTIME
                     DTMO  = DTMO     + DTIME
                     TIMO  = TIMO     + DTIME

                  ENDIF

               ENDIF

            ENDDO

         ENDDO

C        Print timings for this d-batch.
C        -------------------------------

         IF ((IPRINT.GE.INFO) .OR. LOCDBG) THEN
            TDBAT = SECOND() - TDBAT
            WRITE(LUPRI,'(/,5X,A,A,I6,A)')
     &      SECNAM,': d-batch number',IBATD,':'
            WRITE(LUPRI,'(5X,A,I10,A,I10)')
     &      'Memory used: ',KLAST,'. Available: ',LWORK
            WRITE(LUPRI,'(5X,A,I6)')
     &      '#d:',NUMD
            WRITE(LUPRI,'(5X,A,I6,A,I2)')
     &      'First d: ',IOFD1(ISDBEG),' symmetry:',ISDBEG
            WRITE(LUPRI,'(5X,A,I6,A,I2)')
     &      'Last  d: ',IOFD1(ISDEND)+LVIRD(ISDEND)-1,
     &      ' symmetry:',ISDEND
            WRITE(LUPRI,'(5X,A,F10.2,A)')
     &      'Time used reading left amplitudes     : ',DTMR,' seconds.'
            WRITE(LUPRI,'(5X,A,F10.2,A)')
     &      'Time used calculating occupied density: ',DTMO,' seconds.'
            WRITE(LUPRI,'(5X,A,F10.2,A)')
     &      'Time used calculating virtual  density: ',DTMV,' seconds.'
            WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &      'Total time for this d-batch           : ',TDBAT,' seconds.'
            CALL FLSHFO(LUPRI)
         ENDIF

C        Go to next d-batch.
C        -------------------

         GO TO 100

C     Exit point for d-batch loop.
C     ----------------------------

  200 CONTINUE

C     Calculate de-excitation density block.
C     NB!!! XR2AM is replaced with 2CME here!
C     =======================================

      TIMD = SECOND()
      CALL CC_CYITCME(XR2AM,WORK,LWORK,ISYMXR,KLST1)
      CALL CC_CYICNI(DENJB,XL1AM,XR2AM,WORK,LWORK,ISYMXL,ISYMXR,KLST2)
      TIMD = SECOND() - TIMD

      MEM = MAX(MXMEMT,KLST1,KLST2)

C     Debug: print.
C     -------------

      IF (LOCDBG) THEN
         XL1NRM = DDOT(NT1AM(ISYMXL),XL1AM,1,XL1AM,1)
         XR2NRM = ZERO
         CALL CC_CYINRM(XR2AM,ISYMXR,XR2NRM)
         XJKNRM = DDOT(NMATIJ(ISYDEN),DENJK,1,DENJK,1)
         XJBNRM = DDOT(NT1AM(ISYDEN),DENJB,1,DENJB,1)
         XCBNRM = DDOT(NMATAB(ISYDEN),DENCB,1,DENCB,1)
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) SECNAM,': Exiting: Debug info:'
         WRITE(LUPRI,*) 'Symmetry of left  amplitudes: ',ISYMXL
         WRITE(LUPRI,*) 'Symmetry of right amplitudes: ',ISYMXR
         WRITE(LUPRI,*) 'Symmetry of density matrix  : ',ISYDEN
         WRITE(LUPRI,*) ' Sq-norm of left singles amplitudes: ',XL1NRM
         WRITE(LUPRI,*) ' 2--norm of left singles amplitudes: ',
     &                  DSQRT(XL1NRM)
         WRITE(LUPRI,*) ' Sq-norm of right 2CME   amplitudes: ',XR2NRM,
     &                  ' (packed)'
         WRITE(LUPRI,*) ' 2--norm of right 2CME   amplitudes: ',
     &                  DSQRT(XR2NRM),' (packed)'
         WRITE(LUPRI,*) ' Sq-norm of DENJK                  : ',XJKNRM
         WRITE(LUPRI,*) ' 2--norm of DENJK                  : ',
     &                  DSQRT(XJKNRM)
         WRITE(LUPRI,*) ' Sq-norm of DENJB                  : ',XJBNRM
         WRITE(LUPRI,*) ' 2--norm of DENJB                  : ',
     &                  DSQRT(XJBNRM)
         WRITE(LUPRI,*) ' Sq-norm of DENCB                  : ',XCBNRM
         WRITE(LUPRI,*) ' 2--norm of DENCB                  : ',
     &                  DSQRT(XCBNRM)
         CALL FLSHFO(LUPRI)
      ENDIF

C     Print.
C     ------

      IF ((IPRINT.GE.INFO) .OR. LOCDBG) THEN
         TIMT = SECOND() - TIMT
         CALL HEADER('Info from '//SECNAM,-1)
         WRITE(LUPRI,'(5X,A,I6,A,/)')
     &   'Calculation done in',IBATD,' d-batches.'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Memory available: ',LWORK,
     &   'Memory used     : ',MEM
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used reading left amplitudes     : ',TIMR,' seconds.'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used calculating occupied density: ',TIMO,' seconds.'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used calculating virtual  density: ',TIMV,' seconds.'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used calculating de-exc.  density: ',TIMD,' seconds.'
         WRITE(LUPRI,'(5X,A)')
     &   '-----------------------------------------------------------'
         WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &   'Time used in total                    : ',TIMT,' seconds.'
         CALL FLSHFO(LUPRI)
      ENDIF

      RETURN
      END
C  /* Deck cc_chor1dbgao */
      SUBROUTINE CC_CHOR1DBGAO(IPDEN,NPDEN,LISTR,LISTD,MODEL,WORK,LWORK,
     &                         STPDBG)
C
C     Thomas Bondo Pedersen, April 2003.
C
C     Purpose: Calculate Eta(A)*T(B) doubles contribution to the linear
C              response function using the right first-order density.
C              For debugging.
C
#include <implicit.h>
      INTEGER       IPDEN(3,NPDEN)
      CHARACTER*(*) LISTR,LISTD
      CHARACTER*10  MODEL
      DIMENSION     WORK(LWORK)
      LOGICAL       STPDBG
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>
#include <ccr1rsp.h>

      CHARACTER*13 SECNAM
      PARAMETER (SECNAM = 'CC_CHOR1DBGAO')

      LOGICAL LONLYAO
      DATA LONLYAO /.TRUE./

      CHARACTER*8 LABA,LABB

      PARAMETER (THRDFF = 1.0D-14)

      INTEGER ILSTSYM

C     Skip if no densities.
C     ---------------------

      IF (NPDEN .LE. 0) RETURN

C     Print header.
C     -------------

      CALL AROUND('Output from '//SECNAM)
      WRITE(LUPRI,'(5X,A,I10,1X,A,1X,A,A,/)')
     & 'Testing ',NPDEN,MODEL,LISTD(1:3),' AO densities.'
      WRITE(LUPRI,'(5X,A,A,/,5X,A,A)')
     & 'LABELA   Sym.A  LABELB   Sym.B    Frequency',
     & '           Contribution',
     & '-------------------------------------------',
     & '-----------------------'
      CALL FLSHFO(LUPRI)

C     Start loop over densities on disk.
C     ----------------------------------

      DO IDEN = 1,NPDEN

         ILSTDN = IPDEN(3,IDEN)
         ISYDEN = ILSTSYM(LISTD,ILSTDN)
         LABB   = LRTLBL(ILSTDN)
         FREQ   = FRQLRT(ILSTDN)

         ISYMB = ISYDEN

C        Allocation.
C        -----------

         KDENB = 1
         KEND1 = KDENB + N2BST(ISYDEN)
         LWRK1 = LWORK - KEND1 + 1

         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient memory in '//SECNAM//' [1]')
         ENDIF

C        Read density.
C        -------------

         IOPT = 1
         CALL CC_RDRSPD(LISTD,ILSTDN,ISYDEN,IOPT,MODEL,.FALSE.,
     &                  WORK(KDENB))

C        Loop over operators A.
C        ----------------------

         DO IOPA = 1,NLRTLBL

            ISYMA = ISYLRT(IOPA)
            FREQA = FRQLRT(IOPA)

            DIFF = DABS(FREQA + FREQ)

            IF ((ISYMA.EQ.ISYMB) .AND. (DIFF.LE.THRDFF)) THEN

               LABA  = LRTLBL(IOPA)

               KINTA = KEND1
               KEND2 = KINTA + N2BST(ISYMA)
               LWRK2 = LWORK - KEND2 + 1

               IF (LWRK2 .LT. 0) THEN
                  CALL QUIT('Insufficient memory in '//SECNAM//' [2]')
               ENDIF

C              Read AO integrals.
C              ------------------

               IERR = -1
               CALL CCPRPAO(LABA,LONLYAO,WORK(KINTA),IRREP,IREAL,IERR,
     &                      WORK(KEND2),LWRK2)

C              Check.
C              ------

               IF (IERR .GT. 0) THEN
                  WRITE(LUPRI,'(//,5X,A,A)')
     &            SECNAM,': I/O error ocurred in CCPRPAO'
                  WRITE(LUPRI,'(5X,A,I10,A,/,5X,A,2X,A8,/,5X,A,I10)')
     &            'Operator number  : ',IOPA,' (in LRTLBL list)',
     &            'Operator label   : ',LABA,
     &            'Operator symmetry: ',ISYMA
                  WRITE(LUPRI,'(5X,A,I10,/)')
     &            'CCPRPAO returned IERR = ',IERR
                  CALL QUIT('I/O error in CCPRPAO ('//SECNAM//')')
               ELSE IF (IERR .LT. 0) THEN
                  WRITE(LUPRI,'(//,5X,A,A)')
     &          SECNAM,': CCPRPAO failed to determine operator symmetry'
                  WRITE(LUPRI,'(5X,A,I10,A,/,5X,A,2X,A8,/,5X,A,I10)')
     &            'Operator number  : ',IOPA,' (in LRTLBL list)',
     &            'Operator label   : ',LABA,
     &            'Operator symmetry: ',ISYMA
                  WRITE(LUPRI,'(5X,A,I10,/)')
     &            'CCPRPAO returned IERR = ',IERR
                  CALL QUIT('Irrep. error in CCPRPAO ('//SECNAM//')')
               ENDIF

               IF (IRREP .NE. ISYDEN) THEN
                  WRITE(LUPRI,'(//,5X,A,A)')
     &            'Symmetry mismatch in ',SECNAM
                  WRITE(LUPRI,'(5X,A,I10,A,/,5X,A,2X,A8)')
     &            'Operator number  : ',IOPA,' (in LRTLBL list)',
     &            'Operator label   : ',LABA
                  WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &            'CCPRPAO returned: ',IRREP,
     &            '- expected      : ',ISYMA
                  CALL QUIT('Symmetry mismatch in '//SECNAM)
               ENDIF

C              Calculate contribution and print.
C              ---------------------------------

               CON = DDOT(N2BST(ISYMA),WORK(KDENB),1,WORK(KINTA),1)
               WRITE(LUPRI,1) LABA,ISYMA,LABB,ISYMB,FREQ,CON
    1          FORMAT(5X,A8,3X,I2,3X,A8,3X,I2,2X,F12.5,1X,1P,D22.15)
               CALL FLSHFO(LUPRI)

            ENDIF

         ENDDO

      ENDDO

C     Print end of table.
C     -------------------

      WRITE(LUPRI,'(5X,A,A,/)')
     & '--------------------------------------------',
     & '----------------------'
      CALL FLSHFO(LUPRI)

C     Stop if requested.
C     ------------------

      IF (STPDBG) THEN
         WRITE(LUPRI,*) SECNAM,': Stopping as requested!'
         CALL QUIT('STOP in CC_CHOR1DBGAO as requested')
      ENDIF

      RETURN
      END
C  /* Deck cc_chor1dbgmo */
      SUBROUTINE CC_CHOR1DBGMO(IPDEN,NPDEN,LISTR,LISTD,MODEL,WORK,LWORK,
     &                         FILDEN,DELDEN,STPDBG)
C
C     Thomas Bondo Pedersen, April 2003.
C
C     Purpose: Calculate Eta(A)*T(B) doubles contribution to the linear
C              response function using the right first-order density.
C              For debugging.
C
#include <implicit.h>
      INTEGER       IPDEN(3,NPDEN)
      CHARACTER*(*) LISTR,LISTD
      CHARACTER*10  MODEL,FILDEN
      DIMENSION     WORK(LWORK)
      LOGICAL       DELDEN,STPDBG
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>
#include <ccr1rsp.h>
#include <dummy.h>

      CHARACTER*13 SECNAM
      PARAMETER (SECNAM = 'CC_CHOR1DBGMO')

      LOGICAL LONLYAO
      DATA LONLYAO /.TRUE./

      CHARACTER*8 LABA,LABB

      PARAMETER (THRDFF = 1.0D-14)

      INTEGER ADR

      INTEGER ILSTSYM

C     Skip if no densities.
C     ---------------------

      IF (NPDEN .LE. 0) RETURN

C     Open FILDEN.
C     ------------

      LUDEN = 97
      CALL WOPEN2(LUDEN,FILDEN,64,0)

C     Allocation.
C     -----------

      KLAMDP = 1
      KLAMDH = KLAMDP + NGLMDT(1)
      KT1AM  = KLAMDH + NGLMDT(1)
      KEND   = KT1AM  + NT1AM(1)
      LWRK   = LWORK  - KEND + 1

      KEND0 = KT1AM

      IF (LWRK .LT. 0) THEN
         CALL QUIT('Insufficient memory in '//SECNAM//' [0]')
      ENDIF

C     Calculate Lambda matrices.
C     --------------------------

      IOPT = 1
      CALL CC_RDRSP('R0 ',1,1,IOPT,MODEL,WORK(KT1AM),DUMMY)
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),
     &            WORK(KEND),LWRK)

C     Print header.
C     -------------

      CALL AROUND('Output from '//SECNAM)
      WRITE(LUPRI,'(5X,A,I10,1X,A,1X,A,A,/)')
     & 'Testing ',NPDEN,MODEL,LISTD(1:3),' MO densities.'
      WRITE(LUPRI,'(5X,A,A,/,5X,A,A)')
     & 'LABELA   Sym.A  LABELB   Sym.B    Frequency',
     & '           Contribution',
     & '-------------------------------------------',
     & '-----------------------'
      CALL FLSHFO(LUPRI)

C     Start loop over densities on disk.
C     ----------------------------------

      ADR = 1
      DO IDEN = 1,NPDEN

         ILSTDN = IPDEN(3,IDEN)
         ISYDEN = ILSTSYM(LISTD,ILSTDN)
         LABB   = LRTLBL(ILSTDN)
         FREQ   = FRQLRT(ILSTDN)

         ISYMB = ISYDEN

C        Allocation.
C        -----------

         LENDEN = NMATIJ(ISYDEN) + NT1AM(ISYDEN) + NMATAB(ISYDEN)

         KDNIJ = KEND0
         KDNIA = KDNIJ + NMATIJ(ISYDEN)
         KDNAB = KDNIA + NT1AM(ISYDEN)
         KEND1 = KDNAB + NMATAB(ISYDEN)
         LWRK1 = LWORK - KEND1 + 1

         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient memory in '//SECNAM//' [1]')
         ENDIF

         CALL GETWA2(LUDEN,FILDEN,WORK(KDNIJ),ADR,LENDEN)
         ADR = ADR + LENDEN

C        Loop over operators A.
C        ----------------------

         DO IOPA = 1,NLRTLBL

            ISYMA = ISYLRT(IOPA)
            FREQA = FRQLRT(IOPA)

            DIFF = DABS(FREQA + FREQ)

            IF ((ISYMA.EQ.ISYMB) .AND. (DIFF.LE.THRDFF)) THEN

               LABA  = LRTLBL(IOPA)

               KINTA = KEND1 
               KAIJ  = KINTA + N2BST(ISYMA)
               KAAI  = KAIJ  + NMATIJ(ISYMA)
               KAIA  = KAAI  + NT1AM(ISYMA)
               KAAB  = KAIA  + NT1AM(ISYMA)
               KEND2 = KAAB  + NMATAB(ISYMA)
               LWRK2 = LWORK - KEND2 + 1

               IF (LWRK2 .LT. 0) THEN
                  CALL QUIT('Insufficient memory in '//SECNAM//' [2]')
               ENDIF

C              Read AO integrals.
C              ------------------
               
               IERR = -1
               CALL CCPRPAO(LABA,LONLYAO,WORK(KINTA),IRREP,IREAL,IERR,
     &                      WORK(KEND2),LWRK2)

C              Check.
C              ------

               IF (IERR .GT. 0) THEN
                  WRITE(LUPRI,'(//,5X,A,A)')
     &            SECNAM,': I/O error ocurred in CCPRPAO'
                  WRITE(LUPRI,'(5X,A,I10,A,/,5X,A,2X,A8,/,5X,A,I10)')
     &            'Operator number  : ',IOPA,' (in LRTLBL list)',
     &            'Operator label   : ',LABA,
     &            'Operator symmetry: ',ISYMA
                  WRITE(LUPRI,'(5X,A,I10,/)')
     &            'CCPRPAO returned IERR = ',IERR
                  CALL QUIT('I/O error in CCPRPAO ('//SECNAM//')')
               ELSE IF (IERR .LT. 0) THEN
                  WRITE(LUPRI,'(//,5X,A,A)')
     &          SECNAM,': CCPRPAO failed to determine operator symmetry'
                  WRITE(LUPRI,'(5X,A,I10,A,/,5X,A,2X,A8,/,5X,A,I10)')
     &            'Operator number  : ',IOPA,' (in LRTLBL list)',
     &            'Operator label   : ',LABA,
     &            'Operator symmetry: ',ISYMA
                  WRITE(LUPRI,'(5X,A,I10,/)')
     &            'CCPRPAO returned IERR = ',IERR
                  CALL QUIT('Irrep. error in CCPRPAO ('//SECNAM//')')
               ENDIF

               IF (IRREP .NE. ISYDEN) THEN
                  WRITE(LUPRI,'(//,5X,A,A)')
     &            'Symmetry mismatch in ',SECNAM
                  WRITE(LUPRI,'(5X,A,I10,A,/,5X,A,2X,A8)')
     &            'Operator number  : ',IOPA,' (in LRTLBL list)',
     &            'Operator label   : ',LABA
                  WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &            'CCPRPAO returned: ',IRREP,
     &            '- expected      : ',ISYMA
                  CALL QUIT('Symmetry mismatch in '//SECNAM)
               ENDIF

C              Transform to MO basis.
C              ----------------------

               CALL CC_CHOTRFOV(WORK(KLAMDP),WORK(KLAMDH),WORK(KINTA),
     &                          WORK(KAIJ),WORK(KAAB),WORK(KEND2),LWRK2,
     &                          1,1,ISYMA)
               CALL CC_CHOTRFED(WORK(KLAMDP),WORK(KLAMDH),WORK(KINTA),
     &                          WORK(KAAI),WORK(KAIA),WORK(KEND2),LWRK2,
     &                          1,1,ISYMA)

C              Calculate contribution and print.
C              ---------------------------------

               CON = DDOT(NMATIJ(ISYMA),WORK(KDNIJ),1,WORK(KAIJ),1)
     &             + DDOT(NT1AM(ISYMA),WORK(KDNIA),1,WORK(KAIA),1)
     &             + DDOT(NMATAB(ISYMA),WORK(KDNAB),1,WORK(KAAB),1)
               WRITE(LUPRI,1) LABA,ISYMA,LABB,ISYMB,FREQ,CON
    1          FORMAT(5X,A8,3X,I2,3X,A8,3X,I2,2X,F12.5,1X,1P,D22.15)
               CALL FLSHFO(LUPRI)

            ENDIF

         ENDDO

      ENDDO

C     Print end of table.
C     -------------------
      
      WRITE(LUPRI,'(5X,A,A,/)')
     & '--------------------------------------------',
     & '----------------------'
      CALL FLSHFO(LUPRI)

C     Close FILDEN.
C     -------------

      IF (DELDEN) THEN
         CALL WCLOSE2(LUDEN,FILDEN,'DELETE')
      ELSE
         CALL WCLOSE2(LUDEN,FILDEN,'KEEP')
      ENDIF

C     Stop if requested.
C     ------------------

      IF (STPDBG) THEN
         WRITE(LUPRI,*) SECNAM,': Stopping as requested!'
         CALL QUIT(' STOP AS REQUESTED ')
      ENDIF

      RETURN
      END
C  /* Deck cc_chodini */
      SUBROUTINE CC_CHODINI(IPDEN,NPDEN,LISTR,WORK,LWORK,FILDEN)
C
C     Thomas Bondo Pedersen, April 2003.
C
C     Purpose: Initialize densities on disk.
C
#include <implicit.h>
      INTEGER IPDEN(3,NPDEN)
      CHARACTER*(*) LISTR
      CHARACTER*10  FILDEN
      DIMENSION WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>
#include <ccr1rsp.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CC_CHODINI')

      INTEGER ADR
      INTEGER   ILSTSYM

      IF (NPDEN .LE. 0) RETURN

      LUDEN = 97
      CALL WOPEN2(LUDEN,FILDEN,64,0)

      ADR    = 1
      MAXDEN = -999
      DO IVEC = 1,NPDEN

         ILSTNR = IPDEN(2,IVEC)
         ISYVEC = ILSTSYM(LISTR,ILSTNR)
         LENDEN = NT1AM(ISYVEC) + NMATIJ(ISYVEC) + NMATAB(ISYVEC)

         IF (LENDEN .GT. MAXDEN) THEN
            MAXDEN = LENDEN
            IF (LWORK .LT. LENDEN) THEN
               CALL QUIT('Insufficient memory in '//SECNAM)
            ENDIF
            CALL DZERO(WORK,LENDEN)
         ENDIF

         CALL PUTWA2(LUDEN,FILDEN,WORK,ADR,LENDEN)
         ADR = ADR + LENDEN

      ENDDO

      CALL WCLOSE2(LUDEN,FILDEN,'KEEP')

      RETURN
      END
C  /* Deck cc_chor1d1dbg */
      SUBROUTINE CC_CHOR1D1DBG(XL2AM,XR2AM,WORK,LWORK,ISYMXL,ISYMXR)
C
C     Thomas Bondo Pedersen, April 2003.
C
C     Purpose: Debug CC_CHOR1D1.
C
#include <implicit.h>
      DIMENSION XL2AM(*), XR2AM(*), WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      CHARACTER*13 SECNAM
      PARAMETER (SECNAM = 'CC_CHOR1D1DBG')

      INTEGER AI

      ISYDEN = MULD2H(ISYMXL,ISYMXR)

      KDNJK = 1
      KDNCB = KDNJK + NMATIJ(ISYDEN)
      KEND1 = KDNCB + NMATAB(ISYDEN)
      LWRK1 = LWORK - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

      CALL DZERO(WORK(KDNJK),NMATIJ(ISYDEN))
      CALL DZERO(WORK(KDNCB),NMATAB(ISYDEN))

      DO ISYMAI = 1,NSYM

         ISYMBJ = MULD2H(ISYMAI,ISYMXR)
         ISYMCJ = MULD2H(ISYMAI,ISYMXL)
         ISYMBK = ISYMCJ

         DO AI = 1,NT1AM(ISYMAI)

            DO ISYMJ = 1,NSYM

               ISYMB = MULD2H(ISYMJ,ISYMBJ)
               ISYMC = MULD2H(ISYMJ,ISYMCJ)

               KOFFL = IT2SQ(ISYMCJ,ISYMAI)
     &               + NT1AM(ISYMCJ)*(AI - 1)
     &               + IT1AM(ISYMC,ISYMJ) + 1
               KOFFR = IT2SQ(ISYMBJ,ISYMAI)
     &               + NT1AM(ISYMBJ)*(AI - 1)
     &               + IT1AM(ISYMB,ISYMJ) + 1
               KOFFD = KDNCB + IMATAB(ISYMC,ISYMB)

               NTOTB = MAX(NVIR(ISYMB),1)
               NTOTC = MAX(NVIR(ISYMC),1)

               CALL DGEMM('N','T',NVIR(ISYMC),NVIR(ISYMB),NRHF(ISYMJ),
     &                    1.0D0,XL2AM(KOFFL),NTOTC,XR2AM(KOFFR),NTOTB,
     &                    1.0D0,WORK(KOFFD),NTOTC)

            ENDDO

            DO ISYMK = 1,NSYM

               ISYMB = MULD2H(ISYMK,ISYMBK)
               ISYMJ = MULD2H(ISYMB,ISYMBJ)

               KOFFR = IT2SQ(ISYMBJ,ISYMAI)
     &               + NT1AM(ISYMBJ)*(AI - 1)
     &               + IT1AM(ISYMB,ISYMJ) + 1
               KOFFL = IT2SQ(ISYMBK,ISYMAI)
     &               + NT1AM(ISYMBK)*(AI - 1)
     &               + IT1AM(ISYMB,ISYMK) + 1
               KOFFD = KDNJK + IMATIJ(ISYMJ,ISYMK)

               NTOTB = MAX(NVIR(ISYMB),1)
               NTOTJ = MAX(NRHF(ISYMJ),1)

               CALL DGEMM('T','N',NRHF(ISYMJ),NRHF(ISYMK),NVIR(ISYMB),
     &                    -1.0D0,XR2AM(KOFFR),NTOTB,XL2AM(KOFFL),NTOTB,
     &                    1.0D0,WORK(KOFFD),NTOTJ)

            ENDDO


         ENDDO

      ENDDO

      XJKNRM = DDOT(NMATIJ(ISYDEN),WORK(KDNJK),1,WORK(KDNJK),1)
      XCBNRM = DDOT(NMATAB(ISYDEN),WORK(KDNCB),1,WORK(KDNCB),1)

      CALL AROUND(SECNAM//': Occupied density')
      WRITE(LUPRI,*) 'Sq-norm: ',XJKNRM
      WRITE(LUPRI,*) '2--norm: ',DSQRT(XJKNRM)
c     CALL NOCC_PRT(WORK(KDNJK),ISYDEN,'IJ  ')

      CALL AROUND(SECNAM//': Virtual density')
      WRITE(LUPRI,*) 'Sq-norm: ',XCBNRM
      WRITE(LUPRI,*) '2--norm: ',DSQRT(XCBNRM)
c     CALL NOCC_PRT(WORK(KDNCB),ISYDEN,'AB  ')

      RETURN
      END
C  /* Deck cc_chor1ds */
      SUBROUTINE CC_CHOR1DS(IPDEN,NPDEN,LISTR,WORK,LWORK,MODEL,FILDEN)
C
C     Thomas Bondo Pedersen, April 2003.
C
C     Purpose: Calculate singles/0th order density contributions to
C              the 1st order right density.
C
C     Formulas:
C
C     D(jk) = - sum(b) R(bj) * D0(bk)
C     D(jb) =   sum(i) D0(ji) * R(bi)
C             - sum(a) R(aj) * D0(ab)
C     D(cb) =   sum(i) D0(ci) * R(bi)
C
C     where D0 is the 0th order density which must be available on disk
C     (D0(ai) = tbar(ai)).
C
C     The MO blocks are stored consecutively on file FILDEN (via long
C     integer crayio). Thus, the ordering of this file is (for each amplitude
C     in LISTR as specified in IPDEN(2,*)):
C
C        D(jk),D(jb),D(cb).
C
#include <implicit.h>
      INTEGER IPDEN(3,NPDEN)
      CHARACTER*(*) LISTR
      CHARACTER*10  MODEL, FILDEN
      DIMENSION WORK(LWORK)
#include <maxorb.h>
#include <ccdeco.h>
#include <ccisvi.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <dccorb.h>
#include <dccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>
#include <dummy.h>
#include <ciarc.h>
#include <ccr1rsp.h>
#include <chocc2.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CC_CHOR1DS')

      INTEGER ADR

      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0)

      INTEGER ILSTSYM

C     Return if nothing to do.
C     ------------------------

      IF (NPDEN .LE. 0) RETURN

C     Allocation.
C     -----------

      KDOCC = 1
      KDVIR = KDOCC + NMATIJ(1)
      KDEXC = KDVIR + NMATAB(1)
      KEND1 = KDEXC + NT1AM(1)
      LWRK1 = LWORK - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: D0'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need (more than): ',KEND1-1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Read 0th order density blocks.
C     ------------------------------

      LUDN0 = 96
      CALL WOPEN2(LUDN0,FDEN00,64,0)
      ADR = 1
      CALL GETWA2(LUDN0,FDEN00,WORK(KDOCC),ADR,NMATIJ(1))
      ADR = ADR + NMATIJ(1) + NT1AM(1)
      CALL GETWA2(LUDN0,FDEN00,WORK(KDVIR),ADR,NMATAB(1))
      CALL WCLOSE2(LUDN0,FDEN00,'KEEP')
      IOPT = 1
      CALL CC_RDRSP('L0 ',1,1,IOPT,MODEL,WORK(KDEXC),DUMMY)

C     Open FILDEN.
C     ------------

      LUDEN = 97
      CALL WOPEN2(LUDEN,FILDEN,64,0)

C     Initialize ADR on FILDEN.
C     -------------------------

      ADR = 1

C     Loop over amplitudes in LISTR.
C     ------------------------------

      DO IVEC = 1,NPDEN

C        Get list info.
C        --------------

         ILSTNR = IPDEN(2,IVEC)
         ISYVEC = ILSTSYM(LISTR,ILSTNR)
         ISYDEN = ISYVEC

C        Allocation.
C        -----------

         KR1AM = KEND1
         KDNJK = KR1AM + NT1AM(ISYVEC)
         KDNJB = KDNJK + NMATIJ(ISYDEN)
         KDNCB = KDNJB + NT1AM(ISYDEN)
         KEND2 = KDNCB + NMATAB(ISYDEN)
         LWRK2 = LWORK - KEND2 + 1

         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Insufficient memory in ',SECNAM,' - allocation: D1'
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need     : ',KEND2-1,
     &      'Available: ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Read singles LISTR amplitudes.
C        ------------------------------

         IOPT = 1
         CALL CC_RDRSP(LISTR,ILSTNR,ISYVEC,IOPT,MODEL,WORK(KR1AM),DUMMY)

C        Calculate occupied block.
C        -------------------------

         DO ISYMK = 1,NSYM

            ISYMJ = MULD2H(ISYMK,ISYDEN)
            ISYMB = MULD2H(ISYMJ,ISYVEC)

            KOFFR = KR1AM + IT1AM(ISYMB,ISYMJ)
            KOFF0 = KDEXC + IT1AM(ISYMB,ISYMK)
            KOFFD = KDNJK + IMATIJ(ISYMJ,ISYMK)

            NTOTB = MAX(NVIR(ISYMB),1)
            NTOTJ = MAX(NRHF(ISYMJ),1)

            CALL DGEMM('T','N',NRHF(ISYMJ),NRHF(ISYMK),NVIR(ISYMB),
     &                 XMONE,WORK(KOFFR),NTOTB,WORK(KOFF0),NTOTB,
     &                 ZERO,WORK(KOFFD),NTOTJ)

         ENDDO

C        Calculate virtual block.
C        ------------------------

         DO ISYMI = 1,NSYM

            ISYMB = MULD2H(ISYMI,ISYVEC)
            ISYMC = MULD2H(ISYMB,ISYDEN)

            KOFF0 = KDEXC + IT1AM(ISYMC,ISYMI)
            KOFFR = KR1AM + IT1AM(ISYMB,ISYMI)
            KOFFD = KDNCB + IMATAB(ISYMC,ISYMB)

            NTOTB = MAX(NVIR(ISYMB),1)
            NTOTC = MAX(NVIR(ISYMC),1)

            CALL DGEMM('N','T',NVIR(ISYMC),NVIR(ISYMB),NRHF(ISYMI),
     &                 ONE,WORK(KOFF0),NTOTC,WORK(KOFFR),NTOTB,
     &                 ZERO,WORK(KOFFD),NTOTC)

         ENDDO

C        Calculate de-excitation block.
C        ------------------------------

         DO ISYMI = 1,NSYM

            ISYMB = MULD2H(ISYMI,ISYVEC)
            ISYMJ = MULD2H(ISYMB,ISYDEN)

            KOFFR = KR1AM + IT1AM(ISYMB,ISYMI)
            KOFF0 = KDOCC + IMATIJ(ISYMJ,ISYMI)
            KOFFD = KDNJB + IT1AM(ISYMB,ISYMJ)

            NTOTB = MAX(NVIR(ISYMB),1)
            NTOTJ = MAX(NRHF(ISYMJ),1)

            CALL DGEMM('N','T',NVIR(ISYMB),NRHF(ISYMJ),NRHF(ISYMI),
     &                 ONE,WORK(KOFFR),NTOTB,WORK(KOFF0),NTOTJ,
     &                 ZERO,WORK(KOFFD),NTOTB)

         ENDDO

         DO ISYMJ = 1,NSYM

            ISYMA = MULD2H(ISYMJ,ISYVEC)
            ISYMB = MULD2H(ISYMJ,ISYDEN)

            KOFF0 = KDVIR + IMATAB(ISYMA,ISYMB)
            KOFFR = KR1AM + IT1AM(ISYMA,ISYMJ)
            KOFFD = KDNJB + IT1AM(ISYMB,ISYMJ)

            NTOTA = MAX(NVIR(ISYMA),1)
            NTOTB = MAX(NVIR(ISYMB),1)

            CALL DGEMM('T','N',NVIR(ISYMB),NRHF(ISYMJ),NVIR(ISYMA),
     &                 XMONE,WORK(KOFF0),NTOTA,WORK(KOFFR),NTOTA,
     &                 ONE,WORK(KOFFD),NTOTB)

         ENDDO

C        Write density to disk.
C        ----------------------

         LDEN = NMATIJ(ISYDEN) + NT1AM(ISYDEN) + NMATAB(ISYDEN)
         CALL PUTWA2(LUDEN,FILDEN,WORK(KDNJK),ADR,LDEN)

C        Update ADR.
C        -----------

         ADR = ADR + LDEN

      ENDDO

C     Close FILDEN.
C     -------------

      CALL WCLOSE2(LUDEN,FILDEN,'KEEP')

      RETURN
      END
