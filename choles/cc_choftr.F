C  /* Deck cc_chofmat */
      SUBROUTINE CC_CHOFMAT(IFTRAN,NFTRAN,LISTL,LISTR,IOPTRES,
     &                      FILFMA,IFDOTS,FCONS,MXVEC,WORK,LWORK)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Drive the Cholesky CC2 F matrix transformations.
C              The interface is identical to that of the conventional
C              CC_FMATRIX routine by Christof Haettig, although this
C              routine is *far* less flexible. See CC_FMATRIX/CC_FMATNEW
C              for a detailed description of the input.
C
#include <implicit.h>
      INTEGER IFTRAN(3,NFTRAN), IFDOTS(MXVEC,NFTRAN)
      CHARACTER*(*) LISTL, LISTR, FILFMA
      DIMENSION FCONS(MXVEC,NFTRAN), WORK(LWORK)
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CC_CHOFMAT')

C     Return if nothing to do.
C     ------------------------

      IF (NFTRAN .LE. 0) RETURN

C     Check that this is Cholesky CC2.
C     --------------------------------

      IF (.NOT. CHOINT) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   SECNAM,': integrals must be decomposed!'
         CALL QUIT('Integrals not decomposed in '//SECNAM)
      ENDIF

      IF (.NOT. CC2) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   SECNAM,': model is not CC2!'
         CALL QUIT('Only CC2 implemented in '//SECNAM)
      ENDIF

C     Call the appropriate routine.
C     -----------------------------

      IF (LISTL(1:2) .EQ. 'L0') THEN
         CALL CC_CHOFTR0R(IFTRAN,NFTRAN,LISTR,IOPTRES,FILFMA,
     &                    IFDOTS,FCONS,MXVEC,WORK,LWORK)
      ELSE
         WRITE(LUPRI,'(//,5X,A,A,A,A)')
     &   SECNAM,': LISTL = ',LISTL(1:3),' not implemented.'
         CALL QUIT('List '//LISTL(1:3)//' not implemented in '//SECNAM)
      ENDIF

      RETURN
      END
C  /* Deck cc_choftr0r */
      SUBROUTINE CC_CHOFTR0R(IFTRAN,NFTRAN,LISTR,IOPTRES,FILFMA,
     &                       IFDOTS,FCONS,MXVEC,WORK,LWORK)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Calculate Cholesky CC2 F matrix transformations,
C              singles only. Result vectors are stored on disk
C              according to FILFMA (via CC_WRRSP) for IOPTRES=3.
C              The partial result array F12*R2 is also stored on disk.
C              For IOPTRES=5, the following dot products are calculated:
C                 d = (F*R)*R'
C                   = (F11*R1 + F12*R2)*R1'
C                   + (F21*R1)*R2'
C                   = (F11*R1 + F12*R2)*R1'
C                   + (F12*R2')*R1
C              and returned in FCONS (as in the conventional routine). Note
C              that the F transformed vectors as well as the F12*R2 partial
C              result *must* be available on disk for this option.
C
C     TODO/FIXME: dot products should also be available "on-the-fly"....
C
#include <implicit.h>
      INTEGER IFTRAN(3,NFTRAN), IFDOTS(MXVEC,NFTRAN)
      CHARACTER*(*) LISTR, FILFMA
      DIMENSION FCONS(MXVEC,NFTRAN), WORK(LWORK)
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>
#include <dummy.h>
#include <chocc2.h>
#include <ccr1rsp.h>

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CC_CHOFTR0R')

      LOGICAL LONLYAO
      DATA LONLYAO /.TRUE./

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      PARAMETER (INFO = 3)
      PARAMETER (ITYP1 = 1, ITYP2 = 2, ITYP3 = 3, ITYP4 = 4)
      PARAMETER (ITYP5 = 5, ITYP6 = 6, ITYP7 = 7, ITYP8 = 8)
      PARAMETER (ITYP9 = 9, ITYP10 = 10, ITYP11 = 11)
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)

      CHARACTER*3  FLXFMA
      CHARACTER*8  LABEL
      CHARACTER*10 MODEL

      INTEGER ILSTSYM

C     Start timing.
C     -------------

      TIMTOT = SECOND()
      TIMMOL = ZERO
      TIMMOR = ZERO
      TIMYIR = ZERO
      TIMYIL = ZERO
      TIMCON = ZERO

C     Redirect for IOPTRES=5.
C     -----------------------

      IF (IOPTRES .EQ. 5) THEN
         CALL DZERO(FCONS,MXVEC*NFTRAN)
         DO ITRAN = 1,NFTRAN
            IFTRAN(3,ITRAN) = IFTRAN(2,ITRAN)
         ENDDO
         CALL CC_CHOFTR5(IFTRAN,NFTRAN,LISTR,IOPTRES,FILFMA,
     &                   IFDOTS,FCONS,MXVEC,WORK,LWORK)
         RETURN
      ELSE IF (IOPTRES .NE. 3) THEN
         WRITE(LUPRI,'(//,5X,A,A,I4,A)')
     &   SECNAM,': IOPTRES =',IOPTRES,' not recognized.'
         CALL QUIT('Illegal IOPTRES in '//SECNAM)
      ENDIF

C     Check that LISTR has been implemented.
C     --------------------------------------

      IF (LISTR(1:2) .NE. 'R1') THEN
         WRITE(LUPRI,'(//,5X,A,A,A,A)')
     &   SECNAM,': LISTR = ',LISTR(1:3),' not implemented.'
         CALL QUIT('List '//LISTR(1:3)//' not implemented in '//SECNAM)
      ENDIF

C     Set model.
C     (It should have been checked that this is CC2 in calling routine.)
C     ------------------------------------------------------------------

      MODEL = 'CC2       '

C     Set file list for F12*R2 contribution.
C     NB! this will cause trouble if FILFMA contains a 3-character identifier!
C     ------------------------------------------------------------------------

      FLXFMA(1:1) = 'X'
      FLXFMA(2:2) = FILFMA(1:1)
      FLXFMA(3:3) = FILFMA(2:2)

C     Print.
C     ------

      IF (IPRINT .GE. INFO) THEN
         CALL AROUND('Output from '//SECNAM)
         WRITE(LUPRI,'(5X,A,A,A,I10)')
     &   'Number of ',MODEL,' F-matrix transformations requested:',
     &   NFTRAN
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
            DO IVEC = 1,NFTRAN
               ILSTNR = IFTRAN(2,IVEC)
               ISYMR1 = ILSTSYM(LISTR,ILSTNR)
               LABEL  = LRTLBL(ILSTNR)
               FREQ   = FRQLRT(ILSTNR)
               WRITE(LUPRI,'(8X,I6,4X,I4,6X,A8,2X,F12.5)')
     &         ILSTNR,ISYMR1,LABEL,FREQ
            ENDDO
            WRITE(LUPRI,'(8X,A,/)')
     &      '------------------------------------------'
         ELSE
            WRITE(LUPRI,'(5X,A)')
     &      'Sorry! This type has not yet been implemented!'
            CALL QUIT('Error in '//SECNAM)
         ENDIF
         CALL FLSHFO(LUPRI)
      ENDIF

C     Allocation.
C     -----------

      KFOCKD = 1
      KEND0  = KFOCKD + NORBTS
      LWRK0  = LWORK  - KEND0 + 1

      IF (LWRK0 .LT. 0) THEN
         CALL QUIT('Insufficient memory in '//SECNAM//' [0]')
      ENDIF

C     Read canonical orbital energies.
C     --------------------------------

      CALL CHO_RDSIR(DUMMY,DUMMY,WORK(KFOCKD),DUMMY,WORK(KEND0),LWRK0,
     &               .FALSE.,.TRUE.,.FALSE.)

C     Make sure that the 0th order doubles amplitudes have
C     been Cholesky decomposed.
C     ----------------------------------------------------

      IF (.NOT. CHOT2C) THEN
         CALL CC_CHOCIM1(WORK(KFOCKD),DUMMY,DUMMY,WORK(KEND0),LWRK0,0,0)
      ENDIF

C     Allocation.
C     -----------

      KFIA   = KEND0
      KLAMDP = KFIA   + NT1AM(1)
      KLAMDH = KLAMDP + NGLMDT(1)
      KTBAR  = KLAMDH + NGLMDT(1)
      KLAM2  = KTBAR  + NT1AM(1)
      KEND1  = KLAM2  + NGLMDT(1)
      LWRK1  = LWORK  - KEND1 + 1

      KT1AM = KTBAR

      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient memory in '//SECNAM//' [1]')
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
     &            WORK(KEND1),LWRK1)

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
     &   SECNAM,' after reading TBAR vector'
         WRITE(LUPRI,*) '   Square norm TBAR: ',T1NM2
         WRITE(LUPRI,*) '   2-norm      TBAR: ',T1NRM
         CALL FLSHFO(LUPRI)
      ENDIF

C     Calculate Ltilde(ia,J) Cholesky vectors.
C     ----------------------------------------

      TIMMOL = SECOND()

      CALL CC_CHOTG(WORK(KLAMDP),WORK(KLAMDH),WORK(KTBAR),
     &              WORK(KEND1),LWRK1,1,1,-1)

C     Add: L(ia,J) + Ltilde(ia,J), store in ITYP8 files.
C     --------------------------------------------------

      CALL CC_CHOADD(ITYP1,ITYP8,1,NUMCHO,WORK(KEND1),LWRK1,ONE,ONE)

      TIMMOL = SECOND() - TIMMOL

C     Print.
C     ------

      IF (IPRINT .GE. INFO) THEN
         WRITE(LUPRI,'(5X,A,/,5X,A)')
     &   'List                 Time (seconds)',
     &   ' #    L Y Interm.  R Y Interm.  R Ch. Prep.  Contributions'
         WRITE(LUPRI,'(5X,A)')
     &   '----------------------------------------------------------'
         CALL FLSHFO(LUPRI)
      ENDIF

C     Start loop over vectors to transform.
C     =====================================

      DO IR = 1,NFTRAN

C        Time this transformation.
C        -------------------------

         TIMVEC = SECOND()
         DTMMOR = ZERO
         DTMYIR = ZERO
         DTMYIL = ZERO
         DTMCON = ZERO

C        Get vector/operator info.
C        -------------------------

         ILSTNR = IFTRAN(2,IR)
         ISYMR1 = ILSTSYM(LISTR,ILSTNR)
         LABEL  = LRTLBL(ILSTNR)
         FREQ   = FRQLRT(ILSTNR)

C        Read AO property integrals corresponding to operator LABEL.
C        Check info.
C        -----------------------------------------------------------

         KRES  = KEND1
         KRESP = KRES  + NT1AM(ISYMR1)
         KRESI = KRESP + NT1AM(ISYMR1)
         KXIJ  = KRESI + NT1AM(ISYMR1)
         KXAB  = KXIJ  + NMATIJ(ISYMR1)
         KXAO  = KXAB  + NMATAB(ISYMR1)
         KEND2 = KXAO  + N2BST(ISYMR1)
         LWRK2 = LWORK - KEND2 + 1

         IF (LWRK2 .LT. 0) THEN
            CALL QUIT('Insufficient memory in '//SECNAM//' [2]')
         ENDIF

         IERR = -1
         CALL CCPRPAO(LABEL,LONLYAO,WORK(KXAO),IRREP,IREAL,IERR,
     &                WORK(KEND2),LWRK2)

         IF (IERR .GT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      SECNAM,': I/O error ocurred in CCPRPAO'
            WRITE(LUPRI,'(5X,A,I10,A,/,5X,A,2X,A8,/,5X,A,I10)')
     &      'Operator number  : ',ILSTNR,' (in LRTLBL list)',
     &      'Operator label   : ',LABEL,
     &      'Operator symmetry: ',ISYMR1
            WRITE(LUPRI,'(5X,A,I10,/)')
     &      'CCPRPAO returned IERR = ',IERR
            CALL QUIT('I/O error in CCPRPAO ('//SECNAM//')')
         ELSE IF (IERR .LT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      SECNAM,': CCPRPAO failed to determine operator symmetry'
            WRITE(LUPRI,'(5X,A,I10,A,/,5X,A,2X,A8,/,5X,A,I10)')
     &      'Operator number  : ',ILSTNR,' (in LRTLBL list)',
     &      'Operator label   : ',LABEL,
     &      'Operator symmetry: ',ISYMR1
            WRITE(LUPRI,'(5X,A,I10,/)')
     &      'CCPRPAO returned IERR = ',IERR
            CALL QUIT('Irrep. error in CCPRPAO ('//SECNAM//')')
         ENDIF

         IF (IRREP .NE. ISYMR1) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Symmetry mismatch in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,A,/,5X,A,2X,A8)')
     &      'Operator number  : ',ILSTNR,' (in LRTLBL list)',
     &      'Operator label   : ',LABEL
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'CCPRPAO returned: ',IRREP,
     &      '- expected      : ',ISYMR1
            CALL QUIT('Symmetry mismatch in '//SECNAM)
         ENDIF

C        Transform integrals: occupied and virtual blocks only.
C        ------------------------------------------------------

         CALL CC_CHOTRFOV(WORK(KLAMDP),WORK(KLAMDH),WORK(KXAO),
     &                    WORK(KXIJ),WORK(KXAB),WORK(KEND2),LWRK2,
     &                    1,1,ISYMR1)

C        Reset memory: AO integrals no longer needed.
C        --------------------------------------------

         KEND2 = KXAO
         LWRK2 = LWORK - KEND2 + 1

C        Transform Cholesky vectors corresponding to 0th order
C        amplitudes, incorporating the minus sign. Store result
C        in ITYP7 files.
C        ------------------------------------------------------

         DTIME = SECOND()
         CALL CC_CHOTRPRP(WORK(KXIJ),WORK(KXAB),WORK(KEND2),LWRK2,
     &                    ISYMR1,ITYP6,ITYP7,1,NCHMOC)
         DTIME = SECOND() - DTIME
         TIMMOR = TIMMOR  + DTIME
         DTMMOR = DTMMOR  + DTIME

C        Allocation.
C        -----------

         KR1AM = KXIJ
         KLAM2 = KR1AM + NT1AM(ISYMR1)
         KEND3 = KLAM2 + NGLMDT(ISYMR1)
         LWRK3 = LWORK - KEND3 + 1

         IF (LWRK3 .LT. 0) THEN
            CALL QUIT('Insufficient memory in '//SECNAM//' [3]')
         ENDIF

C        Read R1 amplitudes.
C        -------------------

         IOPT = 1
         CALL CC_RDRSP(LISTR,ILSTNR,ISYMR1,IOPT,MODEL,WORK(KR1AM),DUMMY)

         IF (LOCDBG) THEN
            R1NM2 = DDOT(NT1AM(ISYMR1),WORK(KR1AM),1,WORK(KR1AM),1)
            R1NRM = DSQRT(R1NM2)
            WRITE(LUPRI,*)
     &      SECNAM,' after reading R vector number: ',ILSTNR
            WRITE(LUPRI,*) '   symmetry        : ',ISYMR1
            WRITE(LUPRI,*) '   assoc. operator : ',LABEL 
            WRITE(LUPRI,*) '   assoc. frequency: ',FREQ
            WRITE(LUPRI,*) '   Square norm R1AM: ',R1NM2
            WRITE(LUPRI,*) '   2-norm      R1AM: ',R1NRM
            CALL FLSHFO(LUPRI)
         ENDIF

C        Calculate perturbed lambda matrix.
C        ----------------------------------

         DTIME = SECOND()

         CALL CC_LAMPT(WORK(KLAMDP),WORK(KLAMDH),WORK(KR1AM),
     &                 WORK(KLAM2),1,ISYMR1,1)

C        Transform Cholesky vectors: Lbar(ai,J), Lbar(ij,J). Store in
C        ITYP9 and ITYP10 files, respectively.
C        ------------------------------------------------------------

         CALL CC_CHOTG2(WORK(KLAMDP),WORK(KLAMDH),WORK(KLAM2),
     &                  WORK(KEND3),LWRK3,1,ISYMR1,ITYP9,ITYP10)

         DTIME = SECOND() - DTIME
         TIMMOR = TIMMOR  + DTIME
         DTMMOR = DTMMOR  + DTIME

C        Allocation.
C        -----------

         KCBAR = KEND3
         KEND4 = KCBAR + NT1AM(ISYMR1)
         LWRK4 = LWORK - KEND4 + 1

         IF (LWRK4 .LT. 0) THEN
            CALL QUIT('Insufficient memory in '//SECNAM//' [4]')
         ENDIF

C        Calculate Y3 and Cbar intermediates:
C        Y3(ai,J) = sum(bj) [2*R(ai,bj)-R(aj,bi)] * L(jb,J)
C        Cbar(ai) = sum(bj) [2*R(ai,bj)-R(aj,bi)] * TBAR(bj).
C        Y3 intermediates stored in files identified by ITYPY=4.
C        -------------------------------------------------------

         DTIME = SECOND()
         CALL CC_CYIFTRR(LISTR,WORK(KFOCKD),WORK(KTBAR),WORK(KCBAR),
     &                   WORK(KEND4),LWRK4,ISYMR1,1,FREQ,.TRUE.)
         DTIME = SECOND() - DTIME
         TIMYIR = TIMYIR  + DTIME
         DTMYIR = DTMYIR  + DTIME

c     cbnrm = dsqrt(ddot(nt1am(isymr1),work(kcbar),1,work(kcbar),1))
c     write(lupri,*) SECNAM,': 2-norm of Cbar interm.: ',cbnrm

C        Calculate Ebar intermediates and the reference contribution,
C        stored in KRREF. Save the Y3-contributions separately. Delete
C        Y3 files after use.
C        -------------------------------------------------------------

         KRREF = KEND4
         KEIJ  = KRREF + NT1AM(ISYMR1)
         KEAB  = KEIJ  + 2*NMATIJ(ISYMR1)
         KEND5 = KEAB  + 2*NMATAB(ISYMR1)
         LWRK5 = LWORK - KEND5 + 1

         IF (LWRK5 .LT. 0) THEN
            CALL QUIT('Insufficient memory in '//SECNAM//' [5]')
         ENDIF

         DTIME  = SECOND()

         ITYPY3 = 4
         ISYMY3 = ISYMR1
         CALL CC_CHOEBAR(WORK(KLAMDP),WORK(KLAMDH),WORK(KFIA),
     &                   WORK(KR1AM),WORK(KRREF),WORK(KEIJ),WORK(KEAB),
     &                   WORK(KEND5),LWRK5,ISYMR1,ITYPY3,ISYMY3,.TRUE.,
     &                   .TRUE.,.TRUE.)

C        Calculate contributions from Ebar intermediates, and
C        add reference contribution.
C        ----------------------------------------------------

         CALL CC_CHOATR0M1(WORK(KRES),WORK(KTBAR),WORK(KEIJ),WORK(KEAB),
     &                     1,ISYMR1,1)
         CALL DAXPY(NT1AM(ISYMR1),ONE,WORK(KRREF),1,WORK(KRES),1)

c     rfnm2 = ddot(nt1am(isymr1),work(krref),1,work(krref),1)
c     write(lupri,*) SECNAM,': Sq-norm of ref. contr.: ',rfnm2

C        Calculate Ebar contributions to F12*R2 partial result.
C        ------------------------------------------------------

         KOFF1 = KEIJ + NMATIJ(ISYMR1)
         KOFF2 = KEAB + NMATAB(ISYMR1)
         CALL CC_CHOATR0M1(WORK(KRESP),WORK(KTBAR),
     &                     WORK(KOFF1),WORK(KOFF2),
     &                     1,ISYMR1,1)

         DTIME  = SECOND() - DTIME
         TIMCON = TIMCON   + DTIME
         DTMCON = DTMCON   + DTIME

C        Calculate Y1 and Y2 intermediates, the latter only in first
C        call.
C        Y1(ai,J) = sum(bj) tbar(ai,bj) * Lbar(bj,J), store in ITYPY=2
C        Y2(ai,J) = sum(bj) tbar(ai,bj) * L(bj,J),    store in ITYPY=33.
C        Use work space from KRREF, overwriting
C        the Ebar intermediates and the ref. contr. 2*Fbar(ia).
C        ---------------------------------------------------------------

         KEND6 = KRREF
         LWRK6 = LWORK - KEND6 + 1

         DTIME = SECOND()
         CALL CC_CYIFTRL('L0 ',WORK(KFOCKD),WORK(KFIA),WORK(KTBAR),
     &                   WORK(KEND6),LWRK6,1,ISYMR1,DUMMY)
         DTIME = SECOND() - DTIME
         TIMYIL = TIMYIL  + DTIME
         DTMYIL = DTMYIL  + DTIME

C        Calculate GJ and H contributions.
C        ---------------------------------

         DTIME = SECOND()

         CALL CC_CHOFTRGJH(WORK(KLAMDP),WORK(KLAMDH),WORK(KLAM2),
     &                     WORK(KRES),WORK(KTBAR),
     &                     WORK(KEND6),LWRK6,ISYMR1,1)

C        Calculate I term.
C        -----------------

         CALL DZERO(WORK(KRESI),NT1AM(ISYMR1))
         CALL CC_CHOFTRI(WORK(KRESI),WORK(KCBAR),WORK(KEND6),LWRK6,
     &                   ISYMR1,ITYP1,1,NUMCHO)

C        Add I term to results.
C        ----------------------

         CALL DAXPY(NT1AM(ISYMR1),ONE,WORK(KRESI),1,WORK(KRES),1)
         CALL DAXPY(NT1AM(ISYMR1),ONE,WORK(KRESI),1,WORK(KRESP),1)

         DTIME = SECOND() - DTIME
         TIMCON = TIMCON  + DTIME
         DTMCON = DTMCON  + DTIME

C        Write results to disk.
C        ----------------------

         KEND = KRESP  + NT1AM(ISYMR1)
         LWRK = LWORK - KEND + 1

         KDUM = -99999

         IF (LOCDBG) THEN
            R1NM2 = DDOT(NT1AM(ISYMR1),WORK(KR1AM),1,WORK(KR1AM),1)
            R1NRM = DSQRT(R1NM2)
            S1NM2 = DDOT(NT1AM(ISYMR1),WORK(KRES),1,WORK(KRES),1)
            S1NRM = DSQRT(S1NM2)
            P1NM2 = DDOT(NT1AM(ISYMR1),WORK(KRESP),1,WORK(KRESP),1)
            P1NRM = DSQRT(P1NM2)
            WRITE(LUPRI,*)
     &      SECNAM,' before writing result vector number: ',
     &      IFTRAN(3,IR)
            WRITE(LUPRI,*) '   to file list    : ',FILFMA
            WRITE(LUPRI,*) '   Symmetry        : ',ISYMR1
            WRITE(LUPRI,*) '   R vector number : ',ILSTNR
            WRITE(LUPRI,*) '   Sq. norm, R vec.: ',R1NM2
            WRITE(LUPRI,*) '   2-norm,   R vec.: ',R1NRM
            WRITE(LUPRI,*) '   Sq. norm, result: ',S1NM2
            WRITE(LUPRI,*) '   2-norm,   result: ',S1NRM
            WRITE(LUPRI,*) '   Sq. norm, p. res: ',P1NM2
            WRITE(LUPRI,*) '   2-norm,   p. res: ',P1NRM
            CALL FLSHFO(LUPRI)
c           CALL HEADER('Result Vector',-1)
c           CALL NOCC_PRT(WORK(KRES),ISYMR1,'AI  ')
c           CALL HEADER('Partial Result Vector',-1)
c           CALL NOCC_PRT(WORK(KRESP),ISYMR1,'AI  ')
         ENDIF

         IOPT = 1
         CALL CC_WRRSP(FILFMA,IFTRAN(3,IR),ISYMR1,IOPT,MODEL,
     &                 WORK(KDUM),WORK(KRES),DUMMY,WORK(KEND),LWRK)

         IOPT = 1
         CALL CC_WRRSP(FLXFMA,IFTRAN(3,IR),ISYMR1,IOPT,MODEL,
     &                 WORK(KDUM),WORK(KRESP),DUMMY,WORK(KEND),LWRK)

C        Print timings.
C        --------------

         IF (IPRINT .GE. INFO) THEN
            WRITE(LUPRI,'(I9,2X,F11.2,2X,F11.2,2X,F11.2,4X,F11.2)')
     &      ILSTNR,DTMYIL,DTMYIR,DTMMOR,DTMCON
            CALL FLSHFO(LUPRI)
         ENDIF

      ENDDO

C     Print.
C     ------

      IF (IPRINT .GE. INFO) THEN
         WRITE(LUPRI,'(5X,A,/)')
     &   '----------------------------------------------------------'
         TIMTOT = SECOND() - TIMTOT
         CALL HEADER('Timing of '//SECNAM,-1)
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for left  Cholesky vector preparation: ',TIMMOL,
     &   ' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for right Cholesky vector preparation: ',TIMMOR,
     &   ' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for left  Y intermediates            : ',TIMYIL,
     &   ' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for right Y intermediates            : ',TIMYIR,
     &   ' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for calculating contributions        : ',TIMCON,
     &   ' seconds'
         WRITE(LUPRI,'(5X,A,A)')
     &   '-----------------------------------------------------------',
     &   '--------'
         WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &   'Time used in total                             : ',TIMTOT,
     &   ' seconds'
         CALL FLSHFO(LUPRI)
      ENDIF

      RETURN
      END
C  /* Deck cc_cyiftrr */
      SUBROUTINE CC_CYIFTRR(LISTR,FOCKD,X1AM,CIM,WORK,LWORK,
     &                      ISYMR,ISYMX,FREQ,INICIM)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Calculate Y intermediates from the "right-hand vectors"
C              specified through LIST for F matrix transformation,
C              calculating the C intermediate as a byproduct.
C
C     Formulas:
C
C        Y3(ai,J) = sum(bj) [2*R(ai,bj) - R(aj,bi)] * L(jb,J)
C
C        C(ai)    = C(ai) + sum(bj) [2*R(ai,bj) - R(aj,bi)] * X1AM(bj)
C
C     where the "right-hand vectors" are identified through LISTR.
C
C     IF (INICIM): C intermediates will be initialized here.
C
C     Implemented lists:
C
C        LISTR(1:2) = 'R1' -- 1st order amplitudes.
C
C
C     NOTE: all Cholesky vectors needed must be available on disk.
C           No checking is done here.
C
#include <implicit.h>
      CHARACTER*(*) LISTR
      DIMENSION FOCKD(*), X1AM(*), CIM(*), WORK(LWORK)
      LOGICAL INICIM
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>
#include <chocc2.h>
#include <dummy.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CC_CYIFTRR')

      PARAMETER (MXTYP = 2)
      DIMENSION SCD(MXTYP)
      INTEGER JTYP1(MXTYP), JTYP2(MXTYP), ISYMP1(MXTYP), ISYMP2(MXTYP)
      INTEGER NCHP12(8,MXTYP)
      LOGICAL DELP1(MXTYP), DELP2(MXTYP)

      PARAMETER (MXTPZ = 1)
      INTEGER JTYPZ(MXTPZ), ISYMZ(MXTPZ), JTYPY(MXTPZ)
      INTEGER NCHOZ(8,MXTPZ)

      PARAMETER (MXNUMW = 1)
      INTEGER ISYMW(MXNUMW)

      IF (LISTR(1:2) .EQ. 'R1') THEN

C        First order amplitudes.
C        =======================

C        Scaling factors.
C        ----------------

         FAC1 = 1.0D0
         FAC2 = 0.0D0
         SCDG = 1.0D0

C        Incorporate denominators and set up 2CME.
C        -----------------------------------------

         IOPTDN = 1
         IOPTCE = 1

C        Two Cholesky vector contractions needed (NUMP12).
C        No "extra-terms" (NUMUV).
C        One set of Y intermediates (NUMZ).
C        One C intermediate (NUMW).
C        ---------------------------------------------------

         NUMP12 = 2
         NUMUV  = 0
         NUMZ   = 1
         NUMW   = 1

C        First Cholesky vector contraction:
C        JTYP=6: M(ai,K) from 0th order amplitude decomp. [repr. -t]
C        JTYP=7: transformed M(ai,K) [incorporating the minus sign]
C                stored in dummy files.
C        Keep 6, delete 7 after use.
C        --------------------------------------------------------------

         JTYP1(1)  = 6
         JTYP2(1)  = 7
         ISYMP1(1) = 1
         ISYMP2(1) = ISYMR
         DO ISYCHO = 1,NSYM
            NCHP12(ISYCHO,1) = NCHMOC(ISYCHO)
         ENDDO
         DELP1(1) = .FALSE.
         DELP2(1) = .TRUE.
         SCD(1)   = 1.0D0

C        Second Cholesky vector contraction:
C        JTYP=3: L(ai,J), stored permanently from ground state calculation.
C        JTYP=9: Lbar(ai,J), "perturbed" Cholesky vectors.
C        Keep both after use.
C        ------------------------------------------------------------------

         JTYP1(2)  = 3
         JTYP2(2)  = 9
         ISYMP1(2) = 1
         ISYMP2(2) = ISYMR
         DO ISYCHO = 1,NSYM
            NCHP12(ISYCHO,2) = NUMCHO(ISYCHO)
         ENDDO
         DELP1(2) = .FALSE.
         DELP2(2) = .FALSE.
         SCD(2)   = 1.0D0

C        Y3 intermediate calculation:
C        JTYPZ=1: L(ia,J), stored permanently from 1st iter. of ground state
C                 calculation (MP2).
C        JTYPY=4: Points to Y3-files in CC_CYI I/O routines.
C        -------------------------------------------------------------------

         JTYPZ(1) = 1
         ISYMZ(1) = 1
         JTYPY(1) = 4
         DO ISYCHO = 1,NSYM
            NCHOZ(ISYCHO,1) = NUMCHO(ISYCHO)
         ENDDO

C        C intermediate calculation:
C        ---------------------------

         ISYMW(1) = ISYMX

         IF (INICIM) THEN
            ISYCIM = MULD2H(ISYMR,ISYMX)
            CALL DZERO(CIM,NT1AM(ISYCIM))
         ENDIF

C        Calculate intermediates.
C        ------------------------

         ISAVP  = IPRINT
         IPRINT = -1

         CALL CC_CYI(JTYP1,ISYMP1,JTYP2,ISYMP2,NCHP12,FAC1,NUMP12,
     &               DUMMY,IDUMMY,DUMMY,IDUMMY,FAC2,NUMUV,
     &               SCD,SCDG,IOPTDN,FOCKD,FREQ,IOPTCE,
     &               JTYPZ,ISYMZ,NCHOZ,NUMZ,JTYPY,
     &               X1AM,ISYMW,NUMW,CIM,
     &               WORK,LWORK,R2NRM,R2CNM,DELP1,DELP2)

         IPRINT = ISAVP

c-tbp:
c     write(lupri,*)
c     write(lupri,*) SECNAM,': Packed amplitude norm:'
c     write(lupri,*) 'Sq-norm: ',R2NRM
c     write(lupri,*) '2--norm: ',DSQRT(R2NRM)
c     write(lupri,*)

      ELSE

C        Error.
C        ------

         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   SECNAM,': unknown LISTR: ',LISTR
         CALL QUIT('List '//LISTR(1:3)//' not implemented in '//SECNAM)

      ENDIF

      RETURN
      END
C  /* Deck cc_cyiftrl */
      SUBROUTINE CC_CYIFTRL(LISTX,FOCKD,FIA,X1AM,WORK,LWORK,
     &                      ISYMX,ISYMR,FREQ)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Calculate Y intermediates from the "left-hand vectors"
C              specified through LISTV for F matrix transformation.
C
C     Formulas:
C
C        Y1(ai,J) = sum(bj) X(ai,bj) * Lbar(bj,J)
C
C        Y2(ai,J) = sum(bj) X(ai,bj) * L(bj,J)
C
C     where the "left-hand vectors" are identified through LISTX.
C     ISYMR is the symmetry of the "right-hand vector" and hence
C     of Lbar; ISYMX is the symmetry of the "left-hand vector".
C
C     Implemented lists:
C
C        LISTX(1:2) = 'L0' -- 0th order multipliers.
C                             Y2 is not calculated if DSKFY2 = .FALSE.
C                             [flag stored in chocc2.h].
C                             FREQ is ignored (zero is used).
C
C
C     NOTE: all Cholesky vectors needed must be available on disk.
C           No checking is done here.
C
#include <implicit.h>
      CHARACTER*(*) LISTX
      DIMENSION FOCKD(*), WORK(LWORK)
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>
#include <chocc2.h>
#include <dummy.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CC_CYIFTRL')

      PARAMETER (MXTYP = 1)
      DIMENSION SCD(MXTYP)
      INTEGER JTYP1(MXTYP), JTYP2(MXTYP), ISYMP1(MXTYP), ISYMP2(MXTYP)
      INTEGER NCHP12(8,MXTYP)
      LOGICAL DELP1(MXTYP), DELP2(MXTYP)

      PARAMETER (MXTPZ = 2)
      INTEGER JTYPZ(MXTPZ), ISYMZ(MXTPZ), JTYPY(MXTPZ)
      INTEGER NCHOZ(8,MXTPZ)

      PARAMETER (MXNMUV = 1)
      INTEGER ISYMU(MXNMUV), ISYMV(MXNMUV)

      IF (LISTX(1:2) .EQ. 'L0') THEN

C        0th order multipliers.
C        ======================

         IF (ISYMX .NE. 1) THEN
            CALL QUIT('Symmetry mismatch for list '//LISTX(1:2)//' in '
     &                //SECNAM)
         ENDIF

C        Scaling factors.
C        ----------------

         FAC1 = 1.0D0
         FAC2 = 1.0D0
         SCDG = 1.0D0

C        Incorporate denominators and set up 2CME (inherent to 'L0' vector).
C        -------------------------------------------------------------------

         IOPTDN = 1
         IOPTCE = 1

C        One Cholesky vector contraction needed (NUMP12).
C        One "extra-term" (NUMUV).
C        Two set of Y intermediates (NUMZ), one if DSKFY2 is set.
C        No additional intermediates (NUMW).
C        --------------------------------------------------------

         NUMP12 = 1
         NUMUV  = 1
         IF (DSKFY2) THEN
            NUMZ = 1
         ELSE
            NUMZ = 2
         ENDIF
         NUMW = 0

C        Cholesky vector contraction:
C        JTYP=1: L(ia,J).
C        JTYP=8: [L(ia,J) + Ltilde(ia,J)]
C                stored in the same files as used in LH transformation!
C        Keep both files.
C        --------------------------------------------------------------

         JTYP1(1)  = 1
         JTYP2(1)  = 8
         ISYMP1(1) = 1
         ISYMP2(1) = 1
         DO ISYCHO = 1,NSYM
            NCHP12(ISYCHO,1) = NUMCHO(ISYCHO)
         ENDDO
         DELP1(1) = .FALSE.
         DELP2(1) = .FALSE.
         SCD(1)   = 1.0D0

C        Additional term in amplitude construction: X1AM(ai) * FIA(bj).
C        --------------------------------------------------------------

         ISYMU(1) = 1
         ISYMV(1) = 1

C        Y1 intermediate calculation:
C        JTYPZ=9: Lbar(ai,J).
C        JTYPY=2: Points to Y1-files in CC_CYI I/O routines.
C        -------------------------------------------------------------------

         JTYPZ(1) = 9
         ISYMZ(1) = ISYMR
         JTYPY(1) = 2
         DO ISYCHO = 1,NSYM
            NCHOZ(ISYCHO,1) = NUMCHO(ISYCHO)
         ENDDO

C        Y2 intermediate calculation (if any).
C        JTYPZ=3: L(ai,J)
C        JTYPY=33: Points to Y2-files in CC_CYI I/O routines.
C        ----------------------------------------------------

         IF (.NOT. DSKFY2) THEN
            JTYPZ(2) = 3
            ISYMZ(2) = 1
            JTYPY(2) = 33
            DO ISYCHO = 1,NSYM
               NCHOZ(ISYCHO,2) = NUMCHO(ISYCHO)
            ENDDO
         ENDIF

C        Calculate intermediates.
C        Frequency is zero for 'L0'; disregard input FREQ.
C        -------------------------------------------------

         ISAVP  = IPRINT
         IPRINT = -1

         FRQL0 = 0.0D0
         CALL CC_CYI(JTYP1,ISYMP1,JTYP2,ISYMP2,NCHP12,FAC1,NUMP12,
     &               X1AM,ISYMU,FIA,ISYMV,FAC2,NUMUV,
     &               SCD,SCDG,IOPTDN,FOCKD,FRQL0,IOPTCE,
     &               JTYPZ,ISYMZ,NCHOZ,NUMZ,JTYPY,
     &               DUMMY,IDUMMY,NUMW,DUMMY,
     &               WORK,LWORK,X2NRM,X2CNM,DELP1,DELP2)

         IPRINT = ISAVP

C        Set DSKFY2 flag.
C        ----------------

         DSKFY2 = .TRUE.

      ELSE

C        Error.
C        ------

         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   SECNAM,': unknown LISTX: ',LISTX
         CALL QUIT('List '//LISTX(1:3)//' not implemented in '//SECNAM)

      ENDIF

      RETURN
      END
C  /* Deck cc_choftri */
      SUBROUTINE CC_CHOFTRI(SIGMA1,CIM,WORK,LWORK,ISYCIM,ISYVEC,ITYVEC,
     &                      NUMVEC)
C
C     Thomas Bondo Pedersen, April 2003.
C
C     Purpose: Calculate the I-term directly in MO basis,
C
C              SIGMA1(ai) = SIGMA1(ai)
C                         + 2 * sum(J)  L(ai,J) * sum(bj) L(bj,J) * CIM(bj)
C                         -     sum(Jj) L(aj,J) * sum(b)  L(bi,J) * CIM(bj)
C
C              The symmetry and type of the Cholesky vectors are ISYVEC and
C              ITYVEC, respectively, and the number of vectors is NUMVEC.
C
#include <implicit.h>
      DIMENSION SIGMA1(*), CIM(*), WORK(LWORK)
      INTEGER   NUMVEC(8)
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CC_CHOFTRI')

      PARAMETER (IOPEN = -1, IKEEP = 1)
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)

C     Start Cholesky symmetry loop.
C     -----------------------------

      DO ISYCHO = 1,NSYM

C        If nothing to do, skip this symmetry.
C        -------------------------------------

         ISYMAI = MULD2H(ISYCHO,ISYVEC)

         IF (NUMVEC(ISYCHO) .LE. 0) GO TO 1000
         IF (NT1AM(ISYMAI)  .LE. 0) GO TO 1000

C        Allocation.
C        -----------

         ISYMIJ = MULD2H(ISYMAI,ISYCIM)
         LSCR   = NRHF(1)*NRHF(ISYMIJ)
         DO ISYMI = 2,NSYM
            ISYMJ = MULD2H(ISYMI,ISYMIJ)
            LSCR  = MAX(LSCR,NRHF(ISYMI)*NRHF(ISYMJ))
         ENDDO

         KSCR = 1
         KEND = KSCR  + LSCR
         LWRK = LWORK - KEND + 1

         IF (LWRK .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need (more than): ',KEND-1,
     &      'Available       : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Open file.
C        ----------

         CALL CHO_MOP(IOPEN,ITYVEC,ISYCHO,LUVEC,1,ISYVEC)

C        Set up batch.
C        -------------

         MINMEM = NT1AM(ISYMAI)
         IF (ISYMAI .EQ. ISYCIM) MINMEM = MINMEM + 1
         NVEC = MIN(LWRK/MINMEM,NUMVEC(ISYCHO))
         IF (NVEC .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory for batch in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,I10,/,5X,A,I10,/)')
     &      'Minimum memory needed     : ',MINMEM,
     &      'Memory available for batch: ',LWRK,
     &      'Memory available in total : ',LWORK,
     &      'Cholesky symmetry         : ',ISYCHO
            CALL QUIT('Insufficient memory for batch in '//SECNAM)
         ENDIF
         NBATCH = (NUMVEC(ISYCHO) - 1)/NVEC + 1

C        Start batch loop.
C        -----------------

         DO IBATCH = 1,NBATCH

            NUMV = NVEC
            IF (IBATCH .EQ. NBATCH) THEN
               NUMV = NUMVEC(ISYCHO) - NVEC*(NBATCH - 1)
            ENDIF
            JVEC1 = NVEC*(IBATCH - 1) + 1

C           Finish allocation.
C           ------------------

            KCHOL = KEND
            KUVEC = KCHOL + NT1AM(ISYMAI)*NUMV

C           Read vectors.
C           -------------

            CALL CHO_MOREAD(WORK(KCHOL),NT1AM(ISYMAI),NUMV,JVEC1,LUVEC)

C           Calculate Coulomb part.           
C           -----------------------

            IF (ISYMAI .EQ. ISYCIM) THEN
               CALL DGEMV('T',NT1AM(ISYMAI),NUMV,
     &                    ONE,WORK(KCHOL),NT1AM(ISYMAI),CIM,1,
     &                    ZERO,WORK(KUVEC),1)
               CALL DGEMV('N',NT1AM(ISYMAI),NUMV,
     &                    TWO,WORK(KCHOL),NT1AM(ISYMAI),WORK(KUVEC),1,
     &                    ONE,SIGMA1,1)
            ENDIF

C           Calculate exchange part.
C           ------------------------

            DO IVEC = 1,NUMV
               DO ISYMI = 1,NSYM

                  ISYMB = MULD2H(ISYMI,ISYMAI)
                  ISYMJ = MULD2H(ISYMB,ISYCIM)
                  ISYMA = MULD2H(ISYMJ,ISYMAI)

                  NTOTA = MAX(NVIR(ISYMA),1)
                  NTOTB = MAX(NVIR(ISYMB),1)
                  NTOTI = MAX(NRHF(ISYMI),1)

                  KOFF1 = KCHOL + NT1AM(ISYMAI)*(IVEC - 1)
     &                  + IT1AM(ISYMB,ISYMI)
                  KOFF2 = IT1AM(ISYMB,ISYMJ) + 1
                  CALL DGEMM('T','N',
     &                       NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMB),
     &                       ONE,WORK(KOFF1),NTOTB,CIM(KOFF2),NTOTB,
     &                       ZERO,WORK(KSCR),NTOTI)

                  KOFF1 = KCHOL + NT1AM(ISYMAI)*(IVEC - 1) 
     &                  + IT1AM(ISYMA,ISYMJ)
                  KOFF2 = IT1AM(ISYMA,ISYMI) + 1
                  CALL DGEMM('N','T',
     &                       NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMJ),
     &                       XMONE,WORK(KOFF1),NTOTA,WORK(KSCR),NTOTI,
     &                       ONE,SIGMA1(KOFF2),NTOTA)

               ENDDO
            ENDDO

         ENDDO

C        Close file.
C        -----------

         CALL CHO_MOP(IKEEP,ITYVEC,ISYCHO,LUVEC,1,ISYVEC)

C        Escape point for empty symmetry.
C        --------------------------------

 1000    CONTINUE

      ENDDO

      RETURN
      END
C  /* Deck cc_choftrgjh */
      SUBROUTINE CC_CHOFTRGJH(XLAMDP,XLAMDH,XLAMD2,SIGMA1,X1AM,
     &                        WORK,LWORK,ISYML2,ISYMX)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Calculate the GJ and H terms for F-matrix transformations.
C
C     Formula:
C
C     SIGMA1(ai) = SIGMA1(ai)
C                - sum(Jj) Y1(aj,J) * L(ij,J)
C                - sum(Jj) Y2(aj,J) * Lbar(ij,J)
C                + sum(g) LamH(g,a)
C                 * sum(Jd) L(gd,J)
C                   * [ sum(b) LamP(d,b) * {Y1(bi,J) - sum(j) X(bj)*Lbar(ij,J)}
C                     + sum(b) Lam2(d,b) * {Y2(bi,J) - sum(j) X(bj)*L(ij,J)}
C                     + 2 LamP(d,i) * sum(bj) Lbar(bj,J)*X(bj)
C                     ]
C
C     where g,d refer to AO basis, and L(gd,J) = L(dg,J) is used.
C
C     Note: X1AM is the "left-hand" vector (singles part), and XLAMD2 is
C           the perturbed Lambda matrix (Lambda-bar) constructed from the
C           "right-hand" vector (singles part). Thus, ISYML2 is interpreted
C           as the symmetry of the "right-hand" vector, ISYMTR.
C
#include <implicit.h>
      DIMENSION XLAMDP(*),XLAMDH(*),XLAMD2(*),SIGMA1(*),X1AM(*)
      DIMENSION WORK(LWORK)
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      INTEGER IOFFL(8), IOFFZ(8), IOFY1(8), IOFY2(8)

      CHARACTER*12 SECNAM
      PARAMETER (SECNAM = 'CC_CHOFTRGJH')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      PARAMETER (IOPTR = 2)
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (ITYPY1 = 2, ITYPY2 = 33)

c...debug:
c     logical doy1, doy2, docbar
c     parameter (doy1 = .true., doy2 = .true., docbar = .true.)
c     write(lupri,*) SECNAM,': doy1   = ',doy1
c     write(lupri,*) SECNAM,': doy2   = ',doy2
c     write(lupri,*) SECNAM,': docbar = ',docbar

C     Set symmetries.
C     ---------------

      ISYMTR = ISYML2
      ISYCIM = MULD2H(ISYMTR,ISYMX)
      ISYMY1 = ISYCIM
      ISYMY2 = ISYMX

      IF (LOCDBG) THEN
         S1NRM = DSQRT(DDOT(NT1AM(ISYMY1),SIGMA1,1,SIGMA1,1))
         X1NRM = DSQRT(DDOT(NT1AM(ISYMX),X1AM,1,X1AM,1))
         WRITE(LUPRI,*) 'DEBUG INFO FROM ',SECNAM,':'
         WRITE(LUPRI,*) '   Symmetry of left  vector : ',ISYMX
         WRITE(LUPRI,*) '   Symmetry of right vector : ',ISYMTR
         WRITE(LUPRI,*) '   Symmetry of Y2 interm.   : ',ISYMY2
         WRITE(LUPRI,*) '   Symmetry of Y1 interm.   : ',ISYMY1
         WRITE(LUPRI,*) '   Norm of result vector    : ',S1NRM
         WRITE(LUPRI,*) '   Norm of left singles vec.: ',X1NRM
      ENDIF

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

C     Allocation.
C     -----------

      KSIGMA = KEND0
      KEND1  = KSIGMA + NT1AO(ISYMY1)
      LWRK1  = LWORK  - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   'Insufficient memory in ',SECNAM
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need (more than): ',KEND1-1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Initialize AO GJ term.
C     ----------------------

      CALL DZERO(WORK(KSIGMA),NT1AO(ISYMY1))

C     Start Cholesky symmetry loop.
C     -----------------------------

      DO ISYCHO = 1,NSYM

         IF (NUMCHO(ISYCHO) .LE. 0) GO TO 1000

         ISYBI1 = MULD2H(ISYCHO,ISYMY1)
         ISYBI2 = MULD2H(ISYCHO,ISYMY2)
         ISYAJ1 = ISYBI1
         ISYAJ2 = ISYBI2
         ISYMDI = ISYBI1
         ISYIJ1 = ISYCHO
         ISYIJ2 = MULD2H(ISYCHO,ISYMTR)
         ISYMKI = MULD2H(ISYCHO,ISYCIM)
         ISYBJ1 = ISYCHO
         ISYBJ2 = ISYIJ2

C        Allocation for one AO Cholesky vector.
C        --------------------------------------

         LREAD = MEMRD(1,ISYCHO,IOPTR)
         MAXKI = 0
         DO ISYMI = 1,NSYM
            ISYMK = MULD2H(ISYMI,ISYMKI)
            NTEST = NRHF(ISYMK)*NRHF(ISYMI)
            MAXKI = MAX(MAXKI,NTEST)
         ENDDO
         LSCR = MAX(LREAD,MAXKI)

         KCHOL = KEND1
         KREAD = KCHOL + N2BST(ISYCHO)
         KEND2 = KREAD + LSCR
         LWRK2 = LWORK - KEND2 + 1

         IF (LWRK2 .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need (more than): ',KEND2-1,
     &      'Available       : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Set up batch.
C        SCRL: L(ij,#J),Lbar(ij,#J),Y1(b#J,i),Y2(b#J,i), <-- simultaneously
C              Lbar(bj,#J),
C        SCRZ: Y1(bi,#J),
C              Y2(bi,#J),
C              Z(d#J,i).
C        ------------------------------------------------------------------

         LENS = NMATIJ(ISYIJ1) + NMATIJ(ISYIJ2)
     &        + NT1AM(ISYBI1)  + NT1AM(ISYBI2)

         LENL = MAX(LENS,NT1AM(ISYBJ2),N2BST(ISYCHO))
         LENZ = MAX(NT1AM(ISYBI1),NT1AM(ISYBI2),NT1AO(ISYMDI))

         MINMEM = LENL + LENZ
         IF (ISYBJ1 .EQ. ISYCIM) MINMEM = MINMEM + 1
         IF (MINMEM .GT. 0) THEN
            NVEC = MIN(LWRK2/MINMEM,NUMCHO(ISYCHO))
         ELSE IF (MINMEM .EQ. 0) THEN
            GO TO 1000
         ELSE
            WRITE(LUPRI,'(//,5X,A,A,A,I10,/)')
     &      'MINMEM is negative in ',SECNAM,': ',MINMEM
            CALL QUIT('Error in '//SECNAM)
         ENDIF

         IF (NVEC .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory for Cholesky batch in ',SECNAM
            WRITE(LUPRI,'(5X,A,I2,A)')
     &      '(Cholesky symmetry:',ISYCHO,')'
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,I10,/)')
     &      'Minimum memory required for batch: ',MINMEM,
     &      'Memory available for batch       : ',LWRK2,
     &      'Total memory available           : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

         NBATCH = (NUMCHO(ISYCHO) - 1)/NVEC + 1

C        Open Y intermediate files.
C        --------------------------

         CALL CC_CYIOP(-1,ISYCHO,ITYPY1)
         CALL CC_CYIOP(-1,ISYCHO,ITYPY2)

C        Open MO Cholesky vector files.
C        ------------------------------

         IF (ISYBJ2 .EQ. ISYMX) CALL CHO_MOP(-1,9,ISYCHO,LUBJ2,1,ISYMTR)
         CALL CHO_MOP(-1,4,ISYCHO,LUIJ1,1,1)
         CALL CHO_MOP(-1,10,ISYCHO,LUIJ2,1,ISYMTR)

C        Start batch loop.
C        -----------------

         DO IBATCH = 1,NBATCH

            NUMV = NVEC
            IF (IBATCH .EQ. NBATCH) THEN
               NUMV = NUMCHO(ISYCHO) - NVEC*(NBATCH - 1)
            ENDIF
            JVEC1 = NVEC*(IBATCH - 1) + 1

C           Set up local index arrays.
C           --------------------------

            ICOUNL = 0
            ICOUNZ = 0
            DO ISYMD = 1,NSYM
               ISYMI = MULD2H(ISYMD,ISYMDI)
               ISYMG = MULD2H(ISYMD,ISYCHO)
               IOFFL(ISYMD) = ICOUNL
               IOFFZ(ISYMD) = ICOUNZ
               LENDJ  = NBAS(ISYMD)*NUMV
               ICOUNL = ICOUNL + NBAS(ISYMG)*LENDJ
               ICOUNZ = ICOUNZ + LENDJ*NRHF(ISYMI)
            ENDDO

            ICOUY1 = 0
            ICOUY2 = 0
            DO ISYMB = 1,NSYM
               ISYMI1 = MULD2H(ISYMB,ISYBI1)
               ISYMI2 = MULD2H(ISYMB,ISYBI2)
               IOFY1(ISYMB) = ICOUY1
               IOFY2(ISYMB) = ICOUY2
               LENBJ  = NVIR(ISYMB)*NUMV
               ICOUY1 = ICOUY1 + LENBJ*NRHF(ISYMI1)
               ICOUY2 = ICOUY2 + LENBJ*NRHF(ISYMI2)
            ENDDO

C           Complete allocation.
C           --------------------

            KSCRL = KEND2
            KSCRZ = KSCRL + LENL*NUMV
            KUVEC = KSCRZ + LENZ*NUMV

            KIJ1 = KSCRL
            KIJ2 = KIJ1 + NMATIJ(ISYIJ1)*NUMV
            KY1  = KIJ2 + NMATIJ(ISYIJ2)*NUMV
            KY2  = KY1  + NT1AM(ISYBI1)*NUMV

C           Read Y1(bi,#J) into SCRZ.
C           -------------------------

c     if (doy1) then

            CALL CC_CYIRDF(WORK(KSCRZ),NUMV,JVEC1,ISYCHO,ISYMY1,ITYPY1)

c     endif

C           Read L(ij,#J)    into KIJ1 (in SCRL).
C           Read Lbar(ij,#J) into KIJ2 (in SCRL).
C           -------------------------------------

            CALL CHO_MOREAD(WORK(KIJ1),NMATIJ(ISYIJ1),NUMV,JVEC1,LUIJ1)
            CALL CHO_MOREAD(WORK(KIJ2),NMATIJ(ISYIJ2),NUMV,JVEC1,LUIJ2)

c     if (doy1) then

C           Calculate H term 1.
C           -------------------

            DO IVEC = 1,NUMV
               DO ISYMJ = 1,NSYM

                  ISYMI = MULD2H(ISYMJ,ISYIJ1)
                  ISYMA = MULD2H(ISYMJ,ISYAJ1)

                  KOFF1 = KSCRZ + NT1AM(ISYAJ1)*(IVEC - 1)
     &                  + IT1AM(ISYMA,ISYMJ)
                  KOFF2 = KIJ1 + NMATIJ(ISYIJ1)*(IVEC - 1)
     &                  + IMATIJ(ISYMI,ISYMJ)
                  KOFF3 = IT1AM(ISYMA,ISYMI) + 1

                  NTOTA = MAX(NVIR(ISYMA),1)
                  NTOTI = MAX(NRHF(ISYMI),1)

                  CALL DGEMM('N','T',
     &                       NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMJ),
     &                       XMONE,WORK(KOFF1),NTOTA,WORK(KOFF2),NTOTI,
     &                       ONE,SIGMA1(KOFF3),NTOTA)

               ENDDO
            ENDDO

C           Subtract sum(j) X(bj) * Lbar(ij,#J) from Y1(bi,#J).
C           ---------------------------------------------------

            DO IVEC = 1,NUMV
               DO ISYMJ = 1,NSYM

                  ISYMI = MULD2H(ISYMJ,ISYIJ2)
                  ISYMB = MULD2H(ISYMJ,ISYMX)

                  KOFF1 = IT1AM(ISYMB,ISYMJ) + 1
                  KOFF2 = KIJ2 + NMATIJ(ISYIJ2)*(IVEC - 1)
     &                  + IMATIJ(ISYMI,ISYMJ)
                  KOFF3 = KSCRZ + NT1AM(ISYBI1)*(IVEC - 1)
     &                  + IT1AM(ISYMB,ISYMI)

                  NTOTB = MAX(NVIR(ISYMB),1)
                  NTOTI = MAX(NRHF(ISYMI),1)

                  CALL DGEMM('N','T',
     &                       NVIR(ISYMB),NRHF(ISYMI),NRHF(ISYMJ),
     &                       XMONE,X1AM(KOFF1),NTOTB,WORK(KOFF2),NTOTI,
     &                       ONE,WORK(KOFF3),NTOTB)

               ENDDO
            ENDDO

C           Resort Y1 intermediate as Y1(b#J,i), store in KY1 (in SCRL).
C           ------------------------------------------------------------

            DO IVEC = 1,NUMV
               DO ISYMI = 1,NSYM

                  ISYMB = MULD2H(ISYMI,ISYBI1)
                  LENBJ = NVIR(ISYMB)*NUMV

                  DO I = 1,NRHF(ISYMI)

                     KOFF1 = KSCRZ + NT1AM(ISYBI1)*(IVEC - 1)
     &                     + IT1AM(ISYMB,ISYMI) + NVIR(ISYMB)*(I - 1)
                     KOFF2 = KY1 + IOFY1(ISYMB) + LENBJ*(I - 1)
     &                     + NVIR(ISYMB)*(IVEC - 1)

                     CALL DCOPY(NVIR(ISYMB),WORK(KOFF1),1,WORK(KOFF2),1)

                  ENDDO

               ENDDO
            ENDDO

c     endif

C           Read Y2(bi,#J) into SCRZ.
C           -------------------------

c     if (doy2) then

            CALL CC_CYIRDF(WORK(KSCRZ),NUMV,JVEC1,ISYCHO,ISYMY2,ITYPY2)

C           Calculate H term 2.
C           -------------------

            DO IVEC = 1,NUMV
               DO ISYMJ = 1,NSYM

                  ISYMI = MULD2H(ISYMJ,ISYIJ2)
                  ISYMA = MULD2H(ISYMJ,ISYAJ2)

                  KOFF1 = KSCRZ + NT1AM(ISYAJ2)*(IVEC - 1)
     &                  + IT1AM(ISYMA,ISYMJ)
                  KOFF2 = KIJ2 + NMATIJ(ISYIJ2)*(IVEC - 1)
     &                  + IMATIJ(ISYMI,ISYMJ)
                  KOFF3 = IT1AM(ISYMA,ISYMI) + 1

                  NTOTA = MAX(NVIR(ISYMA),1)
                  NTOTI = MAX(NRHF(ISYMI),1)

                  CALL DGEMM('N','T',
     &                       NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMJ),
     &                       XMONE,WORK(KOFF1),NTOTA,WORK(KOFF2),NTOTI,
     &                       ONE,SIGMA1(KOFF3),NTOTA)

               ENDDO
            ENDDO

C           Subtract sum(j) X(bj) * L(ij,#J) from Y2(bi,#J).
C           ------------------------------------------------

            DO IVEC = 1,NUMV
               DO ISYMJ = 1,NSYM

                  ISYMI = MULD2H(ISYMJ,ISYIJ1)
                  ISYMB = MULD2H(ISYMJ,ISYMX)

                  KOFF1 = IT1AM(ISYMB,ISYMJ) + 1
                  KOFF2 = KIJ1 + NMATIJ(ISYIJ1)*(IVEC - 1)
     &                  + IMATIJ(ISYMI,ISYMJ)
                  KOFF3 = KSCRZ + NT1AM(ISYBI2)*(IVEC - 1)
     &                  + IT1AM(ISYMB,ISYMI)

                  NTOTB = MAX(NVIR(ISYMB),1)
                  NTOTI = MAX(NRHF(ISYMI),1)

                  CALL DGEMM('N','T',
     &                       NVIR(ISYMB),NRHF(ISYMI),NRHF(ISYMJ),
     &                       XMONE,X1AM(KOFF1),NTOTB,WORK(KOFF2),NTOTI,
     &                       ONE,WORK(KOFF3),NTOTB)

               ENDDO
            ENDDO

C           Resort Y2 intermediate as Y2(b#J,i), store in KY2 (in SCRL).
C           ------------------------------------------------------------

            DO IVEC = 1,NUMV
               DO ISYMI = 1,NSYM

                  ISYMB = MULD2H(ISYMI,ISYBI2)
                  LENBJ = NVIR(ISYMB)*NUMV

                  DO I = 1,NRHF(ISYMI)

                     KOFF1 = KSCRZ + NT1AM(ISYBI2)*(IVEC - 1)
     &                     + IT1AM(ISYMB,ISYMI) + NVIR(ISYMB)*(I - 1)
                     KOFF2 = KY2 + IOFY2(ISYMB) + LENBJ*(I - 1)
     &                     + NVIR(ISYMB)*(IVEC - 1)

                     CALL DCOPY(NVIR(ISYMB),WORK(KOFF1),1,WORK(KOFF2),1)

                  ENDDO

               ENDDO
            ENDDO

c     endif

C           Backtransform the Y intermediates to get Z(d#J,i).
C           First Y1, then Y2.
C           --------------------------------------------------

c     if (doy1) then

            DO ISYMB = 1,NSYM

               ISYMI = MULD2H(ISYMB,ISYBI1)
               ISYMD = ISYMB

               KOFF1 = IGLMVI(ISYMD,ISYMB) + 1
               KOFF2 = KY1   + IOFY1(ISYMB)
               KOFF3 = KSCRZ + IOFFZ(ISYMD)

               NTOTD = MAX(NBAS(ISYMD),1)
               NTOTB = MAX(NVIR(ISYMB),1)

               CALL DGEMM('N','N',
     &                    NBAS(ISYMD),NUMV*NRHF(ISYMI),NVIR(ISYMB),
     &                    ONE,XLAMDP(KOFF1),NTOTD,WORK(KOFF2),NTOTB,
     &                    ZERO,WORK(KOFF3),NTOTD)

            ENDDO

c     else
c        call dzero(work(kscrz),nt1ao(isymdi)*numv)
c     endif

c     if (doy2) then

            DO ISYMB = 1,NSYM

               ISYMI = MULD2H(ISYMB,ISYBI2)
               ISYMD = MULD2H(ISYMB,ISYML2)

               KOFF1 = IGLMVI(ISYMD,ISYMB) + 1
               KOFF2 = KY2   + IOFY2(ISYMB)
               KOFF3 = KSCRZ + IOFFZ(ISYMD)

               NTOTD = MAX(NBAS(ISYMD),1)
               NTOTB = MAX(NVIR(ISYMB),1)

               CALL DGEMM('N','N',
     &                    NBAS(ISYMD),NUMV*NRHF(ISYMI),NVIR(ISYMB),
     &                    ONE,XLAMD2(KOFF1),NTOTD,WORK(KOFF2),NTOTB,
     &                    ONE,WORK(KOFF3),NTOTD)

            ENDDO

c     endif

c     if (docbar) then

            IF (ISYBJ2 .EQ. ISYMX) THEN

C              Read Lbar(bj,#J) into SCRL.
C              ---------------------------

               CALL CHO_MOREAD(WORK(KSCRL),NT1AM(ISYBJ2),NUMV,JVEC1,
     &                         LUBJ2)

C              Calculate U(#J) = 2 * sum(bj) Lbar(bj,#J) * X(bj).
C              --------------------------------------------------

               NTOBJ = MAX(NT1AM(ISYBJ2),1)
               CALL DGEMV('T',NT1AM(ISYBJ2),NUMV,
     &                    ONE,WORK(KSCRL),NTOBJ,X1AM,1,
     &                    ZERO,WORK(KUVEC),1)
               CALL DSCAL(NUMV,TWO,WORK(KUVEC),1)

C              Set up Z contribution:
C              Z(d#J,i) = Z(d#J,i) + LamP(d,i) * U(#J)
C              ---------------------------------------

               DO ISYMD = 1,NSYM
                  ISYMI = ISYMD    ! true, since ISYBJ1 .EQ. ISYCIM
                  DO I = 1,NRHF(ISYMI)
                     KOFF1 = IGLMRH(ISYMD,ISYMI)
     &                     + NBAS(ISYMD)*(I - 1) + 1
                     DO IVEC = 1,NUMV
                        FAC = WORK(KUVEC+IVEC-1)
                        KOFF2 = KSCRZ + IOFFZ(ISYMD)
     &                        + NBAS(ISYMD)*NUMV*(I - 1)
     &                        + NBAS(ISYMD)*(IVEC - 1)
                        CALL DAXPY(NBAS(ISYMD),FAC,XLAMDP(KOFF1),1,
     &                                             WORK(KOFF2),1)
                     ENDDO
                  ENDDO
               ENDDO

            ENDIF

c     else
c        if (isybj2 .eq. isymx) then
c           call dzero(work(kuvec),numv)
c        endif
c     endif

C           Read AO Cholesky vectors and store as L(g,d#J) in SCRL.
C           -------------------------------------------------------

            DO IVEC = 1,NUMV

               JVEC = JVEC1 + IVEC - 1
               CALL CHO_READN(WORK(KCHOL),JVEC,1,WORK(KIND1),IDUM2,
     &                        ISYCHO,IOPTR,WORK(KREAD),LREAD)

               DO ISYMD = 1,NSYM

                  ISYMG = MULD2H(ISYMD,ISYCHO)

                  LENGD = NBAS(ISYMG)*NBAS(ISYMD)

                  KOFF1 = KCHOL + IAODIS(ISYMG,ISYMD)
                  KOFF2 = KSCRL + IOFFL(ISYMD) + LENGD*(IVEC - 1)

                  CALL DCOPY(LENGD,WORK(KOFF1),1,WORK(KOFF2),1)

               ENDDO

            ENDDO

C           Calculate GJ term in AO basis.
C           ------------------------------

            DO ISYMD = 1,NSYM

               ISYMG = MULD2H(ISYMD,ISYCHO)
               ISYMI = MULD2H(ISYMD,ISYMDI)

               KOFF1 = KSCRL  + IOFFL(ISYMD)
               KOFF2 = KSCRZ  + IOFFZ(ISYMD)
               KOFF3 = KSIGMA + IT1AO(ISYMG,ISYMI)

               LENDJ = NBAS(ISYMD)*NUMV

               NTOTG = MAX(NBAS(ISYMG),1)
               NTODJ = MAX(LENDJ,1)

               CALL DGEMM('N','N',NBAS(ISYMG),NRHF(ISYMI),LENDJ,
     &                    ONE,WORK(KOFF1),NTOTG,WORK(KOFF2),NTODJ,
     &                    ONE,WORK(KOFF3),NTOTG)

            ENDDO

         ENDDO

C        Close Y intermediate files.
C        ---------------------------

         CALL CC_CYIOP(1,ISYCHO,ITYPY1)
         CALL CC_CYIOP(1,ISYCHO,ITYPY2)

C        Close MO Cholesky vector files.
C        -------------------------------

         IF (ISYBJ1.EQ.ISYCIM) CALL CHO_MOP(1,9,ISYCHO,LUBJ2,1,ISYMTR)
         CALL CHO_MOP(1,4,ISYCHO,LUIJ1,1,1)
         CALL CHO_MOP(1,10,ISYCHO,LUIJ2,1,ISYMTR)

C        Escape point for empty symmetry.
C        --------------------------------

 1000    CONTINUE

      ENDDO

C     Transform the GJ term to MO basis.
C     ----------------------------------

      DO ISYMI = 1,NSYM

         ISYMA = MULD2H(ISYMI,ISYMY1)
         ISYMG = ISYMA

         NTOTG = MAX(NBAS(ISYMG),1)
         NTOTA = MAX(NVIR(ISYMA),1)

         KOFF1 = IGLMVI(ISYMG,ISYMA) + 1
         KOFF2 = KSIGMA + IT1AO(ISYMG,ISYMI)
         KOFF3 = IT1AM(ISYMA,ISYMI)  + 1

         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMG),
     &              ONE,XLAMDH(KOFF1),NTOTG,WORK(KOFF2),NTOTG,
     &              ONE,SIGMA1(KOFF3),NTOTA)

      ENDDO

      RETURN
      END
C  /* Deck cc_choftr5 */
      SUBROUTINE CC_CHOFTR5(IFTRAN,NFTRAN,LISTR,IOPTRES,FILFMA,
     &                      IFDOTS,FCONS,MXVEC,WORK,LWORK)
C
C     Thomas Bondo Pedersen, April 2003.
C     
C     Purpose: Calculate dot products of Cholesky CC2 F matrix
C              transformed vectors and vectors specified through
C              list type FILFMA with indices IFDOTS. Results are
C              added into FCONS. NOTE that IFTRAN(3,*) is needed
C              to address the F-transformed vectors on disk in contrast
C              to the conventional code, where the F-transformed vectors
C              are calculated on-the-fly.
C
C     Formula (F22=0 for CC2):
C
C        d = (F*R)*R'
C          = (F11*R1 + F12*R2)*R1' + (F21*R1)*R2'
C          = (F11*R1 + F12*R2)*R1' + (F12*R2')*R1
C
C     where primed vectors are specified through FILFMA and unprimed ones
C     through LISTR with indices stored in IFDOTS(*) and IFTRAN(3,*),
C     respectively.
C
#include <implicit.h>
      INTEGER IFTRAN(3,NFTRAN), IFDOTS(MXVEC,NFTRAN)
      CHARACTER*(*) LISTR, FILFMA
      DIMENSION FCONS(MXVEC,NFTRAN), WORK(LWORK)
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>
#include <dummy.h>
#include <chocc2.h>
#include <ccr1rsp.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CC_CHOFTR5')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      CHARACTER*3  FILRES, FLXRES
      CHARACTER*8  LABEL
      CHARACTER*10 MODEL

      INTEGER ILSTSYM

C     Build file lists for F transformed vectors.
C     -------------------------------------------

      IF ((LISTR(1:2).EQ. 'R1') .AND. (FILFMA(1:2).EQ.'R1')) THEN
         FILRES(1:3) = 'F1 '
         FLXRES(1:3) = 'XF1'
      ELSE
         WRITE(LUPRI,'(//,5X,A,A,A,A,A,A)')
     &   SECNAM,
     &   ': LISTR = ',LISTR(1:3),
     &   '; FILFMA = ',FILFMA(1:3),
     &   ' not implemented.'
         CALL QUIT('List combination not implemented in '//SECNAM)
      ENDIF

C     Set MODEL.
C     ----------

      IF (CC2) THEN
         MODEL = 'CC2       '
      ELSE
         WRITE(LUPRI,'(//,5X,A,A)')
     &   'CC2 is the only possible model in ',SECNAM
         CALL QUIT('Model not implemented in '//SECNAM)
      ENDIF

      DO ITRAN = 1,NFTRAN

C        Get R-vector/operator info.
C        ---------------------------

         ILSTNR = IFTRAN(2,ITRAN)
         ISYMR1 = ILSTSYM(LISTR,ILSTNR)
         LABEL  = LRTLBL(ILSTNR)
         FREQ   = FRQLRT(ILSTNR)

C        Allocation.
C        -----------

         KFTR1 = 1
         KEND1 = KFTR1 + NT1AM(ISYMR1)
         LWRK1 = LWORK - KEND1 + 1

         KR1AM = KFTR1

         IF (LWRK1 .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Insufficient memory in ',SECNAM,' - allocation: FR'
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need     : ',KEND1-1,
     &      'Available: ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Read singles F-transformed vector.
C        ----------------------------------

         IOPT = 1
         CALL CC_RDRSP(FILRES,IFTRAN(3,ITRAN),ISYMR1,IOPT,MODEL,
     &                 WORK(KFTR1),DUMMY)

C        Calculate requested dot-products, singles contributions.
C        --------------------------------------------------------

         IOPT = 1
         CALL CCDOTRSP(IFDOTS,FCONS,IOPT,FILFMA,ITRAN,NFTRAN,MXVEC,
     &                 WORK(KFTR1),DUMMY,ISYMR1,WORK(KEND1),LWRK1)

C        Read singles vector.
C        --------------------

         IOPT = 1
         CALL CC_RDRSP(LISTR,ILSTNR,ISYMR1,IOPT,MODEL,WORK(KR1AM),DUMMY)

C        Allocation.
C        -----------

         KXF1  = KEND1
         KEND2 = KXF1  + NT1AM(ISYMR1)
         LWRK2 = LWORK - KEND2 + 1

         IF (LWRK2 .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Insufficient memory in ',SECNAM,' - allocation: XFR'
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need     : ',KEND2-1,
     &      'Available: ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Calculate requested dot-products, doubles contributions.
C        --------------------------------------------------------

         IVEC = 1
         DO WHILE ((IFDOTS(IVEC,ITRAN).NE.0) .AND. (IVEC.LE.MXVEC))

            IXF1   = IFDOTS(IVEC,ITRAN)
            ISYVEC = ILSTSYM(FLXRES,IXF1)

            IF (ISYVEC .NE. ISYMR1) THEN
               WRITE(LUPRI,*) 'Symmetry mismatch in ',SECNAM,':'
               WRITE(LUPRI,*) 'IVEC  : ',IVEC
               WRITE(LUPRI,*) 'ITRAN : ',ITRAN
               WRITE(LUPRI,*) 'IXF1  : ',IXF1
               WRITE(LUPRI,*) 'ISYVEC: ',ISYVEC
               WRITE(LUPRI,*) 'ISYMR1: ',ISYMR1
               WRITE(LUPRI,*) 'FLXRES: ',FLXRES
               CALL QUIT('Symmetry mismatch in '//SECNAM)
            ENDIF

            IOPT = 1
            CALL CC_RDRSP(FLXRES,IXF1,ISYMR1,IOPT,MODEL,
     &                    WORK(KXF1),DUMMY)

            CON1 = FCONS(IVEC,ITRAN)
            CON2 = DDOT(NT1AM(ISYMR1),WORK(KR1AM),1,WORK(KXF1),1)
            FCONS(IVEC,ITRAN) = FCONS(IVEC,ITRAN) + CON2

            IF (LOCDBG) THEN
               WRITE (LUPRI,*) SECNAM,': IVEC,ITRAN:',IVEC,ITRAN
               WRITE (LUPRI,*) SECNAM,': DOTPROD,CON1,CON2:',
     &                         FCONS(IVEC,ITRAN),CON1,CON2
            ENDIF

c-tbp:
c     write(lupri,*) 'FCONS(',IVEC,',',ITRAN,') is:'
c     write(lupri,*) 'F * R1(',LABEL,',',FREQ,') * ',
c    &                 'R1(',LRTLBL(IXF1),',',FRQLRT(IXF1),')'

            IVEC = IVEC + 1

         ENDDO

      ENDDO

      RETURN
      END
