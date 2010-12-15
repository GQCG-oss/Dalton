C  /* Deck cc_choeim */
      SUBROUTINE CC_CHOEIM(FOCKD,T1AM,WORK,LWORK,DELY)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Calculate global E intermediates for Cholesky CC2 and
C              store them on disk.
C
C     IF (DELY): the Y intermediate files are deleted after use.
C
#include <implicit.h>
      DIMENSION FOCKD(*), T1AM(*), WORK(LWORK)
      LOGICAL DELY
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_CHOEIM')

      PARAMETER (INFO = 2)

C     Start timing.
C     -------------

      TIMT = SECOND()

C     Calculate Lambda matrices.
C     --------------------------

      TIML = SECOND()

      KLAMDP = 1
      KLAMDH = KLAMDP + NGLMDT(1)
      KEIJ   = KLAMDH + NGLMDT(1)
      KEAB   = KEIJ   + NMATIJ(1)
      KEND1  = KEAB   + NMATAB(1)
      LWRK1  = LWORK  - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   'Insufficient memory in ',SECNAM
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need (more than): ',KEND1-1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),T1AM,WORK(KEND1),LWRK1)

C     Calculate E intermediates.
C     --------------------------

      CALL CC_CHOEIM0(FOCKD,T1AM,WORK(KLAMDP),WORK(KLAMDH),
     &                WORK(KEIJ),WORK(KEAB),WORK(KEND1),LWRK1,DELY)

      TIML = SECOND() - TIML

C     Write E intermediates to disk.
C     ------------------------------

      TIMW = SECOND()
      CALL CHO_IMOP(-1,1,LUE,1)
      CALL CHO_MOWRITE(WORK(KEIJ),NMATIJ(1),1,1,LUE)
      CALL CHO_IMOP(1,1,LUE,1)
      CALL CHO_IMOP(-1,2,LUE,1)
      CALL CHO_MOWRITE(WORK(KEAB),NMATAB(1),1,1,LUE)
      CALL CHO_IMOP(1,2,LUE,1)
      TIMW = SECOND() - TIMW

      IF (LOCDBG) THEN
         EIJNRM = DSQRT(DDOT(NMATIJ(1),WORK(KEIJ),1,WORK(KEIJ),1))
         EABNRM = DSQRT(DDOT(NMATAB(1),WORK(KEAB),1,WORK(KEAB),1))
         WRITE(LUPRI,'(A,A,1P,D22.15)')
     &   SECNAM,': norm of EIJ: ',EIJNRM
         WRITE(LUPRI,'(A,A,1P,D22.15)')
     &   SECNAM,': norm of EAB: ',EABNRM
      ENDIF

C     Print.
C     ------

      IF (IPRINT .GE. INFO) THEN
         TIMT = SECOND() - TIMT
         CALL HEADER('Timing of '//SECNAM,-1)
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for calculating E intermediates: ',TIML,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for I/O of the  E intermediates: ',TIMW,' seconds'
         WRITE(LUPRI,'(5X,A)')
     &  '-------------------------------------------------------------'
         WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &   'Total time                                 ',TIMT,' seconds'
      ENDIF

      RETURN
      END
C  /* Deck cc_choeim0 */
      SUBROUTINE CC_CHOEIM0(FOCKD,T1AM,XLAMDP,XLAMDH,EIJ,EAB,WORK,LWORK,
     &                      DELY)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Calculate tot. sym. (ground state) E intermediates,
C
C              E(ab) = F(ab) - sum(ckj) [2*t(aj,ck) - t(ak,cj)] * (kc|jb)
C
C              E(ij) = F(ij) + sum(ckb) [2*t(bi,ck) - t(bk,ci)] * (kc|jb)
C
C              using Y-intermediates from disk for the T2-dependent terms
C              and
C
C              F(ij) = delta(ij) e(i)
C                    + sum(alpha,beta) LamP(alpha,i) W(alpha,beta) LamH(beta,j)
C
C              F(ab) = delta(ab) e(a)
C                    + sum(alpha,beta) LamP(alpha,a) W(alpha,beta) LamH(beta,b)
C
C              where
C
C              W(alpha,beta) =
C                sum(ck) [2(alpha beta|kc)-(alpha c|k beta)] * t(ck)
C
C              is evaluated as a generalized CC Fock matrix.
C
C     IF (DELY): the Y intermediate files are deleted after use.
C
#include <implicit.h>
      DIMENSION FOCKD(*), T1AM(*), XLAMDP(*), XLAMDH(*)
      DIMENSION EIJ(*), EAB(*), WORK(LWORK)
      LOGICAL DELY
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>

      INTEGER AA

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CC_CHOEIM0')

      PARAMETER (INFO = 3)
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0)

C     Start timing.
C     -------------

      TIMT = SECOND()

C     Set symmetries.
C     ---------------

      ISIDE = 0
      ISYMT = 1
      ISYMY = ISYMT

C     Initialize.
C     -----------

      TIMI = SECOND()

      CALL DZERO(EIJ,NMATIJ(1))
      CALL DZERO(EAB,NMATAB(1))

      DO ISYM = 1,NSYM
         DO I = 1,NRHF(ISYM)
            KI = IRHF(ISYM) + I
            II = IMATIJ(ISYM,ISYM) + NRHF(ISYM)*(I - 1) + I
            EIJ(II) = FOCKD(KI)
         ENDDO
         DO A = 1,NVIR(ISYM)
            KA = IVIR(ISYM) + A
            AA = IMATAB(ISYM,ISYM) + NVIR(ISYM)*(A - 1) + A
            EAB(AA) = FOCKD(KA)
         ENDDO
      ENDDO

      TIMI = SECOND() - TIMI

C     Calculate contributions from Y intermediates.
C     ---------------------------------------------

      TIMY = SECOND()
      CALL CC_CHOEY(EIJ,EAB,WORK,LWORK,ISYMY,ISIDE,1,1,NUMCHO,DELY)
      TIMY = SECOND() - TIMY

C     Allocation.
C     -----------

      KLAMH = 1
      KWMAT = KLAMH + NGLMRH(ISYMT)
      KEND1 = KWMAT + N2BST(ISYMT)
      LWRK1 = LWORK - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   'Insufficient memory in ',SECNAM
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need (more than): ',KEND1-1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Calculate "perturbed" occupied Lambda-hole matrix.
C     --------------------------------------------------

      TIMB = SECOND()
      DO ISYMK = 1,NSYM

         ISYMC = MULD2H(ISYMK,ISYMT)
         ISYMB = ISYMC

         KOFFH = IGLMVI(ISYMB,ISYMC) + 1
         KOFFT = IT1AM(ISYMC,ISYMK)  + 1
         KOFF3 = KLAMH + IGLMRH(ISYMB,ISYMK)

         NTOTB = MAX(NBAS(ISYMB),1)
         NTOTC = MAX(NVIR(ISYMC),1)

         CALL DGEMM('N','N',NBAS(ISYMB),NRHF(ISYMK),NVIR(ISYMC),
     &              ONE,XLAMDH(KOFFH),NTOTB,T1AM(KOFFT),NTOTC,
     &              ZERO,WORK(KOFF3),NTOTB)

      ENDDO
      TIMB = SECOND() - TIMB

C     Calculate W intermediate; essentially a "perturbed" Fock matrix:
C     W(ab) = sum(ck) [2(ab|kc)-(ac|kb)] * t(ck), where a=alpha,b=beta.
C     -----------------------------------------------------------------

      TIMF = SECOND()
      CALL DZERO(WORK(KWMAT),N2BST(ISYMT))
      CALL CC_CHOFCKP(XLAMDP,WORK(KLAMH),WORK(KWMAT),WORK(KEND1),LWRK1,
     &                1,ISYMT,1)
      TIMF = SECOND() - TIMF

C     Transform to MO basis.
C     ----------------------

      TITR = SECOND()
      CALL CC_CHOETRF(XLAMDP,XLAMDH,EIJ,EAB,WORK(KWMAT),
     &                WORK(KEND1),LWRK1,1,1,ISYMT)
      TITR = SECOND() - TITR

C     Debug: Check that we get the right F(ia) block
C            by transforming W.
C     ----------------------------------------------

      IF (LOCDBG) THEN

         KTRGT = KEND1
         KSCR1 = KTRGT + NT1AM(1)
         KEND2 = KSCR1 + NT1AM(1)
         LWRK2 = LWORK - KEND2 + 1

         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,'(/,5X,A,A,/)')
     &      SECNAM,': Insuf. mem. for Fock matrix test!'
            GO TO 999
         ENDIF

         CALL ONEL_OP(-1,3,LUFIA)
         CALL CHO_MOREAD(WORK(KTRGT),NT1AM(1),1,1,LUFIA)
         CALL ONEL_OP(1,3,LUFIA)

         DO ISYMI = 1,NSYM

            ISYMD = ISYMI
            ISYMA = ISYMI
            ISYMG = ISYMD

            KSCR2 = KEND2
            KEND3 = KSCR2 + NBAS(ISYMD)*NRHF(ISYMI)
            LWRK3 = LWORK - KEND3 + 1

            IF (LWRK3 .LT. 0) THEN
               WRITE(LUPRI,'(/,5X,A,A,/)')
     &         SECNAM,': Insuf. mem. for Fock matrix test!'
               GO TO 999
            ENDIF

            NTOTG = MAX(NBAS(ISYMG),1)
            NTOTD = MAX(NBAS(ISYMD),1)
            NTOTA = MAX(NVIR(ISYMA),1)

            KOFFW = KWMAT + IAODIS(ISYMG,ISYMD)
            KOFFP = IGLMRH(ISYMG,ISYMI) + 1
            KOFFH = IGLMVI(ISYMD,ISYMA) + 1
            KOFFF = KSCR1 + IT1AM(ISYMA,ISYMI)

            CALL DGEMM('T','N',NBAS(ISYMD),NRHF(ISYMI),NBAS(ISYMG),
     &                 ONE,WORK(KOFFW),NTOTG,XLAMDP(KOFFP),NTOTG,
     &                 ZERO,WORK(KSCR2),NTOTD)
            CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMD),
     &                 ONE,XLAMDH(KOFFH),NTOTD,WORK(KSCR2),NTOTD,
     &                 ZERO,WORK(KOFFF),NTOTA)

         ENDDO

         CALL HEADER('Checking F(ia) in '//SECNAM,-1)

         TRGNRM = DSQRT(DDOT(NT1AM(1),WORK(KTRGT),1,WORK(KTRGT),1))
         TSTNRM = DSQRT(DDOT(NT1AM(1),WORK(KSCR1),1,WORK(KSCR1),1))
         CALL DAXPY(NT1AM(1),XMONE,WORK(KTRGT),1,WORK(KSCR1),1)
         DIFNRM = DDOT(NT1AM(1),WORK(KSCR1),1,WORK(KSCR1),1)
         DIFRMS = DSQRT(DIFNRM/NT1AM(1))
         DIFNRM = DSQRT(DIFNRM)

         WRITE(LUPRI,'(5X,A,1P,D22.15)')
     &   'Target F(ia) norm : ',TRGNRM
         WRITE(LUPRI,'(5X,A,1P,D22.15)')
     &   'Local  F(ia) norm : ',TSTNRM
         WRITE(LUPRI,'(5X,A,1P,D22.15)')
     &   'Difference        : ',TRGNRM-TSTNRM
         WRITE(LUPRI,'(5X,A,1P,D22.15)')
     &   'Norm of difference: ',DIFNRM
         WRITE(LUPRI,'(5X,A,1P,D22.15,/)')
     &   'RMS error         : ',DIFRMS

c         CALL HEADER('Difference F(ia) Matrix',-1)
c         CALL NOCC_PRT(WORK(KSCR1),1,'AI  ')

  999    CONTINUE

      ENDIF

C     Print.
C     ------

      IF (IPRINT .GT. INFO) THEN
         TIMT = SECOND() - TIMT
         CALL AROUND('Output from '//SECNAM)
         WRITE(LUPRI,'(10X,A,/)')
     &   '- E intermediates generated.'
         WRITE(LUPRI,'(10X,A,F10.2,A)')
     &   'Time used for initializations: ',TIMI,' seconds'
         WRITE(LUPRI,'(10X,A,F10.2,A)')
     &   'Time used for Y contributions: ',TIMY,' seconds'
         WRITE(LUPRI,'(10X,A,F10.2,A)')
     &   'Time used for T1 back transf.: ',TIMB,' seconds'
         WRITE(LUPRI,'(10X,A,F10.2,A)')
     &   'Time used for AO Fock matrix : ',TIMF,' seconds'
         WRITE(LUPRI,'(10X,A,F10.2,A)')
     &   'Time used for MO transform.  : ',TITR,' seconds'
         WRITE(LUPRI,'(10X,A)')
     &   '-------------------------------------------------'
         WRITE(LUPRI,'(10X,A,F10.2,A,/)')
     &   'Total time                   : ',TIMT,' seconds'
      ENDIF

      RETURN
      END
C  /* Deck cc_chofckp */
      SUBROUTINE CC_CHOFCKP(XLAMDP,XLAMDH,FOCK,WORK,LWORK,ISYMP,ISYMH,
     &                      NUMF)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Calculate NUMF perturbed Fock matrices (in AO basis) from
C              NUMF perturbed Lambda-hole matrices (and one Lambda-particle
C              matrix). The results are accumulated in the FOCK array.
C
C     Only occupied parts of the Lambda matrices should be transferred to
C     this routine.
C
#include <implicit.h>
      DIMENSION XLAMDP(*), XLAMDH(*), FOCK(*), WORK(LWORK)
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      INTEGER IOFF1(8), IOFF2(8)

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CC_CHOFCKP')

      PARAMETER (IOPTR = 2)
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)

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

C     Set symmetry of FOCK.
C     ---------------------

      ISYMF = MULD2H(ISYMP,ISYMH)

C     Set symmetries.
C     ---------------

      ISYCH1 = ISYMH
      ISYCH2 = ISYMP

C     Start Cholesky symmetry loop.
C     -----------------------------

      DO ISYCHO = 1,NSYM

         IF (NUMCHO(ISYCHO) .LE. 0) GO TO 1000
         IF (N2BST(ISYCHO)  .LE. 0) GO TO 1000

         ISYMAK = MULD2H(ISYCHO,ISYCH1)
         ISYMBK = MULD2H(ISYCHO,ISYCH2)

C        Allocation.
C        -----------

         LREAD = MEMRD(1,ISYCHO,IOPTR)

         KCHOL = KEND0
         KREAD = KCHOL + N2BST(ISYCHO)
         KEND1 = KREAD + LREAD
         LWRK1 = LWORK - KEND1 + 1

         IF (LWRK1 .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Insufficient memory in ',SECNAM,
     &      ' - allocation: AO Cholesky'
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need (more than): ',KEND1-1,
     &      'Available       : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Set up batch.
C        -------------

         MINMEM = NUMF*NT1AO(ISYMAK) + NT1AO(ISYMBK)
         IF (MINMEM .GT. 0) THEN
            NVEC = MIN(LWRK1/MINMEM,NUMCHO(ISYCHO))
         ELSE IF (MINMEM .EQ. 0) THEN
            GO TO 1000
         ELSE
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Fatal batch error in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/)')
     &      'Minimum memory required is negative: ',MINMEM
            CALL QUIT('Batch error in '//SECNAM)
         ENDIF

         IF (NVEC .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory for Cholesky batch in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,I10,/)')
     &      'Minimum memory required: ',MINMEM,
     &      'Available for batch    : ',LWRK1,
     &      'Available in total     : ',LWORK
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

            JVEC1 = NVEC*(IBATCH - 1) + 1

C           Complete allocation.
C           --------------------

            KCHO1 = KEND1
            KCHO2 = KCHO1 + NT1AO(ISYMAK)*NUMF*NUMV

C           Set up local index arrays.
C           --------------------------

            CALL IZERO(IOFF1,NSYM)
            CALL IZERO(IOFF2,NSYM)

            ICOUN1 = 0
            ICOUN2 = 0
            DO ISYMK = 1,NSYM
               ISYMA  = MULD2H(ISYMK,ISYMAK)
               IOFF1(ISYMK) = ICOUN1
               ICOUN1 = ICOUN1 + NBAS(ISYMA)*NRHF(ISYMK)*NUMV
               ISYMB = MULD2H(ISYMK,ISYMBK)
               IOFF2(ISYMK) = ICOUN2
               ICOUN2 = ICOUN2 + NBAS(ISYMB)*NRHF(ISYMK)*NUMV
            ENDDO

            DO IVEC = 1,NUMV

C              Read one vector, full square.
C              -----------------------------

               JVEC = JVEC1 + IVEC - 1
               CALL CHO_READN(WORK(KCHOL),JVEC,1,WORK(KIND1),IDUM2,
     &                        ISYCHO,IOPTR,WORK(KREAD),LREAD)

C              Transform.
C              ----------

               DO IF = 1,NUMF
                  DO ISYMK = 1,NSYM

                     ISYMA = MULD2H(ISYMK,ISYMAK)
                     ISYMB = MULD2H(ISYMK,ISYMH)

                     NAK = NBAS(ISYMA)*NRHF(ISYMK)

                     KOFFC = KCHOL + IAODIS(ISYMB,ISYMA)
                     KOFFH = NGLMRH(ISYMH)*(IF - 1)
     &                     + IGLMRH(ISYMB,ISYMK) + 1
                     KOFF1 = KCHO1
     &                     + NT1AO(ISYMAK)*NUMV*(IF - 1)
     &                     + IOFF1(ISYMK) + NAK*(IVEC - 1)
                  
                     NTOTA = MAX(NBAS(ISYMA),1)
                     NTOTB = MAX(NBAS(ISYMB),1)

                     CALL DGEMM('T','N',
     &                        NBAS(ISYMA),NRHF(ISYMK),NBAS(ISYMB),
     &                        ONE,WORK(KOFFC),NTOTB,XLAMDH(KOFFH),NTOTB,
     &                        ZERO,WORK(KOFF1),NTOTA)

                  ENDDO
               ENDDO

               DO ISYMK = 1,NSYM

                  ISYMB = MULD2H(ISYMK,ISYMBK)
                  ISYMA = MULD2H(ISYMK,ISYMP)

                  NBK = NBAS(ISYMB)*NRHF(ISYMK)

                  KOFFC = KCHOL + IAODIS(ISYMA,ISYMB)
                  KOFFP = IGLMRH(ISYMA,ISYMK) + 1
                  KOFF2 = KCHO2 + IOFF2(ISYMK) + NBK*(IVEC - 1)

                  NTOTA = MAX(NBAS(ISYMA),1)
                  NTOTB = MAX(NBAS(ISYMB),1)

                  CALL DGEMM('T','N',
     &                       NBAS(ISYMB),NRHF(ISYMK),NBAS(ISYMA),
     &                       ONE,WORK(KOFFC),NTOTA,XLAMDP(KOFFP),NTOTA,
     &                       ZERO,WORK(KOFF2),NTOTB)

               ENDDO

C              Coulomb part.
C              -------------

               IF (ISYCHO .EQ. ISYMF) THEN

                  DO IF = 1,NUMF
                     FAC = ZERO
                     DO ISYMK = 1,NSYM
                        ISYMB = MULD2H(ISYMK,ISYMBK)
                        NBK   = NBAS(ISYMB)*NRHF(ISYMK)
                        KOFF2 = KCHO2 + IOFF2(ISYMK) + NBK*(IVEC - 1)
                        KOFFH = NGLMRH(ISYMH)*(IF - 1)
     &                        + IGLMRH(ISYMB,ISYMK) + 1
                        FAC = FAC
     &                      + DDOT(NBK,WORK(KOFF2),1,XLAMDH(KOFFH),1)
                     ENDDO
                     FAC   = TWO*FAC
                     KOFFF = N2BST(ISYMF)*(IF - 1) + 1
                     CALL DAXPY(N2BST(ISYMF),FAC,WORK(KCHOL),1,
     &                                           FOCK(KOFFF),1)
                  ENDDO

               ENDIF

            ENDDO

C           Calculate exchange part.
C           ------------------------

            DO IF = 1,NUMF
               DO ISYMK = 1,NSYM

                  ISYMA  = MULD2H(ISYMK,ISYMAK)
                  ISYMB  = MULD2H(ISYMK,ISYMBK)

                  NKJ   = NRHF(ISYMK)*NUMV
                  NTOTA = MAX(NBAS(ISYMA),1)
                  NTOTB = MAX(NBAS(ISYMB),1)

                  KOFF1 = KCHO1 + NT1AO(ISYMAK)*NUMV*(IF - 1)
     &                  + IOFF1(ISYMK)
                  KOFF2 = KCHO2 + IOFF2(ISYMK)
                  KOFFF = N2BST(ISYMF)*(IF - 1) + IAODIS(ISYMA,ISYMB)
     &                  + 1

                  CALL DGEMM('N','T',NBAS(ISYMA),NBAS(ISYMB),NKJ,
     &                       XMONE,WORK(KOFF1),NTOTA,WORK(KOFF2),NTOTB,
     &                       ONE,FOCK(KOFFF),NTOTA)

               ENDDO
            ENDDO

         ENDDO

C        Escape point for empty symmetry block.
C        --------------------------------------

 1000    CONTINUE

      ENDDO

      RETURN
      END
C  /* Deck cc_choetrf */
      SUBROUTINE CC_CHOETRF(XLAMDP,XLAMDH,EIJ,EAB,WMAT,WORK,LWORK,
     &                      ISYMP,ISYMH,ISYMW)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Transform WMAT to MO basis, adding into EIJ and EAB.
C
#include <implicit.h>
      DIMENSION XLAMDP(*), XLAMDH(*), EIJ(*), EAB(*), WMAT(*)
      DIMENSION WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CC_CHOETRF')

      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)

      DO ISYMD = 1,NSYM

         ISYMG = MULD2H(ISYMD,ISYMW)
         ISYMB = MULD2H(ISYMD,ISYMH)
         ISYMA = MULD2H(ISYMG,ISYMP)
         ISYMJ = ISYMB
         ISYMI = ISYMA

C        Allocation.
C        -----------

         KSCR1 = 1
         KEND1 = KSCR1 + NBAS(ISYMG)*MAX(NVIR(ISYMB),NRHF(ISYMJ))
         LWRK1 = LWORK - KEND1 + 1

         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A,I1)')
     &      'Insufficient memory in ',SECNAM,' - symmetry block ',ISYMD
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need     : ',KEND1-1,
     &      'Available: ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

         NTOTG = MAX(NBAS(ISYMG),1)
         NTOTD = MAX(NBAS(ISYMD),1)

         KOFFW = IAODIS(ISYMG,ISYMD) + 1

C        Occupied part.
C        --------------

         NTOTI = MAX(NRHF(ISYMI),1)

         KOFFH = IGLMRH(ISYMD,ISYMJ) + 1
         KOFFP = IGLMRH(ISYMG,ISYMI) + 1
         KOFFE = IMATIJ(ISYMI,ISYMJ) + 1

         CALL DGEMM('N','N',NBAS(ISYMG),NRHF(ISYMJ),NBAS(ISYMD),
     &              ONE,WMAT(KOFFW),NTOTG,XLAMDH(KOFFH),NTOTD,
     &              ZERO,WORK(KSCR1),NTOTG)
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NBAS(ISYMG),
     &              ONE,XLAMDP(KOFFP),NTOTG,WORK(KSCR1),NTOTG,
     &              ONE,EIJ(KOFFE),NTOTI)

C        Virtual part.
C        -------------

         NTOTA = MAX(NVIR(ISYMA),1)

         KOFFH = IGLMVI(ISYMD,ISYMB) + 1
         KOFFP = IGLMVI(ISYMG,ISYMA) + 1
         KOFFE = IMATAB(ISYMA,ISYMB) + 1

         CALL DGEMM('N','N',NBAS(ISYMG),NVIR(ISYMB),NBAS(ISYMD),
     &              ONE,WMAT(KOFFW),NTOTG,XLAMDH(KOFFH),NTOTD,
     &              ZERO,WORK(KSCR1),NTOTG)
         CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMB),NBAS(ISYMG),
     &              ONE,XLAMDP(KOFFP),NTOTG,WORK(KSCR1),NTOTG,
     &              ONE,EAB(KOFFE),NTOTA)

      ENDDO

      RETURN
      END
C  /* Deck cc_choey */
      SUBROUTINE CC_CHOEY(EIJ,EAB,WORK,LWORK,ISYMY,ISIDE,ISYCHL,ITYP,
     &                    NUMCHO,DELY)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Calculate contribution from Y intermediates to E intermediates.
C              ISIDE determines which Y intermediate file to use in the
C              CC_CYI* I/O routines in cc_choio.F; the Cholesky vectors
C              are of symmetry ISYCHL and read via CHO_MOREAD according to the
C              Cholesky type specifier ITYP (see CHO_MOP). Note that NUMCHO
C              (array containing # Cholesky vectors in each symmetry) must
C              be passed as an argument for flexibility.
C
C     Contributions:
C
C        EIJ(ij) = EIJ(ij) + sum(Jc) L(ic,J) * Y(cj,J)
C        EAB(ab) = EAB(ab) - sum(Jk) Y(ak,J) * L(kb,J)
C
C     IF (DELY): the Y intermediate files are deleted after use.
C
#include <implicit.h>
      DIMENSION EIJ(*), EAB(*), WORK(LWORK)
      INTEGER NUMCHO(8)
      LOGICAL DELY
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      CHARACTER*8 SECNAM
      PARAMETER (SECNAM = 'CC_CHOEY')

      INTEGER IOFFY(8), IOFFC(8)

      PARAMETER (IOPEN = -1, IDEL = 0, IKEEP = 1)
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0)

      IF ((ISIDE.EQ.0) .AND. (ISYMY.NE.1)) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Error in ',SECNAM,':'
         WRITE(LUPRI,'(5X,A,I2,/,5X,A,I2,/)')
     &   'Y intermediates must be tot. sym. for ISIDE =',ISIDE,
     &   'Symmetry specified:',ISYMY
         CALL QUIT('Error in '//SECNAM)
      ENDIF

C     Start Cholesky symmetry loop.
C     -----------------------------

      DO ISYCHO = 1,NSYM

         ISYMAK = MULD2H(ISYCHO,ISYMY)
         ISYMBK = MULD2H(ISYCHO,ISYCHL)

         ISYMCJ = ISYMAK
         ISYMCI = ISYMBK

C        Cycle if nothing to do.
C        -----------------------

         IF (NUMCHO(ISYCHO) .LE. 0) GO TO 1000

C        Allocate space for one full Cholesky/Y vector.
C        ----------------------------------------------

         KYIM  = 1
         KCHOL = KYIM  + NT1AM(ISYMAK)
         KEND1 = KCHOL + NT1AM(ISYMBK)
         LWRK1 = LWORK - KEND1 + 1

         IF (LWRK1 .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Insufficient memory in ',SECNAM,' - allocation: Full'
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need (more than): ',KEND1-1,
     &      'Available       : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Set up Cholesky batch.
C        ----------------------

         MINMEM = NT1AM(ISYMAK) + NT1AM(ISYMBK)
         IF (MINMEM .GT. 0) THEN
            NVEC = MIN(LWRK1/MINMEM,NUMCHO(ISYCHO))
         ELSE IF (MINMEM .EQ. 0) THEN
            GO TO 1000
         ELSE
            WRITE(LUPRI,'(//,5X,A,A,I10)')
     &      SECNAM,': MINMEM is negative: ',MINMEM
            CALL QUIT('Error in '//SECNAM)
         ENDIF

         IF (NVEC .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A,I2,A)')
     &      'Insufficient memory for Cholesky batch in ',SECNAM,
     &      '(symmetry',ISYCHO,')'
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,I10,/)')
     &      'Minimum memory required for batch: ',MINMEM,
     &      'Memory available for batch       : ',LWRK1,
     &      'Total memory available           : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

         NBATCH = (NUMCHO(ISYCHO) - 1)/NVEC + 1

C        Open files for this symmetry.
C        -----------------------------

         CALL CHO_MOP(IOPEN,ITYP,ISYCHO,LUCHMO,1,ISYCHL)
         CALL CC_CYIOP(IOPEN,ISYCHO,ISIDE)

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

            KYAK  = KEND1
            KCHBK = KYAK  + NT1AM(ISYMAK)*NUMV

C           Set up local index arrays IOFFY and IOFFC
C           for addressing vec(c,kJ).
C           -----------------------------------------

            CALL IZERO(IOFFY,NSYM)
            CALL IZERO(IOFFC,NSYM)

            ICOUNY = 0
            ICOUNC = 0
            DO ISYMK = 1,NSYM
               ISYMA = MULD2H(ISYMK,ISYMAK)
               ISYMB = MULD2H(ISYMK,ISYMBK)
               IOFFY(ISYMK) = ICOUNY
               IOFFC(ISYMK) = ICOUNC
               ICOUNY = ICOUNY + NVIR(ISYMA)*NRHF(ISYMK)*NUMV
               ICOUNC = ICOUNC + NVIR(ISYMB)*NRHF(ISYMK)*NUMV
            ENDDO

C           Read vectors.
C           Calculate contribution to EIJ.
C           Resort as vec(a,kJ).
C           -------------------------------

            DO IVEC = 1,NUMV

               JVEC = JVEC1 + IVEC - 1
               CALL CC_CYIRDF(WORK(KYIM),1,JVEC,ISYCHO,ISYMY,ISIDE)
               CALL CHO_MOREAD(WORK(KCHOL),NT1AM(ISYMBK),1,JVEC,LUCHMO)

               DO ISYMJ = 1,NSYM

                  ISYMC = MULD2H(ISYMJ,ISYMCJ)
                  ISYMI = MULD2H(ISYMC,ISYMCI)

                  KOFFC = KCHOL + IT1AM(ISYMC,ISYMI)
                  KOFFY = KYIM  + IT1AM(ISYMC,ISYMJ)
                  KOFFE = IMATIJ(ISYMI,ISYMJ) + 1

                  NTOTI = MAX(NRHF(ISYMI),1)
                  NTOTC = MAX(NVIR(ISYMC),1)

                  CALL DGEMM('T','N',
     &                       NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMC),
     &                       ONE,WORK(KOFFC),NTOTC,WORK(KOFFY),NTOTC,
     &                       ONE,EIJ(KOFFE),NTOTI)

               ENDDO

               DO ISYMK = 1,NSYM

                  ISYMA = MULD2H(ISYMK,ISYMAK)
                  NAK   = NVIR(ISYMA)*NRHF(ISYMK)
                  KOFF1 = KYIM + IT1AM(ISYMA,ISYMK)
                  KOFF2 = KYAK + IOFFY(ISYMK) + NAK*(IVEC - 1)
                  CALL DCOPY(NAK,WORK(KOFF1),1,WORK(KOFF2),1)

                  ISYMB = MULD2H(ISYMK,ISYMBK)
                  NBK   = NVIR(ISYMB)*NRHF(ISYMK)
                  KOFF1 = KCHOL + IT1AM(ISYMB,ISYMK)
                  KOFF2 = KCHBK + IOFFC(ISYMK) + NBK*(IVEC - 1)
                  CALL DCOPY(NBK,WORK(KOFF1),1,WORK(KOFF2),1)
 
               ENDDO

            ENDDO

C           Calculate contribution to EAB.
C           ------------------------------

            DO ISYMK = 1,NSYM

               ISYMB = MULD2H(ISYMK,ISYMBK)
               ISYMA = MULD2H(ISYMK,ISYMAK)

               NTOTA = MAX(NVIR(ISYMA),1)
               NTOTB = MAX(NVIR(ISYMB),1)

               KOFFY = KYAK  + IOFFY(ISYMK)
               KOFFC = KCHBK + IOFFC(ISYMK)
               KOFFE = IMATAB(ISYMA,ISYMB) + 1

               CALL DGEMM('N','T',
     &                    NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMK)*NUMV,
     &                    XMONE,WORK(KOFFY),NTOTA,WORK(KOFFC),NTOTB,
     &                    ONE,EAB(KOFFE),NTOTA)

            ENDDO

         ENDDO

C        Close files for this symmetry.
C        Delete if requested.
C        ------------------------------

         CALL CHO_MOP(IKEEP,ITYP,ISYCHO,LUCHMO,1,ISYCHL)
         IF (DELY) THEN
            CALL CC_CYIOP(IDEL,ISYCHO,ISIDE)
         ELSE
            CALL CC_CYIOP(IKEEP,ISYCHO,ISIDE)
         ENDIF

C        Escape point for empty symmetry.
C        --------------------------------

 1000    CONTINUE

      ENDDO

      RETURN
      END
C  /* Deck cc_chotg */
      SUBROUTINE CC_CHOTG(XLAMDP,XLAMDH,X1AM,WORK,LWORK,ISYLAM,ISYMX,
     &                    ISIDE)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Transform AO Cholesky vectors to generalized MO basis
C              [L(ai,J)] for use in Jacobian transformations.
C
C     ISIDE = -1: Left-hand Jacobian transformations.
C     ISIDE =  1: Right-hand Jacobian transformations.
C
#include <implicit.h>
      DIMENSION XLAMDP(*), XLAMDH(*), X1AM(*), WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      CHARACTER*8 SECNAM
      PARAMETER (SECNAM = 'CC_CHOTG')

      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0)

C     Allocation: perturbed Lambda matrix.
C     ------------------------------------

      ISYLM2 = MULD2H(ISYLAM,ISYMX)

      KLAMD2 = 1
      KEND1  = KLAMD2 + NGLMDT(ISYLM2)
      LWRK1  = LWORK  - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   'Insufficient memory in ',SECNAM
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need (more than): ',KEND1-1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Calculate perturbed Lambda matrix.
C     ----------------------------------

      CALL CC_LAMPT(XLAMDP,XLAMDH,X1AM,WORK(KLAMD2),ISYLAM,ISYMX,ISIDE)

C     Transform Cholesky vectors.
C     ---------------------------

      IF (ISIDE .EQ. -1) THEN
         ITYP = 8
         CALL CC_CHOTG1(XLAMDP,XLAMDH,WORK(KLAMD2),WORK(KEND1),LWRK1,
     &                  ISYLAM,ISYLM2,ITYP)
      ELSE IF (ISIDE .EQ. 1) THEN
         ITYP = 9
         CALL CC_CHOTG1(XLAMDH,XLAMDP,WORK(KLAMD2),WORK(KEND1),LWRK1,
     &                  ISYLAM,ISYLM2,ITYP)
      ELSE
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Error in ',SECNAM,':'
         WRITE(LUPRI,'(5X,A,I10)')
     &   'ISIDE must be -1 or +1; ISIDE =',ISIDE
      ENDIF

      RETURN
      END
C  /* Deck cc_chotg1 */
      SUBROUTINE CC_CHOTG1(XLAMDP,XLAMDH,XLAMD2,WORK,LWORK,
     &                      ISYLAM,ISYLM2,ITYP)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Calculate perturbed MO Cholesky vectors L(ia,J), stored
C              as ai.
C              The result vectors are stored on disk as specified
C              by ITYP which is passed on to CHO_MOP for file opening.
C
C     Formula:
C
C     L(ai,J) = sum(g,d) [XLAMD2(gi)*XLAMDH(da) + XLAMDP(gi)*XLAMD2(da)]
C
C                       * L(gd,J)
C
C     where g,d are AO and a,i MO indices.
C
C     Examples of usage (from CC_CHOTG):
C       - left-hand Jacobian transformation: call this routine with
C         the arguments "as is" and ITYP=8.
C       - right-hand Jacobian transformation: call this routine with
C         XLAMDP and XLAMDH interchanged and ITYP=9.
C
#include <implicit.h>
      DIMENSION XLAMDP(*), XLAMDH(*), XLAMD2(*), WORK(LWORK)
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_CHOTG1')

      PARAMETER (IOPEN = -1, IKEEP = 1, IOPTR = 2)
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)

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

C     Get result symmetry.
C     --------------------

      ISYCH = MULD2H(ISYLAM,ISYLM2)

C     Start Cholesky symmetry loop.
C     -----------------------------

      DO ISYCHO = 1,NSYM

         ISYMAI = MULD2H(ISYCHO,ISYCH)

         IF (NUMCHO(ISYCHO) .LE. 0) GO TO 1000
         IF (NT1AM(ISYMAI)  .LE. 0) GO TO 1000

C        Open Cholesky file for this symmetry.
C        -------------------------------------

         CALL CHO_MOP(IOPEN,ITYP,ISYCHO,LUCHAI,1,ISYCH)

C        Calculate size of scratch space needed for read
C        and first half-transformations.
C        -----------------------------------------------

         MAXDI = 0
         DO ISYMI = 1,NSYM
            ISYMG2 = MULD2H(ISYMI,ISYLM2)
            ISYMD2 = MULD2H(ISYMG2,ISYCHO)
            ISYMG  = MULD2H(ISYMI,ISYLAM)
            ISYMD  = MULD2H(ISYMG,ISYCHO)
            ND     = MAX(NBAS(ISYMD2),NBAS(ISYMD))
            NEEDDI = ND*NRHF(ISYMI)
            MAXDI  = MAX(MAXDI,NEEDDI)
         ENDDO
         LREAD = MEMRD(1,ISYCHO,IOPTR)
         LSCR1 = MAX(LREAD,MAXDI)

C        Allocation.
C        -----------

         KCHOL = KEND0
         KSCR1 = KCHOL + N2BST(ISYCHO)
         KEND1 = KSCR1 + LSCR1
         LWRK1 = LWORK - KEND1 + 1

         IF (LWRK1 .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Insufficient memory in ',SECNAM,' - allocation: Cholesky'
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need (more than): ',KEND1-1,
     &      'Available       : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Set up batch.
C        -------------

         MINMEM = NT1AM(ISYMAI)
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

            KCHAI = KEND1

C           Loop over vectors in this batch.
C           --------------------------------

            DO IVEC = 1,NUMV

               JVEC = JVEC1 + IVEC - 1
               CALL CHO_READN(WORK(KCHOL),JVEC,1,WORK(KIND1),IDUM2,
     &                        ISYCHO,IOPTR,WORK(KSCR1),LSCR1)

C              Perform transformation with perturbed i-index.
C              ----------------------------------------------

               DO ISYMI = 1,NSYM

                  ISYMG = MULD2H(ISYMI,ISYLM2)
                  ISYMD = MULD2H(ISYMG,ISYCHO)
                  ISYMA = MULD2H(ISYMD,ISYLAM)

                  NTOTG = MAX(NBAS(ISYMG),1)
                  NTOTD = MAX(NBAS(ISYMD),1)
                  NTOTA = MAX(NVIR(ISYMA),1)

                  KOFFC = KCHOL + IAODIS(ISYMG,ISYMD)
                  KOFP2 = IGLMRH(ISYMG,ISYMI) + 1
                  KOFFH = IGLMVI(ISYMD,ISYMA) + 1
                  KOFFR = KCHAI + NT1AM(ISYMAI)*(IVEC - 1)
     &                  + IT1AM(ISYMA,ISYMI)

                  CALL DGEMM('T','N',
     &                       NBAS(ISYMD),NRHF(ISYMI),NBAS(ISYMG),
     &                       ONE,WORK(KOFFC),NTOTG,XLAMD2(KOFP2),NTOTG,
     &                       ZERO,WORK(KSCR1),NTOTD)
                  CALL DGEMM('T','N',
     &                       NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMD),
     &                       ONE,XLAMDH(KOFFH),NTOTD,WORK(KSCR1),NTOTD,
     &                       ZERO,WORK(KOFFR),NTOTA)

               ENDDO

C              Perform transformation with perturbed a-index.
C              ----------------------------------------------

               DO ISYMI = 1,NSYM

                  ISYMG = MULD2H(ISYMI,ISYLAM)
                  ISYMD = MULD2H(ISYMG,ISYCHO)
                  ISYMA = MULD2H(ISYMD,ISYLM2)

                  NTOTG = MAX(NBAS(ISYMG),1)
                  NTOTD = MAX(NBAS(ISYMD),1)
                  NTOTA = MAX(NVIR(ISYMA),1)

                  KOFFC = KCHOL + IAODIS(ISYMG,ISYMD)
                  KOFFP = IGLMRH(ISYMG,ISYMI) + 1
                  KOFH2 = IGLMVI(ISYMD,ISYMA) + 1
                  KOFFR = KCHAI + NT1AM(ISYMAI)*(IVEC - 1)
     &                  + IT1AM(ISYMA,ISYMI)

                  CALL DGEMM('T','N',
     &                       NBAS(ISYMD),NRHF(ISYMI),NBAS(ISYMG),
     &                       ONE,WORK(KOFFC),NTOTG,XLAMDP(KOFFP),NTOTG,
     &                       ZERO,WORK(KSCR1),NTOTD)
                  CALL DGEMM('T','N',
     &                       NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMD),
     &                       ONE,XLAMD2(KOFH2),NTOTD,WORK(KSCR1),NTOTD,
     &                       ONE,WORK(KOFFR),NTOTA)

               ENDDO

            ENDDO

C           Write this batch to disk.
C           -------------------------

            CALL CHO_MOWRITE(WORK(KCHAI),NT1AM(ISYMAI),NUMV,JVEC1,
     &                       LUCHAI)

         ENDDO

C        Close Cholesky file for this symmetry.
C        --------------------------------------

         CALL CHO_MOP(IKEEP,ITYP,ISYCHO,LUCHAI,1,ISYCH)

C        Escape point for empty symmetry.
C        --------------------------------

 1000    CONTINUE

      ENDDO

      RETURN
      END
C  /* Deck cc_chotg2 */
      SUBROUTINE CC_CHOTG2(XLAMDP,XLAMDH,XLAMD2,WORK,LWORK,
     &                      ISYLAM,ISYLM2,ITYPAI,ITYPJI)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Calculate perturbed MO Cholesky vectors L(ai,J) and L(ji,J).
C              The result vectors are stored on disk as specified
C              by ITYPAI and ITYPJI which is passed on to CHO_MOP for
C              file opening.
C
C     Formulas:
C
C     L(ai,J) = sum(g,d) [XLAMD2(ga)*XLAMDH(di) + XLAMDP(ga)*XLAMD2(di)]
C                       * L(gd,J)
C     L(ji,J) = sum(g,d) [XLAMDP(gj)*XLAMD2(di)]
C                       * L(gd,J)
C
C     where g,d are AO and a,i,j MO indices.
C
C     Note: the ordering of Lambda-p/h differs from cc_chotg1: this routine
C           is written as a "right-hand" type transformation (Lbar), while
C           cc_chotg1 is written as a "left-hand" type (Ltilde). It should
C           be easily fixed by interchanging p/h matrices, though.
C
#include <implicit.h>
      DIMENSION XLAMDP(*), XLAMDH(*), XLAMD2(*), WORK(LWORK)
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_CHOTG2')

      PARAMETER (IOPEN = -1, IKEEP = 1, IOPTR = 2)
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)

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

C     Get result symmetry.
C     --------------------

      ISYCH = MULD2H(ISYLAM,ISYLM2)

C     Start Cholesky symmetry loop.
C     -----------------------------

      DO ISYCHO = 1,NSYM

         ISYMAI = MULD2H(ISYCHO,ISYCH)
         ISYMJI = ISYMAI

         IF (NUMCHO(ISYCHO) .LE. 0) GO TO 1000

C        Calculate size of scratch space needed for read
C        and first half-transformations.
C        -----------------------------------------------

         MAXGI = 0
         DO ISYMI = 1,NSYM
            ISYMD2 = MULD2H(ISYMI,ISYLM2)
            ISYMG2 = MULD2H(ISYMD2,ISYCHO)
            ISYMD  = MULD2H(ISYMI,ISYLAM)
            ISYMG  = MULD2H(ISYMD,ISYCHO)
            NG     = MAX(NBAS(ISYMG2),NBAS(ISYMG))
            NEEDGI = NG*NRHF(ISYMI)
            MAXGI  = MAX(MAXGI,NEEDGI)
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
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need (more than): ',KEND1-1,
     &      'Available       : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Set up batch.
C        -------------

         MINMEM = NT1AM(ISYMAI) + NMATIJ(ISYMJI)
         IF (MINMEM .GT. 0) THEN
            NVEC = MIN(LWRK1/MINMEM,NUMCHO(ISYCHO))
         ELSE IF (MINMEM .EQ. 0) THEN
            GO TO 1000
         ELSE
            WRITE(LUPRI,'(//,5X,A,A,A,I10)')
     &      'MINMEM negative in ',SECNAM,': ',MINMEM
            CALL QUIT('Error in '//SECNAM)
         ENDIF

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

C        Open Cholesky file for this symmetry.
C        -------------------------------------

         CALL CHO_MOP(IOPEN,ITYPAI,ISYCHO,LUCHAI,1,ISYCH)
         CALL CHO_MOP(IOPEN,ITYPJI,ISYCHO,LUCHJI,1,ISYCH)

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

            KCHAI = KEND1
            KCHJI = KCHAI + NT1AM(ISYMAI)*NUMV

C           Loop over vectors in this batch.
C           --------------------------------

            DO IVEC = 1,NUMV

               JVEC = JVEC1 + IVEC - 1
               CALL CHO_READN(WORK(KCHOL),JVEC,1,WORK(KIND1),IDUM2,
     &                        ISYCHO,IOPTR,WORK(KSCR1),LSCR1)

C              Perform transformation with perturbed i-index.
C              ----------------------------------------------

               DO ISYMI = 1,NSYM

                  ISYMD = MULD2H(ISYMI,ISYLM2)
                  ISYMG = MULD2H(ISYMD,ISYCHO)
                  ISYMA = MULD2H(ISYMG,ISYLAM)
                  ISYMJ = ISYMA

                  NTOTG = MAX(NBAS(ISYMG),1)
                  NTOTD = MAX(NBAS(ISYMD),1)
                  NTOTA = MAX(NVIR(ISYMA),1)
                  NTOTJ = MAX(NRHF(ISYMJ),1)

                  KOFFC = KCHOL + IAODIS(ISYMD,ISYMG)
                  KOFF2 = IGLMRH(ISYMD,ISYMI) + 1
                  CALL DGEMM('T','N',
     &                       NBAS(ISYMG),NRHF(ISYMI),NBAS(ISYMD),
     &                       ONE,WORK(KOFFC),NTOTD,XLAMD2(KOFF2),NTOTD,
     &                       ZERO,WORK(KSCR1),NTOTG)

                  KOFFP = IGLMVI(ISYMG,ISYMA) + 1
                  KOFFR = KCHAI + NT1AM(ISYMAI)*(IVEC - 1)
     &                  + IT1AM(ISYMA,ISYMI)
                  CALL DGEMM('T','N',
     &                       NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMG),
     &                       ONE,XLAMDP(KOFFP),NTOTG,WORK(KSCR1),NTOTG,
     &                       ZERO,WORK(KOFFR),NTOTA)

                  KOFFP = IGLMRH(ISYMG,ISYMJ) + 1
                  KOFFR = KCHJI + NMATIJ(ISYMJI)*(IVEC - 1)
     &                  + IMATIJ(ISYMJ,ISYMI)
                  CALL DGEMM('T','N',
     &                       NRHF(ISYMJ),NRHF(ISYMI),NBAS(ISYMG),
     &                       ONE,XLAMDP(KOFFP),NTOTG,WORK(KSCR1),NTOTG,
     &                       ZERO,WORK(KOFFR),NTOTJ)

               ENDDO

C              Perform transformation with perturbed a-index.
C              (No perturbed j-index!!)
C              ----------------------------------------------

               DO ISYMI = 1,NSYM

                  ISYMD = MULD2H(ISYMI,ISYLAM)
                  ISYMG = MULD2H(ISYMD,ISYCHO)
                  ISYMA = MULD2H(ISYMG,ISYLM2)

                  NTOTG = MAX(NBAS(ISYMG),1)
                  NTOTD = MAX(NBAS(ISYMD),1)
                  NTOTA = MAX(NVIR(ISYMA),1)

                  KOFFC = KCHOL + IAODIS(ISYMD,ISYMG)
                  KOFFH = IGLMRH(ISYMD,ISYMI) + 1
                  KOFF2 = IGLMVI(ISYMG,ISYMA) + 1
                  KOFFR = KCHAI + NT1AM(ISYMAI)*(IVEC - 1)
     &                  + IT1AM(ISYMA,ISYMI)

                  CALL DGEMM('T','N',
     &                       NBAS(ISYMG),NRHF(ISYMI),NBAS(ISYMD),
     &                       ONE,WORK(KOFFC),NTOTD,XLAMDH(KOFFH),NTOTD,
     &                       ZERO,WORK(KSCR1),NTOTG)
                  CALL DGEMM('T','N',
     &                       NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMG),
     &                       ONE,XLAMD2(KOFF2),NTOTG,WORK(KSCR1),NTOTG,
     &                       ONE,WORK(KOFFR),NTOTA)

               ENDDO

            ENDDO

C           Write this batch to disk.
C           -------------------------

            CALL CHO_MOWRITE(WORK(KCHAI),NT1AM(ISYMAI),NUMV,JVEC1,
     &                       LUCHAI)
            CALL CHO_MOWRITE(WORK(KCHJI),NMATIJ(ISYMJI),NUMV,JVEC1,
     &                       LUCHJI)

         ENDDO

C        Close Cholesky file for this symmetry.
C        --------------------------------------

         CALL CHO_MOP(IKEEP,ITYPAI,ISYCHO,LUCHAI,1,ISYCH)
         CALL CHO_MOP(IKEEP,ITYPJI,ISYCHO,LUCHJI,1,ISYCH)

C        Escape point for empty symmetry.
C        --------------------------------

 1000    CONTINUE

      ENDDO

      RETURN
      END
C  /* Deck cc_lampt */
      SUBROUTINE CC_LAMPT(XLAMDP,XLAMDH,X1AM,GLAMDA,ISYLAM,ISYMX,ISIDE)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Calculate general lambda matrix according to ISIDE.
C
C              ISIDE = -1 (LH transformation):
C                 GLAMDA(al,p) = [GLAMDP(al,i) ; GLAMDH(al,a)]
C                 GLAMDP(al,i) =   sum(b) XLAMDP(al,b) * X1AM(bi)
C                 GLAMDH(al,a) = - sum(j) XLAMDH(al,j) * X1AM(aj)
C
C              ISIDE = 1 (RH transformation):
C                 GLAMDA(al,p) = [GLAMDH(al,i) ; GLAMDP(al,a)]
C                 GLAMDH(al,i) =   sum(b) XLAMDH(al,b) * X1AM(bi)
C                 GLAMDP(al,a) = - sum(j) XLAMDP(al,j) * X1AM(aj)
C
C     This routine serves the same purpose as CCLR_LAMTRA by Ove Christiansen.
C     The major difference is the output format, and that this routine uses
C     explicit code for left and right perturbed Lambda matrices
C     (ISIDE=-1 or +1). Moreover, this routines is slightly more general, as
C     the unperturbed Lambda matrices need not be tot. sym. (perhaps useless).
C
#include <implicit.h>
      DIMENSION XLAMDP(*), XLAMDH(*),  X1AM(*), GLAMDA(*)
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      CHARACTER*8 SECNAM
      PARAMETER (SECNAM = 'CC_LAMPT')

      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0)

C     Calculate according to ISIDE.
C     -----------------------------

      IF (ISIDE .EQ. -1) THEN

C        Occupied part.
C        --------------

         DO ISYMI = 1,NSYM

            ISYMB = MULD2H(ISYMI,ISYMX)
            ISYMA = MULD2H(ISYMB,ISYLAM)

            KOFF1 = IGLMVI(ISYMA,ISYMB) + 1
            KOFF2 = IT1AM(ISYMB,ISYMI)  + 1
            KOFF3 = IGLMRH(ISYMA,ISYMI) + 1

            NTOTA = MAX(NBAS(ISYMA),1)
            NTOTB = MAX(NVIR(ISYMB),1)

            CALL DGEMM('N','N',NBAS(ISYMA),NRHF(ISYMI),NVIR(ISYMB),
     &                 ONE,XLAMDP(KOFF1),NTOTA,X1AM(KOFF2),NTOTB,
     &                 ZERO,GLAMDA(KOFF3),NTOTA)

         ENDDO

C        Virtual part.
C        -------------

         DO ISYMA = 1,NSYM

            ISYMJ = MULD2H(ISYMA,ISYMX)
            ISYMG = MULD2H(ISYMJ,ISYLAM)

            KOFF1 = IGLMRH(ISYMG,ISYMJ) + 1
            KOFF2 = IT1AM(ISYMA,ISYMJ)  + 1
            KOFF3 = IGLMVI(ISYMG,ISYMA) + 1

            NTOTG = MAX(NBAS(ISYMG),1)
            NTOTA = MAX(NVIR(ISYMA),1)

            CALL DGEMM('N','T',NBAS(ISYMG),NVIR(ISYMA),NRHF(ISYMJ),
     &                 XMONE,XLAMDH(KOFF1),NTOTG,X1AM(KOFF2),NTOTA,
     &                 ZERO,GLAMDA(KOFF3),NTOTG)

         ENDDO

      ELSE IF (ISIDE .EQ. 1) THEN

C        Occupied part.
C        --------------

         DO ISYMI = 1,NSYM

            ISYMB = MULD2H(ISYMI,ISYMX)
            ISYMA = MULD2H(ISYMB,ISYLAM)

            KOFF1 = IGLMVI(ISYMA,ISYMB) + 1
            KOFF2 = IT1AM(ISYMB,ISYMI)  + 1
            KOFF3 = IGLMRH(ISYMA,ISYMI) + 1

            NTOTA = MAX(NBAS(ISYMA),1)
            NTOTB = MAX(NVIR(ISYMB),1)

            CALL DGEMM('N','N',NBAS(ISYMA),NRHF(ISYMI),NVIR(ISYMB),
     &                 ONE,XLAMDH(KOFF1),NTOTA,X1AM(KOFF2),NTOTB,
     &                 ZERO,GLAMDA(KOFF3),NTOTA)

         ENDDO

C        Virtual part.
C        -------------

         DO ISYMA = 1,NSYM

            ISYMJ = MULD2H(ISYMA,ISYMX)
            ISYMG = MULD2H(ISYMJ,ISYLAM)

            KOFF1 = IGLMRH(ISYMG,ISYMJ) + 1
            KOFF2 = IT1AM(ISYMA,ISYMJ)  + 1
            KOFF3 = IGLMVI(ISYMG,ISYMA) + 1

            NTOTG = MAX(NBAS(ISYMG),1)
            NTOTA = MAX(NVIR(ISYMA),1)

            CALL DGEMM('N','T',NBAS(ISYMG),NVIR(ISYMA),NRHF(ISYMJ),
     &                 XMONE,XLAMDP(KOFF1),NTOTG,X1AM(KOFF2),NTOTA,
     &                 ZERO,GLAMDA(KOFF3),NTOTG)

         ENDDO

      ELSE

         WRITE(LUPRI,'(//,5X,A,A)')
     &   'Error in ',SECNAM
         WRITE(LUPRI,'(5X,A,I10,/)')
     &   'ISIDE not recognized: ',ISIDE
         CALL QUIT('Error in '//SECNAM)

      ENDIF

      RETURN
      END
C  /* Deck cc_chocim */
      SUBROUTINE CC_CHOCIM(FOCKD,X1AM,WORK,LWORK,ISYMX,NUMX)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Calculate NUMX C intermediates for use in left-hand
C              Jacobian transformations. The intermediates are written
C              to disk. Use CC_CHOCIM1 for in-core storage.
C
C     Formula:
C
C        C(ai) = sum(bj) [2*t(ai,bj - t(aj,bi)] * X1AM(bj)
C
C     where the doubles are 0th order amplitudes. These amplitudes are
C     never constructed explicitly, but represented by their Cholesky vectors.
C     Thus, if the energy calculation employed decomposed amplitudes
C     [CHOT2 = .TRUE.] the corresponding Cholesky vectors are used here
C     and the orbital denominators [FOCKD] are never used. Otherwise, the
C     amplitudes are decomposed at this stage, storing the resulting vectors
C     for later use.
C
#include <implicit.h>
      DIMENSION FOCKD(*), X1AM(*), WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>
#include <chocc2.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_CHOCIM')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

C     Allocation.
C     -----------

      KCIM  = 1
      KEND1 = KCIM  + NT1AM(ISYMX)*NUMX
      LWRK1 = LWORK - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   'Insufficient memory in ',SECNAM
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need (more than): ',KEND1-1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Calculate the C intermediates.
C     ------------------------------

      CALL DZERO(WORK(KCIM),NT1AM(ISYMX)*NUMX)
      CALL CC_CHOCIM1(FOCKD,X1AM,WORK(KCIM),WORK(KEND1),LWRK1,ISYMX,
     &                NUMX)

C     Write the C intermediates to disk.
C     ----------------------------------

      IF (LOCDBG) THEN
         WRITE(LUPRI,'(/,A,A)')
     &   SECNAM,': writing C intermediates to disk:'
         DO IX = 1,NUMX
            KOFF = KCIM + NT1AM(ISYMX)*(IX - 1)
            XNRM = DSQRT(DDOT(NT1AM(ISYMX),WORK(KOFF),1,WORK(KOFF),1))
            WRITE(LUPRI,'(A,I4,A,I1,A,1P,D22.15)')
     &      'Norm of CIM number',ix,' of sym. ',isymx,': ',XNRM
         ENDDO
         WRITE(LUPRI,*)
      ENDIF

      CALL CHO_IMOP(-1,3,LUCHOC,ISYMX)
      CALL CHO_MOWRITE(WORK(KCIM),NT1AM(ISYMX),NUMX,1,LUCHOC)
      CALL CHO_IMOP(1,3,LUCHOC,ISYMX)

      RETURN
      END
C  /* Deck cc_chocim1 */
      SUBROUTINE CC_CHOCIM1(FOCKD,X1AM,CIM,WORK,LWORK,ISYMX,NUMX)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Calculate NUMX C intermediates for use in left-hand
C              Jacobian transformations. Decompose amplitudes (if
C              necessary) before calculating intermediates.
C
C     Notes: 
C        - The CIM array must be initialized on input.
C        - Specifying NUMX <= 0 will produce the amplitude decomposition
C          only!
C
#include <implicit.h>
      DIMENSION FOCKD(*), X1AM(*), CIM(*), WORK(LWORK)
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>
#include <chocc2.h>
#include <iratdef.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CC_CHOCIM1')

      PARAMETER (ZERO = 0.0D0)

C     Decompose 0th order amplitudes.
C     ===============================

      IF (.NOT. CHOT2C) THEN

         IF (CHOT2) THEN

C           Already done. Copy information.
C           -------------------------------

            IF ((THRCCC.GT.ZERO) .OR. (SPACCC.GT.ZERO)) THEN
               WRITE(LUPRI,'(//,5X,A,A,/,5X,A,A,A)')
     &         '*** NOTICE: Using decomposed CC2 doubles amplitudes ',
     &         'from ground state',
     &         '            optimization in ',SECNAM,'.'
               WRITE(LUPRI,'(5X,A)')
     &         '            User specifications of decomposition are ',
     &         'disregarded!'
               WRITE(LUPRI,'(5X,A)')
     &         '            Decomposition specification used:'
               WRITE(LUPRI,'(5X,A,1P,D15.6)')
     &         '            Threshold  : ',THRCC2
               WRITE(LUPRI,'(5X,A,1P,D15.6)')
     &         '            Span factor: ',SPACC2
               WRITE(LUPRI,'(5X,A,I10)')
     &         '            Max. number of diagonals qualified: ',MXDEC2
               WRITE(LUPRI,'(5X,A,I10,//)')
     &         '            Max. number of prev. vectors read : ',NCHRD2
            ENDIF

            MXDECC = MXDEC2
            NCHRDC = NCHRD2
            THRCCC = THRCC2
            SPACCC = SPACC2
            DO ISYCHO = 1,NSYM
               NCHMOC(ISYCHO) = NCHMO2(ISYCHO)
            ENDDO

         ELSE

C           Decompose.
C           ----------

            IF (THRCCC .LE. ZERO) THEN
               THRCCC = THRCOM
            ENDIF
            IF (SPACCC .LE. ZERO) THEN
               SPACCC = SPAN
            ENDIF
            MXDECC = MXDEC2
            NCHRDC = NCHRD2
            THZCCC = THZCC2

            IMAT = 3

            KDIANL = 1
            KCOLUM = KDIANL + 3*NSYM
            KINDIA = KCOLUM + (NSYM - 1)/IRAT   + 1
            KEND1  = KINDIA + (MXDECC - 1)/IRAT + 1
            LWRK1  = LWORK  - KEND1

            IF (LWRK1 .LE. 0) THEN
               WRITE(LUPRI,'(//,5X,A,A,A)')
     &         'Insufficient memory in ',SECNAM,' - allocation: decomp.'
               WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &         'Need (more than): ',KEND1,
     &         'Available       : ',LWORK
               CALL QUIT('Insufficient memory in '//SECNAM)
            ENDIF

            CALL CC_DECMO(WORK(KDIANL),NCHMOC,WORK(KCOLUM),WORK(KINDIA),
     &                    THRCCC,SPACCC,
     &                    FCC2S,.FALSE.,MXDECC,NCHRDC,NT1AM,IPRINT,IMAT,
     &                    NSYM,THZCCC,FOCKD,WORK(KEND1),LWRK1)

         ENDIF

C        Set CHOT2C flag.
C        ----------------

         CHOT2C = .TRUE.

      ENDIF

C     Construct the C intermediates.
C     ==============================

      IF (NUMX .GT. 0) THEN
         CALL CC_CHOCIM2(X1AM,CIM,WORK,LWORK,ISYMX,NUMX)
      ENDIF

      RETURN
      END
C  /* Deck cc_chocim2 */
      SUBROUTINE CC_CHOCIM2(X1AM,CIM,WORK,LWORK,ISYMX,NUMX)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Calculate NUMX C intermediates from decomposed T2 amplitudes.
C
C     Note: All contributions appear to have the wrong sign, but this is
C           in fact correct, since the decomposition represents, not T2,
C           but -T2 (for positive definiteness).
C
C     Formula:
C
C        CIM(ai) = sum(bj) [2*t(ai,bj)-t(aj,bi)] * X1AM(bj)
C
#include <implicit.h>
      DIMENSION X1AM(*), CIM(*), WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
#include <chocc2.h>
#include <priunit.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CC_CHOCIM2')
 
      INTEGER IOFF1(8), IOFF2(8)

      PARAMETER (IOPEN = -1, IKEEP = 1, ITYP = 6)
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, XMTWO = -2.0D0)

      IF (NUMX .GT. 0) THEN

         DO ISYCHO = 1,NSYM

            IF (NCHMOC(ISYCHO) .LE. 0) GO TO 1000
            IF (NT1AM(ISYCHO)  .LE. 0) GO TO 1000

            ISYMIJ = MULD2H(ISYCHO,ISYMX)

C           Open Cholesky file.
C           -------------------

            CALL CHO_MOP(IOPEN,ITYP,ISYCHO,LUCHMO,1,1)

C           Allocate space for 1 vector.
C           ----------------------------

            KCHOL = 1
            IF (ISYCHO .EQ. ISYMX) THEN
               KCOUL = KCHOL + NT1AM(ISYCHO)
               KEND1 = KCOUL + NUMX
            ELSE
               KEND1 = KCHOL + NT1AM(ISYCHO)
            ENDIF
            LWRK1 = LWORK - KEND1 + 1

            IF (LWRK1 .LE. 0) THEN
               WRITE(LUPRI,'(//,5X,A,A,A)')
     &         'Insufficient memory in ',SECNAM,' - allocation: Chol.'
               NEED = KEND1 - 1 + NT1AM(ISYCHO) + NUMX*NMATIJ(ISYMIJ)
               WRITE(LUPRI,'(5X,A,I10,A,I1,A,/,5X,A,I10,/)')
     &         'Minimum memory required: ',NEED,
     &         ' (Cholesky symmetry ',ISYCHO,')',
     &         'Total memory available : ',LWORK
               CALL QUIT('Insufficient memory in '//SECNAM)
            ENDIF

C           Set up batch.
C           -------------

            MINMEM = NT1AM(ISYCHO) + NUMX*NMATIJ(ISYMIJ)
            NVEC   = MIN(LWRK1/MINMEM,NCHMOC(ISYCHO))

            IF (NVEC .LE. 0) THEN
               WRITE(LUPRI,'(//,5X,A,A)')
     &         'Insufficient memory for Cholesky batch in ',SECNAM
               WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,I10,/)')
     &         'Minimum memory required: ',MINMEM,
     &         'Available for batching : ',LWRK1,
     &         'Total memory available : ',LWORK
               CALL QUIT('Insufficient memory in '//SECNAM)
            ENDIF

            NBATCH = (NCHMOC(ISYCHO) - 1)/NVEC + 1

C           Start batch loop.
C           -----------------

            DO IBATCH = 1,NBATCH

               NUMV = NVEC
               IF (IBATCH .EQ. NBATCH) THEN
                  NUMV = NCHMOC(ISYCHO) - NVEC*(NBATCH - 1)
               ENDIF

               JVEC1 = NVEC*(IBATCH - 1) + 1

C              Set up local index arrays.
C              --------------------------

               ICOUN1 = 0
               ICOUN2 = 0
               DO ISYMJ = 1,NSYM
                  ISYMA = MULD2H(ISYMJ,ISYCHO)
                  ISYMI = MULD2H(ISYMJ,ISYMIJ)
                  IOFF1(ISYMJ) = ICOUN1
                  IOFF2(ISYMJ) = ICOUN2
                  ICOUN1 = ICOUN1 + NVIR(ISYMA)*NRHF(ISYMJ)*NUMV
                  ICOUN2 = ICOUN2 + NRHF(ISYMI)*NRHF(ISYMJ)*NUMV
               ENDDO
               NIJVEC = ICOUN2

C              Complete allocation.
C              --------------------

               KCHAJ = KEND1
               KCHIJ = KCHAJ + NT1AM(ISYCHO)*NUMV

C              Read, one vector at a time.
C              Calculate Coulomb part if symmetry fits.
C              Sort and transform for exchange part.
C              ----------------------------------------

               DO IVEC = 1,NUMV

                  JVEC = JVEC1 + IVEC - 1
                  CALL CHO_MOREAD(WORK(KCHOL),NT1AM(ISYCHO),1,JVEC,
     &                            LUCHMO)

                  IF (ISYCHO .EQ. ISYMX) THEN

                     NTOBJ = MAX(NT1AM(ISYMX),1)
                     CALL DGEMV('T',NT1AM(ISYMX),NUMX,
     &                          ONE,X1AM,NTOBJ,WORK(KCHOL),1,
     &                          ZERO,WORK(KCOUL),1)

                     DO IX = 1,NUMX
                        FAC   = XMTWO*WORK(KCOUL+IX-1)
                        KOFFC = NT1AM(ISYMX)*(IX - 1) + 1
                        CALL DAXPY(NT1AM(ISYMX),FAC,WORK(KCHOL),1,
     &                                              CIM(KOFFC),1)
                     ENDDO

                  ENDIF

                  DO ISYMJ = 1,NSYM

                     ISYMA = MULD2H(ISYMJ,ISYCHO)

                     NAJ = NVIR(ISYMA)*NRHF(ISYMJ)

                     KOFF1 = KCHOL + IT1AM(ISYMA,ISYMJ)
                     KOFF2 = KCHAJ + IOFF1(ISYMJ) + NAJ*(IVEC - 1)

                     CALL DCOPY(NAJ,WORK(KOFF1),1,WORK(KOFF2),1)

                  ENDDO

                  DO IX = 1,NUMX
                     DO ISYMJ = 1,NSYM

                        ISYMB = MULD2H(ISYMJ,ISYMX)
                        ISYMI = MULD2H(ISYMB,ISYCHO)

                        NIJ = NRHF(ISYMI)*NRHF(ISYMJ)

                        KOFF1 = KCHOL + IT1AM(ISYMB,ISYMI)
                        KOFF2 = NT1AM(ISYMX)*(IX - 1)
     &                        + IT1AM(ISYMB,ISYMJ) + 1
                        KOFF3 = KCHIJ + NIJVEC*(IX - 1)
     &                        + IOFF2(ISYMJ) + NIJ*(IVEC - 1)

                        NTB = MAX(NVIR(ISYMB),1)
                        NTI = MAX(NRHF(ISYMI),1)

                        CALL DGEMM('T','N',
     &                             NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMB),
     &                             ONE,WORK(KOFF1),NTB,X1AM(KOFF2),NTB,
     &                             ZERO,WORK(KOFF3),NTI)

                     ENDDO
                  ENDDO

               ENDDO

C              Contract to give final exchange contribution.
C              ---------------------------------------------

               DO IX = 1,NUMX
                  DO ISYMJ = 1,NSYM

                     ISYMA = MULD2H(ISYMJ,ISYCHO)
                     ISYMI = MULD2H(ISYMJ,ISYMIJ)

                     KOFF1 = KCHAJ + IOFF1(ISYMJ)
                     KOFF2 = KCHIJ + NIJVEC*(IX - 1) + IOFF2(ISYMJ)
                     KOFF3 = NT1AM(ISYMX)*(IX - 1) + IT1AM(ISYMA,ISYMI)
     &                     + 1

                     NTA = MAX(NVIR(ISYMA),1)
                     NTI = MAX(NRHF(ISYMI),1)

                     NJV = NRHF(ISYMJ)*NUMV

                     CALL DGEMM('N','T',NVIR(ISYMA),NRHF(ISYMI),NJV,
     &                          ONE,WORK(KOFF1),NTA,WORK(KOFF2),NTI,
     &                          ONE,CIM(KOFF3),NTA)

                  ENDDO
               ENDDO

            ENDDO

C           Close Cholesky file.
C           --------------------

            CALL CHO_MOP(IKEEP,ITYP,ISYCHO,LUCHMO,1,1)

C           Escape point for empty symmetry.
C           --------------------------------

 1000       CONTINUE

         ENDDO

      ENDIF

      RETURN
      END
C  /* Deck cc_chobcktr */
      SUBROUTINE CC_CHOBCKTR(XLAMDH,X1AM,XLAMD2,ISYMX,NUMX) 
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Backtransform singles amplitudes to obtain "perturbed"
C              occupied Lambda-hole matrix,
C
C              XLAMD2(g,i) = sum(a) XLAMDH(g,a) * X1AM(a,i)
C
C              for NUMX amplitude arrays of symmetry ISYMX. XLAMD2 is
C              ordered according to IGLMRH(*,*) and stored consecutively.
C
#include <implicit.h>
      DIMENSION XLAMDH(*), X1AM(*), XLAMD2(*)
#include <ccorb.h>
#include <ccsdsym.h>

      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)

      DO IX = 1,NUMX
         DO ISYMI = 1,NSYM

            ISYMA = MULD2H(ISYMI,ISYMX)
            ISYMG = ISYMA

            KOFF1 = IGLMVI(ISYMG,ISYMA)    + 1
            KOFF2 = NT1AM(ISYMX)*(IX - 1)  + IT1AM(ISYMA,ISYMI)  + 1
            KOFF3 = NGLMRH(ISYMX)*(IX - 1) + IGLMRH(ISYMG,ISYMI) + 1

            NTOTG = MAX(NBAS(ISYMG),1)
            NTOTA = MAX(NVIR(ISYMA),1)
         
            CALL DGEMM('N','N',NBAS(ISYMG),NRHF(ISYMI),NVIR(ISYMA),
     &                 ONE,XLAMDH(KOFF1),NTOTG,X1AM(KOFF2),NTOTA,
     &                 ZERO,XLAMD2(KOFF3),NTOTG)

         ENDDO
      ENDDO

      RETURN
      END
C  /* Deck cc_chotrprp */
      SUBROUTINE CC_CHOTRPRP(XIJ,XAB,WORK,LWORK,ISYMX,ITYP1,ITYP2,
     &                       ISYCH1,NUMCHO)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Transform Cholesky vectors with property integrals
C
C              L2(ai,J) = sum(j) L1(aj,J) * XIJ(ji)
C                       - sum(b) XAB(ab)  * L1(bi,J)
C
C     ISYMX : Symmetry of property integrals.
C     ITYP1 : Type specifier for L1 as recognized in CHO_MOP.
C     ITYP2 : Type specifier for L2 as recognized in CHO_MOP.
C     ISYCH1: Symmetry of L1.
C     NUMCHO: Array containing the number of vectors (of L1 and L2)
C
C     The L1 vectors must be available on disk as specified through ITYP1,
C     and the L2 vectors will are written to disk as specified through ITYP2.
C
#include <implicit.h>
      DIMENSION XIJ(*), XAB(*), WORK(LWORK)
      INTEGER NUMCHO(8)
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CC_CHOTRPRP')

      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0)

C     Set result symmetry.
C     --------------------

      ISYCH2 = MULD2H(ISYCH1,ISYMX)

C     Start Cholesky symmetry loop.
C     -----------------------------

      DO ISYCHO = 1,NSYM

         ISYMAI = MULD2H(ISYCHO,ISYCH2)
         ISYMBI = MULD2H(ISYCHO,ISYCH1)
         ISYMAJ = ISYMBI

         IF (NUMCHO(ISYCHO) .LE. 0) GO TO 1000
         IF (NT1AM(ISYMAI)  .LE. 0) GO TO 1000
         IF (NT1AM(ISYMBI)  .LE. 0) GO TO 1000

C        Set up batch. Half for original vectors, half for result vectors.
C        -----------------------------------------------------------------

         MINMEM = NT1AM(ISYMBI) + NT1AM(ISYMAI)
         NVEC   = MIN(LWORK/MINMEM,NUMCHO(ISYCHO))

         IF (NVEC .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory for Cholesky batch in ',SECNAM
            WRITE(LUPRI,'(5X,A,I1)')
     &      'Cholesky symmetry: ',ISYCHO
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Minimum memory required: ',MINMEM,
     &      'Memory available       : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

         NBATCH = (NUMCHO(ISYCHO) - 1)/NVEC + 1

C        Open files.
C        -----------

         CALL CHO_MOP(-1,ITYP1,ISYCHO,LUCH1,1,ISYCH1)
         CALL CHO_MOP(-1,ITYP2,ISYCHO,LUCH2,1,ISYCH2)

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

            KCHO1 = 1
            KCHO2 = KCHO1 + NT1AM(ISYMBI)*NUMV

C           Read original vectors.
C           ----------------------

            CALL CHO_MOREAD(WORK(KCHO1),NT1AM(ISYMBI),NUMV,JVEC1,LUCH1)

C           Occupied transformation.
C           ------------------------

            DO IVEC = 1,NUMV
               DO ISYMI = 1,NSYM

                  ISYMJ = MULD2H(ISYMI,ISYMX)
                  ISYMA = MULD2H(ISYMI,ISYMAI)

                  NTOTJ = MAX(NRHF(ISYMJ),1)
                  NTOTA = MAX(NVIR(ISYMA),1)

                  KOFF1 = KCHO1 + NT1AM(ISYMAJ)*(IVEC - 1)
     &                  + IT1AM(ISYMA,ISYMJ)
                  KOFFX = IMATIJ(ISYMJ,ISYMI) + 1
                  KOFF2 = KCHO2 + NT1AM(ISYMAI)*(IVEC - 1)
     &                  + IT1AM(ISYMA,ISYMI)

                  CALL DGEMM('N','N',
     &                       NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMJ),
     &                       ONE,WORK(KOFF1),NTOTA,XIJ(KOFFX),NTOTJ,
     &                       ZERO,WORK(KOFF2),NTOTA)

               ENDDO
            ENDDO

C           Virtual transformation.
C           -----------------------

            DO IVEC = 1,NUMV
               DO ISYMI = 1,NSYM

                  ISYMB = MULD2H(ISYMI,ISYMBI)
                  ISYMA = MULD2H(ISYMI,ISYMAI)

                  NTOTA = MAX(NVIR(ISYMA),1)
                  NTOTB = MAX(NVIR(ISYMB),1)

                  KOFFX = IMATAB(ISYMA,ISYMB) + 1
                  KOFF1 = KCHO1 + NT1AM(ISYMBI)*(IVEC - 1)
     &                  + IT1AM(ISYMB,ISYMI)
                  KOFF2 = KCHO2 + NT1AM(ISYMAI)*(IVEC - 1)
     &                  + IT1AM(ISYMA,ISYMI)

                  CALL DGEMM('N','N',
     &                       NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMB),
     &                       XMONE,XAB(KOFFX),NTOTA,WORK(KOFF1),NTOTB,
     &                       ONE,WORK(KOFF2),NTOTA)

               ENDDO
            ENDDO

C           Write transformed vectors.
C           --------------------------

            CALL CHO_MOWRITE(WORK(KCHO2),NT1AM(ISYMAI),NUMV,JVEC1,LUCH2)

         ENDDO

C        Close files.
C        ------------

         CALL CHO_MOP(1,ITYP1,ISYCHO,LUCH1,1,ISYCH1)
         CALL CHO_MOP(1,ITYP2,ISYCHO,LUCH2,1,ISYCH2)

C        Escape point for empty symmetry.
C        --------------------------------

 1000    CONTINUE

      ENDDO

      RETURN
      END
C  /* Deck cc_chotrfov */
      SUBROUTINE CC_CHOTRFOV(XLAMDP,XLAMDH,XAO,XIJ,XAB,
     &                       WORK,LWORK,ISYMP,ISYMH,ISYMX)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Transform AO one-electron integrals to MO blocks
C              ij [occupied] and ab [virtual].
C
#include <implicit.h>
      DIMENSION XLAMDP(*), XLAMDH(*), XAO(*), XIJ(*), XAB(*)
      DIMENSION WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CC_CHOTRFOV')

      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)

C     Check memory.
C     -------------

      NEED = 0
      DO ISYMD = 1,NSYM
         ISYMG = MULD2H(ISYMD,ISYMX)
         ISYMI = MULD2H(ISYMD,ISYMH)
         ISYMB = ISYMI
         NTST1 = NBAS(ISYMG)*MAX(NRHF(ISYMI),NVIR(ISYMB))
         NEED  = MAX(NEED,NTST1)
      ENDDO

      IF (NEED .GT. LWORK) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   'Insufficient memory in ',SECNAM
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need     : ',NEED,
     &   'Available: ',LWORK
         CALL QUIT('Insufficient meory in '//SECNAM)
      ENDIF

C     Calculate XIJ.
C     --------------

      DO ISYMD = 1,NSYM

         ISYMG = MULD2H(ISYMD,ISYMX)
         ISYMI = MULD2H(ISYMD,ISYMH)
         ISYMJ = MULD2H(ISYMG,ISYMP)

         NTOTG = MAX(NBAS(ISYMG),1)
         NTOTD = MAX(NBAS(ISYMD),1)
         NTOTJ = MAX(NRHF(ISYMJ),1)

         KOFFX = IAODIS(ISYMG,ISYMD) + 1

         KOFFH = IGLMRH(ISYMD,ISYMI) + 1
         CALL DGEMM('N','N',NBAS(ISYMG),NRHF(ISYMI),NBAS(ISYMD),
     &              ONE,XAO(KOFFX),NTOTG,XLAMDH(KOFFH),NTOTD,
     &              ZERO,WORK,NTOTG)

         KOFFP = IGLMRH(ISYMG,ISYMJ) + 1
         KOFJI = IMATIJ(ISYMJ,ISYMI) + 1
         CALL DGEMM('T','N',NRHF(ISYMJ),NRHF(ISYMI),NBAS(ISYMG),
     &              ONE,XLAMDP(KOFFP),NTOTG,WORK,NTOTG,
     &              ZERO,XIJ(KOFJI),NTOTJ)

      ENDDO

C     Calculate XAB.
C     --------------

      DO ISYMD = 1,NSYM

         ISYMG = MULD2H(ISYMD,ISYMX)
         ISYMB = MULD2H(ISYMD,ISYMH)
         ISYMA = MULD2H(ISYMG,ISYMP)

         NTOTG = MAX(NBAS(ISYMG),1)
         NTOTD = MAX(NBAS(ISYMD),1)
         NTOTA = MAX(NVIR(ISYMA),1)

         KOFFX = IAODIS(ISYMG,ISYMD) + 1

         KOFFH = IGLMVI(ISYMD,ISYMB) + 1
         CALL DGEMM('N','N',NBAS(ISYMG),NVIR(ISYMB),NBAS(ISYMD),
     &              ONE,XAO(KOFFX),NTOTG,XLAMDH(KOFFH),NTOTD,
     &              ZERO,WORK,NTOTG)

         KOFFP = IGLMVI(ISYMG,ISYMA) + 1
         KOFAB = IMATAB(ISYMA,ISYMB) + 1
         CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMB),NBAS(ISYMG),
     &              ONE,XLAMDP(KOFFP),NTOTG,WORK,NTOTG,
     &              ZERO,XAB(KOFAB),NTOTA)

      ENDDO

      RETURN
      END
C  /* Deck cc_chotrfed */
      SUBROUTINE CC_CHOTRFED(XLAMDP,XLAMDH,XAO,XAI,XIA,
     &                       WORK,LWORK,ISYMP,ISYMH,ISYMX)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Transform AO one-electron integrals to MO blocks
C              ai [excitation] and ia [de-excitation, stored as ai]
C
#include <implicit.h>
      DIMENSION XLAMDP(*), XLAMDH(*), XAO(*), XAI(*), XIA(*)
      DIMENSION WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CC_CHOTRFED')

      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)

C     Check memory.
C     -------------

      NEED = 0
      DO ISYMD = 1,NSYM
         ISYMG = MULD2H(ISYMD,ISYMX)
         ISYMI = MULD2H(ISYMD,ISYMH)
         ISYMJ = MULD2H(ISYMG,ISYMP)
         NTST1 = NBAS(ISYMG)*NRHF(ISYMI)
         NTST2 = NBAS(ISYMD)*NRHF(ISYMJ)
         NEED  = MAX(NEED,NTST1,NTST2)
      ENDDO

      IF (NEED .GT. LWORK) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   'Insufficient memory in ',SECNAM
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need     : ',NEED,
     &   'Available: ',LWORK
         CALL QUIT('Insufficient meory in '//SECNAM)
      ENDIF

C     Calculate XAI.
C     --------------

      DO ISYMD = 1,NSYM

         ISYMG = MULD2H(ISYMD,ISYMX)
         ISYMI = MULD2H(ISYMD,ISYMH)
         ISYMA = MULD2H(ISYMG,ISYMP)

         NTOTG = MAX(NBAS(ISYMG),1)
         NTOTD = MAX(NBAS(ISYMD),1)
         NTOTA = MAX(NVIR(ISYMA),1)

         KOFFX = IAODIS(ISYMG,ISYMD) + 1

         KOFFH = IGLMRH(ISYMD,ISYMI) + 1
         CALL DGEMM('N','N',NBAS(ISYMG),NRHF(ISYMI),NBAS(ISYMD),
     &              ONE,XAO(KOFFX),NTOTG,XLAMDH(KOFFH),NTOTD,
     &              ZERO,WORK,NTOTG)

         KOFFP = IGLMVI(ISYMG,ISYMA) + 1
         KOFAI = IT1AM(ISYMA,ISYMI)  + 1
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMG),
     &              ONE,XLAMDP(KOFFP),NTOTG,WORK,NTOTG,
     &              ZERO,XAI(KOFAI),NTOTA)

      ENDDO

C     Calculate XIA.
C     --------------

      DO ISYMD = 1,NSYM

         ISYMG = MULD2H(ISYMD,ISYMX)
         ISYMA = MULD2H(ISYMD,ISYMH)
         ISYMI = MULD2H(ISYMG,ISYMP)

         NTOTG = MAX(NBAS(ISYMG),1)
         NTOTD = MAX(NBAS(ISYMD),1)
         NTOTA = MAX(NVIR(ISYMA),1)

         KOFFX = IAODIS(ISYMG,ISYMD) + 1

         KOFFP = IGLMRH(ISYMG,ISYMI) + 1
         CALL DGEMM('T','N',NBAS(ISYMD),NRHF(ISYMI),NBAS(ISYMG),
     &              ONE,XAO(KOFFX),NTOTG,XLAMDP(KOFFP),NTOTG,
     &              ZERO,WORK,NTOTD)

         KOFFH = IGLMVI(ISYMD,ISYMA) + 1
         KOFAI = IT1AM(ISYMA,ISYMI)  + 1
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMD),
     &              ONE,XLAMDH(KOFFH),NTOTD,WORK,NTOTD,
     &              ZERO,XIA(KOFAI),NTOTA)

      ENDDO

      RETURN
      END
C  /* Deck cc_chotrfmo */
      SUBROUTINE CC_CHOTRFMO(XLAMDP,XLAMDH,XAO,XIJ,XAI,XIA,XAB,
     &                       WORK,LWORK,ISYMP,ISYMH,ISYMX)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Transform the AO property integrals in XAO to MO
C              basis.
C
#include <implicit.h>
      DIMENSION XLAMDP(*), XLAMDH(*), XAO(*), XIJ(*), XAI(*), XIA(*)
      DIMENSION XAB(*), WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CC_CHOTRFMO')

      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)

C     Check memory.
C     -------------

      NEED = 0
      DO ISYMD = 1,NSYM
         ISYMG = MULD2H(ISYMD,ISYMX)
         ISYMI = MULD2H(ISYMD,ISYMH)
         ISYMB = ISYMI
         NTST1 = NBAS(ISYMG)*MAX(NRHF(ISYMI),NVIR(ISYMB))
         NEED  = MAX(NEED,NTST1)
      ENDDO

      IF (NEED .GT. LWORK) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   'Insufficient memory in ',SECNAM
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need     : ',NEED,
     &   'Available: ',LWORK
         CALL QUIT('Insufficient meory in '//SECNAM)
      ENDIF

C     Calculate XIJ and XAI.
C     ----------------------

      DO ISYMD = 1,NSYM

         ISYMG = MULD2H(ISYMD,ISYMX)
         ISYMI = MULD2H(ISYMD,ISYMH)
         ISYMJ = MULD2H(ISYMG,ISYMP)
         ISYMA = ISYMJ

         NTOTG = MAX(NBAS(ISYMG),1)
         NTOTD = MAX(NBAS(ISYMD),1)
         NTOTJ = MAX(NRHF(ISYMJ),1)
         NTOTA = MAX(NVIR(ISYMA),1)

         KOFFX = IAODIS(ISYMG,ISYMD) + 1

         KOFFH = IGLMRH(ISYMD,ISYMI) + 1
         CALL DGEMM('N','N',NBAS(ISYMG),NRHF(ISYMI),NBAS(ISYMD),
     &              ONE,XAO(KOFFX),NTOTG,XLAMDH(KOFFH),NTOTD,
     &              ZERO,WORK,NTOTG)

         KOFFP = IGLMRH(ISYMG,ISYMJ) + 1
         KOFJI = IMATIJ(ISYMJ,ISYMI) + 1
         CALL DGEMM('T','N',NRHF(ISYMJ),NRHF(ISYMI),NBAS(ISYMG),
     &              ONE,XLAMDP(KOFFP),NTOTG,WORK,NTOTG,
     &              ZERO,XIJ(KOFJI),NTOTJ)

         KOFFP = IGLMVI(ISYMG,ISYMA) + 1
         KOFAI = IT1AM(ISYMA,ISYMI)  + 1
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMG),
     &              ONE,XLAMDP(KOFFP),NTOTG,WORK,NTOTG,
     &              ZERO,XAI(KOFAI),NTOTA)

      ENDDO

C     Calculate XIA and XAB.
C     ----------------------

      DO ISYMD = 1,NSYM

         ISYMG = MULD2H(ISYMD,ISYMX)
         ISYMB = MULD2H(ISYMD,ISYMH)
         ISYMI = MULD2H(ISYMG,ISYMP)
         ISYMA = ISYMI

         NTOTG = MAX(NBAS(ISYMG),1)
         NTOTD = MAX(NBAS(ISYMD),1)
         NTOTA = MAX(NVIR(ISYMA),1)
         NTOTB = MAX(NVIR(ISYMB),1)

         KOFFX = IAODIS(ISYMG,ISYMD) + 1

         KOFFH = IGLMVI(ISYMD,ISYMB) + 1
         CALL DGEMM('N','N',NBAS(ISYMG),NVIR(ISYMB),NBAS(ISYMD),
     &              ONE,XAO(KOFFX),NTOTG,XLAMDH(KOFFH),NTOTD,
     &              ZERO,WORK,NTOTG)

         KOFFP = IGLMRH(ISYMG,ISYMI) + 1
         KOFBI = IT1AM(ISYMB,ISYMI)  + 1
         CALL DGEMM('T','N',NVIR(ISYMB),NRHF(ISYMI),NBAS(ISYMG),
     &              ONE,WORK,NTOTG,XLAMDP(KOFFP),NTOTG,
     &              ZERO,XIA(KOFBI),NTOTB)

         KOFFP = IGLMVI(ISYMG,ISYMA) + 1
         KOFAB = IMATAB(ISYMA,ISYMB) + 1
         CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMB),NBAS(ISYMG),
     &              ONE,XLAMDP(KOFFP),NTOTG,WORK,NTOTG,
     &              ZERO,XAB(KOFAB),NTOTA)

      ENDDO

      RETURN
      END
C  /* Deck cc_choadd */
      SUBROUTINE CC_CHOADD(ITYP1,ITYP2,ISYM,NUMCHO,WORK,LWORK,
     &                     ALPHA,BETA)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Add the vectors specified through ITYP1 to those specified
C              by ITYP2 and store in the same file as the latter! ISYM is
C              the overall symmetry of the vectors.
C
C     Formula:
C
C        L2(ai,J) <-- ALPHA*L1(ai,J) + BETA*L2(ai,J)
C
#include <implicit.h>
      INTEGER NUMCHO(8)
      DIMENSION WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_CHOADD')

      LOGICAL L1CON, L2SCL

      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)

C     Set up calculation. Return if nothing to do.
C     --------------------------------------------

      IF (ALPHA .EQ. ZERO) THEN
         L1CON = .FALSE.
      ELSE
         L1CON = .TRUE.
      ENDIF

      IF (BETA .EQ. ONE) THEN
         L2SCL = .FALSE.
      ELSE
         L2SCL = .TRUE.
      ENDIF

      IF ((.NOT.L1CON) .AND. (.NOT.L2SCL)) RETURN

C     Add in loop over symmetries.
C     ----------------------------

      DO ISYCHO = 1,NSYM

         ISYMAI = MULD2H(ISYCHO,ISYM)

         IF (NUMCHO(ISYCHO) .LE. 0) GO TO 1000
         IF (NT1AM(ISYMAI)  .LE. 0) GO TO 1000

C        Open files.
C        -----------

         IF (L1CON) CALL CHO_MOP(-1,ITYP1,ISYCHO,LUCHO1,1,ISYM)
         CALL CHO_MOP(-1,ITYP2,ISYCHO,LUCHO2,1,ISYM)

C        Set up batch.
C        -------------

         IF (L1CON) THEN
            MINMEM = 2*NT1AM(ISYMAI)
         ELSE
            MINMEM = NT1AM(ISYMAI) 
         ENDIF
         NVEC = MIN(LWORK/MINMEM,NUMCHO(ISYCHO))

         IF (NVEC .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need at least: ',MINMEM,
     &      'Available    : ',LWORK
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

            JVEC1 = NVEC*(IBATCH - 1) + 1

C           Complete allocation.
C           --------------------

            LAIV = NT1AM(ISYMAI)*NUMV

            KCHO2 = 1
            KCHO1 = KCHO2 + LAIV

C           Read and scale L2 vectors.
C           --------------------------

            CALL CHO_MOREAD(WORK(KCHO2),NT1AM(ISYMAI),NUMV,JVEC1,LUCHO2)
            IF (L2SCL) CALL DSCAL(LAIV,BETA,WORK(KCHO2),1)

C           Read L1 vectors and add.
C           ------------------------

            IF (L1CON) THEN
               CALL CHO_MOREAD(WORK(KCHO1),NT1AM(ISYMAI),NUMV,JVEC1,
     &                         LUCHO1)
               CALL DAXPY(LAIV,ALPHA,WORK(KCHO1),1,WORK(KCHO2),1)
            ENDIF

C           Write result to disk.
C           ---------------------

            CALL CHO_MOWRITE(WORK(KCHO2),NT1AM(ISYMAI),NUMV,JVEC1,
     &                       LUCHO2)

         ENDDO

C        Close files.
C        ------------

         IF (L1CON) CALL CHO_MOP(1,ITYP1,ISYCHO,LUCHO1,1,ISYM)
         CALL CHO_MOP(1,ITYP2,ISYCHO,LUCHO2,1,ISYM)

C        Escape point for empty symmetry.
C        --------------------------------

 1000    CONTINUE

      ENDDO

      RETURN
      END
C  /* Deck cc_choebar */
      SUBROUTINE CC_CHOEBAR(XLAMDP,XLAMDH,FIA,X1AM,SIGMA,EIJ,EAB,
     &                      WORK,LWORK,ISYMX,ITYPY,ISYMY,
     &                      DELY,INCLRF,SVYCON)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Calculate Ebar intermediates for F-matrix transformation
C              and, if requested (flag INCLRF), the reference contribution.
C              The Y intermediates are read from disk according to specifier
C              ITYPY. If DELY = .TRUE., the Y intermediate files will be
C              deleted after use. If SVYCON = .TRUE., copies of the Y
C              contributions are saved at the end of EIJ and EAB (making these
C              effectively 2-dimensional).
C
C     Formulas:
C
C     E(ij) = sum(gd) LamP(gi) * W(gd) * LamH(dj)
C           + sum(j)  FIA(ib) * X(bj)
C           + sum(bJ) L(ib,J) * Y(bj,J)
C
C     E(ab) = sum(gd) LamP(ga) * W(gd) * LamH(db)
C           - sum(j)  X(aj) * FIA(jb)
C           - sum(jJ) Y(aj,J) * L(jb,J)
C
C     SIGMA(ai) = 2 * sum(gd) LamP(gi) * W(gd) * LamH(da)
C
C     where
C
C     W(gd) = sum(ck) [2*(gd|kc)-(gc|kd)] * X(ck)
C
C     which is calculated through a backtransformation of the c-index
C     of X1AM.
C
#include <implicit.h>
      DIMENSION XLAMDP(*), XLAMDH(*), FIA(*), X1AM(*), SIGMA(*)
      DIMENSION EIJ(*), EAB(*), WORK(LWORK)
      LOGICAL DELY, INCLRF, SVYCON
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CC_CHOEBAR')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)

C     Check symmetries.
C     -----------------

      IF (ISYMY .NE. ISYMX) THEN
         WRITE(LUPRI,'(5X,A,A,A)')
     &   'Symmetry mismatch in ',SECNAM,':'
         WRITE(LUPRI,'(5X,A,I2,A,I2)')
     &   'ISYMX =',ISYMX,' != ISYMY =',ISYMY
         CALL QUIT('Symmetry mismatch in '//SECNAM)
      ENDIF

      ISYMW = ISYMX

      IF (LOCDBG) THEN
         X1NM2 = DDOT(NT1AM(ISYMX),X1AM,1,X1AM,1)
         X1NRM = DSQRT(X1NM2)
         WRITE(LUPRI,*) '   Entering ',SECNAM
         WRITE(LUPRI,*) '   ISYMX [R vector ]: ',ISYMX
         WRITE(LUPRI,*) '   ISYMY [Y interm.]: ',ISYMY
         WRITE(LUPRI,*) '   ITYPY [Y interm.]: ',ITYPY
         WRITE(LUPRI,*) '   DELY             : ',DELY
         WRITE(LUPRI,*) '   INCLRF           : ',INCLRF
         WRITE(LUPRI,*) '   SVYCON           : ',SVYCON
         WRITE(LUPRI,*) '   Sq. norm, R vec. : ',X1NM2
         WRITE(LUPRI,*) '   2---norm, R vec. : ',X1NRM
         CALL FLSHFO(LUPRI)
      ENDIF

C     Initialize EIJ and EAB.
C     -----------------------

      CALL DZERO(EIJ,NMATIJ(ISYMX))
      CALL DZERO(EAB,NMATAB(ISYMX))

C     Calculate contributions from Y intermediates.
C     ---------------------------------------------

      CALL CC_CHOEY(EIJ,EAB,WORK,LWORK,ISYMY,ITYPY,1,1,NUMCHO,DELY)

C     If requested, save copies of these contributions at the end of
C     result arrays.
C     --------------------------------------------------------------

      IF (SVYCON) THEN
         KOFF1 = 1
         KOFF2 = NMATIJ(ISYMX) + 1
         CALL DCOPY(NMATIJ(ISYMX),EIJ(KOFF1),1,EIJ(KOFF2),1)
         KOFF2 = NMATAB(ISYMX) + 1
         CALL DCOPY(NMATAB(ISYMX),EAB(KOFF1),1,EAB(KOFF2),1)
      ENDIF

c     xxnrm = dsqrt(ddot(nmatij(isymx),eij,1,eij,1))
c     yynrm = dsqrt(ddot(nmatab(isymx),eab,1,eab,1))
c     write(lupri,*) '...2-norm of Y3 con. to Eij.: ',xxnrm
c     write(lupri,*) '...2-norm of Y3 con. to Eab.: ',yynrm

C     Calculate contributions from FIA.
C     ---------------------------------

      DO ISYMJ = 1,NSYM

         ISYMB = MULD2H(ISYMJ,ISYMX)
         ISYMI = ISYMB

         NTOTB = MAX(NVIR(ISYMB),1)
         NTOTI = MAX(NRHF(ISYMI),1)

         KOFFF = IT1AM(ISYMB,ISYMI)  + 1
         KOFFX = IT1AM(ISYMB,ISYMJ)  + 1
         KOFFE = IMATIJ(ISYMI,ISYMJ) + 1

         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMB),
     &              ONE,FIA(KOFFF),NTOTB,X1AM(KOFFX),NTOTB,
     &              ONE,EIJ(KOFFE),NTOTI)

      ENDDO

      DO ISYMJ = 1,NSYM

         ISYMA = MULD2H(ISYMJ,ISYMX)
         ISYMB = ISYMJ

         NTOTA = MAX(NVIR(ISYMA),1)
         NTOTB = MAX(NVIR(ISYMB),1)

         KOFFX = IT1AM(ISYMA,ISYMJ)  + 1
         KOFFF = IT1AM(ISYMB,ISYMJ)  + 1
         KOFFE = IMATAB(ISYMA,ISYMB) + 1

         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMJ),
     &              XMONE,X1AM(KOFFX),NTOTA,FIA(KOFFF),NTOTB,
     &              ONE,EAB(KOFFE),NTOTA)

      ENDDO

C     Allocation.
C     -----------

      KWMAT = 1
      KLAMH = KWMAT + N2BST(ISYMW)
      KEND1 = KLAMH + NGLMRH(ISYMX)
      LWRK1 = LWORK - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation W'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need     : ',KEND1-1,
     &   'Available: ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Backtransform X1AM.
C     -------------------

      DO ISYMK = 1,NSYM

         ISYMC = MULD2H(ISYMK,ISYMX)
         ISYMA = ISYMC

         NTOTA = MAX(NBAS(ISYMA),1)
         NTOTC = MAX(NVIR(ISYMC),1)

         KOFFH = IGLMVI(ISYMA,ISYMC) + 1
         KOFFX = IT1AM(ISYMC,ISYMK)  + 1
         KOFF2 = KLAMH + IGLMRH(ISYMA,ISYMK)

         CALL DGEMM('N','N',NBAS(ISYMA),NRHF(ISYMK),NVIR(ISYMC),
     &              ONE,XLAMDH(KOFFH),NTOTA,X1AM(KOFFX),NTOTC,
     &              ZERO,WORK(KOFF2),NTOTA)

      ENDDO

C     Calculate W.
C     ------------

      CALL DZERO(WORK(KWMAT),N2BST(ISYMW))
      CALL CC_CHOFCKP(XLAMDP,WORK(KLAMH),WORK(KWMAT),WORK(KEND1),LWRK1,
     &                1,ISYMX,1)

C     Reset memory.
C     -------------

      KEND1 = KLAMH
      LWRK1 = LWORK - KEND1 + 1

C     Transform W to MO basis, adding into EIJ and EAB.
C     -------------------------------------------------

      CALL CC_CHOETRF(XLAMDP,XLAMDH,EIJ,EAB,WORK(KWMAT),
     &                WORK(KEND1),LWRK1,1,1,ISYMW)

C     If requested, calculate reference contribution.
C     -----------------------------------------------

      IF (INCLRF) THEN
         DO ISYMD = 1,NSYM

            ISYMG = MULD2H(ISYMD,ISYMW)
            ISYMA = ISYMD
            ISYMI = ISYMG

            NEED = NBAS(ISYMD)*NRHF(ISYMI)

            IF (NEED .GT. LWRK1) THEN
               CALL QUIT('Insufficient memory in '//SECNAM//' (Ref.)')
            ENDIF

            NTOTG = MAX(NBAS(ISYMG),1)
            NTOTD = MAX(NBAS(ISYMD),1)
            NTOTA = MAX(NVIR(ISYMA),1)

            KOFFW = KWMAT + IAODIS(ISYMG,ISYMD)
            KOFFP = IGLMRH(ISYMG,ISYMI) + 1

            CALL DGEMM('T','N',NBAS(ISYMD),NRHF(ISYMI),NBAS(ISYMG),
     &                 ONE,WORK(KOFFW),NTOTG,XLAMDP(KOFFP),NTOTG,
     &                 ZERO,WORK(KEND1),NTOTD)

            KOFFH = IGLMVI(ISYMD,ISYMA) + 1
            KOFFS = IT1AM(ISYMA,ISYMI)  + 1

            CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMD),
     &                 ONE,XLAMDH(KOFFH),NTOTD,WORK(KEND1),NTOTD,
     &                 ZERO,SIGMA(KOFFS),NTOTA)

         ENDDO
         CALL DSCAL(NT1AM(ISYMX),TWO,SIGMA,1)
      ENDIF

      RETURN
      END
