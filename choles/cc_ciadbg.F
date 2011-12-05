C  /* Deck cc_ciadbg */
      SUBROUTINE CC_CIADBG(WORK,LWORK,IDBG)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Test cc_cia.
C
C     IDBG = 0: Do all tests.
C          = 1: Check (ia|jb) integrals.
C          = 2: Check 'artificial' integrals constructed from (ia|jb).
C
C     WARNING: This routine requires *lots* of memory and, furthermore,
C              produces a *lot* of print. Since it calls CC_CIA hundreds,
C              even thousands, of times (the more symmetry, the more calls),
C              it will also spend a *lot* of CPU time. So, use it only for
C              small debug cases....
C
#include <implicit.h>
      DIMENSION WORK(LWORK)
#include <priunit.h>
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_CIADBG')

      LOGICAL DOALL

C     Print header.
C     -------------

      CALL AROUND(SECNAM//': TESTING CC_CIA')

C     Test that this is a Cholesky run.
C     ---------------------------------

      IF (.NOT. CHOINT) THEN
         WRITE(LUPRI,'(5X,A,A,A,//)')
     &   'This is not a Cholesky run! ',
     &   SECNAM,' returns without further processing.'
         RETURN
      ENDIF

C     Initializations.
C     ----------------

      DOALL = IDBG .EQ. 0
      NERR  = 0
      NCALL = 0

C     Calculate tot. sym. (ia|jb) integrals.
C     --------------------------------------

      IF ((IDBG.EQ.1) .OR. DOALL) THEN
         CALL CC_CIADBG1(WORK,LWORK,NERR,NCALL)
      ENDIF

C     Calculate artificial non-tot. sym. integrals.
C     ---------------------------------------------

      IF ((IDBG.EQ.2) .OR. DOALL) THEN
         CALL CC_CIADBG2(WORK,LWORK,NERR,NCALL)
         DO ISYM = 1,NSYM
            CALL CHO_MOP(-1,7,ISYM,LUDUM,1,NSYM)
            CALL CHO_MOP(0,7,ISYM,LUDUM,1,NSYM)
         ENDDO
      ENDIF

C     Print total number of errors.
C     -----------------------------

      CALL HEADER('Final Result from '//SECNAM,-1)
      WRITE(LUPRI,'(5X,A,I10)')
     & 'Total number of calls to subroutine CC_CIA:',NCALL
      WRITE(LUPRI,'(5X,A,I10)')
     & 'Total number of errors detected in  CC_CIA:',NERR
      IF (NERR .EQ. 0) THEN
         WRITE(LUPRI,'(5X,A)')
     &   '- everything seems fine!'
      ELSE
         WRITE(LUPRI,'(5X,A)')
     &   '- you need to do some debugging! See above for detailed info.'
      ENDIF

C     Print exit message.
C     -------------------

      CALL AROUND('END OF '//SECNAM)

      RETURN
      END
C  /* Deck cc_ciadbg1 */
      SUBROUTINE CC_CIADBG1(WORK,LWORK,NERR,NCALT)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Test cc_cia for (ia|jb) integral calculation.
C
C     This routine tests only for tot. sym. integral generation,
C     trying to invoke cc_cia in all conceivable ways.
C
C     NERR is *updated* with the number of errors detected.
C
#include <implicit.h>
      DIMENSION WORK(LWORK)
#include <priunit.h>
#include <maxorb.h>
#include <ccisvi.h>
#include <ccdeco.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <dccsdsym.h>
#include <ccorb.h>
#include <ciarc.h>
#include <chomp2.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CC_CIADBG1')

      INTEGER LSYM(8)

      LOGICAL CHMOSV, MP2SSV

      PARAMETER (TINY = 1.0D-14, TINY2 = 1.0D-15)
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (HALF = 0.5D0)

C     Print start message.
C     --------------------

      CALL HEADER('Starting '//SECNAM,-1)
      WRITE(LUPRI,'(5X,A,/)')
     & 'Testing CC_CIA assembly of (ia|jb) integrals.'

      XLW = ONE*LWORK
      IF (TWO*XT2SQ(1) .GT. XLW) THEN
         WRITE(LUPRI,'(5X,A,/)')
     &   'Oooops! Insufficient memory for this test routine....exiting!'
         RETURN
      ENDIF 

C     Set print level in CC_CIA.
C     --------------------------

C     IPRCIA = IPRLVL + 1
      IPRCIA = IPRLVL - 1

C     Turn off calculation of min. mem. requirement.
C     Turn on initialization of integral array in CC_CIA.
C     ---------------------------------------------------

      GETMNM = .FALSE.
      INXINT = .TRUE.

C     Initialize counters.
C     --------------------

      IERR  = 0
      NCALL = 0

C     Initialize LSYM (for opening files).
C     ------------------------------------

      DO ISYM = 1,NSYM
         LSYM(ISYM) = ISYM
      ENDDO

C     Get the MO vectors L(ia,J).
C     ---------------------------

      LCMO  = NLAMDS
      IF (SKIPTR) LCMO = 0
      KEMO  = 1
      KTDIA = KEMO  + NORBTS
      KCMO  = KTDIA + 1
      KEND0 = KCMO  + LCMO
      LWRK0 = LWORK - KEND0 + 1

      IF (LWRK0 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,
     &   ' - Allocation: MO'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need     : ',KEND0-1,
     &   'Available: ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

      CALL CHO_RDSIR(DUM1,DUM2,WORK(KEMO),WORK(KCMO),WORK(KEND0),LWRK0,
     &               .FALSE.,.TRUE.,(.NOT.SKIPTR))

      CALL CHO_TRFAI(WORK(KCMO),WORK(KEND0),LWRK0)

C     Reset memory; CMO no longer needed.
C     -----------------------------------

      KEND0 = KCMO
      LWRK0 = LWORK - KCMO + 1

C     Allocate space for target (ia|jb) integrals.
C     --------------------------------------------

      KXINT = KEND0
      KEND1 = KXINT + NT2SQ(1)
      LWRK1 = LWORK - KEND1 + 1

      IF (LWRK1. LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,
     &   ' - Allocation: MO'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need     : ',KEND1-1,
     &   'Available: ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Calculate integrals through the MP2 module.
C     -------------------------------------------

      CHMOSV = CHOMO
      MP2SSV = MP2SAV
      CHOMO  = .FALSE.
      MP2SAV = .FALSE.
      CALL MP2CHO_2(WORK(KTDIA),WORK(KEMO),WORK(KXINT),
     &              WORK(KEND1),LWRK1,EMP2)
      CHOMO  = CHMOSV
      MP2SAV = MP2SSV

      WRITE(LUPRI,'(/,5X,A,F20.10,/)')
     & 'MP2 energy correction obtained from target integrals: ',EMP2

C     ============================
C     START TESTING CC_CIA MODULE.
C     ============================

C     Diagonal-scaling-or-not-loop
C     ----------------------------

      IDOSCL = -1
  111 IDOSCL = IDOSCL + 1
      IF (IDOSCL .GT. 1) GO TO 112
      IF (IDOSCL .EQ. 1) THEN
         DO ISYMBJ = 1,NSYM
            KOFFX = KXINT + IT2SQ(ISYMBJ,ISYMBJ)
            INCRX = NT1AM(ISYMBJ) + 1
            CALL DSCAL(NT1AM(ISYMBJ),HALF,WORK(KOFFX),INCRX)
         ENDDO
      ENDIF

C     Test 1.1: L1=L2 option, full calculation.
C     -----------------------------------------

      KXCIA = KEND1
      KEND2 = KXCIA + NT2SQ(1)
      LWRK2 = LWORK - KEND2 + 1

      IF (LWRK2. LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,
     &   ' - Allocation: Test 1'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need     : ',KEND2-1,
     &   'Available: ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

      DO ICLCIA = 1,2
         NPASS  = 0
         DIAFAC = ONE
         IF (IDOSCL .EQ. 1) DIAFAC = HALF
         CALL CC_CIADBGRC1(1,1,.TRUE.,.FALSE.,.FALSE.,.FALSE.)
         ITYP1 = 1
         ITYP2 = 1
         IF (ICLCIA .EQ. 2) THEN
            LW2S = LWRK2
            NEED = 0
            NVEC = NUMCHO(1)
            DO ISYM = 1,NSYM
               NEEDS = NT1AM(ISYM)
               NEED  = MAX(NEED,NEEDS)
               NVEC  = MIN(NVEC,NUMCHO(ISYM))
            ENDDO
            MULM  = MIN(NVEC,5)
            NEED  = MULM*(NEED + 1)
            LWRK2 = MIN(LWRK2,NEED)
         ENDIF
         CALL CC_CIA(WORK(KXCIA),WORK(KEND2),LWRK2,DIAFAC,KLAST,NPASS)
         NCALL = NCALL + 1
         IF (ICLCIA .EQ. 2) THEN
            LWRK2 = LW2S
         ENDIF

         CALL HEADER('Test 1.1 Result',-1)
         IF (IDOSCL .EQ. 1) THEN
            WRITE(LUPRI,'(5X,A,1P,D20.10,A)')
     &      '(diagonal scaled with ',DIAFAC,')'
         ENDIF
         IF (ICLCIA .EQ. 2) THEN
            WRITE(LUPRI,'(5X,A)')
     &      '(Cholesky batching was forced in CC_CIA)'
         ENDIF
         CALL CC_DIADBGANL(WORK(KXINT),WORK(KXCIA),1,TRGNRM,CIANRM,
     &                     DIFNRM,DIFRMS,RELNRM,IDIVZ)
         IF ((IDIVZ.EQ.1) .OR. (IDIVZ.EQ.12)) THEN
            WRITE(LUPRI,'(5X,A)')
     &      '- the difference RMS was assigned to avoid division with 0'
         ENDIF
         IF ((IDIVZ.EQ.2) .OR. (IDIVZ.EQ.12)) THEN
            WRITE(LUPRI,'(5X,A)')
     &     '- the rel. diff. norm was assigned to avoid division with 0'
         ENDIF
       WRITE(LUPRI,'(5X,A,1P,D20.10,/,5X,A,1P,D20.10,/,5X,A,1P,D20.10)')
     &    'Norm of (ia|jb) [target]: ',TRGNRM,
     &    'Norm of (ia|jb) [cc_cia]: ',CIANRM,
     &    'Difference              : ',TRGNRM-CIANRM
         WRITE(LUPRI,*)
       WRITE(LUPRI,'(5X,A,1P,D20.10,/,5X,A,1P,D20.10,/,5X,A,1P,D20.10)')
     &    'Difference norm         : ',DIFNRM,
     &    'Difference RMS          : ',DIFRMS,
     &    'Relative difference norm: ',RELNRM
         IF ((RELNRM.GT.TINY) .OR. (DIFRMS.GT.TINY)) THEN
            WRITE(LUPRI,'(5X,A,/)')
     &      '- That counts as an error!'
            IERR = IERR + 1
         ENDIF
      ENDDO

C     Test 1.4: L1 = L2 option, batch calculation,
C                               no symmetrization.
C     --------------------------------------------

      DO ICLCIA = 1,2
         DIFNM2 = ZERO
         CIANM2 = ZERO
         NPASS  = 0
         DIAFAC = ONE
         IF (IDOSCL .EQ. 1) DIAFAC = HALF
         ISYCH1 = 1
         ISYCH2 = 1
         ISYMX  = MULD2H(ISYCH1,ISYCH2)
         LIDEN  = .TRUE.
         SYMTRZ = .FALSE.
         CIAMIO = .FALSE.
         ITYP1  = 1
         ITYP2  = 1
         DO ISYM = 1,NSYM
            NTOVEC(ISYM) = NUMCHO(ISYM)
         ENDDO
         IF (NSYM .EQ. 1) THEN
            NVIRB = MAX(NVIRT/2,1)
            NVIRA = MAX(NVIRT/3,1)
         ELSE
            NVIRB = MAX(NVIRT/NSYM,1)
            NVIRA = MAX(NVIRT/(NSYM+1),1)
         ENDIF
         NBATB = (NVIRT - 1)/NVIRB + 1
         NBATA = (NVIRT - 1)/NVIRA + 1
         DO IBATB = 1,NBATB
            NUMB = NVIRB
            IF (IBATB .EQ. NBATB) THEN
               NUMB = NVIRT - NVIRB*(NBATB - 1)
            ENDIF
            IB1 = NVIRB*(IBATB - 1) + 1
            CALL CC_CIADBGSET(IB1,NUMB,IOFB1,LVIRB,NX1AMB,IX1AMB)
            DO IBATA = 1,NBATA
               NUMA = NVIRA
               IF (IBATA .EQ. NBATA) THEN
                  NUMA = NVIRT - NVIRA*(NBATA - 1)
               ENDIF
               IA1 = NVIRA*(IBATA - 1) + 1
               CALL CC_CIADBGSET(IA1,NUMA,IOFA1,LVIRA,NX1AMA,IX1AMA)
               CALL IZERO(IX2SQ,64)
               DO ISYM = 1,NSYM
                  ICOUNT = 0
                  DO ISYMBJ = 1,NSYM
                     ISYMAI = MULD2H(ISYMBJ,ISYM)
                     IX2SQ(ISYMAI,ISYMBJ) = ICOUNT
                     ICOUNT = ICOUNT + NX1AMA(ISYMAI)*NX1AMB(ISYMBJ)
                  ENDDO
                  IF (ISYM .EQ. ISYMX) NX2SQ = ICOUNT
               ENDDO
               IF (ICLCIA .GT. 1) THEN
                  LW2S = LWRK2
                  NEED = 0
                  NVEC = NUMCHO(1)
                  DO ISYCHO = 1,NSYM
                     NEEDS = NT1AM(ISYCHO)
     &                     + NX1AMA(ISYCHO) + NX1AMB(ISYCHO)
                     NEED  = MAX(NEED,NEEDS)
                     NVEC  = MIN(NVEC,NUMCHO(ISYCHO))
                  ENDDO
                  MULM  = MIN(NVEC,5)
                  NEED  = MULM*NEED
                  LWRK2 = MIN(LWRK2,NEED)
               ENDIF
               CALL CC_CIA(WORK(KXCIA),WORK(KEND2),LWRK2,DIAFAC,KLAST,
     &                     NPASS)
               NCALL = NCALL + 1
               IF (ICLCIA .GT. 1) THEN
                  LWRK2 = LW2S
               ENDIF
               CALL CC_CIADBGANL2(WORK(KXINT),WORK(KXCIA),ISYMX,
     &                            CIANM2,DIFNM2,NXCIA,JERR,TINY2)
               CALL HEADER('Test 1.4 Partial Result',-1)
               IF (IDOSCL .EQ. 1) THEN
                  WRITE(LUPRI,'(5X,A,1P,D20.10,A)')
     &            '(diagonal scaled with ',DIAFAC,')'
               ENDIF
               IA2 = IA1 + NUMA - 1
               IB2 = IB1 + NUMB - 1
               IF (ICLCIA .GT. 1) THEN
                  WRITE(LUPRI,'(5X,A)')
     &            '(Cholesky batching was forced in CC_CIA)'
               ENDIF
               WRITE(LUPRI,'(5X,A,I4,A,I4,/,5X,A,I4,A,I4)')
     &         'This is a-batch number',IBATA,' of',NBATA,
     &         '--- for b-batch number',IBATB,' of',NBATB
               WRITE(LUPRI,'(5X,A,I4,1X,I1,A,I4,1X,I1)')
     &         'a-index [a,sym.a] runs from',IA1,ISVI(IA1),' to',
     &         IA2,ISVI(IA2)
               WRITE(LUPRI,'(5X,A,I4,1X,I1,A,I4,1X,I1)')
     &         'b-index [b,sym.b] runs from',IB1,ISVI(IB1),' to',
     &         IB2,ISVI(IB2)
               WRITE(LUPRI,'(5X,A,I10)')
     &         'Discrepancies in this batch: ',JERR
               IF (JERR .NE. 0) THEN
                  IERR = IERR + 1
                  WRITE(LUPRI,'(5X,A)')
     &            '- That counts as an error!'
               ENDIF
               IF (NXCIA .NE. NX2SQ) THEN
                  WRITE(LUPRI,'(5X,A)')
     &            'Something is terribly wrong: NXCIA != NX2SQ'
                  WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     &            'NXCIA = ',NXCIA,
     &            'NX2SQ = ',NX2SQ
                  WRITE(LUPRI,'(5X,A)')
     &            '- That counts as an error!'
                  IERR = IERR + 1
               ENDIF
            ENDDO
         ENDDO
         CALL HEADER('Test 1.4 Global Result',-1)
         IF (IDOSCL .EQ. 1) THEN
            WRITE(LUPRI,'(5X,A,1P,D20.10,A)')
     &      '(diagonal scaled with ',DIAFAC,')'
         ENDIF
         IF (ICLCIA .GT. 1) THEN
            WRITE(LUPRI,'(5X,A)')
     &      '(Cholesky batching was forced in CC_CIA)'
         ENDIF
         TRGNRM = DSQRT(DDOT(NT2SQ(ISYMX),WORK(KXINT),1,WORK(KXINT),1))
         CIANRM = DSQRT(CIANM2)
         DIFNRM = DSQRT(DIFNM2)
         DIFRMS = DSQRT(DIFNM2/XT2SQ(1))
         RELNRM = DIFNRM/TRGNRM
       WRITE(LUPRI,'(5X,A,1P,D20.10,/,5X,A,1P,D20.10,/,5X,A,1P,D20.10)')
     &   'Norm of (ia|jb) [target]: ',TRGNRM,
     &   'Norm of (ia|jb) [cc_cia]: ',CIANRM,
     &   'Difference              : ',TRGNRM-CIANRM
         WRITE(LUPRI,*)
       WRITE(LUPRI,'(5X,A,1P,D20.10,/,5X,A,1P,D20.10,/,5X,A,1P,D20.10)')
     &   'Difference norm         : ',DIFNRM,
     &   'Difference RMS          : ',DIFRMS,
     &   'Relative difference norm: ',RELNRM
         IF ((RELNRM.GT.TINY) .OR. (DIFRMS.GT.TINY)) THEN
            WRITE(LUPRI,'(5X,A,/)')
     &      '- That counts as an error!'
            IERR = IERR + 1
         ENDIF
      ENDDO

C     Back loop for diag. scaling.
C     ----------------------------

      GO TO 111

C     End diag. scaling loop.
C     -----------------------

  112 CONTINUE

C     Update global error and call counters.
C     --------------------------------------

      NERR  = NERR + IERR
      NCALT = NCALT + NCALL

C     Print end message.
C     ------------------

      CALL HEADER('End of '//SECNAM,-1)
      WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     & 'Total number of calls to CC_CIA in this routine:',NCALL,
     & 'Accumulated number of calls                    :',NCALT
      WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     & 'Number of errors detected for (ia|jb) assembly:',IERR,
     & 'Accumulated number of errors                  :',NERR

      RETURN
      END
C  /* Deck cc_ciadbganl2 */
      SUBROUTINE CC_CIADBGANL2(XINT,XCIA,ISYMX,CIANM2,DIFNM2,NTST,NERR,
     &                         TOL)
C
C     Analyze for errors.
C
#include <implicit.h>
      DIMENSION XINT(*), XCIA(*)
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ciarc.h>

      INTEGER AI, BJ, AIBJ

      PARAMETER (ZERO = 0.0D0)

      NTST = 0
      NERR = 0

      DO ISYMBJ = 1,NSYM

         ISYMAI = MULD2H(ISYMBJ,ISYMX)
         IF (NT1AM(ISYMBJ)  .LE. 0) GO TO 999
         IF (NT1AM(ISYMAI)  .LE. 0) GO TO 999
         IF (NX1AMB(ISYMBJ) .LE. 0) GO TO 999
         IF (NX1AMA(ISYMAI) .LE. 0) GO TO 999

         DO ISYMJ = 1,NSYM

            ISYMB = MULD2H(ISYMJ,ISYMBJ)
            IF (LVIRB(ISYMB) .LE. 0) GO TO 998

            DO J = 1,NRHF(ISYMJ)
               DO LB = 1,LVIRB(ISYMB)

                  LBJ = IX1AMB(ISYMB,ISYMJ)
     &                + LVIRB(ISYMB)*(J - 1) + LB
                  B   = IOFB1(ISYMB) + LB - 1
                  BJ  = IT1AM(ISYMB,ISYMJ)
     &                + NVIR(ISYMB)*(J - 1) + B

                  DO ISYMI = 1,NSYM

                     ISYMA = MULD2H(ISYMI,ISYMAI)
                     IF (LVIRA(ISYMA) .LE. 0) GO TO 997

                     DO I = 1,NRHF(ISYMI)
                        DO LA = 1,LVIRA(ISYMA)

                           LAI = IX1AMA(ISYMA,ISYMI)
     &                         + LVIRA(ISYMA)*(I - 1) + LA
                           A   = IOFA1(ISYMA) + LA - 1
                           AI  = IT1AM(ISYMA,ISYMI)
     &                         + NVIR(ISYMA)*(I - 1) + A

                           LAIBJ = IX2SQ(ISYMAI,ISYMBJ)
     &                           + NX1AMA(ISYMAI)*(LBJ - 1) + LAI
                           AIBJ  = IT2SQ(ISYMAI,ISYMBJ)
     &                           + NT1AM(ISYMAI)*(BJ - 1) + AI

                           DIFF   = XINT(AIBJ) - XCIA(LAIBJ)
                           CIANM2 = CIANM2 + XCIA(LAIBJ)*XCIA(LAIBJ)
                           DIFF2  = DIFF*DIFF
                           DIFNM2 = DIFNM2 + DIFF2

                           IF (ABS(DIFF) .GT. TOL) THEN
                              NERR = NERR + 1
c                             WRITE(LUPRI,'(1X,A,1P,D20.10)')
c     &                       'Target value       : ',XINT(AIBJ)
c                             WRITE(LUPRI,'(1X,A,1P,D20.10)')
c     &                       'CIA    value       : ',XCIA(LAIBJ)
c                             WRITE(LUPRI,'(1X,A,1P,D20.10,/)')
c     &                       'Absolute difference: ',ABS(DIFF)
                           ENDIF

                           NTST = NTST + 1

                        ENDDO
                     ENDDO

  997                CONTINUE

                  ENDDO

               ENDDO
            ENDDO

  998       CONTINUE

         ENDDO

  999    CONTINUE

      ENDDO

      RETURN
      END
C  /* Deck cc_ciadbgset */
      SUBROUTINE CC_CIADBGSET(IB1,NUMB,IOFB1,LVIRB,NX1AMB,IX1AMB)
C
C     Set index arrays.
C
#include <implicit.h>
      INTEGER IOFB1(8), LVIRB(8), NX1AMB(8), IX1AMB(8,8)
#include <maxorb.h>
#include <ccisvi.h>
#include <ccsdsym.h>
#include <ccorb.h>

      IB2 = IB1 + NUMB - 1

      ISYMB1 = ISVI(IB1)
      ISYMB2 = ISVI(IB2)

      CALL IZERO(IOFB1,NSYM)
      CALL IZERO(LVIRB,NSYM)
      CALL IZERO(NX1AMB,NSYM)
      CALL IZERO(IX1AMB,64)

      DO IB = IB1,IB2
         ISYMB = ISVI(IB)
         LVIRB(ISYMB) = LVIRB(ISYMB) + 1
      ENDDO

      IBEG = IB1
      DO ISYMB = ISYMB1,ISYMB2
         IOFB1(ISYMB) = IBEG + NRHFT - IVIR(ISYMB)
         IBEG = IBEG + LVIRB(ISYMB)
      ENDDO

      DO ISYMBJ = 1,NSYM
         ICOUNT = 0
         DO ISYMJ = 1,NSYM
            ISYMB = MULD2H(ISYMJ,ISYMBJ)
            IX1AMB(ISYMB,ISYMJ) = ICOUNT
            ICOUNT = ICOUNT + LVIRB(ISYMB)*NRHF(ISYMJ)
         ENDDO
         NX1AMB(ISYMBJ) = ICOUNT
      ENDDO

      RETURN
      END
C  /* Deck cc_ciadbganl */
      SUBROUTINE CC_DIADBGANL(XINT,XCIA,ISYMX,
     &                        XINNRM,CIANRM,DIFNRM,DIFRMS,RELNRM,IDIVZ)
C
#include <implicit.h>
#include <ccsdsym.h>
#include <dccsdsym.h>
#include <ccorb.h>

      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0)
      PARAMETER (BIG = 1.0D10)

      IDIVZ = 0

      XINNRM = DSQRT(DDOT(NT2SQ(ISYMX),XINT,1,XINT,1))
      CIANRM = DSQRT(DDOT(NT2SQ(ISYMX),XCIA,1,XCIA,1))

      CALL DSCAL(NT2SQ(ISYMX),XMONE,XCIA,1)
      CALL DAXPY(NT2SQ(ISYMX),ONE,XINT,1,XCIA,1)

      DIFNM2 = DDOT(NT2SQ(ISYMX),XCIA,1,XCIA,1)
      DIFNRM = DSQRT(DIFNM2)
      IF (XT2SQ(1) .GT. ZERO) THEN
         DIFRMS = DSQRT(DIFNM2/XT2SQ(1))
      ELSE
         IDIVZ  = 1
         DIFRMS = BIG
      ENDIF
      IF (XINNRM .GT. ZERO) THEN
         RELNRM = DIFNRM/XINNRM
      ELSE
         IF (IDIVZ .EQ. 1) THEN
            IDIVZ = 12
         ELSE
            IDIVZ = 2
         ENDIF
         RELNRM = BIG
      ENDIF

      RETURN
      END
C  /* Deck cc_ciadbgrc1 */
      SUBROUTINE CC_CIADBGRC1(ISY1,ISY2,LLID,LSMT,LMIO,LLID2)
C
C     Set up full calculation.
C
#include <implicit.h>
#include <maxorb.h>
#include <ccdeco.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ciarc.h>

      LOGICAL LLID, LSMT, LMIO, LLID2

C     Process input.
C     --------------

      ISYCH1 = ISY1
      ISYCH2 = ISY2
      LIDEN  = LLID
      SYMTRZ = LSMT
      CIAMIO = LMIO

C     Set up the rest of ciarc.h
C     --------------------------

      DO ISYMA = 1,NSYM
         IOFA1(ISYMA) = 1
         LVIRA(ISYMA) = NVIR(ISYMA)
         IOFB1(ISYMA) = 1
         LVIRB(ISYMA) = NVIR(ISYMA)
      ENDDO
      DO ISYMAI = 1,NSYM
         NX1AMA(ISYMAI) = NT1AM(ISYMAI)
         NX1AMB(ISYMAI) = NT1AM(ISYMAI)
      ENDDO
      DO ISYMI = 1,NSYM
         DO ISYMA = 1,NSYM
            IX1AMA(ISYMA,ISYMI) = IT1AM(ISYMA,ISYMI)
            IX1AMB(ISYMA,ISYMI) = IT1AM(ISYMA,ISYMI)
         ENDDO
      ENDDO

      DO ISYMBJ = 1,NSYM
         DO ISYMAI = 1,NSYM
            IX2SQ(ISYMAI,ISYMBJ) = IT2SQ(ISYMAI,ISYMBJ)
         ENDDO
      ENDDO
      ISYMX = MULD2H(ISYCH1,ISYCH2)
      NX2SQ = NT2SQ(ISYMX)

      DO ISYCHO = 1,NSYM
         NTOVEC(ISYCHO) = NUMCHO(ISYCHO)
      ENDDO

      RETURN
      END
C  /* Deck cc_ciadbg2 */
      SUBROUTINE CC_CIADBG2(WORK,LWORK,NERR,NCALT)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Test cc_cia for non-tot. sym. integral calculation.
C
C     NERR is *updated* with the number of errors detected.
C
#include <implicit.h>
      DIMENSION WORK(LWORK)
#include <priunit.h>
#include <maxorb.h>
#include <ccisvi.h>
#include <ccdeco.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <dccsdsym.h>
#include <ccorb.h>
#include <ciarc.h>
#include <chomp2.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CC_CIADBG2')

      INTEGER LSYM(8)

      LOGICAL CHMOSV, MP2SSV

      PARAMETER (TINY = 1.0D-14, TINY2 = 1.0D-15)
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (HALF = 0.5D0)

      PARAMETER (IOPTR = 2)

C     Print start message.
C     --------------------

      CALL HEADER('Starting '//SECNAM,-1)
      WRITE(LUPRI,'(5X,A,/)')
     & 'Testing CC_CIA assembly of non-tot. symmetric integrals.'

      XLW = ONE*LWORK
      XMX = ZERO
      DO ISYM = 1,NSYM
         IF (XT2SQ(ISYM) .GT. XMX) XMX = XT2SQ(ISYM)
      ENDDO
      IF (TWO*XMX .GT. XLW) THEN
         WRITE(LUPRI,'(5X,A,/)')
     &   'Oooops! Insufficient memory for this test routine....exiting!'
         RETURN
      ENDIF

C     Set print level in CC_CIA.
C     --------------------------

C     IPRCIA = IPRLVL + 1
      IPRCIA = IPRLVL - 1

C     Turn off calculation of min. mem. requirement.
C     Turn on initialization of integral array in CC_CIA.
C     ---------------------------------------------------

      GETMNM = .FALSE.
      INXINT = .TRUE.

C     Initialize counters.
C     --------------------

      IERR  = 0
      NCALL = 0

C     Initialize LSYM (for opening files).
C     ------------------------------------

      DO ISYM = 1,NSYM
         LSYM(ISYM) = ISYM
      ENDDO

C     Read reduce index array.
C     ------------------------

      KIND1 = 1
      CALL CC_GETIND1(WORK(KIND1),LWORK,LIND1)
      KEND = KIND1 + LIND1
      LWRK = LWORK - KEND + 1

      IF (LWRK .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: index'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need (more than): ',KEND-1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Get the MO vectors L(ia,J).
C     ---------------------------

      KCMO  = KEND
      KEND0 = KCMO  + NLAMDS
      LWRK0 = LWORK - KEND0 + 1

      IF (LWRK0 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,
     &   ' - Allocation: MO'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need     : ',KEND0-1,
     &   'Available: ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

      CALL CHO_RDSIR(DUM1,DUM2,DUM3,WORK(KCMO),WORK(KEND0),LWRK0,
     &               .FALSE.,.FALSE.,.TRUE.)

      CALL CHO_TRFAI(WORK(KCMO),WORK(KEND0),LWRK0)

C     Start loop over 'perturbation' symmetries.
C     ------------------------------------------

      DO ISYMX = 1,NSYM

C        Set up a 'perturbation' matrix of sym. ISYMX
C        by transforming the first AO Cholesky vector
C        of that symmetry.
C        --------------------------------------------

         IF (NUMCHO(ISYMX) .LE. 0) THEN
            WRITE(LUPRI,'(5X,A,A,I1,/,5X,A)')
     &      SECNAM,': NOTICE: No Cholesky vectors in symmetry ',ISYMX,
     &      '- hence no test of CC_CIA in this symmetry!'
            GO TO 999
         ENDIF

         LSCR1 = 0
         DO ISYMD = 1,NSYM
            ISYMG = MULD2H(ISYMD,ISYMX)
            NEED1 = NBAS(ISYMG)*NRHF(ISYMD)
            NEED2 = NBAS(ISYMG)*NVIR(ISYMD)
            LSCR1 = MAX(NEED,NEED1,NEED2)
         ENDDO

         KMATJK = KEND0
         KMATCB = KMATJK + NMATIJ(ISYMX)
         KCHOL  = KMATCB + NMATAB(ISYMX)
         KSCR1  = KCHOL  + N2BST(ISYMX)
         KEND1  = KSCR1  + LSCR1
         LWRK1  = LWORK  - KEND1 + 1

         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Insufficient memory in ',SECNAM,
     &      ' - Allocation: pert. mat.'
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need     : ',KEND1-1,
     &      'Available: ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

         CALL CHO_READN(WORK(KCHOL),1,1,WORK(KIND1),IDUM2,
     &                  ISYMX,IOPTR,WORK(KEND1),LWRK1)

         DO ISYMD = 1,NSYM

            ISYMG = MULD2H(ISYMD,ISYMX)
            ISYMK = ISYMD
            ISYMJ = ISYMG
            ISYMB = ISYMD
            ISYMC = ISYMG

            NTOTG = MAX(NBAS(ISYMG),1)
            NTOTD = MAX(NBAS(ISYMD),1)
            NTOTJ = MAX(NRHF(ISYMJ),1)
            NTOTC = MAX(NVIR(ISYMC),1)

            KOFF1 = KCHOL + IAODIS(ISYMG,ISYMD)

            KOFF2 = KCMO   + IGLMRH(ISYMD,ISYMK)
            KOFF3 = KCMO   + IGLMRH(ISYMG,ISYMJ)
            KOFF4 = KMATJK + IMATIJ(ISYMJ,ISYMK)
            CALL DGEMM('N','N',NBAS(ISYMG),NRHF(ISYMK),NBAS(ISYMD),
     &                 ONE,WORK(KOFF1),NTOTG,WORK(KOFF2),NTOTD,
     &                 ZERO,WORK(KSCR1),NTOTG)
            CALL DGEMM('T','N',NRHF(ISYMJ),NRHF(ISYMK),NBAS(ISYMG),
     &                 ONE,WORK(KOFF3),NTOTG,WORK(KSCR1),NTOTG,
     &                 ZERO,WORK(KOFF4),NTOTJ)

            KOFF2 = KCMO   + IGLMVI(ISYMD,ISYMB)
            KOFF3 = KCMO   + IGLMVI(ISYMG,ISYMC)
            KOFF4 = KMATCB + IMATAB(ISYMC,ISYMB)
            CALL DGEMM('N','N',NBAS(ISYMG),NVIR(ISYMB),NBAS(ISYMD),
     &                 ONE,WORK(KOFF1),NTOTG,WORK(KOFF2),NTOTD,
     &                 ZERO,WORK(KSCR1),NTOTG)
            CALL DGEMM('T','N',NVIR(ISYMC),NVIR(ISYMB),NBAS(ISYMG),
     &                 ONE,WORK(KOFF3),NTOTG,WORK(KSCR1),NTOTG,
     &                 ZERO,WORK(KOFF4),NTOTC)

         ENDDO

C        Reset memory.
C        -------------

         KEND1 = KCHOL
         LWRK1 = LWORK - KEND1 + 1

C        Transform all Cholesky vectors L(ia,J) with this
C        perturbation matrix.
C        ------------------------------------------------

         DO ISYCHO = 1,NSYM

            ISYMAI = MULD2H(ISYCHO,ISYMX)

C           Open Cholesky files.
C           --------------------

            CALL CHO_MOP(-1,1,ISYCHO,LUCHMO,1,1)
            CALL CHO_MOP(-1,7,ISYCHO,LUCHPT,1,ISYMX)

C           Allocation.
C           -----------

            KCHMO = KEND1
            KCHPT = KCHMO + NT1AM(ISYCHO)
            KEND2 = KCHPT + NT1AM(ISYMAI)
            LWRK2 = LWORK - KEND2 + 1

            IF (LWRK2 .LT. 0) THEN
               WRITE(LUPRI,'(5X,A,A)')
     &         'Insufficient memory for transforming L(ai,J) in ',SECNAM
               WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,9X,I1)')
     &         'Need             : ',KEND2-1,
     &         'Available        : ',LWORK,
     &         'Cholesky symmetry: ',ISYCHO
               CALL QUIT('Insufficient memory in '//SECNAM)
            ENDIF

C           Loop through vectors.
C           ---------------------

            DO JVEC = 1,NUMCHO(ISYCHO)

C              Read unperturbed L(ia,J).
C              -------------------------

               CALL CHO_MOREAD(WORK(KCHMO),NT1AM(ISYCHO),1,JVEC,LUCHMO)

C              Calculate perturbed vectors.
C              ----------------------------

               DO ISYMI = 1,NSYM

                  ISYMA = MULD2H(ISYMI,ISYMAI)
                  NTOTA = MAX(NVIR(ISYMA),1)
                  KOFF3 = KCHPT  + IT1AM(ISYMA,ISYMI)

C                 Part 1: virtual transformation.
C                 -------------------------------

                  ISYMB = MULD2H(ISYMI,ISYCHO)

                  KOFF1 = KMATCB + IMATAB(ISYMB,ISYMA)
                  KOFF2 = KCHMO  + IT1AM(ISYMB,ISYMI)

                  NTOTB = MAX(NVIR(ISYMB),1)

                  CALL DGEMM('T','N',
     &                       NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMB),
     &                       ONE,WORK(KOFF1),NTOTB,WORK(KOFF2),NTOTB,
     &                       ZERO,WORK(KOFF3),NTOTA)

C                 Part 2: occupied transformation.
C                 --------------------------------

                  ISYMJ = MULD2H(ISYMI,ISYMX)

                  KOFF1 = KMATJK + IMATIJ(ISYMI,ISYMJ)
                  KOFF2 = KCHMO  + IT1AM(ISYMA,ISYMJ)

                  NTOTI = MAX(NRHF(ISYMI),1)

                  CALL DGEMM('N','T',
     &                       NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMJ),
     &                       XMONE,WORK(KOFF1),NTOTA,WORK(KOFF2),NTOTI,
     &                       ONE,WORK(KOFF3),NTOTA)

               ENDDO

C              Write perturbed vector.
C              -----------------------

               CALL CHO_MOWRITE(WORK(KCHPT),NT1AM(ISYMAI),1,JVEC,LUCHPT)

            ENDDO

C           Close Cholesky files.
C           ---------------------

            CALL CHO_MOP(1,7,ISYCHO,LUCHPT,1,ISYMX)
            CALL CHO_MOP(1,1,ISYCHO,LUCHMO,1,1)

         ENDDO

C        ============================
C        START TESTING CC_CIA MODULE.
C        ============================

C        Calculate target integrals: tot. sym. integrals
C        generated from perturbed Cholesky vectors.
C        -----------------------------------------------

         KXINT = KEND0
         KEND1 = KXINT + NT2SQ(1)
         LWRK1 = LWORK - KEND1 + 1

         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,'(5X,A,A)')
     &      'Insufficient memory in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     &      'Need     : ',KEND1-1,
     &      'Available: ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

         DO ISYCHO = 1,NSYM

            ISYMAI = MULD2H(ISYCHO,ISYMX)

            KCHOL = KEND1
            KEND2 = KCHOL + NT1AM(ISYMAI)*NUMCHO(ISYCHO)
            LWRK2 = LWORK - KEND2 + 1

            IF ((LWRK2.LT.0) .OR. (LWRK2.GT.LWORK)) THEN
               WRITE(LUPRI,'(5X,A,A)')
     &         'Insufficient memory in ',SECNAM
               WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     &         'Need     : ',KEND2-1,
     &         'Available: ',LWORK
               IF (LWRK2 .GT. LWORK) THEN
                  WRITE(LUPRI,'(5X,A)')
     &            '(Calculation seems far too large!)'
               ENDIF
               CALL QUIT('Insufficient memory in '//SECNAM)
            ENDIF

            CALL CHO_MOP(-1,7,ISYCHO,LUCHPT,1,ISYMX)
            CALL CHO_MOREAD(WORK(KCHOL),NT1AM(ISYMAI),NUMCHO(ISYCHO),
     &                      1,LUCHPT)
            CALL CHO_MOP(1,7,ISYCHO,LUCHPT,1,ISYMX)

            NTOAI = MAX(NT1AM(ISYMAI),1)
            KOFFX = KXINT + IT2SQ(ISYMAI,ISYMAI)
            CALL DGEMM('N','T',
     &                 NT1AM(ISYMAI),NT1AM(ISYMAI),NUMCHO(ISYCHO),
     &                 ONE,WORK(KCHOL),NTOAI,WORK(KCHOL),NTOAI,
     &                 ZERO,WORK(KOFFX),NTOAI)

         ENDDO

         DO IDOSCL = 0,1

C           Scale target diagonal.
C           ----------------------

            IF (IDOSCL .EQ. 1) THEN
               DO ISYMAI = 1,NSYM
                  INCRX = NT1AM(ISYMAI) + 1
                  KOFFX = KXINT + IT2SQ(ISYMAI,ISYMAI)
                  CALL DSCAL(NT1AM(ISYMAI),HALF,WORK(KOFFX),INCRX)
               ENDDO
            ENDIF

C           Test 2.1: L1 = L2 option, full calculation.
C           -------------------------------------------

            KXCIA = KEND1
            KEND2 = KXCIA + NT2SQ(1)
            LWRK2 = LWORK - KEND2 + 1

            IF (LWRK2. LT. 0) THEN
               WRITE(LUPRI,'(//,5X,A,A,A)')
     &         'Insufficient memory in ',SECNAM,
     &         ' - Allocation: Test 2.1'
               WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &         'Need     : ',KEND2-1,
     &         'Available: ',LWORK
               CALL QUIT('Insufficient memory in '//SECNAM)
            ENDIF

            DO ICLCIA = 1,2
               NPASS  = 0
               DIAFAC = ONE
               IF (IDOSCL .EQ. 1) DIAFAC = HALF
               CALL CC_CIADBGRC1(ISYMX,ISYMX,.TRUE.,.FALSE.,.FALSE.,
     &                           .FALSE.)
               ITYP1 = 7
               ITYP2 = 7
               IF (ICLCIA .EQ. 2) THEN
                  LW2S = LWRK2
                  NEED = 0
                  NVEC = NUMCHO(1)
                  DO ISYM = 1,NSYM
                     ISYMOV = MULD2H(ISYM,ISYMX)
                     NEEDS  = NT1AM(ISYMOV)
                     NEED   = MAX(NEED,NEEDS)
                     NVEC   = MIN(NVEC,NUMCHO(ISYM))
                  ENDDO
                  MULM  = MIN(NVEC,5)
                  NEED  = MULM*(NEED + 1)
                  LWRK2 = MIN(LWRK2,NEED)
               ENDIF
               CALL CC_CIA(WORK(KXCIA),WORK(KEND2),LWRK2,DIAFAC,KLAST,
     &                     NPASS)
               NCALL = NCALL + 1
               IF (ICLCIA .EQ. 2) THEN
                  LWRK2 = LW2S
               ENDIF

               CALL HEADER('Test 2.1 Result',-1)
               IF (IDOSCL .EQ. 1) THEN
                  WRITE(LUPRI,'(5X,A,1P,D20.10,A)')
     &            '(diagonal scaled with ',DIAFAC,')'
               ENDIF
               IF (ICLCIA .GT. 1) THEN
                  WRITE(LUPRI,'(5X,A)')
     &            '(Cholesky batching was forced in CC_CIA)'
               ENDIF
               CALL CC_DIADBGANL(WORK(KXINT),WORK(KXCIA),1,TRGNRM,
     &                           CIANRM,DIFNRM,DIFRMS,RELNRM,IDIVZ)
               IF ((IDIVZ.EQ.1) .OR. (IDIVZ.EQ.12)) THEN
                  WRITE(LUPRI,'(5X,A)')
     &      '- the difference RMS was assigned to avoid division with 0'
               ENDIF
               IF ((IDIVZ.EQ.2) .OR. (IDIVZ.EQ.12)) THEN
                  WRITE(LUPRI,'(5X,A)')
     &     '- the rel. diff. norm was assigned to avoid division with 0'
               ENDIF
       WRITE(LUPRI,'(5X,A,1P,D20.10,/,5X,A,1P,D20.10,/,5X,A,1P,D20.10)')
     &          'Norm of (ia|jb) [target]: ',TRGNRM,
     &          'Norm of (ia|jb) [cc_cia]: ',CIANRM,
     &          'Difference              : ',TRGNRM-CIANRM
               WRITE(LUPRI,*)
       WRITE(LUPRI,'(5X,A,1P,D20.10,/,5X,A,1P,D20.10,/,5X,A,1P,D20.10)')
     &          'Difference norm         : ',DIFNRM,
     &          'Difference RMS          : ',DIFRMS,
     &          'Relative difference norm: ',RELNRM
               IF ((RELNRM.GT.TINY) .OR. (DIFRMS.GT.TINY)) THEN
                  WRITE(LUPRI,'(5X,A,/)')
     &            '- That counts as an error!'
                  IERR = IERR + 1
               ENDIF
            ENDDO

C           Test 2.4: L1 = L2 option, batch calculation,
C                                     no symmetrization.
C           --------------------------------------------

            DO ICLCIA = 1,2
               DIFNM2 = ZERO
               CIANM2 = ZERO
               NPASS  = 0
               DIAFAC = ONE
               IF (IDOSCL .EQ. 1) DIAFAC = HALF
               ISYCH1 = ISYMX
               ISYCH2 = ISYMX
               ISYINT = MULD2H(ISYCH1,ISYCH2)
               LIDEN  = .TRUE.
               SYMTRZ = .FALSE.
               CIAMIO = .FALSE.
               ITYP1  = 7
               ITYP2  = 7
               DO ISYM = 1,NSYM
                  NTOVEC(ISYM) = NUMCHO(ISYM)
               ENDDO
               IF (NSYM .EQ. 1) THEN
                  NVIRB = MAX(NVIRT/2,1)
                  NVIRA = MAX(NVIRT/3,1)
               ELSE
                  NVIRB = MAX(NVIRT/NSYM,1)
                  NVIRA = MAX(NVIRT/(NSYM+1),1)
               ENDIF
               NBATB = (NVIRT - 1)/NVIRB + 1
               NBATA = (NVIRT - 1)/NVIRA + 1
               DO IBATB = 1,NBATB
                  NUMB = NVIRB
                  IF (IBATB .EQ. NBATB) THEN
                     NUMB = NVIRT - NVIRB*(NBATB - 1)
                  ENDIF
                  IB1 = NVIRB*(IBATB - 1) + 1
                  CALL CC_CIADBGSET(IB1,NUMB,IOFB1,LVIRB,NX1AMB,IX1AMB)
                  DO IBATA = 1,NBATA
                     NUMA = NVIRA
                     IF (IBATA .EQ. NBATA) THEN
                        NUMA = NVIRT - NVIRA*(NBATA - 1)
                     ENDIF
                     IA1 = NVIRA*(IBATA - 1) + 1
                     CALL CC_CIADBGSET(IA1,NUMA,IOFA1,LVIRA,
     &                                 NX1AMA,IX1AMA)
                     CALL IZERO(IX2SQ,64)
                     DO ISYM = 1,NSYM
                        ICOUNT = 0
                        DO ISYMBJ = 1,NSYM
                           ISYMAI = MULD2H(ISYMBJ,ISYM)
                           IX2SQ(ISYMAI,ISYMBJ) = ICOUNT
                           ICOUNT = ICOUNT
     &                            + NX1AMA(ISYMAI)*NX1AMB(ISYMBJ)
                        ENDDO
                        IF (ISYM .EQ. ISYINT) NX2SQ = ICOUNT
                     ENDDO
                     IF (ICLCIA .GT. 1) THEN
                        LW2S = LWRK2
                        NEED = 0
                        NVEC = NUMCHO(1)
                        DO ISYCHO = 1,NSYM
                           NEEDS = NT1AM(ISYCHO)
     &                           + NX1AMA(ISYCHO) + NX1AMB(ISYCHO)
                           NEED  = MAX(NEED,NEEDS)
                           NVEC  = MIN(NVEC,NUMCHO(ISYCHO))
                        ENDDO
                        MULM  = MIN(NVEC,5)
                        NEED  = MULM*NEED
                        LWRK2 = MIN(LWRK2,NEED)
                     ENDIF
                     CALL CC_CIA(WORK(KXCIA),WORK(KEND2),LWRK2,DIAFAC,
     &                           KLAST,NPASS)
                     NCALL = NCALL + 1
                     IF (ICLCIA .GT. 1) THEN
                        LWRK2 = LW2S
                     ENDIF
                     CALL CC_CIADBGANL2(WORK(KXINT),WORK(KXCIA),ISYINT,
     &                                  CIANM2,DIFNM2,NXCIA,JERR,TINY2)
                     CALL HEADER('Test 2.4 Partial Result',-1)
                     IF (IDOSCL .EQ. 1) THEN
                        WRITE(LUPRI,'(5X,A,1P,D20.10,A)')
     &                  '(diagonal scaled with ',DIAFAC,')'
                     ENDIF
                     IA2 = IA1 + NUMA - 1
                     IB2 = IB1 + NUMB - 1
                     IF (ICLCIA .GT. 1) THEN
                        WRITE(LUPRI,'(5X,A)')
     &                  '(Cholesky batching was forced in CC_CIA)'
                     ENDIF
                     WRITE(LUPRI,'(5X,A,I4,A,I4,/,5X,A,I4,A,I4)')
     &               'This is a-batch number',IBATA,' of',NBATA,
     &               '--- for b-batch number',IBATB,' of',NBATB
                     WRITE(LUPRI,'(5X,A,I4,1X,I1,A,I4,1X,I1)')
     &               'a-index [a,sym.a] runs from',IA1,ISVI(IA1),' to',
     &               IA2,ISVI(IA2)
                     WRITE(LUPRI,'(5X,A,I4,1X,I1,A,I4,1X,I1)')
     &               'b-index [b,sym.b] runs from',IB1,ISVI(IB1),' to',
     &               IB2,ISVI(IB2)
                     WRITE(LUPRI,'(5X,A,I10)')
     &               'Discrepancies in this batch: ',JERR
                     IF (JERR .NE. 0) THEN
                        IERR = IERR + 1
                        WRITE(LUPRI,'(5X,A)')
     &                  '- That counts as an error!'
                     ENDIF
                     IF (NXCIA .NE. NX2SQ) THEN
                        WRITE(LUPRI,'(5X,A)')
     &                  'Something is terribly wrong: NXCIA != NX2SQ'
                        WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     &                  'NXCIA = ',NXCIA,
     &                  'NX2SQ = ',NX2SQ
                        WRITE(LUPRI,'(5X,A)')
     &                  '- That counts as an error!'
                        IERR = IERR + 1
                     ENDIF
                  ENDDO
               ENDDO
               CALL HEADER('Test 2.4 Global Result',-1)
               IF (IDOSCL .EQ. 1) THEN
                  WRITE(LUPRI,'(5X,A,1P,D20.10,A)')
     &            '(diagonal scaled with ',DIAFAC,')'
               ENDIF
               IF (ICLCIA .GT. 1) THEN
                  WRITE(LUPRI,'(5X,A)')
     &            '(Cholesky batching was forced in CC_CIA)'
               ENDIF
         TRGNRM = DSQRT(DDOT(NT2SQ(ISYINT),WORK(KXINT),1,WORK(KXINT),1))
               CIANRM = DSQRT(CIANM2)
               DIFNRM = DSQRT(DIFNM2)
               DIFRMS = DSQRT(DIFNM2/XT2SQ(ISYINT))
               RELNRM = DIFNRM/TRGNRM
       WRITE(LUPRI,'(5X,A,1P,D20.10,/,5X,A,1P,D20.10,/,5X,A,1P,D20.10)')
     &         'Norm of (ia|jb) [target]: ',TRGNRM,
     &         'Norm of (ia|jb) [cc_cia]: ',CIANRM,
     &         'Difference              : ',TRGNRM-CIANRM
               WRITE(LUPRI,*)
       WRITE(LUPRI,'(5X,A,1P,D20.10,/,5X,A,1P,D20.10,/,5X,A,1P,D20.10)')
     &         'Difference norm         : ',DIFNRM,
     &         'Difference RMS          : ',DIFRMS,
     &         'Relative difference norm: ',RELNRM
               IF ((RELNRM.GT.TINY) .OR. (DIFRMS.GT.TINY)) THEN
                  WRITE(LUPRI,'(5X,A,/)')
     &            '- That counts as an error!'
                  IERR = IERR + 1
               ENDIF
            ENDDO

         ENDDO

C        Calculate new non-tot. sym. target integrals
C        of type L1 * L2, with L1 != L2. Integrals
C        are first calculated non-symmetrically.
C        --------------------------------------------

         DIAFAC = ONE

         ISYCH1 = 1
         ISYCH2 = ISYMX
         ISYINT = MULD2H(ISYCH1,ISYCH2)

         KXINT = KEND0
         KEND1 = KXINT + NT2SQ(ISYINT)
         LWRK1 = LWORK - KEND1 + 1

         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,'(5X,A,A)')
     &      'Insufficient memory in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     &      'Need     : ',KEND1-1,
     &      'Available: ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

         DO ISYCHO = 1,NSYM

            ISYMAI = MULD2H(ISYCHO,ISYCH1)
            ISYMBJ = MULD2H(ISYCHO,ISYCH2)

            KCHO1 = KEND1
            KCHO2 = KCHO1 + NT1AM(ISYMAI)*NUMCHO(ISYCHO)
            KEND2 = KCHO2 + NT1AM(ISYMBJ)*NUMCHO(ISYCHO)
            LWRK2 = LWORK - KEND2 + 1

            IF ((LWRK2.LT.0) .OR. (LWRK2.GT.LWORK)) THEN
               WRITE(LUPRI,'(5X,A,A)')
     &         'Insufficient memory in ',SECNAM
               WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     &         'Need     : ',KEND2-1,
     &         'Available: ',LWORK
               IF (LWRK2 .GT. LWORK) THEN
                  WRITE(LUPRI,'(5X,A)')
     &            '(Calculation seems far too large!)'
               ENDIF
               CALL QUIT('Insufficient memory in '//SECNAM)
            ENDIF

            CALL CHO_MOP(-1,1,ISYCHO,LUCHAI,1,ISYCH1)
            CALL CHO_MOP(-1,7,ISYCHO,LUCHBJ,1,ISYCH2)
            CALL CHO_MOREAD(WORK(KCHO1),NT1AM(ISYMAI),NUMCHO(ISYCHO),
     &                      1,LUCHAI)
            CALL CHO_MOREAD(WORK(KCHO2),NT1AM(ISYMBJ),NUMCHO(ISYCHO),
     &                      1,LUCHBJ)
            CALL CHO_MOP(1,7,ISYCHO,LUCHBJ,1,ISYCH2)
            CALL CHO_MOP(1,7,ISYCHO,LUCHAI,1,ISYCH1)

            NTOAI = MAX(NT1AM(ISYMAI),1)
            NTOBJ = MAX(NT1AM(ISYMBJ),1)
            KOFFX = KXINT + IT2SQ(ISYMAI,ISYMBJ)
            CALL DGEMM('N','T',
     &                 NT1AM(ISYMAI),NT1AM(ISYMBJ),NUMCHO(ISYCHO),
     &                 ONE,WORK(KCHO1),NTOAI,WORK(KCHO2),NTOBJ,
     &                 ZERO,WORK(KOFFX),NTOAI)

         ENDDO

         DO ISMTZ = 1,2

C           Symmetrize target integrals.
C           ----------------------------

            IF (ISMTZ .EQ. 2) THEN
               DO ISYMBJ = 1,NSYM
                  ISYMAI = MULD2H(ISYMBJ,ISYINT)
                  IF (ISYMAI .EQ. ISYMBJ) THEN
                     DO LBJ = 1,NT1AM(ISYMBJ)
                        DO LAI = 1,LBJ
                           LAIBJ = KXINT + IT2SQ(ISYMAI,ISYMBJ)
     &                           + NT1AM(ISYMAI)*(LBJ - 1) + LAI - 1
                           LBJAI = KXINT + IT2SQ(ISYMBJ,ISYMAI)
     &                           + NT1AM(ISYMBJ)*(LAI - 1) + LBJ - 1
                           WORK(LAIBJ) = WORK(LAIBJ) + WORK(LBJAI)
                           WORK(LBJAI) = WORK(LAIBJ)
                        ENDDO
                     ENDDO
                  ELSE IF (ISYMAI .LT. ISYMBJ) THEN
                     DO LBJ = 1,NT1AM(ISYMBJ)
                        DO LAI = 1,NT1AM(ISYMAI)
                           LAIBJ = KXINT + IT2SQ(ISYMAI,ISYMBJ)
     &                           + NT1AM(ISYMAI)*(LBJ - 1) + LAI - 1
                           LBJAI = KXINT + IT2SQ(ISYMBJ,ISYMAI)
     &                           + NT1AM(ISYMBJ)*(LAI - 1) + LBJ - 1
                           WORK(LAIBJ) = WORK(LAIBJ) + WORK(LBJAI)
                           WORK(LBJAI) = WORK(LAIBJ)
                        ENDDO
                     ENDDO
                  ENDIF
               ENDDO
            ENDIF

C           Test 2.8: L1 != L2 option, full calculation.
C           --------------------------------------------

            KXCIA = KEND1
            KEND2 = KXCIA + NT2SQ(ISYINT)
            LWRK2 = LWORK - KEND2 + 1

            IF (LWRK2. LT. 0) THEN
               WRITE(LUPRI,'(//,5X,A,A,A)')
     &         'Insufficient memory in ',SECNAM,
     &         ' - Allocation: Test 2.8'
               WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &         'Need     : ',KEND2-1,
     &         'Available: ',LWORK
               CALL QUIT('Insufficient memory in '//SECNAM)
            ENDIF

            DO ICLCIA = 1,2
               NPASS = 0
               IF (ISMTZ .EQ. 2) THEN
                  CALL CC_CIADBGRC1(ISYCH1,ISYCH2,.FALSE.,.TRUE.,
     &                              .FALSE.,.FALSE.)
               ELSE
                  CALL CC_CIADBGRC1(ISYCH1,ISYCH2,.FALSE.,.FALSE.,
     &                              .FALSE.,.FALSE.)
               ENDIF
               ITYP1 = 1
               ITYP2 = 7
               IF (ICLCIA .EQ. 2) THEN
                  LW2S = LWRK2
                  NEED = 0
                  NVEC = NUMCHO(1)
                  DO ISYM = 1,NSYM
                     ISYMAI = MULD2H(ISYM,ISYCH1)
                     ISYMBJ = MULD2H(ISYM,ISYCH2)
                     NEEDS  = NT1AM(ISYMAI) + NT1AM(ISYMBJ)
                     NEED   = MAX(NEED,NEEDS)
                     NVEC   = MIN(NVEC,NUMCHO(ISYM))
                  ENDDO
                  MULM  = MIN(NVEC,5)
                  NEED  = MULM*(NEED + 1)
                  LWRK2 = MIN(LWRK2,NEED)
               ENDIF
               CALL CC_CIA(WORK(KXCIA),WORK(KEND2),LWRK2,DIAFAC,KLAST,
     &                     NPASS)
               NCALL = NCALL + 1
               IF (ICLCIA .EQ. 2) THEN
                  LWRK2 = LW2S
               ENDIF

               CALL HEADER('Test 2.8 Result',-1)
               IF (ISMTZ .EQ. 2) THEN
                  WRITE(LUPRI,'(5X,A)')
     &            '(Integrals symmetrized)'
               ENDIF
               IF (ICLCIA .GT. 1) THEN
                  WRITE(LUPRI,'(5X,A)')
     &            '(Cholesky batching was forced in CC_CIA)'
               ENDIF
               CALL CC_DIADBGANL(WORK(KXINT),WORK(KXCIA),ISYINT,TRGNRM,
     &                           CIANRM,DIFNRM,DIFRMS,RELNRM,IDIVZ)
               IF ((IDIVZ.EQ.1) .OR. (IDIVZ.EQ.12)) THEN
                  WRITE(LUPRI,'(5X,A)')
     &      '- the difference RMS was assigned to avoid division with 0'
               ENDIF
               IF ((IDIVZ.EQ.2) .OR. (IDIVZ.EQ.12)) THEN
                  WRITE(LUPRI,'(5X,A)')
     &     '- the rel. diff. norm was assigned to avoid division with 0'
               ENDIF
       WRITE(LUPRI,'(5X,A,1P,D20.10,/,5X,A,1P,D20.10,/,5X,A,1P,D20.10)')
     &          'Norm of (ia|jb) [target]: ',TRGNRM,
     &          'Norm of (ia|jb) [cc_cia]: ',CIANRM,
     &          'Difference              : ',TRGNRM-CIANRM
               WRITE(LUPRI,*)
       WRITE(LUPRI,'(5X,A,1P,D20.10,/,5X,A,1P,D20.10,/,5X,A,1P,D20.10)')
     &          'Difference norm         : ',DIFNRM,
     &          'Difference RMS          : ',DIFRMS,
     &          'Relative difference norm: ',RELNRM
               IF ((RELNRM.GT.TINY) .OR. (DIFRMS.GT.TINY)) THEN
                  WRITE(LUPRI,'(5X,A,/)')
     &            '- That counts as an error!'
                  IERR = IERR + 1
               ENDIF
            ENDDO

C           Test 2.9: L1 != L2 option, full calculation.
C                                      L1 <-> L2 swap,
C                     (I.e. test only for ISMTZ = 2)
C           --------------------------------------------

            IF (ISMTZ .EQ. 2) THEN

               DO ICLCIA = 1,2
                  NPASS = 0
                  ISY1SV = ISYCH1
                  ISY2SV = ISYCH2
                  CALL CC_CIADBGRC1(ISY2SV,ISY1SV,.FALSE.,.TRUE.,
     &                              .FALSE.,.FALSE.)
                  ITYP1 = 7
                  ITYP2 = 1
                  IF (ICLCIA .EQ. 2) THEN
                     LW2S = LWRK2
                     NEED = 0
                     NVEC = NUMCHO(1)
                     DO ISYM = 1,NSYM
                        ISYMAI = MULD2H(ISYM,ISY2SV)
                        ISYMBJ = MULD2H(ISYM,ISY1SV)
                        NEEDS  = NT1AM(ISYMAI) + NT1AM(ISYMBJ)
                        NEED   = MAX(NEED,NEEDS)
                        NVEC   = MIN(NVEC,NUMCHO(ISYM))
                     ENDDO
                     MULM  = MIN(NVEC,5)
                     NEED  = MULM*(NEED + 1)
                     LWRK2 = MIN(LWRK2,NEED)
                  ENDIF
                  CALL CC_CIA(WORK(KXCIA),WORK(KEND2),LWRK2,DIAFAC,
     &                        KLAST,NPASS)
                  NCALL = NCALL + 1
                  IF (ICLCIA .EQ. 2) THEN
                     LWRK2 = LW2S
                  ENDIF
                  ISYCH1 = ISY1SV
                  ISYCH2 = ISY2SV

                  CALL HEADER('Test 2.9 Result',-1)
                  WRITE(LUPRI,'(5X,A)')
     &            '(Integrals symmetrized, L2 <-> L1 swap)'
                  IF (ICLCIA .EQ. 2) THEN
                     WRITE(LUPRI,'(5X,A)')
     &               '(Cholesky batching was forced in CC_CIA)'
                  ENDIF
                  CALL CC_DIADBGANL(WORK(KXINT),WORK(KXCIA),ISYINT,
     &                              TRGNRM,CIANRM,DIFNRM,DIFRMS,RELNRM,
     &                              IDIVZ)
                  IF ((IDIVZ.EQ.1) .OR. (IDIVZ.EQ.12)) THEN
                     WRITE(LUPRI,'(5X,A)')
     &      '- the difference RMS was assigned to avoid division with 0'
                  ENDIF
                  IF ((IDIVZ.EQ.2) .OR. (IDIVZ.EQ.12)) THEN
                     WRITE(LUPRI,'(5X,A)')
     &     '- the rel. diff. norm was assigned to avoid division with 0'
                  ENDIF
       WRITE(LUPRI,'(5X,A,1P,D20.10,/,5X,A,1P,D20.10,/,5X,A,1P,D20.10)')
     &            'Norm of (ia|jb) [target]: ',TRGNRM,
     &            'Norm of (ia|jb) [cc_cia]: ',CIANRM,
     &            'Difference              : ',TRGNRM-CIANRM
                  WRITE(LUPRI,*)
       WRITE(LUPRI,'(5X,A,1P,D20.10,/,5X,A,1P,D20.10,/,5X,A,1P,D20.10)')
     &            'Difference norm         : ',DIFNRM,
     &            'Difference RMS          : ',DIFRMS,
     &            'Relative difference norm: ',RELNRM
                  IF ((RELNRM.GT.TINY) .OR. (DIFRMS.GT.TINY)) THEN
                     WRITE(LUPRI,'(5X,A,/)')
     &               '- That counts as an error!'
                     IERR = IERR + 1
                  ENDIF
               ENDDO

            ENDIF

C           Test 2.10: L1 != L2 option, batch calculation.
C           ----------------------------------------------

            DO ICLCIA = 1,2
               DIFNM2 = ZERO
               CIANM2 = ZERO
               NPASS  = 0
               LIDEN  = .FALSE.
               CIAMIO = .FALSE.
               SYMTRZ = .FALSE.
               IF (ISMTZ .EQ. 2) SYMTRZ = .TRUE.
               DO ISYM = 1,NSYM
                  NTOVEC(ISYM) = NUMCHO(ISYM)
               ENDDO
               ITYP1 = 1
               ITYP2 = 7
               IF (NSYM .EQ. 1) THEN
                  NVIRB = MAX(NVIRT/2,1)
                  NVIRA = MAX(NVIRT/3,1)
               ELSE     
                  NVIRB = MAX(NVIRT/NSYM,1)
                  NVIRA = MAX(NVIRT/(NSYM+1),1)
               ENDIF
               NBATB = (NVIRT - 1)/NVIRB + 1
               NBATA = (NVIRT - 1)/NVIRA + 1
               DO IBATB = 1,NBATB
                  NUMB = NVIRB
                  IF (IBATB .EQ. NBATB) THEN
                     NUMB = NVIRT - NVIRB*(NBATB - 1)
                  ENDIF
                  IB1 = NVIRB*(IBATB - 1) + 1
                  CALL CC_CIADBGSET(IB1,NUMB,IOFB1,LVIRB,NX1AMB,IX1AMB)
                  DO IBATA = 1,NBATA
                     NUMA = NVIRA
                     IF (IBATA .EQ. NBATA) THEN
                        NUMA = NVIRT - NVIRA*(NBATA - 1)
                     ENDIF
                     IA1 = NVIRA*(IBATA - 1) + 1
                     CALL CC_CIADBGSET(IA1,NUMA,IOFA1,LVIRA,
     &                                 NX1AMA,IX1AMA)
                     CALL IZERO(IX2SQ,64)
                     DO ISYM = 1,NSYM
                        ICOUNT = 0
                        DO ISYMBJ = 1,NSYM
                           ISYMAI = MULD2H(ISYMBJ,ISYM)
                           IX2SQ(ISYMAI,ISYMBJ) = ICOUNT
                           ICOUNT = ICOUNT
     &                            + NX1AMA(ISYMAI)*NX1AMB(ISYMBJ)
                        ENDDO
                        IF (ISYM .EQ. ISYINT) NX2SQ = ICOUNT
                     ENDDO
                     IF (ICLCIA .GT. 1) THEN
                        LW2S = LWRK2
                        NEED = 0
                        NVEC = NUMCHO(1)
                        DO ISYCHO = 1,NSYM
                           ISYMAI = MULD2H(ISYCHO,ISYCH1)
                           ISYMBJ = MULD2H(ISYCHO,ISYCH2)
                           NEEDS = NT1AM(ISYMAI) + NT1AM(ISYMBJ)
     &                           + NX1AMA(ISYMAI) + NX1AMB(ISYMBJ)
                           NEED  = MAX(NEED,NEEDS)
                           NVEC  = MIN(NVEC,NUMCHO(ISYCHO))
                        ENDDO
                        MULM  = MIN(NVEC,5)
                        NEED  = MULM*NEED
                        LWRK2 = MIN(LWRK2,NEED)
                     ENDIF
                     CALL CC_CIA(WORK(KXCIA),WORK(KEND2),LWRK2,DIAFAC,
     &                           KLAST,NPASS)
                     NCALL = NCALL + 1
                     IF (ICLCIA .GT. 1) THEN
                        LWRK2 = LW2S
                     ENDIF
                     CALL CC_CIADBGANL2(WORK(KXINT),WORK(KXCIA),ISYINT,
     &                                  CIANM2,DIFNM2,NXCIA,JERR,TINY2)
                     CALL HEADER('Test 2.10 Partial Result',-1)
                     IA2 = IA1 + NUMA - 1
                     IB2 = IB1 + NUMB - 1
                     IF (ISMTZ .EQ. 2) THEN
                        WRITE(LUPRI,'(5X,A)')
     &                  '(Integrals symmetrized)'
                     ENDIF
                     IF (ICLCIA .GT. 1) THEN
                        WRITE(LUPRI,'(5X,A)')
     &                  '(Cholesky batching was forced in CC_CIA)'
                     ENDIF
                     WRITE(LUPRI,'(5X,A,I4,A,I4,/,5X,A,I4,A,I4)')
     &               'This is a-batch number',IBATA,' of',NBATA,
     &               '--- for b-batch number',IBATB,' of',NBATB
                     WRITE(LUPRI,'(5X,A,I4,1X,I1,A,I4,1X,I1)')
     &               'a-index [a,sym.a] runs from',IA1,ISVI(IA1),' to',
     &               IA2,ISVI(IA2)
                     WRITE(LUPRI,'(5X,A,I4,1X,I1,A,I4,1X,I1)')
     &               'b-index [b,sym.b] runs from',IB1,ISVI(IB1),' to',
     &               IB2,ISVI(IB2)
                     WRITE(LUPRI,'(5X,A,I10)')
     &               'Discrepancies in this batch: ',JERR
                     IF (JERR .NE. 0) THEN
                        IERR = IERR + 1
                        WRITE(LUPRI,'(5X,A)')
     &                  '- That counts as an error!'
                     ENDIF
                     IF (NXCIA .NE. NX2SQ) THEN
                        WRITE(LUPRI,'(5X,A)')
     &                  'Something is terribly wrong: NXCIA != NX2SQ'
                        WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     &                  'NXCIA = ',NXCIA,
     &                  'NX2SQ = ',NX2SQ
                        WRITE(LUPRI,'(5X,A)')
     &                  '- That counts as an error!'
                        IERR = IERR + 1
                     ENDIF
                  ENDDO
               ENDDO
               CALL HEADER('Test 2.10 Global Result',-1)
               IF (ISMTZ .EQ. 2) THEN
                  WRITE(LUPRI,'(5X,A)')
     &            '(Integrals symmetrized)'
               ENDIF
               IF (ICLCIA .GT. 1) THEN
                  WRITE(LUPRI,'(5X,A)')
     &            '(Cholesky batching was forced in CC_CIA)'
               ENDIF
         TRGNRM = DSQRT(DDOT(NT2SQ(ISYINT),WORK(KXINT),1,WORK(KXINT),1))
               CIANRM = DSQRT(CIANM2)
               DIFNRM = DSQRT(DIFNM2)
               DIFRMS = DSQRT(DIFNM2/XT2SQ(ISYINT))
               RELNRM = DIFNRM/TRGNRM
       WRITE(LUPRI,'(5X,A,1P,D20.10,/,5X,A,1P,D20.10,/,5X,A,1P,D20.10)')
     &         'Norm of (ia|jb) [target]: ',TRGNRM,
     &         'Norm of (ia|jb) [cc_cia]: ',CIANRM,
     &         'Difference              : ',TRGNRM-CIANRM
               WRITE(LUPRI,*)
       WRITE(LUPRI,'(5X,A,1P,D20.10,/,5X,A,1P,D20.10,/,5X,A,1P,D20.10)')
     &         'Difference norm         : ',DIFNRM,
     &         'Difference RMS          : ',DIFRMS,
     &         'Relative difference norm: ',RELNRM
               IF ((RELNRM.GT.TINY) .OR. (DIFRMS.GT.TINY)) THEN
                  WRITE(LUPRI,'(5X,A,/)')
     &            '- That counts as an error!'
                  IERR = IERR + 1
               ENDIF
            ENDDO

C           Test 2.11: As 2.9, but with min. I/O symmetrization.
C                      (I.e. test only for ISMTZ = 2)
C           ----------------------------------------------------

            IF (ISMTZ .EQ. 2) THEN

               SYMTRZ = .TRUE.

               DO ICLCIA = 1,2
                  DIFNM2 = ZERO
                  CIANM2 = ZERO
                  NPASS  = 0
                  LIDEN  = .FALSE.
                  CIAMIO = .TRUE.
                  DO ISYM = 1,NSYM
                     NTOVEC(ISYM) = NUMCHO(ISYM)
                  ENDDO
                  ITYP1 = 1
                  ITYP2 = 7
                  IF (NSYM .EQ. 1) THEN
                     NVIRB = MAX(NVIRT/2,1)
                     NVIRA = MAX(NVIRT/3,1)
                  ELSE
                     NVIRB = MAX(NVIRT/NSYM,1)
                     NVIRA = MAX(NVIRT/(NSYM+1),1)
                  ENDIF
                  NBATB = (NVIRT - 1)/NVIRB + 1
                  NBATA = (NVIRT - 1)/NVIRA + 1
                  DO IBATB = 1,NBATB
                     NUMB = NVIRB
                     IF (IBATB .EQ. NBATB) THEN
                        NUMB = NVIRT - NVIRB*(NBATB - 1)
                     ENDIF
                     IB1 = NVIRB*(IBATB - 1) + 1
                     CALL CC_CIADBGSET(IB1,NUMB,IOFB1,LVIRB,
     &                                 NX1AMB,IX1AMB)
                     DO IBATA = 1,NBATA
                        NUMA = NVIRA
                        IF (IBATA .EQ. NBATA) THEN
                           NUMA = NVIRT - NVIRA*(NBATA - 1)
                        ENDIF
                        IA1 = NVIRA*(IBATA - 1) + 1
                        CALL CC_CIADBGSET(IA1,NUMA,IOFA1,LVIRA,
     &                                    NX1AMA,IX1AMA)
                        CALL IZERO(IX2SQ,64)
                        DO ISYM = 1,NSYM
                           ICOUNT = 0
                           DO ISYMBJ = 1,NSYM
                              ISYMAI = MULD2H(ISYMBJ,ISYM)
                              IX2SQ(ISYMAI,ISYMBJ) = ICOUNT
                              ICOUNT = ICOUNT
     &                               + NX1AMA(ISYMAI)*NX1AMB(ISYMBJ)
                           ENDDO
                           IF (ISYM .EQ. ISYINT) NX2SQ = ICOUNT
                        ENDDO
                        IF (ICLCIA .GT. 1) THEN
                           LW2S = LWRK2
                           NEED = 0
                           NVEC = NUMCHO(1)
                           DO ISYCHO = 1,NSYM
                              ISYMAI = MULD2H(ISYCHO,ISYCH1)
                              ISYMBJ = MULD2H(ISYCHO,ISYCH2)
                              NEEDS = NT1AM(ISYMAI) + NT1AM(ISYMBJ)
     &                              + NX1AMA(ISYMAI) + NX1AMB(ISYMBJ)
                              NEED  = MAX(NEED,NEEDS)
                              NVEC  = MIN(NVEC,NUMCHO(ISYCHO))
                           ENDDO
                           MULM  = MIN(NVEC,5)
                           NEED  = MULM*NEED
                           LWRK2 = MIN(LWRK2,NEED)
                        ENDIF
                        CALL CC_CIA(WORK(KXCIA),WORK(KEND2),LWRK2,
     &                              DIAFAC,KLAST,NPASS)
                        NCALL = NCALL + 1
                        IF (ICLCIA .GT. 1) THEN
                           LWRK2 = LW2S
                        ENDIF
                        CALL CC_CIADBGANL2(WORK(KXINT),WORK(KXCIA),
     &                                     ISYINT,CIANM2,DIFNM2,NXCIA,
     &                                     JERR,TINY2)
                        CALL HEADER('Test 2.11 Partial Result',-1)
                        IA2 = IA1 + NUMA - 1
                        IB2 = IB1 + NUMB - 1
                        WRITE(LUPRI,'(5X,A)')
     &                  '(Integrals symmetrized with minimum I/O)'
                        IF (ICLCIA .GT. 1) THEN
                           WRITE(LUPRI,'(5X,A)')
     &                     '(Cholesky batching was forced in CC_CIA)'
                        ENDIF
                        WRITE(LUPRI,'(5X,A,I4,A,I4,/,5X,A,I4,A,I4)')
     &                  'This is a-batch number',IBATA,' of',NBATA,
     &                  '--- for b-batch number',IBATB,' of',NBATB
                        WRITE(LUPRI,'(5X,A,I4,1X,I1,A,I4,1X,I1)')
     &                'a-index [a,sym.a] runs from',IA1,ISVI(IA1),' to',
     &                  IA2,ISVI(IA2)
                        WRITE(LUPRI,'(5X,A,I4,1X,I1,A,I4,1X,I1)')
     &                'b-index [b,sym.b] runs from',IB1,ISVI(IB1),' to',
     &                  IB2,ISVI(IB2)
                        WRITE(LUPRI,'(5X,A,I10)')
     &                  'Discrepancies in this batch: ',JERR
                        IF (JERR .NE. 0) THEN
                           IERR = IERR + 1
                           WRITE(LUPRI,'(5X,A)')
     &                     '- That counts as an error!'
                        ENDIF
                        IF (NXCIA .NE. NX2SQ) THEN
                           WRITE(LUPRI,'(5X,A)')
     &                     'Something is terribly wrong: NXCIA != NX2SQ'
                           WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     &                     'NXCIA = ',NXCIA,
     &                     'NX2SQ = ',NX2SQ
                           WRITE(LUPRI,'(5X,A)')
     &                     '- That counts as an error!'
                           IERR = IERR + 1
                        ENDIF
                     ENDDO
                  ENDDO
                  CALL HEADER('Test 2.11 Global Result',-1)
                  WRITE(LUPRI,'(5X,A)')
     &            '(Integrals symmetrized with minimum I/O)'
                  IF (ICLCIA .GT. 1) THEN
                     WRITE(LUPRI,'(5X,A)')
     &               '(Cholesky batching was forced in CC_CIA)'
                  ENDIF
         TRGNRM = DSQRT(DDOT(NT2SQ(ISYINT),WORK(KXINT),1,WORK(KXINT),1))
                  CIANRM = DSQRT(CIANM2)
                  DIFNRM = DSQRT(DIFNM2)
                  DIFRMS = DSQRT(DIFNM2/XT2SQ(ISYINT))
                  RELNRM = DIFNRM/TRGNRM
       WRITE(LUPRI,'(5X,A,1P,D20.10,/,5X,A,1P,D20.10,/,5X,A,1P,D20.10)')
     &            'Norm of (ia|jb) [target]: ',TRGNRM,
     &            'Norm of (ia|jb) [cc_cia]: ',CIANRM,
     &            'Difference              : ',TRGNRM-CIANRM
                  WRITE(LUPRI,*)
       WRITE(LUPRI,'(5X,A,1P,D20.10,/,5X,A,1P,D20.10,/,5X,A,1P,D20.10)')
     &            'Difference norm         : ',DIFNRM,
     &            'Difference RMS          : ',DIFRMS,
     &            'Relative difference norm: ',RELNRM
                  IF ((RELNRM.GT.TINY) .OR. (DIFRMS.GT.TINY)) THEN
                     WRITE(LUPRI,'(5X,A,/)')
     &               '- That counts as an error!'
                     IERR = IERR + 1
                  ENDIF
               ENDDO

            ENDIF

         ENDDO

C        Escape point when no vectors of sym. ISYMX.
C        -------------------------------------------

  999    CONTINUE

      ENDDO

C     Update global error and call counters.
C     --------------------------------------

      NERR  = NERR + IERR
      NCALT = NCALT + NCALL

C     Print end message.
C     ------------------

      CALL HEADER('End of '//SECNAM,-1)
      WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     & 'Total number of calls to CC_CIA in this routine:',NCALL,
     & 'Accumulated number of calls                    :',NCALT
      WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     & 'Number of errors detected for non-tot. sym. assembly:',IERR,
     & 'Accumulated number of errors                        :',NERR

      RETURN
      END
