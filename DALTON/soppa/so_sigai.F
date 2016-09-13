C  /* Deck so_sigai */
      SUBROUTINE SO_SIGAI(SIGAI1,LSIGAI1,SIGAI2,LSIGAI2,T2M1,LT2M1,
     &                    XINT,LXINT,BTR1E,LBTR1E,BTR1D,LBTR1D,
     &                    BTJ1E,LBTJ1E,BTJ1D,LBTJ1D,CMO,LCMO,
     &                    ISYMDEL,ISYDIS,ISYMAI,ISYMTR,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak and Henrik Koch, September 1995
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C     Rasmus Faber 2016 : Reordered code, which should result in
C                         improved memory access.
C
C     PURPOSE: Calculate SIGMA1(ALFA,I) and SIGMA2(ALFA,I)
C                                     ~ ~
C              which are defined as ( k c | alfa delta ) x
C              T2M1( c k i delta) with summation of repeated indices.
C              See also eq. (34) and (35).
C
#include "implicit.h"
#include "priunit.h"
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION SIGAI1(LSIGAI1), SIGAI2(LSIGAI2), T2M1(LT2M1)
      DIMENSION XINT(LXINT),     BTR1E(LBTR1E),   BTR1D(LBTR1D)
      DIMENSION BTJ1E(LBTJ1E),   BTJ1D(LBTJ1D),   CMO(LCMO)
      DIMENSION WORK(LWORK)
C
      INTEGER   ALFA
C
#include "ccorb.h"
#include "soppinf.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_SIGAI')
C
C-----------------------------------------------------------
C     Test if T2M1 contains any elements. If not there is no
C     contribution.
C-----------------------------------------------------------
C
      IF (LT2M1 .EQ. 0) THEN
         CALL QEXIT('SO_SIGAI')
         RETURN
      END IF
C
      DO 100 ISALFA = 1,NSYM
C
         ISBEGA = MULD2H(ISALFA,ISYDIS)
         ISYMCK = MULD2H(ISBEGA,ISYMTR)
         ISYMI  = MULD2H(ISALFA,ISYMAI)
C
         LSCR1  = N2BST(ISBEGA)
         LSCR2 = NT1AM(ISYMCK)
         LSCR3 = LSCR2
C
         KSCR1 = 1
         KSCR2 = KSCR1 + LSCR1
         KSCR3 = KSCR2 + LSCR2
         KEND1 = KSCR3 + LSCR3
         LWORK1  = LWORK - KEND1
         KEND1B = KEND1
C
         CALL SO_MEMMAX ('SO_SIGAI.1',LWORK1)
         IF (LWORK1 .LT. 0) CALL STOPIT('SO_SIGAI.1',' ',KEND1,LWORK)
C
         DO 110 ALFA = 1,NBAS(ISALFA)
C
            KOFF1 = IDSAOG(ISALFA,ISYDIS) + NNBST(ISBEGA)*(ALFA-1) + 1


C
C---------------------------------------------------------------------
C           Get a squared set of ( beta gamma | alfa delta ) for given
C           alfa and delta.
C---------------------------------------------------------------------
C
            CALL CCSD_SYMSQ(XINT(KOFF1),ISBEGA,WORK(KSCR1))
C--------------------------------------------------------
C           for each beta gamma symmetry block, generate
C             ( c beta | .. ) * x ( beta, k)
C           - ( beta k | .. ) * x ( c, beta )
C--------------------------------------------------------
            DO 120 ISYMK = 1,NSYM
C
               ISBETA = MULD2H(ISYMK,ISYMTR)
               ISGAM  = MULD2H(ISBETA,ISBEGA)
               ISYMC  = ISGAM
C
               LSCR4  = NBAS(ISGAM)*NRHF(ISYMK)
C
               KPOS2 = KSCR2 + IT1AM(ISYMC,ISYMK)
               KPOS3 = KSCR3 + IT1AM(ISYMC,ISYMK)

               KSCR4   = KEND1B
               KEND2   = KSCR4 + LSCR4
               LWORK2  = LWORK - KEND2
C
               CALL SO_MEMMAX ('SO_SIGAI.2',LWORK2)
               IF (LWORK2 .LT. 0)
     &             CALL STOPIT('SO_SIGAI.2',' ',KEND2,LWORK)
C
               NTOTGM = MAX(NBAS(ISGAM),1)
               NTOTBE = MAX(NBAS(ISBETA),1)
               NTOTC  = MAX(NVIR(ISYMC),1)
C
               KOFF2  = KSCR1 + IAODIS(ISBETA,ISGAM)
               KOFF3  = IT1AO(ISBETA,ISYMK) + 1
               KOFF4  = ILMVIR(ISGAM) + 1
C
C-------------------------------------------------------------------
C                             ~
C              Generate two ( k c | alfa delta ) in KSCR2 and KSCR3.
C-------------------------------------------------------------------
C
               CALL DGEMM('T','N',NBAS(ISGAM),NRHF(ISYMK),
     &                    NBAS(ISBETA),ONE,WORK(KOFF2),NTOTBE,
     &                    BTR1D(KOFF3),NTOTBE,ZERO,WORK(KSCR4),NTOTGM)
C
               CALL DGEMM('T','N',NVIR(ISYMC),NRHF(ISYMK),
     &                    NBAS(ISGAM),ONE,CMO(KOFF4),NTOTGM,
     &                    WORK(KSCR4),NTOTGM,ZERO,WORK(KPOS2),NTOTC)
C
               CALL DGEMM('T','N',NBAS(ISGAM),NRHF(ISYMK),
     &                    NBAS(ISBETA),-ONE,WORK(KOFF2),NTOTBE,
     &                    BTR1E(KOFF3),NTOTBE,ZERO,WORK(KSCR4),NTOTGM)
C
               CALL DGEMM('T','N',NVIR(ISYMC),NRHF(ISYMK),
     &                    NBAS(ISGAM),ONE,CMO(KOFF4),NTOTGM,
     &                    WORK(KSCR4),NTOTGM,ZERO,WORK(KPOS3),NTOTC)
C
               ISBETA = ISYMK
               ISGAM  = MULD2H(ISBETA,ISBEGA)
C
               IF ( ISYMC .NE. MULD2H(ISGAM,ISYMTR) ) THEN
                  WRITE(LUPRI,*)' ERROR in SO_SIGAI: '
                  WRITE(LUPRI,*)' ISYMC .NE.  MULD2H(ISGAM,ISYMTR)'
                  CALL QUIT(' ERROR in SO_SIGAI: ')
               END IF
C
               LSCR4  = NBAS(ISGAM)*NRHF(ISYMK)
C
               KSCR4   = KEND1B
               KEND3   = KSCR4 + LSCR4
               LWORK3  = LWORK - KEND3
C
               CALL SO_MEMMAX ('SO_SIGAI.3',LWORK3)
               IF (LWORK3 .LT. 0)
     &             CALL STOPIT('SO_SIGAI.3',' ',KEND3,LWORK)
C
               NTOTGM = MAX(NBAS(ISGAM),1)
               NTOTBE = MAX(NBAS(ISBETA),1)
               NTOTC  = MAX(NVIR(ISYMC),1)
C
               KOFF5  = KSCR1 + IAODIS(ISBETA,ISGAM)
               KOFF6  = ILMRHF(ISBETA) + 1
               KOFF7  = IMATAV(ISGAM,ISYMC) + 1
C
C----------------------------------------------------------------
C                               ~
C              Generate two ( k c | alfa delta ) and add to KSCR2
C              and KSCR3.
C----------------------------------------------------------------
C
               CALL DGEMM('T','N',NBAS(ISGAM),NRHF(ISYMK),
     &                    NBAS(ISBETA),-ONE,WORK(KOFF5),NTOTBE,
     &                    CMO(KOFF6),NTOTBE,ZERO,WORK(KSCR4),NTOTGM)
C
               CALL DGEMM('T','N',NVIR(ISYMC),NRHF(ISYMK),
     &                    NBAS(ISGAM),ONE,BTJ1D(KOFF7),NTOTGM,
     &                    WORK(KSCR4),NTOTGM,ONE,WORK(KPOS2),NTOTC)
C
               CALL DGEMM('T','N',NVIR(ISYMC),NRHF(ISYMK),
     &                    NBAS(ISGAM),-ONE,BTJ1E(KOFF7),NTOTGM,
     &                    WORK(KSCR4),NTOTGM,ONE,WORK(KPOS3),NTOTC)
C
  120       CONTINUE
C
C---------------------------------------------------------------
C                                            ~ ~
C              Generate Sigma(alfa,i) from ( k c | alfa delta ).
C              Store this as (i, alfa), to avoid strided access
C              to output array.
C---------------------------------------------------------------
C
            KOFF10 = IT1AO(ISALFA,ISYMI) + (ALFA-1)*NRHF(ISYMI)+1
C
            NDIM   = NT1AM(ISYMCK)
            NTOT   = MAX(NDIM,1)
            KOFFK  = IT2BCD(ISYMCK,ISYMI) + 1
C
            CALL DGEMV('T',NDIM,NRHF(ISYMI),ONE,
     &                 T2M1(KOFFK),NTOT,WORK(KSCR2),1,ONE,
     &                 SIGAI1(KOFF10),1)
C
            CALL DGEMV('T',NDIM,NRHF(ISYMI),ONE,
     &                 T2M1(KOFFK),NTOT,WORK(KSCR3),1,ONE,
     &                 SIGAI2(KOFF10),1)
C

C
  110    CONTINUE
C
  100 CONTINUE
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_SIGAI')
C
      RETURN
      END
