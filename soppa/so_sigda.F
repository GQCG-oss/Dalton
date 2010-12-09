C  /* Deck so_sigda */
      SUBROUTINE SO_SIGDA(SIGDA1,LSIGDA1,SIGDA2,LSIGDA2,T2MP,LT2MP,
     &                    DSRHF,LDSRHF,BTR1E,LBTR1E,BTR1D,LBTR1D,
     &                    BTJ1E,LBTJ1E,BTJ1D,LBTJ1D,CMO,LCMO,
     &                    IDEL,ISYMDEL,ISYDIS,ISYRES,ISYMTR,WORK,
     &                    LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, November 1995
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C
C     PURPOSE: Calculate SIGDA1(DELTA,A) and SIGDA2(DELTA,A)
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION SIGDA1(LSIGDA1), SIGDA2(LSIGDA2)
      DIMENSION T2MP(LT2MP),     BTR1E(LBTR1E),  BTR1D(LBTR1D)
      DIMENSION BTJ1E(LBTJ1E),   BTJ1D(LBTJ1D),  DSRHF(LDSRHF)
      DIMENSION CMO(LCMO)
      DIMENSION WORK(LWORK)
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_SIGDA')
C
      ISYMI = ISYMDEL
      ISYMA = MULD2H(ISYMI,ISYRES)
C
      DO 100 ISYML = 1,NSYM
C
         ISGAM  = ISYML
	 ISALBE = MULD2H(ISGAM,ISYDIS)
         ISYMLI = MULD2H(ISYML,ISYMI)
	 ISDEL  = MULD2H(ISGAM,ISALBE)
         if (isdel .ne. isymdel) then
            print*,'Keld: Something is WRONG in SO_SIGDA'
            stop
         end if
C
         LSCR1  = N2BST(ISALBE)
C
         KSCR1   = 1
         KEND1   = KSCR1 + LSCR1
         LWORK1  = LWORK - KEND1
C
         CALL SO_MEMMAX ('SO_SIGDA.1',LWORK1)
         IF (LWORK1 .LT. 0) CALL STOPIT('SO_SIGDA.1',' ',KEND1,LWORK)
C
         DO 110 L = 1,NRHF(ISYML)
C
            KOFF1 = IDSRHF(ISALBE,ISYML) + NNBST(ISALBE)*(L-1) + 1
C
C-----------------------------------------------------------------------
C           Get a squared set of ( alfa beta | l delta ) for given l and
C           delta.
C-----------------------------------------------------------------------
C
            CALL CCSD_SYMSQ(DSRHF(KOFF1),ISALBE,WORK(KSCR1))
C
            DO 120 ISYMK = 1,NSYM
C
               ISALFA = MULD2H(ISYMK,ISYMTR)
               ISBETA = MULD2H(ISALFA,ISALBE)
               ISYMC  = ISBETA
C
               LSCR2  = NRHF(ISYMK)*NVIR(ISYMC)
               LSCR3  = NRHF(ISYMK)*NVIR(ISYMC)
               LSCR4  = NVIR(ISYMC)*NRHF(ISYMK)*NVIR(ISYMA)
               LSCR5  = NRHF(ISYMK)*NBAS(ISBETA)
C
               KSCR2   = KEND1
               KSCR3   = KSCR2 + LSCR2
               KSCR4   = KSCR3 + LSCR3
               KSCR5   = KSCR4 + LSCR4
               KEND2   = KSCR5 + LSCR5
               LWORK2  = LWORK - KEND2
C
               CALL SO_MEMMAX ('SO_SIGDA.2',LWORK2)
               IF (LWORK2 .LT. 0) 
     &             CALL STOPIT('SO_SIGDA.2',' ',KEND2,LWORK)
C
               NTOTK  = MAX(NRHF(ISYMK),1)
               NTOTC  = MAX(NVIR(ISYMC),1)
               NTOTAL = MAX(NBAS(ISALFA),1)
               NTOTBE = MAX(NBAS(ISBETA),1)
C
               KOFF2  = IT1AO(ISALFA,ISYMK) + 1
               KOFF3  = KSCR1 + IAODIS(ISALFA,ISBETA)
               KOFF4  = ILMVIR(ISBETA) + 1
C
C----------------------------------------------------------------
C                             ~
C              Generate two ( k c | l delta ) in KSCR2 and KSCR3.
C----------------------------------------------------------------
C
               CALL DGEMM('T','N',NRHF(ISYMK),NBAS(ISBETA),
     &                    NBAS(ISALFA),ONE,BTR1D(KOFF2),NTOTAL,
     &                    WORK(KOFF3),NTOTAL,ZERO,WORK(KSCR5),NTOTK)
C
               CALL DGEMM('N','N',NRHF(ISYMK),NVIR(ISYMC),
     &                    NBAS(ISBETA),ONE,WORK(KSCR5),NTOTK,
     &                    CMO(KOFF4),NTOTBE,ZERO,WORK(KSCR2),NTOTK)
C
C
               CALL DGEMM('T','N',NRHF(ISYMK),NBAS(ISBETA),
     &                    NBAS(ISALFA),-ONE,BTR1E(KOFF2),NTOTAL,
     &                    WORK(KOFF3),NTOTAL,ZERO,WORK(KSCR5),NTOTK)
C
               CALL DGEMM('N','N',NRHF(ISYMK),NVIR(ISYMC),
     &                    NBAS(ISBETA),ONE,WORK(KSCR5),NTOTK,
     &                    CMO(KOFF4),NTOTBE,ZERO,WORK(KSCR3),NTOTK)
C
               ISALFA = ISYMK
               ISBETA = MULD2H(ISYMK,ISYMLI)
C
               IF ( ISYMC .NE. MULD2H(ISBETA,ISYMTR) ) THEN
                  WRITE(LUPRI,*)' ERROR in SO_SIGDA: '
                  WRITE(LUPRI,*)' ISYMC .NE.  MULD2H(ISBETA,ISYMTR)'
                  CALL QUIT(' ERROR in SO_SIGDA: ')
               END IF
C
               LSCR5  = NRHF(ISYMK)*NBAS(ISBETA)
C
               KSCR5   = KSCR4 + LSCR4
               KEND3   = KSCR5 + LSCR5
               LWORK3  = LWORK - KEND3
C
               CALL SO_MEMMAX ('SO_SIGDA.3',LWORK3)
               IF (LWORK3 .LT. 0) 
     &             CALL STOPIT('SO_SIGDA.3',' ',KEND3,LWORK)
C
               NTOTAL = MAX(NBAS(ISALFA),1)
               NTOTBE = MAX(NBAS(ISBETA),1)
C
               KOFF5  = KSCR1 + IAODIS(ISALFA,ISBETA)
	       KOFF6  = ILMRHF(ISALFA) + 1
	       KOFF7  = IMATAV(ISBETA,ISYMC) + 1
C
C-----------------------------------------------------------------
C                               ~
C              Generate two ( k c | l delta ) and add to KSCR2 and
C              KSCR3.
C-----------------------------------------------------------------
C
               CALL DGEMM('T','N',NRHF(ISYMK),NBAS(ISBETA),
     &                    NBAS(ISALFA),-ONE,CMO(KOFF6),NTOTAL,
     &                    WORK(KOFF5),NTOTAL,ZERO,WORK(KSCR5),NTOTK)
C
               CALL DGEMM('N','N',NRHF(ISYMK),NVIR(ISYMC),
     &                    NBAS(ISBETA),ONE,WORK(KSCR5),NTOTK,
     &                    BTJ1D(KOFF7),NTOTBE,ONE,WORK(KSCR2),NTOTK)
C
C
               CALL DGEMM('T','N',NRHF(ISYMK),NBAS(ISBETA),
     &                    NBAS(ISALFA),ONE,CMO(KOFF6),NTOTAL,
     &                    WORK(KOFF5),NTOTAL,ZERO,WORK(KSCR5),NTOTK)
C
               CALL DGEMM('N','N',NRHF(ISYMK),NVIR(ISYMC),
     &                    NBAS(ISBETA),ONE,WORK(KSCR5),NTOTK,
     &                    BTJ1E(KOFF7),NTOTBE,ONE,WORK(KSCR3),NTOTK)
C
C-----------------------------------------------------------------------
C              Extract T2MP amplitudes for given J and given symmetries
C              of C, K and A.
C-----------------------------------------------------------------------
C
               CALL SO_SQAIT2(WORK(KSCR4),LSCR4,T2MP,LT2MP,ISYMC,ISYMK,
     &                        ISYMA,ISYML,L)
C
C-----------------------------------------------------------------
C                                                 ~ ~
C              Generate two Sigma(delta,a) from ( k c | l delta ).
C-----------------------------------------------------------------
C
               KOFF8  = IMATAV(ISDEL,ISYMA) + IDEL - IBAS(ISYMDEL)
C
               NDIM   = NVIR(ISYMC)*NRHF(ISYMK)
               NTOT   = MAX(NVIR(ISYMA),1)
C
               CALL DGEMV('N',NVIR(ISYMA),NDIM,ONE,
     &                    WORK(KSCR4),NTOT,WORK(KSCR2),1,ONE,
     &                    SIGDA1(KOFF8),NBAS(ISDEL))
C
               CALL DGEMV('N',NVIR(ISYMA),NDIM,ONE,
     &                    WORK(KSCR4),NTOT,WORK(KSCR3),1,ONE,
     &                    SIGDA2(KOFF8),NBAS(ISDEL))
C
  120       CONTINUE
C
  110    CONTINUE
C
  100 CONTINUE
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_SIGDA')
C
      RETURN
      END
