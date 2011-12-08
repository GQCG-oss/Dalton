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
!  /* Deck cc3_t3int3 */
      SUBROUTINE CC3_T3INT3(XINT,XLAMDP,XLAMDH,C1AM,ISYMTR,
     *                      WORK,LWORK,IDEL,ISYDEL,LU3SRT,FN3SRT,
     *                      LU3SRT2,FN3SRT2,LU3SRT3,FN3SRT3,
     *                      LUCKJD,FNCKJD,LUCKJD2,FNCKJD2,
     *                      LUCKJD3,FNCKJD3)
!
!     Symmetry by Henrik Koch and Poul Joergensen. 13-Jan-1995
!     Purpose: Calculate integrals used in CC3 model T3 amplitudes.
!
!     Kasper Hald, January 2001.
!     Transformed to use in triplet case.
!
!
      IMPLICIT NONE
!
      INTEGER ISYMTR, LWORK, IDEL, ISYDEL
      INTEGER KEND1, LWRK1, KOFF1, KOFF2, KOFF3
      INTEGER ISYCKD, ISYCKJ, KXCKD, KXCKJ, KXLAM1, KXLAM2
      INTEGER ISYML, ISYMC, ISYMA, NBASA, NVIRC, ISYMK
      INTEGER ISYMD, NVIRD, ID, LENGTH, IOFF
      INTEGER KXCKD2, KXCKJ2, KXCKD3, KXCKJ3
      INTEGER LU3SRT, LU3SRT2, LU3SRT3, LUCKJD, LUCKJD2, LUCKJD3
!
#if defined (SYS_CRAY)
      REAL XINT(*),WORK(LWORK), XDOT, DDOT
      REAL XLAMDP(*), XLAMDH(*), C1AM(*), ZERO, ONE, XMONE
#else
      DOUBLE PRECISION XINT(*),WORK(LWORK), XDOT, DDOT
      DOUBLE PRECISION XLAMDP(*), XLAMDH(*), C1AM(*), ZERO, ONE, XMONE
#endif
!
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, XMONE = -1.0D0)
!
#include <priunit.h>
#include <ccinftap.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
!
      CHARACTER*(*) FN3SRT, FN3SRT2, FN3SRT3, FNCKJD, FNCKJD2, FNCKJD3
!
      CALL QENTER('CC3_T3INT3')
!
!
!-----------------------------------------
!        Calculation of first integral.
!-----------------------------------------
!
      ISYCKD = MULD2H(ISYDEL,ISYMOP)
      ISYCKJ = ISYCKD
!
!---------------------------------
!        Allocation of work space.
!---------------------------------
!
      KXCKD = 1
      KXCKJ = KXCKD + NCKATR(ISYCKD)
      KEND1 = KXCKJ + NCKI(ISYCKJ)
      LWRK1 = LWORK - KEND1
!
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient core in CC3_T3INT3')
      ENDIF
!
!----------------------------------------
!        Calculate transformed integrals.
!----------------------------------------
!
      CALL DZERO(WORK(KXCKD),NCKATR(ISYCKD))
      CALL DZERO(WORK(KXCKJ),NCKI(ISYCKJ))
!
      CALL CC3_CKD1(XINT,WORK(KXCKD),WORK(KXCKJ),XLAMDP(NT1AOX+1),
     *              XLAMDH,XLAMDP(NT1AOX+1),XLAMDH,1,1,1,1,
     *              WORK(KEND1),LWRK1,IDEL,ISYDEL)
!
!--------------------------------
!     Write to disk (ck|d alpha).
!--------------------------------
!
      ID     = IDEL - IBAS(ISYDEL)
!
      LENGTH = NCKATR(ISYCKD)
!
      IOFF = ICKDAO(ISYCKD,ISYDEL) + NCKATR(ISYCKD)*(ID - 1) + 1
!
      IF (LENGTH .GT. 0) THEN
         CALL PUTWA2(LU3SRT,FN3SRT,WORK(KXCKD),IOFF,LENGTH)
      ENDIF
!
      LENGTH = NCKI(ISYCKJ)
!
      IOFF  = ICKID(ISYCKJ,ISYDEL) + NCKI(ISYCKJ)*(ID - 1) + 1
!
      IF (LENGTH .GT. 0) THEN
         CALL PUTWA2(LUCKJD,FNCKJD,WORK(KXCKJ),IOFF,LENGTH)
      ENDIF
!
!-----------------------------------
!     Print the norm if desired.
!-----------------------------------
!
      IF (IPRINT .GT. 55) THEN
         XDOT = DDOT(NCKATR(ISYCKD),WORK(KXCKD),1,WORK(KXCKD),1)
         WRITE(LUPRI,*) 'Norm of CKD integral ',XDOT
         XDOT = DDOT(NCKI(ISYCKJ),WORK(KXCKJ),1,WORK(KXCKJ),1)
         WRITE(LUPRI,*) 'Norm of CKJ integral ',XDOT
      ENDIF
!
!------------------------------------------
!     Calculate second and third integrals.
!------------------------------------------
!
      ISYCKD = MULD2H(ISYMTR,MULD2H(ISYDEL,ISYMOP))
      ISYCKJ = ISYCKD
!
!---------------------------------
!        Allocation of work space.
!---------------------------------
!
      KXCKD  = 1
      KXCKJ  = KXCKD  + NCKATR(ISYCKD)
      KXCKD2 = KXCKJ  + NCKI(ISYCKJ)
      KXCKJ2 = KXCKD2 + NCKATR(ISYCKD)
      KXLAM1 = KXCKJ2 + NCKI(ISYCKJ)
      KXLAM2 = KXLAM1 + NMATAV(ISYMTR)
      KEND1  = KXLAM2 + NT1AO(ISYMTR)
      LWRK1  = LWORK  - KEND1
!
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient core in CC3_T3INT3')
      ENDIF
!
!------------------------------------------
!        Calculate transformation matrices.
!------------------------------------------
!
      DO ISYML = 1,NSYM
!
         ISYMC = MULD2H(ISYML,ISYMTR)
         ISYMA = ISYML
!
         KOFF1 = ILMRHF(ISYML) + 1
         KOFF2 = IT1AM(ISYMC,ISYML) + 1
         KOFF3 = KXLAM1 + IMATAV(ISYMA,ISYMC)
!
         NBASA = MAX(NBAS(ISYMA),1)
         NVIRC = MAX(NVIR(ISYMC),1)
!
         CALL DGEMM('N','T',NBAS(ISYMA),NVIR(ISYMC),NRHF(ISYML),
     *              -ONE,XLAMDP(KOFF1),NBASA,C1AM(KOFF2),NVIRC,
     *              ZERO,WORK(KOFF3),NBASA)
!
      ENDDO
!
      DO ISYMK = 1,NSYM
!
         ISYMD = MULD2H(ISYMK,ISYMTR)
         ISYMA = ISYMD
!
         KOFF1 = ILMVIR(ISYMD) + 1
         KOFF2 = IT1AM(ISYMD,ISYMK) + 1
         KOFF3 = KXLAM2 + IT1AO(ISYMA,ISYMK)
!
         NBASA = MAX(NBAS(ISYMA),1)
         NVIRD = MAX(NVIR(ISYMD),1)
!
         CALL DGEMM('N','N',NBAS(ISYMA),NRHF(ISYMK),NVIR(ISYMD),
     *              ONE,XLAMDH(KOFF1),NBASA,C1AM(KOFF2),NVIRD,
     *               ZERO,WORK(KOFF3),NBASA)
!
      ENDDO
!
!----------------------------------------
!        Calculate transformed integrals.
!        Have two different kinds. The
!        g^(1) and g^(2) will be calculated
!        from these 2 kinds later
!----------------------------------------
!
      CALL DZERO(WORK(KXCKD),NCKATR(ISYCKD))
      CALL DZERO(WORK(KXCKJ),NCKI(ISYCKJ))
      CALL DZERO(WORK(KXCKD2),NCKATR(ISYCKD))
      CALL DZERO(WORK(KXCKJ2),NCKI(ISYCKJ))
!
!--------------------------------------------
!        Calculate g(c-bar,k,b,del) and
!                  g(c-bar,k,del,j)
!--------------------------------------------
!
      CALL CC3_CKD1(XINT,WORK(KXCKD),WORK(KXCKJ),
     *              WORK(KXLAM1),XLAMDH,XLAMDP(NT1AOX+1),
     *              XLAMDH,ISYMTR,1,1,1,
     *              WORK(KEND1),LWRK1,IDEL,ISYDEL)
!
!--------------------------------------------
!        Calculate g(c,k-bar,b,del) and
!                  g(c,k-bar,del,j)
!--------------------------------------------
!
      CALL CC3_CKD1(XINT,WORK(KXCKD),WORK(KXCKJ),
     *              XLAMDP(NT1AOX+1),WORK(KXLAM2),XLAMDP(NT1AOX+1),
     *              XLAMDH,1,ISYMTR,1,1,
     *              WORK(KEND1),LWRK1,IDEL,ISYDEL)
!
!--------------------------------------------
!        Calculate g(c,k,b-bar,del) and
!                  g(c,k,del,j-bar)
!--------------------------------------------
!
      CALL CC3_CKD1(XINT,WORK(KXCKD2),WORK(KXCKJ2),
     *              XLAMDP(NT1AOX+1),XLAMDH,WORK(KXLAM1),
     *              WORK(KXLAM2),1,1,ISYMTR,ISYMTR,
     *              WORK(KEND1),LWRK1,IDEL,ISYDEL)
!
!
!--------------------
!     Write to disk.
!--------------------
!
      ID     = IDEL - IBAS(ISYDEL)
!
      LENGTH = NCKATR(ISYCKD)
!
      IOFF = ICKDAO(ISYCKD,ISYDEL) + NCKATR(ISYCKD)*(ID - 1) + 1
!
      IF (LENGTH .GT. 0) THEN
         CALL PUTWA2(LU3SRT2,FN3SRT2,WORK(KXCKD),IOFF,LENGTH)
         CALL PUTWA2(LU3SRT3,FN3SRT3,WORK(KXCKD2),IOFF,LENGTH)
      ENDIF
!
      LENGTH = NCKI(ISYCKJ)
!
      IOFF  = ICKID(ISYCKJ,ISYDEL) + NCKI(ISYCKJ)*(ID - 1) + 1
!
      IF (LENGTH .GT. 0) THEN
         CALL PUTWA2(LUCKJD2,FNCKJD2,WORK(KXCKJ),IOFF,LENGTH)
         CALL PUTWA2(LUCKJD3,FNCKJD3,WORK(KXCKJ2),IOFF,LENGTH)
      ENDIF
!
!-----------------------------------
!     Print the norm if desired.
!-----------------------------------
!
      IF (IPRINT .GT. 55) THEN
         XDOT = DDOT(NCKATR(ISYCKD),WORK(KXCKD),1,WORK(KXCKD),1)
         WRITE(LUPRI,*) 'Norm of CKD-2 integral ',XDOT
         XDOT = DDOT(NCKATR(ISYCKD),WORK(KXCKD2),1,WORK(KXCKD2),1)
         WRITE(LUPRI,*) 'Norm of CKD-3 integral ',XDOT
         XDOT = DDOT(NCKI(ISYCKJ),WORK(KXCKJ),1,WORK(KXCKJ),1)
         WRITE(LUPRI,*) 'Norm of CKJ-2 integral ',XDOT
         XDOT = DDOT(NCKI(ISYCKJ),WORK(KXCKJ2),1,WORK(KXCKJ2),1)
         WRITE(LUPRI,*) 'Norm of CKJ-3 integral ',XDOT
      ENDIF
!
      CALL QEXIT('CC3_T3INT3')
!
      RETURN
      END
!  /* Deck cc3_ckjsor */
      SUBROUTINE CC3_CKJSOR(XCKJ,WORK,LWORK,ISYCKJ)
!
!     Written by K. Hald, Jan 2001.
!
!     Purpose: Resort integral distribution from
!              X(ck,j) to X(cj,k).
!
!
      IMPLICIT NONE
!
      INTEGER LWORK, ISYCKJ, ISYMJ, ISYMCK, ISYMK, ISYMC, ISYMCJ
      INTEGER NCKJ, NCJK
!
#if defined (SYS_CRAY)
      REAL XCKD(*),WORK(LWORK)
#else
      DOUBLE PRECISION XCKJ(*), WORK(LWORK)
#endif
!
#include <ccorb.h>
#include <ccsdsym.h>
!
      CALL QENTER('CC3_CKJSOR')
!
      IF (LWORK .LT. NCKI(ISYCKJ))
     *   CALL QUIT('Insufficient work space in CC3_CKJSOR')
!
      DO 100 ISYMJ = 1,NSYM
!
         ISYMCK = MULD2H(ISYMJ,ISYCKJ)
!
         DO 110 J = 1,NRHF(ISYMJ)
!
            DO 120 ISYMK = 1,NSYM
!
               ISYMC  = MULD2H(ISYMK,ISYMCK)
               ISYMCJ = MULD2H(ISYMC,ISYMJ)
!
               DO 130 K = 1,NRHF(ISYMK)
!
                  DO 140 C = 1,NVIR(ISYMC)
!
                     NCKJ = ICKI(ISYMCK,ISYMJ)
     *                    + NT1AM(ISYMCK)*(J - 1)
     *                    + IT1AM(ISYMC,ISYMK) 
     *                    + NVIR(ISYMC)*(K - 1) + C
!
                     NCJK = ICKI(ISYMCJ,ISYMK)
     *                    + NT1AM(ISYMCJ)*(K - 1)
     *                    + IT1AM(ISYMC,ISYMJ) 
     *                    + NVIR(ISYMC)*(J - 1) + C
!
                     WORK(NCJK) = XCKJ(NCKJ)
!
  140             CONTINUE
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
!
      CALL DCOPY(NCKI(ISYCKJ),WORK,1,XCKJ,1)
!
      CALL QEXIT('CC3_CKJSOR')
!
      RETURN
      END
