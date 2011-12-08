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
C  /* Deck cc3_t3int */
      SUBROUTINE CC3_T3INT(XINT,XLAMDP,XLAMDH,C1AM,ISYMTR,
     *                     WORK,LWORK,IDEL,ISYDEL,IOPT,
     *                     LUOUT1,FNOUT1,LUOUT2,FNOUT2)
C
C     Symmetry by Henrik Koch and Poul Joergensen. 13-Jan-1995
C
C     Purpose: Calculate integrals used in CC3 model T3 amplitudes.
C
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION XINT(*),WORK(LWORK)
      DIMENSION XLAMDP(*),XLAMDH(*),C1AM(*)
#include <ccinftap.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CHARACTER*(*) FNOUT1, FNOUT2
C
      CALL QENTER('CC3_T3INT')
C
      IF (IOPT .EQ. 1) THEN
C
         ISYCKD = MULD2H(ISYDEL,ISYMOP)
         ISYCKJ = ISYCKD
C
C---------------------------------
C        Allocation of work space.
C---------------------------------
C
         KXCKD = 1
         KXCKJ = KXCKD + NCKATR(ISYCKD)
         KEND1 = KXCKJ + NCKI(ISYCKJ)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient core in CC3_CKDINT')
         ENDIF
C
C----------------------------------------
C        Calculate transformed integrals.
C----------------------------------------
C
         CALL DZERO(WORK(KXCKD),NCKATR(ISYCKD))
         CALL DZERO(WORK(KXCKJ),NCKI(ISYCKJ))
C
         CALL CC3_CKD1(XINT,WORK(KXCKD),WORK(KXCKJ),XLAMDP(NT1AOX+1),
     *                 XLAMDH,XLAMDP(NT1AOX+1),XLAMDH,1,1,1,1,
     *                 WORK(KEND1),LWRK1,IDEL,ISYDEL)
C
      ELSEIF (IOPT .EQ. 2) THEN
C
         ISYCKD = MULD2H(ISYMTR,MULD2H(ISYDEL,ISYMOP))
         ISYCKJ = ISYCKD
C
C---------------------------------
C        Allocation of work space.
C---------------------------------
C
         KXCKD  = 1
         KXCKJ  = KXCKD  + NCKATR(ISYCKD)
         KXLAM1 = KXCKJ  + NCKI(ISYCKJ)
         KXLAM2 = KXLAM1 + NMATAV(ISYMTR)
         KEND1  = KXLAM2 + NT1AO(ISYMTR)
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient core in CC3_CKDINT')
         ENDIF
C
C------------------------------------------
C        Calculate transformation matrices.
C------------------------------------------
C
         DO 100 ISYML = 1,NSYM
C
            ISYMC = MULD2H(ISYML,ISYMTR)
            ISYMA = ISYML
C
            KOFF1 = ILMRHF(ISYML) + 1
            KOFF2 = IT1AM(ISYMC,ISYML) + 1
            KOFF3 = KXLAM1 + IMATAV(ISYMA,ISYMC)
C
            NBASA = MAX(NBAS(ISYMA),1)
            NVIRC = MAX(NVIR(ISYMC),1)
C
            CALL DGEMM('N','T',NBAS(ISYMA),NVIR(ISYMC),NRHF(ISYML),
     *                 -ONE,XLAMDP(KOFF1),NBASA,C1AM(KOFF2),NVIRC,
     *                 ZERO,WORK(KOFF3),NBASA)
C
  100    CONTINUE
C
         DO 110 ISYMK = 1,NSYM
C
            ISYMD = MULD2H(ISYMK,ISYMTR)
            ISYMA = ISYMD
C
            KOFF1 = ILMVIR(ISYMD) + 1
            KOFF2 = IT1AM(ISYMD,ISYMK) + 1
            KOFF3 = KXLAM2 + IT1AO(ISYMA,ISYMK)
C
            NBASA = MAX(NBAS(ISYMA),1)
            NVIRD = MAX(NVIR(ISYMD),1)
C
            CALL DGEMM('N','N',NBAS(ISYMA),NRHF(ISYMK),NVIR(ISYMD),
     *                 ONE,XLAMDH(KOFF1),NBASA,C1AM(KOFF2),NVIRD,
     *                  ZERO,WORK(KOFF3),NBASA)
C
  110    CONTINUE
C
C----------------------------------------
C        Calculate transformed integrals.
C----------------------------------------
C
         CALL DZERO(WORK(KXCKD),NCKATR(ISYCKD))
         CALL DZERO(WORK(KXCKJ),NCKI(ISYCKJ))
C
         CALL CC3_CKD1(XINT,WORK(KXCKD),WORK(KXCKJ),
     *                 WORK(KXLAM1),XLAMDH,XLAMDP(NT1AOX+1),
     *                 XLAMDH,ISYMTR,1,1,1,
     *                 WORK(KEND1),LWRK1,IDEL,ISYDEL)
C
         CALL CC3_CKD1(XINT,WORK(KXCKD),WORK(KXCKJ),
     *                 XLAMDP(NT1AOX+1),WORK(KXLAM2),XLAMDP(NT1AOX+1),
     *                 XLAMDH,1,ISYMTR,1,1,
     *                 WORK(KEND1),LWRK1,IDEL,ISYDEL)
C
         CALL CC3_CKD1(XINT,WORK(KXCKD),WORK(KXCKJ),
     *                 XLAMDP(NT1AOX+1),XLAMDH,WORK(KXLAM1),
     *                 WORK(KXLAM2),1,1,ISYMTR,ISYMTR,
     *                 WORK(KEND1),LWRK1,IDEL,ISYDEL)
C
      ELSE
         CALL QUIT('Incorrect specification of IOPT in CC3_CKDINT')
      ENDIF
C
C--------------------------------
C     Write to disk (ck|d alpha).
C--------------------------------
C
      ID     = IDEL - IBAS(ISYDEL)
C
      LENGTH = NCKATR(ISYCKD)
C
      IOFF = ICKDAO(ISYCKD,ISYDEL) + NCKATR(ISYCKD)*(ID - 1) + 1
C
      IF (LENGTH .GT. 0) THEN
         CALL PUTWA2(LUOUT1,FNOUT1,WORK(KXCKD),IOFF,LENGTH)
      ENDIF
C
      LENGTH = NCKI(ISYCKJ)
C
      IOFF  = ICKID(ISYCKJ,ISYDEL) + NCKI(ISYCKJ)*(ID - 1) + 1
C
      IF (LENGTH .GT. 0) THEN
         CALL PUTWA2(LUOUT2,FNOUT2,WORK(KXCKJ),IOFF,LENGTH)
      ENDIF
C
      CALL QEXIT('CC3_T3INT')
      RETURN
      END
C  /* Deck cc3_ckd1 */
      SUBROUTINE CC3_CKD1(XINT,XCKD,XCKJ,XLAMD1,XLAMD2,XLAMD3,
     *                    XLAMD4,ISYM1,ISYM2,ISYM3,ISYM4,
     *                    WORK,LWORK,IDEL,ISYDEL)
C
C     Symmetry by Henrik Koch and Poul Joergensen. 13-Jan-1995
C
C     Purpose: Calculate integrals used in CC3 model T3 amplitudes.
C
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION XINT(*),XCKD(*),XCKJ(*),WORK(LWORK)
      DIMENSION XLAMD1(*),XLAMD2(*),XLAMD3(*),XLAMD4(*)
#include <ccorb.h>
#include <ccsdsym.h>
C
      ISYDIS = MULD2H(ISYDEL,ISYMOP)

      DO 100  ISYMG = 1,NSYM
C
         ISYMAB = MULD2H(ISYMG, ISYDIS)
         ISYMAK = MULD2H(ISYMAB,ISYM2)
         ISYMCK = MULD2H(ISYMAK,ISYM1)
         ISYMD  = MULD2H(ISYM3, ISYMG)
         ISYMJ  = MULD2H(ISYM4, ISYMG)
C
C----------------------------
C        Allocate work space.
C----------------------------
C
         KINT1 = 1
         KSCR1 = KINT1 + NT1AM(ISYMCK)*NBAS(ISYMG)
         KSCR2 = KSCR1 + N2BST(ISYMAB)
         KEND1 = KSCR2 + NT1AO(ISYMAK)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient space in CC3_BCDINT')
         ENDIF
C
         DO 110 G = 1,NBAS(ISYMG)
C
            KOFF1 = IDSAOG(ISYMG,ISYDIS) + NNBST(ISYMAB)*(G-1) + 1
            CALL CCSD_SYMSQ(XINT(KOFF1),ISYMAB,WORK(KSCR1))
C
            DO 120 ISYMK = 1,NSYM
C
               ISYMB = MULD2H(ISYMK,ISYM2)
               ISYMA = MULD2H(ISYMB,ISYMAB)
               ISYMC = MULD2H(ISYMA,ISYM1)
C
               KOFF2 = KSCR1 + IAODIS(ISYMA,ISYMB)
               KOFF3 = IT1AO(ISYMB,ISYMK) + 1
               KOFF4 = KSCR2 + IT1AO(ISYMA,ISYMK)
C
               NBASA = MAX(NBAS(ISYMA),1)
               NBASB = MAX(NBAS(ISYMB),1)
C
               CALL DGEMM('N','N',NBAS(ISYMA),NRHF(ISYMK),NBAS(ISYMB),
     *                    ONE,WORK(KOFF2),NBASA,XLAMD2(KOFF3),NBASB,
     *                    ZERO,WORK(KOFF4),NBASA)
C
  120       CONTINUE
C
            DO 130 ISYMK = 1,NSYM
C
               ISYMA = MULD2H(ISYMK,ISYMAK)
               ISYMC = MULD2H(ISYMA,ISYM1)
C
               KOFF2 = IMATAV(ISYMA,ISYMC) + 1
               KOFF3 = KSCR2 + IT1AO(ISYMA,ISYMK)
               KOFF4 = KINT1 + NT1AM(ISYMCK)*(G - 1)+IT1AM(ISYMC,ISYMK)
C
               NBASA = MAX(NBAS(ISYMA),1)
               NVIRC = MAX(NVIR(ISYMC),1)
C
               CALL DGEMM('T','N',NVIR(ISYMC),NRHF(ISYMK),NBAS(ISYMA),
     *                    ONE,XLAMD1(KOFF2),NBASA,WORK(KOFF3),NBASA,
     *                    ZERO,WORK(KOFF4),NVIRC)
C
  130       CONTINUE
C
  110    CONTINUE
C
         NTOTCK = MAX(NT1AM(ISYMCK),1)
         NBASG  = MAX(NBAS(ISYMG),1)
C
         KOFF1 = IMATAV(ISYMG,ISYMD)  + 1
         KOFF2 = ICKATR(ISYMCK,ISYMD) + 1
C
         CALL DGEMM('N','N',NT1AM(ISYMCK),NVIR(ISYMD),NBAS(ISYMG),
     *              ONE,WORK(KINT1),NTOTCK,XLAMD3(KOFF1),NBASG,
     *              ONE,XCKD(KOFF2),NTOTCK)
C
         KOFF1 = IT1AO(ISYMG,ISYMJ) + 1
         KOFF2 = ICKI(ISYMCK,ISYMJ) + 1
C
         CALL DGEMM('N','N',NT1AM(ISYMCK),NRHF(ISYMJ),NBAS(ISYMG),
     *              ONE,WORK(KINT1),NTOTCK,XLAMD4(KOFF1),NBASG,
     *              ONE,XCKJ(KOFF2),NTOTCK)
C
  100 CONTINUE
C
      RETURN
      END
C  /* Deck cc3_ckdsor */
      SUBROUTINE CC3_CKDSOR(XCKD,WORK,LWORK,ISYCKD)
C
C     Symmetry by Henrik Koch and Poul Joergensen. 13-Jan-1995
C
C     Purpose: Resort integral distribution from
C              X(ck,d) to X(dk,c).
C
C
#include <implicit.h>
      DIMENSION XCKD(*),WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
C
      IF (LWORK .LT. NCKATR(ISYCKD))
     *   CALL QUIT('Insufficient work space in CC3_CKDSOR')
C
      DO 100 ISYMD = 1,NSYM
C
         ISYMCK = MULD2H(ISYMD,ISYCKD)
C
         DO 110 D = 1,NVIR(ISYMD)
C
            DO 120 ISYMK = 1,NSYM
C
               ISYMDK = MULD2H(ISYMD,ISYMK)
               ISYMC  = MULD2H(ISYMK,ISYMCK)
C
               DO 130 K = 1,NRHF(ISYMK)
C
                  DO 140 C = 1,NVIR(ISYMC)
C
                     NCKD = ICKATR(ISYMCK,ISYMD)
     *                    + NT1AM(ISYMCK)*(D - 1)
     *                    + IT1AM(ISYMC,ISYMK)
     *                    + NVIR(ISYMC)*(K - 1) + C
C
                     NDKC = ICKATR(ISYMDK,ISYMC)
     *                    + NT1AM(ISYMDK)*(C - 1)
     *                    + IT1AM(ISYMD,ISYMK)
     *                    + NVIR(ISYMD)*(K - 1) + D
C
                       WORK(NDKC) = XCKD(NCKD)
C
  140             CONTINUE
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL DCOPY(NCKATR(ISYCKD),WORK,1,XCKD,1)
C
      RETURN
      END
