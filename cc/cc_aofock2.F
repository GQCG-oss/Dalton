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
C  /* Deck cc_aofock2 */
      SUBROUTINE CC_AOFOCK2(XINT,DENSIT,DENSPK,FOCK,WORK,LWORK,IDEL,
     *                      ISYDIS,ISYMD,ISYDEN,SQRINT)
C
C     Purpose: Calculate the two electron contribution to the
C              AO-fock matrix using matrix vector routines.
C
C     Written by Asger Halkier and Henrik Koch 27-4-95.
C
C     Debugged Ove Christiansen august 1995
C
C     Christof Haettig summer 1998: 
C         arbitrary point group symmetry of the integrals
C
C     Sonia Coriani autumn 1999:
C         Distinguish between packed and squared
C         integral distributions (for London integrals)
C         SQRINT = .FALSE., packed integral distribution  and
C                           Coulomb part calculated with packed density
C         SQRINT = .TRUE.,  squared integral distribution and
C                           Coulomb part calculated with squared density
C
C     Obs: It can be done as F(g>=d) = G(a>=b) I(a>=b,g,d) where
C          G(a>=b) = D(a,b) + D(b,a), the diagonal properly scaled
C
#include <implicit.h>
#include <maxorb.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, HALF = 0.5D0)
      DIMENSION XINT(*),DENSIT(*), DENSPK(*)
      DIMENSION FOCK(*),WORK(LWORK)
      LOGICAL SQRINT
#include <ccorb.h>
#include <symsq.h>
#include <ccsdsym.h>
C
      ISYHOP = MULD2H(ISYMD,ISYDIS)
C
C--------------------------------------------
C     start loop over symmetry blocks:
C--------------------------------------------
C
      DO 100 ISYMG = 1,NSYM
C
         IF (NBAS(ISYMG) .EQ. 0) GOTO 100
C
         D   = IDEL - IBAS(ISYMD)
C
         ISYMAB = MULD2H(ISYMG,ISYDIS)

C
C---------------------------------------------------------
C        calculate coulomb contribution:
C        For SQRINT=.FALSE. use packed form of density and integrals
C        For SQRINT=.TRUE.  use squared form of density and integrals
C---------------------------------------------------------
C
         IF (ISYMAB .EQ. ISYDEN) THEN
C
            KGD = IAODIS(ISYMG,ISYMD) + NBAS(ISYMG)*(D - 1) + 1
C
            IF (.NOT.SQRINT) THEN
               KOFF1 = IDSAOG(ISYMG,ISYDIS) + 1
               NTOBST = MAX(NNBST(ISYMAB),1)
               CALL DGEMV('T',NNBST(ISYMAB),NBAS(ISYMG),
     *                           TWO,XINT(KOFF1),NTOBST,
     *                           DENSPK,1,ONE,FOCK(KGD),1)
            ELSE IF (SQRINT) THEN
               KOFF1  = IDSAOGSQ(ISYMG,ISYDIS) + 1
               NTOBST = MAX(N2BST(ISYMAB),1)
               CALL DGEMV('T',N2BST(ISYMAB),NBAS(ISYMG),
     *                        TWO,XINT(KOFF1),NTOBST,
     *                        DENSIT,1,ONE,FOCK(KGD),1)
            END IF
C
         END IF
C
C-----------------------------------------------------------------
C        calculate exchange contribution in a loop over g:
C        For SQRINT=.FALSE. use packed form of density and integrals
C        For SQRINT=.TRUE.  use squared form of density and integrals
C-----------------------------------------------------------------
C
         IF (LWORK .LT. N2BST(ISYMAB)) THEN
            CALL QUIT('Insufficient work space in CC_AOFOCK2')
         ENDIF
C
         ISYMA = MULD2H(ISYMD,ISYDEN)
         ISYMB = MULD2H(ISYMA,ISYMAB)
C
         KAD = IAODIS(ISYMA,ISYMD) + NBAS(ISYMA)*(D - 1) + 1
         NTOTA = MAX(NBAS(ISYMA),1)
         NTOTG = MAX(NBAS(ISYMG),1)
C
         DO G = 1, NBAS(ISYMG)
C
            KGB = IAODIS(ISYMG,ISYMB) + G

            IF (.NOT.SQRINT) THEN

               KOFF1 = IDSAOG(ISYMG,ISYDIS) + NNBST(ISYMAB)*(G - 1) + 1

               CALL CCSD_SYMSQ(XINT(KOFF1),ISYMAB,WORK)
C
               KAB = IAODIS(ISYMA,ISYMB) + 1
               CALL DGEMV('T',NBAS(ISYMA),NBAS(ISYMB),-ONE,WORK(KAB),
     *                 NTOTA,DENSIT(KAD),1,ONE,FOCK(KGB),NTOTG)

            ELSE IF (SQRINT) THEN

               KOFF1 = IDSAOGSQ(ISYMG,ISYDIS) + N2BST(ISYMAB)*(G - 1) + 
     &                 IAODIS(ISYMA,ISYMB) + 1

               CALL DGEMV('T',NBAS(ISYMA),NBAS(ISYMB),-ONE,XINT(KOFF1),
     *                 NTOTA,DENSIT(KAD),1,ONE,FOCK(KGB),NTOTG)

            END IF

C
         END DO
C
  100 CONTINUE 
C
      RETURN
      END
C  /* Deck cc_dnspk */
      SUBROUTINE CC_DNSPK(DENSQ,DENSPK,ISYDEN)
C
C Purpose: construct lower triangular packed density matrix from
C          a full square density matrix
C
C written by Christof Haettig, July 1998
C
#include <implicit.h>
#include <maxorb.h>
#include <ccorb.h>
#include <symsq.h>
#include <ccsdsym.h>
C
      DIMENSION DENSQ(*), DENSPK(*)
C
C--------------------------------------------
C        construct triangular density matrix:
C--------------------------------------------
C
      CALL DZERO(DENSPK,NNBST(ISYDEN))
      DO ISYMB = 1, NSYM
        ISYMA = MULD2H(ISYDEN,ISYMB)
        DO B = 1, NBAS(ISYMB)
          DO A = 1, NBAS(ISYMA)
            KABSQ = IAODIS(ISYMA,ISYMB) + NBAS(ISYMA)*(B-1) + A
            KABPK = IADRPK( I2BST(ISYDEN) + KABSQ )
            DENSPK(KABPK) = DENSPK(KABPK) + DENSQ(KABSQ)
          END DO
        END DO
      END DO

      RETURN
      END 
