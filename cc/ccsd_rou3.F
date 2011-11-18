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
!  /* Deck cc_bf3 */
      SUBROUTINE CC_BF3(XINT,OMEGA2,XLAMD1,ISYML1,XLAMD2,
     *                  ISYML2,XLAMD3,ISYML3,
     *                  SCRM,ISYMM1,SCRM2,ISYMM2,WORK,LWORK,
     *                  IDEL,ISYMD,IOPT)
!
!     Written by Henrik Koch 3-Jan-1994
!     Symmetry by Henrik Koch and Alfredo Sanchez. 18-July-1994
!     Generalized by Asger Halkier and Henrik Koch 19/9 - 1995
!     to handle left-hand-side transformation contribution as well.
!     Righthand generalizations and debugging Ove Christiansen 23-9-1995
!
!     Ove Christiansen 24-9-1996: Generalization for calculating
!           terms similar to B and F-terms in the transformation
!           of vectors with the F-matrix.
!
!     Kasper Hald and Christof Haettig 22-2-1999
!     Generalized to calculate the BF-term for the triplet case.
!
!     Purpose: Calculate B-term and F-term in the orthonormal basis.
!
!     IOPT equals one for energy-calculations and two or three for
!     response calculations (2 for left trans. and 3 for right trans.)
!     IOPT eq. 4 for F*vector contributions.
!     IOPT equals 5  (Tilde)rhoBF(-)
!
!
!     XLAMD1 is always a true lamda matrix whereas XLAMD2
!     is an AO transformed trialvector in the case af a
!     response calculation.
!
!
!     24-9-1996:
!
!     IF (IOPT .EQ. 1):
!                       scrm is ordinary t2: tci,j(delta)
!                       XLAMD1 and XLAMD2 is ordinary lamda Hole
!                       matrices.
!                       (XLAMD1(gam,i)*XLAMD1(del,j))
!
!     IF (IOPT .EQ. 2/3)
!                       scrm is left/right vector transformed
!                       to tci,j(delta): vector general symmetry
!                       lambda particle/hole matrix is tot.sym.
!                       XLAMD1 is ordinary lambda particle/hole matrix.
!                       XLAMD2 is transformed (barred)
!                       lambda particle/hole matrix.
!                       (XLAMD1(gam,i)*XLAMD2(del,j)
!                       +XLAMD2(gam,i)*XLAMD1(del,j))
!
!     IF (IOPT .EQ. 4)
!                       scrm is left/right vector transformed
!                       to tci,j(delta): vector general symmetry
!                       lambda particle matrix is transformed.
!
!                       scrm2 is left/right vector transformed
!                       to tci,j(delta): vector general symmetry
!                       lambda particle matrix is tot.sym. ordinary
!                       lambda particle matrixes.
!
!                       Total transformed vector to be contracted
!                       with integrals is therefore
!
!                       XLAMD1 is an ordinary lambda particle matrix.
!                       XLAMD2 is a double transformed
!                              lambda particle matrix.
!                              (both R1 and L1)
!                       XLAMD3 is R1-transformed lambda particle matrix.
!
!                       sum(gam,del)(T(gam-bar,i,j,del)+T(gam,i,j,
!                                                            del-bar)
!                                   +3(XLAMD1(del,j)*XLAMD2(gam,i)
!                                     +XLAMD1(del,j)*XLAMD2(gam,i)))
!
!     22/2 - 1999
!
!   IF (IOPT .EQ. 5)  
!                     Triplet minus-intermediate.
!                     SCRM is right vector transformed
!                     to t(ci,j)(delta)
!                     Lambda particle/hole matrix is tot. sym.
!                     XLAMD1 is ordinary lambda particle/hole matr.
!                     XLAMD2 is the transformed lambda part./hole matr.
!                      XLAMD1(gam,i)*XLAMD2(del,j)
!                     +XLAMD2(gam,i)*XLAMD1(del,j)
!
!   IF (IOPT .EQ. 6)
!                     Same as IOPT .EQ. 3 except for the fact
!                     that the product of PLUS and PLUS
!                     (See the routine) is zero. (The T-amplitudes
!                     for the (+)triplet case are antisymmetric 
!                     with respect to the interchange of i and j )
!
!     The symmetry input to this routine is somewhat redundant but
!     hopefully logical and flexible:
!     Isymm1 is symmetry of SCRM
!     Isymm2 is symmetry of SCRM2
!     Isyml1 is symmetry of XLAMD1
!     Isyml2 is symmetry of XLAMD2
!     Isyml3 is symmetry of XLAMD3
!
      IMPLICIT NONE
!
      INTEGER LWORK, ISYMGD, ISYMM1, ISYML1, ISYML2
      INTEGER KMGD, KEND1, LWRK1, IDEL, ISYMD
      INTEGER ISYMJ, ISYMCI, ISYMI, ISYMC, ISYMG, ISYMGI, ISYMM2
      INTEGER NVIRC, NBASG
      INTEGER KOFF1, KOFF2, KOFF3, IOPT, ISYML3
!
#if defined (SYS_CRAY)
      REAL ZERO, HALF, ONE, TWO
      REAL XINT(*), OMEGA2(*), XLAMD1(*), XLAMD2(*), XLAMD3(*)
      REAL SCRM(*), SCRM2(*), WORK(LWORK)
#else
      DOUBLE PRECISION ZERO, HALF, ONE, TWO
      DOUBLE PRECISION XINT(*), OMEGA2(*), XLAMD1(*), XLAMD2(*)
      DOUBLE PRECISION XLAMD3(*), SCRM(*), SCRM2(*), WORK(LWORK)
#endif
!
      PARAMETER(ZERO= 0.0D00, HALF= 0.5D00, ONE= 1.0D00, TWO= 2.0D00)
!
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
!
      CALL QENTER('CC_BF3')
!
!------------------------
!     Dynamic allocation.
!------------------------
!
      ISYMGD = MULD2H(ISYMM1,ISYML1)
!
      KMGD   = 1
      KEND1  = KMGD   + NT2BGD(ISYMGD)
      LWRK1  = LWORK  - KEND1
!
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Need : ',KEND1,'Available : ',LWORK
         CALL QUIT('Insufficient space in CC_BF3')
      ENDIF
!
!-----------------------------
!     Prepare the data arrays.
!-----------------------------
!
      DO 100 ISYMJ = 1,NSYM
!
         ISYMCI = MULD2H(ISYMJ,ISYMM1)
!
         DO 110 ISYMI = 1,NSYM
!
            ISYMC  = MULD2H(ISYMI,ISYMCI)
            ISYMG  = MULD2H(ISYMC,ISYML1)
            ISYMGI = MULD2H(ISYMG,ISYMI)
!
            NVIRC = MAX(NVIR(ISYMC),1)
            NBASG = MAX(NBAS(ISYMG),1)
!
            KOFF1 = IGLMVI(ISYMG,ISYMC) + 1
!
            DO 120 J = 1,NRHF(ISYMJ)
!
               KOFF2 = IT2BCD(ISYMCI,ISYMJ) + IT1AM(ISYMC,ISYMI)
     *               + NT1AM(ISYMCI)*(J - 1) + 1
               KOFF3 = IT2BGD(ISYMGI,ISYMJ) + IT1AO(ISYMG,ISYMI)
     *               + NT1AO(ISYMGI)*(J - 1) + 1
!
               CALL DGEMM('N','N',NBAS(ISYMG),NRHF(ISYMI),NVIR(ISYMC),
     *                    ONE,XLAMD1(KOFF1),NBASG,SCRM(KOFF2),NVIRC,
     *                    ZERO,WORK(KOFF3),NBASG)
!
  120       CONTINUE
!
  110    CONTINUE
!
  100 CONTINUE
!
!---------------------------------------------------------
!     Calculate extra contribution to T2 double AO transf.
!     if F-matrix transformation.
!---------------------------------------------------------
!
      IF (IOPT .EQ. 4) THEN
!
         IF (MULD2H(ISYML3,ISYMM2).NE.ISYMGD) THEN
            CALL QUIT('CC_BF: Symmetry mismatch')
         ENDIF
         DO 200 ISYMJ = 1,NSYM
!
            ISYMCI = MULD2H(ISYMJ,ISYMM2)
!
            DO 210 ISYMI = 1,NSYM
!
               ISYMC  = MULD2H(ISYMI,ISYMCI)
               ISYMG  = MULD2H(ISYMC,ISYML3)
               ISYMGI = MULD2H(ISYMG,ISYMI)
!
               NVIRC = MAX(NVIR(ISYMC),1)
               NBASG = MAX(NBAS(ISYMG),1)
!
               KOFF1 = IGLMVI(ISYMG,ISYMC) + 1
!
                  DO 220 J = 1,NRHF(ISYMJ)
!
                  KOFF2 = IT2BCD(ISYMCI,ISYMJ) + IT1AM(ISYMC,ISYMI)
     *                  + NT1AM(ISYMCI)*(J - 1) + 1
                  KOFF3 = IT2BGD(ISYMGI,ISYMJ) + IT1AO(ISYMG,ISYMI)
     *                  + NT1AO(ISYMGI)*(J - 1) + 1
!
                  CALL DGEMM('N','N',NBAS(ISYMG),NRHF(ISYMI),NVIR(ISYMC)
     *                      ,ONE,XLAMD3(KOFF1),NBASG,SCRM2(KOFF2),NVIRC,
     *                       ONE,WORK(KOFF3),NBASG)
!
  220          CONTINUE
!
  210       CONTINUE
!
  200    CONTINUE
!
      ENDIF
!
!--------------------------------
!     Calculate the contribution.
!--------------------------------
!
      CALL CC_BF31(XINT,OMEGA2,WORK(KMGD),ISYMGD,XLAMD1,ISYML1,
     *             XLAMD2,ISYML2,WORK(KEND1),LWRK1,
     *             IDEL,ISYMD,IOPT)
!
      CALL QEXIT('CC_BF3')
!
      RETURN
      END
C  /* Deck cc_bf3_1 */
      SUBROUTINE CC_BF31(XINT,OMEGA2,XMGD,ISYMGD,XLAMD1,ISYML1,
     *                   XLAMD2,ISYML2,WORK,LWORK,
     *                   IDEL,ISYMD,IOPT)
!
!     Written by Henrik Koch 3-Jan-1994
!
!     Purpose: Calculate B-term.
!
!     See CC_BF for more info.
!
      IMPLICIT NONE
!
#include <priunit.h>
#include <iratdef.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
!
      INTEGER LWORK, INDEX, ISYDIS, ISYMD, ISYRES, ISYMGD, ISYCH
      INTEGER ISYML2, ISYML1, ISYMIJ, ISYMAB, ISYMG, IDEL, KSCRAB
      INTEGER KINDV1, KINDV2, KEND1, LWRK1, NSIZE, IMAXG, NMAXG
      INTEGER NBATCH, IBATCH, NUMG, IG1, IG2, KINTP, KINTM
      INTEGER KT2MP, KT2MM, KEND2, LWRK2, ISYMI, ISYMJ, ISYMGI
      INTEGER ISYMGJ, NTOTI, NGIJ, NGJI, NIJ, NGIJPM, KOFFP, KOFFM
      INTEGER IOPT, KOFF1, KOFF2, ISHELP, KOFF, NUMGM, NTOTAB
!
#if defined (SYS_CRAY)
      REAL ZERO, HALF, ONE, FOURTH, TWO, THREE
      REAL XTWO, XHALF, XONE, FACT
      REAL XINT(*), OMEGA2(*), XMGD(*), XLAMD1(*), XLAMD2(*)
      REAL WORK(LWORK)
#else
      DOUBLE PRECISION ZERO, HALF, ONE, FOURTH, TWO
      DOUBLE PRECISION THREE, XTWO, XHALF, XONE, FACT
      DOUBLE PRECISION XINT(*), OMEGA2(*), XMGD(*), XLAMD1(*)
      DOUBLE PRECISION XLAMD2(*), WORK(LWORK)
#endif
      PARAMETER(ZERO = 0.0D00, HALF = 0.5D00, ONE = 1.0D00)
      PARAMETER(FOURTH = 0.25D00, TWO = 2.0D00, THREE = 3.0D00)
      PARAMETER(XTWO = -2.0D00, XHALF= -0.5D00, XONE= -1.0D00)
!
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
!
      CALL QENTER('CC_BF31')
!
      ISYDIS = MULD2H(ISYMOP,ISYMD)
      ISYRES = MULD2H(ISYDIS,ISYMGD)
      ISYCH  = MULD2H(ISYML2,ISYMD)
!
      IF (ISYML1 .NE. 1)
     &     CALL QUIT('CC_BF3: Symmetry of XLAMD1 must be 1')
      IF (ISYML2 .NE. MULD2H(ISYMGD,ISYMD))
     *            CALL QUIT('Symmetry mismatch in CC_BF3_1')
!
!================================
!     Calculate the contribution.
!================================
!
      DO 100 ISYMIJ = 1 , NSYM
C
         ISYMAB = MULD2H(ISYMIJ,ISYRES)
         ISYMG  = MULD2H(ISYMAB,ISYDIS)
         D      = IDEL - IBAS(ISYMD)
C
         KSCRAB = 1
         KINDV1 = KSCRAB + N2BST(ISYMAB)
         KINDV2 = KINDV1 + (NNBST(ISYMAB) - 1)/IRAT + 1
         KEND1  = KINDV2 + (NNBST(ISYMAB) - 1)/IRAT + 1
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient space in CC_BF3_1')
         ENDIF
C
C--------------------------------
C        Calculate index vectors.
C--------------------------------
C
         CALL CCSD_INDEX(WORK(KINDV1),WORK(KINDV2),ISYMAB)
C
C------------------------------
C        Work space allocation.
C------------------------------
C
         NSIZE  = 2*(NNBST(ISYMAB) + NMIJP(ISYMIJ))
C
         IF ((NNBST(ISYMAB) .EQ. 0) .OR.
     *       (NMIJP(ISYMIJ) .EQ. 0)) GOTO 100
C
         IF (ISYMG .EQ. ISYMD) THEN
            IMAXG = D
         ELSE IF (ISYMG .LT. ISYMD) THEN
            IMAXG = NBAS(ISYMG)
         ELSE
            GOTO 100
         ENDIF
C
         IF (IMAXG.EQ.0) GOTO 100
C
         IF (LWRK1.LT.NSIZE) THEN
           CALL QUIT('Insufficient memory in CC_BF1.')
         END IF
C
         NMAXG  = MIN(IMAXG,LWRK1/NSIZE)
         NBATCH = (IMAXG - 1)/NMAXG + 1
C
         DO 110 IBATCH = 1 , NBATCH
C
            NUMG = NMAXG
            IF (IBATCH .EQ. NBATCH) THEN
               NUMG = IMAXG - NMAXG*(NBATCH - 1)
            ENDIF
C
            IG1 = NMAXG*(IBATCH - 1) + 1
            IG2 = NMAXG*(IBATCH - 1) + NUMG
C
            KINTP = KEND1
            KINTM = KINTP + NNBST(ISYMAB)*NUMG
            KT2MP = KINTM + NNBST(ISYMAB)*NUMG
            KT2MM = KT2MP + NUMG*NMIJP(ISYMIJ)
            KEND2 = KT2MM + NUMG*NMIJP(ISYMIJ)
            LWRK2 = LWORK - KEND2
C
            IF (LWRK2 .LT. 0) THEN
               CALL QUIT('Insufficient space in CC_BF31')
            ENDIF
C
C-----------------------------------
C           Construct T2MP and T2MM.
C-----------------------------------
C
            DO 200 ISYMJ = 1 , NSYM
C
               ISYMI  = MULD2H(ISYMJ,ISYMIJ)
               ISYMGI = MULD2H(ISYMI,ISYMG)
               ISYMGJ = MULD2H(ISYMJ,ISYMG)
C
               IF (ISYMI .GT. ISYMJ) GOTO 200
C
               NTOTI = NRHF(ISYMI)
C
               DO 210 J = 1 , NRHF(ISYMJ)
C
                  IF (ISYMI .EQ. ISYMJ) NTOTI = J
C
                  DO 220 I = 1,NTOTI
C
                     NGIJ = IT2BGD(ISYMGI,ISYMJ)
     *                    + NT1AO(ISYMGI)*(J - 1)
     *                    + IT1AO(ISYMG,ISYMI)
     *                    + NBAS(ISYMG)*(I - 1) + IG1
C
                     NGJI = IT2BGD(ISYMGJ,ISYMI)
     *                    + NT1AO(ISYMGJ)*(I - 1)
     *                    + IT1AO(ISYMG,ISYMJ)
     *                    + NBAS(ISYMG)*(J - 1) + IG1
C
                     IF (ISYMI .EQ. ISYMJ) THEN
                        NIJ = IMIJP(ISYMI,ISYMJ) + INDEX(I,J)
                     ELSE
                        NIJ = IMIJP(ISYMI,ISYMJ)
     *                      + NRHF(ISYMI)*(J - 1) + I
                     ENDIF
C
                     NGIJPM = NUMG*(NIJ - 1)
C
                     KOFFP = KT2MP + NGIJPM
                     KOFFM = KT2MM + NGIJPM
C
C
                        IF (IOPT .NE. 6) THEN
                           CALL DCOPY(NUMG,XMGD(NGIJ),1,WORK(KOFFP),1)
                        ENDIF
!
                        CALL DCOPY(NUMG,XMGD(NGIJ),1,WORK(KOFFM),1)
C
                        IF (IOPT .NE. 6) THEN
                         CALL DAXPY(NUMG,ONE,XMGD(NGJI),1,WORK(KOFFP),1)
                        ENDIF
!
                        CALL DAXPY(NUMG,-ONE,XMGD(NGJI),1,WORK(KOFFM),1)
!
C
C-------------------------------------------------
C                    Add the F-term contributions.
C-------------------------------------------------
C
                     FACT = ONE
C
                     IF ((IOPT .EQ. 2) .OR. (IOPT .EQ. 4)) THEN
                        FACT = THREE
                     ENDIF
C
                     IF ((ISYMJ .EQ. ISYCH).AND.(ISYMI .EQ. ISYMG)) THEN
C
                        KOFF1 = IGLMRH(ISYMD,ISYMJ)
     &                        + NBAS(ISYMD)*(J - 1) + D
                        KOFF2 = ILMRHF(ISYMI) + NBAS(ISYMG)*(I - 1) +IG1
C
                        IF (IOPT .EQ. 5) THEN
!
                           CALL DAXPY(NUMG,XHALF*XLAMD2(KOFF1),
     *                                XLAMD1(KOFF2),1,WORK(KOFFP),1)
                           CALL DAXPY(NUMG,XHALF*XLAMD2(KOFF1),
     &                                XLAMD1(KOFF2),1,WORK(KOFFM),1)
                        ELSE
!
                           IF (IOPT .NE. 6) THEN
                            CALL DAXPY(NUMG,XLAMD2(KOFF1),XLAMD1(KOFF2),
     &                                 1,WORK(KOFFP),1)
                           ENDIF
!
                           CALL DAXPY(NUMG,FACT*XLAMD2(KOFF1),
     &                                XLAMD1(KOFF2),1,WORK(KOFFM),1)
!
                        ENDIF
C
                     ENDIF
C
                     IF ((ISYMI .EQ. ISYCH).AND.(ISYMJ .EQ. ISYMG)) THEN
C
                        KOFF1 = IGLMRH(ISYMD,ISYMI)
     &                        + NBAS(ISYMD)*(I - 1) + D
                        KOFF2 = ILMRHF(ISYMJ) + NBAS(ISYMG)*(J - 1) +IG1
C
                        IF (IOPT .EQ. 5) THEN
!
                           CALL DAXPY(NUMG,XHALF*XLAMD2(KOFF1),
     *                                XLAMD1(KOFF2),1,WORK(KOFFP),1)
                           CALL DAXPY(NUMG,HALF*XLAMD2(KOFF1),
     *                                XLAMD1(KOFF2),1,WORK(KOFFM),1)
!
                        ELSE
!
                           IF (IOPT .NE. 6) THEN
                            CALL DAXPY(NUMG,XLAMD2(KOFF1),XLAMD1(KOFF2),
     *                                 1,WORK(KOFFP),1)
                           ENDIF
!
                           CALL DAXPY(NUMG,-FACT*XLAMD2(KOFF1),
     *                                XLAMD1(KOFF2),1,WORK(KOFFM),1)
C
                        ENDIF
!
                     ENDIF
C
C---------------------------------------------------------------------
C                    For response calculation add permuted terms.
C---------------------------------------------------------------------
C
                     IF (IOPT .GE. 2) THEN
C
                        ISHELP = MULD2H(ISYMG,ISYML2)
C
                        IF ((IOPT .EQ. 2) .OR. (IOPT .EQ. 4)) THEN
                           FACT = THREE
                        ENDIF
C
                        IF ((ISYMJ .EQ. ISYMD) .AND.
     &                      (ISYMI .EQ. ISHELP)) THEN
C
                           KOFF1 = ILMRHF(ISYMJ)
     &                           + NBAS(ISYMD)*(J - 1) + D
                           KOFF2 = IGLMRH(ISYMG,ISYMI)
     &                           + NBAS(ISYMG)*(I - 1) +IG1
C
                           IF (IOPT .EQ. 5) THEN
C
                              CALL DAXPY(NUMG,HALF*XLAMD1(KOFF1),
     &                                XLAMD2(KOFF2),1,WORK(KOFFP),1)
                              CALL DAXPY(NUMG,HALF*XLAMD1(KOFF1),
     &                                XLAMD2(KOFF2),1,WORK(KOFFM),1)
C
                           ELSE
C
                              IF (IOPT .NE. 6) THEN
                                 CALL DAXPY(NUMG,XLAMD1(KOFF1),
     &                                   XLAMD2(KOFF2),1,WORK(KOFFP),1)
                              ENDIF
!
                              CALL DAXPY(NUMG,FACT*XLAMD1(KOFF1),
     &                                XLAMD2(KOFF2),1,WORK(KOFFM),1)
C
                           ENDIF
C
                        ENDIF
C
                        IF ((ISYMI .EQ. ISYMD) .AND.
     &                      (ISYMJ .EQ. ISHELP)) THEN
C
                           KOFF1 = ILMRHF(ISYMI)
     &                           + NBAS(ISYMD)*(I - 1) + D
                           KOFF2 = IGLMRH(ISYMG,ISYMJ)
     &                           + NBAS(ISYMG)*(J - 1) + IG1
C
                           IF (IOPT .EQ. 5) THEN
C
                              CALL DAXPY(NUMG,HALF*XLAMD1(KOFF1),
     &                                XLAMD2(KOFF2),1,WORK(KOFFP),1)
                              CALL DAXPY(NUMG,XHALF*XLAMD1(KOFF1),
     &                                XLAMD2(KOFF2),1,WORK(KOFFM),1)
C
                           ELSE
C
                              IF (IOPT .NE. 6) THEN
                                 CALL DAXPY(NUMG,XLAMD1(KOFF1),
     &                                   XLAMD2(KOFF2),1,WORK(KOFFP),1)
                              ENDIF
!
                              CALL DAXPY(NUMG,-FACT*XLAMD1(KOFF1),
     &                                XLAMD2(KOFF2),1,WORK(KOFFM),1)
C
                           ENDIF
C
                        ENDIF
C
                     ENDIF
C
  220          CONTINUE
C
  210       CONTINUE
C
  200     CONTINUE
C
C-----------------------------------
C           Construct INTP and INTM.
C-----------------------------------
C
            CALL CCRHS_IPM(XINT,WORK(KINTP),WORK(KINTM),WORK(KSCRAB),
     *                     WORK(KINDV1),WORK(KINDV2),ISYMAB,ISYMG,
     *                     NUMG,IG1,IG2)
C
C-------------------------------
C           Scale the diagonals.
C-------------------------------
C
            IF ((ISYMG .EQ. ISYMD) .AND. (IBATCH .EQ. NBATCH)) THEN
               KOFF = KINTP + NNBST(ISYMAB)*(NUMG - 1)
               CALL DSCAL(NNBST(ISYMAB),HALF,WORK(KOFF),1)
            ENDIF
C
C----------------------------------------
C           Add the B-term contributions.
C----------------------------------------
C
            NUMGM  = MAX(NUMG,1)
            NTOTAB = MAX(NNBST(ISYMAB),1)
               
C
            IF (.NOT. (IOPT .EQ. 5)) THEN
!
              IF (IOPT .NE. 6) THEN
                KOFF1 = IT2ORT(ISYMAB,ISYMIJ) + 1
!
                CALL DGEMM('N','N',NNBST(ISYMAB),NMIJP(ISYMIJ),NUMG,
     *                     ONE,WORK(KINTP),NTOTAB,WORK(KT2MP),NUMGM,
     *                     ONE,OMEGA2(KOFF1),NTOTAB)
              ENDIF
!
              IF (IOPT .NE. 6) THEN
                KOFF2 = NT2ORT(ISYRES) + IT2ORT(ISYMAB,ISYMIJ) + 1
              ELSE
                KOFF2 =IT2ORT(ISYMAB,ISYMIJ) + 1
              ENDIF
!
              CALL DGEMM('N','N',NNBST(ISYMAB),NMIJP(ISYMIJ),NUMG,
     *                   ONE,WORK(KINTM),NTOTAB,WORK(KT2MM),NUMGM,
     *                   ONE,OMEGA2(KOFF2),NTOTAB)
!
            ELSE
            KOFF1 = IT2ORT(ISYMAB,ISYMIJ) + 1
!
            CALL DGEMM('N','N',NNBST(ISYMAB),NMIJP(ISYMIJ),NUMG,
     *                 ONE,WORK(KINTM),NTOTAB,WORK(KT2MP),NUMGM,
     *                 ONE,OMEGA2(KOFF1),NTOTAB)
!
            KOFF2 = NT2ORT(ISYRES) + IT2ORT(ISYMAB,ISYMIJ) + 1
!
            CALL DGEMM('N','N',NNBST(ISYMAB),NMIJP(ISYMIJ),NUMG,
     *                 ONE,WORK(KINTP),NTOTAB,WORK(KT2MM),NUMGM,
     *                 ONE,OMEGA2(KOFF2),NTOTAB)
!
            END IF
  110    CONTINUE
!
  100 CONTINUE
!
      CALL QEXIT('CC_BF31')
!
      RETURN
      END
C  /* Deck ccrhs_d3 */
      SUBROUTINE CCRHS_D3(XINT,DSRHF,OMEGA2,T2AM,ISYMT2,
     *                    XLAMDP,XLAMIP,XLAMDH,
     *                    XLAMPC,ISYMPC,XLAMHC,ISYMHC,
     *                    SCRM,WORK,LWORK,IDEL,ISYMD,FACTD,ICON,
     *                    LUD,DFIL,IV)
!
!     Written by Henrik Koch 9-Jan-1994
!
!     Generalisation for CCLR by Ove Christiansen august-september 1995
!     (right transformation) and september 1996 (F-matrix).
!
!     Generalisation to calculate the D-intermediates for the 
!     triplet case by Kasper Hald 17-2-1999
!
!     Purpose: Calculate D-term.
!
      IMPLICIT NONE
!
      INTEGER LWORK, ISYDIS, ISYAIK, ISYMPC, ISYMT2, KSCR1, KSCR2
      INTEGER KSCR3, KEND1, LWRK1, ISYMD, INDEX, ISYMHC, LUD, IV 
      INTEGER KOFF1, ISYML, ISYMA, ISYMG
      INTEGER NBASA, NBASG, NVIRD, KSCR11, KEND2, LWRK2, KOFF2
      INTEGER KOFF3, KOFF5, KOFF6, NRHFK, ISYMAI, NTOTDL
      INTEGER IOFF, IDEL, IERR, ICON
!
#if defined (SYS_CRAY)
      REAL ONE, TWO, FACTD
      REAL XINT(*), DSRHF(*), OMEGA2(*), WORK(LWORK)
      REAL XLAMDP(*), XLAMIP(*), XLAMDH(*), SCRM(*)
      REAL XLAMPC(*), XLAMHC(*), T2AM(*)
#else
      DOUBLE PRECISION ONE, TWO, FACTD
      DOUBLE PRECISION XINT(*), DSRHF(*), OMEGA2(*), WORK(LWORK)
      DOUBLE PRECISION XLAMDP(*), XLAMIP(*), XLAMDH(*), SCRM(*)
      DOUBLE PRECISION XLAMPC(*), XLAMHC(*), T2AM(*)
#endif
!
      PARAMETER (ONE = 1.0D00, TWO = 2.0D00)
      CHARACTER DFIL*(*)
!
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
!
      CALL QENTER('CCRHS_D3')
!
      ISYDIS = MULD2H(ISYMD,ISYMOP)
      ISYAIK = MULD2H(ISYDIS,ISYMPC)
      IF (ISYMT2 .NE. ISYMPC) CALL QUIT('Symmetry Mismatch in CCRHS_D3')
C
C------------------------
C     Dynamic allocation.
C------------------------
C
      KSCR1  = 1
      KSCR2  = KSCR1  + MAX(NT2BCD(ISYAIK),NT2BCD(ISYDIS))
      KSCR3  = KSCR2  + NT2BGD(ISYDIS)
      IF (ICON .EQ. 2) THEN
         KEND1  = KSCR3  + NT2BGD(ISYMD)
      ELSE IF (ICON .EQ.5) THEN
         KEND1  = KSCR3 + MAX(NT2BGD(ISYMD),NT2BCD(ISYAIK),
     *                        NT2BCD(ISYDIS))
      ELSE
         KEND1  = KSCR3  + NT2BGD(ISYAIK)
      ENDIF
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Need : ',KEND1,'Available : ',LWORK
         CALL QUIT('Insufficient space in CCRHS_D3')
      ENDIF
C
C--------------------------------
C     Calculate the contribution.
C--------------------------------
C
      CALL CCRHS_D31(XINT,DSRHF,OMEGA2,T2AM,ISYMT2,
     *                SCRM,WORK(KSCR1),WORK(KSCR2),
     *                WORK(KSCR3),XLAMDP,XLAMIP,
     *                XLAMDH,XLAMPC,ISYMPC,XLAMHC,ISYMHC,
     *                WORK(KEND1),LWRK1,ISYDIS,IDEL,
     *                ISYMD,FACTD,ICON,LUD,DFIL,IV)
C
      CALL QEXIT('CCRHS_D3')
C
      RETURN
      END
      SUBROUTINE CCRHS_D31(XINT,DSRHF,OMEGA2,T2AM,ISYMT2,
     *                      SCRM,SCR1,SCR2,SCR3,
     *                      XLAMDP,XLAMIP,XLAMDH,XLAMPC,ISYMPC,XLAMHC,
     *                      ISYMHC,WORK,LWORK,ISYDIS,IDEL,ISYDEL,FACTD,
     *                      ICON,LUD,DFIL,IV)
C
C     Written by Henrik Koch 3-Jan-1994
C
C     Modification by Ove Christiansen 25-7-1995 to allow for a
C     general factor (FACTD). NB: Assumes DUMCD.
C     - calculate intermediates for CCLR.
C
C     29-9-1995 (17-9-1996 F-matrix.) Ove Christiansen:
C
C                 1. If Icon = 2 both contributions are calculated,
C                    for total sym. case. Otherwise
C                    a.ICON = 1 only the integral Laikc(del)
C                               = La-bar,i,k,c + La,i-bar,k,c
C                      for Jacobian right transformation
C                    b.ICON = 3
C                          La-bar,i,k,c + La,i-bar,k,c + Tx*Int
C                      for F-matrix times vector.
C
C                 2. Allow for general transformation matrix for
C                    alpha to a(XLAMPC) and for i (XLAMHC).
C                    (the extra i transformation introduces new
C                     blocks which is only entered when icon = 1 or 3)
C
C                 3. If icon diff. from 2 (we have linear response)
C                    The D intermediate is stored according to
C                    the number of simultaneous trial vector
C                    given by IV. This is ensured using IT2DLR.
C
!     17-2-1999 Kasper Hald:
!
!     IF ICON = 4 then the triplet intermediate:
!     
!     g(a-bar,i,l,c) + g(a,i-bar,l,c) is calculated
!
!     IF ICON = 5 then the triplet intermediate:
!
!     g(aikc) + sum(dl)t(ai,dl)L(kcld) - sum(dl)t(di,al)g(ldkc)
!
!     ICON 6: g(a-bar,i,l,c) - g(a,i-bar,l,c)
!
!     ICON 4, ICON 5 and ICON 6 assumes DUMPCD
!
!     Purpose: Calculate D-term.
!
      IMPLICIT NONE
!
#include <priunit.h>
#include <maxorb.h>
#include <ccorb.h>
#include <symsq.h>
#include <ccsdsym.h>
#include <ccsdio.h>
!
      INTEGER LWORK, ICON, ISYMK, ISYMAG, ISYMDL, KSCR10, KEND1
      INTEGER LWRK1, KOFF1, ISYML, ISYMD, ISYMA, ISYMG
      INTEGER NBASA, NBASG, NVIRD, KSCR11, KEND2, LWRK2, KOFF2
      INTEGER KOFF3, KOFF5, KOFF6, INDEX, ISYAIK, ISYDIS, ISYMPC
      INTEGER NRHFK, ISYMAI, NTOTDL, IOFF, IERR, ISYMBG
      INTEGER ISYMI, ISYMB, NBASB, KSCR12, NAI, KOFF7, KOFF8
      INTEGER ISYMHC, ISALIK, ISYALG, ISYALI, NT1AOM, ISYMAL
      INTEGER NBASAL, KOFF4, MAI, ISYMBJ, ISYDEL, ISYMJ, NVIRB
      INTEGER NTOTBJ, NBJ, NAIBJ, MALI, IV, IDEL, ISYM5, IOPT5
      INTEGER ISYMT2, LUD
!
#if defined (SYS_CRAY)
      REAL ZERO, ONE, HALF, XMHALF, TWO, XMONE, FACTD
      REAL XINT(*), OMEGA2(*), T2AM(*), DSRHF(*), SCRM(*)
      REAL SCR1(*), SCR2(*), SCR3(*), XLAMDP(*), XLAMIP(*), XLAMDH(*)
      REAL XLAMPC(*), XLAMHC(*), WORK(LWORK)
#else
      DOUBLE PRECISION ZERO, ONE, HALF, XMHALF, TWO, XMONE, FACTD
      DOUBLE PRECISION XINT(*), OMEGA2(*), T2AM(*), DSRHF(*), SCRM(*)
      DOUBLE PRECISION SCR1(*), SCR2(*), SCR3(*), XLAMDP(*), XLAMIP(*)
      DOUBLE PRECISION XLAMDH(*), XLAMPC(*), XLAMHC(*), WORK(LWORK)
#endif
      PARAMETER(ZERO=0.0D00,ONE=1.0D00,HALF=0.5D00,XMHALF=-0.5D00)
      PARAMETER(TWO=2.0D00, XMONE=-1.0D00)
      CHARACTER DFIL*(*)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      CALL QENTER('CCRHS_D31')
C
      ISYAIK = MULD2H(ISYDIS,ISYMPC)
C
C-------------------------------------------------------
C     Calculate the integrals K(k,dl) = (k d | l delta).
C-------------------------------------------------------
C
      IF ((ICON .EQ. 2) .OR. (ICON .EQ. 3) .OR. (ICON .EQ. 5)) THEN
C
         DO 100 ISYMK = 1,NSYM
C
            ISYMAG = MULD2H(ISYMK,ISYDIS)
C
            DO 110 K = 1,NRHF(ISYMK)
C
               ISYMDL = MULD2H(ISYMK,ISYDIS)
C
               KSCR10 = 1
               KEND1  = KSCR10 + N2BST(ISYMAG)
               LWRK1  = LWORK  - KEND1
C
               IF (LWRK1 .LT. 0) THEN
                  CALL QUIT('Not enough space for '//
     &                 'allocation in CCRHS_D31(1)')
               END IF
C
               KOFF1 = IDSRHF(ISYMAG,ISYMK) + NNBST(ISYMAG)*(K-1) + 1
               CALL CCSD_SYMSQ(DSRHF(KOFF1),ISYMAG,WORK(KSCR10))
C
               DO 120 ISYML = 1,NSYM
C
                  ISYMD = MULD2H(ISYML,ISYMDL)
                  ISYMA = ISYML
                  ISYMG = ISYMD
C
                  NBASA = MAX(NBAS(ISYMA),1)
                  NBASG = MAX(NBAS(ISYMG),1)
                  NVIRD = MAX(NVIR(ISYMD),1)
C
                  KSCR11 = KEND1
                  KEND2  = KSCR11 + NBAS(ISYMG)*NRHF(ISYML)
                  LWRK2  = LWORK  - KEND2
C
                  IF (LWRK2 .LT. 0) THEN
                     CALL QUIT('Not enough space for '//
     &                    'allocation in CCRHS_D31')
                  END IF
C
                  KOFF2 = KSCR10 + IAODIS(ISYMA,ISYMG)
                  KOFF3 = ILMRHF(ISYML) + 1
C
                  CALL DGEMM('T','N',NBAS(ISYMG),NRHF(ISYML),
     *                       NBAS(ISYMA),ONE,WORK(KOFF2),NBASA,
     *                       XLAMDP(KOFF3),NBASA,
     *                       ZERO,WORK(KSCR11),NBASG)
C
                  KOFF5 = ILMVIR(ISYMD) + 1
                  KOFF6 = IT2BCD(ISYMDL,ISYMK) + NT1AM(ISYMDL)*(K - 1)
     *                  + IT1AM(ISYMD,ISYML) + 1
C
                  CALL DGEMM('T','N',NVIR(ISYMD),NRHF(ISYML),
     *                       NBAS(ISYMG),ONE,XLAMDH(KOFF5),NBASG,
     *                       WORK(KSCR11),NBASG,
     *                       ZERO,SCR1(KOFF6),NVIRD)
C
  120          CONTINUE
C
  110       CONTINUE
C
  100    CONTINUE
!
!----------------------------------------------------
!      For ICON = 5 calculate the last part
!      of the intermediate. (t2(di,al)*g(ldkc))
!----------------------------------------------------
!
         IF (ICON .EQ. 5) THEN
!
!---------------------------
!      Transpose T2.
!---------------------------
!
           ISYM5 = ISYMT2
           CALL CCSD_T2TP(T2AM,WORK,LWORK,ISYM5)
!
           IF (LWORK .LT. NT2BCD(ISYDIS)) THEN
              CALL QUIT('Not enough space in CCRHS_D3 (IOPT = 5)')
           ENDIF
!
!-----------------------------
!      Calculate the cont.
!-----------------------------
!
           DO 123 ISYMK = 1,NSYM
!
              ISYMDL = MULD2H(ISYMK,ISYDIS)
!
              NRHFK = MAX(NRHF(ISYMK),1)
!
              DO 126 K = 1,NRHF(ISYMK)
!
                 KOFF1 = IT2BCD(ISYMDL,ISYMK)+NT1AM(ISYMDL)*(K-1)+1
                 KOFF2 = IT2BCT(ISYMK,ISYMDL) + K
!
                 CALL DCOPY(NT1AM(ISYMDL),SCR1(KOFF1),1,WORK(KOFF2),
     *                      NRHFK)
!
  126         CONTINUE
!
  123      CONTINUE
!
           CALL DCOPY(NT2BCD(ISYDIS),WORK,1,SCR3,1)
!
           IF (LWORK .LT. NT2BCD(ISYAIK)) THEN
              CALL QUIT('Insufficient work space in CCRHS_D31')
           ENDIF
!
           DO ISYMK = 1,NSYM
!
              ISYMDL = MULD2H(ISYMK,ISYDIS)
              ISYMAI = MULD2H(ISYAIK,ISYMK)
!
              NRHFK  = MAX(NRHF(ISYMK),1)
              NTOTDL = MAX(NT1AM(ISYMDL),1)
!
              KOFF1  = IT2BCT(ISYMK,ISYMDL) + 1
              KOFF2  = IT2SQ(ISYMDL,ISYMAI) + 1
              KOFF3  = IT2BCT(ISYMK,ISYMAI) + 1
!
           CALL DGEMM('N','N',NRHF(ISYMK),NT1AM(ISYMAI),NT1AM(ISYMDL),
     *                ONE,SCR3(KOFF1),NRHFK,T2AM(KOFF2),NTOTDL,ZERO,
     *                WORK(KOFF3),NRHFK)
!
           ENDDO
!
           CALL DCOPY(NT2BCD(ISYAIK),WORK,1,SCR3,1)
!
!----------------------------------------
!     Transpose T2 (back).
!----------------------------------------
!
           ISYM5 = ISYMT2
           CALL CCSD_T2TP(T2AM,WORK,LWORK,ISYM5)
!
         ENDIF
!
!--------------------------------------------
!        Transpose integral array.
!--------------------------------------------
!
         CALL CC_MTCME(SCR1,WORK,LWORK,ISYDIS,1)
C
         IF (LWORK .LT. NT2BCD(ISYDIS)) THEN
            CALL QUIT('Not enough space for allocation in CCRHS_D31(3)')
         END IF
C
         DO 130 ISYMK = 1,NSYM
C
            ISYMDL = MULD2H(ISYMK,ISYDIS)
C
            NRHFK = MAX(NRHF(ISYMK),1)
C
            DO 140 K = 1,NRHF(ISYMK)
C
               KOFF1 = IT2BCD(ISYMDL,ISYMK) + NT1AM(ISYMDL)*(K - 1) + 1
               KOFF2 = IT2BCT(ISYMK,ISYMDL) + K
C
               CALL DCOPY(NT1AM(ISYMDL),SCR1(KOFF1),1,WORK(KOFF2),NRHFK)
C
  140       CONTINUE
C
  130    CONTINUE
C
         CALL DCOPY(NT2BCD(ISYDIS),WORK,1,SCR1,1)
C
C-----------------------------------------
C        Calculate the first contribution.
C        sum(2*t(ai,dl)-t(di,al))*L(ldkc)
C-----------------------------------------
C
         IF (LWORK .LT. NT2BCD(ISYAIK)) THEN
            CALL QUIT('Insufficient work space in CCRHS_D31')
         ENDIF
C
         DO 200 ISYMK = 1,NSYM
C
            ISYMDL = MULD2H(ISYMK,ISYDIS)
            ISYMAI = MULD2H(ISYAIK,ISYMK)
C
            NRHFK  = MAX(NRHF(ISYMK),1)
            NTOTDL = MAX(NT1AM(ISYMDL),1)
C
            KOFF1  = IT2BCT(ISYMK,ISYMDL) + 1
            KOFF2  = IT2SQ(ISYMDL,ISYMAI) + 1
            KOFF3  = IT2BCT(ISYMK,ISYMAI) + 1
C
            CALL DGEMM('N','N',NRHF(ISYMK),NT1AM(ISYMAI),NT1AM(ISYMDL),
     *                 ONE,SCR1(KOFF1),NRHFK,T2AM(KOFF2),NTOTDL,ZERO,
     *                 WORK(KOFF3),NRHFK)
C
  200    CONTINUE
C
         CALL DCOPY(NT2BCD(ISYAIK),WORK,1,SCR1,1)
!
      ENDIF
!
      IF (ICON .EQ. 5) THEN
!
         CALL DAXPY(NT2BCD(ISYAIK),XMONE,SCR3,1,SCR1,1)
!
      ENDIF
!
C----------------------------------------------------------
C     Calculate the integrals K(k,ai) = (k i | alfa delta).
C----------------------------------------------------------
C
      DO 300 ISYMA = 1,NSYM
C
         ISYMBG = MULD2H(ISYMA,ISYDIS)
C
         KSCR10 = 1
         KEND1  = KSCR10 + N2BST(ISYMBG)
         LWRK1  = LWORK  - KEND1
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Not enough space for allocation in CCRHS_D31')
         END IF
C
         DO 310 A = 1,NBAS(ISYMA)
C
            KOFF1 = IDSAOG(ISYMA,ISYDIS) + NNBST(ISYMBG)*(A - 1) + 1
            CALL CCSD_SYMSQ(XINT(KOFF1),ISYMBG,WORK(KSCR10))
C
            DO 320 ISYMG = 1,NSYM
C
               ISYMI  = ISYMG
               ISYMB  = MULD2H(ISYMG,ISYMBG)
               ISYMK  = ISYMB
               ISYMAI = MULD2H(ISYMA,ISYMI)
C
               NBASB = MAX(NBAS(ISYMB),1)
               NBASG = MAX(NBAS(ISYMG),1)
               NRHFK = MAX(NRHF(ISYMK),1)
C
               KSCR11 = KEND1
               KSCR12 = KSCR11 + NRHF(ISYMK)*NBAS(ISYMG)
               KEND2  = KSCR12 + NRHF(ISYMK)*NRHF(ISYMI)
               LWRK2  = LWORK  - KEND2
               IF (LWRK2 .LT. 0) THEN
                  CALL QUIT('Not enough space for '//
     &                 'allocation in CCRHS_D1')
               END IF
C
               KOFF2 = ILMRHF(ISYMK) + 1
               KOFF3 = KSCR10 + IAODIS(ISYMB,ISYMG)
C
               CALL DGEMM('T','N',NRHF(ISYMK),NBAS(ISYMG),NBAS(ISYMB),
     *                    ONE,XLAMDP(KOFF2),NBASB,WORK(KOFF3),NBASB,
     *                    ZERO,WORK(KSCR11),NRHFK)
C
               KOFF5 = ILMRHF(ISYMI) + 1
C
               CALL DGEMM('N','N',NRHF(ISYMK),NRHF(ISYMI),NBAS(ISYMG),
     *                    ONE,WORK(KSCR11),NRHFK,XLAMDH(KOFF5),NBASG,
     *                    ZERO,WORK(KSCR12),NRHFK)
C
               DO 330 I = 1,NRHF(ISYMI)
C
                  NAI = IT1AO(ISYMA,ISYMI) + NBAS(ISYMA)*(I-1) + A
C
                  KOFF8 = IT2BGT(ISYMK,ISYMAI)
     *                  + NRHF(ISYMK)*(NAI - 1) + 1
                  KOFF7 = KSCR12 + NRHF(ISYMK)*(I - 1)
C
                  CALL DCOPY(NRHF(ISYMK),WORK(KOFF7),1,SCR2(KOFF8),1)
C
  330          CONTINUE
C
C-------------------------------------------------------
C              In 2C1 linear transformation extra  cont.
C-------------------------------------------------------
C
               IF ((ICON .EQ. 1) .OR. (ICON.EQ.3)) THEN
C
                  ISYMI  = MULD2H(ISYMG,ISYMHC)
                  ISYMAI = MULD2H(ISYMA,ISYMI)
C
                  KEND2  = KSCR12 + NRHF(ISYMK)*NRHF(ISYMI)
                  LWRK2  = LWORK  - KEND2
                  IF (LWRK2 .LT. 0) THEN
                     CALL QUIT('Not enough space for '//
     &                    'allocation in CCRHS_D1')
                  END IF
C
                  KOFF5 = IGLMRH(ISYMG,ISYMI) + 1
C
                  CALL DGEMM('N','N',NRHF(ISYMK),NRHF(ISYMI),
     *                       NBAS(ISYMG),ONE,WORK(KSCR11),NRHFK,
     *                       XLAMHC(KOFF5),NBASG,
     *                       ZERO,WORK(KSCR12),NRHFK)
C
                  DO 331 I = 1,NRHF(ISYMI)
C
                     NAI = IT1AO(ISYMA,ISYMI) + NBAS(ISYMA)*(I-1) + A
C
                     KOFF8 = IT2BGT(ISYMK,ISYMAI)
     *                     + NRHF(ISYMK)*(NAI - 1) + 1
                     KOFF7 = KSCR12 + NRHF(ISYMK)*(I - 1)
C
                     CALL DCOPY(NRHF(ISYMK),WORK(KOFF7),1,SCR3(KOFF8),1)
C
  331             CONTINUE
C
               ENDIF
C
  320       CONTINUE
C
  310    CONTINUE
C
  300 CONTINUE
C
      IF ((ICON .EQ. 4) .OR. (ICON .EQ. 5) .OR. (ICON .EQ. 6)) THEN
!
         CALL DSCAL(NT2BGD(ISYDIS),ZERO,SCR2,1)
C
         ISALIK = MULD2H(ISYDIS,ISYMHC)
C
         CALL DSCAL(NT2BGD(ISALIK),ZERO,SCR3,1)
!
      ELSE
!
         CALL DSCAL(NT2BGD(ISYDIS),-ONE,SCR2,1)
C
         ISALIK = MULD2H(ISYDIS,ISYMHC)
C
         CALL DSCAL(NT2BGD(ISALIK),-ONE,SCR3,1)
C
      ENDIF
!
!---------------------------------------------------------------
!     For ICON = 4 and ICON = 6 the real calculations begins here
!---------------------------------------------------------------
!
      DO 340 ISYMK = 1,NSYM
C
         ISYALG = MULD2H(ISYMK,ISYDIS)
         ISYALI = MULD2H(ISYMHC,ISYALG)
         NT1AOM = MAX(NT1AO(ISYALG),NT1AO(ISYALI))
C
         KSCR10 = 1
         KSCR11 = KSCR10 + N2BST(ISYALG)
         KEND1  = KSCR11 + NT1AOM
         LWRK1  = LWORK  - KEND1
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient space for allocation in CCRHS_D31')
         END IF
C
         DO 350 K = 1,NRHF(ISYMK)
C
            KOFF1 = IDSRHF(ISYALG,ISYMK) + NNBST(ISYALG)*(K - 1) + 1
            CALL CCSD_SYMSQ(DSRHF(KOFF1),ISYALG,WORK(KSCR10))
C
            ISYALI = ISYALG
            CALL DZERO(WORK(KSCR11),NT1AO(ISYALI))
C
C------------------------------
C           Usual contribution.
C------------------------------
C
            DO 360 ISYMI = 1,NSYM
C
               ISYMAL = MULD2H(ISYMI,ISYALI)
               ISYMG  = ISYMI
C
               NBASAL = MAX(NBAS(ISYMAL),1)
               NBASG = MAX(NBAS(ISYMG),1)
C
               KOFF2 = KSCR10 + IAODIS(ISYMAL,ISYMG)
               KOFF3 = ILMRHF(ISYMI) + 1
               KOFF4 = KSCR11 + IT1AO(ISYMAL,ISYMI)
C
               CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMI),NBAS(ISYMG),
     *                    ONE,WORK(KOFF2),NBASAL,XLAMDH(KOFF3),NBASG,
     *                    ZERO,WORK(KOFF4),NBASAL)
C
  360       CONTINUE
C
            NRHFK = MAX(NRHF(ISYMK),1)
            KOFF5 = IT2BGT(ISYMK,ISYALI) + K
C
            IF ((ICON .EQ. 4) .OR. (ICON .EQ. 5) .OR. (ICON .EQ. 6)
     *                                         ) THEN 
               CALL DAXPY(NT1AO(ISYALI),ONE,WORK(KSCR11),1,SCR2(KOFF5),
     *                    NRHFK)
!
            ELSE
               CALL DAXPY(NT1AO(ISYALI),TWO,WORK(KSCR11),1,SCR2(KOFF5),
     *                    NRHFK)
            ENDIF
C
C----------------------------------------------------
C           In 2C1 linear tronsformation extra  cont.
C----------------------------------------------------
C
            IF ((ICON .EQ. 3) .OR. (ICON .EQ. 1) 
     &           .OR. (ICON .EQ. 4) .OR. (ICON .EQ. 6)) THEN
C
               ISYALI = MULD2H(ISYALG,ISYMHC)
C
               CALL DZERO(WORK(KSCR11),NT1AO(ISYALI))
C
               DO 361 ISYMI = 1,NSYM
C
                  ISYMAL = MULD2H(ISYMI,ISYALI)
                  ISYMG  = MULD2H(ISYMI,ISYMHC)
C
                  NBASAL = MAX(NBAS(ISYMAL),1)
                  NBASG  = MAX(NBAS(ISYMG),1)
C
                  KOFF2 = KSCR10 + IAODIS(ISYMAL,ISYMG)
                  KOFF3 = IGLMRH(ISYMG,ISYMI) + 1
                  KOFF4 = KSCR11 + IT1AO(ISYMAL,ISYMI)
C
                  CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMI),
     &                       NBAS(ISYMG),ONE,WORK(KOFF2),NBASAL,
     &                       XLAMHC(KOFF3),NBASG,
     &                       ZERO,WORK(KOFF4),NBASAL)
C
  361          CONTINUE
C
               NRHFK = MAX(NRHF(ISYMK),1)
               KOFF5 = IT2BGT(ISYMK,ISYALI) + K
C
               IF (ICON .EQ. 4 ) THEN
               CALL DAXPY(NT1AO(ISYALI),ONE,WORK(KSCR11),1,
     &                    SCR3(KOFF5),NRHFK)
               ELSE IF (ICON .EQ. 6) THEN
               CALL DAXPY(NT1AO(ISYALI),XMONE,WORK(KSCR11),1,
     &                    SCR3(KOFF5),NRHFK)
               ELSE
               CALL DAXPY(NT1AO(ISYALI),TWO,WORK(KSCR11),1,
     &                    SCR3(KOFF5),NRHFK)
               ENDIF
C
            ENDIF
C
  350    CONTINUE
C
  340 CONTINUE
C
      IF (.NOT. DUMPCD) THEN
C
C-----------------------------------------
C     Back transformation to the AO basis.
C-----------------------------------------
C
      DO 400 ISYMAI = 1,NSYM
C
         ISYMK = MULD2H(ISYMAI,ISYDIS)
C
         NRHFK = MAX(NRHF(ISYMK),1)
C
         DO 410 ISYMI = 1,NSYM
C
            ISYMA = MULD2H(ISYMI,ISYMAI)
C
            NBASA = MAX(NBAS(ISYMA),1)
C
            DO 420 I = 1,NRHF(ISYMI)
C
               NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + 1
               MAI = IT1AO(ISYMA,ISYMI) + NBAS(ISYMA)*(I - 1) + 1
C
               KOFF1 = IT2BCT(ISYMK,ISYMAI) + NRHF(ISYMK)*(NAI - 1) + 1
               KOFF2 = ILMVIR(ISYMA) + 1
               KOFF3 = IT2BGT(ISYMK,ISYMAI) + NRHF(ISYMK)*(MAI - 1) + 1
C
               CALL DGEMM('N','T',NRHF(ISYMK),NBAS(ISYMA),NVIR(ISYMA),
     *                    HALF,SCR1(KOFF1),NRHFK,XLAMIP(KOFF2),NBASA,
     *                    ONE,SCR2(KOFF3),NRHFK)
C
  420       CONTINUE
C
  410    CONTINUE
C
  400 CONTINUE
C
C
      DO 500 ISYMK = 1,NSYM
C
         ISYMBJ = MULD2H(ISYMK,ISYDEL)
C
         DO 510 K = 1,NRHF(ISYMK)
C
            DO 520 ISYMJ = 1,NSYM
C
               ISYMB = MULD2H(ISYMJ,ISYMBJ)
C
               NBASB = MAX(NBAS(ISYMB),1)
               NVIRB = MAX(NVIR(ISYMB),1)
C
               KOFF1 = ILMVIR(ISYMB) + 1
               KOFF2 = IT2BCD(ISYMBJ,ISYMK) + NT1AM(ISYMBJ)*(K - 1)
     *               + IT1AM(ISYMB,ISYMJ) + 1
               KOFF3 = IT2BGD(ISYMBJ,ISYMK) + NT1AO(ISYMBJ)*(K - 1)
     *               + IT1AO(ISYMB,ISYMJ) + 1
C
               CALL DGEMM('N','N',NBAS(ISYMB),NRHF(ISYMJ),NVIR(ISYMB),
     *                    ONE,XLAMIP(KOFF1),NBASB,SCRM(KOFF2),NVIRB,
     *                    ZERO,SCR3(KOFF3),NBASB)
C
  520       CONTINUE
C
  510    CONTINUE
C
  500 CONTINUE
C
C---------------------------------------
C     Calculate the second contribution.
C---------------------------------------
C
      DO 600 ISYMAI = 1,NSYM
C
         ISYMK  = MULD2H(ISYMAI,ISYDIS)
         ISYMBJ = MULD2H(ISYMK,ISYDEL)
C
         IF (NRHF(ISYMK) .EQ. 0) GOTO 600
C
         IF (LWORK .LT. NT1AO(ISYMBJ)) THEN
            CALL QUIT('Insufficient work space in CCRHS_D31')
         ENDIF
C
         NTOTBJ = MAX(NT1AO(ISYMBJ),1)
         NRHFK  = MAX(NRHF(ISYMK),1)
C
         IF (.NOT. OMEGSQ) THEN
C
            DO 610 NAI = 1,NT1AO(ISYMAI)
C
               KOFF1 = IT2BGD(ISYMBJ,ISYMK) + 1
               KOFF2 = IT2BGT(ISYMK,ISYMAI)
     *               + NRHF(ISYMK)*(NAI - 1) + 1
C
               CALL DGEMV('N',NT1AO(ISYMBJ),NRHF(ISYMK),ONE,
     *                    SCR3(KOFF1),NTOTBJ,SCR2(KOFF2),1,
     *                    ZERO,WORK,1)
C
               IF (ISYMAI .EQ. ISYMBJ) THEN
                  WORK(NAI) = TWO*WORK(NAI)
               ENDIF
C
               DO 620 NBJ = 1,NT1AO(ISYMBJ)
                  NAIBJ = IT2AO(ISYMAI,ISYMBJ) + INDEX(NAI,NBJ)
                  OMEGA2(NAIBJ) = OMEGA2(NAIBJ) + HALF*WORK(NBJ)
  620          CONTINUE
C
  610       CONTINUE
C
         ELSE
C
            KOFF1 = IT2BGD(ISYMBJ,ISYMK)  + 1
            KOFF2 = IT2BGT(ISYMK,ISYMAI)  + 1
            KOFF3 = IT2AOS(ISYMBJ,ISYMAI) + 1
C
            CALL DGEMM('N','N',NT1AO(ISYMBJ),NT1AO(ISYMAI),NRHF(ISYMK),
     *                 HALF,SCR3(KOFF1),NTOTBJ,SCR2(KOFF2),NRHFK,
     *                 ONE,OMEGA2(KOFF3),NT1AO(ISYMBJ))
C
         ENDIF
C
  600 CONTINUE
C
      GOTO 999
C
C-------------------
C     I/O algorithm.
C-------------------
C
      ENDIF
C
C----------------------------------------------------------------------
C  Transform the alpha index of K(k,ai) to a.
C  for 2C1 transformation this means lamdpc is a C1 transformed lambda
C----------------------------------------------------------------------
C
      ISYAIK = MULD2H(ISYDIS,ISYMPC)
C
      DO 710 ISYMAI = 1,NSYM
C
         ISYMK = MULD2H(ISYMAI,ISYAIK)
         NRHFK = MAX(NRHF(ISYMK),1)
C
         DO 720 ISYMI = 1,NSYM
C
            ISYMA  = MULD2H(ISYMI,ISYMAI)
            ISYMAL = MULD2H(ISYMPC,ISYMA)
            ISYALI = MULD2H(ISYMAL,ISYMI)
            NBASAL = MAX(NBAS(ISYMAL),1)
C
            DO 730 I = 1,NRHF(ISYMI)
C
               NAI   = IT1AM(ISYMA,ISYMI)   + NVIR(ISYMA)*(I - 1) + 1
               MALI  = IT1AO(ISYMAL,ISYMI)  + NBAS(ISYMAL)*(I - 1) + 1
C
               KOFF1 = IT2BGT(ISYMK,ISYALI) + NRHF(ISYMK)*(MALI - 1) + 1
               KOFF2 = IGLMVI(ISYMAL,ISYMA) + 1
               KOFF3 = IT2BCT(ISYMK,ISYMAI) + NRHF(ISYMK)*(NAI - 1) + 1
C
               CALL DGEMM('N','N',NRHF(ISYMK),NVIR(ISYMA),NBAS(ISYMAL),
     *                    ONE,SCR2(KOFF1),NRHFK,XLAMPC(KOFF2),NBASAL,
     *                    FACTD ,SCR1(KOFF3),NRHFK)
!
  730       CONTINUE
  720    CONTINUE
  710 CONTINUE
!
!-----------------------------------------------
!     Transform the alpha index of K(k,ai) to a.
!     I is C1 transformed.
!-----------------------------------------------
!
      IF ((ICON .EQ. 3) .OR. (ICON .EQ. 1) 
     *    .OR. (ICON .EQ. 4) .OR. (ICON .EQ. 6)) THEN
C
         ISYAIK = MULD2H(ISYDIS,ISYMHC)
C
         DO 750 ISYMAI = 1,NSYM
C
            ISYMK = MULD2H(ISYMAI,ISYAIK)
            NRHFK = MAX(NRHF(ISYMK),1)
C
            DO 760 ISYMI = 1,NSYM
C
               ISYMA = MULD2H(ISYMI,ISYMAI)
               ISYMAL= ISYMA
               ISYALI= MULD2H(ISYMAL,ISYMI)
               NBASAL = MAX(NBAS(ISYMAL),1)
C
               DO 770 I = 1,NRHF(ISYMI)
C
                  NAI = IT1AM(ISYMA,ISYMI)
     *                + NVIR(ISYMA)*(I - 1) + 1
                  MALI = IT1AO(ISYMAL,ISYMI)
     *                 + NBAS(ISYMAL)*(I - 1) + 1
C
                  KOFF1 = IT2BGT(ISYMK,ISYALI)
     *                  + NRHF(ISYMK)*(MALI - 1) + 1
                  KOFF2 = ILMVIR(ISYMA) + 1
                  KOFF3 = IT2BCT(ISYMK,ISYMAI)
     *                  + NRHF(ISYMK)*(NAI - 1) + 1
C
                  CALL DGEMM('N','N',NRHF(ISYMK),NVIR(ISYMA),
     *                       NBAS(ISYMAL),ONE,SCR3(KOFF1),NRHFK,
     *                       XLAMDP(KOFF2),NBASAL,
     *                       ONE,SCR1(KOFF3),NRHFK)
C
  770          CONTINUE
  760       CONTINUE
  750    CONTINUE
C
      ENDIF
C
C---------------------------------------
C     Dump to disk the new contribution.
C---------------------------------------
C
C
      IF ((ICON .EQ. 2 ) .OR. (ICON .EQ. 5)) THEN
         IOFF = IT2DEL(IDEL) + 1
      ELSE
         IOFF = IT2DLR(IDEL,IV) + 1
      ENDIF
C
      IF (NT2BCD(ISYAIK) .GT. 0) THEN
         CALL PUTWA2(LUD,DFIL,SCR1,IOFF,NT2BCD(ISYAIK))
      ENDIF
C
  999 CONTINUE
C
      CALL QEXIT('CCRHS_D31')
C
      RETURN
      END
C  /* Deck ccrhs_e */
      SUBROUTINE CCRHS_E3(OMEGA2,OM2CONT,T2AM,EMAT1,EMAT2,WORK,LWORK,
     *                    ISYMTR,ISYMIM,OMEGA22,ANTISYM)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     Written by Henrik Koch & Ove Christiansen 20-Jan-1994
C     Symmetry 3-aug
C     Contraction of EI intermediates with double excitaion amplitudes.
C     It is assumed that the fock matrix is included. OC 13-1-1995
C
!     Generalized to the triplet case by Kasper Hald. march-1999
!     IF OM2CONT = .TRUE. => Symmetric permutation of (ai,bj)
!     is calculated in OMEGA2. IF ANTISYM = .TRUE. =>
!     antisymmetric prem. of (ai,bj) is calculated in OMEGA22.
!
C     Purpose: Calculate E-terms
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
      IMPLICIT NONE
!
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
!
      INTEGER LWORK, INDEX, ISYAIBJ, ISYMTR, ISYMIM, ISYMAI, ISYMCJ
      INTEGER ISYMBJ, NAI, ISYMJ, ISYMC, ISYMB, NVIRB, NVIRC
      INTEGER KOFF1, KOFF2, KOFF3, NBJ, NAIBJ, ISYMBK, ISYMK, NRHFK
!
#if defined (SYS_CRAY)
      REAL ZERO, ONE, TWO
      REAL EMAT1(*), EMAT2(*), T2AM(*), OMEGA2(*), OMEGA22(*)
      REAL WORK(LWORK)
#else
      DOUBLE PRECISION ZERO, ONE, TWO
      DOUBLE PRECISION EMAT1(*), EMAT2(*), T2AM(*), OMEGA2(*)
      DOUBLE PRECISION OMEGA22(*), WORK(LWORK)
#endif
!
      LOGICAL ANTISYM, OM2CONT
!
      PARAMETER(ZERO=0.0D00,ONE=1.0D00,TWO=2.0D00)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      CALL QENTER('CCRHS_E3')
C
C--------------------------------------------------------------
C     Contract and accumulate the first intermediate in OMEGA2.
C--------------------------------------------------------------
C
      ISYAIBJ = MULD2H(ISYMTR,ISYMIM)
C
      DO 300 ISYMAI = 1,NSYM
C
         ISYMCJ = MULD2H(ISYMAI,ISYMTR)
         ISYMBJ = MULD2H(ISYMAI,ISYAIBJ)
C
         IF (LWORK .LT. NT1AM(ISYMBJ)) THEN
            CALL QUIT('Insufficient space for allocation in CCRHS_E1')
         END IF
C
         DO 310 NAI = 1,NT1AM(ISYMAI)
C
            CALL DZERO(WORK,NT1AM(ISYMBJ))
C
            DO 320 ISYMJ = 1,NSYM
C
               ISYMC  = MULD2H(ISYMJ,ISYMCJ)
               ISYMB  = MULD2H(ISYMJ,ISYMBJ)
C
               NVIRB = MAX(NVIR(ISYMB),1)
               NVIRC = MAX(NVIR(ISYMC),1)
C
               KOFF1 = IMATAB(ISYMB,ISYMC) + 1
               KOFF2 = IT2SQ(ISYMCJ,ISYMAI) + NT1AM(ISYMCJ)*(NAI - 1)
     *                  + IT1AM(ISYMC,ISYMJ) + 1
               KOFF3 = IT1AM(ISYMB,ISYMJ) + 1
C
               CALL DGEMM('N','N',NVIR(ISYMB),NRHF(ISYMJ),
     *                    NVIR(ISYMC),ONE,EMAT1(KOFF1),NVIRB,
     *                    T2AM(KOFF2),NVIRC,
     *                    ONE,WORK(KOFF3),NVIRB)
  320          CONTINUE
C
            IF (OM2CONT) THEN
               IF (ISYMAI .EQ. ISYMBJ ) THEN
C
                 WORK(NAI) = TWO*WORK(NAI)
                 DO 330 NBJ = 1,NT1AM(ISYMBJ)
                    NAIBJ = IT2AM(ISYMAI,ISYMBJ) + INDEX(NAI,NBJ)
                    OMEGA2(NAIBJ) = OMEGA2(NAIBJ) + WORK(NBJ)
  330          CONTINUE
C
               ENDIF
C
               IF (ISYMAI .LT. ISYMBJ) THEN
C
                  DO 340 NBJ = 1,NT1AM(ISYMBJ)
                    NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                     + NT1AM(ISYMAI)*(NBJ - 1) + NAI
                     OMEGA2(NAIBJ) = OMEGA2(NAIBJ) + WORK(NBJ)
  340             CONTINUE
C
               ENDIF
C
               IF (ISYMBJ .LT. ISYMAI) THEN
C
                  DO 350 NBJ = 1,NT1AM(ISYMBJ)
                     NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                     + NT1AM(ISYMBJ)*(NAI - 1) + NBJ
                     OMEGA2(NAIBJ) = OMEGA2(NAIBJ) + WORK(NBJ)
  350             CONTINUE
C
               ENDIF
!
           ENDIF
           IF (ANTISYM) THEN
               IF (ISYMAI .EQ. ISYMBJ) THEN
!
                IF (NAI .GE. NT1AM(ISYMBJ)) THEN
!      
                  DO NBJ = 1, NT1AM(ISYMBJ)
                      NAIBJ = IT2AM(ISYMAI, ISYMBJ) + INDEX(NAI,NBJ)
                     OMEGA22(NAIBJ) = OMEGA22(NAIBJ) + WORK(NBJ)
                  ENDDO
                ELSE
!
                  DO NBJ = 1, NAI - 1
                     NAIBJ = IT2AM(ISYMAI, ISYMBJ) + INDEX(NAI,NBJ)
                     OMEGA22(NAIBJ) = OMEGA22(NAIBJ) + WORK(NBJ)
                  ENDDO 
                  DO NBJ = NAI + 1, NT1AM(ISYMBJ)
                     NAIBJ = IT2AM(ISYMAI, ISYMBJ) + INDEX(NAI,NBJ)
                     OMEGA22(NAIBJ) = OMEGA22(NAIBJ) - WORK(NBJ)
                  ENDDO
                ENDIF
               ENDIF
!
               IF (ISYMAI .LT. ISYMBJ) THEN
                  DO NBJ = 1, NT1AM(ISYMBJ)
                     NAIBJ = IT2AM(ISYMAI, ISYMBJ)
     *                     + NT1AM(ISYMAI)*(NBJ - 1) + NAI
                     OMEGA22(NAIBJ) = OMEGA22(NAIBJ) - WORK(NBJ)
                  ENDDO
               ENDIF
!
               IF (ISYMAI .GT. ISYMBJ) THEN
                  DO NBJ = 1, NT1AM(ISYMBJ)
                     NAIBJ = IT2AM(ISYMAI, ISYMBJ)
     *                     + NT1AM(ISYMBJ)*(NAI - 1) + NBJ
                     OMEGA22(NAIBJ) = OMEGA22(NAIBJ) + WORK(NBJ)
                  ENDDO
               ENDIF
           ENDIF
  310    CONTINUE
  300 CONTINUE
C
C-----------------------------------------------------
C     Contract and accumulate the second intermediate.
C-----------------------------------------------------
C
C
      DO 400 ISYMAI = 1,NSYM
C
         ISYMBK = MULD2H(ISYMAI,ISYMTR)
         ISYMBJ = MULD2H(ISYMAI,ISYAIBJ)
C
         IF (LWORK .LT. NT1AM(ISYMBJ)) THEN
            CALL QUIT('Insufficient space for allocation in CCRHS_E1')
         END IF
C
         DO 410 NAI = 1,NT1AM(ISYMAI)
C
            CALL DZERO(WORK,NT1AM(ISYMBJ))
C
            DO 420 ISYMB = 1,NSYM
C
               ISYMJ  = MULD2H(ISYMB,ISYMBJ)
               ISYMK  = MULD2H(ISYMJ,ISYMIM)
C
               NVIRB = MAX(NVIR(ISYMB),1)
               NRHFK = MAX(NRHF(ISYMK),1)
C
               KOFF1 = IT2SQ(ISYMBK,ISYMAI) + NT1AM(ISYMBK)*(NAI - 1)
     *               + IT1AM(ISYMB,ISYMK) + 1
               KOFF2 = IMATIJ(ISYMK,ISYMJ) + 1
               KOFF3 = IT1AM(ISYMB,ISYMJ) + 1
C
               CALL DGEMM('N','N',NVIR(ISYMB),NRHF(ISYMJ),
     *                    NRHF(ISYMK),ONE,T2AM(KOFF1),NVIRB,
     *                    EMAT2(KOFF2),NRHFK,
     *                    ONE,WORK(KOFF3),NVIRB)
  420       CONTINUE
C
C
            IF (OM2CONT) THEN
               IF (ISYMAI .EQ. ISYMBJ ) THEN
C
                  WORK(NAI) = TWO*WORK(NAI)
C
                  DO 430 NBJ = 1,NT1AM(ISYMBJ)
                     NAIBJ = IT2AM(ISYMAI,ISYMBJ) + INDEX(NAI,NBJ)
                     OMEGA2(NAIBJ) = OMEGA2(NAIBJ) - WORK(NBJ)
  430             CONTINUE
C
               ENDIF
C
               IF (ISYMAI .LT. ISYMBJ) THEN
C
                  DO 440 NBJ = 1,NT1AM(ISYMBJ)
                    NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                     + NT1AM(ISYMAI)*(NBJ - 1) + NAI
                     OMEGA2(NAIBJ) = OMEGA2(NAIBJ) - WORK(NBJ)
  440             CONTINUE
C
               ENDIF
C
               IF (ISYMBJ .LT. ISYMAI) THEN
C
                  DO 450 NBJ = 1,NT1AM(ISYMBJ)
                     NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                     + NT1AM(ISYMBJ)*(NAI - 1) + NBJ
                     OMEGA2(NAIBJ) = OMEGA2(NAIBJ) - WORK(NBJ)
  450             CONTINUE
C
               ENDIF
            ENDIF
            IF (ANTISYM) THEN
!
               IF (ISYMAI .EQ. ISYMBJ) THEN
!
                  DO NBJ = 1, NAI - 1
                     NAIBJ = IT2AM(ISYMAI, ISYMBJ) + INDEX(NAI,NBJ)
                     OMEGA22(NAIBJ) = OMEGA22(NAIBJ) - WORK(NBJ)
                  ENDDO
                  DO NBJ = NAI + 1, NT1AM(ISYMBJ)
                     NAIBJ = IT2AM(ISYMAI,ISYMBJ) + INDEX(NAI,NBJ)
                     OMEGA22(NAIBJ) = OMEGA22(NAIBJ) + WORK(NBJ)
                  ENDDO
                ENDIF
!
                IF (ISYMAI .LT. ISYMBJ) THEN
                   DO NBJ = 1, NT1AM(ISYMBJ)
                      NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                      + NT1AM(ISYMAI)*(NBJ - 1) + NAI
                      OMEGA22(NAIBJ) = OMEGA22(NAIBJ) + WORK(NBJ)
                   ENDDO
                ENDIF
!
                IF (ISYMAI .GT. ISYMBJ) THEN
                   DO NBJ = 1, NT1AM(ISYMBJ)
                      NAIBJ = IT2AM(ISYMAI, ISYMBJ)
     *                      + NT1AM(ISYMBJ)*(NAI - 1) + NBJ
                      OMEGA22(NAIBJ) = OMEGA22(NAIBJ) - WORK(NBJ)
                   ENDDO
                ENDIF
            ENDIF
  410    CONTINUE
  400 CONTINUE
C
      CALL QEXIT('CCRHS_E3')
C
      RETURN
      END
      SUBROUTINE CCRHS_C3(XINT,DSRHF,OMEGA2,T2AM,ISYMT2,
     *                    XLAMDP,XLAMIP,XLAMDH,
     *                    XLAMPC,ISYMPC,XLAMHC,ISYMHC,SCRM,WORK,LWORK,
     *                    IDEL,ISYMD,FACTC,ICON,LUC,CFIL,IV)
C
C     Written by Henrik Koch 3-Jan-1994
C     Symmetry by Henrik Koch and Alfredo Sanchez. 27-July-1994
C     Generalisation for CCLR by Ove Christiansen august-september 1995
C     (right transformation) and september 1996 (F-matrix).
!     Generalisation to the CCLR triplet case by Kasper Hald
!     11-march-1999.
C
C     Purpose: Calculate C-term.
C
      IMPLICIT NONE
#include <priunit.h>
#include <maxorb.h>
#include <ccorb.h>
#include <symsq.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <ccsdinp.h>
!
      INTEGER LWORK, ISYDIS, ISYMD, ISYAIK, ISYMPC, KSCR1, KSCR2
      INTEGER KSCR3, ICON, KEND1, LWRK1, ISYMT2, ISYMHC, IDEL, LUC
      INTEGER IV
!
#if defined (SYS_CRAY)
      REAL FACTC
      REAL XINT(*), DSRHF(*), OMEGA2(*), XLAMDH(*), WORK(LWORK)
      REAL XLAMDP(*), XLAMIP(*), SCRM(*), XLAMPC(*), XLAMHC(*)
      REAL T2AM(*)
#else
      DOUBLE PRECISION FACTC
      DOUBLE PRECISION XINT(*), DSRHF(*), OMEGA2(*), XLAMDH(*)
      DOUBLE PRECISION WORK(LWORK), XLAMDP(*), XLAMIP(*), SCRM(*)
      DOUBLE PRECISION XLAMPC(*), XLAMHC(*), T2AM(*)
#endif
C
      CHARACTER CFIL*(*)
C
      CALL QENTER('CCRHS_C3')
C
      ISYDIS = MULD2H(ISYMD,ISYMOP)
      ISYAIK = MULD2H(ISYDIS,ISYMPC)
C
C--------------------------------------
C     Dynamic allocation of work space.
C--------------------------------------
C
      KSCR1 = 1
      KSCR2 = KSCR1 + MAX(NT2BCD(ISYAIK),NT2BCD(ISYDIS))
      KSCR3 = KSCR2 + NT2BGD(ISYDIS)
      IF ((ICON .EQ. 2) .OR. (ICON .EQ. 5)) THEN
         KEND1  = KSCR3  + NT2BGD(ISYMD)
      ELSE
         KEND1  = KSCR3  + NT2BGD(ISYAIK)
      ENDIF
      LWRK1 = LWORK - KEND1
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space for allocation in CCRHS_C3')
      END IF
C
C--------------------------------------
C     Transpose the cluster amplitudes.
C--------------------------------------
C
      IF ((ICON .EQ. 2) .OR. (ICON .EQ. 3)) THEN
         IF (.NOT. T2TCOR) THEN
            CALL CCSD_T2TP(T2AM,WORK(KEND1),LWRK1,ISYMT2)
         ENDIF
         IF (.NOT. DUMPCD) CALL CCSD_T2MTP(SCRM,WORK(KEND1),LWRK1,ISYMD)
      ENDIF
C
C--------------------------------
C     Calculate the contribution.
C--------------------------------
C
         CALL CCRHS_C31(XINT,DSRHF,OMEGA2,T2AM,ISYMT2,SCRM,
     *                   WORK(KSCR1),WORK(KSCR2),WORK(KSCR3),XLAMDP,
     *                   XLAMIP,XLAMDH,XLAMPC,ISYMPC,XLAMHC,ISYMHC,
     *                   WORK(KEND1),LWRK1,
     *                   ISYDIS,IDEL,ISYMD,FACTC,ICON,LUC,CFIL,IV)
C
C--------------------------------------
C     Transpose the cluster amplitudes.
C--------------------------------------
C
      IF ((ICON .EQ. 2) .OR. (ICON .EQ. 3)) THEN 
         IF (.NOT. T2TCOR) THEN
            CALL CCSD_T2TP(T2AM,WORK(KEND1),LWRK1,ISYMT2)
         ENDIF
         IF (.NOT. DUMPCD) CALL CCSD_T2MTP(SCRM,WORK(KEND1),LWRK1,ISYMD)
      ENDIF
C
      CALL QEXIT('CCRHS_C3')
C
      RETURN
      END
      SUBROUTINE CCRHS_C31(XINT,DSRHF,OMEGA2,T2AM,ISYMT2,SCRM,SCR1,
     *                      SCR2,SCR3,XLAMDP,XLAMIP,XLAMDH,
     *                      XLAMPC,ISYMPC,XLAMHC,ISYMHC,WORK,
     *                      LWORK,ISYDIS,IDEL,ISYDEL,FACTC,ICON,LUC,
     *                      CFIL,IV)
C
C     Written by Henrik Koch 3-Jan-1994
C     Symmetry by Henrik Koch and Alfredo Sanchez. 27-July-1994
C
C     modification by Ove Christiansen 25-7-1995 to allow for a
C     general factor (FACTC) ( assumes DUMCD )
C     and - calculate intermediates for CCLR.
C
C     modification by Ove Christiansen 17-9-1996 for calculating
C     local C-intermediate for F-matrix transformation.
!
!     Modification by Kasper Hald 15-2-1999 for calculating the
!     local C-intermediate for triplet energy calculations.
!
C     Thus:
C
C     Modification to calculate terms in 2C1 right transformation in 
C     CCLR :
C
C                     1. if icon = 2 both contributions are calculated,
C                        otherwise if ICON =1:only the integral 
C                                       TILDE[(ki | ac)]
C                              = (k i-bar | a c) + (k i | a-bar c)
C
C                         3: (k i-bar | a c) + (k i | a-bar c)
C                              + FACTC*Sum(xT*int)
!                                where xT may be non total symmetric.
!
!                         4: (k i-bar | a c) - (k i | a-bar c)
!
C                      2. Allow for general transformation matrix for
C                            alpha to a(XLAMPC) and for i (XLAMHC).
C                            (the extra i transformation introduces new
C                             blocks which is only entered when
C                             icon =1 or 3)
C
C                      3. If icon diff. from 2 (we have linear response)
C                            The C intermediate is stored according to
C                            the number of simultaneous trial vector
C                            given by IV. This is ensured using IT2DLR.
C
!                      4. ICON = 4 is used for the triplet case.
!    
!
C     Thus in energy calc: icon = 2,fact = 1/2
C     For right transformation:
C         icon=1,fact=anything, iv = current vector being transformed
C     For F-matrix transformation:
C         icon=3,fact=1.0, NB - not implemented several vectors yet.
!     For Triplet calculation:
!         icon=4,fact=anything, iv = current vector being transformed
!
C     Purpose: Calculate C-term intermediate.
C
      IMPLICIT NONE
#include <priunit.h>
#include <maxorb.h>
#include <ccorb.h>
#include <symsq.h>
#include <ccsdsym.h>
#include <ccsdio.h>
!
      INTEGER LWORK, ISYMHC, INDEX, ISYAIK, ISYDIS, ISYMPC, ISAIK2
      INTEGER ISYMT2, ICON, ISYML, ISYMAG, KSCR10, KEND1, LWRK1
      INTEGER KOFF1, ISYMDL, ISYMD, ISYMK, ISYMA, ISYMG, NBASA
      INTEGER NBASG, NRHFK, KSCR11, KEND2, LWRK2, KOFF2, NDL, KOFF3
      INTEGER KOFF5, KOFF6, ISYMAI, NTOTDL, IOFF, IDEL, IERR
      INTEGER KOFF7, KOFF8, ISYMBG, ISYMI, ISYMB, NBASB, KSCR12
      INTEGER NAI, MAI, ISYMBJ, ISYDEL, ISYMJ, NVIRB, NTOTBJ
      INTEGER ISYMAJ, ISYMBI, NAJ, NBI, NBJ, NAIBJ, NAJBI, KOFF
      INTEGER NBIAJ, ISYMAL, ISYALI, NBASAL, MALI, IV, LUC
!
#if defined (SYS_CRAY)
      REAL ZERO, ONE, HALF, XMHALF, TWO, XMONE, FACTC
      REAL XINT(*), OMEGA2(*), T2AM(*), DSRHF(*), SCRM(*)
      REAL SCR1(*), SCR2(*), SCR3(*), XLAMDP(*), XLAMIP(*), XLAMDH(*)
      REAL XLAMPC(*), XLAMHC(ISYMHC), WORK(LWORK)
      REAL FACTOR1
#else
      DOUBLE PRECISION ZERO, ONE, HALF, XMHALF, TWO, XMONE, FACTC 
      DOUBLE PRECISION XINT(*), OMEGA2(*), T2AM(*), DSRHF(*), SCRM(*)
      DOUBLE PRECISION SCR1(*), SCR2(*), SCR3(*), XLAMDP(*), XLAMIP(*)
      DOUBLE PRECISION XLAMDH(*), XLAMPC(*), XLAMHC(ISYMHC)
      DOUBLE PRECISION WORK(LWORK), FACTOR1
#endif
!
      PARAMETER (ZERO=0.0D00,ONE=1.0D00,HALF=0.5D00,XMHALF=-0.5D00)
      PARAMETER (TWO=2.0D00,XMONE= -1.0D00)
      CHARACTER CFIL*(*)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      CALL QENTER('CCRHS_C31')
C
      ISYAIK = MULD2H(ISYDIS,ISYMPC)
      ISAIK2 = MULD2H(ISYDIS,ISYMT2)
      IF (ISYAIK .NE. ISAIK2) THEN
          CALL QUIT('Symmetry mismatch in CCRHS_C3')
      ENDIF
C
C-------------------------------------------------------
C     Calculate the integrals K(k,dl) = (k d | l delta).
C-------------------------------------------------------
C
      IF ((ICON .EQ. 2) .OR. (ICON .EQ. 3)) THEN
C
         DO 100 ISYML = 1,NSYM
C
            ISYMAG = MULD2H(ISYML,ISYDIS)
C
            DO 110 L = 1,NRHF(ISYML)
C
               KSCR10 = 1
               KEND1  = KSCR10 + N2BST(ISYMAG)
               LWRK1  = LWORK  - KEND1
               IF (LWRK1 .LT. 0) THEN
                  CALL QUIT('Not enough space for '//
     &                 'allocation in CCRHS_C31(1)')
               END IF
C
               KOFF1 = IDSRHF(ISYMAG,ISYML) + NNBST(ISYMAG)*(L-1) + 1
               CALL CCSD_SYMSQ(DSRHF(KOFF1),ISYMAG,WORK(KSCR10))
C
               DO 120 ISYMDL = 1,NSYM
C
                  ISYMD = MULD2H(ISYML,ISYMDL)
                  ISYMK = MULD2H(ISYMDL,ISYDIS)
                  ISYMA = ISYMK
                  ISYMG = ISYMD
C
                  NBASA = MAX(NBAS(ISYMA),1)
                  NBASG = MAX(NBAS(ISYMG),1)
                  NRHFK = MAX(NRHF(ISYMK),1)
C
                  KSCR11 = KEND1
                  KEND2  = KSCR11 + NRHF(ISYMK)*NBAS(ISYMG)
                  LWRK2  = LWORK  - KEND2
                  IF (LWRK2 .LT. 0) THEN
                     CALL QUIT('Not enough space for '//
     *                      'allocation in CCRHS_C31 (2)')
                  END IF
C
                  KOFF2 = ILMRHF(ISYMK) + 1
                  KOFF3 = IAODIS(ISYMA,ISYMG) + 1
C
                  CALL DGEMM('T','N',NRHF(ISYMK),NBAS(ISYMG),
     *                       NBAS(ISYMA),ONE,XLAMDP(KOFF2),NBASA,
     *                       WORK(KOFF3),NBASA,
     *                       ZERO,WORK(KSCR11),NRHFK)
C
                  NDL   = IT1AM(ISYMD,ISYML)
     *                  + NVIR(ISYMD)*(L - 1) + 1
                  KOFF5 = ILMVIR(ISYMD) + 1
                  KOFF6 = IT2BCT(ISYMK,ISYMDL)
     *                  + NRHF(ISYMK)*(NDL - 1) + 1
C
                  CALL DGEMM('N','N',NRHF(ISYMK),NVIR(ISYMD),
     *                       NBAS(ISYMG),ONE,WORK(KSCR11),NRHFK,
     *                       XLAMDH(KOFF5),NBASG,
     *                       ZERO,SCR1(KOFF6),NRHFK)
C
  120          CONTINUE
C
  110       CONTINUE
C
  100    CONTINUE
C
C-----------------------------------------
C        Calculate the first contribution.
C        Sum(dl)T(al,di)*I(lckd)
C-----------------------------------------
C
         IF (LWORK .LT. NT2BCD(ISYAIK)) THEN
            CALL QUIT('Insufficient work space in CCRHS_C31 (3)')
         ENDIF
C
         DO 200 ISYMK  = 1,NSYM
C
            ISYMAI = MULD2H(ISYAIK,ISYMK)
            ISYMDL = MULD2H(ISYDIS,ISYMK)
C
            NRHFK  = MAX(NRHF(ISYMK),1)
            NTOTDL = MAX(NT1AM(ISYMDL),1)
C
            KOFF1 = IT2BCT(ISYMK,ISYMDL) + 1
            KOFF2 = IT2SQ(ISYMDL,ISYMAI) + 1
            KOFF3 = IT2BCT(ISYMK,ISYMAI) + 1
C
            CALL DGEMM('N','N',NRHF(ISYMK),NT1AM(ISYMAI),NT1AM(ISYMDL),
     *                 ONE,SCR1(KOFF1),NRHFK,T2AM(KOFF2),NTOTDL,ZERO,
     *                 WORK(KOFF3),NRHFK)
C
  200    CONTINUE
!
            CALL DCOPY(NT2BCD(ISYAIK),WORK,1,SCR1,1)
!
      ENDIF
!
C----------------------------------------------------------
C     Calculate the integrals K(k,ai) = (k i | alfa delta).
C----------------------------------------------------------
C
      DO 300 ISYMA = 1,NSYM
C
         ISYMBG = MULD2H(ISYMA,ISYDIS)
C
         KSCR10 = 1
         KEND1  = KSCR10 + N2BST(ISYMBG)
         LWRK1  = LWORK  - KEND1
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT(
     *           'Not enough space for allocation in CCRHS_C31 (4)')
         END IF
C
         DO 310 A = 1,NBAS(ISYMA)
C
            KOFF1 = IDSAOG(ISYMA,ISYDIS) + NNBST(ISYMBG)*(A - 1) + 1
            CALL CCSD_SYMSQ(XINT(KOFF1),ISYMBG,WORK(KSCR10))
C
            DO 320 ISYMG = 1,NSYM
C
               ISYMI  = ISYMG
               ISYMB  = MULD2H(ISYMG,ISYMBG)
               ISYMK  = ISYMB
               ISYMAI = MULD2H(ISYMA,ISYMI)
C
               NBASB = MAX(NBAS(ISYMB),1)
               NBASG = MAX(NBAS(ISYMG),1)
               NRHFK = MAX(NRHF(ISYMK),1)
C
               KSCR11 = KEND1
               KSCR12 = KSCR11 + NRHF(ISYMK)*NBAS(ISYMG)
               KEND2  = KSCR12 + NRHF(ISYMK)*NRHF(ISYMI)
               LWRK2  = LWORK  - KEND2
               IF (LWRK2 .LT. 0) THEN
                  CALL QUIT('Not enough space for '//
     &                 'allocation in CCRHS_C31(5)')
               END IF
C
               KOFF2 = ILMRHF(ISYMK) + 1
               KOFF3 = KSCR10 + IAODIS(ISYMB,ISYMG)
C
               CALL DGEMM('T','N',NRHF(ISYMK),NBAS(ISYMG),NBAS(ISYMB),
     *                    ONE,XLAMDP(KOFF2),NBASB,WORK(KOFF3),NBASB,
     *                    ZERO,WORK(KSCR11),NRHFK)
C
               KOFF5 = ILMRHF(ISYMI) + 1
C
                   CALL DGEMM('N','N',NRHF(ISYMK),NRHF(ISYMI),
     &                        NBAS(ISYMG),ONE,WORK(KSCR11),NRHFK,
     &                        XLAMDH(KOFF5),NBASG,
     &                        ZERO,WORK(KSCR12),NRHFK)
C
               DO 330 I = 1,NRHF(ISYMI)
C
                  NAI = IT1AO(ISYMA,ISYMI) + NBAS(ISYMA)*(I-1) + A
C
                  KOFF8 = IT2BGT(ISYMK,ISYMAI)
     *                  + NRHF(ISYMK)*(NAI - 1) + 1
                  KOFF7 = KSCR12 + NRHF(ISYMK)*(I - 1)
C
                  CALL DCOPY(NRHF(ISYMK),WORK(KOFF7),1,SCR2(KOFF8),1)
C
  330          CONTINUE
C
C
C-------------------------------------------------------
C              In 2C1 linear transformation extra  cont.
C-------------------------------------------------------
C
               IF ((ICON .EQ. 3) .OR. (ICON .EQ. 1) .OR. 
     &              (ICON .EQ. 4)) THEN
C
                  ISYMI  = MULD2H(ISYMG,ISYMHC)
                  ISYMAI = MULD2H(ISYMA,ISYMI)
C
                  KEND2  = KSCR12 + NRHF(ISYMK)*NRHF(ISYMI)
                  LWRK2  = LWORK  - KEND2
                  IF (LWRK2 .LT. 0) THEN
                     CALL QUIT('Not enough space for '//
     &                    'allocation in CCRHS_C31')
                  END IF
C
                  KOFF5 = IGLMRH(ISYMG,ISYMI) + 1
C
                  IF (ICON .EQ. 4) THEN
                     FACTOR1 = XMONE
                  ELSE
                     FACTOR1 = ONE
                  ENDIF
!
                  CALL DGEMM('N','N',NRHF(ISYMK),NRHF(ISYMI),
     *                       NBAS(ISYMG),FACTOR1,WORK(KSCR11),NRHFK,
     *                       XLAMHC(KOFF5),NBASG,
     *                       ZERO,WORK(KSCR12),NRHFK)
C
                  DO 331 I = 1,NRHF(ISYMI)
C
                     NAI = IT1AO(ISYMA,ISYMI) + NBAS(ISYMA)*(I-1) + A
C
                     KOFF8 = IT2BGT(ISYMK,ISYMAI)
     *                     + NRHF(ISYMK)*(NAI - 1) + 1
                     KOFF7 = KSCR12 + NRHF(ISYMK)*(I - 1)
C
                     CALL DCOPY(NRHF(ISYMK),WORK(KOFF7),1,SCR3(KOFF8),1)
C
  331             CONTINUE
C
               ENDIF
C
  320       CONTINUE
C
  310    CONTINUE
C
  300 CONTINUE
                  
C
      IF (.NOT. DUMPCD) THEN
C
C-----------------------------------------
C     Back transformation to the AO basis.
C-----------------------------------------
C
      DO 400 ISYMAI = 1,NSYM
C
         ISYMK = MULD2H(ISYMAI,ISYDIS)
C
         NRHFK = MAX(NRHF(ISYMK),1)
C
         DO 410 ISYMI = 1,NSYM
C
            ISYMA = MULD2H(ISYMI,ISYMAI)
C
            NBASA = MAX(NBAS(ISYMA),1)
C
            DO 420 I = 1,NRHF(ISYMI)
C
               NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + 1
               MAI = IT1AO(ISYMA,ISYMI) + NBAS(ISYMA)*(I - 1) + 1
C
               KOFF1 = IT2BCT(ISYMK,ISYMAI) + NRHF(ISYMK)*(NAI - 1) + 1
               KOFF2 = ILMVIR(ISYMA) + 1
               KOFF3 = IT2BGT(ISYMK,ISYMAI) + NRHF(ISYMK)*(MAI - 1) + 1
C
               CALL DGEMM('N','T',NRHF(ISYMK),NBAS(ISYMA),NVIR(ISYMA),
     *                    XMHALF,SCR1(KOFF1),NRHFK,XLAMIP(KOFF2),NBASA,
     *                    ONE,SCR2(KOFF3),NRHFK)
C
  420       CONTINUE
C
  410    CONTINUE
C
  400 CONTINUE
C
C
      DO 500 ISYMK = 1,NSYM
C
         ISYMBJ = MULD2H(ISYMK,ISYDEL)
C
         DO 510 K = 1,NRHF(ISYMK)
C
            DO 520 ISYMJ = 1,NSYM
C
               ISYMB = MULD2H(ISYMJ,ISYMBJ)
C
               NBASB = MAX(NBAS(ISYMB),1)
               NVIRB = MAX(NVIR(ISYMB),1)
C
               KOFF1 = ILMVIR(ISYMB) + 1
               KOFF2 = IT2BCD(ISYMBJ,ISYMK) + NT1AM(ISYMBJ)*(K - 1)
     *               + IT1AM(ISYMB,ISYMJ) + 1
               KOFF3 = IT2BGD(ISYMBJ,ISYMK) + NT1AO(ISYMBJ)*(K - 1)
     *               + IT1AO(ISYMB,ISYMJ) + 1
C
               CALL DGEMM('N','N',NBAS(ISYMB),NRHF(ISYMJ),NVIR(ISYMB),
     *                    ONE,XLAMIP(KOFF1),NBASB,SCRM(KOFF2),NVIRB,
     *                    ZERO,SCR3(KOFF3),NBASB)
C
  520       CONTINUE
C
  510    CONTINUE
C
  500 CONTINUE
C
C---------------------------------------
C     Calculate the second contribution.
C
C     Alfredo will introduce the batching over ai before the
C     end of august 1994.
C---------------------------------------
C
      DO 600 ISYMAI = 1,NSYM
C
         ISYMK  = MULD2H(ISYMAI,ISYDIS)
         ISYMBJ = MULD2H(ISYMK,ISYDEL)
C
         IF (NRHF(ISYMK) .EQ. 0) GOTO 600
C
         IF (LWORK .LT. NT1AO(ISYMBJ)) THEN
            CALL QUIT('Insufficient work space in CCRHS_C1')
         ENDIF
C
         NTOTBJ = MAX(NT1AO(ISYMBJ),1)
C
         DO 610 ISYMI = 1,NSYM
C
            ISYMA = MULD2H(ISYMI,ISYMAI)
C
            DO 620 I = 1,NRHF(ISYMI)
C
               DO 630 A = 1,NBAS(ISYMA)
C
                  NAI = IT1AO(ISYMA,ISYMI) + NBAS(ISYMA)*(I-1) + A
C
                  KOFF1 = IT2BGD(ISYMBJ,ISYMK) + 1
                  KOFF2 = IT2BGT(ISYMK,ISYMAI)
     *                  + NRHF(ISYMK)*(NAI - 1) + 1
C
                  CALL DGEMV('N',NT1AO(ISYMBJ),NRHF(ISYMK),ONE,
     *                       SCR3(KOFF1),NTOTBJ,SCR2(KOFF2),1,
     *                       ZERO,WORK,1)
C
                  IF (.NOT. OMEGSQ) THEN
C
C
                  IF (ISYMAI .EQ. ISYMBJ) THEN
                     WORK(NAI) = TWO*WORK(NAI)
                  ENDIF
C
                  DO 640 ISYMJ = 1,NSYM
C
                     ISYMB  = MULD2H(ISYMJ,ISYMBJ)
                     ISYMAJ = MULD2H(ISYMA,ISYMJ)
                     ISYMBI = MULD2H(ISYMB,ISYMI)
C
                     DO 650 J = 1,NRHF(ISYMJ)
C
                        NAJ = IT1AO(ISYMA,ISYMJ) + NBAS(ISYMA)*(J-1) + A
C
                        DO 660 B = 1,NBAS(ISYMB)
C
                           NBI = IT1AO(ISYMB,ISYMI)
     *                         + NBAS(ISYMB)*(I-1) + B
                           NBJ = IT1AO(ISYMB,ISYMJ)
     *                         + NBAS(ISYMB)*(J-1) + B
C
                           NAIBJ = IT2AO(ISYMAI,ISYMBJ) + INDEX(NAI,NBJ)
                           NAJBI = IT2AO(ISYMAJ,ISYMBI) + INDEX(NAJ,NBI)
C
                           OMEGA2(NAIBJ) = OMEGA2(NAIBJ)-HALF*WORK(NBJ)
                           OMEGA2(NAJBI) = OMEGA2(NAJBI)-WORK(NBJ)
C
  660                   CONTINUE
  650                CONTINUE
  640             CONTINUE
C
C
                  ELSE
C
C
                  KOFF = IT2AOS(ISYMBJ,ISYMAI)
     *                 + NT1AO(ISYMBJ)*(NAI - 1) + 1
                  CALL DAXPY(NT1AO(ISYMBJ),-HALF,WORK,1,OMEGA2(KOFF),1)
C
                  DO 740 ISYMJ = 1,NSYM
C
                     ISYMB  = MULD2H(ISYMJ,ISYMBJ)
                     ISYMAJ = MULD2H(ISYMA,ISYMJ)
                     ISYMBI = MULD2H(ISYMB,ISYMI)
C
                     NBI = IT1AO(ISYMB,ISYMI) + NBAS(ISYMB)*(I-1) + 1
C
                     DO 750 J = 1,NRHF(ISYMJ)
C
                        NAJ = IT1AO(ISYMA,ISYMJ) + NBAS(ISYMA)*(J-1) + A
                        NBJ = IT1AO(ISYMB,ISYMJ) + NBAS(ISYMB)*(J-1) + 1
C
                        NBIAJ = IT2AOS(ISYMBI,ISYMAJ)
     *                        + NT1AO(ISYMBI)*(NAJ - 1) + NBI
C
                        CALL DAXPY(NBAS(ISYMB),-ONE,WORK(NBJ),1,
     *                             OMEGA2(NBIAJ),1)
C
  750                CONTINUE
  740             CONTINUE
C
                  ENDIF
C
  630          CONTINUE
  620       CONTINUE
C
  610    CONTINUE
  600 CONTINUE
C
      GOTO 999
C
C-------------------
C     I/O algorithm.
C-------------------
C
      ENDIF
C
C-----------------------------------------------
C     Transform the alpha index of K(k,ai) to a.
C-----------------------------------------------
C
      ISYAIK = MULD2H(ISYDIS,ISYMPC)
C
      IF ((ICON .EQ. 1) .OR. (ICON .EQ. 4)) THEN
          CALL DZERO(SCR1,NT2BCD(ISYAIK))
      ENDIF
C
      DO 810 ISYMAI = 1,NSYM
C
         ISYMK = MULD2H(ISYMAI,ISYAIK)
         NRHFK = MAX(NRHF(ISYMK),1)
C
         DO 820 ISYMI = 1,NSYM
C
            ISYMA = MULD2H(ISYMI,ISYMAI)
            ISYMAL= MULD2H(ISYMPC,ISYMA)
            ISYALI= MULD2H(ISYMAL,ISYMI)
            NBASAL = MAX(NBAS(ISYMAL),1)
C
            DO 830 I = 1,NRHF(ISYMI)
C
               NAI  = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + 1
               MALI = IT1AO(ISYMAL,ISYMI) + NBAS(ISYMAL)*(I - 1) + 1
C
               KOFF1 = IT2BGT(ISYMK,ISYALI) + NRHF(ISYMK)*(MALI- 1) + 1
               KOFF2 = IGLMVI(ISYMAL,ISYMA) + 1
               KOFF3 = IT2BCT(ISYMK,ISYMAI) + NRHF(ISYMK)*(NAI - 1) + 1
C
               CALL DGEMM('N','N',NRHF(ISYMK),NVIR(ISYMA),NBAS(ISYMAL),
     *                    ONE,SCR2(KOFF1),NRHFK,XLAMPC(KOFF2),NBASAL,
     *                    FACTC,SCR1(KOFF3),NRHFK)
C
  830       CONTINUE
  820    CONTINUE
  810 CONTINUE
C
C-----------------------------------------------
C     Transform the alpha index of K(k,ai) to a.
C     I is C1 transformed.
C-----------------------------------------------
C
      IF ((ICON .EQ. 3) .OR. (ICON .EQ. 1) .OR. (ICON .EQ. 4)) THEN
C
         ISYAIK = MULD2H(ISYDIS,ISYMHC)
C
         DO 850 ISYMAI = 1,NSYM
C
            ISYMK = MULD2H(ISYMAI,ISYAIK)
            NRHFK = MAX(NRHF(ISYMK),1)
C
            DO 860 ISYMI = 1,NSYM
C
               ISYMA = MULD2H(ISYMI,ISYMAI)
               ISYMAL= ISYMA
               ISYALI= MULD2H(ISYMAL,ISYMI)
               NBASAL = MAX(NBAS(ISYMAL),1)
C
               DO 870 I = 1,NRHF(ISYMI)
C
                  NAI = IT1AM(ISYMA,ISYMI)
     *                + NVIR(ISYMA)*(I - 1) + 1
                  MALI = IT1AO(ISYMAL,ISYMI)
     *                 + NBAS(ISYMAL)*(I - 1) + 1
C
                  KOFF1 = IT2BGT(ISYMK,ISYALI)
     *                  + NRHF(ISYMK)*(MALI - 1) + 1
                  KOFF2 = ILMVIR(ISYMA) + 1
                  KOFF3 = IT2BCT(ISYMK,ISYMAI)
     *                  + NRHF(ISYMK)*(NAI - 1) + 1
C
                  CALL DGEMM('N','N',NRHF(ISYMK),NVIR(ISYMA),
     *                       NBAS(ISYMAL),ONE,SCR3(KOFF1),NRHFK,
     *                       XLAMDP(KOFF2),NBASAL,
     *                       ONE,SCR1(KOFF3),NRHFK)
C
  870          CONTINUE
  860       CONTINUE
  850    CONTINUE
C
      ENDIF
C
C---------------------------------------------------------
C     Dump to disk the new contribution.
C---------------------------------------------------------
C
      IF ( ICON .EQ. 2 ) THEN
C
         IOFF = IT2DEL(IDEL) + 1
C
      ELSE
C
         IOFF = IT2DLR(IDEL,IV) + 1
C
      ENDIF
C
      IF (NT2BCD(ISYAIK) .GT. 0) THEN
         CALL PUTWA2(LUC,CFIL,SCR1,IOFF,NT2BCD(ISYAIK))
      ENDIF
C
  999 CONTINUE
C
      CALL QEXIT('CCRHS_C31')
C
      RETURN
      END
C  /* Deck ccrhs_t2tr */
      SUBROUTINE CCRHS3_T2TR(T2AM,WORK,LWORK,ISYM,IOPT)
C
C     Alfredo Sanchez and Henrik Koch 30-July 1994
C
!     19-03-99: Kasper Hald
!     Generalized to calculate only the last term i.e.
!     only exchange (IOPT = 2) 
!
!     IOPT = 1 : The normal 2T2 - T2
!     IOPT = 2 : Only exchange
!
C     PURPOSE:
C             Calculate two coulomb minus exchange of t2 amplitudes,
C             Calculate only minus the exchange term.
C             The amplitudes are assumed to be a square matrix.
C
      IMPLICIT NONE
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
!
      INTEGER LWORK, ISYMJ, ISYMB, ISYMBJ, ISYMAI, ISYM, NBJ, ISYMI
      INTEGER ISYMA, ISYMAJ, ISYMBI, KSCR1, KSCR2, KEND1, LWRK1
      INTEGER NRHFI, NBI, NAIBJ, NAJBI, IOPT, KOFF 
! 
#if defined (SYS_CRAY) 
      REAL ONE, TWO, ZERO, XMONE
      REAL T2AM(*), WORK(LWORK)
#else
      DOUBLE PRECISION ONE, TWO, ZERO, XMONE
      DOUBLE PRECISION T2AM(*), WORK(LWORK)
#endif
      PARAMETER (ONE = 1.0D00, TWO = 2.0D00, ZERO = 0.0D00)
      PARAMETER (XMONE = -1.0D00)
C
      CALL QENTER('CCRHS3_T2TR')
C
C----------------------------------------------------------
C     Calculate two coulomb minus exchange of t2-amplitude,
C     or minus exchange. 
C----------------------------------------------------------
C
      DO 100 ISYMJ = 1,NSYM
C
         DO 110 J = 1,NRHF(ISYMJ)
C
            DO 120 ISYMB = 1,NSYM
C
               ISYMBJ = MULD2H(ISYMB,ISYMJ)
               ISYMAI = MULD2H(ISYMBJ,ISYM)
C
               DO 130 B = 1,NVIR(ISYMB)
C
                  NBJ = IT1AM(ISYMB,ISYMJ)+NVIR(ISYMB)*(J-1)+B
C
                  DO 140 ISYMI = 1,ISYMJ
C
                     ISYMA  = MULD2H(ISYMI,ISYMAI)
                     ISYMAJ = MULD2H(ISYMA,ISYMJ)
                     ISYMBI = MULD2H(ISYMB,ISYMI)
C
                     KSCR1 = 1
                     KSCR2 = KSCR1 + NVIR(ISYMA)
                     KEND1 = KSCR2 + NVIR(ISYMA)
                     LWRK1 = LWORK - KEND1
                     IF (LWRK1 .LT. 0) THEN
                        CALL QUIT('Insufficient space in CCRHS3_T2TR')
                     ENDIF
C
                     IF (ISYMI .EQ. ISYMJ) THEN
                        NRHFI = J - 1
                     ELSE
                        NRHFI = NRHF(ISYMI)
                     END IF
C
                     DO 150 I = 1,NRHFI
C
                        NBI = IT1AM(ISYMB,ISYMI)+NVIR(ISYMB)*(I-1)+B
C
                        NAIBJ = IT2SQ(ISYMAI,ISYMBJ)
     *                        + NT1AM(ISYMAI)*(NBJ-1)
     *                        + IT1AM(ISYMA,ISYMI)+NVIR(ISYMA)*(I-1)+1
C
                        NAJBI = IT2SQ(ISYMAJ,ISYMBI)
     *                        + NT1AM(ISYMAJ)*(NBI-1)
     *                        + IT1AM(ISYMA,ISYMJ)+NVIR(ISYMA)*(J-1)+1
C
!
                           CALL DCOPY(NVIR(ISYMA),T2AM(NAIBJ),1,
     *                                WORK(KSCR1),1)
                           CALL DCOPY(NVIR(ISYMA),T2AM(NAJBI),1,
     *                                WORK(KSCR2),1)
C
                        IF (IOPT .EQ. 1) THEN
                           CALL DSCAL(NVIR(ISYMA),TWO,T2AM(NAIBJ),1)
                           CALL DSCAL(NVIR(ISYMA),TWO,T2AM(NAJBI),1)
                        ELSE IF (IOPT .EQ. 2) THEN
                           CALL DSCAL(NVIR(ISYMA),ZERO,T2AM(NAIBJ),1)
                           CALL DSCAL(NVIR(ISYMA),ZERO,T2AM(NAJBI),1)
                        ELSE
                           CALL QUIT('IOPT Mismatch in CCRHS3_T2TR ')
                        ENDIF
!
                        CALL DAXPY(NVIR(ISYMA),-ONE,WORK(KSCR2),1,
     *                             T2AM(NAIBJ),1)
                        CALL DAXPY(NVIR(ISYMA),-ONE,WORK(KSCR1),1,
     *                             T2AM(NAJBI),1)
C
  150               CONTINUE
C
  140             CONTINUE
C
  130          CONTINUE
C
  120       CONTINUE
C
  110    CONTINUE
C
  100 CONTINUE
C
      IF (IPRCC .GT. 20) THEN
         IF (IOPT .EQ. 1) THEN
            CALL AROUND('Two coulomb minus exchamge of t2am')
         ELSE IF (IOPT .EQ. 2) THEN
            CALL AROUND('The minus exchange of the T2AM')
         ENDIF
         DO 200 ISYMBJ = 1,NSYM
            ISYMAI = MULD2H(ISYMBJ,ISYM)
            KOFF = IT2SQ(ISYMAI,ISYMBJ) + 1
            WRITE(LUPRI,*)
            WRITE(LUPRI,*) 'Symmetry block:',ISYMBJ
            CALL OUTPUT(T2AM(KOFF),1,NT1AM(ISYMAI),1,NT1AM(ISYMBJ),
     *                  NT1AM(ISYMAI),NT1AM(ISYMBJ),1,LUPRI)
  200    CONTINUE
      END IF
C
      CALL QEXIT('CCRHS3_T2TR')
C
      RETURN
      END
C  /* Deck ccrhs_cio */
      SUBROUTINE CCRHS3_CIO(OMEGA2,T2AM,XLAMDH,WORK,LWORK,
     *                      ISYVEC,ISYCIM,LUC,CFIL,IV,IOPT,FACCN,
     *                      NORMALCONT,FACCP,PIJCONT,ANTISYM)
!
!     asm 17-aug-1994
!
!     Ove Christiansen 30-7-1995: modified to account for general
!                                 non. total symmetric vectors (ISYVEC)
!                                 and
!                                 intermediates (ISYCIM). LUC and CFIL
!                                 is used to control from which file
!                                 the intermediate is obtained.
!
!                                 if iopt = 1 the C intermediate is 
!                                 assumed to be as in energy clac.
!
!                                 if iopt ne. 1 we use the intermediate
!                                    on luc with address given according
!                                    to
!                                    transformed vector nr iv (iv is not
!                                    1 if several vectros are trans.
!                                    simultaneously.)
!
!
!     Kasper Hald 22-3-1999       Modified to calculate the triplet
!                                 contributions. Use ANTISYM and FACCN,
!                                 NORMALCONT, FACCP, PIJCONT
!
!                                 in energy calc: iv=1,iopt=1
!                                 FACCN = -HALF and FACCP = -1
!                                 NORMALCONT = .TRUE.
!                                 PIJCONT = .TRUE.
!
!     PURPOSE:
!             Calculate the C-term making I/O
!
!     N.B. The diagonal is set to zero since the diagonal elements
!          are identical zero in the triplet case.
!
      IMPLICIT NONE
!
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <maxorb.h>
#include <ccsdio.h>
!
      INTEGER LWORK, INDEX, ISAIBJ, ISYVEC, ISYCIM, ISYMAI, ISYMBJ
      INTEGER ISYMCK, ISYMDK, NT1AI, LENAI, LENMIN, NDISAI, NBATAI
      INTEGER ILSTAI, IBATAI, IFSTAI, KSCR1, KSCR2, KEND, LWRK1
      INTEGER KOFF1, ISYDEL, ISYMK, IDELTA, ID, IOPT, IOFF, IV
      INTEGER LEN, IERR, NAI, KAI, KOFF2, KOFF3, KOFF4, KOFF5, KOFF6
      INTEGER ISYMC, NBASD, NVIRC, NT1BJ, NT1CK, KOFF7, KOFF8
      INTEGER ISYMI, ISYMA, ISYMJ, ISYMB, ISYMAJ, ISYMBI, NAJ
      INTEGER CCRHSCOUNT, LUC
!
#if defined (SYS_CRAY)
      REAL ZERO, HALF, ONE, TWO, FACCN, FACCP
      REAL OMEGA2(*), T2AM(*), XLAMDH(*), WORK(LWORK)
#else
      DOUBLE PRECISION ZERO, HALF, ONE, TWO, FACCN, FACCP
      DOUBLE PRECISION OMEGA2(*), T2AM(*), XLAMDH(*), WORK(LWORK)
#endif
!
      PARAMETER (ZERO= 0.0D00, HALF= 0.5D00, ONE= 1.0D00, TWO= 2.0D00)
      CHARACTER CFIL*(*)
      LOGICAL ANTISYM, NORMALCONT, PIJCONT
!
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
!
      CALL QENTER('CCRHS3_CIO')
!
      IF (OMEGSQ) THEN
         WRITE(LUPRI,*) 'I/O in C-term not implemented for '//
     &        'square Omega2'
         CALL QUIT('OMEGSQ = .TRUE.  in CCRHS3_CIO')
      END IF
!
      ISAIBJ = MULD2H(ISYVEC,ISYCIM)
!
      DO 100 ISYMAI = 1,NSYM
!
         IF (NT1AM(ISYMAI) .EQ. 0) GOTO 100
!
         ISYMBJ = MULD2H(ISYMAI,ISAIBJ)
         ISYMCK = MULD2H(ISYVEC,ISYMBJ)
         ISYMDK = MULD2H(ISYCIM,ISYMAI)
!
!------------------------
!        Batch structure.
!------------------------
!
         NT1AI = NT1AM(ISYMAI)
!
         LENAI  = NT1AO(ISYMDK)
         LENMIN = 2*LENAI
         IF (LENMIN .EQ. 0) GOTO 100
!
         NDISAI = LWORK / LENMIN
         IF (NDISAI .LT. 1) THEN
            CALL QUIT('Insufficient space for '//
     &           'allocation in CCRHS3_CIO-1')
         END IF
         NDISAI = MIN(NDISAI,NT1AI)
!
         NBATAI = (NT1AI - 1) / NDISAI + 1
!
!--------------------------
!        Loop over batches.
!--------------------------
!
         ILSTAI = 0
         DO 110 IBATAI = 1,NBATAI
!
            IFSTAI = ILSTAI + 1
            ILSTAI = ILSTAI + NDISAI
            IF (ILSTAI .GT. NT1AI) THEN
               ILSTAI = NT1AI
               NDISAI = ILSTAI - IFSTAI + 1
            END IF
!
!-----------------------------
!           Memory allocation.
!-----------------------------
!
            KSCR1 = 1
            KSCR2 = KSCR1 + NDISAI*NT1AO(ISYMDK)
            KEND  = KSCR2 + NDISAI*NT1AO(ISYMDK)
            LWRK1 = LWORK - KEND
!
            IF (LWRK1 .LT. 0) THEN
               CALL QUIT('Insufficient space for '//
     &              'allocation in CCRHS3_CIO-2')
            END IF
!
!----------------------------------
!           Construct P(del k,#ai).
!----------------------------------
!
            KOFF1 = KSCR1
            DO 120 ISYDEL = 1,NSYM
!
               ISYMK = MULD2H(ISYDEL,ISYMDK)
!
               DO 130 IDELTA = 1,NBAS(ISYDEL)
!
                  ID = IDELTA + IBAS(ISYDEL)
!
                  IF (IOPT .EQ. 1 ) THEN
                     IOFF = IT2DEL(ID) + IT2BCT(ISYMK,ISYMAI)
     *                    + NRHF(ISYMK)*(IFSTAI - 1) + 1
                  ELSE
                     IOFF = IT2DLR(ID,IV) + IT2BCT(ISYMK,ISYMAI)
     *                    + NRHF(ISYMK)*(IFSTAI - 1) + 1
                  ENDIF
!
                  LEN  = NDISAI*NRHF(ISYMK)
!
                  IF (LEN .GT. 0) THEN
                     CALL GETWA2(LUC,CFIL,WORK(KOFF1),IOFF,LEN)
                  ENDIF
!
                  DO 140 NAI = IFSTAI,ILSTAI
!
                     KAI = NAI - IFSTAI + 1
!
                     KOFF2 = KOFF1 + NRHF(ISYMK)*(KAI - 1)
                     KOFF3 = KSCR2 + NT1AO(ISYMDK)*(KAI - 1)
     *                     + IT1AO(ISYDEL,ISYMK) + IDELTA - 1
!
                     CALL DCOPY(NRHF(ISYMK),WORK(KOFF2),1,WORK(KOFF3),
     *                          NBAS(ISYDEL))
!
  140             CONTINUE
!
                  KOFF1 = KOFF1 + LEN
!
  130          CONTINUE
  120       CONTINUE
!
!-----------------------------------------
!              Transform delta index to c.
!-----------------------------------------
!
            DO 150 NAI = IFSTAI,ILSTAI
!
               KAI = NAI - IFSTAI + 1
!
               DO 160 ISYMC = 1,NSYM
!
                  ISYDEL = ISYMC
                  ISYMK  = MULD2H(ISYMC,ISYMCK)
!
                  NBASD = MAX(NBAS(ISYDEL),1)
                  NVIRC = MAX(NVIR(ISYMC),1)
!
                  KOFF4 = ILMVIR(ISYDEL) + 1
                  KOFF5 = KSCR2 + NT1AO(ISYMDK)*(KAI - 1)
     *                  + IT1AO(ISYDEL,ISYMK)
                  KOFF6 = KSCR1 + NT1AM(ISYMCK)*(KAI - 1)
     *                  + IT1AM(ISYMC,ISYMK)
!
                  CALL DGEMM('T','N',NVIR(ISYMC),NRHF(ISYMK),
     *                       NBAS(ISYDEL),ONE,XLAMDH(KOFF4),NBASD,
     *                       WORK(KOFF5),NBASD,ZERO,WORK(KOFF6),
     *                       NVIRC)
!
  160          CONTINUE
  150       CONTINUE
!
!--------------------------------------------
!           Contract P(ck,#ai) with T(bj,ck).
!--------------------------------------------
!
            NT1BJ = MAX(NT1AM(ISYMBJ),1)
            NT1CK = MAX(NT1AM(ISYMCK),1)
!
            KOFF7 = IT2SQ(ISYMBJ,ISYMCK) + 1
!
            CALL DGEMM('N','N',NT1AM(ISYMBJ),NDISAI,NT1AM(ISYMCK),
     *                 ONE,T2AM(KOFF7),NT1BJ,WORK(KSCR1),NT1CK,
     *                 ZERO,WORK(KSCR2),NT1BJ)
!
!--------------------------------------------------
!           Scale the diagonal with zero if antisym
!           since the diagonal is then identical zero
!           If not antisym scale the diagonal with
!           two.
!--------------------------------------------------
!
               IF (ISYMBJ .EQ. ISYMAI) THEN
!
                  DO 170 NAI = IFSTAI,ILSTAI
                    KOFF8 = KSCR2 + NT1AM(ISYMBJ)*(NAI-IFSTAI) + NAI - 1
                    IF (ANTISYM) THEN
                    WORK(KOFF8) = ZERO * WORK(KOFF8)
                    ELSE
                    WORK(KOFF8) = TWO * WORK(KOFF8)
                    ENDIF
  170             CONTINUE
!
               END IF
!
!-----------------------------------------------
!           Add the result to the packed omega2.
!-----------------------------------------------
!
            DO 180 ISYMI = 1,NSYM
!
               ISYMA = MULD2H(ISYMI,ISYMAI)
!
               DO 190 I = 1,NRHF(ISYMI)
!
                  DO 200 A = 1,NVIR(ISYMA)
!
                     NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A
                     IF ((NAI .LT. IFSTAI) .OR. (NAI .GT. ILSTAI))
     *                  GOTO 200
!
                     DO 210 ISYMJ = 1,NSYM
!
                        ISYMB  = MULD2H(ISYMJ,ISYMBJ)
                        ISYMAJ = MULD2H(ISYMA,ISYMJ)
                        ISYMBI = MULD2H(ISYMB,ISYMI)
!
                        DO 220 J = 1,NRHF(ISYMJ)
!
                           NAJ = IT1AM(ISYMA,ISYMJ)
     *                         + NVIR(ISYMA)*(J-1) + A
!
                           CALL CC_PUTC3(WORK(KSCR2),OMEGA2,ISYMAI,
     *                                  ISYMAJ,ISYMBI,ISYMBJ,ISYMB,
     *                                  ISYMI,ISYMJ,NAI,NAJ,I,J,
     *                                  IFSTAI,FACCN,NORMALCONT,
     *                                  FACCP,PIJCONT,ANTISYM)
!
  220                   CONTINUE
  210                CONTINUE
  200             CONTINUE
  190          CONTINUE
  180       CONTINUE
!
  110    CONTINUE
  100 CONTINUE
!
      CALL QEXIT('CCRHS3_CIO')
!
      RETURN
      END
!  /* Deck cc_putc */
      SUBROUTINE CC_PUTC3(SCR2,OMEGA2,ISYMAI,ISYMAJ,ISYMBI,ISYMBJ,
     *                   ISYMB,ISYMI,ISYMJ,NAI,NAJ,I,J,IFSTAI,FACCN,
     *                   NORMALCONT,FACCP,PIJCONT,ANTISYM)
!
!     Ove Christiansen 30-10-1995: Put in C contribution in omega vector
!                                  avoid troble on cray with optimizatio
!
!     Kasper Hald Spring 1999 : Generalized to triplet excitation 
!                               energies
!.
!
      IMPLICIT NONE
!
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <maxorb.h>
#include <ccsdio.h>
!
      INTEGER INDEX, ISYMAI, ISYMBJ, ISYMB, NBJ, ISYMJ, KOFF9, NAI
      INTEGER IFSTAI, NAIBJ, ISYMAJ, ISYMBI, NBI, NAJBI, ISYMI, NAJ
! 
#if defined (SYS_CRAY) 
      REAL ZERO, HALF, ONE, TWO, FACCN, FACCP
      REAL SCR2(*), OMEGA2(*)
#else
      DOUBLE PRECISION ZERO, HALF, ONE, TWO, FACCN, FACCP
      DOUBLE PRECISION SCR2(*), OMEGA2(*)
#endif
!
      LOGICAL NORMALCONT, PIJCONT, ANTISYM
      PARAMETER (ZERO= 0.0D00, HALF= 0.5D00, ONE= 1.0D00, TWO= 2.0D00)
!
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
!
      CALL QENTER('CC_PUTC3')
!
!------------------------------
!     Symmetric cont.
!------------------------------
!
      IF (.NOT. ANTISYM) THEN
!
        IF (NORMALCONT) THEN
!
          IF ( ISYMAI .EQ. ISYMBJ ) THEN
!
            DO B = 1,NVIR(ISYMB)
!
               NBJ = IT1AM(ISYMB,ISYMJ)
     *             + NVIR(ISYMB)*(J-1) + B
               KOFF9 = NT1AM(ISYMBJ)*(NAI - IFSTAI) + NBJ
               NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *               + INDEX(NAI,NBJ)
               OMEGA2(NAIBJ) = OMEGA2(NAIBJ)
     *                    + FACCN * SCR2(KOFF9)
C
            ENDDO
C
           ENDIF
C
           IF ( ISYMAI .LT. ISYMBJ ) THEN
C
            DO B = 1,NVIR(ISYMB)
C
               NBJ = IT1AM(ISYMB,ISYMJ)
     *             + NVIR(ISYMB)*(J-1) + B
               KOFF9 = NT1AM(ISYMBJ)*(NAI - IFSTAI) + NBJ
               NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *               + NT1AM(ISYMAI)*(NBJ - 1) + NAI
               OMEGA2(NAIBJ) = OMEGA2(NAIBJ)
     *                    + FACCN * SCR2(KOFF9)
C
            ENDDO
C
           ENDIF
C
           IF ( ISYMBJ .LT. ISYMAI ) THEN
C
            DO B = 1,NVIR(ISYMB)
C
               NBJ = IT1AM(ISYMB,ISYMJ)
     *             + NVIR(ISYMB)*(J-1) + B
               KOFF9 = NT1AM(ISYMBJ)*(NAI - IFSTAI) + NBJ
               NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *               + NT1AM(ISYMBJ)*(NAI - 1) + NBJ
               OMEGA2(NAIBJ) = OMEGA2(NAIBJ)
     *                    + FACCN * SCR2(KOFF9)
C
            ENDDO
C
           ENDIF
!
      ENDIF
C
      IF (PIJCONT) THEN
!
         IF (ISYMAJ .EQ. ISYMBI) THEN
C
            DO B = 1,NVIR(ISYMB)
C
               NBI = IT1AM(ISYMB,ISYMI)
     *             + NVIR(ISYMB)*(I-1) + B
               NBJ = IT1AM(ISYMB,ISYMJ)
     *             + NVIR(ISYMB)*(J-1) + B
               KOFF9 = NT1AM(ISYMBJ)*(NAI - IFSTAI) + NBJ
               NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *               + INDEX(NAJ,NBI)
               OMEGA2(NAJBI) = OMEGA2(NAJBI) + FACCP*SCR2(KOFF9)
C
            ENDDO
C
          ENDIF
C
        IF (ISYMAJ .LT. ISYMBI) THEN
C
           DO B = 1,NVIR(ISYMB)
C
              NBI = IT1AM(ISYMB,ISYMI)
     *            + NVIR(ISYMB)*(I-1) + B
              NBJ = IT1AM(ISYMB,ISYMJ)
     *            + NVIR(ISYMB)*(J-1) + B
              KOFF9 = NT1AM(ISYMBJ)*(NAI - IFSTAI) + NBJ
              NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *              + NT1AM(ISYMAJ)*(NBI - 1)
     *              + NAJ
              OMEGA2(NAJBI) = OMEGA2(NAJBI) + FACCP*SCR2(KOFF9)
C
           ENDDO
C
         ENDIF
C
        IF (ISYMBI .LT. ISYMAJ) THEN
C
           DO B = 1,NVIR(ISYMB)
C
              NBI = IT1AM(ISYMB,ISYMI)
     *            + NVIR(ISYMB)*(I-1) + B
              NBJ = IT1AM(ISYMB,ISYMJ)
     *            + NVIR(ISYMB)*(J-1) + B
              KOFF9 = NT1AM(ISYMBJ)*(NAI - IFSTAI) + NBJ
              NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *              + NT1AM(ISYMBI)*(NAJ - 1)
     *              + NBI
              OMEGA2(NAJBI) = OMEGA2(NAJBI) + FACCP*SCR2(KOFF9)
C
           ENDDO
C
         ENDIF
C
        ENDIF
!
!-------------------------
!     ANTISYM cont.
!-------------------------
!
      ELSE   
!
      IF (NORMALCONT) THEN
!
          IF ( ISYMAI .EQ. ISYMBJ ) THEN
!
            DO B = 1,NVIR(ISYMB)
!
               NBI = IT1AM(ISYMB,ISYMI)
     *             + NVIR(ISYMB)*(I-1) + B
               NBJ = IT1AM(ISYMB,ISYMJ)
     *             + NVIR(ISYMB)*(J-1) + B
               KOFF9 = NT1AM(ISYMBJ)*(NAI - IFSTAI) + NBJ
               NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *               + INDEX(NAI,NBJ)
               IF (NAI .GT. NBJ) THEN
                  OMEGA2(NAIBJ) = OMEGA2(NAIBJ)
     *                       + FACCN * SCR2(KOFF9)
               ELSE IF (NAI .LT. NBJ) THEN
                  OMEGA2(NAIBJ) = OMEGA2(NAIBJ)
     *                       - FACCN * SCR2(KOFF9)
               ENDIF
!
            ENDDO
!
         ENDIF
!
         IF ( ISYMAI .LT. ISYMBJ ) THEN
C
            DO B = 1,NVIR(ISYMB)
C
               NBJ = IT1AM(ISYMB,ISYMJ)
     *             + NVIR(ISYMB)*(J-1) + B
               KOFF9 = NT1AM(ISYMBJ)*(NAI - IFSTAI) + NBJ
               NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *               + NT1AM(ISYMAI)*(NBJ - 1) + NAI
               OMEGA2(NAIBJ) = OMEGA2(NAIBJ)
     *                    - FACCN * SCR2(KOFF9)
C
            ENDDO
C
           ENDIF
C
           IF ( ISYMAI .GT. ISYMBJ ) THEN
C
            DO B = 1,NVIR(ISYMB)
C
               NBJ = IT1AM(ISYMB,ISYMJ)
     *             + NVIR(ISYMB)*(J-1) + B
               KOFF9 = NT1AM(ISYMBJ)*(NAI - IFSTAI) + NBJ
               NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *               + NT1AM(ISYMBJ)*(NAI - 1) + NBJ
               OMEGA2(NAIBJ) = OMEGA2(NAIBJ)
     *                    + FACCN * SCR2(KOFF9)
C
            ENDDO
C
           ENDIF
!
        ENDIF
!
        IF (PIJCONT) THEN
!
          IF (ISYMAJ .EQ. ISYMBI) THEN
C
            DO B = 1,NVIR(ISYMB)
C
               NBI = IT1AM(ISYMB,ISYMI)
     *             + NVIR(ISYMB)*(I-1) + B
               NBJ = IT1AM(ISYMB,ISYMJ)
     *             + NVIR(ISYMB)*(J-1) + B
               KOFF9 = NT1AM(ISYMBJ)*(NAI - IFSTAI) + NBJ
               NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *               + INDEX(NAJ,NBI)
!
               IF (NAJ .GT. NBI) THEN
                  OMEGA2(NAJBI) = OMEGA2(NAJBI) + FACCP*SCR2(KOFF9)
               ELSE IF (NAJ .LT. NBI) THEN
                  OMEGA2(NAJBI) = OMEGA2(NAJBI) - FACCP*SCR2(KOFF9)
               ENDIF
!
            ENDDO
!
          ENDIF
C
          IF (ISYMAJ .GT. ISYMBI) THEN
C
             DO B = 1,NVIR(ISYMB)
C
                NBI = IT1AM(ISYMB,ISYMI)
     *              + NVIR(ISYMB)*(I-1) + B
                NBJ = IT1AM(ISYMB,ISYMJ)
     *              + NVIR(ISYMB)*(J-1) + B
                KOFF9 = NT1AM(ISYMBJ)*(NAI - IFSTAI) + NBJ
                NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *                + NT1AM(ISYMBI)*(NAJ - 1)
     *                + NBI
!
                OMEGA2(NAJBI) = OMEGA2(NAJBI) + FACCP*SCR2(KOFF9)
C
             ENDDO
C
          ENDIF
C
          IF (ISYMAJ .LT. ISYMBI) THEN
C
             DO B = 1,NVIR(ISYMB)
C
                NBI = IT1AM(ISYMB,ISYMI)
     *              + NVIR(ISYMB)*(I-1) + B
                NBJ = IT1AM(ISYMB,ISYMJ)
     *              + NVIR(ISYMB)*(J-1) + B
                KOFF9 = NT1AM(ISYMBJ)*(NAI - IFSTAI) + NBJ
                NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *                + NT1AM(ISYMAJ)*(NBI - 1)
     *                + NAJ
!
               OMEGA2(NAJBI) = OMEGA2(NAJBI) - FACCP*SCR2(KOFF9)
C
             ENDDO
C
          ENDIF
C
        ENDIF
!
      ENDIF
!
      CALL QEXIT('CC_PUTC3')
!
      RETURN
      END
C  /* Deck ccrhs_dio */
      SUBROUTINE CCRHS3_DIO(OMEGA2,OM2CONT,T2AM,XLAMDH,WORK,LWORK,
     *                      ISYVEC,ISYDIM,LUD,DFIL,IV,IOPT,FACD,
     *                      ANTISYM,OMEGA22,FACD2)
C
C     asm 20-aug-1994
C
C     Ove Christiansen 30-7-1995: Modified to account for general
C                                 non. total symmetric vectors (ISYVEC)
C                                 and intermediates (ISYDIM). 
C                                 LUD and DFIL is
C                                 used to control from which file the
C                                 intermediate is obtained.
C
C                                 if iopt = 1 the D intermediate 
C                                    is assumed
C                                    to be as in energy calc.
C
C                                 if iopt ne. 1 we use the intermediate
C                                    on luc with address given 
C                                    according to transformed 
C                                    vector nr iv (iv is not 1
C                                    if several vectors are transformed
C                                    simultaneously.)
C
C                                 in energy calc: iv=1,iopt=1
C
!     Kasper Hald 22-3-1999: Generalized to the triplet case and for a 
!                            general factor FACD. Have also introduced
!                            the inputvar. OMEGA22, OM2CONT and
!                            FACD2 so can calculate both
!                            the symmetric and antisymmetric of a
!                            given D-term.
C     PURPOSE:
C             Calculate the D-term making I/O
C
      IMPLICIT NONE
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <maxorb.h>
#include <ccsdio.h>
!
      INTEGER LWORK, INDEX, ISYVEC, ISAIBJ, ISYMAI, ISYMBJ, ISYMCK
      INTEGER ISYMDK, NT1AI, LENAI, LENMIN, NDISAI, NBATAI
      INTEGER ILSTAI, IBATAI, IFSTAI, KSCR1, KSCR2, KEND, LWRK1, KOFF1
      INTEGER ISYDEL, ISYMK, IDELTA, ID, IOPT, IOFF, IV, LEN, IERR
      INTEGER NAI, KAI, KOFF2, KOFF3, ISYDIM, ISYMC, NBASD, NVIRC
      INTEGER KOFF4, KOFF5, KOFF6, NT1BJ, NT1CK, KOFF7, KOFF8, LUD
!
#if defined (SYS_CRAY)
      REAL ZERO, HALF, ONE, TWO, FACD, FACD2
      REAL OMEGA2(*), T2AM(*), XLAMDH(*), WORK(LWORK), OMEGA22(*)
#else
      DOUBLE PRECISION ZERO, HALF, ONE, TWO, FACD, FACD2
      DOUBLE PRECISION OMEGA2(*), T2AM(*), XLAMDH(*), WORK(LWORK)
      DOUBLE PRECISION OMEGA22(*)
#endif
!
      PARAMETER (ZERO= 0.0D00, HALF= 0.5D00, ONE= 1.0D00, TWO= 2.0D00)
      CHARACTER DFIL*(*)
!
      LOGICAL ANTISYM, OM2CONT, LOGIC1
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      CALL QENTER('CCRHS3_DIO')
C
      IF (OMEGSQ) THEN
         WRITE(LUPRI,*) 'I/O in D-term not implemented '//
     &        'for square Omega2'
         CALL QUIT('OMEGSQ = .TRUE.  in CCRHS3_DIO')
      END IF
C
      ISAIBJ = MULD2H(ISYVEC,ISYDIM)
C
      DO 100 ISYMAI = 1,NSYM
C
         IF (NT1AM(ISYMAI) .EQ. 0) GOTO 100
C
C
         ISYMBJ = MULD2H(ISYMAI,ISAIBJ)
         ISYMCK = MULD2H(ISYVEC,ISYMBJ)
         ISYMDK = MULD2H(ISYDIM,ISYMAI)
C
C------------------------
C        Batch structure.
C------------------------
C
         NT1AI = NT1AM(ISYMAI)
C
         LENAI  = NT1AO(ISYMDK)
         LENMIN = 2*LENAI
         IF (LENMIN .EQ. 0) GOTO 100
C
         NDISAI = LWORK / LENMIN
         IF (NDISAI .LT. 1) THEN
            CALL QUIT('Insufficient space for allocation in CCRHS3_DIO')
         END IF
         NDISAI = MIN(NDISAI,NT1AI)
C
         NBATAI = (NT1AI - 1) / NDISAI + 1
C
C--------------------------
C        Loop over batches.
C--------------------------
C
         ILSTAI = 0
         DO 110 IBATAI = 1,NBATAI
C
            IFSTAI = ILSTAI + 1
            ILSTAI = ILSTAI + NDISAI
            IF (ILSTAI .GT. NT1AI) THEN
               ILSTAI = NT1AI
               NDISAI = ILSTAI - IFSTAI + 1
            END IF
C
C-----------------------------
C           Memory allocation.
C-----------------------------
C
            KSCR1 = 1
            KSCR2 = KSCR1 + NDISAI*NT1AO(ISYMDK)
            KEND  = KSCR2 + NDISAI*NT1AO(ISYMDK)
            LWRK1 = LWORK - KEND
C
            IF (LWRK1 .LT. 0) THEN
               CALL QUIT('Insufficient space for '//
     &              'allocation in CCRHS_DIO')
            END IF
C
C----------------------------------
C           Construct P(del k,#ai).
C----------------------------------
C
            KOFF1 = KSCR1
            DO 120 ISYDEL = 1,NSYM
C
               ISYMK = MULD2H(ISYDEL,ISYMDK)
C
               DO 130 IDELTA = 1,NBAS(ISYDEL)
C
                  ID = IDELTA + IBAS(ISYDEL)
C
!------------------------------------
!           IOPT part.
!------------------------------------
!
                  IF (IOPT .EQ. 1 ) THEN
                     IOFF = IT2DEL(ID) + IT2BCT(ISYMK,ISYMAI)
     *                    + NRHF(ISYMK)*(IFSTAI - 1) + 1
                  ELSE
                     IOFF = IT2DLR(ID,IV) + IT2BCT(ISYMK,ISYMAI)
     *                    + NRHF(ISYMK)*(IFSTAI - 1) + 1
                  ENDIF
C
                  LEN  = NDISAI*NRHF(ISYMK)
C
                  IF (LEN .GT. 0) THEN
                     CALL GETWA2(LUD,DFIL,WORK(KOFF1),IOFF,LEN)
                  ENDIF
C
                  DO 140 NAI = IFSTAI,ILSTAI
C
                     KAI = NAI - IFSTAI + 1
C
                     KOFF2 = KOFF1 + NRHF(ISYMK)*(KAI - 1)
                     KOFF3 = KSCR2 + NT1AO(ISYMDK)*(KAI - 1)
     *                     + IT1AO(ISYDEL,ISYMK) + IDELTA - 1
C
                     CALL DCOPY(NRHF(ISYMK),WORK(KOFF2),1,WORK(KOFF3),
     *                          NBAS(ISYDEL))
C
  140             CONTINUE
C
                  KOFF1 = KOFF1 + LEN
C
  130          CONTINUE
  120       CONTINUE
C
C--------------------------------------
C           Transform delta index to c.
C--------------------------------------
C
            DO 150 NAI = IFSTAI,ILSTAI
C
               KAI = NAI - IFSTAI + 1
C
               DO 160 ISYMC = 1,NSYM
C
                  ISYDEL = ISYMC
                  ISYMK  = MULD2H(ISYMC,ISYMCK)
C
                  NBASD = MAX(NBAS(ISYDEL),1)
                  NVIRC = MAX(NVIR(ISYMC),1)
C
                  KOFF4 = ILMVIR(ISYDEL) + 1
                  KOFF5 = KSCR2 + NT1AO(ISYMDK)*(KAI - 1)
     *                  + IT1AO(ISYDEL,ISYMK)
                  KOFF6 = KSCR1 + NT1AM(ISYMCK)*(KAI - 1)
     *                  + IT1AM(ISYMC,ISYMK)
C
                  CALL DGEMM('T','N',NVIR(ISYMC),NRHF(ISYMK),
     *                       NBAS(ISYDEL),ONE,XLAMDH(KOFF4),NBASD,
     *                       WORK(KOFF5),NBASD,ZERO,WORK(KOFF6),
     *                       NVIRC)
C
  160          CONTINUE
  150       CONTINUE
C
C--------------------------------------------
C           Contract P(ck,#ai) with T(bj,ck).
C--------------------------------------------
C
            NT1BJ = MAX(NT1AM(ISYMBJ),1)
            NT1CK = MAX(NT1AM(ISYMCK),1)
C
            KOFF7 = IT2SQ(ISYMBJ,ISYMCK) + 1
C
            CALL DGEMM('N','N',NT1AM(ISYMBJ),NDISAI,NT1AM(ISYMCK),
     *                 ONE,T2AM(KOFF7),NT1BJ,WORK(KSCR1),NT1CK,
     *                 ZERO,WORK(KSCR2),NT1BJ)
C
C------------------------------
C           Scale the diagonal.
C------------------------------
C
            IF (OM2CONT) THEN
!
              IF (ISYMBJ .EQ. ISYMAI) THEN
C
                 DO NAI = IFSTAI,ILSTAI
                    KOFF8 = KSCR2 + NT1AM(ISYMBJ)*(NAI-IFSTAI) + NAI - 1
                    WORK(KOFF8) = TWO * WORK(KOFF8)
                 ENDDO
C
              END IF
C     
C-----------------------------------------------
C           Add the result to the packed omega2.
!           This term is SYMMETRIC in (ai,bj)
C-----------------------------------------------
!
              LOGIC1 = .FALSE.
!
              DO 180 NAI = IFSTAI,ILSTAI
!
                 CALL CC_PUTD3(WORK(KSCR2),OMEGA2,ISYMAI,ISYMBJ,NAI,
     *                         IFSTAI,FACD,LOGIC1)
  180         CONTINUE
!
            ENDIF
!
!------------------------------------
!           Scale the diagonal
!------------------------------------
!
            IF (ANTISYM) THEN
!
              IF (ISYMBJ .EQ. ISYMAI) THEN
C
                 DO 190 NAI = IFSTAI,ILSTAI
                    KOFF8 = KSCR2 + NT1AM(ISYMBJ)*(NAI-IFSTAI) + NAI - 1
                    WORK(KOFF8) = ZERO * WORK(KOFF8)
  190          CONTINUE
C
              END IF
C    
C--------------------------------------------------
C           Add the result to the packed omega2.
!           This term is ANTISYMMETRIC in (ai,bj)
C--------------------------------------------------
!
              LOGIC1 = .TRUE.
!
              DO NAI = IFSTAI,ILSTAI
!
                 CALL CC_PUTD3(WORK(KSCR2),OMEGA22,ISYMAI,ISYMBJ,NAI,
     *                         IFSTAI,FACD2,LOGIC1)
              ENDDO
!
            ENDIF
!
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('CCRHS3_DIO')
C
      RETURN
      END
C  /* Deck cc_putd3 */
      SUBROUTINE CC_PUTD3(SCR2,OMEGA2,ISYMAI,ISYMBJ,NAI,IFSTAI,FACD,
     *                    ANTISYM)
C
C     Ove Christiansen 30-10-1995: Put in D contribution in omega vector
C                                  avoid troble on cray 
C                                  with optimization.
C
!     Kasper Hald 22-3-1999: Generalized to the triplet case with
!                            ANTISYM and FACD.
      IMPLICIT NONE
!
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <maxorb.h>
#include <ccsdio.h>
!
      INTEGER INDEX, ISYMAI, ISYMBJ, NBJ, KOFF9, NAI, IFSTAI, NAIBJ
!
#if defined (SYS_CRAY)
      REAL ZERO, HALF, ONE, TWO, FACD
      REAL SCR2(*), OMEGA2(*)
#else
      DOUBLE PRECISION ZERO, HALF, ONE, TWO, FACD
      DOUBLE PRECISION SCR2(*), OMEGA2(*)
#endif
      PARAMETER (ZERO= 0.0D00, HALF= 0.5D00, ONE= 1.0D00, TWO= 2.0D00)
      LOGICAL ANTISYM
!
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
!
      CALL QENTER('CC_PUTD3')
!
      IF (.NOT. ANTISYM) THEN
!
         IF ( ISYMAI .EQ. ISYMBJ) THEN
            DO 190 NBJ = 1,NT1AM(ISYMBJ)
!
               KOFF9 = NT1AM(ISYMBJ)*(NAI-IFSTAI)+ NBJ
               NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *               + INDEX(NAI,NBJ)
!
               OMEGA2(NAIBJ) = OMEGA2(NAIBJ) + FACD*SCR2(KOFF9)
!
  190    CONTINUE
!
         ENDIF
!
         IF ( ISYMAI .LT. ISYMBJ) THEN
            DO 200 NBJ = 1,NT1AM(ISYMBJ)
!
               KOFF9 = NT1AM(ISYMBJ)*(NAI-IFSTAI)+ NBJ
               NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *               + NT1AM(ISYMAI)*(NBJ - 1) + NAI
               OMEGA2(NAIBJ) = OMEGA2(NAIBJ) + FACD*SCR2(KOFF9)
!
  200    CONTINUE
!
         ENDIF
!
         IF (ISYMBJ .LT. ISYMAI) THEN
            DO 210 NBJ = 1,NT1AM(ISYMBJ)
!
               KOFF9 = NT1AM(ISYMBJ)*(NAI-IFSTAI)+ NBJ
               NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *               + NT1AM(ISYMBJ)*(NAI - 1) + NBJ
               OMEGA2(NAIBJ) = OMEGA2(NAIBJ) + FACD*SCR2(KOFF9)
!
  210    CONTINUE
!
         ENDIF
!
!-----------------------
!     ANTISYM Cont.
!-----------------------
!
      ELSE
!
         IF ( ISYMAI .EQ. ISYMBJ) THEN
!
             DO NBJ  = 1,NT1AM(ISYMBJ)
               KOFF9 = NT1AM(ISYMBJ)*(NAI-IFSTAI)+ NBJ
               NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *               + INDEX(NAI,NBJ)
!
               IF (NAI .LT. NBJ) THEN
                  OMEGA2(NAIBJ) = OMEGA2(NAIBJ) - FACD*SCR2(KOFF9)
               ELSE IF (NAI .GT. NBJ) THEN
                  OMEGA2(NAIBJ) = OMEGA2(NAIBJ) + FACD*SCR2(KOFF9)
               ENDIF
!
             ENDDO
!
         ENDIF
!
      IF ( ISYMAI .LT. ISYMBJ) THEN
!
         DO NBJ = 1,NT1AM(ISYMBJ)
!
            KOFF9 = NT1AM(ISYMBJ)*(NAI-IFSTAI)+ NBJ
            NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *            + NT1AM(ISYMAI)*(NBJ - 1) + NAI
            OMEGA2(NAIBJ) = OMEGA2(NAIBJ) - FACD*SCR2(KOFF9)
!
         ENDDO
!
      ENDIF
!
      IF (ISYMBJ .LT. ISYMAI) THEN
         DO NBJ = 1,NT1AM(ISYMBJ)
!
            KOFF9 = NT1AM(ISYMBJ)*(NAI-IFSTAI)+ NBJ
            NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *            + NT1AM(ISYMBJ)*(NAI - 1) + NBJ
            OMEGA2(NAIBJ) = OMEGA2(NAIBJ) + FACD*SCR2(KOFF9)
!
         ENDDO
!
      ENDIF
!
      ENDIF
!
      CALL QEXIT('CC_PUTD3')
!
      RETURN
      END
C  /* Deck ccrhs_h3 */
      SUBROUTINE CCRHS_H3(DSRHF,OMEGA1,XLAMDP,XLAMDH,SCRM,
     *                   WORK,LWORK,ISYDIS,ISYDEL,ISYMTR,FACH)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     Written by Henrik Koch & Ove Christiansen 19-Jan-1994
C     Generalized to do linear transformation parts by
C     OC 30-1-1995
!     Generalized to a general factor FACH
!     Kasper Hald 25-3-99
!
C     Purpose: Calculate H-term.
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
      IMPLICIT NONE
!
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
!
      INTEGER LWORK, ISYDIS, ISYDEL, ISYMTR
!
#if defined (SYS_CRAY)
      REAL FACH
      REAL DSRHF(*), OMEGA1(*), XLAMDH(*), WORK(LWORK)
      REAL XLAMDP(*), SCRM(*)
#else
      DOUBLE PRECISION FACH
      DOUBLE PRECISION DSRHF(*), OMEGA1(*), XLAMDH(*), WORK(LWORK)
      DOUBLE PRECISION XLAMDP(*), SCRM(*)
#endif
C
      CALL QENTER('CCRHS_H3')
C
C--------------------------------
C     Calculate the contribution.
C--------------------------------
C
      CALL CCRHS_H31(DSRHF,OMEGA1,SCRM,XLAMDP,XLAMDH,WORK,LWORK,
     *              ISYDIS,ISYDEL,ISYMTR,FACH)
C
      CALL QEXIT('CCRHS_H3')
C
      RETURN
      END
      SUBROUTINE CCRHS_H31(DSRHF,OMEGA1,SCRM,XLAMDP,XLAMDH,WORK,LWORK,
     *                    ISYDIS,ISYDEL,ISYMTR,FACH)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     Written by Henrik Koch & Ove Christiansen 19-Jan-1994
C     Generalized to do linear transformation parts by
C     OC 30-1-1995
C
C     Purpose: Calculate H-term.
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
      IMPLICIT NONE
!
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
!
      INTEGER LWORK, INDEX, ISYAKL, ISYMTR, ISYDEL, ISYML, ISYMGB
      INTEGER ISYMAK, ISYMKI, KSCR1, KEND1, LWRK1, KOFF1, ISYMI
      INTEGER ISYMB, ISYMG, ISYMK, ISYMA, KSCR2, KSCR3, KEND2
      INTEGER LWRK2, NBASG, NBASB, NRHFK, NVIRA, KOFF2, KOFF3
      INTEGER KOFF4, KOFF5, KOFF6, ISYDIS
!
#if defined (SYS_CRAY)
      REAL ZERO, ONE, FACH
      REAL DSRHF(*),OMEGA1(*),SCRM(*)
      REAL XLAMDP(*),XLAMDH(*),WORK(LWORK)
#else
      DOUBLE PRECISION ZERO, ONE, FACH
      DOUBLE PRECISION DSRHF(*),OMEGA1(*),SCRM(*)
      DOUBLE PRECISION XLAMDP(*),XLAMDH(*),WORK(LWORK)
#endif
      PARAMETER(ZERO=0.0D00,ONE=1.0D00)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      CALL QENTER('CCRHS_H31')
C
C--------------------------------------
C     Calculate contribution.
C--------------------------------------
C
      ISYAKL = MULD2H(ISYMTR,ISYDEL)
C
      DO 100 ISYML = 1,NSYM
C
         ISYMGB = MULD2H(ISYML,ISYDIS)
         ISYMAK = MULD2H(ISYML,ISYAKL)
         ISYMKI = ISYMGB
C
         KSCR1 = 1
         KEND1 = KSCR1 + N2BST(ISYMGB)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Need : ',KEND1,'Available : ',LWORK
            CALL QUIT('Insufficient space in CCRHS_H1')
         ENDIF
         DO 110 L = 1,NRHF(ISYML)
C
            KOFF1 = IDSRHF(ISYMGB,ISYML) + NNBST(ISYMGB)*(L - 1) + 1
            CALL CCSD_SYMSQ(DSRHF(KOFF1),ISYMGB,WORK(KSCR1))
C
            DO 120 ISYMI = 1,NSYM
C
               ISYMB = ISYMI
               ISYMG = MULD2H(ISYMB,ISYMGB)
               ISYMK = ISYMG
               ISYMA = MULD2H(ISYMK,ISYMAK)
C
               KSCR2 = KEND1
               KSCR3 = KSCR2 + NBAS(ISYMG)*NRHF(ISYMI)
               KEND2 = KSCR3 + NRHF(ISYMK)*NRHF(ISYMI)
               LWRK2 = LWORK - KEND2
C
               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Need : ',KEND2,'Available : ',LWORK
                  CALL QUIT('Insufficient space in CCRHS_H1')
               ENDIF
C
               NBASG = MAX(NBAS(ISYMG),1)
               NBASB = MAX(NBAS(ISYMB),1)
               NRHFK = MAX(NRHF(ISYMK),1)
               NVIRA = MAX(NVIR(ISYMA),1)
C
               KOFF2 = KSCR1 + IAODIS(ISYMG,ISYMB)
               KOFF3 = ILMRHF(ISYMI) + 1
C
               CALL DGEMM('N','N',NBAS(ISYMG),NRHF(ISYMI),NBAS(ISYMB),
     *                    ONE,WORK(KOFF2),NBASG,XLAMDH(KOFF3),NBASB,
     *                    ZERO,WORK(KSCR2),NBASG)
C
               KOFF4 = ILMRHF(ISYMK) + 1
C
               CALL DGEMM('T','N',NRHF(ISYMK),NRHF(ISYMI),NBAS(ISYMG),
     *                    ONE,XLAMDP(KOFF4),NBASG,WORK(KSCR2),NBASG,
     *                    ZERO,WORK(KSCR3),NRHFK)
C
               KOFF5 = IT2BCD(ISYMAK,ISYML) + NT1AM(ISYMAK)*(L - 1)
     *               + IT1AM(ISYMA,ISYMK) + 1
               KOFF6 = IT1AM(ISYMA,ISYMI) + 1
C
               CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMK),
     *                    -FACH,SCRM(KOFF5),NVIRA,WORK(KSCR3),NRHFK,
     *                    ONE,OMEGA1(KOFF6),NVIRA)
C
  120       CONTINUE
C
  110    CONTINUE
C
  100 CONTINUE
C
      CALL QEXIT('CCRHS_H31')
C
      RETURN
      END
C  /* Deck ccrhs_g */
      SUBROUTINE CCRHS_G3(DSRHF,OMEGA1,XLAMP1,ISYMP1,XLAMH1,ISYMH1,SCRM,
     *                   WORK,LWORK,ISYDIS,ISYDEL,ISYMTR,FACG)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     Written by Henrik Koch & Ove Christiansen 19-Jan-1994
C     Generalized to calculated term of linear transformation
C     and handle different transformations on integral indices by 
C     OC 30-1-1995
C
C     G(a,i) = sum(cdk)[t(ci,dk)*Lackd]
C     G(a,i)for fixed del = sum(ck)[M(ci,k)*L(alfa gamma k]
C
C     XLAMP1 is the transformation matrix for a ; XLAMP or a oneindex 
C     transformed
C     XLAMH1 is the transformation matrix for c ; XLAMH or a oneindex 
C     transformed.
C     DSRHF is the (alfa gamma | k) array for a given delta.
C
C     not implemented yet with DSRHF and SCRM index transformed.
C
C     tested for energy with symmetry: ordinary XLAM matrices  
C     - OC 10-2-1995
C     tested for linear transformation without symmetry.  
C     - OC spring 95
C
C     Kasper Hald 25-3-1999 : Generalized to a general factor FACG
C
C     Purpose: Calculate G-term.
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
      IMPLICIT NONE
!
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
!
      INTEGER LWORK, ISYINT, ISYMH1, ISYALI, KSCR1, KEND1
      INTEGER LWRK1, ISYMTR, ISYMP1, ISYDIS, ISYDEL
!
#if defined (SYS_CRAY)
      REAL FACG
      REAL DSRHF(*), OMEGA1(*), XLAMP1(*)
      REAL WORK(LWORK), XLAMH1(*), SCRM(*)
#else
      DOUBLE PRECISION FACG
      DOUBLE PRECISION DSRHF(*), OMEGA1(*), XLAMP1(*)
      DOUBLE PRECISION WORK(LWORK), XLAMH1(*), SCRM(*)
#endif
C
      CALL QENTER('CCRHS_G3')
C
C------------------------
C     Dynamic allocation.
C------------------------
C
      ISYINT = MULD2H(ISYMH1,ISYMOP)
      ISYALI = MULD2H(ISYINT,ISYMTR)
C
      KSCR1  = 1
      KEND1  = KSCR1  + NT1AO(ISYALI)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Need : ',KEND1,'Available : ',LWORK
         CALL QUIT('Insufficient space in CCRHS_G')
      ENDIF
C
C--------------------------------
C     Calculate the contribution.
C--------------------------------
C
      CALL CCRHS_G31(DSRHF,OMEGA1,SCRM,XLAMP1,ISYMP1,XLAMH1,ISYMH1,
     *              WORK(KSCR1),WORK(KEND1),LWRK1,ISYDIS,ISYDEL,ISYMTR,
     *              FACG)
C
C
      CALL QEXIT('CCRHS_G3')
C
      RETURN
      END
      SUBROUTINE CCRHS_G31(DSRHF,OMEGA1,SCRM,XLAMP1,ISYMP1,XLAMH1,
     *             ISYMH1,SCR1,WORK,LWORK,ISYDIS,ISYDEL,ISYMTR,FACG)
C
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     Written by Henrik Koch & Ove Christiansen 19-Jan-1994
C     Generalized to calculated term of linear transformation
C     by OC 30-1-1995
!     FACG by KH 25-3-99
!
C     Purpose: Calculate G-term.
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
      IMPLICIT NONE
!
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
!
      INTEGER LWORK, INDEX, ISYINT, ISYMH1, ISYALI, ISYMTR, ISYMAI
      INTEGER ISYMP1, ISYDEL, ISYCIK, ISYMK, ISYDIS, ISYMAG
      INTEGER ISYMCI, ISYMGI, KSCR10, KEND1, LWRK1, KOFF1, ISYMI
      INTEGER ISYMG, ISYMA, ISYMC, NBASG, NBASA, NVIRC, KSCR11
      INTEGER KEND2, LWRK2, KOFF2, KOFF3, KOFF4, KOFF6, ISYMAL, NVIRA
!
#if defined (SYS_CRAY)
      REAL ZERO, ONE, TWO, FACG
      REAL DSRHF(*), OMEGA1(*), SCRM(*), SCR1(*)
      REAL XLAMP1(*), XLAMH1(*), WORK(LWORK)
#else
      DOUBLE PRECISION ZERO, ONE, TWO, FACG
      DOUBLE PRECISION DSRHF(*), OMEGA1(*), SCRM(*), SCR1(*)
      DOUBLE PRECISION XLAMP1(*), XLAMH1(*), WORK(LWORK)
#endif
!
      PARAMETER(ZERO=0.0D00,ONE=1.0D00)
      PARAMETER(TWO=2.0D00)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      CALL QENTER('CCRHS_G31')
C
      ISYINT = MULD2H(ISYMH1,ISYMOP)
      ISYALI = MULD2H(ISYINT,ISYMTR)
      ISYMAI = MULD2H(ISYALI,ISYMP1)
      ISYCIK = MULD2H(ISYMTR,ISYDEL)
C
      CALL DZERO(SCR1,NT1AO(ISYMAI))
C
      DO 100 ISYMK = 1,NSYM
C
         ISYMAG = MULD2H(ISYMK,ISYDIS)
         ISYMCI = MULD2H(ISYMK,ISYCIK)
         ISYMGI = MULD2H(ISYALI,ISYMAG)
C
         KSCR10 = 1
         KEND1  = KSCR10 + N2BST(ISYMAG)
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Need : ',KEND1,'Available : ',LWORK
            CALL QUIT('Insufficient space in CCRHS_G1')
         ENDIF
C
         DO 110 K = 1,NRHF(ISYMK)
C
            KOFF1 = IDSRHF(ISYMAG,ISYMK) + NNBST(ISYMAG)*(K - 1) + 1
            CALL CCSD_SYMSQ(DSRHF(KOFF1),ISYMAG,WORK(KSCR10))
C
            DO 120 ISYMI = 1,NSYM
C
               ISYMG = MULD2H(ISYMI,ISYMGI)
               ISYMA = MULD2H(ISYMG,ISYMAG)
               ISYMC = ISYMG
C
               NBASG = MAX(NBAS(ISYMG),1)
               NBASA = MAX(NBAS(ISYMA),1)
               NVIRC = MAX(NVIR(ISYMC),1)
C
               KSCR11 = KEND1
               KEND2  = KSCR11 + NBAS(ISYMG)*NRHF(ISYMI)
               LWRK2  = LWORK  - KEND2
C
               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Need : ',KEND2,'Available : ',LWORK
                  CALL QUIT('Insufficient space in CCRHS_G1')
               ENDIF
C
               KOFF2 = IGLMVI(ISYMG,ISYMC) + 1
               KOFF3 = IT2BCD(ISYMCI,ISYMK) + NT1AM(ISYMCI)*(K - 1)
     *               + IT1AM(ISYMC,ISYMI) + 1
C
               CALL DGEMM('N','N',NBAS(ISYMG),NRHF(ISYMI),NVIR(ISYMC),
     *                    ONE,XLAMH1(KOFF2),NBASG,SCRM(KOFF3),NVIRC,
     *                    ZERO,WORK(KSCR11),NBASG)
C
               KOFF4 = KSCR10 + IAODIS(ISYMA,ISYMG)
               KOFF6 = IT1AO(ISYMA,ISYMI) + 1
C
               CALL DGEMM('N','N',NBAS(ISYMA),NRHF(ISYMI),NBAS(ISYMG),
     *                    ONE,WORK(KOFF4),NBASA,WORK(KSCR11),NBASG,
     *                    ONE,SCR1(KOFF6),NBASA)
C
  120       CONTINUE
C
  110    CONTINUE
C
  100 CONTINUE
C
C----------------------------------------------
C     Accumulation into OMEGA1 in the MO basis.
C----------------------------------------------
C
      DO 200 ISYMI = 1,NSYM
C
         ISYMA = MULD2H(ISYMI,ISYMAI)
         ISYMAL= MULD2H(ISYMI,ISYALI)
C
         NBASA = MAX(NBAS(ISYMA),1)
         NVIRA = MAX(NVIR(ISYMA),1)
C
         KOFF1 = IGLMVI(ISYMAL,ISYMA) + 1
         KOFF2 = IT1AO(ISYMA,ISYMI) + 1
         KOFF3 = IT1AM(ISYMA,ISYMI) + 1
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMA),FACG,
     *              XLAMP1(KOFF1),NBASA,SCR1(KOFF2),NBASA,ONE,
     *              OMEGA1(KOFF3),NVIRA)
C
  200 CONTINUE
C
      CALL QEXIT('CCRHS_G31')
C
      RETURN
      END
C  /* Deck cc_rvec3 */
      SUBROUTINE CC_RVEC3(LU,FIL,LLEN,LEN,NR,IDISP,VEC)
!
!     Kasper Hald : April 1. 1999 (And that's not even a joke)
!
!     The routine reads LEN elements from the file FIL with
!     logical unit number LU. The address is given by the length
!     of each file multiplied with (NR - 1). In the triplet 
!     case we store different length files so a displacement IDISP
!     can be given.
!
      IMPLICIT NONE
#if defined (SYS_CRAY)
      REAL VEC(*)
#else
      DOUBLE PRECISION VEC(*)
#endif
      CHARACTER FIL*(*)
      INTEGER LU, LLEN, LEN, NR, IDISP, IOFF, IERR
!
      CALL QENTER('CC_RVEC3')
!
      IOFF = 1 + LLEN*(NR-1) + IDISP
!
      IF (LEN .GT. 0) THEN
         CALL GETWA2(LU,FIL,VEC,IOFF,LEN)
      ENDIF
!
      CALL QEXIT('CC_RVEC3')
!
      RETURN
      END
C  /* Deck cc_wvec3 */
      SUBROUTINE CC_WVEC3(LU,FIL,LLEN,LEN,NR,IDISP,VEC)
!
!     Kasper Hald April 1. 1999 (And that's not even a joke)
!
!     Writes the vector VEC to the file FIL with logical unit number 
!     LU. The address is calculated from LLEN, NR and the displacement
!     IDISP.
!
      IMPLICIT NONE
!
#if defined (SYS_CRAY)
      REAL VEC(*)
#else
      DOUBLE PRECISION VEC(*)
#endif
      CHARACTER FIL*(*)
      INTEGER LU, LLEN, LEN, NR, IDISP, IOFF, IERR
!
      CALL QENTER('CC_WVEC3')
!
      IOFF = 1 + LLEN*(NR-1) + IDISP
!
      IF (LEN .GT. 0) THEN
         CALL PUTWA2(LU,FIL,VEC,IOFF,LEN)
      ENDIF
!
      CALL QEXIT('CC_WVEC3')
!
      RETURN
      END
C  /* Deck ccrhs_t2bt */
      SUBROUTINE CCRHS3_T2BT(T2AM,WORK,LWORK,ISYM,IOPT)
C
C     Alfredo Sanchez and Henrik Koch 30-July 1994
C
C     Kasper Hald 17-may 1999
C
C     Backtransform -exchange
C
C     PURPOSE:
C             Back transform t2 amplitudes.
C             The amplitudes are assumed to be a square matrix.
C
C             IOPT = 1 : 2T2 - T2 -> T2
C             IOPT = 2 : 0   - T2  (Pure exchange) -> T2
C
      IMPLICIT NONE
C
      INTEGER LWORK, ISYMJ, ISYMB, ISYMBJ, ISYMAI, NBJ, ISYMI, ISYMA
      INTEGER ISYMBI, KSCR1, ISYM, ISYMAJ, KSCR2, KEND1, LWRK1
      INTEGER NRHFI, NBI, NAIBJ, NAJBI, IOPT, KOFF
C
#if defined (SYS_CRAY)
      REAL ONETHD, TWOTHD, ONEMINUS
      REAL T2AM(*), WORK(LWORK)
#else
      DOUBLE PRECISION ONETHD, TWOTHD, ONEMINUS
      DOUBLE PRECISION T2AM(*), WORK(LWORK)
#endif
      PARAMETER(ONETHD = 1.0D00/3.0D00,TWOTHD = 2.0D00/3.0D00)
      PARAMETER(ONEMINUS = -1.0D00)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('CCRHS3_T2BT')
C
C----------------------------------
C     Back transform t2-amplitudes.
C----------------------------------
C
      DO 100 ISYMJ = 1,NSYM
C
         DO 110 J = 1,NRHF(ISYMJ)
C
            DO 120 ISYMB = 1,NSYM
C
               ISYMBJ = MULD2H(ISYMB,ISYMJ)
               ISYMAI = MULD2H(ISYMBJ,ISYM)
C
               DO 130 B = 1,NVIR(ISYMB)
C
                  NBJ = IT1AM(ISYMB,ISYMJ)+NVIR(ISYMB)*(J-1)+B
C
                  DO 140 ISYMI = 1,ISYMJ
C
                     ISYMA  = MULD2H(ISYMI,ISYMAI)
                     ISYMAJ = MULD2H(ISYMA,ISYMJ)
                     ISYMBI = MULD2H(ISYMB,ISYMI)
C
                     KSCR1 = 1
                     IF (IOPT .EQ. 1) THEN
                        KSCR2 = KSCR1 + NVIR(ISYMA)
                        KEND1 = KSCR2 + NVIR(ISYMA)
                     ELSE IF (IOPT .EQ. 2) THEN
                        KEND1 = KSCR1 + NVIR(ISYMA)
                     ENDIF
                     LWRK1 = LWORK - KEND1
                     IF (LWRK1 .LT. 0) THEN
                        CALL QUIT('Insufficient space in CCRHS3_T2TR')
                     ENDIF
C
                     IF (ISYMI .EQ. ISYMJ) THEN
                        NRHFI = J - 1
                     ELSE
                        NRHFI = NRHF(ISYMI)
                     END IF
C
                     DO 150 I = 1,NRHFI
C
                        NBI = IT1AM(ISYMB,ISYMI)+NVIR(ISYMB)*(I-1)+B
C
                        NAIBJ = IT2SQ(ISYMAI,ISYMBJ)
     *                        + NT1AM(ISYMAI)*(NBJ-1)
     *                        + IT1AM(ISYMA,ISYMI)+NVIR(ISYMA)*(I-1)+1
C
                        NAJBI = IT2SQ(ISYMAJ,ISYMBI)
     *                        + NT1AM(ISYMAJ)*(NBI-1)
     *                        + IT1AM(ISYMA,ISYMJ)+NVIR(ISYMA)*(J-1)+1
C
                        IF (IOPT .EQ. 1) THEN
!
                        CALL DCOPY(NVIR(ISYMA),T2AM(NAIBJ),1,
     *                             WORK(KSCR1),1)
                        CALL DCOPY(NVIR(ISYMA),T2AM(NAJBI),1,
     *                             WORK(KSCR2),1)
C
                        CALL DSCAL(NVIR(ISYMA),TWOTHD,T2AM(NAIBJ),1)
                        CALL DSCAL(NVIR(ISYMA),TWOTHD,T2AM(NAJBI),1)
C
                        CALL DAXPY(NVIR(ISYMA),ONETHD,WORK(KSCR2),1,
     *                             T2AM(NAIBJ),1)
                        CALL DAXPY(NVIR(ISYMA),ONETHD,WORK(KSCR1),1,
     *                             T2AM(NAJBI),1)
C
                        ELSE IF (IOPT .EQ. 2) THEN
C
                        CALL DSCAL(NVIR(ISYMA),ONEMINUS,T2AM(NAIBJ),1)
                        CALL DSCAL(NVIR(ISYMA),ONEMINUS,T2AM(NAJBI),1)
C
                        CALL DCOPY(NVIR(ISYMA),T2AM(NAIBJ),1,
     *                             WORK(KSCR1),1)
                        CALL DCOPY(NVIR(ISYMA),T2AM(NAJBI),1,
     *                             T2AM(NAIBJ),1)
                        CALL DCOPY(NVIR(ISYMA),WORK(KSCR1),1,
     *                             T2AM(NAJBI),1)
C
                        ELSE
                           CALL QUIT('IOPT mismatch in CCRHS3_T2BT')
                        ENDIF
C
  150               CONTINUE
C
  140             CONTINUE
C
  130          CONTINUE
C
  120       CONTINUE
C
  110    CONTINUE
C
  100 CONTINUE
C
      IF (IPRCC .GT. 20) THEN
         CALL AROUND('Back-transformed t2am')
         DO 200 ISYMBJ = 1,NSYM
            ISYMAI = MULD2H(ISYMBJ,ISYM)
            KOFF = IT2SQ(ISYMAI,ISYMBJ) + 1
            WRITE(LUPRI,*)
            WRITE(LUPRI,*) 'Symmetry block:',ISYMBJ
            CALL OUTPUT(T2AM(KOFF),1,NT1AM(ISYMAI),1,NT1AM(ISYMBJ),
     *                  NT1AM(ISYMAI),NT1AM(ISYMBJ),1,LUPRI)
  200    CONTINUE
      END IF
C
      CALL QEXIT('CCRHS3_T2BT')
C
      RETURN
      END
C  /* Deck ccrhs3_cd */
      SUBROUTINE CCRHS3_CD(LUD,DFIL,LUC,CFIL,IDEL,WORK,LWORK,
     *                     LUCD,CDFIL,ISYMD,ISYMPC)
!
!     Written by Kasper Hald.
!
!
!     Purpose : Calculate (3)D - (1)C and write to disk.
!
!
      IMPLICIT NONE
!
#include <priunit.h>
#include <ccsdsym.h>
#include <maxorb.h>
#include <ccorb.h>
#include <ccsdio.h>
!
      INTEGER LWORK,LUC,LUCD,LUD,IDEL,ISYMTR,ISYMD, ISYMPC
      INTEGER IERRCD, IERRC, IERRD, KSCR1, KSCR2, KEND1, LWRK1
      INTEGER ISYAIK, ISYDIS, IOFF
!
#if defined (SYS_CRAY)
      REAL XMONE, WORK(LWORK)
#else
      DOUBLE PRECISION XMONE, WORK(LWORK)
#endif
!
      PARAMETER(XMONE= -1.0D00)
!
      CHARACTER*8 CFIL,DFIL,CDFIL
!
      CALL QENTER('CCRHS3_CD')
!
      ISYDIS = MULD2H(ISYMD,ISYMOP)
      ISYAIK = MULD2H(ISYDIS,ISYMPC)
!
!--------------------------
!     Allocation.
!--------------------------
!
      KSCR1 = 1
      KSCR2 = KSCR1 + NT2BCD(ISYAIK)
      KEND1 = KSCR2 + NT2BCD(ISYAIK)
      LWRK1 = LWORK - KEND1
!
      IF (LWRK1 .LE. 0 ) THEN
         CALL QUIT('Too little workspace in CCRHS3_CD ')
      ENDIF
!
      IOFF = IT2DEL(IDEL) + 1
!
      IF (NT2BCD(ISYAIK) .GT. 0) THEN
         CALL GETWA2(LUD,DFIL,WORK(KSCR1),IOFF,NT2BCD(ISYAIK))
         CALL GETWA2(LUC,CFIL,WORK(KSCR2),IOFF,NT2BCD(ISYAIK))
      ENDIF
!
!-------------------------------------
!     Calculate the contribution.
!-------------------------------------
!
         CALL DAXPY(NT2BCD(ISYAIK),XMONE,WORK(KSCR2),1,WORK(KSCR1),1)
!
!--------------------------------------
!     Save the new intermediate.
!--------------------------------------
!
      IF (NT2BCD(ISYAIK) .GT. 0) THEN
         CALL PUTWA2(LUCD,CDFIL,WORK(KSCR1),IOFF,NT2BCD(ISYAIK))
      ENDIF
!
      CALL QEXIT('CCRHS3_CD')
!
      RETURN
      END
C  /* Deck ccrhs3_prcd */
      SUBROUTINE CCRHS3_PRCD(LUD,DFIL,IDEL,WORK,LWORK,ISYMD,ISYMPC)
!
!     Written by Kasper Hald
!
!     Purpose : Prints the different elements of the
!               C or D intermediates.
!
!
      IMPLICIT NONE
!
#include <priunit.h>
#include <ccsdsym.h>
#include <maxorb.h>
#include <ccorb.h>
#include <ccsdio.h>
!
      INTEGER LWORK, LUD, IOFF, KSCR1, KEND1, LWRK1, ISYMPC
      INTEGER ISYMD, ISYDIS, ISYAIK, IDEL, IERRD
!
#if defined (SYS_CRAY)
      REAL WORK(LWORK)
#else
      DOUBLE PRECISION WORK(LWORK)
#endif
!
      CHARACTER*8 DFIL
!    
      CALL QENTER('CCRHS3_PRCD')
!
      ISYDIS = MULD2H(ISYMD,ISYMOP)
      ISYAIK = MULD2H(ISYDIS,ISYMPC)
!
!--------------------------------
!     Allocation.
!--------------------------------
!
      KSCR1 = 1
      KEND1 = KSCR1 +NT2BCD(ISYAIK)
      LWRK1 = LWORK - KEND1
!
      IF (LWRK1 .LE. 0) THEN
          CALL QUIT('Too little workspace in CCRHS3_PRCD ')
      ENDIF
!
      IOFF = IT2DEL(IDEL) + 1
!
      IF (NT2BCD(ISYAIK) .GT. 0) THEN
         CALL GETWA2(LUD,DFIL,WORK(KSCR1),IOFF,NT2BCD(ISYAIK))
      ENDIF
!
!---------------------------------------
!     Print the C/D intermediate.
!---------------------------------------
!
      IF (NT2BCD(ISYAIK) .GT. 0) THEN
         WRITE(LUPRI,*) ' The elements of ',DFIL
!
         DO I=1,NT2BCD(ISYAIK)
!
            WRITE(LUPRI,*) 'Element : ',WORK(KSCR1+I-1)
!
         ENDDO
!
      ENDIF
!
      CALL QEXIT('CCRHS3_PRCD')
!
      RETURN
      END
C  /* Deck ccrhs3_ei */
      SUBROUTINE CCRHS3_EI(DSRHF,EMAT1,EMAT2,T2AM,SCRM,XLAMDP,
     *                     XLAMDH,WORK,LWORK,IDEL,ISYMD,ISYDIS,
     *                     ISYMTR,FACE1,FACE2)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     Written by Henrik Koch 12-Jan-1994
C     Symmetry 2-aug
C     Modified slightly by Ove Christiansen 31-1-95 for 
C     construction of linear transformation intermediates. 
C     ISYMTR = SYM OF T2-VEC
!     Kasper Hald : General factor for E1 and E2 (FACE1 and FACE2)
C
C     Purpose: Calculate E-intermediates.
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
      IMPLICIT NONE
!
      INTEGER LWORK, KSCR1, KSCR2, KSCR3, KEND1, LWRK1, ISYDIS
      INTEGER IDEL, ISYMTR, ISYMD
!
#if defined (SYS_CRAY)
      REAL ONE, TWO, FACE1, FACE2
      REAL WORK(LWORK), XLAMDP(*), XLAMDH(*)
      REAL EMAT1(*), EMAT2(*), DSRHF(*), T2AM(*), SCRM(*)
#else 
      DOUBLE PRECISION ONE, TWO,FACE1, FACE2
      DOUBLE PRECISION WORK(LWORK), XLAMDP(*), XLAMDH(*)
      DOUBLE PRECISION EMAT1(*), EMAT2(*), DSRHF(*), T2AM(*)
      DOUBLE PRECISION SCRM(*)
#endif
!
      PARAMETER (ONE = 1.0D00, TWO = 2.0D00)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
!
      CALL QENTER('CCRHS3_EI')
!
!------------------------
!     Dynamic allocation.
!------------------------
!
      KSCR1  = 1
      KSCR2  = KSCR1  + NT2BCD(ISYDIS)
      KSCR3  = KSCR2  + NT2BGD(ISYDIS)
      KEND1  = KSCR3  + NT2BGD(ISYDIS)
      LWRK1  = LWORK  - KEND1
!
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Need : ',KEND1,'Available : ',LWORK
         CALL QUIT('Insufficient space in CCRHS3_EI')
      ENDIF
!
!--------------------------------
!     Calculate the contribution.
!--------------------------------
!
      CALL CCRHS3_EI1(DSRHF,EMAT1,EMAT2,T2AM,SCRM,
     *               WORK(KSCR1),WORK(KSCR2),WORK(KSCR3),
     *               XLAMDP,XLAMDH,WORK(KEND1),LWRK1,IDEL,
     *               ISYMD,ISYDIS,ISYMTR,FACE1,FACE2)
!
      CALL QEXIT('CCRHS3_EI')
!
      RETURN
      END
      SUBROUTINE CCRHS3_EI1(DSRHF,EMAT1,EMAT2,T2AM,SCRM,SCR1,SCR2,
     *                     SCR3,XLAMDP,XLAMDH,WORK,LWORK,IDEL,
     *                     ISYMD,ISYDIS,ISYMTR,FACE1,FACE2)
!
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     Written by Henrik Koch 12-Jan-1994
C     Symmetry 2-aug
C     Modified slightly by Ove Christiansen 31-1-95 for 
C     construction of linear transformation intermediates. 
C     ISYMTR = SYM OF T2-VEC
!     Kasper Hald : General factor for E1 and E2  (FACE1 and FACE2)
C
C     Purpose: Calculate E-intermediates.
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
      IMPLICIT NONE
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
!
      INTEGER LWORK, IDEL, ISYMD, ISYDIS, ISYMTR, KBM
      INTEGER ISYMJ, ISYMDJ, ISYMEM, ISYMGM, ISYME, NVIRE
      INTEGER ISYMK, NT1GM, NRHFK, KOFF1, KOFF2, KOFF3, KOFF4
      INTEGER KOFF5, KOFF6, ISYMBM, ISYMB, NT1DL, ISYMM, ISYMAG
      INTEGER ISYMDL, ISYMGL, KSCR1, KEND1, LWRK1, ISYML
      INTEGER ISYMD1, ISYMA, ISYMG, NBASA, NBASG, NVIRD, INDEX
      
!
#if defined (SYS_CRAY) 
      REAL ZERO, HALF, ONE, TWO, FACE1, FACE2
      REAL WORK(LWORK), XLAMDP(*),XLAMDH(*), DSRHF(*)
      REAL EMAT1(*), EMAT2(*), T2AM(*), SCRM(*), SCR1(*)
      REAL SCR2(*), SCR3(*)
#else
      DOUBLE PRECISION ZERO, HALF, ONE, TWO, FACE1, FACE2
      DOUBLE PRECISION WORK(LWORK), XLAMDP(*), XLAMDH(*)
      DOUBLE PRECISION DSRHF(*), EMAT1(*), EMAT2(*), T2AM(*)
      DOUBLE PRECISION SCRM(*), SCR1(*), SCR2(*), SCR3(*)
#endif
!
      PARAMETER(ZERO=0.0D00,HALF=0.5D00,ONE=1.0D00,TWO=2.0D00)
!
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
!
      CALL QENTER('CCRHS3_EI1')
!
!===================================
!     First intermediate I(b,delta).
!===================================
!
!-------------------------------------------------------
!     Construct the integrals I(dl,m) = (l d | m delta).
!-------------------------------------------------------
!
      DO 100 ISYMM = 1,NSYM
!
         ISYMAG = MULD2H(ISYMM,ISYDIS)
         ISYMDL = ISYMAG
         ISYMGL = ISYMAG
!
         DO 110 M = 1,NRHF(ISYMM)
!
            KSCR1 = 1
            KEND1 = KSCR1 + N2BST(ISYMAG)
            LWRK1 = LWORK - KEND1
            IF (LWRK1. LT. 0) THEN
               CALL QUIT('Insufficient core in CCRHS_EI1')
            END IF
!
            KOFF1 = IDSRHF(ISYMAG,ISYMM)+NNBST(ISYMAG)*(M-1)+1
            CALL CCSD_SYMSQ(DSRHF(KOFF1),ISYMAG,WORK(KSCR1))
!
            DO 120 ISYML = 1,NSYM
!
               ISYMD1 = MULD2H(ISYML,ISYMDL)
               ISYMA  = ISYML
               ISYMG  = ISYMD1
!
               NBASA = MAX(NBAS(ISYMA),1)
               NBASG = MAX(NBAS(ISYMG),1)
               NVIRD = MAX(NVIR(ISYMD1),1)
!
               KOFF2 = KSCR1 + IAODIS(ISYMA,ISYMG)
               KOFF3 = ILMRHF(ISYML) + 1
               KOFF4 = IT2BGD(ISYMGL,ISYMM) + NT1AO(ISYMGL)*(M - 1)
     *               + IT1AO(ISYMG,ISYML) + 1
!
               CALL DGEMM('T','N',NBAS(ISYMG),NRHF(ISYML),
     *                    NBAS(ISYMA),ONE,WORK(KOFF2),NBASA,
     *                    XLAMDP(KOFF3),NBASA,ZERO,SCR2(KOFF4),
     *                    NBASG)
!
               KOFF5 = ILMVIR(ISYMD1) + 1
               KOFF6 = IT2BCD(ISYMDL,ISYMM) + NT1AM(ISYMDL)*(M - 1)
     *               + IT1AM(ISYMD1,ISYML) + 1
!
               CALL DGEMM('T','N',NVIR(ISYMD1),NRHF(ISYML),
     *                    NBAS(ISYMG),ONE,XLAMDH(KOFF5),NBASG,
     *                    SCR2(KOFF4),NBASG,ZERO,SCR1(KOFF6),NVIRD)
!
  120       CONTINUE
!
  110    CONTINUE
!
  100 CONTINUE
!
!-------------------------------------------------------
!     Contract the integrals in SCR1 with t2 amplitudes.
!-------------------------------------------------------
!
      DO 200 ISYMM = 1,NSYM
!
         ISYMDL = MULD2H(ISYMM,ISYDIS)
         ISYMBM = MULD2H(ISYMDL,ISYMTR)
         ISYMB  = MULD2H(ISYMBM,ISYMM)
!
         DO 210 M = 1,NRHF(ISYMM)
!
            NT1DL = MAX(NT1AM(ISYMDL),1)
!
            KBM   = IT1AM(ISYMB,ISYMM) + NVIR(ISYMB)*(M - 1) + 1
            KOFF1 = IT2SQ(ISYMDL,ISYMBM)
     *            + NT1AM(ISYMDL)*(KBM - 1) + 1
            KOFF2 = IT2BCD(ISYMDL,ISYMM)
     *            + NT1AM(ISYMDL)*(M - 1) + 1
            KOFF3 = IEMAT1(ISYMB,ISYMD)
     *            + (IDEL - IBAS(ISYMD) - 1)*NVIR(ISYMB) + 1
!
            CALL DGEMV('T',NT1AM(ISYMDL),NVIR(ISYMB),FACE1,
     *                 T2AM(KOFF1),NT1DL,SCR1(KOFF2),1,ONE,
     *                 EMAT1(KOFF3),1)
!
  210    CONTINUE
!
  200 CONTINUE
!
!================================
!     Second intermediate I(k,j).
!================================
!
!-------------------------------------------
!     Transform the SCRM amplitudes to SCR3.
!-------------------------------------------
!
      DO 300 ISYMJ = 1,NSYM
!
         ISYMDJ = MULD2H(ISYMD,ISYMJ)
         ISYMEM = MULD2H(ISYMDJ,ISYMTR)
         ISYMGM = ISYMEM
!
         DO 310 J = 1,NRHF(ISYMJ)
!
            DO 320 ISYMM = 1,NSYM
!
               ISYME = MULD2H(ISYMM,ISYMEM)
               ISYMG = ISYME
!
               NVIRE = MAX(NVIR(ISYME),1)
               NBASG = MAX(NBAS(ISYMG),1)
!
               KOFF1 = ILMVIR(ISYME) + 1
               KOFF2 = IT2BCD(ISYMEM,ISYMJ) + NT1AM(ISYMEM)*(J - 1)
     *               + IT1AM(ISYME,ISYMM) + 1
               KOFF3 = IT2BGD(ISYMGM,ISYMJ) + NT1AO(ISYMGM)*(J - 1)
     *               + IT1AO(ISYMG,ISYMM) + 1
!
               CALL DGEMM('N','N',NBAS(ISYMG),NRHF(ISYMM),
     *                    NVIR(ISYME),ONE,XLAMDH(KOFF1),
     *                    NBASG,SCRM(KOFF2),NVIRE,ZERO,
     *                    SCR3(KOFF3),NBASG)
!
  320       CONTINUE
  310    CONTINUE
  300 CONTINUE
!
!----------------------------------------------------------------
!     Contract the integrals in SCR2 with the amplitudes in SCR3.
!----------------------------------------------------------------
!
      DO 400 ISYMJ = 1,NSYM
!
         ISYMDJ = MULD2H(ISYMD,ISYMJ)
         ISYMEM = MULD2H(ISYMDJ,ISYMTR)
         ISYMGM = ISYMEM
         ISYMK  = MULD2H(ISYMGM,ISYDIS)
!
         NT1GM = MAX(NT1AO(ISYMGM),1)
         NRHFK = MAX(NRHF(ISYMK),1)
!
         KOFF1 = IT2BGD(ISYMGM,ISYMK) + 1
         KOFF2 = IT2BGD(ISYMGM,ISYMJ) + 1
         KOFF3 = IMATIJ(ISYMK,ISYMJ) + 1
!
         CALL DGEMM('T','N',NRHF(ISYMK),NRHF(ISYMJ),NT1AO(ISYMGM),
     *              FACE2,SCR2(KOFF1),NT1GM,SCR3(KOFF2),NT1GM,
     *              ONE,EMAT2(KOFF3),NRHFK)
!
  400 CONTINUE
!
      CALL QEXIT('CCRHS3_EI1')
!
      RETURN
      END
C  /* Deck cc_aofock3 */
      SUBROUTINE CC_AOFOCK3(XINT,DENSIT,FOCK,WORK,LWORK,IDEL,
     *                      ISYMD,ISYDEN)
C
C     Written by Asger Halkier and Henrik Koch 27-4-95.
C
C     Debugged Ove Christiansen august 1995
C
C     Purpose: Calculate the two electron contribution to the
C              AO-fock matrix using matrix vector routines.
C
C     Obs: It can be done as F(g>=d) = G(a>=b) I(a>=b,g,d) where
C          G(a>=b) = D(a,b) + D(b,a), the diagonal properly scaled
C
      IMPLICIT NONE
      INTEGER ISYDIS, ISYMD, ISYMG, ISYMAB, NDISTG, NBATCH
      INTEGER IBATCH, NUMG, IG1, IG2, KOFF2, IG, KOFF1, ISYMA
      INTEGER ISYDEN, ISYMB, KAD, IDEL, KGB, NTOTA, NTOTG, LWORK
#include <priunit.h>
#include <maxorb.h>
#if defined (SYS_CRAY)
      REAL ZERO, ONE, TWO
      REAL XINT(*), DENSIT(*), FOCK(*), WORK(LWORK)
#else
      DOUBLE PRECISION ZERO, ONE, TWO
      DOUBLE PRECISION XINT(*), DENSIT(*), FOCK(*), WORK(LWORK)
#endif
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
#include <ccorb.h>
#include <symsq.h>
#include <ccsdsym.h>
C
      CALL QENTER('CC_AOFOCK3')
C
      ISYDIS = MULD2H(ISYMD,ISYMOP)
C
      DO 100 ISYMG = 1,NSYM
C
         IF (NBAS(ISYMG) .EQ. 0) GOTO 100
C
         ISYMAB = MULD2H(ISYMG,ISYDIS)
C
         NDISTG = MIN(LWORK/N2BST(ISYMAB),NBAS(ISYMG))
C
         IF (NDISTG .LT. 1) THEN
            CALL QUIT('Insufficient work space in CC_AOFOCK1')
         ENDIF
C
         NBATCH = (NBAS(ISYMG) - 1)/NDISTG + 1
C
C-------------------------------------
C        Start the loops over batches.
C-------------------------------------
C
         DO 110 IBATCH = 1,NBATCH
C
            NUMG = NDISTG
            IF (IBATCH .EQ. NBATCH) THEN
               NUMG = NBAS(ISYMG) - NDISTG*(NBATCH - 1)
            ENDIF
C
            IG1 = NDISTG*(IBATCH - 1) + 1
            IG2 = NDISTG*(IBATCH - 1) + NUMG
C
            KOFF2 = 1
            DO 120 IG = IG1,IG2
C
               KOFF1 = IDSAOG(ISYMG,ISYDIS)
     *               + (IG - 1)*NNBST(ISYMAB) + 1
C
               CALL CCSD_SYMSQ(XINT(KOFF1),ISYMAB,WORK(KOFF2))
C
               KOFF2 = KOFF2 + N2BST(ISYMAB)
C
  120       CONTINUE
C
            ISYMA = MULD2H(ISYMD,ISYDEN)
            ISYMB = MULD2H(ISYMA,ISYMAB)
C
            KAD = IAODIS(ISYMA,ISYMD)
     *          + NBAS(ISYMA)*(IDEL - IBAS(ISYMD) - 1) + 1
C
            DO 130 IG = IG1,IG2
C
               KOFF1 = (IG - IG1)*N2BST(ISYMAB)
     *               + IAODIS(ISYMA,ISYMB) + 1
               KGB   = IAODIS(ISYMG,ISYMB) + IG
C
               NTOTA = MAX(NBAS(ISYMA),1)
               NTOTG = MAX(NBAS(ISYMG),1)
C
               CALL DGEMV('T',NBAS(ISYMA),NBAS(ISYMB),-ONE,WORK(KOFF1),
     *                    NTOTA,DENSIT(KAD),1,ONE,FOCK(KGB),NTOTG)
C
  130       CONTINUE
C
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('CC_AOFOCK3')
C
      RETURN
      END
C  /* Deck cc_mofcon3 */
      SUBROUTINE CC_MOFCON3(XINT,OMEGA2,XLAMDP,XLAMDH,XLAMPC,XLAMHC,
     *                      WORK,LWORK,IDEL,ISYMD,ISYMTR,IOPT,
     *                      ANTISYM)
C
C     Written by Asger Halkier and Henrik Koch 3-5-95.
C     Debugged By Ove Christiansen 25-7-1995
C     ANTISYM flag introduced K. Hald & C. Haettig August 99
C
C     Purpose: To calculate the F-term's contribution to the
C              vector function using matrix vector routines.
C              Special version adapted for triplet case.
C
C     N.B. This routine assumes AO-symmetric integrals, and can therefor
C          not be used directly for calculations with London-orbitals!!!
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
      PARAMETER(ZERO = 0.0D0,ONE = 1.0D0,XMONE=-1.0D0,TWO = 2.0D0)
      DIMENSION XINT(*),OMEGA2(*)
      DIMENSION XLAMPC(*),XLAMHC(*),XLAMDH(*),XLAMDP(*)
      DIMENSION WORK(LWORK)
      LOGICAL ANTISYM
#include <ccorb.h>
#include <symsq.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CC_MOFCON3')
C
      ISYDIS = MULD2H(ISYMD,ISYMOP)
C
      DO 100 ISYMG = 1,NSYM
C
         IF (NBAS(ISYMG) .EQ. 0) GOTO 100
C
         ISALBE = MULD2H(ISYMG,ISYDIS)
         ISYMAI = MULD2H(ISALBE,ISYMTR)
         ISYMJ  = ISYMG
C
C-----------------------------------------
C        Dynamic allocation of work space.
C-----------------------------------------
C
         KSCR1 = 1
         KSCR2 = KSCR1 + NNBST(ISALBE)*NRHF(ISYMJ)
         KSCR3 = KSCR2 + N2BST(ISALBE)
         KSCR4 = KSCR3 + NT1AM(ISYMAI)
         KEND1 = KSCR4 + NT1AM(ISYMAI)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Lwrk1 = ',LWRK1
            CALL QUIT('Insufficient work space area in CC_MOFCON')
         ENDIF
C
C--------------------------------
C        Do first transformation.
C--------------------------------
C
         KOFF1 = IDSAOG(ISYMG,ISYDIS) + 1
         KOFF2 = ILMRHF(ISYMJ) + 1
C
         NTALBE = MAX(NNBST(ISALBE),1)
         NTOTG  = MAX(NBAS(ISYMG),1)
C
         CALL DGEMM('N','N',NNBST(ISALBE),NRHF(ISYMJ),NBAS(ISYMG),
     *              ONE,XINT(KOFF1),NTALBE,XLAMDH(KOFF2),NTOTG,
     *              ZERO,WORK(KSCR1),NTALBE)
C
C-----------------------------------
C        Last index transformations.
C-----------------------------------
C
         DO 110 J = 1,NRHF(ISYMJ)
C
            KOFF1 = KSCR1 + NNBST(ISALBE)*(J - 1)
C
            CALL CCSD_SYMSQ(WORK(KOFF1),ISALBE,WORK(KSCR2))
C
            DO 120 ISYMI = 1,NSYM
C
               ISYMBE = ISYMI
               ISYMAL = MULD2H(ISYMBE,ISALBE)
               ISYMA  = MULD2H(ISYMAL,ISYMTR)
C
               IF (LWRK1 .LT. NBAS(ISYMAL)*NRHF(ISYMI)) THEN
                  CALL QUIT('Insufficient space for '//
     &                 '2. trf. in CC_MOFCON')
               ENDIF
C
               KOFF2 = KSCR2 + IAODIS(ISYMAL,ISYMBE)
               KOFF3 = ILMRHF(ISYMI) + 1
               KOFF4 = IGLMVI(ISYMAL,ISYMA) + 1
               KOFF5 = KSCR3 + IT1AM(ISYMA,ISYMI)
C
               NTOTAL = MAX(NBAS(ISYMAL),1)
               NTOTBE = MAX(NBAS(ISYMBE),1)
               NTOTA  = MAX(NVIR(ISYMA),1)
C
               CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMI),NBAS(ISYMBE),
     *                    ONE,WORK(KOFF2),NTOTAL,XLAMDH(KOFF3),NTOTBE,
     *                    ZERO,WORK(KEND1),NTOTAL)
C
               CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMAL),
     *                    ONE,XLAMPC(KOFF4),NTOTAL,WORK(KEND1),NTOTAL,
     *                    ZERO,WORK(KOFF5),NTOTA)
C
               IF (IOPT .EQ. 2) THEN
C
                  ISYMBE = MULD2H(ISYMI,ISYMTR)
                  ISYMAL = MULD2H(ISYMBE,ISALBE)
                  ISYMA  = ISYMAL
C
                  IF (LWRK1 .LT. NBAS(ISYMAL)*NRHF(ISYMI)) THEN
                     CALL QUIT('Insufficient space for '//
     &                    '2. trf. in CC_MOFCON')
                  ENDIF
C
                  KOFF2 = KSCR2 + IAODIS(ISYMAL,ISYMBE)
                  KOFF3 = IGLMRH(ISYMBE,ISYMI) + 1
                  KOFF4 = ILMVIR(ISYMA) + 1
                  KOFF5 = KSCR3 + IT1AM(ISYMA,ISYMI)
C
                  NTOTAL = MAX(NBAS(ISYMAL),1)
                  NTOTBE = MAX(NBAS(ISYMBE),1)
                  NTOTA  = MAX(NVIR(ISYMA),1)
C
                  CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMI),
     *                       NBAS(ISYMBE),ONE,WORK(KOFF2),NTOTAL,
     *                       XLAMHC(KOFF3),NTOTBE,ZERO,WORK(KEND1),
     *                       NTOTAL)
C
                  CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),
     *                       NBAS(ISYMAL),ONE,XLAMDP(KOFF4),NTOTAL,
     *                       WORK(KEND1),NTOTAL,ONE,WORK(KOFF5),NTOTA)
C
               ENDIF
C

  120       CONTINUE
C
C--------------------------------------------------
C           Storing the result in the omega2-array.
C--------------------------------------------------
C
            ISYMB  = ISYMD
            ISYMBJ = MULD2H(ISYMB,ISYMJ)
C
            DO 130 B = 1,NVIR(ISYMB)
C
               NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
               NDB = ILMVIR(ISYMB) + NBAS(ISYMD)*(B - 1)
     *             + IDEL - IBAS(ISYMD)
C
               CALL DZERO(WORK(KSCR4),NT1AM(ISYMAI))
C
               XLB  = XLAMDP(NDB)
C
               CALL DAXPY(NT1AM(ISYMAI),XLB,WORK(KSCR3),1,WORK(KSCR4),1)
C
               IF (ISYMBJ .EQ. ISYMAI) THEN
C
                  NTOTAI = NBJ
C
                  IF (.NOT. ANTISYM) THEN
!
                  IF (IOPT .EQ. 2) THEN
                     NTOTAI = NT1AM(ISYMAI)
                     WORK(KSCR4+NBJ-1) = TWO*WORK(KSCR4+NBJ-1)
                  ENDIF
C
                  DO 140 NAI = 1,NTOTAI
C
                     NAIBJ = IT2AM(ISYMAI,ISYMBJ) + INDEX(NAI,NBJ)
C
                     OMEGA2(NAIBJ) = OMEGA2(NAIBJ) + WORK(KSCR4+NAI-1)
C
  140             CONTINUE
                  ELSE
                  IF (IOPT .EQ. 1) CALL QUIT( 
     *               'IOPT .EQ. 1 .AND. ANTISYM in MOFCON3 not legal')
!
                  DO NAI = 1,NBJ-1
                     NAIBJ = IT2AM(ISYMAI,ISYMBJ) + INDEX(NAI,NBJ)
                     OMEGA2(NAIBJ) = OMEGA2(NAIBJ) + WORK(KSCR4+NAI-1)
                  ENDDO
C
                  DO NAI = NBJ+1,NT1AM(ISYMAI)
                     NAIBJ = IT2AM(ISYMAI,ISYMBJ) + INDEX(NAI,NBJ)
                     OMEGA2(NAIBJ) = OMEGA2(NAIBJ) - WORK(KSCR4+NAI-1)
                  ENDDO
!
                  ENDIF
C
               ENDIF
C
               IF (ISYMAI .LT. ISYMBJ) THEN
C
                 IF (.NOT. ANTISYM) THEN
                  DO NAI = 1,NT1AM(ISYMAI)
                     NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                     + NT1AM(ISYMAI)*(NBJ - 1) + NAI
                     OMEGA2(NAIBJ) = OMEGA2(NAIBJ) + WORK(KSCR4+NAI-1)
                  END DO
                 ELSE
                  DO NAI = 1,NT1AM(ISYMAI)
                     NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                     + NT1AM(ISYMAI)*(NBJ - 1) + NAI
                     OMEGA2(NAIBJ) = OMEGA2(NAIBJ) + WORK(KSCR4+NAI-1)
                  END DO
                 END IF
C
               ENDIF
C
               IF ((ISYMBJ .LT. ISYMAI) .AND. (IOPT .EQ. 2)) THEN
C
                 IF (.NOT.ANTISYM) THEN
                  DO NAI = 1,NT1AM(ISYMAI)
                     NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                     + NT1AM(ISYMBJ)*(NAI - 1) + NBJ
                     OMEGA2(NAIBJ) = OMEGA2(NAIBJ) + WORK(KSCR4+NAI-1)
                  END DO
                 ELSE
                  DO NAI = 1,NT1AM(ISYMAI)
                     NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                     + NT1AM(ISYMBJ)*(NAI - 1) + NBJ
                     OMEGA2(NAIBJ) = OMEGA2(NAIBJ) - WORK(KSCR4+NAI-1)
                  END DO
                 END IF
C
               ENDIF
C
  130       CONTINUE
C
  110    CONTINUE
C
  100 CONTINUE
C
      CALL QEXIT('CC_MOFCON3')
C
      RETURN
      END
C  /* Deck cc_pram3 */
      SUBROUTINE CC_PRAM3(CAM1,CAMP,CAMM,PT1,PTP,PTM,ISYMTR,LGRS)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     30-5-1995 Ove Christiansen
C     05-8-1999 Kasper Hald & Christof Haettig : adapted for triplet
C
C     Purpose: Writes out vector:
C              %T1 and %T2 and ||T1||/||T2||
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include <implicit.h>
C
      PARAMETER (TWO = 2.0D0, THPRT = 1.0D-9)
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
C
C
      LOGICAL LGRS
      DIMENSION CAM1(*), CAMP(*), CAMM(*)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CC_PRAM3')
C
C------------------------
C     Add up the vectors.
C------------------------
C
      C1NOSQ = DDOT(NT1AM(ISYMTR),CAM1,1,CAM1,1)
      IF (.NOT. CCS) THEN
        C2PNOSQ = DDOT(NT2AM(ISYMTR),CAMP,1,CAMP,1)
        C2MNOSQ = DDOT(NT2AM(ISYMTR),CAMM,1,CAMM,1)
      ELSE
        C2PNOSQ = 0.0D0
        C2MNOSQ = 0.0D0
      END IF
C
      CNOSQ  = C1NOSQ + C2PNOSQ + C2MNOSQ
C
      IF (CNOSQ .EQ. 0.0D0) THEN
         PT1 = 0.0D0
         PTP = 0.0D0
         PTM = 0.0D0
      ELSE
         PT1 = (C1NOSQ /CNOSQ)*100.0D0
         PTP = (C2PNOSQ/CNOSQ)*100.0D0
         PTM = (C2MNOSQ/CNOSQ)*100.0D0
      END IF
C
      IF (.NOT. CCS .AND. CNOSQ .NE. 0.0D0) THEN
        WRITE(LUPRI,'(//5X,A)')
     *     'CC_PRAM:Overall Contribution of the Different Components'
        WRITE(LUPRI,'(5X,A//)')
     *     '--------------------------------------------------------'
        WRITE(LUPRI,'(/5X,A,5X,F10.4,A)')
     *     'Single Excitation Contribution      : ', PT1,' %'
        WRITE(LUPRI,'(/5X,A,5X,F10.4,A,F10.4,A)')
     *     'Double Excitation Contribution (+/-): ',
     *     PTP,' %    /',PTM,' % '
        WRITE(LUPRI,'(/5X,A,5X,F10.4)')
     *     '||T1||/||T2||                       : ',
     *      SQRT(C1NOSQ/(C2PNOSQ+C2MNOSQ))
        IF (LGRS) WRITE(LUPRI,'(/5X,A,5X,F10.4)')
     *     'Tau1 diagnostic                     : ',
     *      SQRT(C1NOSQ/(TWO*DFLOAT(NRHFT)))
      END IF
C
      WRITE(LUPRI,'(/5X,A,5X,F10.4)')
     *  'Norm of Total Amplitude Vector : ',SQRT(CNOSQ)
C
      CALL FLSHFO(LUPRI)
C
C----------------------------------------------
C     Initialize threshold etc from Printlevel.
C----------------------------------------------
C
      NL = MAX(1,2*IPRINT)
      CNOSQ = MAX(CNOSQ,THPRT)
      THR1 = SQRT(C1NOSQ/CNOSQ)/NL
      THR2 = SQRT((C2PNOSQ+C2MNOSQ)/CNOSQ)/NL
      THR1 = MAX(THR1,1.0D-03)
      THR2 = MAX(THR2,1.0D-03)
      SUMOFP = 0.0D00
C
C---------------------------------------
C     Loop through single excitation part.
C---------------------------------------
C
      WRITE(LUPRI,'(//A)')
     *     ' +=============================================='
     *    //'===============================+'
      WRITE(LUPRI,'(1X,A)')
     *     '| symmetry|  orbital index  |   Excitation Numbers'
     *     //'             |   Amplitude  |'
      WRITE(LUPRI,'(1X,A)')
     *     '|  Index  |   a   b   i   j |      NAI      NBJ |'
     *     //'     NAIBJ    |              |'
      WRITE(LUPRI,'(A)')
     *     ' +=============================================='
     *    //'===============================+'
C
      ISYMAI = MULD2H(ISYMTR,ISYMOP)
C
  1   CONTINUE
      N1 = 0
C
      DO 100 ISYMA = 1,NSYM
C
         ISYMI = MULD2H(ISYMAI,ISYMA)
C
         DO 110 I = 1,NRHF(ISYMI)
C
            MI = IORB(ISYMI) + I
C
            DO 120 A=1,NVIR(ISYMA)
C
               NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
C
               MA = IORB(ISYMA) + NRHF(ISYMA) +  A
C
               IF (ABS(CAM1(NAI)) .GT. THR1 ) THEN
C
                  WRITE(LUPRI,9990) ISYMA,ISYMI,A,I,NAI,CAM1(NAI)
C
                  N1 = N1 + 1
                  SUMOFP = SUMOFP + CAM1(NAI)*CAM1(NAI)
C
               ENDIF
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      IF ((N1.LT.1).AND.(SQRT(C1NOSQ/CNOSQ).GT.1.0D-3)) THEN
         THR1 = THR1/5.0D0
         GOTO 1
      ENDIF
C
      CALL FLSHFO(LUPRI)
C
C--------------------------------------------
C     Loop through Doublee excitation vector.
C     If not ccs or ccp2
C--------------------------------------------
C
      IF (.NOT. ( CCS .OR. CCP2 )) THEN
C
      WRITE(LUPRI,'(A)')
     *     ' +----------------------------------------------'
     *    //'-------------------------------+'
C

 2    CONTINUE
      N2 = 0
C
      DO 200 ISYMAI = 1,NSYM
C
         ISYMBJ = MULD2H(ISYMAI,ISYMTR)
C
         DO 210 ISYMJ = 1,NSYM
C
            ISYMB = MULD2H(ISYMJ,ISYMBJ)
C
            DO 220 ISYMI = 1,NSYM
C
               ISYMA = MULD2H(ISYMI,ISYMAI)
C
               DO 230 J = 1,NRHF(ISYMJ)
C
                  MJ = IORB(ISYMJ) + J
C
                  DO 240 B = 1,NVIR(ISYMB)
C
                     NBJ = IT1AM(ISYMB,ISYMJ)
     *                   + NVIR(ISYMB)*(J - 1) + B
C
                     MB = IORB(ISYMB) + NRHF(ISYMB) + B
C
                     DO 250 I = 1,NRHF(ISYMI)
C
                        MI = IORB(ISYMI) + I
C
                        DO 260 A = 1,NVIR(ISYMA)
C
                           NAI = IT1AM(ISYMA,ISYMI)
     *                         + NVIR(ISYMA)*(I - 1) + A
C
                           MA = IORB(ISYMA) + NRHF(ISYMA) +  A
C
                           IF (((ISYMAI.EQ.ISYMBJ).AND.
     *                         (NAI .LT. NBJ)).OR.(ISYMAI.LT.ISYMBJ))
     *                          GOTO 260
C
                           IF (ISYMAI.EQ.ISYMBJ) THEN
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                             + INDEX(NAI,NBJ)
                           ELSE
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                            + NT1AM(ISYMAI)*(NBJ-1) + NAI
                           ENDIF
C
                           IF (ABS(CAMP(NAIBJ)) .GT. THR2 ) THEN
                              WRITE(LUPRI,9991) ISYMA,ISYMB,ISYMI,ISYMJ,
     *                                      A,B,I,J,NAI,NBJ,NAIBJ,
     *                                      CAMP(NAIBJ)
                              N2 = N2 + 1
                              SUMOFP = SUMOFP + CAMP(NAIBJ)*CAMP(NAIBJ)
                           ENDIF
C
                           IF (ABS(CAMM(NAIBJ)) .GT. THR2 ) THEN
                              WRITE(LUPRI,9992) ISYMA,ISYMB,ISYMI,ISYMJ,
     *                                      A,B,I,J,NAI,NBJ,NAIBJ,
     *                                      CAMM(NAIBJ)
                              N2 = N2 + 1
                              SUMOFP = SUMOFP + CAMM(NAIBJ)*CAMM(NAIBJ)
                           ENDIF
C
  260                   CONTINUE
  250                CONTINUE
  240             CONTINUE
  230          CONTINUE
  220       CONTINUE
  210    CONTINUE
  200 CONTINUE
C
      IF ((N2.LT.1).AND.(SQRT((C2PNOSQ+C2MNOSQ)/CNOSQ).GT.1.0D-3)) THEN
         THR2 = THR2/5D00
         GOTO 2
      ENDIF
C
      ENDIF
C
      WRITE(LUPRI,'(A)')
     *     ' +=============================================='
     *    //'===============================+'
C
      WRITE(LUPRI,'(//10X,A,8X,F10.4)')
     *     'Norm of Printed Amplitude Vector : ',SQRT(SUMOFP)
      WRITE(LUPRI,'(//10X,A43,1X,F8.6)')
     *     'Printed all single excitations greater than',THR1
      IF (.NOT. (CCS.OR.CCP2)) THEN
         WRITE(LUPRI,'(//10X,A43,1X,F8.6)')
     *     'Printed all double excitations greater than',THR2
      ENDIF
C
 9990 FORMAT(1X,'| ',I1,3X,I1,2X,' | ',I3,5X,I3,4X,' | ',I8,9x,
     *       ' | ',12x,' | ',1x, F10.6,'  |')
 9991 FORMAT(1X,'| ',I1,1X,I1,1X,I1,1X,I1,' | ',
     *       I3,1X,I3,1X,I3,1X,I3,' | ',
     *       I8,1x,I8,' | (+)',I9,' | ',1x,F10.6,'  |')
 9992 FORMAT(1X,'| ',I1,1X,I1,1X,I1,1X,I1,' | ',
     *       I3,1X,I3,1X,I3,1X,I3,' | ',
     *       I8,1x,I8,' | (-)',I9,' | ',1x,F10.6,'  |')
C
      CALL QEXIT('CC_PRAM3')
C
      RETURN
      END
C  /* Deck ccrhs3_ij */
      SUBROUTINE CCRHS3_IJ(OMEGA2,WORK,LWORK,ISYVEC)
!
!     Written by Kasper Hald and Poul Joergensen
!     Spring 1999.
!
!     Purpose : Calculate Omega(aibj) - Omega(ajbi)
!
!     N.B. It is assumed that omega will be in packed form.
!
!
      IMPLICIT NONE
!
#include <priunit.h>
#include <maxorb.h>
#include <ccorb.h>
#include <symsq.h>
#include <ccsdsym.h>
#include <cclr.h>
#include <ccsdio.h>
!
      INTEGER LWORK, LWRK1, KEND1, KSCR1, ISYMJ, ISYMI, ISYMA
      INTEGER ISYMB, ISYMAI, ISYMAJ, ISYMBI, ISYMBJ
      INTEGER NAI, NAJ, NBI, NBJ, NAIBJ, NAJBI, NTOTA,ISYVEC
      INTEGER INDEX
      INTEGER MA, MB, MI, MJ
!
#if defined (SYS_CRAY)
      REAL WORK(LWORK), OMEGA2(*)
      REAL ZERO
#else
      DOUBLE PRECISION WORK(LWORK), OMEGA2(*)
      DOUBLE PRECISION ZERO
#endif
!
      PARAMETER(ZERO = 0.0D00)
!
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
!
      CALL QENTER('CCRHS3_IJ')
!
!-----------------------------------
!     Allocation of workspace.
!-----------------------------------
!
      KSCR1 = 1
      KEND1 = KSCR1 + NT2AM(ISYVEC)
      LWRK1 = LWORK - KEND1
!
      IF (LWRK1 .LE. 0) THEN
         CALL QUIT('Too little workspace in CCRHS3_IJ ')
      ENDIF
!
!------------------------------------------
!     Copy OMEGA to workspace.
!------------------------------------------
!
      CALL DCOPY(NT2AM(ISYVEC),OMEGA2,1,WORK(KSCR1),1)
!
!------------------------------------------
!     Calculate the contribution.
!------------------------------------------
!
      DO 100 ISYMBJ = 1,NSYM
         ISYMAI = MULD2H(ISYMBJ,ISYVEC)
!
         IF (ISYMAI .LE. ISYMBJ) THEN
!
         DO 110 ISYMI = 1,NSYM
!
            ISYMA = MULD2H(ISYMAI,ISYMI)
!
            DO 120 ISYMJ = 1,NSYM
!
               ISYMB  = MULD2H(ISYMBJ,ISYMJ)
               ISYMBI = MULD2H(ISYMB,ISYMI)
               ISYMAJ = MULD2H(ISYMA,ISYMJ)
!
               DO 130 I = 1,NRHF(ISYMI)
               MI = IORB(ISYMI) + I
                  DO 140 J = 1,NRHF(ISYMJ)
                  MJ = IORB(ISYMJ) + J
!
                        DO 150 A = 1,NVIR(ISYMA)
                        MA = IORB(ISYMA) + NRHF(ISYMA) + A
                           NAJ   = IT1AM(ISYMA,ISYMJ)
     *                           + NVIR(ISYMA)*(J-1) + A
                           NAI   = IT1AM(ISYMA,ISYMI)
     *                           + NVIR(ISYMA)*(I-1) + A
!
                           DO 160 B = 1,NVIR(ISYMB)
                           MB = IORB(ISYMB) + NRHF(ISYMB) + B
!
                              NBI   = IT1AM(ISYMB,ISYMI)
     *                              + NVIR(ISYMB)*(I-1) + B
                              NBJ   = IT1AM(ISYMB,ISYMJ)
     *                              + NVIR(ISYMB)*(J-1) + B
!
                                 IF (ISYMAI .EQ. ISYMBJ) THEN
!
                                  IF ((NAI .LE. NBJ) .AND.
     *                                (MA .LE. MB)) THEN
!
                                    NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                                    + INDEX(NAI,NBJ)
!
                                    IF (ISYMAJ .EQ. ISYMBI) THEN
                                       NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *                                       + INDEX(NAJ,NBI)
                                    ELSEIF (ISYMAJ .LT. ISYMBI) THEN
                                       NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *                                     + NT1AM(ISYMAJ)*(NBI-1)+NAJ
                                    ELSEIF (ISYMAJ .GT. ISYMBI) THEN
                                       NAJBI = IT2AM(ISYMBI,ISYMAJ)
     *                                     + NT1AM(ISYMBI)*(NAJ-1)+NBI
                                    ENDIF
!
                                       OMEGA2(NAIBJ) = OMEGA2(NAIBJ)
     *                                               - WORK(NAJBI)
                                       OMEGA2(NAJBI) = ZERO
!
                                  ENDIF
!
                                 ELSEIF (ISYMAI .LT. ISYMBJ) THEN
!
                                   IF (((MA .LE. MB) .AND.
     *                                  (MI .LE. MJ)) .OR.
     *                                 ((MA .GE. MB) .AND.
     *                                  (MI .GE. MJ))) THEN
!
                                     NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                               + NT1AM(ISYMAI)*(NBJ-1)+NAI
!
                                     IF (ISYMAJ .EQ. ISYMBI) THEN
                                        NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *                                        + INDEX(NAJ,NBI)
                                     ELSEIF (ISYMAJ .LT. ISYMBI) THEN
                                        NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *                                  + NT1AM(ISYMAJ)*(NBI-1)+NAJ
                                     ELSEIF (ISYMAJ .GT. ISYMBI) THEN
                                        NAJBI = IT2AM(ISYMBI,ISYMAJ)
     *                                  + NT1AM(ISYMBI)*(NAJ-1)+NBI
                                     ENDIF
!
                                        OMEGA2(NAIBJ) = OMEGA2(NAIBJ)
     *                                                -WORK(NAJBI)
                                        OMEGA2(NAJBI) = ZERO
!
                                   ENDIF
!
                                 ENDIF
!
  160                      CONTINUE
  150                   CONTINUE
  140                CONTINUE
  130             CONTINUE
  120          CONTINUE
  110    CONTINUE
         ENDIF
  100 CONTINUE
!
      CALL QEXIT('CCRHS3_IJ')
!
      RETURN
      END
C  /* Deck ccrhs3_r2ij */
      SUBROUTINE CCRHS3_R2IJ(C2AM,WORK,LWORK,ISYVEC)
!
!     Written by Kasper Hald.
!     Spring 1999.
!
!     Purpose : Take the (+)R(ab,ij) vector
!               for ai<bj AND i<j and "square" it up
!               to include all terms ai<bj i.e.
!               a lower triangular matrix.
!
!     N.B. It is assumed that omega will be in packed form.
!
!
      IMPLICIT NONE
!
#include <priunit.h>
#include <maxorb.h>
#include <ccorb.h>
#include <symsq.h>
#include <ccsdsym.h>
#include <cclr.h>
#include <ccsdio.h>
!
      INTEGER LWORK, LWRK1, KEND1, KSCR1, ISYMJ, ISYMI, ISYMA
      INTEGER ISYMB, ISYMAI, ISYMAJ, ISYMBI, ISYMBJ, ISYVEC
      INTEGER NAI, NAJ, NBI, NBJ, NAIBJ, NAJBI, NTOTA
      INTEGER INDEX
      INTEGER MA, MB, MI, MJ
!
#if defined (SYS_CRAY)
      REAL WORK(LWORK), C2AM(*)
      REAL ZERO
#else
      DOUBLE PRECISION WORK(LWORK), C2AM(*)
      DOUBLE PRECISION ZERO
#endif
!
      PARAMETER(ZERO = 0.0D00)
!
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
!
      CALL QENTER('CCRHS3_R2IJ')
!
!-----------------------------------
!     Allocation of workspace.
!-----------------------------------
!
      KSCR1 = 1
      KEND1 = KSCR1 + NT2AM(ISYVEC)
      LWRK1 = LWORK - KEND1
!
      IF (LWRK1 .LE. 0) THEN
         CALL QUIT('Too little workspace in CCRHS3_R2IJ ')
      ENDIF
!
!------------------------------------------
!     Copy OMEGA to workspace.
!------------------------------------------
!
      CALL DCOPY(NT2AM(ISYVEC),C2AM,1,WORK(KSCR1),1)
!
!------------------------------------------
!     Calculate the contribution.
!------------------------------------------
!
      DO 100 ISYMBJ = 1,NSYM
         ISYMAI = MULD2H(ISYMBJ,ISYMTR)
!
         IF (ISYMAI .LE. ISYMBJ) THEN
!
            DO 110 ISYMI = 1,NSYM
!
               ISYMA = MULD2H(ISYMAI,ISYMI)
!
               DO 120 ISYMJ = 1,NSYM
!
                  ISYMB  = MULD2H(ISYMBJ,ISYMJ)
                  ISYMAJ = MULD2H(ISYMA,ISYMJ)
                  ISYMBI = MULD2H(ISYMB,ISYMI)
!
                  DO 130 I = 1,NRHF(ISYMI)
                  MI = IORB(ISYMI) + I
                     DO 140 J = 1,NRHF(ISYMJ)
                     MJ = IORB(ISYMJ) + J
!
                     IF (MI .NE. MJ) THEN
                        DO 150 A = 1,NVIR(ISYMA)
                        MA = IORB(ISYMA) + NRHF(ISYMA) + A
                           NAJ   = IT1AM(ISYMA,ISYMJ)
     *                           + NVIR(ISYMA)*(J-1) + A
                           NAI   = IT1AM(ISYMA,ISYMI)
     *                           + NVIR(ISYMA)*(I-1) + A
!
                           DO 160 B = 1,NVIR(ISYMB)
                           MB = IORB(ISYMB) + NRHF(ISYMB) + B
!
                           IF (MA .NE. MB) THEN
!
                              NBI   = IT1AM(ISYMB,ISYMI)
     *                              + NVIR(ISYMB)*(I-1) + B
                              NBJ   = IT1AM(ISYMB,ISYMJ)
     *                              + NVIR(ISYMB)*(J-1) + B
!
                                 IF (ISYMAI .EQ. ISYMBJ) THEN
!
                                  IF ((NAI .LT. NBJ) .AND.
     *                                (MA .LT. MB)) THEN
!
                                    NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                                    + INDEX(NAI,NBJ)
!
                                    IF (ISYMAJ .EQ. ISYMBI) THEN
                                       NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *                                       + INDEX(NAJ,NBI)
                                    ELSEIF (ISYMAJ .LT. ISYMBI) THEN
                                       NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *                                     + NT1AM(ISYMAJ)*(NBI-1)+NAJ
                                    ELSEIF (ISYMAJ .GT. ISYMBI) THEN
                                       NAJBI = IT2AM(ISYMBI,ISYMAJ)
     *                                     + NT1AM(ISYMBI)*(NAJ-1)+NBI
                                    ENDIF
!
                                    C2AM(NAJBI) = - WORK(NAIBJ)
                                    C2AM(NAIBJ) =   WORK(NAIBJ)
!
                                  ENDIF
!
                                 ELSEIF (ISYMAI .LT. ISYMBJ) THEN
!
                                   IF (((MA .LT. MB) .AND. 
     *                                  (MI .LT. MJ)) .OR.
     *                                 ((MA .GT. MB) .AND.
     *                                  (MI .GT. MJ))) THEN
!
                                          NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                                  + NT1AM(ISYMAI)*(NBJ-1) + NAI
!
                                       IF (ISYMAJ .EQ. ISYMBI) THEN
                                          NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *                                          + INDEX(NAJ,NBI)
                                       ELSEIF (ISYMAJ .LT. ISYMBI) THEN
                                          NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *                                    + NT1AM(ISYMAJ)*(NBI-1)+NAJ
                                       ELSEIF (ISYMAJ .GT. ISYMBI) THEN
                                          NAJBI = IT2AM(ISYMBI,ISYMAJ)
     *                                    + NT1AM(ISYMBI)*(NAJ-1)+NBI
                                       ENDIF
!
                                       C2AM(NAJBI) = - WORK(NAIBJ)
                                       C2AM(NAIBJ) =   WORK(NAIBJ)
!
                                   ENDIF
!
                                 ENDIF
!
                           ENDIF
  160                      CONTINUE
  150                   CONTINUE
                     ENDIF
  140                CONTINUE
  130             CONTINUE
  120          CONTINUE
  110    CONTINUE
       ENDIF
  100 CONTINUE
!
      CALL QEXIT('CCRHS3_R2IJ')
!
      RETURN
      END
      SUBROUTINE CCRHS_A3(OMEGA2,T2AM,GAMMA,WORK,LWORK,ISYGAM,ISYVEC,
     *                    IOPT,ANTISYM)
C
C     Written by Henrik Koch & Ove Christiansen 20-Jan-1994
C
C     Generalised to non. total sym gamma (isygam) og non. tot. sym
C     double excitation vector (isyvec) Ove Christiansen 29-7-1995
C
C     Generalised to handle left hand side contribution (IOPT 2) as
C     well as usual contributions (IOPT 1) by Asger Halkier 22/11-95.
C
C     Introduced the ANTISYM logical to calculate either the 
C     symmetric or the antisymmetric square up of GAMMA
C
C     Purpose: Calculate A-term.
C
      IMPLICIT NONE
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
!
      INTEGER LWORK, ISYGAM, ISYVEC, IOPT
      INTEGER ISAIBJ, ISYMLJ, ISYMKI, KSCR1, KEND1, LWRK1
      INTEGER ISYML, ISYMJ, NLJ, ISYMK, ISYMI, NKI, NKILJ
      INTEGER NSTO, ISYMB
      INTEGER KOFF1, KOFF2, KOFF3, NVIRA, ISYMA, NBL, NAI, NAIBJ
      INTEGER NTOT, ISYMBJ, ISYMAI, ISYMBL, ISYMAK, KSCR2, KEND2
      INTEGER LWRK2, NBJ, NRHFK
      INTEGER INDEX
!
#if defined (SYS_CRAY)
      REAL OMEGA2(*), GAMMA(*), T2AM(*), WORK(LWORK)
      REAL ZERO, ONE, XMONE, FACT
#else
      DOUBLE PRECISION OMEGA2(*), GAMMA(*), T2AM(*), WORK(LWORK)
      DOUBLE PRECISION ZERO, ONE, XMONE, FACT
#endif
      PARAMETER(ZERO=0.0D00, ONE=1.0D00, XMONE = -1.0D00)
      LOGICAL ANTISYM
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      CALL QENTER('CCRHS_A3')
C
C----------------------------
C     Calculate contribution.
C----------------------------
C
      ISAIBJ = MULD2H(ISYGAM,ISYVEC)
C
      DO 100 ISYMLJ = 1,NSYM
C
         ISYMKI = MULD2H(ISYMLJ,ISYGAM)
C
         KSCR1 = 1
         KEND1 = KSCR1 + NMATIJ(ISYMKI)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient space for allocation in CCRHS_A3')
         END IF
C
         DO 110 ISYMJ = 1,NSYM
C
            ISYML = MULD2H(ISYMJ,ISYMLJ)
C
            DO 120 J = 1,NRHF(ISYMJ)
C
               DO 130 L = 1,NRHF(ISYML)
C
                  IF (IOPT .EQ. 1) THEN
C
                     NLJ = IMATIJ(ISYML,ISYMJ)
     *                   + NRHF(ISYML)*(J - 1) + L
C
                  ELSE IF (IOPT .EQ. 2) THEN
C
                     NLJ = IMATIJ(ISYMJ,ISYML)
     *                   + NRHF(ISYMJ)*(L - 1) + J
C
                  ENDIF
C
                  DO 140 ISYMK = 1,NSYM
C
                     ISYMI = MULD2H(ISYMK,ISYMKI)
C
                     DO 150 I = 1,NRHF(ISYMI)
C
                        DO 160 K = 1,NRHF(ISYMK)
C
                           IF (IOPT .EQ. 1) THEN
C
                              NKI = IMATIJ(ISYMK,ISYMI)
     *                            + NRHF(ISYMK)*(I - 1) + K
C
                           ELSE IF (IOPT .EQ. 2) THEN
C
                              NKI = IMATIJ(ISYMI,ISYMK)
     *                            + NRHF(ISYMI)*(K - 1) + I
C
                           ENDIF
C
                           IF (ISYMKI .EQ. ISYMLJ) THEN
                              NKILJ = IGAMMA(ISYMKI,ISYMLJ)
     *                              + INDEX(NKI,NLJ)
                              FACT = ONE
                              IF (NKI .EQ. NLJ) FACT = ZERO
                              IF (NKI .LT. NLJ) FACT = XMONE
                           ELSE
                              IF (ISYMKI .LT. ISYMLJ) THEN
                                 NKILJ = IGAMMA(ISYMKI,ISYMLJ)
     *                                 + NMATIJ(ISYMKI)*(NLJ - 1) + NKI
                                 FACT  = XMONE
                              ELSE
                                 NKILJ = IGAMMA(ISYMLJ,ISYMKI)
     *                                 + NMATIJ(ISYMLJ)*(NKI - 1) + NLJ
                                 FACT  = ONE
                              ENDIF
                           ENDIF
C
                           IF (.NOT. ANTISYM) FACT = ONE
!
                           NSTO = IMATIJ(ISYMK,ISYMI)
     *                          + NRHF(ISYMK)*(I - 1) + K
C
                           WORK(KSCR1 + NSTO - 1) = FACT * GAMMA(NKILJ)
C
  160                   CONTINUE
  150                CONTINUE
  140             CONTINUE
C
                  DO 170 ISYMB = 1,NSYM
C
                     ISYMBJ = MULD2H(ISYMB,ISYMJ)
                     ISYMAI = MULD2H(ISYMBJ,ISAIBJ)
                     ISYMBL = MULD2H(ISYMB,ISYML)
                     ISYMAK = MULD2H(ISYVEC,ISYMBL)
C
                     KSCR2 = KEND1
                     KEND2 = KSCR2 + NT1AM(ISYMAI)
                     LWRK2 = LWORK - KEND2
C
                     IF (LWRK2 .LT. 0) THEN
                        CALL QUIT('Insufficient space in CCRHS_A3')
                     END IF
C
                     IF (ISYMAI .GT. ISYMBJ) GOTO 170
C
                     DO 180 B = 1,NVIR(ISYMB)
C
                        NBJ = IT1AM(ISYMB,ISYMJ)
     *                      + NVIR(ISYMB)*(J - 1) + B
                        NBL = IT1AM(ISYMB,ISYML)
     *                      + NVIR(ISYMB)*(L - 1) + B
C
                        CALL DZERO(WORK(KSCR2),NT1AM(ISYMAI))
C
                        DO 190 ISYMI = 1,NSYM
C
                           ISYMK = MULD2H(ISYMI,ISYMKI)
                           ISYMA = MULD2H(ISYMK,ISYMAK)
C
                           NVIRA = MAX(NVIR(ISYMA),1)
                           NRHFK = MAX(NRHF(ISYMK),1)
C
                           KOFF1 = IT2SQ(ISYMAK,ISYMBL)
     *                           + NT1AM(ISYMAK)*(NBL - 1)
     *                           + IT1AM(ISYMA,ISYMK) + 1
                           KOFF2 = KSCR1 + IMATIJ(ISYMK,ISYMI)
                           KOFF3 = KSCR2 + IT1AM(ISYMA,ISYMI)
C
                           CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),
     *                                NRHF(ISYMK),ONE,T2AM(KOFF1),
     *                                NVIRA,WORK(KOFF2),NRHFK,ZERO,
     *                                WORK(KOFF3),NVIRA)
C
  190                   CONTINUE
C
                        IF (ISYMAI .EQ. ISYMBJ) THEN
                           NTOT = NBJ
                        ELSE
                           NTOT = NT1AM(ISYMAI)
                        ENDIF
C
                        DO 200 NAI = 1,NTOT
C
                           IF (ISYMAI .EQ. ISYMBJ) THEN
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                              + INDEX(NAI,NBJ)
                           ELSE
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                              + NT1AM(ISYMAI)*(NBJ - 1) + NAI
                           ENDIF
C
                           OMEGA2(NAIBJ) = OMEGA2(NAIBJ)
     *                                   + WORK(KSCR2 + NAI - 1)
C
  200                   CONTINUE
C
  180                CONTINUE
  170             CONTINUE
C
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('CCRHS_A3')
C
      RETURN
      END
C  /* Deck cc_t2sq3 */
      SUBROUTINE CC_T2SQ3(T2AM,T2SQ,ISYM)
!
!--------------------------------------------------------
!     Kasper Hald 8/3-1999 to squareup a 
!     antisymmetric matrix as in the triplet case.
!
!     Based on CC_T2SQ by Henrik Koch, Alfredo Sanchez
!     and Ove Christiansen.
!--------------------------------------------------------
!
      IMPLICIT NONE
#if defined (SYS_CRAY)
      REAL T2AM(*), T2SQ(*)
#else 
      DOUBLE PRECISION T2AM(*), T2SQ(*)
#endif
      INTEGER ISYM, KOFF1, KOFF2, ISYMBJ, KOFF, ISYMAI, NAMP, NAI
      INTEGER NBJ
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
!
      CALL QENTER('CC_T2SQ3')
!
      IF (ISYM.EQ.1) THEN
         KOFF1 = 1
         KOFF2 = 1
         DO 100 ISYMBJ = 1,NSYM
            IF (NT1AM(ISYMBJ) .GT. 0) THEN
               CALL SQMATR3(NT1AM(ISYMBJ),T2AM(KOFF1),T2SQ(KOFF2))
               KOFF1 = KOFF1 + NT1AM(ISYMBJ)*(NT1AM(ISYMBJ)+1)/2
               KOFF2 = KOFF2 + NT1AM(ISYMBJ)*NT1AM(ISYMBJ)
            ENDIF
  100    CONTINUE
!
      ELSE
!
         KOFF = 1
         DO 200 ISYMBJ = 1,NSYM
            ISYMAI = MULD2H(ISYM,ISYMBJ)
!
            IF (ISYMBJ.GT.ISYMAI) THEN
!
               NAMP = NT1AM(ISYMAI)*NT1AM(ISYMBJ)
!
               IF (NAMP .GT. 0) THEN
                  KOFF1 = IT2SQ(ISYMAI,ISYMBJ) + 1
                  CALL DCOPY(NAMP,T2AM(KOFF),1,T2SQ(KOFF1),1)
                  NAI = MAX(NT1AM(ISYMAI),1)
                  NBJ = MAX(NT1AM(ISYMBJ),1)
                  KOFF2 = IT2SQ(ISYMBJ,ISYMAI) + 1
                  CALL TRM3(T2AM(KOFF),NAI,NT1AM(ISYMAI),NT1AM(ISYMBJ),
     *                        T2SQ(KOFF2),NBJ)
                  KOFF = KOFF + NAMP
!
               ENDIF
!
            ENDIF
!
  200    CONTINUE
!
      ENDIF
!
      CALL QEXIT('CC_T2SQ3')
!
      RETURN
      END
!  /* Deck trm3 */
      SUBROUTINE TRM3(A,LDA,M,N,B,LDB)
!
!---------------------------------------------------------------
!
!     Transpose a matrix A with dimension m,n 
!     in array with logical dim. lda.
!     and put minus the result into B with logical dim. ldb.
!
!     Kasper Hald 8/3 - 1999
!
!     Based on TRM by Ove Christiansen.
!---------------------------------------------------------------
!
      IMPLICIT NONE
#include <priunit.h>
!
      INTEGER LDA, LDB, M, N, I
#if defined (SYS_CRAY)
      REAL A(LDA,*), B(LDB,*)
      REAL XMONE
#else
      DOUBLE PRECISION A(LDA,*), B(LDB,*)
      DOUBLE PRECISION XMONE
#endif
      PARAMETER(XMONE = -1.0D00)
!
      CALL QENTER('TRM3')
!
      DO 100 I = 1, N
!
         CALL DSCAL(M,XMONE,A(1,I),1)
         CALL DCOPY(M,A(1,I),1,B(I,1),LDB)
         CALL DSCAL(M,XMONE,A(1,I),1)
!
 100  CONTINUE
!
      CALL QEXIT('TRM3')
!
      RETURN
      END
C  /* Deck sqmatr3 */
      SUBROUTINE SQMATR3(NDIM,PKMAT,SQMAT)
!
!-----------------------------------------------------
!     Written by Kasper Hald 8/3-1999 
!
!     This subroutine squares up the packed 
!     triplet matrix for the totalsymmetric case.
!
!     Based on SQMATR by Henrik Koch.
!-----------------------------------------------------
!
      IMPLICIT NONE
#include <priunit.h>
      INTEGER I, J, NDIM, IJ
#if defined (SYS_CRAY)
      REAL PKMAT(*), SQMAT(NDIM,NDIM)
      REAL ZERO,XMONE
#else
      DOUBLE PRECISION PKMAT(*), SQMAT(NDIM,NDIM)
      DOUBLE PRECISION ZERO,XMONE
#endif
!
      PARAMETER(XMONE = -1.0D00)
! 
      CALL QENTER('SQMATR3')
!
      DO 100 I = 1,NDIM
         DO 110 J = 1,I
!
               IJ = I*(I-1)/2 + J
               SQMAT(I,J) = XMONE * PKMAT(IJ)
               SQMAT(J,I) = PKMAT(IJ)
!
  110    CONTINUE
  100 CONTINUE
!
      CALL QEXIT('SQMATR3')
!
      RETURN
      END
C  /* Deck cc_t2motrip */
      SUBROUTINE CC_T2MOTRIP(RHO1,CTR2,ISYMC2,OMEGA2,RHO2,GAMMA,
     *                       XLAMDP,XLAMPC,ISYMPC,WORK,LWORK,ISYMBF,
     *                       ICON,RHO22,RHO22CONT,ANTISYM)
C
C     Henrik Koch and Alfredo Sanchez.       15-July-1994
C
C     Transform the Omega2 vector from the AO basis to the MO
C     basis.
C
C     Ove Christiansen 4-8-1995:
C
C     Generalizations for CC response.
C
C        1.ISYMBF is the symmetry of the BF (ali,bej) vector.
C        2.Transform with a non total symmetric lambda matrix.
C          (one with sym 1 and one with sym isympc)
C
C        note that if newgam is true gamma is the gamma vector on return
C        with the same symmetry as the input BF. 
C        (transformed with xlamdp)
C
C        if newgam is false the gamma intermediate is not returned.
C
C        ICON is 2 for response to calculat a-tild,ibj and ai,b-tilde,j
C
C        NB these changes are only carried through completely and
C        tested for omegor
C
C     Asger Halkier 2/11-1995:
C
C        For ICON equal to 3 the contraction of the (ali,bej) vector with
C        the trialvector CTR2 (i.e the LT21BF-term) is calculated and
C        stored in RHO1!
C
C     Ove Christiansen 4-10-1996:
C 
C        For use in F-matrix generalize ICON .EQ. 3 section
C
!     Kasper Hald and Christof Haettig. 12-3-1999
!
!        If ANTISYM then rho is calculated as 
!        INTP*KT2MM + INTM*KT2MP
!
!        For ICON .EQ. 1 AND antisym then we will get 
!        Lambda(al a) * Lambda(be b) * rho(ANTISYM)
!
!        To ONLY calculate the new GAMMA ICON=4
!        The calculated gamma will be added to the excisting gamma 
!
!        For ICON .EQ. 5 we calculate 
!       (Lambda(bar)(be b)Lambda(al a) - Lambda(bar)(al a)Lambda(be b))
!       * rho(symmetric) and store it in RHO2.
!       (Lambda(bar)(be b)Lambda(al a) + Lambda(bar)(al a)Lambda(be b))
!       * rho(symmetric) and store it in RHO22.
!
!        For ICON .EQ. 6 : Here KT2MP is identical to zero (C2+ in the
!        triplet case) so we only calculate INTM*KT2MM
!        
C     NOTE: Linear response options only valid and debugged for OMEGOR!
C
      IMPLICIT NONE
#include <priunit.h>
#include <maxorb.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <symsq.h>
#include <cclr.h>
!
      INTEGER INDEX, ISYMBF, ISYMPC, ISYMO1, ISYMO2, ISYMC2, ICON
      INTEGER ISYMJ, ISYMI, ISYMIJ, ISALBE, ISYMAB, ISYBE, ISYAL
      INTEGER ISYALI, ISYBEJ, ISYMA, NVA, NRA, ISYMB, NVB, NRB
      INTEGER KSCR1, KSCR2, KSCR3, KSCR4, KSCR5, KEND1, LWRK1
      INTEGER LWORK, NAI, NBJ, NAB, NAIBJ, NBJAI, NIJ, NABP
      INTEGER NABIJP, NABIJM, ISYMAI, ISYMBJ, NBASA, NBASB
      INTEGER NVIRA, KOFF1, KOFF2, ISYMK, ISYMC, ISYMD, ISYDI
      INTEGER ISYCJ, LENGTH, NTOTAL, NTOTBE, NTOTK, NCJ
      INTEGER NDICJ, NCK, NRHFA1, ISYML, ISYMKI, ISYMLJ, NLJ
      INTEGER NKL, NKI, NKILJ
!
#if defined (SYS_CRAY)
      REAL ZERO, HALF, ONE, TWO, FAC, FAC1, FAC2, FACT
      REAL RHO1(*), CTR2(*), OMEGA2(*), RHO2(*), GAMMA(*)
      REAL GAMMA(*), XLAMDP(*), WORK(*), XLAMPC(*), RHO22(*)
#else
      DOUBLE PRECISION ZERO, HALF, ONE, TWO, FAC, FAC1, FAC2, FACT
      DOUBLE PRECISION RHO1(*), CTR2(*), OMEGA2(*), RHO2(*), GAMMA(*)
      DOUBLE PRECISION XLAMDP(*), WORK(*), XLAMPC(*), RHO22(*)
#endif
      PARAMETER (ZERO= 0.0D00, HALF= 0.5D00, ONE= 1.0D00, TWO= 2.0D00)
!
      LOGICAL ANTISYM, RHO22CONT
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      CALL QENTER('CC_T2MOTRIP')
C
      ISYMO2 = MULD2H(ISYMBF,ISYMPC)
      ISYMO1 = MULD2H(ISYMO2,ISYMC2)
C
      IF (ICON .NE. 3) THEN
         CALL DZERO(RHO2,NT2AM(ISYMO2))
      ENDIF
C
      DO 100 ISYMJ = 1,NSYM
         DO 110 ISYMI = 1,NSYM
C
            ISYMIJ = MULD2H(ISYMI,ISYMJ)
            ISALBE = MULD2H(ISYMIJ,ISYMBF)
            ISYMAB = MULD2H(ISYMIJ,ISYMO2)
C
            DO 120 ISYBE = 1,NSYM
C
               ISYAL  = MULD2H(ISYBE,ISALBE)
               ISYALI = MULD2H(ISYAL,ISYMI)
               ISYBEJ = MULD2H(ISYBE,ISYMJ)
C
C-----------------------------------------------
C              Dynamic allocation of work space.
C-----------------------------------------------
C
               ISYMA = MULD2H(ISYAL,ISYMPC)
               NVA = MAX(NVIR(ISYMA),NVIR(ISYAL))
               NRA = MAX(NRHF(ISYMA),NRHF(ISYAL))
               ISYMB = MULD2H(ISYBE,ISYMPC)
               NVB = MAX(NVIR(ISYMB),NVIR(ISYBE),NRHF(ISYBE))
               NRB = MAX(NRHF(ISYMB),NRHF(ISYBE))
C
               KSCR1 = 1
               KSCR2 = KSCR1 + NBAS(ISYAL)*NBAS(ISYBE)
               KSCR3 = KSCR2 + NBAS(ISYAL)*NVB
               IF (NEWGAM) THEN
                  KSCR4 = KSCR3 + NVA*NVB
                  KSCR5 = KSCR4 + NBAS(ISYAL)*NRB
                  KEND1 = KSCR5 + NRA*NRB
               ELSE
                  KEND1 = KSCR3 + NVA*NVB
               END IF
               LWRK1 = LWORK - KEND1
C
               IF (LWRK1 .LT. 0) THEN
                  CALL QUIT('Not enough space in CC_T2MOTRIP')
               END IF
C
               DO 130 J = 1,NRHF(ISYMJ)
                  DO 140 I = 1,NRHF(ISYMI)
C
C------------------------------------------
C                    Squareup the AB block.
C------------------------------------------
C
                     IF ((.NOT. OMEGSQ) .AND. (.NOT. OMEGOR)) THEN
C
                     DO 150 B = 1,NBAS(ISYBE)
                        NBJ   = IT1AO(ISYBE,ISYMJ)
     *                        + NBAS(ISYBE)*(J-1) + B
                        DO 155 A = 1,NBAS(ISYAL)
C
                           NAI   = IT1AO(ISYAL,ISYMI)
     *                           + NBAS(ISYAL)*(I-1) + A
                           NAB   = KSCR1 + NBAS(ISYAL)*(B - 1) + A - 1
C
                           IF (ISYMO2 .EQ. 1) THEN
                              NAIBJ = IT2AO(ISYALI,ISYBEJ)
     *                              + INDEX(NAI,NBJ)
                           ELSEIF (ISYALI .LT. ISYBEJ) THEN
                              NAIBJ = IT2AO(ISYALI,ISYBEJ)
     *                              + NT1AO(ISYALI)*(NBJ - 1) + NAI
                           ELSEIF (ISYALI .GT. ISYBEJ) THEN
                              NAIBJ = IT2AO(ISYALI,ISYBEJ)
     *                              + NT1AO(ISYBEJ)*(NAI - 1) + NBJ
                           ENDIF
C
                           WORK(NAB) = OMEGA2(NAIBJ)
C
  155                   CONTINUE
  150                CONTINUE
C
                     ENDIF
C
                     IF (OMEGSQ) THEN
C
                     DO 160 B = 1,NBAS(ISYBE)
                        NBJ   = IT1AO(ISYBE,ISYMJ)
     *                        + NBAS(ISYBE)*(J-1) + B
                        DO 165 A = 1,NBAS(ISYAL)
C
                           NAI   = IT1AO(ISYAL,ISYMI)
     *                           + NBAS(ISYAL)*(I-1) + A
                           NAB   = KSCR1 + NBAS(ISYAL)*(B - 1) + A - 1
C
                           NAIBJ = IT2AOS(ISYALI,ISYBEJ)
     *                           + NT1AO(ISYALI)*(NBJ - 1) + NAI
                           NBJAI = IT2AOS(ISYBEJ,ISYALI)
     *                           + NT1AO(ISYBEJ)*(NAI - 1) + NBJ
C
                           WORK(NAB) = OMEGA2(NAIBJ) + OMEGA2(NBJAI)
C
  165                   CONTINUE
  160                CONTINUE
C
                     ENDIF
C
                     IF (OMEGOR) THEN
!
                     IF (.NOT. ANTISYM) THEN
!
                     IF (ISYMI .EQ. ISYMJ) THEN
                        NIJ = IMIJP(ISYMI,ISYMJ) + INDEX(I,J)
                        FAC1 = ONE
                        IF (I .GT. J) FAC1 = -ONE
                     ELSE IF (ISYMI .LT. ISYMJ) THEN
                        NIJ = IMIJP(ISYMI,ISYMJ)
     *                      + NRHF(ISYMI)*(J - 1) + I
                        FAC1 = ONE
                     ELSE
                        NIJ = IMIJP(ISYMI,ISYMJ)
     *                      + NRHF(ISYMJ)*(I - 1) + J
                        FAC1 = -ONE
                     ENDIF
C
                     DO 166 B = 1,NBAS(ISYBE)
                        DO 167 A = 1,NBAS(ISYAL)
C
                           IF (ISYAL .EQ. ISYBE) THEN
                              NABP = IAODPK(ISYAL,ISYBE)
     *                             + INDEX(A,B)
                              FAC2 = ONE
                              IF (A .GT. B) FAC2 = -ONE
                           ELSE IF (ISYAL .LT. ISYBE) THEN
                              NABP = IAODPK(ISYAL,ISYBE)
     *                             + NBAS(ISYAL)*(B - 1) + A
                              FAC2 = ONE
                           ELSE
                              NABP = IAODPK(ISYAL,ISYBE)
     *                             + NBAS(ISYBE)*(A - 1) + B
                              FAC2 = -ONE
                           ENDIF
C
                           IF (ICON .NE. 6) THEN
                           NABIJP = IT2ORT(ISALBE,ISYMIJ)
     *                            + NNBST(ISALBE)*(NIJ - 1) + NABP
C
                           NABIJM = NT2ORT(ISYMBF)
     *                            + IT2ORT(ISALBE,ISYMIJ)
     *                            + NNBST(ISALBE)*(NIJ - 1) + NABP
!
                           ELSE
                           NABIJM = IT2ORT(ISALBE,ISYMIJ)
     *                            + NNBST(ISALBE)*(NIJ - 1) + NABP
!
                           ENDIF
C
                           NAB   = KSCR1 + NBAS(ISYAL)*(B - 1) + A - 1
C
                           FAC = FAC1*FAC2
C
                           IF (ICON .NE. 6) THEN
                           WORK(NAB) =
     *                       HALF*(OMEGA2(NABIJP) + FAC*OMEGA2(NABIJM))
                           ELSE
                           WORK(NAB) = HALF*FAC*OMEGA2(NABIJM)
                           ENDIF
C
  167                   CONTINUE
  166                CONTINUE
C
                     ELSE
!
                        IF (ISYMI .EQ. ISYMJ) THEN
                           NIJ = IMIJP(ISYMI,ISYMJ) + INDEX(I,J)
                           FAC1 = ONE
                           IF (I .GT. J) FAC1 = -ONE
!
                        ELSE IF (ISYMI .LT. ISYMJ) THEN
                           NIJ = IMIJP(ISYMI,ISYMJ)
     *                         + NRHF(ISYMI)*(J - 1) + I
                           FAC1 = ONE
                        ELSE
                           NIJ = IMIJP(ISYMI,ISYMJ)
     *                         + NRHF(ISYMJ)*(I - 1) + J
                           FAC1 = -ONE
                        ENDIF
!
                     DO 168 B = 1,NBAS(ISYBE)
                        DO 169 A = 1,NBAS(ISYAL)
!
                           IF (ISYAL .EQ. ISYBE) THEN
                              NABP = IAODPK(ISYAL,ISYBE)
     *                             + INDEX(A,B)
                              FAC2 = ONE
                              IF (A .GT. B) FAC2 = -ONE
                           ELSE IF (ISYAL .LT. ISYBE) THEN
                              NABP = IAODPK(ISYAL,ISYBE)
     *                             + NBAS(ISYAL)*(B - 1) + A
                              FAC2 = ONE
                           ELSE
                              NABP = IAODPK(ISYAL,ISYBE)
     *                             + NBAS(ISYBE)*(A - 1) + B
                              FAC2 = -ONE
                           ENDIF
!
                           NABIJP = IT2ORT(ISALBE,ISYMIJ)
     *                            + NNBST(ISALBE)*(NIJ - 1) + NABP
!
                           NABIJM = NT2ORT(ISYMBF)
     *                            + IT2ORT(ISALBE,ISYMIJ)
     *                            + NNBST(ISALBE)*(NIJ - 1) + NABP
!
                           NAB   = KSCR1 + NBAS(ISYAL)*(B - 1) + A - 1
!
!
                           IF ((ISYAL .EQ. ISYBE) .AND. 
     *                         (ISYMI. EQ. ISYMJ) .AND. (A .EQ. B)
     *                         .AND. (I .EQ. J)) THEN
!
                               WORK(NAB) = ZERO
!
                           ELSE
!
                              WORK(NAB) =
     *                          HALF*(OMEGA2(NABIJP)*FAC2
     *                              + FAC1*OMEGA2(NABIJM))
                           ENDIF
!
  169                   CONTINUE
  168                CONTINUE
!
                     ENDIF
C
C------------------------------------------------------------
C                    Transform the AB block to virtual space.
C------------------------------------------------------------
C
!                 IF (.NOT. (ICON .EQ. 4)) THEN 
!                  
                  IF (ICON .NE. 3) THEN
C
                     ISYMA = MULD2H(ISYAL,ISYMPC)
                     ISYMB = ISYBE
                     ISYMAI = MULD2H(ISYMA,ISYMI)
                     ISYMBJ = MULD2H(ISYMB,ISYMJ)
C
                     NBASA = MAX(NBAS(ISYAL),1)
                     NBASB = MAX(NBAS(ISYBE),1)
                     NVIRA = MAX(NVIR(ISYMA),1)
C
                     KOFF1 = ILMVIR(ISYBE) + 1
C
                     CALL DGEMM('N','N',NBAS(ISYAL),NVIR(ISYMB),
     *                          NBAS(ISYBE),ONE,WORK(KSCR1),NBASA,
     *                          XLAMDP(KOFF1),NBASB,ZERO,WORK(KSCR2),
     *                          NBASA)
C
                     KOFF2 = IGLMVI(ISYAL,ISYMA) + 1
C
                     CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMB),
     *                          NBAS(ISYAL),ONE,XLAMPC(KOFF2),NBASA,
     *                          WORK(KSCR2),NBASA,ZERO,WORK(KSCR3),
     *                          NVIRA)
C
C--------------------------------------------
C                    Store the omega2 vector.
C--------------------------------------------
C
                     DO 170 B = 1,NVIR(ISYMB)
                        NBJ   = IT1AM(ISYMB,ISYMJ)
     *                        + NVIR(ISYMB)*(J-1) + B
                        DO 180 A = 1,NVIR(ISYMA)
C
                           NAI   = IT1AM(ISYMA,ISYMI)
     *                           + NVIR(ISYMA)*(I-1) + A
                           NAB   = KSCR3 + NVIR(ISYMA)*(B - 1) + A - 1
C
                           IF (ISYMAI .EQ. ISYMBJ) THEN
C
                              IF (NAI .GT. NBJ) GOTO 180
C
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                              + INDEX(NAI,NBJ)
                           ELSEIF (ISYMAI .LT. ISYMBJ) THEN
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                              + NT1AM(ISYMAI)*(NBJ - 1) + NAI
                           ELSEIF (ISYMAI .GT. ISYMBJ) THEN
                              GOTO 180
c                             NAIBJ = IT2AM(ISYMAI,ISYMBJ)
c    *                              + NT1AM(ISYMBJ)*(NAI - 1) + NBJ
                           ENDIF
C
                           RHO2(NAIBJ) = RHO2(NAIBJ) + WORK(NAB)
!
                           IF ((ICON .EQ.5) .AND. (RHO22CONT)) THEN
                              RHO22(NAIBJ) = RHO22(NAIBJ) + WORK(NAB)
                           ENDIF
C
  180                   CONTINUE
  170                CONTINUE
C
                  ENDIF
C
C--------------------------------------
C                    CCLR contribution.
C--------------------------------------
C
                     IF ((ICON .EQ. 2) .OR. (ICON .EQ. 5)) THEN
C
                        CALL DZERO(WORK(KSCR2),NVA*NVB)
                        ISYMA = ISYAL
                        ISYMB = MULD2H(ISYBE,ISYMPC)
                        ISYMAI = MULD2H(ISYMA,ISYMI)
                        ISYMBJ = MULD2H(ISYMB,ISYMJ)
C
                        NBASA = MAX(NBAS(ISYAL),1)
                        NBASB = MAX(NBAS(ISYBE),1)
                        NVIRA = MAX(NVIR(ISYMA),1)
C
                        IF ((ICON .EQ. 5)) THEN
                            FACT = -ONE
                        ELSE
                            FACT = ONE
                        ENDIF
C
                        KOFF1 = IGLMVI(ISYBE,ISYMB) + 1
C
                        CALL DGEMM('N','N',NBAS(ISYAL),NVIR(ISYMB),
     *                             NBAS(ISYBE),FACT,WORK(KSCR1),NBASA,
     *                             XLAMPC(KOFF1),NBASB,ZERO,WORK(KSCR2),
     *                             NBASA)
C
                        KOFF2 = ILMVIR(ISYAL) + 1
C
                        CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMB),
     *                             NBAS(ISYAL),ONE,XLAMDP(KOFF2),NBASA,
     *                             WORK(KSCR2),NBASA,ZERO,WORK(KSCR3),
     *                             NVIRA)
C
C--------------------------------------------
C                    Store the omega2 vector.
C--------------------------------------------
C
                     DO 181 B = 1,NVIR(ISYMB)
                        NBJ   = IT1AM(ISYMB,ISYMJ)
     *                        + NVIR(ISYMB)*(J-1) + B
                        DO 182 A = 1,NVIR(ISYMA)
C
                           NAI   = IT1AM(ISYMA,ISYMI)
     *                           + NVIR(ISYMA)*(I-1) + A
C
                           IF (ISYMAI .EQ. ISYMBJ) THEN
                              IF (NAI .GT. NBJ ) GOTO 182
                              NAIBJ = IT2AM(ISYALI,ISYBEJ)
     *                              + INDEX(NAI,NBJ)
                           ELSEIF (ISYMAI .LT. ISYMBJ) THEN
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                              + NT1AM(ISYMAI)*(NBJ - 1) + NAI
                           ELSEIF (ISYMAI .GT. ISYMBJ) THEN
                              GOTO 182
c                             NAIBJ = IT2AM(ISYMAI,ISYMBJ)
c    *                              + NT1AM(ISYMBJ)*(NAI - 1) + NBJ
                           ENDIF
C
                           NAB  = KSCR3+ NVIR(ISYMA)*(B - 1) + A - 1
                           RHO2(NAIBJ) = RHO2(NAIBJ) + WORK(NAB)
!
                           IF ((ICON .EQ. 5) .AND. (RHO22CONT)) THEN
                             RHO22(NAIBJ) = RHO22(NAIBJ) - WORK(NAB)
                           ENDIF
C
  182                   CONTINUE
  181                CONTINUE
C
                     ENDIF
C
C============================================================
C                    Section for calculating the LT21BF-term.
C============================================================
C
                     IF (ICON .EQ. 3) THEN
C
                        ISYMK = ISYBE
                        ISYMD = MULD2H(ISYAL,ISYMPC)
                        ISYMC = MULD2H(ISYMK,ISYMO1)
                        ISYDI = MULD2H(ISYMD,ISYMI)
                        ISYCJ = MULD2H(ISYMC,ISYMJ)
C
                        LENGTH = NBAS(ISYAL)*NRHF(ISYMK)
C
                        CALL DZERO(WORK(KSCR2),LENGTH)
C
C----------------------------------------------------------
C                       Transform the AO-block to MO-basis.
C----------------------------------------------------------
C
                        KOFF1  = ILMRHF(ISYMK) + 1
C
                        NTOTAL = MAX(NBAS(ISYAL),1)
                        NTOTBE = MAX(NBAS(ISYBE),1)
C
                        CALL DGEMM('N','N',NBAS(ISYAL),NRHF(ISYMK),
     *                             NBAS(ISYBE),ONE,WORK(KSCR1),NTOTAL,
     *                             XLAMDP(KOFF1),NTOTBE,ZERO,
     *                             WORK(KSCR2),NTOTAL)
C
                        KOFF2  = IGLMVI(ISYAL,ISYMD) + 1
C
                        NTOTAL = MAX(NBAS(ISYAL),1)
                        NTOTK  = MAX(NRHF(ISYMK),1)
C
                        CALL DGEMM('T','N',NRHF(ISYMK),NVIR(ISYMD),
     *                             NBAS(ISYAL),ONE,WORK(KSCR2),NTOTAL,
     *                             XLAMPC(KOFF2),NTOTAL,ZERO,
     *                             WORK(KSCR3),NTOTK)
C
C-----------------------------------------------------------------
C                       Contraction with CTR2 & storage in result.
C-----------------------------------------------------------------
C
                        DO 47 C = 1,NVIR(ISYMC)
C
                           NCJ   = IT1AM(ISYMC,ISYMJ)
     *                           + NVIR(ISYMC)*(J - 1) + C
                           NDICJ = IT2SQ(ISYDI,ISYCJ)
     *                           + NT1AM(ISYDI)*(NCJ - 1)
     *                           + IT1AM(ISYMD,ISYMI)
     *                           + NVIR(ISYMD)*(I - 1) + 1
                           NCK   = IT1AM(ISYMC,ISYMK) + C
C
                           CALL DGEMV('N',NRHF(ISYMK),NVIR(ISYMD),
     *                                -ONE,WORK(KSCR3),NTOTK,
     *                                CTR2(NDICJ),1,ONE,RHO1(NCK),
     *                                NVIR(ISYMC))
C
  47                    CONTINUE
C
                     ENDIF
C
                  ENDIF 
!
C-------------------------------------------------------------
C                    Transform the AB block to occupied space.
C-------------------------------------------------------------
C
                     IF (.NOT. NEWGAM) GOTO 999
C
                     NBASA = MAX(NBAS(ISYAL),1)
                     NBASB = MAX(NBAS(ISYBE),1)
                     NRHFA1 = MAX(NRHF(ISYAL),1)
C
                     KOFF1 = ILMRHF(ISYBE) + 1
C
                     CALL DGEMM('N','N',NBAS(ISYAL),NRHF(ISYBE),
     *                          NBAS(ISYBE),ONE,WORK(KSCR1),NBASA,
     *                          XLAMDP(KOFF1),NBASB,ZERO,WORK(KSCR4),
     *                          NBASA)
C
                     KOFF2 = ILMRHF(ISYAL) + 1
C
                     CALL DGEMM('T','N',NRHF(ISYAL),NRHF(ISYBE),
     *                          NBAS(ISYAL),ONE,XLAMDP(KOFF2),NBASA,
     *                          WORK(KSCR4),NBASA,ZERO,WORK(KSCR5),
     *                          NRHFA1)
C
C-------------------------------------------
C                    Store the gamma matrix.
C-------------------------------------------
C
                     ISYMK = ISYAL
                     ISYML = ISYBE
C
                     ISYMKI = MULD2H(ISYMK,ISYMI)
                     ISYMLJ = MULD2H(ISYML,ISYMJ)
C
                     DO 190 L = 1,NRHF(ISYML)
C
                        NLJ = IMATIJ(ISYML,ISYMJ)
     *                      + NRHF(ISYML)*(J - 1) + L
C
                        DO 200 K = 1,NRHF(ISYMK)
C
                           NKL = KSCR5 + NRHF(ISYMK)*(L - 1) + K - 1
C
                           NKI = IMATIJ(ISYMK,ISYMI)
     *                         + NRHF(ISYMK)*(I - 1) + K
C
                           IF (ISYMKI .EQ. ISYMLJ) THEN
                              NKILJ = IGAMMA(ISYMKI,ISYMLJ)
     *                              + INDEX(NKI,NLJ)
                              IF (NKI .LE. NLJ) THEN
                                 GAMMA(NKILJ) = GAMMA(NKILJ)
     *                                        + WORK(NKL)
                              ENDIF
                           ELSE IF (ISYMKI .LT. ISYMLJ) THEN
                              NKILJ = IGAMMA(ISYMKI,ISYMLJ)
     *                              + NMATIJ(ISYMKI)*(NLJ - 1) + NKI
                              GAMMA(NKILJ) = GAMMA(NKILJ) + WORK(NKL)
                           ENDIF
C
  200                   CONTINUE
  190                CONTINUE
  999                CONTINUE
  140             CONTINUE
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('CC_T2MOTRIP')
C
      RETURN
      END
C  /* Deck hescompa */
      SUBROUTINE HESCOMPA(REDHES1,REDHES2,NREDH,NCOMPO,COMTHRES)
!
!     Written by Kasper Hald 3/2-2000
!
!     Compares 2 Hessian matrices.
!     This routine is only meant to be used for
!     comparisons of "small" arrays, since it scales
!     as N^8
!
!     REDHES1 and REDHES2 are the 2 Hessians (Surprise?)
!     NREDH are the number of vectors/components that are 
!     important.
!     NCOMPO are the no. of components per vector (greater or equal
!     to NREDH)
!     COMTHRES is the threshold that you compare against
!
      IMPLICIT NONE
!
#include <priunit.h>
!
      INTEGER NREDH, NCOMPO, I, J, K, L, KOFF1, KOFF2
!
#if defined (SYS_CRAY)
      REAL REDHES1(*), REDHES2(*), COMTHRES
      REAL DIFF1, DIFF2
#else
      DOUBLE PRECISION REDHES1(*), REDHES2(*), COMTHRES
      DOUBLE PRECISION DIFF1, DIFF2
#endif
!
      CALL QENTER('HESCOMPA')
!
      WRITE(LUPRI,*)'                         '
      WRITE(LUPRI,*)'THRESHOLD FOR HESCOMPA = ',COMTHRES
      WRITE(LUPRI,*)'                         '
!
      DO I=1,NREDH
!
        DO J=1,NREDH
!
          KOFF1 = (I-1)*NCOMPO + J
          KOFF2 = (J-1)*NCOMPO + I
!
C         IF ((ABS(REDHES1(KOFF1)) .LT. COMTHRES) .AND.
C    *        (ABS(REDHES2(KOFF1)) .GT. COMTHRES)) THEN
C             WRITE(LUPRI,*)'Diff.1 for (J,I) = ',J,I
C             WRITE(LUPRI,*)'REDHES1(KOFF1) = ',REDHES1(KOFF1)
C             WRITE(LUPRI,*)'REDHES2(KOFF1) = ',REDHES2(KOFF1)
C         ENDIF
          IF ((ABS(REDHES1(KOFF1)) .GT. COMTHRES) .AND.
     *        (ABS(REDHES2(KOFF1)) .LT. COMTHRES)) THEN
              WRITE(LUPRI,*)'Diff.2 for (J,I) =',J,I
              WRITE(LUPRI,*)'REDHES1(KOFF1) = ',REDHES1(KOFF1)
              WRITE(LUPRI,*)'REDHES2(KOFF1) = ',REDHES2(KOFF1)
              WRITE(LUPRI,*)'REDHES1(KOFF2) = ',REDHES1(KOFF2)
              WRITE(LUPRI,*)'REDHES2(KOFF2) = ',REDHES2(KOFF2)
          ENDIF
!
        ENDDO
!
      ENDDO
!
!     DO I=1,NREDH
!
!       DO J=1,NCOMPO
!
!         DO K=1, NREDH
!
!           DO L=1, NCOMPO
!
!              KOFF1 = (I-1)*NCOMPO + J
!              KOFF2 = (K-1)*NCOMPO + L
!              DIFF1 = REDHES1(KOFF1) - REDHES1(KOFF2)
!              DIFF2 = REDHES2(KOFF1) - REDHES2(KOFF2)
!
!              IF ((DIFF1 .LT. COMTHRES) .AND. 
!    *             (DIFF2 .GT. COMTHRES)) THEN
!                 WRITE(LUPRI,*)'COMPARISON ERROR FOR ELEMENTS NO. (',
!    *                      I,',',J,') and NO. (',K,',',L,') '
!                 WRITE(LUPRI,*)'Difference for Array1 = ',DIFF1
!                 WRITE(LUPRI,*)'Difference for Array2 = ',DIFF2
!              ENDIF
!              IF ((DIFF1 .GT. COMTHRES) .AND. 
!    *             (DIFF2 .LT. COMTHRES)) THEN
!                 WRITE(LUPRI,*)'COMPARISON ERROR FOR ELEMENTS NO. (',
!    *                      I,',',J,') and NO. (',K,',',L,') '
!                 WRITE(LUPRI,*)'Difference for Array1 = ',DIFF1
!                 WRITE(LUPRI,*)'Difference for Array2 = ',DIFF2
!              ENDIF
!
!           ENDDO
!
!         ENDDO
!
!       ENDDO
!
!     ENDDO
!
      CALL QEXIT('HESCOMPA')
!
      RETURN
      END
