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
C  /* Deck ccpt_etars_1e */
      SUBROUTINE CCPT_ETARS_1E(ETAIJ,ETAAB,
     *                        XINTIJ,XINTAI,XINTIA,XINTAB,
     *                        DIA,WORK,LWORK,ISYM)
C
C     Written by S. Coriani 21/1-2002
C
C     Version: 1.0
C
C     Purpose: To set up the one-electron contribution to the 
C              right hand side of the equation for
C              zeta-kappa-0_ij (ETAIJ) and zeta-kappa-0_ab (ETAAB)
C              from MO-integrals (XIN*) and (T) 
C              contribution to CCSD(T) density (D_ia)
C              ISYM is the symmetry of both the density and the
C              integrals!
C
C     Based on CC2_ETIJ/CC2_ETAB by A. Halkier & S. Coriani
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION ETAIJ(*), ETAAB(*)
      DIMENSION XINTIJ(*), XINTAI(*), XINTIA(*), XINTAB(*)
      DIMENSION DIA(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
C
      CALL QENTER('CCPT_ETARS_1E')
C
      DO 100 ISYMI = 1,NSYM
C
C----------------------------------------------------------------
C        Calculate terms to eta_ij.
C----------------------------------------------------------------
C
         ISYMJ  = ISYMI
         ISYMC  = MULD2H(ISYMI,ISYM)
C
         KOFFRE = IMATIJ(ISYMI,ISYMJ) + 1
C
         NTOTRE = MAX(NRHF(ISYMI),1)
         NTOTC  = MAX(NVIR(ISYMC),1)
C
         KOFF5  = IT1AM(ISYMC,ISYMI) + 1
         KOFF6  = IT1AM(ISYMC,ISYMJ) + 1
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMC),ONE,
     *              XINTIA(KOFF5),NTOTC,DIA(KOFF6),NTOTC,ONE,
     *              ETAIJ(KOFFRE),NTOTRE)
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMC),-ONE,
     *              DIA(KOFF5),NTOTC,XINTIA(KOFF6),NTOTC,ONE,
     *              ETAIJ(KOFFRE),NTOTRE)
C
  100 CONTINUE
C
      DO 101 ISYMA = 1,NSYM
C
C----------------------------------------------------------------
C        Calculate terms to eta_ab.
C----------------------------------------------------------------
C
         ISYMB  = ISYMA
         ISYMK  = MULD2H(ISYMA,ISYM)
         ISYMC  = MULD2H(ISYMA,ISYM)
C
         KOFFRE = IMATAB(ISYMA,ISYMB) + 1
C
         NTOTRE = MAX(NVIR(ISYMA),1)
         NTOTA  = MAX(NVIR(ISYMA),1)
         NTOTB  = MAX(NVIR(ISYMB),1)
C
         KOFF1  = IT1AM(ISYMA,ISYMK) + 1
         KOFF2  = IT1AM(ISYMB,ISYMK) + 1
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMK),ONE,
     *              XINTIA(KOFF1),NTOTA,DIA(KOFF2),NTOTB,ONE,
     *              ETAAB(KOFFRE),NTOTRE)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMK),-ONE,
     *              DIA(KOFF1),NTOTA,XINTIA(KOFF2),NTOTB,ONE,
     *              ETAAB(KOFFRE),NTOTRE)
C
  101 CONTINUE
C
      CALL QEXIT('CCPT_ETARS_1E')
C
      RETURN
      END
C
C------------------------------------------------------------------------
C  /* Deck ccpt_etaai_1e */
      SUBROUTINE ccpt_etaai_1e(ETAKA,XINTIJ,XINTAI,XINTIA,XINTAB,
     *                    DIA,WORK,LWORK,ISYM)
C
C     Written by Sonia Coriani 22/2 - 2002
C
C     Version: 1.0
C
C     Purpose: To set up the right hand side of the equation for
C              zeta-kappa-0 (ETAKA) from MO-integrals (XI*) and (T) 
C              one-electron density contribution D_ia to the Coupled 
C              Cluster densities 
C              ISYM is the symmetry of both the density and the 
C              integrals!
C     Based on CCDENZK0 by A. Halkier
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION ETAKA(*), XINTIJ(*), XINTAI(*), XINTIA(*), XINTAB(*)
      DIMENSION DIA(*)
      DIMENSION WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('ccpt_etaai_1e')
C
      DO 100 ISYMA = 1,NSYM
C
         ISYMI  = ISYMA
         ISYMB  = MULD2H(ISYMA,ISYM)
         ISYMJ  = MULD2H(ISYMA,ISYM)
C
         KOFFRE = IT1AM(ISYMA,ISYMI)  + 1
C
         NTOTRE = MAX(NVIR(ISYMA),1)
         NTOTB  = MAX(NVIR(ISYMB),1)
         NTOTJ  = MAX(NRHF(ISYMJ),1)
C
         KOFF1  = IMATAB(ISYMA,ISYMB) + 1
         KOFF2  = IT1AM(ISYMB,ISYMI)  + 1
C
         KOFF5  = IT1AM(ISYMA,ISYMJ)  + 1
         KOFF6  = IMATIJ(ISYMJ,ISYMI) + 1
C
C-------------------------------------------------------
C        Calculate terms originating from [H(t1),E(ia)].
C-------------------------------------------------------
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMB),
     *              ONE,XINTAB(KOFF1),NTOTRE,DIA(KOFF2),NTOTB,
     *              ONE,ETAKA(KOFFRE),NTOTRE)
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMJ),
     *              -ONE,DIA(KOFF5),NTOTRE,XINTIJ(KOFF6),NTOTJ,
     *              ONE,ETAKA(KOFFRE),NTOTRE)
C
  100 CONTINUE
C
      CALL QEXIT('ccpt_etaai_1e')
C
      RETURN
      END
C
C----------------------------------------------------------------
C  /* Deck ccpt_etars_2e */
      SUBROUTINE CCPT_ETARS_2E(ETAIJ,ETAAB,
     &                         XINTIJ,XINTAI,XINTIA,XINTAB,
     &                         DSIJ,DAI,DIA,DSAB,WORK,LWORK,ISYM)
C
C     Written by S. Coriani 11/2-2002
C
C     Version: 1.0
C
C     Purpose: To set up the two-electron contribution to the 
C              right hand side of the equation for
C              zeta-kappa-0_ij (ETAIJ) and zeta-kappa-0_ab (ETAAB)
C              from MO-integrals (XIN*) and (T) 
C              contribution to CCSD(T) densities (d_pq;gamma,delta)
C              ISYM is the symmetry of both the density and the
C              integrals!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION ETAIJ(*), ETAAB(*)
      DIMENSION XINTIJ(*), XINTAI(*), XINTIA(*), XINTAB(*)
      DIMENSION DIA(*),DAI(*),DSIJ(*),DSAB(*),WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
C
      CALL QENTER('CCPT_ETARS_2E')
C
C--------------------------------------------
C Two-electron density contribution to eta_ij
C--------------------------------------------
C
      DO 100 ISYMI = 1,NSYM
C
         ISYMJ  = ISYMI
         ISYMK  = MULD2H(ISYMI,ISYM)
         ISYMC  = MULD2H(ISYMI,ISYM)
C
         KOFFIJ = IMATIJ(ISYMI,ISYMJ) + 1
C
         NTOTI  = MAX(NRHF(ISYMI),1)
         NTOTJ  = MAX(NRHF(ISYMJ),1)
         NTOTK  = MAX(NRHF(ISYMK),1)
         NTOTC  = MAX(NVIR(ISYMC),1)        
C
C----------------------------------------------------------------
C        Calculate sum_k terms to eta_ij.
C----------------------------------------------------------------
C
         KOFF1  = IMATIJ(ISYMK,ISYMI) + 1
         KOFF2  = IMATIJ(ISYMK,ISYMJ) + 1
         KOFF3  = IMATIJ(ISYMI,ISYMK) + 1
         KOFF4  = IMATIJ(ISYMJ,ISYMK) + 1
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NRHF(ISYMK),ONE,
     *              XINTIJ(KOFF1),NTOTK,DSIJ(KOFF2),NTOTK,ONE,
     *              ETAIJ(KOFFIJ),NTOTI)
C
         CALL DGEMM('N','T',NRHF(ISYMI),NRHF(ISYMJ),NRHF(ISYMK),ONE,
     *              XINTIJ(KOFF3),NTOTI,DSIJ(KOFF4),NTOTJ,ONE,
     *              ETAIJ(KOFFIJ),NTOTI)
C
         CALL DGEMM('N','T',NRHF(ISYMI),NRHF(ISYMJ),NRHF(ISYMK),-ONE,
     *              DSIJ(KOFF3),NTOTI,XINTIJ(KOFF4),NTOTJ,ONE,
     *              ETAIJ(KOFFIJ),NTOTI)
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NRHF(ISYMK),-ONE,
     *              DSIJ(KOFF1),NTOTK,XINTIJ(KOFF2),NTOTK,ONE,
     *              ETAIJ(KOFFIJ),NTOTI)
C
C----------------------------------------------------------------
C        Calculate sum_c terms to eta_ij.
C----------------------------------------------------------------
C
         KOFF5  = IT1AM(ISYMC,ISYMI) + 1
         KOFF6  = IT1AM(ISYMC,ISYMJ) + 1
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMC),ONE,
     *              XINTAI(KOFF5),NTOTC,DAI(KOFF6),NTOTC,ONE,
     *              ETAIJ(KOFFIJ),NTOTI)
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMC),ONE,
     *              XINTIA(KOFF5),NTOTC,DIA(KOFF6),NTOTC,ONE,
     *              ETAIJ(KOFFIJ),NTOTI)
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMC),-ONE,
     *              DIA(KOFF5),NTOTC,XINTIA(KOFF6),NTOTC,ONE,
     *              ETAIJ(KOFFIJ),NTOTI)
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMC),-ONE,
     *              DAI(KOFF5),NTOTC,XINTAI(KOFF6),NTOTC,ONE,
     *              ETAIJ(KOFFIJ),NTOTI)

  100 CONTINUE
C
C--------------------------------------------
C Two-electron density contribution to eta_ab
C--------------------------------------------
C
      DO 101 ISYMA = 1,NSYM
C
C----------------------------------------------------------------
C        Calculate sum_k terms to eta_ab.
C----------------------------------------------------------------
C
         ISYMB  = ISYMA
         ISYMK  = MULD2H(ISYMA,ISYM)
         ISYMC  = MULD2H(ISYMA,ISYM)
C
         KOFFAB = IMATAB(ISYMA,ISYMB) + 1
C
         NTOTA  = MAX(NVIR(ISYMA),1)
         NTOTB  = MAX(NVIR(ISYMB),1)
         NTOTC  = MAX(NVIR(ISYMC),1)
         NTOTK  = MAX(NRHF(ISYMK),1)
C
         KOFF1  = IT1AM(ISYMA,ISYMK) + 1
         KOFF2  = IT1AM(ISYMB,ISYMK) + 1
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMK),ONE,
     *              XINTIA(KOFF1),NTOTA,DIA(KOFF2),NTOTB,ONE,
     *              ETAAB(KOFFAB),NTOTA)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMK),-ONE,
     *              DAI(KOFF1),NTOTA,XINTAI(KOFF2),NTOTB,ONE,
     *              ETAAB(KOFFAB),NTOTA)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMK),-ONE,
     *              DIA(KOFF1),NTOTA,XINTIA(KOFF2),NTOTB,ONE,
     *              ETAAB(KOFFAB),NTOTA)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMK),ONE,
     *              XINTAI(KOFF1),NTOTA,DAI(KOFF2),NTOTB,ONE,
     *              ETAAB(KOFFAB),NTOTA)
C
C----------------------------------------------------------------
C        Calculate sum_c terms to eta_ab.
C----------------------------------------------------------------
C
         KOFF3  = IMATAB(ISYMC,ISYMA) + 1
         KOFF4  = IMATAB(ISYMC,ISYMB) + 1
         KOFF5  = IMATAB(ISYMA,ISYMC) + 1
         KOFF6  = IMATAB(ISYMB,ISYMC) + 1
C
         CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMB),NVIR(ISYMC),ONE,
     *              XINTAB(KOFF3),NTOTC,DSAB(KOFF4),NTOTC,ONE,
     *              ETAAB(KOFFAB),NTOTA)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NVIR(ISYMC),-ONE,
     *              DSAB(KOFF5),NTOTA,XINTAB(KOFF6),NTOTB,ONE,
     *              ETAAB(KOFFAB),NTOTA)
C
         CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMB),NVIR(ISYMC),-ONE,
     *              DSAB(KOFF3),NTOTC,XINTAB(KOFF4),NTOTC,ONE,
     *              ETAAB(KOFFAB),NTOTA)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NVIR(ISYMC),ONE,
     *              XINTAB(KOFF5),NTOTA,DSAB(KOFF6),NTOTB,ONE,
     *              ETAAB(KOFFAB),NTOTA)
C
  101 CONTINUE
C
      CALL QEXIT('CCPT_ETARS_2E')
C
      RETURN
      END
C------------------------------------------------------------------------
C  /* Deck ccpt_etaai_2e */
      SUBROUTINE ccpt_etaai_2e(ETAKA,XINTIJ,XINTAI,XINTIA,XINTAB,
     &                         DSIJ,DAI,DIA,DSAB,WORK,LWORK,ISYM)
C
C     Written by Sonia Coriani 08/2 - 2002. Based on CCDENZK0
C
C     Version: 1.0
C
C     Purpose: To set up the right hand side of the equation for
C              zeta-kappa-0 (ETAKA) from MO-integrals (XI*) and 
C              pure (T) two-electron density contribution 
C              d_pq(gamma,delta) to the Coupled Cluster densities 
C              ISYM is the symmetry of both the density and the 
C              integrals!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION ETAKA(*), XINTIJ(*), XINTAI(*), XINTIA(*), XINTAB(*)
      DIMENSION DAI(*), DIA(*), DSIJ(*), DSAB(*)
      DIMENSION WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('ccpt_etaai_2e')
C
      DO 100 ISYMA = 1,NSYM
C
         ISYMI  = ISYMA
         ISYMB  = MULD2H(ISYMA,ISYM)
         ISYMJ  = MULD2H(ISYMA,ISYM)
C
         KOFFAI = IT1AM(ISYMA,ISYMI)  + 1
C
         NTOTA  = MAX(NVIR(ISYMA),1)
         NTOTB  = MAX(NVIR(ISYMB),1)
         NTOTJ  = MAX(NRHF(ISYMJ),1)
         NTOTI  = MAX(NRHF(ISYMI),1)
C
         KOFF1  = IMATAB(ISYMB,ISYMA) + 1
         KOFF2  = IT1AM(ISYMB,ISYMI)  + 1

         ! sum_b d_bi g_ba = sum_b (g_ba)^T d_bi

         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMB),
     *              ONE,XINTAB(KOFF1),NTOTB,DAI(KOFF2),NTOTB,ONE,
     *              ETAKA(KOFFAI),NTOTA)

         ! - sum_b ds_ab g_bi

         KOFF3  = IMATAB(ISYMA,ISYMB) + 1
         KOFF4  = IT1AM(ISYMB,ISYMI)  + 1
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMB),
     *              -ONE,DSAB(KOFF3),NTOTA,XINTAI(KOFF4),NTOTB,ONE,
     *              ETAKA(KOFFAI),NTOTA)

         ! - sum_j g_ja(aj) d_ji ?????????????????????????????????


         KOFF5  = IT1AM(ISYMA,ISYMJ)  + 1
         KOFF6  = IMATIJ(ISYMJ,ISYMI) + 1
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMJ),
     *              ONE,XINTIA(KOFF5),NTOTA,DSIJ(KOFF6),NTOTJ,ONE,
     *              ETAKA(KOFFAI),NTOTA)

         ! - sum_j d_aj g_ij^T 

         KOFF7  = IT1AM(ISYMA,ISYMJ)  + 1
         KOFF8  = IMATIJ(ISYMI,ISYMJ) + 1
C
         CALL DGEMM('N','T',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMJ),
     *              -ONE,DAI(KOFF7),NTOTA,XINTIJ(KOFF8),NTOTI,ONE,
     *              ETAKA(KOFFAI),NTOTA)

         ! - sum_b d_ba^T g_bi

         KOFF9  = IMATAB(ISYMB,ISYMA) + 1
         KOFF10 = IT1AM(ISYMB,ISYMI)  + 1
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMB),
     *              -ONE,DSAB(KOFF9),NTOTB,XINTAI(KOFF10),NTOTB,
     *               ONE,ETAKA(KOFFAI),NTOTA)

         ! sum_b g_ab d_ib(bi)

         KOFF11 = IMATAB(ISYMA,ISYMB) + 1
         KOFF12 = IT1AM(ISYMB,ISYMI)  + 1

         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMB),
     *              ONE,XINTAB(KOFF11),NTOTA,DIA(KOFF12),NTOTB,
     *              ONE,ETAKA(KOFFAI),NTOTA)

         ! sum_j d_ja(aj) g_ji 

         KOFF13 = IT1AM(ISYMA,ISYMJ)  + 1
         KOFF14 = IMATIJ(ISYMJ,ISYMI) + 1
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMJ),
     *              -ONE,DIA(KOFF13),NTOTA,XINTIJ(KOFF14),NTOTJ,
     *              ONE,ETAKA(KOFFAI),NTOTA)

         ! sum_j g_aj d_ij^T

         KOFF15 = IT1AM(ISYMA,ISYMJ)  + 1
         KOFF16 = IMATIJ(ISYMI,ISYMJ) + 1
C
         CALL DGEMM('N','T',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMJ),
     *              ONE,XINTAI(KOFF15),NTOTA,DSIJ(KOFF16),NTOTI,
     *              ONE,ETAKA(KOFFAI),NTOTA)
C

  100 CONTINUE
C
      CALL QEXIT('ccpt_etaai_2e')
C
      RETURN
      END
