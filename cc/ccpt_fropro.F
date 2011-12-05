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
C  /* Deck ccpt_fckd1 */
      SUBROUTINE CCPT_FCKD1(FOCK,ISYFCK,DEN1E,ISYDEN,
     &                      ETAB,ETIJ,ETAI,
     &                      LETAB,LETIJ,LETAI,
     &                      ISYETA,WORK,LWORK)
*
******************************************************************
*
*    PURPOSE:  calculate the frozen core corrections to the rhs 
*              eta^(T)_ij, eta^(T)_ab and eta^(T)_ai of kappa-bar
*              originating from the D^(T)
*              density contracted with the 2-e part of the 
*              "frozen-core" Fock matrix
*      
*    Input : F_pq (XFOCK,ISYFCK), D^(T)_ia (XDEN1E, ISYDEN)
*    Output: eta^(T)_pq (ab,ij,ai, isyeta)
*    OBS: FACTOR 2 TO BE SETTLED
*    Sonia Coriani 02/05-2002
*     
******************************************************************
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <dummy.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, HALF= 0.5D0)
      DIMENSION FOCK(*),DEN1E(*),WORK(LWORK)
      DIMENSION ETAB(*),ETIJ(*),ETAI(*)
      LOGICAL LETIJ,LETAB,LETAI
* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      KSTART = 1
*
*-----------------------------------------
*     (1-Pab) \sum_k D_kb(bk) F_ka(ka)
*-----------------------------------------
*
      IF (LETAB) THEN

        KFCKIA = KSTART       
        DO ISYMA = 1, NSYM
           ISYMI = MULD2H(ISYMA,ISYFCK)
           DO A = 1,NVIR(ISYMA)
             KOFFPQ = IFCVIR(ISYMI,ISYMA) + NORB(ISYMI)*(A - 1) + 1
             KOFFIA = IT1AMT(ISYMI,ISYMA) + NRHF(ISYMI)*(A - 1) + KFCKIA
             CALL DCOPY(NRHF(ISYMI),FOCK(KOFFPQ),1,WORK(KOFFIA),1)
            END DO
        END DO

        DO 100 ISYMB = 1,NSYM
*
           ISYMA  = MULD2H(ISYETA,ISYMB)
           ISYMK  = MULD2H(ISYFCK,ISYMA)
*
           NRHFK  = MAX(NRHF(ISYMK),1)
           NVIRB  = MAX(NVIR(ISYMB),1)
           NVIRA  = MAX(NVIR(ISYMA),1)
*
           KOFF1  = IT1AM(ISYMB,ISYMK)  + 1  
           KOFF2  = IT1AMT(ISYMK,ISYMA)  + KFCKIA     !offset F_ka 
           KOFF3  = IMATAB(ISYMA,ISYMB)  + 1         

           CALL DGEMM('T','T',NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMK),
     &              HALF,WORK(KOFF2),NRHFK,DEN1E(KOFF1),NVIRB,
     &              ONE,ETAB(KOFF3),NVIRA)

           ISYMK  = MULD2H(ISYDEN,ISYMA)
           NRHFK  = MAX(NRHF(ISYMK),1)
           KOFF1  = IT1AM(ISYMA,ISYMK)  + 1  
           KOFF2  = IT1AMT(ISYMK,ISYMB)  + KFCKIA     !offset F_kb 
           KOFF3  = IMATAB(ISYMA,ISYMB)  + 1         

           CALL DGEMM('N','N',NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMK),
     &              -HALF,DEN1E(KOFF1),NVIRA,WORK(KOFF2),NRHFK,
     &              ONE,ETAB(KOFF3),NVIRA)

  100   CONTINUE

      END IF
*
*---------------------------------------
*     (1-P_ij) \sum_c D_jc(cj) F_ci = eta_ij
*---------------------------------------
*
      IF (LETIJ) THEN

        KFCKCI = KSTART       
        DO ISYMI = 1, NSYM
           ISYMC = MULD2H(ISYMI,ISYFCK)
           DO I = 1,NRHF(ISYMI)
             KOFFPQ = IFCRHF(ISYMC,ISYMI) + NORB(ISYMC)*(I - 1) 
     &                                    + NRHF(ISYMC) + 1
             KOFFCI = IT1AM(ISYMC,ISYMI) + NVIR(ISYMC)*(I - 1) + KFCKCI
             CALL DCOPY(NVIR(ISYMC),FOCK(KOFFPQ),1,WORK(KOFFCI),1)
           END DO
        END DO
         
        DO 200 ISYMJ = 1,NSYM

           ISYMI = MULD2H(ISYETA,ISYMJ)
           ISYMC = MULD2H(ISYDEN,ISYMJ)
*
           NRHFI = MAX(NRHF(ISYMI),1)
           NVIRC = MAX(NVIR(ISYMC),1)
*
           KOFF1 = IT1AM(ISYMC,ISYMI)  + KFCKCI  
           KOFF2 = IT1AM(ISYMC,ISYMJ)  + 1
           KOFF3 = IMATIJ(ISYMI,ISYMJ) + 1 
*
           CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMC),
     &              HALF,WORK(KOFF1),NVIRC,DEN1E(KOFF2),NVIRC,
     &              ONE,ETIJ(KOFF3),NRHFI)
*
           ISYMC = MULD2H(ISYFCK,ISYMJ)
           NVIRC = MAX(NVIR(ISYMC),1)
           KOFF1 = IT1AM(ISYMC,ISYMI) + 1
           KOFF2 = IT1AM(ISYMC,ISYMJ) + KFCKCI 
           KOFF3 = IMATIJ(ISYMI,ISYMJ) + 1 
*
           CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMC),
     &              -HALF,DEN1E(KOFF1),NVIRC,WORK(KOFF2),NVIRC,
     &              ONE,ETIJ(KOFF3),NRHFI)
  200    CONTINUE

      END IF
*
*---------------------------------------
* - \sum_l D_la F_li + \sum_c D_ic F_ca 
*---------------------------------------
*
      IF (LETAI) THEN
*
*       Extract F_kj in memory out of F_pq 
*
        KFCKIJ = KSTART

        DO ISYMJ = 1,NSYM
           ISYMK = MULD2H(ISYMJ,ISYFCK)
           DO J = 1,NRHF(ISYMJ)
              KOFF1 = IFCRHF(ISYMK,ISYMJ) + NORB(ISYMK)*(J - 1) + 1
              KOFF2 = IMATIJ(ISYMK,ISYMJ) + NRHF(ISYMK)*(J - 1) + KFCKIJ
              CALL DCOPY(NRHF(ISYMK),FOCK(KOFF1),1,WORK(KOFF2),1)
           END DO
        END DO
 
        if (.true.) then
        DO 300 ISYMI = 1,NSYM
*
           ISYMA  = MULD2H(ISYETA,ISYMI)
           ISYML  = MULD2H(ISYFCK,ISYMI)
*
           NRHFL = MAX(NRHF(ISYML),1)
           NVIRA = MAX(NVIR(ISYMA),1)
*
           KOFF1  = IT1AM(ISYMA,ISYML)  + 1
           KOFF2  = IMATIJ(ISYML,ISYMI) + KFCKIJ 
           KOFF3  = IT1AM(ISYMA,ISYMI)  + 1   
*
           CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYML),
     &              -HALF,DEN1E(KOFF1),NVIRA,WORK(KOFF2),NRHFL,
     &              ONE,ETAI(KOFF3),NVIRA)
*
  300   CONTINUE
        end if
*
*     Extract F_bc in memory out of F_pq 
*
        KFCKAB = KSTART
        DO ISYMC = 1, NSYM
           ISYMB = MULD2H(ISYMC,ISYFCK)
           DO C = 1,NVIR(ISYMC)
              KOFF1 = IFCVIR(ISYMB,ISYMC) + NORB(ISYMB)*(C - 1)
     &                                    + NRHF(ISYMB) + 1
              KOFF2 = IMATAB(ISYMB,ISYMC) + NVIR(ISYMB)*(C - 1) + KFCKAB
              CALL DCOPY(NVIR(ISYMB),FOCK(KOFF1),1,WORK(KOFF2),1)
           END DO
        END DO

*     + \sum_c  F_ca D_ic (ci)

        if (.true.) then
        DO 400 ISYMI = 1,NSYM
*
           ISYMA  = MULD2H(ISYETA,ISYMI)
           ISYMC  = MULD2H(ISYFCK,ISYMA)
*
           NVIRA = MAX(NVIR(ISYMA),1)
           NVIRC = MAX(NVIR(ISYMC),1)
*
           KOFF1  = IT1AM(ISYMC,ISYMI)  + 1
           KOFF2  = IMATAB(ISYMC,ISYMA) + KFCKAB
           KOFF3  = IT1AM(ISYMA,ISYMI)  + 1   
*
           CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMC),
     &              HALF,WORK(KOFF2),NVIRC,DEN1E(KOFF1),NVIRC,
     &              ONE,ETAI(KOFF3),NVIRA)
*
  400    CONTINUE
         end if

      END IF

      RETURN
      END
*------------------------------------------------------------------

C  /* Deck ccpt_etaif */
      SUBROUTINE CCPT_ETAIF(ETAAI,DSAB,DAI,DIA,DSIJ,
     &                      XFIJ,XFJI,XFAI,XFIA,XFII,
     &                      WORK,LWORK,ISYM)
C
C     Written by Sonia Coriani 30/04 - 2002
C
C     Version: 1.0
C
C     Purpose: To calculate eta^(T)_aI-0. ISYM is the symmetry of the
C              integrals and densities. Densities and integrals refer
C              d_pq and g_pq for given gamma and delta
C
#include <implicit.h>
#include <priunit.h>      
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION ETAAI(*), DSAB(*), DAI(*), DIA(*), DSIJ(*)
      DIMENSION XFIJ(*), XFJI(*)
      DIMENSION XFAI(*), XFIA(*), XFII(*), WORK(LWORK)
C
      DO 100 ISYMA = 1,NSYM
C
         ISYMI = ISYMA
         ISYMC = MULD2H(ISYMA,ISYM)
         ISYMK = MULD2H(ISYMA,ISYM)
C
         IF (NRHFFR(ISYMI) .EQ. 0) GOTO 100
C
         NTOTA = MAX(NVIR(ISYMA),1)
         NTOTC = MAX(NVIR(ISYMC),1)
         NTOTK = MAX(NRHF(ISYMK),1)
         NTOTI = MAX(NRHFFR(ISYMI),1)
C
C        !-sum_c dS_ac g_cI
C
         KOFF1 = IMATAB(ISYMA,ISYMC) + 1
         KOFF2 = IT1FRO(ISYMC,ISYMI) + 1
         KOFFR = IT1FRO(ISYMA,ISYMI) + 1
         CALL DGEMM('N','N',NVIR(ISYMA),NRHFFR(ISYMI),NVIR(ISYMC),
     *              -ONE,DSAB(KOFF1),NTOTA,XFAI(KOFF2),NTOTC,ONE,
     *               ETAAI(KOFFR),NTOTA)
C
C        !-sum_c dS_ca g_cI
C
         KOFF1 = IMATAB(ISYMC,ISYMA) + 1
         KOFF2 = IT1FRO(ISYMC,ISYMI) + 1
         KOFFR = IT1FRO(ISYMA,ISYMI) + 1
         CALL DGEMM('T','N',NVIR(ISYMA),NRHFFR(ISYMI),NVIR(ISYMC),
     *              -ONE,DSAB(KOFF1),NTOTC,XFAI(KOFF2),NTOTC,ONE,
     *               ETAAI(KOFFR),NTOTA)
C
C        !-sum_k d_ak g_kI
C
         KOFF3 = IT1AM(ISYMA,ISYMK)  + 1
         KOFF4 = ICOFRO(ISYMK,ISYMI) + 1
         KOFFR = IT1FRO(ISYMA,ISYMI) + 1
         CALL DGEMM('N','N',NVIR(ISYMA),NRHFFR(ISYMI),NRHF(ISYMK),
     *              -ONE,DAI(KOFF3),NTOTA,XFIJ(KOFF4),NTOTK,ONE,
     *               ETAAI(KOFFR),NTOTA)
C
C        !-sum_k d_ka(ak) g_kI
C
         KOFF3 = IT1AM(ISYMA,ISYMK)  + 1
         KOFF4 = ICOFRO(ISYMK,ISYMI) + 1
         KOFFR = IT1FRO(ISYMA,ISYMI) + 1
         CALL DGEMM('N','N',NVIR(ISYMA),NRHFFR(ISYMI),NRHF(ISYMK),
     *              -ONE,DIA(KOFF3),NTOTA,XFIJ(KOFF4),NTOTK,ONE,
     *               ETAAI(KOFFR),NTOTA)
 
  100 CONTINUE
C
      RETURN
      END

*------------------------------------------------------------------

C  /* Deck ccpt_etijf */
      SUBROUTINE CCPT_ETIJF(ETAIJ,DIJ,DAB,DAI,DIA,
     &                      XFIJ,XFJI,XFAI,XFIA,XFII,
     &                      WORK,LWORK,ISYM)
C
C     Written by Sonia Coriani 29/4 - 2002
C
C     Version: 1.0
C
C     Purpose: To calculate the only-frozen-integral contributions to 
C              eta^(T)_iJ-0. ISYM is the symmetry of the
C              integrals and densities 
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION ETAIJ(*), DIJ(*), DAB(*), DAI(*), DIA(*)
      DIMENSION XFIJ(*), XFJI(*)
      DIMENSION XFAI(*), XFIA(*), XFII(*), WORK(LWORK)
C
      DO 100 ISYMI = 1,NSYM
C
         ISYMJ = ISYMI
         ISYMC = MULD2H(ISYMI,ISYM)
         ISYMK = ISYMC
C
         IF (NRHFFR(ISYMJ) .EQ. 0 ) GOTO 100
C
         NTOTI = MAX(NRHF(ISYMI),1)
         NTOTC = MAX(NVIR(ISYMC),1)
         NTOTK = MAX(NRHF(ISYMK),1)
C
C -sum_c d_ic(ci) g_cJ = -sum_c d_ic g_Jc
C
         KOFF1 = IT1AM(ISYMC,ISYMI)  + 1
         KOFF2 = IT1FRO(ISYMC,ISYMJ) + 1
         KOFFR = ICOFRO(ISYMI,ISYMJ) + 1
         CALL DGEMM('T','N',NRHF(ISYMI),NRHFFR(ISYMJ),NVIR(ISYMC),
     *              -ONE,DIA(KOFF1),NTOTC,XFIA(KOFF2),NTOTC,ONE,
     *               ETAIJ(KOFFR),NTOTI)
C
C -sum_c d_ci g_cJ
C
         KOFF1 = IT1AM(ISYMC,ISYMI)  + 1
         KOFF2 = IT1FRO(ISYMC,ISYMJ) + 1
         KOFFR = ICOFRO(ISYMI,ISYMJ) + 1
         CALL DGEMM('T','N',NRHF(ISYMI),NRHFFR(ISYMJ),NVIR(ISYMC),
     *              -ONE,DAI(KOFF1),NTOTC,XFAI(KOFF2),NTOTC,ONE,
     *               ETAIJ(KOFFR),NTOTI)
C
C -sum_a d_ik g_kJ
C
         KOFF3 = IMATIJ(ISYMI,ISYMK) + 1
         KOFF4 = ICOFRO(ISYMK,ISYMJ) + 1
         KOFFR = ICOFRO(ISYMI,ISYMJ) + 1
         CALL DGEMM('N','N',NRHF(ISYMI),NRHFFR(ISYMJ),NRHF(ISYMK),
     *              -ONE,DIJ(KOFF3),NTOTI,XFIJ(KOFF4),NTOTK,ONE,
     *               ETAIJ(KOFFR),NTOTI)
C
C -sum_a d_ki g_kJ
C
         KOFF4 = IMATIJ(ISYMK,ISYMI) + 1
         KOFF5 = ICOFRO(ISYMK,ISYMJ) + 1
         KOFFR = ICOFRO(ISYMI,ISYMJ) + 1
         CALL DGEMM('T','N',NRHF(ISYMI),NRHFFR(ISYMJ),NRHF(ISYMK),
     *              -ONE,DIJ(KOFF4),NTOTK,XFIJ(KOFF5),NTOTK,ONE,
     *               ETAIJ(KOFFR),NTOTI)
C
  100 CONTINUE
C
      RETURN
      END
C
C--------------------------------------------------
C  /* Deck ccptetijf */
      SUBROUTINE CCPTETIJF(ETAIFJ,D1PTIA,CMO,LUN1,WORK,LWORK)
C
C     Written by Sonia Coriani 6/5 - 2002.
C
C     Purpose: calculate the contributions to eta-iJ from D^(T)_lc
C     and integral intermediates on disc.
C     OBS: FACTOR 2 TO BE SETTLED!!!!!
C

      IMPLICIT NONE
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>

      INTEGER LUN1, LWORK
      INTEGER IOFF, ISYMD, ISYJKI, ISYMA, ISYMJK, ISYMI, ISYDEN
      INTEGER KETAJK, KEND1, KSCR1
      INTEGER LWRK1, NTOT, KINIM, KINIM2, NTOTD, NTOTA
      INTEGER NJK, NTOJKI, KOFF1, JK
      INTEGER LWRK3,KEND3, LWRK2, KEND2
      integer icount, icfcco, isyjkii, isym, isym1
      CHARACTER NAME1*8
      dimension icfcco(8,8)

#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL CMO(*), D1PTIA(*), ETAIFJ(*)
      REAL ZERO, ONE, TWO, HALF
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION ZERO, ONE, TWO, HALF
      DOUBLE PRECISION CMO(*), D1PTIA(*), ETAIFJ(*)
#endif

      PARAMETER (ZERO = 0.0D0, TWO = 2.0D0, ONE = 1.0D0, HALF = 0.5D0)
C
      LOGICAL LOCDBG
      PARAMETER(LOCDBG = .FALSE. )

      NAME1 = 'CCFRO1IN'
C
      IF (LOCDBG) THEN
         WRITE(LUPRI,*) 'Warning: I am inside CCPTETIJF '
      END IF
C
C Density is assumed totalsymmetric
C

      ISYDEN = 1
C
C----------------------------------------
C  Construct temporary offset index array 
C----------------------------------------
C
      icount = 0
      do isym = 1, nsym
         icfcco(isym,isym) = icount 
         !icount = icount + nofroo(isym)*nrhf(isym)
         icount = icount + nofroo(isym)*0
      end do

   
C--------------------------------
C
      KETAJK = 1
      KEND1  = KETAJK + NCOFRO(1)
      LWRK1  = LWORK - KEND1
      CALL DZERO(WORK(KETAJK),NCOFRO(1))

      DO 100 ISYMD = 1,NSYM
C
         IF (NBAS(ISYMD) .EQ. 0) GOTO 100
C
         ISYJKI = ISYMD
         ISYMA  = ISYMD
         ISYMI  = MULD2H(ISYDEN,ISYMA)
         ISYMJK = MULD2H(ISYJKI,ISYMI)
C
C----------------------------------
C        Work space allocation one.
C----------------------------------
C
         KSCR1 = KEND1
         KINIM = KSCR1 + NBAS(ISYMD)*NOFROO(ISYJKI)
         KEND2 = KINIM + NVIR(ISYMA)*NOFROO(ISYJKI)
         LWRK2 = LWORK - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
            CALL QUIT('Insufficient memory in CC_ETACOR')
         ENDIF
C
         CALL DZERO(WORK(KSCR1),KEND2)
C
C----------------------------------------------
C        Read integral intermediates from disc.
C----------------------------------------------
C
         NTOT = NOFROO(ISYJKI)*NBAS(ISYMD)
         IOFF = IOFOAO(ISYJKI,ISYMD) + 1
C
         IF (NTOT .GT. 0) THEN
           CALL GETWA2(LUN1,NAME1,WORK(KSCR1),IOFF,NTOT)
         END IF
C
C--------------------------------------
C        Transform AO index to virtual.
C--------------------------------------
C
         KOFF1  = ILMVIR(ISYMD) + 1

         NTOJKI = MAX(NOFROO(ISYJKI),1)
         NTOTD  = MAX(NBAS(ISYMD),1)
C
         !Integrals contain a factor 2 too much!!!
         CALL DGEMM('N','N',NOFROO(ISYJKI),NVIR(ISYMA),NBAS(ISYMD),
     *              HALF,WORK(KSCR1),NTOJKI,CMO(KOFF1),NTOTD,ZERO,
     *              WORK(KINIM),NTOJKI)

C
C----------------------------------------------------------------
C        WORK(KINIM) -> WORK(KINIM + NVIR(ISYMA)*NOFROO(ISYJKI))
C        now contains X_jKia
C----------------------------------------------------------------
C
         KINIM2  = KEND2
         KEND3   = KINIM2 + NOFROO(ISYJKI)*NRHF(ISYMI)
         LWRK3   = LWORK - KEND3
C
C-----------------------------------------------------------------
C        Transform virtual index to occupied with D_ia (ai).
C        WORK(KINIM2) -> WORK(KINIM2 + NRHF(ISYMI)*NOFROO(ISYJKI))
C        contains X_jKii
C-----------------------------------------------------------------
C
         KOFF1  = IT1AM(ISYMA,ISYMI) + 1
         NTOJKI = MAX(NOFROO(ISYJKI),1)
         NTOTA  = MAX(NVIR(ISYMA),1)
C
         CALL DGEMM('N','N',NOFROO(ISYJKI),NRHF(ISYMI),NVIR(ISYMA),
     *              ONE,WORK(KINIM),NTOJKI,D1PTIA(KOFF1),NTOTA,ZERO,
     *              WORK(KINIM2),NTOJKI)
C
C----------------------------------------------------------
C        Realize \sum_i X_jKii = eta_jK
C----------------------------------------------------------
C
         DO I = 1, NRHF(ISYMI)

            NJK = KINIM2 !+ ICFCCO(ISYJKI,ISYMI)
     &                   + NOFROO(ISYJKI)*(I-1)
     &                   + IOFROO(ISYMJK,ISYMI)
     &                   + NCOFRO(ISYMJK)*(I-1) 

           CALL DAXPY(NCOFRO(ISYMJK), ONE, WORK(NJK),1, 
     &                    WORK(KETAJK), 1)
         END DO
C----------------------------------------------------------

  100 CONTINUE

      if (.true.) then
      CALL DAXPY(NCOFRO(1),ONE,WORK(KETAJK),1,ETAIFJ,1)
      end if
C
      RETURN
      END
C---------------------------------
C  /* Deck ccpt_zkfcb */
      SUBROUTINE CCPT_ZKFCB(ETAIF,D1PTIA,WORK,LWORK)
C
C     Written by Sonia Coriani Halkier 2002.
C
C     To calculate D^(T)*P_integrals contributions to eta_aI
C
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0)
      PARAMETER (FOUR = 4.0D0)
      DIMENSION ETAIF(*), D1PTIA(*), WORK(LWORK)
      
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE. )
      
C
C
C---------------------------
C     Work space allocation.
C---------------------------
C
      KINTFR1 = 1
      KINTFR2 = KINTFR1 + NT2FRO(1)
      KINTFR3 = KINTFR2 + NT2FRO(1)
      KEND1   = KINTFR3 + NT2FRO(1)
      LWRK1   = LWORK   - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient work space for allocation '//
     &             'in CCPT_ZKFCB')
      ENDIF
C
C--------------------------------------
C     Read integrals (cJ|dk) from disk.
C--------------------------------------
C
      LUCJDK = -1
      CALL GPOPEN(LUCJDK,'INCJDK','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUCJDK)
      READ(LUCJDK) (WORK(KINTFR2+I-1), I = 1,NT2FRO(1))
      CALL GPCLOSE(LUCJDK,'KEEP')
C
      IF (LOCDBG) THEN
         XFRNOR = DDOT(NT2FRO(1),WORK(KINTFR2),1,WORK(KINTFR2),1)
         WRITE(LUPRI,*) 'Norm of integrals (cJdk):', XFRNOR
      ENDIF

      CALL DCOPY(NT2FRO(1),WORK(KINTFR2),1,WORK(KINTFR1),1)
      CALL DSCAL(NT2FRO(1),FOUR,WORK(KINTFR1),1)
C
C--------------------------------------
C     Read integrals (cd|kJ) from disk.
C--------------------------------------
C
      LUCDKJ = -1
      CALL GPOPEN(LUCDKJ,'INCDKJ','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUCDKJ)
      READ(LUCDKJ) (WORK(KINTFR3+I-1), I = 1,NT2FRO(1))
      CALL GPCLOSE(LUCDKJ,'KEEP')

      IF (LOCDBG) THEN
         XFRNOR = DDOT(NT2FRO(1),WORK(KINTFR3),1,WORK(KINTFR3),1)
         WRITE(LUPRI,*) 'Norm of integrals (cdkJ):', XFRNOR
      ENDIF
C
C------------------------------------------------------------------
C     Construct integral intermediate X_cJ,dk = 4 I_cJ,dk - I_ck,dJ 
C                                             = 4 I_cJ,dk - I_dJ,ck
C     (Sort of Coulomb - Exchange)
C------------------------------------------------------------------
C
      ! Gli integrali CDKJ potrei aggiungerli a blocchi di NVIR(ISYMC
      ! Oppure combinare insieme i DJ e i CD (in fondo CDKJ e DCKJ 
      ! sono la stessa cosa)
      !
      DO ISYMK  = 1, NSYM
         ISYMD  = ISYMK
         ISYMDK = MULD2H(ISYMD,ISYMK)
         ISYMCJ = ISYMDK
         DO ISYMJ = 1, NSYM
            ISYMKJ = MULD2H(ISYMK,ISYMJ)
            ISYMCD = ISYMKJ
            ISYMC  = ISYMJ
            ISYMCK = MULD2H(ISYMC,ISYMK)
            ISYMDJ = MULD2H(ISYMD,ISYMJ)

         DO K = 1, NRHF(ISYMK)
         DO D = 1, NVIR(ISYMD)
            DO J = 1, NRHFFR(ISYMJ)
               DO C = 1, NVIR(ISYMC)

!                  write(lupri,*) 'K, D, J, C : ', K, D, J, C

                  NDK   = IT1AM(ISYMD,ISYMK) + NVIR(ISYMD)*(K-1) + D
                  NCK   = IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K-1) + C
                  NKJ   = ICOFRO(ISYMK,ISYMJ)+ NRHF(ISYMK)*(J-1) + K

!                  write(lupri,*) 'NDK, NCK, NKJ : ', NDK, NCK, NKJ

                  KCJDK = IT2FRO(ISYMCJ,ISYMDK) 
     &                    + NT1FRO(ISYMCJ)*(NDK-1) + IT1FRO(ISYMC,ISYMJ)
     &                    + NVIR(ISYMC)*(J-1) + C
                  KDJCK = IT2FRO(ISYMDJ,ISYMCK) 
     &                    + NT1FRO(ISYMDJ)*(NCK-1) + IT1FRO(ISYMD,ISYMJ)
     &                    + NVIR(ISYMD)*(J-1) + D

                  KCDKJ = ICDKFR(ISYMCD,ISYMKJ) 
     &                    + NMATAB(ISYMCD)*(NKJ-1) + IMATAB(ISYMC,ISYMD)
     &                    + NVIR(ISYMC)*(D-1) + C

                  WORK(KINTFR1+KCJDK-1) =  WORK(KINTFR1+KCJDK-1)
     &                                   - WORK(KINTFR2+KDJCK-1)
     &                                   - WORK(KINTFR3+KCDKJ-1)
               END DO
            END DO
         END DO
         END DO
      END DO
      END DO
C
C----------------------------------------------
C     Contract integrals with density D_kd (dk)
C     sum_dk I_cJ,dk D_dk
C     As D is totalsym, isymdk = 1
C     As a consequence  isymcj = 1
C----------------------------------------------
C
      KETAI = KEND1
      KEND1 = KETAI + NT1FRO(1)
      LWRK1 = LWORK  - KEND1
      CALL DZERO(WORK(KETAI),NT1FRO(1))

      ISYMCJ = 1
      ISYMDK = ISYMCJ
C
      KOFF1 = IT2FRO(ISYMCJ,ISYMDK) + KINTFR1
      KOFF2 = 1
C
         NTOTCJ = MAX(NT1FRO(ISYMCJ),1)
C
         CALL DGEMV('N',NT1FRO(ISYMCJ),NT1AM(ISYMDK),
     *               ONE,WORK(KOFF1),NTOTCJ,D1PTIA(KOFF2),1,
     *               ONE,WORK(KETAI),1)
C

      CALL DAXPY(NT1FRO(1),ONE,WORK(KETAI),1,ETAIF(1),1)
C
      RETURN
      END
C
C--------------------
C  /* Deck ccpt_frin */
      SUBROUTINE CCPT_FRIN(XFRIN,XINT,CMO,WORK,LWORK,IDEL,ISYDEL)
C
C     Written by Sonia Coriani Spring 2002
C
C     To calculate the integrals (cd|kJ) needed for frozen-core gradients.
C     AO integrals are assumed totally symmetric.
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION XFRIN(*), XINT(*), CMO(*), WORK(LWORK)
C
C
C---------------------------------------
C     Initial symmetries and allocation.
C---------------------------------------
C
      ISYDIS = ISYDEL
      ISYMJ  = ISYDEL
C
      KJVEC  = 1
      KEND1  = KJVEC + NRHFFR(ISYMJ)
      LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient work space for allocation '//
     &             'in CCPT_FRIN')
      ENDIF
C
      CALL DZERO(WORK(KJVEC),NRHFFR(ISYMJ))
C
C-----------------------------------
C     Copy vector out of CMO matrix.
C-----------------------------------
C
      KOFF1 = ILRHSI(ISYMJ) + IDEL - IBAS(ISYDEL)
C
C Capire bene quanti ne deve saltare !!
C
      CALL DCOPY(NRHFFR(ISYMJ),CMO(KOFF1),NBAS(ISYDEL),WORK(KJVEC),1)
C
C--------------------------------------------------------
C     Outer symmetry loop and next work space allocation.
C--------------------------------------------------------
C
      DO 100 ISYMK = 1,NSYM
C
         ISYMG  = ISYMK
         ISYMKJ = MULD2H(ISYMJ,ISYMK)
         ISYMCD = ISYMKJ
         ISALBE = ISYMCD 
C
         KINAOK = KEND1
         KSCR1  = KINAOK + NNBST(ISALBE)*NRHF(ISYMK)
         KSCR2  = KSCR1  + N2BST(ISALBE)
         KEND2  = KSCR2  + NMATAB(ISYMCD)
         LWRK2  = LWORK  - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
            CALL QUIT('Insufficient space for allocation in CCPT_FRIN')
         ENDIF
C
         CALL DZERO(WORK(KINAOK),NNBST(ISALBE)*NRHF(ISYMK))
C
C-------------------------------------------
C        Transform gamma-index to correlated
C-------------------------------------------
C
         KOFF2  = IDSAOG(ISYMG,ISYDIS) + 1
         KOFF3  = ILRHSI(ISYMG) + NBAS(ISYMG)*NRHFFR(ISYMK) + 1
C
         NTOTAB = MAX(NNBST(ISALBE),1)
         NTOTG  = MAX(NBAS(ISYMG),1)
C
         CALL DGEMM('N','N',NNBST(ISALBE),NRHF(ISYMK),NBAS(ISYMG),ONE,
     *              XINT(KOFF2),NTOTAB,CMO(KOFF3),NTOTG,ZERO,
     *              WORK(KINAOK),NTOTAB)
C
         DO 110 K = 1,NRHF(ISYMK)
C
            CALL DZERO(WORK(KSCR1),N2BST(ISALBE))
            CALL DZERO(WORK(KSCR2),NMATAB(ISYMCD))
C
C-----------------------------------
C           Square up the integrals.
C-----------------------------------
C
            KOFF4 = KINAOK + NNBST(ISALBE)*(K - 1)
            CALL CCSD_SYMSQ(WORK(KOFF4),ISALBE,WORK(KSCR1))
C
C---------------------------------------------------------------
C           Inner symmetry loop and final work space allocation.
C---------------------------------------------------------------
C
            DO 120 ISYMAL = 1,NSYM
C
               ISYMC  = ISYMAL
               ISYMBE = MULD2H(ISALBE,ISYMAL)
               ISYMD  = ISYMBE
C
               KSCR3 = KEND2
               KEND3 = KSCR3 + NBAS(ISYMAL)*NVIR(ISYMD)
               LWRK3 = LWORK - KEND3
C
               IF (LWRK3 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND3
                  CALL QUIT('Insufficient space for allocation '//
     &                      'in CC_FRCOIN')
               ENDIF
C
C------------------------------------------------
C              Construct the integrals (cd|kdel).
C------------------------------------------------
C
               KOFF5  = KSCR1 + IAODIS(ISYMAL,ISYMBE)
               KOFF6  = ILVISI(ISYMBE) + 1
C
               NTOTAL = MAX(NBAS(ISYMAL),1)
               NTOTBE = MAX(NBAS(ISYMBE),1)
C
               CALL DGEMM('N','N',NBAS(ISYMAL),NVIR(ISYMD),
     *                    NBAS(ISYMBE),ONE,WORK(KOFF5),NTOTAL,
     *                    CMO(KOFF6),NTOTBE,ZERO,WORK(KSCR3),NTOTAL)
C
               KOFF7  = ILVISI(ISYMAL) + 1
               KOFF8  = KSCR2 + IMATAB(ISYMC,ISYMD)
C
               NTOTAL = MAX(NBAS(ISYMAL),1)
               NTOTC  = MAX(NVIR(ISYMC),1)
C
               CALL DGEMM('T','N',NVIR(ISYMC),NVIR(ISYMD),
     *                    NBAS(ISYMAL),ONE,CMO(KOFF7),NTOTAL,
     *                    WORK(KSCR3),NTOTAL,ZERO,WORK(KOFF8),NTOTC)
C
  120       CONTINUE
C
C----------------------------------------------------------------
C           Final scaling with CMO element and storage in result.
C----------------------------------------------------------------
C
            DO 130 J = 1,NRHFFR(ISYMJ)
C
               NKJ    = ICOFRO(ISYMK,ISYMJ) + NRHF(ISYMK)*(J - 1) + K
               KOFF9  = KJVEC + J - 1
               KOFF10 = ICDKFR(ISYMCD,ISYMKJ)
     *                       + NMATAB(ISYMCD)*(NKJ - 1) + 1
C
               CALL DAXPY(NMATAB(ISYMCD),WORK(KOFF9),WORK(KSCR2),1,
     *                    XFRIN(KOFF10),1)
C
  130       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C-------------------------------------------------------------------
C  /* Deck ccfckmfrl */
      SUBROUTINE CCFCKMFRL(ETAAIF,ETAIJF,ISYETA,
     &                     FOCKAO,ISYFAO,
     &                     D1PTIA,ISYDEN,
     &                     CMOF,WORK,LWORK)
C
C     Written by Sonia Coriani
C     Purpose: Calculate special MO Fock Matrix block F_mI
C     \sum_\a\b CMO_\am F_\a\b CMO_\bI
C     and contract with D^(T)_ma to get eta_aI
C
#include <implicit.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
      PARAMETER (ONE = 1.0D0, HALF= 0.5D0, ZERO = 0.0D0)
      DIMENSION ETAAIF(*),ETAIJF(*) 
      DIMENSION FOCKAO(*), D1PTIA(*), CMOF(*),WORK(LWORK)
C

      DO 100 ISYMI = 1,NSYM
C
         ISYMA = MULD2H(ISYETA,ISYMI)

         ISYMBE = ISYMI
         ISYMAL = MULD2H(ISYMBE,ISYFAO)
         ISYMM  = ISYMI
         ISYMMI = MULD2H(ISYMI,ISYMM)

         KSCR1  = 1
         KFCKMI = KSCR1  + NBAS(ISYMAL)*NRHFFR(ISYMI)
         KEND2  = KFCKMI + NCOFRO(ISYMMI)
         LWRK2  = LWORK - KEND2

         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
            CALL QUIT('Insufficient memory for allocation '//
     &                   'in CCFCKMFRI')
         ENDIF
C
         CALL DZERO(WORK(KSCR1),NBAS(ISYMAL)*NRHFFR(ISYMI))
         CALL DZERO(WORK(KFCKMI),NCOFRO(ISYMMI))
C
C----------------------------------------
C        Do first partial transformation.
C----------------------------------------
C
         NTOTBE = MAX(NBAS(ISYMBE),1)
         NTOTAL = MAX(NBAS(ISYMAL),1)
C
         KOFF1 = IAODIS(ISYMAL,ISYMBE) + 1
         KOFF2 = ILRHSI(ISYMBE) + 1
C
         CALL DGEMM('N','N',NBAS(ISYMAL),NRHFFR(ISYMI),NBAS(ISYMBE),
     *              ONE,FOCKAO(KOFF1),NTOTAL,CMOF(KOFF2),NTOTBE,
     *              ZERO,WORK(KSCR1),NTOTAL)
C
C-----------------------------------------
C        Do second partial transformation.
C-----------------------------------------
C
         NTOTM  = MAX(NRHF(ISYMM),1)
         KOFRES = KFCKMI + ICOFRO(ISYMM,ISYMI) 
         KOFF3  = ILRHSI(ISYMAL) + NBAS(ISYMAL)*NRHFFR(ISYMI) + 1
C
         CALL DGEMM('T','N',NRHF(ISYMM),NRHFFR(ISYMI),NBAS(ISYMAL),ONE,
     *              CMOF(KOFF3),NTOTAL,WORK(KSCR1),NTOTAL,ONE,
     *              WORK(KOFRES),NTOTM)
C
C-----------------------------------------
C        Do contraction and add to result
C-----------------------------------------
C

         NTOTM  = MAX(NRHF(ISYMM),1)
         NTOTA  = MAX(NVIR(ISYMA),1)
         KOFF4  = IT1AM(ISYMA,ISYMM) + 1
         KOFF5  = KFCKMI + ICOFRO(ISYMM,ISYMI) 
         KOFRE1 = IT1FRO(ISYMA,ISYMI) + 1
C
         if (.true.) then
         CALL DGEMM('N','N',NVIR(ISYMA),NRHFFR(ISYMI),NRHF(ISYMM),-HALF,
     *              D1PTIA(KOFF4),NTOTA,WORK(KOFF5),NTOTM,ONE,
     *              ETAAIF(KOFRE1),NTOTA)
         end if
C
  100 CONTINUE
C
      DO 110 ISYMJ = 1,NSYM
C
         ISYMI = MULD2H(ISYETA,ISYMJ)

         ISYMBE = ISYMJ
         ISYMAL = MULD2H(ISYMBE,ISYFAO)
         ISYMC  = ISYMJ
         ISYMCJ = MULD2H(ISYMC,ISYMJ)

         KSCR1  = 1
         KFCKCJ = KSCR1  + NBAS(ISYMAL)*NRHFFR(ISYMJ)
         KEND2  = KFCKCJ + NT1FRO(ISYMCJ)
         LWRK2  = LWORK  - KEND2

         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
            CALL QUIT('Insufficient memory for allocation '//
     &                   'in CCFCKMFRI')
         ENDIF
C
         CALL DZERO(WORK(KSCR1),NBAS(ISYMAL)*NRHFFR(ISYMJ))
         CALL DZERO(WORK(KFCKCJ),NT1FRO(ISYMCJ))
C
C----------------------------------------
C        Do first partial transformation.
C        \sum_beta F_\alpha,\beta C_\beta,J
C----------------------------------------
C
         NTOTBE = MAX(NBAS(ISYMBE),1)
         NTOTAL = MAX(NBAS(ISYMAL),1)
C
         KOFF1 = IAODIS(ISYMAL,ISYMBE) + 1
         KOFF2 = ILRHSI(ISYMBE) + 1
C
         CALL DGEMM('N','N',NBAS(ISYMAL),NRHFFR(ISYMJ),NBAS(ISYMBE),
     *              ONE,FOCKAO(KOFF1),NTOTAL,CMOF(KOFF2),NTOTBE,
     *              ZERO,WORK(KSCR1),NTOTAL)
C
C-----------------------------------------
C        Do second partial transformation.
C        \sum_alpha F_\alpha,J C_\alpha,c
C-----------------------------------------
C
         NTOTC  = MAX(NVIR(ISYMC),1)
         KOFRES = KFCKCJ + IT1FRO(ISYMC,ISYMJ) 
         KOFF3  = ILVISI(ISYMAL) + 1
C
         CALL DGEMM('T','N',NVIR(ISYMC),NRHFFR(ISYMJ),NBAS(ISYMAL),ONE,
     *              CMOF(KOFF3),NTOTAL,WORK(KSCR1),NTOTAL,ONE,
     *              WORK(KOFRES),NTOTC)
C
C-----------------------------------------
C        Do contraction and add to result
C        - sum_c D_ic (ci) F_cJ
C-----------------------------------------
C

         NTOTI  = MAX(NRHF(ISYMI),1)
         NTOTC  = MAX(NVIR(ISYMC),1)
         KOFF4  = IT1AM(ISYMC,ISYMI) + 1
         KOFF5  = KFCKCJ + IT1FRO(ISYMC,ISYMJ) 
         KOFRE  = ICOFRO(ISYMI,ISYMJ) + 1
C
         if (.true.) then
         CALL DGEMM('T','N',NRHF(ISYMI),NRHFFR(ISYMJ),NVIR(ISYMC),-HALF,
     *              D1PTIA(KOFF4),NTOTC,WORK(KOFF5),NTOTC,ONE,
     *              ETAIJF(KOFRE),NTOTI)
         end if
C
  110 CONTINUE

      RETURN
      END
C
C----------------------------------------------------------------
C  /* Deck ccpt_fd2bl */
      SUBROUTINE CCPT_FD2BL(D2II,D2IJ,D2JI,D2AI,D2IA,D1IA,
     *                    CMO,CMOF,WORK,LWORK,IDEL,ISYMD,
     *                    G,ISYMG)
C
C     Written by Sonia Coriani Halkier 14/5/2002
C
C     Purpose: To calculate the contributions to d(T)(pq,gam;del) 
C              where at least one of the two indices p & q is frozen.
C
#include <implicit.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION D2II(*), D2IJ(*), D2JI(*), D2AI(*), D2IA(*)
      DIMENSION D1IA(*), CMO(*), CMOF(*)
      DIMENSION WORK(LWORK)
C
      CALL QENTER('CCPT_FD2BL')
C
C-------------------------------
C     Work space allocation one.
C-------------------------------
C
      ISYMA = ISYMD
      ISYML = ISYMD
      ISYMI = ISYMA
      ISYMJ = ISYMG
C
      KVECI = 1
      KVECA = KVECI + NRHF(ISYMI)
      KVECL = KVECA + NVIR(ISYMA)
      KEND1 = KVECL + NRHF(ISYML)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient work space in CCPT_FD2BL')
      ENDIF
C
      CALL DZERO(WORK(KVECI),KEND1)
C
      KOFF9  = IGLMVI(ISYMD,ISYMA) + IDEL - IBAS(ISYMD)
      KOFF10 = IGLMRH(ISYMD,ISYML) + IDEL - IBAS(ISYMD)
      CALL DCOPY(NVIR(ISYMA),CMO(KOFF9),NBAS(ISYMD),WORK(KVECA),1)
      CALL DCOPY(NRHF(ISYML),CMO(KOFF10),NBAS(ISYMD),WORK(KVECL),1)
C
C---------------------------------------------------------------
C     Calculate intermediate vector V_i = \sum_a D_ia C_a;delta.
C---------------------------------------------------------------
C
      KOFF11 = IT1AM(ISYMA,ISYMI)  + 1
      KOFF12 = IMATIJ(ISYMI,ISYML) + 1
C
      NTOTA  = MAX(NVIR(ISYMA),1)
      NTOTI  = MAX(NRHF(ISYMI),1)
C
      CALL DGEMV('T',NVIR(ISYMA),NRHF(ISYMI),ONE,D1IA(KOFF11),NTOTA,
     *           WORK(KVECA),1,ZERO,WORK(KVECI),1)
C
C--------------------------------------------
C     Calculate correlated-frozen block (iJ).
C--------------------------------------------
C
      DO 130 J = 1,NRHFFR(ISYMJ)
C
         KOFF13 = ILRHSI(ISYMG) + NBAS(ISYMG)*(J - 1) + G
         KOFF14 = ICOFRO(ISYMI,ISYMJ) + NRHF(ISYMI)*(J - 1) + 1
         CALL DAXPY(NRHF(ISYMI),-CMOF(KOFF13),WORK(KVECI),1,
     *              D2IJ(KOFF14),1)
C
  130 CONTINUE
C
C----------------------------------------------------------
C     Add contribution to frozen-frozen block (II) from V_a.
C----------------------------------------------------------
C
      IF (ISYMG .EQ. ISYMD) THEN
C
         ISYMK = ISYMD
         KVECK = KVECI
C
         KOFF15 = IGLMRH(ISYMG,ISYMK) + G
         FAC = DDOT(NRHF(ISYMK),WORK(KVECK),1,CMO(KOFF15),
     *              NBAS(ISYMG))
C
         DO 140 ISYMM = 1,NSYM
            DO 150 M = 1,NRHFFR(ISYMM)
               KOFF16 = IFROFR(ISYMM,ISYMM) + NRHFFR(ISYMM)*(M - 1) + M
               D2II(KOFF16) = D2II(KOFF16) + TWO*FAC
  150       CONTINUE
  140    CONTINUE
      ENDIF
C
C----------------------------------------------------------------------
C     Calculate Hartree-Fock like contributions to frozen-frozen block.
C     IJ and II
C----------------------------------------------------------------------
C
      if (.false.) then
      IF (ISYMG .EQ. ISYMD) THEN
C
         KOFF17 = ILRHSI(ISYMG) + G
         KOFF18 = ILRHSI(ISYMD) + IDEL - IBAS(ISYMD)
         FAC = TWO*DDOT(NRHFFR(ISYMG),CMOF(KOFF17),NBAS(ISYMG),
     *              CMOF(KOFF18),NBAS(ISYMD))
C
         DO 160 ISYMI = 1,NSYM
            DO 170 I = 1,NRHFFR(ISYMI)
               KOFF19 = IFROFR(ISYMI,ISYMI) + NRHFFR(ISYMI)*(I - 1) + I
               D2II(KOFF19) = D2II(KOFF19) + TWO*FAC
  170       CONTINUE
  160    CONTINUE
      ENDIF
C
      ISYMI = ISYMD
      ISYMJ = ISYMG
      D     = IDEL - IBAS(ISYMD)
C
      DO 180 J = 1,NRHFFR(ISYMJ)
         DO 190 I = 1,NRHFFR(ISYMI)
            KOFF20 = IFROFR(ISYMI,ISYMJ) + NRHFFR(ISYMI)*(J - 1) + I
            KOFF21 = ILRHSI(ISYMD) + NBAS(ISYMD)*(I - 1) + D
            KOFF22 = ILRHSI(ISYMG) + NBAS(ISYMG)*(J - 1) + G
            D2II(KOFF20) = D2II(KOFF20) - TWO*CMOF(KOFF21)*CMOF(KOFF22)
  190    CONTINUE
  180 CONTINUE
      end if
C
C---------------------------------
C     Work space allocation two.
C---------------------------------
C
      ISYMB = ISYMG
      ISYMJ = ISYMG
      ISYMA = ISYMB
      ISYMI = ISYMD
C
      KVECA = 1
      KVECB = KVECA + NVIR(ISYMA)
      KVECJ = KVECB + NVIR(ISYMB)
      KEND1 = KVECJ + NRHF(ISYMJ)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient work space in CCPT_FD2BL')
      ENDIF
C
      CALL DZERO(WORK(KVECI),KEND1)
C
      KOFF23 = IGLMVI(ISYMG,ISYMB) + G
      KOFF24 = IGLMRH(ISYMG,ISYMJ) + G
      CALL DCOPY(NVIR(ISYMB),CMO(KOFF23),NBAS(ISYMG),WORK(KVECB),1)
      CALL DCOPY(NRHF(ISYMJ),CMO(KOFF24),NBAS(ISYMG),WORK(KVECJ),1)
C
C--------------------------------------
C     Calculate intermediate vector Wa.
C--------------------------------------
C
      KOFF25 = IMATAB(ISYMB,ISYMA) + 1
      KOFF26 = IT1AM(ISYMA,ISYMJ)  + 1
C
      NTOTB  = MAX(NVIR(ISYMB),1)
      NTOTA  = MAX(NVIR(ISYMA),1)
C
      CALL DGEMV('N',NVIR(ISYMA),NRHF(ISYMJ),ONE,D1IA(KOFF26),NTOTA,
     *           WORK(KVECJ),1,ONE,WORK(KVECA),1)
C
C----------------------------------------
C     Calculate frozen-virtual block (Ia).
C--------------------------------------------
C
      DO 200 I = 1,NRHFFR(ISYMI)
C
         KOFF27 = ILRHSI(ISYMD) + NBAS(ISYMD)*(I - 1)
     *          + IDEL - IBAS(ISYMD)
         KOFF28 = IT1FRO(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + 1
         CALL DAXPY(NVIR(ISYMA),-CMOF(KOFF27),WORK(KVECA),1,
     *              D2IA(KOFF28),1)
C
  200 CONTINUE
C
      CALL QEXIT('CCPT_FD2BL')
C
      RETURN
      END
C-------------------------------------------------
C  /* Deck ccpt_d2gaf */
      SUBROUTINE CCPT_D2GAF(D2GIA,D1IA,
     *                      CMOF,IDEL,ISYMD,G,ISYMG)
C
C     Written by Sonia Coriani Halkier 12/5 - 2002
C
C     Purpose: To calculate the contributions to d(T)(pq,gam;del) where
C              gamma has been backtransformed from a frozen index.
C
#include <implicit.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION D2GIA(*), D1IA(*), CMOF(*)
C
      CALL QENTER('CCPT_D2GAF')
C
      IF (ISYMG .EQ. ISYMD) THEN
C
         ISYML = ISYMD
C
         ND = ILRHSI(ISYMD) + IDEL - IBAS(ISYMD)
         NG = ILRHSI(ISYMG) + G
C
         FACT = TWO*DDOT(NRHFFR(ISYML),CMOF(ND),NBAS(ISYMD),CMOF(NG),
     *               NBAS(ISYMG))
C
         CALL DAXPY(NT1AMX,FACT,D1IA,1,D2GIA,1)
C 
      ENDIF
C
      CALL QEXIT('CCPT_D2GAF')
C
      RETURN
      END
C----------------------------
C  /* Deck cc_etfro1e */
      SUBROUTINE CC_ETFRO1E(ETAFI,ETIFJ,DIA,
     *                      XFIJ,XFJI,XFAI, XFIA,
     *                      WORK,LWORK,ISYM)
C
C     Written by Sonia Coriani 2002
C
C     Version: 1.0
C
C     Purpose: To calculate D^(T)_ia contributions to eta-iJ
C              and eta_aI. ISYM is the symmetry of the
C              integrals and densities.
C              IOPT controls the one
C              and two electron contributions.
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, HALF = 0.5D0)
      DIMENSION ETAFI(*), ETIFJ(*), DIA(*)
      DIMENSION XFIJ(*), XFJI(*)
      DIMENSION XFAI(*), XFIA(*), WORK(LWORK)
C
!      IF (LETIJ) THEN
         DO 100 ISYMI = 1,NSYM
C
            ISYMC = MULD2H(ISYMI,ISYM)
            ISYMJ = MULD2H(ISYMC,ISYM)
            IF (NRHFFR(ISYMJ) .EQ. 0 ) GOTO 100
C
            KOFFR = ICOFRO(ISYMI,ISYMJ) + 1
            KOFF1 = IT1AM(ISYMC,ISYMI)  + 1
            KOFF2 = IT1FRO(ISYMC,ISYMJ) + 1
C
            NTOTI = MAX(NRHF(ISYMI),1)
            NTOTC = MAX(NVIR(ISYMC),1)
C
            CALL DGEMM('T','N',NRHF(ISYMI),NRHFFR(ISYMJ),NVIR(ISYMC),
     *              -one,DIA(KOFF1),NTOTC,XFIA(KOFF2),NTOTC,ONE,
     *               ETIFJ(KOFFR),NTOTI)
C
  100    CONTINUE
!      END IF
C
!      IF (LETAI) THEN
         DO 101 ISYMI = 1,NSYM
C
            ISYMK = MULD2H(ISYMI,ISYM)
            ISYMA = MULD2H(ISYMK,ISYM)
            IF (NRHFFR(ISYMI) .EQ. 0 ) GOTO 101
C
            KOFFR = IT1FRO(ISYMA,ISYMI) + 1
            KOFF1 = IT1AM(ISYMA,ISYMK)  + 1
            KOFF2 = ICOFRO(ISYMK,ISYMI) + 1
C
            NTOTK = MAX(NRHF(ISYMK),1)
            NTOTA = MAX(NVIR(ISYMA),1)
C
            CALL DGEMM('N','N',NVIR(ISYMA),NRHFFR(ISYMI),NRHF(ISYMK),
     *              -one,DIA(KOFF1),NTOTA,XFJI(KOFF2),NTOTK,ONE,
     *               ETAFI(KOFFR),NTOTA)
C
  101    CONTINUE
!      END IF
      RETURN
      END
