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
C  /* Deck cc_gam */
      SUBROUTINE CC_GAM(DSRHF,GAMMA,XLAMDP,XLAMDH,
     *                  XLAMPC,XLAMHC,ISYMPC,SCRM,ISYMCM,
     *                  WORK,LWORK,IDEL,ISYMD,IOPT)
C
C     Written by Henrik Koch 3-Jan-1994
C     Symmetry by Henrik Koch and Alfredo Sanchez. 21-July-1994
C
C     Purpose: Calculate the gamma intermediate.
C
C     Ove Christiansen 18-9-1996:
C              General symmetry for F-matrix construction.
C
#include <implicit.h>
      DIMENSION DSRHF(*),GAMMA(*),SCRM(*)
      DIMENSION WORK(LWORK)
      DIMENSION XLAMDP(*),XLAMDH(*),XLAMPC(*),XLAMHC(*)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
C------------------------
C     Dynamic allocation.
C------------------------
C
      ISYMJB = MULD2H(ISYMD,ISYMPC)
      KLAMDA = 1
      KLAMDC = KLAMDA + NRHF(ISYMD)
      KEND1  = KLAMDC + NRHF(ISYMJB)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Need : ',KEND1,'Available : ',LWORK
         CALL QUIT('Insufficient space in CC_GAM')
      ENDIF
C
C---------------------------------------
C     Copy XLAMDH vector for given IDEL.
C---------------------------------------
C
      KOFF1 = ILMRHF(ISYMD) + IDEL - IBAS(ISYMD)
      CALL DCOPY(NRHF(ISYMD),XLAMDH(KOFF1),NBAS(ISYMD),WORK(KLAMDA),1)
C
C--------------------------------
C     Calculate the contribution.
C--------------------------------
C
      ISYDIS = MULD2H(ISYMD,ISYMOP)
C
      DO 100 ISYML = 1,NSYM
C
         ISYMAG = MULD2H(ISYML,ISYDIS)
         ISYMKI = MULD2H(ISYMAG,ISYMPC)
C
C---------------------------
C        Dynamic allocation.
C---------------------------
C
         KSCR1  = KEND1
         KSCR2  = KSCR1  + N2BST(ISYMAG)
         KSCR3  = KSCR2  + NT1AO(ISYMAG)
         KSCR4  = KSCR3  + NT1AM(ISYMAG)
         KSCR5  = KSCR4  + NMATIJ(ISYMAG)
C
         IF (IOPT .EQ. 1) THEN
            KEND2  = KSCR5  
         ELSE 
            KEND2  = KSCR5  + NMATIJ(ISYMKI)
         ENDIF
C
         LWRK2  = LWORK  - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Need : ',KEND2,'Available : ',LWORK
            CALL QUIT('Insufficient space in CC_GAM')
         ENDIF
C
         CALL CC_GAM1(DSRHF,GAMMA,SCRM,ISYMCM,WORK(KLAMDA),WORK(KLAMDC),
     *                XLAMDP,XLAMDH,XLAMPC,XLAMHC,ISYMPC,
     *                WORK(KSCR1),WORK(KSCR2),WORK(KSCR3),WORK(KSCR4),
     *                WORK(KSCR5),WORK(KEND2),LWRK2,ISYMD,IDEL,ISYML,
     *                ISYMAG,IOPT)
C
  100 CONTINUE
C
      RETURN
      END
      SUBROUTINE CC_GAM1(DSRHF,GAMMA,SCRM,ISYMCM,XLAM,XLAMC,
     *                   XLAMDP,XLAMDH,XLAMPC,XLAMHC,ISYMPC,
     *                   SCR1,SCR2,SCR3,SCR4,SCR5,WORK,
     *                   LWORK,ISYMD,IDEL,ISYML,ISYMAG,IOPT)
C
C     Written by Henrik Koch 3-Jan-1994
C
C     Generalized by Ove Christiansen 18-9-1996 for 
C     calculation of Gamma-intermediate for F matrix.
C     F-matrix: IOPT = 2 and ISYMCM is symmetry of SCRM
C               and ISYMPC is symmetry of XLAMPC and XLAMHC.
C               (They should be the same)
C
C     Purpose: Calculate the gamma intermediate.
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION DSRHF(*),GAMMA(*),SCRM(*),XLAM(*),XLAMC(*)
      DIMENSION SCR1(*),SCR2(*),SCR3(*),SCR4(*),SCR5(*),WORK(*)
      DIMENSION XLAMDP(*),XLAMDH(*),XLAMPC(*),XLAMHC(*)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      ISYMJB = MULD2H(ISYMD,ISYMPC)
      ISKILJ = MULD2H(ISYMCM,ISYMOP)
      IF (ISYMPC .NE. ISYMCM) CALL QUIT('Symmetry mismatch in CC_GAM1')
C
      ISYMKC = ISYMAG
C
      DO 100 L = 1,NRHF(ISYML)
C
         KOFF1 = IDSRHF(ISYMAG,ISYML) + NNBST(ISYMAG)*(L - 1) + 1
C
         CALL CCSD_SYMSQ(DSRHF(KOFF1),ISYMAG,SCR1)
C
         DO 110 ISYMG = 1,NSYM
C
            ISYMA = MULD2H(ISYMG,ISYMAG)
            ISYMK = ISYMA
            ISYMC = ISYMG
            ISYMI = ISYMG
C
            NBASA = MAX(NBAS(ISYMA),1)
            NBASG = MAX(NBAS(ISYMG),1)
            NRHFK = MAX(NRHF(ISYMK),1)
            NVIRC = MAX(NVIR(ISYMC),1)
C
            KOFF2 = ILMRHF(ISYMK) + 1
            KOFF3 = IAODIS(ISYMA,ISYMG) + 1
            KOFF4 = IT1AOT(ISYMK,ISYMG) + 1
C
            CALL DGEMM('T','N',NRHF(ISYMK),NBAS(ISYMG),NBAS(ISYMA),
     *                 ONE,XLAMDP(KOFF2),NBASA,SCR1(KOFF3),NBASA,
     *                 ZERO,SCR2(KOFF4),NRHFK)
C
            KOFF5 = ILMVIR(ISYMC) + 1
            KOFF6 = IT1AMT(ISYMK,ISYMC) + 1
C
            CALL DGEMM('N','N',NRHF(ISYMK),NVIR(ISYMC),NBAS(ISYMG),
     *                 ONE,SCR2(KOFF4),NRHFK,XLAMDH(KOFF5),NBASG,
     *                 ZERO,SCR3(KOFF6),NRHFK)
C
            KOFF7 = ILMRHF(ISYMI) + 1
            KOFF8 = IMATIJ(ISYMK,ISYMI) + 1
C
            CALL DGEMM('N','N',NRHF(ISYMK),NRHF(ISYMI),NBAS(ISYMG),
     *                 ONE,SCR2(KOFF4),NRHFK,XLAMDH(KOFF7),NBASG,
     *                 ZERO,SCR4(KOFF8),NRHFK)
C
            IF (IOPT .EQ. 2) THEN
C
C-------------------------------------------------------------
C              Only for calculating Ik,i-bar,l delta integral.
C-------------------------------------------------------------
C
               ISYMA = MULD2H(ISYMG,ISYMAG)
               ISYMK = ISYMA
               ISYMI = MULD2H(ISYMG,ISYMPC)
C
               KOFF7 = IGLMRH(ISYMG,ISYMI) + 1
               KOFF8 = IMATIJ(ISYMK,ISYMI) + 1
C
               CALL DGEMM('N','N',NRHF(ISYMK),NRHF(ISYMI),NBAS(ISYMG),
     *                    ONE,SCR2(KOFF4),NRHFK,XLAMHC(KOFF7),NBASG,
     *                    ZERO,SCR5(KOFF8),NRHFK)
C
            ENDIF
C
  110    CONTINUE
C
         DO 120 ISYMJ = 1,NSYM
C
            ISYMLJ = MULD2H(ISYML,ISYMJ)
            ISYMKI = MULD2H(ISYMLJ,ISKILJ)
            ISYDVI = ISYMD
            ISYDVJ = MULD2H(ISYDVI,ISYMJ)
            ISYMCI = MULD2H(ISYMCM,ISYDVJ)
C
            IF (ISYMKI .GT. ISYMLJ) GOTO 120
C
            KLC = IGLMRH(ISYMD,ISYMJ) + IDEL - IBAS(ISYMD)
            CALL DCOPY(NRHF(ISYMJB),XLAMHC(KLC),NBAS(ISYMD),
     *                 XLAMC,1)
C
            KSCR5 = 1
            KEND1 = KSCR5 + NMATIJ(ISYMKI)
C
            DO 130 J = 1,NRHF(ISYMJ)
C
               DO 140 ISYMI = 1,NSYM
C
                  ISYMC = MULD2H(ISYMI,ISYMCI)
                  ISYMK = MULD2H(ISYMI,ISYMKI)
C
                  NVIRC = MAX(NVIR(ISYMC),1)
                  NRHFK = MAX(NRHF(ISYMK),1)
C
                  KOFF2 = IT1AMT(ISYMK,ISYMC) + 1
                  KOFF3 = IT2BCD(ISYMCI,ISYMJ)
     *                  + NT1AM(ISYMCI)*(J - 1)
     *                  + IT1AM(ISYMC,ISYMI) + 1
                  KOFF4 = KSCR5 + IMATIJ(ISYMK,ISYMI)
C
                  CALL DGEMM('N','N',NRHF(ISYMK),NRHF(ISYMI),
     *                       NVIR(ISYMC),ONE,SCR3(KOFF2),NRHFK,
     *                       SCRM(KOFF3),NVIRC,ZERO,WORK(KOFF4),NRHFK)
C
  140          CONTINUE
C
               IF ( IOPT .EQ. 2 ) THEN 
                  IF (ISYMJ .EQ. ISYMD) THEN
                     CALL DAXPY(NMATIJ(ISYMKI),XLAM(J),SCR5,1,
     *                          WORK(KSCR5),1)
                  ENDIF
                  IF (MULD2H(ISYMJ,ISYMD).EQ.ISYMPC) THEN
                     CALL DAXPY(NMATIJ(ISYMKI),XLAMC(J),SCR4,1,
     *                          WORK(KSCR5),1)
                  ENDIF
               ELSE
                  IF (ISYMJ .EQ. ISYMD) THEN
                     CALL DAXPY(NMATIJ(ISYMKI),XLAM(J),SCR4,1,
     *                          WORK(KSCR5),1)
                  ENDIF
               ENDIF
C
               NLJ = IMATIJ(ISYML,ISYMJ) + NRHF(ISYML)*(J - 1) + L
C
               IF (ISKILJ .EQ. 1) THEN
                  KKILJ = IGAMMA(ISYMKI,ISYMLJ) + NLJ*(NLJ-1)/2
                  DO 150 NKI = 1,NLJ
C
                     KOFF = KSCR5 + NKI - 1
                     NKILJ = KKILJ + NKI
                     GAMMA(NKILJ) = GAMMA(NKILJ) + WORK(KOFF)
C
  150             CONTINUE
               ELSE
                  KKILJ = IGAMMA(ISYMKI,ISYMLJ)
     *                  + NMATIJ(ISYMKI)*(NLJ - 1)
                  DO 160 NKI = 1,NMATIJ(ISYMKI)
C
                     KOFF = KSCR5 + NKI - 1
                     NKILJ = KKILJ + NKI
                     GAMMA(NKILJ) = GAMMA(NKILJ) + WORK(KOFF)
C
  160             CONTINUE
               END IF
C
  130       CONTINUE
  120    CONTINUE
C
         IF (IOPT .EQ. 2) THEN
C
C------------------------------------
C        Extra F-matrix term section.
C------------------------------------
C
         ENDIF
C
  100 CONTINUE
C
      RETURN
      END
C  /* Deck cc_zwvi */
      SUBROUTINE CC_ZWVI(ZINT,CTR2,ISYMC2 ,TINT,ISYTIN,
     *                   WORK,LWORK,IOPT)
C
C     Written by Asger Halkier 26/10 - 1995
C
C     Version: 1.0
C
C     Purpose: To calculate the intermediates entering some of the
C              terms in the 2.1-block.
C
C     IOPT equals 2 if transposition of the occupied indices of
C     TINT intermediate is needed, and 1 if not!
C
C     Ove Christiansen 1-10-1996:
C              general symmetry of ctr2 (isymc2)
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION ZINT(*), CTR2(*), TINT(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
C
C-----------------------------
C     Initialize result array.
C-----------------------------
C
      ISYZIN = MULD2H(ISYMC2,ISYTIN)
C
      CALL DZERO(ZINT,NT2BCD(ISYZIN))
C
C-----------------------------------------------------
C     Transpose occupied indices of TINT if requested.
C-----------------------------------------------------
C
      IF (IOPT .EQ. 2) THEN
C
         CALL CCLT_P21I(TINT,ISYTIN,WORK,LWORK,
     &                  IT2BCD,NT2BCD,IT1AM,NT1AM,NVIR)
C
      ENDIF
C
C------------------------
C     Do the calculation.
C------------------------
C
      DO 100 ISYMDK = 1,NSYM
C
         ISYMEI = MULD2H(ISYMDK,ISYMC2)
         ISYMJ  = MULD2H(ISYMDK,ISYTIN)
C
         KOFF1  = IT2SQ(ISYMEI,ISYMDK) + 1
         KOFF2  = IT2BCD(ISYMDK,ISYMJ) + 1
         KOFF3  = IT2BCD(ISYMEI,ISYMJ) + 1
C
         NTOTEI = MAX(NT1AM(ISYMEI),1)
         NTOTDK = MAX(NT1AM(ISYMDK),1)
C
         CALL DGEMM('N','N',NT1AM(ISYMEI),NRHF(ISYMJ),NT1AM(ISYMDK),
     *              ONE,CTR2(KOFF1),NTOTEI,TINT(KOFF2),NTOTDK,ZERO,
     *              ZINT(KOFF3),NTOTEI)
C
  100 CONTINUE
C
C-------------------------------------------
C     Restore TINT intermediate if necessary
C-------------------------------------------
C
      IF (IOPT .EQ. 2) THEN
C
         CALL CCLT_P21I(TINT,ISYTIN,WORK,LWORK,
     &                  IT2BCD,NT2BCD,IT1AM,NT1AM,NVIR)
C
      ENDIF
C
      RETURN
      END
C  /* Deck cc_ti */
      SUBROUTINE CC_TI(TINT,ISYTIN,T2AM,ISYMT2,XLAMDH,ISYMLH,
     *                 WORK,LWORK,IDEL,ISYMD)
C
C     Written by Asger Halkier 26/10 - 1995
C
C     Version: 1.0
C
C     Purpose: To calculate the T-intermediate entering some of the
C              terms in the 2.1-block.
C
C     Ove Christiansen 30-9-1996: Generalised symmetry of T2AM: ISYMT2
C              Still it is assumed that XLAMDH is total symmetric.
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      INTEGER ISYTIN, ISYMT2, ISYMLH, IDEL, ISYMD
      DIMENSION TINT(*), T2AM(*), XLAMDH(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      CALL QENTER('CC_TI')
C
C-----------------------------
C     Initialize result array.
C-----------------------------
C
      ISYMF  = MULD2H(ISYMD,ISYMLH)
      IF (ISYTIN.NE.MULD2H(ISYMT2,ISYMF)) THEN
        WRITE(LUPRI,*) 'ISYTIN,ISYMT2,ISYMLH,ISYMD:',
     &                  ISYTIN,ISYMT2,ISYMLH,ISYMD
        CALL QUIT('SYMMETRY MISMATCH IN CC_TI')
      END IF
C
      CALL DZERO (TINT,NT2BCD(ISYTIN))
C
C----------------------------------
C     Work space allocation one.
C----------------------------------
C
      KLAHF = 1
      KEND1 = KLAHF + NVIR(ISYMF)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient work space in CC_TI')
      ENDIF
C
C------------------------------------------
C     Copy vector out of lamda hole matrix.
C------------------------------------------
C
      KOFF1  = IGLMVI(ISYMD,ISYMF) + IDEL - IBAS(ISYMD)
C
      CALL DCOPY(NVIR(ISYMF),XLAMDH(KOFF1),NBAS(ISYMD),WORK(KLAHF),1)
C
      DO 100 ISYMDK = 1,NSYM
C
         ISYMFJ = MULD2H(ISYMDK,ISYMT2)
         ISYMJ  = MULD2H(ISYMFJ,ISYMF)
C
         IF (NRHF(ISYMJ) .EQ. 0) GOTO 100
C
C----------------------------------
C        Work space allocation two.
C----------------------------------
C
         KT2SM = KEND1
         KEND2 = KT2SM + NT1AM(ISYMDK)*NVIR(ISYMF)
         LWRK2 = LWORK - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            CALL QUIT('Insufficient work space in CC_TI')
         ENDIF
C
         DO 110 J = 1,NRHF(ISYMJ)
C
C---------------------------------------------
C           Copy submatrix out of packed T2AM.
C---------------------------------------------
C
            DO 120 F = 1,NVIR(ISYMF)
C
               NFJ = IT1AM(ISYMF,ISYMJ) + NVIR(ISYMF)*(J - 1) + F
C
               IF (ISYMDK .EQ. ISYMFJ) THEN
C
                  DO 130 NDK = 1,NT1AM(ISYMDK)
C
                     NDKFJ = IT2AM(ISYMDK,ISYMFJ) + INDEX(NDK,NFJ)
                     NDKF  = KT2SM + NT1AM(ISYMDK)*(F - 1) + NDK - 1
C
                     WORK(NDKF) = T2AM(NDKFJ)
C
  130             CONTINUE
C
               ELSE IF (ISYMDK .LT. ISYMFJ) THEN
C
                  NDKFJ = IT2AM(ISYMDK,ISYMFJ)
     *                  + NT1AM(ISYMDK)*(NFJ - 1) + 1
                  NDKF  = KT2SM + NT1AM(ISYMDK)*(F - 1)
C
                  CALL DCOPY(NT1AM(ISYMDK),T2AM(NDKFJ),1,WORK(NDKF),1)
C
               ELSE IF (ISYMDK .GT. ISYMFJ) THEN
C
                  NDKFJ = IT2AM(ISYMFJ,ISYMDK) + NFJ
                  NDKF  = KT2SM + NT1AM(ISYMDK)*(F - 1)
C
                  CALL DCOPY(NT1AM(ISYMDK),T2AM(NDKFJ),NT1AM(ISYMFJ),
     *                       WORK(NDKF),1)
C
               ENDIF
C
  120       CONTINUE
C
C-----------------------------------------------------
C           Contraction of T2AM-submatrix with XLAMDH.
C-----------------------------------------------------
C
            KOFF2  = IT2BCD(ISYMDK,ISYMJ) + NT1AM(ISYMDK)*(J - 1) + 1
C
            NTOTDK = MAX(NT1AM(ISYMDK),1)
C
            CALL DGEMV('N',NT1AM(ISYMDK),NVIR(ISYMF),ONE,WORK(KT2SM),
     *                 NTOTDK,WORK(KLAHF),1,ZERO,TINT(KOFF2),1)
C
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('CC_TI')
C
      RETURN
      END
C  /* Deck cc_21a */
      SUBROUTINE CC_21A(RHO1,DSRHF,YMAT,ISYMY,YDEN,ISYMYD,
     *                  XLAMDH,XLAMDP,ISYMLP,WORK,
     *                  LWORK,IDEL,ISYMD)
C
C     Written by Asger Halkier 4/10 - 1995.
C
C     Version: 1.0
C
C     Purpose: To calculate the 21A contribution to rho1!
C
C     Ove Christiansen 1-10-1996:
C              Generalization to general symmetry YMAT and YDEN
C              for F-matrix.
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION RHO1(*), DSRHF(*), YMAT(*), YDEN(*), XLAMDH(*),
     *          XLAMDP(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
C
      ISYRES = MULD2H(ISYMYD,ISYMOP)
      IF (MULD2H(ISYMY,ISYMLP).NE.ISYRES) THEN
         WRITE(LUPRI,*) ' Symmetry mismatch in CC_21A'
         CALL QUIT(' Symmetry mismatch in CC_21A')
      ENDIF
      ISYINT = MULD2H(ISYMD,ISYMOP)
C
C================================
C     Calculate the coulomb part.
C================================
C
      ISYMA  = ISYMD
      ISYMI  = MULD2H(ISYMA,ISYRES)
      ISALBE = MULD2H(ISYMI,ISYINT)
C
C----------------------------------------------------------------------
C     Work space allocation 1 & options for contraction with integrals.
C----------------------------------------------------------------------
C
      KINT   = 1
      KVECIS = KINT   + N2BST(ISALBE)
      KENSMA = KVECIS + NRHF(ISYMI)
      KVECIB = KINT   + NRHF(ISYMI)*N2BST(ISALBE)
      KENBIG = KVECIB + NRHF(ISYMI)
      LWRSMA = LWORK  - KENSMA
      LWRBIG = LWORK  - KENBIG
C
      IF (LWRSMA .LT. 0) THEN
         CALL QUIT('Insufficient work space in CC_21A')
      ENDIF
C
      IF (LWRBIG .LE. 0) THEN
C
         DO 110 I = 1,NRHF(ISYMI)
C
            KOFF1 = IDSRHF(ISALBE,ISYMI) + NNBST(ISALBE)*(I - 1) + 1
C
            CALL CCSD_SYMSQ(DSRHF(KOFF1),ISALBE,WORK(KINT))
C
            KOFF2 = KVECIS + I - 1
C
            WORK(KOFF2) = DDOT(N2BST(ISALBE),YDEN,1,
     *                      WORK(KINT),1)
C
  110    CONTINUE
C
      ELSE
C
         DO 120 I = 1,NRHF(ISYMI)
C
            KOFF1 = IDSRHF(ISALBE,ISYMI) + NNBST(ISALBE)*(I - 1) + 1
            KOFF2 = KINT + N2BST(ISALBE)*(I - 1)
C
            CALL CCSD_SYMSQ(DSRHF(KOFF1),ISALBE,WORK(KOFF2))
C
  120    CONTINUE
C
         NTALBE = MAX(N2BST(ISALBE),1)
C
         CALL DGEMV('T',N2BST(ISALBE),NRHF(ISYMI),TWO,WORK(KINT),NTALBE,
     *              YDEN,1,ZERO,WORK(KVECIB),1)
C
         CALL DCOPY(NRHF(ISYMI),WORK(KVECIB),1,WORK(KVECIS),1)
C
      ENDIF
C
C-------------------------------------------------
C     Scale with XLAMDH-element and add to result.
C-------------------------------------------------
C
      DO 130 A = 1,NVIR(ISYMA)
C
         KOFF1 = ILMVIR(ISYMA) + NBAS(ISYMD)*(A - 1)
     *         + IDEL - IBAS(ISYMD)
         KOFF2 = IT1AM(ISYMA,ISYMI) + A
C
         CALL DAXPY(NRHF(ISYMI),XLAMDH(KOFF1),WORK(KVECIS),1,
     *              RHO1(KOFF2),NVIR(ISYMA))
C
  130 CONTINUE
C
C=================================
C     Calculate the exchange part.
C=================================
C
      ISYME  = ISYMD
      ISYMF  = MULD2H(ISYME,ISYMY) 
      ISYMAL = MULD2H(ISYMF,ISYMLP)
      ISYBEI = MULD2H(ISYMAL,ISYINT)
C
      DO 140 ISYMI = 1,NSYM
C
         ISYMBE = MULD2H(ISYMI,ISYBEI)
         ISALBE = MULD2H(ISYMBE,ISYMAL)
         ISYMA  = MULD2H(ISYMI,ISYRES)
C
C--------------------------------
C        Work space allocation 2.
C--------------------------------
C
         KINT  = 1
         KSCR1 = KINT  + N2BST(ISALBE)
         KSCR2 = KSCR1 + NBAS(ISYMAL)*NVIR(ISYMA)
         KSCR3 = KSCR2 + NVIR(ISYMA)*NVIR(ISYMF)
         KEND1 = KSCR3 + NVIR(ISYMA)*NVIR(ISYME)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient work space in CC_21A')
         ENDIF
C
         DO 150 I = 1,NRHF(ISYMI)
C
            KOFFPI = IDSRHF(ISALBE,ISYMI) + NNBST(ISALBE)*(I - 1) + 1
C
            CALL CCSD_SYMSQ(DSRHF(KOFFPI),ISALBE,WORK(KINT))
C
C-------------------------------------------
C           Calculate intermediate matrices.
C-------------------------------------------
C
            KOFFUI = KINT + IAODIS(ISYMAL,ISYMBE)
            KOFF1  = ILMVIR(ISYMA) + 1
C
            NTOTAL = MAX(NBAS(ISYMAL),1)
            NTOTBE = MAX(NBAS(ISYMBE),1)
C
            CALL DGEMM('N','N',NBAS(ISYMAL),NVIR(ISYMA),NBAS(ISYMBE),
     *                 ONE,WORK(KOFFUI),NTOTAL,XLAMDH(KOFF1),NTOTBE,
     *                 ZERO,WORK(KSCR1),NTOTAL)
C
            KOFF2  = IGLMVI(ISYMAL,ISYMF) + 1
C
            NTOTAL = MAX(NBAS(ISYMAL),1)
            NTOTF  = MAX(NVIR(ISYMF),1)
C
            CALL DGEMM('T','N',NVIR(ISYMF),NVIR(ISYMA),NBAS(ISYMAL),
     *                 ONE,XLAMDP(KOFF2),NTOTAL,WORK(KSCR1),NTOTAL,
     *                 ZERO,WORK(KSCR2),NTOTF)
C
            KOFF3 = IMATAB(ISYME,ISYMF) + 1
C
            NTOTE = MAX(NVIR(ISYME),1)
            NTOTF = MAX(NVIR(ISYMF),1)
C
            CALL DGEMM('N','N',NVIR(ISYME),NVIR(ISYMA),NVIR(ISYMF),
     *                 ONE,YMAT(KOFF3),NTOTE,WORK(KSCR2),NTOTF,
     *                 ZERO,WORK(KSCR3),NTOTE)
C
C-------------------------------------------
C           Add contribution to RHO1 vector.
C-------------------------------------------
C
            KOFF4 = ILMVIR(ISYME) + IDEL - IBAS(ISYMD)
            KOFF5 = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + 1
C
            NTOTE = MAX(NVIR(ISYME),1)
C
            CALL DGEMV('T',NVIR(ISYME),NVIR(ISYMA),-ONE,WORK(KSCR3),
     *                 NTOTE,XLAMDH(KOFF4),NBAS(ISYMD),ONE,
     *                 RHO1(KOFF5),1)
C
  150    CONTINUE
C
  140 CONTINUE
C
      RETURN
      END
C  /* Deck cc_yd */
      SUBROUTINE CC_YD(YDEN,YMAT,ISYMY,XLAMDH,XLAMDP,ISYMLP,
     *                 WORK,LWORK)
C
C     Written by Asger Halkier 8/12 - 1995.
C
C     Version: 1.0
C
C     Purpose: To transform the Y-matrix to AO-basis!
C
C     Ove Christiansen 1-10-1996:
C              General symmetry of YMAT (ISYMY) and
C              XLAMDP (ISYMLP)
C              XLAMDH is assumed to have the symmetry ISYMOP
C              (it is a hole-virtuel index.)
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION YMAT(*), YDEN(*), XLAMDH(*), XLAMDP(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
C
      ISYMYD = MULD2H(ISYMLP,ISYMY)
C
C------------------------------------
C     Transform Y-matrix to AO-basis.
C------------------------------------
C
      DO 100 ISYMAL = 1,NSYM
C
         ISYMF  = MULD2H(ISYMAL,ISYMLP)
         ISYMBE = MULD2H(ISYMAL,ISYMYD)
         ISYME  = MULD2H(ISYMBE,ISYMOP)
C
         LWRKCH = LWORK - NBAS(ISYMBE)*NVIR(ISYMF)
C
         IF (LWRKCH .LT. 0) THEN
            CALL QUIT('Insufficient work space in CC_21A')
         ENDIF
C
         KOFF1 = ILMVIR(ISYME) + 1
         KOFF2 = IMATAB(ISYME,ISYMF) + 1
C
         NTOTBE = MAX(NBAS(ISYMBE),1)
         NTOTE  = MAX(NVIR(ISYME),1)
C
         CALL DGEMM('N','N',NBAS(ISYMBE),NVIR(ISYMF),NVIR(ISYME),
     *              ONE,XLAMDH(KOFF1),NTOTBE,YMAT(KOFF2),NTOTE,
     *              ZERO,WORK,NTOTBE)
C
         KOFF1 = IGLMVI(ISYMAL,ISYMF) + 1
         KOFF2 = IAODIS(ISYMAL,ISYMBE) + 1
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTBE = MAX(NBAS(ISYMBE),1)
C
         CALL DGEMM('N','T',NBAS(ISYMAL),NBAS(ISYMBE),NVIR(ISYMF),
     *              ONE,XLAMDP(KOFF1),NTOTAL,WORK,NTOTBE,
     *              ZERO,YDEN(KOFF2),NTOTAL)
C
  100 CONTINUE
C
      RETURN
      END
C  /* Deck cc_21h */
      SUBROUTINE CC_21H(RHO1,ISYRHO,VINT,WINT,ZINT,ISYVWZ,X3OINT,
     *                  ISYINT,WORK,LWORK,ISYMD)
C
C     Written by Asger Halkier 30/10 - 1995
C
C     Version: 1.0
C
C     Purpose: To calculate the 21H contribution to RHO1!
C
C     Ove Christiansen 2-10-1996:
C              Generalisation to general symmetry of intermediates
C              as well as general symmetry of integrals.
C              ISYVWZ is symmetry of V,W, and Z intermediates.
C              ISYINT is symmetry of integrals.
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION RHO1(*), VINT(*), WINT(*), ZINT(*), X3OINT(*),
     &          WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
C
      ISYJLI = MULD2H(ISYINT,ISYMD)
      ISYRES = MULD2H(ISYJLI,ISYVWZ)
      IF (ISYRES .NE. ISYRHO) THEN
         CALL QUIT('Symmetry mismatch in CC_21H')
      ENDIF
C
C-------------------------------
C     Work space allocation one.
C-------------------------------
C
      KVZINT = 1
      KEND1  = KVZINT + NT2BCD(ISYVWZ)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient work space in CC_21H')
      ENDIF
C
C---------------------------------------------------
C     Resort and add together V- and Z-intermediate.
C---------------------------------------------------
C
      CALL CCLT_S21I(WORK(KVZINT),VINT,ZINT,ISYVWZ,ONE)
C
      DO 100 ISYMA = 1,NSYM
C
         ISYMI  = MULD2H(ISYMA,ISYRES)
         ISYMJL = MULD2H(ISYMI,ISYJLI)
C
C--------------------------------------
C        Calculate the VZ-contribution.
C--------------------------------------
C
         KOFF1  = KVZINT + IT2AIJ(ISYMA,ISYMJL)
         KOFF2  = IMAIJK(ISYMJL,ISYMI) + 1
         KOFF3  = IT1AM(ISYMA,ISYMI) + 1
C
         NTOTA  = MAX(NVIR(ISYMA),1)
         NTOTJL = MAX(NMATIJ(ISYMJL),1)
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NMATIJ(ISYMJL),
     &              ONE,WORK(KOFF1),NTOTA,X3OINT(KOFF2),NTOTJL,ONE,
     &              RHO1(KOFF3),NTOTA)
C
  100 CONTINUE
C
C-------------------------------
C     Work space allocation two.
C-------------------------------
C
      KWZINT = 1
      KI3OTR = KWZINT + NT2BCD(ISYVWZ)
      KEND2  = KI3OTR + NMAIJK(ISYJLI)
      LWRK2  = LWORK  - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         CALL QUIT('Insufficient work space in CC_21H')
      ENDIF
C
C---------------------------------------------------------
C     Prepare intermediates and integrals for contraction.
C---------------------------------------------------------
C
      CALL CCLT_S21I(WORK(KWZINT),WINT,ZINT,ISYVWZ,-TWO)
C
      CALL CCLT_PI3O(WORK(KI3OTR),X3OINT,ISYJLI)
C
      DO 110 ISYMA = 1,NSYM
C
         ISYMI  = MULD2H(ISYMA,ISYRES)
         ISYMJL = MULD2H(ISYMI,ISYJLI)
C
C----------------------------------
C        Calculate WZ-contribution.
C----------------------------------
C
         KOFF1  = KWZINT + IT2AIJ(ISYMA,ISYMJL)
         KOFF2  = KI3OTR + IMAIJK(ISYMJL,ISYMI)
         KOFF3  = IT1AM(ISYMA,ISYMI) + 1
C
         NTOTA  = MAX(NVIR(ISYMA),1)
         NTOTJL = MAX(NMATIJ(ISYMJL),1)
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NMATIJ(ISYMJL),
     &              ONE,WORK(KOFF1),NTOTA,WORK(KOFF2),NTOTJL,ONE,
     &              RHO1(KOFF3),NTOTA)
C
  110 CONTINUE
C
      RETURN
      END
C  /* Deck cc_21g */
      SUBROUTINE CC_21G(RHO1,XMINT,ISYMM,XLAMDH,WORK,LWORK,
     *                  ISYINT,LUO3,O3FIL)
C
C     Written by Asger Halkier 30/10 - 1995
C
C     Version: 1.0
C
C     Purpose: To calculate the 21G contribution to RHO1.
C
C     Ove Christiansen 3-10-1996:
C        General symmetries for F-matrix
C        isyint is symmetry of integrals.
C
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION RHO1(*), XMINT(*), XLAMDH(*), WORK(LWORK)
      CHARACTER O3FIL*(*)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
C
      ISYRES = MULD2H(ISYMM,ISYINT)
C
      DO 100 ISYMA = 1,NSYM
C
         ISYMD  = ISYMA
         ISYMI  = MULD2H(ISYMA,ISYRES)
         ISYJKL = MULD2H(ISYMD,ISYINT)
C
         IF (NBAS(ISYMD) .EQ. 0) GOTO 100
C
C----------------------------------
C        Work space allocation one.
C----------------------------------
C
         KMOINT = 1
         KAOINT = KMOINT + NVIR(ISYMA)*NMAIJK(ISYJKL)
         KEND1  = KAOINT + NBAS(ISYMD)*NMAIJK(ISYJKL)
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient work space in CC_21G')
         ENDIF
C
C-------------------------------------------
C        Read integrals (jk|ldel) from disc.
C-------------------------------------------
C
         NTOT = NBAS(ISYMD)*NMAIJK(ISYJKL)
         IOFF = I3ODEL(ISYJKL,ISYMD) + 1
C
         CALL GETWA2(LUO3,O3FIL,WORK(KAOINT),IOFF,NTOT)
C
C-----------------------------------------------------
C        Transform AO integral index to virtual space.
C-----------------------------------------------------
C
         KOFF1  = ILMVIR(ISYMA) + 1
C
         NTOJKL = MAX(NMAIJK(ISYJKL),1)
         NTOTD  = MAX(NBAS(ISYMD),1)
C
         CALL DGEMM('N','N',NMAIJK(ISYJKL),NVIR(ISYMA),NBAS(ISYMD),
     &              ONE,WORK(KAOINT),NTOJKL,XLAMDH(KOFF1),NTOTD,ZERO,
     &              WORK(KMOINT),NTOJKL)
C
C------------------------------------------------------------
C        Contraction with M-intermediate & storage in result.
C------------------------------------------------------------
C
         KOFF2  = I3ORHF(ISYJKL,ISYMI) + 1
         KOFF3  = IT1AM(ISYMA,ISYMI) + 1
C
         NTOJKL = MAX(NMAIJK(ISYJKL),1)
         NTOTA  = MAX(NVIR(ISYMA),1)
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NMAIJK(ISYJKL),
     &              ONE,WORK(KMOINT),NTOJKL,XMINT(KOFF2),NTOJKL,ONE,
     &              RHO1(KOFF3),NTOTA)
C
  100 CONTINUE
C
      RETURN
      END
