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
C  /* Deck cc2_den */
      SUBROUTINE CC2_DEN(ETAAI,ZKDIA,WORK,LWORK,IOPT)
C
C     Written by Asger Halkier & Sonia Coriani 12/1 - 2000
C
C     Version: 1.0
C
C     PURPOSE:
C
C     IOPT=1: To calculate the pp, ab, & ij parts of kappa-bar-0
C             that do not need the solution of any coupled equations
C             in the case of (relaxed) CC2 calculations.
C             ZKDIA holds all the blocks pq in the following order:
C             ij, ab, ai, ia; and these are stored full blocks after
C             each other.
C
C     IOPT=2: To calculate the complete right-hand-side for the
C             kappa-bar-0 ai block (eta_ai) - this includes both the
C             "direct" terms calculated from densities and t1-transformed
C             integrals, as well as the indirect terms coming from
C             kappa-bar-0_ab and kappa-bar-0_ij. The result is stored in
C             ETAAI.
C
C             IMPORTANT NOTE: NO FROZEN CORE SO FAR!!!
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <ccinftap.h>
#include <maxash.h>
#include <maxorb.h>
#include <mxcent.h>
#include <aovec.h>
#include <iratdef.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION ETAAI(*), ZKDIA(*), WORK(LWORK)
      DIMENSION INDEXA(MXCORB) 
#include <ccorb.h>
CCN !#include <infind.h> not compatible with R12!
#include <ccisao.h>
#include <r12int.h>
#include <blocks.h>
#include <ccfield.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <distcl.h>
#include <cbieri.h>
#include <eribuf.h>
#include <cclr.h>
#include <ccfro.h>
C
      CHARACTER MODEL*10
      LOGICAL LOCDBG
      PARAMETER (LOCDBG=.FALSE.)
C
      CALL QENTER('CC2_DEN')
C
      TIMETO = SECOND()
      TIMDEN = ZERO
      TIMHE2 = ZERO
      TIMRES = ZERO
      TIRDAO = ZERO
C
      IF ((IPRINT .GT. 3) .AND. (IOPT .EQ. 1)) THEN
         CALL HEADER('Calculating diagonal blocks of zeta-kappa-0',-1)
      ENDIF
      IF ((IPRINT .GT. 3) .AND. (IOPT .EQ. 2)) THEN
         CALL HEADER('Calculating right-hand-side for kappa-bar_ai',-1)
      ENDIF
C
C------------------------------------------------------------------
C     Both t-vectors and tbar-vectors (zeta) are totally symmetric.
C------------------------------------------------------------------
C
      ISYMTR = 1
      ISYMOP = 1
C
C-------------------------------
C     Work space allocation one.
C-------------------------------
C
      KCMO   = 1
      KETAIJ = KCMO   + NLAMDS
      KETAAB = KETAIJ + NMATIJ(1)
      KT2AM  = KETAAB + NMATAB(1)
      KZ2AM  = KT2AM  + NT2AMX
      KLAMDP = KZ2AM  + NT2SQ(1)
      KLAMDH = KLAMDP + NLAMDT
      KT1AM  = KLAMDH + NLAMDT
      KZ1AM  = KT1AM  + NT1AMX
      KXMAT  = KZ1AM  + NT1AMX
      KYMAT  = KXMAT  + NMATIJ(1)
      KEND1  = KYMAT  + NMATAB(1)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT(
     *         'Insufficient memory for initial allocation in CC2_DEN')
      ENDIF
C
C----------------------------------------------------------
C     Read MO-coefficients from interface file and reorder.
C----------------------------------------------------------
C
      LUSIFC = -1
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',
     *            IDUMMY,.FALSE.)
      REWIND LUSIFC
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KCMO+I-1), I=1,NLAMDS)
      CALL GPCLOSE (LUSIFC,'KEEP')
C
C
      CALL CMO_REORDER(WORK(KCMO),WORK(KEND1),LWRK1)
C
C----------------------------------------
C     Read zero'th order zeta amplitudes.
C----------------------------------------
C
      IOPTRW   = 3
      CALL CC_RDRSP('L0',0,1,IOPTRW,MODEL,WORK(KZ1AM),WORK(KZ2AM))
C
C--------------------------------
C     Square up zeta2 amplitudes.
C--------------------------------
C
      CALL DCOPY(NT2AMX,WORK(KZ2AM),1,WORK(KT2AM),1)
      CALL CC_T2SQ(WORK(KT2AM),WORK(KZ2AM),1)
C
C-------------------------------------------
C     Read zero'th order cluster amplitudes.
C-------------------------------------------
C
      IOPTRW = 3
      CALL CC_RDRSP('R0',0,1,IOPTRW,MODEL,WORK(KT1AM),WORK(KT2AM))
C
C----------------------------------
C     Calculate the lambda matrices.
C----------------------------------
C
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),WORK(KEND1),
     *            LWRK1)
C
C--------------------------------------------------------
C     Calculate X-intermediate of tbar- and t-amplitudes.
C--------------------------------------------------------
C
      IF (IOPT .EQ. 1) THEN
         CALL CC_XI(WORK(KXMAT),WORK(KZ2AM),ISYMTR,WORK(KT2AM),ISYMOP,
     *                WORK(KEND1),LWRK1)
C
C--------------------------------------------------------
C     Calculate Y-intermediate of tbar- and t-amplitudes.
C--------------------------------------------------------
C
         CALL CC_YI(WORK(KYMAT),WORK(KZ2AM),ISYMTR,WORK(KT2AM),ISYMOP,
     *              WORK(KEND1),LWRK1)
C
C--------------------------------------------------------------------------
C     Calculate the diagonal elements ZK0(ii) = -X(ii) and ZK0(aa) = Y(aa).
C--------------------------------------------------------------------------
C
         DO 105 ISYMI = 1,NSYM
            DO 115 I = 1,NRHF(ISYMI)
C
               NII = IMATIJ(ISYMI,ISYMI) + NRHF(ISYMI)*(I - 1) + I
C
               ZKDIA(NII) = -WORK(KXMAT + NII - 1)
C
  115       CONTINUE
  105    CONTINUE
C
         DO 125 ISYMA = 1,NSYM
            DO 135 A = 1,NVIR(ISYMA)
C
               NAA = IMATAB(ISYMA,ISYMA) + NVIR(ISYMA)*(A - 1) + A
C
               ZKDIA(NMATIJ(1) + NAA) = WORK(KYMAT + NAA - 1)
C
  135       CONTINUE
  125    CONTINUE
      ENDIF
C
C-----------------------------------------------
C     Set up 2C-E of cluster amplitudes and save
C     in KT2AM, as we only need T(2c-e) below.
C-----------------------------------------------
C
      ISYOPE = 1
      IOPTTCME = 1
      CALL CCSD_TCMEPK(WORK(KT2AM),1.0D0,ISYOPE,IOPTTCME)
      KT2AMT = KT2AM                  !for safety
C
C-------------------------------
C     Work space allocation one.
C     Note that D(ai) = ZETA(ai)
C     and both D(ia) and h(ia) 
C     are stored transposed!
C-------------------------------
C
      KONEAI = KZ1AM
      KONEAB = KONEAI + NT1AMX
      KONEIJ = KONEAB + NMATAB(1)
      KONEIA = KONEIJ + NMATIJ(1)
      KONINT = KONEIA + NT1AMX
      KINTAI = KONINT + N2BST(ISYMOP)
      KINTAB = KINTAI + NT1AMX
      KINTIJ = KINTAB + NMATAB(1)
      KINTIA = KINTIJ + NMATIJ(1)
      KINABT = KINTIA + NT1AMX
      KINIJT = KINABT + NMATAB(1)
      KD1ABT = KINIJT + NMATIJ(1)
      KD1IJT = KD1ABT + NMATAB(1)
      KEND1  = KD1IJT + NMATIJ(1)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient memory for allocation 1 CC2_DEN')
      ENDIF
C
C
C------------------------------------------------------
C     Initialize remaining one electron density arrays.
C------------------------------------------------------
C
      CALL DZERO(WORK(KONEAB),NMATAB(1))
      CALL DZERO(WORK(KONEIJ),NMATIJ(1))
      CALL DZERO(WORK(KONEIA),NT1AMX)
C
C--------------------------------------------------------
C     Construct remaining blocks of one electron density.
C--------------------------------------------------------
C
      CALL DZERO(WORK(KINTIJ),NMATIJ(1))
      CALL DIJGEN(WORK(KONEIJ),WORK(KINTIJ))
      CALL DIAGEN(WORK(KONEIA),WORK(KT2AMT),WORK(KONEAI))
C
      IF (LOCDBG) THEN
         DIJNO = DDOT(NMATIJ(1),WORK(KONEIJ),1,WORK(KONEIJ),1)
         DAINO = DDOT(NT1AMX,WORK(KONEAI),1,WORK(KONEAI),1)
         DIANO = DDOT(NT1AMX,WORK(KONEIA),1,WORK(KONEIA),1)
         DABNO = DDOT(NMATAB(1),WORK(KONEAB),1,WORK(KONEAB),1)
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,*) 'Norm of one electron densities in MO-basis:'
         WRITE(LUPRI,*) DIJNO, DAINO, DIANO, DABNO
      ENDIF
C
C---------------------------------
C     Read one-electron integrals.
C---------------------------------
C
      CALL CCRHS_ONEAO(WORK(KONINT),WORK(KEND1),LWRK1)
C
C--------------------------------------------------
C     Transform one electron integrals to MO-basis.
C--------------------------------------------------
C
      ISYM = 1
      CALL CCDINTMO(WORK(KINTIJ),WORK(KINTIA),WORK(KINTAB),
     *              WORK(KINTAI),WORK(KONINT),WORK(KLAMDP),
     *              WORK(KLAMDH),WORK(KEND1),LWRK1,ISYM)
C
      IF (LOCDBG) THEN
C
         HIJNO = DDOT(NMATIJ(1),WORK(KINTIJ),1,WORK(KINTIJ),1)
         HAINO = DDOT(NT1AMX,WORK(KINTAI),1,WORK(KINTAI),1)
         HIANO = DDOT(NT1AMX,WORK(KINTIA),1,WORK(KINTIA),1)
         HABNO = DDOT(NMATAB(1),WORK(KINTAB),1,WORK(KINTAB),1)
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,*) 'Norm of one electron integrals in MO-basis:'
         WRITE(LUPRI,*) HIJNO, HAINO, HIANO, HABNO
C
      ENDIF
C
C-----------------------------------------------
C     Set up transposed integrals and densities.
C-----------------------------------------------
C
      CALL DCOPY(NMATAB(1),WORK(KINTAB),1,WORK(KINABT),1)
      CALL DCOPY(NMATIJ(1),WORK(KINTIJ),1,WORK(KINIJT),1)
      CALL DCOPY(NMATAB(1),WORK(KONEAB),1,WORK(KD1ABT),1)
      CALL DCOPY(NMATIJ(1),WORK(KONEIJ),1,WORK(KD1IJT),1)
C
      CALL CC_EITR(WORK(KINABT),WORK(KINIJT),WORK(KEND1),LWRK1,1)
      CALL CC_EITR(WORK(KD1ABT),WORK(KD1IJT),WORK(KEND1),LWRK1,1)
C
      IF (IOPT .EQ. 1) THEN
C
C---------------------------------------------------------------------
C     Calculate direct one electron contribution to eta-ij and eta-ab.
C---------------------------------------------------------------------
C
         CALL DZERO(WORK(KETAIJ),NMATIJ(1))
         CALL DZERO(WORK(KETAAB),NMATAB(1))
C
         ISYM = 1
         CALL CC2_ETIJ(WORK(KETAIJ),WORK(KINTIJ),WORK(KINTAI),
     *                 WORK(KINTIA),WORK(KINTAB),WORK(KONEIJ),
     *                 WORK(KONEAI),WORK(KONEIA),WORK(KONEAB),
     *                 WORK(KT1AM),WORK(KEND1),LWRK1,ISYM)
C
         CALL CC2_ETAB(WORK(KETAAB),WORK(KINTIJ),WORK(KINTAI),
     *                 WORK(KINTIA),WORK(KINTAB),WORK(KONEIJ),
     *                 WORK(KONEAI),WORK(KONEIA),WORK(KONEAB),
     *                 WORK(KT1AM),WORK(KEND1),LWRK1,ISYM)
C
      ELSE IF (IOPT .EQ. 2) THEN
C
C----------------------------------------------------------
C     Calculate direct one electron contribution to eta-ai.
C----------------------------------------------------------
C
         ISYM = 1
         CALL CCDENZK0(ETAAI,WORK(KINTIJ),WORK(KINTAI),WORK(KINTIA),
     *                 WORK(KINTAB),WORK(KONEIJ),WORK(KONEAI),
     *                 WORK(KONEIA),WORK(KONEAB),WORK(KINIJT),
     *                 WORK(KINABT),WORK(KD1IJT),WORK(KD1ABT),
     *                 WORK(KT1AM),WORK(KEND1),LWRK1,ISYM)
C
      ENDIF
C
C     TIMONE = SECOND() - TIMONE
C
      KEND1 = KONINT
      LWRK1 = LWORK - KEND1
C
C---------------------------------------------------------------
C     Start the loop for the 2-electron integrals and densities.
C---------------------------------------------------------------
C
      KENDS2 = KEND1
      LWRKS2 = LWRK1
C
      IF (DIRECT) THEN
         IF (HERDIR) THEN
            CALL HERDI1(WORK(KEND1),LWRK1,IPRERI)
         ELSE
            KCCFB1 = KEND1
            KINDXB = KCCFB1 + MXPRIM*MXCONT
            KEND1  = KINDXB + (8*MXSHEL*MXCONT + 1)/IRAT
            LWRK1  = LWORK  - KEND1
            CALL ERIDI1(KODCL1,KODCL2,KODBC1,KODBC2,KRDBC1,KRDBC2,
     *                  KODPP1,KODPP2,KRDPP1,KRDPP2,
     *                  KFREE,LFREE,KEND1,WORK(KCCFB1),WORK(KINDXB),
     *                  WORK(KEND1),LWRK1,IPRERI)
            KEND1 = KFREE
            LWRK1 = LFREE
         ENDIF
         NTOSYM = 1
      ELSE
         NTOSYM = NSYM
      ENDIF
C
      KENDSV = KEND1
      LWRKSV = LWRK1
C
      ICDEL1 = 0
      DO 100 ISYMD1 = 1,NTOSYM
C
         IF (DIRECT) THEN
            IF (HERDIR) THEN
               NTOT = MAXSHL
            ELSE
               NTOT = MXCALL
            ENDIF
         ELSE
            NTOT = NBAS(ISYMD1)
         ENDIF
C
         DO 110 ILLL = 1,NTOT
C
C---------------------------------------------
C           If direct calculate the integrals.
C---------------------------------------------
C
            IF (DIRECT) THEN
C
               KEND1 = KENDSV
               LWRK1 = LWRKSV
C
               DTIME  = SECOND()
               IF (HERDIR) THEN
                  CALL HERDI2(WORK(KEND1),LWRK1,INDEXA,ILLL,NUMDIS,
     &                        IPRINT)
               ELSE
                  CALL ERIDI2(ILLL,INDEXA,NUMDIS,0,0,
     *                        WORK(KODCL1),WORK(KODCL2),
     *                        WORK(KODBC1),WORK(KODBC2),
     *                        WORK(KRDBC1),WORK(KRDBC2),
     *                        WORK(KODPP1),WORK(KODPP2),
     *                        WORK(KRDPP1),WORK(KRDPP2),
     *                        WORK(KCCFB1),WORK(KINDXB),
     *                        WORK(KEND1), LWRK1,IPRERI)
               ENDIF
               DTIME  = SECOND() - DTIME
               TIMHE2 = TIMHE2   + DTIME
C
               KRECNR = KEND1
               KEND1  = KRECNR + (NBUFX(0) - 1)/IRAT + 1
               LWRK1  = LWORK  - KEND1
               IF (LWRK1 .LT. 0) THEN
                  CALL QUIT('Insufficient core in CC2_DEN')
               END IF
C
            ELSE
               NUMDIS = 1
            ENDIF
C
C-----------------------------------------------------
C           Loop over number of distributions in disk.
C-----------------------------------------------------
C
            DO 120 IDEL2 = 1,NUMDIS
C
               IF (DIRECT) THEN
                  IDEL  = INDEXA(IDEL2)
                  ISYMD = ISAO(IDEL)
                  IF (NOAUXB) THEN
                    IDUM = 1
                    CALL IJKAUX(IDEL,IDUM,IDUM,IDUM)
                  END IF
               ELSE
                  IDEL  = IBAS(ISYMD1) + ILLL
                  ISYMD = ISYMD1
               ENDIF
C
C----------------------------------------
C              Work space allocation two.
C----------------------------------------
C
               ISYDEN = ISYMD
C
               KD2IJG = KEND1
               KD2AIG = KD2IJG + ND2IJG(ISYDEN)
               KD2IAG = KD2AIG + ND2AIG(ISYDEN)
               KD2ABG = KD2IAG + ND2AIG(ISYDEN)
               KEND2  = KD2ABG + ND2ABG(ISYDEN)
               LWRK2  = LWORK  - KEND2
C
               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:',KEND2
                  CALL QUIT(
     *                'Insufficient space for allocation 2 in CC2_DEN')
               ENDIF
C
C-------------------------------------------------------
C              Initialize 4 two electron density arrays.
C-------------------------------------------------------
C
               CALL DZERO(WORK(KD2IJG),ND2IJG(ISYDEN))
               CALL DZERO(WORK(KD2AIG),ND2AIG(ISYDEN))
               CALL DZERO(WORK(KD2IAG),ND2AIG(ISYDEN))
               CALL DZERO(WORK(KD2ABG),ND2ABG(ISYDEN))
C
C-------------------------------------------------------------------
C              Calculate the two electron density d(pq,gamma;delta).
C-------------------------------------------------------------------
C
               AUTIME = SECOND()
C
               CALL CC_DEN2(WORK(KD2IJG),WORK(KD2AIG),WORK(KD2IAG),
     *                      WORK(KD2ABG),WORK(KZ2AM),WORK(KT2AM),
     *                      WORK(KT2AMT),WORK(KEND2),WORK(KEND2),
     *                      WORK(KEND2),WORK(KONEAB),WORK(KONEAI),
     *                      WORK(KONEIA),WORK(KEND2),WORK(KLAMDH),1,
     *                      WORK(KLAMDP),1,WORK(KEND2),LWRK2,IDEL,ISYMD)
C
               AUTIME = SECOND() - AUTIME
               TIMDEN = TIMDEN + AUTIME
C
C------------------------------------------
C              Work space allocation three.
C------------------------------------------
C
               ISYDIS = MULD2H(ISYMD,ISYMOP)
C
               KXINT  = KEND2
               KEND3  = KXINT  + NDISAO(ISYDIS)
               LWRK3  = LWORK  - KEND3
C
               IF (LWRK3 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:',KEND3
                  CALL QUIT(
     *         'Insufficient space for allocation 3 in CC2_DEN')
               ENDIF
C
C--------------------------------------------
C              Read AO integral distribution.
C--------------------------------------------
C
               AUTIME = SECOND()
               CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND3),LWRK3,
     *                     WORK(KRECNR),DIRECT)
               AUTIME = SECOND() - AUTIME
               TIRDAO = TIRDAO + AUTIME
C
               IF (IOPT .EQ. 2) THEN
C
C--------------------------------------------------------------
C              Calculate correction terms to eta_ai originating
C              from kbar_ij and kbar_ab.
C--------------------------------------------------------------
C
                  KDSRHF = KEND3
                  K3OINT = KDSRHF + NDSRHF(ISYMD)
                  KEND3  = K3OINT + NMAIJK(ISYDIS)
                  LWRK3  = LWORK  - KEND3
C
                  IF (LWRK3 .LT. 0) THEN
                     WRITE(LUPRI,*) 'Need : ',KEND3,'Available : ',LWORK
                     CALL QUIT(
     *                   'Insufficient core for integrals in CC2_DEN')
                  ENDIF
C
C--------------------------------------------------------------------
C              Transform one index in the integral batch to occupied.
C--------------------------------------------------------------------
C
                  CALL CCTRBT(WORK(KXINT),WORK(KDSRHF),WORK(KCMO),
     *                        ISYMOP,WORK(KEND3),LWRK3,ISYDIS)
C
C-------------------------------------------------------------------
C              Calculate integral batch with three occupied indices.
C              Since LUDUMLOCAL < 0 the integrals are not written 
C              to disk.
C-------------------------------------------------------------------
C
                  LUDUMLOCAL = -10
                  CALL CC_INT3O(WORK(K3OINT),WORK(KDSRHF),WORK(KCMO),
     *                         ISYMOP,WORK(KCMO),WORK(KEND3),LWRK3,
     *                         IDEL,ISYMD,LUDUMLOCAL,'DUMMY')
C
C----------------------------------------------------
C              Calculate the actual correction terms.
C----------------------------------------------------
C
                  CALL MP2_YTV(ETAAI,ZKDIA(NMATIJ(1)+1),WORK(KDSRHF),
     *                      WORK(KCMO),WORK(KEND3),LWRK3,IDEL,ISYMD)
C
                  CALL MP2_NXY(ETAAI,ZKDIA(1),ZKDIA(NMATIJ(1)+1),
     *                         WORK(K3OINT),WORK(KDSRHF),WORK(KCMO),
     *                         WORK(KEND3),LWRK3,IDEL,ISYMD)
C
                  CALL MP2_XTO(ETAAI,ZKDIA(1),WORK(K3OINT),
     *                         WORK(KCMO),WORK(KEND3),LWRK3,IDEL,ISYMD)
C
               ENDIF
C
C------------------------------------------------------
C              Start loop over second AO-index (gamma).
C------------------------------------------------------
C
               AUTIME = SECOND()
C
               DO 130 ISYMG = 1,NSYM
C
                  ISYMPQ = MULD2H(ISYMG,ISYDEN)
C
                  DO 140 G = 1,NBAS(ISYMG)
C
                     KD2GIJ = KD2IJG + ID2IJG(ISYMPQ,ISYMG)
     *                      + NMATIJ(ISYMPQ)*(G - 1) 
                     KD2GAI = KD2AIG + ID2AIG(ISYMPQ,ISYMG)
     *                      + NT1AM(ISYMPQ)*(G - 1)
                     KD2GAB = KD2ABG + ID2ABG(ISYMPQ,ISYMG)
     *                      + NMATAB(ISYMPQ)*(G - 1)
                     KD2GIA = KD2IAG + ID2AIG(ISYMPQ,ISYMG)
     *                      + NT1AM(ISYMPQ)*(G - 1)
C
C-----------------------------------------------
C                    Work space allocation four.
C-----------------------------------------------
C
                     KINTAO = KEND3
                     KD2AOB = KINTAO + N2BST(ISYMPQ)
                     KEND4  = KD2AOB + N2BST(ISYMPQ)
                     LWRK4  = LWORK  - KEND4
C
                     IF (LWRK4 .LT. 0) THEN
                        WRITE(LUPRI,*) 'Available:', LWORK
                        WRITE(LUPRI,*) 'Needed:', KEND4
                        CALL QUIT('Insufficient work space in CC2_DEN')
                     ENDIF
C
                     CALL DZERO(WORK(KD2AOB),N2BST(ISYMPQ))
C
C-------------------------------------------------------
C                    Square up AO-integral distribution.
C-------------------------------------------------------
C
                     KOFFIN = KXINT + IDSAOG(ISYMG,ISYDIS) 
     *                      + NNBST(ISYMPQ)*(G - 1) 
C
                     CALL CCSD_SYMSQ(WORK(KOFFIN),ISYMPQ,
     *                               WORK(KINTAO))
C
C-----------------------------------------------
C                    Work space allocation five.
C-----------------------------------------------
C
                     KIJINT = KEND4
                     KAIINT = KIJINT + NMATIJ(ISYMPQ)
                     KIAINT = KAIINT + NT1AM(ISYMPQ)
                     KABINT = KIAINT + NT1AM(ISYMPQ)
                     KABTIN = KABINT + NMATAB(ISYMPQ)
                     KIJTIN = KABTIN + NMATAB(ISYMPQ)
                     KD2TAB = KIJTIN + NMATIJ(ISYMPQ)
                     KD2TIJ = KD2TAB + NMATAB(ISYMPQ)
                     KEND5  = KD2TIJ + NMATIJ(ISYMPQ)
                     LWRK5  = LWORK  - KEND5
C
                     IF (LWRK5 .LT. 0) THEN
                        WRITE(LUPRI,*) 'Available:', LWORK
                        WRITE(LUPRI,*) 'Needed:', KEND5
                        CALL QUIT('Insufficient work space in CC2_DEN')
                     ENDIF
C
C-------------------------------------------------------------
C                    Transform 2 integral indices to MO-basis.
C-------------------------------------------------------------
C
                     ISYM = ISYMPQ
                     CALL CCDINTMO(WORK(KIJINT),WORK(KIAINT),
     *                             WORK(KABINT),WORK(KAIINT),
     *                             WORK(KINTAO),WORK(KLAMDP),
     *                             WORK(KLAMDH),WORK(KEND5),
     *                             LWRK5,ISYM)
C
C--------------------------------------------------------------
C                    Set up transposed integrals and densities.
C--------------------------------------------------------------
C
                     CALL DCOPY(NMATAB(ISYMPQ),WORK(KABINT),1,
     *                          WORK(KABTIN),1)
                     CALL DCOPY(NMATIJ(ISYMPQ),WORK(KIJINT),1,
     *                          WORK(KIJTIN),1)
                     CALL DCOPY(NMATAB(ISYMPQ),WORK(KD2GAB),1,
     *                          WORK(KD2TAB),1)
                     CALL DCOPY(NMATIJ(ISYMPQ),WORK(KD2GIJ),1,
     *                          WORK(KD2TIJ),1)
C
                     CALL CC_EITR(WORK(KABTIN),WORK(KIJTIN),
     *                            WORK(KEND5),LWRK5,ISYMPQ)
                     CALL CC_EITR(WORK(KD2TAB),WORK(KD2TIJ),
     *                            WORK(KEND5),LWRK5,ISYMPQ)
C
                     IF (IOPT .EQ. 1) THEN
C
                        ISYM = ISYMPQ
                        CALL CC2_ETIJ(WORK(KETAIJ),WORK(KIJINT),
     *                                WORK(KAIINT),WORK(KIAINT),
     *                                WORK(KABINT),WORK(KD2GIJ),
     *                                WORK(KD2GAI),WORK(KD2GIA),
     *                                WORK(KD2GAB),WORK(KT1AM),
     *                                WORK(KEND5),LWRK5,ISYM)
C
                        CALL CC2_ETAB(WORK(KETAAB),WORK(KIJINT),
     *                                WORK(KAIINT),WORK(KIAINT),
     *                                WORK(KABINT),WORK(KD2GIJ),
     *                                WORK(KD2GAI),WORK(KD2GIA),
     *                                WORK(KD2GAB),WORK(KT1AM),
     *                                WORK(KEND5),LWRK5,ISYM)
C
                     ELSE IF (IOPT .EQ. 2) THEN
C
C-----------------------------------------------------------------
C                    Calculate direct 2 e- contribution to eta_ai.
C-----------------------------------------------------------------
C
                        ISYM = ISYMPQ
                        CALL CCDENZK0(ETAAI,WORK(KIJINT),WORK(KAIINT),
     *                                WORK(KIAINT),WORK(KABINT),
     *                                WORK(KD2GIJ),WORK(KD2GAI),
     *                                WORK(KD2GIA),WORK(KD2GAB),
     *                                WORK(KIJTIN),WORK(KABTIN),
     *                                WORK(KD2TIJ),WORK(KD2TAB),
     *                                WORK(KT1AM),WORK(KEND5),LWRK5,
     *                                ISYM)
C
                     ENDIF
C
  140             CONTINUE
  130          CONTINUE
C
               AUTIME = SECOND() - AUTIME
               TIMRES = TIMRES + AUTIME
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
C-----------------------
C     Regain work space.
C-----------------------
C
      KEND1 = KENDS2
      LWRK1 = LWRKS2
C
C--------------------------------------------------
C     Calculate the diagonal blocks of kappa-bar-0.
C--------------------------------------------------
C
c     IF ((FROIMP) .AND. (IOPT .EQ. 1)) THEN
c        CALL CCSD_ZKBLO(ZKDIA,WORK(KEND1),LWRK1)
c     ENDIF
C
C-----------------------------------------------------------
C     Calculate the correlated-frozen blocks of kappa-bar-0.
C-----------------------------------------------------------
C
c     IF (FROIMP) THEN
c        KOFFIJ = 2*NT1AMX + NMATIJ(1) + NMATAB(1) + 2*NT1FRO(1) + 1
c        CALL CC_ZKFCB(ZKDIA(KOFFIJ),WORK(KEND1),LWRK1)
C
C----------------------------------------------------------------
C        Calculate correction terms from correlated-frozen block.
C----------------------------------------------------------------
C
c        KOFFAI = 2*NT1AMX + NMATIJ(1) + NMATAB(1) + 1
c        CALL CC_ETACOR(ETAAI,ZKDIA(KOFFAI),ZKDIA(KOFFIJ),
c    *                  WORK(KCMOF),WORK(KEND1),LWRK1)
C
C---------------------
C     Reorder results.
C---------------------
C
c        CALL CC_ETARE(ETAAI,ZKDIA(KOFFAI),WORK(KEND1),LWRK1)
c     ENDIF
C
C---------------------------------
C     Write out eta-ai and eta-aI.
C---------------------------------
C
      IF (LOCDBG) THEN
C
         CALL AROUND('Eta-kappa-bar-0-ai vector exiting CC2_DEN')
C
         DO 20 ISYM = 1,NSYM
C
            WRITE(LUPRI,*) ' '
            WRITE(LUPRI,444) 'Sub-symmetry block number:', ISYM
            WRITE(LUPRI,555) '--------------------------'
  444       FORMAT(3X,A26,2X,I1)
  555       FORMAT(3X,A25)
C
            KOFF = IALLAI(ISYM,ISYM) + 1
            CALL OUTPUT(ETAAI(KOFF),1,NVIR(ISYM),1,NRHFS(ISYM),
     *                  NVIR(ISYM),NRHFS(ISYM),1,LUPRI)
C
            IF ((NVIR(ISYM) .EQ. 0) .OR. (NRHFS(ISYM) .EQ. 0)) THEN
               WRITE(LUPRI,*) 'This sub-symmetry is empty'
            ENDIF
C
  20     CONTINUE
      ENDIF
C
      IF (IPRINT .GT. 5) THEN
         ETAKAN = DDOT(NALLAI(1),ETAAI,1,ETAAI,1)
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,*) 'Norm of occupied-virtual block:', ETAKAN
      ENDIF
C
C-------------------------------------------
C     Write out frozen block of kappa-bar-0.
C-------------------------------------------
C
c     IF ((IPRINT .GT. 9) .AND. (FROIMP) .AND. (IOPT .EQ. 1)) THEN
c        CALL AROUND('Zeta-kappa-0 correlated-frozen block')
c        KOFRES = 2*NT1AMX + NMATIJ(1) + NMATAB(1) + 2*NT1FRO(1) + 1
c        ZKAPF1 = DDOT(NCOFRO(1),ZKDIA(KOFRES),1,
c    *                 ZKDIA(KOFRES),1)
c        ZKAPFR = ZKAPF1**0.5
c        WRITE(LUPRI,*) 'Norm of correlated-frozen block:', ZKAPFR
c     ENDIF
C
C-----------------------
C     Write out timings.
C-----------------------
C
  99  TIMETO = SECOND() - TIMETO
C
      IF (IPRINT .GT. 3) THEN
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,*) 
     *         'Requested density dependent quantities calculated'
         WRITE(LUPRI,*) 'Total time used in CC2_DEN:', TIMETO
      ENDIF
      IF (IPRINT .GT. 5) THEN
         WRITE(LUPRI,'(A,F15.2)')
     &      ' Time used for setting up 2 e- density      :', TIMDEN
     &     ,' Time used for contraction with integrals   :', TIMRES
     &     ,' Time used for reading 2 e- AO-integrals    :', TIRDAO
     &     ,' Time used for calculating 2 e- AO-integrals:', TIMHE2 
C    &     ,' Time used for 1 e- density & intermediates :', TIMONE
      ENDIF
C
C----------------------------------------------
C     Calculate the ZK0(ab) and ZK0(ij) blocks.
C----------------------------------------------
C
      IF (IOPT .EQ. 1) THEN
         CALL MP2_ZKBLO(ZKDIA,WORK(KETAIJ),WORK(KETAAB),WORK(KEND1),
     *                  LWRK1)
      ENDIF
C
C--------------------------------------------------- 
C     Calculate frozen core occupied blocks ZK0(iJ).
C--------------------------------------------------- 
C
      IF (FROIMP) THEN
          CALL QUIT('NO FROZEN CORE FOR RELAXED CC2 YET')
c         KOFRES = NMATIJ(1) + NMATAB(1) + 2*NT1AMX + 2*NT1FRO(1) + 1
c         CALL MP2_ZKFCB(ZKDIA(KOFRES),WORK(KZ2AM),WORK(KEND1),LWRK1)
      ENDIF
C
C------------------------------------------------
C     Write out timings and results if requested.
C------------------------------------------------
C
      IF ((IPRINT .GT. 3) .AND. (IOPT .EQ. 1)) THEN
C
         CALL AROUND('eta_ij and eta_ab blocks')
         ETAIJN = DDOT(NMATIJ(1),WORK(KETAIJ),1,WORK(KETAIJ),1)
         ETAABN = DDOT(NMATAB(1),WORK(KETAAB),1,WORK(KETAAB),1)
         ETAIJN = ETAIJN**0.5
         ETAABN = ETAABN**0.5
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,*) 'Norm of occupied-occupied block:', ETAIJN
         WRITE(LUPRI,*) 'Norm of virtual-virtual block  :', ETAABN
C
         CALL AROUND('Zeta-kappa-0 diagonal blocks')
         ZKAPI1 = DDOT(NMATIJ(1),ZKDIA(1),1,ZKDIA(1),1)
         ZKAPA1 = DDOT(NMATAB(1),ZKDIA(NMATIJ(1)+1),1,
     *                 ZKDIA(NMATIJ(1)+1),1)
         ZKAPIJ = ZKAPI1**0.5
         ZKAPAB = ZKAPA1**0.5
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,*) 'Norm of occupied-occupied block:', ZKAPIJ
         WRITE(LUPRI,*) 'Norm of virtual-virtual block  :', ZKAPAB
C
c         IF (FROIMP) THEN
c            ZKAPF1 = DDOT(NCOFRO(1),ZKDIA(KOFRES),1,
c     *                    ZKDIA(KOFRES),1)
c            ZKAPFR = ZKAPF1**0.5
c         WRITE(LUPRI,*) 'Norm of frozen-core-occupied block:', ZKAPFR
c         ENDIF
C
         IF (IPRINT .GT. 50) THEN
            DO 240 ISYM = 1,NSYM
               WRITE(LUPRI,*) ' '
               WRITE(LUPRI,*) 'Symmetry block:', ISYM
               KIJ = IMATIJ(ISYM,ISYM) + 1
               KAB = IMATAB(ISYM,ISYM) + 1 + NMATIJ(1)
               CALL AROUND('occ-occ block')
               CALL OUTPUT(ZKDIA(KIJ),1,NRHF(ISYM),1,NRHF(ISYM),
     *                     NRHF(ISYM),NRHF(ISYM),1,LUPRI)
               CALL AROUND('vir-vir block')
               CALL OUTPUT(ZKDIA(KAB),1,NVIR(ISYM),1,NVIR(ISYM),
     *                     NVIR(ISYM),NVIR(ISYM),1,LUPRI)
  240       CONTINUE
         ENDIF
      ENDIF
C
      TIMETO = SECOND() - TIMETO
C
      IF (IPRINT .GT. 3) THEN
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,*) 'Diagonal blocks of Zeta-kappa-0 calculated'
         WRITE(LUPRI,*) 'Total time used in CC2_DEN:', TIMETO
      ENDIF
C
      CALL QEXIT('CC2_DEN')
C
      RETURN
      END
C  /* Deck cc2_etij */
      SUBROUTINE CC2_ETIJ(ETAIJ,XINTIJ,XINTAI,XINTIA,XINTAB,DIJ,DAI,
     *                    DIA,DAB,T1AM,WORK,LWORK,ISYM)
C
C     Written by A. Halkier & S. Coriani 14/1-2000
C
C     Version: 1.0
C
C     Purpose: To set up the right hand side of the equation for
C              zeta-kappa-0_ij (ETAIJ) from MO-integrals (XIN*), CC2
C              densities (D*) and t1-amplitudes (T1AM)!
C              Note that due to the symmetry in the formulas, this
C              routine is able to handle both the one- and the two-
C              electron contributions!
C              ISYM is the symmetry of both the density and the
C              integrals!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION ETAIJ(*), XINTIJ(*), XINTAI(*), XINTIA(*), XINTAB(*)
      DIMENSION DIJ(*), DAI(*), DIA(*), DAB(*), T1AM(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
C
      CALL QENTER('CC2_ETIJ')
C
      DO 100 ISYMI = 1,NSYM
C
C----------------------------------------------------------------
C        Calculate direct terms - those without T1AM - to eta_ij.
C----------------------------------------------------------------
C
         ISYMJ  = ISYMI
         ISYMK  = MULD2H(ISYMI,ISYM)
         ISYMC  = MULD2H(ISYMI,ISYM)
C
         KOFFRE = IMATIJ(ISYMI,ISYMJ) + 1
C
         NTOTRE = MAX(NRHF(ISYMI),1)
         NTOTI  = MAX(NRHF(ISYMI),1)
         NTOTJ  = MAX(NRHF(ISYMJ),1)
         NTOTK  = MAX(NRHF(ISYMK),1)
         NTOTC  = MAX(NVIR(ISYMC),1)
C
         KOFF1  = IMATIJ(ISYMK,ISYMI) + 1
         KOFF2  = IMATIJ(ISYMK,ISYMJ) + 1
         KOFF3  = IMATIJ(ISYMI,ISYMK) + 1
         KOFF4  = IMATIJ(ISYMJ,ISYMK) + 1
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NRHF(ISYMK),ONE,
     *              XINTIJ(KOFF1),NTOTK,DIJ(KOFF2),NTOTK,ONE,
     *              ETAIJ(KOFFRE),NTOTRE)
C
         CALL DGEMM('N','T',NRHF(ISYMI),NRHF(ISYMJ),NRHF(ISYMK),ONE,
     *              XINTIJ(KOFF3),NTOTI,DIJ(KOFF4),NTOTJ,ONE,
     *              ETAIJ(KOFFRE),NTOTRE)
C
         CALL DGEMM('N','T',NRHF(ISYMI),NRHF(ISYMJ),NRHF(ISYMK),-ONE,
     *              DIJ(KOFF3),NTOTI,XINTIJ(KOFF4),NTOTJ,ONE,
     *              ETAIJ(KOFFRE),NTOTRE)
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NRHF(ISYMK),-ONE,
     *              DIJ(KOFF1),NTOTK,XINTIJ(KOFF2),NTOTK,ONE,
     *              ETAIJ(KOFFRE),NTOTRE)
C
         KOFF5  = IT1AM(ISYMC,ISYMI) + 1
         KOFF6  = IT1AM(ISYMC,ISYMJ) + 1
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMC),ONE,
     *              XINTAI(KOFF5),NTOTC,DAI(KOFF6),NTOTC,ONE,
     *              ETAIJ(KOFFRE),NTOTRE)
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMC),ONE,
     *              XINTIA(KOFF5),NTOTC,DIA(KOFF6),NTOTC,ONE,
     *              ETAIJ(KOFFRE),NTOTRE)
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMC),-ONE,
     *              DIA(KOFF5),NTOTC,XINTIA(KOFF6),NTOTC,ONE,
     *              ETAIJ(KOFFRE),NTOTRE)
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMC),-ONE,
     *              DAI(KOFF5),NTOTC,XINTAI(KOFF6),NTOTC,ONE,
     *              ETAIJ(KOFFRE),NTOTRE)
C
  100 CONTINUE
C
      DO 110 ISYMI = 1,NSYM
C
C---------------------------------------------------------------------------
C     Calculate terms involving T1AM and intermediate looping over occupied.
C     Note that the intermediate of integrals and densities is used for both
C     contributions as ISYMI=ISYMJ.
C---------------------------------------------------------------------------
C
         ISYMJ  = ISYMI
         ISYMB  = ISYMJ
         ISYMK  = MULD2H(ISYMB,ISYM)
C
         KOFFRE = IMATIJ(ISYMI,ISYMJ) + 1
C
         NTOTRE = MAX(NRHF(ISYMI),1)
         NTOTI  = MAX(NRHF(ISYMI),1)
         NTOTJ  = MAX(NRHF(ISYMJ),1)
         NTOTK  = MAX(NRHF(ISYMK),1)
         NTOTB  = MAX(NVIR(ISYMB),1)
C
C----------------------------------
C        Work space allocation one.
C----------------------------------
C
         KSCRHD = 1
         KEND1  = KSCRHD + NVIR(ISYMB)*NRHF(ISYMI)
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT(
     *         'Insufficient memory for allocation (1) in CC2_ETIJ')
         ENDIF
C
         CALL DZERO(WORK(KSCRHD),NVIR(ISYMB)*NRHF(ISYMI))
C
         KOFF1  = IT1AM(ISYMB,ISYMK)  + 1
         KOFF2  = IMATIJ(ISYMK,ISYMI) + 1
         KOFF3  = IMATIJ(ISYMI,ISYMK) + 1
         KOFF4  = IT1AM(ISYMB,ISYMJ)  + 1
C
         CALL DGEMM('N','N',NVIR(ISYMB),NRHF(ISYMI),NRHF(ISYMK),ONE,
     *              XINTIA(KOFF1),NTOTB,DIJ(KOFF2),NTOTK,ZERO,
     *              WORK(KSCRHD),NTOTB)
C
         CALL DGEMM('N','T',NVIR(ISYMB),NRHF(ISYMI),NRHF(ISYMK),-ONE,
     *              DAI(KOFF1),NTOTB,XINTIJ(KOFF3),NTOTI,ONE,
     *              WORK(KSCRHD),NTOTB)
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMB),ONE,
     *              WORK(KSCRHD),NTOTB,T1AM(KOFF4),NTOTB,ONE,
     *              ETAIJ(KOFFRE),NTOTRE)
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMB),-ONE,
     *              T1AM(KOFF4),NTOTB,WORK(KSCRHD),NTOTB,ONE,
     *              ETAIJ(KOFFRE),NTOTRE)
C
  110 CONTINUE
C
      DO 120 ISYMI = 1,NSYM
C
C---------------------------------------------------------------------------
C     Calculate terms involving T1AM and intermediate looping over virtual.
C     Note that the intermediate of integrals and densities is used for both
C     contributions as ISYMI=ISYMJ.
C---------------------------------------------------------------------------
C
         ISYMJ  = ISYMI
         ISYMB  = ISYMJ
         ISYMC  = MULD2H(ISYMB,ISYM)
C
         KOFFRE = IMATIJ(ISYMI,ISYMJ) + 1
C
         NTOTRE = MAX(NRHF(ISYMI),1)
         NTOTI  = MAX(NRHF(ISYMI),1)
         NTOTJ  = MAX(NRHF(ISYMJ),1)
         NTOTC  = MAX(NVIR(ISYMC),1)
         NTOTB  = MAX(NVIR(ISYMB),1)
C
C----------------------------------
C        Work space allocation two.
C----------------------------------
C
         KSCRHD = 1
         KEND1  = KSCRHD + NVIR(ISYMB)*NRHF(ISYMI)
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT(
     *         'Insufficient memory for allocation (2) in CC2_ETIJ')
         ENDIF
C
         CALL DZERO(WORK(KSCRHD),NVIR(ISYMB)*NRHF(ISYMI))
C
         KOFF1  = IMATAB(ISYMC,ISYMB) + 1
         KOFF2  = IT1AM(ISYMC,ISYMI)  + 1
         KOFF3  = IMATAB(ISYMB,ISYMC) + 1
         KOFF4  = IT1AM(ISYMB,ISYMJ)  + 1
C
         CALL DGEMM('T','N',NVIR(ISYMB),NRHF(ISYMI),NVIR(ISYMC),ONE,
     *              XINTAB(KOFF1),NTOTC,DAI(KOFF2),NTOTC,ZERO,
     *              WORK(KSCRHD),NTOTB)
C
         CALL DGEMM('N','N',NVIR(ISYMB),NRHF(ISYMI),NVIR(ISYMC),-ONE,
     *              DAB(KOFF3),NTOTB,XINTIA(KOFF2),NTOTC,ONE,
     *              WORK(KSCRHD),NTOTB)
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMB),ONE,
     *              WORK(KSCRHD),NTOTB,T1AM(KOFF4),NTOTB,ONE,
     *              ETAIJ(KOFFRE),NTOTRE)
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMB),-ONE,
     *              T1AM(KOFF4),NTOTB,WORK(KSCRHD),NTOTB,ONE,
     *              ETAIJ(KOFFRE),NTOTRE)
C
  120 CONTINUE
C
      CALL QEXIT('CC2_ETIJ')
C
      RETURN
      END
C  /* Deck cc2_etab */
      SUBROUTINE CC2_ETAB(ETAAB,XINTIJ,XINTAI,XINTIA,XINTAB,DIJ,DAI,
     *                    DIA,DAB,T1AM,WORK,LWORK,ISYM)
C
C     Written by A. Halkier & S. Coriani 15/1-2000
C
C     Version: 1.0
C
C     Purpose: To set up the right hand side of the equation for
C              zeta-kappa-0_ab (ETAAB) from MO-integrals (XIN*), CC2
C              densities (D*) and t1-amplitudes (T1AM)!
C              Note that due to the symmetry in the formulas, this
C              routine is able to handle both the one- and the two-
C              electron contributions!
C              ISYM is the symmetry of both the density and the
C              integrals!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION ETAAB(*), XINTIJ(*), XINTAI(*), XINTIA(*), XINTAB(*)
      DIMENSION DIJ(*), DAI(*), DIA(*), DAB(*), T1AM(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
C
      CALL QENTER('CC2_ETAB')
C
      DO 100 ISYMA = 1,NSYM
C
C----------------------------------------------------------------
C        Calculate direct terms - those without T1AM - to eta_ab.
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
         NTOTC  = MAX(NVIR(ISYMC),1)
         NTOTK  = MAX(NRHF(ISYMK),1)
C
         KOFF1  = IT1AM(ISYMA,ISYMK) + 1
         KOFF2  = IT1AM(ISYMB,ISYMK) + 1
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMK),ONE,
     *              XINTIA(KOFF1),NTOTA,DIA(KOFF2),NTOTB,ONE,
     *              ETAAB(KOFFRE),NTOTRE)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMK),-ONE,
     *              DAI(KOFF1),NTOTA,XINTAI(KOFF2),NTOTB,ONE,
     *              ETAAB(KOFFRE),NTOTRE)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMK),-ONE,
     *              DIA(KOFF1),NTOTA,XINTIA(KOFF2),NTOTB,ONE,
     *              ETAAB(KOFFRE),NTOTRE)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMK),ONE,
     *              XINTAI(KOFF1),NTOTA,DAI(KOFF2),NTOTB,ONE,
     *              ETAAB(KOFFRE),NTOTRE)
C
         KOFF3  = IMATAB(ISYMC,ISYMA) + 1
         KOFF4  = IMATAB(ISYMC,ISYMB) + 1
         KOFF5  = IMATAB(ISYMA,ISYMC) + 1
         KOFF6  = IMATAB(ISYMB,ISYMC) + 1
C
         CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMB),NVIR(ISYMC),ONE,
     *              XINTAB(KOFF3),NTOTC,DAB(KOFF4),NTOTC,ONE,
     *              ETAAB(KOFFRE),NTOTRE)
C 
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NVIR(ISYMC),-ONE,
     *              DAB(KOFF5),NTOTA,XINTAB(KOFF6),NTOTB,ONE,
     *              ETAAB(KOFFRE),NTOTRE)
C
         CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMB),NVIR(ISYMC),-ONE,
     *              DAB(KOFF3),NTOTC,XINTAB(KOFF4),NTOTC,ONE,
     *              ETAAB(KOFFRE),NTOTRE)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NVIR(ISYMC),ONE,
     *              XINTAB(KOFF5),NTOTA,DAB(KOFF6),NTOTB,ONE,
     *              ETAAB(KOFFRE),NTOTRE)
C
  100 CONTINUE
C
      DO 110 ISYMA = 1,NSYM
C
C---------------------------------------------------------------------------
C     Calculate terms involving T1AM and intermediate looping over occupied.
C     Note that the intermediate of integrals and densities is used for both
C     contributions as ISYMA=ISYMB.
C---------------------------------------------------------------------------
C
         ISYMB  = ISYMA
         ISYMJ  = ISYMB
         ISYMK  = MULD2H(ISYMA,ISYM)
C
         KOFFRE = IMATAB(ISYMA,ISYMB) + 1
C
         NTOTRE = MAX(NVIR(ISYMA),1)
         NTOTA  = MAX(NVIR(ISYMA),1)
         NTOTB  = MAX(NVIR(ISYMB),1)
         NTOTK  = MAX(NRHF(ISYMK),1)
         NTOTJ  = MAX(NRHF(ISYMJ),1)
C
C----------------------------------
C        Work space allocation one.
C----------------------------------
C
         KSCRHD = 1
         KEND1  = KSCRHD + NVIR(ISYMA)*NRHF(ISYMJ)
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT(
     *         'Insufficient memory for allocation (1) in CC2_ETAB')
         ENDIF
C
         CALL DZERO(WORK(KSCRHD),NVIR(ISYMA)*NRHF(ISYMJ))
C
         KOFF1 = IT1AM(ISYMA,ISYMK)  + 1
         KOFF2 = IMATIJ(ISYMK,ISYMJ) + 1
         KOFF3 = IMATIJ(ISYMJ,ISYMK) + 1
         KOFF4 = IT1AM(ISYMB,ISYMJ)  + 1
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMJ),NRHF(ISYMK),ONE,
     *              XINTIA(KOFF1),NTOTA,DIJ(KOFF2),NTOTK,ZERO,
     *              WORK(KSCRHD),NTOTA)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NRHF(ISYMJ),NRHF(ISYMK),-ONE,
     *              DAI(KOFF1),NTOTA,XINTIJ(KOFF3),NTOTJ,ONE,
     *              WORK(KSCRHD),NTOTA)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMJ),ONE,
     *              WORK(KSCRHD),NTOTA,T1AM(KOFF4),NTOTB,ONE,
     *              ETAAB(KOFFRE),NTOTRE)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMJ),-ONE,
     *              T1AM(KOFF4),NTOTA,WORK(KSCRHD),NTOTB,ONE,
     *              ETAAB(KOFFRE),NTOTRE)
C
  110 CONTINUE
C
      DO 120 ISYMA = 1,NSYM
C
C---------------------------------------------------------------------------
C     Calculate terms involving T1AM and intermediate looping over virtual.
C     Note that the intermediate of integrals and densities is used for both
C     contributions as ISYMA=ISYMB.
C---------------------------------------------------------------------------
C
         ISYMB  = ISYMA
         ISYMJ  = ISYMB
         ISYMC  = MULD2H(ISYMA,ISYM)
C
         KOFFRE = IMATAB(ISYMA,ISYMB) + 1
C
         NTOTRE = MAX(NVIR(ISYMA),1)
         NTOTA  = MAX(NVIR(ISYMA),1)
         NTOTB  = MAX(NVIR(ISYMB),1)
         NTOTC  = MAX(NVIR(ISYMC),1)
C
C----------------------------------
C        Work space allocation one.
C----------------------------------
C
         KSCRHD = 1
         KEND1  = KSCRHD + NVIR(ISYMA)*NRHF(ISYMJ)
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT(
     *         'Insufficient memory for allocation (2) in CC2_ETAB')
         ENDIF
C
         CALL DZERO(WORK(KSCRHD),NVIR(ISYMA)*NRHF(ISYMJ))
C
         KOFF1 = IMATAB(ISYMC,ISYMA) + 1
         KOFF2 = IT1AM(ISYMC,ISYMJ)  + 1
         KOFF3 = IMATAB(ISYMA,ISYMC) + 1
         KOFF4 = IT1AM(ISYMB,ISYMJ)  + 1
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMJ),NVIR(ISYMC),ONE,
     *              XINTAB(KOFF1),NTOTC,DAI(KOFF2),NTOTC,ZERO,
     *              WORK(KSCRHD),NTOTA)
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMJ),NVIR(ISYMC),-ONE,
     *              DAB(KOFF3),NTOTA,XINTIA(KOFF2),NTOTC,ONE,
     *              WORK(KSCRHD),NTOTA)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMJ),ONE,
     *              WORK(KSCRHD),NTOTA,T1AM(KOFF4),NTOTB,ONE,
     *              ETAAB(KOFFRE),NTOTRE)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMJ),-ONE,
     *              T1AM(KOFF4),NTOTA,WORK(KSCRHD),NTOTB,ONE,
     *              ETAAB(KOFFRE),NTOTRE)
C
  120 CONTINUE
C
      CALL QEXIT('CC2_ETAB')
C
      RETURN
      END
