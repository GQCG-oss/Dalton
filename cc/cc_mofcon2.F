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
      SUBROUTINE CC_MOFCON2(XINT,OMEGA2,
     *                      XLAMPA,XLAMHA,XLAMPB,XLAMHB,XLAMPC,XLAMHD,
     *                      ISYMXLA,ISYMXLB,ISYMXLC,ISYMXLD,
     *                      WORK,LWORK,IDEL,ISYMD,ISYOMEG,ISYHOP,IOPT)
C
C  Written by Asger Halkier and Henrik Koch 3-5-95.
C
C  Debugged By Ove Christiansen 25-7-1995
C  Generalized for cubic response calculations by C.H. in January 1997
C  Generalized for non-total symmetric integrals by C.H. in July 1998
C
C  Purpose: To calculate the F-term's contribution to the vector 
C           function and its derivatives using matrix vector routines.
C
C
C    IOPT=0: no symmetrization, returns:
C            F_{aibj} =   (a^A i^B | b^C j^D)
C
C    IOPT=1: symmetrization in (ai) <-> (bj), returns:
C            F_{aibj} =   (a^A i^B | b^C j^D) + (b^A j^B | a^C i^D)
C
C    IOPT=2: symmetrization in A <-> B, returns:
C            F_{aibj} =   (a^A i^B | b^C j^D) + (a^B i^A | b^C j^D)
C
C    IOPT=3: symmetrization in (ai) <-> (bj) and in A <-> B, returns:
C            F_{aibj} =   (a^A i^B | b^C j^D) + (b^A j^B | a^C i^D) 
C                       + (a^B i^A | b^C j^D) + (b^B j^A | a^C i^D)
C
C    N.B.: IOPT=0,2 assumes XLAMPA = XLAMPC and XLAMHB = XLAMHD
C
C
C  Symmetries:    ISYOMEG  --  OMEGA2
C                 ISYMXLA  --  XLAMPA, XLAMHA
C                 ISYMXLB  --  XLAMPB, XLAMHB
C                 ISYMXLC  --  XLAMPC
C                 ISYMXLD  --  XLAMHD
C                 ISYMD    --  AO IDEL
C                 ISYHOP   --  symmetry of integrals (all 4 indeces)
C
C
C  N.B. This routine assumes AO-symmetric integrals, and can therefor
C       not be used directly for calculations with London-orbitals!!!
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
      PARAMETER(ZERO = 0.0D0,ONE = 1.0D0,XMONE=-1.0D0,TWO = 2.0D0)
      DIMENSION XINT(*),OMEGA2(*)
      DIMENSION XLAMPA(*),XLAMHA(*),XLAMPB(*),XLAMHB(*),
     *          XLAMPC(*),XLAMHD(*)
      DIMENSION WORK(LWORK)
#include <ccorb.h>
#include <symsq.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      ISYDIS = MULD2H(ISYMD,ISYHOP)
      ISYXAB = MULD2H(ISYMXLA,ISYMXLB)
      ISYXCD = MULD2H(ISYMXLC,ISYMXLD)

      IF (MULD2H(ISYXAB,ISYXCD) .NE. MULD2H(ISYOMEG,ISYHOP)) THEN
        CALL QUIT('SYMMETRY MISMATCH IN CC_MOFCON2.')
      END IF
C
      IF (IOPT.LT.0 .OR. IOPT.GT.3) THEN
        CALL QUIT('CC_MOFCON2 called with an illegal value for IOPT.')
      END IF
      IF (IOPT.EQ.0 .OR. IOPT.EQ.2) THEN
        IF (ISYMXLA.NE.ISYMXLC .OR. ISYMXLB.NE.ISYMXLD) THEN
          CALL QUIT('CC_MOFCON2 called with inconsitent symmetries.')
        END IF
      END IF
C
      DO 100 ISYMG = 1,NSYM
C
         IF (NBAS(ISYMG) .EQ. 0) GOTO 100
C
         ISALBE = MULD2H(ISYMG,ISYDIS)
         ISYMAI = MULD2H(ISALBE,ISYXAB)
         ISYMJ  = MULD2H(ISYMG,ISYMXLD)
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
         KOFF2 = IGLMRH(ISYMG,ISYMJ) + 1
C
         NTALBE = MAX(NNBST(ISALBE),1)
         NTOTG  = MAX(NBAS(ISYMG),1)
C
         CALL DGEMM('N','N',NNBST(ISALBE),NRHF(ISYMJ),NBAS(ISYMG),
     *              ONE,XINT(KOFF1),NTALBE,XLAMHD(KOFF2),NTOTG,
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
               ISYMBE = MULD2H(ISYMI,ISYMXLB)
               ISYMAL = MULD2H(ISYMBE,ISALBE)
               ISYMA  = MULD2H(ISYMAL,ISYMXLA)
C
               IF (LWRK1 .LT. NBAS(ISYMAL)*NRHF(ISYMI)) THEN
                  CALL QUIT('Insufficient space for 2. trf. '//
     &                      'in CC_MOFCON')
               ENDIF
C
               KOFF2 = KSCR2 + IAODIS(ISYMAL,ISYMBE)
               KOFF3 = IGLMRH(ISYMBE,ISYMI) + 1
               KOFF4 = IGLMVI(ISYMAL,ISYMA) + 1
               KOFF5 = KSCR3 + IT1AM(ISYMA,ISYMI)
C
               NTOTAL = MAX(NBAS(ISYMAL),1)
               NTOTBE = MAX(NBAS(ISYMBE),1)
               NTOTA  = MAX(NVIR(ISYMA),1)
C
               CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMI),NBAS(ISYMBE),
     *                    ONE,WORK(KOFF2),NTOTAL,XLAMHB(KOFF3),NTOTBE,
     *                    ZERO,WORK(KEND1),NTOTAL)
C
               CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMAL),
     *                    ONE,XLAMPA(KOFF4),NTOTAL,WORK(KEND1),NTOTAL,
     *                    ZERO,WORK(KOFF5),NTOTA)
C
C-----------------------
C Symmetrize in A <-> B.
C-----------------------
               IF (IOPT .EQ. 2 .OR. IOPT .EQ. 3) THEN
C
                  ISYMBE = MULD2H(ISYMI,ISYMXLA)
                  ISYMAL = MULD2H(ISYMBE,ISALBE)
                  ISYMA  = MULD2H(ISYMAL,ISYMXLB)
C
                  IF (LWRK1 .LT. NBAS(ISYMAL)*NRHF(ISYMI)) THEN
                     CALL QUIT('Insufficient space for 2. trf. '//
     &                         'in CC_MOFCON')
                  ENDIF
C
                  KOFF2 = KSCR2 + IAODIS(ISYMAL,ISYMBE)
                  KOFF3 = IGLMRH(ISYMBE,ISYMI) + 1
                  KOFF4 = IGLMVI(ISYMAL,ISYMA) + 1
                  KOFF5 = KSCR3 + IT1AM(ISYMA,ISYMI)
C
                  NTOTAL = MAX(NBAS(ISYMAL),1)
                  NTOTBE = MAX(NBAS(ISYMBE),1)
                  NTOTA  = MAX(NVIR(ISYMA),1)
C
                  CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMI),
     *                       NBAS(ISYMBE),ONE,WORK(KOFF2),NTOTAL,
     *                       XLAMHA(KOFF3),NTOTBE,ZERO,WORK(KEND1),
     *                       NTOTAL)
C
                  CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),
     *                       NBAS(ISYMAL),ONE,XLAMPB(KOFF4),NTOTAL,
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
            ISYMB  = MULD2H(ISYMD,ISYMXLC)
            ISYMBJ = MULD2H(ISYMB,ISYMJ)
C
            DO 130 B = 1,NVIR(ISYMB)
C
               NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
               NDB = IGLMVI(ISYMD,ISYMB) + NBAS(ISYMD)*(B - 1)
     *             + IDEL - IBAS(ISYMD)
C
               CALL DZERO(WORK(KSCR4),NT1AM(ISYMAI))
C
               XLB  = XLAMPC(NDB)
C
               CALL DAXPY(NT1AM(ISYMAI),XLB,WORK(KSCR3),1,WORK(KSCR4),1)
C
               IF (ISYMBJ .EQ. ISYMAI) THEN
C
                  NTOTAI = NBJ
C
                  IF (IOPT .EQ. 1 .OR. IOPT .EQ. 3) THEN
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
C
               ENDIF
C
               IF (ISYMAI .LT. ISYMBJ) THEN
C
                  DO 150 NAI = 1,NT1AM(ISYMAI)
C
                     NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                     + NT1AM(ISYMAI)*(NBJ - 1) + NAI
C
                     OMEGA2(NAIBJ) = OMEGA2(NAIBJ) + WORK(KSCR4+NAI-1)
C
  150             CONTINUE
C
               ENDIF
C
               IF ((ISYMBJ.LT.ISYMAI).AND.(IOPT.EQ.1.OR.IOPT.EQ.3))THEN
C
                  DO 160 NAI = 1,NT1AM(ISYMAI)
C
                     NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                     + NT1AM(ISYMBJ)*(NAI - 1) + NBJ
C
                     OMEGA2(NAIBJ) = OMEGA2(NAIBJ) + WORK(KSCR4+NAI-1)
C
  160             CONTINUE
C
               ENDIF
C
  130       CONTINUE
C
  110    CONTINUE
C
  100 CONTINUE
C
      RETURN
      END
