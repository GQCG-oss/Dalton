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
C  /* Deck t2am_reorder */
      SUBROUTINE T2AM_REORDER(T2AMOLD,T2AMNEW,IPRINT)
C
C     Stephan P.A. Sauer.                    20-Jun-1995
C
C     Reorder the T2 amplitudes for SOPPA.
C     Outer loop over occupied orbitals, inner loop over virtual
C     orbitals.
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION T2AMOLD(*),T2AMNEW(*)
#include <inforb.h>
#include <ccsdsym.h>
C
      DIMENSION JB(8)
C
C
      DO 100 ISYM = 1,NSYM
         JB(ISYM) = 0
  100 CONTINUE
C
      IF (IPRINT .GT. 20) THEN
         CALL AROUND('for SOPPA : two coulomb minus exchange of t2am')
      END IF
C
C
      JBNOFF = 0
      DO 200 ISYMJ = 1,NSYM
C
         DO 210 J = 1,NRHF(ISYMJ)
C
            DO 220 ISYMB = 1,NSYM
C
               ISYMBJ = MULD2H(ISYMB,ISYMJ)
C
               DO 230 B = 1,NVIR(ISYMB)
                  JB(ISYMBJ) = JB(ISYMBJ) + 1
                  JBNEW  = JBNOFF + 1
                  JBOLD  = IT2SQ(ISYMBJ,ISYMBJ)
     *                     + (JB(ISYMBJ)-1) * NT1AM(ISYMBJ) + 1
C
                  CALL DCOPY(NT1AM(ISYMBJ),T2AMOLD(JBOLD),1,
     *                       T2AMNEW(JBNEW),1)
                  JBNOFF = JBNOFF + NT1AM(ISYMBJ)
C
  230          CONTINUE
C
               IF (IPRINT .GT. 20 .AND. NVIR(ISYMB) .GT.0) THEN
                  WRITE (LUPRI,*)
                  WRITE (LUPRI,*) 'Symmetry block:',ISYMBJ
                  JBSTART = JBNOFF + 1 - NVIR(ISYMB)*NT1AM(ISYMBJ)
                  CALL OUTPUT(T2AMNEW(JBSTART),1,NT1AM(ISYMBJ),1,
     *                        NVIR(ISYMB),NT1AM(ISYMBJ),
     *                        NVIR(ISYMB),1,6)
               END IF
  220       CONTINUE
C
  210    CONTINUE
C
  200 CONTINUE
C
      RETURN
      END

C  /* Deck soppa_density */
      SUBROUTINE SOPPA_DENSITY(DONE,T1AM,T2AM,T2TP,IPRINT)
C
C     Stephan P.A. Sauer.                    21-Jun-1995
C
C     Calculates the second-order one electron density matrix.
C     The density matrix is stored as the full square matrix.
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER(TWO=2.0D0)
      DIMENSION DONE(*),T1AM(*),T2AM(*),T2TP(*)
#include <inforb.h>
#include <ccsdsym.h>
C
C-------------------------------
C     Initialize density matrix.
C-------------------------------
C
      CALL DZERO(DONE,N2ORBX)
cSPAS
c      KOFF   = 0
C
c      DO 30 ISYMJ = 1,NSYM
C
c         ISYMI = ISYMJ
C
c         DO 31 J = 1,NRHF(ISYMJ)
C
c            NJJ = KOFF + NORBT*(J - 1) + J
C
c            DONE(NJJ) = DONE(NJJ) + TWO
C
c  31     CONTINUE
C
c         KOFF = KOFF + (NORBT + 1)*NORB(ISYMJ)
C
c  30  CONTINUE
CKeinSPASmehr
C
C==================================================
C     Adds the contribution from the T1 amplitudes.
C==================================================
C
      KOFF = 0
C
      DO 200 ISYMI = 1,NSYM
C
         ISYMA = ISYMI
C
         DO 210 I = 1,NRHF(ISYMI)
C
            KOFFA = IT1AM(ISYMI,ISYMA) + NVIR(ISYMA)*(I-1)
C
            DO 220 A = 1,NVIR(ISYMA)
C
               IA = NRHF(ISYMI) + A
C
               NIAONE = KOFF  + NORBT*(I - 1) + IA 
               NAI    = KOFFA + A
               DONE(NIAONE) = T1AM(NAI)
C
  220       CONTINUE
C
  210    CONTINUE
C
         DO 230 A = 1,NVIR(ISYMA)
C
            IA = NRHF(ISYMI) + A
C
            DO 240 I = 1,NRHF(ISYMI)
C
               NIAONE = KOFF + NORBT*(IA - 1) + I
               NAI    = IT1AM(ISYMI,ISYMA) + NVIR(ISYMA)*(I-1) + A
               DONE(NIAONE) = T1AM(NAI)
C
  240       CONTINUE
C
  230    CONTINUE
C
         KOFF = KOFF + (NORBT + 1)*NORB(ISYMI)
C
  200 CONTINUE
C
C========================================================
C     Calculates the contribution from the T2 amplitudes.
C========================================================
C
C
C-----------------------------------------
C     Occupied part of the density matrix.
C-----------------------------------------
C
      KOFF   = 0
C
      DO 300 ISYMJ = 1,NSYM
C
         ISYMI = ISYMJ
C
         DO 310 J = 1,NRHF(ISYMJ)
C
C           NJJ = KOFF + NORB(ISYMJ)*(J - 1) + J
            NJJ = KOFF + NORBT*(J - 1) + J
C
            DONE(NJJ) = DONE(NJJ) + TWO
C
  310    CONTINUE
C
         DO 320 J = 1,NRHF(ISYMJ)
C
            DO 330 I = 1,NRHF(ISYMI)
C
               DO 340 ISYMC = 1,NSYM
C
                  ISYMCI = MULD2H(ISYMC,ISYMI)
                  ISYMCJ = MULD2H(ISYMC,ISYMJ)
C
                  ISYMDK = ISYMCI
C
C                 NIJ = KOFF + NORB(ISYMI)*(J - 1) + I
                  NIJ = KOFF + NORBT*(J - 1) + I
                  NCI = IT1AM(ISYMC,ISYMI) + NVIR(ISYMC)*(I - 1) + 1
                  NCJ = IT1AM(ISYMC,ISYMJ) + NVIR(ISYMC)*(J - 1) + 1
C
                  KOFF1 = IT2SQ(ISYMDK,ISYMCI)
     *                  + NT1AM(ISYMDK)*(NCI - 1) + 1
                  KOFF2 = IT2SQ(ISYMDK,ISYMCJ)
     *                  + NT1AM(ISYMDK)*(NCJ - 1) + 1
C
                  NTOT = NT1AM(ISYMDK)*NVIR(ISYMC)
C
                  DONE(NIJ) = DONE(NIJ)
     *                      - TWO*DDOT(NTOT,T2AM(KOFF1),1,T2TP(KOFF2),1)
C
  340          CONTINUE
  330       CONTINUE
  320    CONTINUE
C
C        KOFF = KOFF + NORB(ISYMJ)*NORB(ISYMJ)
         KOFF = KOFF + (NORBT + 1)*NORB(ISYMJ)
C
  300 CONTINUE
C
C-----------------------------
C     Virtual part of density.
C-----------------------------
C
C
      DO 400 ISYML = 1,NSYM
C
         DO 410 L = 1,NRHF(ISYML)
C
            KOFF = 0
C
            DO 420 ISYMB = 1,NSYM
C
               ISYMA = ISYMB
C
               ISYMBL = MULD2H(ISYMB,ISYML)
               ISYMAL = ISYMBL
               ISYMCK = ISYMBL
C
               DO 430 B = 1,NVIR(ISYMB)
C
                  IB = NRHF(ISYMB) + B
C
                  NBL = IT1AM(ISYMB,ISYML) + NVIR(ISYMB)*(L - 1) + B
C
                  KOFF1 = IT2SQ(ISYMCK,ISYMBL)
     *                  + NT1AM(ISYMCK)*(NBL - 1) + 1
C
                  DO 440 A = 1,NVIR(ISYMA)
C
                     IA = NRHF(ISYMA) + A
C
C                    NAB = KOFF + NORB(ISYMB)*(IB -1) + IA
                     NAB = KOFF + NORBT*(IB -1) + IA
C
                     NAL = IT1AM(ISYMA,ISYML) + NVIR(ISYMA)*(L - 1) + A
C
                     KOFF2 = IT2SQ(ISYMCK,ISYMAL)
     *                     + NT1AM(ISYMCK)*(NAL - 1) + 1
C
                     DONE(NAB) = DONE(NAB)
     *                         + TWO*DDOT(NT1AM(ISYMCK),T2AM(KOFF1),1,
     *                                                  T2TP(KOFF2),1)
C
  440             CONTINUE
C
  430          CONTINUE
C
C              KOFF = KOFF + NORB(ISYMB)*NORB(ISYMB)
               KOFF = KOFF + (NORBT + 1)*NORB(ISYMB)
C
  420       CONTINUE
C
  410    CONTINUE
C
  400 CONTINUE
C
C-------------------
C     Print section.
C-------------------
C
      IF (IPRINT .GT. 20) THEN
         CALL AROUND('for SOPPA : One electron density matrix')
C
         KOFF = 1
         KORB = 0
C
         DO 500 ISYM = 1,NSYM
C
            WRITE (LUPRI,*)
            WRITE (LUPRI,*) 'Symmetry ',ISYM,' : block of density'
C
C           CALL OUTPUT(DONE(KOFF),1,NORB(ISYM),1,NORB(ISYM),
C    *                  NORB(ISYM),NORB(ISYM),1,6)
            CALL OUTPUT(DONE(KOFF),KORB+1,KORB+NORB(ISYM),1,NORB(ISYM),
     *                  NORBT,NORBT,1,6)
C
C           KOFF = KOFF + NORB(ISYM)*NORB(ISYM)
            KOFF = KOFF + NORBT*NORB(ISYM)
            KORB = KORB + NORB(ISYM)
C
  500    CONTINUE
C
      END IF
C
      RETURN
      END
