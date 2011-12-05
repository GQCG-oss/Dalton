C
C  /* Deck so_bextract */
      SUBROUTINE SO_BEXTRACT(EXTBE,LEXTBE,EXTBD,LEXTBD,
     &                       TR2E,LTR2E,TR2D,LTR2D,
     &                       ISYMA,ISYMI,ISYMB,ISYMJ,J)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, March 1996
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C
C     PURPOSE: Extract 2*b(a i b j) - b(b i a j) for given j
C              in EXTBE (and EXTBD) for excitations (and de-excitations)
C              respectively.
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION EXTBE(LEXTBE), EXTBD(LEXTBD), TR2E(LTR2E), TR2D(LTR2D)
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
C------------------------------
C     Statement function INDEX.
C------------------------------
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_BEXTRACT')
C
C----------------------------------------------------------------
C     Extract 2*b(a i b j) for excitations and de-excitations and
C     put in EXTBE and EXTBD.
C----------------------------------------------------------------
C
      ISYMAI = MULD2H(ISYMA,ISYMI)
      ISYMBJ = MULD2H(ISYMB,ISYMJ)
C
      DO 100 B = 1,NVIR(ISYMB)
C
         NBJ    = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB) * (J - 1) + B
C
         KOFF1 = NVIR(ISYMA) * NRHF(ISYMI) * (B - 1)
C
         IF ( ISYMAI .EQ. ISYMBJ ) THEN
C
            DO 110 I = 1,NRHF(ISYMI)
C
               KOFF2 = KOFF1 + NVIR(ISYMA) * (I - 1)
               NAI1  = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA) * (I - 1)
C
               DO 111 A = 1,NVIR(ISYMA)
C
                  NAI         = NAI1 + A
                  NAIBJ       = IT2AM(ISYMAI,ISYMBJ)
     &                        + INDEX(NAI,NBJ)
                  NAIB        = KOFF2 + A
C
                  EXTBE(NAIB) = TWO * TR2E(NAIBJ)
                  EXTBD(NAIB) = TWO * TR2D(NAIBJ)
C
  111          CONTINUE
C
  110       CONTINUE
C
         ELSE IF ( ISYMAI .LT. ISYMBJ ) THEN
C
            DO 120 I = 1,NRHF(ISYMI)
C
               KOFF2 = KOFF1 + NVIR(ISYMA) * (I - 1)
               NAI1  = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA) * (I - 1)
C
               DO 121 A = 1,NVIR(ISYMA)
C
                  NAI         = NAI1 + A
                  NAIBJ       = IT2AM(ISYMAI,ISYMBJ)
     &                        + NT1AM(ISYMAI) * (NBJ - 1) + NAI
                  NAIB        = KOFF2 + A
C
                  EXTBE(NAIB) =  TWO * TR2E(NAIBJ)
                  EXTBD(NAIB) =  TWO * TR2D(NAIBJ)
C
  121          CONTINUE
C
  120       CONTINUE
C
         ELSE IF ( ISYMAI .GT. ISYMBJ ) THEN
C
            DO 130 I = 1,NRHF(ISYMI)
C
               KOFF2 = KOFF1 + NVIR(ISYMA) * (I - 1)
               NAI1  = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA) * (I - 1)
C
               DO 131 A = 1,NVIR(ISYMA)
C
                  NAI         = NAI1 + A
                  NAIBJ       = IT2AM(ISYMBJ,ISYMAI)
     &                        + NT1AM(ISYMBJ) * (NAI - 1) + NBJ
                  NAIB        = KOFF2 + A
C
                  EXTBE(NAIB) = TWO * TR2E(NAIBJ)
                  EXTBD(NAIB) = TWO * TR2D(NAIBJ)
C
  131          CONTINUE
C
  130       CONTINUE
C
         END IF
C
  100 CONTINUE
C
C--------------------------------------------------------------
C     Extract b(b i a j) for excitations and de-excitations and
C     subtract from EXTBE and EXTBD to get the final result.
C--------------------------------------------------------------
C
      ISYMBI = MULD2H(ISYMB,ISYMI)
      ISYMAJ = MULD2H(ISYMA,ISYMJ)
C
      NORBAI = NVIR(ISYMA) * NRHF(ISYMI)
C
      DO 200 A = 1,NVIR(ISYMA)
C
         NAJ    = IT1AM(ISYMA,ISYMJ) + NVIR(ISYMA) * (J - 1) + A
C
         IF ( ISYMBI .EQ. ISYMAJ ) THEN
C
            DO 210 I = 1,NRHF(ISYMI)
C
               KOFF1 = NVIR(ISYMA) * (I - 1) + A
               NBI1  = IT1AM(ISYMB,ISYMI) + NVIR(ISYMB) * (I - 1)
C
C---------------------------------------------------------
C              The following if-construnction ensures that
C              if (i .eq. j) .and. (a .eq. b) the term is
C              NOT subtracted freom EXTBE and EXTBD.
C---------------------------------------------------------
C
               IF ( (ISYMI .EQ. ISYMJ) .AND. (I .EQ. J) ) THEN
C
                  DO 211 B = 1,A-1
C
                     NBI         = NBI1 + B
                     NBIAJ       = IT2AM(ISYMBI,ISYMAJ)
     &                           + INDEX(NBI,NAJ)
                     NAIB        = KOFF1 + NORBAI * (B - 1)
C
                     EXTBE(NAIB) = EXTBE(NAIB) - TR2E(NBIAJ)
                     EXTBD(NAIB) = EXTBD(NAIB) - TR2D(NBIAJ)
C
  211             CONTINUE
C
                  DO 212 B = A+1,NVIR(ISYMB)
C
                     NBI         = NBI1 + B
                     NBIAJ       = IT2AM(ISYMBI,ISYMAJ)
     &                           + INDEX(NBI,NAJ)
                     NAIB        = KOFF1 + NORBAI * (B - 1)
C
                     EXTBE(NAIB) = EXTBE(NAIB) - TR2E(NBIAJ)
                     EXTBD(NAIB) = EXTBD(NAIB) - TR2D(NBIAJ)
C
  212             CONTINUE
C
               ELSE
C
                  DO 213 B = 1,NVIR(ISYMB)
C
                     NBI         = NBI1 + B
                     NBIAJ       = IT2AM(ISYMBI,ISYMAJ)
     &                           + INDEX(NBI,NAJ)
                     NAIB        = KOFF1 + NORBAI * (B - 1)
C
                     EXTBE(NAIB) = EXTBE(NAIB) - TR2E(NBIAJ)
                     EXTBD(NAIB) = EXTBD(NAIB) - TR2D(NBIAJ)
C
  213             CONTINUE
C
               END IF
C
  210       CONTINUE
C
         ELSE IF ( ISYMBI .LT. ISYMAJ ) THEN
C
            DO 220 I = 1,NRHF(ISYMI)
C
               KOFF1 = NVIR(ISYMA) * (I - 1) + A
               NBI1  = IT1AM(ISYMB,ISYMI) + NVIR(ISYMB) * (I - 1)
C
               DO 221 B = 1,NVIR(ISYMB)
C
                  NBI         = NBI1 + B
                  NBIAJ       = IT2AM(ISYMBI,ISYMAJ)
     &                        + NT1AM(ISYMBI) * (NAJ - 1) + NBI
                  NAIB        = KOFF1 + NORBAI * (B - 1)
C
                  EXTBE(NAIB) =  EXTBE(NAIB) - TR2E(NBIAJ)
                  EXTBD(NAIB) =  EXTBD(NAIB) - TR2D(NBIAJ)
C
  221          CONTINUE
C
  220       CONTINUE
C
         ELSE IF ( ISYMBI .GT. ISYMAJ ) THEN
C
            DO 230 I = 1,NRHF(ISYMI)
C
               KOFF1 = NVIR(ISYMA) * (I - 1) + A
               NBI1  = IT1AM(ISYMB,ISYMI) + NVIR(ISYMB) * (I - 1)
C
               DO 231 B = 1,NVIR(ISYMB)
C
                  NBI         = NBI1 + B
                  NBIAJ       = IT2AM(ISYMAJ,ISYMBI)
     &                        + NT1AM(ISYMAJ) * (NBI - 1) + NAJ
                  NAIB        = KOFF1 + NORBAI * (B - 1)
C
                  EXTBE(NAIB) = EXTBE(NAIB) - TR2E(NBIAJ)
                  EXTBD(NAIB) = EXTBD(NAIB) - TR2D(NBIAJ)
C
  231          CONTINUE
C
  230       CONTINUE
C
         END IF
C
  200 CONTINUE
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_BEXTRACT')
C
      RETURN
      END
