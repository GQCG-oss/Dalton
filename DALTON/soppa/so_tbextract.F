C
C  /* Deck so_tbextt */
      SUBROUTINE SO_TBEXTRACT(ITYPE,EXTBEA,LEXTBEA,EXTBEB,LEXTBEB,
     &                        EXTBDA,LEXTBDA,EXTBDB,LEXTBDB,
     &                        TR2E,LTR2E,TR2D,LTR2D,
     &                        ISYMA,ISYMI,ISYMB,ISYMJ,J)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Stephan P. A. Sauer: 25.11.2010: based on the SO_TBEXTRACT routine
C
C     Pi 29.03.16: Copied from dalton20aosoppa
C
C     PURPOSE: Extract b(a i b j) +/- b(b i a j) for given j
C              in EXTBE (and EXTBD) for excitations (and de-excitations)
C              respectively for the three tiplet trial vectors.
C
#include "implicit.h"
#include "priunit.h"
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION EXTBEA(LEXTBEA), EXTBEB(LEXTBEB)
      DIMENSION EXTBDA(LEXTBDA), EXTBDB(LEXTBDB)
      DIMENSION TR2E(LTR2E), TR2D(LTR2D)
C
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "infsop.h"
#include "soppinf.h"
C
C------------------------------
C     Statement function INDEX.
C------------------------------
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
      INDLT(I,J) = (MAX(I,J) - 1)*(MAX(I,J) - 2)/2 + MIN(I,J)
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_TBEXTRACT')
C
C-------------------------
C     Calculate constants.
C-------------------------
C
      SQ2  = DSQRT(TWO)
      SQ2I = ONE / SQ2
C
C-----------------------------------------------------------------
C     Clean the arrays for the squared triplet 2p-2h trial vectors
C-----------------------------------------------------------------
C
      CALL DZERO(EXTBEA,LEXTBEA)
      CALL DZERO(EXTBDA,LEXTBDA)
      IF (ITYPE .NE. 1) THEN
         CALL DZERO(EXTBEB,LEXTBEB)
         CALL DZERO(EXTBDB,LEXTBDB)
      END IF
C
C----------------------------------------------------------------
C     Extract b(a i b j) and b(b i a j) for the three triplet
C     trial vectors for excitations and de-excitations
C
C     In order to combine the four summations a "fully squared" version
C     $\underline{^{T1,E/D}b}~^{2p2h}_{aibj}$ of the T1 2p2h trial vector
C     $^{T1,E/D}b^{2p2h}_{aibj}$ is defined as
C $C\underline{^{T1,E/D}b}~^{2p2h}_{aibj} = ^{T1,E/D}b^{2p2h}_{aibj}$
C   for $a > b$  and $i > j$
C\underline{^{T1,E/D}b}~^{2p2h}_{aibj} &=& - ^{T1,E/D}b^{2p2h}_{biaj}
C~~~~~\textrm{for}~ a < b  ~\textrm{and} ~i > j $
C\underline{^{T1,E/D}b}~^{2p2h}_{aibj} &=& - ^{T1,E/D}b^{2p2h}_{ajbi}
C~~~~~\textrm{for}~ a > b  ~\textrm{and} ~i < j $
C\underline{^{T1,E/D}b}~^{2p2h}_{aibj} &=& ^{T1,E/D}b^{2p2h}_{bjai}
C~~~~~~~\textrm{for}~ a < b  ~\textrm{and} ~i < j $
C\underline{^{T1,E/D}b}~^{2p2h}_{aibj} &=& 0
C~~~~~~~~~~~~~~~~~~~~~\textrm{for}~ a = b  ~\textrm{or} ~i = j $
C
C     with the symmetry properties
C     $
C\underline{^{T1,E/D}b}~^{2p2h}_{aibj} &=&
C\underline{^{T1,E/D}b}~^{2p2h}_{bjai}
C\label{eq.square.T1b2p2h.sym} \\
C\underline{^{T1,E/D}b}~^{2p2h}_{aibj} &=& -
C\underline{^{T1,E/D}b}~^{2p2h}_{ajbi}  \nonumber \\
C\underline{^{T1,E/D}b}~^{2p2h}_{aibj} &=& -
C\underline{^{T1,E/D}b}~^{2p2h}_{biaj}  $
C
C
C     and a "fully squared" version  $\underline{^{T2,E/D}b}~^{2p2h}_{aibj}$
C     of the T2 2p2h trial vector $^{T2,E/D}b^{2p2h}_{aibj}$ is defined as
C\underline{^{T2,E/D}b}~^{2p2h}_{aibj} &=& ^{T2,E/D}b^{2p2h}_{aibj}
C~~~~~~~\textrm{for}~ a > b  ~\textrm{and} ~i > j \label{eq.square.T2b2p2h} \\
C\underline{^{T2,E/D}b}~^{2p2h}_{aibj} &=& 2~ ^{T2,E/D}b^{2p2h}_{aibj}
C~~~~~~~\textrm{for}~ a = b  ~\textrm{and} ~i > j \nonumber \\
C\underline{^{T2,E/D}b}~^{2p2h}_{aibj} &=& ^{T2,E/D}b^{2p2h}_{biaj}
C~~~~~\textrm{for}~ a < b  ~\textrm{and} ~i > j \nonumber \\
C\underline{^{T2,E/D}b}~^{2p2h}_{aibj} &=& - ^{T2,E/D}b^{2p2h}_{ajbi}
C~~~~~\textrm{for}~ a > b  ~\textrm{and} ~i < j \nonumber \\
C\underline{^{T2,E/D}b}~^{2p2h}_{aibj} &=& - 2~^{T2,E/D}b^{2p2h}_{ajbi}
C~~~~~\textrm{for}~ a = b  ~\textrm{and} ~i < j \nonumber \\
C\underline{^{T2,E/D}b}~^{2p2h}_{aibj} &=& - ^{T2,E/D}b^{2p2h}_{bjai}
C~~~~~~~\textrm{for}~ a < b  ~\textrm{and} ~i < j \nonumber \\
C\underline{^{T2,E/D}b}~^{2p2h}_{aibj} &=& 0
C~~~~~~~~~~~~~~~~~~~~~\textrm{for}~ i = j \nonumber
C
Cwith the symmetry properties
C\underline{^{T2,E/D}b}~^{2p2h}_{aibj} &=& -
C\underline{^{T2,E/D}b}~^{2p2h}_{bjai}
C\label{eq.square.T2b2p2h.sym} \\
C\underline{^{T2,E/D}b}~^{2p2h}_{aibj} &=& -
C\underline{^{T2,E/D}b}~^{2p2h}_{ajbi}  \nonumber \\
C\underline{^{T2,E/D}b}~^{2p2h}_{aibj} &=&
C\underline{^{T2,E/D}b}~^{2p2h}_{biaj}  \nonumber
C
C
C     and a "fully squared" version  $\underline{^{T3,E/D}b}~^{2p2h}_{aibj}$
C     of the T3 2p2h trial vector $^{T3,E/D}b^{2p2h}_{aibj}$ is defined as
C\underline{^{T3,E/D}b}~^{2p2h}_{aibj} &=& ^{T3,E/D}b^{2p2h}_{aibj}
C~~~~~~~\textrm{for}~ a > b  ~\textrm{and} ~i > j \label{eq.square.T3b2p2h} \\
C\underline{^{T3,E/D}b}~^{2p2h}_{aibj} &=& 2~ ^{T3,E/D}b^{2p2h}_{aibj}
C~~~~~\textrm{for}~ a > b  ~\textrm{and} ~i = j \nonumber \\
C\underline{^{T3,E/D}b}~^{2p2h}_{aibj} &=& ^{T3,E/D}b^{2p2h}_{ajbi}
C~~~~~~~\textrm{for}~ a > b  ~\textrm{and} ~i < j \nonumber \\
C\underline{^{T3,E/D}b}~^{2p2h}_{aibj} &=& - ^{T3,E/D}b^{2p2h}_{biaj}
C~~~~~\textrm{for}~ a < b  ~\textrm{and} ~i > j \nonumber \\
C\underline{^{T3,E/D}b}~^{2p2h}_{aibj} &=& - 2~^{T3,E/D}b^{2p2h}_{ajbi}
C~~~\textrm{for}~ a < b  ~\textrm{and} ~i = j \nonumber \\
C\underline{^{T3,E/D}b}~^{2p2h}_{aibj} &=& - ^{T3,E/D}b^{2p2h}_{bjai}
C~~~~~\textrm{for}~ a < b  ~\textrm{and} ~i < j \nonumber \\
C\underline{^{T3,E/D}b}~^{2p2h}_{aibj} &=& 0
C~~~~~~~~~~~~~~~~~~~~~\textrm{for}~ a = b \nonumber
C
Cwith the symmetry properties
C\underline{^{T3,E/D}b}~^{2p2h}_{aibj} &=& -
C\underline{^{T3,E/D}b}~^{2p2h}_{bjai}
C\label{eq.square.T3b2p2h.sym} \\
C\underline{^{T3,E/D}b}~^{2p2h}_{aibj} &=&
C\underline{^{T3,E/D}b}~^{2p2h}_{ajbi}  \nonumber \\
C\underline{^{T3,E/D}b}~^{2p2h}_{aibj} &=& -
C\underline{^{T3,E/D}b}~^{2p2h}_{biaj}  \nonumber
C
C     and put the following combinations of the "fully squared" 2p2h trial
Cvector $\underline{^{T1/2/3,E/D}b}~^{2p2h}_{aibj}$
C\underline{^{T1+,E/D}\widetilde{b}}~^{2p2h}_{aibj} &=&
C~\underline{^{T1,E/D}b}~^{2p2h}_{aibj}
C+  ~\underline{^{T1,E/D}b}~^{2p2h}_{biaj} \\
C\underline{^{T2+,E/D}\widetilde{b}}~^{2p2h}_{aibj} &=&
C~\underline{^{T2,E/D}b}~^{2p2h}_{aibj}
C+  ~\underline{^{T2,E/D}b}~^{2p2h}_{biaj} \\
C\underline{^{T3+,E/D}\widetilde{b}}~^{2p2h}_{aibj} &=&
C~\underline{^{T3,E/D}b}~^{2p2h}_{aibj}
C+  ~\underline{^{T3,E/D}b}~^{2p2h}_{biaj} \\
C   in EXTBEA and EXTBDA
C
C   and
C\underline{^{T1+,E/D}\widetilde{b}}~^{2p2h}_{aibj} &=&
C~\underline{^{T1,E/D}b}~^{2p2h}_{aibj}
C+  ~\underline{^{T1,E/D}b}~^{2p2h}_{biaj} \\
C   $\underline{^{T2-,E/D}\widetilde{b}}~^{2p2h}_{aibj} =
C                \underline{^{T2,E/D}b}~^{2p2h}_{aibj}
C              -  \underline{^{T2,E/D}b}~^{2p2h}_{biaj}$
C\underline{^{T3-,E/D}\widetilde{b}}~^{2p2h}_{aibj} &=&
C~\underline{^{T3,E/D}b}~^{2p2h}_{aibj}
C-  ~\underline{^{T3,E/D}b}~^{2p2h}_{biaj}
C     in EXTBEB and EXTBDB.
C----------------------------------------------------------------
C
      ISYMAI = MULD2H(ISYMA,ISYMI)
      ISYMBJ = MULD2H(ISYMB,ISYMJ)
      ISYMAB = MULD2H(ISYMA,ISYMB)
      ISYMIJ = MULD2H(ISYMI,ISYMJ)
C
      IIJT   = ITOO(MAX(ISYMI,ISYMJ),MIN(ISYMI,ISYMJ))
      IIJS   = ISOO(MAX(ISYMI,ISYMJ),MIN(ISYMI,ISYMJ))
      IIJ1   = IIJT + J
      IIJ3   = IIJS + J
C
      IABT   = ITVV(MAX(ISYMA,ISYMB),MIN(ISYMA,ISYMB))
      IABS   = ISVV(MAX(ISYMA,ISYMB),MIN(ISYMA,ISYMB))
C
      DO 100 B = 1,NVIR(ISYMB)
C
         KOFF1 = NVIR(ISYMA) * NRHF(ISYMI) * (B - 1)
C
         IAB1   = IABT + B
         IAB2   = IABS + B
C
         DO 110 I = 1,NRHF(ISYMI)
C
            KOFF2 = KOFF1 + NVIR(ISYMA) * (I - 1)
C
            IF (ITYPE .EQ. 1) THEN
C
               IF (ISYMI .EQ. ISYMJ) THEN
                  NIJ1  = IIJT + INDLT(I,J)
                  IF (I.EQ.J) THEN
                     SGNIJ1 = ZERO
                  ELSE IF (I.GT.J) THEN
                     SGNIJ1 = SQ2
                  ELSE
                     SGNIJ1 = -SQ2
                  END IF
               ELSE IF (ISYMI .GT. ISYMJ) THEN
                  KJOFF = NRHF(ISYMJ) * (I - 1)
                  NIJ1  = IIJ1 + KJOFF
                  SGNIJ1 = SQ2
               ELSE IF (ISYMI .LT. ISYMJ) THEN
                  KJOFF = NRHF(ISYMI) * (J - 1)
                  NIJ1  = IIJT + I + KJOFF
                  SGNIJ1 = -SQ2
               END IF
C
               IABIJ1 = IT2AMT1(ISYMIJ,ISYMAB)
     &                  + NTVV(ISYMAB) * (NIJ1 - 1)
C
               DO 120 A = 1,NVIR(ISYMA)
C
                  IF (ISYMA .EQ. ISYMB) THEN
                     NAB1  = IABT + INDLT(A,B)
                     IF(A.EQ.B) THEN
                        SGN1 = ZERO
                     ELSE IF(A.GT.B) THEN
                        SGN1 = SGNIJ1
                     ELSE
                        SGN1 = -SGNIJ1
                     END IF
                  ELSE IF (ISYMA .GT. ISYMB) THEN
                     KBOFF = NVIR(ISYMB) * (A - 1)
                     NAB1  = IAB1 + KBOFF
                     SGN1 = SGNIJ1
                  ELSE IF (ISYMA .LT. ISYMB) THEN
                     KBOFF = NVIR(ISYMA) * (B - 1)
                     NAB1  = IABT + A + KBOFF
                     SGN1 = -SGNIJ1
                  END IF
C
                  NABIJ1 = IABIJ1 + NAB1
C
                  NAIB   = KOFF2 + A
C
                  EXTBEA(NAIB) = EXTBEA(NAIB) + TWO*SGN1 * TR2E(NABIJ1)
                  EXTBDA(NAIB) = EXTBDA(NAIB) + TWO*SGN1 * TR2D(NABIJ1)
C
 120           CONTINUE
C
            ELSE IF (ITYPE .EQ. 2) THEN
C
               IF (ISYMI .EQ. ISYMJ) THEN
                  NIJ2  = IIJT + INDLT(I,J)
                  IF (I.EQ.J) THEN
                     SGNIJ2 = ZERO
                  ELSE IF (I.GT.J) THEN
                     SGNIJ2 = ONE
                  ELSE
                     SGNIJ2 = -ONE
                  END IF
               ELSE IF (ISYMI .GT. ISYMJ) THEN
                  KJOFF = NRHF(ISYMJ) * (I - 1)
                  NIJ2  = IIJ1 + KJOFF
                  SGNIJ2 = ONE
               ELSE IF (ISYMI .LT. ISYMJ) THEN
                  KJOFF = NRHF(ISYMI) * (J - 1)
                  NIJ2  = IIJT + I + KJOFF
                  SGNIJ2 = -ONE
               END IF
C
               IABIJ2 = IT2AMT2(ISYMIJ,ISYMAB)
     &                  + NSVV(ISYMAB) * (NIJ2 - 1)
C
               DO 130 A = 1,NVIR(ISYMA)
C
                  IF (ISYMA .EQ. ISYMB) THEN
                     NAB2  = IABS + INDEX(A,B)
                     IF (A.EQ.B) THEN
                        SGN2 = SQ2*SGNIJ2
                     ELSE IF (A.GT.B) THEN
                        SGN2 = SGNIJ2
                     ELSE
                        SGN2 = SGNIJ2
                     END IF
                  ELSE IF (ISYMA .GT. ISYMB) THEN
                     KBOFF = NVIR(ISYMB) * (A - 1)
                     NAB2  = IAB2 + KBOFF
                     SGN2 = SGNIJ2
                  ELSE IF (ISYMA .LT. ISYMB) THEN
                     KBOFF = NVIR(ISYMA) * (B - 1)
                     NAB2  = IABS + A + KBOFF
                     SGN2 = SGNIJ2
                  END IF
C
                  NABIJ2 = IABIJ2 + NAB2
                  NAIB   = KOFF2 + A
C
                  EXTBEA(NAIB) = EXTBEA(NAIB) - TWO*SGN2 * TR2E(NABIJ2)
                  EXTBDA(NAIB) = EXTBDA(NAIB) - TWO*SGN2 * TR2D(NABIJ2)
                  EXTBEB(NAIB) = EXTBEB(NAIB) - TWO*SGN2 * TR2E(NABIJ2)
                  EXTBDB(NAIB) = EXTBDB(NAIB) - TWO*SGN2 * TR2D(NABIJ2)
C
 130           CONTINUE
C
            ELSE
C
               IF (ISYMI .EQ. ISYMJ) THEN
                   NIJ3  = IIJS + INDEX(I,J)
                  IF (I.EQ.J) THEN
                     SGNIJ3 = SQ2
                  ELSE IF (I.GT.J) THEN
                     SGNIJ3 = ONE
                  ELSE
                     SGNIJ3 = ONE
                  END IF
               ELSE IF (ISYMI .GT. ISYMJ) THEN
                  KJOFF = NRHF(ISYMJ) * (I - 1)
                  NIJ3  = IIJ3 + KJOFF
                  SGNIJ3 = ONE
               ELSE IF (ISYMI .LT. ISYMJ) THEN
                  KJOFF = NRHF(ISYMI) * (J - 1)
                  NIJ3  = IIJS + I + KJOFF
                  SGNIJ3 = ONE
               END IF
C
               IABIJ3 = IT2AMT3(ISYMIJ,ISYMAB)
     &                  + NTVV(ISYMAB) * (NIJ3 - 1)
C
               DO 140 A = 1,NVIR(ISYMA)
C
                  IF (ISYMA .EQ. ISYMB) THEN
                     NAB3  = IABT + INDLT(A,B)
                     IF(A.EQ.B) THEN
                        SGN3 = ZERO
                     ELSE IF(A.GT.B) THEN
                        SGN3 = SGNIJ3
                     ELSE
                        SGN3 = -SGNIJ3
                     END IF
                  ELSE IF (ISYMA .GT. ISYMB) THEN
                     KBOFF = NVIR(ISYMB) * (A - 1)
                     NAB3  = IAB1 + KBOFF
                     SGN3 = SGNIJ3
                  ELSE IF (ISYMA .LT. ISYMB) THEN
                     KBOFF = NVIR(ISYMA) * (B - 1)
                     NAB3  = IABT + A + KBOFF
                     SGN3 = -SGNIJ3
                  END IF
C
                  NABIJ3 = IABIJ3 + NAB3
                  NAIB   = KOFF2 + A
C
                  EXTBEA(NAIB) = EXTBEA(NAIB) - TWO*SGN3 * TR2E(NABIJ3)
                  EXTBDA(NAIB) = EXTBDA(NAIB) - TWO*SGN3 * TR2D(NABIJ3)
                  EXTBEB(NAIB) = EXTBEB(NAIB) - TWO*SGN3 * TR2E(NABIJ3)
                  EXTBDB(NAIB) = EXTBDB(NAIB) - TWO*SGN3 * TR2D(NABIJ3)
C
 140           CONTINUE
C
            END IF
C
 110     CONTINUE
C
 100  CONTINUE
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_TBEXTRACT')
C
      RETURN
      END
