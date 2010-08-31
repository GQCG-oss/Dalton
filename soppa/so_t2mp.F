C  /* Deck so_t2mp */
      SUBROUTINE SO_T2MP(T2AM,LT2AM)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, November 1995
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C
C     PURPOSE: Calculate MP2-amplitudes from regular T2-amplitudes.
C              (MP-amplitudes are 4 * Coulomb - 2 * Exchange of
C              T2-amplitudes). The MP2-amplitudes are triangular
C              symmetry packed like the T2-amplitudes and stored
C              in place.
C
C    Alternatively use subroutine: CCSD_TCMEPK
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (FOUR = 4.0D0)
C
      DIMENSION T2AM(LT2AM)
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
      CALL QENTER('SO_T2MP')
C
C---------------------------------------------------
C     Loop over T2-amplitudes for a =< b and i =< j.
C---------------------------------------------------
C
      DO 100 ISYMJ = 1,NSYM
C
         DO 200 ISYMB = 1,NSYM
C
            ISYMBJ = MULD2H(ISYMB,ISYMJ)
            ISYMAI = ISYMBJ
C
            DO 300 J = 1,NRHF(ISYMJ)
C
               DO 400 B = 1,NVIR(ISYMB)
C
                  NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J-1) + B
C
                  DO 500 ISYMI = 1,NSYM
C
                     ISYMA  = MULD2H(ISYMI,ISYMAI)
                     ISYMAJ = MULD2H(ISYMA,ISYMJ)
                     ISYMBI = MULD2H(ISYMB,ISYMI)
C
C---------------------------------------
C                    Ensure that i =< j.
C---------------------------------------
C
                     IF ( ISYMI .LT. ISYMJ) THEN
                        LI = NRHF(ISYMI)
                     ELSE IF ( ISYMI .EQ. ISYMJ) THEN
                        LI = J
                     ELSE
                        LI = 0
                     END IF
C
                     DO 600 I = 1,LI
C
                        NBI   = IT1AM(ISYMB,ISYMI)
     &                        + NVIR(ISYMB)*(I-1) + B
C
C------------------------------------------
C                       Ensure that a =< b.
C------------------------------------------
C
                        IF ( ISYMA .LT. ISYMB ) THEN
                           LA = NVIR(ISYMA)
                        ELSE IF ( ISYMA .EQ. ISYMB ) THEN
                           LA = B
                        ELSE
                           LA = 0
                        ENDIF
C
                        DO 700 A = 1,LA
C
                           NAI   = IT1AM(ISYMA,ISYMI)
     &                           + NVIR(ISYMA)*(I-1) + A
                           NAJ   = IT1AM(ISYMA,ISYMJ)
     &                           + NVIR(ISYMA)*(J-1) + A
                           NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     &                           + INDEX(NAI,NBJ)
                           NAJBI = IT2AM(ISYMAJ,ISYMBI)
     &                           + INDEX(NAJ,NBI)
C
                           T2MP1 = TWO * T2AM(NAIBJ)
     &                           - ONE * T2AM(NAJBI)
                           T2MP2 = TWO * T2AM(NAJBI)
     &                           - ONE * T2AM(NAIBJ)
C
                           T2AM(NAIBJ) = T2MP1
                           T2AM(NAJBI) = T2MP2
C
  700                   CONTINUE
C
  600                CONTINUE
C
  500             CONTINUE
C
  400          CONTINUE
C
  300       CONTINUE
C
  200    CONTINUE
C
  100 CONTINUE
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_T2MP')
C
      RETURN
      END
