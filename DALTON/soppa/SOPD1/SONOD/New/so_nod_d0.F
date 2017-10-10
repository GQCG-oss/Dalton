C
C  /* Deck so_nod_d0 */
      SUBROUTINE SO_NOD_D0(FOCKD,LFOCKD,DMZO)
C
C     This is a routine for proving the atomic integral direct
C     SOPPA-D1 programm.
C
C     Mats Simmermacher & Stephan P. A. Sauer, November 2014
C
C     PURPOSE: Nody code to calculate the the A(0,1) matrix element.
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION FOCKD(LFOCKD)
      REAL*8 DMZO
      DIMENSION DMZO(10000,10000)
      REAL*8 TEMDM
      INTEGER KOFFA, KOFFB
      INTEGER KOFFI, KOFFJ
      INTEGER NAI, NBJ, NAIBJ
      INTEGER NCK, NDL, NCKDL
      INTEGER X, Y
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
C-----------------------------------------------------------------------
C     Statement function INDEX.
C-----------------------------------------------------------------------
C
      INDEX(X,Y) = MAX(X,Y)*(MAX(X,Y) - 3)/2 + X + Y
C
C-----------------------------------------------------------------------
C     Generate the SOPPA-D(0) matrix.
C-----------------------------------------------------------------------
C
      CALL DZERO(DMZO,100000000)
C
      DO 100 I = 1,NRHF(1)
         DO 200 A = 1,NVIR(1)
            DO 300 J = 1,NRHF(1)
               DO 400 B = 1,NVIR(1)
                  DO 500 C = 1,NVIR(1)
                     DO 600 K = 1,NRHF(1)
                        DO 700 D = 1,NVIR(1)
                           DO 800 L = 1,NRHF(1)
C
                              KOFFA = IVIR(1) + A
                              KOFFB = IVIR(1) + B
                              KOFFI = IRHF(1) + I
                              KOFFJ = IRHF(1) + J
                              NAI = NVIR(1) * (I-1) + A
                              NBJ = NVIR(1) * (J-1) + B
                              NCK = NVIR(1) * (K-1) + C
                              NDL = NVIR(1) * (L-1) + D
C
C-----------------------------------------------------------------------
C
      IF (NAI .GE. NBJ .AND. NCK .GE. NDL) THEN
C
         NAIBJ = INDEX(NAI,NBJ)
         NCKDL = INDEX(NCK,NDL)
C
         TEMDM = DMZO(NAIBJ,NCKDL)
C
C-----------------------------------------------------------------------
C
         IF (A .EQ. C .AND. B .EQ. D .AND. I .EQ. K .AND. J .EQ. L) THEN
            TEMDM = TEMDM + FOCKD(KOFFA) + FOCKD(KOFFB) - FOCKD(KOFFI)
     &                    - FOCKD(KOFFJ)
         END IF
C
         IF (A .EQ. D .AND. B .EQ. C .AND. I .EQ. L .AND. J .EQ. K) THEN
            TEMDM = TEMDM + FOCKD(KOFFA) + FOCKD(KOFFB) - FOCKD(KOFFI)
     &                    - FOCKD(KOFFJ)
         END IF
C
         IF (A .EQ. B .AND. I .EQ. J) THEN
            TEMDM = HALF * TEMDM
         END IF
C
C-----------------------------------------------------------------------
C
         DMZO(NAIBJ,NCKDL) = TEMDM
C
      END IF
C
C-----------------------------------------------------------------------
C
  800                      CONTINUE
  700                   CONTINUE
  600                CONTINUE
  500             CONTINUE
  400          CONTINUE
  300       CONTINUE
  200    CONTINUE
  100 CONTINUE  
C
C-----------------------------------------------------------------------
C
      RETURN
      END
