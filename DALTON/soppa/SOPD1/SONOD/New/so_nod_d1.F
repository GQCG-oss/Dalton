C
C  /* Deck so_nod_d1 */
      SUBROUTINE SO_NOD_D1(RINT,D1MPA)
C
C     This is a routine for proving the atomic integral direct
C     SOPPA-D1 programm.
C
C     Mats Simmermacher & Stephan P. A. Sauer, November 2014
C
C     PURPOSE: Nody code to calculate the the D1 matrix element.
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION RINT(10,10,10,10)
      REAL*8 D1MPA
      DIMENSION D1MPA(10000,10000)
      REAL*8 TEMDM
      REAL*8 PRE
      INTEGER NA, NB, NC, ND
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
C     Generate the SOPPA-D1 matrix.
C-----------------------------------------------------------------------
C
      CALL DZERO(D1MPA,100000000)
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
                           NA = NRHF(1) + A
                           NB = NRHF(1) + B
                           NC = NRHF(1) + C
                           ND = NRHF(1) + D
                           NAI = NVIR(1) * (I-1) + A
                           NBJ = NVIR(1) * (J-1) + B
                           NCK = NVIR(1) * (K-1) + C
                           NDL = NVIR(1) * (L-1) + D
C-----------------------------------------------------------------------
C
      IF (NAI .GE. NBJ .AND. NCK .GE. NDL) THEN
C
         NAIBJ = INDEX(NAI,NBJ)
         NCKDL = INDEX(NCK,NDL)
C
         TEMDM =  D1MPA(NAIBJ,NCKDL)
C
C-----------------------------------------------------------------------
C     Prefactor of the matrix element.
C-----------------------------------------------------------------------
C
      IF (A .EQ. B .AND. I .EQ. J) THEN
         PRE = HALF            
      ELSE
         PRE = ONE
      END IF     
C
C-----------------------------------------------------------------------
C     First part of the Soppa-D1 matrix.
C-----------------------------------------------------------------------
C
      IF (A .EQ. C .AND. I .EQ. K) THEN
         TEMDM = TEMDM + PRE * TWO * RINT(NB,J,L,ND)
      END IF
C
      IF (A .EQ. D .AND. I .EQ. L) THEN
         TEMDM = TEMDM + PRE * TWO * RINT(NB,J,K,NC)
      END IF
C
      IF (B .EQ. C .AND. J .EQ. K) THEN
         TEMDM = TEMDM + PRE * TWO * RINT(NA,I,L,ND)
      END IF
C
      IF (B .EQ. D .AND. J .EQ. L) THEN
         TEMDM = TEMDM + PRE * TWO * RINT(NA,I,K,NC)
      END IF
C
C-----------------------------------------------------------------------
C
      IF (A .EQ. C .AND. I .EQ. L) THEN
         TEMDM = TEMDM - PRE * RINT(NB,J,K,ND)
      END IF
C
      IF (A .EQ. D .AND. I .EQ. K) THEN
         TEMDM = TEMDM - PRE * RINT(NB,J,L,NC)
      END IF
C
      IF (B .EQ. C .AND. J .EQ. L) THEN
         TEMDM = TEMDM - PRE * RINT(NA,I,K,ND)
      END IF
C
      IF (B .EQ. D .AND. J .EQ. K) THEN
         TEMDM = TEMDM - PRE * RINT(NA,I,L,NC)
      END IF
C
C-----------------------------------------------------------------------
C     Second part of the Soppa-D1 matrix.
C-----------------------------------------------------------------------
C
      IF (A .EQ. C .AND. I .EQ. K) THEN
         TEMDM = TEMDM - PRE * RINT(NB,ND,L,J)
      END IF
C
      IF (A .EQ. D .AND. I .EQ. L) THEN
         TEMDM = TEMDM - PRE * RINT(NB,NC,K,J)
      END IF
C
      IF (B .EQ. C .AND. J .EQ. K) THEN
         TEMDM = TEMDM - PRE * RINT(NA,ND,L,I)
      END IF
C
      IF (B .EQ. D .AND. J .EQ. L) THEN
         TEMDM = TEMDM - PRE * RINT(NA,NC,K,I)
      END IF
C
C-----------------------------------------------------------------------
C
      IF (A .EQ. C .AND. J .EQ. L) THEN
         TEMDM = TEMDM - PRE * RINT(NB,ND,K,I)
      END IF
C
      IF (A .EQ. D .AND. J .EQ. K) THEN
         TEMDM = TEMDM - PRE * RINT(NB,NC,L,I)
      END IF
C
      IF (B .EQ. C .AND. I .EQ. L) THEN
         TEMDM = TEMDM - PRE * RINT(NA,ND,K,J)
      END IF
C
      IF (B .EQ. D .AND. I .EQ. K) THEN
         TEMDM = TEMDM - PRE * RINT(NA,NC,L,J)
      END IF
C
C-----------------------------------------------------------------------
C     Third part of the Soppa-D1 matrix.
C-----------------------------------------------------------------------
C
      IF (A .EQ. C .AND. B .EQ. D) THEN
         TEMDM = TEMDM + PRE * RINT(K,I,L,J)
      END IF
C
      IF (A .EQ. D .AND. B .EQ. C) THEN
         TEMDM = TEMDM + PRE * RINT(K,J,L,I)
      END IF
C
C-----------------------------------------------------------------------
C     Fourth part of the Soppa-D1 matrix.
C-----------------------------------------------------------------------
C
      IF (I .EQ. K .AND. J .EQ. L) THEN
         TEMDM = TEMDM + PRE * RINT(NA,NC,NB,ND)
      END IF
C
      IF (I .EQ. L .AND. J .EQ. K) THEN
         TEMDM = TEMDM + PRE * RINT(NA,ND,NB,NC)
      END IF
C
C-----------------------------------------------------------------------
C
      D1MPA(NAIBJ,NCKDL) = TEMDM
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
