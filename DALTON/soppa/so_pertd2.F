      SUBROUTINE SO_PERTD2(ISYMTR,FACTOR,
     &                     RDENSIJ,LRDENSIJ,RDENSAB,LRDENSAB,
     &                     T2AM,LT2AM,SOLVEC2,LSOLVEC2)
C
C  Important statements         
      IMPLICIT NONE
C
C Symmetry-offsets in amplitudes.
#include "ccorb.h"
#include "ccsdsym.h"
#include "soppinf.h"
C      
C  Input
      INTEGER, INTENT(IN) ::
     &                ISYMTR, ! Symmetry of solvec2, RDENS
     &                LRDENSIJ, LRDENSAB, LT2AM, LSOLVEC2 
      DOUBLE PRECISION, INTENT(IN) ::
     &      FACTOR,           ! Scaling factor for contributions   
     &      T2AM(LT2AM),      ! Doubles amplitude (paired basis)
     &      SOLVEC2(LSOLVEC2) ! Doubles part of solution
C
C  Output
      DOUBLE PRECISION, INTENT(INOUT) :: 
     &      RDENSIJ(LRDENSIJ),! IJ part of response density
     &      RDENSAB(LRDENSAB) ! AB part of response density
C
C  Locals 
      INTEGER :: INDEX
C     Symmetry of orbital pair      
      INTEGER :: ISYMAB, ISYMBI, ISYMBK, ISYMAI, ISYMAJ, ISYMCJ, ISYMIJ
C     Symmetry of orbital
      INTEGER :: ISYMA, ISYMB, ISYMC, ISYMI, ISYMJ, ISYMK
C     Number of orbitals
      INTEGER :: NOCCI, NOCCJ, NOCCK, NVIRB, NVIRA, NVIRC
C     Double indices      
      INTEGER :: NIJ, NAI, NAJ, NBK, NAB, NBI, NCJ
C     Quadruple indices
      INTEGER :: NAIBK, NAJBK, NAICJ, NBICJ, NAIAI, NAJAI, NAIBI, NBIBI
      INTEGER :: IOFFAB
      REAL*8  :: TMP      
C            
C------------------------------
C     Statement function INDEX.
C------------------------------
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
C-------------------------------------
C     Doubles contribution to RDENSIJ.
C-------------------------------------
C     D_{ij} += - sum_{abk} x^{ab}_{ik} * T^{ab}_{jk}
C               -  sum_{a}  x^{aa}_{ii} * T^{aa}_{ji}
C  
C  Second term only relevant for totally symmetric x  
C      
C  First implementation only, we probably want to use squared 
C  amplitudes and solution vectors
C     
      ISYMIJ = ISYMTR 
      DO ISYMJ = 1, NSYM
C        
         ISYMI = MULD2H(ISYMJ,ISYMIJ)
C         
         DO ISYMA = 1,NSYM
            ISYMAI = MULD2H(ISYMA,ISYMI)
            ISYMAJ = MULD2H(ISYMA,ISYMJ)
            ISYMBK = ISYMAJ ! T is totally symmetric
C
            NVIRA = NVIR(ISYMA)
            NOCCI = NRHF(ISYMI)
            NOCCJ = NRHF(ISYMJ)
C            
            DO ISYMK = 1, NSYM
C               ISYMA = MULD2H(ISYMJ,ISYMAJ)
C               ISYMI = MULD2H(ISYMA,ISYMAI)
               ISYMB = MULD2H(ISYMK,ISYMBK)
C               
               NOCCK = NRHF(ISYMK)
               NVIRB = NVIR(ISYMB)            
C               
               DO I = 1, NRHF(ISYMI)
                  DO J = 1, NRHF(ISYMJ)
C
                     NIJ = IIJDEN(ISYMJ,ISYMI) + NOCCI*(J-1) + I
                     TMP = 0.0D0
C                  
                     DO K =  1, NOCCK
                        DO A = 1, NVIRA
                           NAI = IT1AM(ISYMA,ISYMI) + NVIRA*(I-1) + A
                           NAJ = IT1AM(ISYMA,ISYMJ) + NVIRA*(J-1) + A
C
                           DO B = 1, NVIRB
                              NBK = IT1AM(ISYMB,ISYMK) + NVIRB*(K-1) + B
                              NAIBK = QINDEX(ISYMAI,NAI,ISYMBK,NBK) 
                              NAJBK = QINDEX(ISYMAJ,NAJ,ISYMBK,NBK) 
                              TMP = TMP + SOLVEC2(NAIBK)*T2AM(NAJBK)
                           END DO ! B
                        END DO ! A
                     END DO ! K 
C
C     For totally symmetric perturbations, there is an aditional
C     contribution
C     Note for improvement: Make two versions of the above for
C     ISYMTR = 1 and != 1, then icorperate this into the above inner
C     loop in the first case.                     
                     IF (ISYMIJ.EQ.1.AND.(ISYMK.EQ.ISYMI) ) THEN
                        DO A = 1, NVIRA
                           NAI = IT1AM(ISYMA,ISYMI) + NVIRA*(I-1) + A
                           NAJ = IT1AM(ISYMA,ISYMJ) + NVIRA*(J-1) + A
                           NAIAI = IT2AM(ISYMAI,ISYMAI) + NAI*(NAI+1)/2 
                           NAJAI = IT2AM(ISYMAJ,ISYMAI) + INDEX(NAI,NAJ)
                           TMP = TMP + SOLVEC2(NAIAI)*T2AM(NAJAI)
                        END DO ! A
                     END IF ! Totally symmetric
                     RDENSIJ(NIJ) = RDENSIJ(NIJ) - TMP*FACTOR

                  END DO ! J
               END DO ! I
            END DO ! ISYMK
         END DO ! ISYMA
      END DO ! ISYMJ  
C
C-------------------------------------
C     Doubles contribution to RDENSAB.
C-------------------------------------
C     D_{ab} += sum_{cij} x^{bc}_{ij} * T^{ac}_{ij} 
C              + sum_i x^{bb}_{ii} * T^{ab}_{ii}
C
C  First implementation only, we probably want to use squared 
C  amplitudes and solution vectors
C      
      ISYMAB = ISYMTR   
      DO ISYMB = 1, NSYM
         ISYMA = MULD2H(ISYMB,ISYMAB)
         IOFFAB= IABDEN(ISYMB,ISYMA)
C
         DO ISYMI = 1, NSYM
            ISYMAI = MULD2H(ISYMA,ISYMI)
            ISYMBI = MULD2H(ISYMB,ISYMI)
            ISYMCJ = ISYMAI ! Amplitudes totally symmetric
            NVIRA = NVIR(ISYMA)
            NVIRB = NVIR(ISYMB)
            NOCCI = NRHF(ISYMI)
C            
            DO ISYMJ = 1, NSYM
               ISYMC = MULD2H(ISYMJ,ISYMCJ)
               NVIRC = NVIR(ISYMC)
               NOCCJ = NRHF(ISYMJ)
C               
               DO B = 1, NVIRB
                  DO A = 1, NVIRA
                     TMP = 0.0D0
                     NAB = IOFFAB + NVIRA*(B-1) + A
C                     
                     DO I = 1, NOCCI
                        NAI = IT1AM(ISYMA,ISYMI) + NVIRA*(I-1) + A
                        NBI = IT1AM(ISYMB,ISYMI) + NVIRB*(I-1) + B
C                       
                        DO J = 1, NOCCJ
                           DO C = 1, NVIRC
                              NCJ = IT1AM(ISYMC,ISYMJ) + NVIRC*(J-1) +C
                              
                              NAICJ = QINDEX(ISYMAI,NAI,ISYMCJ,NCJ) 
                              NBICJ = QINDEX(ISYMBI,NBI,ISYMCJ,NCJ)
                              TMP = TMP + SOLVEC2(NBICJ)*T2AM(NAICJ)
                           END DO ! C
                        END DO ! J
                     END DO ! I
C
C     For totally symmetric perturbations, there is an aditional
C     contribution
C     Note for improvement: Make two versions of the above for
C     ISYMTR = 1 and != 1, then incorperate this into the above inner
C     loop in the first case.                     
                     IF (ISYMAB.EQ.1.AND.(ISYMC.EQ.ISYMB) ) THEN
                        DO I = 1, NOCCI
                           NAI = IT1AM(ISYMA,ISYMI) + NVIRA*(I-1) + A
                           NBI = IT1AM(ISYMB,ISYMI) + NVIRA*(I-1) + B
                           NAIBI = IT2AM(ISYMAI,ISYMBI) + INDEX(NAI,NBI)
                           NBIBI = IT2AM(ISYMBI,ISYMBI) + NBI*(NBI+1)/2
                           TMP = TMP + SOLVEC2(NBIBI)*T2AM(NAIBI)
                        END DO ! A
                     END IF ! Totally symmetric
C
                     RDENSAB(NAB) = RDENSAB(NAB) + FACTOR*TMP
                  END DO ! A 
               END DO ! B
            END DO ! ISYMJ
         END DO ! ISYMB
      END DO !ISYMAI
C
      CONTAINS 
         PURE INTEGER FUNCTION QINDEX(ISYM1,I1,ISYM2,I2) 
C        Get the quadruple index from two doubles indicies
         INTEGER, INTENT(IN) :: ISYM1, ISYM2, I1, I2
         IF ( ISYM1 .EQ. ISYM2) THEN
C           Samme symmetry -> triangular distribution            
            QINDEX = IT2AM(ISYM1,ISYM2) + INDEX(I1,I2)
         ELSE IF (ISYM1 .LT. ISYM2) THEN 
C           Square, I1 inner 
            QINDEX = IT2AM(ISYM1,ISYM2) + NT1AM(ISYM1)*(I2-1) + I1
         ELSE 
C           Square, I2 inner
            QINDEX = IT2AM(ISYM2,ISYM1) + NT1AM(ISYM2)*(I1-1) + I2
         END IF

         RETURN
         END FUNCTION
      END 
