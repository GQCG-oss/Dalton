C
      SUBROUTINE SO_PERTD1(MODEL,ISYMTR,FACTOR,
     &                     RDENSIJ,LRDENSIJ,RDENSAB,LRDENSAB,
     &                     RDENSAI,LRDENSAI,
     &                     SOLVEC1,LSOLVEC1,
     &                     DENSIJ,LDENSIJ,
     &                     DENSAB,LDENSAB,DENSAI,LDENSAI,WORK,LWORK)
C
C     Rasmus Faber, November 2015
C
C     PURPOSE: Calculating the Response/perturbed density matrix for
C              corresponding to the given solution/eigen vector.
C              This should ease the calculation of the final property
C
C     INPUT:
C        MODEL    Determine which terms to include
C        ISYMTR   Symmetry of the trial-vector
C        FACTOR   Factor with which to scale all contributions,
C                 use to correct for the fact that the solution vector
C                 is normalized in SO_OPTVEC
C        SOLVEC1  One-particle solution vector
C        SOLVEC2  Two-particle solution vector (IA)
C        T2AM     T2 amplitudes im complementary basis
C        DENSAB
C        DENSIJ   Second order densities (IA)
C        DENSAI
C
C     OUTPUT:
C        RDENSAB
C        RDENSIJ  Response densities, must be initialized outside (IA)
C        RDENSAI
C
C     Note for the following: DENS?? arrays contain the density with a
C        factor of one half, therefore the difference between comments
C        and code
C
      implicit none
C Get MULD2H (group multiplication table)
#include "ccorb.h"
C Symmetry-offsets in amplitudes
#include "ccsdsym.h"
C Symmetry-offsets in densities
#include "soppinf.h"
C
C  Arguments
C
      CHARACTER*5, INTENT(IN) :: MODEL
C
      INTEGER,INTENT(IN) ::  ISYMTR,
     &                       LDENSIJ, LDENSAB, LDENSAI, LWORK,
     &                       LSOLVEC1,
     &                       LRDENSAB, LRDENSIJ, LRDENSAI

      DOUBLE PRECISION,INTENT(IN) ::
     &                            SOLVEC1(LSOLVEC1),
     &                            DENSIJ(LDENSIJ), DENSAB(LDENSAB),
     &                            DENSAI(LDENSAI), FACTOR
C
      DOUBLE PRECISION, INTENT(INOUT) ::  RDENSAB(LRDENSAB),
     &                                    RDENSIJ(LRDENSIJ),
     &                                    RDENSAI(LRDENSAI),
     &                                    WORK(LWORK)
C

C
C     Local variables
      DOUBLE PRECISION THISFACT
      LOGICAL VALID
C
      INTEGER :: ISYMA, ISYMB, ISYMI, ISYMJ   ! Irrep of index
      INTEGER :: NVIRA, NVIRB, NOCCI          ! Size of orbital groups
      INTEGER :: KOFFSOL, KOFFDEN, KOFFOUT    ! Offsets in arrays
      DOUBLE PRECISION, PARAMETER :: ONE = 1.0D0
      DOUBLE PRECISION, PARAMETER :: SQRT2 = SQRT(2.0D0)
C
      CALL QENTER('SO_PERTD1')
C
C     Check arguments here?
C
C---------------------------------------------------
C     Calculate AI part of response density
C---------------------------------------------------
C     Please note that though the terms are in AI order,
C     actually the equations are for the the IA part of
C     perturbed density is calculated, meaning that we have
C     a factor of -1 on AI term for imagninary properties
C
C     Term one:
C     D_{ai}* <= \sqrt{2} *x_{ai}
C
      THISFACT = FACTOR*SQRT2
      CALL DAXPY(LRDENSAI,THISFACT,SOLVEC1,1,RDENSAI,1)
C
C     This is all for RPA
      IF (MODEL.EQ.'AORPA') GOTO 1010
C
C     Term two:
C                   1
C     D_{ai}* <= -------- \sum_j x_{aj} * \rho_{ij}
C                \sqrt{2}
      THISFACT = FACTOR*SQRT2 ! \rho is stored as half?
      DO ISYMI = 1, NSYM
C
         ISYMA = MULD2H(ISYMI,ISYMTR)
         ISYMJ = ISYMI ! \rho is totally symetric
         KOFFSOL = IT1AM(ISYMA,ISYMJ) + 1
         KOFFOUT = IT1AM(ISYMA,ISYMI) + 1
C         KOFFDEN = IIJDEN(ISYMI,ISYMJ) + 1
         KOFFDEN = IIJDEN(ISYMJ,ISYMI) + 1
         NVIRA = MAX(NVIR(ISYMA),1)
         NOCCI = MAX(NRHF(ISYMI),1)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMJ),
     &              THISFACT,SOLVEC1(KOFFSOL),NVIRA,DENSIJ(KOFFDEN),
     &              NOCCI,ONE,RDENSAI(KOFFOUT),NVIRA)
      END DO
C
C     Term three:
C                  - 1
C     D_{ai}* <= -------- \sum_b x_{bi} * \rho_{ba}
C                \sqrt{2}
      THISFACT = -FACTOR*SQRT2 ! \rho is stored as half?
      DO ISYMI = 1, NSYM
C
         ISYMA = MULD2H(ISYMI,ISYMTR)
         ISYMB = ISYMA
         NVIRA = MAX(NVIR(ISYMA),1)
         NOCCI = MAX(NRHF(ISYMI),1)
         NVIRB = MAX(NVIR(ISYMB),1)
         KOFFSOL = IT1AM(ISYMB,ISYMI) + 1
         KOFFOUT = IT1AM(ISYMA,ISYMI) + 1
C         KOFFDEN = IABDEN(ISYMB,ISYMA) + 1
         KOFFDEN = IABDEN(ISYMA,ISYMB) + 1
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMB),
     &              THISFACT,DENSAB(KOFFDEN),NVIRB,
     &              SOLVEC1(KOFFSOL),NVIRB,ONE,
     &              RDENSAI(KOFFOUT),NVIRA)
      END DO
C
C--------------------------------------------------------
C     Calculate one-electron IJ part of response density.
C--------------------------------------------------------
C                 - 1
C     D_{ij} <= -------- \sum_b x_{bi} * \rho_{bj}
C               \sqrt{2}
      THISFACT = -FACTOR*SQRT2
      DO ISYMI = 1,NSYM
         ISYMB = MULD2H(ISYMI,ISYMTR)
         ISYMJ = ISYMB ! \rho totally symmetric
         NOCCI = MAX(NRHF(ISYMI),1)
         NVIRB = MAX(NVIR(ISYMB),1)
         KOFFSOL = IT1AM(ISYMB,ISYMI) + 1
C         KOFFOUT = IIJDEN(ISYMI,ISYMJ) + 1
         KOFFOUT = IIJDEN(ISYMJ,ISYMI) + 1
         KOFFDEN = IAIDEN(ISYMB,ISYMJ) + 1
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMB),
     &              THISFACT,SOLVEC1(KOFFSOL),NVIRB,
     &              DENSAI(KOFFDEN),NVIRB,ONE,
     &              RDENSIJ(KOFFOUT),NOCCI)
      END DO
C
C--------------------------------------------------------
C     Calculate one-electron AB part of response density.
C--------------------------------------------------------
C                  1
C     D_{ab} <= -------- \sum_i x_{bi} * \rho_{ai}
C               \sqrt{2}
C
      THISFACT = FACTOR*SQRT2
      KOFFOUT = 1
      DO ISYMB = 1, NSYM
         ISYMI = MULD2H(ISYMB,ISYMTR)
         ISYMA = ISYMI
         NVIRA = MAX(NVIR(ISYMA),1)
         NVIRB = MAX(NVIR(ISYMB),1)
         KOFFSOL = IT1AM(ISYMB,ISYMI) + 1
C         KOFFDEN = IAIDEN(ISYMA,ISYMI) + 1
         KOFFDEN = IT1AM(ISYMA,ISYMI) + 1
C         KOFFOUT = IABDEN(ISYMA,ISYMB)

         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMI),
     &              THISFACT,DENSAI(KOFFDEN),NVIRA,
     &              SOLVEC1(KOFFSOL),NVIRB,ONE,
     &              RDENSAB(KOFFOUT),NVIRA)
         KOFFOUT = KOFFOUT + NVIR(ISYMA)*NVIR(ISYMB)
      END DO
C
C     All continue here
1010  CALL QEXIT('SO_PERTD1')

      RETURN
      END SUBROUTINE
