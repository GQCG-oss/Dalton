      SUBROUTINE SO_PERTD2(ISYMTR,FACTOR,
     &                     RDENSIJ,LRDENSIJ,RDENSAB,LRDENSAB,
     &                     T2AM,LT2AM,SOLVEC2,LSOLVEC2,
     &                     WORK,LWORK)
C
C  Important statement
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
     &                LRDENSIJ, LRDENSAB, LT2AM, LSOLVEC2, LWORK
      DOUBLE PRECISION, INTENT(IN) ::
     &      FACTOR,           ! Scaling factor for contributions
     &      T2AM(LT2AM),      ! Doubles amplitude (paired basis)
     &      SOLVEC2(LSOLVEC2) ! Doubles part of solution
C
C  Output
      DOUBLE PRECISION, INTENT(INOUT) ::
     &      RDENSIJ(LRDENSIJ),! IJ part of response density
     &      RDENSAB(LRDENSAB),! AB part of response density
     &      WORK(LWORK)
C
C  Locals
C     Symmetry of orbital pair
      INTEGER :: ISYMBK, ISYMAI, ISYMAJ, ISYMCJ
C     Symmetry of orbital
      INTEGER :: ISYMA, ISYMB, ISYMC, ISYMI, ISYMJ, ISYM
C     OFFSETS
      INTEGER :: IOFFAB, IOFFT, IOFFX
C     Memory Locations
      INTEGER :: KT2SQ, KX2SQ, KEND1
      INTEGER :: IPOSD
C     Lengths
      INTEGER :: LT2SQ, LX2SQ, LD
      DOUBLE PRECISION, PARAMETER ::  ONE = 1.0D0, ONEM = -1.0D0
C
C     Loop over symmetry of right index pair
C
      DO ISYMAJ = 1, NSYM
         ISYMBK = ISYMAJ ! T is totally symmetric
         ISYMAI = MULD2H(ISYMBK,ISYMTR) ! X2 has symmetry ISYMTR
C
C     Allocate space for this block
C
         LT2SQ = NT1AM(ISYMAJ)**2
         LX2SQ = NT1AM(ISYMAI)*NT1AM(ISYMBK)
         KT2SQ = 1
         KX2SQ = KT2SQ + LT2SQ
         KEND1 = KX2SQ + LX2SQ
         IF (KEND1.GT.LWORK) CALL STOPIT('SO_PERTD2.1',' ',KEND1,LWORK)
C        Check size...
C        Square up the blocks
C        T2( (ai)<(bj) ) -> T2( bj , ai)
         CALL SQUARE_TOTSYM(WORK(KT2SQ),T2AM,ISYMAJ)
C        X2( (ai)<(bj) ) -> X2( bj , ai), Diagonal terms scaled by two
         CALL SQUAREX(WORK(KX2SQ),SOLVEC2,ISYMAI,ISYMTR)
C
C-------------------------------------
C     Doubles contribution to RDENSIJ.
C-------------------------------------
C     D_{ij} += - sum_{abk} x^{ab}_{ik} * T^{ab}_{jk}
C               -  sum_{a}  x^{aa}_{ii} * T^{aa}_{ji}
C
C  Second term is taken care off by scaling x^{aa}_{ii} by 2 in
C  the above transpositions
         DO ISYMJ = 1, NSYM
            ISYMA = MULD2H(ISYMJ,ISYMAJ)
            ISYMI = MULD2H(ISYMJ,ISYMTR)
            IOFFT = NT1AM(ISYMBK)*IT1AM(ISYMA,ISYMJ)
            IOFFX = NT1AM(ISYMBK)*IT1AM(ISYMA,ISYMI)
            LD = MAX(NT1AM(ISYMBK)*NVIR(ISYMA),1)
            IPOSD = IIJDEN(ISYMJ,ISYMI) + 1
            CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),
     &                 NT1AM(ISYMBK)*NVIR(ISYMA),ONEM*FACTOR,
     &                 WORK(KX2SQ+IOFFX),LD,
     &                 WORK(KT2SQ+IOFFT),LD,
     &                 ONE,RDENSIJ(IPOSD),MAX(1,NRHF(ISYMI)))
         END DO
C
C-------------------------------------
C     Doubles contribution to RDENSAB.
C-------------------------------------
C     D_{ab} += sum_{cij} x^{bc}_{ij} * T^{ac}_{ij}
C              + sum_i x^{bb}_{ii} * T^{ab}_{ii}
C
C  Second term is taken care off by scaling x^{aa}_{ii} by 2 in
C  the above transpositions
         DO ISYMB = 1, NSYM
            ISYMI = MULD2H(ISYMB,ISYMAI) ! SYM
            ISYMA = MULD2H(ISYMB,ISYMTR)
            ISYMCJ = ISYMBK
            IOFFT = NT1AM(ISYMCJ)*IT1AM(ISYMA,ISYMI)
            IOFFX = NT1AM(ISYMCJ)*IT1AM(ISYMB,ISYMI)
            IPOSD = IABDEN(ISYMB,ISYMA)+1
            DO I = 1, NRHF(ISYMI)
               CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMB),
     &                    NT1AM(ISYMCJ),FACTOR,
     &                    WORK(KT2SQ+IOFFT),MAX(1,NT1AM(ISYMCJ)),
     &                    WORK(KX2SQ+IOFFX),MAX(1,NT1AM(ISYMCJ)),
     &                    ONE,
     &                    RDENSAB(IPOSD),MAX(1,NVIR(ISYMA)))
               IOFFT = IOFFT + NT1AM(ISYMCJ)*NVIR(ISYMA)
               IOFFX = IOFFX + NT1AM(ISYMCJ)*NVIR(ISYMB)
            END DO
         END DO

      ENDDO
C
      CONTAINS
         PURE SUBROUTINE SQUARE_TOTSYM(T2SQ,T2PACK,ISYM)
C            Should work
C           Square up the amplitudes => For totally symmetric amplitudes
C
C           T2( (ai)<(bj) ) -> T( b,j,a,i)
            INTEGER, INTENT(IN)          :: ISYM
            DOUBLE PRECISION, INTENT(IN) :: T2PACK(*)
            DOUBLE PRECISION, INTENT(OUT) :: T2SQ(*)
C
            INTEGER :: NAI, NBJ
            INTEGER :: IOFFOUT, IOFFPACKED, IDXPACKED,IOFFINP
C
            IOFFOUT = 0
            IOFFINP = IT2AM(ISYM,ISYM)
            DO NAI = 1, NT1AM(ISYM)

               IOFFOUT = (NAI-1)* NT1AM(ISYM)
C              Terms (bj) < (ai)
               IOFFPACKED = IOFFINP + NAI*(NAI-1)/2
               DO NBJ = 1, NAI
                  T2SQ(IOFFOUT+NBJ) = T2PACK(IOFFPACKED+NBJ)
               END DO
C              Terms (ai) < (bj)
               DO NBJ = NAI+1, NT1AM(ISYM)
                  IDXPACKED = IOFFINP + NBJ*(NBJ-1)/2 + NAI
                  T2SQ(IOFFOUT+NBJ) = T2PACK(IDXPACKED)
               ENDDO
C               IOFFOUT = IOFFOUT + NT1AM(ISYM)
            END DO

         END SUBROUTINE

         SUBROUTINE SQUAREX(T2SQ,T2PACK,ISYMR,ISYMT)
C
C           Square up the amplitudes => For totally symmetric amplitudes
C
C           X2( (ai)<(bj) ) -> X( b,j,a,i)
C           For totally symmetric X, scale X( a,i,a,i)
C           With a factor of two
            INTEGER, INTENT(IN)          :: ISYMR,ISYMT
            DOUBLE PRECISION, INTENT(IN) :: T2PACK(*)
            DOUBLE PRECISION, INTENT(OUT) :: T2SQ(*)
C
            INTEGER :: NAI, NBJ, NSIZE
            INTEGER :: IOFFOUT, IOFFPACKED, IDXPACKED, IOFFINP
            INTEGER :: ISYML
C
            ISYML = MULD2H(ISYMR,ISYMT)
            IOFFINP = IT2AM(ISYML,ISYMR)
            IF ( ISYMT .EQ. 1 ) THEN ! Totally symmetric, do the
                                 ! Same as above
               IOFFOUT = 0
               DO NAI = 1, NT1AM(ISYMR)

C                 Terms (bj) < (ai)
                  IOFFPACKED = IOFFINP + NAI*(NAI-1)/2
                  DO NBJ = 1, NAI -1
                     T2SQ(IOFFOUT+NBJ) = T2PACK(IOFFPACKED+NBJ)
                  END DO
C                 Term (ai) = (bj), scale factor of two
                  T2SQ(IOFFOUT+NAI) = 2.0D0*T2PACK(IOFFPACKED+NAI)
C                 Terms (ai) < (bj)
                  DO NBJ = NAI +1, NT1AM(ISYMR)
                     IDXPACKED = IOFFINP + NBJ*(NBJ-1)/2 + NAI
                     T2SQ(IOFFOUT+NBJ) = T2PACK(IDXPACKED)
                  ENDDO
                  IOFFOUT = IOFFOUT + NT1AM(ISYMR)
               END DO
            ELSE IF ( ISYMR .GT. ISYML ) THEN
C              All correctly organized, just need to copy it over
               NSIZE = NT1AM(ISYML)*NT1AM(ISYMR)
               CALL DCOPY(NSIZE,T2PACK(IOFFINP+1),1,T2SQ,1)
            ELSE
C              We need to transpose the array
C              This would probably benefit from from blocking
C              MKL also has the nonstandard mkl_domatcopy to do this
               DO NAI = 1, NT1AM(ISYMR)
                  IOFFOUT = (NAI-1)*NT1AM(ISYML)
                  DO NBJ = 1, NT1AM(ISYML)
                     T2SQ(IOFFOUT+NBJ) =
     &                     T2PACK(IOFFINP+(NBJ-1)*NT1AM(ISYMR)+NAI)
                  END DO
               END DO

            END IF

         END SUBROUTINE

      END
