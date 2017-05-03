      SUBROUTINE SO_PERTD2(ISYMTR,FACTOR,
     &                     RDENSIJ,LRDENSIJ,RDENSAB,LRDENSAB,
     &                     T2AM,LT2AM,SOLVEC2,LSOLVEC2,
     &                     WORK,LWORK)
C
      use so_info, only : sop_dp
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
         IF (TRIPLET) THEN
            CALL SQUAREXT(WORK(KX2SQ),SOLVEC2,ISYMAI,ISYMTR)
         ELSE
C        X2( (ai)<(bj) ) -> X2( bj , ai), Diagonal terms scaled by two
            CALL SQUAREX(WORK(KX2SQ),SOLVEC2,ISYMAI,ISYMTR)
         ENDIF
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
C           Square up the amplitudes => For general amplitudes
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
C
            END IF
C
         END SUBROUTINE

         SUBROUTINE SQUAREXT(X2SQ,X2PACK,ISYMR,ISYMT)
            real(sop_dp),intent(out) :: X2SQ(*)
            real(sop_dp),intent(in) :: x2pack(*)
            integer, intent(in) :: isymr, isymt

            real(sop_dp),parameter :: factd =  SQRT(2.0D0),
     &                                fact1 =  -SQRT(2.0D0),
     &                                fact23 =  1.D0,
     &                                zero = 0.0d0

            integer :: isymbj, isymai, isymaj, isymbi, isymij, isymab,
     &                 isyma, isymb, isymi, isymj
            integer :: sab1, sab2, nvira, nrhfi, nvirb, nrhfj
            integer :: nai, naj, nbi, nbj, naibj, najbi, nbiaj, nbjai,
     &                 nbibj, nbjbi, najbj, nbjaj, nbjbj,
     &                 nabij1, nabij2, nabij3, nbbij2, nabjj3
            integer :: ioff1, ioffj1, ioffbj1
            integer :: ioff2, ioffj2, ioffbj2
            integer :: ioff3, ioffj3, ioffbj3
            real(sop_dp) :: f, fij, fab

            isymbj = isymr
            isymai = muld2h(isymr,isymt)
C
C           Calculate the intermediate
C           ~           _
C           x(aibj) = -/2 x1(abij) + x2(abij) + x3(abij)
C
C           Note that x2 changes sign when permuting i and j,
C           x3 when permuting a and b, and x1 change sign on both
C           these permutations.
C

            if (isymt .eq. 1 ) then
               do isymj = 1, nsym
                  isymb = muld2h(isymj,isymbj)

               ! Ensure isymi <= isymj
                  do isymi = 1, isymj
                     isyma = muld2h(isymi,isymai)
                     ! and isuma <= isymb
C                     if (isyma.gt.isymb) cycle
C
                     isymaj = muld2h(isyma,isymj)
                     isymbi = muld2h(isymb,isymi)
                     isymab = muld2h(isyma,isymb)
                     isymij = muld2h(isymi,isymj)
C
                     sab1 = ntvv(isymab)
                     sab2 = nsvv(isymab)
                     ioff1 = it2amt1(isymij,isymab) +
     &                     sab1*itoo(isymj,isymi) + itvv(isymb,isyma)
                     ioff2 = it2amt2(isymij,isymab) +
     &                     sab2*itoo(isymj,isymi) + isvv(isymb,isyma)
                     ioff3 = it2amt3(isymij,isymab) +
     &                     sab1*isoo(isymj,isymi) + itvv(isymb,isyma)
C
                  if (isymi.eq.isymj) then ! isyma == isymb
                     do j = 1, nrhf(isymj)
                        nrhfi = j - 1
                        ioffj1 = ioff1 + ((j-1)*(j-2)/2)*sab1
                        ioffj2 = ioff2 + ((j-1)*(j-2)/2)*sab2
                        ioffj3 = ioff3 + (j*(j-1)/2)*sab1
                        do b = 1, nvir(isymb)
                           nbj = pair_pos(isymb,b,isymj,j)
                           nvira = b - 1
                           ioffbj1 = ioffj1 + (b-1)*(b-2)/2
                           ioffbj2 = ioffj2 + b*(b-1)/2
                           ioffbj3 = ioffj3 + (b-1)*(b-2)/2

                           do i = 1, nrhfi
                              nbi = pair_pos(isymb,b,isymi,i)
                              do a = 1, nvira
                                 nai = pair_pos(isyma,a,isymi,i)
                                 naj = pair_pos(isyma,a,isymj,j)
                                 nabij1 = ioffbj1+(i-1)*sab1+a
                                 nabij2 = ioffbj2+(i-1)*sab2+a
                                 nabij3 = ioffbj3+(i-1)*sab1+a
                                 naibj = quad_pos(isymai,nai,isymbj,nbj)
                                 najbi = quad_pos(isymaj,naj,isymbi,nbi)
                                 nbiaj = quad_pos(isymbi,nbi,isymaj,naj)
                                 nbjai = quad_pos(isymbj,nbj,isymai,nai)
                                 x2sq(naibj) =  fact1*x2pack(nabij1)
     &                      + fact23*( x2pack(nabij2) + x2pack(nabij3) )
                                 x2sq(najbi) = -fact1*x2pack(nabij1)
     &                      - fact23*( x2pack(nabij2) - x2pack(nabij3) )
                                 x2sq(nbiaj) = -fact1*x2pack(nabij1)
     &                      + fact23*( x2pack(nabij2) - x2pack(nabij3) )
                                 x2sq(nbjai) =  fact1*x2pack(nabij1)
     &                      - fact23*( x2pack(nabij2) + x2pack(nabij3) )
                              end do ! a
C                          Handle a == b, i<j
                              nbbij2 = ioffbj2+(i-1)*sab2+b
                              nbjbi = quad_pos(isymbj,nbj,isymbi,nbi)
                              nbibj = quad_pos(isymbi,nbi,isymbj,nbj)
                              x2sq(nbibj) = factd*x2pack(nbbij2)
                              x2sq(nbjbi) = - factd*x2pack(nbbij2)
                           end do ! i
C                       Handle i = j, a < b
                           do a = 1, nvira
                              naj = pair_pos(isyma,a,isymj,j)
                              nabjj3 = ioffbj3+(j-1)*sab1+a
                              najbj = quad_pos(isymaj,naj,isymbj,nbj)
                              nbjaj = quad_pos(isymbj,nbj,isymaj,naj)
                              x2sq(najbj) = factd*x2pack(nabjj3)
                              x2sq(nbjaj) = -factd*x2pack(nabjj3)
                           end do
                           ! i == j, a == b
C                           if (isyma.eq.isymb) then
                           nbjbj = quad_pos(isymbj,nbj,isymbj,nbj)
                           x2sq(nbjbj) = zero
C                           end if

                        end do ! b
                     end do ! j
                  else if (isyma .lt. isymb ) then
                     do j = 1, nrhf(isymj)
                        nrhfi = nrhf(isymi)
                        ioffj1 = ioff1 + (j-1)*nrhfi*sab1
                        ioffj2 = ioff2 + (j-1)*nrhfi*sab2
                        ioffj3 = ioff3 + (j-1)*nrhfi*sab1
                        do b = 1, nvir(isymb)
                           nbj = pair_pos(isymb,b,isymj,j)
                           nvira = nvir(isyma)
                           ioffbj1 = ioffj1 + (b-1)*nvira
                           ioffbj2 = ioffj2 + (b-1)*nvira
                           ioffbj3 = ioffj3 + (b-1)*nvira
                           do i = 1, nrhfi
                              nbi = pair_pos(isymb,b,isymi,i)
                              do a = 1, nvira
                                 nai = pair_pos(isyma,a,isymi,i)
                                 naj = pair_pos(isyma,a,isymj,j)
                                 nabij1 = ioffbj1+(i-1)*sab1+a
                                 nabij2 = ioffbj2+(i-1)*sab2+a
                                 nabij3 = ioffbj3+(i-1)*sab1+a
                                 naibj = quad_pos(isymai,nai,isymbj,nbj)
                                 nbjai = quad_pos(isymbj,nbj,isymai,nai)
                                 x2sq(naibj) =  fact1*x2pack(nabij1)
     &                      + fact23*( x2pack(nabij2) + x2pack(nabij3) )
                                 x2sq(nbjai) =  fact1*x2pack(nabij1)
     &                      - fact23*( x2pack(nabij2) + x2pack(nabij3) )
                              end do ! a
                          end do ! i
                       end do ! b
                    end do ! j
                  else ! isyma > isymb -> a > b !
                     ioff1 = it2amt1(isymij,isymab) +
     &                     sab1*itoo(isymj,isymi) + itvv(isyma,isymb)
                     ioff2 = it2amt2(isymij,isymab) +
     &                     sab2*itoo(isymj,isymi) + isvv(isyma,isymb)
                     ioff3 = it2amt3(isymij,isymab) +
     &                     sab1*isoo(isymj,isymi) + itvv(isyma,isymb)
                     nvirb = nvir(isymb)
                     do j = 1, nrhf(isymj)
                        nrhfi = nrhf(isymi)
                        ioffj1 = ioff1 + (j-1)*nrhfi*sab1
                        ioffj2 = ioff2 + (j-1)*nrhfi*sab2
                        ioffj3 = ioff3 + (j-1)*nrhfi*sab1
                        do b = 1, nvir(isymb)
                           nbj = pair_pos(isymb,b,isymj,j)
                           nvira = nvir(isyma)
                           ioffbj1 = ioffj1 + b
                           ioffbj2 = ioffj2 + b
                           ioffbj3 = ioffj3 + b
                           do i = 1, nrhfi
                              nbi = pair_pos(isymb,b,isymi,i)
                              do a = 1, nvira
                                 nai = pair_pos(isyma,a,isymi,i)
                                 naj = pair_pos(isyma,a,isymj,j)
                                 nabij1 = ioffbj1+(i-1)*sab1+(a-1)*nvirb
                                 nabij2 = ioffbj2+(i-1)*sab2+(a-1)*nvirb
                                 nabij3 = ioffbj3+(i-1)*sab1+(a-1)*nvirb
                                 naibj = quad_pos(isymai,nai,isymbj,nbj)
                                 nbjai = quad_pos(isymbj,nbj,isymai,nai)
                                 x2sq(naibj) = -fact1*x2pack(nabij1)
     &                      + fact23*( x2pack(nabij2) - x2pack(nabij3) )
                                 x2sq(nbjai) = -fact1*x2pack(nabij1)
     &                      - fact23*( x2pack(nabij2) - x2pack(nabij3) )
                              end do ! a1
                          end do ! i
                       end do ! b
                    end do ! j
                  end if

                  end do ! loop isymi
               end do ! loop isymj
C
            else ! isymai != isymbj
C
               do isymj = 1, nsym
                  isymb = muld2h(isymj,isymbj)

                  do isymi = 1, nsym
                     isyma = muld2h(isymi,isymai)
                     ! and isuma <= isymb
C                     if (isyma.gt.isymb) cycle
C
                     isymaj = muld2h(isyma,isymj)
                     isymbi = muld2h(isymb,isymi)
                     isymab = muld2h(isyma,isymb)
                     isymij = muld2h(isymi,isymj)
C
                     sab1 = ntvv(isymab)
                     sab2 = nsvv(isymab)
                     ioff1 = it2amt1(isymij,isymab) +
     &                     sab1*itoo(isymj,isymi) + itvv(isymb,isyma)
                     ioff2 = it2amt2(isymij,isymab) +
     &                     sab2*itoo(isymj,isymi) + isvv(isymb,isyma)
                     ioff3 = it2amt3(isymij,isymab) +
     &                     sab1*isoo(isymj,isymi) + itvv(isymb,isyma)

                  if ((isymi .eq. isymj)) then ! isyma != isymb
                     if (isyma .lt. isymb ) then
                        nvira = nvir(isyma)
                        nvirb = 1
                        f = one
                     else
                        nvira = 1
                        nvirb = nvir(isymb)
                        f = onem
                     end if
                     do j = 1, nrhf(isymj)
                        nrhfi = j - 1
                        ioffj1 = ioff1 + ((j-1)*(j-2)/2)*sab1
                        ioffj2 = ioff2 + ((j-1)*(j-2)/2)*sab2
                        ioffj3 = ioff3 + (j*(j-1)/2)*sab1
                        do b = 1, nvir(isymb)
                           nbj = pair_pos(isymb,b,isymj,j)
                           ioffbj1 = ioffj1 + (b-1)*nvira
                           ioffbj2 = ioffj2 + (b-1)*nvira
                           ioffbj3 = ioffj3 + (b-1)*nvira
                           do i = 1, nrhfi
                              nbi = pair_pos(isymb,b,isymi,i)
                              do a = 1, nvir(isyma)
                                 nai = pair_pos(isyma,a,isymi,i)
                                 naj = pair_pos(isyma,a,isymj,j)
                                 nabij1 = ioffbj1 + (i-1)*sab1+
     &                                    (a-1)*nvirb + 1
                                 nabij2 = ioffbj2 + (i-1)*sab2+
     &                                    (a-1)*nvirb + 1
                                 nabij3 = ioffbj3 + (i-1)*sab1+
     &                                    (a-1)*nvirb + 1
                                 naibj = quad_pos(isymai,nai,isymbj,nbj)
                                 najbi = quad_pos(isymaj,naj,isymbi,nbi)
                                 x2sq(naibj) =  f*fact1*x2pack(nabij1)
     &                    + fact23*( x2pack(nabij2) + f*x2pack(nabij3) )
                                 x2sq(najbi) = -f*fact1*x2pack(nabij1)
     &                    - fact23*( x2pack(nabij2) - f*x2pack(nabij3) )
                              end do ! a
                           end do ! i
C                       Handle i = j, a != b
                           do a = 1, nvir(isyma)
                              naj = pair_pos(isyma,a,isymj,j)
                              nabjj3 = ioffbj3+(j-1)*sab1+(a-1)*nvirb+1
                              najbj = quad_pos(isymaj,naj,isymbj,nbj)
                              x2sq(najbj) =  f*factd*x2pack(nabjj3)
                           end do
                        end do ! b
                     end do ! j
                  else if (isyma .eq. isymb ) then ! isymi != isymj
                     if (isymi .lt. isymj ) then
                        nrhfi = nrhf(isymi)
                        nrhfj = 1
                        fij = one
                     else
                        nrhfi = 1
                        nrhfj = nrhf(isymj)
                        fij = onem
                     end if
                     do j = 1, nrhf(isymj)
                        ioffj1 = ioff1 + (j-1)*nrhfi*sab1
                        ioffj2 = ioff2 + (j-1)*nrhfi*sab2
                        ioffj3 = ioff3 + (j-1)*nrhfi*sab1
                        do b = 1, nvir(isymb)
                           nbj = pair_pos(isymb,b,isymj,j)
                           nvira = b - 1
                           ioffbj1 = ioffj1 + (b-1)*(b-2)/2
                           ioffbj2 = ioffj2 + b*(b-1)/2
                           ioffbj3 = ioffj3 + (b-1)*(b-2)/2

                           do i = 1, nrhf(isymi)
                              nbi = pair_pos(isymb,b,isymi,i)
                              do a = 1, nvira
                                 nai = pair_pos(isyma,a,isymi,i)
                                 naj = pair_pos(isyma,a,isymj,j)
                                 nabij1 = ioffbj1+(i-1)*nrhfj*sab1+a
                                 nabij2 = ioffbj2+(i-1)*nrhfj*sab2+a
                                 nabij3 = ioffbj3+(i-1)*nrhfj*sab1+a
                                 naibj = quad_pos(isymai,nai,isymbj,nbj)
                                 nbiaj = quad_pos(isymbi,nbi,isymaj,naj)
                                 x2sq(naibj) =
     &                               fij*fact1*x2pack(nabij1)
     &                              + fact23*( fij*x2pack(nabij2) +
     &                                         x2pack(nabij3) )
                                 x2sq(nbiaj) =
     &                              -fij*fact1*x2pack(nabij1)
     &                              + fact23*( fij*x2pack(nabij2) -
     &                                         x2pack(nabij3) )
                              end do ! a
C                          Handle a == b, i<j
                              nbbij2 = ioffbj2+(i-1)*nrhfj*sab2+b
                              nbibj = quad_pos(isymbi,nbi,isymbj,nbj)
                              x2sq(nbibj) = fij*factd*x2pack(nbbij2)
                           end do ! i
                        end do ! b
                     end do ! j

                  else ! a != b, i != j
                     if (isyma .lt. isymb ) then
                        nvira = nvir(isyma)
                        nvirb = 1
                        fab = one
                     else
                        nvira = 1
                        nvirb = nvir(isymb)
                        fab = onem
                     end if
                     if (isymi .lt. isymj ) then
                        nrhfi = nrhf(isymi)
                        nrhfj = 1
                        fij = one
                     else
                        nrhfi = 1
                        nrhfj = nrhf(isymj)
                        fij = onem
                     end if
                     do j = 1, nrhf(isymj)
                        ioffj1 = ioff1 + (j-1)*nrhfi*sab1
                        ioffj2 = ioff2 + (j-1)*nrhfi*sab2
                        ioffj3 = ioff3 + (j-1)*nrhfi*sab1
                        do b = 1, nvir(isymb)
                           nbj = pair_pos(isymb,b,isymj,j)
                           ioffbj1 = ioffj1 + (b-1)*nvira
                           ioffbj2 = ioffj2 + (b-1)*nvira
                           ioffbj3 = ioffj3 + (b-1)*nvira
                           do i = 1, nrhf(isymi)
                              do a = 1, nvir(isyma)
                                 nai = pair_pos(isyma,a,isymi,i)
                                 nabij1 = ioffbj1 +
     &                                    (i-1)*nrhfj*sab1+
     &                                    (a-1)*nvirb + 1
                                 nabij2 = ioffbj2 +
     &                                    (i-1)*nrhfj*sab2+
     &                                    (a-1)*nvirb + 1
                                 nabij3 = ioffbj3 +
     &                                    (i-1)*nrhfj*sab1+
     &                                    (a-1)*nvirb + 1
                                 naibj = quad_pos(isymai,nai,isymbj,nbj)
                                 x2sq(naibj) =
     &                                fab*fij*fact1*x2pack(nabij1)
     &                                + fact23*( fij*x2pack(nabij2) +
     &                                           fab*x2pack(nabij3) )
                              end do ! a
                           end do ! i
                        end do ! b
                     end do ! j

                  end if
                  end do ! isymi
               end do ! isymj
            end if
         END SUBROUTINE

         pure function pair_pos(isyma,na,isymi,ni)
            integer :: pair_pos
            integer, intent(in) :: na, ni, isyma, isymi

            pair_pos = it1am(isyma,isymi) + nvir(isyma)*(ni-1) + na
            return
         end function

         pure function quad_pos(isymai,nai,isymbj,nbj)
            integer :: quad_pos
            integer, intent(in) :: isymai, isymbj, nai, nbj
            quad_pos = nt1am(isymai)*(nbj-1) + nai
            return
         end function

C
      END
