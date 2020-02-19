C
C  /* Deck so_tmltr */
      SUBROUTINE SO_TMLTR(T2AM,SCAL,ISYOPE)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Andrea Ligabue, January 2004
C
C     Developed starting frm CCSD_TCMEPKX
C     Henrik Koch and Alfredo Sanchez.                Dec 1994
C     Made workable for non-symmetric T2AM, Keld Bak, Dec 1996
C
C     Purpose: calculate the left T_aibj starting from the right ones
C     (I need to call that subroutine with SCAL = HALF)
C
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, THREE = 3.0D0)
      PARAMETER (FOUR = 4.0D0)
C
      DIMENSION T2AM(*)
C
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_TMLTR')
C
      FAC = TWO/THREE
C
      DO 100 ISYMIJ = 1,NSYM
C
         ISYMAB = MULD2H(ISYMIJ,ISYOPE)
C
         DO 110 ISYMJ = 1,NSYM
C
            ISYMI = MULD2H(ISYMJ,ISYMIJ)
C
            IF (ISYMI .GT. ISYMJ) GOTO 110
C
            DO 120 ISYMB = 1,NSYM
C
               ISYMA = MULD2H(ISYMB,ISYMAB)
C
               IF (ISYMA .GT. ISYMB) GOTO 120
C
               ISYMAI = MULD2H(ISYMA,ISYMI)
               ISYMBJ = MULD2H(ISYMB,ISYMJ)
               ISYMBI = MULD2H(ISYMB,ISYMI)
               ISYMAJ = MULD2H(ISYMA,ISYMJ)
C
               DO 130 J = 1,NRHF(ISYMJ)
C
                  IF (ISYMI .EQ. ISYMJ) THEN
                     NRHFI =  J
                  ELSE
                     NRHFI = NRHF(ISYMI)
                  ENDIF
C
               IF ( ISYMAI .EQ. ISYMBJ ) THEN
C
                  DO 140 I = 1,NRHFI
C
                     DO 150 B = 1,NVIR(ISYMB)
C
                        IF (ISYMB .EQ. ISYMA) THEN
                           NVIRA = B
                        ELSE
                           NVIRA = NVIR(ISYMA)
                        ENDIF
C
                        NBI = IT1AM(ISYMB,ISYMI)
     *                      + NVIR(ISYMB)*(I - 1) + B
                        NBJ = IT1AM(ISYMB,ISYMJ)
     *                      + NVIR(ISYMB)*(J - 1) + B
C
                        DO 160 A = 1,NVIRA
C
                           NAI = IT1AM(ISYMA,ISYMI)
     *                         + NVIR(ISYMA)*(I - 1) + A
                           NAJ = IT1AM(ISYMA,ISYMJ)
     *                         + NVIR(ISYMA)*(J - 1) + A
C
                           NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                           + INDEX(NAI,NBJ)
C
                           NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *                           + INDEX(NAJ,NBI)
C
                           XAIBJ = FAC*(TWO*T2AM(NAIBJ)+T2AM(NAJBI))
                           XAJBI = FAC*(TWO*T2AM(NAJBI)+T2AM(NAIBJ))
C
                           T2AM(NAIBJ) = XAIBJ
                           T2AM(NAJBI) = XAJBI
C
  160                   CONTINUE
  150                CONTINUE
  140             CONTINUE
C
               ELSE IF ((ISYMAI.LT.ISYMBJ).AND.(ISYMAJ.LT.ISYMBI)) THEN
C
                  DO 240 I = 1,NRHFI
C
                     DO 250 B = 1,NVIR(ISYMB)
C
                        IF (ISYMB .EQ. ISYMA) THEN
                           NVIRA = B
                        ELSE
                           NVIRA = NVIR(ISYMA)
                        ENDIF
C
                        NBI = IT1AM(ISYMB,ISYMI)
     *                      + NVIR(ISYMB)*(I - 1) + B
                        NBJ = IT1AM(ISYMB,ISYMJ)
     *                      + NVIR(ISYMB)*(J - 1) + B
C
                        DO 260 A = 1,NVIRA
C
                           NAI = IT1AM(ISYMA,ISYMI)
     *                         + NVIR(ISYMA)*(I - 1) + A
                           NAJ = IT1AM(ISYMA,ISYMJ)
     *                         + NVIR(ISYMA)*(J - 1) + A
C
                           NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                           + NT1AM(ISYMAI) * (NBJ - 1) + NAI
C
                           NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *                           + NT1AM(ISYMAJ) * (NBI - 1) + NAJ
C
                           XAIBJ = FAC*(TWO*T2AM(NAIBJ)+T2AM(NAJBI))
                           XAJBI = FAC*(TWO*T2AM(NAJBI)+T2AM(NAIBJ))
C
                           T2AM(NAIBJ) = XAIBJ
                           T2AM(NAJBI) = XAJBI
C
  260                   CONTINUE
  250                CONTINUE
  240             CONTINUE
C
               ELSE IF ((ISYMAI.GT.ISYMBJ).AND.(ISYMAJ.GT.ISYMBI)) THEN
C
                  DO 340 I = 1,NRHFI
C
                     DO 350 B = 1,NVIR(ISYMB)
C
                        IF (ISYMB .EQ. ISYMA) THEN
                           NVIRA = B
                        ELSE
                           NVIRA = NVIR(ISYMA)
                        ENDIF
C
                        NBI = IT1AM(ISYMB,ISYMI)
     *                      + NVIR(ISYMB)*(I - 1) + B
                        NBJ = IT1AM(ISYMB,ISYMJ)
     *                      + NVIR(ISYMB)*(J - 1) + B
C
                        DO 360 A = 1,NVIRA
C
                           NAI = IT1AM(ISYMA,ISYMI)
     *                         + NVIR(ISYMA)*(I - 1) + A
                           NAJ = IT1AM(ISYMA,ISYMJ)
     *                         + NVIR(ISYMA)*(J - 1) + A
C
                           NAIBJ = IT2AM(ISYMBJ,ISYMAI)
     *                           + NT1AM(ISYMBJ) * (NAI - 1) + NBJ
C
                           NAJBI = IT2AM(ISYMBI,ISYMAJ)
     *                           + NT1AM(ISYMBI) * (NAJ - 1) + NBI
C
                           XAIBJ = FAC*(TWO*T2AM(NAIBJ) + T2AM(NAJBI))
                           XAJBI = FAC*(TWO*T2AM(NAJBI) + T2AM(NAIBJ))
C
                           T2AM(NAIBJ) = XAIBJ
                           T2AM(NAJBI) = XAJBI
C
  360                   CONTINUE
  350                CONTINUE
  340             CONTINUE
C
               ELSE IF ((ISYMAI.LT.ISYMBJ).AND.(ISYMAJ.GT.ISYMBI)) THEN
C
                  DO 440 I = 1,NRHFI
C
                     DO 450 B = 1,NVIR(ISYMB)
C
                        IF (ISYMB .EQ. ISYMA) THEN
                           NVIRA = B
                        ELSE
                           NVIRA = NVIR(ISYMA)
                        ENDIF
C
                        NBI = IT1AM(ISYMB,ISYMI)
     *                      + NVIR(ISYMB)*(I - 1) + B
                        NBJ = IT1AM(ISYMB,ISYMJ)
     *                      + NVIR(ISYMB)*(J - 1) + B
C
                        DO 460 A = 1,NVIRA
C
                           NAI = IT1AM(ISYMA,ISYMI)
     *                         + NVIR(ISYMA)*(I - 1) + A
                           NAJ = IT1AM(ISYMA,ISYMJ)
     *                         + NVIR(ISYMA)*(J - 1) + A
C
                           NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                           + NT1AM(ISYMAI) * (NBJ - 1) + NAI
C
                           NAJBI = IT2AM(ISYMBI,ISYMAJ)
     *                           + NT1AM(ISYMBI) * (NAJ - 1) + NBI
C
                           XAIBJ = FAC*(TWO*T2AM(NAIBJ) + T2AM(NAJBI))
                           XAJBI = FAC*(TWO*T2AM(NAJBI) + T2AM(NAIBJ))
C
                           T2AM(NAIBJ) = XAIBJ
                           T2AM(NAJBI) = XAJBI
C
  460                   CONTINUE
  450                CONTINUE
  440             CONTINUE
C
               ELSE IF ((ISYMAI.GT.ISYMBJ).AND.(ISYMAJ.LT.ISYMBI)) THEN
C
                  DO 540 I = 1,NRHFI
C
                     DO 550 B = 1,NVIR(ISYMB)
C
                        IF (ISYMB .EQ. ISYMA) THEN
                           NVIRA = B
                        ELSE
                           NVIRA = NVIR(ISYMA)
                        ENDIF
C
                        NBI = IT1AM(ISYMB,ISYMI)
     *                      + NVIR(ISYMB)*(I - 1) + B
                        NBJ = IT1AM(ISYMB,ISYMJ)
     *                      + NVIR(ISYMB)*(J - 1) + B
C
                        DO 560 A = 1,NVIRA
C
                           NAI = IT1AM(ISYMA,ISYMI)
     *                         + NVIR(ISYMA)*(I - 1) + A
                           NAJ = IT1AM(ISYMA,ISYMJ)
     *                         + NVIR(ISYMA)*(J - 1) + A
C
                           NAIBJ = IT2AM(ISYMBJ,ISYMAI)
     *                           + NT1AM(ISYMBJ) * (NAI - 1) + NBJ
C
                           NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *                           + NT1AM(ISYMAJ) * (NBI - 1) + NAJ
C
                           XAIBJ = FAC*(TWO*T2AM(NAIBJ) + T2AM(NAJBI))
                           XAJBI = FAC*(TWO*T2AM(NAJBI) + T2AM(NAIBJ))
C
                           T2AM(NAIBJ) = XAIBJ
                           T2AM(NAJBI) = XAJBI
C
  560                   CONTINUE
  550                CONTINUE
  540             CONTINUE
C
               END IF
C
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
C---------------------------------------
C     Scale diagonal elements of result.
C---------------------------------------
C
      IF ((ISYOPE .NE. 1).OR.(SCAL.EQ.1.0D0)) GOTO 1000
C
      DO 600 ISYMAI = 1,NSYM
         DO 610 NAI = 1,NT1AM(ISYMAI)
            NAIAI = IT2AM(ISYMAI,ISYMAI) + INDEX(NAI,NAI)
            T2AM(NAIAI) = SCAL*T2AM(NAIAI)
  610    CONTINUE
  600 CONTINUE
C
C-----------------------
C     Remove from trace.
C-----------------------
C
 1000 CALL QEXIT('SO_TMLTR')
C
      RETURN
      END
C

      SUBROUTINE SO_TRANTRIP(X2SQ,X2PACK,ISYMTR)

      use so_info, only: sop_dp
      implicit none

#include "ccorb.h"
#include "ccsdsym.h"
#include "soppinf.h"

      real(sop_dp),intent(out) :: X2SQ(*)
      real(sop_dp),intent(in) :: x2pack(*)
      integer, intent(in) :: isymtr

      character(len=*), parameter :: myname = 'SO_TRANTRIP'
      integer :: isymbj, isymai
      integer :: KOUT

      CALL QENTER(myname)

      KOUT = 1
      do isymbj = 1, nsym
         isymai = muld2h(isymbj,isymtr)
C
C        Calculate the intermediate
C           ~         _
C        x(aibj) = - /2 x1(abij) - x2(abij) - x3(abij)
C
C        Note that x2 changes sign when permuting i and j,
C        x3 when permuting a and b, and x1 change sign on both
C        these permutations.
C
         CALL SQUAREXT(X2SQ(KOUT),X2PACK,ISYMBJ,ISYMTR)
         KOUT = KOUT + NT1AM(ISYMAI)*NT1AM(ISYMBJ)

      end do

      IF ( KOUT-1.NE. NT2SQ(ISYMTR) ) THEN
         PRINT *, 'WARNING', KOUT, NT2SQ(ISYMTR)
      END IF

      CALL QEXIT(myname)

      contains

         SUBROUTINE SQUAREXT(X2SQ,X2PACK,ISYMR,ISYMT)
            real(sop_dp),intent(out) :: X2SQ(*)
            real(sop_dp),intent(in) :: x2pack(*)
            integer, intent(in) :: isymr, isymt

            real(sop_dp),parameter :: factd = -SQRT(2.0D0),
     &                                fact1 = -SQRT(2.0D0),
     &                                fact23 = -1.D0,
     &                                zero = 0.0d0,
     &                                one = 1.0D0,
     &                                onem = -1.0D0

            integer :: isymbj, isymai, isymaj, isymbi, isymij, isymab,
     &                 isyma, isymb, isymi, isymj
            integer :: sab1, sab2, nvira, nrhfi, nvirb, nrhfj
            integer :: nai, naj, nbi, nbj, naibj, najbi, nbiaj, nbjai,
     &                 nbibj, nbjbi, najbj, nbjaj, nbjbj,
     &                 nabij1, nabij2, nabij3, nbbij2, nabjj3
            integer :: i, j, b, a
            integer :: ioff1, ioffj1, ioffbj1
            integer :: ioff2, ioffj2, ioffbj2
            integer :: ioff3, ioffj3, ioffbj3
            real(sop_dp) :: f, fij, fab

            isymbj = isymr
            isymai = muld2h(isymr,isymt)
C
C           Calculate the intermediate
C           ~            _
C           x(aibj) = - /2 x1(abij) - x2(abij) - x3(abij)
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

      END SUBROUTINE


