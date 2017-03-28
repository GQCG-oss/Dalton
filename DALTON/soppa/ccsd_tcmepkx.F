C  /* Deck ccsd_tcmepkx */
      SUBROUTINE CCSD_TCMEPKX(T2AM,SCAL,ISYOPE)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Henrik Koch and Alfredo Sanchez.                Dec 1994
C     Made workable for non-symmetric T2AM, Keld Bak, Dec 1996
C
C     Calculate in place two coulomb minus exchange of t2 amplitudes.
C
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
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
      CALL QENTER('CCSD_TCMEPKX')
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
                           XAIBJ = TWO*T2AM(NAIBJ) - T2AM(NAJBI)
                           XAJBI = TWO*T2AM(NAJBI) - T2AM(NAIBJ)
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
                           XAIBJ = TWO*T2AM(NAIBJ) - T2AM(NAJBI)
                           XAJBI = TWO*T2AM(NAJBI) - T2AM(NAIBJ)
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
                           XAIBJ = TWO*T2AM(NAIBJ) - T2AM(NAJBI)
                           XAJBI = TWO*T2AM(NAJBI) - T2AM(NAIBJ)
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
                           XAIBJ = TWO*T2AM(NAIBJ) - T2AM(NAJBI)
                           XAJBI = TWO*T2AM(NAJBI) - T2AM(NAIBJ)
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
                           XAIBJ = TWO*T2AM(NAIBJ) - T2AM(NAJBI)
                           XAJBI = TWO*T2AM(NAJBI) - T2AM(NAIBJ)
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
      IF (ISYOPE .NE. 1) GOTO 1000
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
 1000 CALL QEXIT('CCSD_TCMEPKX')
C
      RETURN
      END

      FUNCTION SO_DOUBLES_LENGTH(X2AM,ISYMX2)
C
C     Returns the square norm of a singlet adapted doubles vector
C
C     This is given as
C                           ~
C     0.5 * sum_{(ai)<(bj)} x(aibj)*x(aibj) + sum_{ai} x(aiai)**2
C          ~
C     with x(aibj) = 2*x(aibj) - x(ajbi)
C     Inserting this in the above expression followed by some
C     manipulations of the first sum leads to
C
C     sum_{a<b,i<j} (x(aibj)**2 + x(ajbi)**2 - x(aibj) * x(ajbi) )
C     + 0.5*sum_{a<b,i} x(aibi)**2 + 0.5*sum{a,i<j} x(aiaj)**2
C     + sum_{ai} x(aiai)**2
C
C     which is implemented in the present routine.

      use so_info, only: sop_dp
      implicit none
C
      real(sop_dp),parameter :: ZERO = 0.0D0, HALF = 0.5D0, two = 2.0D0
C
      real(sop_dp), intent(in) :: X2AM(*)
      integer, intent(in)      :: isymx2
      real(sop_dp)             :: so_doubles_length
C
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
C
      INTEGER :: ISYMAB, ISYMIJ, ISYMAI, ISYMAJ, ISYMBJ, ISYMBI,
     &           ISYMA, ISYMB, ISYMI, ISYMJ
      INTEGER :: NAI, NBJ, NAJ, NBI, NAIAI, NAIBJ, NAJBI,
     &           NAJBJ, NBIBJ, NBJBJ
      INTEGER :: NRHFI, NVIRA, INDEX
      real(sop_dp) :: XAIBJ, XAJBI
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      SO_DOUBLES_LENGTH = ZERO

      IF (ISYMX2 .eq. 1) THEN

         do isymbj = 1, nsym
            isymai = isymbj

            do isymj = 1, nsym
               isymb = muld2h(isymj,isymbj)

               ! Ensure isymi <= isymj
               do isymi = 1, isymj
                  isyma = muld2h(isymi,isymai)
                  if (isyma.gt.isymb) cycle

                  isymaj = muld2h(isyma,isymj)
                  isymbi = muld2h(isymb,isymi)

                  do J = 1, nrhf(isymj)
                     ! Ensure i<j
                     if (isymi .eq. isymj) then
                        nrhfi = j-1
                     else
                        nrhfi = nrhf(isymi)
                     endif

                     do b = 1, nvir(isymb)
                        ! Ensure a<b
                        if (isyma.eq.isymb) then
                           nvira = b-1
                        else
                           nvira = nvir(isyma)
                        endif
                        nbj = pair_pos(isymb,b,isymj,j)
C
                        ! Do a<b
                        do i = 1, nrhfi
                           nbi = pair_pos(isymb,b,isymi,i)
                           do a = 1, nvira
                              nai = pair_pos(isyma,a,isymi,i)
                              naj = pair_pos(isyma,a,isymj,j)
                              naibj = it2am(isymai,isymbj) +
     &                                nbj*(nbj-1)/2 + nai
                              najbi = it2am(isymbi,isymaj) +
     &                                naj*(naj-1)/2 + nbi
                              xaibj = x2am(naibj)
                              xajbi = x2am(najbi)
                              SO_DOUBLES_LENGTH = SO_DOUBLES_LENGTH
     &                                          + XAIBJ*XAIBJ
     &                                          + XAJBI*XAJBI
     &                                          - XAJBI*XAIBJ
                           end do ! loop a
C
                           ! handle explicitly a == b, i<j
                           if (isyma.eq.isymb) then
                              nbibj = it2am(isymbi,isymbj) +
     &                                nbj*(nbj-1)/2 + nbi
                              so_doubles_length = so_doubles_length
     &                                          + half*(x2am(nbibj)**2)
                           end if
                        end do ! loop i

                        ! handle explicitly a <= b, i == j
                        if (isymi.eq.isymj) then
                           ! a<b
                           do a = 1, nvira
                              naj = pair_pos(isyma,a,isymj,j)
                              najbj = it2am(isymaj,isymbj) +
     &                                nbj*(nbj-1)/2 + naj
                              so_doubles_length = so_doubles_length
     &                                          + half*(x2am(najbj)**2)
                           end do
                           ! a==b, i==j
                           nbjbj = it2am(isymbj,isymbj) +
     &                             nbj*(nbj+1)/2
                           so_doubles_length = so_doubles_length
     &                                       + x2am(nbjbj)**2
                        end if
                     end do ! loop b
                  end do ! loop j
               end do ! loop isymi
            end do ! loop isymj
         end do ! loop isymbj
      else ! not totally symmetric case
C
         do isymbj = 1, nsym
            isymai = muld2h(isymbj,isymx2)

            do isymj = 1, nsym
               isymb = muld2h(isymj,isymbj)

               do isymi = 1, isymj
                  isyma = muld2h(isymi,isymai)
                  if (isyma.gt.isymb) cycle

                  isymaj = muld2h(isyma,isymj)
                  isymbi = muld2h(isymb,isymi)

                  do J = 1, nrhf(isymj)
                     ! Ensure i<j
                     if (isymi .eq. isymj) then
                        nrhfi = j-1
                     else
                        nrhfi = nrhf(isymi)
                     endif

                     do b = 1, nvir(isymb)
                        ! Ensure a<b
                        if (isyma.eq.isymb) then
                           nvira = b-1
                        else
                           nvira = nvir(isyma)
                        endif
                        nbj = pair_pos(isymb,b,isymj,j)
C
                        ! Do a<b
                        do i = 1, nrhfi
                           nbi = pair_pos(isymb,b,isymi,i)
                           do a = 1, nvira
                              nai = pair_pos(isyma,a,isymi,i)
                              naj = pair_pos(isyma,a,isymj,j)
                              naibj = quad_pos(isymai,nai,isymbj,nbj)
                              najbi = quad_pos(isymaj,naj,isymbi,nbi)
                              xaibj = x2am(naibj)
                              xajbi = x2am(najbi)
                              SO_DOUBLES_LENGTH = SO_DOUBLES_LENGTH
     &                                          + XAIBJ*XAIBJ
     &                                          + XAJBI*XAJBI
     &                                          - XAJBI*XAIBJ
                           end do ! loop a
C
                           ! handle explicitly a == b, i<j
                           if (isyma.eq.isymb) then
                              nbibj = quad_pos(isymbi,nbi,isymbj,nbj)
                              so_doubles_length = so_doubles_length
     &                                          + half*(x2am(nbibj)**2)
                           end if
                        end do ! loop i

                        ! handle explicitly a < b, i == j
                        if (isymi.eq.isymj) then
                           ! a<b
                           do a = 1, nvira
                              naj = pair_pos(isyma,a,isymj,j)
                              najbj = quad_pos(isymaj,naj,isymbj,nbj)
                              so_doubles_length = so_doubles_length
     &                                          + half*(x2am(najbj)**2)
                           end do
                        end if
                     end do ! loop b
                  end do ! loop j
               end do ! loop isymi
            end do ! loop isymj
         end do ! loop isymbj
C
      end if

!      SO_DOUBLES_LENGTH = 0.5D0*SO_DOUBLES_LENGTH

      contains
         pure function pair_pos(isyma,na,isymi,ni)
            integer :: pair_pos
            integer, intent(in) :: na, ni, isyma, isymi

            pair_pos = it1am(isyma,isymi) + nvir(isyma)*(ni-1) + na
            return
         end function

         pure function quad_pos(isymai,nai,isymbj,nbj)
            integer :: quad_pos
            integer, intent(in) :: isymai, isymbj, nai, nbj
            if (isymai .lt. isymbj) then
               quad_pos = IT2AM(ISYMAI,ISYMBJ)
     &                  + nt1am(isymai)*(nbj-1) + nai
            else
               quad_pos = it2am(isymbj,isymai)
     &                  + nt1am(isymbj)*(nai-1) + nbj
            end if
            return
         end function

      END
C
      FUNCTION SO_DOUBLES_PRODUCT(X2AM,Y2AM,ISYMX2)
C
C     Returns the inner product of two singlet adapted doubles vectors
C
C     This is given as
C                           ~
C     0.5 * sum_{(ai)<(bj)} y(aibj)*x(aibj) + sum_{ai} y(aiai)*x(aiai)
C          ~
C     with y(aibj) = 2*y(aibj) - y(ajbi)
C     Inserting this in the above expression followed by some
C     manipulations of the first sum leads to
C                          ~                 ~
C     0.5 * sum_{a<b,i<j} (y(aibj)*x(aibj) + y(ajbi) * x(ajbi) )
C     + 0.5*sum_{a<b,i} y(aibi)*x(aibi) + 0.5*sum{a,i<j} y(aiaj)*x(aiaj)
C     + sum_{ai} y(aiai)*x(aiai)
C
C     which is implemented in the present routine.

      use so_info, only: sop_dp
      implicit none
C
      real(sop_dp),parameter :: ZERO = 0.0D0, HALF = 0.5D0, two = 2.0D0
C
      real(sop_dp), intent(in) :: X2AM(*), Y2AM(*)
      integer, intent(in)      :: isymx2
      real(sop_dp)             :: so_doubles_product
C
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
C
      INTEGER :: ISYMAB, ISYMIJ, ISYMAI, ISYMAJ, ISYMBJ, ISYMBI,
     &           ISYMA, ISYMB, ISYMI, ISYMJ
      INTEGER :: NAI, NBJ, NAJ, NBI, NAIAI, NAIBJ, NAJBI,
     &           NAJBJ, NBIBJ, NBJBJ
      INTEGER :: NRHFI, NVIRA, INDEX
      real(sop_dp) :: XAIBJ, XAJBI, YAIBJ, YAJBI
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      SO_DOUBLES_PRODUCT = ZERO

      IF (ISYMX2 .eq. 1) THEN

         do isymbj = 1, nsym
            isymai = isymbj

            do isymj = 1, nsym
               isymb = muld2h(isymj,isymbj)

               ! Ensure isymi <= isymj
               do isymi = 1, isymj
                  isyma = muld2h(isymi,isymai)
                  if (isyma.gt.isymb) cycle

                  isymaj = muld2h(isyma,isymj)
                  isymbi = muld2h(isymb,isymi)

                  do J = 1, nrhf(isymj)
                     ! Ensure i<j
                     if (isymi .eq. isymj) then
                        nrhfi = j-1
                     else
                        nrhfi = nrhf(isymi)
                     endif

                     do b = 1, nvir(isymb)
                        ! Ensure a<b
                        if (isyma.eq.isymb) then
                           nvira = b-1
                        else
                           nvira = nvir(isyma)
                        endif
                        nbj = pair_pos(isymb,b,isymj,j)
C
                        ! Do a<b
                        do i = 1, nrhfi
                           nbi = pair_pos(isymb,b,isymi,i)
                           do a = 1, nvira
                              nai = pair_pos(isyma,a,isymi,i)
                              naj = pair_pos(isyma,a,isymj,j)
                              naibj = it2am(isymai,isymbj) +
     &                                nbj*(nbj-1)/2 + nai
                              najbi = it2am(isymbi,isymaj) +
     &                                naj*(naj-1)/2 + nbi
                              xaibj = x2am(naibj)
                              xajbi = x2am(najbi)
                              yaibj = y2am(naibj)
                              yajbi = y2am(najbi)
                              SO_DOUBLES_PRODUCT = SO_DOUBLES_PRODUCT
     &                                       + (TWO*YAIBJ-YAJBI)*XAIBJ
     &                                       + (TWO*YAJBI-YAIBJ)*XAJBI
                           end do ! loop a
C
                           ! handle explicitly a == b, i<j
                           if (isyma.eq.isymb) then
                              nbibj = it2am(isymbi,isymbj) +
     &                                nbj*(nbj-1)/2 + nbi
                              so_doubles_product = so_doubles_product
     &                                      + x2am(nbibj)*y2am(nbibj)
                           end if
                        end do ! loop i

                        ! handle explicitly a <= b, i == j
                        if (isymi.eq.isymj) then
                           ! a<b
                           do a = 1, nvira
                              naj = pair_pos(isyma,a,isymj,j)
                              najbj = it2am(isymaj,isymbj) +
     &                                nbj*(nbj-1)/2 + naj
                              so_doubles_product = so_doubles_product
     &                                   + x2am(najbj)*y2am(najbj)
                           end do
                           ! a==b, i==j
                           nbjbj = it2am(isymbj,isymbj) +
     &                             nbj*(nbj+1)/2
                           so_doubles_product = so_doubles_product
     &                                    + two*x2am(nbjbj)*y2am(nbjbj)
                        end if
                     end do ! loop b
                  end do ! loop j
               end do ! loop isymi
            end do ! loop isymj
         end do ! loop isymbj
      else ! not totally symmetric case
C
         do isymbj = 1, nsym
            isymai = muld2h(isymbj,isymx2)

            do isymj = 1, nsym
               isymb = muld2h(isymj,isymbj)

               do isymi = 1, isymj
                  isyma = muld2h(isymi,isymai)
                  if (isyma.gt.isymb) cycle

                  isymaj = muld2h(isyma,isymj)
                  isymbi = muld2h(isymb,isymi)

                  do J = 1, nrhf(isymj)
                     ! Ensure i<j
                     if (isymi .eq. isymj) then
                        nrhfi = j-1
                     else
                        nrhfi = nrhf(isymi)
                     endif

                     do b = 1, nvir(isymb)
                        ! Ensure a<b
                        if (isyma.eq.isymb) then
                           nvira = b-1
                        else
                           nvira = nvir(isyma)
                        endif
                        nbj = pair_pos(isymb,b,isymj,j)
C
                        ! Do a<b
                        do i = 1, nrhfi
                           nbi = pair_pos(isymb,b,isymi,i)
                           do a = 1, nvira
                              nai = pair_pos(isyma,a,isymi,i)
                              naj = pair_pos(isyma,a,isymj,j)
                              naibj = quad_pos(isymai,nai,isymbj,nbj)
                              najbi = quad_pos(isymaj,naj,isymbi,nbi)
                              xaibj = x2am(naibj)
                              xajbi = x2am(najbi)
                              yaibj = y2am(naibj)
                              yajbi = y2am(najbi)
                              SO_DOUBLES_PRODUCT = SO_DOUBLES_PRODUCT
     &                                       + (TWO*YAIBJ-YAJBI)*XAIBJ
     &                                       + (TWO*YAJBI-YAIBJ)*XAJBI
                           end do ! loop a
C
                           ! handle explicitly a == b, i<j
                           if (isyma.eq.isymb) then
                              nbibj = quad_pos(isymbi,nbi,isymbj,nbj)
                              so_doubles_product = so_doubles_product
     &                                       + x2am(nbibj)*y2am(nbibj)
                           end if
                        end do ! loop i

                        ! handle explicitly a < b, i == j
                        if (isymi.eq.isymj) then
                           ! a<b
                           do a = 1, nvira
                              naj = pair_pos(isyma,a,isymj,j)
                              najbj = quad_pos(isymaj,naj,isymbj,nbj)
                              so_doubles_product = so_doubles_product
     &                                       + x2am(najbj)*y2am(najbj)
                           end do
                        end if
                     end do ! loop b
                  end do ! loop j
               end do ! loop isymi
            end do ! loop isymj
         end do ! loop isymbj
C
      end if
C
C     In the above calculations the factor of one half has been ignored
C     apply it now.
C
      SO_DOUBLES_PRODUCT = 0.5D0*SO_DOUBLES_PRODUCT

      contains
         pure function pair_pos(isyma,na,isymi,ni)
            integer :: pair_pos
            integer, intent(in) :: na, ni, isyma, isymi

            pair_pos = it1am(isyma,isymi) + nvir(isyma)*(ni-1) + na
            return
         end function

         pure function quad_pos(isymai,nai,isymbj,nbj)
            integer :: quad_pos
            integer, intent(in) :: isymai, isymbj, nai, nbj
            if (isymai .lt. isymbj) then
               quad_pos = IT2AM(ISYMAI,ISYMBJ)
     &                  + nt1am(isymai)*(nbj-1) + nai
            else
               quad_pos = it2am(isymbj,isymai)
     &                  + nt1am(isymbj)*(nai-1) + nbj
            end if
            return
         end function

      END
