C
C  /* Deck so_firgp */
      SUBROUTINE SO_FIRGP(GPVC2,LGPVC2,T2AM,LT2AM,PRP1,LPRP1,PR1IJ,
     &                    LPR1IJ,PR1AB,LPR1AB,ISYMTR,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, July 1997
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Calculate 2p-2h part of gradient property vectors.
C              There is only a first order contribution.
C
#include "implicit.h"
C
#include "ccsdsym.h"
#include "ccorb.h"
#include "soppinf.h"
C
      PARAMETER   (ONE = 1.0D0, TWO = 2.0D0, HALF = 0.5D0)
      DIMENSION   GPVC2(LGPVC2), T2AM(LT2AM)
      DIMENSION   PRP1(LPRP1),   PR1IJ(LPR1IJ), PR1AB(LPR1AB)
      DIMENSION   WORK(LWORK)
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_FIRGP')
C
C---------------------------------------------------
C     Repack the one-particle property mo-integrals.
C---------------------------------------------------
C
      CALL SO_RPPRP1(PRP1,LPRP1,PR1IJ,LPR1IJ,PR1AB,LPR1AB,ISYMTR)
C
C------------------------------
C     Allocation of work space.
C------------------------------
C
      LT2SQ  = NT2SQ(1)
C
      KT2SQ  = 1
      KEND1  = KT2SQ + LT2SQ
      LWORK1 = LWORK - KEND1
C
      IF (LWORK1 .LT. 0) CALL STOPIT('SO_FIRGP',' ',KEND1,LWORK)
C
C
      IF (TRIPLET) THEN
C 
C     We need to change basis on lt2am
C        
         lt2amh = lt2am
         KT2CP = kend1
         kend2 = kt2cp + lt2amh
         lwork2 = lwork - kend2 
         IF (LWORK2 .LT. 0) CALL STOPIT('SO_FIRGP',' ',KEND2,LWORK)
         CALL DCOPY(LT2AMH,T2AM,1,WORK(kt2cp),1)
         call SO_TMLTR(WORK(kt2cp),ONE,1)
         CALL CC_T2SQ(work(kt2cp),work(kt2sq),1)
C
         CALL SO_GP2T(GPVC2,WORK(KT2SQ),PR1AB,PR1IJ,WORK(KEND1),LWORK1,
     &                ISYMTR)         


      ELSE
C
C---------------------------------
C     Square up the t2-amplitudes.
C---------------------------------
C
         CALL CC_T2SQ(T2AM,WORK(KT2SQ),1)
C
C-----------------------------------------------------------------------
C     Contract T2-amplitudes with one-particle property integrals and
C     scale with a half to calculate the 2p-2h gradient property vector.
C-----------------------------------------------------------------------
C
         CALL CCRHS_E(GPVC2,WORK(KT2SQ),PR1AB,PR1IJ,WORK(KEND1),LWORK1,
     &                1,ISYMTR)
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      END IF
      CALL QEXIT('SO_FIRGP')
C
      RETURN
      END

      SUBROUTINE SO_GP2T(GPVC2,T2SQ,PR1AB,PR1IJ,WORK,LWORK,ISYMTR)
     &               

      use so_info, only: sop_dp
      implicit none
#include "priunit.h"
#include "ccsdsym.h"
#include "ccorb.h"
#include "soppinf.h"
      character(len=*), parameter :: myname = 'SO_GPT2'
C     Arguments
      integer, intent(in) :: isymtr, lwork
      real(sop_dp),intent(out) :: GPVC2(N2P2HOP(ISYMTR))
      real(sop_dp),intent(in) :: T2SQ(NT2SQ(1)), PR1AB(*), PR1IJ(*)
      real(sop_dp),intent(inout) :: work(lwork)

      real(sop_dp), parameter :: one = 1.0D0, onem = -1.0D0, 
     &                           zero = 0.0D0,
     &                           fact1 = -ONE/SQRT(2.0D0),
     &                           fact2 = -ONE/(2.0D0), 
     &                           fact3 = -ONE/(2.0D0),
     %                           fact2d= -ONE/SQRT(2.0D0),
     &                           fact3d= -ONE/SQRT(2.0D0)

      integer :: ktemp, ltemp, kend
      ! Irrep numbers
      integer :: isymab, isymai, isymaj, isymbj, isymbi, isymij,
     &           isymei, isymam,
     &           isyma, isymb, isyme, isymm, isymi, isymj
      ! Pair numbers 
      integer :: nai, naj, nbi, nbj
      ! Quad numbers
      integer :: naibj, najbi, nbiaj, nbjai, najbj, nbjaj, nbibj, nbjbi
      ! Memory positions
      integer :: kofft, koffp, koffy, koff1, koff2, koff3
      integer :: nvira, nvire, nrhfi, nrhfm

      ltemp = NT2SQ(ISYMTR) 
      ktemp = 1
      kend = ktemp + ltemp
      IF (KEND.gt.lwork) CALL STOPIT('myname',' ',KEND,LWORK)
     
C      print *, 'get_gp' 

      do isymbj = 1,nsym
         isymei = isymbj
         isymai = muld2h(isymbj,isymtr)
         isymam = isymbj

         do nbj = 1, nt1am(isymbj)
            do isymi = 1, nsym
               isyma = muld2h(isymai,isymi)
               isyme = muld2h(isymei,isymi)
               nvire = max(1,nvir(isyme))
               nvira = max(1,nvir(isyma))
               kofft = IT2SQ(isymei,isymbj) + nt1am(isymei)*(nbj-1)
     &               + it1am(isyme,isymi) + 1
               koffp = imatab(isyme,isyma) + 1
               koffy = it2sq(isymai,isymbj) + nt1am(isymai)*(nbj-1)
     &               + it1am(isyma,isymi) + ktemp
C
C              \sum P_{ae}*t(eibj) => y(aibj)
C              Calculated instead as
C           
C              (-) \sum P_{ea}*t(eibj)              

               CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYME),
     &                    ONE,pr1ab(koffp),nvire,t2sq(kofft),nvire,
     &                    zero,work(koffy),nvira)

               isymm = muld2h(isyma,isymam)
               nrhfm = max(1,nrhf(isymm))
               kofft = IT2SQ(isymam,isymbj) + nt1am(isymam)*(nbj-1)
     &               + it1am(isyma,isymm) + 1
               koffp = imatij(isymm,isymi) + 1
C
C              - \sum t(ambj)*P_{mi} => y(aibj)
C               
               CALL DGEMM('N','N',nvir(isyma),nrhf(isymi),nrhf(isymm),
     &                    onem,t2sq(kofft),nvira,pr1ij(koffp),nrhfm,
     &                    one,work(koffy),nvira)
            end do ! isymi
         end do ! nbj
      end do ! isymbj

C      print *, 'First half done'
C
C     Reshuffle y(aibj) on to the triplet gp2 vector
C
C     {T1}g(abij) = y(aibj) - y(ajbi) - y(biaj) + y(bjai)
C     {T2}g(abij) = y(aibj) - y(ajbi) + y(biaj) - y(bjai)
C     {T3}g(abij) = y(aibj) + y(ajbi) - y(biaj) - y(bjai)
C
      koff1 = 0
      koff2 = NT2AMT1(isymtr)
      koff3 = koff2 + nt2amt2(isymtr)
      do isymij = 1, nsym
         isymab = muld2h(isymtr,isymij)
         do isymj = 1, nsym
            isymi = muld2h(isymj,isymij)
            if (isymi.gt.isymj) cycle
            do j = 1, nrhf(isymj)
               if (isymi.eq.isymj) then
                  nrhfi = j-1
               else 
                  nrhfi = nrhf(isymi)
               endif
               ! loop i<j
               do i = 1, nrhfi
                  do isymb = 1, nsym
                     isyma = muld2h(isymab,isymb)
                     if (isyma .gt. isymb) cycle
                     isymai = muld2h(isyma,isymi)
                     isymbj = muld2h(isymb,isymj)
                     isymaj = muld2h(isyma,isymj)
                     isymbi = muld2h(isymb,isymi)
                     if(isyma.eq.isymb) then
                        do b = 1, nvir(isymb)
                           nbi = pair_pos(isymb,b,isymi,i)
                           nbj = pair_pos(isymb,b,isymj,j)
                           do a = 1, b-1
                              koff1 = koff1+1
                              koff2 = koff2+1
                              koff3 = koff3+1
!                              print *, 2, koff1, koff2, koff3
                              nai = pair_pos(isyma,a,isymi,i)
                              naj = pair_pos(isyma,a,isymj,j)
                              naibj = quad_pos(isymai,nai,isymbj,nbj)
                              najbi = quad_pos(isymaj,naj,isymbi,nbi)
                              nbiaj = quad_pos(isymbi,nbi,isymaj,naj)
                              nbjai = quad_pos(isymbj,nbj,isymai,nai)

                              gpvc2(koff1) = fact1*(
     &                                       work(naibj) - work(najbi)
     &                                      -work(nbiaj) + work(nbjai) )
                              gpvc2(koff2) = fact2*(
     &                                       work(naibj) - work(najbi)
     &                                      +work(nbiaj) - work(nbjai) )
                              gpvc2(koff3) = fact3*(
     &                                       work(naibj) + work(najbi)
     &                                      -work(nbiaj) - work(nbjai) )
                           end do ! a
                           ! handle a==b (case 2 only)
                           koff2 = koff2 + 1
!                           print *, 1, koff1, koff2, koff3
                           nbibj = quad_pos(isymbi,nbi,isymbj,nbj)
                           nbjbi = quad_pos(isymbj,nbj,isymbi,nbi)
C                           print *, work(nbibj), work(nbjbi)
                           gpvc2(koff2) = fact2d*
     &                                    (work(nbibj)-work(nbjbi))
                        end do ! b
                     else ! isyma .lt. isymb
                        do b = 1, nvir(isymb)
                           nbi = pair_pos(isymb,b,isymi,i)
                           nbj = pair_pos(isymb,b,isymj,j)
                           do a = 1, nvir(isyma)
                              koff1 = koff1+1
                              koff2 = koff2+1
                              koff3 = koff3+1
!                              print *, 3, koff1, koff2, koff3
                              nai = pair_pos(isyma,a,isymi,i)
                              naj = pair_pos(isyma,a,isymj,j)
                              naibj = quad_pos(isymai,nai,isymbj,nbj)
                              najbi = quad_pos(isymaj,naj,isymbi,nbi)
                              nbiaj = quad_pos(isymbi,nbi,isymaj,naj)
                              nbjai = quad_pos(isymbj,nbj,isymai,nai)

                              gpvc2(koff1) = fact1*(
     &                                       work(naibj) - work(najbi)
     &                                      -work(nbiaj) + work(nbjai) )
                              gpvc2(koff2) = fact2*(
     &                                       work(naibj) - work(najbi)
     &                                      +work(nbiaj) - work(nbjai) )
                              gpvc2(koff3) = fact3*(
     &                                       work(naibj) + work(najbi)
     &                                      -work(nbiaj) - work(nbjai) )
                           end do ! a
                        end do ! b
                     end if ! isyma/isymb 
                  end do ! isymb
               end do ! i
               ! handle i == j (case 3 only)
               if (isymi.eq.isymj) then
                  do isymb = 1, nsym
                     isyma = muld2h(isymab,isymb)
                     isymbj = muld2h(isymj,isymb)
                     isymaj = muld2h(isymj,isyma)
                     if (isyma.gt.isymb) cycle
                     do b = 1, nvir(isymb)
                        nbj = pair_pos(isymb,b,isymj,j)
                        if (isyma.eq.isymb) then
                           nvira = b - 1
                        else
                           nvira = nvir(isyma)
                        end if
                        do a = 1, nvira
                           koff3 = koff3+1
!                           print *, 4, koff1, koff2, koff3
                           naj = pair_pos(isyma,a,isymj,j)
                           najbj = quad_pos(isymaj,naj,isymbj,nbj)
                           nbjaj = quad_pos(isymbj,nbj,isymaj,naj)
                           gpvc2(koff3) = fact3d
     &                                    *(work(najbj)-work(nbjaj))
                        end do ! a
                     end do ! b
                  end do ! isymb
               end if ! isymi .eq. isymj
            end do ! j
         end do ! isymj
      end do ! isymij   
C
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
            quad_pos = IT2SQ(ISYMAI,ISYMBJ)
     &               + nt1am(isymai)*(nbj-1) + nai
            return
         end function



      end subroutine
