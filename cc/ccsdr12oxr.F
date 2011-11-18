!
!...   Copyright (c) 2011 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2011 (2011), see http://daltonprogram.org"
!...
!...   This source code is provided under a written licence and may be
!...   used, copied, transmitted, or stored only in accord with that
!...   written licence.
!...
!...   In particular, no part of the source code or compiled modules may
!...   be distributed outside the research group of the licence holder.
!...   This means also that persons (e.g. post-docs) leaving the research
!...   group of the licence holder may not take any part of Dalton,
!...   including modified files, with him/her, unless that person has
!...   obtained his/her own licence.
!...
!...   For further information, including how to get a licence, see:
!...      http://daltonprogram.org
!
!
*======================================================================*
      subroutine ccsdr12oxr(omegpk,ombjpk,isymom,work,lwork)
*----------------------------------------------------------------------*
* Purpose: contract the intermediate result vector Omega(bj,p'k) with
*          r12 integrals to the final R12 result vector:
*
*          Omega^jk_mn += P^jk_mn sum_ap' r12(am,p'n) Omega(aj,p'k)
*
* C. Haettig, spring 2006
*----------------------------------------------------------------------*
      implicit none
#include <priunit.h>
#include <dummy.h>     
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
#include <ccr12int.h>

      logical test, locdbg
      parameter ( test = .false. )
      parameter ( locdbg = .false. )

      real*8  zero,half,one,two
      parameter ( zero=0.0d0, half=0.5d0, one=1.0d0, two=2.0d0 )
      character*8 framnp
      parameter(framnp='R12RMNAP')


* input:
      integer lwork, isymom
      real*8  omegpk(*), ombjpk(*), work(*)

* local:
      logical ldum, lproj
      integer idum, lunit, komegsq, kend1, lwrk1, 
     &        isymp, krbmn, kend2, lwrk2, iadr, isymk, isympk, idxpk,
     &        isyjmn, isymj, idxjk, isymmn, isymb, nvirb, nr12mn,
     &        kbjpk, koffr, kmnjk, isym, isym1, isym2, iramn(8,8),
     &        nramn(8), iramnq(8,8), nramnq(8), isymbj, isymjk
      integer isymm,isymn,isymmj,isymnk,idxmj,idxnk,idxmn,kmnj,idxmjnk
      real*8  fac, ddot

      integer index
      index(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j

*----------------------------------------------------------------------*
* initializations:
*----------------------------------------------------------------------*
CCN     do isym = 1, nsym
CCN       nramn(isym) = 0
CCN       nramnq(isym) = 0
CCN       do isym2 = 1, nsym
CCN         isym1 = muld2h(isym,isym2)
CCN         iramn(isym1,isym2) = nramn(isym)
CCN         nramn(isym) = nramn(isym) + nvir(isym1)*nmatkl(isym2)
CCN         iramnq(isym1,isym2) = nramnq(isym)
CCN         nramnq(isym) = nramnq(isym) + nramn(isym1)*norb2(isym2)
CCN       end do
CCN     end do

CCN   !has to be done in two rounds:
      do isym = 1, nsym
        nramn(isym) = 0
        do isym2 = 1, nsym
          isym1 = muld2h(isym,isym2)
          iramn(isym1,isym2) = nramn(isym)
          nramn(isym) = nramn(isym) + nvir(isym1)*nmatkl(isym2)
        end do
      end do
C
      do isym = 1, nsym
        nramnq(isym) = 0
        do isym2 = 1, nsym
          isym1 = muld2h(isym,isym2)
          iramnq(isym1,isym2) = nramnq(isym)
          nramnq(isym) = nramnq(isym) + nramn(isym1)*norb2(isym2)
        end do
      end do

*----------------------------------------------------------------------*
* allocate work space for final R12 result vector in unpacked 
* rectangular form and initialize with present vector on omegpk:
*----------------------------------------------------------------------*

      komegsq = 1
      kend1   = komegsq + ntr12sq(isymom)
      lwrk1   = lwork - kend1
      if (lwrk1.lt.0) call quit('Insufficient memory in ccsdr12oxr')


      if (test) then
        fac = zero
        call dzero(work(komegsq),ntr12sq(isymom))
      else
        fac = one

        ! scale with 1/2 to compensate factor 2 of symmetrization below
        call dscal(ntr12am(isymom),half,omegpk,1)

        ! unpack to rectangular matrix O(mn,jk)
        call ccr12unpck2(omegpk,isymom,work(komegsq),'N',1)
      end if

*----------------------------------------------------------------------*
* loop over orthogonal auxiliary basis functions:
*----------------------------------------------------------------------*
      ! open file with r12 integrals
      lunit = -1
      call wopen2(lunit,framnp,64,0)

      do isymp = 1, nsym
        do p = 1, norb2(isymp)

          krbmn = kend1
          kend2 = krbmn + nramn(isymp)
          lwrk2 = lwork - kend2
          if (lwrk2.lt.0) call quit('Insufficient memory in ccsdr12oxr')

          ! read (bm|r12|p'n) integrals for this p' sorted as r(b,mn)
          iadr = iramnq(isymp,isymp) + nramn(isymp)*(p-1) + 1
          call getwa2(lunit,framnp,work(krbmn),iadr,nramn(isymp))

          if (locdbg) then
          write(lupri,*) "isymp, p', iadr: ", isymp, p, iadr
          write(lupri,*) "Norm^2 of r(b,mn)^p': ", 
     &      ddot(nramn(isymp),work(krbmn),1,work(krbmn),1)
          end if 

          do isymk = 1, nsym
            isympk = muld2h(isymp,isymk)
            do k = 1, nrhf(isymk)
              idxpk = ig1am(isymp,isymk)+norb2(isymp)*(k-1)+p

              isyjmn = muld2h(isymom,isymk)
              do isymj = 1, nsym
                idxjk = imatij(isymj,isymk) + nrhf(isymj)*(k-1) + 1

                isymmn = muld2h(isyjmn,isymj)
                isymb  = muld2h(isymmn,isymp)
                isymbj = muld2h(isymb,isymj)
                isymjk = muld2h(isymk,isymj)

                nvirb  = max(nvir(isymb),1)
                nr12mn = max(nmatkl(isymmn),1)

                kbjpk = itg2sq(isymbj,isympk) + 
     &             nt1am(isymbj)*(idxpk-1) + it1am(isymb,isymj) + 1
                koffr = krbmn + iramn(isymb,isymmn)
                kmnjk = komegsq + itr12sq(isymmn,isymjk) +
     &                    nmatkl(isymmn)*(idxjk-1)

                if (locdbg) then
                write(lupri,*) 'isymk, k: ', isymk, k
                write(lupri,*) 'Norm^2 of ombjpk(kbjpk): ',
     &            ddot(nvir(isymb)*nrhf(isymj),ombjpk(kbjpk),1,
     &                                         ombjpk(kbjpk),1) 
                end if

                call dgemm('T','N',nmatkl(isymmn),nrhf(isymj),
     &                             nvir(isymb),
     &                     one,work(koffr),nvirb,ombjpk(kbjpk),nvirb,
     &                     fac,work(kmnjk),nr12mn)

                if (locdbg) then
                write(lupri,*) 'Norm^2 of omegsq: ',
     &            ddot(ntr12sq(isymom),work(komegsq),1,
     &                                 work(komegsq),1)
                call cc_prsqr12(work(komegsq),isymom,'N',1,.false.)
                write(lupri,*) 'r-ints used:'
                call output(work(koffr),
     &                      1,nvir(isymb),1,nmatkl(isymmn),
     &                      nvir(isymb),nmatkl(isymmn),1,lupri)
                write(lupri,*) 'omega used:'
                call output(ombjpk(kbjpk),
     &                      1,nvir(isymb),1,nrhf(isymj),
     &                      nvir(isymb),nrhf(isymj),1,lupri)
                end if
              
c               -------------------------------------------------
c               only for testing... to be removed later...
c               ... sort directly into the packed result array...
c               -------------------------------------------------
                if (test) then
                 do isymm = 1, nsym
                  isymn  = muld2h(isymmn,isymm)
                  isymmj = muld2h(isymm,isymj)
                  isymnk = muld2h(isymn,isymk)
               
               
                  do j = 1, nrhf(isymj)
                   do m = 1, nrhfb(isymm)
                    do n = 1, nrhfb(isymn)
                      idxmj=imatki(isymm,isymj)+nrhfb(isymm)*(j-1)+m
                      idxnk=imatki(isymn,isymk)+nrhfb(isymn)*(k-1)+n
                      idxmn=imatkl(isymm,isymn)+nrhfb(isymm)*(n-1)+m
                     
                      kmnj = kmnjk-1 + nmatkl(isymmn)*(j-1) + idxmn

                      if (isymmj.eq.isymnk) then
                        idxmjnk = itr12am(isymmj,isymnk) +
     &                               index(idxmj,idxnk)
                      else if (isymmj.lt.isymnk) then
                        idxmjnk = itr12am(isymmj,isymnk) +
     &                               nmatki(isymmj)*(idxnk-1) + idxmj
                      else if (isymnk.lt.isymmj) then
                        idxmjnk = itr12am(isymmj,isymnk) +
     &                               nmatki(isymnk)*(idxmj-1) + idxnk
                      else
                        call quit('Error in ccsdr12oxr')
                      end if

                      if (idxmj.eq.idxnk) work(kmnj) = two * work(kmnj) 

                      omegpk(idxmjnk) = omegpk(idxmjnk) + work(kmnj)

                    end do
                   end do
                  end do
               
                 end do
                end if ! (test)

              end do

            end do
          end do

        end do
      end do

      call wclose2(lunit,framnp,'KEEP')
*----------------------------------------------------------------------*
* end loop over orthogonal auxiliary basis functions p'...
* resort final result vector back to triangular storage scheme:
*----------------------------------------------------------------------*

      if (.not.test) then
        ! symmetrize and pack result vector to triangular matrix 
        !    O(mj,nk) <-- O(mn,jk) + O(nm,kj)
        lproj = .true.
        call ccr12pck2(omegpk,isymom,lproj,work(komegsq),'N',1)
      end if

      if (locdbg) then
        call cc_prpr12(omegpk,isymom,1,.false.)
      end if

      return
      end 
*----------------------------------------------------------------------*
*                     END OF SUBROUTINE CCSDR12OXR                     *
*======================================================================*
