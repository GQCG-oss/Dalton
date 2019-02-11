!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
!
*======================================================================*
      subroutine cc_fckdela(ibasd,isydel,fock,isyfck,xcou,xexc,ifckvao)
*----------------------------------------------------------------------*
*  Purpose: update Fock matrix with one AO and one virtual index
*           using precomputed partially transformed integrals
*  C. Haettig, spring 2006
*----------------------------------------------------------------------*
      implicit none
#include "ccsdsym.h"
#include "ccorb.h"

      real*8  one, two
      parameter ( one=1.0d0, two=2.0d0 )

* input:
      integer ibasd, isydel, isyfck, ifckvao(8,8)
      real*8  fock(*), xcou(*), xexc(*)

* local:
      integer isyma, kofff, isymi, koffx, isymai

      isyma = muld2h(isyfck,isydel)
      kofff = ifckvao(isyma,isydel) + nvir(isyma)*(ibasd-1) + 1

      do isymi = 1, nsym
        isymai = muld2h(isyma,isymi)
        do i = 1, nrhf(isymi)

          ! address of X^del(1,i,i)
          koffx = it2bcd(isymai,isymi) + nt1am(isymai)*(i-1) + 
     &              it1am(isyma,isymi) + nvir(isyma)*(i-1)   + 1

          call daxpy(nvir(isyma), two,xcou(koffx),1,fock(kofff),1)
          call daxpy(nvir(isyma),-one,xexc(koffx),1,fock(kofff),1)

        end do
      end do

      return
      end 
*======================================================================*
