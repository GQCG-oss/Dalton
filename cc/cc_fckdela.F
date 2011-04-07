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
      subroutine cc_fckdela(ibasd,isydel,fock,isyfck,xcou,xexc,ifckvao)
*----------------------------------------------------------------------*
*  Purpose: update Fock matrix with one AO and one virtual index
*           using precomputed partially transformed integrals
*  C. Haettig, spring 2006
*----------------------------------------------------------------------*
      implicit none
#include <ccsdsym.h>
#include <ccorb.h>

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
