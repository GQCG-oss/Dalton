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

#if defined (SYS_CRAY)
      real one, two
#else
      double precision one, two
#endif
      parameter ( one=1.0d0, two=2.0d0 )

* input:
      integer ibasd, isydel, isyfck, ifckvao(8,8)
#if defined (SYS_CRAY)
      real fock(*), xcou(*), xexc(*)
#else
      double precision fock(*), xcou(*), xexc(*)
#endif

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
