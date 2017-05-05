C
C  /* Deck so_nod4 */
      subroutine so_nod4(rint,work,lwork)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, March 1996
C
C     PURPOSE: Nody code to calculate C matrix elements.
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      dimension work(lwork)
      dimension rint(10,10,10,10)
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      fac = one / dsqrt(two)
C
      do 3000 a = 1,nvir(1)
         ia = nrhf(1) + a
         do 3001 i = 1,nrhf(1)
            do 3002 b = 1,nvir(1)
               ib = nrhf(1) + b
               do 3003 j = 1,nrhf(1)
                  do 3004 c = 1,nvir(1)
                     ic = nrhf(1) + c
                     if (ia.eq.ic) then
                        fac3 = one
                     else
                        fac3 = zero
                     endif
                     if (ic.eq.ib) then
                        fac4 = one
                     else
                        fac4 = zero
                     endif
                     do 3005 k = 1,nrhf(1)
                        if (j.eq.k) then
                           fac2 = one
                        else
                           fac2 = zero
                        endif
                        if (i.eq.k) then
                           fac1 = one
                        else
                           fac1 = zero
                        endif

                        cmat = fac * (  fac1 * (two * rint(ia,ic,ib,j)
     &                                              - rint(ia,j,ib,ic))
     &                                + fac2 * (two * rint(ia,i,ib,ic)
     &                                              - rint(ia,ic,ib,i))
     &                                - fac3 * (two * rint(k,i,ib,j)
     &                                              - rint(k,j,ib,i))
     &                                - fac4 * (two * rint(k,j,ia,i)
     &                                              - rint(k,i,ia,j)) )
C
                        write(6,9998) 'Keld: g(',ia,j,k,i,') =',
     &                                 fac*rint(ia,j,k,i)
                        write(6,9998) 'Keld: g(',k,j,ib,i,') =',
     &                                 fac*rint(k,j,ib,i)
                        write(6,9999) 'Keld: C(',k,ic,',',j,ib,i,ia,
     &                                ') =', cmat
 3005                continue
 3004             continue
 3003          continue
 3002       continue
 3001    continue
 3000 continue
C
 9998 format(1x,a8,i1,i1,i1,i1,a3,f14.8)
 9999 format(1x,a8,i1,i1,a1,i1,i1,i1,i1,a2,f14.8)
      return
      end
