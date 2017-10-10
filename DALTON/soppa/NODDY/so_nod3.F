C
C  /* Deck so_nod3 */
      subroutine so_nod3(rint,work,lwork)
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
      do 3000 k = 1,nrhf(1)
         do 3001 c = 1,nvir(1)
            ic = nrhf(1) + c
            do 3002 j = 1,nrhf(1)
               if (j.eq.k) then
                  fac2 = one
               else
                  fac2 = zero
               endif
               do 3003 b = 1,nvir(1)
                  ib = nrhf(1) + b
                  if (ic.eq.ib) then
                     fac4 = one
                  else
                     fac4 = zero
                  endif
                  do 3004 i = 1,nrhf(1)
                     if (i.eq.k) then
                        fac1 = one
                     else
                        fac1 = zero
                     endif
                     do 3005 a = 1,nvir(1)
                        ia = nrhf(1) + a
                        if (ia.eq.ic) then
                           fac3 = one
                        else
                           fac3 = zero
                        endif
                        if ( (ia.eq.ib) .and.  (i.eq.j) ) then
                           fac = one / dsqrt(two)
                        else
                           fac = dsqrt(two)
                        endif
                        cmat = fac * (   fac1 * rint(ia,ic,ib,j)
     &                                 + fac2 * rint(ia,i,ib,ic)
     &                                 - fac3 * rint(k,i,ib,j)
     &                                 - fac4 * rint(ia,i,k,j) )
                        write(6,9999) 'Keld: C(',ia,i,ib,j,',',ic,k,
     &                                ') =', cmat
 3005                continue
 3004             continue
 3003          continue
 3002       continue
 3001    continue
 3000 continue
C
 9999 format(1x,a8,i1,i1,i1,i1,a1,i1,i1,a2,f14.8)
      return
      end
