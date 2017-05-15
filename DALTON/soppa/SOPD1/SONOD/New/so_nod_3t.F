C
C  /* Deck so_nod */
      subroutine so_nod_3t(xint,lxint,xlamb,lxlamb,delta)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, November 1995
C     Mats Simmermacher & Stephan P. A. Sauer, December 2014
C
C     PURPOSE: Calculation of the one and three times transformed
C              two-electron integrals (alpha beta | c delta) and
C              (a b | c delta).
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      dimension xint(lxint),      xlamb(lxlamb)
      dimension oint(10,10,10),   aint(10,10,10)
      integer   alfa,  beta, gamma, delta
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      index(i,j) = max(i,j)*(max(i,j) - 3)/2 + i + j
C
      call dzero(aint,1000)
      DO 100 gamma = 1,nbas(1)
         l = index(nbas(1),nbas(1))*(gamma-1)
         DO 200 beta = 1,nbas(1)
            DO 300 alfa = 1,nbas(1)
               m = l + index(alfa,beta)
               aint(alfa,beta,gamma) = xint(m)
  300       continue
  200    continue
  100 continue
C
      call dzero(oint,1000)
      do 400 c = 1,nbas(1)
         do 500 gamma = 1,nbas(1)
            do 600 beta = 1,nbas(1)
               do 700 alfa = 1,nbas(1)
                  koff = (c-1)*nbas(1) + gamma
                  oint(alfa,beta,c) = oint(alfa,beta,c)
     &                              + aint(alfa,beta,gamma)*xlamb(koff)
  700          continue
  600       continue
  500    continue
  400 continue
C
      write (LUPRI,'(A)')
     &       '  (alpha beta | c delta) integrals calculated with NOD3T'
      do 101 c = 1,nbas(1)
         do 102 beta = 1,nbas(1)
            do 103 alfa = 1,nbas(1)
               if (oint(alfa,beta,c) .NE. ZERO) then
               write (LUPRI,'(4(I2,1X),2X,E14.7)')
     &                delta,c,beta,alfa,oint(alfa,beta,c)
               end if
  103       continue
  102    continue
  101 continue
C
      call dzero(aint,1000)
      do 800 c = 1,nbas(1)
         do 900 b = 1,nbas(1)
            do 1000 beta = 1,nbas(1)
               do 1100 alfa = 1,nbas(1)
                  koff = (b-1)*nbas(1) + beta
                  aint(alfa,b,c) = aint(alfa,b,c)
     &                              + oint(alfa,beta,c)*xlamb(koff)
 1100          continue
 1000       continue
  900    continue
  800 continue
C
      call dzero(oint,1000)
      do 1200 c = 1,nbas(1)
         do 1300 b = 1,nbas(1)
            do 1400 a = 1,nbas(1)
               do 1500 alfa = 1,nbas(1)
                  koff = (a-1)*nbas(1) + alfa
                  oint(a,b,c)   = oint(a,b,c)
     &                          + aint(alfa,b,c)*xlamb(koff)
 1500          continue
 1400       continue
 1300    continue
 1200 continue
C
      write (LUPRI,'(A)')
     &       '  (a b | c delta) integrals calculated with NOD3T'
      do 104 c = 1,nbas(1)
         do 105 b = 1,nbas(1)
            do 106 a = 1,nbas(1)
               if (oint(a,b,c) .NE. ZERO) then
               write (LUPRI,'(4(I2,1X),2X,E14.7)')
     &                delta,c,b,a,oint(a,b,c)
               end if
  106       continue
  105    continue
  104 continue
C
      write (LUPRI,'(A)')
     &       '  c_(delta d) MO coefficients calculated with NOD3T'
       do 107 d = 1,nbas(1)
          koff2 = (d-1)*nbas(1) + delta
             if (xlamb(koff2) .NE. ZERO) then
             write (LUPRI,'(2(I2,1X),2X,E14.7)')
     &              delta,d,xlamb(koff2)
             end if
  107 continue
C
      return
      end
