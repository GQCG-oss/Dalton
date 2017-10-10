C
C  /* Deck so_nod */
      subroutine so_nod(xint,lxint,xlamb,lxlamb,delta,rint)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, November 1995
C
C     PURPOSE: Nody code to calculate third and fourth term
C              of eq. (34).
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      dimension xint(lxint),      xlamb(lxlamb)
      dimension oint(10,10,10),   aint(10,10,10)
      dimension rint(10,10,10,10)
      integer   alfa,  beta, gamma, delta
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      index(i,j) = max(i,j)*(max(i,j) - 3)/2 + i + j
C
c      print*,'KELD delta = ',delta
c      print*,'KELD lxint = ',lxint
c      print*,'KELD nbas(1) = ',nbas(1)
c      print*,'KELD xint = ',(xint(i), i=1,lxint)
c
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
      do 400 gamma = 1,nbas(1)
         do 500 beta = 1,nbas(1)
            do 600 a = 1,nbas(1)
               do 700 alfa = 1,nbas(1)
                  koff = (a-1)*nbas(1) + alfa
                  oint(a,beta,gamma) = oint(a,beta,gamma)
     &                              + aint(alfa,beta,gamma)*xlamb(koff)
c      print*,'KELD xlamb = ',xlamb(koff)
  700          continue
  600       continue
  500    continue
  400 continue
C
      call dzero(aint,1000)
      do 800 gamma = 1,nbas(1)
         do 900 a = 1,nbas(1)
            do 1000 b = 1,nbas(1)
               do 1100 beta = 1,nbas(1)
                  koff = (b-1)*nbas(1) + beta
                  aint(a,b,gamma) = aint(a,b,gamma)
     &                              + oint(a,beta,gamma)*xlamb(koff)
c      print*,'KELD xlamb = ',xlamb(koff)
 1100          continue
 1000       continue
  900    continue
  800 continue
C
      call dzero(oint,1000)
      do 1200 a = 1,nbas(1)
         do 1300 b = 1,nbas(1)
            do 1400 c = 1,nbas(1)
               do 1500 gamma = 1,nbas(1)
                  koff = (c-1)*nbas(1) + gamma
                  oint(a,b,c)   = oint(a,b,c)
     &                          + aint(a,b,gamma)*xlamb(koff)
c      print*,'KELD xlamb = ',xlamb(koff)
 1500          continue
 1400       continue
 1300    continue
 1200 continue
C
      do 1600 a = 1,nbas(1)
         do 1700 b = 1,nbas(1)
            do 1800 c = 1,nbas(1)
               do 1900 d = 1,nbas(1)
                  koff2 = (d-1)*nbas(1) + delta
                  rint(a,b,c,d) = rint(a,b,c,d)
     &                          + oint(a,b,c)*xlamb(koff2)
c      print*,'KELD xlamb = ',xlamb(koff2)
 1900          continue
 1800       continue
 1700    continue
 1600 continue
C
      return
      end
