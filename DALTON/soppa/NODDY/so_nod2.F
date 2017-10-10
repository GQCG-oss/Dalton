C
C  /* Deck so_nod2 */
      subroutine so_nod2(tr1e,ltr1e,tr1d,ltr1d,t2mp,lt2mp,
     &                   rint,work,lwork)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, November 1995
C
C     PURPOSE: Nody code to calculate first, second, third and fourth
C              term of eq. (34).
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION TR1E(LTR1E), TR1D(LTR1D)
      DIMENSION T2MP(LT2MP)
      DIMENSION WORK(LWORK)
      dimension rint(10,10,10,10)
      dimension g1(10,10,10,10)
      dimension g2(10,10,10,10)
      dimension g3(10,10,10,10)
      dimension g4(10,10,10,10)
      dimension g5(10,10,10,10)
      dimension g6(10,10,10,10)
      dimension g7(10,10,10,10)
      dimension g8(10,10,10,10)
      dimension uai1(10,10)
      dimension uai2(10,10)
      dimension uai3(10,10)
      dimension uai4(10,10)
      dimension uai5(10,10)
      dimension uai6(10,10)
      dimension uai7(10,10)
      dimension uai8(10,10)
      dimension fai(10,10)
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      index(i,j) = max(i,j)*(max(i,j) - 3)/2 + i + j
C
      call dzero(g1,10000)
      call dzero(g2,10000)
      call dzero(g3,10000)
      call dzero(g4,10000)
      call dzero(g5,10000)
      call dzero(g6,10000)
      call dzero(g7,10000)
      call dzero(g8,10000)
      call dzero(uai1,100)
      call dzero(uai2,100)
      call dzero(uai3,100)
      call dzero(uai4,100)
      call dzero(uai5,100)
      call dzero(uai6,100)
      call dzero(uai7,100)
      call dzero(uai8,100)
      call dzero(fai,100)
C
      do 101 a = 1,nvir(1)
         do 102 d = 1,nvir(1)
            do 103 k = 1,nrhf(1)
               do 104 c = 1,nvir(1)
C
                  do 105 j = 1,nrhf(1)
                     id  = nrhf(1) + d
                     ic  = nrhf(1) + c
                     ia  = nrhf(1) + a
                     naj = (j-1)*nvir(1) + a
                     njd = (j-1)*nvir(1) + d
                     g1(a,d,k,c) = g1(a,d,k,c)
     &                           + rint(j,id,k,ic)*tr1e(naj)
     &                           - rint(ia,j,k,ic)*tr1d(njd)
C
                     g5(d,a,c,k) = g5(d,a,c,k)
     &                           + rint(j,ia,ic,k)*tr1e(njd)
     &                           - rint(id,j,ic,k)*tr1d(naj)
  105             continue
C
  104          continue
  103       continue
  102    continue
  101 continue
C
      do 201 l = 1,nrhf(1)
         do 202 i = 1,nrhf(1)
            do 203 k = 1,nrhf(1)
               do 204 c = 1,nvir(1)
C
                  do 205 a = 1,nvir(1)
                     ia  = nrhf(1) + a
                     ic  = nrhf(1) + c
                     nla = (l-1)*nvir(1) + a
                     nai = (i-1)*nvir(1) + a
                     g2(l,i,k,c) = g2(l,i,k,c)
     &                           + rint(ia,i,k,ic)*tr1d(nla)
     &                           - rint(l,ia,k,ic)*tr1e(nai)
C
                     g6(i,l,c,k) = g6(i,l,c,k)
     &                           + rint(ia,l,k,ic)*tr1d(nai)
     &                           - rint(i,ia,k,ic)*tr1e(nla)
C
  205             continue
  204          continue
  203       continue
  202    continue
  201 continue
C
      do 301 k = 1,nrhf(1)
         do 302 c = 1,nvir(1)
            do 303 a = 1,nvir(1)
               do 304 d = 1,nvir(1)
C
                  do 305 b = 1,nvir(1)
                     ib  = nrhf(1) + b
                     ic  = nrhf(1) + c
                     ia  = nrhf(1) + a
                     id  = nrhf(1) + d
                     nkb = (k-1)*nvir(1) + b
                     g3(k,c,a,d) = g3(k,c,a,d)
     &                           + rint(ib,ic,ia,id)*tr1d(nkb)
C
                     g7(c,k,d,a) = g7(c,k,d,a)
     &                           - rint(ic,ib,id,ia)*tr1e(nkb)
C
  305             continue
C
                  do 306 j = 1,nrhf(1)
                     ia  = nrhf(1) + a
                     id  = nrhf(1) + d
                     njc = (j-1)*nvir(1) + c
                     g3(k,c,a,d) = g3(k,c,a,d)
     &                           - rint(k,j,ia,id)*tr1d(njc)
C
                     g7(c,k,d,a) = g7(c,k,d,a)
     &                           + rint(j,k,id,ia)*tr1e(njc)
C
  306             continue
C
  304          continue
  303       continue
  302    continue
  301 continue
C
      do 401 k = 1,nrhf(1)
         do 402 c = 1,nvir(1)
            do 403 l = 1,nrhf(1)
               do 404 i = 1,nrhf(1)
C
                  do 405 b = 1,nvir(1)
                     ib  = nrhf(1) + b
                     ic  = nrhf(1) + c
                     nkb = (k-1)*nvir(1) + b
                     g4(k,c,l,i) = g4(k,c,l,i)
     &                           + rint(ib,ic,l,i)*tr1d(nkb)
C
                     g8(c,k,i,l) = g8(c,k,i,l)
     &                           - rint(ic,ib,i,l)*tr1e(nkb)
C
  405             continue
C
                  do 406 j = 1,nrhf(1)
                     njc = (j-1)*nvir(1) + c
                     g4(k,c,l,i) = g4(k,c,l,i)
     &                           - rint(k,j,l,i)*tr1d(njc)
C
                     g8(c,k,i,l) = g8(c,k,i,l)
     &                           + rint(j,k,i,l)*tr1e(njc)
C
  406             continue
C
  404          continue
  403       continue
  402    continue
  401 continue
C
      do 1001 a = 1,nvir(1)
         do 1002 i = 1,nrhf(1)
            do 1003 k = 1,nrhf(1)
               do 1004 c = 1,nvir(1)
                  nck = nvir(1)*(k-1) + c
C
                  do 1005 d = 1,nvir(1)
C
                     ndi       = nvir(1)*(i-1) + d
                     nckdi     = index(nck,ndi)
                     uai1(a,i) = uai1(a,i) + g1(a,d,k,c)*t2mp(nckdi)
                     uai3(a,i) = uai3(a,i) + g3(k,c,a,d)*t2mp(nckdi)
                     uai5(a,i) = uai5(a,i) + g5(d,a,c,k)*t2mp(nckdi)
                     uai7(a,i) = uai7(a,i) + g7(c,k,d,a)*t2mp(nckdi)
C
 1005             continue
C
                  do 1006 l = 1,nrhf(1)
C
                     nal       = nvir(1)*(l-1) + a
                     nckal     = index(nck,nal)
                     uai2(a,i) = uai2(a,i) - g2(l,i,k,c)*t2mp(nckal)
                     uai4(a,i) = uai4(a,i) - g4(k,c,l,i)*t2mp(nckal)
                     uai6(a,i) = uai6(a,i) - g6(i,l,c,k)*t2mp(nckal)
                     uai8(a,i) = uai8(a,i) - g8(c,k,i,l)*t2mp(nckal)
C
 1006             continue
C
 1004          continue
 1003       continue
 1002    continue
 1001 continue
C
C
      do 2001 a = 1,nvir(1)
         ia = nrhf(1) + a
         do 2002 i = 1,nrhf(1)
            do 2003 b = 1,nvir(1)
               ib = nrhf(1) + b
               do 2004 j = 1,nrhf(1)
                  x1 = two*rint(ia,i,ib,j) - rint(ia,j,ib,i)
                  x2 = two*rint(ia,i,j,ib) - rint(ia,ib,j,i)
                  njb = (j-1)*nvir(1) + b
                  fai(a,i) = fai(a,i) + tr1e(njb)*x1 - tr1d(njb)*x2
 2004          continue
 2003       continue
 2002    continue
 2001 continue
C
      print*,'Keld NOD: fai = ',((fai(a,i),a=1,nvir(1)),i=1,nrhf(1))
C
      return
      end
