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
*====================================================================*
      subroutine cc_r12prepccsd(work,lwork)
c--------------------------------------------------------------------
c      purpose: prepare V intermediates for the CCSD(R12) model
c               using ansatz 1 that do not depend on cluster ampl.
c
c      H. Fliegl, C. Haettig summer 2004
c--------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
#include <ccr12int.h>
#include <dummy.h>
#include <ccsdinp.h>

      logical locdbg,lvajkl,lvabkl,lv,lvijkl,lexist
      parameter (locdbg = .false.)

      integer lwork,ibasx(8),isym,kend0,lwrk0
      integer klamdp,klamdh,kend1,kt1am,kend2,lwrk2
      integer isymtr,kvabkl,kend3,lwrk3,kvajkl,
     &        isymc,isymab,isymkl,lunit,idum
      integer ioptbas,iopt

#if defined (SYS_CRAY)
      real zero,one,work(*),ddummy,timrdao,timfr12,timintr12,ddot
#else
      double precision zero,one,work(*),ddummy,timrdao,timfr12,timintr12
      double precision ddot
#endif
      parameter (zero = 0.0D0, one = 1.0D0)

      call qenter('prepccsd')
      if (locdbg) then
        write(lupri,*)'entered cc_r12prepccsd'
      end if
c
      kend0 = 1
c
c     test if file fvabkl already exists, if yes exit:
      inquire(file=fvabkl,exist=lexist)
      if (lexist.and.ccrstr) then
        write(lupri,*) 'Restart: Found V_(alpha beta)^(kl) on disk'
        call qexit('prepccsd')
        return
      end if
c
c     get CMO coefficients:
c
      klamdp = kend0
      klamdh = klamdp + nlamdt
      kend1  = klamdh + nlamdt
      kt1am  = kend1
      kend2  = kt1am + nt1amx
      lwrk2  = lwork - kend2
      if(lwrk2.lt.0) then
        call quit('insufficient work space in prepccsd')
      end if
      call dzero(work(kt1am),nt1amx)
      call lammat(work(klamdp),work(klamdh),work(kt1am),
     &            work(kend2),lwrk2)

c     get V^kl_alpha_beta
      kvabkl = kend2
      kend2  = kvabkl + nvabkl(1)
      lwrk2  = lwork - kend2 
      if(lwrk2.lt.0) then
        call quit('insufficient work space in prepccsd')
      end if
      ! initialize V_(alpha beta)^(kl)
      iopt = 0
      call cc_r12mkvamkl0(work(kvabkl),nvabkl(1),iopt,ddummy,idummy,
     &                    work(kend2),lwrk2)
     
      if (locdbg) then
        write(lupri,*)'norm^2 Sak,bl =', 
     &                 ddot(nvabkl(1),work(kvabkl),1,work(kvabkl),1)
        write(lupri,*)'Sak*Sbl'
        do isymab = 1, nsym
          isymkl = isymab
          call output(work(kvabkl+ivabkl(isymab,isymkl)),1,
     &              n2bst(isymab),1,nmatkl(isymkl),n2bst(isymab),
     &              nmatkl(isymkl),1,lupri)
        end do
      end if

      if (locdbg) then
c       test if V is ok, transform all indices to occupied
        isymc = 1
        call cc_r12vtest(work(kvabkl),work(klamdp),isymc,
     &                   work(kend2),lwrk2) 
      end if
c
      isymtr = 1
      lvabkl  = .true.
      ioptbas = 2
      call cc_mofconr12(work(klamdh),1,ddummy,ddummy,ddummy,
     &                  isymtr,ddummy,ddummy,ddummy,work(kvabkl),
     &                  .false.,.false.,lvabkl,ioptbas,
     &                  timrdao,timfr12,timintr12,
     &                  idummy,idummy,idummy,idummy,idummy,idummy,
     &                  work(kend2),lwrk2)

      if (locdbg) then
        write(lupri,*)'norm^2 Vabkl: ', 
     &                 ddot(nvabkl(1),work(kvabkl),1,work(kvabkl),1)
        write(lupri,*)'V^ab_kl in prepccsdr12'
        do isymab = 1, nsym
          isymkl = isymab
          call output(work(kvabkl+ivabkl(isymab,isymkl)),
     &                1,n2bst(isymab),1,nmatkl(isymkl),
     &                n2bst(isymab),nmatkl(isymkl),1,lupri)
        end do
      end if

c     symmetrize V^ kl_alpha_beta and write it on file
      call cc_r12_symv(work(kvabkl),work(kend2),lwrk2)
      lunit = -1
      call gpopen(lunit,fvabkl,'unknown',' ','unformatted',idum,.false.)
      write(lunit)(work(kvabkl-1+i),i=1,nvabkl(1))
      call gpclose(lunit,'KEEP')

      if (locdbg) then
        write(lupri,*)'norm^2 Vabkl sym',
     &                 ddot(nvabkl(1),work(kvabkl),1,work(kvabkl),1)
        write(lupri,*)'symmetrized V^ab_kl in prepccsdr12'
        do isymab = 1, nsym
          isymkl = isymab
          call output(work(kvabkl+ivabkl(isymab,isymkl)),
     &                1,n2bst(isymab),1,nmatkl(isymkl),
     &                n2bst(isymab),nmatkl(isymkl),1,lupri)
        end do
c       test if V is ok, transform all indices to occupied
        isymc = 1
        call cc_r12vtest(work(kvabkl),work(klamdp),isymc,
     &                   work(kend2),lwrk2) 
      end if

c     calculate V^kl_alphaj and write it on file
      isymc = 1
      kvajkl = kend2
      kend3  = kvajkl + nvajkl(isymc)
      lwrk3  = lwork - kend3
      if(lwrk3.lt.0) then
        call quit('insufficient work space in prepccsd')
      end if
      lv = .false.
      lvajkl = .true.
      lvijkl = .false.
      call cc_r12mkvtf(work(kvabkl),work(kvajkl),dummy,
     &                 work(klamdp),isymc,lv,lvijkl,lvajkl,fvajkl,
     &                 work(kend3),lwrk3)

c     calculate V^ab_kl = \sum_albe L^h_ala L^h_beb V^albe_kl and
c     write it on file
      isymc = 1
      iopt  = 1
      call cc_r12mkvirt(work(kvabkl),work(klamdh),isymc,
     &                  work(klamdh),isymc,fvcdkl,iopt,
     &                  work(kend2),lwrk2)

      call qexit('prepccsd')
      return
      end 
*====================================================================*
      subroutine cc_r12mkvtf(vabkl,vajkl,vijkl,cmo,isymc,
     &                       lv,lvijkl,lvajkl,filvajkl,work,lwork)
c--------------------------------------------------------------------
c     purpose: transform V^alpha_beta _kl to V^ij_kl to test V with
c              MP2-R12 
c
c     lv       read V^(alpha beta)_kl from file
c     lvijkl   return vijkl (result is ADDED to vijkl)
c     lvajkl   write vajkl to file filvajkl   
c
c     Note: when lvijkl and lvajkl are both .FALSE. only vajkl
c           will be returned (in memory) 
c            
c
c     H. Fliegl, C. Haettig, summer 2004 
c     modified C. Neiss, autumn 2005             
c--------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
#include <ccr12int.h>

      logical lres,locdbg,lv,lvajkl,lvijkl
      parameter (locdbg = .false.)   
      integer isymc,lwork,isymkl,isymab,isymk,
     &        isyml,isyma,isymb,isymaj,isymj,
     &        koffc,ntota,ntotb,idxkl,
     &        kvabkl,kvajkl
      integer luvajkl,lunit,idummy
      character*(*) filvajkl
#if defined (SYS_CRAY)
      real vabkl(*),vajkl(*),vijkl(*),cmo(*),work(*),one,zero,ddot
#else
      double precision vabkl(*),vajkl(*),vijkl(*),cmo(*),work(*),one,
     &                 zero,ddot
#endif
      parameter (one = 1.0d0, zero = 0.0d0)

      call qenter('r12mkvtf')
      
      if (lv) then
c       read in Vabkl
        lunit = -1
        call gpopen(lunit,fvabkl,'unknown',' ','unformatted',
     &            idummy,.false.)     
        read (lunit)(vabkl(i),i=1,nvabkl(1))
        call gpclose(lunit,'KEEP')
      end if

      call dzero(vajkl,nvajkl(isymc))
      do isymkl = 1, nsym
        isymab = isymkl
        do isyma = 1, nsym
          isymb = muld2h(isymab,isyma)
          isymj = muld2h(isymc,isymb)
          isymaj = muld2h(isyma,isymj)        
          do isymk = 1, nsym
            isyml = muld2h(isymkl,isymk)
            do k = 1, nrhfb(isymk)
              do l = 1, nrhfb(isyml)
                idxkl = imatkl(isymk,isyml)+nrhfb(isymk)*(l-1)+k
                kvabkl = ivabkl(isymab,isymkl)+n2bst(isymab)*(idxkl-1)+
     &                   iaodis(isyma,isymb)+1
                kvajkl = ivajkl(isymaj,isymkl)+nt1ao(isymaj)*(idxkl-1)+
     &                   it1ao(isyma,isymj)+1
           
                koffc  = iglmrh(isymb,isymj)+1
                ntotb = max(mbas1(isymb),1)
                ntota = max(mbas1(isyma),1)
                call dgemm('N','N',mbas1(isyma),nrhf(isymj),
     &                     mbas1(isymb),one,vabkl(kvabkl),ntota,
     &                     cmo(koffc),ntotb,one,vajkl(kvajkl),ntota)    
              end do
            end do
          end do
        end do
      end do

      if (locdbg) then
        write(lupri,*)'norm^2 vajkl', 
     &                 ddot(nvajkl(isymc),vajkl,1,vajkl,1)
c       write(lupri,*)'vajkl'
c       do isymkl = 1, nsym
c         isymaj = muld2h(isymc,isymkl)
c         call output(vajkl(1+ivajkl(isymaj,isymkl)),1,
c    &                nt1ao(isymaj),1,nmatkl(isymkl),
c    &                nt1ao(isymaj),nmatkl(isymkl),1,lupri)
c       end do
      end if
c
      if (lvajkl) then       
c       -----------------------------
c        write V(alpha j,kl) on file
c       -----------------------------
        luvajkl = -1
        call gpopen(luvajkl,filvajkl,'unknown',' ','unformatted',
     &       idummy,.false.)
        rewind(luvajkl)
        write(luvajkl) (vajkl(i), i = 1,nvajkl(isymc))
        call gpclose(luvajkl,'KEEP')
      end if
c
      if (lvijkl) then
        lres = .true.
        call cc_r12mkvijkl(vajkl,isymc,cmo,isymc,work,lwork,
     &                  lres,one,vijkl)
     
        if (locdbg) then
          write(lupri,*)'in mkvtf: norm^2 vijkl:', 
     &     ddot(ntr12sq(1),vijkl,1,vijkl,1)
        end if
      end if
c
      call qexit('r12mkvtf')
      end
*====================================================================*
      subroutine cc_r12mkvabkl(vabkl,xint,idel,isymd,isydis,ibastyp,
     &                         ibasx,work,lwork)
c--------------------------------------------------------------------
c     purpose: make V^kl_alpha_beta for CCSD(12) model with ansatz 1
c
c     H. Fliegl, C. Haettig, summer 2004
c
c     C. Neiss, 05.12.2005: adapted for CABS
c--------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <maxorb.h>
#include <ccorb.h>
#include <ccisao.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <r12int.h>
#include <ccr12int.h>

      logical locdbg,lauxd
      parameter(locdbg =.false.)
     
      integer nr1orb(8),nr1bas(8),nr1xbas(8),nr2bas,n2bst1(8)
      integer ir1xbas(8,8),ir1orb(8,8),ir1bas(8,8),ir2bas(8,8),
     &        iaodis1(8,8),ir2xbas(8,8),irgkl(8,8),nrgkl(8),
     &        nalphaj(8),ialphaj(8,8)
      integer nr1xorb(8),ir1xorb(8,8),nrxgkl(8),irxgkl(8,8)
      integer idel,isymd,isydis,ibastyp,lwork,lwrk1,kend1,krgkl,
     &        ibasx(8),isym,isymg,isymab,isymkl,koffr,koffg,koffv,
     &        ntotg,ntotab,isym1,isym2,icount1,ngab(8),igab(8,8)
      integer kscr1,kscr2,isymb,isymga,isymgab,koff1,
     &        idxab,idxga,idxgab,isyma,kend2,lwrk2,isymgkl

#if defined (SYS_CRAY)
      real work(*),xint(*),vabkl(*),one,two,factor,ddot
#else
      double precision work(*),xint(*),vabkl(*),one,two,factor,
     &                 ddot
#endif
      parameter (one = 1.0d0, two = 2.0d0)

      call qenter('mkvabkl')

      call cc_r12offset(nr1orb,nr1xorb,nr1bas,nr1xbas,nr2bas,
     &     nrgkl,nrxgkl,n2bst1,ir1orb,ir1xorb,ir1bas,ir1xbas,
     &     ir2bas,ir2xbas,irgkl,irxgkl,iaodis1,nalphaj,ialphaj)
c
      do isym = 1, nsym
        ngab(isym) = 0
        icount1    = 0
        do isym2 = 1, nsym
          isym1 = muld2h(isym,isym2)
          ngab(isym) = ngab(isym) + mbas1(isym1)*n2bst(isym2)
          igab(isym1,isym2) = icount1
          icount1 = icount1 + mbas1(isym1)*n2bst(isym2)
        end do
      end do
c
      krgkl = 1
      kend1 = krgkl + nrgkl(isymd) 
      lwrk1 = lwork - kend1
      if (lwrk1.lt.0) then
        call quit('insufficient work space in mkvabkl')
      end if

      if (locdbg) then
        write(lupri,*) 'ibastyp,isymd,idel:',ibastyp,isymd,idel
        write(lupri,*) 'ibasx:',ibasx
      end if

c     get r12 integrals
      if (ibastyp.eq.1) then
        lauxd = .false.
        call cc_r12getrint(work(krgkl),idel,isymd,nr1bas,ir1bas,
     &     nr2bas,ir2bas,nrgkl,irgkl,ir1xbas,ir2xbas,
     &     nrhfb,nmatkl,imatkl,ibasx,lauxd,.false.,
     &     fnback,work(kend1),lwrk1)
        factor = one
        if (r12cbs) factor = -one
      else
        lauxd = .true.
        call cc_r12getrint(work(krgkl),idel,isymd,nr1bas,ir1bas,
     &     nr2bas,ir2bas,nrgkl,irgkl,ir1xbas,ir2xbas,
     &     nrhfb,nmatkl,imatkl,ibasx,lauxd,.false.,
     &     fnback,work(kend1),lwrk1)
        factor = - two
      end if
c
      isymgab = isydis
c
      kscr2  = kend1
      kend1  = kscr2 + ngab(isymgab)
      lwrk1  = lwork - kend1
      if (lwrk1.lt.0) then
        call quit('insufficient work space in mkvabkl')
      end if
c
      do isymb = 1, nsym
        isymga = muld2h(isydis,isymb)
        
c       dynamic work space allocation
        kscr1  = kend1
        kend2  = kscr1 + n2bst(isymga)
        lwrk2  = lwork - kend2
        if (lwrk2.lt.0) then
          call quit('insufficient work space in mkvabkl')
        end if
c       
        do b = 1, mbas1(isymb)
c         pack triangular matrix to quadratic matrix
          koff1 = idsaog(isymb,isydis) + 1 + nnbst(isymga)*(b-1) 
          call ccsd_symsq(xint(koff1),isymga,work(kscr1))
c       
          do isyma = 1, nsym 
           isymab = muld2h(isymb,isyma)
           isymg  = muld2h(isymga,isyma)
           do a = 1, mbas1(isyma) 
            idxab = iaodis(isyma,isymb)+mbas1(isyma)*(b-1)+a
            do g = 1, mbas1(isymg)
              idxga = iaodis(isymg,isyma)+mbas1(isymg)*(a-1)+g
              idxgab = igab(isymg,isymab)+mbas1(isymg)*(idxab-1)+g 
              work(kscr2-1+idxgab) = work(kscr1-1+idxga)
            end do
           end do
          end do
c       
        end do !b
      end do ! isymb

      if (locdbg) then
        write(lupri,*)'norm^2 xint:',
     &    ddot(ndisao(isydis),xint,1,xint,1)
        write(lupri,*)'norm^2 work(kscr2): ', 
     &    ddot(ngab(isymgab),work(kscr2),1,work(kscr2),1)
        write(lupri,*)'norm^2 work(krgkl)',
     &    ddot(nrgkl(isymd),work(krgkl),1,work(krgkl),1)
      end if

      do isymg = 1, nsym
        isymab = muld2h(isymgab,isymg)
        isymkl = isymab
        isymgkl = muld2h(isymg,isymkl)

        koffr  = krgkl + irgkl(isymg,isymkl)
        koffg  = kscr2 + igab(isymg,isymab) 
        koffv  = ivabkl(isymab,isymkl)+1
        
        ntotab = max(n2bst(isymab),1)
        ntotg  = max(mbas1(isymg),1)
        
        if (locdbg) then
          write(lupri,*)'norm^2 work(koffg)', 
     &    ddot(mbas1(isymg)*n2bst(isymab),work(koffg),1,work(koffg),1)
          write(lupri,*)'norm^2 work(koffr)',
     &    ddot(mbas1(isymg)*nmatkl(isymkl),work(koffr),1,work(koffr),1)
        end if

        call dgemm('T','N',n2bst(isymab),nmatkl(isymkl),mbas1(isymg),
     &             factor,work(koffg),ntotg,work(koffr),ntotg,one,
     &             vabkl(koffv),ntotab)
      end do

      call qexit('mkvabkl')
      end
*====================================================================*
      subroutine cc_r12vtest(vabkl,cmo,isymc,work,lwork)
c--------------------------------------------------------------------
c     purpose: test if V^ab_kl is ok
c
c     H. Fliegl, C. Haettig, summer 2004
c--------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
#include <ccr12int.h>

      logical lv,lvajkl,lvijkl
      integer kvijkl,kvajkl,isymc,
     &        kend1,lwrk1,lwork

#if defined (SYS_CRAY)
      real vabkl(*),work(*),cmo(*),ddot
#else
      double precision vabkl(*),work(*),cmo(*),ddot
#endif
      call qenter('r12vtest')

      kvijkl = 1
      kvajkl = kvijkl + nrhftb*nrhftb*nrhftb*nrhftb 
      kend1  = kvajkl + nvajkl(isymc)
      lwrk1  = lwork - kend1
      if (lwrk1.lt.0) then
        call quit('insufficient work space in r12vtest')
      end if
      call dzero(work(kvijkl),nrhftb*nrhftb*nrhftb*nrhftb)
      lv = .false.
      lvajkl = .false.
      lvijkl = .true.
      call cc_r12mkvtf(vabkl,work(kvajkl),work(kvijkl),
     &                 cmo,isymc,lv,lvijkl,lvajkl,fvajkl,
     &                 work(kend1),lwrk1)
      write(lupri,*)'Norm^2 Vijkl in r12vtest: ',
     & ddot(nrhftb*nrhftb*nrhftb*nrhftb,work(kvijkl),1,work(kvijkl),1)
      call output(work(kvijkl),1,nrhftb*nrhftb,1,nrhftb*nrhftb,
     &            nrhftb*nrhftb,nrhftb*nrhftb,1,lupri)
      write(lupri,*)'Vijkl in r12vtest in triangular format:'
      call ccr12pck2(work(kend1),isymc,.FALSE.,work(kvijkl),'T',1)
      call cc_prpr12(work(kend1),isymc,1,.FALSE.)
        
      call qexit('r12vtest')
      end  
*====================================================================*
      subroutine cc_r12_symv(vabkl,work,lwork)
c--------------------------------------------------------------------
c     purpose: symmetrize V^kl_ab 
c              P^ab_kl * V^kl_ab = 1/2 * (V^kl_ab + V^lk_ba)
c
c     H. Fliegl, C. Haettig, summer 2004
c--------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>

      integer lwork,kvsym,kend1,lwrk1,isyma,isymb,
     &        isymab,isymkl,isyml,isymk,idxkl,idxlk,idxab,idxba,idxabkl,
     &        idxbalk
#if defined (SYS_CRAY)
      real vabkl(*),work(*),half
#else
      double precision vabkl(*),work(*),half
#endif
      parameter (half = 0.5d0)
      call qenter('r12symv')

      kvsym = 1
      kend1 = kvsym + nvabkl(1)
      lwrk1 = lwork - kend1
      if (lwrk1.lt.0) then
        call quit('insufficient work space in r12symv')
      end if

      do isyma = 1, nsym
        do isymb = 1, nsym
          isymab = muld2h(isyma,isymb)
          isymkl = isymab
          do isymk = 1, nsym
            isyml = muld2h(isymkl,isymk)
            do k = 1, nrhfb(isymk)
              do l = 1, nrhfb(isyml)
                idxkl = imatkl(isymk,isyml)+nrhfb(isymk)*(l-1)+k
                idxlk = imatkl(isyml,isymk)+nrhfb(isyml)*(k-1)+l 
                do a = 1, mbas1(isyma)
                  do b = 1, mbas1(isymb)
                    idxab = iaodis(isyma,isymb)+mbas1(isyma)*(b-1)+a
                    idxba = iaodis(isymb,isyma)+mbas1(isymb)*(a-1)+b
                    idxabkl = ivabkl(isymab,isymkl)+
     &                        n2bst(isymab)*(idxkl-1)+idxab
                    idxbalk = ivabkl(isymab,isymkl)+
     &                        n2bst(isymab)*(idxlk-1)+idxba
                    work(kvsym-1+idxabkl) = 
     &                        half*(vabkl(idxabkl)+vabkl(idxbalk)) 
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
      
      call dcopy(nvabkl(1),work(kvsym),1,vabkl,1)

      call qexit('r12symv')
      return
      end
*====================================================================*
      SUBROUTINE CCRHS_BP(OMEGA2,ISYMTR,IOPTOM,IAMP,FC12AM,LUFC12,IFILE,
     &                    LISTR,IDLSTR,BASSCL2,WORK,LWORK)
c--------------------------------------------------------------------
c     purpose: calculate \sum_mn C^ij_mn * V^mn_ab
c
c     IOPTOM = 0: dimension of OMEGA2 = 2*NT2ORT(ISYMTR)
c            = 1: dimension of OMEGA2 = NT2AO(ISYMTR)
c
c     H. Fliegl, C. Haettig, summer 2004
c     modified C. Neiss, autumn 2005
c--------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <dummy.h>
#include <r12int.h>
#include <ccr12int.h>
#include <ccsdinp.h>


      logical locdbg
      parameter (locdbg = .false.)
      integer lwork,kend1,lwrk1,
     &        isymv,idum,kvabkl,lunit,krpck,isymc,ntamp
      integer isymab,iopt,iamp,lufc12,ifile,isymtr,isymij,idlstr,ioptom
      character*10 model
      character*8 fc12am
      character*3 listr
#if defined (SYS_CRAY)
      real work(*),one,two,factor,zero,omega2(*),
     &     four,ddot,x,basscl2
#else
      double precision work(*),one,two,
     &                 factor,zero,omega2(*),four,ddot,x,basscl2
#endif
      parameter (one = 1.0d0, two = 2.0d0, zero = 0.0d0, four = 4.0d0)

      call qenter('ccrhs_bp')

      if (locdbg.and.(ioptom.eq.0)) then
       write(lupri,*)'Omega 2 entering ccrhs_bp'
       do isymab = 1, nsym
        isymij = muld2h(isymab,isymtr)
        write(lupri,*) 'isymab, isymij: ', isymab,isymij
        write(lupri,*) 'Omega+ ='
        call output(omega2(it2ort(isymab,isymij)+1),1,nnbst(isymab),1,
     &              nmijp(isymij),nnbst(isymab),nmijp(isymij),1,lupri)
        write(lupri,*)'Omega- ='
        call output(omega2(it2ort(isymab,isymij)+1+nt2ort(isymtr)),1,
     &              nnbst(isymab),1,nmijp(isymij),nnbst(isymab),
     &              nmijp(isymij),1,lupri)

       end do
      end if

C     call cc_r12offset(nr1orb,nr1xorb,nr1bas,nr1xbas,nr2bas,
C    & nrgkl,nrxgkl,n2bst1,ir1orb,ir1xorb,ir1bas,ir1xbas,ir2bas,ir2xbas,
C    & irgkl,irxgkl,iaodis1,nalphaj,ialphaj)
c
      krpck  = 1
      kvabkl = krpck + ntr12am(isymtr)
      kend1  = kvabkl + nvabkl(1)
      lwrk1  = lwork - kend1
      if (lwrk1 .lt.0) then
         call quit('Insufficient work space in ccrhs_bp')
      end if

c     read V 
      lunit = -1
      call gpopen(lunit,fvabkl,'unknown',' ','unformatted',idum,.false.)
      read(lunit)(work(kvabkl-1+i),i=1,nvabkl(1))
      call gpclose(lunit,'KEEP')
        
c     read r12 amplitudes 
      if (iamp.eq.0) then 
        iopt = 32
        call cc_rdrsp('R0 ',0,1,iopt,model,dummy,work(krpck))
      else if (iamp.eq.1) then
        call cc_rvec(lufc12,fc12am,ntr12am(isymtr),ntr12am(isymtr),
     &               ifile,work(krpck))
        call cclr_diasclr12(work(krpck),basscl2,isymtr)
      else if (iamp.eq.2) then
        iopt = 32
        call cc_rdrsp(listr,idlstr,isymtr,iopt,model,dummy,work(krpck))
        call cclr_diasclr12(work(krpck),basscl2,isymtr)
      else
        call quit('Unknown value of IAMP in CCRHS_BP')
      end if

c     calculate c*V
      isymv = 1
      isymc = isymtr
      factor = one
      call cc_r12mkcv(omega2,ioptom,work(kvabkl),work(krpck),
     &                isymv,isymc,factor,work(kend1),lwrk1)

      if (locdbg) then
        write(lupri,*)'in ccrhs_bp:'
        write(lupri,*)'norm^2 work(kvabkl):', 
     &   ddot(nvabkl(1),work(kvabkl),1,work(kvabkl),1)
        write(lupri,*)'norm^2 work(krpck):', 
     &   ddot(ntr12am(isymtr),work(krpck),1,work(krpck),1)
        if (ioptom.eq.0) then
          write(lupri,*)'norm^2 omega2', 
     &      ddot(2*nt2ort(isymtr),omega2,1,omega2,1)
        else if (ioptom.eq.1) then
          write(lupri,*)'norm^2 omega2',
     &      ddot(nt2ao(isymtr),omega2,1,omega2,1)
        end if
      end if

      if (locdbg.and.(ioptom.eq.0)) then
       write(lupri,*)'Omega 2 leaving ccrhs_bp'
       do isymab = 1, nsym
        isymij = muld2h(isymab,isymtr)
        write(lupri,*) 'isymab, isymij: ', isymab,isymij
        write(lupri,*) 'Omega+ ='
        call output(omega2(it2ort(isymab,isymij)+1),1,nnbst(isymab),1,
     &              nmijp(isymij),nnbst(isymab),nmijp(isymij),1,lupri)
        write(lupri,*)'Omega- ='
        call output(omega2(it2ort(isymab,isymij)+1+nt2ort(isymtr)),1,
     &              nnbst(isymab),1,nmijp(isymij),nnbst(isymab),
     &              nmijp(isymij),1,lupri)
       end do
      end if
 
      call qexit('ccrhs_bp')
      end
*====================================================================*
      subroutine cc_r12mkcv(omega2,ioptom,vabkl,c2am,isymv,
     &                      isymc,factor,work,lwork)
c--------------------------------------------------------------------
c     purpose: calculate \sum_mn C^ij_mn * V^mn_ab
c
c     H. Fliegl, C. Haettig, summer 2004
c
c     modified C. Neiss, autumn 2005
c
c     IOPTOM = 0: sort Omega2 in singlet/triplet format 
c                 (dim: 2*NT2ORT(ISYMOM)) 
c            = 1: sort Omega2 in triangular matrix format
c                 (dim: NT2AO(ISYMOM))
c--------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <maxorb.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <symsq.h>
#include <r12int.h>
#include <ccr12int.h>
 
      logical locdbg
      parameter (locdbg = .false.)

      integer lwork,isymv,isymc,isymom,isymb,
     &        isyma,isymab,isymkl,isymij,isymj,isymi,isymai,isymbj,
     &        kvabkl,kc2am,kom,kend1,lwrk1,ntoti,ntota,idxbj,idxai,
     &        koffai,maxai,naibj,index,kom2,kend2,lwrk2,idxij,idxab,
     &        idxabijp,idxabijm,idxbi,idxaj,idxaibj,idxbiaj,isymbi,
     &        isymaj,nab,kend3,lwrk3,maxj,maxb,ioptom 
      integer isymlj,idxljlj,idxlj,isymk,isyml,idxki,idxkilj
#if defined (SYS_CRAY)
      real work(*),vabkl(*),one,two,c2am(*),factor,zero
      real half,omega2(*),ddot,x1,xv,xc,dummy
#else
      double precision work(*),vabkl(*),one,two,c2am(*),
     &                 factor,zero,half,omega2(*),ddot,
     &                 x1,xv,xc,dummy
#endif
      parameter (one = 1.0d0, two = 2.0d0, zero = 0.0d0, half = 0.5d0)

      index(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j

      call qenter('r12mkcv')
      if (.false.) then
c       small test if gamma term in ccrhs is exact V^ij_kl matrix
c       Caution!!!! d_ki, d_lj have to be put to one !!!!
        call dzero(c2am,ntr12am(1))
        do isymk = 1, nsym
          isymi = isymk
          do isyml = 1, nsym
            isymj = isyml
            do i = 1, nrhf(isymi)              
              k = i
              idxki = imatki(isymk,isymi)+nrhfb(isymk)*(i-1)+k
              do j = 1, nrhf(isymj)
                l = j
                idxlj = imatki(isyml,isymj)+nrhfb(isyml)*(j-1)+l
                idxkilj = itr12am(1,1)+index(idxki,idxlj) 
                c2am(idxkilj) = one
              end do
            end do
          end do
        end do
        write(lupri,*)'c2am:'
        call outpak(c2am,nmatij(1),1,lupri) 
        call dzero(omega2,2*nt2ort(1))
      end if

      isymom = muld2h(isymv,isymc)

      kend1 = 1

      if (ioptom.eq.0) then
        kom2 = 1
        kend1 = kom2 + nt2ao(isymom)
        lwrk1 = lwork - kend1
        if (lwrk1.lt.0) then
          call quit('insufficient work space in r12mkcv')
        end if
        call dzero(work(kom2),nt2ao(isymom))
      else if (ioptom.ne.1) then
        call quit('Unknown value of "IOPTOM" in CC_R12MKCV')
      end if

      if (locdbg) then
        x1 = 0.0d0
        xv = 0.0d0
        xc = 0.0d0
        write(lupri,*)'norm^2(c2am):',ddot(ntr12am(isymc),c2am,1,c2am,1)
        write(lupri,*)'isymom,nt2ao:',isymom,nt2ao(isymom)
        if (ioptom.eq.0) then
          write(lupri,*) 'norm^2 work(kom2)',
     &                   ddot(nt2ao(isymom),work(kom2),1,work(kom2),1)
        else if (ioptom.eq.1) then
          write(lupri,*) 'norm^2 omega2',
     &                   ddot(nt2ao(isymom),omega2,1,omega2,1)
        end if
      end if

      do isymb = 1, nsym
        do isyma = 1, nsym
          isymab = muld2h(isyma,isymb)
          isymkl = muld2h(isymv,isymab)
          isymij = muld2h(isymkl,isymc)

c         dynamic allocation of work space            
          kvabkl = kend1
          kend2  = kvabkl + mbas1(isyma)*nmatkl(isymkl)
          lwrk2  = lwork - kend2
          if (lwrk2.lt.0) then
            call quit('insufficient work space in r12mkcv')
          end if

c
          do b = 1, mbas1(isymb)
            call cc_r12sortv(work(kvabkl),vabkl,b,isymb,isymv,isyma)
               if (locdbg) then
                 xv = xv + ddot(mbas1(isyma)*nmatkl(isymkl),
     &                          work(kvabkl),1,work(kvabkl),1)
               end if
            do isymj = 1, nsym
              isymi  = muld2h(isymij,isymj)
              isymai = muld2h(isyma,isymi)
              isymbj = muld2h(isymb,isymj)
              kc2am  = kend2
              kom    = kc2am + nrhf(isymi)*nmatkl(isymkl)
              kend3  = kom + mbas1(isyma)*nrhf(isymi)
              lwrk3  = lwork - kend3
              if (lwrk3.lt.0) then
                call quit('insufficient work space in r12mkcv')
              end if
              do j = 1, nrhf(isymj)
                call cc_r12sortc(work(kc2am),c2am,j,isymc,isymj,isymi)
                if (locdbg) then
                  xc = xc + ddot(nrhf(isymi)*nmatkl(isymkl),
     &                           work(kc2am),1,work(kc2am),1)
                end if
c
                ntoti = max(nrhf(isymi),1)
                ntota = max(mbas1(isyma),1)
                call dgemm('N','T',mbas1(isyma),nrhf(isymi),
     &                     nmatkl(isymkl),one,work(kvabkl),ntota,
     &                     work(kc2am),ntoti,zero,work(kom),ntota)
                if (locdbg) then
                  x1 = x1 + ddot(mbas1(isyma)*nrhf(isymi),
     &                           work(kom),1,work(kom),1)
                end if

c               update Omega_alpha_i,beta_j which is stored as a triangular
c               matrix with bj.ge.ai 
                idxbj = it1ao(isymb,isymj)+mbas1(isymb)*(j-1)+b
                koffai = it1ao(isyma,isymi)
                if (isymai.eq.isymbj) then
                  maxai  = min(koffai+mbas1(isyma)*nrhf(isymi),idxbj)
                  do idxai = koffai+1, maxai
                    naibj = it2ao(isymai,isymbj)+index(idxai,idxbj)
                    if (ioptom.eq.0) then    
                      work(kom2-1+naibj) = work(kom2-1+naibj) +
     &                                  factor*work(kom-1+idxai-koffai)
                    else if (ioptom.eq.1) then
                      omega2(naibj) = omega2(naibj) + 
     &                                  factor*work(kom-1+idxai-koffai)
                    end if
                  end do
                else if (isymai.lt.isymbj) then
                  maxai = koffai+mbas1(isyma)*nrhf(isymi)
                  do idxai = koffai+1, maxai
                    naibj = it2ao(isymai,isymbj)+
     &                      nt1ao(isymai)*(idxbj-1) + idxai
                    if (ioptom.eq.0) then
                      work(kom2-1+naibj) = work(kom2-1+naibj) +
     &                                  factor*work(kom-1+idxai-koffai)
                    else if (ioptom.eq.1) then
                      omega2(naibj) = omega2(naibj) +
     &                                  factor*work(kom-1+idxai-koffai)
                    end if
                  end do
                end if 
              end do
            end do
            
          end do
        end do
      end do


      if (locdbg) then
        write(lupri,*)'xv =', xv
        write(lupri,*)'xc =', xc
        write(lupri,*)'x1 =', x1
        write(lupri,*)'in mkcv:'
        if (ioptom.eq.0) then
          write(lupri,*)'norm^2 work(kom2):', 
     &                  ddot(nt2ao(isymom),work(kom2),1,work(kom2),1) 
          call cc_prpao(dummy,work(kom2),isymom,0,1)
          write(lupri,*)'norm^2 omega2: ', 
     &                  ddot(2*nt2ort(isymom),omega2,1,omega2,1)
        else if (ioptom.eq.1) then
          write(lupri,*)'norm^2 omega2:',
     &                  ddot(nt2ao(isymom),omega2,1,omega2,1)
          call cc_prpao(dummy,omega2,isymom,0,1)
        end if
      end if
c
      if (ioptom.eq.0) then
c     sort Omega as Omega_ab,_ij triplet and singlet ..... 
c     note Omega is a triangular matrix with a.ge.b and i.ge.j
      do isymi = 1, nsym
        do isymj = 1, nsym
          if (isymj.le.isymi) then
          isymij = muld2h(isymi,isymj)
          isymab = muld2h(isymij,isymom)
          do isyma = 1, nsym
           isymb  = muld2h(isymab,isyma)
           if (isymb.le.isyma) then 

            isymai = muld2h(isyma,isymi)
            isymbj = muld2h(isymb,isymj) 
            isymbi = muld2h(isymb,isymi)
            isymaj = muld2h(isyma,isymj)
C           if (isymom.ne.(muld2h(isymai,isymbj))) 
C    *            call quit('Symmetry error in cc_r12mkcv')

            do i = 1, nrhf(isymi)

              maxj = nrhf(isymj)
              if (isymj.eq.isymi) maxj = i
              do j = 1, maxj
c               
                if (isymi.eq.isymj) then
                  idxij = imijp(isymi,isymj)+index(i,j)
                else if (isymi.lt.isymj) then
                  idxij = imijp(isymi,isymj)+nrhf(isymi)*(j-1)+i
                else
                  idxij = imijp(isymj,isymi)+nrhf(isymj)*(i-1)+j
                end if
 
                if (locdbg) then
                  write(lupri,'(a,11i5)') 
     &            'isymai,isymbj,ab,a,b,ij,i,j,i,j,idxij:',
     &             isymai,isymbj,isymab,isyma,isymb,
     &             isymij,isymi,isymj,i,j,idxij
                end if
c
                do a = 1, mbas1(isyma)
                  maxb = mbas1(isymb)
                  if (isymb.eq.isyma) maxb = a
                  do b = 1, maxb
c
                    if (isyma.eq.isymb) then
                      idxab = iaodpk(isyma,isymb)+index(a,b) 
                    else if (isyma.lt.isymb) then
                      idxab = iaodpk(isyma,isymb)+mbas1(isyma)*(b-1)+a
                    else
                      idxab = iaodpk(isymb,isyma)+mbas1(isymb)*(a-1)+b
                    end if
c  
                    idxabijp = it2ort(isymab,isymij)+
     &                         nnbst(isymab)*(idxij-1)+idxab
                    idxabijm = nt2ort(isymom)+it2ort(isymab,isymij)+
     &                         nnbst(isymab)*(idxij-1)+idxab 
c                   
                    idxai = it1ao(isyma,isymi)+mbas1(isyma)*(i-1)+a
                    idxbj = it1ao(isymb,isymj)+mbas1(isymb)*(j-1)+b
                    idxbi = it1ao(isymb,isymi)+mbas1(isymb)*(i-1)+b
                    idxaj = it1ao(isyma,isymj)+mbas1(isyma)*(j-1)+a
c
                    if (isymai.eq.isymbj) then
                      idxaibj = it2ao(isymai,isymbj)+index(idxai,idxbj)
                      idxbiaj = it2ao(isymbi,isymaj)+index(idxbi,idxaj)

                      omega2(idxabijp) = omega2(idxabijp) +
     &                 (work(kom2-1+idxaibj)+work(kom2-1+idxbiaj))
                      omega2(idxabijm) = omega2(idxabijm)+
     &                 (work(kom2-1+idxaibj)-work(kom2-1+idxbiaj))
       
                    else 

                      if (isymbj.lt.isymai) then
                        idxaibj = it2ao(isymbj,isymai)+
     &                            nt1ao(isymbj)*(idxai-1)+idxbj
                      else
                        idxaibj = it2ao(isymai,isymbj)+
     &                            nt1ao(isymai)*(idxbj-1)+idxai
                      end if

                      if (isymbi.lt.isymaj) then
                        idxbiaj = it2ao(isymbi,isymaj)+
     &                            nt1ao(isymbi)*(idxaj-1)+idxbi
                      else
                        idxbiaj = it2ao(isymaj,isymbi)+
     &                            nt1ao(isymaj)*(idxbi-1)+idxaj
                      end if

                      omega2(idxabijp) = omega2(idxabijp)+
     &                  (work(kom2-1+idxaibj)+work(kom2-1+idxbiaj))
                      omega2(idxabijm) = omega2(idxabijm)+
     &                  (work(kom2-1+idxaibj)-work(kom2-1+idxbiaj))

                    end if
c
                    if (.false.) then
                      write(lupri,'(a,8i5,2g20.10)')
     &           'a,b,ai,bj,iaibj,ibiaj,iabijp,iabijm,O(aibj),O(biaj):',
     &                a,b,idxai,idxbj,idxaibj,idxbiaj,idxabijp,idxabijm,
     &                work(kom2-1+idxaibj),work(kom2-1+idxbiaj)
                      write(lupri,*)'OM(+)', omega2(idxabijp)
                      write(lupri,*)'OM(-)', omega2(idxabijm)
                    end if

                  end do
                end do
c
              end do
            end do
           end if
          end do
          end if
        end do
      end do
c
      end if     
 
      if (.false.) then
c     if (locdbg.and.(ioptom.eq.0)) then
        write(lupri,*)'norm^2 omega2+: ', 
     &   ddot(nt2ort(isymom),omega2,1,omega2,1)
        write(lupri,*)'norm^2 omega2-: ', 
     &   ddot(nt2ort(isymom),omega2(nt2ort(isymom)+1),1,
     &                       omega2(nt2ort(isymom)+1),1)
      end if

      call qexit('r12mkcv')
      end
*====================================================================*
      subroutine cc_r12sortv(vs,vabkl,b,isymb,isymv,isyma)
c--------------------------------------------------------------------
c     purpose: sort V as V^beta_alpha,_mn with fixed beta
c              expect V stored as a square matrix
c
c     H. Fliegl, C. Haettig, summer 2004
c--------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>

      integer isymb,isymv,isyma,index,isymamn,isymmn,isymm,isymn,
     &        isymam,isymbn, idxmn,idxamn
      integer idxab,idxabmn,isymab
#if defined (SYS_CRAY)
      real vs(*),vabkl(*)
#else
      double precision vs(*),vabkl(*)
#endif
      index(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j

      call qenter('sortv')

      isymab = muld2h(isyma,isymb)
      isymmn  = muld2h(isymv,isymab)
      do idxmn = 1, nmatkl(isymmn)
        do a = 1, mbas1(isyma)
          idxab = iaodis(isyma,isymb)+mbas1(isyma)*(b-1)+a 
          idxamn = mbas1(isyma)*(idxmn-1)+a 
          idxabmn = ivabkl(isymab,isymmn)+
     &              n2bst(isymab)*(idxmn-1)+idxab
          vs(idxamn) = vabkl(idxabmn)
        end do
      end do

      call qexit('sortv')
      end
*====================================================================*
      subroutine cc_r12mkvirt(vabkl,xlamd1,isymc1,xlamd2,isymc2,
     &                        filename,icon,work,lwork)
c--------------------------------------------------------------------
c     purpose: calculate V^ab_kl = 
c              \sum_{\alpha \beta} \Lambda1_{\alpha a} 
c              \Lambda2_{\beta b} V^{\alpha \beta}_kl
c
c     icon = 2: calculate V^ab_kl =
c               \sum_{\alpha \beta} [ \Lambda1_{\alpha a}
c               \Lambda2_{\beta b} V^{\alpha \beta}_kl +
c               Lambda2_{\alpha a} \Lambda1_{\beta b} 
c               V^{\alpha \beta}_kl]
c
c     H.Fliegl, C.Haettig, summer 2004
c
c     generalized for two transformation matrices: icon = 2
c     C. Neiss, Nov. 2005
c--------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
#include <ccr12int.h>

      logical locdbg
      parameter (locdbg = .false.)
      integer lwork,icount1
      integer kvabekl,kvcdkl,isym,idum,isym1,isym2,
     &        kend1,lwrk1,kend2,lwrk2,kscr1,kscr2,kend3,lwrk3
      integer isymkl,isymalbe,isymal,isymbe,isyma,isymab,idxkl,
     &        isymk,isyml,isymak,isymbl,idxab,idxakbl,idxak,idxbl,
     &        kvabkl,koffc,koffr,ntotal,ntota,ntotbe,isymb,isymabe,
     &        idxblak
      integer lunit,icon,index
      integer nck(8),ick(8,8),nvckdl(8),ivckdl(8,8),isymv,isymc1,isymc2
      character*(*) filename
#if defined (SYS_CRAY)
      real vabkl(*),xlamd1(*),xlamd2(*),work(*),one,zero,
     &     ddot,x2,x3,x1
#else
      double precision vabkl(*),xlamd1(*),xlamd2(*),work(*),one,zero,
     &                 ddot,x2,x3,x1
#endif
      parameter(one = 1.0d0 , zero = 0.0d0)
      index(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j

      call qenter('mkvirt')

c     calculate some offsets and dimensions
      do isym = 1, nsym
        nck(isym) = 0
        icount1   = 0
        do isym2 = 1, nsym
          isym1 = muld2h(isym,isym2)
          nck(isym) = nck(isym) + nvir(isym1)*nrhfb(isym2)
          ick(isym1,isym2) = icount1
          icount1 = icount1 + nvir(isym1)*nrhfb(isym2)
        end do
      end do
      do isym = 1, nsym
        nvckdl(isym) = 0
        icount1      = 0
        do isym2 = 1, nsym
          isym1 = muld2h(isym,isym2)
          if (isym2.gt.isym1) then
            nvckdl(isym) = nvckdl(isym) + nck(isym1)*nck(isym2)
            ivckdl(isym1,isym2) = icount1
            ivckdl(isym2,isym1) = icount1
            icount1 = icount1 + nck(isym1)*nck(isym2)
          else if (isym2.eq.isym1) then
            nvckdl(isym) = nvckdl(isym) + nck(isym1)*(nck(isym2)+1)/2 
            ivckdl(isym1,isym2) = icount1
            icount1 = icount1 + nck(isym1)*(nck(isym2)+1)/2
          end if
        end do
      end do

      isymv = muld2h(isymc1,isymc2)

      kvcdkl  = 1
      kend1   = kvcdkl + nvckdl(isymv)
      lwrk1   = lwork - kend1
      if (lwrk1.lt.0) then
        call quit('insufficient work space in mkvirt')
      end if

c     transform to virtual
      call dzero(work(kvcdkl),nvckdl(isymv))      
      if (locdbg) then
        x2 = zero
        x3 = zero
      end if

      do isymkl = 1, nsym
        isymalbe = isymkl
        do isymk = 1, nsym
          isyml = muld2h(isymkl,isymk)
  
          do k = 1, nrhfb(isymk)
            do l = 1, nrhfb(isyml)
              idxkl = imatkl(isymk,isyml)+nrhfb(isymk)*(l-1)+k
              
              do isymal = 1, nsym
                isyma = muld2h(isymc1,isymal)
                isymbe = muld2h(isymalbe,isymal)
                isymb  = muld2h(isymc2,isymbe)
                isymabe = muld2h(isyma,isymbe)
  
                kvabekl = kend1
                kscr1   = kvabekl + nvir(isyma)*nbas(isymbe)
                kend2   = kscr1 + nvir(isyma)*nvir(isymb)
                lwrk2   = lwork - kend2
                if (lwrk2.lt.0) then
                 call quit('insufficient work space in mkvirt')
                end if
  
                kvabkl = ivabkl(isymalbe,isymkl)+
     &                   n2bst(isymalbe)*(idxkl-1)+
     &                   iaodis(isymal,isymbe)+1
                koffc  = iglmvi(isymal,isyma)+1
               
                ntotal = max(nbas(isymal),1)
                ntota  = max(nvir(isyma),1) 
  
                call dgemm('T','N',nvir(isyma),nbas(isymbe),
     &                    nbas(isymal),one,xlamd1(koffc),ntotal,
     &                    vabkl(kvabkl),ntotal,zero,work(kvabekl),
     &                    ntota)
                
                koffc  = iglmvi(isymbe,isymb)+1
                ntotbe = max(nbas(isymbe),1)
                ntota  = max(nvir(isyma),1) 
              
                call dgemm('N','N',nvir(isyma),nvir(isymb),nbas(isymbe),
     &                one,work(kvabekl),ntota,xlamd2(koffc),ntotbe,
     &                zero,work(kscr1),ntota)
  
                if (locdbg) then
                  write(lupri,*) 'idxkl, norm^2(kscr1): ',
     &                            idxkl, ddot(nvir(isyma)*nvir(isymb),
     &                                    work(kscr1),1,work(kscr1),1)
                  x2 = x2 + ddot(nvir(isyma)*nvir(isymb),
     &                    work(kscr1),1,work(kscr1),1)
                end if
  
C               sort result as Vakbl and pack as triangular matrix
                isymak = muld2h(isyma,isymk)
                isymbl = muld2h(isymb,isyml)
  
                do b = 1, nvir(isymb)
                  idxbl = ick(isymb,isyml)+nvir(isymb)*(l-1)+b
                  do a = 1, nvir(isyma)
                    idxab = nvir(isyma)*(b-1)+a
                    idxak = ick(isyma,isymk)+nvir(isyma)*(k-1)+a
  
                    if (isymak.eq.isymbl) then 
                      if (idxak.le.idxbl) then
                        idxakbl = ivckdl(isymak,isymbl)+
     &                            index(idxak,idxbl)
                        work(kvcdkl-1+idxakbl) = work(kvcdkl-1+idxakbl)+
     &                                           work(kscr1-1+idxab)
                      end if
                      if ((idxbl.le.idxak).and.(icon.eq.2)) then
                        idxblak = ivckdl(isymbl,isymak)+
     &                            index(idxbl,idxak) 
                        work(kvcdkl-1+idxblak) = work(kvcdkl-1+idxblak)+
     &                                           work(kscr1-1+idxab)
                      end if
                    else if (isymak.lt.isymbl) then
                      idxakbl = ivckdl(isymak,isymbl)+
     &                          nck(isymak)*(idxbl-1)+idxak
                      work(kvcdkl-1+idxakbl) = work(kvcdkl-1+idxakbl) +
     &                                         work(kscr1-1+idxab)
                    else if ((isymbl.lt.isymak).and.(icon.eq.2)) then
                      idxblak = ivckdl(isymbl,isymak)+
     &                          nck(isymbl)*(idxak-1)+idxbl
                      work(kvcdkl-1+idxblak) = work(kvcdkl-1+idxblak) +
     &                                         work(kscr1-1+idxab)
                    end if
  
                    if (locdbg) then
                      if (idxak.eq.idxbl) then
C                       sum diagonal elements
                        x3 = x3 + ddot(1,work(kvcdkl-1+idxakbl),1,
     &                            work(kvcdkl-1+idxakbl),1)
                      end if
                    end if

                  end do
                end do

              end do
            end do
          end do
  
        end do
      end do

      if (locdbg) then
        write(lupri,*)'x2 =', x2
        write(lupri,*)'x3 =', x3
        write(lupri,*)'nvckdl(1) =', nvckdl(1)

c       write(lupri,*)'norm^2 Vabkl =', 
c    &               ddot(nvckdl(isymv),work(kvcdkl),1,work(kvcdkl),1)
        write(lupri,*)'norm^2 Vabkl =', 
     &  2.0d0*ddot(nvckdl(isymv),work(kvcdkl),1,work(kvcdkl),1)-x3
        write(lupri,*)'norm^2 Vabkl triang =', 
     &               ddot(nvckdl(isymv),work(kvcdkl),1,work(kvcdkl),1)

        if (isymv.eq.1) then
          call dreinorm(work(kvcdkl),nvckdl(1),nck(1),nsym,x1)
          write(lupri,*)'dreinorm =', x1
        end if
        do isymak = 1, nsym
          isymbl = isymak
          call outpkb(work(kvcdkl+ivckdl(isymak,isymbl)),nck(isymak),
     &              1,1,lupri)
        end do
      end if

c     write Vabkl on file
      lunit = -1
      call gpopen(lunit,filename,'unknown',' ','unformatted',idum,
     &            .false.)
      rewind(lunit)
      write(lunit)(work(kvcdkl-1+i),i=1,nvckdl(isymv))
      call gpclose(lunit,'KEEP')

      call qexit('mkvirt')
      end
c====================================================================*
      subroutine dreinorm(a,na,ix,nsym,x)
c--------------------------------------------------------------------
c     purpose: calculate from triangular input matrix the norm
c              of the corresponding quadratic matrix
c
c     W. Klopper, summer 2004
c--------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension a(*),ix(*)
      x=2d0*ddot(na,a,1,a,1)
      ij=0
      do isym=1,nsym
       do k=1,ix(isym)
        ij=ij+k
        x=x-a(ij)*a(ij)
       enddo
      enddo
      return
      end
c====================================================================*
      subroutine ccrhs_bpp(omega,t2am,isymt2am,lt2pcked,
     &                     filev,isymvint,work,lwork)
c--------------------------------------------------------------------
c     purpose: calculate B''-Term and update Omega_ijkl
c              \sum_ab V^ab_kl * t^ij_ab
c
c     H. Fliegl, C. Haettig, summer 2004
c
c     modified C. Neiss, autumn 2005:
c       - general symmetry of V
c       - reduced memory requirement
c--------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
#include <ccr12int.h>

      logical locdbg,lt2pcked
      parameter (locdbg = .false.)

      integer isymt2am,lwork,kvpk,lwrk1,kend1,iopt,lunit,
     &        isymvint,idum 
      integer icount1,isym1,isym2,isym,
     &        nck(8),ick(8,8),nvckdl(8),ivckdl(8,8)
      character*(*) filev
#if defined (SYS_CRAY)
      real omega(*),t2am(*),work(*),factor,ddot,one
#else
      double precision omega(*),t2am(*),work(*),factor,ddot,one
#endif
      parameter(one = 1.0d0)

      call qenter('ccrhs_bpp')

c     calculate some offsets
      do isym = 1, nsym
        nck(isym) = 0
        icount1   = 0
        do isym2 = 1, nsym
          isym1 = muld2h(isym,isym2)
          nck(isym) = nck(isym) + nvir(isym1)*nrhfb(isym2)
          ick(isym1,isym2) = icount1
          icount1 = icount1 + nvir(isym1)*nrhfb(isym2)
        end do
      end do
      do isym = 1, nsym
        nvckdl(isym) = 0
        icount1      = 0
        do isym2 = 1, nsym
          isym1 = muld2h(isym,isym2)
          if (isym2.gt.isym1) then
            nvckdl(isym) = nvckdl(isym) + nck(isym1)*nck(isym2)
            ivckdl(isym1,isym2) = icount1
            icount1 = icount1 + nck(isym1)*nck(isym2)
          else if (isym2.eq.isym1) then
            nvckdl(isym) = nvckdl(isym) + nck(isym1)*(nck(isym2)+1)/2
            ivckdl(isym1,isym2) = icount1
            icount1 = icount1 + nck(isym1)*(nck(isym2)+1)/2
          end if
        end do
      end do

      kvpk  = 1 
      kend1 = kvpk  + max(nvckdl(isymvint),nt2am(isymt2am))
      lwrk1 = lwork - kend1
      if (lwrk1.lt.0) then
        call quit('insufficient work space in ccrhs_bpp')
      end if

      if (locdbg) then
        write(lupri,*)'Norm^2 Omega^ij_kl on entry in ccrhs_bpp',
     &                 ddot(ntr12sq(muld2h(isymvint,isymt2am)),
     &                      omega,1,omega,1)
        write(lupri,*) 'Omega on entry in ccrhs_bpp:'
        call cc_prsqr12(omega,muld2h(isymvint,isymt2am),'T',1,.FALSE.)
      end if

      if (.not.lt2pcked) then
c     pack T2 amplitudes as a triangular matrix and put on T2AM
        iopt = 0
        call cc_t2pk(work(kvpk),t2am,isymt2am,iopt)
        if (nt2sq(isymt2am).lt.nt2am(isymt2am)) 
     &    call quit('Internal Error in CCRHS_BPP')
        call dcopy (nt2am(isymt2am),work(kvpk),1,t2am,1)
      end if

      if (locdbg) then
        write(lupri,*)'Norm^2 T2 ampl. in ccrhs_bpp',
     &                 ddot(nt2am(isymt2am),t2am,1,t2am,1)
        write(lupri,*) 'T2 ampl. in ccrhs_bpp:'
        call cc_prp(one,t2am,isymt2am,0,1)
      end if

c     read V^ab_kl
      lunit = -1
      call gpopen(lunit,filev,'old',' ','unformatted',idum,.false.)
      read(lunit)(work(kvpk-1+i),i=1,nvckdl(isymvint))
      call gpclose(lunit,'KEEP')

c     calculate Omega_ijkl
      factor = one
      call cc_r12mi2(omega,work(kvpk),t2am,isymvint,isymt2am,
     &               factor,work(kend1),lwrk1)

      if (locdbg) then
        write(lupri,*)'Norm^2 Omega^ij_kl on exit in ccrhs_bpp',
     &                 ddot(ntr12sq(muld2h(isymvint,isymt2am)),
     &                      omega,1,omega,1)
        write(lupri,*) 'Omega on exit in ccrhs_bpp:'
        call cc_prsqr12(omega,muld2h(isymvint,isymt2am),'T',1,.FALSE.)
      end if

      if (.not.lt2pcked) then
c     restore T2 amplitudes as quadratic matrix
        call dcopy(nt2am(isymt2am),t2am,1,work(kvpk),1)
        call cc_t2sq(work(kvpk),t2am,isymt2am)
      end if

      call qexit('ccrhs_bpp')
      end
c=====================================================================*

c====================================================================*
      subroutine cc_r12mktbjpk(tbjpk,work,lwork)
c--------------------------------------------------------------------
c     purpose: calculate t(bj,p'k) =  
c              \sum_mn c_mn^jk * r_bp'^mn
c
c     C. Neiss, spring 2006
c--------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <dummy.h>
#include <r12int.h>
#include <ccr12int.h>

      logical locdbg
      parameter (locdbg = .false.)

      integer isymc,isymp,isymb,isymj,isymk,isymbmn,isymmn,isymjmn,
     &        isympk,isymbj,isym,isym1,isym2
      integer nramn(8),iramn(8,8),nramnq(8),iramnq(8,8)
      integer kcamp,krbmn,kcjmn,idxbj,idxpk
      integer lwork,kend1,lwrk1,kend2,lwrk2,kend3,lwrk3
      integer iopt,lunit,iadr,len,koff,koffr
      character*10 model
      character*8 framnp
#if defined (SYS_CRAY)
      real work(*), tbjpk(*), ddot, one, zero 
#else
      double precision work(*), tbjpk(*), ddot, one, zero
#endif
      parameter(zero=0.0D0, one=1.0D0)
      parameter(framnp='R12RMNAP')

      call qenter('cc_r12mktbjpk')
C
      if (locdbg) then
        write(lupri,*) 'Entered cc_r12mktbjpk'
      end if
C
C     calc. some dimensions/offsets:
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
C      
C     initialisation:
      isymc = 1
      call dzero(tbjpk,ntg2sq(isymc))
C
C     allocate memory:
      kcamp = 1
      kend1 = kcamp + ntr12am(isymc)
      lwrk1 = lwork - kend1
      if (lwrk1.lt.0) then
        call quit('Insufficient memory in CC_R12MKTBJPK')
      end if
C
C     read r12 amlitudes:
      iopt = 32
      call cc_rdrsp('R0 ',0,isymc,iopt,model,dummy,work(kcamp))
C
C     open file with r12-integrals r_(b,mn)^p'
      lunit = -1
      call wopen2(lunit,framnp,64,0)
C
C     loop over contractions:
      do isymp = 1, nsym
        isymbmn = isymp
        do p = 1, norb2(isymp)
          len = nramn(isymbmn)
          krbmn = kend1
          kend2 = krbmn + len
          lwrk2 = lwork - kend2
          if (lwrk2.lt.0) then
            call quit('Insufficient memory in CC_R12MKTBJPK')
          end if
          iadr = iramnq(isymbmn,isymp) +
     &           nramn(isymbmn)*(p-1) + 1
          call getwa2(lunit,framnp,work(krbmn),iadr,len)
C
C         if (locdbg) then
C           write(lupri,*) 'isymp, p, iadr= ',isymp,p,iadr
C           write(lupri,*) 'Norm^2 (r-int):', ddot(len,work(krbmn),1,
C    &                                     work(krbmn),1)
C           do isym2 = 1, nsym
C             isym1 = muld2h(isymbmn,isym2)
C             write(lupri,*) 'Symmetry block: ',isym1,isym2
C             call output(work(krbmn+iramn(isym1,isym2)),1,nvir(isym1),
C    &                    1,nmatkl(isym2),nvir(isym1),
C    &                    nmatkl(isym2),1,lupri)
C           end do
C         end if
C
          do isymk = 1, nsym
            isymjmn = muld2h(isymc,isymk)
            isympk  = muld2h(isymp,isymk)
            do isymj = 1, nsym
              isymmn = muld2h(isymjmn,isymj)
              isymb  = muld2h(isymmn,isymbmn)
              isymbj = muld2h(isymb,isymj)
              kcjmn = kend2
              kend3 = kcjmn + nrhfa(isymj)*nmatkl(isymmn)
              lwrk3 = lwork - kend3
              if (lwrk3.lt.0) then
                call quit('Insufficient memory in CC_R12MKTBJPK')
              end if
              do k = 1, nrhfa(isymk)
                call cc_r12sortc(work(kcjmn),work(kcamp),k,isymc,isymk,
     &                           isymj)

                idxbj = it1am(isymb,isymj) + 1
                idxpk = ig1am(isymp,isymk)+norb2(isymp)*(k-1)+p
                koffr = krbmn + iramn(isymb,isymmn)
                koff = itg2sq(isymbj,isympk) +
     &                 nt1am(isymbj)*(idxpk-1) + idxbj
                call dgemm('N','T',nvir(isymb),nrhfa(isymj),
     &                     nmatkl(isymmn),one,
     &                     work(koffr),max(nvir(isymb),1),
     &                     work(kcjmn),max(nrhfa(isymj),1),
     &                     zero,tbjpk(koff),max(nvir(isymb),1))


              end do
            end do
          end do
        end do
      end do
C
      if (locdbg) then
        write(lupri,*) 'Norm^2(c*r):',ddot(ntg2sq(isymc),tbjpk,1,
     &                                     tbjpk,1)
        do isym2 = 1, nsym
          isym1 = muld2h(isymc,isym2)
          koff = 1 + itg2sq(isym1,isym2)
          write(lupri,*) 'Symmetry block:',isym1,isym2
          call output(tbjpk(koff),1,nt1am(isym1),1,ng1am(isym2),
     &                nt1am(isym1),ng1am(isym2),1,lupri)
        end do 
      end if
C
C     close file
      call wclose2(lunit,framnp,'KEEP')
C
      if (locdbg) then
        write(lupri,*) 'Leaving cc_r12mktbjpk'
      end if
C
      call qexit('cc_r12mktbjpk')
      return
      end
c=====================================================================*

c====================================================================*
      subroutine ccrhs_e1pim(e1pim,cmox,ilamdx,xlamdh,work,lwork)
c--------------------------------------------------------------------
c     purpose: add T1-independent Fock-contribution to 
c              E1'_(a delta) intermediate 
c              and transform into MO basis
c
c     E1'_(a p') = Sum_(delta) [ CMOX_(delta p') * E1'_(a delta) + 
c                       Sum_(gamma) XLAMDH_(gamma a) * F_(gamma delta)]
c     (delta runs over MO- and aux-basis)
c
c     C. Neiss, spring 2006
c--------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <dummy.h>
#include <r12int.h>
#include <ccr12int.h>

      logical locdbg,lidxtf1
      parameter (locdbg = .false.)

      integer lwork,kend1,lwrk1,lunit,koff1,koff2,koff3
      integer kfgdp,isymfck,isymh,isyme1p,isymd,isyma,isymp
      integer isym,isym1,isym2,nbas2(8),ngdp(8),igdp(8,8),
     &        nadp(8),iadp(8,8),ilamdx(8,8),napp(8),iapp(8,8)
#if defined (SYS_CRAY)
      real work(*), e1pim(*), cmox(*), xlamdh(*), one, zero
#else
      double precision work(*), e1pim(*), cmox(*), xlamdh(*), one, zero
#endif
      parameter (one=1.0D0, zero=0.0D0)
C
      call qenter('ccrhs_e1pim')
C
C     calc. dimensions/offsets:
      do isym = 1, nsym
        nbas2(isym)  = mbas1(isym) + mbas2(isym)
      end do
      do isym = 1, nsym
        ngdp(isym) = 0
        nadp(isym) = 0
        napp(isym) = 0
        do isym2 = 1, nsym
          isym1 = muld2h(isym,isym2)
          igdp(isym1,isym2) = ngdp(isym)
          ngdp(isym) = ngdp(isym) + mbas1(isym1)*nbas2(isym2)
          iadp(isym1,isym2) = nadp(isym)
          nadp(isym) = nadp(isym) + nvir(isym1)*nbas2(isym2)
          iapp(isym1,isym2) = napp(isym)
          napp(isym) = napp(isym) + nvir(isym1)*norb2(isym2)
        end do
      end do
C
C     initialisations:
      isymh   = 1
      isymfck = 1
      isyme1p = 1
C
      kfgdp = 1
      kend1 = kfgdp + ngdp(isymfck)
      lwrk1 = lwork - kend1
      if (lwrk1.lt.0) then
        call quit('Insufficient memory in CCRHS_E1PIM')
      end if
C
C     read F_(gamma delta); delta runs over MO- and aux-basis (NBAS2)
      lunit = -1
      call gpopen(lunit,'R12FOCK','UNKNOWN',' ','UNFORMATTED',idummy,
     &            .false.)
      call readt(lunit,ngdp(isymfck),work(kfgdp))
      call gpclose(lunit,'KEEP')           
C
C     transform first index to virtual and add on E1PIM (which is 
C     half-transformed already)
      lidxtf1 = .TRUE.
      call cc_r12generaltf(work(kfgdp),e1pim,idummy,isymfck,xlamdh,
     &                     isymh,dummy,idummy,lidxtf1,iglmvi,idummy,
     &                     nbas2,nvir,nbas,igdp,idummy,idummy,
     &                     iadp,idummy,work(kend1),lwrk1)
      if (locdbg) then
        write(lupri,*)'E1PIM after addition of Fock-Matrix:'
        do isym2 = 1, nsym
          isym1 = muld2h(isyme1p,isym2)
          write(lupri,*)'Symmetry block:',isym1,isym2
          call output(e1pim(1+iadp(isym1,isym2)),
     &                1,nvir(isym1),1,nbas2(isym2),
     &                nvir(isym1),nbas2(isym2),1,lupri)
        end do
      end if
C
C     transform second index of e1pim to NORB2:
C     (reusing the work space from the Fock-matrix)
      call dzero(work(kfgdp),ngdp(isymfck))
      do isymd = 1, nsym
        isyma = muld2h(isyme1p,isymd)
        isymp = isymd !cmox is total-symmetric
        koff1 = 1 + iadp(isyma,isymd)
        koff2 = 1 + ilamdx(isymd,isymp) + nbas2(isymd)*norb1(isymp)
        koff3 = kfgdp + iapp(isyma,isymp) 
        call dgemm('N','N',nvir(isyma),norb2(isymp),nbas2(isymd),
     &             one,e1pim(koff1),max(nvir(isyma),1),
     &             cmox(koff2),max(nbas2(isymd),1),
     &             zero,work(koff3),max(nvir(isyma),1))
        if (locdbg) then
          write(lupri,*) 'CMOX:',isymd,isymp
          call output(cmox(koff2),1,nbas2(isymd),1,norb2(isymp),
     &                nbas2(isymd),norb2(isymp),1,lupri)
        end if
      end do
C
C     copy result to E1PIM:
      call dcopy(napp(isyme1p),work(kfgdp),1,e1pim,1)
C
      call qexit('ccrhs_e1pim')
      return
      end
c=====================================================================*

