*====================================================================*
      subroutine cc_r12prepccsd(work,lwork)
c--------------------------------------------------------------------
c      purpose: prepare some intermediates for the CCSD(R12) model
c               using ansatz 1
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

      logical locdbg,lvajkl,lvabkl,lv,lvijkl
      parameter (locdbg = .false.)

      integer nr1orb(8),nr1bas(8),nr1xbas(8),nr2bas,n2bst1(8)
      integer ir1xbas(8,8),ir1orb(8,8),ir1bas(8,8),ir2bas(8,8),
     &        iaodis1(8,8),ir2xbas(8,8),irgkl(8,8),nrgkl(8)
      integer nr1xorb(8),ir1xorb(8,8),nrxgkl(8),irxgkl(8,8)
      integer lwork,ibasx(8),isym,kend0,lusifc,nsymx,
     &        norbtsx,nbastx,nlamdsx,nrhfsx(8),norbsx(8),
     &        nalphaj(8),ialphaj(8,8)
      integer kcmox,lwrk0,kscr1,isymp,isymmu,klamdp,klamdh,
     &        kend1,kt1am,kend2,lwrk2,nnbastx,nbasx(8)
      integer isymtr,kvabkl,nvabkl(8),ivabkl(8,8),kend3,
     &        lwrk3,isymc,isymak,isymbl,isymab,isymkl,lunit,idum
      integer kom2,ioptbas,iopt

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

      call cc_r12offset(nr1orb,nr1xorb,nr1bas,nr1xbas,nr2bas,
     & nrgkl,nrxgkl,n2bst1,ir1orb,ir1xorb,ir1bas,ir1xbas,ir2bas,
     & ir2xbas,irgkl,irxgkl,iaodis1,nvabkl,ivabkl,nalphaj,ialphaj)
c
      ibasx(1) = 0
      do isym = 2,nsym
        ibasx(isym) = ibasx(isym-1)+mbas2(isym-1)
      end do
c
      nnbastx = 0
      do isym = 1, nsym
        nbasx(isym) = mbas1(isym) + mbas2(isym)
        nnbastx = nnbastx + nbasx(isym)*(nbasx(isym)+1)/2
      end do
c
      kend0 = 1

c----------------------------------------------------------------------
c     get MO coefficients defined for combined orbital + aux. basis:
c----------------------------------------------------------------------
      lusifc = 0
      call gpopen(lusifc,'SIRIFC','OLD',' ','UNFORMATTED',
     &            idummy,.false.)

c     read dimensions for CMO coefficients for full basis (nlamdsx)
      rewind(lusifc)
      call mollab('FULLBAS ',lusifc,lupri)
      read(lusifc) nsymx,norbtsx,nbastx,nlamdsx,(nrhfsx(i),i=1,nsym),
     &             (norbsx(i),i=1,nsym)


c     allocate work space for CMO coefficients for full basis
      kcmox  = kend0
      kend0  = kcmox + nlamdsx
      lwrk0  = lwork - kend0
      if (lwrk0 .lt.0) then
         call quit('Insufficient work space in prepccsd')
      end if

c     read CMO coefficients and close file
      read(lusifc)
      read(lusifc) (work(kcmox+i-1),i=1,nlamdsx)
      call gpclose(lusifc,'KEEP')

c     test if CMO is ok
      if (.false.) then
        write(lupri,*)'CMO out of prepccsd'
        kscr1 = 1
        do isymp = 1, nsym
          isymmu = isymp
          write(lupri,*) 'symmetry block,nbasx,norbsx:',isymp,
     &       nbasx(isymmu), norbsx(isymp)
          call output(work(kcmox+kscr1-1),1,nbasx(isymmu),
     &              1,norbsx(isymp),nbasx(isymmu),norbsx(isymp),
     &              1,lupri)
          kscr1 = kscr1 + nbasx(isymmu)*norbsx(isymp)
        end do
      end if

c     order CMO coefficients
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
      ! initialize V^kl_alpha_beta
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
        call cc_r12vtest(work(kvabkl),work(klamdp),isymc,ivabkl,
     &                   nvabkl,work(kend2),lwrk2) 
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
      call cc_r12_symv(work(kvabkl),nvabkl,ivabkl,work(kend2),lwrk2)
      lunit = 0
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
        call cc_r12vtest(work(kvabkl),work(klamdp),isymc,ivabkl,
     &                   nvabkl,work(kend2),lwrk2) 
      end if

c     calculate V^kl_alphaj and write it on file
      isymc =1
      lv = .false.
      lvajkl = .true.
      lvijkl = .false.
      call cc_r12mkvtf(work(kvabkl),dummy,work(klamdp),isymc,
     &                ivabkl,nvabkl,lv,lvijkl,lvajkl,fvajkl,
     &                work(kend2),lwrk2)

c     calculate V^ab_kl = \sum_albe L^h_ala L^h_beb V^albe_kl and
c     write it on file
      isymc = 1
      call cc_r12mkvirt(work(kvabkl),work(klamdh),isymc,fvcdkl,
     &                  work(kend2),lwrk2)

      call qexit('prepccsd')
      end 
*====================================================================*
      subroutine cc_r12mkvtf(vabkl,vijkl,cmo,isymc,ivabkl,nvabkl,
     &                       lv,lvijkl,lvajkl,filvajkl,work,lwork)
c--------------------------------------------------------------------
c     purpose: transform V^alpha_beta _kl to V^ij_kl to test V with
c              MP2-R12 
c
c     H. Fliegl, C. Haettig, summer 2004              
c--------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
#include <ccr12int.h>

      logical lres,locdbg,ldum,lprojv,lv,lvajkl,lvijkl
      parameter (locdbg = .false.)   
      integer isymc,lwork,kvjkbl,kend1,lwrk1,isymkl,isymab,isymk,
     &        isyml,isyma,isymb,isymak,isymbl,isymj,idxbl,idxakbl,
     &        koffc,ntota,ntotj,ivabkl(8,8),nvabkl(8),
     &        kvbjkl,isymbj,kvabkl,ntotb,isymjk,koffv
      integer mbas1t,icount1,imbas1(8),isym,isymaj,idxkl,luvajkl,
     &        idummy,kend2,lwrk2,lunit,ivc,
     &        isymv
      character*(*) filvajkl
#if defined (SYS_CRAY)
      real vabkl(*),vijkl(*),cmo(*),work(*),one,zero,ddot
#else
      double precision vabkl(*),vijkl(*),cmo(*),work(*),one,zero,ddot
#endif
      parameter (one = 1.0d0, zero = 0.0d0)

      call qenter('r12mkvtf')

      mbas1t  = 0
      icount1 = 0
      do isym = 1, nsym
        mbas1t = mbas1t + mbas1(isym)
        imbas1(isym) = icount1
        icount1      = icount1 + mbas1(isym)
      end do

      kvbjkl = 1
      kend1 = kvbjkl + nvajkl(1)
      lwrk1 = lwork - kend1
      if (lwrk1.lt.0) then
        call quit('insufficient work space in r12mkvtf')
      end if
      
      if (lv) then
c       read in Vabkl
        lunit = 0
        call gpopen(lunit,fvabkl,'unknown',' ','unformatted',
     &            idummy,.false.)     
        read (lunit)(vabkl(i),i=1,nvabkl(1))
        call gpclose(lunit,'KEEP')
      end if

      call dzero(work(kvbjkl),nvajkl(1))
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
                koffv  = kvbjkl+ivajkl(isymaj,isymkl)+
     &                   nt1ao(isymaj)*(idxkl-1)+it1ao(isyma,isymj)
           
                koffc  = iglmrh(isymb,isymj)+1
                ntotb = max(mbas1(isymb),1)
                ntota = max(mbas1(isyma),1)
                call dgemm('N','N',mbas1(isyma),nrhf(isymj),
     &                     mbas1(isymb),one,vabkl(kvabkl),ntota,
     &                     cmo(koffc),ntotb,one,work(koffv),ntota)    
              end do
            end do
          end do
        end do
      end do

      if (locdbg) then
        write(lupri,*)'norm^2 work(kvbjkl)', 
     &                 ddot(nvajkl(1),work(kvbjkl),1,work(kvbjkl),1)
        write(lupri,*)'work(kvbjkl)'
        do isymkl = 1, nsym
          isymaj = isymkl
          call output(work(kvbjkl+ivajkl(isymaj,isymkl)),1,
     &                nt1ao(isymaj),1,nmatkl(isymkl),
     &                nt1ao(isymaj),nmatkl(isymkl),1,lupri)
        end do
      end if
c
      if (lvajkl) then       
c       -----------------------------
c        write V(alpha j,kl) on file
c       -----------------------------
        luvajkl = 0
        call gpopen(luvajkl,filvajkl,'unknown',' ','unformatted',
     &       idummy,.false.)
        rewind(luvajkl)
        write(luvajkl) (work(kvbjkl+i-1), i = 1,nvajkl(1))
        call gpclose(luvajkl,'KEEP')
      end if
c
      if (lvijkl) then
        call dzero(vijkl,ntr12sq(1))
        lres = .true.
        call cc_r12mkvijkl(work(kvbjkl),1,cmo,isymc,work(kend1),lwrk1,
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
     &        nvabkl(8),ivabkl(8,8),nalphaj(8),ialphaj(8,8)
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
     &     ir2bas,ir2xbas,irgkl,irxgkl,iaodis1,nvabkl,ivabkl,
     &     nalphaj,ialphaj)
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
      subroutine cc_r12vtest(vabkl,cmo,isymc,ivabkl,nvabkl,
     &                       work,lwork)
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
      integer kvijkl,isymc,kend1,lwrk1,lwork,ivabkl(8,8),nvabkl(8)

#if defined (SYS_CRAY)
      real vabkl(*),work(*),cmo(*),ddot
#else
      double precision vabkl(*),work(*),cmo(*),ddot
#endif
      call qenter('r12vtest')

      kvijkl = 1
      kend1  = kvijkl + nrhftb*nrhftb*nrhftb*nrhftb 
      lwrk1  = lwork - kend1
      if (lwrk1.lt.0) then
        call quit('insufficient work space in r12vtest')
      end if
      call dzero(work(kvijkl),nrhftb*nrhftb*nrhftb*nrhftb)
      lv = .false.
      lvajkl = .false.
      lvijkl = .true.
      call cc_r12mkvtf(vabkl,work(kvijkl),cmo,
     &               isymc,ivabkl,nvabkl,
     &               lv,lvijkl,lvajkl,fvajkl,work(kend1),lwrk1)
      write(lupri,*)'Norm^2 Vijkl in r12vtest: ',
     & ddot(nrhftb*nrhftb*nrhftb*nrhftb,work(kvijkl),1,work(kvijkl),1)
      call output(work(kvijkl),1,nrhftb*nrhftb,1,nrhftb*nrhftb,
     &            nrhftb*nrhftb,nrhftb*nrhftb,1,lupri)
        
      call qexit('r12vtest')
      end  
*====================================================================*
      subroutine cc_r12_symv(vabkl,nvabkl,ivabkl,work,lwork)
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

      integer lwork,kvsym,kend1,nvabkl(8),ivabkl(8,8),lwrk1,isyma,isymb,
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
      SUBROUTINE CCRHS_BP(OMEGA2,ISYMTR,IAMP,FC12AM,LUFC12,IFILE,
     &                    BASSCL2,WORK,LWORK)
c--------------------------------------------------------------------
c     purpose: calculate \sum_mn C^ij_mn * V^mn_ab
c
c     H. Fliegl, C. Haettig, summer 2004
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
     &        isymv,nvabkl(8),ivabkl(8,8),idum,nr1orb(8),nr1bas(8),
     &        nr1xbas(8),nr2bas(8),nrgkl(8),n2bst1(8),ir1orb(8,8),
     &        ir1bas(8,8),ir1xbas(8,8),ir2bas(8,8),irgkl(8,8),
     &        iaodis1(8,8),ir2xbas(8,8),kvabkl,lunit,krpck,isymc,ntamp
      integer nr1xorb(8),ir1xorb(8,8),nrxgkl(8),irxgkl(8,8),
     &        nalphaj(8),ialphaj(8,8)
      integer isymab,iopt,iamp,lufc12,ifile,isymtr,isymij
      character*10 model
      character*8 fc12am
#if defined (SYS_CRAY)
      real work(*),one,two,factor,zero,omega2(*),
     &     four,ddot,x,basscl2
#else
      double precision work(*),one,two,
     &                 factor,zero,omega2(*),four,ddot,x,basscl2
#endif
      parameter (one = 1.0d0, two = 2.0d0, zero = 0.0d0, four = 4.0d0)

      call qenter('ccrhs_bp')

      if (locdbg) then
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

      call cc_r12offset(nr1orb,nr1xorb,nr1bas,nr1xbas,nr2bas,
     & nrgkl,nrxgkl,n2bst1,ir1orb,ir1xorb,ir1bas,ir1xbas,ir2bas,ir2xbas,
     & irgkl,irxgkl,iaodis1,nvabkl,ivabkl,nalphaj,ialphaj)
c
      krpck  = 1
      kvabkl = krpck + ntr12am(isymtr)
      kend1  = kvabkl + nvabkl(1)
      lwrk1  = lwork - kend1
      if (lwrk1 .lt.0) then
         call quit('Insufficient work space in ccrhs_bp')
      end if

c     read V 
      lunit = 0
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
      else
        call quit('Unknown value of IAMP in CCRHS_BP')
      end if

c     calculate c*V
      isymv = 1
      isymc = isymtr
      factor = one
      call cc_r12mkcv(omega2,work(kvabkl),work(krpck),
     &                isymv,isymc,factor,nvabkl,ivabkl,
     &                work(kend1),lwrk1)

      if (locdbg) then
        write(lupri,*)'in ccrhs_bp:'
        write(lupri,*)'norm^2 work(kvabkl):', 
     &   ddot(nvabkl(1),work(kvabkl),1,work(kvabkl),1)
        write(lupri,*)'norm^2 work(krpck):', 
     &   ddot(ntr12am(isymtr),work(krpck),1,work(krpck),1)
        write(lupri,*)'norm^2 omega2', 
     &   ddot(2*nt2ort(isymtr),omega2,1,omega2,1)
      end if

      if (locdbg) then
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
      subroutine cc_r12mkcv(omega2,vabkl,c2am,isymv,
     &                      isymc,factor,nvabkl,ivabkl,work,lwork)
c--------------------------------------------------------------------
c     purpose: calculate \sum_mn C^ij_mn * V^mn_ab
c
c     H. Fliegl, C. Haettig, summer 2004 
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

      integer lwork,nvabkl(8),ivabkl(8,8),isymv,isymc,isymom,isymb,
     &        isyma,isymab,isymkl,isymij,isymj,isymi,isymai,isymbj,
     &        kvabkl,kc2am,kom,kend1,lwrk1,ntoti,ntota,idxbj,idxai,
     &        koffai,maxai,naibj,index,kom2,kend2,lwrk2,idxij,idxab,
     &        idxabijp,idxabijm,idxbi,idxaj,idxaibj,idxbiaj,isymbi,
     &        isymaj,nab,kend3,lwrk3,maxj,maxb 
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

      kom2 = 1
      kend1 = kom2 + nt2ao(isymom)
      lwrk1 = lwork - kend1
      if (lwrk1.lt.0) then
        call quit('insufficient work space in r12mkcv')
      end if
           
      call dzero(work(kom2),nt2ao(isymom))

      if (locdbg) then
        x1 = 0.0d0
        xv = 0.0d0
        xc = 0.0d0
        write(lupri,*)'norm^2(c2am):',ddot(ntr12am(isymc),c2am,1,c2am,1)
        write(lupri,*)'isymom,nt2ao:',isymom,nt2ao(isymom)
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
                    work(kom2-1+naibj) = 
     &                                factor*work(kom-1+idxai-koffai)
                  end do
                else if (isymai.lt.isymbj) then
                  maxai = koffai+mbas1(isyma)*nrhf(isymi)
                  do idxai = koffai+1, maxai
                    naibj = it2ao(isymai,isymbj)+
     &                      nt1ao(isymai)*(idxbj-1) + idxai
                    work(kom2-1+naibj) =
     &                                factor*work(kom-1+idxai-koffai)
                  end do
C               else if (isymai.gt.isymbj) then
C                 maxai = koffai+mbas1(isyma)*nrhf(isymi)
C                 do idxai = koffai+1, maxai
C                   naibj = it2ao(isymbj,isymai)+
C    &                      nt1ao(isymbj)*(idxai-1) + idxbj
C                   work(kom2-1+naibj) =
C    &                                factor*work(kom-1+idxai-koffai)
C                 end do
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
        write(lupri,*)'norm^2 work(kom2):', 
     &               ddot(nt2ao(isymom),work(kom2),1,work(kom2),1)    
        call cc_prpao(dummy,work(kom2),isymom,0,1)
        write(lupri,*)'norm^2 omega2: ', 
     &   ddot(2*nt2ort(isymom),omega2,1,omega2,1)
      end if

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
c           isymom = muld2h(isymai,isymbj)
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
                    if (locdbg) then
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
 
      if (locdbg) then
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
c              expect V stored as a lower triangular matrix
c
c     H. Fliegl, C. Haettig, summer 2004
c--------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>

      integer isymb,isymv,isyma,index,isymamn,isymmn,isymm,isymn,
     &        isymam,isymbn, idxmn,idxbn,idxam,idxamn,idxambn
      integer isym,isym1,isym2,irgkl(8,8),icount,icount1,ioff(8,8)
      integer idxab,idxabmn,isymab,icount2,ivabkl(8,8)
#if defined (SYS_CRAY)
      real vs(*),vabkl(*)
#else
      double precision vs(*),vabkl(*)
#endif
      index(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j

      call qenter('sortv')

      do isym = 1, nsym
        icount = 0
        icount1 = 0
        icount2 = 0
        do isym2 = 1, nsym
          isym1 = muld2h(isym,isym2)
          ivabkl(isym1,isym2) = icount2
          icount2 = icount2 + n2bst(isym1)*nmatkl(isym2)
        end do
      end do

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
      subroutine cc_r12mkvirt(vabkl,xlamdh,isymc,filename,work,lwork)
c--------------------------------------------------------------------
c     purpose: calculate V^ab_kl = 
c              \sum_{\alpha \beta} \Lambda^h_{\alpha a} 
c              \Lambda^h_{\beta b} V^{\alpha \beta}_kl
c
c     H.Fliegl, C.Haettig, summer 2004
c--------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
#include <ccr12int.h>

      logical locdbg
      parameter (locdbg = .false.)
      integer isymc,lwork,nvabkl(8),ivabkl(8,8),icount1,isym1,isym2
      integer kvabekl,kvcdkl,kend1,isym,idum,icount2,
     &        lwrk1,isymkl,isymalbe,isymal,isymbe,isyma,isymabe,idxkl,
     &        kvabkl,koffc,ntotal,ntota,ntotbe,isymb,isymab
      integer lunit,isymk,isyml,kend2,lwrk2,kscr1,isymak,isymbl,idxab,
     &        idxakbltr,idxak,idxbl,index
      integer nck(8),ick(8,8),nvckdl(8),ivckdl(8,8)
      character*(*) filename
#if defined (SYS_CRAY)
      real vabkl(*),xlamdh(*),work(*),one,zero,ddot,x2,x3,xx,x1
#else
      double precision vabkl(*),xlamdh(*),work(*),one,zero,ddot,x2,
     &                 x3,xx,x1
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
        nvabkl(isym) = 0
        nvckdl(isym) = 0
        icount1      = 0
        icount2      = 0
        do isym2 = 1, nsym
          isym1 = muld2h(isym,isym2)
          nvabkl(isym) = nvabkl(isym) + n2bst(isym1)*nmatkl(isym2)
          ivabkl(isym1,isym2) = icount1
          icount1 = icount1 + n2bst(isym1)*nmatkl(isym2)
          if (isym2.gt.isym1) then
            nvckdl(isym) = nvckdl(isym) + nck(isym1)*nck(isym2)
            ivckdl(isym1,isym2) = icount2
            icount2 = icount2 + nck(isym1)*nck(isym2)
          else if (isym2.eq.isym1) then
            nvckdl(isym) = nvckdl(isym) + nck(isym1)*(nck(isym2)+1)/2 
            ivckdl(isym1,isym2) = icount2
            icount2 = icount2 + nck(isym1)*(nck(isym2)+1)/2
          end if
        end do
      end do

      kvcdkl  = 1
      kend1   = kvcdkl + nvckdl(1)
      lwrk1   = lwork - kend1
      if (lwrk1.lt.0) then
        call quit('insufficient work space in mkvirt')
      end if

c     transform to virtual
      call dzero(work(kvcdkl),nvckdl(1))      
      if (locdbg) then
        x2 = zero
        x3 = zero
        xx = zero
      end if

      do isymkl = 1, nsym
        isymalbe = isymkl
        isymab   = isymkl
        do isymk = 1, nsym
          isyml = muld2h(isymkl,isymk)

          do k = 1, nrhfb(isymk)
            do l = 1, nrhfb(isyml)
              idxkl = imatkl(isymk,isyml)+nrhfb(isymk)*(l-1)+k
              
              do isymal = 1, nsym
c               isyma = muld2h(isymc,isymal)
                isyma = isymal
                isymbe = muld2h(isymalbe,isymal)
c               isymb  = muld2h(isymc,isymbe)
                isymb  = isymbe
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
     &                    nbas(isymal),one,xlamdh(koffc),ntotal,
     &                    vabkl(kvabkl),ntotal,zero,work(kvabekl),
     &                    ntota)
                
                koffc  = iglmvi(isymbe,isymb)+1
                ntotbe = max(nbas(isymbe),1)
                ntota  = max(nvir(isyma),1) 
              
                call dgemm('N','N',nvir(isyma),nvir(isymb),nbas(isymbe),
     &                one,work(kvabekl),ntota,xlamdh(koffc),ntotbe,zero,
     &                work(kscr1),ntota)

                if (locdbg) then
                  x2 = x2 + ddot(nvir(isyma)*nvir(isymb),
     &                    work(kscr1),1,work(kscr1),1)
                end if

c               sort result as Vakbl and pack as triangular matrix
                isymak = muld2h(isyma,isymk)
                isymbl = muld2h(isymb,isyml)

                if (isymak.eq.isymbl) then
                 do a = 1, nvir(isyma)
                  do b = 1, nvir(isymb)
                    idxab = nvir(isyma)*(b-1)+a
                    idxak = ick(isyma,isymk)+nvir(isyma)*(k-1)+a
                    idxbl = ick(isymb,isyml)+nvir(isymb)*(l-1)+b
                    idxakbltr = ivckdl(isymak,isymbl)+index(idxak,idxbl)
                    work(kvcdkl-1+idxakbltr) = work(kscr1-1+idxab)

                    if (locdbg) then
                      if (idxak.eq.idxbl) then
c                       sum diagonal elements
                        x3 = x3 + ddot(1,work(kvcdkl-1+idxakbltr),1,
     &                            work(kvcdkl-1+idxakbltr),1)
                        
                      end if
                    end if
                  end do
                 end do
                else
                  call quit('symmetry error in mkvirt')
                end if
c
              end do
            end do
          end do

       end do
      end do

      if (locdbg) then
        write(lupri,*)'x2 =', x2
        write(lupri,*)'x3 =', x3
        write(lupri,*)'nvckdl(1) =', nvckdl(1)
        write(lupri,*)'xx =', xx

c       write(lupri,*)'norm^2 Vabkl =', 
c    &               ddot(nvckdl(1),work(kvcdkl),1,work(kvcdkl),1)
        write(lupri,*)'norm^2 Vabkl =', 
     &  2.0d0*ddot(nvckdl(1),work(kvcdkl),1,work(kvcdkl),1)-x3
        write(lupri,*)'norm^2 Vabkl triang =', 
     &               ddot(nvckdl(1),work(kvcdkl),1,work(kvcdkl),1)

        call dreinorm(work(kvcdkl),nvckdl(1),nck(1),nsym,x1)
        write(lupri,*)'dreinorm =', x1
        do isymak = 1, nsym
          isymbl = isymak
          call outpkb(work(kvcdkl+ivckdl(isymak,isymbl)),nck(isymak),
     &              1,1,lupri)
        end do
      end if

c     write Vabkl on file
      lunit = 0
      call gpopen(lunit,filename,'unknown',' ','unformatted',idum,
     &            .false.)
      write(lunit)(work(kvcdkl-1+i),i=1,nvckdl(1))
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
     &                     filev,work,lwork)
c--------------------------------------------------------------------
c     purpose: calculate B''-Term and update Omega_ijkl
c              \sum_ab V^ab_kl * t^ij_ab
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

      logical locdbg,lt2pcked
      parameter (locdbg = .false.)

      integer isymt2am,lwork,kt2pk,kvpk,lwrk1,kend1,iopt,lunit,
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

      isymvint = 1
      kt2pk = 1
      kvpk  = kt2pk + nt2am(isymt2am)
      kend1 = kvpk  + nvckdl(isymvint)
      lwrk1 = lwork - kend1
      if (lwrk1.lt.0) then
        call quit('insufficient work space in ccrhs_bpp')
      end if

      if (.not.lt2pcked) then
c     pack T2 amplitudes as a triangular matrix
        iopt = 0
        call cc_t2pk(work(kt2pk),t2am,isymt2am,iopt)
      else 
        call dcopy(nt2am(isymt2am),t2am,1,work(kt2pk),1)
      end if

c     read V^ab_kl
      lunit = 0
      call gpopen(lunit,filev,'unknown',' ','unformatted',idum,.false.)
      read(lunit)(work(kvpk-1+i),i=1,nvckdl(isymvint))
      call gpclose(lunit,'KEEP')

c     calculate Omega_kilj
      factor = one
      call cc_r12mi2(omega,work(kvpk),work(kt2pk),isymvint,isymt2am,
     &               factor,work(kend1),lwrk1)

      if (locdbg) then
        write(lupri,*)'Norm^2 Omega^ij_kl in ccrhs_bpp',
     &                 ddot(ntr12sq(muld2h(isymvint,isymt2am)),
     &                      omega,1,omega,1)
      end if

      if (.not.lt2pcked) then
c     pack T2 amplitudes back to a quadratic matrix
        call cc_t2sq(work(kt2pk),t2am,isymt2am)
      end if

      call qexit('ccrhs_bpp')
      end
c====================================================================*
