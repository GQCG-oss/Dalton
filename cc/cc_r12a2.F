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
*=====================================================================*
      subroutine cc_r12whtf(xint,idel,isymd,j,isymj,isymbi,lauxbeta,
     &                      lunitr12,filer12,lunitr12_2,filer12_2,
     &                      work,lwork)
c----------------------------------------------------------------------
c     purpose: get half transformed r12 integrals and write them
c              on file first delta then delta' in ccsd_iabj2   
c 
c     xint     array of half transformed 12 integrals 
c              (i beta|r12|j delta) with j, delta fixed
c     idel     number of delta (running over all symmetries; in 
c              contrast, idelta starts with 1 for every symmetry!)
c     isymd    symmetry of delta
c     j        number of j (starts with 1 for every symmetry)
c     isymj    symmetry of j
c     isymbi   symmetry of (i beta)
c     lauxbeta flag to write out also integrals with beta being
c              auxiliary basis function (needed for R12 response)
c     lunitr12, filer12 unit number, file name to be written to
c     lunitr12_2, filer12_2 unit number, file name for integrals with
c              two auxiliary basis functions
c     
c     H. Fliegl, C. Haettig, winter 2003
c     modified C. Neiss, fall 2004: write also integrals with beta 
c                auxiliary basis functions if desired
c----------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
#include <dummy.h>
       
      logical locdbg,lauxbeta
      parameter (locdbg = .false.)
      integer lunitr12, lunitr12_2, lwork, idel, isymd, isymi, isymb
      integer isymbi, isymdj, isymj, ndelj, iadr, ilen
      integer isym, isym1, isym2, icount2, icount3, icount4, icount8
      integer idxbi, idxdj, idxbi2, koff1
      integer nr1bas(8),nr1xbas(8),nr2bas,nr2bas2
      integer ir1bas(8,8), ir1xbas(8,8), ir2bas(8,8)
      integer ir2xbas(8,8),idelta
      integer ir2bas2(8,8),ir2xbas2(8,8)
      integer nrhf1(8)

      character*(*) filer12, filer12_2

#if defined (SYS_CRAY)
      real work(*), xint(*)
#else
      double precision work(*), xint(*)
#endif

      call qenter('r12whtf')
      
      if (locdbg) then
        write(lupri,*) 'Entered CC_R12WHTF'
        write(lupri,*) 'lauxbeta,r12int = ',lauxbeta,r12int
        write(lupri,'(A,I5,2X,A)')
     &    ' LUNITR12,   FILER12  :',lunitr12,filer12,
     &    ' LUNITR12_2, FILER12_2:',lunitr12_2,filer12_2
        write(lupri,*) 'V12INT = ',V12INT
      end if

      if (V12INT) then
        ! use active occupied orbitals
        call icopy(8,nrhfa,1,nrhf1,1)
      else
        ! use r12 orbitals
        call icopy(8,nrhfb,1,nrhf1,1)
      end if

      if (j.gt.nrhf1(isymj)) then
        call qexit('r12whtf')
        return
      end if

      do isym = 1, nsym
        nr1bas(isym)  = 0
        nr1xbas(isym) = 0
        do isym2 = 1, nsym
           isym1 = muld2h(isym,isym2)
           nr1bas(isym) = nr1bas(isym) + mbas1(isym1)*nrhf1(isym2)
           nr1xbas(isym)= nr1xbas(isym)+ mbas2(isym1)*nrhf1(isym2)
        end do
      end do
      nr2bas = 0
      do isym = 1, nsym
         nr2bas = nr2bas + nr1bas(isym)*nr1bas(isym)
      end do
      do isym = 1, nsym
        icount2 = 0
        icount3 = 0
        icount4 = 0
        icount8 = 0
        do isym2 = 1, nsym
          isym1 = muld2h(isym,isym2)
          ir1bas(isym1,isym2) = icount2
          ir1xbas(isym1,isym2)= icount8
          ir2bas(isym1,isym2) = icount3
          ir2xbas(isym1,isym2)= icount4
          icount2 = icount2 + nrhf1(isym2)*mbas1(isym1)
          icount8 = icount8 + nrhf1(isym2)*mbas2(isym1)
          icount3 = icount3 + nr1bas(isym1)*nr1bas(isym2)
          icount4 = icount4 + nr1bas(isym1)*nr1xbas(isym2)
        end do
      end do

      if (lwork.lt.nr1bas(isymbi)) then
        write(lupri,*) 'lwork, nr1bas(isymbi):',lwork, nr1bas(isymbi)
        call quit('Insufficient work space in cc_r12whtf')
      end if

      idelta = idel - ibas(isymd)  

      koff1 = 1
      isymdj = muld2h(isymd,isymj)
      do isymi = 1, nsym
         isymb  = muld2h(isymbi,isymi)
        do i = 1, nrhf1(isymi)
            do b = 1, mbas1(isymb) 
              idxbi   = ir1bas(isymb,isymi) + mbas1(isymb)*(i-1)+b 
              idxbi2  = koff1 -1 + nbas(isymb)*(i-1) +b
              work(idxbi) = xint(idxbi2) 
            end do 
        end do
        if (locdbg) then
          write(lupri,*)'idel, j, isymdj =',idel,j,isymdj
          write(lupri,*) 'symmetry block (beta,i):',isymb,isymi
          call around('integrals before reordering')
          call output(xint(koff1),1,nbas(isymb),1,nrhf(isymi),
     &                nbas(isymb),nrhf(isymi),1,lupri)
          call around('integrals after reordering')
          call output(work(ir1bas(isymb,isymi)+1),1,mbas1(isymb),
     &                1,nrhf1(isymi),mbas1(isymb),nrhf1(isymi),1,lupri)
          write(lupri,*)
        end if 
        koff1 = koff1 + nbas(isymb)*nrhf(isymi)
      end do
      
     
      if (idelta.le.mbas1(isymd)) then
        ndelj = ir1bas(isymd,isymj) + mbas1(isymd)*(j-1) + idelta 
        iadr  = ir2bas(isymbi,isymdj) + nr1bas(isymbi)*(ndelj-1)+1
      else 
        ndelj = ir1xbas(isymd,isymj) + 
     &          mbas2(isymd)*(j-1) + idelta - mbas1(isymd) 
        iadr  = nr2bas + ir2xbas(isymbi,isymdj) + 
     &          nr1bas(isymbi)*(ndelj-1)+1 
      end if
      ilen = nr1bas(isymbi)
      call putwa2(lunitr12,filer12,work,iadr,ilen)

C     write also integrals with beta being auxiliary function on extra
C     file if desired:      
      if (lauxbeta) then
        if (mbsmax.lt.6) then
          write (lupri,*) 'MBSMAX = ',MBSMAX
          call quit("MBSMAX must at least be 6 for "//
     &              "(k beta'|r12|l delta')")
        end if

        ! nr2bas2 = lenghts for nr1bas x nr1xbas over all symmetries 
        nr2bas2 = 0
        do isym = 1, nsym
           nr2bas2 = nr2bas2 + nr1bas(isym)*nr1xbas(isym)
        end do
        ! ir2bas2  = symmetry offsets for nr1xbas x nr1bas
        ! ir2xbas2 = symmetry offsets for nr1xbas x nr1xbas
        do isym = 1, nsym
          icount3 = 0
          icount4 = 0
          do isym2 = 1, nsym
            isym1 = muld2h(isym,isym2)
            ir2bas2(isym1,isym2)  = icount3 
            ir2xbas2(isym1,isym2) = icount4
            icount3 = icount3 + nr1xbas(isym1)*nr1bas(isym2)
            icount4 = icount4 + nr1xbas(isym1)*nr1xbas(isym2)
          end do
        end do

        ! reallocate work space:
        if (lwork.lt.nr1xbas(isymbi)) then
          call quit('Insufficient work space in cc_r12whtf')
        end if   

        koff1 = 1 
        isymdj = muld2h(isymd,isymj)
        do isymi = 1, nsym
          isymb  = muld2h(isymbi,isymi)
          do i = 1, nrhf1(isymi)
            do b = mbas1(isymb) + 1, mbas1(isymb)+mbas2(isymb)
              idxbi   = ir1xbas(isymb,isymi) + mbas2(isymb)*(i-1) + b -
     &                  mbas1(isymb)
              idxbi2  = koff1 -1 + nbas(isymb)*(i-1) + b
              work(idxbi) = xint(idxbi2)
            end do
          end do
          if (locdbg) then
            call around('integrals after reordering')
            call output(work(ir1xbas(isymb,isymi)+1),1,mbas2(isymb),
     &                1,nrhf1(isymi),mbas2(isymb),nrhf1(isymi),1,lupri)
            write(lupri,*)
          end if
          koff1 = koff1 + nbas(isymb)*nrhf(isymi)
        end do

        if (idelta.le.mbas1(isymd)) then
          ndelj = ir1bas(isymd,isymj) + mbas1(isymd)*(j-1) + idelta 
          iadr  = ir2bas2(isymbi,isymdj) + nr1xbas(isymbi)*(ndelj-1)+1
        else
          ndelj = ir1xbas(isymd,isymj) +
     &            mbas2(isymd)*(j-1) + idelta - mbas1(isymd)
          iadr  = nr2bas2 + ir2xbas2(isymbi,isymdj) +
     &            nr1xbas(isymbi)*(ndelj -1) + 1
        end if
        ilen   = nr1xbas(isymbi)
        call putwa2(lunitr12_2,filer12_2,work,iadr,ilen)
      end if

      if (locdbg) then
        write(lupri,*) 'Leaving CC_R12WHTF'
      end if

      call qexit('r12whtf')
      end 
*=====================================================================*
      subroutine cc_r12intf2(vitjtkl,lamdh,isymh,lamdhs,isymv,
     &                       lamdhcs,isymhcs,work,lwork)
c----------------------------------------------------------------------
c     purpose: calculate V^itjt_kl for ansatz 2
c     
c     H. Fliegl, C. Haettig winter 2003
c
c     \sum_MN r^MtNt_kl * g^alpha jt_MN
c----------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
#include <ccr12int.h>
#include <dummy.h>

      logical ldum,lauxd,locdbg
      parameter (locdbg = .false.)

      integer isymv,krmtntkl,kgmtntkl,nt2aox,isymhs,isymhcs       
      integer krhtf, kend1, lwrk1, lwork, idum, idelta, isymd
      integer isymk, isymdk, isymm, isymij, isymkl, isyml, koffv, 
     &        kend2, koffr, koffg, isymnm, ntotmn, ntotij
      integer nr1orb(8), nr1bas(8), nr1xbas(8), nr2bas, n2bst1(8)
      integer ir1xbas(8,8),ibasx(8),isymdm,idxkl,idxlk,iadr1,iadr2
      integer ir1orb(8,8), ir1bas(8,8), ir2bas(8,8), iaodis1(8,8)
      integer ir2xbas(8,8), irgkl(8,8), irxgkl(8,8)
      integer nrgkl(8), nrxgkl(8), lwrk2
      integer idxd,idxdp,isym,icount1,isym1,isym2,kgtf,isyma,koffc,
     &        koff1,koff2,ntota,ntotm,isymn,idxmn,idxdn
      integer iglmrhs(8,8),iglmvis(8,8),nglmds(8),igamsm2(8,8),
     &        icmo(8,8),ncmo(8),imaklm(8,8),nmaklm(8),ngamsm2(8),kgmkl,
     &        imatijm(8,8),nmatijm(8),icount5,ngamsm(8),igamsm(8,8),
     &        irgijs(8,8),nrgijs(8),ir1basm(8,8),nr1basm(8),
     &        ir2basm(8,8),nr2basm,ir1xbasm(8,8),nr1xbasm(8),nvajkls(8),
     &        ir2xbasm(8,8),imatf(8,8),nmatf(8),icount2,ivajkls(8,8)
      integer nr1xorb(8),ir1xorb(8,8),
     &        nalphaj(8),ialphaj(8,8),nmaijm(8),imaijm(8,8)
      integer isymaj,ntotaj,isymh,lunit,krsym,isymg,isymr

      double precision one, zero,two
      parameter (one = 1.0d0, zero = 0.0d0, two = 2.0d0)
#if defined (SYS_CRAY)
      real lamdhs(*), vitjtkl(*), work(*),ddot,lamdhcs(*),
     &     factor, lamdh(*),x1,x2,x3
#else
      double precision lamdhs(*), vitjtkl(*), work(*), ddot,lamdhcs(*),
     &                 factor, lamdh(*),x1,x2,x3
#endif

      call qenter('r12intf2')

      call cc_r12offset(nr1orb,nr1xorb,nr1bas,nr1xbas,nr2bas,
     &     nrgkl,nrxgkl,n2bst1,ir1orb,ir1xorb,ir1bas,ir1xbas,
     &     ir2bas,ir2xbas,irgkl,irxgkl,iaodis1,nalphaj,ialphaj) 
      call cc_r12offs23(iglmrhs,iglmvis,nglmds,icmo,ncmo,
     &             imaijm,nmaijm,imaklm,nmaklm,
     &             imatijm,nmatijm,igamsm,ngamsm,irgijs,nrgijs,
     &             ir1basm,nr1basm,ir2basm,nr2basm,ir1xbasm,
     &             nr1xbasm,ir2xbasm,imatf,nmatf)

      if (locdbg) then
        write(lupri,*)'V before tf2'
        do isymij = 1, nsym
          isymkl = muld2h(isymv,isymij)
          write(lupri,*)'symmetry block (ij,kl)',isymij,isymkl
          call output(vitjtkl(1+itr12sqt(isymij,isymkl)),
     &              1,nmatij(isymij),1,nmatkl(isymkl),
     &              nmatij(isymij),nmatkl(isymkl),1,lupri)
        end do
        write(lupri,*)'norm^2 V before tf2',
     &        ddot(ntr12sq(isymv),vitjtkl,1,vitjtkl,1)
      end if

      ibasx(1) = 0
      do isym = 2, nsym
        ibasx(isym) = ibasx(isym-1) + mbas2(isym-1)
      end do

      do isym = 1, nsym
        icount1 = 0
        icount2 = 0
        do isym2 = 1, nsym
          isym1 = muld2h(isym,isym2)
          igamsm2(isym1,isym2) = icount1
          ivajkls(isym1,isym2) = icount2
          icount1 = icount1 + nmatij(isym1)*nmatijm(isym2)
          icount2 = icount2 + nalphaj(isym1)*nmatijm(isym2)
        end do
        ngamsm2(isym) = icount1
        nvajkls(isym) = icount2
      end do

      isymhs = isymh
      isymg = muld2h(isymh,isymh)
      isymr = muld2h(isymhcs,isymhs)

      krmtntkl = 1
      kgmtntkl = krmtntkl + ngamsm(isymr)
      kend1 = kgmtntkl + ngamsm2(isymg)
      lwrk1   = lwork - kend1
      if (lwrk1.lt.0) then
        write(lupri,*)'lwork,lwrk1,kend1',lwork,lwrk1,kend1
        call quit('Insufficient work space in cc_r12intf2')
      end if

      call dzero(work(krmtntkl),ngamsm(isymr))
      call dzero(work(kgmtntkl),ngamsm2(isymg))

      do isymd = 1, nsym
        krhtf = kend1
        kend2 = krhtf + max(nrgkl(isymd),nrgijs(isymd))
        lwrk2 = lwork - kend2
        if (lwrk2.lt.0) then
          write(lupri,*)'lwrk2, lwork, kend2', lwrk2, lwork, kend2
          call quit('Insufficient work space in cc_r12intf2')
        end if

        do idelta = 1, mbas1(isymd)
c         read (ka|r12|lb) from file and sort as R^d_(ak,l) 
          lauxd = .false.
          idxd = idelta+ibas(isymd)
          call cc_r12getrint(work(krhtf),idxd,isymd,nr1bas,ir1bas,
     &         nr2bas,ir2bas,nrgkl,irgkl,ir1xbas,ir2xbas,
     &         nrhfb,nmatkl,imatkl,
     &         ibasx,lauxd,.false.,frhtf,work(kend2),lwrk2)
          if (locdbg) then
            write(lupri,*)'Rhtf in tf2'
            write(lupri,*)' norm^2 Rhtf in tf2',
     &                ddot(nrgkl(isymd),work(krhtf),1,work(krhtf),1)
          end if

c         case response calculation: lres = .true.
c         get rmtntkl
          call cc_r12generaltf(work(krhtf),work(krmtntkl),idelta,
     &                       isymd,lamdhcs,isymhcs,lamdhs,isymhs,
     &                       .false.,iglmrhs,iglmrhs,nmatkl,nrhfsa,
     &                       nbas,irgkl,imatijm,nmatijm,imaklm,
     &                       igamsm,work(kend2),lwrk2)

c         read (Ka|Lb) from file and sort as I^d_(aK,L)
          lauxd = .false.
          idxd = idelta+ibas(isymd)
          call cc_r12getrints(work(krhtf),idxd,isymd,nr1basm,ir1basm,
     &           nr2basm,ir2basm,nrgijs,irgijs,ir1xbasm,ir2xbasm,
     &           ibasx,imatijm,nmatijm,lauxd,fccgmnab,work(kend2),lwrk2)

c         get gmtntKL
          call cc_r12generaltf(work(krhtf),work(kgmtntkl),idelta,
     &                         isymd,lamdh,isymh,lamdh,isymh,
     &                         .false.,iglmrh,iglmrh,nmatijm,nrhf,
     &                         nbas,irgijs,imatij,nmatij,imatf,
     &                         igamsm2,work(kend2),lwrk2)

        end do
      end do
c
      if (locdbg) then
c       print result
        write(lupri,*)'Rmtntkl',
     &        ddot(ngamsm(isymr),work(krmtntkl),1,work(krmtntkl),1)
        write(lupri,*)'Imtntkl',
     &        ddot(ngamsm2(isymg),work(kgmtntkl),1,work(kgmtntkl),1)
      end if

      x1 = 0.0d0
      x2 = 0.0d0
      x3 = 0.0d0
c
c     make Vitjtkl
c
      factor = one
      do isymnm = 1, nsym
        isymij = isymnm
        isymkl = muld2h(isymv,isymij)

        koffv = 1 + itr12sqt(isymij,isymkl)
        koffr = krmtntkl + igamsm(isymnm,isymkl)
        koffg = kgmtntkl + igamsm2(isymij,isymnm)
      
        ntotmn = max(nmatijm(isymnm),1)
        ntotij = max(nmatij(isymij),1)
 
        call dgemm('N','N',nmatij(isymij),nmatkl(isymkl),
     &           nmatijm(isymnm),factor,work(koffg),ntotij,
     &           work(koffr),ntotmn,one,vitjtkl(koffv),ntotij)

        if (locdbg) then
          x1 = x1 + 
     &    ddot(nmatijm(isymnm)*nmatij(isymij),work(koffg),1,
     &         work(koffg),1)
          x2 = x2 + 
     &    ddot(nmatijm(isymnm)*nmatkl(isymkl),work(koffr),1,
     &         work(koffr),1)
          x3 = x3 +
     &    ddot(nmatij(isymij)*nmatkl(isymkl),vitjtkl(koffv),1,
     &         vitjtkl(koffv),1)
        end if

        end do
c
      if (locdbg) then
        write(lupri,*)'V after tf2'
        do isymij = 1, nsym
          isymkl = muld2h(isymv,isymij)
          write(lupri,*)'symmetry block (ij,kl)',isymij,isymkl
          call output(vitjtkl(1+itr12sqt(isymij,isymkl)),
     &              1,nmatij(isymij),1,nmatkl(isymkl),
     &              nmatij(isymij),nmatkl(isymkl),1,lupri)
        end do
        write(lupri,*)'norm^2 V after tf2',
     &        ddot(ntr12sq(isymv),vitjtkl,1,vitjtkl,1)
        write(lupri,*)'x1 = ', x1
        write(lupri,*)'x2 = ', x2
        write(lupri,*)'x3 = ', x3
      end if 
c
      call qexit('r12intf2')
      end 
*======================================================================*
      subroutine ccr12prep2(work,lwork)
*----------------------------------------------------------------------*
c     Purpose: make some intermediates used for ansatz 2 & 3
c
c     H. Fliegl, C. Haettig, winter 2003
c     adapted for ansatz 3, Christof Haettig, spring 2005
c     adapted for CABS, Christian Neiss, winter 2006
*----------------------------------------------------------------------*
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
#include <ccr12int.h>
#include <dummy.h>
       
#if defined (SYS_CRAY)
      real zero,one,half
#else
      double precision zero,one,half
#endif
      parameter (zero = 0.0D0, one = 1.0D0, half = 0.5d0)
 
      logical noauxbsave, locdbg, lauxd, do_virt
      parameter (locdbg = .false.)

      integer nr1orb(8),nr1bas(8),nr1xbas(8),nr2bas,n2bst1(8)
      integer ir1xbas(8,8),basoff(8)
      integer ir1orb(8,8),ir1bas(8,8),ir2bas(8,8),iaodis1(8,8)
      integer ir2xbas(8,8),irgkl(8,8),nrgkl(8),ibasx(8)
      integer nbasx(8),nnbastx,isym,isymdk,n2bst2(8),isymd,isymk,
     &        icount,kend0,lusifc,norbsx(8),nbastx,nlamdsx,
     &        nrhfsx(8),kcmox,lwrk0,lwork,ksfull
      integer norbtsx,isympq,isymp,isymq,ksjbp,ksbbp,ksmupbp
      integer isymbp,krhtf,ib,isymdp,idelp,idxbpdp,krdpakl
      integer isymakl,kend1,lwrk1,koff1,klamdp,klamdh,kend2,lwrk2,
     &        isymj,koffs,kcmo,ntotc,ntots,ntots2,kgdakl,
     &        kghtf,idxddp,koff2,idel,kgdpakl,idxd,idxdp 
      integer kt1am,isymmu,kscr1,kstep,luvajkl,icoun1
      integer kiajb, idxb,krtot, icount7, kvajkl, kend3,lwrk3
      integer iglmrhs(8,8),iglmvis(8,8),nglmds(8),icmo(8,8),ncmo(8),
     &        imaklm(8,8),nmaklm(8),imatijm(8,8),nmatijm(8),igamsm(8,8),
     &        ngamsm(8),klamdhs,klamdps,irgijs(8,8),nrgijs(8),
     &        ir1basm(8,8),nr1basm(8),ir2basm(8,8),nr2basm,
     &        ir1xbasm(8,8),nr1xbasm(8),ir2xbasm(8,8),imatf(8,8),
     &        nmatf(8),isbbp(8,8),nsbbp(8),isjbp(8,8),nsjbp(8),nbas2(8)
      integer nr1xorb(8),ir1xorb(8,8),nrxgkl(8),irxgkl(8,8),
     &        nalphaj(8),ialphaj(8,8),nmaijm(8),imaijm(8,8),
     &        nrgij(8),irgij(8,8),nt1xao(8),it1xao(8,8),it2xaos(8,8)
#if defined (SYS_CRAY)
      real work(*),ddot,xs
#else
      double precision work(*),ddot,xs
#endif

c----------------------------------------------------------------------
c     some initializations: 
c----------------------------------------------------------------------
      call qenter('r12prep2')
      if(locdbg) write(lupri,*) 'entered ccr12prep2'

      call cc_r12offset(nr1orb,nr1xorb,nr1bas,nr1xbas,nr2bas,
     & nrgkl,nrxgkl,n2bst1,ir1orb,ir1xorb,ir1bas,ir1xbas,ir2bas,
     & ir2xbas,irgkl,irxgkl,iaodis1,nalphaj,ialphaj) 

      call cc_r12offs23(iglmrhs,iglmvis,nglmds,icmo,ncmo,
     &             imaijm,nmaijm,imaklm,nmaklm,
     &             imatijm,nmatijm,igamsm,ngamsm,irgijs,nrgijs,
     &             ir1basm,nr1basm,ir2basm,nr2basm,ir1xbasm,
     &             nr1xbasm,ir2xbasm,imatf,nmatf)
c
      if (r12cbs) then
        do isym = 1, nsym
          nbas2(isym)  = mbas1(isym) + mbas2(isym)
          basoff(isym) = 0
        end do
      else
        call icopy(8,mbas2,1,nbas2,1)
        call icopy(8,mbas1,1,basoff,1)
      end if
c
c     calculate number of auxbasis functions in previous symmetries   
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
c     do isymdk = 1, nsym
c       n2bst2(isymdk) = 0
c       do isymk = 1, nsym
c          isymd = muld2h(isymdk,isymk)
c          n2bst2(isymdk) = n2bst2(isymdk) + mbas2(isymd)*mbas2(isymk)
c       end do
c     end do
c
c     do isymdk = 1, nsym
c       icount = 0
c       do isymk = 1, nsym
c         isymd = muld2h(isymdk,isymk)
c         iaodis2(isymd,isymk) = icount
c         icount = icount + mbas2(isymd)*mbas2(isymk) 
c       end do
c     end do
c
      do isympq = 1, nsym
        nsjbp(isympq)= 0
        nsbbp(isympq)= 0
        do isymp = 1, nsym
          isymq = muld2h(isympq,isymp)
          isjbp(isymp,isymq) = nsjbp(isympq)
          nsjbp(isympq) = nsjbp(isympq)+nrhfsa(isymp)*nbas2(isymq)
          isbbp(isymp,isymq) = nsbbp(isympq)
          nsbbp(isympq) = nsbbp(isympq)+nvir(isymp)*nbas2(isymq)
        end do
      end do
c
      kend0 = 1

c----------------------------------------------------------------------
c     get MO coefficients defined for combined orbital + aux. basis:
c----------------------------------------------------------------------
c     lusifc = -1
c     call gpopen(lusifc,'SIRIFC','OLD',' ','UNFORMATTED',
c    &            idummy,.false.)
c
c     read dimensions for CMO coefficients for full basis (nlamdsx)
c     rewind(lusifc)
c     call mollab('FULLBAS ',lusifc,lupri)
c     read(lusifc) nsymx,norbtsx,nbastx,nlamdsx,(nrhfsx(i),I=1,nsym),
c    &             (norbsx(i),i=1,nsym)
c
c
c     allocate work space for CMO coefficients for full basis
c     kcmox  = kend0
c     kend0  = kcmox + nlamdsx
c     lwrk0  = lwork - kend0
c     if (lwrk0 .lt.0) then
c        call quit('Insufficient work space in ccr12prep2')
c     end if
c
c     read CMO coefficients and close file
c     read(lusifc)
c     read(lusifc) (work(kcmox+i-1),i=1,nlamdsx)
c     call gpclose(lusifc,'KEEP')
c
c     test if CMO is ok
c     if (locdbg) then
c       write(lupri,*)'CMO out of ccr12prep2'
c       kscr1 = 1 
c       do isymp = 1, nsym
c         isymmu = isymp 
c         write(lupri,*) 'symmetry block,nbasx,norbsx:',isymp,
c    &       nbasx(isymmu), norbsx(isymp)
c         call output(work(kcmox+kscr1-1),1,nbasx(isymmu),
c    &              1,norbsx(isymp),nbasx(isymmu),norbsx(isymp),
c    &              1,lupri)
c         kscr1 = kscr1 + nbasx(isymmu)*norbsx(isymp)
c       end do
c     end if

c     order CMO coefficients
      klamdp = kend0
      klamdh = klamdp + nlamdt  
      kend1  = klamdh + nlamdt
      kt1am  = kend1
      kend2  = kt1am + nt1amx
      lwrk2  = lwork - kend2
      if(lwrk2.lt.0) then
        call quit('insufficient work space in cc_r12prep2')
      end if      
      call dzero(work(kt1am),nt1amx) 
      call lammat(work(klamdp),work(klamdh),work(kt1am),
     &            work(kend2),lwrk2)

c-----------------------------------------------------------------------
c     read overlap matrix for full basis from file:
c-----------------------------------------------------------------------

      ksfull  = kend1
      kend1   = ksfull + nnbastx
      lwrk1   = lwork  - kend1
      if (lwrk1 .lt.0) then
         call quit('Insufficient work space in ccr12prep2')
      end if

      noauxbsave = noauxb
      noauxb = .false. ! read integrals for the full basis
      call rdonel('OVERLAP ',.true.,work(ksfull),nnbastx) 
      noauxb = noauxbsave

      if (locdbg) then
        !test if S is ok compare with sirius/r12aux.F !
        write(lupri,*)'in r12prep nnbastx = ', nnbastx
        write(lupri,*)'in r12prep S'
        call outpkb(work(ksfull),nbasx,1,1,lupri)
        write(lupri,*)'norm^2 ksfull:', 
     &                 ddot(nnbastx,work(ksfull),1,work(ksfull),1)
      end if

c     get S_Jb' and for ansatz 3 also S_bb'
      ksjbp   = kend1
      kend1   = ksjbp + nsjbp(1)

      ksbbp = kend1
      if (ianr12.eq.3) then
        kend1 = ksbbp + nsbbp(1)
        do_virt = .true.
      else
        do_virt = .false.
      end if

      klamdhs = kend1
      klamdps = klamdhs + nglmds(1)
      kend1   = klamdps + nglmds(1) 
      lwrk1   = lwork - kend1
      if (lwrk1 .lt.0) then
        call quit('Insufficient work space in ccr12prep2')
      end if

      call lammats(work(klamdps),work(klamdhs),dummy,
     &             1,.true.,.false.,
     &             nglmds,iglmrhs,iglmvis,icmo,work(kend1),lwrk1)       

      call cc_r12mksjbp(work(ksfull),work(ksjbp),work(ksbbp),do_virt,
     &                  work(klamdhs),nbasx,nbas2,iglmrhs,iglmvis,
     &                  isjbp,isbbp,work(kend1),lwrk1)
      
      if (locdbg) then
        write(lupri,*)'norm^2 S_Jbp', 
     &                 ddot(nsjbp(1),work(ksjbp),1,work(ksjbp),1)
      end if

c----------------------------------------------------------------------
c     make I^dp_a,kl 
c----------------------------------------------------------------------
c     calculate new offsets:
      do isymdk = 1, nsym
        icount = 0
        icoun1 = 0
        do isymk = 1, nsym
          isymd = muld2h(isymdk,isymk)
          irgij(isymd,isymk)  = icount
          it1xao(isymd,isymk) = icoun1
          icount = icount + mbas1(isymd)*nmatij(isymk)
          icoun1 = icoun1 + mbas2(isymd)*nrhf(isymk)
        end do
        nrgij(isymdk) = icount
        nt1xao(isymdk) = icoun1
      end do
      do isymdk = 1, nsym
        icount = 0
        do isymk = 1, nsym
          isymd = muld2h(isymdk,isymk)
          it2xaos(isymd,isymk) = icount
          icount = icount + nt1ao(isymd)*nt1xao(isymk)
        end do
      end do
c
      do isymdp = 1, nsym
        isymj = isymdp
        isymd = isymdp
        isymakl = isymd
c
        koffs = kend1
        kghtf  = koffs + mbas1(isymd)*nbas2(isymdp)
        kgdpakl = kghtf + nrgij(isymakl)
        kend2   = kgdpakl + nrgij(isymakl)*nbas2(isymd) 
        lwrk2  = lwork - kend2
        if (lwrk2.lt.0) then
          call quit('Insufficient work space in ccr12prep')
        end if
        
c       get S_ddp
        kcmo = klamdhs + iglmrhs(isymd,isymj)
        koff1 = ksjbp + isjbp(isymdp,isymj)  
        ntotc  = max(mbas1(isymd),1)
        ntots  = max(nrhfsa(isymj),1)
    
        call dgemm('N','N',mbas1(isymd),nbas2(isymdp),
     &             nrhfsa(isymj),one,work(kcmo),ntotc,work(koff1),
     &             ntots,zero,work(koffs),ntotc)

        if (ianr12.eq.3) then
          kcmo  = klamdhs + iglmvis(isymd,isymj)
          koff1 = ksbbp   + isbbp(isymj,isymdp)  
          ntots = max(nvir(isymj),1)
          call dgemm('N','N',mbas1(isymd),nbas2(isymdp),nvir(isymj),
     &               one,work(kcmo),ntotc,work(koff1),ntots,
     &               one,work(koffs),ntotc)
        end if
c    
c       get I^dp_akl
        do idelp = 1, nbas2(isymdp)
          lauxd = .false.
          idxdp = idelp+ibas(isymdp)+basoff(isymdp)
          koff2 = kgdpakl-1 + nrgij(isymakl)*(idelp-1)+1
          call cc_r12getrint(work(koff2),idxdp,isymdp,
     &       nt1ao,it1ao,nt2aos,it2aos,nrgij,irgij,it1xao,it2xaos,
     &       nrhf,nmatij,imatij,
     &       ibasx,lauxd,.false.,fghtf,work(kend2),lwrk2)
          if (locdbg) then
            write(lupri,*)'norm^2 Idpakl:', 
     &      ddot(nrgij(isymd),work(koff2),1,work(koff2),1)
          end if
        end do  
c       add -I*S to I^dp_akl
        do idel = 1, mbas1(isymd)
          lauxd = .false.
          idxd = idel+ibas(isymd)
          call cc_r12getrint(work(kghtf),idxd,isymd,nt1ao,it1ao,
     &     nt2aos,it2aos,nrgij,irgij,it1xao,it2xaos,
     &     nrhf,nmatij,imatij,
     &     ibasx,lauxd,.false.,fghtf,work(kend2),lwrk2)
           if (locdbg) then
             write(lupri,*)'norm^2 Idakl:', 
     &        ddot(nrgij(isymd),work(kghtf),1,work(kghtf),1)
           end if
        
          do idelp = 1, nbas2(isymdp)
            idxddp = mbas1(isymd)*(idelp-1)+idel 
            koff2 = kgdpakl-1 + nrgij(isymakl)*(idelp-1)+1
            call daxpy(nrgij(isymakl),-work(koffs-1+idxddp),
     &               work(kghtf),1,work(koff2),1) 
          end do
        end do
c       save I^d'_akl on file
        do idelp = 1, nbas2(isymdp)
        koff2 = kgdpakl-1 + nrgij(isymakl)*(idelp-1)+1
        lauxd = .false.
        idxdp = idelp + ibas(isymdp) + basoff(isymdp)
c
        if (locdbg) then
          write(lupri,*)'mbas1(isymd)', mbas1(isymd)
          write(lupri,*) "save I^d'akl... isymdp,idelp,idxdp :",
     &                    isymdp,idelp,idxdp
          write(lupri,*)'norm^2 before write Idpakl:', 
     &    ddot(nrgij(isymd),work(koff2),1,work(koff2),1)
          write(lupri,*)"I^d'akl" 
          call output(work(koff2),1,nrgij(isymd),1,nrgij(isymd),
     &                nrgij(isymd),nrgij(isymd),1,lupri)
        end if
c
        call cc_r12putrint(work(koff2),idxdp,isymdp,nt1ao,it1ao,
     &     nt2aos,it2aos,nrgij,irgij,it1xao,it2xaos,
     &     nrhf,nmatij,imatij,
     &     ibasx,lauxd,fgdpakl,work(kend2),lwrk2)
        end do
      end do 
c
c----------------------------------------------------------------------
c     calculate \sum_mn part for V_alphaj,kl
c----------------------------------------------------------------------
c
      kvajkl  = kend1
      kend2 = kvajkl + nvajkl(1) 
      lwrk2  = lwork - kend2
      if (lwrk2.lt.0) then
        call quit('Insufficient work space in prep2')
      end if 
c
c     read V_alphajkl from file
      luvajkl = -1
      call gpopen(luvajkl,fvajkl,'unknown',' ','unformatted',
     &            idummy,.false.)
      rewind(luvajkl)
      read(luvajkl)(work(kvajkl+i-1),i=1,nvajkl(1))
      call gpclose(luvajkl,'KEEP')
  
      if (locdbg) then
        write(lupri,*)'Valphajkl after read in prep2:', 
     &              (work(kvajkl+i-1),i=1,nvajkl(1)) 
      end if
c     
      call cc_r12mkvaj2(work(kvajkl),1,work(klamdh),1,
     &                  work(klamdhs),1,work(kend2),lwrk2)
      if (locdbg) then
        write(lupri,*)'Valphajkl out of prep2:', 
     &                 (work(kvajkl-1+i),i=1,nvajkl(1))
      end if
c     write V_alphajkl on file
      luvajkl = -1
      call gpopen(luvajkL,fvajkl,'unknown',' ','unformatted',
     &            idummy,.false.)
      rewind(luvajkl)
      write(luvajkl)(work(kvajkl+i-1),i=1,nvajkl(1))
      call gpclose(luvajkl,'KEEP')

c     test if V_alphajkl is ok
      if (locdbg) then
        write(lupri,*)'norm^2 Valphajkl in prep2:', 
     &     ddot(nvajkl(1),work(kvajkl),1,work(kvajkl),1)
        call cc_r12mkvijkl(work(kvajkl),1,work(klamdh),1,
     &                     work(kend2),lwrk2,.false.,idummy,idummy)
        write(lupri,*)'Vajkl for ansatz 2 finished in prep2'
      end if

c----------------------------------------------------------------------
c     close file, clean trace and return:
c----------------------------------------------------------------------
      if(locdbg) write(lupri,*)'leaving ccr12prep2'
      call qexit('r12prep2')

      end 
*======================================================================*
      subroutine cc_r12mksjbp(smat,sjbp,sbbp,do_virt,cmo,
     &                        nbasx,nbas2,iglmrhs,iglmvis,
     &                        isjbp,isbbp,work,lwork)
*----------------------------------------------------------------------*
c     purpose: get S_abp symmetry blocked out of S packed as a 
c              triangular matrix and transform to S_Jbp and Sbbp
c     
c     H. Fliegl, C. Haettig, winter 2003
c     added transformation to Sbbp, Christof Haettig, spring 2005
c     adapted for CABS, Christian Neiss, winter 2006
*----------------------------------------------------------------------*
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>

      logical locdbg
      parameter (locdbg = .false.)

      logical do_virt
      integer lwork,isjbp(8,8),iglmrhs(8,8),iglmvis(8,8),isbbp(8,8)
      integer koff1, kcmo, ks, koffs, ntotb1, ntotbx, ntotrhs, ntotvir
      integer nbasx(8), nbas2(8), ndim, isym

#if defined (SYS_CRAY)
      real work(*), smat(*), sjbp(*), sbbp(*), cmo(*)
      real one, zero
#else
      double precision work(*), smat(*), sjbp(*), sbbp(*), cmo(*)
      double precision one, zero
#endif
      parameter(zero = 0.0d0, one = 1.0d0)

      call qenter('r12mksjbp')

      koff1   = 1 
      do isym = 1, nsym

        ndim = nbasx(isym)
        if (lwork.lt.ndim*ndim) then
          call quit('Insufficient work space in cc_r12mksjbp')
        end if

c       --------------------------
c       unpack triangular S matrix
c       --------------------------
        call sqmatr(ndim,smat(koff1),work)           
        koff1 = koff1 + ndim*(ndim+1)/2

c       --------------------------------------
c       transform S(alpha,beta') to S(J,beta')
c       --------------------------------------
        kcmo    = 1 + iglmrhs(isym,isym)
        ks      = 1 + isjbp(isym,isym) 
c       koffs   = 1 + nbasx(isym)*mbas1(isym)
        koffs   = 1 + nbasx(isym)*(nbasx(isym)-nbas2(isym))

        ntotb1  = max(mbas1(isym),1)
        ntotbx  = max(nbasx(isym),1)
        ntotrhs = max(nrhfsa(isym),1)

        call dgemm('T','N',nrhfsa(isym),nbas2(isym),mbas1(isym),
     &             one,cmo(kcmo),ntotb1,work(koffs),ntotbx,
     &             zero,sjbp(ks),ntotrhs) 

c       --------------------------------------
c       transform S(alpha,beta') to S(b,beta')
c       --------------------------------------
        if (do_virt) then
          kcmo   = 1 + iglmvis(isym,isym)
          ks     = 1 + isbbp(isym,isym)
          
          ntotvir = max(nvir(isym),1)

          call dgemm('T','N',nvir(isym),nbas2(isym),mbas1(isym),
     &               one,cmo(kcmo),ntotb1,work(koffs),ntotbx,
     &               zero,sbbp(ks),ntotvir) 
        end if

      end do  

      call qexit('r12mksjbp')
      end 
*=====================================================================*
      subroutine cc_r12mksmupbp(cmo,smupbp,iaodis2,work,lwork)
c--------------------------------------------------------------------- 
c     purpose: make D_mu'b' = \sum_p' C_b'p' * C_mu'p' for ansatz 2
c              for CABS: sum over MO and aux. basis     
c
c     H. Fliegl, C. Haettig, winter 2003
c     C. Neiss, winter 2006: adapted for CABS
c--------------------------------------------------------------------- 
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
      logical locdbg
      parameter (locdbg = .false.)

      integer lwork, isymp, isymq, isympq, isymbp, isym, 
     &        isymmup, kcmo, ntotc, ntots, 
     &        iaodis2(8,8), ks, isympbp 
      integer idxbpp, idxbpp2, ip, ib, koff0, idxp, idxb
      integer nbasx(8),norbx(8),basoff(8),orboff(8)

#if defined (SYS_CRAY) 
      real cmo(*), smupbp(*), work(*)       
#else
      double precision cmo(*), smupbp(*), work(*)       
#endif
      double precision one, zero
      parameter(one = 1.0d0, zero = 0.0d0)   
 
      call qenter('r12mksmupbp')      
c
      if (r12cbs) then
        do isym = 1, nsym
          nbasx(isym) = mbas1(isym)+mbas2(isym)
          norbx(isym) = norb1(isym)+norb2(isym)
          basoff(isym) = 0
          orboff(isym) = 0
        end do
      else
        call icopy(8,mbas2,1,nbasx,1)
        call icopy(8,norb2,1,norbx,1)
        call icopy(8,mbas1,1,basoff,1)
        call icopy(8,norb1,1,orboff,1)
      end if
c
      koff0 = 1
      do isymp = 1, nsym
        isymbp = isymp
        isymmup = isymp
c       isympbp = muld2h(isymp,isymbp) ! eq 1
c       get C_p'b' out of CMO
        if (locdbg) then
          write(lupri,*)'mbas1(isymbp),mbas2(isymbp)',
     &                   mbas1(isymbp),mbas2(isymbp)
          write(lupri,*)'norb1(isymp),norb2(isymp)',
     &                   norb1(isymp),norb2(isymp)
          write(lupri,*)'norb2(isymp)*mbas2(isymbp)',
     &                   norb2(isymp)*mbas2(isymbp)
          write(lupri,*)'nrhfsb(isymp)',nrhfsb(isymp)
          write(lupri,*)'isymbp, isymp', isymbp, isymp
        end if
        do ip=orboff(isymp)+1+nrhfsb(isymp),
     &       (norb1(isymp)+norb2(isymp)+nrhfsb(isymp))
          idxp = ip - nrhfsb(isymp) - orboff(isymp)
          do ib = basoff(isymbp)+1,(mbas1(isymbp)+mbas2(isymbp))
            idxb   = ib - basoff(isymp)
            idxbpp = koff0-1 + (mbas1(isymbp)+mbas2(isymbp))*(ip-1)+ib 
            idxbpp2 = nbasx(isymbp)*(idxp-1)+idxb 
            work(idxbpp2) = cmo(idxbpp)
            if (locdbg) then
              write(lupri,*)'wrk,idxbpp,idxbpp2',idxbpp, idxbpp2
              write(lupri,*)work(idxbpp2)
            end if
          end do
        end do
        koff0 = koff0 + (mbas1(isymbp)+mbas2(isymbp)) *
     &            (nrhfsb(isymp)+norb1(isymp)+norb2(isymp))
        
c       transform to D
        
        ks    = 1+iaodis2(isymbp,isymmup) 
        ntotc = max(nbasx(isymbp),1)

        call dgemm('N','T',nbasx(isymbp),nbasx(isymmup),norbx(isymp),
     &             one,work,ntotc,work,ntotc,zero,
     &             smupbp(ks),ntotc)
        if (locdbg) then
          write(lupri,*)'resorted CMO in cc_r12mksmupbp'
          call output(work,1,nbasx(isymbp),1,
     &         norbx(isymp),nbasx(isymbp),norbx(isymp),1,lupri)
          write(lupri,*)'S_bp,mup in cc_r12mksmupbp'
          call output(smupbp(ks),1,nbasx(isymbp),1,
     &         nbasx(isymmup),nbasx(isymbp),nbasx(isymmup),1,lupri)
        end if
      end do

      call qexit('r12mksmupbp')      
      end
*=====================================================================*
      subroutine cc_r12getkiajb(xkipjq,xkiajb,oneaux,nt2amx,
     &                          nh2amx,ih2am,nh1am,ih1am,norb1,
     &                          nrhf1,nrhfs1)
c---------------------------------------------------------------------
c     purpose: get K^ij_ab out of K^ij_pq stored as a lower 
c              symmetry packed triangular matrix
c
c     H. Fliegl, C. Haettig, spring 2004
c---------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
      logical locdbg,oneaux
      parameter (locdbg = .false.)
      integer nh2amx, ih2am(8,8), nh1am(8), ih1am(8,8), nvircc(8),
     &        nt2amx, nt1am(8), it1am(8,8), it2am(8,8), isympi,
     &        isymi, isym, icount1, icount2, isymj, isymk, ldr12(8),
     &        isymqj, isymq, isymp, nqj, npi, nbj, nai, idxpiqj, 
     &        idxaibj, index, i,j,p,a,b,q, norb1(8),nrhf1(8),nrhfs1(8) 
#if defined (SYS_CRAY)
      real xkiajb(*), xkipjq(*)
#else
      double precision xkiajb(*), xkipjq(*), ddot
#endif
      index(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j

      call qenter('getkiajb')

c     calculate some offsets and dimensions
c     calculate number of active virtual orbitals in CC
c     nvirfr = frozen virtual orbitals        
      do isym = 1, nsym
        nvircc(isym) = norb1(isym) - nrhfsa(isym) - nvirfr(isym)
      end do

      if (locdbg) then
        write(lupri,* ) 'isym nrhfs1 nrhf1 nvircc norb1'
        do isym = 1, nsym
          write(lupri,*) isym,nrhfs1(isym),nrhf1(isym),
     &        nvircc(isym),norb1(isym)
        end do 
      end if

      nt2amx = 0
      do isympi = 1, nsym
        nt1am(isympi) = 0
        do isymi = 1, nsym
          isymp = muld2h(isympi,isymi)
          nt1am(isympi) = nt1am(isympi) + nvircc(isymp)*nrhf1(isymi)
        end do
        nt2amx = nt2amx + nt1am(isympi)*(nt1am(isympi)+1)/2
      end do
     
      do isymk = 1, nsym
        icount1 = 0
        icount2 = 0
        do isymj = 1, nsym
          isymi = muld2h(isymj,isymk)  
          it1am(isymi,isymj) = icount1
          icount1 = icount1 + nrhf1(isymj)*nvircc(isymi)
          if (isymj .gt. isymi) then
             it2am(isymi,isymj) = icount2
             it2am(isymj,isymi) = icount2
             icount2 = icount2 + nt1am(isymi)*nt1am(isymj)
          else if (isymk .eq. 1) then
             it2am(isymi,isymj) = icount2
             icount2 = icount2 + nt1am(isymi)*(nt1am(isymi)+1)/2
          end if
        end do 
      end do

      do isym = 1, nsym
        if (oneaux) then
          ldr12(isym) = norb1(isym)
        else
          ldr12(isym) = nvir(isym)
        end if
      end do
c     
      if (locdbg) then 
        write(lupri,*) 'norm^2(kpiqj):',ddot(nh2amx,xkipjq,1,xkipjq,1)
      end if
c
c    ---------------------------------------------------------------
c     loop over symmetries and indices pick out p,q as virtual MO's 
c     and reorder integrals
c    
c     reorder integrals
c     npi    : pair p,i with i occupied (r12 orbital) and p MO 
c              virtual dimensions as needed for
c              MP2-R12 code
c     nai    : pair a,i dimensions as needed for
c              CC code
c     idxpiqj: quadrupel pi,qj dimensions as 
c              needed for MP2-R12 code
c     idxaibj: quadrupel ai,bj dimensions as 
c              needed for CC code         
c    ---------------------------------------------------------------
      do isympi = 1, nsym
         isymqj = isympi
         do isymp = 1, nsym
          isymi = muld2h(isympi,isymp)
          do isymq =1, nsym
            isymj = muld2h(isymqj,isymq)
c
              do i = 1, nrhf1(isymi)
                do p = nrhfs1(isymp)+1, nrhfs1(isymp)+nvircc(isymp)
                  a = p - nrhfs1(isymp)
                  do j = 1, nrhf1(isymj)
                    do q = nrhfs1(isymq)+1, nrhfs1(isymq)+
     &                                         nvircc(isymq)
                      b = q - nrhfs1(isymq)
                      nqj = ih1am(isymq,isymj)+ ldr12(isymq)*(j-1)+q 
                      npi = ih1am(isymp,isymi)+ ldr12(isymp)*(i-1)+p 
                      nbj = it1am(isymq,isymj)+nvircc(isymq)*(j-1)+b
                      nai = it1am(isymp,isymi)+nvircc(isymp)*(i-1)+a
                      idxpiqj = ih2am(isympi,isymqj)+ index(npi,nqj) 
                      idxaibj = it2am(isympi,isymqj)+ index(nai,nbj) 
                      xkiajb(idxaibj) = xkipjq(idxpiqj)
                    end do
                  end do
                end do
              end do
c
          end do
        end do
      end do

      if (locdbg) then
        write(lupri,*)'xkiajb'
        do isymi =1, nsym
          call outpak(xkiajb(1+it2am(isymi,isymi)),nt1am(isymi),1,lupri)
        end do
        write(lupri,*) 'norm^2(kiajb):',ddot(nt2amx,xkiajb,1,xkiajb,1)
      end if
c
      call qexit('getkiajb')
      end
*======================================================================*
      subroutine ccrhs_epp(omega,t2am,isymt2am,work,lwork,aproxr12,lres,
     &                     lufc2,fc2am,ifile)
*----------------------------------------------------------------------*
c     purpose: add E'' term to Omega^ij_kl (stored in omega) for the
c              ansaetze 2 and 3
c
c     omega^ij_kl <-- omega^ij_kl 
c                       + sum_ab (K^ab_kl - t^ab_kl) * T^ab_ij
c
c     the contrib. of K is skipped for approximation A
c
c     for the approximations ??? is in ansatz 2 added:
c
c     omega^ij_kl <-- omega^ij_kl 
c                       + sum_ab (e_k+e_l-e_i-e_j) r^ab_kl * T^ab_ij
c
c     while in ansatz 3 the following term is added:
c 
c     omega^ij_kl <-- omega^ij_kl 
c                       + sum_ab (e_k+e_l-e_a-e_b) r^ab_ij * T^ab_ij
c
c
c     H. Fliegl, C. Haettig, spring 2004
c     adapted for ansatz 3, Christof Haettig, spring 2005
*----------------------------------------------------------------------*
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
#include <ccr12int.h>
      logical locdbg,lres
      parameter (locdbg = .false.)
      integer lwork,lunit,isymt2am,isymkint,kt2pck,kpck,lwrk1,idum,
     &        kend1,iopt,idummy,kap,lufc2,ifile,kend2,lwrk2
      character*(*) fc2am
#if defined (SYS_CRAY)
      real work(*), t2am(*), omega(*),factor,ddot,two,one
#else
      double precision work(*), t2am(*), omega(*),factor,ddot,two,one
#endif
      parameter(two=2.0D0, one=1.0D0)
      character*(*) aproxr12

      call qenter('rhs_epp')

      kt2pck = 1
      kpck  = kt2pck + max(nt2am(isymt2am),nt2r12(isymt2am),nt2r12(1)) ! if lres isymt2am = isymtr!
      kend1 = kpck   + nt2r12(1) 
      lwrk1 = lwork - kend1

      kend2 = kend1
      lwrk2 = lwork - kend2

      if (lwrk1.lt.0 .or. lwrk2.lt.0) then
        call quit('Insufficient work space in ccrhs_epp')
      end if

c     ---------------------------------------------
c     read t^ab_kl integrals from file
c     ---------------------------------------------
      lunit = -1
      call gpopen(lunit,fkr12,'unknown',' ','unformatted',idum,.false.)
      read(lunit)(work(kpck-1+i),i=1,nt2r12(1))  
      call gpclose(lunit,'KEEP')

c     ----------------------------------------------------
c     for approximation B substract K^ab_kl from t^ab_kl
c     ----------------------------------------------------
      if (aproxr12(1:1).eq.'B') then
        lunit = -1
        call gpopen(lunit,ftr12,'unknown',' ','unformatted',
     &              idum,.false.)
        read(lunit)(work(kt2pck+i-1),i=1,nt2r12(1))  
        call gpclose(lunit,'KEEP')
        
        ! substract (t - K) => work(kpck)
        call daxpy(nt2r12(1),-1.0d0,work(kt2pck),1,work(kpck),1)
      end if

c     ------------------------------------------------------------
c     in ansatz 3 add for approximations A' and B the term
c        (eps_k + eps_l - eps_a - eps_b) * r^ab_kl
c     ------------------------------------------------------------
      if (ianr12.eq.3 .and. 
     &   (aproxr12(1:2).eq."A'".or.aproxr12(1:1).eq."B") ) then
        lunit = -1
        call gpopen(lunit,fr12r12,'unknown',' ','unformatted',
     &              idum,.false.)
        read(lunit)(work(kt2pck+i-1),i=1,nt2r12(1))  
        call gpclose(lunit,'KEEP')
        call cc_t2xe(work(kt2pck),work(kend2),lwrk2)
        call daxpy(nt2r12(1),+1.0d0,work(kt2pck),1,work(kpck),1)
      end if

c     ----------------------------------------------------------
c     get doubles amplitudes R^mn_ab stored as triangular matrix
c     ----------------------------------------------------------
      if (lres) then
        call cc_rvec(lufc2,fc2am,nt2am(isymt2am),nt2am(isymt2am),
     &               ifile,work(kt2pck))
        call cclr_diascl(work(kt2pck),two,isymt2am)
      else
        iopt = 0
        call cc_t2pk(work(kt2pck),t2am,isymt2am,iopt)      
      end if

c     ------------------------------------------------------
c     calculate Omega^ij_kl = \sum_ab k^ab_kl * t^ij_ab
c     ------------------------------------------------------
      if (locdbg) then
        write(lupri,*) 'in ccrhs_epp before cc_r12mi2:'
        write(lupri,*) 'work(kt2pck):'
        call outpak(work(kt2pck),nt1am(1),1,lupri)
      end if
      isymkint = 1
      factor = 1.0d0
      call cc_r12mi2(omega,work(kpck),work(kt2pck),isymkint,isymt2am,
     &               factor,work(kend1),lwrk1)

      if (locdbg) then
        write(lupri,*)'Norm^2 Omega^ij_kl in ccrhs_epp', 
     &                 ddot(ntr12sq(isymt2am),omega,1,omega,1)
      end if

c     ------------------------------------------------------------
c     in ansatz 2 add for approximation A' and B the term
c      omega^ij_kl += sum_ab (e_k+e_l-e_i-e_j) r^ab_kl * c^kl_ij
c     ------------------------------------------------------------
      if (ianr12.eq.2 .and.
     &   (aproxr12(1:2).eq."A'".or.aproxr12(1:1).eq."B") ) then
        kap   = kend1
        kend1 = kap + ntr12sq(isymt2am)         
        lwrk1 = lwork - kend1
        if (lwrk1.lt.0) then
          call quit('Insufficient work space in ccrhs_epp')
        end if

        lunit = -1
        call gpopen(lunit,fr12r12,'unknown',' ','unformatted',
     &              idum,.false.)
        read(lunit)(work(kpck-1+i),i=1,nt2r12(1))  
        call gpclose(lunit,'KEEP')

        ! calculate \sum_ab r^ab_kl*t^ij_ab
        factor = one
        isymkint = 1
        call dzero(work(kap),ntr12sq(isymt2am))
        call cc_r12mi2(work(kap),work(kpck),work(kt2pck),isymkint,
     &                 isymt2am,factor,work(kend1),lwrk1)

        call cc_r12mkediff(work(kap),isymt2am,work(kend1),lwrk1)

        call daxpy(ntr12sq(isymt2am),1.0d0,work(kap),1,omega,1)
      end if 

c     -----------------------------------------------------------------
c     that's it; restore T2 amplitudues, print debug output and return:
c     -----------------------------------------------------------------
      if (.not.lres) then
c       pack T2 amplitudes back to a quadratic matrix
        call cc_t2sq(work(kt2pck),t2am,isymt2am)
      end if

      if (locdbg) then
        write(lupri,*) 'omega^ij_kl after ccrhs_epp contribution:'
        write(lupri,*) 'norm^2 omega^ij,kl',
     &                 ddot(ntr12sq(isymt2am),omega,1,omega,1)
      end if
 
      call qexit('rhs_epp')
      end
*=====================================================================*
      subroutine cc_r12mkediff(r12int,isymr,work,lwork)
c---------------------------------------------------------------------
c     purpose: calculate r12int*(e_k+e_l-e_i-e_j)
c
c     H. Fliegl, C. Haettig, spring 2004
c     modified C. Neiss, september 2005
c---------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <inftap.h>
      integer lwork, kscr1,kend1,lwrk1,isymkl,isymij,keij,kekl,
     &        idxijkl,isymr,ij,kl,idummy,norbtsx
      integer nij,nkl,inmatij(8),inmatkl(8),isym
#if defined (SYS_CRAY)
      real work(*), r12int(*)
#else
      double precision work(*), r12int(*)
#endif
c
      call qenter('r12mkediff')
c 
c-------------------------------------
c     Read canonical orbital energies.
c-------------------------------------
c
      call gpopen(lusifc,'SIRIFC','UNKNOWN',' ','UNFORMATTED',idummy,
     &            .false.)
      rewind lusifc
c
      call mollab('FULLBAS ',lusifc,lupri)
      read (lusifc) idummy,norbtsx
      kscr1 = 1
      kend1 = kscr1 + norbtsx
      lwrk1 = lwork - kend1
      if (lwrk1.lt.0) then
        call quit('Insufficient work space in cc_r12mkediff')
      end if
      read (lusifc) (work(kscr1+i-1), i=1,norbtsx)
      call gpclose(lusifc,'KEEP')
c
      nkl = 0
      nij = 0
      do isym = 1, nsym
        inmatkl(isym) = nkl
        inmatij(isym) = nij
        nkl = nkl + nmatkl(isym)
        nij = nij + nmatij(isym)
      end do
c
      keij = kend1
      kekl = keij + nij
      kend1 = kekl + nkl
      lwrk1 = lwork - kend1
      if (lwrk1.lt.0) then
        call quit('Insufficient work space in cc_r12mkediff')
      end if
c
      call cc_r12epair(work(kekl),nkl,work(kscr1),inmatkl,imatkl,nrhfb)
      call cc_r12epair(work(keij),nij,work(kscr1),inmatij,imatij,nrhf)
c
      do isymkl = 1, nsym
        isymij = muld2h(isymr,isymkl)
        do kl = 1, nmatkl(isymkl)
          do ij = 1, nmatij(isymij)
            idxijkl = itr12sqt(isymij,isymkl)+nmatij(isymij)*(kl-1)+ij
            r12int(idxijkl) = r12int(idxijkl)*
     &      (work(kekl-1+kl+inmatkl(isymkl))-
     &       work(keij-1+ij+inmatij(isymij))) 
          end do
        end do
      end do
      
      call qexit('r12mkediff')
      end 
*=====================================================================*
      subroutine cc_r12mkediff2(c2am,isymc,work,lwork)
c---------------------------------------------------------------------
c     purpose: calculate c12*(e_k+e_l-e_i-e_j)
c
c     H. Fliegl, C. Haettig, spring 2004
c---------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <inftap.h>
       integer lwork, kscr1,kend1,lwrk1, isymc,index,idummy,
     &         kend2,lwrk2,isymk,isymi,isyml,isymj,
     &         mk,mi,mj,ml
       integer isymki,isymlj,keki,kelj,lj,ki,idxkilj
#if defined (SYS_CRAY)
       real work(*), c2am(*),el,ei,ek,ej
#else
       double precision work(*), c2am(*),el,ei,ek,ej
#endif
      index(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j

      call qenter('r12mkediff2')

      kscr1 = 1
      kend1 = kscr1 + norbts
      lwrk1 = lwork - kend1
      if (lwrk1.lt.0) then
        call quit('Insufficient work space in cc_r12mkediff2')
      end if
c-------------------------------------
c     Read canonical orbital energies.
c-------------------------------------
c
      call gpopen(lusifc,'SIRIFC','UNKNOWN',' ','UNFORMATTED',idummy,
     &            .false.)
      rewind lusifc
c
      call mollab('TRCCINT ',lusifc,lupri)
      read (lusifc)
      read (lusifc) (work(kscr1+i-1), i=1,norbts)
      call gpclose(lusifc,'KEEP')
c
      if (froimp .or. froexp)
     *   call ccsd_delfro(work(kscr1),work(kend1),lwrk1)
c
      do isymki = 1, nsym
        isymlj = muld2h(isymc,isymki)
        if (isymki.le.isymlj) then
         keki = kend1
         kelj = keki + nmatki(isymki)
         kend2 = kelj + nmatki(isymlj)
         lwrk2 = lwork - kend2
         if (lwrk2.lt.0) then
           call quit('Insufficient work space in cc_r12mkediff2')
         end if
         do isymk = 1, nsym
           isymi = muld2h(isymki,isymk)
          do isyml = 1, nsym
           isymj = muld2h(isymlj,isyml)
           do k = 1, nrhfb(isymk)
             mk = iorb(isymk)+k
             ek = work(kscr1+mk-1)
             do i = 1, nrhfa(isymi)
               ki = imatki(isymk,isymi)+nrhfb(isymk)*(i-1)+k
               mi = iorb(isymi)+i
               ei = work(kscr1+mi-1)
               work(keki+ki-1) = ek - ei
             end do
           end do
           do l = 1, nrhfb(isyml)
             ml = iorb(isyml) +l
             el = work(kscr1+ml-1)
             do j = 1, nrhfa(isymj)
               lj = imatki(isyml,isymj)+nrhfb(isyml)*(j-1)+l
               mj = iorb(isymj)+j
               ej = work(kscr1+mj-1)
               work(kelj+lj-1) = el - ej
             end do
           end do
          end do
         end do
c
         if (isymki.eq.isymlj) then
           do lj = 1, nmatki(isymlj)
             do ki = 1, lj
               idxkilj = itr12am(isymki,isymlj)+index(ki,lj) 
               c2am(idxkilj) = c2am(idxkilj) *  
     &          ( work(keki-1+ki) + work(kelj-1+lj) )
             end do
           end do
         else
           do lj = 1, nmatki(isymlj)
             do ki = 1, nmatki(isymki)
               idxkilj = itr12am(isymki,isymlj)+nmatki(isymki)*(lj-1)+
     &                    ki
               c2am(idxkilj) = c2am(idxkilj)*
     &          ( work(keki-1+ki) + work(kelj-1+lj) )
             end do
           end do
         end if
        end if
      end do

      call qexit('r12mkediff2')
      return
      end 
*======================================================================*
      subroutine cc_t2xe(t2am,work,lwork)
*----------------------------------------------------------------------*
*     Purpose: multiply vector on t2am in place with orbital 
*              energy differences
*      
*              t2(ak,bl) <-- t2(ak,bl) * (e_k + e_l - e_a - e_b)
*
*              amplitudes are addressed as triangle
*
*     C. Haettig, spring 2005
*     modified for general r12 orbitals: C. Neiss, september 2005
*----------------------------------------------------------------------*
      implicit none
#include <priunit.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>

      logical locdbg, delfro
      parameter (locdbg = .false.)
      parameter (delfro = .true.)

      integer isymbj, isymai, isymj, isymi, isymb, isyma, index,
     &        koffa, koffb, nbj, nai, naibj, lwork, idum
      integer nkl,inmatkl(8),iak(8,8),icount,kekl,kend1,lwrk1,lunit,
     &        koffij,isymij,keps,norbtsx
#if defined (SYS_CRAY)
      real t2am(*), work(lwork)
#else
      double precision t2am(*), work(lwork)
#endif

      index(i,j) = max(i,j)*(max(i,j) - 3)/2 + i + j

      call qenter('cc_t2xe')

      nkl = 0
      do isymij = 1, nsym
        inmatkl(isymij) = nkl
        nkl = nkl + nmatkl(isymij)
      end do
c
      do isymai = 1, nsym
        icount = 0
        do isymi = 1, nsym
          isyma = muld2h(isymai,isymi)
          iak(isyma,isymi) = icount
          icount = icount + nvir(isyma)*nrhfb(isymi)
        end do
      end do
c
      kekl  = 1
      kend1 = kekl + nkl
      lwrk1 = lwork - kend1
      if (lwrk1.lt.0) call quit('Insufficient memory in CC_T2XE')
c
      lunit = -1
      call gpopen(lunit,'SIRIFC','OLD',' ','UNFORMATTED',idum,
     &            .false.)
      rewind lunit
      call mollab('FULLBAS ',lunit,lupri)
      read (lunit) idum,norbtsx
      keps  = kend1
      kend1 = keps + norbtsx
      lwrk1 = lwork - kend1
      if (lwrk1.lt.0) call quit('Insufficient memory in CC_T2XE')
      read (lunit) (work(keps+i-1), i=1,norbtsx)
c     call gpclose(lunit,'KEEP')
c 
      call cc_r12epair(work(kekl),nkl,work(keps),inmatkl,imatkl,nrhfb)
c
      rewind lunit
      call mollab('TRCCINT ',lunit,lupri)
      read (lunit)
      read (lunit) (work(keps+i-1), i=1,norbts)
      call gpclose(lunit,'KEEP')
c
      if (froimp .or. froexp)
     *  call ccsd_delfro(work(keps),work(kend1),lwrk1)
c
      call fock_reorder(work(keps),work(kend1),lwrk1)
c
      do isymbj = 1,nsym
        isymai = isymbj
        do isymj = 1,nsym
          isymb = muld2h(isymj,isymbj)
          do isymi = 1,nsym
            isyma = muld2h(isymi,isymai)
            do j = 1,nrhfb(isymj)
              do b = 1,nvir(isymb)
                koffb = ivir(isymb) + b
                nbj = iak(isymb,isymj)+nvir(isymb)*(j-1)+b
                do i = 1,nrhfb(isymi)
                  do a = 1,nvir(isyma)
                    koffa = ivir(isyma) + a
                    nai = iak(isyma,isymi)+nvir(isyma)*(i-1)+a

                    if (nai .le. nbj) then
                      isymij = muld2h(isymi,isymj)
                      koffij = inmatkl(isymij)+imatkl(isymi,isymj)+
     &                         nrhfb(isymi)*(j-1)+i-1
                      naibj = it2r12(isymai,isymbj) + index(nai,nbj)
                      t2am(naibj) = t2am(naibj) *
     &                   (work(kekl+koffij)-work(keps+koffa-1)-
     &                                      work(keps+koffb-1))
                    end if

                  end do
                end do
              end do
            end do
          end do
        end do
      end do

      call qexit('cc_t2xe')
      return
      end
*=====================================================================*
      subroutine cc_r12mi2(omega,xkint,t2am,isymkint,isymt2am,
     &                     factor,work,lwork)
c---------------------------------------------------------------------
c     purpose: calculate Omega^ij_kl = \sum_ab k^ab_kl * t^ij_ab
c              input: xkint and t2am packed as a triangular matrix
c
c     H. Fliegl, C. Haettig, spring 2004 
c---------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
      logical locdbg
      parameter (locdbg = .false.)
      integer isymkint, isymt2am, lwork, ntoto, ntotk,
     &        isyma, isymij, isymkl, kxint, kt2am, kend1, lwrk1,
     &        isymab, isymint, isymb, kom
#if defined (SYS_CRAY)
      real work(*), xkint(*),t2am(*), omega(*),factor
      real one,ddot,xk,xt
#else
      double precision work(*), xkint(*),t2am(*), omega(*),factor
      double precision one,ddot,xk,xt
#endif
      parameter (one = 1.0d0)

      call qenter('r12mi2')
c
      isymint = muld2h(isymkint,isymt2am) !isymint = isymomg
       
      if (locdbg) then
        xk = 0.0d0
        xt = 0.0d0 
      end if
      do isymb = 1, nsym
        do isyma = 1, nsym
          isymab = muld2h(isyma,isymb)
          isymkl = muld2h(isymkint,isymab)
          isymij = muld2h(isymkl,isymint)
c
          kxint = 1
          kt2am = kxint + nvir(isyma)*nmatkl(isymkl)
          kend1 = kt2am + nvir(isyma)*nmatij(isymij)
          lwrk1 = lwork - kend1
          if (lwrk1.lt.0) then
            call quit('Insufficient work space in cc_r12mi2')
          end if  
c
          do b = 1, nvir(isymb)

C           call cc_r12sortk(work(kxint),xkint,b,isymkint,isymb,
C    &                       isyma,nrhfb,.false.)
            call cc_r12sort(xkint,isymkint,work(kxint),b,isymb,
     &                        1,nvir(isyma),isyma,nvir,nrhfb)
            if (locdbg) then
              xk = xk +
     &        ddot(nvir(isyma)*nmatkl(isymkl),work(kxint),1,
     &             work(kxint),1)
            end if

C           call cc_r12sortk(work(kt2am),t2am,b,isymt2am,isymb,
C    &                       isyma,nrhf,.false.)
            call cc_r12sort(t2am,isymt2am,work(kt2am),b,isymb,
     &                        1,nvir(isyma),isyma,nvir,nrhf)
            if (locdbg) then
              xt = xt +
     &        ddot(nvir(isyma)*nmatij(isymij),work(kt2am),1,
     &             work(kt2am),1)
            end if
c
            ntotk = max(nvir(isyma),1)
            ntoto = max(nmatij(isymij),1)
c
            kom = itr12sqt(isymij,isymkl) +1 
            call dgemm('T','N',nmatij(isymij),nmatkl(isymkl),
     &               nvir(isyma),factor,work(kt2am),ntotk,work(kxint),
     &               ntotk,one,omega(kom),ntoto)
          end do
        end do
      end do

      if (locdbg) then
        write(lupri,*)'xk =', xk
        write(lupri,*)'xt =', xt
      end if
 
      call qexit('r12mi2')
      end
*=====================================================================*
      subroutine cc_r12sortk(xk,xkint,b,isymkint,isymb,isyma,nr12orb,
     &                       lopt)
c---------------------------------------------------------------------
c     purpose: sort k^ab_kl as three index array with fixed b
c              expect that k is stored as a lower triangular 
c              matrix
c              lopt = .false. : xkint stored with it2am/it1am
c              lopt = .true.  : xkint stored with ih2am/ih1am
c              nr12orb = # of k (nrhf or nrhfb)
c              
c     H. Fliegl, C. Haettig, spring 2004
c     modified C. Neiss august 2005
c---------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
      integer isymkint, index, isymakl, isyma, isyml, 
     &        isymk, isymak, isymbl, isymkl, idxbl, idxkl, idxakl,
     &        idxabkl, idxak, isymb
      integer isym,isym1,isym2,icoun1,icoun2,iak(8,8),nak(8),iakbl(8,8),
     &        ipk(8,8),npk(8),ipkql(8,8),nr12orb(8),nkl(8),ikl(8,8),
     &        icoun3,nakbl(8)
#if defined (SYS_CRAY)
      real xkint(*), xk(*), ddot
#else
      double precision xkint(*), xk(*), ddot
#endif
      logical lopt,locdbg
      parameter (locdbg = .false.)

      index(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j
      call qenter('r12sortk')
c
      do isym = 1, nsym
        icoun1 = 0
        icoun2 = 0
        icoun3 = 0
        do isym2 = 1, nsym
          isym1 = muld2h(isym,isym2)
          iak(isym1,isym2) = icoun1
          ipk(isym1,isym2) = icoun2
          ikl(isym1,isym2) = icoun3
          icoun1 = icoun1 + nvir(isym1)*nr12orb(isym2)
          icoun2 = icoun2 + norb1(isym1)*nr12orb(isym2)
          icoun3 = icoun3 + nr12orb(isym1)*nr12orb(isym2)
        end do
        nak(isym) = icoun1
        npk(isym) = icoun2
        nkl(isym) = icoun3
      end do
      do isym = 1, nsym
        icoun1 = 0
        icoun2 = 0
        do isym2 = 1, nsym
          isym1 = muld2h(isym,isym2)
          if (isym2.gt.isym1) then
            iakbl(isym1,isym2) = icoun1
            iakbl(isym2,isym1) = icoun1
            icoun1 = icoun1 + nak(isym1)*nak(isym2)
            ipkql(isym1,isym2) = icoun2
            ipkql(isym2,isym1) = icoun2
            icoun2 = icoun2 + npk(isym1)*npk(isym2)
          else if (isym2.eq.isym1) then
            iakbl(isym1,isym2) = icoun1
            icoun1 = icoun1 + nak(isym1)*(nak(isym2)+1)/2
            ipkql(isym1,isym2) = icoun2
            icoun2 = icoun2 + npk(isym1)*(npk(isym2)+1)/2
          end if
        end do
        nakbl(isym) = icoun1
      end do
c
      if (locdbg) then
        write(lupri,*) 'Entered CC_R12SORTK'
        write(lupri,*)'norm^2(xkint): ',
     &                 ddot(nakbl(isymkint),xkint,1,xkint,1)
C       if (nsym.eq.1) 
C    &    call outpak(xkint,nak(1),1,lupri)
      end if     
c     
      isymakl = muld2h(isymb,isymkint)      
      isymkl = muld2h(isyma,isymakl)
      do isymk = 1, nsym
        isyml = muld2h(isymkl,isymk)
        isymak = muld2h(isyma,isymk)
        isymbl = muld2h(isymb,isyml)
        if (.not.lopt) then
          do l = 1, nr12orb(isyml)
            idxbl = iak(isymb,isyml)+nvir(isymb)*(l-1)+b
            do k = 1, nr12orb(isymk)  
              idxkl  = ikl(isymk,isyml)+nr12orb(isymk)*(l-1)+k
              if (isymak.eq.isymbl) then
                do a = 1, nvir(isyma)  
                 idxakl = nvir(isyma)*(idxkl-1)+a 
                 idxak = iak(isyma,isymk)+nvir(isyma)*(k-1)+a
                 idxabkl = iakbl(isymak,isymbl)+index(idxak,idxbl)
                 xk(idxakl) = xkint(idxabkl)
                end do
              else if (isymak.lt.isymbl) then  
               idxakl = nvir(isyma)*(idxkl-1)+1 
               idxak = iak(isyma,isymk)+nvir(isyma)*(k-1)+1
               idxabkl = iakbl(isymak,isymbl)+
     &                   nak(isymak)*(idxbl-1)+idxak
               call dcopy(nvir(isyma),xkint(idxabkl),1,
     &                    xk(idxakl),1)
              else if (isymak.gt.isymbl) then  
                do a = 1, nvir(isyma)  
                 idxakl = nvir(isyma)*(idxkl-1)+a 
                 idxak = iak(isyma,isymk)+nvir(isyma)*(k-1)+a
                 idxabkl = iakbl(isymbl,isymak)+
     &                     nak(isymbl)*(idxak-1)+idxbl
                 xk(idxakl) = xkint(idxabkl)
                end do
              end if
            end do
          end do
        else 
          do l = 1, nr12orb(isyml)
            idxbl=ipk(isymb,isyml)+norb1(isymb)*(l-1)+nrhfsa(isymb)+b
            do k = 1, nr12orb(isymk)  
              idxkl  = ikl(isymk,isyml)+nr12orb(isymk)*(l-1)+k
              if (isymak.eq.isymbl) then
                do a = 1, nvir(isyma)  
                 idxakl = nvir(isyma)*(idxkl-1)+a 
                 idxak = ipk(isyma,isymk) +
     &                      norb1(isyma)*(k-1)+nrhfsa(isyma)+a
                 idxabkl = ipkql(isymak,isymbl)+index(idxak,idxbl)
                 xk(idxakl) = xkint(idxabkl)
                end do
              else if (isymak.lt.isymbl) then  
               idxakl = nvir(isyma)*(idxkl-1)+1 
               idxak = ipk(isyma,isymk) +
     &                    norb1(isyma)*(k-1)+nrhfsa(isyma)+1
               idxabkl = ipkql(isymak,isymbl)+
     &                   npk(isymak)*(idxbl-1)+idxak
               call dcopy(nvir(isyma),xkint(idxabkl),1,
     &                    xk(idxakl),1)
              else if (isymak.gt.isymbl) then  
                do a = 1, nvir(isyma)  
                 idxakl = nvir(isyma)*(idxkl-1)+a 
                 idxak = ipk(isyma,isymk) +
     &                     norb1(isyma)*(k-1)+nrhfsa(isyma)+a
                 idxabkl = ipkql(isymbl,isymak)+
     &                     npk(isymbl)*(idxak-1)+idxbl
                 xk(idxakl) = xkint(idxabkl)
                end do
              end if
            end do
          end do
        end if
      end do
c
      if (locdbg) then
        write(lupri,*)'norm^2 sorted integrals: ',
     &                 ddot(nvir(isyma)*nkl(isymkl),xk,1,xk,1)
C       if (nsym.eq.1) then
C        write(lupri,'(a,3i5)')
C    &     'integrals for fixed isymb,b,isyma:',isymb,b,isyma
C        call output(xk,1,nvir(isyma),1,nkl(1),
C    &               nvir(isyma),nkl(1),1,lupri)
C       end if
      end if
c
      call qexit('r12sortk')
      end
*======================================================================*
      subroutine ccrhs_eppp(omega2,work,lwork,aproxr12,lres,
     &                      lufc12,fc12am,ifile,isymtr)
*----------------------------------------------------------------------*
c     Purpose: add E''' term to vector function in omega2 for the
c              ansaetze 2 and 3
c
c     omega^ab_ij <-- omega^ab_ij 
c                       + sum_kl (K^ab_kl - t^ab_kl) * c^kl_ij
c
c     the contrib. of K is skipped for the approximations A, A'
c
c     for the approximations ??? is in ansatz 2 added:
c
c     omega^ab_ij <-- omega^ab_ij 
c                       + sum_kl (e_k+e_l-e_i-e_j) r^ab_kl * c^kl_ij
c
c     while in ansatz 3 the following term is added:
c 
c     omega^ab_ij <-- omega^ab_ij 
c                       + sum_kl (e_k+e_l-e_a-e_b) r^ab_ij * c^kl_ij
c
c     H. Fliegl, C. Haettig, spring 2004
c     adapted for ansatz 3, Christof Haettig, spring 2005
*----------------------------------------------------------------------*
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <dummy.h>
#include <r12int.h>
#include <ccr12int.h>

      logical locdbg, lres, needkrint
      parameter (locdbg = .false.)
      character*(*) fc12am
      character*(*) aproxr12
      character*10 model
      integer lwork, lwrk1, kend1, ktint, kcamp, krint
      integer lunit, idum, isymkint, iopt, kend2, lwrk2
      integer ifile, lufc12, isymtr

#if defined (SYS_CRAY)
      real omega2(*), work(*), factor, ddot
#else
      double precision omega2(*), work(*), factor, ddot
#endif

      call qenter('rhs_eppp')

      if (locdbg) then
        write(lupri,*) 'entered ccrhs_eppp... isymtr =',isymtr
      end if

      ! is a second integral array of length nt2r12(1) needed?
      needkrint = (aproxr12(1:1).eq."B") 
      if (ianr12.eq.3 .and. 
     &  (aproxr12(1:2).eq."A'".or.aproxr12(1:1).eq."B") ) then
        needkrint=.true.
      end if

      ktint  = 1
      kcamp  = ktint + nt2r12(1)
      kend1  = kcamp + ntr12am(isymtr)
      lwrk1  = lwork - kend1

      kend2  = kend1
      if (needkrint) then
        krint = kend1 
        kend2 = krint + nt2r12(1)
      end if
      lwrk2 = lwork - kend2

      if (lwrk1.lt.0 .or. lwrk2.lt.0) then
        call quit('Insufficient work space in ccrhs_eppp')
      end if

c     ---------------------------------------------
c     read t^ab_kl integrals from file
c     ---------------------------------------------
      lunit = -1
      call gpopen(lunit,fkr12,'unknown',' ','unformatted',idum,.false.)
      read(lunit)(work(ktint+i-1),i=1,nt2r12(1))  
      call gpclose(lunit,'KEEP')

c     ----------------------------------------------------
c     for approximation B substract K^ab_kl from t^ab_kl
c     ----------------------------------------------------
      if (aproxr12(1:1).eq."B") then
        ! read K^ab_kl from file
        lunit = -1
        call gpopen(lunit,ftr12,'unknown',' ','unformatted',
     &              idum,.false.)
        read(lunit)(work(krint+i-1),i=1,nt2r12(1))  
        call gpclose(lunit,'KEEP')
        
        ! substract (t - K) => work(ktint)
        call daxpy(nt2r12(1),-1.0d0,work(krint),1,work(ktint),1)
      end if

c     ------------------------------------------------------------
c     in ansatz 3 add for approximation A' and B the term
c        (eps_k + eps_l - eps_a - eps_b) * r^ab_kl
c     ------------------------------------------------------------
      if (ianr12.eq.3 .and. 
     &   (aproxr12(1:2).eq."A'".or.aproxr12(1:1).eq."B") ) then
        lunit = -1
        call gpopen(lunit,fr12r12,'unknown',' ','unformatted',
     &              idum,.false.)
        read(lunit)(work(krint+i-1),i=1,nt2r12(1))  
        call gpclose(lunit,'KEEP')

        call cc_t2xe(work(krint),work(kend2),lwrk2)

        call daxpy(nt2r12(1),+1.0d0,work(krint),1,work(ktint),1)
      end if

c     ------------------------------------------------------
c     get R12 amplitudes c^mn_kl stored as triangular matrix
c     ------------------------------------------------------
      if (lres) then
        call cc_rvec(lufc12,fc12am,ntr12am(isymtr),ntr12am(isymtr),
     &               ifile,work(kcamp))
        call cclr_diasclr12(work(kcamp),ketscl,isymtr)
      else
c       read r12 amplitudes from file
        iopt = 32
        call cc_rdrsp('R0 ',0,1,iopt,model,dummy,work(kcamp))
      end if

      if (locdbg) then
        write(lupri,*) 'in ccrhs_eppp norm^2 r12 amplitudes',
     &                 ddot(ntr12am(isymtr),work(kcamp),1,work(kcamp),1)
        write(lupri,*)'Omega_aibj before contribution of ccrhs_eppp:',
     &                  ddot(nt2am(isymtr),omega2,1,omega2,1)
        call flshfo(lupri)
      end if

c     ------------------------------------------------------
c     calculate Omega^ab_ij += - sum_kl k^ab_kl * c^kl_ij
c     ------------------------------------------------------
      isymkint = 1
      factor   = 1.0d0
      call cc_r12mi22(omega2,work(ktint),work(kcamp),isymkint,isymtr,
     &                factor,work(kend1),lwrk1) 

      if (locdbg) then
        write(lupri,*)'norm^2 omega_aibj after mi22:', 
     &               ddot(nt2am(isymtr),omega2,1,omega2,1)
        call flshfo(lupri)
      end if

c     ------------------------------------------------------------
c     in ansatz 2 add for approximation A' and B the term
c      omega^ab_ij += sum_kl (e_k+e_l-e_i-e_j) r^ab_kl * c^kl_ij
c     ------------------------------------------------------------
      if (ianr12.eq.2 .and. 
     &   (aproxr12(1:2).eq."A'".or.aproxr12(1:1).eq."B") ) then
       
        if (locdbg) then
          write(lupri,*) 'R12 amplitudes before cc_r12mkediff2:'
          call cc_prpr12(work(kcamp),isymtr,1,.false.)
        end if 
        ! multiply R12 amplitudes with orbital energy differences:
        !   c^kl_ij  <--  (e_k + el - e_i - e_j) * c^kl_ij
        call cc_r12mkediff2(work(kcamp),isymtr,work(kend1),lwrk1)
        if (locdbg) then
          write(lupri,*) 'R12 amplitudes after cc_r12mkediff2:'
          call cc_prpr12(work(kcamp),isymtr,1,.false.)
        end if

        lunit = -1
        call gpopen(lunit,fr12r12,'unknown',' ','unformatted',
     &              idum,.false.)
        read(lunit)(work(ktint+i-1),i=1,nt2r12(1))  
        call gpclose(lunit,'KEEP')

        isymkint = 1
        factor   = 1.0d0
        call cc_r12mi22(omega2,work(ktint),work(kcamp),isymkint,isymtr,
     &                  factor,work(kend1),lwrk1) 
      end if

c     ----------------------------------------------
c     that's it; print some debug output and return:
c     ----------------------------------------------
      if (locdbg) then
        write(lupri,*)'Omega_aibj after ccrhs_eppp contribution:',
     &               ddot(nt2am(isymtr),omega2,1,omega2,1)
      end if

      call qexit('rhs_eppp')
      end
*=====================================================================*
      subroutine cc_r12mi22(omega2,xkint,c2am,isymkint,isymc,
     &                      factor,work,lwork)
c--------------------------------------------------------------------
c     purpose: make Omega_ab_ij
c
c     H. Fliegl, C. Haettig, spring 2004
c--------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
      logical locdbg
      parameter (locdbg = .false.)
      integer lwork, isymkint, isymc, isymb, isyma, isymab, isymkl,
     &        isymij, isymj, isymi, kxint, kc2am, kend1, isymint
      integer ntoti, ntota, kom, isymai, isymbj, maxai, idxai, idxbj,
     &        naibj, koffai, lwrk1
#if defined (SYS_CRAY)
      real work(*), xkint(*),c2am(*), omega2(*),factor
      real one, zero,xk
#else
      double precision work(*), xkint(*),c2am(*), omega2(*),factor
      double precision one, zero,ddot,x1,x2,x3,xk
#endif
      parameter (zero = 0.0d0, one = 1.0d0)

      integer index
      index(i,j) = max(i,j)*(max(i,j) - 3)/2 + i + j

      call qenter('r12mi22')
      isymint = muld2h(isymkint,isymc)

      if (locdbg) then
        write(lupri,*) 'entered cc_r12mi22'
        x1 = 0.0d0
        xk = 0.0d0
        x2 = 0.0d0
        x3 = 0.0d0
      end if 

      do isymb = 1, nsym
        do isyma = 1, nsym
          isymab = muld2h(isyma,isymb)
          isymkl = muld2h(isymkint,isymab) 
          isymij = muld2h(isymkl,isymc)    
          do isymj = 1, nsym
            isymi = muld2h(isymij,isymj)
            isymai = muld2h(isyma,isymi)
            isymbj = muld2h(isymb,isymj)

            kxint = 1
            kc2am = kxint + nvir(isyma)*nmatkl(isymkl)
            kom   = kc2am + nrhf(isymi)*nmatkl(isymkl)
            kend1 = kom   + nvir(isyma)*nrhf(isymi)
            lwrk1 = lwork - kend1
c  
            if (lwrk1.lt.0) then
                call quit('Insufficient memory in cc_r12mi22')
            end if
c
            do b = 1, nvir(isymb)
C             call cc_r12sortk(work(kxint),xkint,b,isymkint,isymb,
C    &                         isyma,nrhfb,.false.)
              call cc_r12sort(xkint,isymkint,work(kxint),b,isymb,
     &                        1,nvir(isyma),isyma,nvir,nrhfb)  
                  xk = xk + ddot(nvir(isyma)*nmatkl(isymkl),
     &                         work(kxint),1,work(kxint),1)
              do j = 1, nrhf(isymj)
                call cc_r12sortc(work(kc2am),c2am,j,isymc,isymj,isymi)
 
                ntoti = max(nrhf(isymi),1)
                ntota = max(nvir(isyma),1)
                call dgemm('N','T',nvir(isyma),nrhf(isymi),
     &                     nmatkl(isymkl),one,work(kxint),ntota,
     &                     work(kc2am),ntoti,zero,work(kom),ntota)
c
                if (locdbg) then
                  x1 = x1 + ddot(nvir(isyma)*nmatkl(isymkl),
     &                         work(kxint),1,work(kxint),1)
                  x2 = x2 + ddot(nrhf(isymi)*nmatkl(isymkl),
     &                           work(kc2am),1,work(kc2am),1)
                  x3 = x3 + ddot(nvir(isyma)*nrhf(isymi),
     &                           work(kom),1,work(kom),1)
                end if 
 
                ! subtract the result from Omega(ai,bj) which is 
                ! stored as a triangular matrix
                if (isymbj.eq.isymai) then
                  idxbj  = it1am(isymb,isymj) + nvir(isymb)*(j-1) + b
                  koffai = it1am(isyma,isymi)
                  maxai = min(koffai+nvir(isyma)*nrhf(isymi),idxbj)

                  do idxai = koffai+1, maxai
                    naibj = it2am(isymai,isymbj) + index(idxai,idxbj)
                    omega2(naibj)= omega2(naibj)+
     &                            factor*work(kom-1+idxai-koffai)
                  end do
c                 
                else if (isymbj.gt.isymai) then
c                 ! omega not totally symmetric ... only needed for
c                 ! response
                  idxbj  = it1am(isymb,isymj) + nvir(isymb)*(j-1) + b
                  koffai = it1am(isyma,isymi)
                  maxai = koffai+nvir(isyma)*nrhf(isymi) 

                  do idxai = koffai+1, maxai
                    naibj = it2am(isymai,isymbj) + 
     &                      nt1am(isymai)*(idxbj-1) + idxai 
                    omega2(naibj)= omega2(naibj)+
     &                            factor*work(kom-1+idxai-koffai)
                  end do
                end if
              end do
            end do
c           
          end do
        end do
      end do

      if (locdbg) then
        write(lupri,*)'xk =', xk
        write(lupri,*)'x1 =', x1
        write(lupri,*)'x2 =', x2
        write(lupri,*)'x3 =', x3
      end if
 
      call qexit('r12mi22')
      end
*=====================================================================*
      subroutine cc_r12sortc(xc,c2am,j,isymc,isymj,isymi)
c--------------------------------------------------------------------
c     purpose: sort c^kl_ij as c^j_i,kl
c
c     H. Fliegl, C. Haettig, spring 2004
c--------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
      logical locdbg
      parameter ( locdbg = .false. ) 
      integer isymc, index, isymj, isymi, isymikl, isymkl, isymk,
     &        isyml, isymki, isymlj, idxlj, idxkl, idxikl, idxki,
     &        idxkilj  
#if defined (SYS_CRAY)
      real c2am(*), xc(*)
      real ddot
#else
      double precision c2am(*), xc(*)
      double precision ddot
#endif
      index(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j
      call qenter('r12sortc')
      
      if (locdbg) then
       write(lupri,*)'entered cc_r12sortc..' 
       write(lupri,*)'norm^2(c2am):',ddot(ntr12am(isymc),c2am,1,c2am,1)
       call cc_prpr12(c2am,isymc,1,.false.)
      end if 

      isymikl = muld2h(isymj,isymc) 
      isymkl  = muld2h(isymi,isymikl)
      do isymk = 1, nsym
        isyml = muld2h(isymkl,isymk)
        isymki = muld2h(isymi,isymk)
        isymlj = muld2h(isymj,isyml)

        do l = 1, nrhfb(isyml)
          idxlj = imatki(isyml,isymj) + nrhfb(isyml)*(j-1)+l
          do k = 1, nrhfb(isymk)
            idxkl = imatkl(isymk,isyml) + nrhfb(isymk)*(l-1)+k
            if (isymki.eq.isymlj) then 
              do i =1, nrhf(isymi)
                idxikl = nrhf(isymi)*(idxkl-1)+i 
                idxki  = imatki(isymk,isymi)+nrhfb(isymk)*(i-1)+k
                idxkilj = itr12am(isymki,isymlj)+index(idxki,idxlj) 
                xc(idxikl) = c2am(idxkilj)
              end do
            else if (isymki.lt.isymlj) then
              idxikl = nrhf(isymi)*(idxkl-1)+1 
              idxki  = imatki(isymk,isymi)+k
              idxkilj = itr12am(isymki,isymlj)+
     &                  nmatki(isymki)*(idxlj-1)+idxki
              call dcopy(nrhf(isymi),c2am(idxkilj),nrhfb(isymk),
     &                   xc(idxikl),1)
            else if (isymki.gt.isymlj) then
              do i =1, nrhf(isymi)
                idxikl = nrhf(isymi)*(idxkl-1)+i
                idxki  = imatki(isymk,isymi)+nrhfb(isymk)*(i-1)+k
                idxkilj = itr12am(isymlj,isymki)+
     &                    nmatki(isymlj)*(idxki-1)+idxlj
                xc(idxikl) = c2am(idxkilj)
              end do
            end if
          end do
        end do
      end do

      if (locdbg) then
         write(lupri,'(a,3i5)')
     &     'R12 amplitudes for fixed isymj,j,isymi:',isymj,j,isymi
         call output(xc,1,nrhf(isymi),1,nmatkl(isymkl),
     &               nrhf(isymi),nmatkl(isymkl),1,lupri)
      end if

      call qexit('r12sortc')
      end
*=====================================================================*
      subroutine cc_r12batf2(work,lwork)
c--------------------------------------------------------------------
c     purpose: prepare R^d_akl for ansatz 2 to obtain (V^kl_aj)^+
c
c     H. Fliegl, C. Haettig, sping 2004
c     
c     C. Neiss, winter 2006: adapted for CABS
c--------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
#include <ccr12int.h>
       
#if defined (SYS_CRAY)
      real zero,one, work(*),ddot
#else
      double precision zero,one, work(*),ddot
#endif
      parameter (zero = 0.0D0, one = 1.0D0)
 
      logical locdbg,lauxd
      parameter (locdbg = .false.)

      integer nr1orb(8),nr1bas(8),nr1xbas(8),nr2bas,n2bst1(8)
      integer ir1xbas(8,8),nxx2(8),lwork,isymmu,kscr1
      integer ir1orb(8,8),ir1bas(8,8),ir2bas(8,8),iaodis1(8,8)
      integer ir2xbas(8,8),irgkl(8,8),nrgkl(8)
      integer nnbastx,n2bst2(8),isymd,isymk,isym,
     &        kend0,lusifc,nsymx,norbsx(8),nbastx,nlamdsx,
     &        nrhfsx(8),kcmox,lwrk0,ksfull,iaodis2(8,8)
      integer idummy,norbtsx,isympq,isymp,isymq,ksjbp,ksmupbp
      integer kend1,lwrk1,isymbp,isymdp,isymakl,krdpakl,krtot,kend2,
     &        lwrk2,ib,id,idxbp,idelp,idxbpdp,koff1,idxd,
     &        nalphaj(8),ialphaj(8,8),nbasx(8),basoff(8)
      integer nr1xorb(8),ir1xorb(8,8),nrxgkl(8),irxgkl(8,8)

      call qenter('batf2')
      if(locdbg) write(lupri,*) 'entered batf2'

      call cc_r12offset(nr1orb,nr1xorb,nr1bas,nr1xbas,nr2bas,
     & nrgkl,nrxgkl,n2bst1,ir1orb,ir1xorb,ir1bas,ir1xbas,ir2bas,
     & ir2xbas,irgkl,irxgkl,iaodis1,nalphaj,ialphaj)
c
      if (r12cbs) then
        do isym = 1, nsym
          nbasx(isym)  = mbas1(isym)+mbas2(isym)
          basoff(isym) = 0
        end do
      else
        call icopy(8,mbas2,1,nbasx,1)
        call icopy(8,mbas1,1,basoff,1)
      end if
c
      do isympq = 1, nsym
        nxx2(isympq) = 0
        do isymq = 1, nsym
          isymp = muld2h(isympq,isymq)
          iaodis2(isymp,isymq) = nxx2(isympq)
          nxx2(isympq) = nxx2(isympq) + nbasx(isymp)*nbasx(isymq)
        end do
      end do
c
      kend0 = 1

c----------------------------------------------------------------------
c     get MO coefficients defined for combined orbital + aux. basis:
c----------------------------------------------------------------------
      lusifc = -1
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
         call quit('Insufficient work space in batf2')
      end if

c     read CMO coefficients and close file
      read(lusifc)
      read(lusifc) (work(kcmox+i-1),i=1,nlamdsx)
      call gpclose(lusifc,'KEEP')

c     test if CMO is ok
      if (locdbg) then
        write(lupri,*)'CMO out batf2'
        kscr1 = 1 
        do isymp = 1, nsym
          isymmu = isymp 
          write(lupri,*) 'symmetry block,nbas,norbs:',isymp,
     &       nbas(isymmu), norbs(isymp)
          call output(work(kcmox+kscr1-1),1,nbas(isymmu),
     &              1,norbs(isymp),nbas(isymmu),norbs(isymp),
     &              1,lupri)

          kscr1 = kscr1 + nbas(isymmu)*norbs(isymp)
        end do
      end if

c     get S_mu'b'
      ksmupbp = kend0
      kend1   = ksmupbp + nxx2(1) 
      lwrk1   = lwork - kend1
      if (lwrk1 .lt.0) then
       call quit('Insufficient work space in batf2')
      end if
      call cc_r12mksmupbp(work(kcmox),work(ksmupbp),iaodis2, 
     &                    work(kend1),lwrk1)

      if (locdbg) then
        write(lupri,*) "norm^2(S_mu'b'):",ddot(nxx2(1),work(ksmupbp),1,
     &       work(ksmupbp),1)
        do isymp = 1, nsym
          write(lupri,*) 'symmetry block:',isymp 
          call output(work(ksmupbp+iaodis2(isymp,isymp)),
     &                1,nbasx(isymp),1,nbasx(isymp),
     &                nbasx(isymp),nbasx(isymp),1,lupri)
        end do
      end if

c-----------------------------------------------------------------------
c     make R^d'_a,kl 
c-----------------------------------------------------------------------
      do isymdp = 1, nsym
        isymbp = isymdp
        isymakl = isymdp

        krdpakl   = kend1
        krtot   = krdpakl + nrgkl(isymakl)
        kend2   = krtot + nrgkl(isymakl)
        lwrk2   = lwork - kend2
        if(lwrk2.lt.0) then
          call quit('Insufficient works space in batf2 ')
        end if
       
        do id = 1, mbas1(isymdp)+mbas2(isymdp)
          lauxd = .false.
          idxd = id + ibas(isymdp)
          if (id .gt. basoff(isymdp)) then
            !make backtransformation...
            call dzero(work(krdpakl),nrgkl(isymakl))
            idelp = id - basoff(isymbp)
            do ib = 1, nbasx(isymbp)
              lauxd = .false.
              idxbp = ib + ibas(isymbp) + basoff(isymbp)
              call cc_r12getrint(work(krtot),idxbp,isymbp,nr1bas,ir1bas,
     &           nr2bas,ir2bas,nrgkl,irgkl,ir1xbas,ir2xbas,
     &           nrhfb,nmatkl,imatkl,
     &           idummy,lauxd,.false.,frhtf,work(kend2),lwrk2)
              idxbpdp = iaodis2(isymbp,isymdp) +
     &                  nbasx(isymbp)*(idelp-1) + ib
              call daxpy(nrgkl(isymakl),work(ksmupbp-1+idxbpdp),
     &                   work(krtot),1,work(krdpakl),1)
            end do
          else
            !only copy...
            call cc_r12getrint(work(krdpakl),idxd,isymdp,nr1bas,ir1bas,
     &           nr2bas,ir2bas,nrgkl,irgkl,ir1xbas,ir2xbas,
     &           nrhfb,nmatkl,imatkl,
     &           idummy,lauxd,.false.,frhtf,work(kend2),lwrk2)
          end if
          call cc_r12putrint(work(krdpakl),idxd,isymdp,nr1bas,ir1bas,
     &         nr2bas,ir2bas,nrgkl,irgkl,ir1xbas,ir2xbas,
     &         nrhfb,nmatkl,imatkl,
     &         idummy,lauxd,fnback2,work(kend2),lwrk2)
        end do
      end do
c
      call qexit('batf2')
      end
*=====================================================================*
      subroutine cc_r12mkvaj2(vajkl,isymv,cmo,isymc,cmos,isymcs,
     &                        work,lwork)
c---------------------------------------------------------------------
c     purpose: make \sum_mn r_mn^kl*g_aj^mn
c
c     H. Fliegl, C. Haettig, spring 2004
c---------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
#include <ccr12int.h>
#include <dummy.h>
      logical locdbg, lauxd
      parameter (locdbg = .false.) 
      integer isymv, isymc, lwork, lwrk1, lwrk2, kend1, kend2 
      integer krmnkl, kgmakl, krhtf, koffv, koffg, koffr, idxkl, idxd,
     &        ntotaj, ntotmn, idelta, isymd, isymnm, isymaj, isymkl
      integer nr1orb(8), nr1bas(8), nr1xbas(8), nr2bas, n2bst1(8)
      integer ir1xbas(8,8), ibasx(8),kmkl,iadr1,iadr2,isymdm,isymm,ndim
      integer iglmrhs(8,8),iglmvis(8,8),nglmds(8),icmo(8,8),ncmo(8),
     &        imaklm(8,8),nmaklm(8),imatijm(8,8),nmatijm(8),
     &        igamsm(8,8),ngamsm(8),ivajkls(8,8),nvajkls(8)
      integer ir1orb(8,8), ir1bas(8,8), ir2bas(8,8), iaodis1(8,8)
      integer ir2xbas(8,8), irgkl(8,8), nrgkl(8), isyma,isymj
      integer irgijs(8,8),nrgijs(8),ir1basm(8,8),nr1basm(8),
     &        ir2basm(8,8),nr2basm,ir1xbasm(8,8),nr1xbasm(8),
     &        ir2xbasm(8,8),isym,icount1,isym2,isym1,isymg,koffc,
     &        koff1,ntotg,ntotm,imatf(8,8),nmatf(8)
      integer ioffg(8,8),isymn,isymmn,idxmn,kgtf,idxdn,idxlk,isymk,
     &        isyml,isymcs
      integer nr1xorb(8),ir1xorb(8,8),nrxgkl(8),irxgkl(8,8),
     &        nalphaj(8),ialphaj(8,8),nmaijm(8),imaijm(8,8)
      double precision one, zero
      parameter (one = 1.0d0, zero = 0.0d0)
#if defined (SYS_CRAY)
      real work(*), vajkl(*), cmo(*),ddot,ddummy,cmos(*),factor
#else
      double precision work(*),vajkl(*),cmo(*),ddot,ddummy,cmos(*),
     &                 factor
#endif 

      call qenter('r12mkvaj2')
      call cc_r12offset(nr1orb,nr1xorb,nr1bas,nr1xbas,nr2bas,
     &     nrgkl,nrxgkl,n2bst1,ir1orb,ir1xorb,ir1bas,ir1xbas,
     &     ir2bas,ir2xbas,irgkl,irxgkl,iaodis1,nalphaj,ialphaj) 
      call cc_r12offs23(iglmrhs,iglmvis,nglmds,icmo,ncmo,
     &             imaijm,nmaijm,imaklm,nmaklm,
     &             imatijm,nmatijm,igamsm,ngamsm,irgijs,nrgijs,
     &             ir1basm,nr1basm,ir2basm,nr2basm,ir1xbasm,
     &             nr1xbasm,ir2xbasm,imatf,nmatf)
c
      ibasx(1) = 0
      do isym = 2, nsym
        ibasx(isym) = ibasx(isym-1) + mbas2(isym-1)
      end do
c
      do isym = 1, nsym
        icount1 = 0
        do isym2 = 1, nsym
          isym1 = muld2h(isym,isym2)
          ivajkls(isym1,isym2) = icount1
          icount1 = icount1 + nalphaj(isym1)*nmatijm(isym2)
        end do
        nvajkls(isym) = icount1
      end do
c     
      if (locdbg) then
        do isym = 1, nsym
          icount1 = 0
          do isym2 = 1, nsym
            isym1 = muld2h(isym,isym2)
            ioffg(isym1,isym2) = icount1
            icount1 = icount1 + nmatij(isym1)*nmatijm(isym2)
          end do
        end do
      end if
      
      krmnkl = 1
      kgmakl = krmnkl + ngamsm(1)
      kend1  = kgmakl + nvajkls(1) 
      lwrk1  = lwork - kend1
      if (lwrk1.lt.0) then
        call quit('Insufficient work space in mkvaj2')
      end if
 
      call dzero(work(krmnkl),ngamsm(1))
      call dzero(work(kgmakl),nvajkls(1))
c
      do isymd = 1, nsym
        krhtf = kend1
        kmkl  = krhtf + max(nrgkl(isymd),nrgijs(isymd))
        kend2 = kmkl + nmatf(isymd)
        lwrk2 = lwork - kend2
        if (lwrk2.lt.0) then
          call quit('Insufficient work space in mkvaj2')
        end if
        do idelta = 1, mbas1(isymd)
          lauxd = .false.
          idxd = idelta+ibas(isymd)
          call cc_r12getrint(work(krhtf),idxd,isymd,nr1bas,ir1bas,
     &         nr2bas,ir2bas,nrgkl,irgkl,ir1xbas,ir2xbas,
     &         nrhfb,nmatkl,imatkl,
     &         ibasx,lauxd,.false.,frhtf,work(kend2),lwrk2)

c         get rMNkl
          call cc_r12generaltf(work(krhtf),work(krmnkl),idelta,isymd,
     &                         cmos,isymcs,cmos,isymcs,.false.,iglmrhs,
     &                         iglmrhs,nmatkl,nrhfsa,nbas,irgkl,
     &                         imatijm,nmatijm,imaklm,igamsm,
     &                         work(kend2),lwrk2)
c
          lauxd = .false.
          idxd = idelta+ibas(isymd)
c         get g^delta_gamma,KL
          call cc_r12getrints(work(krhtf),idxd,isymd,nr1basm,ir1basm,
     &           nr2basm,ir2basm,nrgijs,irgijs,ir1xbasm,ir2xbasm,
     &           ibasx,imatijm,nmatijm,lauxd,fccgmnab,work(kend2),lwrk2)
c
          if (locdbg) then
            write(lupri,*) 'norm^2(Gmbnd) after read:',idxd,
     &              ddot(nrgijs(isymd),work(krhtf),1,work(krhtf),1)
          end if
c
c        transform g^delta_gamma,KL to g^delta_m,KL
         do isymg = 1, nsym
           isymkl = muld2h(isymg,isymd)
           isymm  = muld2h(isymc,isymg)

c          koffc = iglmrhs(isymg,isymm)+1+nrhffr(isymm)*nbas(isymg)
           koffc = iglmrh(isymg,isymm)+1
           koffg = krhtf + irgijs(isymg,isymkl)
           koff1 = kmkl  + imatf(isymkl,isymm) 

           ntotg = max(nbas(isymg),1)
           ntotm = max(nrhf(isymm),1)

           call dgemm('T','N',nrhf(isymm),nmatijm(isymkl),nbas(isymg),
     &                one,cmo(koffc),ntotg,work(koffg),ntotg,zero,
     &                work(koff1),ntotm)
         end do
c
         do isymkl = 1, nsym
          isymm = muld2h(isymkl,isymd)
          isymdm = muld2h(isymd,isymm)
          do isymk = 1, nsym
            isyml = muld2h(isymkl,isymk)
            do k = 1, nrhfsa(isymk)
              do l = 1, nrhfsa(isyml)
                idxkl = imatijm(isymk,isyml)+nrhfsa(isymk)*(l-1)+k 
                idxlk = imatijm(isyml,isymk)+nrhfsa(isyml)*(k-1)+l
                iadr1 = kgmakl-1+ivajkls(isymdm,isymkl)+
     &                    nalphaj(isymdm)*(idxlk-1)+
     &                    ialphaj(isymd,isymm)+idelta
                iadr2 = kmkl+imatf(isymkl,isymm)+nrhf(isymm)*(idxkl-1)
                call dcopy(nrhf(isymm),work(iadr2),1,
     &                     work(iadr1),mbas1(isymd))
              end do
            end do
          end do
         end do 
c
        end do
      end do
    
      if (locdbg) then
        write(lupri,*)'work(krmnkl)'
c       write(lupri,*)(work(krmnkl-1+i),i=1,ngamsm(1))
        write(lupri,*)'norm^2 work(krmnkl):', 
     &      ddot(ngamsm(1),work(krmnkl),1,work(krmnkl),1)
        write(lupri,*)'work(kgmakl)'
c       write(lupri,*)(work(kgmakl-1+i),i=1,nvajkls(1))
        write(lupri,*)'norm^2 work(kgmakl):', 
     &      ddot(nvajkls(1),work(kgmakl),1,work(kgmakl),1)
      end if
c
c     make Vajkl
c
      factor = one
      do isymnm = 1, nsym
        isymkl = isymnm
        isymaj = muld2h(isymv,isymkl)
        
        koffv = 1 + ivajkl(isymaj,isymkl)  
        koffr = krmnkl + igamsm(isymnm,isymkl)
        koffg = kgmakl + ivajkls(isymaj,isymnm)
        
        ntotmn = max(1,nmatijm(isymnm))
        ntotaj = max(1,nalphaj(isymaj))
        
        call dgemm('N','N',nalphaj(isymaj),nmatkl(isymkl),
     &             nmatijm(isymnm),factor,work(koffg),ntotaj,
     &             work(koffr),ntotmn,one,vajkl(koffv),ntotaj)
        if (locdbg) then
          write(lupri,*)'DEBUG: Vajkl in mkvaj2'
          call output(vajkl(koffv),1,nalphaj(isymaj),1,
     &                nmatkl(isymkl),nalphaj(isymaj),
     &                nmatkl(isymkl),1,lupri)
        end if
      end do
c
      call qexit('r12mkvaj2')
      end
*=====================================================================*
      subroutine ccrhs_hp(omega1,lamdh,isymh,lamdh0,isymh0,work,
     &                    lwork,iamp,isymc,fc12am,lufc12,ifile,iopte)
c---------------------------------------------------------------------
c     purpose: update Omega_ai - = 
c              \sum_alpha Lambda^h_alphaa * Omega_alphai 
c
c     H. Fliegl, C. Haettig, spring 2004
c     adapted for CABS, C. Neiss, winter 2006
c
c     iopte = 1 make E_ab^(R12) intermediate instead of Omega_ai^H'
c               C. Neiss, spring 2006
c---------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
#include <ccr12int.h>
#include <dummy.h>
      integer lwork,isymc,isymh,ifile,lufc12,komai, isymom,
     &        isymh0,kctil,kend1,lwrk1,isymd,isymilk,isymimn,
     &        idelp,kint,kgtf,kmint,kend2,lwrk2,isymkl,
     &        isymi,isymmn,isyma,isymb,koffr,koffm,koffom,
     &        koffom1,ntota,ntotb,koffc,koffg,ntotkl,ntotmn,
     &        koffl,ibasx(8),isym,nbas2(8),basoff(8) 
      integer nr1orb(8),nr1bas(8),nr1xbas(8),nr2bas,n2bst1(8)
      integer ir1xbas(8,8)
      integer ir1orb(8,8),ir1bas(8,8),ir2bas(8,8),iaodis1(8,8)
      integer ir2xbas(8,8),irgkl(8,8),nrgkl(8),nrxgkl(8),irxgkl(8,8)
      integer nr1xorb(8),ir1xorb(8,8) 
      integer idxdp,ntoti,nalphaj(8),ialphaj(8,8),nr1(8)
      integer icount,isym1,isym2,nmaijx(8),imaijx(8,8)
      integer nmaklx(8),imaklx(8,8),imatax(8,8),nmatax(8)
      integer nrgij(8),irgij(8,8),it1xao(8,8),nt1xao(8),it2xaos(8,8)
      integer iamp,iopte
#if defined (SYS_CRAY)
      real work(*), omega1(*), lamdh(*), lamdh0(*), one
      real zero,ddot
#else
      double precision work(*), omega1(*), lamdh(*),lamdh0(*), one
      double precision zero,ddot
#endif
      parameter(one = 1.0d0, zero = 0.0d0)
      character*(*) fc12am
      character cdummy*(8)
      logical locdbg, lauxd
      parameter(locdbg = .false.)

      call qenter('ccrhs_hp')

      isymom = muld2h(isymh,isymc)

      if (locdbg) then
        if (iopte.eq.1) then
          write(lupri,*)'in ccrhs_hp: E1IM on entry '
          call cc_prei(omega1,dummy,isymom,1)
        else
          write(lupri,*)'in ccrhs_hp: omega1 on entry '
          call cc_prp(omega1,dummy,isymom,1,0)
        end if
      end if  

      call cc_r12offset(nr1orb,nr1xorb,nr1bas,nr1xbas,nr2bas,
     & nrgkl,nrxgkl,n2bst1,ir1orb,ir1xorb,ir1bas,ir1xbas,ir2bas,
     & ir2xbas,irgkl,irxgkl,iaodis1,nalphaj,ialphaj) 

      if (r12cbs) then
        do isym = 1, nsym
          nbas2(isym)  = mbas1(isym) + mbas2(isym)
          basoff(isym) = 0
        end do
      else
        call icopy(8,mbas2,1,nbas2,1)
        call icopy(8,mbas1,1,basoff,1)
      end if

      if (iopte.eq.1) then
        call icopy(8,nvir,1,nr1,1)
      else
        call icopy(8,nrhf,1,nr1,1)
      end if

      do isym = 1, nsym
        nmaijx(isym) = 0
        nmaklx(isym) = 0
        nmatax(isym) = 0
        do isym2 = 1,nsym
          isym1 = muld2h(isym,isym2)
          imaijx(isym1,isym2) = nmaijx(isym)
          imaklx(isym1,isym2) = nmaklx(isym)
          imatax(isym1,isym2) = nmatax(isym)
          nmaijx(isym) = nmaijx(isym) + nmatij(isym1)*nr1(isym2)
          nmaklx(isym) = nmaklx(isym) + nmatkl(isym1)*nr1(isym2) 
          nmatax(isym) = nmatax(isym) + mbas1(isym1)*nr1(isym2)
        end do
      end do
      do isym = 1, nsym
        nrgij(isym) = 0
        nt1xao(isym)= 0
        do isym2 = 1, nsym
          isym1 = muld2h(isym,isym2)
          irgij(isym1,isym2)  = nrgij(isym)
          it1xao(isym1,isym2) = nt1xao(isym)
          nrgij(isym) = nrgij(isym) + mbas1(isym1)*nmatij(isym2)
          nt1xao(isym)= nt1xao(isym) + mbas2(isym1)*nrhf(isym2)
        end do
      end do
      do isym = 1, nsym
        icount = 0
        do isym2 = 1, nsym
          isym1 = muld2h(isym,isym2)
          it2xaos(isym1,isym2) = icount
          icount = icount + nt1ao(isym1)*nt1xao(isym2)
        end do
      end do

      komai = 1
      kctil = komai + nmatax(isymom)
      kend1 = kctil + ntr12sq(isymc)
      lwrk1 = lwork - kend1
      if (lwrk1 .lt.0) then
         call quit('Insufficient work space in ccrhs_hp')
      end if

      ibasx(1) = 0
      do isym = 2,nsym
        ibasx(isym) = ibasx(isym-1)+mbas2(isym-1)
      end do         
c
      call dzero(work(komai),nmatax(isymom))
 
c     read r12 amplitudes or trial vector
      call cc_r12getct(work(kctil),isymc,iamp,ketscl,.true.,'T',
     &                 lufc12,fc12am,ifile,cdummy,idummy,work(kend1),
     &                 lwrk1)    

      do isymd = 1, nsym
        isymilk = muld2h(isymh,isymd)
        isymimn = muld2h(isymc,isymilk)
        kint  = kend1
        kgtf  = kint  + max(nrgkl(isymd),nrgij(isymd))
        kmint = kgtf  + nmaijx(isymilk)
        kend2 = kmint + nmaklx(isymimn)
        lwrk2 = lwork - kend2
        if (lwrk2.lt.0) then
           call quit('Insufficient work space in ccrhs_hp')
        end if
c
        do idelp = 1, nbas2(isymd)   
          call dzero(work(kgtf),nmaijx(isymilk))
c         get I^dp_akl
          lauxd = .false.
          idxdp = idelp+ibas(isymd)+basoff(isymd)
          call cc_r12getrint(work(kint),idxdp,isymd,nt1ao,it1ao,
     &         nt2aos,it2aos,nrgij,irgij,it1xao,it2xaos,
     &         nrhf,nmatij,imatij,
     &         ibasx,lauxd,.false.,fgdpakl,work(kend2),lwrk2)
c         if IOPTE.EQ.1 transform I^dp_akl to I^dp_bkl
          if (iopte.eq.1) then
           call cc_r12generaltf(work(kint),work(kgtf),idelp,isymd,
     &                          lamdh,isymh,lamdh,isymh,.true.,iglmvi,
     &                          iglmvi,nmatij,nr1,nbas,irgij,imatij,
     &                          nmatij,imaijx,idummy,work(kend2),lwrk2)
          else
c         transform I^dp_akl to I^dp_ikl
           call cc_r12generaltf(work(kint),work(kgtf),idelp,isymd,
     &                          lamdh,isymh,lamdh,isymh,.true.,iglmrh,
     &                          iglmrh,nmatij,nrhf,nbas,irgij,imatij,
     &                          nmatij,imaijx,idummy,work(kend2),lwrk2)
          end if
c
          if (locdbg) then
            write(lupri,*)'norm^2 coulomb integrals:', 
     &          ddot(nrgij(isymd),work(kint),1,work(kint),1)
            write(lupri,*)'norm^2 kgtf:', 
     &          ddot(nmaijx(isymilk),work(kgtf),1,work(kgtf),1)
          end if

c         transform M^dpi_mn = /sum_kl C^tilde * I^dp_ikl 
          do isymkl = 1, nsym
            isymi  = muld2h(isymilk,isymkl) 
            isymmn = muld2h(isymimn,isymi)
          
            koffc  = kctil + itr12sqt(isymkl,isymmn)
            koffg  = kgtf  + imaijx(isymkl,isymi)
            koffm  = kmint + imaklx(isymmn,isymi)
          
            ntotkl = max(nmatij(isymkl),1)
            ntotmn = max(nmatkl(isymmn),1)
            ntoti  = max(nr1(isymi),1)
          
            call dgemm('T','T',nmatkl(isymmn),nr1(isymi),
     &                nmatij(isymkl),one,work(koffc),ntotkl,
     &                work(koffg),ntoti,zero,work(koffm),ntotmn)
          end do
          if (locdbg) then
            write(lupri,*)'norm^2 kmint:',
     &       ddot(nmaklx(isymimn),work(kmint),1,work(kmint),1)
          end if

c         get R^dp_amn
          lauxd = .false.
          idxdp = idelp+ibas(isymd)+basoff(isymd)
          call cc_r12getrint(work(kint),idxdp,isymd,nr1bas,ir1bas,
     &       nr2bas,ir2bas,nrgkl,irgkl,ir1xbas,ir2xbas,
     &       nrhfb,nmatkl,imatkl,
     &       ibasx,lauxd,.false.,fnback2,work(kend2),lwrk2)

          if (locdbg) then
            write(lupri,*)'norm^2 kint:', 
     &           ddot(nrgkl(isymd),work(kint),1,work(kint),1)
          end if
          
          do isymmn = 1, nsym
            isyma  = muld2h(isymmn,isymd)
            isymi  = muld2h(isymom,isyma) 
          
            koffr  = kint + irgkl(isyma,isymmn)
            koffm  = kmint + imaklx(isymmn,isymi)
            koffom = komai + imatax(isyma,isymi)
          
            ntota  = max(mbas1(isyma),1)
            ntotmn = max(nmatkl(isymmn),1)
          
            call dgemm('N','N',mbas1(isyma),nr1(isymi),
     &                 nmatkl(isymmn),-one,work(koffr),ntota,
     &                 work(koffm),ntotmn,one,work(koffom),ntota)
          end do
        end do
      end do 
      if (locdbg) then
        write(lupri,*)'norm^2 komai:', 
     &     ddot(nmatax(isymom),work(komai),1,work(komai),1)
      end if
c
c     update Omega_ai
      do isyma = 1, nsym
        isymb = muld2h(isymh0,isyma)
        isymi = muld2h(isymom,isyma)  
        koffl = 1 + iglmvi(isyma,isymb)
        koffom = komai + imatax(isyma,isymi) 
        if (iopte.eq.1) then
          koffom1 = 1 + imatab(isymb,isymi)
        else
          koffom1 = 1 + it1am(isymb,isymi)
        end if
          
        ntota = max(mbas1(isyma),1)
        ntotb = max(nvir(isymb),1)
          
        call dgemm('T','N',nvir(isymb),nr1(isymi),mbas1(isyma),one,
     &             lamdh0(koffl),ntota,work(koffom),ntota,one,
     &             omega1(koffom1),ntotb)
      end do
c
      if (locdbg) then
        if (iopte.eq.1) then
          write(lupri,*)'in ccrhs_hp: E1IM on exit '
          call cc_prei(omega1,dummy,isymom,1)
        else
          write(lupri,*)'in ccrhs_hp: omega1 on exit '
          call cc_prp(omega1,dummy,isymom,1,0)
        end if
      end if  

      call qexit('ccrhs_hp')
      end 
*=====================================================================*
      subroutine ccrhs_ip(omega1,t1am,isymt,lamdh,isymh,iamp,isymc,
     &                    fc12am,lufc12,ifile,work,lwork)
c---------------------------------------------------------------------
c     purpose: update Omega_ai = /sum_klm ctilde * M^ak_mn
c
c     H. Fliegl, C. Haettig, spring 2004
c     adapted for CABS, C. Neiss, winter 2006
c---------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
#include <ccr12int.h>  
#include <dummy.h> 
      logical lauxd, locdbg
      parameter (locdbg = .false.)
      integer lwork,ifile,lufc12,isymt,isymh,isymc,isymt1ao,
     &        isymak,icount7,isymk,isyma,kt1ao,kmint,kend1,
     &        lwrk1,isymb,isyml,koffl,kofft,kofftam,ntota,
     &        ntotb,isymd,idelp,kint,kfint,kend2,lwrk2,
     &        ibasx(8),isymkl,idxkl,idxlk,kofft1am,koffgc,
     &        koffge,isymmn,idxmn,koffr,koffm,kctil,isymomg,isym
      integer nr1orb(8),nr1bas(8),nr1xbas(8),nr2bas,n2bst1(8)
      integer ir1xbas(8,8),nbas2(8),basoff(8)
      integer ir1orb(8,8),ir1bas(8,8),ir2bas(8,8),iaodis1(8,8)
      integer ir2xbas(8,8),irgkl(8,8),nrgkl(8),nrxgkl(8),irxgkl(8,8)
      integer nr1xorb(8),ir1xorb(8,8),nalphaj(8),ialphaj(8,8) 
      integer idxdp,nrgij(8),irgij(8,8),
     &        nt1xao(8),it1xao(8,8),it2xaos(8,8),isymdk,icount,icoun1
      integer iamp
#if defined (SYS_CRAY) 
      real omega1(*), work(*), lamdh(*), t1am(*), 
     &     ddot, er12, one, zero
#else 
      double precision omega1(*), work(*), lamdh(*), t1am(*),
     &                ddot, er12, one, zero
#endif 
      parameter (one = 1.0d0, zero = 0.0d0)
      character*(*) fc12am
      character cdummy*(8)

      call qenter('ccrhs_ip')

      call cc_r12offset(nr1orb,nr1xorb,nr1bas,nr1xbas,nr2bas,
     & nrgkl,nrxgkl,n2bst1,ir1orb,ir1xorb,ir1bas,ir1xbas,
     & ir2bas,ir2xbas,irgkl,irxgkl,iaodis1,nalphaj,ialphaj) 

      isymt1ao = muld2h(isymt,isymh)
       
      if (locdbg) then
        write(lupri,*)'in ccrhs_ip: omega1 on entry '
        call cc_prp(omega1,dummy,isymt1ao,1,0)
      end if 
c
      if (r12cbs) then
        do isym = 1, nsym
          nbas2(isym)  = mbas1(isym) + mbas2(isym)
          basoff(isym) = 0
        end do
      else
        call icopy(8,mbas2,1,nbas2,1)
        call icopy(8,mbas1,1,basoff,1)
      end if
c
      ibasx(1) = 0
      do isym = 2,nsym
        ibasx(isym) = ibasx(isym-1)+mbas2(isym-1)
      end do
      do isymdk = 1, nsym
        icount = 0
        icoun1 = 0
        do isymk = 1, nsym
          isymd = muld2h(isymdk,isymk)
          irgij(isymd,isymk)  = icount
          it1xao(isymd,isymk) = icoun1
          icount = icount + mbas1(isymd)*nmatij(isymk)
          icoun1 = icoun1 + mbas2(isymd)*nrhf(isymk)
        end do
        nrgij(isymdk) = icount
        nt1xao(isymdk) = icoun1
      end do
      do isymdk = 1, nsym
        icount = 0
        do isymk = 1, nsym
          isymd = muld2h(isymdk,isymk)
          it2xaos(isymd,isymk) = icount
          icount = icount + nt1ao(isymd)*nt1xao(isymk)
        end do
      end do
c
      kt1ao = 1
      kmint = kt1ao + nt1ao(isymt1ao)
      kend1 = kmint + nvajkl(isymt1ao)
      lwrk1 = lwork - kend1
      if (lwrk1 .lt.0) then
         call quit('Insufficient work space in ccrhs_ip')
      end if 
c     
      call dzero(work(kmint),nvajkl(isymt1ao))
      call dzero(work(kt1ao),nt1ao(isymt1ao))
c
c     calculate t_al = /sum_b t_bl * Lamda^h_ab
      do isymb = 1, nsym
        isyma = muld2h(isymh,isymb)
        isyml = muld2h(isymt,isymb)

        koffl = 1 + iglmvi(isyma,isymb)
        kofft = 1 + it1am(isymb,isyml)
        kofftam = kt1ao + it1ao(isyma,isyml)
        
        ntota = max(mbas1(isyma),1)
        ntotb = max(nvir(isymb),1)

        call dgemm('N','N',mbas1(isyma),nrhf(isyml),nvir(isymb),one,
     &             lamdh(koffl),ntota,t1am(kofft),ntotb,zero,
     &             work(kofftam),ntota)
      end do      
 
      do isymd = 1, nsym
        isymk = muld2h(isymt1ao,isymd)
        do idelp = 1, nbas2(isymd)
          kint  = kend1
          kfint = kint + nrgkl(isymd)
          kend2 = kfint + nrhf(isymk)
          lwrk2 = lwork - kend2
          if (lwrk2.lt.0) then
             call quit('Insufficient work space in ccrhs_hp')
          end if    

c         get I^dp_akl
          lauxd = .false.
          idxdp = idelp+ibas(isymd)+basoff(isymd)
          call cc_r12getrint(work(kint),idxdp,isymd,nt1ao,it1ao,
     &       nt2aos,it2aos,nrgij,irgij,it1xao,it2xaos,
     &       nrhf,nmatij,imatij,
     &       ibasx,lauxd,.false.,fgdpakl,work(kend2),lwrk2)
          if (locdbg) then
            write(lupri,*)'norm^2 kint:',
     &        ddot(nrgij(isymd),work(kint),1,work(kint),1)
          end if

c         get F intermediate      
          call dzero(work(kfint),nrhf(isymk))
          do isyml = 1, nsym
            isyma = muld2h(isymt1ao,isyml)
            isymkl = muld2h(isymk,isyml)
c
            do k = 1, nrhf(isymk)
              do l = 1, nrhf(isyml)
                idxkl = imatij(isymk,isyml)+nrhf(isymk)*(l-1)+k
                idxlk = imatij(isyml,isymk)+nrhf(isyml)*(k-1)+l

                kofft1am = kt1ao + it1ao(isyma,isyml)+
     &                     mbas1(isyma)*(l-1) 
                koffgc   = kint + irgij(isyma,isymkl)+
     &                     mbas1(isyma)*(idxlk-1)  
                koffge   = kint + irgij(isyma,isymkl)+ 
     &                     mbas1(isyma)*(idxkl-1)  

                work(kfint-1+k) = work(kfint-1+k) + 
     &                           2.0d0*ddot(mbas1(isyma),
     &                           work(kofft1am),1,work(koffgc),1) -
     &                           ddot(mbas1(isyma),work(kofft1am),1,
     &                           work(koffge),1)
              end do
            end do
          end do
          if (locdbg) then
            write(lupri,*)'norm^2: kfint', 
     &        ddot(nrhf(isymk),work(kfint),1,work(kfint),1)
          end if

c         get R^dp_amn
          lauxd = .false.
          idxdp = idelp+ibas(isymd)+basoff(isymd)
          call cc_r12getrint(work(kint),idxdp,isymd,nr1bas,ir1bas,
     &       nr2bas,ir2bas,nrgkl,irgkl,ir1xbas,ir2xbas,
     &       nrhfb,nmatkl,imatkl,
     &       ibasx,lauxd,.false.,fnback2,work(kend2),lwrk2)
  
c         get M intermediate
          do isymmn = 1, nsym
            isyma = muld2h(isymd,isymmn)
            isymak = muld2h(isyma,isymk)
            do k = 1, nrhf(isymk)
              do idxmn = 1, nmatkl(isymmn)
                koffr = kint+irgkl(isyma,isymmn)+
     &                  mbas1(isyma)*(idxmn-1)
                koffm = kmint+ivajkl(isymak,isymmn)+
     &                  nt1ao(isymak)*(idxmn-1)+it1ao(isyma,isymk)+
     &                  mbas1(isyma)*(k-1)

                call daxpy(mbas1(isyma),work(kfint-1+k),work(koffr),
     &                     1,work(koffm),1)
              end do
            end do  
          end do
c        
        end do
      end do
      if (locdbg) then
        write(lupri,*) 'norm^2(M):',ddot(nvajkl(isymt1ao),
     &      work(kmint),1,work(kmint),1)
      end if

      kctil = kend1
      kend2 = kctil + ntr12sq(isymc)
      lwrk2 = lwork - kend2
      if (lwrk2.lt.0) then
        call quit('Insufficient work space in ccrhs_ip')
      end if

c     read r12 amplitudes or trial vector
      call cc_r12getct(work(kctil),isymc,iamp,ketscl,.true.,'T',
     &                 lufc12,fc12am,ifile,cdummy,idummy,work(kend2),
     &                 lwrk2)   
  
c     isymomg = muld2h(isymc,muld2h(isymt1ao,isymh))
      isymomg = muld2h(isymh,muld2h(isymt1ao,isymc))
      call ccrhs_gp0(omega1,isymomg,work(kmint),isymt1ao,lamdh,isymh,
     &         work(kctil),isymc,.false.,er12,locdbg,work(kend2),lwrk2)

      if (locdbg) then
        write(lupri,*)'in ccrhs_ip: omega1 on exit '
        call cc_prp(omega1,dummy,isymt1ao,1,0)
      end if  

      call qexit('ccrhs_ip')
      end
*======================================================================*
      subroutine cc_r12offs23(iglmrhs,iglmvis,nglmds,icmo,ncmo,
     &                        imaijm,nmaijm,
     &                        imaklm,nmaklm,imatijm,nmatijm,
     &                        igamsm,ngamsm,irgijs,nrgijs,
     &                        ir1basm,nr1basm,ir2basm,nr2basm,ir1xbasm,
     &                        nr1xbasm,ir2xbasm,imatf,nmatf)
*----------------------------------------------------------------------*
c     purpose: calculate some offsets needet for ansatz 2 & 3 to include
c              all active and inactive orbitals
c
c     H. Fliegl, C. Haettig, spring 2004
*----------------------------------------------------------------------*
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
      integer iglmrhs(8,8),iglmvis(8,8),nglmds(8),isym,
     &        icmo(8,8),ncmo(8),imaklm(8,8),nmaklm(8),
     &        imatijm(8,8),nmatijm(8),ngamsm(8),igamsm(8,8),
     &        isym2,isym1,irgijs(8,8),nrgijs(8),ir1basm(8,8),nr1basm(8),
     &        ir2basm(8,8),nr2basm,ir1xbasm(8,8),nr1xbasm(8),
     &        ir2xbasm(8,8),imatf(8,8),nmatf(8),imaijm(8,8),nmaijm(8),
     &        icount1,icount2,icount2a,icount3,icount4,icount5,
     &        icount6,icount7

      call qenter('ccoffs23')

      do isym = 1, nsym
     
        ! precalculate dimensions of transformation matrices
        ! and the offset for the virtual block
        icount1 = 0
        icount2a= 0
        do isym2 = 1, nsym
          isym1 = muld2h(isym,isym2)

          icmo(isym1,isym2) = icount1

          icount1 = icount1 + nbas(isym1)*norbs(isym2)
          icount2a= icount2a+ nbas(isym1)*nrhfsa(isym2)
        end do
        ncmo(isym)    = icount1
        nglmds(isym)  = icount1

        icount2 = 0
        icount3 = 0
        icount4 = 0
        icount5 = 0
        icount6 = 0
        icount7 = 0
        do isym2 = 1, nsym
          isym1 = muld2h(isym,isym2)
          iglmrhs(isym1,isym2) = icount2
          iglmvis(isym1,isym2) = icount2a
          imaklm(isym1,isym2)  = icount3
          imaijm(isym1,isym2)  = icount7
          imatijm(isym1,isym2) = icount4
          ir1basm(isym1,isym2) = icount5
          ir1xbasm(isym1,isym2) = icount6

          icount2 = icount2 + nbas(isym1)*nrhfsa(isym2)
          icount2a= icount2a+ nbas(isym1)*nvirs(isym2)
          icount3 = icount3 + nmatkl(isym1)*nrhfsa(isym2)
          icount7 = icount7 + nmatij(isym1)*nrhfsa(isym2)
          icount4 = icount4 + nrhfs(isym1)*nrhfsa(isym2)
          icount5 = icount5 + mbas1(isym1)*nrhfsa(isym2)
          icount6 = icount6 + mbas2(isym1)*nrhfsa(isym2)
        end do
        nmaklm(isym)  = icount3
        nmaijm(isym)  = icount7
        nmatijm(isym) = icount4
        nr1basm(isym) = icount5 
        nr1xbasm(isym) = icount6
      end do

      do isym = 1, nsym
        icount7 = 0
        do isym2 = 1, nsym
          isym1 = muld2h(isym,isym2)
          imatf(isym1,isym2)  = icount7
          icount7 = icount7 + nmatijm(isym1)*nrhf(isym2)
        end do
        nmatf(isym)   = icount7
      end do

      do isym = 1, nsym
        icount1 = 0
        icount2 = 0
        icount3 = 0
        icount4 = 0
        do isym2 = 1, nsym
          isym1 = muld2h(isym,isym2)
          igamsm(isym1,isym2)  = icount1
          irgijs(isym1,isym2)  = icount2
          ir2basm(isym1,isym2) = icount3
          ir2xbasm(isym1,isym2) = icount4

          icount1 = icount1 + nmatijm(isym1)*nmatkl(isym2)
          icount2 = icount2 + mbas1(isym1)*nmatijm(isym2)
          icount3 = icount3 + nr1basm(isym1)*nr1basm(isym2)
          icount4 = icount4 + nr1basm(isym1)*nr1xbasm(isym2)
        end do
        ngamsm(isym)  = icount1
        nrgijs(isym)  = icount2
c       nr2basm(isym) = icount3
      end do
c
      nr2basm = 0
      do isym = 1, nsym  
        nr2basm = nr2basm + nr1basm(isym)*nr1basm(isym)
      end do    

      call qexit('ccoffs23')

      end
*=====================================================================*
      subroutine cc_r12mkgmnab(xint,cmos,isymc,idel,isymd,
     &                         iglmrhs,nglmds,lunitr12,filer12,
     &                         work,lwork)
c---------------------------------------------------------------------
c     purpose: make g^MN_alpha_beta intermediate
c              M and N run over all active and inactive occupied
c              molecular orbitals
c
c     H. Fliegl, C. Haettig, spring 2004
c---------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <r12int.h>
      logical locdbg
      parameter (locdbg = .false.)
      integer icount1,isym,isym1,isym2
      integer idelta,idel,isymd,iglmrhs(8,8),nglmds(8),lunitr12,lwork,
     &        isymn,isymdn,isymab,isymbm,koff1,kscr1,kend1,lwrk1,koffg,
     &        koffc,ntotab,ntotg,idxabn,isymm,isymb,isyma,koff5,kscr2,
     &        koff3,koff4,ntota,ntotb,isymc,isymg
      integer nrhftsorb

      character*(*) filer12

#if defined (SYS_CRAY) 
      real work(*), cmos(*), xint(*), one, zero, xnrm, ddot
      real xnrm1, xnrm2
#else
      double precision work(*), cmos(*), xint(*), one, zero, xnrm, ddot
      double precision xnrm1, xnrm2
#endif
      parameter (one = 1.0d0, zero = 0.0d0)

      idelta = idel - ibas(isymd)  

c     skip if idelta not in orbital basis
      if (idelta.gt.mbas1(isymd)) return

      call qenter('r12mkgmnab')

      nrhftsorb = 0
      do isym = 1, nsym
        nrhftsorb = nrhftsorb + nrhfsa(isym)
      end do
c
c     transform first index 
c
      if (locdbg) then
        xnrm = 0.0d0
        xnrm1 = 0.0d0
        xnrm2 = 0.0d0
      end if
c
      do isymg = 1, nsym
        isymn  = muld2h(isymc,isymg) ! isymn = isymg
        isymdn = muld2h(isymd,isymn)
        isymab = isymdn
        isymbm = isymdn

        koff1  = 1
        kscr1  = koff1 + nnbst(isymab)*nrhfsa(isymn)
        kscr2  = kscr1 + nbast*nbast
        kend1  = kscr2 + nbast*nrhftsorb
        lwrk1  = lwork - kend1
        if (lwrk1.lt.0) then
          call quit('Insufficient work space in cc_r12mkgmnab') 
        end if

        koffg  = idsaog(isymg,isymd) + 1
        koffc  = iglmrhs(isymg,isymn) + 1

        ntotab = max(nnbst(isymab),1)
        ntotg  = max(nbas(isymg),1)

        call dgemm('N','N',nnbst(isymab),nrhfsa(isymn),mbas1(isymg),
     &             one,xint(koffg),ntotab,cmos(koffc),ntotg,zero,
     &             work(koff1),ntotab)
c
        if (locdbg) then
          xnrm1 = xnrm1+ddot(nnbst(isymab)*nrhfsa(isymn),
     &                       work(koff1),1,work(koff1),1)
        end if
c
        
        do n = 1, nrhfsa(isymn)
          idxabn = koff1 + nnbst(isymab)*(n-1)
c         unpack integrals as a quadratic matrix
          call ccsd_symsq(work(idxabn),isymab,work(kscr1))
c
          if (locdbg) then
            xnrm2=xnrm2+ddot(n2bst(isymab),work(kscr1),1,work(kscr1),1) 
          end if

c         transform second index to M
          koff5 = kscr2
          do isyma = 1, nsym
            isymm = muld2h(isymc,isyma) !isymm = isyma
            isymb = muld2h(isymab,isyma)

            koff3 = kscr1 + iaodis(isyma,isymb)
            koff4 = iglmrhs(isyma,isymm) + 1

            ntota = max(nbas(isyma),1)
            ntotb = max(nbas(isymb),1)

            call dgemm('T','N',nbas(isymb),nrhfsa(isymm),mbas1(isyma),
     &                 one,work(koff3),ntota,cmos(koff4),ntota,zero,
     &                 work(koff5),ntotb)
c
            if (locdbg) then
            xnrm = xnrm+ddot(nbas(isymb)*nrhfsa(isymm),work(koff5),1,
     &                                                work(koff5),1) 
            end if
c

            koff5 = koff5 + nbas(isymb)*nrhfsa(isymm)
          end do

c         write result on file
          call cc_r12whtfs(work(kscr2),idel,isymd,n,isymn,isymbm,
     &                     lunitr12,filer12,work(kend1),lwrk1)          
        end do
      end do
c
      if (locdbg) then
        write(lupri,*) 'norm^2(cmo)',ddot(nglmds,cmos,1,cmos,1)
        write(lupri,*) 'Gmbnd for idel=',idel,xnrm
        write(lupri,*) 'Gabnd for idel=',idel,xnrm1,xnrm2
      end if
c
      call qexit('r12mkgmnab')
      end
*=====================================================================*
      subroutine cc_r12whtfs(xint,idel,isymd,j,isymj,isymai,
     &                      lunitr12,filer12,work,lwork)
c----------------------------------------------------------------------
c     purpose: get half transformed coulomb integrals and write them
c              on file first delta in ccsd_aibj   
c              note: M,N run over all active and inactive occupied
c                    molecular orbitals 
c
c     H. Fliegl, C. Haettig, spring 2004
c----------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
#include <dummy.h>
       
      logical locdbg
      parameter (locdbg = .false.)
      integer lunitr12, lwork, idel, isymd, isymi, isyma
      integer isymai, isymdj, isymj, ndelj, iadr, ilen
      integer idxai, idxdj, idxai2, koff1, idelta
      integer isym,isym1,isym2,icount1,ng1bas(8),ig2bas(8,8),ig1bas(8,8)

      character*(*) filer12

#if defined (SYS_CRAY)
      real work(*), xint(*)
#else
      double precision work(*), xint(*)
#endif

      call qenter('whtfs')

c     calculate some offsets and dimensions
      do isym = 1, nsym
        ng1bas(isym) = 0
        icount1 = 0
        do isym2 = 1, nsym
          isym1 = muld2h(isym,isym2)
          ng1bas(isym) = ng1bas(isym)+mbas1(isym1)*nrhfsa(isym2)
          ig1bas(isym1,isym2) = icount1
          icount1 = icount1 + mbas1(isym1)*nrhfsa(isym2)
        end do
      end do
c
      do isym = 1, nsym
        icount1 = 0
        do isym2 = 1, nsym
          isym1 = muld2h(isym,isym2)
          ig2bas(isym1,isym2) = icount1
          icount1 = icount1 + ng1bas(isym1)*ng1bas(isym2)
        end do
      end do
c
      if (lwork.lt.ng1bas(isymai)) then
        call quit('Insufficient work space in cc_r12whtfs')
      end if

      idelta = idel - ibas(isymd)  

      koff1 = 1
      isymdj = muld2h(isymd,isymj)
      do isymi = 1, nsym
         isyma  = muld2h(isymai,isymi)
        do i = 1, nrhfsa(isymi)
            do a = 1, mbas1(isyma) 
              idxai   = ig1bas(isyma,isymi) + mbas1(isyma)*(i-1)+a 
              idxai2  = koff1 -1 + nbas(isyma)*(i-1) +a
              work(idxai) = xint(idxai2) 
            end do 
        end do 
        koff1 = koff1 + nbas(isyma)*nrhfsa(isymi)
      end do
      
      if (locdbg) then
        write(lupri,*)'idel=',idel
        write(lupri,*)'in whtfs'
        do isymi = 1, nsym
         isyma  = muld2h(isymai,isymi)
         write(lupri,*) 'symmetry block (a,i):',isyma,isymi
         call output(work(ig1bas(isyma,isymi)+1),1,mbas1(isyma),
     &    1,nrhfsa(isymi),mbas1(isyma),nrhfsa(isymi),1,lupri) 
        end do
      end if
      
      if (idelta.le.mbas1(isymd)) then
        ndelj = ig1bas(isymd,isymj) + mbas1(isymd)*(j-1) + idelta 
        iadr  = ig2bas(isymai,isymdj) + ng1bas(isymai)*(ndelj-1)+1
      else 
        call quit('There is something wrong in cc_r12whtfs'//
     &            ' with the auxiliary basis')
      end if
      ilen = ng1bas(isymai)

      call putwa2(lunitr12,filer12,work,iadr,ilen)

      call qexit('whtfs')
      end 
*=====================================================================*
      subroutine cc_r12getrints(r12bkl,idel,isymd,nr1basm,ir1basm,
     &           nr2basm,ir2basm,nrgkls,irgkls,ir1xbasm,ir2xbasm,
     &           ibasx,imatijm,nmatijm,lauxd,fconst,work,lwork)
c---------------------------------------------------------------------
c     purpose: get g^MN_gamma_delta integrals, needed for ansatz 2 
c
c     H. Fliegl, C. Haettig spring 2004
c---------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>

      logical lauxd,locdbg 
      parameter(locdbg=.false.)

      character*(*) fconst

#if defined (SYS_CRAY)
      real zero,one
#else
      double precision zero,one
#endif
      parameter (zero = 0.0d0, one = 1.0d0)

      integer nr1basm(8),ir1basm(8,8),nr2basm,ir2basm(8,8)
      integer ir2xbasm(8,8),nrgkls(8),imatijm(8,8)
      integer irgkls(8,8),ir1xbasm(8,8),nmatijm(8)
      integer idelta,idel,lwork,lwrk1,ibasx(8)
      integer lur2back,kscr2,kend2,lwrk2,iadr,ndell,ilen
      integer isymb,isymk,isymkl,isyml,idxbk,idxkl,idxbkl
      integer isymd,isymdl,isymbk,kend1 


#if defined (SYS_CRAY)
      real r12bkl(*),work(*),ddot
#else
      double precision r12bkl(*),work(*),ddot
#endif
      
      kend1 = 1

      idelta = idel - ibas(isymd)
      if (lauxd) idelta = idelta - ibasx(isymd)

      lur2back = -1
      call wopen2(lur2back,fconst,64,0)
      
      do isyml = 1, nsym
         isymdl = muld2h(isymd,isyml)
         isymbk = isymdl
         
         kscr2 = kend1 
         kend2 = kscr2 + nr1basm(isymbk)
         lwrk2 = lwork - kend2
         if (lwrk2 .lt.0) then
            write(lupri,*)'lwork',lwork
            write(lupri,*)'lwrk2,kend2,kscr2', lwrk2,kend2,kscr2
            call quit('Insufficient work space in cc_r12getrints')
         end if
c       ------------------------------------------------------------
c        read r^{beta,delta}_{k,l} from file and pack as 
c        R^{delta}(beta,kl) or  R^{delta'}(beta,kl)
c       ------------------------------------------------------------
         do l = 1, nrhfsa(isyml)
            if (idelta.le.mbas1(isymd)) then
               ndell = ir1basm(isymd,isyml) +
     &              mbas1(isymd)*(l-1)+idelta
               iadr = ir2basm(isymbk,isymdl)+
     &              nr1basm(isymbk)*(ndell-1)+1
            else
               ndell = ir1xbasm(isymd,isyml)+
     &              mbas2(isymd)*(l-1)+idelta-mbas1(isymd)
               iadr = nr2basm + ir2xbasm(isymbk,isymdl)+
     &              nr1basm(isymbk)*(ndell-1)+1

               if(locdbg) then
                  write(lupri,*)'ndell,iadr,l,idelta',ndell,iadr,
     &                 l,idelta
               end if

            end if
            ilen   = nr1basm(isymbk)
            call getwa2(lur2back,fconst,work(kscr2),iadr,ilen)
c 
            do isymb =1, nsym
               isymk = muld2h(isymbk,isymb)
               isymkl = muld2h(isymk,isyml)
               do k=1, nrhfsa(isymk)
                  do b=1, mbas1(isymb)
                     idxbk = ir1basm(isymb,isymk) + mbas1(isymb)*(k-1)+b
                     idxkl = imatijm(isymk,isyml) +nrhfsa(isymk)*(l-1)+k
                     idxbkl = irgkls(isymb,isymkl) + 
     &                    mbas1(isymb)*(idxkl-1)+b

                     if(locdbg) then
                        write(lupri,*)'idxbk,idxkl,idxbkl',
     &                       idxbk,idxkl,idxbkl
                     end if

                     r12bkl(idxbkl) = work(kscr2 -1 + idxbk)
                     
                  end do
               end do

               if (locdbg) then
                 write(lupri,'(2a,4i5)')' R^{delta}(beta,kl) for ',
     &            'delta,l,isymb,isymkl',idelta,l,isymb,isymkl
                 idxbkl = irgkls(isymb,isymkl) + 1
                 call output(r12bkl(idxbkl),1,mbas1(isymb),1,
     &             nmatijm(isymkl),mbas1(isymb),nmatijm(isymkl),1,lupri)
               end if

            end do
              
         end do
      end do
      call wclose2(lur2back,fconst,'KEEP')

      if (locdbg) then
        write(lupri,*) 'leaving cc_r12getrints... norm^2(r12bkl):',
     &     ddot(nrgkls(isymd),r12bkl,1,r12bkl,1)
      end if

      end
*=====================================================================*
      subroutine cc_r12fdvint(c1am,c2am,work,lwork,aproxr12,
     &                        fc12am,lufc12,ivec)
c-----------------------------------------------------------------------
c      purpose: Test routine to calculate V bar numerically used for
c               CC2-R12 ansatz 2
c
c      H. Fliegl, C. Haettig, winter 2004
c-----------------------------------------------------------------------
      implicit none
#include <maxorb.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>
#include <dummy.h>
#include <second.h>
#include <ccsdinp.h>
#include <cclr.h>
#include <r12int.h>
#include <ccr12int.h>
      logical omegorsv, testvajtkl, testvbar, testrbar 
      parameter (testvajtkl = .false., testvbar = .false., 
     &           testrbar = .false.)

      integer lwork,kvintp,kvintm,kend1,krho1,krho2,kt1am,
     &        kt2am,lwrk1,iopt,luvijkl,isymij,isymkl,kvajkl,luvajkl
      integer kvsym,kvsing,kvtrip,nrhftria,krp,krm,ngamsm(8),
     &        imatijm(8,8),nmatijm(8),igamsm(8,8),isym,isym1,isym2,
     &        icount1,icount4,krsym,lufc12,krho1p,krho2p,nrho2,ivec
      integer khintm,khintp
#if defined (SYS_CRAY)
      real work(*),c1am(*),c2am(*)
      real xhalf, xmtwo, delta
#else
      double precision work(*),c1am(*),ddot,c2am(*)
      double precision xhalf, xmtwo, delta, deltai
#endif 
      parameter (xhalf = 0.5d00,xmtwo = -2.0d00, 
     &           delta = 1.0d-7 , deltai = 1.0d00/delta)
      character*10 model
      character*(*) aproxr12,fc12am

      integer index
      index(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j
c
      call qenter('r12fdvint')
      do isym = 1, nsym
        icount4 = 0
        do isym2 = 1, nsym
          isym1 = muld2h(isym,isym2)
          imatijm(isym1,isym2) = icount4
          icount4 = icount4 + nrhfsa(isym1)*nrhfsa(isym2)
        end do
        nmatijm(isym) = icount4
      end do
      do isym = 1, nsym
        icount1 = 0
        do isym2 = 1, nsym
          isym1 = muld2h(isym,isym2)
          igamsm(isym1,isym2)  = icount1
          icount1 = icount1 + nmatijm(isym1)*nmatkl(isym2)
        end do
        ngamsm(isym)  = icount1
      end do

      if (omegsq) call quit('cclr_fdvint INCOMPATIBLE WITH omegsq.')
      if (lcor.or.lsec) 
     &   call quit('cclr_fdvint INCOMPATIBLE WITH fcore AND fsecon.')
c
       call around( 'in CCLR_FDVINT: making finite diff. v')
       omegorsv = omegor
       omegor   =.true.
c
c     isymtr = 1 ! in cclr.h 

      krho1 = 1
      krho2 = krho1   + nt1amx
      kt1am = krho2   + max(nt2amx,nt2ao(isymop),2*nt2ort(isymop))
      kt2am = kt1am   + nt1amx
      kvintm = kt2am  + max(nt2sq(isymop),nt2r12(1))
      kvintp = kvintm + ngamsq(isymtr)
      kvajkl = kvintp + ngamsq(isymtr)
      kvsym  = kvajkl + nvajkl(isymtr)
      krm  = kvsym + ngamsq(isymtr)
      krp  = krm + ngamsm(isymtr)
      krsym  = krp + ngamsm(isymtr)
      krho1p  = krsym + ngamsm(isymtr)
      krho2p  = krho1p + nt1amx
      khintm   = krho2p + max(nt2amx,nt2ao(isymop),2*nt2ort(isymop))
      khintp  = khintm + nt1amx
      kend1 = khintp + nt1amx
      lwrk1 = lwork   - kend1
      if (lwrk1.lt.0) then
         write (lupri,*) 'available: lwork   =  ',lwork
         write (lupri,*) 'needed (at least)  =  ',kend1
         call quit('insufficient work space in cclr_fdvint ')
      endif
      nrho2 = max(nt2amx,nt2ao(isymop),2*nt2ort(isymop))

      if (testvajtkl) then
c       calculate V^alpha_j_ kl
        call dzero(work(kt1am),nt1am(isymtr))
        rspim = .true.
        call ccrhsn(work(krho1),work(krho2),work(kt1am),
     &              work(kt2am),
     *              work(kend1),lwrk1,aproxr12)
      
c       read Valpajtildekl
        luvajkl = -1
        call gpopen(luvajkl,fvajtkl,'unknown',' ','unformatted',
     &              idummy,.false.)
        read(luvajkl)(work(kvajkl-1+i),i=1,nvajkl(isymtr))
        call gpclose(luvajkl,'KEEP')

        write(lupri,*)'Vajtkl out of fdvint:', 
     &                (work(kvajkl-1+i),i=1,nvajkl(isymtr))
        end if
      
c     read the CC reference amplitudes From disk.
      iopt = 3
      call cc_rdrsp('R0',0,1,iopt,model,work(kt1am),work(kt2am))
 
      rspim = .false.

c---------------------------------------------
c     Calculate V bar by finite difference.
c---------------------------------------------

c     Compute V(t - 0.5 delta * C)
      call daxpy(nt1amx,-0.5d0*delta,c1am,1,work(kt1am),1)
      call ccrhsn(work(krho1),work(krho2),work(kt1am),
     &            work(kt2am),
     *            work(kend1),lwrk1,aproxr12)

c     Read V(t - 0.5 delta * C) into WORK(KVINTM)
      luvijkl = -1
      call gpopen(luvijkl,fvijkl,'unknown',' ','unformatted',
     &            idummy,.false.)
      read(luvijkl)(work(kvintm-1+i),i=1,ngamsq(isymtr))
      call gpclose(luvijkl,'KEEP')

c     Read HP or IP TERM 'CCR12HTST' or 'CCR12ITST'
c     luvijkl = -1
c     call gpopen(luvijkl,'CCR12ITST','unknown',' ','unformatted',
c    &            idummy,.false.)
c     read(luvijkl)(work(khintm-1+i),i=1,nt1amx)
c     call gpclose(luvijkl,'KEEP')

      if (testrbar) then
c       Read R(t - 0.5 delta * C) into WORK(KRM)
        luvijkl = -1
        call gpopen(luvijkl,'RMTNTKL','unknown',' ','unformatted',
     &              idummy,.false.)
        read(luvijkl)(work(krm-1+i),i=1,ngamsm(isymtr))
        call gpclose(luvijkl,'KEEP')
      end if
c     Compute V(t + 0.5 delta * C)
      call daxpy(nt1amx,delta,c1am,1,work(kt1am),1)
      call ccrhsn(work(krho1p),work(krho2p),work(kt1am),
     &            work(kt2am),
     *            work(kend1),lwrk1,aproxr12)

c     Read V(t + 0.5 delta * C) into WORK(KVINTP)
      luvijkl = -1
      call gpopen(luvijkl,fvijkl,'unknown',' ','unformatted',
     &            idummy,.false.)
      read(luvijkl)(work(kvintp-1+i),i=1,ngamsq(isymtr))
      call gpclose(luvijkl,'KEEP')

c     Read HP or IP TERM
c     luvijkl = -1
c     call gpopen(luvijkl,'CCR12ITST','unknown',' ','unformatted',
c    &            idummy,.false.)
c     read(luvijkl)(work(khintp-1+i),i=1,nt1amx)
c     call gpclose(luvijkl,'KEEP')

      if (testrbar) then
c       Read R(t - 0.5 delta * C) into WORK(KRP)
        luvijkl = -1
        call gpopen(luvijkl,'RMTNTKL','unknown',' ','unformatted',
     &              idummy,.false.)
        read(luvijkl)(work(krp-1+i),i=1,ngamsm(isymtr))
        call gpclose(luvijkl,'KEEP')
      end if

c     Restore t amplitudes
      call daxpy(nt1amx,-0.5d0*delta,c1am,1,work(kt1am),1)

c     Calculate [V(t+0.5*delta*C)-V(t-0.5*delta*C)]/delta
      call daxpy(ngamsq(1),-1.0d00,work(kvintm),1,work(kvintp),1)
      call dscal(ngamsq(1),deltai,work(kvintp),1)

c     Calculate [RHO1(t+0.5*delta*C)-RHO1(t-0.5*delta*C)]/delta
      call daxpy(nt1amx,-1.0d00,work(krho1),1,work(krho1p),1)
      call dscal(nt1amx,deltai,work(krho1p),1)

c     Calculate [RHO2(t+0.5*delta*C)-RHO2(t-0.5*delta*C)]/delta
      call daxpy(nrho2,-1.0d00,work(krho2),1,work(krho2p),1)
      call dscal(nrho2,deltai,work(krho2p),1)

c     Calculate finite diff HP BAR
c     call daxpy(nt1amx,-1.0d00,work(khintm),1,work(khintp),1)
c     call dscal(nt1amx,deltai,work(khintp),1)

c     write(lupri,*)'IBAR numerically'
c     do isymij = 1, nsym
c      isymkl = muld2h(isymij,isymtr)
c      call output(work(khintp+it1am(isymij,isymkl)),1,
c    &      nvir(isymij),1,nrhf(isymkl),nvir(isymij),
c    &      nrhf(isymkl),1,lupri)
c     end do
      

c     Print dOmega/dt_1 Result:
      write(lupri,*)'numerical: RHO1 and RHO2:'
      call cc_prp(work(krho1p),work(krho2p),1,1,1)
      write(lupri,*)'numerical: RHOR12:'
      do isymij = 1, nsym
       isymkl = muld2h(isymij,isymtr)
       call output(work(kvintp+igamsq(isymij,isymkl)),1,
     &      nmatij(isymij),1,nmatij(isymkl),nmatij(isymij),
     &      nmatij(isymkl),1,lupri)
      end do

      if (testrbar) then
c       Calculate [R(t+0.5*delta*C)-R(t-0.5*delta*C)]/delta
        call daxpy(ngamsm(1),-1.0d00,work(krm),1,work(krp),1)
        call dscal(ngamsm(1),deltai,work(krp),1)
      end if

      if (testvbar) then
c     symmetrize V bar 
        call symV(work(kvintp),isymtr,work(kvsym),
     &            nrhf,imatij,igamsq,nmatij,work(kend1),lwrk1)

c     Print the result:
        write(lupri,*)'norm^2 of V bar numerically:', 
     &    ddot(ngamsq(isymtr),work(kvsym),1,work(kvsym),1)
        do isymij = 1, nsym
         isymkl = muld2h(isymij,isymtr)
         write(lupri,*) 'isymij,isymkl:',isymij,isymkl
         write(lupri,*)'CC2-R12 <V bar numerically>'
         call output(work(kvsym+igamsq(isymij,isymkl)),1,
     &        nmatij(isymij),1,nmatij(isymkl),nmatij(isymij),
     &        nmatij(isymkl),1,lupri)
        end do
c       split into triplett and singulett
        nrhftria = nrhftb*(nrhftb+1)/2
        kvsing = kend1
        kvtrip = kvsing + nrhftria*nrhftria
        kend1  = kvtrip + nrhftria*nrhftria
        write(lupri,*)'Vsing and Vtrip out of fdvint'
        call cc_r12vunpack(work(kvsym),isymtr,work(kvsing),
     &                     work(kvtrip),.true.,nrhfb,nrhf)
       end if

       if (testrbar) then
c        Print result for R:
c        symmetrize R:
         call symV(work(krp),isymtr,work(krsym),
     &             nrhfsa,imatijm,igamsm,nmatijm,work(kend1),lwrk1)
         write(lupri,*)'norm^2 of R bar numerically:', 
     &                ddot(ngamsm(isymtr),work(krsym),1,work(krsym),1)
         do isymij = 1, nsym
           isymkl = muld2h(isymij,isymtr)
           write(lupri,*) 'isymij,isymkl:',isymij,isymkl
           write(lupri,*)'CC2-R12 <R bar numerically>'
           call output(work(krsym+igamsm(isymij,isymkl)),1,
     &          nmatijm(isymij),1,nmatij(isymkl),nmatijm(isymij),
     &          nmatij(isymkl),1,lupri)
         end do
       end if

      omegor = omegorsv

      write(lupri,*)'leaving cc_r12fdvint'
      call qexit('r12fdvint')
      end
*=====================================================================*
      subroutine symV(vijkl,isymv,vsym,ndim1,ioff1,ioff2,ndim2,
     &                work,lwork)
c---------------------------------------------------------------------
c     purpose: symmetrize Vijkl = 1/2 (Vijkl + Vjilk)
c
c     H. Fliegl, C. Haettig, winter 2004
c
c     ndim1 --> nrhf / nrhfs
c     ioff1 --> imatij / imatijm
c     ioff2 --> igamsq / igamsm
c     ndim2 --> nmatij / nmatijm
c      
c---------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
    
      integer lwork,isymv,lwrk1,kend1,isymkl,isymij,isymk,
     &        isyml,isymi,isymj,idxkl,idxlk,idxij,idxji,idxijkl,idxjilk
      integer ndim1(8),ioff1(8,8),ioff2(8,8),ndim2(8)

#if defined (SYS_CRAY)
      real vijkl(*),work(*),half,vsym(*)
#else
      double precision vijkl(*),work(*),half,vsym(*)
#endif
      parameter ( half = 0.5d0 )

      call qenter('symV')

      do isymkl = 1, nsym
        isymij = muld2h(isymv,isymkl)  
        do isymk = 1, nsym
          isyml = muld2h(isymkl,isymk)
          do k = 1, nrhfb(isymk)
            do l = 1, nrhfb(isyml)
              idxkl = imatkl(isymk,isyml)+nrhfb(isymk)*(l-1)+k
              idxlk = imatkl(isyml,isymk)+nrhfb(isyml)*(k-1)+l
              do isymi = 1, nsym
                isymj = muld2h(isymij,isymi)
                do i = 1, ndim1(isymi)
                  do j = 1, ndim1(isymj)
                    idxij = ioff1(isymi,isymj)+ndim1(isymi)*(j-1)+i
                    idxji = ioff1(isymj,isymi)+ndim1(isymj)*(i-1)+j
                    idxijkl = ioff2(isymij,isymkl)+
     &                        ndim2(isymij)*(idxkl-1)+idxij
                    idxjilk = ioff2(isymij,isymkl)+
     &                        ndim2(isymij)*(idxlk-1)+idxji
                    vsym(idxijkl)=half*(vijkl(idxijkl)+vijkl(idxjilk))
                  end do
                end do
              end do  
            end do
          end do
        end do
      end do 

      call qexit('symV')
      end 
*=====================================================================*
      subroutine cc_r12stabmat(amat,maxred,len1,len2,len3,work,lwork)
c----------------------------------------------------------------------
c     purpose: get stability matrix and calculate the corresponding
c              eigenvalues
c             
c              e  C^T
c              C  B
c
c     H. Fliegl, C. Haettig
c----------------------------------------------------------------------
      implicit none
#include <priunit.h>

      integer maxred,len1,len2,len3,lwork,i,j,idxi,idxj,ierr
      integer kevr,kevi,kevec,kiv1,kfv1,kend1,lwrk1,matz,nm
#if defined (SYS_CRAY)
      real amat(maxred,maxred),work(*),stabmat(len2+len3,len2+len3)
      real zero,stabsave(len2+len3,len2+len3)
#else
      double precision amat(maxred,maxred),work(*),
     &                 stabmat(len2+len3,len2+len3),zero
      double precision stabsave(len2+len3,len2+len3)
#endif
      parameter(zero = 0.0d0)

      call qenter('stabmat')
c
      nm = len2+len3 

      kevr = 1
      kevi = kevr + nm
      kevec = kevi + nm
      kiv1  = kevec + nm*nm
      kfv1  = kiv1 + nm*nm
      kend1 = kfv1 + nm*nm
      lwrk1 = lwork - kend1
      if (lwrk1.lt.0) then
        call quit('Insufficient work space in stabmat')
      end if

c     get stability matrix out of jacobian 
      do i = len1+1, len1+len2+len3
        idxi = i-len1
        do j = len1+1, len1+len2+len3
          idxj = j-len1
          stabmat(idxi,idxj) = amat(i,j) 
        end do        
      end do

      write(lupri,*)'STABILITY-MATRIX'
      call output(stabmat,1,len2+len3,1,len2+len3,len2+len3,
     &              len2+len3,1,lupri)
c     save matrix
      call dcopy(nm*nm,stabmat,1,stabsave,1)

c     get eigenvalues of stability matrix
      matz = 1 ! get eigenvalues and eigenvectors
      call rg(nm,nm,stabmat,work(kevr),work(kevi),matz,work(kevec),
     &        work(kiv1),work(kfv1),ierr)
c 
      do i = 1, nm
       if (work(kevr-1+i).lt.zero) then
        write(lupri,*)
     &        'WARNING: NEGATIVE EIGENVALUES OF STABILITY-MATRIX'
       end if
      end do
c
      write(lupri,*)'Eigenvalues of STABILITY-MATRIX'
      call output(work(kevr),1,len2+len3,1,1,len2+len3,
     &              len2+len3,1,lupri)
c
      call dcopy(nm*nm,stabsave,1,stabmat,1)
c     set nondiagonal blocks to zero and recalculate the eigenvalues
      do i = len2+1, nm
        do j = 1, len2
          stabmat(i,j) = zero
          stabmat(j,i) = zero
        end do
      end do
      write(lupri,*)'STABILITY-MATRIX without C-MAT'
      call output(stabmat,1,len2+len3,1,len2+len3,len2+len3,
     &              len2+len3,1,lupri)

      call rg(nm,nm,stabmat,work(kevr),work(kevi),matz,work(kevec),
     &        work(kiv1),work(kfv1),ierr)
c 
      do i = 1, nm
       if (work(kevr-1+i).lt.zero) then
        write(lupri,*)
     &         'WARNING: NEGATIVE EIGENVALUES OF STABILITY-MATRIX'
       end if
      end do
c
      write(lupri,*)'Eigenvalues of STABILITY-MATRIX without C-MAT'
      call output(work(kevr),1,len2+len3,1,1,len2+len3,
     &              len2+len3,1,lupri)
c

      call qexit('stabmat')

      end 
*=====================================================================*
*=====================================================================*
      subroutine cc_r12getramnp(xpmqn,framnp,oneaux,norb1,norb2,nrhf1,
     &                          nrhfs1,ih1am,ih2am,work,lwork)
c---------------------------------------------------------------------
c     purpose: get r_a,mn^q' out of r^mn_pq 
c              and put on file
c
c     C. Neiss, spring 2006
c---------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
      logical locdbg,oneaux
      parameter (locdbg = .false.)
      integer ih2am(8,8), ih1am(8,8), nvircc(8),ldr12(8),
     &        norb1(8),norb2(8),nrhf1(8),nrhfs1(8),
     &        nmatkl(8),imatkl(8,8),nramn(8),iramn(8,8),
     &        nramnq(8),iramnq(8,8)
      integer isym,isym1,isym2,isymp,isymq,isympmn,isymmn,isymm,isymn,
     &        isympm,isymqn
      integer lwork,kend1,lwrk1,kramn,lunit,len,iadr,i,j,a,p,q,m,n
      integer nmn,npm,nqn,idxpmqn,idxamn,index
      character*(*) framnp
#if defined (SYS_CRAY)
      real xpmqn(*),work(*),ddot
#else
      double precision xpmqn(*),work(*),ddot
#endif
      index(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j

      call qenter('cc_r12getramnp')

c     calculate some offsets and dimensions
c     calculate number of active virtual orbitals in CC
c     nvirfr = frozen virtual orbitals        
      do isym = 1, nsym
        nvircc(isym) = norb1(isym) - nrhfsa(isym) - nvirfr(isym)
      end do

      if (locdbg) then
        write(lupri,*) 'Entered cc_r12getramnp'
        write(lupri,*) 'isym nrhfs1 nrhf1 nvircc norb1'
        do isym = 1, nsym
          write(lupri,*) isym,nrhfs1(isym),nrhf1(isym),
     &        nvircc(isym),norb1(isym)
        end do 
      end if

      do isym = 1, nsym
        nmatkl(isym) = 0
        do isym2 = 1, nsym
          isym1 = muld2h(isym,isym2)
          imatkl(isym1,isym2) = nmatkl(isym) 
          nmatkl(isym) = nmatkl(isym) + nrhf1(isym1)*nrhf1(isym2)
        end do
      end do

      do isym = 1, nsym
        nramn(isym) = 0
        do isym2 = 1, nsym
          isym1 = muld2h(isym,isym2)
          iramn(isym1,isym2) = nramn(isym)
          nramn(isym) = nramn(isym) + nvircc(isym1)*nmatkl(isym2)
        end do
      end do

      do isym = 1, nsym
        nramnq(isym) = 0
        do isym2 = 1, nsym
          isym1 = muld2h(isym,isym2)
          iramnq(isym1,isym2) = nramnq(isym)
          nramnq(isym) = nramnq(isym) + nramn(isym1)*norb2(isym2)
        end do
      end do
     
      do isym = 1, nsym
        if (oneaux) then
          ldr12(isym) = norb1(isym)
        else
          ldr12(isym) = nvir(isym)
        end if
      end do

c
c    ---------------------------------------------------------------
c     open output file
c    ---------------------------------------------------------------
      lunit = -1
      call wopen2(lunit,framnp,64,0)
c
c    ---------------------------------------------------------------
c     reorder integrals
c    ---------------------------------------------------------------
      do isymq = 1, nsym
        isympmn = isymq
        do q = norb1(isymq)+1, norb1(isymq)+norb2(isymq)
          len = nramn(isympmn)
          kramn = 1
          kend1 = kramn + len
          lwrk1 = lwork - kend1
          if (lwrk1.lt.0) then
            call quit('Insufficient memory in CC_R12GETRAMNP')
          end if
          do isymp = 1, nsym
            isymmn = muld2h(isympmn,isymp)
c
            do isymm = 1, nsym
              isymn = muld2h(isymmn,isymm)
              isympm = muld2h(isymp,isymm)
              isymqn = muld2h(isymq,isymn)
              do n = 1, nrhf1(isymn)
                do m = 1, nrhf1(isymm)
                  do p = nrhfs1(isymp)+1, nrhfs1(isymp)+nvircc(isymp)
                    a = p - nrhfs1(isymp)
                    npm = ih1am(isymp,isymm)+ldr12(isymp)*(m-1)+p
                    nqn = ih1am(isymq,isymn)+ldr12(isymq)*(n-1)+q
                    nmn = imatkl(isymm,isymn)+nrhf1(isymm)*(n-1)+m
                    idxpmqn = ih2am(isympm,isymqn) + index(npm,nqn)
                    idxamn = iramn(isymp,isymmn) + 
     &                       nvircc(isymp)*(nmn-1)+a
                    work(kramn+idxamn-1) = xpmqn(idxpmqn)
                  end do
                end do
              end do
            end do
          end do
          iadr = iramnq(isympmn,isymq) + 
     &           nramn(isympmn)*(q-norb1(isymq)-1) + 1 
          call putwa2(lunit,framnp,work(kramn),iadr,len)
          if (locdbg) then
            write(lupri,*) 'isymq, q, iadr= ',isymq,q,iadr
            write(lupri,*) 'Norm^2:', ddot(len,work(kramn),1,
     &                                     work(kramn),1) 
            do isym2 = 1, nsym
              isym1 = muld2h(isympmn,isym2)
              write(lupri,*) 'Symmetry block: ',isym1,isym2
              call output(work(kramn+iramn(isym1,isym2)),
     &                    1,nvircc(isym1),1,nmatkl(isym2),
     &                    nvircc(isym1),nmatkl(isym2),1,lupri)
            end do
          end if
c
        end do
      end do
c
c    ---------------------------------------------------------------
c     close file
c    ---------------------------------------------------------------
      call wclose2(lunit,framnp,'KEEP')
c
      if (locdbg) then
        write(lupri,*) 'Leaving cc_r12getramnp'
      end if
c
      call qexit('cc_r12getramnp')
      return
      end
*=====================================================================*
