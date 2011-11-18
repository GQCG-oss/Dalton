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
c*=====================================================================*
      subroutine r12mkvamkl(filback,gabjd,tabjd,vajkl,cmo,isymcmo,
     &                      xlamdhs,xlamdps,xint,
     &                      tint,idel,isymd,isymj,isymab,isymg,scr1,
     &                      ibasx,iglmrhs,nglmds,work,lwork)
*----------------------------------------------------------------------*
c     Purpose: calculation of V^{alpha,j}_{k,l} 
c     
c     where:     alpha      AOs in orbital basis
c                j          active occupied MO
c                k,l        r12 MOs
c
c     H. Fliegl, C. Haettig spring 2003
c
c     C. Neiss 2004/2005: adapted for non-total-symm. V (Ansatz 1)
c                         !Ansatz 2 still needs to be checked!
c     isymcmo    Symmetry of cmo-matrix (needed if cmo=lambda_bar)
c
c     adapted for ansatz 3, Christof Haettig, spring 2005.
c
c     adapted for CABS, Christian Neiss, winter 2005/2006
*----------------------------------------------------------------------*
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
#include <ccr12int.h>
#include <dummy.h>

      logical locdbg
      parameter(locdbg=.false.)
#if defined (SYS_CRAY)
      real zero,one,two
#else
      double precision zero,one,two
#endif
      parameter (zero = 0.0d0, one = 1.0d0, two = 2.0d0)

      character*(*) filback
      integer nr1orb(8),nr1bas(8),nr1xbas(8),n2bst1(8),nrgkl(8)
      integer ir1orb(8,8),ir1bas(8,8),ir2bas(8,8),ir1xbas(8,8)
      integer ir2xbas(8,8),iaodis1(8,8),irgkl(8,8),nglmds(8)
      integer igabj(8),igbi(8),igbgi(8),ibasx(8),iglmrhs(8,8)
      integer nr1xorb(8),ir1xorb(8,8),nrxgkl(8),irxgkl(8,8)
      integer isymd,isymj,isymgk,isymg,isymk,isymkl
      integer isymab,isyma,isymb,isymaj,isymdl,isyml,isymbk
      integer isymi,isymgi,isymag,isymbd,isymbj,isymv,isymcmo
      integer idxkl,idxaj,idxabj,idxbk,idxbkl,idxab,idxgi,idxbi
      integer idxbgi,idxlk,idxblk
      integer ilen,iadr,ndell,koff1,idel,idelta 
      integer ngabj,kr2bkl,lwork,kend1,kend0,ngbi,kgab,kgbi
      integer ntota,ngbgi,kgbgi,ntotgi,iansatz
      integer lwrk1,lur2back,kscr2,kend2,lwrk2,lwrk0,krblk
      integer koffg,koffr,koffv,ntotb,ntotaj,nr2bas,kgabj,koffc
      integer idxij,idxji,kcmoa,idxba 
      integer icou2,iglmro(8,8) 
      integer krtf,kgtf,isymm,isymdm,ntotg,ntotv,isymmj,ntotc, idxd
      integer krdpmlk,kgdpmja,isymmk,idxmlk,idxmkl,idxja,
     &        isymdj,idxmja,idxmaj,isymma, kend3,kend4,
     &        lwrk4, igmaj(8), ngmaj, ngmgi, kcmo,kgmaj,ntotm,igmgi(8)
      integer koffrtf,imaklm(8,8),nmaklm(8),isymklm,icount1,isymij
      integer nalphaj(8),ialphaj(8,8)
#if defined (SYS_CRAY)
      real gabjd(*),tabjd(*),vajkl(*),cmo(*),xint(*),work(*)
      real scr1(*), ddot, nbjd(*), tint(*),xlamdhs(*), xlamdps(*) 
      real factor, xrnmr, xnrmg
#else
      double precision gabjd(*),tabjd(*),vajkl(*),cmo(*),xint(*),work(*)
      double precision scr1(*),ddot,tint(*),xlamdhs(*),xlamdps(*)
      double precision factor, xnrmr, xnrmg,x1,x2,x3
#endif
      save xnrmr, xnrmg,x1,x2,x3
      data xnrmr /0.0d0/
      data xnrmg /0.0d0/
      data x1 /0.0d0/
      data x2 /0.0d0/
      data x3 /0.0d0/


      idelta = idel - ibas(isymd)

c    ----------------------------------------------------
c     for ansatz 2 skip loop over deltas in orbital basis 
c    ----------------------------------------------------
      if (ianr12.eq.2 .and. idelta.le.mbas1(isymd) .and. .not.r12cbs) 
c     if (ianr12.eq.2 .and. idelta.le.mbas1(isymd) ) 
     &  return

      ! map ianr12 --> iansatz to get simpler if statements
      ! for ansatz 3 which has common terms with ansatz 2
      if (ianr12.eq.1 .or. ianr12.eq.2) then
        iansatz = ianr12
        if (idelta.le.mbas1(isymd)) then
           factor = +one
           if (R12CBS) factor = -one
        else 
           factor = -one
        end if
      else if (ianr12.eq.3) then
        if (idelta.le.mbas1(isymd)) then
          ! delta is from orbital basis
          ! -> sum_{delta beta}, as in ansatz 1,
          !    but with AO indeces at R projected onto virtuals
          iansatz = 1
          factor  = -one
        else
          ! delta is from auxiliary basis 
          ! -> sum_{delta' m} as in ansatz 2
          iansatz = 2
          factor  = -one
        end if
      else
        call quit('Unknown ansatz in r12mkvamkl.')
      end if


      call qenter('r12mkvamkl')

      if (locdbg .and. iansatz.eq.2) then
        write(lupri,*)'xlamdhs in mkvamkl:', 
     &                ddot(nglmds(1),xlamdhs,1,xlamdhs,1)
        write(lupri,*)'xlamdps in mkvamkl:', 
     &                ddot(nglmds(1),xlamdps,1,xlamdps,1)
      end if

c    ------------------------------------------------
c     calculate some symmetry offsets and dimensions
c    ------------------------------------------------
      call cc_r12offset(nr1orb,nr1xorb,nr1bas,nr1xbas,nr2bas,
     &     nrgkl,nrxgkl,n2bst1,ir1orb,ir1xorb,ir1bas,ir1xbas,ir2bas,
     &     ir2xbas,irgkl,irxgkl,iaodis1,nalphaj,ialphaj)

      do isymaj = 1,nsym
         icou2 = 0
         do isyma = 1,nsym
            isymi = muld2h(isymaj,isyma)
            iglmro(isymi,isyma) = icou2
            icou2 = icou2 + nbas(isymi)*(norb1(isyma)-nrhffr(isyma))
         enddo
      enddo
     
      
      ngabj = 0
      ngmaj = 0
      ngmgi = 0
      do isymb = 1, nsym
         isyma = muld2h(isymb,isymab)
         isymi = muld2h(isymcmo,isyma)
         isymaj = muld2h(isyma,isymj)
         isymgi = muld2h(isymg,isymi)

         igabj(isymaj) = ngabj
         igmaj(isymaj) = ngmaj
         igmgi(isymgi) = ngmgi

         ngabj  = ngabj + mbas1(isymb)*mbas1(isyma)*nrhfa(isymj)
         ngmaj  = ngmaj + nrhfsa(isymb)*mbas1(isyma)*nrhfa(isymj)
         ngmgi  = ngmgi + nrhfsa(isymb)*mbas1(isymg)*nrhfa(isymi)
      end do
      
c    
      do isymklm = 1, nsym
        icount1 = 0
        do isymm = 1, nsym
          isymkl = muld2h(isymklm,isymm)
          imaklm(isymkl,isymm) = icount1
          icount1 = icount1 + nmatkl(isymkl)*nrhfsa(isymm)
        end do
        nmaklm(isymklm) = icount1
      end do
c
      kr2bkl = 1
      kend0  = kr2bkl + nrgkl(isymd)
      if (iansatz.eq.2) then
        krtf = kend0
        kend0 = krtf + nmaklm(isymd)   
      end if 
      lwrk0  = lwork  - kend0
      
      if (iansatz.eq.1) then 
        kgabj = kend0
        kend1 = kgabj + ngabj
      else if (iansatz.eq.2) then
        kgtf  = kend0 
        kend1 = kgtf + ngmaj
      end if  

      lwrk1 = lwork - kend1
      if (lwrk1 .lt.0) then
         call quit('Insufficient work space in r12mkvamkl')
      end if
c
      ! symmetry of V is equal to isymcmo
      isymv = isymcmo
c
      if (locdbg) then
        write(lupri,*)'nalphaj: ',(nalphaj(i),i=1,nsym)
        write(lupri,*)'Vajkl entered mkvamkl'
        write(lupri,*)'norm^2:',ddot(nvajkl(isymv),vajkl,1,vajkl,1)
c       write(lupri,*)(vajkl(i),i=1,nvajkl(isymv))
        do isymkl = 1, nsym
           isymaj = muld2h(isymv,isymkl)
           write(lupri,*) 'Block isymaj, isymkl: ',isymaj,isymkl
           if ((nalphaj(isymaj).eq.0).or.(nmatkl(isymkl).eq.0)) then
            write(lupri,*) 'This block is empty'
           else
            call output(vajkl(1+ivajkl(isymaj,isymkl)),1,
     &                  nalphaj(isymaj),1,nmatkl(isymkl),
     &                  nalphaj(isymaj),nmatkl(isymkl),1,lupri)
           end if
         end do
      end if

c     -------------------------------------
c     get R12 integrals R^{delta}(beta,kl)
c     -------------------------------------
      call cc_r12getrint(work(kr2bkl),idel,isymd,nr1bas,ir1bas,
     &     nr2bas,ir2bas,nrgkl,irgkl,ir1xbas,ir2xbas,
     &     nrhfb,nmatkl,imatkl,
     &     ibasx,.false.,.false.,filback,work(kend1),lwrk1)
      
      if (locdbg) then
        write(lupri,'(a,l5,i5,/,8i7)')
     &           'idel,ibasx:',idel,ibasx
        write(lupri,*)'Rdp_a,kl in mkvamkl',idel,
     &    ddot(nrgkl(isymd),work(kr2bkl),1,work(kr2bkl),1)
c       write(lupri,*) (work(kr2bkl+i-1),i=1,nrgkl(isymd))
        do isymb = 1, nsym
          isymkl = muld2h(isymb,isymd)
          write(lupri,*) 'isymb,isymkl:',isymb,isymkl
          call output(work(kr2bkl+irgkl(isymb,isymkl)),1,mbas1(isymb),
     &                1,nmatkl(isymkl),
     &                mbas1(isymb),nmatkl(isymkl),1,lupri)
        end do
      end if

c     for ansatz 2 transform beta index to all occupied (m)
      if (iansatz.eq.2) then
        do isymb = 1, nsym
          isymm  = isymb
          isymkl = muld2h(isymb,isymd)
          koffc  = 1+iglmrhs(isymb,isymm) 
          koffr  = kr2bkl+irgkl(isymb,isymkl)
          koffrtf = krtf+imaklm(isymkl,isymm)  !for nrhfs(isymm)
      
          ntotb = max(mbas1(isymb),1)
          ntotc = max(nbas(isymb),1)
          ntotm = max(nrhfsa(isymm),1)

          call dgemm('T','N',nrhfsa(isymm),nmatkl(isymkl),mbas1(isymb),
     &               one,xlamdhs(koffc),ntotc,work(koffr),ntotb,zero,
     &               work(koffrtf),ntotm)

          x1 = x1 + ddot(nrhfsa(isymm)*mbas1(isymb),
     &                   xlamdhs(koffc),1,xlamdhs(koffc),1)
          x2 = x2 + ddot(nrgkl(isymd),work(koffr),1,work(koffr),1)
        end do


        if (locdbg) then
          write(lupri,*)'nmaklm(isymd) =', nmaklm(isymd)
          write(lupri,*)'Norm R_m,kl in mkvamkl',idel,
     &      ddot(nmaklm(isymd),work(krtf),1,work(krtf),1)
          x3 = x3 + ddot(nmaklm(isymd),work(krtf),1,work(krtf),1)
          write(lupri,*)'R_m,kl', (work(krtf-1+i),i=1,nmaklm(isymd))
        end if
      end if

c     ----------------------------------------------------------------
c      for ansatz 1: get g_{alpha beta j delta}, i.e. unpack 
c                    triangular matrix g to quadratic matrix (scr1) 
c                    and remove auxiliary basis functions alpha,beta 
c      for ansatz 2: transform beta index to occupied
c     ----------------------------------------------------------------
      do j=1, nrhfa(isymj)
         koff1= nnbst(isymab)*(j-1) +1
c Elena
         if (u21int) then
           call ccsd_symsq(gabjd(koff1),isymab,scr1)
           call ccsd_asymsq(tabjd(koff1),isymab,scr1,1,0)
c Elena
         else
           if (fnvajkl .eq. 'CCR12BIJAL') then
              call ccsd_symsqt(gabjd(koff1),isymab,scr1)
           else
              call  ccsd_symsq(gabjd(koff1),isymab,scr1)
           endif
         end if
         do isyma = 1, nsym
            isymb = muld2h(isymab,isyma)
            isymaj = muld2h(isyma,isymj)

            if (iansatz.eq.1) then
               do a = 1, mbas1(isyma)
                  idxaj = mbas1(isyma)*(j-1)+a
                  do b = 1, mbas1(isymb)
                     koff1 = iaodis(isyma,isymb) + nbas(isyma)*(b-1)+a
                     idxabj = igabj(isymaj)+mbas1(isymb)*(idxaj-1)+b  
                     work(kgabj-1+idxabj) = scr1(koff1)
                  end do
               end do
            else if (iansatz.eq.2) then  
c             m runs now over all occupied orbitals
              isymm = isymb
              koff1 = iaodis(isyma,isymb) + 1
              kcmo  = iglmrhs(isymb,isymm) + 1
              kgmaj = kend1
              lwrk1 = lwork - kend1
              if (lwrk1.lt.nrhfsa(isymm)*mbas1(isyma)) then
                call quit('Insufficient work space in mkvamkl')
              end if
            
              ntotm = max(nrhfsa(isymm),1)
              ntota = max(nbas(isyma),1)
              ntotb = max(nbas(isymb),1)
c             get G^dp_m,alpha j
              call dgemm('T','T',nrhfsa(isymm),mbas1(isyma),
     &                   mbas1(isymb),one,xlamdps(kcmo),ntotb,
     &                   scr1(koff1),ntota,zero,work(kgmaj),ntotm)

              ! sort into array aranged as G(m,aj)
              do a = 1, mbas1(isyma)
                 idxaj = mbas1(isyma)*(j-1)+a
                 do m = 1, nrhfsa(isymm)
                    koff1  = kgmaj-1 + nrhfsa(isymm)*(a-1) + m
                    idxmaj = igmaj(isymaj)+nrhfsa(isymm)*(idxaj-1)+m  
                    work(kgtf-1+idxmaj) = work(koff1)
                 end do
              end do
            else
              write(lupri,*)'ianr12 =', ianr12
              write(lupri,*)'iansatz =', iansatz
            
              call quit('illegal r12 ansatz')
            end if

         end do
      end do

      if (locdbg.and.iansatz.eq.1) then
        write(lupri,*) 'norm^2 orginal coulomb integrals:',
     &   ddot(nrhfa(isymj)*nnbst(isymab),gabjd,1,gabjd,1)
        write(lupri,*) 'resorted coulomb integrals: ... norm^2=',
     &   ddot(ngabj,work(kgabj),1,work(kgabj),1)
        do isymb = 1, nsym
          isyma  = muld2h(isymab,isymb)
          isymaj = muld2h(isyma,isymj)
          write(lupri,*) 'isymj,isyma,isymb:',isymj,isyma,isymb
          call output(work(kgabj+igabj(isymaj)),1,mbas1(isymb),
     &                1,mbas1(isyma)*nrhfa(isymj),
     &                mbas1(isymb),mbas1(isyma)*nrhfa(isymj),1,lupri)
        end do
      else if (locdbg.and.iansatz.eq.2) then 
        write(lupri,*) 'Norm:G^dp_m,gamma_i:',
     &   ddot(ngmaj,work(kgtf),1,work(kgtf),1)
        write(lupri,*) 'G^dp_m,gamma_i:', (work(kgtf-1+i),i=1,ngmaj)
      end if

c     ----------------------------------------------------------------
c      contract G with R and get V 
c 
c      V_(alpha j,kl) =  G^{delta}(beta,alpha j) * R^{delta}(beta,kl) 
c     ----------------------------------------------------------------
      if (iansatz.eq.1) then
        do isymb = 1, nsym
           isyma = muld2h(isymab,isymb)
           isymaj = muld2h(isyma,isymj)
           isymbd = muld2h(isymb,isymd)
           isymkl = isymbd
           if (isymv.ne.muld2h(isymkl,isymaj)) 
     &       call quit('Symmetry error in R12MKVAMKL')
           
           koffg = kgabj + igabj(isymaj)
           koffr = kr2bkl + irgkl(isymb,isymkl)
           koffv = 1 + ivajkl(isymaj,isymkl) + ialphaj(isyma,isymj)
           
           ntotb  = max(1,mbas1(isymb))
           ntotaj = max(1,nalphaj(isymaj))
celena
C          if (idelta.le.mbas1(isymd)) then 
C             if (R12CBS) then
C                factor = - one
C             else
C                factor = one
C             endif
C          else 
C             factor = - one
C          end if

           if (locdbg) then
             write(lupri,*)'d: V'
             call output(vajkl(koffv),1,mbas1(isyma)*nrhfa(isymj),1,
     &          nmatkl(isymkl),nalphaj(isymaj),
     &          nmatkl(isymkl),1,lupri)
           end if

           call dgemm('T','N',mbas1(isyma)*nrhfa(isymj),
     &          nmatkl(isymkl),mbas1(isymb),factor,work(koffg),ntotb,
     &          work(koffr),ntotb,one,vajkl(koffv),ntotaj)
          
           if (locdbg) then
             write(lupri,*)'d: G'
             call output(work(koffg),1,mbas1(isymb),
     &                1,mbas1(isyma)*nrhfa(isymj),
     &                mbas1(isymb),mbas1(isyma)*nrhfa(isymj),1,lupri)
             write(lupri,*)'d: R'
             call output(work(koffr),1,mbas1(isymb),
     &                1,nmatkl(isymkl),
     &                mbas1(isymb),nmatkl(isymkl),1,lupri)
             write(lupri,*)'d: V = V + G*R'
             call output(vajkl(koffv),1,mbas1(isyma)*nrhfa(isymj),1,
     &          nmatkl(isymkl),nalphaj(isymaj),
     &          nmatkl(isymkl),1,lupri)
           end if
        end do
      else if (iansatz.eq.2) then
c       ----------------------------------------------
c       contract G^dp_m,alpha_j*R^dp_m,kl to obtain V
c       ----------------------------------------------
        do isymm = 1, nsym
c         m runs over all occupied orbitals
          isymmj = muld2h(isymm,isymj)
          isymdm = muld2h(isymd,isymm)
          isymkl = isymdm
          isymaj = isymdm
          isyma = muld2h(isymaj,isymj)

          koffg = kgtf + igmaj(isymaj)
          koffr = krtf + imaklm(isymkl,isymm)
          koffv = 1 + ivajkl(isymaj,isymkl) + ialphaj(isyma,isymj)

          ntotg = max(nrhfsa(isymm),1)
          ntotaj = max(nalphaj(isymaj),1)

          call dgemm('T','N',mbas1(isyma)*nrhfa(isymj),nmatkl(isymkl),
     &        nrhfsa(isymm),factor,work(koffg),ntotg,work(koffr),ntotg,
     &        one,vajkl(koffv),ntotaj)
        end do 
      else
        call quit('illegal case.')
      end if
      if (locdbg.and.iansatz.eq.2) then
        write(lupri,*)'Vajkl ansatz 2 auxbas first part:'
        write(lupri,*)'norm^2:',ddot(nvajkl(isymv),vajkl,1,vajkl,1)
        do isymkl = 1, nsym
           isymaj = muld2h(isymv,isymkl)
           write(lupri,*) 'Block isymaj, isymkl: ',isymaj,isymkl
           if ((nalphaj(isymaj).eq.0).or.(nmatkl(isymkl).eq.0)) then
            write(lupri,*) 'This block is empty'
           else
            call output(vajkl(1+ivajkl(isymaj,isymkl)),1,
     &                  nalphaj(isymaj),1,nmatkl(isymkl),
     &                  nalphaj(isymaj),nmatkl(isymkl),1,lupri)
           end if
        end do
      end if

c     ---------------------------------------------
c      case: delta element of auxiliary basis
c     ---------------------------------------------
      if ((idelta.gt.mbas1(isymd)) .or. (r12cbs.and.(iansatz.eq.2))
     &      .and. mbas1(isymg).gt.0) then
         ngbgi = 0
         do isymb = 1, nsym  
            isyma = muld2h(isymab,isymb)
            isymi = muld2h(isymcmo,isyma)
            isymgi = muld2h(isymg,isymi)
            igbgi(isymgi) = ngbgi
            ngbgi  = ngbgi + mbas1(isymb)*mbas1(isymg)*nrhfa(isymi)
         end do
         kgab = kend0
         kgbi = kgab + nbast*nbast
         kgbgi = kgbi + nbast*nrhft
         krblk = kgbgi + ngbgi
         kend1 = krblk + nrgkl(isymd)

         if (iansatz.eq.2) then
           kgtf    = kend1
           krdpmlk = kgtf + ngmgi 
           kend1   = krdpmlk + nmaklm(isymd)
         end if

         lwrk1 = lwork - kend1
         if (lwrk1.lt.0) then
           call quit('Insufficient work space in r12mkvamkl')
         end if

c        --------------------------------------------------------------
c         get g_{alpha beta gamma delta'} and unpack as quadratic
c         matrix (scr1) then contract g with cmo to obtain 
c         g_{i beta gamma delta'} and pack as G^{delta'}(beta,gamma i)
c         
c         G^{delta'}(beta,gamma i) =
c                 sum_{alpha} C_{alpha i} * g_{alpha beta gamma delta'}
c        --------------------------------------------------------------
         do g =1, mbas1(isymg)
            koff1 = nnbst(isymab)*(g-1) +1
c Elena
            if (u21int) then
               call ccsd_symsq(xint(koff1),isymab,scr1)
               call ccsd_asymsq(tint(koff1),isymab,scr1,1,0)
c Elena
            else 
               if (fnvajkl .eq. 'CCR12BIJAL') then
                  call ccsd_symsqt(xint(koff1),isymab,scr1)
                else
                  call  ccsd_symsq(xint(koff1),isymab,scr1)
                endif
            end if 
            do isyma = 1, nsym
               isymb  = muld2h(isymab,isyma)
               isymi  = muld2h(isymcmo,isyma)
               isymag = muld2h(isyma,isymg)
               isymgi = muld2h(isymg,isymi)

               do a = 1, mbas1(isyma)
                  do b = 1, mbas1(isymb)
                     koff1 = iaodis(isyma,isymb)+nbas(isyma)*(b-1)+a
                     idxab = mbas1(isymb)*(a-1)+b
                     work(kgab-1+idxab) = scr1(koff1)
c                     write(*,*) 'vamkl', koff1,idxab
                  end do
               end do

               if (fnvajkl .eq. 'CCR12QIJAL') then
                  koffc = 1 + iglmro(isyma,isymi)
               else
                  koffc = 1 + iglmrh(isyma,isymi)
               endif
               ntotb = max(mbas1(isymb),1)
               ntota = max(nbas(isyma),1)
                 
               call dgemm('N','N',mbas1(isymb),
     &              nrhfa(isymi),mbas1(isyma),one,work(kgab),ntotb,
     &              cmo(koffc),ntota,zero,work(kgbi),ntotb)

                 
               do i =1, nrhfa(isymi)
                  idxgi = mbas1(isymg)*(i-1)+g
                  do b = 1, mbas1(isymb)
                     idxbi = mbas1(isymb)*(i-1)+b
                     idxbgi = igbgi(isymgi)+mbas1(isymb)*(idxgi-1)+b
                     work(kgbgi-1+idxbgi) = work(kgbi-1+idxbi)
                  end do
               end do 
            end do
         end do
    
         if (iansatz.eq.2) then
c          G^{delta'}(m,gamma i) = sum_beta G^{delta'}(beta,gamma i) 
c                                     * C_{beta m}
c          m runs over all occupied orbitals
           do isymgi = 1, nsym
             isymi = muld2h(isymgi,isymg)
             isymm = muld2h(isymgi,isymd)
             isymb = isymm
           
             koffg = kgbgi + igbgi(isymgi) 
             koffc = iglmrhs(isymb,isymm) + 1
             koff1 = kgtf + igmgi(isymgi)
           
             ntotb = max(mbas1(isymb),1)
             ntotc = max(nbas(isymb),1)
             ntotm = max(nrhfsa(isymm),1)
           
             call dgemm('T','N',nrhfsa(isymm),
     &                mbas1(isymg)*nrhfa(isymi),
     &                mbas1(isymb),one,xlamdps(koffc),ntotc,work(koffg),
     &                ntotb,zero,work(koff1),ntotm)
           end do
         end if
c  
         if (locdbg.and.iansatz.eq.1) then
           write(lupri,*) 'resorted coulomb integrals (second part):',
     &       ddot(ngbgi,work(kgbgi),1,work(kgbgi),1)
           do isymi = 1, nsym
             isymb  = muld2h(isymab,isymi)
             isymgi = muld2h(isymg,isymi)
             write(lupri,*) 'isymi,isymg,isymb:',isymi,isymg,isymb
             call output(work(kgbgi+igbgi(isymgi)),1,mbas1(isymb),
     &                   1,mbas1(isymg)*nrhfa(isymi),
     &                   mbas1(isymb),mbas1(isymg)*nrhfa(isymi),
     &                   1,lupri)
           end do
         else if (locdbg.and.iansatz.eq.2) then
            write(lupri,*)'resorted coulomb integrals (second part):',
     &       ddot(ngbgi,work(kgbgi),1,work(kgbgi),1)
            write(lupri,*)'Norm: G^dp_m,gammai in mkvamkl',
     &         ddot(ngmgi,work(kgtf),1,work(kgtf),1)
            write(lupri,*)'G^dp_m,gammai in mkvamkl',
     &                (work(kgtf-1+i),i=1,ngmgi)
            write(lupri,*)'Norm Vajkl ansatz 2:',
     &          ddot(nvajkl(isymv),vajkl,1,vajkl,1)
         end if
c
         if (iansatz.eq.1) then
c          -------------------------------------------------
c          pack R^{delta'}(beta,kl) as R^{delta'}(beta,lk)           
c          -------------------------------------------------
           do isyml = 1, nsym
              isymdl = muld2h(isymd,isyml)
              isymbk = isymdl
              do l=1,nrhfb(isyml)
                 do isymb = 1, nsym
                    isymk  = muld2h(isymbk,isymb)
                    isymkl = muld2h(isymk,isyml)
                    do k=1,nrhfb(isymk)
                       do b=1, mbas1(isymb)
                          idxkl=imatkl(isymk,isyml)+
     &                            nrhfb(isymk)*(l-1)+k
                          idxlk=imatkl(isyml,isymk)+
     &                            nrhfb(isyml)*(k-1)+l
                          idxbkl = irgkl(isymb,isymkl) + 
     &                         mbas1(isymb)*(idxkl-1)+b
                          idxblk = irgkl(isymb,isymkl) + 
     &                         mbas1(isymb)*(idxlk-1)+b
                          work(krblk-1+idxblk) = work(kr2bkl-1+idxbkl)
                       end do
                    end do
                 end do
              end do
           end do
           if (locdbg) then
             write(lupri,*) 'resorted R12 integrals:',
     &         ddot(nrgkl(isymd),work(krblk),1,work(krblk),1)
             do isymb = 1, nsym
              isymkl = muld2h(isymd,isymb)
              write(lupri,*) 'isymb,isymd,isymkl:',isymb,isymd,isymkl
              call output(work(krblk+irgkl(isymb,isymkl)),1,
     &                    mbas1(isymb),1,nmatkl(isymkl),
     &                    mbas1(isymb),nmatkl(isymkl),1,lupri)
             end do
           end if
         else if (iansatz.eq.2) then
c          ---------------------------
c          pack R^dp_m,kl as R^dp_m,lk 
c          ---------------------------
c          m runs over all occupied
           do isyml = 1, nsym
             isymdl = muld2h(isymd,isyml)
             isymmk = isymdl
             do l = 1, nrhfb(isyml)
               do isymm = 1, nsym
                 isymk = muld2h(isymmk,isymm)
                 isymkl = muld2h(isymk,isyml)
                 do k = 1, nrhfb(isymk)
                   do m = 1, nrhfsa(isymm)
                     idxkl=imatkl(isymk,isyml)+nrhfb(isymk)*(l-1)+k
                     idxlk=imatkl(isyml,isymk)+nrhfb(isyml)*(k-1)+l
                     idxmkl =imaklm(isymkl,isymm)+
     &                       nrhfsa(isymm)*(idxkl-1)+m
                     idxmlk =imaklm(isymkl,isymm)+
     &                       nrhfsa(isymm)*(idxlk-1)+m
                     work(krdpmlk+idxmlk-1) = work(krtf-1+idxmkl)       
                   end do 
                 end do
               end do 
             end do
           end do

           if (locdbg) then
             write(lupri,*) 'norm^2 of resorted R^dp_m,lk:',idel,
     &        ddot(nmaklm(isymd),work(krdpmlk),1,work(krdpmlk),1)
             write(lupri,*)'R^dp_m,lk', 
     &                      (work(krdpmlk-1+i),i=1,nmaklm(isymd))
           end if
         end if


         if (iansatz.eq.1) then
c          -------------------------------------------------------------
c           contract  G^{delta'}(beta,gamma i) with R^{delta'}(beta,lk) 
c           to obtain V(gamma i, kl)
c          -------------------------------------------------------------
           do isymb = 1, nsym
             isymbd = muld2h(isymb,isymd)
             isymkl = isymbd
             isymgi = muld2h(isymv,isymbd)
             isymi  = muld2h(isymgi,isymg)
c            isymgi = muld2h(isymd,isymb)
c            isymkl = muld2h(isymd,isymb)
c            isymi  = muld2h(isymgi,isymg)
          
             koffg = kgbgi + igbgi(isymgi)
             koffr = krblk + irgkl(isymb,isymkl)
             koffv = 1 + ivajkl(isymgi,isymkl) + ialphaj(isymg,isymi)
             
             ntotb  = max(1,mbas1(isymb))
             ntotgi = max(1,nalphaj(isymgi))
          
             call dgemm('T','N',mbas1(isymg)*nrhfa(isymi),
     &            nmatkl(isymkl),
     &            mbas1(isymb),-one,work(koffg),ntotb,
     &            work(koffr),ntotb,one,vajkl(koffv),ntotgi)
             if (locdbg) then 
               write(lupri,*)'DEBUG Vajkl(koffv) auxbas'
               call output(vajkl(koffv),1,mbas1(isymg)*nrhfa(isymi),1,
     &           nmatkl(isymkl),nalphaj(isymgi),nmatkl(isymkl),
     &           1,lupri)
             end if 
           end do

        else if (iansatz.eq.2) then 
c          -----------------------------------------------------
c          contract G^{delta'}(m,gamma i) with  R^{delta'}(m,lk) 
c          -----------------------------------------------------
c          m runs over all occupied
           do isymm = 1, nsym
             isymgi = muld2h(isymd,isymm)
             isymkl = muld2h(isymd,isymm)
             isymi  = muld2h(isymgi,isymg)

             koffg = kgtf + igmgi(isymgi)
             koffr = krdpmlk + imaklm(isymkl,isymm)
             koffv = 1 + ivajkl(isymgi,isymkl) + ialphaj(isymg,isymi)

             ntotm  = max(nrhfsa(isymm),1)
             ntotgi = max(nalphaj(isymgi),1)
             
             call dgemm('T','N',mbas1(isymg)*nrhfa(isymi),
     &           nmatkl(isymkl),
     &           nrhfsa(isymm),-one,work(koffg),ntotm,work(koffr),ntotm,
     &           one,vajkl(koffv),ntotgi)

c
             xnrmr = xnrmr + ddot(nrhfsa(isymm)*nmatkl(isymkl),
     &                           work(koffr),1,work(koffr),1)
             xnrmg = xnrmg + ddot(nrhfsa(isymm)*mbas1(isymg)*
     &                           nrhfa(isymi),
     &                           work(koffg),1,work(koffg),1)
c
             if (locdbg) then 
               write(lupri,*)'Vajkl(koffv) ansatz 2 auxbas'
               call output(vajkl(koffv),1,mbas1(isymg)*nrhfa(isymi),1,
     &                    nmatkl(isymkl),nalphaj(isymgi),
     &                    nmatkl(isymkl),1,lupri)
             end if 
           end do 
        end if
      end if

      if (locdbg) then
        write(lupri,*)'Vajkl leaving mkvamkl'
        write(lupri,*)'norm^2:',ddot(nvajkl(isymv),vajkl,1,vajkl,1)
c       write(lupri,*)(vajkl(i),i=1,nvajkl(isymv))
        do isymkl = 1, nsym
           isymaj = muld2h(isymv,isymkl)
           write(lupri,*) 'Block isymaj, isymkl: ',isymaj,isymkl
           if ((nalphaj(isymaj).eq.0).or.(nmatkl(isymkl).eq.0)) then
            write(lupri,*) 'This block is empty'
           else
            call output(vajkl(1+ivajkl(isymaj,isymkl)),1,
     &                  nalphaj(isymaj),1,nmatkl(isymkl),
     &                  nalphaj(isymaj),nmatkl(isymkl),1,lupri)
           end if
         end do
        write(lupri,*)'xnrmr, xnrmg:',xnrmr,xnrmg
        write(lupri,*)'x1, x2, x3:',x1, x2, x3
      end if

      call qexit('r12mkvamkl')
      end                       
*====================================================================*
      subroutine cc_r12mkvamkl0(vxykl,length,iopt,cmo,isycmo,
     &                          work,lwork) 
c--------------------------------------------------------------------
c     purpose: initialize R12 V-intermediates 
c  
c     H. Fliegl, C. Haettig spring 2003
c     modified for CCSD(R12) model, summer 2004
c
c     C. Neiss august/september 2005:
c     generalized and restructured for arbitrary r12-orbitals,
c     non-totalsymm. Vajkl, Vijkl
c     Note that the meaning of variable "cmo" has changed!
c     iopt = 0      initialize V_(alpha beta)^(kl)
c     iopt = 1      initialize V_(alpha j)^(kl)
c     iopt = 2      initialize V_(ij)^(kl) 
c--------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
#include <ccr12int.h>
#include <dummy.h>
 
      logical locdbg
      parameter (locdbg = .false.)
#if defined (SYS_CRAY)
      real zero,one
#else
      double precision zero,one
#endif
      parameter (zero = 0.0d0, one = 1.0d0)

      integer nr1orb(8),nr1bas(8),nr1xbas(8),n2bst1(8),nrgkl(8)
      integer nr2bas(8)
      integer ir1orb(8,8),ir1bas(8,8),ir2bas(8,8),ir1xbas(8,8)
      integer ir2xbas(8,8),iaodis1(8,8),irgkl(8,8)
      integer nr1xorb(8),ir1xorb(8,8),nrxgkl(8),irxgkl(8,8)
      integer igabj(8),nalphaj(8),ialphaj(8,8)
      integer isymi,isymj,isymk,isymkl,isycmo
      integer isymab,isyma,isymb,isymaj,isyml
      integer idxab,idxak,idxkl,idxlk,idxlj,idxajkl,koffc
      integer ksab,ksak,kski,kend1,lwrk1,lwork,koffsak,koffski
      integer idxaj,nnbast
      integer idxabkl,idxbl,lusak,luvajkl,idxbj
      integer iopt,length,isymij,idxij,idxki,idxijkl
      integer lr12f12htf
      
#if defined (SYS_CRAY)
      real vxykl(length),cmo(*),work(*), ddot
#else
      double precision vxykl(length),cmo(*),work(*), ddot
#endif
      
      call qenter('r12mkvamkl0')

      call cc_r12offset(nr1orb,nr1xorb,nr1bas,nr1xbas,nr2bas,
     &     nrgkl,nrxgkl,n2bst1,ir1orb,ir1xorb,ir1bas,ir1xbas,ir2bas,
     &     ir2xbas,irgkl,irxgkl,iaodis1,nalphaj,ialphaj)


      nnbast = 0
      do i = 1,nsym
         nnbast = nnbast + nbas(i)*(nbas(i)+1)/2
      end do
c
      kend1 = 1
      if (.not.R12EOR) then
        ksak = kend1
        kend1 = ksak + nr1bas(1)  
       if (iopt.ne.0) then
         kski = kend1
         kend1 = kski + nmatki(1)
       end if
      end if

      lwrk1 = lwork - kend1
      if (lwrk1.lt.0) then
         write(lupri,*)'lwrk1,lwork',lwrk1,lwork
         call quit('Insufficient work space in cc_r12mkvamkl0.')
      end if

c    ------------------------------
c     read overlap matrix from file
c    ------------------------------
      if (.not. R12EOR) then
      lusak = -1
      call gpopen(lusak,fnsak,'UNKNOWN',' ','UNFORMATTED',
     &            idummy,.false.)
      rewind(lusak)
      read(lusak) (work(ksak+i-1), i=1, nr1bas(1))
      call gpclose(lusak,'KEEP')

c     -------------------------------------------
c     transform S_{alpha k} with C_{alpha i} 
c     to obtain S_{ki}
c     -------------------------------------------
      if (iopt.ne.0) then   
        do isymk = 1, nsym
          isyma = isymk
          isymi = muld2h(isyma,isycmo)
          koffc = 1 + iglmrh(isyma,isymi)
          koffsak = ksak + ir1bas(isyma,isymk)
          koffski = kski + imatki(isymk,isymi)
          call dgemm('T','N',nrhfb(isymk),nrhfa(isymi),
     &               mbas1(isyma),one,work(koffsak),
     &               max(1,mbas1(isyma)),cmo(koffc),
     &               max(1,nbas(isyma)),zero,work(koffski),
     &               max(1,nrhfb(isymk)))  
        end do
      end if
c
      end if
c  
      call dzero(vxykl,length)
c
      if (R12EOR) then
        call cc_r12getvxykl(vxykl,length,cmo,isycmo,nr1bas,
     &                      ir1bas,ir2bas,nalphaj,ialphaj,
     &                      lr12f12htf,fr12f12htf,
     &                      iopt,work(kend1),lwrk1)
      else
c
      if (iopt.eq.0) then
c       ----------------------------------------------
c       initialize V(ak,bl) with S(alpha,k)*S(beta,l)
c       ----------------------------------------------
        do isymkl = 1, nsym
          isymab = isymkl
          do isymk = 1, nsym
            isyml = muld2h(isymkl,isymk) 
            isyma = isymk
            isymb = isyml
            do k = 1, nrhfb(isymk)
              do l = 1, nrhfb(isyml)
                idxkl = imatkl(isymk,isyml)+nrhfb(isymk)*(l-1)+k
                do a = 1, mbas1(isyma)
                  idxak = ir1bas(isyma,isymk)+mbas1(isyma)*(k-1)+a
                  do b = 1, mbas1(isymb)
                    idxab = iaodis(isyma,isymb)+mbas1(isyma)*(b-1)+a
                    idxbl = ir1bas(isymb,isyml)+mbas1(isymb)*(l-1)+b
                    idxabkl = ivabkl(isymab,isymkl)+
     &                        n2bst(isymab)*(idxkl-1)+idxab
                    vxykl(idxabkl) = 
     &                        work(ksak-1+idxak)*work(ksak-1+idxbl)
                  end do
                end do
              end do
            end do
          end do
        end do 
c
      else if (iopt.eq.1) then
c       -------------------------------------------------
c       initialize V(aj,kl) with S(alpha,k)*S(j,l)
c       -------------------------------------------------
        do isymkl = 1, nsym
          isymaj = muld2h(isymkl,isycmo)
          do isymk = 1, nsym
            isyml = muld2h(isymkl,isymk)
            isyma = isymk
            isymj = muld2h(isyml,isycmo)

            do k = 1, nrhfb(isymk)
              do l = 1, nrhfb(isyml)
                idxkl = imatkl(isymk,isyml)+nrhfb(isymk)*(l-1)+k
                do a = 1, mbas1(isyma)
                  idxak = ir1bas(isyma,isymk)+mbas1(isyma)*(k-1)+a
                  do j = 1, nrhfa(isymj)
                    idxaj = ialphaj(isyma,isymj)+mbas1(isyma)*(j-1)+a
                    idxlj = imatki(isyml,isymj)+nrhfb(isyml)*(j-1)+l
                    idxajkl = ivajkl(isymaj,isymkl)+
     &                        nalphaj(isymaj)*(idxkl-1)+idxaj
                    vxykl(idxajkl) =
     &                        work(ksak-1+idxak)*work(kski-1+idxlj)

                  end do
                end do
              end do
            end do
          end do
        end do
c
      else if (iopt.eq.2) then
c       -------------------------------------------------
c       initialize V(ij,kl) with S(i,k)*S(j,l)
c       -------------------------------------------------
        do isymkl = 1, nsym
          isymij = muld2h(isymkl,isycmo)
          do isymk = 1, nsym
            isyml = muld2h(isymkl,isymk)
            isymi = muld2h(isymk,isycmo) 
            isymj = muld2h(isyml,isycmo)
            do k = 1, nrhfb(isymk)
              do l = 1, nrhfb(isyml)
                idxkl = imatkl(isymk,isyml)+nrhfb(isymk)*(l-1)+k
                do i = 1, nrhfa(isymi)
                  idxki = imatki(isymk,isymi)+nrhfb(isymk)*(i-1)+k
                  do j = 1, nrhfa(isymj)
                    idxij = imatij(isymi,isymj)+nrhfa(isymi)*(j-1)+i
                    idxlj = imatki(isyml,isymj)+nrhfb(isyml)*(j-1)+l
                    idxijkl = itr12sqt(isymij,isymkl)+
     &                        nmatij(isymij)*(idxkl-1)+idxij
                    vxykl(idxijkl) =
     &                        work(kski-1+idxki)*work(kski-1+idxlj)
                  end do
                end do
              end do
            end do
          end do
        end do
c
      else 
        call quit('Unknown value for IOPT in CC_R12MKVAMKL0')
      end if
c
      end if !(R12EOR)
c
      if (locdbg) then
        write(lupri,*) 'in cc_r12mkvamkl0: iopt, norm^2(vxykl):', iopt,
     &  ddot(length,vxykl,1,vxykl,1)
        if (iopt.eq.0) then
         do isymab = 1, nsym
          isymkl = isymab
          write(lupri,*) 'symmetry block ',isymab,isymkl
          call output(vxykl(1+ivabkl(isymab,isymkl)),1,n2bst(isymab),
     &         1,nmatkl(isymkl),n2bst(isymab),nmatkl(isymkl),1,lupri)
         end do
c
        else if (iopt.eq.1) then
         do isymkl = 1, nsym
          isymaj  = muld2h(isycmo,isymkl)
          write(lupri,*) 'symmetry block ',isymaj,isymkl
          call output(vxykl(1+ivajkl(isymaj,isymkl)),1,nalphaj(isymaj),
     &         1,nmatkl(isymkl),nalphaj(isymaj),nmatkl(isymkl),1,lupri)
         end do
c
        else if (iopt.eq.2) then
         do isymij = 1, nsym
          isymkl = isymij
          write(lupri,*)'symmetry block (ij,kl)',isymij,isymkl
          call output(vxykl(1+itr12sqt(isymij,isymkl)),1,
     &         nmatij(isymij),1,nmatkl(isymkl),nmatij(isymij),
     &         nmatkl(isymkl),1,lupri)
         end do
        end if
      end if
c     
      call qexit('r12mkvamkl0')       
      end                       
*=====================================================================*
      subroutine cc_r12mkvijkl(vajkl,isymv,cmo,isycmo,work,lwork,lres,
     &                         fac,vijkl)
c---------------------------------------------------------------------
c     purpose: calculate V(ij,kl) to test if V(alpha j,kl)
c              if correct; this V(ij,kl) must be exactly the 
c              same as the MP2-R12 V(ij,kl)
c
c     H. Fliegl, C. Haettig spring 2003
c     C. Neiss 29-04-2005: result is always added on vijkl (if lres=.T.)
c                          with a factor FAC
c---------------------------------------------------------------------
      implicit none 
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
#if defined (SYS_CRAY)
      real zero,one
#else
      double precision zero,one
#endif
      parameter (zero = 0.0d0, one = 1.0d0)

      logical lres, locdbg
      parameter (locdbg = .false.)

      integer nr1orb(8),nr1bas(8),nr1xbas(8),n2bst1(8),nrgkl(8)
      integer nr2bas(8),nalphaj(8),ialphaj(8,8)
      integer ir1orb(8,8),ir1bas(8,8),ir2bas(8,8),ir1xbas(8,8)
      integer ir2xbas(8,8),iaodis1(8,8),irgkl(8,8)
      integer nr1xorb(8),ir1xorb(8,8),nrxgkl(8),irxgkl(8,8)
      integer igabj(8),isycmo,isymv
      integer isyma,isymk,isymi,isymj,isyml,isymaj,isymkl
      integer idxkl,nrhftria,isymij,idxijkl
      integer lwork,ntota,ntoti,ntota1
      integer kvajkl,koffc,kscr1,kvijkl,ksing,ktrip,kend1 
      integer ik,il,ii,ij,idxij,iklij,lwrk1
#if defined (SYS_CRAY)
      real vajkl(*),cmo(*),work(*),vijkl(*),ddot,fac
#else
      double precision vajkl(*),cmo(*),work(*),vijkl(*),ddot,fac
#endif

      call qenter('r12mkvijkl')

      call cc_r12offset(nr1orb,nr1xorb,nr1bas,nr1xbas,nr2bas,
     &     nrgkl,nrxgkl,n2bst1,ir1orb,ir1xorb,ir1bas,ir1xbas,ir2bas,
     &     ir2xbas,irgkl,irxgkl,iaodis1,nalphaj,ialphaj)

      if (locdbg) then
        write(lupri,*) 'in cc_r12mkvijkl: vajkl',
     &     ddot(nvajkl(isymv),vajkl,1,vajkl,1)
        do isymkl = 1, nsym 
          isymaj  = muld2h(isymv,isymkl)
          write(lupri,*) 'symmetry block ',isymaj,isymkl 
          call output(vajkl(1+ivajkl(isymaj,isymkl)),1,nalphaj(isymaj),
     &        1,nmatkl(isymkl),nalphaj(isymaj),nmatkl(isymkl),1,lupri)
        end do
      end if

      nrhftria = nrhftb*(nrhftb+1)/2

      kscr1  = 1
      kend1 = kscr1  + nrhftb*nrhftb
      if (.not. lres) then
        kvijkl = kend1 
        ksing  = kvijkl + nrhftb*nrhftb*nrhftb*nrhftb
        ktrip  = ksing  + nrhftria * nrhftria
        kend1  = ktrip  + nrhftria * nrhftria
      end if

      lwrk1  = lwork -kend1
      if (lwrk1 .lt.0) then
         call quit('Insufficient work space in r12mkvijkl')
      end if

      if (.not. lres) then
        call dzero(work(kvijkl),nrhftb*nrhftb*nrhftb*nrhftb)
      end if 

      do isymkl = 1, nsym
         isymaj = muld2h(isymv,isymkl)
         do isyml =1, nsym
            isymk = muld2h(isymkl,isyml)
            do k = 1, nrhfb(isymk)
               ik = irhfb(isymk) + k
               do l = 1, nrhfb(isyml)
                  il = irhfb(isyml) + l
                  idxkl=imatkl(isymk,isyml)+nrhfb(isymk)*(l-1)+k
                  do isyma =1, nsym
                     isymj = muld2h(isymaj,isyma)
                     isymi = muld2h(isycmo,isyma)

                     koffc = iglmrh(isyma,isymi)+1
                     kvajkl = ivajkl(isymaj,isymkl)+nalphaj(isymaj)*
     &                       (idxkl -1) +ialphaj(isyma,isymj)+1
                     ntota1 = max(1,mbas1(isyma))
                     ntota = max(1,nbas(isyma))
                     ntoti = max(1,nrhfa(isymi))

                     call dgemm('T','N',nrhfa(isymi),
     &                    nrhfa(isymj),mbas1(isyma),one,cmo(koffc),
     &                    ntota,vajkl(kvajkl),ntota1,zero,work(kscr1),
     &                    ntoti)
c                    ------------------------------------------
c                     case lres (response) pack V(it,jt,kl) as
c                     in vector function
c                    ------------------------------------------
                     if (lres) then
                       isymij  = muld2h(isymi,isymj)
                       idxijkl = itr12sqt(isymij,isymkl)+ 
     &                 nmatij(isymij)*(idxkl-1)+imatij(isymi,isymj)+1
                       call daxpy(nrhfa(isymi)*nrhfa(isymj),fac,
     &                            work(kscr1),1,vijkl(idxijkl),1)
                     else
                       do j = 1, nrhfa(isymj)
                         ij = irhfa(isymj) + j
                         do i = 1, nrhfa(isymi)
                           ii = irhfa(isymi) + i
                           idxij = nrhfa(isymi)*(j-1) + i
                           iklij=nrhftb*(nrhftb*(nrhftb*(il-1)+
     &                                  ik-1)+ij-1)+ii
                           work(kvijkl-1+iklij) = work(kscr1-1+idxij)
                         end do
                       end do
                     end if
c
                  end do
               end do
            end do
         end do
      end do

      call qexit('r12mkvijkl') 
      end 
*=====================================================================*
      subroutine cc_r12mkvijkl2(vijkl,vsing,vtrip,nr12t,nrhftria)
c---------------------------------------------------------------------
c     purpose: make V for s = 0,1 to compare with MP2-R12
c
c     H. Fliegl, C. Haettig spring 2003
c--------------------------------------------------------------------
      implicit none 
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>

      integer ij,kl,nrhftria,nr12t
#if defined (SYS_CRAY)
      real vijkl(nr12t,nr12t,nr12t,nr12t),ff
      real vsing(nrhftria,nrhftria),vtrip(nrhftria,nrhftria)
#else
      double precision vijkl(nr12t,nr12t,nr12t,nr12t),ff
      double precision vsing(nrhftria,nrhftria)
      double precision vtrip(nrhftria,nrhftria)
#endif

      ff = 1.0d0/sqrt(2.0d0)

      ij = 0
      do i=1,nr12t
         do j=1,i
            ij = ij + 1
            kl = 0
            do k=1,nr12t
               do l=1, k
                  kl = kl + 1
                  vsing(kl,ij) = vijkl(i,j,k,l)+ vijkl(i,j,l,k)
                  if (i.eq.j) vsing(kl,ij) = ff * vsing(kl,ij)
                  if (k.eq.l) vsing(kl,ij) = ff * vsing(kl,ij)
                  vtrip(kl,ij) = vijkl(i,j,k,l)- vijkl(i,j,l,k)
                  if (i.eq.j) vtrip(kl,ij) = ff * vtrip(kl,ij)
                  if (k.eq.l) vtrip(kl,ij) = ff * vtrip(kl,ij)
               end do
            end do
         end do
      end do
      write(lupri,*)'Vsing'
      call output(vsing,1,nrhftria,1,nrhftria,nrhftria,nrhftria,1,lupri)
      write(lupri,*)'Vtrip'
      call output(vtrip,1,nrhftria,1,nrhftria,nrhftria,nrhftria,1,lupri)

      end 
*=====================================================================*
      subroutine cc_r12getrint(r12bkl,idel,isymd,nr1bas,ir1bas,
     &     nr2bas,ir2bas,nrgkl,irgkl,ir1xbas,ir2xbas,
     &     nrhf1,nmatij1,imatij1,
     &     ibasx,lauxd,lauxbeta,fconst,work,lwork)
c---------------------------------------------------------------------
c     purpose: get R12 integrals 
c
c     H. Fliegl, C. Haettig spring 2003
c     
c     added possilbility to read r12_(beta,kl) with beta being 
c     auxiliary function:  C. Neiss  winter 2004/2005
c
c     nr1bas     dimension (mbas1 OR mbas2) x nrhfb
c     ir1bas     symmetry offset of mbas1 x nrhfb
c     ir1xbas    symmetry offset of mbas2 x nrhfb
c     nr2bas     length of nr1bas OR nr1xbas over all symmetries
c     ir2bas     symmetry offset nr1bas x nr1bas OR nr1xbas x nr1bas
c     nrgkl      dimension of three index array (result)
c     irgkl      symmetry offset of three index array
c     ir2xbas    symmetry offset nr1bas x nr1xbas OR nr1xbas x nr1xbas
c     lauxbeta   beta is auxiliary function
c     fconst     file which contains integrals 
c
c     Note:   you have to specify different files and dimensions 
c             depending on "lauxbeta"! (see also subroutine r12batf
c             for dimensions and cc_r12xvint for an example of use)
c
c---------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>

      logical lauxd,locdbg,lauxbeta 
      parameter(locdbg=.false.)

      character*(*) fconst

#if defined (SYS_CRAY)
      real zero,one
#else
      double precision zero,one
#endif
      parameter (zero = 0.0d0, one = 1.0d0)

      integer nr1bas(8),ir1bas(8,8),nr2bas,ir2bas(8,8)
      integer ir2xbas(8,8),nrgkl(8)
      integer irgkl(8,8),ir1xbas(8,8)
      integer nrhf1(8), nmatij1(8), imatij1(8,8)
      integer idelta,idel,lwork,lwrk1,ibasx(8)
      integer lur2back,kscr2,kend2,lwrk2,iadr,ndell,ilen
      integer isymb,isymk,isymkl,isyml,idxbk,idxkl,idxbkl
      integer isymd,isymdl,isymbk,kend1 


#if defined (SYS_CRAY)
      real r12bkl(*),work(*),ddot
#else
      double precision r12bkl(*),work(*),ddot
#endif
      
      call qenter('cc_r12getrint')

      kend1 = 1

      if (locdbg) then
        do i = 1, nsym
          write(lupri,*) 'nr1bas(',i,') = ',nr1bas(i)
        end do
      end if

      idelta = idel - ibas(isymd)
      if (lauxd) idelta = idelta - ibasx(isymd)

      lur2back = -1
      call wopen2(lur2back,fconst,64,0)
      if (locdbg) then
        write(lupri,*) 'integrals read from file ',fconst
      end if
      
      do isyml = 1, nsym
         isymdl = muld2h(isymd,isyml)
         isymbk = isymdl
         
         kscr2 = kend1 
         kend2 = kscr2 + nr1bas(isymbk)
         lwrk2 = lwork - kend2
         if (lwrk2 .lt. 0) then
           write(lupri,*)'lwork',lwork
           write(lupri,*)'lwrk2,kend2,kscr2', lwrk2,kend2,kscr2
           call quit('Insufficient work space in cc_r12getrint')
         end if
c       ------------------------------------------------------------
c        read r^{beta,delta}_{k,l} from file and pack as 
c        R^{delta}(beta,kl) or  R^{delta'}(beta,kl)
c       ------------------------------------------------------------
         do l = 1, nrhf1(isyml)
           if (idelta.le.mbas1(isymd)) then
             ! delta in orbital basis
             ndell = ir1bas(isymd,isyml) +
     &               mbas1(isymd)*(l-1)+idelta
             iadr  = ir2bas(isymbk,isymdl)+
     &               nr1bas(isymbk)*(ndell-1)+1
           else 
             ! delta in auxiliary basis
             ndell = ir1xbas(isymd,isyml)+
     &               mbas2(isymd)*(l-1)+idelta-mbas1(isymd)
             iadr  = nr2bas + ir2xbas(isymbk,isymdl)+
     &               nr1bas(isymbk)*(ndell-1)+1
           end if
           if(locdbg) then
             write(lupri,*)'mbas1,nr1bas:',mbas1(isymd),nr1bas(isymbk)
             write(lupri,*)'mbas2,nr2bas:',mbas2(isymd),nr2bas
             write(lupri,*)'ndell,iadr,l,idelta',ndell,iadr,
     &              l,idelta
           end if
           ilen  = nr1bas(isymbk)
           call getwa2(lur2back,fconst,work(kscr2),iadr,ilen)
c 
           do isymb =1, nsym
             isymk = muld2h(isymbk,isymb)
             isymkl = muld2h(isymk,isyml)
             do k=1, nrhf1(isymk)
               if (.not. lauxbeta) then 
                 do b=1, mbas1(isymb)
                   idxbk = ir1bas(isymb,isymk) + mbas1(isymb)*(k-1)+b
                   idxkl = imatij1(isymk,isyml) + nrhf1(isymk)*(l-1)+k
                   idxbkl = irgkl(isymb,isymkl) + 
     &                      mbas1(isymb)*(idxkl-1)+b
                   if(locdbg) then
                     write(lupri,*)'idxbk,idxkl,idxbkl',
     &                       idxbk,idxkl,idxbkl
                   end if
                   r12bkl(idxbkl) = work(kscr2 -1 + idxbk)
                 end do
               else 
                 do b=1, mbas2(isymb)
                   idxbk = ir1xbas(isymb,isymk) + mbas2(isymb)*(k-1)+b
                   idxkl = imatij1(isymk,isyml) + nrhf1(isymk)*(l-1)+k
                   idxbkl = irgkl(isymb,isymkl) + 
     &                      mbas2(isymb)*(idxkl-1)+b 
                   if(locdbg) then
                     write(lupri,*)'idxbk,idxkl,idxbkl',
     &                       idxbk,idxkl,idxbkl
                   end if
                   r12bkl(idxbkl) = work(kscr2 -1 + idxbk)
                 end do
               end if
             end do

             if (locdbg) then
               write(lupri,'(2a,4i5)')' R^{delta}(beta,kl) for ',
     &           'delta,l,isymb,isymkl',idelta,l,isymb,isymkl
               idxbkl = irgkl(isymb,isymkl) + 1
               if (.not. lauxbeta) then
                 call output(r12bkl(idxbkl),1,mbas1(isymb),1,
     &                       nmatij1(isymkl),mbas1(isymb),
     &                       nmatij1(isymkl),1,lupri)
               else 
                 call output(r12bkl(idxbkl),1,mbas2(isymb),1,
     &                       nmatij1(isymkl),mbas2(isymb),
     &                       nmatij1(isymkl),1,lupri) 
               end if 
             end if

           end do
        end do
      end do
      call wclose2(lur2back,fconst,'KEEP')

      if (locdbg) then
        write(lupri,*) 'leaving cc_r12getrint... norm^2(r12bkl):',
     &         ddot(nrgkl(isymd),r12bkl,1,r12bkl,1)
      end if

      call qexit('cc_r12getrint')

      end
*=====================================================================*
      subroutine cc_r12putrint(r12bkl,idel,isymd,nr1bas,ir1bas,
     &     nr2bas,ir2bas,nrgkl,irgkl,ir1xbas,ir2xbas,
     &     nrhf1,nmatij1,imatij1,
     &     ibasx,lauxd,fconst,work,lwork)
c---------------------------------------------------------------------
c     purpose: write R12 integrals on file  
c
c     H. Fliegl, C. Haettig winter 2003
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

      integer nr1bas(8),ir1bas(8,8),nr2bas,ir2bas(8,8)
      integer ir2xbas(8,8),nrgkl(8)
      integer irgkl(8,8),ir1xbas(8,8)
      integer nrhf1(8), nmatij1(8), imatij1(8,8)
      integer idelta,idel,lwork,lwrk1,ibasx(8)
      integer lur2back,kscr2,kend2,lwrk2,iadr,ndell,ilen
      integer isymb,isymk,isymkl,isyml,idxbk,idxkl,idxbkl
      integer isymd,isymdl,isymbk,kend1 


#if defined (SYS_CRAY)
      real r12bkl(*),work(*),ddot
#else
      double precision r12bkl(*),work(*),ddot
#endif
      call qenter('cc_r12putrint') 

      kend1 = 1

      idelta = idel - ibas(isymd)
      if (lauxd) idelta = idelta - ibasx(isymd)

      lur2back = -1
      call wopen2(lur2back,fconst,64,0)
      
      do isyml = 1, nsym
         isymdl = muld2h(isymd,isyml)
         isymbk = isymdl
         
         kscr2 = kend1 
         kend2 = kscr2 + nr1bas(isymbk)
         lwrk2 = lwork - kend2
         if (lwrk2 .lt.0) then
            write(lupri,*)'lwork',lwork
            write(lupri,*)'lwrk2,kend2,kscr2', lwrk2,kend2,kscr2
            call quit('Insufficient work space in cc_r12putrint')
         end if

c       ------------------------------------------------------------
c       pack R^{delta}(beta,kl) or  R^{delta'}(beta,kl) as 
c       r^{beta,delta}_{k,l} and write to file 
c       ------------------------------------------------------------
         do l = 1, nrhf1(isyml)

            do isymb =1, nsym
               isymk = muld2h(isymbk,isymb)
               isymkl = muld2h(isymk,isyml)
               do k=1, nrhf1(isymk)
                  do b=1, mbas1(isymb)
                     idxbk = ir1bas(isymb,isymk)+mbas1(isymb)*(k-1)+b
                     idxkl = imatij1(isymk,isyml)+nrhf1(isymk)*(l-1)+k
                     idxbkl = irgkl(isymb,isymkl) + 
     &                    mbas1(isymb)*(idxkl-1)+b
                     work(kscr2 -1 + idxbk) = r12bkl(idxbkl)
                  end do
               end do
            end do

            if (idelta.le.mbas1(isymd)) then
               ndell = ir1bas(isymd,isyml) +
     &              mbas1(isymd)*(l-1)+idelta
               iadr = ir2bas(isymbk,isymdl)+
     &              nr1bas(isymbk)*(ndell-1)+1
            else
               ndell = ir1xbas(isymd,isyml)+
     &              mbas2(isymd)*(l-1)+idelta-mbas1(isymd)
               iadr = nr2bas + ir2xbas(isymbk,isymdl)+
     &              nr1bas(isymbk)*(ndell-1)+1
            end if
            ilen   = nr1bas(isymbk)
            call putwa2(lur2back,fconst,work(kscr2),iadr,ilen)
         end do
      end do
      call wclose2(lur2back,fconst,'KEEP')

      call qexit('cc_r12putrint')    
      end
*=====================================================================*
      subroutine cc_r12mkvijkl0(vijkl)
c---------------------------------------------------------------------
c     purpose: initialize Vijkl 
c
c     H. Fliegl, C. Haettig spring 2003
c---------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
#if defined (SYS_CRAY)
      real zero,one
      real vijkl(*)
#else
      double precision zero,one
      double precision vijkl(*)
#endif
      parameter (zero = 0.0d0, one = 1.0d0)

      integer isymi,isymj,isymij,idxij,idxkl,idxijkl,icoun1
      integer nr12orb(8),nrhforb(8),nij(8),iij(8,8)

      call quit('CC_R12MKVIJKL0 is obsolete! Use CC_R12MKVAMKL0!')
      
C     do isymij = 1, nsym
C       icoun1 = 0
C       do isymj = 1, nsym
C         isymi = muld2h(isymj,isymij)
C         iij(isymi,isymj) = icoun1
C         icoun1 = icoun1 + nrhforb(isymi)*nrhforb(isymj)
C       end do
C       nij(isymij) = icoun1
C     end do

C     call dzero(vijkl,ntr12sq(1))
C     
C     do isymi = 1, nsym
C       do isymj = 1, nsym
C        isymij = muld2h(isymi,isymj)
C        do i = 1, nrhforb(isymi)
C           do j = 1, nrhforb(isymj)
C              idxij  = iij(isymi,isymj) + nrhforb(isymi)*(j-1)+i
C              idxkl  = imatkl(isymi,isymj) + nr12orb(isymi)*(j-1)+i
C              idxijkl= itr12sqt(isymij,isymij)+
C    &                  nij(isymij)*(idxkl-1)+idxij
C              vijkl(idxijkl) = one
C           end do
C        end do
C      end do
C     end do
      end
*=====================================================================*
      subroutine ccrhs_gp(omega1,lambp,work,lwork,iamp,isymc, 
     &                    fc12am,lufc12,ifile)
*---------------------------------------------------------------------*
c     purpose: calculate V^{at j}_{kl}
c
c     H. Fliegl, C. Haettig spring 2003 
c---------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <dummy.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h> 
#include <ccr12int.h> 
      logical ldum,etest,locdbg
      parameter(etest = .false., locdbg = .false.)
#if defined (SYS_CRAY)
      real zero,one,two
#else
      double precision zero,one,two
#endif
      parameter (zero = 0.0d0, one = 1.0d0, two = 2.0d0)
      character*(*) fc12am 
      character     cdummy*(8) 
      integer lwork,luvajkl,kvajkl,lwrk1,kend1
      integer idum,lu43
      integer kctil,isymmu,isymkl, isymij, isymk, isyml,
     &     isymaj, isymi, isymj, idxkl,
     &     isymmuj,koffl,koffv,kvmujkl,ntotmu,ntotaj,ntotmuj,
     &     isyma,ifile,isymc,isyomg,lufc12
      integer icount7,isymak,kvint,kend2,lwrk2,ntota,ntoti
      integer koffo,koffc,index
      integer iamp
#if defined (SYS_CRAY)
      real work(*),lambp(*),omega1(*),er12,ddot
#else
      double precision work(*),lambp(*),omega1(*),er12,ddot
#endif

      index(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j

      call qenter('ccrhs_gp')
      
      if (locdbg) then
        write(lupri,*)'in ccrhs_gp: omega1 on entry '
        call cc_prp(omega1,dummy,1,1,0)
      end if
               
      do isymak = 1, nsym
         nvajkl(isymak) = 0
         icount7 = 0
         do isymk = 1, nsym
            isyma = muld2h(isymak,isymk)
            ivajkl(isyma,isymk) = icount7
            nvajkl(isymak) = nvajkl(isymak)+nt1ao(isyma)*nmatkl(isymk)
            icount7 = icount7 + nt1ao(isyma)*nmatkl(isymk) 
         end do
      end do

      kctil = 1
      kvajkl = kctil + ntr12sq(isymc)
      kend1 = kvajkl + nvajkl(1)
      lwrk1 = lwork - kend1
      if (lwrk1 .lt.0) then
         call quit('Insufficient work space in ccrhs_gp')
      end if

c     --------------------
c      read V(alpha j,kl)
c     --------------------
      luvajkl = -1
      call gpopen(luvajkl,fvajkl,'unknown',' ','unformatted',
     &     idummy,.false.)
      rewind(luvajkl)
      read(luvajkl) (work(kvajkl+i-1), i = 1,nvajkl(1))
      call gpclose(luvajkl,'KEEP')

      if(locdbg) then
        write(lupri,*) 'norm^2 Valphajkl in ccrhs_gp:', 
     &                  ddot(nvajkl(1),work(kvajkl),1,work(kvajkl),1)
      end if
c     -----------------------------------
c     read r12 amplitudes or trial vector
c     -----------------------------------
      call cc_r12getct(work(kctil),isymc,iamp,ketscl,.true.,'T',
     &              lufc12,fc12am,ifile,cdummy,idummy,work(kend1),lwrk1)

c     -------------------
c     get V^aj_kl*Ctilde
c     -------------------
      isyomg = isymc ! isyomg is symmetry of result vector !
      call ccrhs_gp0(omega1,isyomg,work(kvajkl),1,lambp,1,
     &             work(kctil),isymc,etest,er12,locdbg,
     &             work(kend1),lwrk1)

      if (etest) then
         write(lupri,*)'E-R12 = ', er12
      end if
      if (locdbg) then
        write(lupri,*)'in ccrhs_gp: omega1 on exit '
        call cc_prp(omega1,dummy,1,1,0)
      end if
         
      call qexit('ccrhs_gp')
      end 
*=====================================================================*
      subroutine ccrhs_eintp(eim,lambp,work,lwork,iamp,isymc,
     &                        fc12am,lufc12,ifile,listr,idlstr)
*---------------------------------------------------------------------*
c     purpose: calculate R12 contribution to E^(R12) intermediate  
c
c     H. Fliegl, C. Haettig summer 2003 
c     generalized by C. Neiss autumn 2004
c
c     iamp      flag which indicates from where the amplitudes are to 
c		be read from; three possible values:
c		0 	ground state cluster amplitudes
c		1	amplitudes from trial vector adressed by 
c			fc12am, lufc12 and ifile
c		2	response amplitudes adressed by listr, idlstr
c     isymc     symmetry of the amplitude vecor 
c---------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <dummy.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
#include <ccr12int.h>

      logical locdbg
      parameter(locdbg = .false.)
#if defined (SYS_CRAY)
      real zero,one,two
#else
      double precision zero,one,two
#endif
      parameter (zero = 0.0d0, one = 1.0d0, two = 2.0d0)

      character*(*) fc12am, listr
      integer lwork,luvajkl,kvajkl,lwrk1,kend1
      integer idum,lu43,iamp
      integer kctil,isymkl,isymij,isymmu,isyml,isymk,
     &     idxij,idxji,idxkl,idxklij,idxklji,
     &     isymmuj,koffl,koffv,kvmujkl,ntotmu,ntotmuj,
     &     isymi,isymj,isyma,isymnj,isymn,ntotn,isymc
      integer icount7,isymak,kvint,kend2,lwrk2,ntoti
      integer koffo,koffc,index,ifile,lufc12,idlstr,idxnj
#if defined (SYS_CRAY)
      real work(*),lambp(*),eim(*),ddot, vnrm, vcprod
#else
      double precision work(*),lambp(*),eim(*),ddot, vnrm, vcprod
#endif

      index(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j

      call qenter('ccrhs_eintp')
               
      kctil = 1
      kvajkl = kctil + ntr12sq(isymc)
      kend1 = kvajkl + nvajkl(1)
      lwrk1 = lwork - kend1

      if (lwrk1 .lt.0) then
         call quit('Insufficient work space in ccrhs_gp')
      end if

c     --------------------
c      read V(alpha j,kl)
c     --------------------
        luvajkl = -1
        call gpopen(luvajkl,fvajkl,'old',' ','unformatted',
     &       idummy,.false.)
        rewind(luvajkl)
        read(luvajkl) (work(kvajkl+i-1), i = 1,nvajkl(1))
        call gpclose(luvajkl,'KEEP')

C     ---------------------------------------------------
C     read R12 amplitudes and form (2*c_ij^kl - c_ji^kl)
C     ---------------------------------------------------
       call cc_r12getct(work(kctil),isymc,iamp,ketscl,.true.,'T',
     &                  lufc12,fc12am,ifile,listr,idlstr,
     &                  work(kend1),lwrk1)
C
      if (locdbg) then
        write(lupri,*) 'Norm^2(R12 ampl.)=', ddot(ntr12sq(isymc),
     &      work(kctil),1,work(kctil),1)
        write(lupri,*) 'symmetry of R12 ampl:',isymc
        write(lupri,*) 'Norm^2(EIM) at start of ccrhs_eintp:',
     &     ddot(nmatij(isymc),eim,1,eim,1)
      end if
C

C     ------------------------------------
C     calculate E(R12)-intermediate
C     ------------------------------------

      do isymkl = 1, nsym
         isymnj  = isymkl
         isymmuj = isymkl
         isymij = muld2h(isymkl,isymc)
CCN      uncomment the following lines if locdbg = .TRUE. 
CCN            do idxkl = 1, nmatkl(isymkl)
            if (locdbg) vnrm = zero
            if (locdbg) vcprod = zero
         do isymn =1, nsym
            isymmu = isymn
            isymj  = muld2h(isymmuj,isymmu)
            isymi  = muld2h(isymij,isymj) 

           if (muld2h(isymn,isymi) .NE. isymc) then
             call quit('Symmetry mismatch in ccrhs_eintp!!')
           end if
  
            kvint = kend1
            kend2 = kvint + nrhf(isymn)*nrhf(isymj)
            lwrk2 = lwork - kend2
            if (lwrk2 .lt.0) then
               call quit('Insufficient work space in ccrhs_eint')
            end if
CCN         comment the following line if locdbg = .TRUE.
            do idxkl = 1, nmatkl(isymkl)                      
c            ---------------------------------------
c             transform Vmujkl with lambdap to Vnjkl
c            ---------------------------------------
               kvmujkl = kvajkl+ivajkl(isymmuj,isymkl)+
     &              nt1ao(isymmuj)*(idxkl-1)+it1ao(isymmu,isymj)
               ntotmu = max(1,nbas(isymmu))
               koffl = iglmrh(isymmu,isymn)+1
               ntotn = max(1,nrhf(isymn)) 

               call dgemm('T','N',nrhf(isymn),
     &              nrhf(isymj),nbas(isymmu),one,lambp(koffl),
     &              ntotmu,work(kvmujkl),ntotmu,zero,
     &              work(kvint),ntotn)
               if (locdbg) then
                 vnrm=vnrm+ddot(nrhf(isymn)*nrhf(isymj),work(kvint),1,
     &                                                  work(kvint),1)
                 koffc = kctil + itr12sqt(isymnj,isymkl)+
     &              nmatij(isymnj)*(idxkl-1) + imatij(isymn,isymj)
                 vcprod=vcprod+ddot(nrhf(isymn)*nrhf(isymj),
     &                   work(kvint),1,work(koffc),1)
                 if (idxkl.eq.2 .and. isymkl.eq.1) then
                    do idxnj = 1, nrhf(isymj)*nrhf(isymn)
                       write(lupri,*) work(kvint-1+idxnj), 
     &                                work(koffc-1+idxnj)
                    end do
                 end if     
               end if
c             ---------------------------- 
c              contract with R12 amplitudes 
c             ---------------------------               
               ntoti = max(1,nrhf(isymi)) 
               koffo = 1+imatij(isymn,isymi)
               koffc = kctil + itr12sqt(isymij,isymkl)+
     &              nmatij(isymij)*(idxkl-1) + imatij(isymi,isymj)

               call dgemm('N','T',nrhf(isymn),nrhf(isymi),nrhf(isymj),
     &                one,work(kvint),ntotn,work(koffc),ntoti,
     &                one,eim(koffo),ntotn)

            end do
            if (locdbg) then
             koffc = kctil + itr12sqt(isymij,isymkl)+
     &             nmatij(isymij)*(idxkl-1) 
             write(lupri,*) 'isymkl,idxkl,Norm^2(Vnjkl),Norm^2(cijkl):',
     &          isymkl, idxkl, vnrm,ddot(nmatij(isymij),
     &                        work(koffc),1,work(koffc),1) ,vcprod
            end if
         end do
      end do
      if (locdbg) then
        write(lupri,*) 'Norm^2(EIM) at end of ccrhs_eintp:',
     &     ddot(nmatij(isymc),eim,1,eim,1)
      end if

      call qexit('ccrhs_eintp')
      end 
*=====================================================================*
      subroutine cc_r12getct(ctilde,isymc,iamp,basscl2,do_2cme,trans,
     &                     lufc12,fc12am,ifile,listr,idlstr,work,lwork)
c---------------------------------------------------------------------
c     purpose: get ctilde = (2*c^ij_kl - c^ji_kl) OR c^ij_kl
c              format is: square matrix, upper index pair (ij) leading
c
c     H. Fliegl, C. Haettig, spring 2004
c     generalized by C. Neiss  2004/2005
c
c     ctilde:	result
c     isymc:	symmetry of result
c     iamp:     type of amplitudes used (see CCRHS_EINTP)
c     basscl2:  scale factor for diagonal elements for unpacking trial 
c               vectors into square matrix format (see CCRHS_EP)
c     do_2cme:	flag whether to calculate (2*c^ij_kl - c^ji_kl) or not
c     trans:    flag wheter lower (r12) indices (kl) are leading
C                  (-> 'N'), or upper (occ.) indices (ij) are 
C                  leading (-> 'T') in ctilde
c     lufc12:	unit number of trial vector file
c     fc12am:	file name of trial vector file
c     ifile:	index of trial vector
c     listr:	type of response vector (see CC_WRRSP.F)
c     idlstr:	index of response vector
c---------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <dummy.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h> 
#include <r12int.h> 
#include <ccr12int.h> 

      integer lwork, kend1, kpck, lwrk1, lu43, idum, 
     &        ifile, lufc12, isymc, iopt
      integer iamp,idlstr,lunit
#if defined (SYS_CRAY) 
      real work(*), ctilde(*), two, ddot, basscl2
#else
      double precision work(*), ctilde(*), two, ddot, basscl2
#endif
      parameter (two = 2.0d0)
      logical do_2cme, ldum, locdbg
      parameter (locdbg = .false.)
      character*(*) fc12am
      character     listr*(3), trans*(1)
      character     model*(10)

      CALL QENTER('cc_r12getct')

      kpck  = 1
      kend1 = kpck + ntr12am(isymc)
      lwrk1 = lwork - kend1
      if (lwrk1.lt.0) then
         write(lupri,*) 'lwork, lwrk1: ',lwork, lwrk1
         call quit('Insufficient work space in cc_r12getct')
      end if

      if (iamp .EQ. 1) then
c       --------------------
c       read R12 trial vector 
c       -------------------
        call cc_rvec(lufc12,fc12am,ntr12am(isymc),ntr12am(isymc),
     *               ifile,work(kpck))
        call cclr_diasclr12(work(kpck),basscl2,isymc)
        iopt = 1
        call ccr12unpck2(work(kpck),isymc,ctilde,trans,iopt)
      else if (iamp .EQ. 2) then
C       -------------------------
C       read response amplitudes 
C       -------------------------
        iopt = 32
        call cc_rdrsp(listr,idlstr,isymc,iopt,model,dummy,work(kpck))
        call cclr_diasclr12(work(kpck),basscl2,isymc)
        iopt = 1
        call ccr12unpck2(work(kpck),isymc,ctilde,trans,iopt)
      else if (iamp .eq. 0) then
c       --------------------------------
c       read R12 ground state amplitudes 
c       --------------------------------
        if (isymc.ne.1) call quit('Symmetry error in CC_R12GETCT')
        iopt = 32
        call cc_rdrsp('R0 ',0,1,iopt,model,dummy,work(kpck))
        iopt = 1
        call ccr12unpck2(work(kpck),isymc,ctilde,trans,iopt)
      else
        call quit('Illegal value for iamp in cc_r12getct!')
      end if

      if (do_2cme) then
        call cc_r12tcmesq(ctilde,isymc,trans,.false.)
      end if

      if (locdbg) then
        write(lupri,*) 'resorted R12 amplitudes (ctilde):'
        write(lupri,*) 'norm^2:',ddot(ntr12sq(isymc),
     &         ctilde,1,ctilde,1)
        call cc_prsqr12(ctilde,isymc,trans,1,.FALSE.)
      end if
      CALL QEXIT('cc_r12getct')
      return
      end
c=====================================================================*
      subroutine ccrhs_gp0(omega1,isyomg,vajkl,isymv,lambp,isymp,
     &                     ctilde,isymc,etest,er12,locdbg,work,lwork)
c---------------------------------------------------------------------
c     purpose: transform V with Lamdap to obtain V^aj_kl and 
c              get Omega_ai as V*Ctilde
c
c     H. Fliegl, C. Haettig, spring 2004
c---------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <dummy.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>

      integer lwork,isyomg,isymc,isymkl,isymaj,isymmuj,isyma,
     &        isymmu,isymj,isymi,isymij,kvint,kend1,kend2,
     &        lwrk2,idxkl,koffl,kvmujkl,ntotmu,ntota,
     &        ntoti,koffo,koffc,isymv,isymp
#if defined (SYS_CRAY)
      real work(*),er12, vajkl(*), ctilde(*), lambp(*), 
     &                 omega1(*), zero, one, ddot
#else 
      double precision work(*),er12, vajkl(*), ctilde(*), lambp(*), 
     &                 omega1(*), zero, one, ddot
#endif
      parameter (zero=0.0d0, one=1.0d0)
      logical etest, locdbg

      call qenter('ccrhs_gp0')

      kend1 = 1

      if (isyomg.ne.muld2h(isymc,isymv)) then
        write(lupri,*) 'isyomg,isymc,isymv: ',isyomg,isymc,isymv
        call quit('Symmetry mismatch in ccrhs_gp0')
      end if
      if (isymp.ne.1) then
        write(lupri,*)'isymp.ne.1'
        call quit('Symmetry mismatch in ccrhs_gp0')
      end if

c     -----------------------------
c     transform Vmujkl with lambdap 
c     -----------------------------
      if (etest) er12 = 0.0d0
      do isymkl = 1, nsym
         isymaj  = muld2h(isymkl,isymv)
         isymmuj = isymaj
         do isyma =1, nsym
            isymmu = muld2h(isymp,isyma)
            isymj  = muld2h(isymmuj,isymmu) 
            isymij = muld2h(isymc,isymkl)
            isymi = muld2h(isymij,isymj)

            kvint = kend1
            kend2 = kvint + nvir(isyma)*nrhf(isymj)
            if (etest) then 
               kend2 = kvint+max(nvir(isyma),nrhf(isymi))*nrhf(isymj)
            end if

            lwrk2 = lwork - kend2
            if (lwrk2 .lt.0) then
               call quit('Insufficient work space in ccrhs_gp0')
            end if
                         
            do idxkl = 1, nmatkl(isymkl)
               koffl = iglmvi(isymmu,isyma)+1
               kvmujkl = 1+ivajkl(isymmuj,isymkl)+
     &              nt1ao(isymmuj)*(idxkl-1)+it1ao(isymmu,isymj)
               ntotmu = max(1,nbas(isymmu))
               ntota = max(1,nvir(isyma)) 

               call dgemm('T','N',nvir(isyma),
     &              nrhf(isymj),nbas(isymmu),one,lambp(koffl),
     &              ntotmu,vajkl(kvmujkl),ntotmu,zero,
     &              work(kvint),ntota)
               
               if (locdbg) then
                 write(lupri,*)'in ccrhs_gp0: V*Lambdap'
                 call output(work(kvint),1,nvir(isyma),1,
     &              nrhf(isymj),nvir(isyma),
     &              nrhf(isymj),1,lupri)
                 write(lupri,*)'in ccrhs_gp0: Lambdap'
                 call output(lambp(koffl),1,nbas(isymmu),1,nvir(isyma),
     &                       nbas(isymmu),nvir(isyma),1,lupri)
               end if
c              ------------------------------
c               make ctilde*V(atilde,m,k,l) 
c              ------------------------------

               ntoti = max(1,nrhf(isymi)) 
               koffo = 1 + it1am(isyma,isymi)
               koffc = 1 + itr12sqt(isymij,isymkl)+
     &              nmatij(isymij)*(idxkl-1) + imatij(isymi,isymj)
                 call dgemm('N','T',nvir(isyma),nrhf(isymi),nrhf(isymj),
     &                one,work(kvint),ntota,ctilde(koffc),ntoti,
     &                one,omega1(koffo),ntota)

                 if (locdbg) then
                   write(lupri,*)'in ccrhs_gp0: ctilde',idxkl
                   write(lupri,*)'isymij,isymkl,isymi,isymj:',
     &                        isymij,isymkl,isymi,isymj
                   call output(ctilde(koffc),1,nrhf(isymi),1,
     &                    nrhf(isymj),nrhf(isymi),nrhf(isymj),1,lupri)
                   write(lupri,*)'in ccrhs_gp0: omega1 + V*ctilde'
                   call output(omega1(koffo),1,nvir(isyma),1,
     &                nrhf(isymi),nvir(isyma),
     &                nrhf(isymi),1,lupri)
                 end if

c               ------------------------------------------------------
c               make Vijkl and contract ct*Vijkl to obtain MP2-R12 
c               energy in the first iteration (when t1 = 0)  
c               this is a control if the scaling factors of the R12 
c               amplitudes ct is correct  
c              -------------------------------------------------------
               if (etest) then
                  koffl = iglmrh(isymmu,isymi)+1
                  ntoti = max(1,nrhf(isymi)) 
                  call dgemm('T','N',nrhf(isymi),
     &                 nrhf(isymj),nbas(isymmu),one,lambp(koffl),
     &                 ntotmu,vajkl(kvmujkl),ntotmu,zero,
     &                 work(kvint),ntoti)
                  
                   er12 =er12 + ddot(nrhf(isymi)*nrhf(isymj),
     &                 work(kvint),1,ctilde(koffc),1)
               end if
               
            end do
         end do
      end do
      call qexit('ccrhs_gp0')
      return
      end 
c=====================================================================*

*=====================================================================*
      subroutine cc_r12tcmesq(MATRIX,ISYM,TRANS,LINV)
c---------------------------------------------------------------------
c     purpose: calculate ctilde = (2*c^ij_kl - c^ji_kl)
c
c     MATRIX   input and output, dimension ntr12sq(ISYM)
c     LINV     do the inverse operation  
c     TRANS    r12 index pair leading: 'N'
c              occ. index pair leading: 'T'      
c     
c     Christian Neiss  Feb. 2005 
c---------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <dummy.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>

      logical locdbg,linv
      parameter(locdbg = .false.)

      character*1 trans
      integer isym
      integer isymij,isymkl,isymi,isymj,isymk,isyml
      integer idxij,idxkl,idxji,idxklij,idxklji

#if defined (SYS_CRAY)
      real matrix(*),two,xklij,xklji,fac,one
#else
      double precision matrix(*),two,xklij,xklji,fac,one
#endif
      parameter (one = 1.0D0, two = 2.0d0)

      call qenter('cc_r12tcmesq')
      if (locdbg) then
        write(lupri,*) 'Entered CC_R12TCMESQ'

        call around('Input in CC_R12TCMESQ')
        call cc_prsqr12(matrix,isym,trans,1,.false.)
      end if

      fac = -one  
      if (LINV) fac = one

      do isymkl = 1, nsym
        isymij = muld2h(isymkl,isym)
        do isymk = 1, nsym
          isyml = muld2h(isymk,isymkl)
          do k = 1, nrhfb(isymk)
            do l = 1, nrhfb(isyml)
              idxkl = imatkl(isymk,isyml)+nrhfb(isymk)*(l-1)+k
              do isymi = 1, nsym
                isymj = muld2h(isymi,isymij)
                do i = 1, nrhf(isymi)
                  do j = 1, nrhf(isymj)
                    idxij = imatij(isymi,isymj)+nrhf(isymi)*(j-1)+i
                    idxji = imatij(isymj,isymi)+nrhf(isymj)*(i-1)+j
                    if (trans.eq.'T') then
                      idxklij = itr12sqt(isymij,isymkl)+
     &                          nmatij(isymij)*(idxkl-1)+idxij
                      idxklji = itr12sqt(isymij,isymkl)+
     &                          nmatij(isymij)*(idxkl-1)+idxji
                    else if (trans.eq.'N') then
                      idxklij = itr12sq(isymkl,isymij)+
     &                          nmatkl(isymkl)*(idxij-1)+idxkl
                      idxklji = itr12sq(isymkl,isymij)+
     &                          nmatkl(isymkl)*(idxji-1)+idxkl
                    else
                      call quit('Illegal value for "TRANS" in '//
     &                          'CC_R12TCMESQ') 
                    end if
                    if (idxij.le.idxji) then
                      xklij = two*matrix(idxklij) + fac*matrix(idxklji)
                      xklji = two*matrix(idxklji) + fac*matrix(idxklij)
                      matrix(idxklij) = xklij
                      matrix(idxklji) = xklji
                    end if
                  end do
                end do
              end do
            end do
          end do
        end do
      end do

      if (LINV) then
        call dscal(ntr12sq(isym),1.0D0/3.0D0,matrix,1)
      end if
                    
      if (locdbg) then
        call around ('Result in CC_R12TCMESQ')
        call cc_prsqr12(matrix,isym,trans,1,.false.)

        write(lupri,*) 'Leaving CC_R12TCMESQ'
      end if
      call qexit('cc_r12tcmesq')
      return
      end
*======================================================================
       subroutine cc_r12tcmepk(MATRIX,ISYM,LINV)
C----------------------------------------------------------------------
c     purpose: calculate ctilde = (2*c^ij_kl - c^ji_kl)
c
c     MATRIX   input and output, dimension ntr12am(ISYM)
c     LINV     do the inverse operation
c
c     based on CCSD_TCMEPK
c     Christian Neiss  April 2005
C----------------------------------------------------------------------

      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>

      integer isym,isymi,isymj,isymk,isyml
      integer isymij,isymkl,isymki,isymlj,isymli,isymkj
      integer nrhfk,nrhfi,nli,nlj,nki,nkj,nkilj,nkjli 
      logical linv

      logical locdbg
      PARAMETER (locdbg = .FALSE.)      

#if defined (SYS_CRAY)
      real precision matrix(*), xkilj, xkjli, fac, one, two
#else
      double precision matrix(*), xkilj, xkjli, fac, one, two
#endif
      integer index 
      parameter (one = 1.0D0, two = 2.0D0)

      index(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j  
      
      call qenter('cc_r12tcmepk')

      if (locdbg) then
        write(lupri,*) 'Entered CC_R12TCMEPK'
        call around('Input in CC_R12TCMEPK')
        call cc_prpr12(matrix,isym,1,.false.)
      end if

      fac = -one
      if (LINV) fac = one
C
      DO 100 ISYMIJ = 1,NSYM
C
         ISYMKL = MULD2H(ISYMIJ,ISYM)
C
         DO 110 ISYMJ = 1,NSYM
C
            ISYMI = MULD2H(ISYMJ,ISYMIJ)
C
            IF (ISYMI .GT. ISYMJ) GOTO 110
C
            DO 120 ISYML = 1,NSYM
C
               ISYMK = MULD2H(ISYML,ISYMKL)
C
               IF (ISYMK .GT. ISYML) GOTO 120
C
               ISYMKI = MULD2H(ISYMK,ISYMI)
               ISYMLJ = MULD2H(ISYML,ISYMJ)
               ISYMLI = MULD2H(ISYML,ISYMI)
               ISYMKJ = MULD2H(ISYMK,ISYMJ)
C
               DO 130 J = 1,NRHF(ISYMJ)
C
                  IF (ISYMI .EQ. ISYMJ) THEN
                     NRHFI =  J
                  ELSE
                     NRHFI = NRHF(ISYMI)
                  ENDIF
C
                  DO 140 I = 1,NRHFI
C
                     DO 150 L = 1,NRHFB(ISYML)
C
                        IF (ISYML .EQ. ISYMK) THEN
                           NRHFK = L
                        ELSE
                           NRHFK = NRHFB(ISYMK)
                        ENDIF
C
                        NLI = IMATKI(ISYML,ISYMI)
     *                      + NRHFB(ISYML)*(I - 1) + L
                        NLJ = IMATKI(ISYML,ISYMJ)
     *                      + NRHFB(ISYML)*(J - 1) + L
C
                        DO 160 K = 1,NRHFK
C
                           NKI = IMATKI(ISYMK,ISYMI)
     *                         + NRHFB(ISYMK)*(I - 1) + K
                           NKJ = IMATKI(ISYMK,ISYMJ)
     *                         + NRHFB(ISYMK)*(J - 1) + K
C
                         IF (ISYMKI.EQ.ISYMLJ) THEN
                           NKILJ = ITR12AM(ISYMKI,ISYMLJ)
     *                           + INDEX(NKI,NLJ)
                         ELSE IF (ISYMKI.LT.ISYMLJ) THEN
                           NKILJ = ITR12AM(ISYMKI,ISYMLJ)
     *                           + NMATKI(ISYMKI)*(NLJ-1)+NKI
                         ELSE IF (ISYMLJ.LT.ISYMKI) THEN
                           NKILJ = ITR12AM(ISYMKI,ISYMLJ)
     *                           + NMATKI(ISYMLJ)*(NKI-1)+NLJ
                         END IF
C
                         IF (ISYMKJ.EQ.ISYMLI) THEN
                           NKJLI = ITR12AM(ISYMKJ,ISYMLI)
     *                           + INDEX(NKJ,NLI)
                         ELSE IF (ISYMKJ.LT.ISYMLI) THEN
                           NKJLI = ITR12AM(ISYMKJ,ISYMLI)
     *                           + NMATKI(ISYMKJ)*(NLI-1)+NKJ
                         ELSE IF (ISYMLI.LT.ISYMKJ) THEN
                           NKJLI = ITR12AM(ISYMKJ,ISYMLI)
     *                           + NMATKI(ISYMLI)*(NKJ-1)+NLI
                         END IF
C
                         XKILJ = TWO*MATRIX(NKILJ) + FAC*MATRIX(NKJLI)
                         XKJLI = TWO*MATRIX(NKJLI) + FAC*MATRIX(NKILJ)
C
                         MATRIX(NKILJ) = XKILJ
                         MATRIX(NKJLI) = XKJLI
C
  160                   CONTINUE
  150                CONTINUE
  140             CONTINUE
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C 
      if (LINV) then
        call dscal(ntr12am(isym),1.0D0/3.0D0,matrix,1)
      end if
C    
      if (locdbg) then
        call around('Result in CC_R12TCMEPK')
        call cc_prpr12(matrix,isym,1,.false.)
        write(lupri,*) 'Leaving CC_R12TCMEPK'
      end if
      call qexit('cc_r12tcmepk')
      return

      end 
*======================================================================*
      subroutine cc_r12mkSak(work,lwork) 
c----------------------------------------------------------------------
c     purpose: calculate overlap S_ak
c              a        AO index
c              k        r12 index 
c  
c     C. Neiss  august 2005
c----------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
#include <ccr12int.h>
#include <dummy.h>
 
      logical locdbg
      parameter (locdbg = .false.)
#if defined (SYS_CRAY)
      real zero,one
#else
      double precision zero,one
#endif
      parameter (zero = 0.0d0, one = 1.0d0)

      integer nr1orb(8),nr1bas(8),nr1xbas(8),n2bst1(8),nrgkl(8)
      integer nr2bas(8)
      integer ir1orb(8,8),ir1bas(8,8),ir2bas(8,8),ir1xbas(8,8)
      integer ir2xbas(8,8),iaodis1(8,8),irgkl(8,8)
      integer nr1xorb(8),ir1xorb(8,8),nrxgkl(8),irxgkl(8,8)
      integer igabj(8),nalphaj(8),ialphaj(8,8)
      integer isymab,isyma,isymb,isymk,isym0
      integer idxab,idxak,koffc,ntotb,ntota
      integer koverlp,ksab,ksak,kend1,lwrk1,lwork,koff1,koffs
      integer nnbast,nlamdsx,idum
      integer kcmo,lusifc,lusak
#if defined (SYS_CRAY)
      real work(*)
#else
      double precision work(*)
#endif
      
      call qenter('cc_r12mkSak')

      call cc_r12offset(nr1orb,nr1xorb,nr1bas,nr1xbas,nr2bas,
     &     nrgkl,nrxgkl,n2bst1,ir1orb,ir1xorb,ir1bas,ir1xbas,ir2bas,
     &     ir2xbas,irgkl,irxgkl,iaodis1,nalphaj,ialphaj)

      nnbast = 0
      do i = 1,nsym
         nnbast = nnbast + nbas(i)*(nbas(i)+1)/2
      end do

      koverlp = 1
      ksab = koverlp + n2bst(1)
      ksak = ksab + n2bst1(1)
      kcmo = ksak + nr1bas(1)  
      kend1 = kcmo + nlamds
      lwrk1 = lwork - kend1
      if (lwrk1.lt.nnbast) then
         write(lupri,*)'lwrk1,lwork',lwrk1,lwork
         call quit('Insufficient work space in cc_r12mkSak.')
      end if

c    --------------------------------------------------------------
c     read overlap matrix from file and unpack as quadratic matrix
c    --------------------------------------------------------------
      call rdonel('OVERLAP ',.true.,work(kend1),nnbast)
      isym0 = 1
      call ccsd_symsq(work(kend1),isym0,work(koverlp))

c     -------------------------------------------
c     read MO-coefficients from interface file
c     -------------------------------------------
      LUSIFC = -1
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
      CALL MOLLAB('FULLBAS ',LUSIFC,LUPRI)
      READ (LUSIFC) IDUM,IDUM,IDUM,NLAMDSX
      IF (NLAMDSX.NE.NLAMDS) THEN
        CALL QUIT('CC_R12MKSAK assumes MP2-R12 environment!')
      END IF 
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KCMO+I-1), I=1,NLAMDS)
      CALL GPCLOSE(LUSIFC,'KEEP')

      CALL CMO_REORDER(WORK(KCMO),WORK(KEND1),LWRK1)

      do isymk = 1, nsym
        isyma = isymk
        isymb = isymk
        do b = 1, mbas1(isymb)
          do a = 1, mbas1(isyma)
             koff1 = iaodis(isymb,isyma)+nbas(isyma)*(b-1)+a
             idxab = mbas1(isyma)*(b-1)+a
             work(ksab-1+idxab) = work(koverlp -1 + koff1)
          end do
        end do
c       -------------------------------------------
c       transform S_{alpha beta} with C_{beta k} 
c       to obtain S_{alpha k}
c       -------------------------------------------            
        koffc = kcmo + iglmrh(isymb,isymk) 
        koffs = ksak + ir1bas(isyma,isymk)
        ntotb = max(nbas(isymb),1)
        ntota = max(mbas1(isyma),1)
         
        call dgemm('N','N',mbas1(isyma),
     &       nrhfb(isymk),mbas1(isymb),one,work(ksab),ntota,
     &       work(koffc),ntotb,zero,work(koffs),ntota)
      end do
  
      if (locdbg) then
        write(lupri,*)'Sak in CC_R12MKSAK ='
        do isyma = 1, nsym
          isymk = isyma
          call output(work(ksak+ir1bas(isyma,isymk)),1,mbas1(isyma),
     &               1,nrhfb(isymk),
     &               mbas1(isyma),nrhfb(isymk),1,lupri)
        end do
      end if

c     ----------------------------------
c     write result to file
c     ----------------------------------
      lusak = -1
      call gpopen(lusak,fnsak,'UNKNOWN',' ','UNFORMATTED',
     &            idummy,.false.)
      rewind(lusak)
      write(lusak) (work(ksak+i-1), i=1, nr1bas(1))
      call gpclose(lusak,'KEEP')

      call qexit('cc_r12mkSak')
      return
      end       
*======================================================================*
      subroutine cc_r12getvxykl(vxykl,length,cmo,isycmo,
     &                          nr1bas,ir1bas,ir2bas,
     &                          nalphaj,ialphaj,
     &                          lunit,filename,iopt,work,lwork)
c---------------------------------------------------------------------
c     purpose: get Vabkl for initialisation routine cc_r12mkvamkl0
c              needed for extension of new correlation factor f12
c
c     H. Fliegl, C. Haettig, autumn 2005
c---------------------------------------------------------------------
      implicit none
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
#include <dummy.h>
#include <priunit.h>

#if defined (SYS_CRAY)
      real vxykl(*),work(*),ddot,one,zero,cmo(*)
#else
      double precision vxykl(*),work(*),ddot,one,zero,cmo(*)
#endif
      parameter (one = 1.0d0, zero = 0.0d0)

      character*8 filename
      integer lunit,isyma,isymb,isymk,isyml,isymab,isymkl,isymak,
     &        nbetl,iadr,ilen,idxkl,idxak,kend2,kscr2,kend1,lwrk2 
      integer lwork,isymbl,idxab,idxabkl,iopt,idxbl,idxbi,idxil
      integer nr1bas(8),ir1bas(8,8),ir2bas(8,8)
      integer idxbjlk,kend3,lwrk3,ntota,koffv,koffc,length
      integer isymj,kscr3,isycmo,isymbj,nalphaj(8),ialphaj(8,8)
      integer ntota1,koffvblak,idxbj,isymi,isymij
      integer kscr4,kend4,lwrk4,idxjilk
      integer ntotk,idxkj,idxlk,isymkj,idxij


      call qenter('getvxykl') 

      lunit = -1
      call wopen2(lunit,filename,64,0)

      kend1 = 1

      call dzero(vxykl,length)

      do isymb = 1, nsym
       do b = 1, mbas1(isymb)
        do isyml = 1, nsym
         isymbl = muld2h(isymb,isyml)
         isymak = isymbl

         kscr2 = kend1
         kend2 = kscr2 + nr1bas(isymak)
         lwrk2 = lwork - kend2
         if (lwrk2.lt.0) then
            call quit('Insufficient work space in getvxykl')
         end if

         do l = 1, nrhfb(isyml)
          nbetl = ir1bas(isymb,isyml) + mbas1(isymb)*(l-1)+b
          iadr  = ir2bas(isymak,isymbl)+nr1bas(isymak)*(nbetl-1)+1
          ilen  = nr1bas(isymak)
          call getwa2(lunit,filename,work(kscr2),iadr,ilen)

          do isyma = 1, nsym
           isymab = muld2h(isyma,isymb)
           isymk  = muld2h(isymak,isyma)
           isymkl = muld2h(isymk,isyml)

c          isymbj = muld2h(isymkl,isycmo)
c          isymj  = muld2h(isycmo,isyma)
c          isymkj = muld2h(isymk,isymj)

c          isymi = muld2h(isycmo,isymb)
c          isymij = muld2h(isymi,isymj)
        
           if (iopt.eq.0) then
c            calculate Valpha_beta_kl

             do k = 1, nrhfb(isymk)
               idxkl = imatkl(isymk,isyml)+nrhfb(isymk)*(l-1)+k
               do a = 1, mbas1(isyma)
                 idxak = ir1bas(isyma,isymk) + mbas1(isyma)*(k-1)+a
                 idxab = iaodis(isyma,isymb) + mbas1(isyma)*(b-1)+a
                 idxabkl = ivabkl(isymab,isymkl)+
     &                        n2bst(isymab)*(idxkl-1)+idxab
                 vxykl(idxabkl) = work(kscr2 -1 + idxak)
               end do
             end do
           else if ((iopt.eq.1).or.(iopt.eq.2)) then
             isymbj = muld2h(isymkl,isycmo)
             isymj  = muld2h(isycmo,isyma)
             isymkj = muld2h(isymk,isymj)
             
             isymi = muld2h(isycmo,isymb)
             isymij = muld2h(isymi,isymj)

c            calculate Valpha_j_kl

             kscr3 = kend2
             kend3 = kscr3 + nmatki(isymkj)
             lwrk3 = lwork - kend3
             if (lwrk3.lt.0) then
               call quit('insufficient work space in getvxykl')
             end if

             ntota = max(nbas(isyma),1) 
             ntota1 = max(mbas1(isyma),1) 
             ntotk = max(nrhfb(isymk),1)

             koffc = iglmrh(isyma,isymj)+1
             koffvblak = kscr2 + ir1bas(isyma,isymk)
             koffv = kscr3 + imatki(isymk,isymj) 

             call dgemm('T','N',nrhfb(isymk),nrhfa(isymj),mbas1(isyma),
     &                  one,work(koffvblak),ntota1,cmo(koffc),ntota,
     &                  zero,work(koffv),ntotk)
             
             if (iopt.eq.2) then
c              calculate Vij_kl
              
               do i = 1, nrhfa(isymi)
                 idxbi = iglmrh(isymb,isymi)+mbas1(isymb)*(i-1)+b
                 do k = 1, nrhfb(isymk)
                   idxlk = imatkl(isyml,isymk)+nrhfb(isyml)*(k-1)+l
                   do j = 1, nrhfa(isymj)
                     idxij = imatij(isymi,isymj)+nrhfa(isymi)*(j-1)+i
                     idxkj = imatki(isymk,isymj)+nrhfb(isymk)*(j-1)+k
                     idxjilk = itr12sqt(isymij,isymkl)+
     &                         nmatij(isymij)*(idxlk-1)+idxij

                     vxykl(idxjilk) = vxykl(idxjilk)+
     &                               work(kscr3-1+idxkj)*cmo(idxbi)

                   end do
                 end do
               end do

             else
c              sort Vblkj as Vbjlk             
               do j = 1, nrhfa(isymj)
                 idxbj = ialphaj(isymb,isymj)+mbas1(isymb)*(j-1)+b
                 do k = 1, nrhfb(isymk)
                   idxlk = imatkl(isyml,isymk)+nrhfb(isyml)*(k-1)+l
                   idxkj = imatki(isymk,isymj)+nrhfb(isymk)*(j-1)+k
                   idxbjlk = ivajkl(isymbj,isymkl)+
     &                       nalphaj(isymbj)*(idxlk-1)+idxbj
                   vxykl(idxbjlk) = work(kscr3 -1 +idxkj)
                 end do
               end do
c
             end if
c            
           else
            call quit('This iopt is not implemented yet')
           end if
c
          end do

         end do
        end do
       end do
      end do
       
      call wclose2(lunit,filename,'KEEP')

      call qexit('getvxykl') 
      end
c====================================================================*
