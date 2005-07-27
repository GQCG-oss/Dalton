c*=====================================================================*
      subroutine r12mkvamkl(filback,gabjd,tabjd,vajkl,cmo,isymcmo,
     &                      xlamdhs,xlamdps,xint,
     &                      tint,idel,isymd,isymj,isymab,isymg,scr1,
     &                      ibasx,iglmrhs,nglmds,lauxd,work,lwork)
*----------------------------------------------------------------------*
c     Purpose: calculation of V^{alpha,j}_{k,l} 
c     
c     where:     alpha      AOs in orbital basis
c                j,k,l      active occupied MOs
c
c     H. Fliegl, C. Haettig spring 2003
c
c     C. Neiss 2004/2005: adapted for non-total-symm. V (Ansatz 1)
c                         !Ansatz 2 still needs to be checked!
c     isymcmo    Symmetry of cmo-matrix (needed if cmo=lambda_bar)
c
c     adapted for ansatz 3, Christof Haettig, spring 2005.
*----------------------------------------------------------------------*
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
#include <ccr12int.h>
#include <dummy.h>

      logical locdbg,lauxd
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
      integer koffrtf,imaijm(8,8),nmaijm(8),isymijm,icount1,isymij
      integer nvabkl(8),ivabkl(8,8)
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
      if (lauxd) idelta = idelta -ibasx(isymd)

c    ----------------------------------------------------
c     for ansatz 2 skip loop over deltas in orbital basis 
c    ----------------------------------------------------
      if (ianr12.eq.2.and. idelta.le.mbas1(isymd)) return
      
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
     &     ir2xbas,irgkl,irxgkl,iaodis1,nvabkl,ivabkl)

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

         ngabj  = ngabj + mbas1(isymb)*mbas1(isyma)*nrhf(isymj)
         ngmaj  = ngmaj + nrhfs(isymb)*mbas1(isyma)*nrhf(isymj)
         ngmgi  = ngmgi + nrhfs(isymb)*mbas1(isymg)*nrhf(isymi)
      end do
      
c    
      do isymijm = 1, nsym
        icount1 = 0
        do isymm = 1, nsym
          isymij = muld2h(isymijm,isymm)
          imaijm(isymij,isymm) = icount1
          icount1 = icount1 + nmatij(isymij)*nrhfs(isymm)
        end do
        nmaijm(isymijm) = icount1
      end do
      
      kr2bkl = 1
      kend0  = kr2bkl + nrgkl(isymd)
      if (iansatz.eq.2) then
        krtf = kend0
        kend0 = krtf + nmaijm(isymd)   
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
        write(lupri,*)'Vajkl entered mkvamkl'
        write(lupri,*)'norm^2:',ddot(nvajkl(isymv),vajkl,1,vajkl,1)
c       write(lupri,*)(vajkl(i),i=1,nvajkl(isymv))
        do isymkl = 1, nsym
           isymaj = muld2h(isymv,isymkl)
           write(lupri,*) 'Block isymaj, isymkl: ',isymaj,isymkl
           if ((nt1ao(isymaj).eq.0).or.(nmatij(isymkl).eq.0)) then
            write(lupri,*) 'This block is empty'
           else
            call output(vajkl(1+ivajkl(isymaj,isymkl)),1,
     &                  nt1ao(isymaj),1,nmatij(isymkl),nt1ao(isymaj),
     &                  nmatij(isymkl),1,lupri)
           end if
         end do
      end if

c     -------------------------------------
c     get R12 integrals R^{delta}(beta,kl)
c     -------------------------------------
      call cc_r12getrint(work(kr2bkl),idel,isymd,nr1bas,ir1bas,
     &     nr2bas,ir2bas,nrgkl,irgkl,ir1xbas,ir2xbas,
     &     ibasx,lauxd,.false.,filback,work(kend1),lwrk1)
      
      if (locdbg) then
        write(lupri,'(a,l5,i5,/,8i7)')
     &           'lauxd,idel,ibasx:',lauxd,idel,ibasx
        write(lupri,*)'Rdp_a,kl in mkvamkl',idel,
     &    ddot(nrgkl(isymd),work(kr2bkl),1,work(kr2bkl),1)
c       write(lupri,*) (work(kr2bkl+i-1),i=1,nrgkl(isymd))
        do isymb = 1, nsym
          isymkl = muld2h(isymb,isymd)
          write(lupri,*) 'isymb,isymkl:',isymb,isymkl
          call output(work(kr2bkl+irgkl(isymb,isymkl)),1,mbas1(isymb),
     &                1,nmatij(isymkl),
     &                mbas1(isymb),nmatij(isymkl),1,lupri)
        end do
      end if

c     for ansatz 2 transform beta index to all occupied (m)
      if (iansatz.eq.2) then
        do isymb = 1, nsym
          isymm  = isymb
          isymkl = muld2h(isymb,isymd)
          koffc  = 1+iglmrhs(isymb,isymm) 
          koffr  = kr2bkl+irgkl(isymb,isymkl)
          koffrtf = krtf+imaijm(isymkl,isymm)  !for nrhfs(isymm)
      
          ntotb = max(mbas1(isymb),1)
          ntotc = max(nbas(isymb),1)
          ntotm = max(nrhfs(isymm),1)

          call dgemm('T','N',nrhfs(isymm),nmatij(isymkl),mbas1(isymb),
     &               one,xlamdhs(koffc),ntotc,work(koffr),ntotb,zero,
     &               work(koffrtf),ntotm)

          x1 = x1 + ddot(nrhfs(isymm)*mbas1(isymb),
     &                   xlamdhs(koffc),1,xlamdhs(koffc),1)
          x2 = x2 + ddot(nrgkl(isymd),work(koffr),1,work(koffr),1)
        end do


        if (locdbg) then
          write(lupri,*)'nmaijm(isymd) =', nmaijm(isymd)
          write(lupri,*)'Norm R_m,kl in mkvamkl',idel,
     &      ddot(nmaijm(isymd),work(krtf),1,work(krtf),1)
          x3 = x3 + ddot(nmaijm(isymd),work(krtf),1,work(krtf),1)
          write(lupri,*)'R_m,kl', (work(krtf-1+i),i=1,nmaijm(isymd))
        end if
      end if

c     ----------------------------------------------------------------
c      for ansatz 1: get g_{alpha beta j delta}, i.e. unpack 
c                    triangular matrix g to quadratic matrix (scr1) 
c                    and remove auxiliary basis functions alpha,beta 
c      for ansatz 2: transform beta index to occupied
c     ----------------------------------------------------------------
      do j=1, nrhf(isymj)
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
              if (lwrk1.lt.nrhfs(isymm)*mbas1(isyma)) then
                call quit('Insufficient work space in mkvamkl')
              end if
            
              ntotm = max(nrhfs(isymm),1)
              ntota = max(nbas(isyma),1)
              ntotb = max(nbas(isymb),1)
c             get G^dp_m,alpha j
              call dgemm('T','T',nrhfs(isymm),mbas1(isyma),mbas1(isymb),
     &                   one,xlamdps(kcmo),ntotb,scr1(koff1),ntota,zero,
     &                   work(kgmaj),ntotm)

              ! sort into array aranged as G(m,aj)
              do a = 1, mbas1(isyma)
                 idxaj = mbas1(isyma)*(j-1)+a
                 do m = 1, nrhfs(isymm)
                    koff1  = kgmaj-1 + nrhfs(isymm)*(a-1) + m
                    idxmaj = igmaj(isymaj)+nrhfs(isymm)*(idxaj-1)+m  
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
     &   ddot(nrhf(isymj)*nnbst(isymab),gabjd,1,gabjd,1)
        write(lupri,*) 'resorted coulomb integrals: ... norm^2=',
     &   ddot(ngabj,work(kgabj),1,work(kgabj),1)
        do isymb = 1, nsym
          isyma  = muld2h(isymab,isymb)
          isymaj = muld2h(isyma,isymj)
          write(lupri,*) 'isymj,isyma,isymb:',isymj,isyma,isymb
          call output(work(kgabj+igabj(isymaj)),1,mbas1(isymb),
     &                1,mbas1(isyma)*nrhf(isymj),
     &                mbas1(isymb),mbas1(isyma)*nrhf(isymj),1,lupri)
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
           koffv = 1 + ivajkl(isymaj,isymkl) + ir1bas(isyma,isymj)
           
           ntotb  = max(1,mbas1(isymb))
           ntotaj = max(1,nr1bas(isymaj))
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
             call output(vajkl(koffv),1,mbas1(isyma)*nrhf(isymj),1,
     &          nmatij(isymkl),mbas1(isyma)*nrhf(isymj),nmatij(isymkl),
     &          1,lupri)
           end if

           call dgemm('T','N',mbas1(isyma)*nrhf(isymj),nmatij(isymkl),
     &          mbas1(isymb),factor,work(koffg),ntotb,
     &          work(koffr),ntotb,one,vajkl(koffv),ntotaj)
          
           if (locdbg) then
             write(lupri,*)'d: G'
             call output(work(koffg),1,mbas1(isymb),
     &                1,mbas1(isyma)*nrhf(isymj),
     &                mbas1(isymb),mbas1(isyma)*nrhf(isymj),1,lupri)
             write(lupri,*)'d: R'
             call output(work(koffr),1,mbas1(isymb),
     &                1,nmatij(isymkl),
     &                mbas1(isymb),nmatij(isymkl),1,lupri)
             write(lupri,*)'d: V = V + G*R'
             call output(vajkl(koffv),1,mbas1(isyma)*nrhf(isymj),1,
     &          nmatij(isymkl),mbas1(isyma)*nrhf(isymj),nmatij(isymkl),
     &          1,lupri)
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
          koffr = krtf + imaijm(isymkl,isymm)
          koffv = 1 + ivajkl(isymaj,isymkl) + ir1bas(isyma,isymj)

          ntotg = max(nrhfs(isymm),1)
          ntotaj = max(nr1bas(isymaj),1)

          call dgemm('T','N',mbas1(isyma)*nrhf(isymj),nmatij(isymkl),
     &        nrhfs(isymm),factor,work(koffg),ntotg,work(koffr),ntotg,
     &        one,vajkl(koffv),ntotaj)
        end do 
      else
        call quit('illegal case.')
      end if
      if (locdbg.and.iansatz.eq.2) then
        write(lupri,*)'Vajkl ansatz 2 auxbas first part:'
        call output(vajkl,1,nvajkl(1),1,nvajkl(1),nvajkl(1),
     &              nvajkl(1),1,lupri) 
      end if

c     ---------------------------------------------
c      case: delta element of auxiliary basis
c      calculate some new offsets
c     ---------------------------------------------
      if (idelta.gt.mbas1(isymd) .and. mbas1(isymg).gt.0) then
         ngbgi = 0
         do isymb = 1, nsym  
            isyma = muld2h(isymab,isymb)
            isymi = muld2h(isymcmo,isyma)
            isymgi = muld2h(isymg,isymi)
            igbgi(isymgi) = ngbgi
            ngbgi  = ngbgi + mbas1(isymb)*mbas1(isymg)*nrhf(isymi)
         end do
         kgab = kend0
         kgbi = kgab + nbast*nbast
         kgbgi = kgbi + nbast*nrhft
         krblk = kgbgi + ngbgi
         kend1 = krblk + nrgkl(isymd)

         if (iansatz.eq.2) then
           kgtf    = kend1
           krdpmlk = kgtf + ngmgi 
           kend1   = krdpmlk + nmaijm(isymd)
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
     &              nrhf(isymi),mbas1(isyma),one,work(kgab),ntotb,
     &              cmo(koffc),ntota,zero,work(kgbi),ntotb)

 

                 
               do i =1, nrhf(isymi)
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
             ntotm = max(nrhfs(isymm),1)
           
             call dgemm('T','N',nrhfs(isymm),mbas1(isymg)*nrhf(isymi),
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
     &                   1,mbas1(isymg)*nrhf(isymi),
     &                   mbas1(isymb),mbas1(isymg)*nrhf(isymi),1,lupri)
           end do
         else if (locdbg.and.iansatz.eq.2) then
            write(lupri,*)'resorted coulomb integrals (second part):',
     &       ddot(ngbgi,work(kgbgi),1,work(kgbgi),1)
            write(lupri,*)'Norm: G^dp_m,gammai in mkvamkl',
     &         ddot(ngmgi,work(kgtf),1,work(kgtf),1)
            write(lupri,*)'G^dp_m,gammai in mkvamkl',
     &                (work(kgtf-1+i),i=1,ngmgi)
            write(lupri,*)'Norm Vajkl ansatz 2:',
     &          ddot(nvajkl(1),vajkl,1,vajkl,1)
         end if
c
         if (iansatz.eq.1) then
c          -------------------------------------------------
c          pack R^{delta'}(beta,kl) as R^{delta'}(beta,lk)           
c          -------------------------------------------------
           do isyml = 1, nsym
              isymdl = muld2h(isymd,isyml)
              isymbk = isymdl
              do l=1,nrhf(isyml)
                 do isymb = 1, nsym
                    isymk  = muld2h(isymbk,isymb)
                    isymkl = muld2h(isymk,isyml)
                    do k=1,nrhf(isymk)
                       do b=1, mbas1(isymb)
                          idxkl=imatij(isymk,isyml)+nrhf(isymk)*(l-1)+k
                          idxlk=imatij(isyml,isymk)+nrhf(isyml)*(k-1)+l
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
     &                    mbas1(isymb),1,nmatij(isymkl),
     &                    mbas1(isymb),nmatij(isymkl),1,lupri)
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
             do l = 1, nrhf(isyml)
               do isymm = 1, nsym
                 isymk = muld2h(isymmk,isymm)
                 isymkl = muld2h(isymk,isyml)
                 do k = 1, nrhf(isymk)
                   do m = 1, nrhfs(isymm)
                     idxkl=imatij(isymk,isyml)+nrhf(isymk)*(l-1)+k
                     idxlk=imatij(isyml,isymk)+nrhf(isyml)*(k-1)+l
                     idxmkl =imaijm(isymkl,isymm)+
     &                       nrhfs(isymm)*(idxkl-1)+m
                     idxmlk =imaijm(isymkl,isymm)+
     &                       nrhfs(isymm)*(idxlk-1)+m
                     work(krdpmlk+idxmlk-1) = work(krtf-1+idxmkl)       
                   end do 
                 end do
               end do 
             end do
           end do

           if (locdbg) then
             write(lupri,*) 'norm^2 of resorted R^dp_m,lk:',idel,
     &        ddot(nmaijm(isymd),work(krdpmlk),1,work(krdpmlk),1)
             write(lupri,*)'R^dp_m,lk', 
     &                      (work(krdpmlk-1+i),i=1,nmaijm(isymd))
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
             koffv = 1 + ivajkl(isymgi,isymkl) + ir1bas(isymg,isymi)
             
             ntotb  = max(1,mbas1(isymb))
             ntotgi = max(1,nr1bas(isymgi))
          
             call dgemm('T','N',mbas1(isymg)*nrhf(isymi),nmatij(isymkl),
     &            mbas1(isymb),-one,work(koffg),ntotb,
     &            work(koffr),ntotb,one,vajkl(koffv),ntotgi)
             if (locdbg) then 
               write(lupri,*)'DEBUG Vajkl(koffv) auxbas'
               call output(vajkl(koffv),1,mbas1(isymg)*nrhf(isymi),1,
     &           nmatij(isymkl),nr1bas(isymgi),nmatij(isymkl),
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
             koffr = krdpmlk + imaijm(isymkl,isymm)
             koffv = 1 + ivajkl(isymgi,isymkl) + ir1bas(isymg,isymi)

             ntotm  = max(nrhfs(isymm),1)
             ntotgi = max(nr1bas(isymgi),1)
             
             call dgemm('T','N',mbas1(isymg)*nrhf(isymi),nmatij(isymkl),
     &           nrhfs(isymm),-one,work(koffg),ntotm,work(koffr),ntotm,
     &           one,vajkl(koffv),ntotgi)

c
             xnrmr = xnrmr + ddot(nrhfs(isymm)*nmatij(isymkl),
     &                           work(koffr),1,work(koffr),1)
             xnrmg = xnrmg + ddot(nrhfs(isymm)*mbas1(isymg)*nrhf(isymi),
     &                           work(koffg),1,work(koffg),1)
c
             if (locdbg) then 
               write(lupri,*)'Vajkl(koffv) ansatz 2 auxbas'
               call output(vajkl(koffv),1,mbas1(isymg)*nrhf(isymi),1,
     &                    nmatij(isymkl),nr1bas(isymgi),
     &                    nmatij(isymkl),1,lupri)
             end if 
           end do 
        end if
      end if

      if (locdbg) then
        write(lupri,*)'Vajkl leaving mkvamkl'
        write(lupri,*)'norm^2:',ddot(nvajkl(1),vajkl,1,vajkl,1)
c       write(lupri,*)(vajkl(i),i=1,nvajkl(1))
        write(lupri,*)'xnrmr, xnrmg:',xnrmr,xnrmg
        write(lupri,*)'x1, x2, x3:',x1, x2, x3
      end if

      call qexit('r12mkvamkl')
      end                       
*====================================================================*
      subroutine cc_r12mkvamkl0(vajkl,vabkl,cc2r12,cmo,
     &                          work,lwork) 
c--------------------------------------------------------------------
c     purpose: initialize V_{alpha j, kl}
c  
c     H. Fliegl, C. Haettig spring 2003
c     modified for CCSD(R12) model, summer 2004
c--------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
 
      logical cc2r12,locdbg
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
      integer igabj(8)
      integer isymd,isymj,isymgk,isymg,isymk,isymkl
      integer isymab,isyma,isymb,isymaj,isymdl,isyml,isym0
      integer idxab,idxak,idxkl,idxlk,idxajkl,koffc,ntotb,ntota
      integer koverlp,ksab,ksak,kend1,lwrk1,lwork,koff1,koffs
      integer kvajkl,idxaj,nnbast,nvabkl(8),ivabkl(8,8)
      integer idxabkl,idxbl
#if defined (SYS_CRAY)
      real vajkl(*),cmo(*),work(*),vabkl(*)
#else
      double precision vajkl(*),cmo(*),work(*),vabkl(*)
#endif
      
      call qenter('r12mkvamkl0')
     
      call cc_r12offset(nr1orb,nr1xorb,nr1bas,nr1xbas,nr2bas,
     &     nrgkl,nrxgkl,n2bst1,ir1orb,ir1xorb,ir1bas,ir1xbas,ir2bas,
     &     ir2xbas,irgkl,irxgkl,iaodis1,nvabkl,ivabkl)

      nnbast = 0
      do i = 1,nsym
         nnbast = nnbast + nbas(i)*(nbas(i)+1)/2
      end do

      koverlp = 1
      ksab = koverlp + n2bst(1)
      ksak = ksab + n2bst1(1)
      kend1 = ksak + nr1bas(1)  
      lwrk1 = lwork - kend1
      if (lwrk1.lt.nnbast) then
         write(lupri,*)'lwrk1,lwork',lwrk1,lwork
         call quit('Insufficient work space in cc_r12mkvamkl0.')
      end if

c    --------------------------------------------------------------
c     read overlap matrix from file and unpack as quadratic matrix
c    --------------------------------------------------------------
      call rdonel('OVERLAP ',.true.,work(kend1),nnbast)

      isym0 = 1
      call ccsd_symsq(work(kend1),isym0,work(koverlp))

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
c        transform S_{alpha beta} with C_{beta k} 
c        to obtain S_{alpha k}
c       -------------------------------------------            
         koffc = 1 + iglmrh(isymb,isymk)
         koffs = ksak + ir1bas(isyma,isymk)
         ntotb = max(nbas(isymb),1)
         ntota = max(mbas1(isyma),1)
         
         call dgemm('N','N',mbas1(isyma),
     &        nrhf(isymk),mbas1(isymb),one,work(ksab),ntota,
     &        cmo(koffc),ntotb,zero,work(koffs),ntota)
      end do
  
      if (locdbg) then
        write(lupri,*)'Sak ='
        do isyma = 1, nsym
          isymk = isyma
          call output(work(ksak),1,mbas1(isyma),1,nrhf(isymk),
     &               mbas1(isyma),nrhf(isymk),1,lupri)
        end do
      end if
      
      if (cc2r12) then
c       -------------------------------------------------
c       initialize V(aj,kl) with S(alpha,k)*dkron(j,l)
c       -------------------------------------------------
        call dzero(vajkl,nvajkl(1))
        
        do isymkl=1, nsym
           isymaj = isymkl
           do isymk =1, nsym
              isyml = muld2h(isymkl,isymk)
              do k=1, nrhf(isymk)
                 do l=1, nrhf(isyml)
                    idxkl=imatij(isymk,isyml)+nrhf(isymk)*(l-1)+k
                    isymj = isyml
                    isyma = muld2h(isymaj,isymj)   
                    j=l
                    do a=1, mbas1(isyma)
                       idxaj = ir1bas(isyma,isymj)+mbas1(isyma)*(j-1)+a
                       idxak = ir1bas(isyma,isymk)+mbas1(isyma)*(k-1)+a
                       idxajkl = ivajkl(isymaj,isymkl)+nr1bas(isymaj)*
     &                      (idxkl -1) + idxaj 
                       vajkl(idxajkl) = work(ksak-1+idxak)
                    end do
                 end do
              end do
           end do
        end do
c
      else
c     ----------------------------------------------
c     initialize V(ak,bl) with S(alpha,k)*S(beta,l)
c     ----------------------------------------------
        call dzero(vabkl,nvabkl(1))

        do isymkl = 1, nsym
          isymab = isymkl
          do isymk = 1, nsym
            isyml = muld2h(isymkl,isymk) 
            isyma = isymk
            isymb = isyml
            do k = 1, nrhf(isymk)
              do l = 1, nrhf(isyml)
                idxkl = imatij(isymk,isyml)+nrhf(isymk)*(l-1)+k
                do a = 1, mbas1(isyma)
                  idxak = ir1bas(isyma,isymk)+mbas1(isyma)*(k-1)+a
                  do b = 1, mbas1(isymb)
                    idxab = iaodis(isyma,isymb)+mbas1(isyma)*(b-1)+a
                    idxbl = ir1bas(isymb,isyml)+mbas1(isymb)*(l-1)+b
                    idxabkl = ivabkl(isymab,isymkl)+
     &                        n2bst(isymab)*(idxkl-1)+idxab
                    vabkl(idxabkl) = 
     &                        work(ksak-1+idxak)*work(ksak-1+idxbl)
                  end do
                end do
              end do
            end do
          end do
        end do 
 
      end if 
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
      integer nr2bas(8)
      integer ir1orb(8,8),ir1bas(8,8),ir2bas(8,8),ir1xbas(8,8)
      integer ir2xbas(8,8),iaodis1(8,8),irgkl(8,8)
      integer nr1xorb(8),ir1xorb(8,8),nrxgkl(8),irxgkl(8,8)
      integer igabj(8),isycmo,isymv
      integer isyma,isymk,isymi,isymj,isyml,isymaj,isymkl
      integer idxkl,nrhftria,isymij,idxijkl
      integer lwork,ntota,ntoti,ntota1
      integer kvajkl,koffc,kscr1,kvijkl,ksing,ktrip,kend1 
      integer ik,il,ii,ij,idxij,iklij,lwrk1,nvabkl(8),ivabkl(8,8)
#if defined (SYS_CRAY)
      real vajkl(*),cmo(*),work(*),vijkl(*),ddot,fac
#else
      double precision vajkl(*),cmo(*),work(*),vijkl(*),ddot,fac
#endif

      call qenter('r12mkvijkl')
      
      call cc_r12offset(nr1orb,nr1xorb,nr1bas,nr1xbas,nr2bas,
     &     nrgkl,nrxgkl,n2bst1,ir1orb,ir1xorb,ir1bas,ir1xbas,ir2bas,
     &     ir2xbas,irgkl,irxgkl,iaodis1,nvabkl,ivabkl)

      if (locdbg) then
        write(lupri,*) 'in cc_r12mkvijkl: vajkl',
     &     ddot(nvajkl(isymv),vajkl,1,vajkl,1)
        do isymkl = 1, nsym 
          isymaj  = muld2h(isymv,isymkl)
          write(lupri,*) 'symmetry block ',isymaj,isymkl 
          call output(vajkl(1+ivajkl(isymaj,isymkl)),1,nr1bas(isymaj),
     &        1,nmatij(isymkl),nr1bas(isymaj),nmatij(isymkl),1,lupri)
        end do
      end if

      nrhftria = nrhft*(nrhft+1)/2

      kscr1  = 1
      kend1 = kscr1  + nrhft*nrhft
      if (.not. lres) then
        kvijkl = kend1 
        ksing  = kvijkl + nrhft*nrhft*nrhft*nrhft
        ktrip  = ksing  + nrhftria * nrhftria
        kend1  = ktrip  + nrhftria * nrhftria
      end if

      lwrk1  = lwork -kend1
      if (lwrk1 .lt.0) then
         call quit('Insufficient work space in r12mkvijkl')
      end if

      if (.not. lres) then
        call dzero(work(kvijkl),nrhft*nrhft*nrhft*nrhft)
      end if 

      do isymkl = 1, nsym
         isymaj = muld2h(isymv,isymkl)
         do isyml =1, nsym
            isymk = muld2h(isymkl,isyml)
            do k = 1, nrhf(isymk)
               ik = irhf(isymk) + k
               do l = 1, nrhf(isyml)
                  il = irhf(isyml) + l
                  idxkl=imatij(isymk,isyml)+nrhf(isymk)*(l-1)+k
                  do isyma =1, nsym
                     isymj = muld2h(isymaj,isyma)
                     isymi = muld2h(isycmo,isyma)

                     koffc = iglmrh(isyma,isymi)+1
                     kvajkl = ivajkl(isymaj,isymkl)+nr1bas(isymaj)*
     &                       (idxkl -1) +ir1bas(isyma,isymj)+1
                     ntota1 = max(1,mbas1(isyma))
                     ntota = max(1,nbas(isyma))
                     ntoti = max(1,nrhf(isymi))

                     call dgemm('T','N',nrhf(isymi),
     &                    nrhf(isymj),mbas1(isyma),one,cmo(koffc),
     &                    ntota,vajkl(kvajkl),ntota1,zero,work(kscr1),
     &                    ntoti)
c                    ------------------------------------------
c                     case lres (response) pack V(it,jt,kl) as
c                     in vector function
c                    ------------------------------------------
                     if (lres) then
                       isymij  = muld2h(isymi,isymj)
                       idxijkl = igamsq(isymij,isymkl)+ 
     &                 nmatij(isymij)*(idxkl-1)+imatij(isymi,isymj)+1
                       call daxpy(nrhf(isymi)*nrhf(isymj),fac,
     &                            work(kscr1),1,vijkl(idxijkl),1)
                     else
                       do j = 1, nrhf(isymj)
                         ij = irhf(isymj) + j
                         do i = 1, nrhf(isymi)
                           ii = irhf(isymi) + i
                           idxij = nrhf(isymi)*(j-1) + i
                           iklij=nrhft*(nrhft*(nrhft*(il-1)+
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

      if (.not.lres) then
        if (locdbg) then
          write(lupri,*) ' norm^2 Vijkl:',
     &     ddot(nrhft*nrhft*nrhft*nrhft,work(kvijkl),1,work(kvijkl),1)
          call output(work(kvijkl),1,nrhft*nrhft,1,nrhft*nrhft,
     &                   nrhft*nrhft,  nrhft*nrhft,1,lupri)
          call cc_r12mkvijkl2(work(kvijkl),work(ksing),
     &                      work(ktrip),nrhftria)
          write(lupri,*)'norm^2 Vsing:',
     &           ddot(nrhftria*nrhftria,work(ksing),1,work(ksing),1)
          write(lupri,*)'norm^2 Vtrip:',
     &           ddot(nrhftria*nrhftria,work(ktrip),1,work(ktrip),1)
        end if
      end if
      if (lres .and. locdbg) then
        write(lupri,*) ' Vijkl:',ddot(ngamsq(muld2h(isymv,isycmo)),
     &                                vijkl,1,vijkl,1)
        do isymkl = 1, nsym
           isymij = muld2h(muld2h(isymv,isycmo),isymkl)
           write(lupri,*) ' symmetry block:',isymij,isymkl
           call output(vijkl(1+igamsq(isymij,isymkl)),
     &                 1,nmatij(isymij),1,nmatij(isymkl),
     &                   nmatij(isymij),nmatij(isymkl),1,lupri)
        end do
      end if


      call qexit('r12mkvijkl') 
      end 
*=====================================================================*
      subroutine cc_r12mkvijkl2(vijkl,vsing,vtrip,nrhftria)
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

      integer ij,kl,nrhftria
#if defined (SYS_CRAY)
      real vijkl(nrhft,nrhft,nrhft,nrhft),ff
      real vsing(nrhftria,nrhftria),vtrip(nrhftria,nrhftria)
#else
      double precision vijkl(nrhft,nrhft,nrhft,nrhft),ff
      double precision vsing(nrhftria,nrhftria)
      double precision vtrip(nrhftria,nrhftria)
#endif

      ff = 1.0d0/sqrt(2.0d0)

      ij = 0
      do i=1,nrhft
         do j=1,i
            ij = ij + 1
            kl = 0
            do k=1,nrhft
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
     &     ibasx,lauxd,lauxbeta,fconst,work,lwork)
c---------------------------------------------------------------------
c     purpose: get R12 integrals 
c
c     H. Fliegl, C. Haettig spring 2003
c     
c     added possilbility to read r12_(beta,kl) with beta being 
c     auxiliary function:  C. Neiss  winter 2004/2005
c
c     nr1bas     dimension (mbas1 OR mbas2) x nrhf
c     ir1bas     symmetry offset of mbas1 x nrhf
c     ir1xbas    symmetry offset of mbas2 x nrhf
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

      lur2back = 0
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
         do l = 1, nrhf(isyml)
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
             write(lupri,*)'ndell,iadr,l,idelta',ndell,iadr,
     &              l,idelta
           end if
           ilen  = nr1bas(isymbk)
           call getwa2(lur2back,fconst,work(kscr2),iadr,ilen)
c 
           do isymb =1, nsym
             isymk = muld2h(isymbk,isymb)
             isymkl = muld2h(isymk,isyml)
             do k=1, nrhf(isymk)
               if (.not. lauxbeta) then 
                 do b=1, mbas1(isymb)
                   idxbk = ir1bas(isymb,isymk) + mbas1(isymb)*(k-1)+b
                   idxkl = imatij(isymk,isyml) + nrhf(isymk)*(l-1)+k
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
                   idxkl = imatij(isymk,isyml) + nrhf(isymk)*(l-1)+k
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
     &                       nmatij(isymkl),mbas1(isymb),nmatij(isymkl),
     &                       1,lupri)
               else 
                 call output(r12bkl(idxbkl),1,mbas2(isymb),1,
     &                       nmatij(isymkl),mbas2(isymb),nmatij(isymkl),
     &                       1,lupri) 
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
     &     nr2bas,ir2bas,nrgkl,irgkl,ir1xbas,ir2xbas,ibasx,lauxd,
     &     fconst,work,lwork)
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

      lur2back = 0
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
         do l = 1, nrhf(isyml)

            do isymb =1, nsym
               isymk = muld2h(isymbk,isymb)
               isymkl = muld2h(isymk,isyml)
               do k=1, nrhf(isymk)
                  do b=1, mbas1(isymb)
                     idxbk = ir1bas(isymb,isymk) + mbas1(isymb)*(k-1)+b
                     idxkl = imatij(isymk,isyml) + nrhf(isymk)*(l-1)+k
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
#if defined (SYS_CRAY)
      real zero,one
      real vijkl(*)
#else
      double precision zero,one
      double precision vijkl(*)
#endif
      parameter (zero = 0.0d0, one = 1.0d0)

      integer isyml,isymk,isymkl,idxkl,idxklkl
      
      call dzero(vijkl,ngamsq(1))
      
      do isymk = 1, nsym
        do isyml = 1, nsym
         isymkl = muld2h(isyml,isymk)
         do k =1, nrhf(isymk)
            do l = 1, nrhf(isyml)
               idxkl  = imatij(isymk,isyml) + nrhf(isymk)*(l-1)+k
               idxklkl=igamsq(isymkl,isymkl)+
     &              nmatij(isymkl)*(idxkl-1)+idxkl
               vijkl(idxklkl) = one
            end do
         end do
       end do
      end do
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
            nvajkl(isymak) = nvajkl(isymak) + nt1ao(isyma)*nmatij(isymk)
            icount7 = icount7 + nt1ao(isyma)*nmatij(isymk) 
         end do
      end do

      kctil = 1
      kvajkl = kctil + ngamsq(isymc)
      kend1 = kvajkl + nvajkl(1)
      lwrk1 = lwork - kend1
      if (lwrk1 .lt.0) then
         call quit('Insufficient work space in ccrhs_gp')
      end if

c     --------------------
c      read V(alpha j,kl)
c     --------------------
      luvajkl = 0
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
      call cc_r12getct(work(kctil),isymc,iamp,ketscl,.true.,.true.,
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
               
      do isymak = 1, nsym
         nvajkl(isymak) = 0
         icount7 = 0
         do isymk = 1, nsym
            isyma = muld2h(isymak,isymk)
            ivajkl(isyma,isymk) = icount7
            nvajkl(isymak) = nvajkl(isymak) + nt1ao(isyma)*nmatij(isymk)
            icount7 = icount7 + nt1ao(isyma)*nmatij(isymk) 
         end do
      end do

      kctil = 1
      kvajkl = kctil + ngamsq(isymc)
      kend1 = kvajkl + nvajkl(1)
      lwrk1 = lwork - kend1

      if (lwrk1 .lt.0) then
         call quit('Insufficient work space in ccrhs_gp')
      end if

c     --------------------
c      read V(alpha j,kl)
c     --------------------
        luvajkl = 0
        call gpopen(luvajkl,fvajkl,'old',' ','unformatted',
     &       idummy,.false.)
        rewind(luvajkl)
        read(luvajkl) (work(kvajkl+i-1), i = 1,nvajkl(1))
        call gpclose(luvajkl,'KEEP')

C     ---------------------------------------------------
C     read R12 amplitudes and form (2*c_ij^kl - c_ji^kl)
C     ---------------------------------------------------
       call cc_r12getct(work(kctil),isymc,iamp,ketscl,.true.,.true.,
     &                  lufc12,fc12am,ifile,listr,idlstr,
     &                  work(kend1),lwrk1)
C
      if (locdbg) then
        write(lupri,*) 'Norm^2(R12 ampl.)=', ddot(ngamsq(isymc),
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
CCN            do idxkl = 1, nmatij(isymkl)
CCN         if (locdbg) vnrm = zero
CCN         if (locdbg) vcprod = zero
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
            do idxkl = 1, nmatij(isymkl)                      
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
                 koffc = kctil + igamsq(isymnj,isymkl)+
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
               koffc = kctil + igamsq(isymij,isymkl)+
     &              nmatij(isymij)*(idxkl-1) + imatij(isymi,isymj)

               call dgemm('N','T',nrhf(isymn),nrhf(isymi),nrhf(isymj),
     &                one,work(kvint),ntotn,work(koffc),ntoti,
     &                one,eim(koffo),ntotn)

            end do
            if (locdbg) then
             koffc = kctil + igamsq(isymij,isymkl)+
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
      subroutine cc_r12getct(ctilde,isymc,iamp,basscl2,do_2cme,ltrans,
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
c     ltrans:   flag wheter lower (r12) indices (kl) are leading
C                  (-> .FALSE.), or upper (occ.) indices (ij) are 
C                  leading (-> .TRUE.) in ctilde
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

      integer lwork, kend1, kpck, lwrk1, lu43, idum, isymkl, isymij, 
     &        isyml, isymk, idxkl, isymi, isymj, idxij, idxji, idxklij,
     &        idxklji, ifile, lufc12, isymc, iopt
      integer iamp,idlstr,lunit
#if defined (SYS_CRAY) 
      real work(*), ctilde(*), two, ddot, basscl2
#else
      double precision work(*), ctilde(*), two, ddot, basscl2
#endif
      parameter (two = 2.0d0)
      logical do_2cme, ldum, locdbg, ltrans
      parameter (locdbg = .false.)
      character*(*) fc12am
      character     listr*(3)
      character     model*(10)

      CALL QENTER('cc_r12getct')

      kpck  = 1
      kend1 = kpck + ngamma(isymc)
      lwrk1 = lwork - kend1
      if (lwrk1.lt.0) then
         write(lupri,*) 'lwork, lwrk1: ',lwork, lwrk1
         call quit('Insufficient work space in cc_r12getct')
      end if

      if (iamp .EQ. 1) then
c       --------------------
c       read R12 trial vector 
c       -------------------
        call cc_rvec(lufc12,fc12am,ngamma(isymc),ngamma(isymc),
     *               ifile,work(kpck))
        call cclr_diasclr12(work(kpck),basscl2,isymc)
        call ccr12unpck2(work(kpck),isymc,ctilde,ltrans)
      else if (iamp .EQ. 2) then
C       -------------------------
C       read response amplitudes 
C       -------------------------
        iopt = 32
        call cc_rdrsp(listr,idlstr,isymc,iopt,model,dummy,work(kpck))
        call cclr_diasclr12(work(kpck),basscl2,isymc)
        call ccr12unpck2(work(kpck),isymc,ctilde,ltrans)
      else if (iamp .eq. 0) then
c       --------------------------------
c       read R12 ground state amplitudes 
c       --------------------------------
        if (isymc.ne.1) call quit('Symmetry error in CC_R12GETCT')
        iopt = 32
        call cc_rdrsp('R0 ',0,1,iopt,model,dummy,work(kpck))
        call ccr12unpck2(work(kpck),isymc,ctilde,ltrans)
      else
        call quit('Illegal value for iamp in cc_r12getct!')
      end if

      if (do_2cme) then
        call cc_r12tcmesq(ctilde,isymc,.false.)
      end if

      if (locdbg) then
        write(lupri,*) 'resorted R12 amplitudes (ctilde):'
        write(lupri,*) 'norm^2:',ddot(ngamsq(isymc),
     &         ctilde,1,ctilde,1)
        do isymkl = 1, nsym
          isymij = muld2h(isymkl,isymc)
          write(lupri,*) 'symmetry block ',isymkl,isymij
          call output(ctilde(1+igamsq(isymkl,isymij)),
     &       1,nmatij(isymkl),1,nmatij(isymij),
     &       nmatij(isymkl),nmatij(isymij),1,lupri)
        end do
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
                         
            do idxkl = 1, nmatij(isymkl)
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
               end if
c              ------------------------------
c                make ctilde*V(atilde,m,k,l) 
c              ------------------------------

               ntoti = max(1,nrhf(isymi)) 
               koffo = 1+it1am(isyma,isymi)
               koffc = 1 + igamsq(isymij,isymkl)+
     &              nmatij(isymij)*(idxkl-1) + imatij(isymi,isymj)
                 call dgemm('N','T',nvir(isyma),nrhf(isymi),nrhf(isymj),
     &                one,work(kvint),ntota,ctilde(koffc),ntoti,
     &                one,omega1(koffo),ntota)

c                if (locdbg) then
c                  write(lupri,*)'in ccrhs_gp0: ctilde',idxkl
c                  write(lupri,*)'isymij,isymkl,isymi,isymj:',
c    &                        isymij,isymkl,isymi,isymj
c                  call output(ctilde(koffc),1,nrhf(isymi),1,
c    &                    nrhf(isymj),nrhf(isymi),nrhf(isymj),1,lupri)
c                  write(lupri,*)'in ccrhs_gp0: omega1 + V*ctilde'
c                  call output(omega1(koffo),1,nvir(isyma),1,
c    &                nrhf(isymi),nvir(isyma),
c    &                nrhf(isymi),1,lupri)
c                end if

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
      subroutine cc_r12tcmesq(MATRIX,ISYM,LINV)
c---------------------------------------------------------------------
c     purpose: calculate ctilde = (2*c^ij_kl - c^ji_kl)
c
c     MATRIX   input and output, dimension ngamsq(ISYM)
c     LINV     do the inverse operation        
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
        call cc_prsqr12(matrix,isym,1,.false.)
      end if

      fac = -one  
      if (LINV) fac = one

      do isymkl = 1, nsym
        isymij = muld2h(isymkl,isym)
        do isymk = 1, nsym
          isyml = muld2h(isymk,isymkl)
          do k = 1, nrhf(isymk)
            do l = 1, nrhf(isyml)
              idxkl = imatij(isymk,isyml)+nrhf(isymk)*(l-1)+k
              do isymi = 1, nsym
                isymj = muld2h(isymi,isymij)
                do i = 1, nrhf(isymi)
                  do j = 1, nrhf(isymj)
                    idxij = imatij(isymi,isymj)+nrhf(isymi)*(j-1)+i
                    idxji = imatij(isymj,isymi)+nrhf(isymj)*(i-1)+j
                    idxklij = igamsq(isymkl,isymij)+
     &                        nmatij(isymkl)*(idxij-1)+idxkl
                    idxklji = igamsq(isymkl,isymij)+
     &                        nmatij(isymkl)*(idxji-1)+idxkl
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
        call dscal(ngamsq(isym),1.0D0/3.0D0,matrix,1)
      end if
                    
      if (locdbg) then
        call around ('Result in CC_R12TCMESQ')
        call cc_prsqr12(matrix,isym,1,.false.)

        write(lupri,*) 'Leaving CC_R12TCMESQ'
      end if
      call qexit('cc_r12tcmesq')
      return
      end
c=====================================================================*

*=====================================================================*
      subroutine cc_r12tcmesq_old(MATRIX,ISYM,WORK,LWRK)
c---------------------------------------------------------------------
c     purpose: calculate ctilde = (2*c^ij_kl - c^ji_kl)
c
c     MATRIX   input and output, dimension ngamsq(ISYM)
c
c     Christian Neiss  Feb. 2005
c---------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <dummy.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>

      logical locdbg
      parameter(locdbg = .false.)

      integer isym
      integer isymij,isymkl,isymi,isymj,isymk,isyml
      integer idxij,idxkl,idxji,idxklij,idxklji
      integer lwrk,kscr,kend1,lwrk1

#if defined (SYS_CRAY)
      real matrix(*),two,xklij,xklji,work(*)
#else
      double precision matrix(*),two,xklij,xklji,work(*)
#endif
      parameter (two = 2.0d0)

      call qenter('cc_r12tcmesq')
      if (locdbg) then
        write(lupri,*) 'Entered CC_R12TCMESQ'

        call around('Input in CC_R12TCMESQ')
        call cc_prsqr12(matrix,isym,1,.false.)
      end if

      kscr = 1
      kend1 = kscr + ngamsq(isym)
      lwrk1 = lwrk - kend1
      if (lwrk1 .gt. lwrk) then
        write(lupri,*) 'LWRK, LWRK1: ', lwrk, lwrk1
        call quit('Insufficient memory in CC_R12TCMESQ')
      end if


C     ------------------------------------
C     Calculate 2V(ij,kl)-V(ji,kl)
C     ------------------------------------
      do isymkl = 1, nsym
         isymij = muld2h(isymkl,isym)
         do isyml =1, nsym
            isymk = muld2h(isymkl,isyml)
            do k = 1, nrhf(isymk)
               do l = 1, nrhf(isyml)
                  idxkl=imatij(isymk,isyml)+nrhf(isymk)*(l-1)+k
                  do isymi = 1, nsym
                     isymj = muld2h(isymij,isymi)
                     do j = 1, nrhf(isymj)
                        do i = 1, nrhf(isymi)
                           idxij = imatij(isymi,isymj)+
     &                          nrhf(isymi)*(j-1) + i
                           idxji = imatij(isymj,isymi)+
     &                          nrhf(isymj)*(i-1) + j
                           idxklij =  igamsq(isymkl,isymij)+
     &                          nmatij(isymkl)*(idxij-1)+idxkl
                           idxklji =igamsq(isymkl,isymij)+
     &                          nmatij(isymkl)*(idxji-1)+idxkl
                           work(kscr-1+idxklij) =
     &                          two*matrix(idxklij) - matrix(idxklji)
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do

      call dcopy(ngamsq(isym),work(kscr),1,matrix,1)

      if (locdbg) then
        call around ('Result in CC_R12TCMESQ')
        call cc_prsqr12(matrix,isym,1,.false.)

        write(lupri,*) 'Leaving CC_R12TCMESQ'
      end if

      call qexit('cc_r12tcmesq')
      return
      end
c=====================================================================*

*======================================================================
       subroutine cc_r12tcmepk(MATRIX,ISYM,LINV)
C----------------------------------------------------------------------
c     purpose: calculate ctilde = (2*c^ij_kl - c^ji_kl)
c
c     MATRIX   input and output, dimension ngammaISYM)
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

      if (locdbg) write(lupri,*) 'Entered CC_R12TCMEPK'

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
                     DO 150 L = 1,NRHF(ISYML)
C
                        IF (ISYML .EQ. ISYMK) THEN
                           NRHFK = L
                        ELSE
                           NRHFK = NRHF(ISYMK)
                        ENDIF
C
                        NLI = IMATIJ(ISYML,ISYMI)
     *                      + NRHF(ISYML)*(I - 1) + L
                        NLJ = IMATIJ(ISYML,ISYMJ)
     *                      + NRHF(ISYML)*(J - 1) + L
C
                        DO 160 K = 1,NRHFK
C
                           NKI = IMATIJ(ISYMK,ISYMI)
     *                         + NRHF(ISYMK)*(I - 1) + K
                           NKJ = IMATIJ(ISYMK,ISYMJ)
     *                         + NRHF(ISYMK)*(J - 1) + K
C
                         IF (ISYMKI.EQ.ISYMLJ) THEN
                           NKILJ = IGAMMA(ISYMKI,ISYMLJ)
     *                           + INDEX(NKI,NLJ)
                         ELSE IF (ISYMKI.LT.ISYMLJ) THEN
                           NKILJ = IGAMMA(ISYMKI,ISYMLJ)
     *                           + NMATIJ(ISYMKI)*(NLJ-1)+NKI
                         ELSE IF (ISYMLJ.LT.ISYMKI) THEN
                           NKILJ = IGAMMA(ISYMKI,ISYMLJ)
     *                           + NMATIJ(ISYMLJ)*(NKI-1)+NLJ
                         END IF
C
                         IF (ISYMKJ.EQ.ISYMLI) THEN
                           NKJLI = IGAMMA(ISYMKJ,ISYMLI)
     *                           + INDEX(NKJ,NLI)
                         ELSE IF (ISYMKJ.LT.ISYMLI) THEN
                           NKJLI = IGAMMA(ISYMKJ,ISYMLI)
     *                           + NMATIJ(ISYMKJ)*(NLI-1)+NKJ
                         ELSE IF (ISYMLI.LT.ISYMKJ) THEN
                           NKJLI = IGAMMA(ISYMKJ,ISYMLI)
     *                           + NMATIJ(ISYMLI)*(NKJ-1)+NLI
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
        call dscal(ngamma(isym),1.0D0/3.0D0,matrix,1)
      end if
C    
      if (locdbg) write(lupri,*) 'Leaving CC_R12TCMEPK'
      call qexit('cc_r12tcmepk')
      return

      end 
*======================================================================

