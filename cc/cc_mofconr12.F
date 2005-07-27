*======================================================================*
      subroutine cc_mofconr12(xlamdh,isymh,xlamdhs,xlamdps,xlamdhcs,
     &                        isymhcs,vijkl,facterm23,vajkl,vabkl,
     &                        lvijkl,lvajkl,lvabkl,
     &                        ioptbas,timrdao,timfr12,timintr12,
     &                        iglmrhs,nglmds,imaijm,nmaijm,
     &                        work,lwork)
c-----------------------------------------------------------------------
c     purpose: get V_{kl}^{itilde jtilde} contributions 
c              for auxiliary basis if lvijkl=.TRUE.
c              get V_{kl}^{alpha jtilde} contributions
c              for auxiliary basis if lvajkl=.TRUE.
c              get V_{kl}^{atilde btilde} contributions
c              for auxiliary basis if lvabkl=.TRUE.
c
c     H. Fliegl, C. Haettig, W. Klopper spring 2003
c     modified for CCSD(R12), summer 2004
c     modified for CC2-Response C. Neiss, autumn 2004
c
c     ioptbas controls type of basis functions used for outer loop
c     implemented values:
c              1       loop over auxiliary basis only
c              2       loop over orbital basis AND auxiliary basis
c-----------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <dummy.h>     
#include <maxorb.h>
CCN#include <maxash.h>
#include <mxcent.h>
#include <aovec.h>
#include <iratdef.h>
#include <ccorb.h>
#include <ccisao.h>
#include <blocks.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <cbieri.h>
#include <distcl.h>
#include <eribuf.h>
#include <r12int.h>
#include <ccr12int.h>
#include <second.h>
      logical locdbg,lvijkl,lvabkl,lvajkl,temp_direct
      parameter(locdbg =.false.)
#if defined (SYS_CRAY)
      real zero,one,two
#else
      double precision zero,one,two
#endif
      parameter (zero = 0.0d0, one = 1.0d0, two = 2.0d0)

      integer indexa(mxcorb)
      integer isymhcs,iglmrhs(8,8),nglmds(8),imaijm(8,8),nmaijm(8),
     &        lwork,index,kend1,lwrk1,ioff,isym,isymgk,icount2,
     &        icount5,isymk,isymg
      integer kccfb1,kindxb,kfree,lfree,ntosym,
     &        kendsv,lwrksv,isymd1,ntot,illl,numdis,idel2,idel,
     &        kodcl1,kodcl2,kodbc1,kodbc2,krdbc1,krdbc2,kscr1,kscr5,
     &        kodpp1,kodpp2,krdpp1,krdpp2,krecnr,kend2,lwrk2,kscr2
     &        isymd,isydis,kxint,isymd,kgaijd,isymh,
     &        ibasx(8),irgkl(8,8),nrgkl(8),ir1bas(8,8)
      integer ibastyp,ibastypst,ibastypend,idum,ioptbas

#if defined (SYS_CRAY)
      real xlamdh(*),work(*),vijkl(*),vajkl(*)
      real dtime,timrdao,timfr12,timintr12,xlamdhs(*)
      real vabkl(*),xlamdhcs(*),facterm23
#else
      double precision xlamdh(*),work(*),vijkl(*),vajkl(*)
      double precision dtime,timrdao,timfr12,timintr12,ddot,
     &                 xlamdhs(*),xlamdps(*)
      double precision vabkl(*),xlamdhcs(*),facterm23
#endif
      
      character*8 filback
      index(i,j) = max(i,j)*(max(i,j) - 3)/2 + i + j

      call qenter('mofconr12')

      kend1 = 1     
      lwrk1 = lwork - kend1
      if (lwrk1 .lt.0) then
         call quit('Insufficient work space in CC_MOFCONR12')
      end if

c     loop over orbital and auxiliary basis (e.g. for CCSD(R12) model)
c     or loop over auxiliary basis (e.g. for CC2-R12 model)
      if (ioptbas.eq.1) then
        ibastypst = 2
        ibastypend = ibastypst
      else if (ioptbas.eq.2) then
        ibastypst = 1
        ibastypend = 2
      else 
        call quit('Illegal value for ioptbas in CC_MOFCONR12!!')
      end if
c
      do ibastyp = ibastypst, ibastypend 
        if (ibastyp.eq.2) then 
          mbsmax = 5
          loopdp = .true.
          !these integrals are not on file yet -> switch locally to
          !direct mode and calculate them now
          TEMP_DIRECT = DIRECT
          DIRECT = .TRUE.
c       
          ioff    = 0
          ibasx(1) = 0
          do isym = 1,nsym
            if (isym.gt.1) ibasx(isym) = ibasx(isym-1)+mbas2(isym-1)
            do i = 1,mbas1(isym)+mbas2(isym)
              ioff = ioff + 1
              isao(ioff) = isym
            end do
          end do
  
        else if (ibastyp .eq. 1) then 
          mbsmax = 4
          loopdp = .false.

          ioff   = 0
          do isym = 1,nsym
            ibasx(isym) = 0
            do i = 1,nbas(isym)
              ioff = ioff + 1
              isao(ioff) = isym
            end do
          end do
        else
          call quit('Illegal value for "ioptbas" in CC_MOFCONR12!!')
        end if   

        do isymgk = 1, nsym
          nrgkl(isymgk) = 0
          icount2 = 0
          icount5 = 0
          do isymk = 1, nsym
            isymg = muld2h(isymgk,isymk)
            nrgkl(isymgk) = nrgkl(isymgk) + mbas1(isymg)*nmatij(isymk)
            ir1bas(isymg,isymk)  = icount2
            irgkl(isymg,isymk) = icount5
            icount2 = icount2 + nrhf(isymk)*mbas1(isymg)
            icount5 = icount5 + mbas1(isymg)*nmatij(isymk)
          end do
        end do
C  ====================================================
C       start the loop over distributions of integrals.
C  ====================================================
  
        if (locdbg) then
           write(lupri,'(1x,a,i10)') 'lwork = ',lwork
        end if
  
        if (direct) then
           dtime  = second()
           if (herdir) then
              call herdi1(work(kend1),lwrk1,ipreri)
           else
              kccfb1 = kend1
              kindxb = kccfb1 + mxprim*mxcont
              kend1  = kindxb + (8*mxshel*mxcont + 1)/irat
              lwrk1  = lwork  - kend1
              if (lwrk1 .lt.0) 
     &          call quit('Insufficient work space in CC_MOFCONR12')
              call eridi1(kodcl1,kodcl2,kodbc1,kodbc2,krdbc1,krdbc2,
     &                  kodpp1,kodpp2,krdpp1,krdpp2,
     &                  kfree,lfree,kend1,work(kccfb1),work(kindxb),
     &                  work(kend1),lwrk1,ipreri)
              kend1 = kfree
              lwrk1 = lfree
           endif
           timintr12 = timintr12 + ( second() - dtime )
           ntosym = 1
        else
           ntosym = nsym
        endif
  
        kendsv = kend1
        lwrksv = lwrk1
  
        do isymd1 = 1,ntosym
  
           if (direct) then
              if (herdir) then
                 ntot = maxshl
              else
                 ntot = mxcall
              endif
           else
              ntot = nbas(isymd1)
           endif
       
           do illl = 1,ntot
  
c  ---------------------------------------------
c             if direct calculate the integrals.
c  ---------------------------------------------
  
              if (direct) then
                 dtime = second()
  
                 kend1 = kendsv
                 lwrk1 = lwrksv
  
                 if (herdir) then
                    call herdi2(work(kend1),lwrk1,indexa,illl,numdis,
     &                        ipreri)
                 else
                    call eridi2(illl,indexa,numdis,0,0,
     &                        work(kodcl1),work(kodcl2),work(kodbc1),
     &                        work(kodbc2),work(krdbc1),work(krdbc2),
     &                        work(kodpp1),work(kodpp2),work(krdpp1),
     &                        work(krdpp2),work(kccfb1),work(kindxb),
     &                        work(kend1), lwrk1,ipreri)
                 endif
  
                 krecnr = kend1
                 kend1  = krecnr + (nbufx(0) - 1)/irat + 1
                 lwrk1  = lwork  - kend1
                 if (lwrk1 .lt.0) 
     &             call quit('Insufficient work space in CC_MOFCONR12')
                 timintr12 = timintr12 + ( second() - dtime )
              else
                 numdis = 1
              endif
  
c  -----------------------------------------------------
c             loop over number of distributions in disk.
c             loop delta' for R12 contributions
c  -----------------------------------------------------
  
              do idel2 = 1,numdis
  
                 if (direct) then
                    idel  = indexa(idel2)
                    if (noauxb.and.ibastyp.eq.1) then
                       idum = 1
                       call ijkaux(idel,idum,idum,idum)
                    end if
                    isymd = isao(idel)
                 else
                    idel  = ibas(isymd1) + ibasx(isymd1) + illl
                    isymd = isymd1
                 endif
                 
                 if (locdbg) then
                   write(lupri,*) 'in cc_mofconr12:'
                   write(lupri,*) 'idel,isymd:',idel,isymd
                   write(lupri,*) 'idel-ibas-ibasx,mbas1,mbas2:',
     &           idel-ibas(isymd)-ibasx(isymd),mbas1(isymd),mbas2(isymd)
                 end if
  
                 isydis = muld2h(isymd,isymop)
  
C               ------------------------------
C                Work space allocation no. 2.
C               ------------------------------
                 kxint  = kend1
                 if (lvajkl .or. lvijkl) then
                   kgaijd = kxint  + ndisao(isydis)
                   kend2  = kgaijd + nd2ijg(isydis)
                 else
                   kend2 = kxint + ndisao(isydis)
                 end if
                 lwrk2  = lwork  - kend2
c  
                 if (lwrk2 .lt. 0) then
                    write(lupri,*) 'need : ',kend2,'available : ',lwork
                    call quit('insufficient space in cc_mofconr12')
                 endif
  
c               -----------------------------
c                read in batch of integrals.
C               -----------------------------
                 dtime   = second()
                 call ccrdao(work(kxint),idel,idel2,work(kend2),lwrk2,
     &                     work(krecnr),direct)
                 dtime   = second() - dtime
                 timrdao = timrdao  + dtime
  
c                ---------------------------------------
c                  two index transformation for CC2-R12
c                ---------------------------------------
                 dtime = second()
                 if (lvajkl .or. lvijkl) then
                   call cc_r12mofcc2(xlamdh,isymh,xlamdhs,xlamdps,
     &                              isymhcs,work(kxint),vijkl,facterm23,
     &                               lvijkl,vajkl,lvajkl,work(kgaijd),
     &                               iglmrhs,nglmds,imaijm,nmaijm,
     &                               ibasx,nrgkl,irgkl,ir1bas,idel,
     &                               isymd,isydis,xlamdhcs,work(kend2),
     &                               lwrk2)

                 else if (lvabkl) then
c                  calculate ccsd(r12) part
                   call cc_r12mkvabkl(vabkl,work(kxint),idel,isymd,
     &                                isydis,ibastyp,ibasx,
     &                                work(kend2),lwrk2)
                 end if
                 timfr12 = timfr12 + ( second() - dtime )
              end do ! idel2
           end do ! illl
        end do ! isymd1
        if (ibastyp.eq.2) DIRECT = TEMP_DIRECT
      end do ! ibastyp

      mbsmax = 4
      loopdp = .false.

      ioff   = 0
      do isym = 1,nsym
        do i = 1,nbas(isym)
          ioff = ioff + 1
          isao(ioff) = isym
        end do
      end do

      call qexit('mofconr12')
      return      
      end
*=====================================================================*
      subroutine cc_r12mofcc2(xlamdh,isymh,xlamdhs,xlamdps,isymhcs,
     &                        xint,vijkl,facterm23,lvijkl,vajkl,lvajkl,
     &                        xgaijd,iglmrhs,nglmds,imaijm,nmaijm,
     &                        ibasx,nrgkl,irgkl,ir1bas,idel,
     &                        isymd,isydis,xlamdhcs,work,lwork)
c---------------------------------------------------------------------
c     purpose: calculate V^ij_kl (auxiliary basis) for CC2-R12 model
c
c     H. Fliegl, C. Haettig, summer 2004 
c---------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <maxorb.h>
#include <ccorb.h>
CCN#include <ccisao.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <r12int.h>
#include <ccr12int.h>

      logical locdbg,lauxd,lvajkl,lvijkl
      parameter(locdbg =.false.)

      character*8 filback

      integer iglmrhs(8,8),nglmds(8),imaijm(8,8),nmaijm(8),ibasx(8)
      integer nrgkl(8),irgkl(8,8),ir1bas(8,8),lwork,idel,isydis,
     &        isymg,isymj,isalbe,kscr1,kscr2,kend1,lwrk1,koff1,
     &        koff2,ntalbe,ntotg,isymhcs,isymd,isymi,isymbe,isymal,
     &        isyma,kscr5,kend2,lwrk2,koff3,ntotal,ntotbe,isymij,
     &        idxai,idxij,idxaij,krtf,isymm,ntota,ntotm,koffc,
     &        koffg,koffgtf,isymhs,isymh
#if defined (SYS_CRAY)
      real zero,one,two,xlamdh(*),xlamdhs(*),xlamdps(*),
     &     xint(*),vijkl(*),vajkl(*),work(*),xgaijd(*),ddot,xlamdhcs(*),
     &     facterm23
#else
      double precision zero,one,two,xlamdh(*),xlamdhs(*),
     &                 xlamdps(*),xint(*),vijkl(*),vajkl(*),work(*),
     &                 xgaijd(*),ddot,xlamdhcs(*),facterm23
#endif
      parameter (zero = 0.0d0, one = 1.0d0, two = 2.0d0)
 
      call qenter('mofcc2')
c
      if (locdbg) then
        write(lupri,*) 'Entered CC_R12MOFCC2'
        write(lupri,*) 'lwork = ',lwork
        call flshfo(lupri)
      end if
c
      do isymg = 1,nsym
       if (nbas(isymg) .gt. 0) then
         isalbe = muld2h(isymg,isydis)
         isymj  = muld2h(isymh,isymg)

c        ----------------------------------
c        dynamic allocation of work space.
c        -----------------------------------
         kscr1 = 1
         kscr2 = kscr1 + nnbst(isalbe)*nrhf(isymj)
         kend1 = kscr2 + n2bst(isalbe)
         lwrk1 = lwork - kend1
         
         if (lwrk1.lt.0) then
           call quit('insufficient work space in mofcc2')
         end if
         
c        ---------------------------------
c        do first index transformation
c        ---------------------------------             
         koff1 = 1 + idsaog(isymg,isydis) 
         koff2 = iglmrh(isymg,isymj) + 1
c         
         ntalbe = max(nnbst(isalbe),1)
         ntotg  = max(nbas(isymg),1)
c         
         call dgemm('N','N',nnbst(isalbe),nrhf(isymj),
     &     nbas(isymg),one,xint(koff1),ntalbe,
     &     xlamdh(koff2),ntotg,zero,work(kscr1),ntalbe)

         if (locdbg) then
c          write(lupri,*) 'in cc_r12mofcc2: Ints. before first '//
c    &                    'index transformation: '
c          call output(xint(koff1),1,nnbst(isalbe),1,nbas(isymg),
c    &                 nnbst(isalbe),nbas(isymg),1,lupri)
c          write(lupri,*) 'in cc_r12mofcc2: XLAMDH for first '//
c    &                    'index transformation: '
c          call output(xlamdh(koff2),1,nbas(isymg),1,nrhf(isymj),
c    &                 nbas(isymg),nrhf(isymj),1,lupri)
           write(lupri,*) 'in cc_r12mofcc2: Ints. after first '//
     &                    'index transformation: '
           call output(work(kscr1),1,nnbst(isalbe),1,nrhf(isymj),
     &                 nnbst(isalbe),nrhf(isymj),1,lupri)
         end if
 
c------------------------------------------------
c        compute contributions to V(alpha j,kl)
c------------------------------------------------
         if (lvajkl) then
          if (mbas1(isymg).gt.0 .or. nrhf(isymj).gt.0) then
            if (ianr12.eq.1) then
              filback = fnback
            else if (ianr12.eq.2 .or. ianr12.eq.3) then
              filback = frhtf
            end if
            call r12mkvamkl(filback,work(kscr1),work(kscr1),vajkl,
     &           xlamdh,isymh,xlamdhs,xlamdps,xint(koff1),xint(koff1),
     &           idel,isymd,isymj,
     &           isalbe,isymg,work(kscr2),ibasx,iglmrhs,
     &           nglmds,.true.,work(kend1),lwrk1)
          end if
         end if
c------------------------------------------------
c        compute contributions to V(ij,kl)
c------------------------------------------------
         if (lvijkl) then

c          -------------------------------
c          do last index transformations
c          -------------------------------             
           do j = 1,nrhf(isymj)
c          
             koff1 = kscr1 + nnbst(isalbe)*(j - 1)             
             call ccsd_symsq(work(koff1),isalbe,work(kscr2))
c          
             do isymi = 1,nsym
c          
               isymbe = isymi
               isymal = muld2h(isymbe,isalbe)
               isyma  = muld2h(isymal,isymhcs)
c          
               kscr5 = kend1 
               kend2 = kscr5 + nbas(isymal)*nrhf(isymi)
               lwrk2 = lwork - kend2 
               if (lwrk2 .lt. 0) then
                  call quit('insufficient space for 2. trf. '//
     &               'in cc_mofcc2')
               endif
c          
               koff2 = kscr2 + iaodis(isymal,isymbe)
               koff3 = ilmrhf(isymi) + 1
c          
               ntotal = max(nbas(isymal),1)
               ntotbe = max(nbas(isymbe),1)
c          
               call dgemm('N','N',nbas(isymal),nrhf(isymi),
     &                nbas(isymbe),one,work(koff2),ntotal,
     &                xlamdh(koff3),ntotbe,zero,work(kscr5),
     &                ntotal)
c          
c              --------------------------------------------
c              store (alpha i|j delta) as three index array
c              --------------------------------------------  
               isymij = muld2h(isymi,isymj)
               if (ianr12.eq.1) then
                 do i = 1, nrhf(isymi)
                    do a = 1, nbas(isymal)
                      idxai = nbas(isymal)*(i-1)+a
                      idxij = imatij(isymi,isymj)+
     &                        nrhf(isymi)*(j-1)+i
                      idxaij = id2ijg(isymij,isymal)+
     &                        nbas(isymal)*(idxij-1)+a 
                      xgaijd(idxaij) = 
     &                             work(kscr5-1+idxai) 
                    end do
                 end do
               else if (ianr12.eq.2 .or. ianr12.eq.3) then
                 do i = 1, nrhf(isymi)
                   do a = 1, nbas(isymal)
                     idxai = nbas(isymal)*(i-1)+a
                     idxij = imatij(isymi,isymj)+
     &                        nrhf(isymi)*(j-1)+i
                     idxaij = irgkl(isymal,isymij)+
     &                        nbas(isymal)*(idxij-1)+a 
                     xgaijd(idxaij) = 
     &                             work(kscr5-1+idxai) 
                   end do
                 end do
               end if
             end do ! isymi
           end do ! j
         end if ! (lvijkl)
       end if
      end do ! isymg

c     ----------------------
c     add terms for ansatz 2
c     ----------------------
      if (lvijkl) then
       if (ianr12.eq.2 .or. ianr12.eq.3) then
        isymhs = 1
        krtf  = kend1
        kend2 = krtf + nmaijm(isydis) 
        lwrk2 = lwork - kend2
        if (lwrk1.lt.0) then
          call quit('Insufficient work space in mofcc2')
        end if
      
c       CMO_am is equivalent with Lambda^p_am occupied block 
c       should be equivalent with Lambdas^h_am occupied block
        if (locdbg) then
          write(lupri,*)'I_ai,j in mofcc2',
     &           ddot(nd2ijg(isydis),xgaijd,1,xgaijd,1)
c         write(lupri,*)(xgaijd(i),i=1,nd2ijg(isydis))
        end if
c       transform to I_M,ij, M runs over all active and 
c       inactive molecular orbitals
        do isyma = 1, nsym
          isymij = muld2h(isyma,isymd)
c         isymm  = muld2h(isymhcs,isyma)
          isymm  = muld2h(isymhs,isyma)
c          
          ntota  = max(nbas(isyma),1)
          ntotm  = max(nrhfs(isymm),1)
c     
          koffc = 1+iglmrhs(isyma,isymm)
          koffg = 1+irgkl(isyma,isymij)
          koffgtf = krtf+imaijm(isymij,isymm)
c     
          call dgemm('T','N',nrhfs(isymm),nmatij(isymij),
     &                mbas1(isyma),one,xlamdps(koffc),ntota,
     &                xgaijd(koffg),ntota,zero,work(koffgtf),
     &                ntotm)
        end do
      
        if (locdbg) then
          write(lupri,*)'I_mi,j',idel,
     &      ddot(nmaijm(isydis),work(krtf),1,work(krtf),1)
c         write(lupri,*)(work(krtf+i-1),i=1,nmaijm(isydis))
        end if
c     
c       case respose calculation calculate second part for V
        lauxd = .true.
        call cc_r12mkvkl(work(krtf),vijkl,facterm23,xlamdhcs,iglmrhs,
     &                  isymd,isymhcs,idel,ibasx,imaijm,nmaijm,
     &                  lauxd,work(kend2),lwrk2,frhtf)
c
        if (locdbg) then
          write(lupri,*)'norm^2 Vijkl = ',
     &              ddot(ngamsq(isymhcs),vijkl,1,vijkl,1)
        end if
       else ! ansatz 1
        lauxd = .true.
        call cc_r12mkvkl(xgaijd,vijkl,facterm23,xlamdh,iglmrh,
     &                 isymd,isymhcs,idel,ibasx,imaijm,nmaijm,
     &                 lauxd,work(kend1),lwrk1,fnback)
       end if
      end if ! (lvijkl)
c
      if (locdbg) then
        write(lupri,*) 'Leaving CC_R12MOFCC2'
        call flshfo(lupri)
      end if
c
      call qexit('mofcc2')
      end
*=====================================================================*
      subroutine cc_r12mkvkl(gaijd,vijkl,facterm23,xlamdh,ioffc,isymd,
     &                       isymhcs,
     &                       idel,ibasx,imaijm,nmaijm,lauxd,
     &                       work,lwork,filback)
c---------------------------------------------------------------------
c     purpose: update V^{itilde jtilde}_{kl}
c 
c     H. Fliegl, C. Haettig
c---------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <maxorb.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
      logical locdbg,lauxd
      parameter(locdbg = .false.)
#if defined (SYS_CRAY)
      real zero,one,two
#else
      double precision zero,one,two
#endif
      parameter (zero = 0.0d0, one = 1.0d0, two = 2.0d0)

      integer lwork,isymd,idel,ir1xbas(8,8),irgkl(8,8),ir2xbas(8,8)
      integer kr2akl,kend1,lwrk1,icount8,isymgk,isymk,isymg,ntota1
      integer isymij,isyma,isymkl,koffr,koffg,koffv,ntota,ntotaj,
     &        ntotij,icount4,nr1xbas(8),nr1bas(8),nr2bas,nrgkl(8),
     &        ir1bas(8,8), ir2bas(8,8), ibasx(8),imaijm(8,8),nmaijm(8),
     &        icount2,icount3,icount5, idelta,ioffc(8,8)
      integer kend2,lwrk2,krtf,isymm,koff1,koff2,ntotm,isymhcs,isymr

#if defined (SYS_CRAY)
      real gaijd(*),vijkl(*),work(*),xlamdh(*)
      real factor, ddot,facterm23
#else
      double precision gaijd(*),vijkl(*),work(*),xlamdh(*)
      double precision factor, ddot,facterm23
#endif
      double precision xr
      save xr
      data xr /0.0d0/
      character*(*) filback

      call qenter('r12mkvkl')

      if (locdbg) then
        write(lupri,*)'g(aijd) in cc_r12mkvkl (beginning):'
        do isyma = 1, nsym
          isymij = muld2h(isymd,isyma)
          write(lupri,*) 'isymd,isyma,isymij:',isymd,isyma,isymij
          if (ianr12.eq.1) then
            call output(gaijd(1+id2ijg(isyma,isymij)),
     &                  1,nbas(isyma),1,nmatij(isymij),
     &                  nbas(isyma),nmatij(isymij),1,lupri)
          else if (ianr12.eq.2 .or. ianr12.eq.3) then
            isymm = isyma
            call output(gaijd(1+imaijk(isymij,isymm)),
     &                  1,nrhf(isymm),1,nmatij(isymij),
     &                  nrhf(isymm),nmatij(isymij),1,lupri)
          end if
        end do
      end if

      do isymgk = 1, nsym
        nr1bas(isymgk) = 0
        nr1xbas(isymgk) = 0
        do isymk = 1, nsym
           isymg = muld2h(isymgk,isymk)
           nr1bas(isymgk)  = nr1bas(isymgk) +mbas1(isymg)*nrhf(isymk)
           nr1xbas(isymgk) = nr1xbas(isymgk)+mbas2(isymg)*nrhf(isymk)
        end do
      end do

c     -----------------------------------
c     three index array for r12-integrals
c     -----------------------------------
      do isymgk = 1, nsym
        nrgkl(isymgk) = 0  
        do isymk = 1, nsym
          isymg = muld2h(isymgk,isymk)
          nrgkl(isymgk) = nrgkl(isymgk) + mbas1(isymg)*nmatij(isymk) 
         end do
      end do

c     -------------------------------------- 
c     lenghts for nr1bas over all symmetries
c     -------------------------------------- 
      nr2bas = 0
      do isymg = 1, nsym
         nr2bas = nr2bas + nr1bas(isymg)*nr1bas(isymg)
      end do

      do isymgk = 1, nsym
        icount2 = 0
        icount3 = 0
        icount4 = 0
        icount5 = 0
        icount8 = 0
        do isymk = 1, nsym
           isymg = muld2h(isymgk,isymk)
           ir1bas(isymg,isymk)  = icount2 
           ir2bas(isymg,isymk)  = icount3 
           ir2xbas(isymg,isymk) = icount4
           irgkl(isymg,isymk)   = icount5
           ir1xbas(isymg,isymk) = icount8
           icount2 = icount2 + nrhf(isymk)*mbas1(isymg)
           icount3 = icount3 + nr1bas(isymg)*nr1bas(isymk)
           icount4 = icount4 + nt1ao(isymg)*nr1xbas(isymk)
           icount5 = icount5 + mbas1(isymg)*nmatij(isymk) 
           icount8 = icount8 + nrhf(isymk)*mbas2(isymg)
        end do
      end do

      kr2akl = 1
      kend1 = kr2akl + nrgkl(isymd)
      lwrk1 = lwork-kend1

      if (lwrk1 .lt.0) then
         call quit('Insufficient work space in cc_r12mkvkl')
      end if
c    --------------
c     add ansatz 2
c    --------------
      call cc_r12getrint(work(kr2akl),idel,isymd,nr1bas,ir1bas,
     &               nr2bas,ir2bas,nrgkl,irgkl,ir1xbas,ir2xbas,
     &               ibasx,lauxd,.false.,filback,work(kend1),lwrk1)
      
      if (ianr12.eq.2 .or. ianr12.eq.3) then
        isymr = muld2h(isymhcs,isymd)
        krtf = kend1
        kend2 = krtf + nmaijm(isymr)
        lwrk2 = lwork - kend2
        if (lwrk2.lt.0) then
          call quit('Insufficient work space in cc_r12mkvkl')
        end if
c       
        if (locdbg) then
         write(lupri,*)'R_a,kl',
     &     ddot(nrgkl(isymd),work(kr2akl),1,work(kr2akl),1)
c        write(lupri,*)(work(kr2akl+i-1),i=1,nrgkl(isymd))
        end if

        do isyma = 1, nsym
          isymkl = muld2h(isyma,isymd)
          isymm  = muld2h(isyma,isymhcs)
          
          koff1 = 1+ioffc(isyma,isymm)
          koffr = kr2akl+irgkl(isyma,isymkl)
          koff2 = krtf+imaijm(isymkl,isymm)

          ntota = max(nbas(isyma),1)
          ntota1 = max(mbas1(isyma),1)
          ntotm = max(nrhfs(isymm),1)
 
          call dgemm('T','N',nrhfs(isymm),nmatij(isymkl),mbas1(isyma),
     &               one,xlamdh(koff1),ntota,work(koffr),ntota1,
     &               zero,work(koff2),ntotm) 
        end do 

        xr = xr + ddot(nmaijm(isymr),work(krtf),1,work(krtf),1)
        if (locdbg) then
          write(lupri,*)'R_m,kl in mkvkl', idel,
     &   ddot(nmaijm(isymr),work(krtf),1,work(krtf),1)
        end if

      end if
c
      idelta = idel - ibas(isymd)
      if (lauxd) idelta = idelta - ibasx(isymd)

      do isyma = 1, nsym
         
         if (idelta.le.mbas1(isymd)) then
            factor = one
         else 
c           factor = - two
c           if (.not.lmkvajkl) factor = -one
            factor = - facterm23
         end if
c
         if (ianr12.eq.2 .or. ianr12.eq.3) then
           isymm  = muld2h(isymhcs,isyma)
           isymij = muld2h(isymm,isymd)
           isymkl = muld2h(isymhcs,isymij)

           koffr = krtf + imaijm(isymkl,isymm) 
           koffg = 1 + imaijm(isymij,isymm)
           koffv = 1 + igamsq(isymij,isymkl)
         
           ntotm  = max(1,nrhfs(isymm))
           ntotij = max(1,nmatij(isymij))

           call dgemm('T','N',nmatij(isymij),nmatij(isymkl),
     &          nrhfs(isymm),factor,gaijd(koffg),ntotm,
     &          work(koffr),ntotm,one,vijkl(koffv),ntotij)
         else
           isymij = muld2h(isymd,isyma)
           isymkl = muld2h(isymd,isyma)
           isymm  = muld2h(isymhcs,isyma)

           koffr = kr2akl + irgkl(isyma,isymkl)
           koffg = 1 + id2ijg(isymij,isyma)
           koffv = 1 + igamsq(isymij,isymkl)
        
           ntota  = max(1,nbas(isyma))
           ntota1 = max(1,mbas1(isyma))
           ntotij = max(1,nmatij(isymij))

           call dgemm('T','N',nmatij(isymij),nmatij(isymkl),
     &          mbas1(isyma),factor,gaijd(koffg),ntota,
     &          work(koffr),ntota1,one,vijkl(koffv),ntotij)
         end if
c
      end do

      if (locdbg) then
        write(lupri,*)'DEBUG: vijkl, idelta,isymd =',idelta,isymd
        write(lupri,*)'norm^2 = ',ddot(ngamsq(isymhcs),vijkl,1,vijkl,1)
        do isyma = 1, nsym
         isymij = muld2h(isymd,isyma)
         isymkl = muld2h(isymhcs,isymij)
c        isymkl = isymij
         write(lupri,*) 'isymij,isymkl:',isymij,isymkl
         write(lupri,*)'CC2-R12 <V>'
         call output(vijkl(igamsq(isymij,isymkl)+1),1,nmatij(isymij),1,
     &        nmatij(isymkl),nmatij(isymij),nmatij(isymkl),
     &        1,lupri)
        end do
        write(lupri,*)'xr in mkvkl: ', xr
      end if

      call qexit('r12mkvkl')
      end 
*====================================================================*
