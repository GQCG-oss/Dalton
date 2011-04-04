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
*======================================================================*
      subroutine cc_r12vunpack(vijkl,isymv,vsing,vtrip,lprojv,
     &                         nr12orb,nrhforb)
*----------------------------------------------------------------------*
c     purpose: make V^{it jt }_{kl} = P^{kl}_{ij} * V^{it jt }_{kl}
c              and set it equal to the vector function omega for s = 0,1
c              (no symmetry is used!)
c
c     H. Fliegl, C. Haettig spring 2003 
c
c     modified by C. Neiss, summer 2005
c----------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccsdsym.h>     
#include <ccorb.h>

      logical locdbg,lprojv
      parameter(locdbg = .false.)
#if defined (SYS_CRAY)
      real half
#else
      double precision half
#endif
      parameter ( half = 0.5d0 )

      integer idxkl,idxlk,idxij,idxijkl,idxijlk,nrhftria,idxjikl,idxjilk
      integer isymkl,isymij,isyml,isymk,isymi,isymj,idxji
      integer it,jt,kt,lt,kl,ij,klij,isymv,icoun1,icoun2,isym
      integer nrhforb(8),irhforb(8),nij(8),iij(8,8),iijkl(8,8)
      integer nr12orb(8),ir12orb(8),nkl(8),ikl(8,8),nr12t
      integer index
#if defined (SYS_CRAY)
      real vijkl(*),vsing(*)
      real vtrip(*),ff
#else
      double precision vijkl(*),vsing(*)
      double precision vtrip(*),ff
#endif

      index(i,j) = max(i,j)*(max(i,j) - 3)/2 + i + j

      call qenter('cc_r12vunpack')
C
      ff = 1.0d0/sqrt(2.0d0)
C
      nr12t  = 0
      icoun1 = 0
      icoun2 = 0
      do isym = 1, nsym
        irhforb(isym) = icoun1
        ir12orb(isym) = icoun2
        icoun1 = icoun1 + nrhforb(isym)
        icoun2 = icoun2 + nr12orb(isym)
        nr12t  = nr12t  + nr12orb(isym)
      end do

      do isym = 1, nsym
        icoun1 = 0
        icoun2 = 0
        do isymj = 1, nsym
          isymi = muld2h(isymj,isym)
          iij(isymi,isymj) = icoun1
          ikl(isymi,isymj) = icoun2
          icoun1 = icoun1 + nrhforb(isymi)*nrhforb(isymj)
          icoun2 = icoun2 + nr12orb(isymi)*nr12orb(isymj)
        end do
        nij(isym) = icoun1
        nkl(isym) = icoun2
      end do

      do isym = 1, nsym
        icoun1 = 0
        do isymij = 1, nsym
          isymkl = muld2h(isymij,isym)
          iijkl(isymij,isymkl) = icoun1
          icoun1 = icoun1 + nij(isymij)*nkl(isymkl)
        end do
      end do
C
      nrhftria = nr12t*(nr12t+1)/2
      call dzero(vsing,nrhftria*nrhftria) 
      call dzero(vtrip,nrhftria*nrhftria) 

      do isymkl = 1, nsym
         isymij = muld2h(isymkl,isymv)
         do isyml =1, nsym
            isymk = muld2h(isymkl,isyml)

            do k = 1, nr12orb(isymk)
               kt = ir12orb(isymk) + k
               do l = 1, nr12orb(isyml)
                  lt = ir12orb(isyml) + l
                  if (lt.le.kt) then
                     kl = index(kt,lt)
                     idxkl=ikl(isymk,isyml)+nr12orb(isymk)*(l-1)+k
                     idxlk=ikl(isyml,isymk)+nr12orb(isyml)*(k-1)+l
                     do isymi =1, nsym
                        isymj = muld2h(isymij,isymi)

                        do j = 1, nrhforb(isymj)
                           jt = irhforb(isymj)+j                       
                           do i = 1, nrhforb(isymi)
                              it = irhforb(isymi)+i
                              if (jt.le.it) then
                                 ij = index(it,jt) 
                                 klij = nrhftria*(ij-1) + kl
                                 idxij = iij(isymi,isymj)+
     &                                nrhforb(isymi)*(j-1) + i
                                 idxijkl = iijkl(isymij,isymkl)+
     &                                nij(isymij)*(idxkl-1)+idxij
                                 idxijlk = iijkl(isymij,isymkl)+ 
     &                                nij(isymij)*(idxlk-1)+idxij

                                 vsing(klij)=vijkl(idxijkl)+
     &                                vijkl(idxijlk)
                                 vtrip(klij)=vijkl(idxijkl)-
     &                                vijkl(idxijlk)
c                            -----------------------------------------  
c                             V is not symmetric in ik <-> jl, thus
c                                 use 0.5 (V(ijkl) + V(jilk))
c                            -----------------------------------------
                                 if (lprojv) then
                                   idxji = iij(isymj,isymi)+
     &                                  nrhforb(isymj)*(i-1) + j
                                   idxjilk = iijkl(isymij,isymkl)+
     &                                  nij(isymij)*(idxlk-1)+idxji
                                   idxjikl = iijkl(isymij,isymkl)+
     &                                  nij(isymij)*(idxkl-1)+idxji
                                   vsing(klij) = half*(vsing(klij)
     &                                  + vijkl(idxjikl)+vijkl(idxjilk))
                                   vtrip(klij) = half*(vtrip(klij)
     &                                  + vijkl(idxjilk)-vijkl(idxjikl))
                                 end if

                                 if (it.eq.jt) 
     &                                vsing(klij)=ff*vsing(klij)
                                 if (kt.eq.lt) 
     &                                vsing(klij)=ff*vsing(klij)
                                 if (it.eq.jt) 
     &                                vtrip(klij)=ff*vtrip(klij)
                                 if (kt.eq.lt) 
     &                                vtrip(klij)=ff*vtrip(klij)
                              end if
                           end do
                        end do
                        
                     end do
                  end if

               end do
            end do
         end do
      end do

      if (locdbg) then
         write(lupri,*)'in cc_r12vunpack:'
         write(lupri,*)'Vsing'
         call output(vsing,1,nrhftria,1,nrhftria,nrhftria,
     &        nrhftria,1,lupri)     
         write(lupri,*)'Vtrip'
         call output(vtrip,1,nrhftria,1,nrhftria,nrhftria,
     &        nrhftria,1,lupri)
      end if

      call qexit('cc_r12vunpack')
      end

*=====================================================================* 
      subroutine cc_r12vpack(vijkl,isymv,vsing,vtrip,nr12orb,nrhforb,
     &                       nijkl)
*---------------------------------------------------------------------*
c     purpose: transform omega for s = 0,1 back to omega(ki,lj)
c
c     modified by C. Neiss, summer 2005
c     see CCR12PCK for further comments
c---------------------------------------------------------------------
      implicit none
#include <ccsdsym.h>     
#include <priunit.h>
#include <ccorb.h>
      integer idxkl,idxlk,idxij,idxijkl,idxijlk,nrhftria,idxjikl,idxjilk
      integer isymkl,isymij,isyml,isymk,isymi,isymj,idxji
      integer it,jt,kt,lt,index,kl,ij,klij,isymv,isym,icoun1,icoun2
      integer nr12orb(8),nrhforb(8),irhforb(8),ir12orb(8),
     &        nij(8),nkl(8),iij(8,8),ikl(8,8),nijkl(8),iijkl(8,8),nr12t
#if defined (SYS_CRAY)
      real vijkl(*),vsing(*)
      real vtrip(*),ff
      real ffkl,half
#else
      double precision vijkl(*),vsing(*)
      double precision vtrip(*),ff
      double precision ffkl,half
#endif
      logical locdbg
      parameter(locdbg = .false.)
      parameter (half = 0.5D0)

      index(i,j) = max(i,j)*(max(i,j) - 3)/2 + i + j
      
      call qenter('cc_r12vpack')
C
      nr12t  = 0
      icoun1 = 0
      icoun2 = 0
      do isym = 1, nsym
        irhforb(isym) = icoun1
        ir12orb(isym) = icoun2
        icoun1 = icoun1 + nrhforb(isym)
        icoun2 = icoun2 + nr12orb(isym)
        nr12t  = nr12t  + nr12orb(isym)
      end do

      do isym = 1, nsym
        icoun1 = 0
        icoun2 = 0
        do isymj = 1, nsym
          isymi = muld2h(isymj,isym)
          iij(isymi,isymj) = icoun1
          ikl(isymi,isymj) = icoun2
          icoun1 = icoun1 + nrhforb(isymi)*nrhforb(isymj)
          icoun2 = icoun2 + nr12orb(isymi)*nr12orb(isymj)
        end do
        nij(isym) = icoun1
        nkl(isym) = icoun2
      end do

      do isym = 1, nsym
        icoun1 = 0
        do isymij = 1, nsym
          isymkl = muld2h(isymij,isym)
          iijkl(isymij,isymkl) = icoun1
          icoun1 = icoun1 + nij(isymij)*nkl(isymkl)
        end do
        nijkl(isym) = icoun1
      end do
C
      nrhftria = nr12t*(nr12t+1)/2
      call dzero(vijkl,nijkl(isymv))

      do isymkl = 1, nsym
         isymij = muld2h(isymkl,isymv)
         do isyml =1, nsym
            isymk = muld2h(isymkl,isyml)

            do k = 1, nr12orb(isymk)
               kt = ir12orb(isymk) + k
               do l = 1, nr12orb(isyml)
                  lt = ir12orb(isyml) + l

                  if (lt.eq.kt) then
                     ffkl = sqrt(2d0)
                  else
                     ffkl = 1d0
                  end if

                  if (lt.le.kt) then
                     kl = index(kt,lt)
                     idxkl=ikl(isymk,isyml)+nr12orb(isymk)*(l-1)+k
                     idxlk=ikl(isyml,isymk)+nr12orb(isyml)*(k-1)+l

                     do isymi =1, nsym
                        isymj = muld2h(isymij,isymi)

                        do j = 1, nrhforb(isymj)
                           jt = irhforb(isymj)+j                       
                           do i = 1, nrhforb(isymi)
                              it = irhforb(isymi)+i

                              if (jt.eq.it) then
                                 ff = ffkl * sqrt(2d0)
                              else
                                 ff = ffkl
                              end if

                              if (jt.le.it) then
                                 ij = index(it,jt) 
                                 klij = nrhftria*(ij-1)+kl
                                 idxij = iij(isymi,isymj)+
     &                                nrhforb(isymi)*(j-1) + i
                                 idxijkl = iijkl(isymij,isymkl)+
     &                                nij(isymij)*(idxkl-1)+idxij
                                 idxijlk = iijkl(isymij,isymkl)+ 
     &                                nij(isymij)*(idxlk-1)+idxij
                                 idxji = iij(isymj,isymi)+
     &                                nrhforb(isymj)*(i-1) + j
                                 idxjilk = iijkl(isymij,isymkl)+
     &                                nij(isymij)*(idxlk-1)+idxji
                                 idxjikl = iijkl(isymij,isymkl)+
     &                                nij(isymij)*(idxkl-1)+idxji
                                 vijkl(idxijkl) = ff*half*(vsing(klij)
     &                                          +       vtrip(klij))
                                 vijkl(idxijlk) = ff*half*(vsing(klij)
     &                                          -       vtrip(klij))
                                 vijkl(idxjilk) = vijkl(idxijkl)
                                 vijkl(idxjikl) = vijkl(idxijlk)
                              end if

                           end do
                        end do

                     end do

                  end if

               end do
            end do

         end do
      end do

      if (locdbg) then
        write(lupri,*) 'Result in CC_R12VPACK:'
        do isymij = 1, nsym
          isymkl = muld2h(isymij,isymv)
          write(lupri,*) 'Symmetry block (ij,kl):',
     &             isymij,isymkl
          if (nkl(isymkl).eq.0 .or. nij(isymij).eq.0) then
            write(lupri,*) 'This symmetry is empty'
          else
            call output(vijkl(iijkl(isymij,isymkl)+1),1,nij(isymij),
     &                  1,nkl(isymkl),nij(isymij),nkl(isymkl),1,lupri)
          end if
        end do
      end if

      call qexit('cc_r12vpack')
      end

*=====================================================================*
      subroutine ccrhs_ep_old(vijkl,isymv,lcont,cont,work,lwork,
     &                    iamp,fc12am,lufc12,frho12,lufr12,
     &                    ifile,aproxr12,basscl1,basscl2)
*---------------------------------------------------------------------*
c     purpose: make the CC2-R12 vector function omega(ki,lj)
c
c     H. Fliegl, (WK/UniKA/02-05-2003)
c
c     C. Neiss   April 2005:
c     adapted for left hand transformations
c     basscl1 = brascl for right hand transf. 
c             = 0.5*ketscl for left hand transf.
c     basscl2 = ketscl for right hand transf.
c             = 2*brascl for left hand transf.
c
c     NOTE: Does NOT work with .R12ORB
c---------------------------------------------------------------------
      implicit none
#include <ccsdsym.h>
#include <priunit.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <dummy.h>
#include <r12int.h>
#include <ccr12int.h>

      logical lprojv,ldum,lcont
      integer lwork,isymv,nrhftria,ksing,ktrip,kend1,lwrk1     
      integer kbsing,kbtrip,kcsing,kctrip,kevl,kxsing,kxtrip,ij
      integer idum,ijcsvc,ijctvc,ijvsvc,ijvtvc,n2,lunit,ifile,kpck
      integer lufc12,lufr12,ian,iap,kpckv,kpckc,kbmat,kxint
      integer iamp,iopt
      integer nkilj(8),nijkl(8)
      character*(*) frho12, fc12am
      character*3 aprox
      character*10 model
      character*(*) aproxr12
#if defined (SYS_CRAY)
      real vijkl(*),work(*), cont(2), ddot, basscl1, basscl2
#else
      double precision vijkl(*),work(*), cont(2), ddot, basscl1,
     &                 basscl2
#endif
      logical locdbg
      parameter (locdbg = .false.)

      call qenter('ccrhs_ep')

      if (locdbg) write(lupri,*) 'Entered CCRHS_EP'

      nrhftria = nrhftb*(nrhftb+1)/2
      n2 = nrhftria * nrhftria

      ksing  = 1
      ktrip  = ksing  + n2
      kbsing = ktrip  + n2
      kbtrip = kbsing + n2
      kxsing = kbtrip + n2
      kxtrip = kxsing + n2
      kcsing = kxtrip + n2
      kctrip = kcsing + n2
      kevl   = kctrip + n2
      kpckc  = kevl   + nrhftria 
      kbmat  = kpckc  + ntr12am(isymv)
      kxint  = kbmat  + nr12r12p(1)
      kend1  = kxint  + nr12r12p(1)

      if (iamp .EQ. 1) then
        kpck  = kend1
        kend1 = kpck + ntr12am(isymv)
      end if
      if (lcont) then
        kpckv = kend1
        kend1 = kpckv + ntr12am(isymv)
      end if

      lwrk1  = lwork  - kend1
      if (lwrk1 .lt. 0) then
         call quit('Insufficient work space in ccrhs_ep')
      end if

      call dzero(work(kcsing),n2)
      call dzero(work(kctrip),n2)
     
c     ----------------------------
c     get V and omega for s = 0,1
c     ----------------------------
      lprojv = .true.
      call cc_r12vunpack(vijkl,isymv,work(ksing),work(ktrip),lprojv,
     &                   nrhfb,nrhf)
      if (lcont) then
        call ccr12pck(work(kpckv),isymv,work(ksing),work(ktrip),nrhfb,
     &                nrhf,nkilj)
      end if
      if (locdbg) then
c        call around('Vtilde integrals (singlet):')
c        call output(work(ksing),1,nrhftria,1,nrhftria,
c    &             nrhftria,nrhftria,1,lupri)
c        call around('Vtilde integrals (triplet):')
c        call output(work(ktrip),1,nrhftria,1,nrhftria,
c    &             nrhftria,nrhftria,1,lupri)
         call cc_r12vpack(vijkl,isymv,work(ksing),work(ktrip),nrhfb,
     &                    nrhf,nijkl)
         call around('Vtilde on entry (after 0.5*P_kl^ij):')
         call cc_prsqr12(vijkl,isymv,'T',1,.false.)
      end if 

c     ----------------------
c     read orbital energies 
c     ---------------------
      lunit = -1
      ldum = .false.
      call gpopen(lunit,fccr12e,'old',' ','formatted',
     &                   idum,ldum)
      read(lunit,'(4e30.20)') (work(kevl+ij), ij = 0, nrhftria-1)
      call gpclose(lunit,'KEEP')

c     ----------------     
c     read B matrices 
c     ----------------
      lunit = -1
      call gpopen(lunit,fccr12b,'old',' ','unformatted',
     &                   idum,ldum)
 8888 read(lunit) ian, iap, aprox
      read(lunit) (work(kbmat+i), i=0, nr12r12p(1)-1 )
      if ((ian.ne.ianr12).or.(iap.ne.iapr12)) goto 8888
      call gpclose(lunit,'KEEP')
      if (aproxr12(1:3).ne.aprox) then
        write(lupri,*)'aproxr12 =', aproxr12(1:3)
        write(lupri,*)'aprox    =', aprox(1:3)
        call quit('>>> Warning: inconsistent R12 approximation')
      end if
      call ccr12unpck(work(kbmat),1,work(kbsing),work(kbtrip),nrhfb,
     &                nrhfb)
c     -------------------------
c     read R12 overlap matrices
c     ------------------------- 
      lunit = -1
      call gpopen(lunit,fccr12x,'old',' ','unformatted',
     &                   idum,ldum)
 9999 read(lunit) ian
      read(lunit) (work(kxint+i), i=0, nr12r12p(1)-1 )
      if (ian.ne.ianr12) goto 9999
      call gpclose(lunit,'KEEP')
      call ccr12unpck(work(kxint),1,work(kxsing),work(kxtrip),nrhfb,
     &                nrhfb)

      if (iamp .EQ. 1) then
c       --------------------
c       read trial vector R2
c       --------------------
        call cc_rvec(lufc12,fc12am,ntr12am(isymv),ntr12am(isymv),
     *               ifile,work(kpck))
        if (lcont) then
          call dcopy(ntr12am(isymv),work(kpck),1,work(kpckc),1)
        end if
c       if (locdbg) then
c         call cc_prpr12(work(kpck),isymv,1,.TRUE.)
c       end if
        call cclr_diasclr12(work(kpck),basscl2,isymv)
        call ccr12unpck(work(kpck),isymv,work(kcsing),work(kctrip),
     &                  nrhfb,nrhf)
      else 
c       --------------------
c       read R12 amplitudes
c       --------------------
c       lunit = -1 
c       call gpopen(lunit,fccr12c,'unknown',' ','unformatted',
c    &                   idum,ldum)
c       read(lunit) (work(kpckc+i), i = 0, ntr12am(1)-1 )
c       call gpclose(lunit,'KEEP')
        if (isymv.ne.1) call quit('Symmetry error in CCRHS_EP')
        iopt = 32
        call cc_rdrsp('R0 ',0,1,iopt,model,dummy,work(kpckc))
        call ccr12unpck(work(kpckc),isymv,work(kcsing),work(kctrip),
     &                  nrhfb,nrhf)
      end if 

      if (lcont) then 
        cont(1) = ddot(ntr12am(isymv),work(kpckc),1,work(kpckv),1)
      end if

c     if (locdbg) then
c        call around('r12 vector function (singlet)')
c        call output(work(ksing),1,nrhftria,1,nrhftria,
c    &             nrhftria,nrhftria,1,lupri)
c        call around('r12 vector function (triplet)')
c        call output(work(ktrip),1,nrhftria,1,nrhftria,
c    &             nrhftria,nrhftria,1,lupri)
c        call around('r12 amplitudes (singlet)')
c        call output(work(kcsing),1,nrhftria,1,nrhftria,
c    &             nrhftria,nrhftria,1,lupri)
c        call around('r12 amplitudes (triplet)')
c        call output(work(kctrip),1,nrhftria,1,nrhftria,
c    &             nrhftria,nrhftria,1,lupri)
c        call around('b matrix (singlet)')
c        call output(work(kbsing),1,nrhftria,1,nrhftria,
c    &             nrhftria,nrhftria,1,lupri)
c        call around('b matrix (triplet)')
c        call output(work(kbtrip),1,nrhftria,1,nrhftria,
c    &             nrhftria,nrhftria,1,lupri)
c        call around('x matrix (singlet)')
c        call output(work(kxsing),1,nrhftria,1,nrhftria,
c    &             nrhftria,nrhftria,1,lupri)
c        call around('x matrix (triplet)')
c        call output(work(kxtrip),1,nrhftria,1,nrhftria,
c    &             nrhftria,nrhftria,1,lupri)
c        call around('orbital energies')
c        call outpak(work(kevl),nrhft,1,lupri)
c     end if
      if ( debug ) then
         write(lupri,*) 'ccrhs_ep> norm^2(omega_12 singlet):',
     &   ddot(n2,work(ksing),1,work(ksing),1)
         write(lupri,*) 'ccrhs_ep> norm^2(omega_12 triplet):',
     &   ddot(n2,work(ktrip),1,work(ktrip),1)
         if (iamp.EQ.1) write(lupri,*) 'ccrhs_ep> norm^2(C12 packed):',
     &   ddot(ntr12am(isymv),work(kpck),1,work(kpck),1)
         write(lupri,*) 'ccrhs_ep> norm^2(C12 singlet):',
     &   ddot(n2,work(kcsing),1,work(kcsing),1)
         write(lupri,*) 'ccrhs_ep> norm^2(C12 triplet):',
     &   ddot(n2,work(kctrip),1,work(kctrip),1)
      end if

c     ------------------------------------------------- 
c     loop over pairs and add B*c to omega for s = 0,1
c     -------------------------------------------------
      ijcsvc = kcsing
      ijctvc = kctrip
      ijvsvc = ksing
      ijvtvc = ktrip
      do ij = 1, nrhftria
          call ccrhs_ep0(ij,nrhftria,work(kevl),
     &                   work(ijvsvc),work(ijvtvc),
     &                   work(ijcsvc),work(ijctvc),
     &                   work(kbsing),work(kbtrip),
     &                   work(kxsing),work(kxtrip))
        ijcsvc = ijcsvc + nrhftria
        ijctvc = ijctvc + nrhftria
        ijvsvc = ijvsvc + nrhftria
        ijvtvc = ijvtvc + nrhftria
      end do

c     --------------
c     print results:
c     --------------
      if (locdbg) then
c        call around('r12 vector function (singlet)')
c        call output(work(ksing),1,nrhftria,1,nrhftria,
c    *               nrhftria,nrhftria,1,lupri)
c        call around('r12 vector function (triplet)')
c        call output(work(ktrip),1,nrhftria,1,nrhftria,
c    *               nrhftria,nrhftria,1,lupri)
         call cc_r12vpack(vijkl,isymv,work(ksing),work(ktrip),nrhfb,
     &                    nrhf,nijkl)
         call around('OmegaR12 on exit')
         call cc_prsqr12(vijkl,isymv,'T',1,.false.)
         write(lupri,*) 'ccrhs_ep> norm^2(OmegaR12):',
     &   ddot(ntr12sq(isymv),vijkl,1,vijkl,1)
      end if
      if (locdbg) then
         write(lupri,*) 'ccrhs_ep> norm^2(omega_12 singlet):',
     &   ddot(n2,work(ksing),1,work(ksing),1)
         write(lupri,*) 'ccrhs_ep> norm^2(omega_12 triplet):',
     &   ddot(n2,work(ktrip),1,work(ktrip),1)
      end if
 
      if (iamp .EQ. 1) then
c       ------------------------------------------------------------
c       for jacobian right transformations scale diagonal to 
c       transform to pseudo-orthonormal basis
c       ------------------------------------------------------------
        call ccr12pck(work(kpck),isymv,work(ksing),work(ktrip),nrhfb,
     &                nrhf,nkilj)
        call cclr_diasclr12(work(kpck),basscl1,isymv)
        call cc_wvec(lufr12,frho12,ntr12am(isymv),ntr12am(isymv),
     *               ifile,work(kpck))
        if (locdbg) then
          call around('Jacobian transformed vector in CCRHS_EP')
          call cc_prpr12(work(kpck),isymv,1,.false.)
        end if 
      else
c       ------------------------------------------------------------
c       write result vector as ground state vector function to file:
c       ------------------------------------------------------------
        call ccr12pck(work(kpckc),isymv,work(ksing),work(ktrip),nrhfb,
     &                nrhf,nkilj)
        lunit = -1 
        call gpopen(lunit,'CC_OMEGAR12','unknown',' ','unformatted',
     &                   idum,ldum)
        write(lunit) (work(kpckc+i), i=0, ntr12am(isymv)-1 )
        call gpclose(lunit,'KEEP')
      end if
c
      if (lcont) then 
        call daxpy(ntr12am(isymv),-1.0d0,work(kpckv),1,work(kpck),1)
        cont(2) = ddot(ntr12am(isymv),work(kpckc),1,work(kpck),1)
      end if

      if (locdbg) write(lupri,*) 'Leaving CCRHS_EP'

      call qexit('ccrhs_ep')
      end
*=====================================================================* 
      subroutine ccrhs_ep0(ij,n,e,vsing,vtrip,csing,ctrip,
     &     bsing,btrip,xsing,xtrip)
*---------------------------------------------------------------------*
c     purpose: make the CC2-R12 vector function 
c     omega = sum_{m.le.n} B^{ij}(kl,mn)*c(mn,ij) for s = 0,1
c
c     H. Fliegl, (WK/UniKA/02-05-2003)
c---------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <r12int.h>
      integer n,ij,kl,mn
#if defined (SYS_CRAY)
      real e(n),vsing(n),vtrip(n),csing(n),ctrip(n),
     &                 bsing(n,n),btrip(n,n),xsing(n,n),xtrip(n,n), 
     &                 b,half,two,fact
#else
      double precision e(n),vsing(n),vtrip(n),csing(n),ctrip(n),
     &                 bsing(n,n),btrip(n,n),xsing(n,n),xtrip(n,n), 
     &                 b,half,two,fact
#endif
      parameter (half = 0.5d0, two = 2d0)

      call qenter('ccrhs_ep0')

        if (ianr12.eq.1) then
c         if (iapr12.le.2) then
          if (iapr12.le.2.or.iapr12.eq.4.or.iapr12.eq.6) then
            fact = 0.0D0
          else
            fact = half
          end if
        else if (ianr12.eq.2 .or. ianr12.eq.3) then
          if (iapr12.eq.2.or.iapr12.eq.5) then
            fact = 0.0D0
          else
            fact = half
          end if
        end if
c
      do kl = 1, n 
         do mn = 1, n
            b = - bsing(kl,mn) 
     &   + 
     &           fact * (e(kl) + e(mn) - two * e(ij)) * xsing(kl,mn)
            vsing(kl) = vsing(kl) + b * csing(mn)
            b = - btrip(kl,mn) 
     &    + 
     &           fact * (e(kl) + e(mn) - two * e(ij)) * xtrip(kl,mn)
            vtrip(kl) = vtrip(kl) + b * ctrip(mn)
         end do
      end do

      call qexit('ccrhs_ep0')
      end

*=====================================================================*
      subroutine ccrhs_ep(vijkl,isymv,lcont,cont,work,lwork,
     &                    iamp,fc12am,lufc12,frho12,lufr12,
     &                    ifile,aproxr12,basscl1,basscl2)
*---------------------------------------------------------------------*
c     purpose: make the CC2-R12 vector function omega(ki,lj):
c              take OmegaR12 = Vijkl and add OmegaR12_E' contribution,
c              write final OmegaR12 on disk.
c              alternative subroutine for old ccrhs_ep without 
c              using singlet/triplet matrices in matrix 
c              multiplication 
c
c     vijkl   = OmegaR12  
c     isymv   = symmetry of result
c     basscl1 = brascl for right hand transf. 
c             = 0.5*ketscl for left hand transf.
c     basscl2 = ketscl for right hand transf.
c             = 2*brascl for left hand transf.
c
c     Christian Neiss,  April. 2005,  based on old ccrhs_ep
c---------------------------------------------------------------------
      implicit none
#include <ccsdsym.h>
#include <priunit.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <dummy.h>
#include <r12int.h>
#include <ccr12int.h>

      logical ldum,lcont,lbmat
      integer lwork,isymv,kend1,lwrk1,kend2,lwrk2,kevl,kekl,keij     
      integer iamp,idum,lunit,ifile,ian,iap,isym,iopt
      integer lufc12,lufr12,kpck,kpckv,ktr12,kxintsq
      integer ktr12sq,kbmatsq,koffv,kofftr12,koffb,kbscr
      integer isymi,isymj,isymij,isymmn,isymkl,idxij
      integer nkl,nij,inmatkl(8),inmatij(8),norbtsx
      character*(*) frho12, fc12am
      character*3 cdum
      character*10 model
      character*3 aproxr12,aprox
#if defined (SYS_CRAY)
      real vijkl(*),work(*), cont(2), ddot, basscl1, basscl2, one,
     &     factor, half
#else
      double precision vijkl(*),work(*), cont(2), ddot, basscl1,
     &                 basscl2, one, factor, half
#endif
      logical locdbg
      parameter (locdbg = .false.)
      parameter (one = 1.0D0, half = 0.5D0)

      call qenter('ccrhs_ep')

      if (locdbg) write(lupri,*) 'Entered CCRHS_EP'

      ktr12   = 1 
      ktr12sq = ktr12   + ntr12am(isymv)
      kbmatsq = ktr12sq + ntr12sq(isymv)
      kpck    = kbmatsq + nr12r12sq(1)
      kend1   = kpck    + max(nr12r12p(1),ntr12am(isymv))

      if (lcont) then
        kpckv = kend1
        kend1 = kpckv   + ntr12am(isymv)
      end if

      lwrk1  = lwork  - kend1
      if (lwrk1 .lt. 0) then
         call quit('Insufficient work space in ccrhs_ep')
      end if

c     ---------------
c     get V (= Omega)
c     ---------------
      !apply 0.5*P_kl^ij:
      call cc_r12pklij(vijkl,isymv,'T',work(kend1),lwrk1)
      call dscal(ntr12sq(isymv),0.5D0,vijkl,1)
      if (lcont) then
        iopt = 1
        call ccr12pck2(work(kpckv),isymv,.false.,vijkl,'T',iopt)
      end if
      if (locdbg) then
         call around('OmegaR12 on entry (after 0.5*P_kl^ij):')
         call cc_prsqr12(vijkl,isymv,'T',1,.false.)
         write(lupri,*) 'ccrhs_ep> norm^2(OmegaR12):',
     &   ddot(ntr12sq(isymv),vijkl,1,vijkl,1)
      end if 

c     -----------------------------------
c     read R12 amplitudes or trial vector
c     -----------------------------------
      call cc_r12getct(work(ktr12sq),isymv,iamp,basscl2,.false.,'T',
     &                 lufc12,fc12am,ifile,cdum,idum,work(kend1),lwrk1)
      if (locdbg) then
        call around('R12 amplitudes')
        call cc_prsqr12(work(ktr12sq),isymv,'T',1,.false.)
      end if

      if (lcont) then
        iopt = 1
        call ccr12pck2(work(ktr12),isymv,.false.,work(ktr12sq),'T',
     &                 iopt) 
        cont(1) = ddot(ntr12am(isymv),work(ktr12),1,work(kpckv),1)
      end if

c     ----------------
c     read B matrices
c     ----------------
      lunit = -1
      call gpopen(lunit,fccr12b,'old',' ','unformatted',
     &                   idum,ldum)
 8888 read(lunit) ian, iap, aprox
      read(lunit) (work(kpck+i), i=0, nr12r12p(1)-1 )
      if ((ian.ne.ianr12).or.(iap.ne.iapr12)) goto 8888
      call gpclose(lunit,'KEEP')
      if (aproxr12(1:3).ne.aprox) then
        write(lupri,*)'aproxr12 =', aproxr12(1:3)
        write(lupri,*)'aprox    =', aprox(1:3)
        call quit('>>> Warning: inconsistent R12 approximation')
      end if
      iopt = 2
      call ccr12unpck2(work(kpck),1,work(kbmatsq),'N',iopt)
      call dscal(nr12r12sq(1),-one,work(kbmatsq),1)

      if (locdbg) then
        call around('B-Matrix read from file')
        do isymkl = 1, nsym
          isymmn = isymkl
          call output(work(kbmatsq+ir12r12sq(isymkl,isymmn)),1,
     &                nmatkl(isymkl),1,nmatkl(isymmn),nmatkl(isymkl),
     &                nmatkl(isymmn),1,lupri)
        end do 
      end if

c     ----------------------------------------------------
c     for some cases we are nearly done, since B does not
c     depend on (ij)
c     ----------------------------------------------------
      if (ianr12.eq.1) then
        if (iapr12.le.2 .or. iapr12.eq.4 .or. iapr12.eq.6) then
          lbmat = .false.
        else
          lbmat = .true.
          factor = half
        end if
      else
        if (iapr12.eq.2 .or. iapr12.eq.5) then
          lbmat = .false. 
        else
          lbmat = .true.
          factor = half
        end if
      end if

      if (.not. lbmat) then
        !matrix multiplication:
        call cc_r12xi2b(vijkl,'T',work(kbmatsq),1,'N',work(ktr12sq),
     &                  isymv,'T',one)
      else
c       ----------------------------------------------- 
c       in other cases: 
c        - read R12 overlap matrix, orbital energies
c        - loop over pairs (ij) and add B*c to OmegaR12
c       -----------------------------------------------
        !calculate # of indices (kl) / (ij) and offset over all symmetries:
        nkl = 0
        nij = 0
        do isym = 1, nsym
          inmatkl(isym) = nkl
          inmatij(isym) = nij
          nkl = nkl + nmatkl(isym)
          nij = nij + nmatij(isym)
        end do
c
        kbscr = kend1
        kxintsq = kbscr + nr12r12sq(1)
c       kevl  = kxintsq + nr12r12sq(1)
c       kekl  = kevl + norbts
        kekl  = kxintsq + nr12r12sq(1)
        keij  = kekl + nkl
        keij  = kekl + nkl
        kend2 = keij + nij
        lwrk2 = lwork - kend2
        if (lwrk2 .lt. 0) then
          call quit('Insufficient work space in ccrhs_ep')
        end if
c
        lunit = -1
        call gpopen(lunit,fccr12x,'old',' ','unformatted',
     &              idum,ldum)
 9999   read(lunit) ian
        read(lunit) (work(kpck+i), i=0, nr12r12p(1)-1 )
        if (ian.ne.ianr12) goto 9999
        call gpclose(lunit,'KEEP')
        iopt = 2
        call ccr12unpck2(work(kpck),1,work(kxintsq),'N',iopt)
        if (locdbg) then
          call around('R12 Overlap-Matrix read from file')
          do isymkl = 1, nsym
            isymmn = isymkl
            call output(work(kxintsq+ir12r12sq(isymkl,isymmn)),1,
     &                  nmatkl(isymkl),1,nmatkl(isymmn),nmatkl(isymkl),
     &                  nmatkl(isymmn),1,lupri)
          end do
        end if 
c
c       ------------------------------------
c       Read orbital energies.
c       -------------------------------------
c
        lunit = -1
        call gpopen(lunit,'SIRIFC','OLD',' ','UNFORMATTED',idummy,
     &              .false.)
        rewind lunit
        call mollab('FULLBAS ',lunit,lupri)
        read (lunit) idum,norbtsx 
        kevl = kend2
        kend2 = kevl + norbtsx
        lwrk2 = lwork - kend2
        if (lwrk2 .lt. 0) then
          call quit('Insufficient work space in ccrhs_ep')
        end if
        read (lunit) (work(kevl+i-1), i=1,norbtsx)
        call gpclose(lunit,'KEEP')
c
        call cc_r12epair(work(kekl),nkl,work(kevl),inmatkl,imatkl,nrhfb)
        call cc_r12epair(work(keij),nij,work(kevl),inmatij,imatij,nrhf)
c
        do isymij = 1, nsym
          isymmn = muld2h(isymv,isymij)
          isymkl = isymmn
c         if (locdbg) then
c           write(lupri,*) 'in loop over idxij: isymij = ',isymij
c         end if
          do idxij = 1, nmatij(isymij)
            koffb    = ir12r12sq(isymkl,isymmn)
            call cc_r12buildbmat(work(kbscr+koffb),isymkl,isymmn,
     &                           isymij,idxij,work(kbmatsq+koffb),
     &                           work(kxintsq+koffb),work(kekl),inmatkl,
     &                           work(keij),inmatij,factor)
            kofftr12 = ktr12sq + itr12sqt(isymij,isymmn) + idxij-1 
            koffv    = itr12sqt(isymij,isymkl) + idxij
            if (locdbg) then
              write(lupri,*) 'idxij = ',idxij
              write(lupri,*) 'R12 amplitudes part used:'
              do i = 1, nmatkl(isymmn)
                write(lupri,*) kofftr12+(i-1)*nmatij(isymij),
     &                       work(kofftr12+(i-1)*nmatij(isymij))
              end do
              write(lupri,*) 'Omega part used:'
              do i = 1, nmatkl(isymkl)
                write(lupri,*) koffv+(i-1)*nmatij(isymij),
     &                       vijkl(koffv+(i-1)*nmatij(isymij))
              end do
            end if
            call dgemv('N',nmatkl(isymkl),nmatkl(isymmn),
     &                 one,work(kbscr+koffb),max(nmatkl(isymkl),1),
     &                 work(kofftr12),max(nmatij(isymij),1),
     &                 one,vijkl(koffv),max(nmatij(isymij),1))
          end do
        end do
      end if

c     --------------
c     print results:
c     --------------
      if (locdbg) then
         call around('OmegaR12 on exit')
         call cc_prsqr12(vijkl,isymv,'T',1,.false.)
         write(lupri,*) 'ccrhs_ep> norm^2(OmegaR12):',
     &   ddot(ntr12sq(isymv),vijkl,1,vijkl,1)
      end if
 
      if (iamp .EQ. 1) then
c       ------------------------------------------------------------
c       for jacobian right transformations scale diagonal to 
c       transform to pseudo-orthonormal basis
c       ------------------------------------------------------------
        iopt = 1
        call ccr12pck2(work(kpck),isymv,.false.,vijkl,'T',iopt)
        call cclr_diasclr12(work(kpck),basscl1,isymv)
        call cc_wvec(lufr12,frho12,ntr12am(isymv),ntr12am(isymv),
     *               ifile,work(kpck))
        if (locdbg) then
          call around('Jacobian transformed vector in CCRHS_EP')
          call cc_prpr12(work(kpck),isymv,1,.false.)
        end if 
      else
c       ------------------------------------------------------------
c       write result vector as ground state vector function to file:
c       ------------------------------------------------------------
        iopt = 1
        call ccr12pck2(work(ktr12),isymv,.false.,vijkl,'T',iopt)
        lunit = -1 
        call gpopen(lunit,'CC_OMEGAR12','unknown',' ','unformatted',
     &                   idum,ldum)
        write(lunit) (work(ktr12+i), i=0, ntr12am(isymv)-1 )
        call gpclose(lunit,'KEEP')
      end if
c
      if (lcont) then 
        call daxpy(ntr12am(isymv),-one,work(kpckv),1,work(kpck),1)
        cont(2) = ddot(ntr12am(isymv),work(ktr12),1,work(kpck),1)
      end if

      if (locdbg) write(lupri,*) 'Leaving CCRHS_EP'

      call qexit('ccrhs_ep')
      end
*=====================================================================* 

*=====================================================================*
      subroutine cc_r12nxtam_old(omeg12,isym,tamp12,lcceq,
     &                       er12,work,lwork)
c----------------------------------------------------------------------
c     purpose: get new start R12  amplitudes for each iteration
c
c     H. Fliegl, C. Haettig, (WK/UniKA/08-05-2003)
c     modified by C. Neiss  April 2005, September 2005
c----------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <dummy.h>
#include <r12int.h>
#include <ccr12int.h>
#include <ccsdinp.h>
#include <second.h>

      logical locdbg,ldum,lcceq 
      parameter (locdbg = .false.)

      integer kbsing,kbtrip,lwork,nrhftria,n2,kxsing,kxtrip,kevl
      integer kend1,lwrk1,lunit,idum,ijv,ij,kvsing,kvtrip
      integer komgsing,komgtrip,kcsing,kctrip
      integer kscr,kvint,lu43,kepsij,kepsoff
      integer kbb,kuu,kpp,kqq,krr,ian,iap,iopt,isym
      integer nkilj(8),nij,isymij,inmatij(8),it,jt,isymi,isymj,norbtsx
#if defined (SYS_CRAY)
      real tamp12(*),omeg12(*),er12,work(*)
      real ensing,entrip,delta,xf,fs,ft,ws,wt,ddot,timnxtam
#else
      double precision tamp12(*),omeg12(*),er12,work(*)
      double precision ensing,entrip,delta,xf,fs,ft,ws,wt,ddot,timnxtam
#endif
      character*3 aprox
      character*10 model
      integer index
c
      index(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j
c
      call qenter('r12nxtam')
      if (locdbg) then
        write(lupri,*) 'Entered CC_R12NXTAM'
        call flshfo(lupri)
      end if
c
      timnxtam = second()
c
      nrhftria = nrhftb*(nrhftb+1)/2
      n2 = nrhftria * nrhftria
c
      nij = 0
      do isymij = 1, nsym
        inmatij(isymij) = nij
        nij = nij + nmatij(isymij)
      end do
c
      kend1  = 1
c     
      komgsing = kend1
      komgtrip = komgsing + n2
      kcsing = komgtrip + n2
      kctrip = kcsing + n2
      kbsing = kctrip + n2
      kbtrip = kbsing + n2
      kxsing = kbtrip + n2
      kxtrip = kxsing + n2
      kend1  = kxtrip + n2
     
      if (lcceq) then
        kvsing = kend1
        kvtrip = kvsing + n2
        kend1  = kvtrip + n2
      end if

      kevl   = kend1
      kbb    = kevl   + nrhftria 
      kuu    = kbb    + n2 
      kpp    = kuu    + n2
      kqq    = kpp    + n2
      krr    = kqq    + n2
      kend1  = krr    + n2
c
      kepsij = kend1
      kend1  = kepsij + nij
c
      kscr   = kend1
      kend1  = kscr + nr12r12p(1)
c
      lwrk1  = lwork  - kend1

      if (lwrk1 .lt. 0) then
         call quit('Insufficient work space in cc_r12nxtam')
      end if

c     -------------
c     test symmetry
c     -------------
      if (lcceq .and. (isym.ne.1)) then
        call quit('Symmetry error in CC_R12NXTAM')
      end if 

c     ----------------------
c     read orbital energies 
c     ---------------------
      lunit = -1
      ldum = .false.
      call gpopen(lunit,fccr12e,'old',' ','formatted',
     &                   idum,ldum)
      read(lunit,'(4e30.20)') (work(kevl+ij), ij = 0, nrhftria-1)
      call gpclose(lunit,'KEEP')
      if (locdbg) then
        write(lupri,*) 'Read orbital energies done:'
        call outpak(work(kevl),nrhftria,1,lupri)
        call flshfo(lupri)
      end if

      if (lcceq) then
c       ---------------    
c       read V matrices 
c       ---------------
        lunit = -1
        call gpopen(lunit,fccr12v,'old',' ','unformatted',
     &                     idum,ldum)
 7777   read(lunit) ian
        read(lunit) (work(kscr+i), i=0, ntr12am(1)-1)
        if (ian.ne.ianr12) goto 7777
        call gpclose(lunit,'KEEP')
        if (locdbg) then
          write(lupri,*) 'Read V-intermediate done'
          call flshfo(lupri)
        end if
        call ccr12unpck(work(kscr),1,work(kvsing),work(kvtrip),nrhfb,
     &                  nrhf)
        if (locdbg) then
          write(lupri,*) 'Unpack V-intermediate done'
          write(lupri,*) 'Singlet V-intermediate:'
          call output(work(kvsing),1,nrhftria,1,nrhftria,nrhftria,
     &              nrhftria,1,lupri)
          write(lupri,*) 'Triplet V-intermediate:'
          call output(work(kvtrip),1,nrhftria,1,nrhftria,nrhftria,
     &              nrhftria,1,lupri)
          call flshfo(lupri)
        end if
CCN
C        write(lupri,*) 'V-Matrix in CC_R12NXTAM:'
C        call ccr12unpck2(work(kscr),1,work(kend1),'T',1)
CCN
      end if

c     ---------------    
c     read B matrices
c     --------------- 
      lunit = -1
      call gpopen(lunit,fccr12b,'old',' ','unformatted',
     &                   idum,ldum)
 8888 read(lunit) ian, iap, aprox
      read(lunit) (work(kscr+i), i=0, nr12r12p(1)-1)
      if ((ian.ne.ianr12).or.(iap.ne.iapr12)) goto 8888
      call gpclose(lunit,'KEEP')
      call ccr12unpck(work(kscr),1,work(kbsing),work(kbtrip),nrhfb,
     &                nrhfb)

c     ----------------------------
c     read R12 overlap matrices X 
c     ----------------------------
      lunit = -1
      call gpopen(lunit,fccr12x,'old',' ','unformatted',
     &                   idum,ldum)
 9999 read(lunit) ian
      read(lunit) (work(kscr+i), i=0, nr12r12p(1)-1)
      if (ian.ne.ianr12) goto 9999
      call gpclose(lunit,'KEEP')
      call ccr12unpck(work(kscr),1,work(kxsing),work(kxtrip),nrhfb,
     &                nrhfb)

      if (lcceq) then
c       --------------------
c       read R12 amplitudes
c       --------------------
        iopt =32
        call cc_rdrsp('R0 ',0,1,iopt,model,dummy,tamp12)

        if (iprint .ge. 5) then
          write(lupri,529) 'Norm^2 of r12am     in NXTAM:',
     *                      ddot(ntr12am(1),tamp12,1,tamp12,1)
        end if

c       -------------------------------
c       read vector function omega_R12
c       -------------------------------
        lunit = -1 
        call gpopen(lunit,'CC_OMEGAR12','old',' ','unformatted',
     &              idum,ldum)
        read(lunit) (omeg12(i), i=1, ntr12am(1))
        call gpclose(lunit,'KEEP')   
chf
        write(lupri,529)'Norm^2 of Omegar12  in NXTAM:', 
     &     ddot(ntr12am(1),omeg12,1,omeg12,1)
chf
      end if
      !repack omega_R12
      call ccr12unpck(omeg12,isym,work(komgsing),work(komgtrip),nrhfb,
     &                nrhf)

      if (locdbg) then
        write(lupri,*) 'Singlet Omega-R12:'
        call output(work(komgsing),1,nrhftria,1,nrhftria,nrhftria,
     &              nrhftria,1,lupri)
        write(lupri,*) 'Triplet Omega-R12:'
        call output(work(komgtrip),1,nrhftria,1,nrhftria,nrhftria,
     &              nrhftria,1,lupri)
        call around('b in nxtam matrix (singlet)')
        call output(work(kbsing),1,nrhftria,1,nrhftria,
     &       nrhftria,nrhftria,1,lupri)
        call around('b in nxtam  matrix (triplet)')
        call output(work(kbtrip),1,nrhftria,1,nrhftria,
     &       nrhftria,nrhftria,1,lupri)
        call around('x in nxtam matrix (singlet)')
        call output(work(kxsing),1,nrhftria,1,nrhftria,
     &       nrhftria,nrhftria,1,lupri)
        call around('x in nxtam matrix (triplet)')
        call output(work(kxtrip),1,nrhftria,1,nrhftria,
     &       nrhftria,nrhftria,1,lupri)
        call flshfo(lupri)
      end if

 529  format(7x,a,d24.10)

      if (locdbg) then
         write(lupri,529) 'norm^2 of omega r12 singlet:',
     *        ddot(n2,work(komgsing),1,work(komgsing),1)
         write(lupri,529) 'norm^2 of omega r12 triplet:',
     *        ddot(n2,work(komgtrip),1,work(komgtrip),1)
      end if
         
      call dzero(work(kcsing),n2)
      call dzero(work(kctrip),n2)

      if (locdbg) write(lupri,*) 'R12SVD,R12DIA:',R12SVD,R12DIA
c     ------------------------------------------------------ 
c      get CC2-R12 contributions to the ground state energy
c
c     xf  =0.5d0 for approximation A'
c     xf = 0 for for approximation A
c     ------------------------------------------------------
      xf = 0.0d0
      delta = 0d0
      if (ianr12.eq.1) then
        if (iapr12.le.2.or.iapr12.eq.4.or.iapr12.eq.6) then
          xf = 0
        else
          xf = 0.5d0
        end if
      else if (ianr12.eq.2 .or. ianr12.eq.3) then
c       if (iapr12.le.2) then
        if (iapr12.eq.2.or.iapr12.eq.5) then
          xf = 0
        else
          xf = 0.5d0
        end if
      endif
c     -----------------------------------------------
c     Get orbital energy pairs for OCCUPIED orbitals:
c     -----------------------------------------------
      lunit = -1
      call gpopen(lunit,'SIRIFC','OLD',' ','UNFORMATTED',idum,
     &            .false.)
      rewind lunit
      call mollab('FULLBAS ',lunit,lupri)
      read (lunit) idum,norbtsx
      kend1 = kscr + norbtsx
      lwrk1 = lwork - kend1
      if (lwrk1.lt.0) then
        call quit('Memory allocation error in CC_R12NXTAM')
      end if
      read (lunit) (work(kscr+i-1), i=1, norbtsx)
      call gpclose(lunit,'KEEP')
c
      call cc_r12epair(work(kepsij),nij,work(kscr),inmatij,imatij,nrhf)
c
      do isymj = 1, nsym
       do isymi = 1, nsym
        do j = 1, nrhf(isymj)
         jt = irhf(isymj) + j
         do i = 1, nrhf(isymi)
          it = irhf(isymi) + i
          if (it.le.jt) then
           ij = index(it,jt)
           ijv = nrhftria*(ij-1)
           isymij = muld2h(isymi,isymj)
           kepsoff = inmatij(isymij)+imatij(isymi,isymj)+
     &                 nrhf(isymi)*(j-1)+i-1
           if (R12SVD) then
            CALL R12MP2(work(komgsing+ijv),work(komgtrip+ijv),
     &                  work(komgsing+ijv),work(komgtrip+ijv),
     &                  work(kxsing),work(kxtrip),
     &                  work(kepsij+kepsoff),nrhftb,nrhftria,
     &                  fs,ft,work(kbb),work(kuu),work(kpp),
     &                  work(kqq),
     &                  work(kbsing),work(kbtrip),
     &                  work(kevl),xf,work(krr),ij,ws,wt,
     &                  delta,work(kcsing+ijv),work(kctrip+ijv))
           else if (R12DIA) then
            CALL MP2R12(work(komgsing+ijv),work(komgtrip+ijv),
     &                  work(komgsing+ijv),work(komgtrip+ijv),
     &                  work(kxsing),work(kxtrip),
     &                  work(kepsij+kepsoff),nrhftb,nrhftria,
     &                  fs,ft,work(kbb),work(kuu),work(kpp),
     &                  work(kqq),
     &                  work(kbsing),work(kbtrip),
     &                  work(kevl),xf,work(krr),ij,ws,wt,
     &                  delta,work(kcsing+ijv),work(kctrip+ijv))
           end if
          end if
         end do
        end do
       end do
      end do

c     ------------------------------------------------
c     copy update for amplitudes into omgsing/omgtrip
c     ------------------------------------------------
      call dcopy(n2,work(kcsing),1,work(komgsing),1)
      call dcopy(n2,work(kctrip),1,work(komgtrip),1)
      if (locdbg) then
        write(lupri,*) 'Update for Singlet R12 amplitudes:'
        call output(work(kcsing),1,nrhftria,1,nrhftria,nrhftria,
     &              nrhftria,1,lupri)
        write(lupri,*) 'Update for Triplet R12 amplitudes:'
        call output(work(kctrip),1,nrhftria,1,nrhftria,nrhftria,
     &              nrhftria,1,lupri)
      end if

      if (lcceq) then
        call ccr12unpck(tamp12,1,work(kcsing),work(kctrip),nrhfb,
     &                  nrhf)

        if (locdbg) then
          write(lupri,*) 'Singlet R12 amplitudes:'
          call output(work(kcsing),1,nrhftria,1,nrhftria,nrhftria,
     &                nrhftria,1,lupri)
          write(lupri,*) 'Triplet R12 amplitudes:'
          call output(work(kctrip),1,nrhftria,1,nrhftria,nrhftria,
     &                nrhftria,1,lupri)
        end if

c       ----------------------------------------------------
c       add old amplitudes to update to get new amplitudes
c       ----------------------------------------------------
        call daxpy(n2,1.0d0,work(kcsing),1,work(komgsing),1)
        call daxpy(n2,1.0d0,work(kctrip),1,work(komgtrip),1)

        ensing = ddot(n2,work(kvsing),1,work(komgsing),1)
        entrip = 3.0D0*ddot(n2,work(kvtrip),1,work(komgtrip),1)
        er12 = ensing + entrip

        if (locdbg) then
          write(lupri,*) 'Singlet contribution to R12 energy:',ensing
          write(lupri,*) 'Triplet contribution to R12 energy:',entrip
        end if

      end if ! lcceq

      if (locdbg) then
        write(lupri,*) 'Updated Singlet R12 amplitudes:'
        call output(work(komgsing),1,nrhftria,1,nrhftria,nrhftria,
     &              nrhftria,1,lupri)
        write(lupri,*) 'Updated Triplet R12 amplitudes:'
        call output(work(komgtrip),1,nrhftria,1,nrhftria,nrhftria,
     &              nrhftria,1,lupri)
      end if

      call ccr12pck(omeg12,isym,work(komgsing),work(komgtrip),nrhfb,
     &              nrhf,nkilj)

      if (locdbg) write(lupri,*) 'Leaving CC_R12NXTAM'
c
      timnxtam = second() - timnxtam
c     write(lupri,111) 'CC_R12NXTAM', timnxtam
  111 FORMAT(7x,'Time used in',2x,A12,2x,': ',f10.2,' seconds')
c
      call qexit('r12nxtam')
      end
*=====================================================================*
      subroutine cc_r12metric_old(isymr,basscl1,basscl2,work,lwork,
     &                        fc2am,lufc2,fc12am,lufc12,fs12am,lufs12,
     &                        fs2am,lufs2,ifile)
*---------------------------------------------------------------------*
c     purpose: make the CC2-R12 S*R contibution for excitation 
c              energies
c
c     H. Fliegl, C. Haettig summer 2003
c
c     nondiagonal elements added for ansatz 2, spring 2005
c     --> S*R written on fs2am
c
c     C. Neiss, 05.04.2005: adapted for left hand transformation
c     basscl1 = brascl for right hand transf.
c             = 0.5*ketscl for left hand transf.
c     basscl2 = ketscl for right hand transf.
c             = 2*brascl for left hand transf.
c---------------------------------------------------------------------
      implicit none
#include <ccsdsym.h>
#include <priunit.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <r12int.h>
#include <ccr12int.h>

      logical ldum,locdbg,lprojv
      parameter (locdbg = .false.)

      character*(*) fc12am,fs12am,fc2am,fs2am

      integer lwork,isymr,nrhftria,kend1,lwrk1,lufc12,lufs12,ifile     
      integer lufc2,lufs2
      integer krsing,krtrip,kxsing,kxtrip,ij,ksing,ktrip,kpck,kxint
      integer idum,n2,lunit,ijrsvc,ijrtvc,ijsvc,ijtvc
      integer ian, iap,ks2,krabkl,krmnab,ksr
      integer nkilj(8)
#if defined (SYS_CRAY)
      real work(*), ddot,two,one,half,basscl1,basscl2
#else
      double precision work(*), ddot,two,one,half,basscl1,basscl2
#endif
      parameter (one = 1.0d0, two = 2.0d0, half = 0.5d0)

      call qenter('cc_r12metric')
      if (locdbg) then
        write(lupri,*) 'Entered CC_R12METRIC'
      end if

      nrhftria = nrhftb*(nrhftb+1)/2
      n2 = nrhftria * nrhftria

      ksing  = 1
      ktrip  = ksing  + n2
      kxsing = ktrip  + n2
      kxtrip = kxsing + n2
      krsing = kxtrip + n2
      krtrip = krsing + n2
      kxint  = krtrip + n2
      kend1  = kxint  + nr12r12p(1)

      kpck  = kend1
      kend1 = kpck + ntr12am(isymr)
      if (ianr12.eq.2) then
        krabkl = kend1
        krmnab = krabkl + nt2r12(1)
        ksr    = krmnab + nt2am(isymr)
        kend1  = ksr + ntr12sq(isymr)
      end if

      lwrk1  = lwork  - kend1

      if (lwrk1 .lt. 0) then
         write(lupri,*) 'lwork, lwrk1: ',lwork, lwrk1
         call quit('Insufficient work space in cc_r12metric')
      end if

      call dzero(work(krsing),n2)
      call dzero(work(krtrip),n2)
     
      call dzero(work(krsing),n2)
      call dzero(work(krtrip),n2)

c     -------------------------
c     read R12 overlap matrices
c     ------------------------- 
      lunit = -1
      call gpopen(lunit,fccr12x,'unknown',' ','unformatted',
     &                   idum,ldum)
 9999 read(lunit) ian
      read(lunit) (work(kxint+i), i=0, nr12r12p(1)-1 )
      if (ian.ne.ianr12) goto 9999
      call gpclose(lunit,'KEEP')
      call ccr12unpck(work(kxint),1,work(kxsing),work(kxtrip),nrhfb,
     &                nrhfb)

c     ------------------------
c     read try vector R^ij_kl
c     ------------------------
      lunit = -1
      call cc_rvec(lufc12,fc12am,ntr12am(isymr),ntr12am(isymr),
     *             ifile,work(kpck))
      call cclr_diasclr12(work(kpck),basscl2,isymr)
      call ccr12unpck(work(kpck),isymr,work(krsing),work(krtrip),
     &                nrhfb,nrhf)

      if (locdbg) then
         write(lupri,*) 'fc12am,lufc12,ifile:',fc12am,lufc12,ifile
         write(lupri,*) 'cc_r12metric> norm^2(R12 packed):',
     &   ddot(ntr12am(isymr),work(kpck),1,work(kpck),1)
         write(lupri,*) 'cc_r12metric> norm^2(R12 singlet):',
     &   ddot(n2,work(krsing),1,work(krsing),1)
         write(lupri,*) 'cc_r12metric> norm^2(R12 triplet):',
     &   ddot(n2,work(krtrip),1,work(krtrip),1)
      end if
      if (locdbg) then
         call around('r12 try amplitudes (singlet)')
         call output(work(krsing),1,nrhftria,1,nrhftria,
     &             nrhftria,nrhftria,1,lupri)
         call around('r12 try amplitudes (triplet)')
         call output(work(krtrip),1,nrhftria,1,nrhftria,
     &             nrhftria,nrhftria,1,lupri)
         call around('x matrix (singlet)')
         call output(work(kxsing),1,nrhftria,1,nrhftria,
     &             nrhftria,nrhftria,1,lupri)
         call around('x matrix (triplet)')
         call output(work(kxtrip),1,nrhftria,1,nrhftria,
     &             nrhftria,nrhftria,1,lupri)
      end if

c     ------------------------------------------------- 
c     loop over pairs and make S*R for s = 0,1
c     -------------------------------------------------
      call dzero(work(ksing),n2)
      call dzero(work(ktrip),n2)

      if (ianr12.eq.2) then
c       read r^ab_kl integrals
        lunit = -1
        call gpopen(lunit,fr12r12,'unknown',' ','unformatted',
     &              idum,.false.)
        read(lunit)(work(krabkl-1+i),i=1,nt2r12(1))
        call gpclose(lunit,'KEEP')
c       read R^mn_ab trial vector (stored as triangular matrix)
        call cc_rvec(lufc2,fc2am,nt2am(isymr),nt2am(isymr),
     &               ifile,work(krmnab))
        call cclr_diascl(work(krmnab),two,isymr)
c       calculate \sum_ab r^ab_kl * R^mn_ab
        call dzero(work(ksr),ntr12sq(isymr))
        call cc_r12mi2(work(ksr),work(krabkl),work(krmnab),1,isymr,one,
     &                 work(kend1),lwrk1)
c       get singlet and triplet sr contribution 
c       and write it on work(ksing) and work(ktrip)
        lprojv = .true.
        call cc_r12vunpack(work(ksr),isymr,work(ksing),work(ktrip),
     &                     lprojv,nrhfb,nrhf)
      end if

      ijrsvc = krsing
      ijrtvc = krtrip
      ijsvc = ksing
      ijtvc = ktrip
      do ij = 1, nrhftria
          call cc_r12metric0(ij,nrhftria,
     &                   work(ijsvc),work(ijtvc),
     &                   work(ijrsvc),work(ijrtvc),
     &                   work(kxsing),work(kxtrip))
        ijrsvc = ijrsvc + nrhftria
        ijrtvc = ijrtvc + nrhftria
        ijsvc = ijsvc + nrhftria
        ijtvc = ijtvc + nrhftria
      end do
c     -----------------------
c     print and pack results 
c     -----------------------
      if (locdbg) then
         call around('r12 S*R (singlet)')
         call output(work(ksing),1,nrhftria,1,nrhftria,
     *               nrhftria,nrhftria,1,lupri)
         call around('r12 S*R (triplet)')
         call output(work(ktrip),1,nrhftria,1,nrhftria,
     *               nrhftria,nrhftria,1,lupri)
      end if
      if ( locdbg ) then
         write(lupri,*) 'cc_r12metric> norm^2(S * R12 singlet):',
     &   ddot(n2,work(ksing),1,work(ksing),1)
         write(lupri,*) 'cc_r12metric> norm^2(S * R12 triplet):',
     &   ddot(n2,work(ktrip),1,work(ktrip),1)
      end if
 
c     ------------------
c     write S*R on file
c     ------------------
        call ccr12pck(work(kpck),isymr,work(ksing),work(ktrip),nrhfb,
     &                nrhf,nkilj)
        call cclr_diasclr12(work(kpck),basscl1,isymr)
        if (locdbg) then
          call around('R12 S*R after packing')
          call cc_prpr12(work(kpck),isymr,1,.false.)
        end if
        call cc_wvec(lufs12,fs12am,ntr12am(isymr),ntr12am(isymr),
     *               ifile,work(kpck))
c
        if (ianr12.eq.2) then
c         get \sum_kl r^kl_ab * R^ij_kl + R^ij_ab and put it on file
          call cc_r12getsr22(fc2am,lufc2,fc12am,lufc12,fs2am,
     &                       lufs2,isymr,ifile,work(kend1),lwrk1)
        end if
c
      if (locdbg.and.ianr12.eq.2) then
        call cc_rvec(lufs12,fs12am,ntr12am(isymr),ntr12am(isymr),ifile,
     &                work(kpck))
        write(lupri,*)'FS12AM'
        call cc_prpr12(work(kpck),isymr,1,.false.)
      end if

      if (locdbg) write(lupri,*) 'Leaving CC_R12METRIC'
      call qexit('cc_r12metric')
      end
*=====================================================================* 
      subroutine cc_r12metric0(ij,n,srsing,srtrip,
     &                         rsing,rtrip,xsing,xtrip) 
*---------------------------------------------------------------------*
c     purpose: make the CC2-R12 S*R part for s = 0,1  
c
c     H. Fliegl, C. Haettig summer 2003 
c---------------------------------------------------------------------
      implicit none
#include <priunit.h>
      integer n,ij,kl,mn
#if defined (SYS_CRAY)
      real rsing(n),rtrip(n),xsing(n,n),xtrip(n,n),srsing(n),srtrip(n) 
     &                 
#else
      double precision rsing(n),rtrip(n),xsing(n,n),
     & xtrip(n,n),srsing(n),srtrip(n) 
#endif
      call qenter('cc_r12metric0')

      do kl = 1, n 
         do mn = 1, n
            srsing(kl) = srsing(kl)+ xsing(kl,mn)*rsing(mn)
            srtrip(kl) = srtrip(kl)+ xtrip(kl,mn)*rtrip(mn) 
         end do
      end do

      call qexit('cc_r12metric0')
      end
*=====================================================================*
      subroutine cc_r12metric(isymr,basscl1,basscl2,work,lwork,
     &                        fc2am,lufc2,fc12am,lufc12,fs12am,lufs12,
     &                        fs2am,lufs2,ifile,lnoread,vec)
*---------------------------------------------------------------------*
c     purpose: make the CC-R12 S*R contibution for excitation
c              energies
c              alternative subroutine for cc_r12metric_old without 
c              using singlet/triplet matrices in matrix 
c              multiplication 
c
c     Christian Neiss,  Feb. 2005,  based on old cc_r12metric
c
c     H. Fliegl, spring 2005: nondiagonal elements added for ansatz 2
c     --> S*R written on fs2am
c
c     C. Neiss, 05.04.2005: adapted for left hand transformation
c     basscl1 = brascl for right hand transf. 
c             = 0.5*ketscl for left hand transf.
c     basscl2 = ketscl for right hand transf. 
c             = 2*brascl for left hand transf.
c
c     C. Neiss, 01.07.2004: possibility to get input vector via
c     call-list, not via file: 
c     vec      input vector
c     dimension(vec): nt1am(isymr)+nt2am(isymr)+ntr12am(isymr)
c---------------------------------------------------------------------
      implicit none
#include <ccsdsym.h>
#include <priunit.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <r12int.h>
#include <ccr12int.h>

      logical locdbg
      parameter (locdbg=.FALSE.)

      logical ldum,lnoread
      character*(*) fc12am,fs12am,fc2am,fs2am
      integer lwork,isymr,kend1,lwrk1
      integer lufc12,lufs12,lufc2,lufs2
      integer kxint,ij,kpck,kxintsq,kr12sq,kres,krabkl,krmnab
      integer idum,lunit,ifile,ian,iopt
#if defined (SYS_CRAY)
      real work(*), ddot, one, half, two, basscl1, basscl2, vec(*)
#else
      double precision work(*), ddot, one, half, two, basscl1, basscl2,
     &                 vec(*)
#endif
      parameter (one = 1.0D0, half = 0.5D0, two = 2.0D0)

      call qenter('cc_r12metric')
      if (locdbg) then
        write(lupri,*) 'Entered CC_R12METRIC'
      end if

      kres   = 1 
      kxint  = kres + ntr12sq(isymr)
      kxintsq = kxint + nr12r12p(1)
      kr12sq = kxintsq + nr12r12sq(1)
      kpck   = kr12sq + ntr12sq(isymr)
      kend1   = kpck + ntr12am(isymr)

      if (ianr12.eq.2) then
        krabkl = kend1
        krmnab = krabkl + nt2r12(1)
        kend1  = krmnab + nt2am(isymr)
      end if 

      lwrk1  = lwork  - kend1

      if (lwrk1 .lt. 0) then
         write(lupri,*) 'lwork, lwrk1: ',lwork, lwrk1 
         call quit('Insufficient work space in CC_R12METRIC')
      end if

c     --------------------------------------------------------
c     read R12 overlap matrices and pack to symmetrized square
c     -------------------------------------------------------- 
      lunit = -1
      call gpopen(lunit,fccr12x,'unknown',' ','unformatted',
     &                   idum,ldum)
 9999 read(lunit) ian
      read(lunit) (work(kxint+i), i=0, nr12r12p(1)-1 )
      if (ian.ne.ianr12) goto 9999
      call gpclose(lunit,'KEEP')
      iopt = 2
      call ccr12unpck2(work(kxint),1,work(kxintsq),'T',iopt)

c     -------------------------------------------------
c     read trial vector and pack to symmetrized square
c     -------------------------------------------------
      if (.not. lnoread) then
        call cc_rvec(lufc12,fc12am,ntr12am(isymr),ntr12am(isymr),
     *               ifile,work(kpck))
      else
        call dcopy(ntr12am(isymr),vec(nt1am(isymr)+nt2am(isymr)+1),
     &             1,work(kpck),1)
      end if
      call cclr_diasclr12(work(kpck),basscl2,isymr)
      iopt = 1
      call ccr12unpck2(work(kpck),isymr,work(kr12sq),'T',iopt)

      if (locdbg) then
         write(lupri,*) 'fc12am,lufc12,ifile:',fc12am,lufc12,ifile
         write(lupri,*) 'cc_r12metric> norm^2(R12 packed triangle):',
     &   ddot(ntr12am(isymr),work(kpck),1,work(kpck),1)
         write(lupri,*) 'cc_r12metric> norm^2(R12 packed square):',
     &   ddot(ntr12sq(isymr),work(kr12sq),1,work(kr12sq),1)
      end if
      if (locdbg) then
         call around('R12 trial amplitudes')
         call cc_prsqr12(work(kr12sq),isymr,'T',1,.false.)
      end if

c     ------------- 
c     calculate S*R 
c     -------------
      call dzero(work(kres),ntr12sq(isymr))

      ! conv. doubles contributions for Ansatz 2:
      if (ianr12.eq.2) then
c       read r^ab_kl integrals
        lunit = -1
        call gpopen(lunit,fr12r12,'unknown',' ','unformatted',
     &              idum,.false.)
        read(lunit)(work(krabkl-1+i),i=1,nt2r12(1))
        call gpclose(lunit,'KEEP')
c       read R^mn_ab trial vector (stored as triangular matrix)
        if (.not. lnoread) then
          call cc_rvec(lufc2,fc2am,nt2am(isymr),nt2am(isymr),
     &                 ifile,work(krmnab))
        else
          call dcopy(nt2am(isymr),vec(nt1am(isymr)+1),1,
     &               work(krmnab),1)
        end if
        call cclr_diascl(work(krmnab),two,isymr)
c       calculate \sum_ab r^ab_kl * R^mn_ab
        call cc_r12mi2(work(kres),work(krabkl),work(krmnab),1,isymr,one,
     &                 work(kend1),lwrk1)
c       apply projection operator 0.5*P_kl^ij: 
        call cc_r12pklij(work(kres),isymr,'T',work(kend1),lwrk1)
        call dscal(ntr12sq(isymr),half,work(kres),1)
      end if

      ! add the r12 doubles contribution:
      call cc_r12xi2b(work(kres),'T',work(kxintsq),1,'T',
     &                work(kr12sq),isymr,'T',one)

c     -----------------------
c     print and pack results 
c     -----------------------
      if (locdbg) then
         call around('R12 S*R before packing')
         call cc_prsqr12(work(kres),isymr,'T',1,.false.) 
      end if
      if ( debug ) then
         write(lupri,*) 'cc_r12metric> norm^2(S * R12):',
     &   ddot(ntr12am(isymr),work(kres),1,work(kres),1)
      end if
      iopt = 1
      call ccr12pck2(work(kpck),isymr,.false.,work(kres),'T',iopt)
      call cclr_diasclr12(work(kpck),basscl1,isymr)
      if (locdbg) then
        call around('R12 S*R after packing')
        call cc_prpr12(work(kpck),isymr,1,.false.)
      end if
 
c     ------------------
c     write S*R on file
c     ------------------
      call cc_wvec(lufs12,fs12am,ntr12am(isymr),ntr12am(isymr),
     *             ifile,work(kpck))
c
        if (ianr12.eq.2) then
c         get \sum_kl r^kl_ab * R^ij_kl + R^ij_ab and put it on file
          call cc_r12getsr22(fc2am,lufc2,fc12am,lufc12,fs2am,
     &                       lufs2,isymr,ifile,lnoread,vec,
     &                       work(kend1),lwrk1)
        end if
c
      if (locdbg) write(lupri,*) 'Leaving CC_R12METRIC'
      call qexit('cc_r12metric')
      end
*=====================================================================*

*=====================================================================*
      subroutine cc_r12getsr22(fc2am,lufc2,fc12am,lufc12,
     &                         fs2am,lufs2,isymt,ifile,lnoread,
     &                         vec,work,lwork)
*---------------------------------------------------------------------*
c     purpose: get S_mu2,nu2' * R^ij_kl = \sum_kl r^kl_ab * R^ij_kl
c              needed for ansatz 2
c
c     H. Fliegl, C. Haettig spring 2005
c     C. Neiss, 07.07.2005: modified, see cc_r12metric
*---------------------------------------------------------------------*
      implicit none
#include <priunit.h>
#include <ccsdsym.h>
c
#include <ccsdinp.h>
#include <ccorb.h>
c
#include <r12int.h>
#include <ccr12int.h>
      logical locdbg,lnoread
      parameter (locdbg = .false.)

      character*(*) fc2am,fc12am,fs2am

      integer lufc2,lufc12,lufs2,ifile,isymt,lwork,kend1,lwrk1
      integer krabkl,krmnab,krijkl,lunit,idummy
#if defined (SYS_CRAY)
      real work(*),two,one,half,factor,vec(*)
#else
      double precision work(*),two,one,half,factor,vec(*)
#endif
      parameter (one = 1.0d0, two = 2.0d0, half = 0.5d0)

      call qenter('getsr22')

      krabkl = 1
      krmnab = krabkl + nt2r12(1)
      krijkl = krmnab + nt2am(isymt)
      kend1  = krijkl + ntr12am(isymt)
      lwrk1  = lwork - kend1
      if (lwrk1.lt.0) then
        call quit('isufficient work space in getrs22')
      end if

c     read r^ab_kl integrals
      lunit = -1
      call gpopen(lunit,fr12r12,'unknown',' ','unformatted',
     &            idummy,.false.)
      read(lunit)(work(krabkl-1+i),i=1,nt2r12(1))
      call gpclose(lunit,'KEEP')

      if (.not. lnoread) then
c       read R^mn_ab trial vector (stored as triangular matrix)
        call cc_rvec(lufc2,fc2am,nt2am(isymt),nt2am(isymt),
     &               ifile,work(krmnab))
c       read R^ij_kl
        call cc_rvec(lufc12,fc12am,ntr12am(isymt),ntr12am(isymt),
     *               ifile,work(krijkl))
      else 
        call dcopy(nt2am(isymt),vec(nt1am(isymt)+1),1,
     &             work(krmnab),1)
        call dcopy(ntr12am(isymt),vec(nt1am(isymt)+nt2am(isymt)+1),
     &             1,work(krijkl),1)
      end if
      call cclr_diascl(work(krmnab),two,isymt)
      call cclr_diasclr12(work(krijkl),ketscl,isymt)

c     get \sum_kl r^kl_ab * R^ij_kl and add it to work(krmnab)
      factor = one
      call cc_r12mi22(work(krmnab),work(krabkl),work(krijkl),1,isymt,
     &                factor,work(kend1),lwrk1)
c     call dscal(nt2am(isymt),half,work(krmnab),isymt)
      call cclr_diascl(work(krmnab),half,isymt)
c     write result on file fs2am
      call cc_wvec(lufs2,fs2am,nt2am(isymt),nt2am(isymt),ifile,
     &             work(krmnab)) 
      if (locdbg) then
        call cc_rvec(lufs2,fs2am,nt2am(isymt),nt2am(isymt),ifile,
     &               work(krmnab))
        write(lupri,*)'FS2AM'
        write(lupri,*)'nt2am(isymt):', nt2am(isymt)
        write(lupri,*)'ntr12am(isymt):',ntr12am(isymt)
        call outpak(work(krmnab),nt2am(isymt),1,lupri)
      end if

      call qexit('getsr22')
      end
*=====================================================================* 

*=====================================================================*
      subroutine cc_r12buildbmat(bmatsq,isymkl,isymmn,
     &                           isymij,idxij,bmatsq0,
     &                           xintsq,epskl,inmatkl,epsij,inmatij,
     &                           factor)
*---------------------------------------------------------------------*
c     purpose: Build up R12 B-matrix B_(kl,mn)^(ij) for 
c              fixed index pair (ij)
c
c     bmatsq   B-matrix block of dimension nmatkl(isymkl)*nmatkl(isymmn)
c              analog: bmatsq0, xintsq
c     isymij   symmetry of pair (ij)
c     idxij    sorted pair index (ij)
c
c     Christian Neiss,  April 2005, restructured September 2005
c---------------------------------------------------------------------
      implicit none
#include <ccsdsym.h>
#include <priunit.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <r12int.h>
#include <ccr12int.h>

      logical locdbg
      parameter (locdbg=.false.)

      integer idxij,idxkl,idxmn,idxklmn
      integer isymij,isymkl,isymmn
      integer inmatkl(8),inmatij(8)

#if defined (SYS_CRAY)
      real bmatsq(*),bmatsq0(*),xintsq(*),epskl(*),epsij(*),factor,two
#else
      double precision bmatsq(*),bmatsq0(*),xintsq(*),epskl(*),
     &                 epsij(*),factor,two
#endif
      parameter (two = 2.0D0)

      call qenter('cc_r12buildbmat')
      if (locdbg) then
        write(lupri,*)'Entered CC_R12BUILDBMAT'
      end if

c     --------------
c     build B-matrix
c     --------------
      do idxmn = 1, nmatkl(isymmn)
        do idxkl = 1, nmatkl(isymkl)
          idxklmn =  nmatkl(isymkl)*(idxmn-1) + idxkl
c         write(lupri,*) 'idxklmn: ',idxklmn
          bmatsq(idxklmn) = bmatsq0(idxklmn) +
     &                   factor*(epskl(inmatkl(isymkl)+idxkl)+
     &                           epskl(inmatkl(isymmn)+idxmn)-
     &                   two*epsij(inmatij(isymij)+idxij))*
     &                   xintsq(idxklmn)
        end do
      end do

      if (locdbg) then
        call around('B-Matrix in CC_R12BUILDBMAT')
        write(lupri,*) 'idxij, epsilon(ij): ',idxij,
     &                 epsij(inmatkl(isymij)+idxij)
        call output(bmatsq,1,nmatkl(isymkl),1,nmatkl(isymmn),
     &              nmatkl(isymkl),nmatkl(isymmn),1,lupri)
      end if

      if (locdbg) then
        write(lupri,*)'Leaving CC_R12BUILDBMAT'
      end if

      call qexit('cc_r12buildbmat')
      return
      end 
*=====================================================================*

*=====================================================================*
      subroutine cc_r12nxtam(omeg12,isymom,tamp12,lcceq,
     &                       er12,work,lwork)
c----------------------------------------------------------------------
c     purpose: get new start R12  amplitudes for each iteration
c              by inversion of the B matrix.
c              substitute for old cc_r12nxtam: working without 
c              singlet/triplet format, similiar to conventional
c              part of vector function
c
c     C. Neiss,  december 2005
c----------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
#include <ccr12int.h>
#include <ccsdinp.h>
#include <iratdef.h>
#include <second.h>

      logical locdbg,ldum,lcceq,etest 
      parameter (locdbg = .false., etest = .false.)

      integer kbmatsq,komegsq,kend1,lwrk1,kend2,lwrk2
      integer kxintsq,kscr,keij,kekl,kevl,kpvt,kzvec,koffom,koffb
      integer lwork,ian,iap,iopt,isym,lunit,idum
      integer isymij,isymmn,isymkl,idxij,isymom
      integer nij,nkl,inmatij(8),inmatkl(8),norbtsx
#if defined (SYS_CRAY)
      real tamp12(*),omeg12(*),er12,work(*),ddot,one,half,zero,factor,
     &     rcond,thrzero,timnxtam
#else
      double precision tamp12(*),omeg12(*),er12,work(*),ddot,one,half,
     &                 zero,factor,rcond,thrzero,timnxtam
#endif
      character*3 aprox
      character*10 model
      parameter (one = 1.0D0, half = 0.5D0, zero = 0.0D0)
      parameter (thrzero = 1.0D-12)

      call qenter('r12nxtam')
      if (locdbg) then
        write(lupri,*) 'Entered CC_R12NXTAM'
        call flshfo(lupri)
      end if

      timnxtam = second()
  
      kbmatsq = 1
      komegsq = kbmatsq+ nr12r12sq(1)
      kend1   = komegsq+ ntr12sq(isymom)
     
      kscr   = kend1
      kend1  = kscr + max(nr12r12sq(1),ntr12am(1))

      lwrk1  = lwork  - kend1

      if (lwrk1 .lt. 0) then
         call quit('Insufficient work space in cc_r12nxtam')
      end if

c     -------------
c     test symmetry
c     -------------
      if (lcceq .and. (isymom.ne.1)) then
        call quit('Symmetry error in CC_R12NXTAM')
      end if
c
c     ----------------
c     read B matrices
c     ----------------
      lunit = -1
      call gpopen(lunit,fccr12b,'old',' ','unformatted',
     &                   idum,ldum)
 8881 read(lunit) ian, iap, aprox
      read(lunit) (work(kscr+i), i=0, nr12r12p(1)-1 )
      if ((ian.ne.ianr12).or.(iap.ne.iapr12)) goto 8881
      call gpclose(lunit,'KEEP')
      iopt = 2
      call ccr12unpck2(work(kscr),1,work(kbmatsq),'N',iopt)
c     call dscal(nr12r12sq(1),-one,work(kbmatsq),1)      
c 
      if (lcceq) then
c       --------------------
c       read R12 amplitudes
c       --------------------
        iopt =32
        call cc_rdrsp('R0 ',0,1,iopt,model,idum,tamp12)
c
        if (iprint .ge. 5) then
          write(lupri,529) 'Norm^2 of r12am     in NXTAM:',
     *                      ddot(ntr12am(1),tamp12,1,tamp12,1)
        end if
c
c       -------------------------------
c       read vector function omega_R12
c       -------------------------------
        lunit = -1 
        call gpopen(lunit,'CC_OMEGAR12','old',' ','unformatted',
     &              idum,ldum)
        read(lunit) (omeg12(i), i=1, ntr12am(1))
        call gpclose(lunit,'KEEP')  
c 
        write(lupri,529)'Norm^2 of Omegar12  in NXTAM:', 
     &     ddot(ntr12am(1),omeg12,1,omeg12,1)
c
      end if
c
      !repack Omega as Omega_(kl,ij)
      iopt = 1
      call ccr12unpck2(omeg12,isymom,work(komegsq),'N',iopt)

 529  format(7x,a,d24.10)
c
      if (ianr12.eq.1) then
        if (iapr12.le.2 .or. iapr12.eq.4 .or. iapr12.eq.6) then
          factor = zero
        else
          factor = -half
        end if
      else
        if (iapr12.eq.2 .or. iapr12.eq.5) then
          factor = zero
        else
          factor = -half
        end if
      end if
c
c     -----------------------------------------------
c      - read R12 overlap matrix, orbital energies
c      - loop over pairs (ij)
c     -----------------------------------------------
c
      !calculate # of indices (kl) / (ij) and offset over all symmetries:
      nkl = 0
      nij = 0
      do isym = 1, nsym
        inmatkl(isym) = nkl
        inmatij(isym) = nij
        nkl = nkl + nmatkl(isym)
        nij = nij + nmatij(isym)
      end do
c
      kxintsq = kend1
      kekl  = kxintsq + nr12r12sq(1)
      keij  = kekl + nkl
      keij  = kekl + nkl
      kend1 = keij + nij
      lwrk1 = lwork - kend1
      if (lwrk1 .lt. 0) then
        call quit('Insufficient work space in cc_r12nxtam')
      end if
c
      lunit = -1
      call gpopen(lunit,fccr12x,'old',' ','unformatted',
     &            idum,ldum)
 9991 read(lunit) ian
      read(lunit) (work(kscr+i), i=0, nr12r12p(1)-1 )
      if (ian.ne.ianr12) goto 9991
      call gpclose(lunit,'KEEP')
      iopt = 2
      call ccr12unpck2(work(kscr),1,work(kxintsq),'N',iopt)
c
      !read orbital energies
      lunit = -1
      call gpopen(lunit,'SIRIFC','OLD',' ','UNFORMATTED',idum,
     &            .false.)
      rewind lunit
      call mollab('FULLBAS ',lunit,lupri)
      read (lunit) idum,norbtsx
      kevl = kend1
      kend2 = kevl + norbtsx
      lwrk2 = lwork - kend2
      if (lwrk2 .lt. 0) then
        call quit('Insufficient work space in cc_r12nxtam')
      end if
      read (lunit) (work(kevl+i-1), i=1,norbtsx)
      call gpclose(lunit,'KEEP')
c
      !calculate pair orbital energies:
      call cc_r12epair(work(kekl),nkl,work(kevl),inmatkl,imatkl,nrhfb)
      call cc_r12epair(work(keij),nij,work(kevl),inmatij,imatij,nrhf)
c
c     -------------------------------
c     calculate update for amplitudes
c     -------------------------------
c
      do isymij = 1, nsym
        isymmn = muld2h(isymom,isymij)
        isymkl = isymmn
        if (nmatkl(isymkl).gt.0) then
        do idxij = 1, nmatij(isymij)
c         call dzero(work(kscr),nmatkl(isymkl)*nmatkl(isymmn))
          koffb = ir12r12sq(isymkl,isymmn) 
          call cc_r12buildbmat(work(kscr),isymkl,isymmn,
     &                         isymij,idxij,work(kbmatsq+koffb),
     &                         work(kxintsq+koffb),work(kekl),inmatkl,
     &                         work(keij),inmatij,factor)
c
c         invert B-matrix and contract with omega:
c
          kpvt  = kend1
          kzvec = kpvt + nmatkl(isymkl)/irat + 1
          kend2 = kzvec + nmatkl(isymkl)
          lwrk2 = lwork - kend2
          if (lwrk2 .lt. 0) then
            call quit('Insufficient work space in cc_r12nxtam')
          end if
c
          call dgeco(work(kscr),nmatkl(isymkl),nmatkl(isymkl),
     &               work(kpvt),rcond,work(kzvec))
          if (rcond.lt.thrzero) then
            call quit('B-Matrix is singular in CC_R12NXTAM')
          end if
c
          koffom = komegsq + itr12sq(isymkl,isymij) +
     &             nmatkl(isymkl)*(idxij-1)

c         write(lupri,*) 'Omega before dgesl:'
c         call output(work(koffom),1,nmatkl(isymkl),1,1,
c    &                nmatkl(isymkl),1,1,lupri)

          call dgesl(work(kscr),nmatkl(isymkl),nmatkl(isymkl),
     &               work(kpvt),work(koffom),0)

c         write(lupri,*) 'Omega after dgesl:'
c         call output(work(koffom),1,nmatkl(isymkl),1,1,
c    &                nmatkl(isymkl),1,1,lupri)
c
        end do
        end if
      end do
c
      iopt = 1
      call ccr12pck2(omeg12,isymom,.false.,work(komegsq),'N',iopt)
c
      if (lcceq) then
c       ----------------------------------------------------
c       add old amplitudes to update to get new amplitudes
c       ----------------------------------------------------
        call daxpy(ntr12am(isymom),one,tamp12,1,omeg12,1)
c
        if (locdbg) then
          write(lupri,*) 'R12 amplitudes:'
          call cc_prpr12(tamp12,isymom,1,.false.)
          write(lupri,*) 'Updated R12 amplitudes:'
          call cc_prpr12(omeg12,isymom,1,.false.)
        end if
c
        if (etest) then
c         -------------    
c         read V matrix 
c         -------------
          lunit = -1
          call gpopen(lunit,fccr12v,'old',' ','unformatted',
     &                idum,ldum)
 7771     read(lunit) ian
          read(lunit) (work(kscr+i), i=0, ntr12am(1)-1)
          if (ian.ne.ianr12) goto 7771
          call gpclose(lunit,'KEEP')
c
          call cc_r12tcmepk(work(kscr),1,.false.)
          call cclr_diasclr12(work(kscr),half,isymom)
          er12 = 2.0D0*ddot(ntr12am(isymom),omeg12,1,work(kscr),1)
c
          write(lupri,*) 'R12 contribution to energy: ',er12
        end if
c
      end if ! lcceq
c
      if (locdbg) write(lupri,*) 'Leaving CC_R12NXTAM'
c
      timnxtam = second() - timnxtam
c     write(lupri,111) 'CC_R12NXTAM', timnxtam
  111 FORMAT(7x,'Time used in',2x,A12,2x,': ',f10.2,' seconds')
c
      call qexit('r12nxtam')
      end
*=====================================================================*

*=====================================================================*
      subroutine cc_r12epair(evlout,length,evlin,inmat,imat,nrhf1)
c----------------------------------------------------------------------
c     purpose: get pair orbital energies
c     
c     it is assumed that the R12 orbital energies are leading 
c     (incl. frozen orbitals, as for "FULLBAS" in SIRIFC)
c
c     C. Neiss september 2005
c----------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
#include <ccsdinp.h>

      logical locdbg
      parameter (locdbg = .false.)

      integer kend1,lwrk1,length,kscr1
      integer lwork,idummy,icount,isym
      integer isymij,isymi,isymj,mi,mj,idxij 
      integer nrhf1(8),inmat(8),imat(8,8),ioff(8)
#if defined (SYS_CRAY)
      real evlout(length),evlin(*),ei,ej
#else
      double precision evlout(length),evlin(*),ei,ej
#endif
c
      call qenter('cc_r12epair')
      if (locdbg) then
        write(lupri,*) 'Entered CC_R12EPAIR'
        call flshfo(lupri)
      end if
c
      icount = 0
      do isym = 1, nsym
        icount = icount + nrhffr(isym)
        ioff(isym) = icount 
        icount = icount + nrhfb(isym) + norb1(isym) + norb2(isym)
      end do
c
      call dzero(evlout,length)
c
      do isymij = 1, nsym
        do isymi = 1, nsym
          isymj = muld2h(isymi,isymij)
          do i = 1, nrhf1(isymi)
c           mi = iorb(isymi)+i
            mi = ioff(isymi)+i
            ei = evlin(mi)
            do j = 1, nrhf1(isymj)
              idxij = inmat(isymij) + imat(isymi,isymj) + 
     &                nrhf1(isymi)*(j-1)+i
c             mj = iorb(isymj)+j
              mj = ioff(isymj)+j
              ej = evlin(mj)
              evlout(idxij) = ei + ej
            end do
          end do
        end do
      end do 
c
      if (locdbg) then
        call around('CC_R12EPAIR: pair orbital energies')
        do isymij = 1, nsym
          do isymi = 1, nsym
            isymj = muld2h(isymij,isymi)
            write(lupri,*) 'Symmetry block ',isymi, isymj
            call output(evlout(1+inmat(isymij)+imat(isymi,isymj)),
     &                  1,nrhf1(isymi),1,nrhf1(isymj),
     &                  nrhf1(isymi),nrhf1(isymj),1,lupri)
          end do
        end do
      end if
c
      if (locdbg) write(lupri,*) 'Leaving CC_R12EPAIR'

      call qexit('cc_r12epair')
      end
*=====================================================================*

