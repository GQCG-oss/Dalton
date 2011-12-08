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
#include <single.h>

*=====================================================================*
      subroutine CCLHTR_FP(RHO1,XLAMDH,WORK,LWORK,ISYMC,FC12AM,LUFC12,
     &           IFILE)
*---------------------------------------------------------------------*
C     purpose: Calculate Sum{klj} c^{ij}_{kl}*V^{a jt}_{kl} =: rho_Fp
C
C     C. Neiss, C. Hättig, summer 2004
*---------------------------------------------------------------------*
      implicit none
#include <r12int.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccr12int.h>
#include <dummy.h>

      integer LWORK,ISYMC,LUFC12,IFILE,isymak,isymk,isyma,
     &        kvec,kvajkl,kend1,lwrk1,luvajkl
      character*(*) FC12AM
      character     cdummy*(8)
      logical locdbg
      parameter(locdbg=.false.)

#if defined (SYS_CRAY)
      real RHO1(*),XLAMDH(*),WORK(LWORK),ddot
#else
      double precision RHO1(*),XLAMDH(*),WORK(LWORK),ddot
#endif

      call qenter('cclhtr_fp')

C     do isymak = 1, nsym
C        nvajkl(isymak) = 0
C        do isymk = 1, nsym
C           isyma = muld2h(isymak,isymk)
C           nvajkl(isymak) = nvajkl(isymak) + nt1ao(isyma)*nmatkl(isymk)
C        end do
C     end do

      kvec = 1
      kvajkl = kvec + ntr12sq(isymc)
      kend1 = kvajkl + nvajkl(1)
      lwrk1 = lwork - kend1

      if (lwrk1 .lt.0) then
         call quit('Insufficient work space in cclhtr_fp')
      end if

C     -----------------------
C     read V^{alpha jt}_{kl}
C     -----------------------
      luvajkl=-1
      call gpopen(luvajkl,'CCR12VAJTKL','old',' ','unformatted',idummy,
     &            .false.)
      rewind(luvajkl)
      read(luvajkl) (work(kvajkl+i-1), i=1,nvajkl(1))
      call gpclose(luvajkl,'KEEP')
      if(locdbg) then
        write(lupri,*) 'norm^2 V^{alpha j~}_{kl} in cclhtr_fp:',
     &                  ddot(nvajkl(1),work(kvajkl),1,work(kvajkl),1)
      end if

c     ---------------------
c     read r12 trial vector
c     ---------------------
      call cc_r12getct(work(kvec),isymc,1,2.0D0*brascl,.false.,'T',
     &                 lufc12,fc12am,ifile,cdummy,idummy,work(kend1),
     &                 lwrk1)

c     --------------
c     get V^aj_kl*C
c     --------------
      call ccrhs_gp0(rho1,isymc,work(kvajkl),1,xlamdh,1,
     &               work(kvec),isymc,.false.,dummy,locdbg,
     &               work(kend1),lwrk1)

      if (locdbg) then
        write(lupri,*)'in cclhtr_fp: rho1 on exit '
        call cc_prp(rho1,dummy,1,1,0)
      end if

      call qexit('cclhtr_fp')
      end


*=====================================================================*
      subroutine CCLHTR_GP(CTR1,ISYCTR,XLAMDP,ISYLAM,RHOR12,ISYRES,
     &                    WORK,LWORK)
*---------------------------------------------------------------------*
C     purpose: Calculate Sum_{ai mu} [C_{ai}*Lambda^p_{a mu}*
C                (2V(mu j,kl)-V(mu j,lk)] and add to RHOR12 
C
C     C. Neiss, C. Hättig, summer 2004
*---------------------------------------------------------------------*
      implicit none
#include <r12int.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccr12int.h>
#include <ccsdinp.h>

      logical locdbg
      parameter (locdbg = .false.)

      integer isymmu,isymi,isyma,isymmua,kofflp,kofflpb,koffctr,
     &        ntotmu,ntota,isyctr,isyres,isylam
      integer klamdpb,kend1,lwrk1,lwork

#if defined (SYS_CRAY)
      real ctr1(*),XLAMDP(*),rhor12(*),WORK(LWORK),zero,one,ddot
#else
      double precision ctr1(*),XLAMDP(*),rhor12(*),WORK(LWORK),zero,one,
     &                 ddot
#endif

      parameter (one = 1.0d0, zero = 0.0d0)

      call qenter('cclhtr_gp')
      if (locdbg) then
        write(lupri,*) 'Entered CCLHTR_GP'
      end if

      klamdpb = 1
      kend1 = klamdpb + nglmrh(isyres)
      lwrk1 = lwork - kend1
      if (lwrk1 .lt. 0) then
        call quit('Insufficient work space in cclhtr_gp')
      end if

C     Test if symmetry is consistent
      if (isyres .ne. muld2h(isyctr,isylam)) then
        call quit('Symmetry mismatch in CCLHTR_GP!')
      end if

C     ---------------------------------
C     Calculate Lambda^p*C = Lambda^pb
C     ---------------------------------

      do isymi = 1,nsym
         isyma = muld2h(isymi,isyctr)
         isymmua = isylam
         isymmu = muld2h(isyma,isymmua)
         kofflp  = iglmvi(isymmu,isyma) + 1
         kofflpb = iglmrh(isymmu,isymi) + 1
         koffctr = it1am(isyma,isymi) + 1
         ntotmu = max(1,nbas(isymmu))
         ntota = max(1,nvir(isyma))
         CALL DGEMM('N','N',nbas(isymmu),nrhf(isymi),nvir(isyma),
     &        one,xlamdp(kofflp),ntotmu,ctr1(koffctr),ntota,
     &        zero,work(klamdpb-1+kofflpb),ntotmu)
      end do
      if (locdbg) then
         write(lupri,*) 'isyctr, isylam, isyres: ',isyctr,isylam,isyres
         write(lupri,*) 'norm^2(ctr1):',
     &    ddot(nt1am(isyctr),ctr1,1,ctr1,1)
         write(lupri,*) 'norm^2(lambda^p*C):',
     &    ddot(nglmrh(isyres),work(klamdpb),1,work(klamdpb),1)
      end if

c     ----------------------------
c     make the rest 
c     ----------------------------
      call cc_r12eta0(rhor12,work(klamdpb),isyres,work(kend1),lwrk1)

      if (locdbg) then
         call around('rhor12 in CCLHTR_GP')
         call cc_prsqr12(rhor12,isyres,'T',1,.false.)
         write(lupri,*) 'norm^2(rhor12):',
     &     ddot(ntr12sq(isyres),rhor12,1,rhor12,1)

         write(lupri,*) 'Leaving CCLHTR_GP'
      end if

      call qexit('cclhtr_gp')
      return
      end
*=====================================================================*

*=====================================================================*
      subroutine CC_R12CV(XGAMMA,C12AMP,ISYMC,VINT,ISYMV,IOPTRES,
     &                    WORK,LWORK)
*---------------------------------------------------------------------*
C     purpose: Calculate Sum_{mn} (V\dagger)^mn_ij * c^kl_mn
C              add to XGAMMA
C
C     IOPTRES = 1: Dim. of XGAMMA = NGAMMA(ISYGAM)
C             = 2: Dim. of XGAMMA = N3ORHF(ISYGAM)
C
C     C. Neiss october/november 2005
*---------------------------------------------------------------------*
      implicit none
#include <r12int.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccr12int.h>
#include <ccsdinp.h>

      integer isygam,isymc,isymv,isymij,isymmn,isymkl
      integer lwork,kend1,lwrk1,koff1,koff2,koff3,kscr,iopt,
     &        idxik,idxjl,kcamp,kvint,ioptres,isym1,isym2

#if defined (SYS_CRAY)
      real XGAMMA(*),C12AMP(*),VINT(*),WORK(LWORK),zero,one
#else
      double precision XGAMMA(*),C12AMP(*),VINT(*),WORK(LWORK),zero,one
#endif
      logical locdbg
      parameter (locdbg = .FALSE.)
      parameter (one = 1.0d0, zero = 0.0d0)
C
      call qenter('cc_r12cv')
      if (locdbg) then
        write(lupri,*) 'Entered CC_R12CV'
      end if
C
      kend1 = 1
C
      isygam = muld2h(isymc,isymv)
C 
      if (locdbg) then
        write(Lupri,*) 'Gamma on entry:'
        if (ioptres.eq.1) then
          call cc_prpr12(xgamma,isygam,1,.FALSE.)
        else if (ioptres.eq.2) then
          do isym2 = 1, nsym
            isym1 = muld2h(isygam,isym2)
            write(lupri,*)'Symmetry block ', isym1, isym2
            call output(xgamma(1+i3orhf(isym1,isym2)),
     &                  1,nmaijk(isym1),1,nrhf(isym2),
     &                  nmaijk(isym1),nrhf(isym2),1,lupri)
          end do
        end if
        write(lupri,*) 
        write(Lupri,*) 'R12 amplitudes on entry:'
        call cc_prpr12(C12AMP,isymc,1,.FALSE.)
        write(lupri,*)
        write(Lupri,*) 'V-intermediate on entry:'
        call cc_prpr12(VINT,isymv,1,.FALSE.)
        write(lupri,*)
      end if
C
      !unpack amplitudes and V-interm. to square:
      kcamp = kend1
      kvint = kcamp + max(ntr12sq(isymc),ngamma(isygam))
      kscr  = kvint + ntr12sq(isymv)
      kend1 = kscr  + ngamsq(isygam)
      lwrk1 = lwork - kend1
      if (lwrk1.lt.0) then
        call quit('insufficient work space in CC_R12CV')
      end if
C
      call cclr_diasclr12(c12amp,ketscl,isymc)
      iopt = 1
      call ccr12unpck2(c12amp,isymc,work(kcamp),'T',iopt)
      call ccr12unpck2(vint,isymv,work(kvint),'T',iopt)
      if (locdbg) then
        write(lupri,*) 'R12 amplitudes (squared):'
        call cc_prsqr12(work(kcamp),isymc,'T',1,.FALSE.)
        write(lupri,*) 'V-intermediate (squared):'
        call cc_prsqr12(work(kvint),isymv,'T',1,.FALSE.)
      end if
C
      !contract :
      do isymij = 1, nsym
        isymmn = muld2h(isymij,isymv)
        isymkl = muld2h(isymmn,isymc)
        koff1 = kcamp + itr12sqt(isymkl,isymmn)
        koff2 = kvint + itr12sqt(isymij,isymmn)
        koff3 = kscr + igamsq(isymij,isymkl)
        call dgemm('N','T',nmatij(isymij),nmatij(isymkl),
     &             nmatkl(isymmn),one,work(koff2),max(1,nmatij(isymij)),
     &             work(koff1),max(1,nmatij(isymkl)),zero,
     &             work(koff3),max(1,nmatij(isymij)))
      end do
C
      if (ioptres.eq.1) then
        !pack to triangle and add to Gamma:
        iopt = 0
        call ccr12pck2(work(kcamp),isygam,.false.,work(kscr),'T',iopt)
        call daxpy(ngamma(isygam),one,work(kcamp),1,xgamma,1)
      else if (ioptres.eq.2) then
        iopt = 0
        call cc_sort4o2(xgamma,isygam,work(kscr),iopt,.TRUE.)
      end if
C
      if (locdbg) then
        write(Lupri,*) 'Gamma on exit:'
        write(lupri,*) 'before reorder:'
        call cc_prsqr12(work(kscr),isygam,'T',1,.FALSE.)
        write(lupri,*) 'after reorder:'
        if (ioptres.eq.1) then
          call cc_prpr12(xgamma,isygam,1,.FALSE.)
        else if (ioptres.eq.2) then
          do isym2 = 1, nsym
            isym1 = muld2h(isygam,isym2)
            write(lupri,*)'Symmetry block ', isym1, isym2
            call output(xgamma(1+i3orhf(isym1,isym2)),
     &                  1,nmaijk(isym1),1,nrhf(isym2),
     &                  nmaijk(isym1),nrhf(isym2),1,lupri)
          end do
        end if
      end if
C
      if (locdbg) then
        write(lupri,*) 'Leaving CC_R12CV'
      end if      
      call qexit('cc_r12cv')
      return
      end
*=====================================================================*

*=====================================================================*
      subroutine CCLHTR_AP(RHOR12,VINT,XMINT,ISYMTR,WORK,LWORK)
*---------------------------------------------------------------------*
C     purpose: Calculate Sum_{mn} M^ij_mn * (V\dagger)^kl_mn
C              add to R12-part of rho
C              i,j,m,n: occ. indices
C              k,l    : r12 indices
C
C     C. Neiss october 2005
*---------------------------------------------------------------------*
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
C
C     INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      integer lwork,kend1,lwrk1,kvint,kmint,iopt,isymtr
      integer isymij,isymmn,isymkl,koffv,koffm,koffr,isym1,isym2
#if defined (SYS_CRAY)
      real RHOR12(*),VINT(*),XMINT(*),WORK(LWORK),zero,one,two
#else
      double precision RHOR12(*),VINT(*),XMINT(*),WORK(LWORK),zero,one,
     &                 two
#endif
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      logical locdbg
      parameter (locdbg = .FALSE.)
C
      call qenter('cclhtr_ap')
      if (locdbg) then
        write(lupri,*) 'Entered CCLHTR_AP'
      end if

      kvint = 1
      kmint = kvint + ntr12sq(1)
      kend1 = kmint + ngamsq(isymtr)
      lwrk1 = lwork - kend1
      if (lwrk1.lt.0) call quit('Insufficient memory in CCLHTR_AP')

      !pack V-intermediate to square format:
      iopt = 1
      call ccr12unpck2(vint,1,work(kvint),'T',iopt)
      if (locdbg) then
        write(lupri,*) 'V-intermediate in CCLHTR_AP:'
        call cc_prsqr12(work(kvint),1,'T',1,.false.)
      end if

      !reorder M-intermediate to square format:
      iopt = 1
      call cc_sort4o2(xmint,isymtr,work(kmint),iopt,.FALSE.)
      if (locdbg) then
        write(lupri,*) 'M-intermediate in CCLHTR_AP before reorder:'
        do isym2 = 1, nsym
          isym1 = muld2h(isymtr,isym2)
          write(lupri,*)'Symmetry block ', isym1, isym2
          call output(xmint,1,nmaijk(isym1),1,nrhf(isym2),
     &                nmaijk(isym1),nrhf(isym2),1,lupri)
        end do
C
        write(lupri,*) 'M-intermediate in CCLHTR_AP after reorder:'
        call cc_prsqr12(work(kmint),isymtr,'T',1,.false.)
      end if

      !transform:
      do isymkl = 1, nsym
        isymmn = isymkl
        isymij = muld2h(isymmn,isymtr)
        koffv = kvint + itr12sqt(isymmn,isymkl)
        koffm = kmint + igamsq(isymij,isymmn)
        koffr = 1 + itr12sqt(isymij,isymkl)
        call dgemm('N','N',nmatij(isymij),nmatkl(isymkl),nmatij(isymmn),
     &             one,work(koffm),max(1,nmatij(isymij)),
     &             work(koffv),max(1,nmatij(isymmn)),one,
     &             rhor12(koffr),max(1,nmatij(isymij)))
      end do
C
      if (locdbg) then
        write(lupri,*) 'Leaving CCLHTR_AP'
      end if
      call qexit('cclhtr_ap')
      RETURN
      END
*=====================================================================*

*=====================================================================*
      subroutine CCLHTR_EP(RHOR12,VINT,XMAT,ISYMTR,WORK,LWORK)
*---------------------------------------------------------------------*
C     purpose: Calculate 2C-E [-2 * Sum_{n} (V\dagger)^kl_in * X_nj]
C              add to RHOR12
C              k,l:   r12 indices
C              i,j,n: occ. indices
C
C     C. Neiss october 2005
*---------------------------------------------------------------------*
      implicit none
#include <r12int.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccr12int.h>
#include <ccsdinp.h>

      integer isymkl,isymin,isymij,isymnj,isymi,isymj,isymn,isymk,isyml
      integer idxkl,koffx,kvinkl,kvijkl,isymtr 
      integer lwork,kend1,lwrk1,iopt,kvint0

#if defined (SYS_CRAY)
      real RHOR12(*),VINT(*),XMAT(*),WORK(LWORK),zero,one,two
#else
      double precision RHOR12(*),VINT(*),XMAT(*),WORK(LWORK),zero,one,
     &                 two
#endif
      logical locdbg
      parameter (locdbg = .FALSE.)
      parameter (one = 1.0d0, zero = 0.0d0, two = 2.0d0)

      call qenter('cclhtr_ep')
      if (locdbg) then
        write(lupri,*) 'Entered CCLHTR_EP'
      end if

      kvint0 = 1
      kend1  = kvint0 + ntr12sq(1)
      lwrk1  = lwork - kend1
      if (lwrk1.lt.0) call quit('Insufficient memory in CCLHTR_EP')

      !pack V-intermediate to square format and make 2C-E combination:
      iopt = 1
      call ccr12unpck2(vint,1,work(kvint0),'T',iopt)
      call cc_r12tcmesq(work(kvint0),1,'T',.FALSE.)

      !transform with XMAT and add to rho:
      do isymkl = 1, nsym
        isymin = isymkl
        do isymn = 1, nsym
          isymi = muld2h(isymin,isymn)
          isymj = muld2h(isymtr,isymn)
          isymnj = isymtr
          isymij = muld2h(isymi,isymj)
          do isymk = 1, nsym
            isyml = muld2h(isymkl,isymk)
            do k = 1, nrhfb(isymk)
              do l = 1, nrhfb(isyml)
                idxkl = imatkl(isymk,isyml)+nrhfb(isymk)*(l-1)+k
                kvinkl = kvint0 + itr12sqt(isymin,isymkl) + 
     &                   nmatij(isymin)*(idxkl-1)+
     &                   imatij(isymi,isymn) 
                kvijkl = 1 + itr12sqt(isymij,isymkl) +
     &                   nmatij(isymij)*(idxkl-1)+
     &                   imatij(isymi,isymj)
                koffx  = 1 + imatij(isymn,isymj) 
                call dgemm('N','N',nrhf(isymi),nrhf(isymj),nrhf(isymn),
     &                     -two,work(kvinkl),max(1,nrhf(isymi)),
     &                     xmat(koffx),max(1,nrhf(isymn)),
     &                     one,rhor12(kvijkl),max(1,nrhf(isymi)))
              end do
            end do
          end do
        end do
      end do 

      call qexit('cclhtr_ep')
      if (locdbg) then
        write(lupri,*) 'Leaving CCLHTR_EP'
      end if

      return
      end
*======================================================================*





