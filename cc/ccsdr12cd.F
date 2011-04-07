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
      subroutine ccsdr12cd(ccsdr12,
     &                     omegar12,isymom,t2am,isymt2,isymim,
     &                     fniadj,luiadj,fnijda,luijda,it2del,
     &                     xlamdah,isymlam,
     &                     cmox,ilamdx,
     &                     work,lwork)
*----------------------------------------------------------------------*
*  Purpose: compute the C and D contributions to the R12 result vector
*
* C. Haettig, C. Neiss, spring 2006
*----------------------------------------------------------------------*
      implicit none
#include <priunit.h>
#include <dummy.h>     
#include <ccorb.h>     
#include <ccsdsym.h>     
#include <r12int.h>     
#include <ccr12int.h>     

      logical locdbg
      parameter (locdbg = .FALSE.)
      integer isym0
      parameter ( isym0 = 1)
      real*8  zero,one
      parameter ( zero=0.0d0, one=1.0d0 )

* input:
      logical ccsdr12
      character fnijda*8, fniadj*8
      integer isymlam, ilamdx(8,8), luijda, luiadj, isymom, isymt2,
     &        it2del(*), lwork, isymim
      real*8  xlamdah(*),omegar12(*),t2am(*),cmox(*),work(*)

* local:
      logical lcbs, ldum
      character cdummy*8
      integer isym, nbas2(8), iviraop(8,8),
     &        imatpb(8,8), nviraop(8), nmatpb(8), 
     &        kend1, lwrk1, kemat1, kfckvao, kfckvaop, kfgdp, kend2,
     &        lwrk2, lunit, isymp, isydel, isygam, isymc, kscr1, kscr2,
     &        kend3, lwrk3, koffv, len1, len2, koffl, koffgdp, nbasg,
     &        nvirc, koffe1, koffcx, nbasd, norbp, kcdint, iopt, 
     &        komegpk, isym1, isym2, igdp(8,8), ngdp(8), iopttcme,
     &        ioptr12, idum, mtotbas
      real*8  fac, ddot

*----------------------------------------------------------------------*
* precompute non-standard symmetry offsets and dimensions:
*----------------------------------------------------------------------*
      mtotbas = 0
      do  isym = 1, nsym
        ! total number of basis functions in primary + aux. basis
        nbas2(isym) = mbas1(isym) + mbas2(isym)
        mtotbas = mtotbas + nbas2(isym)
      end do

      do  isym = 1, nsym
        nviraop(isym) = 0
        ngdp(isym)    = 0
        nmatpb(isym)  = 0
        do isym2 = 1, nsym
           isym1 = muld2h(isym2,isym)
           iviraop(isym1,isym2) = nviraop(isym)
           nviraop(isym) = nviraop(isym) + nvir(isym1)*mbas2(isym2)
           imatpb(isym1,isym2) = nmatpb(isym)
           nmatpb(isym) = nmatpb(isym) + norb2(isym1)*nvir(isym2)
           igdp(isym1,isym2) = ngdp(isym)
           ngdp(isym) = ngdp(isym) + mbas1(isym1)*nbas2(isym2)
        end do
      end do

      kend1 = 1
      lwrk1 = lwork

      if (isymom.ne.muld2h(isymim,isymt2)) then
        call quit('Symmetry mismatch in ccsdr12cd!')
      end if

      ! shift it2del by one (offset vs. start address)
      do i = 1, mtotbas
        it2del(i) = it2del(i) + 1
c       write(lupri,*) 'i,it2del(i):',i,it2del(i)
      end do

*----------------------------------------------------------------------*
* transform the leading index of the E1 intermediate (for CCSD(R12)
* identical with the correlation contribution to the Fock-hat matrix) 
* from the AO to the orthogonal complementary basis:
*----------------------------------------------------------------------*
      ! allocate work space for E(p',c) intermediate
      kemat1 = kend1
      kend1  = kemat1 + nmatpb(isymim)
      lwrk1  = lwork  - kend1
      if (lwrk1 .lt. 0) call quit('Insufficient memory in CCSDR12CD 1')

      ! allocate work space for the half-transformed and AO 
      ! intermediates needed in the following section:
      kfckvao  = kend1
      kfckvaop = kfckvao  + nemat1(isymim)
      kfgdp    = kfckvaop + nviraop(isymim)
      kend2    = kfgdp    + ngdp(isym0)
      lwrk2    = lwork    - kend2
      if (lwrk2 .lt. 0) call quit('Insufficient memory in CCSDR12CD 2')

      ! read the Fhat(c,delta) and Fhat(c,delta-p) matrices
      lunit = -1
      call gpopen(lunit,'CCFHATADEL','UNKNOWN',' ','UNFORMATTED',
     &            idummy,.false.)
      read(lunit) (work(kfckvao-1+i), i=1,nemat1(isymim))
      read(lunit) (work(kfckvaop-1+i),i=1,nviraop(isymim))
      call gpclose(lunit,'KEEP')

C     if (locdbg) then
C       write(lupri,*) 'Fhat(c,delta) and Fhat(c,delta-p):'
C       do isym2 = 1, nsym
C         isym1 = muld2h(isym2,isymim)
C         write(lupri,*) 'Symmetry block ',isym1,isym2
C         call output(work(kfckvao+iemat1(isym1,isym2)),
C    &                1,nvir(isym1),1,mbas1(isym2),
C    &                nvir(isym1),mbas1(isym2),1,lupri)
C         call output(work(kfckvaop+iviraop(isym1,isym2)),
C    &                1,nvir(isym1),1,mbas2(isym2),
C    &                nvir(isym1),mbas2(isym2),1,lupri)
C       end do
C     end if

      ! read the F^val(gamma,delta-p) matrix, which contains the
      ! two-electron contribution of a Fock matrix computed from
      ! the SCF valence electron density (i.e. w/o core orbitals)
      lunit = -1
      call gpopen(lunit,'R12FOCK','UNKNOWN',' ','UNFORMATTED',idummy,
     &            .false.)
      call readt(lunit,ngdp(1),work(kfgdp))
      call gpclose(lunit,'KEEP')
 

c     -------------------------------------------------------------
c     transform leading index to virtual and substract it from
c     the corresponding matrix computed from the CC Lambda density,
c     and then transform the remaining AO index to the orthogonal
c     complementary basis to get the matrix E1(p',c):
c     -------------------------------------------------------------
      if (isymlam.ne.isymim) call quit('symmetry mismatch in ccsdr12cd')
      do isymp = 1, nsym
         isydel = isymp
         isygam = isydel
         isymc  = muld2h(isymlam,isygam)

         len1  = nvir(isymc) * mbas1(isydel)
         len2  = nvir(isymc) * mbas2(isydel)

         kscr1 = kend2
         kscr2 = kscr1 + len1
         kend3 = kscr2 + len2
         lwrk3 = lwork - kend3 
         if (lwrk3 .lt. 0) call quit('Insufficient core in CCSDR12CD 3')

         ! .........................................................
         ! initialize scratch array with the half-transformed E1
         ! or fock matrix computed from the Lambda density, with
         ! the primary AOs followed immediately by the aux. AOs
         ! .........................................................
         koffv = kfckvao + iemat1(isymc,isydel)
         call dcopy(len1,work(koffv),1,work(kscr1),1)
         koffv = kfckvaop + iviraop(isymc,isydel)
         call dcopy(len2,work(koffv),1,work(kscr2),1)

         ! .........................................................
         ! substract the contribution from the SCF density matrix
         ! with the leading index transformed to the virtual basis
         ! .........................................................
         koffl   = iglmvi(isygam,isymc) + 1
         koffgdp = kfgdp + igdp(isygam,isydel)

         nbasg = max(mbas1(isygam),1)
         nvirc = max(nvir(isymc),1)

         call dgemm('T','N',nvir(isymc),nbas2(isydel),mbas1(isygam),
     &              -one,xlamdah(koffl),nbasg,work(koffgdp),nbasg,
     &               one,work(kscr1),nvirc)

C        if (locdbg) then
C          write(lupri,*) 'in CCSDR12CD: Norm^2(FGDP) = ',
C    &      ddot(mbas1(isygam)*nbas2(isydel),work(koffgdp),1,
C    &           work(koffgdp),1)
C          call output(work(koffgdp),1,mbas1(isygam),1,nbas2(isydel),
C    &                 mbas1(isygam),nbas2(isydel),1,lupri)
C          write(lupri,*) 'in CCSDR12CD: Norm^2(SCR) = ',
C    &      ddot(nvir(isymc)*nbas2(isydel),work(kscr1),1,work(kscr1),1)
C          call output(work(kscr1),1,nvir(isymc),1,nbas2(isydel),
C    &                 nvir(isymc),nbas2(isydel),1,lupri) 
C        end if

         ! .........................................................
         ! finally transform the outer index to the orthogonal
         ! complementary basis and store at work(kemat1)
         ! .........................................................
         koffe1 = kemat1 + imatpb(isymp,isymc)
         koffcx = 1 + ilamdx(isydel,isymp) + nbas2(isydel)*norb1(isymp)

         nbasd  = max(nbas2(isydel),1)
         norbp  = max(norb2(isymp),1)

         call dgemm('T','T',norb2(isymp),nvir(isymc),nbas2(isydel),
     &              one,cmox(koffcx),nbasd,work(kscr1),nvirc,
     &              zero,work(koffe1),norbp)
        
      end do

      if (locdbg) then
        write(lupri,*) 'in CCSDR12CD: Final Norm^2(EMAT1) = ',
     &     ddot(nmatpb(isymim),work(kemat1),1,work(kemat1),1)
      end if

*----------------------------------------------------------------------*
* allocate work space for C & D intermediates:
*----------------------------------------------------------------------*
      kcdint = kend1
      kend1  = kcdint + ntg2sq(isymim)
      lwrk1  = lwork  - kend1
      if (lwrk1 .lt. 0) call quit('Insufficient memory in CCSDR12CD 4')

*----------------------------------------------------------------------*
* initialize result vector:
*----------------------------------------------------------------------*
      call dzero(omegar12,ntg2sq(isymom))

*----------------------------------------------------------------------*
* transform C intermediate to orthonormal complementary basis
*----------------------------------------------------------------------*
      iopt = 1
      lcbs = .true.
      call dzero(work(kcdint),ntg2sq(isymim))
      call cc_iajb2(work(kcdint), isymim, iopt, .false., .false., lcbs,
     &              luijda, fnijda, it2del, cmox, isym0,
     &              idummy, cdummy, idummy, dummy, idummy,
     &              work(kend1), lwrk1                               )

      if (locdbg) then
        write(lupri,*) 'Norm^2 of C-Int. after CC_IAJB2: ',
     &    ddot(ntg2sq(isymim),work(kcdint),1,work(kcdint),1)
      end if

*----------------------------------------------------------------------*
* add E1(p',c) intermediate to the k=1 diagonal of the C intermediate:
*----------------------------------------------------------------------*
      fac = -0.5d0
      call cc_cdbar2(work(kcdint),work(kemat1),fac,.true.,isymim)
      if (locdbg) then
        write(lupri,*) 'Norm^2 of C-Int. after CC_CDBAR2: ',
     &    ddot(ntg2sq(isymim),work(kcdint),1,work(kcdint),1)
      end if

*----------------------------------------------------------------------*
* calculate the C term and add to result vector:
*----------------------------------------------------------------------*
      ioptr12 = 1
      call cc_cd('C', +1, ioptr12,
     &           omegar12, isymom, t2am, isymt2,
     &           work(kcdint),  isymim, work(kend1), lwrk1       )

      if (locdbg) then
        write(lupri,*) 'Norm^2 of OmegaR12-Int. after CC_CD(C): ',
     &    ddot(ntg2sq(isymom),omegar12,1,omegar12,1)
      end if

*----------------------------------------------------------------------*
* transform D intermediate to orthonormal complementary basis
*----------------------------------------------------------------------*
      iopt = 1
      lcbs = .true.
      call dzero(work(kcdint),ntg2sq(isymim))
      call cc_iajb2(work(kcdint), isymim, iopt, .false., .false., lcbs,
     &              luiadj, fniadj, it2del, cmox, isym0,
     &              idummy, cdummy, idummy, dummy, idummy,
     &              work(kend1), lwrk1                              )

*----------------------------------------------------------------------*
* add E1(p',c) intermediate to the k=1 diagonal of the D intermediate:
*----------------------------------------------------------------------*
      fac = 0.5d0
      call cc_cdbar2(work(kcdint),work(kemat1),fac,.true.,isymim)

*----------------------------------------------------------------------*
* calculate the D term and add to result vector:
*----------------------------------------------------------------------*
      ! form in place 2C-E combination of T2
      iopttcme = 1
      call ccsd_tcmepk(t2am,one,isymt2,iopttcme)

      ioptr12 = 1
      call cc_cd('D', +1, ioptr12,
     &           omegar12, isymom, t2am, isymt2,
     &           work(kcdint),  isymim, work(kend1), lwrk1       )

      if (locdbg) then
        write(lupri,*) 'Norm^2 of OmegaR12-Int. after CC_CD(D): ',
     &    ddot(ntg2sq(isymom),omegar12,1,omegar12,1)
      end if

*----------------------------------------------------------------------*
* contract Omega(ai,p'j) with r12 integrals to Omega(xi,yj):
*----------------------------------------------------------------------*
      komegpk = 1
      kend1   = komegpk + ntr12am(isymom)
      lwrk1   = lwork - kend1
      if (lwrk1 .lt. 0) call quit('Insufficient memory in CCSDR12CD 5')

      ! read present result vector from file:
      lunit = -1 
      call gpopen(lunit,'CC_OMEGAR12','unknown',' ','unformatted',
     &                 idum,ldum)
      read(lunit) (work(komegpk+i), i=0, ntr12am(isymom)-1 )
      call gpclose(lunit,'KEEP')

      if (locdbg) then
        write(lupri,*) 'in CCSDR12CD: Norm^2 of OMEGAR12 before '//
     &                 'CD contributions: ',
     &   ddot(ntr12am(isymom),work(komegpk),1,work(komegpk),1)
      end if
     
      call ccsdr12oxr(work(komegpk),omegar12,isymom,work(kend1),lwrk1)
 
      if (locdbg) then
        write(lupri,*) 'in CCSDR12CD: Norm^2 of OMEGAR12 after '//
     &                 'CD contributions: ',
     &   ddot(ntr12am(isymom),work(komegpk),1,work(komegpk),1)
      end if

      ! write updated result vector back to file:
      lunit = -1 
      call gpopen(lunit,'CC_OMEGAR12','unknown',' ','unformatted',
     &                 idum,ldum)
      write(lunit) (work(komegpk+i), i=0, ntr12am(isymom)-1 )
      call gpclose(lunit,'KEEP')

*----------------------------------------------------------------------*
* restore it2del
*----------------------------------------------------------------------*
      ! shift it2del by one (offset vs. start address)
      do i = 1, mtotbas
        it2del(i) = it2del(i) - 1
      end do

      return
      end 
*----------------------------------------------------------------------*
*                     END OF SUBROUTINE CCSDR12CD                      *
*======================================================================*
