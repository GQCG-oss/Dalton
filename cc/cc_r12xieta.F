C
C...   Copyright (c) 2001 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 1.2
C...   (2001), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, A.A. Auer, K.L. Bak, V. Bakken,
C...   O. Christiansen, S. Coriani, P. Dahle, E. K. Dalskov,
C...   T. Enevoldsen, B. Fernandez, C. Haettig, K. Hald, A. Halkier,
C...   H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, R. Kobayashi,
C...   H. Koch, K. V. Mikkelsen, P. Norman, M. J. Packer,
C...   T. B. Pedersen, T. A. Ruden, A. Sanchez, T. Saue, S. P. A. Sauer,
C...   B. Schimmelpfennig, K. O. Sylvester-Hvid, P. R. Taylor,
C...   and O. Vahtras"
C...
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html
C
#include <single.h>
C

*=======================================================================
      subroutine cc_r12vxint(work,lwork,ltstxint)
*-----------------------------------------------------------------------
* Purpose: Calculate VX-intermediate for CC-R12 response and store 
*          on file
*
*     filer12    filename for half transformed r12 integrals with  
*                no or one auxliary basis function
*     filer12_2  filename for half transformed r12 integrals with
*                two auxiliary basis functions
*     ltstxint   only for testing:
*                if true: ignore V (i.e. set V to unit matrix),
*                and compute ONLY that contribution, which involves
*                only sums over the AO-basis (NOT the aux-basis).
*                The result should be the same as the corresponding
*                part of the R12 overlap matrix.
*
* by  C. Neiss      winter 2004/2005
*-----------------------------------------------------------------------

      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <r12int.h>
#include <ccr12int.h>
#include <ccrspprp.h>
#include <second.h>

C-----------------------------------------------------------------------
C variables/parameter:
C-----------------------------------------------------------------------
      logical locdbg,lauxdelta,lauxd,mkvamkl,lhtf,lauxbeta,ltstxint
      parameter (locdbg = .false.)
      character*8 label1,lab123(3)
      character*8 filer12,filer12_2
      integer nr1orb(8),nr1xorb(8),nr1bas(8),nr1xbas(8),nr2bas,n2bst1(8)
      integer ir1xbas(8,8),ir1xorb(8,8),ibasx(8)
      integer ir1orb(8,8),ir1bas(8,8),ir2bas(8,8),iaodis1(8,8)
      integer ir2xbas(8,8),irgkl(8,8),nrgkl(8),nvabkl(8),ivabkl(8,8) 
      integer nrxgkl(8),irxgkl(8,8)
      integer nr2bas2,ir2bas2(8,8),ir2xbas2(8,8)
      integer iaodis2(8,8),nallbas,nbasxtd(8),norbxtd(8),nlamdsxtd
      integer iaomo(8,8),naomo(8),n2bstxtd(8),iaomox(8,8)
      integer icount1,icount2,icount3,icount4
      integer nsymx,norbtsx,nbastx,nlamdsx,nrhfsx(8),norbsx(8)
      integer isymb,isymd,isymbd,isymp,isymkappa,isymbkl,isymkl,isymmn,
     &        isymmk,isymnl
      integer lwork,kend1,kend2,kend3,kend4,kend5,
     &        lwrk1,lwrk2,lwrk3,lwrk4,lwrk5
      integer kcmo,kcmox,k12back,kr12bkl,kr12bkl2,krvbkl,krvbkl2
      integer kprpao,kprpmo,kvxintsq,kvxint,koff1,koff2
      integer koffv,koffr,koffrv,koffrv2,koffc,koffres,koffx
      integer isympt,imatrix,ierr,iprp,lusifc,idel,idelta,idum
      integer luvxint
      integer kappa,beta,ienddelta
c     integer nrhftria,n2,kvxints,kvxintt

#if defined (SYS_CRAY)
      real zero,one,work(*),fac,ddot,time
#else
      double precision zero,one,work(*),fac,ddot,time
#endif
      parameter (zero = 0.0D0, one = 1.0D0)
      data lab123/'********','********','********'/

C-----------------------------------------------------------------------
C initialisation
C-----------------------------------------------------------------------
      call qenter ('cc_r12vxint')
      if (locdbg) then
        write(lupri,*) 'Entered CC_R12VXINT'
        write(lupri,*) 'memory available: ',lwork,' double words'
        call flshfo(lupri)
      end if
      time = second()

      filer12 = frhtf
      filer12_2 = frhtf2
 
C     additional dimensions/offsets:
      call cc_r12offset(nr1orb,nr1xorb,nr1bas,nr1xbas,nr2bas,
     &                  nrgkl,nrxgkl,n2bst1,ir1orb,ir1xorb,ir1bas,
     &                  ir1xbas,ir2bas,ir2xbas,irgkl,irxgkl,iaodis1,
     &                  nvabkl,ivabkl)
     
      nallbas = 0
      nlamdsxtd = 0
      do i = 1, nsym
        nbasxtd(i) = mbas1(i) + mbas2(i)
        norbxtd(i) = norb1(i) + norb2(i)
        nallbas = nallbas + nbasxtd(i)
        nlamdsxtd = nlamdsxtd + nbasxtd(i)*norbxtd(i) 
      end do

      ibasx(1) = 0
      do i = 2, nsym
        ibasx(i) = ibasx(i-1) + mbas2(i-1)
      end do

      do isymbd = 1, nsym
        icount1 = 0
        icount2 = 0
        icount3 = 0
        icount4 = 0
        do isymd = 1, nsym
          isymb = muld2h(isymbd,isymd)
          iaodis2(isymb,isymd) = icount1
          iaomo(isymb,isymd)  = icount2
          ir2bas2(isymb,isymd) = icount3
          ir2xbas2(isymb,isymd) = icount4
          icount1 = icount1 + nbasxtd(isymb)*nbasxtd(isymd)
          icount2 = icount2 + nbasxtd(isymb)*norbxtd(isymd)
          icount3 = icount3 + nr1xbas(isymb)*nr1bas(isymd)
          icount4 = icount4 + nr1xbas(isymb)*nr1xbas(isymd)
        end do
        n2bstxtd(isymbd) = icount1
        naomo(isymbd) = icount2
      end do

      nr2bas2 = 0
      do i = 1, nsym
        nr2bas2 = nr2bas2 + nr1bas(i)*nr1xbas(i)
      end do

C     if locdbg do some printout of variables:
      if (locdbg) then
        do i = 1, nsym
        write(lupri,*) 'nbas(',i,') = ',nbas(i)
        write(lupri,*) 'nbasxtd(',i,') = ',nbasxtd(i)
        write(lupri,*) 'mbas1(',i,') = ',mbas1(i)
        write(lupri,*) 'ibas(',i,') = ',ibas(i)
        write(lupri,*) 'mbas2(',i,') = ',mbas2(i)
        write(lupri,*) 'ibasx(',i,') = ',ibasx(i)
        write(lupri,*) 'nr1bas(',i,') = ',nr1bas(i)
        write(lupri,*) 'nr1xbas(',i,') = ',nr1xbas(i)
        write(lupri,*) 'nrhf(',i,') = ',nrhf(i)
        write(lupri,*) 'nvir(',i,') = ',nvir(i)
        write(lupri,*) 'norb(',i,') = ',norb(i)
        write(lupri,*) 'norbxtd(',i,') = ',norbxtd(i)  
        write(lupri,*) 'norb1(',i,') = ',norb1(i)
        write(lupri,*) 'norb2(',i,') = ',norb2(i)
        write(lupri,*) 'nmatij(',i,') = ',nmatij(i)
        write(lupri,*) 'nr1orb(',i,') = ',nr1orb(i)
        write(lupri,*) 'nr1xorb(',i,') = ',nr1xorb(i)
        write(lupri,*) 'nrgkl(',i,') = ',nrgkl(i)
        write(lupri,*) 'nrxgkl(',i,') = ',nrxgkl(i)
        write(lupri,*) 'ngamma(',i,') = ',ngamma(i)
        write(lupri,*) 'ngamsq(',i,') = ',ngamsq(i)
        end do
        do i = 1, nsym
          do j = 1, nsym
          write(lupri,*) 'irxgkl(',i,',',j,') = ',irxgkl(i,j)
          write(lupri,*) 'imatav(',i,',',j,') = ',imatav(i,j)
          write(lupri,*) 'iaodis(',i,',',j,') = ',iaodis(i,j)
          write(lupri,*) 'iaodis2(',i,',',j,') = ',iaodis2(i,j)
          end do
        end do
        write(lupri,*) 'n2basx = ',n2basx
        write(lupri,*) 'nlamdt = ',nlamdt
        write(lupri,*) 'nlamdsxtd =  ',nlamdsxtd
        call FLSHFO(LUPRI)
      end if

C-----------------------------------------------------------------------
C open output file
C-----------------------------------------------------------------------
      luvxint = 0
      call gpopen(luvxint,'CCR12VXINT','UNKNOWN',' ','UNFORMATTED',
     &            idum,.false.)

C-----------------------------------------------------------------------
C read in MO coefficients
C-----------------------------------------------------------------------
      lusifc = 0
      call gpopen(lusifc,'SIRIFC','OLD',' ','UNFORMATTED',idum,.false.)
      ! read dimensions for CMO coefficients for full basis (nlamdsx)
      rewind(lusifc)
      call mollab('FULLBAS ',lusifc,lupri)
      read(lusifc) nsymx,norbtsx,nbastx,nlamdsx,(nrhfsx(i),i=1,nsym),
     &             (norbsx(i),i=1,nsym)
C     allocate memory for MO coefficients:
      kcmo  = 1
      kend1 = kcmo + nlamdsxtd
      kcmox = kend1 
      kend2 = kcmox + nlamdsx
      lwrk2 = lwork - kend2
      if (lwrk2 .lt. 0) then
        call quit ('Insufficient work memory in CC_R12VXINT at #1')
      end if
      call dzero(work(kcmo),nlamdsxtd)
      call dzero(work(kcmox),nlamdsx)
      read(lusifc)
      read(lusifc) (work(kcmox+i-1),i=1,nlamdsx)
      call gpclose(lusifc,'KEEP')
      if (locdbg) then
        write(lupri,*) 'nsymx, norbtsx, nbastx, nlamdsx: ',
     &                  nsymx, norbtsx, nbastx, nlamdsx
        do i = 1, nsym
          write(lupri,*) 'nrhfsx(',i,') = ', nrhfsx(i)
          write(lupri,*) 'norbsx(',i,') = ', norbsx(i)
        end do
        do isymbd = 1, nsym
          icount2 = 0
          do isymd = 1, nsym
            isymb = muld2h(isymbd,isymd)
            iaomox(isymb,isymd)  = icount2
            icount2 = icount2 + nbasxtd(isymb)*norbsx(isymd)
          end do
        end do
   
        call around('MO-coefficient matrix incl. redundant orbitals')
        do i = 1, nsym
          koffc = iaomox(i,i) + kcmox
          write(lupri,*) 'Symmetry number:', i
          if (norbsx(i) .eq. 0) then
            write(lupri,*) 'This symmetry is empty'
          else
            call output(work(koffc),1,nbasxtd(i),1,norbsx(i),
     &                  nbasxtd(i),norbsx(i),1,lupri)
          end if
        end do
      end if

C     delete redundant occupied orbitals from CMO matrix:
      koff1 = kcmox 
      koff2 = kcmo
      do i = 1, nsym 
        koff1 = koff1 + nbasxtd(i)*nrhfsx(i) 
        call dcopy(nbasxtd(i)*norbxtd(i),work(koff1),1,work(koff2),1)
        koff1 = koff1 + nbasxtd(i)*norbxtd(i)
        koff2 = koff2 + nbasxtd(i)*norbxtd(i)
      end do
      if (locdbg) then
        call around('MO-coefficient matrix')
        do i = 1, nsym
          koffc = iaomo(i,i) + kcmo
          write(lupri,*) 'Symmetry number:', i
          if (norbxtd(i) .eq. 0) then
            write(lupri,*) 'This symmetry is empty'
          else
            call output(work(koffc),1,nbasxtd(i),1,norbxtd(i),
     &                  nbasxtd(i),norbxtd(i),1,lupri)
          end if
        end do
      end if
C
C-----------------------------------------------------------------------
C start loop over operators V
C-----------------------------------------------------------------------
      IF (LOCDBG) THEN
        CALL AROUND('LABELS on Operator-List:')
        DO IPRP = 1, NPRLBL
          WRITE(LUPRI,*) PRPLBL(IPRP)
        END DO
      END IF

      DO IPRP = 1, NPRLBL
        LABEL1 = PRPLBL(IPRP)
        IF (LABEL1(1:5).NE.'HAM0 ') THEN
C
C-----------------------------------------------------------------------
C read in V integrals
C-----------------------------------------------------------------------
C       allocate memory for V integrals:
        !since we do not know the symmetry of V in advance, it is more
        !memory allocated than needed in general!
        kprpao = kend1
        kend2 = kprpao + nallbas*nallbas
        lwrk2 = lwork - kend2
        if (locdbg) write(lupri,*) 'lwrk2 at #2 = ',lwrk2 
        if (lwrk2 .lt. 0) then
          call quit ('Insufficient work memory in CC_R12VXINT at #2')
        end if
        !just to be on the safe side:
        call dzero(work(kprpao),nallbas*nallbas)

        !symmetry of V is determined here
        CALL CCPRPAO(LABEL1,.FALSE.,WORK(KPRPAO),ISYMPT,IMATRIX,IERR,
     &               WORK(KEND2),LWRK2)

        if (ierr .gt. 0 ) then
          call quit('property not found on file')
        end if
        
        if (isympt .eq. 0) then
          write(lupri,*) 'WARNING: Symmetry for operator could not be'//
     &          ' determined!'
          write(lupri,*) 'WARNING: Will use IRREP = 1'
          isympt = 1
        end if

        if (locdbg) then
          call around('Original V integrals')
          write(lupri,*) 'LABEL = ',label1, 'ISYMPT = ',isympt
          do isymkappa = 1, nsym
            isymb = muld2h(isympt,isymkappa)
            koffv = iaodis2(isymb,isymkappa) + kprpao
            write(lupri,*) 'Symmetry block:', isymb, isymkappa
            if (nbasxtd(isymkappa) .eq. 0) then
              write(lupri,*) 'This symmetry is empty'
            else
              call output(work(koffv),1,nbasxtd(isymb),1,
     &                    nbasxtd(isymkappa),nbasxtd(isymb),
     &                    nbasxtd(isymkappa),1,lupri)
            end if
          end do
        end if

C-----------------------------------------------------------------------
C transform one index to MO
C-----------------------------------------------------------------------
        kprpmo = kend2
        kend3  = kprpmo + n2bstxtd(isympt)
        lwrk3  = lwork - kend3
        if (locdbg) write(lupri,*) 'lwrk3 at #3 = ',lwrk3
        if (lwrk3 .lt. 0) then
          call quit ('Insufficient work memory in CC_R12VXINT at #3')
        end if
        
        do isymkappa = 1, nsym
          isymb = muld2h(isympt,isymkappa)
          isymp = isymkappa  
          koffv = iaodis2(isymb,isymkappa) + kprpao
          koffc = iaomo(isymkappa,isymp) + kcmo
          koffres = iaomo(isymb,isymp) + kprpmo
          call dgemm('N','N',nbasxtd(isymb),norbxtd(isymp),
     &               nbasxtd(isymkappa),one,work(koffv),
     &               max(1,nbasxtd(isymb)),work(koffc),
     &               max(1,nbasxtd(isymkappa)),
     &               zero,work(koffres),max(1,nbasxtd(isymb))) 
c          if (locdbg) then
c            call around('Half transformed V integrals')
c            write(lupri,*) 'symmetry block: ',isymb, isymkappa
c            call output(work(koffres),1,nbasxtd(isymb),1,
c     &                  norbxtd(isymp),nbasxtd(isymb),norbxtd(isymp),
c     &                  1,lupri)
c            call around('MO coefficient sub-matrix')
c            call output(work(koffc),1,nbasxtd(isymkappa),1,
c     &                  norbxtd(isymp),nbasxtd(isymkappa),
c     &                  norbxtd(isymp),1,lupri)
c          end if   

C-----------------------------------------------------------------------
C backtransform to AO => one index in contravariant basis
C-----------------------------------------------------------------------
          call dgemm('N','T',nbasxtd(isymb),nbasxtd(isymkappa),
     &               norbxtd(isymp),
     &               one,work(koffres),max(1,nbasxtd(isymb)),
     &               work(koffc),max(1,nbasxtd(isymkappa)),
     &               zero,work(koffv),max(1,nbasxtd(isymb)))

        end do

        if (ltstxint) then
          !V should be the (-1)-unit matrix
          isympt = 1
          call dzero(work(kprpao),nallbas*nallbas)
          do isymkappa = 1, nsym
            isymb = isymkappa
            do kappa = 1, nbasxtd(isymkappa)
              beta = kappa
              koffv = iaodis2(isymb,isymkappa) + kprpao-1 +  
     &                nbasxtd(isymb)*(kappa-1) + beta
              work(koffv) = -one
            end do
          end do
        end if

        if (locdbg) then
          call around('Half back transformed V integrals')
          write(lupri,*) 'LABEL = ',label1
          do isymkappa = 1, nsym
            isymb = muld2h(isympt,isymkappa)
            koffv = iaodis2(isymb,isymkappa) + kprpao
            write(lupri,*) 'Symmetry block:', isymb, isymkappa
            if (nbasxtd(isymkappa) .eq. 0) then
              write(lupri,*) 'This symmetry is empty'
            else
              call output(work(koffv),1,nbasxtd(isymb),1,
     &                    nbasxtd(isymkappa),nbasxtd(isymb),
     &                    nbasxtd(isymkappa),1,lupri)
            end if
          end do
        end if

C-----------------------------------------------------------------------
C allocate memory for result arrays (square matrix and lower triangle)
C-----------------------------------------------------------------------
        kvxintsq = kend2
        kvxint   = kvxintsq + ngamsq(isympt)
        kend3 = kvxint + ngamma(isympt)
        lwrk3 = lwork - kend3
        if (locdbg) write(lupri,*) 'lwrk3 at #4 = ',lwrk3
        if (lwrk3 .lt. 0) then
          call quit('Insufficient memory in CC_R12VXINT at #4')
        end if
C       zero out result array:
        call dzero(work(kvxintsq),ngamsq(isympt))
        call dzero(work(kvxint),ngamma(isympt))

C-----------------------------------------------------------------------
C start loop over delta/delta':
C-----------------------------------------------------------------------
        do isymd = 1, nsym
 
C         allocate memory for "(r*V)" integrals:
          krvbkl  = kend3
          krvbkl2 = krvbkl + nrgkl(muld2h(isympt,isymd)) 
          kend4   = krvbkl2 + nrxgkl(muld2h(isympt,isymd))
          lwrk4   = lwork - kend4
          if (locdbg) write(lupri,*) 'lwrk4 at #5 = ',lwrk4
          if (lwrk4 .lt. 0) then
            call quit('Insufficient work memory in CC_R12VXINT at #5')
          end if
c          call dzero(work(krvbkl),nrgkl(isymd))
c          call dzero(work(krvbkl2),nrxgkl(isymd)) 

C         allocate memory for r integrals:
          kr12bkl  = kend4
          kr12bkl2 = kr12bkl + nrgkl(isymd) 
          kend5    = kr12bkl2 + nrxgkl(isymd)
          lwrk5    = lwork - kend5
          if (locdbg) write(lupri,*) 'lwrk5 at #6 = ',lwrk5
          if (lwrk5 .lt. 0) then
            call quit('Insufficient work memory in CC_R12VXINT at #6')
          end if
c          call dzero(work(kr12bkl),nrgkl(isymd))
c          call dzero(work(kr12bkl2),nrxgkl(isymd)) 

          ienddelta = nbasxtd(isymd)
          if (ltstxint) ienddelta = mbas1(isymd)
          do idelta = 1, ienddelta
            !is delta auxiliary function?
            if (idelta .gt. mbas1(isymd)) then
              lauxdelta = .true.
            else
              lauxdelta = .false.
            end if
C           determine idel:
            idel = idelta + ibas(isymd) + ibasx(isymd)

C-----------------------------------------------------------------------
C read in r integrals
C-----------------------------------------------------------------------
            lauxd = .TRUE.

            lauxbeta = .FALSE.
            ! read r^(delta)_(kappa,mn) or r^(delta')_(kappa,mn)
            call cc_r12getrint(work(kr12bkl),idel,isymd,nr1bas,ir1bas,
     &                       nr2bas,ir2bas,nrgkl,irgkl,
     &                       ir1xbas,ir2xbas,ibasx,lauxd,lauxbeta,
     &                       filer12,work(kend5),lwrk5) 
            if (.not. ltstxint) then
            ! read r^(delta)_(kappa',mn) or r^(delta')_(kappa',mn)
            lauxbeta = .TRUE.
            call cc_r12getrint(work(kr12bkl2),idel,isymd,nr1xbas,
     &                          ir1bas,nr2bas2,ir2bas2,nrxgkl,irxgkl,
     &                          ir1xbas,ir2xbas2,ibasx,lauxd,
     &                          lauxbeta,filer12_2,work(kend5),lwrk5)
            end if

            if (locdbg) then
              call around('Half transformed r integrals')
              write(lupri,*) 'idelta, idel, isymd: ',idelta,idel,isymd
              do isymkappa = 1, nsym
                isymmn = muld2h(isymd,isymkappa)
                write(lupri,*) 'Symmetry kappa, mn:',isymkappa,isymmn
                write(lupri,*) 'kappa is AO basis function:'
                koffr = irgkl(isymkappa,isymmn) + kr12bkl
                call output(work(koffr),1,mbas1(isymkappa),1,
     &                      nmatij(isymmn),mbas1(isymkappa),
     &                      nmatij(isymmn),1,lupri)
                if (min(mbas1(isymkappa),nmatij(isymmn)).eq.0) then
                  write(lupri,*) 'This symmetry block is empty'
                end if
                if (.not. ltstxint) then
                write(lupri,*) 'kappa is auxiliary basis function:'
                koffr = irxgkl(isymkappa,isymmn) + kr12bkl2
                call output(work(koffr),1,mbas2(isymkappa),1,
     &                      nmatij(isymmn),mbas2(isymkappa),
     &                      nmatij(isymmn),1,lupri)
                if (min(mbas2(isymkappa),nmatij(isymmn)).eq.0) then
                  write(lupri,*) 'This symmetry block is empty'
                end if
                end if
              end do
            end if


C-----------------------------------------------------------------------
C contract r with V
C-----------------------------------------------------------------------
            if (locdbg) then
              call around('"V*r" integrals')
              write(lupri,*) 'idelta, idel: ',idelta,idel
            end if
            do isymb = 1, nsym
              isymkappa  = muld2h(isympt,isymb)
              isymmn = muld2h(isymd,isymkappa)
              if (lauxdelta) then
                fac = one
              else 
                fac = -one
              end if
              koffv = iaodis2(isymb,isymkappa) + kprpao
              koffr = irgkl(isymkappa,isymmn) + kr12bkl 
              koffrv = irgkl(isymb,isymmn) + krvbkl
              call dgemm('N','N', mbas1(isymb),nmatij(isymmn),
     &                 mbas1(isymkappa),fac,work(koffv),
     &                 max(1,nbasxtd(isymb)),
     &                 work(koffr),max(1,mbas1(isymkappa)),zero,
     &                 work(koffrv),max(1,mbas1(isymb)))
c              if (locdbg) then
c                call output(work(koffrv),1,mbas1(isymb),1,
c     &                      nmatij(isymmn),mbas1(isymb),
c     &                      nmatij(isymmn),1,lupri)
c              end if
              koffv = iaodis2(isymb,isymkappa) +
     &              nbasxtd(isymb)*mbas1(isymkappa) 
     &              + kprpao
              koffr = irxgkl(isymkappa,isymmn) + kr12bkl2
              if (.not. ltstxint) then
              call dgemm('N','N',mbas1(isymb),nmatij(isymmn),
     &                 mbas2(isymkappa),-fac,work(koffv),
     &                 max(1,nbasxtd(isymb)),
     &                 work(koffr),max(1,mbas2(isymkappa)),one,
     &                 work(koffrv),max(1,mbas1(isymb)))
              koffv = iaodis2(isymb,isymkappa) + mbas1(isymb) + kprpao
              koffr = irgkl(isymkappa,isymmn) + kr12bkl
              koffrv2 = irxgkl(isymb,isymmn) + krvbkl2
              call dgemm('N','N',mbas2(isymb),nmatij(isymmn),
     &                 mbas1(isymkappa),-fac,work(koffv),
     &                 max(1,nbasxtd(isymb)),
     &                 work(koffr),max(1,mbas1(isymkappa)),zero,
     &                 work(koffrv2),max(1,mbas2(isymb)))
c              if (locdbg) then
c                call output(work(koffrv2),1,mbas2(isymb),1,
c     &                      nmatij(isymmn),mbas2(isymb),
c     &                      nmatij(isymmn),1,lupri)
c              end if
              koffv = iaodis2(isymb,isymkappa) +
     &              nbasxtd(isymb)*mbas1(isymkappa) + mbas1(isymb) +
     &              kprpao
              koffr = irxgkl(isymkappa,isymmn) + kr12bkl2
              call dgemm('N','N',mbas2(isymb),nmatij(isymmn),
     &                 mbas2(isymkappa),fac,work(koffv),
     &                 max(1,nbasxtd(isymb)),
     &                 work(koffr),max(1,mbas2(isymkappa)),one,
     &                 work(koffrv2),max(1,mbas2(isymb)))
              end if
              if (locdbg) then
                write(lupri,*) 'Symmetry beta, mn:',isymb,isymmn
                write(lupri,*) 'beta is AO basis function:'
                call output(work(koffrv),1,mbas1(isymb),1,
     &                      nmatij(isymmn),mbas1(isymb),
     &                      nmatij(isymmn),1,lupri)
                if (min(mbas1(isymb),nmatij(isymmn)).eq.0) then
                  write(lupri,*) 'This symmetry block is empty'
                end if
                if (.not. ltstxint) then
                write(lupri,*) 'beta is auxiliary basis function:'
                call output(work(koffrv2),1,mbas2(isymb),1,
     &                      nmatij(isymmn),mbas2(isymb),
     &                      nmatij(isymmn),1,lupri)
                if (min(mbas2(isymb),nmatij(isymmn)).eq.0) then
                  write(lupri,*) 'This symmetry block is empty'
                end if
                end if
              end if
            end do

C-----------------------------------------------------------------------
C read in R integrals (we use the memory where the r integrals were 
C stored, since these are not longer needed)
C-----------------------------------------------------------------------
            lauxd = .TRUE.

            lauxbeta = .FALSE.
            ! read R^(delta)_(beta,kl) or R^(delta')_(beta,kl)
            call cc_r12getrint(work(kr12bkl),idel,isymd,nr1bas,ir1bas,
     &                       nr2bas,ir2bas,nrgkl,irgkl,
     &                       ir1xbas,ir2xbas,ibasx,lauxd,lauxbeta,
     &                       fnback,work(kend5),lwrk5)
            if (.not. ltstxint) then
            lauxbeta = .TRUE.
            ! read R^(delta)_(beta',kl) or R^(delta')_(beta',kl)
            call cc_r12getrint(work(kr12bkl2),idel,isymd,nr1xbas,
     &                         ir1bas,nr2bas2,ir2bas2,nrxgkl,irxgkl,
     &                         ir1xbas,ir2xbas2,ibasx,lauxd,lauxbeta,
     &                         fnback2,work(kend5),lwrk5)
            end if

            if (locdbg) then
              call around('Half back transformed R integrals')
              write(lupri,*) 'idelta, idel: ',idelta,idel
              do isymb = 1, nsym
                isymkl = muld2h(isymd,isymb)
                write(lupri,*) 'Symmetry beta, kl:',isymb,isymkl
                write(lupri,*) 'beta is AO basis function:'
                koffr = irgkl(isymb,isymkl) + kr12bkl
                call output(work(koffr),1,mbas1(isymb),1,
     &                      nmatij(isymkl),mbas1(isymb),
     &                      nmatij(isymkl),1,lupri)
                if (min(mbas1(isymb),nmatij(isymkl)).eq.0) then
                  write(lupri,*) 'This symmetry block is empty'
                end if
                if (.not. ltstxint) then
                write(lupri,*) 'beta is auxiliary basis function:'
                koffr = irxgkl(isymb,isymkl) + kr12bkl2
                call output(work(koffr),1,mbas2(isymb),1,
     &                      nmatij(isymkl),mbas2(isymb),
     &                      nmatij(isymkl),1,lupri)
                if (min(mbas2(isymb),nmatij(isymkl)).eq.0) then
                  write(lupri,*) 'This symmetry block is empty'
                end if
                end if
              end do
            end if

C----------------------------------------------------------------------
C contract and add up over all deltas
C----------------------------------------------------------------------
            if (locdbg) then
              call around('VXINT before packing')
            end if
            do isymb = 1, nsym
              isymbd = muld2h(isymd,isymb)
              isymmn = muld2h(isymbd,isympt)
              isymkl = muld2h(isymbd,1)
              koffr  = irgkl(isymb,isymkl) + kr12bkl
              koffrv = irgkl(isymb,isymmn) + krvbkl
              koffx  = igamsq(isymkl,isymmn) + kvxintsq
c              call output(work(koffr),1,mbas1(isymb),1,nmatij(isymkl),
c     &                    mbas1(isymb),nmatij(isymkl),1,lupri)
c              call output(work(koffrv),1,mbas1(isymb),1,nmatij(isymmn),
c     &                    mbas1(isymb),nmatij(isymmn),1,lupri)
              call dgemm('T','N',nmatij(isymkl),nmatij(isymmn),
     &                 mbas1(isymb),one,work(koffr),max(1,mbas1(isymb)),
     &                 work(koffrv),max(1,mbas1(isymb)),one,work(koffx),
     &                 max(1,nmatij(isymkl)))
              koffr  = irxgkl(isymb,isymkl) + kr12bkl2
              koffrv2 = irxgkl(isymb,isymmn) + krvbkl2
c              call output(work(koffr),1,mbas2(isymb),1,nmatij(isymkl),
c     &                    mbas2(isymb),nmatij(isymkl),1,lupri)
c              call output(work(koffrv2),1,mbas2(isymb),1,nmatij(isymmn),
c     &                    mbas2(isymb),nmatij(isymmn),1,lupri)
              if (.not. ltstxint) then
              call dgemm('T','N',nmatij(isymkl),nmatij(isymmn),
     &                 mbas2(isymb),one,work(koffr),max(1,mbas2(isymb)),
     &                 work(koffrv2),max(1,mbas2(isymb)),one,
     &                 work(koffx),max(1,nmatij(isymkl)))
              end if
              if (locdbg) then
                write(lupri,*) 'idelta, idel, isymd: ',idelta,idel,isymd
                write(lupri,*) 'isymkl, isymmn :',isymkl,isymmn
                write(lupri,*) 'ngamsq(isympt),igamsq(isymkl,isymmn): ',
     &                          ngamsq(isympt),igamsq(isymkl,isymmn)
                call output(work(koffx),1,nmatij(isymkl),1,
     &                    nmatij(isymmn),nmatij(isymkl),nmatij(isymmn),
     &                    1,lupri)
              end if
            end do
             

C----------------------------------------------------------------------
C end of loops
C----------------------------------------------------------------------
          end do !idelta
        end do !isymd
 
C-----------------------------------------------------------------------
C resort result into a symmetry packed triangular matrix and
C apply the projection operator P_(mn)^(kl):
C-----------------------------------------------------------------------
        call ccr12pck2(work(kvxint),isympt,.TRUE.,work(kvxintsq),
     &                 .false.)

        if (locdbg) then
          call around('final VXINT before packing')
          do isymmn = 1, nsym
            isymkl = muld2h(isympt,isymmn)
            koffx = igamsq(isymkl,isymmn) + kvxintsq
            write(lupri,*) 'Symmetry block:', isymkl, isymmn
            if (nmatij(isymkl).eq.0 .or. nmatij(isymmn).eq.0) then
              write(lupri,*) 'This symmetry is empty'
            else
              call output(work(koffx),1,nmatij(isymkl),1,
     &                    nmatij(isymmn),nmatij(isymkl),nmatij(isymmn),
     &                    1,lupri)
            end if
          end do
          write(lupri,*) 
          write(lupri,*) 'Norm^2 over symmetries: ',
     &           ddot(ngamsq(isympt),work(kvxintsq),1,work(kvxintsq),1)
          call around('final VXINT after packing')
          do isymnl = 1, nsym
            isymmk = muld2h(isympt,isymnl)
            koffx = igamma(isymmk,isymnl) + kvxint 
            if (isymmk .lt. isymnl) then 
              write(lupri,*) 'Symmetry block: ', isymmk, isymnl
              if (nmatij(isymmk).eq.0 .or. nmatij(isymnl).eq.0) then
                write(lupri,*) 'This symmetry is empty' 
              else
                call output(work(koffx),1,nmatij(isymmk),1,
     &                      nmatij(isymnl),nmatij(isymmk),
     &                      nmatij(isymnl),1,lupri)
              end if
            else if (isymmk .eq. isymnl) then
              write(lupri,*) 'Symmetry block: ', isymmk, isymnl
              if (nmatij(isymmk) .eq. 0) then
                write(lupri,*) 'This symmetry is empty'
              else 
                call outpak(work(koffx),nmatij(isymmk),1,lupri)
              end if
            end if
          end do
        end if

C-----------------------------------------------------------------------
C write result on file for future use
C-----------------------------------------------------------------------
        !place date and time in lab123
        call getdat(lab123(2),lab123(3))

        write(luvxint) lab123, label1
        write(luvxint) isympt
        write(luvxint) (work(kvxint-1+i), i=1, ngamma(isympt))

        if (locdbg) then
          write(lupri,*)
          write(lupri,*) 'ngamma(isympt) = ',ngamma(isympt)
          write(lupri,*) label1, isympt
          write(lupri,*) 'norm:',ddot(ngamma(isympt),work(kvxint),1,
     &                                               work(kvxint),1)
          write(lupri,*) (work(kvxint-1+i), i=1, ngamma(isympt))
        end if

C-----------------------------------------------------------------------
C end loop over operators V
C-----------------------------------------------------------------------
        end if
      end do !operator V

C-----------------------------------------------------------------------
C close output file
C-----------------------------------------------------------------------
      call gpclose(luvxint,'KEEP')

C-----------------------------------------------------------------------
C end subroutine
C-----------------------------------------------------------------------
      time = second() - time
      WRITE(LUPRI,'(1X,A)')
     &     'Computation of CC-R12 VXINT response intermediates done'
         WRITE(LUPRI,'(/1X,A,F7.2,A)')
     &     '>>> Time used for VXINT is ',time,' seconds'

      if (locdbg) then
        write(lupri,*) 'Leaving CC_R12VXINT'
        call flshfo(lupri)
      end if
     
      call qexit('cc_r12vxint')

      return 
      end

*=======================================================================

*=====================================================================*
      subroutine cc_r12xi(XIR12,ISYMXI,TR12,ISYMC,XINT,VXINT,ISYMV,
     &                    PRPAO,CMO,LAMH,TRANS,WORK,LWORK)
C----------------------------------------------------------------------
C     purpose: calculate Xi vector for CCR12 response
C
C     XIR12    result array of dimension ngamsq(isymxi)
C     ISYMXI   symmetry of Xi
C     TR12     R12 amplitudes of dimension ngamsq(isymc)
C     ISYMC    symmetry of amplitudes
C     XINT     X-intermediate (R12 overlap integrals) of dim. ngamsq(1)
C     VXINT    VX-intermediate of dimension ngamsq(isymv)
C     ISYMV    symmetry of perturbation
C     PRPAO    one-electron perturbation operator integrals in AO-basis
C     LAMH     Lambda^h matrix
C     CMO      MO coefficient matrix
C     TRANS    use normal ('N') or transposed ('T') matrix V_m^jtilde 
C              when constructing vc-intermediate
C
C     Christian Neiss 2005 
C----------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
#include <dummy.h>

      logical locdbg 
      parameter(locdbg = .false.)

#if defined (SYS_CRAY)
      real zero,one
#else
      double precision zero,one
#endif
      parameter (one = 1.0d0)

      integer isymxi,isymv,isymc
      integer kvcint
      integer lwork,kend1,lwrk1
      character*1 trans !allowed are 'N' and 'T'

#if defined (SYS_CRAY)
      real xir12(*),tr12(*),xint(*),vxint(*),prpao(*),lamh(*),
     &     cmo(*),work(*)
#else
      double precision xir12(*),tr12(*),xint(*),vxint(*),prpao(*),
     &                 lamh(*),cmo(*),work(*)
#endif

      call qenter('cc_r12xi')

      if (locdbg) then
        write(lupri,*) 'Entered CC_R12XI'
        write(lupri,*) 'memory available: ',lwork,' double words'
      end if

      !check variable trans:
      if (.not.((trans.eq.'T').or.(trans.eq.'N'))) then
        write(lupri,*) 'TRANS = ',trans
        call quit('Forbidden value of TRANS in CC_R12XI')
      end if

      !check symmetry:
      if (isymxi.ne.muld2h(isymc,isymv)) then
        call quit('Symmetry error in CC_R12XI')
      end if 

C-----------------------------------------------------------------------
C allocate memory for intermediates needed until end
C-----------------------------------------------------------------------
      kvcint = 1                      !VC-intermediate
      kend1 = kvcint + ngamsq(isymxi)
      lwrk1 = lwork - kend1
      if (lwrk1 .lt. 0) then
        call quit('Insufficient memory in CC_R12XI')
      end if

C-----------------------------------------------------------------------
C calculate Sum_m (c_{kl}^{im}*V_m^jtilde + c_{kl}^{mj}*V_m^itilde)
C or        Sum_m (c_{kl}^{im}*V_j^mtilde + c_{kl}^{mj}*V_j^mtilde)
C and apply P_kl^ij:
C-----------------------------------------------------------------------
      call cc_r12xi2a(work(kvcint),isymxi,TR12,ISYMC,PRPAO,ISYMV,
     &                CMO,LAMH,1,TRANS,WORK(KEND1),LWRK1)
   
C-----------------------------------------------------------------------
C finally calculate Xi
C-----------------------------------------------------------------------
      !zero out Xi:
      call dzero(xir12,ngamsq(isymxi))

      call cc_r12xi2b(xir12,vxint,isymv,'N',
     &               tr12,isymc,'N',one)
      call cc_r12xi2b(xir12,xint,1,'N',
     &               work(kvcint),isymxi,'N',-one)

      if (locdbg) then
        call around('Result in CC_R12XI')
        call cc_prsqr12(xir12,isymxi,1,.false.)
      end if

      if (locdbg) write(lupri,*) 'Leaving CC_R12XI'
      call qexit('cc_r12xi')        
      return 
      end

*=======================================================================

*=====================================================================*
      subroutine cc_r12xi2a(VCINT,ISYMVC,TR12,ISYMC,PRPAO,ISYMV,
     &                      CMO,LAMH,ISYMH,TRANS,WORK,LWORK)
C----------------------------------------------------------------------
C     calculate Sum_m (c_{kl}^{im}*V_m^jtilde + c_{kl}^{mj}*V_m^itilde)
C     or        Sum_m (c_{kl}^{im}*V_j^mtilde + c_{kl}^{mj}*V_j^mtilde)
C     and apply P_kl^ij
C
C     VCINT    result array of dimension ngamsq(isymvc)
C     ISYMVC   symmetry of result
C     TR12     R12 amplitudes of dimension ngamsq(isymc)
C     ISYMC    symmetry of amplitudes
C     PRPAO    one-electron perturbation operator integrals in AO-basis
C     ISYMV    symmetry of perturbation
C     LAMH     Lambda^h matrix
C     ISYMH    Symmetry of Lambda^h matrix
C     CMO      MO coefficient matrix
C     TRANS    use normal ('N') or transposed ('T') matrix V_m^jtilde 
C              when constructing vc-intermediate
C
C     Christian Neiss 2005 
C----------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
#include <dummy.h>

      logical locdbg 
      parameter(locdbg = .false.)

#if defined (SYS_CRAY)
      real zero,one
#else
      double precision zero,one
#endif
      parameter (zero = 0.0d0, one = 1.0d0)

      integer isymvc,isymv,isymij,isymkl,isymc,isymh,isym1,isym2,
     &        isymi,isymj,isymm,isymim,isymmj
      integer kscr1,kvint,koffc,koffv,
     &        koffvc,kofftr12
      integer lwork,kend1,kend2,lwrk1,lwrk2
      character*1 trans !allowed are 'N' and 'T'

#if defined (SYS_CRAY)
      real vcint(*),tr12(*),prpao(*),lamh(*),cmo(*),work(*)
#else
      double precision vcint(*),tr12(*),prpao(*),lamh(*),cmo(*),work(*)
#endif

      call qenter('cc_r12xi2a')

      if (locdbg) then
        write(lupri,*) 'Entered CC_R12XI2A'
        write(lupri,*) 'memory available: ',lwork,' double words'
      end if

      !check variable trans:
      if (.not.((trans.eq.'T').or.(trans.eq.'N'))) then
        write(lupri,*) 'TRANS = ',trans
        call quit('Forbidden value of TRANS in CC_R12XI2A')
      end if

      !check symmetry:
      if (isymvc.ne.muld2h(isymc,muld2h(isymh,isymv))) then
        call quit('Symmetry error in CC_R12XI2A')
      end if
 
C-----------------------------------------------------------------------
C print out R12 amplitudes (debug)
C-----------------------------------------------------------------------
      if (locdbg) then
        call around('R12 amplitudes in CC_R12XI after unpacking')
        call cc_prsqr12(tr12,isymc,1,.false.)

c       !only for testing if ccr12unpck2 works right!
c       !repack amplitudes -> should be the same as read in!
c       kscr1 = kend1
c       kend2 = kscr1 + ngamma(isymc)
c       call ccr12pck2(work(kscr1),isymc,.false.,tr12,.false.)
c       call around('R12 amplitudes in CC_R12XI after repacking')
c       call cc_prpr12(work(kscr1),isymc,1,.false.)
      end if
 
C-----------------------------------------------------------------------
C transform perturbation integrals in MO-basis (occ. orbitals only)
C-----------------------------------------------------------------------
      kvint = 1
      kend1 = kvint + nmatij(muld2h(isymh,isymv))
      kscr1 = kend1
      kend2 = kscr1 + nt1ao(isymv)
      lwrk2 = lwork - kend2
      if (lwrk2 .lt. 0) then
        call quit('Insufficient memory in CC_R12XI2A')
      end if 
      do isym2 = 1, nsym
        isym1 = muld2h(isymv,isym2)
        koffv = 1 + iaodis(isym1,isym2)
        koffc = 1 + iglmrh(isym1,isym1)
        !first transformation
        call  dgemm('T','N',nrhf(isym1),nbas(isym2),nbas(isym1),
     &              one,cmo(koffc),max(1,nbas(isym1)),prpao(koffv),
     &              max(1,nbas(isym1)),zero,work(kscr1),
     &              max(1,nrhf(isym1)))
c       call output(work(kscr1),1,nrhf(isym1),1,nbas(isym2),nrhf(isym1),
c    &              nbas(isym2),1,lupri)
        !second transformation
        koffv = kvint + imatij(isym1,muld2h(isymh,isym2))
        koffc = 1 + iglmrh(isym2,muld2h(isymh,isym2))
        call dgemm('N','N',nrhf(isym1),nrhf(muld2h(isymh,isym2)),
     &             nbas(isym2),one,work(kscr1),max(1,nrhf(isym1)),
     &             lamh(koffc),max(1,nbas(isym2)),zero,work(koffv),
     &             max(1,nrhf(isym1)))
      end do

       
      if (locdbg) then
        call around('Original V integrals:')
        do isym2 = 1, nsym
          isym1 = muld2h(isymv,isym2)
          koffv = 1 + iaodis(isym1,isym2)
          write(lupri,*) 'Symmetry block: ',isym1, isym2
          if (nbas(isym1).eq.0 .or. nbas(isym2).eq.0) then
            write(lupri,*) 'This symmetry is empty'
          else
            call output(prpao(koffv),1,nbas(isym1),1,
     &                  nbas(isym2),nbas(isym1),nbas(isym2),1,lupri)
          end if
        end do
        call around('CMO matrix (occ. part):')
        do isym1 = 1, nsym 
          koffc = 1 + iglmrh(isym1,isym1)
          write(lupri,*) 'Symmetry block: ',isym1, isym1
          if (nbas(isym1).eq.0 .or. nrhf(isym1).eq.0) then
            write(lupri,*) 'This symmetry is empty'
          else
            call output(cmo(koffc),1,nbas(isym1),1,nrhf(isym1),
     &                  nbas(isym1),nrhf(isym1),1,lupri)
          end if
        end do
        call around('Lambda^h matrix (occ. part):')
        do isym2 = 1, nsym
          isym1 = muld2h(isymh,isym2)
          koffc = 1 + iglmrh(isym1,isym2)
          write(lupri,*) 'Symmetry block: ',isym1, isym2
          if (nbas(isym1).eq.0 .or. nrhf(isym2).eq.0) then
            write(lupri,*) 'This symmetry is empty'
          else
            call output(lamh(koffc),1,nbas(isym1),1,nrhf(isym2),
     &                  nbas(isym1),nrhf(isym2),1,lupri)
          end if
        end do
        call around('Transformed V integrals:')
        do isym1 = 1, nsym
          isym2 = muld2h(muld2h(isymh,isymv),isym1)
          koffv = kvint + imatij(isym1,isym2)
          write(lupri,*) 'Symmetry block: ',isym1, isym2
          if (nrhf(isym1).eq.0 .or. nrhf(isym2).eq.0) then
            write(lupri,*) 'This symmetry is empty'
          else
            call output(work(koffv),1,nrhf(isym1),1,nrhf(isym2),
     &                  nrhf(isym1),nrhf(isym2),1,lupri)
          end if
        end do
        call flshfo(lupri)
      end if

C-----------------------------------------------------------------------
C calculate Sum_m (c_{kl}^{im}*V_m^jtilde + c_{kl}^{mj}*V_m^itilde)
C or        Sum_m (c_{kl}^{im}*V_j^mtilde + c_{kl}^{mj}*V_j^mtilde)
C-----------------------------------------------------------------------
      do isymkl = 1, nsym
        isymim  = muld2h(isymkl,isymc)
        isymmj  = muld2h(isymh,isymv)
        isymij  = muld2h(isymim,isymmj)
        do isymm = 1, nsym
          isymi = muld2h(isymim,isymm)
          isymj = muld2h(isymmj,isymm)
          kofftr12 = 1 + igamsq(isymkl,isymim) + 
     &               nmatij(isymkl)*imatij(isymi,isymm)
c         koffv = kvint + imatij(isymm,isymj)
          koffvc = 1 + igamsq(isymkl,isymij) + 
     &             nmatij(isymkl)*imatij(isymi,isymj)
c         write(lupri,*) 'symmetry isymkl, isymi, isymm, isymj: ',
c    &                    isymkl,isymi,isymm,isymj
c         write(lupri,*) 'igamsq(isymkl,isymij): ',igamsq(isymkl,isymij)
c         write(lupri,*) 'nmatij(isymkl), nrhf(isymi), nrhf(isymm)',
c    &                    nmatij(isymkl), nrhf(isymi), nrhf(isymm)
c         call output(tr12(kofftr12),1,nmatij(isymkl)*nrhf(isymi),1,
c    &                nrhf(isymm),nmatij(isymkl)*nrhf(isymi),
c    &                nrhf(isymm),1,lupri)
          if (trans.eq.'N') then
            koffv = kvint + imatij(isymm,isymj)
            call dgemm('N','N',nmatij(isymkl)*nrhf(isymi),nrhf(isymj),
     &               nrhf(isymm),one,tr12(kofftr12),
     &               max(1,nmatij(isymkl)*nrhf(isymi)),
     &               work(koffv),max(1,nrhf(isymm)),zero,
     &               vcint(koffvc),max(1,nmatij(isymkl)*nrhf(isymi))) 
          else if (trans.eq.'T') then
            koffv = kvint + imatij(isymj,isymm)
            call dgemm('N','T',nmatij(isymkl)*nrhf(isymi),nrhf(isymj),
     &               nrhf(isymm),one,tr12(kofftr12),
     &               max(1,nmatij(isymkl)*nrhf(isymi)),
     &               work(koffv),max(1,nrhf(isymj)),zero,
     &               vcint(koffvc),max(1,nmatij(isymkl)*nrhf(isymi)))
          else
            call quit('Forbidden value of TRANS in CC_R12XI2A')
          end if
        end do
      end do

      if (locdbg) then
        call around('"c*V" before P_{kl}^{ij}')
        do isymkl = 1, nsym
          isymij = muld2h(isymkl,isymvc)
          koffvc = 1 + igamsq(isymkl,isymij)
          write(lupri,*) 'Symmetry block isymkl,isymij: ',isymkl,isymij
          if (nmatij(isymkl).eq.0 .or. nmatij(isymij).eq.0) then
            write(lupri,*) 'This symmetry is empty'
          else
            call output(vcint(koffvc),1,nmatij(isymkl),1,nmatij(isymij),
     &                  nmatij(isymkl),nmatij(isymij),1,lupri)
          end if
        end do
      end if

      !apply P_{kl}^{ij}
      call cc_r12pklij(vcint,isymvc,work,lwork)

      if (locdbg) then
        call around('VC-intermediate')
        do isymkl = 1, nsym
          isymij = muld2h(isymkl,isymvc)
          koffvc = 1 + igamsq(isymkl,isymij)
          write(lupri,*) 'Symmetry block isymkl,isymij: ',isymkl,isymij
          if (nmatij(isymkl).eq.0 .or. nmatij(isymij).eq.0) then
            write(lupri,*) 'This symmetry is empty'
          else
            call output(vcint(koffvc),1,nmatij(isymkl),1,nmatij(isymij),
     &                  nmatij(isymkl),nmatij(isymij),1,lupri)
          end if
        end do
      end if
    
      if (locdbg) write(lupri,*) 'Leaving CC_R12XI2A'
      call qexit('cc_r12xi2a')        
      return 
      end

*=======================================================================

*=====================================================================*
      subroutine cc_r12xi2b(RESULT,MAT1,ISYM1,TR1,MAT2,ISYM2,TR2,FAC)
C----------------------------------------------------------------------
C     purpose: calculate FAC*Sum_kl X_{mn}^{kl}*c_{kl}^{ij} and add it 
C              to a given array (RESULT)
C
C     RESULT   result array of dimension ngamsq(muld2h(isym1,isym2))
C     MAT1     Matrix 1 [X]
C     ISYM1    Symmetry of MAT1
C     MAT2     Matrix 2 [c]
C     ISYM2    Symmetry of MAT2
C     FAC      Factor
C     TR1,TR2  Transpose or not?
C
C     Christian Neiss 2005
C----------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>

      logical locdbg
      parameter(locdbg = .false.)

#if defined (SYS_CRAY)
      real zero,one
#else
      double precision zero,one
#endif
      parameter (zero = 0.0d0, one = 1.0d0)

      integer isym1,isym2,isymmn,isymij,isymkl
      integer koff1,koff2,koffr
      character*1 tr1,tr2

#if defined (SYS_CRAY)
      real result(*),mat1(*),mat2(*),fac
#else
      double precision result(*),mat1(*),mat2(*),fac 
#endif

      call qenter('cc_r12xi2b')

      if (locdbg) then
        write(lupri,*) 'Entered CC_R12XI2B'
      end if

      !check TR1, TR2
      IF (((TR1.NE.'N').AND.(TR1.NE.'T')).OR.
     &    ((TR2.NE.'N').AND.(TR2.NE.'T'))) THEN
        WRITE(LUPRI,*) 'TR1, TR2: ',TR1,TR2
        CALL QUIT('Illegal value for TR1 or TR2 in CC_R12XI2B')
      END IF

      if (locdbg) then
        call around('MAT1')
        do isymmn = 1, nsym
          isymkl = muld2h(isymmn,isym1)
          koff1 = 1 + igamsq(isymmn,isymkl)
          write(lupri,*) 'Symmetry block ',isymmn,isymkl
          if (nmatij(isymmn).eq.0 .or. nmatij(isymkl).eq.0) then
            write(lupri,*) 'This symmetry is empty'
          else
            call output(mat1(koff1),1,nmatij(isymmn),1,nmatij(isymkl),
     &                  nmatij(isymmn),nmatij(isymkl),1,lupri)
          end if
        end do
        call around('MAT2')
        do isymkl = 1, nsym
          isymij = muld2h(isymkl,isym2)
          koff2 = 1 + igamsq(isymkl,isymij)
          write(lupri,*) 'Symmetry block ',isymkl,isymij
          if (nmatij(isymkl).eq.0 .or. nmatij(isymij).eq.0) then
            write(lupri,*) 'This symmetry is empty'
          else
            call output(mat2(koff2),1,nmatij(isymkl),1,nmatij(isymij),
     &                  nmatij(isymkl),nmatij(isymij),1,lupri)
          end if
        end do
        write(lupri,*) 'TR1, TR2: ',TR1,TR2
        call flshfo(lupri)
      end if


      do isymmn = 1, nsym
        isymkl = muld2h(isymmn,isym1)
        isymij = muld2h(isymkl,isym2)
        if ((tr1.eq.'N').and.(tr2.eq.'N')) then
        koff1 = 1 + igamsq(isymmn,isymkl)
        koff2 = 1 + igamsq(isymkl,isymij)
        koffr = 1 + igamsq(isymmn,isymij)
        call dgemm(tr1,tr2,nmatij(isymmn),nmatij(isymij),
     &             nmatij(isymkl),fac,mat1(koff1),max(1,nmatij(isymmn)),
     &             mat2(koff2),max(1,nmatij(isymkl)),one,
     &             result(koffr),max(1,nmatij(isymmn)))
        else if ((tr1.eq.'T').and.(tr2.eq.'T')) then
        koff1 = 1 + igamsq(isymkl,isymmn)
        koff2 = 1 + igamsq(isymij,isymkl)
        koffr = 1 + igamsq(isymmn,isymij)
        call dgemm(tr1,tr2,nmatij(isymmn),nmatij(isymij),
     &             nmatij(isymkl),fac,mat1(koff1),max(1,nmatij(isymkl)),
     &             mat2(koff2),max(1,nmatij(isymij)),one,
     &             result(koffr),max(1,nmatij(isymmn)))
        else if ((tr1.eq.'T').and.(tr2.eq.'N')) then
        koff1 = 1 + igamsq(isymkl,isymmn)
        koff2 = 1 + igamsq(isymkl,isymij)
        koffr = 1 + igamsq(isymmn,isymij)
        call dgemm(tr1,tr2,nmatij(isymmn),nmatij(isymij),
     &             nmatij(isymkl),fac,mat1(koff1),max(1,nmatij(isymkl)),
     &             mat2(koff2),max(1,nmatij(isymkl)),one,
     &             result(koffr),max(1,nmatij(isymmn)))
        else if ((tr1.eq.'N').and.(tr2.eq.'T')) then
        koff1 = 1 + igamsq(isymmn,isymkl)
        koff2 = 1 + igamsq(isymij,isymkl)
        koffr = 1 + igamsq(isymmn,isymij)
        call dgemm(tr1,tr2,nmatij(isymmn),nmatij(isymij),
     &             nmatij(isymkl),fac,mat1(koff1),max(1,nmatij(isymmn)),
     &             mat2(koff2),max(1,nmatij(isymij)),one,
     &             result(koffr),max(1,nmatij(isymmn)))
        end if
      end do

      if (locdbg) then
        call around('Result of CC_R12XI2B')
        do isymmn = 1, nsym
          isymij = muld2h(isymmn,muld2h(isym1,isym2))
          koffr = 1 + igamsq(isymmn,isymij)
          write(lupri,*) 'Symmetry block ',isymmn,isymij
          if (nmatij(isymmn).eq.0 .or. nmatij(isymij).eq.0) then
            write(lupri,*) 'This symmetry is empty'
          else
            call output(result(koffr),1,nmatij(isymmn),1,nmatij(isymij),
     &                  nmatij(isymmn),nmatij(isymij),1,lupri)
          end if
        end do
      end if

      if (locdbg) write(lupri,*) 'Leaving CC_R12XI2B'
      call qexit('cc_r12xi2b')
      return
      end

*=======================================================================

*=====================================================================*
      subroutine cc_r12rdvxint(matrix,work,lwork,ff,isympt,labpt)
C----------------------------------------------------------------------
C     purpose: read in VXINT, calculate FF*VXINT(isympt) and 
C              add it to a given array (MATRIX)
C
C     MATRIX   array of dimension ngamsq(isympt)
C     ISYMPT   symmetry of perturbation
C     LABPT    label of perturbation
C     FF       Factor (Fieldstrength)
C     
C     Christian Neiss,  Feb. 2005, based on CC_ONEP
C----------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <dummy.h>

      logical locdbg
      parameter(locdbg = .false.)

#if defined (SYS_CRAY)
      real zero,one
#else
      double precision zero,one
#endif
      parameter (zero = 0.0d0, one = 1.0d0)

      character*8 labpt
      integer isympt,isymv
      integer kvxint,kvxintsq,luvxint
      integer lwork,kend1,lwrk1,idum

#if defined (SYS_CRAY)
      real matrix(*),work(*),ff
#else
      double precision matrix(*),work(*),ff
#endif

      call qenter('cc_r12rdvxint')

      if (locdbg) then
        write(lupri,*) 'Entered CC_R12RDVXINT'
      end if

      kvxint   = 1 
      kvxintsq = kvxint + ngamma(isympt)
      kend1    = kvxintsq + ngamsq(isympt)
      lwrk1 = lwork - kend1
      if (lwrk1 .lt. 0) then
        call quit('Insufficient memory in CC_R12RDVXINT')
      end if
 
C-----------------------------------------------------------------------
C read in VXINT
C-----------------------------------------------------------------------
      luvxint = 0
      call gpopen(luvxint,'CCR12VXINT','OLD',' ','UNFORMATTED',
     &        idum,.false.)
      rewind(luvxint)
      call mollab(labpt,luvxint,lupri)
      read(luvxint) isymv
      if (isymv .ne. isympt) then
        write(lupri,*) 'LABEL, ISYMV, ISYMPT: ',LABPT, ISYMV, ISYMPT
        call quit('Symmetry mismatch when reading CCR12VXINT in '//
     &            'CC_R12RDVXINT')
      end if
      read(luvxint) (work(kvxint+i-1),i=1,ngamma(isympt))
      call gpclose(luvxint,'KEEP')
      !unpack to square:
      call ccr12unpck2(work(kvxint),isympt,work(kvxintsq),.false.)

      if (locdbg) then
        write(lupri,*) 'LABEL = ',labpt
        call around('VXINT in CC_R12RDVXINT')
        call cc_prsqr12(work(kvxintsq),isympt,1,.false.)
      end if

C-----------------------------------------------------------------------
C multiply with FF and add to MATRIX 
C-----------------------------------------------------------------------
      call daxpy(ngamsq(isympt),ff,work(kvxintsq),1,matrix,1) 

      if (locdbg) write(lupri,*) 'Leaving CC_R12RDVXINT'
      call qexit('cc_r12rdvxint')
      return
      end

*=======================================================================

*=====================================================================*
      subroutine cc_r12pklij(MATRIX,ISYM,WORK,LWORK)
C----------------------------------------------------------------------
C     purpose: calculate P_{kl}^{ij}X_{kl}^{ij}
C                       = X_{kl}^{ij} + X_{lk}^{ji} 
C
C     MATRIX   Matrix of dimension ngamsq(ISYM) with {kl} and {ij} as
C              packed pair indices; MATRIX is input and output!
C
C     Christian Neiss,  Feb. 2005
C----------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>

      logical locdbg
      parameter(locdbg = .false.)

#if defined (SYS_CRAY)
      real zero,one
#else
      double precision zero,one
#endif
      parameter (zero = 0.0d0, one = 1.0d0)

      integer isym,isymkl,isymij,isymk,isyml,isymi,isymj
      integer idxkl,idxlk,idxij,idxji,idxklij,idxlkji
      integer kscr1,kscr2,koff
      integer lwork,kend1,lwrk1

#if defined (SYS_CRAY)
      real work(*),matrix(*),dnrm2
#else
      double precision work(*),matrix(*),dnrm2
#endif


      call qenter('cc_r12pklij')

      if (locdbg) then
        write(lupri,*) 'Entered CC_R12PKLIJ'
        write(lupri,*) 'memory available: ',lwork,' double words'
      end if

      !only for testing needed, see below!
      if (locdbg) then
        kscr1 = 1
        kscr2 = kscr1 + ngamma(isym)
        kend1 = kscr2 + ngamsq(isym)
        lwrk1 = lwork - kend1
        if (lwrk1 .lt. 0) then
          call quit('Insufficient memory in CC_R12PKLIJ')
        end if

        call dcopy(ngamsq(isym),matrix,1,work(kscr2),1)
      end if
      
C calculate X_{lk}^{ji}      
      do isymkl = 1, nsym
        isymij = muld2h(isym,isymkl)
        do isymk = 1, nsym
          isyml = muld2h(isymkl,isymk) 
          do isymi = 1, nsym
            isymj = muld2h(isymij,isymi)
            do k = 1, nrhf(isymk)
              do l = 1, nrhf(isyml)
                idxkl = imatij(isymk,isyml) + nrhf(isymk)*(l-1) + k
                idxlk = imatij(isyml,isymk) + nrhf(isyml)*(k-1) + l
                do i = 1, nrhf(isymi)
                  do j = 1, nrhf(isymj)
                    idxij = imatij(isymi,isymj) + nrhf(isymi)*(j-1)+i
                    idxji = imatij(isymj,isymi) + nrhf(isymj)*(i-1)+j
                    idxklij = igamsq(isymkl,isymij) + 
     &                        nmatij(isymkl)*(idxij-1) + idxkl
                    idxlkji = igamsq(isymkl,isymij) + 
     &                        nmatij(isymkl)*(idxji-1) + idxlk
                    if (idxklij.le.idxlkji) then
                      matrix(idxklij)=matrix(idxklij)+matrix(idxlkji)
                      matrix(idxlkji)=matrix(idxklij)
                    end if
                  end do
                end do
              end do
            end do 
          end do
        end do
      end do
 
      if (locdbg) then
        call around('Result in CC_R12PKLIJ')
        call cc_prsqr12(matrix,isym,1,.false.)
      end if
     
      !Test: calling these two routines should give the same as before!
      if (locdbg) then
        call ccr12pck2(work(kscr1),isym,.TRUE.,work(kscr2),.false.)
        call ccr12unpck2(work(kscr1),isym,work(kscr2),.false.)
        call around('alternatively calculated result in CC_R12PKLIJ')
        call cc_prsqr12(work(kscr2),isym,1,.false.)

        call daxpy(ngamsq(isym),-one,matrix,1,work(kscr2),1)
        write(lupri,*) 'Norm^2 of difference: ',
     &    dnrm2(ngamsq(isym),work(kscr2),1)
      end if     

      if (locdbg) write(lupri,*) 'Leaving CC_R12PKLIJ'
      call qexit('cc_r12pklij')
      return
      end

*=======================================================================

*=====================================================================*
      subroutine cc_r12tstform(WORK,LWORK)
C----------------------------------------------------------------------
C     purpose: test the routines for reordering of amplitudes etc.
C              pure test routine!
C
C     Christian Neiss,  Feb. 2005
C----------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <dummy.h>
#include <r12int.h>
#include <ccr12int.h>

#if defined (SYS_CRAY)
      real zero,one,work(*)
#else
      double precision zero,one,work(*)
#endif
      parameter (zero = 0.0d0, one = 1.0d0)

      integer lwork, lwrk0, kend0
      integer KTR12SQ,KTR12PK1,KTR12PK2,KTAMP12S,KTAMP12T
      integer NRHFTRIA     
      integer lunit,idum,ij,n2
      logical ldum


      call qenter('cc_r12tstform')


      nrhftria = nrhft*(nrhft+1)/2 
      n2 = nrhftria*nrhftria

      KTR12SQ = 1
      KEND0   = KTR12SQ + NGAMSQ(1)

      KTR12PK1 = KEND0
      KTR12PK2 = KTR12PK1 + NGAMMA(1)
      KEND0    = KTR12PK2 + NGAMMA(1)

      KTAMP12S = KEND0
      KTAMP12T = KTAMP12S + NRHFTRIA*NRHFTRIA
      KEND0    = KTAMP12T + NRHFTRIA*NRHFTRIA

      LWRK0 = LWORK - KEND0
      if (lwrk0.lt.0) then
        call quit('Not enough memory in CC_R12TSTFORM')
      end if

  
C     method 1:
      CALL CC_R12GETCT(WORK(KTR12SQ),1,0,ketscl,.FALSE.,.TRUE.,
     &                 DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,
     &                 WORK(KEND0),LWRK0)
      call cc_prsqr12(WORK(KTR12SQ),1,1,.true.)
      CALL CCR12PCK2(WORK(KTR12PK1),1,.FALSE.,WORK(KTR12SQ),.true.)
      call cc_prpr12(WORK(KTR12PK1),1,1,.true.)


C     method 2:
C     CALL GPOPEN(lunit,fccr12c,'old',' ','formatted',idum,ldum)
C     REWIND(LUNIT)
C     READ(LUNIT,'(4E30.20)') (WORK(KTAMP12S+IJ), IJ = 0, N2-1)
C     READ(LUNIT,'(4E30.20)') (WORK(KTAMP12T+IJ), IJ = 0, N2-1)
C     CALL GPCLOSE(LUNIT,'KEEP')
C     call around('R12 singlet part')
C     call output(work(KTAMP12S),1,nrhftria,1,nrhftria,
C    &            nrhftria,nrhftria,1,lupri)
C     call around('R12 triplet part')
C     call output(work(KTAMP12T),1,nrhftria,1,nrhftria,
C    &            nrhftria,nrhftria,1,lupri)
C     CALL CCR12PCK(WORK(KTR12PK2),1,WORK(KTAMP12S),WORK(KTAMP12T))
C     call cc_prpr12(WORK(KTR12PK2),1,1,.true.)

C     reconstruct singlet/triplet format from method 1:
      call ccr12unpck(WORK(KTR12PK1),1,WORK(KTAMP12S),WORK(KTAMP12T))
      call around('R12 singlet part')
      call output(work(KTAMP12S),1,nrhftria,1,nrhftria,
     &            nrhftria,nrhftria,1,lupri)
      call around('R12 triplet part')
      call output(work(KTAMP12T),1,nrhftria,1,nrhftria,
     &            nrhftria,nrhftria,1,lupri)

C     resonstruct square matrix format from method 2:
      call ccr12unpck2(WORK(KTR12PK2),1,WORK(KTR12SQ),.true.)
      call cc_prsqr12(WORK(KTR12SQ),1,1,.true.)

C     test P_kl^ij:
      ! Method 1:
      call cc_r12pklij(WORK(KTR12SQ),1,work(kend0),lwrk0)
      call cc_prsqr12(WORK(KTR12SQ),1,1,.true.)
      ! Method 2:
      call cc_r12vunpack(WORK(KTR12SQ),1,work(KTAMP12S),work(KTAMP12T),
     &                   nrhftria,.TRUE.)
      call cc_prsqr12(WORK(KTR12SQ),1,1,.true.)


      call quit('TEST DONE')
      call qexit('cc_r12tstform')
      return
      end

*=======================================================================

*=======================================================================
      subroutine cc_r12eta0(etar12sq,cmo,isyres,work,lwork)
C-----------------------------------------------------------------------
C  purpose: r12 contribution to eta^(0) 
C
C           etar12sq   r12 part of eta^(0) vector of dim. ngamsq(isyres)
C           cmo        MO coefficient or lambda matrix
C           isyres     symmetry of result (for eta = 1)
C
C  Note that the result is ADDED to the input!
C
C  Christian Neiss  Mar. 2005
C-----------------------------------------------------------------------

      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
#include <ccr12int.h>
#include <dummy.h>

      logical locdbg
      parameter (locdbg = .false.)

      integer isyres,isymak,isymk,isyma
      integer kvajkl,kvijkl,kend1,lwrk1
      integer lwork,luvajkl
      
#if defined (SYS_CRAY)
      real cmo(*),etar12sq(*),WORK(LWORK),one,two,ddot
#else
      double precision cmo(*),etar12sq(*),WORK(LWORK),one,two,ddot
#endif

      parameter (one = 1.0D0, two = 2.0D0)

      call qenter('cc_r12eta0')
      if (locdbg) then
        write(lupri,*) 'Entered CC_R12ETA0'
        write(lupri,*) 'memory available: ',lwork,' double words'
      end if

C     do isymak = 1, nsym
C       nvajkl(isymak) = 0
C       do isymk = 1, nsym
C         isyma = muld2h(isymak,isymk)
C         nvajkl(isymak) = nvajkl(isymak) + nt1ao(isyma)*nmatij(isymk)
C       end do
C     end do
 
      kvajkl = 1
      kvijkl = kvajkl + nvajkl(1)
      kend1 = kvijkl + ngamsq(isyres)
      lwrk1 = lwork - kend1
      if (lwrk1 .lt. 0) then
        call quit('Insufficient work space in cc_r12eta0')
      end if

C-----------------------------------------------------------------------
C     Read in V(alpha j,kl)
C-----------------------------------------------------------------------
      luvajkl = 0
      call gpopen(luvajkl,fvajkl,'old',' ','unformatted',
     &     idummy,.false.)
      rewind(luvajkl)
      read(luvajkl) (work(kvajkl+i-1), i = 1,nvajkl(1))
      call gpclose(luvajkl,'KEEP')
      if (locdbg) then
         write(lupri,*) 'fvajkl: ', fvajkl
         write(lupri,*) 'norm^2(V(alpha j,kl)):',
     &    ddot(nvajkl(1),work(kvajkl),1,work(kvajkl),1)
      end if      

C-----------------------------------------------------------------------
C     Calculate CMO*V
C-----------------------------------------------------------------------
      call dzero(work(kvijkl),ngamsq(isyres))
      call cc_r12mkvijkl(work(kvajkl),1,cmo,isyres,work(kend1),
     &                   lwrk1,.true.,one,work(kvijkl))
      if (locdbg) then
         write(lupri,*) 'norm^2(vijkl) after cc_r12mkvijkl:',
     &    ddot(ngamsq(isyres),work(kvijkl),1,work(kvijkl),1)
      end if

C-----------------------------------------------------------------------
C     Make 2*Coulomb-Exchange: 2V(ij,kl)-V(ji,kl)
C-----------------------------------------------------------------------
      call cc_r12tcmesq(work(kvijkl),isyres,.false.)

      call daxpy(ngamsq(isyres),two,work(kvijkl),1,etar12sq,1)

      if (locdbg) write(lupri,*) 'Leaving CC_R12ETA0'
      call qexit('cc_r12eta0')
      return
      end
*=======================================================================

*=======================================================================
      subroutine cc_r12prop(propr12,labelh,aproxr12,work,lwork)
C-----------------------------------------------------------------------
C  purpose: r12 contribution to expectation value of property <label> 
C
C           propr12     r12 part of property
C           labelh      label of property
C           aproxr12    r12 approximation
C
C  Christian Neiss  Mar. 2005
C-----------------------------------------------------------------------

      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <r12int.h>
#include <ccr12int.h>
#include <dummy.h>

      logical locdbg
      parameter (locdbg = .false.)

      character*8 labelh
      CHARACTER*3 APROXR12, LISTL, filxi
      CHARACTER*10 MODEL
      integer idxxi,idlstl,lenmod,iopt,isymxi,isymv,isym,idx 
      integer kxir12,kzetar12,kzeta1,kzeta2,kzetar12sq,kxir12sq 
      integer ian,luxint,ierr
      integer ktr12,ktr12sq,kxint,kxintsq,kvxint,kvxintsq,
     &        kprpao,klamp0,klamh0,kt1amp0
      integer lwork,lwrk1,kend1,lwrk2,kend2

      ! external function:
      integer IRHSR1,ILSTSYM       

#if defined (SYS_CRAY)
      real work(*),propr12,ddot,one,diacont
#else
      double precision work(*),propr12,ddot,one,diacont
#endif

      parameter(one = 1.0D0)

      call qenter('cc_r12prop')
      if (locdbg) then
        write(lupri,*) 'Entered CC_R12PROP'
        write(lupri,*) 'memory available: ',lwork,' double words'
      end if

      ! only total symmetric Xi
      isymxi = 1

      kxir12 = 1
      kzetar12 = kxir12 + ngamma(isymxi)
      kzetar12sq = kzetar12 + ngamma(1)
      kend1  = kzetar12sq + ngamsq(1)
      lwrk1 = lwork - kend1
      if (locdbg) then
        kzeta1 = kend1
        kzeta2 = kzeta1 + nt1amx
        kend1 = kzeta2 + nt2amx
      end if
      if (lwrk1 .lt. 0) then
          call quit ('Insufficient work memory in CC_R12PROP')
      end if

C-----------------------------------------------------------------------
C Read Xi vector from file if it exists...
C-----------------------------------------------------------------------
      ! check if this is correct !
c      idxxi = indprpcc(labelh)
c      idxxi  = IRHSR1(labelh,.FALSE.,0.0D0,isym)
c      if (idxxi .ne. -1) then
c        if (isym .ne. isymxi) then
c          call quit('Symmetry error in CC_R12PROP')
c        end if
c        iopt = 32
c        filxi  = 'R1 '
c        call cc_rdrsp(filxi,idxxi,isymxi,iopt,model,dummy,work(kxir12))
c      else
C-----------------------------------------------------------------------
C ...otherwise calculate Xi vector
C-----------------------------------------------------------------------
C       write(lupri,*) 'Will calculate R12 part for ',labelh
        call flshfo(lupri)
        ! allocate memory
        ktr12   = kend1
        ktr12sq = ktr12 + ngamma(1)
        kxir12sq = ktr12sq + ngamsq(1)
        kxint   = kxir12sq + ngamsq(1)
        kxintsq = kxint + ngamma(1)
        kvxintsq  = kxintsq + ngamsq(1)
        kprpao  = kvxintsq + ngamsq(isymxi)
        kt1amp0 = kprpao + n2bst(isymxi)
        klamp0  = kt1amp0 + NT1AMX
        klamh0  = klamp0 + NLAMDT
        kend2   = klamh0 + NLAMDT
        LWRK2   = LWORK - KEND2
        IF (LWRK2 .LT. 0) THEN
          CALL QUIT('Insufficient work space in CC_R12PROP')
        END IF

        ! read R12 amplitudes from disk
        iopt = 32
        call cc_rdrsp('R0 ',0,1,iopt,model,dummy,work(ktr12))
        call ccr12unpck2(work(ktr12),1,work(ktr12sq),.false.)

        ! read X-integrals (R12 overlap matrix)
        luxint = 0
        call gpopen(luxint,fccr12x,'old',' ','unformatted',idummy,
     &              .false.)
        rewind(luxint)
 9999   read(luxint) ian
        read(luxint) (work(kxint+i), i=0, ngamma(1)-1 )
        if (ian.ne.ianr12) goto 9999
        call gpclose(luxint,'KEEP')
        call ccr12unpck2(work(kxint),1,work(kxintsq),.false.)

        ! read in VXINT
        call dzero(work(kvxintsq),ngamsq(isymxi))
        call cc_r12rdvxint(work(kvxintsq),work(kend2),lwrk2,one,
     &                     isymxi,labelh)

        ! read in V (perturbation operator) matrix in AO-basis
        call ccprpao(labelh,.TRUE.,work(kprpao),isymv,isym,
     &               ierr,work(kend2),lwrk2)
        IF ((IERR.GT.0) .OR. (IERR.EQ.0 .AND. isymv.NE.isymxi)) THEN
          CALL QUIT('CC_XIETA1: error while reading operator '//LABELH)
        ELSE IF (IERR.LT.0) THEN
          CALL DZERO(work(kprpao),N2BST(isymxi))
        END IF

        ! generate Lambda matrices
        iopt = 1
        call cc_rdrsp('R0 ',0,1,iopt,model,work(kt1amp0),dummy)
        call lammat(work(klamp0),work(klamh0),work(kt1amp0),
     &              work(kend2),lwrk2)

        ! calulate R12 part of Xi
        call cc_r12xi(work(kxir12sq),isymxi,work(ktr12sq),1,
     &                work(kxintsq),work(kvxintsq),isymxi,work(kprpao),
     &                work(klamp0),work(klamh0),'N',work(kend2),lwrk2)

        ! pack Xi to triangular format
        call ccr12pck2(work(kxir12),isymxi,.false.,work(kxir12sq),
     &                 .false.)
        call cclr_diasclr12(work(kxir12),brascl,isymxi)

c      end if

      if (locdbg) then
        call around('Xi R12 part in CC_R12PROP')
        call cc_prpr12(work(kxir12),isymxi,1,.FALSE.)
        call FLSHFO(LUPRI)
      end if

C-----------------------------------------------------------------------
C Read Lagrangian multipliers from file
C-----------------------------------------------------------------------
      listl = 'L0 '
      idlstl = 0
      iopt = 32
      call cc_rdrsp(listl,idlstl,1,iopt,model,dummy,work(kzetar12))

C     if (locdbg) then
C       call around('Lagrangian multipliers R12 part as read from file')
C       call cc_prpr12(work(kzetar12),1,1,.FALSE.)
C     end if

      if (locdbg) then
        !Read conventional multipliers (for comparison)
        IOPT = 3
        CALL CC_RDRSP(listl,idlstl,1,IOPT,MODEL,WORK(kzeta1),
     &                WORK(kzeta2))
 
        call around('Lagrangian multipliers')
        call cc_prp(work(kzeta1),work(kzeta2),1,1,1)
        call cc_prpr12(work(kzetar12),1,1,.TRUE.)
      end if

C-----------------------------------------------------------------------
C Calculate dot product
C-----------------------------------------------------------------------
      propr12 = ddot(ngamma(1),work(kzetar12),1,work(kxir12),1)

      if (locdbg) then
        diacont = 0.0D0
        do isym = 1, nsym
          do i = 1, nmatij(isym)
            idx = igamma(isym,isym) + i*(i+1)/2 - 1
            diacont = diacont + work(kzetar12+idx)*work(kxir12+idx)
          end do
        end do 
        write(lupri,*) 'R12 contribution analysis for operator ',
     &                 labelh 
        write(lupri,*) 'propr12 = ', propr12
        write(lupri,*) 'Diagonal contribution = ',diacont
        write(lupri,*) 'Off-diagonal cont. = ',propr12-diacont
      end if
      
      if (locdbg) write(lupri,*) 'Leaving CC_R12PROP'
      call qexit('cc_r12prop')
      return
      end
*=======================================================================

*=======================================================================
      subroutine cc_r12etaa(ETAA,ISYRES,CTR12,ISYCTR,TR12,ISYT12,XINT,
     &                      PRPAO,ISYMV,LAMDP,LAMDH,LAO,WORK,LWORK)
C-----------------------------------------------------------------------
C  purpose: r12 contribution to Eta{A} singles part 
C
C  ETAA     singles part of Eta{A}
C  ISYRES   symmetry of Eta{A}
C  CTR12    R12 doubles Lagr. multipliers (ngamsq(isyctr))
C  ISYCTR   symmetry of Lagr. multipliers 
C  TR12     R12 doubles amplitudes / trial vector (ngamsq(isyt12))
C  ISYT12   R12 amplitudes symmetry
C  XINT     R12 overlap matrix (ngamsq(1))
C  PRPAO    Perturbation operator integrals in AO basis
C  ISYMV    Symmtetry of perturbation operator
C  LAMDP    CMO matrix used for occ. index (assumed symm. 1) 
C  LAMDH    CMO matrix used for vir. index (assumed symm. 1)
C  LAO      flag: compute Eta{A}_(alpha i), i.e. "vitual" index in AO
C                 basis (skip second transformation of PRPAO)
C
C  Christian Neiss  April 2005
C-----------------------------------------------------------------------

      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <r12int.h>
#include <ccr12int.h>
#include <dummy.h>

      logical locdbg, lao
      parameter (locdbg = .false.)

      integer isyres,isyctr,isym1,isym2,isymkl,isymmj,isymmi,isymm,
     &        isymi,isymj,isyma,isyt12,isymv
      integer lwork,kend1,kend2,lwrk1,lwrk2,kvint,kscr1,kscr2
      integer koffv,koffc,koffscr1,koffscr2,koffres,koffctr

#if defined (SYS_CRAY)
      real work(*),etaa(*),ctr12(*),tr12(*),xint(*),prpao(*),
     &     lamdp(*),lamdh(*),zero,one
#else
      double precision work(*),etaa(*),ctr12(*),tr12(*),xint(*),
     &                 prpao(*),lamdp(*),lamdh(*),zero,one
#endif

      parameter(zero = 0.0D0, one = 1.0D0)

      call qenter('cc_r12etaa')
      if (locdbg) then
        write(lupri,*) 'Entered CC_R12ETAA'
        write(lupri,*) 'memory available: ',lwork,' double words'
      end if

      !check symmetry:
      if (isyres .ne. muld2h(isyctr,muld2h(isyt12,isymv))) then
        call quit('Symmetry error in CC_R12ETAA')
      end if

      if (locdbg) then
        call around('Eta{O} singles before R12 contribution')
        call cc_prp(etaa,dummy,isyres,1,0)
      end if

C--------------------
C allocate memory
C--------------------
      kvint = 1
      if (lao) then
        kend1 = kvint + nt1ao(isymv)
      else
        kend1 = kvint + nt1am(isymv)
      end if
      lwrk1 = lwork - kend1
      if (lwrk1 .lt. 0) then
        call quit('Insufficient work space in CC_R12ETAA')
      end if

C-----------------------------------------------------
C transform perturbation integrals in MO-basis: V(a,j)
C-----------------------------------------------------
      kscr1 = kend1
      kend2 = kscr1 + nt1ao(isymv)
      lwrk2 = lwork - kend2
      if (lwrk2 .lt. 0) then
        call quit('Insufficient work space in CC_R12ETAA')
      end if

      do isym2 = 1, nsym
        isym1 = muld2h(isymv,isym2)
        koffv = 1 + iaodis(isym1,isym2)
        koffscr1 = kscr1 + it1ao(isym1,isym2)
        koffc = 1 + iglmrh(isym2,isym2)
        !first transformation
        call  dgemm('N','N',nbas(isym1),nrhf(isym2),nbas(isym2),
     &              one,prpao(koffv),max(1,nbas(isym1)),lamdp(koffc),
     &              max(1,nbas(isym2)),zero,work(koffscr1),
     &              max(1,nbas(isym1)))
        !second transformation
        if (lao) then
          call dcopy(nt1ao(isymv),work(koffscr1),1,work(kvint+
     &               it1ao(isym1,isym2)),1)
        else
          koffv = kvint + it1am(isym1,isym2)
          koffc = 1 + iglmvi(isym1,isym1)
          !Note that the virtual index is leading in V(a,j)!
          call dgemm('T','N',nvir(isym1),nrhf(isym2),nbas(isym1),
     &               one,lamdh(koffc),max(1,nbas(isym1)),work(koffscr1),
     &               max(1,nbas(isym1)),zero,work(koffv),
     &               max(1,nvir(isym1)))
        end if
      end do

      if (locdbg ) then
        call around('Lambda^p matrix, all active orbitals:')
        do isym2 = 1, nsym
          isym1 = muld2h(isymv,isym2)
          write(lupri,*) 'Symmetry block: ',isym1, isym2
          call output(lamdp(1+iglmrh(isym1,isym2)),1,nbas(isym1),1,
     &                norb(isym2),nbas(isym1),norb(isym2),1,lupri)
        end do
        call around('Lambda^p matrix, occ. part:')
        do isym2 = 1, nsym
          isym1 = muld2h(isymv,isym2)
          write(lupri,*) 'Symmetry block: ',isym1, isym2
          call output(lamdp(1+iglmrh(isym1,isym2)),1,nbas(isym1),1,
     &                nrhf(isym2),nbas(isym1),nrhf(isym2),1,lupri)
        end do
        call around('Lambda^p matrix, virt. part:')
        do isym2 = 1, nsym
          isym1 = muld2h(isymv,isym2)
          write(lupri,*) 'Symmetry block: ',isym1, isym2
          call output(lamdp(1+iglmvi(isym1,isym2)),1,nbas(isym1),1,
     &                nvir(isym2),nbas(isym1),nvir(isym2),1,lupri)
        end do
        call around('Original V integrals:')
        do isym2 = 1, nsym
          isym1 = muld2h(isymv,isym2)
          write(lupri,*) 'Symmetry block: ',isym1, isym2
          call output(prpao(1+iaodis(isym1,isym2)),1,nbas(isym1),1,
     &                nbas(isym2),nbas(isym1),nbas(isym2),1,lupri)
        end do
        call around('Half transformed V integrals:')
        do isym2 = 1, nsym
          isym1 = muld2h(isymv,isym2)
          write(lupri,*) 'Symmetry block: ',isym1, isym2
          call output(work(kscr1+it1ao(isym1,isym2)),1,nbas(isym1),
     &                1,nrhf(isym2),nbas(isym1),nrhf(isym2),1,lupri)
        end do
      end if

      if (locdbg.and.(.not.lao)) then
        call around('Perturbation operator integrals V(a,j)')
        call cc_prp(work(kvint),dummy,isymv,1,0)
      end if

      if (locdbg) then
        call around('R12 overlap matrix X')
        call cc_prsqr12(xint,1,1,.false.)
        call around('R12 ground state amplitudes')
        call cc_prsqr12(tr12,isyt12,1,.false.)
        call around('R12 Lagrangian multipliers')
        call cc_prsqr12(ctr12,isyctr,1,.false.)
      end if

C----------------------
C make contractions
C----------------------
      kscr1 = kend1
      kscr2 = kscr1 + ngamsq(isyt12)
      kend2 = kscr2 + nmatij(muld2h(isyctr,isyt12))
      lwrk2 = lwork - kend2
      if (lwrk2 .lt. 0) then
        call quit('Insufficient work space in CC_R12ETAA')
      end if

      ! first contraction over k'l':
      call dzero(work(kscr1),ngamsq(isyt12))
      call cc_r12xi2b(work(kscr1),xint,1,'N',tr12,isyt12,'N',one)

      if (locdbg) then
        call around('c*X')
        call cc_prsqr12(work(kscr1),isyt12,1,.false.)
      end if

      ! second contraction over klm:
      call dzero(work(kscr2),nmatij(muld2h(isyctr,isyt12)))
      do isymkl = 1, nsym
        isymmj = muld2h(isymkl,isyt12)
        isymmi = muld2h(isymkl,isyctr)
        do isymm = 1, nsym
          isymj = muld2h(isymmj,isymm)
          isymi = muld2h(isymmi,isymm)
          koffscr1 = kscr1 + igamsq(isymkl,isymmj) + 
     &               nmatij(isymkl)*imatij(isymm,isymj) 
          koffctr  = 1 + igamsq(isymkl,isymmi) + 
     &               nmatij(isymkl)*imatij(isymm,isymi)
          koffscr2 = kscr2 + imatij(isymj,isymi)
C         write(lupri,*) 'symmetry isymkl, isymm: ',
C    &                    isymkl,isymm
C         write(lupri,*) 'SCR1:'
C         call output(work(koffscr1),1,nmatij(isymkl)*nrhf(isymm),1,
C    &                nrhf(isymj),nmatij(isymkl)*nrhf(isymm),
C    &                nrhf(isymj),1,lupri)
C         write(lupri,*) 'CTR:'
C         call output(ctr12(koffctr),1,nmatij(isymkl)*nrhf(isymm),1,
C    &                nrhf(isymi),nmatij(isymkl)*nrhf(isymm),
C    &                nrhf(isymi),1,lupri)
          call dgemm('T','N',nrhf(isymj),nrhf(isymi),
     &               nmatij(isymkl)*nrhf(isymm),one,work(koffscr1),
     &               max(1,nmatij(isymkl)*nrhf(isymm)),ctr12(koffctr),
     &               max(1,nmatij(isymkl)*nrhf(isymm)),one,
     &               work(koffscr2),max(1,nrhf(isymj)))
C         write(lupri,*)'SCR2: imatij(isymj,isymi)=',imatij(isymj,isymi)
C         call output(work(koffscr2),1,nrhf(isymj),1,
C    &                nrhf(isymi),nrhf(isymj),nrhf(isymi),1,lupri)
        end do
      end do

      if (locdbg) then
        call around('ctr*c*X')
        do isymi = 1, nsym
          isymj = muld2h(isymi,muld2h(isyctr,isyt12))
          koffscr2 = kscr2 + imatij(isymj,isymi)
          write(lupri,*) 'Symmetry block: ',isymj, isymi
          if (nrhf(isymj).eq.0 .or. nrhf(isymi).eq.0) then
            write(lupri,*) 'This symmetry is empty'
          else
            call output(work(koffscr2),1,nrhf(isymj),1,nrhf(isymi),
     &                  nrhf(isymj),nrhf(isymi),1,lupri)
          end if
        end do
      end if

      ! last contraction over j:
      if (lao) then
        do isymj = 1, nsym
          isyma = muld2h(isymv,isymj)
          isymi = muld2h(isymj,muld2h(isyctr,isyt12))
          koffscr2 = kscr2 + imatij(isymj,isymi)
          koffv    = kvint + it1ao(isyma,isymj)
          koffres  = 1 + it1ao(isyma,isymi)
c         call output(work(koffv),1,nbas(isyma),1,nrhf(isymj),
c    &                nbas(isyma),nrhf(isymj),1,lupri)
c         call output(work(koffscr2),1,nrhf(isymj),1,nrhf(isymi),
c    &                nrhf(isymj),nrhf(isymi),1,lupri)
          call dgemm('N','N',nbas(isyma),nrhf(isymi),nrhf(isymj),-one,
     &               work(koffv),max(1,nbas(isyma)),work(koffscr2),
     &               max(1,nrhf(isymj)),one,etaa(koffres),
     &               max(1,nbas(isyma)))
        end do
      else
        do isymj = 1, nsym
          isyma = muld2h(isymv,isymj)
          isymi = muld2h(isymj,muld2h(isyctr,isyt12))
          koffscr2 = kscr2 + imatij(isymj,isymi)
          koffv    = kvint + it1am(isyma,isymj)
          koffres  = 1 + it1am(isyma,isymi)
          call dgemm('N','N',nvir(isyma),nrhf(isymi),nrhf(isymj),-one,
     &               work(koffv),max(1,nvir(isyma)),work(koffscr2),
     &               max(1,nrhf(isymj)),one,etaa(koffres),
     &               max(1,nvir(isyma))) 
        end do
      end if

      if (locdbg) then
        if (lao) then
          call around('Eta{O}_(alpha j)')
          do isym2 = 1, nsym
            isym1 = muld2h(isyres,isym2)
            write(lupri,*) 'Symmetry block: ',isym1, isym2
            call output(etaa(1+it1ao(isym1,isym2)),1,nbas(isym1),1,
     &                nrhf(isym2),nbas(isym1),nrhf(isym2),1,lupri)
          end do    
        else
          call around('Eta{O} singles after R12 contribution')
          call cc_prp(etaa,dummy,isyres,1,0)
        end if
      end if

      if (locdbg) write(lupri,*) 'Leaving CC_R12ETAA'
      call qexit('cc_r12etaa')
      return
      end
*=======================================================================


