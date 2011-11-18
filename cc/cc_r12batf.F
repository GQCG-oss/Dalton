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
      subroutine r12batf(r2am,fback,fback2,lauxbeta,projvir,work,lwork)
c----------------------------------------------------------------------
c   Purpose: back transforms two indices of the integrals
c            (kp|r_12|lq) integrals to the contravariant AO basis
c            and store them on file
c
c   (k gamma|r_12|l delta) = sum_pq C_gamma,p C_delta,q (kp|r_12|lq)
c
c   where:  k,l      active occupied / r12 MOs
c           p,q      general MOs (orbital + auxiliary basis)
c           gamma    AO in orbital basis
c           delta    AO in orbital or auxiliary basis
c
c           r2am     contains (kp|r_12|lq) stored as symmetry
c                    packed triangular matrix using it2am,nt1am,it1am
c           fback    file name for back transformed integrals
c           fback2   file name for back transformed integrals with
c                    two indices from the auxiliary basis
c           lauxbeta flag whether to calculate also integrals with
c                    p transformed to the auxiliary basis
c           projvir  include only virtual orbitals in transformation
c
c   Note: this routine is used in the R12 environment with a
c         number of dimension defined differently form the
c         definitions used in the Coupled Cluster program:
c  
c         mbas1   AOs in orbital basis
c         mbas2   AOs in auxiliary basis
c         norb1   orbital basis functions per irrep
c         norb2   auxiliary basis functions per irrep
c         nvir    general MOs per irrep ( = norb1 + norb2 )
c         nrhf    active occupied / r12 MOs
c         nrhfa   active occupied MOs
c         nvircc  active virtual MOs
c         nbas    mbas1 + mbas2
c
c  Heike Fliegl, Christof Haettig  spring 2003
c  adapted for max. two auxiliary functions: Christian Neiss, fall 2004
c  introduced projvir option: Christof Haettig, spring 2005
c  adapted for CABS, Ansatz 3: Christian Neiss, winter 2006
c----------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>

      integer isymr2
      parameter (isymr2 = 1)

      logical locdbg
      parameter (locdbg = .false.)

      character*(*) fback,fback2
      logical lauxbeta, projvir, modifyinttyp1

#if defined (SYS_CRAY)
      real zero,one
#else
      double precision zero,one
#endif
      parameter (zero = 0.0D0, one = 1.0D0)

      integer nr1orb(8),nrgkl(8),irgkl(8,8)
      integer ir1orb(8,8),ir1bas(8,8),nr1bas(8),ir1xbas(8,8),nr2xbas(8) 
      integer nr1xbas(8),ir2bas(8,8),ir2xbas(8,8),n2bst1(8),iaodis1(8,8)
      integer nr1xorb(8),ir1xorb(8,8),nrxgkl(8),irxgkl(8,8)
      integer nr2bas2,ir2bas2(8,8),ir2xbas2(8,8),nvircc(8)
      integer isymql,lur2back,isymp,isympk,isymk,isyml
      integer lwork,isymq,isymd,idelta,idel,isym,isym1,isym2
      integer icount3,icount4
      integer kcmo,lwrk0,lusifc,idummy
      integer nr2bas,kcmobas,kcmoaux,lwrk1,koff1,ntotg1
      integer kr2sm,kend2,kr2gk,kr2pk,lwrk2,istartq,iendq,ntotdk,isymg
      integer koffl,kofft,ntotp,ntotg,ndell,iadr,len
      integer kend0,kend1,kend3,lwrk3
      integer iendtyp,inttyp,istartdelta,ienddelta,istartgamma
      integer istartp,iendp,norbp,norbq,nbasg,lur2back2
      integer nalphaj(8),ialphaj(8,8)
celena
      integer isymdk,kr2sm2,kr2gk2,kr2pk2,icount1
      integer koffs 
celena
#if defined (SYS_CRAY)
      real r2am(*), work(*), dnrm2, ddot
#else
      double precision r2am(*), work(*), dnrm2, ddot
#endif

c-----------------------------------------------------------------------
c     some initializations: 
c         set label for trace routines
c         calculate a number of symmetry offsets and dimensions
c-----------------------------------------------------------------------
      call qenter('r12batf')
      if(locdbg) write(lupri,*)'entered r12batf'

      call cc_r12offset(nr1orb,nr1xorb,nr1bas,nr1xbas,nr2bas,
     & nrgkl,nrxgkl,n2bst1,ir1orb,ir1xorb,ir1bas,ir1xbas,ir2bas,ir2xbas,
     & irgkl,irxgkl,iaodis1,nalphaj,ialphaj)

      ! nr2bas2 = lenghts for nr1bas x nr1xbas over all symmetries 
      ! nvircc  = number of active virtual orbitals in a symmetry
      nr2bas2 = 0
      do isym = 1, nsym
         nr2bas2 = nr2bas2 + nr1bas(isym)*nr1xbas(isym)
         nvircc(isym) = norb1(isym) - nrhfsa(isym) - nvirfr(isym)
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

celena
c      do isymdk = 1, nsym
c        nr2orb(isymdk) = 0
c        icount1 = 0
c        do isymk = 1, nsym
c           isymd = muld2h(isymdk,isymk)
c           ir2orb(isymd,isymk) = icount1
c           icount1 = icount1 + nrhf(isymk)*norb2(isymd)
c           nr2orb(isymdk) = nr2orb(isymdk) + norb2(isymd)*nrhf(isymk)
c        end do
c      end do
celena
c  nr2orb:=nr1xorb
c  ir2orb:=ir1xorb


      if (locdbg) then
        do i = 1, nsym
          write(lupri,*) 'nrhf(',i,') = ',nrhf(i)
          write(lupri,*) 'norb1(',i,') = ',norb1(i)
          write(lupri,*) 'norb2(',i,') = ',norb2(i)
          write(lupri,*) 'nvir(',i,') = ',nvir(i)
          write(lupri,*) 'mbas1(',i,') = ',mbas1(i)
          write(lupri,*) 'mbas2(',i,') = ',mbas2(i)
          write(lupri,*) 'nbas(',i,') = ',nbas(i)
        end do
      end if

c-----------------------------------------------------------------------
c     get MO coefficients (defined for combined orbital + aux. basis):
c-----------------------------------------------------------------------
      kcmo  = 1
      kend0 = kcmo + nlamds
      lwrk0 = lwork - kend0
      if (lwrk0 .lt.0) then
         call quit('Insufficient work space in R12BATF')
      end if
      
      lusifc = -1
      call gpopen(lusifc,'SIRIFC','OLD',' ','UNFORMATTED',
     &            idummy,.false.)
      rewind(lusifc)
      call mollab('FULLBAS ',lusifc,lupri)
      read(lusifc)
      read(lusifc)
      read(lusifc) (work(kcmo+i-1),i=1,nlamds)
      call gpclose(lusifc,'KEEP')
      !Reorder C_MO coefficients such that first all occ. orbitals
      !in all symmetries, then all virtual(here in MP2-R12: = ALL orbitals!!)
      !orbitals in all symmetries; by this, the redundant occ. 
      !orbitals are grouped in one single block!
      call cmo_reorder(work(kcmo),work(kend0),lwrk0)

C-----------------------------------------------------------------------
C     Loop over types of integrals inttyp
C-----------------------------------------------------------------------
      if (projvir) then
        iendtyp = 1   ! integrals of type (kp|r12|lq)
      else if (.not. lauxbeta) then
        iendtyp = 2   ! integrals of type (kp|r12|lq) and (kp|r12|lq')
        if (mbsmax.lt.5) then
          write (lupri,*) 'MBSMAX = ',MBSMAX 
          call quit("MBSMAX must at least be 5 for (kp|r12|lq')")
        end if
      else if (lauxbeta) then
        iendtyp = 4   ! additionally (kp'|r12|lq) and (kp'|r12|lq')
        if (mbsmax.lt.6) then
          write (lupri,*) 'MBSMAX = ',MBSMAX 
          call quit("MBSMAX must at least be 6 for (kp'|r12|lq')")
        end if
      end if

      modifyinttyp1 = ianr12.eq.1 .and. r12cbs

      if (locdbg) then
        write (lupri,*) 'MBSMAX = ', MBSMAX
      end if

      do inttyp = 1, iendtyp
       
c-----------------------------------------------------------------------
c     open file for back transformed r12 integrals
c-----------------------------------------------------------------------
        if ((inttyp .eq. 1) .or. (inttyp .eq. 2)) then
          lur2back = -1
          call wopen2(lur2back,fback,64,0)
        else if ((inttyp .eq. 3) .or. (inttyp .eq. 4)) then
          lur2back2 = -1
          call wopen2(lur2back2,fback2,64,0)
        end if
 
C-----------------------------------------------------------------------
C       Loop over symmetries of delta
C-----------------------------------------------------------------------
        do isymd = 1, nsym
          isymq = isymd 

C-----------------------------------------------------------------------
C         Set start and end values for orbital spaces (part1)
C-----------------------------------------------------------------------
          if (modifyinttyp1 .and. (inttyp .eq. 1)) then
            istartdelta = 1
            ienddelta   = mbas1(isymd)
            istartq   = 1
            iendq     = norb1(isymq) + norb2(isymq)
            norbq     = norb1(isymq) + norb2(isymq)
          elseif ((inttyp .eq. 1) .or. (inttyp .eq. 3)) then
            istartdelta = 1
            ienddelta   = mbas1(isymd)
           if (projvir) then
            istartq   = nrhfsa(isymq) + 1
            iendq     = nrhfsa(isymq) + nvircc(isymq)
            norbq     = nvircc(isymq)
           else
            istartq   = 1
            iendq     = norb1(isymq)
            norbq     = norb1(isymq)
           end if
          else if ((inttyp .eq. 2) .or. (inttyp .eq. 4)) then
            istartdelta = mbas1(isymd) + 1
            ienddelta  = mbas1(isymd) + mbas2(isymd)
            istartq   = norb1(isymq) + 1
            iendq     = norb1(isymq) + norb2(isymq)
            norbq     = norb2(isymq)
          end if

C-----------------------------------------------------------------------
C         Start loop over delta
C-----------------------------------------------------------------------
          do idelta = istartdelta, ienddelta
            idel = ibas(isymd) + idelta

            if (locdbg) write(lupri,*) 'isymd, idelta, idel: ',
     &                                  isymd, idelta, idel

c           ------------------------------------------------------------
c           get row vector out of transformation matrix C_delta,q
c           for fixed delta:
c           ------------------------------------------------------------
            kcmobas = kend0
            if (modifyinttyp1 .and. inttyp .eq. 1) then
              kend1 = kcmobas + norb1(isymd) + norb2(isymd)  
              koff1 = kcmo - 1 + iglmvi(isymd,isymq) + idelta
            elseif ((inttyp.eq.1) .or. (inttyp.eq.3)) then
              if (projvir) then
                kend1   = kcmobas + norbq
                koff1 = kcmo - 1 + iglmvi(isymd,isymq) +
     &                  nbas(isymd)*(istartq-1) + idelta
              else
                kend1 = kcmobas + norb1(isymd)
                koff1 = kcmo - 1 + iglmvi(isymd,isymq) + idelta
              end if
            else
              kend1 = kcmobas + norb2(isymd)
              koff1 = kcmo - 1 + iglmvi(isymd,isymq) +
     &                nbas(isymd)*norb1(isymq) + idelta
            end if
            lwrk1   = lwork - kend1
            if (lwrk1 .lt. 0) then
              call quit('Insufficient work space in R12BATF')
            end if
            call dcopy(norbq,work(koff1),nbas(isymd),
     &                    work(kcmobas),1)
            if (locdbg) then
              write(lupri,*) 'C_MO vector for fixed delta:'
              write(lupri,*) 'idelta, C_MO: ',idelta
              call output(work(kcmobas),1,1,1,kend1-kcmobas,
     &                    1,kend1-kcmobas,1,lupri)
            end if

c           ------------------------------------------------------------
c           loop over symmetry blocks of the pair index (pk) ->
c           determines also symmetry of pair (ql) and of index l
c           ------------------------------------------------------------
            do isympk = 1, nsym
              isymql = muld2h(isympk,isymr2)
              isyml  = muld2h(isymql,isymq)

              if (nrhf(isyml) .ne. 0) then

C               --------------------------------------------------------
C               allocate memory for:
C               kr2gk  : R(gamma,k)    matrix AO * occupied
C               --------------------------------------------------------
                kr2gk = kend1
                if ((inttyp .eq. 1) .or. (inttyp .eq. 2)) then
                  if (modifyinttyp1 .and. (inttyp .eq. 1)) then
                     kr2gk2= kr2gk + nr1bas(isympk)
                     kend2 = kr2gk2+ nr1bas(isympk)
                  else 
                     kend2 = kr2gk + nr1bas(isympk)
                  endif
                else if ((inttyp .eq. 3) .or. (inttyp .eq. 4)) then
                  kend2 = kr2gk + nr1xbas(isympk)
                end if
                lwrk2 = lwork - kend2
                if (lwrk2 .lt. 0) then
                  call quit('Insufficient work in R12BATF')
                end if  

C               --------------------------------------------------------
C               Loop over occupied orbitals l
C               --------------------------------------------------------
                do l = 1, nrhf(isyml)

C                 ------------------------------------------------------
C                 Loop over symmetries of gamma
C                 ------------------------------------------------------
                  do isymg = 1, nsym
                    isymp = isymg
                    isymk = muld2h(isympk,isymp)
                 
C                   ----------------------------------------------------
C                   Set start and end values for orbital spaces (part 2)
C                   ----------------------------------------------------
                    if ((inttyp .eq. 1) .or. (inttyp .eq. 2)) then
                     istartgamma = 1
                     nbasg       = mbas1(isymg)
                     if (projvir) then
                      istartp = nrhfsa(isymp) + 1
                      iendp   = nrhfsa(isymp) + nvircc(isymp)
                      norbp   = nvircc(isymp)
                     else
                      istartp = 1
                      iendp   = norb1(isymp)
                      norbp   = norb1(isymp)
                     end if
                    else if ((inttyp .eq. 3) .or. (inttyp .eq. 4)) then
                      istartgamma = mbas1(isymg) + 1
                      nbasg       = mbas2(isymg)
                      istartp = norb1(isymp) + 1
                      iendp   = norb1(isymp) + norb2(isymp)
                      norbp   = norb2(isymp)
                    end if
                 
c                   ----------------------------------------------------
c                   allocate work space for:
c                   kr2sm  : r(pk,q)  three-index array batch of r2am
c                   kr2pk  : R(pk)    matrix MO * occupied
c                   ----------------------------------------------------
                    kr2sm = kend2
                    if (modifyinttyp1 .and. inttyp .eq. 1) then
                       kr2sm2= kr2sm + norb1(isymp)*nrhf(isymk)*
     &                          (norb1(isymq)+norb2(isymq))
                       kr2pk = kr2sm2+ norb2(isymp)*nrhf(isymk)*
     &                           norb1(isymq)
                       kr2pk2= kr2pk + norb1(isymp)*nrhf(isymk)
                       kend3 = kr2pk2+ norb2(isymp)*nrhf(isymk)
                    else 
                       kr2pk = kr2sm + norbp * nrhf(isymk) * norbq
                       kend3 = kr2pk + norbp * nrhf(isymk) 
                    endif
                    lwrk3 = lwork - kend3
                    if (lwrk3 .lt. 0) then
                     call quit('Insufficient work space in R12BATF')
                    end if

                    call cc_r12getrpkq(fback,work(kr2sm),l,isyml,
     &                                 istartq,iendq,isymq,istartp,
     &                                 iendp,isymp,r2am,
     &                                 isympk,isymql,nr1orb,ir1orb)
                    if (modifyinttyp1 .and. (inttyp .eq.1)) then
                       call cc_r12getrpkq(fback,work(kr2sm2),l,isyml,
     &                                    1,norb1(isymq),isymq,
     &                                    norb1(isymp)+1,
     &                                    norb1(isymp) +norb2(isymp),
     &                                    isymp,r2am,
     &                                    isympk,isymql,nr1orb,ir1orb)
                    endif
                    if (locdbg) then
                      write(lupri,*) 'r^(l)_(pkq):'
                      call output(work(kr2sm),1,norbp*nrhf(isymk),1,
     &                            norbq,norbp*nrhf(isymk),norbq,
     &                            1,lupri)
                      if (modifyinttyp1 .and. (inttyp .eq.1)) then
                        write(lupri,*) 'r^(l)_(pkq) 2:'
                        call output(work(kr2sm2),1,norb2(isymp)*
     &                              nrhf(isymk),1,norb1(isymq),
     &                              norb2(isymp)*nrhf(isymk),
     &                              norb1(isymq),1,lupri)
                      end if 
                    end if

c                   ----------------------------------------------------
c                   transform q to delta (contravariant AO basis)
c                   ----------------------------------------------------
                    ntotdk = max(norbp*nrhf(isymk),1)
                    !zero out result (dgemv may give wrong result)!
                    call dzero(work(kr2pk),norbp*nrhf(isymk))
                    call dgemv('N',norbp*nrhf(isymk),
     &                          norbq,one,
     &                          work(kr2sm),ntotdk,
     &                          work(kcmobas),1,
     &                          zero,work(kr2pk),1)
                    if (modifyinttyp1 .and. (inttyp .eq.1)) then
                       ntotdk = max(norb2(isymp)*nrhf(isymk),1)
                       !zero out result (dgemv may give wrong result)!
                       call dzero(work(kr2pk2),norb2(isymp)*nrhf(isymk))
                       call dgemv('N',norb2(isymp)*nrhf(isymk),
     &                            norb1(isymq),one,
     &                             work(kr2sm2),ntotdk,
     &                             work(kcmobas),1,
     &                             zero,work(kr2pk2),1)
                    endif

c                   ----------------------------------------------------
c                   transform p to gamma (contravariant AO basis)
c                   ----------------------------------------------------
                    koffl = kcmo + iglmvi(isymg,isymp)
     &                       + nbas(isymg)*(istartp-1) + istartgamma-1
                    if ((inttyp .eq. 1) .or. (inttyp .eq. 2)) then
                      kofft = kr2gk + ir1bas(isymg,isymk)
                    else if ((inttyp .eq. 3) .or. (inttyp .eq. 4)) then
                      kofft = kr2gk + ir1xbas(isymg,isymk)
                    end if

                    ntotp  = max(norbp,1)
                    ntotg  = max(nbas(isymg),1)
                    ntotg1 = max(nbasg,1)
                    call dgemm('N','N',nbasg,nrhf(isymk),norbp,
     &                         one,work(koffl),ntotg,
     &                             work(kr2pk),ntotp,
     &                         zero,work(kofft),ntotg1) 
                    if (modifyinttyp1 .and. (inttyp .eq. 1)) then
                       koffl=kcmo + norb1(isymp)*nbas(isymg)
     &                            + iglmvi(isymg,isymp)
                       koffs=kr2gk2 + ir1bas(isymg,isymk)
                       ntotp=max(norb2(isymp),1)
                       call dgemm('N','N',mbas1(isymg),nrhf(isymk),
     &                       norb2(isymp),one,work(koffl),ntotg,
     &                       work(kr2pk2),ntotp,zero,work(koffs),ntotg1)
                       call daxpy(mbas1(isymg)*nrhf(isymk),one,
     &                            work(koffs),1,work(kofft),1)
                    endif

                  end do ! loop over symmetries of gamma

c                 ----------------------------------------------------
c                 store R(gamma k, l, delta) as R_(gamma k),(delta l)
c                 on file:
c                 (note: isympk = isymgk, isymql = isymdl)
c                 ----------------------------------------------------
                  if (inttyp .eq. 1) then
                    ndell = ir1bas(isymd,isyml) + 
     &                      mbas1(isymd)*(l-1) + idelta
                    iadr  = ir2bas(isympk,isymql) + 
     &                      nr1bas(isympk)*(ndell -1) + 1
                    len   = nr1bas(isympk)
                    call putwa2(lur2back,fback,work(kr2gk),
     &                          iadr,len)
                  else if (inttyp .eq. 2) then
                    ndell = ir1xbas(isymd,isyml) +
     &                      mbas2(isymd)*(l-1) + idelta - mbas1(isymd)
                    iadr  = nr2bas + ir2xbas(isympk,isymql) +
     &                      nr1bas(isympk)*(ndell -1) + 1
                    len   = nr1bas(isympk)
                    call putwa2(lur2back,fback,work(kr2gk),
     &                          iadr,len)
                  else if (inttyp .eq. 3) then
                    ndell = ir1bas(isymd,isyml) +
     &                      mbas1(isymd)*(l-1) + idelta
                    iadr  = ir2bas2(isympk,isymql) +
     &                      nr1xbas(isympk)*(ndell -1) + 1
                    len   = nr1xbas(isympk)
                    call putwa2(lur2back2,fback2,work(kr2gk),
     &                          iadr,len)
                  else if (inttyp .eq. 4) then 
                    ndell = ir1xbas(isymd,isyml) +
     &                      mbas2(isymd)*(l-1) + idelta - mbas1(isymd)
                    iadr  = nr2bas2 + ir2xbas2(isympk,isymql) +
     &                      nr1xbas(isympk)*(ndell -1) + 1  
                    len   = nr1xbas(isympk)
                    call putwa2(lur2back2,fback2,work(kr2gk),
     &                          iadr,len)
                  end if
                  if (locdbg) then
                    write(lupri,'(a,2i5,g20.10)') 
     &                'iadr,len,norm^2(R)=',iadr,len,
     &                 ddot(len,work(kr2gk),1,work(kr2gk),1)
                    write(lupri,*) 'R^(l delta)_(gamma k):'
                    call output(work(kr2gk),1,len,1,1,
     &                          len,1,1,lupri)
                    write(lupri,*)'nr2bas,ir2xbas',nr2bas,
     &                   ir2xbas(isympk,isymql)
                    write(lupri,*)'ndell',ndell
                    write(lupri,*)'ir1xbas',ir1xbas(isymd,isyml)
                    write(lupri,*)'mbas2',mbas2(isymd)
                    write(lupri,*)
                  end if

                end do ! loop over l
              end if
            end do ! loop over symmetry blocks (pk)
          end do ! loop over delta
        end do ! loop over symmetries of delta

c-----------------------------------------------------------------------
c      close file, clean trace and return:
c-----------------------------------------------------------------------
        if ((inttyp .eq. 1) .or. (inttyp .eq. 2)) then
          call wclose2(lur2back,fback,'KEEP')
        else 
          call wclose2(lur2back2,fback2,'KEEP')
        end if

      end do ! loop over integral types

      if(locdbg) write(lupri,*)'leaving cc_r12batf'
      call qexit('r12batf')

      return
      end 
*=====================================================================*
      subroutine cc_r12getrpkq(fback,rpkq,l,isyml,istartq,iendq,isymq,
     &                         istartp,iendp,isymp,
     &                         r2am,isympk,isymql,nr1orb,ir1orb)
c---------------------------------------------------------------------- 
c     Purpose:  collect a three index batch r^l_{pk,q} with fixed l
c               and q and symmetry of p, ranging from istartq to iendq 
c               and from istartp to iendp out of the
c               symmetry packed lower triangular matrix r_{pk,ql} 
c               stored in r2am
c
c     rpkq = r^l_{pk,q} with :   p MOs ranging from istartp ... iendp 
c                                k occupied / r12 MOs
c                                q MOs ranging from istartq ... iendq
c
c     r2am = r_{pk,ql}  with :  p,q general MOs (orbital + aux. basis)
c                               k,l occupied / r12 MOs
c
c     H. Fliegl, C. Haettig  spring 2003
c     modified C. Neiss      autumn 2004
c----------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
!for first order properties, app. B (Elena):
#include <ccr12int.h>

      logical locdbg
      character*(*) fback !for first order properties, app. B (Elena)
      parameter (locdbg=.false.)
      integer ir1orb(8,8),nr1orb(8),npk,npk2,npkq,isymk,isymp
      integer npkql,istartq,iendq,nql,istartp,iendp,norbp
      integer isymq,isyml,isympk,isymql, iqloff
      integer index
#if defined (SYS_CRAY)
      real rpkq(*), r2am(*)
#else
      double precision rpkq(*), r2am(*)
#endif
      index(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j
c
      if (locdbg) then
        write(lupri,*) 'ONEAUX = ',ONEAUX
      end if

c
c     number of p orbitals
      norbp = iendp - istartp + 1
c
c     symmetry of k
      isymk = muld2h(isymp,isympk)      

      iqloff = nh1am(isymql) * (nh1am(isymql)+1) / 2

      if (locdbg) then
         write(lupri,*) 'istartq,iendq,istartp,iendp,l: ',
     &                   istartq,iendq,istartp,iendp,l
         write(lupri,*)
         write(lupri,*) 'Input: r_{pk,ql}'
         call outpak(r2am,nbast*nrhf(isymk),1,lupri)
      end if

c     ------------------------------------------------
c     loop over q and index k
c     ------------------------------------------------
      do q = istartq, iendq

        if (.not.oneaux) then
          nql = it1am(isymq,isyml) + nvir(isymq)*(l-1) + q
        else if (q. le. norb1(isymq)) then
          nql = ih1am(isymq,isyml) + norb1(isymq)*(l-1) + q
        else
          nql = ig1am(isymq,isyml) + norb2(isymq)*(l-1) + q-norb1(isymq)
        end if

        do k = 1, nrhf(isymk)

          if (locdbg) then
            write(lupri,*)'nql,isymk,k',nql,isymk,k
          end if

c           ---------------------------------------------------------
c           reorder integrals, the following indeces are used:
c            npk   : pair p,k with k occupied, p gen. MO (orb.+aux.)
c            npk2  : pair p,k with k occupied, p MO (orbital only)
c            npkql : quadrupel pk,ql as needed for r2am
c            npkq  : tripel    pk,q  as needed for rpkq
c           ---------------------------------------------------------
            if (isympk.eq.isymql) then 
              do p = istartp, iendp
                if (oneaux) then
                  npk = ih1am(isymp,isymk)+norb1(isymp)*(k-1)+p
                  if (q. le. norb1(isymq)) then
                    npkql = ih2am(isymql,isympk)+index(nql,npk)
                  else
                    npkql = ih2am(isymql,isympk)+iqloff
     *                    + nh1am(isympk)*(nql-1)+npk
                  end if
                else
                  npk   = it1am(isymp,isymk) + nvir(isymp)*(k-1)+p
                  if (fback.eq. fq12back) then
                     npkql = it2sq(isymql,isympk)+ nt1am(isymql)
     &                       *(npk-1)+nql
                  elseif (fback.eq. fu12back) then
                     npkql = it2sq(isympk,isymql)+ nt1am(isympk)
     &                      *(nql-1)+npk
                  else
                     npkql = it2am(isympk,isymql)+index(npk,nql)
                  endif
                end if
                npk2  = norbp*(k-1)+(p-istartp+1)
                npkq  = norbp*nrhf(isymk)*(q-istartq)+npk2
                
                if (locdbg) then
                   write(lupri,'(a,5i5)')'p,npk,npk2,npkql,npkq',
     &                  p,npk,npk2,npkql,npkq
                   call flshfo(lupri)
                end if

                rpkq(npkq) = r2am(npkql)
              end do
            else if (isympk.lt.isymql) then 
              if (oneaux) 
     &          call quit('isympk.ne.isymql has not been implemented')
              do p = istartp, iendp
                npk   = it1am(isymp,isymk) + nvir(isymp)*(k-1)+p 
                npk2  = norbp*(k-1)+(p-istartp+1)
                npkql = it2am(isympk,isymql)+ nt1am(isympk)*(nql-1)+npk
                npkq  = norbp*nrhf(isymk)*(q-istartq)+npk2
                rpkq(npkq) = r2am(npkql)
              end do
            else if (isympk.gt.isymql) then 
              if (oneaux) 
     &          call quit('isympk.ne.isymql has not been implemented')
              do p = istartp, iendp
                npk   = it1am(isymp,isymk) + nvir(isymp)*(k-1)+p 
                npk2  = norbp*(k-1)+(p-istartp+1)
                npkql = it2am(isymql,isympk)+nt1am(isymql)*(npk-1)+nql
                npkq  = norbp*nrhf(isymk)*(q-istartq)+npk2
                rpkq(npkq) = r2am(npkql)
              end do
            end if
        end do
      end do

      end 
*=====================================================================*
      subroutine cc_r12offset(nr1orb,nr1xorb,nr1bas,nr1xbas,nr2bas,
     & nrgkl,nrxgkl,n2bst1,ir1orb,ir1xorb,ir1bas,ir1xbas,ir2bas,ir2xbas,
     & irgkl,irxgkl,iaodis1,nalphaj,ialphaj)
c----------------------------------------------------------------------
c     set some symmetry offsets and dimensions needed only when 
c     working in the R12 environment and the standard CC offset 
c     cannot be used because some dimensions were overwritten
c
c     H. Fliegl, C. Haettig spring 2003
c
c     adapted for more general R12-MOs:
c     C. Neiss, summer 2005
c----------------------------------------------------------------------
      implicit none 
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
      integer nr1orb(8),nr1xorb(8),nr1bas(8),nr1xbas(8),nr2bas,n2bst1(8)
      integer ir1xbas(8,8),ir1xorb(8,8)
      integer ir1orb(8,8),ir1bas(8,8),ir2bas(8,8),iaodis1(8,8)
      integer ir2xbas(8,8),irgkl(8,8),nrgkl(8)
      integer nrxgkl(8),irxgkl(8,8)
      integer isymdk,isymk,isymd,isym,isym1,isym2
      integer icount1,icount2,icount3,icount4,icount5,icount6
      integer icount7,icount8,icount9,icount10,icount11
      integer nalphaj(8),ialphaj(8,8)
      logical locdbg
      parameter (locdbg = .FALSE.)

      call qenter('cc_r12offset')
      if (locdbg) then
        write(LUPRI,*) 'Entered CC_R12OFFSET'
        call flshfo(LUPRI)
      end if

      do isymdk = 1, nsym
        nr1orb(isymdk)  = 0
        nr1xorb(isymdk) = 0
        nr1bas(isymdk)  = 0
        nr1xbas(isymdk) = 0
        n2bst1(isymdk)  = 0
        nalphaj(isymdk) = 0
        do isymk = 1, nsym
           isymd = muld2h(isymdk,isymk)
           nr1orb(isymdk) = nr1orb(isymdk) + norb1(isymd)*nrhfb(isymk)
           nr1xorb(isymdk)= nr1xorb(isymdk)+ norb2(isymd)*nrhfb(isymk)
           nr1bas(isymdk) = nr1bas(isymdk) + mbas1(isymd)*nrhfb(isymk)
           nr1xbas(isymdk)= nr1xbas(isymdk)+ mbas2(isymd)*nrhfb(isymk)
           n2bst1(isymdk) = n2bst1(isymdk) + mbas1(isymd)*mbas1(isymk)
           nalphaj(isymdk)= nalphaj(isymdk)+ mbas1(isymd)*nrhfa(isymk)
        end do
        if (locdbg) then
          write (lupri,*) 'nrhf, nrhfa, nrhfb:',
     &                     nrhf(isymdk), nrhfa(isymdk), nrhfb(isymdk)
          write (lupri,*) 'nr1orb, nr1xorb, nr1bas, nr1xbas: ',
     &                     nr1orb(isymdk),nr1xorb(isymdk),
     &                     nr1bas(isymdk),nr1xbas(isymdk),n2bst1(isymdk)
        end if
      end do

c     ------------------------------------
c     three index array for r12-integrals
c     ------------------------------------
      do isymdk = 1, nsym
        nrgkl(isymdk) = 0  
        nrxgkl(isymdk) = 0
        nvajkl(isymdk) = 0 
        do isymk = 1, nsym
          isymd = muld2h(isymdk,isymk)
          nrgkl(isymdk) = nrgkl(isymdk) + mbas1(isymd)*nmatkl(isymk) 
          nrxgkl(isymdk) = nrxgkl(isymdk) + mbas2(isymd)*nmatkl(isymk)
          nvajkl(isymdk) = nvajkl(isymdk) + nalphaj(isymd)*nmatkl(isymk)
         end do
      end do

c
c     do isym = 1, nsym
c       nvabkl(isym) = 0
c       icount1      = 0
c       do isym2 = 1, nsym
c         isym1 = muld2h(isym,isym2)
c         nvabkl(isym) = nvabkl(isym) + n2bst(isym1)*nmatkl(isym2)
c         ivabkl(isym1,isym2) = icount1
c         icount1 = icount1 + n2bst(isym1)*nmatkl(isym2)
c       end do
c     end do

c     ---------------------------------------
c     lenghts for nr1bas over all symmetries 
c     ---------------------------------------
      nr2bas = 0
      do isym = 1, nsym
         nr2bas = nr2bas + nr1bas(isym)*nr1bas(isym)
      end do

c     --------------------------------------------
c     calculate now the offsets for r12 integrals
c     --------------------------------------------
      do isymdk = 1, nsym
        icount1 = 0
        icount2 = 0
        icount3 = 0
        icount4 = 0
        icount5 = 0
        icount6 = 0
        icount7 = 0
        icount8 = 0
        icount9 = 0
        icount10 = 0
        icount11 = 0
        do isymk = 1, nsym
          isymd = muld2h(isymdk,isymk)
          ir1orb(isymd,isymk) = icount1 !occ+vir
          ir1xorb(isymd,isymk)= icount9 
          ir1bas(isymd,isymk) = icount2 !AO
          ir1xbas(isymd,isymk)= icount8
          ir2bas(isymd,isymk) = icount3 !AO+aux
          ir2xbas(isymd,isymk)= icount4 !aux
          irgkl(isymd,isymk)  = icount5
          irxgkl(isymd,isymk) = icount10
          iaodis1(isymd,isymk)= icount6
          ivajkl(isymd,isymk) = icount7
          ialphaj(isymd,isymk)= icount11
          icount1 = icount1 + nrhfb(isymk)*norb1(isymd)
          icount9 = icount9 + nrhfb(isymk)*norb2(isymd)
          icount2 = icount2 + nrhfb(isymk)*mbas1(isymd)
          icount8 = icount8 + nrhfb(isymk)*mbas2(isymd)
          icount3 = icount3 + nr1bas(isymd)*nr1bas(isymk)
          icount4 = icount4 + nr1bas(isymd)*nr1xbas(isymk)
          icount5 = icount5 + mbas1(isymd)*nmatkl(isymk) 
          icount6 = icount6 + mbas1(isymd)*mbas1(isymk) 
          icount7 = icount7 + nalphaj(isymd)*nmatkl(isymk) 
          icount10 = icount10 + mbas2(isymd)*nmatkl(isymk)
          icount11 = icount11 + nrhfa(isymk)*mbas1(isymd)
        end do
      end do
      if (locdbg) then
        write(LUPRI,*) 'Leaving CC_R12OFFSET'
        call flshfo(LUPRI)
      end if
      call qexit('cc_r12offset')
      return
      end 
*======================================================================*
