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
C
*=====================================================================*
      subroutine ccexnorm(isym,imult,threig,work,lwork)
*---------------------------------------------------------------------*
*
*    Purpose: normalize left and right eigenvectors such that
*
*                  (RE^f|RE^f) = 1
*                  (LE^i|RE^f) = delta_if
*
*    Written by Christof Haettig, November 2002, Friedrichstal
*
*=====================================================================*
      implicit none  
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccexci.h>
#include <ccexcinf.h>
#include <codata.h>
#include <dummy.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccnoddy.h>

      logical locdbg
      parameter (locdbg = .false.)
      integer mxovlp
      parameter (mxovlp = maxexci*maxexci)
      character fndeg*10
      parameter (fndeg='CCDEGEN_LE')
 
      integer isym, imult, lwork, istartend(maxexci,2)

#if defined (SYS_CRAY) 
      real threig
      real work(lwork), overlap(mxovlp)
      real ovlpslr(mxovlp), ovlpdlr(mxovlp), 
     &     ovlptlr(mxovlp), ovlpsrr(mxovlp), 
     &     ovlpdrr(mxovlp), ovlptrr(mxovlp)
      real seigr(maxexci), seigi(maxexci)
      real srr, drr, trr, xnrmrr, slr, dlr, tlr, xnrmlr
      real zero, one, xhund, tlovlp, xnrm, fac, xerror
      real tiny
#else 
      double precision threig
      double precision work(lwork), overlap(mxovlp)
      double precision ovlpslr(mxovlp), ovlpdlr(mxovlp), 
     &                 ovlptlr(mxovlp), ovlpsrr(mxovlp), 
     &                 ovlpdrr(mxovlp), ovlptrr(mxovlp)
      double precision seigr(maxexci), seigi(maxexci)
      double precision srr, drr, trr, xnrmrr, slr, dlr, tlr, xnrmlr
      double precision zero, one, xhund, tlovlp, xnrm, fac, xerror
      double precision tiny, xnrmr, xnrml
#endif 
      parameter ( one = 1.0d0, xhund = 100.0d0, tlovlp = 0.5d0,
     &            zero = 0.0d0, tiny = 1.0d-10 )

      character listl*3, listr*3, model*10
      logical triplet, degeneracy
      integer novlp, ilstnrst, ilstnrnd, ileft, ncctot, iadrij,
     &        kl1am, kl2am, kr1am, kr2am, kwrk1, lwrk1, iopt, ipass,
     &        iright, iovlp, istate,
     &        nccvar, nccvars, nccvard, nst, ioff, k1amp, k2amp,
     &        ilstnr, nleft, idxr, idxl, ndeg, matz, ierr, ludeg

      integer ilstsym
#if defined (SYS_CRAY) 
      real ddot
#else
      double precision ddot
#endif 

      call qenter('ccexnorm')

*---------------------------------------------------------------------*
* some initializations:
*---------------------------------------------------------------------*
      if      (imult.eq.1) then
         triplet = .false.
      else if (imult.eq.3) then
         triplet = .true.
      else
         call quit('Illegal multiplicity in CCEXNORM.')
      end if

      if (triplet .and. nccexci(isym,imult).gt.0) then
         call quit('Triplet not yet implemented in CCEXNORM.')
      end if

      ! set index of first and last vector in symmetry class
      if (triplet) then
        ilstnrst = itrofe(isym) + 1
        ilstnrnd = itrofe(isym) + nccexci(isym,imult)
      else
        ilstnrst = isyofe(isym) + 1
        ilstnrnd = isyofe(isym) + nccexci(isym,imult)
      end if

      ! work space allocation for two singles (and doubles) vectors
      nccvars = nt1am(isym)
      nccvard = nt2am(isym)
      if (triplet) nccvard = nccvard + nt2ama(isym)
      ncctot  = nccvars 
      if (.not.ccs) ncctot = ncctot + nccvard

      kl1am = 1
      kwrk1 = kl1am + nccvars
      if (.not.ccs) then
        kl2am = kl1am + nccvars
        kwrk1 = kl2am + nccvard
      end if

      kr1am = kwrk1
      kwrk1 = kr1am + nccvars
      if (.not.ccs) then
        kr2am = kr1am + nccvars
        kwrk1 = kr2am + nccvard
      end if

      lwrk1 = lwork - kwrk1
      if (lwrk1.lt.0) call quit('out of memory in CCEXNORM (1)')

*---------------------------------------------------------------------*
* check for degeneracies within this symmetry class and try to get
* the correct matching of left and right eigenvectors
*---------------------------------------------------------------------*
      ileft = ilstnrst
      do while (ileft .lt. ilstnrnd)

        ! check if state is degenerate
        ndeg = 1 
        do while (abs(eigval(ileft)-eigval(ileft+ndeg)).lt.threig) 
           ndeg = ndeg + 1
        end do

        if (ndeg.gt.1) then
          write(lupri,'(//1x,a)') 'Found a degeneracy between:'
          write(lupri,'(1x,a)')   '---------------------------'
          write(lupri,'(1x,a)')   'State  Excitation Energy'

          ludeg = -1
          call wopen2(ludeg,fndeg,64,0) 

c         ---------------------------------------------------
c         set up (LE|RE) overlap matrix, which is used as the
c         transformation matrix:
c         ---------------------------------------------------
          iadrij = 0
          do idxl = ileft, ileft-1+ndeg
            write(lupri,'(1x,i3,f15.8)') idxl,eigval(idxl)

            iopt = 3
            call cc_rdrsp('LE ',idxl,isym,iopt,model,
     &                    work(kl1am),work(kl2am))
            xnrml = dsqrt(ddot(ncctot,work(kl1am),1,work(kl1am),1))

            iadrij = (idxl-ileft)*ncctot + 1
            call putwa2(ludeg,fndeg,work(kl1am),iadrij,ncctot)

            xnrm = zero
            do idxr = ileft, ileft-1+ndeg
              iopt  = 3
              call cc_rdrsp('RE ',idxr,isym,iopt,model,
     &                      work(kr1am),work(kr2am))
              xnrmr = dsqrt(ddot(ncctot,work(kr1am),1,work(kr1am),1))
              iovlp = (idxl-ileft)*ndeg+(idxr-ileft+1)
              ovlpslr(iovlp)=ddot(ncctot,work(kl1am),1,work(kr1am),1)/
     &                         (xnrml*xnrmr)
              xnrm = xnrm + ovlpslr(iovlp)*ovlpslr(iovlp)
            end do
            fac   = one/dsqrt(xnrm)
            iovlp = (idxl-ileft)*ndeg+1
            call dscal(ndeg,fac,ovlpslr(iovlp),1)
          end do

          ! check if matrix is orthogonal:
          call dgemm('n','t',ndeg,ndeg,ndeg,
     &               one,ovlpslr,ndeg,ovlpslr,ndeg,
     &               zero,ovlpdlr,ndeg)
          xerror = zero
          do idxl = ileft, ileft-1+ndeg
            do idxr = ileft, ileft-1+ndeg
              iovlp = (idxl-ileft)*ndeg+(idxr-ileft+1)
              if (idxl.eq.idxr) xerror=xerror+abs(ovlpdlr(iovlp)-one)
              if (idxl.ne.idxr) xerror=xerror+abs(ovlpdlr(iovlp)-zero)
            end do
          end do

          ! some print out
          write(lupri,'(//,1x,a)') 
     &      'transformation matrix computed from (LE|RE) overlap:'
          call output(ovlpslr,1,ndeg,1,ndeg,ndeg,ndeg,1,lupri)
          write(lupri,'(/,1x,a,f20.10)') 
     &      'deviation from orthogonality:',xerror

c         -----------------------------------------------------------
c         transform the left eigenvector to match right eigenvectors:
c                     LE(j)' = sum_i LE(i) (LE(i)|RE(j))
c         -----------------------------------------------------------
          do idxl = ileft, ileft-1+ndeg
            call dzero(work(kl1am),ncctot)
            do idxr = ileft, ileft-1+ndeg
              iovlp = (idxr-ileft)*ndeg+(idxl-ileft+1)
              if (abs(ovlpslr(iovlp)).gt.tiny) then
                iadrij = (idxr-ileft)*ncctot + 1
                call getwa2(ludeg,fndeg,work(kr1am),iadrij,ncctot)
                call daxpy(ncctot,ovlpslr(iovlp),
     &                     work(kr1am),1,work(kl1am),1)
              end if
            end do
            iopt = 3
            call cc_wrrsp('LE ',idxl,isym,iopt,model,dummy,
     &                    work(kl1am),work(kl2am),work(kwrk1),lwrk1)
          end do

          call wclose2(ludeg,fndeg,'DELETE')

        end if

        ileft = ileft + ndeg
      end do
      
*---------------------------------------------------------------------*
* build list of overlap matrix elements to be computed:
*---------------------------------------------------------------------*
      novlp = 0
      do ileft = ilstnrst, ilstnrnd
        if (sqrovlp) then
          istartend(ileft,1) = ilstnrst
          istartend(ileft,2) = ilstnrnd
          novlp              = novlp + nccexci(isym,imult)
        else
          istartend(ileft,1) = ileft
          istartend(ileft,2) = ileft
          novlp              = novlp + 1
        end if
      end do
    
      if (locdbg) then
        write(lupri,*) 'isym,imult,nccexci,novlp:',
     &                  isym,imult,nccexci(isym,imult),novlp
      end if

*     ----------------------------------------------------------------
*     compute overlap matrices: singles and doubles part is done here,
*     for the triples part we call ccovlpt(_noddy)
*     ----------------------------------------------------------------
      do ipass = 1, 2
       if (ipass.eq.1) then
         listl = 'RE '
         listr = 'RE '
       else
         listl = 'LE '
         listr = 'RE '
       end if

       call dzero(ovlpslr,novlp)
       call dzero(ovlpdlr,novlp)
       call dzero(ovlptlr,novlp)

       iovlp = 0
       do ileft = ilstnrst, ilstnrnd
        iopt    = 3
        call cc_rdrsp(listl,ileft,isym,iopt,model,
     &                work(kl1am),work(kl2am))

        do iright = istartend(ileft,1), istartend(ileft,2)
          iovlp = iovlp + 1
          iopt  = 3
          call cc_rdrsp(listr,iright,isym,iopt,model,
     &                  work(kr1am),work(kr2am))
          ovlpslr(iovlp) = ddot(nccvars,work(kl1am),1,work(kr1am),1)
          if (.not.ccs) ovlpdlr(iovlp) = 
     &                     ddot(nccvard,work(kl2am),1,work(kr2am),1)
          if (locdbg) then
            write(lupri,*) 'listl,listr,ileft,iright,iovlp:',
     &                      listl,listr,ileft,iright,iovlp
            write(lupri,*) 'norm^2 (S|S), (D|D):',
     &          ovlpslr(iovlp),ovlpdlr(iovlp)
          end if
        end do
       end do
       if (iovlp.ne.novlp) call quit('error in CCEXNORM.')
  
       if (ccsdt) then
         if (ipass.eq.1 .and. (.not.ccsdtrenrm)) then
           continue 
         else
           if (locdbg) write(lupri,*)'entering overlap calc'
           IF (NODDY_OVLP) THEN 
              if (locdbg) write(lupri,*)'noddy overlap'
              call ccovlp_noddy(listl,listr,1,ovlptlr,
     &                          ilstnrst,ilstnrnd,istartend,
     &                          1,1,1,1,dummy,work,lwork)
           ELSE 
              if (listl.eq.'RE ') call quit('error in CCEXNORM: '//
     &           '(RE|RE) overlap for not yet available for triples')
              if (locdbg) write(lupri,*)'non-noddy overlap'
              call cc3_lr_ovlp(listl,listr,1,ovlptlr,
     &                         ilstnrst,ilstnrnd,istartend,
     &                         1,1,1,1,dummy,work,lwork)
           END IF
         end if
       end if

       if (ipass.eq.1) then
         call dcopy(novlp,ovlpslr,1,ovlpsrr,1)
         call dcopy(novlp,ovlpdlr,1,ovlpdrr,1)
         call dcopy(novlp,ovlptlr,1,ovlptrr,1)
       end if

      end do

*     --------------------------------
*     print computed overlap matrices:
*     --------------------------------
      if (sqrovlp) then
        do ipass = 1, 2
         write(lupri,*)
         if (ipass.eq.1) then
           write(lupri,*) 'Overlap of right eigenvectors (RE|RE):'
           write(lupri,*) '======================================'
           call dcopy(novlp,    ovlpsrr,1,overlap,1)
           call daxpy(novlp,one,ovlpdrr,1,overlap,1)
           call daxpy(novlp,one,ovlptrr,1,overlap,1)
         else
           write(lupri,*) 'Overlap of eigenvectors (LE|RE):'
           write(lupri,*) '================================'
           call dcopy(novlp,    ovlpslr,1,overlap,1)
           call daxpy(novlp,one,ovlpdlr,1,overlap,1)
           call daxpy(novlp,one,ovlptlr,1,overlap,1)
         end if

         iovlp = 1
         nst = nccexci(isym,imult)
         write(lupri,*) 'overlap matrices for isym,imult:',isym,imult
         write(lupri,*) 'single excitation contribution:'
         if (ipass.eq.1) then 
           call output(ovlpsrr(iovlp),1,nst,1,nst,nst,nst,1,lupri)
         else
           call output(ovlpslr(iovlp),1,nst,1,nst,nst,nst,1,lupri)
         end if
        
         if (.not.ccs) then
          write(lupri,*) 'double excitation contribution:'
          if (ipass.eq.1) then 
            call output(ovlpdrr(iovlp),1,nst,1,nst,nst,nst,1,lupri)
          else
            call output(ovlpdlr(iovlp),1,nst,1,nst,nst,nst,1,lupri)
          end if
         end if
        
         if (ccsdt) then
          write(lupri,*) 'triple excitation contribution:'
          if (ipass.eq.1) then 
            call output(ovlptrr(iovlp),1,nst,1,nst,nst,nst,1,lupri)
          else
            call output(ovlptlr(iovlp),1,nst,1,nst,nst,nst,1,lupri)
          end if
         end if

         write(lupri,*) 'total overlap matrix:'
         call output(overlap(iovlp),1,nst,1,nst,nst,nst,1,lupri)

         iovlp = iovlp + nst*nst
        end do
        write(lupri,*)
      end if

*     ----------------------------------------------------------------
*     print short summary of contributions to norm and renormalize:
*     ----------------------------------------------------------------
      do ipass = 1, 2

       if (ipass.eq.1) then
         write(lupri,'(///,"+",70("="),"+")')
         write(lupri,'(a,/,a,/,"+",70("-"),"+",/,a)')
     &  '| contributions of excitation levels to right eigenvectors'//
     &      ' as computed |',
     &  '| from the (RE|RE) self overlap                           '//
     &      '             |',
     &  '| sym | multi | state |    E / eV  | % singles |'//
     &  ' % doubles | % triples |'
         write(lupri,'("+",70("="),"+")')
       else
         write(lupri,'(///,"+",70("="),"+")')
         write(lupri,'(a,/,a,/,"+",70("-"),"+",/,a)')
     &  '| contributions of excitation levels to excited states as '//
     &      'computed     |',
     &  '| from the (LE|RE) overlap                                '//
     &      '             |',
     &  '| sym | multi | state |    E / eV  | % singles |'//
     &  ' % doubles | % triples |'
         write(lupri,'("+",70("="),"+")')
       end if

       ioff = 0
       do istate = 1, nccexci(isym,imult)

        nccvar  = nt1am(isym)
        nccvard = nt2am(isym)
        if (imult.eq.3) nccvard = nccvard + nt2ama(isym)
        if (.not.(ccs.or.cis)) nccvar  = nccvar + nccvard
        
        k1amp = 1
        kwrk1 = k1amp + nt1am(isym)
        if (.not.ccs) then
          k2amp = k1amp + nt1am(isym)
          kwrk1 = k2amp + nccvard
        end if
        lwrk1 = lwork - kwrk1
        if (lwrk1.lt.0) call quit('out of memory in CCEXNORM (2)')

        if (sqrovlp) then
          iovlp = ioff + nccexci(isym,imult)*(istate-1)+istate
        else
          iovlp = ioff + istate 
        end if
       
        if (imult.eq.3) then
          ilstnr = itrofe(isym) + istate
        else
          ilstnr = isyofe(isym) + istate
        end if
          
        srr = ovlpsrr(iovlp)
        drr = ovlpdrr(iovlp)
        trr = ovlptrr(iovlp)
        xnrmrr = srr + drr + trr

        slr = ovlpslr(iovlp)
        dlr = ovlpdlr(iovlp)
        tlr = ovlptlr(iovlp)
        xnrmlr = slr + dlr + tlr
       
        if (ipass.eq.1) then
          write(lupri,'("| ",i4,"|",2(i4,3x,"|"),f11.5,'//
     &                '" |",3(1x,f8.4,2x,"|"))')
     &     isym,imult,istate,eigval(ilstnr)*xtev,
     &     xhund*srr/xnrmrr,xhund*drr/xnrmrr,xhund*trr/xnrmrr
        else
          write(lupri,'("| ",i4,"|",2(i4,3x,"|"),f11.5,'//
     &                '" |",3(1x,f8.4,2x,"|"))')
     &     isym,imult,istate,eigval(ilstnr)*xtev,
     &     xhund*slr/xnrmlr,xhund*dlr/xnrmlr,xhund*tlr/xnrmlr

c         ---------------------------------------------------
c         renormalize right eigenvectors such that (RE|RE)=1:
c         ---------------------------------------------------
          iopt = 3
          call cc_rdrsp('RE ',ilstnr,isym,iopt,model,
     &                  work(k1amp),work(k2amp))
          
          call dscal(nccvar,one/dsqrt(xnrmrr),work(k1amp),1)
          
          iopt = 3
          call cc_wrrsp('RE ',ilstnr,isym,iopt,model,dummy,
     &                  work(k1amp),work(k2amp),work(kwrk1),lwrk1)
          
c         --------------------------------------------------
c         renormalize left eigenvectors such that (LE|RE)=1:
c         (but first check if this  can be done)
c         --------------------------------------------------
          if (abs(xnrmlr).lt.tlovlp) then
           write(lupri,*)'SERIOUS WARNING !!!!!!!!!!!!!!!!!!!!'
           write(lupri,*)'Overlap is less than tolerance.'
           write(lupri,*)' - You may encounter serious problems.'
           write(lupri,*)'Symmetry    : ',isym
           write(lupri,*)'Multiplicity: ',imult
           write(lupri,*)'State       : ',istate
           write(lupri,*)'Overlap     : ',xnrmlr
           write(lupri,*)'Tolerance   : ',tlovlp
           write(lupri,*)'Left eigenvector will not be normalized!!'
           write(lupri,*)'(Transition) properties will be wrong!!!!'
          else
            iopt = 3
            call cc_rdrsp('LE ',ilstnr,isym,iopt,model,
     &                    work(k1amp),work(k2amp))
            
            call dscal(nccvar,dsqrt(xnrmrr)/xnrmlr,work(k1amp),1)
            
            iopt = 3
            call cc_wrrsp('LE ',ilstnr,isym,iopt,model,dummy,
     &                    work(k1amp),work(k2amp),work(kwrk1),lwrk1)
C FIX OF CCSDR(3) BUG. SAVE UNIT NORM IN APPROPRIATE VECTOR
            xnorm(ilstnr) = one
          endif

        end if
       end do

       write(lupri,'("+",70("="),"+",//)')

      end do

      call qexit('ccexnorm')

      return
      end

*---------------------------------------------------------------------*
*              end of subroutine ccexnorm                             *
*---------------------------------------------------------------------*
*=====================================================================*
      subroutine ccovlp_noddy(listl,listr,imode,
     &                        overlap,ilstnrst,ilstnrnd,istartend,
     &                        ndtran,mxvec,idtran,iddots,dotcon,
     &                        work,lwork)
*---------------------------------------------------------------------*
*
*    Purpose: compute overlap (L|R) between two triples vectors
*
*             imode=1  use ilstnrst,ilstnrnd, and istartend to
*                      set up the loops over left and right vectors
*                      result is returned in overlap
*
*             imode=2  use ndtran,idtran,iddots to set up the loops
*                      over left and right vectors and return result
*                      in dotcon
*
*    Written by Christof Haettig, November 2002, Friedrichstal
*
*=====================================================================*
      implicit none  
#include <priunit.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccfield.h>
#include <ccorb.h>
#include <inftap.h> 
#include <dummy.h> 
#include <ccsdsym.h>
#include <ccexci.h>

      logical locdbg
      parameter (locdbg=.false.)

      integer isym0
      parameter (isym0 = 1)
 
      character listl*3, listr*3
      integer lwork, ilstnrst,ilstnrnd, ndtran, mxvec, imode
      integer istartend(maxexci,2)
      integer idtran(ndtran), iddots(mxvec,ndtran)

#if defined (SYS_CRAY) 
      real work(lwork), overlap(*)
      real tcon, sixth, freql, freqr, xcon(3)
      real dotcon(mxvec,ndtran)
#else 
      double precision work(lwork), overlap(*)
      double precision tcon, sixth, freql, freqr, xcon(3)
      double precision dotcon(mxvec,ndtran)
#endif 
      parameter( sixth = 1.0d0/6.0d0 ) 

      character model*10
      logical do_intxy, do_intt0, do_ints0
      integer klamp0, klamh0, kfock0, kfockd, kt1amp0, kl3am, kr3am,
     &        kxiajb, kyiajb, kint1t, kint2t, kwrk1, lwrk1, nt3amp,
     &        ileft, isyml, imull, iright, isymr, imulr, iopt,
     &        kscr1,kint1sc,kint2sc,klampc,klamhc,kfockc,kwrk2,lwrk2,
     &        lufock,iovlp, kint1s, kint2s, isymd, illl, idel, isydis,
     &        kxint, idxai, idxbj, idxck, idiag, idxt3, kdum,
     &        istart, iend, istartr, iendr, ivec, itran


      integer ilstsym
#if defined (SYS_CRAY) 
      real ddot
#else
      double precision ddot
#endif 

      kdum = 1
      call qenter('ccov-nod')

*---------------------------------------------------------------------*
*     Initializations and work space allocation:
*---------------------------------------------------------------------*
      nt3amp = nt1amx * nt1amx * nt1amx

      ! check if XIAJB and/or YIAJB integrals needed
      if (listl.eq.'LE ' .or. listl.eq.'L1 ' .or.
     &    listr.eq.'LE ' .or. listr.eq.'L1 ' .or.        
     &    listl.eq.'X1B' .or. listr.eq.'X1B'     ) then
        do_intxy = .true.
      end if

      ! check if XINT1T and XINT2T integrals needed
      if (listl.eq.'LE ' .or. listl.eq.'L1 ' .or.
     &    listr.eq.'LE ' .or. listr.eq.'L1 ' .or.       
     &    listl.eq.'X1B' .or. listr.eq.'X1B'     ) then
        do_intt0 = .true.
      end if

      ! check if XINT1S and XINT2S integrals needed
      if (listl.eq.'RE ' .or. listl.eq.'R1 ' .or.
     &    listr.eq.'RE ' .or. listr.eq.'R1 ' .or.
     &    listl.eq.'ER1' .or. listl.eq.'R2 ' .or.
     &    listr.eq.'ER1' .or. listr.eq.'R2 ' .or.
     &    listl.eq.'EO1' .or. listl.eq.'O2 ' .or.
     &    listr.eq.'EO1' .or. listr.eq.'O2 '     ) then
        do_ints0 = .true.
      end if

      klamp0 = 1
      klamh0 = klamp0 + nlamdt
      kfock0 = klamh0 + nlamdt
      kfockd = kfock0 + norbt*norbt
      kt1amp0= kfockd + nt1amx
      kwrk1  = kt1amp0+ nt1amx

      kl3am  = kwrk1
      kwrk1  = kl3am  + nt3amp

      kr3am  = kwrk1
      kwrk1  = kr3am  + nt3amp

      if (do_intxy) then
        kxiajb = kwrk1 
        kyiajb = kxiajb + nt1amx*nt1amx 
        kwrk1  = kyiajb + nt1amx*nt1amx
      end if
        
      if (do_intt0) then
        kint1t = kwrk1  
        kint2t = kint1t + nt1amx*nvirt*nvirt
        kwrk1  = kint2t + nt1amx*nrhft*nrhft
      end if

      if (do_ints0) then
        kint1s = kwrk1  
        kint2s = kint1s + nt1amx*nvirt*nvirt
        kwrk1  = kint2s + nt1amx*nrhft*nrhft
      end if

      lwrk1 = lwork - kwrk1
      if (lwrk1 .lt. 0) then
         call quit('Insufficient space in ccovlp_noddy')
      endif 

*---------------------------------------------------------------------*
*     Read SCF orbital energies from file:
*---------------------------------------------------------------------*
      call gpopen(lusifc,'SIRIFC','OLD',' ','UNFORMATTED',idummy,
     &            .false.)
      rewind lusifc
      call mollab('TRCCINT ',lusifc,lupri)
      read (lusifc)
      read (lusifc) (work(kfockd+i-1), i=1,norbt)
      call gpclose(lusifc,'KEEP')
 
*---------------------------------------------------------------------*
*     Get zeroth-order Lambda matrices:
*---------------------------------------------------------------------*
      iopt   = 1
      call cc_rdrsp('R0',0,isym0,iopt,model,work(kt1amp0),dummy)
 
      call lammat(work(klamp0),work(klamh0),work(kt1amp0),
     &            work(kwrk1),lwrk1)
 
*---------------------------------------------------------------------*
*     Read zeroth-order AO Fock matrix from file:
*---------------------------------------------------------------------*
      lufock = -1
      call gpopen(lufock,'CC_FCKH','OLD',' ','UNFORMATTED',idummy,
     &            .false.)
      rewind(lufock)
      read(lufock) (work(kfock0-1+i),i=1,n2bst(isym0))
      call gpclose(lufock,'KEEP')
 
      call cc_fckmo(work(kfock0),work(klamp0),work(klamh0),
     &              work(kwrk1),lwrk1,isym0,isym0,isym0)
                                                          
*---------------------------------------------------------------------*
*     Precompute integrals which do not depend on eigenvectors:
*---------------------------------------------------------------------*
      if (do_intxy) then
         call dzero(work(kxiajb),nt1amx*nt1amx)
         call dzero(work(kyiajb),nt1amx*nt1amx)
      end if
      if (do_intt0) then
         call dzero(work(kint1t),nt1amx*nvirt*nvirt)
         call dzero(work(kint2t),nt1amx*nrhft*nrhft)
      end if
      if (do_ints0) then
        call dzero(work(kint1s),nt1amx*nvirt*nvirt)
        call dzero(work(kint2s),nt1amx*nrhft*nrhft)
      end if

      do isymd = 1, nsym
         do illl = 1,nbas(isymd)
            idel   = ibas(isymd) + illl
            isydis = muld2h(isymd,isym0)
 
C           ----------------------------
C           Work space allocation no. 2.
C           ----------------------------
            kxint  = kwrk1
            kwrk2  = kxint + ndisao(isydis)
            lwrk2  = lwork - kwrk2
            if (lwrk2 .lt. 0) then
               write(lupri,*) 'Need : ',kwrk2,'Available : ',lwork
               call quit('Insufficient space in ccovlp_noddy')
            endif
 
C           ---------------------------
C           Read in batch of integrals.
C           ---------------------------
            call ccrdao(work(kxint),idel,1,work(kwrk2),lwrk2,
     *                  work(kdum),direct)
 
C           ----------------------------------
C           Calculate integrals needed in CC3:
C           ----------------------------------
            if (do_intxy) then
               call cc3_tran2(work(kxiajb),work(kyiajb),
     &                        work(klamp0),work(klamh0),
     &                        work(kxint),idel)
            end if
 
            if (do_intt0) then
               call ccsdt_tran1(work(kint1t),work(kint2t),
     &                          work(klamp0),work(klamh0),
     &                          work(kxint),idel)
            end if

            if (do_ints0) then
              call ccsdt_tran3(work(kint1s),work(kint2s),
     &                         work(klamp0),work(klamh0),
     &                         work(kxint),idel)
            end if
 
         end do   
      end do  

      if (locdbg) then
       if (do_intxy) write(lupri,*) 'norm^2(xiajb):',
     &        ddot(nt1amx*nt1amx,work(kxiajb),1,work(kxiajb),1)
       if (do_intxy) write(lupri,*) 'norm^2(yiajb):',
     &        ddot(nt1amx*nt1amx,work(kyiajb),1,work(kyiajb),1)
       if (do_intt0) write(lupri,*) 'norm^2(xint1t):',
     &        ddot(nt1amx*nvirt*nvirt,work(kint1t),1,work(kint1t),1)
       if (do_intt0) write(lupri,*) 'norm^2(xint2t):',
     &        ddot(nt1amx*nrhft*nrhft,work(kint2t),1,work(kint2t),1)
       if (do_ints0) write(lupri,*) 'norm^2(xint1s):',
     &        ddot(nt1amx*nvirt*nvirt,work(kint1s),1,work(kint1s),1)
       if (do_ints0) write(lupri,*) 'norm^2(xint2s):',
     &        ddot(nt1amx*nrhft*nrhft,work(kint2s),1,work(kint2s),1)
      end if

*---------------------------------------------------------------------*
*     Loop over "left" vectors, pre-compute them and calculate
*     overlap matrix contributions in loop over "right" vectors
*---------------------------------------------------------------------*
      iovlp = 0

      if      (imode.eq.1) then
        istart = ilstnrst
        iend   = ilstnrnd
      else if (imode.eq.2) then
        istart = 1
        iend   = ndtran
      else
        call quit('illegal imode in ccovlp_noddy.')
      end if

      do itran = istart, iend
        if      (imode.eq.1) then
           ileft = itran
        else if (imode.eq.2) then
           ileft = idtran(itran)
        else
          call quit('illegal imode in ccovlp_noddy.')
        end if

        isyml   = ilstsym(listl,ileft)
        imull   = imulte(ileft)

        if (locdbg) then
          write(lupri,*) 'ileft,isyml,imull:',ileft,isyml,imull
          write(lupri,*) 'istart,iend:',istartend(ileft,1), 
     &                                  istartend(ileft,2)
        end if

        if      (listl.eq.'LE ' .or. listl.eq.'L1 ' .or.
     &           listl.eq.'X1B' .or. listl.eq.'BE ' .or.
     &           listl.eq.'N2 ' .or. listl.eq.'M1 '     ) then
          kscr1 = kwrk1
          kwrk2 = kscr1 + nt1amx
          lwrk2 = lwork - kwrk2

          if (locdbg) then
            write(lupri,*) 'call ccsdt_tbar31_noddy for:',listl,ileft
          end if

          call ccsdt_tbar31_noddy(work(kl3am),freql,listl,ileft,
     &                          work(klamp0),work(klamh0),
     &                          work(kfock0),work(kfockd),work(kscr1),
     &                          work(kxiajb),work(kint1t),work(kint2t),
     &                          work(kwrk2),lwrk2)  


        else if (listl.eq.'RE ' .or. listl.eq.'R1 ' .or.
     &           listl.eq.'RC '                         ) then
          kint1sc = kwrk1
          kint2sc = kint1sc + nt1amx*nvirt*nvirt
          klampc  = kint2sc + nt1amx*nrhft*nrhft
          klamhc  = klampc  + nlamdt
          kfockc  = klamhc  + nlamdt
          kwrk2   = klamhc  + norbt*norbt
          lwrk2   = lwork - kwrk2

          if (locdbg) then
            write(lupri,*) 'call ccsdt_t31_noddy for:',listl,ileft
          end if

          call ccsdt_t31_noddy(work(kl3am),listl,ileft,
     &                         freql,.false.,
     &                         .false.,work(kint1s),work(kint2s),
     &                         .false.,dummy,dummy,
     &                         .false.,dummy,dummy,
     &                         work(kint1sc),work(kint2sc),
     &                         work(klampc),work(klamhc),work(kfockc),
     &                         work(klamp0),work(klamh0),work(kfock0),
     &                         dummy,work(kfockd),
     &                         work(kwrk2),lwrk2)

        else if (listl.eq.'O2 ' .or. listl.eq.'R2 ' .or.
     &           listl.eq.'EO1' .or. listl.eq.'ER1'     ) then
          kscr1 = kwrk1
          kwrk2 = kscr1 + nt1amx
          lwrk2 = lwork - kwrk2

          if (locdbg) then
            write(lupri,*) 'call ccsdt_t32_noddy for:',listl,ileft
          end if

          call ccsdt_t32_noddy(work(kl3am),listl,ileft,freql,
     &                         work(kint1s),work(kint2s),
     &                         work(klamp0),work(klamh0),work(kfock0),
     &                         work(kfockd),
     &                         dummy,dummy,
     &                         work(kscr1),work(kwrk2),lwrk2)

        end if


        if      (imode.eq.1) then
          istartr = istartend(ileft,1)
          iendr   = istartend(ileft,2)
        else if (imode.eq.2) then
          istartr = 1
          iendr   = 1
          do while (iddots(iendr,itran).ne.0 .and. iendr.le.mxvec)
            iendr = iendr + 1
          end do
          iendr = iendr - 1
        else
          call quit('illegal imode in ccovlp_noddy.')
        end if

        do ivec = istartr, iendr
          if      (imode.eq.1) then
             iright = ivec
          else if (imode.eq.2) then
             iright = iddots(ivec,itran)
          else
            call quit('illegal imode in ccovlp_noddy.')
          end if
          
          isymr   = ilstsym(listr,iright)
          imulr   = imulte(iright)

          if (locdbg) then
            write(lupri,*) 'iright,isymr,imulr:',iright,isymr,imulr
          end if

          if (isymr.eq.isyml .and. imulr.eq.imull) then

            if      (listr.eq.'LE ' .or. listr.eq.'L1 ' .or.
     &               listr.eq.'X1B' .or. listr.eq.'BE ' .or.
     &               listr.eq.'N2 ' .or. listr.eq.'M1 '     ) then
              kscr1 = kwrk1
              kwrk2 = kscr1 + nt1amx
              lwrk2 = lwork - kwrk2
            
              if (locdbg) then
               write(lupri,*)'call ccsdt_tbar31_noddy for:',listr,iright
              end if

              call ccsdt_tbar31_noddy(work(kr3am),freqr,listr,iright,
     &                          work(klamp0),work(klamh0),
     &                          work(kfock0),work(kfockd),work(kscr1),
     &                          work(kxiajb),work(kint1t),work(kint2t),
     &                          work(kwrk2),lwrk2)  
            
              do_intxy = .false.
              do_intt0 = .false.
            
            else if (listr(1:2).eq.'RE ' .or. listr(1:3).eq.'R1 ' .or.
     &               listr(1:2).eq.'RC '                         ) then
              kint1sc = kwrk1
              kint2sc = kint1sc + nt1amx*nvirt*nvirt
              klampc  = kint2sc + nt1amx*nrhft*nrhft
              klamhc  = klampc  + nlamdt
              kfockc  = klamhc  + nlamdt
              kwrk2   = klamhc  + norbt*norbt
              lwrk2   = lwork - kwrk2
            
              if (locdbg) then
                write(lupri,*) 'call ccsdt_t31_noddy for:',listr,iright
              end if

              call ccsdt_t31_noddy(work(kr3am),listr,iright,
     &                          freqr,.false.,
     &                          .false.,work(kint1s),work(kint2s),
     &                          .false.,dummy,dummy,
     &                          .false.,dummy,dummy,
     &                                  work(kint1sc),work(kint2sc),
     &                          work(klampc),work(klamhc),work(kfockc),
     &                          work(klamp0),work(klamh0),work(kfock0),
     &                          dummy,work(kfockd),
     &                          work(kwrk2),lwrk2)
            
            else if (listr.eq.'O2 ' .or. listr.eq.'R2 ' .or.
     &               listr.eq.'EO1' .or. listr.eq.'ER1'     ) then

              kscr1 = kwrk1
              kwrk2 = kscr1 + nt1amx
              lwrk2 = lwork - kwrk2
            
              if (locdbg) then
                write(lupri,*) 'call ccsdt_t32_noddy for:',listr,iright
              end if
            
              call ccsdt_t32_noddy(work(kr3am),listr,iright,freqr,
     &                          work(kint1s),work(kint2s),
     &                          work(klamp0),work(klamh0),work(kfock0),
     &                          work(kfockd),dummy,dummy,
     &                          work(kscr1),work(kwrk2),lwrk2)

            end if

            tcon = - sixth * ddot(nt3amp,work(kl3am),1,work(kr3am),1)
            
            if      (imode.eq.1) then
              iovlp = iovlp + 1
              overlap(iovlp) = overlap(iovlp) + tcon
            else if (imode.eq.2) then
              dotcon(ivec,itran) = dotcon(ivec,itran) + tcon
            else
              call quit('illegal imode in ccovlp_noddy.')
            end if

            if (locdbg) then
              write(lupri,*) 'iovlp:',iovlp
              write(lupri,*) 'xcon:',xcon
              write(lupri,*) 'ccovlp_noddy> norm^2(R3):',
     &          sixth*ddot(nt3amp,work(kr3am),1,work(kr3am),1)
              write(lupri,*) 'ccovlp_noddy> norm^2(L3):',
     &          sixth*ddot(nt3amp,work(kl3am),1,work(kl3am),1)
              if (imode.eq.1)
     &          write(lupri,*) 'tcon,overlap:',tcon,overlap(iovlp)
              if (imode.eq.2)
     &          write(lupri,*) 'tcon,dotcon:',tcon,dotcon(ivec,itran)
            end if

          end if
           
        end do
      end do
  
      call qexit('ccov-nod')

      return
      end

*---------------------------------------------------------------------*
*              end of subroutine ccovlp_noddy                         *
*---------------------------------------------------------------------*
C  /* Deck cc3_lr_ovlp */
      subroutine cc3_lr_ovlp(listl,listr,imode,overlap,
     &                       ilstnrst,ilstnrnd,istartend,
     &                       ndtran,mxvec,idtran,iddots,dotcon,
     &                       work,lwork)     
*---------------------------------------------------------------------*
*
*    Purpose: compute overlap (L|R) = (EBAR3!E3) between triples part 
*             of left and right eigenvectors 
*
*=====================================================================*
C
      IMPLICIT NONE
C
#include <ccexci.h>
#include <priunit.h>
C
      INTEGER ilstnrst,ilstnrnd,istartend(maxexci,2),lwork
      INTEGER LU3VI,LU3VI2,LU3FOP,LU3FOP2,LUCKJD,LUTOC,LUDKBC
      INTEGER iovlp,ileft,isyml,imull,iright,isymr,imulr,
     &        imode,ndtran,mxvec
      integer idtran(ndtran), iddots(mxvec,ndtran)
C
      CHARACTER*(*) listl,listr
      CHARACTER*8 FN3VI2,FNTOC
      CHARACTER*6 FNCKJD, FN3VI, FN3FOP2
      CHARACTER*5 FN3FOP
      CHARACTER*4 FNDKBC
C
      PARAMETER (FN3VI  = 'CC3_VI', FN3VI2  = 'CC3_VI12', 
     *           FN3FOP = 'PTFOP' , FN3FOP2 = 'PTFOP2',
     *           FNCKJD = 'CKJDEL', FNTOC   = 'CCSDT_OC',
     *           FNDKBC = 'DKBC')
C
      LOGICAL locdbg
      PARAMETER (locdbg = .false.)
C
#if defined (SYS_CRAY) 
      REAL overlap(*),work(lwork)
      REAL tcon
      REAL d0,dsixth
      REAL dotcon(mxvec,ndtran)
#else 
      DOUBLE PRECISION overlap(*),work(lwork)
      DOUBLE PRECISION tcon
      DOUBLE PRECISION d0,dsixth
      DOUBLE PRECISION dotcon(mxvec,ndtran)
#endif 
      PARAMETER (d0 = 0.0d0)
      PARAMETER (dsixth = 1.0d0/6.0d0)
C
      integer ilstsym
C
      call qenter('cc3_lr_ovlp')
C
*---------------------------------------------------------------------*
*     Initializations and work space allocation:
*---------------------------------------------------------------------*
      if (imode.ne.1) then
        call quit('illegal imode in CC3_LR_OVLP.')
      end if

      if (listl(1:2).ne.'LE ' ) then
         WRITE(LUPRI,*) 'LISTL:',LISTL
         CALL QUIT('Unknown/Illegal list in CC3_LR_OVLP.')
      end if
!
      if ( listr(1:2).ne.'RE ') then
         WRITE(LUPRI,*) 'LISTR:',LISTR
         CALL QUIT('Unknown/Illegal list in CC3_LR_OVLP.')
      end if
C
C-------------------
C     Open the files
C-------------------
C

      LU3VI   = -1
      LU3VI2  = -1
      LU3FOP  = -1
      LU3FOP2 = -1
      LUCKJD  = -1
      LUTOC   = -1
      LUDKBC  = -1
C
      CALL WOPEN2(LU3VI,FN3VI,64,0)
      CALL WOPEN2(LU3VI2,FN3VI2,64,0)
      CALL WOPEN2(LU3FOP,FN3FOP,64,0)
      CALL WOPEN2(LU3FOP2,FN3FOP2,64,0)
      CALL WOPEN2(LUCKJD,FNCKJD,64,0)
      CALL WOPEN2(LUTOC,FNTOC,64,0)
      CALL WOPEN2(LUDKBC,FNDKBC,64,0)

*---------------------------------------------------------------------*
*     Loop over "left" vectors, pre-compute them and calculate
*     overlap matrix contributions in loop over "right" vectors
*---------------------------------------------------------------------*
      iovlp = 0
      do ileft = ilstnrst, ilstnrnd

        isyml   = ilstsym(listl,ileft)
        imull   = imulte(ileft)

        if (locdbg) then
          write(lupri,*) 'ileft,isyml,imull:',ileft,isyml,imull
          write(lupri,*) 'istart,iend:',istartend(ileft,1), 
     &                                  istartend(ileft,2)
        end if
C
        do iright = istartend(ileft,1), istartend(ileft,2)
          isymr   = ilstsym(listr,iright)
          imulr   = imulte(iright)

          if (locdbg) then
            write(lupri,*) 'iright,isymr,imulr:',iright,isymr,imulr
          end if

          if (isymr.eq.isyml .and. imulr.eq.imull) then
 
            tcon = 0.0d0
            call eb3_e3_ovlp(listl,ileft,isyml,listr,iright,isymr,
     *                       lu3vi,fn3vi,lu3vi2,fn3vi2,lu3fop,fn3fop,
     *                       lu3fop2,fn3fop2,luckjd,fnckjd,lutoc,fntoc,
     *                       ludkbc,fndkbc,tcon,work,lwork)
 

            iovlp = iovlp + 1
            overlap(iovlp) = overlap(iovlp) + dsixth*tcon

            if (locdbg) then
              write(lupri,*) 'iovlp:',iovlp
              write(lupri,*) 'tcon,overlap:',tcon,overlap(iovlp)
            end if

          end if
           
        end do
      end do
  
C
C-------------------
C     Close the files
C-------------------
C
      CALL WCLOSE2(LU3VI,FN3VI,'KEEP')
      CALL WCLOSE2(LU3VI2,FN3VI2,'KEEP')
      CALL WCLOSE2(LU3FOP,FN3FOP,'KEEP')
      CALL WCLOSE2(LU3FOP2,FN3FOP2,'KEEP')
      CALL WCLOSE2(LUCKJD,FNCKJD,'KEEP')
      CALL WCLOSE2(LUTOC,FNTOC,'KEEP')
      CALL WCLOSE2(LUDKBC,FNDKBC,'KEEP')
C
      call qexit('cc3_lr_ovlp')
C
      return
      end
*---------------------------------------------------------------------*
C  /* Deck eb3_e3_ovlp */
*=====================================================================*
      SUBROUTINE EB3_E3_OVLP(LISTEB,IDLSTEB,ISYMEB,LISTE,IDLSTE,ISYME, 
     *                       LU3VI,FN3VI,LU3VI2,FN3VI2,LU3FOP,FN3FOP,
     *                       LU3FOP2,FN3FOP2,LUCKJD,FNCKJD,LUTOC,FNTOC,
     *                       LUDKBC,FNDKBC,TCON,WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: compute overlap (L|R) = (EBAR3!E3) between triples part 
*             of left and right eigenvectors 
*
*    Poul Jorgensen and Filip Pawlowski, 18-Dec-2002, Aarhus
*
*    fixed for the memory leakage, 07-Oct-2003, Aarhus, FP.
*=====================================================================*
C
      IMPLICIT NONE 
C
#include <priunit.h>
#include <ccorb.h>
#include <ccexci.h>
#include <ccsdsym.h>
#include <dummy.h>
#include <iratdef.h>
#include <cc3t3d.h>
#include <ccinftap.h>
C
      INTEGER IDLSTEB,ISYMEB,IDLSTE,ISYME,ISYM0,ISINTETR,LWORK
      INTEGER IOPT,IOFF
      INTEGER ISYMC,ISCKB2,ISCKB2E,ISYCKB
      INTEGER ISYMB,ISYMBC,ISCKIJ,ISWEBMAT,ISWEMAT,ISCKD2
      INTEGER ISCKD2E,ISYALJBY,ISYALJCY
      INTEGER ISYML,ISYMCL,ISAIBJ,ISYMJ,ISYMBJ,ISYMAI,ISYAIL
      INTEGER KOFF1,NBJ,IADR
      INTEGER ISAIBJL,ISBJAIL,ISBJAI
      INTEGER KFOCKD,KT2TP,KEB1,KEB2TP,KE1,KE2TP,KLAMP0
      INTEGER KLAMH0,KFOCK0,KFCKCK,KEND1,LWRK1
      INTEGER KTROC0,KTROC0Y,KINTOC
      INTEGER KXIAJB,KTROCCG0,KTROCCL0
      INTEGER KTRVI0,KTRVI0Y,KVGDKCB0,KVGDKBC0,KVLDKCB0,KVLDKBC0
      INTEGER KINTVI
      INTEGER KINDSQWEB,KINDSQWE,KINDEXBY,KINDEXCY,KDIAGEB,KDIAGE
      INTEGER KWEBMAT,KWEMAT,KTMAT,KTRVI8,KTRVI8Y
      INTEGER KEBAIBJ,KEBJAI,KEAIBJ,KEND2,LWRK2
      INTEGER LUCKJDR,LU3SRTR,LUDELDR,LUDKBCR,LUWEB,LUWE
      INTEGER LU3VI,LU3VI2,LU3FOP,LU3FOP2,LUCKJD,LUTOC,LUDKBC
      INTEGER LENGTH,LENSQWEB,LENSQWE
COMMENT COMENT
      INTEGER kx3am
COMMENT COMENT

C
      CHARACTER*(*) LISTEB,LISTE
      CHARACTER CDUMMY*1 
      CHARACTER MODEL*10
C
      CHARACTER*(*) FN3VI,FN3VI2,FN3FOP,FN3FOP2,FNCKJD,FNTOC,FNDKBC
C
      CHARACTER*12 FN3SRTR, FNCKJDR, FNDELDR, FNDKBCR
      CHARACTER*11 FNWEB
      CHARACTER*10 FNWE
C
      PARAMETER (CDUMMY = ' ')
C
      PARAMETER(FN3SRTR  = 'CCSDT_FBMAT1', FNCKJDR  = 'CCSDT_FBMAT2',
     *          FNDELDR  = 'CCSDT_FBMAT3', FNDKBCR  = 'CCSDT_FBMAT4')
C
      PARAMETER(FNWEB = 'CC3_WEB_TMP', FNWE = 'CC3_WE_TMP')
C
#if defined (SYS_CRAY)
      REAL TCON,FREQEB,FREQE,XDIANORM,XOFFDIA,D3,D6,TWO
      REAL WORK(LWORK)
      REAL DDOT
      REAL THRSH
#else
      DOUBLE PRECISION TCON,FREQEB,FREQE,XDIANORM,XOFFDIA,D3,D6,TWO
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION DDOT
      DOUBLE PRECISION THRSH
#endif
C
      PARAMETER (D3 = 3.0d0)
      PARAMETER (D6 = 6.0d0)
      PARAMETER (TWO = 2.0d0)
      PARAMETER (THRSH = 1.0d-3)
C
      CALL QENTER('EB3_E3_OVLP')
C
C------------------------------------
C     Initial test of operator labels 
C------------------------------------
C
      if (LISTEB(1:2).ne.'LE ' ) then
         WRITE(LUPRI,*) 'LISTEB:',LISTEB
         CALL QUIT('Unknown/Illegal list in EB3_E3_OVLP.')
      end if
C
      if ( LISTE(1:2).ne.'RE ') then
         WRITE(LUPRI,*) 'LISTE:',LISTE
         CALL QUIT('Unknown/Illegal list in EB3_E3_OVLP.')
      end if
C
C--------------------------------
C     Open temporary files
C--------------------------------
C
      LU3SRTR  = -1
      LUCKJDR  = -1
      LUDELDR  = -1
      LUDKBCR  = -1
      LUWEB    = -1
      LUWE     = -1
C
      CALL WOPEN2(LU3SRTR,FN3SRTR,64,0)
      CALL WOPEN2(LUCKJDR,FNCKJDR,64,0)
      CALL WOPEN2(LUDELDR,FNDELDR,64,0)
      CALL WOPEN2(LUDKBCR,FNDKBCR,64,0)
      CALL WOPEN2(LUWEB,FNWEB,64,0)
      CALL WOPEN2(LUWE,FNWE,64,0)
C
C-------------------------
C     Some initializations
C-------------------------
C
      ISYM0 = 1 ! symmetry of H_0
      ISINTETR = MULD2H(ISYM0,ISYME)
C
      FREQEB = -EIGVAL(IDLSTEB)
      FREQE  = -EIGVAL(IDLSTE)
C
C----------------------------------
C     Check eigenvalues consistency
C----------------------------------
C
      IF (ABS(FREQEB - FREQE) .GT. THRSH) THEN
         WRITE(LUPRI,*)'Left eigenvalue = ', FREQEB
         WRITE(LUPRI,*)' not equal '
         WRITE(LUPRI,*)'right eigenvalue = ', FREQE
         CALL QUIT('Inconsistency in eigenvalues specification in'
     *              //' EB3_E3_OVLP')
      END IF
C
C--------------
C     Allcation
C--------------
C
      KFOCKD = 1
      KT2TP  = KFOCKD + NORBTS
      KEB1   = KT2TP  + NT2SQ(ISYM0)
      KEB2TP = KEB1   +  NT1AM(ISYMEB)
      KE1    = KEB2TP + NT2SQ(ISYMEB)
      KE2TP  = KE1    +  NT1AM(ISYME)
      KLAMP0 = KE2TP  + NT2SQ(ISYME)
      KLAMH0 = KLAMP0 + NLAMDT
      KFOCK0 = KLAMH0 + NLAMDT
      KFCKCK = KFOCK0 + N2BST(ISYM0)
      KEND1  = KFCKCK + NT1AM(ISYM0)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in EB3_E3_OVLP')
      END IF
C
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT If we want to sum the T3 amplitudes
COMMENT COMMENT COMMENT
      if (.false.) then
         kx3am  = kend1
         kend1 = kx3am + nrhft*nrhft*nrhft*nvirt*nvirt*nvirt
         call dzero(work(kx3am),nrhft*nrhft*nrhft*nvirt*nvirt*nvirt)
         lwrk1 = lwork - kend1
         if (lwrk1 .lt. 0) then
            write(lupri,*) 'Memory available : ',lwork
            write(lupri,*) 'Memory needed    : ',kend1
            call quit('Insufficient space (kx3am) in EB3_E3_OVLP')
         END IF
      endif
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT

C
C--------------------------------------------------------------
C     Read in t2 , square it and reorder
C--------------------------------------------------------------
C
      IOPT = 2
      CALL GET_T1_T2(IOPT,.FALSE.,DUMMY,WORK(KT2TP),'R0',0,ISYM0,
     *               WORK(KEND1),LWRK1)
C
C--------------------------------------------------------------
C     Getting the singles and doubles part of left eigenvectors
C--------------------------------------------------------------
C
      IOPT = 3 
      CALL GET_T1_T2(IOPT,.FALSE.,WORK(KEB1),WORK(KEB2TP),LISTEB,
     *               IDLSTEB,ISYMEB,WORK(KEND1),LWRK1) 
C
C--------------------------------------------------------------
C     Getting the singles and doubles part of right eigenvectors
C--------------------------------------------------------------
C
      IOPT = 3
      CALL GET_T1_T2(IOPT,.TRUE.,WORK(KE1),WORK(KE2TP),LISTE,IDLSTE,
     *               ISYME,WORK(KEND1),LWRK1)
C     IOPT = 3
C     CALL CC_RDRSP(LISTE,IDLSTE,ISYME,IOPT,MODEL,
C    *              WORK(KE1),WORK(KE2TP))
C     CALL CCLR_DIASCL(WORK(KE2TP),TWO,ISYME)
C     CALL CC_T2SQ(WORK(KE2TP),WORK(KEND1),ISYME)
C     CALL CC3_T2TP(WORK(KE2TP),WORK(KEND1),ISYME)

C     write(lupri,*)'KE2TP in real code'
C     call output(WORK(KE2TP),1,nrhft*nrhft*nvirt,1,nvirt,
C    *            nrhft*nrhft*nvirt,nvirt,1,lupri)
C
C--------------------------------------------------------------
C     Getting the Lambda^0 matrices
C--------------------------------------------------------------
C
      CALL GET_LAMBDA0(WORK(KLAMP0),WORK(KLAMH0),WORK(KEND1),LWRK1)
C
C-----------------------------------------------------------------
C     Read in Fock matrix in AO basis (from file) and transform to 
C     Lambda_0 basis.
C-----------------------------------------------------------------
C
      CALL GET_FOCK0(WORK(KFOCK0),WORK(KLAMP0),WORK(KLAMH0),WORK(KEND1),
     *               LWRK1)
C
C--------------------------------------------------------------
C     Sort the Fock matrix to get F(ck) block
C--------------------------------------------------------------
C
      CALL SORT_FOCKCK(WORK(KFCKCK),WORK(KFOCK0),ISYM0)
C
C--------------------------------------------------------------
C     Get orbital energies
C--------------------------------------------------------------
C
      CALL GET_ORBEN(WORK(KFOCKD),WORK(KEND1),LWRK1)
C
C------------------------------
C     Read in integrals for WE3
C------------------------------
C
C---------------------------------------------
C     Alocation for occupied integrals for WE3
C---------------------------------------------
C
      KTROC0 = KEND1
      KEND1  = KTROC0 + NTRAOC(ISYM0)
      LWRK1  = LWORK  - KEND1
C
      KTROC0Y = KEND1
      KEND1   = KTROC0Y + NTRAOC(ISINTETR)
C
      KINTOC = KEND1
      KEND1  = KINTOC + MAX(NTOTOC(ISYM0),NTOTOC(ISINTETR))
      LWRK1  = LWORK  - KEND1
C
C--------------------------------------------------------------
C     Construction of one-index transformed occupied and virtual
C     integrals needed for WE3
C--------------------------------------------------------------
C
C     Construct one-index transformed occupied integrals (FNCKJDR)
C
      CALL CC3_BARINT(WORK(KE1),ISYME,WORK(KLAMP0),
     *                WORK(KLAMH0),WORK(KEND1),LWRK1,
     *                LU3SRTR,FN3SRTR,LUCKJDR,FNCKJDR)
C
      CALL CC3_SORT1(WORK(KEND1),LWRK1,2,ISYME,LU3SRTR,FN3SRTR,
     *               LUDELDR,FNDELDR,IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *               IDUMMY,CDUMMY)
C
C     Virt integrals for WE3 dumped to the file
C
      CALL CC3_SINT(WORK(KLAMH0),WORK(KEND1),LWRK1,ISYME,
     *              LUDELDR,FNDELDR,LUDKBCR,FNDKBCR)
C
C---------------------------------
C     Construct occupied integrals
C---------------------------------
C
      IOFF = 1
      IF (NTOTOC(ISYM0) .GT. 0) THEN
         CALL GETWA2(LUCKJD,FNCKJD,WORK(KINTOC),IOFF,NTOTOC(ISYM0))
      ENDIF
C
C----------------------------------------------------------------------
C     Transform (ai|j delta) integrals to (ai|j k) and sort as (ij,k,a)
C----------------------------------------------------------------------
C
      CALL CC3_TROCC(WORK(KINTOC),WORK(KTROC0),WORK(KLAMP0),
     *               WORK(KEND1),LWRK1,ISYM0)
C     write(lupri,*)'(ij,k,a) occupied int in real'
C     call output(WORK(KTROC0),1,nrhft*nrhft*nrhft,1,nvirt,
C    *            nrhft*nrhft*nrhft,nvirt,1,lupri)
C
C------------------------------------------
C     B transformed Occupied integrals.
C-----------------------------------------
C
      IOFF = 1
      IF (NTOTOC(ISYME) .GT. 0) THEN
         CALL GETWA2(LUCKJDR,FNCKJDR,WORK(KINTOC),IOFF,
     *               NTOTOC(ISYME))
      ENDIF
C
      CALL CC3_TROCC(WORK(KINTOC),WORK(KTROC0Y),WORK(KLAMP0),
     *               WORK(KEND1),LWRK1,ISYME)
C
C----------------
C     Allocation
C----------------
C
      KXIAJB  = KEND1
      KEND1   = KXIAJB  + NT2AM(ISYM0)
C
      KTROCCG0   = KEND1
      KTROCCL0   = KTROCCG0   + NTRAOC(ISYM0)
      KEND1      = KTROCCL0   + NTRAOC(ISYM0)
      LWRK1      = LWORK      - KEND1
C
C------------------------
C     Construct L(ia,jb).
C------------------------
C
      LENGTH = IRAT*NT2AM(ISYM0)
C
      REWIND(LUIAJB)
      CALL READI(LUIAJB,LENGTH,WORK(KXIAJB))
C
      CALL CCSD_TCMEPK(WORK(KXIAJB),1.0D0,ISYM0,1)
C
C------------------------
C     Construct occupied integrals for WEBMAT
C-----------------------
C
      IOFF = 1
      IF (NTOTOC(ISYM0) .GT. 0) THEN
         CALL GETWA2(LUTOC,FNTOC,WORK(KINTOC),IOFF,NTOTOC(ISYM0))
      ENDIF
C
C----------------------------------------------------------------------
C     Transform (ia|j delta) integrals to (ia|j k-) and sort  G(j,i,k-,a)
C----------------------------------------------------------------------
C
      CALL CCLR_TROCC(WORK(KINTOC),WORK(KTROCCG0),
     *                 WORK(KLAMH0),ISYM0,
     *                 WORK(KEND1),LWRK1)
C
C----------------------------------------------------------------------
C integrals g(ia|j k-) sorted as G(j,i,k-,a) resorted as G(k-,i,j,a)
C                       on input KTROCCG   on output     KTROCCG
C           L(ia|j k-)                                   KTROCCL
C----------------------------------------------------------------------
C
      CALL CC3_SRTOCC(WORK(KTROCCG0),WORK(KTROCCL0),ISYM0)
C
C----------------------------
C     Loop over C
C----------------------------
C
      DO ISYMC = 1,NSYM
C
C           -------------------------
C           Memory allocation for WE3
C           -------------------------
C
         ISCKB2  = MULD2H(ISYM0,ISYMC)
         ISCKB2E = MULD2H(ISINTETR,ISYMC)
C
         KTRVI0  = KEND1
         KEND1   = KTRVI0  + NCKATR(ISCKB2)
         LWRK1   = LWORK  - KEND1
C
         KTRVI0Y  = KEND1
         KEND1    = KTRVI0Y  + NCKATR(ISCKB2E)
         LWRK1    = LWORK    - KEND1
C
C           ----------------------------
C           Memory allocation for WEBAR3
C           ----------------------------
C
         ISYCKB  = MULD2H(ISYM0,ISYMC)
C
         KVGDKCB0  = KEND1
         KVGDKBC0  = KVGDKCB0  + NCKATR(ISYCKB)
         KVLDKCB0  = KVGDKBC0  + NCKATR(ISYCKB)
         KVLDKBC0  = KVLDKCB0  + NCKATR(ISYCKB)
         KEND1     = KVLDKBC0  + NCKATR(ISYCKB)
         LWRK1     = LWORK     - KEND1
C
         KINTVI = KEND1
         KEND1  = KINTVI + NCKA(ISYCKB)
         LWRK1  = LWORK  - KEND1
C
         DO C = 1,NVIR(ISYMC)
C
C-----------------------------------------------
C           Read virtual integrals used WE3
C-----------------------------------------------
C
            IOFF = ICKBD(ISCKB2,ISYMC) + NCKATR(ISCKB2)*(C - 1) + 1
            IF (NCKATR(ISCKB2) .GT. 0) THEN
               CALL GETWA2(LUDKBC,FNDKBC,WORK(KTRVI0),IOFF,
     &                        NCKATR(ISCKB2))
            ENDIF
C
C-----------------------------------------------------------------------
C        Read B transformed virtual integrals used for WE3 
C-----------------------------------------------------------------------
C
            IOFF = ICKBD(ISCKB2E,ISYMC) + NCKATR(ISCKB2E)*(C - 1) + 1
            IF (NCKATR(ISCKB2) .GT. 0) THEN
               CALL GETWA2(LUDKBCR,FNDKBCR,WORK(KTRVI0Y),IOFF,
     &                     NCKATR(ISCKB2E))
            ENDIF
C
C---------------------------------------
C           VIRTUAL INTEGRALS FOR WEBAR3
C---------------------------------------
C
C           Read in g(c1^h k1^p | delta b2^h) integrals
C
            IOFF = ICKAD(ISYCKB,ISYMC) + NCKA(ISYCKB)*(C - 1) + 1
            IF (NCKA(ISYCKB) .GT. 0) THEN
C
               CALL GETWA2(LU3VI,FN3VI,WORK(KINTVI),IOFF,
     &                        NCKA(ISYCKB))
            ENDIF
C
C          Transform g(c1^h k1^p | delta b2^h) integrals
C          to g(c1^h k1^p | d2^p b2^h) with lamp0
C
            CALL CCLR_TRVIR(WORK(KINTVI),WORK(KVGDKCB0),
     *                      WORK(KLAMP0),ISYM0,
     *                      ISYMC,C,ISYM0,WORK(KEND1),LWRK1)
C
C          Sort g(c1^h k1^p | d2^p b2^h) as 
C          g(d2^p k1^p | c1^h b2^h)
C
            CALL CCSDT_SRVIR3(WORK(KVGDKCB0),WORK(KEND1),ISYMC,C,ISYM0)
C
C           Read in g(b2^h k1^p | delta c1^h) integrals and 
C           transform and sort --- see above
C
            IOFF = ICKAD(ISYCKB,ISYMC) + NCKA(ISYCKB)*(C - 1) + 1
            IF (NCKA(ISYCKB) .GT. 0) THEN
C
               CALL GETWA2(LU3VI2,FN3VI2,WORK(KINTVI),IOFF,
     &                        NCKA(ISYCKB))
            ENDIF
C
            CALL CCLR_TRVIR(WORK(KINTVI),WORK(KVGDKBC0),
     *                          WORK(KLAMP0),ISYM0,
     *                          ISYMC,C,ISYM0,WORK(KEND1),LWRK1)
C
            CALL CCSDT_SRVIR3(WORK(KVGDKBC0),WORK(KEND1),ISYMC,C,ISYM0)
C
C           Read in L(c1^h k1^p | delta b2^h) integrals
C
            IOFF = ICKAD(ISYCKB,ISYMC) + NCKA(ISYCKB)*(C - 1) + 1
            IF (NCKA(ISYCKB) .GT. 0) THEN
C
               CALL GETWA2(LU3FOP,FN3FOP,WORK(KINTVI),IOFF,
     &                        NCKA(ISYCKB))
            ENDIF
C
C          Transform L(c1^h k1^p | delta b2^h) integrals
C          to L(c1^h k1^p | d2^p b2^h) 
C
            CALL CCLR_TRVIR(WORK(KINTVI),WORK(KVLDKCB0),
     *                          WORK(KLAMP0),ISYM0,
     *                          ISYMC,C,ISYM0,WORK(KEND1),LWRK1)
C
C          Sort L(c1^h k1^p | d2^p b2^h) as 
C          L(d2^p k1^p | c1^h b2^h)
C
            CALL CCSDT_SRVIR3(WORK(KVLDKCB0),WORK(KEND1),ISYMC,C,ISYM0)

C
C           Read in L(b2^h k1^p | delta c1^h) integrals and 
C           transform and sort --- see above
C
            IOFF = ICKAD(ISYCKB,ISYMC) + NCKA(ISYCKB)*(C - 1) + 1
            IF (NCKA(ISYCKB) .GT. 0) THEN
C
               CALL GETWA2(LU3FOP2,FN3FOP2,WORK(KINTVI),IOFF,
     &                           NCKA(ISYCKB))
C
            END IF
C
            CALL CCLR_TRVIR(WORK(KINTVI),WORK(KVLDKBC0),
     *                          WORK(KLAMP0),ISYM0,
     *                          ISYMC,C,ISYM0,WORK(KEND1),LWRK1)
C
            CALL CCSDT_SRVIR3(WORK(KVLDKBC0),WORK(KEND1),ISYMC,C,ISYM0)
C
C----------------------------
C           Loop over B
C----------------------------
C
            DO ISYMB = 1,NSYM
C
               ISYMBC  = MULD2H(ISYMB,ISYMC)
               ISCKIJ  = MULD2H(ISYMBC,ISYM0)
               ISWEBMAT = MULD2H(ISCKIJ,ISYMEB)
               ISWEMAT = MULD2H(ISCKIJ,ISYME)
               ISCKD2  = MULD2H(ISYM0,ISYMB)
               ISCKD2E = MULD2H(ISINTETR,ISYMB)
               ISYALJBY  = MULD2H(ISYMB,ISYMEB)
               ISYALJCY = MULD2H(ISYMC,ISYMEB)
C
               KINDSQWEB = KEND1
               KINDSQWE  = KINDSQWEB  + (6*NCKIJ(ISWEBMAT) - 1)/IRAT + 1
               KINDEXBY  = KINDSQWE  + (6*NCKIJ(ISWEMAT) - 1)/IRAT + 1
               KINDEXCY  = KINDEXBY  + (NCKI(ISYALJBY)  - 1)/IRAT + 1
               KEND2     = KINDEXCY  + (NCKI(ISYALJCY) - 1)/IRAT + 1
C
               KDIAGEB = KEND2 
               KDIAGE  = KDIAGEB + NCKIJ(ISWEBMAT)
               KWEBMAT = KDIAGE  + NCKIJ(ISWEMAT)
               KWEMAT  = KWEBMAT + NCKIJ(ISWEBMAT)
               KTMAT   = KWEMAT  + NCKIJ(ISWEMAT)
               KEND2   = KTMAT   + MAX(NCKIJ(ISWEBMAT),NCKIJ(ISWEMAT))
C
               KTRVI8  = KEND2
               KEND2   = KTRVI8  + NCKATR(ISCKD2)
C
               KTRVI8Y = KEND2
               KEND2   = KTRVI8Y + NCKATR(ISCKD2E)
               LWRK2   = LWORK - KEND2
C
C---------------------------------------------
C           Construct part of the diagonal.
C---------------------------------------------
C
               CALL CC3_DIAG(WORK(KDIAGEB),WORK(KFOCKD),ISWEBMAT)
               CALL CC3_DIAG(WORK(KDIAGE),WORK(KFOCKD),ISWEMAT)
C
C              -----------------------
C              Construct index arrays.
C              -----------------------
C
               CALL CC3_INDEX(WORK(KINDEXBY),ISYALJBY)
               CALL CC3_INDEX(WORK(KINDEXCY),ISYALJCY)
               LENSQWEB  = NCKIJ(ISWEBMAT)
               CALL CC3_INDSQ(WORK(KINDSQWEB),LENSQWEB,ISWEBMAT)
               LENSQWE  = NCKIJ(ISWEMAT)
               CALL CC3_INDSQ(WORK(KINDSQWE),LENSQWE,ISWEMAT)


               DO B = 1,NVIR(ISYMB)
C
C--------------------------------------------
C                 Initialize WEBMAT and WEMAT
C--------------------------------------------
C
                  CALL DZERO(WORK(KWEMAT),NCKIJ(ISWEMAT))
                  CALL DZERO(WORK(KWEBMAT),NCKIJ(ISWEBMAT))
C
C-----------------------------------------------
C           Read virtual integrals used in WE3
C-----------------------------------------------
C
                  IOFF = ICKBD(ISCKD2,ISYMB) +
     &                   NCKATR(ISCKD2)*(B - 1) + 1
                  IF (NCKATR(ISCKD2) .GT. 0) THEN
                     CALL GETWA2(LUDKBC,FNDKBC,WORK(KTRVI8),IOFF,
     &                           NCKATR(ISCKD2))
                  ENDIF
C--------------------------------------------------------------------
C           Read B transformed virtual integrals used in WE3
C--------------------------------------------------------------------

                  IOFF = ICKBD(ISCKD2E,ISYMB) +
     &                   NCKATR(ISCKD2E)*(B - 1) + 1
                  IF (NCKATR(ISCKD2E) .GT. 0) THEN
                     CALL GETWA2(LUDKBCR,FNDKBCR,WORK(KTRVI8Y),IOFF,
     &                           NCKATR(ISCKD2E))
                  ENDIF
C
C--------------------------------------------------------------
C    calculate     <L2Y|[H^,tau3]|HF>  - construction of WEBAR3
C--------------------------------------------------------------
C              
                  CALL WBARBD_TMAT(WORK(KEB2TP),ISYMEB,WORK(KWEBMAT),
     *                             WORK(KTMAT),ISWEBMAT,WORK(KFCKCK),
     *                             ISYM0,WORK(KVLDKBC0),
     *                             WORK(KVLDKCB0),
     *                             WORK(KVGDKBC0),WORK(KVGDKCB0),
     *                             WORK(KTROCCL0),WORK(KTROCCG0),
     *                             ISYM0,WORK(KEND2),LWRK2,
     *                             WORK(KINDEXBY),WORK(KINDEXCY),
     *                             WORK(KINDSQWEB),LENSQWEB,ISYMB,B,
     *                             ISYMC,C)
C
C--------------------------------------------------------------
C    calculate     <L1Y|[H^,tau3]|HF>  - construction of WEBAR3
C--------------------------------------------------------------
C
                  CALL WBARBD_L1(WORK(KEB1),ISYMEB,WORK(KTMAT),
     *                           WORK(KXIAJB),ISYM0,
     *                           WORK(KWEBMAT),WORK(KEND2),LWRK2,
     *                           WORK(KINDSQWEB),LENSQWEB,ISYMB,B,
     *                           ISYMC,C)
C
C------------------------------------------------
C     Divide by the energy difference and
C     remove the forbidden elements
C------------------------------------------------
C
C
                  CALL WBD_DIA(B,ISYMB,C,ISYMC,-FREQEB,
     *                         ISWEBMAT,WORK(KWEBMAT),
     *                         WORK(KDIAGEB),WORK(KFOCKD))
C
                  CALL T3_FORBIDDEN(WORK(KWEBMAT),ISYMEB,
     *                              ISYMB,B,ISYMC,C)
COMMENT COMMENT COMMENT
C      call sum_pt3(work(kwebmat),isymb,b,isymc,c,
C    *             ISWEBMAT,work(kx3am),4)
COMMENT COMMENT COMMENT
C
C--------------------------------------------------------------
C                   -<mu3 | [[H,T1^Y],T2] | HF> -<mu3| [H,T2^Y] |HF> 
C
C                        - construction of WE3
C--------------------------------------------------------------
C
C                 <mu3 | [[H,E1],T2] | HF>
C
                  CALL WBD_GROUND(WORK(KT2TP),ISYM0,WORK(KTMAT),
     *                            WORK(KTRVI0Y),WORK(KTRVI8Y),
     *                            WORK(KTROC0Y),ISINTETR,WORK(KWEMAT),
     *                            WORK(KEND2),LWRK2,
     *                            WORK(KINDSQWE),LENSQWE,
     *                            ISYMB,B,ISYMC,C)
C
C                  <mu3| [H,E2] |HF> 
C
                  CALL WBD_GROUND(WORK(KE2TP),ISYME,WORK(KTMAT),
     *                            WORK(KTRVI0),WORK(KTRVI8),
     *                            WORK(KTROC0),ISYM0,WORK(KWEMAT),
     *                            WORK(KEND2),LWRK2,
     *                            WORK(KINDSQWE),LENSQWE,
     *                            ISYMB,B,ISYMC,C)
C
C------------------------------------------------
C     Divide by the energy difference and
C     remove the forbidden elements
C------------------------------------------------
C
C
                  CALL WBD_DIA(B,ISYMB,C,ISYMC,-FREQE,
     *                         ISWEMAT,WORK(KWEMAT),
     *                         WORK(KDIAGE),WORK(KFOCKD))
C
                  CALL T3_FORBIDDEN(WORK(KWEMAT),ISYME,
     *                              ISYMB,B,ISYMC,C)
COMMENT COMMENT COMMENT
C      call sum_pt3(work(kwemat),isymb,b,isymc,c,
C    *             ISWEMAT,work(kx3am),4)
COMMENT COMMENT COMMENT

C
C----------------------------------------------------------------
C                 Calculate sum_(aikj) WEB^BC(aikj) * WE^BC(aikj)
C----------------------------------------------------------------
C
                  IF ( (NCKIJ(ISWEMAT) .GT. 0) 
     *                 .AND. (ISWEMAT .EQ. ISWEBMAT)) THEN 
C
                     XDIANORM = DDOT(NCKIJ(ISWEMAT),WORK(KWEMAT),1,
     *                               WORK(KWEBMAT),1)
C                    
                     TCON = TCON + D3*XDIANORM
C
                  END IF
C
C-----------------------------------------------------------------------
C                 Construct TMAT^AC(bjik) = WE^AC(bkji) + WE^AC(bjki)
C-----------------------------------------------------------------------
C
C                 CALL WE_TILDE(WORK(KWEBMAT),ISWEBMAT,WORK(KTMAT),
C    *                          WORK(KINDSQWEB),LENSQWEB)
C
C                 -------------------------------------
C                 Write WEBMAT as EB^C(ai,bj,l) to disc
C                 -------------------------------------
                  CALL WRITE_T3_DL(LUWEB,FNWEB,WORK(KWEBMAT),ISYMEB,
     *                            ISYMC,ISYMB,B)
C
C-----------------------------------------------------------------------
C                 Construct TMAT^AC(bjik) = WE^AC(bkji) + WE^AC(bjki)
C-----------------------------------------------------------------------
C
C                 CALL WE_TILDE(WORK(KWEMAT),ISWEMAT,WORK(KTMAT),
C    *                          WORK(KINDSQWE),LENSQWE)
C
C                 -------------------------------------
C                 Write WTMAT as E^C(bj,ai,l) to disc
C                 -------------------------------------
                  CALL WRITE_T3_DL(LUWE,FNWE,WORK(KWEMAT),ISYME,ISYMC,
     *                             ISYMB,B)
C
               END DO ! B
            END DO    ! ISYMB
C
C-----------------------------------------------------
C          Calculate WEB^BC(aijk) TMAT^AC(bjik)
C-----------------------------------------------------
C
            ISAIBJL = MULD2H(ISYMEB,ISYMC)
            ISBJAIL = MULD2H(ISYME,ISYMC)
            DO ISYML = 1, NSYM
C
               ISAIBJ = MULD2H(ISAIBJL,ISYML)
               ISBJAI = MULD2H(ISBJAIL,ISYML)
C
C------------------------
C              Allocation 
C------------------------
C
cfp memory problem fix.
               !reuse KEND1 since stuff from B-loop is not needed any more
               KEBAIBJ = KEND1
               KEBJAI  = KEBAIBJ + NT2SQ(ISAIBJ)
               KEAIBJ  = KEBJAI  + NT2SQ(ISBJAI)
               KEND2   = KEAIBJ  + NT2SQ(ISBJAI)
               LWRK2   = LWORK   - KEND2
C
               IF (LWRK2 .LT. 0) THEN 
                  WRITE(LUPRI,*)'Memory available: ', LWORK
                  WRITE(LUPRI,*)'Memory needed   : ', KEND2
                  CALL QUIT('Insufficient memory in EB3_E3_OVLP (2)')
               END IF
C
               DO L = 1, NRHF(ISYML)
C
C                 --------------------------------------------
C                 Read WEB^C(ai,bj,l) from file:
C                 --------------------------------------------
C
                  IADR = ISWTL(ISAIBJ,ISYML) + NT2SQ(ISAIBJ)*(L-1) + 1
                  CALL GETWA2(LUWEB,FNWEB,WORK(KEBAIBJ),IADR,
     *                        NT2SQ(ISAIBJ))
C
C                 write(lupri,*)'WEBAIBJ after reading,C,L',C,L
C                 call output(WORK(KEBAIBJ),1,NT1AM(1),1,NT1AM(1),
C    *                        NT1AM(1),NT1AM(1),1,lupri)
C
C                 --------------------------------------------
C                 Read WE^C(bj,ai,l) from file:
C                 --------------------------------------------
C
                  IADR = ISWTL(ISBJAI,ISYML) + NT2SQ(ISBJAI)*(L-1) + 1
                  CALL GETWA2(LUWE,FNWE,WORK(KEBJAI),IADR,
     *                        NT2SQ(ISBJAI))
C
C                 write(lupri,*)'WEBJAI after reading,C,L',C,L
C                 call output(WORK(KEBJAI),1,NT1AM(1),1,NT1AM(1),
C    *                        NT1AM(1),NT1AM(1),1,lupri)
C
C---------------------------------------------------------
C                 Transpose WE^C(bj,ai,l) to WE^C(ai,bj,l)
C---------------------------------------------------------
C
                  CALL TRANS_AIBJ_BJAI(WORK(KEBJAI),WORK(KEAIBJ),ISBJAI)
C
C                 write(lupri,*)'WEAIBJ after transformation ,C,L',C,L
C                 call output(WORK(KEAIBJ),1,NT1AM(1),1,NT1AM(1),
C    *                        NT1AM(1),NT1AM(1),1,lupri)
C
C-----------------------------------------------------------
C                 Multiply WEB^C(ai,bj,l) with WE^C(ai,bj,l)
C-----------------------------------------------------------
C
                  IF ( (NT2SQ(ISAIBJ) .GT. 0) 
     *                  .AND. (ISAIBJ .EQ. ISBJAI) ) THEN
C
                     XOFFDIA = DDOT(NT2SQ(ISAIBJ),WORK(KEBAIBJ),1,
     *                            WORK(KEAIBJ),1)
C
                     TCON = TCON + D6*XOFFDIA
C
                  END IF
C
               END DO ! L
            END DO    ! ISYML
C
         END DO       ! C
      END DO          ! ISYMC
C
COMMENT COMMENT
C      write(lupri,*) 'WEBMAT in EB3_E3_OVLP : '
C      call print_pt3(work(kx3am),ISYMEB,4)
C      write(lupri,*) 'WEMAT in EB3_E3_OVLP : '
C      call print_pt3(work(kx3am),ISYME,4)
COMMENT COMMENT
C
C
C--------------------------------
C     Close the temporary files 
C--------------------------------
C
      CALL WCLOSE2(LU3SRTR,FN3SRTR,'DELETE')
      CALL WCLOSE2(LUCKJDR,FNCKJDR,'DELETE')
      CALL WCLOSE2(LUDELDR,FNDELDR,'DELETE')
      CALL WCLOSE2(LUDKBCR,FNDKBCR,'DELETE')
      CALL WCLOSE2(LUWEB,FNWEB,'DELETE')
      CALL WCLOSE2(LUWE,FNWE,'DELETE')
C
C--------
C     End.
C--------
C
      CALL QEXIT('EB3_E3_OVLP')
C
      RETURN
      END
C
C
C  /* Deck we_tilde */
      SUBROUTINE WE_TILDE(WEMAT,ISWEMAT,TMAT,INDSQ,LENSQ)
C
C-----------------------------------------------------------------------
C     Construct TMAT^AC(bjik) = WE^AC(bkji) + WE^AC(bjki)
C-----------------------------------------------------------------------
C
      IMPLICIT NONE
#include <priunit.h>
C
      INTEGER ISWEMAT,LENSQ,INDSQ(LENSQ,6)
      INTEGER I
C
#if defined (SYS_CRAY)
      REAL WEMAT(*),TMAT(*)
#else
      DOUBLE PRECISION WEMAT(*),TMAT(*)
#endif
C
      CALL QENTER('WE_TILDE')

      DO I = 1, LENSQ
          TMAT(I) = WEMAT(INDSQ(I,3))
      ENDDO
C
      CALL QEXIT('WE_TILDE')
C
      RETURN
      END
C  /* Deck trans_aibj_bjai */
      SUBROUTINE TRANS_AIBJ_BJAI(AIBJ,BJAI,ISAIBJ)
C******************************************************************
C     Transpose:
C
C     t2am(bj,ai) = t2tp(ai,bj) 
C
C     Poul Jorgensen and Filip Pawlowski, 18-Dec-2002, Aarhus
C
C******************************************************************
C
      IMPLICIT NONE 
C
#include <ccorb.h>
#include <ccsdsym.h>
C
      INTEGER ISAIBJ,ISYMB,ISYAIJ,ISYMJ,ISYMBJ,ISYMAI,ISYMI,ISYMA
      INTEGER NBJ,NAI,KOFF1,KOFF2
C
#if defined (SYS_CRAY)
      REAL AIBJ(*),BJAI(*)
#else
      DOUBLE PRECISION AIBJ(*),BJAI(*)
#endif
C
      CALL QENTER('TRANS_AIBJ_BJAI')
C
      DO ISYMB = 1,NSYM
C
         ISYAIJ = MULD2H(ISYMB,ISAIBJ)
C
         DO ISYMJ = 1,NSYM
C
            ISYMBJ = MULD2H(ISYMB,ISYMJ)
            ISYMAI = MULD2H(ISYMBJ,ISAIBJ)
C
            DO ISYMI = 1,NSYM
               ISYMA = MULD2H(ISYMAI,ISYMI)
C
               DO J = 1,NRHF(ISYMJ)
C
                  DO B = 1,NVIR(ISYMB)
C
                     NBJ   = IT1AM(ISYMB,ISYMJ)
     *                     + NVIR(ISYMB)*(J - 1) + B
C
                     DO I = 1,NRHF(ISYMI)
C
                        DO A = 1,NVIR(ISYMA)
C
                           NAI   = IT1AM(ISYMA,ISYMI)
     *                           + NVIR(ISYMA)*(I - 1) + A
C
                           KOFF1 = IT2SQ(ISYMAI,ISYMBJ)
     *                           + NT1AM(ISYMAI)*(NBJ - 1) 
     *                           + IT1AM(ISYMA,ISYMI) 
     *                           + NVIR(ISYMA)*(I-1)
     *                           + A
C
                           KOFF2 = IT2SQ(ISYMBJ,ISYMAI)
     *                           + NT1AM(ISYMBJ)*(NAI - 1) 
     *                           + IT1AM(ISYMB,ISYMJ) 
     *                           + NVIR(ISYMB)*(J-1)
     *                           + B
C
                           BJAI(KOFF2) = AIBJ(KOFF1)
C
                        END DO ! A
                     END DO    ! I
                  END DO       ! B
               END DO          ! J
            END DO             ! ISYMI
         END DO                ! ISYMJ
      END DO                   ! ISYMB
C
      CALL QEXIT('TRANS_AIBJ_BJAI')
C
      RETURN
      END
! -- end of ccexnorm.F --
