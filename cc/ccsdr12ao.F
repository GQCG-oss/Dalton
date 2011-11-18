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
      subroutine ccsdr12ao(ccsdr12,
     &                     t2am,xlambdap,xlambdah,
     &                     fniadj,luiadj,fnijda,luijda,
     &                     cpfil,lucp,dpfil,ludp,e1pim,
     &                     timintr12,timrdao,timtrbt,
     &                     timc,timd,timt2tr,timt2bt,
     &                     work,lwork)
*----------------------------------------------------------------------*
*  Purpose: compute contributions to CCSDR12 requiring the calculation
*           of integrals with auxiliary basis function indices
*
* C. Haettig, C. Neiss, spring 2006
*----------------------------------------------------------------------*
      implicit none
#include <priunit.h>
#include <dummy.h>     
#include <maxorb.h>
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
#include <ccsdio.h>
#include <second.h>
      logical locdbg
      parameter (locdbg = .false.)

      integer isym0
      parameter (isym0 = 1)

      real*8  zero,half,one,two,xmhalf,xmone
      parameter (zero = 0.0d0, half = 0.5d0, one = 1.0d0, two = 2.0d0)
      parameter (xmhalf = -0.5d0, xmone= -1.0d0 )

* input:
      logical ccsdr12
      integer lwork, lucp, ludp, luiadj, luijda
      character*(*) cpfil, dpfil, fniadj, fnijda
      real*8  work(*), xlambdap(*), xlambdah(*), t2am(*), e1pim(*),
     &     timintr12, timrdao, timtrbt, timc, timd, timt2tr, timt2bt

* local
      logical lauxg, temp_direct
      integer indexa(mxcorb)
      integer ioff, ibasx(8), kend1, lwrk1, kendsv, lwrksv,
     &        icdel1, isymd1, ntot, illl, numdis, idel2, idel, isymd,
     &        isydis, irecord, leniaj, isygam, isyalbe, igam, iadr,
     &        icon, iv, isym, nviraop(8), iviraop(8,8), isym1, isym2,
     &        kxint, kxiadj, kxijda, kend3, lwrk3, koffg, kdsrhf,
     &        kt2amt, kend2, lwrk2, kfckvaop, isyfck, lunit, ioptr12
      integer kccfb1,kindxb,kfree,lfree,ntosym,icount,ibasd, iopte,
     &        kodcl1,kodcl2,kodbc1,kodbc2,krdbc1,krdbc2,
     &        kodpp1,kodpp2,krdpp1,krdpp2,krecnr
      real*8  dtime, factc, factd, ddot

      call qenter('ccsdr12ao')
      if (.not.dumpcd) call quit('CCSDR12AO requires DUMPCD=.TRUE.')

      ! switch to integrals with delta index from auxiliary basis:
      mbsmax = 5
      loopdp = .true.

      ! integrals with auxiliary basis function are not on file 
      ! -> switch locally to direct mode to calculate them 
      temp_direct = direct
      direct = .true.

      ! substract aux. functions when calculating index of g index
      ! within irrep from index running over both basis sets and
      ! all symmetries
      ! (needed since isao is overwritten)
      lauxg = .true.

      ioff    = 0
      ibasx(1) = 0
      do isym = 1,nsym
        if (isym.gt.1) ibasx(isym) = ibasx(isym-1)+mbas2(isym-1)
        if(ioff+mbas1(isym)+mbas2(isym).gt.mxcorb)call quit('CCSDR12AO')
        do i = 1,mbas1(isym)+mbas2(isym)
          ioff = ioff + 1
          isao(ioff) = isym
        end do
      end do

      kend1 = 1
      lwrk1 = lwork

*----------------------------------------------------------------------*
*     initialization of symmetry arrays and allocation of work space
*     for Fhat(del',a)
*----------------------------------------------------------------------*
      do isym = 1, nsym
        icount = 0
        do isym2 = 1, nsym
          isym1 = muld2h(isym2,isym)
          iviraop(isym1,isym2) = icount
          icount = icount + nvir(isym1)*mbas2(isym2)
        end do
        nviraop(isym) = icount
      end do

      isyfck = 1 ! symmetry of the Fhat matrix that will be computed
  
      kfckvaop = kend1
      kend1   = kfckvaop + nviraop(isyfck)
      lwrk1   = lwork   - kend1
      if (lwrk1 .lt. 0) call quit('Insufficient core in CCSDR12AO')

      call dzero(work(kfckvaop),nviraop(isyfck))

*----------------------------------------------------------------------*
*     prepare cluster amplutides with transposed occupied indices:
*----------------------------------------------------------------------*
      if ((.not. direct) .and. t2tcor) then
         kt2amt = kend1
         kend1  = kt2amt + nt2sq(1)
         lwrk1  = lwork  - kend1
         if (lwrk1 .lt. 0) call quit('Insufficient core in CCSDR12AO')
         call dcopy(nt2sq(1),t2am,1,work(kt2amt),1)
         call ccsd_t2tp(work(kt2amt),work(kend1),lwrk1,isym0)
      end if

*----------------------------------------------------------------------*
*     intialize integral program
*----------------------------------------------------------------------*
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
     &        call quit('Insufficient work space in CC_MOFCONR12')
            call eridi1(kodcl1,kodcl2,kodbc1,kodbc2,krdbc1,krdbc2,
     &                kodpp1,kodpp2,krdpp1,krdpp2,
     &                kfree,lfree,kend1,work(kccfb1),work(kindxb),
     &                work(kend1),lwrk1,ipreri)
            kend1 = kfree
            lwrk1 = lfree
         endif
         timintr12 = timintr12 + ( second() - dtime )
         ntosym = 1
      else
         ntosym = nsym
      endif
  
*----------------------------------------------------------------------*
*     start the loop over integral distributions:
*----------------------------------------------------------------------*

      kendsv = kend1
      lwrksv = lwrk1

      icdel1 = 0
      do isym = 1, nsym
        icdel1 = icdel1  + nt2bcd(isym)*mbas1(isym)
      end do  
  
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
  
          if (direct) then
             dtime = second()
             kend1 = kendsv
             lwrk1 = lwrksv
             if (herdir) then
                call herdi2(work(kend1),lwrk1,indexa,illl,numdis,
     &                      ipreri)
             else
                call eridi2(illl,indexa,numdis,0,0,
     &                    work(kodcl1),work(kodcl2),work(kodbc1),
     &                    work(kodbc2),work(krdbc1),work(krdbc2),
     &                    work(kodpp1),work(kodpp2),work(krdpp1),
     &                    work(krdpp2),work(kccfb1),work(kindxb),
     &                    work(kend1), lwrk1,ipreri)
             endif
  
             krecnr = kend1
             kend1  = krecnr + (nbufx(0) - 1)/irat + 1
             lwrk1  = lwork  - kend1
             if (lwrk1 .lt.0) 
     &         call quit('Insufficient work space in CCSDR12AO')
             timintr12 = timintr12 + ( second() - dtime )

             if (t2tcor) then
                kt2amt = kend1
                kend1  = kt2amt + nt2sq(1)
                lwrk1  = lwork  - kend1
                if (lwrk1.lt.0) call quit('Insuff. core in CCSDR12AO')
                call dcopy(nt2sq(1),t2am,1,work(kt2amt),1)
                call ccsd_t2tp(work(kt2amt),work(kend1),lwrk1,isym0)
             end if
          else
             numdis = 1
          endif
  
c         ------------------------------------------------------
c         loop over AOs delta included in the distribution
c         calculated in the above call to the integrals program:
c         ------------------------------------------------------
          do idel2 = 1,numdis
  
            if (direct) then
               idel  = indexa(idel2)
               isymd = isao(idel)
            else
               idel  = ibas(isymd1) + ibasx(isymd1) + illl
               isymd = isymd1
            endif
            
            isydis = isymd

C           ------------------------------------------------------------
C           record number and start addresses for intermediates on file:
C           append the records for auxiliary basis functions after those
C           where delta is primary basis function
C           ------------------------------------------------------------
            irecord = mbas1t + idel - mbas1(isymd) - ibas(isymd)
            if (locdbg) 
     &        write(lupri,*) 'in CCSDR12AO: irecord = ',irecord 
            it2del(irecord) = icdel1
            icdel1 = icdel1 + nt2bcd(isydis)
  
C           --------------------------------------------------
C           allocate memory for 3-index batch of AO integrals
C           and read batch into memory:
C           --------------------------------------------------
            kxint  = kend1
            kend2  = kxint + ndisao(isydis)
            lwrk2  = lwork  - kend2
            if (lwrk2 .lt. 0) 
     &        call quit('Insufficient work space in CCSDR12AO')
  
            dtime   = second()
            call ccrdao(work(kxint),idel,idel2,work(kend2),lwrk2,
     &                  work(krecnr),direct)
            dtime   = second() - dtime
            timrdao = timrdao  + dtime
  
*----------------------------------------------------------------------*
* compute integrals (ia|delta' j) and (ij|delta' a), where i is obtained
* by transformation with lambda-particle and j and a by transformation
* with lambda-hole. the integtrals are saved on disk.
*----------------------------------------------------------------------*
            if (ccsdr12 .and. (ianr12.eq.2.or.ianr12.eq.3)) then
               leniaj = nt2bcd(isydis)

               kxiadj = kend2
               kxijda = kxiadj + leniaj
               kend3  = kxijda + leniaj
               lwrk3  = lwork  - kend3
               if (lwrk3 .lt. 0) 
     &           call quit('Insufficient space in CCSDR12AO')

               call dzero(work(kxiadj),leniaj)
               call dzero(work(kxijda),leniaj)
CCN
C              write(lupri,*) 'in CCSDR12AO:'
C              write(lupri,*) 'Norm^2 of kxint = ',
C    &           ddot(ndisao(isydis),work(kxint),1,work(kxint),1)
CCN

               do isygam = 1, nsym
                  isyalbe = muld2h(isydis,isygam)
               do g = 1, mbas1(isygam)
                  igam = g + ibas(isygam) + ibasx(isygam) 
              
                  koffg = kxint + idsaog(isygam,isydis)
     &                    + nnbst(isyalbe)*(g-1)

                  call cc_iajb( work(koffg), isyalbe, dummy, isym0, 
     &                          idel, igam, lauxg, ibasx,
     &                          dummy, work(kxiadj), work(kxijda),
     &                          dummy, dummy, dummy, 
     &                          xlambdap, xlambdah, isym0,
     &                          dummy, dummy, isym0,
     &                          xlambdap, xlambdah, isym0,
     &                          dummy, dummy, isym0,
     &                          work(kend3), lwrk3,   3,    
     &                          .false., .false.,  .true.,   
     &                          .false., .false.,  0      )
               end do
               end do

c              ------------------------------------
c              update Fhat_{del a}:
c              ------------------------------------
               ibasd = idel - mbas1(isymd) - ibas(isymd) - ibasx(isymd)
               if (locdbg) then
                 write(lupri,*) 'in CCSDR12AO: isymd, ibasd = ',
     &                                         isymd, ibasd
               end if

               call cc_fckdela(ibasd,isymd,work(kfckvaop),isyfck,
     &                         work(kxijda),work(kxiadj),iviraop)

c              ------------------------------------
c              transform (ia|del j) to L(ia|del j):
c              ------------------------------------
               call dscal(leniaj, two,work(kxiadj),1)
               call daxpy(leniaj,-one,work(kxijda),1,work(kxiadj),1)
               
c              --------------------------------------------
c              write 3-index transformed integrals to disk:
c              --------------------------------------------
               iadr = it2del(irecord) + 1
               call putwa2(luiadj,fniadj,work(kxiadj),iadr,leniaj)
               call putwa2(luijda,fnijda,work(kxijda),iadr,leniaj)

               if (locdbg) then
                 write(lupri,*) 'in CCSDR12AO: iadr = ', iadr
                 write(lupri,*) 'in CCSDR12AO: Norm^2 of XIADJ: ',
     &              ddot(leniaj,work(kxiadj),1,work(kxiadj),1)
                 write(lupri,*) 'in CCSDR12AO: Norm^2 of XIJDA: ',
     &              ddot(leniaj,work(kxijda),1,work(kxijda),1)
               end if

            end if

*----------------------------------------------------------------------*
* Compute the C' and D' intermediates with delta' an auxiliary
* basis function using CCRHS_C and CCRHS_D. These routines require
* the amplitudes as full square matrix on t2am
*----------------------------------------------------------------------*
            if (ccsdr12 .and. (ianr12.eq.2.or.ianr12.eq.3)) then
 
C              -------------------------------------------------------
C              transform gamma index in the integral batch to occupied
C              using the lambda-partical = CMO coefficients:
C              -------------------------------------------------------
               kdsrhf = kend2
               kend3  = kdsrhf + ndsrhf(isymd)
               lwrk3  = lwork  - kend3
               if (lwrk3 .lt. 0) 
     &           call quit('Insufficient space in CCSDR12AO')

               dtime   = second()
               call cctrbt(work(kxint),work(kdsrhf),xlambdap,
     *                     isym0,work(kend3),lwrk3,isydis)
               dtime   = second() - dtime
               timtrbt = timtrbt + dtime
 

C              -----------------------------
C              calculate the C intermediate:
C              -----------------------------
               dtime = second()
               factc = xmone
               icon  = 2
               ioptr12 = 0 !only calculate ONE C-intermediate
               iopte = 1
               iv    = 1
               if (.not. t2tcor) then
                  call ccrhs_c(work(kxint),work(kdsrhf),dummy,
     *                         t2am,isym0,xlambdap,dummy,
     *                         xlambdah,xlambdap,isym0,
     *                         xlambdap,isym0,
     *                         dummy,e1pim,work(kend3),lwrk3,
     *                         irecord,isymd,factc,icon,ioptr12,iopte,
     *                         lucp,cpfil,idummy,dummy,iv)
               else
                  call ccrhs_c(work(kxint),work(kdsrhf),dummy,
     *                         work(kt2amt),isym0,
     *                         xlambdap,dummy,
     *                         xlambdah,xlambdap,isym0,
     *                         xlambdap,isym0,
     *                         dummy,e1pim,work(kend3),lwrk3,
     *                         irecord,isymd,factc,icon,ioptr12,iopte,
     *                         lucp,cpfil,idummy,dummy,iv)
               end if
               dtime   = second() - dtime
               timc    = timc     + dtime
 
c              -------------------------
c              transform T2 to 2T2 - T2.
c              -------------------------
               dtime   = second()
               if (t2tcor) then
                  call dscal(nt2sq(1),two,t2am,1)
                  call daxpy(nt2sq(1),-one,work(kt2amt),1,t2am,1)
               else
                  call ccrhs_t2tr(t2am,work(kend3),lwrk3,isym0)
               end if
               dtime   = second() - dtime
               timt2tr = timt2tr  + dtime
 
c              -----------------------------
c              calculate the D intermediate:
c              -----------------------------
               dtime = second()
               factd = one
               icon  = 2
               ioptr12 = 0
               iopte = 1
               iv    = 1
               call ccrhs_d(work(kxint),work(kdsrhf),dummy,
     *                      t2am,isym0,
     *                      xlambdap,dummy,
     *                      xlambdah,xlambdap,isym0,
     *                      xlambdah,isym0,
     *                      dummy,e1pim,work(kend3),lwrk3,
     *                      irecord,isymd,factd,icon,ioptr12,iopte,
     *                      ludp,dpfil,idummy,dummy,iv)
               dtime   = second() - dtime
               timd    = timd     + dtime
 

c              -------------------------
c              restore T2 from 2T2 - T2:
c              -------------------------
               dtime   = second()
               if (t2tcor) then
                  call daxpy(nt2sq(1),one,work(kt2amt),1,t2am,1)
                  call dscal(nt2sq(1),half,t2am,1)
               else
                  call ccrhs_t2bt(t2am,work(kend3),lwrk3,isym0)
               end if
               dtime   = second() - dtime
               timt2bt = timt2bt  + dtime

            end if

          end do ! idel2
        end do ! illl
      end do ! isymd1

*----------------------------------------------------------------------*
*     end of the loop over integral distributions...
*     restore default settings for the integral evaluation:
*----------------------------------------------------------------------*

      mbsmax = 4
      loopdp = .false.
      direct = temp_direct

      ioff   = 0
      do isym = 1,nsym
        do i = 1,nbas(isym)
          ioff = ioff + 1
          isao(ioff) = isym
        end do
      end do

*----------------------------------------------------------------------*
*     save Fhat(a,del') on file:
*----------------------------------------------------------------------*
      lunit = -1
      call gpopen(lunit,'CCFHATADEL','UNKNOWN',' ','UNFORMATTED',
     &            idummy,.false.)
      read(lunit)
      write(lunit) (work(kfckvaop-1+i),i=1,nviraop(isyfck))
      call gpclose(lunit,'KEEP')

*----------------------------------------------------------------------*
*     return
*----------------------------------------------------------------------*
      call qexit('ccsdr12ao')
      return
      end
*----------------------------------------------------------------------*
*                     END OF SUBROUTINE CCSDR12AO                      *
*======================================================================*
