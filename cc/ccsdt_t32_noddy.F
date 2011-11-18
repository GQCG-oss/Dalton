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
      subroutine ccsdt_t32_noddy(tab3am,listab,idlstab,freqab,
     &                           xint1s0,xint2s0,
     &                           xlamdp0,xlamdh0,fock0,fockd,
     &                           fieldao,field,
     &                           scr1,work,lwork)
*---------------------------------------------------------------------*
*    Purpose: compute triples part of second-order response amplitudes
*            
*    Input:  listab, idlstab
*            xlamdp0, xlamdh0, fock0, fockd
*            xint1s, xint2s
*            work, lwork
*
*    Output: tab3am, freqab
*
*    Written by Christof Haettig, Mai 2003, based on ccsdt_t31_noddy.
*---------------------------------------------------------------------*
      implicit none
#include <priunit.h>
#include <ccsdinp.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <cco2rsp.h>
#include <ccr2rsp.h>
#include <ccer1rsp.h>
#include <ccexci.h>
#include <ccfield.h> 
#include <dummy.h> 


      logical locdbg, print_t3
      parameter (locdbg=.false., print_t3=.false.)

      integer isym0
      parameter ( isym0 = 1 )

      character*3 listab
      integer idlstab, lwork

#if defined (SYS_CRAY) 
      real tab3am(*), freqab, freqa, freqb
      real scr1(*), work(*), xlamdp0(*), xlamdh0(*)
      real fock0(*), fockd(*), xint1s0(*), xint2s0(*)
      real field(*), fieldao(*)
      real two, one, ddot
#else
      double precision tab3am(*), freqab, freqa, freqb
      double precision scr1(*), work(*), xlamdp0(*), xlamdh0(*)
      double precision fock0(*), fockd(*), xint1s0(*), xint2s0(*)
      double precision field(*), fieldao(*)
      double precision two, one, ddot
#endif
      parameter ( one = 1.0d0, two = 2.0d0 )

      character lista*3, listb*3, labela*8, labelb*8, model*10
      logical lorxa, lorxb, rhs_only
      integer isyma, isymb, idlsta, idlstb, ierr, irrep, kend1, lwrk1,
     &        kta1am, kta2am, klampa, klamha, kfocka_ao, kfocka,
     &        ktb1am, ktb2am, klampb, klamhb, kfockb_ao, kfockb,
     &        kint1sa, kint2sa, kint1sb, kint2sb, kint1sab, kint2sab,
     &        kta3am, ktb3am, kt03am, iopt, ktab1am, ktab2am, isymab,
     &        kt02am, kfockab, kfckbuf

* external functions:
      integer ir1tamp, ilstsym

      call qenter('CCT32NOD')
*---------------------------------------------------------------------*
*     Do some initializations and checks:
*---------------------------------------------------------------------*
      if   (listab(1:3).eq.'O2 ') then
        ! we compute the rhs vectors for second-order amplitude resp.
        rhs_only = .true.

        lista  = 'R1 '
        labela = lblo2(idlstab,1)
        isyma  = isyo2(idlstab,1)
        freqa  = frqo2(idlstab,1)
        lorxa  = lorxo2(idlstab,1)
        idlsta = ir1tamp(labela,lorxa,freqa,isyma)

        listb  = 'R1 '
        labelb = lblo2(idlstab,2)
        isymb  = isyo2(idlstab,2)
        freqb  = frqo2(idlstab,2)
        lorxb  = lorxo2(idlstab,2)
        idlstb = ir1tamp(labelb,lorxb,freqb,isymb)

      else if   (listab(1:3).eq.'R2 ') then
        ! we compute the second-order amplitude response vectors
        rhs_only = .false.

        lista  = 'R1 '
        labela = lblr2t(idlstab,1)
        isyma  = isyr2t(idlstab,1)
        freqa  = frqr2t(idlstab,1)
        lorxa  = lorxr2t(idlstab,1)
        idlsta = ir1tamp(labela,lorxa,freqa,isyma)

        listb  = 'R1 '
        labelb = lblr2t(idlstab,2)
        isymb  = isyr2t(idlstab,2)
        freqb  = frqr2t(idlstab,2)
        lorxb  = lorxr2t(idlstab,2)
        idlstb = ir1tamp(labelb,lorxb,freqb,isymb)

      else if (listab(1:3).eq.'EO1' .or. listab(1:3).eq.'ER1') then
        ! ER1 - first-order response of excited states
        ! EO1 - rhs vectors for ER1 equations
        rhs_only = listab(1:3) .eq. 'EO1'

        lista  = 'R1 '
        labela = lbler1(idlstab)
        isyma  = isyoer1(idlstab)
        freqa  = frqer1(idlstab)
        lorxa  = lorxer1(idlstab)
        idlsta = ir1tamp(labela,lorxa,freqa,isyma)

        listb  = 'RE '
        labelb = '-- XX --'
        isymb  = isyser1(idlstab)
        freqb  = eiger1(idlstab)
        lorxb  = .false.
        idlstb = ister1(idlstab)

      else
        call quit('Unknown LISTAB in CCSDT_T32_NODDY.')
      end if

      isymab = ilstsym(listab,idlstab)
      freqab = freqa + freqb
 

      if (lorxa.or.lorxb) then
        call quit('Orbital relaxation not allowed in CCSDT_T32_NODDY.')
      end if

      if (listb(1:3).ne.'R1 ' .and. listb(1:3).NE.'RE ') then
        call quit('Unknown LISTB in CCSDT_T32_NODDY.')
      end if

      if (lista(1:3).ne.'R1 ' .and. lista(1:3).NE.'RE ') then
        call quit('Unknown LISTA in CCSDT_T32_NODDY.')
      end if

      call dzero(tab3am,nt1amx*nt1amx*nt1amx)

*---------------------------------------------------------------------*
*     Memory allocation:
*---------------------------------------------------------------------*
      kend1   = 1

      kta1am  = kend1
      klampa  = kta1am + nt1amx
      klamha  = klampa + nlamdt
      kend1   = klamha + nlamdt

      ktb1am  = kend1
      klampb  = ktb1am + nt1amx
      klamhb  = klampb + nlamdt
      kend1   = klamhb + nlamdt

      kt02am  = kend1
      kta2am  = kt02am + nt1amx*nt1amx
      ktb2am  = kta2am + nt1amx*nt1amx
      kend1   = ktb2am + nt1amx*nt1amx

      if (lista.eq.'R1 ') then
        kfocka_ao = kend1
        kfocka    = kfocka_ao + nbast*nbast
        kend1     = kfocka    + nbast*nbast
      else
        kfocka_ao = -999 999
        kfocka    = -999 999
      end if

      if (listb.eq.'R1 ') then
        kfockb_ao = kend1
        kfockb    = kfockb_ao + nbast*nbast
        kend1     = kfockb    + nbast*nbast
      else
        kfockb_ao = -999 999
        kfockb    = -999 999
      end if

      if (lista.eq.'R1 '.or.listb.eq.'R1 ') then
        kfockab = kend1
        kfckbuf = kfockab + nbast*nbast
        kend1   = kfckbuf + nbast*nbast
      else
        kfockab = -999 999
        kfckbuf = -999 999
      end if

      kint1sa = kend1
      kint2sa = kint1sa + nt1amx*nvirt*nvirt
      kend1   = kint2sa + nrhft*nrhft*nt1amx
 
      kint1sb = kend1
      kint2sb = kint1sb + nt1amx*nvirt*nvirt
      kend1   = kint2sb + nrhft*nrhft*nt1amx
 
      kint1sab = kend1
      kint2sab = kint1sab + nt1amx*nvirt*nvirt
      kend1    = kint2sab + nrhft*nrhft*nt1amx

      kta3am = kend1
      kend1  = kta3am + nt1amx*nt1amx*nt1amx
      ktb3am = kta3am
      kt03am = kta3am

      lwrk1  = lwork  - kend1
      if (lwrk1 .lt. 0) then
         call quit('Insufficient space in CCSDT_T32_NODDY')
      endif

*---------------------------------------------------------------------*
*     Read zeroth-order amplitudes and response vectors into memory
*---------------------------------------------------------------------*
      if (lwrk1.lt.nt2amx) then
        call quit('Insufficient space in CCSDT_T32_NODDY.')
      end if

      iopt = 2
      call cc_rdrsp('R0 ',0,isym0,iopt,model,dummy,work(kend1))  
      call cc_t2sq(work(kend1),work(kt02am),isym0)

      iopt = 3
      call cc_rdrsp(lista,idlsta,isyma,iopt,model,
     &              work(kta1am),work(kend1))  
      call cclr_diascl(work(kend1),two,isyma)
      call cc_t2sq(work(kend1),work(kta2am),isyma)

      iopt = 3
      call cc_rdrsp(listb,idlstb,isymb,iopt,model,
     &              work(ktb1am),work(kend1))  
      call cclr_diascl(work(kend1),two,isymb)
      call cc_t2sq(work(kend1),work(ktb2am),isymb)

*---------------------------------------------------------------------*
*     Get property matrices A and B and the matrix A^B+B^A:
*---------------------------------------------------------------------*
      if (lista(1:3).eq.'R1 ') then
        ! read property integrals from file:
        call ccprpao(labela,.TRUE.,work(kfocka_ao),irrep,isyma,ierr,
     &               work(kend1),lwrk1)
        if ((ierr.gt.0) .or. (ierr.eq.0 .and. irrep.ne.isyma)) then
          call quit('CCSDT_T32_NODDY: error reading operator '//LABELA)
        else if (ierr.lt.0) then
          call dzero(work(kfocka_ao),n2bst(isyma))
        end if
        call dcopy(nbast*nbast,work(kfocka_ao),1,work(kfocka),1)
 
        ! transform property integrals to Lambda-MO basis
        call cc_fckmo(work(kfocka),xlamdp0,xlamdh0,
     &                work(kend1),lwrk1,isyma,1,1)
      end if

      if (listb(1:3).eq.'R1 ') then
        ! read property integrals from file:
        call ccprpao(labelb,.TRUE.,work(kfockb_ao),irrep,isymb,ierr,
     &               work(kend1),lwrk1)
        if ((ierr.gt.0) .or. (ierr.eq.0 .and. irrep.ne.isymb)) then
          call quit('CCSDT_T32_NODDY: error reading operator '//LABELB)
        else if (ierr.lt.0) then
          call dzero(work(kfockb_ao),n2bst(isymb))
        end if
        call dcopy(nbast*nbast,work(kfockb_ao),1,work(kfockb),1)
 
        ! transform property integrals to Lambda-MO basis
        call cc_fckmo(work(kfockb),xlamdp0,xlamdh0,
     &                work(kend1),lwrk1,isymb,1,1)
      end if


*---------------------------------------------------------------------*
*     Compute contributions
*
*            <mu_3| [B,T^A_3] + [[B,T^A_2],T^0_2] |HF>
*
*     Compute corrections to triples vector T^A_3, and corresponding 
*     lambda matrices and the XINT1SA,XINT2SA integrals and set FREQA,
*---------------------------------------------------------------------*
      if (listb(1:3).eq.'R1 ') then

        if (lista(1:3).eq.'R1 ' .or. lista(1:3).eq.'RE ' .or.
     &      lista(1:3).eq.'RC '                              ) then
          call ccsdt_t31_noddy(work(kta3am),lista,idlsta,freqa,.false.,
     &                         .false.,xint1s0,xint2s0,
     &                         .false.,dummy,dummy,
     &                         .false.,dummy,dummy,
     &                         work(kint1sa),work(kint2sa),
     &                         work(klampa),work(klamha),work(kfocka),
     &                         xlamdp0,xlamdh0,fock0,dummy,fockd,
     &                         work(kend1),lwrk1)
          call dscal(nt1amx*nt1amx*nt1amx,-1.0d0,work(kta3am),1)
        else
          call quit('Unknown LISTA in CCSDT_T32_NODDY.')
        end if

c     write(lupri,*) 'norm^2(tab3arm) before cont. of t^B_3:',
c    &   ddot(nt1amx**3,tab3am,1,tab3am,1)
        call ccsdt_xksi3_2(tab3am,work(kfockb),work(kta3am))
c     write(lupri,*) 'norm^2(tab3arm) after cont. of t^B_3-1:',
c    &   ddot(nt1amx**3,tab3am,1,tab3am,1)

        call ccsdt_xksi3_1(tab3am,work(kfockb),
     &                     work(kt02am),work(kta2am),one)
        call ccsdt_xksi3_1(tab3am,work(kfockb),
     &                     work(kta2am),work(kt02am),one)
c     write(lupri,*) 'norm^2(tab3arm) after cont. of t^B_3-2:',
c    &   ddot(nt1amx**3,tab3am,1,tab3am,1)

      else 
        ! we don't need the triples vector, but still we need in
        ! the following the response Lambda matrices and the
        ! one-index transformed integrals

        call cclr_lamtra(xlamdp0,work(klampa),xlamdh0,work(klamha),
     &                   work(kta1am),isyma)

        call ccsdt_ints1_noddy(.true.,work(kint1sa),work(kint2sa),
     &                         .false.,dummy,dummy,
     &                         xlamdp0,xlamdh0,
     &                         work(klampa),work(klamha),
     &                         work(kend1),lwrk1) 

      end if

      if (locdbg) then
        write(lupri,*) 'norm^2(tab3am) after cont. of t^B_3:',
     &     ddot(nt1amx**3,tab3am,1,tab3am,1)
      end if
 
*---------------------------------------------------------------------*
*     Compute the contributions 
*
*            <mu_3| [A,T^B_3] + [[A,T^B_2],T^0_2] |HF>
*
*     Compute corrections to triples vector T^B_3, and corresponding 
*     lambda matrices and the XINT1SB,XINT2SB integrals and set FREQB:
*---------------------------------------------------------------------*
      if (lista(1:3).eq.'R1 ') then

        if (listb(1:3).eq.'R1 ' .or. listb(1:3).eq.'RE ' .or.
     &      listb(1:3).eq.'RC '                              ) then
          call ccsdt_t31_noddy(work(ktb3am),listb,idlstb,freqb,.false.,
     &                         .false.,xint1s0,xint2s0,
     &                         .false.,dummy,dummy,
     &                         .false.,dummy,dummy,
     &                         work(kint1sb),work(kint2sb),
     &                         work(klampb),work(klamhb),work(kfockb),
     &                         xlamdp0,xlamdh0,fock0,dummy,fockd,
     &                         work(kend1),lwrk1)
          call dscal(nt1amx*nt1amx*nt1amx,-1.0d0,work(ktb3am),1)
        else
          call quit('Unknown LISTB in CCSDT_T32_NODDY.')
        end if

        call ccsdt_xksi3_2(tab3am,work(kfocka),work(ktb3am))

        call ccsdt_xksi3_1(tab3am,work(kfocka),
     &                     work(kt02am),work(ktb2am),one)
        call ccsdt_xksi3_1(tab3am,work(kfocka),
     &                     work(ktb2am),work(kt02am),one)

      else 
        ! we don't need the triples vector, but still we need in
        ! the following the response Lambda matrices and the
        ! one-index transformed integrals

        call cclr_lamtra(xlamdp0,work(klampb),xlamdh0,work(klamhb),
     &                   work(ktb1am),isymb)

        call ccsdt_ints1_noddy(.true.,work(kint1sb),work(kint2sb),
     &                         .false.,dummy,dummy,
     &                         xlamdp0,xlamdh0,
     &                         work(klampb),work(klamhb),
     &                         work(kend1),lwrk1) 

      end if

      if (locdbg) then
        write(lupri,*) 'norm^2(tab3arm) after cont. of t^A_3:',
     &     ddot(nt1amx**3,tab3am,1,tab3am,1)
      end if
 
*---------------------------------------------------------------------*
*     Compute contributions:
*
*           <mu_3| [ ([A,T^B_1] + [B,T^A_1]) , T^0_3 ] |HF>
*       
*     Compute the matrix with one-index transformed property 
*     integrals FOCKAB = [A,T^B_1] + [B,T^A_1]
*---------------------------------------------------------------------*
      if (lista(1:3).eq.'R1 ' .or. listb(1:3).eq.'R1 ') then
        call dzero(work(kfockab),nbast*nbast)

        if (lista(1:3).eq.'R1 ') then
          ! add [A,T^B_1]
          call dcopy(nbast*nbast,work(kfocka_ao),1,work(kfckbuf),1)
          call cc_fckmo(work(kfckbuf),work(klampb),xlamdh0,
     &                  work(kend1),lwrk1,isyma,isymb,isym0)
          call daxpy(nbast*nbast,one,work(kfckbuf),1,work(kfockab),1)
          
          call dcopy(nbast*nbast,work(kfocka_ao),1,work(kfckbuf),1)
          call cc_fckmo(work(kfckbuf),xlamdp0,work(klamhb),
     &                  work(kend1),lwrk1,isyma,isym0,isymb)
          call daxpy(nbast*nbast,one,work(kfckbuf),1,work(kfockab),1)
        end if

        if (listb(1:3).eq.'R1 ') then
          ! add [B,T^A_1]
          call dcopy(nbast*nbast,work(kfockb_ao),1,work(kfckbuf),1)
          call cc_fckmo(work(kfckbuf),work(klampa),xlamdh0,
     &                  work(kend1),lwrk1,isymb,isyma,isym0)
          call daxpy(nbast*nbast,one,work(kfckbuf),1,work(kfockab),1)
          
          call dcopy(nbast*nbast,work(kfockb_ao),1,work(kfckbuf),1)
          call cc_fckmo(work(kfckbuf),xlamdp0,work(klamha),
     &                  work(kend1),lwrk1,isymb,isym0,isyma)
          call daxpy(nbast*nbast,one,work(kfckbuf),1,work(kfockab),1)
        end if

        if (nonhf .and. (nfield.gt.0) .and.
     &      (lwrk1.lt.nt1amx*nt1amx*nt1amx)) then
           call quit('Insufficient space in CCSDT_T32_NODDY')
        endif

        call dzero(work(kt03am),nt1amx*nt1amx*nt1amx) 
        call ccsdt_t03am(work(kt03am),xint1s0,xint2s0,
     &                   work(kt02am),scr1,fockd,
     &                   field,work(kend1))

        call ccsdt_xksi3_2(tab3am,work(kfockab),work(kt03am))
      end if

      if (locdbg) then
        call ccsdt_clean_t3(tab3am,nt1amx,nvirt,nrhft)
        write(lupri,*) 'norm^2(tab3am) after A{O} contrib.:',
     &     ddot(nt1amx**3,tab3am,1,tab3am,1)
        call print_pt3_noddy(tab3am)
      end if
 
*---------------------------------------------------------------------*
*     Compute double one-index transformed integrals and evaluate
*     the B matrix contributions
*
*       <mu_3| [H^AB,T^0_2] + [H^A,T^B_2] + [H^B,T^A_2] |HF>
*
*---------------------------------------------------------------------*
      call ccsdt_ints2_noddy(work(kint1sab),work(kint2sab),
     &                       xlamdp0,xlamdh0,
     &                       work(klampb),work(klamhb),
     &                       work(klampa),work(klamha),
     &                       work(kend1),lwrk1) 

      if (locdbg) then
        write(lupri,*) 'norm^2(tab3am) before b3am:',
     &     ddot(nt1amx**3,tab3am,1,tab3am,1)
      end if
 
      ! note: this routine overwrites work(kt02am)
      call ccsdt_b3am(tab3am,
     &                work(kint1sab),work(kint2sab),fockd,
     &                lista,idlsta,work(kint1sa),work(kint2sa),
     &                listb,idlstb,work(kint1sb),work(kint2sb),
     &                scr1,work(kt02am),work(kend1),lwrk1) 
        if (nonhf .and. (nfield.gt.0))
     &    call quit('Finite Field unfinished in CCSDT_T32_NODDY.')
cch
c       call dscal(nt1amx*nt1amx*nt1amx,-1.0d0,tab3am,1)
cch
c     call dscal(nt1amx**3,0.5d0,tab3am,1)

      if (locdbg) then
        call ccsdt_clean_t3(tab3am,nt1amx,nvirt,nrhft)
        write(lupri,*) 'norm^2(tab3am) after b3am:',
     &     ddot(nt1amx**3,tab3am,1,tab3am,1)
      end if
 
*---------------------------------------------------------------------*
*     if solution vector requested add contribution from jacobian:
*
*        <mu_3| [[H,T^AB_1],T^0_2] + [H,T^AB_2] |HF>
*
*     and solve the triples equations:
*      
*
*---------------------------------------------------------------------*
      if (.not. rhs_only) then
        ktab1am = ktb1am
        ktab2am = ktb2am 

        if (lwrk1.lt.max(nt2amx,nt1amx*nt1amx*nt1amx)) then
          call quit('Insufficient space in CCSDT_T32_NODDY.')
        end if

        iopt = 2
        call cc_rdrsp('R0 ',0,isym0,iopt,model,dummy,work(kend1))  
        call cc_t2sq(work(kend1),work(kt02am),isym0)

        iopt = 3
        call cc_rdrsp(listab,idlstab,isymab,iopt,model,
     &                work(ktab1am),work(kend1))  
        call cclr_diascl(work(kend1),two,isymab)
        call cc_t2sq(work(kend1),work(ktab2am),isymab)

        ! note: here the doubly one-index transformed integrals
        !       on kint1sab and kint2sab are overwritten by the
        !       one t^ab one-index transformed integrals
        call ccsdt_a3am(tab3am,work(ktab1am),work(ktab2am),
     &                  freqab,isymab,work(kt02am),work(kt03am),
     &                  xint1s0,xint2s0,
     &                  work(kint1sab),work(kint2sab),
     &                  fockd,xlamdp0,xlamdh0,
     &                  fieldao,field,scr1,work(kend1),lwrk1)
        call ccsdt_clean_t3(tab3am,nt1amx,nvirt,nrhft)
     
        if (locdbg) then
          write(lupri,*) 'norm^2(tab3am) after a3am:',
     &       ddot(nt1amx**3,tab3am,1,tab3am,1)
        end if
 
        call ccsdt_3am(tab3am,freqab,scr1,fockd,
     &                 nonhf,field,.false.,work(kend1))
 
        call dscal(nt1amx*nt1amx*nt1amx,-1.0d0,tab3am,1)
 
      end if

      call ccsdt_clean_t3(tab3am,nt1amx,nvirt,nrhft)

      if (print_t3) then
        write(lupri,*)'CCSDT_T32_AM> vector type:',listab
        write(lupri,*)'CCSDT_T32_AM> list,idlst:',listab,idlstab
        write(lupri,*)'CCSDT_T32_AM> freq:',freqab
        call print_pt3_noddy(tab3am)
      end if

      call qexit('CCT32NOD')
      return
      end 
*---------------------------------------------------------------------*
*                    end of subroutine ccsdt_t32_noddy
*---------------------------------------------------------------------*
