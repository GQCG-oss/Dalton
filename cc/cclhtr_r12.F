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

*=====================================================================*
      subroutine CCLHTR_FP(RHO1,XLAMDH,WORK,LWORK,ISYMC,FC12AM,LUFC12,
     &           IFILE)
*---------------------------------------------------------------------*
C     purpose: Calculate Sum{klj}c^{ij}_{kl}*V^{a jt}_{kl} =: rho_Fp
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

      do isymak = 1, nsym
         nvajkl(isymak) = 0
         do isymk = 1, nsym
            isyma = muld2h(isymak,isymk)
            nvajkl(isymak) = nvajkl(isymak) + nt1ao(isyma)*nmatij(isymk)
         end do
      end do

      kvec = 1
      kvajkl = kvec + ngamsq(isymc)
      kend1 = kvajkl + nvajkl(1)
      lwrk1 = lwork - kend1

      if (lwrk1 .lt.0) then
         call quit('Insufficient work space in cclhtr_fp')
      end if

C     -----------------------
C     read V^{alpha jt}_{kl}
C     -----------------------
      luvajkl=0
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
      call cc_r12getct(work(kvec),isymc,1,2.0D0*brascl,.false.,.true.,
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
         call cc_prsqr12(rhor12,isyres,1,.false.)
         write(lupri,*) 'norm^2(rhor12):',
     &     ddot(ngamsq(isyres),rhor12,1,rhor12,1)

         write(lupri,*) 'Leaving CCLHTR_GP'
      end if

      call qexit('cclhtr_gp')
      return
      end
*=====================================================================*
