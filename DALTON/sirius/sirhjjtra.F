#define HJJTRA_DEBUG 1
!
!...   Copyright (c) 2014 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2014 (2014), see http://daltonprogram.org"
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

!  FILE: DALTON/sirius/sirhjjtra.F
!  (c) Copyright Hans Joergen Aa. Jensen, hjj@sdu.dk (2014)

      subroutine SIR_INTOPEN_HJJTRA()
!
!  Written by Hans Joergen Aa. Jensen October 2014
!
      implicit none
#include "priunit.h"

      write(LUPRI,*) 'SIR_INTOPEN_HJJTRA called'

      end
! ===================================================================
!     sirhjjtra.F section 1: calculate MO integrals
! ===================================================================
      subroutine HJJTRA_CTL(ITRSIR, CMO, WORK, LWORK)
!
!  Written by Hans Joergen Aa. Jensen October 2014
!
      implicit none
      integer ITRSIR, LWORK
      real*8  CMO(*), WORK(LWORK)
#include "priunit.h"

! Used from include files:
!  gnrinf.h : PARCAL
!  inforb.h : NNORBX, NORBT, ???
!  inftra.h : IPRTRA
!  inftap.h : LUINTM, LBINTM
#include "gnrinf.h"
#include "inforb.h"
#include "inftra.h"
#include "inftap.h"

      integer JTRSIR, NTRLVL, MX_NDMAT
      integer KFREE, LFREE, KFRSAV, KUCMO, KDMAT, KFMAT
      integer ICDTRA(NNORBX), ITRTYP(NORBT) ! automatic arrays
      character*8 TABLE(2),LAB123(3)
      data        TABLE /'MO2ELHJJ','END INTM'/
      data        LAB123/'********','********','********'/

C --- start of executable code

      call QENTER('HJJTRA_CTL')

      KFREE  = 1
      LFREE  = LWORK
      KFRSAV = KFREE
     
      ! call SIR_INTOPEN_HJJTRA()

!     translate sirius integral level code to HJJTRA and NEWTRA integral level code
      JTRSIR = abs(ITRSIR)
      call TRALVL(JTRSIR, NTRLVL) 
#ifdef HJJTRA_DEBUG
      write(LUPRI,*) 'HJJTRA_CTL called'
      write(LUPRI,*) 'ITRSIR, NTRLVL, LWORK:',ITRSIR, NTRLVL, LWORK
#endif

!     call sirntra.F(TRASET) in order to be able to call sirntra.F(TRALIM)
!     which sets :
!        ICDTRA: index array for C,D distributions (**|CD) to calculate
!        ITRTYP(1:NORBT) = number of integral indices in which this orbital
!                          enters (i.e. 0,1,2,3, or 4)
      call TRASET(-1,LWORK)
#ifdef HJJTRA_DEBUG
      IPRTRA = MAX(IPRTRA,5)
#endif
      call TRALIM(NTRLVL,ICDTRA,ITRTYP)

      CALL REWSPL(LUINTM)
      READ  (LUINTM) ! info not used for HJJTRA
      WRITE (LUINTM) LBINTM,ITRSIR,NSYM,NORB,NBAS
      WRITE (LUINTM) CMO(1:NCMOT)
      ! WRITE (LUINTM) MISH(1:8), MASH(1:8)
      WRITE (LUINTM) NISH(1:8), NASH(1:8)

      call GETDAT(LAB123(2),LAB123(3))   ! place date in LAB123(2) and time in LAB123(3)
      WRITE (LUINTM) LAB123,TABLE(1)

!     if (PARCAL) then
!        CALL QUIT('HJJTRA_CTL parcal not finished yet')
!        ! CALL HJJTRA_PARCAL(ICDTRA, ITRTYP, CMO, WORK, LWORK)
!     else ! not PARCAL

         MX_NDMAT = 10 ! to be optimized TODO
         call MEMGET2('REAL','UCMO',  KUCMO,NBAS*NORBT,
     &      WORK,KFREE,LFREE)
         call MEMGET2('REAL','DMAT',  KDMAT,MX_NDMAT*N2BASX,
     &      WORK,KFREE,LFREE)
         call MEMGET2('REAL','FMAT',  KFMAT,MX_NDMAT*N2BASX,
     &      WORK,KFREE,LFREE)

         call UPKCMO(CMO,WORK(KUCMO))
         call HJJTRA_SEQUENTIAL(ICDTRA, ITRTYP, WORK(KUCMO),
     &      MX_NDMAT, WORK(KDMAT), WORK(KFMAT), WORK(KFREE), LFREE)

         call MEMREL('after HJJTRA_SEQUENTIAL',WORK,
     &      1,KFRSAV,KFREE,LFREE)

!     end if ! if (PARCAL) then ... else ...

      WRITE (LUINTM) LAB123,TABLE(2)
      CALL REWSPL(LUINTM)

      call QEXIT('HJJTRA_CTL')
      end

      subroutine HJJTRA_SEQUENTIAL(
     &     ICDTRA, ITRTYP, UCMO, MX_NDMAT, DMAT, FMAT, WORK, LWORK)
!
!  Written by Hans Joergen Aa. Jensen October 2014
!
      implicit none
      integer ICDTRA(NNORBX), ITRTYP(NORBT), MX_NDMAT, LWORK
      real*8  DMAT(N2BASX,MX_NDMAT), FMAT(N2BASX,MX_NDMAT)
      real*8  UCMO(NBAST,NORBT), WORK(LWORK)
#include "priunit.h"

! Used from include files:
!  infinp.h : DIRFCK
!  inforb.h : N2BASX, NORBT,NNORBX, ...
!  infind.h : ISMO(:)
#include "maxorb.h"
#include "maxash.h"
#include "infinp.h"
#include "inforb.h"
#include "infind.h"

      integer ICD, IC, ID, NDMAT
      integer ICD_DMAT(2,MX_NDMAT), IFCTYP(MX_NDMAT), ISYMDM(MX_NDMAT) ! automatic arrays

      write(LUPRI,*) 'HJJTRA_SEQUENTIAL called'
      write(LUPRI,*) 'LWORK:',LWORK

      IFCTYP(1:MX_NDMAT) = 11 ! symmetric density matrix, only Coulomb

      NDMAT = 0
      ICD = 0
      do IC = 1, NORBT
         do ID = 1, IC
            ICD = ICD + 1
         if (ICDTRA(ICD) .le. 0) cycle
            NDMAT = NDMAT + 1
            ISYMDM(NDMAT) = MULD2H(ISMO(IC), ISMO(ID))
            ICD_DMAT(1,NDMAT) = IC
            ICD_DMAT(2,NDMAT) = ID
            call HJJTRA_CD_DMAT(IC,ID,UCMO,DMAT(1,NDMAT))
            if (NDMAT .eq. MX_NDMAT .or. ICD .eq. NNORBX) then
               call SIRFCK(FMAT,DMAT,NDMAT,ISYMDM,IFCTYP,DIRFCK,
     &                     WORK,LWORK)
               call HJJTRA_WR_LUINTM(UCMO,NDMAT,FMAT,DMAT,
     &                     ICD_DMAT,ISYMDM)
               NDMAT = 0
            end if
         end do
      end do

      if (NDMAT .ne. 0) CALL QUIT('Programming error, NDMAT .ne. 0')
      
      end
! ===================================================================
      subroutine HJJTRA_CD_DMAT(IC,ID,UCMO,DMAT)
!
!  Written by Hans Joergen Aa. Jensen October 2014
!
      implicit none
      integer IC, ID
      real*8  UCMO(NBAST,NORBT), DMAT(NBAST,NBAST)
#include "priunit.h"

! Used from include files:
!  inforb.h : NORBT,NBAAT
#include "inforb.h"

      integer JA, JB
      real*8  UCMO_BC, UCMO_BD

      do JB = 1, NBAST
         UCMO_BC = UCMO(JB,IC)
         UCMO_BD = UCMO(JB,ID)
         do JA = 1,NBAST
            DMAT(JA,JB) = 0.5D0*
     &         (UCMO(JA,IC)*UCMO_BD + UCMO(JA,ID)*UCMO_BC) ! symmetrize
         end do
      end do

      call DGETSI(NBAST,DMAT,DMAT) ! symmetrize

      end
! ===================================================================
      subroutine HJJTRA_WR_LUINTM(UCMO,NDMAT,FMAT,DMAT,ICD_DMAT,ISYMDM)
!
!  Written by Hans Joergen Aa. Jensen October 2014
!
      implicit none
      integer NDMAT, ICD_DMAT(2,NDMAT), ISYMDM(NDMAT)
      real*8  UCMO(NCMOT), DMAT(N2BASX,NDMAT), FMAT(N2BASX,NDMAT)
#include "priunit.h"

! Used from include files:
!  inforb.h : NCMOT,N2BASX, ...
!  inftap.h : LUINTM
#include "inforb.h"
#include "inftap.h"

      integer IDMAT, IC, ID

      do IDMAT = 1, NDMAT
         IC = ICD_DMAT(1,IDMAT)
         ID = ICD_DMAT(2,IDMAT)
      end do

      end

! ===================================================================
!     sirhjjtra.F section 2: read MO integrals
! ===================================================================
      subroutine HJJTRA_NXTH2M(IC,ID,H2CD,NEEDTP,IDIST)
!
!  Written by Hans Joergen Aa. Jensen October 2014
!  This version is interface routine for .HJJTRA integral transformation.
!
! Purpose:
!    Read next Mulliken two-electron integral distribution (**|cd)
!    where (cd) distribution is needed according to NEEDTP(ITYPCD)
!
! Usage:
! 1) Set IDIST = 0 before first call.
!    Do NOT change IDIST in calling routine
!    until last distribution has been read (signalled by IDIST .eq. -1)
! 2) if IDIST<0, then IC, ID must be set to correspond to the distribution
!    wanted in H2CD
!

      implicit none
#include "priunit.h"

      integer  IC, ID, NEEDTP(*), IDIST
      real*8   H2CD(NORBT,NORBT)

! Used from include files:
!  inforb.h : NORBT, ???
#include "inforb.h"

      write(LUPRI,*) 'HJJTRA_NXTH2M called, IC, ID, IDIST =',
     &   IC, ID, IDIST

      if (IDIST .le. 0) then
         if (IC .gt. 0 .and. IC .le. NORBT .and.
     &       ID .gt. 0 .and. ID .le. NORBT) then
            write (LUPRI,*) 'HJJTRA_NXTH2M: invalid input IC, ID',IC,ID
            call quit('HJJTRA_NXTH2M: invalid input IC, ID')
         end if
      else if (IDIST .eq. 0) then
      else
      end if

      end

      subroutine HJJTRA_NXTH2D(IC,ID,H2CD,NEEDTP,IDIST)
!
!  Written by Hans Joergen Aa. Jensen October 2014
!  This version is interface routine for .HJJTRA integral transformation.
!  
! Purpose:
!    Read next Dirac two-electron integral distribution <**|cd> = (*c|*d)
!    where (cd) distribution is needed according to NEEDTP(ITYPCD)
!
!

      implicit none
#include "priunit.h"

      integer  IC, ID, NEEDTP(*), IDIST
      real*8   H2CD(NORBT,NORBT)

! Used from include files:
!  inforb.h : NORBT, ???
#include "inforb.h"

      call quit('ERROR: HJJTRA_NXTH2D is not implemented')

      end
! -- end of DALTON/sirius/sirhjjtra.F --
