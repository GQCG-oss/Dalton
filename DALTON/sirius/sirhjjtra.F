!#define HJJTRA_DEBUG 1
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

!ifdef HJJTRA_DEBUG
      write(LUPRI,*) 'SIR_INTOPEN_HJJTRA called'
!#endif

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
      integer KFREE, LFREE, KFRSAV, KUCMO, KDMAT, KFMAT, I
      integer ICDTRA(NNORBX), ITRTYP(NORBT) ! automatic arrays
      character*8 TABLE(2),LAB123(3)
      data        TABLE /'MO2INTHJ','END INTM'/
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
      IF (JTRSIR .GE. 1 .AND. JTRSIR .LE. 4) NTRLVL = NTRLVL + 1
      ! we need third order, because USEDRC not defined for HJJTRA
!#ifdef HJJTRA_DEBUG
      write(LUPRI,*) 'HJJTRA_CTL called'
      write(LUPRI,*) 'ITRSIR, NTRLVL, LWORK:',ITRSIR, NTRLVL, LWORK
!#endif

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
      WRITE (LUINTM) (0, I=1,8) ! info not used for HJJTRA
      WRITE (LUINTM) LBINTM,ITRSIR,NSYM,NORB,NBAS
      WRITE (LUINTM) CMO(1:NCMOT)
      ! WRITE (LUINTM) MISH(1:8), MASH(1:8)
      WRITE (LUINTM) NISH(1:8), NASH(1:8)

      call GETDAT(LAB123(2),LAB123(3))   ! place date in LAB123(2) and time in LAB123(3)
      WRITE (LUINTM) LAB123,TABLE(1)
      WRITE (LUINTM) NNORBX, NCDTRA, ICDTRA

!     if (PARCAL) then
!        CALL QUIT('HJJTRA_CTL parcal not finished yet')
!        ! CALL HJJTRA_PARCAL(ICDTRA, ITRTYP, CMO, WORK, LWORK)
!     else ! not PARCAL

         call MEMGET2('REAL','UCMO',  KUCMO,NBAST*NORBT,
     &      WORK,KFREE,LFREE)

         MX_NDMAT = LFREE/(4*N2BASX) ! use half of free memory for DMAT and FMAT
         MX_NDMAT = MIN(NNORBX,MX_NDMAT)
         write(lupri,*) 'MX_NDMAT =',MX_NDMAT,NNORBX
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
      !implicit none
#include "implicit.h"
      integer ICDTRA(NNORBX), ITRTYP(NORBT), MX_NDMAT, LWORK
      real*8  DMAT(N2BASX,MX_NDMAT), FMAT(N2BASX,MX_NDMAT)
      real*8  UCMO(NBAST,NORBT), WORK(LWORK)
#include "priunit.h"

! Used from include files:
!  infinp.h : DIRFCK
!  inforb.h : N2BASX, NORBT,NNORBX, ...
!  infind.h : ISMO(:)
!  inftap.h : LUMINT
#include "maxorb.h"
#include "maxash.h"
#include "infinp.h"
#include "inforb.h"
#include "infind.h"
#include "inftap.h"

!DEBUG
#include "infpri.h"

#include "iratdef.h"

      integer ICD, IC, ID, NDMAT,  LBUF
      integer ICD_REC(MX_NDMAT), IFCTYP(MX_NDMAT), ISYMDM(MX_NDMAT) ! automatic arrays
      logical FEXIST, OLDDX

      write(LUPRI,*) 'HJJTRA_SEQUENTIAL called'
#ifdef HJJTRA_DEBUG
      write(LUPRI,*) 'LWORK:',LWORK
      write(LUPRI,*) 'Test output of UCMO:'
      call OUTPUT(UCMO,1,NBAST,1,NORBT,NBAST,NORBT,-1,LUPRI)
      IPRFCK = 19
#endif

!
!...  Open direct access file for final MO integrals
!     (the ICDTRA array which links a C,D to a record has been saved on LUINTM)
!
      CALL GPINQ('MO2INT_HJJTRA','EXIST',FEXIST)
      IF (FEXIST) THEN
         IF (LUMINT .LE. 0) CALL
     &      GPOPEN(LUMINT,'MO2INT_HJJTRA','OLD','DIRECT',' ',LBUF,OLDDX)
         CALL GPCLOSE(LUMINT,'DELETE')
      END IF
      LUMINT = -1
      LBUF = IRAT*NNORBT
      CALL GPOPEN(LUMINT,'MO2INT_HJJTRA','NEW',
     &   'DIRECT','UNFORMATTED',LBUF,OLDDX)

      IFCTYP(1:MX_NDMAT) = 11 ! symmetric density matrix, only Coulomb

      NDMAT = 0
      ICD   = 0
      do IC = 1, NORBT
         do ID = 1, IC
            ICD = ICD + 1
         if (ICDTRA(ICD) .le. 0) cycle
            NDMAT = NDMAT + 1
            ISYMDM(NDMAT)  = MULD2H(ISMO(IC), ISMO(ID))
            ICD_REC(NDMAT) = ICDTRA(ICD)
            call HJJTRA_CD_DMAT(IC,ID,UCMO,DMAT(1,NDMAT))
#ifdef HJJTRA_DEBUG
            write(lupri,*) 'NDMAT, IC, ID, CD_DMAT:',
     &         NDMAT,IC,ID,ICD,ICDTRA(ICD)
            call output(DMAT(1,NDMAT),1,NBAST,1,NBAST,
     &         NBAST,NBAST,-1,LUPRI)
#endif
            if (NDMAT .eq. MX_NDMAT) then
#ifdef HJJTRA_DEBUG
         write (lupri,*) 'calling SIRFCK, NDMAT',NDMAT
#endif
               call DZERO(FMAT,NDMAT*N2BASX)
               call SIRFCK(FMAT,DMAT,NDMAT,ISYMDM,IFCTYP,DIRFCK,
     &                     WORK,LWORK)
               call HJJTRA_WR_LUMINT(UCMO,NDMAT,FMAT,DMAT,
     &                     ICD_REC,ISYMDM)
               NDMAT = 0
            end if
         end do
      end do

      if (NDMAT .gt. 0) then
#ifdef HJJTRA_DEBUG
         write (lupri,*) 'calling SIRFCK, NDMAT',NDMAT
#endif
         call DZERO(FMAT,NDMAT*N2BASX)
         call SIRFCK(FMAT,DMAT,NDMAT,ISYMDM,IFCTYP,DIRFCK,
     &               WORK,LWORK)
         call HJJTRA_WR_LUMINT(UCMO,NDMAT,FMAT,DMAT,
     &               ICD_REC,ISYMDM)
      end if
      
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
      subroutine HJJTRA_WR_LUMINT(UCMO,NDMAT,FMAT,DMAT,ICD_REC,ISYMDM)
!
!  Written by Hans Joergen Aa. Jensen October 2014
!
      implicit none
      integer NDMAT, ICD_REC(NDMAT), ISYMDM(NDMAT), IOS
      real*8  UCMO(NCMOT), DMAT(N2BASX,NDMAT), FMAT(N2BASX,NDMAT)
      real*8  H2CDPK(NNORBT)  ! automatic array
#include "priunit.h"

! Used from include files:
!  inforb.h : NCMOT,N2BASX, ...
!  inftap.h : LUMINT
#include "inforb.h"
#include "inftap.h"

      integer IDMAT, IC, ID

      do IDMAT = 1, NDMAT
         ! step 1: transform to MO basis, using DMAT(:,1) and FMAT(:,1) as work arrays
         call DGEMM('T','N',NORBT,NBAST,NBAST,1.0D0,
     &                 UCMO,NBAST,
     &                 FMAT(1,IDMAT),NBAST, 0.0D0,
     &                 DMAT,NORBT)
         call DGEMM('N','N',NORBT,NORBT,NBAST,1.0D0,
     &                 DMAT,NORBT,
     &                 UCMO,NBAST, 0.0D0,
     &                 FMAT,NORBT)
         ! step 2: pack
         call DGETSP(NORBT,FMAT,DMAT)
         call TRDPAK(DMAT,H2CDPK,NORB,IORB,NORBT, ISYMDM(IDMAT),1)
#ifdef HJJTRA_DEBUG
         write (LUPRI,*) 'IDMAT, rec, sym',
     &      IDMAT,ICD_REC(IDMAT),ISYMDM(IDMAT)
         write (LUPRI,*) 'Unpacked H2CD'
         call output(FMAT,1,NORBT,1,NORBT,NORBT,NORBT,-1,LUPRI)
         if (idmat == 1) then
         write (LUPRI,*) 'Packed H2CDPK'
         call outpak(H2CDPK,NORBT,-1,LUPRI)
         end if
#endif
         ! step 3: write to file
         write(LUMINT, rec=ICD_REC(IDMAT), iostat=IOS) H2CDPK
         if (IOS .ne. 0) then
            call QENTER('HJJTRA_WR_LUMINT')
            write(lupri,*) 'IDMAT, rec, IOS',IDMAT,ICD_REC(IDMAT),IOS
            write(lupri,*) 'LUMINT, NNORBT/X',LUMINT,NNORBT,NNORBX
            call QUIT('IOS .ne. 0')
         end if
      end do

      end

! ===================================================================
!     sirhjjtra.F section 2: read MO integrals
! ===================================================================
      subroutine HJJTRA_NXTH2M(IC,ID,H2CD,NEEDTP,WORK,IDIST)
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

      integer  :: IC, ID, NEEDTP(*), IDIST
      real*8   :: H2CD(NORBT,NORBT), WORK(*)

! Used from include files:
!  inforb.h : NORBT, ???
!  infind.h : ISMO(:)
!  inftap.h : LUINTM, LUMINT
#include "maxorb.h"
#include "maxash.h"
#include "inforb.h"
#include "infind.h"
#include "inftap.h"

      integer  :: MMORBX, ICD, ICD_1, ICDSYM
      !integer, allocatable :: ICDTRA(:)
      integer  :: ICDTRA(10000)
      integer  :: NCDTRA, ICD_SAVE
      save    NCDTRA, ICD_SAVE

#ifdef HJJTRA_DEBUG
      write(LUPRI,*) 'HJJTRA_NXTH2M called, IC, ID, IDIST =',
     &   IC, ID, IDIST
      call flshfo(LUPRI)
#endif

      if (IDIST .le. 0) then
         !allocate(ICDTRA(1:NNORBX))
         if (LUINTM .le. 0)
     &      call quit('LUINTM.le.0 in HJJTRA_NXTH2M')
         if (LUMINT .le. 0)
     &      call quit('LUMINT.le.0 in HJJTRA_NXTH2M')
         rewind (LUINTM)
         call MOLLAB('MO2INTHJ',LUINTM,LUPRI)
         read (LUINTM) MMORBX, NCDTRA, ICDTRA(1:MMORBX)
         if (NNORBX.ne.MMORBX)
     &      call QUIT('NNORBX.ne.MMORBX in HJJTRA_NXTH2M')
         ICD_SAVE = 0
      end if

      if (IDIST .lt. 0) then

         if (IC .gt. 0 .and. IC .le. NORBT .and.
     &       ID .gt. 0 .and. ID .le. NORBT) then
            write (LUPRI,*) 'HJJTRA_NXTH2M: invalid input IC, ID',IC,ID
            call quit('HJJTRA_NXTH2M: invalid input IC, ID')
         end if
         if (IC .ge. ID) then
            ICD = IC*(IC-1)/2 + ID
         else
            ICD = ID*(ID-1)/2 + IC
         end if
         IDIST = ICDTRA(ICD)
         if (ICD_SAVE .le. 0)
     &      call QUIT('IC,ID record no. .le. 0 in HJJTRA_NXTH2M')
         !deallocate(ICDTRA)

      else ! IDIST .ge. 0

         if (IDIST .eq. NCDTRA) then
            ! finished 
            IDIST = -1
            RETURN
         end if

         ICD_1 = ICD_SAVE + 1
         do ICD = ICD_1,NNORBX
            if (ICDTRA(ICD) .gt. 0) then
               IDIST    = ICDTRA(ICD)
               ICD_SAVE = ICD
               exit
            end if
         end do
         do IC = 1, NORBT
            ICD = IC*(IC+1)/2
            if (ICD .lt. ICD_SAVE) cycle
            ID = ICD_SAVE - ICD + IC
            exit
         end do

      end if

         ICDSYM = MULD2H( ISMO(IC), ISMO(ID) )
         read ( LUMINT, rec=IDIST ) WORK(1:NNORBT)
         call TRDPAK(WORK(1+NNORBT),WORK,NORB,IORB,NORBT,
     &      ICDSYM,-1)
         call DSPTGE(NORBT,WORK(1+NNORBT),H2CD)

#ifdef HJJTRA_DEBUG
         write (LUPRI,*) 'HJJTRA_NXTH2M, rec',IDIST
         write (LUPRI,*) 'Unpacked H2CD'
         call output(H2CD,1,NORBT,1,NORBT,NORBT,NORBT,-1,LUPRI)
         if (idist == 1) then
         write (LUPRI,*) 'Packed H2CDPK'
         call outpak(WORK,NORBT,-1,LUPRI)
         end if
#endif

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
