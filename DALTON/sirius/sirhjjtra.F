#define HJJTRA_DEBUG 2
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

#if HJJTRA_DEBUG > 0
      write(LUPRI,*) 'SIR_INTOPEN_HJJTRA called'
#endif

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
!  inftra.h : NCDTRA,IPRTRA
!  inftap.h : LUINTM, LBINTM
#include "gnrinf.h"
#include "inforb.h"
#include "inftra.h"
#include "inftap.h"

      integer JTRSIR, NTRLVL, MX_NDMAT
      integer KFREE, LFREE, KFRSAV, KUCMO, KDMAT, KFMAT, I
      integer, allocatable :: ICDTRA(:), ITRTYP(:)
      character*8 TABLE(2),LAB123(3)
      data        TABLE /'MO2INTHJ','END INTM'/
      data        LAB123/'********','********','********'/

C --- start of executable code

      call QENTER('HJJTRA_CTL')

      allocate(ICDTRA(NNORBX))
      allocate(ITRTYP(NORBT))

      KFREE  = 1
      LFREE  = LWORK
      KFRSAV = KFREE
     
      ! call SIR_INTOPEN_HJJTRA()

!     translate sirius integral level code to HJJTRA and NEWTRA integral level code
      JTRSIR = abs(ITRSIR)
      call TRALVL(JTRSIR, NTRLVL) 
      IF (JTRSIR .GE. 1 .AND. JTRSIR .LE. 4) NTRLVL = NTRLVL + 1
      ! we need third order, because USEDRC not defined for HJJTRA
#if HJJTRA_DEBUG > 0
      write(LUPRI,*) 'DEBUG: HJJTRA_CTL called, debug level', HJJTRA_DEBUG
      write(LUPRI,*) 'DEBUG: ITRSIR, NTRLVL, LWORK:',ITRSIR,NTRLVL,LWORK
#endif

!     call TRASET(sirntra.F) in order to be able to call TRALIM(sirntra.F)
!     which sets :
!        ICDTRA: index array for C,D distributions (**|CD) to calculate
!        ITRTYP(1:NORBT) = number of integral indices in which this orbital
!                          enters (i.e. 0,1,2,3, or 4)
      call TRASET(-1,LWORK)
#if HJJTRA_DEBUG > 0
      IPRTRA = MAX(IPRTRA,HJJTRA_DEBUG)
      write(lupri,*) 'DEBUG: IPRTRA set to',IPRTRA
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

#ifdef VAR_MPI
      if (HJJTRA_TYPE .eq. 1 .AND. PARCAL) then

         CALL HJJTRA_DISTRIBUTED_MASTER(ICDTRA, ITRTYP,
     &      CMO, WORK, LWORK)

      else ! not PARCAL
#endif

         call MEMGET2('REAL','UCMO',  KUCMO,NBAST*NORBT,
     &      WORK,KFREE,LFREE)

         MX_NDMAT = LFREE/(4*N2BASX) ! use half of free memory for DMAT and FMAT
         MX_NDMAT = MIN(NCDTRA,MX_NDMAT)
         write(lupri,*) 'MX_NDMAT, NNORBX =',MX_NDMAT,NNORBX
         call MEMGET2('REAL','DMAT',  KDMAT,MX_NDMAT*N2BASX,
     &      WORK,KFREE,LFREE)
         call MEMGET2('REAL','FMAT',  KFMAT,MX_NDMAT*N2BASX,
     &      WORK,KFREE,LFREE)

         call UPKCMO(CMO,WORK(KUCMO))
         call HJJTRA_FCK_DISTRIBUTED(ICDTRA, ITRTYP, WORK(KUCMO),
     &      MX_NDMAT, WORK(KDMAT), WORK(KFMAT), WORK(KFREE), LFREE)

         call MEMREL('after HJJTRA_FCK_DISTRIBUTED',WORK,
     &      1,KFRSAV,KFREE,LFREE)

#ifdef VAR_MPI
      end if !  if (HJJTRA_TYPE .eq. 1 .AND. PARCAL) then ... else ...
#endif

      WRITE (LUINTM) LAB123,TABLE(2)
      CALL REWSPL(LUINTM)

      deallocate ( ICDTRA )
      deallocate ( ITRTYP )
      call QEXIT('HJJTRA_CTL')
      return
      end

      subroutine HJJTRA_FCK_DISTRIBUTED(
     &     ICDTRA, ITRTYP, UCMO, MX_NDMAT, DMAT, FMAT, WORK, LWORK)
!
!  Written by Hans Joergen Aa. Jensen October 2014
!
      ! implicit none
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
!  dftcom.h : HFXFAC
#include "maxorb.h"
#include "maxash.h"
#include "infinp.h"
#include "inforb.h"
#include "infind.h"
#include "inftap.h"
#include "dftcom.h"
#include "cbihrs.h"

!  infpri.h : IPRTRA -- TODO : infpri.h not OK for implicit none yet
#include "infpri.h"

#include "iratdef.h"

      integer ICD, IC, ID, NDMAT,  LBUF
      integer, allocatable :: ICD_REC(:), IFCTYP(:), ISYMDM(:)
      logical FEXIST, OLDDX, TIMING
      real*8  TCPU0,TCPU1,TCPU2,TWAL0,TWAL1,TWAL2
      real*8  TCPU_FCK,TCPU_WR,TWAL_FCK,TWAL_WR

      logical NOSSUP_SAVE
      real*8  HFXFAC_SAVE

      allocate ( ICD_REC(MX_NDMAT) )
      allocate ( IFCTYP(MX_NDMAT) )
      allocate ( ISYMDM(MX_NDMAT) )

      TIMING = IPRTRA .ge. 0 ! DEBUG
      IF (TIMING) THEN
         call GETTIM(TCPU0,TWAL0)
         TCPU_FCK = 0.0D0
         TCPU_WR  = 0.0D0
         TWAL_FCK = 0.0D0
         TWAL_WR  = 0.0D0
      END IF
#if HJJTRA_DEBUG > 0
      write(LUPRI,*) 'HJJTRA_FCK_DISTRIBUTED called'
      write(LUPRI,*) 'LWORK:',LWORK
#if HJJTRA_DEBUG > 5
      write(LUPRI,*) 'Test output of UCMO:'
      call OUTPUT(UCMO,1,NBAST,1,NORBT,NBAST,NORBT,-1,LUPRI)
      IPRFCK = 19
#endif
#endif

!     for RDSUPM calls
      HFXFAC_SAVE = HFXFAC
      HFXFAC      = 0.0D0
      NOSSUP_SAVE = NOSSUP
      NOSSUP      = .true.
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
#if HJJTRA_DEBUG > 5
            write(lupri,*) 'NDMAT, IC, ID, CD_DMAT:',
     &         NDMAT,IC,ID,ICD,ICDTRA(ICD)
            call output(DMAT(1,NDMAT),1,NBAST,1,NBAST,
     &         NBAST,NBAST,-1,LUPRI)
#endif
            if (NDMAT .eq. MX_NDMAT) then
               call DZERO(FMAT,NDMAT*N2BASX)
               if (TIMING) call GETTIM(TCPU1,TWAL1)
               if (LUSUPM .ne. -1) then
#if HJJTRA_DEBUG > 1
         write (lupri,*) 'calling RDSUPM, NDMAT',NDMAT
#endif
                  call RDSUPM(NDMAT,FMAT,DMAT,ISYMDM,WORK,LWORK)
               else
#if HJJTRA_DEBUG > 1
         write (lupri,*) 'calling SIRFCK, NDMAT',NDMAT
#endif
                  call SIRFCK(FMAT,DMAT,NDMAT,ISYMDM,IFCTYP,DIRFCK,
     &                     WORK,LWORK)
               end if
               if (TIMING) then
                  call GETTIM(TCPU2,TWAL2)
#if HJJTRA_DEBUG > 1
         write (lupri,*) 'CPU and wall time',TCPU2-TCPU1,TWAL2-TWAL1
#endif
                  TCPU_FCK = TCPU_FCK + TCPU2 - TCPU1
                  TWAL_FCK = TWAL_FCK + TWAL2 - TWAL1
               end if
               call HJJTRA_WR_LUMINT(UCMO,NDMAT,FMAT,DMAT,
     &                     ICD_REC,ISYMDM)
               if (TIMING) then
                  call GETTIM(TCPU1,TWAL1)
                  TCPU_WR  = TCPU_WR  + TCPU1 - TCPU2
                  TWAL_WR  = TWAL_WR  + TWAL1 - TWAL2
               end if
               NDMAT = 0
            end if
         end do
      end do

      if (NDMAT .gt. 0) then
#if HJJTRA_DEBUG > 0
         write (lupri,*) 'calling SIRFCK, NDMAT',NDMAT
#endif
         call DZERO(FMAT,NDMAT*N2BASX)
         if (TIMING) call GETTIM(TCPU1,TWAL1)
         if (LUSUPM .ne. -1) then
#if HJJTRA_DEBUG > 1
            write (lupri,*) 'calling RDSUPM, NDMAT',NDMAT
#endif
            call RDSUPM(NDMAT,FMAT,DMAT,ISYMDM,WORK,LWORK)
         else
#if HJJTRA_DEBUG > 1
            write (lupri,*) 'calling SIRFCK, NDMAT',NDMAT
#endif
            call SIRFCK(FMAT,DMAT,NDMAT,ISYMDM,IFCTYP,DIRFCK,
     &                  WORK,LWORK)
         end if
         if (TIMING) then
            call GETTIM(TCPU2,TWAL2)
            TCPU_FCK = TCPU_FCK + TCPU2 - TCPU1
            TWAL_FCK = TWAL_FCK + TWAL2 - TWAL1
         end if
         call HJJTRA_WR_LUMINT(UCMO,NDMAT,FMAT,DMAT,
     &               ICD_REC,ISYMDM)
         if (TIMING) then
            call GETTIM(TCPU1,TWAL1)
            TCPU_WR  = TCPU_WR  + TCPU1 - TCPU2
            TWAL_WR  = TWAL_WR  + TWAL1 - TWAL2
         end if
      end if
      
      if (TIMING) then
         call GETTIM(TCPU1,TWAL1)
         write(lupri,*) 'CPU and WALL time (s) AO Fock',
     &      TCPU_FCK,TWAL_FCK
         write(lupri,*) 'CPU and WALL time (s) AO->MO ',
     &      TCPU_WR,TWAL_WR
         TCPU2 = TCPU1 - TCPU0 - TCPU_FCK - TCPU_WR
         TWAL2 = TWAL1 - TWAL0 - TWAL_FCK - TWAL_WR
         write(lupri,*) 'CPU and WALL time (s) rest   ',
     &      TCPU2,TWAL2
      end if

      HFXFAC = HFXFAC_SAVE
      NOSSUP = NOSSUP_SAVE

      deallocate ( ICD_REC )
      deallocate ( IFCTYP )
      deallocate ( ISYMDM )

      return
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
!  inforb.h : NORBT,NBAST
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
#if HJJTRA_DEBUG > 5
         write (LUPRI,*) 'IDMAT, rec, sym',
     &      IDMAT,ICD_REC(IDMAT),ISYMDM(IDMAT)
         write (LUPRI,*) 'Unpacked H2CD'
         call OUTPUT(FMAT,1,NORBT,1,NORBT,NORBT,NORBT,-1,LUPRI)
         if (idmat == 1) then
         write (LUPRI,*) 'Packed H2CDPK'
         call OUTPKB(H2CDPK,NORB,NSYM,-1,LUPRI)
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
!     sirhjjtra.F section 2: HJJTRA type 1 routines (if VAR_MPI)
! ===================================================================
#ifdef VAR_MPI
      subroutine HJJTRA_DISTRIBUTED_MASTER(ICDTRA, ITRTYP, CMO,
     &   WORK, LWORK)
!
!  Written by Hans Joergen Aa. Jensen October 2014
!
      implicit none
      integer ICDTRA(NNORBX), ITRTYP(NORBT), LWORK
      real*8  CMO(NCMOT), WORK(LWORK)
#include "priunit.h"

! Used from include files:
!  inforb.h : N2BASX, NORBT,NNORBX, ...
!  infpar.h, mpif.h : for MPI
!  iprtyp.h : defined parallel calculation types  
#include "maxorb.h"
#include "inforb.h"
#include "infpar.h"
#include "mpif.h"
#include "iprtyp.h"

      integer IPRTYP, IPRINT

!     start nodes on this task
!     1) make sure basis set information etc. is transferred
      IPRTYP = CALL_HJJTRA_DISTRIBUTED
      write(lupri,*) 'IPRTYP=',IPRTYP
      call MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
      IPRINT = 29
      call MPIXBCAST(IPRINT,1,'INTEGER',MASTER)

!     transfer information to nodes; we assume basis set information
!     already has been transferred in a Fock matrix call

      write (*,*) 'HJJTRA_DISTRIBUTED_MASTER has started slaves'
      write (lupri,*) 'HJJTRA_DISTRIBUTED_MASTER has started slaves'
      call flush(lupri)
      call sleep(5)
      call quit('more programming to do ...')
      
      return
      end
! ===================================================================
      subroutine HJJTRA_DISTRIBUTED_NODE(WORK,LWORK,IPRINT)
!
!  Written by Hans Joergen Aa. Jensen October 2014
!
      implicit none
      integer LWORK, IPRINT
      real*8  WORK(LWORK)
#include "priunit.h"

! Used from include files:
!  infpar.h, mpif.h : for MPI
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"

      write (lupri,*)
     &'HJJTRA_DISTRIBUTED_NODE has started, MYNYM =',MYNUM
      call flush(lupri)
      call sleep(2)
      call quit('more programming to do ...')

      return
      end
#endif   /* VAR_MPI */
! ===================================================================
!     sirhjjtra.F section 3: read MO integrals
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

      integer  :: IC, ID, NEEDTP(-4:6), IDIST
      real*8   :: H2CD(NORBT,NORBT), WORK(*)

! Used from include files:
!  inforb.h : NORBT, ???
!  infind.h : ISMO(:)
!  inftra.h : NCDTRA
!  inftap.h : LUINTM, LUMINT
#include "maxorb.h"
#include "maxash.h"
#include "inforb.h"
#include "infind.h"
#include "inftra.h"
#include "inftap.h"
#include "orbtypdef.h"

      integer       :: MMORBX, ICD, ICD_1, ICDSYM
      integer       :: ITYPC, ITYPD, ITYPCD
      integer, parameter :: LEN_ICDTRA = 160000
      integer, save :: ICDTRA(LEN_ICDTRA), ICD_SAVE

#if HJJTRA_DEBUG > 3
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
         if (MMORBX .gt. LEN_ICDTRA) then
            write (LUPRI,*) 'NNORBX>LEN_ICDTRA',MMORBX,LEN_ICDTRA
            call QUIT('NNORBX > LEN_ICDTRA')
         end if
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
         if (IDIST .le. 0) then
            call QUIT('IC,ID record no. .le. 0 in HJJTRA_NXTH2M')
         end if
         !deallocate(ICDTRA)

      else ! IDIST .ge. 0

  100    continue
         if (IDIST .eq. NCDTRA) then
            ! finished 
            IDIST = -1
            !deallocate(ICDTRA)
            return
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
         ITYPC  = IOBTYP(IC)
         ITYPD  = IOBTYP(ID)
         ITYPCD = IDBTYP(ITYPC,ITYPD)
         if ( NEEDTP(ITYPCD) .eq. 0 ) go to 100

      end if

         ICDSYM = MULD2H( ISMO(IC), ISMO(ID) )
         read ( LUMINT, rec=IDIST ) WORK(1:NNORBT)
         call TRDPAK(WORK(1+NNORBT),WORK,NORB,IORB,NORBT,
     &      ICDSYM,-1)
         call DSPTGE(NORBT,WORK(1+NNORBT),H2CD)

#if HJJTRA_DEBUG > 9
         write (LUPRI,*) 'HJJTRA_NXTH2M, rec',IDIST
         write (LUPRI,*) 'Unpacked H2CD'
         call output(H2CD,1,NORBT,1,NORBT,NORBT,NORBT,-1,LUPRI)
         if (idist == 1) then
         write (LUPRI,*) 'Packed H2CDPK'
         call outpkb(WORK,NORB,NSYM,-1,LUPRI)
         end if
#endif

      return
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
