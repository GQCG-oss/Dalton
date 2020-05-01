#define SRDFT_PAR_DEBUG -1
!
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2016 (2015), see http://daltonprogram.org"
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
!  FILE: DALTON/srdft/srdft_par.F
!  (c) Copyright Hans Joergen Aa. Jensen, hjj@sdu.dk (2017)


C*****************************************************************************
#ifdef VAR_MPI
      SUBROUTINE SRDFT_PAR_MASTER(ND_SIM,DMAT,EDFT,
     &                 DOERG,DO_MOLGRAD,DOATR,TRIPLET,
     &                 WORK,LWORK,IPRINT_in)
C*****************************************************************************
C     srdft_par.F master and slave routines
C
C  Written by Hans Joergen Aa. Jensen 2017
C
C  Master driver for parallel integration over srdft functional.
C
C  Parameter list same as for subroutine SRDFT,
C  see SRDFT for an explanation of the parameters.
C
C*****************************************************************************
      implicit none
      logical :: DOERG,DO_MOLGRAD,DOATR,TRIPLET
      integer :: ND_SIM, LWORK, IPRINT_in
      real*8 :: DMAT(N2BASX,*),EDFT(3),
     &          WORK(LWORK)

#include "priunit.h"

!  mpif.h   : for MPI
#include "mpif.h"


! Used from include files:
!  gnrinf.h : CHIVAL
!  infinp.h : /LOGINP/ for DOHFSRDFT, SRHYBR, ?
!  infvar.h : JWOPSY
!  inforb.h : N2BASX, ?
!  infpar.h : MYNUM, MASTER
!  dftcom.h : /DFTCOM/ transferred to slaves
!  iprtyp.h : PAR_SRDFT_INTEGRATION

#include "maxorb.h"
#include "maxash.h"
#include "gnrinf.h"
#include "infinp.h"
#include "infvar.h"
#include "inforb.h"
#include "infpar.h"
#include "dftcom.h"
#include "iprtyp.h"

C     variables for SRCPBERI
      real*8   mu_0, mu_0_m6, mu_0_m7, mu_0_m4
      common /cb_srcpbe_RI/ mu_0, mu_0_m6, mu_0_m7, mu_0_m4

      integer IPRTYP, IPRINT, IPRINT_slaves, IWHO, NWHO
      integer KFREE, LFREE
      integer NDMAT, NEXCMAT

      integer(mpi_integer_kind) :: bytesize_mpi, ierr_mpi

      call qenter('SRDFT_PAR_MASTER')

      IPRINT = max(SRDFT_PAR_DEBUG,iprint_in)

      IF (ND_SIM .NE. 1) THEN ! ND_SIM > 1 not implemented yet
         WRITE(LUPRI,*) 'SRDFT_PAR_MASTER: illegal ND_SIM ',ND_SIM
         CALL QUIT('SRDFT_PAR_MASTER ERROR: illegal ND_SIM value')
      END IF

!     start nodes on the distributed 2-electron integral transformation task

      IPRTYP = PAR_SRDFT_INTEGRATION ! get code number from iprtyp.h
      call MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
      IPRINT_slaves = IPRINT ! TODO set print level in input?
      call MPIXBCAST(IPRINT_slaves,1,'INTEGER',MASTER)

!     Now nodes are in subroutine SRDFT_PAR_NODE.

!     1) Make sure basis set information is known to nodes
      KFREE = 1
      LFREE = LWORK
      CALL HER_sendreceive_molinfo(IPRINT_slaves,WORK,KFREE,LFREE)

!     2) Make sure general information in gnrinf.h is known to nodes
!        (need at least CHIVAL and ERFEXP(:))
      call getbytespan(GRADML, GNRINFlast, bytesize_mpi)
      call mpi_bcast(GRADML, bytesize_mpi, mpi_byte,
     &   MASTER, mpi_comm_world, ierr_mpi)

!     3) Make sure dft information in dftcom.h is known to nodes
      call getbytespan(HFXFAC, DFTCOMlast, bytesize_mpi)
      call mpi_bcast(HFXFAC, bytesize_mpi, mpi_byte,
     &   MASTER, mpi_comm_world, ierr_mpi)

!     4) Make sure /LOGINP/ in infinp.h is known to nodes
!        (need at least DOHFSRDFT and SRHYBR)
      call getbytespan(FLAG(1), LOGINPlast, bytesize_mpi)
      call mpi_bcast(FLAG(1), bytesize_mpi, mpi_byte,
     &   MASTER, mpi_comm_world, ierr_mpi)

!     5) Transfer input control variables
      CALL MPI_BCAST(DOERG,1,my_MPI_LOGICAL,
     &     MASTER,MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(DO_MOLGRAD,1,my_MPI_LOGICAL,
     &     MASTER,MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(DOATR,1,my_MPI_LOGICAL,
     &     MASTER,MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(TRIPLET,1,my_MPI_LOGICAL,
     &     MASTER,MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(JWOPSY,1,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(mu_0,1,MPI_DOUBLE_PRECISION,
     &     MASTER,MPI_COMM_WORLD,ierr_mpi)

!     6) Transfer density matrices
      IF (DFT_SPINDNS) THEN
         NEXCMAT = 2
         NDMAT   = 3
      ELSE
         NEXCMAT = 1
         NDMAT   = 1
      END IF
      IF (DOATR) NDMAT = NDMAT + 1
      CALL MPI_BCAST(NDMAT,1,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(DMAT,NDMAT*N2BASX,MPI_DOUBLE_PRECISION,
     &     MASTER,MPI_COMM_WORLD,ierr_mpi)

      IF (iprint .gt. 5) then
      write (*,*) 'SRDFT_PAR_MASTER has started slaves, lupri', lupri
      write (lupri,*) 'SRDFT_PAR_MASTER has started slaves'
      write (lupri,*) 'IPRTYP,IPRINT_slaves=',IPRTYP,IPRINT_slaves
      flush(lupri)
      end if

      call qexit('SRDFT_PAR_MASTER')
      return
      end subroutine SRDFT_PAR_MASTER


C*****************************************************************************
      subroutine SRDFT_PAR_NODE(WORK,LWORK,IPRINT_in)
C*****************************************************************************
C
C  Written by Hans Joergen Aa. Jensen 2017
C
C  Node (slave) driver for parallel integration over srdft functional.
C
C*****************************************************************************
      implicit none
      integer :: LWORK, IPRINT_in
      real*8 :: WORK(LWORK)

!  mpif.h   : for MPI
#include "mpif.h"

      logical :: DOERG,DO_MOLGRAD,DOATR,TRIPLET
      integer :: ND_SIM, NDMAT, NEXCMAT, NDMAT_in, IPRINT
      integer :: KFREE, LFREE, KDMAT, KEXCMAT
      real*8 :: EDFT(3)

      logical, save :: first_call
      DATA      first_call/.true./

#include "priunit.h"

! Used from include files:
!  gnrinf.h : CHIVAL
!  infinp.h : /LOGINP/ for DOHFSRDFT, SRHYBR, ?
!  infvar.h : JWOPSY
!  inforb.h : N2BASX, ?
!  infpar.h : MYNUM, MASTER
!  dftcom.h : /DFTCOM/ transferred to slaves

#include "maxorb.h"
#include "maxash.h"
#include "gnrinf.h"
#include "infinp.h"
#include "infvar.h"
#include "inforb.h"
#include "infpar.h"
#include "dftcom.h"

C     variables for SRCPBERI
      real*8   mu_0, mu_0_m6, mu_0_m7, mu_0_m4
      common /cb_srcpbe_RI/ mu_0, mu_0_m6, mu_0_m7, mu_0_m4

      integer(mpi_integer_kind) :: bytesize_mpi, ierr_mpi

      call qenter('SRDFT_PAR_NODE')
      IPRINT = max(SRDFT_PAR_DEBUG,iprint_in)
      if (iprint .ge. 2) then
         write(lupri,*)
     &      'Hello from SRDFT_PAR_NODE mynum, lwork, iprint, first_call'
     &      ,MYNUM,LWORK,IPRINT,first_call
         call qdump(lupri)
      end if

!     ==================
!     Retrieve information from master
!     ==================

!     1) Make sure basis set information is known to nodes
      KFREE = 1
      LFREE = LWORK
      CALL HER_sendreceive_molinfo(IPRINT,WORK,KFREE,LFREE)
      CALL SETORB ! define many inforb.h variables

!     2) Make sure general information in gnrinf.h is known to nodes
!        (need at least CHIVAL and ERFEXP(:))
      call getbytespan(GRADML, GNRINFlast, bytesize_mpi)
      call mpi_bcast(GRADML, bytesize_mpi, mpi_byte,
     &   MASTER, mpi_comm_world, ierr_mpi)

!     3) Make sure dft information in dftcom.h is known to nodes
      call getbytespan(HFXFAC, DFTCOMlast, bytesize_mpi)
      call mpi_bcast(HFXFAC, bytesize_mpi, mpi_byte,
     &   MASTER, mpi_comm_world, ierr_mpi)
      ! master does not know if the grid file
      ! DALTON.QUAD has been generated for this node or not
      if (first_call .or. NEWGEO .or. NEWBAS) then
         DFTGRID_DONE_OLD = .FALSE.
      else
         DFTGRID_DONE_OLD = .TRUE.
      end if

!     4) Make sure /LOGINP/ in infinp.h is known to nodes
!        (need at least DOHFSRDFT and SRHYBR)
      call getbytespan(FLAG(1), LOGINPlast, bytesize_mpi)
      call mpi_bcast(FLAG(1), bytesize_mpi, mpi_byte,
     &   MASTER, mpi_comm_world, ierr_mpi)

!     5) Transfer input control variables
      CALL MPI_BCAST(DOERG,1,my_MPI_LOGICAL,
     &     MASTER,MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(DO_MOLGRAD,1,my_MPI_LOGICAL,
     &     MASTER,MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(DOATR,1,my_MPI_LOGICAL,
     &     MASTER,MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(TRIPLET,1,my_MPI_LOGICAL,
     &     MASTER,MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(JWOPSY,1,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(mu_0,1,MPI_DOUBLE_PRECISION,
     &     MASTER,MPI_COMM_WORLD,ierr_mpi)

!     6) Receive density matrices from master
      IF (DFT_SPINDNS) THEN
         NEXCMAT = 2
         NDMAT   = 3
      ELSE
         NEXCMAT = 1
         NDMAT   = 1
      END IF
      IF (DOATR) NDMAT = NDMAT + 1
      CALL MPI_BCAST(NDMAT_in,1,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,ierr_mpi)
      IF (NDMAT .NE. NDMAT_in) THEN
         WRITE(LUERR,*) 'ERROR: MYNUM, NDMAT, NDMAT_in:',
     &   MYNUM, NDMAT,NDMAT_in
         WRITE(LUPRI,*) 'ERROR: NDMAT, NDMAT_in:',NDMAT,NDMAT_in
         CALL QUIT('NDMAT error')
      END IF

      KFREE = 1
      LFREE = LWORK
      call MEMGET2('REAL','DMAT', KDMAT,NDMAT*N2BASX,
     &     WORK,KFREE,LFREE)
      call MEMGET2('REAL','EXCMAT', KEXCMAT,NEXCMAT*N2BASX,
     &     WORK,KFREE,LFREE)

      CALL MPI_BCAST(WORK(KDMAT),NDMAT*N2BASX,MPI_DOUBLE_PRECISION,
     &     MASTER,MPI_COMM_WORLD,ierr_mpi)

!     7) A check
#if SRDFT_PAR_DEBUG > 1
      write (lupri,*)
     &'SRDFT_PAR_NODE has started, MYNYM =',MYNUM
      write (lupri,*) 'DOERG,DO_MOLGRAD,DOATR',DOERG,DO_MOLGRAD,DOATR
      write (lupri,*) 'NBAST, N2BASX, NDMAT, NEXCMAT =',
     &   NBAST, N2BASX, NDMAT, NEXCMAT
      write (lupri,*) 'KDMAT,KEXCMAT,KFREE,LFREE',
     &                 KDMAT,KEXCMAT,KFREE,LFREE
      flush(lupri)
#endif

      IF (N2BASX .LE. 0) THEN
         CALL QUIT(
     &   'ERROR: Basis set information not transferred to slaves.')
      END IF

!     ==================
!     Do work
!     and transfer results to master in SRDFT
!     ==================

      ND_SIM = 1 
      EDFT(1:3)  = 0.0D0
      CALL DZERO(WORK(KEXCMAT),NEXCMAT*N2BASX)
      CALL SRDFT(
     &      ND_SIM,WORK(KEXCMAT),WORK(KDMAT),EDFT,
     &      DOERG,DO_MOLGRAD,DOATR,TRIPLET,
     &      WORK(KFREE),LFREE,IPRINT)
!     SUBROUTINE SRDFT(ND_SIM,EXCMAT,DMAT,EDFT(1:3),DOERG,DO_MOLGRAD,DOATR,
!    &                 WORK,LWORK,IPRINT)

      if (iprint .ge. 5) then
         write(lupri,*) 'Bye from SRDFT_PAR_NODE'
         flush(lupri)
      end if
      first_call = .false.
      call qexit('SRDFT_PAR_NODE')
      return
      end subroutine SRDFT_PAR_NODE
#else
      subroutine dummy_srdft_par
      ! to quench loader warnings for no symbols
      end
#endif   /* VAR_MPI */

! -- end of DALTON/sirius/srdft_par.F --
