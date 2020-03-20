!  Copyright (C) 2018 Andre Severo Pereira Gomes, Christoph Jacob, Lucas Visscher and collaborators
!
!  This file is part of Embed, a program implementing the Frozen Density Embedding (FDE) framework
! 
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at http://mozilla.org/MPL/2.0/.
!

module fde_evaluators_dalton

! embed modules
   use fde_types
   use fde_cfg
   use fde_data
   use fde_io

!! dalton-specific modules
!   use numoper_integrals
!   use fde_dalton_matrices_integration
!   use electrostatic_potential
!   use density_eval
!   use ao_eval
!   use interface_ao
!   use interface_mo_specific

! these are the interfaces to the dirac-specific code
   public fde_calculate_elpot

! below we have the routines that actually calculate
   public fde_get_density
   interface fde_get_density
      module procedure fde_dalton_get_density_from_dfcoef
   end interface

   public fde_get_density_matrix

   public fde_qccode_data_interface
   interface fde_qccode_data_interface
      module procedure fde_dalton_data_interface
   end interface

! the matrix for the embedding potential, or its contributions to
! different fock matrices can be calculated in different ways; 
! we take this into account by defining an interface name to the
! code-specific routines
 
   public fde_calculate_emb_pot_mat
   interface fde_calculate_emb_pot_mat
      module procedure fde_dalton_embpot_via_integrator
      module procedure fde_dalton_embpot_via_oneelectron
   end interface


   public fde_calculate_emb_linrsp_mat
   interface fde_calculate_emb_linrsp_mat
      module procedure fde_dalton_linrsp_via_integrator
   end interface


   private 

   integer, save :: nba, nba_orb, nba_aux

   integer, parameter :: max_irreps = 8
   integer, save :: nr_irreps
   integer, save :: nr_bas_irrep(max_irreps)

   integer, save :: nr_type_so
   real(kind=8), save :: dft_hri 

   logical, save :: is_dalton_initialized = .false.

   contains

!   ----------------------------------------------------------------------------
      subroutine fde_calculate_elpot(level,grid,gf)
!   ----------------------------------------------------------------------------
         type(grid_function), intent(inout) :: gf
         type(fde_grid), intent(in)         :: grid
         character(len=4)                   :: level
         real(kind=8), allocatable          :: dmat(:)

! dmat is allocated and set in the routine below
         call fde_dalton_get_dmat_from_dummy(dmat,level)
! we calculate both the nuclear and electrostatic potentials in one call
! if needed, we can split this up later 
         call fde_dalton_calculate_elpot(dmat,grid%r,gf%elpot,gf%nucpot)

         deallocate(dmat)

      end subroutine fde_calculate_elpot

!   ----------------------------------------------------------------------------
      subroutine fde_dalton_calculate_elpot(dmat,r,vc,vnuc)
!   ----------------------------------------------------------------------------
         real(kind=8), intent(in)  :: dmat(:)
         real(kind=8), intent(in)  :: r(:,:)
         real(kind=8), intent(inout) :: vc(:)
         real(kind=8), intent(in) :: vnuc(:)
         integer :: i
         integer :: np
         integer :: irrep = 1
         logical :: nuc_part = .false.    
         logical :: ele_part = .true.    
         
         np = size(r,2)
!
! calculating the hartree potential first
! 
         write (*,*) 'Generating the Hartree potential over the grid'
         write (*,*) '... write me !'
!
! now the nuclear potential
!
         write (*,*) 'Generating the nuclear potential over the grid'
         write (*,*) '... write me !'
         call quit('fde error:  Hartree and nuclear potentials over the grid not implemented')
!
! assemble the final electrostatic potential
!
         vc   =  vnuc + vc 

      end subroutine fde_dalton_calculate_elpot
      
!   ----------------------------------------------------------------------------
      subroutine fde_dalton_get_dmat_from_dummy(dmat,level)
!   ----------------------------------------------------------------------------
      character(len=4) :: level
      real(kind=8), allocatable, intent(out) :: dmat(:)

         write (*,*) '... write me !'
         call quit('fde_dalton_get_dmat_from_dummy called, but not implented')

      end subroutine fde_dalton_get_dmat_from_dummy 


!   ----------------------------------------------------------------------------
      subroutine fde_dalton_get_density_from_dfcoef(level,grid,gf)
!   ----------------------------------------------------------------------------
      type(grid_function), intent(inout) :: gf
      type(fde_grid), intent(in)         :: grid
      character(len=4), intent(in)       :: level

      INTEGER :: i
      INTEGER :: NDIMENSION
      INTEGER :: NDER=0, DOLND=0

      REAL(KIND=8),ALLOCATABLE :: GAO(:), GAO1(:),GAO2(:),BUF(:),NCNT(:)
      REAL(KIND=8),ALLOCATABLE :: GAB1(:)

      real(kind=8), allocatable :: dmat(:)

      REAL(KIND=8)::dummy,TOLS,TOLOG
      REAL(KIND=8)::gexp=0.0d0,ncentc=1.0d0,factor=1.0d0

      NDIMENSION=nba_orb

!     DMAT does into account the occupation numbers for scf and mp2 ...
      call fde_dalton_get_dmat_from_dummy(dmat,level)
#if 0

hjaaj 19-Mar-2020:
   this code is not used (fde_dalton_get_dmat_from_dummy quits ...)
   and it causes compiler warnings (some compilers)
  
!     so we correct it for scf and mp2...
      dmat = 2.0d0*dmat

      allocate(gao (   ndimension))
      allocate(gao1( 3*ndimension))
      allocate(gao2( 6*ndimension))
      allocate(gab1( 3*ndimension))
      allocate(ncnt(   ndimension))
      allocate(buf ( 4*ndimension))

!     Initialize the matrices to zero.

      GAO  = 0.0d0
      GAO1 = 0.0d0
      GAO2 = 0.0d0
      BUF  = 0.0d0
      NCNT = 0.0d0
      GAB1 = 0.0d0

      WRITE(*,*) 'Calculating the density and its derivatives.'
      WRITE(*,*) '- not implemented, write me!'
      call quit('fde_dalton_get_density_from_dfcoef not implemented')

      do i = 1, size(grid%r,2)

#ifdef FIXME
        CALL GETSOS(GAO,GAO1,GAO2,GAB1,NCNT,  &
                grid%r(1,i), grid%r(2,i), grid%r(3,i), &
                    BUF,NDIMENSION,2,DOLND,0)
        CALL GETRHO(gf%n(i),0,GAO,DMAT,BUF)
        CALL GETGRHO(gf%gn(1:3,i),0,GAO,GAO1,DMAT,BUF)
        CALL GETLRHO(gf%hn(1,i),gf%hn(1,i),gf%hn(1,i), &
                     gf%hn(4,i),gf%hn(5,i),gf%hn(6,i), &
                     GAO,GAO1,GAO2,DMAT,BUF)
#endif

      END DO

      deallocate (buf)
      deallocate (gao)
      deallocate (gao1)
      deallocate (gao2)
      deallocate (gab1)
      deallocate (ncnt) 
      deallocate (dmat)
#endif

      end subroutine fde_dalton_get_density_from_dfcoef



!   ----------------------------------------------------------------------------
      subroutine fde_dalton_embpot_via_oneelectron(mat_dim,fmat)
!   ----------------------------------------------------------------------------
         integer, intent(in) :: mat_dim
         real(kind=8), target  :: fmat(*)  ! fdemat(nnbasx)

         integer        :: p, q, nbuffer, lambda, pq
         integer        :: i ,irrep_offt, jpkds_offt, ibas_offt
         integer        :: iprdfe
         real(kind=8)   :: gamma_pq, ep_pq_update 
         real(kind=8), allocatable :: gao(:), buffer(:)
         integer, allocatable :: ncnt(:)

!        write (*,*) 'dimension: ',mat_dim,'  bla: ',fmat(1:5)
         call setupsos(0,.false.,idum1,idum2)
         call fde_dalton_data_interface

         iprdfe       = 1
         ndimension   = nba_orb
         nbuffer      = nr_irreps*nr_irreps*ndimension*nr_type_so 

         allocate(gao(ndimension*nr_type_so))
         allocate(ncnt(ndimension))
         allocate(buffer(nbuffer))

         nder   = 0
         irepop = 0

         do i = 1, size(fde_grid_sv%r,2)  ! begin numerical integration
!
!       getsos_fde is a simplification of getsos for our case. if there was
!       no symmetry, getaos could have been used directly... 
!
            call fde_getsos(gao, ncnt, fde_grid_sv%r(:,i), buffer, nbuffer,    &
                            ndimension,.false.,.false.,dft_hri, iprdfe)

!         in order to treat the symmetry blocking with packed storage, 
!         i (aspg) came up with the following (remember that the matrix
!         is stored as an array with the upper/lower triangular of the
!         first irrep coming before the second and so on):
!         
!         a. "irrep offset" (irrep_offt) will tell us where in the array a 
!             given irrep begins 
!         b. "packed storage offset" (jpkds_offt)) will tell us the number 
!            of elements from the packed storage within a given irrep that
!            we've looped over so far 
!         c. "basis offset" (ibas_offt) will tell us the starting index of 
!            the basis set in the gao array for the different irreps 
!
             irrep_offt = 0 
             jpkds_offt = 0
             ibas_offt  = 0

             do lambda = 1, nr_irreps

                do p = 1, nr_bas_irrep(lambda)
                   do q = p, nr_bas_irrep(lambda)    ! loop over the unique elements

                      gamma_pq      = gao(ibas_offt+p)*gao(ibas_offt+q)
                      ep_pq_update  = fde_grid_sv%w(i)*gamma_pq*fde_static_vemb(i)

!              we do the indexing for packed array storage [ij=i+j(j-1)/2] (see
!              e.g the lapack documentation), but taking care of the offsetting
!               due to the symmetry blocking
                     pq = irrep_offt + p + q*(q-1)/2
                     jpkds_offt =  jpkds_offt + 1

                     fmat(pq)  = fmat(pq) + ep_pq_update 
                  enddo ! loop over q
               enddo ! loop over p

               irrep_offt = irrep_offt + jpkds_offt 
               ibas_offt  = ibas_offt  + nr_bas_irrep(lambda)
               jpkds_offt = 0

            enddo ! loop over symmetry representations
         enddo ! loop over points

!        call flush()

         deallocate(gao)
         deallocate(ncnt)
         deallocate(buffer)

      end subroutine fde_dalton_embpot_via_oneelectron

!   ----------------------------------------------------------------------------
      subroutine fde_dalton_embpot_via_integrator(mat_dim,dmat,fmat)
!   ----------------------------------------------------------------------------
         integer, intent(in) :: mat_dim
         real(kind=8), target, intent(in) :: dmat(*)
         real(kind=8), target  :: fmat(*)

         write(*,*) 'write me!'
         call quit('fde error: embedded potential via integrator not implemented')

      end subroutine fde_dalton_embpot_via_integrator

!   ----------------------------------------------------------------------------
      subroutine fde_dalton_linrsp_via_integrator(mat_dim,dmat0,ndmat,dmat,fmat)
!   ----------------------------------------------------------------------------
         integer, intent(in)  :: mat_dim
         integer, intent(in)  :: ndmat
         real(kind=8), target, intent(in) :: dmat0(*)
         real(kind=8), target :: dmat(*)
         real(kind=8), target :: fmat(*)

         write(*,*) 'write me!'
         call quit('fde error: linear response via integrator not implemented')

      end subroutine fde_dalton_linrsp_via_integrator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! routines to set/recover values for certain dirac-specific variable
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   ----------------------------------------------------------------------------
      subroutine fde_dalton_data_interface
!   ----------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "nuclei.h" 
#include "dftcom.h"
#include "ccom.h"
#include "inforb.h"
#include "dftinf.h" 
#include "symmet.h"


         nba_orb       = nbast
         nba           = nba_orb
         mat_dim_as_1d = nnbasx

         nr_irreps     = maxrep + 1
         nr_bas_irrep  = 0
         do i = 1, nr_irreps
            nr_bas_irrep(i)  = naos(i) 
         enddo

         nr_type_so   = ntypso
         dft_hri      = DFTHRI

      end subroutine fde_dalton_data_interface


      SUBROUTINE FDE_GETSOS(GSO,NCNT,COR,WORK,LWORK,NBAST,              &
                        DOLND,DOGGA,DFTHRI,IPRINT)
!
!     original implementation of getsos: T. Helgaker feb 01
!     adaptation to fde case: A. Gomes, apr 07, review nov 12
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
!
      LOGICAL DOLND, DOGGA
      DIMENSION GSO(NBAST*NTYPSO), WORK(LWORK), NCNT(NBAST),COR(3)
!
#include "dftinf.h"
#include "symmet.h"
!
      IF (MAXREP.EQ.0) THEN
         CALL DFTAOS(GSO(KSO0),GSO(KSO1),GSO(KSO2),GSO(KSOB),GSO(KSOB1), &
     &               NCNT,COR(1),COR(2),COR(3),NBAST,                    &
     &               DOLND,DOGGA,DFTHRI,IPRINT)
      ELSE
         KGAO = 1
         KLST = KGAO + NTYPSO*NBAST
         IF (KLST.GT.LWORK) CALL STOPIT('DFTAOS','LWORK',KLST,LWORK)
         CALL DFTAOS(WORK(KSO0),WORK(KSO1),WORK(KSO2),WORK(KSOB),        &
     &               WORK(KSOB1),NCNT,COR(1),COR(2),COR(3),              &
     &               NBAST,DOLND,DOGGA,DFTHRI,IPRINT)
         CALL DFTSOS(WORK(KSO0),GSO,NBAST,NTYPSO,NCNT,IPRINT)
      END IF

      END SUBROUTINE FDE_GETSOS

end module fde_evaluators_dalton
