module dirac_gen1int_interface

!  radovan: this file is only used in DIRAC
!           but i keep it here for interface refactoring

   use gen1int_api

   implicit none

   public get_1el_integrals

   private

contains

   subroutine get_1el_integrals(io_viewer)

      integer, intent(in) :: io_viewer

      character(20), parameter :: PROP_NAME               = INT_OVERLAP
      integer,       parameter :: ORDER_MOM               = 0
      integer,       parameter :: ORDER_GEO_TOTAL         = 0
      integer,       parameter :: MAX_NUM_CENT            = 0
      logical,       parameter :: write_integrals_to_file = .false.
      type(prop_comp_t)        :: prop_operator
      integer                  :: num_ao
      integer, parameter       :: nr_dens_mat = 1
      integer                  :: ierr

! uses \var(MXCENT)
#include "mxcent.h"
! number of atomic centers
#include "nuclei.h"
! uses \var(MXCORB)
#include "maxorb.h"
! uses \var(MXQN)
#include "maxaqn.h"
! uses \var(MAXREP)
#include "symmet.h"

      integer                   :: num_prop         !number of property integral matrices
      integer                   :: kind_prop        !kind of property integral matrices
      logical                   :: triangular       !integral matrices are triangularized or squared
      logical                   :: symmetric        !integral matrices are symmetric or anti-symmetric
      integer                   :: num_redunt_geo   !number of all redundant total geometric derivatives
      integer                   :: num_opt_derv     !number of operators including derivatives
      type(matrix), allocatable :: val_ints(:)      !integral matrices
      integer                   :: imat             !incremental recorder over matrices

      ! gets the number of atomic orbitals
      call Gen1IntAPIGetNumAO(num_ao=num_ao)

!FIXME: implements subroutines in gen1int_host.F90 to return the number of property integrals
      ! gets the number of property integrals and total geometric derivatives
      num_prop = 1
      kind_prop = SYMM_INT_MAT
      num_redunt_geo = (3*NUCDEP)**ORDER_GEO_TOTAL
      num_opt_derv = num_redunt_geo*num_prop

      ! number and addresses of expectation values
      ! size and addresses of referenced integrals

      select case (kind_prop)
         case(SYMM_INT_MAT, ANTI_INT_MAT)
            triangular = .true.
            symmetric  = (kind_prop == SYMM_INT_MAT)
         case default
            triangular = .false.
            symmetric  = .false.
      end select

      ! allocates integral matrices
      allocate(val_ints(num_opt_derv))

      do imat = 1, num_opt_derv
         call MatCreate(A=val_ints(imat),      &
                        num_row=num_ao,        &
                        info_mat=ierr,         &
                        triangular=triangular, &
                        symmetric=symmetric)
         if (ierr /= 0) then
            print *, 'failed to create integral matrices'
            stop 1
         end if
      end do

      ! computes the integrals (using \fn(gen1int_host_get_expt) for expectation values);
      ! the operator will be created in \fn(Gen1IntAPIPropCreate) of gen1int_api.F90 (please
      ! update the operator there)
      call gen1int_host_get_int(NON_LAO, trim(PROP_NAME),      &
                                ORDER_MOM,                     &  !multipole moments
                                0, 0, 0,                       &  !magnetic derivatives
                                0, 0, 0,                       &  !derivatives w.r.t. total RAM
                                0, 0,                          &  !partial geometric derivatives
                                MAX_NUM_CENT, ORDER_GEO_TOTAL, &  !total geometric derivatives
                                0, (/0/), REDUNDANT_GEO,       &  !total geometric derivatives
                                .false., .false., .false.,     &  !not implemented
                                num_opt_derv, val_ints,        &
                                write_integrals_to_file,       &
                                io_viewer, 0)

      ! frees spaces of matrices
      do imat = 1, num_opt_derv
         print *, 'raboof'
         call MatView(A=val_ints(imat), io_viewer=io_viewer)
         call MatDestroy(A=val_ints(imat))
      end do
      deallocate(val_ints)

   end subroutine

end module
