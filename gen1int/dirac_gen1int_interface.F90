module dirac_gen1int_interface

!  radovan: this file is only used in DIRAC
!           but i keep it here for interface refactoring

   use dalton_shell

   implicit none

   public get_1el_integrals

   private

contains

   subroutine get_1el_integrals(io_viewer)

      integer, intent(in) :: io_viewer

      character(20), parameter :: PROP_NAME               = INT_OVERLAP
      logical,       parameter :: LONDON_AO               = .false.
      integer,       parameter :: ORDER_MOM               = 0
      integer,       parameter :: ORDER_GEO_TOTAL         = 0
      integer,       parameter :: MAX_NUM_CENT            = 0
      logical,       parameter :: write_integrals_to_file = .false.
      type(one_prop_t)         :: prop_operator
      integer                  :: num_ao
      integer, parameter       :: nr_dens_mat = 1
      integer                  :: ierr

! origins
#include "orgcom.h"
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
      integer                   :: dim_unique_geo   !dimension of all unique total geometric derivatives
      integer                   :: num_opt_derv     !number of operators including derivatives
      type(matrix), allocatable :: val_ints(:)      !integral matrices
      integer                   :: imat             !incremental recorder over matrices

      ! gets the number of atomic orbitals
      call DaltonShellGetNumAO(num_ao=num_ao)

      ! initializes the information of one-electron property integrals
      call OnePropCreate(prop_name=trim(PROP_NAME),       &
                         one_prop=prop_operator,          &
                         info_prop=ierr,                  &
                         num_prop=num_prop,               &
                         kind_prop=kind_prop,             &
                         coord_nuclei=CORD(:,1:NUCDEP),   &
                         charge_nuclei=-CHARGE(1:NUCDEP), &
                         dipole_origin=DIPORG,            &
                         order_mom=ORDER_MOM)
      if (ierr /= 0) then
         print *, 'failed to create '//trim(PROP_NAME)
         stop 1
      end if

      dim_unique_geo = (3*NUCDEP)**ORDER_GEO_TOTAL
      num_opt_derv = dim_unique_geo*num_prop
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

      ! computes the integrals and/or expectation values

      call DaltonShellIntegral(comp_bra=1,                       &
                               comp_ket=1,                       &
                               one_prop=prop_operator,           &
                               london_ao=LONDON_AO,              &
                               num_ints=num_opt_derv,            &
                               val_ints=val_ints,                &
                               wrt_ints=write_integrals_to_file, &
                               num_dens=nr_dens_mat,             &
                               io_viewer=io_viewer,              &
                               level_print=0)

      call DaltonShellIntegral(comp_bra=2,                       &
                               comp_ket=2,                       &
                               one_prop=prop_operator,           &
                               london_ao=LONDON_AO,              &
                               num_ints=num_opt_derv,            &
                               val_ints=val_ints,                &
                               wrt_ints=write_integrals_to_file, &
                               num_dens=nr_dens_mat,             &
                               io_viewer=io_viewer,              &
                               level_print=0)

      ! frees space taken by the information of one-electron property integrals
      call OnePropDestroy(one_prop=prop_operator)

      do imat = 1, num_opt_derv
         print *, 'raboof'
         call MatView(A=val_ints(imat), io_viewer=io_viewer)
         call MatDestroy(A=val_ints(imat))
      end do
      deallocate(val_ints)

   end subroutine

end module
