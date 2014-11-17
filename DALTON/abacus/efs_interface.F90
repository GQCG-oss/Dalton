#include "iprtyp.h"

module efs_interface

   use iso_c_binding

   implicit none

   !MAXPRD: max primitives per block
   !MXCONT: max number of contractions
   !MAXCONT = MAXPRD = 35

   type, bind(c) :: gto
      integer(c_int) :: angular_momentum
      integer(c_int) :: nr_primitives
      integer(c_int) :: nr_contractions
      real(c_double) :: primitive_exp(35)
      real(c_double) :: contract_coef(35, 35)
   end type

   type, bind(c) :: point
      real(c_double) :: x
      real(c_double) :: y
      real(c_double) :: z
   end type

   ! this is not used to interface with C
   type :: atom_type
      integer     :: nr_centers
      integer     :: nr_gtos
      type(point), allocatable :: centers(:)
      type(gto)  , allocatable :: gtos(:)
   end type

   logical :: use_efs = .false.

   logical :: efs_interface_is_initialized = .false.

   integer :: const_add_atomtype   = efs_add_atomtype_work
   integer :: const_init           = efs_init_work
   integer :: const_generate_basis = efs_generate_basis_work
   integer :: const_init_2efock    = efs_init_2efock_work
   integer :: const_fock_update    = efs_fock_update_work

end module


subroutine EFSprint(length, string)
    integer, intent(in) :: length
    character(length), intent(in) ::string

#include "priunit.h"

    write(LUPRI,*) string
end subroutine


subroutine EFSerror(length, string)
    integer, intent(in) :: length
    character(length), intent(in) :: string

    call quit(string)
end subroutine


