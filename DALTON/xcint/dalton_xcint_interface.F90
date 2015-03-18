#ifdef ENABLE_XCINT
module dalton_xcint_interface

   use, intrinsic :: iso_c_binding
   use xcint_fortran_interface

   implicit none

   public interface_xc_init
   public interface_xc_finalize
   public is_ks_calculation
   public openrsp_set_functional

   private

   ! if false the interface will refuse to be accessed
   logical :: is_initialized = .false.

   ! true if this is a KS calculation
   logical :: this_is_ks_calculation = .false.

   ! FIXME
   real(8) :: openrsp_cfg_radint = 1.0d-12
   integer :: openrsp_cfg_angmin = 86
   integer :: openrsp_cfg_angint = 302
   logical :: openrsp_cfg_use_xcint_grid = .false.

contains

   subroutine check_if_initialized()
      if (.not. is_initialized) then
         print *, 'ERROR: you try to access interface_xc'
         print *, '       but this interface is not initialized'
         stop 1
      end if
   end subroutine


   subroutine interface_xc_finalize()
      is_initialized = .false.
      this_is_ks_calculation = .false.
   end subroutine


   logical function is_ks_calculation()
      call check_if_initialized()
      is_ks_calculation = this_is_ks_calculation
   end function


   integer(c_int) function fortran_stdout_function(string) bind(c)

   ! input
      ! string to print (c_null_char-terminated)
      character(kind=c_char, len=1), intent(in) :: string(*)

   ! local
      integer(c_int) :: i
      integer(c_int) :: string_len
      character(80) :: line

   ! returns
      ! 0 upon success

#include "priunit.h"

      i = 1
      do while (.true.)
         if (string(i) == c_null_char) then
            string_len = i - 1
            exit
         end if
         i = i + 1
      end do
      do i = 1, 80
         line(i:i) = ' '
      end do
      if (string_len < 81) then
         ! at the end we remove newline
         do i = 1, string_len-1
            line(i:i) = string(i)
         end do
      else
         do i = 1, 80
            line(i:i) = string(i)
         end do
      end if
      write(lupri, '(a80)') line
      fortran_stdout_function = 0

   end function


   subroutine interface_xc_init()

   ! local
      real(c_double), allocatable :: primitive_exp(:)
      real(c_double), allocatable :: contraction_coef(:)
      real(c_double), allocatable :: center_xyz(:)
      integer(c_int), allocatable :: shell_num_primitives(:)
      integer(c_int), allocatable :: l_quantum_num(:)
      integer(c_int), allocatable :: shell_center(:)
      integer(c_int), allocatable :: center_element(:)
      integer(c_int)              :: basis_type
      integer(c_int)              :: num_shells
      integer(c_int)              :: num_centers
      integer(c_int)              :: ierr
      integer                     :: i
      integer                     :: j
      integer                     :: icount
      integer                     :: iprim
      integer                     :: ishell
      integer                     :: iround
      integer                     :: n
      integer                     :: ixyz
      integer                     :: icenter

#include "aovec.h"
#include "mxcent.h"
#include "nuclei.h"
#include "maxorb.h"
#include "shells.h"
#include "primit.h"
#ifdef VAR_MPI
#include "mpif.h"
#endif

      if (is_initialized) return

      num_shells  = kmax
      num_centers = nucind
      allocate(center_xyz(3*num_centers))
      n = 1
      do icenter = 1, num_centers
         do ixyz = 1, 3
            center_xyz(n) = cord(ixyz, icenter)
            n = n + 1
         end do
      end do

      allocate(center_element(num_centers))
      do i = 1, num_centers
         center_element(i) = nint(charge(i))
      end do

      allocate(shell_num_primitives(num_shells))
      n = 0
      do iround = 1, 2
         if (iround == 2) then
            allocate(primitive_exp(n))
            allocate(contraction_coef(n))
            n = 0
         end if
         do ishell = 1, num_shells
            i = jstrt(ishell) + 1
            j = jstrt(ishell) + nuco(ishell)
            icount = 0
            do iprim = i, j
               if (dabs(priccf(iprim, numcf(ishell))) > tiny(0.0d0)) then
                  icount = icount + 1
                  shell_num_primitives(ishell) = icount
                  n = n + 1
                  if (iround == 2) then
                     contraction_coef(n) = priccf(iprim, numcf(ishell))
                     primitive_exp(n)    = priexp(iprim)
                  end if
               end if
            end do
         end do
      end do

      allocate(l_quantum_num(num_shells))
      allocate(shell_center(num_shells))
      do ishell = 1, num_shells
         l_quantum_num(ishell) = nhkt(ishell) - 1
         shell_center(ishell)  = ncent(ishell)
      end do

      ! we default to spherical basis
      basis_type = XCINT_BASIS_SPHERICAL
      do ishell = 1, num_shells
         if ((l_quantum_num(ishell) > 1) .and. .not. sphr(ishell)) then
            ! basis is cartesian
            basis_type = XCINT_BASIS_CARTESIAN
            exit
         end if
      end do

#ifdef VAR_MPI
      call xcint_set_mpi_comm(MPI_COMM_WORLD)
#endif

      call xcint_set_basis(basis_type,           &
                           num_centers,          &
                           center_xyz,           &
                           center_element,       &
                           num_shells,           &
                           shell_center,         &
                           l_quantum_num,        &
                           shell_num_primitives, &
                           primitive_exp,        &
                           contraction_coef)

      if (openrsp_cfg_use_xcint_grid) then
         ! this generates XCint's internal grid
         ierr = xcint_generate_grid(openrsp_cfg_radint,   &
                                    openrsp_cfg_angmin,   &
                                    openrsp_cfg_angint,   &
                                    num_centers,          &
                                    center_xyz,           &
                                    center_element,       &
                                    num_shells,           &
                                    shell_center,         &
                                    l_quantum_num,        &
                                    shell_num_primitives, &
                                    primitive_exp)
         if (ierr /= 0) then
            print *, 'OpenRSP ERROR: problem in xcint_generate_grid'
            stop 1
         end if
      end if

      deallocate(primitive_exp)
      deallocate(contraction_coef)
      deallocate(center_xyz)
      deallocate(shell_num_primitives)
      deallocate(l_quantum_num)
      deallocate(shell_center)
      deallocate(center_element)

      is_initialized = .true.

   end subroutine


   subroutine openrsp_set_functional(line, hfx, mu, beta)

   ! input
      character(*), intent(in) :: line

   ! output
      real(c_double), intent(out) :: hfx
      real(c_double), intent(out) :: mu
      real(c_double), intent(out) :: beta

   ! local
      integer        :: ierr
      type(c_funptr) :: stdout_function

      stdout_function = c_funloc(fortran_stdout_function)
      call xcint_set_stdout_function(stdout_function)
      call xcint_set_stderr_function(stdout_function)

      call xcint_print_splash()

      ierr = xcint_set_functional(line//C_NULL_CHAR, hfx, mu, beta)
      if (ierr /= 0) then
         print *, 'OpenRSP ERROR: problem in xcint_set_functional'
         stop 1
      end if

      this_is_ks_calculation = .true.

   end subroutine

end module


logical function is_dalton_ks_calculation()
   use dalton_xcint_interface, only: is_ks_calculation
   is_dalton_ks_calculation = is_ks_calculation()
end function


subroutine xcint_wakeup_workers()
#ifdef VAR_MPI

   ! local
      integer :: num_proc
      integer :: iprint
      integer :: ierr

#include "mpif.h"
#include "iprtyp.h"
#include "maxorb.h"
#include "infpar.h"

      call MPI_Comm_size(MPI_COMM_WORLD, num_proc, ierr)
      if (num_proc > 1) then
         CALL MPIXBCAST(XCINT_MPI_WAKEUP_SIGNAL, 1,'INTEGER', MASTER)
         iprint = 0
         CALL MPIXBCAST(iprint, 1,'INTEGER', MASTER)
      end if
#endif
end subroutine

#else
subroutine empty()
end subroutine
#endif /* ENABLE_XCINT */
