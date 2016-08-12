!> @file
!> configure PCM input DALTON-side
module pcm_config

use, intrinsic :: iso_c_binding
use pcmsolver
use pcm_write, only: init_host_writer

implicit none

private

! if false the interface will refuse to be accessed
logical :: is_initialized = .false.

public pcm_initialize
public pcm_finalize

type, public :: pcm_configuration
  ! Polarizable continuum model calculation is turned off by default
  logical :: do_pcm = .false.
  ! Use of separate charges and potentials is turned off by default
  logical :: separate = .false.
  ! Use of vectorized integration routines is turned off by default
  logical :: fast_integration = .false.
  ! Print level is set to 0 by default
  integer :: print_level = 0
  ! Host does not parse PCMSolver input by default
  logical :: host_provides_input = .false.
end type pcm_configuration

type, public :: pcm_data
  integer :: print_unit
  integer :: error_unit
  integer(c_size_t)  :: nr_points
  integer(c_size_t)  :: nr_points_irr
  real(c_double), allocatable :: tess_cent(:)
end type pcm_data

type(pcm_configuration), public :: pcm_cfg
type(pcm_data),          public :: pcm_global
type(c_ptr),             public :: context_

! *PCMSOL section
! cavity specification *PCMSOL section
character(len=8), public, save :: pcmmod_cavity_type = 'gepol  '//c_null_char
integer(4), public, save  :: pcmmod_patch_level = 2
real(8), public, save  :: pcmmod_coarsity = 0.5
real(8), public, save  :: pcmmod_cavity_area = 0.3
real(8), public, save  :: pcmmod_min_distance = 0.1
integer(4), public, save  :: pcmmod_der_order = 4
logical, public, save :: pcmmod_scaling = .true.
character(len=8), public, save :: pcmmod_radii_set = 'bondi  '//c_null_char
character(len=20), public, save :: pcmmod_restart_name = 'cavity.npz         '//c_null_char
real(8), public, save  :: pcmmod_min_radius = 100.0
! solver specification *PCMSOL section
character(len=7), public, save :: pcmmod_solver_type = 'iefpcm'//c_null_char
character(len=16), public, save :: pcmmod_solvent = '               '//c_null_char
character(len=11), public, save :: pcmmod_equation_type = 'secondkind'//c_null_char
real(8), public, save  :: pcmmod_correction = 0.0
real(8), public, save  :: pcmmod_probe_radius = 1.0
! green specification *PCMSOL section
character(len=7), public, save :: pcmmod_inside_type = 'vacuum'//c_null_char
character(len=22), public, save :: pcmmod_outside_type = 'uniformdielectric    '//c_null_char
real(8), public, save :: pcmmod_outside_epsilon = 1.0

contains

!> \brief initializes interface to PCM
!> \author R. Di Remigio
!> \date 2015
!> \param print_unit global print unit
!> \param err_unit global error unit
subroutine pcm_initialize(print_unit, err_unit)

  integer :: print_unit, err_unit

  type(pcm_data) :: pcm_tmp
  integer(c_int) :: nr_nuclei
  real(c_double), allocatable :: charges(:)
  real(c_double), allocatable :: centers(:)
  integer(c_int) :: symmetry_info(4)
  type(PCMInput) :: host_input

  if (.not. pcmsolver_is_compatible_library()) then
      call quit('Error: incompatible version of PCMSolver library.')
  end if

  nr_nuclei = collect_nctot()
  allocate(charges(nr_nuclei))
  charges = 0.0_c_double
  allocate(centers(3*nr_nuclei))
  centers = 0.0_c_double
  call collect_atoms(charges, centers)
  call collect_symmetry_info(symmetry_info)
  if (pcm_cfg%host_provides_input) then
      host_input = pcmsolver_input()
      context_ = pcmsolver_new(PCMSOLVER_READER_HOST, nr_nuclei, charges, centers, symmetry_info, host_input)
  else
      context_ = pcmsolver_new(PCMSOLVER_READER_OWN, nr_nuclei, charges, centers, symmetry_info, host_input)
  end if

  deallocate(charges)
  deallocate(centers)

  pcm_tmp%nr_points = pcmsolver_get_cavity_size(context_)
  pcm_tmp%nr_points_irr = pcmsolver_get_irreducible_cavity_size(context_)

  ! Time initialization of tess_cent
  allocate(pcm_tmp%tess_cent(3*pcm_tmp%nr_points))
  pcm_tmp%tess_cent = 0.0d0
  call pcmsolver_get_centers(context_, pcm_tmp%tess_cent)

  pcm_tmp%print_unit = print_unit
  pcm_tmp%error_unit = err_unit

  pcm_global = pcm_tmp

  is_initialized = .true.

  call pcmsolver_print(context_)

end subroutine pcm_initialize

!> \brief finalizes interface to PCM-SCF
!> \author R. Di Remigio
!> \date 2015
!>
!> This subroutine finalizes various global objects
subroutine pcm_finalize()

   if (.not. is_initialized) then
       call quit('Error: PCM was never initialized.')
   else
       ! Free the memory taken from the free store both in Fortran and in C++
       deallocate(pcm_global%tess_cent)
       !call pcmsolver_delete(context_)
   end if

   is_initialized = .false.

end subroutine pcm_finalize

!> \brief sets PCMSolver input parameters from LSDALTON input
!> \author R. Di Remigio
!> \date 2014
!> \param cavity struct holding cavity parameters
! Performs syntactic checks on PCMSolver input and fills the data structures
! holding input data
function pcmsolver_input() result(host_input)

  type(PCMInput) :: host_input

  character(kind=c_char, len=1) :: cavity_type(8)
  character(kind=c_char, len=1) :: radii_set(8)
  character(kind=c_char, len=1) :: restart_name(20)
  character(kind=c_char, len=1) :: solver_type(7)
  character(kind=c_char, len=1) :: solvent(16)
  character(kind=c_char, len=1) :: equation_type(11)
  character(kind=c_char, len=1) :: inside_type(7)
  character(kind=c_char, len=1) :: outside_type(22)

  call pcmsolver_f2c_string(pcmmod_cavity_type, cavity_type, 8_c_int)
  host_input%cavity_type  = cavity_type
  host_input%patch_level  = int(pcmmod_patch_level, kind=c_int)
  host_input%coarsity     = pcmmod_coarsity
  host_input%area         = pcmmod_cavity_area
  host_input%min_distance = pcmmod_min_distance
  host_input%der_order    = int(pcmmod_der_order, kind=c_int)
  host_input%scaling      = pcmmod_scaling
  call pcmsolver_f2c_string(pcmmod_radii_set, radii_set, 8_c_int)
  host_input%radii_set    = radii_set
  call pcmsolver_f2c_string(pcmmod_restart_name, restart_name, 20_c_int)
  host_input%restart_name = restart_name
  host_input%min_radius   = pcmmod_min_radius

  call pcmsolver_f2c_string(pcmmod_solver_type, solver_type, 7_c_int)
  host_input%solver_type   = solver_type
  call pcmsolver_f2c_string(pcmmod_solvent, solvent, 16_c_int)
  host_input%solvent       = solvent
  call pcmsolver_f2c_string(pcmmod_equation_type, equation_type, 11_c_int)
  host_input%equation_type = equation_type
  host_input%correction    = pcmmod_correction
  host_input%probe_radius  = pcmmod_probe_radius

  call pcmsolver_f2c_string(pcmmod_inside_type, inside_type, 7_c_int)
  host_input%inside_type     = inside_type
  host_input%outside_epsilon = pcmmod_outside_epsilon
  call pcmsolver_f2c_string(pcmmod_outside_type, outside_type, 7_c_int)
  host_input%outside_type    = outside_type

end function pcmsolver_input

function collect_nctot() result(nr_nuclei)

#include "mxcent.h"
#include "nuclei.h"

  integer(4) :: nr_nuclei

  nr_nuclei = nucdep

end function collect_nctot

subroutine collect_atoms(atomic_charges, atomic_centers)

#include "mxcent.h"
#include "nuclei.h"

  real(8), intent(out) :: atomic_charges(*)
  real(8), intent(out) :: atomic_centers(3,*)

  integer :: i, j, k

  ! Get coordinates
  call getacord(atomic_centers)
  ! Get charges
  i = 0
  do j = 1, nucind
    do k = 1, nucdeg(j)
      i = i + 1
      atomic_charges(i) = charge(j)
    enddo
  enddo

end subroutine collect_atoms

subroutine collect_symmetry_info(symmetry_info)

#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "symmet.h"

  integer(4), intent(inout) :: symmetry_info(4)

  symmetry_info = (/pcm_igen(1), pcm_igen(2), pcm_igen(3), pcm_igen(4)/)

end subroutine collect_symmetry_info

end module pcm_config

!> \brief print relevant setting of both LSDALTON and PCMSolver
!> \author R. Di Remigio
!> \date 2014
!> \param print_unit the printing unit to be used
!> \param pcm_input PCM input section
subroutine report_after_pcm_input(print_unit, pcm_cfg)

  use pcm_config, only: pcm_configuration
  use pcm_write, only: init_host_writer

  integer, optional, intent(in) :: print_unit
  type(pcm_configuration), intent(in) :: pcm_cfg

  if (present(print_unit)) then
    ! Initialize host writer
    call init_host_writer(print_unit)
    write(print_unit, *)
    write(print_unit, *) &
      ' ===== Polarizable Continuum Model calculation set-up ====='
    write(print_unit, *) &
      '* Polarizable Continuum Model using PCMSolver external module:'
    write(print_unit, *) &
      '  1: Converged potentials and charges at tesserae representative points written on file.'

    if (pcm_cfg%separate) then
      write(print_unit, *) &
        '  2: Separate potentials and apparent charges in nuclear and electronic.'
    else
      write(print_unit, *) &
        '  2: Use total potentials and apparent charges.'
    end if

    if (pcm_cfg%print_level > 5 .and. pcm_cfg%print_level < 10) then
      write(print_unit, *) &
        '  3: Print potentials at tesserae representative points.'
    else if (pcm_cfg%print_level > 10) then
      write(print_unit, *) &
        '  3: Print potentials and charges at tesserae representative points.'
    else
      write(print_unit, *) &
        '  3: Do not print potentials and charges.'
    end if
  else
    write(*, *)
    write(*, *) &
      ' ===== Polarizable Continuum Model calculation set-up ====='
    write(*, *) &
      '* Polarizable Continuum Model using PCMSolver external module:'
    write(*, *) &
      '  1: Converged potentials and charges at tesserae representative points written on file.'

    if (pcm_cfg%separate) then
      write(*, *) &
        '  2: Separate potentials and apparent charges in nuclear and electronic.'
    else
      write(*, *) &
        '  2: Use total potentials and apparent charges.'
    end if

    if (pcm_cfg%print_level > 5 .and. pcm_cfg%print_level < 10) then
      write(*, *) &
        '  3: Print potentials at tesserae representative points.'
    else if (pcm_cfg%print_level > 10) then
      write(*, *) &
        '  3: Print potentials and charges at tesserae representative points.'
    else
      write(*, *) &
        '  3: Do not print potentials and charges.'
    end if
  end if

end subroutine report_after_pcm_input

subroutine read_input_pcm(word, kw_section)

  use pcm_config, only: pcm_configuration, pcm_cfg
  use input_reader

  implicit none

  !    --------------------------------------------------------------------------
  character(kw_length), intent(in) :: word
  character(kw_length), intent(in) :: kw_section
  !    --------------------------------------------------------------------------

  call reset_available_kw_list()

  if (kw_matches(word, '.FASTIN')) then
    ! Use "new" charge-attraction integrals evaluation subroutines
    pcm_cfg%fast_integration = .true.
  end if

  if (kw_matches(word, '.SEPARA')) then
    ! Split potentials and polarization charges in nuclear and electronic
    pcm_cfg%separate = .true.
  end if

  if (kw_matches(word, '.PRINT ')) then
    call kw_read(word, pcm_cfg%print_level)
  end if

  call check_whether_kw_found(word, kw_section)

end subroutine

subroutine getacord(coora)
!*****************************************************************************
!
! getacord : Make list atomic coordinates
!
!  Written oct.2001 by Jesper Kielberg Pedersen
!  Copied here from DIRAC by Roberto Di Remigio, February 2012
!
!*****************************************************************************
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "nuclei.h"
#include "symmet.h"
#include "pgroup.h"
  real(8), intent(out) :: coora(3,*)

  integer              :: jatom, icent, mulcnt, isymop
  !
  !     Make the full matrix of cartesian coordinates from CORD(NUCIND)
  !
  jatom = 0
  do icent = 1, nucind
    mulcnt = istbnu(icent)
    if (mult(mulcnt) .eq. 1) then
      jatom = jatom + 1
      coora(1,jatom) = cord(1,icent)
      coora(2,jatom) = cord(2,icent)
      coora(3,jatom) = cord(3,icent)
    else
      do isymop = 0, maxopr
        if (iand(isymop,mulcnt) .eq. 0) then
          jatom = jatom + 1
          coora(1,jatom) = pt(iand(isymax(1,1),isymop))*cord(1,icent)
          coora(2,jatom) = pt(iand(isymax(2,1),isymop))*cord(2,icent)
          coora(3,jatom) = pt(iand(isymax(3,1),isymop))*cord(3,icent)
        end if
      enddo
    end if
  enddo

end subroutine

