!
! DALTON-side interface routines for the Polarizable Continuum Model
! Luca Frediani, Roberto Di Remigio 2011-2013
!
! We divide the interface routines into:
!   1. cavity-related routines;
!   2. solver-related routines.
!

!
!                   CAVITY-RELATED ROUTINES
!
! We shall provide the following data to the module:
!   1. the number of nuclei in the molecule. NUCLEI variable;
!   2. the atomic number of each nucleus. CHARGES vector;
!   3. the coordinates of each nucleus. CENTERS matrix.
! PCMSolver will provide to cavity formation depending
! on the selected mode:
!   1. Mode = Implicit means that the spheres are
!      centered according to the CENTERS matrix with
!      radii taken from the internal library;
!   2. Mode = Atoms means that we specify sphere centers
!      and radii for certain atoms in the molecule. All the
!      other spheres are created as in the Implicit mode.
!


module pcm_scf

  use iso_c_binding, only: c_null_char
  use pcmsolver
  use pcm_config
  use pcm_integrals, only: get_nuclear_mep, get_electronic_mep, get_mep
  use pcm_write, only: pcm_write_file, pcm_write_file_separate

  implicit none

  public pcm_scf_driver
  public pcm_mo_fock

  private

  real(8)              :: pcm_energy
  ! A counter for the number of SCF iterations
  integer, save        :: scf_iteration_counter = 1

contains

  function pcm_scf_driver(density_matrix, fock_matrix, work, lwork) result(pol_ene)

    real(8),    intent(in) :: density_matrix(*)
    real(8),   intent(out) :: fock_matrix(*)
    real(8), intent(inout) :: work(*)
    integer, intent(inout) :: lwork
    real(8)                :: pol_ene

    character(7) :: mep_name, asc_name
    real(8), allocatable :: asc(:)
    integer              :: kfree, lfree

    kfree = 1
    lfree = lwork - kfree

    if(lfree .le. 0) then
      call quit('Not enough mem in pcm_scf_driver')
    end if

    mep_name = 'TotMEP'//c_null_char
    asc_name = 'TotASC'//c_null_char
    call compute_mep_asc(density_matrix, work, lfree)
    pol_ene = pcmsolver_compute_polarization_energy(context_, mep_name, asc_name)
    pcm_energy = pol_ene
    allocate(asc(pcm_global%nr_points))
    asc = 0.0d0
    call pcmsolver_get_surface_function(context_, pcm_global%nr_points, asc, asc_name)
    call get_electronic_mep(pcm_global%nr_points, pcm_global%nr_points_irr, pcm_global%tess_cent, &
      asc, fock_matrix, work(kfree), lfree, .true.)
    deallocate(asc)
    scf_iteration_counter = scf_iteration_counter + 1

  end function pcm_scf_driver

  subroutine pcm_mo_fock(fock_matrix, mo_coeff, work, lwork)

#include "inforb.h"

    real(8),   intent(out) :: fock_matrix(*)
    real(8),    intent(in) :: mo_coeff(*)
    real(8), intent(inout) :: work(*)
    integer, intent(inout) :: lwork

    real(8), allocatable :: asc(:)
    real(8), allocatable :: ao_oper(:), ao_oper_packed(:)
    integer              :: kfree, lfree
    integer              :: isym, jcmo, j1ao, j1mo
    character(7)         :: asc_name

    kfree = 1
    lfree = lwork - kfree

    if(lfree .le. 0) then
      call quit('Not enough mem in pcm_mo_fock')
    end if

    asc_name = 'TotASC'//c_null_char
    allocate(asc(pcm_global%nr_points))
    asc = 0.0d0
    allocate(ao_oper(nnbasx))
    ao_oper = 0.0d0
    call pcmsolver_get_surface_function(context_, pcm_global%nr_points, asc, asc_name)
    call get_electronic_mep(pcm_global%nr_points, pcm_global%nr_points_irr, pcm_global%tess_cent, &
                            asc, ao_oper, work(kfree), lfree, .true.)
    deallocate(asc)

    ! Now we have the AO basis operator in ao_oper. Transform to MO basis
    allocate(ao_oper_packed(nnbasx))
    ao_oper_packed = 0.0d0
    call pksym1(ao_oper, ao_oper_packed, nbas, nsym, +1)

    deallocate(ao_oper)

    do isym = 1, nsym
      jcmo = icmo(isym) + 1
      j1ao = 1 + iibas(isym)
      j1mo = 1 + iiorb(isym)
      call uthu(ao_oper_packed(j1ao),fock_matrix(j1mo),mo_coeff(jcmo),work(kfree),nbas(isym),norb(isym))
    end do

    deallocate(ao_oper_packed)

  end subroutine pcm_mo_fock

  real(8) function get_pcm_energy()

    get_pcm_energy = pcm_energy

  end function get_pcm_energy

  subroutine compute_mep_asc(density_matrix, work, lfree)
    !
    ! Calculate the molecular electrostatic potential and
    ! the apparent surface charge at the cavity points.
    !
    ! The user can control via the DALTON input the following:
    !    * switch between separate and total evaluation of the
    !      nuclear and electronic parts;
    !    * switch between point-by-point and vectorized
    !      charge attraction integrals evaluation subroutines.
    !
    real(8), intent(in)    :: density_matrix(*)
    real(8), intent(inout) :: work(*)
    integer                :: lfree

    ! Local variables
    real(8), allocatable :: mep(:)
    real(8), allocatable :: asc(:)
    real(8), allocatable :: nuc_pot(:), nuc_pol_chg(:)
    real(8), allocatable :: ele_pot(:), ele_pol_chg(:)
    character(7)         :: potName, chgName
    character(7)         :: potName1, chgName1, potName2, chgName2
    integer              :: kfree, i, irrep

    kfree   = 1

    allocate(mep(pcm_global%nr_points))
    mep = 0.0d0
    allocate(asc(pcm_global%nr_points))
    asc = 0.0d0
    ! The totally symmetric irrep
    irrep = 0

    if (.not.(pcm_cfg%separate)) then
      potName = 'TotMEP'//c_null_char
      chgName = 'TotASC'//c_null_char
      ! Calculate the (total) Molecular Electrostatic Potential
      call get_mep(pcm_global%nr_points, pcm_global%nr_points_irr, pcm_global%tess_cent, mep, density_matrix, work(kfree), lfree)
      ! Set a cavity surface function with the MEP
      call pcmsolver_set_surface_function(context_, pcm_global%nr_points, mep, potName)
      ! Compute polarization charges and set the proper surface function
      call pcmsolver_compute_asc(context_, potName, chgName, irrep)
      ! Get polarization charges @tesserae centers
      call pcmsolver_get_surface_function(context_, pcm_global%nr_points, asc, chgName)

      ! Print some information
      if (pcm_cfg%print_level > 5) then
        write(pcm_global%print_unit, '(20X, A, 6X, I6)') "MEP and ASC at iteration", scf_iteration_counter
        write(pcm_global%print_unit, '(A, T27, A, T62, A)') "Finite element #", "Total MEP", "Total ASC"
        do i = 1, pcm_global%nr_points
          write(pcm_global%print_unit, '(I6, 2(20X, F15.12))') i, mep(i), asc(i)
        end do
      end if

      ! Write to file MEP and ASC
      call pcm_write_file(pcm_global%nr_points, mep, asc)
    else
      ! Allocation
      allocate(nuc_pot(pcm_global%nr_points))
      nuc_pot = 0.0d0
      allocate(nuc_pol_chg(pcm_global%nr_points))
      nuc_pol_chg = 0.0d0
      allocate(ele_pot(pcm_global%nr_points))
      ele_pot = 0.0d0
      allocate(ele_pol_chg(pcm_global%nr_points))
      ele_pol_chg = 0.0d0

      potName1 = 'NucMEP'//c_null_char
      chgName1 = 'NucASC'//c_null_char
      call get_nuclear_mep(pcm_global%nr_points, pcm_global%nr_points_irr, pcm_global%tess_cent, nuc_pot)
      call pcmsolver_set_surface_function(context_, pcm_global%nr_points, nuc_pot, potName1)
      call pcmsolver_compute_asc(context_, potName1, chgName1, irrep)
      call pcmsolver_get_surface_function(context_, pcm_global%nr_points, nuc_pol_chg, chgName1)

      potName2 = 'EleMEP'//c_null_char
      chgName2 = 'EleASC'//c_null_char
      call get_electronic_mep(pcm_global%nr_points, pcm_global%nr_points_irr, pcm_global%tess_cent, ele_pot, &
        density_matrix, work(kfree), lfree, .false.)
      call pcmsolver_set_surface_function(context_, pcm_global%nr_points, ele_pot, potName2)
      call pcmsolver_compute_asc(context_, potName2, chgName2, irrep)
      call pcmsolver_get_surface_function(context_, pcm_global%nr_points, ele_pol_chg, chgName2)

      ! Print some information
      if (pcm_cfg%print_level > 5) then
        write(pcm_global%print_unit, '(60X, A, 6X, I6)') "MEP and ASC at iteration", scf_iteration_counter
        write(pcm_global%print_unit, '(A, T27, A, T62, A, T97, A, T132, A)') "Finite element #", &
          "Nuclear MEP", "Nuclear ASC", "Electronic MEP", "Electronic ASC"
        do i = 1, pcm_global%nr_points
          write(pcm_global%print_unit, '(I6, 4(20X, F15.12))') i, nuc_pot(i), nuc_pol_chg(i), ele_pot(i), ele_pol_chg(i)
        end do
      end if

      ! Obtain vector of total MEP
      potName  = 'TotMEP'//c_null_char
      mep(:) = nuc_pot(:) + ele_pot(:)
      call pcmsolver_set_surface_function(context_, pcm_global%nr_points, mep, potName)

      ! Obtain vector of total polarization charges
      chgName  = 'TotASC'//c_null_char
      asc(:) = nuc_pol_chg(:) + ele_pol_chg(:)
      call pcmsolver_set_surface_function(context_, pcm_global%nr_points, asc, chgName)

      ! Write to file MEP and ASC
      call pcm_write_file_separate(pcm_global%nr_points, nuc_pot, nuc_pol_chg, ele_pot, ele_pol_chg)

      deallocate(nuc_pot)
      deallocate(nuc_pol_chg)
      deallocate(ele_pot)
      deallocate(ele_pol_chg)
    end if
    deallocate(mep)
    deallocate(asc)

  end subroutine compute_mep_asc

end module
