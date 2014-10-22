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
  
   use iso_c_binding

   implicit none 

   public pcm_scf_initialize
   public pcm_scf_finalize

!   public collect_nctot
!   public collect_atoms
   public pcm_energy_driver
   public pcm_oper_ao_driver
!   public pcm_lintra_driver

   private

!  if false the interface will refuse to be accessed
   logical :: is_initialized = .false.

   real(c_double), allocatable :: tess_cent(:, :)
   real(c_double)              :: pcm_energy
   integer(c_int)              :: nr_points
   integer(c_int)              :: nr_points_irr
! A (maybe clumsy) way of passing LUPRI without using common blocks   
   integer                     :: global_print_unit
! A counter for the number of SCF iterations
   integer, save               :: scf_iteration_counter = 1

   contains 
      
      subroutine pcm_scf_initialize(print_unit)                              
      
      use pcmmod_cfg, only: pcmmod_host_provides_input
     
      integer, intent(in) :: print_unit
      
      integer :: host_provides_input = 0

      global_print_unit = print_unit 

      if (pcmmod_host_provides_input) host_provides_input = 1 

      call set_up_pcm(host_provides_input)
      call print_pcm

      call get_cavity_size(nr_points, nr_points_irr)

      allocate(tess_cent(3, nr_points))
      tess_cent = 0.0d0
      call get_tesserae(tess_cent)

      pcm_energy = 0.0d0
              
      is_initialized = .true.
                                                                 
      end subroutine
                                                                    
      subroutine pcm_scf_finalize()

      if (.not. is_initialized) then
         print *, 'Error: pcm_scf was never initialized.'
         stop 1
      end if
! Free the memory taken from the free store both in Fortran and in C++
      deallocate(tess_cent)

      call tear_down_pcm
                                                                    
      is_initialized = .false.
                                                                    
      end subroutine
                                                                    
      subroutine check_if_interface_is_initialized()

      if (.not. is_initialized) then
         print *, 'Error: pcm_scf is not initialized.'
         stop 1
      end if

      end subroutine

      subroutine collect_nctot(nr_nuclei) bind(c, name='collect_nctot')

#include "mxcent.h"
#include "nuclei.h"

      integer(c_int), intent(out) :: nr_nuclei

      nr_nuclei = nucdep

      end subroutine collect_nctot
      
      subroutine collect_atoms(atomic_charges, atomic_centers) bind(c, name='collect_atoms')
  
      use pcm_utils, only: getacord

#include "mxcent.h"
#include "nuclei.h"
      
      real(c_double), intent(out) :: atomic_charges(*)
      real(c_double), intent(out) :: atomic_centers(3,*)
      
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

      subroutine set_point_group(nr_gen, gen1, gen2, gen3) bind(c, name='set_point_group')

#include "mxcent.h"                     
#include "maxorb.h"
#include "maxaqn.h"
#include "symmet.h"

      integer(c_int), intent(inout) :: nr_gen 
      integer(c_int), intent(inout) :: gen1, gen2, gen3

      nr_gen = pcm_igen(1) 
      gen1   = pcm_igen(2) 
      gen2   = pcm_igen(3)
      gen3   = pcm_igen(4)

      end subroutine set_point_group

      subroutine pcm_energy_driver(density_matrix, pol_ene, work, lfree)

      real(8)        :: density_matrix(*)
      real(c_double) :: pol_ene
      real(8)        :: work(*)
      integer        :: lfree

! Make sure that the interface is initialized first
      call check_if_interface_is_initialized

! OK, now compute MEP and ASC
      call compute_mep_asc(density_matrix, work, lfree)

! pcm_energy is the polarization energy:
! U_pol = 0.5 * (U_NN + U_Ne + U_eN + U_ee)
      call compute_polarization_energy(pol_ene)

! Now make the value of the polarization energy known throughout the module
      pcm_energy = pol_ene

      end subroutine pcm_energy_driver
      
      subroutine pcm_oper_ao_driver(oper, charge_name, work, lwork)
!
! Calculate exp values of potentials on tesserae
! Input: symmetry packed Density matrix in AO basis
!        cavity points
! Output: expectation values of electrostatic potential on tesserae
!
      use pcm_integrals, only: get_electronic_mep

      real(8), intent(out)        :: oper(*)
      real(8)                     :: work(*)
      character(*), intent(in)    :: charge_name
      integer                     :: lwork

      real(c_double), allocatable :: asc(:)
      integer                     :: kfree, lfree

      call check_if_interface_is_initialized

      kfree = 1
      lfree = lwork - kfree
                                                                                  
      if(lfree .le. 0) then 
              call quit('Not enough mem in pcm_oper_ao_driver')
      end if
       
! Here we need the ASC, we can either get it from C++ (as it 
! has been saved as a surface function) or from Fortran asking
! for the mep and asc to be recomputed.
! RDR, 220813 !WARNING! what happens in the second case if we use this
! subroutine as a driver also for the linear response part??
      allocate(asc(nr_points))
      asc = 0.0d0
      call get_surface_function(nr_points, asc, charge_name)
      call get_electronic_mep(nr_points, nr_points_irr, tess_cent, asc, oper, work(kfree), lfree, .true.)
      deallocate(asc)

      scf_iteration_counter = scf_iteration_counter + 1

      end subroutine pcm_oper_ao_driver

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
      use pcm_integrals, only: get_nuclear_mep, get_electronic_mep, get_mep
      use pcm_write, only: pcm_write_file, pcm_write_file_separate
      use pcmmod_cfg

      real(8), intent(in)    :: density_matrix(*)
      real(8), intent(inout) :: work(*)
      integer                :: lfree

! Local variables
      real(c_double), allocatable :: mep(:)
      real(c_double), allocatable :: asc(:)
      real(c_double), allocatable :: nuc_pot(:), nuc_pol_chg(:)
      real(c_double), allocatable :: ele_pot(:), ele_pol_chg(:)
      character(7)                :: potName, chgName
      character(7)                :: potName1, chgName1, potName2, chgName2
      integer                     :: kfree, i, irrep
      
      kfree   = 1
      
      allocate(mep(nr_points))
      mep = 0.0d0
      allocate(asc(nr_points))
      asc = 0.0d0
      ! The totally symmetric irrep
      irrep = 0

      if (.not.(pcmmod_separate)) then
         potName = 'TotMEP'//char(0) 
         chgName = 'TotASC'//char(0)
! Calculate the (total) Molecular Electrostatic Potential
         call get_mep(nr_points, nr_points_irr, tess_cent, mep, density_matrix, work(kfree), lfree)
! Set a cavity surface function with the MEP
         call set_surface_function(nr_points, mep, potName)
! Compute polarization charges and set the proper surface function
         call compute_asc(potName, chgName, irrep)
! Get polarization charges @tesserae centers
         call get_surface_function(nr_points, asc, chgName)

! Print some information
         if (pcmmod_print > 5) then
            write(global_print_unit, '(20X, A, 6X, I6)') "MEP and ASC at iteration", scf_iteration_counter
            write(global_print_unit, '(A, T27, A, T62, A)') "Finite element #", "Total MEP", "Total ASC"
            do i = 1, nr_points
              write(global_print_unit, '(I6, 2(20X, F15.12))') i, mep(i), asc(i)
            end do
         end if

! Write to file MEP and ASC
         call pcm_write_file(nr_points, mep, asc)
      else
! Allocation
         allocate(nuc_pot(nr_points))
         nuc_pot = 0.0d0
         allocate(nuc_pol_chg(nr_points))
         nuc_pol_chg = 0.0d0
         allocate(ele_pot(nr_points))
         ele_pot = 0.0d0
         allocate(ele_pol_chg(nr_points))
         ele_pol_chg = 0.0d0
      
         potName1 = 'NucMEP'//char(0)
         chgName1 = 'NucASC'//char(0)
         call get_nuclear_mep(nr_points, nr_points_irr, tess_cent, nuc_pot)
         call set_surface_function(nr_points, nuc_pot, potName1)
         call compute_asc(potName1, chgName1, irrep)
         call get_surface_function(nr_points, nuc_pol_chg, chgName1)

         potName2 = 'EleMEP'//char(0)
         chgName2 = 'EleASC'//char(0)
         call get_electronic_mep(nr_points, nr_points_irr, tess_cent, ele_pot, density_matrix, work(kfree), lfree, .false.)
         call set_surface_function(nr_points, ele_pot, potName2)
         call compute_asc(potName2, chgName2, irrep)
         call get_surface_function(nr_points, ele_pol_chg, chgName2)

! Print some information
        if (pcmmod_print > 5) then
           write(global_print_unit, '(60X, A, 6X, I6)') "MEP and ASC at iteration", scf_iteration_counter
           write(global_print_unit, '(A, T27, A, T62, A, T97, A, T132, A)') "Finite element #", &
           "Nuclear MEP", "Nuclear ASC", "Electronic MEP", "Electronic ASC"
           do i = 1, nr_points
             write(global_print_unit, '(I6, 4(20X, F15.12))') i, nuc_pot(i), nuc_pol_chg(i), ele_pot(i), ele_pol_chg(i)
           end do
        end if

! Obtain vector of total MEP
        potName  = 'TotMEP'//char(0)
        mep(:) = nuc_pot(:) + ele_pot(:)
        call set_surface_function(nr_points, mep, potName)

! Obtain vector of total polarization charges 
        chgName  = 'TotASC'//char(0)
        asc(:) = nuc_pol_chg(:) + ele_pol_chg(:)
        call set_surface_function(nr_points, asc, chgName)

! Write to file MEP and ASC
        call pcm_write_file_separate(nr_points, nuc_pot, nuc_pol_chg, ele_pot, ele_pol_chg)
 
        deallocate(nuc_pot)
        deallocate(nuc_pol_chg)
        deallocate(ele_pot)
        deallocate(ele_pol_chg)
      end if
      deallocate(mep)
      deallocate(asc)

      end subroutine compute_mep_asc
                                                                    
end module
