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


module pcm_interface
  
   use iso_c_binding

   implicit none 

   public pcm_interface_initialize
   public pcm_interface_finalize

   public collect_nctot
   public collect_atoms
   public pcm_energy_driver
   public pcm_oper_ao_driver
!   public pcm_lintra_driver

   private

!  if false the interface will refuse to be accessed
   logical :: is_initialized = .false.

   real(c_double), allocatable :: tess_cent(:, :)
   real(c_double)              :: pcm_energy
   integer(c_int)              :: nr_points

! A (maybe clumsy) way of passing LUPRI without using common blocks   
   integer                     :: global_print_unit

   contains 
      
      subroutine pcm_interface_initialize(print_unit)                              
     
      integer, intent(in) :: print_unit

      global_print_unit = print_unit 

      call init_pcm
      call print_pcm

      call get_cavity_size(nr_points)

      allocate(tess_cent(3, nr_points))
      tess_cent = 0.0d0
      call get_tess_centers(tess_cent)

      pcm_energy = 0.0d0
              
      is_initialized = .true.
                                                                 
      end subroutine
                                                                    
      subroutine pcm_interface_finalize()

      if (.not. is_initialized) then
         print *, 'Error: pcm_interface was never initialized.'
         stop 1
      end if
! Free the memory taken from the free store both in Fortran and in C++
      deallocate(tess_cent)

      call tear_down_pcm
                                                                    
      is_initialized = .false.
                                                                    
      end subroutine
                                                                    
      subroutine check_if_interface_is_initialized()

      if (.not. is_initialized) then
         print *, 'Error: pcm_interface is not initialized.'
         stop 1
      end if

      end subroutine

      subroutine collect_nctot(nr_nuclei) bind(c, name='collect_nctot_')

#include "mxcent.h"
#include "nuclei.h"

      integer(c_int), intent(out) :: nr_nuclei

      nr_nuclei = nucdep

      end subroutine collect_nctot
      
      subroutine collect_atoms(atomic_charges, atomic_centers) bind(c, name='collect_atoms_')
  
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

      subroutine pcm_energy_driver(density_matrix, pol_ene, work, lfree)

      use pcm_write, only: pcm_write_file_separate
      use pcmmod_cfg

      real(8)      :: density_matrix(*)
      real(8)      :: pol_ene
      real(8)      :: work(*)
      integer      :: lfree

! Make sure that the interface is initialized first
      call check_if_interface_is_initialized
! OK, now compute MEP and ASC
      call compute_mep_asc(density_matrix, work, lfree)

! pcm_energy is the polarization energy:
! U_pol = 0.5 * (U_NN + U_Ne + U_eN + U_ee)
      call comp_pol_ene_pcm(pol_ene)

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
      call get_electronic_mep(nr_points, tess_cent, asc, oper, work(kfree), lfree, .true.)
      deallocate(asc)

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
      integer                     :: kfree, i

      
      kfree   = 1
      
      allocate(mep(nr_points))
      mep = 0.0d0
      allocate(asc(nr_points))
      asc = 0.0d0

      if (.not.(pcmmod_separate)) then
         potName = 'TotMEP'//char(0) 
         chgName = 'TotASC'//char(0)
! Calculate the (total) Molecular Electrostatic Potential
         call get_mep(nr_points, tess_cent, mep, density_matrix, work(kfree), lfree)
! Set a cavity surface function with the MEP
         call set_surface_function(nr_points, mep, potName)
! Compute polarization charges and set the proper surface function
         call comp_chg_pcm(potName, chgName)
! Get polarization charges @tesserae centers
         call get_surface_function(nr_points, asc, chgName)

! Print some information
         if (pcmmod_print > 5) then
            do i = 1, nr_points
              write(global_print_unit, *) "MEP @point", i, mep(i)
              if (pcmmod_print > 10) then
                 write(global_print_unit, *) "ASC @point", i, asc(i)
              end if
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
         call get_nuclear_mep(nr_points, tess_cent, nuc_pot)
         call set_surface_function(nr_points, nuc_pot, potName1)
         call comp_chg_pcm(potName1, chgName1)
         call get_surface_function(nr_points, nuc_pol_chg, chgName1)

         potName2 = 'EleMEP'//char(0)
         chgName2 = 'EleASC'//char(0)
         call get_electronic_mep(nr_points, tess_cent, ele_pot, density_matrix, work(kfree), lfree, .false.)
         call set_surface_function(nr_points, ele_pot, potName2)
         call comp_chg_pcm(potName2, chgName2)
         call get_surface_function(nr_points, ele_pol_chg, chgName2)

! Print some information
        if (pcmmod_print > 5) then
           do i = 1, nr_points
             write(global_print_unit, *) "NMEP @point", i, nuc_pot(i)
             write(global_print_unit, *) "EMEP @point", i, ele_pot(i)
             if (pcmmod_print > 10) then
                write(global_print_unit, *) "NASC @point", i, nuc_pol_chg(i)
                write(global_print_unit, *) "EASC @point", i, ele_pol_chg(i)
             end if
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