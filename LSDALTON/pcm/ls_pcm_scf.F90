!> @file
!> LSDALTON-side interface routines for the Polarizable Continuum Model
module ls_pcm_scf
  
use iso_c_binding
use typedeftype, only: lssetting 

implicit none 

public ls_pcm_scf_initialize
public ls_pcm_scf_finalize

public ls_pcm_energy_driver
public ls_pcm_oper_ao_driver

public get_pcm_energy

private

! if false the interface will refuse to be accessed
logical :: is_initialized = .false.

real(c_double), allocatable :: tess_cent(:, :)
real(c_double)              :: pcm_energy
integer(c_int)              :: nr_points
! A (maybe clumsy) way of passing LUPRI and LUERR  
integer                     :: global_print_unit
integer                     :: global_error_unit
type(lssetting)             :: integral_settings
! A counter for the number of SCF iterations
integer, save               :: scf_iteration_counter = 1

contains 
      
!> \brief initializes PCM-SCF interface  
!> \author R. Di Remigio
!> \date 2014
!> \param print_unit the LSDALTON output print unit
!>
!> This subroutine initializes various global objects internal to this
!> module.
subroutine ls_pcm_scf_initialize(setting, print_unit, err_unit)                              
     
   use ls_pcm_config, only: pcmtype, pcm_config

   type(lssetting), intent(in) :: setting
   integer,         intent(in) :: print_unit
   integer,         intent(in) :: err_unit
   ! nr_points_irr contains the number of points that are irreducible
   ! by symmetry. We won't be using symmetry but we need to preserve
   ! consistency with the API!
   integer(c_int)      :: nr_points_irr
   integer :: host_provides_input = 0
   
   global_print_unit = print_unit 
   global_error_unit = err_unit
   integral_settings = setting
   
   if (pcm_config%host_provides_input) host_provides_input = 1 

   call set_up_pcm(host_provides_input)
   call print_pcm
   
   call get_cavity_size(nr_points, nr_points_irr)
   
   allocate(tess_cent(3, nr_points))
   tess_cent = 0.0d0
   call get_tesserae(tess_cent)
   
   pcm_energy = 0.0d0
           
   is_initialized = .true.
                                                           
end subroutine ls_pcm_scf_initialize
                                                              
!> \brief finalizes PCM-SCF interface  
!> \author R. Di Remigio
!> \date 2014
!>
!> This subroutine finalizes various global objects internal to this
!> module.
subroutine ls_pcm_scf_finalize()

   if (.not. is_initialized) then
      call lsquit('Error: ls_pcm_scf was never initialized.', -1)
   end if
   ! Free the memory taken from the free store both in Fortran and in C++
   deallocate(tess_cent)
   
   call tear_down_pcm
                                                                 
   is_initialized = .false.
                                                              
end subroutine ls_pcm_scf_finalize
                                                                    
!> \brief checks if the PCM-SCF interface has been initialized
!> \author R. Di Remigio
!> \date 2014
subroutine check_if_interface_is_initialized()

   if (.not. is_initialized) then
      call lsquit('Error: ls_pcm_scf was never initialized.', -1)
   end if

end subroutine check_if_interface_is_initialized

!> \brief driver for PCM energy contribution calculation  
!> \author R. Di Remigio
!> \date 2014
!> \param density_matrix the current iteration density matrix
!> \param pol_ene the polarization energy
subroutine ls_pcm_energy_driver(density_matrix, pol_ene)

   use matrix_module

   type(matrix),      intent(in) :: density_matrix
   real(c_double), intent(inout) :: pol_ene
   
   ! Make sure that the interface is initialized first
   call check_if_interface_is_initialized
   
   ! OK, now compute MEP and ASC
   call compute_mep_asc(density_matrix)
   
   ! pcm_energy is the polarization energy:
   ! U_pol = 0.5 * (U_NN + U_Ne + U_eN + U_ee)
   call compute_polarization_energy(pol_ene)
   
   ! Now make the value of the polarization energy known throughout the module
   pcm_energy = pol_ene

end subroutine ls_pcm_energy_driver
      
!> \brief driver for PCM Fock matrix contribution calculation  
!> \author R. Di Remigio
!> \date 2014
!> \param oper NxN matrix to allocate PCM Fock matrix contribution
!>
!> This subroutine retrieves the uncontracted potentials at cavity points 
!> \latexonly $V_{\mu\nu}^I$ \endlatexonly and contracts them with the
!> apparent surface charge specified by the charge_name variable.
subroutine ls_pcm_oper_ao_driver(oper)
   
   use integralinterfacemod, only: II_get_ep_ab               
   use matrix_module
   
   type(matrix) :: oper
   
   character(7) :: charge_name
   real(c_double), allocatable :: asc(:)
   
   call check_if_interface_is_initialized
   
   charge_name = 'TotASC'//c_null_char
   allocate(asc(nr_points))
   asc = 0.0d0
   call get_surface_function(nr_points, asc, charge_name)
   call II_get_ep_ab(global_print_unit, &
                     global_error_unit, &
                     integral_settings, &
                     oper,              &
                     nr_points,         &
                     tess_cent,         &
                     asc)

   deallocate(asc)
   
   scf_iteration_counter = scf_iteration_counter + 1

end subroutine ls_pcm_oper_ao_driver

!> \brief handle to the PCM polarization energy  
!> \author R. Di Remigio
!> \date 2014
!>
!> Retrieve the PCM polarization energy contribution. 
real(c_double) function get_pcm_energy()
   
   get_pcm_energy = pcm_energy

end function get_pcm_energy                

!> \brief driver for PCM Fock matrix contribution calculation  
!> \author R. Di Remigio
!> \date 2014
!> \param density_matrix the current iteration density matrix
!>
!> Calculate the molecular electrostatic potential and        
!> the apparent surface charge at the cavity points.
!>                                                            
!> The user can control via the LSDALTON input the following:
!>    * switch between separate and total evaluation of the
!>      nuclear and electronic parts;
subroutine compute_mep_asc(density_matrix)
   
   use ls_pcm_integrals, only: get_nuclear_mep, get_electronic_mep, get_mep
   use ls_pcm_config, only: pcm_config
   use ls_pcm_write, only: ls_pcm_write_file, ls_pcm_write_file_separate
   use matrix_module
   
   type(matrix),    intent(in) :: density_matrix
   
   ! Local variables
   real(c_double), allocatable :: mep(:)
   real(c_double), allocatable :: asc(:)
   real(c_double), allocatable :: nuc_pot(:), nuc_pol_chg(:)
   real(c_double), allocatable :: ele_pot(:), ele_pol_chg(:)
   character(7)                :: potName, chgName
   character(7)                :: potName1, chgName1, potName2, chgName2
   integer                     :: i, irrep
   
   allocate(mep(nr_points))
   mep = 0.0d0
   allocate(asc(nr_points))
   asc = 0.0d0
   ! The totally symmetric irrep
   irrep = 0
   
   SeparateMEPandASC: if (.not.(pcm_config%separate)) then
      potName = 'TotMEP'//char(0) 
      chgName = 'TotASC'//char(0)
      ! Calculate the (total) Molecular Electrostatic Potential
      call get_mep(nr_points, tess_cent, mep, density_matrix, &
                   integral_settings, global_print_unit, global_error_unit)
      ! Set a cavity surface function with the MEP
      call set_surface_function(nr_points, mep, potName)
      ! Compute polarization charges and set the proper surface function
      call compute_asc(potName, chgName, irrep)
      ! Get polarization charges @tesserae centers
      call get_surface_function(nr_points, asc, chgName)
   
      ! Print some information
      PrintoutTotal: if (pcm_config%print_level > 5) then
         write(global_print_unit, '(20X, A, 6X, I6)') "MEP and ASC at iteration", scf_iteration_counter
         write(global_print_unit, '(A, T27, A, T62, A)') "Finite element #", "Total MEP", "Total ASC"
         do i = 1, nr_points
           write(global_print_unit, '(I6, 2(20X, F15.12))') i, mep(i), asc(i)
         end do
      end if PrintoutTotal
   
      ! Write to file MEP and ASC
      call ls_pcm_write_file(nr_points, mep, asc)
   else SeparateMEPandASC
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
      call compute_asc(potName1, chgName1, irrep)
      call get_surface_function(nr_points, nuc_pol_chg, chgName1)
   
      potName2 = 'EleMEP'//char(0)
      chgName2 = 'EleASC'//char(0)
      call get_electronic_mep(nr_points, tess_cent, ele_pot, density_matrix, &
                              integral_settings, global_print_unit, global_error_unit)
      call set_surface_function(nr_points, ele_pot, potName2)
      call compute_asc(potName2, chgName2, irrep)
      call get_surface_function(nr_points, ele_pol_chg, chgName2)
   
      ! Print some information
      PrintoutSeparate: if (pcm_config%print_level > 5) then
         write(global_print_unit, '(60X, A, 6X, I6)') "MEP and ASC at iteration", scf_iteration_counter
         write(global_print_unit, '(A, T27, A, T62, A, T97, A, T132, A)') "Finite element #", &
         "Nuclear MEP", "Nuclear ASC", "Electronic MEP", "Electronic ASC"
         do i = 1, nr_points
           write(global_print_unit, '(I6, 4(20X, F15.12))') i, nuc_pot(i), nuc_pol_chg(i), ele_pot(i), ele_pol_chg(i)
         end do
      end if PrintoutSeparate
   
      ! Obtain vector of total MEP
      potName  = 'TotMEP'//char(0)
      mep(:) = nuc_pot(:) + ele_pot(:)
      call set_surface_function(nr_points, mep, potName)
   
      ! Obtain vector of total polarization charges 
      chgName  = 'TotASC'//char(0)
      asc(:) = nuc_pol_chg(:) + ele_pol_chg(:)
      call set_surface_function(nr_points, asc, chgName)
   
      ! Write to file MEP and ASC
      call ls_pcm_write_file_separate(nr_points, nuc_pot, nuc_pol_chg, ele_pot, ele_pol_chg)
    
      deallocate(nuc_pot)
      deallocate(nuc_pol_chg)
      deallocate(ele_pot)
      deallocate(ele_pol_chg)
   end if SeparateMEPandASC
   deallocate(mep)
   deallocate(asc)

end subroutine compute_mep_asc
                                                                    
end module
