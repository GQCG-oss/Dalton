!> @file
!> thin wrappers for integral evaluation routines needed fo PCM 
module ls_pcm_integrals

use typedeftype, only: lssetting

implicit none

public get_mep
public get_nuclear_mep
public get_electronic_mep

contains

!> \brief calculates molecular electrostatic potential (MEP) at cavity points 
!> \author R. Di Remigio
!> \date 2014
!> \param nr_points the number of cavity points
!> \param centers the cavity points
!> \param mep the molecular electrostatic potential vector
!> \param density a density matrix
!>
!> Retrieves both electronic and nuclear MEP given an electronic density and a list
!> of points.
subroutine get_mep(nr_points, centers, mep, density, setting, lupri, luerr)

   use matrix_module
   use matrix_operations

   integer,         intent(in)  :: nr_points             
   real(8),         intent(in)  :: centers(3, nr_points)
   real(8),         intent(out) :: mep(nr_points)
   type(matrix),    intent(in)  :: density
   type(lssetting), intent(in)  :: setting
   integer,         intent(in)  :: lupri, luerr
   
   call get_electronic_mep(nr_points, centers, mep, density, .false., setting, lupri, luerr)
   call get_nuclear_mep(nr_points, centers, mep)

end subroutine get_mep
             
!> \brief calculates nuclear molecular electrostatic potential (MEP) at cavity points 
!> \author R. Di Remigio
!> \date 2014
!> \param nr_points the number of cavity points
!> \param centers the cavity points
!> \param mep the nuclear molecular electrostatic potential vector
!>
!> Retrieves the nuclear MEP given a list of points.
subroutine get_nuclear_mep(nr_points, centers, mep)

   use molecule_typetype, only: moleculeinfo, atomitem
   use ls_pcm_utils, only: get_molecule 
                                                    
   integer, intent(in)  :: nr_points                           
   real(8), intent(in)  :: centers(3, nr_points)               
   real(8), intent(out) :: mep(nr_points)
   
   type(moleculeinfo)   :: molecule
   real(8)              :: dist                          
   integer              :: i, ipoint
   
   molecule = get_molecule()

   LoopOnAtoms: do i = 1, molecule%nAtoms
       LoopOnPoints: do ipoint = 1, nr_points 
          dist = (molecule%Atom(i)%Center(1) - centers(1, ipoint))**2 + &
                 (molecule%Atom(i)%Center(2) - centers(2, ipoint))**2 + &
                 (molecule%Atom(i)%Center(3) - centers(3, ipoint))**2
          dist = sqrt(dist)
          mep(ipoint) = mep(ipoint) + (molecule%Atom(i)%Charge / dist)
       end do LoopOnPoints
   end do LoopOnAtoms
                                                                 
   !print *, "Debug print of v_nuc at cavity points" 
   !do ipoint = 1, nr_points
   !   print *, "v_nuc(",ipoint,") = ", mep(ipoint)
   !end do
         
end subroutine get_nuclear_mep
         
!> \brief calculates molecular electrostatic potential (MEP) at cavity points 
!> \author R. Di Remigio
!> \date 2014
!> \param nr_points the number of cavity points
!> \param centers the cavity points
!> \param vector MEP or apparent surface charge (ASC) vector
!> \param mat density or Fock matrix
!> \param get_matrix whether to calculate a MEP or the PCM Fock matrix contribution 
!>
!> Driver routine for the calculation of the electronic part of the molecular               
!> electrostatic potential on a certain grid of points {r_i}:
!>     v_el(r_i) = tr(DV_i)
!> tr is the trace operator, D is the density matrix, V^i is the matrix of the 
!> "nuclear attraction" integrals calculated at point r_i of the grid:
!>     V_mu,nu,i =  - <mu|1/|r-r_i||nu>
!> 
!> Written, tested, debugged: R. Di Remigio
!> 
!> RDR 060914 This routine will be used to form both potentials and Fock
!>            matrix contribution for PCM.
!>            matrix is Fock or density matrix, vector is potentials 
!>            or charges vector.
!>            get_matrix logical is present and TRUE:
!>            charges vector as input, Fock matrix contribution as output.
!>            get_matrix logical is absent or is present and FALSE:
!>            density matrix as input, potentials vector as output.                        
subroutine get_electronic_mep(nr_points, centers, vector, mat, get_matrix, setting, lupri, luerr)
   
   use integralinterfacemod, only: II_get_ep_integrals3 
   use matrix_module
   use matrix_operations

   integer,         intent(in)    :: nr_points                            
   real(8),         intent(in)    :: centers(3, nr_points)
   real(8),         intent(out)   :: vector(nr_points) 
   type(matrix)                   :: mat                           
   logical, optional              :: get_matrix 
   type(lssetting), intent(in)    :: setting
   integer,         intent(in)    :: lupri, luerr
   ! Local variables
   logical              :: do_matrix
   integer              :: ipoint
            
   ! Fock matrix contribution or potential calculation?      
   do_matrix = .false. 
   if (present(get_matrix)) then
     if (get_matrix) then
       do_matrix = .true.
     end if
   end if
   
   FockContribution: if (do_matrix) then
      ! Calculate a Fock matrix contribution
   else FockContribution
      ! Just get the electronic MEP here
      call II_get_ep_integrals3(lupri,     &
                                luerr,     &
                                setting,   &
                                centers,   &
                                nr_points, &
                                mat,       &
                                vector)
   end if FockContribution
  
   !print *, "Debug print of v_ele at cavity points"
   !if (do_matrix) then
   !   print *, "Called with get_matrix"
   !else
   !   do ipoint = 1, nr_points
   !     print *, "v_ele(", ipoint,") = ", vector(ipoint)
   !   end do
   !end if

end subroutine get_electronic_mep
      
end module ls_pcm_integrals
