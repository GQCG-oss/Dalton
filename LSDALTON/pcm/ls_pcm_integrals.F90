module ls_pcm_integrals

implicit none

public get_mep
public get_nuclear_mep
public get_electronic_mep

contains

subroutine get_mep(nr_points, centers, mep, density)

   integer, intent(in)  :: nr_points
   real(8), intent(in)  :: centers(3, nr_points)
   real(8), intent(out) :: mep(nr_points)
   real(8)              :: density(*)
   
   call get_electronic_mep(nr_points, centers, mep, density, .false.)
   call get_nuclear_mep(nr_points, centers, mep)

end subroutine get_mep
             
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
         
subroutine get_electronic_mep(nr_points, centers, vector, matrix, get_matrix)
!      
! Driver routine for the calculation of the electronic part of the molecular 
! electrostatic potential on a certain grid of points {r_i}:
!     v_el(r_i) = tr(DV_i)
! tr is the trace operator, D is the density matrix, V^i is the matrix of the 
! "nuclear attraction" integrals calculated at point r_i of the grid:
!     V_mu,nu,i =  - <mu|1/|r-r_i||nu>
! 
! Written, tested, debugged: R. Di Remigio
! 
! RDR 060914 This routine will be used to form both potentials and Fock
!            matrix contribution for PCM.
!            matrix is Fock or density matrix, vector is potentials 
!            or charges vector.
!            get_matrix logical is present and TRUE:
!            charges vector as input, Fock matrix contribution as output.
!            get_matrix logical is absent or is present and FALSE:
!            density matrix as input, potentials vector as output.                        
!                           
   integer, intent(in)  :: nr_points 
   real(8), intent(in)  :: centers(3, nr_points)
   real(8), intent(out) :: vector(nr_points) 
   real(8)              :: matrix(*)                           
   logical, optional    :: get_matrix 
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
