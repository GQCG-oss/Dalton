module ls_pcm_utils

use iso_c_binding
use molecule_typetype, only: moleculeinfo

implicit none

public init_molecule

private

type(moleculeinfo) :: molecule

contains

subroutine init_molecule(molec)

   type(moleculeinfo), intent(in) :: molec

   molecule = molec

end subroutine init_molecule

subroutine extract_coordinates(molec, coordinates)

   use molecule_typetype, only: atomitem
   
   type(moleculeinfo), intent(in) :: molec
   real(c_double),    intent(out) :: coordinates(3, *)
   integer :: i
   
   do i = 1, molec%nAtoms
     coordinates(:, i) = molec%Atom(i)%Center(:)
   enddo

end subroutine extract_coordinates

subroutine collect_nctot(nr_nuclei) bind(c, name='collect_nctot')

   integer(c_int), intent(out) :: nr_nuclei
   
   nr_nuclei = molecule%nAtoms 

end subroutine collect_nctot
      
subroutine collect_atoms(atomic_charges, atomic_centers) bind(c, name='collect_atoms')

   use molecule_typetype, only: atomitem
   
   real(c_double), intent(out) :: atomic_charges(*)
   real(c_double), intent(out) :: atomic_centers(3, *)
   
   integer :: i, j, k 
   
   ! Get coordinates
   call extract_coordinates(molecule, atomic_centers)
   ! Get charges      
   do i = 1, molecule%nAtoms
     atomic_charges(i) = molecule%Atom(i)%Charge
   enddo

end subroutine collect_atoms

subroutine set_point_group(nr_gen, gen1, gen2, gen3) bind(c, name='set_point_group')

   integer(c_int), intent(inout) :: nr_gen 
   integer(c_int), intent(inout) :: gen1, gen2, gen3
   
   nr_gen = 1 
   gen1   = 0
   gen2   = 0
   gen3   = 0

end subroutine set_point_group
 
end module ls_pcm_utils
