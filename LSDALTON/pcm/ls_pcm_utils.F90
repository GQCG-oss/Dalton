module ls_pcm_utils

use iso_c_binding
use molecule_typetype, only: moleculeinfo

implicit none

public init_molecule
public get_molecule
public get_coordinates

private

type(moleculeinfo) :: molecule

contains

subroutine init_molecule(molec)

   type(moleculeinfo), intent(in) :: molec

   molecule = molec

end subroutine init_molecule

type(moleculeinfo) function get_molecule()

   get_molecule = molecule

end function get_molecule

subroutine get_coordinates(coordinates)

   use molecule_typetype, only: atomitem
   
   real(c_double),    intent(out) :: coordinates(3, *)
   integer :: i
   
   do i = 1, molecule%nAtoms
     coordinates(:, i) = molecule%Atom(i)%Center(:)
   enddo

end subroutine get_coordinates

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
   call get_coordinates(atomic_centers)
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
