!> @file
!> Some utilities for the PCM interface LSDALTON-side 
module ls_pcm_utils

use iso_c_binding
use molecule_typetype, only: moleculeinfo

implicit none

public init_molecule
public get_molecule

private

type(moleculeinfo) :: molecule

contains

!> \brief initializes global molecule derived type
!> \author R. Di Remigio
!> \date 2014
!>
!> This subroutines retrieves the molecule derived type as initialized
!> at input reading.
subroutine init_molecule(molec)

   type(moleculeinfo), intent(in) :: molec

   molecule = molec

end subroutine init_molecule

!> \brief handle to the global molecule derived type 
!> \author R. Di Remigio
!> \date 2014
!>
!> This function provides a handle to the molecule derived type within
!> the PCM interface LSDALTON-side.
type(moleculeinfo) function get_molecule()

   get_molecule = molecule

end function get_molecule

!> \brief extract Cartesian coordinates from molecule derived type 
!> \author R. Di Remigio
!> \date 2014
!> \param coordinates 3xnAtoms matrix containing the Cartesian coordinates, in Bohr
!>
!> This subroutine is called internally to extract the matrix of Cartesian coordinates
!> from molecule derived type.
subroutine get_coordinates(coordinates)

   use molecule_typetype, only: atomitem
   
   real(c_double), intent(out) :: coordinates(3, *)
   integer :: i
   
   do i = 1, molecule%nAtoms
     coordinates(:, i) = molecule%Atom(i)%Center(:)
   enddo

end subroutine get_coordinates

!> \brief sets number of atoms 
!> \author R. Di Remigio
!> \date 2014
!> \param nr_nuclei number of atoms 
!>
!> This subroutine is called by PCMSolver to set up the number of atoms 
!> in the current molecule, i.e. to dimension various arrays.
subroutine collect_nctot(nr_nuclei) bind(c, name='collect_nctot')

   integer(c_int), intent(out) :: nr_nuclei
   
   nr_nuclei = molecule%nAtoms 

end subroutine collect_nctot
      
!> \brief sets vectors of atomic charges and geometry matrix 
!> \author R. Di Remigio
!> \date 2014
!> \param atomic_charges vector of atomic charges 
!> \param atomic_centers 3xnAtoms matrix containing the Cartesian coordinates, in Bohr
!>
!> This subroutine is called by PCMSolver to set up the geometry to be used 
!> when generating the cavity. The vector of atomic charges is needed in order
!> to provide the right atomic radii. 
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

!> \brief sets point group information
!> \author R. Di Remigio
!> \date 2014
!> \param nr_gen the number of generators for the point group
!> \param gen1 the first generator
!> \param gen2 the second generator 
!> \param gen3 the third generator
!>
!> This subroutine is called by PCMSolver to set point group information 
!> to be used when generating the cavity. 
!> The generators gen1, gen2, gen3 are specified exploiting the mapping
!> between the integers in the 0-7 range, their bitmap representation
!> and the action of the correspoding symmetry operations on the
!> Cartesian frame.
!> Since we are not exploiting symmetry in LSDALTON simply set the
!> number of generators to zero and the three generators to zero too.
subroutine set_point_group(nr_gen, gen1, gen2, gen3) bind(c, name='set_point_group')

   integer(c_int), intent(inout) :: nr_gen 
   integer(c_int), intent(inout) :: gen1, gen2, gen3
   
   nr_gen = 0 
   gen1   = 0
   gen2   = 0
   gen3   = 0

end subroutine set_point_group
 
end module ls_pcm_utils
