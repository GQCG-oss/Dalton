program pe_unit_tests

!------------------------------------------------------------------------------
! Unit test of subroutine get_surface
         use double_precision
         use polarizable_embedding

         integer :: natoms
         real(dp), dimension(:,:), allocatable :: all_coords, surf_atoms, cart_new
         real(dp), dimension(:,:), allocatable :: all_chg
         real(dp), dimension(20) :: area
         real(dp), dimension(3,20) :: tri_c
         character(len=2) :: elem_label
         logical :: lexist
         integer :: i, j, luin=1
         real(dp), parameter :: aa2au = 1.8897261249935897d0
         real(dp), dimension(3) :: atom_c
         real(dp) :: r_atom
         integer, dimension(:), allocatable :: ksurf

         open (luin, FILE='molecule.xyz', STATUS='OLD')

         read(luin,*) natoms
         allocate(all_coords(3,natoms))
         allocate(all_chg(1,natoms))
         read(luin,*)
         do i = 1, natoms
             read(luin,*) elem_label, all_coords(1,i), all_coords(2,i), all_coords(3,i)
             all_chg(1,i) = elem2charge(elem_label)
         end do

         close(luin)
!        write(luout,*) 'All coords'
!        do i=1,natoms
!            write (luout,*) i, all_coords(:,i)
!        end do
         all_coords = aa2au*all_coords
         allocate( cart_new(3,natoms) )
         allocate( ksurf(natoms) )
         call get_surface(natoms, all_coords, all_chg, surf_atoms, ksurf, cart_new)

!         atom_c = 0.0d0
!         r_atom = 1.0d0
!         call icosahedron(atom_c, r_atom,area,tri_c)

end program unit_test_polarizable_embedding

