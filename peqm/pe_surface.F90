module pe_surface

    implicit none

    private

    public :: get_surface

!------------------------------------------------------------------------------
#ifdef WORK_CODE

subroutine surface_atoms()
! This routine finds the number of neighbors to each atom. Does not work to find 
! surface atoms. Maybe delete later MNP

     real(dp), dimension(:,:), allocatable :: all_coords
     real(dp), dimension(:,:), allocatable :: all_charges
     real(dp), dimension(:,:), allocatable :: surfatm, surfatm2
     integer :: i, j, neighbors, nsa, k
     real(dp) :: vdwrad_i, vdwrad_j, vdw_ij, dist_rij, charge_check
     real(dp), dimension(3) :: rij

     allocate(all_coords(3,qmnucs+nsites))
     allocate(all_charges(1,qmnucs+nsites))
     allocate(surfatm(3,qmnucs+nsites))
     all_coords(:,1:qmnucs) = Rm
     all_coords(:,qmnucs+1:) = Rs
     all_charges(:,1:qmnucs) = Zm
     all_charges(:,qmnucs+1:) = Zs
     if (print_lvl > 1000) then
        write(luout,*) 'Rs = ', Rs
        write(luout,*) 'Rm = ', Rm
        write(luout,*) 'Rs+Rm = ', all_coords
        write(luout,*) 'Zs = ', Zs
        write(luout,*) 'Zm = ', Zm
        write(luout,*) 'Zs+Zm = ', all_charges
     end if
     nsa = 0
     do i = 1, qmnucs + nsites
        neighbors = 0
        charge_check = all_charges(1,i)
        if (charge_check < 1 ) cycle
        vdwrad_i = charge2vdw(all_charges(1,i))
        do j = 1, qmnucs + nsites
           if (i == j ) cycle
           vdwrad_j = charge2vdw(all_charges(1,j))
           rij = (all_coords(:,j) - all_coords(:,i))
!          The number 1.5 is the radius of a single water molecule
!          must be included in definition of neighboring atoms.
           dist_rij = nrm2(rij) + 1.5*aa2au
           vdw_ij = vdwrad_i + vdwrad_j
!DEBUG           write(luout,*) 'vdwrad_i', vdwrad_i
!DEBUG           write(luout,*) 'vdwrad_j', vdwrad_j
!DEBUG           write(luout,*) 'vdwrad_j+i', vdw_ij
!DEBUG           write(luout,*) 'dist_rij', dist_rij
           if (dist_rij <= vdw_ij ) then
              neighbors = neighbors + 1
           end if
        end do
        if (neighbors < 10 ) then
!DEBUG           write(luout,*) 'I found a surface atom... YEPEE'
           nsa = nsa + 1
           surfatm(:,nsa) = all_coords(:,i)
!DEBUG           write(luout,*) 'Checking surfatm', surfatm(:,nsa)
        end if
        write(luout,*) 'Atom nr:',i,'with',neighbors,'neighbors'
!DEBUG        write(luout,*) 'Atom coords', all_coords(:,i)
     end do
     allocate(surfatm2(3,nsa))
     do k = 1, nsa
        surfatm2(:,k) = surfatm(:,k)
     end do
         write(luout,*) 'Coords of all atoms'
         write(luout,*)  all_coords
         write(luout,*) 'Coords of surface atoms2'
         write(luout,*)  surfatm2
         write(luout,*) 'Number of surface atoms = ', nsa
     if (print_lvl > 1000) then
         write(luout,*) 'Coords of all atoms'
         write(luout,*)  all_coords
         write(luout,*) 'Coords of surface atoms'
         write(luout,*)  surfatm
         write(luout,*) 'Coords of surface atoms2'
         write(luout,*)  surfatm2
     end if
end subroutine surface_atoms

subroutine make_pe_cav()
     
     real(dp), dimension(:,:), allocatable :: all_coords, all_coords_new
     real(dp), dimension(:,:), allocatable :: all_charges
     real(dp), dimension(:,:), allocatable :: surf_atoms
     real(dp), dimension(:,:), allocatable :: tes_p
     real(dp), dimension(:), allocatable :: tes_area
     integer, dimension(:), allocatable :: ksurf ! points to the surface atom
!     real(dp), dimension(:,:) :: surface_atoms
     integer :: natoms, charge_check, j, max_tess_p, i, k
     real(dp) :: r_atom
     real(dp), dimension(:,:), allocatable :: tes_p_new, verts 
     real(dp), dimension(:), allocatable :: tes_area_new
     real(dp) :: r_i, r_j, r_i2, r2, r_ij, r_i_check
     real(dp), dimension(3) :: r_dist
     integer :: l, m, n, max_vert_p
!     real(dp) , dimension(:), allocatable :: n
     logical :: inside

!1) Find the surface atoms of a molecule
     allocate(all_coords(3,qmnucs+nsites(0)))
     allocate(all_charges(1,qmnucs+nsites(0)))
     all_coords(:,1:qmnucs) = Rm
     all_coords(:,qmnucs+1:) = Rs
     all_charges(:,1:qmnucs) = Zm
     all_charges(:,qmnucs+1:) = Zs
     natoms = qmnucs + nsites(0)
     write(luout,*) 'Rm', Rm
     write(luout,*) 'Rs', Rs
     write(luout,*) 'Zm', Zm
     write(luout,*) 'Zs', Zs
     write(luout,*) 'all coords', all_coords
!     do j = 1, size(all_charges, dim=2)
!        charge_check = all_charges(1,j)
!        if (charge_check < 1 ) cycle
!        natoms = natoms + 1 
!     end do
!1)  Find all surface atoms
     allocate( surf_atoms(3,qmnucs+nsites(0)) )
     allocate( ksurf(qmnucs+nsites(0)) )
     surf_atoms = 0.0d0
!     if (natoms < 1000) then
!        surf_atoms = all_coords
!        write(luout,*) 'All coords'
!        write(luout,*) all_coords
!        write(luout,*) 'Surf coords'
!        write(luout,*) surf_atoms
!     else
     allocate(all_coords_new(3,qmnucs+nsites(0))) ! cartesian coordinates in the new basis
     call get_surface(natoms,all_coords,all_charges,surf_atoms,all_coords_new,kk)
     write(luout,*) 'After get_surface: ksurf', ksurf
     write(luout,*) 'After get_surface: all_coords_new', all_coords_new
     write(luout,*) 'After get_surface: surf_atoms', surf_atoms
!2)     Do the actual tesselation on each surface atom!
!     end if
     max_tess_p = (qmnucs + nsites(0))*20 ! 20 if icosahedron is used. 
     max_vert_p = (qmnucs + nsites(0))*12 ! 12 if icosahedron is used. 
     allocate( tes_p(3,max_tess_p) ) 
     allocate( verts(3,max_vert_p) ) 
     allocate( tes_area(max_tess_p) ) 
     tes_p = 0.0d0
     tes_area = 0.0d0
     k = 0
     m = 0
     write(luout,*) 'max_tess_p, qmnucs, nsites(0), size(surf_atoms,dim=2)'
     write(luout,*) max_tess_p, qmnucs, nsites(0), size(surf_atoms,dim=2)
     do j = 1, size(surf_atoms,dim=2)
        r_atom = charge2vdw(all_charges(1,ksurf(j)))*1.170d0  ! add_pcm in divide, to compute SAS
        write(luout,*) 'r_atom for surface atom',j,'oF atom',ksurf(j),'is',r_atom
        write(luout,*) 'surf_atom(:,j)', surf_atoms(:,j)
        call icosahedron(surf_atoms(:,j),r_atom,tes_area(k),tes_p(:,k+1:k+20),verts(:,m+1:m+12))
        write(luout,*) 'Tesselation points for atom ', j
        write(luout,'(3F15.10)') tes_p(:,k+1:k+20)
        write(luout,*) 'Vert points for atom ', j
        write(luout,'(3F15.10)') verts(:,m+1:m+12)
        k = k + 20
        m = m + 12
      end do
!3)  Compute overlap between surface atoms and all other atoms
      allocate( tes_p_new(3,max_tess_p) ) 
      allocate( tes_area_new(max_tess_p) ) 
      tes_p_new = 0.0d0
      tes_area_new = 0.0d0
      m = 1
      k = 0
      do j = 1, size(surf_atoms,dim=2)
         r_j = charge2vdw(all_charges(1,ksurf(j)))*1.170d0 + 2.0d0
         do l = k+1, k+20
            n = 0
            do i = 1, qmnucs + nsites(0) 
               if (i == j) cycle ! same atom
               r_i = charge2vdw(all_charges(1,i))*1.170d0 + 2.0d0
               r_i2 = r_i**2
               r_i_check = (r_i + r_i*0.010d0)**2
!               r_dist = all_coords(:,j) - all_coords(:,i)
!               r_ij = nrm2(r_dist)
!               if ( r_ij > (r_i + r_j) ) cycle ! Atoms are to far apart, we do not need to compute the overlap
               write(luout,*) 'l, i, j =', l, i, j
               write(luout,*) 'tes_p(:,l)', tes_p(:,l)
               write(luout,*) 'all_coords(:,i)', all_coords_new(:,i)
               r2 = ( tes_p(1,l) - all_coords_new(1,i) )**2 + ( tes_p(2,l) - all_coords_new(2,i) )**2 + ( tes_p(3,l) - all_coords_new(3,i) )**2
               write(luout,*) 'r2,r_i2', r2, r_i2
               if ( r2 <= r_i2 ) then
                  n = 0 
                  write(luout,*) 'n =', n
                  exit
               else
                  n = n + 1
                  write(luout,*) 'n =', n
               end if 
            end do
!            if (.not. inside) then
!               tes_p_new(:,m) = tes_p(:,l)
!               tes_area_new(m) = tes_area_new(l)
!               write(luout,*) 'm - 1', m
!               m = m + 1
!            end if
            if (n > 0) then
               tes_p_new(1,m) = tes_p(1,l)
               tes_p_new(2,m) = tes_p(2,l)
               tes_p_new(3,m) = tes_p(3,l)
               tes_area_new(m) = tes_area(l)
               write(luout,*) 'm - 1', m
               m = m + 1
            end if
         end do
         k = k + 20
      end do
      write(luout,*) 'Coordinates of remaining points'
      do j = 1, size(tes_p_new,dim=2)
         write(luout,*) tes_p_new(:,j)
      end do
      write(luout,*) 'Area of remaining points'
      do j = 1, size(tes_p_new,dim=2)
         write(luout,*) tes_area_new(j)
      end do

end subroutine make_pe_cav
#endif
!------------------------------------------------------------------------------
subroutine get_surface(natoms,cart_coords,all_charges,surf_atoms,cart_new,kk)
     
     integer, intent(in) :: natoms
     real(dp), dimension(3,natoms), intent(in) :: cart_coords
     real(dp), dimension(1,natoms), intent(in) :: all_charges
     real(dp), dimension(:,:), allocatable, intent(out) :: surf_atoms, cart_new
     integer, dimension(:), allocatable, intent(out) :: kk
     real(dp), dimension(3,natoms) :: sph_coords
! 1) First all coordinates are transformed to spherical coordinates
     call trans_to_sp(natoms,cart_coords, sph_coords)
! 2) All atoms are divided into boxes 
     allocate( surf_atoms(3,natoms) ) 
     allocate( cart_new(3,natoms) ) 
     allocate( kk(natoms) )
     call divide(natoms,sph_coords,all_charges, cart_coords,surf_atoms,kk)
! 3) The cartesian coordinates of all surface atoms are returned in surf_atoms
!    but they are rotated compared to the original cartisian coordinates of all atoms
!    therefor we must rotate them back

end subroutine get_surface
!------------------------------------------------------------------------------

subroutine trans_to_sp(natoms,cart_coords, sph_coords)


     intrinsic :: acos
     integer, intent(in) :: natoms
     real(dp), dimension(3,natoms), intent(in) :: cart_coords
     real(dp), dimension(3,natoms), intent(out) :: sph_coords
     real(dp), dimension(3) :: kv
     real(dp), dimension(3,natoms) :: ki
     real(dp), dimension(3,3) :: inert_eigv
     real(dp), dimension(3) :: inert_eig, u2
     real(dp), dimension(9) :: work
     real(dp), dimension(6) :: cinrtp
     integer :: i, info, imax, min_eign, imin
     real(dp) :: dotpz, dotpx, dotpy, sinphi, rmax, r, min_eigv, temp, bla


! 1) Find the center of volume
     kv = 0.0d0 
     do i = 1, natoms 
        kv(1) = kv(1) + cart_coords(1,i)
        kv(2) = kv(2) + cart_coords(2,i)
        kv(3) = kv(3) + cart_coords(3,i)
     end do
     
     kv(1) = kv(1)/natoms
     kv(2) = kv(2)/natoms
     kv(3) = kv(3)/natoms

! 2) Find the three principal axis centered in the center of volume

     ki = 0.0d0
     rmax = 0.0d0
     imax = 0
     do i = 1, natoms
        ki(1,i) = cart_coords(1,i) - kv(1)
        ki(2,i) = cart_coords(2,i) - kv(2)
        ki(3,i) = cart_coords(3,i) - kv(3)
        r = ki(1,i) ** 2 + ki(2,i) ** 2 + ki(3,i) ** 2
        if (r > rmax) then
           rmax = r
           imax = i
        end if
     end do
     rmax = sqrt(rmax)
     inert_eigv(1:3,3) = ki(1:3,imax) / rmax

! Do a gram schmidt orthogonalization 
     min_eigv = 1.1d0
     do i = 1,3
        if (abs(inert_eigv(i,3)) < min_eigv) then
           imin = i
           min_eigv = abs(inert_eigv(i,3))
        end if
     end do
     
     u2 = 0.0d0 
     u2(imin) = 1.0d0

     ! temp = u2(1)*inert_eigv(1,3) + u2(2)*inert_eigv(2,3) + u2(3)*inert_eigv(3,3)
     temp = dot(u2, inert_eigv(:,3))
     u2 = u2 - temp*inert_eigv(:,3)
     temp = dot(u2,u2)
     inert_eigv(:,1) = u2 / sqrt(temp)
     
     inert_eigv(1,2) = inert_eigv(2,1)*inert_eigv(3,3) - inert_eigv(3,1)*inert_eigv(2,3)
     inert_eigv(2,2) = inert_eigv(3,1)*inert_eigv(1,3) - inert_eigv(1,1)*inert_eigv(1,3)
     inert_eigv(3,2) = inert_eigv(1,1)*inert_eigv(2,3) - inert_eigv(2,1)*inert_eigv(3,3)

     if (print_lvl > 1000) then
        write(luout,*) 'imax, rmax, imin',imax,rmax,imin
        write(luout,*) 'u z',inert_eigv(:,3)
        write(luout,*) 'u x',inert_eigv(:,1)
        write(luout,*) 'u y',inert_eigv(:,2)
     end if 
! 3) Transform to spherical coordinates
     sph_coords = 0.0d0
     do i = 1, natoms
!sph_coords(1,:) = r
        sph_coords(1,i) = sqrt(ki(1,i)*ki(1,i) + ki(2,i)*ki(2,i) + ki(3,i)*ki(3,i)) 
!sph_coords(2,:) = theta
        dotpz = dot(inert_eigv(:,3),ki(:,i)) 
!TODO, find out why acos(1) does not work when 1 is dp
        bla = real(dotpz/(sph_coords(1,i)))
        if (bla < -1.0d0 ) then
           bla = -1.0d0
        end if
        sph_coords(2,i) = acos(real(dotpz/sph_coords(1,i))) 
!sph_coords(3,:) = phi
        dotpx = dot(inert_eigv(:,1),ki(:,i))
        if (sph_coords(2,i) <= zero    ) then
! We are then directly on the z-axis, so sph_coords(3,i) is redundant
           sph_coords(3,i) = 0.0d0
        else   
           bla = real(dotpx/(sph_coords(1,i)*sin(sph_coords(2,i))))
           if (bla < -1.0d0 ) then
              bla = -1.0d0
           end if
           sph_coords(3,i) = acos(bla)
           dotpy = dot(inert_eigv(:,2),ki(:,i))
           sinphi = dotpy/(sph_coords(1,i)*sin(real(sph_coords(2,i))))
           sph_coords(2,i) = sph_coords(2,i) * 180.d0/pi
           sph_coords(3,i) = sph_coords(3,i) * 180.d0/pi
           if (sinphi < 0.0d0 ) then
              sph_coords(3,i) = 360.0d0 - sph_coords(3,i)
           end if
        end if
     end do
     
     if (print_lvl > 1000) then
         write(luout,*) 'Number of atoms =', natoms
         write(luout,*) 'Center of volume:', kv
         write(luout,*) 'Eigenvector 1 =', inert_eigv(:,1)
         write(luout,*) 'Eigenvector 2 =', inert_eigv(:,2)
         write(luout,*) 'Eigenvector 3 =', inert_eigv(:,3)
         do i = 1, natoms
            write(luout,*) 'Spherical and cartesian coordinates of atom',i,' = ', sph_coords(:,i), cart_coords(:,i)
         end do
     end if
     
end subroutine trans_to_sp
!------------------------------------------------------------------------------
subroutine divide(natoms,sph_coords,charges,cart_coords,cart_coords_surface,kk)

     integer, intent(in) :: natoms
     real(dp), dimension(3,natoms), intent(in) :: sph_coords, cart_coords
     real(dp), dimension(1,natoms), intent(in) :: charges
     real(dp), dimension(:,:),   allocatable, intent(out) :: cart_coords_surface
     integer, dimension(:), allocatable, intent(out) :: kk
     real(dp), dimension(3,natoms) :: cart_coords_new
     integer,  dimension(:,:,:), allocatable :: nbox
     integer,  dimension(:,:),   allocatable :: j_surface ! j_surface(3,n_surface)
     real(dp), dimension(:,:),   allocatable :: r_surface ! r of each surface point
     integer,  dimension(:,:),   allocatable :: k_surface ! which atom gives this surface point
     integer,  dimension(:),     allocatable :: i_temp
     integer,  dimension(:,:,:), allocatable :: kbox
     real(dp), dimension(:,:,:), allocatable :: rbox, rsurf
     real(dp), dimension(3) :: kv
     real(dp) :: rmax, ltheta, lphi, temp, dtheta, dphimin, dphimax
     real(dp) :: r_atom, r_pcm, c_pcm
     real(dp) :: r_dis
     real(dp), parameter :: fac_pcm = 1.17D0, add_pcm = 4.0d0 ! standard values: 1.17D0, 4.0d0 (~ 2 water molecules)
     integer :: ntheta, nphi, i, ibox, jbox, nbox_max, j, itemp
     integer :: k, l, ii, jj, n_surface, jmax, imax
     integer :: max_ntheta = 1000
     logical :: sorted
     real(dp) :: dotp, l_start, cosalpha, r_surf
     real(dp), dimension(:), allocatable :: costheta, sintheta
     real(dp), dimension(:), allocatable :: cosphi, sinphi
     real(dp), dimension(3) :: cart_ijl, cart_surf
     integer :: add_atom, lunit, bla
     real(dp) :: tot_area, area

     rmax = maxval(sph_coords(1,:))
     ntheta = 30*int(rmax*0.3d0 + 1.0d0)
     ntheta = min(ntheta,max_ntheta)
     nphi = 2*ntheta
     ltheta = 180d0/dble(ntheta)
     lphi = ltheta
     if (print_lvl > 1000) then
        write(luout,*) 'rmax   =', rmax
        write(luout,*) 'ntheta =', ntheta
        write(luout,*) 'nphi   =', nphi
        write(luout,*) 'ltheta =', ltheta
        write(luout,*) 'lphi   =', lphi
        write(luout,*) 'natoms =', natoms
     end if
     allocate(nbox(2,nphi,ntheta))
     nbox = 0
! Find the number of atoms in each box
     do i = 1, natoms
         if (sph_coords(2,i) < zero    ) then
! Apparently this if statement is necessacy to place "north pole" atom in box 1,1
            ibox = 1
            jbox = 1
            nbox(1,jbox,ibox) = nbox(1,jbox,ibox) + 1
         else
           ibox = int(sph_coords(2,i)/ltheta + 0.99999d0)
           jbox = int(sph_coords(3,i)/lphi + 0.99999d0)
!           write(luout,*) 'atom number =', i
!           write(luout,*) 'nbox, jbox, ibox =', nbox(1,jbox,ibox), jbox, ibox
           nbox(1,jbox,ibox) = nbox(1,jbox,ibox) + 1
         end if
     end do
     nbox_max = 0
! Find the box with the largest number of atoms
! Is used in the allocation of kbox and rbox
     do j = 1, ntheta
        do i = 1, nphi
           if (nbox_max < nbox(1,i,j)) nbox_max = nbox(1,i,j)
        end do
     end do
     allocate(kbox(nphi,ntheta,nbox_max))
     allocate(rbox(nphi,ntheta,nbox_max))
     nbox = 0
! Assign atom number to kbox and spherical distance in rbox
     do i = 1, natoms
         if (sph_coords(2,i) < zero    ) then
            ibox = 1
            jbox = 1
            nbox(1,jbox,ibox) = nbox(1,jbox,ibox) + 1
            kbox(jbox,ibox,nbox(1,jbox,ibox)) = i
            rbox(jbox,ibox,nbox(1,jbox,ibox)) = sph_coords(1,i)
         else
           ibox = int(sph_coords(2,i)/ltheta + 0.99999d0)
           jbox = int(sph_coords(3,i)/lphi + 0.99999d0)
           nbox(1,jbox,ibox) = nbox(1,jbox,ibox) + 1
           kbox(jbox,ibox,nbox(1,jbox,ibox)) = i
           rbox(jbox,ibox,nbox(1,jbox,ibox)) = sph_coords(1,i)
         end if
     end do
! Sort elements in box n from largest to smallest element
! using a bubblesort algorithm
     do j = 1, ntheta
        jj = 0
        rmin = 1000000.0d0
        rmax = 0.0d0
        do i = 1, nphi
           sorted = .false.
           k = 0
           do while ( .not. sorted )
              sorted = .true.
              k = k + 1
              do l = 1, nbox(1,i,j) - k
                 if (rbox(i,j,l) < rbox(i,j,l+1)) then
                    temp  = rbox(i,j,l)
                    itemp = kbox(i,j,l)
                    rbox(i,j,l) = rbox(i,j,l+1)
                    kbox(i,j,l) = kbox(i,j,l+1)
                    rbox(i,j,l+1) = temp
                    kbox(i,j,l+1) = itemp
                    sorted = .false.
                 end if
              end do
           end do  
!#ifdef UNIT_TEST
!           if (print_lvl > 100) write(luout,*) '*** box no.',i,j
!           if (nbox(1,i,j) > 0) then
!             if (print_lvl > 100) then
!              do l = 1, nbox(1,i,j)
!                 write(luout,'(A,2I8,F10.5,F10.2)') 'l,kbox,rbox= ',l, kbox(i,j,l), rbox(i,j,l), charges(1,kbox(i,j,l))
!              end do
!             end if
              rmin = min( rmin, rbox(i,j,1) )
              rmax = max( rmax, rbox(i,j,1) )
              jj = jj + nbox(1,i,j)
!           end if
!#endif
        end do
!#ifdef UNIT_TEST
!        write(luout,*) 'Total number of atoms in section',j,jj, rmin, rmax
        dtheta = rmax*pi/dble(ntheta)
        if ( j > (ntheta/2) ) then
           dphimax   = rmax * sin( pi*dble(j-1)/dble(ntheta) ) * 2.0d0*pi/dble(nphi)
           dphimin   = rmax * sin( pi*dble(j)/dble(ntheta) ) * 2.0d0*pi/dble(nphi)
        else
           dphimax   = rmax * sin( pi*dble(j)/dble(ntheta) ) * 2.0d0*pi/dble(nphi)
           dphimin   = rmax * sin( pi*dble(j-1)/dble(ntheta) ) * 2.0d0*pi/dble(nphi)
        end if
!        write(luout,'(A,3F10.4)') ' -- dtheta, dphi for rmax',dtheta,dphimin,dphimax
!#endif
     end do
! rbox(i,j,1) now holds a (possible) surface atom
! make list of surface atos
! - the first one is the atom with the greatest radial value

     allocate( j_surface(3,natoms) )
     allocate( r_surface(nphi,0:ntheta) )
     allocate( k_surface(nphi,0:ntheta) )
     r_surface = 0.0d0
     k_surface = 0

! North pole atom
     !write(luout,*) 'nbox(2,1,1)', nbox(2,1,1)
     i = 1
     j = 1
     l = 1
     nbox(2,i,j) = l
     n_surface = 1
     j_surface(1,n_surface) = l ! j_surface points to the atom in kbox(i,j,l)
     j_surface(2,n_surface) = i
     j_surface(3,n_surface) = j
     r_atom = rbox(i,j,l)
     c_pcm = charges(1,kbox(i,j,l))
     r_pcm = charge2vdw(c_pcm) * fac_pcm + add_pcm

     allocate( costheta(ntheta) )
     allocate( sintheta(ntheta) )
     allocate( cosphi(nphi) )
     allocate( sinphi(nphi) )
     do jj = 1, ntheta
        costheta(jj) = cos(dble(jj)*pi/dble(ntheta))
        sintheta(jj) = sin(dble(jj)*pi/dble(ntheta))
     end do
     do ii = 1, nphi
        cosphi(ii) = cos(dble(ii)*2.0d0*pi/dble(nphi))
        sinphi(ii) = sin(dble(ii)*2.0d0*pi/dble(nphi))
     end do
     r_surface(1,0) = r_atom + r_pcm
     k_surface(1,0) = n_surface
     !write(luout,*) 'jj, r_surface(jj)',  0,r_surface(1,0), k_surface(1,0)
     do jj = 1,ntheta
        r_dis = (rbox(1,1,1)*costheta(jj))**2 - (rbox(1,1,1)**2 - r_pcm**2)
        !write(luout,*) 'r_dis for "north pole" atom is:', r_dis,jj
        if ( r_dis < 0.0d0 ) then
           exit 
        else
           r_surface(:,jj) = rbox(1,1,1)*costheta(jj) + sqrt(r_dis)
           k_surface(:,jj) = n_surface
           !write(luout,*) 'jj, r_surface(jj)', jj, r_surface(1,jj),k_surface(1,jj)
        end if
     end do
! End north pole atom

     do j = 1, ntheta
     do i = 1, nphi
        l_start = nbox(2,i,j) + 1
        do l = l_start, nbox(1,i,j)
           add_atom = 0
           r_atom = rbox(i,j,l)
           c_pcm = charges(1,kbox(i,j,l))
           r_pcm = charge2vdw(c_pcm) * fac_pcm + add_pcm
! get the cartesian coordinates from atom i,j,l
           cart_ijl(1) = sin( sph_coords(2,kbox(i,j,l))*(pi/180d0) ) * cos( sph_coords(3,kbox(i,j,l))*(pi/180d0) ) 
           cart_ijl(2) = sin( sph_coords(2,kbox(i,j,l))*(pi/180d0) ) * sin( sph_coords(3,kbox(i,j,l))*(pi/180d0) ) 
           cart_ijl(3) = cos( sph_coords(2,kbox(i,j,l))*(pi/180d0) ) 
           !write(luout,*) 'Next atom',l,i,j,r_atom
           !write(luout,*) 'pcm',c_pcm, r_pcm
           !write(luout,*) 'unit vector',cart_ijl(1:3)
           do jj = 1, ntheta
              do ii = 1, nphi
! get the cartesian coordinates from surface point nphi,ntheta
                   cart_surf(1) = sintheta(jj) * cosphi(ii)
                   cart_surf(2) = sintheta(jj) * sinphi(ii)
                   cart_surf(3) = costheta(jj)
                   !write(luout,*) 'surface point',ii,jj
                   !if (j .eq. 1) then
                   !   write(luout,*) 'unit vector  ',cart_surf(1:3)
                   !end if
! Find the angle alpha between cart_ijl and cart_surf
                   cosalpha = dot( cart_ijl(:), cart_surf(:) )  
                   r_dis = (r_atom*cosalpha)**2 - (r_atom**2 - r_pcm**2)
                   !write (luout,*) 'cosalpha, r_dis',cosalpha, r_dis
                   if ( r_dis < 0.0d0 ) then
                      cycle
                   else
                      r_surf = r_atom*cosalpha + sqrt(r_dis)
                      if ( r_surf < r_surface(ii,jj)) then
                         cycle
                      else
                         add_atom = add_atom + 1
                         r_surface(ii,jj) = r_surf
                         k_surface(ii,jj) = n_surface + 1
!                         if (print_lvl > 1000) then
!                             write(luout,*) 'Surface atom nr.', n_surface+1, 'add_atom', add_atom
!                             write(luout,*) ii, jj, r_atom,r_surface(ii,jj)
!                         end if
                      end if
                   end if
              end do
           end do
           if (add_atom > 0 ) then
! Add atom      
              nbox(2,i,j) = l
              n_surface = n_surface + 1
              j_surface(1,n_surface) = l ! j_surface points to the atom in kbox(i,j,l)
              j_surface(2,n_surface) = i
              j_surface(3,n_surface) = j
           end if
        end do
     end do
     end do
     write(luout,*) 'Total number of surface atoms:', n_surface


! Eliminate non-surface atoms
     add_atom = 1 ! north-pole atom
k_loop: do k = 2, n_surface
        l = j_surface(1,k)
        i = j_surface(2,k)
        j = j_surface(3,k)
           r_atom = rbox(i,j,l)
           c_pcm = charges(1,kbox(i,j,l))
           r_pcm = charge2vdw(c_pcm) * fac_pcm + add_pcm
! get the cartesian coordinates from atom i,j,l
           cart_ijl(1) = sin( sph_coords(2,kbox(i,j,l))*(pi/180d0) ) * cos( sph_coords(3,kbox(i,j,l))*(pi/180d0) ) 
           cart_ijl(2) = sin( sph_coords(2,kbox(i,j,l))*(pi/180d0) ) * sin( sph_coords(3,kbox(i,j,l))*(pi/180d0) ) 
           cart_ijl(3) = cos( sph_coords(2,kbox(i,j,l))*(pi/180d0) ) 
           !write(luout,*) 'Next atom',l,i,j,r_atom
           !write(luout,*) 'pcm',c_pcm, r_pcm
           !write(luout,*) 'unit vector',cart_ijl(1:3)
           do jj = 1, ntheta
              do ii = 1, nphi
                 if (r_atom + r_pcm < r_surface(ii,jj)) cycle
! get the cartesian coordinates from surface point nphi,ntheta
                   cart_surf(1) = sintheta(jj) * cosphi(ii)
                   cart_surf(2) = sintheta(jj) * sinphi(ii)
                   cart_surf(3) = costheta(jj)
! Find the angle alpha between cart_ijl and cart_surf
                   cosalpha = dot( cart_ijl(:), cart_surf(:) )  
                   r_dis = (r_atom*cosalpha)**2 - (r_atom**2 - r_pcm**2)
                   if ( r_dis < 0.0d0 ) then
                      cycle
                   else
                      r_surf = r_atom*cosalpha + sqrt(r_dis)
                      if ( r_surf < 0.99999d0*r_surface(ii,jj)) then
                         cycle
                      else
                         add_atom = add_atom + 1
                        ! if (print_lvl > 1000) then
                         !    write(luout,*) 'Surface atom nr.', k, 'add_atom', add_atom
                          !   write(luout,*) ii, jj, r_surf,r_surface(ii,jj)
                        ! end if
                         if (add_atom < k) then
                            j_surface(1,add_atom) = j_surface(1,k)
                            j_surface(2,add_atom) = j_surface(2,k)
                            j_surface(3,add_atom) = j_surface(3,k)
                         end if
                         cycle k_loop
                      end if
                   end if
              end do
           end do
      end do k_loop
      write (luout,*) 'Revised number of surface atoms',add_atom, ' from', n_surface
      allocate ( i_temp(n_surface) )
      i_temp = 0
      i_temp(1) = 1
      do jj = 1, ntheta
         do ii = 1, nphi
           i_temp( k_surface(ii,jj) ) = i_temp( k_surface(ii,jj) ) + 1 
         end do
      end do
!      write(luout,*) 'the real surface atoms:'
!      write(luout,'(20I5)') i_temp(1:n_surface)
      itemp = 0
      do i = 1, n_surface
         if ( i_temp(i) > 0 ) itemp = itemp + 1
      end do
      write(luout,*) 'Revised number of surface atoms from k_surface',itemp

      n_surface = add_atom
          allocate( cart_coords_surface(3,n_surface) )
          allocate( kk(n_surface) )
          write(luout,*) 'Surface atoms'
          do i = 1, n_surface
             kk(i) = kbox(j_surface(2,i), j_surface(3,i), j_surface(1,i) )
             cart_coords_surface(1,i) = sph_coords(1,kk(i)) &
                                      & * sin(sph_coords(2,kk(i))*pi/180.0D0) &
                                      & * cos(sph_coords(3,kk(i))*pi/180.0D0)
             cart_coords_surface(2,i) = sph_coords(1,kk(i)) &
                                      & * sin(sph_coords(2,kk(i))*pi/180.0D0) &
                                      & * sin(sph_coords(3,kk(i))*pi/180.0D0)
             cart_coords_surface(3,i) = sph_coords(1,kk(i)) &
                                      & * cos(sph_coords(2,kk(i))*pi/180.0D0)
          end do
!          write(luout,*) 'Surface points'
          write( luout,*) 0.0d0,0.0d0,r_surface(1,0), area, 1
          i_temp = 0
          bla = 1
             do jj = 1, ntheta
                do ii = 1, nphi
                   cart_surf(1) = r_surface(ii,jj) * sintheta(jj) * cosphi(ii)
                   cart_surf(2) = r_surface(ii,jj) * sintheta(jj) * sinphi(ii)
                   cart_surf(3) = r_surface(ii,jj) * costheta(jj)
                   bla = bla + 1
!                   write(luout,*) cart_surf(1:3), area, bla
                end do
             end do  
          write(luout,*) 'Number of surface point for each surface atom'
          write(luout,'(50I2)') i_temp(1:n_surface)
          close( lunit )

      if (print_lvl > 10 ) then
          write(luout,*) 'All coords'
          do i = 1, natoms
                 cart_coords_new(1,i) = sph_coords(1,i) &
                                      & * sin(sph_coords(2,i)*pi/180.0D0) &
                                      & * cos( sph_coords(3,i)*pi/180.0D0 ) 
                 cart_coords_new(2,i) = sph_coords(1,i) & 
                                      & * sin(sph_coords(2,i)*pi/180.0D0) & 
                                      & * sin( sph_coords(3,i)*pi/180.0D0 ) 
                 cart_coords_new(3,i) = sph_coords(1,i) &
                                      & * cos(sph_coords(2,i)*pi/180.0D0) 
                 write(luout,'(3F15.10)') cart_coords_new(1:3,i)
          end do 
      end if

end subroutine divide

!------------------------------------------------------------------------------

#ifdef WORK_CODE
subroutine icosahedron(atom_c, r_atom,area,tri_c,verts)

      real(dp), dimension(3), intent(in) :: atom_c
      real(dp), intent(in) :: r_atom
      real(dp), dimension(20), intent(out) :: area
      real(dp), dimension(3,20), intent(out) :: tri_c
      real(dp), parameter :: phi = 26.56505d0
      real(dp) :: phia, theta
      real(dp), dimension(3,12), intent(out) :: verts
      integer, dimension(1:3,20) :: faces ! Holds pointers to the atoms that define the faces 
      integer :: i
      real(dp) :: d, e, f, d2, e2, f2
!
      phia = phi * pi/180.0d0

      verts = 0.0d0
      !North pole coordinate
      verts(1,1) = atom_c(1)
      verts(2,1) = atom_c(2)
      verts(3,1) = atom_c(3) + r_atom 

      theta = 0.0d0

      do i = 2, 6
         verts(1,i) = atom_c(1) + r_atom * cos(theta) * cos(phia)
         verts(2,i) = atom_c(2) + r_atom * sin(theta) * cos(phia)
         verts(3,i) = atom_c(3) + r_atom * sin(phia) 
         theta = theta + pi * 72.0d0 / 180.0d0
      end do

      theta = pi * 36.0d0 / 180.0d0
 
      do i = 7, 11
         verts(1,i) = atom_c(1) + r_atom * cos(theta) * cos(-phia)
         verts(2,i) = atom_c(2) + r_atom * sin(theta) * cos(-phia)
         verts(3,i) = atom_c(3) + r_atom * sin(-phia) 
         theta = theta + pi * 72.0d0 / 180.0d0
      end do

      ! South pole coordinate
      verts(1,12) = atom_c(1)
      verts(2,12) = atom_c(2)
      verts(3,12) = atom_c(3) - r_atom 
      
!      write(luout,*) 'Coordinates of a icosahedron'
!      do i = 1, 12
!         write(luout,*) verts(:,i)
!      end do

! calculate the center and area of each triangle 
      faces(1,1) = 1
      faces(2,1) = 2
      faces(3,1) = 3
      faces(1,2) = 1
      faces(2,2) = 3
      faces(3,2) = 4
      faces(1,3) = 1
      faces(2,3) = 4 
      faces(3,3) = 5
      faces(1,4) = 1
      faces(2,4) = 5
      faces(3,4) = 6 
      faces(1,5) = 1
      faces(2,5) = 6
      faces(3,5) = 2 
      faces(1,6) = 12
      faces(2,6) = 7
      faces(3,6) = 8
      faces(1,7) = 12
      faces(2,7) = 8
      faces(3,7) = 9
      faces(1,8) = 12
      faces(2,8) = 9
      faces(3,8) = 10
      faces(1,9) = 12
      faces(2,9) = 10
      faces(3,9) = 11
      faces(1,10) = 12
      faces(2,10) = 11
      faces(3,10) = 7 
      faces(1,11) = 2 
      faces(2,11) = 3 
      faces(3,11) = 7 
      faces(1,12) = 3 
      faces(2,12) = 4 
      faces(3,12) = 8 
      faces(1,13) = 4 
      faces(2,13) = 5 
      faces(3,13) = 9 
      faces(1,14) = 5 
      faces(2,14) = 6 
      faces(3,14) = 10
      faces(1,15) = 6 
      faces(2,15) = 2 
      faces(3,15) = 11
      faces(1,16) = 7 
      faces(2,16) = 8 
      faces(3,16) = 3 
      faces(1,17) = 8 
      faces(2,17) = 9 
      faces(3,17) = 4 
      faces(1,18) = 9 
      faces(2,18) = 10
      faces(3,18) = 5 
      faces(1,19) = 10
      faces(2,19) = 11
      faces(3,19) = 6 
      faces(1,20) = 11
      faces(2,20) = 7 
      faces(3,20) = 2

      do i = 1, 20 ! loop over faces
        !d = y1*z2*1 + z1*1*y3 + 1*y2*z3 - y3*z2*1 - z3*1*y1 - 1*y2*z1
         d = verts(2,faces(1,i))*verts(3,faces(2,i)) + verts(3,faces(1,i))*verts(2,faces(3,i))&
         & + verts(2,faces(2,i))*verts(3,faces(3,i)) - verts(2,faces(3,i))*verts(3,faces(2,i))&
         & - verts(3,faces(3,i))*verts(2,faces(1,i)) - verts(2,faces(2,i))*verts(3,faces(1,i))      
        !e = z1*x2*1 + x1*1*z3 + 1*z2*x3 - z3*x2*1 - x3*1*z1 - 1*z2*x1
         e = verts(3,faces(1,i))*verts(1,faces(2,i)) + verts(1,faces(1,i))*verts(3,faces(3,i))&
         & + verts(3,faces(2,i))*verts(1,faces(3,i)) - verts(3,faces(3,i))*verts(1,faces(2,i))&
         & - verts(1,faces(3,i))*verts(3,faces(1,i)) - verts(3,faces(2,i))*verts(1,faces(1,i))      
        !f = x1*y2*1 + y1*1*x3 + 1*x2*y3 - x3*y2*1 - y3*1*x1 - 1*x2*y1
         f = verts(1,faces(1,i))*verts(2,faces(2,i)) + verts(2,faces(1,i))*verts(1,faces(3,i))&
         & + verts(1,faces(2,i))*verts(2,faces(3,i)) - verts(1,faces(3,i))*verts(2,faces(2,i))&
         & - verts(2,faces(3,i))*verts(1,faces(1,i)) - verts(1,faces(2,i))*verts(2,faces(1,i))      
         d2 = d**2
         e2 = e**2
         f2 = f**2
         area(i) = 0.50d0 * sqrt(d2 + e2 + f2)
         tri_c(1,i) = (1.0d0/3.0d0) * ( verts(1,faces(1,i)) + verts(1,faces(2,i)) + verts(1,faces(3,i)) ) 
         tri_c(2,i) = (1.0d0/3.0d0) * ( verts(2,faces(1,i)) + verts(2,faces(2,i)) + verts(2,faces(3,i)) ) 
         tri_c(3,i) = (1.0d0/3.0d0) * ( verts(3,faces(1,i)) + verts(3,faces(2,i)) + verts(3,faces(3,i)) ) 
      end do
 
!      write(luout,*) 'Triangular area'
!      do i = 1, 20
!         write(luout,*) area(i)
!      end do
!      write(luout,*) 'Triangular center'
!      do i = 1, 20
!         write(luout,*) tri_c(:,i)
!      end do
       
end subroutine icosahedron
#endif
!------------------------------------------------------------------------------

end module pe_surface
