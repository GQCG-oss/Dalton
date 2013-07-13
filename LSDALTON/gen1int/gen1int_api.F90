!> \brief module of API of Gen1Int interface
!> \author Bin Gao
!> \date 2012-05-23
module gen1int_api

  ! Fortran 90 module of Gen1Int library
  use gen1int, Gen1IntAPIPropView => OnePropView,               &
               Gen1IntAPIPropGetNumProp => OnePropGetNumProp,   &
               Gen1IntAPIPropGetSymmetry => OnePropGetSymmetry, &
               Gen1IntAPIPropDestroy => OnePropDestroy,         &
               Gen1IntAPIGeoTreeDestroy => GeomTreeDestroy,     &
               Gen1IntAPIGeoTreeView => GeomTreeView
  use precision

  implicit none

  public :: Gen1IntAPIPropCreate
  public :: Gen1IntAPIPropView
  public :: Gen1IntAPIPropGetNumProp
  public :: Gen1IntAPIPropGetSymmetry
  public :: Gen1IntAPIPropGetIntExpt
  public :: Gen1IntAPIPropDestroy

  public :: Gen1IntAPIGeoTreeCreate
  public :: Gen1IntAPIGeoTreeView
  public :: Gen1IntAPIGeoTreeDestroy

  contains

  !> \brief creates the operator of property integrals
  !> \author Bin Gao
  !> \date 2012-05-23
  !> \param mol_info contains the information of molecule
  !> \param gto_type specifies the type of GTOs, should be either NON_LAO (non London atomic
  !>        orbital), LONDON (London atomic orbital, LAO), or ROT_LAO (rotational LAO), only
  !>        NON_LAO implemented
  !> \param prop_name is the name of property integrals
  !> \param order_mom is the order of multipole moments
  !> \param order_mag_bra is the order of partial magnetic derivatives on bra center, not implemented
  !> \param order_mag_ket is the order of partial magnetic derivatives on ket center, not implemented
  !> \param order_mag_total is the order of total magnetic derivatives, not implemented
  !> \param order_ram_bra is the order of partial derivatives w.r.t. the total rotational
  !>        angular momentum on bra center, not implemented
  !> \param order_ram_ket is the order of partial derivatives w.r.t. the total rotational
  !>        angular momentum on ket center, not implemented
  !> \param order_ram_total is the order of total derivatives w.r.t. the total rotational
  !>        angular momentum, not implemented
  !> \param order_geo_bra is the order of partial geometric derivatives on bra center
  !> \param order_geo_ket is the order of partial geometric derivatives on ket center
  !> \param add_sr is for scalar-relativistic (SR) correction, not implemented
  !> \param add_so is for spin-orbit (SO) correction, not implemented
  !> \param add_london transforms the operator by the LAO type gauge-including projector, not implemented
  !> \param io_viewer is the logical unit number of the viewer
  !> \return one_prop is the operator of property integrals
  subroutine Gen1IntAPIPropCreate(mol_info, gto_type,           &
                                  prop_name, order_mom,         &
                                  order_mag_bra, order_mag_ket, &
                                  order_mag_total,              &
                                  order_ram_bra, order_ram_ket, &
                                  order_ram_total,              &
                                  order_geo_bra, order_geo_ket, &
                                  add_sr, add_so, add_london,   &
                                  one_prop, io_viewer)
    use molecule_type, only: MOLECULEINFO,            &
                             molecule_get_num_atoms,  &
                             molecule_get_atom_coord, &
                             molecule_get_atom_charge
    type(MOLECULEINFO), intent(in) :: mol_info
    integer, optional, intent(in) :: gto_type
    character*(*), intent(in) :: prop_name
    integer, optional, intent(in) :: order_mom
    integer, optional, intent(in) :: order_mag_bra
    integer, optional, intent(in) :: order_mag_ket
    integer, optional, intent(in) :: order_mag_total
    integer, optional, intent(in) :: order_ram_bra
    integer, optional, intent(in) :: order_ram_ket
    integer, optional, intent(in) :: order_ram_total
    integer, optional, intent(in) :: order_geo_bra
    integer, optional, intent(in) :: order_geo_ket
    logical, optional, intent(in) :: add_sr
    logical, optional, intent(in) :: add_so
    logical, optional, intent(in) :: add_london
    type(one_prop_t), intent(inout) :: one_prop
    integer, intent(in) :: io_viewer
!FIXME: gets information of dipole origin from mol_info?
    real(realk), parameter :: DIPORG(3) = (/0.0,0.0,0.0/)  !coordinates of dipole origin
    integer num_atoms                             !number of atoms
    real(realk), allocatable :: coord_atoms(:,:)  !coordinates of atoms
    real(realk), allocatable :: charge_atoms(:)   !charges of atoms
    integer ierr                                  !error information
    select case (trim(prop_name))
    ! one-electron Hamiltonian
    case (INT_ONE_HAMIL)
      ! gets the coordinates and charges of atoms
      call molecule_get_num_atoms(mol_info=mol_info, num_atoms=num_atoms)
      allocate(coord_atoms(3,num_atoms), stat=ierr)
      if (ierr/=0) then
        stop "Gen1IntAPIPropCreate>> failed to allocate coord_atoms!"
      end if
      allocate(charge_atoms(num_atoms), stat=ierr)
      if (ierr/=0) then
        stop "Gen1IntAPIPropCreate>> failed to allocate charge_atoms!"
      end if
      call molecule_get_atom_coord(mol_info=mol_info,       &
                                   num_atoms=num_atoms,     &
                                   coord_atoms=coord_atoms, &
                                   lupri=io_viewer)
      call molecule_get_atom_charge(mol_info=mol_info,         &
                                    num_atoms=num_atoms,       &
                                    charge_atoms=charge_atoms, &
                                    lupri=io_viewer)
      call OnePropCreate(prop_name=INT_ONE_HAMIL,  &
                         one_prop=one_prop,        &
                         info_prop=ierr,           &
                         coord_nuclei=coord_atoms, &
                         charge_nuclei=charge_atoms)
      deallocate(coord_atoms)
      deallocate(charge_atoms)
    ! Cartesian multipole moments
    case (INT_CART_MULTIPOLE)
      call OnePropCreate(prop_name=INT_CART_MULTIPOLE, &
                         one_prop=one_prop,            &
                         info_prop=ierr,               &
                         dipole_origin=DIPORG,         &
                         order_mom=order_mom)
    ! overlap integrals
    case (INT_OVERLAP)
      call OnePropCreate(prop_name=INT_OVERLAP, &
                         one_prop=one_prop,     &
                         info_prop=ierr)
    ! kinetic energy integrals
    case (INT_KIN_ENERGY)
      call OnePropCreate(prop_name=INT_KIN_ENERGY, &
                         one_prop=one_prop,        &
                         info_prop=ierr)
    ! one-electron potential energy integrals
    case (INT_POT_ENERGY)
      ! gets the coordinates and charges of atoms
      call molecule_get_num_atoms(mol_info=mol_info, num_atoms=num_atoms) 
      allocate(coord_atoms(3,num_atoms), stat=ierr)
      if (ierr/=0) then
        stop "Gen1IntAPIPropCreate>> failed to allocate coord_atoms!"
      end if
      allocate(charge_atoms(num_atoms), stat=ierr)
      if (ierr/=0) then
        stop "Gen1IntAPIPropCreate>> failed to allocate charge_atoms!"
      end if
      call molecule_get_atom_coord(mol_info=mol_info,       &
                                   num_atoms=num_atoms,     &
                                   coord_atoms=coord_atoms, &
                                   lupri=io_viewer)
      call molecule_get_atom_charge(mol_info=mol_info,         &
                                    num_atoms=num_atoms,       &
                                    charge_atoms=charge_atoms, &
                                    lupri=io_viewer)
      call OnePropCreate(prop_name=INT_POT_ENERGY, &
                         one_prop=one_prop,        &
                         info_prop=ierr,           &
                         coord_nuclei=coord_atoms, &
                         charge_nuclei=charge_atoms)
      deallocate(coord_atoms)
      deallocate(charge_atoms)
    case default
      write(io_viewer,999) "unknown property "//trim(prop_name)//"!"
      stop
    end select
    if (ierr/=0) then
      write(io_viewer,999) "failed to create operator of "//trim(prop_name)//"!"
      stop
    end if
    ! sets partial geometric derivatives
    if (present(order_geo_bra) .or. present(order_geo_ket))   &
      call OnePropSetPartialGeom(one_prop=one_prop,           &
                                 order_geo_bra=order_geo_bra, &
                                 order_geo_ket=order_geo_ket)
    ! sets magnetic derivatives
    if (present(order_mag_total) .or.                 &
        present(order_mag_bra) .or.                   &
        present(order_mag_ket))                       &
      call OnePropSetMag(one_prop=one_prop,           &
                         order_mag=order_mag_total,   &
                         order_mag_bra=order_mag_bra, &
                         order_mag_ket=order_mag_ket)
    ! sets derivatives w.r.t. total rotational angular momentum
    if (present(order_ram_total) .or.                 &
        present(order_ram_bra) .or.                   &
        present(order_ram_ket))                       &
      call OnePropSetRAM(one_prop=one_prop,           &
                         order_ram=order_ram_total,   &
                         order_ram_bra=order_ram_bra, &
                         order_ram_ket=order_ram_ket)
    ! sets the type of GTOs
    if (present(gto_type)) then
      call OnePropSetGTO(one_prop=one_prop, &
                         gto_type=gto_type, &
                         info_prop=ierr)
      if (ierr/=0) then
        stop "Gen1IntAPIPropCreate>> invalid type of GTOs!"
      end if
    end if
999 format("Gen1IntAPIPropCreate>> ",A)
  end subroutine Gen1IntAPIPropCreate

  !> \brief evaluates the integral matrices and/or expectation values
  !> \author Bin Gao
  !> \date 2012-05-23
  !> \param ao_item_bra contains the AO item on bra center
  !> \param ao_item_ket contains the AO item on ket center
  !> \param same_braket indicates if the AO items are the same on bra and ket centers
  !> \param do_screen indicates if using screen scheme?
  !> \param one_prop contains the information of one-electron property integrals
  !> \param geom_tree contains the information of N-ary tree for total geometric derivatives
  !> \param geom_type is the type of returned total geometric derivatives, UNIQUE_GEO is for unique total
  !>        geometric derivatives, and REDUNDANT_GEO for redundant total geometric derivatives;
  !>        for instance, (xx,xy,yy,xz,yz,zz) are unique while (xx,yx,zx,xy,yy,zy,xz,yz,zz) are
  !>        redundant, note that the "triangular" total geometric derivatives could be obtained
  !>        from the unique total geometric derivatives by giving \var(max_num_cent)=\var(order_geo_total)
  !> \param write_ints indicates if writing integral matrices on file
  !> \param ao_dens contains the AO density matrices
  !> \param write_expt indicates if writing expectation values on file
  !> \param io_viewer is the logical unit number of the viewer
  !> \param level_print is the level of print
  !> \return val_ints contains the integral matrices
  !> \return val_expt contains the expectation values
  !> \note the arrangement of var(val_ints) and \var(val_expt) will be in the order of
  !>       \var(order_mom), \var(order_mag_bra), ..., \var(order_geo_total), and each of
  !>       them is arranged in the order of (xx,xy,yy,xz,yz,zz) or (xx,yx,zx,xy,yy,zy,xz,yz,zz),
  !>       see Gen1Int library manual, for instance Section 2.2;
  !>       \var(val_expt) should be zero by users before calculations
  subroutine Gen1IntAPIPropGetIntExpt(ao_item_bra, ao_item_ket,        &
                                      same_braket, do_screen,          &
                                      one_prop, geom_tree, geom_type,  &
                                      val_ints, write_ints,            &
                                      ao_dens, val_expt, write_expt,   &
                                      io_viewer, level_print)
    use AO_Type, only: AOITEM
    use Matrix_module, only: Matrix
    ! AO sub-shells
    use gen1int_shell
    type(AOITEM), intent(in) :: ao_item_bra
    type(AOITEM), intent(in) :: ao_item_ket
    logical, intent(in) :: same_braket
    logical, intent(in) :: do_screen
    type(one_prop_t), intent(in) :: one_prop
    type(geom_tree_t), optional, intent(inout) :: geom_tree
    integer, optional, intent(in) :: geom_type
    type(Matrix), optional, intent(inout) :: val_ints(:)
    logical, optional, intent(in) :: write_ints
    type(Matrix), optional, intent(in) :: ao_dens(:)
    real(realk), optional, intent(inout) :: val_expt(:,:)
    logical, optional, intent(in) :: write_expt
    integer, intent(in) :: io_viewer
    integer, intent(in) :: level_print
    integer num_prop        !number of property integrals
    integer order_geo       !order of total geometric derivatives
    integer p_geom_type     !type of returned total geometric derivatives
    integer num_return_geo  !number of returned total geometric derivatives
    integer idx_path        !index of current path of N-ary tree
    integer num_paths       !total number of different paths of N-ary tree
    integer ipath           !incremental recorder over different paths
    ! gets the number of property integrals
    call OnePropGetNumProp(one_prop=one_prop, num_prop=num_prop)
    ! sets the type of total geometric derivatives
    if (present(geom_tree)) then
      ! gets the order of total geometric derivatives
      call GeomTreeGetOrder(geom_tree=geom_tree, order_geo=order_geo)
      ! gets the type of total geometric derivatives
      if (order_geo>1) then
        if (present(geom_type)) then
          select case (geom_type)
          case (REDUNDANT_GEO)
            p_geom_type = REDUNDANT_GEO
          case default
            p_geom_type = UNIQUE_GEO
          end select
        ! default are unique total geometric derivatives
        else
          p_geom_type = UNIQUE_GEO
        end if
      ! the first order total geometric derivatives are unique
      else
        p_geom_type = UNIQUE_GEO
      end if
      ! gets the number of returned geometric derivatives
      select case (p_geom_type)
      case (REDUNDANT_GEO)
        call GeomTreeGetNumAtoms(geom_tree=geom_tree, num_atoms=num_return_geo)
        num_return_geo = (3*num_return_geo)**order_geo
      case default
        call GeomTreeGetNumGeo(geom_tree=geom_tree, num_unique_geo=num_return_geo)
      end select
    ! no total geometric derivatives
    else
      p_geom_type = UNIQUE_GEO
      num_return_geo = 1
    end if
    ! checks the sizes of returned arguments
    if (present(val_ints)) then
      if (num_prop*num_return_geo>size(val_ints)) then
        stop "Gen1IntAPIPropGetIntExpt>> inconsistent size of val_ints!"
      end if
    end if
    if (present(ao_dens) .and. present(val_expt)) then
      if (size(ao_dens)/=size(val_expt,2)) then
        stop "Gen1IntAPIPropGetIntExpt>> inconsistent number of AO density matrices!"
      end if
      if (num_prop*num_return_geo/=size(val_expt,1)) then
        stop "Gen1IntAPIPropGetIntExpt>> inconsistent size of val_expt!"
      end if
    end if
    ! dumps information to check
    if (level_print>=5) then
      if (same_braket) then
        write(io_viewer,100) "same AO item on bra and ket centers"
      else
        write(io_viewer,100) "different AO items on bra and ket centers"
      end if
      if (do_screen) then
        write(io_viewer,100) "do screen"
      else
        write(io_viewer,100) "no screen"
      end if
    end if
    if (level_print>=10) then
      ! the information of one-electron property integrals
      call Gen1IntAPIPropView(one_prop=one_prop, io_viewer=io_viewer)
      if (present(geom_tree)) then
        write(io_viewer,100) "evalutes total geometric derivatives"
        call Gen1IntAPIGeoTreeView(geom_tree=geom_tree, io_viewer=io_viewer)
        ! type of total geometric derivatives
        select case (p_geom_type)
        case (REDUNDANT_GEO)
          write(io_viewer,100) "redundant total geometric derivatives required"
        case default
          write(io_viewer,100) "unique total geometric derivatives required"
        end select
      end if
      ! arguments related to integral matrices
      if (present(val_ints)) write(io_viewer,100) "integral matrices returned"
      if (present(write_ints)) then
        if (write_ints) write(io_viewer,100) "integral matrices will be written on file"
      end if
      ! arguments related to expectation values
      if (present(ao_dens)) then
        write(io_viewer,100) "number of AO density matrices", size(ao_dens)
        if (present(val_expt)) write(io_viewer,100) "expectation values returned"
        if (present(write_expt)) then
          if (write_expt) write(io_viewer,100) "expectation values will be written on file"
        end if
      end if
    end if
    ! calculates total geometric derivatives
    if (present(geom_tree)) then
      ! gets the index of current path and total number of different paths
      call GeomPathGetIndex(geom_tree=geom_tree, idx_path=idx_path)
      call GeomTreeGetNumPaths(geom_tree=geom_tree, num_paths=num_paths)
      ! calculates the property integrals of current path
      call Gen1IntShellGetIntExpt(ao_item_bra=ao_item_bra, &
                                  ao_item_ket=ao_item_ket, &
                                  same_braket=same_braket, &
                                  do_screen=do_screen,     &
                                  one_prop=one_prop,       &
                                  geom_tree=geom_tree,     &
                                  geom_type=p_geom_type,   &
                                  val_ints=val_ints,       &
                                  write_ints=write_ints,   &
                                  ao_dens=ao_dens,         &
                                  val_expt=val_expt,       &
                                  write_expt=write_expt,   &
                                  io_viewer=io_viewer)
      ! loops over other paths
      do ipath = idx_path+1, num_paths
        ! generates the differentiated centers and their orders
        call GeomTreeSearch(geom_tree=geom_tree)
        ! dumps the information of current path
        if (level_print>=20) &
          call Gen1IntAPIGeoTreeView(geom_tree=geom_tree, io_viewer=io_viewer)
        ! calculates the property integrals of current path
        call Gen1IntShellGetIntExpt(ao_item_bra=ao_item_bra, &
                                    ao_item_ket=ao_item_ket, &
                                    same_braket=same_braket, &
                                    do_screen=do_screen,     &
                                    one_prop=one_prop,       &
                                    geom_tree=geom_tree,     &
                                    geom_type=p_geom_type,   &
                                    val_ints=val_ints,       &
                                    write_ints=write_ints,   &
                                    ao_dens=ao_dens,         &
                                    val_expt=val_expt,       &
                                    write_expt=write_expt,   &
                                    io_viewer=io_viewer)
      end do
    ! no total geometric derivatives
    else
      call Gen1IntShellGetIntExpt(ao_item_bra=ao_item_bra, &
                                  ao_item_ket=ao_item_ket, &
                                  same_braket=same_braket, &
                                  do_screen=do_screen,     &
                                  one_prop=one_prop,       &
                                  geom_tree=geom_tree,     &
                                  geom_type=p_geom_type,   &
                                  val_ints=val_ints,       &
                                  write_ints=write_ints,   &
                                  ao_dens=ao_dens,         &
                                  val_expt=val_expt,       &
                                  write_expt=write_expt,   &
                                  io_viewer=io_viewer)
    end if
100 format("Gen1IntAPIPropGetIntExpt>> ",A,I4)
  end subroutine Gen1IntAPIPropGetIntExpt

  !> \brief creates N-ary tree for total geometric derivatives
  !> \author Bin Gao
  !> \date 2012-05-23
  !> \param num_atoms is the number of atoms
  !> \param max_num_cent is the maximum number of differentiated centers for total
  !>        geometric derivatives
  !> \param order_geo_total is the order of total geometric derivatives
  !> \param idx_geo_atoms contains the indices of the selected atoms as the differentiated centers
  !> \return geom_tree is the N-ary tree for total geometric derivatives
  subroutine Gen1IntAPIGeoTreeCreate(num_atoms, max_num_cent, order_geo_total, &
                                     idx_geo_atoms, geom_tree)
    integer, intent(in) :: num_atoms
    integer, intent(in) :: max_num_cent
    integer, intent(in) :: order_geo_total
    integer, optional, intent(in) :: idx_geo_atoms(num_atoms)
    type(geom_tree_t), intent(inout) :: geom_tree
    integer ierr  !error information
    if (order_geo_total>=0) then
      if (present(idx_geo_atoms)) then
        call GeomTreeCreate(num_atoms=num_atoms,       &
                            order_geo=order_geo_total, &
                            max_num_cent=max_num_cent, &
                            geom_tree=geom_tree,       &
                            info_geom=ierr)
        if (ierr/=0) then
          stop "Gen1IntAPIGeoTreeCreate>> error occurred when calling GeomTreeCreate!"
        end if
        call GeomTreeSetAtoms(num_atoms=num_atoms,     &
                              idx_atoms=idx_geo_atoms, &
                              geom_tree=geom_tree,     &
                              info_geom=ierr)
        if (ierr/=0) then
          stop "Gen1IntAPIGeoTreeCreate>> error occurred when calling GeomTreeSetAtoms!"
        end if
      else
        call GeomTreeCreate(num_atoms=num_atoms,       &
                            order_geo=order_geo_total, &
                            max_num_cent=max_num_cent, &
                            geom_tree=geom_tree,       &
                            info_geom=ierr)
        if (ierr/=0) then
          stop "Gen1IntAPIGeoTreeCreate>> error occurred when calling GeomTreeCreate!"
        end if
      end if
    else
      stop "Gen1IntAPIGeoTreeCreate>> negative order of geometric derivatives!"
    end if
  end subroutine Gen1IntAPIGeoTreeCreate

end module gen1int_api
