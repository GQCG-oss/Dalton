!FIXME: the default is that we use the so called "grand canonical" basis 
!FIXME: - a modification of the contraction coefficients, 
!FIXME: for the time being we need to deactivate use of this
!FIXME: "grand canonical" basis by using the keyword .NOGCBASIS
!FIXME: In conclusion use the the following in LSDALTON.INP:
!FIXME: *DENSOPT
!FIXME: .START
!FIXME: H1DIAG
!FIXME: .NOGCBASIS
!FIXME: we are able to use the normal AO basis
!FIXME: how to treat the so called "grand canonical" basis?

!> \brief module of Gen1Int interface in host program
!> \author Bin Gao and T. Kjaergaard
!> \date 2012-05-23
module gen1int_host
  use precision
  use molecule_typeType, only: MOLECULEINFO
  use AO_TypeType, only: AOITEM
  use integral_type, only: INTEGRALINPUT
  use gen1int_api
  use gen1int_shell!, only: REDUNDANT_GEO, UNIQUE_GEO
  use matrix_util, only: VerifyMatrices
  use lsmatrix_operations_dense
  implicit none

!FIXME: going to remove if they could be read from disk
  type(MOLECULEINFO), save, private :: mol_info      !molecule information
  type(AOITEM), save, private :: ao_items            !AO items
  type(INTEGRALINPUT), save, private :: int_setting  !integral settings

!FIXME: going to remove if the required information could be read from disk
  public :: gen1int_host_init
  public :: gen1int_host_finalize

  public :: gen1int_host_one_prop
  public :: gen1int_host_test

  private :: gen1int_host_norm_gto

  public :: gen1int_host_get_overlap
  public :: gen1int_host_get_first_geoderiv_overlap
  public :: gen1int_host_get_first_geoderiv_overlap_expval
  public :: gen1int_host_get_second_geoderiv_overlap_expval
  public :: gen1int_host_get_second_geoderiv_overlap
  public :: gen1int_host_get_h1
  public :: gen1int_host_get_first_geoderiv_h1
  public :: gen1int_host_get_first_geoderiv_h1_expval
  public :: gen1int_host_get_second_geoderiv_h1_expval

  contains

  !> \brief initializes the information needed for Gen1Int interface
  !> \author Bin Gao and T. Kjaergaard
  !> \date 2012-05-23
  !> \param ls_setting is the settings from linsca
  !> \param io_viewer is the logical unit number of the viewer
  subroutine gen1int_host_init(ls_setting, io_viewer)
    use TYPEDEFTYPE, only: LSSETTING
    use molecule_type, only: copy_molecule
    use BUILDAOBATCH, only: BUILD_AO
    use Integralinfo, only: init_integral_INPUT
    type(LSSETTING), intent(inout) :: ls_setting
    integer, intent(in) :: io_viewer
    !
    logical :: nofamily
    ! copies the molecule information
    call copy_molecule(oldMolecule=ls_setting%molecule(1)%p, &
                       newMolecule=mol_info,                 &
                       lupri=io_viewer)
    ! sets AO items
    !deactivate use of family basis set
    nofamily = ls_setting%SCHEME%nofamily
    ls_setting%SCHEME%nofamily = .TRUE.
    call BUILD_AO(LUPRI=io_viewer,                         &
                  SCHEME=ls_setting%SCHEME,                &
                  IPRINT=ls_setting%SCHEME%AOPRINT,        &
                  MOLECULE=mol_info,                       &
                  BASISINFO=ls_setting%Basis(1)%p%REGULAR, &
                  AO=ao_items,                             &
                  UNCONTRACTED=.false.,                    &
                  INTNRM=.false.,                          &
                  NORMA=.false.,                           &
                  EXTEND=.false.)

    ! normalizes the basis sets
    call gen1int_host_norm_gto(ao_items=ao_items, io_viewer=io_viewer)
    ! sets integral setting
    call init_integral_INPUT(INTINPUT=int_setting, &
                             SETTING=ls_setting)
    !restore settings
    ls_setting%SCHEME%nofamily = nofamily
  end subroutine gen1int_host_init

  !> \brief terminates Gen1Int interface
  !> \author Bin Gao and T. Kjaergaard
  !> \date 2012-05-23
  !> \param io_viewer is the logical unit number of the viewer
  subroutine gen1int_host_finalize(io_viewer)
    use molecule_type, only: free_Moleculeinfo
    use AO_Type, only: FREE_AOITEM
    integer, intent(in) :: io_viewer
    integer :: iang,ibtch
    call free_Moleculeinfo(Molecule=mol_info)

!    do ibtch = 1, ao_items%nbatches
!     do iang = 1, ao_items%BATCH(ibtch)%nAngmom
!      call lsmat_dense_free(ao_items%BATCH(ibtch)%pCC(iang)%p)
!      deallocate(ao_items%BATCH(ibtch)%pCC(iang)%p)
!      nullify(ao_items%BATCH(ibtch)%pCC(iang)%p)
!     enddo
!    enddo
    call FREE_AOITEM(LUPRI=io_viewer, AO=ao_items)
!FIXME: do we need to clean \var(int_setting)?
  end subroutine gen1int_host_finalize

  !> \brief evaluates the one-electron property integrals and/or expectation values
  !> \author Bin Gao
  !> \date 2012-05-24
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
  !> \param max_num_cent is the maximum number of differentiated centers for total
  !>        geometric derivatives
  !> \param order_geo_total is the order of total geometric derivatives
  !> \param idx_geo_atoms contains the indices of the selected atoms which might be chosen as the
  !>        differentiated centers
  !> \param geom_type is the type of returned total geometric derivatives, UNIQUE_GEO is for unique total
  !>        geometric derivatives, and REDUNDANT_GEO for redundant total geometric derivatives;
  !>        for instance, (xx,xy,yy,xz,yz,zz) are unique while (xx,yx,zx,xy,yy,zy,xz,yz,zz) are
  !>        redundant, note that the "triangular" total geometric derivatives could be obtained
  !>        from the unique total geometric derivatives by giving \var(max_num_cent)=\var(order_geo_total)
  !> \param add_sr is for scalar-relativistic (SR) correction, not implemented
  !> \param add_so is for spin-orbit (SO) correction, not implemented
  !> \param add_london transforms the operator by the LAO type gauge-including projector, not implemented
  !> \param write_ints indicates if writing integral matrices on file
  !> \param ao_dens contains the AO density matrices
  !> \param write_expt indicates if writing expectation values on file
  !> \param io_viewer is the logical unit number of the viewer
  !> \param level_print is the level of print
  !> \return val_ints contains the property integral matrices
  !> \return val_expt contains the expectation values
  !> \note the integrals and expectation values will be in the order of \var(order_mom),
  !>       \var(order_mag_bra), ..., \var(order_geo_total), and each of them is arranged
  !>       in the order of (xx,xy,yy,xz,yz,zz) or (xx,yx,zx,xy,yy,zy,xz,yz,zz), see
  !>       Gen1Int library manual, for instance Section 2.2.
  !>       \var(val_ints) should be created by users before calculations
  !>       \var(val_expt) should be zero by users before calculations
  subroutine gen1int_host_one_prop(gto_type, prop_name, order_mom, &
                                   order_mag_bra, order_mag_ket,   &
                                   order_mag_total,                &
                                   order_ram_bra, order_ram_ket,   &
                                   order_ram_total,                &
                                   order_geo_bra, order_geo_ket,   &
                                   max_num_cent, order_geo_total,  &
                                   idx_geo_atoms, geom_type,       &
                                   add_sr, add_so, add_london,     &
                                   val_ints, write_ints,           &
                                   ao_dens, val_expt, write_expt,  &
                                   io_viewer, level_print)
    use Matrix_module, only: Matrix
    use molecule_type, only: MOLECULEINFO, &
                             molecule_get_num_atoms
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
    integer, optional, intent(in) :: max_num_cent
    integer, optional, intent(in) :: order_geo_total
    integer, optional, intent(in) :: idx_geo_atoms(:)
    integer, optional, intent(in) :: geom_type
    logical, optional, intent(in) :: add_sr
    logical, optional, intent(in) :: add_so
    logical, optional, intent(in) :: add_london
    type(Matrix), optional, intent(inout) :: val_ints(:)
    logical, optional, intent(in) :: write_ints
    type(Matrix), optional, intent(in) :: ao_dens(:)
    real(realk), optional, intent(inout) :: val_expt(:,:)
    logical, optional, intent(in) :: write_expt
    integer, intent(in) :: io_viewer
    integer, intent(in) :: level_print
    logical do_geom              !if calculating total geometric derivatives
    integer num_atoms            !number of atoms
    type(one_prop_t) one_prop    !operator of property integrals with non-zero components
    type(geom_tree_t) geom_tree  !N-ary tree for total geometric derivatives
    integer ierr                 !error information
    real(realk) start_time       !start time
    ! gets the start time
    call xtimer_set(start_time)
    ! creates the operator of property integrals and N-ary tree for total geometric derivatives
    call Gen1IntAPIPropCreate(mol_info=mol_info,               &
                              gto_type=gto_type,               &
                              prop_name=prop_name,             &
                              order_mom=order_mom,             &
                              order_mag_bra=order_mag_bra,     &
                              order_mag_ket=order_mag_ket,     &
                              order_mag_total=order_mag_total, &
                              order_ram_bra=order_ram_bra,     &
                              order_ram_ket=order_ram_ket,     &
                              order_ram_total=order_ram_total, &
                              order_geo_bra=order_geo_bra,     &
                              order_geo_ket=order_geo_ket,     &
                              add_sr=add_sr,                   &
                              add_so=add_so,                   &
                              add_london=add_london,           &
                              one_prop=one_prop,               &
                              io_viewer=io_viewer)
    do_geom = present(max_num_cent) .and. present(order_geo_total)
    if (do_geom) then
      if (present(idx_geo_atoms)) then
        call Gen1IntAPIGeoTreeCreate(num_atoms=size(idx_geo_atoms),   &
                                     max_num_cent=max_num_cent,       &
                                     order_geo_total=order_geo_total, &
                                     idx_geo_atoms=idx_geo_atoms,     &
                                     geom_tree=geom_tree)
      else
        call molecule_get_num_atoms(mol_info=mol_info, num_atoms=num_atoms)
        call Gen1IntAPIGeoTreeCreate(num_atoms=num_atoms,             &
                                     max_num_cent=max_num_cent,       &
                                     order_geo_total=order_geo_total, &
                                     geom_tree=geom_tree)
      end if
    end if
    ! performs calculations
    if (do_geom) then
      call Gen1IntAPIPropGetIntExpt(ao_item_bra=ao_items,            &
                                    ao_item_ket=ao_items,            &
                                    same_braket=.true.,              &
                                    do_screen=int_setting%OD_SCREEN, &
                                    one_prop=one_prop,               &
                                    geom_tree=geom_tree,             &
                                    geom_type=geom_type,             &
                                    val_ints=val_ints,               &
                                    write_ints=write_ints,           &
                                    ao_dens=ao_dens,                 &
                                    val_expt=val_expt,               &
                                    write_expt=write_expt,           &
                                    io_viewer=io_viewer,             &
                                    level_print=level_print)
    else
      call Gen1IntAPIPropGetIntExpt(ao_item_bra=ao_items,            &
                                    ao_item_ket=ao_items,            &
                                    same_braket=.true.,              &
                                    do_screen=int_setting%OD_SCREEN, &
                                    one_prop=one_prop,               &
                                    geom_type=geom_type,             &
                                    val_ints=val_ints,               &
                                    write_ints=write_ints,           &
                                    ao_dens=ao_dens,                 &
                                    val_expt=val_expt,               &
                                    write_expt=write_expt,           &
                                    io_viewer=io_viewer,             &
                                    level_print=level_print)
    end if
    ! free spaces
    call Gen1IntAPIPropDestroy(one_prop=one_prop)
    if (do_geom) call Gen1IntAPIGeoTreeDestroy(geom_tree=geom_tree)
    ! prints the CPU elapsed time
    call xtimer_view(start_time, trim(prop_name)//"@gen1int_host_one_prop", io_viewer)
  end subroutine gen1int_host_one_prop

  !> \brief test suite of Gen1Int interface
  !> \author Bin Gao and T. Kjaergaard
  !> \date 2011-10-09
  !> \param ls_setting is the integral evalualtion settings
  !> \param D contains the AO density matrix
  !> \param io_viewer logical unit number of the standard output
  subroutine gen1int_host_test(ls_setting, D, io_viewer)
    use TYPEDEFTYPE, only: LSSETTING
    use TYPEDEF, only: getNbasis
    use Matrix_module, only: Matrix
    use matrix_operations, only: mat_init,  &
                                 mat_print, &
                                 mat_add,   &
                                 mat_free,  &
                                 mat_trab
    use LSparameters, only: AORegular, &
                                 Contractedinttype
    use IntegralInterfaceMOD, only: II_get_prop,ii_get_kinetic,ii_get_nucel_mat,ii_get_overlap
    type(LSSETTING), intent(inout) :: ls_setting
    type(Matrix), intent(in) :: D
    integer, intent(in) :: io_viewer
    integer num_ao                             !number of atomic orbitals
    integer level_print                        !level of print
    type(Matrix), allocatable :: val_ints(:)   !one-electron property integral matrices from Gen1Int
    real(realk), allocatable :: val_expt(:,:)  !expectation values from Gen1Int
    type(Matrix), pointer :: host_ints(:)      !one-electron property integral matrices from host program
!    type(Matrix) diff_matrix                   !difference matrix
    integer iopt                               !incremental recorder over operators
    integer ierr                               !error information
    real(realk) :: THRESH
    THRESH=1.0E-10
    ! initializes the information needed for Gen1Int interface
    call gen1int_host_init(ls_setting=ls_setting, io_viewer=io_viewer)
    ! gets the number of atomic orbitals
    num_ao = getNbasis(AORegular, Contractedinttype, mol_info, io_viewer)
    ! gets the level of print
    level_print = ls_setting%SCHEME%INTPRINT
    if (level_print>=10) &
      write(io_viewer,100) "number of atomic orbitals", num_ao
    !====================================================================
    ! (1) tests Cartesian multipole moments
    !====================================================================
    allocate(val_ints(3), stat=ierr)
    if (ierr/=0) &
      call lsquit("gen1int_host_test>> failed to allocate val_ints!", io_viewer)
    do iopt = 1, 3
      call mat_init(a=val_ints(iopt), nrow=num_ao, ncol=num_ao, complex=.false.)
    end do
    allocate(val_expt(3,1), stat=ierr)
    if (ierr/=0) &
      call lsquit("gen1int_host_test>> failed to allocate val_expt!", io_viewer)
    val_expt = 0.0
    ! gets results from Gen1Int
    call gen1int_host_one_prop(prop_name=INT_CART_MULTIPOLE, &
                               order_mom=1,                  &
                               val_ints=val_ints,            &
                               ao_dens=(/D/),                &
                               val_expt=val_expt,            &
                               io_viewer=io_viewer,          &
                               level_print=level_print)
    ! gets results from lsdalton
    allocate(host_ints(3), stat=ierr)
    if (ierr/=0) &
      call lsquit("gen1int_host_test>> failed to allocate host_ints!", io_viewer)
    do iopt = 1, 3
      call mat_init(a=host_ints(iopt), nrow=num_ao, ncol=num_ao, complex=.false.)
    end do
    call II_get_prop(io_viewer, io_viewer, ls_setting, host_ints, 3, 'DIPLEN ')
 !   call II_get_integral(io_viewer, io_viewer, ls_setting, host_ints, 3, "DIPLEN ")
    ! checks the results from Gen1Int
!    call mat_init(a=diff_matrix, nrow=num_ao, ncol=num_ao, complex=.false.)
    call VerifyMatrices(MAT1=host_ints(1),MAT2=val_ints(1),STRING='DIPLENX',THR=THRESH,lu=io_viewer)
    call VerifyMatrices(MAT1=host_ints(2),MAT2=val_ints(2),STRING='DIPLENY',THR=THRESH,lu=io_viewer)
    call VerifyMatrices(MAT1=host_ints(3),MAT2=val_ints(3),STRING='DIPLENZ',THR=THRESH,lu=io_viewer)
!     write(io_viewer,100) "dipole length matrix from linsca", iopt
!     call mat_print(a=host_ints(iopt),i_row1=1,i_rown=num_ao,j_col1=1,j_coln=num_ao,lu=io_viewer)
!     write(io_viewer,100) "dipole length matrix from Gen1Int", iopt
!     call mat_print(a=val_ints(iopt),i_row1=1,i_rown=num_ao,j_col1=1,j_coln=num_ao,lu=io_viewer)
!     calculcates the difference matrix
!     call mat_add(alpha=1.0_realk,a=host_ints(iopt),beta=-1.0_realk,b=val_ints(iopt),c=diff_matrix)
!     write(io_viewer,110) "difference matrix", mat_trab(diff_matrix, diff_matrix)
!     call mat_print(a=diff_matrix,i_row1=1,i_rown=num_ao,j_col1=1,j_coln=num_ao,lu=io_viewer)
    do iopt = 1, 3
       call mat_free(a=host_ints(iopt))
       write(io_viewer,110) "DIPLEN expectation value from mat_trab", &
            mat_trab(D,val_ints(iopt))
       write(io_viewer,110) "DIPLEN expectation value from Gen1Int ", &
            val_expt(iopt,1)
       IF(ABS(val_expt(iopt,1)-mat_trab(D,val_ints(iopt))).GT.THRESH)THEN
          call lsquit("gen1int_host_test>> failed diplen expectation val fail!", io_viewer)
       ENDIF
       call mat_free(a=val_ints(iopt))
    end do
    deallocate(val_ints)
    deallocate(val_expt)
    deallocate(host_ints)
!    call mat_free(a=diff_matrix)
    !====================================================================
    ! (2) tests kinetic energy integrals
    !====================================================================

    allocate(val_ints(1), stat=ierr)
    if (ierr/=0) &
      call lsquit("gen1int_host_test>> failed to allocate val_ints!", io_viewer)
    call mat_init(a=val_ints(1), nrow=num_ao, ncol=num_ao, complex=.false.)
    allocate(val_expt(1,1), stat=ierr)
    if (ierr/=0) &
      call lsquit("gen1int_host_test>> failed to allocate val_expt!", io_viewer)
    val_expt = 0.0
    ! gets results from Gen1Int
    call gen1int_host_one_prop(prop_name=INT_KIN_ENERGY, &
                               val_ints=val_ints,        &
                               ao_dens=(/D/),            &
                               val_expt=val_expt,        &
                               io_viewer=io_viewer,      &
                               level_print=level_print)
    ! gets results from linsca
    allocate(host_ints(1), stat=ierr)
    if (ierr/=0) &
      call lsquit("gen1int_host_test>> failed to allocate host_ints!", io_viewer)
    call mat_init(a=host_ints(1), nrow=num_ao, ncol=num_ao, complex=.false.)
    call II_get_kinetic(io_viewer, io_viewer, ls_setting, host_ints(1))
    call VerifyMatrices(MAT1=host_ints(1),MAT2=val_ints(1),STRING='Kinetic',THR=THRESH,lu=io_viewer)

!    write(io_viewer,100) "kinetic energy matrix from linsca"
!    call mat_print(a=host_ints(1),          &
!                   i_row1=1, i_rown=num_ao, &
!                   j_col1=1, j_coln=num_ao, &
!                   lu=io_viewer)
!    write(io_viewer,100) "kinetic energy matrix from Gen1Int"
!    call mat_print(a=val_ints(1),           &
!                   i_row1=1, i_rown=num_ao, & 
!                   j_col1=1, j_coln=num_ao, & 
!                   lu=io_viewer)
!    ! checks the results from Gen1Int
!    call mat_init(a=diff_matrix, nrow=num_ao, ncol=num_ao, complex=.false.)
!    call mat_add(alpha=1.0_realk, &
!                 a=host_ints(1),  &
!                 beta=-1.0_realk, &
!                 b=val_ints(1),   &
!                 c=diff_matrix)
!    write(io_viewer,110) "difference matrix", &
!                         mat_trab(diff_matrix, diff_matrix)
!    call mat_print(a=diff_matrix,           &
!                   i_row1=1, i_rown=num_ao, &
!                   j_col1=1, j_coln=num_ao, &
!                   lu=io_viewer)
    call mat_free(a=host_ints(1))
    write(io_viewer,110) "Kinetic expectation value from mat_trab", &
                         mat_trab(D,val_ints(1))
    write(io_viewer,110) "Kinetic expectation value from Gen1Int ", &
                         val_expt(1,1)
    IF(ABS(val_expt(1,1)-mat_trab(D,val_ints(1))).GT.THRESH)THEN
       call lsquit("gen1int_host_test>> failed Kinetic expectation val fail!", io_viewer)
    ENDIF
    call mat_free(a=val_ints(1))
    deallocate(val_ints)
    deallocate(val_expt)
    deallocate(host_ints)
!    call mat_free(a=diff_matrix)
    !====================================================================
    ! (3) tests overlap integrals
    !====================================================================
    allocate(val_ints(1), stat=ierr)
    if (ierr/=0) &
      call lsquit("gen1int_host_test>> failed to allocate val_ints!", io_viewer)
    call mat_init(a=val_ints(1), nrow=num_ao, ncol=num_ao, complex=.false.)
    allocate(val_expt(1,1), stat=ierr)
    if (ierr/=0) &
      call lsquit("gen1int_host_test>> failed to allocate val_expt!", io_viewer)
    val_expt = 0.0
    ! gets results from Gen1Int
    call gen1int_host_one_prop(prop_name=INT_OVERLAP, &
                               val_ints=val_ints,     &
                               ao_dens=(/D/),         &
                               val_expt=val_expt,     &
                               io_viewer=io_viewer,   &
                               level_print=level_print)
    ! gets results from linsca
    allocate(host_ints(1), stat=ierr)
    if (ierr/=0) &
      call lsquit("gen1int_host_test>> failed to allocate host_ints!", io_viewer)
    call mat_init(a=host_ints(1), nrow=num_ao, ncol=num_ao, complex=.false.)
    call II_get_overlap(io_viewer, io_viewer, ls_setting, host_ints(1))
    call VerifyMatrices(MAT1=host_ints(1),MAT2=val_ints(1),STRING='Overlap',THR=THRESH,lu=io_viewer)
!    write(io_viewer,100) "overlap matrix from linsca"
!    call mat_print(a=host_ints(1),          &
!                   i_row1=1, i_rown=num_ao, &
!                   j_col1=1, j_coln=num_ao, &
!                   lu=io_viewer)
!    write(io_viewer,100) "overlap matrix from Gen1Int"
!    call mat_print(a=val_ints(1),           &
!                   i_row1=1, i_rown=num_ao, &
!                   j_col1=1, j_coln=num_ao, &
!                   lu=io_viewer)
!    ! checks the results from Gen1Int
!    call mat_init(a=diff_matrix, nrow=num_ao, ncol=num_ao, complex=.false.)
!    call mat_add(alpha=1.0_realk, &
!                 a=host_ints(1),  &
!                 beta=-1.0_realk, &
!                 b=val_ints(1),   &
!                 c=diff_matrix)
!    write(io_viewer,110) "difference matrix", &
!                         mat_trab(diff_matrix, diff_matrix)
!    call mat_print(a=diff_matrix,           &
!                   i_row1=1, i_rown=num_ao, &
!                   j_col1=1, j_coln=num_ao, &
!                   lu=io_viewer)
    call mat_free(a=host_ints(1))
    write(io_viewer,110) "Overlap expectation value from mat_trab", &
                         mat_trab(D,val_ints(1))
    write(io_viewer,110) "Overlap expectation value from Gen1Int ", &
                         val_expt(1,1)
    IF(ABS(val_expt(1,1)-mat_trab(D,val_ints(1))).GT.THRESH)THEN
       call lsquit("gen1int_host_test>> failed Overlap expectation val fail!", io_viewer)
    ENDIF
    call mat_free(a=val_ints(1))
    deallocate(val_ints)
    deallocate(val_expt)
    deallocate(host_ints)
!    call mat_free(a=diff_matrix)
    !====================================================================
    ! (4) tests one-electron potential energy integrals
    !====================================================================
    allocate(val_ints(1), stat=ierr)
    if (ierr/=0) &
      call lsquit("gen1int_host_test>> failed to allocate val_ints!", io_viewer)
    call mat_init(a=val_ints(1), nrow=num_ao, ncol=num_ao, complex=.false.)
    allocate(val_expt(1,1), stat=ierr)
    if (ierr/=0) &
      call lsquit("gen1int_host_test>> failed to allocate val_expt!", io_viewer)
    val_expt = 0.0
    ! gets results from Gen1Int
    call gen1int_host_one_prop(prop_name=INT_POT_ENERGY, &
                               val_ints=val_ints,        &
                               ao_dens=(/D/),            &
                               val_expt=val_expt,        &
                               io_viewer=io_viewer,      &
                               level_print=level_print)
    ! gets results from linsca
    allocate(host_ints(1), stat=ierr)
    if (ierr/=0) &
      call lsquit("gen1int_host_test>> failed to allocate host_ints!", io_viewer)
    call mat_init(a=host_ints(1), nrow=num_ao, ncol=num_ao, complex=.false.)
    call II_get_nucel_mat(io_viewer, io_viewer, ls_setting, host_ints(1))
    call VerifyMatrices(MAT1=host_ints(1),MAT2=val_ints(1),STRING='Nucel',THR=THRESH,lu=io_viewer)
!    write(io_viewer,100) "potential energy matrix from linsca"
!    call mat_print(a=host_ints(1),          &
!                   i_row1=1, i_rown=num_ao, &
!                   j_col1=1, j_coln=num_ao, &
!                   lu=io_viewer)
!    write(io_viewer,100) "potential energy matrix from Gen1Int"
!    call mat_print(a=val_ints(1),           &
!                   i_row1=1, i_rown=num_ao, &
!                   j_col1=1, j_coln=num_ao, &
!                   lu=io_viewer)
!    ! checks the results from Gen1Int
!    call mat_init(a=diff_matrix, nrow=num_ao, ncol=num_ao, complex=.false.)
!    call mat_add(alpha=1.0_realk, &
!                 a=host_ints(1),  &
!                 beta=-1.0_realk, &
!                 b=val_ints(1),   &
!                 c=diff_matrix)
!    write(io_viewer,110) "difference matrix", &
!                         mat_trab(diff_matrix, diff_matrix)
!    call mat_print(a=diff_matrix,           &
!                   i_row1=1, i_rown=num_ao, &
!                   j_col1=1, j_coln=num_ao, &
!                   lu=io_viewer)
    call mat_free(a=host_ints(1))
    write(io_viewer,110) "nucel expectation value from mat_trab", &
                         mat_trab(D,val_ints(1))
    write(io_viewer,110) "nucel expectation value from Gen1Int ", &
                         val_expt(1,1)
    IF(ABS(val_expt(1,1)-mat_trab(D,val_ints(1))).GT.THRESH)THEN
       call lsquit("gen1int_host_test>> failed nucel expectation val fail!", io_viewer)
    ENDIF
    call mat_free(a=val_ints(1))
    deallocate(val_ints)
    deallocate(val_expt)
    deallocate(host_ints)
!    call mat_free(a=diff_matrix)
    ! terminates Gen1Int interface
    call gen1int_host_finalize(io_viewer=io_viewer)
100 format("gen1int_host_test>> ",A,I12)
110 format("gen1int_host_test>> ",A,F16.10)
  end subroutine gen1int_host_test

  !> \brief normalizes the contraction coefficients of GTOs
  !> \author Bin Gao and T. Kjaergaard
  !> \date 2011-10-08
  !> \param ao_items contains the information of basis sets
  !> \param io_viewer logical unit number of the standard output file
  subroutine gen1int_host_norm_gto(ao_items, io_viewer)
    use lsmatrix_operations_dense, only: lsmat_dense_init
    type(AOITEM), intent(inout) :: ao_items
    integer, intent(in) :: io_viewer
    integer ibtch                                !incremental recorder over \var(nbatches)
    integer iang                                 !incremental recorder over angular momentum?
    integer num_prim                             !number of primitive GTOs
    integer ang_num                              !angular number
    integer num_contr                            !number of contractions
    real(realk), allocatable :: contr_coef(:,:)  !temporary contraction coefficients for normalization
    integer icontr, iprim                        !incremental recorders over contractions and primitives
    integer ierr                                 !error information
    integer :: iCC,i
    !sanity check
    IF(AO_items%nCC.NE.AO_items%nEXP)THEN
       print*,'AO_items%nCC',AO_items%nCC
       print*,'AO_items%nEXP',AO_items%nEXP
       CALL LSQUIT('AO_items%nCC.NE.AO_items%nEXP',-1)
    ENDIF
    
    do I=1,AO_items%nbatches
       IF(AO_items%BATCH(I)%nANGMOM.GT.1)call lsquit('Error with NOFAMILY',-1)
    enddo
    
    do iCC = 1, ao_items%nCC
       num_prim = ao_items%CC(iCC)%nrow
       num_contr = ao_items%CC(iCC)%ncol
       ang_num = ao_items%angmom(iCC)
       allocate(contr_coef(num_contr,num_prim), stat=ierr)
       if (ierr/=0) then
          call lsquit("gen1int_host_norm_gto>> failed to allocate contr_coef!", &
               io_viewer)
       endif
       do icontr = 1, num_contr
          do iprim = 1, num_prim
             contr_coef(icontr,iprim) &
                  = ao_items%CC(iCC)%elms(iprim+(icontr-1)*num_prim)
          end do
       end do
       ! normalizes SGTOs
       if (ao_items%BATCH(1)%spherical) then
          call norm_contr_sgto(ang_num, num_prim, &
               ao_items%Exponents(iCC)%elms(1:num_prim), &
               num_contr, contr_coef)
          ! normalizes CGTOs
       else
          call norm_contr_cgto(ang_num, num_prim, &
               ao_items%Exponents(iCC)%elms(1:num_prim), &
               num_contr, contr_coef) 
       end if
       ao_items%CC(iCC)%nrow = num_contr
       ao_items%CC(iCC)%ncol = num_prim
       ! gets the normalized contraction coefficients back
       do iprim = 1, num_prim
          do icontr = 1, num_contr
             ao_items%CC(iCC)%elms(icontr+(iprim-1)*num_contr) &
                  = contr_coef(icontr,iprim)
          end do
       end do
       deallocate(contr_coef)
    enddo
  end subroutine gen1int_host_norm_gto

  !> \brief Calc Overlap using Gen1Int
  !> \author Bin Gao and T. Kjaergaard
  !> \date 2011-10-09
  !> \param ls_setting is the integral evalualtion settings
  !> \param S contains the AO overlap matrix
  !> \param io_viewer logical unit number of the standard output
  subroutine gen1int_host_get_overlap(ls_setting, S, io_viewer)
    use TYPEDEFTYPE, only: LSSETTING
    use Matrix_module, only: Matrix
    type(LSSETTING), intent(inout) :: ls_setting
    type(Matrix), intent(inout) :: S(1)
    integer, intent(in) :: io_viewer         
    integer num_ao,level_print                 !number of atomic orbitals and level of print
    integer ierr !error information
    ! initializes the information needed for Gen1Int interface
    call gen1int_host_init(ls_setting=ls_setting, io_viewer=io_viewer)
    ! gets the number of atomic orbitals
    num_ao = S(1)%nrow
    ! gets the level of print
    level_print = ls_setting%SCHEME%INTPRINT
    ! gets results from Gen1Int
    call gen1int_host_one_prop(prop_name=INT_OVERLAP, val_ints=S,&
         & io_viewer=io_viewer, level_print=level_print)
    ! terminates Gen1Int interface
    call gen1int_host_finalize(io_viewer=io_viewer)
  end subroutine gen1int_host_get_overlap

  !> \brief Calc Overlap using Gen1Int
  !> \author Bin Gao and T. Kjaergaard
  !> \date 2011-10-09
  !> \param ls_setting is the integral evalualtion settings
  !> \param S contains the AO overlap matrix
  !> \param io_viewer logical unit number of the standard output
  subroutine gen1int_host_get_first_geoderiv_overlap(ls_setting, Sx, natoms,io_viewer)
    use TYPEDEFTYPE, only: LSSETTING
    use Matrix_module, only: Matrix
    type(LSSETTING), intent(inout) :: ls_setting
    type(Matrix), intent(inout) :: Sx(3*natoms)
    integer, intent(in) :: io_viewer,natoms         
    integer num_ao,level_print                 !number of atomic orbitals and level of print
    integer ierr !error information
    ! initializes the information needed for Gen1Int interface
    call gen1int_host_init(ls_setting=ls_setting, io_viewer=io_viewer)
    ! gets the number of atomic orbitals
    num_ao = Sx(1)%nrow
    ! gets the level of print
    level_print = ls_setting%SCHEME%INTPRINT
    ! gets results from Gen1Int
    ! Overlap S_ab = (a|b) 
    ! first geoderiv overlap S^x = (a^x|b)+(a|b^x)
    call gen1int_host_one_prop(prop_name=INT_OVERLAP, &
         & max_num_cent=1, &
         & order_geo_total=1,&
         & val_ints=Sx,io_viewer=io_viewer, level_print=level_print)
    ! terminates Gen1Int interface
    call gen1int_host_finalize(io_viewer=io_viewer)
  end subroutine gen1int_host_get_first_geoderiv_overlap

  !> \brief Calc Overlap using Gen1Int
  !> \author Bin Gao and T. Kjaergaard
  !> \date 2011-10-09
  !> \param ls_setting is the integral evalualtion settings
  !> \param S contains the AO overlap matrix
  !> \param io_viewer logical unit number of the standard output
  subroutine gen1int_host_get_first_geoderiv_overlap_expval(ls_setting, D, grad, natoms,io_viewer)
    use TYPEDEFTYPE, only: LSSETTING
    use Matrix_module, only: Matrix
    type(LSSETTING), intent(inout) :: ls_setting
    type(Matrix), intent(in) :: D
    real(realk), intent(inout) :: grad(3*natoms,1)
    integer, intent(in) :: io_viewer,natoms         
    integer num_ao,level_print                 !number of atomic orbitals and level of print
    integer ierr !error information
    ! initializes the information needed for Gen1Int interface
    call gen1int_host_init(ls_setting=ls_setting, io_viewer=io_viewer)
    ! gets the number of atomic orbitals
    num_ao = D%nrow
    ! gets the level of print
    level_print = ls_setting%SCHEME%INTPRINT
    ! gets results from Gen1Int
    ! Overlap S_ab = (a|b) 
    ! first geoderiv overlap Tr(D,S^x) = Tr(D_ab*(a^x|b)+(a|b^x))
    call gen1int_host_one_prop(prop_name=INT_OVERLAP, &
         & max_num_cent=1, &
         & order_geo_total=1,&
         & ao_dens=(/D/), val_expt=grad, &
         & io_viewer=io_viewer, level_print=level_print)
    ! terminates Gen1Int interface
    call gen1int_host_finalize(io_viewer=io_viewer)
  end subroutine gen1int_host_get_first_geoderiv_overlap_expval

  !> \brief Calc Overlap using Gen1Int
  !> \author Bin Gao and T. Kjaergaard
  !> \date 2011-10-09
  !> \param ls_setting is the integral evalualtion settings
  !> \param S contains the AO overlap matrix
  !> \param io_viewer logical unit number of the standard output
  subroutine gen1int_host_get_second_geoderiv_overlap_expval(&
       & ls_setting, D, hess, natoms,io_viewer)
    use TYPEDEFTYPE, only: LSSETTING
    use Matrix_module, only: Matrix
    type(LSSETTING), intent(inout) :: ls_setting
    type(Matrix), intent(in) :: D
    real(realk), intent(inout) :: hess(9*natoms*natoms,1)
    integer, intent(in) :: io_viewer,natoms         
    integer num_ao,level_print                 !number of atomic orbitals and level of print
    integer ierr !error information
    ! initializes the information needed for Gen1Int interface
    call gen1int_host_init(ls_setting=ls_setting, io_viewer=io_viewer)
    ! gets the number of atomic orbitals
    num_ao = D%nrow
    ! gets the level of print
    level_print = ls_setting%SCHEME%INTPRINT
    ! gets results from Gen1Int
    ! Overlap S_ab = (a|b) 
    ! second geoderiv overlap Tr(D,S^xy) = Tr(D_ab*((a^xy|b)+(a^y|b^x)+(a^x|b^y)+(a|b^xy)))
    !WARNING full set of contributions return because of geom_type = REDUNDANT_GEO,&
    !we should not use this and instead only get the triangular part (3Natoms*(3Natoms+1)/2)
    call gen1int_host_one_prop(prop_name=INT_OVERLAP, &
         & max_num_cent=2, &
         & order_geo_total=2,&
         & geom_type = REDUNDANT_GEO,&
         & ao_dens=(/D/), val_expt=hess, &
         & io_viewer=io_viewer, level_print=level_print)
    ! terminates Gen1Int interface
    call gen1int_host_finalize(io_viewer=io_viewer)
  end subroutine gen1int_host_get_second_geoderiv_overlap_expval

  !> \brief Calc Overlap using Gen1Int
  !> \author Bin Gao and T. Kjaergaard
  !> \date 2011-10-09
  !> \param ls_setting is the integral evalualtion settings
  !> \param S contains the AO overlap matrix
  !> \param io_viewer logical unit number of the standard output
  subroutine gen1int_host_get_second_geoderiv_overlap(ls_setting, Sxy, natoms,io_viewer)
    use TYPEDEFTYPE, only: LSSETTING
    use Matrix_module, only: Matrix
    integer, intent(in) :: io_viewer,natoms
    type(LSSETTING), intent(inout) :: ls_setting
    type(Matrix), intent(inout) :: Sxy(9*natoms*natoms)
    integer num_ao,level_print                 !number of atomic orbitals and level of print
    integer ierr !error information
    ! initializes the information needed for Gen1Int interface
    call gen1int_host_init(ls_setting=ls_setting, io_viewer=io_viewer)
    ! gets the number of atomic orbitals
    num_ao = Sxy(1)%nrow
    ! gets the level of print
    level_print = ls_setting%SCHEME%INTPRINT
    ! gets results from Gen1Int
    ! Overlap S_ab = (a|b) 
    ! first geoderiv overlap S^x = (a^x|b)+(a|b^x)
    ! second geoderiv overlap S^xy = (a^xy|b)+(a^x|b^y)+(a^y|b^x)+(a|b^xy)
    !WARNING full set of contributions return because of geom_type = REDUNDANT_GEO,&
    !we should not use this and instead only get the triangular part (3Natoms*(3Natoms+1)/2)
    call gen1int_host_one_prop(prop_name=INT_OVERLAP, &
         & max_num_cent=2, &
         & order_geo_total=2,&
         & geom_type = REDUNDANT_GEO,&
         & val_ints=Sxy,io_viewer=io_viewer, level_print=level_print)
    ! terminates Gen1Int interface
    call gen1int_host_finalize(io_viewer=io_viewer)
  end subroutine gen1int_host_get_second_geoderiv_overlap

  !> \brief Calc Overlap using Gen1Int
  !> \author Bin Gao and T. Kjaergaard
  !> \date 2011-10-09
  !> \param ls_setting is the integral evalualtion settings
  !> \param S contains the AO overlap matrix
  !> \param io_viewer logical unit number of the standard output
  subroutine gen1int_host_get_h1(ls_setting, h1, io_viewer)
    use TYPEDEFTYPE, only: LSSETTING
    use Matrix_module, only: Matrix
    type(LSSETTING), intent(inout) :: ls_setting
    type(Matrix), intent(inout) :: h1(1)
    integer, intent(in) :: io_viewer         
    integer num_ao,level_print                 !number of atomic orbitals and level of print
    integer ierr !error information
    ! initializes the information needed for Gen1Int interface
    call gen1int_host_init(ls_setting=ls_setting, io_viewer=io_viewer)
    ! gets the number of atomic orbitals
    num_ao = h1(1)%nrow
    ! gets the level of print
    level_print = ls_setting%SCHEME%INTPRINT
    ! gets results from Gen1Int
    call gen1int_host_one_prop(prop_name=INT_ONE_HAMIL, val_ints=h1,&
         & io_viewer=io_viewer, level_print=level_print)
    ! terminates Gen1Int interface
    call gen1int_host_finalize(io_viewer=io_viewer)
  end subroutine gen1int_host_get_h1

  !> \brief Calc Overlap using Gen1Int
  !> \author Bin Gao and T. Kjaergaard
  !> \date 2011-10-09
  !> \param ls_setting is the integral evalualtion settings
  !> \param S contains the AO overlap matrix
  !> \param io_viewer logical unit number of the standard output
  subroutine gen1int_host_get_first_geoderiv_h1(ls_setting, h1x,natoms,io_viewer)
    use TYPEDEFTYPE, only: LSSETTING
    use Matrix_module, only: Matrix
    integer, intent(in) :: io_viewer,natoms
    type(LSSETTING), intent(inout) :: ls_setting
    type(Matrix), intent(inout) :: h1x(3*natoms)
    integer num_ao,level_print    !number of atomic orbitals and level of print
    integer ierr !error information
    ! initializes the information needed for Gen1Int interface
    call gen1int_host_init(ls_setting=ls_setting, io_viewer=io_viewer)
    ! gets the number of atomic orbitals
    num_ao = h1x(1)%nrow
    ! gets the level of print
    level_print = ls_setting%SCHEME%INTPRINT
    ! gets results from Gen1Int
    call gen1int_host_one_prop(prop_name=INT_ONE_HAMIL,&
         & max_num_cent=1, &
         & order_geo_total=1,&
         & val_ints=h1x,io_viewer=io_viewer, level_print=level_print)
    ! terminates Gen1Int interface
    call gen1int_host_finalize(io_viewer=io_viewer)
  end subroutine gen1int_host_get_first_geoderiv_h1

  !> \brief Calc Overlap using Gen1Int
  !> \author Bin Gao and T. Kjaergaard
  !> \date 2011-10-09
  !> \param ls_setting is the integral evalualtion settings
  !> \param S contains the AO overlap matrix
  !> \param io_viewer logical unit number of the standard output
  subroutine gen1int_host_get_first_geoderiv_h1_expval(&
       & ls_setting, D,h1xgrad, natoms,io_viewer)
    use TYPEDEFTYPE, only: LSSETTING
    use Matrix_module, only: Matrix
    integer, intent(in) :: io_viewer,natoms
    type(LSSETTING), intent(inout) :: ls_setting
    type(Matrix), intent(in) :: D
    real(realk) :: h1xgrad(3*natoms,1)
    integer num_ao,level_print                 !number of atomic orbitals and level of print
    integer ierr !error information
    ! initializes the information needed for Gen1Int interface
    call gen1int_host_init(ls_setting=ls_setting, io_viewer=io_viewer)
    ! gets the number of atomic orbitals
    num_ao = D%nrow
    ! gets the level of print
    level_print = ls_setting%SCHEME%INTPRINT
    ! gets results from Gen1Int
    call gen1int_host_one_prop(prop_name=INT_ONE_HAMIL, &
         & max_num_cent=1, &
         & order_geo_total=1,&
         & ao_dens=(/D/), val_expt=h1xgrad, &
         & io_viewer=io_viewer, level_print=level_print)
    ! terminates Gen1Int interface
    call gen1int_host_finalize(io_viewer=io_viewer)
  end subroutine gen1int_host_get_first_geoderiv_h1_expval

  !> \brief Calc Overlap using Gen1Int
  !> \author Bin Gao and T. Kjaergaard
  !> \date 2011-10-09
  !> \param ls_setting is the integral evalualtion settings
  !> \param S contains the AO overlap matrix
  !> \param io_viewer logical unit number of the standard output
  subroutine gen1int_host_get_second_geoderiv_h1_expval(&
       & ls_setting, D,h1xy, natoms,io_viewer)
    use TYPEDEFTYPE, only: LSSETTING
    use Matrix_module, only: Matrix
    integer, intent(in) :: io_viewer,natoms
    type(LSSETTING), intent(inout) :: ls_setting
    type(Matrix), intent(in) :: D
    real(realk) :: h1xy(3*natoms,3*natoms)
    integer num_ao,level_print                 !number of atomic orbitals and level of print
    integer ierr !error information
    ! initializes the information needed for Gen1Int interface
    call gen1int_host_init(ls_setting=ls_setting, io_viewer=io_viewer)
    ! gets the number of atomic orbitals
    num_ao = D%nrow
    ! gets the level of print
    level_print = ls_setting%SCHEME%INTPRINT
    ! gets results from Gen1Int
    call gen1int_host_one_prop(prop_name=INT_ONE_HAMIL, &
         & max_num_cent=2, &
         & order_geo_total=2,&
         & ao_dens=(/D/), val_expt=h1xy, &
         & io_viewer=io_viewer, level_print=level_print)
    ! terminates Gen1Int interface
    call gen1int_host_finalize(io_viewer=io_viewer)
  end subroutine gen1int_host_get_second_geoderiv_h1_expval

end module gen1int_host
