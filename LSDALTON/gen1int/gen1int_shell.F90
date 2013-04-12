!> \brief calculates integrals/expectation values between AO sub-shells
!> \author Bin Gao and T. Kjaergaard
!> \date 2012-05-23
module gen1int_shell

  ! Fortran 90 module of Gen1Int library
  use gen1int

  use precision
  use OD_Type, only: getODscreening
  use Matrix_module, only: Matrix
!FIXME: implements subroutines to put the values into matrix or get the values from AO density matrices
  !use matrix_operations

  implicit none

  ! types of returned total geometric derivatives, UNIQUE_GEO is for unique total
  ! geometric derivatives, and REDUNDANT_GEO for redundant total geometric derivatives;
  ! for instance, (xx,xy,yy,xz,yz,zz) are unique while (xx,yx,zx,xy,yy,zy,xz,yz,zz) are
  ! redundant, note that the "triangular" total geometric derivatives could be obtained
  ! from the unique total geometric derivatives by giving \var(max_num_cent)=\var(order_geo_total)
  integer, public, parameter :: UNIQUE_GEO = 1
  integer, public, parameter :: REDUNDANT_GEO = 3

  public :: Gen1IntShellGetIntExpt

  private :: gen1int_reorder_p_sgto

  contains

  !> \brief calculates property integrals and/or expectation values
  !> \author Bin Gao and T. Kjaergaard
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
  !> \return val_ints contains the integral matrices
  !> \return val_expt contains the expectation values
  !> \note the arrangement of var(val_ints) and \var(val_expt) will be in the order of
  !>       \var(order_mom), \var(order_mag_bra), ..., \var(order_geo_total), and each of
  !>       them is arranged in the order of (xx,xy,yy,xz,yz,zz) or (xx,yx,zx,xy,yy,zy,xz,yz,zz),
  !>       see Gen1Int library manual, for instance Section 2.2;
  !>       \var(val_expt) should be zero by users before calculations
  subroutine Gen1IntShellGetIntExpt(ao_item_bra, ao_item_ket,        &
                                    same_braket, do_screen,          &
                                    one_prop, geom_tree, geom_type,  &
                                    val_ints, write_ints,            &
                                    ao_dens, val_expt, write_expt,   &
                                    io_viewer)
    use AO_Type, only: AOITEM
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
    real(realk), target, optional, intent(inout) :: val_expt(:,:)
    logical, optional, intent(in) :: write_expt
    integer, intent(in) :: io_viewer
    integer num_prop                                 !number of property integrals
    integer prop_sym                                 !symmetry of property integrals
    integer order_geo                                !order of total geometric derivatives
    logical do_redunt_geo                            !calculates redundant total geometric derivatives
    integer path_num_unique                          !number of unique derivatives of current path
    integer path_num_redunt                          !number of redundnat derivatives of current path
    integer, allocatable :: redunt_list(:,:)         !list addresses of redundant total geometric derivatives
    integer num_redunt_geo                           !number of all redundant total geometric derivatives
    integer path_offset                              !offset of unique derivatives of current path
    integer num_matrices                             !number of integral matrices
    logical p_write_ints                             !if writing integral matrices on file
    logical p_write_expt                             !if writing expectation values on file
    logical do_integral                              !if returning and/or writing integral matrices
    logical do_expectation                           !if calculating or writing expectaion values on file
    integer num_dens                                 !number of AO density matrices
    integer start_expt                               !start address of expectation values
    integer end_expt                                 !end address of expectation values
    real(realk), pointer :: unique_expt(:,:)         !expectation values with unique geometric derivatives
    integer offset_prop                              !offset of property integrals
    integer offset_ao_bra                            !offset of AOs on bra center
    integer offset_ao_ket                            !offset of AOs on ket center
    integer addr_prop                                !address of property integrals
    integer addr_ao_bra                              !address of AOs on bra center
    integer addr_ao_ket                              !address of AOs on ket center
    integer angular_bra                              !angular number on bra center
    integer angular_ket                              !angular number on ket center
    real(realk) coord_bra(3)                         !coordinates of bra center
    real(realk) coord_ket(3)                         !coordinates of ket center
    integer num_prim_bra                             !number of primitive GTOs on bra center
    integer num_prim_ket                             !number of primitive GTOs on ket center
    integer num_contr_bra                            !number of contractions on bra center
    integer num_contr_ket                            !number of contractions on ket center
    integer num_gto_bra                              !number of GTOs on bra center
    integer num_gto_ket                              !number of orbitals on ket center
    integer num_ao_bra                               !number of AOs on bra center
    integer num_ao_ket                               !number of AOs on ket center
    real(realk), allocatable :: contr_ints(:,:,:,:)  !contracted integrals
!FIXME: to remove the powers? they should be defined in AOITEM?
    integer, allocatable :: powers_bra(:,:)          !Cartesian powers on bra center
    integer, allocatable :: powers_ket(:,:)          !Cartesian powers on ket center
    logical screen                                   !if screening
    integer ibtch, jbtch                             !incremental recorders over AO batches
    integer iang, jang                               !incremental recorders over angular momentum?
    integer igto                                     !incremental recorder over GTOs
    integer idens, igeo, iprop, ibra, iket           !incremental recorders
    integer xpow, ypow                               !incremental recorders over xy powers
    integer ierr                                     !error information
!FIXME: assumes SGTOs used, to be removed, should get from AOITEM
    logical, parameter :: SPHER_GTO = .true.
    screen=.FALSE.
    ! gets the number of property integrals
    call OnePropGetNumProp(one_prop=one_prop, num_prop=num_prop)
    ! gets the symmetry of property integrals
    if (same_braket) then
      call OnePropGetSymmetry(one_prop=one_prop, prop_sym=prop_sym)
    else
      prop_sym = SQUARE_INT_MAT
    end if
    ! sets total geometric derivatives
    if (present(geom_tree)) then
      ! gets the order of total geometric derivatives
      call GeomTreeGetOrder(geom_tree=geom_tree, order_geo=order_geo)
      ! sets the type of total geometric derivatives
      if (order_geo>1) then
        if (present(geom_type)) then
          if (geom_type/=UNIQUE_GEO .and. geom_type/=REDUNDANT_GEO) then
            call lsquit("Gen1IntShellGetIntExpt>> invalid type of total "// &
                        "geometric derivatives!", io_viewer)
          else
            do_redunt_geo = geom_type==REDUNDANT_GEO
          end if
        ! default are unique total geometric derivatives
        else
          do_redunt_geo = .false.
        end if
      ! the first order total geometric derivatives are unique
      else
        do_redunt_geo = .false.
      end if
      ! gets the number of unique derivatives of current path
      call GeomPathGetNumUnique(geom_tree=geom_tree, path_num_unique=path_num_unique)
      ! redundant total geometric derivatives
      if (do_redunt_geo) then
        ! gets the number of redundant derivatives of current path
        call GeomPathGetNumRedunt(geom_tree=geom_tree, &
                                  path_num_redunt=path_num_redunt)
        ! sets the list addresses of redundant total geometric derivatives
        allocate(redunt_list(2,path_num_redunt), stat=ierr)
        if (ierr/=0) then
          call lsquit("Gen1IntShellGetIntAve>> failed to allocate redunt_list!", &
                      io_viewer)
        end if
        call GeomPathGetReduntList(geom_tree=geom_tree, redunt_list=redunt_list)
        ! gets the number of all redundant derivatives in the N-ary tree
        call GeomTreeGetNumAtoms(geom_tree=geom_tree, num_atoms=num_redunt_geo)
        num_redunt_geo = (3*num_redunt_geo)**order_geo
        ! sets the offset of unique derivatives of current path as 0, which will only
        ! be used for calculating expecatation values \var(unique_expt)
        path_offset = 0
      else
        path_num_redunt = 1  !not used
        num_redunt_geo = 1   !not used
        ! gets the offset of unique derivatives of current path
        call GeomPathGetOffset(geom_tree=geom_tree, path_offset=path_offset)
      end if
    ! no total geometric derivatives
    else
      do_redunt_geo = .false.
      path_num_unique = 1
      path_num_redunt = 1  !not used
      num_redunt_geo = 1   !not used
      path_offset = 0
    end if
    ! sets the number of integral matrices
    num_matrices = num_prop*path_num_unique
    ! sets the start and end addresses of expectation values
    start_expt = path_offset*num_prop+1
    end_expt = path_offset*num_prop+num_matrices
    ! if writing integral matrices on file
    if (present(write_ints)) then
      p_write_ints = write_ints
    else
      p_write_ints = .false.
    end if
    ! if calculating contracted integrals (only needed for parallel mode)
    do_integral = present(val_ints) .or. p_write_ints
    ! if writing expectation values on file
    if (present(write_expt)) then
      p_write_expt = write_expt
    else
      p_write_expt = .false.
    end if
    ! if calculating expectation values
    do_expectation = present(ao_dens) .and. (present(val_expt) .or. p_write_expt)
    ! sets the information related to expectation values
    if (do_expectation) then
      ! gets the number of AO density matrices
      num_dens = size(ao_dens)
      ! if redundant total geometric derivatives are required or no output argument is
      ! given, we need to allocate memory to save the expectation values with the unique
      ! total geometric derivatives
      if (do_redunt_geo .or. .not.present(val_expt)) then
        allocate(unique_expt(start_expt:end_expt,num_dens), stat=ierr)
        if (ierr/=0) then
          call lsquit("Gen1IntShellGetIntExpt>> failed to allocate unique_expt!", &
                      io_viewer)
        end if
        unique_expt = 0.0
      ! for required unique total geometric derivatives, we just point \var(unique_expt) to them
      else
        unique_expt => val_expt
      end if
    end if
    ! sets the offset of AOs on ket center
    offset_ao_ket = 0
    ! loops over AO sub-shells on ket center
    do jbtch = 1, ao_item_ket%nbatches
      ! gets the number of primitive GTOs
      num_prim_ket = ao_item_ket%BATCH(jbtch)%nPrimitives
      ! gets the coordinates of ket center
      coord_ket(1) = ao_item_ket%BATCH(jbtch)%CENTER(1)
      coord_ket(2) = ao_item_ket%BATCH(jbtch)%CENTER(2)
      coord_ket(3) = ao_item_ket%BATCH(jbtch)%CENTER(3)
      do jang = 1, ao_item_ket%BATCH(jbtch)%nAngmom
        ! gets the angular number
        angular_ket = ao_item_ket%BATCH(jbtch)%Angmom(jang)
        ! gets the number of contractions
        num_contr_ket = ao_item_ket%BATCH(jbtch)%nContracted(jang)
        ! sets the number of atomic orbitals on ket center
        if (SPHER_GTO) then
          num_gto_ket = 2*angular_ket+1
        else
          num_gto_ket = (angular_ket+1)*(angular_ket+2)/2
          ! sets the Cartesian powers on ket center
          allocate(powers_ket(3,num_gto_ket), stat=ierr)
          if (ierr/=0) then
            call lsquit("Gen1IntShellGetIntExpt>> failed to allocate powers_ket!", &
                        io_viewer)
          end if
          ! Dalton's order of CGTOs, for instance, dxx, dxy, dxz, dyy, dyz, dzz
          igto = 0
          do xpow = angular_ket, 0, -1
            do ypow = angular_ket-xpow, 0, -1
              igto = igto+1
              powers_ket(1,igto) = xpow
              powers_ket(2,igto) = ypow
              powers_ket(3,igto) = angular_ket-(xpow+ypow)
            end do
          end do
        end if
        ! sets the number of AOs on ket center
        num_ao_ket = num_gto_ket*num_contr_ket
        ! sets the offset of AOs on bra center
        offset_ao_bra = 0
        ! loops over AO sub-shells on bra center
        do ibtch = 1, ao_item_bra%nbatches
          if (do_screen) then
            call getODscreening(ao_item_bra%BATCH(ibtch), &
                                ao_item_ket%BATCH(jbtch), screen)
          end if
          ! updates the offset of AOs on bra center
          if (screen) then
            do iang = 1, ao_item_bra%BATCH(ibtch)%nAngmom
              ! gets the angular number
              angular_bra = ao_item_bra%BATCH(ibtch)%Angmom(iang)
              ! gets the number of contractions
              num_contr_bra = ao_item_bra%BATCH(ibtch)%nContracted(iang)
              ! sets the number of atomic orbitals on bra center
              if (SPHER_GTO) then
                num_gto_bra = 2*angular_bra+1
              else
                num_gto_bra = (angular_bra+1)*(angular_bra+2)/2
              end if
              offset_ao_bra = offset_ao_bra+num_gto_bra*num_contr_bra
            end do
          else
            ! gets the number of primitive GTOs
            num_prim_bra = ao_item_bra%BATCH(ibtch)%nPrimitives
            ! gets the coordinates of bra center
            coord_bra(1) = ao_item_bra%BATCH(ibtch)%CENTER(1)
            coord_bra(2) = ao_item_bra%BATCH(ibtch)%CENTER(2)
            coord_bra(3) = ao_item_bra%BATCH(ibtch)%CENTER(3)
            do iang = 1, ao_item_bra%BATCH(ibtch)%nAngmom
              ! gets the angular number
              angular_bra = ao_item_bra%BATCH(ibtch)%Angmom(iang)
              ! gets the number of contractions
              num_contr_bra = ao_item_bra%BATCH(ibtch)%nContracted(iang)
              ! sets the number of atomic orbitals on bra center
              if (SPHER_GTO) then
                num_gto_bra = 2*angular_bra+1
              else
                num_gto_bra = (angular_bra+1)*(angular_bra+2)/2
                ! sets the Cartesian powers on bra center
                allocate(powers_bra(3,num_gto_bra), stat=ierr)
                if (ierr/=0) then
                  call lsquit("Gen1IntShellGetIntExpt>> failed to allocate powers_bra!", &
                              io_viewer)
                end if
                ! Dalton's order of CGTOs, for instance, dxx, dxy, dxz, dyy, dyz, dzz
                igto = 0
                do xpow = angular_bra, 0, -1
                  do ypow = angular_bra-xpow, 0, -1
                    igto = igto+1
                    powers_bra(1,igto) = xpow
                    powers_bra(2,igto) = ypow
                    powers_bra(3,igto) = angular_bra-(xpow+ypow)
                  end do
                end do
              end if
              ! sets the number of orbitals on bra center
              num_ao_bra = num_gto_bra*num_contr_bra
!FIXME: repeatedly allocate and deallocate may not be efficient?
              ! allocates memory for contracted integrals between two AO sub-shells
              allocate(contr_ints(num_ao_bra,num_ao_ket,num_prop,path_num_unique), stat=ierr)
              if (ierr/=0) then
                call lsquit("Gen1IntShellGetIntExpt>> failed to allocate contr_ints!", &
                            io_viewer)
              end if
              ! spherical GTOs
              if (SPHER_GTO) then
                ! calls Gen1Int subroutines to evaluate property integrals
                call OnePropGetIntegral(idx_bra=ao_item_bra%BATCH(ibtch)%atom,               &
                                        coord_bra=coord_bra,                                 &
                                        angular_bra=angular_bra,                             &
                                        num_prim_bra=num_prim_bra,                           &
                                        exponent_bra=ao_item_bra%BATCH(ibtch)%               &
                                                     pExponents%elms(1:num_prim_bra),        &
                                        num_contr_bra=num_contr_bra,                         &
                                        contr_coef_bra=ao_item_bra%BATCH(ibtch)%pCC(iang)%   &
                                                       p%elms(1:num_contr_bra*num_prim_bra), &
                                        idx_ket=ao_item_ket%BATCH(jbtch)%atom,               &
                                        coord_ket=coord_ket,                                 &
                                        angular_ket=angular_ket,                             &
                                        num_prim_ket=num_prim_ket,                           &
                                        exponent_ket=ao_item_ket%BATCH(jbtch)%               &
                                                     pExponents%elms(1:num_prim_ket),        &
                                        num_contr_ket=num_contr_ket,                         &
                                        contr_coef_ket=ao_item_ket%BATCH(jbtch)%pCC(jang)%   &
                                                       p%elms(1:num_contr_ket*num_prim_ket), &
                                        spher_gto=SPHER_GTO,                                 &
                                        one_prop=one_prop,                                   &
                                        geom_tree=geom_tree,                                 &
                                        num_gto_bra=num_gto_bra,                             &
                                        num_gto_ket=num_gto_ket,                             &
                                        num_opt=num_matrices,                                &
                                        contr_ints=contr_ints)
                ! reorders p-shell spherical GTOs, since Dalton uses x(+1), y(-1) and z(0),
                ! while Gen1Int uses y(-1), z(0) and x(+1)
                if (angular_bra==1)                                                       &
                  call gen1int_reorder_p_sgto(1, num_contr_bra*num_gto_ket*num_contr_ket, &
                                              num_matrices, contr_ints, io_viewer)
                if (angular_ket==1)                                                     &
                  call gen1int_reorder_p_sgto(num_gto_bra*num_contr_bra, num_contr_ket, &
                                              num_matrices, contr_ints, io_viewer)
              ! Cartesian GTOs
              else
                ! calls Gen1Int subroutines to evaluate property integrals, and reorders
                ! the integrals according to Cartesian powers
                call OnePropGetIntegral(idx_bra=ao_item_bra%BATCH(ibtch)%atom,               &
                                        coord_bra=coord_bra,                                 &
                                        angular_bra=angular_bra,                             &
                                        num_prim_bra=num_prim_bra,                           &
                                        exponent_bra=ao_item_bra%BATCH(ibtch)%               &
                                                     pExponents%elms(1:num_prim_bra),        &
                                        num_contr_bra=num_contr_bra,                         &
                                        contr_coef_bra=ao_item_bra%BATCH(ibtch)%pCC(iang)%   &
                                                       p%elms(1:num_contr_bra*num_prim_bra), &
                                        idx_ket=ao_item_ket%BATCH(jbtch)%atom,               &
                                        coord_ket=coord_ket,                                 &
                                        angular_ket=angular_ket,                             &
                                        num_prim_ket=num_prim_ket,                           &
                                        exponent_ket=ao_item_ket%BATCH(jbtch)%               &
                                                     pExponents%elms(1:num_prim_ket),        &
                                        num_contr_ket=num_contr_ket,                         &
                                        contr_coef_ket=ao_item_ket%BATCH(jbtch)%pCC(jang)%   &
                                                       p%elms(1:num_contr_ket*num_prim_ket), &
                                        spher_gto=SPHER_GTO,                                 &
                                        one_prop=one_prop,                                   &
                                        geom_tree=geom_tree,                                 &
                                        num_gto_bra=num_gto_bra,                             &
                                        num_gto_ket=num_gto_ket,                             &
                                        num_opt=num_matrices,                                &
                                        contr_ints=contr_ints,                               &
                                        powers_bra=powers_bra,                               &
                                        powers_ket=powers_ket)
                deallocate(powers_bra)
              end if
              ! assigns the returned integrals
              if (present(val_ints)) then
                ! returns integral matrices with redundant total geometric derivatives
                if (do_redunt_geo) then
                  do igeo = 1, path_num_redunt
                    offset_prop = num_prop*(redunt_list(2,igeo)-1)
                    do iprop = 1, num_prop
                      addr_prop = offset_prop+iprop
                      do iket = 1, num_ao_ket
!FIXME: it may not be safe to touch %nrow and %elms here, using subroutines from matrix module
                        addr_ao_ket = (iket+offset_ao_ket-1)*val_ints(addr_prop)%nrow
                        do ibra = 1, num_ao_bra
                          val_ints(addr_prop)%elms(ibra+offset_ao_bra+addr_ao_ket) &
                            = contr_ints(ibra,iket,iprop,redunt_list(1,igeo))
                        end do
                      end do
!FIXME: to implement?
                      !if (prop_sym==SYMM_INT_MAT .and. jbtch/=ibtch .and. jang/=iang) then
                      !else if (prop_sym==ANTI_INT_MAT .and. jbtch/=ibtch .and. jang/=iang) then
                      !end if
                    end do
                  end do
                ! returns integrals matrices with unique total geometric derivatives
                else
                  do igeo = 1, path_num_unique
                    offset_prop = num_prop*(path_offset+igeo-1)
                    do iprop = 1, num_prop
                      addr_prop = offset_prop+iprop
                      do iket = 1, num_ao_ket
!FIXME: it may not be safe to touch %nrow and %elms here, using subroutines from matrix module
                        addr_ao_ket = (iket+offset_ao_ket-1)*val_ints(addr_prop)%nrow
                        do ibra = 1, num_ao_bra
                          val_ints(addr_prop)%elms(ibra+offset_ao_bra+addr_ao_ket) &
                            = contr_ints(ibra,iket,iprop,igeo)
                        end do
                      end do
!FIXME: to implement? 
                      !if (prop_sym==SYMM_INT_MAT .and. jbtch/=ibtch .and. jang/=iang) then
                      !else if (prop_sym==ANTI_INT_MAT .and. jbtch/=ibtch .and. jang/=iang) then
                      !end if
                    end do
                  end do
                end if
              end if
!FIXME: to implement
              ! writes the integrals on file
              if (p_write_ints) then
              end if
              ! calculates the expectation values with unique total geometric derivatives
              if (do_expectation) then
                do idens = 1, num_dens
                  do igeo = 1, path_num_unique
                    offset_prop = (path_offset+igeo-1)*num_prop
                    do iprop = 1, num_prop
                      addr_prop = offset_prop+iprop
                      addr_ao_bra = (offset_ao_bra-1)*ao_dens(idens)%nrow
                      do ibra = 1, num_ao_bra
                        addr_ao_bra = addr_ao_bra+ao_dens(idens)%nrow
                        do iket = 1, num_ao_ket
                          unique_expt(addr_prop,idens)         &
                            = unique_expt(addr_prop,idens)     &
                            + contr_ints(ibra,iket,iprop,igeo) &
                            * ao_dens(idens)%elms(iket+offset_ao_ket+addr_ao_bra)
                        end do
                      end do
!FIXME: to implement?
                      !if (prop_sym==SYMM_INT_MAT .and. jbtch/=ibtch .and. jang/=iang) then
                      !else if (prop_sym==ANTI_INT_MAT .and. jbtch/=ibtch .and. jang/=iang) then
                      !end if
                    end do
                  end do
                end do
              end if
!FIXME: repeatedly allocate and deallocate may not be efficient?
              deallocate(contr_ints)
              ! updates the offset of AOs on bra center
              offset_ao_bra = offset_ao_bra+num_ao_bra
            end do !iang
          end if
        end do !ibtch
        if (.not.SPHER_GTO) deallocate(powers_ket)
        ! updates the offset of AOs on ket center
        offset_ao_ket = offset_ao_ket+num_ao_ket
      end do !jang
    end do !jbtch
    ! frees spaces
    if (allocated(redunt_list)) deallocate(redunt_list)
    if (do_expectation) then
!FIXME: to implement
      ! writes expectation values \var(unique_expt) on file
      if (p_write_expt) then
      end if
      ! returns expectation values with redundant total geometric derivatives
      if (present(val_expt)) then
        if (do_redunt_geo) then
          call GeomPathSetReduntExpt(geom_tree=geom_tree,                            &
                                     num_opt=num_prop,                               &
                                     path_num_unique=path_num_unique,                &
                                     num_dens=num_dens,                              &
                                     unique_expt=unique_expt(start_expt:end_expt,:), &
                                     num_redunt_geo=num_redunt_geo,                  &
                                     redunt_expt=val_expt)
          deallocate(unique_expt)
        end if
      else
        deallocate(unique_expt)
      end if
      nullify(unique_expt)
    end if
100 format("Gen1IntShellGetIntExpt>> ",A,I12)
  end subroutine Gen1IntShellGetIntExpt

  !> \brief reorders the p-shell contracted real solid-harmonic Gaussians in Dalton's order
  !> \author Bin Gao
  !> \date 2011-08-03
  !> \param dim_sgto_bra is the dimension of SGTOs on bra center
  !> \param num_contr_ket is the number of contractions of ket center
  !> \param num_opt is the number of operators
  !> \param gen_ints contains the contracted integrals from Gen1Int library
  !> \param io_viewer is the logical unit number of the viewer
  subroutine gen1int_reorder_p_sgto(dim_sgto_bra, num_contr_ket, num_opt, &
                                    gen_ints, io_viewer)
    integer, intent(in) :: dim_sgto_bra
    integer, intent(in) :: num_contr_ket
    integer, intent(in) :: num_opt
    real(realk), intent(inout) :: gen_ints(dim_sgto_bra,3,num_contr_ket,num_opt)
    integer, intent(in) :: io_viewer
    real(realk), allocatable :: tmp_ints(:)  !temporary integrals
    integer icontr, iopt                     !incremental recorders
    integer ierr                             !error information
    ! Dalton's order of SGTOs: px(1), py(-1), pz(0),
    ! while those in Gen1Int is: py(-1), pz(0), px(1)
    allocate(tmp_ints(dim_sgto_bra), stat=ierr)
    if (ierr/=0) then
      call lsquit("gen1int_reorder_p_sgto>> failed to allocate tmp_ints!", &
                  io_viewer)
    end if
    do iopt = 1, num_opt
      do icontr = 1, num_contr_ket
        tmp_ints = gen_ints(:,3,icontr,iopt)
        gen_ints(:,3,icontr,iopt) = gen_ints(:,2,icontr,iopt)
        gen_ints(:,2,icontr,iopt) = gen_ints(:,1,icontr,iopt)
        gen_ints(:,1,icontr,iopt) = tmp_ints
      end do
    end do
    deallocate(tmp_ints)
  end subroutine gen1int_reorder_p_sgto

end module gen1int_shell
