!
!...   Copyright (c) 2011 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2011 (2011), see http://daltonprogram.org"
!...
!...   This source code is provided under a written licence and may be
!...   used, copied, transmitted, or stored only in accord with that
!...   written licence.
!...
!...   In particular, no part of the source code or compiled modules may
!...   be distributed outside the research group of the licence holder.
!...   This means also that persons (e.g. post-docs) leaving the research
!...   group of the licence holder may not take any part of Dalton,
!...   including modified files, with him/her, unless that person has
!...   obtained his/her own licence.
!...
!...   For further information, including how to get a licence, see:
!...      http://daltonprogram.org
!
!
!...  This file takes care the atomic orbital (AO) sub-shell used in Gen1Int interface.
!
!...  2012-03-11, Bin Gao
!...  * adds subroutine Gen1IntShellMO to calculate molecular orbitals at grid points
!
!...  2011-10-02, Bin Gao
!...  * first version

#include "stdout.h"
#include "xkind.h"

!> \brief defines the AO sub-shell used in Gen1Int interface and corresponding subroutines
!> \author Bin Gao
!> \date 2011-10-02
module gen1int_shell

  ! Fortran 90 module of Gen1Int library
  use gen1int

  implicit none

  ! contracted Cartesian or spherical GTOs, or an atomic orbital (AO) sub-shell
  ! here, an AO sub-shell contains the orbitals with the same exponents, and the
  ! same/different contraction coefficients; for instance 2p and 3p sub-shells
  type, public :: sub_shell_t
    private
    ! if spherical GTOs
    logical spher_gto
    ! index of atomic center where this AO sub-shell locates
    integer idx_cent
    ! coordinates of the atomic center
    real(REALK) coord_cent(3)
    ! angular number
    integer ang_num
    ! number of primitive GTOs
    integer num_prim
    ! exponents of primitive GTOs
    real(REALK), allocatable :: exponents(:)
    ! number of contractions
    integer num_contr
    ! contraction coefficients
    real(REALK), allocatable :: contr_coef(:,:)
    ! number of atomic orbitals
    integer num_ao
    !-Dalton uses the same sequence of SGTOs as those in Gen1Int except for p-shell
    !-which could be done by hand coding
    !-! magnetic numbers if spherical GTOs
    !-integer, allocatable :: mag_num(:)
    ! Cartesian powers if Cartesian GTOs
    integer, allocatable :: powers(:,:)
    ! base index of the atomic orbitals in this sub-shell, i.e., the index of the
    ! first atomic orbital in this shell minus 1
    integer base_idx
  end type sub_shell_t

  public :: Gen1IntShellCreate
  public :: Gen1IntShellAttr
  public :: Gen1IntShellGetRangeAO
  public :: Gen1IntShellView
  public :: Gen1IntShellIntegral
  public :: Gen1IntShellMO
  public :: Gen1IntShellDestroy

  contains

  !> \brief creates an AO sub-shell
  !> \author Bin Gao
  !> \date 2011-10-03
  !> \param spher_gto indicates if spherical GTOs
  !> \param idx_cent is the index of atomic center where this AO sub-shell locates
  !> \param coord_cent contains the coordinates of the atomic center
  !> \param ang_num is the angular number
  !> \param num_prim is the number of primitive GTOs
  !> \param exponents contains the exponents of primitive GTOs
  !> \param num_contr is the number of contractions
  !> \param contr_coef contains the contraction coefficients
  !> \param powers contains the Cartesian powers if Cartesian GTOs
  !> \param last_shell is the last sub-shell before this AO sub-shell
  !> \return sub_shell is the initialized AO sub-shell
  subroutine Gen1IntShellCreate(spher_gto, idx_cent, coord_cent, ang_num,   &
                                num_prim, exponents, num_contr, contr_coef, &
                                powers, last_shell, sub_shell)
    logical, intent(in) :: spher_gto
    integer, intent(in) :: idx_cent
    real(REALK), intent(in) :: coord_cent(3)
    integer, intent(in) :: ang_num
    integer, intent(in) :: num_prim
    real(REALK), intent(in) :: exponents(num_prim)
    integer, intent(in) :: num_contr
    real(REALK), intent(in) :: contr_coef(num_contr,num_prim)
    integer, optional, intent(in) :: powers(3,(ang_num+1)*(ang_num+2)/2)
    type(sub_shell_t), optional, intent(in) :: last_shell
    type(sub_shell_t), intent(inout) :: sub_shell
    integer iao         !incremental recorder
    integer xpow, ypow  !incremental recorders over xy powers
    integer ierr        !error information
    sub_shell%spher_gto = spher_gto
    sub_shell%idx_cent = idx_cent
    sub_shell%coord_cent = coord_cent
    sub_shell%ang_num = ang_num
    sub_shell%num_prim = num_prim
    allocate(sub_shell%exponents(num_prim), stat=ierr)
    if (ierr/=0) stop "Gen1IntShellCreate>> failed to allocate exponents!"
    sub_shell%exponents = exponents
    sub_shell%num_contr = num_contr
    allocate(sub_shell%contr_coef(num_contr,num_prim), stat=ierr)
    if (ierr/=0) stop "Gen1IntShellCreate>> failed to allocate contr_coef!"
    sub_shell%contr_coef = contr_coef
    ! spherical GTOs
    if (spher_gto) then
      sub_shell%num_ao = 2*ang_num+1
      !-allocate(sub_shell%mag_num(sub_shell%num_ao), stat=ierr)
      !-if (ierr/=0) stop "Gen1IntShellCreate>> failed to allocate mag_num!"
      !-if (present(mag_num)) then
      !-  sub_shell%mag_num = mag_num
      !-else
      !-  ! Dalton's order of SGTOs
      !-  if (ang_num==1) then
      !-    sub_shell%mag_num(1) = 1   !px
      !-    sub_shell%mag_num(2) = -1  !py
      !-    sub_shell%mag_num(3) = 0   !pz
      !-  else
      !-    do iao = -ang_num, ang_num
      !-      sub_shell%mag_num(ang_num+iao+1) = iao
      !-    end do
      !-  end if
      !-end if
    ! Cartesian GTOs
    else
      sub_shell%num_ao = (ang_num+1)*(ang_num+2)/2
      allocate(sub_shell%powers(3,sub_shell%num_ao), stat=ierr)
      if (ierr/=0) stop "Gen1IntShellCreate>> failed to allocate powers!"
      if (present(powers)) then
        sub_shell%powers = powers
      ! Dalton's order of CGTOs, for instance, dxx, dxy, dxz, dyy, dyz, dzz
      else
        iao = 0
        do xpow = ang_num, 0, -1
          do ypow = ang_num-xpow, 0, -1
            iao = iao+1
            sub_shell%powers(1,iao) = xpow
            sub_shell%powers(2,iao) = ypow
            sub_shell%powers(3,iao) = ang_num-(xpow+ypow)
          end do
        end do
      end if  
    end if
    if (present(last_shell)) then
      sub_shell%base_idx = last_shell%base_idx &
                         + last_shell%num_ao*last_shell%num_contr
    ! this is the first sub-shell
    else
      sub_shell%base_idx = 0
    end if
  end subroutine Gen1IntShellCreate

  !> \brief gets the numbers of atomic orbitals and contractions of an AO sub-shell
  !> \author Bin Gao
  !> \date 2011-10-04
  !> \param sub_shell is the AO sub-shell
  !> \return num_ao is the number of atomic orbitals
  !> \return num_contr is the number of contractions
  subroutine Gen1IntShellAttr(sub_shell, num_ao, num_contr)
    type(sub_shell_t), intent(in) :: sub_shell
    integer, intent(out) :: num_ao
    integer, intent(out) :: num_contr
    num_ao = sub_shell%num_ao
    num_contr = sub_shell%num_contr
  end subroutine Gen1IntShellAttr

  !> \brief gets the indices of the first and last orbtials for the given AO sub-shell
  !> \author Bin Gao
  !> \date 2011-10-06
  !> \param sub_shell is the AO sub-shell
  !> \return idx_first is the index of the first orbital
  !> \return idx_last is the index of the last orbital
  subroutine Gen1IntShellGetRangeAO(sub_shell, idx_first, idx_last)
    type(sub_shell_t), intent(in) :: sub_shell
    integer, intent(out) :: idx_first
    integer, intent(out) :: idx_last
    idx_first = sub_shell%base_idx+1
    idx_last = sub_shell%base_idx+sub_shell%num_ao*sub_shell%num_contr
  end subroutine Gen1IntShellGetRangeAO

  !> \brief visualizes the information of several AO sub-shells
  !> \author Bin Gao
  !> \date 2011-10-03
  !> \param num_shells is the number of AO sub-shells
  !> \param sub_shells are the AO sub-shells
  !> \param io_viewer is the logical unit number of the viewer
  subroutine Gen1IntShellView(num_shells, sub_shells, io_viewer)
    integer, intent(in) :: num_shells
    type(sub_shell_t), intent(in) :: sub_shells(*)
    integer, intent(in) :: io_viewer
    integer ishell, iprim  !incremental recorder
    write(io_viewer,100) "number of AO sub-shells", num_shells
    do ishell = 1, num_shells
      write(io_viewer,100) "sub-shell", ishell
      write(io_viewer,100) "index of atomic center of the sub-shell", &
                         sub_shells(ishell)%idx_cent
      write(io_viewer,110) "coordinates of the atomic center", &
                         sub_shells(ishell)%coord_cent
      write(io_viewer,100) "angular number", sub_shells(ishell)%ang_num
      write(io_viewer,100) "number of primitive GTOs", sub_shells(ishell)%num_prim
      write(io_viewer,100) "number of contractions", sub_shells(ishell)%num_contr
      write(io_viewer,100) "    exponents      coefficients"
      do iprim = 1, sub_shells(ishell)%num_prim
        write(io_viewer,120)                   &
          sub_shells(ishell)%exponents(iprim), &
          sub_shells(ishell)%contr_coef(1:min(4,sub_shells(ishell)%num_contr),iprim)
        write(io_viewer,130) &
          sub_shells(ishell)%contr_coef(5:sub_shells(ishell)%num_contr,iprim)
      end do
      if (sub_shells(ishell)%spher_gto) then
        write(io_viewer,100) "SGTOs used"
        !-write(io_viewer,140) sub_shells(ishell)%mag_num
      else
        write(io_viewer,100) "CGTOs used"
        write(io_viewer,150) sub_shells(ishell)%powers
      end if
      write(io_viewer,100) "base index of the orbitals", sub_shells(ishell)%base_idx
    end do
100 format("Gen1IntShellView>> "A,2I8)
110 format("Gen1IntShellView>> "A,3F16.8)
120 format("Gen1IntShellView>> ",5Es16.8)
130 format("Gen1IntShellView>> ",16X,4Es16.8)
140 format("Gen1IntShellView>> magnetic numbers>> ",8I4)
150 format("Gen1IntShellView>> Cartesian powers>> ",4(3I4,2X))
  end subroutine Gen1IntShellView

  !> \brief calculates property integrals for given AO sub-shells
  !> \author Bin Gao
  !> \date 2011-10-07
  !> \param num_shells is the number of AO sub-shells
  !> \param sub_shells are the AO sub-shells
  !> \param one_prop contains the information of one-electron property integrals
  !> \param london_ao indicates if using London atomic orbitals
  !> \param order_mag_bra is the order of magnetic derivatives on bra center
  !> \param order_mag_ket is the order of magnetic derivatives on ket center
  !> \param order_mag_total is the order of total magnetic derivatives
  !> \param order_ram_bra is the order of derivatives w.r.t. total rotational angular momentum on bra center
  !> \param order_ram_ket is the order of derivatives w.r.t. total rotational angular momentum on ket center
  !> \param order_ram_total is the order of total derivatives w.r.t. total rotational angular momentum
  !> \param order_geo_bra is the order of geometric derivatives with respect to bra center
  !> \param order_geo_ket is the order of geometric derivatives with respect to ket center
  !> \param geom_tree contains the information of n-ary tree for total geometric derivatives
  !> \param num_ints is the number of integral matrices including kinds of derivatives
  !> \param redunt_ints indicates if integral matrices with redundant total geometric derivatives returned
  !> \param wrt_ints indicates if writing integral matrices on file
  !> \param num_dens is the number of AO density matrices
  !> \param ao_dens contains the AO density matrices
  !> \param redunt_expt indicates if expectation values with redundant total geometric derivatives returned
  !> \param wrt_expt indicates if writing expectation values on file
  !> \return val_ints contains the integral matrices
  !> \return val_expt contains the expectation values
  subroutine Gen1IntShellIntegral(num_shells, sub_shells, one_prop, london_ao,   &
                                  order_mag_bra, order_mag_ket, order_mag_total, &
                                  order_ram_bra, order_ram_ket, order_ram_total, &
                                  order_geo_bra, order_geo_ket, geom_tree,       &
                                  num_ints, val_ints, redunt_ints, wrt_ints,     &
                                  num_dens, ao_dens, val_expt, redunt_expt, wrt_expt)
    ! matrix module
    use gen1int_matrix
    integer, intent(in) :: num_shells
    type(sub_shell_t), intent(in) :: sub_shells(*)
    type(one_prop_t), intent(in) :: one_prop
    logical, optional, intent(in) :: london_ao
    integer, optional, intent(in) :: order_mag_bra
    integer, optional, intent(in) :: order_mag_ket
    integer, optional, intent(in) :: order_mag_total
    integer, optional, intent(in) :: order_ram_bra
    integer, optional, intent(in) :: order_ram_ket
    integer, optional, intent(in) :: order_ram_total
    integer, optional, intent(in) :: order_geo_bra
    integer, optional, intent(in) :: order_geo_ket
    type(geom_tree_t), optional, intent(in) :: geom_tree
    integer, intent(in) :: num_ints
    type(matrix), optional, intent(inout) :: val_ints(num_ints)
    logical, optional, intent(in) :: redunt_ints
    logical, optional, intent(in) :: wrt_ints
    integer, intent(in) :: num_dens
    type(matrix), optional, intent(in) :: ao_dens(num_dens)
    real(REALK), optional, intent(inout) :: val_expt(num_ints*num_dens)
    logical, optional, intent(in) :: redunt_expt
    logical, optional, intent(in) :: wrt_expt
    logical p_wrt_ints       !if writing integral matrices on file
    logical p_redunt_ints    !if returning matrices with redundant total geometric derivatives
    integer dim_unique_geo   !dimension of all unique total geometric derivatives
    integer base_unique_geo  !base address of unique geometric derivatives for current path
    integer dim_redunt_geo   !dimension of all redundant total geometric derivatives
    integer order_geo_total  !order of total geometric derivatives
    logical p_wrt_expt       !if writing expectation values on file
    logical evaluate_expt    !if evaluating expectation values
    logical p_redunt_expt    !if returning expectation values with redundant total geometric derivatives
    integer num_prop         !number of integral matrices
    integer kind_prop        !kind of integral matrices
    integer num_unique_geo   !number of unique total geometric derivatives
    integer num_opt          !number of operators including derivatives
    integer ishell, jshell   !incremental recorders over AO sub-shells
    logical spher_gto        !if spherical GTOs
    integer max_num_ao       !maximum number of AOs in sub-shells
    integer max_num_contr    !maximum number of contractions of sub-shells
    real(REALK), allocatable :: contr_ints(:,:,:,:,:,:)  !contracted integrals between two AO sub-shells
    real(REALK), allocatable :: unique_expt(:,:,:)       !expectation values with unique geometric derivatives
    integer, allocatable :: redunt_list(:,:)             !list addresses of redundant total geometric derivatives
    integer num_redunt_geo                               !number of redundant total geometric derivatives
    integer iopt, iprop, igeo, idens                     !incremental recorders
    integer min_row_idx, max_row_idx                     !minimum and maximum indices of rows
    integer min_col_idx, max_col_idx                     !minimum and maximum indices of columns
    integer ierr                                         !error information
    if (present(wrt_ints)) then
      p_wrt_ints = wrt_ints
    else
      p_wrt_ints = .false.
    end if
    ! sets the information of total geometric derivatives
    if (present(geom_tree)) then
      call GeomTreeAttr(geom_tree=geom_tree, num_atoms=dim_redunt_geo, &
                        order_geo=order_geo_total)
      call GeomTreeDimUnique(geom_tree, dim_unique_geo)
      call GeomTreeBaseUnique(geom_tree, base_unique_geo)
      dim_redunt_geo = (3*dim_redunt_geo)**order_geo_total
    else
      order_geo_total = 0
      dim_unique_geo = 1
      base_unique_geo = 0
      dim_redunt_geo = 1
    end if
    ! returns integral matrices with redundant high order (>1) total geometric derivatives
    if ((present(val_ints) .or. p_wrt_ints) .and. &
        present(redunt_ints) .and. present(geom_tree)) then
      p_redunt_ints = redunt_ints .and. order_geo_total>1
    else
      p_redunt_ints = .false.
    end if
    if (present(wrt_expt)) then
      p_wrt_expt = wrt_expt
    else
      p_wrt_expt = .false.
    end if
    evaluate_expt = present(ao_dens) .and. (present(val_expt) .or. p_wrt_expt)
    ! returns expectation values with redundant high order (>1) total geometric derivatives
    if (evaluate_expt .and. present(redunt_expt) .and. present(geom_tree)) then
      p_redunt_expt = redunt_expt .and. order_geo_total>1
    else
      p_redunt_expt = .false.
    end if
    ! gets the number and kind of integral matrices
    call OnePropAttr(one_prop, num_prop, kind_prop)
    ! gets the number of derivatives
    if (present(order_mag_bra))   &
      num_prop = num_prop*(order_mag_bra+1)*(order_mag_bra+2)/2
    if (present(order_mag_ket))   &
       num_prop = num_prop*(order_mag_ket+1)*(order_mag_ket+2)/2
    if (present(order_mag_total)) &
       num_prop = num_prop*(order_mag_total+1)*(order_mag_total+2)/2
    if (present(order_ram_bra))   &
      num_prop = num_prop*(order_ram_bra+1)*(order_ram_bra+2)/2
    if (present(order_ram_ket))   &
      num_prop = num_prop*(order_ram_ket+1)*(order_ram_ket+2)/2
    if (present(order_ram_total)) &
      num_prop = num_prop*(order_ram_total+1)*(order_ram_total+2)/2
    if (present(order_geo_bra))   &
      num_prop = num_prop*(order_geo_bra+1)*(order_geo_bra+2)/2
    if (present(order_geo_ket))   &
      num_prop = num_prop*(order_geo_ket+1)*(order_geo_ket+2)/2
    ! gets the number of unique geometric derivatives
    if (present(geom_tree)) then
      call GeomTreeNumUnique(geom_tree, num_unique_geo)
    else
      num_unique_geo = 1
    end if
    ! number of operators
    num_opt = num_prop*num_unique_geo
    ! allocates the expectation values with unique geometric derivatives
    if (p_redunt_expt) then
      allocate(unique_expt(num_unique_geo,num_prop,num_dens), stat=ierr)
      if (ierr/=0) call QUIT("Gen1IntShellIntegral>> failed to allocate unique_expt!")
      unique_expt = 0.0_REALK
    end if
    select case(kind_prop)
    ! symmetric integral matrices (column- or ket-major, upper and diagonal parts)
    case(SYMM_INT_MAT)
      ! loops over AO sub-shells on ket center
      do jshell = 1, num_shells
        ! sets the minimum and maximum of indices of columns of the integral matrices
        min_col_idx = sub_shells(jshell)%base_idx+1
        max_col_idx = sub_shells(jshell)%base_idx &
                    + sub_shells(jshell)%num_ao*sub_shells(jshell)%num_contr
        ! loops over AO sub-shells on bra center (upper and diagonal parts)
        do ishell = 1, jshell
          spher_gto = sub_shells(ishell)%spher_gto .or. sub_shells(jshell)%spher_gto
          allocate(contr_ints(sub_shells(ishell)%num_ao,    &
                              sub_shells(ishell)%num_contr, &
                              sub_shells(jshell)%num_ao,    &
                              sub_shells(jshell)%num_contr, &
                              num_prop,num_unique_geo), stat=ierr)
          if (ierr/=0) call QUIT("Gen1IntShellIntegral>> failed to allocate contr_ints!")
          ! spherical GTOs
          if (spher_gto) then
            ! calls Gen1Int subroutines to evaluate property integrals
            call OnePropEvaluate(idx_bra=sub_shells(ishell)%idx_cent,          &
                                 coord_bra=sub_shells(ishell)%coord_cent,      &
                                 angular_bra=sub_shells(ishell)%ang_num,       &
                                 num_prim_bra=sub_shells(ishell)%num_prim,     &
                                 exponent_bra=sub_shells(ishell)%exponents,    &
                                 num_contr_bra=sub_shells(ishell)%num_contr,   &
                                 contr_coef_bra=sub_shells(ishell)%contr_coef, &
                                 idx_ket=sub_shells(jshell)%idx_cent,          &
                                 coord_ket=sub_shells(jshell)%coord_cent,      &
                                 angular_ket=sub_shells(jshell)%ang_num,       &
                                 num_prim_ket=sub_shells(jshell)%num_prim,     &
                                 exponent_ket=sub_shells(jshell)%exponents,    &
                                 num_contr_ket=sub_shells(jshell)%num_contr,   &
                                 contr_coef_ket=sub_shells(jshell)%contr_coef, &
                                 spher_gto=spher_gto,                          &
                                 london_ao=london_ao,                          &
                                 one_prop=one_prop,                            &
                                 order_mag_bra=order_mag_bra,                  &
                                 order_mag_ket=order_mag_ket,                  &
                                 order_mag_total=order_mag_total,              &
                                 order_ram_bra=order_ram_bra,                  &
                                 order_ram_ket=order_ram_ket,                  &
                                 order_ram_total=order_ram_total,              &
                                 order_geo_bra=order_geo_bra,                  &
                                 order_geo_ket=order_geo_ket,                  &
                                 geom_tree=geom_tree,                          &
                                 num_gto_bra=sub_shells(ishell)%num_ao,        &
                                 num_gto_ket=sub_shells(jshell)%num_ao,        &
                                 num_opt=num_opt, contr_ints=contr_ints)
            ! reorders p-shell spherical GTOs
            if (sub_shells(ishell)%ang_num==1)                      &
              call reorder_p_sgto(1, sub_shells(ishell)%num_contr   &
                                     *sub_shells(jshell)%num_ao     &
                                     *sub_shells(jshell)%num_contr, &
                                  num_opt, contr_ints)
            if (sub_shells(jshell)%ang_num==1)                   &
              call reorder_p_sgto(sub_shells(ishell)%num_ao      &
                                  *sub_shells(ishell)%num_contr, &
                                  sub_shells(jshell)%num_contr,  &
                                  num_opt, contr_ints)
          ! Cartesian GTOs
          else
            ! calls Gen1Int subroutines to evaluate property integrals, and reorders
            ! the integrals according to Cartesian powers
            call OnePropEvaluate(idx_bra=sub_shells(ishell)%idx_cent,          &
                                 coord_bra=sub_shells(ishell)%coord_cent,      &
                                 angular_bra=sub_shells(ishell)%ang_num,       &
                                 num_prim_bra=sub_shells(ishell)%num_prim,     &
                                 exponent_bra=sub_shells(ishell)%exponents,    &
                                 num_contr_bra=sub_shells(ishell)%num_contr,   &
                                 contr_coef_bra=sub_shells(ishell)%contr_coef, &
                                 idx_ket=sub_shells(jshell)%idx_cent,          &
                                 coord_ket=sub_shells(jshell)%coord_cent,      &
                                 angular_ket=sub_shells(jshell)%ang_num,       &
                                 num_prim_ket=sub_shells(jshell)%num_prim,     &
                                 exponent_ket=sub_shells(jshell)%exponents,    &
                                 num_contr_ket=sub_shells(jshell)%num_contr,   &
                                 contr_coef_ket=sub_shells(jshell)%contr_coef, &
                                 spher_gto=spher_gto,                          &
                                 london_ao=london_ao,                          &
                                 one_prop=one_prop,                            &
                                 order_mag_bra=order_mag_bra,                  &
                                 order_mag_ket=order_mag_ket,                  &
                                 order_mag_total=order_mag_total,              &
                                 order_ram_bra=order_ram_bra,                  &
                                 order_ram_ket=order_ram_ket,                  &
                                 order_ram_total=order_ram_total,              &
                                 order_geo_bra=order_geo_bra,                  &
                                 order_geo_ket=order_geo_ket,                  &
                                 geom_tree=geom_tree,                          &
                                 num_gto_bra=sub_shells(ishell)%num_ao,        &
                                 num_gto_ket=sub_shells(jshell)%num_ao,        &
                                 num_opt=num_opt, contr_ints=contr_ints,       &
                                 powers_bra=sub_shells(ishell)%powers,         &
                                 powers_ket=sub_shells(jshell)%powers)
          end if
          ! sets the minimum and maximum of indices of rows of the integral matrices
          min_row_idx = sub_shells(ishell)%base_idx+1
          max_row_idx = sub_shells(ishell)%base_idx &
                      + sub_shells(ishell)%num_ao*sub_shells(ishell)%num_contr
!FIXME: the returned index order of \var(num_prop)
          ! assigns the returned integrals, and write the integrals on file if required
          if (present(val_ints)) then
            ! returns integral matrices with redundant high order (>1) total geometric derivatives
            if (p_redunt_ints) then
              call GeomTreeNumRedunt(geom_tree, num_redunt_geo)
              allocate(redunt_list(2,num_redunt_geo), stat=ierr)
              if (ierr/=0) call QUIT("Gen1IntShellIntegral>> failed to allocate redunt_list!")
              call GeomTreeReduntList(geom_tree, redunt_list)
              iopt = 0
              do iprop = 1, num_prop
                do igeo = 1, num_redunt_geo
                  call MatSetBlockedValues(val_ints(iopt+redunt_list(2,igeo)), &
                                           min_row_idx, max_row_idx,           &
                                           min_col_idx, max_col_idx,           &
                                           contr_ints(:,:,:,:,iprop,           &
                                                      redunt_list(1,igeo)),    &
                                           .false.)
                  if (ishell/=jshell)                                            &
                    call MatSetBlockedValues(val_ints(iopt+redunt_list(2,igeo)), &
                                             min_row_idx, max_row_idx,           &
                                             min_col_idx, max_col_idx,           &
                                             contr_ints(:,:,:,:,iprop,           &
                                                        redunt_list(1,igeo)),    &
                                             .true.)
                end do
                iopt = iopt+dim_redunt_geo
              end do
              deallocate(redunt_list)
            ! returns integrals matrices with unique total geometric derivatives
            else
              iopt = 0
              do iprop = 1, num_prop
                do igeo = 1, num_unique_geo
                  call MatSetBlockedValues(val_ints(iopt+base_unique_geo+igeo), &
                                           min_row_idx, max_row_idx,            &
                                           min_col_idx, max_col_idx,            &
                                           contr_ints(:,:,:,:,iprop,igeo), .false.)
                  if (ishell/=jshell)                                             &
                    call MatSetBlockedValues(val_ints(iopt+base_unique_geo+igeo), &
                                             min_row_idx, max_row_idx,            &
                                             min_col_idx, max_col_idx,            &
                                             contr_ints(:,:,:,:,iprop,igeo), .true.)
                end do
                iopt = iopt+dim_unique_geo
              end do
            end if
          end if
          ! calculates the expectation values
          if (evaluate_expt) then
            ! returns expectation values with redundant high order (>1) total geometric derivatives
            if (p_redunt_expt) then
              do idens = 1, num_dens
                do igeo = 1, num_unique_geo
                  do iprop = 1, num_prop
                    call MatMultBlockedTrace(ao_dens(idens),                 &
                                             min_row_idx, max_row_idx,       &
                                             min_col_idx, max_col_idx,       &
                                             contr_ints(:,:,:,:,iprop,igeo), &
                                             unique_expt(igeo,iprop,idens), .false.)
                    if (ishell/=jshell)                                        &
                      call MatMultBlockedTrace(ao_dens(idens),                 &
                                               min_row_idx, max_row_idx,       &
                                               min_col_idx, max_col_idx,       &
                                               contr_ints(:,:,:,:,iprop,igeo), &
                                               unique_expt(igeo,iprop,idens), .true.)
                  end do
                end do
              end do
            ! returns expectation values with unique total geometric derivatives
            else
              do idens = 1, num_dens
                do igeo = 1, num_unique_geo
                  do iprop = 1, num_prop
                    iopt = base_unique_geo+igeo+(iprop-1)*dim_unique_geo &
                         + (idens-1)*dim_unique_geo*num_prop
                    call MatMultBlockedTrace(ao_dens(idens),                 &
                                             min_row_idx, max_row_idx,       &
                                             min_col_idx, max_col_idx,       &
                                             contr_ints(:,:,:,:,iprop,igeo), &
                                             val_expt(iopt), .false.)
                    if (ishell/=jshell)                                        &
                      call MatMultBlockedTrace(ao_dens(idens),                 &
                                               min_row_idx, max_row_idx,       &
                                               min_col_idx, max_col_idx,       &
                                               contr_ints(:,:,:,:,iprop,igeo), &
                                               val_expt(iopt), .true.)
                  end do
                end do
              end do
            end if
          end if
          deallocate(contr_ints)
        end do
      end do
    ! anti-symmetric integral matrices (column- or ket-major, upper and diagonal parts)
    case(ANTI_INT_MAT)
      ! loops over AO sub-shells on ket center
      do jshell = 1, num_shells
        ! sets the minimum and maximum of indices of columns of the integral matrices
        min_col_idx = sub_shells(jshell)%base_idx+1
        max_col_idx = sub_shells(jshell)%base_idx & 
                    + sub_shells(jshell)%num_ao*sub_shells(jshell)%num_contr
        ! loops over AO sub-shells on bra center (upper and diagonal parts)
        do ishell = 1, jshell
          spher_gto = sub_shells(ishell)%spher_gto .or. sub_shells(jshell)%spher_gto
          allocate(contr_ints(sub_shells(ishell)%num_ao,    &
                              sub_shells(ishell)%num_contr, &
                              sub_shells(jshell)%num_ao,    &
                              sub_shells(jshell)%num_contr, &
                              num_prop,num_unique_geo), stat=ierr)
          if (ierr/=0) call QUIT("Gen1IntShellIntegral>> failed to allocate contr_ints!")
          ! spherical GTOs
          if (spher_gto) then
            ! calls Gen1Int subroutines to evaluate property integrals
            call OnePropEvaluate(idx_bra=sub_shells(ishell)%idx_cent,          &
                                 coord_bra=sub_shells(ishell)%coord_cent,      &
                                 angular_bra=sub_shells(ishell)%ang_num,       &
                                 num_prim_bra=sub_shells(ishell)%num_prim,     &
                                 exponent_bra=sub_shells(ishell)%exponents,    &
                                 num_contr_bra=sub_shells(ishell)%num_contr,   &
                                 contr_coef_bra=sub_shells(ishell)%contr_coef, &
                                 idx_ket=sub_shells(jshell)%idx_cent,          &
                                 coord_ket=sub_shells(jshell)%coord_cent,      &
                                 angular_ket=sub_shells(jshell)%ang_num,       &
                                 num_prim_ket=sub_shells(jshell)%num_prim,     &
                                 exponent_ket=sub_shells(jshell)%exponents,    &
                                 num_contr_ket=sub_shells(jshell)%num_contr,   &
                                 contr_coef_ket=sub_shells(jshell)%contr_coef, &
                                 spher_gto=spher_gto,                          &
                                 london_ao=london_ao,                          &
                                 one_prop=one_prop,                            &
                                 order_mag_bra=order_mag_bra,                  &
                                 order_mag_ket=order_mag_ket,                  &
                                 order_mag_total=order_mag_total,              &
                                 order_ram_bra=order_ram_bra,                  &
                                 order_ram_ket=order_ram_ket,                  &
                                 order_ram_total=order_ram_total,              &
                                 order_geo_bra=order_geo_bra,                  &
                                 order_geo_ket=order_geo_ket,                  &
                                 geom_tree=geom_tree,                          &
                                 num_gto_bra=sub_shells(ishell)%num_ao,        &
                                 num_gto_ket=sub_shells(jshell)%num_ao,        &
                                 num_opt=num_opt, contr_ints=contr_ints)
            ! reorders p-shell spherical GTOs
            if (sub_shells(ishell)%ang_num==1)                      &
              call reorder_p_sgto(1, sub_shells(ishell)%num_contr   &
                                     *sub_shells(jshell)%num_ao     &
                                     *sub_shells(jshell)%num_contr, &
                                  num_opt, contr_ints)
            if (sub_shells(jshell)%ang_num==1)                   &
              call reorder_p_sgto(sub_shells(ishell)%num_ao      &
                                  *sub_shells(ishell)%num_contr, &
                                  sub_shells(jshell)%num_contr,  &
                                  num_opt, contr_ints)
          ! Cartesian GTOs
          else
            ! calls Gen1Int subroutines to evaluate property integrals, and reorders
            ! the integrals according to Cartesian powers
            call OnePropEvaluate(idx_bra=sub_shells(ishell)%idx_cent,          &
                                 coord_bra=sub_shells(ishell)%coord_cent,      &
                                 angular_bra=sub_shells(ishell)%ang_num,       &
                                 num_prim_bra=sub_shells(ishell)%num_prim,     &
                                 exponent_bra=sub_shells(ishell)%exponents,    &
                                 num_contr_bra=sub_shells(ishell)%num_contr,   &
                                 contr_coef_bra=sub_shells(ishell)%contr_coef, &
                                 idx_ket=sub_shells(jshell)%idx_cent,          &
                                 coord_ket=sub_shells(jshell)%coord_cent,      &
                                 angular_ket=sub_shells(jshell)%ang_num,       &
                                 num_prim_ket=sub_shells(jshell)%num_prim,     &
                                 exponent_ket=sub_shells(jshell)%exponents,    &
                                 num_contr_ket=sub_shells(jshell)%num_contr,   &
                                 contr_coef_ket=sub_shells(jshell)%contr_coef, &
                                 spher_gto=spher_gto,                          &
                                 london_ao=london_ao,                          &
                                 one_prop=one_prop,                            &
                                 order_mag_bra=order_mag_bra,                  &
                                 order_mag_ket=order_mag_ket,                  &
                                 order_mag_total=order_mag_total,              &
                                 order_ram_bra=order_ram_bra,                  &
                                 order_ram_ket=order_ram_ket,                  &
                                 order_ram_total=order_ram_total,              &
                                 order_geo_bra=order_geo_bra,                  &
                                 order_geo_ket=order_geo_ket,                  &
                                 geom_tree=geom_tree,                          &
                                 num_gto_bra=sub_shells(ishell)%num_ao,        &
                                 num_gto_ket=sub_shells(jshell)%num_ao,        &
                                 num_opt=num_opt, contr_ints=contr_ints,       &
                                 powers_bra=sub_shells(ishell)%powers,         &
                                 powers_ket=sub_shells(jshell)%powers)
          end if
          ! sets the minimum and maximum of indices of rows of the integral matrices
          min_row_idx = sub_shells(ishell)%base_idx+1
          max_row_idx = sub_shells(ishell)%base_idx &
                      + sub_shells(ishell)%num_ao*sub_shells(ishell)%num_contr
!FIXME: the returned index order of \var(num_prop)
          ! assigns the returned integrals, and write the integrals on file if required
          if (present(val_ints)) then
            ! returns integral matrices with redundant high order (>1) total geometric derivatives
            if (p_redunt_ints) then
              call GeomTreeNumRedunt(geom_tree, num_redunt_geo)
              allocate(redunt_list(2,num_redunt_geo), stat=ierr)
              if (ierr/=0) call QUIT("Gen1IntShellIntegral>> failed to allocate redunt_list!")
              call GeomTreeReduntList(geom_tree, redunt_list)
              iopt = 0
              do iprop = 1, num_prop
                do igeo = 1, num_redunt_geo
                  call MatSetBlockedValues(val_ints(iopt+redunt_list(2,igeo)), &
                                           min_row_idx, max_row_idx,           &
                                           min_col_idx, max_col_idx,           &
                                           contr_ints(:,:,:,:,iprop,           &
                                                      redunt_list(1,igeo)),    &
                                           .false.)
                  if (ishell/=jshell)                                            &
                    call MatSetBlockedValues(val_ints(iopt+redunt_list(2,igeo)), &
                                             min_row_idx, max_row_idx,           &
                                             min_col_idx, max_col_idx,           &
                                             -contr_ints(:,:,:,:,iprop,          &
                                                         redunt_list(1,igeo)),   &
                                             .true.)
                end do
                iopt = iopt+dim_redunt_geo
              end do
              deallocate(redunt_list)
            ! returns integrals matrices with unique total geometric derivatives
            else
              iopt = 0
              do iprop = 1, num_prop
                do igeo = 1, num_unique_geo
                  call MatSetBlockedValues(val_ints(iopt+base_unique_geo+igeo), &
                                           min_row_idx, max_row_idx,            &
                                           min_col_idx, max_col_idx,            &
                                           contr_ints(:,:,:,:,iprop,igeo), .false.)
                  if (ishell/=jshell)                                             &
                    call MatSetBlockedValues(val_ints(iopt+base_unique_geo+igeo), &
                                             min_row_idx, max_row_idx,            &
                                             min_col_idx, max_col_idx,            &
                                             -contr_ints(:,:,:,:,iprop,igeo), .true.)
                end do
                iopt = iopt+dim_unique_geo
              end do
            end if
          end if
          ! calculates the expectation values
          if (evaluate_expt) then
            ! returns expectation values with redundant high order (>1) total geometric derivatives
            if (p_redunt_expt) then
              do idens = 1, num_dens
                do igeo = 1, num_unique_geo
                  do iprop = 1, num_prop
                    call MatMultBlockedTrace(ao_dens(idens),                 &
                                             min_row_idx, max_row_idx,       &
                                             min_col_idx, max_col_idx,       &
                                             contr_ints(:,:,:,:,iprop,igeo), &
                                             unique_expt(igeo,iprop,idens), .false.)
                    if (ishell/=jshell)                                         &
                      call MatMultBlockedTrace(ao_dens(idens),                  &
                                               min_row_idx, max_row_idx,        &
                                               min_col_idx, max_col_idx,        &
                                               -contr_ints(:,:,:,:,iprop,igeo), &
                                               unique_expt(igeo,iprop,idens), .true.)
                  end do
                end do
              end do
            ! returns expectation values with unique total geometric derivatives
            else
              do idens = 1, num_dens
                do igeo = 1, num_unique_geo
                  do iprop = 1, num_prop
                    iopt = base_unique_geo+igeo+(iprop-1)*dim_unique_geo &
                         + (idens-1)*dim_unique_geo*num_prop
                    call MatMultBlockedTrace(ao_dens(idens),                 &
                                             min_row_idx, max_row_idx,       &
                                             min_col_idx, max_col_idx,       &
                                             contr_ints(:,:,:,:,iprop,igeo), &
                                             val_expt(iopt), .false.)
                    if (ishell/=jshell)                                         &
                      call MatMultBlockedTrace(ao_dens(idens),                  &
                                               min_row_idx, max_row_idx,        &
                                               min_col_idx, max_col_idx,        &
                                               -contr_ints(:,:,:,:,iprop,igeo), &
                                               val_expt(iopt), .true.)
                  end do
                end do
              end do
            end if
          end if
          deallocate(contr_ints)
        end do
      end do
    ! square integral matrices (column- or ket-major)
    case default
      ! loops over AO sub-shells on ket center
      do jshell = 1, num_shells
        ! sets the minimum and maximum of indices of columns of the integral matrices
        min_col_idx = sub_shells(jshell)%base_idx+1
        max_col_idx = sub_shells(jshell)%base_idx & 
                    + sub_shells(jshell)%num_ao*sub_shells(jshell)%num_contr
        ! loops over AO sub-shells on bra center (upper and diagonal parts)
        do ishell = 1, num_shells
          spher_gto = sub_shells(ishell)%spher_gto .or. sub_shells(jshell)%spher_gto
          allocate(contr_ints(sub_shells(ishell)%num_ao,    &
                              sub_shells(ishell)%num_contr, &
                              sub_shells(jshell)%num_ao,    &
                              sub_shells(jshell)%num_contr, &
                              num_prop,num_unique_geo), stat=ierr)
          if (ierr/=0) call QUIT("Gen1IntShellIntegral>> failed to allocate contr_ints!")
          ! spherical GTOs
          if (spher_gto) then
            ! calls Gen1Int subroutines to evaluate property integrals
            call OnePropEvaluate(idx_bra=sub_shells(ishell)%idx_cent,          &
                                 coord_bra=sub_shells(ishell)%coord_cent,      &
                                 angular_bra=sub_shells(ishell)%ang_num,       &
                                 num_prim_bra=sub_shells(ishell)%num_prim,     &
                                 exponent_bra=sub_shells(ishell)%exponents,    &
                                 num_contr_bra=sub_shells(ishell)%num_contr,   &
                                 contr_coef_bra=sub_shells(ishell)%contr_coef, &
                                 idx_ket=sub_shells(jshell)%idx_cent,          &
                                 coord_ket=sub_shells(jshell)%coord_cent,      &
                                 angular_ket=sub_shells(jshell)%ang_num,       &
                                 num_prim_ket=sub_shells(jshell)%num_prim,     &
                                 exponent_ket=sub_shells(jshell)%exponents,    &
                                 num_contr_ket=sub_shells(jshell)%num_contr,   &
                                 contr_coef_ket=sub_shells(jshell)%contr_coef, &
                                 spher_gto=spher_gto,                          &
                                 london_ao=london_ao,                          &
                                 one_prop=one_prop,                            &
                                 order_mag_bra=order_mag_bra,                  &
                                 order_mag_ket=order_mag_ket,                  &
                                 order_mag_total=order_mag_total,              &
                                 order_ram_bra=order_ram_bra,                  &
                                 order_ram_ket=order_ram_ket,                  &
                                 order_ram_total=order_ram_total,              &
                                 order_geo_bra=order_geo_bra,                  &
                                 order_geo_ket=order_geo_ket,                  &
                                 geom_tree=geom_tree,                          &
                                 num_gto_bra=sub_shells(ishell)%num_ao,        &
                                 num_gto_ket=sub_shells(jshell)%num_ao,        &
                                 num_opt=num_opt, contr_ints=contr_ints)
            ! reorders p-shell spherical GTOs
            if (sub_shells(ishell)%ang_num==1)                      &
              call reorder_p_sgto(1, sub_shells(ishell)%num_contr   &
                                     *sub_shells(jshell)%num_ao     &
                                     *sub_shells(jshell)%num_contr, &
                                  num_opt, contr_ints)
            if (sub_shells(jshell)%ang_num==1)                   &
              call reorder_p_sgto(sub_shells(ishell)%num_ao      &
                                  *sub_shells(ishell)%num_contr, &
                                  sub_shells(jshell)%num_contr,  &
                                  num_opt, contr_ints)
          ! Cartesian GTOs
          else
            ! calls Gen1Int subroutines to evaluate property integrals, and reorders
            ! the integrals according to Cartesian powers
            call OnePropEvaluate(idx_bra=sub_shells(ishell)%idx_cent,          &
                                 coord_bra=sub_shells(ishell)%coord_cent,      &
                                 angular_bra=sub_shells(ishell)%ang_num,       &
                                 num_prim_bra=sub_shells(ishell)%num_prim,     &
                                 exponent_bra=sub_shells(ishell)%exponents,    &
                                 num_contr_bra=sub_shells(ishell)%num_contr,   &
                                 contr_coef_bra=sub_shells(ishell)%contr_coef, &
                                 idx_ket=sub_shells(jshell)%idx_cent,          &
                                 coord_ket=sub_shells(jshell)%coord_cent,      &
                                 angular_ket=sub_shells(jshell)%ang_num,       &
                                 num_prim_ket=sub_shells(jshell)%num_prim,     &
                                 exponent_ket=sub_shells(jshell)%exponents,    &
                                 num_contr_ket=sub_shells(jshell)%num_contr,   &
                                 contr_coef_ket=sub_shells(jshell)%contr_coef, &
                                 spher_gto=spher_gto,                          &
                                 london_ao=london_ao,                          &
                                 one_prop=one_prop,                            &
                                 order_mag_bra=order_mag_bra,                  &
                                 order_mag_ket=order_mag_ket,                  &
                                 order_mag_total=order_mag_total,              &
                                 order_ram_bra=order_ram_bra,                  &
                                 order_ram_ket=order_ram_ket,                  &
                                 order_ram_total=order_ram_total,              &
                                 order_geo_bra=order_geo_bra,                  &
                                 order_geo_ket=order_geo_ket,                  &
                                 geom_tree=geom_tree,                          &
                                 num_gto_bra=sub_shells(ishell)%num_ao,        &
                                 num_gto_ket=sub_shells(jshell)%num_ao,        &
                                 num_opt=num_opt, contr_ints=contr_ints,       &
                                 powers_bra=sub_shells(ishell)%powers,         &
                                 powers_ket=sub_shells(jshell)%powers)
          end if
          ! sets the minimum and maximum of indices of rows of the integral matrices
          min_row_idx = sub_shells(ishell)%base_idx+1
          max_row_idx = sub_shells(ishell)%base_idx &
                      + sub_shells(ishell)%num_ao*sub_shells(ishell)%num_contr
!FIXME: the returned index order of \var(num_prop)
          ! assigns the returned integrals, and write the integrals on file if required
          if (present(val_ints)) then
            ! returns integral matrices with redundant high order (>1) total geometric derivatives
            if (p_redunt_ints) then
              call GeomTreeNumRedunt(geom_tree, num_redunt_geo)
              allocate(redunt_list(2,num_redunt_geo), stat=ierr)
              if (ierr/=0) call QUIT("Gen1IntShellIntegral>> failed to allocate redunt_list!")
              call GeomTreeReduntList(geom_tree, redunt_list)
              iopt = 0
              do iprop = 1, num_prop
                do igeo = 1, num_redunt_geo
                  call MatSetBlockedValues(val_ints(iopt+redunt_list(2,igeo)), &
                                           min_row_idx, max_row_idx,           &
                                           min_col_idx, max_col_idx,           &
                                           contr_ints(:,:,:,:,iprop,           &
                                                      redunt_list(1,igeo)),    &
                                           .false.)
                end do
                iopt = iopt+dim_redunt_geo
              end do
              deallocate(redunt_list)
            ! returns integrals matrices with unique total geometric derivatives
            else
              iopt = 0
              do iprop = 1, num_prop
                do igeo = 1, num_unique_geo
                  call MatSetBlockedValues(val_ints(iopt+base_unique_geo+igeo), &
                                           min_row_idx, max_row_idx,            &
                                           min_col_idx, max_col_idx,            &
                                           contr_ints(:,:,:,:,iprop,igeo), .false.)
                end do
                iopt = iopt+dim_unique_geo
              end do
            end if
          end if
          ! calculates the expectation values
          if (evaluate_expt) then
            ! returns expectation values with redundant high order (>1) total geometric derivatives
            if (p_redunt_expt) then
              do idens = 1, num_dens
                do igeo = 1, num_unique_geo
                  do iprop = 1, num_prop
                    call MatMultBlockedTrace(ao_dens(idens),                 &
                                             min_row_idx, max_row_idx,       &
                                             min_col_idx, max_col_idx,       &
                                             contr_ints(:,:,:,:,iprop,igeo), &
                                             unique_expt(igeo,iprop,idens), .false.)
                  end do
                end do
              end do
            ! returns expectation values with unique total geometric derivatives
            else
              do idens = 1, num_dens
                do igeo = 1, num_unique_geo
                  do iprop = 1, num_prop
                    iopt = base_unique_geo+igeo+(iprop-1)*dim_unique_geo &
                         + (idens-1)*dim_unique_geo*num_prop
                    call MatMultBlockedTrace(ao_dens(idens),                 &
                                             min_row_idx, max_row_idx,       &
                                             min_col_idx, max_col_idx,       &
                                             contr_ints(:,:,:,:,iprop,igeo), &
                                             val_expt(iopt), .false.)
                  end do
                end do
              end do
            end if
          end if
          deallocate(contr_ints)
        end do
      end do
    end select
    ! returns expectation values with redundant high order (>1) total geometric derivatives
    if (p_redunt_expt) then
      call GeomTreeReduntExpt(geom_tree, num_unique_geo, num_prop*num_dens, &
                              unique_expt, dim_redunt_geo, val_expt)
      deallocate(unique_expt)
    end if
  end subroutine Gen1IntShellIntegral

  !> \brief calculates molecular orbitals at grid points
  !> \author Bin Gao
  !> \date 2012-03-11
  !> \param num_shells is the number of AO sub-shells
  !> \param sub_shells are the AO sub-shells
  !> \param num_ao is the number of atomic orbitals
  !> \param num_mo is the number of molecular orbitals
  !> \param mo_coef contains the molecular orbital coefficients
  !> \param num_points is the number of grid points
  !> \param grid_points contains the coordinates of grid points
  !> \param london_ao indicates if using London atomic orbitals
  !> \param num_derv is the number of derivatives
  !> \param order_mag is the order of magnetic derivatives
  !> \param order_ram is the order of derivatives w.r.t. total rotational angular momentum
  !> \param order_geo is the order of geometric derivatives
  !> \return val_mo contains the value of molecular orbitals at grid points
  subroutine Gen1IntShellMO(num_shells, sub_shells, num_ao, num_mo, mo_coef, &
                            num_points, grid_points, num_derv, val_mo,       &
                            london_ao, order_mag, order_ram, order_geo)
    integer, intent(in) :: num_shells
    type(sub_shell_t), intent(in) :: sub_shells(*)
    integer, intent(in) :: num_ao
    integer, intent(in) :: num_mo
    real(REALK), intent(in) :: mo_coef(num_ao,num_mo)
    integer, intent(in) :: num_points
    real(REALK), intent(in) :: grid_points(3,num_points)
    integer, intent(in) :: num_derv
    real(REALK), intent(out) :: val_mo(num_points*num_derv,num_mo)
    logical, optional, intent(in) :: london_ao
    integer, optional, intent(in) :: order_mag
    integer, optional, intent(in) :: order_ram
    integer, optional, intent(in) :: order_geo
    integer p_order_mag  !order of magnetic derivatives
    integer p_order_ram  !order of derivatives w.r.t. total rotational angular momentum
    integer p_order_geo  !order of geometric derivatives
    integer num_opt      !number of operators including grid points and derivatives
    integer ishell       !incremental recorders over AO sub-shells
    logical spher_gto    !if spherical GTOs
    real(REALK), allocatable :: contr_value(:,:,:)  !contracted GTOs at grid points
    real(REALK), allocatable :: tmp_value(:,:,:)    !temporary results from Gen1Int
    integer imo, iopt, icontr, iao, jao             !incremental recorders
    integer offset_ao                               !offset of AOs
    integer ierr                                    !error information
    ! checks the number of AOs
    ierr = sub_shells(num_shells)%base_idx &
         + sub_shells(num_shells)%num_ao*sub_shells(num_shells)%num_contr
    if (num_ao/=ierr) stop "Gen1IntShellMO>> incorrect number of AOs!"
    ! gets the number of derivatives
    num_opt = 1
    if (present(order_mag)) then
      p_order_mag = order_mag
      num_opt = num_opt*(order_mag+1)*(order_mag+2)/2
    else
      p_order_mag = 0
    end if
    if (present(order_ram)) then
      p_order_ram = order_ram
      num_opt = num_opt*(order_ram+1)*(order_ram+2)/2
    else
      p_order_ram = 0
    end if
    if (present(order_geo)) then
      p_order_geo = order_geo
      num_opt = num_opt*(order_geo+1)*(order_geo+2)/2
    else
      p_order_geo = 0
    end if
    if (num_opt/=num_derv) stop "Gen1IntShellMO>> incorrect number of derivatives!"
    num_opt = num_opt*num_points
    ! initializes
    val_mo = 0.0_REALK
    ! loops over AO sub-shells
    do ishell = 1, num_shells
      spher_gto = sub_shells(ishell)%spher_gto
      allocate(contr_value(sub_shells(ishell)%num_ao,    &
                           sub_shells(ishell)%num_contr, &
                           num_opt), stat=ierr)
      if (ierr/=0) call QUIT("Gen1IntShellMO>> failed to allocate contr_value!")
      ! spherical GTOs
      if (spher_gto) then
        call contr_sgto_value(sub_shells(ishell)%coord_cent,        &
                              sub_shells(ishell)%ang_num,           &
                              sub_shells(ishell)%num_prim,          &
                              sub_shells(ishell)%exponents,         &
                              sub_shells(ishell)%num_contr,         &
                              sub_shells(ishell)%contr_coef,        &
                              p_order_geo, num_points, grid_points, &
                              sub_shells(ishell)%num_ao, num_derv,  &
                              contr_value)
        ! reorders p-shell spherical GTOs
        if (sub_shells(ishell)%ang_num==1)                     &
          call reorder_p_sgto(1, sub_shells(ishell)%num_contr, &
                              num_opt, contr_value)
      ! Cartesian GTOs
      else
        if (sub_shells(ishell)%ang_num==0) then
          call contr_cgto_value(sub_shells(ishell)%coord_cent,        &
                                sub_shells(ishell)%ang_num,           &
                                sub_shells(ishell)%num_prim,          &
                                sub_shells(ishell)%exponents,         &
                                sub_shells(ishell)%num_contr,         &
                                sub_shells(ishell)%contr_coef,        &
                                p_order_geo, num_points, grid_points, &
                                sub_shells(ishell)%num_ao, num_derv,  &
                                contr_value)
        else
          allocate(tmp_value(sub_shells(ishell)%num_ao,    &
                             sub_shells(ishell)%num_contr, &
                             num_opt), stat=ierr)
          if (ierr/=0) call QUIT("Gen1IntShellMO>> failed to allocate tmp_value!")
          call contr_cgto_value(sub_shells(ishell)%coord_cent,        &
                                sub_shells(ishell)%ang_num,           &
                                sub_shells(ishell)%num_prim,          &
                                sub_shells(ishell)%exponents,         &
                                sub_shells(ishell)%num_contr,         &
                                sub_shells(ishell)%contr_coef,        &
                                p_order_geo, num_points, grid_points, &
                                sub_shells(ishell)%num_ao, num_derv,  &
                                tmp_value)
          ! reorders the results according to Cartesian powers
          call reorder_cgtos(sub_shells(ishell)%ang_num, sub_shells(ishell)%num_ao,      &
                             sub_shells(ishell)%powers, 1, sub_shells(ishell)%num_contr, &
                             num_opt, tmp_value, contr_value)
          deallocate(tmp_value)
        end if
      end if
      ! sets the offset of AO
      offset_ao = sub_shells(ishell)%base_idx
      ! gets the value of MOs at grid points
      do imo = 1, num_mo
        do iopt = 1, num_opt
          jao = offset_ao
          do icontr = 1, sub_shells(ishell)%num_contr
            do iao = 1, sub_shells(ishell)%num_ao
              jao = jao+1
              val_mo(iopt,imo) = val_mo(iopt,imo) &
                               + mo_coef(jao,imo)*contr_value(iao,icontr,iopt)
            end do
          end do
        end do
      end do
      deallocate(contr_value)
    end do
  end subroutine Gen1IntShellMO

  !> \brief frees space taken by the AO sub-shells
  !> \author Bin Gao
  !> \date 2011-10-03
  !> \param num_shells is the number of AO sub-shells
  !> \param sub_shells are the AO sub-shells
  subroutine Gen1IntShellDestroy(num_shells, sub_shells)
    integer, intent(in) :: num_shells
    type(sub_shell_t), intent(inout) :: sub_shells(*)
    integer ishell  !incremental recorder
    do ishell = 1, num_shells
      deallocate(sub_shells(ishell)%exponents)
      deallocate(sub_shells(ishell)%contr_coef)
      if (sub_shells(ishell)%spher_gto) then
        !-deallocate(sub_shells(ishell)%mag_num)
      else
        deallocate(sub_shells(ishell)%powers)
      end if
    end do
  end subroutine Gen1IntShellDestroy

end module gen1int_shell
