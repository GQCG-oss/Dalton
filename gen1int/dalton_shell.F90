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
!...  This file takes care the atomic orbital (AO) sub-shells in Dalton.
!
!...  2012-01-10, Bin Gao
!...  * first version

#include "xkind.h"

!> \brief module of Dalton AO sub-shells
!> \author Bin Gao
!> \date 2012-01-10
module dalton_shell

  ! Fortran 90 module of Gen1Int library
  use gen1int
  ! AO sub-shells
  use gen1int_shell
  ! matrix module
  use gen1int_matrix

  implicit none

  ! if Dalton AO sub-shells is created
  logical, private, save :: shells_created = .false.
  ! number of AO sub-shells from Dalton
  integer, private, save :: num_sub_shells = 0
  ! AO sub-shells from Dalton
  type(sub_shell_t), private, allocatable, save :: sub_shells(:)

  public :: DaltonShellCreate

  public :: DaltonShellView
  public :: DaltonShellGetNumAO
  public :: DaltonShellIntegral
  public :: DaltonShellMO
  public :: DaltonShellDestroy

  contains

  !> \brief initializes the AO sub-shells from Dalton, should be called before any calculation,
  !>        based on \fn(ORBPRO) subroutine, getting the unnormalized contraction coefficients
  !> \author Bin Gao
  !> \date 2010-12-06
  !> \param num_comp is the number of components
  !> \param num_atom_type is the number of atomic types
  !> \param num_sym_atom contains the number of symmetry independent centers of atomic types
  !> \param ang_numbers contains the angular momentum (1=s, 2=p, 3=d, ...)
  !> \param num_cgto contains the number of CGTOs in the AO blocks for an angular momentum
  !> \param num_prim contains the number of uncontracted functions
  !> \param num_contr contains the number of contracted functions
  !> \param exponents contains the exponents of primitive shells
  !> \param ucontr_coefs contains the unnormalized contraction coefficients
  subroutine DaltonShellCreate(num_comp, num_atom_type, KATOM, num_sym_atom, &
                               ang_numbers, NBLCK, KANG, num_cgto, KBLOCK,   &
                               num_prim, num_contr, KPRIM, exponents, ucontr_coefs)
    integer, intent(in) :: num_comp
    integer, intent(in) :: num_atom_type
    integer, intent(in) :: KATOM
    integer, intent(in) :: num_sym_atom(KATOM)
    integer, intent(in) :: ang_numbers(KATOM,num_comp)
    integer, intent(in) :: NBLCK(KATOM,num_comp)
    integer, intent(in) :: KANG
    integer, intent(in) :: num_cgto(KANG,KATOM,num_comp)
    integer, intent(in) :: KBLOCK
    integer, intent(in) :: num_prim(KBLOCK,num_comp)
    integer, intent(in) :: num_contr(KBLOCK,num_comp)
    integer, intent(in) :: KPRIM
    real(REALK), intent(in) :: exponents(KPRIM,KBLOCK,num_comp)
    real(REALK), intent(in) :: ucontr_coefs(KPRIM,KPRIM,KBLOCK,num_comp)
#include "mxcent.h"
#include "maxaqn.h"
#include "ccom.h"
#include "nuclei.h"
    integer icomp      !incremental recorder over components
    integer IDX_CENT   !index of symmetry independent center
    integer IDX_BLOCK  !
    integer ITYP       !incremental recorder over number of atomic types
    integer ICENT      !incremental recorder over number of symmetry independent centers
    integer IANG       !incremental recorder over angular momentum
    integer KBCH       !
    logical spher_gto  !if SGTOs
    integer ishell     !incremental recorder over AO sub-shells
    integer ang_num    !angular number
    real(REALK), allocatable :: contr_coef(:,:)  !contraction coefficients
    integer icontr, iprim                        !incremental recorder over contractions
    integer ierr                                 !error information
    ! backs if the AO sub-shells are created
    if (shells_created) return
    call QENTER("DaltonShellCreate")
    ! gets the number of AO sub-shells and kind of GTOs
    num_sub_shells = 0
    spher_gto = .false.
    ! loops over components of basis sets
    do icomp = 1, num_comp
      ! number of atomic types
      do ITYP = 1, num_atom_type
        ! number of symmetry independent centers of this type
        do ICENT = 1, num_sym_atom(ITYP)
          ! angular momentum 1=s, 2=p, 3=d, etc.
          do IANG = 1, ang_numbers(ITYP,icomp)
            if (num_cgto(IANG, ITYP, icomp) > 0) then
              ! radovan: basis does not have to start with s
              num_sub_shells = num_sub_shells+1
              if (SPH(IANG)) spher_gto = .true.
            end if
          end do
         end do
      end do
    end do
    ! initializes the AO sub-shells
    allocate(sub_shells(num_sub_shells), stat=ierr)
    if (ierr/=0) call QUIT("DaltonShellCreate>> failed to allocate sub_shells!")
    ishell = 0
    ! loops over components of basis sets
    do icomp = 1, num_comp
      IDX_BLOCK = 0
      IDX_CENT = 0
      ! number of atomic types
      do ITYP = 1, num_atom_type
        ! number of symmetry independent centers of this type
        do ICENT = 1, num_sym_atom(ITYP)
          IDX_CENT = IDX_CENT+1
          KBCH = IDX_BLOCK
          ! angular momentum 1=s, 2=p, 3=d, etc.
          do IANG = 1, ang_numbers(ITYP,icomp)

            ! radovan: basis does not have to start with s
            if (num_cgto(IANG, ITYP, icomp) < 1) cycle

            ! next block
            KBCH = KBCH+1
            ! gets the contraction coefficients
            allocate(contr_coef(num_contr(KBCH,icomp),num_prim(KBCH,icomp)), stat=ierr)
            if (ierr/=0) call QUIT("DaltonShellCreate>> failed to allocate contr_coef!")
            do iprim = 1, num_prim(KBCH,icomp)
              do icontr = 1, num_contr(KBCH,icomp)
                contr_coef(icontr,iprim) = ucontr_coefs(iprim,icontr,KBCH,icomp)
              end do
            end do
            ! normalizes the contraction coefficients
            ang_num = IANG-1
            !-if (SPH(IANG)) then
            if (spher_gto) then
              call norm_contr_sgto(ang_num, num_prim(KBCH,icomp),                &
                                   exponents(1:num_prim(KBCH,icomp),KBCH,icomp), &
                                   num_contr(KBCH,icomp), contr_coef)
            else
              call norm_contr_cgto(ang_num, num_prim(KBCH,icomp),                &
                                   exponents(1:num_prim(KBCH,icomp),KBCH,icomp), &
                                   num_contr(KBCH,icomp), contr_coef)
            end if
            !-ISTBNU(IDX_CENT)  !stabiliser: basic sym. op. that do not move center
            ishell = ishell+1
            if (ishell>1) then
              !-call Gen1IntShellCreate(spher_gto=SPH(IANG), idx_cent=IDX_CENT,  &
              call Gen1IntShellCreate(spher_gto=spher_gto, idx_cent=IDX_CENT, &
                     coord_cent=CORD(1:3,IDX_CENT), ang_num=ang_num,          &
                     num_prim=num_prim(KBCH,icomp),                           &
                     exponents=exponents(1:num_prim(KBCH,icomp),KBCH,icomp),  &
                     num_contr=num_contr(KBCH,icomp), contr_coef=contr_coef,  &
                     last_shell=sub_shells(ishell-1), sub_shell=sub_shells(ishell))
            ! sets the first AO sub-shell
            else
              !-call Gen1IntShellCreate(spher_gto=SPH(IANG), idx_cent=IDX_CENT,  &
              call Gen1IntShellCreate(spher_gto=spher_gto, idx_cent=IDX_CENT, &
                     coord_cent=CORD(1:3,IDX_CENT), ang_num=ang_num,          &
                     num_prim=num_prim(KBCH,icomp),                           &
                     exponents=exponents(1:num_prim(KBCH,icomp),KBCH,icomp),  &
                     num_contr=num_contr(KBCH,icomp), contr_coef=contr_coef,  &
                     sub_shell=sub_shells(ishell))
            end if
            deallocate(contr_coef)
            ! skips other CGTOs in this AO block for this angular momentum
            KBCH = KBCH+num_cgto(IANG,ITYP,icomp)-1
          end do
        end do
        IDX_BLOCK = IDX_BLOCK+NBLCK(ITYP,icomp)
      end do
    end do
    shells_created = .true.
    call QEXIT("DaltonShellCreate")
  end subroutine DaltonShellCreate

  !> \brief visualizes the information of Dalton AO sub-shells
  !> \author Bin Gao
  !> \date 2012-03-09
  !> \param io_viewer is the logical unit number of the viewer
  subroutine DaltonShellView(io_viewer)
    integer, intent(in) :: io_viewer
    call Gen1IntShellView(num_shells=num_sub_shells, sub_shells=sub_shells, &
                          io_viewer=io_viewer)
  end subroutine DaltonShellView

  !> \brief gets the number of atomic orbitals from the Dalton AO sub-shells
  !> \author Bin Gao
  !> \date 2012-03-09
  !> \return num_ao is the number of atomic orbitals
  subroutine DaltonShellGetNumAO(num_ao)
    integer, intent(out) :: num_ao
    integer idx_first  !index of the first orbital in the last AO sub-shells
    ! gets the number of atomic orbitals
    call Gen1IntShellGetRangeAO(sub_shell=sub_shells(num_sub_shells), &
                                idx_first=idx_first, idx_last=num_ao)
  end subroutine DaltonShellGetNumAO

  !> \brief evaluates the integral matrices and/or expectation values with Dalton AO sub-shells
  !> \author Bin Gao
  !> \date 2011-01-11
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
  !> \param io_viewer is the IO unit of standard output
  !> \param level_print is the level of print
  !> \return val_ints contains the integral matrices
  !> \return val_expt contains the expectation values
  !> \note \var(val_expt) should be zero by users before calculations
  subroutine DaltonShellIntegral(one_prop, london_ao,                                &
                                 order_mag_bra, order_mag_ket, order_mag_total,      &
                                 order_ram_bra, order_ram_ket, order_ram_total,      &
                                 order_geo_bra, order_geo_ket, geom_tree,            &
                                 num_ints, val_ints, redunt_ints, wrt_ints,          &
                                 num_dens, ao_dens, val_expt, redunt_expt, wrt_expt, &
                                 io_viewer, level_print)
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
    type(geom_tree_t), optional, intent(inout) :: geom_tree
    integer, intent(in) :: num_ints
    type(matrix), optional, intent(inout) :: val_ints(num_ints)
    logical, optional, intent(in) :: redunt_ints
    logical, optional, intent(in) :: wrt_ints
    integer, intent(in) :: num_dens
    type(matrix), optional, intent(in) :: ao_dens(num_dens)
    real(REALK), optional, intent(inout) :: val_expt(num_ints*num_dens)
    logical, optional, intent(in) :: redunt_expt
    logical, optional, intent(in) :: wrt_expt
    integer, intent(in) :: io_viewer
    integer, intent(in) :: level_print
    integer curr_path  !index of current path of n-ary tree
    integer num_paths  !total number of different paths of n-ary tree
    integer ipath      !incremental recorder over different paths
    ! checks if we have the Dalton AO sub-shells
    if (.not.shells_created) &
      stop "DaltonShellIntegral>> Dalton AO sub-shells are not created!"
    ! dumps the information of Dalton AO sub-shells
    if (level_print>=10) call DaltonShellView(io_viewer=io_viewer)
    if (level_print>=5) then
      ! dumps the information of one-electron property integrals
      call OnePropView(one_prop=one_prop, io_viewer=io_viewer)
      if (present(london_ao)) then
        if (london_ao) write(io_viewer,100) "London atomic orbitals used"
      end if
      ! dumps different derivatives
      if (present(order_mag_bra))   &
        write(io_viewer,100) "order of magnetic derivatives on bra center", order_mag_bra
      if (present(order_mag_ket))   &
        write(io_viewer,100) "order of magnetic derivatives on ket center", order_mag_ket
      if (present(order_mag_total)) &
        write(io_viewer,100) "order of total magnetic derivatives", order_mag_total
      if (present(order_ram_bra))   &
        write(io_viewer,100) "order of derivatives w.r.t. total rotational "// &
                             "angular momentum on bra center", order_ram_bra
      if (present(order_ram_ket))   &
        write(io_viewer,100) "order of derivatives w.r.t. total rotational "// &
                             "angular momentum on ket center", order_ram_ket
      if (present(order_ram_total)) &
        write(io_viewer,100) "order of total derivatives w.r.t. total rotational "// &
                             "angular momentum", order_ram_total
      if (present(order_geo_bra))   &
        write(io_viewer,100) "order of geometric derivatives with respect to bra center", &
                             order_geo_bra
      if (present(order_geo_ket))   &
        write(io_viewer,100) "order of geometric derivatives with respect to ket center", &
                             order_geo_ket
      if (present(geom_tree)) then
        write(io_viewer,100) "evalutes total geometric derivatives"
        call GeomTreeView(geom_tree=geom_tree, io_viewer=io_viewer)
      end if
      ! dumps arguments related to returning integral matrices
      if (present(val_ints)) then
        write(io_viewer,100) "integral matrices returned"
        if (present(redunt_ints)) then
          if (redunt_ints)       &
            write(io_viewer,100) &
              "integral matrices with redundant total geometric derivatives returned"
        end if
      end if
      if (present(wrt_ints)) then
        if (wrt_ints) then
          write(io_viewer,100) "integral matrices written on file"
          if (present(redunt_ints)) then
            if (redunt_ints)       &
              write(io_viewer,100) &
                "integral matrices with redundant total geometric derivatives written"
          end if
        end if
      end if
      ! dumps arguments related to expectation values
      if (present(ao_dens)) then
        if (present(val_expt)) then
          write(io_viewer,100) "number of AO density matrices", num_dens
          write(io_viewer,100) "expectation values returned"
          if (present(redunt_expt)) then
            if (redunt_expt)       &
              write(io_viewer,100) &
                "expectation values with redundant total geometric derivatives returned"
          end if
        end if
        if (present(wrt_expt)) then
          if (wrt_expt) then
            write(io_viewer,100) "number of AO density matrices", num_dens
            write(io_viewer,100) "expectation values written on file"
            if (present(redunt_expt)) then
              if (redunt_expt)       &
                write(io_viewer,100) &
                  "expectation values with redundant total geometric derivatives written"
            end if
          end if
        end if
      end if
    end if
    ! calculates total geometric derivatives
    if (present(geom_tree)) then
      ! gets the index and current path and total number of different paths
      call GeomTreeIdxPath(geom_tree=geom_tree, idx_path=curr_path)
      call GeomTreeNumPath(geom_tree=geom_tree, num_paths=num_paths)
      ! calculates the property integrals of current path
      call Gen1IntShellIntegral(num_shells=num_sub_shells,       &
                                sub_shells=sub_shells,           &
                                one_prop=one_prop,               &
                                london_ao=london_ao,             &
                                order_mag_bra=order_mag_bra,     &
                                order_mag_ket=order_mag_ket,     &
                                order_mag_total=order_mag_total, &
                                order_ram_bra=order_ram_bra,     &
                                order_ram_ket=order_ram_ket,     &
                                order_ram_total=order_ram_total, &
                                order_geo_bra=order_geo_bra,     &
                                order_geo_ket=order_geo_ket,     &
                                geom_tree=geom_tree,             &
                                num_ints=num_ints,               &
                                val_ints=val_ints,               &
                                redunt_ints=redunt_ints,         &
                                wrt_ints=wrt_ints,               &
                                num_dens=num_dens,               &
                                ao_dens=ao_dens,                 &
                                val_expt=val_expt,               &
                                redunt_expt=redunt_expt,         &
                                wrt_expt=wrt_expt)
      ! loops over other paths
      do ipath = curr_path+1, num_paths
        ! generates the differentiated centers and their orders
        call GeomTreeSearch(geom_tree=geom_tree)
        ! dumps the information of current path
        if (level_print>=20) &
          call GeomTreeView(geom_tree=geom_tree, io_viewer=io_viewer)
        ! calculates the property integrals of current path
        call Gen1IntShellIntegral(num_shells=num_sub_shells,       &
                                  sub_shells=sub_shells,           &
                                  one_prop=one_prop,               &
                                  london_ao=london_ao,             &
                                  order_mag_bra=order_mag_bra,     &
                                  order_mag_ket=order_mag_ket,     &
                                  order_mag_total=order_mag_total, &
                                  order_ram_bra=order_ram_bra,     &
                                  order_ram_ket=order_ram_ket,     &
                                  order_ram_total=order_ram_total, &
                                  order_geo_bra=order_geo_bra,     &
                                  order_geo_ket=order_geo_ket,     &
                                  geom_tree=geom_tree,             &
                                  num_ints=num_ints,               &
                                  val_ints=val_ints,               &
                                  redunt_ints=redunt_ints,         &
                                  wrt_ints=wrt_ints,               &
                                  num_dens=num_dens,               &
                                  ao_dens=ao_dens,                 &
                                  val_expt=val_expt,               &
                                  redunt_expt=redunt_expt,         &
                                  wrt_expt=wrt_expt)
      end do
    ! no total geometric derivatives
    else
      call Gen1IntShellIntegral(num_shells=num_sub_shells,       &
                                sub_shells=sub_shells,           &
                                one_prop=one_prop,               &
                                london_ao=london_ao,             &
                                order_mag_bra=order_mag_bra,     &
                                order_mag_ket=order_mag_ket,     &
                                order_mag_total=order_mag_total, &
                                order_ram_bra=order_ram_bra,     &
                                order_ram_ket=order_ram_ket,     &
                                order_ram_total=order_ram_total, &
                                order_geo_bra=order_geo_bra,     &
                                order_geo_ket=order_geo_ket,     &
                                num_ints=num_ints,               &
                                val_ints=val_ints,               &
                                redunt_ints=redunt_ints,         &
                                wrt_ints=wrt_ints,               &
                                num_dens=num_dens,               &
                                ao_dens=ao_dens,                 &
                                val_expt=val_expt,               &
                                redunt_expt=redunt_expt,         &
                                wrt_expt=wrt_expt)
    end if
100 format("DaltonShellIntegral>> ",A,I4)
  end subroutine DaltonShellIntegral

  !> \brief calculates molecular orbitals at grid points
  !> \author Bin Gao
  !> \date 2012-03-11
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
  subroutine DaltonShellMO(num_ao, num_mo, mo_coef, num_points, grid_points,  &
                           num_derv, val_mo, london_ao, order_mag, order_ram, &
                           order_geo)
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
    call Gen1IntShellMO(num_shells=num_sub_shells, sub_shells=sub_shells,      &
                        num_ao=num_ao, num_mo=num_mo, mo_coef=mo_coef,         &
                        num_points=num_points, grid_points=grid_points,        &
                        num_derv=num_derv, val_mo=val_mo, london_ao=london_ao, &
                        order_mag=order_mag, order_ram=order_ram, order_geo=order_geo)
  end subroutine DaltonShellMO

  !> \brief frees the space taken by Dalton AO sub-shells after all calculations
  !> \author Bin Gao
  !> \date 2011-10-02
  subroutine DaltonShellDestroy
    ! backs if the AO sub-shells are not created
    if (.not. shells_created) return
    call Gen1IntShellDestroy(num_sub_shells, sub_shells)
    deallocate(sub_shells)
    num_sub_shells = 0
    shells_created = .false.
  end subroutine DaltonShellDestroy

end module dalton_shell
