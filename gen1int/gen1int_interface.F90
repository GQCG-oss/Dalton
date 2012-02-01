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
!...  This file contains the interface of calling Gen1Int library.
!
!...  2012-01-10, Bin Gao
!...  * first version

#include "xkind.h"

!> \brief interface of calling Gen1Int library
!> \author Bin Gao
!> \date 2012-01-10
module gen1int_interface

  ! Fortran 90 module of Gen1Int library
  use gen1int
  ! AO sub-shells
  use gen1int_shell
  ! matrix module
  use gen1int_matrix

  implicit none

  ! if performing test suite
  logical, public, save :: test_gen1int_interface = .false.
  ! if this interface are initialized
  logical, private, save :: init_interface = .false.
  ! number of AO sub-shells from Dalton
  integer, private, save :: num_sub_shells = 0
  ! AO sub-shells from Dalton
  type(sub_shell_t), private, allocatable, save :: sub_shells(:)

  public :: gen1int_ifc_init
  public :: gen1int_ifc_main
  public :: gen1int_ifc_clean
  public :: gen1int_ifc_test

  contains

  !> \brief initializes the data used in Gen1Int interface before any calculation,
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
  subroutine gen1int_ifc_init(num_comp, num_atom_type, KATOM, num_sym_atom, &
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
#ifdef GEN1INT_DEBUG
#include "priunit.h"
#endif
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
    ! backs if already initialized the basis sets
    if (init_interface) return
    call QENTER("gen1int_ifc_init")
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
            num_sub_shells = num_sub_shells+1
            if (SPH(IANG)) spher_gto = .true.
          end do
         end do
      end do
    end do
    ! initializes the AO sub-shells
    allocate(sub_shells(num_sub_shells), stat=ierr)
    if (ierr/=0) call QUIT("gen1int_ifc_init>> failed to allocate sub_shells!")
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
            ! next block
            KBCH = KBCH+1
            ! gets the contraction coefficients
            allocate(contr_coef(num_contr(KBCH,icomp),num_prim(KBCH,icomp)), stat=ierr)
            if (ierr/=0) call QUIT("gen1int_ifc_init>> failed to allocate contr_coef!")
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
#ifdef GEN1INT_DEBUG
    ! visualizes the information of AO sub-shells
    call Gen1IntShellView(num_sub_shells, sub_shells, LUPRI)
#endif
    init_interface = .true.
    call QEXIT("gen1int_ifc_init")
  end subroutine gen1int_ifc_init

  !> \brief main driver of calling Gen1Int interface
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
  !> \param order_geo_total is the order of total geometric derivatives
  !> \param max_num_cent is the maximum number of differentiated centers
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
  subroutine gen1int_ifc_main(one_prop, london_ao,             &
                              order_mag_bra, order_mag_ket,    &
                              order_mag_total,                 &
                              order_ram_bra, order_ram_ket,    &
                              order_ram_total,                 &
                              order_geo_bra, order_geo_ket,    &
                              order_geo_total, max_num_cent,   &
                              num_ints, val_ints, redunt_ints, &
                              wrt_ints, num_dens, ao_dens,     &
                              val_expt, redunt_expt, wrt_expt, &
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
    integer, optional, intent(in) :: order_geo_total
    integer, optional, intent(in) :: max_num_cent
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
    integer p_max_ncent          !maximum number of differentiated centers (private)
    type(geom_tree_t) geom_tree  !n-ary tree for total geometric derivatives
    integer num_paths            !total number of different paths
    integer info_ifc             !information of the interface
    integer ipath                !incremental recorder over different paths
    ! Dalton stuff
#include "mxcent.h"
#include "nuclei.h"
    real(REALK) TIMHER, TEND
    real(REALK) WALHER, WEND
    call QENTER("gen1int_ifc_main")
    call GETTIM(TIMHER, WALHER)
    ! checks if the interface is initialized
    if (.not.init_interface) &
      call QUIT("gen1int_ifc_main>> Gen1Int interface is not initialized!")
    ! dumps to check
    if (level_print>=10) then
      ! visualizes the information of AO sub-shells
      call Gen1IntShellView(num_sub_shells, sub_shells, io_viewer)
    end if
    if (level_print>=5) then
      call OnePropView(one_prop, io_viewer)
      if (present(london_ao)) then
        if (london_ao) write(io_viewer,100) "London atomic orbitals used"
      end if
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
      if (present(order_geo_total) .and. present(max_num_cent)) then
        write(io_viewer,100) "evalutes total geometric derivatives"
        write(io_viewer,100) "maximum number of differentiated centers", max_num_cent
        write(io_viewer,100) "order of total geometric derivatives", order_geo_total
      end if
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
    if (present(order_geo_total)) then
      if (present(max_num_cent)) then
        p_max_ncent = max_num_cent
      else
        p_max_ncent = order_geo_total
      end if
      ! computes the total number of different paths, and generates
      ! the first path (composition of centers)
      call GeomTreeCreate(num_atoms=NUCDEP, order_geo=order_geo_total, &
                          max_ncent=p_max_ncent, geom_tree=geom_tree,  &
                          num_paths=num_paths, info_geom=info_ifc)
      if (info_ifc/=0) &
        call QUIT("gen1int_ifc_main>> error occurred when calling GeomTreeCreate!")
      ! dumps to check
      if (level_print>=10) call GeomTreeView(geom_tree, io_viewer)
      ! calculates the property integrals
      call Gen1IntShellEvaluate(num_shells=num_sub_shells,       &
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
      do ipath = 2, num_paths
        ! generates the differentiated centers and their orders
        call GeomTreeSearch(geom_tree=geom_tree)
        ! dumps to check
        if (level_print>=20) call GeomTreeView(geom_tree, io_viewer)
        ! calculates the property integrals
        call Gen1IntShellEvaluate(num_shells=num_sub_shells,       &
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
      ! frees space taken by a n-ary tree for total geometric derivatives
      call GeomTreeDestroy(geom_tree)
    ! no total geometric derivatives
    else
      call Gen1IntShellEvaluate(num_shells=num_sub_shells,       &
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
    ! Dalton stuff
    call GETTIM(TEND, WEND)
    TIMHER = TEND-TIMHER
    WALHER = WEND-WALHER
    write(io_viewer,"()")
    call TIMTXT(">>>> Total CPU  time used in Gen1Int", TIMHER, io_viewer)
    call TIMTXT(">>>> Total wall time used in Gen1Int", WALHER, io_viewer)
    call HEADER("End of Gen1Int", -1)
    call QEXIT("gen1int_ifc_main")
100 format("gen1int_ifc_main>> ",A,I4)
  end subroutine gen1int_ifc_main

  !> \brief cleans up the data used in Gen1Int interface after all calculations
  !> \author Bin Gao
  !> \date 2011-10-02
  subroutine gen1int_ifc_clean
    if (.not. init_interface) return
    call QENTER("gen1int_ifc_clean")
    ! frees space taken by the AO sub-shells
    call Gen1IntShellDestroy(num_sub_shells, sub_shells)
    deallocate(sub_shells)
    num_sub_shells = 0
    init_interface = .false.
    call QEXIT("gen1int_ifc_clean")
  end subroutine gen1int_ifc_clean

  !> \brief test suite of Gen1Int interface, enabled by adding the following lines
  !>        in DALTON.INP/DIRAC.INP
  !>        **INTEGRAL
  !>        .GENINT
  !> \author Bin Gao
  !> \date 2011-10-04
  !> \param len_work is the length of Dalton/Dirac workspace
  !> \param wrk_space is the Dalton/Dirac workspace
  !> \param io_viewer is the IO unit of standard output
  !> \param level_print is the level of print
  subroutine gen1int_ifc_test(len_work, wrk_space, io_viewer, level_print)
    integer, intent(in) :: len_work
    real(REALK), intent(inout) :: wrk_space(len_work)
    integer, intent(in) :: io_viewer
    integer, intent(in) :: level_print
    real(REALK), parameter :: ERR_THRSH = 10.0_REALK**(-8)    !threshold of error
    real(REALK), parameter :: RATIO_THRSH = 10.0_REALK**(-6)  !threshold of ratio to the referenced result
    integer, parameter :: NUM_TEST = 3                        !number of tests
    character*14, parameter :: PROP_NAME(NUM_TEST) = &        !labels of testing property integrals,
      (/INT_KIN_ENERGY, INT_OVERLAP, INT_POT_ENERGY/)         !see Gen1int library src/gen1int.F90
    character*8, parameter :: HERM_PROP(NUM_TEST) = &         !labels of property integrals,
      (/"KINENERG", "OVERLAP ", "POTENERG"/)                  !see \fn(PR1IN1) in abacus/her1pro.F
    logical, parameter :: LONDON_AO(NUM_TEST) = &             !if using London atomic orbitals
      (/.false., .false., .false./)
    integer, parameter :: ORDER_MOM(NUM_TEST) = &             !order of Cartesian multipole moments
      (/0, 0, 0/)
    integer, parameter :: ORDER_GEO_TOTAL(NUM_TEST) = &       !order of total geometric derivatives
      (/0, 0, 0/)
    integer, parameter :: MAX_NUM_CENT(NUM_TEST) = &          !maximum number of differentiated centers
      (/0, 0, 0/)
    logical, parameter :: WRT_INTS = .false.                  !if writing integrals on file
    logical, parameter :: WRT_EXPT = .false.                  !if writing expectation values on file
    type(one_prop_t) prop_operator                !operator for property integrals
    integer num_orb                               !number of orbitals
    integer, parameter :: NUM_DENS = 1            !number of AO density matrices
    integer dim_sym_dens, dim_sq_dens             !dimensions of AO density matrices
    real(REALK), allocatable :: symmetry_dens(:)  !artifical AO density matrices
    real(REALK), allocatable :: square_dens(:,:)
    integer ierr                                  !error information
    integer itst                                  !incremental recorder of tests
! origins in Dalton
#include "orgcom.h"
! uses \var(MXCENT)
#include "mxcent.h"
! number of atomic centers
#include "nuclei.h"
! uses \var(MXCORB)
#include "maxorb.h"
! uses \var(MXQN)
#include "maxaqn.h"
! uses \var(MAXREP)
#include "symmet.h"
    integer num_prop                              !number of property integral matrices
    integer kind_prop                             !kind of property integral matrices
    logical triangular                            !integral matrices are triangularized or squared
    logical symmetric                             !integral matrices are symmetric or anti-symmetric
    integer dim_unique_geo                        !dimension of all unique total geometric derivatives
    integer num_opt_derv                          !number of operators including derivatives
    type(matrix), allocatable :: val_ints(:)      !integral matrices
    integer imat                                  !incremental recorder over matrices
    integer num_expt                              !number of expectation values
    integer strt_herm_expt, end_herm_expt         !addresses for expectation values from \fn(PR1IN1)
    integer size_int                              !size of property integrals
    integer strt_herm_int, end_herm_int           !addresses for integrals from \fn(PR1IN1)
    ! variables related to \fn(PR1IN1) ...
!FIXME: having problem of calling \fn(PR1IN1) with \var(GET_EXPT)=.true., which gives wrong integrals
    logical, parameter :: GET_EXPT = .false.      !if getting expectation values back
    integer max_typ
    integer, allocatable :: int_rep(:)
    integer lint_ad
    integer, allocatable :: int_adr(:)
    character*8, allocatable :: lb_int(:)
    integer NCOMP
    logical, parameter :: TOFILE = .false.        !if writing integrals on file
    character*6, parameter :: MTFORM = "TRIANG"
    logical, parameter :: DOINT(4) = (/.true., .false., .false., .false./)
    logical, parameter :: PROP_PRINT = .false.    !if printing referenced property integrals
    integer, parameter :: NUM_PQUAD = 40          !number of integration points for DSO integrals
    integer len_free                              !length of free Dalton/Dirac workspace
#if !defined (PRG_DIRAC)
    integer base_free                             !base address of free Dalton/Dirac workspace
#endif
    logical almost_equal                          !indicates if the results from Gen1Int are almost equal
                                                  !to those from \fn(PR1IN1)
    logical test_failed                           !indicator if the test failed
    call QENTER("gen1int_ifc_test")
#if !defined (BUILD_OPENRSP)
    ! test suite of matrix module
    call MatTestSuite(test_failed=test_failed, io_viewer=io_viewer, &
                      level_print=level_print, threshold=ERR_THRSH)
#else
    test_failed = .false.
#endif
    ! checks if the interface is initialized
    if (.not.init_interface) &
      call QUIT("gen1int_ifc_test>> Gen1Int interface is not initialized!")
    ! gets the number of atomic orbitals
    call Gen1IntShellIdx(sub_shell=sub_shells(num_sub_shells), &
                         idx_first=ierr, idx_last=num_orb)
    write(io_viewer,100) "number of orbitals", num_orb
    write(io_viewer,110) "threshold of error", ERR_THRSH
    write(io_viewer,110) "threshold of ratio to the referenced result", RATIO_THRSH
    ! sets artifical AO density matrices
    dim_sym_dens = num_orb*(num_orb+1)/2
    dim_sq_dens = num_orb*num_orb
    allocate(symmetry_dens(dim_sym_dens), stat=ierr)
    if (ierr/=0) call QUIT("gen1int_ifc_test>> failed to allocate symmetry_dens!")
    allocate(square_dens(num_orb,num_orb), stat=ierr)
    if (ierr/=0) call QUIT("gen1int_ifc_test>> failed to allocate square_dens!")
    call random_number(symmetry_dens)
    call random_number(square_dens)
    symmetry_dens = 0.1_REALK*symmetry_dens
    square_dens = 0.1_REALK*square_dens
    ! loops over different tests
    do itst = 1, NUM_TEST
#if defined(PRG_DIRAC)
      if (HERM_PROP(itst)=="POTENERG") cycle
#endif
      ! initializes the information of one-electron property integrals
      call OnePropCreate(prop_name=PROP_NAME(itst), one_prop=prop_operator,      &
                         info_prop=ierr, num_prop=num_prop, kind_prop=kind_prop, &
                         dipole_origin=DIPORG, coord_nuclei=CORD(:,1:NUCDEP),    &
                         charge_nuclei=-CHARGE(1:NUCDEP))
      if (ierr/=0) then
        write(io_viewer,999) "failed to creat "//trim(PROP_NAME(itst))
        test_failed = .true.
        cycle
      end if
!FIXME: creates n-ary tree here
      dim_unique_geo = (3*NUCDEP)**ORDER_GEO_TOTAL(itst)
      num_opt_derv = dim_unique_geo*num_prop
      ! number and addresses of expectation values
      num_expt = num_opt_derv*NUM_DENS
      strt_herm_expt = num_expt+1
      end_herm_expt = 2*num_expt
      ! size and addresses of referenced integrals
      select case(kind_prop)
      case(SYMM_INT_MAT,ANTI_INT_MAT)
        size_int = num_orb*(num_orb+1)*num_opt_derv/2
        triangular = .true.
        symmetric = kind_prop==SYMM_INT_MAT
      case default
        size_int = num_orb*num_orb*num_opt_derv
        triangular = .false.
        symmetric = .false.
      end select
      strt_herm_int = end_herm_expt+1
      end_herm_int = strt_herm_int+size_int-1
      if (end_herm_int>len_work) then
        write(io_viewer,100) "test of "//trim(PROP_NAME(itst))//".Gen1Int failed!"
        write(io_viewer,999) "required memory", end_herm_int
        write(io_viewer,999) "available memory", len_work
        call QUIT("gen1int_ifc_test>> increase workspace!")
      end if
      ! allocates integral matrices
      allocate(val_ints(num_opt_derv), stat=ierr)
      if (ierr/=0) call QUIT("gen1int_ifc_test>> failed to allocate val_ints!")
      do imat = 1, num_opt_derv
        call MatCreate(A=val_ints(imat), num_row=num_orb, info_mat=ierr, &
                       triangular=triangular, symmetric=symmetric)
        if (ierr/=0) call QUIT("gen1int_ifc_test>> failed to creates integral matrices!")
      end do
      ! computes the integrals and/or expectation values
      wrk_space(1:num_expt) = 0.0_REALK
!FIXME: to evaluate expectaion values
      call gen1int_ifc_main(one_prop=prop_operator, london_ao=LONDON_AO(itst), &
                            order_geo_total=ORDER_GEO_TOTAL(itst),             &
                            max_num_cent=MAX_NUM_CENT(itst),                   &
                            num_ints=num_opt_derv, val_ints=val_ints,          &
                            wrt_ints=WRT_INTS, num_dens=NUM_DENS,              &
                            io_viewer=io_viewer, level_print=level_print)
      ! frees space taken by the information of one-electron property integrals
      call OnePropDestroy(one_prop=prop_operator)
      ! gets the referenced results from HERMIT
!FIXME: \var(FORQM3)
      if (trim(PROP_NAME(itst))=="DSO") then
        max_typ = (3*NUCDEP)**2
      else
        max_typ = 3*MXCOOR
      end if
      allocate(int_rep(max_typ), stat=ierr)
      int_rep = 0
      if (ierr/=0) call QUIT("gen1int_ifc_test>> failed to allocate int_rep!")
      if (trim(PROP_NAME(itst))=="ELFGRDC" .or. trim(PROP_NAME(itst))=="ELFGRDS") then
         lint_ad = 9*NUCIND*(MAXREP+1)
      else
         lint_ad = max_typ
      end if
      allocate(int_adr(lint_ad), stat=ierr)
      int_adr = 0
      if (ierr/=0) call QUIT("gen1int_ifc_test>> failed to allocate int_adr!")
      allocate(lb_int(max_typ), stat=ierr)
      if (ierr/=0) call QUIT("gen1int_ifc_test>> failed to allocate lb_int!")
#if !defined (PRG_DIRAC)
      base_free = 1
#endif
      len_free = len_work-end_herm_int
      write(io_viewer,100) "gets the referenced results from HERMIT ..."
      ! not equals to 0 so that \fn(PR1IN1) will copy the results back when first calling it
      NCOMP = -1
      wrk_space(strt_herm_expt:end_herm_expt) = 0.0_REALK
      if (triangular) then
#if defined(PRG_DIRAC)
        call PR1IN1(wrk_space(end_herm_int+1:), len_free, int_rep, int_adr,   &
                    lb_int, HERM_PROP(itst)(1:7), ORDER_MOM(itst), NUM_PQUAD, &
                    triangular, PROP_PRINT, level_print,                      &
                    wrk_space(strt_herm_int:end_herm_int), NCOMP, TOFILE,     &
                    MTFORM, DOINT)
#else
        call PR1IN1(wrk_space(end_herm_int+1:), base_free, len_free, int_rep, &
                    int_adr, lb_int, HERM_PROP(itst)(1:7), ORDER_MOM(itst),   &
                    NUM_PQUAD, triangular, PROP_PRINT, level_print,           &
                    wrk_space(strt_herm_int:end_herm_int),                    &
                    NCOMP, TOFILE, MTFORM, DOINT,                             &
                    wrk_space(strt_herm_expt:end_herm_expt), GET_EXPT, symmetry_dens)
#endif
      else
#if defined(PRG_DIRAC)
        call PR1IN1(wrk_space(end_herm_int+1:), len_free, int_rep, int_adr,   &
                    lb_int, HERM_PROP(itst)(1:7), ORDER_MOM(itst), NUM_PQUAD, &
                    triangular, PROP_PRINT, level_print,                      &
                    wrk_space(strt_herm_int:end_herm_int), NCOMP, TOFILE,     &
                    MTFORM, DOINT)
#else
        call PR1IN1(wrk_space(end_herm_int+1:), base_free, len_free, int_rep, &
                    int_adr, lb_int, HERM_PROP(itst)(1:7), ORDER_MOM(itst),   &
                    NUM_PQUAD, triangular, PROP_PRINT, level_print,           &
                    wrk_space(strt_herm_int:end_herm_int), NCOMP, TOFILE,     &
                    MTFORM, DOINT, wrk_space(strt_herm_expt:end_herm_expt),   &
                    GET_EXPT, square_dens)
#endif
      end if
      deallocate(int_rep)
      deallocate(int_adr)
      deallocate(lb_int)
      ! HERMIT uses different sign for the following integrals
      if (HERM_PROP(itst)=="DPLGRA  " .or. HERM_PROP(itst)=="POTENERG" .or. &
          HERM_PROP(itst)=="NUCSLO  " .or. HERM_PROP(itst)=="PSO     ") then
        wrk_space(strt_herm_int:end_herm_int) = -wrk_space(strt_herm_int:end_herm_int)
        wrk_space(strt_herm_expt:end_herm_expt) = -wrk_space(strt_herm_expt:end_herm_expt)
      end if
      ! checks the results
      write(io_viewer,100) "checks the results of "//trim(PROP_NAME(itst))
      do imat = 1, num_opt_derv
!FIXME for several integral matrices, values=
        call MatArrayAlmostEqual(A=val_ints(imat),                               &
                                 values=wrk_space(strt_herm_int:end_herm_int),   &
                                 io_viewer=io_viewer, almost_equal=almost_equal, &
                                 triangular=triangular, symmetric=symmetric,     &
                                 threshold=ERR_THRSH, ratio_thrsh=RATIO_THRSH)
        if (.not.almost_equal) test_failed = .true.
        call MatDestroy(A=val_ints(imat))
      end do
      deallocate(val_ints)
#if !defined (PRG_DIRAC)
!      if (GET_EXPT) then
!        do ierr = 1, num_expt
!          ! we check the ratio for absolutely greater values
!          if (abs(wrk_space(num_expt+ierr))>ERR_THRSH) then
!            ratio_to_ref = wrk_space(ierr)/wrk_space(num_expt+ierr)
!            if (ratio_to_ref<RATIO_THRSH(1) .or. ratio_to_ref>RATIO_THRSH(2)) then
!              write(io_viewer,997) trim(PROP_NAME(itst)), ierr, wrk_space(num_expt+ierr), &
!                                   wrk_space(ierr)
!              test_failed = .true.
!            end if
!          ! checks the difference for smaller values
!          else
!            if (abs(wrk_space(ierr)-wrk_space(num_expt+ierr))>ERR_THRSH) then
!              write(io_viewer,997) trim(PROP_NAME(itst)), ierr, wrk_space(num_expt+ierr), &
!                                   wrk_space(ierr)
!              test_failed = .true.
!            end if
!          end if
!        end do
!      end if
#endif
      if (test_failed) then
        write(io_viewer,100) "test of "//trim(PROP_NAME(itst))//".Gen1Int failed!"
      else
        write(io_viewer,100) "test of "//trim(PROP_NAME(itst))//".Gen1Int passed!"
      end if
    end do
    deallocate(symmetry_dens)
    deallocate(square_dens)
    ! cleans up the data in Gen1Int interface after all calculations will be
    ! performed in abacus/dalton.F (Dalton), main/dirac.F(Dirac)
    call QEXIT("gen1int_ifc_test")
100 format("gen1int_ifc_test>> ",A,I8)
110 format("gen1int_ifc_test>> ",A,Es16.6)
!995 format("gen1int_ifc_test>> ",A,".EXP(",I6,")>> SQUARE>>",F18.12,", SYM>>",F18.12)
!996 format("gen1int_ifc_test>> ",A,".EXP(",I6,")>> SQUARE>>",F18.12,", ANTI>>",F18.12)
#if !defined (PRG_DIRAC)
997 format("gen1int_ifc_test>> ",A,".EXP(",I6,")>> HERMIT>>",F18.12,", Gen1Int>>",F18.12)
#endif
!998 format("gen1int_ifc_test>> ",A,".INT(",I6,",",I6,")>> HERMIT>>",F18.12, &
!           ", Gen1Int>>",F18.12)
999 format("gen1int_ifc_test>> ",A,I16)
  end subroutine gen1int_ifc_test

end module gen1int_interface
