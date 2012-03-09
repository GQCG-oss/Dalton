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
!...  This file contains the interface for Fortran 77 users.
!
!...  2012-01-10, Bin Gao
!...  * first version

#include "xkind.h"

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
  subroutine gen1int_ifc_init(num_comp, num_atom_type, KATOM, num_sym_atom, &
                              ang_numbers, NBLCK, KANG, num_cgto, KBLOCK,   &
                              num_prim, num_contr, KPRIM, exponents, ucontr_coefs)
    ! module of Dalton AO sub-shells
    use dalton_shell
    implicit none
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
    call QENTER("gen1int_ifc_init")
    call DaltonShellCreate(num_comp, num_atom_type, KATOM, num_sym_atom, &
                           ang_numbers, NBLCK, KANG, num_cgto, KBLOCK,   &
                           num_prim, num_contr, KPRIM, exponents, ucontr_coefs)
    call QEXIT("gen1int_ifc_init")
    return
  end subroutine gen1int_ifc_init

  !> \brief frees the space taken by Dalton AO sub-shells after all calculations
  !> \author Bin Gao
  !> \date 2011-10-02
  subroutine gen1int_ifc_clean
    ! module of Dalton AO sub-shells
    use dalton_shell
    implicit none
    call QENTER("gen1int_ifc_clean")
    call DaltonShellDestroy()
    call QEXIT("gen1int_ifc_clean")
    return
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
    ! module of Dalton AO sub-shells
    use dalton_shell
    implicit none
    integer, intent(in) :: len_work
    real(REALK), intent(inout) :: wrk_space(len_work)
    integer, intent(in) :: io_viewer
    integer, intent(in) :: level_print
    real(REALK), parameter :: ERR_THRSH = 10.0_REALK**(-8)    !threshold of error
    real(REALK), parameter :: RATIO_THRSH = 10.0_REALK**(-6)  !threshold of ratio to the referenced result
    integer, parameter :: NUM_TEST = 4                        !number of tests
    character*20, parameter :: PROP_NAME(NUM_TEST) = &        !labels of testing property integrals,
      (/INT_KIN_ENERGY, INT_OVERLAP, INT_POT_ENERGY, &        !see Gen1int library src/gen1int.F90
        INT_CART_MULTIPOLE/)
    character*8, parameter :: HERM_PROP(NUM_TEST) = &         !labels of property integrals,
      (/"KINENERG", "OVERLAP ", "POTENERG", "CARMOM  "/)      !see \fn(PR1IN1) in abacus/her1pro.F
    logical, parameter :: LONDON_AO(NUM_TEST) = &             !if using London atomic orbitals
      (/.false., .false., .false., .false./)
    integer, parameter :: ORDER_MOM(NUM_TEST) = &             !order of Cartesian multipole moments
      (/0, 0, 0, 1/)
    integer, parameter :: ORDER_GEO_TOTAL(NUM_TEST) = &       !order of total geometric derivatives
      (/0, 0, 0, 0/)
    integer, parameter :: MAX_NUM_CENT(NUM_TEST) = &          !maximum number of differentiated centers
      (/0, 0, 0, 0/)
    logical, parameter :: WRT_INTS = .false.                  !if writing integrals on file
    logical, parameter :: WRT_EXPT = .false.                  !if writing expectation values on file
    type(one_prop_t) prop_operator                !operator for property integrals
    integer num_ao                                !number of orbitals
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
    ! gets the number of atomic orbitals
    call DaltonShellGetNumAO(num_ao=num_ao)
    write(io_viewer,100) "number of orbitals", num_ao
    write(io_viewer,110) "threshold of error", ERR_THRSH
    write(io_viewer,110) "threshold of ratio to the referenced result", RATIO_THRSH
    ! sets artifical AO density matrices
    dim_sym_dens = num_ao*(num_ao+1)/2
    dim_sq_dens = num_ao*num_ao
    allocate(symmetry_dens(dim_sym_dens), stat=ierr)
    if (ierr/=0) call QUIT("gen1int_ifc_test>> failed to allocate symmetry_dens!")
    allocate(square_dens(num_ao,num_ao), stat=ierr)
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
      call OnePropCreate(prop_name=trim(PROP_NAME(itst)), &
                         one_prop=prop_operator,          &
                         info_prop=ierr,                  &
                         num_prop=num_prop,               &
                         kind_prop=kind_prop,             &
                         coord_nuclei=CORD(:,1:NUCDEP),   &
                         charge_nuclei=-CHARGE(1:NUCDEP), &
                         dipole_origin=DIPORG,            &
                         order_mom=ORDER_MOM(itst))
      if (ierr/=0) then
        write(io_viewer,999) "failed to creat "//trim(PROP_NAME(itst))
        test_failed = .true.
        cycle
      else
        test_failed = .false.
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
        size_int = num_ao*(num_ao+1)/2
        triangular = .true.
        symmetric = kind_prop==SYMM_INT_MAT
      case default
        size_int = num_ao*num_ao
        triangular = .false.
        symmetric = .false.
      end select
      strt_herm_int = end_herm_expt+1
      end_herm_int = strt_herm_int+size_int*num_opt_derv-1
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
        call MatCreate(A=val_ints(imat), num_row=num_ao, info_mat=ierr, &
                       triangular=triangular, symmetric=symmetric)
        if (ierr/=0) call QUIT("gen1int_ifc_test>> failed to creates integral matrices!")
      end do
      ! computes the integrals and/or expectation values
      wrk_space(1:num_expt) = 0.0_REALK
!FIXME: to evaluate expectaion values
      call DaltonShellEvaluate(one_prop=prop_operator, london_ao=LONDON_AO(itst), &
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
        end_herm_int = strt_herm_int+size_int-1
        call MatArrayAlmostEqual(A=val_ints(imat),                               &
                                 values=wrk_space(strt_herm_int:end_herm_int),   &
                                 io_viewer=io_viewer, almost_equal=almost_equal, &
                                 triangular=triangular, symmetric=symmetric,     &
                                 threshold=ERR_THRSH, ratio_thrsh=RATIO_THRSH)
        if (.not.almost_equal) test_failed = .true.
        call MatDestroy(A=val_ints(imat))
        strt_herm_int = end_herm_int+1
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
    call QEXIT("gen1int_ifc_test")
    return
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
