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
!...  This file is the test suite of Gen1Int interface.
!
!...  2011-10-03, Bin Gao
!...  * first version

#include "xkind.h"

  !> \brief test suite of Gen1Int interface, enabled by adding the following lines
  !>        in DALTON.INP
  !>        **INTEGRAL
  !>        .GENINT
  !> \author Bin Gao
  !> \date 2011-10-04
  !> \param len_work is the length of Dalton workspace
  !> \param dal_work is the Dalton workspace
  !> \param io_std is the IO unit of standard output
  !> \param level_print is the level of print
  subroutine gen1int_ifc_test(len_work, dal_work, io_std, level_print)
    ! AO sub-shells
    use gen1int_shell
    implicit none
    integer, intent(in) :: len_work
    real(REALK), intent(inout) :: dal_work(len_work)
    integer, intent(in) :: io_std
    integer, intent(in) :: level_print
#ifdef BUILD_GEN1INT
! threshold of error
#include "err_thrsh.h"
    ! number of expectation value tests
    integer, parameter :: NUM_EXP_TEST = 3
    ! property integral names for testing expectation values, see subroutine
    ! \fn(gen1int_shell_prop) in gen1int_shell.F90
    character*8, parameter :: PROP_EXPT(NUM_EXP_TEST) = &
      (/"ANGLON  ", "KINENERG", "PSO     "/)
    ! number of tests
    integer, parameter :: NUM_TEST = 8
    ! property integral names, see subroutine \fn(gen1int_shell_prop) in gen1int_shell.F90
    character*8, parameter :: PROP_NAME(NUM_TEST) =                 &
      (/"ANGLON  ", "DIPLEN  ", "KINENERG", "NSTLON  ", "NSTNOL  ", &
        "OVERLAP ", "POTENERG", "PSO     "/)
    ! property integral names in subroutine \fn(PR1IN1) in abacus/her1pro.F
    character*8, parameter :: DAL_PROP(NUM_TEST) =                  &
      (/"ANGLON  ", "DPLGRA  ", "KINENERG", "NUCSLO  ", "NUCSNLO ", &
        "OVERLAP ", "POTENERG", "PSO     "/)
    ! if London atomic orbitals
    logical, parameter :: IS_LAO(NUM_TEST) = &
      (/.false., .false., .false., .false., .false., .false., .false., .false./)
    ! order of Cartesian multipole moments
    integer, parameter :: ORDER_MOM(NUM_TEST) = &
      (/0, 0, 0, 0, 0, 0, 0, 0/)
    ! maximum number of differentiated centers
    integer, parameter :: MAX_NUM_CENT(NUM_TEST) = &
      (/0, 3, 0, 0, 0, 0, 0, 0/)
    ! order of total geometric derivatives
    integer, parameter :: ORDER_GEO_TOT(NUM_TEST) = &
      (/0, 1, 0, 0, 0, 0, 0, 0/)
    ! if getting integrals back from Gen1Int
    logical, parameter :: GET_INT = .true.
    ! if writing integrals on file
    logical, parameter :: WRT_INT = .false.
!FIXME: having problem of calling \fn(PR1IN1) with \var(GET_EXPT)=.true., which gives wrong integrals
    ! if getting expectation values back
    logical, parameter :: GET_EXPT = .false.
    ! if writing expectation values on file
    logical, parameter :: WRT_EXPT = .false.
    ! number of operators including derivatives
    integer num_opt_derv
    ! kind of integral matrices
    integer kind_int
    ! integrals are triangularized or squared
    logical is_triang
    ! dimensions of AO density matrices
    integer dim_tri_dens, dim_sq_dens
    ! number of AO density matrices
    integer, parameter :: NUM_DENS = 1
    ! AO density matrices
    real(REALK), allocatable :: sym_ao_dens(:)
    real(REALK), allocatable :: anti_ao_dens(:)
    real(REALK), allocatable :: sq_sym_dens(:,:)
    real(REALK), allocatable :: sq_anti_dens(:,:)
    real(REALK), allocatable :: sq_gen_dens(:,:)
    ! variables related to \fn(PR1IN1) ...
    integer max_typ
    integer, allocatable :: int_rep(:)
    !
    integer lint_ad
    integer, allocatable :: int_adr(:)
    !
    character*8, allocatable :: lb_int(:)
    !
    integer :: NCOMP = 0
    ! if writing integrals on file
    logical, parameter :: TOFILE = .false.
    !
    character*6, parameter :: MTFORM = "TRIANG"
    !
    logical, parameter :: DOINT(4) = (/.true., .false., .false., .false./)
    ! if printing referenced property integrals
    logical, parameter :: PROP_PRINT = .false.
    ! number of integration points for diamagnetic spin-orbit integrals
    integer, parameter :: NUM_PQUAD = 40
    ! uses \var(NBAST) (number of orbitals)
#include "inforb.h"
    !-! number of orbitals
    !-integer num_orb
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
    ! number of expectation values
    integer num_expt
    ! addresses for expectation values of square matrices
    integer strt_sq_expt, end_sq_expt
    ! addresses for expectation values from Dalton's routines
    integer strt_dal_expt, end_dal_expt
    ! dimension of integral matrices
    integer dim_int
    ! size of integrals
    integer size_int
    ! addresses for integrals from Gen1Int library
    integer strt_gen_int, end_gen_int
    ! addresses for integrals from Dalton's routines
    integer strt_dal_int, end_dal_int
    ! base address of free Dalton workspace
    integer base_free
    ! length of free Dalton workspace
    integer len_free
    ! ratio of values from Gen1Int to Dalton's routines
    real(REALK) ratio_to_dal
    ! indicator if the test failed
    logical test_failed
    ! incremental recorder of tests
    integer itst
    ! indices of row and column
    integer irow, icol
    ! addresses of results from Gen1Int and Dalton's routines
    integer addr_gen, addr_dal
    ! error information
    integer ierr
    call QENTER("gen1int_ifc_test")
    ! checks if the Gen1Int interface is initialized
    if (.not.shell_init) &
      call QUIT("gen1int_ifc_test>> Gen1Int interface is not properly initialized!")
    !-! gets the number of orbitals
    !-call gen1int_shell_idx_orb(ao_shells(num_ao_shells), ierr, num_orb)
    write(io_std,100) "threshold of error", ERR_THRSH
    write(io_std,110) "threshold of ratio to the referenced result", RATIO_THRSH
    ! sets artifical AO density matrices
    dim_tri_dens = NBAST*(NBAST+1)/2
    dim_sq_dens = NBAST*NBAST
    allocate(sym_ao_dens(dim_tri_dens), stat=ierr)
    if (ierr/=0) call QUIT("gen1int_ifc_test>> failed to allocate sym_ao_dens!")
    allocate(anti_ao_dens(dim_tri_dens), stat=ierr)
    if (ierr/=0) call QUIT("gen1int_ifc_test>> failed to allocate anti_ao_dens!")
    allocate(sq_sym_dens(NBAST,NBAST), stat=ierr)
    if (ierr/=0) call QUIT("gen1int_ifc_test>> failed to allocate sq_sym_dens!")
    allocate(sq_anti_dens(NBAST,NBAST), stat=ierr)
    if (ierr/=0) call QUIT("gen1int_ifc_test>> failed to allocate sq_anti_dens!")
    allocate(sq_gen_dens(NBAST,NBAST), stat=ierr)
    if (ierr/=0) call QUIT("gen1int_ifc_test>> failed to allocate sq_gen_dens!")
    call random_number(sym_ao_dens)
    call random_number(sq_gen_dens)
    ! we uses the upper and diagonal parts in triangular format
    addr_gen = 0
    do icol = 1, NBAST
      do irow = 1, icol-1
        addr_gen = addr_gen+1
        sym_ao_dens(addr_gen) = 0.1_REALK*sym_ao_dens(addr_gen)
        anti_ao_dens(addr_gen) = sym_ao_dens(addr_gen)
        sq_sym_dens(irow,icol) = sym_ao_dens(addr_gen)
        sq_sym_dens(icol,irow) = sym_ao_dens(addr_gen)
        sq_anti_dens(irow,icol) = sym_ao_dens(addr_gen)
        sq_anti_dens(icol,irow) = -sym_ao_dens(addr_gen)
      end do
      addr_gen = addr_gen+1
      sym_ao_dens(addr_gen) = 0.1_REALK*sym_ao_dens(addr_gen)
      anti_ao_dens(addr_gen) = 0.0_REALK
      sq_sym_dens(icol,icol) = sym_ao_dens(addr_gen)
      sq_anti_dens(icol,icol) = 0.0_REALK
    end do
    ! loops over different expectation value tests
    do itst = 1, NUM_EXP_TEST
      ! gets the number of operators including derivatives
      call gen1int_prop_attr(PROP_EXPT(itst), IS_LAO(itst), ORDER_MOM(itst),  &
                             MAX_NUM_CENT(itst), ORDER_GEO_TOT(itst), NUCDEP, &
                             num_opt_derv, kind_int)
      ! number and addresses of expectation values
      num_expt = NUM_DENS*num_opt_derv
      strt_sq_expt = num_expt+1
      end_sq_expt = 2*num_expt
      ! format of integral matrix
      select case(kind_int)
      case(1, -1)
        is_triang = .true.
      case default
        is_triang = .false.
      end select
      ! addresses for integrals
      strt_gen_int = end_sq_expt+1
      end_gen_int = strt_gen_int
      if (end_gen_int>len_work) then
        write(io_std,100) "expectation value test of "//trim(PROP_EXPT(itst))// &
                          ".Gen1Int failed!"
        write(io_std,999) "required memory", end_dal_int
        write(io_std,999) "available memory", len_work
        call QUIT("gen1int_ifc_test>> increase Dalton workspace!")
      end if
      ! computes the expectation values of symmetric AO density matrix in triangular format
      dal_work(1:num_expt) = 0.0_REALK
      call gen1int_ifc_main(PROP_EXPT(itst), IS_LAO(itst), ORDER_MOM(itst), &
                            MAX_NUM_CENT(itst), ORDER_GEO_TOT(itst),        &
                            is_triang, .false., WRT_INT, 1,                 &
                            dal_work(strt_gen_int:end_gen_int),             &
                            1, dim_tri_dens, NUM_DENS, sym_ao_dens, .true., &
                            .false., WRT_EXPT, dal_work(1:num_expt),        &
                            io_std, level_print)
      ! computes the expectation values of symmetric AO density matrix in square format
      dal_work(strt_sq_expt:end_sq_expt) = 0.0_REALK
      call gen1int_ifc_main(PROP_EXPT(itst), IS_LAO(itst), ORDER_MOM(itst),        &
                            MAX_NUM_CENT(itst), ORDER_GEO_TOT(itst),               &
                            is_triang, .false., WRT_INT, 1,                        &
                            dal_work(strt_gen_int:end_gen_int),                    &
                            0, dim_sq_dens, NUM_DENS, sq_sym_dens, .true.,         &
                            .false., WRT_EXPT, dal_work(strt_sq_expt:end_sq_expt), &
                            io_std, level_print)
      if (kind_int==-1) then
        if (any(abs(dal_work(1:num_expt))>ERR_THRSH)) then
          test_failed = .true.
          write(io_std,100) "expectation value of symmetric AO density matrix"
          write(io_std,100) "and anti-symmetric integral matrix should be zero!"
          write(io_std,100) "wrong result using triangular AO density matrix!"
        else
          test_failed = .false.
        end if
        if (any(abs(dal_work(strt_sq_expt:end_sq_expt))>ERR_THRSH)) then
          test_failed = .true.
          write(io_std,100) "expectation value of symmetric AO density matrix"
          write(io_std,100) "and anti-symmetric integral matrix should be zero!"
          write(io_std,100) "wrong result using square AO density matrix!"
        else
          test_failed = .false.
        end if
      else
        test_failed = .false.
      end if
      if (.not.test_failed) then
        do ierr = 1, num_expt
          ! we check the ratio for absolutely greater values
          if (abs(dal_work(num_expt+ierr))>ERR_THRSH) then
            ratio_to_dal = dal_work(ierr)/dal_work(num_expt+ierr)
            if (ratio_to_dal<RATIO_THRSH(1) .or. ratio_to_dal>RATIO_THRSH(2)) then
              write(io_std,995) trim(PROP_EXPT(itst)), ierr, dal_work(num_expt+ierr), &
                                dal_work(ierr)
              test_failed = .true.
            end if
          ! checks the difference for smaller values
          else
            if (abs(dal_work(ierr)-dal_work(num_expt+ierr))>ERR_THRSH) then
              write(io_std,995) trim(PROP_EXPT(itst)), ierr, dal_work(num_expt+ierr), &
                                dal_work(ierr)
              test_failed = .true.
            end if
          end if
        end do
      end if
      if (test_failed) then
        write(io_std,100) "expectation value test of "//trim(PROP_EXPT(itst))// &
                          ".Gen1Int failed!"
      else
        write(io_std,100) "expectation value test of "//trim(PROP_EXPT(itst))// &
                          ".Gen1Int passed!"
      end if
      ! computes the expectation values of anti-symmetric AO density matrix in triangular format
      dal_work(1:num_expt) = 0.0_REALK
      call gen1int_ifc_main(PROP_EXPT(itst), IS_LAO(itst), ORDER_MOM(itst),   &
                            MAX_NUM_CENT(itst), ORDER_GEO_TOT(itst),          &
                            is_triang, .false., WRT_INT, 1,                   &
                            dal_work(strt_gen_int:end_gen_int),               &
                            -1, dim_tri_dens, NUM_DENS, anti_ao_dens, .true., &
                            .false., WRT_EXPT, dal_work(1:num_expt),          &
                            io_std, level_print)
      ! computes the expectation values of anti-symmetric AO density matrix in square format
      dal_work(strt_sq_expt:end_sq_expt) = 0.0_REALK
      call gen1int_ifc_main(PROP_EXPT(itst), IS_LAO(itst), ORDER_MOM(itst),        &
                            MAX_NUM_CENT(itst), ORDER_GEO_TOT(itst),               &
                            is_triang, .false., WRT_INT, 1,                        &
                            dal_work(strt_gen_int:end_gen_int),                    &
                            0, dim_sq_dens, NUM_DENS, sq_anti_dens, .true.,        &
                            .false., WRT_EXPT, dal_work(strt_sq_expt:end_sq_expt), &
                            io_std, level_print)
      if (kind_int==1) then
        if (any(abs(dal_work(1:num_expt))>ERR_THRSH)) then
          test_failed = .true.
          write(io_std,100) "expectation value of anti-symmetric AO density matrix"
          write(io_std,100) "and symmetric integral matrix should be zero!"
          write(io_std,100) "wrong result using triangular AO density matrix!"
        else
          test_failed = .false.
        end if
        if (any(abs(dal_work(strt_sq_expt:end_sq_expt))>ERR_THRSH)) then
          test_failed = .true.
          write(io_std,100) "expectation value of anti-symmetric AO density matrix"
          write(io_std,100) "and symmetric integral matrix should be zero!"
          write(io_std,100) "wrong result using square AO density matrix!"
        else
          test_failed = .false.
        end if
      else
        test_failed = .false.
      end if
      if (.not.test_failed) then
        do ierr = 1, num_expt
          ! we check the ratio for absolutely greater values
          if (abs(dal_work(num_expt+ierr))>ERR_THRSH) then
            ratio_to_dal = dal_work(ierr)/dal_work(num_expt+ierr)
            if (ratio_to_dal<RATIO_THRSH(1) .or. ratio_to_dal>RATIO_THRSH(2)) then
              write(io_std,996) trim(PROP_EXPT(itst)), ierr, dal_work(num_expt+ierr), &
                                dal_work(ierr)
              test_failed = .true.
            end if
          ! checks the difference for smaller values
          else
            if (abs(dal_work(ierr)-dal_work(num_expt+ierr))>ERR_THRSH) then
              write(io_std,996) trim(PROP_EXPT(itst)), ierr, dal_work(num_expt+ierr), &
                                dal_work(ierr)
              test_failed = .true.
            end if
          end if
        end do
      end if
      if (test_failed) then
        write(io_std,100) "expectation value test of "//trim(PROP_EXPT(itst))// &
                          ".Gen1Int failed!"
      else
        write(io_std,100) "expectation value test of "//trim(PROP_EXPT(itst))// &
                          ".Gen1Int passed!"
      end if
    end do
    deallocate(anti_ao_dens)
    deallocate(sq_sym_dens)
    deallocate(sq_anti_dens)
    ! loops over different tests
    do itst = 1, NUM_TEST
      ! gets the number of operators including derivatives
      call gen1int_prop_attr(PROP_NAME(itst), IS_LAO(itst), ORDER_MOM(itst),  &
                             MAX_NUM_CENT(itst), ORDER_GEO_TOT(itst), NUCDEP, &
                             num_opt_derv, kind_int)
      ! number and addresses of expectation values
      num_expt = NUM_DENS*num_opt_derv
      strt_dal_expt = num_expt+1
      end_dal_expt = 2*num_expt
      ! size of integrals
      select case(kind_int)
      case(1, -1)
        !-dim_int = num_orb*(num_orb+1)/2
        dim_int = NBAST*(NBAST+1)/2
        is_triang = .true.
      case default
        !-dim_int = num_orb*num_orb
        dim_int = NBAST*NBAST
        is_triang = .false.
      end select
      size_int = dim_int*num_opt_derv
      ! addresses for integrals and referenced results
      strt_gen_int = end_dal_expt+1
      strt_dal_int = strt_gen_int+size_int
      end_gen_int = strt_dal_int-1
      end_dal_int = strt_dal_int+size_int-1
      if (end_dal_int>len_work) then
        write(io_std,100) "test of "//trim(PROP_NAME(itst))//".Gen1Int failed!"
        write(io_std,999) "required memory", end_dal_int
        write(io_std,999) "available memory", len_work
        call QUIT("gen1int_ifc_test>> increase Dalton workspace!")
      end if
      ! computes the integrals and/or expectation values
      dal_work(1:num_expt) = 0.0_REALK
      if (is_triang) then
        call gen1int_ifc_main(PROP_NAME(itst), IS_LAO(itst), ORDER_MOM(itst),    &
                              MAX_NUM_CENT(itst), ORDER_GEO_TOT(itst),           &
                              is_triang, GET_INT, WRT_INT, dim_int,              &
                              dal_work(strt_gen_int:end_gen_int),                &
                              kind_int, dim_tri_dens, NUM_DENS, sym_ao_dens,     &
                              GET_EXPT, .false., WRT_EXPT, dal_work(1:num_expt), &
                              io_std, level_print)
      else
        call gen1int_ifc_main(PROP_NAME(itst), IS_LAO(itst), ORDER_MOM(itst),    &
                              MAX_NUM_CENT(itst), ORDER_GEO_TOT(itst),           &
                              is_triang, GET_INT, WRT_INT, dim_int,              &
                              dal_work(strt_gen_int:end_gen_int),                &
                              kind_int, dim_sq_dens, NUM_DENS, sq_gen_dens,      &
                              GET_EXPT, .false., WRT_EXPT, dal_work(1:num_expt), &
                              io_std, level_print)
      end if
      ! gets the referenced results from HERMIT
!FIXME: \var(FORQM3)
      if (trim(PROP_NAME(itst))=="DSO") then
        max_typ = (3*NUCDEP)**2
      else
        max_typ = 3*MXCOOR
      end if
      allocate(int_rep(max_typ), stat=ierr)
      if (ierr/=0) call QUIT("gen1int_ifc_test>> failed to allocate int_rep!")
      !
      if (trim(PROP_NAME(itst))=="ELFGRDC" .or. trim(PROP_NAME(itst))=="ELFGRDS") then
         lint_ad = 9*NUCIND*(MAXREP+1)
      else
         lint_ad = max_typ
      end if
      allocate(int_adr(lint_ad), stat=ierr)
      if (ierr/=0) call QUIT("gen1int_ifc_test>> failed to allocate int_adr!")
      !
      allocate(lb_int(max_typ), stat=ierr)
      if (ierr/=0) call QUIT("gen1int_ifc_test>> failed to allocate lb_int!")
      base_free = 1
      len_free = len_work-end_dal_int
      write(io_std,100) "gets the referenced results from HERMIT ..."
      dal_work(strt_dal_expt:end_dal_expt) = 0.0_REALK
      if (is_triang) then
        call PR1IN1(dal_work(end_dal_int+1:), base_free, len_free, int_rep, int_adr,    &
                    lb_int, DAL_PROP(itst)(1:7), ORDER_MOM(itst), NUM_PQUAD, is_triang, &
                    PROP_PRINT, level_print, dal_work(strt_dal_int:end_dal_int),        &
                    NCOMP, TOFILE, MTFORM, DOINT, dal_work(strt_dal_expt:end_dal_expt), &
                    GET_EXPT, sym_ao_dens)
      else
        call PR1IN1(dal_work(end_dal_int+1:), base_free, len_free, int_rep, int_adr,    &
                    lb_int, DAL_PROP(itst)(1:7), ORDER_MOM(itst), NUM_PQUAD, is_triang, &
                    PROP_PRINT, level_print, dal_work(strt_dal_int:end_dal_int),        &
                    NCOMP, TOFILE, MTFORM, DOINT, dal_work(strt_dal_expt:end_dal_expt), &
                    GET_EXPT, sq_gen_dens)
      end if
!FIXME: why the first time calling PR1IN1 gives wrong results??
      if (itst==1) then
        base_free = 1
        len_free = len_work-end_dal_int
        dal_work(strt_dal_expt:end_dal_expt) = 0.0_REALK
        if (is_triang) then
          call PR1IN1(dal_work(end_dal_int+1:), base_free, len_free, int_rep, int_adr,    &
                      lb_int, DAL_PROP(itst)(1:7), ORDER_MOM(itst), NUM_PQUAD, is_triang, &
                      PROP_PRINT, level_print, dal_work(strt_dal_int:end_dal_int),        &
                      NCOMP, TOFILE, MTFORM, DOINT, dal_work(strt_dal_expt:end_dal_expt), &
                      GET_EXPT, sym_ao_dens)
        else
          call PR1IN1(dal_work(end_dal_int+1:), base_free, len_free, int_rep, int_adr,    &
                      lb_int, DAL_PROP(itst)(1:7), ORDER_MOM(itst), NUM_PQUAD, is_triang, &
                      PROP_PRINT, level_print, dal_work(strt_dal_int:end_dal_int),        &
                      NCOMP, TOFILE, MTFORM, DOINT, dal_work(strt_dal_expt:end_dal_expt), &
                      GET_EXPT, sq_gen_dens)
        end if
      end if
      deallocate(int_rep)
      deallocate(int_adr)
      deallocate(lb_int)
      ! Dalton uses different sign for the following integrals
      if (DAL_PROP(itst)=="DPLGRA  " .or. DAL_PROP(itst)=="POTENERG" .or. &
          DAL_PROP(itst)=="NUCSLO  " .or. DAL_PROP(itst)=="PSO     ") then
        dal_work(strt_gen_int:end_gen_int) = -dal_work(strt_gen_int:end_gen_int)
        dal_work(1:num_expt) = -dal_work(1:num_expt)
      end if
      ! checks the results
      write(io_std,100) "checks the results of "//trim(PROP_NAME(itst))
      test_failed = .false.
      addr_gen = end_dal_expt
      addr_dal = end_gen_int
      if (is_triang) then
        do icol = 1, NBAST
          do irow = 1, icol
            addr_gen = addr_gen+1
            addr_dal = addr_dal+1
            ! we check the ratio for absolutely greater values
            if (abs(dal_work(addr_dal))>ERR_THRSH) then
              ratio_to_dal = dal_work(addr_gen)/dal_work(addr_dal)
              if (ratio_to_dal<RATIO_THRSH(1) .or. ratio_to_dal>RATIO_THRSH(2)) then
                write(io_std,998) trim(PROP_NAME(itst)), irow, icol, &
                                  dal_work(addr_dal), dal_work(addr_gen)
                test_failed = .true.
              end if
            ! checks the difference for smaller values
            else
              if (abs(dal_work(addr_gen)-dal_work(addr_dal))>ERR_THRSH) then
                write(io_std,998) trim(PROP_NAME(itst)), irow, icol, &
                                  dal_work(addr_dal), dal_work(addr_gen)
                test_failed = .true.
              end if
            end if
          end do
        end do
      else
        do icol = 1, NBAST
          do irow = 1, NBAST
            addr_gen = addr_gen+1
            addr_dal = addr_dal+1
            ! we check the ratio for absolutely greater values
            if (abs(dal_work(addr_dal))>ERR_THRSH) then
              ratio_to_dal = dal_work(addr_gen)/dal_work(addr_dal)
              if (ratio_to_dal<RATIO_THRSH(1) .or. ratio_to_dal>RATIO_THRSH(2)) then
                write(io_std,998) trim(PROP_NAME(itst)), irow, icol, &
                                  dal_work(addr_dal), dal_work(addr_gen)
                test_failed = .true.
              end if
            ! checks the difference for smaller values
            else
              if (abs(dal_work(addr_gen)-dal_work(addr_dal))>ERR_THRSH) then
                write(io_std,998) trim(PROP_NAME(itst)), irow, icol, &
                                  dal_work(addr_dal), dal_work(addr_gen)
                test_failed = .true.
              end if
            end if
          end do
        end do
      end if
      if (GET_EXPT) then
        do ierr = 1, num_expt
          ! we check the ratio for absolutely greater values
          if (abs(dal_work(num_expt+ierr))>ERR_THRSH) then
            ratio_to_dal = dal_work(ierr)/dal_work(num_expt+ierr)
            if (ratio_to_dal<RATIO_THRSH(1) .or. ratio_to_dal>RATIO_THRSH(2)) then
              write(io_std,997) trim(PROP_NAME(itst)), ierr, dal_work(num_expt+ierr), &
                                dal_work(ierr)
              test_failed = .true.
            end if
          ! checks the difference for smaller values
          else
            if (abs(dal_work(ierr)-dal_work(num_expt+ierr))>ERR_THRSH) then
              write(io_std,997) trim(PROP_NAME(itst)), ierr, dal_work(num_expt+ierr), &
                                dal_work(ierr)
              test_failed = .true.
            end if
          end if
        end do
      end if
      if (test_failed) then
        write(io_std,100) "test of "//trim(PROP_NAME(itst))//".Gen1Int failed!"
      else
        write(io_std,100) "test of "//trim(PROP_NAME(itst))//".Gen1Int passed!"
      end if
    end do
    deallocate(sym_ao_dens)
    deallocate(sq_gen_dens)
    ! cleans up the data in Gen1Int interface after all calculations will be
    ! performed in abacus/dalton.F
    call QEXIT("gen1int_ifc_test")
    return
100 format("gen1int_ifc_test>> ",A,Es16.6)
110 format("gen1int_ifc_test>> ",A,Es16.6,"-->",Es16.6)
995 format("gen1int_ifc_test>> ",A,".EXP(",I6,")>> SQUARE>>",F18.12,", SYM>>",F18.12)
996 format("gen1int_ifc_test>> ",A,".EXP(",I6,")>> SQUARE>>",F18.12,", ANTI>>",F18.12)
997 format("gen1int_ifc_test>> ",A,".EXP(",I6,")>> HERMIT>>",F18.12,", Gen1Int>>",F18.12)
998 format("gen1int_ifc_test>> ",A,".INT(",I6,",",I6,")>> HERMIT>>",F18.12, &
           ", Gen1Int>>",F18.12)
999 format("gen1int_ifc_test>> ",A,I16)
#else
    write(io_std,100) "Gen1Int is not installed! No tests will be performed!"
    return
100 format("gen1int_ifc_test>> ",A)
#endif
  end subroutine gen1int_ifc_test
