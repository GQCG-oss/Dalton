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
  subroutine gen1int_dal_test(len_work, dal_work, io_std, level_print)
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
    ! number of tests
    integer, parameter :: NUM_TEST = 3
    ! property integral names, see subroutine \fn(gen1int_shell_prop) in gen1int_shell.F90,
    ! and variable \var(TABLE) in subroutine \fn(PR1IN1) in abacus/her1pro.F
    character*8, parameter :: PROP_NAME(NUM_TEST) = &
      (/"KINENERG", "OVERLAP ", "POTENERG"/)
    ! if London atomic orbitals
    logical, parameter :: IS_LAO(NUM_TEST) = &
      (/.false., .false., .false./)
    ! order of Cartesian multipole moments
    integer, parameter :: ORDER_MOM(NUM_TEST) = &
      (/0, 0, 0/)
    ! maximum number of differentiated centers
    integer, parameter :: MAX_NUM_CENT(NUM_TEST) = &
      (/0, 0, 0/)
    ! order of total geometric derivatives
    integer, parameter :: ORDER_GEO_TOT(NUM_TEST) = &
      (/0, 0, 0/)
    ! if the integral matrices are symmetric, anti-symmetric or square
    integer, parameter :: KIND_INT(NUM_TEST) = &
      (/1, 1, 1/)
    ! if getting integrals back from Gen1Int
    logical, parameter :: GET_INT(NUM_TEST) = &
      (/.true., .true., .true./)
    ! if writing integrals on file
    logical, parameter :: WRT_INT(NUM_TEST) = &
      (/.false., .false., .false./)
    ! if calculating expectation values
    logical, parameter :: DO_EXPT(NUM_TEST) = &
      (/.false., .false., .false./)
    ! if getting expectation values back
    logical, parameter :: GET_EXPT(NUM_TEST) = &
      (/.false., .false., .false./)
    ! if writing expectation values on file
    logical, parameter :: WRT_EXPT(NUM_TEST) = &
      (/.false., .false., .false./)
    ! number of operators including derivatives
    integer num_opt_derv
    ! number of AO density matrices
    integer, parameter :: NUM_DENS = 1
    ! AO density matrices
    real(REALK), allocatable :: ao_dens(:,:)
    !
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
    ! integrals are triangularized or squared
    logical is_triang
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
    ! error information
    integer ierr
    call QENTER("gen1int_dal_test")
    ! checks if the Gen1Int interface is initialized
    if (.not.shell_init) &
      call QUIT("gen1int_dal_test>> Gen1Int interface is not properly initialized!")
    !-! gets the number of orbitals
    !-call gen1int_shell_idx_orb(ao_shells(num_ao_shells), ierr, num_orb)
    write(io_std,100) "threshold of error", ERR_THRSH
    write(io_std,110) "threshold of ratio to the referenced result", RATIO_THRSH
    ! loops over different tests
    do itst = 1, NUM_TEST
      ! gets the number of operators including derivatives
      call gen1int_num_opt(PROP_NAME(itst), IS_LAO(itst), ORDER_MOM(itst),  &
                           MAX_NUM_CENT(itst), ORDER_GEO_TOT(itst), NUCDEP, &
                           num_opt_derv)
      ! number and addresses of expectation values
      num_expt = NUM_DENS*num_opt_derv
      strt_dal_expt = num_expt+1
      end_dal_expt = 2*num_expt
      ! size of integrals
      select case(KIND_INT(itst))
      case(1, -1)
        !-dim_int = num_orb*(num_orb+1)/2
        dim_int = NBAST*(NBAST+1)/2
        is_triang = .true.
      case default
        !-dim_int = num_orb*num_orb
        dim_int = NBAST*NBAST
        is_triang = .false.
      end select
!FIXME: ao_dens
      allocate(ao_dens(dim_int,NUM_DENS), stat=ierr)
      if (ierr/=0) call QUIT("gen1int_dal_test>> failed to allocate ao_dens!")
      size_int = dim_int*num_opt_derv
      ! addresses for integrals and referenced results
      strt_gen_int = end_dal_expt+1
      strt_dal_int = strt_gen_int+size_int
      end_gen_int = strt_dal_int-1
      end_dal_int = strt_dal_int+size_int-1
      if (end_dal_int>len_work) then
        write(io_std,999) "ID of test", itst
        write(io_std,999) "required memory", end_dal_int
        write(io_std,999) "available memory", len_work
        call QUIT("gen1int_dal_test>> increase Dalton workspace!")
      end if
      ! computes the integrals and/or expectation values
      call gen1int_dal_main(PROP_NAME(itst), IS_LAO(itst), ORDER_MOM(itst),       &
                            MAX_NUM_CENT(itst), ORDER_GEO_TOT(itst),              &
                            KIND_INT(itst), GET_INT(itst), WRT_INT(itst),         &
                            dim_int, dal_work(strt_gen_int:end_gen_int),          &
                            DO_EXPT(itst), NUM_DENS, ao_dens,                     &
                            GET_EXPT(itst), WRT_EXPT(itst), dal_work(1:num_expt), &
                            io_std, level_print)
      ! gets the referenced results from HERMIT
!FIXME: \var(FORQM3)
      if (trim(PROP_NAME(itst))=="DSO") then
        max_typ = (3*NUCDEP)**2
      else
        max_typ = 3*MXCOOR
      end if
      allocate(int_rep(max_typ), stat=ierr)
      if (ierr/=0) call QUIT("gen1int_dal_test>> failed to allocate int_rep!")
      !
      if (trim(PROP_NAME(itst))=="ELFGRDC" .or. trim(PROP_NAME(itst))=="ELFGRDS") then
         lint_ad = 9*NUCIND*(MAXREP+1)
      else
         lint_ad = max_typ
      end if
      allocate(int_adr(lint_ad), stat=ierr)
      if (ierr/=0) call QUIT("gen1int_dal_test>> failed to allocate int_adr!")
      !
      allocate(lb_int(max_typ), stat=ierr)
      if (ierr/=0) call QUIT("gen1int_dal_test>> failed to allocate lb_int!")
      base_free = 1
      len_free = len_work-end_dal_int
      write(io_std,100) "gets the referenced results from HERMIT ..."
      call PR1IN1(dal_work(end_dal_int+1:), base_free, len_free, int_rep, int_adr,     &
                  lb_int, PROP_NAME(itst)(1:7), ORDER_MOM(itst), NUM_PQUAD, is_triang, &
                  PROP_PRINT, level_print, dal_work(strt_dal_int:end_dal_int),         &
                  NCOMP, TOFILE, MTFORM, DOINT, dal_work(strt_dal_expt:end_dal_expt),  &
                  DO_EXPT(itst), ao_dens)
!FIXME: why the first time calling PR1IN1 gives wrong results??
      if (itst==1) then
        base_free = 1
        len_free = len_work-end_dal_int
        call PR1IN1(dal_work(end_dal_int+1:), base_free, len_free, int_rep, int_adr,     &
                    lb_int, PROP_NAME(itst)(1:7), ORDER_MOM(itst), NUM_PQUAD, is_triang, &
                    PROP_PRINT, level_print, dal_work(strt_dal_int:end_dal_int),         &
                    NCOMP, TOFILE, MTFORM, DOINT, dal_work(strt_dal_expt:end_dal_expt),  &
                    DO_EXPT(itst), ao_dens)
      end if
      deallocate(int_rep)
      deallocate(int_adr)
      deallocate(lb_int)
      deallocate(ao_dens)
      ! checks the results
      write(io_std,100) "checks the results of "//trim(PROP_NAME(itst))
      test_failed = .false.
      do ierr = 1, size_int
        ! we check the ratio for absolutely greater values
        if (abs(dal_work(end_gen_int+ierr))>ERR_THRSH) then
          ratio_to_dal = dal_work(end_dal_expt+ierr)/dal_work(end_gen_int+ierr)
          if (ratio_to_dal<RATIO_THRSH(1) .or. ratio_to_dal>RATIO_THRSH(2)) then
            write(io_std,998) PROP_NAME(itst), dal_work(end_gen_int+ierr), &
                              dal_work(end_dal_expt+ierr)
            test_failed = .true.
          end if
        ! checks the difference for smaller values
        else
          if (abs(dal_work(end_dal_expt+ierr)-dal_work(end_gen_int+ierr))>ERR_THRSH) then
            write(io_std,998) PROP_NAME(itst), dal_work(end_gen_int+ierr), &
                              dal_work(end_dal_expt+ierr)
            test_failed = .true.
          end if
        end if
      end do
      if (DO_EXPT(itst)) then
        do ierr = 1, num_expt
          ! we check the ratio for absolutely greater values
          if (abs(dal_work(num_expt+ierr))>ERR_THRSH) then
            ratio_to_dal = dal_work(ierr)/dal_work(num_expt+ierr)
            if (ratio_to_dal<RATIO_THRSH(1) .or. ratio_to_dal>RATIO_THRSH(2)) then
              write(io_std,997) PROP_NAME(itst), dal_work(num_expt+ierr), &
                                dal_work(ierr)
              test_failed = .true.
            end if
          ! checks the difference for smaller values
          else
            if (abs(dal_work(ierr)-dal_work(num_expt+ierr))>ERR_THRSH) then
              write(io_std,997) PROP_NAME(itst), dal_work(num_expt+ierr), &
                                dal_work(ierr)
              test_failed = .true.
            end if
          end if
        end do
      end if
      if (test_failed) then
        write(io_std,100) "test of "//trim(PROP_NAME(itst))//" failed!"
      else
        write(io_std,100) "test of "//trim(PROP_NAME(itst))//" passed!"
      end if
    end do
    ! cleans up the data in Gen1Int interface after all calculations will be
    ! performed in abacus/dalton.F
    call QEXIT("gen1int_dal_test")
    return
100 format("gen1int_dal_test>> ",A,Es16.6)
110 format("gen1int_dal_test>> ",A,Es16.6,"-->",Es16.6)
997 format("gen1int_dal_test>> ",A,"_expt, HERMIT>>",F18.12,", Gen1Int>>",F18.12)
998 format("gen1int_dal_test>> ",A,"_int, HERMIT>>",F18.12,", Gen1Int>>",F18.12)
999 format("gen1int_dal_test>> ",A,I16)
#else
    call QUIT("gen1int_dal_test>> Gen1Int is not installed!")
#endif
  end subroutine gen1int_dal_test
