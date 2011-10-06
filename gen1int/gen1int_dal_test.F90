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
  !> \param io_unit is the IO unit of standard output
  !> \param level_print is the level of print
  subroutine gen1int_dal_test(len_work, dal_work, io_unit, level_print)
    ! AO sub-shells
    use gen1int_shell
    implicit none
    integer, intent(in) :: len_work
    real(REALK), intent(inout) :: dal_work(len_work)
    integer, intent(in) :: io_unit
    integer, intent(in) :: level_print
#ifdef BUILD_GEN1INT
    ! number of tests
    integer, parameter :: NUM_TEST = 1
    ! property integral names
    character*7, parameter :: PROP_NAME(NUM_TEST) = &
      (/"OVERLAP"/)
    ! if London atomic orbitals
    logical, parameter :: IS_LAO(NUM_TEST) = &
      (/.false./)
    ! order of magnetic derivatives on bra center
    integer, parameter :: ORDER_MAG_BRA(NUM_TEST) = &
      (/0/)
    ! order of magnetic derivatives on ket center
    integer, parameter :: ORDER_MAG_KET(NUM_TEST) = &
      (/0/)
    ! order of total magnetic derivatives
    integer, parameter :: ORDER_MAG_TOT(NUM_TEST) = &
      (/0/)
    ! order of derivatives w.r.t. total rotational angular momentum on bra center
    integer, parameter :: ORDER_RAM_BRA(NUM_TEST) = &
      (/0/)
    ! order of derivatives w.r.t. total rotational angular momentum on ket center
    integer, parameter :: ORDER_RAM_KET(NUM_TEST) = &
      (/0/)
    ! order of total derivatives w.r.t. total rotational angular momentum
    integer, parameter :: ORDER_RAM_TOT(NUM_TEST) = &
      (/0/)
    ! order of geometric derivatives on bra center
    integer, parameter :: ORDER_GEO_BRA(NUM_TEST) = &
      (/0/)
    ! order of geometric derivatives on ket center
    integer, parameter :: ORDER_GEO_KET(NUM_TEST) = &
      (/0/)
    ! maximum number of differentiated centers
    integer, parameter :: MAX_NUM_CENT(NUM_TEST) = &
      (/0/)
    ! order of total geometric derivatives
    integer, parameter :: ORDER_GEO_TOT(NUM_TEST) = &
      (/0/)
    ! if the integral matrices are symmetric
    logical, parameter :: SYM_INT(NUM_TEST) = &
      (/.true./)
    ! if getting integrals back from Gen1Int
    logical, parameter :: GET_INT(NUM_TEST) = &
      (/.true./)
    ! if writing integrals on file
    logical, parameter :: WRT_INT(NUM_TEST) = &
      (/.false./)
    ! if calculating expectation values
    logical, parameter :: DO_EXPT(NUM_TEST) = &
      (/.false./)
    ! number of operators including derivatives
    integer, parameter :: NUM_OPT(NUM_TEST) = &
      (/1/)
    ! number of AO density matrices
    integer, parameter :: NUM_DENS = 1
    ! AO density matrices
    real(REALK), allocatable :: ao_dens(:,:)
    ! if getting expectation values back
    logical, parameter :: GET_EXPT(NUM_TEST) = &
      (/.false./)
    ! if writing expectation values on file
    logical, parameter :: WRT_EXPT(NUM_TEST) = &
      (/.false./)
    ! number of orbitals
    integer num_orb
    ! number of expectation values
    integer num_expt
    ! dimension of integral matrices
    integer dim_int
    ! size of integrals
    integer size_int
    ! addresses for integrals
    integer start_int, end_int
    ! addresses for referenced results
    integer start_ref, end_ref
    ! incremental recorder of tests
    integer itst
    ! error information
    integer ierr
    call QENTER("gen1int_dal_test")
    ! initializes the data used in Gen1Int interface before any calculation
    call gen1int_dal_init
    ! gets the number of orbitals
    call gen1int_shell_num_orb(num_ao_shells, ao_shells, num_orb)
    ! loops over different tests
    do itst = 1, NUM_TEST
      ! number of expectation values
      num_expt = NUM_DENS*NUM_OPT(itst)
      ! size of integrals
      if (SYM_INT(itst)) then
        dim_int = num_orb*(num_orb+1)/2
      else
        dim_int = num_orb*num_orb
      end if
!FIXME: ao_dens
      allocate(ao_dens(dim_int,NUM_DENS), stat=ierr)
      if (ierr/=0) call QUIT("gen1int_dal_test>> failed to allocate ao_dens!")
      size_int = dim_int*NUM_OPT(itst)
      ! addresses for integrals and referenced results
      start_int = num_expt+1
      start_ref = start_int+size_int
      end_int = start_ref-1
      end_ref = start_ref+size_int-1
      if (end_ref>len_work) then
        write(io_unit,100) "ID of test", itst
        write(io_unit,100) "required memory", end_ref
        write(io_unit,100) "available memory", len_work
        call QUIT("Increase Dalton workspace!")
      end if
      ! computes the integrals and/or expectation values
      call gen1int_dal_main(PROP_NAME(itst), IS_LAO(itst),                     &
             ORDER_MAG_BRA(itst), ORDER_MAG_KET(itst), ORDER_MAG_TOT(itst),    &
             ORDER_RAM_BRA(itst), ORDER_RAM_KET(itst), ORDER_RAM_TOT(itst),    &
             ORDER_GEO_BRA(itst), ORDER_GEO_KET(itst), MAX_NUM_CENT(itst),     &
             ORDER_GEO_TOT(itst), SYM_INT(itst), GET_INT(itst), WRT_INT(itst), &
             dim_int, dal_work(start_int:end_int), DO_EXPT(itst), NUM_DENS,    &
             ao_dens, GET_EXPT(itst), WRT_EXPT(itst), dal_work(1:num_expt),    &
             io_unit, level_print)
      ! gets the referenced results
!FIXME
      deallocate(ao_dens)
      ! checks the results
!FIXME
      call AROUND('gen1int_dal_test>> '//trim(PROP_NAME(itst)))
      if (SYM_INT(itst)) then
        call OUTPAK(dal_work(start_int:end_int), num_orb, 1, io_unit)
      else
        call OUTPUT(dal_work(start_int:end_int), 1, num_orb, 1, num_orb, &
                    num_orb, num_orb, 1, io_unit)
      end if
    end do
    ! cleans up the data in Gen1Int interface after all calculations will be
    ! performed in abacus/dalton.F
    call QEXIT("gen1int_dal_test")
    return
100 format("gen1int_dal_test>> ",A,I16)
#else
    call QUIT("Gen1Int is not installed!")
#endif
  end subroutine gen1int_dal_test
