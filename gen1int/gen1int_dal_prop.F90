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
!...  This file calculates property integrals and/or expectation values.
!
!...  2011-10-04, Bin Gao
!...  * first version

#include "xkind.h"

  !> \brief calculates property integrals and/or expectation values
  !> \author Bin Gao
  !> \date 2011-10-04
  !> \param prop_name is the name of property integrals to calculated
  !> \param is_lao indicates if using rotational London atomic orbitals
  !> \param order_mom is the order of Cartesian multipole moments, only used for
  !>        integrals "CARMOM"
  !> \param base_geo_derv is the base address of current total geometric derivatives in the
  !>        list of all total geometric derivatives
  !> \param num_cents is the number of geometric differentiated centers
  !> \param idx_cent contains the indices of differentiated centers
  !> \param order_cent contains the orders of differentiated centers
  !> \param kind_int indicates if the kind of integral matrices, 1 for symmetric, -1 for
  !>        anti-symmetric, others for square
  !> \param get_int indicates if getting integrals back
  !> \param wrt_int indicates if writing integrals to file
  !> \param dim_int is the dimension of integral matrices
  !> \param do_expt indicates if calculating expectation values
  !> \param num_dens is the number of AO density matrices
  !> \param ao_dens contains the AO density matrices
  !> \param get_expt indicates if getting expectation values back
  !> \param wrt_expt indicates if writing expectation values to file
  !> \param io_std is the IO unit of standard output
  !> \param level_print is the level of print
  !> \return val_ints contains the calculated integrals
  !> \return val_expt contains the calculated expectation values
  subroutine gen1int_dal_prop(prop_name, is_lao, order_mom, &
                              base_geo_derv, num_cents,     &
                              idx_cent, order_cent,         &
                              kind_int, get_int, wrt_int,   &
                              dim_int, val_ints,            &
                              do_expt, num_dens, ao_dens,   &
                              get_expt, wrt_expt, val_expt, &
                              io_std, level_print)
    ! AO sub-shells
    use gen1int_shell
    implicit none
    character*(*), intent(in) :: prop_name
    logical, intent(in) :: is_lao
    integer, intent(in) :: order_mom
    integer, intent(in) :: base_geo_derv
    integer, intent(in) :: num_cents
    integer, intent(in) :: idx_cent(num_cents)
    integer, intent(in) :: order_cent(num_cents)
    integer, intent(in) :: kind_int
    logical, intent(in) :: get_int
    logical, intent(in) :: wrt_int
    integer, intent(in) :: dim_int
    real(REALK), intent(inout) :: val_ints(dim_int,*)
    logical, intent(in) :: do_expt
    integer, intent(in) :: num_dens
    real(REALK), intent(in) :: ao_dens(dim_int,num_dens)
    logical, intent(in) :: get_expt
    logical, intent(in) :: wrt_expt
    real(REALK), intent(inout) :: val_expt(num_dens,*)
    integer, intent(in) :: io_std
    integer, intent(in) :: level_print
#ifdef BUILD_GEN1INT
    integer num_geo_derv        !number of total geometric derivatives
    integer icent               !incremental recorder over differetiated centers
    type(shell_int_t) prop_int  !contracted property integrals between two AO sub-shells
    integer ishell, jshell      !incremental recorders over sub-shells
    integer num_ao_orb          !total number of atomic orbitals
    call QENTER("gen1int_dal_prop")
    ! computes the number of total geometric derivatives
    if (num_cents>0) then
      num_geo_derv = (order_cent(1)+1)*(order_cent(1)+2)/2
      do icent = 2, num_cents
        num_geo_derv = num_geo_derv*(order_cent(icent)+1)*(order_cent(icent)+2)/2
      end do
    else
      num_geo_derv = 1
    end if
    if (level_print>=10) &
      write(io_std,100) "number of total geometric derivatives", num_geo_derv
    select case(kind_int)
    ! symmetric integral matrices (column- or ket-major, upper and diagonal parts)
    case(1)
      ! loops over AO sub-shells on ket center
      do jshell = 1, num_ao_shells
        ! loops over AO sub-shells on bra center (upper and diagonal parts)
        do ishell = 1, jshell
          ! different property integrals
          call gen1int_shell_prop(prop_name, ao_shells(ishell), ao_shells(jshell), &
                                  is_lao, order_mom, num_cents, idx_cent,          &
                                  order_cent, num_geo_derv, prop_int)
          ! assigns the returned integrals, and write the integrals to file if required
          if (get_int) then
            if (ishell<jshell) then
              call gen1int_shell_int_tri_off(prop_int, base_geo_derv, &
                                             wrt_int, dim_int, val_ints)
            else
              call gen1int_shell_int_tri_diag(prop_int, base_geo_derv, &
                                              wrt_int, dim_int, val_ints)
            end if
          end if
          ! calculates the expectation values
          if (do_expt) then
            if (ishell<jshell) then
            else
            end if
            if (get_expt) then
            end if
            if (wrt_expt) then
            end if
          end if
          ! cleans contracted property integrals
          call gen1int_shell_int_clean(prop_int)
        end do
      end do
    ! anti-symmetric integral matrices (column- or ket-major, upper and diagonal parts)
    case(-1)
      ! loops over AO sub-shells on ket center
      do jshell = 1, num_ao_shells
        ! loops over AO sub-shells on bra center (upper and diagonal parts)
        do ishell = 1, jshell
          ! different property integrals
          call gen1int_shell_prop(prop_name, ao_shells(ishell), ao_shells(jshell), &
                                  is_lao, order_mom, num_cents, idx_cent,          &
                                  order_cent, num_geo_derv, prop_int)
          ! assigns the returned integrals, and write the integrals to file if required
          if (get_int) then
            if (ishell<jshell) then
              call gen1int_shell_int_tri_off(prop_int, base_geo_derv, &
                                             wrt_int, dim_int, val_ints)
            else
              call gen1int_shell_int_tri_diag(prop_int, base_geo_derv, &
                                              wrt_int, dim_int, val_ints)
            end if
          end if
          ! calculates the expectation values, and notice that the integral matrices
          ! are anti-symmetric
          if (do_expt) then
            if (ishell<jshell) then
            else
            end if
            if (get_expt) then
            end if
            if (wrt_expt) then
            end if
          end if
          ! cleans contracted property integrals
          call gen1int_shell_int_clean(prop_int)
        end do
      end do
    ! square integral matrices (column- or ket-major)
    case default
      ! gets the total number of atomic orbitals
      call gen1int_shell_idx_orb(ao_shells(num_ao_shells), icent, num_ao_orb)
      if (level_print>=10) write(io_std,100) "number of atomic orbitals", num_ao_orb
      ! loops over AO sub-shells on ket center
      do jshell = 1, num_ao_shells
        ! loops over AO sub-shells on bra center
        do ishell = 1, num_ao_shells
          ! different property integrals
          call gen1int_shell_prop(prop_name, ao_shells(ishell), ao_shells(jshell), &
                                  is_lao, order_mom, num_cents, idx_cent,          &
                                  order_cent, num_geo_derv, prop_int)
          ! assigns the returned integrals, and write the integrals to file if required
          if (get_int) then
            call gen1int_shell_int_square(prop_int, base_geo_derv, num_ao_orb, &
                                          wrt_int, val_ints)
          end if
          ! calculates the expectation values
          if (do_expt) then
            if (get_expt) then
            end if
            if (wrt_expt) then
            end if
          end if
          ! cleans contracted property integrals
          call gen1int_shell_int_clean(prop_int)
        end do
      end do
    end select
    call QEXIT("gen1int_dal_prop")
    return
100 format("gen1int_dal_prop>> ",A,I10)
#else
    call QUIT("gen1int_dal_prop>> Gen1Int is not installed!")
#endif
  end subroutine gen1int_dal_prop
