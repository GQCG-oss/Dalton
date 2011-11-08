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
  !> \param is_triang indicates if returning integral matrices are in triangular format,
  !>        otherwise square format
  !> \param get_int indicates if getting integrals back
  !> \param wrt_int indicates if writing integrals on file
  !> \param dim_int is the dimension of integral matrices
  !> \param kind_dens indicates if the kind of AO density matrices, 1 for symmetric, -1 for
  !>        anti-symmetric, others for square
  !> \param dim_dens is the dimension of AO density matrices
  !> \param num_dens is the number of AO density matrices
  !> \param ao_dens contains the AO density matrices
  !> \param get_expt indicates if getting expectation values back
  !> \param redunt_expt indicates if getting expectation values of redundant total
  !>        geometric derivatives
  !> \param wrt_expt indicates if writing expectation values on file
  !> \param io_std is the IO unit of standard output
  !> \param level_print is the level of print
  !> \return val_ints contains the calculated integrals
  !> \return val_expt contains the calculated expectation values
  subroutine gen1int_ifc_prop(prop_name, is_lao, order_mom,   &
                              base_geo_derv, num_cents,       &
                              idx_cent, order_cent,           &
                              is_triang, get_int, wrt_int,    &
                              dim_int, val_ints,              &
                              kind_dens, dim_dens, num_dens,  &
                              ao_dens, get_expt, redunt_expt, &
                              wrt_expt, val_expt, io_std, level_print)
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
    logical, intent(in) :: is_triang
    logical, intent(in) :: get_int
    logical, intent(in) :: wrt_int
    integer, intent(in) :: dim_int
    real(REALK), intent(inout) :: val_ints(dim_int,*)
    integer, intent(in) :: kind_dens
    integer, intent(in) :: dim_dens
    integer, intent(in) :: num_dens
    real(REALK), intent(in) :: ao_dens(dim_dens,num_dens)
    logical, intent(in) :: get_expt
    logical, intent(in) :: redunt_expt
    logical, intent(in) :: wrt_expt
    real(REALK), intent(inout) :: val_expt(num_dens,*)
    integer, intent(in) :: io_std
    integer, intent(in) :: level_print
#ifdef BUILD_GEN1INT
    ! uses \var(MXCENT)
#include "mxcent.h"
    ! uses \var(NUCDEP)
#include "nuclei.h"
    integer order_geo           !order of total geometric derivatives
    integer num_opt             !number of operators
    integer kind_int            !kind of integral matrices, 1 for symmetric, -1 for
                                !anti-symmetric, others for square
    integer num_geo_derv        !number of total geometric derivatives
    logical alloc_expt          !if allocating temporary expectation values
    real(REALK), allocatable :: tmp_expt(:,:,:)
                                !temporary expectation values
    integer num_ao_orb          !total number of atomic orbitals
    integer icent               !incremental recorder over differetiated centers
    type(shell_int_t) prop_int  !contracted property integrals between two AO sub-shells
    integer ishell, jshell      !incremental recorders over sub-shells
    integer ierr                !error information
    call QENTER("gen1int_ifc_prop")
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
    ! gets the kind of integral matrices
    call gen1int_prop_attr(prop_name, is_lao, order_mom, 0, 0, NUCDEP, &
                           num_opt, kind_int)
    ! gets the order of total geometric derivatives
    order_geo = sum(order_cent)
    ! only writes the expectation values but without returning, or gets the
    ! expectation values of redundant total geometric derivatives (order>1)
    alloc_expt = (wrt_expt .and. .not.get_expt) .or. (redunt_expt .and. order_geo>1)
    if (alloc_expt) then
      allocate(tmp_expt(num_dens,num_opt,num_geo_derv), stat=ierr)
      if (ierr/=0) call QUIT("gen1int_ifc_prop>> failed to allocate tmp_expt!")
      tmp_expt = 0.0_REALK
    end if
    ! gets the total number of atomic orbitals
    call gen1int_shell_idx_orb(ao_shells(num_ao_shells), icent, num_ao_orb)
    if (level_print>=10) write(io_std,100) "number of atomic orbitals", num_ao_orb
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
          ! assigns the returned integrals, and write the integrals on file if required
          if (get_int) then
            if (ishell<jshell) then
              call gen1int_shell_int_tri_off(prop_int, base_geo_derv, is_triang,   &
                                             .true., num_ao_orb, wrt_int, dim_int, &
                                             val_ints)
            else
              call gen1int_shell_int_tri_diag(prop_int, base_geo_derv, is_triang, &
                                              num_ao_orb, wrt_int, dim_int, val_ints)
            end if
          end if
          ! calculates the expectation values
          if (alloc_expt) then
            if (ishell<jshell) then
              call gen1int_shell_tr_tri_off(prop_int, 0, .true.,             &
                                            num_ao_orb, kind_dens, dim_dens, &
                                            num_dens, ao_dens, tmp_expt)
            else
              call gen1int_shell_tr_tri_diag(prop_int, 0, .true.,             &
                                             num_ao_orb, kind_dens, dim_dens, &
                                             num_dens, ao_dens, tmp_expt)
            end if
          else if (get_expt) then
            if (ishell<jshell) then
              call gen1int_shell_tr_tri_off(prop_int, base_geo_derv, .true., &
                                            num_ao_orb, kind_dens, dim_dens, &
                                            num_dens, ao_dens, val_expt)
            else
              call gen1int_shell_tr_tri_diag(prop_int, base_geo_derv, .true., &
                                             num_ao_orb, kind_dens, dim_dens, &
                                             num_dens, ao_dens, val_expt)
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
          ! assigns the returned integrals, and write the integrals on file if required
          if (get_int) then
            if (ishell<jshell) then
              call gen1int_shell_int_tri_off(prop_int, base_geo_derv, is_triang,    &
                                             .false., num_ao_orb, wrt_int, dim_int, &
                                             val_ints)
            else
              call gen1int_shell_int_tri_diag(prop_int, base_geo_derv, is_triang, &
                                              num_ao_orb, wrt_int, dim_int, val_ints)
            end if
          end if
          ! calculates the expectation values
          if (alloc_expt) then
            if (ishell<jshell) then
              ! notice that the integrals in \var(prop_int) will be changed afterwards
              call gen1int_shell_tr_tri_off(prop_int, 0, .false.,            &
                                            num_ao_orb, kind_dens, dim_dens, &
                                            num_dens, ao_dens, tmp_expt)
            else
              call gen1int_shell_tr_tri_diag(prop_int, 0, .false.,            &
                                             num_ao_orb, kind_dens, dim_dens, &
                                             num_dens, ao_dens, tmp_expt)
            end if
          else if (get_expt) then
            if (ishell<jshell) then
              ! notice that the integrals in \var(prop_int) will be changed afterwards
              call gen1int_shell_tr_tri_off(prop_int, base_geo_derv, .false., &
                                            num_ao_orb, kind_dens, dim_dens,  &
                                            num_dens, ao_dens, val_expt)
            else
              call gen1int_shell_tr_tri_diag(prop_int, base_geo_derv, .false., &
                                             num_ao_orb, kind_dens, dim_dens,  &
                                             num_dens, ao_dens, val_expt)
            end if
          end if
          ! cleans contracted property integrals
          call gen1int_shell_int_clean(prop_int)
        end do
      end do
    ! square integral matrices (column- or ket-major)
    case default
      if (is_triang) then
        call QUIT(trim(prop_name)//" integrals can not return in triangular format!")
      end if
      ! loops over AO sub-shells on ket center
      do jshell = 1, num_ao_shells
        ! loops over AO sub-shells on bra center
        do ishell = 1, num_ao_shells
          ! different property integrals
          call gen1int_shell_prop(prop_name, ao_shells(ishell), ao_shells(jshell), &
                                  is_lao, order_mom, num_cents, idx_cent,          &
                                  order_cent, num_geo_derv, prop_int)
          ! assigns the returned integrals, and write the integrals on file if required
          if (get_int) then
            call gen1int_shell_int_square(prop_int, base_geo_derv, num_ao_orb, &
                                          wrt_int, val_ints)
          end if
          ! calculates the expectation values
          if (alloc_expt) then
            call gen1int_shell_tr_square(prop_int, 0, num_ao_orb, kind_dens, &
                                         dim_dens, num_dens, ao_dens,        &
                                         tmp_expt)
          else if (get_expt) then
            call gen1int_shell_tr_square(prop_int, base_geo_derv, num_ao_orb,    &
                                         kind_dens, dim_dens, num_dens, ao_dens, &
                                         val_expt)
          end if
          ! cleans contracted property integrals
          call gen1int_shell_int_clean(prop_int)
        end do
      end do
    end select
    if (alloc_expt) then
      ! puts the expectation values from Gen1Int into appropriate positions
      if (redunt_expt) then
        call geom_total_redundant_expt(num_cents, idx_cent, order_cent, NUCDEP,  &
                                       num_dens*num_opt, num_geo_derv, tmp_expt, &
                                       (3*NUCDEP)**order_geo, val_expt)
      end if
!FIXME: the path and dimensions of expectation values should already be created
      ! writes \var(tmp_expt) into file
      if (wrt_expt) then

      end if
      deallocate(tmp_expt)
    ! writes parts of \var(val_expt) into file
    else if (wrt_expt) then

    end if
    call QEXIT("gen1int_ifc_prop")
    return
100 format("gen1int_ifc_prop>> ",A,I10)
#else
    call QUIT("gen1int_ifc_prop>> Gen1Int is not installed!")
#endif
  end subroutine gen1int_ifc_prop
