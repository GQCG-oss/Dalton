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
!...  This file calculates property integrals and/or expectation values, in which
!...  the integral matrices are symmetric
!
!...  2011-10-04, Bin Gao
!...  * first version

#include "xkind.h"

  !> \brief calculates property integrals and/or expectation values, in which
  !>        the integral matrices are symmetric
  !> \author Bin Gao
  !> \date 2011-10-04
  !> \param prop_name is the name of property integrals to calculated
  !> \param is_lao indicates if using rotational London atomic orbitals
  !> \param order_mag_bra is the order of magnetic derivatives on bra center
  !> \param order_mag_ket is the order of magnetic derivatives on ket center
  !> \param order_mag_total is the order of total magnetic derivatives
  !> \param order_ram_bra is the order of derivatives w.r.t. total rotational angular momentum on bra center
  !> \param order_ram_ket is the order of derivatives w.r.t. total rotational angular momentum on ket center
  !> \param order_ram_total is the order of total derivatives w.r.t. total rotational angular momentum
  !> \param order_geo_bra is the order of geometric derivatives with respect to bra center
  !> \param order_geo_ket is the order of geometric derivatives with respect to ket center
  !> \param num_cents is the number of geometric differentiated centers
  !> \param idx_cent contains the indices of differentiated centers
  !> \param order_cent contains the orders of differentiated centers
  !> \param sym_int indicates if the integral matrices are symmetric
  !> \param get_int indicates if getting integrals back
  !> \param wrt_int indicates if writing integrals to file
  !> \param dim_int is the dimension of integral matrices
  !> \param do_expt indicates if calculating expectation values
  !> \param num_dens is the number of AO density matrices
  !> \param ao_dens contains the AO density matrices
  !> \param get_expt indicates if getting expectation values back
  !> \param wrt_expt indicates if writing expectation values to file
  !> \param io_unit is the IO unit of standard output
  !> \param level_print is the level of print
  !> \return vals_int contains the calculated integrals
  !> \return vals_expt contains the calculated expectation values
  subroutine gen1int_dal_prop(prop_name, is_lao,                              &
                              order_mag_bra, order_mag_ket, order_mag_total,  &
                              order_ram_bra, order_ram_ket, order_ram_total,  &
                              order_geo_bra, order_geo_ket,                   &
                              num_cents, idx_cent, order_cent,                &
                              sym_int, get_int, wrt_int, dim_int, vals_int,   &
                              do_expt, num_dens, ao_dens, get_expt, wrt_expt, &
                              vals_expt, io_unit, level_print)
    ! AO sub-shells
    use gen1int_shell
    implicit none
    character*(*), intent(in) :: prop_name
    logical, intent(in) :: is_lao
    integer, intent(in) :: order_mag_bra
    integer, intent(in) :: order_mag_ket
    integer, intent(in) :: order_mag_total
    integer, intent(in) :: order_ram_bra
    integer, intent(in) :: order_ram_ket
    integer, intent(in) :: order_ram_total
    integer, intent(in) :: order_geo_bra
    integer, intent(in) :: order_geo_ket
    integer, intent(in) :: num_cents
    integer, intent(in) :: idx_cent(num_cents)
    integer, intent(in) :: order_cent(num_cents)
    logical, intent(in) :: sym_int
    logical, intent(in) :: get_int
    logical, intent(in) :: wrt_int
    integer, intent(in) :: dim_int
    real(REALK), intent(out) :: vals_int(dim_int,*)
    logical, intent(in) :: do_expt
    integer, intent(in) :: num_dens
    real(REALK), intent(in) :: ao_dens(dim_int,num_dens)
    logical, intent(in) :: get_expt
    logical, intent(in) :: wrt_expt
    real(REALK), intent(out) :: vals_expt(num_dens,*)
    integer, intent(in) :: io_unit
    integer, intent(in) :: level_print
#ifdef BUILD_GEN1INT
! common block with origins
#include "orgcom.h"
    integer num_derv                      !number of derivatives
    integer num_opt                       !number of operators including derivatives
    integer num_ao_bra, num_ao_ket        !number of AOs on bra and ket centers
    integer num_contr_bra, num_contr_ket  !number of contractions on bra and ket centers
    integer num_orb_bra, num_orb_ket      !number of orbitals on bra and ket centers for
                                          !the contracted integrals
    real(REALK), allocatable :: contr_ints(:,:,:)
                                          !contracted integrals between sub-shells on bra and ket centers
    integer ishell, jshell                !incremental recorders over sub-shells
    integer ibra, iket, iopt              !incremental recorders over contracted integrals
    integer ival                          !incremental recorder over returned values
    integer ierr                          !error information
    call QENTER("gen1int_dal_prop")
    ! computes the number of derivatives
    num_derv = (order_mag_bra+1)*(order_mag_bra+2)     &
             * (order_mag_ket+1)*(order_mag_ket+2)     &
             * (order_mag_total+1)*(order_mag_total+2) &
             * (order_ram_bra+1)*(order_ram_bra+2)     &
             * (order_ram_ket+1)*(order_ram_ket+2)     &
             * (order_ram_total+1)*(order_ram_total+2) &
             * (order_geo_bra+1)*(order_geo_bra+2)     &
             * (order_geo_ket+1)*(order_geo_ket+2)/256
    do ishell = 1, num_cents
      num_derv = num_derv*(order_cent(ishell)+1)*(order_cent(ishell)+2)/2
    end do
    if (level_print>=10) write(io_unit,100) "number of derivatives", num_derv
    ! symmetric integral matrices
    if (sym_int) then
      ival = 0
      ! loops over AO sub-shells on ket center
      do jshell = 1, num_ao_shells
        ! gets the number of AOs and contractions on ket center
        call gen1int_shell_dims(ao_shells(jshell), num_ao_ket, num_contr_ket)
        num_orb_ket = num_ao_ket*num_contr_ket
        ! loops over AO sub-shells on bra center (lower part of integral matrix)
        do ishell = jshell, num_ao_shells
          ! gets the number of AOs and contractions on bra center
          call gen1int_shell_dims(ao_shells(ishell), num_ao_bra, num_contr_bra)
          num_orb_bra = num_ao_bra*num_contr_bra
          ! different property integrals
          select case(trim(prop_name))
          case("OVERLAP")
            num_opt = num_derv
            allocate(contr_ints(num_orb_bra,num_orb_ket,num_opt), stat=ierr)
            if (ierr/=0) call QUIT("gen1int_dal_prop>> failed to allocate contr_ints!")
            call gen1int_shell_carmom(ao_shells(ishell), ao_shells(jshell), is_lao,  &
                                      0, -1, DIPORG, 1.0_REALK, 0,                   &
                                      order_mag_bra, order_mag_ket, order_mag_total, &
                                      order_ram_bra, order_ram_ket, order_ram_total, &
                                      order_geo_bra, order_geo_ket, 0,               &
                                      num_cents, idx_cent, order_cent,               &
                                      num_ao_bra, num_contr_bra, num_ao_ket,         &
                                      num_contr_ket, num_opt, contr_ints)
          case default
            call QUIT("gen1int_dal_prop>> "//trim(prop_name)//" is not implemented!")
          end select
          ! assigns returned integrals
!FIXME
          if (get_int) then
            if (jshell<ishell) then
              do iopt = 1, num_opt
                do iket = 1, num_orb_ket
                  do ibra = 1, num_orb_bra
                    ival = ival+1
                    vals_int(ival,iopt) = contr_ints(ibra,iket,iopt)
                  end do
                end do
              end do
            else
              do iopt = 1, num_opt
                do iket = 1, num_orb_ket
                  do ibra = iket, num_orb_bra
                    ival = ival+1
                    vals_int(ival,iopt) = contr_ints(ibra,iket,iopt)
                  end do
                end do
              end do
            end if
          end if
          if (wrt_int) then
            if (jshell<ishell) then
            else
            end if
          end if
          if (do_expt) then
            if (jshell<ishell) then
            else
            end if
            if (get_expt) then
            end if
            if (wrt_expt) then
            end if
          end if
          deallocate(contr_ints)
        end do
      end do
    ! non-symmetric integral matrices
    else
      ival = 0
      ! loops over AO sub-shells on ket center
      do jshell = 1, num_ao_shells
        ! gets the number of AOs and contractions on ket center
        call gen1int_shell_dims(ao_shells(jshell), num_ao_ket, num_contr_ket)
        num_orb_ket = num_ao_ket*num_contr_ket
        ! loops over AO sub-shells on bra center
        do ishell = 1, num_ao_shells
          ! gets the number of AOs and contractions on bra center
          call gen1int_shell_dims(ao_shells(ishell), num_ao_bra, num_contr_bra)
          num_orb_bra = num_ao_bra*num_contr_bra
          ! different property integrals
          select case(trim(prop_name))
          case("OVERLAP")
            num_opt = num_derv
            allocate(contr_ints(num_orb_bra,num_orb_ket,num_opt), stat=ierr)
            if (ierr/=0) call QUIT("gen1int_dal_prop>> failed to allocate contr_ints!")
            call gen1int_shell_carmom(ao_shells(ishell), ao_shells(jshell), is_lao,  &
                                      0, -1, DIPORG, 1.0_REALK, 0,                   &
                                      order_mag_bra, order_mag_ket, order_mag_total, &
                                      order_ram_bra, order_ram_ket, order_ram_total, &
                                      order_geo_bra, order_geo_ket, 0,               &
                                      num_cents, idx_cent, order_cent,               &
                                      num_ao_bra, num_contr_bra, num_ao_ket,         &
                                      num_contr_ket, num_opt, contr_ints)
          case default
            call QUIT("gen1int_dal_prop>> "//trim(prop_name)//" is not implemented!")
          end select
          ! assigns returned integrals
!FIXME
          if (get_int) then
            do iopt = 1, num_opt
              do iket = 1, num_orb_ket
                do ibra = 1, num_orb_bra
                  ival = ival+1
                  vals_int(ival,iopt) = contr_ints(ibra,iket,iopt)
                end do
              end do
            end do
          end if
          if (wrt_int) then
          end if
          if (do_expt) then
            if (get_expt) then
            end if
            if (wrt_expt) then
            end if
          end if
          deallocate(contr_ints)
        end do
      end do
    end if
    call QEXIT("gen1int_dal_prop")
    return
100 format("gen1int_dal_prop>> ",A,I10)
#else
    call QUIT("Gen1Int is not installed!")
#endif
  end subroutine gen1int_dal_prop
