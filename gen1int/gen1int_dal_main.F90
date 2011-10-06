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
!...  This file is the main driver of calling Gen1Int interface from Dalton,
!...  in which the integral matrices are symmetric
!
!...  2011-02-13, Bin Gao
!...  * this file does not read information from input anymore, it will only
!...    be called by other subroutines
!
!...  2010-07-28, Bin Gao
!...  * first version

#include "xkind.h"

  !> \brief main driver of calling Gen1Int interface from Dalton, in which
  !>        the integral matrices are symmetric
  !> \author Bin Gao
  !> \date 2010-07-28
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
  !> \param max_num_cent is the maximum number of geometric differentiated centers
  !> \param order_geo_total is the order of total geometric derivatives
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
  subroutine gen1int_dal_main(prop_name, is_lao,                             &
                              order_mag_bra, order_mag_ket, order_mag_total, &
                              order_ram_bra, order_ram_ket, order_ram_total, &
                              order_geo_bra, order_geo_ket,                  &
                              max_num_cent, order_geo_total,                 &
                              sym_int, get_int, wrt_int, dim_int, vals_int,  &
                              do_expt, num_dens, ao_dens, get_expt,          &
                              wrt_expt, vals_expt, io_unit, level_print)
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
    integer, intent(in) :: max_num_cent
    integer, intent(in) :: order_geo_total
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
    ! number of operators, magnetic derivatives, derivatives w.r.t. total
    ! rotational angular momentum, and partial geometric derivatives
    integer num_opt_derv
    ! indices of atomic centers
    integer, allocatable :: idx_cent(:)
    ! order of derivatives of the corresponding atomic centers
    integer, allocatable :: order_cent(:)
    ! selected atom nodes
    integer, allocatable :: idx_node(:)
    ! weights of the selected atom nodes
    integer, allocatable :: wt_node(:)
    ! number of different paths for total geometric derivatives
    integer num_paths
    ! depth of atom to visit
    integer visit_depth
    ! number of geometric derivatives per path
    integer num_geo_cent
    ! start index of geometric derivatives per path
    integer start_idx_geo
    ! end index of geometric derivatives per path
    integer end_idx_geo
    ! start index of operators per path
    integer start_idx_opt
    ! end index of operators per path
    integer end_idx_opt
    ! incremental recorder over different paths of total geometric derivatives
    integer ipath
    ! incremental recorder over centers
    integer icent
    ! error information
    integer ierr
    ! Dalton stuff
#include "mxcent.h"
#include "nuclei.h"
    real(REALK) TIMHER, TEND
    real(REALK) WALHER, WEND
    call QENTER("gen1int_dal_main")
    call GETTIM(TIMHER, WALHER)
    call HEADER("Gen1Int calculates "//trim(prop_name), -1)
    ! checks if the Gen1Int interface is initialized
    call gen1int_dal_init
    ! dumps to check
    if (level_print>=5) then
      if (is_lao) write(io_unit,100) "using London atomic orbitals"
      write(io_unit,100) "order of magnetic derivaives on bra center", order_mag_bra
      write(io_unit,100) "order of magnetic derivaives on ket center", order_mag_ket
      write(io_unit,100) "order of total magnetic derivaives", order_mag_total
      write(io_unit,100) "order of derivatives w.r.t. total rotational "// &
                         "angular momentum on bra center", order_ram_bra
      write(io_unit,100) "order of derivatives w.r.t. total rotational "// &
                         "angular momentum on ket center", order_ram_ket
      write(io_unit,100) "order of total derivatives w.r.t. total "// &
                         "rotational angular momentum", order_ram_total
      write(io_unit,100) "order of geometric derivatives on bra center", order_geo_bra
      write(io_unit,100) "order of geometric derivatives on ket center", order_geo_ket
      write(io_unit,100) "maximum number of differentiated centers", max_num_cent
      write(io_unit,100) "order of total geometric derivatives", order_geo_total
      if (sym_int) then
        write(io_unit,100) "symmetric integral matrices"
      else
        write(io_unit,100) "non-symmetric integral matrices"
      end if
      if (get_int) write(io_unit,100) "return integrals"
      if (wrt_int) write(io_unit,100) "write integrals on file"
      write(io_unit,110) "dimension of integral matrices", dim_int
      if (do_expt) then
        write(io_unit,100) "calculate expectation values"
        write(io_unit,100) "number of AO density matrices", num_dens
        if (get_expt) write(io_unit,100) "return expectation values"
        if (wrt_expt) write(io_unit,100) "write expectation values on file"
      end if
    end if
    ! calculates geometric derivatives
    if (max_num_cent>0 .and. order_geo_total>0) then
      ! number of operators, magnetic derivatives, derivatives w.r.t. total
      ! rotational angular momentum, and partial geometric derivatives
      num_opt_derv = (order_mag_bra+1)*(order_mag_bra+2)     &
                   * (order_mag_ket+1)*(order_mag_ket+2)     &
                   * (order_mag_total+1)*(order_mag_total+2) &
                   * (order_ram_bra+1)*(order_ram_bra+2)     &
                   * (order_ram_ket+1)*(order_ram_ket+2)     &
                   * (order_ram_total+1)*(order_ram_total+2) &
                   * (order_geo_bra+1)*(order_geo_bra+2)     &
                   * (order_geo_ket+1)*(order_geo_ket+2)/256
      select case(trim(prop_name))
      case("OVERLAP")
      case default
        call QUIT("gen1int_dal_main>> "//trim(prop_name)//" is not implemented!")
      end select
      allocate(idx_cent(max_num_cent), stat=ierr)
      if (ierr/=0) call QUIT("gen1int_dal_main>> failed to allocate idx_cent!")
      allocate(order_cent(max_num_cent), stat=ierr)
      if (ierr/=0) call QUIT("gen1int_dal_main>> failed to allocate order_cent!")
      allocate(idx_node(order_geo_total), stat=ierr)
      if (ierr/=0) call QUIT("gen1int_dal_main>> failed to allocate idx_node!")
      allocate(wt_node(order_geo_total), stat=ierr)
      if (ierr/=0) call QUIT("gen1int_dal_main>> failed to allocate wt_node!")
      ! computes the total number of different paths, and generates
      ! the first path (composition of centers)
      call geom_total_tree_init(NUCDEP, order_geo_total, max_num_cent,     &
                                num_paths, visit_depth, idx_node, wt_node, &
                                idx_cent, order_cent, num_geo_cent)
      start_idx_geo = 1
      end_idx_geo = num_geo_cent
      ! dumps to check
      if (level_print>=10)                                                     &
        write(io_unit,110) "number of compositions of differentiated centers", &
                           num_paths
      if (level_print>=20) then
        write(io_unit,120) ipath, start_idx_geo, wt_node(order_geo_total), &
                           (idx_cent(icent),"(",order_cent(icent),         &
                           ")",icent=1,wt_node(order_geo_total))
      end if
      ! calculates the geometric and magnetic derivatives
      start_idx_opt = num_opt_derv*(start_idx_geo-1)+1
      end_idx_opt = num_opt_derv*end_idx_geo
      call gen1int_dal_prop(prop_name, is_lao,                              &
                            order_mag_bra, order_mag_ket, order_mag_total,  &
                            order_ram_bra, order_ram_ket, order_ram_total,  &
                            order_geo_bra, order_geo_ket,                   &
                            wt_node(order_geo_total),                       &
                            idx_cent(1:wt_node(order_geo_total)),           &
                            order_cent(1:wt_node(order_geo_total)),         &
                            sym_int, get_int, wrt_int, dim_int,             &
                            vals_int(:,start_idx_opt:end_idx_opt),          &
                            do_expt, num_dens, ao_dens, get_expt, wrt_expt, &
                            vals_expt(:,start_idx_opt:end_idx_opt),         &
                            io_unit, level_print)
      if (level_print>=10) then
        write(io_unit,110) "current number of geometric derivatives", &
                           end_idx_geo
      end if
      ! loops over other paths
      do ipath = 2, num_paths
        ! generates the differentiated centers and their orders
        call geom_total_tree_search(NUCDEP, order_geo_total, max_num_cent,    &
                                    visit_depth, idx_node, wt_node, idx_cent, &
                                    order_cent, num_geo_cent)
        ! computes the position of current geometric derivatives
        start_idx_geo = end_idx_geo
        end_idx_geo = start_idx_geo+num_geo_cent
        start_idx_geo = start_idx_geo+1
        ! dumps to check
        if (level_print>=20) then
          write(io_unit,120) ipath, start_idx_geo, wt_node(order_geo_total), &
                             (idx_cent(icent),"(",order_cent(icent),         &
                             ")",icent=1,wt_node(order_geo_total))
        end if
        ! calculates the geometric and magnetic derivatives
        start_idx_opt = num_opt_derv*(start_idx_geo-1)+1
        end_idx_opt = num_opt_derv*end_idx_geo
        call gen1int_dal_prop(prop_name, is_lao,                              &
                              order_mag_bra, order_mag_ket, order_mag_total,  &
                              order_ram_bra, order_ram_ket, order_ram_total,  &
                              order_geo_bra, order_geo_ket,                   &
                              wt_node(order_geo_total),                       &
                              idx_cent(1:wt_node(order_geo_total)),           &
                              order_cent(1:wt_node(order_geo_total)),         &
                              sym_int, get_int, wrt_int, dim_int,             &
                              vals_int(:,start_idx_opt:end_idx_opt),          &
                              do_expt, num_dens, ao_dens, get_expt, wrt_expt, &
                              vals_expt(:,start_idx_opt:end_idx_opt),         &
                              io_unit, level_print)
        if (level_print>=10) then
          write(io_unit,110) "current number of geometric derivatives", &
                             end_idx_geo
        end if
      end do
      deallocate(idx_cent)
      deallocate(order_cent)
      deallocate(idx_node)
      deallocate(wt_node)
    ! no geometric derivatives
    else
      call gen1int_dal_prop(prop_name, is_lao,                             &
                            order_mag_bra, order_mag_ket, order_mag_total, & 
                            order_ram_bra, order_ram_ket, order_ram_total, & 
                            order_geo_bra, order_geo_ket,                  & 
                            0, (/0/), (/0/),                               &
                            sym_int, get_int, wrt_int, dim_int, vals_int,  & 
                            do_expt, num_dens, ao_dens, get_expt,          & 
                            wrt_expt, vals_expt, io_unit, level_print)
    end if
    ! Dalton stuff
    call GETTIM(TEND, WEND)
    TIMHER = TEND-TIMHER
    WALHER = WEND-WALHER
    write(io_unit,"()")
    call TIMTXT(">>>> Total CPU  time used in Gen1Int", TIMHER, io_unit)
    call TIMTXT(">>>> Total wall time used in Gen1Int", WALHER, io_unit)
    call HEADER("End of Gen1Int", -1)
    call QEXIT("gen1int_dal_main")
    return
100 format("gen1int_dal_main>> ",A,I4)
110 format("gen1int_dal_main>> ",A,I12)
120 format("gen1int_dal_main>> ","path",I8,", start",I10,", N_{cent}",I3, &
           ", atom(order)",20(I3,A,I2,A))
#else
    call QUIT("Gen1Int is not installed!")
#endif
  end subroutine gen1int_dal_main
