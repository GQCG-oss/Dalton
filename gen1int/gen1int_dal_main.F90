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
!...  This file is the main driver of calling Gen1Int interface from Dalton.
!
!...  2011-02-13, Bin Gao
!...  * this file does not read information from input anymore, it will only
!...    be called by other subroutines
!
!...  2010-07-28, Bin Gao
!...  * first version

#include "xkind.h"

  !> \brief main driver of calling Gen1Int interface from Dalton
  !> \author Bin Gao
  !> \date 2010-07-28
  !> \param prop_name is the name of property integrals to calculated
  !> \param is_lao indicates if using rotational London atomic orbitals
  !> \param order_mom is the order of Cartesian multipole moments, only used for
  !>        integrals "CARMOM"
  !> \param max_num_cent is the maximum number of geometric differentiated centers
  !> \param order_geo_total is the order of total geometric derivatives
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
  subroutine gen1int_dal_main(prop_name, is_lao, order_mom,  &
                              max_num_cent, order_geo_total, &
                              kind_int, get_int, wrt_int,    &
                              dim_int, val_ints,             &
                              do_expt, num_dens, ao_dens,    &
                              get_expt, wrt_expt, val_expt,  &
                              io_std, level_print)
    ! AO sub-shells
    use gen1int_shell
    implicit none
    character*(*), intent(in) :: prop_name
    logical, intent(in) :: is_lao
    integer, intent(in) :: order_mom
    integer, intent(in) :: max_num_cent
    integer, intent(in) :: order_geo_total
    integer, intent(in) :: kind_int
    logical, intent(in) :: get_int
    logical, intent(in) :: wrt_int
    integer, intent(in) :: dim_int
    real(REALK), intent(out) :: val_ints(dim_int,*)
    logical, intent(in) :: do_expt
    integer, intent(in) :: num_dens
    real(REALK), intent(in) :: ao_dens(dim_int,num_dens)
    logical, intent(in) :: get_expt
    logical, intent(in) :: wrt_expt
    real(REALK), intent(out) :: val_expt(num_dens,*)
    integer, intent(in) :: io_std
    integer, intent(in) :: level_print
#ifdef BUILD_GEN1INT
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
    ! number of total geometric derivatives per path
    integer num_geo_cent
    ! base address of total geometric derivatives per path
    integer base_geo_derv
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
    if (.not.shell_init) call QUIT("Gen1Int interface is not properly initialized!")
    ! dumps to check
    if (level_print>=5) then
      if (is_lao) write(io_std,100) "using London atomic orbitals"
      write(io_std,100) "maximum number of differentiated centers", max_num_cent
      write(io_std,100) "order of total geometric derivatives", order_geo_total
      select case(kind_int)
      case(1)
        write(io_std,100) "symmetric integral matrices"
      case(-1)
        write(io_std,100) "anti-symmetric integral matrices"
      case default
        write(io_std,100) "square integral matrices"
      end select
      if (get_int) write(io_std,100) "return integrals"
      if (wrt_int) write(io_std,100) "write integrals on file"
      write(io_std,110) "dimension of integral matrices", dim_int
      if (do_expt) then
        write(io_std,100) "calculate expectation values"
        write(io_std,100) "number of AO density matrices", num_dens
        if (get_expt) write(io_std,100) "return expectation values"
        if (wrt_expt) write(io_std,100) "write expectation values on file"
      end if
    end if
    ! calculates total geometric derivatives
    if (max_num_cent>0 .and. order_geo_total>0) then
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
      base_geo_derv = 0
      call geom_total_tree_init(NUCDEP, order_geo_total, max_num_cent,     &
                                num_paths, visit_depth, idx_node, wt_node, &
                                idx_cent, order_cent, num_geo_cent)
      ! dumps to check
      if (level_print>=10)                                                    &
        write(io_std,110) "number of compositions of differentiated centers", &
                           num_paths
      if (level_print>=20) then
        write(io_std,120) ipath, base_geo_derv, wt_node(order_geo_total), &
                          (idx_cent(icent),"(",order_cent(icent),         &
                          ")",icent=1,wt_node(order_geo_total))
      end if
      ! calculates the property integrals
      call gen1int_dal_prop(prop_name, is_lao, order_mom,            &
                            base_geo_derv, wt_node(order_geo_total), &
                            idx_cent(1:wt_node(order_geo_total)),    &
                            order_cent(1:wt_node(order_geo_total)),  &
                            kind_int, get_int, wrt_int, dim_int,     &
                            val_ints, do_expt, num_dens, ao_dens,    &
                            get_expt, wrt_expt, val_expt, io_std, level_print)
      if (level_print>=10) then
        write(io_std,110) "current number of total geometric derivatives", &
                           base_geo_derv+num_geo_cent
      end if
      ! loops over other paths
      do ipath = 2, num_paths
        ! updates the base address of current total geometric derivatives,
        ! before calling \fn(geom_total_tree_search)!
        base_geo_derv = base_geo_derv+num_geo_cent
        ! generates the differentiated centers and their orders
        call geom_total_tree_search(NUCDEP, order_geo_total, max_num_cent,    &
                                    visit_depth, idx_node, wt_node, idx_cent, &
                                    order_cent, num_geo_cent)
        ! dumps to check
        if (level_print>=20) then
          write(io_std,120) ipath, base_geo_derv, wt_node(order_geo_total), &
                            (idx_cent(icent),"(",order_cent(icent),         &
                            ")",icent=1,wt_node(order_geo_total))
        end if
        ! calculates the property integrals
        call gen1int_dal_prop(prop_name, is_lao, order_mom,            &
                              base_geo_derv, wt_node(order_geo_total), &
                              idx_cent(1:wt_node(order_geo_total)),    &
                              order_cent(1:wt_node(order_geo_total)),  &
                              kind_int, get_int, wrt_int, dim_int,     &
                              val_ints, do_expt, num_dens, ao_dens,    &
                              get_expt, wrt_expt, val_expt, io_std, level_print)
        if (level_print>=10) then
          write(io_std,110) "current number of total geometric derivatives", &
                             base_geo_derv+num_geo_cent
        end if
      end do
      deallocate(idx_cent)
      deallocate(order_cent)
      deallocate(idx_node)
      deallocate(wt_node)
    ! no total geometric derivatives
    else
      call gen1int_dal_prop(prop_name, is_lao, order_mom,         &
                            0, 0, (/0/), (/0/),                   &
                            kind_int, get_int, wrt_int, dim_int,  &
                            val_ints, do_expt, num_dens, ao_dens, &
                            get_expt, wrt_expt, val_expt, io_std, level_print)
    end if
    ! Dalton stuff
    call GETTIM(TEND, WEND)
    TIMHER = TEND-TIMHER
    WALHER = WEND-WALHER
    write(io_std,"()")
    call TIMTXT(">>>> Total CPU  time used in Gen1Int", TIMHER, io_std)
    call TIMTXT(">>>> Total wall time used in Gen1Int", WALHER, io_std)
    call HEADER("End of Gen1Int", -1)
    call QEXIT("gen1int_dal_main")
    return
100 format("gen1int_dal_main>> ",A,I4)
110 format("gen1int_dal_main>> ",A,I12)
120 format("gen1int_dal_main>> ","path",I8,", base",I10,", N_{cent}",I3, &
           ", atom(order)",20(I3,A,I2,A))
#else
    call QUIT("Gen1Int is not installed!")
#endif
  end subroutine gen1int_dal_main
