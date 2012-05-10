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
!...  This file generates the cube file of the electron density and/or molecular
!...  orbitals using Gen1Int library.
!
!...  2012-03-09, Bin Gao
!...  * first version

#include "xkind.h"
#include "max_len_str.h"

!> \brief module of generating cube files
module gen1int_cube

  ! Dalton AO sub-shells
  use dalton_shell

  implicit none

  logical, save, public :: do_density_cube = .false.                !electron density cube file
  logical, save, public :: do_homo_cube = .false.                   !HOMO cube file
  logical, save, public :: do_lumo_cube = .false.                   !LUMO cube file
  logical, save, public :: do_mo_cube = .false.                     !MO cube files
  integer, save, public :: num_cube_mo = 0                          !number of MOs to generate
  integer, save, allocatable, public :: idx_cube_mo(:)              !indices of MOs
  character(MAX_LEN_STR), save, public :: cube_format = "GAUSSIAN"  !format of cube file
  integer, save, public :: num_cube_points = 0                      !number of points in cube file
  real(REALK), save, public :: cube_origin(3) = 0.0_REALK           !origin of cube file
  real(REALK), save, public :: cube_increment(3,3) = 0.0_REALK      !increments of cube file
  integer, save, public :: cube_num_inc(3) = 0                      !number of increments of cube file
  real(REALK), save, allocatable, public :: cube_coord(:,:)         !XYZ coordinates of points in cube file

  public :: Gen1IntCubeCreate

  contains

  !> \brief generates the cube file of the electron density and/or molecular orbitals
  !>        using Gen1Int library
  !> \author Bin Gao
  !> \date 2012-03-09
  !> \param len_work is the length of Dalton/Dirac workspace
  !> \param wrk_space is the Dalton/Dirac workspace
  !> \param io_viewer is the IO unit of standard output
  subroutine Gen1IntCubeCreate(len_work, wrk_space, io_viewer)
    integer, intent(in) :: len_work
    real(REALK), intent(inout) :: wrk_space(len_work)
    integer, intent(in) :: io_viewer
    ! uses MXCENT, etc.
#include "mxcent.h"
    ! uses NUCDEP, CHARGE, CORD
#include "nuclei.h"
    ! uses NCMOT, NNASHX, N2BASX
#include "inforb.h"
    ! uses LUSIFC
#include "inftap.h"
    integer ipoint                                    !incremental recorder over points
    integer ix, iy, iz                                !incremental recorders along XYZ directions
    integer num_ao                                    !number of atomic orbitals from Gen1Int interface
    integer size_ao                                   !size of matrices in atomic orbitals
    type(matrix) ao_dens(1)                           !AO density matrix
    logical found                                     !if found required data from SIRIFC
    type(one_prop_t) prop_operator                    !property operator
    real(REALK), allocatable :: cube_values(:,:,:,:)  !values of points in cube file
    real(REALK), allocatable :: mo_coef(:,:)          !molecular orbital coefficients
    integer start_ao, end_ao                          !start and end addresses of AOs
    integer imo                                       !incremental recorder over MOs
    integer io_cube                                   !IO unit of cube file
    integer ierr                                      !error information
    call QENTER("Gen1IntCubeCreate")
    if (allocated(cube_coord)) then
      ! gets the number of atomic orbitals from Gen1Int interface
      write(io_viewer,"()")
      call DaltonShellGetNumAO(num_ao=num_ao)
      write(io_viewer,100) "number of orbitals from Gen1Int interface:", num_ao
      ! electron density
      if (do_density_cube) then
        ! gets AO density matrix
        size_ao = num_ao*num_ao
        call get_ao_dens(size_ao, wrk_space(1:size_ao), len_work-size_ao, &
                         wrk_space(size_ao+1))
!FIXME: could be triangular and symmetric
        call MatAssociate(work_alpha=wrk_space(1:size_ao), &
                          num_row=num_ao, A=ao_dens(1),    &
                          info_mat=ierr)
        if (ierr/=0) call QUIT("Failed to associate AO density matrix!")
#define LEVEL_PRINT 1
        ! writes matrix to check
        if (LEVEL_PRINT>10) then
          write(io_viewer,"()")
          write(io_viewer,100) "AO density matrix"
          call MatView(A=ao_dens(1), io_viewer=io_viewer)
        end if
        ! initializes the information of property operator
        call OnePropCreate(prop_name=INT_OVERLAP_DIST, one_prop=prop_operator, &
                           info_prop=ierr, grid_points=cube_coord)
        if (ierr/=0) call QUIT("Failed to creat property operator!")
        ! evaluates the electron density at points of cube file
        allocate(cube_values(cube_num_inc(3),cube_num_inc(2),cube_num_inc(1),1), &
                 stat=ierr)
        if (ierr/=0) call QUIT("Failed to allocate cube_values!")
        cube_values = 0.0_REALK  !necessary to zero
        call DaltonShellIntegral(comp_bra=1,                  &
                                 comp_ket=1,                  &
                                 one_prop=prop_operator,      &
                                 num_ints=num_cube_points,    &
                                 num_dens=1, ao_dens=ao_dens, &
                                 val_expt=cube_values,        &
                                 io_viewer=io_viewer,         &
                                 level_print=LEVEL_PRINT)
        call OnePropDestroy(one_prop=prop_operator)
        ! writes cube file
        write(io_viewer,100) "writes cube file of electron density"
        io_cube = -1
        call GPOPEN(io_cube, "density.cube", "unknown", " ", "formatted", &
                    ierr, .false.)
        write(io_cube,"(1X,A)") "molecule density=scf"
        write(io_cube,"(1X,A)") "electron density from total SCF density"
        write(io_cube,"(I5,3F12.6)") NUCDEP, cube_origin
        do ix = 1, 3
          write(io_cube,"(I5,3F12.6)") cube_num_inc(ix), cube_increment(ix,:)
        end do
        do ipoint = 1, NUCDEP
          write(io_cube,"(I5,4F12.6)") int(CHARGE(ipoint)), CHARGE(ipoint), &
                                       CORD(:,ipoint)
        end do
        do iz = 1, cube_num_inc(1)
          do iy = 1, cube_num_inc(2)
            write(io_cube,"(6Es13.5)") cube_values(:,iy,iz,1)
          end do
        end do
        call GPCLOSE(io_cube, "KEEP")
        ! cleans
        deallocate(cube_values)
        call MatNullify(A=ao_dens(1))
      end if
      ! molecular orbitals
      if (do_homo_cube .or. do_lumo_cube .or. do_mo_cube) then
        ! gets molecular orbital coefficienct from SIRIFC file
        if (LUSIFC<=0) &
          call GPOPEN(LUSIFC, "SIRIFC", "OLD", " ", "UNFORMATTED", ierr, .false.)
        rewind(LUSIFC)
        ! reads the molecular orbital coefficients
#if REALK == 4
        call SZERO(wrk_space(1), NCMOT)
#elif REALK == 8
        call DZERO(wrk_space(1), NCMOT)
#else
        call QUIT("Unknown kind of real numbers!")
#endif
#ifdef PRG_DIRAC
        print *, 'error: RD_SIRIFC not available in DIRAC'
        stop 1
#else
        call RD_SIRIFC("CMO", found, wrk_space(1), wrk_space(NCMOT+1), &
                       len_work-NCMOT)
#endif
        if (.not.found) call QUIT("CMO IS NOT FOUND ON SIRIFC!")
        ! gets the required MOs
        if (do_homo_cube) num_cube_mo = num_cube_mo+1
        if (do_lumo_cube) num_cube_mo = num_cube_mo+1
        allocate(mo_coef(num_ao,num_cube_mo), stat=ierr)
        if (ierr/=0) call QUIT("Failed to allocate mo_coef!")
        if (do_mo_cube) then
          do imo = 1, size(idx_cube_mo)
            if (idx_cube_mo(imo)>NCMOT .or. idx_cube_mo(imo)<1) &
              call QUIT("Incorrect MOs!")
            start_ao = (idx_cube_mo(imo)-1)*num_ao
            end_ao = start_ao+num_ao
            start_ao = start_ao+1
            mo_coef(:,imo) = wrk_space(start_ao:end_ao)
          end do
          imo = size(idx_cube_mo)
        else
          imo = 0
        end if
        if (do_homo_cube) then
          start_ao = (NOCCT-1)*num_ao
          end_ao = start_ao+num_ao
          start_ao = start_ao+1
          imo = imo+1
          mo_coef(:,imo) = wrk_space(start_ao:end_ao)
        end if
        if (do_lumo_cube) then
          start_ao = NOCCT*num_ao
          end_ao = start_ao+num_ao
          start_ao = start_ao+1
          imo = imo+1
          mo_coef(:,imo) = wrk_space(start_ao:end_ao)
        end if
        ! evaluates MOs at points in cube file
        allocate(cube_values(cube_num_inc(3),cube_num_inc(2), &
                             cube_num_inc(1),num_cube_mo),    &
                 stat=ierr)
        if (ierr/=0) call QUIT("Failed to allocate cube_values!")
        call DaltonShellMO(num_mo=num_cube_mo, num_ao=num_ao, mo_coef=mo_coef, &
                           num_points=num_cube_points, grid_points=cube_coord, &
                           num_derv=1, val_mo=cube_values)
        deallocate(mo_coef)
        ! resets the number of indices of MOs in cube file
        if (do_mo_cube) then
          num_cube_mo = size(idx_cube_mo)
        else
          num_cube_mo = 0
        end if
        ! writes cube file of MOs
        write(io_viewer,100) "writes cube file of HOMO, LUMO and/or MOs"
        if (do_mo_cube) then
          io_cube = -1
          call GPOPEN(io_cube, "mo.cube", "unknown", " ", "formatted", &
                      ierr, .false.)
          write(io_cube,"(1X,A)") "molecule mo=selected"
          write(io_cube,"(1X,A)") "MO coefficients"
          write(io_cube,"(I5,3F12.6)") -NUCDEP, cube_origin
          do ix = 1, 3
            write(io_cube,"(I5,3F12.6)") cube_num_inc(ix), cube_increment(ix,:)
          end do
          do ipoint = 1, NUCDEP
            write(io_cube,"(I5,4F12.6)") int(CHARGE(ipoint)), CHARGE(ipoint), &
                                         CORD(:,ipoint)
          end do
          write(io_cube,"(10I5)") num_cube_mo, idx_cube_mo
          do iz = 1, cube_num_inc(1)
            do iy = 1, cube_num_inc(2)
!FIXME: the memory access here is awful, might need to change
              write(io_cube,"(6Es13.5)") &
                ((cube_values(ix,iy,iz,imo), imo=1,num_cube_mo), ix=1,cube_num_inc(3))
            end do
          end do
          call GPCLOSE(io_cube, "KEEP")
          imo = num_cube_mo
        else
          imo = 0
        end if
        ! writes cube file of HOMO
        if (do_homo_cube) then
          io_cube = -1
          call GPOPEN(io_cube, "homo.cube", "unknown", " ", "formatted", &
                      ierr, .false.)
          write(io_cube,"(1X,A)") "molecule mo=HOMO"
          write(io_cube,"(1X,A)") "MO coefficients"
          write(io_cube,"(I5,3F12.6)") -NUCDEP, cube_origin
          do ix = 1, 3
            write(io_cube,"(I5,3F12.6)") cube_num_inc(ix), cube_increment(ix,:)
          end do
          do ipoint = 1, NUCDEP
            write(io_cube,"(I5,4F12.6)") int(CHARGE(ipoint)), CHARGE(ipoint), &
                                         CORD(:,ipoint)
          end do
          write(io_cube,"(10I5)") 1, NOCCT
          imo = imo+1
          do iz = 1, cube_num_inc(1)
            do iy = 1, cube_num_inc(2)
              write(io_cube,"(6Es13.5)") cube_values(:,iy,iz,imo)
            end do
          end do
          call GPCLOSE(io_cube, "KEEP")
        end if
        ! writes cube file of LUMO
        if (do_lumo_cube) then
          io_cube = -1
          call GPOPEN(io_cube, "lumo.cube", "unknown", " ", "formatted", &
                      ierr, .false.)
          write(io_cube,"(1X,A)") "molecule mo=LUMO"
          write(io_cube,"(1X,A)") "MO coefficients"
          write(io_cube,"(I5,3F12.6)") -NUCDEP, cube_origin
          do ix = 1, 3
            write(io_cube,"(I5,3F12.6)") cube_num_inc(ix), cube_increment(ix,:)
          end do
          do ipoint = 1, NUCDEP
            write(io_cube,"(I5,4F12.6)") int(CHARGE(ipoint)), CHARGE(ipoint), &
                                         CORD(:,ipoint)
          end do
          write(io_cube,"(10I5)") 1, NOCCT+1
          imo = imo+1
          do iz = 1, cube_num_inc(1)
            do iy = 1, cube_num_inc(2)
              write(io_cube,"(6Es13.5)") cube_values(:,iy,iz,imo)
            end do
          end do
          call GPCLOSE(io_cube, "KEEP")
        end if
        deallocate(cube_values)
      end if
    else
      write(io_viewer,"()")
      write(io_viewer,100) "points of cube file are not available"
      write(io_viewer,100) "no cube file will be generated"
      write(io_viewer,"()")
    end if
    ! cleans
    do_density_cube = .false.
    do_homo_cube = .false.
    do_lumo_cube = .false.
    do_mo_cube = .false.
    num_cube_mo = 0
    if (allocated(idx_cube_mo)) deallocate(idx_cube_mo)
    cube_format = "GAUSSIAN"
    num_cube_points = 0
    cube_origin = 0.0_REALK
    cube_increment = 0.0_REALK
    cube_num_inc = 0
    if (allocated(cube_coord)) deallocate(cube_coord)
    call QEXIT("Gen1IntCubeCreate")
100 format("Gen1IntCubeCreate>> ",A,I8)
  end subroutine Gen1IntCubeCreate

end module gen1int_cube
