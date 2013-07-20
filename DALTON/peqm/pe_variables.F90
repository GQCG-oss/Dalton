!
!   Polarizable Embedding (PE) library
!   Copyright (C) 2013 Jógvan Magnus Haugaard Olsen
!
!   This file is part of the PE library.
!
!   The PE library is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License as
!   published by the Free Software Foundation, either version 3 of the
!   License, or (at your option) any later version.
!
!   The PE library is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with the PE library. If not, see <http://www.gnu.org/licenses/>.
!
!   Contact information:
!
!   Jógvan Magnus Haugaard Olsen
!   E-mail: foeroyingur@gmail.com
!
module pe_variables

    use pe_precision

#if defined(VAR_MPI)
#if defined(VAR_USE_MPIF)
    implicit none
#include "mpif.h"
#else
    use mpi
    implicit none
#endif
#endif

    ! options
    logical, save :: peqm = .false.
    logical, save :: pe_iter = .true.
    logical, save :: pe_border = .false.
    logical, save :: pe_damp = .false.
    logical, save :: pe_gspol = .false.
    logical, save :: pe_nomb = .false.
    logical, save :: pe_polar = .false.
    logical, save :: pe_cube = .false.
    logical, save :: pe_restart = .false.
    logical, save :: pe_verbose = .false.
    logical, save :: pe_debug = .false.

    ! calculation type
    logical, save :: fock = .false.
    logical, save :: energy = .false.
    logical, save :: response = .false.

    ! temporary solution for work array thing
    real(dp), dimension(:), pointer :: work

    ! filenames
    character(len=80) :: potfile = 'POTENTIAL.INP'

    ! MPI stuff
#if defined(VAR_MPI)
    integer, parameter :: comm = MPI_COMM_WORLD
    integer, parameter :: impi = MPI_INTEGER
    integer, parameter :: rmpi = MPI_REAL8
    integer, parameter :: lmpi = MPI_LOGICAL
#endif
    integer, save :: myid, nprocs, ierr
    integer, save :: site_start, site_finish
    integer, save :: cube_start, cube_finish
    logical, save :: synced = .false.
    integer, dimension(:), save, allocatable :: siteloops
    integer, dimension(:), save, allocatable :: cubeloops
    integer, dimension(:), save, allocatable :: poldists, sitedists
    integer, dimension(:), save, allocatable :: cubedists
    integer, dimension(:), save, allocatable :: displs

    ! logical unit for output file (default is stdout)
    integer, save :: luout = 6

    ! constants, thresholds and stuff
    ! 1 bohr = 0.5291772109217 Aa (codata 2010)
    real(dp), parameter :: aa2au = 1.0 / 0.5291772109217
    real(dp), parameter :: aa2au2 = aa2au * aa2au
    real(dp), parameter :: pi = 3.141592653589793
    real(dp), parameter :: zero = 1.0d-8
    integer, save :: scfcycle = 0
    real(dp), save :: thriter = 1.0d-8
    real(dp), save :: damp = 2.1304
    real(dp), save :: Rmin = 2.2
    integer, save :: nredist = 1
    character(len=6), save :: border_type = 'REDIST'
    ! use Cholesky factorization of classical response matrix
    logical, save :: chol = .true.

    ! C. E. Dykstra, J. Comp. Chem., 9 (1988), 476
    ! C^(n)_ij coefficients for calculating T(k) tensor elements
    integer, dimension(:,:,:), allocatable, save :: Cnij

    ! polarizable embedding potential info
    ! ------------------------------------
    ! total number of classical sites
    integer, save :: nsites = 0
    ! number of polarizable sites
    integer, save :: npols = 0
    ! exclusion list length
    integer, save :: lexlst = 0
    ! number of density matrices
    integer :: ndens = 0
    ! number of basis functions in core fragment
    integer, save :: nbas
    ! size of packed matrices
    integer, save :: nnbas
    ! size of full matrices
    integer, save :: n2bas
    ! number of nuclei in core region
    integer, save :: qmnucs = 0

    ! specifies what type of parameters are present
    ! lmul(0): monopoles, lmul(1): dipoles etc.
    logical, dimension(0:5), save :: lmul = .false.
    ! lpol(1): (an)isotropic dipole-dipole polarizabilities
    logical, dimension(1), save :: lpol = .false.
    ! lhypol(1): dipole-dipole-dipole polarizabilities/1st hyperpolarizability
!    logical, dimension(1), save :: lhypol

    ! charges, areas, coordinates, elements and exclusion lists
    ! site nuclear charges
    real(dp), dimension(:,:), allocatable, save :: Zs
    ! core fragment nuclear charges
    real(dp), dimension(:,:), allocatable, save :: Zm
    ! site coordinates
    real(dp), dimension(:,:), allocatable, save :: Rs
    ! core fragment nuclear coordinates
    real(dp), dimension(:,:), allocatable, save :: Rm
    ! site elements
    character(len=2), dimension(:), allocatable, save :: elems
    ! exclusion list
    integer, dimension(:,:), allocatable, save :: exclists

    ! energy contributions
    ! total
    real(dp), dimension(:), allocatable, save :: Epe
    ! electrostatic
    real(dp), dimension(:,:), allocatable, save :: Ees
    ! polarization
    real(dp), dimension(:,:), allocatable, save :: Epol

    ! multipole moments
    ! order of the highest order multipole moment
    integer, save :: mulorder = -1
    ! monopoles, dipoles, quadrupoles, octopoles, etc.
    real(dp), dimension(:,:), allocatable, save :: M0s, M1s, M2s, M3s, M4s, M5s
    ! (hyper)polarizabilities
    ! order of highest order polarizability
    integer, save :: polorder = -1
    ! dipole-dipole polarizabilities
    real(dp), dimension(:,:), allocatable, save :: P1s
    ! .true. if P1 > 0 else .false.
    logical, dimension(:), allocatable, save :: zeroalphas


    ! CUBE stuff
    ! ---------
    ! calculate electric field
    logical, save :: cube_field = .false.
    ! general cube information
    ! number of grid points
    integer, save :: npoints
    ! grid points
    real(dp), dimension(:,:), allocatable, save :: Rp
    ! CUBE file origin and step sizes
    real(dp), dimension(3), save :: origin, step
    ! grid density in x, y and z direction
    integer, save :: xgrid = 6
    integer, save :: ygrid = 6
    integer, save :: zgrid = 6
    ! numberof steps in x, y and z direction
    integer, save :: xsteps
    integer, save :: ysteps
    integer, save :: zsteps
    ! box size relative to molecule size
    real(dp), save :: xsize = 8.0
    real(dp), save :: ysize = 8.0
    real(dp), save :: zsize = 8.0

end module pe_variables
