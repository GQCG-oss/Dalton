!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
#if defined(BUILD_PELIB)
module pelib_interface

    implicit none

    private

    public :: use_pelib, pelib_ifc_gspol
    public :: pelib_ifc_do_mep, pelib_ifc_do_mep_noqm, pelib_ifc_do_cube
    public :: pelib_ifc_do_infld, pelib_ifc_do_lf
    public :: pelib_ifc_activate, pelib_ifc_deactivate
    public :: pelib_ifc_init, pelib_ifc_finalize, pelib_ifc_input_reader
    public :: pelib_ifc_fock, pelib_ifc_energy, pelib_ifc_response, pelib_ifc_london
    public :: pelib_ifc_molgrad, pelib_ifc_infld, pelib_ifc_lf, pelib_ifc_localfield
    public :: pelib_ifc_mep, pelib_ifc_mep_noqm, pelib_ifc_cube
    public :: pelib_ifc_set_mixed, pelib_ifc_mixed
    public :: pelib_ifc_do_savden, pelib_ifc_do_twoints
    public :: pelib_ifc_save_density, pelib_ifc_twoints
    public :: pelib_ifc_get_num_core_nuclei
#if defined(VAR_MPI)
    public :: pelib_ifc_slave
#endif
    ! TODO: update the following interface routines
    public :: pelib_ifc_grad, pelib_ifc_lin, pelib_ifc_lr, pelib_ifc_qro
    public :: pelib_ifc_cro, pelib_ifc_rspmcqr
    public :: pelib_ifc_pecc
    public :: pelib_ifc_transformer, pelib_ifc_qrtransformer
    ! edh for debug purposes
    public :: pelib_ifc_qrtest !called instead from E3INIT

contains

logical function use_pelib()
    use pe_variables, only: peqm
    if (peqm) then
        use_pelib = .true.
    else
        use_pelib = .false.
    end if
end function use_pelib

logical function pelib_ifc_gspol()
    use pe_variables, only: pe_gspol
    if (pe_gspol) then
        pelib_ifc_gspol = .true.
    else
        pelib_ifc_gspol = .false.
    end if
end function pelib_ifc_gspol

subroutine pelib_ifc_set_mixed(do_mixed)
    use pe_variables, only: mixed
    logical :: do_mixed
    if (do_mixed) then
        mixed = .true.
    else
        mixed = .false.
    end if
end subroutine pelib_ifc_set_mixed

logical function pelib_ifc_do_mep()
    use pe_variables, only: pe_mep
    if (pe_mep) then
        pelib_ifc_do_mep = .true.
    else
        pelib_ifc_do_mep = .false.
    end if
end function pelib_ifc_do_mep

logical function pelib_ifc_do_mep_noqm()
    use pe_variables, only: pe_mep, mep_qmcube
    if (pe_mep .and. .not. mep_qmcube) then
        pelib_ifc_do_mep_noqm = .true.
    else
        pelib_ifc_do_mep_noqm = .false.
    end if
end function pelib_ifc_do_mep_noqm

logical function pelib_ifc_do_cube()
    use pe_variables, only: pe_cube
    if (pe_cube) then
        pelib_ifc_do_cube = .true.
    else
        pelib_ifc_do_cube = .false.
    end if
end function pelib_ifc_do_cube

logical function pelib_ifc_do_infld()
    use pe_variables, only: pe_infld
    if (pe_infld) then
        pelib_ifc_do_infld = .true.
    else
        pelib_ifc_do_infld = .false.
    end if
end function pelib_ifc_do_infld

logical function pelib_ifc_do_lf()
    use pe_variables, only: pe_lf
    if (pe_lf) then
        pelib_ifc_do_lf = .true.
    else
        pelib_ifc_do_lf = .false.
    end if
end function pelib_ifc_do_lf

logical function pelib_ifc_do_savden()
    use pe_variables, only: pe_savden
    if (pe_savden) then
        pelib_ifc_do_savden = .true.
    else
        pelib_ifc_do_savden = .false.
    end if
end function pelib_ifc_do_savden

logical function pelib_ifc_do_twoints()
    use pe_variables, only: pe_twoints
    if (pe_twoints) then
        pelib_ifc_do_twoints = .true.
    else
        pelib_ifc_do_twoints = .false.
    end if
end function pelib_ifc_do_twoints

subroutine pelib_ifc_activate()
    use pe_variables, only: peqm
    call qenter('pelib_ifc_activate')
    if (use_pelib()) call quit('PElib already active')
    peqm = .true.
    call qexit('pelib_ifc_activate')
end subroutine pelib_ifc_activate

subroutine pelib_ifc_deactivate()
    use pe_variables, only: peqm
    call qenter('pelib_ifc_deactivate')
    if (.not. use_pelib()) call quit('PElib already deactivated')
    peqm = .false.
    call qexit('pelib_ifc_deactivate')
end subroutine pelib_ifc_deactivate

subroutine pelib_ifc_input_reader(word)
    use polarizable_embedding, only: pe_input_reader
#include "priunit.h"
    character(len=7), intent(inout) :: word
    call qenter('pelib_ifc_input_reader')
    call pe_input_reader(word, lucmd)
    call qexit('pelib_ifc_input_reader')
end subroutine pelib_ifc_input_reader

subroutine pelib_ifc_init()
    use polarizable_embedding, only: pe_init
#include "priunit.h"
#include "mxcent.h"
#include "nuclei.h"
    call qenter('pelib_ifc_init')
    call pe_init(lupri, cord(1:3,1:natoms), charge(1:natoms))
    call qexit('pelib_ifc_init')
end subroutine pelib_ifc_init

subroutine pelib_ifc_finalize()
    use pe_variables, only: peqm
    use polarizable_embedding, only: pe_finalize
    call qenter('pelib_ifc_finalize')
    if (.not. use_pelib()) call quit('PElib not active')
    call pe_finalize()
    call qexit('pelib_ifc_finalize')
end subroutine pelib_ifc_finalize

subroutine pelib_ifc_fock(denmats, fckmats, energy)
    use polarizable_embedding, only: pe_master
#include "inforb.h"
    real*8, dimension(nnbasx), intent(in) :: denmats
    real*8, dimension(nnbasx), intent(out) :: fckmats
    real*8, intent(out) :: energy
    real*8, dimension(1) :: energies
    call qenter('pelib_ifc_fock')
    if (.not. use_pelib()) call quit('PElib not active')
#if defined(VAR_MPI)
    call pelib_ifc_start_slaves(1)
#endif
    call pe_master(runtype='full_fock', &
                   triang=.true., &
                   ndim=nbast, &
                   nmats=1, &
                   denmats=denmats, &
                   fckmats=fckmats, &
                   expvals=energies)
    energy = energies(1)
    call qexit('pelib_ifc_fock')
end subroutine pelib_ifc_fock

subroutine pelib_ifc_mixed(denmats, fckmats, energy)
    use polarizable_embedding, only: pe_master
#include "inforb.h"
    real*8, dimension(2*nnbasx), intent(in) :: denmats
    real*8, dimension(nnbasx), intent(out) :: fckmats
    real*8, intent(out) :: energy
    real*8, dimension(1) :: energies
    call qenter('pelib_ifc_mixed')
    if (.not. use_pelib()) call quit('PElib not active')
#if defined(VAR_MPI)
    call pelib_ifc_start_slaves(1)
#endif
    call pe_master(runtype='full_fock', &
                   triang=.true., &
                   ndim=nbast, &
                   nmats=1, &
                   denmats=denmats, &
                   fckmats=fckmats, &
                   expvals=energies)
    energy = energies(1)
    call qexit('pelib_ifc_mixed')
end subroutine pelib_ifc_mixed

subroutine pelib_ifc_energy(denmats, energy)
    use polarizable_embedding, only: pe_master
#include "inforb.h"
    real*8, dimension(nnbasx), intent(in) :: denmats
    real*8, intent(out), optional :: energy
    real*8, dimension(1) :: energies
    call qenter('pelib_ifc_energy')
    if (.not. use_pelib()) call quit('PElib not active')
#if defined(VAR_MPI)
    call pelib_ifc_start_slaves(2)
#endif
    call pe_master(runtype='print_energy', &
                   triang=.true., &
                   ndim=nbast, &
                   nmats=1, &
                   denmats=denmats, &
                   expvals=energies)
    if (present(energy)) then
        energy = energies(1)
    end if
    call qexit('pelib_ifc_energy')
end subroutine pelib_ifc_energy

subroutine pelib_ifc_molgrad(denmats, molgrad)
    use polarizable_embedding, only: pe_master
#include "inforb.h"
#include "mxcent.h"
#include "nuclei.h"
    real*8, dimension(nnbasx), intent(in) :: denmats
    real*8, dimension(3*natoms), intent(out) :: molgrad
    call qenter('pelib_ifc_molgrad')
    if (.not. use_pelib()) call quit('PElib not active')
#if defined(VAR_MPI)
    call pelib_ifc_start_slaves(5)
#endif
    call pe_master(runtype='molecular_gradient', &
                   triang=.true., &
                   ndim=nbast, &
                   nmats=1, &
                   denmats=denmats, &
                   expvals=molgrad)
        call qexit('pelib_ifc_molgrad')
end subroutine pelib_ifc_molgrad

subroutine pelib_ifc_infld()
    use polarizable_embedding, only: pe_master
    call qenter('pelib_ifc_infld')
    if (.not. use_pelib()) call quit('PElib not active')
    call pe_master(runtype='infld')
    call qexit('pelib_ifc_infld')
end subroutine pelib_ifc_infld

subroutine pelib_ifc_response(denmats, fckmats, nmats)
    use polarizable_embedding, only: pe_master
#include "inforb.h"
    integer, intent(in) :: nmats
    real*8, dimension(nmats*nnbasx), intent(in) :: denmats
    real*8, dimension(nmats*nnbasx), intent(out) :: fckmats
    call qenter('pelib_ifc_response')
    if (.not. use_pelib()) call quit('PElib not active')
#if defined(VAR_MPI)
    call pelib_ifc_start_slaves(3)
#endif
    call pe_master(runtype='dynamic_response', &
                   triang=.true., &
                   ndim=nbast, &
                   nmats=nmats, &
                   denmats=denmats, &
                   fckmats=fckmats)
    call qexit('pelib_ifc_response')
end subroutine pelib_ifc_response

subroutine pelib_ifc_london(fckmats)
    use polarizable_embedding, only: pe_master
#include "inforb.h"
    real*8, dimension(3*n2basx), intent(out) :: fckmats
    integer :: i, j, k, l, m
    real*8, dimension(:), allocatable :: fckmats_packed
    call qenter('pelib_ifc_london')
    if (.not. use_pelib()) call quit('PElib not active')
    allocate(fckmats_packed(3*nnbasx))
#if defined(VAR_MPI)
    call pelib_ifc_start_slaves(4)
#endif
    call pe_master('magnetic_gradient', &
                   triang=.true., &
                   ndim=nbast, &
                   fckmats=fckmats_packed)
    do i = 1, 3
        j = (i - 1) * nnbasx + 1
        k = i * nnbasx
        l = (i - 1) * n2basx + 1
        m = i * n2basx
        call daptge(nbas, fckmats_packed(j:k), fckmats(l:m))
    end do
    deallocate(fckmats_packed)
    call qexit('pelib_ifc_london')
end subroutine pelib_ifc_london

subroutine pelib_ifc_localfield(eefmats)
    use polarizable_embedding, only: pe_master
#include "inforb.h"
    real*8, dimension(:), intent(out) :: eefmats
    integer :: i, ndim
    call qenter('pelib_ifc_localfield')
    if (.not. use_pelib()) call quit('PElib not active')
#if defined(VAR_MPI)
    call pelib_ifc_start_slaves(8)
#endif
    call pe_master(runtype='effdipole', &
                   triang=.true., &
                   ndim=nbast, &
                   fckmats=eefmats)
    call qexit('pelib_ifc_localfield')
end subroutine pelib_ifc_localfield

subroutine pelib_ifc_lf()
    use polarizable_embedding, only: pe_master
#include "priunit.h"
#include "inforb.h"
#include "inftap.h"
#include "orgcom.h"
    real*8, dimension(nnbasx) :: dip
    real*8, dimension(:),allocatable :: fckmats
    integer :: i, j, k
    logical :: lopen
    character*8 :: lblinf(2)
    allocate(fckmats(3*nnbasx))
    call qenter('pelib_ifc_lf')
    if (.not. use_pelib()) call quit('PElib not active')

    call flshfo(lupri)
    lopen = .false.
    dip = 0.0d0
    fckmats = 0.0d0

#if defined(VAR_MPI)
    call pelib_ifc_start_slaves(8)
#endif
    call flshfo(lupri)
    call pe_master(runtype='effdipole', &
                   triang=.true., &
                   ndim=nbast, &
                   fckmats=fckmats)

    if (luprop <= 0) then
        call gpopen(luprop, 'AOPROPER', 'OLD', ' ', 'UNFORMATTED', 0, .false.)
        lopen = .true.
    end if

!   dipole integrals are stored in triangular form
    rewind(luprop)
    call mollb2('XDIPLEN ',lblinf,luprop,LUPRI)
    call readt(luprop, nnbasx, dip)
    dip(:) = dip(:) + fckmats(1:nnbasx)
    call WRTPRO(dip,nnbasx,'XLFDIPLN',lblinf,0)

    rewind(luprop)
    call mollb2('YDIPLEN ',lblinf,luprop,LUPRI)
    call readt(luprop, nnbasx, dip)
    dip(:) = dip(:) + fckmats(nnbasx+1:2*nnbasx)
    call WRTPRO(dip,nnbasx,'YLFDIPLN',lblinf,0)

    rewind(luprop)
    call mollb2('ZDIPLEN ',lblinf,luprop,LUPRI)
    call readt(luprop, nnbasx, dip)
    dip(:) = dip(:) + fckmats(2*nnbasx+1:3*nnbasx)
    call WRTPRO(dip,nnbasx,'ZLFDIPLN',lblinf,0)

    deallocate(fckmats)
    if (lopen) call gpclose(luprop,'KEEP')
    call qexit('pelib_ifc_lf')

end subroutine pelib_ifc_lf

subroutine pelib_ifc_mep(denmats)
  use polarizable_embedding, only: pe_master
  implicit none
#include "inforb.h"
  real*8, dimension(nnbasx), intent(in) :: denmats
  call qenter('pelib_ifc_mep')
#if defined(VAR_MPI)
  call pelib_ifc_start_slaves(6)
#endif
  call pe_master(runtype='mep', &
                 triang=.true., &
                 ndim=nbast, &
                 nmats=1, &
                 denmats=denmats)
  call qexit('pelib_ifc_mep')
end subroutine pelib_ifc_mep

subroutine pelib_ifc_mep_noqm()
    use polarizable_embedding, only: pe_master
    implicit none
    call qenter('pelib_ifc_mep_noqm')
#if defined(VAR_MPI)
    call pelib_ifc_start_slaves(6)
#endif
    call pe_master(runtype='mep', &
                   triang=.true., &
                   ndim=0, &
                   nmats=0)
    call qexit('pelib_ifc_mep_noqm')
end subroutine pelib_ifc_mep_noqm

subroutine pelib_ifc_cube(denmats, idx)
    use polarizable_embedding, only: pe_master
    implicit none
#include "inforb.h"
    real*8, dimension(nnbasx), intent(in) :: denmats
    integer, intent(in) :: idx
    call qenter('pelib_ifc_cube')
#if defined(VAR_MPI)
    call pelib_ifc_start_slaves(7)
#endif
    call pe_master(runtype='cube', &
                   triang=.true., &
                   ndim=nbast, &
                   nmats=1, &
                   denmats=denmats, &
                   idx=idx)
    call qexit('pelib_ifc_cube')
end subroutine pelib_ifc_cube

subroutine pelib_ifc_save_density(ao_denmat, mo_fckmat, mo_coefficients)
    use pde_utils, only: pde_save_density
#include "iprtyp.h"
#include "maxorb.h"
#include "infpar.h"
#include "inforb.h"
    real*8, dimension(:), intent(in) :: ao_denmat
    real*8, dimension(:), intent(in) :: mo_fckmat
    real*8, dimension(nbast,norbt), intent(in) :: mo_coefficients
    real*8, dimension(:), allocatable :: mo_energies
    real*8, dimension(:,:), allocatable :: ew_denmat
    real*8, dimension(:,:), allocatable :: temp
    integer, parameter :: iprtyp = POLARIZABLE_EMBEDDING
    integer, parameter :: runtyp = 9
    integer :: i
    call qenter('pelib_ifc_save_density')
    allocate(temp(nbast,nisht))
    allocate(mo_energies(nisht))
    do i = 1, nisht
        mo_energies(i) = mo_fckmat(i*(i+1)/2)
        temp(:,i) = mo_energies(i) * mo_coefficients(:,i)
    end do
    allocate(ew_denmat(nbast,nbast))
    ew_denmat = matmul(mo_coefficients(:,1:nisht), transpose(temp))
    deallocate(temp)
#if defined(VAR_MPI)
    if (nodtot >= 1) then
        call mpixbcast(iprtyp, 1, 'INTEGER', master)
        call mpixbcast(runtyp, 1, 'INTEGER', master)
    end if
#endif
    call pde_save_density(ao_denmat, ew_denmat, nbast)
    deallocate(ew_denmat)
    deallocate(mo_energies)
    call qexit('pelib_ifc_save_density')
end subroutine pelib_ifc_save_density

subroutine pelib_ifc_twoints(work, lwork)
    use pde_utils, only: pde_twoints, pde_get_fragment_density
#include "inforb.h"
    real*8, dimension(:), intent(inout) :: work
    integer, intent(in) :: lwork
    real*8, dimension(:), allocatable :: overlap
    real*8, dimension(:), allocatable :: core_fckmat
    real*8, dimension(:), allocatable :: packed_frag_denmat
    real*8, dimension(:,:), allocatable :: full_overlap
    real*8, dimension(:,:), allocatable :: full_fckmat
    real*8, dimension(:,:), allocatable :: full_denmat
    real*8, dimension(:,:), allocatable :: frag_denmat
    integer :: i, j, k
    integer :: core_nbast
    integer :: frag_nbast
    integer, dimension(1) :: isymdm, ifctyp
    call qenter('pelib_ifc_twoints')
    call pde_get_fragment_density(packed_frag_denmat, frag_nbast)
    allocate(frag_denmat(frag_nbast,frag_nbast))
    frag_denmat = 0.0d0
    call dunfld(frag_nbast, packed_frag_denmat, frag_denmat)
    core_nbast = nbast - frag_nbast
    allocate(full_denmat(nbast,nbast))
    full_denmat = 0.0d0
    full_denmat(core_nbast+1:nbast,core_nbast+1:nbast) = frag_denmat
    deallocate(frag_denmat)
    ! IFCTYP = +/-XY
    !   X indicates symmetry about diagonal
    !     X = 0 No symmetry
    !     X = 1 Symmetric
    !     X = 2 Anti-symmetric
    !   Y indicates contributions
    !     Y = 0 No contribution
    !     Y = 1 Coulomb
    !     Y = 2 Exchange
    !     Y = 3 Coulomb + Exchange
    !   + sign: alpha + beta matrix (singlet)
    !   - sign: alpha - beta matrix (triplet)
    ! SIRFCK(fckmat, denmat, ?, isymdm, ifctyp, direct, work, nwrk)
    allocate(full_fckmat(nbast,nbast))
    full_fckmat = 0.0d0
    isymdm = 1
    ifctyp = 11
    call sirfck(full_fckmat, full_denmat, 1, isymdm, ifctyp, .true., work(1), lwork)
    deallocate(full_denmat)
    allocate(core_fckmat(core_nbast*(core_nbast+1)/2))
    core_fckmat = 0.0d0
    call dgetsp(core_nbast, full_fckmat(1:core_nbast,1:core_nbast), core_fckmat)
    deallocate(full_fckmat)
    allocate(overlap(nnbast))
    overlap = 0.0d0
    call rdonel('OVERLAP', .true., overlap, nnbast)
    allocate(full_overlap(nbast,nbast))
    full_overlap = 0.0d0
    call dsptge(nbast, overlap, full_overlap)
    deallocate(overlap)
    call pde_twoints(core_fckmat, full_overlap(1:core_nbast,core_nbast+1:nbast), nbast)
    deallocate(full_overlap, core_fckmat)
    call qexit('pelib_ifc_twoints')
end subroutine pelib_ifc_twoints

integer function pelib_ifc_get_num_core_nuclei()
    use pde_utils, only: pde_get_num_core_nuclei
    integer :: num_nuclei
    num_nuclei = pde_get_num_core_nuclei()
    if (num_nuclei <= 0) then
        call quit('Number of core nuclei must be more than zero')
    end if
    pelib_ifc_get_num_core_nuclei = num_nuclei
end function pelib_ifc_get_num_core_nuclei

#if defined(VAR_MPI)
subroutine pelib_ifc_slave(runtype)
    use polarizable_embedding, only: pe_slave
    use pde_utils, only: pde_save_density
    implicit none
#include "inforb.h"
    integer, intent(in) :: runtype
    real*8, dimension(:), allocatable :: ao_denmat
    real*8, dimension(:,:), allocatable :: dummy
    call qenter('pelib_ifc_slave')
    if (runtype == 1) then
        call pe_slave('full_fock')
    else if (runtype == 2) then
        call pe_slave('print_energy')
    else if (runtype == 3) then
        call pe_slave('dynamic_response')
    else if (runtype == 4) then
        call pe_slave('magnetic_gradient')
    else if (runtype == 5) then
        call pe_slave('molecular_gradient')
    else if (runtype == 6) then
        call pe_slave('mep')
    else if (runtype == 7) then
        call pe_slave('cube')
    else if (runtype == 8) then
        call pe_slave('effdipole')
    else if (runtype == 9) then
        allocate(ao_denmat(nnbasx))
        allocate(dummy(1,1))
        call pde_save_density(ao_denmat, dummy, nbast)
        deallocate(ao_denmat, dummy)
    end if
    call qexit('pelib_ifc_slave')
end subroutine pelib_ifc_slave
#endif

subroutine pelib_ifc_grad(cref, cmo, cindx, dv, grd, energy, wrk, nwrk)
!     Written by Erik Donovan Hedegård (edh) and Jogvan Magnus H. Olsen
!                based on PCMGRAD
!
!     Purpose:  calculate (MCSCF) energy and gradient contribution
!               from an embedding potential using the PE library
!
!     Output:
!     grd       MCSCF gradient with PE contribution added
!     energy    total PE energy
!
! Used from common blocks:
!   INFVAR: NCONF,  NWOPT,  NVAR,   NVARH
!   INFORB: NNASHX, NNBASX, NNORBX, etc.
!   INFTAP: LUIT2
    implicit none
#include "priunit.h"
#include "infvar.h"
#include "inforb.h"
#include "inftap.h"

    integer :: nwrk
    real*8 :: energy
    real*8, dimension(*) :: cref, cmo, cindx, dv, grd
    real*8, dimension(nwrk) :: wrk
    character*8 :: star8 = '********'
    character*8 :: solvdi = 'SOLVDIAG'
    character*8 :: eodata = 'EODATA  '

    logical :: fndlab
    integer :: nc4, nw4, i
    real*8 :: solelm, ddot
    real*8 :: tmo, tac, test
    real*8, dimension(1) :: etmp
    real*8, dimension(:), allocatable :: fckmo, fckac
    real*8, dimension(:), allocatable :: pegrd, diape
    real*8, dimension(:), allocatable :: dcao, dvao, fdtao, fckao

    call qenter('pelib_ifc_grad')
    if (.not. use_pelib()) call quit('PElib not active')

    allocate(dcao(n2basx), dvao(n2basx))
    call fckden((nisht > 0), (nasht > 0), dcao, dvao, cmo, dv, wrk, nwrk)
    if (nisht == 0) dcao = 0.0d0
    if (nasht > 0) dcao = dcao + dvao
    deallocate(dvao)
    allocate(fdtao(nnbasx))
    call dgefsp(nbast, dcao, fdtao)
    deallocate(dcao)
    allocate(fckao(nnbasx))
    call pelib_ifc_fock(fdtao, fckao, energy)
    deallocate(fdtao)
    allocate(fckmo(nnorbx))
    call uthu(fckao, fckmo, cmo, wrk, nbast, norbt)
    deallocate(fckao)

    allocate(fckac(nnashx))
    fckac = 0.0d0
    if (nasht > 0) call getac2(fckmo, fckac)
    tmo = solelm(dv, fckac, fckmo, tac)

    allocate(pegrd(nvarh))
    pegrd = 0.0d0
    if (nconf > 1) then
        ! edh: SOLGC calc. < u | Fg | 0 > + < 0 | Fg | 0 > c_u
        call solgc(cref, fckac, tac, pegrd, cindx, wrk, nwrk)
    end if
    if (nwopt > 0) then
        ! edh: SOLGO calc. 2 < 0 | [Ers, Fg] | 0 >
        call solgo(2.0d0, dv, fckmo, pegrd(1+nconf:nvarh))
    end if

    allocate(diape(nvar))
    diape = 0.0d0
    call soldia(tac, fckac, cindx, fckmo, dv, diape, wrk, nwrk)
    diape = - diape
    deallocate(fckmo, fckac)

    !--------------- Orthogonality test ----------------
    test = ddot(nconf, cref, 1, pegrd, 1)
    if (abs(test) > 1.0d-8) then
        nwarn = nwarn + 1
        write(lupri,*) ' --- PE GRADIENT WARNING --- '
        write(lupri,*) ' < CREF | GRAD > =', test
    end if

    ! Add PE gradient contribution to MCSCF gradient
    call daxpy(nvarh, 1.0d0, pegrd, 1, grd, 1)
    deallocate(pegrd)
    if (luit2 > 0) then
        nc4 = max(nconf, 4)
        nw4 = max(nwopt, 4)
        rewind luit2
        if (fndlab(eodata,luit2)) backspace luit2
        write(luit2) star8, star8, star8, solvdi
        if (nconf > 1) call writt(luit2, nc4, diape)
        write(luit2) star8, star8, star8, eodata
    end if

    call qexit('pelib_ifc_grad')

end subroutine pelib_ifc_grad

subroutine pelib_ifc_lin(ncsim, nosim, bcvecs, bovecs, cref, cmo, cindx, dv, dtv,&
                         scvecs, sovecs, orblin, wrk, nwrk)
!
! Written by Erik Donovan Hedegård and Jogvan Magnus H. Olsen
!            after original code by  Hans Joergen Aa. Jensen
!
! Common driver for pe_lnc and pe_lno
!
!   Used from common blocks:
!   INFLIN : NWOPPT,NVARPT
    implicit none
#include "priunit.h"
#include "inflin.h"
#include "infvar.h"
#include "inforb.h"

    logical :: orblin
    integer :: ncsim, nosim, nwrk
    integer, dimension(*) :: cindx
    real*8, dimension(*) :: bcvecs, bovecs, scvecs, sovecs
    real*8, dimension(*) :: cmo, cref, dv, dtv
    real*8, dimension(nwrk) :: wrk

    integer :: nso

    call qenter('pelib_ifc_lin')
    if (.not. use_pelib()) call quit('PElib not active')

    if (ncsim > 0) then
        call pe_lnc(ncsim, bcvecs, cref, cmo, cindx, dv, dtv, scvecs, wrk, nwrk)
    end if

    if (nosim > 0) then
        if (.not. orblin) then
            nso = nvarpt
        else
            nso = nwoppt
        end if
        call pe_lno(nosim, bovecs, cref, cmo, cindx, dv, sovecs, nso, wrk, nwrk)
    end if

    call qexit('pelib_ifc_lin')

end subroutine pelib_ifc_lin

subroutine pe_lnc(ncsim, bcvecs, cref, cmo, cindx, dv, dtv, scvecs, wrk, nwrk)
!
!  Written by Erik Donovan Hedegaard and Jogvan Magnus H. Olsen
!             after original routine by Hans Jørgen Aa. Jensen
!
!  Purpose:  Calculate Hessian contribution from a polarizable
!            embedding potantial to a csf trial vector.
!
!  Used from common blocks:
!    INFORB : NNASHX, NNORBX, NNBASX, etc.
!    INFVAR : NWOPH
!    INFLIN : NCONST, NVARPT, NWOPPT
!
    use pe_variables, only: pe_polar
    implicit none
#include "priunit.h"
#include "dummy.h"
#include "inforb.h"
#include "infvar.h"
#include "inflin.h"
#include "infdim.h"

    integer :: ncsim, nwrk
    integer, dimension(*) :: cindx
    real*8, dimension(*) :: bcvecs, cref, cmo, dv
    real*8, dimension(nnashx,*) :: dtv
    real*8, dimension(nvarpt,*) :: scvecs
    real*8, dimension(nwrk) :: wrk

    logical :: fndlab
    integer :: i, j, jscvec, mwoph
    real*8 :: tfxc, tfyc, tfycac, solelm, energy
    real*8, dimension(:), allocatable :: dcao, dvao, fdtao, fycao
    real*8, dimension(:), allocatable :: dtvao, fdtvaos, fxcaos
    real*8, dimension(:), allocatable :: tfxcacs, fyc, fycac
    real*8, dimension(:,:), allocatable :: fxcs, fxcacs

    call qenter('pe_lnc')
    if (.not. use_pelib()) call quit('PElib not active')

    allocate(fxcs(nnorbx,ncsim))
    allocate(fxcacs(nnashx,ncsim))
    allocate(tfxcacs(ncsim))
    fxcs = 0.0d0
    fxcacs = 0.0d0
    tfxcacs = 0.0d0

    if (pe_polar) then
        ! Fxc = -R<0|Fe|B>Fe in fxcaos
        allocate(dtvao(n2basx))
        allocate(fdtvaos(ncsim*nnbasx))
        do i = 1, ncsim
            j = (i - 1) * nnbasx + 1
            call fckden(.false., .true., dummy, dtvao, cmo, dtv(:,i), wrk, nwrk)
            call dgefsp(nbast, dtvao, fdtvaos(j))
        end do
        deallocate(dtvao)
        allocate(fxcaos(ncsim*nnbasx))
        call pelib_ifc_response(fdtvaos, fxcaos, ncsim)
        deallocate(fdtvaos)
        do i = 1, ncsim
            j = (i - 1) * nnbasx + 1
            call uthu(fxcaos(j), fxcs(:,i), cmo, wrk, nbast, norbt)
            if (nasht > 0) call getac2(fxcs(:,i), fxcacs(:,i))
            tfxc = solelm(dv, fxcacs(:,i), fxcs(:,i), tfxcacs(i))
        end do
        deallocate(fxcaos)
    end if

    ! Fg = Vmul -R<0|Fe|0>Fe in fyc
    allocate(dcao(n2basx), dvao(n2basx))
    call fckden((nisht > 0), (nasht > 0), dcao, dvao, cmo, dv, wrk, nwrk)
    if (nisht == 0) dcao = 0.0d0
    if (nasht > 0) dcao = dcao + dvao
    deallocate(dvao)
    allocate(fdtao(nnbasx))
    call dgefsp(nbast, dcao, fdtao)
    deallocate(dcao)
    allocate(fycao(nnbasx))
    call pelib_ifc_fock(fdtao, fycao, energy)
    deallocate(fdtao)
    allocate(fyc(nnorbx))
    call uthu(fycao, fyc, cmo, wrk, nbast, norbt)
    deallocate(fycao)
    allocate(fycac(nnashx))
    if (nasht > 0) call getac2(fyc, fycac)
    tfyc = solelm(dv, fycac, fyc, tfycac)

    ! ...CSF part of sigma vectors
    call solsc(ncsim, 0, bcvecs, cref, scvecs, fxcacs, fycac, tfxcacs, tfycac,&
               cindx, wrk, nwrk)

    if (nwoppt > 0) then
        mwoph = nwoph
        nwoph = nwoppt
        jscvec = 1 + nconst
        do i = 1, ncsim
            if (pe_polar) then
                call solgo(2.0d0, dv, fxcs(:,i), scvecs(jscvec,i))
            end if
            call solgo(0.0d0, dtv(:,i), fyc, scvecs(jscvec,i))
        end do
        nwoph = mwoph
    end if

    deallocate(fxcacs, fycac, tfxcacs)
    deallocate(fyc, fxcs)

    call qexit('pe_lnc')

end subroutine pe_lnc

subroutine pe_lno(nosim, bovecs, cref, cmo, cindx, dv, sovecs, nso,&
                  wrk, nwrk)
!
!  Written by Erik Donovan Hedegaard and Jogvan Magnus H. Olsen
!             after original code by Hans Jorgen Aa. Jensen
!
!  Purpose:  Calculate Hessian contribution from a
!            PE potential to an orbital trial vector.
!
!  NSVEC     may be NVAR or NWOPT, dependent on LINTRN
!
!  Used from common blocks:
!    INFORB : NNASHX, NNORBX, NNBASX, etc.
!    INFVAR : JWOP
!    INFLIN : NWOPPT, NVARPT, NCONST, NCONRF
!
    use pe_variables, only: pe_polar
    implicit none
#include "priunit.h"
#include "dummy.h"
#include "inforb.h"
#include "infvar.h"
#include "inflin.h"

    integer :: nosim, nso, nwrk
    integer, dimension(*) :: cindx
    real*8, dimension(*) :: cref, cmo, dv
    real*8, dimension(nwrk) :: wrk
    real*8, dimension(nwoppt,*) :: bovecs
    real*8, dimension(nso,*) :: sovecs

    integer :: i, j, jsovec, mwoph, ncolim
    logical :: fulhes, fndlab
    real*8 :: solelm
    real*8 :: txyo
    real*8 :: energy
    real*8, dimension(:), allocatable :: txyoacs
    real*8, dimension(:), allocatable :: ubodcao, ubodvao
    real*8, dimension(:), allocatable :: bodtaos, fxoaos
    real*8, dimension(:), allocatable :: fckmo, fyo, ufyo
    real*8, dimension(:), allocatable :: dcao, dvao, fdtao, fckao
    real*8, dimension(:,:), allocatable :: ubovecs, fxos
    real*8, dimension(:,:), allocatable :: fxyos, fxyoacs


    call qenter('pe_lno')
    if (.not. use_pelib()) call quit('PElib not active')

    allocate(ubovecs(n2orbx,nosim))
    if (nosim > 0) then
        do i = 1, nosim
            call upkwop(nwoppt, jwop, bovecs(:,i), ubovecs(:,i))
        end do
    end if

    ! 1. Calculation of Fxo = R*<0|Fe(k)|O>Fe
    !    Store in fxos
    if (pe_polar) then
        allocate(ubodcao(n2basx), ubodvao(n2basx))
        allocate(bodtaos(nosim*nnbasx))
        do i = 1, nosim
            j = (i - 1) * nnbasx + 1
            call tr1den(cmo, ubovecs(:,i), dv, ubodcao, ubodvao, wrk, nwrk)
            if (nasht > 0) ubodcao = ubodcao + ubodvao
            call dgefsp(nbast, ubodcao, bodtaos(j))
        end do
        deallocate(ubodcao, ubodvao)
        allocate(fxoaos(nosim*nnbasx))
        call pelib_ifc_response(bodtaos, fxoaos, nosim)
        deallocate(bodtaos)
        allocate(fxos(nnorbx,nosim))
        do i = 1, nosim
            j = (i - 1) * nnbasx + 1
            call uthu(fxoaos(j), fxos(:,i), cmo, wrk, nbast, norbt)
        end do
        deallocate(fxoaos)
    end if

    ! 2. Calculation of Fyo = V(k) + R<0|F|0>Fe(k)
    !    Store in fyos
    allocate(dcao(n2basx), dvao(n2basx))
    call fckden((nisht > 0), (nasht > 0), dcao, dvao, cmo, dv, wrk, nwrk)
    if (nisht == 0) dcao = 0.0d0
    if (nasht > 0) dcao = dcao + dvao
    deallocate(dvao)
    allocate(fdtao(nnbasx))
    call dgefsp(nbast, dcao, fdtao)
    deallocate(dcao)
    allocate(fckao(nnbasx))
    call pelib_ifc_fock(fdtao, fckao, energy)
    deallocate(fdtao)
    allocate(fckmo(nnorbx))
    call uthu(fckao, fckmo, cmo, wrk, nbast, norbt)
    deallocate(fckao)
    allocate(fyo(n2orbx))
    call dsptsi(norbt, fckmo, fyo)
    deallocate(fckmo)

    allocate(ufyo(n2orbx), txyoacs(nosim))
    allocate(fxyos(nnorbx,nosim), fxyoacs(nnashx,nosim))
    do i = 1, nosim
        ufyo = 0.0d0
        call tr1uh1(ubovecs(:,i), fyo, ufyo, 1)
        call dgetsp(norbt, ufyo, fxyos(:,i))
        if (pe_polar) then
            call daxpy(nnorbx, 1.0d0, fxos(:,i), 1, fxyos(:,i), 1)
        end if
        if (nasht > 0) then
            call getac2(fxyos(:,i), fxyoacs(:,i))
        end if
        txyo = solelm(dv, fxyoacs(:,i), fxyos(:,i), txyoacs(i))
    end do
    ! 3.   /        <0[Epq,Fxo + Fyo]|0>      \  orbital part
    !      \ 2<0|Fyo + Fxo|mu> - <0|Fyo|0>*c0 /  CSF part
    !     ... CSF part of sigma vectors
    if (lsymrf == lsymst) then
        ncolim = 1
    else
        ncolim = 0
    end if

    ! Determine if full Hessian or only orbital Hessian
    fulhes = (nso == nvarpt)
    if (fulhes) then
        jsovec = 1 + nconst
    else
        jsovec = 1
    end if

    if (fulhes .and. (nconst > ncolim)) then
        call solsc(0, nosim, dummy, cref, sovecs, fxyoacs, dummy, txyoacs,&
                   dummy, cindx, wrk, nwrk)
    end if

    ! ... orbital part of sigma vectors
    mwoph = nwoph
    nwoph = nwoppt
    ! ... tell SOLGO only to use the NWOPPT first JWOP entries
    do i = 1, nosim
        call solgo(2.0d0, dv, fxyos(:,i), sovecs(jsovec,i))
    end do
    nwoph = mwoph

    call qexit('pe_lno')

end subroutine pe_lno

subroutine pelib_ifc_lr(ncsim, nosim, bcvecs, bovecs, cref, cmo, cindx, udv,&
                        dv, udvtr, dvtr, dtv, dtvtr, scvecs, sovecs, wrk, nwrk)
    implicit none
#include "priunit.h"
#include "dummy.h"
#include "wrkrsp.h"
#include "infrsp.h"
#include "inftap.h"

    integer :: ncsim, nosim, nwrk
    integer, dimension(*) :: cindx
    real*8, dimension(*) :: bcvecs, bovecs
    real*8, dimension(*) :: cref, cmo, udv, dv
    real*8, dimension(*) :: udvtr, dvtr, dtv, dtvtr
    real*8, dimension(*) :: scvecs, sovecs
    real*8, dimension(nwrk) :: wrk

    call qenter('pelib_ifc_lr')
    if (.not. use_pelib()) call quit('PElib not active')

    if (ncsim > 0 .and. .not. soppa) then
        call pe_rsplnc(ncsim, bcvecs, cref, cmo, cindx, udv, dv,&
                       udvtr, dvtr, dtv, dtvtr, scvecs, wrk, nwrk)
    end if

    if (nosim > 0) then
        call pe_rsplno(nosim, bovecs, cref, cmo, cindx,&
                       udv, dv, udvtr, dvtr, sovecs, wrk, nwrk)
    end if

    call qexit('pelib_ifc_lr')

end subroutine pelib_ifc_lr

subroutine pe_rsplnc(ncsim, bcvecs, cref, cmo, cindx, udv, dv,&
                     udvtr, dvtr, dtv, dtvtr, scvecs, wrk, nwrk)
    use pe_variables, only: pe_polar
    implicit none
#include "priunit.h"
#include "dummy.h"
#include "infrsp.h"
#include "inftap.h"
#include "wrkrsp.h"
#include "inforb.h"
#include "qrinf.h"
#include "infvar.h"

    integer :: i, j
    integer :: ncsim, nwrk
    integer, dimension(*) :: cindx

    real*8, dimension(*) :: bcvecs, cref, cmo, udv, dv
    real*8, dimension(*) :: udvtr, dvtr
    real*8, dimension(n2ashx,*) :: dtv
    real*8, dimension(n2ashx,*) :: dtvtr
    real*8, dimension(kzyvar,*) :: scvecs
    real*8, dimension(nwrk) :: wrk

    real*8 :: ovlap, solelm, tfpeac, tfpe
    real*8, dimension(:,:), allocatable :: udtv
    real*8, dimension(:,:), allocatable :: fxcs, fxcacs
    real*8, dimension(:), allocatable :: dtvao, fuxcs
    real*8, dimension(:), allocatable :: fdtvaos, fxcaos
    real*8, dimension(:), allocatable :: tfxc, tfxcacs
    real*8, dimension(:), allocatable :: fpe, fupe, fpeac

    logical :: lexist, lopen, locdeb
    logical :: fndlab
    logical :: tdm, norho2

    call qenter('pe_rsplnc')
    if (.not. use_pelib()) call quit('PElib not active')

    locdeb = .false.

    lopen = .false.
    tdm = .true.
    norho2 = .true.

    allocate(fxcs(nnorbx,ncsim))
    allocate(fuxcs(n2orbx))
    allocate(fxcacs(nnashx,ncsim))
    allocate(tfxc(ncsim))
    allocate(tfxcacs(ncsim))
    fxcs    = 0.0d0
    fuxcs   = 0.0d0
    fxcacs  = 0.0d0
    tfxc    = 0.0d0
    tfxcacs = 0.0d0

    ! Fxc = R*(<0(L)|Fe|0> + <0|Fe|0(R)>)Fe
    if (pe_polar .or. .not. trplet) then
        call getref(cref, ncref)
        ! ...Construct <0(L)|...|0> + <0|...|0(R)>
        allocate(udtv(n2ashx,ncsim))
        udtv = 0.0d0
        call rsptdm(ncsim, irefsy, ksymst, ncref, kzconf, cref, bcvecs,&
                    udtv, dummy, 0, 0, tdm, norho2, cindx, wrk, 1, nwrk)
        udtv = - 1.0d0 * udtv

        if ( ncsim > 0 ) then
            allocate(fdtvaos(nnbasx*ncsim))
            fdtvaos = 0.0d0
            allocate(dtvao(n2basx))
            dtvao = 0.0d0
            do i = 1, ncsim
               j = (i - 1) * nnbasx + 1
                call fckden2(.false., .true., dummy, dtvao, cmo,&
                             udtv(:,i), wrk, nwrk)
                call dgefsp(nbast, dtvao, fdtvaos(j))
            end do
            deallocate(udtv, dtvao)
        end if

        allocate(fxcaos(ncsim*nnbasx))
        fxcaos = 0.0d0
        call pelib_ifc_response(fdtvaos, fxcaos, ncsim)
        deallocate(fdtvaos)

        do i = 1, ncsim
            j = (i - 1) * nnbasx + 1
            call uthu(fxcaos(j), fxcs(:,i), cmo, wrk, nbast, norbt)
            if (nasht > 0) call getac2(fxcs(:,i), fxcacs(:,i))
            tfxc = solelm(dv, fxcacs(:,i), fxcs(:,i), tfxcacs(i))
        end do
        deallocate(fxcaos)
    end if

    ! Fg = V - <0|F|0>Fe -unpack into fupe
    if (.not. tdhf) then
        allocate(fpe(nnorbx))
        fpe = 0.0d0
        if (lusifc <= 0) then
            call gpopen(lusifc, 'SIRIFC', 'OLD', ' ', 'UNFORMATTED',&
                        idummy, .false.)
            lopen = .true.
        end if
        rewind(lusifc)
        call mollab('PEFMAT  ', lusifc, lupri)
        call readt(lusifc, nnorbx, fpe)
        if (lopen) call gpclose(lusifc, 'KEEP')
        allocate(fupe(n2orbx), fpeac(nnashx))
        fupe = 0.0d0
        fpeac = 0.0d0
        call dsptsi(norbt, fpe, fupe)
        if (nasht > 0) call getac2(fpe, fpeac)
            tfpe = solelm(dv, fpeac, fpe, tfpeac)
        deallocate(fpe)
    end if

    ! Calculate Fxc(Rxc) and Fg(Ryc) contributions to SCVECS(NVAR,NCSIM)
    !  ... CSF part of sigma vectors
    if (locdeb) then
       write(lupri,*)' Linear transformed configuration vector'
       write(lupri,*)' **** Before slvsc in pe_rsplnc **** '
       call output(scvecs,1,kzyvar,1,ncsim,kzyvar,ncsim,1,lupri)
    endif

    call slvsc(ncsim, 0, nnashx, bcvecs, cref, scvecs, fxcacs,&
               fpeac, tfxcacs, tfpeac, cindx, wrk, nwrk)
    deallocate(fxcacs, tfxcacs, fpeac)

    if (locdeb) then
        write(lupri,*)' Linear transformed configuration vector'
        write(lupri,*)' **** After slvsc in pe_rsplnc **** '
        call output(scvecs,1,kzyvar,1,ncsim,kzyvar,ncsim,1,lupri)
    end if

    ! edh: The triplet will only work for response from a single reference
    ! state, where the term from Fxc is 0. Can be generalized (at least to
    ! dublet) by including Fxc...
    ! ... orbital part of sigma vector(s)
    if (kzwopt .gt. 0) then
        do i = 1,ncsim
            fuxcs = 0.0d0
            call dsptsi(norbt,fxcs(:,i), fuxcs)
            if (trplet) then
               ! zero for singlet reference state
            else
               call slvsor(.true.,.true., 1, udv, scvecs(1,i), fuxcs)
            end if
            if (locdeb) then
                write(lupri,*) ' *** After slvsor in pe_rsplnc *** '
                write(lupri,*) 'Orbital part of lin transf conf vec no ', i
                write(lupri,*) ' Txc contribution'
                call output(scvecs(1,i), 1, kzyvar, 1, 1, kzyvar, 1, 1,&
                            lupri)
            end if
            if (trplet) then
               call slvsor(.false.,.false.,1, dtvtr(1,i),scvecs(1,i),fupe)
            else
               call slvsor(.false., .false., 1, dtv(1,i), scvecs(1,i), fupe)
            end if
            if (locdeb) then
                write(lupri,*) 'Orbital part of lin transf conf vec no ', i
                write(lupri,*)' Tg contribution'
                call output(scvecs(1,i), 1, kzyvar, 1, 1, kzyvar, 1, 1,&
                            lupri)
            end if
        end do
        deallocate(fupe, fuxcs)

        if (locdeb) then
            write(lupri,*)' linear transformed conf. vector'
            write(lupri,*)' *** after slvsor in pe_rsplnc *** '
            call output(scvecs, 1, kzyvar, 1, ncsim, kzyvar, ncsim, 1,&
                        lupri)
        end if
    end if

    if (ncref /= kzconf) then
        call quit('ERROR in pe_rsplnc: ncref /= kzconf')
    end if

    call qexit('pe_rsplnc')

end subroutine pe_rsplnc

subroutine pe_rsplno(nosim, bovecs, cref, cmo, cindx, udv, dv,&
                     udvtr, dvtr, sovecs, wrk, nwrk)
    use pe_variables, only: pe_polar, pe_gspol
    implicit none
#include "priunit.h"
#include "dummy.h"
#include "wrkrsp.h"
#include "inforb.h"
#include "infrsp.h"
#include "inftap.h"

    integer :: nosim, nwrk
    integer, dimension(*) :: cindx
    real*8, dimension(*) :: bovecs
    real*8, dimension(kzyvar,*) :: sovecs
    real*8, dimension(*) :: cref, cmo, udv, dv, udvtr, dvtr
    real*8, dimension(nwrk) :: wrk

    integer :: i, j
    real*8 :: txyo
    real*8 :: ddot, slvqlm
    real*8, dimension(:), allocatable :: dcao, dvao
    real*8, dimension(:), allocatable :: daos, fckaos
    real*8, dimension(:), allocatable :: daotrs
    real*8, dimension(:), allocatable :: evec
    real*8, dimension(:,:), allocatable :: ubovecs, evecs, eacs
    real*8, dimension(:), allocatable :: fpemo,fupemo
    real*8, dimension(:), allocatable :: txyoacs
    real*8, dimension(:), allocatable :: ovlp
    logical :: lexist, lopen

    call qenter('pe_rsplno')
    if (.not. use_pelib()) call quit('PElib not active')

    ! return if no polarization and not MCSCF
    if (tdhf .and. .not. pe_polar) then
        call qexit('pe_rsplno')
        return
    ! no polarization for triplet excitations in closed shell SCF
    else if ((nasht == 0) .and. trplet) then
        !write(lupri,*)'WARNING: Triplet PE-response experimental'
        call qexit('pe_rsplno')
        return
    ! ground state polarization approximation
    else if (pe_gspol) then
        call qexit('pe_rsplno')
        return
    ! triplet response for open shell (and MCSCF) not ready yet
    else if (tdhf .and. (nasht > 0) .and. trplet) then
        !write(lupri,*)'WARNING: Triplet code experimental'
        call quit('ERROR: triplet operators for open shell'//&
                  ' systems not implemented')
    end if

    lopen = .false.

    if (.not. tdhf) then
        ! Read Fg = V - <0|F|0>Fe from file
        allocate(fpemo(nnorbx))
        if (lusifc <= 0) then
            call gpopen(lusifc, 'SIRIFC', 'OLD', ' ', 'UNFORMATTED',&
                        idummy, .false.)
            lopen = .true.
        end if
        rewind(lusifc)
        call mollab('PEFMAT  ', lusifc, lupri)
        call readt(lusifc, nnorbx, fpemo)
        if (lopen) call gpclose(lusifc, 'KEEP')
        allocate(fupemo(n2orbx))
        call dsptsi(norbt, fpemo, fupemo)
        deallocate(fpemo)
    end if

    allocate(ubovecs(n2orbx,nosim))
    call rspzym(nosim, bovecs, ubovecs)

    ubovecs = - ubovecs

    if (.not. trplet) then
       allocate(dcao(n2basx), dvao(n2basx), daos(nosim*nnbasx))
       ! Calculate Fxo = <0|Fe(k)|0>Fe
       do i = 1, nosim
           j = (i - 1) * nnbasx + 1
           call deq27(cmo, ubovecs(:,i), udv, dcao, dvao, wrk, nwrk)
           if (nasht > 0) then
               dcao = dcao + dvao
           end if
           call dgefsp(nbast, dcao, daos(j))
       end do
       deallocate(dcao, dvao)

       allocate(fckaos(nosim*nnbasx))
       call pelib_ifc_response(daos, fckaos, nosim)
       deallocate(daos)
    end if

    allocate(evec(nnorbx))
    allocate(evecs(n2orbx,nosim))
    evecs = 0.0d0
    if (.not. tdhf) then
        allocate(eacs(n2ashx,nosim))
        allocate(txyoacs(nosim))
        eacs = 0.0d0
        txyoacs = 0.0d0
    end if
    do i = 1, nosim
        j = (i - 1) * nnbasx + 1
        if (.not. trplet) then
           call uthu(fckaos(j), evec, cmo, wrk, nbast, norbt)
           call dsptsi(norbt, evec, evecs(:,i))
        end if
        ! Fyo = V(k) - <0|F|0>Fe(k)
        if (.not. tdhf) then
            call onexh1(ubovecs(:,i), fupemo, evecs(:,i))
            call getacq(evecs(:,i), eacs(:,i))
            if (trplet) then
              ! do nothing (txyoacs should be zero for triplet with singlet cref)
            else
              txyo = slvqlm(udv, eacs(:,i), evecs(:,i), txyoacs(i))
            end if
        end if
    end do

    deallocate(evec)
    if (.not. tdhf) then
        deallocate(fupemo)
        ! Special triplet handling is taken care of inside slvsc!
        call slvsc(0, nosim, n2ashx, dummy, cref, sovecs, eacs,&
                   dummy, txyoacs, dummy, cindx, wrk, nwrk)
        deallocate(eacs)
        deallocate(txyoacs)
    end if
    ! Note: triplet and singlet calls to slvsor get identical.
    ! Singlet case: singlet evecs and singlet orbital operator.
    ! Triplet case: triplet evecs and triplet orbital operator.
    ! The two cases give the same expressions for gradient elements.
    call slvsor(.true., .true., nosim, udv, sovecs, evecs)
    deallocate(evecs)

    call qexit('pe_rsplno')

end subroutine pe_rsplno

subroutine pelib_ifc_qro(vecb, vecc, etrs, xindx, zymb, zymc, udv, wrk, nwrk,&
                         kzyva, kzyvb, kzyvc, isyma, isymb, isymc, cmo, mjwop)
    implicit none
#include "inforb.h"
#include "infvar.h"
#include "infrsp.h"
#include "infdim.h"
#include "qrinf.h"

    integer :: kzyva, kzyvb, kzyvc
    integer :: isyma, isymb, isymc
    integer :: nwrk
    real*8, dimension(nwrk) :: wrk
    real*8, dimension(kzyva) :: etrs
    real*8, dimension(kzyvb) :: vecb
    real*8, dimension(kzyvc) :: vecc
    real*8, dimension(ncmot) :: cmo
    real*8, dimension(norbt,norbt) :: zymb, zymc
    real*8, dimension(nashdi,nashdi) :: udv
    real*8, dimension(lcindx) :: xindx
    integer, dimension(2,maxwop,8) :: mjwop

    integer :: i, j, k
    integer :: idum = 1
    real*8, dimension(:), allocatable :: udcao, ufcmo
    real*8, dimension(:), allocatable :: dcaos, fcaos
    real*8, dimension(:), allocatable :: fcmo

    call qenter('pelib_ifc_qro')
    if (.not. use_pelib()) call quit('PElib not active')

    call gtzymt(1, vecb, kzyvb, isymb, zymb, mjwop)
    call gtzymt(1, vecc, kzyvc, isymc, zymc, mjwop)

    allocate(udcao(n2basx))
    allocate(ufcmo(n2orbx))
    allocate(dcaos(4*nnbasx))
    dcaos = 0.0d0

    !  D(1k)
    udcao = 0.0d0
    call cdens1(isymb, cmo, zymb, udcao, wrk, nwrk)
    call dgefsp(nbast, udcao, dcaos(1:nnbasx))

    ! D(1k,2k)
    udcao = 0.0d0
    call cdens2(isymb, isymc, cmo, zymb, zymc, udcao,&
                wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
    call dgefsp(nbast, udcao, dcaos(nnbasx+1:2*nnbasx))

    !  D(2k)
    udcao = 0.0d0
    call cdens1(isymc, cmo, zymc, udcao, wrk, nwrk)
    call dgefsp(nbast, udcao, dcaos(2*nnbasx+1:3*nnbasx))

    !  D(2k,1k)
    udcao = 0.0d0
    call cdens2(isymc, isymb, cmo, zymc, zymb, udcao,&
                wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
    call dgefsp(nbast, udcao, dcaos(3*nnbasx+1:4*nnbasx))

    deallocate(udcao)

    allocate(fcaos(4*nnbasx))
    call pelib_ifc_response(dcaos, fcaos, 4)
    deallocate(dcaos)

    allocate(fcmo(nnorbx))
    ufcmo = 0.0d0

    i = 1
    j = nnbasx
    call uthu(2.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
    wrk(1:n2orbx) = 0.0d0
    call dsptsi(norbt, fcmo, wrk(1:n2orbx))
    call oith1(isymc, zymc, wrk(1:n2orbx), ufcmo, isyma)

    i = i + nnbasx
    j = j + nnbasx
    call uthu(fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
    wrk(1:n2orbx) = 0.0d0
    call dsptsi(norbt, fcmo, wrk(1:n2orbx))
    ufcmo = ufcmo + wrk(1:n2orbx)

    i = i + nnbasx
    j = j + nnbasx
    call uthu(2.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
    wrk(1:n2orbx) = 0.0d0
    call dsptsi(norbt, fcmo, wrk(1:n2orbx))
    call oith1(isymb, zymb, wrk(1:n2orbx), ufcmo, isyma)

    i = i + nnbasx
    j = j + nnbasx
    call uthu(fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
    wrk(1:n2orbx) = 0.0d0
    call dsptsi(norbt, fcmo, wrk(1:n2orbx))
    ufcmo = ufcmo + wrk(1:n2orbx)

    call rsp1gr(1, kzyva, idum, 0, isyma, 0, 1, etrs,&
                wrk, idum, idum, 1.0d0, 1, udv, ufcmo, xindx,&
                mjwop, wrk, nwrk, .true., .false., .false.)

    deallocate(fcaos, fcmo, ufcmo)

    call qexit('pelib_ifc_qro')

end subroutine pelib_ifc_qro

subroutine pelib_ifc_qrtest(vecb, vecc, veca, atest, etrs, xindx, zymb, zymc,&
                            den1, udv, wrk, lfree, kzyva, kzyvb, kzyvc,&
                            isyma, isymb, isymc, cmo, mjwop)

      use pe_variables, only: pe_polar
      use polarizable_embedding, only: pe_master

      implicit none

#include "inforb.h"
#include "infvar.h"
#include "infdim.h"
#include "qrinf.h"
#include "priunit.h"
#include "dummy.h"
#include "inftap.h"
#include "infrsp.h"
#include "wrkrsp.h"

      integer :: kzyva, kzyvb, kzyvc
      integer :: isyma, isymb, isymc, isymbc
      integer :: lfree
      integer :: ilsym, irsym, ncl, ncr, kzvarl, kzvarr
      integer :: isymdn, isymst
      integer :: kcref, nzyvec, nzcvec
      integer :: iprone, nzconf, nzvar
      integer :: n2ash

      real*8 :: ovlap
      real*8 :: fact
      real*8 :: ddot
      real*8 :: e3test_value, e3test_old

      real*8, dimension(*) :: wrk
      real*8, dimension(1) :: tmpwrk
      real*8, dimension(*) :: cmo, xindx

      real*8, dimension(kzyva) :: etrs
      real*8, dimension(kzyvb) :: vecb
      real*8, dimension(kzyvc) :: vecc
      real*8, dimension(kzyva) :: veca

      real*8, dimension(norbt,norbt) :: zymb, zymc
      real*8, dimension(nashdi,nashdi) :: udv, den1
      real*8, dimension(nnashx) :: dv

      integer, dimension(2,maxwop,8) :: mjwop

      real*8, dimension(:), allocatable :: fpe
      real*8, dimension(:), allocatable :: cref
      real*8, dimension(:), allocatable :: dcaos, fcaos
      real*8, dimension(:), allocatable :: udtv, udtvao
      real*8, dimension(:), allocatable :: dvaao, dvbao, dvatr
      real*8, dimension(:), allocatable :: udcao, udvao
      real*8, dimension(:), allocatable :: udcmo, udvmo
      real*8, dimension(:), allocatable :: fcmo

      real*8, dimension(:,:), allocatable :: dva, dvb
      real*8, dimension(:,:), allocatable :: fupe
      real*8, dimension(:,:), allocatable :: fxpeb,fxpec, fx2pe
      real*8, dimension(:,:), allocatable :: fxo1k, fxc1s
      real*8, dimension(:,:), allocatable :: fxo2k, fxc2s
      !real*8, dimension(:,:), allocatable :: fcas2_1, fcas2_2
      !real*8, dimension(:,:), allocatable :: fcas3_1, fcas3_2
      real*8, dimension(:,:), allocatable :: fxo, fxo1k2k, fxo2k1k
      real*8, dimension(:,:), allocatable :: f1kxo2k, f1sxc2k, f2kxo1k, f2sxc1k
      real*8, dimension(:,:), allocatable :: fxo1s2s
      real*8, dimension(:,:), allocatable :: fxc1s2s
      real*8, dimension(:,:), allocatable :: fxc1s2k, fxc1k2s

      logical :: lexist, lopen, lcon, lorb
      logical :: fndlab
      logical :: atest

      lopen = .false.

      call qenter('pe_rspmcqr')

      call gtzymt(1, vecb, kzyvb, isymb, zymb, mjwop)
      call gtzymt(1, vecc, kzyvc, isymc, zymc, mjwop)

      if (atest) then
         e3test_old = ddot(kzyva,veca,1,etrs,1)
         write(lupri,*) 'PE-E3TEST, on entry', e3test_old
      end if
      !-----------------------------------------------------------
      ! Get Fg = Vmul - R*<0|F|>Fe from file
      !-----------------------------------------------------------
      if (.not. tdhf) then
         allocate(fpe(nnorbx))
         if (lusifc <= 0) then
             call gpopen(lusifc, 'SIRIFC', 'OLD', ' ', 'UNFORMATTED',&
                         idummy, .false.)
                 lopen = .true.
         end if
         rewind(lusifc)
         call mollab('PEFMAT  ', lusifc, lupri)
         call readt(lusifc, nnorbx, fpe)
         if (lopen) call gpclose(lusifc, 'KEEP')
         allocate(fupe(norbt,norbt))
         call dsptsi(norbt, fpe, fupe)
         deallocate(fpe)
      end if
      !-----------------------------------------------------------
      ! Density Factory ...
      !-----------------------------------------------------------
      allocate(cref(mzconf(1)))
      call getref(cref, mzconf(1))

      if (pe_polar) then
         allocate(udcao(n2basx))
         allocate(udvao(n2basx))
         if (.not. tdhf) then
            allocate(dcaos(10*nnbasx))
            else
            allocate(dcaos(4*nnbasx))
         end if
         dcaos = 0.0d0
         !  DTX = D_pq(k1) = <0|[k1,Epq]|0>
         allocate(udcmo(n2orbx),udvmo(n2orbx))
         udcmo = 0.0d0
         udvmo = 0.0d0
         call deq27mo(isymb, zymb, udv, udcmo, udvmo, wrk, lfree)
         if (nasht > 0) then
            udcmo = udcmo + udvmo
         end if
         udcao = 0.0d0
         call motoao(udcmo,udcao,cmo,isymb,wrk,lfree)
         call dgefsp(nbast, udcao, dcaos(1:nnbasx))
         ! needed to fit with HF code
         dcaos(1:nnbasx) = 0.5d0*dcaos(1:nnbasx)
         !  DT2X = D_pq(k2,k1) = <0|[k2,[k1,Epq]|0>
         udvmo = 0.0d0
         call oitd1(isymc,zymc,udcmo,udvmo,isymb)
         ! DT2X in udvmo (re-used to save memory)
         udcao = 0.0d0
         isymbc = muld2h(isymb,isymc)
         call motoao(udvmo,udcao,cmo,isymbc,wrk,lfree)
         call dgefsp(nbast, udcao, dcaos(nnbasx+1:2*nnbasx))
         ! needed to fit with HF code
         dcaos(nnbasx+1:2*nnbasx) = 0.5d0*dcaos(nnbasx+1:2*nnbasx)
         !  DTX = D_pq(k2) = <0|[k2,Epq]|0>
         udcmo = 0.0d0
         udvmo = 0.0d0
         call deq27mo(isymc, zymc, udv, udcmo, udvmo,wrk, lfree)
         if (nasht > 0) then
            udcmo = udcmo + udvmo
         end if
         udcao = 0.0d0
         call motoao(udcmo,udcao,cmo,isymb,wrk,lfree)
         call dgefsp(nbast, udcao, dcaos(2*nnbasx+1:3*nnbasx))
         ! needed to fit with HF code
         dcaos(2*nnbasx+1:3*nnbasx) = 0.5d0*dcaos(2*nnbasx+1:3*nnbasx)
         !  DT2X = D_pq(k1,k2) = <0|[k1,[k2,Epq]|0>
         udvmo = 0.0d0
         call oitd1(isymb,zymb,udcmo,udvmo,isymc)
         ! DT2X in udvmo (re-used to save memory)
         udcao = 0.0d0
         isymbc = muld2h(isymc,isymb)
         call motoao(udvmo,udcao,cmo,isymbc,wrk,lfree)
         call dgefsp(nbast, udcao, dcaos(3*nnbasx+1:4*nnbasx))
         ! needed to fit with HF code
         dcaos(3*nnbasx+1:4*nnbasx) = 0.5d0*dcaos(3*nnbasx+1:4*nnbasx)
         deallocate(udcmo,udvmo)
         if (tdhf) then
            write(lupri,*) 'PE-DFT or HF QR detected: Skipping CI dens.'
         end if
         if (.not. tdhf ) then
            write(lupri,*) 'PE-MCSCF QR detected: Constructing CI dens.'
            ! Construct the density matrix <02L|..|0> + <0|..|02R>
            ilsym  = irefsy
            irsym  = muld2h(irefsy,isymc)
            ncl    = mzconf(1)
            ncr    = mzconf(isymc)
            kzvarl = mzconf(1)
            kzvarr = mzyvar(isymc)

            den1 = 0.0d0 ! edh: This is equal to udtv later...
            allocate(udtv(n2ashx), udtvao(n2basx))

            udtv = 0.0d0
            udtvao = 0.0d0
            call rspgdm(1, ilsym, irsym, ncl, ncr, kzvarl, kzvarr,&
                        cref, vecc, ovlap, udtv, dummy, 0 ,0, .true.,&
                        .true., xindx, wrk, 1, lfree, .true.)
            call fckden2(.false.,.true., dummy, udtvao, cmo,&
                         udtv, wrk, lfree)
            call dgefsp(nbast, udtvao, dcaos(4*nnbasx+1:5*nnbasx))
            dcaos(4*nnbasx+1:5*nnbasx) = 1.0d0*dcaos(4*nnbasx+1:5*nnbasx)

            ! Construct the density matrix <01L|..|0> + <0|..|01R>
            ilsym  = irefsy
            irsym  = muld2h(irefsy,isymb)
            ncl    = mzconf(1)
            ncr    = mzconf(isymb)
            kzvarl = mzconf(1)
            kzvarr = mzyvar(isymb)

            udtv = 0.0d0
            udtvao = 0.0d0
            call rspgdm(1, ilsym, irsym, ncl, ncr, kzvarl, kzvarr,&
                        cref, vecb, ovlap, udtv, dummy, 0 ,0, .true.,&
                        .true., xindx, wrk, 1, lfree, .true.)
            call fckden2(.false.,.true., dummy, udtvao, cmo,&
                         udtv, wrk, lfree)
            call dgefsp(nbast, udtvao, dcaos(5*nnbasx+1:6*nnbasx))
            dcaos(5*nnbasx+1:6*nnbasx) = 1.0d0*dcaos(5*nnbasx+1:6*nnbasx)

            if (mzconf(isymb) .gt. 0 .and. mzconf(isymc) .gt. 0) then
               ! Construct <01L|..|02R> + <02L|..|01R> density
               ilsym  = muld2h(irefsy,isymb)
               irsym  = muld2h(irefsy,isymc)
               ncl    = mzconf(isymb)
               ncr    = mzconf(isymc)
               kzvarl = mzyvar(isymb)
               kzvarr = mzyvar(isymc)
               isymdn = muld2h(ilsym,irsym)

               udtv = 0.0d0
               udtvao = 0.0d0
               call rspgdm(1, ilsym, irsym, ncl, ncr, kzvarl, kzvarr,&
                           vecb, vecc, ovlap, udtv, dummy, 0 ,0,&
                          .true., .true., xindx, wrk, 1, lfree,&
                          .false.)
               call fckden2(.false.,.true., dummy, udtvao, cmo,&
                            udtv, wrk, lfree)
               call dgefsp(nbast, udtvao, dcaos(6*nnbasx+1:7*nnbasx))
            end if
            dcaos(6*nnbasx+1:7*nnbasx) = 1.0d0*dcaos(6*nnbasx+1:7*nnbasx)

            ! D_pq = <0|Epq|0>
            udcao = 0.0d0
            udvao = 0.0d0
            call dgetsp(nasht, udv, dv)
            call fckden((nisht>0), (nasht>0), udcao, udvao,&
                        cmo, dv, wrk, lfree)
            if (nisht==0) udcao = 0.0d0
            udcao = udcao + udvao
            call dgefsp(nbast, udcao, dcaos(7*nnbasx+1:8*nnbasx))
            dcaos(7*nnbasx+1:8*nnbasx) = 0.5d0*dcaos(7*nnbasx+1:8*nnbasx)

            ! D_pq(S1,k2) = <01L|[k2,Epq]|0> + <0|[k2,Epq]|01R>
            allocate(dva(norbt,nasht), dvb(norbt,nasht))
            allocate(dvaao(n2basx), dvbao(n2basx), dvatr(n2basx))
            dva    = 0.0d0
            dvb    = 0.0d0
            dvaao  = 0.0d0
            dvbao  = 0.0d0
            dvatr  = 0.0d0
            udtvao = 0.0d0
            call rsptr1(1, udv, zymb, dva, dvb)
            call fckden2(.false.,.true., dummy, dvaao, cmo,&
                         dva, wrk, lfree)
            call fckden2(.false.,.true., dummy, dvbao, cmo,&
                         dvb, wrk, lfree)
            call mtrsp(nbast, nbast, dvaao, nbast, dvatr, nbast)
            udtvao = dvbao - dvatr
            call dgefsp(nbast, udtvao, dcaos(8*nnbasx+1:9*nnbasx))
            dcaos(8*nnbasx+1:9*nnbasx) = 2.0d0*dcaos(8*nnbasx+1:9*nnbasx)

            ! D_pq(k1,S2) = <02L|[k1,Epq]|0> + <0|[k1,Epq]|02R>
            dva    = 0.0d0
            dvb    = 0.0d0
            dvaao  = 0.0d0
            dvbao  = 0.0d0
            dvatr  = 0.0d0
            udtvao = 0.0d0
            call rsptr1(1, udv, zymc, dva, dvb)
            call fckden2(.false.,.true., dummy, dvaao, cmo,&
                         dva, wrk, lfree)
            call fckden2(.false.,.true., dummy, dvbao, cmo,&
                         dvb, wrk, lfree)
            call mtrsp(nbast, nbast, dvaao, nbast, dvatr, nbast)
            udtvao = dvbao - dvatr
            call dgefsp(nbast, udtvao, dcaos(9*nnbasx+1:10*nnbasx))
            dcaos(9*nnbasx+1:10*nnbasx) = 2.0d0*dcaos(9*nnbasx+1:10*nnbasx)

            deallocate(dva, dvb, dvaao, dvbao, dvatr)
            deallocate(udtv,udtvao)

         end if
         deallocate(udcao, udvao)
         write(lupri,*)'entering PElib (pe_master)'
         !-----------------------------------------------------------
         ! Calculate PE response operators in AO basis
         !-----------------------------------------------------------
         if (.not. tdhf) then
            allocate(fcaos(10*nnbasx))
            fcaos = 0.0d0
#if defined(VAR_MPI)
            call pelib_ifc_start_slaves(3)
#endif
            call pe_master(runtype='dynamic_response', &
                           triang=.true., &
                           ndim=nbast, &
                           nmats=10, &
                           denmats=dcaos, &
                           fckmats=fcaos)
         else
            allocate(fcaos(4*nnbasx))
            fcaos = 0.0d0
#if defined(VAR_MPI)
            call pelib_ifc_start_slaves(3)
#endif
            call pe_master(runtype='dynamic_response', &
                           triang=.true., &
                           ndim=nbast, &
                           nmats=4, &
                           denmats=dcaos, &
                           fckmats=fcaos)
         end if
         deallocate(dcaos)
      end if ! pe_polar

      if ( .not. tdhf ) then
         !-----------------------------------------------------------
         !case 1
         !-----------------------------------------------------------
         if ( mzconf(isymb) .eq. 0 .or. mzconf(isymc) .eq. 0 ) return
         !/   <01L| [qj,TB] |02R>  + <02L| [qj,TB] |01R>  \
         !|                       0                       |
         !|   <01L| [qj+,TB] |02R> + <02L| [qj+,TB] |01R> |
         !\                       0                       /

         ! ionstruct <01L|..|02R> + <02L|..|01R> density
         ilsym  = muld2h(irefsy,isymb)
         irsym  = muld2h(irefsy,isymc)
         ncl    = mzconf(isymb)
         ncr    = mzconf(isymc)
         kzvarl = mzyvar(isymb)
         kzvarr = mzyvar(isymc)

         den1 = 0.0d0
         call rspgdm(1, ilsym, irsym, ncl, ncr, kzvarl, kzvarr,&
                     vecb, vecc, ovlap, den1, dummy, 0, 0, .true.,&
                     .true., xindx, wrk, 1, lfree, .false.)

         ! Make the gradient
         isymdn = muld2h(ilsym,irsym)

         if ( mzwopt(isyma) .gt. 0 ) then
            call orbsx(1, isyma, kzyva, etrs, fupe, ovlap,&
                       isymdn, den1, mjwop, 1, lfree)
            write(lupri,*)'atest:', atest
         end if
      endif
      if (atest) then
         e3test_value = ddot(kzyva,veca,1,etrs,1)
         write(lupri,*) 'PE-E3TEST, case 1 ', e3test_value-e3test_old
         e3test_old = e3test_value
      end if
      !-----------------------------------------------------------
      !case 2
      !-----------------------------------------------------------
      if (pe_polar) then
         allocate(fcmo(nnorbx))
         allocate(fxo1k(norbt,norbt))
         fxo1k = 0.0d0
         if (.not. tdhf) then
            allocate(fxc1s(norbt,norbt))
            fxc1s = 0.0d0
         end if
         ! Fxo = R*<0|Fe(1k)|0>Fe
         fcmo = 0.0d0
         call uthu(2.0d0*fcaos(1:nnbasx), fcmo, cmo,&
                   wrk, nbast, norbt)
         call dsptsi(norbt, fcmo, fxo1k)
         if (.not. tdhf) then
            if (mzconf(isymc) .le. 0) return
            !/   0    \
            !| Sj(2)  | * <0| Fxo[1](1k) |0>
            !|   0    |
            !\ Sj(2)* /
            if (isyma .eq. isymc) then
                fact  = 0.0d0
                ovlap = 1.0d0
                call melone(fxo1k, 1, udv, ovlap, fact,&
                            200,'fact for Fxo(1k)')
                nzconf = mzconf(isyma)
                nzvar  = mzvar(isyma)
                call daxpy(nzconf, fact, vecc, 1, etrs, 1)
                call daxpy(nzconf, fact, vecc(nzvar+1), 1, etrs(nzvar+1), 1)
            end if
         end if
         ! For testing the Fxo(1k) term
         if (atest) then
            e3test_value = ddot(kzyva,veca,1,etrs,1)
            write(lupri,*) 'PE-E3TEST, case 2 Fxo(1k)', e3test_value-e3test_old
            e3test_old = e3test_value
         end if
         ! Fxc(1S) = ( R*<01lE|0>+<0|E01R> )Fe
         if (.not. tdhf) then
            ! edh: Should it be 1.0d0 or 2.0d0 ???
            fcmo = 0.0d0
            call uthu(1.0d0*fcaos(4*nnbasx+1:5*nnbasx), fcmo, cmo,&
                      wrk, nbast, norbt)
            call dsptsi(norbt, fcmo, fxc1s)
         end if
         if (.not. tdhf) then
            if (mzconf(isymc) .le. 0) return
            !/   0    \
            !| Sj(2)  | * <0| Fxc(1S) |0>
            !|   0    |
            !\ Sj(2)* /
            if (isyma .eq. isymc) then
                fact = 0.0d0
                ovlap = 1.0d0
                call melone(fxc1s, 1, udv, ovlap, fact,&
                            200,'fact for Fxc(1S) ')
                nzconf = mzconf(isyma)
                nzvar  = mzvar(isyma)
                call daxpy(nzconf, fact, vecc, 1, etrs, 1)
                call daxpy(nzconf, fact, vecc(nzvar+1), 1, etrs(nzvar+1), 1)
            end if
         endif
         if (atest) then
            e3test_value = ddot(kzyva,veca,1,etrs,1)
            write(lupri,*) 'PE-E3TEST case 2 Fxc(1S)', e3test_value-e3test_old
            e3test_old = e3test_value
         end if
         !
         ! second part of case 2
         allocate(fxo2k(norbt,norbt))
         fxo2k = 0.0d0
         if (.not. tdhf) then
            allocate(fxc2s(norbt,norbt))
            fxc2s = 0.0d0
         end if
         ! Fxo(2k) = R*<0|[2k,Epq]|0>Fe
         fcmo = 0.0d0
         call uthu(2.0d0*fcaos(2*nnbasx+1:3*nnbasx),&
                   fcmo, cmo, wrk, nbast, norbt)
         call dsptsi(norbt, fcmo, fxo2k)
         if (.not. tdhf) then
            if (mzconf(isymb) .le. 0) return
            !/   0    \
            !| Sj(1)  | * <0| Fxo(2k) |0>
            !|   0    |
            !\ Sj(1)* /
            if (isyma .eq. isymb) then
                ovlap = 1.0d0
                fact = 0.0d0
                call melone(fxo2k, 1, udv, ovlap, fact,&
                            200,'fact for Fxo(2k)')
                  nzconf = mzconf(isyma)
                  nzvar  = mzvar(isyma)
                  call daxpy(nzconf, fact, vecb, 1, etrs, 1)
                  call daxpy(nzconf, fact, vecb(nzvar+1), 1, etrs(nzvar+1), 1)
            end if
         end if
         if (atest) then
            e3test_value = ddot(kzyva,veca,1,etrs,1)
            write(lupri,*) 'PE-E3TEST, Fxo(2k)', e3test_value-e3test_old
            e3test_old = e3test_value
         end if
         ! Fxc(2S) = ( R*<02LE|0>+<0|E02R> )Fe
         if (.not. tdhf) then
            ! edh: Should it be 1.0d0 or 2.0d0 ???
            fcmo = 0.0d0
            call uthu(1.0d0*fcaos(5*nnbasx+1:6*nnbasx), fcmo, cmo,&
                      wrk, nbast, norbt)
             call dsptsi(norbt, fcmo, fxc2s)
         end if
         if (.not. tdhf) then
            if (mzconf(isymb) .le. 0) return
            !/   0    \
            !| Sj(1)  | * <0| Fxc(2S) |0>
            !|   0    |
            !\ Sj(1)* /
            if (isyma .eq. isymb) then
                ovlap = 1.0d0
                fact = 0.0d0
                call melone(fxc2s, 1, udv, ovlap, fact,&
                            200,'fact for Fxc(2S) ')
                  nzconf = mzconf(isyma)
                  nzvar  = mzvar(isyma)
                  call daxpy(nzconf, fact, vecb, 1, etrs, 1)
                  call daxpy(nzconf, fact, vecb(nzvar+1), 1, etrs(nzvar+1), 1)
            end if
         end if
         if (atest) then
            e3test_value = ddot(kzyva,veca,1,etrs,1)
            write(lupri,*) 'PE-E3TEST, fxc(2S)', e3test_value-e3test_old
            e3test_old = e3test_value
         end if

         deallocate(fcmo)
      endif ! pe_polar
      !-----------------------------------------------------------
      !case 3
      !-----------------------------------------------------------
      allocate(fxpeb(norbt,norbt))
      fxpeb = 0.0d0
      call oith1(isymb,zymb,fupe,fxpeb,1)
      if (.not. tdhf) then
         if (mzconf(isymc) .le. 0) return
         !/   <0| [qj,Fg(1k)] |02R>  + <02L| [qj,Fg(1k)] |0>  \
         !|   <j| Fg(1k) |02R>                                |
         !|   <0| [qj+,Fg(1k)] |02R> + <02L| [qj+,Fg(1k)] |0> |
         !\  -<02L| Fg(1k) |j>                                /
         !Construct the density matrix <02L|..|0> + <0|..|02R>
         ilsym  = irefsy
         irsym  = muld2h(irefsy,isymc)
         ncl    = mzconf(1)
         ncr    = mzconf(isymc)
         kzvarl = mzconf(1)
         kzvarr = mzyvar(isymc)
         den1 = 0.0d0
         call rspgdm(1, ilsym, irsym, ncl, ncr, kzvarl, kzvarr,&
                     cref, vecc, ovlap, den1, dummy, 0 ,0, .true.,&
                     .true., xindx, wrk, 1, lfree, .true.)
         !Make the gradient
         isymdn = muld2h(ilsym,irsym)
         isymst = muld2h(isyma,irefsy)
         if ( isymst .eq. irefsy ) then
            lcon = ( mzconf(isyma) .gt. 1 )
         else
            lcon = ( mzconf(isyma) .gt. 0 )
         end if
         lorb    = ( mzwopt(isyma) .gt. 0 )
         nzyvec = mzyvar(isymc)
         nzcvec = mzconf(isymc)
         call rsp1gr(1, kzyva, idummy, 0 , isyma, 0, isymc, etrs,&
                     vecc, nzyvec, nzcvec, ovlap, isymdn, den1, fxpeb,&
                     xindx, mjwop, wrk(1), lfree, lorb, lcon, .false.)
      end if
      if (atest) then
         e3test_value = ddot(kzyva,veca,1,etrs,1)
         write(lupri,*) 'PE-E3TEST, case 3 Fg(1k)', e3test_value-e3test_old
         e3test_old = e3test_value
      end if
      if (.not. tdhf) then
         if (mzconf(isymc) .le. 0) return
         !/   <0| [qj,Fxo(1k)] |02R>  + <02L| [qj,Fxo(1k)] |0>  \
         !|   <j| Fxo(1k) |02R>                                 |
         !|   <0| [qj+,Fxo(1k)] |02R> + <02L| [qj+,Fxo(1k)] |0> |
         !\  -<02L| Fxo(1k) |j>                                 /
         ! Construct the density matrix <02L|..|0> + <0|..|02R>
         ilsym  = irefsy
         irsym  = muld2h(irefsy,isymc)
         ncl    = mzconf(1)
         ncr    = mzconf(isymc)
         kzvarl = mzconf(1)
         kzvarr = mzyvar(isymc)
         den1 = 0.0d0
         call rspgdm(1, ilsym, irsym, ncl, ncr, kzvarl, kzvarr,&
                     cref, vecc, ovlap, den1, dummy, 0 ,0, .true.,&
                     .true., xindx, wrk, 1, lfree, .true.)
         ! Make the gradient
         isymdn = muld2h(ilsym,irsym)
         isymst = muld2h(isyma,irefsy)
         if ( isymst .eq. irefsy ) then
            lcon = ( mzconf(isyma) .gt. 1 )
         else
            lcon = ( mzconf(isyma) .gt. 0 )
         end if
         lorb    = ( mzwopt(isyma) .gt. 0 )
         nzyvec = mzyvar(isymc)
         nzcvec = mzconf(isymc)
         call rsp1gr(1, kzyva, idummy, 0 , isyma, 0, isymc, etrs,&
                     vecc, nzyvec, nzcvec, ovlap, isymdn, den1, fxo1k,&
                     xindx, mjwop, wrk(1), lfree, lorb, lcon, .false.)
      end if
      if (atest) then
         e3test_value = ddot(kzyva,veca,1,etrs,1)
         write(lupri,*) 'PE-E3TEST, case 3 fxo1k', e3test_value-e3test_old
         e3test_old = e3test_value
      end if
      if (.not. tdhf) then
         if (mzconf(isymc) .le. 0) return
         !/   <0| [qj,Fxc(1S)] |02R>  + <02L| [qj,Fxc(1S)] |0>  \
         !|   <j| Fxc(1S) |02R>                                 |
         !|   <0| [qj+,Fxc(1S)] |02R> + <02L| [qj+,Fxc(1S)] |0> |
         !\  -<02L| Fxc(1S) |j>                                 /
         ! Construct the density matrix <02L|..|0> + <0|..|02R>
         ilsym  = irefsy
         irsym  = muld2h(irefsy,isymc)
         ncl    = mzconf(1)
         ncr    = mzconf(isymc)
         kzvarl = mzconf(1)
         kzvarr = mzyvar(isymc)
         den1 = 0.0d0
         call rspgdm(1, ilsym, irsym, ncl, ncr, kzvarl, kzvarr,&
                     cref, vecc, ovlap, den1, dummy, 0 ,0, .true.,&
                     .true., xindx, wrk, 1, lfree, .true.)
         !Make the gradient
         isymdn = muld2h(ilsym,irsym)
         isymst = muld2h(isyma,irefsy)
         if ( isymst .eq. irefsy ) then
            lcon = ( mzconf(isyma) .gt. 1 )
         else
            lcon = ( mzconf(isyma) .gt. 0 )
         end if
         lorb    = ( mzwopt(isyma) .gt. 0 )
         nzyvec = mzyvar(isymc)
         nzcvec = mzconf(isymc)
         call rsp1gr(1, kzyva, idummy, 0 , isyma, 0, isymc, etrs,&
                     vecc, nzyvec, nzcvec, ovlap, isymdn, den1, fxc1s,&
                     xindx, mjwop, wrk(1), lfree, lorb, lcon, .false.)
      end if
      if (atest) then
         e3test_value = ddot(kzyva,veca,1,etrs,1)
         write(lupri,*) 'PE-E3TEST, case 3 fxc1s', e3test_value-e3test_old
         e3test_old = e3test_value
      end if
      ! second part
      allocate(fxpec(norbt,norbt))
      fxpec = 0.0d0
      call oith1(isymc,zymc,fupe,fxpec,1)
      if (.not. tdhf) then
         if (mzconf(isymb) .le. 0) return
         !/   <0| [qj,F(2k)] |01R>  + <01L| [qj,Fg(2k)] |0>   \
         !|   <j|Fg(2k) |01R>                                 |
         !|   <0| [qj+,Fg(2k)] |01R> + <01L| [qj+,Fg(2k)] |0> |
         !\  -<01L| Fg(2k) |j>                                /
         ! Construct the density matrix <01L|..|0> + <0|..|01R>
         ilsym  = irefsy
         irsym  = muld2h(irefsy,isymb)
         ncl    = mzconf(1)
         ncr    = mzconf(isymb)
         kzvarl = mzconf(1)
         kzvarr = mzyvar(isymb)
         den1 = 0.0d0
         call rspgdm(1, ilsym, irsym, ncl, ncr, kzvarl, kzvarr,&
                     cref, vecb, ovlap, den1, dummy, 0 ,0, .true.,&
                     .true., xindx, wrk, 1, lfree, .true.)
         ! Make the gradient
         isymdn = muld2h(ilsym,irsym)
         isymst = muld2h(isyma,irefsy)
         if ( isymst .eq. irefsy ) then
            lcon = ( mzconf(isyma) .gt. 1 )
         else
            lcon = ( mzconf(isyma) .gt. 0 )
         end if
         lorb    = ( mzwopt(isyma) .gt. 0 )
         nzyvec = mzyvar(isymb)
         nzcvec = mzconf(isymb)
         call rsp1gr(1, kzyva, idummy, 0 , isyma, 0, isymb, etrs,&
                     vecb, nzyvec, nzcvec, ovlap, isymdn, den1, fxpec,&
                     xindx, mjwop, wrk(1), lfree, lorb, lcon, .false.)
      end if
      if (atest) then
         e3test_value = ddot(kzyva,veca,1,etrs,1)
         write(lupri,*) 'PE-E3TEST, case 3 Fg(2k)', e3test_value-e3test_old
         e3test_old = e3test_value
      end if
      if (.not. tdhf) then
         if (mzconf(isymb) .le. 0) return
         !/   <0| [qj,Fxo(2k)] |01R>  + <01L| [qj,Fxo(2k)] |0>  \
         !|   <j| Fxo(2k) |01R>                                 |
         !|   <0| [qj+,Fxo(2k)] |01R> + <01L| [qj+,Fxo(2k)] |0> |
         !\  -<01L| Fxo(2k) |j>                                 /
         !Construct the density matrix <01L|..|0> + <0|..|01R>
         ilsym  = irefsy
         irsym  = muld2h(irefsy,isymb)
         ncl    = mzconf(1)
         ncr    = mzconf(isymb)
         kzvarl = mzconf(1)
         kzvarr = mzyvar(isymb)
         den1 = 0.0d0
         call rspgdm(1, ilsym, irsym, ncl, ncr, kzvarl, kzvarr,&
                     cref, vecb, ovlap, den1, dummy, 0 ,0, .true.,&
                     .true., xindx, wrk, 1, lfree, .true.)
         !Make the gradient
         isymdn = muld2h(ilsym,irsym)
         isymst = muld2h(isyma,irefsy)
         if ( isymst .eq. irefsy ) then
            lcon = ( mzconf(isyma) .gt. 1 )
         else
            lcon = ( mzconf(isyma) .gt. 0 )
         end if
         lorb    = ( mzwopt(isyma) .gt. 0 )
         nzyvec = mzyvar(isymb)
         nzcvec = mzconf(isymb)
         call rsp1gr(1, kzyva, idummy, 0 , isyma, 0, isymb, etrs,&
                     vecb, nzyvec, nzcvec, ovlap, isymdn, den1, fxo2k,&
                     xindx, mjwop, wrk(1), lfree, lorb, lcon, .false.)
      end if
      if (atest) then
         e3test_value = ddot(kzyva,veca,1,etrs,1)
         write(lupri,*) 'PE-E3TEST, case fxo2k', e3test_value-e3test_old
         e3test_old = e3test_value
      end if
      if (.not. tdhf) then
         if (mzconf(isymb) .le. 0) return
         !/   <0| [qj,Fxc(2S)] |01R>  + <01L| [qj,Fxc(2S)] |0>  \
         !|   <j| Fxc(2S) |01R>                                 |
         !|   <0| [qj+,Fxc(2S)] |01R> + <01L| [qj+,Fxc(2S)] |0> |
         !\  -<01L| Fxc(2S) |j>                                 /
         ! Construct the density matrix <01L|..|0> + <0|..|01R>
         ilsym  = irefsy
         irsym  = muld2h(irefsy,isymb)
         ncl    = mzconf(1)
         ncr    = mzconf(isymb)
         kzvarl = mzconf(1)
         kzvarr = mzyvar(isymb)
         den1 = 0.0d0
         call rspgdm(1, ilsym, irsym, ncl, ncr, kzvarl, kzvarr,&
                     cref, vecb, ovlap, den1, dummy, 0 ,0, .true.,&
                     .true., xindx, wrk, 1, lfree, .true.)
         !Make the gradient
         isymdn = muld2h(ilsym,irsym)
         isymst = muld2h(isyma,irefsy)
         if ( isymst .eq. irefsy ) then
            lcon = ( mzconf(isyma) .gt. 1 )
         else
            lcon = ( mzconf(isyma) .gt. 0 )
         end if
         lorb    = ( mzwopt(isyma) .gt. 0 )
         nzyvec = mzyvar(isymb)
         nzcvec = mzconf(isymb)
         call rsp1gr(1, kzyva, idummy, 0 , isyma, 0, isymb, etrs,&
                     vecb, nzyvec, nzcvec, ovlap, isymdn, den1, fxc2s,&
                     xindx, mjwop, wrk(1), lfree, lorb, lcon, .false.)
            if (atest) then
               e3test_value = ddot(kzyva,veca,1,etrs,1)
               write(lupri,*) 'PE-E3TEST, case fxc2s', e3test_value-e3test_old
               e3test_old = e3test_value
            end if
      end if
!     !-----------------------------------------------------------
!     !case 4
!     !-----------------------------------------------------------
!     !fx2pe = 0.5*Fg(1k,2k) + 0.5*Fg(2k,1k) + ...
      allocate(fx2pe(norbt,norbt))
      fx2pe = 0.0d0
      if (.not. tdhf) then
         call oith1(isymc, zymc, fxpeb, fx2pe, isymb)
         call oith1(isymb, zymb, fxpec, fx2pe, isymc)
         deallocate(fxpeb,fxpec)
      end if
      fx2pe = 0.5d0 * fx2pe
      !/ <0| [qj, Fg(1k,2k) + Fg(2k,1k)] |0>  \
      !| <j| Fg(1k,2k) + Fg(2k,1k) |0>        |*0.5
      !| <0| [qj+,Fg(1k,2k) + Fg(2k,1k) ] |0> |
      !\-<0| Fg(1k,2k) + Fg(2k,1k) |j>        /
      isymdn = 1
      ovlap  = 1.0d0
      isymst = muld2h(isyma, irefsy)
      if ( isymst .eq. irefsy ) then
         lcon = ( mzconf(isyma) .gt. 1 )
      else
         lcon = ( mzconf(isyma) .gt. 0 )
      end if
      lorb   = ( mzwopt(isyma) .gt. 0 )
      nzyvec = mzconf(1)
      nzcvec = mzconf(1)
      call rsp1gr(1 ,kzyva, idummy,0, isyma, 0, irefsy, etrs,&
                  cref, nzyvec, nzcvec, ovlap, isymdn, udv, fx2pe,&
                  xindx, mjwop, wrk(1), lfree, lorb, lcon, .true.)
      deallocate(fx2pe)
      if (atest) then
         e3test_value = ddot(kzyva,veca,1,etrs,1)
         write(lupri,*) 'PE-E3TEST, case 4 Fg(1k,2k) + F(2k,1k) ', e3test_value-e3test_old
         write(lupri,*) 'PE-E3TEST, total  ', e3test_value
      end if
      if (pe_polar) then
          allocate(f1kxo2k(norbt,norbt))
          f1kxo2k = 0.0d0
          call oith1(isymc,zymc,fxo1k,f1kxo2k,isymb)
          !/ <0| [qj, Fxo[1k](2k)] |0> \
          !| <j| Fxo[1k](2k) |0>       |
          !| <0| [qj+,Fxo[1k](2k)] |0> |
          !\ -<0| Fxo[1k](2k) |j>      /
          isymdn = 1
          ovlap  = 1.0d0
          isymst = muld2h(isyma, irefsy)
          if ( isymst .eq. irefsy ) then
             lcon = ( mzconf(isyma) .gt. 1 )
          else
             lcon = ( mzconf(isyma) .gt. 0 )
          end if
          lorb   = ( mzwopt(isyma) .gt. 0 )
          nzyvec = mzconf(1)
          nzcvec = mzconf(1)
          call rsp1gr(1 ,kzyva, idummy,0, isyma, 0, irefsy, etrs,&
                      cref, nzyvec, nzcvec, ovlap, isymdn, udv, f1kxo2k,&
                      xindx, mjwop, wrk(1), lfree, lorb, lcon, .true.)
          deallocate(f1kxo2k)
          if (atest) then
             e3test_value = ddot(kzyva,veca,1,etrs,1)
             write(lupri,*) 'PE-E3TEST, case f1kxo2k ', e3test_value-e3test_old
             write(lupri,*) 'PE-E3TEST, total  ', e3test_value
          end if
          allocate(f1sxc2k(norbt,norbt))
          f1sxc2k = 0.0d0
          call oith1(isymc,zymc,fxc1s,f1sxc2k,isymb)
          !/ <0| [qj ,Fxc[1s](2k)] |0>  \
          !| <j| Fxc[1s](2k) |0>        |
          !| <0| [qj+,Fxc[1s](2k)] |0>  |
          !\ -<0| Fxc[1s](2k) |j>       /
          isymdn = 1
          ovlap  = 1.0d0
          isymst = muld2h(isyma, irefsy)
          if ( isymst .eq. irefsy ) then
             lcon = ( mzconf(isyma) .gt. 1 )
          else
             lcon = ( mzconf(isyma) .gt. 0 )
          end if
          lorb   = ( mzwopt(isyma) .gt. 0 )
          nzyvec = mzconf(1)
          nzcvec = mzconf(1)
          call rsp1gr(1 ,kzyva, idummy,0, isyma, 0, irefsy, etrs,&
                      cref, nzyvec, nzcvec, ovlap, isymdn, udv, f1sxc2k,&
                      xindx, mjwop, wrk(1), lfree, lorb, lcon, .true.)
          deallocate(f1sxc2k)
          if (atest) then
             e3test_value = ddot(kzyva,veca,1,etrs,1)
             write(lupri,*) 'PE-E3TEST, case 4 f1sxc2k ', e3test_value-e3test_old
             write(lupri,*) 'PE-E3TEST, total  ', e3test_value
          end if
          !construct case2 and add to gradient...
          allocate(f2kxo1k(norbt,norbt))
          f2kxo1k = 0.0d0
          call oith1(isymb,zymb,fxo2k,f2kxo1k,isymc)
          !/ <0| [qj ,Fxo[2k](1k)] |0> \
          !| <j| Fxo[2k](1k) |0>       |
          !| <0| [qj+,Fxo[2k](1k)] |0> |
          !\ -<0| Fxo[2k](1k) |j>      /
          isymdn = 1
          ovlap  = 1.0d0
          isymst = muld2h(isyma, irefsy)
          if ( isymst .eq. irefsy ) then
             lcon = ( mzconf(isyma) .gt. 1 )
          else
             lcon = ( mzconf(isyma) .gt. 0 )
          end if
          lorb   = ( mzwopt(isyma) .gt. 0 )
          nzyvec = mzconf(1)
          nzcvec = mzconf(1)
          call rsp1gr(1 ,kzyva, idummy,0, isyma, 0, irefsy, etrs,&
                      cref, nzyvec, nzcvec, ovlap, isymdn, udv, f2kxo1k,&
                      xindx, mjwop, wrk(1), lfree, lorb, lcon, .true.)
          deallocate(fxo2k,f2kxo1k)
          if (atest) then
             e3test_value = ddot(kzyva,veca,1,etrs,1)
             write(lupri,*) 'PE-E3TEST, case 4 f2kxo1k ', e3test_value-e3test_old
             write(lupri,*) 'PE-E3TEST, total  ', e3test_value
          end if
          allocate(f2sxc1k(norbt,norbt))
          f2sxc1k = 0.0d0
          call oith1(isymb,zymb,fxc2s,f2sxc1k,isymc)
          !/ <0| [qj ,Fxc[2s](1k)] |0> \
          !| <j| Fxc[2s](1k) |0>       |
          !| <0| [qj+,Fxc[2s](1k)] |0> |
          !\ -<0| Fxc[2s](1k) |j>      /
          isymdn = 1
          ovlap  = 1.0d0
          isymst = muld2h(isyma, irefsy)
          if ( isymst .eq. irefsy ) then
             lcon = ( mzconf(isyma) .gt. 1 )
          else
             lcon = ( mzconf(isyma) .gt. 0 )
          end if
          lorb   = ( mzwopt(isyma) .gt. 0 )
          nzyvec = mzconf(1)
          nzcvec = mzconf(1)
          call rsp1gr(1 ,kzyva, idummy,0, isyma, 0, irefsy, etrs,&
                      cref, nzyvec, nzcvec, ovlap, isymdn, udv, f2sxc1k,&
                      xindx, mjwop, wrk(1), lfree, lorb, lcon, .true.)
          deallocate(fxc2s,f2sxc1k)
          if (atest) then
             e3test_value = ddot(kzyva,veca,1,etrs,1)
             write(lupri,*) 'PE-E3TEST, case 4 f2sxc1k ', e3test_value-e3test_old
             write(lupri,*) 'PE-E3TEST, total  ', e3test_value
          end if
          ! + ( S(1)S*(2) + S(2)S*(1) ) * fxo  + ...
          if (.not. tdhf ) then
             if ((isymb .eq. isymc) .and. (mzconf(isymb) .gt. 0)) then
                allocate(fxo(norbt,norbt))
                fcmo = 0.0d0
                fxo  = 0.0d0
                fact = 0.0d0
                ! edh: was changed from 0.25 to 1.0 (check this factor!)
                call uthu(1.0d0*fcaos(7*nnbasx+1:8*nnbasx), fcmo, cmo,&
                          wrk, nbast, norbt)
                call dsptsi(norbt, fcmo, fxo)
                allocate(fxo1s2s(norbt,norbt))
                nzconf = mzconf(isymb)
                nzvar  = mzvar(isymb)
                fact   = ddot(nzconf, vecb, 1, vecc(nzvar+1), 1) + &
                         ddot(nzconf, vecc, 1, vecb(nzvar+1), 1)
                call daxpy(n2orbx, fact, fxo, 1, fxo1s2s, 1)
                deallocate(fxo)
             end if
          end if
          !/ <0| [qj ,TE] |0> \
          !| <j| TE |0>       |
          !| <0| [qj+,TE] |0> |
          !\ -<0| TE |j>      /
          isymdn = 1
          ovlap  = 1.0d0
          isymst = muld2h(isyma, irefsy)
          if ( isymst .eq. irefsy ) then
             lcon = ( mzconf(isyma) .gt. 1 )
          else
             lcon = ( mzconf(isyma) .gt. 0 )
          end if
          lorb   = ( mzwopt(isyma) .gt. 0 )
          nzyvec = mzconf(1)
          nzcvec = mzconf(1)
          call rsp1gr(1 ,kzyva, idummy,0, isyma, 0, irefsy, etrs,&
                      cref, nzyvec, nzcvec, ovlap, isymdn, udv, fxo1s2s,&
                      xindx, mjwop, wrk(1), lfree, lorb, lcon, .true.)
               deallocate(fxo1s2s)
               if (atest) then
                  e3test_value = ddot(kzyva,veca,1,etrs,1)
                  write(lupri,*) 'PE-E3TEST, case 4 fxo1s2s ', e3test_value-e3test_old
                  write(lupri,*) 'PE-E3TEST, total  ', e3test_value
               end if
          if (mzwopt(isymb).gt.0 .and. mzwopt(isymc).gt.0) then
             allocate(fxo1k2k(norbt,norbt))
             fxo1k2k = 0.0d0
             !fxo(1k,2k) = <0|Fe(1k,2k)|0>Fe
             fcmo = 0.0d0
             call uthu(1.0d0*fcaos(nnbasx+1:2*nnbasx),&
                       fcmo, cmo, wrk, nbast, norbt)
             call dsptsi(norbt, fcmo, fxo1k2k)
             !/ <0| [qj ,Fxo(1k,2k)] |0> \
             !| <j| Fxo(1k,2k) |0>       |
             !| <0| [qj+,Fxo(1k,2k)] |0> |
             !\ -<0| Fxo(1k,2k) |j>      /
             isymdn = 1
             ovlap  = 1.0d0
             isymst = muld2h(isyma, irefsy)
             if ( isymst .eq. irefsy ) then
                lcon = ( mzconf(isyma) .gt. 1 )
             else
                lcon = ( mzconf(isyma) .gt. 0 )
             end if
             lorb   = ( mzwopt(isyma) .gt. 0 )
             nzyvec = mzconf(1)
             nzcvec = mzconf(1)
             call rsp1gr(1 ,kzyva, idummy,0, isyma, 0, irefsy, etrs,&
                         cref, nzyvec, nzcvec, ovlap, isymdn, udv, fxo1k2k,&
                         xindx, mjwop, wrk(1), lfree, lorb, lcon, .true.)
             deallocate(fxo1k2k)
             if (atest) then
                e3test_value = ddot(kzyva,veca,1,etrs,1)
                write(lupri,*) 'PE-E3TEST, case fxo1k2k ', e3test_value-e3test_old
                write(lupri,*) 'PE-E3TEST, total  ', e3test_value
             end if
             allocate(fxo2k1k(norbt,norbt))
             fxo2k1k = 0.0d0
             !fxo(2k,1k) = R*<0|Fe(2k,1k)|0>Fe
             fcmo    = 0.0d0
             call uthu(1.0d0*fcaos(3*nnbasx+1:4*nnbasx),&
                       fcmo, cmo, wrk, nbast, norbt)
             call dsptsi(norbt, fcmo, fxo2k1k)
             !/ <0| [qj ,Fxo(2k,1k)] |0> \
             !| <j| Fxo(2k,1k) |0>       |
             !| <0| [qj+,Fxo(2k,1k)] |0> |
             !\ -<0| Fxo(2k,1k) |j>      /
             isymdn = 1
             ovlap  = 1.0d0
             isymst = muld2h(isyma, irefsy)
             if ( isymst .eq. irefsy ) then
                lcon = ( mzconf(isyma) .gt. 1 )
             else
                lcon = ( mzconf(isyma) .gt. 0 )
             end if
             lorb   = ( mzwopt(isyma) .gt. 0 )
             nzyvec = mzconf(1)
             nzcvec = mzconf(1)
             call rsp1gr(1, kzyva, idummy,0, isyma, 0, irefsy, etrs,&
                         cref, nzyvec, nzcvec, ovlap, isymdn, udv, fxo2k1k,&
                         xindx, mjwop, wrk(1), lfree, lorb, lcon, .true.)
             deallocate(fxo2k1k)
             if (atest) then
                e3test_value = ddot(kzyva,veca,1,etrs,1)
                write(lupri,*) 'PE-E3TEST, case 4 fxo2k1k ', e3test_value-e3test_old
                write(lupri,*) 'PE-E3TEST, total  ', e3test_value
             end if
          end if
          if (.not. tdhf) then
             ! fxc(1s2s) = R* ( <01L|..|02R> + <02L|..|01R> )Fe
             fcmo    = 0.0d0
             allocate(fxc1s2s(norbt,norbt))
             fxc1s2s = 0.0d0
             call uthu(0.125d0*fcaos(6*nnbasx+1:7*nnbasx), fcmo, cmo,&
                       wrk, nbast, norbt)
             call dsptsi(norbt, fcmo, fxc1s2s)
             !/ <0| [qj ,Fxc[1S,2S]] |0> \
             !| <j| Fxc[1S,2S] |0>       |
             !| <0| [qj+,Fxc[1S,2S]] |0> |
             !\ -<0| Fxc[1S,2S] |j>      /
             isymdn = 1
             ovlap  = 1.0d0
             isymst = muld2h(isyma, irefsy)
             if ( isymst .eq. irefsy ) then
                lcon = ( mzconf(isyma) .gt. 1 )
             else
                lcon = ( mzconf(isyma) .gt. 0 )
             end if
             lorb   = ( mzwopt(isyma) .gt. 0 )
             nzyvec = mzconf(1)
             nzcvec = mzconf(1)
             call rsp1gr(1 ,kzyva, idummy,0, isyma, 0, irefsy, etrs,&
                         cref, nzyvec, nzcvec, ovlap, isymdn, udv, fxc1s2s,&
                         xindx, mjwop, wrk(1), lfree, lorb, lcon, .true.)
             deallocate(fxc1s2s)
             if (atest) then
                e3test_value = ddot(kzyva,veca,1,etrs,1)
                write(lupri,*) 'PE-E3TEST, case 4 fxc1s2s', e3test_value-e3test_old
                write(lupri,*) 'PE-E3TEST, total  ', e3test_value
             end if
          end if
          ! fxc(1s2k) = ( <01L|[k2,Epq]|0> + <0|[k2,Epq]|01R> )Fe
          allocate(fxc1s2k(norbt,norbt))
          fcmo    = 0.0d0
          fxc1s2k = 0.0d0
          call uthu(2.0d0*fcaos(8*nnbasx+1:9*nnbasx), fcmo, cmo,&
                    wrk, nbast, norbt)
          call dsptsi(norbt, fcmo, fxc1s2k)
          !/ <0| [qj ,Fxc[2k,1S]] |0> \
          !| <j| Fxc[2k,1S] |0>       |
          !| <0| [qj+,Fxc[2k,1S]] |0> |
          !\ -<0| Fxc[2k,1S] |j>      /
          isymdn = 1
          ovlap  = 1.0d0
          isymst = muld2h(isyma, irefsy)
          if ( isymst .eq. irefsy ) then
             lcon = ( mzconf(isyma) .gt. 1 )
          else
             lcon = ( mzconf(isyma) .gt. 0 )
          end if
          lorb   = ( mzwopt(isyma) .gt. 0 )
          nzyvec = mzconf(1)
          nzcvec = mzconf(1)
          call rsp1gr(1 ,kzyva, idummy,0, isyma, 0, irefsy, etrs,&
                      cref, nzyvec, nzcvec, ovlap, isymdn, udv, fxc1s2k,&
                      xindx, mjwop, wrk(1), lfree, lorb, lcon, .true.)
          deallocate(fxc1s2k)
          if (atest) then
             e3test_value = ddot(kzyva,veca,1,etrs,1)
             write(lupri,*) 'PE-E3TEST, case 4 fxc1s2k ', e3test_value-e3test_old
             write(lupri,*) 'PE-E3TEST, total  ', e3test_value
          end if
          !fxc(1k2s) = ( <02L|[k1,Epq]|0> + <0|[k1,Epq]|02R> )Fe
          allocate(fxc1k2s(norbt,norbt))
          fcmo = 0.0d0
          fxc1k2s = 0.0d0
          call uthu(2.0d0*fcaos(9*nnbasx+1:10*nnbasx), fcmo, cmo,&
                    wrk, nbast, norbt)
          call dsptsi(norbt, fcmo, fxc1k2s)
          !/ <0| [qj ,Fxc[1k,2S]] |0> \
          !| <j| Fxc[1k,2S] |0>       |
          !| <0| [qj+,Fxc[1k,2S]] |0> |
          !\ -<0| Fxc[1S,2k] |j>      /
          isymdn = 1
          ovlap  = 1.0d0
          isymst = muld2h(isyma, irefsy)
          if ( isymst .eq. irefsy ) then
             lcon = ( mzconf(isyma) .gt. 1 )
          else
             lcon = ( mzconf(isyma) .gt. 0 )
          end if
          lorb   = ( mzwopt(isyma) .gt. 0 )
          nzyvec = mzconf(1)
          nzcvec = mzconf(1)
          call rsp1gr(1 ,kzyva, idummy,0, isyma, 0, irefsy, etrs,&
                      cref, nzyvec, nzcvec, ovlap, isymdn, udv, fxc1k2s,&
                      xindx, mjwop, wrk(1), lfree, lorb, lcon, .true.)
          deallocate(fxc1k2s)
          if (atest) then
             e3test_value = ddot(kzyva,veca,1,etrs,1)
             write(lupri,*) 'PE-E3TEST, case 4 fxc1k2s ', e3test_value-e3test_old
             write(lupri,*) 'PE-E3TEST, total  ', e3test_value
          end if
      end if ! pe_polar
      deallocate(cref)
      deallocate(fupe)

      call qexit('pe_rspmcqr')

end subroutine pelib_ifc_qrtest

subroutine pelib_ifc_rspmcqr(vecb, vecc, veca, atest, etrs, xindx, zymb, zymc,&
                             den1, udv, wrk, lfree, kzyva, kzyvb, kzyvc,&
                             isyma, isymb, isymc, cmo, mjwop)

      use pe_variables, only: pe_polar
      use polarizable_embedding, only: pe_master

      implicit none

#include "inforb.h"
#include "infvar.h"
#include "infdim.h"
#include "qrinf.h"
#include "priunit.h"
#include "dummy.h"
#include "inftap.h"
#include "infrsp.h"
#include "wrkrsp.h"

      integer :: kzyva, kzyvb, kzyvc
      integer :: isyma, isymb, isymc, isymbc
      integer :: lfree
      integer :: ilsym, irsym, ncl, ncr, kzvarl, kzvarr
      integer :: isymdn, isymst
      integer :: kcref, nzyvec, nzcvec
      integer :: iprone, nzconf, nzvar
      integer :: n2ash

      real*8 :: ovlap
      real*8 :: fact
      real*8 :: ddot
      real*8 :: e3test_value, e3test_old

      real*8, dimension(*) :: wrk
      real*8, dimension(1) :: tmpwrk
      real*8, dimension(*) :: cmo, xindx

      real*8, dimension(kzyva) :: etrs
      real*8, dimension(kzyvb) :: vecb
      real*8, dimension(kzyvc) :: vecc
      real*8, dimension(kzyva) :: veca

      real*8, dimension(norbt,norbt) :: zymb, zymc
      real*8, dimension(nashdi,nashdi) :: udv, den1
      real*8, dimension(nnashx) :: dv

      integer, dimension(2,maxwop,8) :: mjwop

      real*8, dimension(:), allocatable :: fpe
      real*8, dimension(:), allocatable :: cref
      real*8, dimension(:), allocatable :: dcaos, fcaos
      real*8, dimension(:), allocatable :: udtv, udtvao
      real*8, dimension(:), allocatable :: dvaao, dvbao, dvatr
      real*8, dimension(:), allocatable :: udcao, udvao
      real*8, dimension(:), allocatable :: udcmo, udvmo
      real*8, dimension(:), allocatable :: fcmo

      real*8, dimension(:,:), allocatable :: dva, dvb
      real*8, dimension(:,:), allocatable :: fupe
      real*8, dimension(:,:), allocatable :: fxpeb,fxpec, fx2pe
      real*8, dimension(:,:), allocatable :: fxo1k, fxc1s
      real*8, dimension(:,:), allocatable :: fxo2k, fxc2s
      real*8, dimension(:,:), allocatable :: fcas2_1, fcas2_2
      real*8, dimension(:,:), allocatable :: fcas3_1, fcas3_2
      real*8, dimension(:,:), allocatable :: fxo, fxo1k2k, fxo2k1k
      real*8, dimension(:,:), allocatable :: fxc1s2s
      real*8, dimension(:,:), allocatable :: fxc1s2k, fxc1k2s

      logical :: lexist, lopen, lcon, lorb
      logical :: fndlab
      logical :: atest

      lopen = .false.

      call qenter('pe_rspmcqr')

      call gtzymt(1, vecb, kzyvb, isymb, zymb, mjwop)
      call gtzymt(1, vecc, kzyvc, isymc, zymc, mjwop)

      if (atest) then
         e3test_old = ddot(kzyva,veca,1,etrs,1)
         write(lupri,*) 'PE-E3TEST, on entry', e3test_old
      end if

!      CALL GTZYMT(NSIM,VEC1,KZYV1,ISYMV1,ZYM1,MJWOP)
!      CALL GTZYMT(NSIM,VEC2,KZYV2,ISYMV2,ZYM2,MJWOP)

      !-----------------------------------------------------------
      ! Get Fg = Vmul - R*<0|F|>Fe from file
      !-----------------------------------------------------------
      if (.not. tdhf) then
         allocate(fpe(nnorbx))
         if (lusifc <= 0) then
             call gpopen(lusifc, 'SIRIFC', 'OLD', ' ', 'UNFORMATTED',&
                         idummy, .false.)
                 lopen = .true.
         end if
         rewind(lusifc)
         call mollab('PEFMAT  ', lusifc, lupri)
         call readt(lusifc, nnorbx, fpe)
         if (lopen) call gpclose(lusifc, 'KEEP')
         allocate(fupe(norbt,norbt))
         call dsptsi(norbt, fpe, fupe)
         deallocate(fpe)
      end if
      !-----------------------------------------------------------
      ! Density Factory ...
      !-----------------------------------------------------------

      allocate(cref(mzconf(1)))
      call getref(cref, mzconf(1))

      if (pe_polar) then
         allocate(udcao(n2basx))
         allocate(udvao(n2basx))
         if (.not. tdhf) then
            allocate(dcaos(10*nnbasx))
            else
            allocate(dcaos(4*nnbasx))
         end if
         dcaos = 0.0d0

         !  DTX = D_pq(k1) = <0|[k1,Epq]|0>
         allocate(udcmo(n2orbx),udvmo(n2orbx))
         udcmo = 0.0d0
         udvmo = 0.0d0
         call deq27mo(isymb, zymb, udv, udcmo, udvmo, wrk, lfree)
         if (nasht > 0) then
            udcmo = udcmo + udvmo
         end if
         udcao = 0.0d0
         call motoao(udcmo,udcao,cmo,isymb,wrk,lfree)
         call dgefsp(nbast, udcao, dcaos(1:nnbasx))
         ! needed to fit with HF code
         dcaos(1:nnbasx) = 0.5d0*dcaos(1:nnbasx)

         !  DT2X = D_pq(k2,k1) = <0|[k2,[k1,Epq]|0>
         udvmo = 0.0d0
         call oitd1(isymc,zymc,udcmo,udvmo,isymb)
         ! DT2X in udvmo (re-used to save memory)
         udcao = 0.0d0
         isymbc = muld2h(isymb,isymc)
         call motoao(udvmo,udcao,cmo,isymbc,wrk,lfree)
         call dgefsp(nbast, udcao, dcaos(nnbasx+1:2*nnbasx))
         ! needed to fit with HF code
         dcaos(nnbasx+1:2*nnbasx) = 0.5d0*dcaos(nnbasx+1:2*nnbasx)

         !  DTX = D_pq(k2) = <0|[k2,Epq]|0>
         udcmo = 0.0d0
         udvmo = 0.0d0
         call deq27mo(isymc, zymc, udv, udcmo, udvmo,wrk, lfree)
         if (nasht > 0) then
            udcmo = udcmo + udvmo
         end if
         udcao = 0.0d0
         call motoao(udcmo,udcao,cmo,isymb,wrk,lfree)
         call dgefsp(nbast, udcao, dcaos(2*nnbasx+1:3*nnbasx))
         ! needed to fit with HF code
         dcaos(2*nnbasx+1:3*nnbasx) = 0.5d0*dcaos(2*nnbasx+1:3*nnbasx)

         !  DT2X = D_pq(k1,k2) = <0|[k1,[k2,Epq]|0>
         udvmo = 0.0d0
         call oitd1(isymb,zymb,udcmo,udvmo,isymc)
         ! DT2X in udvmo (re-used to save memory)
         udcao = 0.0d0
         isymbc = muld2h(isymc,isymb)
         call motoao(udvmo,udcao,cmo,isymbc,wrk,lfree)
         call dgefsp(nbast, udcao, dcaos(3*nnbasx+1:4*nnbasx))
         ! needed to fit with HF code
         dcaos(3*nnbasx+1:4*nnbasx) = 0.5d0*dcaos(3*nnbasx+1:4*nnbasx)
         deallocate(udcmo,udvmo)

         if (tdhf) then
            write(lupri,*) 'PE-DFT or HF QR detected: Skipping CI dens.'
         end if
         if (.not. tdhf ) then
            write(lupri,*) 'PE-MCSCF QR detected: Constructing CI dens.'

            ! Construct the density matrix <02L|..|0> + <0|..|02R>
            ilsym  = irefsy
            irsym  = muld2h(irefsy,isymc)
            ncl    = mzconf(1)
            ncr    = mzconf(isymc)
            kzvarl = mzconf(1)
            kzvarr = mzyvar(isymc)

            den1 = 0.0d0 ! edh: This is equal to udtv later...
            allocate(udtv(n2ashx), udtvao(n2basx))

            udtv = 0.0d0
            udtvao = 0.0d0
            call rspgdm(1, ilsym, irsym, ncl, ncr, kzvarl, kzvarr,&
                        cref, vecc, ovlap, udtv, dummy, 0 ,0, .true.,&
                        .true., xindx, wrk, 1, lfree, .true.)
            call fckden2(.false.,.true., dummy, udtvao, cmo,&
                         udtv, wrk, lfree)
            call dgefsp(nbast, udtvao, dcaos(4*nnbasx+1:5*nnbasx))
            dcaos(4*nnbasx+1:5*nnbasx) = 1.0d0*dcaos(4*nnbasx+1:5*nnbasx)

            ! Construct the density matrix <01L|..|0> + <0|..|01R>
            ilsym  = irefsy
            irsym  = muld2h(irefsy,isymb)
            ncl    = mzconf(1)
            ncr    = mzconf(isymb)
            kzvarl = mzconf(1)
            kzvarr = mzyvar(isymb)

            udtv = 0.0d0
            udtvao = 0.0d0
            call rspgdm(1, ilsym, irsym, ncl, ncr, kzvarl, kzvarr,&
                        cref, vecb, ovlap, udtv, dummy, 0 ,0, .true.,&
                        .true., xindx, wrk, 1, lfree, .true.)
            call fckden2(.false.,.true., dummy, udtvao, cmo,&
                         udtv, wrk, lfree)
            call dgefsp(nbast, udtvao, dcaos(5*nnbasx+1:6*nnbasx))
            dcaos(5*nnbasx+1:6*nnbasx) = 1.0d0*dcaos(5*nnbasx+1:6*nnbasx)

            if (mzconf(isymb) .gt. 0 .and. mzconf(isymc) .gt. 0) then

               ! Construct <01L|..|02R> + <02L|..|01R> density
               ilsym  = muld2h(irefsy,isymb)
               irsym  = muld2h(irefsy,isymc)
               ncl    = mzconf(isymb)
               ncr    = mzconf(isymc)
               kzvarl = mzyvar(isymb)
               kzvarr = mzyvar(isymc)
               isymdn = muld2h(ilsym,irsym)

               udtv = 0.0d0
               udtvao = 0.0d0
               call rspgdm(1, ilsym, irsym, ncl, ncr, kzvarl, kzvarr,&
                           vecb, vecc, ovlap, udtv, dummy, 0 ,0,&
                          .true., .true., xindx, wrk, 1, lfree,&
                          .false.)
               call fckden2(.false.,.true., dummy, udtvao, cmo,&
                            udtv, wrk, lfree)
               call dgefsp(nbast, udtvao, dcaos(6*nnbasx+1:7*nnbasx))
            end if
            dcaos(6*nnbasx+1:7*nnbasx) = 1.0d0*dcaos(6*nnbasx+1:7*nnbasx)

            ! D_pq = <0|Epq|0>
            udcao = 0.0d0
            udvao = 0.0d0
            call dgetsp(nasht, udv, dv)
            call fckden((nisht>0), (nasht>0), udcao, udvao,&
                        cmo, dv, wrk, lfree)
            if (nisht==0) udcao = 0.0d0
            udcao = udcao + udvao
            call dgefsp(nbast, udcao, dcaos(7*nnbasx+1:8*nnbasx))
            dcaos(7*nnbasx+1:8*nnbasx) = 0.5d0*dcaos(7*nnbasx+1:8*nnbasx)

            ! D_pq(S1,k2) = <01L|[k2,Epq]|0> + <0|[k2,Epq]|01R>
            allocate(dva(norbt,nasht), dvb(norbt,nasht))
            allocate(dvaao(n2basx), dvbao(n2basx), dvatr(n2basx))
            dva    = 0.0d0
            dvb    = 0.0d0
            dvaao  = 0.0d0
            dvbao  = 0.0d0
            dvatr  = 0.0d0
            udtvao = 0.0d0
            call rsptr1(1, udv, zymb, dva, dvb)
            call fckden2(.false.,.true., dummy, dvaao, cmo,&
                         dva, wrk, lfree)
            call fckden2(.false.,.true., dummy, dvbao, cmo,&
                         dvb, wrk, lfree)
            call mtrsp(nbast, nbast, dvaao, nbast, dvatr, nbast)
            udtvao = dvbao - dvatr
            call dgefsp(nbast, udtvao, dcaos(8*nnbasx+1:9*nnbasx))
            dcaos(8*nnbasx+1:9*nnbasx) = 2.0d0*dcaos(8*nnbasx+1:9*nnbasx)

            ! D_pq(k1,S2) = <02L|[k1,Epq]|0> + <0|[k1,Epq]|02R>
            dva    = 0.0d0
            dvb    = 0.0d0
            dvaao  = 0.0d0
            dvbao  = 0.0d0
            dvatr  = 0.0d0
            udtvao = 0.0d0
            call rsptr1(1, udv, zymc, dva, dvb)
            call fckden2(.false.,.true., dummy, dvaao, cmo,&
                         dva, wrk, lfree)
            call fckden2(.false.,.true., dummy, dvbao, cmo,&
                         dvb, wrk, lfree)
            call mtrsp(nbast, nbast, dvaao, nbast, dvatr, nbast)
            udtvao = dvbao - dvatr
            call dgefsp(nbast, udtvao, dcaos(9*nnbasx+1:10*nnbasx))
            dcaos(9*nnbasx+1:10*nnbasx) = 2.0d0*dcaos(9*nnbasx+1:10*nnbasx)

            deallocate(dva, dvb, dvaao, dvbao, dvatr)
            deallocate(udtv,udtvao)

         end if
         deallocate(udcao, udvao)
         !-----------------------------------------------------------
         ! Calculate PE response operators in AO basis
         !-----------------------------------------------------------

         if (.not. tdhf) then
            allocate(fcaos(10*nnbasx))
            fcaos = 0.0d0
#if defined(VAR_MPI)
            call pelib_ifc_start_slaves(3)
#endif
            call pe_master(runtype='dynamic_response', &
                           triang=.true., &
                          ndim=nbast, &
                          nmats=10, &
                          denmats=dcaos, &
                          fckmats=fcaos)
         else
            allocate(fcaos(4*nnbasx))
            fcaos = 0.0d0
#if defined(VAR_MPI)
            call pelib_ifc_start_slaves(3)
#endif
            call pe_master(runtype='dynamic_response', &
                           triang=.true., &
                           ndim=nbast, &
                           nmats=4, &
                           denmats=dcaos, &
                           fckmats=fcaos)
         end if
         deallocate(dcaos)

      end if ! pe_polar

      if ( .not. tdhf ) then
         !-----------------------------------------------------------
         !case 1
         !-----------------------------------------------------------
         if ( mzconf(isymb) .eq. 0 .or. mzconf(isymc) .eq. 0 ) return

         !/   <01L| [qj,TB] |02R>  + <02L| [qj,TB] |01R>  \
         !|                       0                       |
         !|   <01L| [qj+,TB] |02R> + <02L| [qj+,TB] |01R> |
         !\                       0                       /

         ! ionstruct <01L|..|02R> + <02L|..|01R> density
         ilsym  = muld2h(irefsy,isymb)
         irsym  = muld2h(irefsy,isymc)
         ncl    = mzconf(isymb)
         ncr    = mzconf(isymc)
         kzvarl = mzyvar(isymb)
         kzvarr = mzyvar(isymc)

         den1 = 0.0d0
         call rspgdm(1, ilsym, irsym, ncl, ncr, kzvarl, kzvarr,&
                     vecb, vecc, ovlap, den1, dummy, 0, 0, .true.,&
                     .true., xindx, wrk, 1, lfree, .false.)

         ! Make the gradient
         isymdn = muld2h(ilsym,irsym)

         if ( mzwopt(isyma) .gt. 0 ) then
            call orbsx(1, isyma, kzyva, etrs, fupe, ovlap,&
                       isymdn, den1, mjwop, 1, lfree)
            if (atest) then
               e3test_value = ddot(kzyva,veca,1,etrs,1)
               write(lupri,*) 'PE-E3TEST, case 1 ', e3test_value-e3test_old
               e3test_old = e3test_value
            end if
        end if
      end if
      !-----------------------------------------------------------
      !case 2
      !-----------------------------------------------------------

      if (pe_polar) then

         allocate(fcmo(nnorbx))
         allocate(fxo1k(norbt,norbt))
         fxo1k = 0.0d0
         if (.not. tdhf) then
            allocate(fxc1s(norbt,norbt))
            fxc1s = 0.0d0
         end if

         ! Fxo = R*<0|Fe(1k)|0>Fe
         fcmo = 0.0d0
         call uthu(2.0d0*fcaos(1:nnbasx), fcmo, cmo,&
                   wrk, nbast, norbt)
         call dsptsi(norbt, fcmo, fxo1k)

         if (.not. tdhf) then
            ! Fxc(1S) = ( R*<01lE|0>+<0|E01R> )Fe
            ! edh: Should it be 1.0d0 or 2.0d0 ???
            fcmo = 0.0d0
            call uthu(1.0d0*fcaos(4*nnbasx+1:5*nnbasx), fcmo, cmo,&
                      wrk, nbast, norbt)
            call dsptsi(norbt, fcmo, fxc1s)
         end if

         ! fcas2_1 = Fa[1](1k)
         allocate(fcas2_1(norbt,norbt))
         fcas2_1 = 0.0d0
         if (.not. tdhf) then
            fcas2_1 = fxo1k + fxc1s
            deallocate(fxc1s)
         else
            fcas2_1 = fxo1k
         end if
         deallocate(fxo1k)

         if (.not. tdhf) then
            if (mzconf(isymc) .le. 0) return
            !/   0    \
            !| Sj(2)  | * <0| Fa[1](1k) |0>
            !|   0    |
            !\ Sj(2)* /
            ! edh should fact not be zeroed? here and below for fcas2_2?
            if (isyma .eq. isymc) then
                ovlap = 1.0d0
                call melone(fcas2_1, 1, udv, ovlap, fact,&
                            200, 'fact for Fxo(1k) + Fxc(1S) ')
                nzconf = mzconf(isyma)
                nzvar  = mzvar(isyma)
                call daxpy(nzconf, fact, vecc, 1, etrs, 1)
                call daxpy(nzconf, fact, vecc(nzvar+1), 1, etrs(nzvar+1), 1)
               if (atest) then
                  e3test_value = ddot(kzyva,veca,1,etrs,1)
                  write(lupri,*) 'PE-E3TEST, case 2a', e3test_value-e3test_old
                  e3test_old = e3test_value
               end if
            end if
         end if

         allocate(fxo2k(norbt,norbt))
         fxo2k = 0.0d0
         if (.not. tdhf) then
            allocate(fxc2s(norbt,norbt))
            fxc2s = 0.0d0
         end if

         ! Fxo(2k) = R*<0|[2k,Epq]|0>Fe
         fcmo = 0.0d0
         call uthu(2.0d0*fcaos(2*nnbasx+1:3*nnbasx),&
                   fcmo, cmo, wrk, nbast, norbt)
         call dsptsi(norbt, fcmo, fxo2k)

         if (.not. tdhf) then
            ! Fxc(1S) = ( R*<01lE|0>+<0|E01R> )Fe
            ! edh: Should it be 1.0d0 or 2.0d0 ???
            fcmo = 0.0d0
            call uthu(1.0d0*fcaos(5*nnbasx+1:6*nnbasx), fcmo, cmo,&
                      wrk, nbast, norbt)
             call dsptsi(norbt, fcmo, fxc2s)
         end if

         ! fcas2_2 = Fa[1](2k)
         allocate(fcas2_2(norbt,norbt))
         fcas2_2 = 0.0d0
         if (.not. tdhf) then
            fcas2_2 = fxo2k + fxc2s
            deallocate(fxc2s)
         else
            fcas2_2 = fxo2k
         end if
         deallocate(fxo2k)

         if (.not. tdhf) then
            if (mzconf(isymb) .le. 0) return
            !/   0    \
            !| Sj(1)  | * <0| Fa[1](2k) |0>
            !|   0    |
            !\ Sj(1)* /
            if (isyma .eq. isymb) then
                ovlap = 1.0d0
                ! edh should fact not be zeroed??
                call melone(fcas2_2, 1, udv, ovlap, fact,&
                            200, 'fact for Fxo(2k) + Fxc(2S) ')
                  nzconf = mzconf(isyma)
                  nzvar  = mzvar(isyma)
                  call daxpy(nzconf, fact, vecb, 1, etrs, 1)
                  call daxpy(nzconf, fact, vecb(nzvar+1), 1, etrs(nzvar+1), 1)
               if (atest) then
                  e3test_value = ddot(kzyva,veca,1,etrs,1)
                  write(lupri,*) 'PE-E3TEST, case 2b', e3test_value-e3test_old
                  e3test_old = e3test_value
               end if
            end if
         end if
      end if ! pe_polar
!     !-----------------------------------------------------------
!     !case 3
!     !-----------------------------------------------------------
      allocate(fxpeb(norbt,norbt))
      fxpeb = 0.0d0
      call oith1(isymb,zymb,fupe,fxpeb,1)

      ! fcas3_1 = Fg(k1) + F[a1]
      allocate(fcas3_1(norbt,norbt))
      if (pe_polar) then
          fcas3_1 = fxpeb + fcas2_1
      else
          fcas3_1 = fxpeb
      end if

      if (.not. tdhf) then

         if (mzconf(isymc) .le. 0) return

        !/   <0| [qj,TD1] |02R>  + <02L| [qj,TD1] |0>  \
        !|   <j| TD1 |02R>                             |
        !|   <0| [qj+,TD1] |02R> + <02L| [qj+,TD1] |0> |
        !\  -<02L| TD1 |j>                             /

        ! 1a. Construct the density matrix <02L|..|0> + <0|..|02R>
         ilsym  = irefsy
         irsym  = muld2h(irefsy,isymc)
         ncl    = mzconf(1)
         ncr    = mzconf(isymc)
         kzvarl = mzconf(1)
         kzvarr = mzyvar(isymc)

         den1 = 0.0d0
         call rspgdm(1, ilsym, irsym, ncl, ncr, kzvarl, kzvarr,&
                     cref, vecc, ovlap, den1, dummy, 0 ,0, .true.,&
                     .true., xindx, wrk, 1, lfree, .true.)

!        1b. Make the gradient
         isymdn = muld2h(ilsym,irsym)
         isymst = muld2h(isyma,irefsy)
         if ( isymst .eq. irefsy ) then
            lcon = ( mzconf(isyma) .gt. 1 )
         else
            lcon = ( mzconf(isyma) .gt. 0 )
         end if
         lorb    = ( mzwopt(isyma) .gt. 0 )
         nzyvec = mzyvar(isymc)
         nzcvec = mzconf(isymc)

         call rsp1gr(1, kzyva, idummy, 0 , isyma, 0, isymc, etrs,&
                     vecc, nzyvec, nzcvec, ovlap, isymdn, den1, fcas3_1,&
                     xindx, mjwop, wrk(1), lfree, lorb, lcon, .false.)
         if (atest) then
            e3test_value = ddot(kzyva,veca,1,etrs,1)
            write(lupri,*) 'PE-E3TEST, case 3a', e3test_value-e3test_old
            e3test_old = e3test_value
         end if
      end if
      deallocate(fcas3_1)

      allocate(fxpec(norbt,norbt))
      fxpec = 0.0d0
      call oith1(isymc,zymc,fupe,fxpec,1)

      ! fcas3_2 = Fg(2k) + Fa[1]
      allocate(fcas3_2(norbt,norbt))
      if (pe_polar) then
          fcas3_2 = fxpec + fcas2_2
      else
          fcas3_2 = fxpec
      end if

      if (.not. tdhf) then

         if (mzconf(isymb) .le. 0) return

        !/   <0| [qj,TD2] |01R>  + <01L| [qj,TD2] |0>  \
        !|   <j| TD2 |01R>                             |
        !|   <0| [qj+,TD2] |01R> + <01L| [qj+,TD2] |0> |
        !\  -<01L| TD2 |j>                             /

        ! 2a. Construct the density matrix <01L|..|0> + <0|..|01R>
         ilsym  = irefsy
         irsym  = muld2h(irefsy,isymb)
         ncl    = mzconf(1)
         ncr    = mzconf(isymb)
         kzvarl = mzconf(1)
         kzvarr = mzyvar(isymb)

         den1 = 0.0d0
         call rspgdm(1, ilsym, irsym, ncl, ncr, kzvarl, kzvarr,&
                     cref, vecb, ovlap, den1, dummy, 0 ,0, .true.,&
                     .true., xindx, wrk, 1, lfree, .true.)

        ! 2b. Make the gradient
         isymdn = muld2h(ilsym,irsym)
         isymst = muld2h(isyma,irefsy)
         if ( isymst .eq. irefsy ) then
            lcon = ( mzconf(isyma) .gt. 1 )
         else
            lcon = ( mzconf(isyma) .gt. 0 )
         end if
         lorb    = ( mzwopt(isyma) .gt. 0 )
         nzyvec = mzyvar(isymb)
         nzcvec = mzconf(isymb)

         call rsp1gr(1, kzyva, idummy, 0 , isyma, 0, isymb, etrs,&
                     vecb, nzyvec, nzcvec, ovlap, isymdn, den1, fcas3_2,&
                     xindx, mjwop, wrk(1), lfree, lorb, lcon, .false.)
            if (atest) then
               e3test_value = ddot(kzyva,veca,1,etrs,1)
               write(lupri,*) 'PE-E3TEST, case 3b', e3test_value-e3test_old
               e3test_old = e3test_value
            end if
      end if
      deallocate(fcas3_2)
!     !-----------------------------------------------------------
!     !case 4
!     !-----------------------------------------------------------

!     !fx2pe = 0.5*Fg(1k,2k) + 0.5*Fg(2k,1k) + ...
      allocate(fx2pe(norbt,norbt))
      fx2pe = 0.0d0
      if (.not. tdhf) then
         call oith1(isymc, zymc, fxpeb, fx2pe, isymb)
         call oith1(isymb, zymb, fxpec, fx2pe, isymc)
         deallocate(fxpeb,fxpec)
      end if
      fx2pe = 0.5d0 * fx2pe

      ! ... + fcas2_1(2k) + fcas2_2(k1)
      if (pe_polar) then
          call oith1(isymc,zymc,fcas2_1, fx2pe,  isymb)
          call oith1(isymb,zymb,fcas2_2, fx2pe, isymc)
          deallocate(fcas2_1,fcas2_2)

          if (.not. tdhf ) then
          ! + ( S(1)S*(2) + S(2)S*(1) ) * fxo  + ...
             if ((isymb .eq. isymc) .and. (mzconf(isymb) .gt. 0)) then
                allocate(fxo(norbt,norbt))
                fcmo = 0.0d0
                fxo = 0.0d0
               ! fact = 0.0d0??
               ! edh: was changed from 0.25 to 1.0 (check this factor!)
                call uthu(1.0d0*fcaos(7*nnbasx+1:8*nnbasx), fcmo, cmo,&
                          wrk, nbast, norbt)
                call dsptsi(norbt, fcmo, fxo)
                nzconf = mzconf(isymb)
                nzvar  = mzvar(isymb)
                fact   = ddot(nzconf, vecb, 1, vecc(nzvar+1), 1) + &
                         ddot(nzconf, vecc, 1, vecb(nzvar+1), 1)
                call daxpy(n2orbx, fact, fxo, 1, fx2pe, 1)
                deallocate(fxo)
             end if
          end if

          allocate(fxo1k2k(norbt,norbt), fxo2k1k(norbt,norbt))
          fxo1k2k = 0.0d0
          fxo2k1k = 0.0d0
          if (.not. tdhf) then
             allocate(fxc1s2s(norbt,norbt))
             allocate(fxc1k2s(norbt,norbt), fxc1s2k(norbt,norbt))
             fxc1s2s = 0.0d0
             fxc1s2k = 0.0d0
             fxc1k2s = 0.0d0
          end if

          ! ... + fxo(1k,2k) + fxo(2k,1k)
          if (mzwopt(isymb).gt.0 .and. mzwopt(isymc).gt.0) then
             !fxo(1k,2k) = <0|Fe(1k,2k)|0>Fe
             fcmo = 0.0d0
             call uthu(1.0d0*fcaos(nnbasx+1:2*nnbasx),&
                       fcmo, cmo, wrk, nbast, norbt)
             call dsptsi(norbt, fcmo, fxo1k2k)

             !fxo(2k,1k) = R*<0|Fe(2k,1k)|0>Fe
             fcmo    = 0.0d0
             call uthu(1.0d0*fcaos(3*nnbasx+1:4*nnbasx),&
                       fcmo, cmo, wrk, nbast, norbt)
             call dsptsi(norbt, fcmo, fxo2k1k)
          end if

          ! + ... fxc(1s2s) + fxc(2s1s) + 2fxc(1s2k) + 2fxc(1k2S)
          if (.not. tdhf) then
             ! fxc(1s2s) = R* ( <01L|..|02R> + <02L|..|01R> )Fe
             fcmo    = 0.0d0
             call uthu(0.125d0*fcaos(6*nnbasx+1:7*nnbasx), fcmo, cmo,&
                       wrk, nbast, norbt)
             call dsptsi(norbt, fcmo, fxc1s2s)
             ! edh should fcmo be re-zeroed or is it done in UTHU??
             ! + ... fxc(1s2k) + fxc(1k2S)
             ! fxc(1s2k) = ( <01L|[k2,Epq]|0> + <0|[k2,Epq]|01R> )Fe
             call uthu(2.0d0*fcaos(8*nnbasx+1:9*nnbasx), fcmo, cmo,&
                       wrk, nbast, norbt)
             call dsptsi(norbt, fcmo, fxc1s2k)

             ! fxc(1k2s) = ( <02L|[k1,Epq]|0> + <0|[k1,Epq]|02R> )Fe
             call uthu(2.0d0*fcaos(9*nnbasx+1:10*nnbasx), fcmo, cmo,&
                       wrk, nbast, norbt)
             call dsptsi(norbt, fcmo, fxc1k2s)
          end if
          deallocate(fcmo)

          if (.not. tdhf) then
             fx2pe = fx2pe + fxo1k2k + fxo2k1k &
                   + fxc1s2s + fxc1s2k + fxc1k2s
             deallocate(fxc1s2s, fxc1s2k, fxc1k2s)
          else
             fx2pe = fx2pe + fxo1k2k + fxo2k1k
          end if
          deallocate(fxo1k2k, fxo2k1k)
      end if

      !/ <0| [qj ,TE] |0> \
      !| <j| TE |0>       |
      !| <0| [qj+,TE] |0> |
      !\ -<0| TE |j>      /

      isymdn = 1
      ovlap  = 1.0d0
      isymst = muld2h(isyma, irefsy)
      if ( isymst .eq. irefsy ) then
         lcon = ( mzconf(isyma) .gt. 1 )
      else
         lcon = ( mzconf(isyma) .gt. 0 )
      end if
      lorb   = ( mzwopt(isyma) .gt. 0 )
      nzyvec = mzconf(1)
      nzcvec = mzconf(1)

      call rsp1gr(1 ,kzyva, idummy,0, isyma, 0, irefsy, etrs,&
                  cref, nzyvec, nzcvec, ovlap, isymdn, udv, fx2pe,&
                  xindx, mjwop, wrk(1), lfree, lorb, lcon, .true.)
      if (atest) then
         e3test_value = ddot(kzyva,veca,1,etrs,1)
         write(lupri,*) 'PE-E3TEST, case 4 ', e3test_value-e3test_old
         write(lupri,*) 'PE-E3TEST, total  ', e3test_value
      end if
      deallocate(cref)
      deallocate(fupe)

      call qexit('pe_rspmcqr')

end subroutine pelib_ifc_rspmcqr

subroutine pelib_ifc_cro(vecb, vecc, vecd, etrs, xindx, zymb, zymc, zymd, udv,&
                         wrk, nwrk, kzyva, kzyvb, kzyvc, kzyvd, isyma, isymb,&
                         isymc, isymd, cmo,mjwop)
    implicit none
#include "inforb.h"
#include "infvar.h"
#include "infrsp.h"
#include "infdim.h"
#include "qrinf.h"

    integer :: kzyva, kzyvb, kzyvc, kzyvd
    integer :: isyma, isymb, isymc, isymd
    integer :: nwrk
    real*8, dimension(nwrk) :: wrk
    real*8, dimension(kzyva) :: etrs
    real*8, dimension(kzyvb) :: vecb
    real*8, dimension(kzyvc) :: vecc
    real*8, dimension(kzyvd) :: vecd
    real*8, dimension(ncmot) :: cmo
    real*8, dimension(norbt,norbt) :: zymb, zymc, zymd
    real*8, dimension(nashdi,nashdi) :: udv
    real*8, dimension(lcindx) :: xindx
    integer, dimension(2,maxwop,8) :: mjwop

    integer :: i, j, k
    integer :: idum = 1
    real*8, dimension(:), allocatable :: udcao, ufcmo
    real*8, dimension(:), allocatable :: dcaos, fcaos
    real*8, dimension(:), allocatable :: fcmo

    call qenter('pelib_ifc_cro')
    if (.not. use_pelib()) call quit('PElib not active')

    call gtzymt(1, vecb, kzyvb, isymb, zymb, mjwop)
    call gtzymt(1, vecc, kzyvc, isymc, zymc, mjwop)
    call gtzymt(1, vecd, kzyvd, isymd, zymd, mjwop)

    allocate(udcao(n2basx))
    allocate(ufcmo(n2orbx))
    allocate(dcaos(15*nnbasx))
    dcaos = 0.0d0

    udcao = 0.0d0
    call cdens1(isymb, cmo, zymb, udcao, wrk, nwrk)
    call dgefsp(nbast, udcao, dcaos(1:nnbasx))
    udcao = 0.0d0
    call cdens1(isymc, cmo, zymc, udcao, wrk, nwrk)
    call dgefsp(nbast, udcao, dcaos(nnbasx+1:2*nnbasx))
    udcao = 0.0d0
    call cdens1(isymd, cmo, zymd, udcao, wrk, nwrk)
    call dgefsp(nbast, udcao, dcaos(2*nnbasx+1:3*nnbasx))

    udcao = 0.0d0
    call cdens2(isymb, isymc, cmo, zymb, zymc, udcao,&
                wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
    call dgefsp(nbast, udcao, dcaos(3*nnbasx+1:4*nnbasx))
    udcao = 0.0d0
    call cdens2(isymc, isymb, cmo, zymc, zymb, udcao,&
                wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
    call dgefsp(nbast, udcao, dcaos(4*nnbasx+1:5*nnbasx))
    udcao = 0.0d0
    call cdens2(isymb, isymd, cmo, zymb, zymd, udcao,&
                wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
    call dgefsp(nbast, udcao, dcaos(5*nnbasx+1:6*nnbasx))
    udcao = 0.0d0
    call cdens2(isymd, isymb, cmo, zymd, zymb, udcao,&
                wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
    call dgefsp(nbast, udcao, dcaos(6*nnbasx+1:7*nnbasx))
    udcao = 0.0d0
    call cdens2(isymc, isymd, cmo, zymc, zymd, udcao,&
                wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
    call dgefsp(nbast, udcao, dcaos(7*nnbasx+1:8*nnbasx))
    udcao = 0.0d0
    call cdens2(isymd, isymc, cmo, zymd, zymc, udcao,&
                wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
    call dgefsp(nbast, udcao, dcaos(8*nnbasx+1:9*nnbasx))

    udcao = 0.0d0
    call cdens3(isymb, isymc, isymd, cmo, zymb, zymc, zymd, udcao,&
                wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
    call dgefsp(nbast, udcao, dcaos(9*nnbasx+1:10*nnbasx))
    udcao = 0.0d0
    call cdens3(isymd, isymb, isymc, cmo, zymd, zymb, zymc, udcao,&
                wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
    call dgefsp(nbast, udcao, dcaos(10*nnbasx+1:11*nnbasx))
    udcao = 0.0d0
    call cdens3(isymc, isymd, isymb, cmo, zymc, zymd, zymb, udcao,&
                wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
    call dgefsp(nbast, udcao, dcaos(11*nnbasx+1:12*nnbasx))
    udcao = 0.0d0
    call cdens3(isymb, isymd, isymc, cmo, zymb, zymd, zymc, udcao,&
                wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
    call dgefsp(nbast, udcao, dcaos(12*nnbasx+1:13*nnbasx))
    udcao = 0.0d0
    call cdens3(isymc, isymb, isymd, cmo, zymc, zymb, zymd, udcao,&
                wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
    call dgefsp(nbast, udcao, dcaos(13*nnbasx+1:14*nnbasx))
    udcao = 0.0d0
    call cdens3(isymd, isymc, isymb, cmo, zymd, zymc, zymb, udcao,&
                wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
    call dgefsp(nbast, udcao, dcaos(14*nnbasx+1:15*nnbasx))

    deallocate(udcao)

    allocate(fcaos(15*nnbasx))
    call pelib_ifc_response(dcaos, fcaos, 15)
    deallocate(dcaos)

    allocate(fcmo(nnorbx))
    ufcmo = 0.0d0

    i = 1
    j = nnbasx
    call uthu(1.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
    wrk(1:2*n2orbx) = 0.0d0
    call dsptsi(norbt, fcmo, wrk(1:n2orbx))
    call oith1(isymc, zymc, wrk(1:n2orbx), wrk(n2orbx+1:2*n2orbx), isyma)
    call oith1(isymd, zymd, wrk(n2orbx+1:2*n2orbx), ufcmo, isyma)
    i = i + nnbasx
    j = j + nnbasx
    call uthu(1.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
    wrk(1:2*n2orbx) = 0.0d0
    call dsptsi(norbt, fcmo, wrk(1:n2orbx))
    call oith1(isymb, zymb, wrk(1:n2orbx), wrk(n2orbx+1:2*n2orbx), isyma)
    call oith1(isymd, zymd, wrk(n2orbx+1:2*n2orbx), ufcmo, isyma)
    i = i + nnbasx
    j = j + nnbasx
    call uthu(1.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
    wrk(1:2*n2orbx) = 0.0d0
    call dsptsi(norbt, fcmo, wrk(1:n2orbx))
    call oith1(isymb, zymb, wrk(1:n2orbx), wrk(n2orbx+1:2*n2orbx), isyma)
    call oith1(isymc, zymc, wrk(n2orbx+1:2*n2orbx), ufcmo, isyma)
    i = 1
    j = nnbasx
    call uthu(1.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
    wrk(1:2*n2orbx) = 0.0d0
    call dsptsi(norbt, fcmo, wrk(1:n2orbx))
    call oith1(isymd, zymd, wrk(1:n2orbx), wrk(n2orbx+1:2*n2orbx), isyma)
    call oith1(isymc, zymc, wrk(n2orbx+1:2*n2orbx), ufcmo, isyma)
    i = i + nnbasx
    j = j + nnbasx
    call uthu(1.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
    wrk(1:2*n2orbx) = 0.0d0
    call dsptsi(norbt, fcmo, wrk(1:n2orbx))
    call oith1(isymd, zymd, wrk(1:n2orbx), wrk(n2orbx+1:2*n2orbx), isyma)
    call oith1(isymb, zymb, wrk(n2orbx+1:2*n2orbx), ufcmo, isyma)
    i = i + nnbasx
    j = j + nnbasx
    call uthu(1.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
    wrk(1:2*n2orbx) = 0.0d0
    call dsptsi(norbt, fcmo, wrk(1:n2orbx))
    call oith1(isymc, zymc, wrk(1:n2orbx), wrk(n2orbx+1:2*n2orbx), isyma)
    call oith1(isymb, zymb, wrk(n2orbx+1:2*n2orbx), ufcmo, isyma)

    i = i + nnbasx
    j = j + nnbasx
    call uthu(1.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
    wrk(1:n2orbx) = 0.0d0
    call dsptsi(norbt, fcmo, wrk(1:n2orbx))
    call oith1(isymd, zymd, wrk(1:n2orbx), ufcmo, isyma)
    i = i + nnbasx
    j = j + nnbasx
    call uthu(1.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
    wrk(1:n2orbx) = 0.0d0
    call dsptsi(norbt, fcmo, wrk(1:n2orbx))
    call oith1(isymd, zymd, wrk(1:n2orbx), ufcmo, isyma)
    i = i + nnbasx
    j = j + nnbasx
    call uthu(1.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
    wrk(1:n2orbx) = 0.0d0
    call dsptsi(norbt, fcmo, wrk(1:n2orbx))
    call oith1(isymc, zymc, wrk(1:n2orbx), ufcmo, isyma)
    i = i + nnbasx
    j = j + nnbasx
    call uthu(1.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
    wrk(1:n2orbx) = 0.0d0
    call dsptsi(norbt, fcmo, wrk(1:n2orbx))
    call oith1(isymc, zymc, wrk(1:n2orbx), ufcmo, isyma)
    i = i + nnbasx
    j = j + nnbasx
    call uthu(1.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
    wrk(1:n2orbx) = 0.0d0
    call dsptsi(norbt, fcmo, wrk(1:n2orbx))
    call oith1(isymb, zymb, wrk(1:n2orbx), ufcmo, isyma)
    i = i + nnbasx
    j = j + nnbasx
    call uthu(1.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
    wrk(1:n2orbx) = 0.0d0
    call dsptsi(norbt, fcmo, wrk(1:n2orbx))
    call oith1(isymb, zymb, wrk(1:n2orbx), ufcmo, isyma)

    i = i + nnbasx
    j = j + nnbasx
    call uthu(1.0d0/3.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
    wrk(1:n2orbx) = 0.0d0
    call dsptsi(norbt, fcmo, wrk(1:n2orbx))
    ufcmo = ufcmo + wrk(1:n2orbx)
    i = i + nnbasx
    j = j + nnbasx
    call uthu(1.0d0/3.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
    wrk(1:n2orbx) = 0.0d0
    call dsptsi(norbt, fcmo, wrk(1:n2orbx))
    ufcmo = ufcmo + wrk(1:n2orbx)
    i = i + nnbasx
    j = j + nnbasx
    call uthu(1.0d0/3.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
    wrk(1:n2orbx) = 0.0d0
    call dsptsi(norbt, fcmo, wrk(1:n2orbx))
    ufcmo = ufcmo + wrk(1:n2orbx)
    i = i + nnbasx
    j = j + nnbasx
    call uthu(1.0d0/3.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
    wrk(1:n2orbx) = 0.0d0
    call dsptsi(norbt, fcmo, wrk(1:n2orbx))
    ufcmo = ufcmo + wrk(1:n2orbx)
    i = i + nnbasx
    j = j + nnbasx
    call uthu(1.0d0/3.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
    wrk(1:n2orbx) = 0.0d0
    call dsptsi(norbt, fcmo, wrk(1:n2orbx))
    ufcmo = ufcmo + wrk(1:n2orbx)
    i = i + nnbasx
    j = j + nnbasx
    call uthu(1.0d0/3.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
    wrk(1:n2orbx) = 0.0d0
    call dsptsi(norbt, fcmo, wrk(1:n2orbx))
    ufcmo = ufcmo + wrk(1:n2orbx)

    call rsp1gr(1, kzyva, idum, 0, isyma, 0, 1, etrs,&
                wrk, idum, idum, 1.0d0, 1, udv, ufcmo, xindx,&
                mjwop, wrk, nwrk, .true., .false., .false.)

    deallocate(fcaos, fcmo, ufcmo)

    call qexit('pelib_ifc_cro')

end subroutine pelib_ifc_cro

subroutine pelib_ifc_pecc(aoden, aodencc, converged, t_or_tbar)
    implicit none
#include "priunit.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccorb.h"
#include "qm3.h"
#include "ccfdgeo.h"
#include "ccslvinf.h"
#include "ccinftap.h"

    real*8, parameter :: zero = 0.0d0
    integer, intent(in) :: t_or_tbar
    logical :: converged
    real*8, intent(in) :: aoden(n2bst(isymop)), aodencc(n2bst(isymop))
    real*8 :: energy, ecmcu, eccgrs1
    real*8, allocatable :: denmat(:), fockmat(:)
    character*10 :: model, label
    character*6 :: modelpri

    call qenter('PELIB_IFC_PECC')

    write (lupri,'(1X,A,/)') ' '
    write (lupri,'(1X,A,/)') &
       '******************************************************************'
    write (lupri,'(1X,A,/)') &
       '**** Output from Polarizable Embedding Coupled Cluster module ****'
    write (lupri,'(1X,A,/)') &
       '******************************************************************'

    model = 'CCSD'
    if (ccsd) then
        call around('The Coupled Cluster model is CCSD')
        model = 'CCSD'
        modelpri = 'CCSD'
    end if
    if ((cc2) .and. (.not. mcc2)) then
        call around('The Coupled Cluster model is CC2')
        model = 'CC2'
        modelpri = ' CC2'
    end if
    if (.not. (ccsd .or. cc2)) then
        call quit('PECC only implemented for CCSD and CC2 parametrizations')
    end if
    allocate(denmat(nnbasx), fockmat(nnbasx))
    fockmat = 0.0d0
    call dgefsp(nbas, aoden, denmat)
    if (hffld) then
        call pelib_ifc_energy(denmat, energy)
    else
        call pelib_ifc_fock(denmat, fockmat, energy)
        call pelib_ifc_energy(denmat)
        call around('Writing the Fock matrix to FOCKMAT file')
        call put_to_file('FOCKMAT', nnbasx, fockmat)
    end if
    deallocate(denmat, fockmat)
    ecmcu = eccgrs + energy
    if (abs(ecmcu-eccpr).lt.cvgesol) lslecvg = .true.
    write(lupri,*) 'E(PECC) contribution in iteration', iccslit, ': ', energy
    eccpr = ecmcu
    converged = .false.
    if (lslecvg .and. lsltcvg .and. lsllcvg) then
        converged = .true.
        if (loiter) then
            write(lupri,'(1X,A,I3,A)') 'Maximum inner iterations for '// &
                                       't set to', mxtinit, 'in each outer iteration.'
            write(lupri,'(1X,A,I3,A)') 'Maximum inner iterations for '//&
                                       't-bar set to', mxlinit, 'in each outer iteration.'
        end if
        write(lupri,*) 'PECC equations are converged in ', iccslit, &
                       ' outer iterations'
        write(lupri,*) 'PECC equations are converged in ', nslvinit, &
                       ' inner iterations'
        write(lures,'(12X,A4,A,F20.10)') modelpri, ' Total  energy:             ',&
                                         ecmcu
        write(lures,'(12X,A4,A,F20.10)') modelpri, ' E(PE-CC)     :             ',&
                                         energy
        eccgrs1 = ecmcu
        label = 'PE-'//modelpri//' '
        call wripro(eccgrs1,model,0,label,label,label,label,zero,zero,zero,1,0,0,0)
        label = 'E(PE-CC) '
        call wripro(energy,model,0,label,label,label,label,zero,zero,zero,1,0,0,0)
    else
        iccslit = iccslit + 1
        if (iccslit .gt. mxccslit) then
            write(lupri,*) 'Maximum number of PE-CC iterations:', mxccslit, &
                           'is reached!'
            call quit('Maximum number of PE-CC iterations reached')
        end if
    end if

    call qexit('PELIB_IFC_PECC')

end subroutine pelib_ifc_pecc

subroutine pelib_ifc_transformer(rho1,rho2,ctr1,ctr2,model,isymtr,lr,work,lwork)
    implicit none

#include "mxcent.h"
#include "qmmm.h"
#include "qm3.h"
#include "ccsdsym.h"
#include "priunit.h"
#include "ccsdinp.h"
#include "ccslvinf.h"
#include "ccorb.h"

    real*8, parameter :: zero = 0.0d0, one = 1.0d0, half = 0.5d0, two = 2.0d0
    integer, parameter :: izero = 0
    integer :: lwork, ndim, idldum, isydum, idlino, isymtr
    character*1, intent(in) :: lr
    character*2 :: list
    real*8, dimension(lwork) :: work
    real*8 :: rho1(nt1am(isymtr)), rho2(nt2am(isymtr)), ctr1(nt1am(isymtr)),&
              ctr2(nt2am(isymtr)), ddot, rho1n, rho2n, tal1, tal2, fact
    real*8, allocatable :: denmats(:), gmat(:), eta(:), denmattemp(:), gmattemp(:)
    character*2 :: lisdum
    character*10 :: model
    logical, parameter :: locdeb = .false.
    character*8 :: label

    call qenter('PELIB_IFC_TRANSFORMER')

    if (ipqmmm .gt. 10) then
        write(lupri,*) 'PECC TRANSFORMER : ISYMTR : ', isymtr, ' and LR : ', lr
    end if

    if (ccs) call quit('PECC TRANSFORMER not implemented for CCS')

    if (discex) call quit('PECC TRANSFORMER not implemented for DISCEX')
    if (hffld) call quit('HFFLD should not be here')
    if (ccfixf) then
        call qexit('PELIB_IFC_TRANSFORMER')
        return
    end if

    if (lwork .lt. 0) then
        write(lupri,*) 'Available LWORK: ', lwork
        call quit('Too little work in PECC TRANSFORMER')
    end if

    allocate(denmats(nnbasx),gmat(n2bst(isymtr)),eta(nt1am(isymtr)+nt2am(isymtr)), &
             denmattemp(n2bst(isymtr)),gmattemp(nnbasx))
    eta = zero
    denmattemp = zero
    gmat = zero
    if ((ipqmmm .gt. 10) .or. (locdeb)) then
       rho1n = ddot(nt1am(isymtr),rho1,1,rho1,1)
       rho2n = ddot(nt2am(isymtr),rho2,1,rho2,1)
       write(lupri,*) 'Norm of RHO1 in PECC TRANSFORMER on input: ', rho1n
       write(lupri,*) 'Norm of RHO2 in PECC TRANSFORMER on input: ', rho2n
       rho1n = ddot(nt1am(isymtr),ctr1,1,ctr1,1)
       rho2n = ddot(nt2am(isymtr),ctr2,1,ctr2,1)
       write(lupri,*) 'Norm of C1AM in PECC TRANSFORMER on input: ', rho1n
       write(lupri,*) 'Norm of C2AM in PECC TRANSFORMER on input: ', rho2n
    end if
    call ccmm_d1ao(denmattemp,ctr1,ctr2,trim(model),lr,lisdum,idldum, &
                   isydum,work,lwork)
    call dgefsp(nbas,denmattemp,denmats)
    call pelib_ifc_response(denmats,gmattemp,1)
    call dsptsi(nbas,gmattemp,gmat)
    if ((lr .eq. 'L') .or. (lr .eq. 'F')) then
        label = 'GIVE INT'
        list = 'L0'
        idlino = 1
        call cc_etac(isymtr,label,eta,list,idlino,0,gmat,work,lwork)
    else if ((lr .eq. 'R') .or. (lr .eq. 'P')) then
        label = 'GIVE INT'
        call cc_xksi(eta,label,isymtr,0,gmat,work,lwork)
        if (lr .eq. 'R') then
            call cclr_diascl(eta(nt1am(isymtr)+1),two,isymtr)
        end if
    end if

    if ((locdeb) .or. (ipqmmm .gt. 14)) then
        tal1 = ddot(nt1am(isymtr),eta,1,eta,1)
        tal2 = ddot(nt2am(isymtr),eta(nt1am(isymtr)+1),1,eta(nt1am(isymtr)+1),1)
        write(lupri,*) 'Printing TRANSFORMATION contribution. &
                        Norm2 of singles: ', tal1, 'Norm2 of doubles: ', tal2
    end if

    fact=one

    call daxpy(nt1am(isymtr),fact,eta,1,rho1,1)
    call daxpy(nt2am(isymtr),fact,eta(nt1am(isymtr)+1),1,rho2,1)

    if ((locdeb) .or. (ipqmmm .gt. 14)) then
        tal1 = ddot(nt1am(isymtr),rho1,1,rho1,1)
        tal2 = ddot(nt2am(isymtr),rho2,1,rho2,1)
        write(lupri,*) 'Printing RHO: &
                        Norm2 of singles: ', tal1, 'Norm2 of doubles: ', tal2
    end if
    deallocate(denmats,gmat,eta,denmattemp,gmattemp)

    call qexit('PELIB_IFC_TRANSFORMER')

end subroutine pelib_ifc_transformer

subroutine pelib_ifc_qrtransformer(rho1,rho2,isyres,listb,idlstb,isymtb, &
                                   listc,idlstc,isymtc,model,rsptyp, &
                                   work,lwork)
    implicit none

#include "mxcent.h"
#include "qmmm.h"
#include "qm3.h"
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccslvinf.h"
#include "ccsdsym.h"
#include "priunit.h"

    real*8, parameter :: zero = 0.0d0, half = 0.5d0, one = 1.0d0, two = 2.0d0
    integer, parameter :: izero = 0
    integer :: lwork, ndim, isycb
    real*8 :: work(lwork), rho1(*), rho2(*), ddot,rho1n, rho2n, tal1, tal2, &
              factks, fact
    real*8, allocatable :: b1am(:), b2am(:), c1am(:), c2am(:), denmattemp(:), &
                           denmats(:), eta1(:), eta2(:), eta3(:), gbmat(:), &
                           gbmattemp(:)
    character*(*) :: listb, listc
    integer :: idlstb, isymtb, idlstc, isymtc, isyres, iopt
    logical, parameter :: locdeb = .false.
    logical :: lsame
    character*5 :: model, moddum
    character*8 :: label
    character*1 :: dentyp, rsptyp
    character*2 :: listl

    call qenter('PELIB_IFC_QRTRANSFORMER')

    listl = 'L0'

    if (ccs) call quit('PECC QR TRANSFORMER: Not implemented for CCS')
    if (discex) call quit('PECC QR TRANSFORMER: Not implemented for DISCEX')
    if ((rsptyp .ne. 'F') .and. (rsptyp .ne. 'G') .and. (rsptyp .ne. 'B') .and. &
        (rsptyp .ne. 'K')) then
        write(lupri,*) 'Response flag in QRTRANSFORMER is : ', rsptyp
        write(lupri,*) 'Response flag in QRTRANSFORMER must be F, G, B or K!'
        call quit('Wrong flag in PECC QR TRANSFORMER')
    end if

!   if B=C, just multiply by 2 for second contribution
    if (rsptyp .eq. 'F') then
        factks = one
        lsame = .false.
    else
        lsame = ((listc .eq. listb) .and. (idlstc .eq. idlstb))
        if (lsame) then
            factks = two
        else
            factks = one
        end if
    end if

    if ((locdeb) .or. (ipqmmm .gt. 10)) then
        write(lupri,*) 'PECC QR TRANSFORMER: RSPTYP = ', rsptyp
        write(lupri,*) 'PECC QR TRANSFORMER: ISYRES = ', isyres
        write(lupri,*) 'PECC QR TRANSFORMER: ISYMTB = ', isymtb
        write(lupri,*) 'PECC QR TRANSFORMER: ISYMTC = ', isymtc
        write(lupri,*) 'PECC QR TRANSFORMER: LISTB = ', listb
        write(lupri,*) 'PECC QR TRANSFORMER: LISTC = ', listc
        write(lupri,*) 'PECC QR TRANSFORMER: IDLISTB = ', idlstb
        write(lupri,*) 'PECC QR TRANSFORMER: IDLISTC = ', idlstc
        write(lupri,*) 'PECC QR TRANSFORMER: LSAME = ', lsame
        call flshfo(lupri)
    end if

    isycb = muld2h(isymtb,isymtc)
    if (isycb .ne. isyres) then
        call quit('Symmetry problem in PECC QR TRANSFORMER')
    end if
    if (lwork .lt. 0) then
        write(lupri,*) 'Available memory: ', lwork
        call quit('Too little work in PECC QR TRANSFORMER (1).')
    end if
    if ((.not. lsame) .and. (rsptyp .ne. 'K')) then
        allocate(b1am(nt1am(isymtb)),b2am(nt2am(isymtb)), &
                 c1am(nt1am(isymtc)),c2am(nt2am(isymtc)), &
                 denmattemp(n2bst(isymtb)+n2bst(isymtc)+n2bst(isycb)), &
                 gbmat(n2bst(isymtb)+n2bst(isymtc)+n2bst(isycb)), &
                 denmats(3*nnbasx),gbmattemp(3*nnbasx), &
                 eta1(nt1am(isycb)+nt2am(isycb)), &
                 eta2(nt1am(isycb)+nt2am(isycb)), &
                 eta3(nt1am(isycb)+nt2am(isycb)))
        ndim = 3
        denmattemp = zero
        gbmattemp = zero
        eta1 = zero
        eta2 = zero
        eta3 = zero
        fact = one
    else if ((lsame) .and. (rsptyp .ne. 'K')) then
        allocate(b1am(nt1am(isymtb)),b2am(nt2am(isymtb)), &
                 denmattemp(n2bst(isymtb)+n2bst(isymtc)), &
                 denmats(2*nnbasx),gbmattemp(2*nnbasx),&
                 gbmat(n2bst(isymtb)+n2bst(isymtc)), &
                 eta1(nt1am(isycb)+nt2am(isycb)), &
                 eta2(nt1am(isycb)+nt2am(isycb)))
        ndim = 2
        denmattemp = zero
        gbmattemp = zero
        eta1 = zero
        eta2 = zero
        fact = one
    else if ((.not.lsame) .and. ((rsptyp .eq. 'K'))) then
        allocate(b1am(nt1am(isymtb)),b2am(nt2am(isymtb)), &
                 c1am(nt1am(isymtc)),c2am(nt2am(isymtc)), &
                 denmattemp(n2bst(isymtb)+n2bst(isymtc)), &
                 denmats(2*nnbasx),gbmattemp(2*nnbasx), &
                 gbmat(n2bst(isymtb)+n2bst(isymtc)), &
                 eta1(nt1am(isycb)+nt2am(isycb)), &
                 eta2(nt1am(isycb)+nt2am(isycb)))
        ndim = 2
        denmattemp = zero
        gbmattemp = zero
        eta1 = zero
        eta2 = zero
        fact = one
    else if ((lsame) .and. ((rsptyp .eq. 'K'))) then
        allocate(b1am(nt1am(isymtb)),b2am(nt2am(isymtb)), &
                 denmattemp(n2bst(isymtb)),denmats(nnbasx), &
                 gbmat(n2bst(isycb)),gbmattemp(nnbasx), &
                 eta1(nt1am(isycb)+nt2am(isycb)))
        ndim = 1
        denmattemp = zero
        gbmattemp = zero
        eta1 = zero
        fact = one
    end if

    iopt = 3
    call cc_rdrsp(listb,idlstb,isymtb,iopt,model,b1am,b2am)
    if (.not. lsame) then
        call cc_rdrsp(listc,idlstc,isymtc,iopt,model,c1am,c2am)
    end if

    if ((locdeb) .or. (ipqmmm .gt. 10)) then
        rho1n = ddot(nt1am(isycb),rho1,1,rho1,1)
        rho2n = ddot(nt2am(isycb),rho2,1,rho2,1)
        write(lupri,*) 'Norm of RHO1 in PECC QM TRANSFORMER on input (1): ', &
                       rho1n
        write(lupri,*) 'Norm of RHO2 in PECC QM TRANSFORMER on input (1): ', &
                       rho2n
        rho1n = ddot(nt1am(isycb),b1am,1,b1am,1)
        rho2n = ddot(nt2am(isycb),b2am,1,b2am,1)
        write(lupri,*) 'Norm of B1AM in PECC QM TRANSFORMER on input (1): ', &
                       rho1n
        write(lupri,*) 'Norm of B2AM in PECC QM TRANSFORMER on input (1): ', &
                       rho2n
        if (.not. lsame) then
            rho1n = ddot(nt1am(isycb),c1am,1,c1am,1)
            rho2n = ddot(nt2am(isycb),c2am,1,c2am,1)
            write(lupri,*) 'Norm of C1AM in PECC QM TRANSFORMER on input (1): ', &
                           rho1n
            write(lupri,*) 'Norm of C2AM in PECC QM TRANSFORMER on input (1): ', &
                           rho2n
        end if
    end if

!   check kind of response (for F start with left density)
    if ((rsptyp .eq. 'G') .or. (rsptyp .eq. 'B')) then
        dentyp = 'R'
    else if ((rsptyp .eq. 'K') .or. (rsptyp .eq. 'F')) then
        dentyp = 'L'
    end if
    call ccmm_d1ao(denmattemp,b1am,b2am,model,dentyp,listc,idlstc,isymtc,&
                   work,lwork)
    if (.not. lsame) then
        if (rsptyp .eq. 'F') dentyp = 'R'
        call ccmm_d1ao(denmattemp(n2bst(isymtb)+1),c1am,c2am,model,dentyp, &
                        listb,idlstb,isymtb,work,lwork)
        if ((rsptyp .eq. 'G').or.(rsptyp .eq. 'B')) then
            call ccmm_d2ao(denmattemp(n2bst(isymtb)+n2bst(isymtc)+1),isyres, &
                           listb,idlstb,isymtb,listc,idlstc,isymtc,model, &
                           work,lwork)
        else if (rsptyp .eq. 'F') then
            dentyp = 'Q'
            call ccmm_d1ao(denmattemp(n2bst(isymtb)+n2bst(isymtc)+1),c1am,c2am, &
                           model,dentyp,listb,idlstb,isymtb,work,lwork)
        end if
    else if (lsame) then
        if ((rsptyp .eq. 'G') .or. (rsptyp .eq. 'B')) then
            call ccmm_d2ao(denmattemp(n2bst(isymtb)+1),isyres,listb,idlstb, &
                           isymtb,listc,idlstc,isymtc,model,work,lwork)
        end if
    end if
!   construct effective G operator
    call dgefsp(nbas,denmattemp,denmats)
    if (.not. lsame) then
        call dgefsp(nbas,denmattemp(n2bst(isymtb)+1),denmats(nnbasx+1))
        if (rsptyp .eq. 'B') then
            call dgefsp(nbas,denmattemp(n2bst(isymtb)+n2bst(isymtc)+1), &
                        denmats(2*nnbasx+1))
        else if ((rsptyp .eq. 'F').or.(rsptyp.eq.'G')) then
            call dgefsp(nbas,denmattemp(n2bst(isymtb)+n2bst(isymtc)+1), &
                        denmats(2*nnbasx+1))
        end if
    else if (lsame) then
        if ((rsptyp .eq. 'G') .or. (rsptyp .eq. 'B')) then
            call dgefsp(nbas,denmattemp(n2bst(isymtb)+1),denmats(nnbasx+1))
        end if
    end if
    call pelib_ifc_response(denmats,gbmattemp,ndim)
    call dsptsi(nbas,gbmattemp,gbmat)
    if (.not. lsame) then
        call dsptsi(nbas,gbmattemp(nnbasx+1),gbmat(n2bst(isymtb)+1))
        if ((rsptyp .eq. 'B')) then
            call dsptsi(nbas,gbmattemp(2*nnbasx+1),gbmat(n2bst(isymtb)+ &
                        n2bst(isymtc)+1))
        else if ((rsptyp .eq. 'F').or.(rsptyp.eq.'G')) then
            call dsptsi(nbas,gbmattemp(2*nnbasx+1),gbmat(n2bst(isymtb)+ &
                        n2bst(isymtc)+1))
        end if
    else if (lsame) then
        if ((rsptyp .eq. 'G') .or. (rsptyp .eq. 'B')) then
            call dsptsi(nbas,gbmattemp(nnbasx+1),gbmat(n2bst(isymtb)+1))
        end if
    end if
    if ((locdeb) .or. (ipqmmm .gt. 14)) then
        tal1 = ddot(n2bst(isymtb),gbmat,1,gbmat,1)
        write(lupri,*) 'Print Norm2 GBMAT (1): ', tal1
        if (.not. lsame) then
            tal1 = ddot(n2bst(isymtc),gbmat(n2bst(isymtb)+1),1, &
                        gbmat(n2bst(isymtb)+1),1)
            write(lupri,*) 'Print Norm2 GBMAT (2): ', tal1
            if ((rsptyp.eq.'B').or.(rsptyp.eq.'G').or.(rsptyp.eq.'F')) then
                tal1 = ddot(n2bst(isymtb),gbmat(n2bst(isymtb)+n2bst(isymtc)+1), &
                            1,gbmat(n2bst(isymtb)+n2bst(isymtc)+1),1)
                write(lupri,*) 'Print Norm2 GBCMAT (3): ', tal1
            end if
        else if (lsame) then
            if ((rsptyp .eq. 'G') .or. (rsptyp .eq. 'B')) then
                tal1 = ddot(n2bst(isymtc),gbmat(n2bst(isymtb)+1),1, &
                            gbmat(n2bst(isymtb)+1),1)
                write(lupri,*) 'Print Norm2 GBMAT (2): ', tal1
            end if
        end if
    end if
    label = 'GIVE INT'
    if ((rsptyp .eq. 'G') .or. (rsptyp .eq. 'F')) then
        call cclr_fa(label,isymtb,listc,idlstc,listl,0,gbmat, &
                     work,lwork)
        call daxpy(nt1am(isycb)+nt2am(isycb),1.0d0,work,1,eta1,1)
    else if (rsptyp .eq. 'B') then
        call cccr_aa(label,isymtb,listc,idlstc,gbmat,work,lwork)
        call daxpy(nt1am(isycb)+nt2am(isycb),1.0d0,work,1,eta1,1)
        call cclr_diascl(eta1(nt1am(isycb)+1),two,isycb)
    else if (rsptyp .eq. 'K') then
        call cc_etac(isymtb,label,eta1,listc,idlstc,0,gbmat,work,lwork)
    end if
    if (locdeb .or.(ipqmmm .gt. 14)) then
        tal1 = ddot(nt1am(isycb),eta1,1,eta1,1)
        tal2 = ddot(nt2am(isycb),eta1(nt1am(isycb)+1),1,eta1(nt1am(isycb)+1),1)
        write(lupri,*) 'Printing transformed G^B*C contribution.'
        write(lupri,*) 'Norm2 of singles : ', tal1
        write(lupri,*) 'Norm2 of doubles : ', tal2
    end if
    if (.not. lsame) then
        if (rsptyp .eq. 'G') then
            call cclr_fa(label,isymtc,listb,idlstb,listl,0,gbmat(n2bst(isymtb)+1), &
                         work,lwork)
            call daxpy(nt1am(isycb)+nt2am(isycb),1.0d0,work,1,eta2,1)
        else if (rsptyp .eq. 'B') then
            call cccr_aa(label,isymtc,listb,idlstb,gbmat(n2bst(isymtb)+1),work,lwork)
            call daxpy(nt1am(isycb)+nt2am(isycb),1.0d0,work,1,eta2,1)
            call cclr_diascl(eta2(nt1am(isycb)+1),two,isycb)
        else if ((rsptyp .eq. 'K') .or. (rsptyp .eq. 'F')) then
            call cc_etac(isymtc,label,eta2,listb,idlstb,0,gbmat(n2bst(isymtb)+1), &
                         work,lwork)
        end if
        if (locdeb .or.(ipqmmm .gt. 14)) then
            tal1 = ddot(nt1am(isycb),eta2,1,eta2,1)
            tal2 = ddot(nt2am(isycb),eta2(nt1am(isycb)+1),1,eta2(nt1am(isycb)+1),1)
            write(lupri,*) 'Printing transformed G^C*B contribution.'
            write(lupri,*) 'Norm2 of singles : ', tal1
            write(lupri,*) 'Norm2 of doubles : ', tal2
        end if
        if ((rsptyp .eq. 'G').or.(rsptyp .eq. 'F')) then
            call cc_etac(isycb,label,eta3,listl,0,0,gbmat(n2bst(isymtb)+ &
                         n2bst(isymtc)+1),work,lwork)
        else if ((rsptyp .eq. 'B')) then
            call cc_xksi(eta3,label,isycb,0,gbmat(n2bst(isymtb)+n2bst(isymtc)+1), &
                         work,lwork)
            call cclr_diascl(eta3(nt1am(isycb)+1),two,isycb)
        end if
        if (locdeb .or.(ipqmmm .gt. 14)) then
            if ((rsptyp .eq. 'G').or.(rsptyp .eq. 'B').or.(rsptyp.eq.'F')) then
                tal1 = ddot(nt1am(isycb),eta3,1,eta3,1)
                tal2 = ddot(nt2am(isycb),eta3(nt1am(isycb)+1),1,eta3(nt1am(isycb)+1),1)
                write(lupri,*) 'Printing transformed G^{BC} contribution.'
                write(lupri,*) 'Norm2 of singles : ', tal1
                write(lupri,*) 'Norm2 of doubles : ', tal2
            end if
        end if
    else if (lsame) then
        if (rsptyp .eq. 'G') then
            call cc_etac(isycb,label,eta2,listl,0,0,gbmat(n2bst(isymtb)+1), &
                         work,lwork)
        else if ((rsptyp .eq. 'B')) then
            call cc_xksi(eta2, label, isycb, 0, gbmat(n2bst(isymtb)+1), work,lwork)
            call cclr_diascl(eta2(nt1am(isycb)+1),two,isycb)
        end if
        if (locdeb .or.(ipqmmm .gt. 14)) then
            if ((rsptyp .eq. 'G').or.(rsptyp .eq. 'B')) then
                tal1 = ddot(nt1am(isycb),eta2,1,eta2,1)
                tal2 = ddot(nt2am(isycb),eta2(nt1am(isycb)+1),1,eta2(nt1am(isycb)+1),1)
                write(lupri,*) 'Printing transformed G^{BC} contribution.'
                write(lupri,*) 'Norm2 of singles : ', tal1
                write(lupri,*) 'Norm2 of doubles : ', tal2
            end if
        end if
    end if

    if ((locdeb) .or. (ipqmmm .gt. 14)) then
        tal1 = ddot(nt1am(isycb),eta1,1,eta1,1)
        tal2 = ddot(nt2am(isycb),eta1(nt1am(isycb)+1),1,eta1(nt1am(isycb)+1),1)
        write(lupri,*) 'Printing transformed G^B*C contribution.'
        write(lupri,*) 'Norm2 of singles: ', tal1
        write(lupri,*) 'Norm2 of doubles: ', tal2
    end if
    call daxpy(nt1am(isycb),factks*fact,eta1,1,rho1,1)
    call daxpy(nt2am(isycb),factks*fact,eta1(nt1am(isycb)+1),1,rho2,1)
    if (.not. lsame) then
        call daxpy(nt1am(isycb),one*fact,eta2,1,rho1,1)
        call daxpy(nt2am(isycb),one*fact,eta2(nt1am(isycb)+1),1,rho2,1)
        if ((rsptyp.eq.'G').or.(rsptyp.eq.'B').or.(rsptyp.eq.'F')) then
            call daxpy(nt1am(isycb),one*fact,eta3,1,rho1,1)
            call daxpy(nt2am(isycb),one*fact,eta3(nt1am(isycb)+1),1,rho2,1)
        end if
    else if (lsame) then
        if ((rsptyp .eq. 'G') .or. (rsptyp .eq. 'B')) then
            call daxpy(nt1am(isycb),one*fact,eta2,1,rho1,1)
            call daxpy(nt2am(isycb),one*fact,eta2(nt1am(isycb)+1),1,rho2,1)
        end if
    end if

    if ((locdeb) .or. (ipqmmm .gt. 14)) then
        tal1 = ddot(nt1am(isycb),rho1,1,rho1,1)
        tal2 = ddot(nt2am(isycb),rho2,1,rho2,1)
        write(lupri,*) 'Printing RHO.'
        write(lupri,*) 'Norm2 of singles: ', tal1
        write(lupri,*) 'Norm2 of doubles: ', tal2
    end if

    deallocate(denmattemp,denmats,gbmat,gbmattemp,b1am,b2am,eta1)
    if (allocated(c1am)) deallocate(c1am)
    if (allocated(c2am)) deallocate(c2am)
    if (allocated(eta2)) deallocate(eta2)
    if (allocated(eta3)) deallocate(eta3)

    call qexit('PELIB_IFC_QRTRANSFORMER')

end subroutine pelib_ifc_qrtransformer

end module pelib_interface

subroutine pelib_ifc_start_slaves(runtyp)
    use pelib_interface, only: use_pelib
    integer :: runtyp
#include "iprtyp.h"
#include "maxorb.h"
#include "infpar.h"
    integer, parameter :: iprtyp = POLARIZABLE_EMBEDDING
    call qenter('pelib_ifc_start_slaves')
    if (nodtot >= 1) then
        call mpixbcast(iprtyp, 1, 'INTEGER', master)
        call mpixbcast(runtyp, 1, 'INTEGER', master)
    end if
    call qexit('pelib_ifc_start_slaves')
end subroutine pelib_ifc_start_slaves

#else

module pelib_interface

    implicit none

    private

    public :: use_pelib, pelib_ifc_gspol
    public :: pelib_ifc_do_mep, pelib_ifc_do_mep_noqm, pelib_ifc_do_cube
    public :: pelib_ifc_do_infld, pelib_ifc_do_lf
    public :: pelib_ifc_activate, pelib_ifc_deactivate
    public :: pelib_ifc_init, pelib_ifc_finalize, pelib_ifc_input_reader
    public :: pelib_ifc_fock, pelib_ifc_energy, pelib_ifc_response, pelib_ifc_london
    public :: pelib_ifc_molgrad, pelib_ifc_infld, pelib_ifc_lf, pelib_ifc_localfield
    public :: pelib_ifc_mep, pelib_ifc_mep_noqm, pelib_ifc_cube
    public :: pelib_ifc_set_mixed, pelib_ifc_mixed
    public :: pelib_ifc_do_savden, pelib_ifc_do_twoints
    public :: pelib_ifc_save_density, pelib_ifc_twoints
    public :: pelib_ifc_get_num_fragment_nuclei
#if defined(VAR_MPI)
    public :: pelib_ifc_slave
#endif
    public :: pelib_ifc_grad, pelib_ifc_lin, pelib_ifc_lr, pelib_ifc_qro
    public :: pelib_ifc_cro, pelib_ifc_rspmcqr
    public :: pelib_ifc_pecc
    public :: pelib_ifc_transformer, pelib_ifc_qrtransformer
    public :: pelib_ifc_qrtest

contains

logical function use_pelib()
    use_pelib = .false.
end function use_pelib

logical function pelib_ifc_gspol()
    call quit('using dummy PElib interface routines')
end function pelib_ifc_gspol

subroutine pelib_ifc_set_mixed(do_mixed)
    logical :: do_mixed
    call quit('using dummy PElib interface routines')
end subroutine pelib_ifc_set_mixed

logical function pelib_ifc_do_mep()
    pelib_ifc_do_mep = .false.
end function pelib_ifc_do_mep

logical function pelib_ifc_do_mep_noqm()
    pelib_ifc_do_mep_noqm = .false.
end function pelib_ifc_do_mep_noqm

logical function pelib_ifc_do_cube()
    call quit('using dummy PElib interface routines')
end function pelib_ifc_do_cube

logical function pelib_ifc_do_infld()
    call quit('using dummy PElib interface routines')
end function pelib_ifc_do_infld

logical function pelib_ifc_do_lf()
    call quit('using dummy PElib interface routines')
end function pelib_ifc_do_lf

logical function pelib_ifc_do_savden()
    call quit('using dummy PElib interface routines')
end function pelib_ifc_do_savden

logical function pelib_ifc_do_twoints()
    call quit('using dummy PElib interface routines')
end function pelib_ifc_do_twoints

subroutine pelib_ifc_activate()
    call qenter('pelib_ifc_activate')
    call quit('PElib not compiled, please use -DENABLE_PELIB=ON to enable PElib')
    call qexit('pelib_ifc_activate')
end subroutine pelib_ifc_activate

subroutine pelib_ifc_deactivate()
    call qenter('pelib_ifc_deactivate')
    call quit('using dummy PElib interface routines')
    call qexit('pelib_ifc_deactivate')
end subroutine pelib_ifc_deactivate

subroutine pelib_ifc_input_reader(word)
    character(len=7), intent(in) :: word
    call qenter('pelib_ifc_input_reader')
    call quit('using dummy PElib interface routines')
    call qexit('pelib_ifc_input_reader')
end subroutine pelib_ifc_input_reader

subroutine pelib_ifc_init()
    call qenter('pelib_ifc_init')
    call quit('using dummy PElib interface routines')
    call qexit('pelib_ifc_init')
end subroutine pelib_ifc_init

subroutine pelib_ifc_finalize()
    call qenter('pelib_ifc_finalize')
    call quit('using dummy PElib interface routines')
    call qexit('pelib_ifc_finalize')
end subroutine pelib_ifc_finalize

subroutine pelib_ifc_fock(denmats, fckmats, energy)
    real*8, dimension(*), intent(in) :: denmats
    real*8, dimension(*), intent(out) :: fckmats
    real*8, intent(out) :: energy
    call qenter('pelib_ifc_fock')
    call quit('using dummy PElib interface routines')
    call qexit('pelib_ifc_fock')
end subroutine pelib_ifc_fock

subroutine pelib_ifc_mixed(denmats, fckmats, energy)
    real*8, dimension(*), intent(in) :: denmats
    real*8, dimension(*), intent(in) :: fckmats
    real*8, intent(in) :: energy
    call qenter('pelib_ifc_mixed')
    call quit('using dummy PElib interface routines')
    call qexit('pelib_ifc_mixed')
end subroutine pelib_ifc_mixed

subroutine pelib_ifc_energy(denmats)
    real*8, dimension(*), intent(in) :: denmats
    call qenter('pelib_ifc_energy')
    call quit('using dummy PElib interface routines')
    call qexit('pelib_ifc_energy')
end subroutine pelib_ifc_energy

subroutine pelib_ifc_molgrad(denmats, molgrad)
    real*8, dimension(*), intent(in) :: denmats
    real*8, dimension(*), intent(out) :: molgrad
    call qenter('pelib_ifc_molgrad')
    call quit('using dummy PElib interface routines')
    call qexit('pelib_ifc_molgrad')
end subroutine pelib_ifc_molgrad

subroutine pelib_ifc_infld()
    call qenter('pelib_ifc_infld')
    call quit('using dummy PElib interface routines')
    call qexit('pelib_ifc_infld')
end subroutine pelib_ifc_infld

subroutine pelib_ifc_response(denmats, fckmats, nmats)
    integer, intent(in) :: nmats
    real*8, dimension(*), intent(in) :: denmats
    real*8, dimension(*), intent(out) :: fckmats
    call qenter('pelib_ifc_response')
    call quit('using dummy PElib interface routines')
    call qexit('pelib_ifc_response')
end subroutine pelib_ifc_response

subroutine pelib_ifc_london(fckmats)
    real*8, dimension(*), intent(out) :: fckmats
    call qenter('pelib_ifc_london')
    call quit('using dummy PElib interface routines')
    call qexit('pelib_ifc_london')
end subroutine pelib_ifc_london

subroutine pelib_ifc_localfield(eefmats)
    real*8, dimension(:), intent(in) :: eefmats
    call qenter('pelib_ifc_localfield')
    call quit('using dummy PElib interface routines')
    call qexit('pelib_ifc_localfield')
end subroutine pelib_ifc_localfield

subroutine pelib_ifc_lf()
    call qenter('pelib_ifc_lf')
    call quit('using dummy PElib interface routines')
    call qexit('pelib_ifc_lf')
end subroutine pelib_ifc_lf

subroutine pelib_ifc_mep(denmats)
    real*8, dimension(*), intent(in) :: denmats
    call qenter('pelib_ifc_mep')
    call quit('using dummy PElib interface routines')
    call qexit('pelib_ifc_mep')
end subroutine pelib_ifc_mep

subroutine pelib_ifc_mep_noqm()
    call qenter('pelib_ifc_mep_noqm')
    call quit('using dummy PElib interface routines')
    call qexit('pelib_ifc_mep_noqm')
end subroutine pelib_ifc_mep_noqm

subroutine pelib_ifc_cube(denmats, idx)
    real*8, dimension(*), intent(in) :: denmats
    integer, intent(in) :: idx
    call qenter('pelib_ifc_cube')
    call quit('using dummy PElib interface routines')
    call qexit('pelib_ifc_cube')
end subroutine pelib_ifc_cube

subroutine pelib_ifc_save_density(ao_denmat, mo_fckmat, mo_coefficients)
    real*8, dimension(*), intent(in) :: ao_denmat
    real*8, dimension(*), intent(in) :: mo_fckmat
    real*8, dimension(*), intent(in) :: mo_coefficients
    call qenter('pelib_ifc_save_density')
    call quit('using dummy PElib interface routines')
    call qexit('pelib_ifc_save_density')
end subroutine pelib_ifc_save_density

subroutine pelib_ifc_twoints(work, lwork)
    real*8, dimension(*), intent(in) :: work
    integer, intent(in) :: lwork
    call qenter('pelib_ifc_twoints')
    call quit('using dummy PElib interface routines')
    call qexit('pelib_ifc_twoints')
end subroutine pelib_ifc_twoints

integer function pelib_ifc_get_num_fragment_nuclei()
    use pde_utils, only: pde_get_num_fragment_nuclei
    call quit('using dummy PElib interface routines')
end function pelib_ifc_get_num_fragment_nuclei

#if defined(VAR_MPI)
subroutine pelib_ifc_slave(runtype)
    integer, intent(in) :: runtype
    call qenter('pelib_ifc_slave')
    call quit('using dummy PElib interface routines')
    call qexit('pelib_ifc_slave')
end subroutine pelib_ifc_slave
#endif

subroutine pelib_ifc_grad(cref, cmo, cindx, dv, grd, energy, wrk, nwrk)
    integer :: nwrk
    real*8 :: energy
    real*8, dimension(*) :: cref, cmo, cindx, dv, grd
    real*8, dimension(*) :: wrk
    call qenter('pelib_ifc_grad')
    call quit('using dummy PElib interface routines')
    call qexit('pelib_ifc_grad')
end subroutine pelib_ifc_grad

subroutine pelib_ifc_lin(ncsim, nosim, bcvecs, bovecs, cref, cmo, cindx, dv, dtv,&
                         scvecs, sovecs, orblin, wrk, nwrk)
    logical :: orblin
    integer :: ncsim, nosim, nwrk
    integer, dimension(*) :: cindx
    real*8, dimension(*) :: bcvecs, bovecs, scvecs, sovecs
    real*8, dimension(*) :: cmo, cref, dv, dtv
    real*8, dimension(*) :: wrk
    call qenter('pelib_ifc_lin')
    call quit('using dummy PElib interface routines')
    call qexit('pelib_ifc_lin')
end subroutine pelib_ifc_lin

subroutine pelib_ifc_lr(ncsim, nosim, bcvecs, bovecs, cref, cmo, cindx, udv,&
                        dv, udvtr, dvtr, dtv, dtvtr, scvecs, sovecs, wrk, nwrk)
    integer :: ncsim, nosim, nwrk
    integer, dimension(*) :: cindx
    real*8, dimension(*) :: bcvecs, bovecs
    real*8, dimension(*) :: cref, cmo, udv, dv
    real*8, dimension(*) :: udvtr, dvtr, dtv, dtvtr
    real*8, dimension(*) :: scvecs, sovecs
    real*8, dimension(nwrk) :: wrk
    call qenter('pelib_ifc_lr')
    call quit('using dummy PElib interface routines')
    call qexit('pelib_ifc_lr')
end subroutine pelib_ifc_lr

subroutine pelib_ifc_qro(vecb, vecc, etrs, xindx, zymb, zymc, udv, wrk, nwrk,&
                         kzyva, kzyvb, kzyvc, isyma, isymb, isymc, cmo, mjwop)
    integer :: kzyva, kzyvb, kzyvc
    integer :: isyma, isymb, isymc
    integer :: nwrk
    real*8, dimension(*) :: wrk
    real*8, dimension(*) :: etrs
    real*8, dimension(*) :: vecb
    real*8, dimension(*) :: vecc
    real*8, dimension(*) :: cmo
    real*8, dimension(*) :: zymb, zymc
    real*8, dimension(*) :: udv
    real*8, dimension(*) :: xindx
    integer, dimension(*) :: mjwop
    call qenter('pelib_ifc_qro')
    call quit('using dummy PElib interface routines')
    call qexit('pelib_ifc_qro')
end subroutine pelib_ifc_qro

subroutine pelib_ifc_cro(vecb, vecc, vecd, etrs, xindx, zymb, zymc, zymd, udv,&
                         wrk, nwrk, kzyva, kzyvb, kzyvc, kzyvd, isyma, isymb,&
                         isymc, isymd, cmo,mjwop)
    integer :: kzyva, kzyvb, kzyvc, kzyvd
    integer :: isyma, isymb, isymc, isymd
    integer :: nwrk
    real*8, dimension(*) :: wrk
    real*8, dimension(*) :: etrs
    real*8, dimension(*) :: vecb
    real*8, dimension(*) :: vecc
    real*8, dimension(*) :: vecd
    real*8, dimension(*) :: cmo
    real*8, dimension(*) :: zymb, zymc, zymd
    real*8, dimension(*) :: udv
    real*8, dimension(*) :: xindx
    integer, dimension(*) :: mjwop
    call qenter('pelib_ifc_cro')
    call quit('using dummy PElib interface routines')
    call qexit('pelib_ifc_cro')
end subroutine pelib_ifc_cro

subroutine pelib_ifc_pecc(aoden,converged,work,lwork)
    use pelib_interface, only: pelib_ifc_fock
    integer :: lwork
    real*8 :: work(lwork), aoden(*)
    logical :: converged
    call qenter('PELIB_IFC_PECC')
    call quit('using dummy PElib interface routines')
    call qexit('PELIB_IFC_PECC')
end subroutine pelib_ifc_pecc

subroutine pelib_ifc_transformer(rho1,rho2,ctr1,ctr2,model,isymtr,lr,work,lwork)
    integer :: lwork, isymtr
    character*1, intent(in) :: lr
    real*8, dimension(*) :: work
    real*8 :: rho1(*), rho2(*), ctr1(*), ctr2(*)
    character*10 :: model
    call qenter('PELIB_IFC_TRANSFORMER')
    call quit('using dummy PElib interface routines')
    call qexit('PELIB_IFC_TRANSFORMER')
end subroutine pelib_ifc_transformer

subroutine pelib_ifc_qrtransformer(rho1,rho2,isyres,listb,idlstb,isymtb, &
   &                               listc,idlstc,isymtc,model,rsptyp, &
   &                               work,lwork)
    integer :: lwork
    real*8 :: work(*), rho1(*), rho2(*)
    character*(*) :: listb, listc
    integer :: idlstb, isymtb, idlstc, isymtc
    character*5 :: model
    character*1 :: rsptyp
    call qenter('PELIB_IFC_QRTRANSFORMER')
    call quit('using dummy PElib interface routines')
    call qexit('PELIB_IFC_QRTRANSFORMER')
end subroutine pelib_ifc_qrtransformer

end module pelib_interface

subroutine pelib_ifc_start_slaves(runtyp)
    integer :: runtyp
    call qenter('pelib_ifc_start_slaves')
    call quit('using dummy PElib interface routines')
    call qexit('pelib_ifc_start_slaves')
end subroutine pelib_ifc_start_slaves

subroutine pelib_ifc_pecc(aoden,converged,work,lwork)
    use pelib_interface, only: pelib_ifc_fock
    implicit none
    integer :: lwork
    real*8 :: work(lwork), aoden(*)
    logical :: converged
    call qenter('PELIB_IFC_PECC')
    call quit('using dummy PElib interface routines')
    call qexit('PELIB_IFC_PECC')
end subroutine pelib_ifc_pecc

#endif
