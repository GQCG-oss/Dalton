!
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2016 (2015), see http://daltonprogram.org"
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
#if defined(BUILD_PELIB)
module pelib_interface

    implicit none

    private

    public :: use_pelib, pelib_ifc_gspol
    public :: pelib_ifc_dolf
    public :: pelib_ifc_activate, pelib_ifc_deactivate
    public :: pelib_ifc_init, pelib_ifc_finalize, pelib_ifc_input_reader
    public :: pelib_ifc_fock, pelib_ifc_energy, pelib_ifc_response, pelib_ifc_london
    public :: pelib_ifc_molgrad, pelib_ifc_lf
#if defined(VAR_MPI)
    public :: pelib_ifc_slave
#endif
    ! TODO: update the following interface routines
    public :: pelib_ifc_grad, pelib_ifc_lin, pelib_ifc_lr, pelib_ifc_qro

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

logical function pelib_ifc_dolf()
  use pe_variables, only: pe_lf
  if (pe_lf) then
      pelib_ifc_dolf = .true.
  else
      pelib_ifc_dolf = .false.
  end if
end function pelib_ifc_dolf

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
    if (.not. use_pelib()) call quit('PElib not active')
    call pe_input_reader(word, lucmd)
    call qexit('pelib_ifc_input_reader')
end subroutine pelib_ifc_input_reader

subroutine pelib_ifc_init()
    use pe_variables, only: peqm
    use polarizable_embedding, only: pe_init
#include "priunit.h"
#include "mxcent.h"
#include "nuclei.h"
    call qenter('pelib_ifc_init')
    if (.not. use_pelib()) call quit('PElib not active')
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
                   triang=.true.,       &
                   ndim=nbast,          &
                   nmats=1,             &
                   denmats=denmats,     &
                   fckmats=fckmats,     &
                   expvals=energies)
    energy = energies(1)
    call qexit('pelib_ifc_fock')
end subroutine pelib_ifc_fock

subroutine pelib_ifc_energy(denmats)
    use polarizable_embedding, only: pe_master
#include "inforb.h"
    real*8, dimension(nnbasx), intent(in) :: denmats
    call qenter('pelib_ifc_energy')
    if (.not. use_pelib()) call quit('PElib not active')
#if defined(VAR_MPI)
    call pelib_ifc_start_slaves(2)
#endif
    call pe_master(runtype='print_energy', &
                   triang=.true.,          &
                   ndim=nbast,             &
                   nmats=1,                &
                   denmats=denmats)
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
                   triang=.true.,                &
                   ndim=nbast,                   &
                   nmats=1,                      &
                   denmats=denmats,              &
                   expvals=molgrad)
        call qexit('pelib_ifc_molgrad')
end subroutine pelib_ifc_molgrad

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
                   triang=.true.,              &
                   ndim=nbast,                 &
                   nmats=nmats,                &
                   denmats=denmats,            &
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
                   triang=.true.,       &
                   ndim=nbast,          &
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
                   triang=.true.,       &
                   ndim=nbast,          &
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

#if defined(VAR_MPI)
subroutine pelib_ifc_slave(runtype)
    use polarizable_embedding, only: pe_slave
    integer, intent(in) :: runtype
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
    else if (runtype == 8) then
        call pe_slave('effdipole')
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
                 & scvecs, sovecs, orblin, wrk, nwrk)
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
              & cindx, wrk, nwrk)

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
                 & wrk, nwrk)
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
                  & dummy, cindx, wrk, nwrk)
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
                       & dv, udvtr, dvtr, dtv, dtvtr, scvecs, sovecs, wrk, nwrk)
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
                      & udvtr, dvtr, dtv, dtvtr, scvecs, wrk, nwrk)
    end if

    if (nosim > 0) then
        call pe_rsplno(nosim, bovecs, cref, cmo, cindx,&
                      & udv, dv, udvtr, dvtr, sovecs, wrk, nwrk)
    end if

    call qexit('pelib_ifc_lr')

end subroutine pelib_ifc_lr

subroutine pe_rsplnc(ncsim, bcvecs, cref, cmo, cindx, udv, dv,&
                    & udvtr, dvtr, dtv, dtvtr, scvecs, wrk, nwrk)
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
                   & udtv, dummy, 0, 0, tdm, norho2, cindx, wrk, 1, nwrk)
        udtv = - 1.0d0 * udtv

        if ( ncsim > 0 ) then
            allocate(fdtvaos(nnbasx*ncsim))
            fdtvaos = 0.0d0
            allocate(dtvao(n2basx))
            dtvao = 0.0d0
            do i = 1, ncsim
               j = (i - 1) * nnbasx + 1
                call fckden2(.false., .true., dummy, dtvao, cmo,&
                            & udtv(:,i), wrk, nwrk)
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
                       & idummy, .false.)
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
              & fpeac, tfxcacs, tfpeac, cindx, wrk, nwrk)
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
                           & lupri)
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
                           & lupri)
            end if
        end do
        deallocate(fupe, fuxcs)

        if (locdeb) then
            write(lupri,*)' linear transformed conf. vector'
            write(lupri,*)' *** after slvsor in pe_rsplnc *** '
            call output(scvecs, 1, kzyvar, 1, ncsim, kzyvar, ncsim, 1,&
                       & lupri)
        end if
    end if

    if (ncref /= kzconf) then
        call quit('ERROR in pe_rsplnc: ncref /= kzconf')
    end if

    call qexit('pe_rsplnc')

end subroutine pe_rsplnc

subroutine pe_rsplno(nosim, bovecs, cref, cmo, cindx, udv, dv,&
                    & udvtr, dvtr, sovecs, wrk, nwrk)
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
                 & ' systems not implemented')
    end if

    lopen = .false.

    if (.not. tdhf) then
        ! Read Fg = V - <0|F|0>Fe from file
        allocate(fpemo(nnorbx))
        if (lusifc <= 0) then
            call gpopen(lusifc, 'SIRIFC', 'OLD', ' ', 'UNFORMATTED',&
                       & idummy, .false.)
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
                  & dummy, txyoacs, dummy, cindx, wrk, nwrk)
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
                    & kzyva, kzyvb, kzyvc, isyma, isymb, isymc, cmo, mjwop)
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
               & wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
    call dgefsp(nbast, udcao, dcaos(nnbasx+1:2*nnbasx))

    !  D(2k)
    udcao = 0.0d0
    call cdens1(isymc, cmo, zymc, udcao, wrk, nwrk)
    call dgefsp(nbast, udcao, dcaos(2*nnbasx+1:3*nnbasx))

    !  D(2k,1k)
    udcao = 0.0d0
    call cdens2(isymc, isymb, cmo, zymc, zymb, udcao,&
               & wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
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
               & wrk, idum, idum, 1.0d0, 1, udv, ufcmo, xindx,&
               & mjwop, wrk, nwrk, .true., .false., .false.)

    deallocate(fcaos, fcmo, ufcmo)

    call qexit('pelib_ifc_qro')

end subroutine pelib_ifc_qro

end module pelib_interface

subroutine pelib_ifc_start_slaves(runtyp)
    use pelib_interface, only: use_pelib
    integer :: runtyp
#include "iprtyp.h"
#include "maxorb.h"
#include "infpar.h"
    integer, parameter :: iprtyp = POLARIZABLE_EMBEDDING
    call qenter('pelib_ifc_start_slaves')
    if (.not. use_pelib()) call quit('PElib not active')
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
    public :: pelib_ifc_activate, pelib_ifc_deactivate
    public :: pelib_ifc_init, pelib_ifc_finalize, pelib_ifc_input_reader
    public :: pelib_ifc_fock, pelib_ifc_energy, pelib_ifc_response, pelib_ifc_london
    public :: pelib_ifc_molgrad
#if defined(VAR_MPI)
    public :: pelib_ifc_slave
#endif
    ! TODO: update the following interface routines
    public :: pelib_ifc_grad, pelib_ifc_lin, pelib_ifc_lr, pelib_ifc_qro
    public :: pelib_ifc_cro

contains

logical function use_pelib()
    use_pelib = .false.
end function use_pelib

logical function pelib_ifc_gspol()
    call quit('using dummy PElib interface routines')
end function pelib_ifc_gspol

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
                 & scvecs, sovecs, orblin, wrk, nwrk)
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
                       & dv, udvtr, dvtr, dtv, dtvtr, scvecs, sovecs, wrk, nwrk)
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
                    & kzyva, kzyvb, kzyvc, isyma, isymb, isymc, cmo, mjwop)
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
                    & wrk, nwrk, kzyva, kzyvb, kzyvc, kzyvd, isyma, isymb,&
                    & isymc, isymd, cmo,mjwop)
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

end module pelib_interface

subroutine pelib_ifc_start_slaves(runtyp)
    integer :: runtyp
    call qenter('pelib_ifc_start_slaves')
    call quit('using dummy PElib interface routines')
    call qexit('pelib_ifc_start_slaves')
end subroutine pelib_ifc_start_slaves

#endif
