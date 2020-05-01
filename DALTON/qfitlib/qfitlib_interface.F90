!------------------------------------------------------------------------------
!> @brief Interface between DALTON and QFITLIB
!!
!! @author Casper Steinmann
#if defined (BUILD_QFITLIB)
module qfitlib_interface
  implicit none
  private :: qfitlib_ifc_compute_dipole

  public :: qfitlib_ifc_input_reader
  public :: qfitlib_ifc_initialize
  public :: qfitlib_ifc_finalize
  public :: qfitlib_ifc_information
  public :: qfitlib_ifc_fit
  public :: qfitlib_ifc_results

  public :: qfitlib_ifc_sirius_fit
  public :: qfitlib_ifc_response_fit

#if defined(VAR_MPI)
  public :: qfitlib_ifc_slave ! entry-point for parallel nodes
#endif

contains

!------------------------------------------------------------------------------
!> @brief DALTON input reader interface
!!
!! @author Casper Steinmann
subroutine qfitlib_ifc_input_reader(word)
  use qfit_input_readers, only: dalton_input
#include "priunit.h"
  character(len=7), intent(inout) :: word
  call qenter('qfitlib_ifc_input_reader')
  call dalton_input(word, lucmd, lupri)
  call qexit('qfitlib_ifc_input_reader')
end subroutine qfitlib_ifc_input_reader


!------------------------------------------------------------------------------
!> @brief Interface subroutine to initialize QFITLIB from DALTON
!!
!! Initializes QFITLIB with geometrical and molecular data.
!!
!! @note that constraints on the total dipole (and/or) quadrupole should
!! be set just before the fit is calculated (see helper routines calling
!! qfit_fit)
subroutine qfitlib_ifc_initialize
  use qfit, only : qfit_initialize
! mxcent used in other includes
! nuclei.h : cord, charge, nucdep
! orgcom.h : cmxyz
#include "mxcent.h"
#include "nuclei.h"
#include "orgcom.h"
  integer :: mol_charge
  call molchr(mol_charge)
  call qenter('qfitlib_ifc_initialize')
  call qfit_initialize(cord(1:3,1:nucdep), charge(1:nucdep), mol_charge, RCM=cmxyz)
  call qexit('qfitlib_ifc_initialize')
end subroutine qfitlib_ifc_initialize


!------------------------------------------------------------------------------
!> @brief Interface subroutine to finalize QFITLIB from DALTON
subroutine qfitlib_ifc_finalize
  use qfit, only : qfit_finalize
  call qenter('qfitlib_ifc_finalize')
  call qfit_finalize
  call qexit('qfitlib_ifc_finalize')
end subroutine qfitlib_ifc_finalize

!------------------------------------------------------------------------------
!> @brief Interface subroutine to print QFITLIB information from DALTON
subroutine qfitlib_ifc_information
  use qfit, only : qfit_print_info
  call qenter('qfitlib_ifc_information')
  CALL HEADER('Charge and moment fitting (QFITLIB) settings',-1)
  call qfit_print_info
  call qexit('qfitlib_ifc_information')
end subroutine qfitlib_ifc_information

!------------------------------------------------------------------------------
!> @brief Interface subroutine to compute a fit based on an input dens
subroutine qfitlib_ifc_fit(dens)
  use qfit, only : qfit_fit
#include "iprtyp.h"
#include "maxorb.h"
#include "inforb.h"
#include "infpar.h"
  integer, parameter :: iprtyp = QFIT_WORK
  integer, parameter :: dummy = 0
  real*8, dimension(n2basx), intent(in) :: dens
  call qenter('q_ifc_fit')
#if defined(VAR_MPI)
  call mpixbcast(iprtyp, 1, 'INTEGER', master)
  call mpixbcast(dummy, 1, 'INTEGER', master)
  call mpixbcast(n2basx, 1, 'INTEGER', master)
  call mpixbcast(dens, n2basx, 'DOUBLE', master)
#endif
  call qfit_fit(dens)
  call qexit('q_ifc_fit')
end subroutine qfitlib_ifc_fit

!------------------------------------------------------------------------------
!> @brief Interface subroutine to compute a fit based on an input density from SIRIUS
!!
!! The density is formed from the CMOs and from the density we compute the total dipole.
subroutine qfitlib_ifc_sirius_fit(cmo, wrk, nwrk)
  use qfit, only : qfit_set_dipole
#include "mxcent.h"
#include "maxorb.h"
#include "infinp.h"
#include "inforb.h"
#include "moldip.h"
! parameter variables
  integer :: nwrk
  real*8, dimension(*), intent(in) :: cmo
  real*8, dimension(nwrk), intent(in) :: wrk
! local variables
  real*8, dimension(:), allocatable :: dv, dcao, dvao, dens
  integer :: dummy, i, ii
  call qenter('qfitlib_ifc_sirius_fit')
  call qfitlib_ifc_information
  allocate(dcao(n2basx))
  allocate(dvao(n2basx))
  allocate(dv(nnashx))
  if (nasht == 1) then
      dv(1) = 1.0d0
  else if (hsrohf) then
    do i = 1, nasht
      ii = i*(i+1)/2
      dv(ii) = 1.0d0
    end do
  endif
  call fckden((nisht > 0), (nasht > 0), dcao, dvao, cmo, dv, wrk, nwrk)
  if (nisht == 0) dcao = 0.0d0
  if (nasht > 0) dcao = dcao + dvao ! only add active shells if they are present
  deallocate(dv)

  call qfitlib_ifc_compute_dipole(dcao)  ! stores dipole in DIP0 in "moldip.h"
  deallocate(dvao)

  call qfit_set_dipole(dip0)
  call qfitlib_ifc_fit(dcao)
  deallocate(dcao)
  call qexit('qfitlib_ifc_sirius_fit')
end subroutine qfitlib_ifc_sirius_fit

!------------------------------------------------------------------------------
subroutine qfitlib_ifc_response_fit(cmo, wop, udv, wrk, nwrk)
  use qfit, only : qfit_set_transition_dipole
#include "mxcent.h"
#include "inforb.h"
#include "dipole.h"
! parameter variables
  integer :: nwrk
  real*8, dimension(*), intent(in) :: cmo, wop, udv
  real*8, dimension(nwrk), intent(inout) :: wrk
! local variables
  real*8, dimension(:), allocatable :: ubove, dcao, dvao
  call qenter('qfitlib_ifc_response_fit')
  allocate(ubove(n2orbx))
  call rspzym(1, wop, ubove)
  ubove = -ubove

  allocate(dcao(n2basx))
  allocate(dvao(n2basx))
  call deq27(cmo, ubove, udv, dcao, dvao, wrk, nwrk)
  if (nasht > 0) dcao = dcao + dvao ! only add active shells if they are present
  deallocate(dvao)
  deallocate(ubove)

  call qfitlib_ifc_compute_dipole(dcao)    ! stores total dipole in DIP0 in "moldip.h"
                                           ! electric dipole stored in dipme
  call qfit_set_transition_dipole(dipme)
  call qfitlib_ifc_fit(dcao)

  deallocate(dcao)
  call qexit('qfitlib_ifc_response_fit')
end subroutine qfitlib_ifc_response_fit

!------------------------------------------------------------------------------
!> @brief Computes the dipole moment from any square AO density
!!
!! The total dipole moment \f$\mu=\mu_e + \mu_n\f$ from the
!! electron density and nuclei, respectively. The electronic
!! contribution is computed for a cartesian direction \f$\alpha\f$
!! \f[
!!   \mu_\alpha = - \mathrm{Tr}\,\mathbf{D}\mathbf{V}
!! \f]
!! where \f$V\f$ are the dipole integrals in cartesian direction \f$\alpha\f$
!! @param dcao square AO density matrix
subroutine qfitlib_ifc_compute_dipole(dcao)
#include "mxcent.h"
#include "dipole.h"
#include "inftap.h"
#include "inforb.h"
#include "priunit.h"
! parameter variables
  real*8, dimension(:) :: dcao
! local variables
  logical :: lopen
  real*8 :: dummy
  real*8, dimension(:), allocatable :: dipints
  real*8, dimension(:), allocatable :: dens
  real*8, dimension(:), allocatable :: dvao ! temporary for SP storage of dcao
  real*8 :: ddot ! function ddot
  call qenter('qfitlib_ifc_compute_dipole')

  allocate(dens(nnbasx))
  allocate(dvao(n2basx))
  dvao = 0.0d0
  dens = 0.0d0
  ! first convert to packed symmetric storage
  call dgefsp(nbast,dcao,dvao)       ! store in some other (ASP?) format
  call pksym1(dvao,dens,nbas,nsym,1) ! store as upper triangular matrix in dens
  deallocate(dvao)

  ! we compute elements we need for the total dipole moment
  ! nuclear contribution first (dipmn)
  dummy = 0.0d0
  call dipnuc(dummy,dummy,0,.false.)

  ! then electronic contribution (dipme)
  allocate(dipints(nnbasx))
  if (luprop <= 0) then
    call gpopen(luprop, 'AOPROPER', 'OLD', 'SEQUENTIAL', 'UNFORMATTED', 0, .false.)
    lopen = .true.
  end if
  call mollab('XDIPLEN ', luprop, lupri)
  call readt(luprop, nnbasx, dipints)
  dipme(1) = -ddot(nnbasx, dens, 1, dipints, 1)

  call mollab('YDIPLEN ', luprop, lupri)
  call readt(luprop, nnbasx, dipints)
  dipme(2) = -ddot(nnbasx, dens, 1, dipints, 1)

  call mollab('ZDIPLEN ', luprop, lupri)
  call readt(luprop, nnbasx, dipints)
  dipme(3) = -ddot(nnbasx, dens, 1, dipints, 1)

  ! compute total dipole and update all common blocks
  call dp0sum

  if (lopen) then
    call gpclose(luprop,'KEEP')
  end if
  deallocate(dipints)
  deallocate(dens)
  call qexit('qfitlib_ifc_compute_dipole')
end subroutine qfitlib_ifc_compute_dipole

subroutine qfitlib_ifc_results
      use qfit, only : qfit_get_results
      use qfit_variables, only : qfit_multipole_rank
#include "mxcent.h"
#include "nuclei.h"
#include "priunit.h"
  real*8, dimension(:), allocatable :: charges
  real*8, dimension(:), allocatable :: dipoles
  real*8, dimension(:), allocatable :: quadrupoles
  real*8 :: rmsd_value
  integer :: m, k
  call qenter('qfitlib_ifc_results')
  CALL HEADER('Potential fitted multipole moments (QFITLIB)',-1)

  allocate( charges( nucdep ) )
  allocate( dipoles( 3*nucdep ) )
  allocate( quadrupoles(6*nucdep) )
  call qfit_get_results( charges, dipoles, quadrupoles, rmsd_value )

  write(lupri,'(a)') " Charges:"
  write(lupri,'(a,6a,a12)') "  ", "      ", "        Q   "
  do m = 1, size(charges)
      write(lupri,'(a,a6,f12.6)') '@ ',namdep(m), charges(m)
  enddo
  write(lupri,'(/a,f12.6)') '  Sum = total charge:', sum(charges)

  if (qfit_multipole_rank >= 1) then
      write(lupri,'(/a)') " Dipoles:"
      write(lupri,'(a,6a,3a12)') "  ", "      ", "        X   ", &
       "        Y   ", "        Z   "
      do m = 1, size(charges)
          write(lupri,'(a,a6,3f12.6)') '@ ',namdep(m), (dipoles((m-1)*3+k),k=1,3)
      enddo
  endif
  if (qfit_multipole_rank >= 2) then
      write(lupri,'(/a)') " Quadrupoles:"
      write(lupri,'(a,6a,6a12)') "  ", "      ", "       XX   ", &
       "       XY   ", "       XZ   ", "       YY   ", "       YZ  ", &
       "       ZZ  "
       do m = 1, size(charges)
           write(lupri,'(a,a6,6f12.6)') '@ ',namdep(m), (quadrupoles((m-1)*6+k),k=1,6)
       enddo
  endif
  write(lupri,'(/a,f12.6)') '@ RMSD of fit', rmsd_value

  !call qfit_finalize
  deallocate( quadrupoles )
  deallocate( dipoles )
  deallocate( charges )
  call qexit('qfitlib_ifc_results')
end subroutine qfitlib_ifc_results

#if defined(VAR_MPI)
subroutine qfitlib_ifc_slave(lupri, iprint)
  use qfit, only : qfit_fit
#include "maxorb.h"
#include "infpar.h"
! parameter values
  integer :: lupri, iprint
! local variables
  integer :: dens_size
  real*8, dimension(:), allocatable :: dens
  call qenter('qfitlib_ifc_slave')
  ! sync density so we can call qfit_fit
  call mpixbcast(dens_size, 1, 'INTEGER', master)
  allocate(dens(dens_size))
  call mpixbcast(dens, dens_size, 'DOUBLE', master)
  call qfit_fit(dens, lupri=lupri)
  deallocate(dens)
  flush(lupri)
  call qexit('qfitlib_ifc_slave')
end subroutine qfitlib_ifc_slave
#endif

end module qfitlib_interface

#else

module qfitlib_interface
  implicit none

  private :: qfitlib_ifc_compute_dipole

  public :: qfitlib_ifc_input_reader
  public :: qfitlib_ifc_initialize
  public :: qfitlib_ifc_finalize
  public :: qfitlib_ifc_information
  public :: qfitlib_ifc_fit
  public :: qfitlib_ifc_results

  public :: qfitlib_ifc_sirius_fit
  public :: qfitlib_ifc_response_fit

#if defined(VAR_MPI)
  public :: qfitlib_ifc_slave ! entry-point for parallel nodes
#endif

contains

subroutine qfitlib_ifc_input_reader(word)
  character(len=7), intent(inout) :: word
  call qenter('qfitlib_ifc_input_reader')
  call quit('QFITLIB not enabled - using QFITLIB dummy routines.')
  call qexit('qfitlib_ifc_input_reader')
end subroutine qfitlib_ifc_input_reader

subroutine qfitlib_ifc_initialize
  call qenter('qfitlib_ifc_initialize')
  call quit('QFITLIB not enabled - using QFITLIB dummy routines.')
  call qexit('qfitlib_ifc_initialize')
end subroutine qfitlib_ifc_initialize

subroutine qfitlib_ifc_finalize
  call qenter('qfitlib_ifc_finalize')
  call quit('QFITLIB not enabled - using QFITLIB dummy routines.')
  call qexit('qfitlib_ifc_finalize')
end subroutine qfitlib_ifc_finalize

subroutine qfitlib_ifc_information
  call qenter('qfitlib_ifc_information')
  call quit('QFITLIB not enabled - using QFITLIB dummy routines.')
  call qexit('qfitlib_ifc_information')
end subroutine qfitlib_ifc_information

subroutine qfitlib_ifc_fit(dens)
  real*8, dimension(:), intent(in) :: dens
  call qenter('qfitlib_ifc_fit')
  call quit('QFITLIB not enabled - using QFITLIB dummy routines.')
  call qexit('qfitlib_ifc_fit')
end subroutine qfitlib_ifc_fit

!------------------------------------------------------------------------------
subroutine qfitlib_ifc_sirius_fit(cmo, wrk, nwrk)
  integer :: nwrk
  real*8, dimension(*), intent(in) :: cmo
  real*8, dimension(nwrk), intent(in) :: wrk
  call qenter('qfitlib_ifc_sirius_fit')
  call quit('QFITLIB not enabled - using QFITLIB dummy routines.')
  call qexit('qfitlib_ifc_sirius_fit')
end subroutine qfitlib_ifc_sirius_fit

!------------------------------------------------------------------------------
subroutine qfitlib_ifc_response_fit(cmo, wop, udv, wrk, nwrk)
  integer :: nwrk
  real*8, dimension(*), intent(in) :: cmo, wop, udv
  real*8, dimension(nwrk), intent(in) :: wrk
  call qenter('qfitlib_ifc_response_fit')
  call quit('QFITLIB not enabled - using QFITLIB dummy routines.')
  call qexit('qfitlib_ifc_response_fit')
end subroutine qfitlib_ifc_response_fit

!------------------------------------------------------------------------------
subroutine qfitlib_ifc_compute_dipole(dcao)
  real*8, dimension(:) :: dcao
  call qenter('qfitlib_ifc_compute_dipole')
  call quit('QFITLIB not enabled - using QFITLIB dummy routines.')
  call qexit('qfitlib_ifc_compute_dipole')
end subroutine qfitlib_ifc_compute_dipole

subroutine qfitlib_ifc_results
  call qenter('qfitlib_ifc_results')
  call quit('QFITLIB not enabled - using QFITLIB dummy routines.')
  call qexit('qfitlib_ifc_results')
end subroutine qfitlib_ifc_results

#if defined(VAR_MPI)
subroutine qfitlib_ifc_slave(lupri, iprint)
  integer :: lupri, iprint
  call qenter('qfitlib_ifc_slave')
  call quit('QFITLIB not enabled - using QFITLIB dummy routines.')
  call qexit('qfitlib_ifc_slave')
end subroutine qfitlib_ifc_slave
#endif

end module qfitlib_interface

#endif
