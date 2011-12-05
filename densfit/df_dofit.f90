! 
!...   Copyright (c) 2004 by the authors of Dirac (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   Dirac, a relativistic ab initio electronic structure program,
!...   release DIRAC04.0 (2004),
!...   written by H. J. Aa. Jensen, T. Saue, and L. Visscher
!...   with contributions from V. Bakken, E. Eliav, T. Enevoldsen, T. Fleig,
!...   O. Fossgaard, T. Helgaker, J. Henriksson, J. K. Laerdahl, C. V. Larsen,
!...   P. Norman, J. Olsen, M. Pernpointner, J. K. Pedersen, K. Ruud,
!...   P. Salek, J. N. P. van Stralen, J. Thyssen, O. Visser, and T. Winther
!...   (http://dirac.chem.sdu.dk).
!...
!...   For a suitable BibTEX entry, see:
!...      http://dirac.chem.sdu.dk/doc/reference.shtml
!...
!...   This source code is provided under a written licence and may be
!...   used, copied, transmitted, or stored only in accord with that
!...   written licence.
!...
!...   In particular, no part of the source code or compiled modules may
!...   be distributed outside the research group of the licence holder.
!...   This means also that persons (e.g. post-docs) leaving the research
!...   group of the licence holder may not take any part of Dirac,
!...   including modified files, with him/her, unless that person has
!...   obtained his/her own licence.
!...
!...   For questions concerning this copyright write to:
!...      dirac-admin@dirac.chem.sdu.dk
!...
!...   For information on how to get a licence see:
!...      http://dirac.chem.sdu.dk
! 
module Varlengths

! ----------------------
! define KIND parameters
! ----------------------

integer, parameter :: KINT  = 4
integer, parameter :: KREAL = 8

end module Varlengths

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

module densfit_common

use varlengths

! ---------------------------------
! define data structure for centers
! ---------------------------------

type center_definition
  integer(KINT) :: n_basis, n_fit
end type

end module densfit_common

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
 
subroutine FitPairDensity ( index_A, index_B, &
                            Density_AB, ld_DM, &
                            fit_coefficients_A, fit_coefficients_B, &
                            ierr)

use Varlengths
use densfit_common

!
! ==== Description ===>
!
! Calculate fit coefficients given the MO-coefficients for the atom pair.
! Written by L. Visscher, november 2003.
!
! ==== Dummy arguments ===>
!
! General information
! Index_X : Identification of nucleus X so that we can find its characteristics
!           by calling the appropriate routine for the code that we're using. What we
!           in any case need to know are the number of basis and fit functions on the center.
!
! Information about the density
! Density_AB : density matrix block for atom pair AB
! ld_DM : leading dimension of density matrix
!
! Result
! fit_coefficients_X on input  : value of the fit coefficients without pair density contribution
!                    on output : value of the fit coefficients with pair density contribution
! ierr                         : error code (0:OK, 1:error computing weight matrix, 2:error in fit)

  integer(KINT), intent(in)  :: index_A, index_B, ld_DM
  real(KREAL), intent(in)    :: Density_AB(ld_dM,*)
  real(KREAL), intent(inout) :: fit_coefficients_A(*), fit_coefficients_B(*)
  integer(KINT), intent(out) :: ierr
!
! ==== Local variables ===>
!
  integer(KINT)                 :: n_fit, n_fit_A, n_fit_B
  real(KREAL), allocatable      :: fit_coefficients(:), weight_matrix(:,:), weight_vector(:)
  type(center_definition)       :: definition_center_A, definition_center_B
!
! ==== Executable code ===>
!
  ierr = 0
!
! Get information about the centers
!
  call GetCenterInfo (index_A, definition_center_A)
  call GetCenterInfo (index_B, definition_center_B)
!
! All fit functions from A and B may contribute
!
  n_fit_A = definition_center_A % n_fit
  n_fit_B = definition_center_B % n_fit
  n_fit = n_fit_A + n_fit_B
! 
  allocate ( weight_matrix(n_fit, n_fit) )
  allocate ( weight_vector(n_fit) )
!
! Compute weight matrix and vector for the density fit. This may be based on Coulomb (V) or
! overlap (S) matrix or possibly something else. This is decided inside GetWeights,
! we only need to know that it's positive definite here.
!
  call GetWeights (definition_center_A, definition_center_B, &
                   Density_AB, ld_DM, &
                   weight_matrix, weight_vector, n_fit, ierr)
!
! Do the density fit
!
  allocate ( fit_coefficients(n_fit) )
  call ComputeFitCoefficients (n_fit, weight_matrix, n_fit, &
                               weight_vector, fit_coefficients, info)
! 
  if (info.eq.0) then
    fit_coefficients_A(1:n_fit_A) = fit_coefficients_A(1:n_fit_A) + fit_coefficients(1:n_fit_A)
    fit_coefficients_B(1:n_fit_B) = fit_coefficients_B(1:n_fit_B) + fit_coefficients(n_fit_A+1:n_fit)
  else
    print*,"error in density fit routine ComputeFitCoefficients",info
    ierr = 2
  end if
!
  deallocate (fit_coefficients)
  deallocate (weight_vector)
  deallocate (weight_matrix)
  return
!
end subroutine FitPairDensity
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine ComputeFitCoefficients (n_fit, weight_matrix, ld_weight_matrix, &
                                   weight_vector,  fit_coefficients, info)
! ==== Description ===>
!
! Calculate fit coefficients for density fit.
! Written by L. Visscher, november 2003.
!
! ==== Dummy arguments ===>
!

use Varlengths
use densfit_common

  integer(KINT), intent(in)  :: n_fit, ld_weight_matrix
  integer(KINT), intent(out) :: info
  real(KREAL), intent(in)    :: weight_matrix(ld_weight_matrix, n_fit), &
                                weight_vector(n_fit)
  real(KREAL), intent(out)   :: fit_coefficients(n_fit)
!
! ==== Local variables ===>
!
  real(KREAL), allocatable :: work_matrix(:,:)
!
! ==== Implementation based on the simple LAPACK driver available in ATLAS ===>
!
  allocate ( work_matrix(n_fit, n_fit) )
!
! -----------------------------------------------------
! Our matrix should be positive definite and symmetric.
! We call a LAPACK driver that solves Ax=B for this case.
! Since it will overwrite its input data we need to make 
! a copy first. If the procedure becomes
! numerically unstable we might want to change to the
! dposvx driver that gives more flexibility and error control, 
! but requires more memory and is not included in ATLAS.
!
! See also http://www.netlib.org/lapack
! -----------------------------------------------------
!
! Copy the data to be used in dposv
!
  work_matrix(1:n_fit,1:n_fit) = weight_matrix(1:n_fit,1:n_fit)
  fit_coefficients(1:n_fit) = weight_vector(1:n_fit)
!
! Perform the fit
!
!  call dposv ( 'L', n_fit, 1, work_matrix, n_fit, &
!              fit_coefficients, n_fit, info )
!
  deallocate ( work_matrix )
! 
  return
end subroutine ComputeFitCoefficients
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
 
subroutine GetCenterInfo (index, definition_center)

use Varlengths
use densfit_common

!
! ==== Description ===>
!
! Get relevant data about an expansion center.
! This is a dummy routine, used to test the code !!!!!!
! Written by L. Visscher, november 2003.
!
! ==== Dummy arguments ===>
!
! Index : Identification of the nucleus.
! definition_center  : Structure with all the relevant information
!
  integer(KINT), intent(in)                  :: index
  type(center_definition), intent(out)       :: definition_center
!
  definition_center % n_fit = 2
  return
end subroutine GetCenterInfo
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
 
subroutine GetWeights (definition_center_A, definition_center_B, &
                       Density_AB, ld_DM, &
                       weight_matrix, weight_vector, n_fit, ierr)

use Varlengths
use densfit_common

!
! ==== Description ===>
!
! Get the matrices needed for the density fitting.
! This is a dummy routine, used to test the code !!!!!!
!
! Written by L. Visscher, november 2003.
!
! ==== Dummy arguments ===>
!
! Index : Identification of the nucleus.
! definition_center  : Structure with all the relevant information
!
  integer(KINT), intent(out)                :: ierr
  type(center_definition), intent(in)       :: definition_center_A, definition_center_B
  real(KREAL), intent(in)                   :: Density_AB(ld_DM,*)
  real(KREAL), intent(out)                  :: weight_matrix(n_fit,n_fit), weight_vector(n_fit)
!
  ierr = 0
  do i = 1, n_fit
     weight_vector(i) = i * 0.1D0 
  end do
  do i = 1, n_fit
     do j = 1, n_fit
        if (i.eq.j) then
           weight_matrix(i,j) = 1.d0
        else
           weight_matrix(i,j) = 0.d0
        endif
     end do
  end do
!
  return
end subroutine GetWeights
 
