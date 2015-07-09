
!***************************************************************************
!*                                                                         *
!*  MAIN ROUTINES TO COMPUTE GRADIENT,HESSIAN ON TRIAL VECTOR AND PRECOND  *
!*                                                                         *
!* Implemented for                                                         *
!* - Pipek-Mezey with external power m, Lowdin charge only                 *
!* - Generalized Pipek-Mezey, external power m, Lowdin charge              *
!*                                                                         *
!* Lowdin orbitals may be optimized using orbspread scheme                 *
!*                                                                         *
!* NOTE: For Pipek-Mezey --> .PML m=2, Lowdin pop.analysis                 *
!* NOTE: For Pipek-Mezey --> .PMM m=2, Mulliken pop.analysis               *
!*                                                                         *
!* Author: Ida-Marie Hoeyvik                                                *
!***************************************************************************


module loc_types
use matrix_module
use typedef
use typedeftype

TYPE PMitem
! ******************************************
! *         Essential Settings              *
! ******************************************
! Exponent m
real(realk) :: m
! Use Pipek-Mezey scheme with Lowdin charge
logical :: PipekMezeyLowdin
! Use Pipek-Mezey scheme with Lowdin charge
logical :: PipekMezeyMull
! use linesearch in optimization
logical :: linesearch
! if using precond
logical :: precond


! ******************************************
! * Used for construction of LT and/or grad *
! ******************************************
! SU = S_sqrt for Lowdin orbitals
! SU = S for Mulliken based Pipek measure  
type(matrix) :: SU
!SU*Coeff
type(matrix) :: SC
!CMOS
type(matrix),pointer :: CMO
! number og atoms in molecule
integer :: natoms
! number of basis functs to use
integer :: nbas
! number of orbitals in mol (CMO%ncol)
integer :: norb
! d_i^m
real(realk),pointer :: di(:)
! d_i^(m+1)
real(realk),pointer :: d_m(:)
! d_i^(m+2)
real(realk),pointer :: d_m2(:)
! Pos(A); start pos. basisfunct on atom A  
integer,pointer :: Pos(:)
! nbasA(A); # of basisfunct on atom A
integer,pointer :: nbasA(:)
! Q(i,A) --> Q_{ii}^A
real(realk), pointer :: Q(:,:)
! Qkappa(i,A)-->[Q^A Kappa]_{ii}
real(realk),pointer :: Qkappa(:,:)


! ******************************************
! *       Used in localization scheme      *
! ******************************************
! current function value
real(realk) :: funcVal
! value indicating degree of localization
real(realk) :: loc_degree
! Matrix of diagonal Hessian elements
type(matrix),pointer :: P


! ******************************************
! *       Useful output etc.               *
! ******************************************
! Array that stores sorted lowdin charges
real(realk),pointer :: SortedLow(:,:)
! extra print statements (debug purposes)
logical :: orb_debug

END TYPE PMitem

contains

subroutine loc_types_dummy()

end subroutine loc_types_dummy

end module loc_types


