!> @file
!> Contains the module soeo_util
!
!> brief Routines used by soeo
!> author C. Nygaard
!> date 2010 - 2011
!
!> Routines written in 2009 and 2010 and collected in this module 2011.02.02
module soeo_util

use precision
use matrix_module
use Matrix_operations
use Matrix_operations_aux
use soeo_typedef
use soeo_matop
use matrix_util
use dal_interface
use memory_handling

contains
!  soeo_Dmo_init
!  soeo_oldtheta_init
!  soeo_getnew_theta
!  soeo_getnew_Dmo
!  soeo_getnew_Dao
!  soeo_getnew_C
!  soeo_get_nfirst
!  soeo_get_nsecond
!  soeo_get_gradient
!  soeo_linear_transform
!  soeo_get_GmoL
!  soeo_get_Epred
!  soeo_transformplus
!  soeo_gradplus


!> \brief Finds the MO density from the AO density
!> \author C. Nygaard
!> \date 2010
!> \param S The overlap matrix
!> \param C The MO orbital coefficients
!> \param Dao The density matrix in the AO basis
!> \param Dmo The density matrix in the MO basis
!> \param unres True if alculation is unrestricted
!
!> Dmo = C^T S Dao S C
!
!=======================================================================
subroutine soeo_Dmo_init (S, C, Dao, Dmo, unres)
implicit none
!I/O
type(matrix), intent(inout)  :: Dmo
type(matrix), intent(in)     :: S, C, Dao
logical, intent(in)          :: unres
!Other
integer                      :: Nbast, i, j
real(realk), pointer         :: tmp(:)

if (unres) then
  call mem_alloc (tmp,2)
else
  call mem_alloc (tmp,1)
endif

Nbast = Dmo%nrow
call mat_zero (Dmo)

call util_AO_to_MO_2 (S, C, Dao, Dmo, .false.) !Dmo = C^T S Dao S C

!Cleaning up so that zero elements are really zero
do i=1,Nbast
  do j=1,Nbast
    call mat_get_ab_elms (Dmo, i, j, tmp)
    if (tmp(1) < 1.0E-5_realk) then
      tmp(1) = 0.0E0_realk
    endif
    if (unres) then
      if (tmp(2) < 1.0E-5_realk) then
        tmp(2) = 0.0E0_realk
      endif
    endif
    call mat_create_ab_elms (i, j, tmp, Dmo)
  enddo
enddo

call mem_dealloc (tmp)

end subroutine soeo_Dmo_init
!=======================================================================

!> \brief Initializes the vector oldtheta from Dmo
!> \author C. Nygaard
!> \date 2010
!> \param Dmo The density matrix in the MO basis
!> \param oldtheta The starting point for the occupation angles
!> \param unres True if calculation is unrestricted
!
!>  oldtheta(i) = acos(sqrt(Dmo(i,i)))
!>  because
!>  Dmo(i,i) = cos^2(oldtheta(i))
!
!=======================================================================
subroutine soeo_oldtheta_init (Dmo, oldtheta, unres)

implicit none

!I/O:
type(matrix), intent(in)     :: Dmo
type(matrix), intent(inout)  :: oldtheta
logical, intent(in)          :: unres
!Other:
integer                      :: i, Nbast
real(realk), pointer         :: tmp(:)

Nbast = Dmo%nrow

if (unres) then
  call mem_alloc (tmp,2)
else
  call mem_alloc (tmp,1)
endif

do i=1,Nbast
  call mat_get_ab_elms (Dmo, i, i, tmp)

  !clean numbers
  ! acos(n) = NaN if n > 1 or n < 0
  if (tmp(1) > 1.0E0_realk) then 
    tmp(1) = 1.0E0_realk
  endif
  if (tmp(1) < 0.0E0_realk) then
    tmp(1) = 0.0E0_realk
  endif
  tmp(1) = sqrt(tmp(1))
  tmp(1) = acos(tmp(1))

  if (unres) then !beta part
    if (tmp(2) > 1.0E0_realk) then 
      tmp(2) = 1.0E0_realk
    endif
    if (tmp(2) < 0.0E0_realk) then
      tmp(2) = 0.0E0_realk
    endif
    tmp(2) = sqrt(tmp(2))
    tmp(2) = acos(tmp(2))
  endif

  call mat_create_ab_elms (i, 1, tmp, oldtheta)
enddo

call mem_dealloc (tmp)

end subroutine soeo_oldtheta_init
!=======================================================================


!> \brief Updates the occupation angles
!> \author C. Nygaard
!> \date 2010
!> \param deltatheta The change in occupation angles
!> \param space Information about the active space
!> \param oldtheta Starting point for the occupation numbers
!> \param unres True if calculation is unrestricted
!
!> theta = oldtheta + deltatheta
!
!=======================================================================
subroutine soeo_getnew_theta (deltatheta, space, oldtheta, unres)
implicit none
!I/O
type(matrix), intent(inout)      :: oldtheta
type(matrix), intent(in)         :: deltatheta
type(soeoItem_space), intent(in) :: space
logical, intent(in)              :: unres


!Other
integer                      :: i
real(realk), pointer         :: old(:), new(:)

if (unres) then
  call mem_alloc (old,2)
  call mem_alloc (new,2)
else
  call mem_alloc (old,1)
  call mem_alloc (new,1)
endif

!only change theta in the active space
do i=space%Nocc+1,space%Nocc+space%Nact
  call mat_get_ab_elms (oldtheta, i, 1, old)
  call mat_get_ab_elms (deltatheta, i-space%Nocc, 1, new)
  old = old + new
  call mat_create_ab_elms (i, 1, old, oldtheta)
enddo

call mem_dealloc (old)
call mem_dealloc (new)

end subroutine soeo_getnew_theta
!=======================================================================

!> \brief Updates the MO density matrix
!> \author C. Nygaard
!> \date 2010
!> \param theta The occupation angles
!> \param Dmo The density matrix in the MO basis
!> \param unres True if calculation is unrestricted
!
!>  Dmo(i,j) = delta_ij cos^2 (oldtheta(i))
!
!=======================================================================
subroutine soeo_getnew_Dmo (theta, Dmo, unres)

implicit none

!I/O
type(matrix), intent(inout)  :: theta
type(matrix), intent(inout)  :: Dmo
logical, intent(in)          :: unres
!Other
integer                      :: Nbast
integer                      :: nact, s, i, j
real(realk), pointer         :: n(:), t(:)
real(realk)                  :: thresh, dn

thresh = 1.0E-5_realk

call mat_zero (Dmo)
if (unres) then
  s = 2
else
  s = 1
endif
call mem_alloc (n,s)
call mem_alloc (t,s)

Nbast = Dmo%nrow

dn = 0.0E0_realk
nact = 0
do j = 1, s
  do i=1,Nbast
  
    call mat_get_ab_elms (theta, i, 1, t)
    n = cos(t) * cos(t)
  
    !cleaning
    if (n(s) < thresh) then
       dn = dn + n(s)
       n(s) = 0.0E0_realk
       t(s) = dacos(n(s))
    elseif (1.0E0_realk-n(s) < thresh) then
       dn = dn - (1.0E0_realk-n(s))
       n(s) = 1.0E0_realk
       t(s) = dacos(n(s))
    else
       nact = nact + 1
    endif
    call mat_create_ab_elms (i, i, n, Dmo)
    call mat_create_ab_elms (i, 1, t, theta)
  enddo


  !finish cleaning
  dn = dn / nact
  do i=1,Nbast
    call mat_get_ab_elms (Dmo, i, i, n)
    if ((n(s) > thresh) .and. (1.0E0_realk-n(s) > thresh)) then
      n(s) = n(s) + dn
      t(s) = dacos(dsqrt(n(s)))
      call mat_create_ab_elms (i, i, n, Dmo)
      call mat_create_ab_elms (i, 1, t, theta)
    endif
  enddo

enddo

call mem_dealloc (n)
call mem_dealloc (t)

end subroutine soeo_getnew_Dmo
!=======================================================================



!> \brief Finds the new MO-coefficients
!> \author C. Nygaard
!> \date Sep. 1. 2010
!> \param C The MO orbital coefficients (inout)
!> \param Kin The MO orbital rotation matrix
!
!> C(i+1) = C(i)exp(-K)
!>        = C(i)(1 - K + 1/2 K^2 - 1/6 K^3 + ...)
!=======================================================================
subroutine soeo_getnew_C (Kin, C)

implicit none

type(matrix), intent(in)    :: Kin
type(matrix), intent(inout) :: C
type(matrix)                :: K, Cstep, tmp
integer                     :: i
real(realk)                 :: fac, n, thresh, stepnorm

call mat_init (K, Kin%nrow, Kin%ncol)
call mat_init (Cstep, C%nrow, C%ncol)
call mat_init (tmp, C%nrow, C%ncol)

call mat_assign(K,Kin)
call mat_assign(Cstep,C)
i = 1
thresh = 1.0E-15_realk
do
  !The factor in front:
  !---------------------------------------------------------------------
  if (i == 1) then
    fac = 1.0E0_realk
  else
    fac = fac * i
  endif
  i = i + 1
  n = 1.0E0_realk / fac
  !---------------------------------------------------------------------

  !This terms product:
  !---------------------------------------------------------------------
  call mat_mul (Cstep, K, 'n', 'n', -1.0E0_realk, 0.0E0_realk, tmp)
  call mat_assign(Cstep,tmp)
  !---------------------------------------------------------------------

  !Add the term:
  !---------------------------------------------------------------------
  call mat_daxpy (n, Cstep, C)
  !---------------------------------------------------------------------

  !Convergence check:
  !---------------------------------------------------------------------
  call mat_scal (n, tmp)
  stepnorm = mat_sqnorm2 (tmp)
  stepnorm = sqrt(stepnorm)
  if (stepnorm < thresh) exit
  !---------------------------------------------------------------------
enddo

call mat_free (K)
call mat_free (Cstep)
call mat_free (tmp)

end subroutine soeo_getnew_C
!=======================================================================



!> \brief Calculates first derivative of Dmo wrt the occupation angles
!> \author C. Nygaard
!> \date 2010
!> \param oldtheta point where derivative is calculated
!> \param p Which occupation angle to derive wrt
!> \param part 'a' for alpha-part and 'b' for beta-part
!> \param nfirst The derivative
!> \param unres True if calculation is unrestricted
!
!>  nfirst(p)_ij = -2 * delta_ij * delta_ip
!>                    * sin(oldtheta(p)) * cos(oldtheta(p))
!
!=======================================================================
subroutine soeo_get_nfirst (oldtheta, p, part, nfirst, unres)

implicit none

!I/O:
integer, intent(in)          :: p
type(matrix), intent(in)     :: oldtheta
type(matrix), intent(inout)  :: nfirst
character, intent(in)        :: part !part = 'a' or 'A' for alpha
                                     !     = 'b' or 'B' for beta
logical, intent(in)          :: unres

!Other:
real(realk), pointer         :: tmp1(:), tmp2(:)

if (unres) then
  call mem_alloc (tmp1,2)
  call mem_alloc (tmp2,2)
else
  call mem_alloc (tmp1,1)
  call mem_alloc (tmp2,1)
endif

call mat_zero (nfirst)

call mat_get_ab_elms (oldtheta, p, 1, tmp1)

tmp2 = 0.0E0_realk
if (part=='a' .or. part=='A') then !alpha-part of the matrix - elms(:)
  tmp2(1) = -2.0E0_realk * sin(tmp1(1)) * cos(tmp1(1))
elseif (part=='b' .or. part=='B') then !beta-part of the matrix - elmsb(:)
  if (unres) then
    tmp2(2) = -2.0E0_realk * sin(tmp1(2)) * cos(tmp1(2))
  else
    call lsquit ('Beta-part of nfirst doesnt exist in restricted calculation',-1)
  endif
else !part = alpha or beta, otherways you're stupid
  call lsquit ('Wrong argument (part) in soeo_get_nfirst',-1)
endif

call mat_create_ab_elms (p, p, tmp2, nfirst)

call mem_dealloc (tmp1)
call mem_dealloc (tmp2)

end subroutine soeo_get_nfirst
!=======================================================================



!> \brief Calculates first derivative of Dmo wrt the occupation angles
!> \author C. Nygaard
!> \date 2010
!> \param oldtheta point where derivative is calculated
!> \param p Which occupation angle to derive wrt
!> \param part 'a' for alpha-part and 'b' for beta-part
!> \param nfirst The derivative
!> \param unres True if calculation is unrestricted
!
!>  nfirst(p)_ij = -2 * delta_ij * delta_ip
!>                    * sin(oldtheta(p)) * cos(oldtheta(p))
!
!=======================================================================
subroutine soeo_get_nfirst2 (space, oldtheta, p, nfirst, unres)

implicit none

!I/O:
integer, intent(in)              :: p
type(matrix), intent(in)         :: oldtheta
type(soeoItem_space), intent(in) :: space
type(matrix), intent(inout)      :: nfirst
logical, intent(in)              :: unres

!Other:
real(realk), pointer             :: tmp1(:), tmp2(:)
integer                          :: elm, part

if (p<=space%Nact) then
  elm = space%Nocc+p
  part = 1
else
  elm = space%Nocc+p-space%Nact
  part = 2
endif

if (unres) then
  call mem_alloc (tmp1,2)
  call mem_alloc (tmp2,2)
else
  call mem_alloc (tmp1,1)
  call mem_alloc (tmp2,1)
endif

call mat_zero (nfirst)

call mat_get_ab_elms (oldtheta, elm, 1, tmp1)

tmp2 = 0.0E0_realk
if (part==1) then !alpha-part of the matrix - elms(:)
  tmp2(1) = -2.0E0_realk * sin(tmp1(1)) * cos(tmp1(1))
elseif (part==2) then !beta-part of the matrix - elmsb(:)
  if (unres) then
    tmp2(2) = -2.0E0_realk * sin(tmp1(2)) * cos(tmp1(2))
  else
    call lsquit ('Beta-part of nfirst doesnt exist in restricted calculation',-1)
  endif
endif

call mat_create_ab_elms (elm, elm, tmp2, nfirst)

call mem_dealloc (tmp1)
call mem_dealloc (tmp2)

end subroutine soeo_get_nfirst2
!=======================================================================



!> \brief Calculates the second derivative of Dmo wrt the occupation angles
!> \author C. Nygaard
!> \date 2010
!> \param oldtheta Where the derivative should be calculated
!> \param p Which occupation angle to make first derivation wrt
!> \param q Which occupation angle so make second edrivation wrt
!> \param partp 'a' if p is in the alpha-part, 'b' for beta-part
!> \param partq 'a' if q is in the alpha-part, 'b' for beta-part
!> \param nsecond The derivative
!> \param unres True if calculation is unrestricted
!
!>  nsecond(i,j)_pq = -2 * delta_ij * delta_ip * delta_pq
!>                       * (cos^2 (oldtheta(q)) - sin^2 (oldtheta(q)))
!
!=======================================================================
subroutine soeo_get_nsecond (oldtheta, p, q, partp, partq, nsecond, unres)

implicit none

!Input/output
integer, intent(in)          :: p, q
type(matrix), intent(in)     :: oldtheta
type(matrix), intent(inout)  :: nsecond
character, intent(in)        :: partp, partq !part = 'a' or 'A' for alpha
                                             !     = 'b' or 'B' for beta
logical, intent(in)          :: unres

!Other:
real(realk), pointer         :: tmp1(:), tmp2(:)

if (unres) then
  call mem_alloc (tmp1,2)
  call mem_alloc (tmp2,2)
else
  call mem_alloc (tmp1,1)
  call mem_alloc (tmp2,1)
endif

call mat_zero (nsecond)

if (partp == partq) then
  if (p == q) then
    call mat_get_ab_elms (oldtheta, q, 1, tmp1)

    tmp2 = 0.0E0_realk
    if (partq == 'a' .or. partq == 'A') then !alpha part
      tmp2(1) = -2.0E0_realk * cos(tmp1(1)) * cos(tmp1(1))
      tmp2(1) = tmp2(1) + 2.0E0_realk * sin(tmp1(1)) * sin(tmp1(1))
    elseif (partq == 'b' .or. partq == 'B') then !beta part
      if (unres) then
        tmp2(2) = -2.0E0_realk * cos(tmp1(2)) * cos(tmp1(2))
        tmp2(2) = tmp2(2) + 2.0E0_realk * sin(tmp1(2)) * sin(tmp1(2))
      else
        call lsquit ('Beta-part of nsecond doesnt exist in restricted calculation',-1)
      endif
    else
      call lsquit ('Wrong argument (part) in soeo_get_nsecond',-1)
    endif

    call mat_create_ab_elms (p, p, tmp2, nsecond)
  endif
endif

call mem_dealloc (tmp1)
call mem_dealloc (tmp2)

end subroutine soeo_get_nsecond
!=======================================================================



!> \brief Calculates the second derivative of Dmo wrt the occupation angles
!> \author C. Nygaard
!> \date 2010
!> \param oldtheta Where the derivative should be calculated
!> \param p Which occupation angle to make first derivation wrt
!> \param q Which occupation angle so make second edrivation wrt
!> \param partp 'a' if p is in the alpha-part, 'b' for beta-part
!> \param partq 'a' if q is in the alpha-part, 'b' for beta-part
!> \param nsecond The derivative
!> \param unres True if calculation is unrestricted
!
!>  nsecond(i,j)_pq = -2 * delta_ij * delta_ip * delta_pq
!>                       * (cos^2 (oldtheta(q)) - sin^2 (oldtheta(q)))
!
!=======================================================================
subroutine soeo_get_nsecond2 (space, oldtheta, p, q, nsecond, unres)

implicit none

!Input/output
type(soeoItem_space), intent(in) :: space
integer, intent(in)              :: p, q
type(matrix), intent(in)         :: oldtheta
type(matrix), intent(inout)      :: nsecond
logical, intent(in)              :: unres

!Other:
real(realk), pointer         :: tmp1(:), tmp2(:)
integer                      :: elm, part

if (unres) then
  call mem_alloc (tmp1,2)
  call mem_alloc (tmp2,2)
else
  call mem_alloc (tmp1,1)
  call mem_alloc (tmp2,1)
endif

call mat_zero (nsecond)

if (p == q) then
  if (p<=space%Nact) then
    elm = space%Nocc+p
    part = 1
  else
    elm = space%Nocc+p-space%Nact
    part = 2
  endif

  call mat_get_ab_elms (oldtheta, elm, 1, tmp1)

  tmp2 = 0.0E0_realk
  if (part == 1) then !alpha part
    tmp2(1) = -2.0E0_realk * cos(tmp1(1)) * cos(tmp1(1))
    tmp2(1) = tmp2(1) + 2.0E0_realk * sin(tmp1(1)) * sin(tmp1(1))
  elseif (part == 2) then !beta part
    if (unres) then
      tmp2(2) = -2.0E0_realk * cos(tmp1(2)) * cos(tmp1(2))
      tmp2(2) = tmp2(2) + 2.0E0_realk * sin(tmp1(2)) * sin(tmp1(2))
    else
      call lsquit ('Beta-part of nsecond doesnt exist in restricted calculation',-1)
    endif
  endif

  call mat_create_ab_elms (elm, elm, tmp2, nsecond)
endif

call mem_dealloc (tmp1)
call mem_dealloc (tmp2)

end subroutine soeo_get_nsecond2
!=======================================================================



!> \brief Calculates left-hand-side (lhs) of the SOEO equatios
!> \author C. Nygaard
!> \date 2010
!> \param soeo Contains all matrices used in soeo
!> \param grad_m Matrix part of the lhs
!> \param grad_v Vector part of the lhs
!
!>  grad_m    = -4*[Fmo, Dmo]
!>  grad_v(p) = 2*Tr(Fmo*nfirst(p))
!
!=======================================================================
subroutine soeo_get_gradient (soeo, grad_m, grad_v)

implicit none

!Input/output:
type(soeoItem), intent(in)   :: soeo
type(matrix), intent(inout)  :: grad_m, grad_v
!Other:
integer                      :: i
real(realk), pointer         :: tmp(:)

if (soeo%cfg_unres) then
  call mem_alloc (tmp,2)
else
  call mem_alloc (tmp,1)
endif

call mat_zero (grad_m)
call mat_zero (grad_v)

!The matrix-part:
!-----------------------------------------------------------------------
call commutator (-2.0E0_realk, soeo%mats%Fmo, soeo%mats%Dmo, "n", "n", grad_m)
if (.not. soeo%cfg_unres) then
  call mat_scal (2.0E0_realk, grad_m)
endif
!-----------------------------------------------------------------------

!The vector-part:
!-----------------------------------------------------------------------
do i=1,soeo%space%Nact
  tmp(1) = mat_TrAB (soeo%mats%Fmo, soeo%mats%nfirst(i))
  if (soeo%cfg_unres) then
    tmp(2) = mat_TrAB (soeo%mats%Fmo, soeo%mats%nfirst(i+soeo%space%Nact))
  else
    tmp = 2.0E0_realk*tmp
  endif
  call mat_create_ab_elms (i, 1, tmp, grad_v)
enddo
!-----------------------------------------------------------------------

call mem_dealloc (tmp)

end subroutine soeo_get_gradient
!=======================================================================

!> \brief Same as soeo_get_gradient, but for properties
!> \author C. Nygaard
!> \date Nov. 20 2012
!> \param soeo Structure containing all matrices and info from soeo
!> \param Vao First order perturbation matrix, AO-basis
!> \param pgm Matrix-part of the lhs
!> \param pgv Vector-part of the lhs
!=======================================================================
subroutine soeo_get_propgrad (soeo, Vao, pgm, pgv)

implicit none

!Input/output:
type(soeoItem), intent(in)   :: soeo
type(matrix), intent(in)     :: Vao
type(matrix), intent(inout)  :: pgm, pgv
!Other:
type(matrix)                 :: Vmo
integer                      :: i
real(realk), pointer         :: tmp(:)

if (soeo%cfg_unres) then
  call mem_alloc (tmp,2)
else
  call mem_alloc (tmp,1)
endif

call mat_init (Vmo, soeo%space%Nbast, soeo%space%Nbast)
call mat_zero (Vmo)
call util_AO_to_MO_2 (soeo%mats%S, soeo%mats%C, Vao, Vmo, .true.) !Vmo=CT Vao C

call mat_zero (pgm)
call mat_zero (pgv)

!The matrix-part:
!-----------------------------------------------------------------------
call commutator (-2.0E0_realk, Vmo, soeo%mats%Dmo, "n", "n", pgm)
if (.not. soeo%cfg_unres) then
  call mat_scal (2.0E0_realk, pgm)
endif
!-----------------------------------------------------------------------

!The vector-part:
!-----------------------------------------------------------------------
do i=1,soeo%space%Nact
  tmp(1) = mat_TrAB (Vmo, soeo%mats%nfirst(i))
  if (soeo%cfg_unres) then
    tmp(2) = mat_TrAB (Vmo, soeo%mats%nfirst(i+soeo%space%Nact))
  else
    tmp = 2.0E0_realk*tmp
  endif
  call mat_create_ab_elms (i, 1, tmp, pgv)
enddo
!-----------------------------------------------------------------------

write (soeo%lupri, *) 'pgm'
call mat_print (pgm, 1, soeo%space%Nbast, 1, soeo%space%Nbast, soeo%lupri)
write (soeo%lupri, *) 'pgv'
call mat_print (pgv, 1, soeo%space%Nact, 1, 1, soeo%lupri)

call mem_dealloc (tmp)
call mat_free (Vmo)

end subroutine soeo_get_propgrad
!=======================================================================

!> \brief Calculates right-hand-side (rhs) of the SOEO equatios
!> \author C. Nygaard
!> \date 2010
!> \param soeo Contains all matrices used in soeo
!> \param b_m Matrix part of trial vector
!> \param b_v Vector part of trial vector
!> \param sigma_m Matrix part of the rhs
!> \param sigma_v Vector part of the rhs
!
!>  sigma_m    = -2*[[b_m,Fmo],Dmo] - 2*[Fmo,[Dmo,b_m]]
!>               - 4*sum_(r)(b_v(r)*[Fmo,nfirst(r)])
!>               - 8*[gmo*L,Dmo]
!>  sigma_v(p) = 2sum_(r)(b_v(r)*Tr(Fmo*nsecond(p,r)))
!>               - 2Tr([Fmo,b_m]nfirst(p))
!>               + 4*Tr(gmo*L*nfirst(p))
!>
!>  where L = [Dmo,b_m] + sum_(r)(b_v(r)*nfirst(r))
!
!> (nfirst(i) and nsecond(i,j) are matrices, not scalars)
!=======================================================================
subroutine soeo_linear_transform (soeo, b_m, b_v, sigma_m, sigma_v)

implicit none

!I/O
type(soeoItem), intent(in)   :: soeo
type(matrix), intent(in)     :: b_m, b_v
type(matrix), intent(inout)  :: sigma_m, sigma_v

!Other
integer                      :: Nbast, Nact, ndim, i, j
type(matrix)                 :: L, GmoL, commDb, sumn
!Temporary stuff:
type(matrix)                 :: tmp_m
real(realk)                  :: val1, val2
real(realk), pointer         :: tmp(:), tmp1(:)

!Initializations:
Nbast = soeo%space%Nbast
Nact  = soeo%space%Nact

if (soeo%cfg_unres) then
  call mem_alloc (tmp,2)
  call mem_alloc (tmp1,2)
else
  call mem_alloc (tmp,1)
  call mem_alloc (tmp1,1)
endif

call mat_init (L, Nbast, Nbast)
call mat_init (GmoL, Nbast, Nbast)
call mat_init (commDb, Nbast, Nbast)
call mat_init (sumn, Nbast, Nbast)
call mat_init (tmp_m, Nbast, Nbast)

ndim = size(soeo%mats%nfirst)

!Some usefull matrices
!-----------------------------------------------------------------------
call commutator (1.0E0_realk, soeo%mats%Dmo, b_m, "n", "n", commDb) ! [D,b_m]

call mat_zero (sumn)                                          ! sum_(r)(b_v(r) nfirst(r))
do i=1,Nact                                                   !
  call mat_get_ab_elms (b_v, i, 1, tmp)                       !
  call mat_ab_daxpy (tmp, soeo%mats%nfirst(i), sumn)          !
  if (soeo%cfg_unres) then                                    ! for unrestricted
    call mat_ab_daxpy (tmp, soeo%mats%nfirst(i+Nact), sumn)   !  the sum is also
  endif                                                       !  over spin
enddo                                                         !

!Creation of L and G(L)
call mat_add (1.0E0_realk, commDb, 1.0E0_realk, sumn, L)
call soeo_get_GmoL (L, soeo%mats%S, soeo%mats%C, soeo%mats%Dao, GmoL, soeo%lupri, soeo%do_dft)
!-----------------------------------------------------------------------

!!DEBUGGING
!!
!write (soeo%lupri, *)
!write (soeo%lupri, *) 'Input in soeo_linear_transform:'
!write (soeo%lupri, *) 'b_m ='
!call mat_print (b_m, 1, soeo%space%Nbast, 1, soeo%space%Nbast, soeo%lupri)
!write (soeo%lupri, *) 'b_v ='
!call mat_print (b_v, 1, soeo%space%Nact, 1, 1, soeo%lupri)
!write (soeo%lupri, *)
!write (soeo%lupri, *) 'Matrices in use'
!write (soeo%lupri, *) 'Dmo ='
!call mat_print (soeo%mats%Dmo, 1, soeo%space%Nbast, 1, soeo%space%Nbast, soeo%lupri)
!write (soeo%lupri, *) 'sumn ='
!call mat_print (sumn, 1, soeo%space%Nbast, 1, soeo%space%Nbast, soeo%lupri)
!write (soeo%lupri, *) 'commDb ='
!call mat_print (commDb, 1, soeo%space%Nbast, 1, soeo%space%Nbast, soeo%lupri)
!write (soeo%lupri, *) 'L ='
!call mat_print (L, 1, soeo%space%Nbast, 1, soeo%space%Nbast, soeo%lupri)
!write (soeo%lupri, *) 'GmoL ='
!call mat_print (GmoL, 1, soeo%space%Nbast, 1, soeo%space%Nbast, soeo%lupri)
!write (soeo%lupri, *)
!!
!!END DEBUGGING


!The matrix-part:
!-----------------------------------------------------------------------
!First term: sigma_m = 2*[[b_m,Fmo],Dmo]
! unres => *0.5
call commutator (1.0E0_realk, b_m, soeo%mats%Fmo, "n", "n", tmp_m)
call commutator (1.0E0_realk, tmp_m, soeo%mats%Dmo, "n", "n", sigma_m)
if (.not. soeo%cfg_unres) then
  call mat_scal (2.0E0_realk, sigma_m)
endif

!Second term: sigma_m = sigma_m + 2*[Fmo,[Dmo,b_m]]
! unres => *0.5
call commutator (1.0E0_realk, soeo%mats%Fmo, commDb, "n", "n", tmp_m)
if (.not. soeo%cfg_unres) then
  call mat_scal (2.0E0_realk, tmp_m)
endif
call mat_daxpy (1.0E0_realk, tmp_m, sigma_m)

!Third term: sigma_m = sigma_m + 4*sum_(r)(b_v(r)*[Fmo,nfirst(r)])
!                    = sigma_m + 4*[Fmo,sumn]
! unres => *0.5
call commutator (1.0E0_realk, soeo%mats%Fmo, sumn, "n", "n", tmp_m)
if (.not. soeo%cfg_unres) then
  call mat_scal (2.0E0_realk, tmp_m)
endif
call mat_daxpy (2.0E0_realk, tmp_m, sigma_m)

!Fourth term: sigma_m = sigma_m + 8*[GmoL,Dmo]
! unres => *0.25
call commutator (1.0E0_realk, GmoL, soeo%mats%Dmo, "n", "n", tmp_m)
if (.not. soeo%cfg_unres) then
  call mat_scal (4.0E0_realk, tmp_m)
endif
call mat_daxpy (2.0E0_realk, tmp_m, sigma_m)

call mat_scal (-1.0E0_realk, sigma_m)
!-----------------------------------------------------------------------

!The vector-part:
!-----------------------------------------------------------------------
call commutator (1.0E0_realk, soeo%mats%Fmo, b_m, "n", "n", tmp_m) ! [F,b_m]

!write (soeo%lupri, *)
!write (soeo%lupri, *) 'DEBUG LINEAR TRANSFORM'
!write (soeo%lupri, *) 'Nbast =', Nbast
!write (soeo%lupri, *) 'Nocc =', soeo%space%Nocc
!write (soeo%lupri, *) 'Nact =', Nact
!write (soeo%lupri, *) 'oldtheta ='
!call mat_print (soeo%mats%oldtheta, 1, Nbast, 1, 1, soeo%lupri)
!write (soeo%lupri, *) 'b_v ='
!call mat_print (b_v, 1, Nact, 1, 1, soeo%lupri)
!write (soeo%lupri, *) '[F,b_m] ='
!call mat_print (tmp_m, 1, Nbast, 1, Nbast, soeo%lupri)

do i=1,ndim

!write (soeo%lupri, *) '------------------------------------'
!write (soeo%lupri, *) 'i =', i

  tmp = 0.0E0_realk
  call mat_zero (sumn)                                  ! sumn = sum_(r)(b_v(r) nsecond(p,r))
  if (i <= Nact) then                                   !      = b_v(p)nsecond(p,p)
    call mat_get_ab_elms (b_v, i, 1, tmp)               !
  else                                                  !
    call mat_get_ab_elms (b_v, i-Nact, 1, tmp)          ! (also sum over spin)
  endif                                                 !
  call mat_ab_daxpy (tmp, soeo%mats%nsecond(i,i), sumn) !

!write (soeo%lupri, *) 'b_v(i) =', tmp
!write (soeo%lupri, *) 'nsecond(i,i) ='
!call mat_print (soeo%mats%nsecond(i,i), 1, Nbast, 1, Nbast, soeo%lupri)
!write (soeo%lupri, *) 'sumn ='
!call mat_print (sumn, 1, Nbast, 1, Nbast, soeo%lupri)

  !First term: sigma_v(i) = 2*sum_(r)(b_v(r)*Tr(Fmo*nsecond(p,r)))
  ! unres => *0.5
  val1 = mat_TrAB (soeo%mats%Fmo, sumn)
  if (.not. soeo%cfg_unres) then
    val1 = 2.0E0_realk * val1
  endif
  if (i <= Nact) then
    tmp(1) = val1
  else
    tmp(2) = val1
  endif

!write (soeo%lupri, *) 'First term =', val1, tmp

  !Second term: sigma_v(i) = sigma_v(i) - 2*Tr([Fmo,b_m]*nfirst(i))
  val1 = mat_TrAB (tmp_m, soeo%mats%nfirst(i))
  if (.not. soeo%cfg_unres) then
    val1 = 2.0E0_realk * val1
  endif
  if (i <= Nact) then
    tmp(1) = tmp(1) - val1
  else
    tmp(2) = tmp(2) - val1
  endif

!write (soeo%lupri, *) 'Second term =', val1, tmp

  !Third term: sigma_v(i) = sigma_v(i) + 4*Tr(GmoL*nfirst(i))
  val1 = mat_TrAB (GmoL, soeo%mats%nfirst(i))
  if (.not. soeo%cfg_unres) then
    val1 = 4.0E0_realk * val1
  endif
  if (i <= Nact) then
    tmp(1) = tmp(1) + val1
  else
    tmp(2) = tmp(2) + val1
  endif

!write (soeo%lupri, *) 'Third term =', val1, tmp

  !Element is put into sigma_v in the right place
  if (i <= Nact) then
    call mat_get_ab_elms (sigma_v, i, 1, tmp1)
    tmp1(1) = tmp(1)
    call mat_create_ab_elms (i, 1, tmp1, sigma_v)
  else
    call mat_get_ab_elms (sigma_v, i-Nact, 1, tmp1)
    tmp1(2) = tmp(2)
    call mat_create_ab_elms (i-Nact, 1, tmp1, sigma_v)
  endif
enddo
!-----------------------------------------------------------------------
!write (soeo%lupri, *) 'sigma_v ='
!call mat_print (sigma_v, 1, Nact, 1, 1, soeo%lupri)
!!stop 'velociraptor'

!Finalizations:
call mat_free (L)
call mat_free (GmoL)
call mat_free (commDb)
call mat_free (sumn)
call mat_free (tmp_m)

call mem_dealloc (tmp)
call mem_dealloc (tmp1)

end subroutine soeo_linear_transform
!=======================================================================



!> \brief Gets G(L) in the MO basis (see di_get_GAOL in dalton_interface.f90)
!> \author C. Nygaard
!> \date 2010
!> \param Lmo The matrix L in the MO basis
!> \param S The overlap matrix
!> \param C The MO orbital coefficients
!> \param Dao The density matrix in the AO basis
!> \param GmoL The wanted matrix G(L)
!> \param lupri logical unit number for standard output
!=======================================================================
subroutine soeo_get_GmoL (Lmo, S, C, Dao, GmoL, lupri, do_dft)

implicit none

!I/O
type(matrix), intent(in)    :: Lmo, S, C, Dao
integer, intent(in)         :: lupri
type(matrix), intent(inout) :: GmoL !Intent(out)
!Other
type(matrix)                :: Lao, GaoL, tmp
integer                     :: Nbast
logical, intent(in) :: do_dft
!do_dft = .false.

Nbast = Lmo%nrow

call mat_init (Lao , Nbast, Nbast)
call mat_init (GaoL, Nbast, Nbast)
call mat_init (tmp , Nbast, Nbast)

!write (lupri, *) 'C ='
!call mat_print (C, 1, Nbast, 1, Nbast, lupri)
!write (lupri, *) 'Lmo ='
!call mat_print (Lmo, 1, Nbast, 1, Nbast, lupri)

!Transform L to the AO-basis:
!-----------------------------------------------------------------------
call util_MO_to_AO_2 (S, C, Lmo, Lao, .false.)
!-----------------------------------------------------------------------

!write (lupri, *) 'Lao ='
!call mat_print (Lao, 1, Nbast, 1, Nbast, lupri)

!Get G(L) in the AO-basis:
!-----------------------------------------------------------------------
call di_get_GAOL_lsdalton(Lao, GaoL)
if (do_dft) then
  !dft only implemented for restricted
  call di_get_sigma_xc_cont_lsdalton(Lao,Dao,tmp)
  call mat_daxpy (0.5E0_realk, tmp, GaoL)
endif
!-----------------------------------------------------------------------

!write (lupri, *) 'GaoL ='
!call mat_print (GaoL, 1, Nbast, 1, Nbast, lupri)

!Transform G(L) back to the MO-basis:
!-----------------------------------------------------------------------
call util_AO_to_MO_2 (S, C, GaoL, GmoL, .true.)
!-----------------------------------------------------------------------

!write (lupri, *) 'GmoL ='
!call mat_print (GmoL, 1, Nbast, 1, Nbast, lupri)

call mat_free (Lao )
call mat_free (GaoL)
call mat_free (tmp )

end subroutine soeo_get_GmoL
!=======================================================================



!> \brief Finds the predicted energy change
!> \author C. Nygaard
!> \date 2010-07-01
!> \param soeo Contains all matrices used in soeo
!> \param K The orbital rotation matrix
!> \param deltatheta The change in occupation angles
!
!> Epred = (K,deltatheta)^T g + 1/2 * (K,deltatheta)^T H(K,deltatheta)
!
!=======================================================================
subroutine soeo_get_Epred (soeo, K, deltatheta)

implicit none

type(soeoItem), intent(inout) :: soeo
type(matrix), intent(in)      :: K, deltatheta

integer                       :: Nbast, Nact
type(matrix)                  :: sigma_m, sigma_v
type(matrix)                  :: grad_m, grad_v

Nbast = soeo%space%Nbast
Nact  = soeo%space%Nact

call mat_init (sigma_m, Nbast, Nbast)
call mat_init (sigma_v, Nact, 1)
call mat_init (grad_m, Nbast, Nbast)
call mat_init (grad_v, Nact, 1)

soeo%dEpred = 0.0E0_realk

call soeo_linear_transform (soeo, K, deltatheta, sigma_m, sigma_v) !H(v)

call soeo_get_gradient (soeo, grad_m, grad_v) !g

!write (soeo%lupri, *) 'v ='
!call mat_print (K, 1, Nbast, 1, Nbast, soeo%lupri)
!call mat_print (deltatheta, 1, Nact, 1, 1, soeo%lupri)
!write (soeo%lupri, *) 'g ='
!call mat_print (grad_m, 1, Nbast, 1, Nbast, soeo%lupri)
!call mat_print (grad_v, 1, Nact, 1, 1, soeo%lupri)
!write (soeo%lupri, *) 'v*g =', soeo_dotproduct (K, deltatheta, grad_m, grad_v)
!write (soeo%lupri, *) 'H(v) ='
!call mat_print (sigma_m, 1, Nbast, 1, Nbast, soeo%lupri)
!call mat_print (sigma_v, 1, Nact, 1, 1, soeo%lupri)
!write (soeo%lupri, *) 'v*H(v) =', soeo_dotproduct (K, deltatheta, sigma_m, sigma_v)

call soeo_daxpy (0.5E0_realk, sigma_m, sigma_v, grad_m, grad_v) !g + 1/2 H(v)

!write (soeo%lupri, *) 'g + 1/2 H(v) ='
!call mat_print (grad_m, 1, Nbast, 1, Nbast, soeo%lupri)
!call mat_print (grad_v, 1, Nact, 1, 1, soeo%lupri)

soeo%dEpred = soeo_dotproduct (K, deltatheta, grad_m, grad_v) !v*(g + 1/2 H(v))

!write (soeo%lupri, *) 'v*(g + 1/2 H(v)) =', soeo%dEpred

call mat_free (sigma_m)
call mat_free (sigma_v)
call mat_free (grad_m)
call mat_free (grad_v)

end subroutine soeo_get_Epred
!=======================================================================


!> \brief Calculates extra term on linear transform when Nelec is restricted
!> \author C. Nygaard
!> \date Feb. 18. 2011
!> \param space Information about the active space
!> \param theta The occupation angles
!> \param b The trial vector for change in occupation angles
!> \param sigmaplus The extra term to be calculated
!
!> Has only a vector part
!>   sigmaplus_i = -2*(cos^2(theta_i) - sin^2(theta_i))b_i
!=======================================================================
subroutine soeo_transformplus (space, unres, grandcan, theta, b, sigmaplus)

implicit none

type(soeoItem_space), intent(in) :: space
logical, intent(in)              :: unres, grandcan
type(matrix), intent(in)         :: theta, b
type(matrix), intent(inout)      :: sigmaplus !intent(out)

integer                          :: i
real(realk), pointer             :: thetaval(:), sigmaval(:), bi(:)

if (grandcan) then

  call mat_zero (sigmaplus)

else

  if (unres) then
    call mem_alloc (thetaval,2)
    call mem_alloc (sigmaval,2)
    call mem_alloc (bi,2)
  else
    call mem_alloc (thetaval,1)
    call mem_alloc (sigmaval,1)
    call mem_alloc (bi,1)
  endif

  call mat_zero (sigmaplus) 
  do i=1,space%Nact
    call mat_get_ab_elms (theta, space%Nocc+i, 1, thetaval)
    call mat_get_ab_elms (b, i, 1, bi)
    sigmaval(1) = cos(thetaval(1))*cos(thetaval(1))
    sigmaval(1) = sigmaval(1) - sin(thetaval(1))*sin(thetaval(1))
    sigmaval(1) = -2.0E0_realk*sigmaval(1)*bi(1)
    if (unres) then
      sigmaval(2) = cos(thetaval(2))*cos(thetaval(2))
      sigmaval(2) = sigmaval(2) - sin(thetaval(2))*sin(thetaval(2))
      sigmaval(2) = -2.0E0_realk*sigmaval(2)*bi(2)
    endif
    call mat_create_ab_elms (i, 1, sigmaval, sigmaplus)
  enddo

  call mem_dealloc (thetaval)
  call mem_dealloc (sigmaval)
  call mem_dealloc (bi)

!  if (.not. unres) then
!    call mat_scal (2.0E0_realk, sigmaplus)
!  endif

endif

end subroutine soeo_transformplus
!=======================================================================



!> \brief Calculates extra term on gradient when Nelec is restricted
!> \author C. Nygaard
!> \date Feb. 18. 2011
!> \param theta The occupation angles
!> \param b The trial vector for change in occupation angles
!> \param gradplus The extra term to be calculated
!
!> Has only a vector part
!>   gradplus_i = -2*cos(theta_i)sin(theta_i)
!=======================================================================
subroutine soeo_gradplus (space, unres, theta, gradplus)

implicit none

type(soeoItem_space), intent(in) :: space
logical, intent(in)              :: unres
type(matrix), intent(in)         :: theta
type(matrix), intent(inout)      :: gradplus !intent(out)

integer                          :: i
real(realk), pointer             :: thetaval(:), gradval(:)
real(realk)                      :: full, empty
real(realk), parameter           :: pi = 3.1415926553897932

if (unres) then
  call mem_alloc (thetaval,2)
  call mem_alloc (gradval,2)
else
  call mem_alloc (thetaval,1)
  call mem_alloc (gradval,1)
endif

!full = 1.0E-5_realk !if theta < full the orbital is fully occupied
!empty = pi/2.0E0_realk - 1.0E-5_realk !if theta > empty the orbital is empty

do i=1,space%Nact
  call mat_get_ab_elms (theta, space%Nocc+i, 1, thetaval)
!  if (thetaval(1) > full .and. thetaval(1) < empty) then
    gradval(1) = -2.0E0_realk*cos(thetaval(1))*sin(thetaval(1))
!  else
!    gradval(1) = 0.0E0_realk
!  endif
  if (unres) then
!    if (thetaval(2) > full .and. thetaval(2) < empty) then
      gradval(2) = -2.0E0_realk*cos(thetaval(2))*sin(thetaval(2))
!    else
!      gradval(2) = 0.0E0_realk
!    endif
  endif
  call mat_create_ab_elms (i, 1, gradval, gradplus)
enddo

!if (.not. unres) then
!  call mat_scal (2.0E0_realk, gradplus)
!endif

call mem_dealloc (thetaval)
call mem_dealloc (gradval)

end subroutine soeo_gradplus
!=======================================================================



!> \brief Determines the total number of electrons in the system
!> \author C. Nygaard
!> \date Mar. 11 2011
!> \param Nocc Number of fully occupied orbitals (location of active space)
!> \param t The old occupation angles
!> \param dt The correction to the occupation angles
!> \param unres True if calculation is unrestricted
!=======================================================================
function soeo_numberof_electrons (Nocc, t, dt, unres)

implicit none

real(realk) soeo_numberof_electrons
integer, intent(in)      :: Nocc
type(matrix), intent(in) :: t, dt
logical, intent(in)      :: unres
integer                  :: Nbast, Nact, i
real(realk), pointer     :: tmp1(:), tmp2(:)
real(realk)              :: N

Nbast = t%nrow
Nact = dt%nrow

if (unres) then
  call mem_alloc (tmp1,2)
  call mem_alloc (tmp2,2)
else
  call mem_alloc (tmp1,1)
  call mem_alloc (tmp2,1)
endif

N = 0.0E0_realk
do i=1,Nbast
  call mat_get_ab_elms (t, i, 1, tmp1)
  if (i>Nocc .and. i<=Nocc+Nact) then
    call mat_get_ab_elms (dt, i-Nocc, 1, tmp2)
    N = N + cos(tmp1(1)+tmp2(1))*cos(tmp1(1)+tmp2(1))
    if (unres) then
      N = N + cos(tmp1(2)+tmp2(2))*cos(tmp1(2)+tmp2(2))
    endif
  else
    N = N + cos(tmp1(1))*cos(tmp1(1))
    if (unres) then
      N = N + cos(tmp1(2))*cos(tmp1(2))
    endif
  endif
enddo

soeo_numberof_electrons = N

call mem_dealloc (tmp1)
call mem_dealloc (tmp2)

end function soeo_numberof_electrons
!=======================================================================

end module soeo_util
