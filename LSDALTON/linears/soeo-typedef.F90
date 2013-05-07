!> @file
!> Contains the soeo_typedef module
!
!> brief Structure with information for SOEO
!> author C. Nygaard
!> date 2010.07.07
module soeo_typedef

use matrix_operations

implicit none

!Contains:
!  type soeoItem_mats
!  type soeoItem_old
!  type soeoItem_space
!  type soeoItem_settings
!  type soeoItem
!  type soeoItem_input
!
!  soeoinp_set_default_config
!  soeoItem_init
!  soeoItem_init_space
!  soeoItem_init_settings
!  soeoItem_init_mats
!  soeoItem_init_old
!  soeoItem_free
!  soeoItem_free_mats
!  soeoItem_free_old



!=======================================================================
type soeoItem_mats

!AO matrices:
!> Density matrix in AO basis
type(matrix)              :: Dao
!> Fock matrix in the AO basis
type(matrix)              :: Fao
!> Overlap matrix
type(matrix)              :: S
!> MO coefficient matrix
type(matrix)              :: C

!MO matrices:
!> Density matrix in the MO basis
type(matrix)              :: Dmo
!> Fock matrix in the MO basis
type(matrix)              :: Fmo
!> First derivatives of Dmo wrt occupations
type(matrix), pointer :: nfirst(:)
!> Second derivatives of Dmo wrt occupations
type(matrix), pointer :: nsecond(:,:)
!> Occupations Dmo(i,i) = cos^2(theta(i))
type(matrix)              :: oldtheta

end type soeoItem_mats
!=======================================================================



!=======================================================================
type soeoItem_old

!> Occupations
type(matrix)              :: oldtheta
!> Density in MO basis
type(matrix)              :: Dmo
!> Density in AO basis
type(matrix)              :: Dao
!> Fock matrix in MO basis
type(matrix)              :: Fmo
!> Fock matrix in AO basis
type(matrix)              :: Fao
!> MO coefficient matrix
type(matrix)              :: C
!> First derivatives of Dmo wrt occupations
type(matrix), pointer :: nfirst(:)
!> Second derivatives of Dmo wrt occupations
type(matrix), pointer :: nsecond(:,:)

end type soeoItem_old
!=======================================================================



!=======================================================================
type soeoItem_space

!> Size of all space (number of basis functions)
integer                   :: Nbast
!> Size of occupied space
integer                   :: Nocc
!> Size of active space
integer                   :: Nact
!> Number of electrons
integer                   :: Nelec

end type soeoItem_space
!=======================================================================



!=======================================================================
type soeoItem_settings

!> Maximum macro-iterations
integer                   :: macromaxiter
!> Convergence threshold for macro-iterations
real(realk)               :: macrothresh
!> Maximum micro-iterations
integer                   :: micromaxiter
!> Convergence threshold for micro-iterations
real(realk)               :: microthresh
!> Trust radius
real(realk)               :: trust

end type soeoItem_settings
!=======================================================================



!=======================================================================
type soeoItem

!> Logical unit number for DALTON.OUT
integer                 :: lupri

!> Matrices used in soeo
type(soeoItem_mats)     :: mats

!> Old matrices (matrices from last iteration is saved in case of rejection of the step)
type(soeoItem_old)      :: old

!> Information about the active space
type(soeoItem_space)    :: space

!> Settings (thresholds, maximum number of iterations and trust radius)
type(soeoItem_settings) :: settings

!> Is calculation unrestricted?
logical                 :: cfg_unres
!> Grand canonican ensemble? (no restriction on total number of electrons)
logical                 :: cfg_grandcan
!> Is it a dft-calculation?
logical                 :: do_dft
!> Is this a test?
logical                 :: test
!> Print level
integer                 :: prnt

!Energy:
!> The total energy
real(realk)             :: Etotal
!> The predicted energy change
real(realk)             :: dEpred
!> The actual energy change
real(realk)             :: dE

!> The number of macroiterations done
integer                 :: iter

end type soeoItem
!=======================================================================



!=======================================================================
type soeoItem_input

!> Logical unit number for DALTON.OUT
integer     :: lupri

!> Should soeo be used?
logical     :: cfg_soeo
!> Is calculation unrestricted?
logical     :: cfg_unres
!> Grand canonican ensemble? (no restriction on total number of electrons)
logical     :: cfg_grandcan
!> Should matrices be saved after optimization?
logical     :: cfg_save
!> Should soeo-calculation be started from matrices saved in soeosave.out ?
logical     :: cfg_restart
!> Is it a dft_calculation?
logical     :: do_dft
!> Is this a test?
logical     :: test
!> Print level 
integer     :: prnt
!> Calculate dipole moment?
logical     :: cfg_dipole

!> Is active space given as input?
logical     :: spaceinput
!> Size of occupied space
integer     :: Nocc
!> Size of active space
integer     :: Nact

!> Are occupations given in input?
logical     :: occsinput
!> Number of fully occupied orbitals
integer     :: Nfullocc
!> Number of fractionally occupied orbitals
integer     :: Nfracocc
!> Fractional occupation numbers
real(realk) :: fracoccs(100) !dimension Nfracocc
!> Number of fully occupied beta orbitals
integer     :: Nfulloccb
!> Number of fractionally occupied beta orbitals
integer     :: Nfracoccb
!> Fractional occupation numbers in beta orbitals
real(realk) :: fracoccsb(100) !dimension Nfracoccb

!> Maximum macro-iterations
integer     :: macromaxiter
!> Convergence threshold for macro-iterations
real(realk) :: macrothresh
!> Maximum micro-iterations
integer     :: micromaxiter
!> Convergence threshold for micro-iterations
real(realk) :: microthresh
!> Trust radius
real(realk) :: trust

end type soeoItem_input
!=======================================================================



Contains



!> \brief Sets default options for soeo
!> \author C. Nygaard
!> \date Feb 2. 2011
!> \param soeoinp Contains options for soeo
!> Active space will be defined later, if it is not defined in the input
!=======================================================================
subroutine soeoinp_set_default_config (soeoinp)

implicit none

type(soeoItem_input) :: soeoinp

soeoinp%cfg_soeo     = .false.
soeoinp%cfg_unres    = .false.
soeoinp%cfg_grandcan = .false.
soeoinp%cfg_save     = .false.
soeoinp%cfg_restart  = .false.
soeoinp%do_dft       = .false.
soeoinp%test         = .false.
soeoinp%prnt         = 1
soeoinp%cfg_dipole   = .false.

soeoinp%spaceinput   = .false.
soeoinp%occsinput    = .false.
soeoinp%Nfracocc     = 0
soeoinp%Nfracoccb    = 0

soeoinp%macromaxiter = 200
soeoinp%macrothresh  = 1.0E-5_realk
soeoinp%micromaxiter = 200
soeoinp%microthresh  = 1.0E-8_realk
soeoinp%trust        = 0.5E0_realk

end subroutine soeoinp_set_default_config
!=======================================================================



!> \brief Initializes the structure soeo
!> \author C. Nygaard
!> \date Feb 2. 2011
!> \param soeo Structure containing all information needed in soeo
!> \param soeoinp Contains information given as input
!=======================================================================
subroutine soeoItem_init (Nbast, Nelec, soeo, soeoinp)

implicit none

integer, intent(in)              :: Nbast, Nelec
type(soeoItem), intent(inout)    :: soeo
type(soeoItem_input), intent(in) :: soeoinp

soeo%cfg_unres    = soeoinp%cfg_unres
soeo%cfg_grandcan = soeoinp%cfg_grandcan
soeo%do_dft       = soeoinp%do_dft
soeo%test         = soeoinp%test
soeo%lupri        = soeoinp%lupri
soeo%prnt         = soeoinp%prnt
soeo%iter         = 0

call soeoItem_init_space (Nbast, Nelec, soeo%space, soeoinp)
call soeoItem_init_settings (soeo%settings, soeoinp)
call soeoItem_init_mats (soeo%mats, soeo%space, soeo%cfg_unres)
call soeoItem_init_old (soeo%old, soeo%space, soeo%cfg_unres)

end subroutine soeoItem_init
!=======================================================================



!> \brief Defines the active space in soeo (default: all space is active)
!> \author C. Nygaard
!> \date Feb 2. 2011
!> \param space Information about the active space
!> \param input Information from input
!=======================================================================
subroutine soeoItem_init_space (Nbast, Nelec, space, input)

implicit none

integer, intent(in)                 :: Nbast, Nelec
type(soeoItem_space), intent(inout) :: space
type(soeoItem_input), intent(in) :: input

space%Nbast = Nbast
space%Nelec = Nelec
if (input%spaceinput) then !use input
  space%Nocc = input%Nocc
  space%Nact = input%Nact
  if (space%Nocc+space%Nact>space%Nbast) then
    !negative number of virtual orbitals...
    print *, 'Nocc =', space%Nocc
    print *, 'Nact =', space%Nact
    print *, 'Nbast =', space%Nbast
    call lsquit ('Nocc + Nact > Nbast in SOEO',-1)
  endif
else !default: all space is active
  space%Nocc = 0
  space%Nact = space%Nbast
endif

end subroutine soeoItem_init_space
!=======================================================================



!> \brief Initializes the settings for a soeo calculation
!> \author C. Nygaard
!> \date Feb 2. 2011
!> \param settings The settings to be initialized
!> \param input Settings from input
!=======================================================================
subroutine soeoItem_init_settings (settings, input)

implicit none

type(soeoItem_settings), intent(inout) :: settings
type(soeoItem_input), intent(in)       :: input

settings%macromaxiter = input%macromaxiter
settings%macrothresh  = input%macrothresh
settings%micromaxiter = input%micromaxiter
settings%microthresh  = input%microthresh
settings%trust        = input%trust

end subroutine soeoItem_init_settings
!=======================================================================



!> \brief Initializes the matrices used in the soeo iterations
!> \author C. Nygaard
!> \date Feb 2. 2011
!> \param mats Structure containing the matrices to be initialized
!> \param space Contains information about the active space
!> \param unres True if calculation is unrestricted
!=======================================================================
subroutine soeoItem_init_mats (mats, space, unres)

implicit none

type(soeoItem_mats), intent(inout) :: mats
type(soeoItem_space), intent(in)   :: space
logical, intent(in)                :: unres

integer                            :: ndim, i, j

call mat_init (mats%Dao     , space%Nbast, space%Nbast)
call mat_init (mats%Fao     , space%Nbast, space%Nbast)
call mat_init (mats%S       , space%Nbast, space%Nbast)
call mat_init (mats%C       , space%Nbast, space%Nbast)
call mat_init (mats%Dmo     , space%Nbast, space%Nbast)
call mat_init (mats%Fmo     , space%Nbast, space%Nbast)
call mat_init (mats%oldtheta, space%Nbast, 1          )
if (unres) then
  ndim = 2.0E0_realk*space%Nact
else
  ndim = space%Nact
endif
allocate (mats%nfirst(ndim), mats%nsecond(ndim,ndim))
!call mem_alloc (mats%nfirst,ndim)
!call mem_alloc (mats%nsecond,ndim,ndim)
do i=1,ndim
  call mat_init(mats%nfirst(i), space%Nbast, space%Nbast)
  do j=1,ndim
    call mat_init (mats%nsecond(i,j), space%Nbast, space%Nbast)
  enddo
enddo

end subroutine soeoItem_init_mats
!=======================================================================



!> \brief Initializes the old matrices used in the soeo iterations (for stepbacks)
!> \author C. Nygaard
!> \date Feb 2. 2011
!> \param old Structure containing the matrices to be initialized
!> \param space Contains information about the active space
!> \param unres True if calculation is unrestricted
!=======================================================================
subroutine soeoItem_init_old (old, space, unres)

implicit none

type(soeoItem_old), intent(inout)  :: old
type(soeoItem_space), intent(in)   :: space
logical, intent(in)                :: unres

integer                            :: ndim, i, j

call mat_init (old%Dao     , space%Nbast, space%Nbast)
call mat_init (old%Fao     , space%Nbast, space%Nbast)
call mat_init (old%C       , space%Nbast, space%Nbast)
call mat_init (old%Dmo     , space%Nbast, space%Nbast)
call mat_init (old%Fmo     , space%Nbast, space%Nbast)
call mat_init (old%oldtheta, space%Nbast, 1          )
if (unres) then
  ndim = 2.0E0_realk*space%Nact
else
  ndim = space%Nact
endif
allocate (old%nfirst(ndim), old%nsecond(ndim,ndim))
!call mem_alloc (old%nfirst,ndim)
!call mem_alloc (old%nsecond,ndim,ndim)
do i=1,ndim
  call mat_init(old%nfirst(i), space%Nbast, space%Nbast)
  do j=1,ndim
    call mat_init (old%nsecond(i,j), space%Nbast, space%Nbast)
  enddo
enddo

end subroutine soeoItem_init_old
!=======================================================================



!> \brief Frees the matrices initialized in soeoItem_init
!> \author C. Nygaard
!> \date Feb 2. 2011
!> \param soeo Structure containing all the matrices
!=======================================================================
subroutine soeoItem_free (soeo)

implicit none

type(soeoItem), intent(inout) :: soeo

call soeoItem_free_mats (soeo%mats, soeo%space, soeo%cfg_unres)
call soeoItem_free_old (soeo%old, soeo%space, soeo%cfg_unres)

end subroutine soeoItem_free
!=======================================================================



!> \brief Frees the matrices initialized in soeoItem_init_mats
!> \author C. Nygaard
!> \date Feb 2. 2011
!> \param mats Structure containing the matrices
!> \param space Information about the dimensions of the matrices
!> \param unres True if calculation is unrestricted
!=======================================================================
subroutine soeoItem_free_mats (mats, space, unres)

implicit none

type(soeoItem_mats), intent(inout) :: mats
type(soeoItem_space), intent(in)   :: space
logical, intent(in)                :: unres

integer                            :: ndim, i, j

call mat_free (mats%Dao     )
call mat_free (mats%Fao     )
call mat_free (mats%S       )
call mat_free (mats%C       )
call mat_free (mats%Dmo     )
call mat_free (mats%Fmo     )
call mat_free (mats%oldtheta)
if (unres) then
  ndim = 2*space%Nact
else
  ndim = space%Nact
endif
do i=1,ndim
  call mat_free(mats%nfirst(i))
  do j=1,ndim
    call mat_free (mats%nsecond(i,j))
  enddo
enddo
deallocate (mats%nfirst, mats%nsecond)
!call mem_dealloc (mats%nfirst)
!call mem_dealloc (mats%nsecond)

end subroutine soeoItem_free_mats
!=======================================================================



!> \brief Frees matrices initialized in soeoItem_init_old
!> \author C. Nygaard
!> \date Feb 2. 2011
!> \param old Structure containing the matrices
!> \param space Contains information about the active space
!> \param unres True if calculation is unrestricted
!=======================================================================
subroutine soeoItem_free_old (old, space, unres)

implicit none

type(soeoItem_old), intent(inout)  :: old
type(soeoItem_space), intent(in)   :: space
logical, intent(in)                :: unres

integer                            :: ndim, i, j

call mat_free (old%Dao     )
call mat_free (old%Fao     )
call mat_free (old%C       )
call mat_free (old%Dmo     )
call mat_free (old%Fmo     )
call mat_free (old%oldtheta)
if (unres) then
  ndim = 2.0E0_realk*space%Nact
else
  ndim = space%Nact
endif
do i=1,ndim
  call mat_free(old%nfirst(i))
  do j=1,ndim
    call mat_free (old%nsecond(i,j))
  enddo
enddo
deallocate (old%nfirst, old%nsecond)
!call mem_dealloc (old%nfirst)
!call mem_dealloc (old%nsecond)

end subroutine soeoItem_free_old
!=======================================================================

end module soeo_typedef
