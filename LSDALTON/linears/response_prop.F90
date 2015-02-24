!> @fileddd
!> Contains property wrappers to the response_driver_module in response_driver.f90

!> \author Thomas Kjaergaard and Kasper Kristensen
!> \date 2010-03
module response_wrapper_module
#ifdef VAR_RSP
  use precision
  use memory_handling, only: mem_alloc, mem_dealloc
  use files, only: lsopen, lsclose
  use response_driver_module, only: RspFunc, get_first_order_property, &
       & initialize_rspfunc, linear_response, quadratic_response_2nplus1,&
       & cubic_response, quadratic_response_nplus1, &
       & transition_moment_density_matrix, print_response_func,&
       & rspfunc_or_transmoment
  use matrix_module, only: matrix
  use lsdalton_rsp_contribs
  use lstiming
  !warning matrix_operations and matrix_defop do not mix
  use matrix_operations,only: mat_write_to_disk
  use lsdalton_matrix_defop
  use rspsolver, only: rsp_molcfg
  use lsdalton_rsp_equations, only: rsp_eq_sol, pert_dens, pert_fock, rsp_eq_sol_empty
  use response_wrapper_type_module, only: ALPHAinputItem, &
       & BETAinputItem, GAMMAinputItem, TPAinputItem, ESGinputItem, &
       & DTPAinputItem, RSPSOLVERinputitem, MCDinputItem, ESDinputItem
  use decompMod, only: decompItem
  use TYPEDEFTYPE, only: LSSETTING
  use complexsolver, only: rsp_complex_init, rsp_complex_solver
  use complexsymsolver, only: rsp_sym_complex_init, rsp_sym_complex_solver


  public MCDresponse_driver, ALPHAresponse_driver, BETAresponse_driver, &
       & GAMMAresponse_driver, OPAresponse_driver, &
       & Calculate_and_store_transition_density_matrices, &
       & free_transition_density_matrices, TPAresponse_driver, &
       & ESGresponse_driver, Get_dipole_moment, ESDresponse_driver, &
       & DTPAresponse_driver, NMRshieldresponse_driver, &
       & get_excitation_energies, dipolemomentmatrix_driver
  private

Contains



  !> \brief Calculate dipole moment.
  !> \author Kasper Kristensen                                                                
  !> \date 2010-08
  subroutine Get_dipole_moment(molcfg,F,D,S,doPrint,DipoleMoment)
    implicit none
    !> Info on molecule needed by solver and integral programs
    type(rsp_molcfg),intent(inout)    :: molcfg
    !> Unperturbed Fock matrix
    type(matrix), intent(in)          :: F
    !> Unperturbed density matrix
    type(matrix), intent(in)          :: D
    !> Unperturbed overlap matrix
    type(matrix), intent(in)          :: S
    !> Print out dipole or not?
    logical, intent(in) :: doPrint
    !> Molecular dipole moment vector
    real(realk), intent(inout) :: DipoleMoment(3)
    complex(realk) :: DipoleMoment_complex(3)
    integer :: i
    real(realk) :: au_to_debye, au_to_SI, normDipole
    integer, pointer :: lupri

    lupri => molcfg%lupri

    ! Zero dipole
    DipoleMoment(:) = 0E0_realk
    DipoleMoment_complex(:) = (0E0_realk,0E0_realk)

    ! Get dipole, we need it to be complex here to use general framework.
    call Get_first_order_property(molcfg,F,D,S,DipoleMoment_complex,3,'EL  ')

    ! By convention the dipole moment equals MINUS the energy derivative.
    ! Get_first_order_property calls rsp_oneave, where the sign convention
    ! is always based on energy derivatives.
    normDipole=0E0_realk
    do i=1,3
       DipoleMoment(i) = -real(DipoleMoment_complex(i))
       normDipole = normDipole + DipoleMoment(i)**2
    end do
    normDipole = sqrt(normDipole)

    if(doPrint) then
       call PrintDipoleMoment(DipoleMoment,lupri)
    end if

  end subroutine Get_dipole_moment

  subroutine PrintDipoleMoment(DipoleMoment,lupri)
    implicit none
    real(realk), intent(in) :: DipoleMoment(3)
    integer,intent(in) :: lupri
    integer :: i
    real(realk) :: au_to_debye, au_to_SI, normDipole
    normDipole=0E0_realk
    do i=1,3
       normDipole = normDipole + DipoleMoment(i)**2
    end do
    normDipole = sqrt(normDipole)         
    au_to_debye=2.54175
    au_to_SI=8.47835
    write(lupri,*)
    write(lupri,*)
    write(lupri,*)
    write(lupri,'(6X,A)') '                 Permanent dipole moment'
    write(lupri,'(6X,A)') '                 -----------------------'
    write(lupri,'(12X,A)') '   au              Debye           10**-30 C m'
    write(lupri,'(5X,3g18.6)') normDipole, au_to_debye*normDipole, au_to_SI*normDipole
    write(lupri,*)
    write(lupri,*)
    write(lupri,'(6X,A)') '                 Dipole moment components'
    write(lupri,'(6X,A)') '                 ------------------------'
    write(lupri,'(12X,A)') '   au              Debye           10**-30 C m'
    write(lupri,'(1X,A,3X,3g18.6)') 'x', DipoleMoment(1), &
         & au_to_debye*DipoleMoment(1), au_to_SI*DipoleMoment(1)
    write(lupri,'(1X,A,3X,3g18.6)') 'y', DipoleMoment(2), &
         & au_to_debye*DipoleMoment(2), au_to_SI*DipoleMoment(2)
    write(lupri,'(1X,A,3X,3g18.6)') 'z', DipoleMoment(3), &
         & au_to_debye*DipoleMoment(3), au_to_SI*DipoleMoment(3)
    
    write(lupri,*)
    write(lupri,*)
    write(lupri,*)
  end subroutine PrintDipoleMoment

  !> \brief Driver for calculating polarizabilites for frequencies defined in input.
  !> If no frequencies are specified the frequency is set to zero.
  !> \author Kasper Kristensen                                                                
  !> \date 2010-08
  subroutine ALPHAresponse_driver(molcfg,F,D,S,alphainput)
    implicit none
    !> Info on molecule needed by solver and integral programs
    type(rsp_molcfg),intent(inout)    :: molcfg
    !> Unperturbed Fock matrix
    type(matrix), intent(in)          :: F
    !> Unperturbed density matrix
    type(matrix), intent(in)          :: D
    !> Unperturbed overlap matrix
    type(matrix), intent(in)          :: S
    !> Contains alpha input
    type(ALPHAinputItem),intent(inout)       :: ALPHAinput
    !> Dimension of response vector
    integer :: response_size
    !> Vector for collecting response results          
    complex(realk), allocatable :: rsp_results_vec(:,:)
    type(RspFunc), allocatable :: MyRspFunc(:)
    complex(realk), allocatable :: alpha(:,:,:), isotropic(:)
    integer :: i,j, k, counter, nfreq
    integer, pointer :: lupri

    lupri => molcfg%lupri

    ! 1. Determine number of different frequencies and set frequencies in response function type
    ! ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

    ! Check that input is consistent
    ! (If both real and imaginary frequencies have been specified
    !  the numbers of real and imaginary frequencies must be equal.)
    if(alphainput%real_frequencies_in_input &
         & .and. alphainput%imag_frequencies_in_input) then
       if(alphainput%nfreq /= alphainput%nimfreq) then
          write(lupri,*) 'Error in ALPHAresponse_driver: The number of real and imaginary frequencies &
               & specified in the input file must be equal!'
          CALL lsQUIT('The number of real and imaginary frequencies &
               & specified in the input file must be equal!',lupri)
       end if
    end if

    nfreq = alphainput%nfreq
    allocate(MyRspFunc(nfreq))
    allocate(isotropic(nfreq))

    ! Initialize response functions
    do i=1,nfreq
       call initialize_RspFunc(MyRspFunc(i))
    end do

    ! If no frequencies have been defined, only calculate static polarizability (freq=0).
    ! Note: The frequencies are complex and we set the frequency omegaB for operator B based on input,
    ! whereas for operator A, omegaA=-omegaB.
    SetFrequencies: do i=1,nfreq

       ! Real frequencies
       if(alphainput%real_frequencies_in_input) then
          MyRspFunc(i)%freq(2)=(1E0_realk,0E0_realk)*alphainput%bfreq(i)
       else
          MyRspFunc(i)%freq(2)=(0E0_realk,0E0_realk)
       end if

       ! Imaginary frequencies
       if(alphainput%imag_frequencies_in_input) then
          MyRspFunc(i)%freq(2)=MyRspFunc(i)%freq(2) + (0E0_realk,1E0_realk)*alphainput%imbfreq(i)
       end if

       MyRspFunc(i)%freq(1) = - MyRspFunc(i)%freq(2) ! omegaA = -omegaB

    end do SetFrequencies


    ! 2. Set parameters other than frequencies
    ! ''''''''''''''''''''''''''''''''''''''''

    ! Allocate response result vector (3x3=9) and alpha tensor (3,3) for each frequency.
    response_size = 3*3
    allocate(rsp_results_vec(nfreq,response_size))
    allocate(alpha(nfreq,3,3))

    SetOtherParameters: do i=1,nfreq
       ! Order of response function=2 (linear response)
       MyRspFunc(i)%order = 2

       ! Operators = (EL,EL)
       MyRspFunc(i)%code(1) = 'EL  '
       MyRspFunc(i)%code(2) = 'EL  '

       ! Calculate all 9 tensor components (xx,xy,..., zz)
       MyRspFunc(i)%first_comp(1) = 1
       MyRspFunc(i)%first_comp(2) = 1
       MyRspFunc(i)%dims(1) = 3
       MyRspFunc(i)%dims(2) = 3


       ! Zero response function and polarizabilities
       rsp_results_vec(i,:)=(0E0_realk,0E0_realk)
       alpha(i,:,:)=(0E0_realk,0E0_realk)

    end do SetOtherParameters



    ! 3. Calculate linear response functions
    ! ''''''''''''''''''''''''''''''''''''''

    ! Linear response function <<EL;EL>> for each frequency.
    CalculateLinearResponseFunc: do i=1,nfreq

       ! Calculate <<EL;EL>> response functions for all 9 tensor components
       call linear_response(molcfg,F,D,S,MyRspFunc(i), &
            & rsp_results_vec(i,1:response_size),response_size)

       ! Alpha equals MINUS the linear response function and the linear response
       ! function results are arranged in the vector rsp_results_vec as follows:
       ! (1,1) (1,2) (1,3) (2,1) (2,2) (2,3) (3,1) (3,2) (3,3).
       ! Thus, for each frequency (i) alpha may be calculated in the following way:
       counter=0
       do j=1,3
          do k=1,3
             counter=counter+1
             alpha(i,j,k) = -rsp_results_vec(i,counter)
          end do
       end do

       ! Isotropically averaged polarizability
       Isotropic(i) = ( alpha(i,1,1) + alpha(i,2,2) + alpha(i,3,3) )/3E0_realk

    end do CalculateLinearResponseFunc


    ! 4. Print alpha tensor results
    ! '''''''''''''''''''''''''''''

    write(lupri,*) 
    write(lupri,*) 
    write(lupri,*) 
    write(lupri,'(1X,A)') '*************************************************************'
    write(lupri,'(1X,A)') '*          POLARIZABILITY TENSOR RESULTS (in a.u.)          *'
    write(lupri,'(1X,A)') '*************************************************************'
    write(lupri,*) 
    write(lupri,*)
    write(lupri,*) 
    write(lupri,*) 
    PrintFrequencyLoop: do i=1,nfreq

       if(alphainput%imag_frequencies_in_input) then ! print complex frequency
          write(lupri,'(1X,A,g12.5,A,g12.5,A)')  'Frequency = ', &
               & real(MyRspFunc(i)%freq(2)), " + ", imag(MyRspFunc(i)%freq(2)), "i"
          write(lupri,*) '===================================='
       else ! print real frequency
          write(lupri,'(1X,A,g12.5)')  'Frequency = ', &
               & real(MyRspFunc(i)%freq(2))
          write(lupri,*) '=================='
       end if

       write(lupri,*) 
       if(alphainput%imag_frequencies_in_input) then ! complex frequencies
          write(lupri,*) 
          write(lupri,*) 'Real part of polarizability tensor'
          write(lupri,*) "''''''''''''''''''''''''''''''''''"
       end if

       write(lupri,*) "                Ex               Ey               Ez"
       write(lupri,'(1X,A, 3g18.8)') "Ex    ", real(alpha(i,1,:))
       write(lupri,'(1X,A, 3g18.8)') "Ey    ", real(alpha(i,2,:))
       write(lupri,'(1X,A, 3g18.8)') "Ez    ", real(alpha(i,3,:))
       write(lupri,*) 
       if(alphainput%imag_frequencies_in_input) then ! complex frequencies
          write(lupri,'(1X,A,g18.8)') 'Isotropic real polarizability = ', real(isotropic(i))
       else
          write(lupri,'(1X,A,g18.8)') 'Isotropic polarizability = ', real(isotropic(i))
       end if


       if(alphainput%imag_frequencies_in_input) then ! complex frequencies
          write(lupri,*) 
          write(lupri,*) 
          write(lupri,*) 'Imaginary part of polarizability tensor'
          write(lupri,*) "'''''''''''''''''''''''''''''''''''''''"
          write(lupri,*) "                Ex               Ey               Ez"
          write(lupri,'(1X,A, 3g18.8)') "Ex    ", imag(alpha(i,1,:))
          write(lupri,'(1X,A, 3g18.8)') "Ey    ", imag(alpha(i,2,:))
          write(lupri,'(1X,A, 3g18.8)') "Ez    ", imag(alpha(i,3,:))
          write(lupri,*) 
          write(lupri,'(1X,A,g18.8)') 'Isotropic imaginary polarizability = ', imag(isotropic(i))
       end if

       write(lupri,*) 
       write(lupri,*) 
       write(lupri,*) 
       write(lupri,*) 
       write(lupri,*) 

    end do PrintFrequencyLoop


    ! Deallocate stuff
    deallocate(rsp_results_vec)
    deallocate(isotropic)
    deallocate(MyRspFunc)
    deallocate(alpha)

    if(alphainput%real_frequencies_in_input) deallocate(alphainput%BFREQ)
    if(alphainput%imag_frequencies_in_input) deallocate(alphainput%IMBFREQ)


  end subroutine ALPHAresponse_driver





  !> \brief Driver for calculating 1st hyperpolarizabilites for frequencies defined in input.
  !>  If a given frequency is not specified in the input, it is set to zero.
  !> \author Kasper Kristensen                                                                
  !> \date 2010-08
  subroutine BETAresponse_driver(molcfg,F,D,S,betainput,DipoleMoment)
    implicit none
    !> Info on molecule needed by solver and integral programs
    type(rsp_molcfg),intent(inout)    :: molcfg
    !> Unperturbed Fock matrix
    type(matrix), intent(in)          :: F
    !> Unperturbed density matrix
    type(matrix), intent(in)          :: D
    !> Unperturbed overlap matrix
    type(matrix), intent(in)          :: S
    !> Contains beta input
    type(BETAinputItem),intent(inout)       :: BETAinput
    !> Dipole moment
    real(realk) :: DipoleMoment(3)
    !> Dimension of response vector
    integer :: response_size
    !> Vector for collecting response results          
    complex(realk), allocatable :: rsp_results_vec(:,:)
    type(RspFunc), allocatable :: MyRspFunc(:)
    complex(realk), allocatable :: beta(:,:,:,:)
    integer :: i,j, k, l, counter, nfreq
    character(len=1) :: comp(3)
    real(realk) :: DipoleNorm, DipoleUnit(3)
    complex(realk), allocatable :: BetaPara(:), BetaPerp(:)
    logical :: DipoleIsZero, ImagFreq
    integer, pointer :: lupri

    lupri => molcfg%lupri

    ! 1. Determine number of different frequencies and set frequencies in response function type
    ! ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

    ! Check that input is consistent.
    ! Code quite ugly, but it probably pays off to do this check...
    ! If present in input - number of real and imaginary frequencies must be identical.


    ! BFREQ and IMBFREQ
    if(betainput%real_bfrequencies_in_input &
         & .and. betainput%imag_bfrequencies_in_input) then
       if(betainput%nbfreq /= betainput%nimbfreq) then
          write(lupri,*) 'Error in BETAresponse_driver: The number of real and imaginary frequencies &
               & specified for operator B in the input file must be equal!'
          CALL lsQUIT('The number of real and imaginary frequencies  for operator B&
               & specified in the input file must be equal!',lupri)
       end if
    end if

    ! CFREQ and IMCFREQ
    if(betainput%real_cfrequencies_in_input &
         & .and. betainput%imag_cfrequencies_in_input) then
       if(betainput%ncfreq /= betainput%nimcfreq) then
          write(lupri,*) 'Error in BETAresponse_driver: The number of real and imaginary frequencies &
               & specified for operator C in the input file must be equal!'
          CALL lsQUIT('The number of real and imaginary frequencies  for operator C&
               & specified in the input file must be equal!',lupri)
       end if
    end if

    ! Also - if specified - the number of frequencies for operator B
    ! and operator C must be the same
    if(betainput%real_bfrequencies_in_input &
         & .and. betainput%real_cfrequencies_in_input) then
       if(betainput%nbfreq /= betainput%ncfreq) then
          write(lupri,*) 'Error in BETAresponse_driver: The number of frequencies specfied &
               &for operator B and C must be equal!'
          CALL lsQUIT('Error in BETAresponse_driver: The number of frequencies specfied &
               &for operator B and C must be equal!',lupri)
       end if
    end if

    ! Set number of frequencies equal to the number of frequencies for operator B or C
    ! (in case only one of them has been specified by the input)
    nfreq = max(betainput%nbfreq,betainput%ncfreq)
    allocate(MyRspFunc(nfreq))
    allocate(BetaPara(nfreq))
    allocate(BetaPerp(nfreq))

    ! Initialize response functions
    do i=1,nfreq
       call initialize_RspFunc(MyRspFunc(i))
    end do

    ! If no frequencies have been defined, only calculate static polarizability 
    ! (omegaA = omegaB = omegaC = 0).
    SetFrequencies: do i=1,nfreq

       ! Real frequencies operator B
       if(betainput%real_bfrequencies_in_input) then
          MyRspFunc(i)%freq(2)=(1E0_realk,0E0_realk)*betainput%bfreq(i)
       else
          MyRspFunc(i)%freq(2)=(0E0_realk,0E0_realk)
       end if

       ! Imaginary frequencies operator B
       if(betainput%imag_bfrequencies_in_input) then
          MyRspFunc(i)%freq(2)=MyRspFunc(i)%freq(2) + (0E0_realk,1E0_realk)*betainput%imbfreq(i)
       end if


       ! Real frequencies operator C
       if(betainput%real_cfrequencies_in_input) then
          MyRspFunc(i)%freq(3)=(1E0_realk,0E0_realk)*betainput%cfreq(i)
       else
          MyRspFunc(i)%freq(3)=(0E0_realk,0E0_realk)
       end if

       ! Imaginary frequencies operator C
       if(betainput%imag_cfrequencies_in_input) then
          MyRspFunc(i)%freq(3)=MyRspFunc(i)%freq(3) + (0E0_realk,1E0_realk)*betainput%imcfreq(i)
       end if

       ! omegaA = -omegaB-omegaC
       MyRspFunc(i)%freq(1) = - MyRspFunc(i)%freq(2) - MyRspFunc(i)%freq(3) 

    end do SetFrequencies


    ! Any imaginary frequencies?
    if(betainput%imag_bfrequencies_in_input .or. betainput%imag_cfrequencies_in_input) then
       ImagFreq=.true.
    else
       ImagFreq=.false. 
    end if


    ! 2. Set parameters other than frequencies
    ! ''''''''''''''''''''''''''''''''''''''''

    ! Allocate response result vector (3x3x3=27) and beta tensor (3,3,3) for each frequency.
    response_size = 3*3*3
    allocate(rsp_results_vec(nfreq,response_size))
    allocate(beta(nfreq,3,3,3))

    SetOtherParameters: do i=1,nfreq
       ! Order of response function=3 (quadratic response)
       MyRspFunc(i)%order = 3

       ! Operators = (EL,EL,EL)
       MyRspFunc(i)%code(1) = 'EL  '
       MyRspFunc(i)%code(2) = 'EL  '
       MyRspFunc(i)%code(3) = 'EL  '

       ! Calculate all 27 tensor components (xxx,xxy,..., zzz) by default
       MyRspFunc(i)%first_comp(1) = 1
       MyRspFunc(i)%first_comp(2) = 1
       MyRspFunc(i)%first_comp(3) = 1
       MyRspFunc(i)%dims(1) = 3
       MyRspFunc(i)%dims(2) = 3
       MyRspFunc(i)%dims(3) = 3


       ! Zero response function and polarizabilities
       rsp_results_vec(i,:)=(0E0_realk,0E0_realk)
       beta(i,:,:,:)=(0E0_realk,0E0_realk)

    end do SetOtherParameters


    ! 3. Calculate quadratic response functions
    ! '''''''''''''''''''''''''''''''''''''''''

    ! Quadratic response function <<EL;EL,EL>> for each frequency.
    CalculateQuadResponseFunc: do i=1,nfreq

       ! Calculate <<EL;EL,EL>> response functions for all 27 tensor components
       ! using 2n+1 rule
       call quadratic_response_2nplus1(molcfg,F,D,S,MyRspFunc(i), &
            & rsp_results_vec(i,1:response_size),response_size)

       ! Beta equals MINUS the quadratic response function, and the quadratic response
       ! function results are arranged in the vector rsp_results_vec as follows:
       !> (1,1,1) (1,1,2) (1,1,3) (1,2,1) (1,2,2) (1,2,3) (1,3,1) (1,3,2) (1,3,3)
       !> (2,1,1) (2,1,2) (2,1,3) (2,2,1) (2,2,2) (2,2,3) (2,3,1) (2,3,2) (2,3,3)
       !> (3,1,1) (3,1,2) (3,1,3) (3,2,1) (3,2,2) (3,2,3) (3,3,1) (3,3,2) (3,3,3)

       ! Thus, for each frequency i beta may be calculated in the following way:
       counter=0
       do j=1,3
          do k=1,3
             do l=1,3
                counter=counter+1
                beta(i,j,k,l) = -rsp_results_vec(i,counter)
             end do
          end do
       end do

    end do CalculateQuadResponseFunc



    ! 4. Calculate unit vector in dipole moment direction
    ! '''''''''''''''''''''''''''''''''''''''''''''''''''

    DipoleNorm = DipoleMoment(1)**2+DipoleMoment(2)**2+DipoleMoment(3)**2
    DipoleNorm = sqrt(DipoleNorm)

    ! Check if dipole is zero
    if(DipoleNorm < 1e-9_realk) then
       DipoleIsZero = .true.
    else
       DipoleIsZero = .false.
    end if

    ! Dipole moment unit vector
    if(.not. DipoleIsZero) then
       do j=1,3
          DipoleUnit(j) = DipoleMoment(j)/DipoleNorm
       end do
    end if


    ! 5. Get isotropically averaged beta parallel and beta perpendicular
    ! ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

    ! For each frequency i, get averaged betas (see subroutine Get_averaged_betas for definitions).
    ! Only do this if dipole is nonzero!
    if(.not. DipoleIsZero) then
       do i=1,nfreq
          call Get_averaged_betas(BetaPara(i),BetaPerp(i),beta(i,1:3,1:3,1:3),DipoleUnit)
       end do
    end if


    ! 6. Print beta tensor results
    ! ''''''''''''''''''''''''''''

    ! Character component vector
    comp(1) = 'X'
    comp(2) = 'Y'
    comp(3) = 'Z'

    write(lupri,*) 
    write(lupri,*) 
    write(lupri,*) 
    write(lupri,'(1X,A)') '********************************************************************'
    write(lupri,'(1X,A)') '*        FIRST HYPERPOLARIZABILITY TENSOR RESULTS (in a.u.)        *'
    write(lupri,'(1X,A)') '********************************************************************'
    write(lupri,*) 
    write(lupri,*)

    ! Only print isotropically averaged betas if dipole is nonzero
    if(.not. DipoleIsZero) then
       write(lupri,*) "The isotropically averaged parallel and perpendicular beta values"
       write(lupri,*)  "are calculated in the following way:"
       write(lupri,*) 
       write(lupri,*) "BetaParallel      = 1/5 * SUM_{i,j} &
            &( beta_{jii} + beta_{iji} + beta_{iij} ) * DipoleUnit(j)"
       write(lupri,*) "BetaPerpendicular = 1/5 * SUM_{i,j} &
            &( 2beta_{jii} - 3beta_{iji} + 2beta_{iij} ) * DipoleUnit(j)"
       write(lupri,*) 
       write(lupri,*) "where DipoleUnit is a unit vector in the direction of the permanent dipole moment,"
       write(lupri,*) "and the i and j summations run independently over x,y,z."
       write(lupri,*) 
       write(lupri,*) 
    end if
    write(lupri,*) 
    write(lupri,*) 
    PrintFrequencyLoop: do i=1,nfreq

       if(ImagFreq) then  
          write(lupri,'(1X,A,g12.5,A,g12.5,A)')  'Frequency B = ', &
               & real(MyRspFunc(i)%freq(2)), " + ", imag(MyRspFunc(i)%freq(2)), "i"
          write(lupri,'(1X,A,g12.5,A,g12.5,A)')  'Frequency C = ', &
               & real(MyRspFunc(i)%freq(3)), " + ", imag(MyRspFunc(i)%freq(3)), "i"
          write(lupri,*) '======================================'
       else
          write(lupri,'(1X,A,g12.5)')  'Frequency B = ', &
               & real(MyRspFunc(i)%freq(2))
          write(lupri,'(1X,A,g12.5)')  'Frequency C = ', &
               & real(MyRspFunc(i)%freq(3))
          write(lupri,*) '========================'
       end if

       write(lupri,*)

       if(ImagFreq) then
          write(lupri,*) "Components        Real part"
          write(lupri,*) "'''''''''''''''''''''''''''"
       else
          write(lupri,*) "Components        1st hyp. "
          write(lupri,*) "''''''''''''''''''''''''''"
       end if

       j_loopREAL: do j=1,3
          k_loopREAL: do k=1,3
             l_loopREAL: do l=1,3
                write(lupri,'(4X,3A,6X,g18.8)') comp(j),comp(k),comp(l),real(beta(i,j,k,l))
             end do l_loopREAL
          end do k_loopREAL
       end do j_loopREAL

       ! Only print isotropically averaged betas if dipole is nonzero
       if(.not. DipoleIsZero) then
          write(lupri,*) 
          write(lupri,'(a,g18.8)') 'BetaParallel      =', real(BetaPara(i))
          write(lupri,'(a,g18.8)') 'BetaPerpendicular =', real(BetaPerp(i))
       end if

       if(ImagFreq) then
          write(lupri,*) 
          write(lupri,*) 
          write(lupri,*) 'Components         Imag part'
          write(lupri,*) "'''''''''''''''''''''''''''''"

          j_loopIMAG: do j=1,3
             k_loopIMAG: do k=1,3
                l_loopIMAG: do l=1,3
                   write(lupri,'(4X,3A,6X,g18.8)') comp(j),comp(k),comp(l),imag(beta(i,j,k,l))
                end do l_loopIMAG
             end do k_loopIMAG
          end do j_loopIMAG

          ! Only print isotropically averaged betas if dipole is nonzero
          if(.not. DipoleIsZero) then
             write(lupri,*) 
             write(lupri,'(a,g18.8)') 'BetaParallel      =', imag(BetaPara(i))
             write(lupri,'(a,g18.8)') 'BetaPerpendicular =', imag(BetaPerp(i))
          end if
       end if



       write(lupri,*) 
       write(lupri,*) 
       write(lupri,*) 
       write(lupri,*) 
       write(lupri,*) 

    end do PrintFrequencyLoop


    deallocate(MyRspFunc)
    deallocate(BetaPara)
    deallocate(BetaPerp)
    deallocate(rsp_results_vec)
    deallocate(beta)

    if(betainput%real_bfrequencies_in_input) deallocate(betainput%BFREQ)
    if(betainput%imag_bfrequencies_in_input) deallocate(betainput%IMBFREQ)
    if(betainput%real_cfrequencies_in_input) deallocate(betainput%CFREQ)
    if(betainput%imag_cfrequencies_in_input) deallocate(betainput%IMCFREQ)


  end subroutine BETAresponse_driver


  !> \brief Get parallel and perpendicular averaged 1st hyperpolarizabilities.
  !> (See definition on http://folk.uio.no/michalj/node21.html or inside code).
  !> \author Kasper Kristensen                                                                
  !> \date 2010-08
  subroutine Get_averaged_betas(BetaPara,BetaPerp,beta,DipoleUnit)
    implicit none
    !> Parallel avaraged beta
    complex(realk), intent(inout) :: BetaPara
    !> Perpendicular avaraged beta
    complex(realk), intent(inout) :: BetaPerp
    !> Full beta tensor (3,3,3)
    complex(realk), intent(in)    :: beta(3,3,3)
    !> Dipole vector normalized to unity
    real(realk), intent(in)       :: DipoleUnit(3)
    integer :: i,j

    ! Beta parallel and beta perpendiular are for example defined on:
    ! http://folk.uio.no/michalj/node21.html

    ! Beta parallel:
    ! BetaPara = 1/5 * SUM_{i,j} ( beta_{jii} + beta_{iji} + beta_{iij} ) * DipoleUnit(j)
    ! where the i and j summations run independently over x,y,z components.
    BetaPara=(0E0_realk,0E0_realk)
    do i=1,3
       do j=1,3
          BetaPara = BetaPara + &
               & ( beta(j,i,i)+beta(i,j,i)+beta(i,i,j) )*DipoleUnit(j)
       end do
    end do
    BetaPara = (1E0_realk/5E0_realk)*BetaPara


    ! Beta perpendicular:
    ! BetaPerp = 1/5 * SUM_{i,j} ( 2*beta_{jii} - 3*beta_{iji} + 2*beta_{iij} ) * DipoleUnit(j)
    ! where the i and j summations run independently over x,y,z components.
    BetaPerp=(0E0_realk,0E0_realk)
    do i=1,3
       do j=1,3
          BetaPerp = BetaPerp + &
               & ( 2E0_realk*beta(j,i,i) - 3E0_realk*beta(i,j,i) + 2E0_realk*beta(i,i,j) )*DipoleUnit(j)
       end do
    end do
    BetaPerp = (1E0_realk/5E0_realk)*BetaPerp



  end subroutine Get_averaged_betas




  !> \brief Driver for calculating 2nd hyperpolarizabilites for frequencies defined in input.
  !>  If a given frequency is not specified in the input, it is set to zero.
  !> \author Kasper Kristensen                                                                
  !> \date 2010-08
  subroutine GAMMAresponse_driver(molcfg,F,D,S,gammainput)
    implicit none
    !> Info on molecule needed by solver and integral programs
    type(rsp_molcfg),intent(inout)    :: molcfg
    !> Unperturbed Fock matrix
    type(matrix), intent(in)          :: F
    !> Unperturbed density matrix
    type(matrix), intent(in)          :: D
    !> Unperturbed overlap matrix
    type(matrix), intent(in)          :: S
    !> Contains gamma input
    type(GAMMAinputItem),intent(inout)   :: GAMMAinput
    !> Dimension of response vector
    integer :: response_size
    !> Vector for collecting response results          
    complex(realk), allocatable :: rsp_results_vec(:,:)
    type(RspFunc), allocatable :: MyRspFunc(:)
    complex(realk), allocatable :: gamma(:,:,:,:,:)
    integer :: i, j, k, l, m, counter, nfreq
    character(len=1) :: comp(3)
    complex(realk), allocatable :: GammaPara(:), GammaPerp(:)
    integer, pointer :: lupri
    logical :: ImagFreq

    lupri => molcfg%lupri

    ! 1. Determine number of different frequencies and put frequencies into response function type
    ! ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

    ! Check that input is consistent.
    ! Code quite ugly, but it probably pays off to do this check...

    ! If present in input - number of real and imaginary frequencies must be identical.

    ! BFREQ and IMBFREQ
    if(gammainput%real_bfrequencies_in_input &
         & .and. gammainput%imag_bfrequencies_in_input) then
       if(gammainput%nbfreq /= gammainput%nimbfreq) then
          write(lupri,*) 'Error in GAMMAresponse_driver: The number of real and imaginary frequencies &
               & specified for operator B in the input file must be equal!'
          CALL lsQUIT('The number of real and imaginary frequencies  for operator B&
               & specified in the input file must be equal!',lupri)
       end if
    end if

    ! CFREQ and IMCFREQ
    if(gammainput%real_cfrequencies_in_input &
         & .and. gammainput%imag_cfrequencies_in_input) then
       if(gammainput%ncfreq /= gammainput%nimcfreq) then
          write(lupri,*) 'Error in GAMMAresponse_driver: The number of real and imaginary frequencies &
               & specified for operator C in the input file must be equal!'
          CALL lsQUIT('The number of real and imaginary frequencies  for operator C&
               & specified in the input file must be equal!',lupri)
       end if
    end if

    ! DFREQ and IMDFREQ
    if(gammainput%real_dfrequencies_in_input &
         & .and. gammainput%imag_dfrequencies_in_input) then
       if(gammainput%ndfreq /= gammainput%nimdfreq) then
          write(lupri,*) 'Error in GAMMAresponse_driver: The number of real and imaginary frequencies &
               & specified for operator D in the input file must be equal!'
          CALL lsQUIT('The number of real and imaginary frequencies  for operator D&
               & specified in the input file must be equal!',lupri)
       end if
    end if


    ! Also - if specified - the number of frequencies for operator B
    ! and operator C must be the same
    if(gammainput%real_bfrequencies_in_input &
         & .and. gammainput%real_cfrequencies_in_input) then
       if(gammainput%nbfreq /= gammainput%ncfreq) then
          write(lupri,*) 'Error in GAMMAresponse_driver: The number of frequencies specfied &
               &for operator B and C must be equal!'
          CALL lsQUIT('Error in GAMMAresponse_driver: The number of frequencies specfied &
               &for operator B and C must be equal!',lupri)
       end if
    end if

    ! Number of frequencies for operator B and operator D must be the same
    if(gammainput%real_bfrequencies_in_input &
         & .and. gammainput%real_dfrequencies_in_input) then
       if(gammainput%nbfreq /= gammainput%ndfreq) then
          write(lupri,*) 'Error in GAMMAresponse_driver: The number of frequencies specfied &
               &for operator B and D must be equal!'
          CALL lsQUIT('Error in GAMMAresponse_driver: The number of frequencies specfied &
               &for operator B and D must be equal!',lupri)
       end if
    end if


    ! Set number of frequencies equal to the number of frequencies for operator B or C or D
    ! (in case not all of them has been specified by the input)
    nfreq = max(gammainput%nbfreq,gammainput%ncfreq,gammainput%ndfreq)
    allocate(MyRspFunc(nfreq))
    allocate(GammaPara(nfreq))
    allocate(GammaPerp(nfreq))


    ! Initialize response functions
    do i=1,nfreq
       call initialize_RspFunc(MyRspFunc(i))
    end do

    ! If no frequencies have been defined, only calculate static polarizability 
    ! (omegaA = omegaB = omegaC = omegaD = 0).
    SetFrequencies: do i=1,nfreq

       ! Real frequencies operator B
       if(gammainput%real_bfrequencies_in_input) then
          MyRspFunc(i)%freq(2)=(1E0_realk,0E0_realk)*gammainput%bfreq(i)
       else
          MyRspFunc(i)%freq(2)=(0E0_realk,0E0_realk)
       end if

       ! Imaginary frequencies operator B
       if(gammainput%imag_bfrequencies_in_input) then
          MyRspFunc(i)%freq(2)=MyRspFunc(i)%freq(2) + (0E0_realk,1E0_realk)*gammainput%imbfreq(i)
       end if


       ! Real frequencies operator C
       if(gammainput%real_cfrequencies_in_input) then
          MyRspFunc(i)%freq(3)=(1E0_realk,0E0_realk)*gammainput%cfreq(i)
       else
          MyRspFunc(i)%freq(3)=(0E0_realk,0E0_realk)
       end if

       ! Imaginary frequencies operator C
       if(gammainput%imag_cfrequencies_in_input) then
          MyRspFunc(i)%freq(3)=MyRspFunc(i)%freq(3) + (0E0_realk,1E0_realk)*gammainput%imcfreq(i)
       end if

       ! Real frequencies operator D
       if(gammainput%real_dfrequencies_in_input) then
          MyRspFunc(i)%freq(4)=(1E0_realk,0E0_realk)*gammainput%dfreq(i)
       else
          MyRspFunc(i)%freq(4)=(0E0_realk,0E0_realk)
       end if

       ! Imaginary frequencies operator D
       if(gammainput%imag_dfrequencies_in_input) then
          MyRspFunc(i)%freq(4)=MyRspFunc(i)%freq(4) + (0E0_realk,1E0_realk)*gammainput%imdfreq(i)
       end if


       ! omegaA = -omegaB -omegaC -omegaD
       MyRspFunc(i)%freq(1) = - MyRspFunc(i)%freq(2) &
            & - MyRspFunc(i)%freq(3) - MyRspFunc(i)%freq(4) 

    end do SetFrequencies


    ! Any imaginary frequencies?
    if(gammainput%imag_bfrequencies_in_input .or. &
         &  gammainput%imag_cfrequencies_in_input .or. &
         &  gammainput%imag_dfrequencies_in_input) then
       ImagFreq=.true.
    else
       ImagFreq=.false. 
    end if


    ! 2. Set parameters other than frequencies
    ! ''''''''''''''''''''''''''''''''''''''''

    ! Allocate response result vector (3x3x3x3=81) and gamma tensor (3,3,3,3) for each frequency.
    response_size = 3*3*3*3
    allocate(rsp_results_vec(nfreq,response_size))
    allocate(gamma(nfreq,3,3,3,3))

    SetOtherParameters: do i=1,nfreq
       ! Order of response function=4 (cubic response)
       MyRspFunc(i)%order = 4

       ! Operators = (EL,EL,EL,EL)
       MyRspFunc(i)%code(1) = 'EL  '
       MyRspFunc(i)%code(2) = 'EL  '
       MyRspFunc(i)%code(3) = 'EL  '
       MyRspFunc(i)%code(4) = 'EL  '

       ! Calculate all 81 tensor components (xxxx,xxxy,..., zzzz) by default
       MyRspFunc(i)%first_comp(1) = 1
       MyRspFunc(i)%first_comp(2) = 1
       MyRspFunc(i)%first_comp(3) = 1
       MyRspFunc(i)%first_comp(4) = 1
       MyRspFunc(i)%dims(1) = 3
       MyRspFunc(i)%dims(2) = 3
       MyRspFunc(i)%dims(3) = 3
       MyRspFunc(i)%dims(4) = 3

       ! Zero response function and polarizabilities
       rsp_results_vec(i,:)=(0E0_realk,0E0_realk)
       gamma(i,:,:,:,:)=(0E0_realk,0E0_realk)

    end do SetOtherParameters


    ! 3. Calculate cubic response functions
    ! '''''''''''''''''''''''''''''''''''''

    ! Cubic response function <<EL;EL,EL,EL>> for each frequency.
    CalculateCubicResponseFunc: do i=1,nfreq

       ! Calculate <<EL;EL,EL,EL>> response functions for all 81 tensor components
       ! using 2n+1 rule
       call cubic_response(molcfg,F,D,S,MyRspFunc(i), &
            & rsp_results_vec(i,1:response_size),response_size)

       ! Gamma equals MINUS the cubic response function, and the cubic response
       ! function results are arranged in the vector rsp_results_vec as follows:
       !> (1,1,1,1) (1,1,1,2) (1,1,1,3) (1,1,2,1) (1,1,2,2) (1,1,2,3) (1,1,3,1) (1,1,3,2) (1,1,3,3)    
       !> (1,2,1,1) (1,2,1,2) (1,2,1,3) (1,2,2,1) (1,2,2,2) (1,2,2,3) (1,2,3,1) (1,2,3,2) (1,2,3,3)    
       !> (1,3,1,1) (1,3,1,2) (1,3,1,3) (1,3,2,1) (1,3,2,2) (1,3,2,3) (1,3,3,1) (1,3,3,2) (1,3,3,3)    
       !> (2,1,1,1) (2,1,1,2) (2,1,1,3) (2,1,2,1) (2,1,2,2) (2,1,2,3) (2,1,3,1) (2,1,3,2) (2,1,3,3)    
       !> (2,2,1,1) (2,2,1,2) (2,2,1,3) (2,2,2,1) (2,2,2,2) (2,2,2,3) (2,2,3,1) (2,2,3,2) (2,2,3,3)    
       !> (2,3,1,1) (2,3,1,2) (2,3,1,3) (2,3,2,1) (2,3,2,2) (2,3,2,3) (2,3,3,1) (2,3,3,2) (2,3,3,3)    
       !> (3,1,1,1) (3,1,1,2) (3,1,1,3) (3,1,2,1) (3,1,2,2) (3,1,2,3) (3,1,3,1) (3,1,3,2) (3,1,3,3)    
       !> (3,2,1,1) (3,2,1,2) (3,2,1,3) (3,2,2,1) (3,2,2,2) (3,2,2,3) (3,2,3,1) (3,2,3,2) (3,2,3,3)    
       !> (3,3,1,1) (3,3,1,2) (3,3,1,3) (3,3,2,1) (3,3,2,2) (3,3,2,3) (3,3,3,1) (3,3,3,2) (3,3,3,3)

       ! Thus, for each frequency (i) gamma may be calculated in the following way:
       counter=0
       do j=1,3
          do k=1,3
             do l=1,3
                do m=1,3
                   counter=counter+1
                   gamma(i,j,k,l,m) = -rsp_results_vec(i,counter)
                end do
             end do
          end do
       end do

    end do CalculateCubicResponseFunc




    ! 4. Get isotropically averaged gamma parallel and gamma perpendicular
    ! ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

    ! For each frequency i, get averaged gammas (see subroutine Get_averaged_gammas for definitions).
    do i=1,nfreq
       call Get_averaged_gammas(GammaPara(i),GammaPerp(i),gamma(i,1:3,1:3,1:3,1:3))
    end do


    ! 5. Print gamma tensor results
    ! '''''''''''''''''''''''''''''

    ! Character component vector
    comp(1) = 'X'
    comp(2) = 'Y'
    comp(3) = 'Z'


    write(lupri,*) 
    write(lupri,*) 
    write(lupri,*) 
    write(lupri,'(1X,A)') '*********************************************************************'
    write(lupri,'(1X,A)') '*        SECOND HYPERPOLARIZABILITY TENSOR RESULTS (in a.u.)        *'
    write(lupri,'(1X,A)') '*********************************************************************'
    write(lupri,*) 
    write(lupri,*)
    write(lupri,*) "The isotropically averaged parallel and perpendicular gamma values"
    write(lupri,*)  "are calculated in the following way:"
    write(lupri,*) 
    write(lupri,*) "GammaParallel      = 1/15 * SUM_{i,j} &
         &( gamma_{jiij} + gamma_{jiji} + gamma_{jjii} )"
    write(lupri,*) "GammaPerpendicular = 1/15 * SUM_{i,j} &
         &( 2*gamma_{jiij} - gamma_{jjii} )"
    write(lupri,*) 
    write(lupri,*) "where the i and j summations run independently over x,y,z."
    write(lupri,*) 
    write(lupri,*) 
    write(lupri,*) 
    write(lupri,*) 
    PrintFrequencyLoop: do i=1,nfreq

       if(ImagFreq) then
          write(lupri,'(1X,A,g12.5,A,g12.5,A)')  'Frequency B = ', &
               & real(MyRspFunc(i)%freq(2)), " + ", imag(MyRspFunc(i)%freq(2)), "i"
          write(lupri,'(1X,A,g12.5,A,g12.5,A)')  'Frequency C = ', &
               & real(MyRspFunc(i)%freq(3)), " + ", imag(MyRspFunc(i)%freq(3)), "i"
          write(lupri,'(1X,A,g12.5,A,g12.5,A)')  'Frequency D = ', &
               & real(MyRspFunc(i)%freq(4)), " + ", imag(MyRspFunc(i)%freq(4)), "i"
          write(lupri,*) '======================================'
       else
          write(lupri,'(1X,A,g12.5)')  'Frequency B = ', &
               & real(MyRspFunc(i)%freq(2))
          write(lupri,'(1X,A,g12.5)')  'Frequency C = ', &
               & real(MyRspFunc(i)%freq(3))
          write(lupri,'(1X,A,g12.5)')  'Frequency D = ', &
               & real(MyRspFunc(i)%freq(4))
          write(lupri,*) '========================'
       end if


       write(lupri,*) 
       if(ImagFreq) then
          write(lupri,*) "Components         Real part"
          write(lupri,*) "''''''''''''''''''''''''''''"
       else
          write(lupri,*) "Components         2nd hyp."
          write(lupri,*) "'''''''''''''''''''''''''''"
       end if

       j_loopREAL: do j=1,3
          k_loopREAL: do k=1,3
             l_loopREAL: do l=1,3
                m_loopREAL: do m=1,3
                   write(lupri,'(4X,4A,6X,g18.8)') comp(j),comp(k),&
                        &comp(l),comp(m),real(gamma(i,j,k,l,m))
                end do m_loopREAL
             end do l_loopREAL
          end do k_loopREAL
       end do j_loopREAL

       write(lupri,*) 
       write(lupri,'(a,g18.8)') 'GammaParallel      =', real(GammaPara(i))
       write(lupri,'(a,g18.8)') 'GammaPerpendicular =', real(GammaPerp(i))

       if(ImagFreq) then
          write(lupri,*) 
          write(lupri,*) 
          write(lupri,*) 'Components         Imag part'
          write(lupri,*) "'''''''''''''''''''''''''''"

          j_loopIMAG: do j=1,3
             k_loopIMAG: do k=1,3
                l_loopIMAG: do l=1,3
                   m_loopIMAG: do m=1,3
                      write(lupri,'(4X,4A,6X,g18.8)') comp(j),comp(k),&
                           &comp(l),comp(m),imag(gamma(i,j,k,l,m))
                   end do m_loopIMAG
                end do l_loopIMAG
             end do k_loopIMAG
          end do j_loopIMAG

          write(lupri,*) 
          write(lupri,'(a,g18.8)') 'GammaParallel      =', imag(GammaPara(i))
          write(lupri,'(a,g18.8)') 'GammaPerpendicular =', imag(GammaPerp(i))

       end if

       write(lupri,*) 
       write(lupri,*) 
       write(lupri,*) 
       write(lupri,*) 
       write(lupri,*) 

    end do PrintFrequencyLoop


    deallocate(MyRspFunc)
    deallocate(GammaPara)
    deallocate(GammaPerp)
    deallocate(rsp_results_vec)
    deallocate(gamma)

    if(gammainput%real_bfrequencies_in_input) deallocate(gammainput%BFREQ)
    if(gammainput%imag_bfrequencies_in_input) deallocate(gammainput%IMBFREQ)
    if(gammainput%real_cfrequencies_in_input) deallocate(gammainput%CFREQ)
    if(gammainput%imag_cfrequencies_in_input) deallocate(gammainput%IMCFREQ)
    if(gammainput%real_dfrequencies_in_input) deallocate(gammainput%DFREQ)
    if(gammainput%imag_dfrequencies_in_input) deallocate(gammainput%IMDFREQ)


  end subroutine GAMMAresponse_driver



  !> \brief Get parallel and perpendicular averaged 2nd hyperpolarizabilities.
  !> (See definition on http://folk.uio.no/michalj/node21.html or inside code).
  !> \author Kasper Kristensen                                                                
  !> \date 2010-08
  subroutine Get_averaged_gammas(GammaPara,GammaPerp,gamma)
    implicit none
    !> Parallel avaraged gamma
    complex(realk), intent(inout) :: GammaPara
    !> Perpendicular avaraged gamma
    complex(realk), intent(inout) :: GammaPerp
    !> Full gamma tensor (3,3,3,3)
    complex(realk), intent(in)    :: gamma(3,3,3,3)
    integer :: i,j

    ! Gamma parallel and gamma perpendiular are for example defined on:
    ! http://folk.uio.no/michalj/node21.html

    ! Gamma parallel:
    ! GammaPara = 1/15 * SUM_{i,j} ( gamma_{jiij} + gamma_{jiji} + gamma_{jjii} )
    ! where the i and j summations run independently over x,y,z components.
    GammaPara=(0E0_realk,0E0_realk)
    do i=1,3
       do j=1,3
          GammaPara = GammaPara + &
               & gamma(j,i,i,j)+gamma(j,i,j,i)+gamma(j,j,i,i)
       end do
    end do
    GammaPara = (1E0_realk/15E0_realk)*GammaPara


    ! Gamma perpendicular:
    ! GammaPerp = 1/15 * SUM_{i,j} ( 2*gamma_{jiij} - gamma_{jjii}) 
    ! where the i and j summations run independently over x,y,z components.
    GammaPerp=(0E0_realk,0E0_realk)
    do i=1,3
       do j=1,3
          GammaPerp = GammaPerp + &
               & 2E0_realk*gamma(j,i,i,j) - gamma(j,j,i,i)
       end do
    end do
    GammaPerp = (1E0_realk/15E0_realk)*GammaPerp



  end subroutine Get_averaged_gammas



  !> \brief Driver for calculating transition density matrices and storing them
  !> in rsp_eq_sol(i)%D in rsp_equations.f90 for the number of excited states
  !> defined by cfg_rsp_nexcit in decomp.
  !> \author Kasper Kristensen                                                                
  !> \date 2010-08
  subroutine calculate_and_store_transition_density_matrices(molcfg,F,D,S)
    implicit none
    !> Info on molecule needed by solver and integral programs
    type(rsp_molcfg),intent(inout)    :: molcfg
    !> Unperturbed Fock matrix
    type(matrix), intent(in)          :: F
    !> Unperturbed density matrix
    type(matrix), intent(in)          :: D
    !> Unperturbed overlap matrix
    type(matrix), intent(in)          :: S
    integer :: nexci,i


    ! 1. Set number of excitation energies
    ! ''''''''''''''''''''''''''''''''''''
    ! (Defined by .NEXCIT in input)
    nexci = molcfg%decomp%cfg_rsp_nexcit


    ! 2. Determine transition density matrices for the nexci excitations
    ! ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

    ! The transition density matrices are calculated as:
    !
    ! Dx(i) = D*S*X(i) - X(i)*S*D
    !
    ! for (i=1,nexci), where X(i) is a solution to eigenvalue problem                              
    !
    ! E2*X(i) = ExEnergy(i)*S2*X(i)
    !
    ! for each excitation vector X(i).
    ! Then Dx(i) is saved in rsp_eq_sol(i)%D in rsp_equations.f90.
    ! At the same time the corresponding excitation energies are
    ! stored in rsp_eq_sol(i)%fld(1)%freq as complex frequencies
    ! where the imaginary component is zero
    ! When this has been done the subroutines in response_driver.f90 automatically
    ! reads the transition density matrices from rsp_eq_sol(i)%D when needed. 
    ! (See subroutine transition_moment_density_matrix in response_driver.f90 for more details)
    call transition_moment_density_matrix(molcfg,F,D,S,nexci)

  end subroutine Calculate_and_store_transition_density_matrices

  !> \brief Driver for freeing the transition density matrices in
  !> rsp_eq_sol(i)%D in rsp_equations.f90 for the number of excited states
  !> defined by cfg_rsp_nexcit in decomp.
  !> \author Kasper Kristensen
  !> \date 2010-08
  subroutine free_transition_density_matrices(molcfg)
    implicit none
    type(rsp_molcfg),target,intent(in) :: molcfg
    integer :: i

    if (allocated(rsp_eq_sol)) then
       ! Run through rsp_eq_sol in rsp_equations.f90, freeing all 'EXCI' densities.
       do i=1,size(rsp_eq_sol)
          if (.not.associated(rsp_eq_sol(i)%mol,molcfg)) cycle
          if (size(rsp_eq_sol(i)%fld)/=1) cycle
          if (rsp_eq_sol(i)%fld(1)%label/='EXCI') cycle
          nullify(rsp_eq_sol(i)%mol)
          deallocate(rsp_eq_sol(i)%fld)
          rsp_eq_sol(i)%D = 0
          rsp_eq_sol(i)%refc = 0
       end do
    end if

  end subroutine Free_transition_density_matrices

  !> \brief Get the excitation energies
  !> \author Thomas Kjaergaard
  !> \date 2012-12
  subroutine get_excitation_energies(ExcitE,nexcit)
    implicit none
    !> number of Excitation Energies
    integer :: nexcit
    !> Excitation Energies
    real(realk) :: ExcitE(nexcit)
    !
    integer :: i
    do i=1,nexcit
       ExcitE(i) = real(rsp_eq_sol(i)%fld(1)%freq)
    enddo
  end subroutine Get_excitation_energies

  !> \brief Driver for calculating one-photon absorption transition momemnts.
  !> \author Kasper Kristensen                                                                
  !> \date 2010-08
  subroutine OPAresponse_driver(molcfg,F,D,S)
    implicit none
    !> Info on molecule needed by solver and integral programs
    type(rsp_molcfg),intent(inout)    :: molcfg
    !> Unperturbed Fock matrix
    type(matrix), intent(in)          :: F
    !> Unperturbed density matrix
    type(matrix), intent(in)          :: D
    !> Unperturbed overlap matrix
    type(matrix), intent(in)          :: S
    integer :: nexci
    !> Vector for collecting response results          
    complex(realk), allocatable :: rsp_results_vec(:,:)
    real(realk), allocatable :: OscillatorStrength(:), ExEnergies(:), trans_moments(:,:)
    type(RspFunc), allocatable :: MyRspFunc(:) 
    integer :: i,j, nELcomponents
    integer, pointer :: lupri

    lupri => molcfg%lupri

    ! 1. Set number of excitation energies and allocate response functions
    ! ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    ! (Defined by .NEXCIT in input or set to 1 by default)
    nexci = molcfg%decomp%cfg_rsp_nexcit

    allocate(MyRspFunc(nexci))
    allocate(OscillatorStrength(nexci))
    allocate(ExEnergies(nexci))
    ! 3 transition moments (x,y,z) for each excitation energy
    allocate(trans_moments(nexci,3))
    allocate(rsp_results_vec(nexci,3))



    ! We assume now that the transition density matrices has been calculated using
    ! the Calculate_and_store_transition_density_matrices subroutine.
    ! Therefore these matrices are stored in rsp_eq_sol(i)%D in rsp_equations.f90,
    ! while the corresponding excitation energies are
    ! stored in rsp_eq_sol(i)%fld(1)%freq as complex frequencies
    ! where the imaginary component is zero.


    ! 3. Setup response functions (=transition moments)
    ! '''''''''''''''''''''''''''''''''''''''''''''''''

    SetupTransMomentStructure: do i=1,nexci

       ! Initialize response function
       call initialize_RspFunc(MyRspFunc(i))

       ! Order of response function=2 (linear response)
       MyRspFunc(i)%order = 2

       ! Operators = (EL,EXCI): Transition moment <0|EL|n>
       MyRspFunc(i)%code(1) = 'EL  '
       MyRspFunc(i)%code(2) = 'EXCI'

       ! Calculate all 3 transition moments, <0|ELx|n>, <0|ELy|n>, <0|ELz|n>.
       ! Note that the number of components for the second (EXCI) operator is only 1,
       ! and the number of response functions is nexci.
       ! This is because EACH excitation energy requires its own response function.
       MyRspFunc(i)%first_comp(1) = 1
       MyRspFunc(i)%first_comp(2) = 1
       MyRspFunc(i)%dims(1) = 3
       MyRspFunc(i)%dims(2) = 1

       ! We are interested in a transition moment.
       MyRspFunc(i)%trans_moment = .true.

       ! The transition moment corresponds to a SINGLE residue,
       ! therefore the order of the transition moment is 1.
       MyRspFunc(i)%trans_moment_order = 1

       ! For MyRspFunc(i) we calculate transition moments for the i'th excited state.
       MyRspFunc(i)%ExNumber1 = i

       ! The frequecy of the response function equals the excitation energy
       ! which is stored in rsp_eq_sol(i)%fld(1)%freq as a complex frequency
       ! where the imaginary component is zero (see point 2 above).
       ! Note: The excitation energy is the frequency for the SECOND operator (EXCI)
       MyRspFunc(i)%freq(2) = rsp_eq_sol(i)%fld(1)%freq
       ! Frequncy 1 = -Frequency 2
       MyRspFunc(i)%freq(1) = - MyRspFunc(i)%freq(2) 
       ! For later convience we also have a vector of REAL excitation energies
       ExEnergies(i) = real(MyRspFunc(i)%freq(2))

       ! Zero response function results
       rsp_results_vec(i,:)=(0E0_realk,0E0_realk)

    end do SetupTransMomentStructure


    ! Calculate transition momemnts for all 3 (x,y,z) components of electric dipole.
    nELcomponents = 3



    ! 4. Calculate one-photon transition moments
    ! ''''''''''''''''''''''''''''''''''''''''''

    ! Linear response function <<EL;EXCI>> for EACH excitation energy.
    CalculateTransMoment: do i=1,nexci


       ! Calculate <<EL;EXCI>> response functions for all 3 (x,y,z) components of EL.
       ! The linear_response subroutine now reads the transition density matrices from
       ! rsp_eq_sol(i)%D in rsp_equations.f90.
       call linear_response(molcfg,F,D,S,MyRspFunc(i), &
            & rsp_results_vec(i,1:nELcomponents),nELcomponents)

       ! The transition moments are now stored in rsp_results_vec as
       ! complex numbers with zero imaginary part. We are only interested in the real part.
       ! Oscillator strength for excited state i:
       ! strength = 2/3 * Excitation energy * SUM_{j=x,y,z} <0|ELj|i>
       OscillatorStrength(i) = 0E0_realk
       do j=1,nELcomponents
          trans_moments(i,j) = real(rsp_results_vec(i,j))
          OscillatorStrength(i) = OscillatorStrength(i) + trans_moments(i,j)**2
       end do
       OscillatorStrength(i) = (2E0_realk/3E0_realk)*ExEnergies(i)*OscillatorStrength(i)

    end do CalculateTransMoment



    ! 5. Print transition moments and excitation energies
    ! '''''''''''''''''''''''''''''''''''''''''''''''''''
    call print_OPA_results(lupri,nexci,ExEnergies, trans_moments, OscillatorStrength)


    deallocate(MyRspFunc)
    deallocate(OscillatorStrength)
    deallocate(ExEnergies)
    deallocate(trans_moments)
    deallocate(rsp_results_vec)


  end subroutine OPAresponse_driver


  !> \brief Routine for printing one-photon absorption results.
  !> \author Kasper Kristensen                                                                
  !> \date 2010-08
  subroutine print_OPA_results(lupri,nexci,ExEnergies, trans_moments, OscillatorStrength)
    implicit none
    !> File unit number for LSDALTON.OUT
    integer, intent(in) :: lupri
    !> Number of excitation energies
    integer, intent(in) :: nexci
    !> Vector containing excitation energies
    real(realk), intent(in) :: ExEnergies(nexci)
    !> Vector containing transition moments
    real(realk), intent(in) :: trans_moments(nexci,3)
    !> Vector containing oscillator strengths
    real(realk), intent(in) :: Oscillatorstrength(nexci)
    integer :: i

    write(lupri,*) 
    write(lupri,*) 
    write(lupri,*) 
    write(lupri,*) 
    write(lupri,'(2X,A)') '********************************************************&
         &**********************'
    write(lupri,'(2X,A)') '*                   ONE-PHOTON ABSORPTION RESULTS (in a.u.)&
         &                  *'
    write(lupri,'(2X,A)') '********************************************************&
         &**********************'
    write(lupri,*)
    write(lupri,*) 
    write(lupri,*) 
    write(lupri,*) '     Excitation              Transition Dipole Moments      &
         &         Oscillator'
    write(lupri,*) '      Energies            x               y               z  &
         &         Strengths'
    write(lupri,*) '===================================================================&
         &============'

    PrintExcitationLoop: do i=1,nexci

       write(lupri,'(5f16.8)') ExEnergies(i), trans_moments(i,1:3), OscillatorStrength(i)

    end do PrintExcitationLoop

    write(lupri,*) 
    write(lupri,*)
    write(lupri,*) 
    write(lupri,*) 


  end subroutine print_OPA_results





  !> \brief Driver for calculating two-photon absorption transition moments.
  !> TPA is determined from a residue of the quadratic response function.
  !> See e.g. JCP 82, 3235 (1985), Eqs. (3.18-3.20).
  !> \author Kasper Kristensen                                                                
  !> \date 2010-08
  subroutine TPAresponse_driver(molcfg,F,D,S,TPAinput)
    implicit none
    !> Info on molecule needed by solver and integral programs
    type(rsp_molcfg),intent(inout)    :: molcfg
    !> Unperturbed Fock matrix
    type(matrix), intent(in)          :: F
    !> Unperturbed density matrix
    type(matrix), intent(in)          :: D
    !> Unperturbed overlap matrix
    type(matrix), intent(in)          :: S
    !> Contains TPA input
    type(TPAinputItem),intent(inout)   :: TPAinput
    integer :: nexci, state, size_rsp
    logical :: AllStates
    logical, allocatable :: doState(:)
    !> Vector for collecting response results          
    complex(realk), allocatable :: rsp_results_vec(:,:)
    real(realk), allocatable :: ExEnergies(:), TPA_trans_moments(:,:,:), &
         & tpa_lin(:), tpa_circ(:)
    type(RspFunc), allocatable :: MyRspFunc(:) 
    integer :: i,j,k,l, counter, nELcomponents
    character(len=1) :: comp(3)
    integer, pointer :: lupri

    lupri => molcfg%lupri
    comp(1) = 'X'
    comp(2) = 'Y'
    comp(3) = 'Z'



    ! 1. Set number of excited states
    ! '''''''''''''''''''''''''''''''
    ! (Defined by .NEXCIT in input or set to 1 by default)

    nexci = molcfg%decomp%cfg_rsp_nexcit
    if(TPAinput%specific_states_in_input) then
       AllStates = .false.  ! only calculate TPA for specific states
    else 
       AllStates = .true.   ! calculate TPA for ALL excited states from 1 to nexci
    end if



    ! 2. Determine which particular excited states to calculate TPA for
    ! '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

    ! Logical for keeping track of which states to address
    allocate(doState(nexci))

    doState(:)=.false.
    if(AllStates) then ! TPA for all excited states from 1 to nexci
       size_rsp=nexci
       doState(:) = .true.
    else ! Only states specified by input
       size_rsp=TPAinput%tpa_nexci
       do i=1,TPAinput%tpa_nexci
          state=TPAinput%ExStates(i)
          doState(state) = .true.
       end do
    end if


    allocate(MyRspFunc(size_rsp))
    allocate(ExEnergies(size_rsp))
    allocate(tpa_lin(size_rsp))
    allocate(tpa_circ(size_rsp))
    ! 3x3=9 transition moments (xx,xy,...,zz) for each excitation energy
    allocate(TPA_trans_moments(size_rsp,3,3))
    allocate(rsp_results_vec(size_rsp,9))



    ! We assume now that the transition density matrices for state 1 to nexci
    ! has been calculated using
    ! the Calculate_and_store_transition_density_matrices subroutine.
    ! Therefore these matrices are stored in rsp_eq_sol(i)%D in rsp_equations.f90,
    ! while the corresponding excitation energies are
    ! stored in rsp_eq_sol(i)%fld(1)%freq as complex frequencies
    ! where the imaginary component is zero.
    ! (It is checked in configuration.f90 that the TPA input is consistent).


    ! 3. Setup response functions (=TPA transition moments)
    ! '''''''''''''''''''''''''''''''''''''''''''''''''''''

    ! The TPA transition moment may be formulated as a single residue
    ! of a quadratic response function, see e.g. JCP 82, 3235 (1985), Eqs. (3.18-3.20).
    ! In our framework this corresponds to a response function of the
    ! form <<EL;EL,EXCI>>, where the frequency of the EXCI operator
    ! equals the corresponding excitation energy and the frequencies
    ! for the EL operators both equal -1/2*excitation energy.
    ! (Assuming that the two photons have the same frequency).


    ! Use j as counter in case not all excited states are requested
    j = 0
    SetupTPAResponseFunc: do i=1,nexci

       StateRequested: if(doState(i)) then
          j=j+1

          ! Initialize response function
          call initialize_RspFunc(MyRspFunc(j))

          ! Order of response function=3 (quadratic response function residue)
          MyRspFunc(j)%order = 3

          ! Operators = (EL,EL,EXCI) yield TPA transition moment (see above)
          MyRspFunc(j)%code(1) = 'EL  '
          MyRspFunc(j)%code(2) = 'EL  '
          MyRspFunc(j)%code(3) = 'EXCI'

          ! Calculate all 9 transition moments (xx, xy, xz, ... , zz),
          ! i.e. 3 components for the EL-operators, whereas
          ! the number of components for the third (EXCI) operator is only 1.
          ! Note that EACH excitation energy requires its own response function.
          MyRspFunc(j)%first_comp(1) = 1
          MyRspFunc(j)%first_comp(2) = 1
          MyRspFunc(j)%first_comp(3) = 1
          MyRspFunc(j)%dims(1) = 3
          MyRspFunc(j)%dims(2) = 3
          MyRspFunc(j)%dims(3) = 1

          ! We are interested in a transition moment.
          MyRspFunc(j)%trans_moment = .true.

          ! The transition moment corresponds to a SINGLE residue,
          ! therefore the order of the transition moment is 1.
          MyRspFunc(j)%trans_moment_order = 1

          ! If all excited states are requested MyRspFunc(j) simply refers to
          ! excited state number j, whereas if specific excited states
          ! have been defined in the input, then MyRspFunc(j) refers to the
          ! excited state stored in TPAinput%ExStates(j).
          if(AllStates) then
             MyRspFunc(j)%ExNumber1 = j
          else
             MyRspFunc(j)%ExNumber1 = TPAinput%ExStates(j)
          end if

          ! As descibed above:
          ! Frequency 1 = -1/2*excitation energy
          ! Frequency 2 = -1/2*excitation energy
          ! Frequency 3 = excitation energy

          ! The excitation energy,
          ! is stored in rsp_eq_sol(X)%fld(1)%freq where
          ! X refers to the position of the excited state now
          ! stored in MyRspFunc(j)%ExNumber1.
          ! Note that the frequencies stored in rsp_eq_sol(X)%fld(1)%freq are complex
          ! frequencies where the imaginary component is zero.
          MyRspFunc(j)%freq(3) = rsp_eq_sol(MyRspFunc(j)%ExNumber1)%fld(1)%freq
          MyRspFunc(j)%freq(2) = -0.5E0_realk*MyRspFunc(j)%freq(3)
          MyRspFunc(j)%freq(1) = - MyRspFunc(j)%freq(2) - MyRspFunc(j)%freq(3)
          ! For later convience we also have a vector of REAL excitation energies
          ExEnergies(j) = real(MyRspFunc(j)%freq(3))

          ! Zero response function
          rsp_results_vec(j,:)=(0E0_realk,0E0_realk)
       end if StateRequested

    end do SetupTPAResponseFunc


    ! Calculate transition moments for all 3x3 (x,y,z) X (x,y,z) components of electric dipoles.
    nELcomponents = 3*3



    ! 4. Calculate two-photon transition moments
    ! ''''''''''''''''''''''''''''''''''''''''''

    ! Quadratic TPA response function <<EL;EL,EXCI>> for each excited state given in doState.
    j=0
    CalculateTransMoment: do i=1,nexci

       CalculateState: if(doState(i)) then
          j=j+1

          ! Calculate <<EL;EL,EXCI>> response functions for all 3*3 components of EL.
          ! The quadratic_response subroutine now reads the transition density matrices from
          ! rsp_eq_sol(j)%D in rsp_equations.f90.
          call quadratic_response_2nplus1(molcfg,F,D,S,MyRspFunc(j), &
               & rsp_results_vec(j,1:nELcomponents),nELcomponents)

          ! The transition moments are now stored in rsp_results_vec as
          ! complex numbers with zero imaginary part in the order:
          ! (x,x) (x,y) (x,z) (y,x) (y,y) (y,z) (z,x) (z,y) (z,z)
          ! For each excitation energy j, we now store the results in the
          ! TPA_trans_moments array as real numbers.
          counter=0
          do k=1,3
             do l=1,3
                counter=counter+1
                TPA_trans_moments(j,k,l) = real(rsp_results_vec(j,counter))
             end do
          end do

          ! Get isotropically averaged TPA for linearly and for circularly polarized light.
          call Get_averaged_TPA(TPA_trans_moments(j,1:3,1:3),tpa_lin(j),tpa_circ(j))


       end if CalculateState

    end do CalculateTransMoment




    ! 5. Print TPA transition moments and excitation energies
    ! '''''''''''''''''''''''''''''''''''''''''''''''''''''''

    write(lupri,*) 
    write(lupri,*) 
    write(lupri,*) 
    write(lupri,*) 
    write(lupri,'(2X,A)') '*************************************************************&
         &**********************'
    write(lupri,'(2X,A)') '*                     TWO-PHOTON ABSORPTION RESULTS (in a.u.)&
         &                     *'
    write(lupri,'(2X,A)') '*************************************************************&
         &**********************'
    write(lupri,*)
    write(lupri,*) 
    write(lupri,*) 

    write(lupri,*) 'Isotropically averaged TPA'
    write(lupri,*) '**************************'
    write(lupri,*) 
    write(lupri,*) 'tpa_averaged =  1/30 [ F*tpaF + G*tpaG + H*tpaH ]'
    write(lupri,*) "tpaF = sum_{j,k=x,y,z} S(j,j)*S(k,k)"
    write(lupri,*) "tpaG = sum_{j,k=x,y,z} S(j,k)*S(j,k)"
    write(lupri,*) "tpaH = sum_{j,k=x,y,z} S(j,k)*S(k,j)"

    write(lupri,*) 
    write(lupri,*) "Linearly polarized light  :  F=G=H=2"
    write(lupri,*) "Circularly polarized light:  F=-2; G=H=3"

    write(lupri,*) 
    write(lupri,*) 
    write(lupri,*) ' State No.   Exc. Energy        Linear pol.       Circ. pol. '
    write(lupri,*) '============================================================='


    j=0
    PrintAveragedTPA: do i=1,nexci

       StateWasCalculated1: if(doState(i)) then
          j=j+1


          write(lupri,'(I7,3X,3g18.8)') MyRspFunc(j)%ExNumber1, ExEnergies(j), &
               & tpa_lin(j), tpa_circ(j)


       end if StateWasCalculated1

    end do PrintAveragedTPA

    write(lupri,*) 
    write(lupri,*) 
    write(lupri,*) 
    write(lupri,*) 



    j=0
    PrintTPAtensor: do i=1,nexci

       StateWasCalculated2: if(doState(i)) then
          j=j+1

          write(lupri,'(1X,A,I6)') ' Two-photon transition tensor S for state no. ', &
               & MyRspFunc(j)%ExNumber1
          write(lupri,'(1X,A,I6)') '======================================================'
          write(lupri,*) "                Ex                Ey                Ez"
          write(lupri,'(1X,A, 3g18.8)') "Ex    ", TPA_trans_moments(j,1,1:3)
          write(lupri,'(1X,A, 3g18.8)') "Ey    ", TPA_trans_moments(j,2,1:3)
          write(lupri,'(1X,A, 3g18.8)') "Ez    ", TPA_trans_moments(j,3,1:3)
          write(lupri,*) 
          write(lupri,*) 
          write(lupri,*) 
          write(lupri,*) 

       end if StateWasCalculated2

    end do PrintTPATensor




    deallocate(MyRspFunc)
    deallocate(ExEnergies)
    deallocate(TPA_trans_moments)
    deallocate(rsp_results_vec)
    deallocate(doState)
    deallocate(tpa_lin)
    deallocate(tpa_circ)
    if(tpainput%specific_states_in_input) deallocate(tpainput%ExStates)


  end subroutine TPAresponse_driver



  !> \brief Driver for calculating isotropically averaged 
  !> two-photon absorption cross sections for linearly
  !> and circularly polarized laser light from
  !> TPA transition moment components,
  !> see e.g. Journal of Chemical Physics, 53, 29 (1970).
  !> \author Kasper Kristensen                                                                
  !> \date 2010-08
  subroutine Get_averaged_tpa(tpa,tpa_lin,tpa_circ)
    implicit none
    !> TPA transition dipole moment tensor
    real(realk), intent(in)  :: tpa(3,3)
    !> Absorption cross section for linearly polarized laster light
    real(realk), intent(inout)  :: tpa_lin
    !> Absorption cross section for circularly polarized laster light
    real(realk), intent(inout)  ::  tpa_circ
    real(realk)  :: F, G, H,tpaF, tpaG,tpaH
    integer      :: j,k


    ! The isotropically averaged TPA spectrum is given by:
    !
    ! tpa_av =  1/30 [ F*tpaF + G*tpaG + H*tpaH ]    (*)
    !
    ! tpaF = sum_{j,k} tpa(j,j)*tpa(k,k)
    ! tpaG = sum_{j,k} tpa(j,k)*tpa(j,k)
    ! tpaH = sum_{j,k} tpa(j,k)*tpa(k,j)
    !
    ! where j and k runs independently over the X,Y,Z components, 
    ! and where F and G are numbers which depend on the polarization of the photons,
    ! see Journal of Chemical Physics, 53, 29 (1970).
    ! 
    tpa_lin=0E0_realk
    tpa_circ=0E0_realk


    ! 1. Calculate tpaF, tpaG, and tpaH
    ! '''''''''''''''''''''''''''''''''

    tpaF = 0E0_realk
    do j=1,3
       do k=1,3
          tpaF = tpaF + tpa(j,j)*tpa(k,k)
       end do
    end do

    tpaG = 0E0_realk
    do j=1,3
       do k=1,3
          tpaG = tpaG + tpa(j,k)*tpa(j,k)
       end do
    end do

    tpaH = 0E0_realk
    do j=1,3
       do k=1,3
          tpaH = tpaH + tpa(j,k)*tpa(k,j)
       end do
    end do


    ! 2. Linearly polarized light
    ! '''''''''''''''''''''''''''

    ! According to Journal of Chemical Physics, 53, 29 (1970)
    ! for linearly polarized light we have:
    F=2; G=2; H=2
    tpa_lin = (1E0_realk/30E0_realk) * (F*tpaF + G*tpaG + H*tpaH)


    ! 3. Circularly polarized light
    ! '''''''''''''''''''''''''''''

    ! According to Journal of Chemical Physics, 53, 29 (1970)
    ! for circularly polarized light we have:
    F=-2; G=3; H=3
    tpa_circ = (1E0_realk/30E0_realk) * (F*tpaF + G*tpaG + H*tpaH)


  end subroutine Get_averaged_tpa



  !> \brief Driver for calculating damped TPA for frequencies specified in input.
  !> Damped TPA is obtained from the complex cubic response function but with some terms omitted.
  !> Therefore we do not use the general response driver, but instead calls the solver
  !> explicitly inside the DTPAresponse_driver.
  !> More details are given inside the subroutine and in JCP, 214104 (2011).
  !> \author Kasper Kristensen                                                                
  !> \date 2010-09
  subroutine DTPAresponse_driver(molcfg,F,D,S,dtpainput)
    ! The damped TPA function for the (a,b,c,d)th
    ! component is given by:
    !
    ! tpa_abcd = - Im [tr(RHS_ab * Xcd)]
    !
    ! where RHS_ab is the second order damped right-hand side matrix for the electric dipole
    ! operators mu_a and mu_b, and Xcd is the second order damped response parameter matrix
    ! obtained by solving:
    !
    ! [E2 - (wcd+i*gamma)S2] Xcd = RHS_cd
    !
    ! where gamma is the damping factor specified in the input.
    ! Note that RHS_ab, RHS_cd and Xcd are all COMPLEX.
    !
    ! We only consider TPA where the two photons have the same energy,
    ! which means that the frequencies are related as follows:
    !
    ! wa = wb = -w;      wc = wd = w;    wcd = 2*w
    !
    ! where w is the (positive) optical frequency of interest.
    ! Due to the frequency relations RHS_ab and RHS_cd are related by (for a=c and b=d):
    !
    ! RHS_ab = - (RHS_cd)^dagger
    !
    ! Therefore, assuming a=c and b=d, we can write tpa_abcd in the form:
    !
    ! tpa_abcd = - Im [tr(RHS_ab * Xcd)]
    !          =   Im[ tr( RHS_cd^dagger * Xcd) ]
    !          =   Im[ tr{ ( Re(RHS_cd)^T - i Im(RHS_cd)^T ) ( Re(Xcd) + i Im(Xcd) ) } ]
    !          = -tr{ Im(RHS_cd)^T Re(Xcd) } + tr{ Re(RHS_cd)^T Im(Xcd) }
    !          = -dot_product[ Im(RHS_cd) Re(Xcd) ] + dot_product[ Re(RHS_cd) Im(Xcd) ]    (**)
    !
    ! The damped one-photon spectrum is also calculated because
    ! the damped first-order equations are solved in the TPA procedure anyway.

    implicit none
    !> Info on molecule needed by solver and integral programs
    type(rsp_molcfg),intent(inout)    :: molcfg
    !> Unperturbed Fock matrix
    type(matrix), intent(in)          :: F
    !> Unperturbed density matrix
    type(matrix), intent(in)          :: D
    !> Unperturbed overlap matrix
    type(matrix), intent(in)          :: S
    !> Contains damped TPA input
    type(DTPAinputItem),intent(inout)       :: DTPAinput
    type(matrix)               :: RHS(6), Xreal(1), Ximag(1), D1(3)
    integer                    :: I,J,K,nfreq
    real(realk)                :: tmp1, tmp2, normRHS,w(1),gamma_save
    real(realk), allocatable   :: freqs(:), one_photon(:,:),tpa(:,:,:)
    integer, pointer           :: lupri
    type(decompItem), pointer  :: decomp


    ! Initialize stuff and sanity checks
    ! **********************************
    lupri => molcfg%lupri
    decomp => molcfg%decomp

    write(lupri,*) 
    write(lupri,*) 'Starting damped TPA calculation...'
    write(lupri,*) 

    ! Sanity check: Meaningful number of frequencies
    if(dtpainput%nfreq < 1) then
       write(lupri,*) 'Error in DTPAresponse_driver: The number frequencies is smaller than 1!'
       write(lupri,*) 'One-photon frequencies must be &
            &specified in the input by the .OPFREQ keyword.'
       write(lupri,*) 'Note that the two photons ALWAYS have the same frequency.'
       write(lupri,*) 'For example, to calculate three points in the damped TPA spectrum'
       write(lupri,*) 'where each one-photon frequency is 0.07, 0.1 and 0.15, respectively,'
       write(lupri,*) '(meaning that excitation energies of 0.14, 0.2, and 0.3 are probed)'
       write(lupri,*) 'the following response input is used:'
       write(lupri,*) '**RESPONSE'
       write(lupri,*) '*DAMPED_TPA'
       write(lupri,*) '.OPFREQ'
       write(lupri,*) '3'
       write(lupri,*) '0.07  0.1  0.15'
       CALL lsQUIT('The number frequencies is smaller than one! See output file for example &
            &of input',lupri)
    end if

    ! Set frequencies
    nfreq = dtpainput%nfreq
    allocate(freqs(nfreq))
    allocate(one_photon(nfreq,3))
    allocate(tpa(nfreq,6,6))
    freqs(1:nfreq) = dtpainput%freq(1:nfreq)
    ! Saving gamma stored in solver
    gamma_save=molcfg%solver%rsp_gamma
    molcfg%solver%rsp_gamma=dtpainput%gamma

    ! From now on the "gamma" stored in molcfg%solver%rsp_gamma
    ! is used as the damping parameter in ALL equations!
    ! Therefore all our frequencies are considered to be real
    ! and then i*gamma is added to them inside the solver to make them complex.
    ! At the end we set molcfg%solver%rsp_gamma=gamma_save again.

    if(dtpainput%gamma_specified) then
       write(lupri,*) 'Damping parameter gamma is set to:', dtpainput%gamma
    else
       write(lupri,*) 'Damping parameter gamma set to default value:', dtpainput%gamma
       write(lupri,*) 'Note: Gamma value can be changed by .GAMMA keyword &
            &  under *DAMPED_TPA.'
    end if


    ! Initializing complex solver
    if (molcfg%solver%rsp_cmplxnew .or. molcfg%solver%rsp_cpp) then
       call rsp_sym_complex_init(1, 1, 1, 1, 1)
    else 
       call rsp_complex_init(1, 1, 1, 1, 1)
    endif
    ! Initialize solution matrices
    Xreal(1) = 0*D
    call mat_ensure_alloc(Xreal(1), only_alloc=.true.)
    Ximag(1) = 0*D
    call mat_ensure_alloc(Ximag(1), only_alloc=.true.)


    ! Loop over frequencies and components and calculate damped TPA
    ! *************************************************************


    ! Note: We use the following component/number identification:
    ! 1: XX,   2: YY,   3: ZZ,   4: XY=YX,   5: XZ=ZX,   6: YZ=ZY    (*)

    ! Also note: The indices of tpa are as follows:
    ! First index: Frequency index
    ! Second index: Index for the first two TPA components according to (*) 
    ! Third index: Index for the last two TPA components according to (*) 
    tpa(:,:,:) = 0E0_realk


    ! Loop over frequencies specified in input
    FrequencyLoop: do I = 1,nfreq

       Write(LUPRI,'( "Starting damped TPA calculation for frequency ", I3, " of", &
            & I3, " frequencies: ", F9.6 )') I, nfreq, freqs(I)

       ! Complex second order RHS is constructed.
       call build_2nd_order_RHS_for_DTPA(molcfg,F,D,S,RHS,D1,freqs(I),lupri)
       ! RHS now generally contains six matrices corresponding to the index convention in (*).

       ! Complex second order RHS is constructed.
       ! Calculate one-photon absorption for frequency I for the components specified in input.
       call calculate_one_photon_absorption(molcfg,D,D1,one_photon(I,:))

       ! Complex second order RHS is constructed.
       ! 2nd order wcd frequency to be used for solving complex 2nd order eq.
       w(1) = 2E0_realk*freqs(I)


       ! Determine isotropically averaged damped TPA.
       ! ********************************************

       ! According to the description in the subroutine get_averaged_dtpa
       ! 12 different components are needed for the averaging:
       !
       ! XXXX,  YYYY,  ZZZZ,  XXYY,  YYXX,  XXZZ,  ZZXX,  YYZZ,  ZZYY,  XYXY,  XZXZ,  YZYZ
       ! (1 1)  (2 2)  (3 3)  (1 2)  (2 1)  (1 3)  (3 1)  (2 3)  (3 2)  (4 4)  (5 5)  (6 6)
       !
       ! where the corresponding indices according to (*) are given below. 
       ! The remaining components, e.g. XXYZ (1 6), are not calculated.
       ! Note the symmetry of the tpa tensor for the following permuations:
       ! tpa(J,K,L,M) = tpa(K,J,M,L)
       ! and the lack of symmetry for the permutations:
       ! tpa(J,J,K,K) /= tpa(K,K,J,J)
       ! because tpa is complex.


       !call print_mat(F,label='F')


       ! Calculating damped TPA according to (**) for the first 9 components:
       ! XXXX, XXYY, XXZZ, YYXX, YYYY, YYZZ, ZZXX, ZZYY, ZZZZ (in that particular order)
       ! '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
       ! Complex second order RHS is constructed.
       do K = 1,3

          normRHS = dot(RHS(K),RHS(K))
          normRHS = sqrt(normRHS)

          ! Problems with solver if normRHS = 0.
          if (normRHS > 1e-9_realk) then

             ! Complex second order RHS is constructed.
             ! Solve 2nd order equations for Kth component
             if (molcfg%solver%rsp_cmplxnew .or. molcfg%solver%rsp_cpp) then
                call rsp_sym_complex_solver(molcfg, F, D, S, 1, &
                     & (/mat_get_part(RHS(K), imag=.false.)/),  &
                     & w, 1, Xreal(1:1), Ximag(1:1), .true.,    &
                     & gdi=(/mat_get_part(RHS(K), imag=.true.)/))
             else
                call rsp_complex_solver(molcfg, F, D, S, 1,    &
                     & (/mat_get_part(RHS(K), imag=.false.)/), &
                     & w, 1, Xreal(1:1), Ximag(1:1), .true.,   &
                     & gdi=(/mat_get_part(RHS(K), imag=.true.)/))
             endif

             ! Complex second order RHS is constructed.
             do J =1,3
                ! Calculate the damped TPA strength according to (**).
                ! Multiply by 2 for closed shell systems, must be generalized to include open shell systems.
                tmp1 = -2E0_realk*dot(mat_get_part(RHS(J), imag=.true. ), Xreal(1))  !im(RHS)
                tmp2 =  2E0_realk*dot(mat_get_part(RHS(J), imag=.false.), Ximag(1))  !re(RHS)
                tpa(I,J,K) = tmp1 + tmp2
             enddo
          else ! Set tpa=0 if normRHS for component K is zero.
             do J=1,3
                tpa(I,J,K) = 0E0_realk
             enddo
          endif

       enddo


       ! Calculating damped TPA according to (**) above for the components:
       ! XYXY, XZXZ, YZYZ (in that order)
       ! ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
       do K=4,6

          normRHS = dot(RHS(K),RHS(K))
          normRHS = sqrt(normRHS)

          ! Problems with solver if RHS = 0
          if (normRHS > 1e-9_realk) then

             if (molcfg%solver%rsp_cmplxnew .or. molcfg%solver%rsp_cpp) then
                call rsp_sym_complex_solver(molcfg, F, D, S, 1, &
                     & (/mat_get_part(RHS(K), imag=.false.)/),  &  ! Real part of RHS
                     & w(1), 1, Xreal(1), Ximag(1), .true.,     &
                     & gdi=(/mat_get_part(RHS(K), imag=.true.)/))  ! Imaginary part of RHS
             else
                call rsp_complex_solver(molcfg, F, D, S, 1,    &
                     & (/mat_get_part(RHS(K), imag=.false.)/), &  ! Real part of RHS
                     & w(1), 1, Xreal(1), Ximag(1), .true.,    &
                     & gdi=(/mat_get_part(RHS(K), imag=.true.)/))  ! Imaginary part of RHS
             endif

             ! Calculate the damped TPA strength according to (**).
             ! Multiply by 2 for closed shell systems, must be generalized to include open shell systems.
             tmp1 = -2E0_realk*dot(mat_get_part(RHS(K), imag=.true. ), Xreal(1))  ! im(RHS)
             tmp2 =  2E0_realk*dot(mat_get_part(RHS(K), imag=.false.), Ximag(1))  ! re(RHS)
             tpa(I,K,K) = tmp1 + tmp2
          else
             tpa(I,K,K) = 0E0_realk
          endif

       enddo

    enddo FrequencyLoop


    ! Print the calculated damped TPA spectrum and the one-photon absorption spectrum
    ! *******************************************************************************
    ! The isotropical averaging is also carried out in the print_tpa routine
    call print_dtpa(lupri,nfreq,freqs,tpa,one_photon)


    ! Restore gamma
    molcfg%solver%rsp_gamma = gamma_save


    ! Free stuff
    ! **********
    D1(:)=0
    RHS(:)=0
    Xreal(1)=0; Ximag(1)=0
    deallocate(dtpainput%freq)
    deallocate(freqs)
    deallocate(tpa)
    deallocate(one_photon)


  end subroutine DTPARESPONSE_DRIVER



  !> \brief Build 6 (XX,XY,XZ,YY,YZ,ZZ) complex second order RHS matrices required for damped TPA calculation. 
  !> More details are given inside the subroutine build_2nd_order_RHS_for_DTPA.
  !> \author Kasper Kristensen                                                                
  !> \date 2010-09
  subroutine build_2nd_order_RHS_for_DTPA(molcfg,F,D,S,RHS,D1,freq,lupri)

    ! According to "Journal of Chemical Physics, 129, 214108 (2008)"
    ! (and after projecting out the non-redundant occupied/virtual parts)
    ! the second order RHS with respect to perturbations c and d become:
    !
    ! RHS_cd = Projector [  FPcd*D*S - S*D*FPcd
    !                     + Fc*Dd*S - S*Dd*Fc
    !                     + Fd*Dd*S - S*Dc*Fc ]
    !
    ! FPcd : second order (cd) derivative of the Fock matrix, where ONLY the first order parameters are kept.
    ! Fc/Fd: Total first order derivative of the Fock matrix with respect to perturbation c/d.
    ! Dc/Dd: Total first order derivative of the density matrix with respect to perturbation c/d.
    ! The projector acts by projecting out the PQ and QP components to
    ! avoid redundancies (P=DS; Q=1-P):
    !
    ! Projector(M) = P^dagger*M*Q +  Q^dagger*M*P
    !
    ! Or, equivalently, by the two operations:
    !
    !    1. M = M*D*S-S*D*M
    !    2. M = M*D*S-S*D*M
    !
    ! Outputs:
    ! RHS: Complex second order density matrix (6 components)
    ! D1(3): Complex first order perturbed density matrix for frequency "freq".

    implicit none
    !> structure containing the molecule, integral and solver settings
    type(rsp_molcfg), intent(inout)        :: molcfg
    !> Unperturbed Fock matrix
    type(matrix), intent(in)               :: F
    !> Unperturbed density matrix
    type(matrix), intent(in)               :: D
    !> Unperturbed overlap matrix
    type(matrix), intent(in)               :: S
    !> Complex second order right-hand side matrix for electric field perturbations
    type(matrix), intent(inout)            :: RHS(6)
    !> First order complex density matrix for electric field perturbations
    type(matrix), intent(inout)            ::  D1(3)
    !> Frequency for operators in second order RHS (they both have the same frequency)
    real(realk),intent(in)                 :: freq
    !> File unit number for LSDALTON.OUT
    integer, intent(inout)                 :: lupri
    type(matrix)                           :: F1(3),DS,SD
    integer                                :: J,K,comp
    complex(realk)                         :: freq_complex


    ! Real part of complex frequency is input frequency 
    ! Note: Here the frequency has only a real part.
    ! Inside the solver i*gamma is added to the frequency to make it real
    freq_complex = freq*(1E0_realk,0E0_realk) 

    ! Calculate perturbed first order density matrix D1 and Fock matrix F1.
    call pert_dens(molcfg, S, (/'EL'/), (/3/), (/D/), (/F/), D1, F1, &
         & comp=(/1/), freq=(/freq_complex/))


    DS = D*S
    SD = S*D

    do comp=1,6
       ! Get first-order indices (J,K) corresponding to 
       ! second order index "comp" according to the ordering:
       ! 1: XX,   2: YY,   3: ZZ,   4: XY,   5: XZ,   6: YZ   
       ! E.g. for comp=5 we get J=1 and K=3
       call determine_J_and_K(J,K,comp,lupri) 
       call calculate_component_of_RHS_for_DTPA(molcfg, D1(J), D1(K), F1(J), F1(K), &
            & J, K, RHS(comp),D, S, DS, SD)
    enddo


    ! Free matrices
    do J=1,3
       F1(J)=0
    enddo

    DS=0  
    SD=0

  end subroutine build_2nd_order_RHS_for_DTPA



  !> \brief Get first-order indices (J,K) corresponding to second order index "comp" according to the ordering:
  !> <br>
  !> 1: XX,   2: YY,   3: ZZ,   4: XY,   5: XZ,   6: YZ   
  !> <br>
  !> E.g. for comp=5 we get J=1 and K=3
  !> \author Kasper Kristensen                                                                
  !> \date 2010-09
  subroutine determine_J_and_K(J,K,comp,lupri)

    !> Combined second order component (XX,XY,XZ,YY,YZ,ZZ)
    integer, intent(in)       :: comp
    !> First first-order component (X,Y,Z)
    integer, intent(inout)    :: J
    !> Second first-order component (X,Y,Z)
    integer, intent(inout)    :: K
    !> File unit for LSDALTON.OUT
    integer, intent(in)       :: lupri

    select case(comp)

    case(1) ! XX-component
       J = 1; K = 1
    case(2) ! YY-component
       J = 2; K = 2
    case(3) ! ZZ-component
       J = 3; K = 3
    case(4) ! XY-component = YX-component
       J = 1; K = 2
    case(5) ! XZ-component = ZX-component
       J = 1; K = 3
    case(6) ! YZ-component = ZY-component
       J = 2; K = 3

    case default
       call LSQUIT('determine_J_and_K: Index for comp must be between 1 and 6!',lupri)
    end select

  end subroutine determine_J_and_K



  !> \brief Calculates a particular component of the second order RHS matrix 
  !> required for damped TPA as described
  !> in subroutine build_2nd_order_RHS_for_DTPA.
  !> \author Kasper Kristensen                                                                
  !> \date 2010-09
  subroutine calculate_component_of_RHS_for_DTPA(molcfg, Dj, Dk, Fj, Fk, J, K, &
       & RHS,D,S,DS,SD)

    implicit none
    !> structure containing the molecule, integral and solver settings
    type(rsp_molcfg), intent(inout) :: molcfg
    !> Unperturbed density matrix
    type(matrix), intent(in)      :: D
    !> Unperturbed overlap matrix
    type(matrix), intent(in)      :: S
    !> First-order perturbed density matrix (component j in build_2nd_order_RHS_for_DTPA)
    type(matrix), intent(in)      :: Dj
    !> First-order perturbed density matrix (component k in build_2nd_order_RHS_for_DTPA)
    type(matrix), intent(in)      :: Dk
    !> First-order perturbed Fock matrix (component j in build_2nd_order_RHS_for_DTPA)
    type(matrix), intent(in)      :: Fj
    !> First-order perturbed Fock matrix (component k in build_2nd_order_RHS_for_DTPA)
    type(matrix), intent(in)      :: Fk
    !> D*S  (to avoid calculating this matrix for all components)
    type(matrix), intent(in)      :: DS
    !> S*D  (to avoid calculating this matrix for all components)
    type(matrix), intent(in)      :: SD
    !> Component (j,k) of RHS matrix in build_2nd_order_RHS_for_DTPA
    type(matrix),intent(inout)    :: RHS
    !> First component (X,Y, or Z) or second order RHS
    integer, intent(in)           :: J
    !> Second component (X,Y, or Z) or second order RHS
    integer, intent(in)           :: K

    type(matrix)                  :: FPcd(1,1), DPcd


    ! Calculate the second order "particular" density matrix, 
    ! where the second order response parameters are omitted.
    ! It corresponds to a projected form of the differentiated idempotency relation,
    ! see "Journal of Chemical Physics, 129, 214108 (2008)".
    if (J == K) then 
       DPcd = -2E0_realk*Dj*S*Dj
    else
       DPcd = -Dj*S*Dk - Dk*S*Dj
    endif
    DPcd = - DPcd + DPcd*SD  +DS*DPcd

    ! Second order particular Fock matrix FPcd containing only first order response parameters.
    ! (Comment: The terms containing Dj and Dk only contribute for DFT calculations,
    !           whereas the terms containing DPcd contribute for both HF and DFT.)
    call pert_fock(molcfg, (/'EL', 'EL'/), (/1,1/), (/D, Dj, Dk, DPcd /), FPcd, &
         & comp=(/J,K/))

    if(J == K) then ! Fewer matrix multiplications if J=K
       RHS= FPcd(1,1)*DS - SD*FPcd(1,1) + &
            & 2E0_realk*(Fj*Dj*S - S*Dj*Fj)
    else
       RHS = FPcd(1,1)*DS - SD*FPcd(1,1) + &
            & Fj*Dk*S - S*Dk*Fj + Fk*Dj*S - S*Dj*Fk
    endif

    ! Projecting RHS as described in subroutine build_2nd_order_RHS_for_DTPA.
    RHS = RHS*DS-SD*RHS
    RHS = RHS*DS-SD*RHS

    ! Free matrices
    Dpcd=0; FPcd(1,1)=0

  end subroutine calculate_component_of_RHS_for_DTPA


  !> \brief Calculate one-photon absorption with perturbed density matrix as input.
  !> Only used in damped TPA routine to supplement the damped TPA spectrum
  !> by a damped one-photon absorption spectrum at "no" extra cost.
  !> (We need to have special treatment instead of just calling the linear response driver 
  !> to avoid calculating the first-order density matrices twice.)
  !> \author Kasper Kristensen                                                                
  !> \date 2010-09
  subroutine calculate_one_photon_absorption(molcfg,D,D1,one_photon)

    implicit none
    !> structure containing the molecule, integral and solver settings
    type(rsp_molcfg), intent(inout) :: molcfg
    !> Unperturbed density matrix
    type(matrix), intent(in)   :: D
    !> First-order perturbed density matrices for electric field perturbation
    type(matrix), intent(in)   :: D1(3)
    !> One-photon absorption for the frequency associated with D1
    real(realk), intent(inout) :: one_photon(3)
    complex(realk)             :: linear_response_func(1)
    integer                    :: J

    ! One-photon absorption = - Im<< mu_a ; mu_b >> 
    ! (See "Journal of Chemical Physics, 131, 044112 (2009)")

    one_photon(:) = 0E0_realk

    do J=1,3
       ! Zero before rsp_oneint call, otherwise the components are added
       linear_response_func(1) = (0E0_realk,0E0_realk)
       call rsp_oneave(molcfg, D, (/'EL'/), (/D1(J)/), (/1/), &
            & linear_response_func(1:1), comp=(/J/))
       one_photon(J) = -imag(linear_response_func(1))
    enddo

  end subroutine calculate_one_photon_absorption


  !> \brief Print two-photon absorption spectra AND one-photon absorption spectra.
  !> The isotropical averaging of TPA and OPA components is also carried out here.
  !> \author Kasper Kristensen                                                                
  !> \date 2010-09
  subroutine print_dtpa(lupri,nfreq,freqs,tpa,one_photon)
    ! 
    ! If the isotropically averaged spectra are requested the individual components are also printed. 
    implicit none
    !> File unit for LSDALTON.OUT
    integer, intent(in)        :: lupri
    !> Number of different frequencies
    integer, intent(in)        :: nfreq
    !> One-photon frequencies
    real(realk), intent(in)    :: freqs(nfreq)
    !> Damped TPA strength components using labelling described inside subroutine
    real(realk),intent(in)     :: tpa(nfreq,6,6) 
    !> Damped OPA strength components
    real(realk),intent(in)     :: one_photon(nfreq,3)
    integer                    :: J,K
    real(realk)                :: tpa_lin(nfreq), tpa_circ(nfreq), op_av(nfreq)
    character*1                :: one_photon_label(3)
    character*2                :: labels(6)


    ! Labels according to the ordering:
    ! 1: XX,   2: YY,   3: ZZ,   4: XY=YX,   5: XZ=ZX,   6: YZ=ZY   
    labels(1) = 'XX'; labels(2) = 'YY'; labels(3) = 'ZZ'
    labels(4) = 'XY'; labels(5) = 'XZ'; labels(6) = 'YZ'

    ! One-photon labels
    one_photon_label(1) = 'X'; one_photon_label(2) = 'Y'; one_photon_label(3) = 'Z'

    ! Perform isotropic average for linear and circular polarizations.
    call get_averaged_dtpa(nfreq,tpa,tpa_lin,'linear')
    call get_averaged_dtpa(nfreq,tpa,tpa_circ,'circular')
    ! Calculate isotropically averaged one-photon absorption
    call get_averaged_opa(nfreq,one_photon, op_av)


    write(lupri,*) 
    write(lupri,*) 
    write(lupri,*) 
    write(lupri,*) 
    write(lupri,'(2X,A)') '*************************************************************&
         &**********************'
    write(lupri,'(2X,A)') '*                DAMPED TWO-PHOTON ABSORPTION RESULTS (in a.u.)&
         &                   *'
    write(lupri,'(2X,A)') '*************************************************************&
         &**********************'
    write(lupri,*)
    write(lupri,*) 
    write(lupri,*) '*** See Kristensen et al., J. Chem. Phys. 134, 214104 (2011) ***'
    write(lupri,*) 'In particular, the damped TPA components given below correspond to Eq. (108).'
    write(lupri,*) 'See also Eq. (48) for the corresponding expression in exact theory along with'
    write(lupri,*) 'a physical interpretation of damped TPA spectra.'
    write(lupri,*)
    write(lupri,*) 


    write(lupri,*) 'Isotropically averaged TPA'
    write(lupri,*) '**************************'
    write(lupri,*) 
    write(lupri,*) 'tpa_averaged =  1/30 [ F*tpaF + G*tpaG + H*tpaH ]'
    write(lupri,*) "tpaF = sum_{j,k=x,y,z} S(j,j)*S(k,k)"
    write(lupri,*) "tpaG = sum_{j,k=x,y,z} S(j,k)*S(j,k)"
    write(lupri,*) "tpaH = sum_{j,k=x,y,z} S(j,k)*S(k,j)"

    write(lupri,*) 
    write(lupri,*) "Linearly polarized light  :  F=G=H=2"
    write(lupri,*) "Circularly polarized light:  F=-2; G=H=3"

    write(lupri,*) 
    write(lupri,*) 
    write(lupri,*) '  Frequency            Linear pol.       Circ. pol.'
    write(lupri,*) '==================================================='
    call print_freq_vs_2strengths(lupri,nfreq,freqs,tpa_lin,tpa_circ)

    write(lupri,*) 
    write(lupri,*) 'Individual TPA components'
    write(lupri,*) "========================="
    write(lupri,*) 

    ! XXXX, XXYY, XXZZ, YYXX, YYYY, YYZZ, ZZXX, ZZYY, ZZZZ components
    do J=1,3
       do K=1,3
          write(lupri,*) '  Frequency      TPA component: ', labels(J), labels(K) 
          call print_freq_vs_strength(lupri,nfreq,freqs,tpa(:,J,K))
       enddo
    enddo

    ! XYXY, XZXZ, YZYZ components
    do J=4,6
       write(lupri,*) '  Frequency      TPA component: ', labels(J), labels(J) 
       call print_freq_vs_strength(lupri,nfreq,freqs,tpa(:,J,J))
    enddo

    ! One-photon spectrum
    write(lupri,*) 
    write(lupri,*) 
    write(lupri,*) 
    write(lupri,*) 
    write(lupri,'(2X,A)') '*************************************************************&
         &**********************'
    write(lupri,'(2X,A)') '*                DAMPED ONE-PHOTON ABSORPTION RESULTS (in a.u.)&
         &                   *'
    write(lupri,'(2X,A)') '*************************************************************&
         &**********************'
    write(lupri,*)
    write(lupri,*) 

    write(lupri,*) 
    write(lupri,*) 
    write(lupri,*) 'ISOTROPICALLY AVERAGED ONE-PHOTON ABSORPTION'
    write(lupri,*) "********************************************"
    write(lupri,*) '  Frequency       One-photon absorption'
    call print_freq_vs_strength(lupri,nfreq,freqs,op_av)

    write(lupri,*) 
    write(lupri,*) 'Individual OPA components'
    write(lupri,*) "========================="
    write(lupri,*) 

    do J = 1,3
       write(lupri,*) '  Frequency      OPA component: ', one_photon_label(J)
       call print_freq_vs_strength(lupri,nfreq,freqs,one_photon(:,J))
    enddo

    write(lupri,*)
    write(lupri,*)

  end subroutine print_dtpa


  !> \brief Carry out isotropical averaging of damped TPA components for
  !> either linearly or circularly polarized light.
  !> (Since the component labelling is special for damped TPA we separate this completely
  !>  from the averaging of standard TPA spectra although the underlying
  !>  averaging equations are of course the same.)
  !> More details inside subroutine.
  !> \author Kasper Kristensen                                                                
  !> \date 2010-09
  subroutine get_averaged_dtpa(nfreq,tpa,tpa_av,poltype)

    implicit none
    !> Number of frequencies
    integer, intent(in) :: nfreq
    !> Damped TPA components for each frequency using labelling described inside subroutine
    real(realk),intent(in)   :: tpa(nfreq,6,6)
    !> Isotropically averaged damped TPA for eah frequency
    real(realk),intent(inout)   :: tpa_av(nfreq)
    !> Which polarization type (linear or circular?)
    character(*), intent(in)    :: poltype
    real(realk)                 :: F, GH, tpaF, tpaGH
    integer                     :: I, J, K


    ! The isotropically averaged TPA spectrum is given by:
    !
    ! tpa_av =  1/30 [ F*tpaF + GH*tpaGH ]    (***)
    !
    ! tpaF = sum_{J,K} tpa(J,J,K,K)
    ! tpaF = sum_{J,K} tpa(J,K,J,K)
    !
    ! where J and K runs independently over the X,Y,Z components, 
    ! and where F and GH are numbers which depend on the polarization of the photons,
    ! see "Journal of Chemical Physics, 53, 29 (1970)".
    ! In fact, our 'F' corresponds to F in this article and 'GH' corresponds to G+H.
    ! 
    ! Two different polarization types are implemented:
    !
    ! poltype = "linear": Two linearly polarized photons with parallel propagation.
    ! poltype="circular": Two circularly polarized photons with parallel propagation.
    !
    ! According to (***) only 12 different components are needed for the averaging:
    ! XXXX,  YYYY,  ZZZZ,  XXYY,  YYXX,  XXZZ,  ZZXX,  YYZZ,  ZZYY,  XYXY,  XZXZ,  YZYZ
    ! (1 1)  (2 2)  (3 3)  (1 2)  (2 1)  (1 3)  (3 1)  (2 3)  (3 2)  (4 4)  (5 5)  (6 6)


    if(poltype == 'linear') then ! See article reference above
       F = 2.0E0_realk
       GH = 4.0E0_realk
    elseif(poltype == 'circular') then
       F = -2.0E0_realk
       GH = 6.0E0_realk
    else
       call LSquit('get_averaged_dtpa: Requested polarization case not implemented!',-1)
    endif

    tpa_av(:) = 0E0_realk

    do I = 1,nfreq

       ! Setting tpaF and tpaGH terms in (***) to zero.
       tpaF = 0E0_realk; tpaGH = 0E0_realk

       ! Loop over (JJ,KK) components (F-term):
       ! XXXX, XXYY, XXZZ, YYXX, YYYY, YYZZ, ZZXX, ZZYY, ZZZZ.
       do J=1,3
          do K=1,3
             tpaF = tpaF + tpa(I,J,K)
          enddo
       enddo

       ! Loop over (JK,JK) components (GH-term):
       ! XXXX, YYYY, ZZZZ
       ! XYXY=YXYX, XZXZ=ZXZX, YZYZ=ZYZY
       do J=1,6
          if(J<4) then ! XXXX, YYYY, ZZZZ
             tpaGH = tpaGH + tpa(I,J,J)
          else     ! XYXY=YXYX, XZXZ=ZXZX, YZYZ=ZYZY
             ! Multiply by 2 since XYXY=YXYX, XZXZ=ZXZX, YZYZ=ZYZY.
             tpaGH = tpaGH + 2E0_realk*tpa(I,J,J)
          endif
       enddo

       tpa_av(I) = (1E0_realk/30E0_realk) * (F*tpaF + GH*tpaGH)

    enddo

  end subroutine get_averaged_dtpa


  !> \brief Carry out isotropical averaging of one-photon absorption spectra
  !> \author Kasper Kristensen                                                                
  !> \date 2010-09
  subroutine get_averaged_opa(nfreq,one_photon, op_av)

    implicit none
    !> Number of frequencies
    integer, intent(in) :: nfreq
    ! XYZ components of one-photon absorption strengths
    real(realk), intent(in) ::  one_photon(nfreq,3)
    !> Isotropically averaged one-photon absorption strength
    real(realk)      :: op_av(nfreq)
    integer          :: I,J


    ! Isotropically averaged one-photon absorption:
    ! op_av = 1/3 * [ one_photon(x) + one_photon(y) + one_photon(z) ]

    op_av(:) = 0E0_realk

    do I = 1,nfreq

       do J = 1,3
          op_av(I) = op_av(I) + one_photon(I,J)
       enddo
       op_av(I) = (1E0_realk/3E0_realk) * op_av(I)

    enddo

  end subroutine get_averaged_opa



  !> \brief Prints a list of numbers: (frequency, absorption strength)
  !> \author Kasper Kristensen                                                                
  !> \date 2010-09
  subroutine print_freq_vs_strength(lupri,nfreq,freqs,abs_strength)
    implicit none
    !> File unit for LSDALTON.OUT
    integer,intent(in) :: lupri
    !> Number of frequencies
    integer, intent(in) :: nfreq
    !> frequencies
    real(realk),intent(in) :: freqs(nfreq)
    !> Absorption strengths
    real(realk), intent(in)  :: abs_strength(nfreq)
    integer         :: I


    do I =1,nfreq
       Write(LUPRI,'(g16.7, "   ", g16.7)') freqs(I), abs_strength(I)
    enddo

    write(lupri,*) 

  end subroutine print_freq_vs_strength


  !> \brief Prints a list of numbers: (frequency, absorption strength 1, absorption strength 2)
  !> Abs. strength 1 and 2 may e.g. refer to linearly and circularly polarized light.
  !> \author Kasper Kristensen                                                                
  !> \date 2010-09
  subroutine print_freq_vs_2strengths(lupri,nfreq,freqs,abs_strength1,abs_strength2)
    implicit none
    !> File unit for LSDALTON.OUT
    integer,intent(in) :: lupri
    !> Number of frequencies
    integer, intent(in) :: nfreq
    !> frequencies
    real(realk),intent(in) :: freqs(nfreq)
    !> Absorption strength 1
    real(realk), intent(in)  :: abs_strength1(nfreq)
    !> Absorption strength 2
    real(realk), intent(in)  :: abs_strength2(nfreq)
    integer         :: I


    do I =1,nfreq
       Write(LUPRI,'(g16.7, "   ", g16.7, "   ", g16.7)') freqs(I), abs_strength1(I), abs_strength2(I)
    enddo

    write(lupri,*) 

  end subroutine print_freq_vs_2strengths





  !> \brief Driver for calculating excited state gradient (ESG).
  !> ESG is determined by subtracting a part of the double residue of the quadratic response function
  !> from the ground state gradient.
  !> See e.g. JCP 82, 3235 (1985), Eq. (3.21) with mu_a --> dH/dR and m=q.
  !> \author Kasper Kristensen                                                                
  !> \date 2010-08
  subroutine ESGresponse_driver(molcfg,F,D,S,ESGinput,ESGgrad)
    implicit none
    !> Info on molecule needed by solver and integral programs
    type(rsp_molcfg), intent(inout)   :: molcfg
    !> Unperturbed Fock matrix
    type(matrix), intent(in)          :: F
    !> Unperturbed density matrix
    type(matrix), intent(in)          :: D
    !> Unperturbed overlap matrix
    type(matrix), intent(in)          :: S
    !> Contains ESG input
    type(ESGinputItem),intent(inout)   :: ESGinput
    !> desired ESG output
    real(realk),optional   :: ESGgrad(3,molcfg%setting%molecule(1)%p%natoms)
    integer :: nexci, state, size_rsp, natoms
    integer, pointer :: lupri
    logical :: AllStates
    logical, pointer :: doState(:)
    complex(realk), allocatable :: rsp_results_vec(:,:), GSG(:)
    real(realk), pointer :: ExEnergies(:), ESG(:,:,:)
    type(RspFunc), allocatable :: MyRspFunc(:) 
    integer :: i,j,k,l, counter, totdim


    ! Number of atoms and file unit number for LSDALTON.OUT
    natoms = molcfg%setting%molecule(1)%p%natoms
    lupri => molcfg%lupri




    ! 1. Set number of excited states
    ! '''''''''''''''''''''''''''''''
    ! (Defined by .NEXCIT in input or set to 1 by default)

    nexci = molcfg%decomp%cfg_rsp_nexcit
    if(ESGinput%specific_states_in_input) then
       AllStates = .false.  ! only calculate ESG for specific states
    else 
       AllStates = .true.   ! calculate ESG for ALL excited states from 1 to nexci
    end if



    ! 2. Determine which particular excited states to calculate ESG for
    ! '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

    ! Logical for keeping track of which states to address
    call mem_alloc(doState,nexci)

    doState(:)=.false.
    if(AllStates) then ! ESG for all excited states from 1 to nexci
       size_rsp=nexci
       doState(:) = .true.
    else ! Only states specified by input
       size_rsp=ESGinput%esg_nexci
       do i=1,ESGinput%esg_nexci
          state=ESGinput%ExStates(i)
          doState(state) = .true.
       end do
    end if


    allocate(MyRspFunc(size_rsp))
    call mem_alloc(ExEnergies,size_rsp)
    ! 3*natoms coordinates describing ESG for each excitation energy
    totdim=3*natoms
    call mem_alloc(ESG,size_rsp,natoms,3)
    allocate(rsp_results_vec(size_rsp,totdim))
    allocate(GSG(totdim))





    ! 3. Setup response functions
    ! '''''''''''''''''''''''''''

    ! We assume now that the transition density matrices for state 1 to nexci
    ! has been calculated using
    ! the Calculate_and_store_transition_density_matrices subroutine.
    ! Therefore these matrices are stored in rsp_eq_sol(i)%D in rsp_equations.f90,
    ! while the corresponding excitation energies are
    ! stored in rsp_eq_sol(i)%fld(1)%freq as complex frequencies
    ! where the imaginary component is zero.
    ! (It is checked in configuration.f90 that the ESG input is consistent).


    ! The ESG may be formulated in terms of a double residue
    ! of a quadratic response function, 
    ! see e.g. JCP 82, 3235 (1985), Eq. (3.21) with mu_a --> dH/dR, mu_b=mu_c=B, and m=q:
    !
    ! Double residue = - |<0|B|m>|^2 * { <m| dH/dR |m> - <0| dH/dR |0> }
    !
    ! Thus the excited state gradient <m| dH/dR |m> becomes:
    ! 
    ! <m| dH/dR |m> = <0| dH/dR |0> - { Double residue / |<0|B|m>|^2 } 
    !
    ! We do not care about the B operator since in a practical calculation
    ! we factorize the double residue of <<A;B,C>> as:
    ! 
    ! Double residue = <<A; EXCI,EXCI>> * <<B;EXCI>> * <<C;EXCI>>
    ! 
    ! where <<B; EXCI>> = <0|B|m>
    ! 
    ! Thus, we may calculate the excited state gradient as:
    !
    ! <m| dH/dR |m> = <0| dH/dR |0> - <<dH/dR;EXCI,EXCI>>   (*)
    !
    ! (Recall that perturbation-dependent basis sets are built into the 
    !  theory - and into the code - when we call the response routines
    !  in response_driver. Therefore we do not need to worry about this
    !  and we can use the exact theory response expression directly.)


    ! Use j as counter in case not all excited states are requested
    j = 0
    SetupESGResponseFunc: do i=1,nexci

       StateRequested: if(doState(i)) then
          j=j+1

          ! Initialize response function
          call initialize_RspFunc(MyRspFunc(j))

          ! Order of response function=3 (quadratic response function double residue)
          MyRspFunc(j)%order = 3

          ! Operators = (GEO,EXCI,EXCI) yield ESG (see above)
          MyRspFunc(j)%code(1) = 'GEO '
          MyRspFunc(j)%code(2) = 'EXCI'
          MyRspFunc(j)%code(3) = 'EXCI'

          ! Calculate all totdim=3*natoms components of the ESG.
          ! I.e. we have 3*natoms components for the GEO operator, whereas
          ! the number of components for the second and third (EXCI) operators is only 1.
          ! Note that EACH excitation energy requires its own response function.
          MyRspFunc(j)%first_comp(1) = 1
          MyRspFunc(j)%first_comp(2) = 1
          MyRspFunc(j)%first_comp(3) = 1
          MyRspFunc(j)%dims(1) = totdim
          MyRspFunc(j)%dims(2) = 1
          MyRspFunc(j)%dims(3) = 1

          ! We are interested in a double residue (transition moment).
          MyRspFunc(j)%trans_moment = .true.

          ! The transition moment corresponds to a DOUBLE residue,
          ! therefore the order of the transition moment is 2.
          MyRspFunc(j)%trans_moment_order = 2


          ! If all excited states are requested MyRspFunc(j) simply refers to
          ! excited state number j, whereas if specific excited states
          ! have been defined in the input, then MyRspFunc(j) refers to the
          ! excited state number stored in ESGinput%ExStates(j).
          ! Both excited states in the double residue are the same 
          ! since q=m in JCP 82, 3235 (1985), Eq. (3.21): ExNumber1=ExNumber2
          if(AllStates) then
             MyRspFunc(j)%ExNumber1 = j
             MyRspFunc(j)%ExNumber2 = j
          else
             MyRspFunc(j)%ExNumber1 = ESGinput%ExStates(j)
             MyRspFunc(j)%ExNumber2 = ESGinput%ExStates(j)
          end if

          ! Accoding to JCP 82, 3235 (1985), Eq. (3.21) the residue of interest has
          ! Frequency 2 = -excitation energy
          ! Frequency 3 = excitation energy
          ! Therefore Frequency 1 = 0

          ! The excitation energy
          ! is stored in rsp_eq_sol(X)%fld(1)%freq where
          ! X refers to the position of the excited state now
          ! stored in MyRspFunc(j)%ExNumber1.
          ! Note that the frequencies stored in rsp_eq_sol(X)%fld(1)%freq are complex
          ! frequencies where the imaginary component is zero.
          MyRspFunc(j)%freq(3) =  rsp_eq_sol(MyRspFunc(j)%ExNumber1)%fld(1)%freq
          MyRspFunc(j)%freq(2) = -rsp_eq_sol(MyRspFunc(j)%ExNumber1)%fld(1)%freq
          MyRspFunc(j)%freq(1) = 0E0_realk
          ! For later convience we also have a vector of REAL excitation energies
          ExEnergies(j) = real(MyRspFunc(j)%freq(3))

          ! Zero response function
          rsp_results_vec(j,:)=(0E0_realk,0E0_realk)

       end if StateRequested

    end do SetupESGResponseFunc




    ! 4. Calculate excited state gradient
    ! '''''''''''''''''''''''''''''''''''

    ! According to Eq (*) above we need the ground state gradient (GSG) to calculate the
    ! excited state gradient using a response approach.
    ! The 3*natoms gradient elements are stored in the GSG vector as complex
    ! entities with zero imaginary part (to keep the general structure of the response code).
    call get_first_order_property(molcfg,F,D,S,GSG,totdim,'GEO ')



    ! Response function  <<GEO;EXCI,EXCI>> [see Eq. (*) above] for each excited state given in doState.
    j=0
    CalculateTransMoment: do i=1,nexci

       CalculateState: if(doState(i)) then
          j=j+1
          ! Calculate <<GEO;EXCI,EXCI>> response functions for all 3*natoms components of GEO.
          ! The quadratic_response subroutine now reads the transition density matrices from
          ! rsp_eq_sol(j)%D in rsp_equations.f90.
          ! Note: We use the n+1 rule here!
          call quadratic_response_nplus1(molcfg,F,D,S,MyRspFunc(j), &
               & rsp_results_vec(j,1:totdim),totdim)

          ! The transition moments are now stored in rsp_results_vec as
          ! complex numbers with zero imaginary part in the order:
          ! Atom1(x), Atom1(y), Atom1(z), Atom2(x), ...
          ! For each excitation energy j, we calculate the ESG using Eq. (*) above.
          ! The ESG array is a REAL array and has the structure:
          ! (excitation,atom,xyz-component)
          counter=0
          do k=1,natoms
             do l=1,3
                counter=counter+1
                ESG(j,k,l) = real(GSG(counter)) - real(rsp_results_vec(j,counter))
             end do
          end do
       end if CalculateState

    end do CalculateTransMoment

    IF(present(ESGgrad))THEN
       IF(j.GT.1)call lsquit('ESG opt error j.gt.1',-1)
       do k=1,natoms
          do l=1,3
             ESGgrad(l,k) = ESG(1,k,l)
          end do
       end do
    ENDIF

    ! 5. Print ESG and excitation energies
    ! ''''''''''''''''''''''''''''''''''''

    write(lupri,*) 
    write(lupri,*) 
    write(lupri,*) 
    write(lupri,*) 
    write(lupri,'(2X,A)') '*************************************************************&
         &**********************'
    write(lupri,'(2X,A)') '*                     EXCITED STATE GRADIENT RESULTS (in a.u.)&
         &                    *'
    write(lupri,'(2X,A)') '*************************************************************&
         &**********************'
    write(lupri,*)
    write(lupri,*) 
    write(lupri,*) 


    j=0
    PrintESG: do i=1,nexci

       StateWasCalculated: if(doState(i)) then
          j=j+1

          write(lupri,*) ' State No.   Exc. Energy'
          write(lupri,*) ' ======================='
          write(lupri,'(I7,3X,g18.8)') MyRspFunc(j)%ExNumber1, ExEnergies(j)

          write(lupri,*)
          write(lupri,'(2X,A,17X,A,17X,A,17X,A)') 'Atom', 'x', 'y', 'z'

          counter=0
          do k=1,natoms
             write(lupri,'(3X,A,7X,3g18.8)') molcfg%setting%molecule(1)%p%ATOM(k)%Name,&
                  &ESG(j,k,1:3)
          end do

          write(lupri,*) 
          write(lupri,*) 
          write(lupri,*) 
          write(lupri,*) 


       end if StateWasCalculated

    end do PrintESG

    write(lupri,*) 
    write(lupri,*) 
    write(lupri,*) 
    write(lupri,*) 


    call mem_dealloc(doState)
    deallocate(MyRspFunc)
    call mem_dealloc(ExEnergies)
    call mem_dealloc(ESG)
    deallocate(rsp_results_vec)
    deallocate(GSG)

  end subroutine ESGresponse_driver






  !> \brief Driver for calculating excited state dipole moment (ESD).
  !> ESD is determined by subtracting a part of the double residue of the quadratic response function
  !> from the ground state dipole moment.
  !> See e.g. JCP 82, 3235 (1985), Eq. (3.21) with m=q.
  !> \author Kasper Kristensen                                                                
  !> \date 2010-08
  subroutine ESDresponse_driver(molcfg,F,D,S,ESDinput,DipoleMoment)
    implicit none
    !> Info on molecule needed by solver and integral programs
    type(rsp_molcfg), intent(inout)   :: molcfg
    !> Unperturbed Fock matrix
    type(matrix), intent(in)          :: F
    !> Unperturbed density matrix
    type(matrix), intent(in)          :: D
    !> Unperturbed overlap matrix
    type(matrix), intent(in)          :: S
    !> Contains ESD input
    type(ESDinputItem),intent(inout)   :: ESDinput
    !> Dipole moment
    real(realk) :: DipoleMoment(3)
    integer :: nexci, state, size_rsp
    integer, pointer :: lupri
    logical :: AllStates
    logical, allocatable :: doState(:)
    complex(realk), allocatable :: rsp_results_vec(:,:)
    real(realk), allocatable :: ExEnergies(:), ESD(:,:)
    type(RspFunc), allocatable :: MyRspFunc(:) 
    integer :: i,j,k


    ! File unit number for LSDALTON.OUT
    lupri => molcfg%lupri




    ! 1. Set number of excited states
    ! '''''''''''''''''''''''''''''''
    ! (Defined by .NEXCIT in input or set to 1 by default)

    nexci = molcfg%decomp%cfg_rsp_nexcit
    if(ESDinput%specific_states_in_input) then
       AllStates = .false.  ! only calculate ESD for specific states
    else 
       AllStates = .true.   ! calculate ESD for ALL excited states from 1 to nexci
    end if



    ! 2. Determine which particular excited states to calculate ESD for
    ! '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

    ! Logical for keeping track of which states to address
    allocate(doState(nexci))

    doState(:)=.false.
    if(AllStates) then ! ESD for all excited states from 1 to nexci
       size_rsp=nexci
       doState(:) = .true.
    else ! Only states specified by input
       size_rsp=ESDinput%esd_nexci
       do i=1,ESDinput%esd_nexci
          state=ESDinput%ExStates(i)
          doState(state) = .true.
       end do
    end if


    allocate(MyRspFunc(size_rsp))
    allocate(ExEnergies(size_rsp))
    ! 3 dipole coordinates (xyz) for each excitation energy
    allocate(ESD(size_rsp,3))
    allocate(rsp_results_vec(size_rsp,3))




    ! 3. Setup response functions
    ! '''''''''''''''''''''''''''

    ! We assume now that the transition density matrices for state 1 to nexci
    ! has been calculated using
    ! the Calculate_and_store_transition_density_matrices subroutine.
    ! Therefore these matrices are stored in rsp_eq_sol(i)%D in rsp_equations.f90,
    ! while the corresponding excitation energies are
    ! stored in rsp_eq_sol(i)%fld(1)%freq as complex frequencies
    ! where the imaginary component is zero.
    ! (It is checked in configuration.f90 that the ESD input is consistent).


    ! The ESD may be formulated in terms of a double residue
    ! of a quadratic response function, 
    ! see e.g. JCP 82, 3235 (1985), Eq. (3.21) with mu_b=mu_c=B and m=q:
    !
    ! Double residue = - |<0|B|m>|^2 * { <m| mu_a |m> - <0| mu_a |0> }
    !
    ! Thus the excited state expectation value <m| mu_a |m> becomes:
    ! 
    ! <m| mu_a |m> = <0| mu_a |0> - { Double residue / |<0|B|m>|^2 } 
    !
    ! We do not care about the B operator since in a practical calculation
    ! we factorize the double residue of <<A;B,C>> as:
    ! 
    ! Double residue = <<A; EXCI,EXCI>> * <<B;EXCI>> * <<C;EXCI>>
    ! 
    ! where <<B; EXCI>> = <<C; EXCI>> = <0|B|m>
    ! 
    ! Thus, we may calculate the expectation value for the mu_a operator in state |m> as:
    !
    ! <m| mu_a |m> = <0| mu_a |0> - <<mu_a;EXCI,EXCI>>  
    !
    ! Now the tricky sign conventions. Our formalism is based on energy derivatives.
    ! The dipole is by definition MINUS the energy derivative with respect to an electric field.
    ! Therefore, each dipole property that we calculate by using the 'EL  ' operator
    ! will have opposite sign to the corresponding dipole property.
    ! In other words, we make the equivalence:
    !
    ! mu_a <--->  - EL_a
    !
    ! Appylying this to the expression above we see that:
    !
    ! <m| mu |m> = <0| mu |0> + <<EL;EXCI,EXCI>>  (*)
    !
    ! where <<EL;EXCI,EXCI>> is the generalized quadratic response function that 
    ! we get by calling the general quadratic response routine with operators (EL,EXCI,EXCI),
    ! and <0| mu_a |0> is the dipole moment for the ground state which is input to 
    ! the ESD subroutine.
    ! (Note that when the ground state dipole moment was calculated in
    !  subroutine Get_dipole_moment the sign was also changed after calling the
    !  general response framework with operator 'EL  '. Therefore, DipoleMoment(3)
    !  already has the correct sign).


    ! Use j as counter in case not all excited states are requested
    j = 0
    SetupESDResponseFunc: do i=1,nexci

       StateRequested: if(doState(i)) then
          j=j+1

          ! Initialize response function
          call initialize_RspFunc(MyRspFunc(j))

          ! Order of response function=3 (quadratic response function double residue)
          MyRspFunc(j)%order = 3

          ! Operators = (EL,EXCI,EXCI) yield ESD (see above)
          MyRspFunc(j)%code(1) = 'EL  '
          MyRspFunc(j)%code(2) = 'EXCI'
          MyRspFunc(j)%code(3) = 'EXCI'

          ! Calculate 3 (xyz) components of the ESD, whereas
          ! the number of components for the second and third (EXCI) operators is only 1.
          ! Note that EACH excitation energy requires its own response function.
          MyRspFunc(j)%first_comp(1) = 1
          MyRspFunc(j)%first_comp(2) = 1
          MyRspFunc(j)%first_comp(3) = 1
          MyRspFunc(j)%dims(1) = 3
          MyRspFunc(j)%dims(2) = 1
          MyRspFunc(j)%dims(3) = 1

          ! We are interested in a double residue (transition moment).
          MyRspFunc(j)%trans_moment = .true.

          ! The transition moment corresponds to a DOUBLE residue,
          ! therefore the order of the transition moment is 2.
          MyRspFunc(j)%trans_moment_order = 2


          ! If all excited states are requested MyRspFunc(j) simply refers to
          ! excited state number j, whereas if specific excited states
          ! have been defined in the input, then MyRspFunc(j) refers to the
          ! excited state number stored in ESDinput%ExStates(j).
          ! Both excited states in the double residue are the same 
          ! since q=m in JCP 82, 3235 (1985), Eq. (3.21): ExNumber1=ExNumber2
          if(AllStates) then
             MyRspFunc(j)%ExNumber1 = j
             MyRspFunc(j)%ExNumber2 = j
          else
             MyRspFunc(j)%ExNumber1 = ESDinput%ExStates(j)
             MyRspFunc(j)%ExNumber2 = ESDinput%ExStates(j)
          end if

          ! Accoding to JCP 82, 3235 (1985), Eq. (3.21) the residue of interest has
          ! Frequency 2 = -excitation energy
          ! Frequency 3 = excitation energy
          ! Therefore Frequency 1 = 0

          ! The excitation energy
          ! is stored in rsp_eq_sol(X)%fld(1)%freq where
          ! X refers to the position of the excited state now
          ! stored in MyRspFunc(j)%ExNumber1.
          ! Note that the frequencies stored in rsp_eq_sol(X)%fld(1)%freq are complex
          ! frequencies where the imaginary component is zero.
          MyRspFunc(j)%freq(3) =  rsp_eq_sol(MyRspFunc(j)%ExNumber1)%fld(1)%freq
          MyRspFunc(j)%freq(2) = -rsp_eq_sol(MyRspFunc(j)%ExNumber1)%fld(1)%freq
          MyRspFunc(j)%freq(1) = 0E0_realk
          ! For later convience we also have a vector of REAL excitation energies
          ExEnergies(j) = real(MyRspFunc(j)%freq(3))

          ! Zero response function
          rsp_results_vec(j,:)=(0E0_realk,0E0_realk)

       end if StateRequested

    end do SetupESDResponseFunc



    ! 4. Calculate excited state dipole moment
    ! ''''''''''''''''''''''''''''''''''''''''


    j=0
    CalculateTransMoment: do i=1,nexci

       CalculateState: if(doState(i)) then
          j=j+1

          ! Calculate <<EL;EXCI,EXCI>> response functions for all 3 (xyz) components of EL.
          ! The quadratic_response subroutine now reads the transition density matrices from
          ! rsp_eq_sol(j)%D in rsp_equations.f90.
          ! Note: We use the n+1 rule here!
          call quadratic_response_nplus1(molcfg,F,D,S,MyRspFunc(j), &
               & rsp_results_vec(j,1:3),3)

          ! The three components of <<EL;EXCI,EXCI>> are now stored in rsp_results_vec as
          ! complex numbers with zero imaginary part.
          ! For each excitation energy j, we calculate the ESD using Eq. (*) above.
          ! The ESD array is a REAL array and has the structure:
          ! (excitation, xyz-component)
          do k=1,3
             ESD(j,k) = DipoleMoment(k) + real(rsp_results_vec(j,k))
          end do

       end if CalculateState

    end do CalculateTransMoment




    ! 5. Print ESD and excitation energies
    ! ''''''''''''''''''''''''''''''''''''

    write(lupri,*) 
    write(lupri,*) 
    write(lupri,*) 
    write(lupri,*) 
    write(lupri,'(2X,A)') '*************************************************************&
         &**********************'
    write(lupri,'(2X,A)') '*                      EXCITED STATE DIPOLE RESULTS (in a.u.)&
         &                     *'
    write(lupri,'(2X,A)') '*************************************************************&
         &**********************'
    write(lupri,*)
    write(lupri,*) 
    write(lupri,*) 

    write(lupri,'(1X,A,12X,A,17X,A,17X,A)') ' State No.   Exc. Energy', 'x','y','z'
    write(lupri,*) ' ========================================================&
         &=========================='

    write(lupri,'(7X,A,7X,3g18.8)') '(Ground state)', DipoleMoment(1:3)


    j=0
    PrintESD: do i=1,nexci

       StateWasCalculated: if(doState(i)) then
          j=j+1

          write(lupri,'(I7,3X,4g18.8)') MyRspFunc(j)%ExNumber1, ExEnergies(j), ESD(j,1:3)

       end if StateWasCalculated

    end do PrintESD

    write(lupri,*) 
    write(lupri,*) 
    write(lupri,*) 
    write(lupri,*) 


    deallocate(doState)
    deallocate(MyRspFunc)
    deallocate(ExEnergies)
    deallocate(ESD)
    deallocate(rsp_results_vec)
    if(esdinput%specific_states_in_input) deallocate(esdinput%ExStates)

  end subroutine ESDresponse_driver

 subroutine set_Bterm_RspFunc(MyRspFunc,i,mC,eC,freq,MAGNETIC)
   implicit none
   !> response function structure to be built
   type(RspFunc) :: MyRspFunc
   !> exited state 
   integer,intent(in) :: i
   !> magnetic component
   integer,intent(in) :: mC
   !> electric component
   integer,intent(in) :: eC
   !> the excitation freq
   complex(realk),intent(in) :: freq
   !> string of MAG or MAGO
   character(len=4),intent(in) :: MAGNETIC

   call initialize_RspFunc(MyRspFunc)
   MyRspFunc%order = 3     
   MyRspFunc%code(1) = MAGNETIC
   MyRspFunc%code(2) = 'EL  '
   MyRspFunc%code(3) = 'EXCI'     
   MyRspFunc%first_comp(1) = mC
   MyRspFunc%first_comp(2) = eC  
   MyRspFunc%dims(1) =  1    
   MyRspFunc%dims(2) =  1    
   MyRspFunc%freq(1) =  0E0_realk    
   MyRspFunc%freq(2) = -freq  
   MyRspFunc%freq(3) =  freq        
   MyRspFunc%trans_moment = .TRUE. !single residue
   MyRspFunc%trans_moment_order = 1
   MyRspFunc%ExNumber1 = i
 end subroutine set_Bterm_RspFunc

 Subroutine set_Aterm_RspFunc(MyRspFunc,i1,i2,allowed,freq,MAGNETIC)
   implicit none
   !> response function structure to be built
   type(RspFunc) :: MyRspFunc
   !> exited state 
   integer,intent(in) :: i1,i2
   !> dipole allowed or not
   logical,intent(in) :: allowed(:,:)    
   !> the excitation freq
   complex(realk),intent(in) :: freq
   !> string of MAG or MAGO
   character(len=4),intent(in) :: MAGNETIC
   !> local parameters
   logical :: MAG_ALLOWED(3),allow,yforbidden               
   integer :: dim,j,start

   call initialize_RspFunc(MyRspFunc)
   MAG_ALLOWED = .FALSE.
   IF((allowed(1,i1).AND.allowed(2,i2)).OR.(allowed(1,i2).AND.allowed(2,i1)))&
        & MAG_allowed(3) = .TRUE.
   IF((allowed(3,i1).AND.allowed(1,i2)).OR.(allowed(3,i2).AND.allowed(1,i1)))&
        MAG_allowed(2) = .TRUE.
   IF((allowed(2,i1).AND.allowed(3,i2)).OR.(allowed(2,i2).AND.allowed(3,i1)))&
        MAG_allowed(1) = .TRUE.
   MyRspFunc%order = 3     
   MyRspFunc%code(1) = MAGNETIC
   MyRspFunc%code(2) = 'EXCI'     
   MyRspFunc%code(3) = 'EXCI'     
   dim=0
   yforbidden=.false.
   do j=3,1,-1
      if(MAG_allowed(j))then
         start = j 
         dim = dim+1
         if((j.eq. 1).and.yforbidden)dim = dim+1 !x and z allowed but y forbidden, we need to include y
      endif
      if((j.eq. 2).and.(dim.eq. 1))yforbidden=.true. !z allowed but y forbidden
   enddo
   MyRspFunc%first_comp(1) = start    
   MyRspFunc%dims(1) = dim    
   MyRspFunc%freq(1) = 0E0_realk    
   MyRspFunc%freq(2) = -freq
   MyRspFunc%freq(3) = freq   
   MyRspFunc%trans_moment = .TRUE. !double residue
   MyRspFunc%trans_moment_order = 2 
   MyRspFunc%ExNumber1 = i1
   MyRspFunc%ExNumber2 = i2
 end Subroutine set_Aterm_RspFunc

 subroutine set_verdet_RspFunc(MyRspFunc,mC,eC1,eC2,freq,gamma,MAGNETIC)
   implicit none
   !> response function structure to be built
   type(RspFunc) :: MyRspFunc
   !> magnetic component
   integer,intent(in) :: mC
   !> electric component
   integer,intent(in) :: eC1,eC2
   !> the excitation freq
   real(realk),intent(in) :: freq,gamma
   !> string of MAG or MAGO
   character(len=4),intent(in) :: MAGNETIC

   call initialize_RspFunc(MyRspFunc)
   MyRspFunc%order = 3     
   MyRspFunc%code(1) = MAGNETIC
   MyRspFunc%code(2) = 'EL  '
   MyRspFunc%code(3) = 'EL  '     
   MyRspFunc%first_comp(1) = mC
   MyRspFunc%first_comp(2) = eC1  
   MyRspFunc%first_comp(3) = eC2  
   MyRspFunc%dims(1) =  1    
   MyRspFunc%dims(2) =  1    
   MyRspFunc%dims(3) =  1    
   MyRspFunc%freq(1) =  0E0_realk    
!   cfg_rsp_gamma = gamma           
   MyRspFunc%freq(2) = -freq*(1E0_realk,0E0_realk)-gamma*(0E0_realk,1E0_realk)  
   MyRspFunc%freq(3) =  freq*(1E0_realk,0E0_realk)+gamma*(0E0_realk,1E0_realk)  
   MyRspFunc%trans_moment = .FALSE. 
   MyRspFunc%trans_moment_order = 1
 end subroutine set_verdet_RspFunc

 subroutine set_doubleResidue_RspFunc(MyRspFunc,i,j,eC,freqI,freqJ,OperatorString)
   implicit none
   !> response function structure to be built
   type(RspFunc) :: MyRspFunc
   !> exited state i
   integer,intent(in) :: i
   !> exited state j
   integer,intent(in) :: j
   !> operator component (for instance electric component)
   integer,intent(in) :: eC
   !> the excitation freq for excitated state nr i
   complex(realk),intent(in) :: freqI
   !> the excitation freq for excitated state nr j
   complex(realk),intent(in) :: freqJ
   !> operator string like ( MAG or MAGO or EL )
   character(len=4),intent(in) :: OPERATORSTRING
!   !> logical which determine if the double residue should be app. 
!   logical,intent(in) :: TrueDoubleResidue
   call initialize_RspFunc(MyRspFunc)
   MyRspFunc%order = 3     
   MyRspFunc%code(1) = OPERATORSTRING
   MyRspFunc%code(2) = 'EXCI'
   MyRspFunc%code(3) = 'EXCI'
   MyRspFunc%first_comp(1) = eC
   MyRspFunc%first_comp(2) = 1  
   MyRspFunc%first_comp(3) = 1  
   MyRspFunc%dims(1) =  1
   MyRspFunc%dims(2) =  1    
   MyRspFunc%dims(3) =  1    
   MyRspFunc%freq(1) =  0E0_realk    
   MyRspFunc%freq(2) = -freqI  
   MyRspFunc%freq(3) =  freqJ       
   MyRspFunc%trans_moment = .TRUE.   !double residue
   MyRspFunc%trans_moment_order = 2  !double residue
   MyRspFunc%ExNumber1 = i
   MyRspFunc%ExNumber2 = j
!   MyRspFunc%TrueDoubleResidue = TrueDoubleResidue
 end subroutine set_doubleResidue_RspFunc

 !> \brief main MCD driver
 !> \author Thomas Kjaergaard
 !> \date 2010-03
 !>
 !> Main MCD A and B terms driver. Sets up response input, prints input, 
 !> calculates response function,
 !> and prints out the response results.
 !> Calculates A and B terms of Magnetic Circular Dichroisme (MCD)
 !> based on J. Chem. Theory Comput. 2009, 5, 1997-2020  but reformulated using JCP, 129, 214108 (2008).
 !> <br><br>
 !> It should work for any number of excitations requested but at the moments it calculates 
 !> all up to the excitation specified in input 
 !> perturbation dependent basis sets (e.g. London orbitals) as well as nonlondon basis sets are default
 !> but purely london or purely non london calculations can be performed with the correct input
 !> can be done and they can be combined.
 !> <br><br>
 !> A damped MCD spectrum can also be obtained and as default 
 !> 1. the standard A and B terms are constructed
 !> 2. A simulated spectrum is created from the calculated A and B terms
 !> 3. A damped spectrum is calculated using points a round the peaks determined in 1.
 !> <br><br>
 !> The input structure for response functions is as follows:
 !> <br><br>
 !> **RESPONS
 !> <br>
 !> *QUASIMCD    
 !> <br>
 !> .MCDEXCIT
 !> <br>
 !> (number) (highest excitation to be examined)
 !> <br>
 !> .NO LONDON (optional, turns off the use of london orbitals - pure basis sets are used)
 !> <br>v
 !> .NO NONLONDON (optional, turns off the use of nonlondon orbitals - purely london orbitals are used)
 !> <br>
 !> .NO SIMULATE (optional, turns of the creation of a simulated spectra from calculated A and B terms, and a gnuplot script to plot it)
 !> <br>
 !> .GAUSSIAN (optional, request that the simulated is buildt using gaussians lineshape functions, lorentz type lineshape functions are default)
 !> <br>
 !> .LINESHAPEPARAM (optional, gives the lineshape parameter used to plot the simulated spectra)
 !> <br>
 !> (lineshape parameter) (default is 0.005E0_realk for lorentz, 0.0070851079363103793E0_realk for gaussian)
 !> <br>
 !> .NSTEPS (optional, gives points used in the simulated spectrum)
 !> <br>
 !> (number of points) (default is 5000)
 !> <br>
 !> .NVECFORPEAK (optional, gives how many damped point calculations there is scheduled for each A and B terms)
 !> <br>
 !> (number of points) (default is 10)
 !> <br>
 !> .NO DAMPED (optional, turns of the creation of a damped spectra, and a gnuplot script to plot it)
 !> <br>
 !> .NO ATERM (optional, turns of the calculation of A terms)
 !> <br>
 !> .NO BTERM (optional, turns of the calculation of B terms)
 !> <br>
 !> .DAMPEDXCOOR (optional, schedule points for a damped calculation)
 !> <br>
 !> (number of points) (for instance 3)
 !> <br>
 !> (energy point 1)  (point in units of au. )
 !> <br>
 !> (energy point 2)  (point in units of au. )
 !> <br>
 !> (energy point 3)  (point in units of au. )
 !> <br>
 !> *END QUASIMCD
 !> <br>
 !> $END RESPONS
 !> <br> <br>
 !> Questions regarding this driver, the input structure etc. may be addressed to tkjaergaard@chem.au.dk.
 subroutine MCDresponse_driver(lupri,setting,decomp,solver,F,D,S,MCDinput)
  implicit none
  !> logical unit number of output 
  integer,intent(in)                :: lupri
  !> Info on molecule needed by solver and integral programs
  type(LSSETTING),intent(in),target :: setting
  !> Contains matrices from OAO decomposition of overlap matrix
  type(decompItem),intent(inout),target :: decomp
  !> RSPsolver settings
  type(RSPSOLVERinputitem),intent(inout),target  :: solver
  !> Unperturbed Fock matrix
  type(matrix), intent(in)          :: F
  !> Unperturbed density matrix
  type(matrix), intent(in)          :: D
  !> Unperturbed overlap matrix
  type(matrix), intent(in)          :: S
  !> Contains MCD input
  type(MCDinputItem),intent(in)      :: MCDinput
  !> Total number of response functions and transition moments requested
  integer                           :: n_rspfunc
  !> The position of the highest lying excited state requested
  integer                           :: nexci_max
  !> Stack of all response function and transition moments to be determined
  type(RspFunc), allocatable        :: RspFunc_stack(:),RspFunc_stack2(:)
  !> temporary response function structure
  type(RspFunc)                     :: tmpRspFunc
  !> local parameters
  integer                           :: i,j,k,l,maxdim,start,MCDnexci,dim
  integer                           :: nAterms,oldAterm,oldbterm
  complex(realk), allocatable       :: rsp_results(:,:)
  type(rsp_molcfg) :: molcfg
  integer, target  :: natoms
  !    complex(realk), allocatable, target, save           :: rsp_results(:,:)
  complex(realk) :: freqold,freq
  integer,pointer :: MCDexci1(:),MCDexci2(:)
  complex(realk),pointer :: MCDfreq(:)
  logical :: allow
  logical,pointer :: allowed(:,:)    
  Real(realk),parameter :: PI=3.14159265358979323846E0_realk
  !the Bterm is in standard units given in D**2*muB/(cm^-1) (au : e**3*a0**4/hbar)
  ! D is the unit debye        :   1 D = 0.393430307 (au.) (au : e*a0)
  ! muB is the Bohr magneton   : 1 muB = 0.5 (au.)
  ! and inverse centimeter     : 1cm^⁻1= 219475 (au.)
  ! so the standard units are : D**2*muB/(cm^-1) = 16985.974069713367 (au.)
  ! so to convert from au to standard we inverse this to obtain
  Real(realk),parameter :: BtermConstant=5.88764*1E-05_realk
  !the Aterm is in standard units given in D**2*muB (au : e**3*hbar*a0**2/me)
  ! so the standard units are : D**2*muB = 0.07739370323305711 (au.)
  ! so to convert from au to standard we inverse this to obtain
  Real(realk),parameter :: AtermConstant=12.920947806163007
  Real(realk),parameter :: VerdetConstant=5.88764*1E-05_realk*33.53!0.152123007*1.0E-7_realk
  !hc = 1239.8419300923944 eV*nm  
  !so  lambda (in nm) = hc/E (E in eV)
  !and E (in eV) = hc/lambda (lambda in nm)
  Real(realk),parameter :: hcConstant_eVnm=1239.8419300923944
  character(len=4) :: STRING4
  real(realk) :: eigenvalue(3),eigenvalueL(3),sign,signL,X1,Y1,Z1,X2,Y2,Z2
  real(realk),allocatable :: twophoton(:,:,:),Atermmoment(:,:),MCDBterm(:)
  real(realk),allocatable :: twophotonL(:,:,:),AtermmomentL(:,:),MCDBtermL(:)
  real(realk),allocatable :: Bmoment(:,:),Amoment(:,:,:),MCDAterm(:),MCDAtermL(:),freqA(:),freqB(:)
  integer,allocatable :: BtermRspFunc(:),AtermRspFunc(:)
  integer,allocatable :: BtermRspFuncLondon(:),AtermRspFuncLondon(:)
  !
  INTEGER :: nBterms,LUG,Nstep,lengthX,maxA,maxB,lumcd,lumcd2,nXcoor2,nXcoor3,lumcd3,lumcd4
  real(realk),pointer  :: Y(:),YL(:),GAMMAATERM(:),GAMMABTERM(:)
  real(realk)    :: Range,step,Estart,Eend,FUN,INT,x,DISPLACEMENT,FA
  real(realk) :: gamma,DTHR
  real(realk),allocatable :: PUREB(:),PUREBL(:),PUREA(:),PUREAL(:)
  real(realk),allocatable :: Xcoor(:),COMBINED(:),COMBINEDL(:),Deltavec(:),Xcoor2(:),Xcoor3(:)
  complex(realk),allocatable :: verdet(:),verdetL(:) 
  character(len=30)  :: INPUTFILENAME,OUTPUTFILENAME
  character(len=7)   :: TYPE
  logical            :: input_exists,LORENTZ,unique,degenerateStates
  logical            :: FOUNDINUSERLIST
  real(realk),parameter :: TWOPIM1 = 0.15915494309189535E0_realk
  real(realk),parameter :: hartree=27.21138386E0_realk
  Real(realk), parameter ::eVTOcm= 8065.54E0_realk !eV to cm^-1
  real(realk), parameter :: autocm=hartree*eVTOcm 
  write(lupri,*)'******************************************************'
  write(lupri,*)'**                                                  **'  
  write(lupri,*)'**  MCD Calculation                                 **'
  write(lupri,*)'**                                                  **'  
  write(lupri,*)'******************************************************'
  
  nexci_max = MCDinput%nexci
  degenerateStates = MCDinput%degenerateStates !same as solver%degenerateState
  DTHR = solver%degenerateTHR
  !create config struct to be passed to rsp_contribs / rsp_equations
!/*point to natoms within structure*/
  molcfg = rsp_molcfg(S,setting%MOLECULE(1)%p%Natoms, &
  & decomp%lupri,decomp%luerr,setting,decomp,solver)
  molcfg%zeromat = 0*S
  LORENTZ = .TRUE.
  IF(MCDinput%useinputgamma)then
     GAMMA = MCDinput%Gamma
  else!default
     if(MCDinput%lorentz)then
        GAMMA=0.005E0_realk
     else
        GAMMA=0.0070851079363103793E0_realk
     endif
  endif

  if(MCDinput%nexci.GT. 0)then
   !calc excitation energies and transition densities
   call transition_moment_density_matrix(molcfg,F,D,S,nexci_max)
   call write_transition_density_matrix(nexci_max,decomp%lupri)

   write(lupri,*)'MCD: done transition_moment_density_matrix'
   maxdim = 3
   n_rspfunc = nexci_max 
   allocate(rsp_results(maxdim,n_rspfunc))
   rsp_results(:,:) = (0.0E0_realk, 0.0E0_realk)
   allocate(allowed(3,nexci_max ))

   WRITE(lupri,'(A)')'MCD: TRANSITION MOMENTS'
   WRITE(lupri,'(2X,A5,2X,A8,2X,A8,9X,A1,11X,A1,11X,A1,7X,A8)')&
        &'STATE','Freq(au)','Freq(eV)','X','Y','Z','ALLOWED?'
   !calc transition moments to determine dipole allowed transitions
   do i = 1,nexci_max
      call initialize_RspFunc(tmpRspFunc)
      tmpRspFunc%order = 2     
      tmpRspFunc%code(1) = 'EL  '
      tmpRspFunc%code(2) = 'EXCI'     
      tmpRspFunc%first_comp(1) = 1    
      tmpRspFunc%dims(1) = 3    
      tmpRspFunc%freq(2) = -rsp_eq_sol(i)%fld(1)%freq
      tmpRspFunc%freq(3) =  rsp_eq_sol(i)%fld(1)%freq
      tmpRspFunc%trans_moment = .TRUE. !single residue
      tmpRspFunc%trans_moment_order = 1 
      tmpRspFunc%ExNumber1 = i
      !calc transition moments 
      call linear_response(molcfg,F,D,S,tmpRspFunc, &
           &rsp_results(1:maxdim,i), maxdim)
      do j=1,3 !electricfield components
         IF(ABS(rsp_results(j,i)) .GT. 1.0E-5_realk)THEN
            allowed(j,i) = .TRUE.
         ELSE
            allowed(j,i) = .FALSE.
         ENDIF
      enddo
      STRING4 = ' NO ' 
      IF(ABS(rsp_results(1,i)*rsp_results(1,i)+rsp_results(2,i)*rsp_results(2,i)+&
           &rsp_results(3,i)*rsp_results(3,i)) .GT. 1.0E-5_realk)THEN
         STRING4 = ' YES' 
      ENDIF
      WRITE(lupri,'(I5,2X,F10.6,2X,F10.6,2X,F10.7,2X,F10.7,2X,F10.7,6X,A4)') &
           &i,REAL(rsp_eq_sol(i)%fld(1)%freq),REAL(rsp_eq_sol(i)%fld(1)%freq*hartree),&
           &REAL(rsp_results(1,i)),REAL(rsp_results(2,i)),REAL(rsp_results(3,i)),STRING4
   enddo

   !remove dipole forbidden transitions and find degenerate states
   MCDnexci = 0
   freqold=(0E0_realk,0E0_realk)
   allocate(MCDexci1(nexci_max))
   allocate(MCDexci2(nexci_max))
   allocate(MCDfreq(nexci_max))
   nAterms = 1
   allocate(Amoment(nexci_max,2,3)) !transitionmoments saved for Aterm eval.
   dipoleloop: do i = 1,nexci_max
      freq=rsp_eq_sol(i)%fld(1)%freq
      allow = .FALSE.
      do j=1,3 
         if(allowed(j,i))allow=.TRUE.
      enddo
      if(allow)then
         IF(MCDinput%specific_states_in_input)THEN
            do j=1,MCDinput%nMCDexci
               FOUNDINUSERLIST = .FALSE.
               IF(MCDinput%ExStates(j).EQ.i)THEN
                  FOUNDINUSERLIST = .TRUE.
               ENDIF
            enddo
            IF(.NOT.FOUNDINUSERLIST) CYCLE dipoleloop
         ENDIF
         if(degenerateStates.AND.abs(freq-freqold)<DTHR) then
            !degenerate
            MCDexci2(MCDnexci) = MCDexci2(MCDnexci)+i
            Amoment(nAterms,2,1) = REAL(rsp_results(1,i))
            Amoment(nAterms,2,2) = REAL(rsp_results(2,i))
            Amoment(nAterms,2,3) = REAL(rsp_results(3,i))
            nAterms = nAterms+1
            if(MCDexci2(MCDnexci).GT.i)then
               !degeneracy higher than 2
               call lsquit('MCD: degeneracy higher than 2',lupri)
            endif
            MCDnexci = MCDnexci + 1
            MCDexci1(MCDnexci) = i
            MCDexci2(MCDnexci) = 0
            MCDfreq(MCDnexci) = rsp_eq_sol(i)%fld(1)%freq
         else
            Amoment(nAterms,1,1) = REAL(rsp_results(1,i))
            Amoment(nAterms,1,2) = REAL(rsp_results(2,i))
            Amoment(nAterms,1,3) = REAL(rsp_results(3,i))
            MCDnexci = MCDnexci + 1
            MCDexci1(MCDnexci) = i
            MCDexci2(MCDnexci) = 0
            MCDfreq(MCDnexci) = rsp_eq_sol(i)%fld(1)%freq
         endif
      endif
      freqold=freq
   enddo dipoleloop
   write(lupri,'(A)')' '
   write(lupri,'(A)')'Excited states used for MCD calculation'
   do i = 1,MCDnexci
      if(MCDexci2(i).EQ. 0)then
         write(lupri,'(A,F12.8)')  'Excitation Energies(au):',REAL(rsp_eq_sol(MCDexci1(i))%fld(1)%freq)
      else
         write(lupri,'(A,F12.8,A)')'Excitation Energies(au):',REAL(rsp_eq_sol(MCDexci1(i))%fld(1)%freq),'Degenerate'
      endif
   enddo
  else
     IF(MCDinput%simulate)WRITE(lupri,'(A)')&
          &'WARNING: In order to do a simulated MCD spectrum you must &
          &specify the excited states of intrerest using the .MCDEXCIT keyword'
     MCDnexci=MCDinput%nexci
  endif

   IF(MCDnexci .GT. 0)THEN
      n_rspfunc = 9*nexci_max 
      allocate(BtermRspFuncLondon(n_rspfunc))
      allocate(AtermRspFuncLondon(n_rspfunc))
      allocate(BtermRspFunc(n_rspfunc))
      allocate(AtermRspFunc(n_rspfunc))
      if(MCDinput%london.AND.MCDinput%nolondon)n_rspfunc=n_rspfunc*2
      allocate(RspFunc_stack(n_rspfunc))
      BtermRspFuncLondon = 0 
      AtermRspFuncLondon = 0
      BtermRspFunc = 0 
      AtermRspFunc = 0
      !****************************************************************************
      !* Set up Aterm response input to prepare double resdiue <i|m|j> 
      !* where i,j are dipole allowed degenerate states
      !****************************************************************************
      k=0
      IF(degenerateStates.AND.MCDinput%doAterms)THEN
         if(MCDnexci.EQ. 0)call lsquit('somethings wrong 0 bterms',lupri)
         allocate(freqA(MCDnexci))
         nAterms = 0
         do i = 1,MCDnexci
            if(MCDexci2(i).ne. 0)then !degenerate
               nAterms = nAterms+1
               freqA(nAterms) = MCDfreq(i)   
               if(MCDinput%london)then
                  k=k+1
                  call set_Aterm_RspFunc(RspFunc_stack(k),MCDexci1(i),MCDexci2(i),allowed,MCDfreq(i),'MAG ')
                  AtermRspFuncLondon(k) = nAterms
               endif
               if(MCDinput%nolondon)then
                  k=k+1
                  call set_Aterm_RspFunc(RspFunc_stack(k),MCDexci1(i),MCDexci2(i),allowed,MCDfreq(i),'MAGO')
                  AtermRspFunc(k) = nAterms
               endif
            endif
         enddo
         IF(nAterms.EQ.0)WRITE(lupri,*)'No Aterms found'
      ELSE
         WRITE(lupri,*)'No Aterms is calculated'
         nAterms = 0
      ENDIF
      IF(.NOT.degenerateStates)THEN
         WRITE(lupri,*)'Note that a calculation of Aterms requires degenerate states'
         WRITE(lupri,*)'and therefore the keyword .DEGENERATE'
      ENDIF
      !****************************************************************************
      !* Set up Bterm response input to prepare single resdiue (<i|m|k><k|mu|j>/w..+..)
      !* also called two photon transition moment, where i,j are dipole allowed states
      !****************************************************************************
      IF(MCDinput%doBterms)THEN
         allocate(Bmoment(MCDnexci,3)) !transitionmoments saved for Bterm eval
         allocate(freqB(MCDnexci)) 
         do i = 1,MCDnexci
            do j=1,3
               Bmoment(i,j) = REAL(rsp_results(j,MCDexci1(i)))
            enddo
            freqB(i) = MCDfreq(i)   
            if(MCDinput%london)then
               if(allowed(1,MCDexci1(i)))then   !X dipole allowed - yz and zy two photon transition moment
                  k = k+1 ; BtermRspFuncLondon(k) = i
                  call set_Bterm_RspFunc(RspFunc_stack(k),MCDexci1(i),2,3,MCDfreq(i),'MAG ') ! yz
                  k = k+1; BtermRspFuncLondon(k) = i
                  call set_Bterm_RspFunc(RspFunc_stack(k),MCDexci1(i),3,2,MCDfreq(i),'MAG ') ! zy
               endif
               if(allowed(2,MCDexci1(i)))then   !Y dipole allowed - zx and xz two photon transition moment
                  k = k+1; BtermRspFuncLondon(k) = i
                  call set_Bterm_RspFunc(RspFunc_stack(k),MCDexci1(i),3,1,MCDfreq(i),'MAG ') !zx
                  k = k+1; BtermRspFuncLondon(k) = i
                  call set_Bterm_RspFunc(RspFunc_stack(k),MCDexci1(i),1,3,MCDfreq(i),'MAG ') !xz
               endif
               if(allowed(3,MCDexci1(i)))then   !Z dipole allowed - xy and yx two photon transition moment
                  k = k+1; BtermRspFuncLondon(k) = i
                  call set_Bterm_RspFunc(RspFunc_stack(k),MCDexci1(i),1,2,MCDfreq(i),'MAG ') !xy
                  k = k+1; BtermRspFuncLondon(k) = i
                  call set_Bterm_RspFunc(RspFunc_stack(k),MCDexci1(i),2,1,MCDfreq(i),'MAG ') !yx
               endif
            endif
            if(MCDinput%nolondon)then
               if(allowed(1,MCDexci1(i)))then   !X dipole allowed - yz and zy two photon transition moment
                  k = k+1; BtermRspFunc(k) = i
                  call set_Bterm_RspFunc(RspFunc_stack(k),MCDexci1(i),2,3,MCDfreq(i),'MAGO') ! yz
                  k = k+1; BtermRspFunc(k) = i
                  call set_Bterm_RspFunc(RspFunc_stack(k),MCDexci1(i),3,2,MCDfreq(i),'MAGO') ! zy
               endif
               if(allowed(2,MCDexci1(i)))then   !Y dipole allowed - zx and xz two photon transition moment
                  k = k+1; BtermRspFunc(k) = i
                  call set_Bterm_RspFunc(RspFunc_stack(k),MCDexci1(i),3,1,MCDfreq(i),'MAGO') !zx
                  k = k+1; BtermRspFunc(k) = i
                  call set_Bterm_RspFunc(RspFunc_stack(k),MCDexci1(i),1,3,MCDfreq(i),'MAGO') !xz
               endif
               if(allowed(3,MCDexci1(i)))then   !Z dipole allowed - xy and yx two photon transition moment
                  k = k+1; BtermRspFunc(k) = i
                  call set_Bterm_RspFunc(RspFunc_stack(k),MCDexci1(i),1,2,MCDfreq(i),'MAGO') !xy
                  k = k+1; BtermRspFunc(k) = i
                  call set_Bterm_RspFunc(RspFunc_stack(k),MCDexci1(i),2,1,MCDfreq(i),'MAGO') !yx
               endif
            endif
         enddo
      ENDIF
      IF(MCDinput%doBterms.OR.MCDinput%doAterms.AND.degenerateStates)THEN
         n_RspFunc = k 
         !    Print input for response calculation, this print should ease identification of input errors.
         call print_response_func(RspFunc_stack, n_RspFunc, .false.,lupri)

         !********************************************************************************
         !*              Calculate all response functions.
         !********************************************************************************
         IF(n_RspFunc .GT. 0)THEN
            call rspfunc_or_transmoment(molcfg,F,D,S,RspFunc_stack,n_rspfunc)
            call print_response_func(RspFunc_stack, n_RspFunc, .true.,lupri) 
         ENDIF
         !********************************************************************************
         !*              Print Results
         !********************************************************************************
         IF(n_RspFunc .GT. 0)THEN
            allocate(twophotonL(n_RspFunc,3,3))
            twophotonL=0E0_realk
            allocate(AtermmomentL(n_RspFunc,3))
            AtermmomentL=0E0_realk
            allocate(twophoton(n_RspFunc,3,3))
            twophoton=0E0_realk
            allocate(Atermmoment(n_RspFunc,3))
            Atermmoment=0E0_realk
            do i = 1,n_RspFunc
               if(BtermRspFuncLondon(i) .ne. 0)then     !London bterm
                  twophotonL(BtermRspFuncLondon(i),RspFunc_stack(i)%first_comp(2),RspFunc_stack(i)%first_comp(1))&
                       & = REAL(RspFunc_stack(i)%result(1))
               elseif(BtermRspFunc(i) .ne. 0)then       !bterm
                  twophoton(BtermRspFunc(i),RspFunc_stack(i)%first_comp(2),RspFunc_stack(i)%first_comp(1))&
                       & = REAL(RspFunc_stack(i)%result(1))
               elseif(AtermRspFuncLondon(i) .ne. 0)then !London Aterm
                  k=0
                  do j=RspFunc_stack(i)%first_comp(1),RspFunc_stack(i)%first_comp(1)+RspFunc_stack(i)%dims(1)-1
                     k=k+1
                     AtermmomentL(AtermRspFuncLondon(i),j) = REAL(RspFunc_stack(i)%result(k))
                  enddo
               elseif(AtermRspFunc(i) .ne. 0)then       !Aterm
                  k=0
                  do j=RspFunc_stack(i)%first_comp(1),RspFunc_stack(i)%first_comp(1)+RspFunc_stack(i)%dims(1)-1
                     k=k+1
                     Atermmoment(AtermRspFunc(i),j) = REAL(RspFunc_stack(i)%result(k))
                  enddo
               else
                  call lsquit('neither A or B term?????????? Thomas K',lupri)
               endif
            enddo
         ENDIF
         
         allocate(MCDBterm(MCDnexci))
         allocate(MCDBtermL(MCDnexci))
         do i=1,MCDnexci    
            write(lupri,'(A)')' '
            write(lupri,'(A,F12.8)')'Excitation Energy (a.u.)         :', freqB(i)
            write(lupri,'(A,F12.8)')'Excitation Energy (eV)           :', freqB(i)*hartree
            write(lupri,'(A,F12.8)')'Excitation Energy (nm)           :', hcConstant_eVnm/(freqB(i)*hartree)
            write(lupri,'(A,F16.4)')'Excitation Energy (cm^-1)        :', freqB(i)*hartree*eVTOcm
            write(lupri,'(A)')' '
            write(lupri,'(A,F16.8)')'SR-contribution from XY L:    ',twophotonL(i,1,2)
            write(lupri,'(A,F16.8)')'SR-contribution from XZ L:    ',twophotonL(i,1,3)
            write(lupri,'(A,F16.8)')'SR-contribution from YX L:    ',twophotonL(i,2,1)
            write(lupri,'(A,F16.8)')'SR-contribution from YZ L:    ',twophotonL(i,2,3)
            write(lupri,'(A,F16.8)')'SR-contribution from ZX L:    ',twophotonL(i,3,1)
            write(lupri,'(A,F16.8)')'SR-contribution from ZY L:    ',twophotonL(i,3,2)
            write(lupri,'(A,F16.8)')'BTERM-contribution from X*YZ L:',Bmoment(i,1)*twophotonL(i,2,3)
            write(lupri,'(A,F16.8)')'BTERM-contribution from X*ZY L:',Bmoment(i,1)*twophotonL(i,3,2)
            write(lupri,'(A,F16.8)')'BTERM-contribution from Y*XZ L:',Bmoment(i,2)*twophotonL(i,1,3)
            write(lupri,'(A,F16.8)')'BTERM-contribution from Y*ZX L:',Bmoment(i,2)*twophotonL(i,3,1)
            write(lupri,'(A,F16.8)')'BTERM-contribution from Z*XY L:',Bmoment(i,3)*twophotonL(i,1,2)
            write(lupri,'(A,F16.8)')'BTERM-contribution from Z*YX L:',Bmoment(i,3)*twophotonL(i,2,1)
            MCDBtermL(i)=-Bmoment(i,1)*twophotonL(i,2,3)+Bmoment(i,1)*twophotonL(i,3,2)&
                 &-Bmoment(i,2)*twophotonL(i,3,1)+Bmoment(i,2)*twophotonL(i,1,3)&
                 &-Bmoment(i,3)*twophotonL(i,1,2)+Bmoment(i,3)*twophotonL(i,2,1)
            
            write(lupri,'(A,F16.8)')'LAO    BTERM-result  (a.u.)      :', MCDBtermL(i)
            write(lupri,'(A,F16.8)')'LAO    BTERM-result  (std units) :', BtermConstant*MCDBtermL(i)
            !   IF(MCDexci2(i).NE. 0)then
            !      MCDBtermL(i)=2E0_realk*MCDBtermL(i)
            !      write(lupri,'(A,F16.8)')'this state is degenerate therefore'
            !      write(lupri,'(A,F16.8)')'LAO    BTERM-result  (a.u.)      :', MCDBtermL(i)
            !      write(lupri,'(A,F16.8)')'LAO    BTERM-result  (std units) :', BtermConstant*MCDBtermL(i)
            !   endif
            write(lupri,'(A)')' '
            write(lupri,'(A,F16.8)')'SR-contribution from XY :    ',twophoton(i,1,2)
            write(lupri,'(A,F16.8)')'SR-contribution from XZ :    ',twophoton(i,1,3)
            write(lupri,'(A,F16.8)')'SR-contribution from YX :    ',twophoton(i,2,1)
            write(lupri,'(A,F16.8)')'SR-contribution from YZ :    ',twophoton(i,2,3)
            write(lupri,'(A,F16.8)')'SR-contribution from ZX :    ',twophoton(i,3,1)
            write(lupri,'(A,F16.8)')'SR-contribution from ZY :    ',twophoton(i,3,2)
            write(lupri,'(A,F16.8)')'BTERM-contribution from X*YZ :',Bmoment(i,1)*twophoton(i,2,3)
            write(lupri,'(A,F16.8)')'BTERM-contribution from X*ZY :',Bmoment(i,1)*twophoton(i,3,2)
            write(lupri,'(A,F16.8)')'BTERM-contribution from Y*XZ :',Bmoment(i,2)*twophoton(i,1,3)
            write(lupri,'(A,F16.8)')'BTERM-contribution from Y*ZX :',Bmoment(i,2)*twophoton(i,3,1)
            write(lupri,'(A,F16.8)')'BTERM-contribution from Z*XY :',Bmoment(i,3)*twophoton(i,1,2)
            write(lupri,'(A,F16.8)')'BTERM-contribution from Z*YX :',Bmoment(i,3)*twophoton(i,2,1)
            MCDBterm(i)=-Bmoment(i,1)*twophoton(i,2,3)+Bmoment(i,1)*twophoton(i,3,2)&
                 &-Bmoment(i,2)*twophoton(i,3,1)+Bmoment(i,2)*twophoton(i,1,3)&
                 &-Bmoment(i,3)*twophoton(i,1,2)+Bmoment(i,3)*twophoton(i,2,1)
            write(lupri,'(A,F16.8)')'NONLAO BTERM-result  (a.u.)      :', MCDBterm(i)
            write(lupri,'(A,F16.8)')'NONLAO BTERM-result  (std units) :', BtermConstant*MCDBterm(i)
            !   IF(MCDexci2(i).NE. 0)then
            !      MCDBterm(i)=2E0_realk*MCDBterm(i)
            !      write(lupri,'(A,F16.8)')'this state is degenerate therefore'
            !      write(lupri,'(A,F16.8)')'NONLAO BTERM-result  (a.u.)      :', MCDBterm(i)
            !      write(lupri,'(A,F16.8)')'NONLAO BTERM-result  (std units) :', BtermConstant*MCDBterm(i)
            !   endif
         enddo
         
         allocate(MCDAterm(nAterms))
         allocate(MCDAtermL(nAterms))
         do i=1,nAterms
            write(lupri,'(A,F16.8)')'<i|mx|j>',Atermmoment(i,1)
            write(lupri,'(A,F16.8)')'<i|my|j>',Atermmoment(i,2)
            write(lupri,'(A,F16.8)')'<i|mz|j>',Atermmoment(i,3)
            write(lupri,'(A,F16.8)')'London <i|mx|j>',AtermmomentL(i,1)
            write(lupri,'(A,F16.8)')'London <i|my|j>',AtermmomentL(i,2)
            write(lupri,'(A,F16.8)')'London <i|mz|j>',AtermmomentL(i,3)
            MCDAtermL(i) = 0E0_realk  
            MCDAterm(i) = 0E0_realk  
            do j=1,3
               SIGNL =  1E0_realk
               IF(AtermmomentL(i,j).LT. 0E0_realk) SIGNL = -1E0_realk
               eigenvalueL(j) = ABS(AtermmomentL(i,j))
               SIGN =  1E0_realk
               IF(Atermmoment(i,j).LT. 0E0_realk) SIGN = -1E0_realk
               eigenvalue(j) = ABS(Atermmoment(i,j))
               if(j.eq. 3)then
                  k=1;l=2
               elseif(j.eq. 2)then
                  k=3;l=1
               else !j=1
                  k=2;l=3
               endif
               X1=Amoment(i,1,k); Y1=Amoment(i,1,l); X2=Amoment(i,2,k); Y2=Amoment(i,2,l)          
               !there is 2 contribution of same magnitude so i need a factor 2
               !but I also need to divide with the degeneracy so it fits already. 
               write(lupri,'(a26,F16.8)')'eigenvalue(j)',eigenvalue(j)
               write(lupri,'(a26,F16.8)')'eigenvalueL(j)',eigenvalueL(j)
               write(lupri,'(a26,F16.8)')'transitionmomen',(X1*Y2+X2*Y1)
               write(lupri,'(a26,F16.8)')'X1',X1
               write(lupri,'(a26,F16.8)')'X2',X2
               write(lupri,'(a26,F16.8)')'Y1',Y1
               write(lupri,'(a26,F16.8)')'Y2',Y2
               write(lupri,'(a26,F16.8)')'SIGN',SIGN
               write(lupri,'(a26,F16.8)')'SIGNL',SIGNL
               MCDAterm(i) = MCDAterm(i)-SIGN*(X1*Y2-X2*Y1)*eigenvalue(j)
               MCDAtermL(i) = MCDAtermL(i)-SIGNL*(X1*Y2-X2*Y1)*eigenvalueL(j)
            enddo
            write(lupri,'(A,F12.8)')'Excitation Energy (au)        :', freqA(i)
            write(lupri,'(A,F12.8)')'Excitation Energy (eV)        :', freqA(i)*hartree
            write(lupri,'(A,F12.8)')'Excitation Energy (nm)        :', hcConstant_eVnm/(freqA(i)*hartree)
            write(lupri,'(A,F16.4)')'Excitation Energy (cm^-1)     :', freqA(i)*hartree*eVTOcm
            write(lupri,'(A,F16.8)')'   london A term (au)         :', MCDAtermL(i)
            write(lupri,'(A,F16.8)')'   london A term (std units)  :', AtermConstant*MCDAtermL(i)
            write(lupri,'(A,F16.8)')'nonlondon A term (au)         :', MCDAterm(i)
            write(lupri,'(A,F16.8)')'nonlondon A term (std units)  :', AtermConstant*MCDAterm(i)
         enddo
         !    call print_response_function(decomp,RspFunc_stack, n_RspFunc, .true.)

         deallocate(RspFunc_stack)
         call rsp_eq_sol_empty()
         nBterms = MCDnexci    
      endif
      !#################################################################
      !#   SIMULATED MCD SPECTRA BASED ON CALCULATED A AND B TERMS   
      !#################################################################
      IF(MCDinput%simulate)THEN
         call mem_alloc(GAMMAATERM,nAterms)
         call mem_alloc(GAMMABTERM,nBterms)
         do i=1,nAterms
            GAMMAATERM(i) = GAMMA
         enddo
         do i=1,nBterms
            GAMMABTERM(i) = GAMMA
         enddo
         call simulateSpectra(lupri,nBterms,nAterms,MCDBterm,MCDBtermL,&
              & MCDAterm,MCDAtermL,freqA,freqB,LORENTZ,MCDinput%Nsteps,&
              & GAMMAATERM,GAMMABTERM)
         call mem_dealloc(GAMMAATERM)
         call mem_dealloc(GAMMABTERM)
      endif
   endif
   if(MCDinput%dampedMCD)then
      !#################################################################
      !#    DO THE DAMPED MCD atomic unit SPECTRA CALCULATION   
      !#################################################################
      IF(MCDinput%nXcoor.EQ. 0)THEN
         !BUILD XCOOR BASED ON Standard MCD calculation: 
         !we determine points around the MCD peaks for which to calculate the 
         !damped MCD spectra
         allocate(Deltavec(MCDinput%nVecForPeak))
         do i=1,MCDinput%nVecForPeak-1
            Deltavec(i) = i*(1E0_realk/MCDinput%nVecForPeak)  
         enddo
         Deltavec(MCDinput%nVecForPeak) = 0.998E0_realk 
         allocate(Xcoor2(2*MCDinput%nVecForPeak*MCDnexci))
         Xcoor2 = 0E0_realk
         nXcoor2 = 0
         if(LORENTZ)then
            do i=1,MCDnexci
               do j=1,MCDinput%nVecForPeak
                  nXcoor2 = nXcoor2+1
                  Xcoor2(nXcoor2) = -sqrt(GAMMA*GAMMA*(1E0_realk-DELTAVEC(j))/(DELTAVEC(j)))+freqB(i)
               enddo
               do j=MCDinput%nVecForPeak-1,1,-1
                  nXcoor2 = nXcoor2+1
                  Xcoor2(nXcoor2) = sqrt(GAMMA*GAMMA*(1E0_realk-DELTAVEC(j))/(DELTAVEC(j)))+freqB(i)
               enddo
            enddo
         else
            do i=1,MCDnexci
               do j=1,MCDinput%nVecForPeak
                  nXcoor2 = nXcoor2+1
                  Xcoor2(nXcoor2) = freqB(i)-sqrt(-2*GAMMA*GAMMA*log(DELTAVEC(j)))
               enddo
               do j=MCDinput%nVecForPeak-1,1,-1
                  nXcoor2 = nXcoor2+1
                  Xcoor2(nXcoor2) = freqB(i)+sqrt(-2*GAMMA*GAMMA*log(DELTAVEC(j)))
               enddo
            enddo
         endif
         
         allocate(Xcoor3(nXcoor2))
         nXcoor3 = 0 
         do i=1,nXcoor2
            unique=.TRUE.
            do j=1,nXcoor3
               IF(ABS(Xcoor3(j)-Xcoor2(i)).LT. 1.0E-7_realk)unique=.FALSE.
            enddo
            if(unique)then
               nXcoor3 = nXcoor3+1
               Xcoor3(nXcoor3) = Xcoor2(i)
            endif
         enddo
          call INSERTION_increasing(Xcoor3(1:nXcoor3))
          deallocate(Xcoor2)
          deallocate(deltavec)
      else
         !We use points given in input using the .DAMPEDXCOOR keyword
         nXcoor3 = MCDinput%nXcoor
         allocate(Xcoor3(nXcoor3))
         do I =1,nXcoor3
            Xcoor3(I)=MCDinput%Xcoor(I)
         enddo
      endif
      if(nXcoor3.GT. 0)then
         n_rspfunc = nXcoor3*6
         allocate(RspFunc_stack2(n_rspfunc))
         k=1
         !gamma = 0E0_realk
         do i=1,nXcoor3
            !as both 2 and 3 is EL there is a symmetry connection 
            !<<mag,el,el>>(-omega,omega)
            !<< X, Y, Z>>(-omega,omega) = -<< X, Z, Y>>(-omega,omega)   ??????
            call set_verdet_RspFunc(RspFunc_stack2(k),1,2,3,Xcoor3(i),gamma,'MAG ')
            !   call set_verdet_RspFunc(RspFunc_stack2(k+1),1,3,2,Xcoor3(i),gamma,'MAG ') !not needed ?
            call set_verdet_RspFunc(RspFunc_stack2(k+1),2,3,1,Xcoor3(i),gamma,'MAG ')
            !   call set_verdet_RspFunc(RspFunc_stack2(k+3),2,1,3,Xcoor3(i),gamma,'MAG ') !not needed ?
            call set_verdet_RspFunc(RspFunc_stack2(k+2),3,1,2,Xcoor3(i),gamma,'MAG ')
            !   call set_verdet_RspFunc(RspFunc_stack2(k+5),3,2,1,Xcoor3(i),gamma,'MAG ') !not needed ?
            k=k+3
         enddo
         do i=1,nXcoor3
            call set_verdet_RspFunc(RspFunc_stack2(k),1,2,3,Xcoor3(i),gamma,'MAGO')
            !   call set_verdet_RspFunc(RspFunc_stack2(k+1),1,3,2,Xcoor3(i),gamma,'MAGO') !not needed ?
            call set_verdet_RspFunc(RspFunc_stack2(k+1),2,3,1,Xcoor3(i),gamma,'MAGO')
            !   call set_verdet_RspFunc(RspFunc_stack2(k+3),2,1,3,Xcoor3(i),gamma,'MAGO') !not needed ?
            call set_verdet_RspFunc(RspFunc_stack2(k+2),3,1,2,Xcoor3(i),gamma,'MAGO')
            !   call set_verdet_RspFunc(RspFunc_stack2(k+5),3,2,1,Xcoor3(i),gamma,'MAGO') !not needed ?
            k=k+3
         enddo
         n_rspfunc = k-1
         call rspfunc_or_transmoment(molcfg,F,D,S,RspFunc_stack2,n_rspfunc)
!         call print_response_func(RspFunc_stack2, n_RspFunc,.true.,molcfg%lupri)

         allocate(verdetL(n_rspfunc))
         verdetL=0E0_realk
         allocate(verdet(n_rspfunc))
         verdet=0E0_realk
         k=1
         do i=1,nXcoor3
            verdetL(i)=2*RspFunc_stack2(k)%result(1)+2*RspFunc_stack2(k+1)%result(1)+2*RspFunc_stack2(k+2)%result(1)
            k=k+3
         enddo
         do i=1,nXcoor3
            verdet(i)=2*RspFunc_stack2(k)%result(1)+2*RspFunc_stack2(k+1)%result(1)+2*RspFunc_stack2(k+2)%result(1)
            k=k+3
         enddo
         !#################################################################
         !#    PRINT THE Atomic unit damped MCD SPECTRA   
         !#################################################################
         LUMCD3=-1
         CALL LSOPEN(LUMCD3,'dampedMCDspectraRAU.dat','REPLACE','FORMATTED')
         WRITE(lupri,*)'Damped MCD Output (AU)'
         DO I=1,nXcoor3
            FA =  TWOPIM1*Xcoor3(I)
            WRITE(LUMCD3,*)Xcoor3(I),FA*AIMAG(verdetL(i)),FA*AIMAG(verdet(i))
            WRITE(LUPRI,'(3ES22.10)')Xcoor3(I),FA*AIMAG(verdetL(i)),FA*AIMAG(verdet(i))
         ENDDO
         CALL LSCLOSE(LUMCD3,'KEEP')

         
         LUMCD4=-1
         CALL LSOPEN(LUMCD4,'gnuplot_dampedMCDAU.gnu','REPLACE','FORMATTED')
         WRITE(LUMCD4,'(A)')'reset'
         WRITE(LUMCD4,'(A)')'set title ''MCD Spectra with and without london orbitals obtained using damped theory'' '
         WRITE(LUMCD4,'(A)')'set terminal postscript eps enhanced color'
         WRITE(LUMCD4,'(A)')'set output ''dampedMCDspectraAU.ps'' '
         WRITE(LUMCD4,'(A)')'set style data linespoints'
         WRITE(LUMCD4,'(A)')'#set linestyle  3' 
         WRITE(LUMCD4,'(A)')'#set linestyle  4'
         WRITE(LUMCD4,'(A)')'#set linestyle  7'
         WRITE(LUMCD4,'(A)')'set size 2.0,1.0'
         WRITE(LUMCD4,'(A)')'set style data linespoints'
         WRITE(LUMCD4,'(A,F12.4,A,F12.4,A)')'set xrange [',Xcoor3(1),':',Xcoor3(nXcoor3),']'
         WRITE(LUMCD4,'(A)')'P(x)=0.000'
         WRITE(LUMCD4,'(A)')'plot ''dampedMCDspectraRAU.dat'' using 1:2 title ''london'' w l ls 3,'//Achar(92)
         WRITE(LUMCD4,'(A)')'     ''dampedMCDspectraRAU.dat'' using 1:3 title ''nolondon'' w l ls 4,'//Achar(92)
         WRITE(LUMCD4,'(TL11,A)')'     P(x) w l ls 7'
         CALL LSCLOSE(LUMCD4,'KEEP')

         LUMCD4=-1
         CALL LSOPEN(LUMCD4,'gnuplot_compareMCDAU.gnu','REPLACE','FORMATTED')
         WRITE(LUMCD4,'(A)')'reset'
         WRITE(LUMCD4,'(A)')'set title ''comparison of MCD Spectra obtained from damped theory vs. explicit calculations'' '
         WRITE(LUMCD4,'(A)')'set terminal postscript eps enhanced color'
         WRITE(LUMCD4,'(A)')'set output ''comparedMCDspectraAU.ps'' '
         WRITE(LUMCD4,'(A)')'set style data linespoints'
         WRITE(LUMCD4,'(A)')'#set linestyle  1' 
         WRITE(LUMCD4,'(A)')'#set linestyle  2'
         WRITE(LUMCD4,'(A)')'#set linestyle  3' 
         WRITE(LUMCD4,'(A)')'#set linestyle  4'
         WRITE(LUMCD4,'(A)')'#set linestyle  7'
         WRITE(LUMCD4,'(A)')'set size 2.0,1.0'
         WRITE(LUMCD4,'(A)')'set style data linespoints'
         WRITE(LUMCD4,'(A,F12.4,A,F12.4,A)')'set xrange [',Xcoor3(1),':',Xcoor3(nXcoor3),']'
         WRITE(LUMCD4,'(A)')'P(x)=0.000'
         WRITE(LUMCD4,'(A)')'plot ''MCDspectraAU.dat'' using 1:2 title ''london'' w l ls 1,'//Achar(92)
         WRITE(LUMCD4,'(A)')'     ''MCDspectraAU.dat'' using 1:3 title ''nolondon'' w l ls 2,'//Achar(92)
         WRITE(LUMCD4,'(A)')'     ''dampedMCDspectraRAU.dat'' using 1:2 title ''damped-london'' w l ls 3,'//Achar(92)
         WRITE(LUMCD4,'(A)')'     ''dampedMCDspectraRAU.dat'' using 1:3 title ''damped-nolondon'' w l ls 4,'//Achar(92)
         WRITE(LUMCD4,'(TL11,A)')'     P(x) w l ls 7'
         CALL LSCLOSE(LUMCD4,'KEEP')
         !#################################################################
         !#    PRINT THE Standard unit damped MCD SPECTRA   
         !#################################################################
         LUMCD3=-1
         CALL LSOPEN(LUMCD3,'dampedMCDspectraR.dat','REPLACE','FORMATTED')
         WRITE(lupri,*)'Damped MCD Output '
         DO I=1,nXcoor3
            FA =  TWOPIM1*Xcoor3(I)*verdetConstant
            WRITE(LUMCD3,*)Xcoor3(I)*auTOcm,FA*AIMAG(verdetL(i)),FA*AIMAG(verdet(i))
            WRITE(LUPRI,'(3ES22.10)')Xcoor3(I)*auTOcm,FA*AIMAG(verdetL(i)),FA*AIMAG(verdet(i))
         ENDDO
         CALL LSCLOSE(LUMCD3,'KEEP')

         LUMCD4=-1
         CALL LSOPEN(LUMCD4,'gnuplot_dampedMCD.gnu','REPLACE','FORMATTED')
         WRITE(LUMCD4,'(A)')'reset'
         WRITE(LUMCD4,'(A)')'set title ''MCD Spectra with and without london orbitals obtained using damped theory'' '
         WRITE(LUMCD4,'(A)')'set terminal postscript eps enhanced color'
         WRITE(LUMCD4,'(A)')'set output ''dampedMCDspectra.ps'' '
         WRITE(LUMCD4,'(A)')'set style data linespoints'
         WRITE(LUMCD4,'(A)')'#set linestyle  3' 
         WRITE(LUMCD4,'(A)')'#set linestyle  4'
         WRITE(LUMCD4,'(A)')'#set linestyle  7'
         WRITE(LUMCD4,'(A)')'set size 2.0,1.0'
         WRITE(LUMCD4,'(A)')'set style data linespoints'
         WRITE(LUMCD4,'(A,F9.2,A,F9.2,A)')'set xrange [',Xcoor3(1)*auTOcm,':',Xcoor3(nXcoor3)*auTOcm,']'
         WRITE(LUMCD4,'(A)')'P(x)=0.000'
         WRITE(LUMCD4,'(A)')'plot ''dampedMCDspectraR.dat'' using 1:2 title ''london'' w l ls 3,'//Achar(92)
         WRITE(LUMCD4,'(A)')'     ''dampedMCDspectraR.dat'' using 1:3 title ''nolondon'' w l ls 4,'//Achar(92)
         WRITE(LUMCD4,'(TL11,A)')'     P(x) w l ls 7'
         CALL LSCLOSE(LUMCD4,'KEEP')

         LUMCD4=-1
         CALL LSOPEN(LUMCD4,'gnuplot_compareMCD.gnu','REPLACE','FORMATTED')
         WRITE(LUMCD4,'(A)')'reset'
         WRITE(LUMCD4,'(A)')'set title ''comparison of MCD Spectra obtained from damped theory vs. explicit calculations'' '
         WRITE(LUMCD4,'(A)')'set terminal postscript eps enhanced color'
         WRITE(LUMCD4,'(A)')'set output ''comparedMCDspectra.ps'' '
         WRITE(LUMCD4,'(A)')'set style data linespoints'
         WRITE(LUMCD4,'(A)')'#set linestyle  1' 
         WRITE(LUMCD4,'(A)')'#set linestyle  2'
         WRITE(LUMCD4,'(A)')'#set linestyle  3' 
         WRITE(LUMCD4,'(A)')'#set linestyle  4'
         WRITE(LUMCD4,'(A)')'#set linestyle  7'
         WRITE(LUMCD4,'(A)')'set size 2.0,1.0'
         WRITE(LUMCD4,'(A)')'set style data linespoints'
         WRITE(LUMCD4,'(A,F9.2,A,F9.2,A)')'set xrange [',Xcoor3(1)*auTOcm,':',Xcoor3(nXcoor3)*auTOcm,']'
         WRITE(LUMCD4,'(A)')'P(x)=0.000'
         WRITE(LUMCD4,'(A)')'plot ''MCDspectra.dat'' using 1:2 title ''london'' w l ls 1,'//Achar(92)
         WRITE(LUMCD4,'(A)')'     ''MCDspectra.dat'' using 1:3 title ''nolondon'' w l ls 2,'//Achar(92)
         WRITE(LUMCD4,'(A)')'     ''dampedMCDspectraR.dat'' using 1:2 title ''damped-london'' w l ls 3,'//Achar(92)
         WRITE(LUMCD4,'(A)')'     ''dampedMCDspectraR.dat'' using 1:3 title ''damped-nolondon'' w l ls 4,'//Achar(92)
         WRITE(LUMCD4,'(TL11,A)')'     P(x) w l ls 7'
         CALL LSCLOSE(LUMCD4,'KEEP')
         deallocate(verdetL)
         deallocate(verdet)
         deallocate(Xcoor3)
         deallocate(RspFunc_stack2)
      endif
      !Print the final result
   endif
   IF(MCDnexci.GT. 0)then
      write(lupri,'(A)')'The B terms of Magnetic Circular Dicroism'
      do i=1,MCDnexci    
         write(lupri,'(A)')' '
         write(lupri,'(A,F12.8)')'Excitation Energy (a.u.)         :', freqB(i)
         write(lupri,'(A,F12.8)')'Excitation Energy (eV)           :', freqB(i)*hartree
         write(lupri,'(A,F12.8)')'Excitation Energy (nm)           :', hcConstant_eVnm/(freqB(i)*hartree)
         write(lupri,'(A,F16.4)')'Excitation Energy (cm^-1)        :', freqB(i)*hartree*eVTOcm
         write(lupri,'(A)')' '
         write(lupri,'(A,F16.8)')'LAO    BTERM-result  (a.u.)      :', MCDBtermL(i)
         write(lupri,'(A,F16.8)')'LAO    BTERM-result  (std units) :', BtermConstant*MCDBtermL(i)
         write(lupri,'(A)')' '
         write(lupri,'(A,F16.8)')'NONLAO BTERM-result  (a.u.)      :', MCDBterm(i)
         write(lupri,'(A,F16.8)')'NONLAO BTERM-result  (std units) :', BtermConstant*MCDBterm(i)
      enddo
      
      IF(nAterms.GT.0)write(lupri,'(A)')'The A terms of Magnetic Circular Dicroism'
      do i=1,nAterms
         write(lupri,'(A,F12.8)')'Excitation Energy (a.u.)         :', freqA(i)
         write(lupri,'(A,F12.8)')'Excitation Energy (eV)           :', freqA(i)*hartree
         write(lupri,'(A,F12.8)')'Excitation Energy (nm)           :', hcConstant_eVnm/(freqA(i)*hartree)
         write(lupri,'(A,F16.4)')'Excitation Energy (cm^-1)        :', freqA(i)*hartree*eVTOcm
         write(lupri,'(a26,F16.8)')'   london A term (a.u.)        :',MCDAtermL(i)
         write(lupri,'(a26,F16.8)')'   london A term (std units)   :',AtermConstant*MCDAtermL(i)
         write(lupri,'(a26,F16.8)')'nonlondon A term (a.u.)        :',MCDAterm(i)
         write(lupri,'(a26,F16.8)')'nonlondon A term (std units)   :',AtermConstant*MCDAterm(i)
         write(lupri,'(A)')' '
      enddo
   else
      write(lupri,'(A)')'There are no A or B terms of Magnetic Circular Dicroism'
      write(lupri,'(A)')'As none of the investigated excitations were allowed'
   endif

 end subroutine MCDresponse_driver

SUBROUTINE INSERTION_increasing(a)
IMPLICIT NONE
real(realk), INTENT(INOUT) :: a(:)
INTEGER                :: I,J
real(realk) :: temp

Do I = 2, Size(a)
  temp = a(I)
  DO J = I-1, 1, -1
     IF(temp .LT. a(J)) Then
        a(J+1) = a(J)
     ELSE
        Exit
     ENDIF
  ENDDO
  a(J+1) = temp
End Do

END SUBROUTINE INSERTION_increasing

  !> \brief main NMR nuclear magnetic shielding tensor driver
  !> \author Thomas Kjaergaard
  !> \date 2010-03
  !>
  !> Questions regarding this driver, the input structure etc. may be addressed to tkjaergaard@chem.au.dk.
  !> The result is placed in the tensor NMST which has dimensions (3*natoms)*3
subroutine NMRshieldresponse_driver(molcfg,F,D,S)
implicit none
type(rsp_molcfg), intent(inout) :: molcfg
type(Matrix),intent(in) :: F,D,S
!
complex(8),allocatable       :: NMST(:,:)
type(Matrix)                 :: Dx(3),Fx(3)
real(realk)                  :: Factor
integer                      :: natoms,icoor,jcoor,k,lupri,nbast
Character(len=4),allocatable :: atomName(:)
real(realk)         :: TS,TE
character(len=1)        :: CHRXYZ(-3:3)
DATA CHRXYZ /'z','y','x',' ','X','Y','Z'/
nbast = D%nrow
lupri = molcfg%lupri
natoms = molcfg%natoms
CALL LSTIMER('START ',TS,TE,LUPRI)

WRITE(LUPRI,*) '=========================================================='
WRITE(LUPRI,*) '      NUCLEAR MAGNETIC SHIELDING(NMS) TENSOR RESULTS   '
WRITE(LUPRI,*) '=========================================================='

call pert_dens(molcfg,S,(/'MAG '/),(/3/),(/D/),(/F/),Dx,Fx,freq=(/(0d0,0d0)/))
Fx=0
Factor=53.25135351884770E0_realk!53.2513539566280 !1e6*alpha^2 
allocate(NMST(3*natoms,3))
do icoor=1,3
   do jcoor=1,3*natoms  
      NMST(jcoor,icoor) = 0.0E0_realk
   enddo
enddo
call rsp_oneave(molcfg, S,(/'NUCM'/),(/Dx/),(/3*natoms,3/), NMST)
call rsp_oneave(molcfg, S,(/'NUCM','MAG '/),(/D/),(/3*natoms,3/), NMST)
Dx=0
do icoor=1,3
   do jcoor=1,3*natoms  
      NMST(jcoor,icoor) = factor*NMST(jcoor,icoor)
   enddo
enddo

WRITE(LUPRI,*) " Total shielding tensor"
allocate(atomname(natoms))
do jcoor=1,natoms  
   atomname(jcoor)=molcfg%setting%molecule(1)%p%ATOM(jcoor)%Name
enddo

do jcoor=1,natoms  
   WRITE (LUPRI,'(20X,3(A,13X),/)') 'Bx', 'By', 'Bz'
   WRITE (LUPRI,'(2X,A8,3F15.8)')  atomname(jcoor)//CHRXYZ(1),&
   & (REAL(NMST(3*jCOOR-2,K)),K=1,3)
   WRITE (LUPRI,'(2X,A8,3F15.8)')  atomname(jcoor)//CHRXYZ(2),&
   & (REAL(NMST(3*jCOOR-1,K)),K=1,3)
   WRITE (LUPRI,'(2X,A8,3F15.8)')  atomname(jcoor)//CHRXYZ(3),&
   & (REAL(NMST(3*jCOOR,K)),K=1,3)
enddo
WRITE(LUPRI,*) "      Absolute chemical shift"
WRITE(LUPRI,*) "number of atoms"
WRITE(LUPRI,*) natoms, "                          "
do jcoor=1,natoms  
   WRITE (LUPRI,'(2X,A8,f15.8)')  atomname(jcoor), 1.d0/3.d0*&
   &(REAL(NMST(3*jCOOR-2,1))+REAL(NMST(3*jCOOR-1,2))+REAL(NMST(3*jCOOR,3)))
enddo
deallocate(NMST)
deallocate(atomname)
CALL LSTIMER('*SHIELD',TS,TE,LUPRI)
WRITE(LUPRI,*) " Done with shielding tensor calculation"

end subroutine NMRshieldresponse_driver

subroutine simulateSpectra(lupri,nBterms,nAterms,MCDBterm,MCDBtermL,MCDAterm,&
     &MCDAtermL,freqA,freqB,LORENTZ,Nsteps,GAMMAATERM,GAMMABTERM)
  implicit none           
  !> logical unit number of output 
  integer,intent(in)     :: lupri
  integer,intent(in)     :: Nsteps,nAterms,nBterms
  !> MCD B terms in atomic units
  real(realk),intent(in) :: MCDBterm(nBterms),MCDBtermL(nBterms)
  !> MCD A terms in atomic units
  real(realk),intent(in) :: MCDAterm(nAterms),MCDAtermL(nAterms)
  !> MCD A and B term freqs
  real(realk),intent(in) :: freqA(nAterms),freqB(nBterms)
  real(realk),intent(in) :: GAMMAATERM(nAterms)
  real(realk),intent(in) :: GAMMABTERM(nBterms)
  logical,intent(in)     :: LORENTZ
!
  integer   :: i,j,k,l,LUMCD,LUMCD2,Nstep,lengthX
  real(realk),allocatable  :: Y(:),YL(:)
  real(realk),allocatable :: PUREB(:,:),PUREBL(:,:),PUREA(:,:),PUREAL(:,:)
  real(realk),allocatable :: Xcoor(:),COMBINED(:),COMBINEDL(:),Deltavec(:),Xcoor2(:),Xcoor3(:)
  complex(realk),allocatable :: verdet(:),verdetL(:) 
  real(realk) :: Range,step,Estart,Eend,FUN,INT,x,DISPLACEMENT,FA
  character(len=16) :: filenameMCD
!
  Real(realk),parameter :: PI=3.14159265358979323846E0_realk
  !the Bterm is in standard units given in D**2*muB/(cm^-1) (au : e**3*a0**4/hbar)
  ! D is the unit debye        :   1 D = 0.393430307 (au.) (au : e*a0)
  ! muB is the Bohr magneton   : 1 muB = 0.5 (au.)
  ! and inverse centimeter     : 1cm^⁻1= 219475 (au.)
  ! so the standard units are : D**2*muB/(cm^-1) = 16985.974069713367 (au.)
  ! so to convert from au to standard we inverse this to obtain
  Real(realk),parameter :: BtermConstant=5.88764*1E-05_realk
  !the Aterm is in standard units given in D**2*muB (au : e**3*hbar*a0**2/me)
  !so the standard units are : D**2*muB = 0.07739370323305711 (au.)
  !so to convert from au to standard we inverse this to obtain
  Real(realk),parameter :: AtermConstant=12.920947806163007
  Real(realk),parameter :: VerdetConstant=5.88764*1E-05_realk*33.53!0.152123007*1.0E-7_realk
  !hc = 1239.8419300923944 eV*nm  
  !so  lambda (in nm) = hc/E (E in eV)
  !and E (in eV) = hc/lambda (lambda in nm)
  Real(realk),parameter :: hcConstant_eVnm=1239.8419300923944
  real(realk),parameter :: TWOPIM1 = 0.15915494309189535E0_realk
  real(realk),parameter :: hartree=27.21138386E0_realk
  Real(realk), parameter ::eVTOcm= 8065.54E0_realk !eV to cm^-1
  real(realk), parameter :: autocm=hartree*eVTOcm 

  print*,'freqA',freqA
  print*,'freqB',freqB
  print*,'GAMMAATERM',GAMMAATERM
  print*,'GAMMABTERM',GAMMABTERM
  print*,'MCDBterm',MCDBterm
  print*,'MCDAterm',MCDAterm
  print*,'MCDBtermL',MCDBtermL
  print*,'MCDAtermL',MCDAtermL

  !#################################################################
  !#   Atomic units UNITS
  !#################################################################
  IF(nAterms.NE. 0 .AND. nBterms.NE. 0)then
     Eend=MAX(freqA(nAterms),freqB(nBterms))+0.1E0_realk*MAX(freqA(nAterms),freqB(nBterms))
  ELSE
     IF(nAterms.NE. 0)THEN
        Eend=freqA(nAterms)+0.1E0_realk*freqA(nAterms)
     ELSEIF(nBterms.NE. 0)THEN
        Eend = freqB(nBterms)+0.1E0_realk*freqB(nBterms)
     ELSE
        Eend = -1
     ENDIF
  ENDIF
  IF(nAterms.NE. 0)THEN
     IF(nBterms.NE. 0)THEN
        Estart = MAX(MIN(freqA(1),freqB(1))-0.1E0_realk*MAX(freqA(nAterms),freqB(nBterms)),0E0_realk)
     ELSE
        Estart = MAX(freqA(1)-0.1E0_realk*MAX(freqA(nAterms),freqB(nBterms)),0E0_realk)
     ENDIF
  ELSE
     IF(nBterms.NE. 0)THEN
        Estart = MAX(freqB(1)-0.1E0_realk*freqB(nBterms),0E0_realk)
     ELSE
        Estart = 0
     ENDIF
  ENDIF
  IF(Eend.NE.-1)THEN
     Nstep = Nsteps
     Range = Eend-Estart
     step = Range/Nstep
     x=Estart
     lengthX=0
     DO WHILE(x<Eend)
        x=x+step
        lengthX=lengthX+1
     ENDDO
     Allocate(PUREB(lengthX,nBterms))
     Allocate(PUREBL(lengthX,nBterms))
     PUREB = 0E0_realk
     PUREBL = 0E0_realk
     Allocate(PUREA(lengthX,nAterms))
     Allocate(PUREAL(lengthX,nAterms))
     PUREA = 0E0_realk
     PUREAL = 0E0_realk
     Allocate(Xcoor(lengthX))
     Xcoor = 0E0_realk
     x=Estart
     lengthX=0
     DO WHILE(x<Eend)
        x=x+step
        lengthX=lengthX+1
        XCOOR(lengthX)=x
     ENDDO
     Allocate(COMBINED(lengthX))
     Allocate(COMBINEDL(lengthX))
     COMBINED = 0E0_realk
     COMBINEDL = 0E0_realk
     DISPLACEMENT = 0E0_realk          
     !#################################################################
     !#    THE Atomic unit B TERMS   
     !#################################################################
     IF(nBterms .NE. 0)THEN
        Allocate(Y(nBterms))
        Allocate(YL(nBterms))
        IF(LORENTZ)THEN
           !DETERMINE THE SCALING FACTORS FOR THE BTERMS 
           DO i=1,nBterms
              x=Estart   
              INT=0
              DO WHILE(x<Eend)
                 x=x+step
                 FUN=GAMMABTERM(i)/((x-freqB(i))**2+GAMMABTERM(i)**2)
                 INT=INT+(FUN/x)*step
              ENDDO
              Y(i)=-MCDBterm(i)/INT 
              YL(i)=-MCDBtermL(i)/INT
           ENDDO
           !DO THE LINE SHAPE
           DO i=1,nBterms
              x=Estart
              lengthX=0
              DO WHILE(x<Eend)
                 x=x+step
                 lengthX=lengthX+1
                 PUREB(lengthX,i)=Y(i)*GAMMABTERM(i)/((x-freqB(i))**2+GAMMABTERM(i)**2) 
                 PUREBL(lengthX,i)=YL(i)*GAMMABTERM(i)/((x-freqB(i))**2+GAMMABTERM(i)**2) 
                 COMBINED(lengthX) = COMBINED(lengthX) + PUREB(lengthX,i)  
                 COMBINEDL(lengthX) = COMBINEDL(lengthX) + PUREBL(lengthX,i)  
              ENDDO
           ENDDO
        ELSE !GAUSSIAN
           !DETERMINE THE SCALING FACTORS FOR THE BTERMS 
           DO i=1,nBterms
              x=Estart   
              INT=0
              DO WHILE(x<Eend)
                 x=x+step
                 FUN=1/(gammaBTERM(i)*sqrt(2*PI))*EXP(-((x-freqB(i))**2)/(2*gammaBTERM(i)*gammaBTERM(i))) 
                 INT=INT+(FUN/x)*step
              ENDDO
              !PLOT IN ATOMIC UNIT THEN ELSE
              Y(i)=-MCDBterm(i)/INT 
              YL(i)=-MCDBtermL(i)/INT 
           ENDDO
           !DO THE LINE SHAPE
           DO i=1,nBterms
              x=Estart
              lengthX=0
              DO WHILE(x<Eend)
                 x=x+step
                 lengthX=lengthX+1
                 PUREB(lengthX,i)=Y(i)/(gammaBTERM(i)*sqrt(2*PI))*&
                      & EXP(-((x-freqB(i))**2)/(2*gammaBTERM(i)*gammaBTERM(i))) 
                 PUREBL(lengthX,i)=YL(i)/(gammaBTERM(i)*sqrt(2*PI))*&
                      & EXP(-((x-freqB(i))**2)/(2*gammaBTERM(i)*gammaBTERM(i))) 
                 COMBINED(lengthX) = COMBINED(lengthX) + PUREB(lengthX,i)  
                 COMBINEDL(lengthX) = COMBINEDL(lengthX) + PUREBL(lengthX,i)  
              ENDDO
           ENDDO
        ENDIF
        DeAllocate(Y)
        DeAllocate(YL)
     ENDIF
     !#################################################################
     !#    THE Atomic unit A TERMS   
     !#################################################################
     IF(nAterms .NE. 0)THEN
        Allocate(Y(nAterms))
        Allocate(YL(nAterms))
        IF(LORENTZ)THEN
           !DETERMINE THE SCALING FACTORS FOR THE ATERMS 
           DO i=1,nAterms
              x=Estart   
              INT=0
              DO WHILE(x<Eend)
                 x=x+step
                 FUN=(2E0_realk*GAMMAATERM(i)*(x-freqA(i)))/((x-freqA(i))**2+GAMMAATERM(i)**2)**2
                 INT=INT+(x-freqA(i))*(FUN/x)*step
              ENDDO
              Y(i)=MCDAterm(i)/INT 
              YL(i)=MCDAtermL(i)/INT 
           ENDDO
           !DO THE LINE SHAPE
           DO i=1,nAterms
              x=Estart
              lengthX=0
              DO WHILE(x<Eend)
                 x=x+step
                 lengthX=lengthX+1
                 PUREA(lengthX,i)  =Y(i)*(2E0_realk*GAMMAATERM(i)*(x-freqA(i)))/((x-freqA(i))**2+GAMMAATERM(i)**2)**2
                 PUREAL(lengthX,i)=YL(i)*(2E0_realk*GAMMAATERM(i)*(x-freqA(i)))/((x-freqA(i))**2+GAMMAATERM(i)**2)**2
                 COMBINED(lengthX) = COMBINED(lengthX) + PUREA(lengthX,i)  
                 COMBINEDL(lengthX)= COMBINEDL(lengthX) + PUREAL(lengthX,i)  
              ENDDO
           ENDDO
        ELSE !GAUSSIAN
           !DETERMINE THE SCALING FACTORS FOR THE ATERMS 
           DO i=1,nAterms
              x=Estart   
              INT=0
              DO WHILE(x<Eend)
                 x=x+step
                 FUN=((x-freqA(i))/(gammaATERM(i)*gammaATERM(i)))/(gammaATERM(i)*sqrt(2*PI))*&
                      &EXP(-((x-freqA(i))**2)/(2*gammaATERM(i)*gammaATERM(i))) 
                 INT=INT+(x-freqA(i))*(FUN/x)*step
              ENDDO
              Y(i)=MCDAterm(i)/INT
              YL(i)=MCDAtermL(i)/INT
           ENDDO
           !DO THE LINE SHAPE
           DO i=1,nAterms
              x=Estart
              lengthX=0
              DO WHILE(x<Eend)
                 x=x+step
                 lengthX=lengthX+1
                 PUREA(lengthX,i)=Y(i)*((x-freqA(i))/(gammaATERM(i)*gammaATERM(i)))/(gammaATERM(i)*sqrt(2*PI))*&
                      &EXP(-((x-freqA(i))**2)/(2*gammaATERM(i)*gammaATERM(i))) 
                 PUREAL(lengthX,i)=YL(i)*((x-freqA(i))/(gammaATERM(i)*gammaATERM(i)))/(gammaATERM(i)*sqrt(2*PI))*&
                      &EXP(-((x-freqA(i))**2)/(2*gammaATERM(i)*gammaATERM(i))) 
                 COMBINED(lengthX) = COMBINED(lengthX) + PUREA(lengthX,i)  
                 COMBINEDL(lengthX)= COMBINEDL(lengthX)+ PUREAL(lengthX,i)  
              ENDDO
           ENDDO
        ENDIF
        DeAllocate(Y)
        DeAllocate(YL)
     ENDIF
     !#################################################################
     !#    PRINT THE Atomic unit MCD SPECTRA   
     !#################################################################
     LUMCD=-1
     CALL LSOPEN(LUMCD,'MCDspectraAU.dat','REPLACE','FORMATTED')
     DO I=1,lengthX
        !      WRITE(LUMCD,*)'',Xcoor(I),PUREB(1,I),PUREB(2,I),PUREB(3,I),COMBINED(I)
        !      WRITE(LUMCD,*)'',Xcoor(I),PUREA(1,I),COMBINED(I)
        WRITE(LUMCD,*)Xcoor(I),COMBINEDL(I),COMBINED(I)
     ENDDO
     CALL LSCLOSE(LUMCD,'KEEP')

     DO J=1,nBterms
        LUMCD=-1
        call generateMCDFilename(filenameMCD,'BtermAU',J)
        CALL LSOPEN(LUMCD,TRIM(filenameMCD),'REPLACE','FORMATTED')
        DO I=1,lengthX
           WRITE(LUMCD,*)Xcoor(I),PUREBL(I,J),PUREB(I,J)
        ENDDO
        CALL LSCLOSE(LUMCD,'KEEP')
     ENDDO

     DO J=1,nAterms
        LUMCD=-1
        call generateMCDFilename(filenameMCD,'AtermAU',J)
        CALL LSOPEN(LUMCD,TRIM(filenameMCD),'REPLACE','FORMATTED')
        DO I=1,lengthX
           WRITE(LUMCD,*)Xcoor(I),PUREAL(I,J),PUREA(I,J)
        ENDDO
        CALL LSCLOSE(LUMCD,'KEEP')
     ENDDO

     LUMCD2=-1
     CALL LSOPEN(LUMCD2,'gnuplot_MCDAU.gnu','REPLACE','FORMATTED')
     WRITE(LUMCD2,'(A)')'reset'
     WRITE(LUMCD2,'(A)')'set title '//Achar(39)//'MCD Spectra with london orbitals'//Achar(39)
     WRITE(LUMCD2,'(A)')'set terminal postscript eps enhanced color'
     WRITE(LUMCD2,'(A)')'set output ''simulatedMCDspectraAU.ps'' '
     WRITE(LUMCD2,'(A)')'set style data linespoints'
     WRITE(LUMCD2,'(A)')'#set linestyle  1' 
     WRITE(LUMCD2,'(A)')'#set linestyle  2'
     WRITE(LUMCD2,'(A)')'#set linestyle  7'
     WRITE(LUMCD2,'(A)')'set size 2.0,1.0'
     WRITE(LUMCD2,'(A)')'set style data linespoints'
     WRITE(LUMCD2,'(A,F12.4,A,F12.4,A)')'set xrange [',Xcoor(1),':',Xcoor(lengthX),']'
     WRITE(LUMCD2,'(A)')'P(x)=0.000'
     WRITE(LUMCD2,'(A)')'plot ''MCDspectraAU.dat'' using 1:2 title ''london'' w l ls 1,'//Achar(92)
     WRITE(LUMCD2,'(A)')'     ''MCDspectraAU.dat'' using 1:3 title ''nolondon'' w l ls 2,'//Achar(92)
     WRITE(LUMCD2,'(TL11,A)')'     P(x) w l ls 7'
     CALL LSCLOSE(LUMCD2,'KEEP')

     deAllocate(PUREB)
     deAllocate(PUREBL)
     deAllocate(PUREA)
     deAllocate(PUREAL)
     deAllocate(Xcoor)
     deAllocate(COMBINED)
     deAllocate(COMBINEDL)

  endif

  !#################################################################
  !#   Standard UNITS which means CM^-1 for the freq and 
  !#   molar ellipticity for the ellipticity
  !#################################################################
  IF(nAterms.NE. 0 .AND. nBterms.NE. 0)then
     Eend=MAX(freqA(nAterms)*auTOcm,freqB(nBterms)*auTOcm)+0.1E0_realk*MAX(freqA(nAterms)*auTOcm,freqB(nBterms)*auTOcm)
  ELSE
     IF(nAterms.NE. 0)THEN
        Eend=freqA(nAterms)*auTOcm+0.1E0_realk*freqA(nAterms)*auTOcm
     ELSEIF(nBterms.NE. 0)THEN
        Eend = freqB(nBterms)*auTOcm+0.1E0_realk*freqB(nBterms)*auTOcm
     ELSE
        Eend = -1
     ENDIF
  ENDIF
  IF(nAterms.NE. 0)THEN
     IF(nBterms.NE. 0)THEN
        Estart = MAX(MIN(freqA(1)*auTOcm,freqB(1)*auTOcm)-0.1E0_realk*MAX(freqA(nAterms)*auTOcm,freqB(nBterms)*auTOcm),0E0_realk)
     ELSE
        Estart = MAX(freqA(1)*auTOcm-0.1E0_realk*MAX(freqA(nAterms)*auTOcm,freqB(nBterms)*auTOcm),0E0_realk)
     ENDIF
  ELSE
     IF(nBterms.NE. 0)THEN
        Estart = MAX(freqB(1)*auTOcm-0.1E0_realk*freqB(nBterms)*auTOcm,0E0_realk)
     ELSE
        Estart = 0
     ENDIF
  ENDIF
  IF(Eend.NE.-1)THEN
     Range = Eend-Estart
     step = Range/Nstep
     x=Estart
     lengthX=0
     DO WHILE(x<Eend)
        x=x+step
        lengthX=lengthX+1
     ENDDO
     Allocate(PUREB(lengthX,nBterms))
     Allocate(PUREBL(lengthX,nBterms))
     PUREB = 0E0_realk
     PUREBL = 0E0_realk
     Allocate(PUREA(lengthX,nAterms))
     Allocate(PUREAL(lengthX,nAterms))
     PUREA = 0E0_realk
     PUREAL = 0E0_realk
     Allocate(Xcoor(lengthX))
     Xcoor = 0E0_realk
     x=Estart
     lengthX=0
     DO WHILE(x<Eend)
        x=x+step
        lengthX=lengthX+1
        XCOOR(lengthX)=x
     ENDDO
     Allocate(COMBINED(lengthX))
     Allocate(COMBINEDL(lengthX))
     COMBINED = 0E0_realk
     COMBINEDL = 0E0_realk
     DISPLACEMENT = 0E0_realk
     !#################################################################
     !#    THE B TERMS   
     !#################################################################
     IF(nBterms .NE. 0)THEN
        Allocate(Y(nBterms))
        Allocate(YL(nBterms))
        IF(LORENTZ)THEN
           !DETERMINE THE SCALING FACTORS FOR THE BTERMS 
           DO i=1,nBterms
              x=Estart   
              INT=0
              DO WHILE(x<Eend)
                 x=x+step
                 FUN=GAMMABTERM(i)*auTOcm/((x-freqB(i)*auTOcm)**2+(GAMMABTERM(i)*auTOcm)**2)
                 INT=INT+(FUN/x)*step
              ENDDO
              Y(i)=-BtermConstant*MCDBterm(i)*33.53/INT 
              YL(i)=-BtermConstant*MCDBtermL(i)*33.53/INT
           ENDDO
           !DO THE LINE SHAPE
           DO i=1,nBterms
              x=Estart
              lengthX=0
              DO WHILE(x<Eend)
                 x=x+step
                 lengthX=lengthX+1
                 PUREB(lengthX,i)=Y(i)*GAMMABTERM(i)*auTOcm/((x-freqB(i)*auTOcm)**2+(GAMMABTERM(i)*auTOcm)**2) 
                 PUREBL(lengthX,i)=YL(i)*GAMMABTERM(i)*auTOcm/((x-freqB(i)*auTOcm)**2+(GAMMABTERM(i)*auTOcm)**2) 
                 COMBINED(lengthX) = COMBINED(lengthX) + PUREB(lengthX,i)  
                 COMBINEDL(lengthX) = COMBINEDL(lengthX) + PUREBL(lengthX,i)  
              ENDDO
           ENDDO
        ELSE !gaussian
           !DETERMINE THE SCALING FACTORS FOR THE BTERMS 
           DO i=1,nBterms
              x=Estart   
              INT=0
              DO WHILE(x<Eend)
                 x=x+step
                 FUN=1/(gammaBTERM(i)*auTOcm*sqrt(2*PI))*&
                      & EXP(-((x-freqB(i)*auTOcm)**2)/(2*(gammaBTERM(i)*gammaBTERM(i)*auTOcm*auTOcm))) 
                 INT=INT+(FUN/x)*step
              ENDDO
              !PLOT IN ATOMIC UNIT THEN ELSE
              Y(i)=-BtermConstant*MCDBterm(i)*33.53/INT 
              YL(i)=-BtermConstant*MCDBtermL(i)*33.53/INT 
           ENDDO
           !DO THE LINE SHAPE
           DO i=1,nBterms
              x=Estart
              lengthX=0
              DO WHILE(x<Eend)
                 x=x+step
                 lengthX=lengthX+1
                 PUREB(lengthX,i)=Y(i)/(gammaBTERM(i)*auTOcm*sqrt(2*PI))*&
                      & EXP(-((x-freqB(i)*auTOcm)**2)/(2*gammaBTERM(i)*auTOcm*gammaBTERM(i)*auTOcm)) 
                 PUREBL(lengthX,i)=YL(i)/(gammaBTERM(i)*auTOcm*sqrt(2*PI))*&
                      & EXP(-((x-freqB(i)*auTOcm)**2)/(2*gammaBTERM(i)*auTOcm*gammaBTERM(i)*auTOcm)) 
                 COMBINED(lengthX) = COMBINED(lengthX) + PUREB(lengthX,i)  
                 COMBINEDL(lengthX) = COMBINEDL(lengthX) + PUREBL(lengthX,i)  
              ENDDO
           ENDDO
        ENDIF
        DeAllocate(Y)
        DeAllocate(YL)
     ENDIF
     !#################################################################
     !#    THE A TERMS   
     !#################################################################
     IF(nAterms .NE. 0)THEN
        Allocate(Y(nAterms))
        Allocate(YL(nAterms))
        IF(LORENTZ)THEN
           !DETERMINE THE SCALING FACTORS FOR THE ATERMS 
           DO i=1,nAterms
              x=Estart   
              INT=0
              DO WHILE(x<Eend)
                 x=x+step
                 FUN=(2E0_realk*GAMMAATERM(i)*auTOcm*(x-freqA(i)*auTOcm))/((x-freqA(i)*auTOcm)**2+(GAMMAATERM(i)*auTOcm)**2)**2
                 INT=INT+(x-freqA(i)*auTOcm)*(FUN/x)*step
              ENDDO
              Y(i)=AtermConstant*MCDAterm(i)*33.53/INT 
              YL(i)=AtermConstant*MCDAtermL(i)*33.53/INT 
           ENDDO
           !DO THE LINE SHAPE
           DO i=1,nAterms
              x=Estart
              lengthX=0
              DO WHILE(x<Eend)
                 x=x+step
                 lengthX=lengthX+1
                 PUREA(lengthX,i)  = Y(i)*(2E0_realk*GAMMAATERM(i)*auTOcm*(x-freqA(i)*auTOcm))/&
		 &((x-freqA(i)*auTOcm)**2+(GAMMAATERM(i)*auTOcm)**2)**2
                 PUREAL(lengthX,i) =YL(i)*(2E0_realk*GAMMAATERM(i)*auTOcm*(x-freqA(i)*auTOcm))/&
		 &((x-freqA(i)*auTOcm)**2+(GAMMAATERM(i)*auTOcm)**2)**2
                 COMBINED(lengthX) = COMBINED(lengthX) + PUREA(lengthX,i)  
                 COMBINEDL(lengthX)= COMBINEDL(lengthX) + PUREAL(lengthX,i)  
              ENDDO
           ENDDO
        ELSE !GAUSSIAN
           !DETERMINE THE SCALING FACTORS FOR THE ATERMS 
           DO i=1,nAterms
              x=Estart   
              INT=0
              DO WHILE(x<Eend)
                 x=x+step
                 FUN=((x-freqA(i)*auTOcm)/(gammaATERM(i)*auTOcm*gammaATERM(i)*auTOcm))/(gammaATERM(i)*auTOcm*sqrt(2*PI))*&
                      &EXP(-((x-freqA(i)*auTOcm)**2)/(2*gammaATERM(i)*auTOcm*gammaATERM(i)*auTOcm)) 
                 INT=INT+(x-freqA(i)*auTOcm)*(FUN/x)*step
              ENDDO
              Y(i)=AtermConstant*MCDAterm(i)*33.53/INT
              YL(i)=AtermConstant*MCDAtermL(i)*33.53/INT
           ENDDO
           !DO THE LINE SHAPE
           DO i=1,nAterms
              x=Estart
              lengthX=0
              DO WHILE(x<Eend)
                 x=x+step
                 lengthX=lengthX+1
                 PUREA(lengthX,i) = Y(i)*((x-freqA(i)*auTOcm)/(gammaATERM(i)*auTOcm*gammaATERM(i)*auTOcm))/&
		 &(gammaATERM(i)*auTOcm*sqrt(2*PI))*&
                      &EXP(-((x-freqA(i)*auTOcm)**2)/(2*gammaATERM(i)*auTOcm*gammaATERM(i)*auTOcm)) 
                 PUREAL(lengthX,i)=YL(i)*((x-freqA(i)*auTOcm)/(gammaATERM(i)*auTOcm*gammaATERM(i)*auTOcm))/&
		 &(gammaATERM(i)*auTOcm*sqrt(2*PI))*&
                      &EXP(-((x-freqA(i)*auTOcm)**2)/(2*gammaATERM(i)*auTOcm*gammaATERM(i)*auTOcm)) 
                 COMBINED(lengthX) = COMBINED(lengthX) + PUREA(lengthX,i)  
                 COMBINEDL(lengthX)= COMBINEDL(lengthX)+ PUREAL(lengthX,i)  
              ENDDO
           ENDDO
        ENDIF
        DeAllocate(Y)
        DeAllocate(YL)
     ENDIF
     !#################################################################
     !#    PRINT THE MCD SPECTRA   
     !#################################################################
     LUMCD=-1
     CALL LSOPEN(LUMCD,'MCDspectra.dat','REPLACE','FORMATTED')
     DO I=1,lengthX
        !      WRITE(LUMCD,*)'',Xcoor(I),PUREB(1,I),PUREB(2,I),PUREB(3,I),COMBINED(I)
        !      WRITE(LUMCD,*)'',Xcoor(I),PUREA(1,I),COMBINED(I)
        WRITE(LUMCD,*)Xcoor(I),COMBINEDL(I),COMBINED(I)
     ENDDO
     CALL LSCLOSE(LUMCD,'KEEP')

     DO J=1,nBterms     
        LUMCD=-1
        call generateMCDFilename(filenameMCD,'Bterm',J)
        CALL LSOPEN(LUMCD,TRIM(filenameMCD),'REPLACE','FORMATTED')
        DO I=1,lengthX
           WRITE(LUMCD,*)Xcoor(I),PUREBL(I,J),PUREB(I,J)
        ENDDO
        CALL LSCLOSE(LUMCD,'KEEP')
     ENDDO
     DO J=1,nAterms     
        LUMCD=-1
        call generateMCDFilename(filenameMCD,'Aterm',J)
        CALL LSOPEN(LUMCD,TRIM(filenameMCD),'REPLACE','FORMATTED')
        DO I=1,lengthX
           WRITE(LUMCD,*)Xcoor(I),PUREAL(I,J),PUREA(I,J)
        ENDDO
        CALL LSCLOSE(LUMCD,'KEEP')
     ENDDO

     LUMCD2=-1
     CALL LSOPEN(LUMCD2,'gnuplot_MCD.gnu','REPLACE','FORMATTED')
     WRITE(LUMCD2,'(A)')'reset'
     WRITE(LUMCD2,'(A)')'set title '//Achar(39)//'MCD Spectra with london orbitals'//Achar(39)
     WRITE(LUMCD2,'(A)')'set terminal postscript eps enhanced color'
     WRITE(LUMCD2,'(A)')'set output ''simulatedMCDspectra.ps'' '
     WRITE(LUMCD2,'(A)')'set style data linespoints'
     WRITE(LUMCD2,'(A)')'#set linestyle  1' 
     WRITE(LUMCD2,'(A)')'#set linestyle  2'
     WRITE(LUMCD2,'(A)')'#set linestyle  7'
     WRITE(LUMCD2,'(A)')'set size 2.0,1.0'
     WRITE(LUMCD2,'(A)')'set style data linespoints'
     WRITE(LUMCD2,'(A,F9.2,A,F9.2,A)')'set xrange [',Xcoor(1),':',Xcoor(lengthX),']'
     WRITE(LUMCD2,'(A)')'P(x)=0.000'
     WRITE(LUMCD2,'(A)')'plot ''MCDspectra.dat'' using 1:2 title ''london'' w l ls 1,'//Achar(92)
     WRITE(LUMCD2,'(A)')'     ''MCDspectra.dat'' using 1:3 title ''nolondon'' w l ls 2,'//Achar(92)
     WRITE(LUMCD2,'(TL11,A)')'     P(x) w l ls 7'
     CALL LSCLOSE(LUMCD2,'KEEP')
     deAllocate(PUREB)
     deAllocate(PUREBL)
     deAllocate(PUREA)
     deAllocate(PUREAL)
     deAllocate(Xcoor)
     deAllocate(COMBINED)
     deAllocate(COMBINEDL)
  endif
end subroutine simulateSpectra

subroutine generateMCDFilename(filenameBtermAU,SPEC,I)
  implicit none
  character(len=16) :: filenameBtermAU
  character*(*) :: SPEC !like BtermAU
  integer :: I,J
  do J=1,16
     filenameBtermAU(J:J) = ' '
  enddo
  IF(i.LT.10)THEN
     write(filenameBtermAU,'(A,I1,A)')SPEC,I,'.dat'
  ELSEIF(i.LT.100)THEN
     write(filenameBtermAU,'(A,I2,A)')SPEC,I,'.dat'
  ELSEIF(i.LT.1000)THEN
     write(filenameBtermAU,'(A,I3,A)')SPEC,I,'.dat'
  ELSE
     CALL LSQUIT('nBterms too large in generateMCDFilename',-1)
  ENDIF
end subroutine generateMCDFilename

subroutine write_transition_density_matrix(nexci_max,lupri)
  implicit none
  !> Position of the highest lying excited state to be printet
  integer, intent(in)            :: nexci_max,lupri
  !
  integer :: j,luExciTransDensMat
  real(realk) :: ExEnergies(nexci_max) 
  character(len=19) :: Filename 
  logical :: OnMaster
  OnMaster=.TRUE.
  Do j=1,nexci_max
     ExEnergies(j) = REAL(rsp_eq_sol(j)%fld(1)%freq)
     WRITE(lupri,'(A,F16.8)')&
      & 'Writing transition density matrix to the state with excitation energy',ExEnergies(j)
     Filename = 'ExciTransDensMat'//Char(j/100+48)//Char(mod(j,100)/10+48)//Char(mod(mod(j,100),10)+48)
     WRITE(lupri,'(A,A)')'To Disk with Filename',Filename    
     luExciTransDensMat = -1
     CALL LSOPEN(luExciTransDensMat,Filename,'unknown','UNFORMATTED')
     rewind(luExciTransDensMat)    
     call mat_write_to_disk(luExciTransDensMat,rsp_eq_sol(j)%D,OnMaster)
     CALL LSCLOSE(luExciTransDensMat,'KEEP')
  enddo
end subroutine write_transition_density_matrix

 !> \brief main DipoleMomentMatrix driver
 !> \author Thomas Kjaergaard
 !> \date 2010-03
 !>
 !> Main driver to determine the full dipole moment matrix $\mu_{xy}$ where $\mu_{00}$ is the 
 !> dipole moment of the ground state while $\mu_{11}$ are the state dipole moments of the first 
 !> excited state, and $\mu_{ij}$ it the transition dipole moment between state i and j. 
 !> <br><br>
 !> The input structure for response functions is as follows:
 !> <br><br>
 !> **RESPONS
 !> <br>
 !> *DIPOLEMOMENTMATRIX
 !> <br> <br>
 !> Questions regarding this driver, the input structure etc. 
 !> may be addressed to tkjaergaard@chem.au.dk.
 subroutine DipoleMomentMatrix_driver(lupri,setting,decomp,solver,F,D,S)
  implicit none
  !> logical unit number of output 
  integer,intent(in)                :: lupri
  !> Info on molecule needed by solver and integral programs
  type(LSSETTING),intent(in),target :: setting
  !> Contains matrices from OAO decomposition of overlap matrix
  type(decompItem),intent(inout),target :: decomp
  !> RSPsolver settings
  type(RSPSOLVERinputitem),intent(inout),target  :: solver
  !> Unperturbed Fock matrix
  type(matrix), intent(in)          :: F
  !> Unperturbed density matrix
  type(matrix), intent(in)          :: D
  !> Unperturbed overlap matrix
  type(matrix), intent(in)          :: S
  !> local parameters
  !> Total number of response functions and transition moments requested
  integer                           :: n_rspfunc
  !> The position of the highest lying excited state requested
  integer                           :: nexci_max
  !> Stack of all response function and transition moments to be determined
  type(RspFunc), allocatable        :: RspFunc_stack(:)
  !> temporary response function structure
  type(RspFunc)                     :: tmpRspFunc
  complex(realk), allocatable       :: rsp_results(:,:)
  type(rsp_molcfg) :: molcfg
  real(realk),parameter :: hartree=27.21138386E0_realk
  logical     :: doprint
  real(realk) :: GSDipoleMoment(3)
  complex(realk) :: freqI,freqK
  real(realk),pointer :: DipoleMomentMatrix(:,:,:),excitationFreq(:)
  character(len=12) :: Dir(3)
  integer :: I,J,K,nrsp,KK,maxdim,l
  Dir(1) = 'X coordinate'
  Dir(2) = 'Y coordinate'
  Dir(3) = 'Z coordinate'
  write(lupri,*)'******************************************************'
  write(lupri,*)'**                                                  **'  
  write(lupri,*)'**  DipoleMomentMatrix Calculation                  **'
  write(lupri,*)'**                                                  **'  
  write(lupri,*)'******************************************************'
  
  !create config struct to be passed to rsp_contribs / rsp_equations
!/*point to natoms within structure*/
  molcfg = rsp_molcfg(S,setting%MOLECULE(1)%p%Natoms, &
  & decomp%lupri,decomp%luerr,setting,decomp,solver)
  molcfg%zeromat = 0*S
  !defined by .NEXCIT in input
  nexci_max = molcfg%decomp%cfg_rsp_nexcit
  call mem_alloc(DipoleMomentMatrix,nexci_max+1,nexci_max+1,3)
  doPrint = .TRUE.
! Step 1: Ground state Dipole MOment GSDipoleMoment(1:3)
  call Get_dipole_moment(molcfg,F,D,S,doPrint,GSDipoleMoment)
  do I = 1,3
     DipoleMomentMatrix(1,1,I) = GSDipoleMoment(I)
  enddo

  ! Step 2: Transition moments between the groundstate and the excited states
  
  !calc excitation energies and transition densities
  call transition_moment_density_matrix(molcfg,F,D,S,nexci_max)
  call write_transition_density_matrix(nexci_max,decomp%lupri)

!  write(lupri,*)'DipoleMomentMatrix: done transition_moment_density_matrix'
  maxdim = 3
  n_rspfunc = nexci_max 
  allocate(rsp_results(maxdim,n_rspfunc))
  rsp_results(:,:) = (0.0E0_realk, 0.0E0_realk)
!  WRITE(lupri,'(A)')'DipoleMomentMatrix: Transition moments between the groundstate and the excited states'
!  WRITE(lupri,'(2X,A5,2X,A8,2X,A8,9X,A1,11X,A1,11X,A1,7X)')&
!       &'STATE','Freq(au)','Freq(eV)','X','Y','Z'
  !calc transition moments to determine dipole allowed transitions
  call mem_alloc(excitationFreq,nexci_max)
  do i = 1,nexci_max
     excitationFreq(i) = REAL(rsp_eq_sol(i)%fld(1)%freq)
     call initialize_RspFunc(tmpRspFunc)
     tmpRspFunc%order = 2     
     tmpRspFunc%code(1) = 'EL  '
     tmpRspFunc%code(2) = 'EXCI'     
     tmpRspFunc%first_comp(1) = 1    
     tmpRspFunc%dims(1) = 3    
     tmpRspFunc%freq(2) = -rsp_eq_sol(i)%fld(1)%freq
     tmpRspFunc%freq(3) =  rsp_eq_sol(i)%fld(1)%freq
     tmpRspFunc%trans_moment = .TRUE. !single residue
     tmpRspFunc%trans_moment_order = 1 
     tmpRspFunc%ExNumber1 = i
     !calc transition moments 
     call linear_response(molcfg,F,D,S,tmpRspFunc, &
          &rsp_results(1:maxdim,i), maxdim)
!     WRITE(lupri,'(I5,2X,F10.6,2X,F10.6,2X,F10.7,2X,F10.7,2X,F10.7,6X,A4)') &
!          &i,REAL(rsp_eq_sol(i)%fld(1)%freq),REAL(rsp_eq_sol(i)%fld(1)%freq*hartree),&
!          &REAL(rsp_results(1,i)),REAL(rsp_results(2,i)),REAL(rsp_results(3,i))

     do J = 1,3
        DipoleMomentMatrix(1,1+I,J) = REAL(rsp_results(J,i))
        DipoleMomentMatrix(1+I,1,J) = REAL(rsp_results(J,i))
     enddo
  enddo
  deallocate(rsp_results)

  WRITE(lupri,'(A)')'DipoleMomentMatrix: Transition moments between the groundstate and the excited states'
  WRITE(lupri,'(A)')' '
  WRITE(lupri,'(2X,A5,2X,A8,2X,A8,9X,A1,11X,A1,11X,A1,7X)')&
       &'STATE','Freq(au)','Freq(eV)','X','Y','Z'
  do i = 1,nexci_max
     WRITE(lupri,'(I5,2X,F10.6,2X,F10.6,2X,F10.7,2X,F10.7,2X,F10.7,6X,A4)') &
          &i,excitationFreq(i),excitationFreq(i)*hartree,(DipoleMomentMatrix(1,1+I,J),J=1,3)
  enddo
  WRITE(lupri,'(A)')' '

  ! Step 3: Transition moments between the groundstate and the excited states

!  TrueDoubleResidue = .TRUE.
  nrsp = 0
  do i = 1,nexci_max
     do k = i,nexci_max
        nrsp = nrsp + 3
     enddo
  enddo
  n_rspfunc = nrsp
  allocate(RspFunc_stack(n_rspfunc))
  nrsp = 0
  do i = 1,nexci_max
     freqI = excitationFreq(i)
     do k = i,nexci_max
        freqK = excitationFreq(k)
        call set_doubleResidue_RspFunc(RspFunc_stack(nrsp+1),k,i,1,freqK,freqI,'EL  ')
        call set_doubleResidue_RspFunc(RspFunc_stack(nrsp+2),k,i,2,freqK,freqI,'EL  ')
        call set_doubleResidue_RspFunc(RspFunc_stack(nrsp+3),k,i,3,freqK,freqI,'EL  ')
        nrsp = nrsp + 3
        DipoleMomentMatrix(1+I,1+K,1) = 0.0E0_realk
        DipoleMomentMatrix(1+I,1+K,2) = 0.0E0_realk
        DipoleMomentMatrix(1+I,1+K,3) = 0.0E0_realk
        DipoleMomentMatrix(1+K,1+I,1) = 0.0E0_realk
        DipoleMomentMatrix(1+K,1+I,2) = 0.0E0_realk
        DipoleMomentMatrix(1+K,1+I,3) = 0.0E0_realk
     enddo
  enddo
  n_rspfunc = nrsp
  call rspfunc_or_transmoment(molcfg,F,D,S,RspFunc_stack,n_rspfunc)
  nrsp = 0
  do i = 1,nexci_max
     freqI = excitationFreq(i)
     do k = i,nexci_max
        freqK = excitationFreq(k)
        DipoleMomentMatrix(1+I,1+K,1) = REAL(RspFunc_stack(nrsp+1)%result(1))
        DipoleMomentMatrix(1+I,1+K,2) = REAL(RspFunc_stack(nrsp+2)%result(1))
        DipoleMomentMatrix(1+I,1+K,3) = REAL(RspFunc_stack(nrsp+3)%result(1))
        DipoleMomentMatrix(1+K,1+I,1) = REAL(RspFunc_stack(nrsp+1)%result(1))
        DipoleMomentMatrix(1+K,1+I,2) = REAL(RspFunc_stack(nrsp+2)%result(1))
        DipoleMomentMatrix(1+K,1+I,3) = REAL(RspFunc_stack(nrsp+3)%result(1))
        nrsp = nrsp + 3
     enddo
  enddo
  deallocate(RspFunc_stack)

  WRITE(lupri,'(A)')'DipoleMomentMatrix Summary'
  call PrintDipoleMoment(GSDipoleMoment,lupri)

  WRITE(lupri,'(A)')'DipoleMomentMatrix: Transition moments between the groundstate and the excited states'
  WRITE(lupri,'(A)')' '
  WRITE(lupri,'(2X,A5,2X,A8,2X,A8,9X,A1,11X,A1,11X,A1,7X)')&
       &'STATE','Freq(au)','Freq(eV)','X','Y','Z'
  do i = 1,nexci_max
     WRITE(lupri,'(I5,2X,F10.6,2X,F10.6,2X,F10.7,2X,F10.7,2X,F10.7,6X,A4)') &
          &i,excitationFreq(i),excitationFreq(i)*hartree,(DipoleMomentMatrix(1,1+I,J),J=1,3)
  enddo
  WRITE(lupri,'(A)')' '
  WRITE(lupri,'(A)')'The transition dipole moments between state X and Y for operator XDIPLEN: <X | A - <A> | Y>'
  call output(DipoleMomentMatrix(2:nexci_max+1,2:nexci_max+1,1),1,nexci_max,1,nexci_max,nexci_max,nexci_max,1,lupri)
  WRITE(lupri,'(A)')' '

  WRITE(lupri,'(A)')'The transition dipole moments between state X and Y for operator YDIPLEN: <X | A - <A> | Y>'
  call output(DipoleMomentMatrix(2:nexci_max+1,2:nexci_max+1,2),1,nexci_max,1,nexci_max,nexci_max,nexci_max,1,lupri)
  WRITE(lupri,'(A)')' '

  WRITE(lupri,'(A)')'The transition dipole moments between state X and Y for operator ZDIPLEN: <X | A - <A> | Y>'
  call output(DipoleMomentMatrix(2:nexci_max+1,2:nexci_max+1,3),1,nexci_max,1,nexci_max,nexci_max,nexci_max,1,lupri)
  WRITE(lupri,'(A)')' '
  WRITE(lupri,'(A)')'The Full Dipole Moment Matrix'
  WRITE(lupri,'(A)')' '
  WRITE(lupri,'(A)')'Column 0 indicate the Ground State. Element (0,0) is therefore the ground state dipole moment'
  WRITE(lupri,'(A)')'Elements (0,X) is therefore the transition state dipole moment of the X. excited state'
  WRITE(lupri,'(A)')'Elements (X,Y) is the transition dipole moment between state X and Y'
  WRITE(lupri,'(A)')'Note (X,X) X!=0, denotes <X|A|X> - <A> = <X|A|X> - <0|A|0>'

  do J=1,3
     WRITE(lupri,'(A,A12)')'Dipole Moment Matrix for ',Dir(J)
     WRITE(lupri,'(A)')' '

!     call output(DipoleMomentMatrix(:,:,J),1,1+nexci_max,1,1+nexci_max,1+nexci_max,1+nexci_max,1,lupri)

     DO K = 1,nexci_max+1,4
        IF(K.LE.nexci_max+1-3)THEN
           IF(K.EQ.1)THEN
              WRITE(lupri,'(A,I4,A,I4,A,I4,A,I4)')&
                   &'               Column',K-1,' |   Column',K,'     Column',K+1,'     Column',K+2
              WRITE(lupri,'(A,I4,A2,F15.8,A,F13.8,2F15.8)')&
                      &'    ',0,'  ',DipoleMomentMatrix(1,K,J),' |',(DipoleMomentMatrix(1,K+KK,J),KK=1,3)
              WRITE(lupri,'(A)')'       ---------------------------------------------------------------'
              DO L = 2,nexci_max+1
                 WRITE(lupri,'(A,I4,A2,F15.8,A,F13.8,2F15.8)')&
                      &'    ',L-1,'  ',DipoleMomentMatrix(L,K,J),' |',(DipoleMomentMatrix(L,K+KK,J),KK=1,3)
              ENDDO
           ELSE
              WRITE(lupri,'(A,I4,A,I4,A,I4,A,I4)')&
                   &'               Column',K-1,'     Column',K,'     Column',K+1,'     Column',K+2
              WRITE(lupri,'(A,4F15.8)')'       0  ',(DipoleMomentMatrix(1,K+KK,J),KK=0,3)
              WRITE(lupri,'(A)')'       ---------------------------------------------------------------'
              DO L = 2,nexci_max+1
                 WRITE(lupri,'(A,I4,A2,4F15.8)')&
                      &'    ',L-1,'  ',(DipoleMomentMatrix(L,K+KK,J),KK=0,3)
              ENDDO
           ENDIF
        ELSEIF(K.LE.nexci_max+1-2)THEN
           IF(K.EQ.1)THEN
              WRITE(lupri,'(A,I4,A,I4,A,I4)')&
                   &'               Column',K-1,' |   Column',K,'     Column',K+1
!              WRITE(lupri,'(A,3F15.8)')'       0  ',(DipoleMomentMatrix(1,K+KK,J),KK=0,2)
              WRITE(lupri,'(A,I4,A2,F15.8,A,F13.8,F15.8)')&
                   &'    ',0,'  ',DipoleMomentMatrix(1,K,J),' |',(DipoleMomentMatrix(1,K+KK,J),KK=1,2)
              WRITE(lupri,'(A)')'       ------------------------------------------------'
              DO L = 2,nexci_max+1
                 WRITE(lupri,'(A,I4,A2,F15.8,A,F13.8,F15.8)')&
                      &'    ',L-1,'  ',DipoleMomentMatrix(L,K,J),' |',(DipoleMomentMatrix(L,K+KK,J),KK=1,2)
              ENDDO
           ELSE
              WRITE(lupri,'(A,I4,A,I4,A,I4)')&
                   &'               Column',K-1,'     Column',K,'     Column',K+1
              WRITE(lupri,'(A,3F15.8)')'       0  ',(DipoleMomentMatrix(1,K+KK,J),KK=0,2)
              WRITE(lupri,'(A)')'       ------------------------------------------------'
              DO L = 2,nexci_max+1
                 WRITE(lupri,'(A,I4,A2,3F15.8)')&
                      &'    ',L-1,'  ',(DipoleMomentMatrix(L,K+KK,J),KK=0,2)
              ENDDO
           ENDIF
        ELSEIF(K.LE.nexci_max+1-1)THEN
           IF(K.EQ.1)THEN
              WRITE(lupri,'(A,I4,A,I4)')&
                   &'               Column',K-1,' |   Column',K
!              WRITE(lupri,'(A,2F15.8)')'       0  ',(DipoleMomentMatrix(1,K+KK,J),KK=0,1)
              WRITE(lupri,'(A,I4,A2,F15.8,A,F13.8)')&
                   &'    ',0,'  ',DipoleMomentMatrix(1,K,J),' |',DipoleMomentMatrix(1,K+1,J)
              WRITE(lupri,'(A)')'       ---------------------------------'
              DO L = 2,nexci_max+1
                 WRITE(lupri,'(A,I4,A2,F15.8,A,F13.8)')&
                      &'    ',L-1,'  ',DipoleMomentMatrix(L,K,J),' |',DipoleMomentMatrix(L,K+1,J)
              ENDDO
           ELSE
              WRITE(lupri,'(A,I4,A,I4)')&
                   &'               Column',K-1,'     Column',K
              WRITE(lupri,'(A,2F15.8)')'       0  ',(DipoleMomentMatrix(1,K+KK,J),KK=0,1)
              WRITE(lupri,'(A)')'       ---------------------------------'
              DO L = 2,nexci_max+1
                 WRITE(lupri,'(A,I4,A2,2F15.8)')&
                      &'    ',L-1,'  ',(DipoleMomentMatrix(L,K+KK,J),KK=0,1)
              ENDDO
           ENDIF
        ELSEIF(K.LE.nexci_max+1)THEN
           IF(K.EQ.1)THEN
              WRITE(lupri,'(A,I4)')&
                   &'               Column',K-1
              WRITE(lupri,'(A,1F15.8)')'       0  ',DipoleMomentMatrix(1,K,J)
              WRITE(lupri,'(A)')'       ------------------'
              DO L = 2,nexci_max+1
                 WRITE(lupri,'(A,I4,A2,1F15.8)')&
                      &'    ',L-1,'  ',DipoleMomentMatrix(L,K,J)
              ENDDO
           ELSE
              WRITE(lupri,'(A,I4)')&
                   &'               Column',K-1
              WRITE(lupri,'(A,1F15.8)')'       0  ',DipoleMomentMatrix(1,K,J)
              WRITE(lupri,'(A)')'       ------------------'
              DO L = 2,nexci_max+1
                 WRITE(lupri,'(A,I4,A2,1F15.8)')&
                      &'    ',L-1,'  ',DipoleMomentMatrix(L,K,J)
              ENDDO
           ENDIF
        ENDIF
        WRITE(lupri,'(A)')' '
     ENDDO
  enddo
  WRITE(lupri,'(A)')'DipoleMomentMatrix Summary done'
  call mem_dealloc(DipoleMomentMatrix)
  call mem_dealloc(excitationFreq)

end subroutine DipoleMomentMatrix_driver

#else
CONTAINS
subroutine response_wrapper_dummy_routine()
implicit none

end subroutine response_wrapper_dummy_routine
#endif
end module response_wrapper_module
