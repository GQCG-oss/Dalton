!> @filed
!> Contains response_driver_module.

!> Calculates linear, quadratic, and cubic response functions and their residues.
!> JCP, 129, 214108 (2008).
!> \author Kasper Kristensen
!> \date 2010-01
module response_driver_module

#ifdef VAR_RSP
  use decompMod, only: DecompItem
  !use scf_config, only: AOR_stacksize, cfg_rsp_imfreqs_specified, cfg_rsp_complex
  use precision
  use TYPEDEFTYPE, only: LSSETTING
  use matrix_module
  use lsdalton_matrix_defop
  use lsdalton_rsp_contribs
  use lsdalton_rsp_equations
  use RSPsolver, only: rsp_init, rsp_solver, rsp_molcfg  
  use RSPsymsolver, only: rsp_sym_init, rsp_sym_solver
  use matrix_operations, only: mat_write_to_disk, mat_read_from_disk
  use memory_handling

  implicit none
  !> Maximum order of response functions, current 4 (cubic response)
  integer, parameter :: max_order=4

  !> Type for response functions and transition moments
  type RspFunc

     !> order of response function (2 for linear, 3 for quadratic, 4 for cubic)
     integer                 :: order
     !> code(i): Operator code (see rsp_contribs.F90) for operator (i)
     character*4             :: code(max_order)
     !> first_comp(i): First field component to be calculated for operator (i)
     integer                 :: first_comp(max_order)
     !> dims(i): Total number of field components to be calculated for operator (i)
     integer                 :: dims(max_order)
     !> result: Vector containing response function/transition moment result(s)
     complex(realk),pointer  :: result(:)
     !> freq(i): Frequency for operator (i) (may equal plus/minus the excitation energy
     !> for transition moments)
     complex(realk)          :: freq(max_order)
     !> Is the requested property a transition moment (.true.) or not (.false.)
     logical                 :: trans_moment
     !> Order of the transition moment (redundant for response function), i.e. the
     !> number of poles of the original response function (1/2 if the transition moment
     !> originates from a single/double residue of the response function).
     integer :: trans_moment_order
     !> Position of the excited state at the (first) pole (redundant for response function).
     !> Thus each type RspFunc only contains the transition moment for ONE excited state.
     integer :: ExNumber1
     !> Position of the excited state at the second pole (redundant for response function).
     !> It is only used for transition moments originating from double residues 
     !> of response function.
     integer :: ExNumber2

  end type RspFunc

  !  public response_driver, RspFunc, initialize_rspfunc, transition_moment_density_matrix, &
  !       &linear_response, rspfunc_or_transmoment, print_response_function, &
  !       &quadratic_response_2nplus1, cubic_response, Get_first_order_property, &
  !       &quadratic_response_nplus1
  public RspFunc, initialize_rspfunc, transition_moment_density_matrix, &
       &linear_response, rspfunc_or_transmoment, &
       &quadratic_response_2nplus1, cubic_response, Get_first_order_property, &
       &quadratic_response_nplus1, print_response_func
  private

Contains

  !> Main response driver. Sets up response input, prints input, 
  !> calculates response function,
  !> and prints out the response results.
  !> <br><br>
  !> It should work for any combination of operators, i.e.
  !> we may calculate <<A;B>> and <<A;B,C>> and their residues for
  !> any choice of A,B, and C - also if the basis functions
  !> are perturbation dependent (e.g. London orbitals and geometric dependencies).
  !> <br><br>
  !> The input structure for response functions is as follows:
  !> <br><br>
  !> $RESPONS
  !> <br>
  !> *AORESPONSE
  !> <br>
  !> .RESPONSE-FUNCTION
  !> <br>
  !> (order) (operators)
  !> <br>
  !> .COMPONENTS   (optional)
  !> <br>
  !> (order) (components for the different operators)
  !> <br>
  !> .FREQ (optional)
  !> <br>
  !> (number of real frequencies) (which frequencies)
  !> <br>
  !> .IMFREQ (optional)
  !> <br>
  !> (number of imaginary frequencies) (which frequencies)
  !> <br>
  !> *END AORESPONSE
  !> <br>
  !> $END RESPONS
  !> <br> <br>

  !> This input calculates a set of response functions for the components and frequencies
  !> specified.
  !> If .COMPONENTS is not present, ALL field components for ALL operators are determined.
  !> <br>
  !> NOTE: Only ONE set of frequencies may be specified - e.g. if linear response functions
  !> <<A;B>> for the frequencies 0.1 and 0.2 are requested two .RESPONSE-FUNCTION 
  !> inputs are required!
  !> If no frequencies are specified all frequencies are set to zero.
  !> See pert_table in rsp_contribs.F90 for a list of operator abbreviations.

  !> <br> <br>

  !> Example: <<EL;EL,MAG>> for components x (1:1) for the first operator (EL);
  !> x and y (1:2) for the second operator (EL); x,y, and z (1:3) for the third operator (MAG).
  !> Real frequencies for the second and third operators are 0.1 and 0.2, respectively
  !> (meaning that the frequency of the first operator is set to -0.1-0.2 = -0.3).
  !> Imaginary frequencies are set to zero (default). Input:

  !> <br> <br>

  !> $RESPONS
  !> <br>
  !> *AORESPONSE
  !> <br>
  !> .RESPONSE-FUNCTION
  !> <br>
  !> 3 EL EL MAG
  !> <br>
  !> .COMPONENTS
  !> <br>
  !> 3 1:1 1:2 1:3
  !> <br>
  !> .FREQ
  !> <br>
  !> 2 0.1 0.2
  !> <br>
  !> *END AORESPONSE
  !> <br>
  !> $END RESPONS
  !> <br> <br>

  !> The input structure for transition moments is very similar.
  !> Consider the linear transition moment of the form <0|A|n>. These are identifed from
  !> residues of the linear response function, see Eq. (3.12) i J. Chem. Phys. 82, 3242 (1985).

  !> <br><br>

  !> lim (omega_b - omega_n) <<A;B>> = <0|A|n> * <n|B|0>

  !> <br><br>

  !> In the input structure we use the notation:

  !> <br><br>

  !> <0|A|n> * <n|B|0> = <<A;EXCI>> * <<B;EXCI>>

  !> <br><br>

  !> To exemplify the input structure for calculating transition moments originating from
  !> residues of the linear response function, 
  !> consider the linear transition moment <0|EL|n> = <<EL;EXCI>> 
  !> for the x-component of the electric dipole
  !> and for the first 5 excited states (n=1,5).
  !> The input is then:

  !> <br><br>
  !> $RESPONS
  !> <br>
  !> *AORESPONSE
  !> <br>
  !> .TRANSITION-MOMENT

  !> <br>

  !> 2 EL EXCI

  !> <br>

  !> .COMPONENTS

  !> <br>

  !> 1 1:1

  !> <br>

  !> .NEXCIT

  !> <br>

  !> 5

  !> <br>
  !> *END AORESPONSE
  !> <br>
  !> $END RESPONS
  !> <br> <br>

  !> Note that only ONE component is specified because the "components"
  !> for EXCI is specified by the number of excited states in .NEXCIT.
  !> Also note that transition moments for different excited states are stored
  !> as seperate entries in the response function stack.
  !> <br><br>
  !> For single residues of quadratic response functions (quadratic transition moments)
  !> we choose the residue convention in Eq. (3.18) in
  !> J. Chem. Phys. 82, 3242 (1985) -
  !> HOWEVER with a different sign on omega_b!!!
  !> <br><br>
  !> lim(omega_c - omega_m) <<A;B,C>>_{omega_b, omega_c} = <<A;B,EXCI>> * <<C;EXCI>>
  !> <br><br>
  !> where <<C;EXCI>> = <m|C|0> (see description of linear transition moment above)
  !> and <<A;B,EXCI>> equals the remaining 'quadratic part' of Eq. (3.18).
  !> <br><br>
  !> To calculate <<A;B,EXCI>> where A is an electric dipole operator,
  !> B is a magnetic dipole operator (using London orbitals), 
  !> and omega_b=0 (and thus omega_a = -omega_c = -omega_m)
  !> for the first four excited states
  !> we use the following input:
  !> <br><br>
  !> $RESPONS
  !> <br>
  !> *AORESPONSE
  !> <br>
  !> .TRANSITION-MOMENT

  !> <br>

  !> 3 EL MAG EXCI

  !> <br>

  !> .NEXCIT

  !> <br>

  !> 4

  !> <br>

  !> .FREQ

  !> <br>

  !> 1    0E0_realk
  !> <br>
  !> *END AORESPONSE
  !> <br>
  !> $END RESPONS
  !> <br> <br>

  !> From the example it is seen that for the residue the .FREQ keyword sets set
  !> frequency for the omega_b frequency.
  !> <br> <br>

  !> For double residues
  !> we choose the residue convention in Eq. (3.21) in
  !> J. Chem. Phys. 82, 3242 (1985).
  !> <br><br>
  !> lim(omega_c - omega_m)(omega_b+omega_q) <<A;B,C>>_{omega_b, omega_c} =
  !> <br>
  !> <0|B|q> * (- <q| (A - <0|A|0>)|m> ) <m|C|0> =
  !> <br>
  !> <<B;EXCIq>> * <<A;EXCIq,EXCIm>> * <<C;EXCIm>>
  !> <br><br>
  !> where the minus sign is formally absorbed into <<A;EXCIq,EXCIm>>.
  !> To calculate <<A;EXCIq,EXCIm>> where A is an electric dipole
  !> where m runs from state 1 to 5 and q runs from state 1 to 3 we
  !> use the following input:
  !> <br><br>
  !> $RESPONS
  !> <br>
  !> *AORESPONSE
  !> <br>

  !> .TRANSITION-MOMENT

  !> <br>

  !> 3 EL EXCI EXCI

  !> <br>

  !> .NEXCIT

  !> <br>

  !> 5

  !> <br>

  !> .NEXCIT2

  !> <br>

  !> 3

  !> <br>
  !> *END AORESPONSE
  !> <br>
  !> $END RESPONS
  !> <br> <br>
  !> Finally, I know the input is semi-complicated but it is very difficult to
  !> make simpler when it needs to be so general. Anyways, most specific properties will
  !> soon be accesible through special and much simpler inputs.
  !> Questions regarding this driver, the input structure etc. may be addressed to kasperk@chem.au.dk.
  !> \author Kasper Kristensen
  !> \date 2010-01
  !  subroutine response_driver(setting,decomp,solver,F,D,S,AOR_input)
  !    use TYPEDEF,       only: LSSETTING
  !    use molecule_type, only: MOLECULE_PT, MOLECULEINFO
  !    use response_wrapper_type_module, only: RSPSOLVERinputitem
  !    !> Contains matrices from OAO decomposition of overlap matrix
  !    type(LSSETTING),intent(in),target     :: setting
  !    !> Contains matrices from OAO decomposition of overlap matrix
  !    type(decompItem),intent(inout),target :: decomp
  !    !> RSPsolver settings
  !    type(RSPSOLVERinputitem),intent(inout),target  :: solver
  !    !> Unperturbed Fock matrix
  !    type(matrix), intent(in)          :: F
  !    !> Unperturbed density matrix
  !    type(matrix), intent(in)          :: D
  !    !> Unperturbed overlap matrix
  !    type(matrix), intent(in)          :: S
  !    !> Input from *AORESPONSE to *END AORESPONSE in the dalton input file - but
  !    !> stripped from (potential) comment lines
  !    character(len=40), intent(inout)  :: AOR_input(AOR_stacksize)
  !    !> Info on molecule needed by solver and integral programs
  !    type(rsp_molcfg)                  :: molcfg
  !    !> Total number of response functions and transition moments requested
  !    integer                           :: n_rspfunc
  !    !> The position of the highest lying excited state requested
  !    integer ::  nexci_max
  !    !> Stack of all response function and transition moments to be determined
  !    type(RspFunc), allocatable        :: RspFunc_stack(:)
  !    !> Do any excited states need to be determined
  !    logical                           :: ExEnergies_requested=.false.
  !    integer, target :: natoms
  !    integer :: i
  !    call determine_number_of_rspfuncs_and_exenergies(n_rspfunc,AOR_input, &
  !         & nexci_max, ExEnergies_requested)
  !
  !    allocate(RspFunc_stack(n_rspfunc))
  !
  !    !create config struct to be passed to rsp_contribs / rsp_equations
  !    call init_rsp_molcfg(molcfg,S,setting%MOLECULE(1)%p%Natoms, &
  !                       & decomp%lupri,decomp%luerr,setting,decomp,solver)
  !
  !    ! Set up response input to prepare response calculation(s)
  !    call response_input(molcfg,F,D,S,AOR_input,RspFunc_stack,n_rspfunc, &
  !         & nexci_max, ExEnergies_requested)
  !
  !
  !    ! Print input for response calculation, this print should ease identification of input errors.
  !    call print_response_function(RspFunc_stack,n_RspFunc,.false.,molcfg%lupri)
  !
  !    ! Calculate all response functions.
  !    call rspfunc_or_transmoment(molcfg,F,D,S,RspFunc_stack,n_rspfunc)
  !
  !    ! Print response results
  !    call print_response_function(RspFunc_stack, n_RspFunc,.true.,molcfg%lupri)
  !
  !    deallocate(RspFunc_stack)
  !
  !    call rsp_eq_sol_empty()
  !
  !
  !  end subroutine response_driver


  !> \brief Sets up the response function stack, which contains 
  !> information about all response functions to be determined.
  !> \author Kasper Kristensen
  !> \date 2010-01
  !  subroutine response_input(molcfg,F,D,S,AOR_input,RspFunc_stack,n_rspfunc, &
  !       & nexci_max,ExEnergies_requested)
  !    !> Info on molecule needed by solver and integral programs
  !    type(rsp_molcfg), intent(inout)   :: molcfg
  !    !> Unperturbed Fock matrix
  !    type(matrix), intent(in)          :: F
  !    !> Unperturbed density matrix
  !    type(matrix), intent(in)          :: D
  !    !> Unperturbed overlap matrix
  !    type(matrix), intent(in)          :: S
  !    !> Input from *AORESPONSE to *END AORESPONSE in the dalton input file - but
  !    !> stripped from (potential) comment lines
  !    character(len=40), intent(in)     :: AOR_input(AOR_stacksize)
  !    !> Stack of all response function and transition moments to be determined
  !    type(RspFunc), intent(inout)  :: RspFunc_stack(n_rspfunc)
  !    !> Total number of response functions and transition moments requested
  !    integer, intent(in)               :: n_rspfunc
  !    !> The position of the highest lying excited state requested
  !    integer, intent(in)               :: nexci_max
  !    !> Do any excited states need to be determined?
  !    logical, intent(in)               :: ExEnergies_requested
  !    !> Local parameters.
  !    type(RspFunc)                    :: MyRspFunc
  !    character(len=40)                :: word , component_string, freq_string
  !    integer                          :: i,j, k, nexci, nexci2, order_temp, &
  !         & rsp_counter, num_aux
  !    real(realk)                      :: Ifreq(max_order), Rfreq(max_order)
  !    integer, pointer                 :: lupri
  !
  !    lupri => molcfg%lupri
  !    rsp_counter=0; num_aux=0
  !
  !    ! If transition moment properties are requested we must solve the eigenvalue
  !    ! problem for nexci_max excitation energies.
  !    if(ExEnergies_requested) then
  !       call transition_moment_density_matrix(molcfg,F,D,S,nexci_max)
  !    endif
  !
  !    ! Loop over dalton input and store information about response functions 
  !    ! and transition moments in RspFunc_stack. 
  !    ! For clarity, for each input we store all response information in MyRspFunc
  !    ! and copy the information to the RspFunc_stack at the end of the loop.
  !    do i=1,AOR_stacksize
  !
  !       select case(AOR_input(i))
  !
  !          ! ****************************************************************************
  !          ! *                          RESPONSE FUNCTION INPUT                         *
  !          ! ****************************************************************************
  !
  !       CASE('.RESPONSE-FUNCTION')
  !
  !          ! A response function has been requested in the input, i.e. NOT a transition moment.
  !          MyRspFunc%trans_moment = .false.
  !
  !          ! Set response function parameters to the default values.
  !          call initialize_RspFunc(MyRspFunc)
  !
  !
  !
  !          ! ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
  !          !         Specifying the order of the response function and which operators        !
  !          !                  (not the specific components of the operators)                  !
  !          ! ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
  !
  !
  !          ! Read in the number of operators (order) in response function 
  !          ! (2 for linear, 3 for quadratic, 4 for cubic)
  !          word=AOR_input(i+1)
  !          READ (word(1:1), *) MyRspFunc%order
  !
  !          ! Read in operators in code
  !          ! E.g. for the linear response function <<EL; MAG>>:
  !          ! code(1) = 'EL  '
  !          ! code(2) = 'MAG '
  !          ! code(i) = 'NOOP'    (i=3, max_order)
  !          ! (num_aux is only relevant for properties using the HERMIT integral label, see
  !          !  pert_table in rsp_contribs.F90)
  !          call read_words_from_string(word(2:40),MyRspFunc%order, MyRspFunc%code,num_aux)
  !
  !
  !          ! ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
  !          !  Specifying the specific components for each operators in the response function  !
  !          ! ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
  !
  !          ! The (optional) .COMPONENTS keyword must be specied in 
  !          ! position 2 after .RESPONSE-FUNCTION.
  !          word = AOR_input(i+2)
  !          if(word(1:11) == '.COMPONENTS') then ! Components specified
  !             component_string = AOR_input(i+3)
  !
  !             ! Clean up string of components for easy read in.
  !             call clean_up_string(component_string,':', ' ')
  !
  !             ! Read in order and components of operators (see input in response_driver_module)
  !             READ (component_string, *) order_temp, &
  !                  & (MyRspFunc%first_comp(j), MyRspFunc%dims(j), j=1,MyRspFunc%order)
  !
  !             ! Now MyRspFunc%dims(j) equals that last component requested.
  !             ! The actual dimension equals the last component + 1 minus the first comp.
  !             do j=1, MyRspFunc%order
  !                MyRspFunc%dims(j) = MyRspFunc%dims(j) + 1 - MyRspFunc%first_comp(j)
  !             enddo
  !
  !             if(order_temp /= MyRspFunc%order) CALL lsQUIT('response_input: The order of the&
  !                  & response function must equal the order specified in .COMPONENTS',lupri)
  !
  !          else ! Default: Calculate all components (e.g. x,y, and z for electric dipole)
  !
  !             MyRspFunc%dims(1:MyRspFunc%order) = pert_shape(molcfg,MyRspFunc%code(1:MyRspFunc%order))
  !
  !          endif ! .COMPONENTS
  !
  !
  !
  !          ! ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
  !          !                     Specifying the frequencies for response functions            !
  !          ! ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
  !
  !
  !          Rfreq(:) = 0E0_realk
  !          Ifreq(:) = 0E0_realk
  !          do k=1,2
  !
  !             ! REAL FREQUENCIES:
  !             ! Real frequencies have been specified using .FREQ 
  !             ! in the second or fourth position after .RESPONSE-FUNCTION.
  !
  !             ! Check that we do not address non-defined position in string vector.
  !             if(size(AOR_input)>=i+2*k) then
  !                freq_string =AOR_input(i+2*k) 
  !                ! In case no frequencies are specified but instead another 
  !                ! response function or transition moment -> exit loop.
  !                ! In that case all frequencies are zero by the initialization above.
  !                if(freq_string(1:18) == '.RESPONSE-FUNCTION' .or. &
  !                     freq_string(1:18) == '.TRANSITION-MOMENT') exit
  !
  !                if(freq_string(1:5) == '.FREQ') then
  !                   ! Read in real frequencies for all operators except the first,
  !                   ! e.g. for << GEO; EL, MAG>> we read in frequencies for EL and MAG.
  !                   ! (See input structure in response_driver_module)
  !                   READ (AOR_input(i+2*k+1), *) order_temp, &
  !                        & (Rfreq(j), j=2,MyRspFunc%order)
  !                   ! w1 = -w2 -w3 -...
  !                   Rfreq(1) = -sum(Rfreq(2:MyRspFunc%order))
  !
  !                   ! Check correct real frequency input
  !                   if(order_temp /= MyRspFunc%order-1)  then
  !                      CALL lsQUIT('response_input: The number of real frequencies specified &
  !                           &for a response function must equal the order &
  !                           &of the response function minus one.',lupri)
  !                   endif !QUIT
  !
  !                endif !.FREQ
  !             end if
  !
  !
  !
  !             ! IMAGINARY FREQUENCIES:
  !             ! Imaginary frequencies have been specified using .IMFREQ 
  !             ! in the fourth or sixth position after .RESPONSE-FUNCTION.
  !             ! Check that we do not address non-defined position in string vector.
  !             if(size(AOR_input)>=i+2*k+2) then
  !                freq_string =AOR_input(i+2*k+2)
  !                if(freq_string(1:7) == '.IMFREQ') then
  !                   cfg_rsp_imfreqs_specified = .true.
  !                   ! Same as above for .FREQ
  !                   READ (AOR_input(i+2*k+3), *) order_temp, &
  !                        & (Ifreq(j), j=2,MyRspFunc%order)
  !                   ! KK HACK: FIX ME, set gamma = Ifreq(2)
  !                   molcfg%solver%rsp_gamma = Ifreq(2)
  !                   Ifreq(1) = -sum(Ifreq(2:MyRspFunc%order))
  !
  !                   ! Check correct imaginary frequency input
  !                   if(order_temp /= MyRspFunc%order-1)  then
  !                      CALL lsQUIT('response_input: The number of imaginary frequencies specified &
  !                           &for a response function must equal the order &
  !                           &of the response function minus one.',lupri)
  !                   endif !QUIT
  !
  !                endif !.IMFREQ
  !             end if
  !
  !          enddo
  !
  !          ! Copy real and imaginary frequencies to MyRspFunc.
  !          do k=1,MyRspFunc%order
  !             MyRspFunc%freq(k) = Rfreq(k)*(1E0_realk,0E0_realk) + Ifreq(k)*(0E0_realk,1E0_realk)
  !          enddo
  !
  !
  !
  !          ! Finally, all information about the response function is now stored in
  !          ! MyRspFunc and can be copied to the stack.
  !          rsp_counter = rsp_counter+1
  !          RspFunc_stack(rsp_counter) = MyRspFunc
  !
  !
  !
  !
  !
  !
  !          ! ****************************************************************************
  !          ! *                          TRANSITION MOMENT INPUT                         *
  !          ! ****************************************************************************
  !
  !       Case('.TRANSITION-MOMENT')
  !
  !          ! Much of the code is the same as for .RESPONSE-FUNCTION, 
  !          ! but there are some important differences so we keep it separate here.
  !
  !          call initialize_RspFunc(MyRspFunc)
  !
  !          MyRspFunc%trans_moment=.true.
  !
  !          ! First we must determine how many excited states are present
  !          ! (it will usually be the same as nexci_max, but not necessarilly).
  !          ! This is (maybe) written under .NEXCIT. Otherwise, set value to nexci_max.
  !          nexci=0
  !          do j=1,AOR_stacksize-i
  !             word = AOR_input(i+j)
  !
  !             ! In case no excitations are specified but instead another 
  !             ! response function or transition moment -> exit loop.
  !             if(word(1:18) == '.RESPONSE-FUNCTION' .or. &
  !                  word(1:18) == '.TRANSITION-MOMENT') exit
  !
  !             ! nexci specified in the next line
  !             if(word(1:8) == '.NEXCIT ')  READ (AOR_input(i+j+1), '(I10)') nexci
  !
  !          enddo
  !
  !          ! If nexci has not been specified for this transition moment, set nexci equal nexci_max
  !          if(nexci == 0) nexci=nexci_max
  !
  !
  !          ! ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
  !          !         Specifying the order of the transition moment and which operator(s)      !
  !          !                       (the same for all excited states)                          !
  !          ! ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
  !
  !
  !          ! Read in number of operators (order) in transition moment
  !          word=AOR_input(i+1)
  !          READ (word(1:1), '(I10)') MyRspFunc%order
  !          ! Read in operators in code (see response_driver_module)
  !          ! E.g. for the electric dipole transition moment <<EL;EXCI>> = <0|EL|n>:
  !          ! code(1) = 'EL  '
  !          ! code(2) = 'EXCI'
  !          ! code(i) = 'NOOP'    (i=3, max_order)
  !          call read_words_from_string(word(2:40),MyRspFunc%order, MyRspFunc%code,num_aux)
  !
  !
  !          ! ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
  !          !   Specifying the specific components for each operator in the transition moment  !
  !          ! ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
  !
  !          word = AOR_input(i+2)
  !          if(word(1:11) == '.COMPONENTS') then ! Components specified in the next line
  !             component_string = AOR_input(i+3)
  !
  !             ! Clean up string of components for easy read in.
  !             call clean_up_string(component_string,':', ' ')
  !
  !             ! Read in order and components of operators (see input in response_driver_module)
  !             ! This is only for the operators which are NOT 'EXCI'.
  !             ! For 'EXCI' operator, first_comp=1 and dims=1.
  !             READ (component_string, *) order_temp, &
  !                  & (MyRspFunc%first_comp(j), MyRspFunc%dims(j), j=1,order_temp)
  !
  !             ! Now MyRspFunc%dims(j) equals that last component requsted.
  !             ! The actual dimension equals the last component + 1 minus the first comp.
  !             do j=1, order_temp
  !                MyRspFunc%dims(j) = MyRspFunc%dims(j) + 1 - MyRspFunc%first_comp(j)
  !             enddo
  !
  !
  !          else  ! Default: Calculate all components (e.g. x,y, and z for electric dipole)
  !
  !             MyRspFunc%dims(1:MyRspFunc%order) = pert_shape(molcfg,MyRspFunc%code(1:MyRspFunc%order))
  !
  !          endif
  !
  !          ! ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
  !          !     Specifying the frequencies (excitation energies) for transition moments      !
  !          ! ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
  !
  !          ! Loop over excitaiton energies.
  !          ! If double transition moment, this j-loop is 
  !          ! over the excited states defined by .NEXCIT,
  !          ! whereas the loop over excited states defined by .NEXCIT2 occurs in the k-loop below.
  !
  !          do j=1,nexci ! loop over excited states
  !
  !             ! In all cases the last operator (position=order) is always 'EXCI' and 
  !             ! its frequency equals the excication energy.
  !             MyRspFunc%freq(MyRspFunc%order) = rsp_eq_sol(j)%fld(1)%freq
  !             ! Store the position of the excitation.
  !             MyRspFunc%ExNumber1 = j
  !
  !             ! Linear transition moment: <<A;EXCI>> = <0|A|n>
  !             if(MyRspFunc%order == 2) then
  !                MyRspFunc%trans_moment_order = 1 ! Single transition moment
  !                ! The last frequency (for EXCI) equals the excitation energy stored in rsp_eq_sol
  !                ! The first frequency (for A) then equals minus the excitation energy.
  !                MyRspFunc%freq(1) =-MyRspFunc%freq(2)
  !
  !                ! For each individual excited state the information needed to determine the
  !                ! linear transition moment for excitation number (j) is now stored in
  !                ! MyRspFunc and can be copied to the stack.
  !                rsp_counter = rsp_counter+1
  !                RspFunc_stack(rsp_counter) = MyRspFunc
  !             endif
  !
  !
  !             ! Quadratic transition moment (single or double)
  !             if(MyRspFunc%order == 3) then
  !
  !
  !
  !                ! QUADRATIC SINGLE TRANSITION MOMENT (e.g. <<A; B, EXCI>>)
  !                if(MyRspFunc%code(2) /= 'EXCI' .and. MyRspFunc%code(3) == 'EXCI') then
  !                   MyRspFunc%trans_moment_order = 1 ! Single transition moment
  !
  !                   k=0
  !                   Rfreq(1)=0E0_realk
  !                   Ifreq(1)=0E0_realk
  !
  !                   ! Read in frequency for operator B.
  !                   do k=1,AOR_stacksize-i
  !                      word=AOR_input(i+k)
  !
  !                      ! No frequencies defined for operator B, 
  !                      ! it is then zero by default initialization.
  !                      if(word(1:18) == 'TRANSITION-MOMENT' .or. &
  !                           & word(1:18) == 'RESPONSE-FUNCTION') exit
  !
  !                      ! Read in real and imaginary part of frequency
  !                      if(word(1:5) == '.FREQ') READ (AOR_input(i+k+1), *) order_temp, Rfreq(1)
  !                      if(word(1:7) == '.IMFREQ')  READ (AOR_input(i+k+1), *) order_temp, Ifreq(1)
  !                   enddo ! k-loop over AOR_input
  !
  !                   ! Freq3 = excitation energy (defined above)
  !                   ! Freq2 = input frequency (0 if nothing is specified)
  !                   ! Freq1 = -Freq2 - Freq3
  !                   MyRspFunc%freq(2) = (1E0_realk,0E0_realk)*Rfreq(1) + (0E0_realk,1E0_realk)*Ifreq(1)
  !                   MyRspFunc%freq(1) = - sum(MyRspFunc%freq(2:MyRspFunc%order))
  !
  !                   ! Finally, all information about the quadratic single 
  !                   ! transition moment for excitation number (j) is now stored in
  !                   ! MyRspFunc and can be copied to the stack.
  !                   rsp_counter = rsp_counter+1
  !                   RspFunc_stack(rsp_counter) = MyRspFunc
  !
  !                endif ! Single quadratic residue
  !
  !
  !
  !                ! QUADRATIC DOUBLE TRANSITION MOMENT (e.g. <<A; EXCI, EXCI>>) 
  !                ! Double residue defined as Eq. (3.21) of jcp 82, 3235 (1985):
  !                ! freq3 = (excitation energy m)
  !                ! freq2 = -(excitation energy q) 
  !                ! freq1 = -freq2-freq3
  !                ! We are now inside the j-loop which is over 
  !                ! the excited states for the second EXCI (operator C) in <<A; EXCI, EXCI>>. 
  !                ! We must also loop over the first exci (operator B), i.e. from 1 to nexci2.
  !
  !                if(MyRspFunc%code(2) == 'EXCI' .and. MyRspFunc%code(3) == 'EXCI') then
  !
  !                   MyRspFunc%trans_moment_order = 2 ! Double transition moment
  !
  !                   ! Determine nexci2 (number of excited states for operator B)
  !                   nexci2=0
  !                   do k=1,AOR_stacksize-i
  !                      word=AOR_input(i+k)
  !
  !                      if(word(1:18) == 'TRANSITION-MOMENT' .or. &
  !                           & word(1:18) == 'RESPONSE-FUNCTION') exit
  !
  !                      ! Read in real and imaginary part of frequency
  !                      if(word(1:8) == '.NEXCIT2') READ (AOR_input(i+k+1), *) nexci2
  !                   enddo ! loop: k=1,AOR_stacksize-i
  !
  !                   if(nexci2 ==0)  CALL lsQUIT('response_input: &                    
  !                        &The number of excitations for both poles must be specified &      
  !                        &for double residues using .NEXCIT and .NEXCIT2, respectively.',lupri)
  !
  !                   ! Loop over excited states for the second EXCI
  !                   do k=1,nexci2
  !
  !                      ! The excitation energies are stored in rsp_eq_sol(j)%fld(1)%freq,
  !                      ! see subroutine transition_moment_density_matrix below.
  !                      MyRspFunc%freq(3) = rsp_eq_sol(j)%fld(1)%freq
  !                      MyRspFunc%freq(2) =-rsp_eq_sol(k)%fld(1)%freq
  !                      MyRspFunc%freq(1) = - sum(MyRspFunc%freq(2:MyRspFunc%order))
  !                      MyRspFunc%ExNumber2 = k
  !
  !                      ! Finally, all information about the quadratic double
  !                      ! transition moment is now stored in
  !                      ! MyRspFunc and can be copied to the stack.
  !                      rsp_counter = rsp_counter+1
  !                      RspFunc_stack(rsp_counter) = MyRspFunc
  !
  !                   enddo ! loop: k=1,nexci2
  !
  !                endif ! quadratic double transition moment
  !
  !             endif ! quadratic transition moment (general)
  !
  !
  !             ! Cubic transition moment
  !             CubicTransMoment: if(MyRspFunc%order == 4) then
  !
  !                ! CUBIC SINGLE TRANSITION MOMENT (<<A; B, C, EXCI>>)
  !                if(MyRspFunc%code(3) /= 'EXCI' .and. MyRspFunc%code(4) == 'EXCI') then
  !                   MyRspFunc%trans_moment_order = 1 ! Single transition moment
  !
  !                   ! CUBIC DOUBLE TRANSITION MOMENT (<<A; B, EXCI, EXCI>>)
  !                elseif(MyRspFunc%code(3) == 'EXCI' .and. MyRspFunc%code(4) == 'EXCI') then
  !                   MyRspFunc%trans_moment_order = 2 ! Double transition moment
  !                else 
  !                   CALL lsQUIT('response_input: &
  !                        &For cubic single transition moment the fourth operator must be "EXCI". &
  !                        &E.g. 4  EL EL EL EXCI.',lupri)
  !                end if
  !
  !
  !                SingleCubicMoment: if(MyRspFunc%trans_moment_order == 1) then
  !
  !                   MyRspFunc%trans_moment_order = 1 ! Single transition moment
  !
  !                   k=0
  !                   ! Freq B
  !                   Rfreq(2)=0E0_realk
  !                   Ifreq(2)=0E0_realk
  !
  !                   ! Freq C
  !                   Rfreq(3)=0E0_realk
  !                   Ifreq(3)=0E0_realk
  !
  !                   ! Read in frequency for operator B and C.
  !                   do k=1,AOR_stacksize-i
  !                      word=AOR_input(i+k)
  !
  !                      ! No frequencies defined for operator B, 
  !                      ! it is then zero by default initialization.
  !                      if(word(1:18) == '.TRANSITION-MOMENT' .or. &
  !                           & word(1:18) == '.RESPONSE-FUNCTION') exit
  !
  !                      ! Read in real and imaginary parts of frequency
  !                      if(word(1:5) == '.FREQ') READ (AOR_input(i+k+1), *) order_temp, Rfreq(2), Rfreq(3)
  !                      if(word(1:7) == '.IMFREQ')  READ (AOR_input(i+k+1), *) order_temp, Ifreq(2), Ifreq(3)
  !                   enddo ! k-loop over AOR_input
  !
  !                   ! Freq4 = excitation energy (defined above)
  !                   ! Freq3 = input frequency(3) (0 if nothing is specified)
  !                   ! Freq2 = input frequency(2) (0 if nothing is specified)
  !                   ! Freq1 = -Freq2 - Freq3 - Freq4
  !                   MyRspFunc%freq(3) = (1E0_realk,0E0_realk)*Rfreq(3) + (0E0_realk,1E0_realk)*Ifreq(3)
  !                   MyRspFunc%freq(2) = (1E0_realk,0E0_realk)*Rfreq(2) + (0E0_realk,1E0_realk)*Ifreq(2)
  !                   MyRspFunc%freq(1) = - sum(MyRspFunc%freq(2:MyRspFunc%order))
  !
  !                   ! Finally, all information about the quadratic single 
  !                   ! transition moment for excitation number (j) is now stored in
  !                   ! MyRspFunc and can be copied to the stack.
  !                   rsp_counter = rsp_counter+1
  !                   RspFunc_stack(rsp_counter) = MyRspFunc
  !
  !                end if SingleCubicMoment
  !
  !
  !                DoubleCubicMoment: if(MyRspFunc%trans_moment_order == 2) then
  !
  !
  !                   ! QUADRATIC DOUBLE TRANSITION MOMENT (e.g. <<A; B, EXCI, EXCI>>) 
  !                   ! Double residue defined as Eq. (3.31) of jcp 82, 3235 (1985):
  !                   ! freq4 = (excitation energy m)
  !                   ! freq3 = (excitation energy n)
  !                   ! freq2 = input frequnecy
  !                   ! freq1 = -freq2-freq3-freq4
  !                   ! We are now inside the j-loop which is over 
  !                   ! the excited states for the second EXCI (operator D) in <<A; B, EXCI, EXCI>>. 
  !                   ! We must also loop over the first exci (operator C), i.e. from 1 to nexci2.
  !
  !
  !                   ! Determine nexci2 (number of excited states for operator C)
  !                   ! and frequency for operator B.
  !                   nexci2=0
  !                   Rfreq(2) = 0E0_realk
  !                   Ifreq(2) = 0E0_realk
  !
  !                   do k=1,AOR_stacksize-i
  !                      word=AOR_input(i+k)
  !
  !                      if(word(1:18) == 'TRANSITION-MOMENT' .or. &
  !                           & word(1:18) == 'RESPONSE-FUNCTION') exit
  !
  !                      if(word(1:8) == '.NEXCIT2') READ (AOR_input(i+k+1), *) nexci2
  !
  !                      ! Read in real and imaginary part of frequency for operator B
  !                      if(word(1:5) == '.FREQ') READ (AOR_input(i+k+1), *) order_temp, Rfreq(2)
  !                      if(word(1:7) == '.IMFREQ')  READ (AOR_input(i+k+1), *) order_temp, Ifreq(2)
  !                   enddo ! loop: k=1,AOR_stacksize-i
  !
  !                   if(nexci2 ==0)  CALL lsQUIT('response_input: &                    
  !                        &The number of excitations for both poles must be specified &      
  !                        &for double residues using .NEXCIT and .NEXCIT2, respectively.',lupri)
  !
  !                   ! Set freq2 equal to input frequency.
  !                   MyRspFunc%freq(2) = (1E0_realk,0E0_realk)*Rfreq(2) + (0E0_realk,1E0_realk)*Ifreq(2)
  !
  !                   ! Loop over excited states for the second EXCI
  !                   do k=1,nexci2
  !
  !                      ! The excitation energies are stored in rsp_eq_sol(j)%fld(1)%freq,
  !                      ! see subroutine transition_moment_density_matrix below.
  !                      MyRspFunc%freq(4) = rsp_eq_sol(j)%fld(1)%freq
  !                      MyRspFunc%freq(3) =-rsp_eq_sol(k)%fld(1)%freq
  !                      ! MyRspFunc%freq(2) is set above
  !                      MyRspFunc%freq(1) = - sum(MyRspFunc%freq(2:MyRspFunc%order))
  !                      MyRspFunc%ExNumber2 = k
  !
  !                      ! Finally, all information about the quadratic double
  !                      ! transition moment is now stored in
  !                      ! MyRspFunc and can be copied to the stack.
  !                      rsp_counter = rsp_counter+1
  !                      RspFunc_stack(rsp_counter) = MyRspFunc
  !
  !                   enddo ! loop: k=1,nexci2
  !
  !                end if DoubleCubicMoment
  !
  !
  !
  !
  !             endif CubicTransMoment
  !
  !
  !          enddo ! Loop over excitation energies
  !
  !
  !       end select
  !
  !
  !    enddo ! Loop over lines in AOR_input
  !
  !
  !  End subroutine response_input


  !> \brief Calculates response functions and transition moments according to the input
  !> information in RspFunc_stack and associate the pointer in RspFunc_stack%result
  !> to the response results.
  !> <br> <br>
  !> Details:
  !> <br>
  !> At the end of the subroutine,
  !> all response results have been saved in rsp_results. The response results for the
  !> first property specified in the input is saved in column one of rsp_results,
  !> the results for the second property in column 2 etc. 
  !> The arrangement in the individual columns is best explained by an example.
  !> Consider <<EL;MAG>> where the x,y, and z (1,2,3) components are specified
  !> for both EL and MAG.
  !> This yields 9 entries in the column in the order:
  !> <br> <br>
  !> (1,1) (1,2) (1,3) (2,1) (2,2) (2,3) (3,1) (3,2) (3,3)
  !> <br> <br>
  !> where the first index refers to EL and the second index refers to MAG.
  !> <br>
  !> Finally, RspFunc_stack(i)%result is associated with column (i) of rsp_results.
  !> \author Kasper Kristensen
  !> \date 2010-01
  subroutine rspfunc_or_transmoment(molcfg,F,D,S,RspFunc_stack,n_rspfunc)
    !> Info on molecule and other info needed by solver and integral routines
    type(rsp_molcfg), intent(inout)   :: molcfg
    !> Unperturbed Fock matrix
    type(matrix), intent(in)          :: F
    !> Unperturbed density matrix
    type(matrix), intent(in)          :: D
    !> Unperturbed overlap matrix
    type(matrix), intent(in)          :: S
    !> Stack of all response function and transition moments to be determined
    type(RspFunc), intent(inout)  :: RspFunc_stack(n_rspfunc)
    !> Total number of response functions and transition moments requested
    integer, intent(in)               :: n_rspfunc
    !> Local parameters
    integer                             :: maxdim, i, rsp_dim(n_rspfunc)
    complex(realk), allocatable, target :: rsp_results(:,:)

    maxdim=0

    ! Determine dimensions for the response tensors requested
    ! (e.g. for the response function <<EL;MAG>> with x,y, and z components
    ! for both EL and MAG, the dimension is 3x3=9).
    ! Maxdim is the maximum dimension requested
    do i=1,n_rspfunc
       rsp_dim(i) = product(RspFunc_stack(i)%dims(1:RspFunc_stack(i)%order))
       if(rsp_dim(i) > maxdim) maxdim = rsp_dim(i)
    enddo

    allocate(rsp_results(maxdim,n_rspfunc))

    do i=1,n_rspfunc

       rsp_results(:,:) = (0.0E0_realk, 0.0E0_realk)
       select case(RspFunc_stack(i)%order)

       case(2)
          ! Determine linear response function or transition moment
          call linear_response(molcfg,F,D,S,RspFunc_stack(i), &
               &rsp_results(1:rsp_dim(i),i), rsp_dim(i))

       case(3) 
          ! Determine quadratic response function or transition moment

          if(RspFunc_stack(i)%code(1) == 'GEO ') then !n+1 rule
             call quadratic_response_nplus1(molcfg,F,D,S,RspFunc_stack(i), &
                  &rsp_results(1:rsp_dim(i),i), rsp_dim(i))
          else ! 2n+1 rule
             call quadratic_response_2nplus1(molcfg,F,D,S,RspFunc_stack(i), &
                  &rsp_results(1:rsp_dim(i),i), rsp_dim(i))
          endif

       case(4) 
          ! Determine cubic response function or transition moment
          ! Always use 2n+1 rule for cubic response
          call cubic_response(molcfg,F,D,S,RspFunc_stack(i), &
               &rsp_results(1:rsp_dim(i),i), rsp_dim(i))

       case default
          CALL lsQUIT('rspfunc_or_transmoment: Only linear, quadratic, and cubic response &
               &are currently implemented',molcfg%lupri)
       end select

       ! Associate result-pointer with the response results.
       nullify(RspFunc_stack(i)%result)
       allocate(RspFunc_stack(i)%result(rsp_dim(i)))
       RspFunc_stack(i)%result(1:rsp_dim(i)) = rsp_results(1:rsp_dim(i),i)

    enddo

    deallocate(rsp_results)
  end subroutine rspfunc_or_transmoment


  !> \brief Determines the number of response functions and excitation energies
  !> requested in input.
  !> \author Kasper Kristensen
  !> \date 2010-01
  !  subroutine determine_number_of_rspfuncs_and_exenergies(n_rspfunc, AOR_input, &
  !       & nexci_max, ExEnergies_requested)
  !    !> Total number of response functions and transition moments requested
  !    integer, intent(inout)          :: n_rspfunc
  !    !> Input from *AORESPONSE to *END AORESPONSE in the dalton input file - but
  !    !> stripped from (potential) comment lines
  !    character(len=40), intent(in)   :: AOR_input(AOR_stacksize)
  !    !> The position of the highest lying excited state requested
  !    integer, intent(inout)          :: nexci_max
  !    !> Do any excited state energies need to be determined?
  !    logical, intent(inout)          :: ExEnergies_requested 
  !    !> Local parameters
  !    character(len=40)               :: word1, word2
  !    integer                         :: nexci,nexci2,i,j
  !
  !    n_rspfunc=0
  !    nexci_max=0
  !    nexci=0
  !
  !    ! Determine maximum number of excitation energies nexci_max
  !    ! requested for any transition moment property.
  !    do i=1,AOR_stacksize
  !       word1 = AOR_input(i)
  !       if(word1(1:7) == '.NEXCIT') then
  !          READ(AOR_input(i+1), '(I10)') nexci
  !          if(nexci_max < nexci) nexci_max =nexci
  !       endif
  !    enddo
  !
  !    ! If .NEXCIT has not been declared, default: determine one excitation energy
  !    if(nexci_max == 0) nexci_max =1
  !
  !    ! Count total number of response functions and transition moments
  !    do i=1,AOR_stacksize
  !       word1 = AOR_input(i)
  !
  !       Select case(word1(1:18))
  !
  !       case('.RESPONSE-FUNCTION')
  !          n_rspfunc= n_rspfunc+1
  !       case('.TRANSITION-MOMENT')
  !          ExEnergies_requested = .true.
  !
  !          nexci=nexci_max
  !          nexci2=0
  !          ! Determine the number of excitation energies. 
  !          do j=1,AOR_stacksize-i
  !             word2 = AOR_input(i+j)
  !
  !             if(word2(1:18) == '.TRANSITION-MOMENT' &
  !                  & .or. word2(1:18) == '.RESPONSE-FUNCTION') exit
  !
  !
  !             ! nexci specified in the next line
  !             if(word2(1:8) == '.NEXCIT ') READ (AOR_input(i+j+1:i+j+1), '(I10)') nexci
  !
  !             ! Only for double residues, the number of excitation energies to be determined
  !             ! for the second set of excitation energies
  !             ! [the number of w_q's in Eq. 3.21 in J. Chem. Phys. 82, 3235 (1985)] 
  !             if(word2(1:8) == '.NEXCIT2') READ (AOR_input(i+j+1:i+j+1), '(I10)') nexci2
  !
  !          enddo ! j-loop
  !
  !          if(nexci2 == 0) then ! single transition moment
  !             n_rspfunc=n_rspfunc+nexci
  !          else ! double transition moment
  !             n_rspfunc=n_rspfunc+nexci*nexci2
  !          endif
  !
  !       end select
  !
  !    enddo ! i-loop
  !
  !  end subroutine determine_number_of_rspfuncs_and_exenergies


  !> \brief Initializes all parameters in MyRspFunc.
  !> \author Kasper Kristensen
  !> \date 2010-01
  subroutine initialize_RspFunc(MyRspFunc)
    !> MyRspFunc is set to the default values:
    !> <br> <br>
    !> order=0
    !> <br>
    !> Code for all operators='NOOP'  (No operator)
    !> <br>
    !> first component = 1
    !> <br>
    !> dimension = 1  (total number of components)
    !> <br>
    !> All (complex) frequencies = (0,0)
    !> <br>
    !> trans_moment=.false. (response function, not transition moment)
    !> <br>
    !> trans_moment_order=0 (response function, not transition moment)
    !> <br>
    !> ExNumber1=0 (response function, not transition moment)
    !> <br>
    !> ExNumber2=0 (response function, not transition moment)

    type(RspFunc)    :: MyRspFunc
    Integer           :: i

    MyRspFunc%order=0

    do i=1,max_order
       MyRspFunc%code(i) = 'NOOP'
       MyRspFunc%first_comp(i) = 1
       MyRspFunc%dims(i) = 1
       MyRspFunc%freq(i)=(0E0_realk,0E0_realk)
    enddo

    MyRspFunc%trans_moment = .false.
    MyRspFunc%trans_moment_order = 0
    MyRspFunc%ExNumber1 = 0
    MyRspFunc%ExNumber2 = 0
  end subroutine initialize_RspFunc


  !> \brief The information in the string about the order of the response function
  !> and the codes for the operators are read into 'order' and the 'code'-vector, respectively.
  !> The formal size of the code-vector is max_order. 
  !> The indices not used are labeled 'NOOP'. E.g. if max_order=4 we write the quadr. RF:
  !> <<EL;EL,MAG>> = (/'EL  ', 'EL  ', 'MAG ', 'NOOP'/)
  !> \author Kasper Kristensen
  !> \date 2010-01
  subroutine read_words_from_string(string,order,code,num_aux)
    !> Code for the operators in the response function
    Character(len=4), intent(inout) :: code(max_order)
    !> String specified in dalton input containing order and codes.
    character(len=39), intent(in)   :: string
    !> Order of response funciton.
    integer, intent(in)             :: order
    !> Counter used when perturbations are specified using the HERMIT integral label
    !> (see pert_table in rsp_contribs.F90).
    integer, intent(inout)          :: num_aux
    !> Local variables
    character(len=8)                :: code_temp
    integer                         :: i,j,k


    k=0

    ! Loop over letters in string
    do i=1,38

       ! code_temp contains 8 letters to be able to hold input like XDIPLEN.
       code_temp(1:8) = '        '

       ! If the position in string is NOT an empty space we read in the code for the operator
       ! ('EL' for electric dipole etc., see pert_table in rsp_contribs.F90).
       ! Two string letters of the form ' X' where X is the beginning of the operator code
       ! indicates that a new operator is written in the string.
       if(string(i:i) == ' ' .and. string(i+1:i+1) /= ' ') then
          j=1

          ! Read in operator in code_temp
          do
             if(string(i+j:i+j) == ' ') exit
             code_temp(j:j) = string(i+j:i+j)
             j=j+1
          enddo

          k=k+1

          ! If code_temp contains 4 letters or less it can be directly copied to code(k)
          ! using the nomenclature in pert_table in rsp_contribs.F90.
          if(code_temp(5:8) == '    ') then
             code(k) = code_temp(1:4)


             ! Else, e.g. for input XDIPLEN, we must assign a 4-letter code AUXI (I=0,9)
             ! to prop_auxlab in rsp_contribs.F90.
             ! For example we may have prop_auxlab(2)='XDIPLEN ', then code(k)='AUX2'.
             ! In this case code(k) implicitly refers to XDIPLEN via entry 2
             ! in prop_auxlab.
          else 

             call lsquit('error in read_words_from_string: '// &
                  &'long operator lables are no longer permitted in LSresponse',-1)

          endif ! endif: 4 letter code or "old" code with more letters

       endif ! endif: operator specified in string

    Enddo ! i-loop

    if(k /= order) then
       CALL lsQUIT('read_words_from_string: k /= order of response function. &
            &Input error in specifying the number of operators: &
            &The order of the response function must equal the number of operators - e.g. &
            & 3 EL MAG MAG',-1)
    endif

  end subroutine read_words_from_string


  !> \brief Simply replaces the letter oldletter by newletter in string.
  !> (E.g. replace ':' by ' ').
  !> \author Kasper Kristensen
  !> \date 2010-01
  subroutine clean_up_string(string, oldletter, newletter)
    !> Input string
    character(len=40), intent(inout)   :: string
    !> Old letter to be replaced
    character(len=1),intent(in)        :: oldletter
    !> New letter to replace old letter.
    character(len=1),intent(in)        :: newletter
    integer                            :: i

    do i=1,40
       if(string(i:i) == oldletter) string(i:i) = newletter
    enddo

  end subroutine clean_up_string


  !> \brief This routine calculates
  !> transition density matrices
  !> <br> <br>
  !> Dx(i) = D*S*X(i) - X(i)*S*D
  !> <br> <br>
  !> for (i=1,nexci_max), where X(i) is a solution to eigenvalue problem
  !> <br> <br>
  !> E2*X(i) = ExEnergy(i)*S2*X(i).
  !> <br> <br>
  !> Dx(i) is saved in rsp_eq_sol(i)\%D in rsp_equations.f90.
  !> <br>
  !> ExEnergies(:) is saved and rsp_eq_sol(i)\%fld(1)\%freq in rsp_equations is assigned to ExEnergies(i).
  !> \author Kasper Kristensen
  !> \date 2010-01 
  subroutine transition_moment_density_matrix(molcfg,F,D,S,nexci_max)
    !> Info on molecule needed by solver and integral programs
    type(rsp_molcfg),target,intent(inout) :: molcfg
    !> Unperturbed Fock matrix
    type(matrix), intent(in)       :: F
    !> Unperturbed density matrix
    type(matrix), intent(in)       :: D
    !> Unperturbed overlap matrix
    type(matrix), intent(in)       :: S
    !> Position of the highest lying excited state requested.
    integer, intent(in)            :: nexci_max
    !> Dx will contain the transition densities for nexci_max excited states.
    type(matrix)                   :: Dx(nexci_max)
    !> ExEnergies will contain the excitation energies for nexci_max excited states.
    real(realk)                    :: ExEnergies(nexci_max)
    !> Local parameters
    type(matrix)                   :: B(1)  ! dummy matrix for calling solver 
    integer                        :: nstart, i, j
    type(decompItem), pointer      :: decomp
    integer, pointer               :: lupri

    ! Rsp_solver has decomp as inout, but molcfg is in, so re-point
    lupri => molcfg%lupri

    ! Initializing B and Dx, necessary for calling response solver
    B(1) = 0*D
    call mat_ensure_alloc(B(1))

    do I = 1,nexci_max
       Dx(I) = 0*D
       call mat_ensure_alloc(Dx(I))
    enddo

    ! Set number of start vectors for solving eigenvalue equation.
    if (molcfg%decomp%cfg_startvectors) then
       nstart=molcfg%decomp%cfg_no_of_startvectors
    else ! Default: Number of startvectors is the number of excitation energies.
       nstart = nexci_max
    endif

    if (molcfg%solver%rsp_stdnew) then
       ! Initialize solver parameters.
       call rsp_sym_init(nexci_max,1,nexci_max,nexci_max,nstart)
       ! Calling solver. 
       call rsp_sym_solver(molcfg,D,S,F,.false.,nexci_max,B,ExEnergies,Dx)
    else
       ! Initialize solver parameters.
       call rsp_init(nexci_max,1,nexci_max,nexci_max,nstart)
       ! Calling solver. 
       call rsp_solver(molcfg,D,S,F,.false.,nexci_max,B,ExEnergies,Dx)
    endif
    ! At this point Dx contain the excitation vectors, not the transition densities.

    write(lupri,*)
    write(lupri,*) 'The following excitation energies are used for the response calculation:'
    write(lupri,*)

    do i=1,nexci_max
       write(lupri,*) i, ExEnergies(i)
    enddo

    ! Ensuring rsp_eq_sol in module rsp_equations is large enough
    ! ajt FIXME use subroutine for inserting into rsp_eq_sol
    call rsp_eq_sol_grow(nexci_max)

    j = 1
    Do i=1,nexci_max
       ! Finding vacant slot in rsp_eq_sol
       do while (associated(rsp_eq_sol(j)%mol))
          if (j==size(rsp_eq_sol)) &
               call lsquit('Unexpected error in transition_moment_density_matrix' &
               & // ': rsp_eq_sol is full',lupri)
          j = j+1
       end do
       ! Setting reference to molecule + rsp config
       rsp_eq_sol(j)%mol => molcfg
       ! Setting up field
       allocate(rsp_eq_sol(j)%fld(1))
       rsp_eq_sol(j)%fld(1)%label='EXCI'
       rsp_eq_sol(j)%fld(1)%comp=1
       rsp_eq_sol(j)%fld(1)%ncomp=1
       rsp_eq_sol(j)%fld(1)%freq=ExEnergies(i) * (1E0_realk,0E0_realk)
       ! Transforming Dx, such that it becomes the transition DENSITY.
       rsp_eq_sol(j)%D = D*S*Dx(i) - Dx(i)*S*D
       Dx(i) = 0  ! free
    enddo

    ! Free B(1) matrix
    B(1)=0
  end subroutine transition_moment_density_matrix

  !> \brief Subroutine for calculating first order properties such as dipole and gradient.
  !> \author Kasper Kristensen
  !> \date 2010-08 
  subroutine Get_first_order_property(molcfg,F,D,S,property,dim_property,property_type)
    !ajt PS This line is redundant if there is "implicit none"
    implicit none  !in the encompassing module
    !> Info on molecule and other info needed by solver and integral routines
    type(rsp_molcfg), intent(inout)  :: molcfg
    !> Unperturbed Fock matrix
    type(matrix), intent(in)         :: F
    !> Unperturbed density matrix
    type(matrix), intent(in)         :: D
    !> Unperturbed overlap matrix
    type(matrix), intent(in)         :: S
    !> Length of property vector
    integer, intent(in) :: dim_property
    !> Vector containing calculated property
    complex(realk), intent(inout)    :: property(dim_property)
    !> 4-index code for which property to be determined (e.g. 'EL  ' for dipole)
    character(len=4), intent(in)     :: property_type
    type(RspFunc) :: MyRspFunc
    type(matrix) :: W
    logical                        :: pdbs(1)


    call initialize_RspFunc(MyRspFunc)

    ! Order of "response function" is just one - simple expectation value
    MyRspFunc%order = 1

    ! Which property operator
    MyRspFunc%code(1) = property_type

    ! Calculate components specified by dim_property
    MyRspFunc%first_comp(1) = 1
    MyRspFunc%dims(1) = dim_property

    ! Zero property before calling
    property  = (0E0_realk,0E0_realk)

    ! Determine if the operator use perturbation dependent basis sets
    pdbs=pert_basdep(MyRspFunc%code(1:1))


    ! One-electron contributions
    if(pdbs(1)) then ! Sx*DFD contribution
       W = D*F*D
       call rsp_oneave(molcfg, S, MyRspFunc%code(1:1), (/D/), shape(property), &
            & property,comp= (/MyRspFunc%first_comp(1)/), DFD=(/W/) )
    else
       call rsp_oneave(molcfg, S, MyRspFunc%code(1:1), (/D/), shape(property), &
            & property,comp= (/MyRspFunc%first_comp(1)/) )
    end if


    ! Two-electron contributions
    call rsp_twoave(molcfg, MyRspFunc%code(1:1), (/D/), &
         & shape(property), property, &
         & comp=(/MyRspFunc%first_comp(1)/) )



    ! Free W if used
    if(pdbs(1)) W=0

  end subroutine Get_first_order_property





  !> \brief Calculates linear response function or linear transition moment
  !> according to the response information in MyRspFunc and saves
  !> the results in the rsp_results_vec vector.
  !> <br> <br>
  !> Details: 
  !> <br>
  !> We use the AO density matrix-based response formulation in J. Chem. Phys. 129, 214108,
  !> see in particular Eq. 207.
  !> <br>
  !> rsp_results_vec will contain all response results, possibly for several
  !> field components, but always for only one frequency (or excitation energy).
  !> The arrangement of the results is best explained by an example.
  !> Consider <<EL;MAG>> where the x,y, and z (1,2,3) components are specified
  !> for both EL and MAG.
  !> In this case the 3x3=9 components 
  !> in rsp_results_vec are arranged as follows:
  !> <br> <br>
  !> (1,1) (1,2) (1,3) (2,1) (2,2) (2,3) (3,1) (3,2) (3,3)
  !> <br> <br>
  !> where the first index refers to EL and the second index refers to MAG.
  !> <br> <br>
  !> NOTE! 
  !> <br>
  !> In case you want to calculate a linear transition moment by 
  !> setting up MyRspFunc in a separate subroutine (not in the input),
  !> and calling linear_response.
  !> You must then first call transition_moment_density_matrix in response_driver_module
  !> to determine the transition density matrix, and then call linear_response.
  !> <br>
  !> On the other hand, in a linear response function calculation,
  !> the perturbed density matrix is calculated on the fly by linear_response itself,
  !> and you may obtain a linear response function by a single call to linear_response.
  !> \author Kasper Kristensen
  !> \date 2010-01 
  subroutine linear_response(molcfg,F,D,S,MyRspFunc, rsp_results_vec,response_size)
    implicit none
    !> Info on molecule and other info needed by solver and integral routines
    type(rsp_molcfg), intent(inout)   :: molcfg
    !> Unperturbed Fock matrix
    type(matrix), intent(in)          :: F
    !> Unperturbed density matrix
    type(matrix), intent(in)          :: D
    !> Unperturbed overlap matrix
    type(matrix), intent(in)          :: S
    !> Contains all information about the linear 
    !> response function or transition moment to be determined.
    type(RspFunc), intent(inout)      :: MyRspFunc
    !> Dimension of response vector
    integer, intent(in)               :: response_size
    !> Vector for collecting response results
    complex(realk), intent(inout)     :: rsp_results_vec(response_size)
    !> Local parameters
    type(matrix)                      :: Db(MyRspFunc%dims(2)), &
         & Fb(MyRspFunc%dims(2)), Wb(MyRspFunc%dims(2)), W(1)
    complex(realk)                 :: lin_RspFunc(MyRspFunc%dims(1),MyRspFunc%dims(2))
    logical                        :: pdbs(2)
    integer                        :: i,j, k
    character*4                    :: NOOP(1)
    integer, pointer               :: lupri
    NOOP(1) = 'NOOP'
    lupri => molcfg%lupri

    write(lupri, *) 'Starting linear response function calculation for: ', &
         & '<< ',MyRspFunc%code(1), '; ', MyRspFunc%code(2),'>>'

    write(lupri,*)

    ! Determine if the operators in fields contain perturbation dependent basis sets
    ! E.g. if MyRspFunc%code(1:2)= ('EL', 'GEO') then pdbs=(.false. , .true.).
    pdbs=pert_basdep(MyRspFunc%code(1:2))

    ! Calculate perturbed density matrix for operator 2 (Db), 
    ! only required if a linear response function
    ! is requested. If a linear transition moment is requested the transition density matrix
    ! is already stored in rsp_eq_sol%D in rsp_equations.f90 in position "i",
    ! where "i" is the number of the excitation in question (MyRspFunc%ExNumber1).
    if(MyRspFunc%trans_moment) then
       Db(1) = rsp_eq_sol(MyRspFunc%ExNumber1)%D
    else
       call pert_dens(molcfg, S, (/ MyRspFunc%code(2) /), (/MyRspFunc%dims(2)/), (/D/), (/F/), &
            & Db, Fb, comp=(/MyRspFunc%first_comp(2)/), freq=(/MyRspFunc%freq(2)/) )
    endif

    ! **********************************************************************
    ! *                START LINEAR RESPONSE CALCULATION                   *
    ! **********************************************************************

    ! As described in [jcp, 129, 214108 (2008)] there are four contributions
    ! to the linear response function, see in particular Eq. 207:
    !
    ! <<A;B>> = Tr { E^{0,ab} - Sab W } + 
    !         + Tr { E^{1,a} Db - Sa Wb }
    !
    ! Where all quantities needed are defined after Eq. 207 of that paper.

    lin_RspFunc(:,:) = (0E0_realk,0E0_realk)

    !************************************************************************
    ! E^{0,ab} - S^{ab}*W contributions                                     !
    ! (only contributes for PDBS and never for transition moment)           !
    !************************************************************************
    ! Determine W
    if(  (.not. MyRspFunc%trans_moment) .and. (any(pdbs)) ) then

       ! Calculate zeroth order W
       call get_W_matrix(F,D,S,MyRspFunc,0,pdbs,W)

       ! One-electron contributions to E^{0,ab} - S^{ab}*W
       call rsp_oneave(molcfg, S, MyRspFunc%code(1:2), (/D/), shape(lin_RspFunc), lin_RspFunc, &
            & comp= (/MyRspFunc%first_comp(1), MyRspFunc%first_comp(2)/), &
            & freq= (/MyRspFunc%freq(1), MyRspFunc%freq(2) /), DFD=(/W/) )

       ! Two-electron contributions to E^{0,ab} - S^{ab}*W
       call rsp_twoave(molcfg, MyRspFunc%code(1:2), (/D/), shape(lin_RspFunc), lin_RspFunc, &
            & comp= (/MyRspFunc%first_comp(1), MyRspFunc%first_comp(2)/) )

       ! Free W
       if(isdef(W(1))) W(1)=0

    endif


    !************************************************************************
    ! E^{1,a}*Db - S^{a}*Wb contributions                                   !
    !************************************************************************

    ! Only if perturbation 1 use PDBS (Sa /= 0) we need to determine Wb
    ! using Eq. 213 in jcp, 129, 214108 (2008)
    if(pdbs(1)) then 

       ! If Response function, the perturbed Fock matrix has already been determined above
       ! when calling pert_dens.
       ! However, if transition moment, we need to calculate the transition Fock matrix Fb
       ! manually now to determine Wb below.
       if(MyRspFunc%trans_moment) then
          Fb(1)=0E0_realk*F
          call rsp_twoint(molcfg, NOOP, (/D,Db(1)/), (/1/), Fb(1))
       endif

       ! Determine Wb [Eq. 213 in jcp, 129, 214108 (2008)]
       call get_W_matrix(F,D,S,MyRspFunc,1,pdbs,Wb,ops=(/2/),Dx=Db,Fx=Fb)

       ! One electron contributions to E^{1,a}*Db - S^{a}*Wb
       call rsp_oneave(molcfg, S, (/MyRspFunc%code(1)/), (/Db/), &
            & shape(lin_RspFunc), lin_RspFunc, &
            & comp=(/MyRspFunc%first_comp(1)/), &
            & freq=(/MyRspFunc%freq(1)/), DFD=(/Wb/) )

       ! Two electron contributions to E^{1,a}*Db - S^{a}*Wb 
       call rsp_twoave(molcfg, (/ MyRspFunc%code(1) /), (/D, Db/), &
            & shape(lin_RspFunc), lin_RspFunc, &
            & comp=(/MyRspFunc%first_comp(1)/) )

       if(isdef(Wb(1))) Wb(:)=0

    else ! If not pdbs(1) there is no (Sa Wb)-contribution (since Sa=0) and no two-electron terms
       ! [see Eq. 209 in jcp, 129, 214108 (2008)]

       ! One electron contributions to E^{1,a}*Db 
       call rsp_oneave(molcfg, S, (/MyRspFunc%code(1)/), (/Db/), &
            & shape(lin_RspFunc), lin_RspFunc, &
            & comp=(/MyRspFunc%first_comp(1)/), &
            & freq=(/MyRspFunc%freq(1)/) )

    endif


    ! Done calculating quadratic response function.
    ! Copying results in lin_RspFunc to rsp_results_vec according to the ordering described
    ! at the beginning of this subroutine.
    k=0
    do i=1,MyRspFunc%dims(1)
       do j=1,MyRspFunc%dims(2)
          k=k+1
          rsp_results_vec(k) = lin_RspFunc(i,j)
       enddo
    enddo

    ! Free matrices
    Db(:)=0
    if(isdef(Fb(1))) Fb(:)=0

  end subroutine linear_response


  !> \brief Using the 2n+1 rule (determine Da, Db, and Dc)
  !> this subroutine calculates quadratic response function or transition moments
  !> obtained from single or double residues of quadratic response functions
  !> according to the response information in MyRspFunc and saves
  !> the results in the rsp_results_vec vector.
  !> <br> <br>
  !> Details: 
  !> <br>
  !> We use the AO density matrix-based response formulation in J. Chem. Phys. 129, 214108,
  !> see in particular Eqs. 255-258 for the quadratic response function using the 2n+1 rule.
  !> Single and double residues are found from Eq. (258) using the strategy presented in Sec. IV.H.
  !> <br>
  !> rsp_results_vec will contain all response results, possibly for several
  !> field components, but always for only one set of frequencies (or excitation energies).
  !> The arrangement of the results is best explained by an example.
  !> Consider <<A;B,C>> where the x,y, and z (1,2,3) components are specified
  !> for all operators.
  !> In this case the 3x3x3=27 components 
  !> in rsp_results_vec are arranged as follows (read from left to right):
  !> <br> <br>
  !> (1,1,1) (1,1,2) (1,1,3) (1,2,1) (1,2,2) (1,2,3) (1,3,1) (1,3,2) (1,3,3)
  !> <br>
  !> (2,1,1) (2,1,2) (2,1,3) (2,2,1) (2,2,2) (2,2,3) (2,3,1) (2,3,2) (2,3,3)
  !> <br>
  !> (3,1,1) (3,1,2) (3,1,3) (3,2,1) (3,2,2) (3,2,3) (3,3,1) (3,3,2) (3,3,3)
  !> <br> <br>
  !> where the first, second, and third indices refer to A, B, and C, respectively.
  !> <br> <br>
  !> NOTE! 
  !> <br>
  !> In case you want to calculate a quadratic transition moment by 
  !> setting up MyRspFunc in a separate subroutine (not in the input),
  !> and calling quadratic_response_2nplus1:
  !> You must then first call transition_moment_density_matrix in response_driver_module
  !> to determine the transition density matrix, and then call linear_response.
  !> <br>
  !> On the other hand, in a quadratic response function calculation,
  !> the perturbed density matrix is calculated on the fly by quadratic_response_2nplus1 itself,
  !> and you may obtain a quadratic response function by a single call to this routine.
  !> \author Kasper Kristensen
  !> \date 2010-03
  subroutine quadratic_response_2nplus1(molcfg,F,D,S,MyRspFunc, rsp_results_vec,response_size)
    implicit none
    !> References to molecule input, and solver and integral settings
    type(rsp_molcfg), intent(inout)   :: molcfg
    !> Unperturbed Fock matrix
    type(matrix), intent(in)          :: F
    !> Unperturbed density matrix
    type(matrix), intent(in)          :: D
    !> Unperturbed overlap matrix
    type(matrix), intent(in)          :: S
    !> Contains all information about the quadratic
    !> response function or transition moment to be determined.
    type(RspFunc), intent(inout)      :: MyRspFunc
    !> Dimension of response vector
    integer,intent(in) :: response_size
    !> Vector for collecting response results
    complex(realk), intent(inout)     :: rsp_results_vec(response_size)
    !> Local parameters
    type(matrix)                      :: Da(MyRspFunc%dims(1)), &
         & Fa(MyRspFunc%dims(1)), Wa(MyRspFunc%dims(1)), &
         & Db(MyRspFunc%dims(2)), &
         & Fb(MyRspFunc%dims(2)), Wb(MyRspFunc%dims(2)), &
         & Dc(MyRspFunc%dims(3)), &
         & Fc(MyRspFunc%dims(3)), Wc(MyRspFunc%dims(3)), &
         & Sa(MyRspFunc%dims(1)), Sb(MyRspFunc%dims(2)), Sc(MyRspFunc%dims(3)), &
         & Wbc(1), Ybc(1), Zbc(1), lambda(MyRspFunc%dims(1)), &
         & zeta(MyRspFunc%dims(1)), W(1), zeroD, zeroF
    complex(realk)                 :: quad_RspFunc(MyRspFunc%dims(1),MyRspFunc%dims(2), MyRspFunc%dims(3))
    logical                        :: pdbs(3)
    integer                        :: i,j, k, l, b,c
    character*4                    :: NOOP(1)
    logical :: AeqB, AeqC, BeqC
    integer, pointer               :: lupri
    NOOP(1) = 'NOOP'

    lupri => molcfg%lupri

    write(lupri, *) 'Starting quadratic response function calculation for: ', &
         & '<< ',MyRspFunc%code(1), '; ', MyRspFunc%code(2), ', ', MyRspFunc%code(3), '>>'

    write(lupri,*)

    ! Determine if the operators in fields contain perturbation dependent basis sets
    ! E.g. if MyRspFunc%code(1:3)= ('GEO','EL', 'EL') then pdbs=(.true. , .false., .false.).
    pdbs=pert_basdep(MyRspFunc%code(1:3))

    ! Determine if any of the operators are identical
    AeqB = Identical_operators(MyRspFunc,1,2)
    AeqC = Identical_operators(MyRspFunc,1,3)
    BeqC = Identical_operators(MyRspFunc,2,3)


    ! *******************************************************************************************
    !                    Calculate perturbed density matrices Da, Db, and Dc.                   !
    ! *******************************************************************************************
    ! Da
    call pert_dens(molcfg, S, (/ MyRspFunc%code(1) /), (/MyRspFunc%dims(1)/), (/D/), (/F/), &
         & Da, Fa, comp=(/MyRspFunc%first_comp(1)/), freq=(/MyRspFunc%freq(1)/) )
    ! Db
    if(MyRspFunc%trans_moment .and. MyRspFunc%trans_moment_order ==2) then ! double residue
       ! Db = -Dx^T
       ! where Dx contain the transition moment density matrices.
       ! Due to the transposition we effectively set omega_b = - excitation energy
       ! (see the transition_moment_density_matrix subroutine and the conventions for 
       ! residues in the beginning of the module and Sec. IV.H. of jcp, 129, 214108).
       Db(1) = -trans(rsp_eq_sol(MyRspFunc%ExNumber2)%D)
       Fb(1)=0E0_realk*F
       call rsp_twoint(molcfg, NOOP, (/D,Db(1)/), (/1/), Fb(1))
    else ! if single residue or ordinary quadratic response function
       ! KK quick fix. If A and B refer to the same operators and frequencies,
       ! we simply copy the matrices already calculated for A.

       AeqBorNot: if(AeqB) then
          do i=1,MyRspFunc%dims(1) - MyRspFunc%first_comp(1) + 1
             Db(i)=Da(i)
             Fb(i)=Fa(i)
          end do
       else
          call pert_dens(molcfg, S, (/ MyRspFunc%code(2) /), (/MyRspFunc%dims(2)/), (/D/), (/F/), &
               & Db, Fb, comp=(/MyRspFunc%first_comp(2)/), freq=(/MyRspFunc%freq(2)/) )
       end if AeqBorNot

    endif


    ! Dc matrix
    if(MyRspFunc%trans_moment) then ! If perturbation c is EXCI (single or double residue)
       ! Dc = Dx
       ! because we set omega_c = omega_m
       ! (see the transition_moment_density_matrix subroutine and the conventions for 
       ! residues in the beginning of the module and Sec. IV.H. if jcp, 129, 214108).
       Dc(1) = rsp_eq_sol(MyRspFunc%ExNumber1)%D
       ! Determine Fc
       Fc(1)=0E0_realk*F
       call rsp_twoint(molcfg, NOOP, (/D,Dc(1)/), (/1/), Fc(1))
    else ! ordinary quadratic response function
       ! KK quick fix. If B and C refer to the same operators and frequencies,
       ! we simply copy the matrices already calculated for B.
       ! Should eventually be done creating work list of response functions.

       BeqCorNot: if(BeqC) then
          do i=1,MyRspFunc%dims(2) - MyRspFunc%first_comp(2) + 1
             Dc(i)=Db(i)
             Fc(i)=Fb(i)
          end do
       else
          call pert_dens(molcfg, S, (/ MyRspFunc%code(3) /), (/MyRspFunc%dims(3)/), (/D/), (/F/), &
               & Dc, Fc, comp=(/MyRspFunc%first_comp(3)/), freq=(/MyRspFunc%freq(3)/) )
       end if BeqCorNot
    end if

    ! *******************************************************************************************
    !                    Calculate matrices for perturbation-dependent basis sets               !
    ! *******************************************************************************************

    ! Determine Sa
    ! KK fix: Added "comp=..." in rsp_oneint call. Necessary if the first component is not 1... 
    if(pdbs(1)) then
       call rsp_oneint(molcfg, S, (/MyRspFunc%code(1) /), (/MyRspFunc%dims(1) /), &
            comp=(/MyRspFunc%first_comp(1)/), S=Sa)
    else
       do i=1,MyRspFunc%dims(1)
          Sa(i)=0E0_realk*S
       enddo
    endif

    ! Determine Sb
    if(pdbs(2)) then
       call rsp_oneint(molcfg, S, (/MyRspFunc%code(2) /), (/MyRspFunc%dims(2) /), &
            comp=(/MyRspFunc%first_comp(2)/), S=Sb)
    else
       do i=1,MyRspFunc%dims(2)
          Sb(i)=0E0_realk*S
       enddo
    endif


    ! Determine Sc
    ! Setting Sc=0 for trans_moment ensures that possible terms containing Sc
    ! do not contribute to the transition moment.
    if(pdbs(3) .and. .not. MyRspFunc%trans_moment) then
       call rsp_oneint(molcfg, S, (/MyRspFunc%code(3) /), (/MyRspFunc%dims(3) /), &
            comp=(/MyRspFunc%first_comp(3)/), S=Sc)
    else
       do i=1,MyRspFunc%dims(3)
          Sc(i)=0E0_realk*S
       enddo
    endif



    ! **********************************************************************
    ! **********************************************************************
    ! *                START QUADRATIC RESPONSE CALCULATION                *
    ! **********************************************************************
    ! **********************************************************************

    quad_RspFunc(:,:,:) = (0E0_realk,0E0_realk)

    !************************************************************************
    ! E^{0,abc} - S^{abc}*W contributions                                   !
    ! (only contributes for PDBS and never for transition moment)           !
    !************************************************************************

    ! One-electron contributions to E^{0,abc} - S^{abc}*W
    if(any(pdbs(2:3)) .and. .not. MyRspFunc%trans_moment) then

       ! Calculate zeroth order W
       call get_W_matrix(F,D,S,MyRspFunc,0,pdbs,W)

       call rsp_oneave(molcfg, S, MyRspFunc%code(1:3), (/D/), shape(quad_RspFunc), quad_RspFunc, &
            & comp= (/ MyRspFunc%first_comp(1:3) /), &
            & freq= (/ MyRspFunc%freq(1:3) /), DFD=(/W/) )

       ! Two-electron contributions to E^{0,abc}
       ! only if ALL perturbation use PDBS.
       if(all(pdbs)) then
          call rsp_twoave(molcfg, MyRspFunc%code(1:3), (/D/), shape(quad_RspFunc), quad_RspFunc, &
               & comp= (/ MyRspFunc%first_comp(1:3) /) )
       endif

       if(isdef(W(1))) W(1)=0

    endif
    !************************************************************************
    ! E^{1,bc}*Da - S^{bc}*Wa contributions                                 !
    !************************************************************************
    ! Never for transition moments.

    if( .not. MyRspFunc%trans_moment ) then

       ! Get Wa
       call get_W_matrix(F,D,S,MyRspFunc,1,pdbs,Wa,ops=(/1/),Dx=Da,Fx=Fa)

       ! One-electron contributions to E^{1,bc}*Da - S^{bc}*Wa
       call rsp_oneave(molcfg, S, (/MyRspFunc%code(2:3) /), (/Da/), &
            & shape(quad_RspFunc), quad_RspFunc, &
            & comp=(/MyRspFunc%first_comp(2:3)/), &
            & perm = (/2,3,1 /), &
            & freq=(/MyRspFunc%freq(2:3)/), DFD=(/Wa/) )

       ! Two-electron contributions to E^{1,bc}*Da
       call rsp_twoave(molcfg, (/ MyRspFunc%code(2:3) /), (/D, Da/), &
            & shape(quad_RspFunc), quad_RspFunc, &
            & perm = (/2,3,1 /), &
            & comp=(/MyRspFunc%first_comp(2:3)/) )

       if(isdef(Wa(1))) Wa(:)=0

    endif
    !************************************************************************
    ! E^{1,ac}*Db - S^{ac}*Wb contributions                                 !
    !************************************************************************
    ! Never for transition moments.

    if(.not. MyRspFunc%trans_moment) then

       ! Get Wb
       call get_W_matrix(F,D,S,MyRspFunc,1,pdbs,Wb,ops=(/2/),Dx=Db,Fx=Fb)

       ! One-electron contributions to E^{1,ac}*Db - S^{ac}*Wb
       call rsp_oneave(molcfg, S, (/MyRspFunc%code(1), MyRspFunc%code(3) /), (/Db/), &
            & shape(quad_RspFunc), quad_RspFunc, &
            & comp=(/MyRspFunc%first_comp(1), MyRspFunc%first_comp(3) /), &
            & perm = (/1,3,2 /), &
            & freq=(/MyRspFunc%freq(1), MyRspFunc%freq(3) /), DFD=(/Wb/) )

       ! Two-electron contributions to E^{1,ac}*Db
       call rsp_twoave(molcfg, (/ MyRspFunc%code(1), MyRspFunc%code(3) /), (/D, Db/), &
            & shape(quad_RspFunc), quad_RspFunc, &
            & perm = (/1,3,2 /), &
            & comp=(/MyRspFunc%first_comp(1), MyRspFunc%first_comp(3) /) )

       if(isdef(Wb(1))) Wb(:)=0

    endif
    !************************************************************************
    ! E^{1,ab}*Dc - S^{ab}*Wc contributions                                 !
    !************************************************************************
    ! May contribute to single transition moment

    if(any(pdbs(1:2)) .and. MyRspFunc%trans_moment_order /= 2) then

       ! Get Wc
       call get_W_matrix(F,D,S,MyRspFunc,1,pdbs,Wc,ops=(/3/),Dx=Dc,Fx=Fc)

       ! One-electron contributions to E^{1,ab}*Dc - S^{ab}*Wc
       call rsp_oneave(molcfg, S, (/MyRspFunc%code(1:2) /), (/Dc/), &
            & shape(quad_RspFunc), quad_RspFunc, &
            & comp=(/MyRspFunc%first_comp(1:2)/), &
            & perm = (/1,2,3 /), &
            & freq=(/MyRspFunc%freq(1:2)/), DFD=(/Wc/) )

       ! Two-electron contributions to E^{1,bc}*Da
       call rsp_twoave(molcfg, (/ MyRspFunc%code(1:2) /), (/D, Dc/), &
            & shape(quad_RspFunc), quad_RspFunc, &
            & perm = (/1,2,3 /), &
            & comp=(/MyRspFunc%first_comp(1:2)/) )

       if(isdef(Wc(1))) Wc(:)=0

    endif

    !************************************************************************
    ! E^{2,c} Da Db + E^{2,b} Da Dc + E^{2,a} Db Dc                         !
    !************************************************************************

    zeroD = 0*D
    zeroF = 0*F

    ! E^{2,c} Da Db, never for transition moments
    if(.not. MyRspFunc%trans_moment) then
       call rsp_twoave(molcfg, (/ MyRspFunc%code(3) /), (/D, Da, Db, &
            & ( zeroD, i=1,MyRspFunc%dims(1)*MyRspFunc%dims(2) ) /), &
            & shape(quad_RspFunc), quad_RspFunc, perm=(/3,1,2/), &
            & comp=(/MyRspFunc%first_comp(3)/) )
    endif

    ! E^{2,b} Da Dc, may contribute to single transition moment
    if(MyRspFunc%trans_moment_order /= 2) then
       call rsp_twoave(molcfg, (/ MyRspFunc%code(2) /), (/D, Da, Dc, &
            & ( zeroD, i=1,MyRspFunc%dims(1)*MyRspFunc%dims(3) ) /), &
            & shape(quad_RspFunc), quad_RspFunc, perm=(/2,1,3/), &
            & comp=(/MyRspFunc%first_comp(2)/) )
    endif

    ! E^{2,a} Db Dc, may contribute to single transition moment.
    ! This is the first term that may contribute to double transition moment
    ! since both Db and Dc occur.
    call rsp_twoave(molcfg, (/ MyRspFunc%code(1) /), (/D, Db, Dc, &
         & ( zeroD, i=1,MyRspFunc%dims(2)*MyRspFunc%dims(3) ) /), &
         & shape(quad_RspFunc), quad_RspFunc, perm=(/1,2,3/), &
         & comp=(/MyRspFunc%first_comp(1)/) )
    !************************************************************************
    ! E^{3,0} Da Db Dc (only DFT)                                           !
    !************************************************************************

    call rsp_twoave(molcfg, NOOP, (/D, Da, Db, Dc, &
         & ( zeroD, i=1,MyRspFunc%dims(1)*MyRspFunc%dims(2) ), &
         & ( zeroD, i=1,MyRspFunc%dims(1)*MyRspFunc%dims(3) ), &
         & ( zeroD, i=1,MyRspFunc%dims(2)*MyRspFunc%dims(3) ), &
         & ( zeroD, i=1,MyRspFunc%dims(1)*MyRspFunc%dims(2)*MyRspFunc%dims(3) ) /), &
         & shape(quad_RspFunc), quad_RspFunc )

    !************************************************************************
    ! -Sa W^{bc}_{1'}     -lambda_a * Y^{bc}_{1'}    - zeta_a * Z^{bc}_{1'}
    !************************************************************************

    ! Determine lambda and zeta multipliers
    do i=1,MyRspFunc%dims(1)
       lambda(i) = Da(i)*S*D - D*S*Da(i)
       zeta(i) = Fa(i)*D*S + S*D*Fa(i) - Fa(i) - F*D*Sa(i) - Sa(i)*D*F
    enddo

    ! loop over b and c components, one component at a time to save space.
    do j=1,MyRspFunc%dims(2)

       do k=1,MyRspFunc%dims(3)

          ! Construct Ybc and Zbc
          b = MyRspFunc%first_comp(2)+j-1
          c = MyRspFunc%first_comp(3)+k-1


          call pert_scf_eq(molcfg, S, (/ MyRspFunc%code(2), MyRspFunc%code(3) /), &
               & (/1,1/), (/D,Db(j),Dc(k), zeroD/), &
               & (/F,Fb(j),Fc(k), zeroF/), &
               & Ybc, Zbc, comp = (/b,c/), & 
               & freq=(/MyRspFunc%freq(2), MyRspFunc%freq(3)/) )

          if(pdbs(1)) then
             ! Get Wbc' 
             ! (prime indicates that only first order derivatives should be included in Wbc)
             call get_W_matrix(F,D,S,MyRspFunc,2,pdbs,Wbc,ops=(/2,3/),&
                  Dx=(/Db(j)/), Fx=(/Fb(j)/),Sx=(/Sb(j)/), &
                  Dy=(/Dc(k)/), Fy=(/Fc(k)/),Sy=(/Sc(k)/), &
                  Dxy= (/ ( zeroD ) /), &
                  Fxy=(/ ( zeroF ) /) )
          endif

          ! Calculate contributions.
          do i=1,MyRspFunc%dims(1)

             quad_RspFunc(i,j,k) = quad_RspFunc(i,j,k) -trace(lambda(i),Ybc(1))
             quad_RspFunc(i,j,k) = quad_RspFunc(i,j,k) -trace(zeta(i),Zbc(1))
             ! Only Sa Wbc contribution if perturbation A uses PDBS.
             if(pdbs(1)) then
                quad_RspFunc(i,j,k) = quad_RspFunc(i,j,k) -trace(Sa(i),Wbc(1))
             endif

          enddo

       enddo

    enddo
    ! Done calculating quadratic response function.
    ! Set rsp_results_vec equation to the entries in the quadratic response function
    ! using the index ordering described in the beginning of this subroutine.
    l=0
    do i=1,MyRspFunc%dims(1)
       do j=1,MyRspFunc%dims(2)
          do k=1,MyRspFunc%dims(3)
             l=l+1
             rsp_results_vec(l) = quad_RspFunc(i,j,k)
          enddo
       enddo
    enddo



    ! Free matrices
    Da(:)=0; Fa(:)=0; Db(:)=0; Fb(:)=0; Dc(:)=0; Fc(:)=0
    lambda(:)=0; zeta(:)=0

    if(isdef(Sa(1))) Sa(:)=0
    if(isdef(Sb(1))) Sb(:)=0
    if(isdef(Sc(1))) Sc(:)=0
    if(isdef(Ybc(1))) Ybc(1)=0
    if(isdef(Zbc(1))) Zbc(1)=0
    if(isdef(Wbc(1))) Wbc(1)=0

  end subroutine quadratic_response_2nplus1



  !> \brief Using the n+1 rule (Determine Db,Dc, and Dbc)
  !> this subroutine calculates quadratic response function or transition moments
  !> obtained from single or double residues of quadratic response functions
  !> according to the response information in MyRspFunc and saves
  !> the results in the rsp_results_vec vector.
  !> <br> <br>
  !> Details: 
  !> <br>
  !> We use the AO density matrix-based response formulation in J. Chem. Phys. 129, 214108,
  !> see in particular Eq. 237 for the quadratic response function using the n+1 rule.
  !> Single and double residues are found from Eq. (237) using the strategy presented in Sec. IV.H.
  !> <br>
  !> rsp_results_vec will contain all response results, possibly for several
  !> field components, but always for only one set of frequencies (or excitation energies).
  !> The arrangement of the results is best explained by an example.
  !> Consider <<A;B,C>> where the x,y, and z (1,2,3) components are specified
  !> for all operators.
  !> In this case the 3x3x3=27 components 
  !> in rsp_results_vec are arranged as follows:
  !> <br> <br>
  !> (1,1,1) (1,1,2) (1,1,3) (1,2,1) (1,2,2) (1,2,3) (1,3,1) (1,3,2) (1,3,3)
  !> <br>
  !> (2,1,1) (2,1,2) (2,1,3) (2,2,1) (2,2,2) (2,2,3) (2,3,1) (2,3,2) (2,3,3)
  !> <br>
  !> (3,1,1) (3,1,2) (3,1,3) (3,2,1) (3,2,2) (3,2,3) (3,3,1) (3,3,2) (3,3,3)
  !> <br> <br>
  !> where the first, second, and third indices refer to A, B, and C, respectively.
  !> <br> <br>
  !> NOTE! 
  !> <br>
  !> In case you want to calculate a quadratic transition moment by 
  !> setting up MyRspFunc in a separate subroutine (not in the input),
  !> and calling quadratic_response_nplus1:
  !> You must then first call transition_moment_density_matrix in response_driver_module
  !> to determine the transition density matrix, and then call linear_response.
  !> <br>
  !> On the other hand, in a quadratic response function calculation,
  !> the perturbed density matrix is calculated on the fly by quadratic_response_nplus1 itself,
  !> and you may obtain the a quadratic response function by a single call to this routine.
  !> \author Kasper Kristensen
  !> \date 2010-03
  subroutine quadratic_response_nplus1(molcfg,F,D,S,MyRspFunc, rsp_results_vec,response_size)
    implicit none
    !> References to molecule input, and solver and integral settings
    type(rsp_molcfg), intent(inout)   :: molcfg
    !> Unperturbed Fock matrix
    type(matrix), intent(in)          :: F
    !> Unperturbed density matrix
    type(matrix), intent(in)          :: D
    !> Unperturbed overlap matrix
    type(matrix), intent(in)          :: S
    !> Contains all information about the quadratic 
    !> response function or transition moment to be determined.
    type(RspFunc), intent(inout)      :: MyRspFunc
    !> Dimension of response vector
    integer, intent(in)               :: response_size
    !> Vector for collecting response results
    complex(realk), intent(inout)     :: rsp_results_vec(response_size)
    !> Local parameters
    type(matrix)                      :: Db(MyRspFunc%dims(2)), &
         & Fb(MyRspFunc%dims(2)), Wb(MyRspFunc%dims(2)), &
         & Dc(MyRspFunc%dims(3)), &
         & Fc(MyRspFunc%dims(3)), Wc(MyRspFunc%dims(3)), &
         & Dbc(MyRspFunc%dims(2),MyRspFunc%dims(3)), &
         & Fbc(MyRspFunc%dims(2),MyRspFunc%dims(3)), &
         & Wbc(MyRspFunc%dims(2),MyRspFunc%dims(3)), &
         & Sb(MyRspFunc%dims(2)), Sc(MyRspFunc%dims(3)), W(1)
    complex(realk)          :: quad_RspFunc(MyRspFunc%dims(1),MyRspFunc%dims(2), MyRspFunc%dims(3)), &
         & freqbc
    logical                 :: pdbs(3)
    integer                 :: i,j, k, l, b,c, singularity_component
    character*4             :: NOOP(1)
    integer, pointer        :: lupri
    NOOP(1) = 'NOOP'

    lupri => molcfg%lupri

    write(lupri, *) 'Starting quadratic response function calculation for: ', &
         & '<< ',MyRspFunc%code(1), '; ', MyRspFunc%code(2), ', ', MyRspFunc%code(3), '>>'

    write(lupri,*)

    ! omega_bc = omega_b + omega_c
    freqbc = MyRspFunc%freq(2)+MyRspFunc%freq(3)

    ! Determine if the operators in fields contain perturbation dependent basis sets.
    ! E.g. if MyRspFunc%code(1:2)= ('GEO','EL', 'EL') then pdbs=(.true. , .false., .false.).
    pdbs=pert_basdep(MyRspFunc%code(1:3))




    ! *******************************************************************************************
    !                    Calculate perturbed density matrices Db, Dc, and Dbc.                  !
    ! *******************************************************************************************

    ! Db matrix
    if(MyRspFunc%trans_moment .and. MyRspFunc%trans_moment_order ==2) then ! Double transition moment
       ! Db = -Dx^T
       ! where Dx contain the transition moment density matrices.
       ! Due to the transposition we effectively set omega_b = - excitation energy
       ! (see the transition_moment_density_matrix subroutine and the conventions for 
       ! residues in the beginning of the module and Sec. IV.H. if jcp, 129, 214108).
       Db(1) = -trans(rsp_eq_sol(MyRspFunc%ExNumber2)%D)
       Fb(1)=0E0_realk*F
       call rsp_twoint(molcfg, NOOP, (/D,Db(1)/), (/1/), Fb(1))
    else ! single residue or ordinary quadratic response function
       call pert_dens(molcfg, S, (/ MyRspFunc%code(2) /), (/MyRspFunc%dims(2)/), (/D/), (/F/), &
            & Db, Fb, comp=(/MyRspFunc%first_comp(2)/), freq=(/MyRspFunc%freq(2)/) )
    endif

    ! Dc matrix
    if(MyRspFunc%trans_moment) then     ! If perturbation c is EXCI
       ! Dc = Dx
       ! because we set omega_c = omega_m
       ! (see the transition_moment_density_matrix subroutine and the conventions for 
       ! residues in the beginning of the module and Sec. IV.H. if jcp, 129, 214108).
       Dc(1) = rsp_eq_sol(MyRspFunc%ExNumber1)%D
       ! Determine Fc
       Fc(1)=0E0_realk*F
       call rsp_twoint(molcfg, NOOP, (/D,Dc(1)/), (/1/), Fc(1))
    else ! ordinary quadratic response function
       call pert_dens(molcfg, S, (/ MyRspFunc%code(3) /), (/MyRspFunc%dims(3)/), (/D/), (/F/), &
            & Dc, Fc, comp=(/MyRspFunc%first_comp(3)/), freq=(/MyRspFunc%freq(3)/) )
    endif

    ! Dbc matrix
    call pert_dens(molcfg, S, (/ MyRspFunc%code(2), MyRspFunc%code(3) /), &
         & (/ MyRspFunc%dims(2), MyRspFunc%dims(3) /), (/D,Db,Dc/), (/F,Fb,Fc/), &
         & Dbc, Fbc, &
         & comp=(/MyRspFunc%first_comp(2), MyRspFunc%first_comp(3)/), &
         &  freq=(/MyRspFunc%freq(2), MyRspFunc%freq(3)/) )


    ! *******************************************************************************************
    !                    Calculate matrices for perturbation-dependent basis sets               !
    ! *******************************************************************************************

    ! Determine Sb
    if(pdbs(2)) then
       call rsp_oneint(molcfg, S, (/MyRspFunc%code(2) /), (/MyRspFunc%dims(2) /), &
            comp=(/MyRspFunc%first_comp(2)/), S=Sb)
    else
       do i=1,MyRspFunc%dims(2)
          Sb(i)=0E0_realk*S
       enddo
    endif

    ! Determine Sc
    ! Setting Sc=0 for trans_moment ensures that possible terms containing Sc
    ! do not contribute to the transition moment.
    if(pdbs(3) .and. .not. MyRspFunc%trans_moment) then
       call rsp_oneint(molcfg, S, (/MyRspFunc%code(3) /), (/MyRspFunc%dims(3) /), &
            comp=(/MyRspFunc%first_comp(3)/), S=Sc)
    else
       do i=1,MyRspFunc%dims(3)
          Sc(i)=0E0_realk*S
       enddo
    endif




    ! **********************************************************************
    ! *                START QUADRATIC RESPONSE CALCULATION                *
    ! **********************************************************************

    quad_RspFunc(:,:,:) = (0E0_realk,0E0_realk)


    !************************************************************************
    ! E^{0,abc} - S^{abc}*W contributions                                   !
    ! (only contributes for PDBS and never for transition moment)           !
    !************************************************************************

    ! One-electron contributions to E^{0,abc} - S^{abc}*W
    if(any(pdbs(2:3)) .and. .not. MyRspFunc%trans_moment) then

       ! Calculate zeroth order W
       call get_W_matrix(F,D,S,MyRspFunc,0,pdbs,W)

       call rsp_oneave(molcfg, S, MyRspFunc%code(1:3), (/D/), shape(quad_RspFunc), quad_RspFunc, &
            & comp= (/ MyRspFunc%first_comp(1:3) /), &
            & freq= (/ MyRspFunc%freq(1:3) /), DFD=(/W/) )

       ! Two-electron contributions to E^{0,abc}
       ! only if ALL perturbation use PDBS.
       if(all(pdbs)) then
          call rsp_twoave(molcfg, MyRspFunc%code(1:3), (/D/), shape(quad_RspFunc), quad_RspFunc, &
               & comp= (/ MyRspFunc%first_comp(1:3) /) )
       endif

       if(isdef(W(1))) W(1)=0

    endif


    !************************************************************************
    ! E^{1,ac}*Db - S^{ac}*Wb contributions                                 !
    !************************************************************************

    if( (pdbs(1) .or. pdbs(3)) .and. .not. MyRspFunc%trans_moment) then

       ! Get Wb
       call get_W_matrix(F,D,S,MyRspFunc,1,pdbs,Wb,ops=(/2/),Dx=Db,Fx=Fb)

       ! One-electron contributions to E^{1,ac}*Db - S^{ac}*Wb
       call rsp_oneave(molcfg, S, (/MyRspFunc%code(1), MyRspFunc%code(3) /), (/Db/), &
            & shape(quad_RspFunc), quad_RspFunc, &
            & comp=(/MyRspFunc%first_comp(1), MyRspFunc%first_comp(3) /), &
            & perm = (/1,3,2 /), &
            & freq=(/MyRspFunc%freq(1), MyRspFunc%freq(3) /), DFD=(/Wb/) )

       ! Two-electron contributions to E^{1,ac}*Db
       call rsp_twoave(molcfg, (/ MyRspFunc%code(1), MyRspFunc%code(3) /), (/D, Db/), &
            & shape(quad_RspFunc), quad_RspFunc, &
            & perm = (/1,3,2 /), &
            & comp=(/MyRspFunc%first_comp(1), MyRspFunc%first_comp(3) /) )

       if(isdef(Wb(1))) Wb(:)=0
    endif


    !************************************************************************
    ! E^{1,ab}*Dc - S^{ab}*Wc contributions                                 !
    !************************************************************************

    ! May contribute to single transition moment
    if(any(pdbs(1:2)) .and. MyRspFunc%trans_moment_order /= 2) then

       ! Get Wc
       call get_W_matrix(F,D,S,MyRspFunc,1,pdbs,Wc,ops=(/3/),Dx=Dc,Fx=Fc)

       ! One-electron contributions to E^{1,ab}*Dc - S^{ab}*Wc
       call rsp_oneave(molcfg, S, (/MyRspFunc%code(1:2) /), (/Dc/), &
            & shape(quad_RspFunc), quad_RspFunc, &
            & comp=(/MyRspFunc%first_comp(1:2)/), &
            & perm = (/1,2,3 /), &
            & freq=(/MyRspFunc%freq(1:2)/), DFD=(/Wc/) )

       ! Two-electron contributions to E^{1,ab}*Dc
       call rsp_twoave(molcfg, (/ MyRspFunc%code(1:2) /), (/D, Dc/), &
            & shape(quad_RspFunc), quad_RspFunc, &
            & perm = (/1,2,3 /), &
            & comp=(/MyRspFunc%first_comp(1:2)/) )

       if(isdef(Wc(1))) Wc(:)=0
    endif

    !************************************************************************
    ! E^{1,a} Dbc - Sa Wbc + E^{2,a} Db Dc                                  !
    !************************************************************************

    ! Get Wbc
    call get_W_matrix(F,D,S,MyRspFunc,2,pdbs,Wbc,ops=(/2,3/),Dx=Db,Fx=Fb,Sx=Sb,&
         Dy=Dc,Fy=Fc,Sy=Sc, Dxy=Dbc, Fxy=Fbc)

    ! One-electron part of E^{1,a} Dbc - Sa Wbc
    call rsp_oneave(molcfg, S, (/MyRspFunc%code(1) /), (/Dbc/), &
         & shape(quad_RspFunc), quad_RspFunc, &
         & comp=(/MyRspFunc%first_comp(1)/), &
         & perm = (/1,2,3 /), &
         & freq=(/MyRspFunc%freq(1)/), DFD=(/Wbc/) )

    if(isdef(Wbc(1,1))) Wbc(:,:)=0

    ! Two-electron part of E^{1,a} Dbc and  E^{2,a} Db Dc (the latter is purely two-electron)
    call rsp_twoave(molcfg, (/ MyRspFunc%code(1) /), (/D, Db,Dc,Dbc/), &
         & shape(quad_RspFunc), quad_RspFunc, &
         & perm = (/1,2,3 /), &
         & comp=(/MyRspFunc%first_comp(1)/) )


    ! Done calculating quadratic response function.
    ! Set rsp_results equal to quadratic response function based on the index ordering
    ! described in the beginning of the subroutine.
    l=0
    do i=1,MyRspFunc%dims(1)
       do j=1,MyRspFunc%dims(2)
          do k=1,MyRspFunc%dims(3)
             l=l+1
             rsp_results_vec(l) = quad_RspFunc(i,j,k)
          enddo
       enddo
    enddo


    ! Free matrices
    Db(:)=0; Fb(:)=0; Dc(:)=0; Fc(:)=0; Dbc(:,:)=0; Fbc(:,:)=0

    if(isdef(Sb(1))) Sb(:)=0
    if(isdef(Sc(1))) Sc(:)=0

  end subroutine quadratic_response_nplus1


  !> \brief Calculates the cubic response function using Eq. (259)
  !> in jcp, 129, 214108 (the 2n+1 and 2n+2 rules).
  !> <br>
  !> It is also possible to get transition moments. Single transition moments are
  !> defined according to Eqs. (3.30) and Eq. (3.31), respectively, in jcp 82, 3235 (1985).
  !> <br>
  !> Only Hartree-Fock is possible right now so DFT terms have been commented out.
  !> <br> <br>
  !> The arrangement of the cubic response results is best explained by an example.
  !> Consider <<A;B,C,D>> where the x,y, and z (1,2,3) components are specified
  !> for all operators.
  !> In this case the 3x3x3x3=81 components 
  !> in rsp_results_vec are arranged as follows (read from left to right):
  !> <br> <br>
  !> (1,1,1,1) (1,1,1,2) (1,1,1,3) (1,1,2,1) (1,1,2,2) (1,1,2,3) (1,1,3,1) (1,1,3,2) (1,1,3,3)
  !> <br>
  !> (1,2,1,1) (1,2,1,2) (1,2,1,3) (1,2,2,1) (1,2,2,2) (1,2,2,3) (1,2,3,1) (1,2,3,2) (1,2,3,3)
  !> <br>
  !> (1,3,1,1) (1,3,1,2) (1,3,1,3) (1,3,2,1) (1,3,2,2) (1,3,2,3) (1,3,3,1) (1,3,3,2) (1,3,3,3)
  !> <br>
  !> (2,1,1,1) (2,1,1,2) (2,1,1,3) (2,1,2,1) (2,1,2,2) (2,1,2,3) (2,1,3,1) (2,1,3,2) (2,1,3,3)
  !> <br>
  !> (2,2,1,1) (2,2,1,2) (2,2,1,3) (2,2,2,1) (2,2,2,2) (2,2,2,3) (2,2,3,1) (2,2,3,2) (2,2,3,3)
  !> <br>
  !> (2,3,1,1) (2,3,1,2) (2,3,1,3) (2,3,2,1) (2,3,2,2) (2,3,2,3) (2,3,3,1) (2,3,3,2) (2,3,3,3)
  !> <br>
  !> (3,1,1,1) (3,1,1,2) (3,1,1,3) (3,1,2,1) (3,1,2,2) (3,1,2,3) (3,1,3,1) (3,1,3,2) (3,1,3,3)
  !> <br>
  !> (3,2,1,1) (3,2,1,2) (3,2,1,3) (3,2,2,1) (3,2,2,2) (3,2,2,3) (3,2,3,1) (3,2,3,2) (3,2,3,3)
  !> <br>
  !> (3,3,1,1) (3,3,1,2) (3,3,1,3) (3,3,2,1) (3,3,2,2) (3,3,2,3) (3,3,3,1) (3,3,3,2) (3,3,3,3)
  !> <br> <br>
  !> where the first, second, and third indices refer to A, B, and C, respectively.
  !> <br> <br>
  !> NOTE!!
  !> <br>
  !> In case you want to calculate a cubic transition moment by 
  !> setting up MyRspFunc in a separate subroutine (not in the input),
  !> and calling cubic_response:
  !> You must then first call transition_moment_density_matrix in response_driver_module
  !> to determine the transition density matrix, and then call linear_response.
  !> <br>
  !> On the other hand, in a cubic response function calculation,
  !> the perturbed density matrix is calculated on the fly by cubic_response itself,
  !> and you may obtain cubic response function by a single call to this routine.
  !> \author Kasper Kristensen
  !> \date 2010-03
  subroutine cubic_response(molcfg,F,D,S,MyRspFunc, rsp_results_vec,response_size)

    implicit none
    !> References to molecule input, and solver and integral settings
    type(rsp_molcfg), intent(inout)   :: molcfg
    !> Unperturbed Fock matrix
    type(matrix), intent(in)          :: F
    !> Unperturbed density matrix
    type(matrix), intent(in)          :: D
    !> Unperturbed overlap matrix
    type(matrix), intent(in)          :: S
    !> Contains all information about the cubic
    !> response function or transition moment to be determined.
    type(RspFunc), intent(inout)      :: MyRspFunc
    !> Dimension of response vector containing response function
    integer, intent(in)               :: response_size
    !> Vector for collecting response results
    complex(realk), intent(inout)     :: rsp_results_vec(response_size)
    !> Local parameters (many of the matrices are not used for most properties)
    type(matrix)  :: Da(MyRspFunc%dims(1)), &
         & Fa(MyRspFunc%dims(1)), Wa(MyRspFunc%dims(1)), &
         & Db(MyRspFunc%dims(2)), &
         & Fb(MyRspFunc%dims(2)), Wb(MyRspFunc%dims(2)), &
         & Dc(MyRspFunc%dims(3)), &
         & Fc(MyRspFunc%dims(3)), Wc(MyRspFunc%dims(3)), &
         & Fd(MyRspFunc%dims(4)), Wd(MyRspFunc%dims(4)), &
         & Dd(MyRspFunc%dims(4)), &
         & Dab(MyRspFunc%dims(1),MyRspFunc%dims(2)), &
         & Fab(MyRspFunc%dims(1),MyRspFunc%dims(2)), &
         & Dac(MyRspFunc%dims(1),MyRspFunc%dims(3)), &
         & Fac(MyRspFunc%dims(1),MyRspFunc%dims(3)), &
         & Dad(MyRspFunc%dims(1),MyRspFunc%dims(4)), &
         & Fad(MyRspFunc%dims(1),MyRspFunc%dims(4)), &
         & Sa(MyRspFunc%dims(1)), Sb(MyRspFunc%dims(2)), &
         & Sc(MyRspFunc%dims(3)), Sd(MyRspFunc%dims(4)), &
         & Sab(MyRspFunc%dims(1),MyRspFunc%dims(2)), &
         & Sac(MyRspFunc%dims(1),MyRspFunc%dims(3)), &
         & Sad(MyRspFunc%dims(1),MyRspFunc%dims(4)), &
         & Wab(MyRspFunc%dims(1), MyRspFunc%dims(2)), &
         & Wac(MyRspFunc%dims(1), MyRspFunc%dims(3)), &
         & Wad(MyRspFunc%dims(1), MyRspFunc%dims(4)), &
         & Wbc(MyRspFunc%dims(2), MyRspFunc%dims(3)), &
         & Wbd(MyRspFunc%dims(2), MyRspFunc%dims(4)), &
         & Wcd(MyRspFunc%dims(3), MyRspFunc%dims(4)), &
         & W(1), zeroD, zeroF
    complex(realk)                 :: cubic_RspFunc(MyRspFunc%dims(1),MyRspFunc%dims(2), &
         & MyRspFunc%dims(3), MyRspFunc%dims(4))
    logical                        :: pdbs(4), AeqB, AeqC, AeqD, BeqC, BeqD, CeqD
    integer                        :: i,j, k, l, m
    character(4)                   :: NOOP(1)
    complex(realk) :: freqA,freqB,freqC,freqD
    integer :: nA, nB, nC, nD, first_compA,first_compB, first_compC, first_compD
    ! KK HACK: Quick fix of problem that "MAG, nd > 2 not implemented" in rsp_contribs.F90
    complex(realk), allocatable :: hack_bcad(:,:,:), hack_bdac(:,:,:), &
         & hack_cbad(:,:,:), hack_cdab(:,:,:), &
         & hack_dbac(:,:,:), hack_dcab(:,:,:)
    integer :: count,ia,ib,ic,id
    integer, pointer :: lupri
    NOOP(1) = 'NOOP'

    lupri => molcfg%lupri
    i=0

    ! Initialize stuff
    nA=MyRspFunc%dims(1)
    nB=MyRspFunc%dims(2)
    nC=MyRspFunc%dims(3)
    nD=MyRspFunc%dims(4)
    freqA = MyRspFunc%freq(1)
    freqB = MyRspFunc%freq(2)
    freqC = MyRspFunc%freq(3)
    freqD = MyRspFunc%freq(4)
    first_compA = MyRspFunc%first_comp(1)
    first_compB = MyRspFunc%first_comp(2)
    first_compC = MyRspFunc%first_comp(3)
    first_compD = MyRspFunc%first_comp(4)

    ! KK HACK: Quick fix of problem that "MAG, nd > 2 not implemented" in rsp_contribs.F90
    allocate(hack_bcad(nB,nC,nA*nD))
    allocate(hack_bdac(nB,nD,nA*nC))
    allocate(hack_cbad(nC,nB,nA*nD))
    allocate(hack_cdab(nC,nD,nA*nB))
    allocate(hack_dbac(nD,nB,nA*nC))
    allocate(hack_dcab(nD,nC,nA*nB))
    hack_bcad(:,:,:) = (0E0_realk,0E0_realk)
    hack_bdac(:,:,:) = (0E0_realk,0E0_realk)
    hack_cbad(:,:,:) = (0E0_realk,0E0_realk)
    hack_cdab(:,:,:) = (0E0_realk,0E0_realk)
    hack_dbac(:,:,:) = (0E0_realk,0E0_realk)
    hack_dcab(:,:,:) = (0E0_realk,0E0_realk)

    write(lupri,*)
    write(lupri, *) 'Starting cubic response function calculation for: ', &
         & '<< ',MyRspFunc%code(1), '; ', MyRspFunc%code(2), ', ', &
         & MyRspFunc%code(3), ', ', MyRspFunc%code(4), '>>'
    write(lupri,*)

    ! Determine if the operators in fields contain perturbation dependent basis sets.
    ! E.g. if MyRspFunc%code(1:4)= ('GEO','EL', 'EL', 'MAG') 
    ! then pdbs=(.true. , .false., .false., .true.).
    pdbs=pert_basdep(MyRspFunc%code(1:4))

    ! Check whether any of the perturbations are identical
    ! If so, we avoid re-calculating various matrices and dot products.
    AeqB = Identical_operators(MyRspFunc,1,2)
    AeqC = Identical_operators(MyRspFunc,1,3)
    AeqD = Identical_operators(MyRspFunc,1,4)
    BeqC = Identical_operators(MyRspFunc,2,3)
    BeqD = Identical_operators(MyRspFunc,2,4)
    CeqD = Identical_operators(MyRspFunc,3,4)


    ! *******************************************************************************************
    !                Calculate perturbed density matrices and Fock matrices                     !
    ! *******************************************************************************************
    call cubic_response_dens_and_fock(molcfg,F,D,S,MyRspFunc,&
         & Da,Fa,Db,Fb,Dc,Fc,Fd,Dd,Dab,Fab,Dac,Fac,Dad,Fad)



    ! *******************************************************************************************
    !               Calculate overlap matrices for perturbation-dependent basis sets            !
    ! *******************************************************************************************
    call cubic_response_overlap(molcfg,F,D,S,MyRspFunc,&
         & Sa,Sb,Sc,Sd,Sab,Sac,Sad)




    ! **********************************************************************
    !                    START CUBIC RESPONSE CALCULATION                  *
    ! **********************************************************************

    cubic_RspFunc(:,:,:,:) = (0E0_realk,0E0_realk)


    !************************************************************************
    !                  E^{0,abcd} - S^{abcd}*W contributions                !
    !************************************************************************
    ! Only for PDBS
    if( any(pdbs) .and. (.not. MyRspFunc%trans_moment) ) then

       ! Calculate zeroth order W
       call get_W_matrix(F,D,S,MyRspFunc,0,pdbs,W)

       ! One-electron contributions to E^{0,abcd} - S^{abcd}*W
       call rsp_oneave(molcfg, S, MyRspFunc%code(1:4), (/D/), shape(cubic_RspFunc), cubic_RspFunc, &
            & comp= (/ MyRspFunc%first_comp(1:4) /), &
            & freq= (/ MyRspFunc%freq(1:4) /), DFD=(/W/) )

       ! Two-electron contributions to E^{0,abcd},
       ! only if ALL perturbations use PDBS.
       if(all(pdbs)) then
          call rsp_twoave(molcfg, MyRspFunc%code(1:4), (/D/), shape(cubic_RspFunc), cubic_RspFunc, &
               & comp= (/ MyRspFunc%first_comp(1:4) /) )
       endif

       ! Free W
       if(isdef(W(1))) W(1)=0

    endif


    !************************************************************************
    !               E^{1,bcd}*Da - S^{bcd}*Wa contributions                 !
    !************************************************************************
    ! Only for PDBS and never if pert d = EXCI
    if( any(pdbs(2:4)) .and. MyRspFunc%code(4)/='EXCI' ) then

       ! Get Wa
       call get_W_matrix(F,D,S,MyRspFunc,1,pdbs,Wa,ops=(/1/),Dx=Da,Fx=Fa)

       ! One-electron contributions to E^{1,bcd}*Da - S^{bcd}*Wa
       call rsp_oneave(molcfg, S, (/MyRspFunc%code(2:4) /), (/Da/), &
            & shape(cubic_RspFunc), cubic_RspFunc, &
            & comp=(/first_compB,first_compC,first_compD/), &
            & perm = (/2,3,4,1 /), &   ! Order b-c-d-a in E^{1,bcd}*Da
            & freq=(/freqB,freqC,freqD/), DFD=(/Wa/) )

       ! Two-electron contributions to E^{1,bcd}*Da
       call rsp_twoave(molcfg, (/ MyRspFunc%code(2:4) /), (/D, Da/), &
            & shape(cubic_RspFunc), cubic_RspFunc, &
            & perm = (/2,3,4,1 /), &
            & comp=(/first_compB,first_compC,first_compD/) )

       ! Free Wa
       if(isdef(Wa(1))) Wa(:)=0

    end if


    !************************************************************************
    !               E^{1,acd}*Db - S^{acd}*Wb contributions                 !
    !************************************************************************
    ! Only for PDBS and never if pert d = EXCI
    if( (pdbs(1) .or. any(pdbs(3:4))) .and. MyRspFunc%code(4)/='EXCI' ) then

       ! Get Wb
       call get_W_matrix(F,D,S,MyRspFunc,1,pdbs,Wb,ops=(/2/),Dx=Db,Fx=Fb)

       ! One-electron contributions to E^{1,acd}*Db - S^{acd}*Wb
       call rsp_oneave(molcfg, S, (/MyRspFunc%code(1), MyRspFunc%code(3:4) /), (/Db/), &
            & shape(cubic_RspFunc), cubic_RspFunc, &
            & comp=(/first_compA, first_compC, first_compD /), &
            & perm = (/1,3,4,2 /), &   ! Order a-c-d-b in E^{1,acd}*Db
            & freq=(/freqA, freqC, freqD /), DFD=(/Wb/) )

       ! Two-electron contributions to E^{1,acd}*Db
       call rsp_twoave(molcfg, (/MyRspFunc%code(1), MyRspFunc%code(3:4) /), (/D, Db/), &
            & shape(cubic_RspFunc), cubic_RspFunc, &
            & perm = (/1,3,4,2 /), &
            & comp=(/first_compA, first_compC, first_compD /) )

       ! Free Wb
       if(isdef(Wb(1))) Wb(:)=0

    end if


    !************************************************************************
    !               E^{1,abd}*Dc - S^{abd}*Wc contributions                 !
    !************************************************************************
    ! Only for PDBS and never if pert d = EXCI
    if( (any(pdbs(1:2)) .or. pdbs(4)) .and. MyRspFunc%code(4)/='EXCI' ) then

       ! Get Wc
       call get_W_matrix(F,D,S,MyRspFunc,1,pdbs,Wc,ops=(/3/),Dx=Dc,Fx=Fc)

       ! One-electron contributions to E^{1,abd}*Dc - S^{abd}*Wc
       call rsp_oneave(molcfg, S, (/MyRspFunc%code(1:2), MyRspFunc%code(4) /), (/Dc/), &
            & shape(cubic_RspFunc), cubic_RspFunc, &
            & comp=(/ first_compA, first_compB, first_compD /), &
            & perm = (/1,2,4,3 /), &   ! Order a-b-d-c in E^{1,abd}*Dc
            & freq=(/freqA, freqB, freqD /), DFD= (/Wc/) )

       ! Two-electron contributions to E^{1,abd}*Dc
       call rsp_twoave(molcfg, (/MyRspFunc%code(1:2), MyRspFunc%code(4) /), (/D, Dc/), &
            & shape(cubic_RspFunc), cubic_RspFunc, &
            & perm = (/1,2,4,3 /), &
            & comp=(/first_compA, first_compB, first_compD /) )

       ! Free Wc
       if(isdef(Wc(1))) Wc(:)=0

    end if


    !************************************************************************
    !               E^{1,abc}*Dd - S^{abc}*Wd contributions                 !
    !************************************************************************
    ! Only for PDBS and never for double transition moment
    if(any(pdbs(1:3)) .and. MyRspFunc%trans_moment_order /= 2) then

       ! Get Wd
       call get_W_matrix(F,D,S,MyRspFunc,1,pdbs,Wd,ops=(/4/),Dx=Dd,Fx=Fd)

       ! One-electron contributions to E^{1,abc}*Dd - S^{abc}*Wd
       call rsp_oneave(molcfg, S, (/MyRspFunc%code(1:3) /), (/Dd/), &
            & shape(cubic_RspFunc), cubic_RspFunc, &
            & comp=(/ first_compA, first_compB, first_compC/), &
            & perm = (/1,2,3,4 /), &   ! Order a-b-c-d in E^{1,abc}*Dd
            & freq=(/freqA, freqB, freqC/), DFD=(/Wd/) )

       ! Two-electron contributions to E^{1,abc}*Dd
       call rsp_twoave(molcfg, (/ MyRspFunc%code(1:3) /), (/D, Dd/), &
            & shape(cubic_RspFunc), cubic_RspFunc, &
            & perm = (/1,2,3,4 /), &
            & comp=(/ first_compA, first_compB, first_compC/) )

       ! Free Wd
       if(isdef(Wd(1))) Wd(:)=0

    end if


    !************************************************************************
    !       E^{1,cd}*Dab + E^{2,cd} Da Db - S^{cd}*Wab contributions        !
    !************************************************************************
    ! Only for PDBS or terms quadratic in the field (e.g. magnetic)
    ! And not if pert d = EXCI

    if(MyRspFunc%code(4) /= 'EXCI' ) then

       ! Get Wab (only if c and d use PDBS)
       call get_W_matrix(F,D,S,MyRspFunc,2,pdbs,Wab,ops=(/1,2/),Dx=Da,Fx=Fa,Sx=Sa,&
            Dy=Db,Fy=Fb,Sy=Sb, Dxy=Dab, Fxy=Fab)

       ! One-electron contributions
       call rsp_oneave(molcfg, S, (/MyRspFunc%code(3:4) /), (/Dab/), &
            & shape(cubic_RspFunc), cubic_RspFunc, &
            & comp=(/first_compC, first_compD/), &
            & perm = (/3,4,1,2 /), &   ! Order c-d-a-b in E^{1,cd}*Dab
            & freq=(/freqC,freqD/), DFD=(/Wab/) )


       ! Two-electron contributions
       call rsp_twoave(molcfg, (/ MyRspFunc%code(3:4) /), (/D, Da, Db, Dab/), &
            & shape(cubic_RspFunc), cubic_RspFunc, &
            & perm = (/3,4,1,2 /), &
            & comp=(/first_compC, first_compD/) )

       ! Free Wab
       if(isdef(Wab(1,1))) Wab(:,:)=0

    end if


    !************************************************************************
    !       E^{1,bd}*Dac + E^{2,bd} Da Dc - S^{bd}*Wac contributions        !
    !************************************************************************
    ! Only for PDBS or terms quadratic in the field

    ! And not if pert d = EXCI
    if(MyRspFunc%code(4) /= 'EXCI' ) then

       ! Get Wac
       call get_W_matrix(F,D,S,MyRspFunc,2,pdbs,Wac,ops=(/1,3/),Dx=Da,Fx=Fa,Sx=Sa,&
            Dy=Dc,Fy=Fc,Sy=Sc, Dxy=Dac, Fxy=Fac)

       ! One-electron contributions
       call rsp_oneave(molcfg, S, (/MyRspFunc%code(2), MyRspFunc%code(4) /), (/Dac/), &
            & shape(cubic_RspFunc), cubic_RspFunc, &
            & comp=(/first_compB, first_compD/), &
            & perm = (/2,4,1,3 /), &   ! Order b-d-a-c in E^{1,bd}*Dac
            & freq=(/freqB, freqD/), DFD=(/Wac/) )

       ! Two-electron contributions
       call rsp_twoave(molcfg, (/MyRspFunc%code(2), MyRspFunc%code(4) /), &
            (/D, Da, Dc, Dac/), &
            & shape(cubic_RspFunc), cubic_RspFunc, &
            & perm = (/2,4,1,3 /), &
            & comp=(/first_compB, first_compD/) )

       ! Free Wac
       if(isdef(Wac(1,1))) Wac(:,:)=0

    end if


    !************************************************************************
    !       E^{1,bc}*Dad + E^{2,bc} Da Dd - S^{bc}*Wad contributions        !
    !************************************************************************
    ! Only for PDBS or terms quadratic in the field and never for double transition moment.

    if(MyRspFunc%trans_moment_order /= 2) then

       ! Get Wad
       call get_W_matrix(F,D,S,MyRspFunc,2,pdbs,Wad,ops=(/1,4/),Dx=Da,Fx=Fa,Sx=Sa,&
            Dy=Dd,Fy=Fd,Sy=Sd, Dxy=Dad, Fxy=Fad)

       call rsp_oneave(molcfg, S, (/MyRspFunc%code(2:3) /), (/Dad/), &
            & shape(cubic_RspFunc), cubic_RspFunc, &
            & comp=(/first_compB, first_compC/), &
            & perm = (/2,3,1,4 /), &   ! Order b-c-a-d in E^{1,bc}*Dad
            & freq=(/freqB, freqC/), DFD=(/Wad/) )

       ! Two-electron contributions
       call rsp_twoave(molcfg, (/ MyRspFunc%code(2:3) /), (/D, Da, Dd, Dad/), &
            & shape(cubic_RspFunc), cubic_RspFunc, &
            & perm = (/2,3,1,4 /), &
            & comp=(/first_compB, first_compC/) )

       ! Free Wad
       if(isdef(Wad(1,1))) Wad(:,:)=0

    end if


    !************************************************************************
    !                E^{2,ab} Dc Dd - S^{ab}*Wcd' contributions             !
    !************************************************************************
    ! Only for PDBS or terms quadratic in the field

    zeroD = 0*D
    zeroF = 0*F

    ! Get Wcd' 
    ! (prime indicates that only first order derivatives should be included in Wcd)
    call get_W_matrix(F,D,S,MyRspFunc,2,pdbs,Wcd,ops=(/3,4/),Dx=Dc,Fx=Fc,Sx=Sc,&
         Dy=Dd,Fy=Fd,Sy=Sd, &
         Dxy= (/ ( zeroD, i=1,nC*nD ) /), &
         Fxy=(/ ( zeroD, i=1,nC*nD ) /) )

    ! One-electron contributions
    ! In this case, -Sab*Wcd is the only 'one-electron' contribution
    call rsp_oneave(molcfg, S, (/MyRspFunc%code(1:2) /), &
         & (/ ( zeroD, i=1,nC*nD ) /), &
         & shape(cubic_RspFunc), cubic_RspFunc, &
         & comp=(/first_compA, first_compB/), &
         & perm = (/1,2,3,4 /), &  
         & freq=(/freqA, freqB/), DFD=(/Wcd/) )

    ! Two-electron contributions
    call rsp_twoave(molcfg, (/ MyRspFunc%code(1:2) /), &
         & (/ D, Dc, Dd, (zeroD, i=1,nC*nD) /), &
         & shape(cubic_RspFunc), cubic_RspFunc, &
         & perm = (/1,2,3,4 /), &
         & comp=(/first_compA, first_compB/) )

    ! Free Wcd'
    if(isdef(Wcd(1,1))) Wcd(:,:)=0


    !************************************************************************
    !                E^{2,ac} Db Dd - S^{ac}*Wbd' contributions             !
    !************************************************************************
    ! Only for PDBS or terms quadratic in the field and never for double transition moment

    if(MyRspFunc%trans_moment_order /= 2) then

       ! Get Wbd' 
       ! (prime indicates that only first order derivatives should be included in Wbd)
       call get_W_matrix(F,D,S,MyRspFunc,2,pdbs,Wbd,ops=(/2,4/),Dx=Db,Fx=Fb,Sx=Sb,&
            Dy=Dd,Fy=Fd,Sy=Sd, &
            Dxy= (/ ( zeroD, i=1,nB*nD ) /), &
            Fxy=(/ ( zeroF, i=1,nB*nD ) /) )

       ! One-electron contributions
       ! In this case, -Sac*Wbd is the only 'one-electron' contribution
       call rsp_oneave(molcfg, S, (/MyRspFunc%code(1),MyRspFunc%code(3) /), &
            & (/ ( zeroD, i=1,nB*nD ) /), &
            & shape(cubic_RspFunc), cubic_RspFunc, &
            & comp=(/first_compA,first_compC/), &
            & perm = (/1,3,2,4 /), &  
            & freq=(/freqA,freqC/), DFD=(/Wbd/) )

       ! Two-electron contributions
       call rsp_twoave(molcfg, (/ MyRspFunc%code(1),MyRspFunc%code(3) /), &
            & (/ D, Db, Dd, (zeroD, i=1,nB*nD) /), &
            & shape(cubic_RspFunc), cubic_RspFunc, &
            & perm = (/1,3,2,4 /), &
            & comp=(/first_compA,first_compC/) )

       ! Free Wbd'
       if(isdef(Wbd(1,1))) Wbd(:,:)=0

    end if


    !************************************************************************
    !                E^{2,ad} Db Dc - S^{ad}*Wbc' contributions             !
    !************************************************************************
    ! Only for PDBS or terms quadratic in the field and never for double transition moment

    if(MyRspFunc%code(4) /= 'EXCI') then

       ! Get Wbc' 
       ! (prime indicates that only first order derivatives should be included in Wbc)
       call get_W_matrix(F,D,S,MyRspFunc,2,pdbs,Wbc,ops=(/2,3/),Dx=Db,Fx=Fb,Sx=Sb,&
            Dy=Dc,Fy=Fc,Sy=Sc, &
            Dxy= (/ ( zeroD, i=1,nB*nC ) /), &
            Fxy=(/ ( zeroF, i=1,nB*nC ) /) )

       ! One-electron contributions
       ! In this case, -Sad*Wbc is the only 'one-electron' contribution
       call rsp_oneave(molcfg, S, (/MyRspFunc%code(1),MyRspFunc%code(4) /), &
            & (/ ( zeroD, i=1,nB*nC ) /), &
            & shape(cubic_RspFunc), cubic_RspFunc, &
            & comp=(/first_compA,first_compD/), &
            & perm = (/1,4,2,3 /), &  
            & freq=(/freqA,freqD/), DFD=(/Wbc/) )

       ! Two-electron contributions
       call rsp_twoave(molcfg, (/ MyRspFunc%code(1),MyRspFunc%code(4) /), &
            & (/ D, Db, Dc, (zeroD, i=1,nB*nC) /), &
            & shape(cubic_RspFunc), cubic_RspFunc, &
            & perm = (/1,4,2,3 /), &
            & comp=(/first_compA,first_compD/) )

       ! Free Wbc'
       if(isdef(Wbc(1,1))) Wbc(:,:)=0

    end if


    !************************************************************************
    !             E^{2,b} Dac Dd + E^{2,b} Dad Dc contributions             !
    !************************************************************************
    ! Only if b uses PDBS
    if(pdbs(2)) then

       ! KK HACK: Circumvent problem that "MAG, nd > 2 not implmented" in rsp_twoave.
       ! Temporary solution: For E^{2,b} Dac Dd = E^{2,b} Dd Dac
       ! we consider Dac and Dd as first order matrices!
       ! (And similarly for the other terms with abcd permuted)

       ! E^{2,b} Dd Dac
       call rsp_twoave(molcfg, (/MyRspFunc%code(2) /), &
            & (/D,Dd,Dac, &
            & (zeroD, i=1,nD*nA*nC) /), &
            & shape(hack_bdac), hack_bdac, &
            & comp=(/first_compB /) )

       ! Add hack_bdac to cubic response function.
       count=0
       do ic=1,nC
          do ia=1,nA
             count=count+1
             do id=1,nD
                do ib=1,nB
                   cubic_RspFunc(ia,ib,ic,id) = cubic_RspFunc(ia,ib,ic,id) &
                        + hack_bdac(ib,id,count)
                enddo
             end do
          end do
       end do

       ! E^{2,b} Dc Dad
       call rsp_twoave(molcfg, (/MyRspFunc%code(2) /), &
            & (/D,Dc,Dad, &
            & (zeroD, i=1,nC*nA*nD) /), &
            & shape(hack_bcad), hack_bcad, &
            & comp=(/first_compB /) )

       ! Add hack_bcad to cubic response function.
       count=0
       do id=1,nD
          do ia=1,nA
             count=count+1
             do ic=1,nC
                do ib=1,nB
                   cubic_RspFunc(ia,ib,ic,id) = cubic_RspFunc(ia,ib,ic,id) &
                        + hack_bcad(ib,ic,count)
                enddo
             end do
          end do
       end do

    end if


    !************************************************************************
    !             E^{2,c} Dad Db + E^{2,c} Dab Dd contributions             !
    !************************************************************************
    ! Only if c uses PDBS and not if c=EXCI
    if(pdbs(3) .and. MyRspFunc%trans_moment_order/=2 ) then

       ! E^{2,c} Db Dad
       call rsp_twoave(molcfg, (/MyRspFunc%code(3) /), &
            & (/D,Db,Dad, &
            & (zeroD, i=1,nB*nA*nD) /), &
            & shape(hack_cbad), hack_cbad, &
            & comp=(/first_compC /) )

       ! Add hack_cbad to cubic response function.
       count=0
       do id=1,nD
          do ia=1,nA
             count=count+1
             do ib=1,nB
                do ic=1,nC
                   cubic_RspFunc(ia,ib,ic,id) = cubic_RspFunc(ia,ib,ic,id) &
                        + hack_cbad(ic,ib,count)
                enddo
             end do
          end do
       end do

       ! E^{2,c} Dd Dab
       call rsp_twoave(molcfg, (/MyRspFunc%code(3) /), &
            & (/D,Dd,Dab, &
            & (zeroD, i=1,nD*nA*nB) /), &
            & shape(hack_cdab), hack_cdab, &
            & comp=(/first_compC /) )

       ! Add hack_cdab to cubic response function.
       count=0
       do ib=1,nB
          do ia=1,nA
             count=count+1
             do id=1,nD
                do ic=1,nC
                   cubic_RspFunc(ia,ib,ic,id) = cubic_RspFunc(ia,ib,ic,id) &
                        + hack_cdab(ic,id,count)
                enddo
             end do
          end do
       end do

    end if


    !************************************************************************
    !             E^{2,d} Dab Dc + E^{2,d} Dac Db contributions             !
    !************************************************************************
    ! Only if d uses PDBS and never if pert d = EXCI.
    if(pdbs(4) .and. MyRspFunc%code(4) /= 'EXCI') then

       ! E^{2,d} Dc Dab
       call rsp_twoave(molcfg, (/MyRspFunc%code(4) /), &
            & (/D,Dc,Dab, &
            & (zeroD, i=1,nC*nA*nB) /), &
            & shape(hack_dcab), hack_dcab, &
            & comp=(/first_compD /) )

       ! Add hack_dcab to cubic response function.
       count=0
       do ib=1,nB
          do ia=1,nA
             count=count+1
             do ic=1,nC
                do id=1,nD
                   cubic_RspFunc(ia,ib,ic,id) = cubic_RspFunc(ia,ib,ic,id) &
                        + hack_dcab(id,ic,count)
                enddo
             end do
          end do
       end do

       ! E^{2,d} Db Dac
       call rsp_twoave(molcfg, (/MyRspFunc%code(4) /), &
            & (/D,Db,Dac, &
            & (zeroD, i=1,nB*nA*nC) /), &
            & shape(hack_dbac), hack_dbac, &
            & comp=(/first_compD /) )


       ! Add hack_dbac to cubic response function.
       count=0
       do ic=1,nC
          do ia=1,nA
             count=count+1
             do ib=1,nB
                do id=1,nD
                   cubic_RspFunc(ia,ib,ic,id) = cubic_RspFunc(ia,ib,ic,id) &
                        + hack_dbac(id,ib,count)
                enddo
             end do
          end do
       end do

    end if


    !************************************************************************
    !                         E^{3,a} Db Dc Dd                              !
    !************************************************************************
    ! Only if a uses PDBS and only for DFT
    ! DFT stuff is not working right now and therefore it is commented out.
!!$      if(pdbs(1)) then
!!$
!!$         call rsp_twoave(molcfg, (/MyRspFunc%code(1) /), (/D,Db,Dc,Dd,&
!!$              & (zeroD, i=1,nB*nC), &
!!$              & (zeroD, i=1,nB*nD), &
!!$              & (zeroD, i=1,nC*nD), &
!!$              & (zeroD, i=1,nB*nC*nD) /), &
!!$              & shape(cubic_RspFunc), cubic_RspFunc, &
!!$              & perm = (/1,2,3,4/), &
!!$              & comp=(/first_compA /) )
!!$
!!$      end if



    !************************************************************************
    !                         E^{3,b} Da Dc Dd                              !
    !************************************************************************
    ! Only if b uses PDBS
    ! DFT stuff is not working right now and therefore it is commented out.
!!$      if(pdbs(2)) then
!!$
!!$         call rsp_twoave(molcfg, (/MyRspFunc%code(2) /), (/D,Da,Dc,Dd,&
!!$              & (zeroD, i=1,nA*nC), &
!!$              & (zeroD, i=1,nA*nD), &
!!$              & (zeroD, i=1,nC*nD), &
!!$              & (zeroD, i=1,nA*nC*nD) /), &
!!$              & shape(cubic_RspFunc), cubic_RspFunc, &
!!$              & perm = (/2,1,3,4/), &
!!$              & comp=(/first_compB /) )
!!$
!!$      end if


    !************************************************************************
    !                         E^{3,c} Da Db Dd                              !
    !************************************************************************
    ! Only if c uses PDBS
    ! DFT stuff is not working right now and therefore it is commented out.
!!$      if(pdbs(3) .and. MyRspFunc%trans_moment_order/=2) then
!!$
!!$         call rsp_twoave(molcfg, (/MyRspFunc%code(3) /), (/D,Da,Db,Dd,&
!!$              & (zeroD, i=1,nA*nB), &
!!$              & (zeroD, i=1,nA*nD), &
!!$              & (zeroD, i=1,nB*nD), &
!!$              & (zeroD, i=1,nA*nB*nD) /), &
!!$              & shape(cubic_RspFunc), cubic_RspFunc, &
!!$              & perm = (/3,1,2,4/), &
!!$              & comp=(/first_compC /) )
!!$
!!$      end if


    !************************************************************************
    !                         E^{3,d} Da Db Dc                              !
    !************************************************************************
    ! Only if d uses PDBS
    ! DFT stuff is not working right now and therefore it is commented out.
!!$      if(pdbs(4) .and. MyRspFunc%code(4) /= 'EXCI') then
!!$
!!$         call rsp_twoave(molcfg, (/MyRspFunc%code(4) /), (/D,Da,Db,Dc,&
!!$              & (zeroD, i=1,nA*nB), &
!!$              & (zeroD, i=1,nA*nC), &
!!$              & (zeroD, i=1,nB*nC), &
!!$              & (zeroD, i=1,nA*nB*nC) /), &
!!$              & shape(cubic_RspFunc), cubic_RspFunc, &
!!$              & perm = (/4,1,2,3/), &
!!$              & comp=(/first_compD /) )
!!$
!!$      end if


    !************************************************************************
    !          E^{3,0} Dab Dc Dd  + E^{3,0} Dac Db Dd                       !
    !        + E^{3,0} Dad Db Dc  + E^{4,0} Da Db Dc Dd      (only DFT)     !
    !************************************************************************
    ! DFT stuff is not working right now and therefore it is commented out.
    call rsp_twoave(molcfg, (/NOOP/), (/D,Da,Db,Dc,Dd,Dab,Dac, &
         & ( zeroD, i=1,nB*nC ), Dad, &
         & ( zeroD, i=1,nB*nD ), &
         & ( zeroD, i=1,nC*nD ), &
         & ( zeroD, i=1,nA*nB*nC ), &
         & ( zeroD, i=1,nA*nB*nD ), &
         & ( zeroD, i=1,nA*nC*nD ), &
         & ( zeroD, i=1,nB*nC*nD ), &
         & ( zeroD, i=1,nA*nB*nC*nD ) /), &
         & shape(cubic_RspFunc), cubic_RspFunc, &
         & perm = (/1,2,3,4/) )




    ! ****************************************************************************
    !                              MULTIPLIERS TERMS                             !
    ! ****************************************************************************
    
    call cubic_response_multiplier_contribs(molcfg,F,D,S,MyRspFunc,&
         & Da,Fa,Db,Fb,Dc,Fc,Dd,Fd,Dab,Fab,Dac,Fac,Dad,Fad,Sa,Sb,Sc,Sd,Sab,Sac,Sad,cubic_RspFunc)

    ! Done calculating cubic response function.



    !************************************************************************
    !             Store the results in the rsp_results_vec vector           !
    !************************************************************************
    m=0
    do i=1,nA
       do j=1,nB
          do k=1,nC
             do l=1,nD
                m=m+1
                rsp_results_vec(m) = cubic_RspFunc(i,j,k,l)
             enddo
          enddo
       enddo
    enddo

    ! FREE Stuff

    ! KK HACK: Quick fix of problem that "MAG, nd > 2 not implemented" in rsp_contribs.F90
    deallocate(hack_bcad, hack_bdac, hack_cbad, hack_cdab, hack_dbac, hack_dcab)

    ! Matrices always used
    Da(:)=0; Db(:)=0; Dc(:)=0; Dd(:)=0
    Dab(:,:)=0; Dac(:,:)=0; Dad(:,:)=0
    Fa(:)=0; Fb(:)=0; Fc(:)=0; Fd(:)=0
    Fab(:,:)=0; Fac(:,:)=0; Fad(:,:)=0

    ! Matrices used only for perturbation-dependent basis sets.
    if(isdef(Sa(1))) Sa(:)=0                                                                
    if(isdef(Sb(1))) Sb(:)=0
    if(isdef(Sc(1))) Sc(:)=0
    if(isdef(Sd(1))) Sd(:)=0
    if(isdef(Sab(1,1))) Sab(:,:)=0
    if(isdef(Sac(1,1))) Sac(:,:)=0
    if(isdef(Sad(1,1))) Sad(:,:)=0
  end subroutine cubic_response




  !> \brief Get perturbed density and Fock matrices for cubic response.
  !> \author Kasper Kristensen
  !> \date April 2013
  subroutine cubic_response_dens_and_fock(molcfg,F,D,S,MyRspFunc,&
       & Da,Fa,Db,Fb,Dc,Fc,Fd,Dd,Dab,Fab,Dac,Fac,Dad,Fad)

    implicit none
    !> References to molecule input, and solver and integral settings
    type(rsp_molcfg), intent(inout)   :: molcfg
    !> Unperturbed Fock matrix
    type(matrix), intent(in)          :: F
    !> Unperturbed density matrix
    type(matrix), intent(in)          :: D
    !> Unperturbed overlap matrix
    type(matrix), intent(in)          :: S
    !> Contains all information about the cubic
    !> response function or transition moment to be determined.
    type(RspFunc), intent(inout)      :: MyRspFunc
    !> Perturbed Fock and density matrices w.r.t perturbations a,b,c,d.
    type(matrix),intent(inout)  :: Da(MyRspFunc%dims(1)), &
         & Fa(MyRspFunc%dims(1)), Db(MyRspFunc%dims(2)), &
         & Fb(MyRspFunc%dims(2)), Dc(MyRspFunc%dims(3)), &
         & Fc(MyRspFunc%dims(3)), Fd(MyRspFunc%dims(4)), &
         & Dd(MyRspFunc%dims(4)), &
         & Dab(MyRspFunc%dims(1),MyRspFunc%dims(2)), &
         & Fab(MyRspFunc%dims(1),MyRspFunc%dims(2)), &
         & Dac(MyRspFunc%dims(1),MyRspFunc%dims(3)), &
         & Fac(MyRspFunc%dims(1),MyRspFunc%dims(3)), &
         & Dad(MyRspFunc%dims(1),MyRspFunc%dims(4)), &
         & Fad(MyRspFunc%dims(1),MyRspFunc%dims(4))
    logical                        :: AeqB, AeqC, AeqD, BeqC, BeqD, CeqD
    character(4)                   :: NOOP(1)
    complex(realk) :: freqA,freqB,freqC,freqD
    integer :: nA, nB, nC, nD, first_compA,first_compB, first_compC, first_compD,i,j


    ! Initialize stuff
    NOOP(1) = 'NOOP'
    nA=MyRspFunc%dims(1)
    nB=MyRspFunc%dims(2)
    nC=MyRspFunc%dims(3)
    nD=MyRspFunc%dims(4)
    freqA = MyRspFunc%freq(1)
    freqB = MyRspFunc%freq(2)
    freqC = MyRspFunc%freq(3)
    freqD = MyRspFunc%freq(4)
    first_compA = MyRspFunc%first_comp(1)
    first_compB = MyRspFunc%first_comp(2)
    first_compC = MyRspFunc%first_comp(3)
    first_compD = MyRspFunc%first_comp(4)


    ! Check whether any of the perturbations are identical
    ! If so, we avoid re-calculating various matrices and dot products.
    AeqB = Identical_operators(MyRspFunc,1,2)
    AeqC = Identical_operators(MyRspFunc,1,3)
    AeqD = Identical_operators(MyRspFunc,1,4)
    BeqC = Identical_operators(MyRspFunc,2,3)
    BeqD = Identical_operators(MyRspFunc,2,4)
    CeqD = Identical_operators(MyRspFunc,3,4)


    ! *******************************************************************************************
    !                Calculate perturbed density matrices and Fock matrices                     !
    ! *******************************************************************************************

    ! Da matrix
    call pert_dens(molcfg, S, (/ MyRspFunc%code(1) /), (/nA/), (/D/), (/F/), &
         & Da, Fa, comp=(/first_compA/), freq=(/freqA/) )

    ! Db matrix
    if(AeqB) then
       do i=1,MyRspFunc%dims(1) - MyRspFunc%first_comp(1) + 1
          Db(i)=Da(i)
          Fb(i)=Fa(i)
       end do
    else
       call pert_dens(molcfg, S, (/ MyRspFunc%code(2) /), (/nB/), (/D/), (/F/), &
            & Db, Fb, comp=(/first_compB/), freq=(/freqB/) )
    end if


    ! Dc matrix:
    ! If perturation 3 = 'EXCI' we want the DOUBLE cubic transition moment
    ! and we replace the perturbed Dc matrix by the corresponding transition density matrix.
    if(MyRspFunc%code(3) == 'EXCI') then 
       ! Dc = -Dn^T (transition density matrix)
       ! because we set omega_c = -omega_n
       ! (see the transition_moment_density_matrix subroutine and
       !  Sec. IV.H. if jcp, 129, 214108).
       ! The position of the excited state for operator C is defined by ExNumber2 in MyRspFunc
       Dc(1) = -trans(rsp_eq_sol(MyRspFunc%ExNumber2)%D)
       ! Determine transition Fock matrix Fc => Fn = G(-Dn^T)
       Fc(1)=0E0_realk*F
       call rsp_twoint(molcfg, NOOP, (/D,Dc(1)/), (/1/), Fc(1))
    else ! ordinary cubic response function or single transition moment
       if(AeqC) then
          do i=1,MyRspFunc%dims(1) - MyRspFunc%first_comp(1) + 1
             Dc(i)=Da(i)
             Fc(i)=Fa(i)
          end do
       elseif(BeqC) then
          do i=1,MyRspFunc%dims(2) - MyRspFunc%first_comp(2) + 1
             Dc(i)=Db(i)
             Fc(i)=Fb(i)
          end do
       else
          call pert_dens(molcfg, S, (/ MyRspFunc%code(3) /), (/nC/), (/D/), (/F/), &
               & Dc, Fc, comp=(/first_compC/), freq=(/freqC/) )
       end if
    end if

    ! Dd matrix:
    ! If perturation 4 = 'EXCI' we want the SINGLE or DOUBLE cubic transition moment
    ! and we replace the perturbed Dd matrix by the transition density matrix.
    if(MyRspFunc%code(4) == 'EXCI') then 
       ! Dd = Dm (transition density matrix)
       ! because we set omega_d = omega_m
       ! (see the transition_moment_density_matrix subroutine and
       !  Sec. IV.H. if jcp, 129, 214108).
       ! The position of the excited state is defined by ExNumber1 in MyRspFunc.
       Dd(1) = rsp_eq_sol(MyRspFunc%ExNumber1)%D
       ! Determine transition Fock matrix Fd => Fm = G(Dm)
       Fd(1)=0E0_realk*F
       call rsp_twoint(molcfg, NOOP, (/D,Dd(1)/), (/1/), Fd(1))
    else ! ordinary cubic response function
       if(AeqD) then
          do i=1,MyRspFunc%dims(1) - MyRspFunc%first_comp(1) + 1
             Dd(i)=Da(i)
             Fd(i)=Fa(i)
          end do
       elseif(BeqD) then
          do i=1,MyRspFunc%dims(2) - MyRspFunc%first_comp(2) + 1
             Dd(i)=Db(i)
             Fd(i)=Fb(i)
          end do
       elseif(CeqD) then
          do i=1,MyRspFunc%dims(3) - MyRspFunc%first_comp(3) + 1
             Dd(i)=Dc(i)
             Fd(i)=Fc(i)
          end do
       else
          call pert_dens(molcfg, S, (/ MyRspFunc%code(4) /), (/nD/), (/D/), (/F/), &
               & Dd, Fd, comp=(/first_compD/), freq=(/freqD/) )
       end if
    endif


    ! Dab matrix
    call pert_dens(molcfg, S, (/ MyRspFunc%code(1), MyRspFunc%code(2) /), &
         & (/ nA, nB /), (/D,Da,Db/), (/F,Fa,Fb/), &
         & Dab, Fab, &
         & comp=(/first_compA, first_compB/), &
         &  freq=(/freqA, freqB/) )

    ! Dac matrix
    if(BeqC) then
       do i=1,MyRspFunc%dims(1) - MyRspFunc%first_comp(1) + 1
          do j=1,MyRspFunc%dims(2) - MyRspFunc%first_comp(2) + 1
             Dac(i,j) = Dab(i,j)
             Fac(i,j) = Fab(i,j)
          end do
       end do
    else
       call pert_dens(molcfg, S, (/ MyRspFunc%code(1), MyRspFunc%code(3) /), &
            & (/ nA, nC /), (/D,Da,Dc/), (/F,Fa,Fc/), &
            & Dac, Fac, &
            & comp=(/first_compA, first_compC/), &
            &  freq=(/freqA, freqC/) )
    end if

    ! Dad matrix
    if(BeqD) then
       do i=1,MyRspFunc%dims(1) - MyRspFunc%first_comp(1) + 1
          do j=1,MyRspFunc%dims(2) - MyRspFunc%first_comp(2) + 1
             Dad(i,j) = Dab(i,j)
             Fad(i,j) = Fab(i,j)
          end do
       end do
    elseif(CeqD) then
       do i=1,MyRspFunc%dims(1) - MyRspFunc%first_comp(1) + 1
          do j=1,MyRspFunc%dims(3) - MyRspFunc%first_comp(3) + 1
             Dad(i,j) = Dac(i,j)
             Fad(i,j) = Fac(i,j)
          end do
       end do
    else
       call pert_dens(molcfg, S, (/ MyRspFunc%code(1), MyRspFunc%code(4) /), &
            & (/ nA, nD /), (/D,Da,Dd/), (/F,Fa,Fd/), &
            & Dad, Fad, &
            & comp=(/first_compA, first_compD/), &
            &  freq=(/freqA, freqD/) )
    end if
  end subroutine cubic_response_dens_and_fock



  !> \brief Get perturbed overlap matrices for cubic response.
  !> \author Kasper Kristensen
  !> \date April 2013
  subroutine cubic_response_overlap(molcfg,F,D,S,MyRspFunc,Sa,Sb,Sc,Sd,Sab,Sac,Sad)

    implicit none
    !> References to molecule input, and solver and integral settings
    type(rsp_molcfg), intent(inout)   :: molcfg
    !> Unperturbed Fock matrix
    type(matrix), intent(in)          :: F
    !> Unperturbed density matrix
    type(matrix), intent(in)          :: D
    !> Unperturbed overlap matrix
    type(matrix), intent(in)          :: S
    !> Contains all information about the cubic
    !> response function or transition moment to be determined.
    type(RspFunc), intent(inout)      :: MyRspFunc
    !> Overlap matrices w.r.t. a,b,c,d perturbations.
    type(matrix),intent(inout) :: Sa(MyRspFunc%dims(1)), Sb(MyRspFunc%dims(2)), &
         & Sc(MyRspFunc%dims(3)), Sd(MyRspFunc%dims(4)), &
         & Sab(MyRspFunc%dims(1),MyRspFunc%dims(2)), &
         & Sac(MyRspFunc%dims(1),MyRspFunc%dims(3)), &
         & Sad(MyRspFunc%dims(1),MyRspFunc%dims(4))
   logical                        :: AeqB, AeqC, AeqD, BeqC, BeqD, CeqD,pdbs(4)
    integer :: nA, nB, nC, nD, first_compA,first_compB, first_compC, first_compD,i,j


    ! Initialize stuff
    nA=MyRspFunc%dims(1)
    nB=MyRspFunc%dims(2)
    nC=MyRspFunc%dims(3)
    nD=MyRspFunc%dims(4)
    first_compA = MyRspFunc%first_comp(1)
    first_compB = MyRspFunc%first_comp(2)
    first_compC = MyRspFunc%first_comp(3)
    first_compD = MyRspFunc%first_comp(4)

    ! Determine if the operators in fields contain perturbation dependent basis sets.
    ! E.g. if MyRspFunc%code(1:4)= ('GEO','EL', 'EL', 'MAG') 
    ! then pdbs=(.true. , .false., .false., .true.).
    pdbs=pert_basdep(MyRspFunc%code(1:4))


    ! Check whether any of the perturbations are identical
    ! If so, we avoid re-calculating various matrices and dot products.
    AeqB = Identical_operators(MyRspFunc,1,2)
    AeqC = Identical_operators(MyRspFunc,1,3)
    AeqD = Identical_operators(MyRspFunc,1,4)
    BeqC = Identical_operators(MyRspFunc,2,3)
    BeqD = Identical_operators(MyRspFunc,2,4)
    CeqD = Identical_operators(MyRspFunc,3,4)


    ! *******************************************************************************************
    !               Calculate overlap matrices for perturbation-dependent basis sets            !
    ! *******************************************************************************************

    ! Determine Sa
    if(pdbs(1)) then
       call rsp_oneint(molcfg, S, (/MyRspFunc%code(1) /), (/nA /), & 
            comp=(/first_compA/), S=Sa)
    else
       do i=1,nA
          Sa(i)=0E0_realk*S
       enddo
    endif

    ! Determine Sb
    if(pdbs(2)) then
       call rsp_oneint(molcfg, S, (/MyRspFunc%code(2) /), (/nB /), &
            comp=(/first_compB/), S=Sb)
    else
       do i=1,nB
          Sb(i)=0E0_realk*S
       enddo
    endif

    ! Determine Sc
    if(pdbs(3) .and. MyRspFunc%code(3) /= 'EXCI') then
       call rsp_oneint(molcfg, S, (/MyRspFunc%code(3) /), (/nC /), &
            comp=(/first_compC/), S=Sc)
    else ! If not pdbs or if we take the residue: freqC -> -Excitation energy.
       do i=1,nC
          Sc(i)=0E0_realk*S
       enddo
    endif

    ! Determine Sd
    if(pdbs(4) .and. MyRspFunc%code(4) /= 'EXCI') then
       call rsp_oneint(molcfg, S, (/MyRspFunc%code(4) /), (/nD /), &
            comp=(/first_compD/), S=Sd)
    else ! If not pdbs or if we take the residue: freqD -> Excitation energy.
       do i=1,nD
          Sd(i)=0E0_realk*S
       enddo
    endif


    ! Determine Sab
    if(pdbs(1) .and. pdbs(2)) then
       call rsp_oneint(molcfg, S, (/ MyRspFunc%code(1),MyRspFunc%code(2) /), &
            (/ nA,nB /), &
            comp=(/ first_compA,first_compB /), S=Sab)
    else
       do i=1,nA
          do j=1,nB
             Sab(i,j)=0E0_realk*S
          end do
       enddo
    endif

    ! Determine Sac
    IfSac: if(pdbs(1) .and. pdbs(3) .and. (MyRspFunc%code(3) /= 'EXCI') ) then

       IfBeqC: if(BeqC) then
          do i=1,nA
             do j=1,nC
                Sac(i,j) = Sab(i,j)
             enddo
          enddo
       else
          call rsp_oneint(molcfg, S, (/ MyRspFunc%code(1),MyRspFunc%code(3) /), &
               (/ nA,nC /), &
               comp=(/ first_compA,first_compC /), S=Sac)
       end if IfBeqC

    else
       ! Set Sac=0 if a and c do not both use PDBS
       ! and always if c=EXCI.
       do i=1,nA
          do j=1,nC
             Sac(i,j)=0E0_realk*S
          end do
       enddo
    end if IfSac


    ! Determine Sad
    IfSad: if(pdbs(1) .and. pdbs(4) .and. (MyRspFunc%code(4) /= 'EXCI') ) then

       IfBeqD_or_CeqD: if(BeqD) then
          do i=1,nA
             do j=1,nD
                Sad(i,j) = Sab(i,j)
             enddo
          enddo
       elseif(CeqD) then
          do i=1,nA
             do j=1,nD
                Sad(i,j) = Sac(i,j)
             enddo
          enddo
       else
          call rsp_oneint(molcfg, S, (/ MyRspFunc%code(1),MyRspFunc%code(4) /), &
               (/ nA,nD /), &
               comp=(/ first_compA,first_compD /), S=Sad)
       end if IfBeqD_or_CeqD
    else 
       ! Set Sad=0 if a and d do not both use PDBS
       ! and always if d=EXCI.
       do i=1,nA
          do j=1,nD
             Sad(i,j)=0E0_realk*S
          end do
       enddo
    end if IfSad

  end subroutine cubic_response_overlap





  !> \brief Calculate multiplier contrbutions to cubic response function.
  !> \author Kasper Kristensen
  !> \date April 2013
  subroutine cubic_response_multiplier_contribs(molcfg,F,D,S,MyRspFunc,&
       & Da,Fa,Db,Fb,Dc,Fc,Dd,Fd,Dab,Fab,Dac,Fac,Dad,Fad,Sa,Sb,Sc,Sd,Sab,Sac,Sad,cubic_RspFunc)

    implicit none
    !> References to molecule input, and solver and integral settings
    type(rsp_molcfg), intent(inout)   :: molcfg
    !> Unperturbed Fock matrix
    type(matrix), intent(in)          :: F
    !> Unperturbed density matrix
    type(matrix), intent(in)          :: D
    !> Unperturbed overlap matrix
    type(matrix), intent(in)          :: S
    !> Contains all information about the cubic
    !> response function or transition moment to be determined.
    type(RspFunc), intent(inout)      :: MyRspFunc
    !> Cubic response function where multiplier terms are added (note: cubic_rspfunc is not zeroed here!)
    complex(realk),intent(inout) :: cubic_RspFunc(MyRspFunc%dims(1),MyRspFunc%dims(2),MyRspFunc%dims(3), MyRspFunc%dims(4))

    !> Pertubed density,Fock, and overlap matrices w.r.t. perturbations a,b,c,d.
    type(matrix),intent(in)  :: Da(MyRspFunc%dims(1)), &
         & Fa(MyRspFunc%dims(1)), &
         & Db(MyRspFunc%dims(2)), &
         & Fb(MyRspFunc%dims(2)), &
         & Dc(MyRspFunc%dims(3)), &
         & Fc(MyRspFunc%dims(3)), &
         & Fd(MyRspFunc%dims(4)), &
         & Dd(MyRspFunc%dims(4)), &
         & Dab(MyRspFunc%dims(1),MyRspFunc%dims(2)), &
         & Fab(MyRspFunc%dims(1),MyRspFunc%dims(2)), &
         & Dac(MyRspFunc%dims(1),MyRspFunc%dims(3)), &
         & Fac(MyRspFunc%dims(1),MyRspFunc%dims(3)), &
         & Dad(MyRspFunc%dims(1),MyRspFunc%dims(4)), &
         & Fad(MyRspFunc%dims(1),MyRspFunc%dims(4)), &
         & Sa(MyRspFunc%dims(1)), Sb(MyRspFunc%dims(2)), &
         & Sc(MyRspFunc%dims(3)), Sd(MyRspFunc%dims(4)), &
         & Sab(MyRspFunc%dims(1),MyRspFunc%dims(2)), &
         & Sac(MyRspFunc%dims(1),MyRspFunc%dims(3)), &
         & Sad(MyRspFunc%dims(1),MyRspFunc%dims(4))

    type(matrix) :: Ybc(MyRspFunc%dims(2), MyRspFunc%dims(3)), &
         & Ybd(MyRspFunc%dims(2), MyRspFunc%dims(4)), &
         & Ycd(MyRspFunc%dims(3), MyRspFunc%dims(4)), &
         & Zbc(MyRspFunc%dims(2), MyRspFunc%dims(3)), &
         & Zbd(MyRspFunc%dims(2), MyRspFunc%dims(4)), &
         & Zcd(MyRspFunc%dims(3), MyRspFunc%dims(4)), &
         & Wbcd,Ybcd,Zbcd, &
         & lambdaA(MyRspFunc%dims(1)), &
         & lambdaAB(MyRspFunc%dims(1),MyRspFunc%dims(2)), &
         & lambdaAC(MyRspFunc%dims(1),MyRspFunc%dims(3)), &
         & lambdaAD(MyRspFunc%dims(1),MyRspFunc%dims(4)), &
         & zetaA(MyRspFunc%dims(1)), &
         & zetaAB(MyRspFunc%dims(1),MyRspFunc%dims(2)), &
         & zetaAC(MyRspFunc%dims(1),MyRspFunc%dims(3)), &
         & zetaAD(MyRspFunc%dims(1),MyRspFunc%dims(4)), zeroD, zeroF

    logical                        :: pdbs(4), BeqC, BeqD, CeqD
    integer                        :: i,j, k, l, m
    complex(realk) :: freqA,freqB,freqC,freqD
    integer :: nA, nB, nC, nD, first_compA,first_compB, first_compC, first_compD


    ! Initialize stuff
    nA=MyRspFunc%dims(1)
    nB=MyRspFunc%dims(2)
    nC=MyRspFunc%dims(3)
    nD=MyRspFunc%dims(4)
    freqA = MyRspFunc%freq(1)
    freqB = MyRspFunc%freq(2)
    freqC = MyRspFunc%freq(3)
    freqD = MyRspFunc%freq(4)
    first_compA = MyRspFunc%first_comp(1)
    first_compB = MyRspFunc%first_comp(2)
    first_compC = MyRspFunc%first_comp(3)
    first_compD = MyRspFunc%first_comp(4)
    BeqC = Identical_operators(MyRspFunc,2,3)
    BeqD = Identical_operators(MyRspFunc,2,4)
    CeqD = Identical_operators(MyRspFunc,3,4)
    zeroD = 0*D
    zeroF = 0*F


    ! Below the multipliers are first constructed, and subsequently contracted with
    ! Y,Z, and W matrices.


    ! ****************************************************************************
    !                    Contruct lambda and zeta multipliers                    !
    !                    (Eqs. 220 and 224 in jcp 129, 214108)                   !
    ! ****************************************************************************
    MultiplierLoopA: do i=1,nA

       ! '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
       !                      Zeroth order multipliers                 !
       ! '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
       ! Only if at least one of the bcd perturbations use PDBS
       if(any(pdbs(2:4))) then
          lambdaA(i) = Da(i)*S*D - D*S*Da(i)
          zetaA(i) = Fa(i)*D*S + S*D*Fa(i) - Fa(i) - F*D*Sa(i) - Sa(i)*D*F
       end if


       ! '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
       !                       First order multipliers                   !
       ! '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

       ! lambdaAB and zetaAB
       MultiplierLoopB: do j=1,nB

          lambdaAB(i,j) = Dab(i,j)*S*D - D*S*Dab(i,j) + &
               & Da(i)*Sb(j)*D-D*Sb(j)*Da(i) + Da(i)*S*Db(j) - Db(j)*S*Da(i)

          zetaAB(i,j) = Fab(i,j)*D*S + S*D*Fab(i,j) + Fa(i)*Db(j)*S + S*Db(j)*Fa(i) &
               - Fab(i,j) &
               + Fa(i)*(D*Sb(j)) + (Sb(j)*D)*Fa(i)  &
               -Fb(j)*(D*Sa(i)) - (Sa(i)*D)*Fb(j) - F*(Db(j)*Sa(i)) - (Sa(i)*Db(j))*F &
               - F*(D*Sab(i,j)) - (Sab(i,j)*D)*F &
               + 0.5E0_realk*freqB*(Sb(j)*D*Sa(i) - Sa(i)*D*Sb(j)) &
               + freqB*(S*Db(j)*Sa(i) - Sa(i)*Db(j)*S)


       end do MultiplierLoopB


       ! lambdaAC and zetaAC. Copy if C is identical to B
       MultiplierLoopC: do k=1,nC

          if(BeqC) then
             lambdaAC(i,k) = lambdaAB(i,k)
             zetaAC(i,k) = zetaAB(i,k)
          else
             lambdaAC(i,k) = Dac(i,k)*S*D - D*S*Dac(i,k) + &
                  & Da(i)*Sc(k)*D-D*Sc(k)*Da(i) + Da(i)*S*Dc(k) - Dc(k)*S*Da(i)

             zetaAC(i,k) = Fac(i,k)*D*S + S*D*Fac(i,k) + Fa(i)*Dc(k)*S + S*Dc(k)*Fa(i) &
                  - Fac(i,k) &
                  + Fa(i)*(D*Sc(k)) + (Sc(k)*D)*Fa(i)  &
                  -Fc(k)*(D*Sa(i)) - (Sa(i)*D)*Fc(k) - F*(Dc(k)*Sa(i)) - (Sa(i)*Dc(k))*F &
                  - F*(D*Sac(i,k)) - (Sac(i,k)*D)*F &
                  + 0.5E0_realk*freqC*(Sc(k)*D*Sa(i) - Sa(i)*D*Sc(k)) &
                  + freqC*(S*Dc(k)*Sa(i) - Sa(i)*Dc(k)*S)
          end if

       end do MultiplierLoopC


       ! lambdaAD and zetaAD. Copy if D is identical to B or C
       MultiplierLoopD: do l=1,nD
          if(BeqD) then
             lambdaAD(i,l) = lambdaAB(i,l)
             zetaAD(i,l) = zetaAB(i,l)
          elseif(CeqD) then
             lambdaAD(i,l) = lambdaAC(i,l)
             zetaAD(i,l) = zetaAC(i,l)
          else
             lambdaAD(i,l) = Dad(i,l)*S*D - D*S*Dad(i,l) + &
                  & Da(i)*Sd(l)*D-D*Sd(l)*Da(i) + Da(i)*S*Dd(l) - Dd(l)*S*Da(i)
             zetaAD(i,l) = Fad(i,l)*D*S + S*D*Fad(i,l) + Fa(i)*Dd(l)*S + S*Dd(l)*Fa(i) &
                  - Fad(i,l) &
                  + Fa(i)*(D*Sd(l)) + (Sd(l)*D)*Fa(i)  &
                  -Fd(l)*(D*Sa(i)) - (Sa(i)*D)*Fd(l) - F*(Dd(l)*Sa(i)) - (Sa(i)*Dd(l))*F &
                  - F*(D*Sad(i,l)) - (Sad(i,l)*D)*F &
                  + 0.5E0_realk*freqD*(Sd(l)*D*Sa(i) - Sa(i)*D*Sd(l)) &
                  + freqD*(S*Dd(l)*Sa(i) - Sa(i)*Dd(l)*S)
          end if
       end do MultiplierLoopD


    end do MultiplierLoopA

    !************************************************************************
    !                Construct second order Y,Z and W matrices              !
    !************************************************************************

    ! Ybc and Zbc
    call pert_scf_eq(molcfg, S, (/ MyRspFunc%code(2), MyRspFunc%code(3) /), &
         & (/ nB,nC /), &
         & (/D,Db,Dc, (zeroD, i=1,nB*nC ) /), &
         & (/F,Fb,Fc, (zeroF, i=1,nB*nC ) /), &
         & Ybc, Zbc, &
         & comp = (/first_compB, first_compC/), & 
         & freq=(/freqB, freqC/) )

    ! Ybd and Zbd: Copy Ybc and Zbc if c=d
    if(CeqD) then
       do i=1,nB
          do j=1,nD
             Ybd(i,j) = Ybc(i,j)
             Zbd(i,j) = Zbc(i,j)
          enddo
       enddo
    else
       call pert_scf_eq(molcfg, S, (/ MyRspFunc%code(2), MyRspFunc%code(4) /), &
            & (/ nB,nD /), &
            & (/D,Db,Dd, (zeroD, i=1,nB*nD ) /), &
            & (/F,Fb,Fd, (zeroF, i=1,nB*nD ) /), &
            & Ybd, Zbd, &
            & comp = (/first_compB, first_compD/), & 
            & freq=(/freqB, freqD/) )
    endif


    ! Ycd and Zcd: Copy from Ybc/Zbc or Ybd/Zbd if Ycd/Zcd is identical.
    if(BeqC) then
       do i=1,nC
          do j=1,nD
             Ycd(i,j) = Ybd(i,j)
             Zcd(i,j) = Zbd(i,j)
          enddo
       enddo
    elseif(BeqD) then
       do i=1,nC
          do j=1,nD
             Ycd(i,j) = Ybc(j,i)
             Zcd(i,j) = Zbc(j,i)
          enddo
       enddo
    else
       call pert_scf_eq(molcfg, S, (/ MyRspFunc%code(3), MyRspFunc%code(4) /), &
            & (/ nC,nD /), &
            & (/D,Dc,Dd, (zeroD, i=1,nC*nD )/), &
            & (/F,Fc,Fd, (zeroF, i=1,nC*nD )/), &
            & Ycd, Zcd, &
            & comp = (/first_compC, first_compD/), & 
            & freq=(/freqC, freqD/) )
    endif


    !************************************************************************
    !                     Construct Ybcd, Zbcd, and Wbcd                    !
    !************************************************************************
    j_loop: do j=1,nB
       k_loop: do k=1,nC
          l_loop: do l=1,nD

             ! Construct Ybcd' and Zbcd' (keeping only terms containing first order derivatives)
             ! These terms contribute only if at least one 
             ! of the bcd perturbations use PDBS.
             IfYbcdZbcd: if(any(pdbs(2:4))) then
                Ybcd = Fb(j)*(Dc(k)*Sd(l)) - (Sd(l)*Dc(k))*Fb(j) + &
                     Fb(j)*(Dd(l)*Sc(k)) - (Sc(k)*Dd(l))*Fb(j) &
                     + Fc(k)*(Db(j)*Sd(l)) - (Sd(l)*Db(j))*Fc(k) + &
                     Fc(k)*(Dd(l)*Sb(j)) - (Sb(j)*Dd(l))*Fc(k) &
                     + Fd(l)*(Db(j)*Sc(k)) - (Sc(k)*Db(j))*Fd(l) + &
                     Fd(l)*(Dc(k)*Sb(j)) - (Sb(j)*Dc(k))*Fd(l) &
                     - (freqC + 0.5E0_realk*(freqB+freqD)) * &
                     ( (Sb(j)*Dc(k))*Sd(l) + (Sd(l)*Dc(k))*Sb(j) ) &
                     - (freqB + 0.5E0_realk*(freqC+freqD)) * &
                     ( (Sc(k)*Db(j))*Sd(l) + (Sd(l)*Db(j))*Sc(k) ) &
                     - (freqD + 0.5E0_realk*(freqB+freqC)) * &
                     ( (Sb(j)*Dd(l))*Sc(k) + (Sc(k)*Dd(l))*Sb(j) ) 

                Zbcd = Db(j)*( Sc(k)*Dd(l) + Sd(l)*Dc(k) ) &
                     + Dc(k)*( Sb(j)*Dd(l) + Sd(l)*Db(j) ) &
                     + Dd(l)*( Sb(j)*Dc(k) + Sc(k)*Db(j) )

             end if IfYbcdZbcd

             ! Construct Wbcd' [see Eq. 95 in jcp, 129, 214108 (2008)]
             ! The prime indicates that only first order matrices should be included.
             ! (Since Wbcd' is only needed this one time
             !  we have not included third order W matrices in the get_W_matrix routine.
             !  Instead we construct Wbcd' explicitly here.)
             IfWbcd: if(pdbs(1)) then
                Wbcd = Db(j)*(Fc(k)*Dd(l)+Fd(l)*Dc(k)) + &
                     & Dc(k)*(Fb(j)*Dd(l)+Fd(l)*Db(j)) + &
                     & Dd(l)*(Fb(j)*Dc(k)+Fc(k)*Db(j)) + &
                     & 0.5E0_realk*( &
                     & (freqB-freqD)*(Db(j)*Sc(k)*Dd(l)-Dd(l)*Sc(k)*Db(j)) &
                     & + (freqC-freqD)*(Dc(k)*Sb(j)*Dd(l)-Dd(l)*Sb(j)*Dc(k)) &
                     & + (freqB-freqC)*(Db(j)*Sd(l)*Dc(k)-Dc(k)*Sd(l)*Db(j)) &
                     & )
             end if IfWbcd


             !************************************************************************
             !     Calculate multiplier contributions to cubic response function     !
             !************************************************************************


             i_loop: do i=1,nA

                if(any(pdbs(2:4))) then
                   ! Zeroth order lambda contribution, lambdaA Ybcd
                   cubic_RspFunc(i,j,k,l) = cubic_RspFunc(i,j,k,l) &
                        -trace(lambdaA(i),Ybcd)

                   ! Zeroth order zeta contribution, zetaA Zbcd
                   cubic_RspFunc(i,j,k,l) = cubic_RspFunc(i,j,k,l) &
                        -trace(zetaA(i),Zbcd)
                end if

                ! -Sa Wbcd' contribution
                if(pdbs(1)) then
                   cubic_RspFunc(i,j,k,l) = cubic_RspFunc(i,j,k,l) &
                        -trace(Sa(i),Wbcd)
                endif


                ! For Hartree-Fock, perturbation-independent basis, and
                ! operators that are linear in the perturbing field(s), the
                ! first-order multipliers terms below are the only terms that contribute
                ! to the cubic response function.

                ! first order lambda contributions
                cubic_RspFunc(i,j,k,l) = cubic_RspFunc(i,j,k,l) &
                     -trace(lambdaAB(i,j),Ycd(k,l)) &
                     -trace(lambdaAC(i,k),Ybd(j,l)) &
                     -trace(lambdaAD(i,l),Ybc(j,k))

                ! First order zeta contributions
                cubic_RspFunc(i,j,k,l) = cubic_RspFunc(i,j,k,l) &
                     -trace(zetaAB(i,j),Zcd(k,l)) &
                     -trace(zetaAC(i,k),Zbd(j,l)) &
                     -trace(zetaAD(i,l),Zbc(j,k))

             end do i_loop
          end do l_loop
       end do k_loop
    end do j_loop


    ! Free stuff
    lambdaAB(:,:)=0; lambdaAC(:,:)=0; lambdaAD(:,:)=0
    zetaAB(:,:)=0; zetaAC(:,:)=0; zetaAD(:,:)=0
    Ybc(:,:)=0; Ybd(:,:)=0; Ycd(:,:)=0
    Zbc(:,:)=0; Zbd(:,:)=0; Zcd(:,:)=0

    ! Matrices used only for perturbation-dependent basis sets.
    if(isdef(lambdaA(1))) lambdaA(:)=0 
    if(isdef(zetaA(1))) zetaA(:)=0 
    if(isdef(Wbcd)) Wbcd=0
    if(isdef(Ybcd)) Ybcd=0
    if(isdef(Zbcd)) Zbcd=0

  end subroutine cubic_response_multiplier_contribs




  !> \brief Calculates W matrix required for calculating response functions
  !> using perturbations-dependent basis sets, see Eq. (95) in jcp 129, 214108.
  !> Currently up to second order W matrices are implemented.
  !> For second order matrices the order of operators in W is W(x,y)
  !> \author Kasper Kristensen
  !> \date 2010-01 

  subroutine get_W_matrix(F,D,S,MyRspFunc,order,pdbs,W,ops,Dx,Fx,Sx,&
       Dy,Fy,Sy,Dxy,Fxy)
    implicit none
    !> Unperturbed Fock matrix
    type(matrix), intent(in)          :: F
    !> Unperturbed density matrix
    type(matrix), intent(in)          :: D
    !> Unperturbed overlap matrix
    type(matrix), intent(in)          :: S
    !> Response function information
    type(RspFunc), intent(in)      :: MyRspFunc
    !> Order of W matrix
    integer, intent(in) :: order
    !> Use perturbation-dependent basis sets for operators A,B,C, ...
    logical, intent(in) :: pdbs(MyRspFunc%order)
    !> W matrix to be generated
    type(matrix), intent(inout) :: W(*)
    !> Operator positions for operators
    !> for which W is requested [e.g. for Wbd ops = (/2,4/)]
    integer, intent(in), optional :: ops(:)
    !> First order (x) density, Fock, and overlap matrix
    type(matrix), intent(in), optional :: Dx(*),Fx(*),Sx(*)
    !> First order (y) density, Fock, and overlap matrix
    type(matrix), intent(in), optional :: Dy(*),Fy(*),Sy(*)
    !> Second order (xy) density, Fock, and overlap matrix
    type(matrix), intent(in), optional :: Dxy(*),&
         Fxy(*)
    !>
    ! Local parameters
    integer :: i,j,k
    logical :: pdbs_tmp(MyRspFunc%order)
    complex :: freqX, freqY, freqXY

    ! Initialize stuff
    pdbs_tmp = pdbs
    if(present(ops)) freqX = MyRspFunc%freq(ops(1))

    select case(order)

    case(0) ! W

       ! W = DFD
       W(1) = D*F*D

    case(1) ! Wx

       ! Calculate Wx, where x is the ops(1)'th conponent in MyRspFunc.
       ! Wx is only nonzero if all the other components than x in MyRspFunc use PDBS.
       ! E.g. for the quadratic response function and x=c the contribution is
       ! -Sab*Wc
       ! and we therefore require both a and b (but not necessarily c) to use PDBS.

       pdbs_tmp(ops(1)) = .true.

       do i=1,MyRspFunc%dims(ops(1))
          if(all(pdbs_tmp)) then
             W(i) = Dx(i)*(F+(0.5E0_realk*freqX)*S)*D + &
                  & D*(F-(0.5E0_realk*freqX)*S)*Dx(i) + D*Fx(i)*D
          else
             W(i)=0E0_realk*D
          endif
       enddo

    case(2) ! Wxy

       ! Calculate Wxy, where x and y are the op1'th and op2'th conponents in MyRspFunc.
       ! Wxy is only nonzero if all other components than x and y in MyRspFunc use PDBS.

       pdbs_tmp(ops(1)) =.true.
       pdbs_tmp(ops(2)) =.true.
       freqY= MyRspFunc%freq(ops(2))
       freqXY = freqX+freqY

       k=0
       do j=1,MyRspFunc%dims(ops(2))
          do i=1,MyRspFunc%dims(ops(1))
             k=k+1
             if(all(pdbs_tmp)) then

                W(k) = Dxy(k)*F*D + D*F*Dxy(k) + &
                     & D*Fxy(k)*D +&
                     & 0.5E0_realk*(freqX+freqY)*(Dxy(k)*S*D - D*S*Dxy(k)) +&
                     & Dx(i)*F*Dy(j) + Dy(j)*F*Dx(i) +&
                     & Dx(i)*Fy(j)*D + D*Fy(j)*Dx(i) +&
                     & Dy(j)*Fx(i)*D + D*Fx(i)*Dy(j) +&
                     & 0.5E0_realk*freqX*(Dx(i)*Sy(j)*D - D*Sy(j)*Dx(i)) +&
                     & 0.5E0_realk*freqY*(Dy(j)*Sx(i)*D - D*Sx(i)*Dy(j)) +&
                     & 0.5E0_realk*(freqX-freqY)*(Dx(i)*S*Dy(j) - Dy(j)*S*Dx(i))
             else
                W(k)=0E0_realk*D
             endif
          enddo
       end do

    case default
       CALL lsQUIT('get_W_matrix: Only up to second order W matrices have been implemented.',-1)

    end select
  end subroutine get_W_matrix




  !> \brief Prints out information about all response functions and/or
  !> transition moments requested in input.
  !> If print_results=.true. the final response response results are also printed.
  !> \author Kasper Kristensen
  !> \date 2010-01 
  subroutine print_response_func(RspFunc_stack,n_RspFunc,print_results,lupri)
    !> Response function stack containing information and (possibly) results
    !> for all response requested in input.
    type(RspFunc), intent(in)   :: RspFunc_stack(n_RspFunc)
    !> Number of response functions and/or transition moments.
    integer, intent(in)         :: n_RspFunc
    !> Should the response results be printed (.true.) or not (.false.).
    logical, intent(in)         :: print_results
    !> Unit to print to
    integer, intent(in)         :: lupri
    !> Local parameters
    integer                     :: i,j, k,l,m,jx,kx,lx,mx,num_aux, counter
    character(8) :: prop_auxlab(0:9)


    RspFunc_loop: do i=1,n_RspFunc
       write(lupri,*)
       write(lupri,*)
       write(lupri,*)
       write(lupri, '(A65)') ' ****************************************************************'

       if(print_results) then
          write(lupri, '(A2,21X,A18,I2,21X,A1)') " *"," RESPONSE RESULTS ", i, '*'
       else
          write(lupri, '(A2,22X,A16,I2,22X,A1)') " *"," RESPONSE INPUT ", i, '*'
       endif

       if(RspFunc_stack(i)%trans_moment .and. &
            & (RspFunc_stack(i)%trans_moment_order == 1) ) then
          write(lupri, '(A2,11X,A39,12X,A1)') ' *','Property type: Single transition moment', '*'
       elseif(RspFunc_stack(i)%trans_moment .and. &
            & RspFunc_stack(i)%trans_moment_order == 2) then
          write(lupri, '(A2,11X,A39,12X,A1)') ' *','Property type: Double transition moment', '*'
       else
          write(lupri, '(A2,10X,A40,I1,A1,10X,A1)') ' *','Property type: Response function (order ', &
               &RspFunc_stack(i)%order, ')','*'
       endif

       write(lupri, '(A65)') ' ****************************************************************'

       Operator_loop: do j=1,RspFunc_stack(i)%order 
          write(lupri,*)
          write(lupri, '(A10, I1)') ' OPERATOR ', j
          write(lupri, '(A11)') " ''''''''''"
          ! If AUX property we print out which integral, e.g. XDIPLEN, XANGMOM, ...
          ! (see subroutine read_words_from_string in response_driver_module for more details)
          if(RspFunc_stack(i)%code(j)(1:3) == 'AUX') then
             read (RspFunc_stack(i)%code(j)(4:4), '(I10)') num_aux
             call lsquit('prop_auxlab not defined',-1)
             !                     write(lupri, '(A16,A8,A2,A4,A1)') ' Operator type: ', prop_auxlab(num_aux), &
             !                          & ' (', RspFunc_stack(i)%code(j), ')'
          else
             write(lupri, '(A16,A4)') ' Operator type: ', RspFunc_stack(i)%code(j)
          endif


          ExciInput: if(RspFunc_stack(i)%code(j)(1:4) == 'EXCI') then ! Operator is 'EXCI'

             ! (Convention: The first residue is taken for the last operator and
             !  the second residue is taken for the second-last operator)

             ! Second excitation energy for double residue (second operator, j=order-1)
             if(RspFunc_stack(i)%trans_moment_order==2 .and. j==RspFunc_stack(i)%order-1) then 
                write(lupri, '(A38,F12.7)') ' Real frequency (-excitation energy): ', &
                     real(RspFunc_stack(i)%freq(j))
                write(lupri, '(A40,I3)') ' Excited state number for this operator:', &
                     & RspFunc_stack(i)%ExNumber2
                ! Single residue og first excitation energy for double residue
             else 
                write(lupri, '(A37,F12.7)') ' Real frequency (excitation energy): ', &
                     real(RspFunc_stack(i)%freq(j))
                write(lupri, '(A40,I3)') ' Excited state number for this operator:', &
                     & RspFunc_stack(i)%ExNumber1
             endif


          else ! Normal operator

             write(lupri, '(A17,I5)') ' First component:', RspFunc_stack(i)%first_comp(j)
             write(lupri, '(A17,I5)') ' Last component: ', &
                  RspFunc_stack(i)%dims(j)+RspFunc_stack(i)%first_comp(j)-1


             write(lupri, '(A17,F12.7)') ' Real frequency: ', real(RspFunc_stack(i)%freq(j))
             !                     if(cfg_rsp_imfreqs_specified .or. cfg_rsp_complex) then
             !                        write(lupri, '(A22,F12.7)') ' Imaginary frequency: ', imag(RspFunc_stack(i)%freq(j))
             !                     endif

          end if ExciInput


       end do Operator_loop


       LinearResponse: if(print_results .and. RspFunc_stack(i)%order == 2) then 


          LinTransMoment: if(RspFunc_stack(i)%trans_moment) then ! Transition moment

             write(lupri,*)
             write(lupri, '(A25)') ' LINEAR TRANSITION MOMENT'
             write(lupri, '(A25)') " ''''''''''''''''''''''''"


             write(lupri,'(1X,A)') '    Comp1            Result '

             do j=1,RspFunc_stack(i)%dims(1)
                jx = j + RspFunc_stack(i)%first_comp(1)-1
                write(lupri,'(1X,I7,F21.7)') jx,real(RspFunc_stack(i)%result(j))
             enddo

          else ! Linear response function

             write(lupri,*)
             write(lupri, '(A36)') ' LINEAR RESPONSE RESULTS (REAL PART)'
             write(lupri, '(A36)') " '''''''''''''''''''''''''''''''''''"


             write(lupri,'(1X,A)') '    Comp1   Comp2           Result '

             counter = 0
             do j=1,RspFunc_stack(i)%dims(1)
                do k=1,RspFunc_stack(i)%dims(2)
                   jx = j + RspFunc_stack(i)%first_comp(1)-1
                   kx = k + RspFunc_stack(i)%first_comp(2)-1
                   counter=counter+1
                   write(lupri,'(1X,2I7,F21.7)') jx,kx,real(RspFunc_stack(i)%result(counter))
                enddo
             enddo

             !                     LinearImag: if(cfg_rsp_imfreqs_specified .or. cfg_rsp_complex) then
             !                        write(lupri,*)
             !                        write(lupri, '(A36)') ' LINEAR RESPONSE RESULTS (IMAG PART)'
             !                        write(lupri, '(A36)') " '''''''''''''''''''''''''''''''''''"
             !
             !                        write(lupri,*) '  Comp1   Comp2         Result '
             !
             !                        counter = 0
             !                        do j=1,RspFunc_stack(i)%dims(1)
             !                           do k=1,RspFunc_stack(i)%dims(2)
             !                              jx = j + RspFunc_stack(i)%first_comp(1)-1
             !                              kx = k + RspFunc_stack(i)%first_comp(2)-1
             !                              counter=counter+1
             !                              write(lupri,'(X,2I6,F21.7)') jx,kx,aimag(RspFunc_stack(i)%result(counter))
             !                           enddo
             !                        enddo
             !
             !                     end if LinearImag

          end if LinTransMoment

          write(lupri,*) 

       end if LinearResponse





       QuadraticResponse: if(print_results .and. RspFunc_stack(i)%order == 3) then 


          QuadType: if(RspFunc_stack(i)%trans_moment) then ! Transition moment

             TransMomentOrder: if(RspFunc_stack(i)%trans_moment_order == 1) then ! single transition moment

                write(lupri,*)
                write(lupri, '(A35)') ' QUADRATIC SINGLE TRANSITION MOMENT'
                write(lupri, '(A35)') " ''''''''''''''''''''''''''''''''''"

                write(lupri,'(1X,A)') '    Comp1   Comp2           Result '
                counter = 0
                do j=1,RspFunc_stack(i)%dims(1)
                   do k=1,RspFunc_stack(i)%dims(2)
                      counter=counter+1
                      jx = j + RspFunc_stack(i)%first_comp(1)-1
                      kx = k + RspFunc_stack(i)%first_comp(2)-1
                      write(lupri,'(1X,2I7,F21.7)') jx,kx,real(RspFunc_stack(i)%result(counter))
                   enddo
                enddo

             else ! Double transition moment

                write(lupri,*)
                write(lupri, '(A35)') ' QUADRATIC DOUBLE TRANSITION MOMENT'
                write(lupri, '(A35)') " ''''''''''''''''''''''''''''''''''"

                write(lupri,'(1X,A)') '    Comp1            Result '
                do j=1,RspFunc_stack(i)%dims(1)
                   jx = j + RspFunc_stack(i)%first_comp(1)-1
                   write(lupri,'(1X,I7,F21.7)') jx,real(RspFunc_stack(i)%result(j))
                enddo
             end if TransMomentOrder


          else ! Quadratic response function

             write(lupri,*)
             write(lupri, '(A39)') ' QUADRATIC RESPONSE RESULTS (REAL PART)'
             write(lupri, '(A39)') " ''''''''''''''''''''''''''''''''''''''"

             write(lupri,'(1X,A)') '    Comp1  Comp2  Comp3           Result '

             counter = 0
             do j=1,RspFunc_stack(i)%dims(1)
                do k=1,RspFunc_stack(i)%dims(2)
                   do l=1,RspFunc_stack(i)%dims(3)
                      counter=counter+1
                      jx = j + RspFunc_stack(i)%first_comp(1)-1
                      kx = k + RspFunc_stack(i)%first_comp(2)-1
                      lx = l + RspFunc_stack(i)%first_comp(3)-1
                      write(lupri,'(1X,3I7,F21.7)') jx,kx,lx,real(RspFunc_stack(i)%result(counter))
                   enddo
                enddo
             enddo


             !                     QuadraticImag: if(cfg_rsp_imfreqs_specified .or. cfg_rsp_complex) then
             !
             !                        write(lupri,*)
             !                        write(lupri, '(A39)') ' QUADRATIC RESPONSE RESULTS (IMAG PART)'
             !                        write(lupri, '(A39)') " ''''''''''''''''''''''''''''''''''''''"
             !
             !                        write(lupri,'(1X,A)') '    Comp1   Comp2   Comp3           Result '
             !
             !                        counter = 0
             !                        do j=1,RspFunc_stack(i)%dims(1)
             !                           do k=1,RspFunc_stack(i)%dims(2)
             !                              do l=1,RspFunc_stack(i)%dims(3)
             !                                 counter=counter+1
             !                                 jx = j + RspFunc_stack(i)%first_comp(1)-1
             !                                 kx = k + RspFunc_stack(i)%first_comp(2)-1
             !                                 lx = l + RspFunc_stack(i)%first_comp(3)-1
             !                                 write(lupri,'(X,3I7,F21.7)') jx,kx,lx,&
             !                                      aimag(RspFunc_stack(i)%result(counter))
             !                              enddo
             !                           enddo
             !                        enddo
             !
             !                     end if QuadraticImag

          end if QuadType

          write(lupri,*) 

       end if QuadraticResponse


       CubicResponse: if(print_results .and. RspFunc_stack(i)%order == 4) then 


          write(lupri,*)

          counter = 0

          RspFunc_or_TransMoment: if(RspFunc_stack(i)%trans_moment) then

             if(RspFunc_stack(i)%trans_moment_order == 1) then 

                write(lupri, '(A39)') ' CUBIC SINGLE TRANSITION MOMENT RESULTS'
                write(lupri, '(A39)') " ''''''''''''''''''''''''''''''''''''''"

                ! Single transition moment
                write(lupri,'(1X,A)') '    Comp1  Comp2  Comp3           Result '
                do j=1,RspFunc_stack(i)%dims(1)
                   do k=1,RspFunc_stack(i)%dims(2)
                      do l=1,RspFunc_stack(i)%dims(3)
                         counter=counter+1
                         jx = j + RspFunc_stack(i)%first_comp(1)-1
                         kx = k + RspFunc_stack(i)%first_comp(2)-1
                         lx = l + RspFunc_stack(i)%first_comp(3)-1
                         write(lupri,'(1X,3I7,F21.7)') jx,kx,lx,&
                              real(RspFunc_stack(i)%result(counter))
                      enddo
                   enddo
                enddo

             end if

             if(RspFunc_stack(i)%trans_moment_order == 2) then

                write(lupri, '(A39)') ' CUBIC DOUBLE TRANSITION MOMENT RESULTS'
                write(lupri, '(A39)') " ''''''''''''''''''''''''''''''''''''''"

                ! Double transition moment
                write(lupri,'(1X,A)') '    Comp1  Comp2           Result '
                do j=1,RspFunc_stack(i)%dims(1)
                   do k=1,RspFunc_stack(i)%dims(2)
                      counter=counter+1
                      jx = j + RspFunc_stack(i)%first_comp(1)-1
                      kx = k + RspFunc_stack(i)%first_comp(2)-1
                      write(lupri,'(1X,2I7,F21.7)') jx,kx,&
                           real(RspFunc_stack(i)%result(counter))
                   enddo
                enddo

             end if

          else 

             write(lupri, '(A35)') ' CUBIC RESPONSE RESULTS (REAL PART)'
             write(lupri, '(A35)') " ''''''''''''''''''''''''''''''''''"

             ! Standard response function
             write(lupri,'(1X,A)') '    Comp1  Comp2  Comp3  Comp4           Result '

             do j=1,RspFunc_stack(i)%dims(1)
                do k=1,RspFunc_stack(i)%dims(2)
                   do l=1,RspFunc_stack(i)%dims(3)
                      do m=1,RspFunc_stack(i)%dims(4)
                         counter=counter+1
                         jx = j + RspFunc_stack(i)%first_comp(1)-1
                         kx = k + RspFunc_stack(i)%first_comp(2)-1
                         lx = l + RspFunc_stack(i)%first_comp(3)-1
                         mx = m + RspFunc_stack(i)%first_comp(4)-1
                         write(lupri,'(1X,4I7,F21.7)') jx,kx,lx,mx,&
                              real(RspFunc_stack(i)%result(counter))
                      enddo
                   enddo
                enddo
             enddo


             !                     CubicImag: if(cfg_rsp_imfreqs_specified .or. cfg_rsp_complex) then
             !
             !                        write(lupri,*)
             !                        write(lupri, '(A35)') ' CUBIC RESPONSE RESULTS (IMAG PART)'
             !                        write(lupri, '(A35)') " ''''''''''''''''''''''''''''''''''"
             !
             !                        write(lupri,'(1X,A)') '    Comp1  Comp2  Comp3  Comp4           Result '
             !
             !                        do j=1,RspFunc_stack(i)%dims(1)
             !                           do k=1,RspFunc_stack(i)%dims(2)
             !                              do l=1,RspFunc_stack(i)%dims(3)
             !                                 do m=1,RspFunc_stack(i)%dims(4)
             !                                    counter=counter+1
             !                                    jx = j + RspFunc_stack(i)%first_comp(1)-1
             !                                    kx = k + RspFunc_stack(i)%first_comp(2)-1
             !                                    lx = l + RspFunc_stack(i)%first_comp(3)-1
             !                                    mx = m + RspFunc_stack(i)%first_comp(4)-1
             !                                    write(lupri,'(X,4I7,F21.7)') jx,kx,lx,mx,&
             !                                         aimag(RspFunc_stack(i)%result(counter))
             !                                 enddo
             !                              enddo
             !                           enddo
             !                        enddo
             !
             !                     end if CubicImag

          end if RspFunc_or_TransMoment

          write(lupri,*) 

       end if CubicResponse

    end do RspFunc_loop


    write(lupri,*)
  end subroutine print_response_func



  !> \brief Check whether position A and B
  !> in the response function MyRspFunc are identical - in the sense that 
  !> the OPERATORS on position A and B are identical AND
  !> the FREQUENCIES on position A and B are identical AND
  !> that the same COMPONENTS (e.g. x,y, and z) are requested for both operators.
  !> \author Kasper Kristensen
  !> \date 2010-03
  function Identical_operators(MyRspFunc,A,B) result(identical)
    !> Are the operators identical (true) or not (false)
    logical :: identical
    !> Response function to be investigated
    type(RspFunc), intent(in)      :: MyRspFunc
    !> Positions of the operators to be investigated
    integer, intent(in) :: A,B
    ! Local parameters
    integer :: i,j
    logical :: operators_identical, components_identical, frequencies_identical
    real(realk) :: faR, fbR, faI, fbI

    operators_identical = .false.
    components_identical = .false.
    frequencies_identical = .false.
    identical = .false.

    ! Operators identical?
    if(MyRspFunc%code(A) == MyRspFunc%code(B)) operators_identical = .true.

    ! Components identical?
    if( (MyRspFunc%first_comp(A) == MyRspFunc%first_comp(B)) .and. &
         (MyRspFunc%dims(A) == MyRspFunc%dims(B)) ) then
       components_identical = .true.
    endif

    ! Frequencies identical?
    faR = real(MyRspFunc%freq(A))
    faI = aimag(MyRspFunc%freq(A))
    fbR = real(MyRspFunc%freq(B))
    fbI = aimag(MyRspFunc%freq(B))
    if( ( abs(faR-fbR) < 1e-9 ) .and. ( abs(faI-fbI) < 1e-9 ) ) frequencies_identical = .true.

    ! All true?
    if(operators_identical .and. components_identical &
         .and. frequencies_identical) identical = .true.

    return
  end function Identical_operators
#else
  CONTAINS
  subroutine response_driver_dummy_sub()
    implicit none
  end subroutine response_driver_dummy_sub
#endif

end module response_driver_module
