!> @file
!> Driver for creating PLT files of orbitals, densities, electrostatic potentials etc.
!> which may be plotted e.g. using the Chimera program.
!> \author Kasper Kristensen (based on program by Branislav Jansik)

module plt_driver_module

    use TYPEDEFTYPE
    use ls_util
    use print_moorb_grid_mod
    use matrix_operations
    use IntegralInterfaceMOD,only:II_Get_overlap

contains

  !> \brief Driver for creating PLT files of orbitals, densities, electrostatic potentials etc.
  !> which may be plotted e.g. using the Chimera program.
  !> \author Kasper Kristensen
  !> \date May 2013
  subroutine contruct_plt_file_driver(ls,MyPlt)

    implicit none

    !> Integral info
    type(LsItem),intent(inout) :: ls
    !> Information about which PLT information to calculation (see type pltinfo)
    type(pltinfo),intent(in) :: MyPlt
    integer :: nocc,nrow,ncol,funit, IOS,I,J,natoms
    Type(Matrix) :: InputMat,S,tmpMat,tmpMat2
    real(4), allocatable       :: ATOMXYZ(:,:)


    ! Sanity check 1: Even number of electrons
    nocc=ls%input%molecule%nelectrons/2
    if(nocc*2 /= ls%input%molecule%nelectrons) then
       call lsquit('PLT DRIVER not implemented for unrestricted calculations!',-1)
    end if

    ! Sanity check 2: Density matrix algebra
    if(matrix_type/=mtype_dense) then
       call lsquit('PLT DRIVER only implemented for dense matrices!',-1)
    end if



    ! Store atomic coordinates in simple array
    ! ****************************************
    nATOMS = ls%setting%MOLECULE(1)%p%natoms
    ! Since this is the only place with hardcoded real(4) 
    ! we allow the use of allocate rather than mem_alloc....
    allocate(ATOMXYZ(3,nATOMS))
    DO I=1,nATOMS 
       DO J=1,3
          ATOMXYZ(J,I) = real(ls%setting%molecule(1)%p%atom(I)%CENTER(J),4)
       ENDDO
    ENDDO



    ! read density/orbitals
    ! *********************
    write(6,*) 'PLT driver reads input file...'
    funit = 33
    OPEN(UNIT=funit,FILE=trim(MyPlt%inputfile),STATUS='OLD', &
         & FORM='UNFORMATTED',IOSTAT=IOS)
    READ (funit) nrow,ncol
    write(6,*) 'Matrix size: ', nrow, ncol
    call mat_init(InputMat,nrow,ncol)
    READ(funit) InputMat%elms
    CLOSE(funit,STATUS='KEEP',IOSTAT=IOS)

    ! Safety precaution for GC basis
    IF(ls%setting%integraltransformGC)THEN
       call transform_basis_to_GCAO(ls%input%basis)
       ls%setting%integraltransformGC = .FALSE.
    ENDIF



    ! Which calculation?
    ! ******************

    select case (trim(MyPlt%frmt))

       ! Density
    case('DENS')
       write(ls%lupri,*) 'Writing density distribution plt file...'
       call calculate_density2(trim(MyPlt%outputfile),InputMat,ls,natoms,ATOMXYZ)

       ! Electrostatic potential
    case('EP')
       write(ls%lupri,*) 'Writing electrostatic potential plt file...'
       call calculate_ep3(trim(MyPlt%outputfile),InputMat,ls,natoms,ATOMXYZ)

       ! Orbital 
    case('ORB')

       ! Check whether orbitals are orthogonal
       call mat_init(S,nrow,ncol)
       CALL II_get_overlap(ls%lupri,ls%luerr,ls%setting,S)
       call mat_init(tmpMat,nrow,ncol)
       call mat_init(tmpMat2,nrow,ncol)
       call mat_mul(InputMat,S,'t','n',1E0_realk,0E0_realk,tmpMat)
       call mat_mul(tmpMat,InputMat,'n','n',1E0_realk,0E0_realk,S)
       !S should now be the identity
       call mat_identity(tmpMat)
       call mat_add(1E0_realk,S,-1E0_realk,tmpMat,tmpMat2)
       print*,'mat_sqnorm(C^T*S*C - I) =',mat_sqnorm2(tmpMat2)/tmpMat2%nrow

       ! KK: I think we should just print a warning and not exit here.
       ! Otherwise it is not possible to plot nonorthogonal orbitals...
       IF(ABS(mat_sqnorm2(tmpMat2)/tmpMat2%nrow).GT.1.0E-10_realk)THEN
          write(ls%lupri,*) 'WARNING WARNING WARNING'
          write(ls%lupri,*) 'The input orbitals are not orthogonal!'
          write(ls%lupri,*) 'WARNING WARNING WARNING'
          write(*,*) 'WARNING WARNING WARNING'
          write(*,*) 'The input orbitals are not orthogonal!'
          write(*,*) 'WARNING WARNING WARNING'
       ENDIF
       call mat_free(S)
       call mat_free(tmpMat)
       call mat_free(tmpMat2)

       write(6,*) 'Writing orbital plt file...'
       call calculate_pplt(trim(MyPlt%outputfile),MyPlt%iorb,InputMat,ls,natoms,ATOMXYZ)

       ! Charge distribution
    case('CHARGEDIST')
       write(6,*) 'Writing charge distribution plt file...'
       call calculate_charge(trim(MyPlt%outputfile),MyPlt%iorb,MyPlt%jorb,InputMat,ls,natoms,ATOMXYZ)

    case default
       write(ls%lupri,*) 'PLT driver unknown input format: ', MyPlt%frmt
       call lsquit('PLT driver unknown input format',-1)
    end select


    call mat_free(InputMat)
    deallocate(ATOMXYZ)

  end subroutine contruct_plt_file_driver


  !> \brief Set default configurations for pltinfo type.
  !> \author Kasper Kristensen
  !> \date May 2013
  subroutine pltinfo_set_default_config(MyPlt)
    implicit none
    !> Information about which PLT information to calculation (see type pltinfo)
    type(pltinfo),intent(inout) :: MyPlt
    integer :: l,i

    ! Length characters hardcoded to 80 (see type pltinfo)
    l=80

    ! Make characters blank (just in case)
    do i=1,l
       MyPlt%inputfile(i:i) = ' '
       MyPlt%outputfile(i:i) = ' '
       MyPlt%frmt(i:i) = ' '
    end do

    ! Zero orbital indices
    MyPlt%iorb = 0
    MyPlt%jorb = 0

  end subroutine pltinfo_set_default_config



  !> \brief Read the **DEC or **CC input section in DALTON.INP and set 
  !> configuration structure accordingly.
  !> \author Kasper Kristensen
  !> \date September 2010
  SUBROUTINE config_plt_input(input,output,readword,word,myplt)
    implicit none
    !> Logical for keeping track of when to read
    LOGICAL,intent(inout)                :: READWORD
    !> Logical unit number for DALTON.INP
    integer,intent(in) :: input
    !> Logical unit number for LSDALTON.OUT
    integer,intent(in) :: output
    !> Word read from input
    character(len=80),intent(inout) :: word
    !> PLT info
    type(pltinfo),intent(inout) :: myplt
    logical :: prop_spec, input_spec, output_spec

    ! Check that everything has been set
    prop_spec=.false.
    input_spec = .false.
    output_spec = .false.


    DO

       IF(READWORD) THEN
          READ (input, '(A80)') WORD
          call capitalize_string(word)
          READWORD=.TRUE.
       ENDIF

       IF ((WORD(1:1) .EQ. '!') .OR. (WORD(1:1) .EQ. '#')) CYCLE

       IF(WORD(1:2) .EQ. '**') THEN
          READWORD=.FALSE.
          goto 100
       ENDIF

       IF(WORD(1:13) == '*END OF INPUT') THEN
          goto 100
       END IF



       PLT_INPUT_INFO: SELECT CASE(trim(WORD))


          ! POSSIBLE PROPERTIES TO CALCULATE
          ! ================================
          ! See type pltinfo

          ! Density
       case('.DENS')
          myplt%frmt='DENS'
          prop_spec=.true.

          ! Electrostatic potentian
       case('.EP')
          myplt%frmt='EP'
          prop_spec=.true.

          ! Single orbital
       case('.ORB')
          myplt%frmt='ORB'
          read(input,*) myplt%iorb
          prop_spec=.true.

          ! Charge distribution between two orbitals
       case('.CHARGEDIST')
          myplt%frmt='CHARGEDIST'
          read(input,*) myplt%iorb, myplt%jorb
          prop_spec=.true.


          ! Input file
          ! ==========
       case('.INPUT')
          read(input,*) myplt%inputfile
          input_spec=.true.


          ! Output file
          ! ==========
       case('.OUTPUT')
          read(input,*) myplt%outputfile
          output_spec=.true.


       CASE DEFAULT
          WRITE (output,'(/,3A,/)') ' Keyword "',WORD,&
               & '" not recognized in config_plt_input'
          CALL lsQUIT('Illegal keyword in config_dec_input',output)

       END SELECT PLT_INPUT_INFO

    ENDDO

    100 continue

    ! Sanity checks
    if(.not. prop_spec) then
       call lsquit('Error in **PLOT input: Property not specified! Possible inputs: &
            & .DENS  .EP  .ORB  .CHARGEDIST ',-1) 
    end if
    if(.not. input_spec) then
       call lsquit('Error in **PLOT input: Input file not specified using .INPUT',-1)
    end if
    if(.not. output_spec) then
       call lsquit('Error in **PLOT input: Output file not specified using .OUTPUT',-1)
    end if


  END SUBROUTINE config_plt_input



end module plt_driver_module
