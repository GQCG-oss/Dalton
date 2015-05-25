!> @file
!> Driver for creating PLT files of orbitals, densities, electrostatic potentials etc.
!> which may be plotted e.g. using the Chimera program.
!> \author Kasper Kristensen (based on program by Branislav Jansik)

module plt_driver_module

  use TYPEDEFTYPE
  use ls_util
  use matrix_operations
  use fundamental
  use BUILDAOBATCH
  use integralinterfaceMod
  use matrix_module
  use matrix_operations
  use typedef
  use typedefTYPE
  use grid_utilities_module
  use davidson_settings
  use IntegralInterfaceMOD!,only:II_Get_overlap
  private
  public :: plt_wrapper, pltinfo_set_default_config, config_plt_input, &
       & config_pltgrid_input, construct_plt_file_driver, &
       & plt_test_suite, calculate_charge, calculate_pplt, calculate_ep,&
       & calculate_density, make_orbitalplot_file
contains

  !> \brief Wrapper for creating PLT files of orbitals, densities, electrostatic potentials etc.
  !> which may be plotted e.g. using the Chimera program.
  !> \author Kasper Kristensen
  !> \date May 2013
  subroutine plt_wrapper(ls,MyPlt)

    implicit none

    !> Integral info
    type(LsItem),intent(inout) :: ls
    !> Information about which PLT information to calculation (see type pltinfo)
    type(pltinfo),intent(inout) :: MyPlt
    integer :: nrow,ncol,nocc,funit, IOS
    integer(kind=long) :: nrow8,ncol8
    Type(Matrix) :: InputMat
    logical :: file_exist


    ! Sanity check 1: Even number of electrons
    nocc=ls%input%molecule%nelectrons/2
    if(nocc*2 /= ls%input%molecule%nelectrons) then
       call lsquit('PLT DRIVER not implemented for unrestricted calculations!',-1)
    end if

    ! Sanity check 2: Density matrix algebra
    if(matrix_type/=mtype_dense) then
       call lsquit('PLT DRIVER only implemented for dense matrices!',-1)
    end if

    ! Sanity check 3: Input file exists
    inquire(file=trim(MyPlt%inputfile),exist=file_exist)
    if(.not. file_exist) then
       write(*,*) 'Input filename: ', trim(MyPlt%inputfile)
       call lsquit('PLT DRIVER: Input file does not exist!',-1)
    end if



    ! read density/orbitals
    ! *********************

    write(ls%lupri,*) 'PLT driver reads input file...'
    funit = 33
    OPEN(UNIT=funit,FILE=trim(MyPlt%inputfile),STATUS='OLD', &
         & FORM='UNFORMATTED',IOSTAT=IOS)
    READ (funit) nrow8,ncol8
    nrow=nrow8
    ncol=ncol8
    write(ls%lupri,*) 'Matrix size: ', nrow, ncol
    call mat_init(InputMat,nrow,ncol)
    READ(funit) InputMat%elms
    CLOSE(funit,STATUS='KEEP',IOSTAT=IOS)

    ! Safety precaution for GC basis
    IF(ls%setting%integraltransformGC)THEN
       call transform_basis_to_GCAO(ls%input%basis)
       ls%setting%integraltransformGC = .FALSE.
    ENDIF

    ! Which calculation?
    call construct_plt_file_driver(MyPlt,InputMat,ls)


    call mat_free(InputMat)

  end subroutine plt_wrapper


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

    ! Zero gridbox info
    MyPlt%nX=0
    MyPlt%nY=0
    MyPlt%nZ=0
    MyPlt%nGRIDPOINTS=0
    MyPlt%deltax=0.0
    MyPlt%deltay=0.0
    MyPlt%deltaz=0.0
    MyPlt%X1=0.0
    MyPlt%Y1=0.0
    MyPlt%Z1=0.0
    MyPlt%buffer=0.0
    MyPlt%gridbox_defined=.false.
    MyPlt%manual_gridbox=.false.
    MyPlt%test=.false.

  end subroutine pltinfo_set_default_config



  !> \brief Read the **PLT section and set the PLT structure accordingly.
  !> Note: The gridbox information of the plt structure is not set here but in config_pltgrid_input.
  !> \author Kasper Kristensen
  !> \date May 2013
  SUBROUTINE config_plt_input(input,output,readword,word,myplt)
    implicit none
    !> Logical for keeping track of when to read
    LOGICAL,intent(inout)                :: READWORD
    !> Logical unit number for LSDALTON.INP
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

          ! Run test case
          ! =============
       case('.TESTPLT')
          myplt%test=.true.


       CASE DEFAULT
          WRITE (output,'(/,3A,/)') ' Keyword "',WORD,&
               & '" not recognized in config_plt_input'
          CALL lsQUIT('Illegal keyword in config_dec_input',output)

       END SELECT PLT_INPUT_INFO

    ENDDO

100 continue

    ! Sanity checks
    if(.not. prop_spec) then
       call lsquit('Error in **PLT input: Property not specified! Possible inputs: &
            & .DENS  .EP  .ORB  .CHARGEDIST ',-1) 
    end if
    if(.not. input_spec) then
       call lsquit('Error in **PLT input: Input file not specified using .INPUT',-1)
    end if
    if(.not. output_spec) then
       call lsquit('Error in **PLT input: Output file not specified using .OUTPUT',-1)
    end if


  END SUBROUTINE config_plt_input



  !> \brief Read the **PLTGRID section and set gridbox related information of the PLT structure
  !> (nX,nY,nZ,nGRIDPOINTS,deltax,deltay,deltaz,X1,Y1,Z1,gridbox_defined) accordingly.
  !> \author Kasper Kristensen
  !> \date May 2013
  SUBROUTINE config_pltgrid_input(input,output,readword,word,myplt)
    implicit none
    !> Logical for keeping track of when to read
    LOGICAL,intent(inout)                :: READWORD
    !> Logical unit number for LSDALTON.INP
    integer,intent(in) :: input
    !> Logical unit number for LSDALTON.OUT
    integer,intent(in) :: output
    !> Word read from input
    character(len=80),intent(inout) :: word
    !> PLT info
    type(pltinfo),intent(inout) :: myplt

    ! Grid box has been defined by input
    MyPlt%gridbox_defined=.true.

    ! There are two possible allowed input forms:
    !
    ! 1. Manually defined gridbox (see type pltinfo for definition of these quantities):
    ! **PLTGRID
    ! .MANUAL
    ! X1 Y1 Z1
    ! deltax deltay deltaz
    ! nX nY nZ
    !
    ! 2. Molecule-adapted gridbox (see DETERMINE_GRIDBOX)
    ! **PLTGRID
    ! .MOLECULE
    ! buffer 
    ! deltax deltay deltaz

    ! First word should be either .MANUAL or .MOLECULE
    READ (input, '(A80)') WORD
    if(trim(word)=='.MANUAL') then
       MyPlt%manual_gridbox=.true.   ! manually defined gridbox
    elseif(trim(word)=='.MOLECULE') then
       MyPlt%manual_gridbox=.false.
    else
       call lsquit('Line following **PLTGRID must be either .MANUAL or .MOLECULE. Please&
            & consult manual!',-1)
    end if

    if(MyPlt%manual_gridbox) then
       read(input,*) MyPlt%X1, MyPlt%Y1, MyPlt%Z1
       read(input,*) MyPlt%deltax, MyPlt%deltay, MyPlt%deltaz
       read(input,*) MyPlt%nX, MyPlt%nY, MyPlt%nZ
       ! Number of gridpoints
       MyPlt%nGRIDPOINTS=MyPlt%nX*MyPlt%nY*MyPlt%nZ

       write(output,*) 
       write(output,*) 'Manual gridbox parameters are:'
       write(output,'(a,3g16.5)') 'X1,Y1,Z1: ', MyPlt%X1, MyPlt%Y1, MyPlt%Z1
       write(output,'(a,3g16.5)') 'deltax,deltay,deltaz: ', MyPlt%deltax, MyPlt%deltay, MyPlt%deltaz
       write(output,'(a,3i16)') 'nX,nY,nZ: ', MyPlt%nX, MyPlt%nY, MyPlt%nZ 
       write(output,'(a,i16)') 'Number of gridpoints: ', MyPlt%nGRIDPOINTS
       write(output,*) 
    else
       read(input,*) MyPlt%buffer
       read(input,*) MyPlt%deltax, MyPlt%deltay, MyPlt%deltaz
       write(output,*) 
       write(output,*) 'Molecule specific gridbox parameters are:'
       write(output,'(a,g16.5)') 'buffer: ', MyPlt%buffer
       write(output,'(a,3g16.5)') 'deltax,deltay,deltaz: ', MyPlt%deltax, MyPlt%deltay, MyPlt%deltaz
       write(output,*) 

    end if

  END SUBROUTINE config_pltgrid_input



  !> \brief Driver for creating PLT files of orbitals, densities, electrostatic potentials etc.
  !> which may be plotted e.g. using the Chimera program.
  !> \author Kasper Kristensen
  !> \date May 2013
  subroutine construct_plt_file_driver(MyPlt,InputMat,ls)
    implicit none
    !> Information about which PLT information to calculate (see type pltinfo)
    type(pltinfo),intent(inout) :: MyPlt
    !> Matrix containing orbitals or density
    type(matrix),intent(in)        :: InputMat
    !> Integral info
    type(lsitem),intent(inout)     :: ls
    integer                        :: I, J, nATOMS,nrow,ncol
    real(4), allocatable       :: ATOMXYZ(:,:)
    type(matrix) :: S,tmpMat,tmpMat2

    nrow = InputMat%nrow
    ncol=InputMat%ncol

    ! Store atomic coordinates in simple array
    ! (We accept simple allocation without mem_alloc here because it is hardcoded real(4))
    nATOMS = ls%setting%MOLECULE(1)%p%natoms
    allocate(ATOMXYZ(3,nATOMS))
    DO I=1,nATOMS 
       DO J=1,3  
          ATOMXYZ(J,I) = real(ls%setting%molecule(1)%p%atom(I)%CENTER(J),4)
       ENDDO
    ENDDO


    ! CHOICE OF GRIDBOX
    ! ******************
    if(MyPLt%gridbox_defined) then ! gridbox was defined by input

       if(MyPlt%manual_gridbox) then
          ! Manually defined gridbox, all parameters were set by input
          write(ls%lupri,*) 'Using PLT gridbox defined manually by input'
       else
          ! Molecule-specific gridbox
          ! deltax,deltay,deltaz, and buffer defined by input.
          write(ls%lupri,*) 'Using molecule-specific PLT gridbox defined by input'
          call DETERMINE_GRIDBOX(natoms,ATOMXYZ,MyPlt)
       end if

    else ! gridbox not defined by input, use default molecule-specific gridbox

       write(ls%lupri,*) 'Using molecule-specific PLT gridbox with default parameters'
       ! Distance of 0.3 a.u. between grid points
       MyPlt%deltax = 0.3_4
       MyPlt%deltay = MyPlt%deltax
       MyPlt%deltaz = MyPlt%deltax
       ! Buffer zone of 6 a.u. around molecule (see details in DETERMINE_GRIDBOX)
       MyPlt%buffer = 6.0_4
       call DETERMINE_GRIDBOX(natoms,ATOMXYZ,MyPlt)
    end if
    write(ls%lupri,*) 'Gridbox parameters are:'
    write(ls%lupri,'(a,3g20.10)') 'X1,Y1,Z1: ', MyPlt%X1, MyPlt%Y1, MyPlt%Z1
    write(ls%lupri,'(a,3g20.10)') 'deltax,deltay,deltaz: ', MyPlt%deltax, MyPlt%deltay, MyPlt%deltaz
    write(ls%lupri,'(a,3i16)') 'nX,nY,nZ: ', MyPlt%nX, MyPlt%nY, MyPlt%nZ 
    write(ls%lupri,'(a,i16)') 'Number of gridpoints: ', MyPlt%nGRIDPOINTS


    ! Sanity check
    if(MyPlt%nGRIDPOINTS==0) then
       call lsquit('Zero gridpoints: Cannot construct PLT file!',-1)
    end if

    ! Which calculation?
    ! ******************

    select case (myplt%frmt)

       ! Density
    case('DENS')
       write(ls%lupri,*) 'Writing density distribution plt file...'
       call calculate_density(trim(MyPlt%outputfile),InputMat,ls,natoms,ATOMXYZ,myplt)

       ! Electrostatic potential
    case('EP')
       write(ls%lupri,*) 'Writing electrostatic potential plt file...'
       call calculate_ep(trim(MyPlt%outputfile),InputMat,ls,natoms,ATOMXYZ,myplt)

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
       call calculate_pplt(trim(MyPlt%outputfile),MyPlt%iorb,InputMat,ls,natoms,ATOMXYZ,myplt)

       ! Charge distribution
    case('CHARGEDIST')
       write(6,*) 'Writing charge distribution plt file...'
       call calculate_charge(trim(MyPlt%outputfile),Myplt%iorb,myplt%jorb,InputMat,&
            & ls,natoms,ATOMXYZ,myplt)

    case default
       write(ls%lupri,*) 'PLT driver unknown input format: ', myplt%frmt
       call lsquit('PLT driver unknown input format',-1)
    end select

    deallocate(ATOMXYZ)

    ! Tets calculate grid points
    if(MyPlt%test) then
       call plt_test_suite(MyPlt,ls%lupri)
    end if

  end subroutine construct_plt_file_driver

  
  !> \brief Test PLT driver.
  !> \author Kasper Kristensen
  !> \date May 2013
  subroutine plt_test_suite(MyPlt,lupri)
    implicit none
    !> Information about which PLT information to calculation (see type pltinfo)
    type(pltinfo),intent(inout) :: MyPlt
    !> LSDALTON.OUT unit number
    integer,intent(in) :: lupri
    integer :: iunit
    integer(4) :: II,TD,reclen,nX,nY,nZ,i
    real(4) :: Z1,Zn,Y1,Yn,X1,Xn
    real(4),pointer :: values(:)

    write(lupri,*) 'Testing PLT driver...'

    !> Read in file generated by PLT driver
    IUNIT=111
    reclen=(MyPlt%nGRIDPOINTS + 11)
    reclen=reclen*4
    allocate(values(MyPlt%ngridpoints))

    OPEN (IUNIT,FILE=trim(MyPlt%outputfile),STATUS='OLD',FORM='UNFORMATTED',&
         & ACCESS='DIRECT',RECL=reclen)
    READ(IUNIT,REC=1) II,TD,nZ,nY,nX,Z1,Zn,Y1,Yn,X1,Xn,values
    CLOSE(IUNIT,STATUS='KEEP')

    ! Check that number of points in X,Y,Z directions are consistent with PLT structure
    if(MyPlt%nX/=nX .or. MyPlt%nY/=nY .or. MyPlt%nZ/=nZ) then
       print *, 'Plt structure : ', MyPlt%nX, MyPlt%nY, MyPlt%nZ
       print *, 'Read from file: ', nX, nY, nZ
       call lsquit('plt_test_suite: Something wrong number of grid points!',-1)
    end if

    ! Calculate sum of values and min/max values
    write(lupri,*) 'PLT TEST: Sum of values = ', sum(values)
    write(lupri,*) 'PLT TEST: Minimum value = ', minval(values)
    write(lupri,*) 'PLT TEST: Maximum value = ', maxval(values)
    deallocate(values)


  end subroutine plt_test_suite

  !> \brief Calculates charge distribution between two orbitals at grid points.
  !> KK: Whoever is responsible, please document this properly!
  subroutine calculate_charge(filename,iorb,jorb,dCMO,ls,natoms,ATOMXYZ,MyPlt)
    implicit none
    type(matrix)                   :: dCMO
    type(lsitem)                   :: ls
    Character*(*)                  :: filename
    integer,intent(in)             :: natoms
    integer                        :: nX,nY,nZ
    real(4), intent(in)            :: ATOMXYZ(3,nATOMS)
    !> PLT info, including grid box info
    type(pltinfo),intent(in) :: MyPlt
    !Local
    integer                        :: I,J,P,Q,Xg,Yg,Zg,nGRIDPOINTS
    integer                        :: orbnr, nORBITALS,nBASIS,iorb,jorb,iunit,ig
    real(4)                        :: X,Y,Z,X1,Y1,Z1,Xgrid, Ygrid,Zgrid,deltax, deltay, deltaz 
    real(4), allocatable           :: moorb(:),CMO(:,:),GAO(:),tmp(:)
    real(4)                        :: fZ1,fY1,fX1,fZn,fYn,fXn
    integer                        :: II,TD, ijk
    integer                        :: reclen
    real(4)                        :: br
    double precision, external     :: SDOT
    br=real(bohr_to_angstrom,4)

    nORBITALS = dCMO%ncol
    nBASIS    = dCMO%nrow

    ! Copy grid box info from MyPlt structure for easier overview
    nX = MyPlt%nX
    nY = MyPlt%nY
    nZ = MyPlt%nZ
    deltax = MyPlt%deltax
    deltay = MyPlt%deltay
    deltaz = MyPlt%deltaz
    X1 = MyPlt%X1
    Y1 = MyPlt%Y1
    Z1 = MyPlt%Z1
    nGRIDPOINTS = MyPlt%nGRIDPOINTS

    allocate(moorb(nGRIDPOINTS), GAO(nBASIS),CMO(nBASIS,nORBITALS),tmp(nBASIS*nORBITALS)) !MO orbitals
    moorb = 0.0; ijk = 0

    tmp=real(dCMO%elms,4); call SCOPY(nBASIS*nORBITALS,tmp,1,CMO,1)
    deallocate(tmp)

    DO Zg = 1,nZ 
       Zgrid = Z1 + deltaz*(Zg-1) 
       DO Yg = 1,nY 
          Ygrid = Y1 + deltay*(Yg-1) 
          DO Xg = 1,nX
             Xgrid = X1 + deltax*(Xg-1) 



             GAO = 0.0_4
             orbnr=1

             DO I=1,natoms
                ! X-, Y-, and Z-distances from gridpoint to atom I
                X = Xgrid - ATOMXYZ(1,I)
                Y = Ygrid - ATOMXYZ(2,I)
                Z = Zgrid - ATOMXYZ(3,I)
                call determine_orbitals(I,GAO,orbnr, nBASIS, X,Y,Z,ls)
             ENDDO


             ijk = ijk + 1

             !moorb(ijk) = real(SDOT(nBASIS,CMO(:,iorb),1,GAO,1) * SDOT(nBASIS,CMO(:,jorb),1,GAO,1))
             moorb(ijk) = dot_product(CMO(:,iorb),GAO) * dot_product(CMO(:,jorb),GAO)


          ENDDO
       ENDDO
    ENDDO

    deallocate(GAO,CMO)


    !Write plt file

    II=3; TD=200
    fZ1=Z1*br; fZn = fZ1 + br*deltaz*(nZ-1);
    fY1=Y1*br; fYn = fY1 + br*deltay*(nY-1);
    fX1=X1*br; fXn = fX1 + br*deltax*(nX-1);

    IUNIT=111
    reclen=(nGRIDPOINTS + 11)
    reclen=reclen*4

    OPEN (IUNIT,FILE=trim(filename),STATUS='REPLACE',FORM='UNFORMATTED',&
         & ACCESS='DIRECT',RECL=reclen)

    WRITE(IUNIT,REC=1) int(II,4),int(TD,4),int(nZ,4),int(nY,4),int(nX,4),fZ1,fZn,fY1,fYn,fX1,fXn,moorb

    CLOSE(IUNIT,STATUS='KEEP')

    deallocate(moorb)


  end subroutine calculate_charge


  !> \brief Calculates orbital value at gridpoints.
  !> KK: Whoever is responsible, please document this properly!
  subroutine calculate_pplt(filename,iorb,dCMO,ls,natoms,ATOMXYZ,myplt)
    implicit none
    type(matrix)                   :: dCMO
    type(lsitem)                   :: ls
    Character*(*)                  :: filename
    integer,intent(in)             :: natoms
    integer                        :: nX,nY,nZ
    real(4), intent(in)            :: ATOMXYZ(3,nATOMS)
    !> PLT info, including grid box info
    type(pltinfo),intent(in) :: MyPlt

    !Local
    integer                        :: I,J,P,Q,Xg,Yg,Zg,nGRIDPOINTS
    integer                        :: orbnr, nORBITALS,nBASIS,iorb,iunit,ig
    real(4)                        :: X,Y,Z,X1,Y1,Z1,Xgrid, Ygrid,Zgrid,deltax, deltay, deltaz 
    real(4),     pointer           :: GAO(:)
    real(4), allocatable           :: moorb(:),CMO(:)
    real(4)                        :: fZ1,fY1,fX1,fZn,fYn,fXn
    integer                        :: II,TD, ijk
    integer                        :: reclen
    real(4)                        :: br
    double precision, external     :: SDOT
    br=real(bohr_to_angstrom,4)

    nORBITALS = dCMO%ncol
    nBASIS    = dCMO%nrow

    ! Copy grid box info from MyPlt structure for easier overview
    nX = MyPlt%nX
    nY = MyPlt%nY
    nZ = MyPlt%nZ
    deltax = MyPlt%deltax
    deltay = MyPlt%deltay
    deltaz = MyPlt%deltaz
    X1 = MyPlt%X1
    Y1 = MyPlt%Y1
    Z1 = MyPlt%Z1
    nGRIDPOINTS = MyPlt%nGRIDPOINTS

    allocate(moorb(nGRIDPOINTS),CMO(nBASIS)) !MO orbitals
    moorb = 0; ijk = 0

    CMO=real((dCMO%elms(1+(iorb-1)*nBASIS:iorb*nBASIS)),4)

    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(Zg,Zgrid,Yg,Ygrid,Xg,Xgrid,orbnr,I,X,Y,Z,ijk,GAO)

    allocate(GAO(nBASIS))

    !$OMP DO SCHEDULE(STATIC)
    DO Zg = 1,nZ 
       Zgrid = Z1 + deltaz*(Zg-1) 
       DO Yg = 1,nY 
          Ygrid = Y1 + deltay*(Yg-1) 
          DO Xg = 1,nX
             Xgrid = X1 + deltax*(Xg-1) 



             GAO = 0.0_4
             orbnr=1

             DO I=1,natoms
                ! X-, Y-, and Z-distances from gridpoint to atom I
                X = Xgrid - ATOMXYZ(1,I)
                Y = Ygrid - ATOMXYZ(2,I)
                Z = Zgrid - ATOMXYZ(3,I)
                call determine_orbitals(I,GAO,orbnr, nBASIS, X,Y,Z,ls)
             ENDDO

             ijk = Xg + (Yg-1)*nX + (Zg-1)*nX*nY

             !moorb(ijk) = real(SDOT(nBASIS,CMO,1,GAO,1))
             moorb(ijk) = dot_product(CMO,GAO)


          ENDDO
       ENDDO
    ENDDO
    !$OMP END DO

    deallocate(GAO)
    !$OMP END PARALLEL



    !Write plt file

    II=3; TD=200
    fZ1=Z1*br; fZn = fZ1 + br*deltaz*(nZ-1);
    fY1=Y1*br; fYn = fY1 + br*deltay*(nY-1);
    fX1=X1*br; fXn = fX1 + br*deltax*(nX-1);

    IUNIT=111
    reclen=(nGRIDPOINTS + 11)
    reclen=reclen*4

    OPEN (IUNIT,FILE=trim(filename),STATUS='REPLACE',FORM='UNFORMATTED',&
         & ACCESS='DIRECT',RECL=reclen)

    WRITE(IUNIT,REC=1) int(II,4),int(TD,4),int(nZ,4),int(nY,4),int(nX,4),fZ1,fZn,fY1,fYn,fX1,fXn,moorb

    CLOSE(IUNIT,STATUS='KEEP')

    deallocate(moorb,CMO)


  end subroutine calculate_pplt


  !> \brief Calculates electrostatic potential at gridpoints.
  !> KK: Whoever is responsible, please document this properly!
  subroutine calculate_ep(filename,D,ls,natoms,ATOMXYZ,MyPlt)
    implicit none
    type(matrix)                   :: D
    type(matrix)                   :: dEPint
    type(lsitem)                   :: ls
    TYPE(LSSETTING)                :: psetting
    Character*(*)                  :: filename
    integer,intent(in)             :: natoms
    integer                        :: nX,nY,nZ
    real(4), intent(in)            :: ATOMXYZ(3,nATOMS)
    !> PLT info, including grid box info
    type(pltinfo),intent(in) :: MyPlt
    !Local
    real(realk), allocatable       :: R(:,:), emoorb(:),moorb1(:)
    integer                        :: I,Q,Xg,Yg,Zg,nGRIDPOINTS
    integer                        :: orbnr, nORBITALS,nBASIS,iunit,ig
    real(realk)                    :: X,Y,Z,X1,Y1,Z1,Xgrid, Ygrid,Zgrid,deltax, deltay, deltaz 
    real(4), allocatable           :: moorb(:)
    real(realk)                    :: fZ1,fY1,fX1,fZn,fYn,fXn,epnuc
    integer                        :: II,TD, ijk, icount
    integer                        :: reclen
    real(realk)                    :: br
    br=real(bohr_to_angstrom,4)

    nORBITALS = D%ncol
    nBASIS    = D%nrow

    ! Copy grid box info from MyPlt structure for easier overview
    nX = MyPlt%nX
    nY = MyPlt%nY
    nZ = MyPlt%nZ
    deltax = MyPlt%deltax
    deltay = MyPlt%deltay
    deltaz = MyPlt%deltaz
    X1 = MyPlt%X1
    Y1 = MyPlt%Y1
    Z1 = MyPlt%Z1
    nGRIDPOINTS = MyPlt%nGRIDPOINTS

    allocate(moorb1(nGRIDPOINTS),emoorb(nGRIDPOINTS),R(3,nGRIDPOINTS)) !MO orbitals
    emoorb = 0_realk;

    write(6,*) 'Computing ', nGRIDPOINTS, ' gridpoints'

    !   Nuclear contribution

    !$OMP PARALLEL DO SCHEDULE(STATIC)  DEFAULT(SHARED) PRIVATE(Zg,Zgrid,Yg,Ygrid,Xg,Xgrid,epnuc,I,X,Y,Z,Q,ijk)
    DO Zg = 1,nZ 
       Zgrid = Z1 + deltaz*(Zg-1) 
       DO Yg = 1,nY 
          Ygrid = Y1 + deltay*(Yg-1) 
          DO Xg = 1,nX
             Xgrid = X1 + deltax*(Xg-1) 


             ijk = Xg + (Yg-1)*nX + (Zg-1)*nX*nY

             epnuc = 0_realk

             DO I=1,natoms
                ! X-, Y-, and Z-distances from gridpoint to atom I
                X = Xgrid - ATOMXYZ(1,I)
                Y = Ygrid - ATOMXYZ(2,I)
                Z = Zgrid - ATOMXYZ(3,I)
                Q = ls%setting%molecule(1)%p%atom(I)%Charge
                epnuc = epnuc + (Q/SQRT(X*X + Y*Y + Z*Z))
             ENDDO

             R(1,ijk) = Xgrid
             R(2,ijk) = Ygrid
             R(3,ijk) = Zgrid

             moorb1(ijk) = epnuc 

          ENDDO
       ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    !   Electronic contribution
    call II_get_ep_integrals3(ls%lupri,ls%luerr,ls%setting,R,nGRIDPOINTS,D,emoorb)
    !   Add the two
    deallocate(R)

    allocate(moorb(nGRIDPOINTS)) !MO orbitals

    moorb = real(moorb1 - emoorb,4)

    deallocate(emoorb)
    deallocate(moorb1)

    !Write plt file

    II=3; TD=200
    fZ1=Z1*br; fZn = fZ1 + br*deltaz*(nZ-1);
    fY1=Y1*br; fYn = fY1 + br*deltay*(nY-1);
    fX1=X1*br; fXn = fX1 + br*deltax*(nX-1);

    IUNIT=111
    reclen=(nGRIDPOINTS + 11)
    reclen=reclen*4

    OPEN (IUNIT,FILE=trim(filename),STATUS='REPLACE',FORM='UNFORMATTED',&
         & ACCESS='DIRECT',RECL=reclen)


    WRITE(IUNIT,REC=1) int(II,4),int(TD,4),int(nZ,4),int(nY,4),int(nX,4),&
         & real(fZ1,4),real(fZn,4),real(fY1,4),real(fYn,4),real(fX1,4),&
         & real(fXn,4),moorb
    CLOSE(IUNIT,STATUS='KEEP')

    deallocate(moorb)


  end subroutine calculate_ep


  !> \brief Calculates electron density at gridpoints.
  !> KK: Whoever is responsible, please document this properly!
  subroutine calculate_density(filename,dD,ls,natoms,ATOMXYZ,myplt)
    implicit none
    type(matrix)                   :: dD
    type(lsitem)                   :: ls
    Character*(*)                  :: filename
    integer,intent(in)             :: natoms
    integer                        :: nX,nY,nZ
    real(4), intent(in)            :: ATOMXYZ(3,nATOMS)
    !> PLT info, including grid box info
    type(pltinfo),intent(in) :: MyPlt
    !Local
    integer                        :: I,J,P,Q,Xg,Yg,Zg,nGRIDPOINTS
    integer                        :: orbnr,nBASIS,iorb,iunit,ig
    real(4)                        :: X,Y,Z,X1,Y1,Z1,Xgrid, Ygrid,Zgrid,deltax, deltay, deltaz 
    real(4), allocatable           :: GAO(:,:),D(:),tmp(:,:)
    real(4), allocatable           :: moorb(:),XgridB(:),YgridB(:),ZgridB(:)
    real(4)                        :: fZ1,fY1,fX1,fZn,fYn,fXn
    integer                        :: II,TD, ijk, klm, nLeftover
    integer                        :: reclen
    real(4)                        :: br
    double precision, external     :: SDOT
    br=real(bohr_to_angstrom,4)

    nBASIS    = dD%nrow

    ! Copy grid box info from MyPlt structure for easier overview
    nX = MyPlt%nX
    nY = MyPlt%nY
    nZ = MyPlt%nZ
    deltax = MyPlt%deltax
    deltay = MyPlt%deltay
    deltaz = MyPlt%deltaz
    X1 = MyPlt%X1
    Y1 = MyPlt%Y1
    Z1 = MyPlt%Z1
    nGRIDPOINTS = MyPlt%nGRIDPOINTS

    write(*,*)  'calculate_density2 computing',nGRIDPOINTS, ' gridpoints' 

    !allocate
    allocate(moorb(nGRIDPOINTS), D(nBASIS*nBASIS))
    allocate(XgridB(nBASIS),YgridB(nBASIS),ZgridB(nBASIS))
    allocate(GAO(nBASIS,nBASIS),tmp(nBASIS,nBASIS))

    moorb = 0e0; ijk = 0; klm = 0;
    D = real(dD%elms,4)


    DO Zg = 1,nZ 
       Zgrid = Z1 + deltaz*(Zg-1) 
       DO Yg = 1,nY 
          Ygrid = Y1 + deltay*(Yg-1) 
          DO Xg = 1,nX
             Xgrid = X1 + deltax*(Xg-1) 


             ijk = ijk + 1
             klm = klm + 1

             XgridB(klm)=Xgrid
             YgridB(klm)=Ygrid
             ZgridB(klm)=Zgrid


             IF (mod(ijk,nBASIS).eq. 0) THEN


                GAO = 0E0_4

                !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(klm,I,orbnr,X,Y,Z)
                DO klm=1,nBASIS
                   orbnr=1

                   DO I=1,natoms
                      ! X-, Y-, and Z-distances from gridpoint to atom I
                      X = XgridB(klm) - ATOMXYZ(1,I)
                      Y = YgridB(klm) - ATOMXYZ(2,I)
                      Z = ZgridB(klm) - ATOMXYZ(3,I)
                      call determine_orbitals(I,GAO(:,klm),orbnr, nBASIS, X,Y,Z,ls)
                   ENDDO
                ENDDO
                !$OMP END PARALLEL DO


                CALL SSYMM('L','L',nBASIS,nBASIS,1.0,D,nBASIS,GAO,nBASIS,0.0,tmp,nBASIS)


                !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(klm)
                DO klm=1,nBASIS
                   moorb(ijk-nBASIS+klm) = 2.0*dot_product(GAO(:,klm),tmp(:,klm))
                ENDDO
                !$OMP END PARALLEL DO

                klm = 0

             ENDIF


          ENDDO
       ENDDO
    ENDDO

    !Compute leftover gridpoints
    nLeftover = mod(nGRIDPOINTS,nBASIS)
    IF (nLeftover.ne. 0) THEN

       GAO=0E0_4

       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(klm,I,orbnr,X,Y,Z)
       DO klm=1,nLeftover
          orbnr=1

          DO I=1,natoms
             ! X-, Y-, and Z-distances from gridpoint to atom I
             X = XgridB(klm) - ATOMXYZ(1,I)
             Y = YgridB(klm) - ATOMXYZ(2,I)
             Z = ZgridB(klm) - ATOMXYZ(3,I)
             call determine_orbitals(I,GAO(:,klm),orbnr, nBASIS, X,Y,Z,ls)
          ENDDO
       ENDDO
       !$OMP END PARALLEL DO

       CALL SSYMM('L','L',nBASIS,nLeftover,1.0,D,nBASIS,GAO,nBASIS,0.0,tmp,nBASIS)

       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(klm)
       DO klm=1,nLeftover
          moorb(ijk-nLeftover+klm) = 2.0*dot_product(GAO(:,klm),tmp(:,klm))
       ENDDO
       !$OMP END PARALLEL DO

    ENDIF


    deallocate(GAO,tmp,XgridB,YgridB,ZgridB,D)

    !Write plt file

    II=3; TD=200
    fZ1=Z1*br; fZn = fZ1 + br*deltaz*(nZ-1);
    fY1=Y1*br; fYn = fY1 + br*deltay*(nY-1);
    fX1=X1*br; fXn = fX1 + br*deltax*(nX-1);

    IUNIT=111
    reclen=(nGRIDPOINTS + 11)
    reclen=reclen*4

    OPEN (IUNIT,FILE=trim(filename),STATUS='REPLACE',FORM='UNFORMATTED',&
         & ACCESS='DIRECT',RECL=reclen)

    WRITE(IUNIT,REC=1) int(II,4),int(TD,4),int(nZ,4),int(nY,4),int(nX,4),fZ1,fZn,fY1,fYn,fX1,fXn,moorb

    CLOSE(IUNIT,STATUS='KEEP')

    deallocate(moorb)


  end subroutine calculate_density


  !> \brief Make PLT files for specific orbitals
  !> Ida, please document properly!
  subroutine make_orbitalplot_file(CMO,CFG,ls,MyPlt)
    implicit none
    !> file containing localized orbitals
    type(matrix) :: CMO 
    type(RedSpaceItem) :: CFG
    type(lsitem) :: ls
    !> Info for PLT generation
    type(pltinfo),intent(inout) :: MyPlt
    integer :: nocc


    IF(ls%setting%integraltransformGC)THEN
       call transform_basis_to_GCAO(ls%input%basis)
       ls%setting%integraltransformGC = .FALSE.
    ENDIF

    nocc=int(ls%input%molecule%nelectrons/2)
    MyPlt%frmt = 'ORB'

    select case (trim(CFG%plt_orbital))
    case('LEASTL')
       write(ls%lupri,*) 'Writing plt file for the least local occupied and virtual orbitals...'
       MyPlt%iorb=CFG%leastl_occ
       Myplt%outputfile = 'LeastLocal_occ.plt'
       if (MyPlt%iorb>0) call construct_plt_file_driver(MyPlt,CMO,ls)

       MyPlt%iorb=CFG%leastl_virt
       Myplt%outputfile = 'LeastLocal_virt.plt'
       if (MyPlt%iorb>0) call construct_plt_file_driver(MyPlt,CMO,ls)
    case('MOSTL')
       write(ls%lupri,*) 'Writing plt file for the most local occupied and virtual orbitals...'
       MyPlt%iorb=CFG%mostl_occ
       Myplt%outputfile = 'MostLocal_occ.plt'
       if (MyPlt%iorb>0) call construct_plt_file_driver(MyPlt,CMO,ls)

       Myplt%outputfile = 'MostLocal_virt.plt'
       MyPlt%iorb=CFG%mostl_virt
       if (MyPlt%iorb>0) call construct_plt_file_driver(MyPlt,CMO,ls)
    case('ALL')
       write(ls%lupri,*) 'Writing plt file for the least local occupied and virtual orbitals...'
       MyPlt%iorb=CFG%leastl_occ
       Myplt%outputfile = 'LeastLocal_occ.plt'
       if (MyPlt%iorb>0) call construct_plt_file_driver(MyPlt,CMO,ls)

       MyPlt%iorb=CFG%leastl_virt
       Myplt%outputfile = 'LeastLocal_virt.plt'
       if (MyPlt%iorb>0) call construct_plt_file_driver(MyPlt,CMO,ls)

       write(ls%lupri,*) 'Writing plt file for the most local occupied and virtual orbitals...'
       MyPlt%iorb=CFG%mostl_occ
       Myplt%outputfile = 'MostLocal_occ.plt'
       if (MyPlt%iorb>0) call construct_plt_file_driver(MyPlt,CMO,ls)

       MyPlt%iorb=CFG%mostl_virt
       Myplt%outputfile = 'MostLocal_virt.plt'
       if (MyPlt%iorb>0) call construct_plt_file_driver(MyPlt,CMO,ls)

    case default
       write(ls%lupri,*) '========== SOMETHING WRONG WHEN MAKING .plt FILE ================'
       write(ls%lupri,*) 'None of the valid options chosen. Use the options'
       write(ls%lupri,*) ' MOSTL,LEASTL or ALL below keyword .ORBITAL PLOT'
       write(ls%lupri,*) ' If you want to plot other than the above orbitals, consider using'
       write(ls%lupri,*) ' **PLT option where you may specify any orbital index.'
       write(ls%lupri,*) '================================================================='

       write(*,*) '========== SOMETHING WRONG WHEN MAKING .plt FILE ================'
       write(*,*) 'None of the valid options chosen. Use the options'
       write(*,*) ' MOSTL,LEASTL or ALL below keyword .ORBITAL PLOT'
       write(*,*) ' If you want to plot other than the above orbitals, consider using'
       write(*,*) ' **PLT option where you may specify any orbital index.'
       write(*,*) '================================================================='
    end select



  end subroutine make_orbitalplot_file


end module plt_driver_module
