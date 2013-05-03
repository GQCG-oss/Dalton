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
  use IntegralInterfaceMOD,only:II_Get_overlap

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
    type(pltinfo),intent(in) :: MyPlt
    integer :: nocc,nrow,ncol,funit, IOS
    Type(Matrix) :: InputMat


    ! Sanity check 1: Even number of electrons
    nocc=ls%input%molecule%nelectrons/2
    if(nocc*2 /= ls%input%molecule%nelectrons) then
       call lsquit('PLT DRIVER not implemented for unrestricted calculations!',-1)
    end if

    ! Sanity check 2: Density matrix algebra
    if(matrix_type/=mtype_dense) then
       call lsquit('PLT DRIVER only implemented for dense matrices!',-1)
    end if


    ! read density/orbitals
    ! *********************
    write(ls%lupri,*) 'PLT driver reads input file...'
    funit = 33
    OPEN(UNIT=funit,FILE=trim(MyPlt%inputfile),STATUS='OLD', &
         & FORM='UNFORMATTED',IOSTAT=IOS)
    READ (funit) nrow,ncol
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
    call construct_plt_file_driver(trim(MyPlt%outputfile),trim(MyPlt%frmt),InputMat,&
         & ls,MyPlt%iorb,MyPlt%jorb)


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


  !> \brief Driver for creating PLT files of orbitals, densities, electrostatic potentials etc.
  !> which may be plotted e.g. using the Chimera program.
  !> \author Kasper Kristensen
  !> \date May 2013
  subroutine construct_plt_file_driver(filename,frmt,InputMat,ls,iorb,jorb)
    implicit none
    !> Matrix containing orbitals or density
    type(matrix),intent(in)        :: InputMat
    !> Integral info
    type(lsitem),intent(inout)     :: ls
    !> Filename for output PLT file
    Character*(*),intent(in)       :: filename
    !> What do we want to calculate (see type pltinfo)
    Character*(*),intent(in)       :: frmt
    !> Orbital index to plot or first index in charge distribution (see type pltinfo)
    integer,intent(in)              :: iorb
    !> Second index in charge distribution (see type pltinfo)
    integer,intent(in)              :: jorb
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


    ! Which calculation?
    ! ******************

    select case (frmt)

       ! Density
    case('DENS')
       write(ls%lupri,*) 'Writing density distribution plt file...'
       call calculate_density2(trim(filename),InputMat,ls,natoms,ATOMXYZ)

       ! Electrostatic potential
    case('EP')
       write(ls%lupri,*) 'Writing electrostatic potential plt file...'
       call calculate_ep3(trim(filename),InputMat,ls,natoms,ATOMXYZ)

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
       call calculate_pplt(trim(filename),iorb,InputMat,ls,natoms,ATOMXYZ)

       ! Charge distribution
    case('CHARGEDIST')
       write(6,*) 'Writing charge distribution plt file...'
       call calculate_charge(trim(filename),iorb,jorb,InputMat,ls,natoms,ATOMXYZ)

    case default
       write(ls%lupri,*) 'PLT driver unknown input format: ', frmt
       call lsquit('PLT driver unknown input format',-1)
    end select

    deallocate(ATOMXYZ)
  end subroutine construct_plt_file_driver


  subroutine print_xyz(filename,ls)
    implicit none
    type(lsitem)                   :: ls
    Character*(*)                  :: filename
    integer                        :: nATOMS
    integer                        :: I,J, IUNIT, ANUM
    real(4)                        :: ATOMXYZ(3)
    real(4)                        :: br
    br=real(bohr_to_angstrom,4)
    IUNIT=111
    nATOMS = ls%setting%MOLECULE(1)%p%natoms

    OPEN (IUNIT,FILE=trim(filename),STATUS='REPLACE',FORM='FORMATTED')

    WRITE(IUNIT,'(I5)') nATOMS

    WRITE(IUNIT,*)

    DO I=1,nATOMS !Atom coordinates        
       ANUM=NINT(ls%setting%molecule(1)%p%atom(I)%Charge)
       DO J=1,3                                                     
          ATOMXYZ(J) = ls%setting%molecule(1)%p%atom(I)%CENTER(J)
       ENDDO

       WRITE(IUNIT,'(I3,F10.5,F10.5,F10.5)') ANUM, (ATOMXYZ(J)*br,J=1,3)

    ENDDO

    CLOSE(IUNIT,STATUS='KEEP')

  end subroutine print_xyz

  subroutine calculate_daltongrid(filename,nORBITALS,dCMO,ls,natoms,ATOMXYZ)
    implicit none
    type(matrix)                   :: dCMO
    type(lsitem)                   :: ls
    Character*(*)                  :: filename
    integer,intent(in)             :: natoms
    integer                        :: nX,nY,nZ
    real(4), intent(in)            :: ATOMXYZ(3,nATOMS)
    !Local
    integer                        :: I,J,P,Q,Xg,Yg,Zg,nGRIDPOINTS
    integer                        :: orbnr, nORBITALS,nBASIS,iorb,iunit,ig
    real(4)                        :: X,Y,Z,X1,Y1,Z1,Xgrid, Ygrid,Zgrid,deltax, deltay, deltaz 
    real(4), allocatable           :: GAO(:,:),CMO(:)
    real(4), allocatable           :: moorb(:,:)
    real(4)                        :: fZ1,fY1,fX1,fZn,fYn,fXn
    integer                        :: II,TD, ijk,irec, nLeftover
    integer                        :: reclen
    real(4)                        :: br
    integer, parameter             :: nGRS = 2048
    br=real(bohr_to_angstrom,4)

    nBASIS    = dCMO%nrow
    ! Determine grid-box. It extends from
    ! X1/Y1/Z1 to Xn/Yn/Zn in the x-/y-/z-directions,
    ! and the number of gridpoints in these directions are nX/nY/nZ.
    ! Thus, the total number of gridpoints is nGRIDPOINT=nX*nY*nZ.
    call DETERMINE_GRIDBOX(X1,nX,Y1,nY,Z1,nZ,deltax,deltay,deltaz,&
         &ATOMXYZ,natoms,nGRIDPOINTS)


    !Open daltongrid file and write header
    fZ1=Z1*br; fZn = fZ1 + br*deltaz*(nZ-1);
    fY1=Y1*br; fYn = fY1 + br*deltay*(nY-1);
    fX1=X1*br; fXn = fX1 + br*deltax*(nX-1);

    IUNIT=111

    reclen=nGRS

#ifdef GFORTRAN
    reclen=reclen*4
#endif

    OPEN (IUNIT,FILE=trim(filename),STATUS='UNKNOWN',FORM='UNFORMATTED',&
         & ACCESS='DIRECT',RECL=reclen)

    irec=1
    WRITE(IUNIT,REC=irec) nZ,nY,nX,fZ1,fZn,fY1,fYn,fX1,fXn,nORBITALS,nGRIDPOINTS,nGRS

    !allocate
    allocate(moorb(nGRS,nORBITALS), GAO(nBASIS,nGRS),CMO(nORBITALS*nBASIS)) !MO orbitals

    moorb = 0e0; ijk = 0
    CMO = real(dCMO%elms,4)

    DO Zg = 1,nZ 
       Zgrid = Z1 + deltaz*(Zg-1) 
       DO Yg = 1,nY 
          Ygrid = Y1 + deltay*(Yg-1) 
          DO Xg = 1,nX
             Xgrid = X1 + deltax*(Xg-1) 



             orbnr=1
             ijk = ijk + 1

             GAO(:,ijk) = 0.0_4

             DO I=1,natoms
                ! X-, Y-, and Z-distances from gridpoint to atom I
                X = Xgrid - ATOMXYZ(1,I)
                Y = Ygrid - ATOMXYZ(2,I)
                Z = Zgrid - ATOMXYZ(3,I)
                call determine_orbitals(I,GAO(:,ijk),orbnr, nBASIS, X,Y,Z,ls)
             ENDDO



             !Compute and write nGRD grid points             
             IF (mod(ijk,nGRS).eq. 0) THEN

                call SGEMM('t','n',nGRS,nORBITALS,nBASIS,1.0,GAO,nBASIS,CMO,nBASIS,0.0,moorb,nGRS)

                DO iorb=1,nORBITALS
                   irec = irec + 1
                   WRITE(IUNIT,REC=irec) moorb(:,iorb)
                ENDDO

                ijk = 0
             ENDIF

          ENDDO
       ENDDO
    ENDDO


    !Write leftover gridpoints
    nLeftover = mod(nGRIDPOINTS,nGRS)
    IF (nLeftover.ne. 0) THEN

       call SGEMM('t','n',nGRS,nORBITALS,nBASIS,1.0,GAO,nBASIS,CMO,nBASIS,0.0,moorb,nGRS)

       DO iorb=1,nORBITALS
          irec = irec + 1
          WRITE(IUNIT,REC=irec) moorb(:,iorb)
       ENDDO
    ENDIF


    CLOSE(IUNIT,STATUS='KEEP')

    deallocate(moorb,GAO,CMO)


  end subroutine calculate_daltongrid



  subroutine calculate_charge(filename,iorb,jorb,dCMO,ls,natoms,ATOMXYZ)
    implicit none
    type(matrix)                   :: dCMO
    type(lsitem)                   :: ls
    Character*(*)                  :: filename
    integer,intent(in)             :: natoms
    integer                        :: nX,nY,nZ
    real(4), intent(in)            :: ATOMXYZ(3,nATOMS)
    !Local
    integer                        :: I,J,P,Q,Xg,Yg,Zg,nGRIDPOINTS
    integer                        :: orbnr, nORBITALS,nBASIS,iorb,jorb,iunit,ig
    real(4)                        :: X,Y,Z,X1,Y1,Z1,Xgrid, Ygrid,Zgrid,deltax, deltay, deltaz 
    real(4), allocatable           :: moorb(:),CMO(:,:),GAO(:),tmp(:)
    real(4)                        :: fZ1,fY1,fX1,fZn,fYn,fXn, SDOT
    integer                        :: II,TD, ijk
    integer                        :: reclen
    real(4)                       :: br
    br=real(bohr_to_angstrom,4)

    nORBITALS = dCMO%ncol
    nBASIS    = dCMO%nrow
    ! Determine grid-box. It extends from
    ! X1/Y1/Z1 to Xn/Yn/Zn in the x-/y-/z-directions,
    ! and the number of gridpoints in these directions are nX/nY/nZ.
    ! Thus, the total number of gridpoints is nGRIDPOINT=nX*nY*nZ.
    call DETERMINE_GRIDBOX(X1,nX,Y1,nY,Z1,nZ,deltax,deltay,deltaz,&
         &ATOMXYZ,natoms,nGRIDPOINTS)

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

             moorb(ijk) = SDOT(nBASIS,CMO(:,iorb),1,GAO,1) * SDOT(nBASIS,CMO(:,jorb),1,GAO,1)


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
#ifdef GFORTRAN
    reclen=reclen*4
#endif

    OPEN (IUNIT,FILE=trim(filename),STATUS='REPLACE',FORM='UNFORMATTED',&
         & ACCESS='DIRECT',RECL=reclen)

    WRITE(IUNIT,REC=1) int(II,4),int(TD,4),int(nZ,4),int(nY,4),int(nX,4),fZ1,fZn,fY1,fYn,fX1,fXn,moorb

    CLOSE(IUNIT,STATUS='KEEP')

    deallocate(moorb)


  end subroutine calculate_charge


  subroutine calculate_pplt(filename,iorb,dCMO,ls,natoms,ATOMXYZ)
    implicit none
    type(matrix)                   :: dCMO
    type(lsitem)                   :: ls
    Character*(*)                  :: filename
    integer,intent(in)             :: natoms
    integer                        :: nX,nY,nZ
    real(4), intent(in)            :: ATOMXYZ(3,nATOMS)
    !Local
    integer                        :: I,J,P,Q,Xg,Yg,Zg,nGRIDPOINTS
    integer                        :: orbnr, nORBITALS,nBASIS,iorb,iunit,ig
    real(4)                        :: X,Y,Z,X1,Y1,Z1,Xgrid, Ygrid,Zgrid,deltax, deltay, deltaz 
    real(4),     pointer           :: GAO(:)
    real(4), allocatable           :: moorb(:),CMO(:)
    real(4)                        :: fZ1,fY1,fX1,fZn,fYn,fXn,SDOT
    integer                        :: II,TD, ijk
    integer                        :: reclen
    real(4)                       :: br
    br=real(bohr_to_angstrom,4)

    nORBITALS = dCMO%ncol
    nBASIS    = dCMO%nrow
    ! Determine grid-box. It extends from
    ! X1/Y1/Z1 to Xn/Yn/Zn in the x-/y-/z-directions,
    ! and the number of gridpoints in these directions are nX/nY/nZ.
    ! Thus, the total number of gridpoints is nGRIDPOINT=nX*nY*nZ.
    call DETERMINE_GRIDBOX(X1,nX,Y1,nY,Z1,nZ,deltax,deltay,deltaz,&
         &ATOMXYZ,natoms,nGRIDPOINTS)

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

             moorb(ijk) = SDOT(nBASIS,CMO,1,GAO,1)


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

#if 0
  ! subroutine calculate_ep now obsolete
  subroutine calculate_ep(filename,dD,ls,natoms,ATOMXYZ)
    implicit none
    type(matrix)                   :: dD
    type(matrix)                   :: dEPint
    type(lsitem)                   :: ls
    TYPE(LSSETTING)                :: psetting
    Character*(*)                  :: filename
    integer,intent(in)             :: natoms
    integer                        :: nX,nY,nZ
    real(4), intent(in)            :: ATOMXYZ(3,nATOMS)
    !Local
    real(realk)                    :: R(3)
    integer                        :: I,Q,Xg,Yg,Zg,nGRIDPOINTS
    integer                        :: orbnr, nORBITALS,nBASIS,iunit,ig
    real(4)                        :: X,Y,Z,X1,Y1,Z1,Xgrid, Ygrid,Zgrid,deltax, deltay, deltaz 
    real(4), allocatable           :: moorb(:),D(:),EPint(:)
    real(4)                        :: fZ1,fY1,fX1,fZn,fYn,fXn,epnuc,SDOT
    integer                        :: II,TD, ijk, icount
    integer                        :: reclen
    real(4)                       :: br
    br=real(bohr_to_angstrom,4)

    nORBITALS = dD%ncol
    nBASIS    = dD%nrow
    ! Determine grid-box. It extends from
    ! X1/Y1/Z1 to Xn/Yn/Zn in the x-/y-/z-directions,
    ! and the number of gridpoints in these directions are nX/nY/nZ.
    ! Thus, the total number of gridpoints is nGRIDPOINT=nX*nY*nZ.
    call DETERMINE_GRIDBOX(X1,nX,Y1,nY,Z1,nZ,deltax,deltay,deltaz,&
         &ATOMXYZ,natoms,nGRIDPOINTS)

    allocate(moorb(nGRIDPOINTS),D(nORBITALS*nBASIS)) !MO orbitals
    moorb = 0_4; ijk = 0

    D=real(dD%elms,4)


    call mem_TurnONThread_Memory()
    ls%setting%scheme%noOMP = .TRUE. 


    write(6,*) 'Computing ', nGRIDPOINTS, ' gridpoints'


    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(psetting,Zg,Zgrid,Yg,Ygrid,Xg,Xgrid,epnuc,I,X,Y,Z,Q,R,dEPint,EPint,ijk)

    call init_threadmemvar()

    call mat_init(dEPint,nBASIS,nBASIS)
    allocate(EPint(nBASIS*nBASIS))

    !$OMP CRITICAL
    call copy_setting(psetting,ls%setting,ls%lupri)
    !$OMP END CRITICAL


    !$OMP DO SCHEDULE(STATIC)
    DO Zg = 1,nZ 
       Zgrid = Z1 + deltaz*(Zg-1) 
       DO Yg = 1,nY 
          Ygrid = Y1 + deltay*(Yg-1) 
          DO Xg = 1,nX
             Xgrid = X1 + deltax*(Xg-1) 


             epnuc = 0_4

             DO I=1,natoms
                ! X-, Y-, and Z-distances from gridpoint to atom I
                X = Xgrid - ATOMXYZ(1,I)
                Y = Ygrid - ATOMXYZ(2,I)
                Z = Zgrid - ATOMXYZ(3,I)
                Q = real(ls%setting%molecule(1)%p%atom(I)%Charge,4)
                epnuc = epnuc + (Q/SQRT(X*X + Y*Y + Z*Z))
             ENDDO

             R = real((/Xgrid, Ygrid, Zgrid/),realk)
             call II_get_ep_integrals(ls%lupri,ls%luerr,psetting,dEPint,R) 
             EPint=real(dEPint%elms,4)

             ijk = Xg + (Yg-1)*nX + (Zg-1)*nX*nY

             moorb(ijk) = epnuc - 2.0*SDOT(nBASIS*nBASIS,D,1,EPint,1)

          ENDDO
       ENDDO
    ENDDO
    !$OMP END DO

    !$OMP CRITICAL
    call typedef_free_setting(psetting)
    !$OMP END CRITICAL

    deallocate(EPint)
    call mat_free(dEPint)

    call collect_thread_memory()

    !$OMP END PARALLEL

    call mem_TurnOffThread_Memory()
    ls%setting%scheme%noOMP = .FALSE. 


    !Write plt file

    II=3; TD=200
    fZ1=Z1*br; fZn = fZ1 + br*deltaz*(nZ-1);
    fY1=Y1*br; fYn = fY1 + br*deltay*(nY-1);
    fX1=X1*br; fXn = fX1 + br*deltax*(nX-1);

    IUNIT=111
    reclen=(nGRIDPOINTS + 11)
#ifdef GFORTRAN
    reclen=reclen*4
#endif

    OPEN (IUNIT,FILE=trim(filename),STATUS='REPLACE',FORM='UNFORMATTED',&
         & ACCESS='DIRECT',RECL=reclen)

    WRITE(IUNIT,REC=1) int(II,4),int(TD,4),int(nZ,4),int(nY,4),int(nX,4),fZ1,fZn,fY1,fYn,fX1,fXn,moorb

    CLOSE(IUNIT,STATUS='KEEP')

    deallocate(moorb,D)


  end subroutine calculate_ep
#endif

  subroutine calculate_ep3(filename,D,ls,natoms,ATOMXYZ)
    implicit none
    type(matrix)                   :: D
    type(matrix)                   :: dEPint
    type(lsitem)                   :: ls
    TYPE(LSSETTING)                :: psetting
    Character*(*)                  :: filename
    integer,intent(in)             :: natoms
    integer                        :: nX,nY,nZ
    real(4), intent(in)            :: ATOMXYZ(3,nATOMS)
    !Local
    real(realk), allocatable       :: R(:,:), emoorb(:)
    integer                        :: I,Q,Xg,Yg,Zg,nGRIDPOINTS
    integer                        :: orbnr, nORBITALS,nBASIS,iunit,ig
    real(4)                        :: X,Y,Z,X1,Y1,Z1,Xgrid, Ygrid,Zgrid,deltax, deltay, deltaz 
    real(4), allocatable           :: moorb(:)
    real(4)                        :: fZ1,fY1,fX1,fZn,fYn,fXn,epnuc
    integer                        :: II,TD, ijk, icount
    integer                        :: reclen
    real(4)                       :: br
    br=real(bohr_to_angstrom,4)

    nORBITALS = D%ncol
    nBASIS    = D%nrow
    ! Determine grid-box. It extends from
    ! X1/Y1/Z1 to Xn/Yn/Zn in the x-/y-/z-directions,
    ! and the number of gridpoints in these directions are nX/nY/nZ.
    ! Thus, the total number of gridpoints is nGRIDPOINT=nX*nY*nZ.
    call DETERMINE_GRIDBOX(X1,nX,Y1,nY,Z1,nZ,deltax,deltay,deltaz,&
         &ATOMXYZ,natoms,nGRIDPOINTS)

    allocate(moorb(nGRIDPOINTS),emoorb(nGRIDPOINTS),R(3,nGRIDPOINTS)) !MO orbitals
    moorb = 0_4; emoorb = 0_realk;

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

             epnuc = 0_4

             DO I=1,natoms
                ! X-, Y-, and Z-distances from gridpoint to atom I
                X = Xgrid - ATOMXYZ(1,I)
                Y = Ygrid - ATOMXYZ(2,I)
                Z = Zgrid - ATOMXYZ(3,I)
                Q = real(ls%setting%molecule(1)%p%atom(I)%Charge,4)
                epnuc = epnuc + (Q/SQRT(X*X + Y*Y + Z*Z))
             ENDDO

             R(:,ijk) = real((/Xgrid, Ygrid, Zgrid/),realk)

             moorb(ijk) = epnuc 

          ENDDO
       ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    !   Electronic contribution

    call II_get_ep_integrals3(ls%lupri,ls%luerr,ls%setting,R,nGRIDPOINTS,D,emoorb)

    !   Add the two

    deallocate(R)

    moorb = moorb - real(emoorb,4)

    deallocate(emoorb)


    !Write plt file

    II=3; TD=200
    fZ1=Z1*br; fZn = fZ1 + br*deltaz*(nZ-1);
    fY1=Y1*br; fYn = fY1 + br*deltay*(nY-1);
    fX1=X1*br; fXn = fX1 + br*deltax*(nX-1);

    IUNIT=111
    reclen=(nGRIDPOINTS + 11)
#ifdef GFORTRAN
    reclen=reclen*4
#endif

    OPEN (IUNIT,FILE=trim(filename),STATUS='REPLACE',FORM='UNFORMATTED',&
         & ACCESS='DIRECT',RECL=reclen)


    WRITE(IUNIT,REC=1) int(II,4),int(TD,4),int(nZ,4),int(nY,4),int(nX,4),fZ1,fZn,fY1,fYn,fX1,fXn,moorb

    CLOSE(IUNIT,STATUS='KEEP')

    deallocate(moorb)


  end subroutine calculate_ep3



  subroutine calculate_density2(filename,dD,ls,natoms,ATOMXYZ)
    implicit none
    type(matrix)                   :: dD
    type(lsitem)                   :: ls
    Character*(*)                  :: filename
    integer,intent(in)             :: natoms
    integer                        :: nX,nY,nZ
    real(4), intent(in)            :: ATOMXYZ(3,nATOMS)
    !Local
    integer                        :: I,J,P,Q,Xg,Yg,Zg,nGRIDPOINTS
    integer                        :: orbnr,nBASIS,iorb,iunit,ig
    real(4)                        :: X,Y,Z,X1,Y1,Z1,Xgrid, Ygrid,Zgrid,deltax, deltay, deltaz 
    real(4), allocatable           :: GAO(:,:),D(:),tmp(:,:)
    real(4), allocatable           :: moorb(:),XgridB(:),YgridB(:),ZgridB(:)
    real(4)                        :: fZ1,fY1,fX1,fZn,fYn,fXn,SDOT
    integer                        :: II,TD, ijk, klm, nLeftover
    integer                        :: reclen
    real(4)                        :: br
    br=real(bohr_to_angstrom,4)

    nBASIS    = dD%nrow
    ! Determine grid-box. It extends from
    ! X1/Y1/Z1 to Xn/Yn/Zn in the x-/y-/z-directions,
    ! and the number of gridpoints in these directions are nX/nY/nZ.
    ! Thus, the total number of gridpoints is nGRIDPOINT=nX*nY*nZ.
    call DETERMINE_GRIDBOX(X1,nX,Y1,nY,Z1,nZ,deltax,deltay,deltaz,&
         &ATOMXYZ,natoms,nGRIDPOINTS)

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
                   moorb(ijk-nBASIS+klm) = 2.0*SDOT(nBASIS,GAO(:,klm),1,tmp(:,klm),1)
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
          moorb(ijk-nLeftover+klm) = 2.0*SDOT(nBASIS,GAO(:,klm),1,tmp(:,klm),1)
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
    !removed the #ifdef VAR_GFORTRAN as this is also needed for ifort. 
    reclen=reclen*4

    OPEN (IUNIT,FILE=trim(filename),STATUS='REPLACE',FORM='UNFORMATTED',&
         & ACCESS='DIRECT',RECL=reclen)

    WRITE(IUNIT,REC=1) int(II,4),int(TD,4),int(nZ,4),int(nY,4),int(nX,4),fZ1,fZn,fY1,fYn,fX1,fXn,moorb

    CLOSE(IUNIT,STATUS='KEEP')

    deallocate(moorb)


  end subroutine calculate_density2


  subroutine make_orbitalplot_file(CMO,CFG,ls)
    implicit none
    !> file containing localized orbitals
    type(matrix) :: CMO 
    type(RedSpaceItem) :: CFG
    type(lsitem) :: ls
    integer :: indx1,indx2,nocc


    IF(ls%setting%integraltransformGC)THEN
       call transform_basis_to_GCAO(ls%input%basis)
       ls%setting%integraltransformGC = .FALSE.
    ENDIF

    nocc=int(ls%input%molecule%nelectrons/2)

    select case (trim(CFG%plt_orbital))
    case('HOMO')
       indx1 = nocc
       write(ls%lupri,*) 'Writing plt file for HOMO...'
       call construct_plt_file_driver('HOMO.plt','ORB',CMO,ls,indx1,-1)
    case('LUMO')
       indx1 = nocc+1
       write(ls%lupri,*) 'Writing plt file for LUMO...'
       call construct_plt_file_driver('LUMO.plt','ORB',CMO,ls,indx1,-1)
    case('LEASTL')
       write(ls%lupri,*) 'Writing plt file for the least local occupied and virtual orbitals...'
       indx1=CFG%leastl_occ
       indx2=CFG%leastl_virt
       if (indx1>0) call construct_plt_file_driver('LeastLocal_occ.plt','ORB',CMO,ls,indx1,-1)
       if (indx2>0) call construct_plt_file_driver('Leastlocal_virt.plt','ORB',CMO,ls,indx2,-1)
    case('MOSTL')
       write(ls%lupri,*) 'Writing plt file for the most local occupied and virtual orbitals...'
       indx1=CFG%mostl_occ
       indx2=CFG%mostl_virt
       if (indx1>0) call construct_plt_file_driver('MostLocal_occ.plt','ORB',CMO,ls,indx1,-1)
       if (indx2>0) call construct_plt_file_driver('Mostlocal_virt.plt','ORB',CMO,ls,indx2,-1)
    case('ALL')
       write(ls%lupri,*) 'Writing plt file for the least local occupied and virtual orbitals...'
       indx1=CFG%leastl_occ
       indx2=CFG%leastl_virt
       if (indx1>0) call construct_plt_file_driver('LeastLocal_occ.plt','ORB',CMO,ls,indx1,-1)
       if (indx2>0) call construct_plt_file_driver('Leastlocal_virt.plt','ORB',CMO,ls,indx2,-1)
       write(ls%lupri,*) 'Writing plt file for the most local occupied and virtual orbitals...'
       indx1=CFG%mostl_occ
       indx2=CFG%mostl_virt
       if (indx1>0) call construct_plt_file_driver('MostLocal_occ.plt','ORB',CMO,ls,indx1,-1)
       if (indx2>0) call construct_plt_file_driver('Mostlocal_virt.plt','ORB',CMO,ls,indx2,-1)
    case default
       write(ls%lupri,*) '========== SOMETHING WRONG WHEN MAKING .plt FILE ================'
       write(ls%lupri,*) 'None of the valid options chosen. Use the options'
       write(ls%lupri,*) ' HOMO,LUMO,MOSTL,LEASTL or ALL below keyword .ORBITAL PLOT'
       write(ls%lupri,*) ' If you want to plot other than the above orbitals, consider using'
       write(ls%lupri,*) ' **PLOT option where you may specify any orbital index.'
       write(ls%lupri,*) '================================================================='

       write(*,*) '========== SOMETHING WRONG WHEN MAKING .plt FILE ================'
       write(*,*) 'None of the valid options chosen. Use the options'
       write(*,*) ' HOMO,LUMO,MOSTL,LEASTL or ALL below keyword .ORBITAL PLOT'
       write(*,*) ' If you want to plot other than the above orbitals, consider using'
       write(*,*) ' **PLOT option where you may specify any orbital index.'
       write(*,*) '================================================================='
    end select



  end subroutine make_orbitalplot_file


end module plt_driver_module
