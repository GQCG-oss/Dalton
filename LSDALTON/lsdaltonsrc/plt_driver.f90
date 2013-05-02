!> @file
!> Driver for creating PLT files of orbitals, densities, electrostatic potentials etc.
!> which may be plotted e.g. using the Chimera program.
!> \author Kasper Kristensen (based on program by Branislav Jansik)

module plt_driver_module

contains

!> Driver for creating PLT files of orbitals, densities, electrostatic potentials etc.
!> which may be plotted e.g. using the Chimera program.
  subroutine contruct_plt_file_driver(ls,inputfile,outputfile,frmt,iorb,jorb)

    use print_moorb_grid_mod
    use matrix_operations
    use init_lsdalton_mod
    use IntegralInterfaceMOD,only:II_Get_overlap

    !> Integral info
    type(LsItem) :: ls
    !> Input file containing orbitals (e.g. lcm_orbitals.u) or density (e.g. dens.restart)
    character*(80) :: inputfile
    !> PLT output file to which information is written
    character*(80) :: outputfile
    !> What to calculate:
    !> frmt = 'DENS': Calculate electron density at grid points (inputfile=density matrix)
    !> frmt = 'EP'  : Calculate electrostatic potential at grid points (inputfile=density matrix)
    !> frmt = 'ORB' : Calculate specific molecular orbitals at grid points (inputfile=orbital matrix)
    !> frmt = 'CHARGEDIST': Calculate charge distribution between two orbitals (input file=orbital mat)
    character*(80) :: frmt
    !> Index for which orbital to plot (only used for ORB and CHARGEDIST)
    integer,intent(in) :: iorb
    !> Index for second orbital in charge distribution (only used for CHARGEDIST)
    integer,intent(in) :: jorb
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
    OPEN(UNIT=funit,FILE=trim(inputfile),STATUS='OLD', &
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

    select case (trim(frmt))

       ! Density
    case('DENS')
       write(ls%lupri,*) 'Writing density distribution plt file...'
       call calculate_density2(trim(outputfile),InputMat,ls,natoms,ATOMXYZ)

       ! Electrostatic potential
    case('EP')
       write(ls%lupri,*) 'Writing electrostatic potential plt file...'
       call calculate_ep3(trim(outputfile),InputMat,ls,natoms,ATOMXYZ)

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
       call calculate_pplt(trim(outputfile),iorb,InputMat,ls,natoms,ATOMXYZ)

       ! Charge distribution
    case('CHPLT')
       write(6,*) 'Writing charge distribution plt file...'
       call calculate_charge(trim(outputfile),iorb,jorb,InputMat,ls,natoms,ATOMXYZ)

    case default
       write(ls%lupri,*) 'PLT driver unknown input format: ', frmt
       call lsquit('PLT driver unknown input format',-1)
    end select


    call mat_free(InputMat)
    deallocate(ATOMXYZ)

  end subroutine contruct_plt_file_driver


end module plt_driver_module
