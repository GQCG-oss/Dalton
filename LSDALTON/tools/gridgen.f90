program gridgen
  use print_moorb_grid_mod
  use init_lsdalton_mod
  use IntegralInterfaceMOD,only:II_Get_overlap
  !use typedef
  !use matrix_operations
  integer :: nocc, nrow,ncol,lun, IOS, m(2), iorb, jorb
  Type(Matrix) :: CMO,S,tmpMat,tmpMat2
  type(LsItem) :: ls
  character*(80) :: inputfile
  character*(80) :: outputfile
  character*(80) :: frmt
  character*(80) :: arg

  matrix_type=mtype_dense

  call getarg(1,inputfile)
  call getarg(2,outputfile)
  call getarg(3,frmt)


  write(6,*) 'Getting lsitem from input files DALTON.INP and MOLECULE.INP'

  call init_globalmemvar  !initialize the global memory counters
  call set_matrix_default !initialize global matrix counters
  call lstmem_init
  call MatrixmemBuf_init()
  call init_timers !initialize timers
  ! MPI initialization
  call lsmpi_init
  call get_lsitem_from_input(ls)


  !    write(6,*) 'Reading lsitem...'
  !    call read_lsitem_from_disk(ls)

  ls%lupri=6
  ls%luerr=6

  nocc=int(ls%input%molecule%nelectrons/2)
  write(6,*) 'nocc= ', nocc

  ! read orbitals
  write(6,*) 'Reading orbitals...'
  LUN=33
  OPEN(UNIT=LUN,FILE=trim(inputfile),STATUS='OLD', &
       &     FORM='UNFORMATTED',IOSTAT=IOS)
  READ (LUN) nrow,ncol

  write(6,*) 'Orbital matrix size: ', nrow, ncol
  call mat_init(Cmo,nrow,ncol)

  READ (LUN) CMO%elms
  CLOSE(LUN,STATUS='KEEP',IOSTAT=IOS)
  write(6,*) 'Done'

  IF(ls%setting%integraltransformGC)THEN
     call transform_basis_to_GCAO(ls%input%basis)
     ls%setting%integraltransformGC = .FALSE.
  ENDIF


  select case (trim(frmt))
  case('CUBE')
     write(6,*) 'Writing cube file...'
     call print_moorb(trim(outputfile),'CUBE',CMO,ls,-1,-1)
  case('DPLT')
     write(6,*) 'Writing dplt file...'
     call getarg(4,arg)
     read(arg,*) iorb
     call print_moorb(trim(outputfile),'DPLT',CMO,ls,iorb,-1)
  case('DENS')
     write(6,*) 'Writing density distribution plt file...'
     call print_moorb(trim(outputfile),'DENS',CMO,ls,-1,-1)
  case('EP')
     write(6,*) 'Writing electrostatic potential plt file...'
     call print_moorb(trim(outputfile),'EP',CMO,ls,-1,-1)
  case('ORB')

     call mat_init(S,nrow,ncol)
     CALL II_get_overlap(ls%lupri,ls%luerr,ls%setting,S)
     call mat_init(tmpMat,nrow,ncol)
     call mat_init(tmpMat2,nrow,ncol)
     call mat_mul(CMO,S,'t','n',1E0_realk,0E0_realk,tmpMat)
     call mat_mul(tmpMat,CMO,'n','n',1E0_realk,0E0_realk,S)
     !S should now be the identity
     call mat_identity(tmpMat)
     call mat_add(1E0_realk,S,-1E0_realk,tmpMat,tmpMat2)
     print*,'mat_sqnorm(C^T*S*C - I) =',mat_sqnorm2(tmpMat2)/tmpMat2%nrow
     IF(ABS(mat_sqnorm2(tmpMat2)/tmpMat2%nrow).GT.1.0E-15_realk)THEN
        stop 'gridgen error C^T*S*C = I Error.'       
     ENDIF
     call mat_free(S)
     call mat_free(tmpMat)
     call mat_free(tmpMat2)

     write(6,*) 'Writing plt file...'
     call getarg(4,arg)
     read(arg,*) iorb
     call print_moorb(trim(outputfile),'ORB',CMO,ls,iorb,-1)
  case('CHARGEDIST')
     write(6,*) 'Writing charge distribution plt file...'
     call getarg(4,arg)
     read(arg,*) iorb
     call getarg(5,arg)
     read(arg,*) jorb
     call print_moorb(trim(outputfile),'CHARGEDIST',CMO,ls,iorb,jorb)
  case default
     print *, 'Input: ', trim(frmt)
     stop 'Gridgen: Input keyword not recognized!'
  end select


  write(6,*) 'Done'

  call mat_free(Cmo)

  call MatrixmemBuf_free()
  call lstmem_free
  stop 'gridgen done'

endprogram gridgen
