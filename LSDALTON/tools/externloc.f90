program externloc
!use print_moorb_grid_mod
use typedefType, only: lsitem
use matrix_module
use matrix_operations
use davidson_settings
use kurtosis
use init_lsdalton_mod, only: init_lsdalton_and_get_lsitem, get_lsitem_from_input
use memory_handling
use lstiming
use lstensorMem, only: lstmem_init, lstmem_free
integer :: nocc, nrow,ncol,lun, IOS, m(2)
Type(Matrix) :: CMO,cmoocc
integer :: nb, nr
real(realk),allocatable  :: cmomat(:,:)
type(LsItem) :: ls
real(realk)  :: mx(2)
type(RedSpaceItem) :: CFG
character*(80) :: inputfile
character*(80) :: outputfile
character*(80) :: arg
character*(80) :: loctype
character*(80) :: loctype1


 call davidson_default_OrbLoc(CFG)

 matrix_type=mtype_dense

 call getarg(1,inputfile)
 call getarg(2,outputfile)
 call getarg(3,arg)
 read(arg,*) m(1)
 call getarg(4,arg)
 read(arg,*) m(2)
 call getarg(5,loctype1)


    loctype  = trim(loctype1)

 
    if (loctype .eq. 'PRINTINFO') then
        CFG%PRINT_INFO=.true. 
    elseif (loctype .eq. 'CLL') then
        CFG%PM_input%ChargeLocLowdin=.true.
        CFG%PM=.true.
    elseif (loctype .eq. 'KURT' ) then
        CFG%PFM=.true.
	CFG%precond=.true.
        CFG%PFM_input%crossterms=.false.
	if (m(1) .ne. 0) then
  	  CFG%PFM_input%m=m(1)
	elseif(m(2).ne. 0) then
	  CFG%PFM_input%m=m(2)
	end if
	CFG%PM=.false.
	CFG%orbspread=.false.
    elseif (loctype .eq. 'PML') then
        CFG%PM_input%PipekMezeyLowdin=.true.
        CFG%PM=.true.   
	CFG%PM_input%Precond = .true.
	CFG%Precond= .true.
    elseif (loctype .eq. 'PMM') then
        CFG%PM_input%PipekMezeyMull=.true.
        CFG%PM=.true.   
	CFG%PM_input%Precond = .false.
	CFG%Precond= .false.
    elseif (loctype .eq. 'CLM') then
        CFG%PM_input%ChargeLocMulliken=.true.
        CFG%PM=.true.   
    elseif (loctype .eq. 'MLO') then
        CFG%orbspread =.true.
    end if


    write(6,*) 'Getting lsitem from input files DALTON.INP and MOLECULE.INP'
    call init_globalmemvar !initialize the global memory counters
    call set_matrix_default !initialize global matrix counters
    call lstmem_init
    call MatrixmemBuf_init()
    call init_timers 
    call get_lsitem_from_input(ls)
!    write(6,*) 'Reading lsitem...'
!    call read_lsitem_from_disk(ls)

    ls%lupri=6
  
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

    write(6,*)
    write(6,*)  '  LOCALIZATION INFO'
    write(6,*)  '*********************************'
    write(6,*)  'Localization:   ',loctype
    write(6,*)  'Exponent, occ:  ', m(1)
    write(6,*)  'Exponent, virt: ', m(2)
    write(6,*)  'Orbital file:   ', trim(inputfile) 
    write(6,*)  '**********************************'
    write(6,*)
! nb=cmo%nrow
! nr=5
! allocate(cmomat(nb,nr))
! call mat_retrieve_block(cmo,cmomat,nb,nr,1,1) 
! call mat_init(cmoocc,nb,nr)
! call mat_set_from_full(cmomat,1.0d0,cmoocc)
! call kurtosis_driver(ls,cmoocc,nb,nr)
!
    if ((m(1).ne.0).or.(m(2).ne.0)) then
       call optimloc(CMO,nocc,m,ls,CFG)

       write(6,*) 'Writing localized orbitals...'
       OPEN(UNIT=LUN,FILE=trim(outputfile),STATUS='REPLACE', &
      &     FORM='UNFORMATTED',IOSTAT=IOS)

       call mat_dense_write_to_disk(lun,CMO)
       CLOSE(LUN,STATUS='KEEP',IOSTAT=IOS)
       write(6,*) 'Done'
    endif
    call lstmem_free
    call MatrixmemBuf_free()

    stop 'externloc done'

endprogram
