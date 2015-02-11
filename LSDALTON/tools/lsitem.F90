program make_lsitem
  use configuration
  use daltoninfo
  use IIDFTINT
  use IIDFTD, only: II_DFT_DISP

implicit none
  TYPE(lsitem),target :: ls
  integer             :: nbast !,lu_error,lu_output
!  CHARACTER(len=9)    :: BASISLABEL
  LOGICAL             :: isdft,restart_from_dens,isgga, do_decomp!,isminimum
  integer             :: matmul1, matmul2, matmultot,idum,ldum,restart_lun, lun
  integer             :: nshell,natoms,I,J,bast,idummy, lupri, luerr
  type(configItem)    :: config
  logical             :: doDFT,mem_monitor
  real(realk) :: DUMMY(1,1)


  mem_monitor = .false. !Mostly for memory debugging
  call init_globalmemvar !initialize the global memory counters
! Read LSDALTON.INP FILE,MOLECULE.INP AND BASISSET FILES TO SET UP
! THE  DALTONINPUT STRUCTURE
  LUPRI=6
  LUERR=6

  call mat_no_of_matmuls(no_of_matmuls)
  ! READ *LINSCA SECTION UNDER **WAVEFUNCTION AND INIT SCF_CONFIG
  call config_set_default_config(config)
  ! Setting default starting-guess to ATOMS instead of HUCKEL
  config%opt%cfg_start_guess = 'ATOMS'

  call config_read_input(config,lupri,luerr)
  ls%input%dalton = config%integral

  doDFT = config%opt%calctype.EQ.config%opt%dftcalc
  call ls_init(ls,lupri,luerr,nbast,config%integral,doDFT,config%doDEC,.true.)
  call set_final_config_and_print(lupri,config,ls,nbast)

  ! eventually empirical dispersion correction in case of dft
  CALL II_DFT_DISP(LS%SETTING,0,DUMMY,LUPRI)

  call mat_pass_info(LUPRI,config%opt%info_matop,mem_monitor)

#if defined(VAR_INT64)
  if (bit_size(nbast)==64) then
     write(lupri,*) 'Using 64-bit integers!'
  endif
#endif

 !Initialize grand canonical electronic configuration data table

 call EcData_init 

  if (config%decomp%cfg_gcbasis) then
       call trilevel_basis(config%opt,ls)
  endif

   call write_lsitem_to_disk(ls)

 

end program make_lsitem
