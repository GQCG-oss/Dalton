module localization_input
use davidson_settings
use configurationType, only: configitem
use ls_util
use kurtosis!, only: KurtosisItem
use precision
use matrix_module!, only: matrix
use loc_types!, only: OrbitalLoc
use decompMod!, only: orbspread_data, decompitem
use queue_module!, only: modFIFO



CONTAINS
 !> \brief Read the **LOCALIZE ORBITALS input section in LSDALTON.INP and set 
 !> configuration structure accordingly.
 !> \author Ida-Marie Hoeyvik
 !> \date May 2013

subroutine orbitalloc_input(input,output,config,READWORD,word)
implicit none
 type(ConfigItem) :: config
 !> Logical for keeping track of when to read
 LOGICAL,intent(inout) :: READWORD
 !> Logical unit number for LSDALTON.INP
 integer,intent(in) :: input
 !> Logical unit number for LSDALTON.OUT
 integer,intent(in) :: output
 !> Word read from input
 character(len=80),intent(inout) :: word

    DO

       IF(READWORD) THEN
          READ (input, '(A80)') WORD
          call capitalize_string(word)
          READWORD=.TRUE.
       ENDIF

       IF ((WORD(1:1) .EQ. '!') .OR. (WORD(1:1) .EQ. '#')) CYCLE

       IF(WORD(1:2) .EQ. '**') THEN
          READWORD=.FALSE.
          EXIT
       ENDIF

       IF(WORD(1:13) == '*END OF INPUT') THEN
          EXIT
       END IF

       ORBLOC_INPUT_INFO: SELECT CASE(WORD)
       ! Localization function
       case('.PSM')
           config%decomp%cfg_mlo = .true.
           config%davidOrbLoc%orbspread=.true.
           config%davidOrbLoc%linesearch=.true.
           READ(input,*) config%decomp%cfg_mlo_m(1), config%decomp%cfg_mlo_m(2)
       case('.PFM')
           config%decomp%cfg_mlo = .true.
           config%davidOrbLoc%PFM = .true.
           READ(input,*) config%decomp%cfg_mlo_m(1), config%decomp%cfg_mlo_m(2)
           config%davidOrbLoc%PFM_input%crossterms=.true.
       case('.PIPEKMEZEY')
           config%davidOrbLoc%PM_input%PipekMezeyLowdin=.true.
           config%davidOrbLoc%PM=.true.
           config%decomp%cfg_mlo = .true.
           config%davidOrbLoc%NOL2OPT = .true.
           config%decomp%cfg_mlo_m(1)=2
           config%decomp%cfg_mlo_m(2)=2
       case('.OCCPIPEKMEZEY')
           config%davidOrbLoc%PM_input%PipekMezeyLowdin=.true.
           config%davidOrbLoc%PM=.true.
           config%decomp%cfg_mlo = .true.
           config%davidOrbLoc%NOL2OPT = .true.
           config%decomp%cfg_mlo_m(1)=2
           config%decomp%cfg_mlo_m(2)=0
       case('.VIRTPIPEKMEZEY')
           config%davidOrbLoc%PM_input%PipekMezeyLowdin=.true.
           config%davidOrbLoc%PM=.true.
           config%decomp%cfg_mlo = .true.
           config%davidOrbLoc%NOL2OPT = .true.
           config%decomp%cfg_mlo_m(1)=0
           config%decomp%cfg_mlo_m(2)=2
       CASE('.PIPEKM(MULL)')
           config%davidOrbLoc%PM_input%PipekMezeyMull=.true.
           config%davidOrbLoc%PM=.true.
           config%decomp%cfg_mlo = .true.
           config%davidOrbLoc%NOL2OPT = .true.
           config%decomp%cfg_mlo_m(1) = 2
           config%decomp%cfg_mlo_m(2) = 2
       CASE('.OCCPIPEKM(MULL)')
           config%davidOrbLoc%PM_input%PipekMezeyMull=.true.
           config%davidOrbLoc%PM=.true.
           config%decomp%cfg_mlo = .true.
           config%davidOrbLoc%NOL2OPT = .true.
           config%decomp%cfg_mlo_m(1) = 2
           config%decomp%cfg_mlo_m(2) = 0
       CASE('.VIRTPIPEKM(MULL)')
           config%davidOrbLoc%PM_input%PipekMezeyMull=.true.
           config%davidOrbLoc%PM=.true.
           config%decomp%cfg_mlo = .true.
           config%davidOrbLoc%NOL2OPT = .true.
           config%decomp%cfg_mlo_m(1) = 0
           config%decomp%cfg_mlo_m(2) = 2

       ! SPECIFIC SETTINGS
       CASE('.NO LEVEL2 LOCALIZATION')
           config%decomp%cfg_mlo = .true.
           config%davidOrbLoc%NOL2OPT = .true.
       CASE('.ONLY LOC')
           config%davidOrbLoc%OnlyLocalize=.true.
       CASE('.NOPRECOND')
           config%davidOrbLoc%precond=.false.
           config%davidSCF%precond=.false.
           config%davidOrbLoc%PM_input%precond=.false.
       CASE('.MACRO IT')
           READ(input,*)  config%davidOrbLoc%max_macroit
       CASE('.NOT_QUIT_AFTER_10IT')
           config%davidOrbLoc%quit_after_10it = .false.
       CASE('.LOOSE MICRO THRESH')
           config%davidOrbLoc%conv_thresh= 0.1_realk
           config%davidOrbLoc%global_conv_thresh= 0.1_realk
           config%davidOrbLoc%local_conv_thresh= 0.01_realk
       CASE('.RESTART LOCALIZATION')
           config%davidOrbLoc%orbloc_restart = .true.
           config%davidOrbLoc%OnlyLocalize=.true.
       CASE('.SAVE ITERATIONS')
           READ(input,*) config%davidOrbLoc%orbital_save_interval
       ! TESTING AND DEBUG OPTIONS
       CASE('.TEST PFM')
           config%davidOrbLoc%PFM_input%TESTCASE = .true.
           config%decomp%cfg_mlo = .true.
       CASE('.ORBLOC DEBUG')
           config%davidOrbLoc%orb_debug = .true.
           config%davidOrbLoc%PM_input%orb_debug = .true.


       ! PRINTING AND PLOTTING OPTIONS 
       CASE('.ORBITAL LOCALITY')
           config%davidOrbLoc%all_orb_locality=.true.
       CASE('.ORBITAL PLOT')
           config%davidOrbLoc%make_orb_plot=.true.
           READ(input,*) config%davidOrbLoc%plt_orbital

       CASE DEFAULT
          WRITE (output,'(/,3A,/)') ' Keyword "',WORD,&
               & '" not recognized in orbitalloc_input'
          CALL lsQUIT('Illegal keyword in orbitalloc_input',output)

       END SELECT ORBLOC_INPUT_INFO

    ENDDO

end subroutine orbitalloc_input


end module

