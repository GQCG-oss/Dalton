module LS_optimizer_mod
  use precision
  use matrix_module!, only: matrix
  use TYPEDEFTYPE!, only: lsitem
  use configurationType!, only: configitem
  use optimization_input!, only: opt_setting, MXCOOR, MXCENT
  use memory_handling!, only: mem_alloc, mem_dealloc
  use ls_util!, only: lsheader, ls_print_gradient
  use files!, only: lsopen, lsclose
  use molecule_module!, only: print_geometry
  use lstiming!, only: lstimer
  use WRITEMOLEFILE!, only: write_molecule_output
  use Energy_and_deriv!, only: get_energy, get_gradient, get_num_grad
  use Fundamental  
#ifdef VAR_DEC
  use DEC_settings_mod
#endif
  use dec_typedef_module
  use lsdalton_rsp_mod,only: LS_RSP_EQ_SOL_EMPTY
!===========================================================!
!    The main driver for geometry optimization in LSDALTON  !
!===========================================================!
! Written by Vladimir Rybkin in 11/2010
!
  private
  public :: LS_RUNOPT
CONTAINS
  SUBROUTINE LS_RUNOPT(E,config,H1,F,D,S,CMO,ls)
    Implicit none
    !  All these general entities needed to get energy and gradient
    Type(lsitem)                    :: ls   ! General information,used only to get E and gradient
    Type(Matrix), intent(inout)     :: F(1),D(1),S ! Fock,density,overlap matrices
    Type(Matrix), intent(inout)     :: H1   ! One electron matrix
    Type(Matrix), intent(inout)     :: CMO       ! Orbitals
    Type(ConfigItem), intent(inout) :: Config ! General information
    Real(realk),intent(inout)       :: E(1)   ! Energy
    !
    INTEGER :: iOpt, ipre
    Real(realk) :: GradThr, ThrStep, ThGradMax, ThStepMax, trsave
    ! Special preoptimization step currently only possible for DEC calculations
    config%noDecEnergy = .TRUE.
    IF (config%optinfo%dynopt) THEN
      IF (.NOT.DECinfo%doDec) THEN
        CALL LSQUIT('.DYNOPT currently just available for DEC-type calculations',ls%lupri)
      ENDIF
      trsave    = config%optinfo%TrustRad
      GradThr   = config%optinfo%GradThr
      ThrStep   = config%optinfo%ThrStep
      ThGradMax = config%optinfo%ThGradMax
      ThStepMax = config%optinfo%ThStepMax
      ipre=1
      DO 
        config%optinfo%dynamicThreshold = ipre.NE.nFOTs
        IF (DECinfo%doDec) THEN
           DECinfo%FOT = DECinfo%GeoFOTs(Ipre)
           DECinfo%FOTlevel=ipre
        ENDIF
        WRITE(ls%lupri,'(A,I2)') '*** Starting dynamical optimization step number',ipre
        ! Reset trust radius to input value for each dynamical optimization step
        config%optinfo%TrustRad  = trsave
        config%optinfo%ItrNmr    = 0
        config%optinfo%GradThr   = GradThr
        config%optinfo%ThrStep   = ThrStep
        config%optinfo%ThGradMax = ThGradMax
        config%optinfo%ThStepMax = ThStepMax
        CALL LS_FLSHFO(ls%lupri)
        CALL LS_RUNOP1(E,config,H1,F,D,S,CMO,ls)
        IF (config%optinfo%dynamicConvergence) EXIT
        IF (ipre==nFOTs) EXIT
        Ipre = Ipre +1
      ENDDO
      IF (config%optinfo%dynamicConvergence) THEN
        WRITE(ls%lupri,'(A)') '*** Dynamical optimization converged'
      ELSE
        WRITE(ls%lupri,'(A)') '*** Dynamical optimization did not converge!'
      ENDIF
      !Reset FOT level
      DECinfo%FOTlevel=1
      DECinfo%FOT = DECinfo%GeoFOTs(1)
    ! Deafult geoemetry optimization
    ELSE
      CALL LS_RUNOP1(E,config,H1,F,D,S,CMO,ls)
    ENDIF
    config%noDecEnergy = .FALSE.
    Call mem_dealloc(config%optinfo%IConstr)
  END SUBROUTINE LS_RUNOPT

  SUBROUTINE LS_RUNOP1(E,config,H1,F,D,S,CMO,ls)
    Implicit none
    !  All these general entities needed to get energy and gradient
    Type(lsitem) :: ls   ! General information,used only to get E and gradient
    Type(Matrix), intent(inout) :: F(1),D(1),S ! Fock,density,overlap matrices
    Type(Matrix), intent(inout) :: H1   ! One electron matrix
    Type(Matrix), intent(inout) :: CMO       ! Orbitals
    Type(ConfigItem), intent(inout) :: Config ! General information
    !
    LOGICAL ls_minend, INDXOK, STATPO, ACTIVE
    LOGICAL REJGEO, TRU, FAL, TMPLOG, NEWSTP, NEWBMT
    CHARACTER TMPLIN*80, WRDRSP*7
    Integer :: lupri, luerr   ! File units
    Real(realk) :: E(1),Eerr   ! Energy
    Real(realk), pointer ::  EGRAD(:), CSTEP(:)
    Real(realk), pointer ::  GRDOLD(:), GRDMAT(:,:)
    Real(realk), pointer ::  STPMAT(:,:), HESOLD(:,:)
    Real(realk), pointer ::  GRDARR(:,:), STPARR(:,:)
    ! Reduced Hessian
    Real(realk), pointer :: Red_Hess(:,:)
    Real(realk), pointer :: Red_Cart_Hess(:,:)
    !
    !     The array geinfo contains optimization information for each
    !     iteration. The first index is the iteration, the second gives
    !     the property:   1  -  Energy
    !                     2  -  Gradient norm
    !                     3  -  Index of Hessian
    !                           (a negative index indicates symmetry break)
    !                     4  -  Step length
    !                     5  -  Trust radius
    !                     6  -  # rejected steps
    !
    Real(realk), pointer :: GEINFO(:,:)
    !
    Real(realk), pointer :: WILBMT(:,:), BMTRAN(:,:)
    Real(realk), pointer :: HESINT(:,:), VECMOD(:)
    Real(realk) :: TE,TS ! CPU time
    Integer, PARAMETER :: IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12
    !
    Real(realk), pointer :: KEHESS(:),KALHES(:)
    Real(realk), pointer :: KBMINV(:),KPJINM(:),KEVEC(:), KCONMT(:)
    Real(realk), pointer :: KTEMP1(:),KTEMP2(:),KTEMP3(:),KTEMP4(:), &
         &   KTEMP5(:),KTEMP6(:),KTEMP7(:),KTEMP8(:),KTEMP9(:)
    Real(realk), pointer :: TEMP1(:,:),TEMP3(:,:),TEMP2(:,:)
    Integer,pointer :: ITEMP2(:,:)
    ! 
    CHARACTER(len=70) :: WORD
    LOGICAL :: ReadWord
    Integer :: luinfo,fileStatus
    Integer :: NAtoms,i,j,IPrint,IWOFF,NCRDHS,JI,IREJ,NCRD
    Real(realk) :: THRLDP,THRIND,TOLST,EMOD,ERGDIF
    Integer :: NRedint

    lupri = config%lupri
    luerr = config%luerr
    ! Quit statement for failed print level
    If ((config%optinfo%RedInt .OR. config%optinfo%DelInt) .AND. config%optinfo%IPrint .GE. 15) then
       Call LSQuit('The print level in **OPTIMI is too high for optimization &
            & redundant internals. Try to reduce it', lupri)
    Endif
    !
    ! First, grabbing number of atoms
    NAtoms = config%Molecule%nAtoms
    ! Second, determining numbers of various coordinates
    ! needed for memory allocation
    MXCENT = NAtoms
    MXCOOR = 3*NAtoms
    config%optinfo%ICartCoord = 3*NAtoms
    config%optinfo%NTempMat = 6*config%optinfo%ICartCoord
    config%optinfo%NCoordTot = config%optinfo%ICartCoord
    config%optinfo%energy = E(1)
    ! Knowing number of cartesians,allocating cartesian coordinates vector
    Call mem_alloc(config%optinfo%Coordinates,3,MXCENT)
    Call ls_dzero(config%optinfo%Coordinates,3*MXCENT)
    !
    ! Third, grabbing the initial coordinates
    !
    Do i = 1, NAtoms
       config%optinfo%Coordinates(:,i) = config%Molecule%Atom(i)%Center(:)  
    Enddo
    ! Do the force modification of energy if asked
    If (config%optinfo%FMPES) call FM_energy(E(1),config%optinfo)
    config%optinfo%energy = E(1)
    ! Fourth, finding the number of internals if needed
    If (config%optinfo%RedInt .OR. config%optinfo%DelInt .OR. &
     &  config%optinfo%InmdHess.OR. config%optinfo%InrdHess) then  
       ! Allocate some memory for finding number of redundant internals
       Call mem_alloc(TEMP1,NAtoms,8)
       Call mem_alloc(ITEMP2,NAtoms,NAtoms)
       ! We allocate INTCRD with memory in excess, once we get
       ! the exact number of redundant internals we will allocate
       ! as much memory as needed
       Call mem_alloc(config%optinfo%INTCRD,NAtoms*3*8,6)
       ! Find redundant internals
       IPrint = config%optinfo%IPrint   ! To avoid a printout if LS_FNDRED
       config%optinfo%IPrint = 0
       Call LS_FNDRED(TEMP1,ITEMP2,NAtoms,Config%Molecule,lupri,config%optinfo,.TRUE.)
       config%optinfo%IPrint = IPrint   ! Print level set back
       ! Maximum number of coordinates set
       MXCOOR = MAX(config%optinfo%ICartCoord,config%optinfo%NIntCoord)
       ! Deallocate memory
       Call mem_dealloc(TEMP1)
       Call mem_dealloc(ITEMP2)
       Call mem_dealloc(config%optinfo%INTCRD)
       ! If .Findre no optimization is carried out, but only the redundant internals are found
       If (config%optinfo%FindRe) then
          Call lsheader(lupri,'Redundant internals found')
          Call mem_dealloc(config%optinfo%Coordinates)
          Return
       Endif
    Endif

    ! Allocating gradient and Hessian
    Call mem_alloc(config%optinfo%GradMol,config%optinfo%ICartCoord)
    Call mem_alloc(config%optinfo%HessMol,config%optinfo%ICartCoord,config%optinfo%ICartCoord)
    !Initializing them
    Call ls_dzero(config%optinfo%GradMol,config%optinfo%ICartCoord)
    Call ls_dzero(config%optinfo%HessMol,config%optinfo%ICartCoord*config%optinfo%ICartCoord)
    ! Allocating more
    Call mem_alloc(config%optinfo%STPDIA,MXCOOR)
    Call mem_alloc(config%optinfo%STPSYM,MXCOOR)
    Call mem_alloc(config%optinfo%GRDDIA,MXCOOR)
    Call mem_alloc(config%optinfo%EVAL,  MXCOOR)
    Call mem_alloc(config%optinfo%EVALOL,MXCOOR)
    Call mem_alloc(GRDMAT,MXCOOR,25)
    Call mem_alloc(STPMAT,MXCOOR,25)
    Call mem_alloc(GRDARR,MXCOOR,25)
    Call mem_alloc(STPARR,MXCOOR,25)
    Call mem_alloc(EGRAD,MXCOOR)
    Call mem_alloc(GRDOLD,MXCOOR)
    If (config%optinfo%RedInt .OR. config%optinfo%DelInt .OR. &
     &  config%optinfo%InrdHess .OR. config%optinfo%InmdHess) then 
       Call mem_alloc(config%optinfo%GRDINT,MXCOOR)
       Call mem_alloc(config%optinfo%STPINT,MXCOOR)
       Call mem_alloc(config%optinfo%CoordInt,MXCOOR)
       Call mem_alloc(config%optinfo%INTCRD,8*MXCOOR,6)
       Call mem_alloc(HESINT,MXCOOR,MXCOOR)
       Call mem_alloc(HESOLD,MXCOOR,MXCOOR)
    else
       Call mem_alloc(HESINT,MXCOOR,MXCOOR)
       Call mem_alloc(HESOLD,MXCOOR,MXCOOR)
    Endif
    Call mem_alloc(BMTRAN,MXCOOR,MXCOOR)
    Call mem_alloc(KBMINV,MXCOOR*MXCOOR)
    Call mem_alloc(WILBMT,MXCOOR,MXCOOR)
    Call mem_alloc(VECMOD,MXCOOR)
    !
    !     Allocate GEINFO 
    !
    call mem_alloc(GEINFO,config%optinfo%MaxIter+1,6,.TRUE.,.FALSE.)
    call LSTIMER('START ',TS,TE,lupri)
    !
    !     Initialization of variables.
    !
    THRLDP = 1.0E-4_realk
    THRIND = 5.0E-4_realk
    TOLST  = 1.0E-5_realk
    call ls_DZERO(GEINFO,(config%optinfo%MaxIter+2)*6)
    call ls_DZERO(GRDARR,25*MXCOOR)
    call ls_DZERO(STPARR,25*MXCOOR)
    ACTIVE = .FALSE.
    config%optinfo%KEPTIT = 0
    GEINFO(0,5) = config%optinfo%TrustRad
    IWOFF = 0
    !
    !     Perform preoptimization if requested.
    !
    IF (config%optinfo%DoPre) call lsquit('DOPRE not an option in lsdalton',lupri)
    !
    !     Allocating some memory
    !
    call mem_alloc(KEHESS,MXCOOR*MXCOOR)
    call mem_alloc(KALHES,MXCOOR*MXCOOR)
    !
    !     Calculate gradient and Hessian for second order method and
    !     first order method with initial Hessian.
    !
    If (.NOT. config%optinfo%Findre .AND. .NOT. config%optinfo%simple_scan) then
       !
       !     First order methods only require the energy and the gradient.
       !
       call Obtain_Gradient(E(1),Eerr,lupri,NAtoms,S,F(1),D(1),ls,config,CMO,config%optinfo)
       config%optinfo%energy = E(1)
    Endif
    !
    !     Make config%optinfo%VRLM-file of initial geometry if requested.
    !
    IF (config%optinfo%VRLM) CALL LSQUIT('No VRLM implemented!',LUPRI)
    !
    !     Save initial geometry and energy to MOLDEN file
    !
    !      IF (MOLDEN) call LSQUIT('No MOLDEN implemented!',LUPRI)

    !
    !     We allocate more 
    !
    Call mem_alloc(CSTEP,MXCOOR)
    Call mem_alloc(KPJINM,MXCOOR*MXCOOR)
    Call mem_alloc(KEVEC,MXCOOR*MXCOOR)
    Call mem_alloc(KCONMT,MXCOOR*MXCOOR)
    Call mem_alloc(KTEMP1,MXCOOR*MXCOOR)
    Call mem_alloc(KTEMP2,MXCOOR*MXCOOR)
    Call mem_alloc(KTEMP3,MXCOOR*MXCOOR)
    Call mem_alloc(KTEMP4,MXCOOR*MXCOOR)
    Call mem_alloc(KTEMP5,MXCOOR*MXCOOR)
    Call mem_alloc(KTEMP6,MXCOOR*MXCOOR)
    Call mem_alloc(KTEMP7,MXCOOR*MXCOOR)
    Call mem_alloc(KTEMP8,MXCOOR*MXCOOR)
    Call mem_alloc(KTEMP9,MXCOOR*MXCOOR)
    !     Set cartesian step equal to zero
    call ls_dzero(CSTEP,MXCOOR)
    !
    !     Check if redundant internal coordinates should be used.
    !

    IF (config%optinfo%DelInt .OR. config%optinfo%RedInt .OR. &
         &   config%optinfo%InrdHess .OR. config%optinfo%InmdHess) THEN
       call ls_INIRED(config%optinfo%NIntCoord,MXCOOR,WILBMT,BMTRAN,KBMINV, &
            &        KPJINM,KTEMP1,KTEMP2,KTEMP3,KTEMP4,KTEMP5,KTEMP6, &
            Config%Molecule,NAtoms,lupri,config%optinfo)
    END IF
!
! Branching for doing PES scan
!    
If (config%optinfo%simple_scan) then
      call mem_alloc(config%optinfo%Scan_info,config%optinfo%scanMaxIter+1,2)
      Call Merge_geometry(config%optinfo,config%optinfo%ScanCoord)
      Call Rotate_coordinates(MXCENT,config%optinfo%Coordinates,config%optinfo%Atoms_to_move(1,1),&
           & config%optinfo%Atoms_to_move(2,1))
      Call PES_scan(ls,config,config%optinfo,H1,F,D,S,CMO,E,NAtoms,lupri,luerr)
      ! Final printout
      call Print_scan(config%optinfo,config%optinfo%scanMaxIter,lupri)
      ! Deallocate
      call mem_dealloc(config%optinfo%Atoms_to_move)
      call mem_dealloc(config%optinfo%Scan_info)
Else
!
!  Entering normal optimizer
!
    
    IF (config%optinfo%DelInt .OR. config%optinfo%RedInt) THEN
       NCRDHS = config%optinfo%NIntCoord
    ELSE
       NCRDHS = config%optinfo%ICartCoord
    END IF
    IF (config%optinfo%RatFun) NCRDHS = NCRDHS + 1
    !     Numeric Hessian in reduced space if asked
    ! Internal
    If (config%optinfo%RedSpa) then
       call mem_alloc(Red_Hess,config%optinfo%Hess_dim,config%optinfo%Hess_dim)
       call Num_Int_Hess(ls,config,F,D,S,H1,CMO,config%optinfo,Red_Hess, &
            &config%optinfo%Hess_dim,config%optinfo%NIntCoord,MXCOOR,&
            &NAtoms,KBMINV,WILBMT,lupri,luerr)
    Endif
    ! Cartesian
    If (config%optinfo%CartRS) then
       call mem_alloc(Red_Cart_Hess,config%optinfo%Red_Atoms*3,config%optinfo%Red_Atoms*3) 
       call Num_Cart_Hess(ls,config,F,D,S,H1,CMO,config%optinfo%Red_Atoms, &
            &NAtoms,Red_Cart_Hess,config%optinfo,lupri,luerr)
       call mem_dealloc(Red_Cart_Hess)
    Endif
    !
    !     Initialize Hessian if first order method is used.
    !
7   CONTINUE
    IF (.NOT. config%optinfo%Newton .AND. .NOT. config%optinfo%simple_scan) THEN
       If (config%optinfo%Redint .OR. config%optinfo%Delint) then
          call ls_INIHES(config%Molecule,config%optinfo%NIntCoord, &
               &  config%optinfo%NIntCoord,MXCOOR,GRDOLD,HESOLD, &
               &  KTEMP1,KTEMP2,KTEMP3,KTEMP4,WILBMT,BMTRAN, & 
               &  KBMINV,HESINT,lupri,config%optinfo)
       ELSE
          call ls_INIHES(config%Molecule,MXCOOR, &
               &  config%optinfo%NIntCoord,MXCOOR,GRDOLD,HESOLD, &
               &  KTEMP1,KTEMP2,KTEMP3,KTEMP4,WILBMT,BMTRAN, & 
               &  KBMINV,HESINT,lupri,config%optinfo)

       ENDIF
    ENDIF
    !
    !     For REDSPA: build the full Hessian      
    !
    IF (config%optinfo%RedSpa) THEN
       call Add_Redspa(config%optinfo,MXCOOR,config%optinfo%Hess_Dim,Red_Hess,HESINT,lupri)
       call mem_dealloc(Red_Hess)
       ! Transform Hessian to Cartesian frame
       call ls_DZERO(config%optinfo%HessMol,MXCOOR*MXCOOR)
       ! Some memory
       call mem_alloc(TEMP1,MXCENT,8)
       TEMP1 = 0E0_realk
       call mem_alloc(TEMP2,MXCOOR,MXCOOR)
       call mem_alloc(TEMP3,MXCOOR,MXCOOR)
       call ls_HQ2HX(config%Molecule,MXCOOR,MXCOOR,TEMP1,TEMP2,TEMP3,HESINT, &
            &     config%optinfo%GRDINT,config%optinfo%HessMol,WILBMT,BMTRAN, &
            config%optinfo)
       ! Free memory
       call mem_dealloc(TEMP1)
       call mem_dealloc(TEMP2)
       call mem_dealloc(TEMP3)
    ENDIF
    !
    call ls_DZERO(EGRAD,MXCOOR)
    call ls_DZERO(KEHESS,MXCOOR*MXCOOR)
    call ls_DZERO(KALHES,MXCOOR*MXCOOR)
!
! If we do simple scan all the following is not needed
!
If (.NOT. config%optinfo%simple_scan) then
    DO i = 1, config%optinfo%ICartCoord
       EGRAD(i) = config%optinfo%GradMol(i)
    ENDDO
    JI = 1
    DO I = 1, config%optinfo%ICartCoord
       DO J = 1, config%optinfo%ICartCoord
          KEHESS(JI) = config%optinfo%HessMol(J,I) 
          KALHES(JI) = config%optinfo%HessMol(J,I) 
          JI = JI + 1
       ENDDO
    ENDDO
    !      Call ls_copyGH(EGRAD,KEHESS,KALHES,config%optinfo%GradMol,config%optinfo%HessMol,config%optinfo%ICartCoord)
    !     Construct projection operator and use it.
    !     Then diagonalize Hessian.
    !
    IF ((config%optinfo%RedInt .OR. config%optinfo%DelInt)) THEN!.AND. (.NOT. config%optinfo%RedSpa)) THEN

       IF (config%optinfo%Newton) call ls_CGHINT(config%Molecule, &
            &        config%optinfo%NIntCoord,MXCOOR,KTEMP1, &
            &        KTEMP2,KTEMP3,KTEMP4,KTEMP5, &
            &        WILBMT,KBMINV,BMTRAN,HESINT,lupri,config%optinfo)
       call ls_PRJINT(config%optinfo%NIntCoord,config%optinfo%NIntCoord,KPJINM,KCONMT, &
            &        HESINT,KTEMP1,KTEMP2,KTEMP3,KTEMP4,lupri,config%optinfo)
       !
       !     Note that the contents of KTEMP7 is passed on
       !     from LINSRC to FNSTIN below.
       !
       IF (config%optinfo%LnSearch .AND. config%optinfo%RatFun .AND. (config%optinfo%ItrNmr .GT. 0)) &
            &        call ls_LINSRC(config%optinfo%NIntCoord,config%optinfo%NIntCoord,config%optinfo%GRDINT,GRDARR(1,1), &
            &        KTEMP7,STPARR(1,1),KTEMP3,KTEMP4, &
            &        ACTIVE,EMOD,lupri,config%optinfo)
       IF (config%optinfo%RatFun .AND. config%optinfo%Saddle) NCRDHS = NCRDHS - 1
       call ls_DIAINT(config%optinfo%NIntCoord,MXCOOR,NCRDHS,KEVEC,KTEMP1, &
            &        KTEMP2,KTEMP3,KTEMP4,THRIND,HESINT,KTEMP5,lupri,config%optinfo)
       IF (config%optinfo%RatFun .AND. config%optinfo%Saddle) NCRDHS = NCRDHS + 1
    ELSE
       !
       !     Note that the contents of KTEMP1 is passed on
       !     from PROJGH to DIAHES below.
       !
       call ls_PROJGH(EGRAD,KEHESS,KALHES,KTEMP1, &
            &   KTEMP2,KTEMP3,KTEMP4,config%optinfo%NCoordTot, &
            &   config%optinfo%NTempMat,lupri,config%optinfo)
       IF (config%optinfo%LnSearch .AND. config%optinfo%RatFun .AND. (config%optinfo%ItrNmr .GT. 0)) &
            &      call ls_LINSRC(config%optinfo%ICartCoord,MXCOOR,EGRAD,GRDARR(1,1),CSTEP, &
            &        STPARR(1,1),KTEMP3,KTEMP4,ACTIVE,EMOD,lupri,config%optinfo)
       call ls_DIAHES(config%optinfo%NIntCoord,MXCOOR,NCRDHS,EGRAD,KEHESS,KALHES, &
            & KTEMP1,THRIND, &
            & KEVEC,KTEMP2,KTEMP3,KTEMP4,config%optinfo%NCoordTot,config%optinfo%NTempMat, &
            & lupri,config%optinfo)
    END IF
!
ENDIF  ! If not simple scan
    GEINFO(0,1) = E(1)
    GEINFO(0,3) = config%optinfo%IndHes*1.0E0_realk
!!!!!! Vladimir:: currently disabled
    !     Write Hessian to file (for 1st order restarts).
    !
    !      IF (.NOT. config%optinfo%NoHessianWrite) &
    !     &     call ls_PNCHES(config%optinfo%NIntCoord,MXCOOR,HESINT,WILBMT,BMTRAN,KTEMP1, &
    !     &     KTEMP2,KTEMP3,KTEMP4,WORK(KWRK2),LWRK2,lupri,config%optinfo)

    !
    !     Determine step, check for convergence, print output and
    !     and update geometry.
    !
    IREJ = 0
755 CONTINUE
    IF (config%optinfo%RedInt .OR. config%optinfo%DelInt) THEN
       call ls_FNSTIN(config%Molecule,config%optinfo%NIntCoord, &
            &        MXCOOR,NCRDHS,HESINT,KEVEC, &
            &        KTEMP1,KTEMP2,KTEMP3,KTEMP4, &
            &        KTEMP5,CSTEP,WILBMT,BMTRAN,KBMINV,GRDARR, &
            &        STPARR,ACTIVE,EMOD,VECMOD,KTEMP7,lupri,config%optinfo)
    ELSE
       call ls_FNDSTP(MXCOOR,MXCOOR,NCRDHS,EGRAD,KEHESS, &
            &        KEVEC,KTEMP1,KTEMP2,KTEMP3, &
            &        KTEMP4,KTEMP5,CSTEP,GRDARR,STPARR, &
            &        ACTIVE,EMOD,VECMOD,lupri,config%optinfo)
    END IF
!
! If we don't do a simple scan then we don't care about convergence
!
If (.NOT. config%optinfo%simple_scan) then 
    config%optinfo%GeConv = ls_minend(config%optinfo%NIntCoord,BMTRAN,KTEMP1,KTEMP2,lupri,config%optinfo)
    !
    !     If there has been a completely failed step, the geometry has
    !     by default not converged.
    !
    IF (ABS(GEINFO(0,6)) .GT. 1.0E-3_realk) config%optinfo%GeConv = .FALSE.
    call ls_PRIALL(config%Molecule,NAtoms,CSTEP,KTEMP1,lupri,config%optinfo)
    !
    !     Save this geometry and energy to MOLDEN file
    !
    NEWSTP = .FALSE.
    NEWBMT = .FALSE.
    !
    !     To allow reinitialization
    !
    config%optinfo%InitHess = .FALSE.
    !
Else  ! simple scan
    config%optinfo%GeConv = .FALSE.
    call ls_PRIALL(config%Molecule,NAtoms,CSTEP,KTEMP1,lupri,config%optinfo)
Endif
!    
    IF (.NOT. config%optinfo%GeConv) then 
       config%optinfo%dynamicChange = .FALSE.
       call Find_Geometry(E,CSTEP,EGRAD,KTEMP1, &
            &     KTEMP2,IREJ,GEINFO,NEWSTP,NEWBMT,NAtoms,lupri, &
            luerr,config,ls,H1,F,D,S,CMO,config%optinfo)
       IF (config%optinfo%dynamicChange) goto 300
       !    Grabbing coordinates
       Do i = 1, NAtoms
          config%optinfo%Coordinates(:,i) = KTEMP1(3*i-2:3*i)
       Enddo
       !    Renovating Config%Molecule
       Do i = 1,NAtoms
          Config%Molecule%Atom(i)%Center(:)=config%optinfo%Coordinates(:,i)
       Enddo
    ENDIF
    !
    IF (NEWSTP) GOTO 755
    GEINFO(0,2) = config%optinfo%GradNorm
    GEINFO(0,4) = config%optinfo%StepNorm
    IF (config%optinfo%ItrNmr .LT. config%optinfo%MaxIter) GEINFO(1,5) = config%optinfo%TrustRad
    IF (ABS(GEINFO(0,6)) .LT. 1.0E-3_realk) THEN
       GEINFO(0,6) = IREJ*1.0E0_realk
    ELSE
       GEINFO(0,6) = -(ABS(GEINFO(0,6))+ABS(IREJ)*1.0E0_realk)
    END IF
    config%optinfo%TotRj = config%optinfo%TotRj + ABS(IREJ)
    !
    !     Determine value of the various coordinates
    !
    IF (config%optinfo%RedInt .AND. (config%optinfo%IPrint .GE. 1)) THEN
       call ATOM_INI(KTEMP1,Config%Molecule,config%optinfo,NAtoms,.TRUE.,lupri)
       call ls_GETINT(NAtoms,config%optinfo%NIntCoord,KTEMP1,config%optinfo%CoordInt,lupri,config%optinfo)
       call lsheader(lupri,'New internal coordinates') 
       call ls_output(config%optinfo%CoordInt,1,1,1,config%optinfo%NIntCoord,1,config%optinfo%NIntCoord,1,LUPRI)
       WRITE(LUPRI,'(//)')
    END IF
    !
    !     If the step has failed
    !
    IF (IREJ .LT. 0) THEN
       GOTO 7
       !      ELSE IF (config%optinfo%RejIni .AND. config%optinfo%RedInt .AND. (config%optinfo%TotRj .GE. 3)) THEN
       !         WRITE(LUPRI,*)'***** NOTE! *****'
       !         WRITE(LUPRI,*)
       !     &        'The number of dihedral angles will be reduced!'
       !         call ls_RREDUN
       !         config%optinfo%TotRj = 0
    END IF
    !
    IF (config%optinfo%GeConv .AND. config%optinfo%DoPre .AND. (.NOT. config%optinfo%FinPre)) THEN
       config%optinfo%KeepHessian = .FALSE.
    END IF
    !
    !     Switch the iteration counter
    !
    !      config%optinfo%ItrNmr = config%optinfo%ItrNmr + 1
    !
    !     DO WHILE-loop that runs until geometry has converged or
    !     maximum number of iterations is reached.
    !
10  CONTINUE
    call LSTIMER('Geom. opt.',TS,TE,lupri)

    IF ((config%optinfo%ItrNmr .LT. config%optinfo%MaxIter) .AND. (.NOT. config%optinfo%GeConv)) THEN
       config%optinfo%ItrNmr = config%optinfo%ItrNmr + 1
       NCRD = config%optinfo%NCoordTot
       IF (config%optinfo%RedInt .OR. config%optinfo%DelInt) NCRD = config%optinfo%NIntCoord
       DO I = 1, NCRD
          config%optinfo%EVALOL(I) = config%optinfo%EVAL(I)
       ENDDO
       ! If we do a scan we don't need to update step and we don't need a gradient
       IF (.NOT. config%optinfo%simple_scan) THEN
          IF (config%optinfo%RedInt .OR. config%optinfo%DelInt) THEN
             call ls_UPGDST(config%optinfo%NIntCoord,config%optinfo%NIntCoord,GRDARR,STPARR, &
               & config%optinfo%GRDINT,config%optinfo%STPINT,lupri,config%optinfo)
          ELSE
             call ls_UPGDST(config%optinfo%ICartCoord,MXCOOR,GRDARR,STPARR,& 
               & EGRAD,config%optinfo%STPSYM,lupri,config%optinfo)
          END IF
       !
       !     We go through the same procedure as for the first iteration,
       !     but here we call Get_Energy. We called the optimizer after 
       !     calculating the first energy in LSDALTON

!!!!! Vladimir:: No Newton yet!     
       !         IF (config%optinfo%Newton) THEN
       !            call ls_GTHESS(EGRAD,KEHESS,KALHES, &
       !     &           EXHER,EXSIR,EXABA,WORK(KWRK1),LWRK1,IWOFF, &
       !     &           WRKDLM,config%optinfo%NCoordTot,config%optinfo)
       !         ELSE
       !            call Get_Energy(E,config,config%optinfo,H1,F,D,S,ls,NAtoms,lupri,luerr)
          IF (.not.DECinfo%doDEC) &
       &  call Obtain_Gradient(E(1),Eerr,lupri,NAtoms,S,F(1),D(1),ls,config,CMO,config%optinfo)

       ENDIF  ! Simple scan

       !         END IF

       !
       !     Numeric Hessian in reduced space if asked
       !            Print *,'REDSPA',config%optinfo%RedSpa
       !            If (config%optinfo%RedSpa) then
       !               call Num_Int_Hess(ls,config,F,D,S,H1,config%optinfo, &
       !               &config%optinfo%Hess_dim,config%optinfo%NIntCoord,MXCOOR,&
       !               &NAtoms,KBMINV,WILBMT,lupri,luerr)
       !            Endif
       !        
       !
       !     If redundant internal coordinates are used, Wilson's B matrix,
       !     its derivative and its inverse must be updated.
       !
       IF (config%optinfo%RedInt .OR. config%optinfo%DelInt) THEN
          call ls_GETWIL(config%Molecule,config%optinfo%NIntCoord,MXCOOR, &
               &           KTEMP1,WILBMT, &
               &           BMTRAN,KTEMP2,lupri,NAtoms,config%optinfo)
          IF (config%optinfo%IPrint .GE. IPRMAX) &
               &           call ls_GETDWL(config%Molecule,config%optinfo%NIntCoord, &
               &                KTEMP1,KTEMP2, &
               &                KTEMP3,WILBMT,lupri,NAtoms,config%optinfo)
          call ls_GTBINV(config%optinfo%NIntCoord,KTEMP1,KTEMP2,KTEMP3, &
               &           KTEMP4,WILBMT,BMTRAN,KBMINV,KPJINM, &
               &           KTEMP5,KTEMP6,lupri,config%optinfo)
       END IF
       !
       IF (config%optinfo%RedInt .OR. config%optinfo%DelInt) THEN
          NCRDHS = config%optinfo%NIntCoord
       ELSE
          NCRDHS = config%optinfo%ICartCoord
       END IF
       IF (config%optinfo%RatFun) NCRDHS = NCRDHS + 1
!
! In case of simple scan we don't care about Hessians
!
If (.NOT. config%optinfo%simple_scan) then
!
       IF (.NOT. config%optinfo%Newton) THEN
          IF (config%optinfo%Rebid) THEN
             If (config%optinfo%Redint .OR. config%optinfo%Delint) then
                call ls_INIHES(config%Molecule,config%optinfo%NIntCoord,config%optinfo%NIntCoord, &
                     &              MXCOOR,GRDOLD,HESOLD, &
                     &              KTEMP1,KTEMP2,KTEMP3,KTEMP4, &
                     &              WILBMT,BMTRAN,KBMINV,HESINT,lupri,config%optinfo)
             Else
                call ls_INIHES(config%Molecule,MXCOOR,config%optinfo%NIntCoord, &
                     &              MXCOOR,GRDOLD,HESOLD, &
                     &              KTEMP1,KTEMP2,KTEMP3,KTEMP4, &
                     &              WILBMT,BMTRAN,KBMINV,HESINT,lupri,config%optinfo)
             Endif
             config%optinfo%Rebid = .FALSE.
          ELSE ! Update Hessian
             If (config%optinfo%RedInt .OR. config%optinfo%DelInt) then ! Internals
                call ls_UPDHES(config%Molecule,config%optinfo%NIntCoord,MXCOOR, &
                     &              GRDOLD,GRDMAT,STPMAT, &
                     &              HESOLD,KTEMP1,KTEMP2,KTEMP3,KTEMP4,KTEMP5,KTEMP6, &
                     &              KTEMP7,KTEMP8,KTEMP9,WILBMT,BMTRAN,KBMINV, &
                     &              HESINT,NINT(ABS(GEINFO(config%optinfo%ItrNmr-1,6))), &
                     &              NINT(ABS(GEINFO(config%optinfo%ItrNmr,6))),lupri,config%optinfo)
             Else   ! Cartesians
                call ls_UPDHES(config%Molecule,MXCOOR,MXCOOR, &
                     &              GRDOLD,GRDMAT,STPMAT, &
                     &              HESOLD,KTEMP1,KTEMP2,KTEMP3,KTEMP4,KTEMP5,KTEMP6, &
                     &              KTEMP7,KTEMP8,KTEMP9,WILBMT,BMTRAN,KBMINV, &
                     &              HESINT,NINT(ABS(GEINFO(config%optinfo%ItrNmr-1,6))), &
                     &              NINT(ABS(GEINFO(config%optinfo%ItrNmr,6))),lupri,config%optinfo)
             Endif
          END IF
          !
          call ls_DZERO(EGRAD,MXCOOR)
          call ls_DZERO(KEHESS,MXCOOR*MXCOOR)
          call ls_DZERO(KALHES,MXCOOR*MXCOOR)
          DO i = 1, config%optinfo%ICartCoord
             EGRAD(i) = config%optinfo%GradMol(i)
          ENDDO
          JI = 1
          DO I = 1, config%optinfo%ICartCoord
             DO J = 1, config%optinfo%ICartCoord
                KEHESS(JI) = config%optinfo%HessMol(J,I) 
                KALHES(JI) = config%optinfo%HessMol(J,I) 
                JI = JI + 1
             ENDDO
          ENDDO
       END IF
!
ENDIF  ! Simple scan
33     CONTINUE
!
! In case of simple scan we don't care about Hessians
!
If (.NOT. config%optinfo%simple_scan) then
!
       IF (config%optinfo%RedInt .OR. config%optinfo%DelInt) THEN
          IF (config%optinfo%Newton) call ls_CGHINT(config%Molecule, &
               &           config%optinfo%NIntCoord,MXCOOR,KTEMP1, &
               &           KTEMP2,KTEMP3,KTEMP4,KTEMP5, &
               &           WILBMT,KBMINV,BMTRAN,HESINT,lupri,config%optinfo)
          call ls_PRJINT(config%optinfo%NIntCoord,config%optinfo%NIntCoord,KPJINM,KCONMT, &
               &           HESINT,KTEMP1,KTEMP2,KTEMP3,KTEMP4,lupri,config%optinfo)
          IF (config%optinfo%LnSearch .AND. config%optinfo%RatFun .AND. (config%optinfo%ItrNmr .GT. 0)) &
               &           call ls_LINSRC(config%optinfo%NIntCoord,config%optinfo%NIntCoord,config%optinfo%GRDINT,GRDARR(1,1), &
               &           KTEMP7,STPARR(1,1),KTEMP3,KTEMP4, &
               &           ACTIVE,EMOD,lupri,config%optinfo)
          IF (config%optinfo%RatFun .AND. config%optinfo%Saddle) NCRDHS = NCRDHS - 1
          call ls_DIAINT(config%optinfo%NIntCoord,MXCOOR,NCRDHS,KEVEC,KTEMP1, &
               &           KTEMP2,KTEMP3,KTEMP4,THRIND,HESINT, &
               &           KTEMP5,lupri,config%optinfo)
          IF (config%optinfo%RatFun .AND. config%optinfo%Saddle) NCRDHS = NCRDHS + 1
       ELSE
          call ls_PROJGH(EGRAD,KEHESS,KALHES,KTEMP1, &
               &      KTEMP2,KTEMP3,KTEMP4,config%optinfo%NCoordTot, &
               &      config%optinfo%NTempMat,lupri,config%optinfo)
          IF (config%optinfo%LnSearch .AND. config%optinfo%RatFun .AND. (config%optinfo%ItrNmr .GT. 0)) &
               &           call ls_LINSRC(config%optinfo%ICartCoord,MXCOOR,EGRAD,GRDARR(1,1),CSTEP, &
               &           STPARR(1,1),KTEMP3,KTEMP4,ACTIVE,EMOD,lupri,config%optinfo)
          call ls_DIAHES(config%optinfo%NIntCoord,MXCOOR,NCRDHS,EGRAD,KEHESS, &
               &           KALHES,KTEMP1,THRIND,KEVEC,KTEMP2,KTEMP3,KTEMP4, &
               &           config%optinfo%NCoordTot,config%optinfo%NTempMat,lupri,config%optinfo)
       END IF

       !
       !     Update information for this iteration
       !
       GEINFO(config%optinfo%ItrNmr,3) = config%optinfo%IndHes*1.0E0_realk
       IREJ = 0
!!!!!! Vladimir:: currently disabled
       !     Write Hessian to file
       !
       !         IF (.NOT. config%optinfo%NoHessianWrite) &
       !     &      call ls_PNCHES(config%optinfo%NIntCoord,MXCOOR,HESINT,WILBMT,BMTRAN,KTEMP1, &
       !     &      KTEMP2,KTEMP3,KTEMP4,WORK(KWRK2),LWRK2, &
       !     &      lupri,config%optinfo)

ENDIF  ! Simple scan

       !
756    CONTINUE
       call ls_dzero(CSTEP,config%optinfo%ICartCoord)
       IF (config%optinfo%RedInt .OR. config%optinfo%DelInt) THEN
          call ls_FNSTIN(config%Molecule,config%optinfo%NIntCoord, &
               &           MXCOOR,NCRDHS,HESINT,KEVEC, &
               &           KTEMP1,KTEMP2,KTEMP3,KTEMP4, &
               &           KTEMP5,CSTEP,WILBMT,BMTRAN,KBMINV,GRDARR, &
               &           STPARR,ACTIVE,EMOD,VECMOD,KTEMP7,lupri,config%optinfo)
       ELSE
          call ls_FNDSTP(MXCOOR,MXCOOR,NCRDHS,EGRAD,KEHESS, &
               &           KEVEC,KTEMP1,KTEMP2,KTEMP3, &
               &           KTEMP4,KTEMP5,CSTEP,GRDARR,STPARR, &
               &           ACTIVE,EMOD,VECMOD,lupri,config%optinfo)
       END IF
!
! If we don't do a simple scan then we don't care about convergence
!
If (.NOT. config%optinfo%simple_scan) then 
       config%optinfo%GeConv = ls_minend(config%optinfo%NIntCoord,BMTRAN,KTEMP1,KTEMP2,lupri,config%optinfo)
       IF (ABS(GEINFO(config%optinfo%ItrNmr,6)) .GT. 1.0E-3_realk) config%optinfo%GeConv = .FALSE.
       call ls_PRIALL(Config%Molecule,NAtoms,CSTEP,KTEMP1,lupri,config%optinfo)
       !
       !
Else  ! simple scan
    config%optinfo%GeConv = .FALSE.
    call ls_PRIALL(Config%Molecule,NAtoms,CSTEP,KTEMP1,lupri,config%optinfo)
Endif
!
       IF (config%optinfo%RedInt .AND. (config%optinfo%IPrint .GE. 1)) THEN
          call ATOM_INI(KTEMP1,Config%Molecule,config%optinfo,NAtoms,.TRUE.,lupri)
          call ls_GETINT(NAtoms,config%optinfo%NIntCoord,KTEMP1,config%optinfo%CoordInt,lupri,config%optinfo)
          call lsheader(lupri,'New internal coordinates')
          call ls_output(config%optinfo%CoordInt,1,1,1,config%optinfo%NIntCoord,1,config%optinfo%NIntCoord,1,LUPRI)
          WRITE(LUPRI,'(//)')
       END IF
!
       NEWSTP = .FALSE.
       IF (.NOT. config%optinfo%GeConv) then 
          config%optinfo%dynamicChange = .FALSE.
          call Find_Geometry(E,CSTEP,EGRAD,KTEMP1,KTEMP2,&
               &        IREJ,GEINFO,NEWSTP,NEWBMT,NAtoms,lupri,luerr, &
               &        config,ls,H1,F,D,S,CMO,config%optinfo)
          IF (config%optinfo%dynamicChange) goto 300
          !  Grabbing coordinates
          Do i = 1, NAtoms
             config%optinfo%Coordinates(:,i) = KTEMP1(3*i-2:3*i)
          Enddo
          !    Renovating Config%Molecule
          Do i = 1,NAtoms
             Config%Molecule%Atom(i)%Center(:)=config%optinfo%Coordinates(:,i)
          Enddo
       ENDIF
       !            
       IF (NEWSTP) GOTO 756
       GEINFO(config%optinfo%ItrNmr,2) = config%optinfo%GradNorm
       GEINFO(config%optinfo%ItrNmr,4) = config%optinfo%StepNorm
       IF (config%optinfo%ItrNmr .LT. config%optinfo%MaxIter) GEINFO(config%optinfo%ItrNmr+1,5) = config%optinfo%TrustRad
       IF (ABS(GEINFO(config%optinfo%ItrNmr,6)) .LT. 1.0E-3_realk) THEN
          GEINFO(config%optinfo%ItrNmr,6) = IREJ*1.0E0_realk
       ELSE
          GEINFO(config%optinfo%ItrNmr,6) = -(ABS(GEINFO(config%optinfo%ItrNmr,6))+ABS(IREJ)*1.0E0_realk)
       END IF
       config%optinfo%TotRj = config%optinfo%TotRj + ABS(IREJ)
       IF (config%optinfo%Rebid) THEN
          call dcopy(config%optinfo%NIntCoord,config%optinfo%STPINT,1,KTEMP7,1)
          call ls_DZERO(config%optinfo%STPINT,config%optinfo%NIntCoord)
          call ls_DZERO(GRDOLD,config%optinfo%NIntCoord)
          DO I = 1, config%optinfo%NIntCoord
             DO J = 1, config%optinfo%NIntCoord
                config%optinfo%STPINT(I) = config%optinfo%STPINT(I) + BMTRAN(I,J)*KTEMP7(J-1)
                GRDOLD(I) = GRDOLD(I) + BMTRAN(I,J)*config%optinfo%GRDINT(J)
             ENDDO
          ENDDO
       END IF
       !
       !     If the step has failed
       !
       IF (IREJ .LT. 0) THEN
          IF (.NOT. config%optinfo%Newton) THEN
             IF (NEWBMT) THEN
                NCRDHS = config%optinfo%NIntCoord
                IF (config%optinfo%RatFun) NCRDHS = NCRDHS + 1
                call ls_GETWIL(config%Molecule,config%optinfo%NIntCoord, &
                     &                 MXCOOR,KTEMP1,WILBMT, &
                     &                 BMTRAN,KTEMP2,lupri,NAtoms,config%optinfo)
                call ls_GTBINV(config%optinfo%NIntCoord,KTEMP1,KTEMP2, &
                     &                 KTEMP3,KTEMP4,WILBMT,BMTRAN, &
                     &                 KBMINV,KPJINM,KTEMP5,KTEMP6,lupri,config%optinfo)
             END IF
             If (config%optinfo%Redint .OR. config%optinfo%Delint) then
                call ls_INIHES(config%Molecule,config%optinfo%NIntCoord,config%optinfo%NIntCoord, &
                     &              MXCOOR,GRDOLD,HESOLD, &
                     &              KTEMP1,KTEMP2,KTEMP3,KTEMP4, &
                     &              WILBMT,BMTRAN,KBMINV,HESINT,lupri,config%optinfo)
             Else
                call ls_INIHES(config%Molecule,MXCOOR,config%optinfo%NIntCoord, &
                     &              MXCOOR,GRDOLD,HESOLD, &
                     &              KTEMP1,KTEMP2,KTEMP3,KTEMP4, &
                     &              WILBMT,BMTRAN,KBMINV,HESINT,lupri,config%optinfo)
             Endif

             call ls_DZERO(EGRAD,MXCOOR)
             call ls_DZERO(KEHESS,MXCOOR*MXCOOR)
             call ls_DZERO(KALHES,MXCOOR*MXCOOR)
             DO i = 1, config%optinfo%ICartCoord
                EGRAD(i) = config%optinfo%GradMol(i)
             ENDDO
             JI = 1
             DO I = 1, config%optinfo%ICartCoord
                DO J = 1, config%optinfo%ICartCoord
                   KEHESS(JI) = config%optinfo%HessMol(J,I) 
                   KALHES(JI) = config%optinfo%HessMol(J,I) 
                   JI = JI + 1
                ENDDO
             ENDDO
             !      Call ls_copyGH(EGRAD,KEHESS,KALHES,config%optinfo%GradMol,config%optinfo%HessMol,config%optinfo%ICartCoord)
          END IF
          GOTO 33
          !         ELSE IF (config%optinfo%RejIni .AND. config%optinfo%RedInt .AND. (config%optinfo%TotRj .GE. 5)) THEN
          !            WRITE(LUPRI,*)'***** NOTE! *****'
          !            WRITE(LUPRI,*)
          !     &           'The number of dihedral angles will be reduced!'
          !            call ls_RREDUN
          !            config%optinfo%TotRj = 0
       END IF
       !
       !     Check if rejected steps should cause reinitialization of Hessian.
       !
       IF ((.NOT. config%optinfo%Newton) .AND. (config%optinfo%RejIni .AND. (IREJ .GE. 1))) THEN
          WRITE(LUPRI,*)
          WRITE(LUPRI,*)'***** NOTE! *****'
          WRITE(LUPRI,*) &
               &           'Due to rejected step, Hessian is reinitialized.'
          WRITE(LUPRI,*)
          If (config%optinfo%Redint .OR. config%optinfo%Delint) then
             call ls_INIHES(config%Molecule,config%optinfo%NIntCoord,config%optinfo%NIntCoord, &
                  &           MXCOOR,GRDOLD,HESOLD, &
                  &           KTEMP1,KTEMP2,KTEMP3,KTEMP4,WILBMT,BMTRAN, &
                  &           KBMINV,HESINT,lupri,config%optinfo)
          Else
             call ls_INIHES(config%Molecule,MXCOOR,config%optinfo%NIntCoord, &
                  &           MXCOOR,GRDOLD,HESOLD, &
                  &           KTEMP1,KTEMP2,KTEMP3,KTEMP4,WILBMT,BMTRAN, &
                  &           KBMINV,HESINT,lupri,config%optinfo)
          Endif
          config%optinfo%TrustRad = GEINFO(0,5)
          GEINFO(config%optinfo%ItrNmr+1,5) = config%optinfo%TrustRad
          config%optinfo%Restart = .TRUE.
       END IF
       !
       !     Check if increase of gradient norm should cause reinitialization
       !     of Hessian. Reinitialization occurs when the norm of the gradient
       !     is larger than the norm of the gradient two iterations earlier.
       !
       IF (.NOT.config%optinfo%Newton .AND. config%optinfo%GradIni .AND. (config%optinfo%ItrNmr .GE. 2)) THEN
          IF (GEINFO(config%optinfo%ItrNmr,2) .GE. GEINFO(config%optinfo%ItrNmr-2,2)) THEN
             WRITE(LUPRI,*)
             WRITE(LUPRI,*)'***** NOTE! *****'
             WRITE(LUPRI,*)'Due to increasing gradient norm,   &
                  &              Hessian is reinitialized.'
             WRITE(LUPRI,*)
             If (config%optinfo%Redint .OR. config%optinfo%Delint) then
                call ls_INIHES(config%Molecule,config%optinfo%NIntCoord,config%optinfo%NIntCoord, &
                     &              MXCOOR,GRDOLD,HESOLD, &
                     &              KTEMP1,KTEMP2,KTEMP3,KTEMP4,WILBMT,BMTRAN, &
                     &              KBMINV,HESINT,lupri,config%optinfo)
             Else
                call ls_INIHES(config%Molecule,MXCOOR,config%optinfo%NIntCoord, &
                     &              MXCOOR,GRDOLD,HESOLD, &
                     &              KTEMP1,KTEMP2,KTEMP3,KTEMP4,WILBMT,BMTRAN, &
                     &              KBMINV,HESINT,lupri,config%optinfo)
             Endif
             config%optinfo%TrustRad = GEINFO(0,5)
             GEINFO(config%optinfo%ItrNmr+1,5) = config%optinfo%TrustRad
          END IF
       END IF
       !
       IF (config%optinfo%GeConv .AND. config%optinfo%DoPre .AND. (.NOT. config%optinfo%FinPre)) THEN
          config%optinfo%KeepHessian = .FALSE.
       END IF
       !
       GOTO 10
       !
       !     Finished case 1: Geometry has converged.
       !
    ELSE IF (config%optinfo%GeConv) THEN
       !
       !     Final results are printed, partially through PRIINF.
       !

       config%optinfo%dynamicConvergence = .TRUE.
       call lsheader(lupri,' End of Optimization ')
       call ls_PRIINF(Config%Molecule,Config%Molecule%NAtoms,GEINFO,lupri,config%optinfo)
       WRITE(LUPRI,*)
       IF (config%optinfo%ConOpt) THEN
          WRITE(LUPRI,*) 'Constrained optimization converged in ', &
               &           config%optinfo%ItrNmr+1, ' iterations!'
          IF (config%optinfo%GradNorm .GT. config%optinfo%GradThr) THEN
             WRITE(LUPRI,*) 'Removing the  &
                  &            constraint(s) might decrease the energy further.'
          ELSE
             WRITE(LUPRI,*) 'A saddle point might have been reached.'
          END IF
       ELSE
          WRITE(LUPRI,*) 'Geometry converged in ', config%optinfo%ItrNmr+1, &
               &           ' iterations!'
       END IF
       IF (config%optinfo%Newton .AND. config%optinfo%Saddle .AND. (config%optinfo%IndHes .NE. 1)) THEN
          WRITE(LUPRI,'(/A/A)') &
               &         ' Please note that Hessian index does not correspond', &
               &         ' to a first order saddle point (transition state).'
       END IF
       E = GEINFO(config%optinfo%ItrNmr,1)
       WRITE(LUPRI,*)
       WRITE(LUPRI,'(A,F14.6,A)') &
            &          ' Energy at final geometry is       : ',E,' a.u.'
       ERGDIF = config%optinfo%energy - GEINFO (0,1)
       WRITE(LUPRI,'(A,F14.6,A)') &
            &          ' Energy change during optimization : ',ERGDIF,' a.u.'
       ERGDIF = ERGDIF * XKJMOL
       WRITE(LUPRI,'(A,F14.3,A)') &
            &        '                                     ',ERGDIF,' kJ/mol'
       WRITE(LUPRI,*)
       IF (config%optinfo%DoPre) THEN
          WRITE(LUPRI,'(A)') ' Preoptimization was performed using'// &
               &           ' the basis set(s):'
          DO I = 1, config%optinfo%Pre-1
             WRITE(LUPRI,'(A,A60)') '     ',config%optinfo%PreText(I)
          ENDDO
          WRITE(LUPRI,*)
       END IF
       IF (config%optinfo%DoSpE) THEN
          E = GEINFO(config%optinfo%ItrNmr+1,1)
          WRITE(LUPRI,'(A,A60)') ' Using the basis ',config%optinfo%SpBText
          WRITE(LUPRI,'(A,F14.6,A)') &
               &          ' single point energy was calculated: ',E,' a.u.'
          WRITE(LUPRI,*)
       END IF
       Do i = 1,NAtoms
          ls%input%Molecule%Atom(i)%Center(:)=Config%Molecule%Atom(i)%Center(:)
       Enddo
       CALL WRITE_MOLECULE_OUTPUT('MOLECULE.OUT',ls%input%Molecule,ls%input%basis,lupri)
       !
       !     Finished case 2: Exceeded maximum number of iterations.
       !
    ELSE
       !     No single point energy has been calculated.
       TMPLOG = config%optinfo%DoSpE
       config%optinfo%DoSpE = .FALSE.
       call lsheader(lupri,'Optimization Control Center')
       call ls_PRIINF(Config%Molecule,Config%Molecule%NAtoms,GEINFO,lupri,config%optinfo)
       config%optinfo%DoSpE = TMPLOG
       WRITE(LUPRI,*)
       WRITE(LUPRI,*) 'Geometry has NOT converged!'
       WRITE(LUPRI,*) 'Maximum number of iterations (', config%optinfo%MaxIter, &
            &                                        ') has been reached and'
       WRITE(LUPRI,*) 'optimization halted. Increase number or ', &
            &                                   'restart from last geometry.'
       IF (config%optinfo%DoSpE) WRITE(LUPRI,*) 'No single point energy has been  &
            &                          calculated.'
       WRITE(LUPRI,*)
    END IF
!
Endif ! Optimization
!
300 continue
    !
    !  Deallocate everything
    !      
    Call mem_dealloc(GEINFO)
    Call mem_dealloc(KEHESS)
    Call mem_dealloc(KALHES)
    Call mem_dealloc(CSTEP)
    Call mem_dealloc(KPJINM)
    Call mem_dealloc(KEVEC)
    Call mem_dealloc(KCONMT)
    Call mem_dealloc(KTEMP1)
    Call mem_dealloc(KTEMP2)
    Call mem_dealloc(KTEMP3)
    Call mem_dealloc(KTEMP4)
    Call mem_dealloc(KTEMP5)
    Call mem_dealloc(KTEMP6)
    Call mem_dealloc(KTEMP7)
    Call mem_dealloc(KTEMP8)
    Call mem_dealloc(KTEMP9)
    Call mem_dealloc(config%optinfo%Coordinates)
    Call mem_dealloc(config%optinfo%GradMol)
    Call mem_dealloc(config%optinfo%HessMol)
    Call mem_dealloc(config%optinfo%STPDIA)
    Call mem_dealloc(config%optinfo%STPSYM)
    Call mem_dealloc(config%optinfo%GRDDIA)
    Call mem_dealloc(config%optinfo%EVAL)
    Call mem_dealloc(config%optinfo%EVALOL)
    Call mem_dealloc(GRDMAT)
    Call mem_dealloc(STPMAT)
    Call mem_dealloc(GRDARR)
    Call mem_dealloc(EGRAD)
    Call mem_dealloc(GRDOLD)
    Call mem_dealloc(STPARR)
    Call mem_dealloc(HESOLD)
    If (config%optinfo%RedInt .OR. config%optinfo%DelInt .OR. &
     &  config%optinfo%InmdHess  .OR. config%optinfo%InrdHess) then 
       Call mem_dealloc(config%optinfo%GRDINT)
       Call mem_dealloc(config%optinfo%STPINT)
       Call mem_dealloc(config%optinfo%CoordInt)
       Call mem_dealloc(config%optinfo%INTCRD)
    Endif
    Call mem_dealloc(HESINT)
    Call mem_dealloc(BMTRAN)
    Call mem_dealloc(KBMINV)
    Call mem_dealloc(WILBMT)
    Call mem_dealloc(VECMOD)
    IF (config%optinfo%NFreeze.GT. 0) call mem_dealloc(config%optinfo%FreezeArray)
    IF (config%optinfo%NAdd.GT. 0) call mem_dealloc(config%optinfo%AddCoordArray)
    IF (config%optinfo%NumPre.GT. 0) call mem_dealloc(config%optinfo%PreText)
    IF (config%optinfo%RedSpa) call mem_dealloc(config%optinfo%Red_Space)
    IF (config%optinfo%CartRS) call mem_dealloc(config%optinfo%Cart_Red_Space)

    !
    RETURN
  END subroutine LS_RUNOP1

!=========================!
! Find_Geometry           !  
!=========================!

  SUBROUTINE Find_Geometry(E,CSTEP,EGRAD,COONEW,COOOLD, &
       &     IREJ,GEINFO,NEWSTP,NEWBMT,NAtoms,lupri,luerr, &
       config,ls,H1,F,D,S,C,optinfo)
    !
    !     If the step is acceptable, the geometry is updated
    !     and written to file.
    !
    implicit none
!    Implicit Real(realk) (A-H,O-Z)
    !
    Type(ConfigItem), intent(inout) :: Config ! General information
    Type(lsitem) :: ls   ! General information,used only to get E and gradient
    Type(Matrix), intent(inout) :: F(1),D(1),S ! Fock,density,overlap matrices
    Type(Matrix), intent(inout) :: H1   ! One electron matrix
    Type(Matrix), intent(inout) :: C       ! Orbitals
    !
    Integer :: lupri,luerr,NAtoms
    TYPE(opt_setting) :: optinfo
    Real(realk) :: CSTEP(MXCOOR), EGRAD(MXCOOR)
    Real(realk) :: COONEW(3,MXCENT), COOOLD(3,MXCENT)
    Integer     :: ICRD(3), IFAILD,IREJ,IJ,J,I,JJ
    Real(realk) :: GEINFO(0:optinfo%MaxIter+1,6)
    Real(realk) :: E(1)
    Real(realk) :: Eerr, Egeodiff, Eerrsave,graddi,fac
    CHARACTER*10 FILENM
    CHARACTER*12 molname
    LOGICAL REJGEO,NEWSTP,NEWBMT
    LOGICAL FAILED
    SAVE FAILED, IFAILD
    DATA FAILED, IFAILD /.FALSE.,0/
    !
    NEWSTP = .FALSE.
    REJGEO = .TRUE.
    Egeodiff = 0.0_realk
    Eerrsave=0.0_realk
    optinfo%energyOld = optinfo%energy
    IJ = 1
    DO J = 1, NAtoms
       DO I = 1, 3
          COOOLD(I,J) = optinfo%Coordinates(I,J)
          COONEW(I,J) = optinfo%Coordinates(I,J) + CSTEP(IJ)
          IJ = IJ + 1
       ENDDO
    ENDDO
    !
    !     Here we start a loop to obtain acceptable step, note that REJGEO
    !     is initially set TRUE to enter the loop.
    !

50  CONTINUE

    IF ((IREJ .LE. optinfo%MaxRej) .AND. REJGEO) THEN
       !
       !     If geometry stabilization has been requested (experimental feature!!!),
       !     all coordinates with a difference less than the limit are set equal.
       !
       IF (optinfo%Stblz .GT. 0) THEN
          IF (optinfo%IPrint .GT. 6) THEN
             call lsheader(lupri,'Non-stabilized geometry')
             Do i = 1,NAtoms
                Config%Molecule%Atom(i)%Center(:)=COONEW(:,i)
             Enddo
             call Print_Geometry(config%Molecule,lupri)
          END IF
          DO J = 1, nAtoms - 1
             ICRD(1) = NINT(COONEW(1,J)*10**optinfo%Stblz)
             ICRD(2) = NINT(COONEW(2,J)*10**optinfo%Stblz)
             ICRD(3) = NINT(COONEW(3,J)*10**optinfo%Stblz)
             DO JJ = J, nAtoms
                DO I = 1, 3
                   IF (ABS(ICRD(I) - NINT(COONEW(I,JJ)*10**optinfo%Stblz)) &
                        &                    .LE. 1) COONEW(I,JJ) = COONEW(I,J)
                ENDDO
             ENDDO
          ENDDO
          IF (optinfo%IPrint .GT. 6) THEN
             call lsheader(lupri,'Stabilized geometry')
             Do i = 1,NAtoms
                Config%Molecule%Atom(i)%Center(:)=COONEW(:,i)
             Enddo
             call Print_Geometry(config%Molecule,lupri)
          END IF
       END IF
       !
       IF (optinfo%IPrint .GT. 2) THEN
          call lsheader(lupri,'New geometry')
          Do i = 1,NAtoms
             Config%Molecule%Atom(i)%Center(:)=COONEW(:,i)
          Enddo
          call Print_Geometry(config%Molecule,lupri)
       END IF
       !
       !     Calculate energy at new geometry, which is compared to predicted
       !     energy(change) in UPTRAD.
       !     The temporary update of optinfo%ItrNmr is in case the molecule input is
       !     provided in the DALTON input file
       !

!!! Vladimir: disabled         optinfo%ItrNmr = optinfo%ItrNmr + 1
       !     Renovating geometry first
       optinfo%Coordinates = COONEW
       Call Update_coordinates(ls,config,optinfo,NAtoms)
       Call Get_Energy(E,Eerr,config,H1,F,D,S,ls,C,NAtoms,lupri,luerr)
       ! Do the force modification of energy if asked
       If (optinfo%FMPES) call FM_energy(E(1),optinfo)
       !
#ifdef VAR_DEC
       IF (DECinfo%dodec) call Obtain_Gradient(E(1),Eerr,lupri,NAtoms,S,F(1),D(1),ls,config,C,config%optinfo)
#else
       call lsquit('DEC requires -DVAR_DEC (-DENABLE_DEC=ON) ',-1)
#endif
       optinfo%energy = E(1)
       GEINFO(optinfo%ItrNmr+1,1) = E(1)
       IF (IREJ .EQ. 0) THEN
          call ls_UPTRAD(REJGEO,lupri,optinfo)
          IF (.NOT. REJGEO) THEN
             FAILED = .FALSE.
             IFAILD = 0
          END IF
          !
          !     After the first failure, we are satisfied if the new energy is below
          !     the last. No comparison with predicted energy is done.
          !
          ! Special case for dec and dynamical thresholding
          if( DECinfo%dodec) then
             write(lupri,'(a,i4,2g16.6,g22.12)') '1. DEC STAT: FOT level, Eerr, Ediff, E', &
                  & DECinfo%FOTlevel,Eerr,abs(optinfo%energyOld-E),E(1)
             if( Eerr.GT.abs(optinfo%energyOld-E(1)) ) then
                Egeodiff = abs(optinfo%energyOld-E(1))
                Eerrsave = Eerr
                optinfo%dynamicChange = .TRUE.
             end if
          end if
          ! Write molecule file for each iteration
          write(molname,'(A9,I3)') 'MOLECULE.',config%optinfo%ItrNmr+1
          DO I=10,11
            IF (molname(I:I).EQ.' ') molname(I:I)='0'
          ENDDO
          CALL WRITE_MOLECULE_OUTPUT(molname,ls%input%Molecule,ls%input%basis,lupri)
       ELSE IF ((optinfo%energy .GT. optinfo%energyOld)) THEN
          IF (FAILED .AND. (IFAILD .LE. optinfo%MaxRej) .AND. &
               &           (ABS(optinfo%energy-optinfo%energyOld) .LT. 1.0E-5_realk)) THEN
             WRITE(LUPRI,'(/A)') 'Trouble determining step,  &
                  &              accepting small energy increase.'
             IFAILD = IFAILD + 1
             REJGEO = .FALSE.
          ELSE
             WRITE(LUPRI,'(/A)') &
                  &              'Step rejected because energy is increasing.'
             WRITE(LUPRI,'(A,F10.5)')' Updated trust radius', optinfo%TrustRad
             REJGEO = .TRUE.
          END IF
       ELSE
          WRITE(LUPRI,'(/A)') 'Acceptable step has been found.'
          REJGEO = .FALSE.
          FAILED = .FALSE.
          ! Special case for dec and dynamical thresholding
          if( DECinfo%dodec) then
             write(lupri,'(a,i4,2g16.6,g22.12)') '2. DEC STAT: FOT level, Eerr, Ediff, E', &
                  & DECinfo%FOTlevel,Eerr,abs(optinfo%energyOld-E),E(1)
             if( Eerr.GT.abs(optinfo%energyOld-E(1)) ) then
                Egeodiff = abs(optinfo%energyOld-E(1))
                Eerrsave = Eerr
                optinfo%dynamicChange = .TRUE.
             end if
          end if

          ! Write molecule file for each iteration
          write(molname,'(A9,I3)') 'MOLECULE.',config%optinfo%ItrNmr
          DO I=10,11
            IF (molname(I:I).EQ.' ') molname(I:I)='0'
          ENDDO
          CALL WRITE_MOLECULE_OUTPUT(molname,ls%input%Molecule,ls%input%basis,lupri)
       END IF
       !
       IF (REJGEO) THEN

          IF(config%doESGopt)then
             !rejected step in excited state geometry optimization 
             !means we need to clean the cache in rsp module.
             call LS_rsp_eq_sol_empty!rsp_eq_sol in module rsp_equations
          endif

          IREJ = IREJ + 1
          !
          !     Line search based on quadratic model
          !
          GRADDI = 0.0E0_realk
          DO I = 1, optinfo%ICartCoord
             GRADDI = GRADDI + optinfo%GRDDIA(I)*optinfo%STPDIA(I)
          ENDDO
          GRADDI = GRADDI/optinfo%StepNorm
 

    ! Allocate gradient first

          !
          IF (optinfo%IPrint .GE. 12) THEN
             call lsheader(lupri,'Line search based on quadratic model')
             WRITE(LUPRI,'(A,F12.6)') &
                  &              ' Energy at last geometry     : ', optinfo%energyOld
             WRITE(LUPRI,'(A,F12.6)') &
                  &              ' Energy at rejected geometry : ', optinfo%energy
             WRITE(LUPRI,'(A,F12.6)') &
                  &              ' Norm of rejected step       : ', optinfo%StepNorm
             WRITE(LUPRI,'(A,F12.6)') &
                  &              ' Norm of gradient            : ', optinfo%GradNorm
             WRITE(LUPRI,'(A,F12.6)') &
                  &              ' Gradient along step         : ', GRADDI
          END IF
          !
          !     The minimum for a quadratic model is calculated with the formula
          !                    -f'(0)
          !     x     =  -------------------
          !      min     2*(f(1)-f(0)-f'(0))
          !
          FAC = -0.5E0_realk*GRADDI/(optinfo%energy-optinfo%energyOld-GRADDI)
          !
          !     If the factor found is very small or very large, we don't trust
          !     it. The factor is replaced by "safer" (but rather atbitrary) numbers.
          !
          IF (FAC .LT. 0.1E0_realk) FAC = 0.25E0_realk
          IF (FAC .GT. 0.9E0_realk) FAC = 0.75E0_realk
          !
          !     We have to update both steps and their norm.
          !
          DO I = 1, optinfo%NIntCoord
             optinfo%STPINT(I) = optinfo%STPINT(I)*FAC
          ENDDO
          DO I = 1, optinfo%ICartCoord
             optinfo%STPDIA(I) = optinfo%STPDIA(I)*FAC
             optinfo%STPSYM(I) = optinfo%STPSYM(I)*FAC
          ENDDO
          optinfo%StepNorm = optinfo%StepNorm*FAC
          !
          !     We also set the trust radius equal to the new norm
          !
          optinfo%TrustRad = optinfo%StepNorm
          !
          WRITE(LUPRI,'(A,F12.6)') &
               &           ' Minimum for quadratic model : ', FAC
          WRITE(LUPRI,'(A,F12.6)') &
               &           ' Norm of new step            : ', optinfo%StepNorm
          !
          !     Finally we construct a new geometry based on the factor found
          !
          DO J = 1, nAtoms
             DO I = 1, 3
                COONEW(I,J)=FAC*COONEW(I,J)+(1.0E0_realk-FAC)*COOOLD(I,J)
             ENDDO
          ENDDO
       END IF

       ! Special case for DEC, no more calculations if dynamic change is requested
       if(DECinfo%dodec .and. optinfo%dynamicChange) then
          GOTO 100
       else
          GOTO 50
       end if
    ELSE IF (REJGEO) THEN
       !
       !     Maximum number of allowed rejections reached
       !
       GEINFO(optinfo%ItrNmr,4) = optinfo%StepNorm
       IF (optinfo%ItrNmr .LT. optinfo%MaxIter) GEINFO(optinfo%ItrNmr+1,5) = optinfo%TrustRad
       GEINFO(optinfo%ItrNmr,6) = IREJ*1.0E0_realk
       !
       !     If redundant internal coordinates are used, we try reducing the
       !     number of dihedral angles to one third the original number (high
       !     redundancy might cause problems). We only allow this once before
       !     we give up (this should be viewed as an emergency solution!).
       !
       IF ((optinfo%RedInt .AND. (.NOT. FAILED)) .AND. (.NOT. optinfo%ConOpt)) THEN
          FAILED = .TRUE.
          IREJ = -IREJ
          GEINFO(optinfo%ItrNmr,6) = 0.0E0_realk
          WRITE(LUPRI,*) 'Maximum number of rejected steps (',optinfo%MaxRej, &
               &           ') reached.'
          WRITE(LUPRI,'(A)') 'No acceptable step found.'
          WRITE(LUPRI,'(/A)')'***** NOTE! *****'
          WRITE(LUPRI,'(A)')'As an emergency solution,  &
               &           the number of dihedral angles will be reduced!'
          call ls_RREDUN(lupri,optinfo)
          NEWBMT = .TRUE.
          optinfo%TrustRad = 0.5E0_realk
          RETURN
       ELSE IF (((.NOT. optinfo%Newton) .AND. (.NOT. FAILED)) .AND. &
            &           (.NOT. optinfo%ConOpt)) THEN
          FAILED = .TRUE.
          IREJ = -IREJ
          GEINFO(optinfo%ItrNmr,6) = 0.0E0_realk
          WRITE(LUPRI,*) 'Maximum number of rejected steps (',optinfo%MaxRej, &
               &           ') reached.'
          WRITE(LUPRI,'(A)') 'No acceptable step found.'
          WRITE(LUPRI,'(/A)')'***** NOTE! *****'
          WRITE(LUPRI,'(A)')'As a last resort,  &
               &           the Hessian is initialized to unity!'
          optinfo%EvLini = 1.0E0_realk
          optinfo%TrustRad = 0.5E0_realk
          RETURN
       ELSEIF(DECinfo%dodec) then
          ! For DEC we have the possibility to tighten FOT and start over
          goto 100

          !
          !     Otherwise we give up...
          !
       ELSE
          call ls_PRIINF(config%Molecule,config%Molecule%natoms,GEINFO,lupri,optinfo)
          WRITE(LUPRI,*) 'Maximum number of rejected steps (',optinfo%MaxRej, &
               &        ') reached.'
          WRITE(LUPRI,'(A)') 'No acceptable step found. Aborting.'
          call lsquit('*** FNDGEO *** No acceptable step found.',lupri)
       END IF
    END IF


100 CONTINUE
    IF (optinfo%dynamicChange) THEN
       write(lupri,*) 
       write(lupri,*) 'UNSTABLE DEC GEOMETRY CONVERGENCE - TRY TO TIGHTEN FOT!'
       write(lupri,*) '======================================================='
       write(lupri,'(1X,a,i5)')    'DECERR: Current FOT level              = ', DECinfo%FOTlevel
       if(DECinfo%FOTlevel<nFOTs) then
          write(lupri,'(1X,a,i5)') 'DECERR: FOT level will be increased to = ', DECinfo%FOTlevel+1
       else
         write(lupri,'(1X,a)') 'DECERR: Stopping geometry optimization because the intrinsic'
         write(lupri,'(1X,a)') 'DECERR: DEC error is too large!'
       end if
       write(lupri,'(1X,a,f20.10)') 'DECERR: Estimated intrinsic DEC error  = ', Eerrsave
       write(lupri,'(1X,a,f20.10)') 'DECERR: Energy diff between geometries = ', Egeodiff
       write(lupri,*) 
    ENDIF

  End subroutine Find_Geometry
!=================!
!  Num_Int_Hess   !
!=================!
  Subroutine Num_Int_Hess(ls,config,F,D,S,H1,C,optinfo,Red_Hess,Hess_dim,&
       NIntCoord,MaxCoor,NAtoms,B_Inv,B,lupri,luerr)
    Implicit none
    Integer :: Hess_dim, NIntCoord,NAtoms,MaxCoor
    Real(realk) :: Aux_Energy(1)
    Type(Opt_Setting) :: optinfo
    Real(realk) Red_Hess(Hess_dim,Hess_dim)
    Real(realk), pointer :: Aux_Int_Grad(:), Aux_Int_Grad_2(:)
    Real(realk), pointer :: Aux_Cart_Grad(:)
    Real(realk), pointer :: Aux_Cart_Coord(:,:)
    Real(realk) :: B_Inv(NIntCoord,MaxCoor)
    Real(realk) :: B(NIntCoord,MaxCoor)
    Type(lsitem),target :: ls   ! General information,used only to get E and gradient
    Type(Matrix), intent(inout) :: F(1),D(1),S     ! Fock,density,overlap matrices
    Type(Matrix), intent(inout),target :: H1 ! One electron matrix
    Type(Matrix), intent(inout) :: C       ! Orbitals
    Type(ConfigItem), intent(inout) :: Config ! General information
    Real(realk) :: delta  ! A finite difference in coordinate
    Real(realk) :: Eerr
    Integer :: i,j,k,l,lupri,luerr
    ! Defining delta
    delta = optinfo%displ
    ! Allocate some memory
    Call mem_alloc(Aux_Int_Grad,NIntCoord)
    If (optinfo%ForBac) then
       Call mem_alloc(Aux_Int_Grad_2,NIntCoord)
    Endif
    Call mem_alloc(Aux_Cart_Coord,3,NAtoms)
    Call mem_alloc(Aux_Cart_Grad,3*NAtoms)
    ! Save the geometry and gradient
    Aux_Cart_Coord = optinfo%Coordinates
    Aux_Cart_Grad = optinfo%GradMol
    ! Calculate Hessian
    Do i=1, Hess_dim
       ! Find the cartesian vector after displacement in internals
       Do j =1, NAtoms
          optinfo%Coordinates(:,j) = optinfo%Coordinates(:,j) + &
               & B_Inv(optinfo%Red_space(i),(3*j-2):3*j)*delta 
       Enddo
       ! Calculate energy and gradient
       Call Update_coordinates(ls,config,optinfo,NAtoms)
       Call Get_Energy(Aux_Energy,Eerr,config,H1,F,D,S,ls,C,NAtoms,lupri,luerr)
       Call obtain_Gradient(Aux_Energy(1),Eerr,lupri,NAtoms,S,F(1),D(1),ls,config,C,optinfo)
       ! Transform Cartesian gradient to internal coordinates
       Call LS_GX2GQ(NIntCoord,optinfo%GradMol,Aux_Int_Grad,B_Inv,optinfo)
       !
       ! Now the backward step if needed
       !
       If (optinfo%ForBac) then
          ! Restore coordinates 
          optinfo%Coordinates = Aux_Cart_Coord
          ! Find the cartesian vector after displacement in internal
          Do j =1, NAtoms
             optinfo%Coordinates(:,j) = optinfo%Coordinates(:,j) - &
                  & B_Inv(optinfo%Red_space(i),(3*j-2):3*j)*delta 
          Enddo
          ! Calculate energy and gradient
          Call Update_coordinates(ls,config,optinfo,NAtoms)
          Call Get_Energy(Aux_Energy,Eerr,config,H1,F,D,S,ls,C,NAtoms,lupri,luerr)
          Call Obtain_Gradient(Aux_Energy(1),Eerr,lupri,NAtoms,S,F(1),D(1),ls,config,C,optinfo)
          ! Transform Cartesian gradient to internal coordinates
          Call LS_GX2GQ(NIntCoord,optinfo%GradMol,Aux_Int_Grad_2,B_Inv,optinfo)
       Endif

       ! Forward finite difference
       If (.NOT. optinfo%ForBac) then
          ! Calculate Hessian in reduced space
          Do k = 1,Hess_dim
             Red_Hess(i,k) = &
                  (Aux_Int_Grad(optinfo%Red_space(k))-optinfo%GRDINT(optinfo%Red_space(k)))/delta
          Enddo
          ! Central
       Else
          ! Calculate Hessian in reduced space
          Do k = 1,Hess_dim
             Red_Hess(i,k) = &
                  &(Aux_Int_Grad(optinfo%Red_space(k))-Aux_Int_Grad_2(optinfo%Red_space(k)))/(2.00E0_realk*delta)
          Enddo
       Endif
       ! Restore coordinates 
       optinfo%Coordinates = Aux_Cart_Coord
    Enddo
    ! Symmetrize the internal reduced Hessian
    Do  i = 1, Hess_dim
       Do j = 1, i
          If (i .NE. j) then
             Red_Hess(i,j) = (Red_Hess(i,j)+Red_Hess(j,i))/2E0_realk
             Red_Hess(j,i)=Red_Hess(i,j)
          Endif
       Enddo
    Enddo
    ! Restore gradient
    optinfo%GradMol = Aux_Cart_Grad
    Print *, 'Reduced Hessian'
    Print *,  Red_Hess
    ! Deallocate memory
    Call mem_dealloc(Aux_Int_Grad)
    If (optinfo%ForBac) then
       Call mem_dealloc(Aux_Int_Grad_2)
    Endif
    Call mem_dealloc(Aux_Cart_Coord)
    Call mem_dealloc(Aux_Cart_Grad)
    !
  End subroutine Num_Int_Hess
!=================!
! Num_Cart_Hess   !
!=================!
  Subroutine Num_Cart_Hess(ls,config,F,D,S,H1,C,Red_Atoms,NAtoms,Red_Cart_Hess,optinfo,lupri,luerr)
    Implicit none
    Integer :: Red_Atoms, NAtoms ! Number of atoms in the reduced space and total number of atoms
    Real(realk) :: Aux_Energy(1)
    Type(Opt_Setting) :: optinfo
    Real(realk) Red_Cart_Hess(3*Red_Atoms,3*Red_Atoms)
    Type(lsitem),target :: ls   ! General information,used only to get E and gradient
    Type(Matrix), intent(inout) :: F(1),D(1),S     ! Fock,density,overlap matrices
    Type(Matrix), intent(inout),target :: H1 ! One electron matrix
    Type(Matrix), intent(inout) :: C       ! Orbitals
    Type(ConfigItem), intent(inout) :: Config ! General information
    Real(realk) :: Eerr
    Real(realk) :: delta  ! A finite difference in coordinate
    Real(realk), pointer :: Save_Cart_Grad(:)
    Real(realk), pointer :: Aux_Cart_Grad(:)
    Real(realk), pointer :: Aux_Cart_Coord(:,:)
    Integer :: i,j,k,l,m,n,lupri,luerr
    ! Defining delta
    delta = optinfo%displ
    ! Allocate some memory
    Call mem_alloc(Aux_Cart_Coord,3,NAtoms)
    Call mem_alloc(Save_Cart_Grad,3*NAtoms)
    If (optinfo%ForBac) Call mem_alloc(Aux_Cart_Grad,3*NAtoms)
    ! Save the geometry and gradient
    Aux_Cart_Coord = optinfo%Coordinates
    Save_Cart_Grad = optinfo%GradMol
    !
    ! Cartesian forward Hessian
    !
    Do i = 1, Red_Atoms
       l = optinfo%Cart_Red_Space(i)
       Do j = 1,3
          ! Shift coordinates
          optinfo%Coordinates(j,l) = optinfo%Coordinates(j,l) + delta
          ! Calculate energy and gradient
          Call Update_coordinates(ls,config,optinfo,NAtoms)
          Call Get_Energy(Aux_Energy,Eerr,config,H1,F,D,S,ls,C,NAtoms,lupri,luerr)
          Call Obtain_Gradient(Aux_Energy(1),Eerr,lupri,NAtoms,S,F(1),D(1),ls,config,C,optinfo)
          ! For central FD
          If (optinfo%ForBac) Aux_Cart_Grad = optinfo%GradMol 
          ! Central FD
          If (optinfo%ForBac) then
             ! Restore coordinates 
             optinfo%Coordinates = Aux_Cart_Coord
             ! Shift coordinates
             optinfo%Coordinates(j,l) = optinfo%Coordinates(j,l) - delta
             ! Calculate energy and gradient
             Call Update_coordinates(ls,config,optinfo,NAtoms)
             Call Get_Energy(Aux_Energy,Eerr,config,H1,F,D,S,ls,C,NAtoms,lupri,luerr)
             Call Obtain_Gradient(Aux_Energy(1),Eerr,lupri,NAtoms,S,F(1),D(1),ls,config,C,optinfo)
          Endif
          ! Forward FD Hessian
          If (.NOT. optinfo%ForBac) then
             Do k = 1, Red_Atoms
                n = optinfo%Cart_Red_Space(k)
                Do m = 1,3
                   Red_Cart_Hess(3*(i-1)+j,3*(k-1)+m)= &
                        & ( optinfo%GradMol(3*(n-1)+m)-Save_Cart_Grad(3*(n-1)+m) )/delta
                Enddo
             Enddo
             ! Central FD Hessian
          Else
             Do k = 1, Red_Atoms
                n = optinfo%Cart_Red_Space(k)
                Do m = 1,3
                   Red_Cart_Hess(3*(i-1)+j,3*(k-1)+m)= &
                        & ( Aux_Cart_Grad(3*(n-1)+m)-optinfo%GradMol(3*(n-1)+m) )/(2E0_realk*delta)
                Enddo
             Enddo
          Endif
          ! Restore coordinates 
          optinfo%Coordinates = Aux_Cart_Coord
       Enddo
    Enddo
    ! Restore gradient
    optinfo%GradMol = Save_Cart_Grad
    ! Symmetrize Hessian
    Do  i = 1, Red_Atoms*3
       Do j = 1, i
          If (i .NE. j) then
             Red_Cart_Hess(i,j) = (Red_Cart_Hess(i,j)+Red_Cart_Hess(j,i))/2E0_realk
             Red_Cart_Hess(j,i) = Red_Cart_Hess(i,j)
          Endif
       Enddo
    Enddo
    Write(*,*)'Reduced Cartesian Hessian'
    Write(*,*) Red_Cart_Hess(1:3,1:3)
    Write(*,*) Red_Cart_Hess(4:6,4:6)
    Write(*,*) Red_Cart_Hess(7:9,7:9)
    Write(*,*) Red_Cart_Hess(10:12,10:12)
    Write(*,*) Red_Cart_Hess(13:15,13:15)
    Write(*,*) Red_Cart_Hess
    ! Deallocate memory
    Call mem_dealloc(Aux_Cart_Coord)
    Call mem_dealloc(Save_Cart_Grad)
    If (optinfo%ForBac) Call mem_dealloc(Aux_Cart_Grad)
    !
  End subroutine Num_Cart_Hess
!==============!
! PES_scan     !
!==============!
! Scans potential energy surface along a bond
!
Subroutine PES_scan(ls,config,optinfo,H1,F,D,S,C,E,NAtoms,lupri,luerr)
Implicit none
Real(realk) :: E(1)
Type(opt_setting) :: optinfo
Integer :: Active_atom, Origin_atom,i,j,NAtoms,lupri,luerr
Real(realk), pointer :: Atom_array(:,:)
Type(lsitem),target :: ls   ! General information,used only to get E and gradient
Type(Matrix), intent(inout) :: F(1),D(1),S     ! Fock,density,overlap matrices
Type(Matrix), intent(inout),target :: H1 ! One electron matrix
    Type(Matrix), intent(inout) :: C       ! Orbitals
Type(ConfigItem), intent(inout) :: Config ! General information
    Real(realk) :: Eerr
! Allocate Atoms_array
Call mem_alloc(Atom_array,Natoms,8)
! Set active and origin atoms
Origin_atom = optinfo%Atoms_to_move(1,1)
Active_atom = optinfo%Atoms_to_move(2,1)
! Reference initial data
optinfo%Scan_info(1,1) = optinfo%coordinates(1,Active_atom)
optinfo%Scan_info(1,2) = E(1)
! Loop over all steps
Do i = 2, optinfo%MaxIter+1
   ! Displace the fragment connected to active atom
   Do j = 1, optinfo%N_to_move(2)
      optinfo%coordinates(1,optinfo%Atoms_to_move(2,j)) = &
      & optinfo%coordinates(1,optinfo%Atoms_to_move(2,j)) + optinfo%scan_step     
   Enddo
   ! Get energy
   Call Update_coordinates(ls,config,optinfo,NAtoms)
   Call Get_Energy(E,Eerr,config,H1,F,D,S,ls,C,NAtoms,lupri,luerr)
   ! Do the force modification of energy if asked
   If (optinfo%FMPES) call FM_energy(E(1),optinfo)
   ! Print Cartesian coordinates
   call lsheader(lupri,'New Cartesian coordinates') 
   call Print_Geometry(config%Molecule,lupri)
   ! Get and print new values of internal coordinates
   call ATOM_INI(Atom_array,config%Molecule,optinfo,NAtoms,.TRUE.,lupri)
   call ls_GETINT(NAtoms,optinfo%NIntCoord,Atom_array,optinfo%CoordInt,lupri,optinfo)
   call lsheader(lupri,'New internal coordinates') 
   call ls_output(optinfo%CoordInt,1,1,1,optinfo%NIntCoord,1,optinfo%NIntCoord,1,LUPRI)
   ! Reference the data
   optinfo%energy = E(1)
   optinfo%Scan_info(i,1) = optinfo%coordinates(1,Active_atom)
   optinfo%Scan_info(i,2) = E(1)
Enddo
! Deallocate Atoms_array
Call mem_dealloc(Atom_array)
!
End subroutine PES_scan
!======================!
! Update_coordinates   !
!======================!
! Copy new coordinates to ls and config
!
Subroutine Update_coordinates(ls,config,optinfo,NAtoms)
!
Implicit none
Integer :: NAtoms,i
Type(opt_setting) :: optinfo
Type(ConfigItem), intent(inout) :: Config ! General information
Type(lsitem),target :: ls   ! General information,used only to get E and gradient
       !
       Do i = 1,NAtoms
          ls%input%Molecule%Atom(i)%Center(:)=optinfo%Coordinates(:,i)
          Config%Molecule%Atom(i)%Center(:)=optinfo%Coordinates(:,i)
       Enddo
       !
End subroutine Update_coordinates
!===================!
! Obtain_gradient   !
!===================!
Subroutine Obtain_gradient(E,Eerr,lupri,NAtoms,S,F,D,ls,config,C,optinfo)
! A brief wrapper to call get_gradient
Implicit none
Type(opt_setting) :: optinfo
Integer :: NAtoms,lupri,i
Type(Matrix), intent(inout),target :: S  ! overlap matrices
Type(Matrix), intent(inout) :: F,D   ! Fock and density matrix
Type(Matrix), intent(inout) :: C       ! Orbitals
Type(lsitem) :: ls
Type(ConfigItem), intent(inout) :: Config ! General information
Real(realk), pointer :: Gradient(:,:)
Real(realk) :: h,Eerr,E ! Energy
Real(realk), pointer :: anaGrad(:,:)
Real(realk) :: direction(3),R_a(3),R_b(3)
logical    :: DEBUG_PAT
! Allocate gradient
Call mem_alloc(Gradient,3,NAtoms)

if( optinfo%doNumGradGeomOpt )then
   ! Calculate numerical gradient
   h = optinfo%findif_mesh
   call get_num_grad(h,lupri,config%luerr,ls,S,F,D,C,config,Gradient)
   DEBUG_PAT = .TRUE.
   IF (DEBUG_PAT) THEN
      Call mem_alloc(anaGrad,3,NAtoms)
      Call Get_Gradient(E,Eerr,lupri,NAtoms,S,F,D,ls,config,C,anaGrad)
      CALL LS_PRINT_GRADIENT(lupri,ls%setting%molecule(1)%p,anaGrad,nAtoms,'Ana grad')
      CALL LS_PRINT_GRADIENT(lupri,ls%setting%molecule(1)%p,Gradient,nAtoms,'Num grad')
      Do i = 1,NAtoms
            anaGrad(:,i) = Gradient(:,i) - anaGrad(:,i)
      Enddo
      write (*,*) "print difference AnaGradient - NumGradient"
      write (lupri,*) "print difference AnaGradient - NumGradient"
      CALL LS_PRINT_GRADIENT(lupri,ls%setting%molecule(1)%p,anaGrad,nAtoms,'Ana-Num grad')

      Call mem_dealloc(anaGrad)
   ENDIF
else
   ! Calculate analytical gradient
   Call Get_Gradient(E,Eerr,lupri,NAtoms,S,F,D,ls,config,C,Gradient)
endif


! Expand gradient to optinfo%GradMol
Do i = 1,NAtoms
   optinfo%GradMol(3*i-2:3*i) = Gradient(:,i)
Enddo
! Deallocate gradient
Call mem_dealloc(Gradient)
!
If (optinfo%FMPES) then
   ! Define the direction
   R_a = optinfo%Coordinates(:,optinfo%Att_atom(1))
   R_b = optinfo%Coordinates(:,optinfo%Att_atom(2))
   direction = (R_b - R_a)/(sqrt(dot_product(R_b-R_a,R_b-R_a)))
   ! Add external force
   optinfo%GradMol(optinfo%Att_atom(1)*3-2:optinfo%Att_atom(1)*3) = &
   optinfo%GradMol(optinfo%Att_atom(1)*3-2:optinfo%Att_atom(1)*3)+ &
   &  direction*optinfo%Ext_force
   optinfo%GradMol(optinfo%Att_atom(2)*3-2:optinfo%Att_atom(2)*3) = &
   optinfo%GradMol(optinfo%Att_atom(2)*3-2:optinfo%Att_atom(2)*3) - &
   &  direction*optinfo%Ext_force
Endif
!
end subroutine Obtain_gradient
!
end module LS_optimizer_mod




