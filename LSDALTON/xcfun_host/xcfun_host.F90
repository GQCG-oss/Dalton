module xcfun_host
  use precision
#ifdef VAR_XCFUN
  use xcfun
#endif
  implicit none
  integer :: XCFunFunctional
  logical :: USEXCFUN
  public :: xcfun_host_init, xcfun_host_free, USEXCFUN, &
       & xcfun_gga_xc_single_eval, xcfun2_gga_xc_single_eval, &
       & xcfun_lda_xc_single_eval, xcfun_gga_unres_xc_single_eval, &
       & xcfun_lda_unres_xc_single_eval, xcfun3_gga_xc_single_eval,&
       & xcfun_meta_xc_single_eval, xcfun2_lda_xc_single_eval
  private
  contains
  subroutine xcfun_host_init(DFTfuncString,hfweight,lupri)
    integer :: lupri
    character(len=80),intent(in)  :: DFTfuncString
    real(realk),intent(out) :: hfweight
    !
    character(len=80),pointer  :: DFTfuncStringSingle(:)
    integer :: Ipos,nStrings,ierr,I,ierrLDA,ierrGGA,ierrMETA,ierrHF
    logical :: GGAkeyString
    real(realk),pointer :: WeightSingle(:)
    GGAkeyString = .FALSE.
    IPOS = INDEX(DFTfuncString,'GGAKEY')
    IF(IPOS.NE.0)GGAkeyString = .TRUE.
    IPOS = INDEX(DFTfuncString,'GGAkey')
    IF(IPOS.NE.0)GGAkeyString = .TRUE.
    IPOS = INDEX(DFTfuncString,'GGAKey')
    IF(IPOS.NE.0)GGAkeyString = .TRUE.

    IF(GGAkeyString)THEN
       IPOS = INDEX(DFTfuncString,'GGAKey')
       call determine_nStrings(DFTfuncString(IPOS+7:80),nStrings)
       allocate(DFTfuncStringSingle(nStrings))
       allocate(WeightSingle(nStrings))
       call trim_strings(DFTfuncString(IPOS+7:80),nStrings,&
            & DFTfuncStringSingle,WeightSingle)
    ELSE
       nStrings = 1
       allocate(DFTfuncStringSingle(nStrings))
       allocate(WeightSingle(nStrings))
       call trim_strings(DFTfuncString,nStrings,&
            & DFTfuncStringSingle,WeightSingle)
       WeightSingle(1) = 1.0E0_realk
    ENDIF
#ifdef VAR_XCFUN
    XCFUNfunctional = xc_new_functional()
    do I=1,nStrings
       ierr = xc_set(XCFUNfunctional,DFTfuncStringSingle(I),WeightSingle(I))
       IF(ierr.NE.0)THEN
          print*,'The functional name:',DFTfuncStringSingle(I),'was not recognized'
          print*,' by the XCFUN program '
          write(lupri,*)'The functional name:',DFTfuncStringSingle(I),'was not recognized'
          write(lupri,*)' by the XCFUN program '
          call lsquit('xcfun_host_init: error',lupri)
       ENDIF
    enddo   
#else
    call lsquit('xcfun not activated -DVAR_XCFUN (can only be done using cmake)',-1)
#endif
    !Test for the different types of functional (LDA,GGA,META)
    !we start assuming that we only need the parameters need for LDA:
    ! \frac{\partial f}{\partial \rho} 
    !if xc_eval_setup gives an error then this assumption is wrong and 
    !we try to assume that we only need the parameters need for GGA:
    ! \frac{\partial f}{\partial \rho} 
    ! \frac{\partial f}{\partial |\nabla \rho|^{2}}
    !if xc_eval_setup gives an error then this assumption is wrong and 
    !it must be a META GGA
    ! \frac{\partial f}{\partial \rho} 
    ! \frac{\partial f}{\partial |\nabla \rho|^{2}}
    ! \frac{\partial f}{\partial \tau}
    !Test if this is a LDA type
#ifdef VAR_XCFUN
    ierrLDA = xc_eval_setup(XCFUNfunctional,XC_N,XC_PARTIAL_DERIVATIVES,1)
    ierrLDA=2
    IF(ierrLDA.NE.0)THEN
       !Test if this is a GGA type
       ierrGGA = xc_eval_setup(XCFUNfunctional,XC_N_GNN,XC_PARTIAL_DERIVATIVES,1)
       IF(ierrGGA.NE.0)THEN
       !Test if this is a Meta type
          ierrMETA = xc_eval_setup(XCFUNfunctional,XC_N_GNN_TAUN,XC_PARTIAL_DERIVATIVES,1)
          IF(ierrMETA.NE.0)THEN
             call lsquit('Error determining the correct XC type.',lupri)
          ELSE
             WRITE(lupri,*)'The Functional chosen is a META type functional'
          ENDIF
       ELSE
          WRITE(lupri,*)'The Functional chosen is a GGA type functional'
       ENDIF
    ELSE
       WRITE(lupri,*)'The Functional chosen is a LDA type functional'
    ENDIF
    ierrHF = xc_get(XCFUNfunctional,'EXX',hfweight)
#else
    call lsquit('xcfun not activated -DVAR_XCFUN (can only be done using cmake)',-1)
#endif
    IF(ierrHF.NE.0)THEN
       call lsquit('Error determining exact (HF) exchange weight.',lupri)
    ENDIF
    IF(ABS(hfweight).GT.1.0E-16_realk)THEN
       WRITE(lupri,*)'The Functional chosen contains a exact exchange contribution'
       WRITE(lupri,*)'with the weight:', hfweight
    ELSE
       WRITE(lupri,*)'The Functional chosen contains no exact exchange contribution'
    ENDIF
    deallocate(DFTfuncStringSingle)
    deallocate(WeightSingle)

  end subroutine xcfun_host_init


  subroutine xcfun3_gga_components_xc_single_eval(XCFUNINPUT,XCFUNOUTPUT)
    implicit none
    REAL(REALK),intent(in) :: XCFUNINPUT(2,1)
    REAL(REALK),intent(inout) :: XCFUNOUTPUT(35,1)
    integer :: ierr
#ifdef VAR_XCFUN
    !rho    = XCFUNINPUT(1,1) 
    !grad_x = XCFUNINPUT(2,1) 
    !grad_y = XCFUNINPUT(3,1) 
    !grad_z = XCFUNINPUT(4,1) 
    ierr = xc_eval_setup(XCFUNfunctional,XC_N_NX_NY_NZ,XC_PARTIAL_DERIVATIVES,3)
    IF(ierr.NE.0) then
       print*,'ierr from xcfun',ierr
       call lsquit('Unexpected error from xcfun',-1)
    endif
    call xc_eval(XCFUNfunctional,1,XCFUNINPUT,XCFUNOUTPUT)
    ! XCFUNOUTPUT(1,1) - Exc
    ! XCFUNOUTPUT(2,1) - d Exc/d rho
    ! XCFUNOUTPUT(3,1) - d Exc/d gx
    ! XCFUNOUTPUT(4,1) - d Exc/d gy
    ! XCFUNOUTPUT(5,1) - d Exc/d gz

    ! XCFUNOUTPUT(6,1) - d^2 Exc/d rho d rho
    ! XCFUNOUTPUT(7,1) - d^2 Exc/d rho d gx
    ! XCFUNOUTPUT(8,1) - d^2 Exc/d rho d gy
    ! XCFUNOUTPUT(9,1) - d^2 Exc/d rho d gz
    ! XCFUNOUTPUT(10,1) - d^2 Exc/d gx d gx
    ! XCFUNOUTPUT(11,1) - d^2 Exc/d gx d gy
    ! XCFUNOUTPUT(12,1) - d^2 Exc/d gx d gz
    ! XCFUNOUTPUT(13,1) - d^2 Exc/d gy d gy
    ! XCFUNOUTPUT(14,1) - d^2 Exc/d gy d gz
    ! XCFUNOUTPUT(15,1) - d^2 Exc/d gz d gz

    ! XCFUNOUTPUT(16,1) - d^3 Exc/d rho d rho d rho
    ! XCFUNOUTPUT(17,1) - d^3 Exc/d rho d rho d gx
    ! XCFUNOUTPUT(18,1) - d^3 Exc/d rho d rho d gy
    ! XCFUNOUTPUT(19,1) - d^3 Exc/d rho d rho d gz
    ! XCFUNOUTPUT(20,1) - d^3 Exc/d rho d gx d gx
    ! XCFUNOUTPUT(21,1) - d^3 Exc/d rho d gx d gy
    ! XCFUNOUTPUT(22,1) - d^3 Exc/d rho d gx d gz
    ! XCFUNOUTPUT(23,1) - d^3 Exc/d rho d gy d gy
    ! XCFUNOUTPUT(24,1) - d^3 Exc/d rho d gy d gz
    ! XCFUNOUTPUT(25,1) - d^3 Exc/d rho d gz d gz
    ! XCFUNOUTPUT(26,1) - d^3 Exc/d gx d gx d gx
    ! XCFUNOUTPUT(27,1) - d^3 Exc/d gx d gx d gy
    ! XCFUNOUTPUT(28,1) - d^3 Exc/d gx d gx d gz
    ! XCFUNOUTPUT(29,1) - d^3 Exc/d gx d gy d gy
    ! XCFUNOUTPUT(30,1) - d^3 Exc/d gx d gy d gz
    ! XCFUNOUTPUT(31,1) - d^3 Exc/d gx d gz d gz
    ! XCFUNOUTPUT(32,1) - d^3 Exc/d gy d gy d gy
    ! XCFUNOUTPUT(33,1) - d^3 Exc/d gy d gy d gz
    ! XCFUNOUTPUT(34,1) - d^3 Exc/d gy d gz d gz
    ! XCFUNOUTPUT(35,1) - d^3 Exc/d gz d gz d gz
#else
    call lsquit('xcfun not activated -DVAR_XCFUN (can only be done using cmake)',-1)
#endif
  end subroutine xcfun3_gga_components_xc_single_eval


  subroutine xcfun2_gga_components_xc_single_eval(XCFUNINPUT,XCFUNOUTPUT)
    implicit none
    REAL(REALK),intent(in) :: XCFUNINPUT(2,1)
    REAL(REALK),intent(inout) :: XCFUNOUTPUT(15,1)
    integer :: ierr
#ifdef VAR_XCFUN
    !rho    = XCFUNINPUT(1,1) 
    !grad_x = XCFUNINPUT(2,1) 
    !grad_y = XCFUNINPUT(3,1) 
    !grad_z = XCFUNINPUT(4,1) 
    ierr = xc_eval_setup(XCFUNfunctional,XC_N_NX_NY_NZ,XC_PARTIAL_DERIVATIVES,3)
    IF(ierr.NE.0) then
       print*,'ierr from xcfun',ierr
       call lsquit('Unexpected error from xcfun',-1)
    endif
    call xc_eval(XCFUNfunctional,1,XCFUNINPUT,XCFUNOUTPUT)
    ! XCFUNOUTPUT(1,1) - Exc
    ! XCFUNOUTPUT(2,1) - d Exc/d rho
    ! XCFUNOUTPUT(3,1) - d Exc/d gx
    ! XCFUNOUTPUT(4,1) - d Exc/d gy
    ! XCFUNOUTPUT(5,1) - d Exc/d gz

    ! XCFUNOUTPUT(6,1) - d^2 Exc/d rho d rho
    ! XCFUNOUTPUT(7,1) - d^2 Exc/d rho d gx
    ! XCFUNOUTPUT(8,1) - d^2 Exc/d rho d gy
    ! XCFUNOUTPUT(9,1) - d^2 Exc/d rho d gz
    ! XCFUNOUTPUT(10,1) - d^2 Exc/d gx d gx
    ! XCFUNOUTPUT(11,1) - d^2 Exc/d gx d gy
    ! XCFUNOUTPUT(12,1) - d^2 Exc/d gx d gz
    ! XCFUNOUTPUT(13,1) - d^2 Exc/d gy d gy
    ! XCFUNOUTPUT(14,1) - d^2 Exc/d gy d gz
    ! XCFUNOUTPUT(15,1) - d^2 Exc/d gz d gz
#else
    call lsquit('xcfun not activated -DVAR_XCFUN (can only be done using cmake)',-1)
#endif
  end subroutine xcfun2_gga_components_xc_single_eval

  subroutine xcfun1_gga_components_xc_single_eval(XCFUNINPUT,XCFUNOUTPUT)
    implicit none
    REAL(REALK),intent(in) :: XCFUNINPUT(2,1)
    REAL(REALK),intent(inout) :: XCFUNOUTPUT(5,1)
    integer :: ierr
#ifdef VAR_XCFUN
    !rho    = XCFUNINPUT(1,1) 
    !grad_x = XCFUNINPUT(2,1) 
    !grad_y = XCFUNINPUT(3,1) 
    !grad_z = XCFUNINPUT(4,1) 
    ierr = xc_eval_setup(XCFUNfunctional,XC_N_NX_NY_NZ,XC_PARTIAL_DERIVATIVES,3)
    IF(ierr.NE.0) then
       print*,'ierr from xcfun',ierr
       call lsquit('Unexpected error from xcfun',-1)
    endif
    call xc_eval(XCFUNfunctional,1,XCFUNINPUT,XCFUNOUTPUT)
    ! XCFUNOUTPUT(1,1) - Exc
    ! XCFUNOUTPUT(2,1) - d Exc/d rho
    ! XCFUNOUTPUT(3,1) - d Exc/d gx
    ! XCFUNOUTPUT(4,1) - d Exc/d gy
    ! XCFUNOUTPUT(5,1) - d Exc/d gz
#else
    call lsquit('xcfun not activated -DVAR_XCFUN (can only be done using cmake)',-1)
#endif
  end subroutine xcfun1_gga_components_xc_single_eval

  subroutine xcfun2_gga_unres_components_xc_single_eval(XCFUNINPUT,XCFUNOUTPUT)
    implicit none
    REAL(REALK),intent(in) :: XCFUNINPUT(2,1)
    REAL(REALK),intent(inout) :: XCFUNOUTPUT(45,1)
    integer :: ierr
#ifdef VAR_XCFUN
    !rho_a   = XCFUNINPUT(1,1) 
    !rho_b   = XCFUNINPUT(2,1) 
    !grada_x = XCFUNINPUT(3,1) 
    !grada_y = XCFUNINPUT(4,1) 
    !grada_z = XCFUNINPUT(5,1) 
    !gradb_x = XCFUNINPUT(6,1) 
    !gradb_y = XCFUNINPUT(7,1) 
    !gradb_z = XCFUNINPUT(8,1) 
    ierr = xc_eval_setup(XCFUNfunctional,XC_N_NX_NY_NZ,XC_PARTIAL_DERIVATIVES,3)
    IF(ierr.NE.0) then
       print*,'ierr from xcfun',ierr
       call lsquit('Unexpected error from xcfun',-1)
    endif
    call xc_eval(XCFUNfunctional,1,XCFUNINPUT,XCFUNOUTPUT)
    ! XCFUNOUTPUT(1,1) - Exc

    ! XCFUNOUTPUT(2,1) - d Exc/d rho_a
    ! XCFUNOUTPUT(3,1) - d Exc/d rho_b
    ! XCFUNOUTPUT(4,1) - d Exc/d gax
    ! XCFUNOUTPUT(5,1) - d Exc/d gay
    ! XCFUNOUTPUT(6,1) - d Exc/d gaz
    ! XCFUNOUTPUT(7,1) - d Exc/d gbx
    ! XCFUNOUTPUT(8,1) - d Exc/d gby
    ! XCFUNOUTPUT(9,1) - d Exc/d gbz

    ! XCFUNOUTPUT(10,1) - d^2 Exc/d rho_a d rho_a
    ! XCFUNOUTPUT(11,1) - d^2 Exc/d rho_a d rho_b
    ! XCFUNOUTPUT(12,1) - d^2 Exc/d rho_a d gax
    ! XCFUNOUTPUT(13,1) - d^2 Exc/d rho_a d gay
    ! XCFUNOUTPUT(14,1) - d^2 Exc/d rho_a d gaz
    ! XCFUNOUTPUT(15,1) - d^2 Exc/d rho_a d gbx
    ! XCFUNOUTPUT(16,1) - d^2 Exc/d rho_a d gby
    ! XCFUNOUTPUT(17,1) - d^2 Exc/d rho_a d gbz

    ! XCFUNOUTPUT(18,1) - d^2 Exc/d rho_b d rho_b
    ! XCFUNOUTPUT(19,1) - d^2 Exc/d rho_b d gax
    ! XCFUNOUTPUT(20,1) - d^2 Exc/d rho_b d gay
    ! XCFUNOUTPUT(21,1) - d^2 Exc/d rho_b d gaz
    ! XCFUNOUTPUT(22,1) - d^2 Exc/d rho_b d gbx
    ! XCFUNOUTPUT(23,1) - d^2 Exc/d rho_b d gby
    ! XCFUNOUTPUT(24,1) - d^2 Exc/d rho_b d gbz

    ! XCFUNOUTPUT(25,1) - d^2 Exc/d gax d gax
    ! XCFUNOUTPUT(26,1) - d^2 Exc/d gax d gay
    ! XCFUNOUTPUT(27,1) - d^2 Exc/d gax d gaz
    ! XCFUNOUTPUT(28,1) - d^2 Exc/d gax d gbx
    ! XCFUNOUTPUT(29,1) - d^2 Exc/d gax d gby
    ! XCFUNOUTPUT(30,1) - d^2 Exc/d gax d gbz

    ! XCFUNOUTPUT(31,1) - d^2 Exc/d gay d gay
    ! XCFUNOUTPUT(32,1) - d^2 Exc/d gay d gaz
    ! XCFUNOUTPUT(33,1) - d^2 Exc/d gay d gbx
    ! XCFUNOUTPUT(34,1) - d^2 Exc/d gay d gby
    ! XCFUNOUTPUT(35,1) - d^2 Exc/d gay d gbz

    ! XCFUNOUTPUT(36,1) - d^2 Exc/d gaz d gaz
    ! XCFUNOUTPUT(37,1) - d^2 Exc/d gaz d gbx
    ! XCFUNOUTPUT(38,1) - d^2 Exc/d gaz d gby
    ! XCFUNOUTPUT(39,1) - d^2 Exc/d gaz d gbz

    ! XCFUNOUTPUT(40,1) - d^2 Exc/d gbx d gbx
    ! XCFUNOUTPUT(41,1) - d^2 Exc/d gbx d gby
    ! XCFUNOUTPUT(42,1) - d^2 Exc/d gbx d gbz

    ! XCFUNOUTPUT(43,1) - d^2 Exc/d gby d gby
    ! XCFUNOUTPUT(44,1) - d^2 Exc/d gby d gbz

    ! XCFUNOUTPUT(45,1) - d^2 Exc/d gbz d gbz
#else
    call lsquit('xcfun not activated -DVAR_XCFUN (can only be done using cmake)',-1)
#endif
  end subroutine xcfun2_gga_unres_components_xc_single_eval


  subroutine xcfun1_gga_unres_components_xc_single_eval(XCFUNINPUT,XCFUNOUTPUT)
    implicit none
    REAL(REALK),intent(in) :: XCFUNINPUT(2,1)
    REAL(REALK),intent(inout) :: XCFUNOUTPUT(9,1)
    integer :: ierr
#ifdef VAR_XCFUN
    !rho_a   = XCFUNINPUT(1,1) 
    !rho_b   = XCFUNINPUT(2,1) 
    !grada_x = XCFUNINPUT(3,1) 
    !grada_y = XCFUNINPUT(4,1) 
    !grada_z = XCFUNINPUT(5,1) 
    !gradb_x = XCFUNINPUT(6,1) 
    !gradb_y = XCFUNINPUT(7,1) 
    !gradb_z = XCFUNINPUT(8,1) 
    ierr = xc_eval_setup(XCFUNfunctional,XC_N_NX_NY_NZ,XC_PARTIAL_DERIVATIVES,3)
    IF(ierr.NE.0) then
       print*,'ierr from xcfun',ierr
       call lsquit('Unexpected error from xcfun',-1)
    endif
    call xc_eval(XCFUNfunctional,1,XCFUNINPUT,XCFUNOUTPUT)
    ! XCFUNOUTPUT(1,1) - Exc

    ! XCFUNOUTPUT(2,1) - d Exc/d rho_a
    ! XCFUNOUTPUT(3,1) - d Exc/d rho_b
    ! XCFUNOUTPUT(4,1) - d Exc/d gax
    ! XCFUNOUTPUT(5,1) - d Exc/d gay
    ! XCFUNOUTPUT(6,1) - d Exc/d gaz
    ! XCFUNOUTPUT(7,1) - d Exc/d gbx
    ! XCFUNOUTPUT(8,1) - d Exc/d gby
    ! XCFUNOUTPUT(9,1) - d Exc/d gbz
#else
    call lsquit('xcfun not activated -DVAR_XCFUN (can only be done using cmake)',-1)
#endif
  end subroutine xcfun1_gga_unres_components_xc_single_eval


  subroutine xcfun_gga_xc_single_eval(XCFUNINPUT,XCFUNOUTPUT)
    implicit none
    REAL(REALK),intent(in) :: XCFUNINPUT(2,1)
    REAL(REALK),intent(inout) :: XCFUNOUTPUT(3,1)
#ifdef VAR_XCFUN
    !rho = XCFUNINPUT(1,1) 
    !gg = |grad|^2 = XCFUNINPUT(2,1) 
    call xc_eval(XCFUNfunctional,1,XCFUNINPUT,XCFUNOUTPUT)
    ! XCFUNOUTPUT(1,1) - Exc
    ! XCFUNOUTPUT(2,1) - d Exc/d rho
    ! XCFUNOUTPUT(3,1) - d Exc/d gg
#else
    call lsquit('xcfun not activated -DVAR_XCFUN (can only be done using cmake)',-1)
#endif
  end subroutine xcfun_gga_xc_single_eval



  subroutine xcfun_meta_xc_single_eval(XCFUNINPUT,XCFUNOUTPUT)
    implicit none
    REAL(REALK),intent(in) :: XCFUNINPUT(3,1)
    REAL(REALK),intent(inout) :: XCFUNOUTPUT(4,1)
#ifdef VAR_XCFUN
    !rho = XCFUNINPUT(1,1)   = \sum_{\mu \nu} D_{\mu \nu} chi_{mu} \times chi_{nu}
    !gg = |grad|^2 = XCFUNINPUT(2,1) = \sum_{\mu \nu} D_{\mu \nu} (\nabla chi_{mu} \times chi_{nu} + chi_{mu} \times \nabla chi_{nu})
    !tau = XCFUNINPUT(3,1)  = \sum_{\mu \nu} D_{\mu \nu} half nabla chi_{mu} \times \nabla chi_{nu}
    call xc_eval(XCFUNfunctional,1,XCFUNINPUT,XCFUNOUTPUT)
    ! XCFUNOUTPUT(1,1) - Exc
    ! XCFUNOUTPUT(2,1) - d Exc/d rho
    ! XCFUNOUTPUT(3,1) - d Exc/d gg
    ! XCFUNOUTPUT(4,1) - d Exc/d tau
#else
    call lsquit('xcfun not activated -DVAR_XCFUN (can only be done using cmake)',-1)
#endif
  end subroutine xcfun_meta_xc_single_eval

  subroutine xcfun2_gga_xc_single_eval(XCFUNINPUT,XCFUNOUTPUT)
    implicit none
    REAL(REALK),intent(in) :: XCFUNINPUT(2,1)
    REAL(REALK),intent(inout) :: XCFUNOUTPUT(6,1)
    !
    integer :: ierrGGA
#ifdef VAR_XCFUN
    !rho = XCFUNINPUT(1,1) 
    !gg = |grad|^2 = XCFUNINPUT(2,1) 
    ierrGGA = xc_eval_setup(XCFUNfunctional,XC_N_GNN,XC_PARTIAL_DERIVATIVES,2)
    call xc_eval(XCFUNfunctional,1,XCFUNINPUT,XCFUNOUTPUT)
    ierrGGA = xc_eval_setup(XCFUNfunctional,XC_N_GNN,XC_PARTIAL_DERIVATIVES,1)
    ! XCFUNOUTPUT(1,1) - Exc
    ! XCFUNOUTPUT(2,1) - d Exc/d rho
    ! XCFUNOUTPUT(3,1) - d Exc/d gg
    ! XCFUNOUTPUT(4,1) - d^2 Exc/d rho^2
    ! XCFUNOUTPUT(5,1) - d^2 Exc/d rho d gg
    ! XCFUNOUTPUT(6,1) - d^2 Exc/d gg^2
#else
    call lsquit('xcfun not activated -DVAR_XCFUN (can only be done using cmake)',-1)
#endif
  end subroutine xcfun2_gga_xc_single_eval

  subroutine xcfun3_gga_xc_single_eval(XCFUNINPUT,XCFUNOUTPUT)
    implicit none
    REAL(REALK),intent(in) :: XCFUNINPUT(2,1)
    REAL(REALK),intent(inout) :: XCFUNOUTPUT(10,1)
    !
    integer :: ierrGGA
#ifdef VAR_XCFUN
    !rho = XCFUNINPUT(1,1) 
    !|grad|^2 = gg = XCFUNINPUT(2,1) 
    ierrGGA = xc_eval_setup(XCFUNfunctional,XC_N_GNN,XC_PARTIAL_DERIVATIVES,3)
    IF(ierrGGA.NE.0) then
       print*,'ierrGGA',ierrGGA
       call lsquit('fun3 too large1',-1)
    endif
    call xc_eval(XCFUNfunctional,1,XCFUNINPUT,XCFUNOUTPUT)
    ierrGGA = xc_eval_setup(XCFUNfunctional,XC_N_GNN,XC_PARTIAL_DERIVATIVES,1)
    IF(ierrGGA.NE.0) call lsquit('fun3 too large2',-1)
    ! XCFUNOUTPUT(1,1) - Exc
    ! XCFUNOUTPUT(2,1) - d Exc/d rho
    ! XCFUNOUTPUT(3,1) - d Exc/d gg
    ! XCFUNOUTPUT(4,1) - d^2 Exc/d rho^2
    ! XCFUNOUTPUT(5,1) - d^2 Exc/d rho d gg
    ! XCFUNOUTPUT(6,1) - d^2 Exc/d gg^2
    ! XCFUNOUTPUT(7,1) - d^3 Exc/d rho^3
    ! XCFUNOUTPUT(8,1) - d^3 Exc/d rho^2 d gg
    ! XCFUNOUTPUT(9,1) - d^3 Exc/d rho gg^2
    ! XCFUNOUTPUT(10,1) -d^3 Exc/d gg^3
#else
    call lsquit('xcfun not activated -DVAR_XCFUN (can only be done using cmake)',-1)
#endif
  end subroutine xcfun3_gga_xc_single_eval

  subroutine xcfun_lda_xc_single_eval(XCFUNINPUT,XCFUNOUTPUT)
    implicit none
    REAL(REALK),intent(in) :: XCFUNINPUT(1,1)
    REAL(REALK),intent(inout) :: XCFUNOUTPUT(2,1)
#ifdef VAR_XCFUN
    !rho = XCFUNINPUT(1,1) 
    call xc_eval(XCFUNfunctional,1,XCFUNINPUT,XCFUNOUTPUT)
    ! XCFUNOUTPUT(1,1) - Exc
    ! XCFUNOUTPUT(2,1) - d Exc/d rho
#else
    call lsquit('xcfun not activated -DVAR_XCFUN (can only be done using cmake)',-1)
#endif
  end subroutine xcfun_lda_xc_single_eval

  subroutine xcfun2_lda_xc_single_eval(XCFUNINPUT,XCFUNOUTPUT)
    implicit none
    REAL(REALK),intent(in) :: XCFUNINPUT(1,1)
    REAL(REALK),intent(inout) :: XCFUNOUTPUT(3,1)
    integer :: ierrLDA
#ifdef VAR_XCFUN
    !rho = XCFUNINPUT(1,1) 
    ierrLDA = xc_eval_setup(XCFUNfunctional,XC_N,XC_PARTIAL_DERIVATIVES,2)
    call xc_eval(XCFUNfunctional,1,XCFUNINPUT,XCFUNOUTPUT)
    ierrLDA = xc_eval_setup(XCFUNfunctional,XC_N,XC_PARTIAL_DERIVATIVES,1)
    ! XCFUNOUTPUT(1,1) - Exc
    ! XCFUNOUTPUT(2,1) - d Exc/d rho
    ! XCFUNOUTPUT(3,1) - d^2 Exc/d rho^2
#else
    call lsquit('xcfun not activated -DVAR_XCFUN (can only be done using cmake)',-1)
#endif
  end subroutine xcfun2_lda_xc_single_eval

  subroutine xcfun_gga_unres_xc_single_eval(XCFUNINPUT,XCFUNOUTPUT)
    implicit none
    REAL(REALK),intent(in) :: XCFUNINPUT(5,1)
    REAL(REALK),intent(inout) :: XCFUNOUTPUT(6,1)
#ifdef VAR_XCFUN
    !na = rho_alpha = XCFUNINPUT(1,1) 
    !nb = rho_beta = XCFUNINPUT(2,1) 
    !gaa = |grad_alpha|^2 = XCFUNINPUT(3,1) 
    !gab = grad_alpha*grad_beta = XCFUNINPUT(4,1) 
    !gbb = |grad_beta|^2 = XCFUNINPUT(5,1) 
    call xc_eval(XCFUNfunctional,1,XCFUNINPUT,XCFUNOUTPUT)
    ! XCFUNOUTPUT(1,1) - Exc
    ! XCFUNOUTPUT(2,1) - d Exc/d na
    ! XCFUNOUTPUT(3,1) - d Exc/d nb
    ! XCFUNOUTPUT(4,1) - d Exc/d gaa
    ! XCFUNOUTPUT(5,1) - d Exc/d gab
    ! XCFUNOUTPUT(6,1) - d Exc/d gbb
#else
    call lsquit('xcfun not activated -DVAR_XCFUN (can only be done using cmake)',-1)
#endif
  end subroutine xcfun_gga_unres_xc_single_eval

  subroutine xcfun_lda_unres_xc_single_eval(XCFUNINPUT,XCFUNOUTPUT)
    implicit none
    REAL(REALK),intent(in) :: XCFUNINPUT(2,1)
    REAL(REALK),intent(inout) :: XCFUNOUTPUT(3,1)
#ifdef VAR_XCFUN
    !na = rho_alpha = XCFUNINPUT(1,1) 
    !nb = rho_beta = XCFUNINPUT(2,1) 
    call xc_eval(XCFUNfunctional,1,XCFUNINPUT,XCFUNOUTPUT)
#else
    call lsquit('xcfun not activated -DVAR_XCFUN (can only be done using cmake)',-1)
#endif
    ! XCFUNOUTPUT(1,1) - Exc
    ! XCFUNOUTPUT(2,1) - d Exc/d na
    ! XCFUNOUTPUT(3,1) - d Exc/d nb
  end subroutine xcfun_lda_unres_xc_single_eval

  subroutine xcfun_host_free()
    implicit none
#ifdef VAR_XCFUN
    call xc_free_functional(XCFUNfunctional)
#else
    call lsquit('xcfun not activated -DVAR_XCFUN (can only be done using cmake)',-1)
#endif
  end subroutine xcfun_host_free

  SUBROUTINE determine_nStrings(string,n)
    implicit none
    character(len=*), intent(in) :: string
    integer, intent(inout) :: n !number of found words on output
    integer :: i,imax,j,k,nspaces
    logical :: inword
    inword = .false.
    imax = LEN(string)
    nspaces = 0
    i = 0  !string index
    j = 0  !number of words found
    do
       i = i + 1
       if (i > imax) then
          n = j
          exit
       endif
       if (string(i:i) == ' ') then
          inword = .false.
          nspaces = nspaces + 1
          cycle
       else
          nspaces = 0
          if (.not. inword) then
             !found beginning of word
             j = j + 1
             k = 1 !word index
             inword = .true.
          else
             k = k + 1
          endif
       endif
    enddo    
  END SUBROUTINE DETERMINE_NSTRINGS

  SUBROUTINE TRIM_STRINGS(string,n,words,WeightSingle)
    implicit none
    character(len=*), intent(in) :: string
    integer, intent(inout) :: n  !max number of wanted words on input
                               !number of found words on output
    character(len=80), intent(inout) :: words(n)
    real(realk),intent(inout) :: WeightSingle(n)
    !
    integer :: i,imax,j,k,nspaces,nweight,jj
    logical :: inword
    character(len=80) :: TMPweight
    do i=1,n
       do j=1,LEN(words(i))
          words(i)(j:j) = ' '
       enddo
    enddo
    inword = .false.
    imax = LEN(string)
    nspaces = 0
    i = 0  !string index
    j = 0  !number of words found
    MAJOR: do
       i = i + 1
       if (i > imax) then
          print*,'Finished searching the string:',string, &
               & 'found ',j,' words out of ',n
          call lsquit('Finished searching the string number mismatch',-1)
       endif
       if (string(i:i) == '='.AND. inword) then
          print*,'found func =',words(j)
          do jj=1,LEN(TMPweight)
             TMPweight(jj:jj) = ' '
          enddo
          inword = .false.
          nspaces = nspaces + 1
          nweight = 0
          MINOR: do
             i = i + 1
             if (i > imax) then
                print*,'Finished searching the string:',string, &
                     & 'found ',j,' words out of ',n
                call lsquit('Finished searching the string number mismatch',-1)
             endif
             if (string(i:i) == ' ') then
                TMPweight(1:nweight) = string(i-nweight+1:i)
                READ(TMPweight,*) WeightSingle(j)
                print*,'found weight =',WeightSingle(j)
                if (j == n) then
                   exit MAJOR !found all words
                endif
                cycle MAJOR
             else
                nweight = nweight + 1                
             endif
          enddo MINOR
       elseif (string(i:i) == ' ') then
          inword = .false.
          nspaces = nspaces + 1
          if (j == n) then
             print*,'found func =',words(j)
             exit MAJOR !found all words
          endif
          cycle MAJOR
       else
          nspaces = 0
          if (.not. inword) then
             !found beginning of word
             j = j + 1
             k = 1 !word index
             inword = .true.
          else
             k = k + 1
          endif
          words(j)(k:k) = string(i:i)
       endif
    enddo MAJOR

  END SUBROUTINE TRIM_STRINGS

end module xcfun_host
